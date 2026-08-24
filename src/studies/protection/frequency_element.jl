struct ROCOFEstimatorSettings
    id::Symbol
    frequency_measurement_id::Symbol
    frequency_channel::Symbol
    sample_period_ticks::Int
    window_intervals::Int
    minimum_frequency_hz::Float64
    maximum_frequency_hz::Float64
    orientation::String
    numerical_provenance::ParameterProvenance

    function ROCOFEstimatorSettings(
        id::Symbol,
        frequency_measurement_id::Symbol;
        frequency_channel::Symbol=:frequency,
        sample_period_ticks::Integer,
        window_intervals::Integer,
        minimum_frequency_hz::Real,
        maximum_frequency_hz::Real,
        orientation::AbstractString,
        numerical_provenance::ParameterProvenance,
    )
        all(name -> !isempty(String(name)), (
            id,
            frequency_measurement_id,
            frequency_channel,
        )) || throw(ArgumentError("ROCOF estimator identities must not be empty"))
        period = Int(sample_period_ticks)
        window = Int(window_intervals)
        period > 0 || throw(ArgumentError(
            "ROCOF sample period must be a positive tick count",
        ))
        window > 0 || throw(ArgumentError(
            "ROCOF window must contain at least one past interval",
        ))
        frequency_bounds = Float64.((minimum_frequency_hz, maximum_frequency_hz))
        all(isfinite, frequency_bounds) &&
            0.0 < frequency_bounds[1] < frequency_bounds[2] ||
            throw(ArgumentError("ROCOF frequency bounds must be finite positive and ordered"))
        orientation_string = String(orientation)
        isempty(strip(orientation_string)) && throw(ArgumentError(
            "ROCOF orientation must not be empty",
        ))
        numerical_provenance.nature === NumericalPolicyParameter ||
            throw(ArgumentError("ROCOF window provenance must be numerical policy"))
        return new(
            id,
            frequency_measurement_id,
            frequency_channel,
            period,
            window,
            frequency_bounds...,
            orientation_string,
            numerical_provenance,
        )
    end
end

mutable struct ROCOFEstimatorState
    source_ticks::Vector{Int}
    frequencies_hz::Vector{Float64}
    last_release_tick::Int
    evaluation_count::Int
    estimate_count::Int
    blocked_count::Int

    ROCOFEstimatorState() = new(Int[], Float64[], -1, 0, 0, 0)
end

struct ROCOFEstimatorSnapshot
    settings_signature_sha256::String
    source_ticks::Vector{Int}
    frequencies_hz::Vector{Float64}
    last_release_tick::Int
    evaluation_count::Int
    estimate_count::Int
    blocked_count::Int
end

function _rocof_estimator_signature(settings::ROCOFEstimatorSettings)
    io = IOBuffer()
    for field in fieldnames(ROCOFEstimatorSettings)
        print(io, field, '=', repr(getfield(settings, field)), '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

function rocof_estimator_snapshot(
    state::ROCOFEstimatorState,
    settings::ROCOFEstimatorSettings,
)
    return ROCOFEstimatorSnapshot(
        _rocof_estimator_signature(settings),
        copy(state.source_ticks),
        copy(state.frequencies_hz),
        state.last_release_tick,
        state.evaluation_count,
        state.estimate_count,
        state.blocked_count,
    )
end

function restore_rocof_estimator_snapshot!(
    state::ROCOFEstimatorState,
    settings::ROCOFEstimatorSettings,
    snapshot::ROCOFEstimatorSnapshot,
)
    snapshot.settings_signature_sha256 == _rocof_estimator_signature(settings) ||
        throw(ArgumentError("ROCOF snapshot settings identity is stale"))
    length(snapshot.source_ticks) == length(snapshot.frequencies_hz) ||
        throw(ArgumentError("ROCOF snapshot history lengths differ"))
    length(snapshot.source_ticks) <= settings.window_intervals + 1 ||
        throw(ArgumentError("ROCOF snapshot exceeds its declared window"))
    issorted(snapshot.source_ticks) && allunique(snapshot.source_ticks) ||
        throw(ArgumentError("ROCOF snapshot source ticks are not strictly ordered"))
    all(isfinite, snapshot.frequencies_hz) || throw(ArgumentError(
        "ROCOF snapshot frequency history must be finite",
    ))
    snapshot.last_release_tick >= -1 && all(>=(0), (
        snapshot.evaluation_count,
        snapshot.estimate_count,
        snapshot.blocked_count,
    )) || throw(ArgumentError("ROCOF snapshot counters are invalid"))
    state.source_ticks = copy(snapshot.source_ticks)
    state.frequencies_hz = copy(snapshot.frequencies_hz)
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.estimate_count = snapshot.estimate_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function _rocof_output(
    settings::ROCOFEstimatorSettings,
    measurement::ProtectionMeasurement,
    value,
    reason,
)
    owner_signature = bytes2hex(sha256(
        _rocof_estimator_signature(settings) *
        measurement.deterministic_signature_sha256,
    ))
    return ProtectionMeasurement(
        settings.id,
        :rocof,
        :rocof,
        "Hz/s",
        settings.orientation,
        ProtectionROCOFStage,
        measurement.source_tick,
        measurement.release_tick,
        measurement.tick_s,
        value,
        reason === nothing ? :valid : :history_unavailable,
        reason,
        owner_signature,
    )
end

"""Estimate a causal past-window slope and invalidate history across missing samples."""
function estimate_rocof!(
    state::ROCOFEstimatorState,
    settings::ROCOFEstimatorSettings,
    measurement::ProtectionMeasurement,
)
    measurement.measurement_id == settings.frequency_measurement_id ||
        throw(ArgumentError("ROCOF frequency owner identity does not match its settings"))
    measurement.channel == settings.frequency_channel || throw(ArgumentError(
        "ROCOF frequency channel identity does not match its settings",
    ))
    measurement.stage === ProtectionFrequencyStage || throw(ArgumentError(
        "ROCOF input must be an accepted A200 frequency-stage measurement",
    ))
    measurement.unit == "Hz" || throw(ArgumentError(
        "ROCOF input must use hertz",
    ))
    measurement.release_tick > state.last_release_tick || throw(ArgumentError(
        "ROCOF measurements must advance in strict release-tick order",
    ))
    state.last_release_tick = measurement.release_tick
    state.evaluation_count += 1
    if measurement.value === nothing
        empty!(state.source_ticks)
        empty!(state.frequencies_hz)
        state.blocked_count += 1
        return _rocof_output(settings, measurement, nothing, measurement.unavailable_reason)
    end
    value = measurement.value
    imaginary_tolerance = 64.0 * eps(Float64) * max(1.0, abs(real(value)))
    abs(imag(value)) <= imaginary_tolerance || throw(ArgumentError(
        "ROCOF frequency input must be a real scalar",
    ))
    frequency = real(value)
    settings.minimum_frequency_hz <= frequency <= settings.maximum_frequency_hz || begin
        empty!(state.source_ticks)
        empty!(state.frequencies_hz)
        state.blocked_count += 1
        return _rocof_output(settings, measurement, nothing, :frequency_out_of_domain)
    end
    if !isempty(state.source_ticks) &&
       measurement.source_tick - last(state.source_ticks) != settings.sample_period_ticks
        empty!(state.source_ticks)
        empty!(state.frequencies_hz)
        state.blocked_count += 1
    end
    push!(state.source_ticks, measurement.source_tick)
    push!(state.frequencies_hz, frequency)
    maximum_count = settings.window_intervals + 1
    length(state.source_ticks) > maximum_count && popfirst!(state.source_ticks)
    length(state.frequencies_hz) > maximum_count && popfirst!(state.frequencies_hz)
    if length(state.source_ticks) < maximum_count
        return _rocof_output(settings, measurement, nothing, :rocof_window_incomplete)
    end
    elapsed_ticks = last(state.source_ticks) - first(state.source_ticks)
    elapsed_ticks == settings.window_intervals * settings.sample_period_ticks ||
        throw(ArgumentError("ROCOF accepted history is not on its declared exact calendar"))
    estimate = (last(state.frequencies_hz) - first(state.frequencies_hz)) /
        (elapsed_ticks * measurement.tick_s)
    state.estimate_count += 1
    return _rocof_output(settings, measurement, estimate, nothing)
end

struct DirectionalRelaySettings
    id::Symbol
    voltage_measurement_id::Symbol
    voltage_channel::Symbol
    current_measurement_id::Symbol
    current_channel::Symbol
    characteristic_angle_rad::Float64
    minimum_polarizing_voltage_v::Float64
    minimum_operating_torque_w::Float64
    forward_polarity::Int
    memory_decay_per_tick::Float64
    maximum_memory_age_ticks::Int
    provenance::ParameterProvenance

    function DirectionalRelaySettings(
        id::Symbol,
        voltage_measurement_id::Symbol,
        voltage_channel::Symbol,
        current_measurement_id::Symbol,
        current_channel::Symbol;
        characteristic_angle_rad::Real,
        minimum_polarizing_voltage_v::Real,
        minimum_operating_torque_w::Real=0.0,
        forward_polarity::Integer=1,
        memory_decay_per_tick::Real=1.0,
        maximum_memory_age_ticks::Integer=0,
        provenance::ParameterProvenance,
    )
        all(name -> !isempty(String(name)), (
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
        )) || throw(ArgumentError("directional relay identities must not be empty"))
        angle, minimum_voltage, minimum_torque, decay = Float64.((
            characteristic_angle_rad,
            minimum_polarizing_voltage_v,
            minimum_operating_torque_w,
            memory_decay_per_tick,
        ))
        all(isfinite, (angle, minimum_voltage, minimum_torque, decay)) ||
            throw(ArgumentError("directional relay settings must be finite"))
        minimum_voltage >= 0.0 && minimum_torque >= 0.0 || throw(ArgumentError(
            "directional relay minima must be nonnegative",
        ))
        0.0 < decay <= 1.0 || throw(ArgumentError(
            "directional relay memory decay must be in (0, 1]",
        ))
        polarity = Int(forward_polarity)
        polarity in (-1, 1) || throw(ArgumentError(
            "directional relay forward polarity must be -1 or 1",
        ))
        maximum_age = Int(maximum_memory_age_ticks)
        maximum_age >= 0 || throw(ArgumentError(
            "directional relay memory age must be nonnegative",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "directional relay setting provenance must be physical",
        ))
        return new(
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
            angle,
            minimum_voltage,
            minimum_torque,
            polarity,
            decay,
            maximum_age,
            provenance,
        )
    end
end

mutable struct DirectionalRelayState
    memory_voltage_v::Union{Nothing,ComplexF64}
    memory_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int

    DirectionalRelayState() = new(nothing, -1, -1, 0, 0)
end

struct DirectionalRelaySnapshot
    settings_signature_sha256::String
    memory_voltage_v::Union{Nothing,ComplexF64}
    memory_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int
end

struct DirectionalRelayDecision
    relay_id::Symbol
    source_tick::Int
    release_tick::Int
    polarizing_voltage_v::Union{Nothing,ComplexF64}
    operating_current_a::Union{Nothing,ComplexF64}
    operating_torque_w::Union{Nothing,Float64}
    forward::Bool
    used_memory_polarization::Bool
    blocked_reason::Union{Nothing,Symbol}
end

function _validate_synchronized_phasors(
    voltage::ProtectionMeasurement,
    current::ProtectionMeasurement,
)
    voltage.source_tick == current.source_tick &&
        voltage.release_tick == current.release_tick &&
        voltage.tick_s == current.tick_s || throw(ArgumentError(
            "protection phasors must share exact source and release clocks",
        ))
    voltage.stage in (ProtectionFundamentalPhasorStage, ProtectionSequencePhasorStage) &&
        current.stage in (ProtectionFundamentalPhasorStage, ProtectionSequencePhasorStage) ||
        throw(ArgumentError("directional and impedance inputs must be RMS phasors"))
    voltage.unit == "V" && current.unit == "A" || throw(ArgumentError(
        "directional and impedance phasors must use volts and amperes",
    ))
    return nothing
end

function evaluate_directional_relay!(
    state::DirectionalRelayState,
    settings::DirectionalRelaySettings,
    voltage::ProtectionMeasurement,
    current::ProtectionMeasurement,
)
    _validate_synchronized_phasors(voltage, current)
    voltage.measurement_id == settings.voltage_measurement_id &&
        voltage.channel == settings.voltage_channel &&
        current.measurement_id == settings.current_measurement_id &&
        current.channel == settings.current_channel || throw(ArgumentError(
            "directional relay input identities do not match their settings",
        ))
    voltage.release_tick > state.last_release_tick || throw(ArgumentError(
        "directional relay measurements must advance in strict release-tick order",
    ))
    state.last_release_tick = voltage.release_tick
    state.evaluation_count += 1
    current.value === nothing && begin
        state.blocked_count += 1
        return DirectionalRelayDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            nothing,
            false,
            false,
            current.unavailable_reason,
        )
    end
    present_voltage = voltage.value
    polarizing_voltage = present_voltage
    used_memory = false
    if present_voltage === nothing || abs(present_voltage) < settings.minimum_polarizing_voltage_v
        age = state.memory_source_tick < 0 ? typemax(Int) :
            voltage.source_tick - state.memory_source_tick
        if state.memory_voltage_v === nothing || age < 0 || age > settings.maximum_memory_age_ticks
            state.blocked_count += 1
            reason = present_voltage === nothing ? voltage.unavailable_reason :
                :polarizing_voltage_below_minimum
            return DirectionalRelayDecision(
                settings.id,
                voltage.source_tick,
                voltage.release_tick,
                nothing,
                current.value,
                nothing,
                false,
                false,
                reason,
            )
        end
        polarizing_voltage = state.memory_voltage_v * settings.memory_decay_per_tick^age
        used_memory = true
    else
        state.memory_voltage_v = present_voltage
        state.memory_source_tick = voltage.source_tick
    end
    torque = settings.forward_polarity * real(
        polarizing_voltage * conj(current.value) * cis(-settings.characteristic_angle_rad),
    )
    return DirectionalRelayDecision(
        settings.id,
        voltage.source_tick,
        voltage.release_tick,
        polarizing_voltage,
        current.value,
        torque,
        torque >= settings.minimum_operating_torque_w,
        used_memory,
        nothing,
    )
end

abstract type AbstractDistanceZoneShape end

struct MhoDistanceZone <: AbstractDistanceZoneShape
    center_ohm::ComplexF64
    radius_ohm::Float64

    function MhoDistanceZone(center_ohm::Number, radius_ohm::Real)
        center = ComplexF64(center_ohm)
        radius = Float64(radius_ohm)
        isfinite(real(center)) && isfinite(imag(center)) &&
            isfinite(radius) && radius > 0.0 || throw(ArgumentError(
            "mho zone center must be finite and radius finite positive",
        ))
        return new(center, radius)
    end
end

struct PolygonDistanceZone <: AbstractDistanceZoneShape
    outward_normals::Vector{ComplexF64}
    limits_ohm::Vector{Float64}

    function PolygonDistanceZone(
        outward_normals::AbstractVector{<:Number},
        limits_ohm::AbstractVector{<:Real},
    )
        normals = ComplexF64.(outward_normals)
        limits = Float64.(limits_ohm)
        length(normals) == length(limits) >= 3 || throw(ArgumentError(
            "polygonal distance zone requires at least three matched half-planes",
        ))
        all(value -> isfinite(real(value)) && isfinite(imag(value)) && !iszero(value), normals) &&
            all(isfinite, limits) || throw(ArgumentError(
            "polygonal distance-zone half-planes must be finite with nonzero normals",
        ))
        return new(normals, limits)
    end
end

struct DistanceRelaySettings{Z<:AbstractDistanceZoneShape}
    id::Symbol
    voltage_measurement_id::Symbol
    voltage_channel::Symbol
    current_measurement_id::Symbol
    current_channel::Symbol
    minimum_loop_current_a::Float64
    zone::Z
    provenance::ParameterProvenance

    function DistanceRelaySettings(
        id::Symbol,
        voltage_measurement_id::Symbol,
        voltage_channel::Symbol,
        current_measurement_id::Symbol,
        current_channel::Symbol,
        zone::Z;
        minimum_loop_current_a::Real,
        provenance::ParameterProvenance,
    ) where {Z<:AbstractDistanceZoneShape}
        all(name -> !isempty(String(name)), (
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
        )) || throw(ArgumentError("distance relay identities must not be empty"))
        minimum_current = Float64(minimum_loop_current_a)
        isfinite(minimum_current) && minimum_current > 0.0 || throw(ArgumentError(
            "distance relay minimum loop current must be finite and positive",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "distance relay setting provenance must be physical",
        ))
        return new{Z}(
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
            minimum_current,
            zone,
            provenance,
        )
    end
end

mutable struct DistanceRelayState
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int

    DistanceRelayState() = new(-1, 0, 0)
end

struct DistanceRelaySnapshot
    settings_signature_sha256::String
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int
end

struct DistanceRelayDecision
    relay_id::Symbol
    source_tick::Int
    release_tick::Int
    apparent_impedance_ohm::Union{Nothing,ComplexF64}
    zone_margin_ohm::Union{Nothing,Float64}
    asserted::Bool
    blocked_reason::Union{Nothing,Symbol}
end

_distance_zone_margin(zone::MhoDistanceZone, impedance::ComplexF64) =
    zone.radius_ohm - abs(impedance - zone.center_ohm)
function _distance_zone_margin(zone::PolygonDistanceZone, impedance::ComplexF64)
    return minimum(
        limit - real(conj(normal) * impedance)
        for (normal, limit) in zip(zone.outward_normals, zone.limits_ohm)
    )
end

function evaluate_distance_relay!(
    state::DistanceRelayState,
    settings::DistanceRelaySettings,
    voltage::ProtectionMeasurement,
    current::ProtectionMeasurement,
)
    _validate_synchronized_phasors(voltage, current)
    voltage.measurement_id == settings.voltage_measurement_id &&
        voltage.channel == settings.voltage_channel &&
        current.measurement_id == settings.current_measurement_id &&
        current.channel == settings.current_channel || throw(ArgumentError(
            "distance relay input identities do not match their settings",
        ))
    voltage.release_tick > state.last_release_tick || throw(ArgumentError(
        "distance relay measurements must advance in strict release-tick order",
    ))
    state.last_release_tick = voltage.release_tick
    state.evaluation_count += 1
    if voltage.value === nothing || current.value === nothing
        state.blocked_count += 1
        reason = voltage.value === nothing ? voltage.unavailable_reason : current.unavailable_reason
        return DistanceRelayDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            false,
            reason,
        )
    end
    if abs(current.value) < settings.minimum_loop_current_a
        state.blocked_count += 1
        return DistanceRelayDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            false,
            :loop_current_below_minimum,
        )
    end
    impedance = voltage.value / current.value
    margin = _distance_zone_margin(settings.zone, impedance)
    return DistanceRelayDecision(
        settings.id,
        voltage.source_tick,
        voltage.release_tick,
        impedance,
        margin,
        margin >= 0.0,
        nothing,
    )
end

@enum DifferentialRestraintMode::UInt8 begin
    HalfSumDifferentialRestraint = 0x01
    MaximumDifferentialRestraint = 0x02
end

struct DifferentialRelaySettings
    id::Symbol
    terminal_measurement_ids::Vector{Symbol}
    terminal_channels::Vector{Symbol}
    compensation::Vector{ComplexF64}
    restraint_mode::DifferentialRestraintMode
    minimum_operate_a::Float64
    initial_bias_a::Float64
    restraint_breakpoints_a::Vector{Float64}
    region_slopes::Vector{Float64}
    provenance::ParameterProvenance

    function DifferentialRelaySettings(
        id::Symbol,
        terminal_measurement_ids::AbstractVector{Symbol},
        terminal_channels::AbstractVector{Symbol};
        compensation::AbstractVector{<:Number},
        restraint_mode::DifferentialRestraintMode=HalfSumDifferentialRestraint,
        minimum_operate_a::Real,
        initial_bias_a::Real=0.0,
        restraint_breakpoints_a::AbstractVector{<:Real}=Float64[],
        region_slopes::AbstractVector{<:Real},
        provenance::ParameterProvenance,
    )
        ids = Symbol.(terminal_measurement_ids)
        channels = Symbol.(terminal_channels)
        factors = ComplexF64.(compensation)
        terminal_count = length(ids)
        terminal_count >= 2 && length(channels) == terminal_count &&
            length(factors) == terminal_count || throw(ArgumentError(
            "differential relay requires at least two matched terminals and compensation factors",
        ))
        all(name -> !isempty(String(name)), vcat([id], ids, channels)) ||
            throw(ArgumentError("differential relay identities must not be empty"))
        all(value -> isfinite(real(value)) && isfinite(imag(value)) && !iszero(value), factors) ||
            throw(ArgumentError("differential compensation factors must be finite and nonzero"))
        minimum_operate = Float64(minimum_operate_a)
        initial_bias = Float64(initial_bias_a)
        breakpoints = Float64.(restraint_breakpoints_a)
        slopes = Float64.(region_slopes)
        isfinite(minimum_operate) && minimum_operate >= 0.0 &&
            isfinite(initial_bias) && initial_bias >= 0.0 || throw(ArgumentError(
            "differential operate minimum and bias must be finite and nonnegative",
        ))
        length(slopes) == length(breakpoints) + 1 && !isempty(slopes) ||
            throw(ArgumentError("differential slopes must define every bias region"))
        all(isfinite, breakpoints) && all(>=(0.0), breakpoints) && issorted(breakpoints) &&
            allunique(breakpoints) || throw(ArgumentError(
            "differential restraint breakpoints must be finite nonnegative and strictly ordered",
        ))
        all(value -> isfinite(value) && value >= 0.0, slopes) || throw(ArgumentError(
            "differential slopes must be finite and nonnegative",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "differential relay setting provenance must be physical",
        ))
        return new(
            id,
            ids,
            channels,
            factors,
            restraint_mode,
            minimum_operate,
            initial_bias,
            breakpoints,
            slopes,
            provenance,
        )
    end
end

mutable struct DifferentialRelayState
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int

    DifferentialRelayState() = new(-1, 0, 0)
end

struct DifferentialRelaySnapshot
    settings_signature_sha256::String
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int
end

struct DifferentialRelayDecision
    relay_id::Symbol
    source_tick::Int
    release_tick::Int
    compensated_currents_a::Union{Nothing,Vector{ComplexF64}}
    operate_current_a::Union{Nothing,Float64}
    restraint_current_a::Union{Nothing,Float64}
    threshold_current_a::Union{Nothing,Float64}
    margin_a::Union{Nothing,Float64}
    asserted::Bool
    blocked_reason::Union{Nothing,Symbol}
end

function _differential_bias_threshold(settings::DifferentialRelaySettings, restraint::Float64)
    threshold = settings.initial_bias_a
    lower = 0.0
    for (region, slope) in pairs(settings.region_slopes)
        upper = region <= length(settings.restraint_breakpoints_a) ?
            settings.restraint_breakpoints_a[region] : restraint
        width = max(0.0, min(restraint, upper) - lower)
        threshold += slope * width
        restraint <= upper && break
        lower = upper
    end
    return max(settings.minimum_operate_a, threshold)
end

function evaluate_differential_relay!(
    state::DifferentialRelayState,
    settings::DifferentialRelaySettings,
    measurements::AbstractVector{ProtectionMeasurement},
)
    length(measurements) == length(settings.terminal_measurement_ids) ||
        throw(DimensionMismatch("differential relay terminal measurement count changed"))
    first_measurement = first(measurements)
    all(measurement -> measurement.source_tick == first_measurement.source_tick &&
        measurement.release_tick == first_measurement.release_tick &&
        measurement.tick_s == first_measurement.tick_s,
        measurements,
    ) || throw(ArgumentError("differential terminal measurements must be exactly synchronized"))
    first_measurement.release_tick > state.last_release_tick || throw(ArgumentError(
        "differential relay measurements must advance in strict release-tick order",
    ))
    for index in eachindex(measurements)
        measurement = measurements[index]
        measurement.measurement_id == settings.terminal_measurement_ids[index] &&
            measurement.channel == settings.terminal_channels[index] ||
            throw(ArgumentError("differential terminal identity does not match its ordered setting"))
        measurement.stage in (ProtectionFundamentalPhasorStage, ProtectionSequencePhasorStage) &&
            measurement.unit == "A" || throw(ArgumentError(
            "differential terminal input must be an RMS current phasor",
        ))
    end
    state.last_release_tick = first_measurement.release_tick
    state.evaluation_count += 1
    unavailable_index = findfirst(measurement -> measurement.value === nothing, measurements)
    if unavailable_index !== nothing
        state.blocked_count += 1
        return DifferentialRelayDecision(
            settings.id,
            first_measurement.source_tick,
            first_measurement.release_tick,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            false,
            measurements[unavailable_index].unavailable_reason,
        )
    end
    compensated = ComplexF64[
        something(measurement.value) * settings.compensation[index]
        for (index, measurement) in pairs(measurements)
    ]
    operate = abs(sum(compensated))
    restraint = settings.restraint_mode === HalfSumDifferentialRestraint ?
        sum(abs, compensated) / 2.0 : maximum(abs, compensated)
    threshold = _differential_bias_threshold(settings, restraint)
    margin = operate - threshold
    return DifferentialRelayDecision(
        settings.id,
        first_measurement.source_tick,
        first_measurement.release_tick,
        compensated,
        operate,
        restraint,
        threshold,
        margin,
        margin >= 0.0,
        nothing,
    )
end

@enum IncrementalWaveDirection::UInt8 begin
    ForwardIncrementalWave = 0x01
    ReverseIncrementalWave = 0x02
end

struct IncrementalWaveSettings
    id::Symbol
    voltage_measurement_id::Symbol
    voltage_channel::Symbol
    current_measurement_id::Symbol
    current_channel::Symbol
    sample_period_ticks::Int
    reference_impedance_ohm::Float64
    direction::IncrementalWaveDirection
    polarity::Int
    threshold_v::Float64
    provenance::ParameterProvenance

    function IncrementalWaveSettings(
        id::Symbol,
        voltage_measurement_id::Symbol,
        voltage_channel::Symbol,
        current_measurement_id::Symbol,
        current_channel::Symbol;
        sample_period_ticks::Integer,
        reference_impedance_ohm::Real,
        direction::IncrementalWaveDirection,
        polarity::Integer=1,
        threshold_v::Real,
        provenance::ParameterProvenance,
    )
        all(name -> !isempty(String(name)), (
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
        )) || throw(ArgumentError("incremental-wave identities must not be empty"))
        period = Int(sample_period_ticks)
        period > 0 || throw(ArgumentError(
            "incremental-wave sample period must be a positive tick count",
        ))
        impedance, threshold = Float64.((reference_impedance_ohm, threshold_v))
        isfinite(impedance) && impedance > 0.0 && isfinite(threshold) && threshold > 0.0 ||
            throw(ArgumentError("incremental-wave impedance and threshold must be finite positive"))
        sign = Int(polarity)
        sign in (-1, 1) || throw(ArgumentError(
            "incremental-wave polarity must be -1 or 1",
        ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "incremental-wave setting provenance must be physical",
        ))
        return new(
            id,
            voltage_measurement_id,
            voltage_channel,
            current_measurement_id,
            current_channel,
            period,
            impedance,
            direction,
            sign,
            threshold,
            provenance,
        )
    end
end

mutable struct IncrementalWaveState
    previous_voltage_v::Union{Nothing,Float64}
    previous_current_a::Union{Nothing,Float64}
    previous_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int

    IncrementalWaveState() = new(nothing, nothing, -1, -1, 0, 0)
end

struct IncrementalWaveSnapshot
    settings_signature_sha256::String
    previous_voltage_v::Union{Nothing,Float64}
    previous_current_a::Union{Nothing,Float64}
    previous_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    blocked_count::Int
end

struct IncrementalWaveDecision
    relay_id::Symbol
    source_tick::Int
    release_tick::Int
    forward_wave_v::Union{Nothing,Float64}
    reverse_wave_v::Union{Nothing,Float64}
    selected_wave_v::Union{Nothing,Float64}
    asserted::Bool
    blocked_reason::Union{Nothing,Symbol}
end

function _protection_element_settings_signature(settings)
    io = IOBuffer()
    for field in fieldnames(typeof(settings))
        print(io, field, '=', repr(getfield(settings, field)), '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

function directional_relay_snapshot(
    state::DirectionalRelayState,
    settings::DirectionalRelaySettings,
)
    return DirectionalRelaySnapshot(
        _protection_element_settings_signature(settings),
        state.memory_voltage_v,
        state.memory_source_tick,
        state.last_release_tick,
        state.evaluation_count,
        state.blocked_count,
    )
end

function restore_directional_relay_snapshot!(
    state::DirectionalRelayState,
    settings::DirectionalRelaySettings,
    snapshot::DirectionalRelaySnapshot,
)
    snapshot.settings_signature_sha256 ==
        _protection_element_settings_signature(settings) || throw(ArgumentError(
        "directional relay snapshot settings identity is stale",
    ))
    memory_valid = snapshot.memory_voltage_v === nothing ?
        snapshot.memory_source_tick == -1 :
        snapshot.memory_source_tick >= 0 &&
            isfinite(real(snapshot.memory_voltage_v)) &&
            isfinite(imag(snapshot.memory_voltage_v))
    memory_valid && snapshot.last_release_tick >= -1 &&
        snapshot.memory_source_tick <= snapshot.last_release_tick &&
        all(>=(0), (snapshot.evaluation_count, snapshot.blocked_count)) &&
        snapshot.blocked_count <= snapshot.evaluation_count || throw(ArgumentError(
        "directional relay snapshot state is invalid",
    ))
    state.memory_voltage_v = snapshot.memory_voltage_v
    state.memory_source_tick = snapshot.memory_source_tick
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function distance_relay_snapshot(
    state::DistanceRelayState,
    settings::DistanceRelaySettings,
)
    return DistanceRelaySnapshot(
        _protection_element_settings_signature(settings),
        state.last_release_tick,
        state.evaluation_count,
        state.blocked_count,
    )
end

function restore_distance_relay_snapshot!(
    state::DistanceRelayState,
    settings::DistanceRelaySettings,
    snapshot::DistanceRelaySnapshot,
)
    snapshot.settings_signature_sha256 ==
        _protection_element_settings_signature(settings) || throw(ArgumentError(
        "distance relay snapshot settings identity is stale",
    ))
    snapshot.last_release_tick >= -1 &&
        all(>=(0), (snapshot.evaluation_count, snapshot.blocked_count)) &&
        snapshot.blocked_count <= snapshot.evaluation_count || throw(ArgumentError(
        "distance relay snapshot state is invalid",
    ))
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function differential_relay_snapshot(
    state::DifferentialRelayState,
    settings::DifferentialRelaySettings,
)
    return DifferentialRelaySnapshot(
        _protection_element_settings_signature(settings),
        state.last_release_tick,
        state.evaluation_count,
        state.blocked_count,
    )
end

function restore_differential_relay_snapshot!(
    state::DifferentialRelayState,
    settings::DifferentialRelaySettings,
    snapshot::DifferentialRelaySnapshot,
)
    snapshot.settings_signature_sha256 ==
        _protection_element_settings_signature(settings) || throw(ArgumentError(
        "differential relay snapshot settings identity is stale",
    ))
    snapshot.last_release_tick >= -1 &&
        all(>=(0), (snapshot.evaluation_count, snapshot.blocked_count)) &&
        snapshot.blocked_count <= snapshot.evaluation_count || throw(ArgumentError(
        "differential relay snapshot state is invalid",
    ))
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function incremental_wave_snapshot(
    state::IncrementalWaveState,
    settings::IncrementalWaveSettings,
)
    return IncrementalWaveSnapshot(
        _protection_element_settings_signature(settings),
        state.previous_voltage_v,
        state.previous_current_a,
        state.previous_source_tick,
        state.last_release_tick,
        state.evaluation_count,
        state.blocked_count,
    )
end

function restore_incremental_wave_snapshot!(
    state::IncrementalWaveState,
    settings::IncrementalWaveSettings,
    snapshot::IncrementalWaveSnapshot,
)
    snapshot.settings_signature_sha256 ==
        _protection_element_settings_signature(settings) || throw(ArgumentError(
        "incremental-wave snapshot settings identity is stale",
    ))
    history_present = snapshot.previous_voltage_v !== nothing &&
        snapshot.previous_current_a !== nothing
    history_absent = snapshot.previous_voltage_v === nothing &&
        snapshot.previous_current_a === nothing
    history_valid = history_absent ? snapshot.previous_source_tick == -1 :
        history_present && snapshot.previous_source_tick >= 0 &&
            isfinite(snapshot.previous_voltage_v) &&
            isfinite(snapshot.previous_current_a)
    history_valid && snapshot.last_release_tick >= -1 &&
        snapshot.previous_source_tick <= snapshot.last_release_tick &&
        all(>=(0), (snapshot.evaluation_count, snapshot.blocked_count)) &&
        snapshot.blocked_count <= snapshot.evaluation_count || throw(ArgumentError(
        "incremental-wave snapshot state is invalid",
    ))
    state.previous_voltage_v = snapshot.previous_voltage_v
    state.previous_current_a = snapshot.previous_current_a
    state.previous_source_tick = snapshot.previous_source_tick
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function _real_instantaneous_value(measurement::ProtectionMeasurement, unit::String)
    measurement.stage === ProtectionInstantaneousStage && measurement.unit == unit ||
        throw(ArgumentError("incremental-wave input stage or unit is incompatible"))
    measurement.value === nothing && return nothing
    value = measurement.value
    tolerance = 64.0 * eps(Float64) * max(1.0, abs(real(value)))
    abs(imag(value)) <= tolerance || throw(ArgumentError(
        "incremental-wave input must be a real instantaneous scalar",
    ))
    return real(value)
end

function evaluate_incremental_wave!(
    state::IncrementalWaveState,
    settings::IncrementalWaveSettings,
    voltage::ProtectionMeasurement,
    current::ProtectionMeasurement,
)
    voltage.measurement_id == settings.voltage_measurement_id &&
        voltage.channel == settings.voltage_channel &&
        current.measurement_id == settings.current_measurement_id &&
        current.channel == settings.current_channel || throw(ArgumentError(
            "incremental-wave input identities do not match their settings",
        ))
    voltage.source_tick == current.source_tick &&
        voltage.release_tick == current.release_tick &&
        voltage.tick_s == current.tick_s || throw(ArgumentError(
            "incremental-wave voltage and current must be exactly synchronized",
        ))
    voltage.release_tick > state.last_release_tick || throw(ArgumentError(
        "incremental-wave measurements must advance in strict release-tick order",
    ))
    voltage_value = _real_instantaneous_value(voltage, "V")
    current_value = _real_instantaneous_value(current, "A")
    state.last_release_tick = voltage.release_tick
    state.evaluation_count += 1
    if voltage_value === nothing || current_value === nothing
        state.previous_voltage_v = nothing
        state.previous_current_a = nothing
        state.previous_source_tick = -1
        state.blocked_count += 1
        reason = voltage_value === nothing ? voltage.unavailable_reason : current.unavailable_reason
        return IncrementalWaveDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            nothing,
            false,
            reason,
        )
    end
    if state.previous_voltage_v === nothing || state.previous_current_a === nothing
        state.previous_voltage_v = voltage_value
        state.previous_current_a = current_value
        state.previous_source_tick = voltage.source_tick
        return IncrementalWaveDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            nothing,
            false,
            :incremental_history_unavailable,
        )
    end
    expected_source_tick = state.previous_source_tick + settings.sample_period_ticks
    if voltage.source_tick != expected_source_tick
        state.previous_voltage_v = voltage_value
        state.previous_current_a = current_value
        state.previous_source_tick = voltage.source_tick
        state.blocked_count += 1
        return IncrementalWaveDecision(
            settings.id,
            voltage.source_tick,
            voltage.release_tick,
            nothing,
            nothing,
            nothing,
            false,
            :incremental_calendar_gap,
        )
    end
    delta_voltage = voltage_value - state.previous_voltage_v
    delta_current = current_value - state.previous_current_a
    forward = delta_voltage + settings.reference_impedance_ohm * delta_current
    reverse = delta_voltage - settings.reference_impedance_ohm * delta_current
    selected = settings.direction === ForwardIncrementalWave ? forward : reverse
    oriented = settings.polarity * selected
    state.previous_voltage_v = voltage_value
    state.previous_current_a = current_value
    state.previous_source_tick = voltage.source_tick
    return IncrementalWaveDecision(
        settings.id,
        voltage.source_tick,
        voltage.release_tick,
        forward,
        reverse,
        selected,
        oriented >= settings.threshold_v,
        nothing,
    )
end

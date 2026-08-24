mutable struct MagnitudeRelayState
    pickup_active::Bool
    operated::Bool
    timer_fraction::Float64
    last_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    operation_count::Int
    blocked_count::Int

    MagnitudeRelayState() = new(false, false, 0.0, -1, -1, 0, 0, 0)
end

struct MagnitudeRelayDecision
    relay_id::Symbol
    source_tick::Int
    release_tick::Int
    measured_value::Union{Nothing,Float64}
    pickup_multiple::Union{Nothing,Float64}
    pickup_active::Bool
    operated::Bool
    timer_fraction::Float64
    blocked_reason::Union{Nothing,Symbol}
    deterministic_signature_sha256::String
end

struct MagnitudeRelayCandidate
    state::MagnitudeRelayState
    decision::MagnitudeRelayDecision
end

struct MagnitudeRelaySnapshot
    settings_signature_sha256::String
    pickup_active::Bool
    operated::Bool
    timer_fraction::Float64
    last_source_tick::Int
    last_release_tick::Int
    evaluation_count::Int
    operation_count::Int
    blocked_count::Int
end

function magnitude_relay_snapshot(
    state::MagnitudeRelayState,
    settings::MagnitudeRelaySettings,
)
    return MagnitudeRelaySnapshot(
        magnitude_relay_signature(settings),
        state.pickup_active,
        state.operated,
        state.timer_fraction,
        state.last_source_tick,
        state.last_release_tick,
        state.evaluation_count,
        state.operation_count,
        state.blocked_count,
    )
end

function restore_magnitude_relay_snapshot!(
    state::MagnitudeRelayState,
    settings::MagnitudeRelaySettings,
    snapshot::MagnitudeRelaySnapshot,
)
    snapshot.settings_signature_sha256 == magnitude_relay_signature(settings) ||
        throw(ArgumentError("magnitude relay snapshot settings identity is stale"))
    isfinite(snapshot.timer_fraction) && 0.0 <= snapshot.timer_fraction <= 1.0 ||
        throw(ArgumentError("magnitude relay snapshot timer fraction is invalid"))
    snapshot.last_source_tick >= -1 &&
        snapshot.last_release_tick >= snapshot.last_source_tick &&
        all(>=(0), (
            snapshot.evaluation_count,
            snapshot.operation_count,
            snapshot.blocked_count,
        )) || throw(ArgumentError("magnitude relay snapshot state is invalid"))
    snapshot.operation_count <= snapshot.evaluation_count &&
        snapshot.blocked_count <= snapshot.evaluation_count ||
        throw(ArgumentError("magnitude relay snapshot counters are inconsistent"))
    state.pickup_active = snapshot.pickup_active
    state.operated = snapshot.operated
    state.timer_fraction = snapshot.timer_fraction
    state.last_source_tick = snapshot.last_source_tick
    state.last_release_tick = snapshot.last_release_tick
    state.evaluation_count = snapshot.evaluation_count
    state.operation_count = snapshot.operation_count
    state.blocked_count = snapshot.blocked_count
    return state
end

function _magnitude_relay_decision_signature(
    settings::MagnitudeRelaySettings,
    measurement::ProtectionMeasurement,
    measured_value,
    pickup_multiple,
    state::MagnitudeRelayState,
    blocked_reason,
)
    io = IOBuffer()
    print(
        io,
        magnitude_relay_signature(settings),
        '|',
        measurement.deterministic_signature_sha256,
        '|',
        repr(measured_value),
        '|',
        repr(pickup_multiple),
        '|',
        state.pickup_active,
        '|',
        state.operated,
        '|',
        repr(state.timer_fraction),
        '|',
        state.evaluation_count,
        '|',
        state.operation_count,
        '|',
        state.blocked_count,
        '|',
        something(blocked_reason, :available),
    )
    return bytes2hex(sha256(take!(io)))
end

function _validate_magnitude_measurement(
    settings::MagnitudeRelaySettings,
    measurement::ProtectionMeasurement,
)
    measurement.measurement_id == settings.measurement_id || throw(ArgumentError(
        "magnitude relay measurement owner identity does not match its settings",
    ))
    measurement.channel == settings.channel || throw(ArgumentError(
        "magnitude relay channel identity does not match its settings",
    ))
    measurement.stage === settings.stage || throw(ArgumentError(
        "magnitude relay measurement stage does not match its settings",
    ))
    measurement.unit == settings.unit || throw(ArgumentError(
        "magnitude relay measurement unit does not match its settings",
    ))
    return nothing
end

function _magnitude_pickup_state(
    settings::MagnitudeRelaySettings,
    was_picked_up::Bool,
    value::Float64,
)
    if settings.direction === OverMagnitude
        threshold = was_picked_up ? settings.pickup * settings.dropout_ratio :
            settings.pickup
        return value >= threshold
    end
    threshold = was_picked_up ? settings.pickup / settings.dropout_ratio :
        settings.pickup
    return value <= threshold
end

function _magnitude_relay_value(
    settings::MagnitudeRelaySettings,
    measurement::ProtectionMeasurement,
)
    value = something(measurement.value)
    settings.value_mode === AbsoluteMagnitudeValue && return abs(value)
    imaginary_tolerance = 64.0 * eps(Float64) * max(1.0, abs(real(value)))
    abs(imag(value)) <= imaginary_tolerance || throw(ArgumentError(
        "signed magnitude relay input must be a real scalar measurement",
    ))
    return settings.orientation_polarity * real(value)
end

function _magnitude_pickup_multiple(
    settings::MagnitudeRelaySettings,
    value::Float64,
)
    settings.direction === OverMagnitude && return value / settings.pickup
    iszero(value) && return Inf
    return settings.pickup / value
end

function _advance_magnitude_timer(
    settings::MagnitudeRelaySettings,
    prior_fraction::Float64,
    pickup_active::Bool,
    pickup_multiple::Float64,
    elapsed_s::Float64,
)
    if !pickup_active
        settings.reset_time_s == 0.0 && return 0.0
        return max(0.0, prior_fraction - elapsed_s / settings.reset_time_s)
    elseif settings.timer_mode === ProtectionTimerInstantaneous
        return 1.0
    elseif settings.timer_mode === ProtectionDefiniteTimer
        return min(1.0, prior_fraction + elapsed_s / settings.definite_time_s)
    end
    pickup_multiple > 1.0 || return prior_fraction
    operate_time_s = settings.time_dial_s * (
        settings.inverse_a / (pickup_multiple^settings.inverse_p - 1.0) +
        settings.inverse_b
    )
    isfinite(operate_time_s) && operate_time_s > 0.0 || throw(ArgumentError(
        "inverse magnitude relay operate time is outside its positive finite domain",
    ))
    return min(1.0, prior_fraction + elapsed_s / operate_time_s)
end

"""Evaluate one causal candidate without mutating the accepted relay state."""
function magnitude_relay_candidate(
    settings::MagnitudeRelaySettings,
    accepted_state::MagnitudeRelayState,
    measurement::ProtectionMeasurement,
)
    _validate_magnitude_measurement(settings, measurement)
    measurement.release_tick > accepted_state.last_release_tick || throw(ArgumentError(
        "magnitude relay measurements must advance in strict release-tick order",
    ))
    state = deepcopy(accepted_state)
    elapsed_s = accepted_state.last_release_tick < 0 ? 0.0 :
        (measurement.release_tick - accepted_state.last_release_tick) * measurement.tick_s
    state.last_source_tick = measurement.source_tick
    state.last_release_tick = measurement.release_tick
    state.evaluation_count += 1
    if measurement.value === nothing
        state.blocked_count += 1
        signature = _magnitude_relay_decision_signature(
            settings,
            measurement,
            nothing,
            nothing,
            state,
            measurement.unavailable_reason,
        )
        return MagnitudeRelayCandidate(
            state,
            MagnitudeRelayDecision(
                settings.id,
                measurement.source_tick,
                measurement.release_tick,
                nothing,
                nothing,
                state.pickup_active,
                state.operated,
                state.timer_fraction,
                measurement.unavailable_reason,
                signature,
            ),
        )
    end
    value = _magnitude_relay_value(settings, measurement)
    pickup_active = _magnitude_pickup_state(
        settings,
        accepted_state.pickup_active,
        value,
    )
    pickup_multiple = _magnitude_pickup_multiple(settings, value)
    timer_fraction = _advance_magnitude_timer(
        settings,
        accepted_state.timer_fraction,
        pickup_active,
        pickup_multiple,
        elapsed_s,
    )
    operated = pickup_active && timer_fraction >= 1.0
    !accepted_state.operated && operated && (state.operation_count += 1)
    state.pickup_active = pickup_active
    state.operated = operated
    state.timer_fraction = timer_fraction
    signature = _magnitude_relay_decision_signature(
        settings,
        measurement,
        value,
        pickup_multiple,
        state,
        nothing,
    )
    return MagnitudeRelayCandidate(
        state,
        MagnitudeRelayDecision(
            settings.id,
            measurement.source_tick,
            measurement.release_tick,
            value,
            pickup_multiple,
            pickup_active,
            operated,
            timer_fraction,
            nothing,
            signature,
        ),
    )
end

function evaluate_magnitude_relay!(
    state::MagnitudeRelayState,
    settings::MagnitudeRelaySettings,
    measurement::ProtectionMeasurement,
)
    candidate = magnitude_relay_candidate(settings, state, measurement)
    accepted = candidate.state
    state.pickup_active = accepted.pickup_active
    state.operated = accepted.operated
    state.timer_fraction = accepted.timer_fraction
    state.last_source_tick = accepted.last_source_tick
    state.last_release_tick = accepted.last_release_tick
    state.evaluation_count = accepted.evaluation_count
    state.operation_count = accepted.operation_count
    state.blocked_count = accepted.blocked_count
    return candidate.decision
end

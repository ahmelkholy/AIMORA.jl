export DisruptiveEffectInsulator,
       LeaderProgressionInsulator,
       disruptive_effect_state,
       leader_progression_state,
       insulation_flashover_margin,
       apply_insulation_flashover!

"""Polarity-specific disruptive-effect insulator with a physical flashover conductance."""
mutable struct DisruptiveEffectInsulator <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    positive_threshold_voltage_v::Float64
    negative_threshold_voltage_v::Float64
    positive_exponent::Float64
    negative_exponent::Float64
    positive_critical_effect::Float64
    negative_critical_effect::Float64
    leakage_conductance_s::Float64
    flashover_conductance_s::Float64
    positive_effect::Float64
    negative_effect::Float64
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    flashed::Bool
    flashover_time_s::Float64
    flashover_polarity::Int
    flashover_count::Int
    provenance::SurgeParameterProvenance
end

function DisruptiveEffectInsulator(
    positive_node::Integer,
    negative_node::Integer;
    positive_threshold_voltage_v::Real,
    negative_threshold_voltage_v::Real,
    positive_exponent::Real,
    negative_exponent::Real,
    positive_critical_effect::Real,
    negative_critical_effect::Real,
    leakage_conductance_s::Real=1.0e-12,
    flashover_conductance_s::Real=1.0e4,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "volt, volt-to-exponent second, and siemens",
        "polarity-specific disruptive-effect insulation model",
    ),
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    positive_threshold = _nonnegative_finite(
        positive_threshold_voltage_v,
        "positive disruptive-effect threshold",
    )
    negative_threshold = _nonnegative_finite(
        negative_threshold_voltage_v,
        "negative disruptive-effect threshold",
    )
    positive_order = _positive_finite(positive_exponent, "positive disruptive-effect exponent")
    negative_order = _positive_finite(negative_exponent, "negative disruptive-effect exponent")
    positive_critical = _positive_finite(
        positive_critical_effect,
        "positive critical disruptive effect",
    )
    negative_critical = _positive_finite(
        negative_critical_effect,
        "negative critical disruptive effect",
    )
    leakage = _nonnegative_finite(leakage_conductance_s, "insulator leakage conductance")
    flashover = _positive_finite(flashover_conductance_s, "flashover conductance")
    flashover > leakage || throw(ArgumentError(
        "flashover conductance must exceed leakage conductance",
    ))
    _require_physical_provenance(provenance, "disruptive-effect insulator")
    return DisruptiveEffectInsulator(
        positive,
        negative,
        positive_threshold,
        negative_threshold,
        positive_order,
        negative_order,
        positive_critical,
        negative_critical,
        leakage,
        flashover,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        false,
        Inf,
        0,
        0,
        provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::DisruptiveEffectInsulator) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::DisruptiveEffectInsulator) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::DisruptiveEffectInsulator) =
    device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::DisruptiveEffectInsulator,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("insulator trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("insulator trial step must be positive"))
    _surge_companion_method(companion_method, "disruptive-effect insulator")
    device.trial_step_s = step_s
    return nothing
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::DisruptiveEffectInsulator,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "insulator voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("insulator evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    conductance = device.flashed ? device.flashover_conductance_s :
        device.leakage_conductance_s
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        conductance * voltage,
        conductance,
    )
end

function _disruptive_effect_increment(device::DisruptiveEffectInsulator, voltage_v::Float64)
    if voltage_v >= 0.0
        stress = max(voltage_v - device.positive_threshold_voltage_v, 0.0)
        return stress^device.positive_exponent, 0.0
    end
    stress = max(-voltage_v - device.negative_threshold_voltage_v, 0.0)
    return 0.0, stress^device.negative_exponent
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::DisruptiveEffectInsulator,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    conductance = device.flashed ? device.flashover_conductance_s :
        device.leakage_conductance_s
    expected = conductance * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted insulator current does not match trial"))
    elapsed = time_s - device.accepted_time_s
    elapsed >= 0.0 || throw(ArgumentError("insulator accepted time must be nondecreasing"))
    if !device.flashed
        old_positive, old_negative = _disruptive_effect_increment(
            device,
            device.accepted_voltage_v,
        )
        new_positive, new_negative = _disruptive_effect_increment(device, voltage)
        device.positive_effect += 0.5 * elapsed * (old_positive + new_positive)
        device.negative_effect += 0.5 * elapsed * (old_negative + new_negative)
    end
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    return nothing
end

function _apply_disruptive_effect_flashover!(
    device::DisruptiveEffectInsulator,
    time_s::Float64,
)
    apply_insulation_flashover!(device, time_s)
end

function NonlinearNetwork.nonlinear_device_event_surfaces(device::DisruptiveEffectInsulator)
    device.flashed && return ()
    return (
        NonlinearNetwork.NonlinearDeviceEventSurface(
            :positive_disruptive_effect_flashover,
            insulator -> insulator.positive_effect - insulator.positive_critical_effect,
            _apply_disruptive_effect_flashover!;
            direction=:rising,
            priority=10,
            topology_invalidating=true,
        ),
        NonlinearNetwork.NonlinearDeviceEventSurface(
            :negative_disruptive_effect_flashover,
            insulator -> insulator.negative_effect - insulator.negative_critical_effect,
            _apply_disruptive_effect_flashover!;
            direction=:rising,
            priority=10,
            topology_invalidating=true,
        ),
    )
end

function apply_insulation_flashover!(device::DisruptiveEffectInsulator, time_s::Real)
    time = _nonnegative_finite(time_s, "insulation flashover time")
    if !device.flashed
        positive_margin = device.positive_effect / device.positive_critical_effect
        negative_margin = device.negative_effect / device.negative_critical_effect
        max(positive_margin, negative_margin) >= 1.0 || throw(ArgumentError(
            "disruptive-effect flashover surface has not been reached",
        ))
        device.flashed = true
        device.flashover_time_s = time
        device.flashover_polarity = positive_margin >= negative_margin ? 1 : -1
        device.flashover_count += 1
    end
    return device
end

disruptive_effect_state(device::DisruptiveEffectInsulator) = (
    positive=device.positive_effect,
    negative=device.negative_effect,
    flashed=device.flashed,
    flashover_time_s=device.flashover_time_s,
    polarity=device.flashover_polarity,
)

function insulation_flashover_margin(device::DisruptiveEffectInsulator)
    return 1.0 - max(
        device.positive_effect / device.positive_critical_effect,
        device.negative_effect / device.negative_critical_effect,
    )
end

"""Polarity-aware engineering leader-progression insulator."""
mutable struct LeaderProgressionInsulator <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    gap_length_m::Float64
    positive_inception_field_v_per_m::Float64
    negative_inception_field_v_per_m::Float64
    positive_velocity_coefficient::Float64
    negative_velocity_coefficient::Float64
    velocity_exponent::Float64
    leakage_conductance_s::Float64
    flashover_conductance_s::Float64
    leader_length_m::Float64
    leader_polarity::Int
    incepted::Bool
    flashed::Bool
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    flashover_time_s::Float64
    provenance::SurgeParameterProvenance
end

function LeaderProgressionInsulator(
    positive_node::Integer,
    negative_node::Integer;
    gap_length_m::Real,
    positive_inception_field_v_per_m::Real,
    negative_inception_field_v_per_m::Real,
    positive_velocity_coefficient::Real,
    negative_velocity_coefficient::Real,
    velocity_exponent::Real,
    leakage_conductance_s::Real=1.0e-12,
    flashover_conductance_s::Real=1.0e4,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "metre, volt per metre, model-specific velocity coefficient, and siemens",
        "bounded polarity-specific engineering leader-progression model",
    ),
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    gap = _positive_finite(gap_length_m, "leader gap length")
    positive_field = _positive_finite(
        positive_inception_field_v_per_m,
        "positive leader inception field",
    )
    negative_field = _positive_finite(
        negative_inception_field_v_per_m,
        "negative leader inception field",
    )
    positive_velocity = _positive_finite(
        positive_velocity_coefficient,
        "positive leader velocity coefficient",
    )
    negative_velocity = _positive_finite(
        negative_velocity_coefficient,
        "negative leader velocity coefficient",
    )
    exponent = _positive_finite(velocity_exponent, "leader velocity exponent")
    leakage = _nonnegative_finite(leakage_conductance_s, "leader insulator leakage")
    flashover = _positive_finite(flashover_conductance_s, "leader flashover conductance")
    flashover > leakage || throw(ArgumentError(
        "leader flashover conductance must exceed leakage",
    ))
    _require_physical_provenance(provenance, "leader-progression insulator")
    return LeaderProgressionInsulator(
        positive, negative, gap, positive_field, negative_field,
        positive_velocity, negative_velocity, exponent, leakage, flashover,
        0.0, 0, false, false, 0.0, 0.0, 0.0, 0.0, Inf, provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::LeaderProgressionInsulator) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::LeaderProgressionInsulator) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::LeaderProgressionInsulator) =
    device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::LeaderProgressionInsulator,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("leader trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("leader trial step must be positive"))
    _surge_companion_method(companion_method, "leader-progression insulator")
    device.trial_step_s = step_s
    return nothing
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::LeaderProgressionInsulator,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "leader voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("leader evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    conductance = device.flashed ? device.flashover_conductance_s :
        device.leakage_conductance_s
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        conductance * voltage,
        conductance,
    )
end

function _leader_velocity_m_per_s(device::LeaderProgressionInsulator, voltage_v::Float64)
    remaining_gap = max(device.gap_length_m - device.leader_length_m, eps(Float64))
    field = abs(voltage_v) / remaining_gap
    polarity = voltage_v >= 0.0 ? 1 : -1
    inception = polarity > 0 ? device.positive_inception_field_v_per_m :
        device.negative_inception_field_v_per_m
    coefficient = polarity > 0 ? device.positive_velocity_coefficient :
        device.negative_velocity_coefficient
    velocity = coefficient * max(field - inception, 0.0)^device.velocity_exponent
    return velocity, polarity, field, inception
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::LeaderProgressionInsulator,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    conductance = device.flashed ? device.flashover_conductance_s :
        device.leakage_conductance_s
    expected = conductance * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted leader insulator current does not match trial"))
    elapsed = time_s - device.accepted_time_s
    elapsed >= 0.0 || throw(ArgumentError("leader accepted time must be nondecreasing"))
    if !device.flashed
        old_velocity, old_polarity, old_field, old_inception = _leader_velocity_m_per_s(
            device,
            device.accepted_voltage_v,
        )
        new_velocity, new_polarity, new_field, new_inception =
            _leader_velocity_m_per_s(device, voltage)
        if max(old_field - old_inception, new_field - new_inception) > 0.0
            device.incepted = true
            device.leader_polarity = new_velocity >= old_velocity ? new_polarity : old_polarity
            device.leader_length_m = min(
                device.gap_length_m,
                device.leader_length_m + 0.5 * elapsed * (old_velocity + new_velocity),
            )
        end
    end
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    return nothing
end

function _apply_leader_flashover!(device::LeaderProgressionInsulator, time_s::Float64)
    apply_insulation_flashover!(device, time_s)
end

function NonlinearNetwork.nonlinear_device_event_surfaces(device::LeaderProgressionInsulator)
    device.flashed && return ()
    return (
        NonlinearNetwork.NonlinearDeviceEventSurface(
            :leader_terminal_flashover,
            insulator -> insulator.leader_length_m - insulator.gap_length_m,
            _apply_leader_flashover!;
            direction=:rising,
            priority=10,
            topology_invalidating=true,
        ),
    )
end

function apply_insulation_flashover!(device::LeaderProgressionInsulator, time_s::Real)
    time = _nonnegative_finite(time_s, "leader flashover time")
    device.leader_length_m >= device.gap_length_m || throw(ArgumentError(
        "leader has not reached the terminal flashover surface",
    ))
    if !device.flashed
        device.flashed = true
        device.flashover_time_s = time
    end
    return device
end

leader_progression_state(device::LeaderProgressionInsulator) = (
    length_m=device.leader_length_m,
    remaining_gap_m=max(device.gap_length_m - device.leader_length_m, 0.0),
    polarity=device.leader_polarity,
    incepted=device.incepted,
    flashed=device.flashed,
    flashover_time_s=device.flashover_time_s,
)

insulation_flashover_margin(device::LeaderProgressionInsulator) =
    1.0 - device.leader_length_m / device.gap_length_m

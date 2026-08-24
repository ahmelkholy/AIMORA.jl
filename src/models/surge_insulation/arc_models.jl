export CombinedArcBranch,
       FaultArcBranch,
       VacuumInterruptionState,
       VacuumInterruptionBranch,
       arc_conductance_s,
       arc_dissipated_energy_j,
       arc_last_transition_time_s,
       ignite_arc!,
       request_vacuum_open!,
       vacuum_interruption_event_values,
       apply_vacuum_chop!,
       apply_vacuum_restrike!

"""Stateful two-terminal combined Cassie-Mayr arc with accepted-state mutation only."""
mutable struct CombinedArcBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    cassie_power_w::Float64
    mayr_power_w::Float64
    cassie_time_constant_s::Float64
    mayr_time_constant_s::Float64
    transition_power_w::Float64
    minimum_conductance_s::Float64
    maximum_conductance_s::Float64
    extinction_current_a::Float64
    accepted_conductance_s::Float64
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    dissipated_energy_j::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    mode::Symbol
    ignition_count::Int
    extinction_count::Int
    last_transition_time_s::Float64
    provenance::SurgeParameterProvenance
end

function CombinedArcBranch(
    positive_node::Integer,
    negative_node::Integer;
    cassie_power_w::Real,
    mayr_power_w::Real,
    cassie_time_constant_s::Real,
    mayr_time_constant_s::Real,
    transition_power_w::Real,
    initial_conductance_s::Real,
    minimum_conductance_s::Real=1.0e-9,
    maximum_conductance_s::Real=1.0e6,
    extinction_current_a::Real=1.0e-3,
    initially_ignited::Bool=true,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "watt, second, siemens, and ampere",
        "bounded positive combined Cassie-Mayr conductance arc",
    ),
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    cassie_power = _positive_finite(cassie_power_w, "Cassie arc power")
    mayr_power = _positive_finite(mayr_power_w, "Mayr arc power")
    cassie_time = _positive_finite(cassie_time_constant_s, "Cassie time constant")
    mayr_time = _positive_finite(mayr_time_constant_s, "Mayr time constant")
    transition_power = _positive_finite(transition_power_w, "arc transition power")
    minimum_conductance = _positive_finite(minimum_conductance_s, "minimum arc conductance")
    maximum_conductance = _positive_finite(maximum_conductance_s, "maximum arc conductance")
    maximum_conductance > minimum_conductance || throw(ArgumentError(
        "maximum arc conductance must exceed minimum conductance",
    ))
    initial_conductance = _positive_finite(initial_conductance_s, "initial arc conductance")
    minimum_conductance <= initial_conductance <= maximum_conductance ||
        throw(ArgumentError("initial arc conductance is outside declared bounds"))
    extinction_current = _nonnegative_finite(extinction_current_a, "arc extinction current")
    _require_physical_provenance(provenance, "combined arc")
    return CombinedArcBranch(
        positive,
        negative,
        cassie_power,
        mayr_power,
        cassie_time,
        mayr_time,
        transition_power,
        minimum_conductance,
        maximum_conductance,
        extinction_current,
        initially_ignited ? initial_conductance : minimum_conductance,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        initially_ignited ? :conducting : :extinguished,
        initially_ignited ? 1 : 0,
        0,
        0.0,
        provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::CombinedArcBranch) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::CombinedArcBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::CombinedArcBranch) = device.provenance

function _combined_arc_trial(device::CombinedArcBranch, voltage_v::Float64)
    device.mode === :conducting || return (
        conductance_s=device.minimum_conductance_s,
        derivative_s_per_v=0.0,
    )
    conductance = device.accepted_conductance_s
    step = device.trial_step_s
    step > 0.0 || return (conductance_s=conductance, derivative_s_per_v=0.0)
    power = conductance * voltage_v * voltage_v
    transition = device.transition_power_w
    weight = power / (power + transition)
    cassie_rate = (power / device.cassie_power_w - 1.0) /
        device.cassie_time_constant_s
    mayr_rate = (power / device.mayr_power_w - 1.0) /
        device.mayr_time_constant_s
    rate = weight * cassie_rate + (1.0 - weight) * mayr_rate
    raw_conductance = conductance * exp(clamp(step * rate, -80.0, 80.0))
    candidate = clamp(
        raw_conductance,
        device.minimum_conductance_s,
        device.maximum_conductance_s,
    )
    if candidate != raw_conductance
        return (conductance_s=candidate, derivative_s_per_v=0.0)
    end
    weight_derivative = transition / (power + transition)^2
    rate_derivative_per_power =
        weight_derivative * (cassie_rate - mayr_rate) +
        weight / (device.cassie_power_w * device.cassie_time_constant_s) +
        (1.0 - weight) / (device.mayr_power_w * device.mayr_time_constant_s)
    derivative = candidate * step * rate_derivative_per_power *
        2.0 * conductance * voltage_v
    return (conductance_s=candidate, derivative_s_per_v=derivative)
end

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::CombinedArcBranch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("arc trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("arc trial step must be positive"))
    _surge_companion_method(companion_method, "combined arc")
    device.trial_step_s = step_s
    return nothing
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::CombinedArcBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 ||
        throw(DimensionMismatch("arc voltage workspace must contain two terminals"))
    isfinite(time_s) || throw(ArgumentError("arc evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    trial = _combined_arc_trial(device, voltage)
    current = trial.conductance_s * voltage
    differential = trial.conductance_s + voltage * trial.derivative_s_per_v
    isfinite(current) && isfinite(differential) && differential >= 0.0 ||
        throw(ArgumentError("combined arc trial is nonfinite or nonpassive"))
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        current,
        differential,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::CombinedArcBranch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 && length(terminal_current_a) >= 2 ||
        throw(DimensionMismatch("arc acceptance requires two terminal values"))
    time_s >= device.accepted_time_s && isfinite(time_s) ||
        throw(ArgumentError("arc accepted time must be finite and nondecreasing"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    trial = _combined_arc_trial(device, voltage)
    expected = trial.conductance_s * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted arc current does not match the converged trial"))
    elapsed = time_s - device.accepted_time_s
    old_power = device.accepted_voltage_v * device.accepted_current_a
    new_power = voltage * current
    min(old_power, new_power) >= -256.0 * eps(Float64) * max(abs(old_power), abs(new_power), 1.0) ||
        throw(ArgumentError("combined arc accepted negative dissipative power"))
    device.dissipated_energy_j += 0.5 * elapsed * max(old_power + new_power, 0.0)
    device.accepted_conductance_s = trial.conductance_s
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    return nothing
end

function _extinguish_combined_arc!(device::CombinedArcBranch, time_s::Float64)
    if device.mode === :conducting
        device.mode = :extinguished
        device.accepted_conductance_s = device.minimum_conductance_s
        device.extinction_count += 1
        device.last_transition_time_s = time_s
    end
    return device
end

function NonlinearNetwork.nonlinear_device_event_surfaces(device::CombinedArcBranch)
    device.mode === :conducting || return ()
    return (
        NonlinearNetwork.NonlinearDeviceEventSurface(
            :arc_current_extinction,
            arc -> abs(arc.accepted_current_a) - arc.extinction_current_a,
            _extinguish_combined_arc!;
            direction=:falling,
            priority=30,
            topology_invalidating=true,
        ),
    )
end

function ignite_arc!(
    device::CombinedArcBranch,
    conductance_s::Real;
    time_s::Real=device.accepted_time_s,
)
    conductance = _positive_finite(conductance_s, "ignition conductance")
    device.minimum_conductance_s <= conductance <= device.maximum_conductance_s ||
        throw(ArgumentError("ignition conductance is outside declared bounds"))
    device.mode = :conducting
    device.accepted_conductance_s = conductance
    device.ignition_count += 1
    device.last_transition_time_s = _nonnegative_finite(time_s, "arc ignition time")
    return device
end

arc_conductance_s(device::CombinedArcBranch) = device.accepted_conductance_s
arc_dissipated_energy_j(device::CombinedArcBranch) = device.dissipated_energy_j
arc_last_transition_time_s(device::CombinedArcBranch) = device.last_transition_time_s

"""Dissipative Mayr-family fault arc with an explicit identity and length boundary."""
mutable struct FaultArcBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    cooling_power_w::Float64
    time_constant_s::Float64
    arc_length_m::Float64
    minimum_conductance_s::Float64
    maximum_conductance_s::Float64
    accepted_conductance_s::Float64
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    dissipated_energy_j::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    provenance::SurgeParameterProvenance
end

function FaultArcBranch(
    positive_node::Integer,
    negative_node::Integer;
    cooling_power_w::Real,
    time_constant_s::Real,
    arc_length_m::Real,
    initial_conductance_s::Real,
    minimum_conductance_s::Real=1.0e-9,
    maximum_conductance_s::Real=1.0e5,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "watt, second, metre, and siemens",
        "bounded positive conductance fault arc",
    ),
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    cooling = _positive_finite(cooling_power_w, "fault arc cooling power")
    time_constant = _positive_finite(time_constant_s, "fault arc time constant")
    length_m = _positive_finite(arc_length_m, "fault arc length")
    minimum = _positive_finite(minimum_conductance_s, "fault arc minimum conductance")
    maximum = _positive_finite(maximum_conductance_s, "fault arc maximum conductance")
    maximum > minimum || throw(ArgumentError("fault arc conductance bounds are invalid"))
    initial = _positive_finite(initial_conductance_s, "fault arc initial conductance")
    minimum <= initial <= maximum || throw(ArgumentError(
        "fault arc initial conductance is outside bounds",
    ))
    _require_physical_provenance(provenance, "fault arc")
    return FaultArcBranch(
        positive, negative, cooling, time_constant, length_m, minimum, maximum,
        initial, 0.0, 0.0, 0.0, 0.0, 0.0, provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::FaultArcBranch) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::FaultArcBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::FaultArcBranch) = device.provenance

function _fault_arc_trial(device::FaultArcBranch, voltage_v::Float64)
    conductance = device.accepted_conductance_s
    step = device.trial_step_s
    power = conductance * voltage_v * voltage_v
    rate = (power / device.cooling_power_w - 1.0) / device.time_constant_s
    raw = conductance * exp(clamp(step * rate, -80.0, 80.0))
    candidate = clamp(raw, device.minimum_conductance_s, device.maximum_conductance_s)
    derivative = candidate == raw ?
        candidate * step * 2.0 * conductance * voltage_v /
        (device.cooling_power_w * device.time_constant_s) : 0.0
    return candidate, derivative
end

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::FaultArcBranch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("fault arc trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("fault arc step must be positive"))
    _surge_companion_method(companion_method, "fault arc")
    device.trial_step_s = step_s
    return nothing
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::FaultArcBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "fault arc voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("fault arc evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    conductance, derivative = _fault_arc_trial(device, voltage)
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        conductance * voltage,
        conductance + voltage * derivative,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::FaultArcBranch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    conductance, _ = _fault_arc_trial(device, voltage)
    expected = conductance * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted fault arc current does not match trial"))
    elapsed = time_s - device.accepted_time_s
    elapsed >= 0.0 || throw(ArgumentError("fault arc time must be nondecreasing"))
    device.dissipated_energy_j += 0.5 * elapsed * max(
        device.accepted_voltage_v * device.accepted_current_a + voltage * current,
        0.0,
    )
    device.accepted_conductance_s = conductance
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    return nothing
end

arc_conductance_s(device::FaultArcBranch) = device.accepted_conductance_s
arc_dissipated_energy_j(device::FaultArcBranch) = device.dissipated_energy_j

"""Explicit vacuum interruption calendar and dielectric-recovery state."""
mutable struct VacuumInterruptionState
    identity::Symbol
    chopping_current_a::Float64
    initial_dielectric_strength_v::Float64
    dielectric_recovery_rate_v_per_s::Float64
    maximum_dielectric_strength_v::Float64
    separation_time_s::Float64
    last_transition_time_s::Float64
    last_chop_time_s::Float64
    last_restrike_time_s::Float64
    state::Symbol
    chop_count::Int
    restrike_count::Int
    provenance::SurgeParameterProvenance
end

function VacuumInterruptionState(
    identity::Symbol;
    chopping_current_a::Real,
    initial_dielectric_strength_v::Real,
    dielectric_recovery_rate_v_per_s::Real,
    maximum_dielectric_strength_v::Real,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "ampere, volt, volt per second",
        "bounded vacuum chopping and linear dielectric recovery model",
    ),
)
    identity == Symbol("") && throw(ArgumentError("vacuum interruption identity must not be empty"))
    chopping = _positive_finite(chopping_current_a, "vacuum chopping current")
    initial_strength = _nonnegative_finite(
        initial_dielectric_strength_v,
        "initial dielectric strength",
    )
    recovery_rate = _positive_finite(
        dielectric_recovery_rate_v_per_s,
        "dielectric recovery rate",
    )
    maximum_strength = _positive_finite(
        maximum_dielectric_strength_v,
        "maximum dielectric strength",
    )
    maximum_strength >= initial_strength || throw(ArgumentError(
        "maximum dielectric strength must not be below its initial value",
    ))
    _require_physical_provenance(provenance, "vacuum interruption")
    return VacuumInterruptionState(
        identity, chopping, initial_strength, recovery_rate, maximum_strength,
        Inf, 0.0, -1.0, -1.0, :closed, 0, 0, provenance,
    )
end

function request_vacuum_open!(state::VacuumInterruptionState, time_s::Real)
    time = _nonnegative_finite(time_s, "vacuum separation time")
    state.state === :closed || throw(ArgumentError(
        "vacuum opening can be requested only from the closed state",
    ))
    state.state = :separating
    state.separation_time_s = time
    state.last_transition_time_s = time
    return state
end

function _vacuum_dielectric_strength_v(state::VacuumInterruptionState, time_s::Float64)
    isfinite(state.separation_time_s) || return 0.0
    elapsed = max(time_s - state.separation_time_s, 0.0)
    return min(
        state.initial_dielectric_strength_v +
        state.dielectric_recovery_rate_v_per_s * elapsed,
        state.maximum_dielectric_strength_v,
    )
end

function vacuum_interruption_event_values(
    state::VacuumInterruptionState,
    current_a::Real,
    recovery_voltage_v::Real,
    time_s::Real,
)
    current = _finite_value(current_a, "vacuum current")
    voltage = _finite_value(recovery_voltage_v, "vacuum recovery voltage")
    time = _nonnegative_finite(time_s, "vacuum event time")
    chop_surface = state.state === :separating ?
        abs(current) - state.chopping_current_a : Inf
    dielectric_strength = _vacuum_dielectric_strength_v(state, time)
    restrike_surface = state.state === :open ?
        abs(voltage) - dielectric_strength : -Inf
    return (
        chop_surface_a=chop_surface,
        restrike_surface_v=restrike_surface,
        dielectric_strength_v=dielectric_strength,
    )
end

function apply_vacuum_chop!(state::VacuumInterruptionState, time_s::Real)
    time = _nonnegative_finite(time_s, "vacuum chop time")
    state.state === :separating || throw(ArgumentError(
        "vacuum chopping requires a separating contact",
    ))
    time >= state.separation_time_s || throw(ArgumentError(
        "vacuum chop cannot precede contact separation",
    ))
    state.state = :open
    state.last_transition_time_s = time
    state.last_chop_time_s = time
    state.chop_count += 1
    return state
end

function apply_vacuum_restrike!(state::VacuumInterruptionState, time_s::Real)
    time = _nonnegative_finite(time_s, "vacuum restrike time")
    state.state === :open || throw(ArgumentError("vacuum restrike requires an open gap"))
    time >= state.last_transition_time_s || throw(ArgumentError(
        "vacuum restrike time must be nondecreasing",
    ))
    state.state = :restruck
    state.last_transition_time_s = time
    state.last_restrike_time_s = time
    state.restrike_count += 1
    return state
end

"""Two-terminal vacuum contact with solver-localized chop and restrike surfaces.

The metallic contact carries current only while closed. After separation, the
associated combined arc owns the material current path. Chopping extinguishes
that exact arc owner and dielectric failure reignites it at the localized time.
"""
mutable struct VacuumInterruptionBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    state::VacuumInterruptionState
    arc_owner::CombinedArcBranch
    closed_conductance_s::Float64
    open_conductance_s::Float64
    reignition_conductance_s::Float64
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    provenance::SurgeParameterProvenance
end

function VacuumInterruptionBranch(
    positive_node::Integer,
    negative_node::Integer,
    state::VacuumInterruptionState,
    arc_owner::CombinedArcBranch;
    closed_conductance_s::Real=1.0e6,
    open_conductance_s::Real=1.0e-9,
    reignition_conductance_s::Real=1.0,
    provenance::SurgeParameterProvenance=state.provenance,
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    (positive, negative) == (arc_owner.positive_node, arc_owner.negative_node) ||
        throw(ArgumentError("vacuum contact and arc terminals must be identical"))
    closed = _positive_finite(closed_conductance_s, "closed vacuum-contact conductance")
    open = _nonnegative_finite(open_conductance_s, "open vacuum-contact conductance")
    closed > open || throw(ArgumentError(
        "closed vacuum-contact conductance must exceed open conductance",
    ))
    reignition = _positive_finite(reignition_conductance_s, "vacuum reignition conductance")
    arc_owner.minimum_conductance_s <= reignition <= arc_owner.maximum_conductance_s ||
        throw(ArgumentError("vacuum reignition conductance is outside arc bounds"))
    _require_physical_provenance(provenance, "vacuum interruption branch")
    return VacuumInterruptionBranch(
        positive,
        negative,
        state,
        arc_owner,
        closed,
        open,
        reignition,
        0.0,
        0.0,
        0.0,
        0.0,
        provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::VacuumInterruptionBranch) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::VacuumInterruptionBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::VacuumInterruptionBranch) =
    device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::VacuumInterruptionBranch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("vacuum-contact trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError(
        "vacuum-contact trial step must be finite and positive",
    ))
    _surge_companion_method(companion_method, "vacuum interruption branch")
    device.trial_step_s = step_s
    return nothing
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::VacuumInterruptionBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "vacuum-contact voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("vacuum-contact evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    conductance = device.state.state === :closed ?
        device.closed_conductance_s : device.open_conductance_s
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        conductance * voltage,
        conductance,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::VacuumInterruptionBranch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 && length(terminal_current_a) >= 2 ||
        throw(DimensionMismatch("vacuum-contact acceptance requires two terminals"))
    time_s >= device.accepted_time_s && isfinite(time_s) || throw(ArgumentError(
        "vacuum-contact accepted time must be finite and nondecreasing",
    ))
    device.accepted_voltage_v = terminal_voltage_v[1] - terminal_voltage_v[2]
    device.accepted_current_a = terminal_current_a[1]
    device.accepted_time_s = time_s
    return nothing
end

function _apply_vacuum_branch_chop!(device::VacuumInterruptionBranch, time_s::Float64)
    apply_vacuum_chop!(device.state, time_s)
    _extinguish_combined_arc!(device.arc_owner, time_s)
    return device
end

function _apply_vacuum_branch_restrike!(device::VacuumInterruptionBranch, time_s::Float64)
    apply_vacuum_restrike!(device.state, time_s)
    ignite_arc!(device.arc_owner, device.reignition_conductance_s; time_s)
    return device
end

function NonlinearNetwork.nonlinear_device_event_surfaces(
    device::VacuumInterruptionBranch,
)
    if device.state.state === :separating
        return (
            NonlinearNetwork.NonlinearDeviceEventSurface(
                :vacuum_current_chop,
                branch -> abs(branch.arc_owner.accepted_current_a) -
                    branch.state.chopping_current_a,
                _apply_vacuum_branch_chop!;
                direction=:falling,
                priority=20,
                topology_invalidating=true,
            ),
        )
    elseif device.state.state === :open
        return (
            NonlinearNetwork.NonlinearDeviceEventSurface(
                :vacuum_dielectric_restrike,
                (branch, time_s) -> abs(branch.accepted_voltage_v) -
                    _vacuum_dielectric_strength_v(branch.state, time_s),
                _apply_vacuum_branch_restrike!;
                direction=:rising,
                priority=10,
                topology_invalidating=true,
            ),
        )
    end
    return ()
end

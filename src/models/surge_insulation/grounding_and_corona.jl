export PositiveRealGroundingModel,
       IonizingGroundBranch,
       DynamicCoronaBranch,
       grounding_admittance_s,
       grounding_impedance_ohm,
       ground_potential_rise_v,
       ionized_radius_m,
       grounding_dissipated_energy_j,
       corona_charge_c,
       corona_dissipated_energy_j

"""Passive rational grounding admittance `G0 + sum(r_k*s/(s+p_k))`."""
struct PositiveRealGroundingModel
    direct_conductance_s::Float64
    pole_rates_per_s::Vector{Float64}
    residue_conductances_s::Vector{Float64}
    minimum_frequency_hz::Float64
    maximum_frequency_hz::Float64
    provenance::SurgeParameterProvenance

    function PositiveRealGroundingModel(
        direct_conductance_s::Real,
        pole_rates_per_s::AbstractVector{<:Real},
        residue_conductances_s::AbstractVector{<:Real};
        minimum_frequency_hz::Real=0.0,
        maximum_frequency_hz::Real,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "siemens, inverse second, and hertz",
            "positive-real rational grounding admittance over a declared band",
        ),
    )
        direct = _positive_finite(direct_conductance_s, "ground direct conductance")
        poles = Float64.(pole_rates_per_s)
        residues = Float64.(residue_conductances_s)
        length(poles) == length(residues) || throw(ArgumentError(
            "ground pole and residue counts must match",
        ))
        all(value -> isfinite(value) && value > 0.0, poles) || throw(ArgumentError(
            "ground pole rates must be finite and positive",
        ))
        all(value -> isfinite(value) && value >= 0.0, residues) || throw(ArgumentError(
            "ground residues must be finite and nonnegative",
        ))
        minimum_frequency = _nonnegative_finite(minimum_frequency_hz, "minimum ground frequency")
        maximum_frequency = _positive_finite(maximum_frequency_hz, "maximum ground frequency")
        maximum_frequency > minimum_frequency || throw(ArgumentError(
            "maximum ground frequency must exceed minimum frequency",
        ))
        _require_physical_provenance(provenance, "grounding model")
        return new(
            direct,
            poles,
            residues,
            minimum_frequency,
            maximum_frequency,
            provenance,
        )
    end
end

function grounding_admittance_s(model::PositiveRealGroundingModel, frequency_hz::Real)
    frequency = _nonnegative_finite(frequency_hz, "grounding evaluation frequency")
    model.minimum_frequency_hz <= frequency <= model.maximum_frequency_hz ||
        throw(DomainError(frequency, "grounding frequency is outside the declared fit band"))
    laplace = complex(0.0, 2.0 * pi * frequency)
    admittance = complex(model.direct_conductance_s, 0.0)
    for (pole, residue) in zip(model.pole_rates_per_s, model.residue_conductances_s)
        admittance += residue * laplace / (laplace + pole)
    end
    real(admittance) >= -256.0 * eps(Float64) * max(abs(admittance), 1.0) ||
        throw(ArgumentError("grounding admittance violates positive-real passivity"))
    return admittance
end

grounding_impedance_ohm(model::PositiveRealGroundingModel, frequency_hz::Real) =
    inv(grounding_admittance_s(model, frequency_hz))

"""Compact dissipative soil-ionization branch with bounded effective radius state."""
mutable struct IonizingGroundBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    node::Int
    linear_conductance_s::Float64
    electrode_radius_m::Float64
    maximum_ionized_radius_m::Float64
    critical_field_v_per_m::Float64
    expansion_rate_m_per_v_s::Float64
    recovery_rate_per_s::Float64
    accepted_ionized_radius_m::Float64
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    dissipated_energy_j::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    provenance::SurgeParameterProvenance
end

function IonizingGroundBranch(
    node::Integer;
    linear_resistance_ohm::Real,
    electrode_radius_m::Real,
    maximum_ionized_radius_m::Real,
    critical_field_v_per_m::Real,
    expansion_rate_m_per_v_s::Real,
    recovery_rate_per_s::Real,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "ohm, metre, volt per metre, metre per volt-second, and inverse second",
        "compact positive dissipative effective-radius soil-ionization model",
    ),
)
    checked_node = Int(node)
    checked_node > 0 || throw(ArgumentError("ionizing ground node must be positive"))
    resistance = _positive_finite(linear_resistance_ohm, "linear grounding resistance")
    radius = _positive_finite(electrode_radius_m, "ground electrode radius")
    maximum_radius = _positive_finite(maximum_ionized_radius_m, "maximum ionized radius")
    maximum_radius >= radius || throw(ArgumentError(
        "maximum ionized radius must be at least the electrode radius",
    ))
    critical_field = _positive_finite(critical_field_v_per_m, "critical soil field")
    expansion_rate = _positive_finite(
        expansion_rate_m_per_v_s,
        "soil ionization expansion rate",
    )
    recovery_rate = _positive_finite(recovery_rate_per_s, "soil recovery rate")
    _require_physical_provenance(provenance, "ionizing ground")
    return IonizingGroundBranch(
        checked_node,
        inv(resistance),
        radius,
        maximum_radius,
        critical_field,
        expansion_rate,
        recovery_rate,
        radius,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::IonizingGroundBranch) = (device.node, 0)
NonlinearNetwork.nonlinear_device_formulation(::IonizingGroundBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::IonizingGroundBranch) = device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::IonizingGroundBranch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("ground trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("ground trial step must be positive"))
    _surge_companion_method(companion_method, "ionizing ground")
    device.trial_step_s = step_s
    return nothing
end

function _ionizing_ground_trial(device::IonizingGroundBranch, voltage_v::Float64)
    radius = device.accepted_ionized_radius_m
    field = abs(voltage_v) / max(radius, device.electrode_radius_m)
    step = device.trial_step_s
    if field > device.critical_field_v_per_m
        raw_radius = radius + step * device.expansion_rate_m_per_v_s *
            (field - device.critical_field_v_per_m)
        derivative_radius = step * device.expansion_rate_m_per_v_s *
            sign(voltage_v) / max(radius, device.electrode_radius_m)
    else
        raw_radius = radius - step * device.recovery_rate_per_s *
            (radius - device.electrode_radius_m)
        derivative_radius = 0.0
    end
    candidate_radius = clamp(
        raw_radius,
        device.electrode_radius_m,
        device.maximum_ionized_radius_m,
    )
    candidate_radius == raw_radius || (derivative_radius = 0.0)
    conductance = device.linear_conductance_s *
        candidate_radius / device.electrode_radius_m
    conductance_derivative = device.linear_conductance_s *
        derivative_radius / device.electrode_radius_m
    return candidate_radius, conductance, conductance_derivative
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::IonizingGroundBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "ionizing ground voltage workspace must contain node and ground",
    ))
    isfinite(time_s) || throw(ArgumentError("ground evaluation time must be finite"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    _, conductance, conductance_derivative = _ionizing_ground_trial(device, voltage)
    current = conductance * voltage
    differential = conductance + conductance_derivative * voltage
    differential >= 0.0 && isfinite(differential) || throw(ArgumentError(
        "ionizing ground trial is nonfinite or nonpassive",
    ))
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        current,
        differential,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::IonizingGroundBranch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    radius, conductance, _ = _ionizing_ground_trial(device, voltage)
    expected = conductance * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted ground current does not match trial"))
    elapsed = time_s - device.accepted_time_s
    elapsed >= 0.0 || throw(ArgumentError("ground accepted time must be nondecreasing"))
    device.dissipated_energy_j += 0.5 * elapsed * max(
        device.accepted_voltage_v * device.accepted_current_a + voltage * current,
        0.0,
    )
    device.accepted_ionized_radius_m = radius
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    return nothing
end

ground_potential_rise_v(device::IonizingGroundBranch) = device.accepted_voltage_v
ionized_radius_m(device::IonizingGroundBranch) = device.accepted_ionized_radius_m
grounding_dissipated_energy_j(device::IonizingGroundBranch) = device.dissipated_energy_j

"""Passive voltage-dependent corona charge companion with explicit hysteresis."""
mutable struct DynamicCoronaBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    base_capacitance_f::Float64
    incremental_capacitance_f_per_v::Float64
    onset_voltage_v::Float64
    extinction_voltage_v::Float64
    loss_conductance_s::Float64
    active::Bool
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    accepted_charge_c::Float64
    dissipated_energy_j::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    provenance::SurgeParameterProvenance
end

function DynamicCoronaBranch(
    positive_node::Integer,
    negative_node::Integer;
    base_capacitance_f::Real,
    incremental_capacitance_f_per_v::Real,
    onset_voltage_v::Real,
    extinction_voltage_v::Real,
    loss_conductance_s::Real,
    provenance::SurgeParameterProvenance=_surge_default_provenance(
        "farad, farad per volt, volt, and siemens",
        "passive bounded dynamic corona charge with explicit hysteresis",
    ),
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    base_capacitance = _positive_finite(base_capacitance_f, "corona base capacitance")
    incremental = _nonnegative_finite(
        incremental_capacitance_f_per_v,
        "corona incremental capacitance",
    )
    onset = _positive_finite(onset_voltage_v, "corona onset voltage")
    extinction = _nonnegative_finite(extinction_voltage_v, "corona extinction voltage")
    extinction <= onset || throw(ArgumentError(
        "corona extinction voltage must not exceed onset voltage",
    ))
    loss = _nonnegative_finite(loss_conductance_s, "corona loss conductance")
    _require_physical_provenance(provenance, "dynamic corona")
    return DynamicCoronaBranch(
        positive, negative, base_capacitance, incremental, onset, extinction, loss,
        false, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::DynamicCoronaBranch) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::DynamicCoronaBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::DynamicCoronaBranch) = device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::DynamicCoronaBranch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("corona trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("corona trial step must be positive"))
    _surge_companion_method(companion_method, "dynamic corona")
    device.trial_step_s = step_s
    return nothing
end

function _corona_charge_and_derivative(device::DynamicCoronaBranch, voltage_v::Float64)
    magnitude = abs(voltage_v)
    active = device.active ? magnitude > device.extinction_voltage_v :
        magnitude >= device.onset_voltage_v
    overvoltage = active ? max(magnitude - device.extinction_voltage_v, 0.0) : 0.0
    charge_magnitude = device.base_capacitance_f * magnitude +
        0.5 * device.incremental_capacitance_f_per_v * overvoltage^2
    differential_capacitance = device.base_capacitance_f +
        device.incremental_capacitance_f_per_v * overvoltage
    return copysign(charge_magnitude, voltage_v), differential_capacitance, active
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::DynamicCoronaBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "corona voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("corona evaluation time must be finite"))
    device.trial_step_s > 0.0 || throw(ArgumentError("corona step must be prepared"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    charge, differential_capacitance, _ = _corona_charge_and_derivative(device, voltage)
    current = (charge - device.accepted_charge_c) / device.trial_step_s +
        device.loss_conductance_s * voltage
    differential = differential_capacitance / device.trial_step_s +
        device.loss_conductance_s
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        current,
        differential,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::DynamicCoronaBranch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    charge, _, active = _corona_charge_and_derivative(device, voltage)
    expected = (charge - device.accepted_charge_c) / device.trial_step_s +
        device.loss_conductance_s * voltage
    abs(current - expected) <= 256.0 * eps(Float64) * max(abs(current), abs(expected), 1.0) ||
        throw(ArgumentError("accepted corona current does not match trial"))
    elapsed = time_s - device.accepted_time_s
    elapsed >= 0.0 || throw(ArgumentError("corona accepted time must be nondecreasing"))
    device.dissipated_energy_j += 0.5 * elapsed * device.loss_conductance_s *
        (device.accepted_voltage_v^2 + voltage^2)
    device.active = active
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_charge_c = charge
    device.accepted_time_s = time_s
    return nothing
end

corona_charge_c(device::DynamicCoronaBranch) = device.accepted_charge_c
corona_dissipated_energy_j(device::DynamicCoronaBranch) = device.dissipated_energy_j

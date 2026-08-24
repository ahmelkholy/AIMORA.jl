export MetalOxideCharacteristic,
       MetalOxideArrester,
       DynamicArresterEquivalent,
       metal_oxide_current_and_derivative,
       arrester_charge_c,
       arrester_absorbed_energy_j,
       arrester_temperature_k,
       arrester_duty_margin

"""Continuous monotone piecewise power-law metal-oxide current characteristic."""
struct MetalOxideCharacteristic
    current_knots_a::Vector{Float64}
    voltage_knots_v::Vector{Float64}
    exponents::Vector{Float64}
    coefficients_a_per_v_alpha::Vector{Float64}
    extrapolation::Symbol
    provenance::SurgeParameterProvenance

    function MetalOxideCharacteristic(
        current_knots_a::AbstractVector{<:Real},
        voltage_knots_v::AbstractVector{<:Real};
        extrapolation::Symbol=:refuse,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "ampere and volt",
            "positive monotone piecewise power-law metal-oxide characteristic",
        ),
    )
        currents = Float64.(current_knots_a)
        voltages = Float64.(voltage_knots_v)
        length(currents) == length(voltages) >= 2 || throw(ArgumentError(
            "metal-oxide characteristic requires at least two aligned knots",
        ))
        all(value -> isfinite(value) && value > 0.0, currents) &&
            all(value -> isfinite(value) && value > 0.0, voltages) ||
            throw(ArgumentError("metal-oxide knots must be finite and positive"))
        all(index -> currents[index] > currents[index - 1], 2:length(currents)) &&
            all(index -> voltages[index] > voltages[index - 1], 2:length(voltages)) ||
            throw(ArgumentError("metal-oxide current and voltage knots must increase strictly"))
        extrapolation in (:power_law, :refuse) || throw(ArgumentError(
            "metal-oxide extrapolation must be :power_law or :refuse",
        ))
        exponents = Float64[]
        coefficients = Float64[]
        for index in 1:(length(currents) - 1)
            exponent = log(currents[index + 1] / currents[index]) /
                log(voltages[index + 1] / voltages[index])
            isfinite(exponent) && exponent > 0.0 || throw(ArgumentError(
                "metal-oxide segment exponent must be finite and positive",
            ))
            coefficient = currents[index] / voltages[index]^exponent
            push!(exponents, exponent)
            push!(coefficients, coefficient)
        end
        _require_physical_provenance(provenance, "metal-oxide characteristic")
        return new(currents, voltages, exponents, coefficients, extrapolation, provenance)
    end
end

function metal_oxide_current_and_derivative(
    characteristic::MetalOxideCharacteristic,
    voltage_v::Real,
)
    voltage = _finite_value(voltage_v, "metal-oxide voltage")
    magnitude = abs(voltage)
    magnitude == 0.0 && return (current_a=0.0, derivative_s=0.0, segment=0)
    knots = characteristic.voltage_knots_v
    if magnitude < first(knots)
        characteristic.extrapolation === :power_law || throw(DomainError(
            voltage,
            "metal-oxide voltage is below the declared characteristic",
        ))
        segment = 1
    elseif magnitude > last(knots)
        characteristic.extrapolation === :power_law || throw(DomainError(
            voltage,
            "metal-oxide voltage is above the declared characteristic",
        ))
        segment = length(characteristic.exponents)
    else
        segment = clamp(searchsortedlast(knots, magnitude), 1, length(characteristic.exponents))
    end
    exponent = characteristic.exponents[segment]
    coefficient = characteristic.coefficients_a_per_v_alpha[segment]
    current_magnitude = coefficient * magnitude^exponent
    derivative = exponent * coefficient * magnitude^(exponent - 1.0)
    return (
        current_a=copysign(current_magnitude, voltage),
        derivative_s=derivative,
        segment=segment,
    )
end

"""Metal-oxide arrester with optional gap, charge, energy, and admitted thermal state."""
mutable struct MetalOxideArrester <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    characteristic::MetalOxideCharacteristic
    gap_sparkover_voltage_v::Float64
    leakage_conductance_s::Float64
    thermal_capacitance_j_per_k::Float64
    thermal_resistance_k_per_w::Float64
    ambient_temperature_k::Float64
    maximum_temperature_k::Float64
    maximum_energy_j::Float64
    conducting::Bool
    failed::Bool
    accepted_voltage_v::Float64
    accepted_current_a::Float64
    charge_c::Float64
    absorbed_energy_j::Float64
    temperature_k::Float64
    peak_current_a::Float64
    peak_residual_voltage_v::Float64
    accepted_time_s::Float64
    trial_step_s::Float64
    sparkover_count::Int
    provenance::SurgeParameterProvenance
end

function MetalOxideArrester(
    positive_node::Integer,
    negative_node::Integer,
    characteristic::MetalOxideCharacteristic;
    gap_sparkover_voltage_v::Real=0.0,
    leakage_conductance_s::Real=1.0e-12,
    thermal_capacitance_j_per_k::Real=1.0,
    thermal_resistance_k_per_w::Real=Inf,
    ambient_temperature_k::Real=293.15,
    maximum_temperature_k::Real=Inf,
    maximum_energy_j::Real=Inf,
    provenance::SurgeParameterProvenance=characteristic.provenance,
)
    positive, negative = _surge_nodes(positive_node, negative_node)
    gap = _nonnegative_finite(gap_sparkover_voltage_v, "arrester gap voltage")
    leakage = _nonnegative_finite(leakage_conductance_s, "arrester leakage conductance")
    thermal_capacitance = _positive_finite(
        thermal_capacitance_j_per_k,
        "arrester thermal capacitance",
    )
    thermal_resistance = Float64(thermal_resistance_k_per_w)
    (isinf(thermal_resistance) && thermal_resistance > 0.0) ||
        (isfinite(thermal_resistance) && thermal_resistance > 0.0) ||
        throw(ArgumentError("arrester thermal resistance must be positive or Inf"))
    ambient = _positive_finite(ambient_temperature_k, "arrester ambient temperature")
    maximum_temperature = Float64(maximum_temperature_k)
    (isinf(maximum_temperature) && maximum_temperature > 0.0) ||
        (isfinite(maximum_temperature) && maximum_temperature >= ambient) ||
        throw(ArgumentError("arrester maximum temperature must be at least ambient or Inf"))
    maximum_energy = Float64(maximum_energy_j)
    (isinf(maximum_energy) && maximum_energy > 0.0) ||
        (isfinite(maximum_energy) && maximum_energy > 0.0) ||
        throw(ArgumentError("arrester maximum energy must be positive or Inf"))
    _require_physical_provenance(provenance, "metal-oxide arrester")
    return MetalOxideArrester(
        positive,
        negative,
        characteristic,
        gap,
        leakage,
        thermal_capacitance,
        thermal_resistance,
        ambient,
        maximum_temperature,
        maximum_energy,
        gap == 0.0,
        false,
        0.0,
        0.0,
        0.0,
        0.0,
        ambient,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        provenance,
    )
end

NonlinearNetwork.nonlinear_terminal_nodes(device::MetalOxideArrester) =
    (device.positive_node, device.negative_node)
NonlinearNetwork.nonlinear_device_formulation(::MetalOxideArrester) =
    NonlinearNetwork.PhysicalConstitutiveCurrent
NonlinearNetwork.nonlinear_device_provenance(device::MetalOxideArrester) = device.provenance

function NonlinearNetwork.prepare_nonlinear_device_step!(
    device::MetalOxideArrester,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) || throw(ArgumentError("arrester trial time must be finite"))
    step_s > 0.0 && isfinite(step_s) || throw(ArgumentError("arrester trial step must be positive"))
    _surge_companion_method(companion_method, "metal-oxide arrester")
    device.failed && throw(ArgumentError("thermally or energetically failed arrester cannot be stamped"))
    device.trial_step_s = step_s
    return nothing
end

function _arrester_current_and_derivative(device::MetalOxideArrester, voltage_v::Float64)
    if device.conducting
        return metal_oxide_current_and_derivative(device.characteristic, voltage_v)
    end
    return (
        current_a=device.leakage_conductance_s * voltage_v,
        derivative_s=device.leakage_conductance_s,
        segment=0,
    )
end

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::MetalOxideArrester,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "arrester voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("arrester evaluation time must be finite"))
    device.failed && throw(ArgumentError("failed arrester has no admitted constitutive stamp"))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    evaluation = _arrester_current_and_derivative(device, voltage)
    voltage * evaluation.current_a >= -256.0 * eps(Float64) *
        max(abs(voltage * evaluation.current_a), 1.0) ||
        throw(ArgumentError("metal-oxide trial violates passive power orientation"))
    return _fill_two_terminal_stamp!(
        terminal_current_a,
        terminal_jacobian_s,
        evaluation.current_a,
        evaluation.derivative_s,
    )
end

function NonlinearNetwork.accept_nonlinear_device_state!(
    device::MetalOxideArrester,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_voltage_v) >= 2 && length(terminal_current_a) >= 2 ||
        throw(DimensionMismatch("arrester acceptance requires two terminal values"))
    time_s >= device.accepted_time_s && isfinite(time_s) || throw(ArgumentError(
        "arrester accepted time must be finite and nondecreasing",
    ))
    voltage = terminal_voltage_v[1] - terminal_voltage_v[2]
    current = terminal_current_a[1]
    evaluation = _arrester_current_and_derivative(device, voltage)
    abs(current - evaluation.current_a) <= 256.0 * eps(Float64) *
        max(abs(current), abs(evaluation.current_a), 1.0) ||
        throw(ArgumentError("accepted arrester current does not match converged trial"))
    elapsed = time_s - device.accepted_time_s
    old_power = device.accepted_voltage_v * device.accepted_current_a
    new_power = voltage * current
    average_power = 0.5 * max(old_power + new_power, 0.0)
    device.charge_c += 0.5 * elapsed * (device.accepted_current_a + current)
    device.absorbed_energy_j += elapsed * average_power
    cooling_power = isinf(device.thermal_resistance_k_per_w) ? 0.0 :
        (device.temperature_k - device.ambient_temperature_k) /
        device.thermal_resistance_k_per_w
    device.temperature_k += elapsed * (average_power - cooling_power) /
        device.thermal_capacitance_j_per_k
    device.temperature_k = max(device.temperature_k, device.ambient_temperature_k)
    device.peak_current_a = max(device.peak_current_a, abs(current))
    device.peak_residual_voltage_v = max(device.peak_residual_voltage_v, abs(voltage))
    device.accepted_voltage_v = voltage
    device.accepted_current_a = current
    device.accepted_time_s = time_s
    device.failed = device.absorbed_energy_j > device.maximum_energy_j ||
        device.temperature_k > device.maximum_temperature_k
    return nothing
end

function _sparkover_arrester!(device::MetalOxideArrester, _time_s::Float64)
    if !device.conducting
        device.conducting = true
        device.sparkover_count += 1
    end
    return device
end

function NonlinearNetwork.nonlinear_device_event_surfaces(device::MetalOxideArrester)
    device.gap_sparkover_voltage_v > 0.0 && !device.conducting && !device.failed || return ()
    return (
        NonlinearNetwork.NonlinearDeviceEventSurface(
            :arrester_gap_sparkover,
            arrester -> abs(arrester.accepted_voltage_v) - arrester.gap_sparkover_voltage_v,
            _sparkover_arrester!;
            direction=:rising,
            priority=20,
            topology_invalidating=true,
        ),
    )
end

arrester_charge_c(device::MetalOxideArrester) = device.charge_c
arrester_absorbed_energy_j(device::MetalOxideArrester) = device.absorbed_energy_j
arrester_temperature_k(device::MetalOxideArrester) = device.temperature_k

function arrester_duty_margin(device::MetalOxideArrester)
    energy_margin = isinf(device.maximum_energy_j) ? Inf :
        device.maximum_energy_j - device.absorbed_energy_j
    temperature_margin = isinf(device.maximum_temperature_k) ? Inf :
        device.maximum_temperature_k - device.temperature_k
    return (
        energy_margin_j=energy_margin,
        temperature_margin_k=temperature_margin,
        failed=device.failed,
    )
end

"""Explicit two-nonlinear-column dynamic arrester topology parameters."""
struct DynamicArresterEquivalent
    first_characteristic::MetalOxideCharacteristic
    second_characteristic::MetalOxideCharacteristic
    series_resistance_ohm::Float64
    series_inductance_h::Float64
    filter_resistance_ohm::Float64
    filter_inductance_h::Float64
    shunt_capacitance_f::Float64
    provenance::SurgeParameterProvenance

    function DynamicArresterEquivalent(
        first_characteristic::MetalOxideCharacteristic,
        second_characteristic::MetalOxideCharacteristic;
        series_resistance_ohm::Real,
        series_inductance_h::Real,
        filter_resistance_ohm::Real,
        filter_inductance_h::Real,
        shunt_capacitance_f::Real,
        provenance::SurgeParameterProvenance=first_characteristic.provenance,
    )
        series_resistance = _nonnegative_finite(series_resistance_ohm, "arrester series resistance")
        series_inductance = _nonnegative_finite(series_inductance_h, "arrester series inductance")
        filter_resistance = _positive_finite(filter_resistance_ohm, "arrester filter resistance")
        filter_inductance = _positive_finite(filter_inductance_h, "arrester filter inductance")
        shunt_capacitance = _nonnegative_finite(shunt_capacitance_f, "arrester shunt capacitance")
        _require_physical_provenance(provenance, "dynamic arrester equivalent")
        return new(
            first_characteristic,
            second_characteristic,
            series_resistance,
            series_inductance,
            filter_resistance,
            filter_inductance,
            shunt_capacitance,
            provenance,
        )
    end
end

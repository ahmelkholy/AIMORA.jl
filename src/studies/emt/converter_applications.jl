using ..ConverterSystems

const _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_OWNER =
    "emt.average_converter_application"
const _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_SCHEMA = Int64(2)

function _average_converter_application_snapshot_fail(code::Symbol, message)
    throw(PortableSnapshots.PortableSnapshotFailure(code, String(message)))
end

mutable struct AverageConverterApplicationRuntime{S}
    study::S
    state::Vector{Float64}
    time_s::Float64
    accepted_step_index::Int
    operating_mode::ConverterSystems.ConverterApplicationOperatingMode
    application_event_cursor::Int
    input_energy_j::Float64
    dissipated_energy_j::Float64
    time_trace_s::Vector{Float64}
    input_voltage_trace_v::Matrix{Float64}
    input_current_trace_a::Matrix{Float64}
    output_voltage_trace_v::Matrix{Float64}
    output_current_trace_a::Matrix{Float64}
    converter_current_trace_a::Matrix{Float64}
    dc_link_voltage_trace_v::Matrix{Float64}
    stage_power_trace_w::Matrix{Float64}
    stored_energy_trace_j::Vector{Float64}
    dissipated_power_trace_w::Vector{Float64}
    energy_residual_trace_w::Vector{Float64}
    control_error_trace::Matrix{Float64}
    operating_mode_trace::Vector{Int64}
    load_delivery_trace::Vector{Int64}
    bypass_synchronized_trace::Vector{Int64}
end

function _application_phase_values(peak, frequency, time_s; phase_shift_rad=0.0)
    angle = 2.0 * pi * frequency * time_s + phase_shift_rad
    return [peak * sin(angle - (phase - 1) * 2.0 * pi / 3.0) for phase in 1:3]
end

function _application_phase_derivative(peak, frequency, time_s; phase_shift_rad=0.0)
    angle = 2.0 * pi * frequency * time_s + phase_shift_rad
    angular_frequency = 2.0 * pi * frequency
    return [angular_frequency * peak * cos(
        angle - (phase - 1) * 2.0 * pi / 3.0,
    ) for phase in 1:3]
end

function _application_three_wire_limit(values, maximum_phase_voltage_v)
    centered = values .- sum(values) / 3.0
    maximum_value = maximum(abs, centered; init=0.0)
    maximum_value <= maximum_phase_voltage_v && return centered
    return centered .* (maximum_phase_voltage_v / maximum_value)
end

function _application_positive_voltage(value, identity)
    isfinite(value) && value > 1.0e-6 || throw(ArgumentError(
        "$identity collapsed outside the admitted positive-voltage domain",
    ))
    return value
end

function _converter_application_same_boundary(time_s, event_time_s, fixed_step_s)
    return isapprox(
        time_s,
        event_time_s;
        atol=16.0 * eps(Float64) * max(abs(event_time_s), fixed_step_s),
        rtol=0.0,
    )
end

function _converter_application_boundary_has_passed(
    time_s,
    event_time_s,
    fixed_step_s,
    event_side,
)
    _converter_application_same_boundary(time_s, event_time_s, fixed_step_s) &&
        return event_side === :right
    return time_s > event_time_s
end


function _uninterruptible_power_supply_bypass_synchronization(parameters)
    transfer_time = parameters.bypass_transfer_time_s
    phase_error = rem2pi(
        parameters.bypass_source_phase_shift_rad + 2.0 * pi *
            (parameters.bypass_source_frequency_hz - parameters.output_frequency_hz) *
            transfer_time,
        RoundNearest,
    )
    voltage_mismatch = abs(
        parameters.bypass_source_phase_voltage_peak_v -
            parameters.output_phase_voltage_peak_v,
    ) / parameters.output_phase_voltage_peak_v
    frequency_mismatch = abs(
        parameters.bypass_source_frequency_hz - parameters.output_frequency_hz,
    )
    accepted = voltage_mismatch <= parameters.maximum_bypass_voltage_mismatch_pu &&
        frequency_mismatch <= parameters.maximum_bypass_frequency_mismatch_hz &&
        abs(phase_error) <= parameters.maximum_bypass_phase_mismatch_rad
    return (; accepted, voltage_mismatch, frequency_mismatch, phase_error)
end

function _uninterruptible_power_supply_primary_source_available(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    loss_has_passed = _converter_application_boundary_has_passed(
        time_s,
        parameters.source_loss_time_s,
        fixed_step_s,
        event_side,
    )
    recovery_has_passed = _converter_application_boundary_has_passed(
        time_s,
        parameters.source_recovery_time_s,
        fixed_step_s,
        event_side,
    )
    return !loss_has_passed || recovery_has_passed
end

function _uninterruptible_power_supply_operating_mode(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    _converter_application_boundary_has_passed(
        time_s,
        parameters.double_conversion_restore_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationRestoredOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.bypass_transfer_time_s,
        fixed_step_s,
        event_side,
    ) && _uninterruptible_power_supply_bypass_synchronization(parameters).accepted &&
        return ConverterSystems.ConverterApplicationSynchronizedBypass
    _converter_application_boundary_has_passed(
        time_s,
        parameters.source_loss_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationStoredEnergySupport
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _conductive_electric_vehicle_charger_operating_mode(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    _converter_application_boundary_has_passed(
        time_s,
        parameters.output_short_clear_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationRestoredOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.output_short_start_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationOutputTerminalShortCircuit
    _converter_application_boundary_has_passed(
        time_s,
        parameters.load_step_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationLoadStepOperation
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _dynamic_voltage_restorer_operating_mode(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    _converter_application_boundary_has_passed(
        time_s,
        parameters.sag_stop_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationRestoredOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.bypass_stop_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationVoltageInjectionOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.bypass_start_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationBypassOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.sag_start_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationVoltageInjectionOperation
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _shunt_active_filter_operating_mode(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    _converter_application_boundary_has_passed(
        time_s,
        parameters.reference_restore_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationRestoredOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.reference_loss_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationReferenceLossOperation
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _solid_state_transformer_operating_mode(
    parameters,
    time_s,
    event_side,
    fixed_step_s,
)
    _converter_application_boundary_has_passed(
        time_s,
        parameters.transformer_side_fault_clear_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationRestoredOperation
    _converter_application_boundary_has_passed(
        time_s,
        parameters.transformer_side_fault_start_time_s,
        fixed_step_s,
        event_side,
    ) && return ConverterSystems.ConverterApplicationTransformerSideFaultOperation
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _average_converter_application_operating_mode(study, time_s, event_side)
    parameters = study.parameters
    parameters isa ConverterSystems.UninterruptiblePowerSupplyParameters &&
        return _uninterruptible_power_supply_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        )
    parameters isa ConverterSystems.ConductiveElectricVehicleChargerParameters &&
        return _conductive_electric_vehicle_charger_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        )
    parameters isa ConverterSystems.DynamicVoltageRestorerParameters &&
        return _dynamic_voltage_restorer_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        )
    parameters isa ConverterSystems.ShuntActiveFilterParameters &&
        return _shunt_active_filter_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        )
    parameters isa ConverterSystems.SolidStateTransformerParameters &&
        return _solid_state_transformer_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        )
    return ConverterSystems.ConverterApplicationNormalOperation
end

function _average_converter_application_event_cursor(study, time_s, event_side)
    parameters = study.parameters
    parameters isa ConverterSystems.DynamicVoltageRestorerParameters && return count(
        event_time -> _converter_application_boundary_has_passed(
            time_s,
            event_time,
            study.specification.timing.fixed_step_s,
            event_side,
        ),
        (
            parameters.sag_start_time_s,
            parameters.bypass_start_time_s,
            parameters.bypass_stop_time_s,
            parameters.sag_stop_time_s,
        ),
    )
    parameters isa ConverterSystems.UninterruptiblePowerSupplyParameters && return count(
        event_time -> _converter_application_boundary_has_passed(
            time_s,
            event_time,
            study.specification.timing.fixed_step_s,
            event_side,
        ),
        (
            parameters.source_loss_time_s,
            parameters.bypass_transfer_time_s,
            parameters.source_recovery_time_s,
            parameters.double_conversion_restore_time_s,
        ),
    )
    parameters isa ConverterSystems.ConductiveElectricVehicleChargerParameters &&
        return count(
            event_time -> _converter_application_boundary_has_passed(
                time_s,
                event_time,
                study.specification.timing.fixed_step_s,
                event_side,
            ),
            (
                parameters.load_step_time_s,
                parameters.output_short_start_time_s,
                parameters.output_short_clear_time_s,
            ),
        )
    parameters isa ConverterSystems.ShuntActiveFilterParameters && return count(
        event_time -> _converter_application_boundary_has_passed(
            time_s,
            event_time,
            study.specification.timing.fixed_step_s,
            event_side,
        ),
        (parameters.reference_loss_time_s, parameters.reference_restore_time_s),
    )
    parameters isa ConverterSystems.SolidStateTransformerParameters && return count(
        event_time -> _converter_application_boundary_has_passed(
            time_s,
            event_time,
            study.specification.timing.fixed_step_s,
            event_side,
        ),
        (
            parameters.transformer_side_fault_start_time_s,
            parameters.transformer_side_fault_clear_time_s,
        ),
    )
    return 0
end

function _application_input_controller(
    grid_voltage,
    input_current,
    dc_link_voltage,
    dc_link_reference_voltage,
    requested_output_power,
    resistance,
    inductance,
    bandwidth,
    dc_gain,
    modulation_limit,
)
    voltage = _application_positive_voltage(dc_link_voltage, "application DC link")
    requested_input_power = requested_output_power +
        dc_gain * (dc_link_reference_voltage - voltage)
    voltage_norm = sum(abs2, grid_voltage)
    voltage_norm > eps(Float64) || throw(ArgumentError(
        "converter application input controller encountered zero grid-voltage norm",
    ))
    reference_current = requested_input_power .* grid_voltage ./ voltage_norm
    requested_derivative = bandwidth .* (reference_current .- input_current)
    unconstrained_converter_voltage = grid_voltage .- resistance .* input_current .-
        inductance .* requested_derivative
    converter_voltage = _application_three_wire_limit(
        unconstrained_converter_voltage,
        0.5 * modulation_limit * voltage,
    )
    derivative = (grid_voltage .- converter_voltage .- resistance .* input_current) ./
        inductance
    return (; derivative, converter_voltage, reference_current)
end

function _application_initial_vector(initial::ConverterSystems.ShuntActiveFilterInitialState)
    return [initial.filter_current_a..., initial.dc_link_voltage_v,
        initial.current_control_integral_a_s...]
end

function _application_initial_vector(initial::ConverterSystems.DynamicVoltageRestorerInitialState)
    return [initial.load_current_a..., initial.dc_link_voltage_v]
end

function _application_initial_vector(initial::ConverterSystems.UninterruptiblePowerSupplyInitialState)
    return [initial.input_filter_current_a..., initial.dc_link_voltage_v,
        initial.output_filter_current_a...]
end

function _application_initial_vector(
    initial::ConverterSystems.ConductiveElectricVehicleChargerInitialState,
)
    return [initial.input_filter_current_a..., initial.dc_link_voltage_v,
        initial.output_inductor_current_a, initial.output_capacitor_voltage_v]
end

function _application_initial_vector(initial::ConverterSystems.SolidStateTransformerInitialState)
    return [initial.input_filter_current_a..., initial.primary_dc_link_voltage_v,
        initial.transformer_current_a, initial.secondary_dc_link_voltage_v,
        initial.output_filter_current_a...]
end

function _application_observation(study, time_s, state,
    parameters::ConverterSystems.ShuntActiveFilterParameters;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    filter_current = state[1:3]
    dc_voltage = _application_positive_voltage(state[4], "shunt-filter DC link")
    control_integral = state[5:7]
    grid_voltage = _application_phase_values(
        parameters.grid_phase_voltage_peak_v,
        parameters.grid_frequency_hz,
        time_s,
    )
    load_fundamental = _application_phase_values(
        parameters.load_fundamental_current_peak_a,
        parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=-parameters.load_phase_shift_rad,
    )
    fifth_harmonic = _application_phase_values(
        parameters.load_fifth_harmonic_current_peak_a,
        5.0 * parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=parameters.load_phase_shift_rad,
    )
    seventh_harmonic = _application_phase_values(
        parameters.load_seventh_harmonic_current_peak_a,
        7.0 * parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=-parameters.load_phase_shift_rad,
    )
    load_current = load_fundamental .+ fifth_harmonic .+ seventh_harmonic
    active_source_current = _application_phase_values(
        parameters.load_fundamental_current_peak_a * cos(parameters.load_phase_shift_rad),
        parameters.grid_frequency_hz,
        time_s,
    )
    mode = isnothing(operating_mode) ?
        _shunt_active_filter_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        ) : operating_mode
    compensation_available = mode !==
        ConverterSystems.ConverterApplicationReferenceLossOperation
    base_reference = compensation_available ?
        load_current .- active_source_current : zeros(3)
    unit_in_phase = _application_phase_values(1.0, parameters.grid_frequency_hz, time_s)
    dc_charge_current = parameters.dc_voltage_control_a_per_v *
        (parameters.dc_link_reference_voltage_v - dc_voltage)
    reference_current = base_reference .- dc_charge_current .* unit_in_phase
    available_base_derivative = _application_phase_derivative(
        parameters.load_fundamental_current_peak_a,
        parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=-parameters.load_phase_shift_rad,
    ) .+ _application_phase_derivative(
        parameters.load_fifth_harmonic_current_peak_a,
        5.0 * parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=parameters.load_phase_shift_rad,
    ) .+ _application_phase_derivative(
        parameters.load_seventh_harmonic_current_peak_a,
        7.0 * parameters.grid_frequency_hz,
        time_s;
        phase_shift_rad=-parameters.load_phase_shift_rad,
    ) .- _application_phase_derivative(
        parameters.load_fundamental_current_peak_a * cos(parameters.load_phase_shift_rad),
        parameters.grid_frequency_hz,
        time_s,
    )
    base_derivative = compensation_available ? available_base_derivative : zeros(3)
    reference_derivative = base_derivative .- dc_charge_current .*
        _application_phase_derivative(1.0, parameters.grid_frequency_hz, time_s)
    error = reference_current .- filter_current
    unconstrained_voltage = grid_voltage .+
        parameters.filter_resistance_ohm .* filter_current .+
        parameters.filter_inductance_h .* reference_derivative .+
        parameters.current_control_resistance_ohm .* error .+
        parameters.current_control_integral_ohm_per_s .* control_integral
    converter_voltage = _application_three_wire_limit(
        unconstrained_voltage,
        0.5 * parameters.modulation_limit * dc_voltage,
    )
    current_derivative = (converter_voltage .- grid_voltage .-
        parameters.filter_resistance_ohm .* filter_current) ./
        parameters.filter_inductance_h
    converter_power = dot(converter_voltage, filter_current)
    dc_derivative = -converter_power /
        (parameters.dc_link_capacitance_f * dc_voltage)
    source_current = load_current .- filter_current
    stored_energy = 0.5 * parameters.filter_inductance_h * sum(abs2, filter_current) +
        0.5 * parameters.dc_link_capacitance_f * dc_voltage^2
    dissipated_power = parameters.filter_resistance_ohm * sum(abs2, filter_current)
    external_input_power = -dot(grid_voltage, filter_current)
    return (
        derivative=[current_derivative..., dc_derivative, error...],
        input_voltage=grid_voltage,
        input_current=-filter_current,
        output_voltage=grid_voltage,
        output_current=source_current,
        converter_current=filter_current,
        dc_link_voltage=[dc_voltage, 0.0],
        stage_power=[converter_power, 0.0, 0.0],
        stored_energy,
        dissipated_power,
        external_input_power,
        control_error=error,
        operating_mode=mode,
        load_delivery=true,
        bypass_synchronized=false,
    )
end

function _application_observation(study, time_s, state,
    parameters::ConverterSystems.DynamicVoltageRestorerParameters;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    load_current = state[1:3]
    dc_voltage = _application_positive_voltage(state[4], "DVR DC link")
    event_side in (:left, :right) || throw(ArgumentError(
        "converter-application event side must be :left or :right",
    ))
    step = study.specification.timing.fixed_step_s
    mode = operating_mode === nothing ?
        _dynamic_voltage_restorer_operating_mode(
            parameters,
            time_s,
            event_side,
            step,
        ) : operating_mode
    mode isa ConverterSystems.ConverterApplicationOperatingMode ||
        throw(ArgumentError("DVR operating mode has the wrong type"))
    same_boundary(event_time_s) = isapprox(
        time_s,
        event_time_s;
        atol=16.0 * eps(Float64) * max(abs(event_time_s), step),
        rtol=0.0,
    )
    at_start = same_boundary(parameters.sag_start_time_s)
    at_stop = same_boundary(parameters.sag_stop_time_s)
    sag_active = if at_start
        event_side === :right
    elseif at_stop
        event_side === :left
    else
        parameters.sag_start_time_s < time_s < parameters.sag_stop_time_s
    end
    retained = sag_active ?
        parameters.sag_retained_voltage_pu : 1.0
    source_voltage = _application_phase_values(
        retained * parameters.source_phase_voltage_peak_v,
        parameters.source_frequency_hz,
        time_s,
    )
    target_voltage = _application_phase_values(
        parameters.target_phase_voltage_peak_v,
        parameters.source_frequency_hz,
        time_s,
    )
    injection_error = target_voltage .- source_voltage
    requested_injection_voltage = _application_three_wire_limit(
        parameters.injection_control_gain .* injection_error,
        parameters.series_transformer_ratio * 0.5 * parameters.modulation_limit *
            dc_voltage,
    )
    injection_voltage = mode === ConverterSystems.ConverterApplicationBypassOperation ?
        zeros(3) : requested_injection_voltage
    load_voltage = source_voltage .+ injection_voltage
    current_derivative = (load_voltage .-
        parameters.load_resistance_ohm .* load_current) ./ parameters.load_inductance_h
    injection_power = dot(injection_voltage, load_current)
    dc_derivative = -injection_power /
        (parameters.dc_link_capacitance_f * dc_voltage)
    stored_energy = 0.5 * parameters.load_inductance_h * sum(abs2, load_current) +
        0.5 * parameters.dc_link_capacitance_f * dc_voltage^2
    dissipated_power = parameters.load_resistance_ohm * sum(abs2, load_current)
    delivered_voltage_pu = sqrt(sum(abs2, load_voltage) /
        (1.5 * parameters.target_phase_voltage_peak_v^2))
    return (
        derivative=[current_derivative..., dc_derivative],
        input_voltage=source_voltage,
        input_current=load_current,
        output_voltage=load_voltage,
        output_current=load_current,
        converter_current=load_current,
        dc_link_voltage=[dc_voltage, 0.0],
        stage_power=[injection_power, 0.0, 0.0],
        stored_energy,
        dissipated_power,
        external_input_power=dot(source_voltage, load_current),
        control_error=target_voltage .- load_voltage,
        operating_mode=mode,
        load_delivery=delivered_voltage_pu >= 0.9,
        bypass_synchronized=false,
    )
end

function _application_observation(study, time_s, state,
    parameters::ConverterSystems.UninterruptiblePowerSupplyParameters;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    input_current = state[1:3]
    dc_voltage = _application_positive_voltage(state[4], "UPS DC link")
    output_current = state[5:7]
    event_side in (:left, :right) || throw(ArgumentError(
        "converter-application event side must be :left or :right",
    ))
    fixed_step = study.specification.timing.fixed_step_s
    mode = operating_mode === nothing ?
        _uninterruptible_power_supply_operating_mode(
            parameters,
            time_s,
            event_side,
            fixed_step,
        ) : operating_mode
    mode isa ConverterSystems.ConverterApplicationOperatingMode ||
        throw(ArgumentError("UPS operating mode has the wrong type"))
    source_available = _uninterruptible_power_supply_primary_source_available(
        parameters,
        time_s,
        event_side,
        fixed_step,
    )
    source_voltage = source_available ? _application_phase_values(
        parameters.source_phase_voltage_peak_v,
        parameters.source_frequency_hz,
        time_s,
    ) : zeros(3)
    requested_output_voltage = _application_phase_values(
        parameters.output_phase_voltage_peak_v,
        parameters.output_frequency_hz,
        time_s,
    )
    inverter_output_voltage = _application_three_wire_limit(
        requested_output_voltage,
        0.5 * parameters.modulation_limit * dc_voltage,
    )
    bypass_synchronization =
        _uninterruptible_power_supply_bypass_synchronization(parameters)
    bypass_voltage = _application_phase_values(
        parameters.bypass_source_phase_voltage_peak_v,
        parameters.bypass_source_frequency_hz,
        time_s;
        phase_shift_rad=parameters.bypass_source_phase_shift_rad,
    )
    bypass_active = mode === ConverterSystems.ConverterApplicationSynchronizedBypass
    output_voltage = bypass_active ? bypass_voltage : inverter_output_voltage
    output_derivative = (output_voltage .-
        parameters.load_resistance_ohm .* output_current) ./ parameters.load_inductance_h
    inverter_power = bypass_active ? 0.0 : dot(output_voltage, output_current)
    input_control = if source_available
        _application_input_controller(
            source_voltage,
            input_current,
            dc_voltage,
            parameters.dc_link_reference_voltage_v,
            max(inverter_power, 0.0),
            parameters.input_filter_resistance_ohm,
            parameters.input_filter_inductance_h,
            parameters.input_current_control_bandwidth_rad_s,
            parameters.dc_voltage_control_w_per_v,
            parameters.modulation_limit,
        )
    else
        decay = -parameters.input_filter_resistance_ohm .* input_current ./
            parameters.input_filter_inductance_h
        (; derivative=decay, converter_voltage=zeros(3), reference_current=zeros(3))
    end
    rectifier_power = source_available ?
        dot(input_control.converter_voltage, input_current) : 0.0
    dc_derivative = (rectifier_power - inverter_power) /
        (parameters.dc_link_capacitance_f * dc_voltage)
    stored_energy = 0.5 * parameters.input_filter_inductance_h *
        sum(abs2, input_current) +
        0.5 * parameters.dc_link_capacitance_f * dc_voltage^2 +
        0.5 * parameters.load_inductance_h * sum(abs2, output_current)
    dissipated_power = parameters.input_filter_resistance_ohm *
        sum(abs2, input_current) +
        parameters.load_resistance_ohm * sum(abs2, output_current)
    bypass_power = bypass_active ? dot(bypass_voltage, output_current) : 0.0
    delivered_voltage_pu = sqrt(sum(abs2, output_voltage) /
        (1.5 * parameters.output_phase_voltage_peak_v^2))
    load_delivery = delivered_voltage_pu >=
        parameters.minimum_load_delivery_voltage_pu
    return (
        derivative=[input_control.derivative..., dc_derivative, output_derivative...],
        input_voltage=source_voltage,
        input_current,
        output_voltage,
        output_current,
        converter_current=input_control.reference_current,
        dc_link_voltage=[dc_voltage, 0.0],
        stage_power=[rectifier_power, inverter_power, bypass_power],
        stored_energy,
        dissipated_power,
        external_input_power=dot(source_voltage, input_current) + bypass_power,
        control_error=[parameters.dc_link_reference_voltage_v - dc_voltage,
            maximum(abs, requested_output_voltage .- output_voltage),
            bypass_synchronization.phase_error],
        operating_mode=mode,
        load_delivery,
        bypass_synchronized=bypass_synchronization.accepted,
    )
end

function _application_observation(study, time_s, state,
    parameters::ConverterSystems.ConductiveElectricVehicleChargerParameters;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    input_current = state[1:3]
    dc_voltage = _application_positive_voltage(state[4], "EV-charger DC link")
    output_inductor_current = state[5]
    output_voltage = max(state[6], 0.0)
    event_side in (:left, :right) || throw(ArgumentError(
        "converter-application event side must be :left or :right",
    ))
    mode = operating_mode === nothing ?
        _conductive_electric_vehicle_charger_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        ) : operating_mode
    mode isa ConverterSystems.ConverterApplicationOperatingMode ||
        throw(ArgumentError("EV-charger operating mode has the wrong type"))
    output_load_resistance =
        mode === ConverterSystems.ConverterApplicationOutputTerminalShortCircuit ?
            parameters.output_short_resistance_ohm :
        mode in (
            ConverterSystems.ConverterApplicationLoadStepOperation,
            ConverterSystems.ConverterApplicationRestoredOperation,
        ) ? parameters.stepped_output_load_resistance_ohm :
        parameters.output_load_resistance_ohm
    source_voltage = _application_phase_values(
        parameters.source_phase_voltage_peak_v,
        parameters.source_frequency_hz,
        time_s,
    )
    load_power = output_voltage^2 / output_load_resistance
    input_control = _application_input_controller(
        source_voltage,
        input_current,
        dc_voltage,
        parameters.dc_link_reference_voltage_v,
        load_power,
        parameters.input_filter_resistance_ohm,
        parameters.input_filter_inductance_h,
        parameters.input_current_control_bandwidth_rad_s,
        parameters.dc_voltage_control_w_per_v,
        parameters.modulation_limit,
    )
    rectifier_power = dot(input_control.converter_voltage, input_current)
    voltage_error = parameters.output_voltage_reference_v - output_voltage
    duty = clamp(
        (output_voltage + parameters.output_voltage_control_a_per_v * voltage_error) /
            dc_voltage,
        0.0,
        parameters.modulation_limit,
    )
    chopper_voltage = duty * dc_voltage
    output_inductor_derivative = (chopper_voltage - output_voltage) /
        parameters.output_inductance_h
    output_capacitor_derivative = (output_inductor_current -
        output_voltage / output_load_resistance) /
        parameters.output_capacitance_f
    chopper_power = chopper_voltage * output_inductor_current
    dc_derivative = (rectifier_power - chopper_power) /
        (parameters.dc_link_capacitance_f * dc_voltage)
    stored_energy = 0.5 * parameters.input_filter_inductance_h *
        sum(abs2, input_current) +
        0.5 * parameters.dc_link_capacitance_f * dc_voltage^2 +
        0.5 * parameters.output_inductance_h * output_inductor_current^2 +
        0.5 * parameters.output_capacitance_f * output_voltage^2
    dissipated_power = parameters.input_filter_resistance_ohm *
        sum(abs2, input_current) + load_power
    return (
        derivative=[input_control.derivative..., dc_derivative,
            output_inductor_derivative, output_capacitor_derivative],
        input_voltage=source_voltage,
        input_current,
        output_voltage=[output_voltage, 0.0, 0.0],
        output_current=[output_voltage / output_load_resistance,
            0.0, 0.0],
        converter_current=[output_inductor_current, 0.0, 0.0],
        dc_link_voltage=[dc_voltage, 0.0],
        stage_power=[rectifier_power, chopper_power, load_power],
        stored_energy,
        dissipated_power,
        external_input_power=dot(source_voltage, input_current),
        control_error=[parameters.dc_link_reference_voltage_v - dc_voltage,
            voltage_error, duty],
        operating_mode=mode,
        load_delivery=mode !==
            ConverterSystems.ConverterApplicationOutputTerminalShortCircuit &&
            output_voltage >= 0.9 * parameters.output_voltage_reference_v,
        bypass_synchronized=false,
    )
end

function _application_observation(study, time_s, state,
    parameters::ConverterSystems.SolidStateTransformerParameters;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    input_current = state[1:3]
    primary_dc_voltage = _application_positive_voltage(state[4], "SST primary DC link")
    transformer_current = state[5]
    secondary_dc_voltage = _application_positive_voltage(state[6], "SST secondary DC link")
    output_current = state[7:9]
    source_voltage = _application_phase_values(
        parameters.source_phase_voltage_peak_v,
        parameters.source_frequency_hz,
        time_s,
    )
    requested_output_voltage = _application_phase_values(
        parameters.output_phase_voltage_peak_v,
        parameters.output_frequency_hz,
        time_s,
    )
    output_voltage = _application_three_wire_limit(
        requested_output_voltage,
        0.5 * parameters.modulation_limit * secondary_dc_voltage,
    )
    output_derivative = (output_voltage .-
        parameters.load_resistance_ohm .* output_current) ./ parameters.load_inductance_h
    inverter_power = dot(output_voltage, output_current)
    mode = isnothing(operating_mode) ?
        _solid_state_transformer_operating_mode(
            parameters,
            time_s,
            event_side,
            study.specification.timing.fixed_step_s,
        ) : operating_mode
    transformer_side_fault_current = mode ===
        ConverterSystems.ConverterApplicationTransformerSideFaultOperation ?
        secondary_dc_voltage / parameters.transformer_side_fault_resistance_ohm : 0.0
    transformer_side_fault_power = secondary_dc_voltage *
        transformer_side_fault_current
    transformer_drive_voltage = clamp(
        secondary_dc_voltage + parameters.secondary_dc_voltage_control_v_per_v *
            (parameters.secondary_dc_link_reference_voltage_v - secondary_dc_voltage),
        0.0,
        parameters.modulation_limit * parameters.transformer_ratio * primary_dc_voltage,
    )
    transformer_derivative = (transformer_drive_voltage - secondary_dc_voltage -
        parameters.transformer_resistance_ohm * transformer_current) /
        parameters.transformer_leakage_inductance_h
    transformer_input_power = transformer_drive_voltage * transformer_current
    transformer_output_power = secondary_dc_voltage * transformer_current
    input_control = _application_input_controller(
        source_voltage,
        input_current,
        primary_dc_voltage,
        parameters.primary_dc_link_reference_voltage_v,
        max(transformer_input_power, inverter_power, 0.0),
        parameters.input_filter_resistance_ohm,
        parameters.input_filter_inductance_h,
        parameters.input_current_control_bandwidth_rad_s,
        parameters.primary_dc_voltage_control_w_per_v,
        parameters.modulation_limit,
    )
    rectifier_power = dot(input_control.converter_voltage, input_current)
    primary_dc_derivative = (rectifier_power - transformer_input_power) /
        (parameters.primary_dc_link_capacitance_f * primary_dc_voltage)
    secondary_dc_derivative = (transformer_output_power - inverter_power -
        transformer_side_fault_power) /
        (parameters.secondary_dc_link_capacitance_f * secondary_dc_voltage)
    stored_energy = 0.5 * parameters.input_filter_inductance_h *
        sum(abs2, input_current) +
        0.5 * parameters.primary_dc_link_capacitance_f * primary_dc_voltage^2 +
        0.5 * parameters.transformer_leakage_inductance_h * transformer_current^2 +
        0.5 * parameters.secondary_dc_link_capacitance_f * secondary_dc_voltage^2 +
        0.5 * parameters.load_inductance_h * sum(abs2, output_current)
    dissipated_power = parameters.input_filter_resistance_ohm *
        sum(abs2, input_current) +
        parameters.transformer_resistance_ohm * transformer_current^2 +
        parameters.load_resistance_ohm * sum(abs2, output_current) +
        transformer_side_fault_power
    return (
        derivative=[input_control.derivative..., primary_dc_derivative,
            transformer_derivative, secondary_dc_derivative, output_derivative...],
        input_voltage=source_voltage,
        input_current,
        output_voltage,
        output_current,
        converter_current=[transformer_current, transformer_current, transformer_current],
        dc_link_voltage=[primary_dc_voltage, secondary_dc_voltage],
        stage_power=[rectifier_power, transformer_output_power, inverter_power],
        stored_energy,
        dissipated_power,
        external_input_power=dot(source_voltage, input_current),
        control_error=[parameters.primary_dc_link_reference_voltage_v - primary_dc_voltage,
            parameters.secondary_dc_link_reference_voltage_v - secondary_dc_voltage,
            maximum(abs, requested_output_voltage .- output_voltage)],
        operating_mode=mode,
        load_delivery=sqrt(sum(abs2, output_voltage) /
            (1.5 * parameters.output_phase_voltage_peak_v^2)) >= 0.9,
        bypass_synchronized=false,
    )
end

function _application_observation(
    study,
    time_s,
    state;
    event_side::Symbol=:right,
    operating_mode=nothing,
)
    if study.parameters isa Union{
        ConverterSystems.ShuntActiveFilterParameters,
        ConverterSystems.UninterruptiblePowerSupplyParameters,
        ConverterSystems.ConductiveElectricVehicleChargerParameters,
        ConverterSystems.DynamicVoltageRestorerParameters,
        ConverterSystems.SolidStateTransformerParameters,
    }
        return _application_observation(
            study,
            time_s,
            state,
            study.parameters;
            event_side,
            operating_mode,
        )
    end
    return _application_observation(study, time_s, state, study.parameters; event_side)
end

function _average_converter_application_events(runtime)
    parameters = runtime.study.parameters
    events = ConverterSystems.ConverterSystemEventRecord[]
    step = runtime.study.specification.timing.fixed_step_s
    accepted(event_time_s) = runtime.time_s >= event_time_s -
        16.0 * eps(Float64) * max(abs(event_time_s), step)
    if parameters isa ConverterSystems.DynamicVoltageRestorerParameters
        if accepted(parameters.sag_start_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.sag_start_time_s,
                :source_voltage_sag_applied,
                :dynamic_voltage_restorer_source,
                true;
                message="accepted source-voltage sag boundary and entered injection mode",
            ))
        end
        if accepted(parameters.bypass_start_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.bypass_start_time_s,
                :dynamic_voltage_restorer_bypass_applied,
                :dynamic_voltage_restorer_load_path,
                true;
                message="accepted physical bypass and removed series injection",
            ))
        end
        if accepted(parameters.bypass_stop_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.bypass_stop_time_s,
                :dynamic_voltage_restorer_bypass_cleared,
                :dynamic_voltage_restorer_load_path,
                true;
                message="cleared physical bypass and restored series injection",
            ))
        end
        if accepted(parameters.sag_stop_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.sag_stop_time_s,
                :source_voltage_sag_cleared,
                :dynamic_voltage_restorer_source,
                true;
                message="accepted source-voltage recovery boundary and left injection mode",
            ))
        end
    elseif parameters isa ConverterSystems.UninterruptiblePowerSupplyParameters
        synchronization =
            _uninterruptible_power_supply_bypass_synchronization(parameters)
        if accepted(parameters.source_loss_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.source_loss_time_s,
                :uninterruptible_power_supply_source_lost,
                :uninterruptible_power_supply_primary_source,
                true;
                message="accepted primary-source loss and entered stored-energy support",
            ))
        end
        if accepted(parameters.bypass_transfer_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.bypass_transfer_time_s,
                synchronization.accepted ?
                    :uninterruptible_power_supply_bypass_transfer :
                    :uninterruptible_power_supply_bypass_refused,
                :uninterruptible_power_supply_bypass_source,
                synchronization.accepted;
                message=synchronization.accepted ?
                    "accepted synchronized bypass transfer" :
                    "refused bypass transfer outside voltage, frequency, or phase synchronization limits",
            ))
        end
        if accepted(parameters.source_recovery_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.source_recovery_time_s,
                :uninterruptible_power_supply_source_recovered,
                :uninterruptible_power_supply_primary_source,
                true;
                message="accepted primary-source recovery while preserving the active load path",
            ))
        end
        if accepted(parameters.double_conversion_restore_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.double_conversion_restore_time_s,
                :uninterruptible_power_supply_double_conversion_restored,
                :uninterruptible_power_supply_load_path,
                true;
                message="accepted return to the double-conversion load path",
            ))
        end
    elseif parameters isa ConverterSystems.ConductiveElectricVehicleChargerParameters
        if accepted(parameters.load_step_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.load_step_time_s,
                :electric_vehicle_charger_load_step_applied,
                :electric_vehicle_charger_dc_terminal,
                true;
                message="accepted the declared DC-terminal load-resistance step",
            ))
        end
        if accepted(parameters.output_short_start_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.output_short_start_time_s,
                :electric_vehicle_charger_output_short_applied,
                :electric_vehicle_charger_dc_terminal,
                true;
                message="accepted the finite-resistance DC-terminal short circuit",
            ))
        end
        if accepted(parameters.output_short_clear_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.output_short_clear_time_s,
                :electric_vehicle_charger_output_short_cleared,
                :electric_vehicle_charger_dc_terminal,
                true;
                message="cleared the DC-terminal short and restored the stepped load",
            ))
        end
    elseif parameters isa ConverterSystems.ShuntActiveFilterParameters
        if accepted(parameters.reference_loss_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.reference_loss_time_s,
                :shunt_active_filter_reference_lost,
                :shunt_active_filter_compensation_reference,
                true;
                message="accepted compensation-reference loss while retaining DC-link regulation",
            ))
        end
        if accepted(parameters.reference_restore_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.reference_restore_time_s,
                :shunt_active_filter_reference_restored,
                :shunt_active_filter_compensation_reference,
                true;
                message="accepted restoration of harmonic-current compensation",
            ))
        end
    elseif parameters isa ConverterSystems.SolidStateTransformerParameters
        if accepted(parameters.transformer_side_fault_start_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.transformer_side_fault_start_time_s,
                :solid_state_transformer_secondary_terminal_fault_applied,
                :solid_state_transformer_medium_frequency_transformer_secondary_terminal,
                true;
                message="accepted finite-resistance transformer-secondary terminal fault",
            ))
        end
        if accepted(parameters.transformer_side_fault_clear_time_s)
            push!(events, ConverterSystems.ConverterSystemEventRecord(
                parameters.transformer_side_fault_clear_time_s,
                :solid_state_transformer_secondary_terminal_fault_cleared,
                :solid_state_transformer_medium_frequency_transformer_secondary_terminal,
                true;
                message="cleared transformer-secondary terminal fault and restored stage operation",
            ))
        end
    end
    return events
end

function prepare_average_converter_application(
    study::ConverterSystems.AverageConverterApplicationStudy,
)
    step = study.specification.timing.fixed_step_s
    sample_count = round(Int, (study.stop_time_s - study.start_time_s) / step) + 1
    initial_mode = _average_converter_application_operating_mode(
        study,
        study.start_time_s,
        :right,
    )
    runtime = AverageConverterApplicationRuntime(
        study,
        _application_initial_vector(study.initial_state),
        study.start_time_s,
        0,
        initial_mode,
        _average_converter_application_event_cursor(
            study,
            study.start_time_s,
            :right,
        ),
        0.0,
        0.0,
        zeros(sample_count),
        zeros(3, sample_count),
        zeros(3, sample_count),
        zeros(3, sample_count),
        zeros(3, sample_count),
        zeros(3, sample_count),
        zeros(2, sample_count),
        zeros(3, sample_count),
        zeros(sample_count),
        zeros(sample_count),
        zeros(sample_count),
        zeros(3, sample_count),
        zeros(Int64, sample_count),
        zeros(Int64, sample_count),
        zeros(Int64, sample_count),
    )
    _record_average_converter_application!(
        runtime,
        1,
        _application_observation(
            study,
            study.start_time_s,
            runtime.state;
            operating_mode=runtime.operating_mode,
        ),
        0.0,
    )
    return runtime
end

function _average_converter_application_study_signature(study)
    return bytes2hex(sha256(join((
        study.specification.signature_sha256,
        string(typeof(study.parameters)),
        repr(study.parameters),
        string(typeof(study.initial_state)),
        repr(study.initial_state),
        repr(study.start_time_s),
        repr(study.stop_time_s),
    ), '\n')))
end

function _average_converter_application_snapshot_field(
    identity,
    family,
    unit,
    value;
    axes=String[],
)
    encoded = value isa AbstractArray ? PortableSnapshots.portable_snapshot_array(
        value;
        unit,
        axes,
    ) : value
    return PortableSnapshots.PortableSnapshotStateField(
        identity,
        _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_OWNER,
        family,
        unit,
        axes,
        encoded,
    )
end

function _average_converter_application_snapshot_invariants(runtime)
    final_step = length(runtime.time_trace_s) - 1
    0 <= runtime.accepted_step_index <= final_step ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application accepted step lies outside its configured horizon",
        )
    expected_time = runtime.study.start_time_s + runtime.accepted_step_index *
        runtime.study.specification.timing.fixed_step_s
    runtime.time_s == expected_time || _average_converter_application_snapshot_fail(
        :backend_state_boundary,
        "converter-application time differs from its accepted fixed-step index",
    )
    expected_mode = _average_converter_application_operating_mode(
        runtime.study,
        runtime.time_s,
        :right,
    )
    runtime.operating_mode === expected_mode ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application operating mode differs from its accepted event calendar",
        )
    expected_event_cursor = _average_converter_application_event_cursor(
        runtime.study,
        runtime.time_s,
        :right,
    )
    runtime.application_event_cursor == expected_event_cursor ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application event cursor differs from its accepted event calendar",
        )
    prefix = 1:(runtime.accepted_step_index + 1)
    values = (
        runtime.state,
        @view(runtime.time_trace_s[prefix]),
        @view(runtime.input_voltage_trace_v[:, prefix]),
        @view(runtime.input_current_trace_a[:, prefix]),
        @view(runtime.output_voltage_trace_v[:, prefix]),
        @view(runtime.output_current_trace_a[:, prefix]),
        @view(runtime.converter_current_trace_a[:, prefix]),
        @view(runtime.dc_link_voltage_trace_v[:, prefix]),
        @view(runtime.stage_power_trace_w[:, prefix]),
        @view(runtime.stored_energy_trace_j[prefix]),
        @view(runtime.dissipated_power_trace_w[prefix]),
        @view(runtime.energy_residual_trace_w[prefix]),
        @view(runtime.control_error_trace[:, prefix]),
        @view(runtime.operating_mode_trace[prefix]),
        @view(runtime.load_delivery_trace[prefix]),
        @view(runtime.bypass_synchronized_trace[prefix]),
    )
    all(array -> all(isfinite, array), values) &&
        all(isfinite, (runtime.time_s, runtime.input_energy_j,
            runtime.dissipated_energy_j)) ||
        _average_converter_application_snapshot_fail(
            :nonfinite_value,
            "converter-application accepted state contains a nonfinite value",
        )
    expected_trace = [
        runtime.study.start_time_s + index *
            runtime.study.specification.timing.fixed_step_s for
        index in 0:runtime.accepted_step_index
    ]
    runtime.time_trace_s[prefix] == expected_trace ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application output cursor does not follow its fixed-step grid",
        )
    return prefix
end

function average_converter_application_snapshot(runtime::AverageConverterApplicationRuntime)
    prefix = _average_converter_application_snapshot_invariants(runtime)
    fields = PortableSnapshots.PortableSnapshotStateField[
        _average_converter_application_snapshot_field(
            "converter_application.schema_version", :checkpoint, "1",
            _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_SCHEMA,
        ),
        _average_converter_application_snapshot_field(
            "converter_application.study_signature", :checkpoint, "1",
            _average_converter_application_study_signature(runtime.study),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.application_id", :checkpoint, "1",
            Int64(Int(runtime.study.specification.selection.application)),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.sample_capacity", :checkpoint, "1",
            Int64(length(runtime.time_trace_s)),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.accepted_step_index", :discrete, "1",
            Int64(runtime.accepted_step_index),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.operating_mode", :discrete, "1",
            Int64(runtime.operating_mode),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.application_event_cursor", :discrete, "1",
            Int64(runtime.application_event_cursor),
        ),
        _average_converter_application_snapshot_field(
            "converter_application.time_s", :continuous, "s", runtime.time_s,
        ),
        _average_converter_application_snapshot_field(
            "converter_application.input_energy_j", :continuous, "J",
            runtime.input_energy_j,
        ),
        _average_converter_application_snapshot_field(
            "converter_application.dissipated_energy_j", :continuous, "J",
            runtime.dissipated_energy_j,
        ),
        _average_converter_application_snapshot_field(
            "converter_application.state", :continuous, "model_coordinate",
            copy(runtime.state); axes=["state_coordinate"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.time_trace_s", :output, "s",
            copy(@view(runtime.time_trace_s[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.input_voltage_trace_v", :output, "V",
            copy(@view(runtime.input_voltage_trace_v[:, prefix]));
            axes=["phase", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.input_current_trace_a", :output, "A",
            copy(@view(runtime.input_current_trace_a[:, prefix]));
            axes=["phase", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.output_voltage_trace_v", :output, "V",
            copy(@view(runtime.output_voltage_trace_v[:, prefix]));
            axes=["phase", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.output_current_trace_a", :output, "A",
            copy(@view(runtime.output_current_trace_a[:, prefix]));
            axes=["phase", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.converter_current_trace_a", :output, "A",
            copy(@view(runtime.converter_current_trace_a[:, prefix]));
            axes=["phase", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.dc_link_voltage_trace_v", :output, "V",
            copy(@view(runtime.dc_link_voltage_trace_v[:, prefix]));
            axes=["dc_link", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.stage_power_trace_w", :output, "W",
            copy(@view(runtime.stage_power_trace_w[:, prefix]));
            axes=["stage", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.stored_energy_trace_j", :output, "J",
            copy(@view(runtime.stored_energy_trace_j[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.dissipated_power_trace_w", :output, "W",
            copy(@view(runtime.dissipated_power_trace_w[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.energy_residual_trace_w", :output, "W",
            copy(@view(runtime.energy_residual_trace_w[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.control_error_trace", :output,
            "model_coordinate", copy(@view(runtime.control_error_trace[:, prefix]));
            axes=["control_coordinate", "sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.operating_mode_trace", :output, "1",
            copy(@view(runtime.operating_mode_trace[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.load_delivery_trace", :output, "1",
            copy(@view(runtime.load_delivery_trace[prefix])); axes=["sample"],
        ),
        _average_converter_application_snapshot_field(
            "converter_application.bypass_synchronized_trace", :output, "1",
            copy(@view(runtime.bypass_synchronized_trace[prefix])); axes=["sample"],
        ),
    ]
    return PortableSnapshots.PortableSnapshotStateInventory(fields)
end

function _average_converter_application_snapshot_value(
    fields,
    identity,
    family,
    unit,
    axes=String[],
)
    field = get(fields, identity, nothing)
    field === nothing && _average_converter_application_snapshot_fail(
        :missing_backend_state,
        "converter-application snapshot field $identity is missing",
    )
    field.owner == _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_OWNER &&
        field.family === family && field.unit == unit && field.axes == axes ||
        _average_converter_application_snapshot_fail(
            :backend_state_schema,
            "converter-application snapshot field $identity has incompatible metadata",
        )
    return field.value
end

function _average_converter_application_snapshot_array(
    fields,
    identity,
    family,
    unit,
    axes,
    expected_shape,
    expected_element_type::Type=Float64,
)
    encoded = _average_converter_application_snapshot_value(
        fields, identity, family, unit, axes,
    )
    encoded isa PortableSnapshots.PortableSnapshotArray ||
        _average_converter_application_snapshot_fail(
            :backend_state_type,
            "converter-application snapshot field $identity is not a portable array",
        )
    values = PortableSnapshots.portable_snapshot_array_values(encoded)
    eltype(values) === expected_element_type && size(values) == expected_shape &&
        all(isfinite, values) ||
        _average_converter_application_snapshot_fail(
            :backend_state_shape,
            "converter-application snapshot field $identity has invalid shape or values",
        )
    return values
end

function _average_converter_application_restore_data(runtime, snapshot)
    fields = Dict(field.identity => field for field in snapshot.fields)
    required_identities = Set((
        "converter_application.schema_version",
        "converter_application.study_signature",
        "converter_application.application_id",
        "converter_application.sample_capacity",
        "converter_application.accepted_step_index",
        "converter_application.operating_mode",
        "converter_application.application_event_cursor",
        "converter_application.time_s",
        "converter_application.input_energy_j",
        "converter_application.dissipated_energy_j",
        "converter_application.state",
        "converter_application.time_trace_s",
        "converter_application.input_voltage_trace_v",
        "converter_application.input_current_trace_a",
        "converter_application.output_voltage_trace_v",
        "converter_application.output_current_trace_a",
        "converter_application.converter_current_trace_a",
        "converter_application.dc_link_voltage_trace_v",
        "converter_application.stage_power_trace_w",
        "converter_application.stored_energy_trace_j",
        "converter_application.dissipated_power_trace_w",
        "converter_application.energy_residual_trace_w",
        "converter_application.control_error_trace",
        "converter_application.operating_mode_trace",
        "converter_application.load_delivery_trace",
        "converter_application.bypass_synchronized_trace",
    ))
    Set(keys(fields)) == required_identities ||
        _average_converter_application_snapshot_fail(
            :backend_state_schema,
            "converter-application snapshot has missing or unknown fields",
        )
    scalar(identity, family, unit) = _average_converter_application_snapshot_value(
        fields, identity, family, unit,
    )
    schema = scalar("converter_application.schema_version", :checkpoint, "1")
    schema === _AVERAGE_CONVERTER_APPLICATION_SNAPSHOT_SCHEMA ||
        _average_converter_application_snapshot_fail(
            :backend_state_schema,
            "converter-application snapshot schema is unsupported",
        )
    study_signature = scalar(
        "converter_application.study_signature", :checkpoint, "1",
    )
    study_signature == _average_converter_application_study_signature(runtime.study) ||
        _average_converter_application_snapshot_fail(
            :backend_inventory_mismatch,
            "converter-application snapshot belongs to a different prepared study",
        )
    application_id = scalar("converter_application.application_id", :checkpoint, "1")
    application_id === Int64(Int(runtime.study.specification.selection.application)) ||
        _average_converter_application_snapshot_fail(
            :backend_inventory_mismatch,
            "converter-application snapshot application identity changed",
        )
    sample_capacity = scalar(
        "converter_application.sample_capacity", :checkpoint, "1",
    )
    sample_capacity === Int64(length(runtime.time_trace_s)) ||
        _average_converter_application_snapshot_fail(
            :backend_state_shape,
            "converter-application snapshot horizon changed",
        )
    accepted_step_index = scalar(
        "converter_application.accepted_step_index", :discrete, "1",
    )
    accepted_step_index isa Int64 && 0 <= accepted_step_index < sample_capacity ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application snapshot accepted step is invalid",
        )
    operating_mode_value = scalar(
        "converter_application.operating_mode", :discrete, "1",
    )
    operating_mode_value isa Int64 && operating_mode_value in Int64.(Int.(instances(
        ConverterSystems.ConverterApplicationOperatingMode,
    ))) || _average_converter_application_snapshot_fail(
        :backend_state_boundary,
        "converter-application snapshot operating mode is invalid",
    )
    operating_mode = ConverterSystems.ConverterApplicationOperatingMode(
        operating_mode_value,
    )
    application_event_cursor = scalar(
        "converter_application.application_event_cursor", :discrete, "1",
    )
    application_event_cursor isa Int64 && application_event_cursor >= 0 ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application snapshot event cursor is invalid",
        )
    prefix_length = Int(accepted_step_index) + 1
    time_s = scalar("converter_application.time_s", :continuous, "s")
    input_energy_j = scalar(
        "converter_application.input_energy_j", :continuous, "J",
    )
    dissipated_energy_j = scalar(
        "converter_application.dissipated_energy_j", :continuous, "J",
    )
    all(value -> value isa Float64 && isfinite(value),
        (time_s, input_energy_j, dissipated_energy_j)) ||
        _average_converter_application_snapshot_fail(
            :nonfinite_value,
            "converter-application snapshot scalar state is invalid",
        )
    expected_time = runtime.study.start_time_s + Int(accepted_step_index) *
        runtime.study.specification.timing.fixed_step_s
    time_s == expected_time || _average_converter_application_snapshot_fail(
        :backend_state_boundary,
        "converter-application snapshot time differs from its accepted step",
    )
    state = _average_converter_application_snapshot_array(
        fields, "converter_application.state", :continuous, "model_coordinate",
        ["state_coordinate"], size(runtime.state),
    )
    time_trace_s = _average_converter_application_snapshot_array(
        fields, "converter_application.time_trace_s", :output, "s", ["sample"],
        (prefix_length,),
    )
    expected_trace = [
        runtime.study.start_time_s + index *
            runtime.study.specification.timing.fixed_step_s for
        index in 0:Int(accepted_step_index)
    ]
    vec(time_trace_s) == expected_trace || _average_converter_application_snapshot_fail(
        :backend_state_boundary,
        "converter-application snapshot output cursor does not follow its fixed-step grid",
    )
    matrix(identity, unit, rows, axis) = _average_converter_application_snapshot_array(
        fields, identity, :output, unit, [axis, "sample"],
        (rows, prefix_length),
    )
    vector(identity, unit) = _average_converter_application_snapshot_array(
        fields, identity, :output, unit, ["sample"], (prefix_length,),
    )
    integer_vector(identity) = _average_converter_application_snapshot_array(
        fields,
        identity,
        :output,
        "1",
        ["sample"],
        (prefix_length,),
        Int64,
    )
    operating_mode_trace = vec(integer_vector(
        "converter_application.operating_mode_trace",
    ))
    all(value -> value in Int64.(Int.(instances(
        ConverterSystems.ConverterApplicationOperatingMode,
    ))), operating_mode_trace) || _average_converter_application_snapshot_fail(
        :backend_state_boundary,
        "converter-application snapshot mode trace contains an invalid mode",
    )
    load_delivery_trace = vec(integer_vector(
        "converter_application.load_delivery_trace",
    ))
    bypass_synchronized_trace = vec(integer_vector(
        "converter_application.bypass_synchronized_trace",
    ))
    all(value -> value in (Int64(0), Int64(1)), load_delivery_trace) &&
        all(value -> value in (Int64(0), Int64(1)), bypass_synchronized_trace) ||
        _average_converter_application_snapshot_fail(
            :backend_state_boundary,
            "converter-application snapshot Boolean state traces are invalid",
        )
    return (
        accepted_step_index=Int(accepted_step_index),
        operating_mode,
        application_event_cursor=Int(application_event_cursor),
        time_s,
        input_energy_j,
        dissipated_energy_j,
        state=vec(state),
        time_trace_s=vec(time_trace_s),
        input_voltage_trace_v=matrix(
            "converter_application.input_voltage_trace_v", "V", 3, "phase"),
        input_current_trace_a=matrix(
            "converter_application.input_current_trace_a", "A", 3, "phase"),
        output_voltage_trace_v=matrix(
            "converter_application.output_voltage_trace_v", "V", 3, "phase"),
        output_current_trace_a=matrix(
            "converter_application.output_current_trace_a", "A", 3, "phase"),
        converter_current_trace_a=matrix(
            "converter_application.converter_current_trace_a", "A", 3, "phase"),
        dc_link_voltage_trace_v=matrix(
            "converter_application.dc_link_voltage_trace_v", "V", 2, "dc_link"),
        stage_power_trace_w=matrix(
            "converter_application.stage_power_trace_w", "W", 3, "stage"),
        stored_energy_trace_j=vector(
            "converter_application.stored_energy_trace_j", "J"),
        dissipated_power_trace_w=vector(
            "converter_application.dissipated_power_trace_w", "W"),
        energy_residual_trace_w=vector(
            "converter_application.energy_residual_trace_w", "W"),
        control_error_trace=matrix(
            "converter_application.control_error_trace", "model_coordinate", 3,
            "control_coordinate"),
        operating_mode_trace,
        load_delivery_trace,
        bypass_synchronized_trace,
    )
end

function _apply_average_converter_application_restore_data!(runtime, data)
    prefix = 1:(data.accepted_step_index + 1)
    runtime.state .= data.state
    runtime.time_s = data.time_s
    runtime.accepted_step_index = data.accepted_step_index
    runtime.operating_mode = data.operating_mode
    runtime.application_event_cursor = data.application_event_cursor
    runtime.input_energy_j = data.input_energy_j
    runtime.dissipated_energy_j = data.dissipated_energy_j
    for trace in (
        runtime.time_trace_s,
        runtime.input_voltage_trace_v,
        runtime.input_current_trace_a,
        runtime.output_voltage_trace_v,
        runtime.output_current_trace_a,
        runtime.converter_current_trace_a,
        runtime.dc_link_voltage_trace_v,
        runtime.stage_power_trace_w,
        runtime.stored_energy_trace_j,
        runtime.dissipated_power_trace_w,
        runtime.energy_residual_trace_w,
        runtime.control_error_trace,
        runtime.operating_mode_trace,
        runtime.load_delivery_trace,
        runtime.bypass_synchronized_trace,
    )
        fill!(trace, 0.0)
    end
    runtime.time_trace_s[prefix] .= data.time_trace_s
    runtime.input_voltage_trace_v[:, prefix] .= data.input_voltage_trace_v
    runtime.input_current_trace_a[:, prefix] .= data.input_current_trace_a
    runtime.output_voltage_trace_v[:, prefix] .= data.output_voltage_trace_v
    runtime.output_current_trace_a[:, prefix] .= data.output_current_trace_a
    runtime.converter_current_trace_a[:, prefix] .= data.converter_current_trace_a
    runtime.dc_link_voltage_trace_v[:, prefix] .= data.dc_link_voltage_trace_v
    runtime.stage_power_trace_w[:, prefix] .= data.stage_power_trace_w
    runtime.stored_energy_trace_j[prefix] .= data.stored_energy_trace_j
    runtime.dissipated_power_trace_w[prefix] .= data.dissipated_power_trace_w
    runtime.energy_residual_trace_w[prefix] .= data.energy_residual_trace_w
    runtime.control_error_trace[:, prefix] .= data.control_error_trace
    runtime.operating_mode_trace[prefix] .= data.operating_mode_trace
    runtime.load_delivery_trace[prefix] .= data.load_delivery_trace
    runtime.bypass_synchronized_trace[prefix] .= data.bypass_synchronized_trace
    return runtime
end

function restore_average_converter_application_snapshot!(
    runtime::AverageConverterApplicationRuntime,
    snapshot::PortableSnapshots.PortableSnapshotStateInventory,
)
    data = _average_converter_application_restore_data(runtime, snapshot)
    probe = prepare_average_converter_application(runtime.study)
    _apply_average_converter_application_restore_data!(probe, data)
    reconstructed = average_converter_application_snapshot(probe)
    reconstructed.signature_sha256 == snapshot.signature_sha256 ||
        _average_converter_application_snapshot_fail(
            :backend_reconstruction,
            "converter-application snapshot changed during isolated reconstruction",
        )
    _apply_average_converter_application_restore_data!(runtime, data)
    return runtime
end

function _record_average_converter_application!(runtime, sample, observation,
    energy_residual_w)
    runtime.time_trace_s[sample] = runtime.time_s
    runtime.input_voltage_trace_v[:, sample] .= observation.input_voltage
    runtime.input_current_trace_a[:, sample] .= observation.input_current
    runtime.output_voltage_trace_v[:, sample] .= observation.output_voltage
    runtime.output_current_trace_a[:, sample] .= observation.output_current
    runtime.converter_current_trace_a[:, sample] .= observation.converter_current
    runtime.dc_link_voltage_trace_v[:, sample] .= observation.dc_link_voltage
    runtime.stage_power_trace_w[:, sample] .= observation.stage_power
    runtime.stored_energy_trace_j[sample] = observation.stored_energy
    runtime.dissipated_power_trace_w[sample] = observation.dissipated_power
    runtime.energy_residual_trace_w[sample] = energy_residual_w
    runtime.control_error_trace[:, sample] .= observation.control_error
    operating_mode = hasproperty(observation, :operating_mode) ?
        observation.operating_mode :
        ConverterSystems.ConverterApplicationNormalOperation
    load_delivery = hasproperty(observation, :load_delivery) ?
        observation.load_delivery : true
    bypass_synchronized = hasproperty(observation, :bypass_synchronized) ?
        observation.bypass_synchronized : false
    runtime.operating_mode_trace[sample] = Int64(operating_mode)
    runtime.load_delivery_trace[sample] = Int64(load_delivery)
    runtime.bypass_synchronized_trace[sample] = Int64(bypass_synchronized)
    return runtime
end

function _advance_average_converter_application!(runtime)
    study = runtime.study
    final_step = length(runtime.time_trace_s) - 1
    runtime.accepted_step_index < final_step || return false
    next_step = runtime.accepted_step_index + 1
    sample = next_step + 1
    step = study.specification.timing.fixed_step_s
    time = runtime.time_s
    state = runtime.state
    first_observation = _application_observation(
        study,
        time,
        state;
        operating_mode=runtime.operating_mode,
    )
    second_state = state .+ 0.5 * step .* first_observation.derivative
    second_observation = _application_observation(study, time + 0.5 * step, second_state)
    third_state = state .+ 0.5 * step .* second_observation.derivative
    third_observation = _application_observation(study, time + 0.5 * step, third_state)
    fourth_state = state .+ step .* third_observation.derivative
    endpoint_time = time + step
    fourth_observation = _application_observation(
        study,
        endpoint_time,
        fourth_state,
        event_side=:left,
    )
    previous_energy = first_observation.stored_energy
    runtime.state .= state .+ step / 6.0 .* (
        first_observation.derivative .+ 2.0 .* second_observation.derivative .+
        2.0 .* third_observation.derivative .+ fourth_observation.derivative
    )
    runtime.time_s = study.start_time_s + (sample - 1) * step
    runtime.operating_mode = _average_converter_application_operating_mode(
        study,
        runtime.time_s,
        :right,
    )
    runtime.application_event_cursor = _average_converter_application_event_cursor(
        study,
        runtime.time_s,
        :right,
    )
    final_observation = _application_observation(
        study,
        runtime.time_s,
        runtime.state;
        operating_mode=runtime.operating_mode,
    )
    average_input_power = (
        first_observation.external_input_power +
        2.0 * second_observation.external_input_power +
        2.0 * third_observation.external_input_power +
        fourth_observation.external_input_power
    ) / 6.0
    average_dissipated_power = (
        first_observation.dissipated_power +
        2.0 * second_observation.dissipated_power +
        2.0 * third_observation.dissipated_power +
        fourth_observation.dissipated_power
    ) / 6.0
    energy_residual = average_input_power - average_dissipated_power -
        (final_observation.stored_energy - previous_energy) / step
    runtime.input_energy_j += step * average_input_power
    runtime.dissipated_energy_j += step * average_dissipated_power
    runtime.accepted_step_index = next_step
    _record_average_converter_application!(
        runtime,
        sample,
        final_observation,
        energy_residual,
    )
    return true
end

function advance_average_converter_application!(
    runtime::AverageConverterApplicationRuntime,
    accepted_step_count::Integer=1,
)
    accepted_step_count >= 0 || throw(ArgumentError(
        "converter-application accepted step count must be nonnegative",
    ))
    final_step = length(runtime.time_trace_s) - 1
    runtime.accepted_step_index + accepted_step_count <= final_step || throw(ArgumentError(
        "converter-application advancement exceeds the configured horizon",
    ))
    for _ in 1:accepted_step_count
        _advance_average_converter_application!(runtime) || throw(ArgumentError(
            "converter-application runtime completed before the requested step count",
        ))
    end
    return runtime
end

function _average_converter_application_result(runtime)
    study = runtime.study
    observation = _application_observation(
        study,
        runtime.time_s,
        runtime.state;
        operating_mode=runtime.operating_mode,
    )
    events = _average_converter_application_events(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(runtime.state),
        repr(observation.stored_energy),
        string(runtime.accepted_step_index),
        string(Int(runtime.operating_mode)),
        string(runtime.application_event_cursor),
        repr(events),
    ), '\n')))
    capacitor_charge = Float64[]
    inductor_flux = Float64[]
    parameters = study.parameters
    if parameters isa ConverterSystems.ShuntActiveFilterParameters
        push!(capacitor_charge,
            parameters.dc_link_capacitance_f * runtime.state[4])
        append!(inductor_flux, parameters.filter_inductance_h .* runtime.state[1:3])
    elseif parameters isa ConverterSystems.DynamicVoltageRestorerParameters
        push!(capacitor_charge,
            parameters.dc_link_capacitance_f * runtime.state[4])
        append!(inductor_flux, parameters.load_inductance_h .* runtime.state[1:3])
    elseif parameters isa ConverterSystems.UninterruptiblePowerSupplyParameters
        push!(capacitor_charge,
            parameters.dc_link_capacitance_f * runtime.state[4])
        append!(inductor_flux,
            parameters.input_filter_inductance_h .* runtime.state[1:3])
        append!(inductor_flux, parameters.load_inductance_h .* runtime.state[5:7])
    elseif parameters isa ConverterSystems.ConductiveElectricVehicleChargerParameters
        append!(capacitor_charge, (
            parameters.dc_link_capacitance_f * runtime.state[4],
            parameters.output_capacitance_f * runtime.state[6],
        ))
        append!(inductor_flux,
            parameters.input_filter_inductance_h .* runtime.state[1:3])
        push!(inductor_flux, parameters.output_inductance_h * runtime.state[5])
    else
        append!(capacitor_charge, (
            parameters.primary_dc_link_capacitance_f * runtime.state[4],
            parameters.secondary_dc_link_capacitance_f * runtime.state[6],
        ))
        append!(inductor_flux,
            parameters.input_filter_inductance_h .* runtime.state[1:3])
        push!(inductor_flux,
            parameters.transformer_leakage_inductance_h * runtime.state[5])
        append!(inductor_flux, parameters.load_inductance_h .* runtime.state[7:9])
    end
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [observation.input_voltage..., observation.output_voltage...],
        [observation.input_current..., observation.output_current...],
        BitVector(),
        BitVector(),
        BitVector(),
        capacitor_charge,
        inductor_flux,
        Float64[],
        copy(runtime.state),
        observation.stored_energy,
        runtime.dissipated_energy_j,
        runtime.accepted_step_index,
        length(events),
        signature,
    )
    prefix = 1:(runtime.accepted_step_index + 1)
    power_scale = max(maximum(abs, @view(runtime.stage_power_trace_w[:, prefix])),
        maximum(abs, @view(runtime.dissipated_power_trace_w[prefix])), 1.0)
    return ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=true,
        status=:ok,
        events,
        maximum_kcl_residual_a=0.0,
        relative_charge_residual=0.0,
        relative_energy_residual=maximum(abs,
            @view(runtime.energy_residual_trace_w[prefix])) /
            power_scale,
        harmonic_metrics=Dict(
            :application_stage_count => Float64(length(study.topologies)),
            :average_value_fidelity => 1.0,
        ),
    )
end

function execute_average_converter_application!(
    runtime::AverageConverterApplicationRuntime,
)
    while _advance_average_converter_application!(runtime)
    end
    return ConverterSystems.AverageConverterApplicationTrace(
        runtime.time_trace_s,
        runtime.input_voltage_trace_v,
        runtime.input_current_trace_a,
        runtime.output_voltage_trace_v,
        runtime.output_current_trace_a,
        runtime.converter_current_trace_a,
        runtime.dc_link_voltage_trace_v,
        runtime.stage_power_trace_w,
        runtime.stored_energy_trace_j,
        runtime.dissipated_power_trace_w,
        runtime.energy_residual_trace_w,
        runtime.control_error_trace,
        ConverterSystems.ConverterApplicationOperatingMode.(
            runtime.operating_mode_trace,
        ),
        BitVector(runtime.load_delivery_trace .== 1),
        BitVector(runtime.bypass_synchronized_trace .== 1),
        _average_converter_application_result(runtime),
    )
end

export AbstractAverageConverterApplicationParameters,
       ShuntActiveFilterParameters,
       DynamicVoltageRestorerParameters,
       UninterruptiblePowerSupplyParameters,
       ConductiveElectricVehicleChargerParameters,
       SolidStateTransformerParameters,
       ShuntActiveFilterInitialState,
       DynamicVoltageRestorerInitialState,
       UninterruptiblePowerSupplyInitialState,
       ConductiveElectricVehicleChargerInitialState,
       SolidStateTransformerInitialState,
       ConverterApplicationOperatingMode,
       ConverterApplicationNormalOperation,
       ConverterApplicationStoredEnergySupport,
       ConverterApplicationSynchronizedBypass,
       ConverterApplicationLoadStepOperation,
       ConverterApplicationOutputTerminalShortCircuit,
       ConverterApplicationVoltageInjectionOperation,
       ConverterApplicationBypassOperation,
       ConverterApplicationReferenceLossOperation,
       ConverterApplicationTransformerSideFaultOperation,
       ConverterApplicationRestoredOperation,
       AverageConverterApplicationStudy,
       AverageConverterApplicationTrace

abstract type AbstractAverageConverterApplicationParameters end

Base.@kwdef struct ShuntActiveFilterParameters <:
                       AbstractAverageConverterApplicationParameters
    grid_phase_voltage_peak_v::Float64
    grid_frequency_hz::Float64
    load_fundamental_current_peak_a::Float64
    load_fifth_harmonic_current_peak_a::Float64
    load_seventh_harmonic_current_peak_a::Float64
    load_phase_shift_rad::Float64 = 0.0
    filter_resistance_ohm::Float64
    filter_inductance_h::Float64
    dc_link_capacitance_f::Float64
    dc_link_reference_voltage_v::Float64
    current_control_resistance_ohm::Float64
    current_control_integral_ohm_per_s::Float64
    dc_voltage_control_a_per_v::Float64
    reference_loss_time_s::Float64
    reference_restore_time_s::Float64
    modulation_limit::Float64 = 0.95
end

Base.@kwdef struct DynamicVoltageRestorerParameters <:
                       AbstractAverageConverterApplicationParameters
    source_phase_voltage_peak_v::Float64
    source_frequency_hz::Float64
    sag_start_time_s::Float64
    sag_stop_time_s::Float64
    sag_retained_voltage_pu::Float64
    bypass_start_time_s::Float64
    bypass_stop_time_s::Float64
    target_phase_voltage_peak_v::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    dc_link_capacitance_f::Float64
    dc_link_reference_voltage_v::Float64
    injection_control_gain::Float64
    series_transformer_ratio::Float64
    modulation_limit::Float64 = 0.95
end

Base.@kwdef struct UninterruptiblePowerSupplyParameters <:
                       AbstractAverageConverterApplicationParameters
    source_phase_voltage_peak_v::Float64
    source_frequency_hz::Float64
    input_filter_resistance_ohm::Float64
    input_filter_inductance_h::Float64
    dc_link_capacitance_f::Float64
    dc_link_reference_voltage_v::Float64
    output_phase_voltage_peak_v::Float64
    output_frequency_hz::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    input_current_control_bandwidth_rad_s::Float64
    dc_voltage_control_w_per_v::Float64
    source_loss_time_s::Float64
    bypass_transfer_time_s::Float64
    source_recovery_time_s::Float64
    double_conversion_restore_time_s::Float64
    bypass_source_phase_voltage_peak_v::Float64
    bypass_source_frequency_hz::Float64
    bypass_source_phase_shift_rad::Float64 = 0.0
    maximum_bypass_voltage_mismatch_pu::Float64 = 0.05
    maximum_bypass_frequency_mismatch_hz::Float64 = 0.1
    maximum_bypass_phase_mismatch_rad::Float64 = 5.0 * pi / 180.0
    minimum_load_delivery_voltage_pu::Float64 = 0.9
    modulation_limit::Float64 = 0.95
end

@enum ConverterApplicationOperatingMode::Int64 begin
    ConverterApplicationNormalOperation = 1
    ConverterApplicationStoredEnergySupport = 2
    ConverterApplicationSynchronizedBypass = 3
    ConverterApplicationLoadStepOperation = 4
    ConverterApplicationOutputTerminalShortCircuit = 5
    ConverterApplicationVoltageInjectionOperation = 6
    ConverterApplicationBypassOperation = 7
    ConverterApplicationReferenceLossOperation = 8
    ConverterApplicationRestoredOperation = 9
    ConverterApplicationTransformerSideFaultOperation = 10
end

Base.@kwdef struct ConductiveElectricVehicleChargerParameters <:
                       AbstractAverageConverterApplicationParameters
    source_phase_voltage_peak_v::Float64
    source_frequency_hz::Float64
    input_filter_resistance_ohm::Float64
    input_filter_inductance_h::Float64
    dc_link_capacitance_f::Float64
    dc_link_reference_voltage_v::Float64
    output_inductance_h::Float64
    output_capacitance_f::Float64
    output_load_resistance_ohm::Float64
    output_voltage_reference_v::Float64
    input_current_control_bandwidth_rad_s::Float64
    dc_voltage_control_w_per_v::Float64
    output_voltage_control_a_per_v::Float64
    load_step_time_s::Float64
    stepped_output_load_resistance_ohm::Float64
    output_short_start_time_s::Float64
    output_short_clear_time_s::Float64
    output_short_resistance_ohm::Float64
    modulation_limit::Float64 = 0.95
end

Base.@kwdef struct SolidStateTransformerParameters <:
                       AbstractAverageConverterApplicationParameters
    source_phase_voltage_peak_v::Float64
    source_frequency_hz::Float64
    input_filter_resistance_ohm::Float64
    input_filter_inductance_h::Float64
    primary_dc_link_capacitance_f::Float64
    primary_dc_link_reference_voltage_v::Float64
    transformer_ratio::Float64
    transformer_resistance_ohm::Float64
    transformer_leakage_inductance_h::Float64
    secondary_dc_link_capacitance_f::Float64
    secondary_dc_link_reference_voltage_v::Float64
    output_phase_voltage_peak_v::Float64
    output_frequency_hz::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    input_current_control_bandwidth_rad_s::Float64
    primary_dc_voltage_control_w_per_v::Float64
    secondary_dc_voltage_control_v_per_v::Float64
    transformer_side_fault_start_time_s::Float64
    transformer_side_fault_clear_time_s::Float64
    transformer_side_fault_resistance_ohm::Float64
    modulation_limit::Float64 = 0.95
end

struct ShuntActiveFilterInitialState
    filter_current_a::NTuple{3,Float64}
    dc_link_voltage_v::Float64
    current_control_integral_a_s::NTuple{3,Float64}

    function ShuntActiveFilterInitialState(
        filter_current_a=(0.0, 0.0, 0.0),
        dc_link_voltage_v=1.0;
        current_control_integral_a_s=(0.0, 0.0, 0.0),
    )
        current = ntuple(index -> Float64(filter_current_a[index]), 3)
        integral = ntuple(index -> Float64(current_control_integral_a_s[index]), 3)
        voltage = Float64(dc_link_voltage_v)
        all(isfinite, (current..., voltage, integral...)) && voltage > 0.0 ||
            throw(ArgumentError("shunt-filter initial state must be finite with positive DC voltage"))
        return new(current, voltage, integral)
    end
end

struct DynamicVoltageRestorerInitialState
    load_current_a::NTuple{3,Float64}
    dc_link_voltage_v::Float64

    function DynamicVoltageRestorerInitialState(
        load_current_a=(0.0, 0.0, 0.0),
        dc_link_voltage_v=1.0,
    )
        current = ntuple(index -> Float64(load_current_a[index]), 3)
        voltage = Float64(dc_link_voltage_v)
        all(isfinite, (current..., voltage)) && voltage > 0.0 ||
            throw(ArgumentError("DVR initial state must be finite with positive DC voltage"))
        return new(current, voltage)
    end
end

struct UninterruptiblePowerSupplyInitialState
    input_filter_current_a::NTuple{3,Float64}
    dc_link_voltage_v::Float64
    output_filter_current_a::NTuple{3,Float64}

    function UninterruptiblePowerSupplyInitialState(
        input_filter_current_a=(0.0, 0.0, 0.0),
        dc_link_voltage_v=1.0,
        output_filter_current_a=(0.0, 0.0, 0.0),
    )
        input = ntuple(index -> Float64(input_filter_current_a[index]), 3)
        output = ntuple(index -> Float64(output_filter_current_a[index]), 3)
        voltage = Float64(dc_link_voltage_v)
        all(isfinite, (input..., voltage, output...)) && voltage > 0.0 ||
            throw(ArgumentError("UPS initial state must be finite with positive DC voltage"))
        return new(input, voltage, output)
    end
end

struct ConductiveElectricVehicleChargerInitialState
    input_filter_current_a::NTuple{3,Float64}
    dc_link_voltage_v::Float64
    output_inductor_current_a::Float64
    output_capacitor_voltage_v::Float64

    function ConductiveElectricVehicleChargerInitialState(
        input_filter_current_a=(0.0, 0.0, 0.0),
        dc_link_voltage_v=1.0,
        output_inductor_current_a=0.0,
        output_capacitor_voltage_v=0.0,
    )
        input = ntuple(index -> Float64(input_filter_current_a[index]), 3)
        values = Float64.((dc_link_voltage_v, output_inductor_current_a,
            output_capacitor_voltage_v))
        all(isfinite, (input..., values...)) && values[1] > 0.0 && values[3] >= 0.0 ||
            throw(ArgumentError("EV-charger initial state must be finite with positive DC-link and nonnegative output voltage"))
        return new(input, values...)
    end
end

struct SolidStateTransformerInitialState
    input_filter_current_a::NTuple{3,Float64}
    primary_dc_link_voltage_v::Float64
    transformer_current_a::Float64
    secondary_dc_link_voltage_v::Float64
    output_filter_current_a::NTuple{3,Float64}

    function SolidStateTransformerInitialState(
        input_filter_current_a=(0.0, 0.0, 0.0),
        primary_dc_link_voltage_v=1.0,
        transformer_current_a=0.0,
        secondary_dc_link_voltage_v=1.0,
        output_filter_current_a=(0.0, 0.0, 0.0),
    )
        input = ntuple(index -> Float64(input_filter_current_a[index]), 3)
        output = ntuple(index -> Float64(output_filter_current_a[index]), 3)
        values = Float64.((primary_dc_link_voltage_v, transformer_current_a,
            secondary_dc_link_voltage_v))
        all(isfinite, (input..., values..., output...)) && values[1] > 0.0 &&
            values[3] > 0.0 || throw(ArgumentError(
                "SST initial state must be finite with positive primary and secondary DC voltages",
            ))
        return new(input, values..., output)
    end
end

function _average_application_parameter_values(parameters)
    return Float64[getfield(parameters, field) for field in fieldnames(typeof(parameters))]
end

function _application_event_is_on_fixed_step_calendar(
    event_time_s,
    start_time_s,
    fixed_step_s,
)
    event_step = (event_time_s - start_time_s) / fixed_step_s
    return isapprox(event_step, round(event_step); atol=1.0e-10, rtol=1.0e-10)
end

function _validate_average_application_parameters(
    parameters,
    start_time_s,
    stop_time_s,
    fixed_step_s,
)
    values = _average_application_parameter_values(parameters)
    all(isfinite, values) || throw(ArgumentError(
        "converter-application parameters must be finite",
    ))
    parameters.modulation_limit > 0.0 && parameters.modulation_limit <= 1.0 ||
        throw(ArgumentError("converter-application modulation limit must lie in (0, 1]"))
    if parameters isa DynamicVoltageRestorerParameters
        0.0 <= parameters.sag_retained_voltage_pu <= 1.0 || throw(ArgumentError(
            "DVR retained sag voltage must lie in [0, 1]",
        ))
        start_time_s <= parameters.sag_start_time_s < parameters.sag_stop_time_s <=
            stop_time_s || throw(ArgumentError(
                "DVR sag interval must lie inside the study horizon",
            ))
        parameters.sag_start_time_s < parameters.bypass_start_time_s <
            parameters.bypass_stop_time_s < parameters.sag_stop_time_s ||
            throw(ArgumentError(
                "DVR bypass interval must lie strictly inside the sag interval",
            ))
        all(
            event_time -> _application_event_is_on_fixed_step_calendar(
                event_time,
                start_time_s,
                fixed_step_s,
            ),
            (
                parameters.sag_start_time_s,
                parameters.bypass_start_time_s,
                parameters.bypass_stop_time_s,
                parameters.sag_stop_time_s,
            ),
        ) || throw(ArgumentError(
            "DVR sag application and clearance must lie on the accepted fixed-step calendar",
        ))
        positive_fields = setdiff(eachindex(values), (3, 4, 5, 6, 7))
        all(index -> values[index] > 0.0, positive_fields) || throw(ArgumentError(
            "DVR physical and control parameters must be positive",
        ))
    elseif parameters isa UninterruptiblePowerSupplyParameters
        event_times = (
            parameters.source_loss_time_s,
            parameters.bypass_transfer_time_s,
            parameters.source_recovery_time_s,
            parameters.double_conversion_restore_time_s,
        )
        start_time_s <= event_times[1] <= event_times[2] < event_times[3] <=
            event_times[4] <= stop_time_s || throw(ArgumentError(
                "UPS source-loss, bypass-transfer, source-recovery, and double-conversion restoration times must be ordered inside the study horizon",
            ))
        all(
            event_time -> _application_event_is_on_fixed_step_calendar(
                event_time,
                start_time_s,
                fixed_step_s,
            ),
            event_times,
        ) || throw(ArgumentError(
            "UPS source and transfer events must lie on the accepted fixed-step calendar",
        ))
        positive_values = (
            parameters.source_phase_voltage_peak_v,
            parameters.source_frequency_hz,
            parameters.input_filter_resistance_ohm,
            parameters.input_filter_inductance_h,
            parameters.dc_link_capacitance_f,
            parameters.dc_link_reference_voltage_v,
            parameters.output_phase_voltage_peak_v,
            parameters.output_frequency_hz,
            parameters.load_resistance_ohm,
            parameters.load_inductance_h,
            parameters.input_current_control_bandwidth_rad_s,
            parameters.dc_voltage_control_w_per_v,
            parameters.bypass_source_phase_voltage_peak_v,
            parameters.bypass_source_frequency_hz,
        )
        all(>(0.0), positive_values) || throw(ArgumentError(
            "UPS physical, source, load, and control parameters must be positive",
        ))
        abs(parameters.bypass_source_phase_shift_rad) <= pi || throw(ArgumentError(
            "UPS bypass-source phase shift must lie in [-pi, pi]",
        ))
        0.0 <= parameters.maximum_bypass_voltage_mismatch_pu <= 1.0 ||
            throw(ArgumentError(
                "UPS bypass voltage-mismatch limit must lie in [0, 1]",
            ))
        0.0 <= parameters.maximum_bypass_frequency_mismatch_hz <=
            parameters.output_frequency_hz || throw(ArgumentError(
                "UPS bypass frequency-mismatch limit must lie inside the output-frequency domain",
            ))
        0.0 <= parameters.maximum_bypass_phase_mismatch_rad <= pi ||
            throw(ArgumentError(
                "UPS bypass phase-mismatch limit must lie in [0, pi]",
            ))
        0.0 < parameters.minimum_load_delivery_voltage_pu <= 1.0 ||
            throw(ArgumentError(
                "UPS minimum load-delivery voltage must lie in (0, 1]",
            ))
    elseif parameters isa ConductiveElectricVehicleChargerParameters
        event_times = (
            parameters.load_step_time_s,
            parameters.output_short_start_time_s,
            parameters.output_short_clear_time_s,
        )
        start_time_s <= event_times[1] < event_times[2] < event_times[3] <=
            stop_time_s || throw(ArgumentError(
                "EV-charger load-step and output-short events must be ordered inside the study horizon",
            ))
        all(
            event_time -> _application_event_is_on_fixed_step_calendar(
                event_time,
                start_time_s,
                fixed_step_s,
            ),
            event_times,
        ) || throw(ArgumentError(
            "EV-charger load events must lie on the accepted fixed-step calendar",
        ))
        all(>(0.0), values) || throw(ArgumentError(
            "EV-charger physical, load, event, and control parameters must be positive",
        ))
        parameters.output_short_resistance_ohm < min(
            parameters.output_load_resistance_ohm,
            parameters.stepped_output_load_resistance_ohm,
        ) || throw(ArgumentError(
            "EV-charger output-short resistance must be below both admitted load resistances",
        ))
    elseif parameters isa ShuntActiveFilterParameters
        start_time_s <= parameters.reference_loss_time_s <
            parameters.reference_restore_time_s <= stop_time_s ||
            throw(ArgumentError(
                "shunt-filter reference-loss interval must lie inside the study horizon",
            ))
        all(
            event_time -> _application_event_is_on_fixed_step_calendar(
                event_time,
                start_time_s,
                fixed_step_s,
            ),
            (parameters.reference_loss_time_s, parameters.reference_restore_time_s),
        ) || throw(ArgumentError(
            "shunt-filter reference events must lie on the accepted fixed-step calendar",
        ))
        all(>(0.0), values[[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]]) &&
            all(>=(0.0), values[[4, 5]]) || throw(ArgumentError(
                "shunt-filter physical, harmonic, and control parameters are invalid",
            ))
    elseif parameters isa SolidStateTransformerParameters
        start_time_s <= parameters.transformer_side_fault_start_time_s <
            parameters.transformer_side_fault_clear_time_s <= stop_time_s ||
            throw(ArgumentError(
                "SST transformer-side fault interval must lie inside the study horizon",
            ))
        all(
            event_time -> _application_event_is_on_fixed_step_calendar(
                event_time,
                start_time_s,
                fixed_step_s,
            ),
            (
                parameters.transformer_side_fault_start_time_s,
                parameters.transformer_side_fault_clear_time_s,
            ),
        ) || throw(ArgumentError(
            "SST transformer-side fault events must lie on the accepted fixed-step calendar",
        ))
        all(>(0.0), values) || throw(ArgumentError(
            "SST physical, fault, load, and control parameters must be positive",
        ))
    else
        all(>(0.0), values) || throw(ArgumentError(
            "converter-application physical and control parameters must be positive",
        ))
    end
    return nothing
end

function _average_application_contract(application)
    application === ShuntActiveHarmonicFilter && return (
        ShuntActiveFilterParameters,
        ShuntActiveFilterInitialState,
        (:two_level_bridge,),
    )
    application === DynamicVoltageRestorer && return (
        DynamicVoltageRestorerParameters,
        DynamicVoltageRestorerInitialState,
        (:two_level_bridge,),
    )
    application === DoubleConversionUninterruptiblePowerSupply && return (
        UninterruptiblePowerSupplyParameters,
        UninterruptiblePowerSupplyInitialState,
        (:two_level_bridge, :two_level_bridge),
    )
    application === ConductiveElectricVehicleCharger && return (
        ConductiveElectricVehicleChargerParameters,
        ConductiveElectricVehicleChargerInitialState,
        (:two_level_bridge, :step_down_chopper),
    )
    application === ThreeStageSolidStateTransformer && return (
        SolidStateTransformerParameters,
        SolidStateTransformerInitialState,
        (:two_level_bridge, :full_bridge, :full_bridge, :two_level_bridge),
    )
    throw(ArgumentError("average converter application requires a released application identity"))
end

struct AverageConverterApplicationStudy{T<:Tuple,P,I}
    specification::ConverterSystemSpecification
    topologies::T
    parameters::P
    initial_state::I
    start_time_s::Float64
    stop_time_s::Float64

    function AverageConverterApplicationStudy(
        specification::ConverterSystemSpecification;
        topologies,
        parameters::AbstractAverageConverterApplicationParameters,
        initial_state,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.fidelity === AverageValue || throw(ArgumentError(
            "converter application study requires explicit average-value fidelity",
        ))
        selection.application !== StandaloneConversion || throw(ArgumentError(
            "converter application study cannot use standalone conversion",
        ))
        selection.family === ThreePhaseTwoLevelBridge && selection.phase_count == 3 ||
            throw(ArgumentError(
                "released converter applications use a three-phase two-level grid-interface owner",
            ))
        parameter_type, initial_type, topology_families =
            _average_application_contract(selection.application)
        parameters isa parameter_type && initial_state isa initial_type ||
            throw(ArgumentError(
                "converter application parameters or initial state do not match the selected application",
            ))
        topology_values = Tuple(topologies)
        length(topology_values) == length(topology_families) && all(index ->
            topology_values[index] isa BridgeTopologyDescriptor &&
            topology_values[index].family === topology_families[index],
            eachindex(topology_families),
        ) || throw(ArgumentError(
            "converter application omits or misorders a canonical stage topology",
        ))
        specification.topology_signatures ==
            Tuple(bridge_topology_signature.(topology_values)) || throw(ArgumentError(
                "converter application specification does not bind its exact stage topologies",
            ))
        !isempty(specification.passive_and_transformer_signatures) ||
            throw(ArgumentError(
                "converter application requires explicit passive and transformer owner identities",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("converter application specification is not ready"))
        start_time, stop_time = Float64.((start_time_s, stop_time_s))
        isfinite(start_time) && isfinite(stop_time) && start_time >= 0.0 &&
            stop_time > start_time || throw(ArgumentError(
                "converter application horizon must be finite and ordered",
            ))
        step_count = (stop_time - start_time) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "converter application horizon must contain integer fixed steps",
            ))
        _validate_average_application_parameters(
            parameters,
            start_time,
            stop_time,
            specification.timing.fixed_step_s,
        )
        return new{typeof(topology_values),typeof(parameters),typeof(initial_state)}(
            specification,
            topology_values,
            parameters,
            initial_state,
            start_time,
            stop_time,
        )
    end
end

struct AverageConverterApplicationTrace
    time_s::Vector{Float64}
    input_voltage_v::Matrix{Float64}
    input_current_a::Matrix{Float64}
    output_voltage_v::Matrix{Float64}
    output_current_a::Matrix{Float64}
    converter_current_a::Matrix{Float64}
    dc_link_voltage_v::Matrix{Float64}
    stage_power_w::Matrix{Float64}
    stored_energy_j::Vector{Float64}
    dissipated_power_w::Vector{Float64}
    energy_residual_w::Vector{Float64}
    control_error::Matrix{Float64}
    operating_mode::Vector{ConverterApplicationOperatingMode}
    load_delivery_state::BitVector
    bypass_synchronized_state::BitVector
    result::ConverterSystemResult
end

export DetailedConverterSemiconductorParameters,
       DetailedChopperSemiconductorParameters,
       DetailedBuckSemiconductorParameters,
       detailed_chopper_semiconductor_signatures,
       detailed_buck_semiconductor_signatures,
       SwitchingChopperConverterStudy,
       SwitchingBuckConverterStudy,
       SwitchingBoostConverterStudy,
       SwitchingInvertingBuckBoostConverterStudy,
       SwitchingChopperConverterTrace,
       SwitchingBuckConverterTrace,
       SwitchingBoostConverterTrace,
       SwitchingInvertingBuckBoostConverterTrace

struct DetailedConverterSemiconductorParameters
    controlled_on_conductance_s::Float64
    controlled_off_conductance_s::Float64
    controlled_forward_voltage_v::Float64
    freewheel_on_conductance_s::Float64
    freewheel_off_conductance_s::Float64
    freewheel_forward_voltage_v::Float64
    controlled_junction_capacitance_f::Float64
    freewheel_junction_capacitance_f::Float64
    junction_potential_v::Float64
    junction_grading_exponent::Float64
    recovered_charge_lifetime_s::Float64
    turn_off_tail_time_s::Float64
    turn_off_tail_cutoff_current_a::Float64
    snubber_resistance_ohm::Float64
    snubber_capacitance_f::Float64
    switching_energy_current_domain_maximum_a::Float64
    switching_energy_blocking_voltage_domain_maximum_v::Float64
    turn_on_energy_at_domain_maximum_j::Float64
    turn_off_energy_at_domain_maximum_j::Float64
    reverse_recovery_energy_at_domain_maximum_j::Float64
    thermal_capacitance_j_per_k::NTuple{2,Float64}
    thermal_resistance_k_per_w::NTuple{2,Float64}
    ambient_temperature_k::Float64
    junction_voltage_limit_v::Float64
    provenance::ParameterProvenance
end

function DetailedConverterSemiconductorParameters(;
    controlled_on_conductance_s::Real=200.0,
    controlled_off_conductance_s::Real=1.0e-6,
    controlled_forward_voltage_v::Real=1.2,
    freewheel_on_conductance_s::Real=200.0,
    freewheel_off_conductance_s::Real=1.0e-6,
    freewheel_forward_voltage_v::Real=0.9,
    controlled_junction_capacitance_f::Real=1.0e-9,
    freewheel_junction_capacitance_f::Real=2.0e-9,
    junction_potential_v::Real=100.0,
    junction_grading_exponent::Real=0.4,
    recovered_charge_lifetime_s::Real=5.0e-6,
    turn_off_tail_time_s::Real=5.0e-6,
    turn_off_tail_cutoff_current_a::Real=1.0e-3,
    snubber_resistance_ohm::Real=100.0e3,
    snubber_capacitance_f::Real=100.0e-12,
    switching_energy_current_domain_maximum_a::Real=2_000.0,
    switching_energy_blocking_voltage_domain_maximum_v::Real=200.0,
    turn_on_energy_at_domain_maximum_j::Real=0.8e-3,
    turn_off_energy_at_domain_maximum_j::Real=1.0e-3,
    reverse_recovery_energy_at_domain_maximum_j::Real=0.4e-3,
    thermal_capacitance_j_per_k=(0.5, 5.0),
    thermal_resistance_k_per_w=(0.2, 0.8),
    ambient_temperature_k::Real=293.15,
    junction_voltage_limit_v::Real=20.0e3,
    provenance::ParameterProvenance,
)
    scalar_values = Float64.((
        controlled_on_conductance_s,
        controlled_off_conductance_s,
        controlled_forward_voltage_v,
        freewheel_on_conductance_s,
        freewheel_off_conductance_s,
        freewheel_forward_voltage_v,
        controlled_junction_capacitance_f,
        freewheel_junction_capacitance_f,
        junction_potential_v,
        junction_grading_exponent,
        recovered_charge_lifetime_s,
        turn_off_tail_time_s,
        turn_off_tail_cutoff_current_a,
        snubber_resistance_ohm,
        snubber_capacitance_f,
        switching_energy_current_domain_maximum_a,
        switching_energy_blocking_voltage_domain_maximum_v,
        turn_on_energy_at_domain_maximum_j,
        turn_off_energy_at_domain_maximum_j,
        reverse_recovery_energy_at_domain_maximum_j,
        ambient_temperature_k,
        junction_voltage_limit_v,
    ))
    all(isfinite, scalar_values) || throw(ArgumentError(
        "detailed buck semiconductor parameters must be finite",
    ))
    all(>(0.0), scalar_values[[1, 4, 7, 8, 9, 11, 12, 14, 15, 16, 17, 21, 22]]) ||
        throw(ArgumentError("detailed buck positive semiconductor parameters must be positive"))
    0.0 <= scalar_values[2] <= scalar_values[1] || throw(ArgumentError(
        "controlled-device off conductance must lie between zero and its on conductance",
    ))
    0.0 <= scalar_values[5] <= scalar_values[4] || throw(ArgumentError(
        "freewheel-device off conductance must lie between zero and its on conductance",
    ))
    all(>=(0.0), scalar_values[[3, 6, 13, 18, 19, 20]]) || throw(ArgumentError(
        "detailed buck voltage, cutoff-current, and event-energy values must be nonnegative",
    ))
    0.0 <= scalar_values[10] < 1.0 || throw(ArgumentError(
        "detailed buck junction grading exponent must lie in [0, 1)",
    ))
    thermal_capacitance = Tuple(Float64.(thermal_capacitance_j_per_k))
    thermal_resistance = Tuple(Float64.(thermal_resistance_k_per_w))
    length(thermal_capacitance) == 2 && length(thermal_resistance) == 2 ||
        throw(DimensionMismatch("detailed buck thermal model requires exactly two Cauer stages"))
    all(value -> isfinite(value) && value > 0.0,
        (thermal_capacitance..., thermal_resistance...)) || throw(ArgumentError(
        "detailed buck Cauer parameters must be finite and positive",
    ))
    provenance.nature === PhysicalModelParameter || throw(ArgumentError(
        "detailed buck semiconductor provenance must describe physical model parameters",
    ))
    return DetailedConverterSemiconductorParameters(
        scalar_values[1:20]...,
        thermal_capacitance,
        thermal_resistance,
        scalar_values[21],
        scalar_values[22],
        provenance,
    )
end

const DetailedChopperSemiconductorParameters = DetailedConverterSemiconductorParameters
const DetailedBuckSemiconductorParameters = DetailedConverterSemiconductorParameters

function detailed_chopper_semiconductor_signatures(
    parameters::DetailedChopperSemiconductorParameters,
)
    body = repr(parameters)
    return (
        bytes2hex(sha256("controlled_igbt\n" * body)),
        bytes2hex(sha256("freewheel_diode\n" * body)),
    )
end


detailed_buck_semiconductor_signatures(
    parameters::DetailedChopperSemiconductorParameters,
) = detailed_chopper_semiconductor_signatures(parameters)

struct SwitchingChopperConverterStudy{T<:BridgeTopologyDescriptor,D}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    inductor_resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
    load_resistance_ohm::Float64
    initial_state::DCDCConverterInitialState
    detailed_semiconductor::D
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingChopperConverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        inductor_resistance_ohm::Real,
        inductance_h::Real,
        capacitance_f::Real,
        load_resistance_ohm::Real,
        initial_state::DCDCConverterInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family in (BuckChopper, BoostChopper, InvertingBuckBoostChopper) ||
            throw(ArgumentError(
            "switching chopper study requires a canonical buck, boost, or inverting buck-boost family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "switching chopper study requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "switching buck study is a standalone converter, not an application composition",
        ))
        expected_topology_family = selection.family === BuckChopper ?
            :step_down_chopper : selection.family === BoostChopper ?
            :step_up_chopper : :inverting_buck_boost
        topology.family === expected_topology_family || throw(ArgumentError(
            "switching chopper study requires the canonical B200 topology for its selected family",
        ))
        topology_signature = bridge_topology_signature(topology)
        specification.topology_signatures == (topology_signature,) || throw(ArgumentError(
            "switching buck specification does not bind the exact supplied B200 topology identity",
        ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("switching chopper converter specification is not ready"))
        initial_state.inductor_current_a >= 0.0 || throw(ArgumentError(
            "switching chopper initial inductor current must follow its nonnegative reference direction",
        ))
        correctly_signed_output = selection.family === InvertingBuckBoostChopper ?
            initial_state.output_voltage_v <= 0.0 : initial_state.output_voltage_v >= 0.0
        correctly_signed_output || throw(ArgumentError(
            "switching chopper initial output voltage has the wrong sign for its selected family",
        ))
        if selection.fidelity === SwitchingDetailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed chopper execution requires typed semiconductor parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                    "switching-detailed chopper specification does not bind its exact D200 device parameters",
                ))
            specification.timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "switching-detailed chopper timestep must provide at least ten points across its fastest recovered-charge or tail-current state",
            ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state chopper execution cannot claim unused detailed semiconductor parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state chopper execution cannot claim detailed device-fidelity signatures",
            ))
        end
        specification.modulation.kind === CarrierSinusoidalPulseWidthModulation ||
            throw(ArgumentError("switching chopper study requires carrier PWM"))
        timing = specification.timing
        timing.dead_time_s == 0.0 && timing.minimum_pulse_s == 0.0 || throw(ArgumentError(
            "the switching-state chopper slice requires zero gate delay and minimum pulse",
        ))
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            inductor_resistance_ohm,
            inductance_h,
            capacitance_f,
            load_resistance_ohm,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) || throw(ArgumentError("switching chopper study parameters must be finite"))
        values[1] > 0.0 || throw(ArgumentError("switching chopper input voltage must be positive"))
        values[2] > 0.0 || throw(ArgumentError("switching chopper source resistance must be positive"))
        values[3] >= 0.0 || throw(ArgumentError("switching chopper inductor resistance must be nonnegative"))
        all(>(0.0), values[4:6]) || throw(ArgumentError(
            "switching chopper L, C, and load resistance must be positive",
        ))
        if selection.fidelity === SwitchingDetailed
            detailed_semiconductor.switching_energy_current_domain_maximum_a >=
                specification.rated_bases.current_a || throw(ArgumentError(
                    "switching-detailed chopper switching-energy current domain must include rated current",
                ))
            detailed_semiconductor.switching_energy_blocking_voltage_domain_maximum_v >=
                max(values[1], abs(initial_state.output_voltage_v)) || throw(ArgumentError(
                    "switching-detailed chopper switching-energy voltage domain must include its maximum initial blocking voltage",
                ))
        end
        values[7] >= 0.0 && values[8] > values[7] || throw(ArgumentError(
            "switching chopper stop time must be later than its nonnegative start time",
        ))
        step_count = (values[8] - values[7]) / timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("switching chopper horizon must contain an integer number of fixed steps"))
        carrier_steps = inv(timing.carrier_frequency_hz * timing.fixed_step_s)
        isapprox(carrier_steps, round(carrier_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("switching chopper carrier period must contain an integer number of fixed steps"))
        rising_edge_step = carrier_steps * specification.modulation.duty / 2.0
        falling_edge_step = carrier_steps * (1.0 - specification.modulation.duty / 2.0)
        all(edge -> isapprox(edge, round(edge); atol=1.0e-10, rtol=1.0e-10),
            (rising_edge_step, falling_edge_step)) || throw(ArgumentError(
            "switching chopper PWM edges must lie exactly on the fixed-step calendar",
        ))
        validate_converter_system_event_calendar(
            specification;
            start_time_s=values[7],
            stop_time_s=values[8],
            allowed_target_valve_indices=(1,),
        )
        return new{typeof(topology),typeof(detailed_semiconductor)}(
            specification,
            topology,
            values[1:6]...,
            initial_state,
            detailed_semiconductor,
            values[7:8]...,
        )
    end
end

const SwitchingBuckConverterStudy = SwitchingChopperConverterStudy
const SwitchingBoostConverterStudy = SwitchingChopperConverterStudy
const SwitchingInvertingBuckBoostConverterStudy = SwitchingChopperConverterStudy

struct SwitchingChopperConverterTrace
    fidelity::ModelFidelity
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    switch_node_voltage_v::Vector{Float64}
    output_voltage_v::Vector{Float64}
    inductor_current_a::Vector{Float64}
    load_current_a::Vector{Float64}
    requested_gate_state::BitVector
    applied_gate_state::BitVector
    controlled_valve_conducting_state::BitVector
    freewheel_diode_conducting_state::BitVector
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Vector{Float64}
    freewheel_junction_temperature_k::Vector{Float64}
    freewheel_recovered_charge_c::Vector{Float64}
    controlled_turn_off_tail_current_a::Vector{Float64}
    controlled_junction_stored_energy_j::Vector{Float64}
    freewheel_junction_stored_energy_j::Vector{Float64}
    result::ConverterSystemResult
end


const SwitchingBuckConverterTrace = SwitchingChopperConverterTrace
const SwitchingBoostConverterTrace = SwitchingChopperConverterTrace
const SwitchingInvertingBuckBoostConverterTrace = SwitchingChopperConverterTrace

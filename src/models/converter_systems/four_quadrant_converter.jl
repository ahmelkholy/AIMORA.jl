export FullBridgeConverterInitialState,
       FourQuadrantConverterInitialState,
       SinglePhaseTwoLevelInverterInitialState,
       AverageFullBridgeConverterStudy,
       AverageFullBridgeConverterTrace,
       AverageFourQuadrantConverterStudy,
       AverageFourQuadrantConverterTrace,
       AverageSinglePhaseTwoLevelInverterStudy,
       AverageSinglePhaseTwoLevelInverterTrace,
       average_four_quadrant_operating_point,
       SwitchingFullBridgeConverterStudy,
       SwitchingFullBridgeConverterTrace,
       SwitchingFourQuadrantConverterStudy,
       SwitchingFourQuadrantConverterTrace,
       SwitchingSinglePhaseTwoLevelInverterStudy,
       SwitchingSinglePhaseTwoLevelInverterTrace

struct FullBridgeConverterInitialState
    load_current_a::Float64

    function FullBridgeConverterInitialState(load_current_a::Real=0.0)
        current = Float64(load_current_a)
        isfinite(current) || throw(ArgumentError(
            "four-quadrant initial load current must be finite",
        ))
        return new(current)
    end
end

function average_four_quadrant_operating_point(
    input_voltage_v::Real,
    duty::Real,
    source_resistance_ohm::Real,
    load_resistance_ohm::Real,
)
    input_voltage, modulation_duty, source_resistance, load_resistance =
        Float64.((input_voltage_v, duty, source_resistance_ohm, load_resistance_ohm))
    input_voltage > 0.0 && 0.0 < modulation_duty < 1.0 &&
        source_resistance > 0.0 && load_resistance > 0.0 || throw(ArgumentError(
            "average four-quadrant equilibrium requires positive voltage and resistances with 0<d<1",
        ))
    average_polarity = 2.0 * modulation_duty - 1.0
    return FullBridgeConverterInitialState(
        average_polarity * input_voltage / (source_resistance + load_resistance),
    )
end

struct AverageFullBridgeConverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::FullBridgeConverterInitialState
    start_time_s::Float64
    stop_time_s::Float64

    function AverageFullBridgeConverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::FullBridgeConverterInitialState,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family in (FourQuadrantChopper, SinglePhaseTwoLevelBridge) ||
            throw(ArgumentError(
            "average full-bridge study requires four-quadrant or single-phase two-level conversion",
        ))
        selection.fidelity === AverageValue || throw(ArgumentError(
            "average four-quadrant study requires explicit AverageValue fidelity",
        ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "average four-quadrant study is a standalone conversion owner",
        ))
        topology.family === :full_bridge || throw(ArgumentError(
            "average four-quadrant study requires the canonical B200 full bridge",
        ))
        specification.topology_signatures == (bridge_topology_signature(topology),) ||
            throw(ArgumentError(
                "average four-quadrant specification does not bind its exact topology",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("average four-quadrant specification is not ready"))
        specification.modulation.kind === CarrierSinusoidalPulseWidthModulation ||
            throw(ArgumentError("average four-quadrant study requires bipolar PWM averaging"))
        selection.family === SinglePhaseTwoLevelBridge &&
            !(0.0 < specification.modulation.modulation_index <= 1.0) &&
            throw(ArgumentError(
                "average single-phase two-level modulation index must lie in (0, 1]",
            ))
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) || throw(ArgumentError(
            "average four-quadrant study parameters must be finite",
        ))
        all(>(0.0), values[1:4]) || throw(ArgumentError(
            "average four-quadrant voltage, resistances, and inductance must be positive",
        ))
        values[5] >= 0.0 && values[6] > values[5] || throw(ArgumentError(
            "average four-quadrant stop time must follow its nonnegative start time",
        ))
        step_count = (values[6] - values[5]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "average four-quadrant horizon must contain an integer number of fixed steps",
            ))
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:4]...,
            initial_state,
            values[5:6]...,
        )
    end
end

struct AverageFullBridgeConverterTrace
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    output_voltage_v::Vector{Float64}
    load_current_a::Vector{Float64}
    stored_energy_j::Vector{Float64}
    source_dissipated_energy_j::Vector{Float64}
    circuit_residual_v::Vector{Float64}
    energy_residual_w::Vector{Float64}
    result::ConverterSystemResult
end

struct SwitchingFullBridgeConverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::FullBridgeConverterInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingFullBridgeConverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::FullBridgeConverterInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family in (FourQuadrantChopper, SinglePhaseTwoLevelBridge) ||
            throw(ArgumentError(
            "switching full-bridge study requires four-quadrant or single-phase two-level conversion",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
            "four-quadrant switching execution requires switching-state or switching-detailed fidelity",
        ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "four-quadrant converter study is a standalone conversion owner",
        ))
        topology.family === :full_bridge || throw(ArgumentError(
            "four-quadrant converter requires the canonical B200 full bridge",
        ))
        specification.topology_signatures == (bridge_topology_signature(topology),) ||
            throw(ArgumentError(
                "four-quadrant specification does not bind its exact full-bridge topology",
            ))
        if selection.fidelity === SwitchingDetailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed four-quadrant execution requires typed semiconductor parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                    "four-quadrant specification does not bind its exact D200 device parameters",
                ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state four-quadrant execution cannot claim detailed semiconductor parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state four-quadrant execution cannot bind detailed device signatures",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("four-quadrant converter specification is not ready"))
        specification.modulation.kind === CarrierSinusoidalPulseWidthModulation ||
            throw(ArgumentError("four-quadrant converter requires bipolar carrier PWM"))
        selection.family === SinglePhaseTwoLevelBridge &&
            !(0.0 < specification.modulation.modulation_index <= 1.0) &&
            throw(ArgumentError(
                "single-phase two-level modulation index must lie in (0, 1]",
            ))
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) || throw(ArgumentError(
            "four-quadrant converter parameters must be finite",
        ))
        all(>(0.0), values[1:4]) || throw(ArgumentError(
            "four-quadrant voltage, source/load resistance, and inductance must be positive",
        ))
        values[5] >= 0.0 && values[6] > values[5] || throw(ArgumentError(
            "four-quadrant stop time must follow its nonnegative start time",
        ))
        timing = specification.timing
        if selection.fidelity === SwitchingDetailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "four-quadrant timestep must resolve recovered-charge and tail-current state",
            ))
        end
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "four-quadrant commutation dead time must span at least one fixed step",
        ))
        dead_time_steps = timing.dead_time_s / timing.fixed_step_s
        isapprox(dead_time_steps, round(dead_time_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "four-quadrant dead time must lie on the fixed-step calendar",
            ))
        step_count = (values[6] - values[5]) / timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "four-quadrant horizon must contain an integer number of fixed steps",
            ))
        carrier_steps = inv(timing.carrier_frequency_hz * timing.fixed_step_s)
        isapprox(carrier_steps, round(carrier_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "four-quadrant carrier period must contain an integer number of fixed steps",
            ))
        if selection.family === FourQuadrantChopper
            rising_edge_step = carrier_steps * specification.modulation.duty / 2.0
            falling_edge_step = carrier_steps *
                (1.0 - specification.modulation.duty / 2.0)
            all(edge -> isapprox(edge, round(edge); atol=1.0e-10, rtol=1.0e-10),
                (rising_edge_step, falling_edge_step)) || throw(ArgumentError(
                    "four-quadrant PWM edges must lie on the fixed-step calendar",
                ))
        end
        if selection.fidelity === SwitchingDetailed
            detailed_semiconductor.switching_energy_current_domain_maximum_a >=
                specification.rated_bases.current_a || throw(ArgumentError(
                    "four-quadrant switching-energy current domain must include rated current",
                ))
            detailed_semiconductor.switching_energy_blocking_voltage_domain_maximum_v >=
                values[1] || throw(ArgumentError(
                    "four-quadrant switching-energy voltage domain must include DC blocking voltage",
                ))
        end
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:4]...,
            initial_state,
            detailed_semiconductor,
            values[5:6]...,
        )
    end
end

struct SwitchingFullBridgeConverterTrace
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    output_voltage_v::Vector{Float64}
    load_current_a::Vector{Float64}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    controlled_conducting_state::BitMatrix
    diode_conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Matrix{Float64}
    diode_junction_temperature_k::Matrix{Float64}
    diode_recovered_charge_c::Matrix{Float64}
    controlled_turn_off_tail_current_a::Matrix{Float64}
    result::ConverterSystemResult
end

const FourQuadrantConverterInitialState = FullBridgeConverterInitialState
const SinglePhaseTwoLevelInverterInitialState = FullBridgeConverterInitialState
const AverageFourQuadrantConverterStudy = AverageFullBridgeConverterStudy
const AverageFourQuadrantConverterTrace = AverageFullBridgeConverterTrace
const AverageSinglePhaseTwoLevelInverterStudy = AverageFullBridgeConverterStudy
const AverageSinglePhaseTwoLevelInverterTrace = AverageFullBridgeConverterTrace
const SwitchingFourQuadrantConverterStudy = SwitchingFullBridgeConverterStudy
const SwitchingFourQuadrantConverterTrace = SwitchingFullBridgeConverterTrace
const SwitchingSinglePhaseTwoLevelInverterStudy = SwitchingFullBridgeConverterStudy
const SwitchingSinglePhaseTwoLevelInverterTrace = SwitchingFullBridgeConverterTrace

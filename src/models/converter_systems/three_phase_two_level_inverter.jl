export ThreePhaseTwoLevelInverterInitialState,
       AverageThreePhaseTwoLevelInverterStudy,
       AverageThreePhaseTwoLevelInverterTrace,
       SwitchingThreePhaseTwoLevelInverterStudy,
       SwitchingThreePhaseTwoLevelInverterTrace

struct ThreePhaseTwoLevelInverterInitialState
    phase_current_a::NTuple{3,Float64}

    function ThreePhaseTwoLevelInverterInitialState(phase_current_a=(0.0, 0.0, 0.0))
        length(phase_current_a) == 3 || throw(ArgumentError(
            "three-phase inverter initial state requires three finite phase currents",
        ))
        current = ntuple(index -> Float64(phase_current_a[index]), 3)
        all(isfinite, current) || throw(ArgumentError(
            "three-phase inverter initial state requires three finite phase currents",
        ))
        abs(sum(current)) <= 1.0e-10 * max(maximum(abs, current), 1.0) ||
            throw(ArgumentError(
                "three-wire inverter initial phase currents must sum to zero",
            ))
        return new(current)
    end
end

function _validate_three_phase_two_level_inverter(
    specification,
    topology,
    expected_fidelities,
)
    selection = specification.selection
    selection.family === ThreePhaseTwoLevelBridge || throw(ArgumentError(
        "three-phase inverter study requires the canonical three-phase two-level family",
    ))
    selection.fidelity in expected_fidelities || throw(ArgumentError(
        "three-phase inverter study received an incompatible fidelity",
    ))
    selection.application === StandaloneConversion || throw(ArgumentError(
        "three-phase inverter study is a standalone conversion owner",
    ))
    topology.family === :two_level_bridge && length(topology.valve_positions) == 6 ||
        throw(ArgumentError(
            "three-phase inverter requires the canonical six-position B200 two-level bridge",
        ))
    specification.topology_signatures == (bridge_topology_signature(topology),) ||
        throw(ArgumentError(
            "three-phase inverter specification does not bind its exact topology",
        ))
    specification.modulation.kind in (
        CarrierSinusoidalPulseWidthModulation,
        SpaceVectorPulseWidthModulation,
    ) || throw(ArgumentError(
        "three-phase inverter requires sinusoidal or space-vector PWM",
    ))
    0.0 < specification.modulation.modulation_index <= 1.0 ||
        throw(ArgumentError("three-phase inverter modulation index must lie in (0, 1]"))
    converter_system_is_ready(converter_system_readiness(specification)) ||
        throw(ArgumentError("three-phase inverter specification is not ready"))
    return nothing
end

struct AverageThreePhaseTwoLevelInverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::ThreePhaseTwoLevelInverterInitialState
    start_time_s::Float64
    stop_time_s::Float64

    function AverageThreePhaseTwoLevelInverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::ThreePhaseTwoLevelInverterInitialState,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        _validate_three_phase_two_level_inverter(
            specification,
            topology,
            (AverageValue,),
        )
        values = Float64.((
            input_voltage_v,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:3]) || throw(ArgumentError(
            "average three-phase inverter voltage, resistance, and inductance must be finite and positive",
        ))
        values[4] >= 0.0 && values[5] > values[4] || throw(ArgumentError(
            "average three-phase inverter stop time must follow its nonnegative start time",
        ))
        step_count = (values[5] - values[4]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "average three-phase inverter horizon must contain integer fixed steps",
            ))
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:3]...,
            initial_state,
            values[4:5]...,
        )
    end
end

struct AverageThreePhaseTwoLevelInverterTrace
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    stored_energy_j::Vector{Float64}
    circuit_residual_v::Matrix{Float64}
    energy_residual_w::Vector{Float64}
    result::ConverterSystemResult
end

struct SwitchingThreePhaseTwoLevelInverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::ThreePhaseTwoLevelInverterInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingThreePhaseTwoLevelInverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::ThreePhaseTwoLevelInverterInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        _validate_three_phase_two_level_inverter(
            specification,
            topology,
            (SwitchingStateEquivalent, SwitchingDetailed),
        )
        detailed = specification.selection.fidelity === SwitchingDetailed
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed three-phase inverter requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                    "three-phase inverter specification does not bind its D200 parameters",
                ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state three-phase inverter cannot claim D200 parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state three-phase inverter cannot bind device-fidelity signatures",
            ))
        end
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:4]) || throw(ArgumentError(
            "switching three-phase inverter electrical parameters must be finite and positive",
        ))
        values[5] >= 0.0 && values[6] > values[5] || throw(ArgumentError(
            "switching three-phase inverter stop time must follow its nonnegative start time",
        ))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "three-phase inverter dead time must span at least one fixed step",
        ))
        all(value -> isapprox(value, round(value); atol=1.0e-10, rtol=1.0e-10), (
            (values[6] - values[5]) / timing.fixed_step_s,
            inv(timing.carrier_frequency_hz * timing.fixed_step_s),
            timing.dead_time_s / timing.fixed_step_s,
        )) || throw(ArgumentError(
            "three-phase inverter horizon, carrier, and dead time must lie on the fixed-step calendar",
        ))
        if detailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "three-phase inverter timestep must resolve recovery and tail state",
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

struct SwitchingThreePhaseTwoLevelInverterTrace
    time_s::Vector{Float64}
    dc_link_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
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

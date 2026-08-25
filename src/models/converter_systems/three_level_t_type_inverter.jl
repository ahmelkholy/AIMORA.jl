export ThreeLevelTTypeInitialState,
       SwitchingThreeLevelTTypeStudy,
       SwitchingThreeLevelTTypeTrace

struct ThreeLevelTTypeInitialState
    phase_current_a::NTuple{3,Float64}
    upper_dc_link_voltage_v::Float64
    lower_dc_link_voltage_v::Float64

    function ThreeLevelTTypeInitialState(
        phase_current_a=(0.0, 0.0, 0.0);
        upper_dc_link_voltage_v::Real,
        lower_dc_link_voltage_v::Real,
    )
        length(phase_current_a) == 3 || throw(ArgumentError(
            "T-type initial state requires three phase currents",
        ))
        current = ntuple(index -> Float64(phase_current_a[index]), 3)
        upper = Float64(upper_dc_link_voltage_v)
        lower = Float64(lower_dc_link_voltage_v)
        all(isfinite, (current..., upper, lower)) && upper > 0.0 && lower > 0.0 ||
            throw(ArgumentError(
                "T-type initial currents must be finite and split-link voltages finite and positive",
            ))
        abs(sum(current)) <= 1.0e-10 * max(maximum(abs, current), 1.0) ||
            throw(ArgumentError("three-wire T-type initial phase currents must sum to zero"))
        return new(current, upper, lower)
    end
end

struct SwitchingThreeLevelTTypeStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topologies::NTuple{3,T}
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    dc_link_capacitance_f::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::ThreeLevelTTypeInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingThreeLevelTTypeStudy(
        specification::ConverterSystemSpecification;
        topologies::NTuple{3,<:BridgeTopologyDescriptor},
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        dc_link_capacitance_f::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::ThreeLevelTTypeInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === ThreeLevelTTypeBridge || throw(ArgumentError(
            "T-type study requires the canonical three-level T-type family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError("T-type execution requires switching-state or switching-detailed fidelity"))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "T-type study is a standalone conversion owner",
        ))
        all(topology -> topology.family === :t_type_leg &&
            length(topology.valve_positions) == 4, topologies) || throw(ArgumentError(
                "T-type study requires three canonical four-position B200 T-type legs",
            ))
        specification.topology_signatures == Tuple(bridge_topology_signature.(topologies)) ||
            throw(ArgumentError("T-type specification does not bind its exact three leg topologies"))
        specification.modulation.kind in (
            CarrierSinusoidalPulseWidthModulation,
            PhaseShiftedCarrierPulseWidthModulation,
        ) || throw(ArgumentError("T-type execution requires level-shifted carrier PWM"))
        0.0 < specification.modulation.modulation_index <= 1.0 || throw(ArgumentError(
            "T-type modulation index must lie in (0, 1]",
        ))
        detailed = selection.fidelity === SwitchingDetailed
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed T-type execution requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError("T-type specification does not bind its D200 parameters"))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state T-type execution cannot claim D200 parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state T-type execution cannot bind device-fidelity signatures",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("T-type specification is not ready"))
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            dc_link_capacitance_f,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:5]) || throw(ArgumentError(
            "T-type source, split-link, and load parameters must be finite and positive",
        ))
        values[6] >= 0.0 && values[7] > values[6] || throw(ArgumentError(
            "T-type stop time must follow its nonnegative start time",
        ))
        isapprox(
            initial_state.upper_dc_link_voltage_v + initial_state.lower_dc_link_voltage_v,
            values[1];
            atol=1.0e-10 * max(values[1], 1.0),
            rtol=1.0e-10,
        ) || throw(ArgumentError("T-type initial split-link voltages must sum to the source voltage"))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "T-type dead time must span at least one fixed step",
        ))
        all(value -> isapprox(value, round(value); atol=1.0e-10, rtol=1.0e-10), (
            (values[7] - values[6]) / timing.fixed_step_s,
            inv(timing.carrier_frequency_hz * timing.fixed_step_s),
            timing.dead_time_s / timing.fixed_step_s,
        )) || throw(ArgumentError(
            "T-type horizon, carrier, and dead time must lie on the fixed-step calendar",
        ))
        if detailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "T-type timestep must resolve recovery and tail state",
            ))
        end
        return new{typeof(topologies[1])}(
            specification,
            topologies,
            values[1:5]...,
            initial_state,
            detailed_semiconductor,
            values[6:7]...,
        )
    end
end

struct SwitchingThreeLevelTTypeTrace
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    upper_dc_link_voltage_v::Vector{Float64}
    lower_dc_link_voltage_v::Vector{Float64}
    midpoint_current_a::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    requested_level::Matrix{Int8}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    controlled_conducting_state::BitMatrix
    antiparallel_diode_conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Matrix{Float64}
    antiparallel_junction_temperature_k::Matrix{Float64}
    antiparallel_recovered_charge_c::Matrix{Float64}
    controlled_turn_off_tail_current_a::Matrix{Float64}
    result::ConverterSystemResult
end

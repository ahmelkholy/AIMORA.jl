export FlyingCapacitorInitialState,
       SwitchingFlyingCapacitorStudy,
       SwitchingFlyingCapacitorTrace

struct FlyingCapacitorInitialState
    phase_current_a::NTuple{3,Float64}
    flying_capacitor_voltage_v::NTuple{3,Float64}

    function FlyingCapacitorInitialState(
        phase_current_a=(0.0, 0.0, 0.0);
        flying_capacitor_voltage_v,
    )
        length(phase_current_a) == 3 && length(flying_capacitor_voltage_v) == 3 ||
            throw(DimensionMismatch(
                "flying-capacitor initial state requires three currents and three voltages",
            ))
        current = ntuple(index -> Float64(phase_current_a[index]), 3)
        voltage = ntuple(index -> Float64(flying_capacitor_voltage_v[index]), 3)
        all(isfinite, (current..., voltage...)) && all(>(0.0), voltage) ||
            throw(ArgumentError(
                "flying-capacitor initial currents must be finite and voltages positive",
            ))
        abs(sum(current)) <= 1.0e-10 * max(maximum(abs, current), 1.0) ||
            throw(ArgumentError(
                "three-wire flying-capacitor initial phase currents must sum to zero",
            ))
        return new(current, voltage)
    end
end

struct SwitchingFlyingCapacitorStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topologies::NTuple{3,T}
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    flying_capacitance_f::Float64
    balance_voltage_tolerance_v::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::FlyingCapacitorInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingFlyingCapacitorStudy(
        specification::ConverterSystemSpecification;
        topologies::NTuple{3,<:BridgeTopologyDescriptor},
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        flying_capacitance_f::Real,
        balance_voltage_tolerance_v::Real=1.0e-3,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::FlyingCapacitorInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === FlyingCapacitorBridge || throw(ArgumentError(
            "flying-capacitor study requires the canonical flying-capacitor family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "flying-capacitor execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "flying-capacitor study is a standalone conversion owner",
        ))
        all(topology -> topology.family === :flying_capacitor_leg &&
            length(topology.valve_positions) == 4 &&
            length(topology.passive_positions) == 1 &&
            only(topology.passive_positions).kind === :capacitor, topologies) ||
            throw(ArgumentError(
                "flying-capacitor study requires three canonical four-position capacitive legs",
            ))
        specification.topology_signatures == Tuple(bridge_topology_signature.(topologies)) ||
            throw(ArgumentError(
                "flying-capacitor specification does not bind its exact leg topologies",
            ))
        specification.modulation.kind in (
            CarrierSinusoidalPulseWidthModulation,
            PhaseShiftedCarrierPulseWidthModulation,
        ) || throw(ArgumentError(
            "flying-capacitor execution requires a carrier modulation",
        ))
        0.0 < specification.modulation.modulation_index <= 1.0 || throw(ArgumentError(
            "flying-capacitor modulation index must lie in (0, 1]",
        ))
        detailed = selection.fidelity === SwitchingDetailed
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed flying-capacitor execution requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                    "flying-capacitor specification does not bind its D200 parameters",
                ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state flying-capacitor execution cannot accept D200 parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state flying-capacitor execution cannot bind D200 signatures",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("flying-capacitor specification is not ready"))
        values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            flying_capacitance_f,
            balance_voltage_tolerance_v,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:3]) && values[4] >= 0.0 &&
            all(>(0.0), values[5:6]) || throw(ArgumentError(
            "flying-capacitor source, capacitor, load, and balance tolerance parameters are invalid",
        ))
        values[7] >= 0.0 && values[8] > values[7] || throw(ArgumentError(
            "flying-capacitor stop time must follow its nonnegative start time",
        ))
        all(voltage -> isapprox(
            voltage,
            0.5 * values[1];
            atol=1.0e-10 * max(values[1], 1.0),
            rtol=1.0e-10,
        ), initial_state.flying_capacitor_voltage_v) || throw(ArgumentError(
            "flying-capacitor initial voltages must equal half the DC source voltage",
        ))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "flying-capacitor dead time must span at least one fixed step",
        ))
        all(value -> isapprox(value, round(value); atol=1.0e-10, rtol=1.0e-10), (
            (values[8] - values[7]) / timing.fixed_step_s,
            inv(timing.carrier_frequency_hz * timing.fixed_step_s),
            timing.dead_time_s / timing.fixed_step_s,
        )) || throw(ArgumentError(
            "flying-capacitor horizon, carrier, and dead time must lie on the fixed-step calendar",
        ))
        if detailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "flying-capacitor timestep must resolve recovery and tail state",
            ))
        end
        return new{typeof(topologies[1])}(
            specification,
            topologies,
            values[1:6]...,
            initial_state,
            detailed_semiconductor,
            values[7:8]...,
        )
    end
end

struct SwitchingFlyingCapacitorTrace
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    flying_capacitor_voltage_v::Matrix{Float64}
    flying_capacitor_current_a::Matrix{Float64}
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

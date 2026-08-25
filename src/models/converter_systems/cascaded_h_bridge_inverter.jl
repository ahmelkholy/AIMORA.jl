export CascadedHBridgeInitialState,
       SwitchingCascadedHBridgeStudy,
       SwitchingCascadedHBridgeTrace

struct CascadedHBridgeInitialState
    phase_current_a::NTuple{3,Float64}
    cell_dc_voltage_v::Matrix{Float64}

    function CascadedHBridgeInitialState(
        phase_current_a=(0.0, 0.0, 0.0);
        cell_dc_voltage_v,
    )
        length(phase_current_a) == 3 || throw(DimensionMismatch(
            "cascaded H-bridge initial state requires three phase currents",
        ))
        current = ntuple(index -> Float64(phase_current_a[index]), 3)
        voltage = Matrix{Float64}(cell_dc_voltage_v)
        length(current) == 3 && size(voltage, 1) == 3 && 2 <= size(voltage, 2) <= 8 ||
            throw(DimensionMismatch(
                "cascaded H-bridge initial state requires three currents and a three-by-two-through-eight cell-voltage matrix",
            ))
        all(isfinite, (current..., voltage...)) && all(>(0.0), voltage) ||
            throw(ArgumentError(
                "cascaded H-bridge initial currents must be finite and cell voltages positive",
            ))
        abs(sum(current)) <= 1.0e-10 * max(maximum(abs, current), 1.0) ||
            throw(ArgumentError(
                "three-wire cascaded H-bridge initial phase currents must sum to zero",
            ))
        return new(current, copy(voltage))
    end
end

function _cascaded_h_bridge_parameter_matrix(values, cell_count, name)
    matrix = if values isa Real
        fill(Float64(values), 3, cell_count)
    elseif values isa AbstractVector || values isa Tuple
        vector = Float64.(values)
        length(vector) == cell_count || throw(DimensionMismatch(
            "cascaded H-bridge $name vector must contain one value per cell",
        ))
        repeat(reshape(vector, 1, :), 3, 1)
    else
        Matrix{Float64}(values)
    end
    size(matrix) == (3, cell_count) || throw(DimensionMismatch(
        "cascaded H-bridge $name must be scalar, cell-count vector, or three-by-cell-count matrix",
    ))
    all(isfinite, matrix) && all(>(0.0), matrix) || throw(ArgumentError(
        "cascaded H-bridge $name values must be finite and positive",
    ))
    return matrix
end

struct SwitchingCascadedHBridgeStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topologies::NTuple{3,T}
    cell_dc_capacitance_f::Matrix{Float64}
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::CascadedHBridgeInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingCascadedHBridgeStudy(
        specification::ConverterSystemSpecification;
        topologies::NTuple{3,<:BridgeTopologyDescriptor},
        cell_dc_capacitance_f,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::CascadedHBridgeInitialState,
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === CascadedHBridge || throw(ArgumentError(
            "cascaded H-bridge study requires the canonical cascaded family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "cascaded H-bridge execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "cascaded H-bridge study is a standalone conversion owner",
        ))
        cell_count = selection.cell_count
        all(topology -> topology.family === :cascaded_h_bridge_phase &&
            length(topology.valve_positions) == 4 * cell_count &&
            length(topology.state_groups) == cell_count &&
            isempty(topology.passive_positions), topologies) || throw(ArgumentError(
            "cascaded H-bridge study requires three canonical equal-cell phase topologies",
        ))
        specification.topology_signatures ==
            Tuple(bridge_topology_signature.(topologies)) || throw(ArgumentError(
            "cascaded H-bridge specification does not bind its exact phase topologies",
        ))
        specification.modulation.kind in (
            CarrierSinusoidalPulseWidthModulation,
            PhaseShiftedCarrierPulseWidthModulation,
            SelectiveHarmonicElimination,
            NearestLevelModulation,
        ) || throw(ArgumentError(
            "cascaded H-bridge execution requires carrier, selective-harmonic, or nearest-level modulation",
        ))
        0.0 < specification.modulation.modulation_index <= 1.0 ||
            throw(ArgumentError(
                "cascaded H-bridge modulation index must lie in (0, 1]",
            ))
        if specification.modulation.kind === SelectiveHarmonicElimination
            angles = specification.modulation.selective_harmonic_angles_rad
            1 <= length(angles) <= cell_count && issorted(angles) &&
                all(angle -> 0.0 < angle < pi / 2.0, angles) || throw(ArgumentError(
                "cascaded H-bridge selective-harmonic angles must be ordered inside (0, pi/2)",
            ))
        elseif !isempty(specification.modulation.selective_harmonic_angles_rad)
            throw(ArgumentError(
                "cascaded H-bridge selective-harmonic angles require selective-harmonic modulation",
            ))
        end
        detailed = selection.fidelity === SwitchingDetailed
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed cascaded H-bridge execution requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                "cascaded H-bridge specification does not bind its complete D200 parameters",
            ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state cascaded H-bridge execution cannot accept D200 parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state cascaded H-bridge execution cannot bind D200 signatures",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("cascaded H-bridge specification is not ready"))
        capacitance = _cascaded_h_bridge_parameter_matrix(
            cell_dc_capacitance_f,
            cell_count,
            "cell capacitance",
        )
        size(initial_state.cell_dc_voltage_v) == (3, cell_count) ||
            throw(DimensionMismatch(
                "cascaded H-bridge initial cell-voltage matrix does not match its selected cell count",
            ))
        values = Float64.((
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && values[1] > 0.0 && values[2] > 0.0 &&
            values[3] >= 0.0 && values[4] > values[3] || throw(ArgumentError(
            "cascaded H-bridge load and horizon domain is invalid",
        ))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "cascaded H-bridge dead time must span at least one fixed step",
        ))
        all(value -> isapprox(value, round(value); atol=1.0e-10, rtol=1.0e-10), (
            (values[4] - values[3]) / timing.fixed_step_s,
            inv(timing.carrier_frequency_hz * timing.fixed_step_s),
            timing.dead_time_s / timing.fixed_step_s,
        )) || throw(ArgumentError(
            "cascaded H-bridge horizon, carrier, and dead time must lie on the fixed-step calendar",
        ))
        if detailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "cascaded H-bridge timestep must resolve recovery and tail state",
            ))
        end
        return new{typeof(topologies[1])}(
            specification,
            topologies,
            capacitance,
            values[1],
            values[2],
            initial_state,
            detailed_semiconductor,
            values[3],
            values[4],
        )
    end
end

struct SwitchingCascadedHBridgeTrace
    time_s::Vector{Float64}
    cell_dc_voltage_v::Array{Float64,3}
    cell_dc_current_a::Array{Float64,3}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    requested_level::Matrix{Int8}
    requested_cell_state::Array{Int8,3}
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

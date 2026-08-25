export MatrixConverterInitialState,
       SwitchingMatrixConverterStudy,
       SwitchingMatrixConverterTrace

struct MatrixConverterInitialState
    output_phase_current_a::NTuple{3,Float64}
    input_for_output::NTuple{3,Int}

    function MatrixConverterInitialState(
        output_phase_current_a=(0.0, 0.0, 0.0);
        input_for_output=(1, 2, 3),
    )
        length(output_phase_current_a) == 3 || throw(DimensionMismatch(
            "matrix-converter initial state requires three output currents",
        ))
        current = ntuple(index -> Float64(output_phase_current_a[index]), 3)
        connection = ntuple(index -> Int(input_for_output[index]), 3)
        all(isfinite, current) || throw(ArgumentError(
            "matrix-converter initial currents must be finite",
        ))
        abs(sum(current)) <= 1.0e-10 * max(maximum(abs, current), 1.0) ||
            throw(ArgumentError(
                "three-wire matrix-converter initial output currents must sum to zero",
            ))
        all(value -> value in 1:3, connection) || throw(ArgumentError(
            "each matrix-converter output must initially select one input in 1:3",
        ))
        return new(current, connection)
    end
end

struct SwitchingMatrixConverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_phase_voltage_peak_v::Float64
    input_frequency_hz::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::MatrixConverterInitialState
    detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingMatrixConverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_phase_voltage_peak_v::Real,
        input_frequency_hz::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::MatrixConverterInitialState=MatrixConverterInitialState(),
        detailed_semiconductor::Union{Nothing,DetailedChopperSemiconductorParameters}=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === ThreePhaseMatrixConverter || throw(ArgumentError(
            "matrix-converter study requires the canonical three-phase matrix family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "matrix-converter execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "matrix-converter study is a standalone conversion owner",
        ))
        topology.family === :matrix_converter &&
            length(topology.valve_positions) == 18 &&
            length(topology.state_groups) == 1 &&
            size(only(topology.state_groups).admitted_states) == (18, 27) &&
            isempty(topology.passive_positions) || throw(ArgumentError(
                "matrix-converter study requires the canonical 18-position 3x3 topology",
            ))
        specification.topology_signatures == (bridge_topology_signature(topology),) ||
            throw(ArgumentError(
                "matrix-converter specification does not bind its exact topology",
            ))
        specification.modulation.kind === MatrixSpaceVectorModulation ||
            throw(ArgumentError(
                "matrix-converter execution requires matrix space-vector modulation",
            ))
        0.0 < specification.modulation.modulation_index <= sqrt(3.0) / 2.0 ||
            throw(ArgumentError(
                "matrix-converter modulation index must lie in (0, sqrt(3)/2]",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("matrix-converter specification is not ready"))
        detailed = selection.fidelity === SwitchingDetailed
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed matrix conversion requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_chopper_semiconductor_signatures(detailed_semiconductor) ||
                throw(ArgumentError(
                    "matrix-converter specification does not bind its complete D200 parameters",
                ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state matrix conversion cannot accept D200 parameters",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state matrix conversion cannot bind D200 signatures",
            ))
        end
        values = Float64.((
            input_phase_voltage_peak_v,
            input_frequency_hz,
            source_resistance_ohm,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:5]) &&
            values[6] >= 0.0 && values[7] > values[6] || throw(ArgumentError(
                "matrix-converter source, load, and horizon domain is invalid",
            ))
        values[1] <= 1.0e6 && 1.0 <= values[2] <= 1.0e3 || throw(ArgumentError(
            "matrix-converter source voltage or frequency is outside the released domain",
        ))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "matrix safe commutation requires at least one fixed step per stage",
        ))
        timing.minimum_pulse_s <= timing.fixed_step_s || throw(ArgumentError(
            "matrix safe-commutation stages require minimum pulse no longer than one fixed step",
        ))
        all(value -> isapprox(value, round(value); atol=1.0e-10, rtol=1.0e-10), (
            (values[7] - values[6]) / timing.fixed_step_s,
            inv(timing.carrier_frequency_hz * timing.fixed_step_s),
            timing.dead_time_s / timing.fixed_step_s,
        )) || throw(ArgumentError(
            "matrix-converter horizon, carrier, and commutation time must lie on the fixed-step calendar",
        ))
        inv(timing.carrier_frequency_hz * timing.fixed_step_s) >= 40.0 ||
            throw(ArgumentError(
                "matrix-converter carrier requires at least forty fixed steps per period",
            ))
        if detailed
            timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "matrix-converter timestep must resolve recovery and tail state",
            ))
        end
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:5]...,
            initial_state,
            detailed_semiconductor,
            values[6:7]...,
        )
    end
end

struct SwitchingMatrixConverterTrace
    time_s::Vector{Float64}
    input_source_voltage_v::Matrix{Float64}
    input_terminal_voltage_v::Matrix{Float64}
    input_current_a::Matrix{Float64}
    output_reference_voltage_v::Matrix{Float64}
    output_phase_voltage_v::Matrix{Float64}
    output_phase_current_a::Matrix{Float64}
    requested_connection::Array{Bool,3}
    applied_connection::Array{Bool,3}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    conducting_state::BitMatrix
    commutation_stage::Matrix{UInt8}
    commutation_direction::Matrix{Int8}
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    junction_temperature_k::Matrix{Float64}
    recovered_charge_c::Matrix{Float64}
    turn_off_tail_current_a::Matrix{Float64}
    result::ConverterSystemResult
end

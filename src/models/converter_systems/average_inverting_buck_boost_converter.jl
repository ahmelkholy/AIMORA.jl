export AverageInvertingBuckBoostConverterInitialState,
       AverageInvertingBuckBoostConverterStudy,
       AverageInvertingBuckBoostConverterTrace

const AverageInvertingBuckBoostConverterInitialState = DCDCConverterInitialState

struct AverageInvertingBuckBoostConverterStudy{T<:BridgeTopologyDescriptor}
    specification::ConverterSystemSpecification
    topology::T
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    inductor_resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
    load_resistance_ohm::Float64
    initial_state::DCDCConverterInitialState
    start_time_s::Float64
    stop_time_s::Float64

    function AverageInvertingBuckBoostConverterStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        inductor_resistance_ohm::Real,
        inductance_h::Real,
        capacitance_f::Real,
        load_resistance_ohm::Real,
        initial_state::DCDCConverterInitialState,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === InvertingBuckBoostChopper || throw(ArgumentError(
            "average inverting buck-boost study requires the canonical inverting buck-boost family",
        ))
        selection.fidelity === AverageValue || throw(ArgumentError(
            "average inverting buck-boost study requires explicit AverageValue fidelity",
        ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "average inverting buck-boost study is a standalone converter",
        ))
        topology.family === :inverting_buck_boost || throw(ArgumentError(
            "average inverting buck-boost study requires its canonical B200 topology",
        ))
        specification.topology_signatures == (bridge_topology_signature(topology),) ||
            throw(ArgumentError(
                "average inverting buck-boost specification does not bind the supplied B200 topology",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError(
                "average inverting buck-boost converter specification is not ready",
            ))
        initial_state.inductor_current_a >= 0.0 && initial_state.output_voltage_v <= 0.0 ||
            throw(ArgumentError(
                "average inverting buck-boost initial current must be nonnegative and output voltage nonpositive",
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
        all(isfinite, values) || throw(ArgumentError(
            "average inverting buck-boost study parameters must be finite",
        ))
        values[1] > 0.0 && values[2] > 0.0 && values[3] >= 0.0 &&
            all(>(0.0), values[4:6]) || throw(ArgumentError(
                "average inverting buck-boost voltage, source resistance, L, C, and load must be positive and inductor resistance nonnegative",
            ))
        values[7] >= 0.0 && values[8] > values[7] || throw(ArgumentError(
            "average inverting buck-boost stop time must follow its nonnegative start time",
        ))
        step_count = (values[8] - values[7]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "average inverting buck-boost horizon must contain an integer number of fixed steps",
            ))
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:6]...,
            initial_state,
            values[7:8]...,
        )
    end
end

struct AverageInvertingBuckBoostConverterTrace
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    output_voltage_v::Vector{Float64}
    inductor_current_a::Vector{Float64}
    load_current_a::Vector{Float64}
    stored_energy_j::Vector{Float64}
    dissipated_energy_j::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    result::ConverterSystemResult
end

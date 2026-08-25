export DCDCConverterInitialState,
       AverageBuckConverterInitialState,
       AverageBuckConverterStudy,
       AverageBuckConverterTrace,
       average_buck_operating_point,
       boost_converter_operating_point,
       inverting_buck_boost_operating_point

struct DCDCConverterInitialState
    inductor_current_a::Float64
    output_voltage_v::Float64

    function DCDCConverterInitialState(inductor_current_a, output_voltage_v)
        values = Float64.((inductor_current_a, output_voltage_v))
        all(isfinite, values) || throw(ArgumentError("DC/DC initial state must be finite"))
        return new(values...)
    end
end

const AverageBuckConverterInitialState = DCDCConverterInitialState

struct AverageBuckConverterStudy{T<:BridgeTopologyDescriptor}
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

    function AverageBuckConverterStudy(
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
        selection.family === BuckChopper || throw(ArgumentError(
            "average buck study requires the canonical buck-chopper family",
        ))
        selection.fidelity === AverageValue || throw(ArgumentError(
            "average buck study requires explicit AverageValue fidelity",
        ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "average buck study is a standalone converter, not an application composition",
        ))
        topology.family === :step_down_chopper || throw(ArgumentError(
            "average buck study requires the canonical B200 step-down chopper topology",
        ))
        topology_signature = bridge_topology_signature(topology)
        specification.topology_signatures == (topology_signature,) || throw(ArgumentError(
            "average buck specification does not bind the exact supplied B200 topology identity",
        ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("average buck converter specification is not ready"))
        initial_state.inductor_current_a >= 0.0 && initial_state.output_voltage_v >= 0.0 ||
            throw(ArgumentError(
                "average buck initial current and output voltage must follow their nonnegative reference directions",
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
        all(isfinite, values) || throw(ArgumentError("average buck study parameters must be finite"))
        values[1] > 0.0 || throw(ArgumentError("average buck input voltage must be positive"))
        values[2] > 0.0 || throw(ArgumentError("average buck source resistance must be positive and finite"))
        values[3] >= 0.0 || throw(ArgumentError("average buck inductor resistance must be nonnegative"))
        all(>(0.0), values[4:6]) || throw(ArgumentError("average buck L, C, and load resistance must be positive"))
        values[7] >= 0.0 && values[8] > values[7] || throw(ArgumentError(
            "average buck stop time must be later than its nonnegative start time",
        ))
        step_count = (values[8] - values[7]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("average buck horizon must contain an integer number of fixed steps"))
        return new{typeof(topology)}(
            specification,
            topology,
            values[1:6]...,
            initial_state,
            values[7:8]...,
        )
    end
end

function average_buck_operating_point(
    input_voltage_v::Real,
    duty::Real,
    source_resistance_ohm::Real,
    inductor_resistance_ohm::Real,
    load_resistance_ohm::Real,
)
    input_voltage = Float64(input_voltage_v)
    d = Float64(duty)
    source_resistance = Float64(source_resistance_ohm)
    inductor_resistance = Float64(inductor_resistance_ohm)
    load_resistance = Float64(load_resistance_ohm)
    input_voltage > 0.0 && 0.0 < d < 1.0 && source_resistance > 0.0 &&
        inductor_resistance >= 0.0 && load_resistance > 0.0 || throw(ArgumentError(
            "average buck equilibrium requires positive voltage, source/load resistance, 0<d<1, and nonnegative inductor resistance",
        ))
    total_series_resistance = source_resistance + inductor_resistance
    output_voltage = d * input_voltage * load_resistance /
        (load_resistance + total_series_resistance)
    inductor_current = output_voltage / load_resistance
    return DCDCConverterInitialState(inductor_current, output_voltage)
end

function boost_converter_operating_point(
    input_voltage_v::Real,
    duty::Real,
    source_resistance_ohm::Real,
    inductor_resistance_ohm::Real,
    load_resistance_ohm::Real,
)
    input_voltage, d, source_resistance, inductor_resistance, load_resistance =
        Float64.((input_voltage_v, duty, source_resistance_ohm,
            inductor_resistance_ohm, load_resistance_ohm))
    input_voltage > 0.0 && 0.0 < d < 1.0 && source_resistance > 0.0 &&
        inductor_resistance >= 0.0 && load_resistance > 0.0 || throw(ArgumentError(
            "boost equilibrium requires positive voltage, source/load resistance, 0<d<1, and nonnegative inductor resistance",
        ))
    transfer_fraction = 1.0 - d
    series_resistance = source_resistance + inductor_resistance
    output_voltage = input_voltage / (
        transfer_fraction + series_resistance / (load_resistance * transfer_fraction)
    )
    inductor_current = output_voltage / (load_resistance * transfer_fraction)
    return DCDCConverterInitialState(inductor_current, output_voltage)
end

function inverting_buck_boost_operating_point(
    input_voltage_v::Real,
    duty::Real,
    source_resistance_ohm::Real,
    inductor_resistance_ohm::Real,
    load_resistance_ohm::Real,
)
    input_voltage, d, source_resistance, inductor_resistance, load_resistance =
        Float64.((input_voltage_v, duty, source_resistance_ohm,
            inductor_resistance_ohm, load_resistance_ohm))
    input_voltage > 0.0 && 0.0 < d < 1.0 && source_resistance > 0.0 &&
        inductor_resistance >= 0.0 && load_resistance > 0.0 || throw(ArgumentError(
            "inverting buck-boost equilibrium requires positive voltage, source/load resistance, 0<d<1, and nonnegative inductor resistance",
        ))
    transfer_fraction = 1.0 - d
    loop_resistance = inductor_resistance + d * source_resistance +
        load_resistance * transfer_fraction^2
    inductor_current = d * input_voltage / loop_resistance
    output_voltage = -load_resistance * transfer_fraction * inductor_current
    return DCDCConverterInitialState(inductor_current, output_voltage)
end

struct AverageBuckConverterTrace
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

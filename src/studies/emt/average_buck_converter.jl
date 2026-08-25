using ..ConverterSystems
using ..Branches
using ..Nodal

mutable struct AverageBuckConverterRuntime{S,N,B,C}
    study::S
    network::N
    inductor::B
    capacitor::C
    time_s::Float64
    dissipated_energy_j::Float64
    accepted_step_index::Int
    time_trace_s::Vector{Float64}
    input_voltage_trace_v::Vector{Float64}
    input_current_trace_a::Vector{Float64}
    output_voltage_trace_v::Vector{Float64}
    inductor_current_trace_a::Vector{Float64}
    load_current_trace_a::Vector{Float64}
    stored_energy_trace_j::Vector{Float64}
    dissipated_energy_trace_j::Vector{Float64}
    kcl_residual_trace_a::Vector{Float64}
    energy_residual_trace_w::Vector{Float64}
end

function _seed_average_buck_network(study::AverageBuckConverterStudy)
    duty = study.specification.modulation.duty
    controlled_voltage = duty * study.input_voltage_v
    initial = study.initial_state
    initial_source_drop = study.source_resistance_ohm * initial.inductor_current_a
    initial_inductor_voltage =
        controlled_voltage - initial_source_drop - initial.output_voltage_v
    initial_capacitor_current =
        initial.inductor_current_a - initial.output_voltage_v / study.load_resistance_ohm
    source = Branches.TwoTerminalTheveninSource(
        1,
        0,
        inv(study.source_resistance_ohm),
        _time_s -> controlled_voltage,
    )
    inductor = Branches.SeriesRLBranch(
        1,
        2,
        study.inductor_resistance_ohm,
        study.inductance_h,
        initial.inductor_current_a,
        initial_inductor_voltage,
        initial.inductor_current_a,
    )
    capacitor = Branches.CapacitorBranch(
        2,
        0,
        study.capacitance_f,
        initial_capacitor_current,
        initial.output_voltage_v,
        initial_capacitor_current,
    )
    load = Branches.ConductanceBranch(2, 0, inv(study.load_resistance_ohm))
    network = Nodal.NodalSystem(2, Any[source, inductor, capacitor, load])
    network.v[1] = initial.output_voltage_v + initial_inductor_voltage
    network.v[2] = initial.output_voltage_v
    return network, inductor, capacitor
end

function _average_buck_stored_energy(study, inductor_current, output_voltage)
    return 0.5 * study.inductance_h * inductor_current^2 +
        0.5 * study.capacitance_f * output_voltage^2
end

function prepare_average_buck_converter(study::AverageBuckConverterStudy)
    network, inductor, capacitor = _seed_average_buck_network(study)
    sample_count = Int(round(
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    )) + 1
    traces = ntuple(_ -> zeros(Float64, sample_count), 10)
    runtime = AverageBuckConverterRuntime(
        study,
        network,
        inductor,
        capacitor,
        study.start_time_s,
        0.0,
        0,
        traces...,
    )
    _record_average_buck_sample!(runtime, 1, 0.0)
    return runtime
end

function _record_average_buck_sample!(runtime, sample_index, energy_residual_w)
    study = runtime.study
    duty = study.specification.modulation.duty
    inductor_current = runtime.inductor.i_last
    output_voltage = runtime.network.v[2]
    source_node_voltage = runtime.network.v[1]
    load_current = output_voltage / study.load_resistance_ohm
    input_current = duty * inductor_current
    capacitor_current = runtime.capacitor.i_last
    kcl_residual = inductor_current - load_current - capacitor_current
    stored_energy = _average_buck_stored_energy(study, inductor_current, output_voltage)
    runtime.time_trace_s[sample_index] = runtime.time_s
    runtime.input_voltage_trace_v[sample_index] = study.input_voltage_v
    runtime.input_current_trace_a[sample_index] = input_current
    runtime.output_voltage_trace_v[sample_index] = output_voltage
    runtime.inductor_current_trace_a[sample_index] = inductor_current
    runtime.load_current_trace_a[sample_index] = load_current
    runtime.stored_energy_trace_j[sample_index] = stored_energy
    runtime.dissipated_energy_trace_j[sample_index] = runtime.dissipated_energy_j
    runtime.kcl_residual_trace_a[sample_index] = kcl_residual
    runtime.energy_residual_trace_w[sample_index] = energy_residual_w
    return source_node_voltage
end

function _advance_average_buck_converter!(runtime, sample_index)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    previous_energy = runtime.stored_energy_trace_j[sample_index - 1]
    endpoint = study.start_time_s + (sample_index - 1) * step
    Nodal.solve_algebraic_state!(runtime.network, endpoint, step)
    Nodal.accept_algebraic_state!(runtime.network, step)
    runtime.time_s = endpoint
    current = runtime.inductor.i_last
    output_voltage = runtime.network.v[2]
    load_power = output_voltage^2 / study.load_resistance_ohm
    loss_power = current^2 * (study.source_resistance_ohm + study.inductor_resistance_ohm)
    input_power = study.specification.modulation.duty * study.input_voltage_v * current
    stored_energy = _average_buck_stored_energy(study, current, output_voltage)
    energy_rate = (stored_energy - previous_energy) / step
    energy_residual = input_power - load_power - loss_power - energy_rate
    runtime.dissipated_energy_j += step * loss_power
    _record_average_buck_sample!(runtime, sample_index, energy_residual)
    return runtime
end

function _average_buck_state(runtime)
    study = runtime.study
    input_current = study.specification.modulation.duty * runtime.inductor.i_last
    output_voltage = runtime.network.v[2]
    stored_energy = _average_buck_stored_energy(study, runtime.inductor.i_last, output_voltage)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(runtime.inductor.i_last),
        repr(output_voltage),
        repr(stored_energy),
        repr(runtime.dissipated_energy_j),
    ), '\n')))
    return ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, output_voltage],
        [input_current, output_voltage / study.load_resistance_ohm],
        BitVector(),
        BitVector(),
        BitVector(),
        [study.capacitance_f * output_voltage],
        [study.inductance_h * runtime.inductor.i_last],
        Float64[],
        [study.specification.modulation.duty],
        stored_energy,
        runtime.dissipated_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
end

function execute_average_buck_converter!(runtime::AverageBuckConverterRuntime)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    state = _average_buck_state(runtime)
    maximum_kcl = maximum(abs, runtime.kcl_residual_trace_a; init=0.0)
    energy_scale = max(
        maximum(abs, runtime.input_voltage_trace_v .* runtime.input_current_trace_a; init=0.0),
        maximum(abs, runtime.output_voltage_trace_v .* runtime.load_current_trace_a; init=0.0),
        1.0,
    )
    relative_energy = maximum(abs, runtime.energy_residual_trace_w; init=0.0) / energy_scale
    result = ConverterSystems.converter_system_result(
        runtime.study.specification,
        state;
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=maximum_kcl,
        relative_charge_residual=maximum_kcl,
        relative_energy_residual=relative_energy,
    )
    return ConverterSystems.AverageBuckConverterTrace(
        runtime.time_trace_s,
        runtime.input_voltage_trace_v,
        runtime.input_current_trace_a,
        runtime.output_voltage_trace_v,
        runtime.inductor_current_trace_a,
        runtime.load_current_trace_a,
        runtime.stored_energy_trace_j,
        runtime.dissipated_energy_trace_j,
        runtime.kcl_residual_trace_a,
        runtime.energy_residual_trace_w,
        result,
    )
end

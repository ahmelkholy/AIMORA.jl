using ..ConverterSystems

mutable struct AverageBoostConverterRuntime{S}
    study::S
    inductor_current_a::Float64
    output_voltage_v::Float64
    time_s::Float64
    dissipated_energy_j::Float64
    accepted_step_index::Int
    transition_matrix::Matrix{Float64}
    source_increment::Vector{Float64}
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

function _average_boost_state_matrix(study)
    transfer_fraction = 1.0 - study.specification.modulation.duty
    series_resistance = study.source_resistance_ohm + study.inductor_resistance_ohm
    return [
        -series_resistance / study.inductance_h -transfer_fraction / study.inductance_h
        transfer_fraction / study.capacitance_f -1.0 / (study.load_resistance_ohm * study.capacitance_f)
    ]
end

function _average_boost_stored_energy(study, inductor_current, output_voltage)
    return 0.5 * study.inductance_h * inductor_current^2 +
        0.5 * study.capacitance_f * output_voltage^2
end

function prepare_average_boost_converter(study::AverageBoostConverterStudy)
    step = study.specification.timing.fixed_step_s
    state_matrix = _average_boost_state_matrix(study)
    identity_matrix = Matrix{Float64}(I, 2, 2)
    left_matrix = identity_matrix - 0.5 * step * state_matrix
    transition_matrix = left_matrix \ (identity_matrix + 0.5 * step * state_matrix)
    source_increment = left_matrix \ [step * study.input_voltage_v / study.inductance_h, 0.0]
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / step,
    ) + 1
    traces = ntuple(_ -> zeros(Float64, sample_count), 10)
    runtime = AverageBoostConverterRuntime(
        study,
        study.initial_state.inductor_current_a,
        study.initial_state.output_voltage_v,
        study.start_time_s,
        0.0,
        0,
        transition_matrix,
        source_increment,
        traces...,
    )
    _record_average_boost_sample!(runtime, 1, 0.0, 0.0)
    return runtime
end

function _record_average_boost_sample!(runtime, sample_index, kcl_residual_a, energy_residual_w)
    study = runtime.study
    current = runtime.inductor_current_a
    voltage = runtime.output_voltage_v
    runtime.time_trace_s[sample_index] = runtime.time_s
    runtime.input_voltage_trace_v[sample_index] = study.input_voltage_v
    runtime.input_current_trace_a[sample_index] = current
    runtime.output_voltage_trace_v[sample_index] = voltage
    runtime.inductor_current_trace_a[sample_index] = current
    runtime.load_current_trace_a[sample_index] = voltage / study.load_resistance_ohm
    runtime.stored_energy_trace_j[sample_index] =
        _average_boost_stored_energy(study, current, voltage)
    runtime.dissipated_energy_trace_j[sample_index] = runtime.dissipated_energy_j
    runtime.kcl_residual_trace_a[sample_index] = kcl_residual_a
    runtime.energy_residual_trace_w[sample_index] = energy_residual_w
    return runtime
end

function _advance_average_boost_converter!(runtime, sample_index)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    transfer_fraction = 1.0 - study.specification.modulation.duty
    series_resistance = study.source_resistance_ohm + study.inductor_resistance_ohm
    previous_current = runtime.inductor_current_a
    previous_voltage = runtime.output_voltage_v
    previous_energy = _average_boost_stored_energy(
        study,
        previous_current,
        previous_voltage,
    )
    transition = runtime.transition_matrix
    runtime.inductor_current_a = transition[1, 1] * previous_current +
        transition[1, 2] * previous_voltage + runtime.source_increment[1]
    runtime.output_voltage_v = transition[2, 1] * previous_current +
        transition[2, 2] * previous_voltage + runtime.source_increment[2]
    runtime.inductor_current_a >= 0.0 && runtime.output_voltage_v >= 0.0 ||
        throw(DomainError(
            (runtime.inductor_current_a, runtime.output_voltage_v),
            "average boost state left its nonnegative continuous-conduction domain",
        ))
    current = runtime.inductor_current_a
    voltage = runtime.output_voltage_v
    runtime.time_s = study.start_time_s + (sample_index - 1) * step
    input_energy = 0.5 * step * study.input_voltage_v *
        (previous_current + current)
    load_energy = 0.5 * step *
        (previous_voltage^2 + voltage^2) / study.load_resistance_ohm
    loss_energy = 0.5 * step * series_resistance *
        (previous_current^2 + current^2)
    stored_energy = _average_boost_stored_energy(study, current, voltage)
    energy_residual = (
        input_energy - load_energy - loss_energy - (stored_energy - previous_energy)
    ) / step
    capacitor_current = study.capacitance_f * (voltage - previous_voltage) / step
    average_inductor_current = 0.5 * (previous_current + current)
    average_output_voltage = 0.5 * (previous_voltage + voltage)
    kcl_residual = transfer_fraction * average_inductor_current -
        average_output_voltage / study.load_resistance_ohm - capacitor_current
    runtime.dissipated_energy_j += loss_energy
    return _record_average_boost_sample!(
        runtime,
        sample_index,
        kcl_residual,
        energy_residual,
    )
end

function _average_boost_result_state(runtime)
    study = runtime.study
    current = runtime.inductor_current_a
    voltage = runtime.output_voltage_v
    stored_energy = _average_boost_stored_energy(study, current, voltage)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(current),
        repr(voltage),
        repr(stored_energy),
        repr(runtime.dissipated_energy_j),
    ), '\n')))
    return ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, voltage],
        [current, voltage / study.load_resistance_ohm],
        BitVector(),
        BitVector(),
        BitVector(),
        [study.capacitance_f * voltage],
        [study.inductance_h * current],
        Float64[],
        [study.specification.modulation.duty],
        stored_energy,
        runtime.dissipated_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
end

function execute_average_boost_converter!(runtime::AverageBoostConverterRuntime)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    maximum_kcl = maximum(abs, runtime.kcl_residual_trace_a; init=0.0)
    energy_scale = max(
        maximum(abs, runtime.input_voltage_trace_v .* runtime.input_current_trace_a; init=0.0),
        maximum(abs, runtime.output_voltage_trace_v .* runtime.load_current_trace_a; init=0.0),
        1.0,
    )
    result = ConverterSystems.converter_system_result(
        runtime.study.specification,
        _average_boost_result_state(runtime);
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=maximum_kcl,
        relative_charge_residual=maximum_kcl,
        relative_energy_residual=
            maximum(abs, runtime.energy_residual_trace_w; init=0.0) / energy_scale,
    )
    return ConverterSystems.AverageBoostConverterTrace(
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

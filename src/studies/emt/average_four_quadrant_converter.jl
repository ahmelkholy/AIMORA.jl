using ..ConverterSystems

mutable struct AverageFullBridgeConverterRuntime{S}
    study::S
    load_current_a::Float64
    time_s::Float64
    source_dissipated_energy_j::Float64
    accepted_step_index::Int
    transition_factor::Float64
    source_gain_a_per_v::Float64
    time_trace_s::Vector{Float64}
    input_voltage_trace_v::Vector{Float64}
    input_current_trace_a::Vector{Float64}
    output_voltage_trace_v::Vector{Float64}
    load_current_trace_a::Vector{Float64}
    stored_energy_trace_j::Vector{Float64}
    source_dissipated_energy_trace_j::Vector{Float64}
    circuit_residual_trace_v::Vector{Float64}
    energy_residual_trace_w::Vector{Float64}
end

function _average_full_bridge_polarity(study, time_s)
    selection = study.specification.selection
    modulation = study.specification.modulation
    selection.family === ConverterSystems.FourQuadrantChopper &&
        return 2.0 * modulation.duty - 1.0
    selection.family === ConverterSystems.SinglePhaseTwoLevelBridge ||
        throw(ArgumentError("average full-bridge runtime received an unsupported family"))
    angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
        modulation.phase_shift_rad
    return modulation.modulation_index * sin(angle)
end

function prepare_average_full_bridge_converter(study::AverageFullBridgeConverterStudy)
    step = study.specification.timing.fixed_step_s
    total_resistance = study.source_resistance_ohm + study.load_resistance_ohm
    rate = -total_resistance / study.load_inductance_h
    left = 1.0 - 0.5 * step * rate
    transition_factor = (1.0 + 0.5 * step * rate) / left
    source_gain = 0.5 * step / study.load_inductance_h / left
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / step,
    ) + 1
    traces = ntuple(_ -> zeros(Float64, sample_count), 9)
    runtime = AverageFullBridgeConverterRuntime(
        study,
        study.initial_state.load_current_a,
        study.start_time_s,
        0.0,
        0,
        transition_factor,
        source_gain,
        traces...,
    )
    _record_average_full_bridge_sample!(runtime, 1, 0.0, 0.0)
    return runtime
end

function _record_average_full_bridge_sample!(
    runtime,
    sample_index,
    circuit_residual_v,
    energy_residual_w,
)
    study = runtime.study
    polarity = _average_full_bridge_polarity(study, runtime.time_s)
    current = runtime.load_current_a
    input_current = polarity * current
    output_voltage = polarity * study.input_voltage_v -
        study.source_resistance_ohm * current
    runtime.time_trace_s[sample_index] = runtime.time_s
    runtime.input_voltage_trace_v[sample_index] = study.input_voltage_v
    runtime.input_current_trace_a[sample_index] = input_current
    runtime.output_voltage_trace_v[sample_index] = output_voltage
    runtime.load_current_trace_a[sample_index] = current
    runtime.stored_energy_trace_j[sample_index] =
        0.5 * study.load_inductance_h * current^2
    runtime.source_dissipated_energy_trace_j[sample_index] =
        runtime.source_dissipated_energy_j
    runtime.circuit_residual_trace_v[sample_index] = circuit_residual_v
    runtime.energy_residual_trace_w[sample_index] = energy_residual_w
    return runtime
end

function _advance_average_full_bridge_converter!(runtime, sample_index)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    previous_polarity = _average_full_bridge_polarity(study, runtime.time_s)
    previous_current = runtime.load_current_a
    previous_energy = 0.5 * study.load_inductance_h * previous_current^2
    endpoint = study.start_time_s + (sample_index - 1) * step
    polarity = _average_full_bridge_polarity(study, endpoint)
    runtime.load_current_a =
        runtime.transition_factor * previous_current +
        runtime.source_gain_a_per_v * study.input_voltage_v *
        (previous_polarity + polarity)
    current = runtime.load_current_a
    runtime.time_s = endpoint
    average_current = 0.5 * (previous_current + current)
    average_polarity = 0.5 * (previous_polarity + polarity)
    current_derivative = (current - previous_current) / step
    average_output_voltage = average_polarity * study.input_voltage_v -
        study.source_resistance_ohm * average_current
    circuit_residual = average_output_voltage -
        study.load_resistance_ohm * average_current -
        study.load_inductance_h * current_derivative
    input_energy = step * average_polarity * study.input_voltage_v * average_current
    load_energy = step * study.load_resistance_ohm * average_current^2
    source_loss_energy = step * study.source_resistance_ohm * average_current^2
    stored_energy = 0.5 * study.load_inductance_h * current^2
    energy_residual = (
        input_energy - load_energy - source_loss_energy -
        (stored_energy - previous_energy)
    ) / step
    runtime.source_dissipated_energy_j += source_loss_energy
    return _record_average_full_bridge_sample!(
        runtime,
        sample_index,
        circuit_residual,
        energy_residual,
    )
end

function _average_full_bridge_result_state(runtime)
    study = runtime.study
    polarity = _average_full_bridge_polarity(study, runtime.time_s)
    current = runtime.load_current_a
    output_voltage = polarity * study.input_voltage_v -
        study.source_resistance_ohm * current
    stored_energy = 0.5 * study.load_inductance_h * current^2
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(current),
        repr(output_voltage),
        repr(stored_energy),
        repr(runtime.source_dissipated_energy_j),
    ), '\n')))
    return ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, output_voltage],
        [polarity * current, current],
        BitVector(),
        BitVector(),
        BitVector(),
        Float64[],
        [study.load_inductance_h * current],
        Float64[],
        [polarity],
        stored_energy,
        runtime.source_dissipated_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
end

function execute_average_full_bridge_converter!(runtime::AverageFullBridgeConverterRuntime)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    voltage_scale = max(
        maximum(abs, runtime.output_voltage_trace_v; init=0.0),
        runtime.study.input_voltage_v,
        1.0,
    )
    energy_scale = max(
        maximum(abs, runtime.input_voltage_trace_v .* runtime.input_current_trace_a; init=0.0),
        maximum(abs, runtime.output_voltage_trace_v .* runtime.load_current_trace_a; init=0.0),
        1.0,
    )
    relative_circuit_residual =
        maximum(abs, runtime.circuit_residual_trace_v; init=0.0) / voltage_scale
    result = ConverterSystems.converter_system_result(
        runtime.study.specification,
        _average_full_bridge_result_state(runtime);
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=0.0,
        relative_charge_residual=relative_circuit_residual,
        relative_energy_residual=
            maximum(abs, runtime.energy_residual_trace_w; init=0.0) / energy_scale,
    )
    return ConverterSystems.AverageFullBridgeConverterTrace(
        runtime.time_trace_s,
        runtime.input_voltage_trace_v,
        runtime.input_current_trace_a,
        runtime.output_voltage_trace_v,
        runtime.load_current_trace_a,
        runtime.stored_energy_trace_j,
        runtime.source_dissipated_energy_trace_j,
        runtime.circuit_residual_trace_v,
        runtime.energy_residual_trace_w,
        result,
    )
end

using ..ConverterSystems

mutable struct AverageDualActiveBridgeRuntime{S}
    study::S
    primary_dc_voltage_v::Float64
    secondary_dc_voltage_v::Float64
    primary_dc_current_a::Float64
    secondary_dc_current_a::Float64
    transferred_power_w::Float64
    source_dissipated_energy_j::Float64
    accepted_step_index::Int
    time_trace_s::Vector{Float64}
    primary_dc_voltage_trace_v::Vector{Float64}
    secondary_dc_voltage_trace_v::Vector{Float64}
    primary_dc_current_trace_a::Vector{Float64}
    secondary_dc_current_trace_a::Vector{Float64}
    transferred_power_trace_w::Vector{Float64}
    source_dissipated_energy_trace_j::Vector{Float64}
    energy_residual_trace_w::Vector{Float64}
end

function _average_dab_port_current(source_voltage_v, source_resistance_ohm, terminal_power_w)
    discriminant = source_voltage_v^2 -
        4.0 * source_resistance_ohm * terminal_power_w
    discriminant >= 0.0 || throw(ArgumentError(
        "average DAB terminal power exceeds its source-resistance limit",
    ))
    denominator = source_voltage_v + sqrt(discriminant)
    return iszero(terminal_power_w) ? 0.0 : 2.0 * terminal_power_w / denominator
end

function prepare_average_dual_active_bridge(study::AverageDualActiveBridgeStudy)
    angular_frequency = 2.0 * pi * study.specification.timing.carrier_frequency_hz
    transferred_power = ConverterSystems.dual_active_bridge_power_w(
        study.primary_dc_voltage_v,
        study.secondary_dc_voltage_v,
        study.transformer_link.turns_ratio,
        study.specification.modulation.phase_shift_rad,
        angular_frequency,
        study.transformer_link.leakage_inductance_h,
    )
    primary_current = _average_dab_port_current(
        study.primary_dc_voltage_v,
        study.primary_source_resistance_ohm,
        transferred_power,
    )
    secondary_current = _average_dab_port_current(
        study.secondary_dc_voltage_v,
        study.secondary_source_resistance_ohm,
        -transferred_power,
    )
    primary_terminal_voltage = study.primary_dc_voltage_v -
        study.primary_source_resistance_ohm * primary_current
    secondary_terminal_voltage = study.secondary_dc_voltage_v -
        study.secondary_source_resistance_ohm * secondary_current
    step = study.specification.timing.fixed_step_s
    sample_count = round(Int, (study.stop_time_s - study.start_time_s) / step) + 1
    traces = ntuple(_ -> zeros(Float64, sample_count), 8)
    runtime = AverageDualActiveBridgeRuntime(
        study,
        primary_terminal_voltage,
        secondary_terminal_voltage,
        primary_current,
        secondary_current,
        transferred_power,
        0.0,
        0,
        traces...,
    )
    _record_average_dab_sample!(runtime, 1, 0.0)
    return runtime
end

function _record_average_dab_sample!(runtime, sample_index, energy_residual_w)
    study = runtime.study
    runtime.time_trace_s[sample_index] = study.start_time_s +
        (sample_index - 1) * study.specification.timing.fixed_step_s
    runtime.primary_dc_voltage_trace_v[sample_index] = runtime.primary_dc_voltage_v
    runtime.secondary_dc_voltage_trace_v[sample_index] = runtime.secondary_dc_voltage_v
    runtime.primary_dc_current_trace_a[sample_index] = runtime.primary_dc_current_a
    runtime.secondary_dc_current_trace_a[sample_index] = runtime.secondary_dc_current_a
    runtime.transferred_power_trace_w[sample_index] = runtime.transferred_power_w
    runtime.source_dissipated_energy_trace_j[sample_index] =
        runtime.source_dissipated_energy_j
    runtime.energy_residual_trace_w[sample_index] = energy_residual_w
    return runtime
end

function _average_dab_result_state(runtime)
    study = runtime.study
    final_time = last(runtime.time_trace_s)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        study.transformer_link.deterministic_signature_sha256,
        repr(final_time),
        repr(runtime.primary_dc_voltage_v),
        repr(runtime.secondary_dc_voltage_v),
        repr(runtime.primary_dc_current_a),
        repr(runtime.secondary_dc_current_a),
        repr(runtime.transferred_power_w),
        repr(runtime.source_dissipated_energy_j),
    ), '\n')))
    return ConverterSystems.ConverterSystemState(
        final_time,
        [runtime.primary_dc_voltage_v, runtime.secondary_dc_voltage_v],
        [runtime.primary_dc_current_a, runtime.secondary_dc_current_a],
        BitVector(),
        BitVector(),
        BitVector(),
        Float64[],
        Float64[],
        Float64[],
        [
            study.specification.modulation.phase_shift_rad,
            runtime.transferred_power_w,
        ],
        0.0,
        runtime.source_dissipated_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
end

function _advance_average_dual_active_bridge!(runtime, sample_index)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    primary_source_power = study.primary_dc_voltage_v * runtime.primary_dc_current_a
    secondary_source_power = study.secondary_dc_voltage_v * runtime.secondary_dc_current_a
    source_loss = study.primary_source_resistance_ohm * runtime.primary_dc_current_a^2 +
        study.secondary_source_resistance_ohm * runtime.secondary_dc_current_a^2
    energy_residual = primary_source_power + secondary_source_power - source_loss
    runtime.source_dissipated_energy_j += step * source_loss
    return _record_average_dab_sample!(runtime, sample_index, energy_residual)
end

function execute_average_dual_active_bridge!(runtime::AverageDualActiveBridgeRuntime)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    study = runtime.study
    primary_source_power = study.primary_dc_voltage_v * runtime.primary_dc_current_a
    secondary_source_power = study.secondary_dc_voltage_v * runtime.secondary_dc_current_a
    source_loss = study.primary_source_resistance_ohm * runtime.primary_dc_current_a^2 +
        study.secondary_source_resistance_ohm * runtime.secondary_dc_current_a^2
    energy_residual = primary_source_power + secondary_source_power - source_loss
    scale = max(
        abs(primary_source_power),
        abs(secondary_source_power),
        abs(source_loss),
        1.0,
    )
    result = ConverterSystems.converter_system_result(
        study.specification,
        _average_dab_result_state(runtime);
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=0.0,
        relative_charge_residual=0.0,
        relative_energy_residual=abs(energy_residual) / scale,
    )
    return ConverterSystems.AverageDualActiveBridgeTrace(
        runtime.time_trace_s,
        runtime.primary_dc_voltage_trace_v,
        runtime.secondary_dc_voltage_trace_v,
        runtime.primary_dc_current_trace_a,
        runtime.secondary_dc_current_trace_a,
        runtime.transferred_power_trace_w,
        runtime.source_dissipated_energy_trace_j,
        runtime.energy_residual_trace_w,
        result,
    )
end

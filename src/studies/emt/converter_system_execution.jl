function _converter_accepted_step_count(value::Integer)
    value >= 0 || throw(ArgumentError(
        "converter accepted-step count must be nonnegative",
    ))
    value <= typemax(Int) || throw(ArgumentError(
        "converter accepted-step count exceeds the host integer range",
    ))
    return Int(value)
end

function _advance_average_converter_runtime!(
    runtime,
    accepted_step_count::Integer,
    advance_step!,
)
    count = _converter_accepted_step_count(accepted_step_count)
    final_step = length(runtime.time_trace_s) - 1
    0 <= runtime.accepted_step_index <= final_step || throw(ArgumentError(
        "average converter accepted-step cursor is outside its configured horizon",
    ))
    count <= final_step - runtime.accepted_step_index || throw(ArgumentError(
        "average converter advancement exceeds its configured horizon",
    ))
    for _ in 1:count
        sample_index = runtime.accepted_step_index + 2
        advance_step!(runtime, sample_index)
        runtime.accepted_step_index += 1
    end
    return runtime
end

function _advance_switching_converter_integrator!(
    integrator,
    accepted_step_count::Integer,
    advance_step!,
)
    count = _converter_accepted_step_count(accepted_step_count)
    integrator.failed && throw(ArgumentError(
        "converter integrator is terminally failed: $(integrator.last_failure)",
    ))
    final_step = length(integrator.recorder.time_s) - 1
    0 <= integrator.accepted_step_index <= final_step || throw(ArgumentError(
        "switching converter accepted-step cursor is outside its configured horizon",
    ))
    count <= final_step - integrator.accepted_step_index || throw(ArgumentError(
        "switching converter advancement exceeds its configured horizon",
    ))
    for _ in 1:count
        advance_step!(integrator) || throw(ArgumentError(
            "switching converter stopped before the requested accepted-step boundary",
        ))
    end
    return integrator
end

advance_prepared_converter_system!(
    runtime::AverageBuckConverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_buck_converter!)

advance_prepared_converter_system!(
    runtime::AverageBoostConverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_boost_converter!)

advance_prepared_converter_system!(
    runtime::AverageInvertingBuckBoostConverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(
    runtime,
    count,
    _advance_average_inverting_buck_boost_converter!,
)

advance_prepared_converter_system!(
    runtime::AverageFullBridgeConverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_full_bridge_converter!)

advance_prepared_converter_system!(
    runtime::AverageThreePhaseTwoLevelInverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_three_phase_inverter!)

advance_prepared_converter_system!(
    runtime::AverageThreeLevelSplitLinkInverterRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_npc!)

advance_prepared_converter_system!(
    runtime::AverageDualActiveBridgeRuntime,
    count::Integer,
) = _advance_average_converter_runtime!(runtime, count, _advance_average_dual_active_bridge!)

advance_prepared_converter_system!(
    integrator::SwitchingChopperConverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_switching_chopper_converter!,
)

advance_prepared_converter_system!(
    integrator::FullBridgeConverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(integrator, count, _advance_full_bridge_converter!)

advance_prepared_converter_system!(
    integrator::InterleavedChopperIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(integrator, count, _advance_interleaved_chopper!)

advance_prepared_converter_system!(
    integrator::DualActiveBridgeIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(integrator, count, _advance_dual_active_bridge!)

advance_prepared_converter_system!(
    integrator::LineCommutatedRectifierIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_line_commutated_rectifier!,
)

advance_prepared_converter_system!(
    integrator::ThreePhaseTwoLevelInverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_three_phase_two_level_inverter!,
)

advance_prepared_converter_system!(
    integrator::ThreeLevelSplitLinkInverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_three_level_split_link_inverter!,
)

advance_prepared_converter_system!(
    integrator::FlyingCapacitorInverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_flying_capacitor_inverter!,
)

advance_prepared_converter_system!(
    integrator::CascadedHBridgeInverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(
    integrator,
    count,
    _advance_cascaded_h_bridge_inverter!,
)

advance_prepared_converter_system!(
    integrator::MatrixConverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(integrator, count, _advance_matrix_converter!)

advance_prepared_converter_system!(
    integrator::CycloconverterIntegrator,
    count::Integer,
) = _advance_switching_converter_integrator!(integrator, count, _advance_cycloconverter!)

advance_prepared_converter_system!(
    runtime::AverageConverterApplicationRuntime,
    count::Integer,
) = advance_average_converter_application!(runtime, count)

function _execute_nonlinear_nodal_step_solver(
    nonlinear_nodal_step_solver,
    context,
    current_injections,
    source_voltage_constraint_result,
)
    applicable(
        nonlinear_nodal_step_solver,
        context,
        current_injections,
        source_voltage_constraint_result,
    ) || throw(ArgumentError(
        "nonlinear nodal step solver must be callable with the EMT context, current injections, and source-voltage constraints",
    ))
    return nonlinear_nodal_step_solver(
        context,
        current_injections,
        source_voltage_constraint_result,
    )
end

function step_with_over16_boundary!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState;
    collect_step_diagnostics::Bool = true,
    over16_kwargs...,
)
    context.step_index <= context.step_count ||
        throw(ArgumentError("fixed-step EMT context is already complete"))
    step_before = context.step_index
    t_before = context.t_s
    raw_over16_kwargs = (; over16_kwargs...)
    nonlinear_nodal_step_solver = get(
        raw_over16_kwargs,
        :nonlinear_nodal_step_solver,
        nothing,
    )
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (nonlinear_nodal_step_solver = nothing,),
    )
    hybrid_substep = get(raw_over16_kwargs, :hybrid_substep, false)
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (hybrid_substep = nothing,),
    )
    hybrid_report_step_s = Float64(get(
        raw_over16_kwargs,
        :hybrid_report_step_s,
        context.dt_s,
    ))
    hybrid_report_step_s > 0.0 && isfinite(hybrid_report_step_s) ||
        throw(ArgumentError("hybrid report step must be finite and positive"))
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (hybrid_report_step_s = nothing,),
    )
    record_presolve_voltage_state =
        get(raw_over16_kwargs, :record_presolve_voltage_state, false)
    solver_over16_kwargs = _without_presolve_trace_config(raw_over16_kwargs)
    lean_empty_float_values =
        haskey(solver_over16_kwargs, :lean_empty_float_values) ?
        solver_over16_kwargs.lean_empty_float_values : Float64[]
    presolve_voltage =
        record_presolve_voltage_state ?
        copy(context.system.v) : lean_empty_float_values
    sparse_node_group_config =
        get(solver_over16_kwargs, :sparse_node_group_config, nothing)
    sparse_node_group_enabled =
        sparse_node_group_config !== nothing &&
        get(sparse_node_group_config, :enabled, true)
    solver_over16_kwargs = _without_sparse_node_group_config(solver_over16_kwargs)
    dense_primary_nonlinear_compensation = get(
        solver_over16_kwargs,
        :dense_primary_nonlinear_compensation,
        false,
    )
    solver_over16_kwargs =
        _without_dense_primary_nonlinear_config(solver_over16_kwargs)
    seed_current_source_values =
        get(solver_over16_kwargs, :seed_current_source_values, false)
    current_source_seed_element_count =
        get(solver_over16_kwargs, :current_source_seed_element_count, nothing)
    current_source_seed_clear_nodes =
        get(solver_over16_kwargs, :current_source_seed_clear_nodes, Int[])
    solver_over16_kwargs = _without_current_source_seed(solver_over16_kwargs)
    source_config = get(solver_over16_kwargs, :source_config, nothing)
    distributed_transposed_line_config =
        get(solver_over16_kwargs, :distributed_transposed_line_config, nothing)
    frequency_dependent_line_config =
        get(solver_over16_kwargs, :frequency_dependent_line_config, nothing)
    nonlinear_current_config =
        get(solver_over16_kwargs, :nonlinear_current_config, nothing)
    nonlinear_current_compensation_enabled =
        nonlinear_current_config !== nothing &&
        get(nonlinear_current_config, :enabled, true)
    reset_current_source_values =
        get(solver_over16_kwargs, :reset_current_source_values, false)
    solver_over16_kwargs =
        _without_distributed_transposed_line_config(solver_over16_kwargs)
    solver_over16_kwargs =
        _without_frequency_dependent_line_config(solver_over16_kwargs)
    solver_over16_kwargs = _without_current_source_reset(solver_over16_kwargs)
    sync_switch_base_admittance =
        get(solver_over16_kwargs, :switch_base_admittance_from_step_context, false)
    source_function_control_values =
        source_config === nothing ?
        lean_empty_float_values :
        get(get(source_config, :kwargs, NamedTuple()), :xtcs_values, Float64[])
    source_function_row_result =
        context.source_function_runtime === nothing ?
        nothing :
        _synchronize_source_function_row_slots!(
            context.source_function_runtime,
            context.t_s,
            source_function_control_values,
        )
    source_voltage_constraint_result =
        _electromagnetic_source_voltage_constraints(
            over16_state,
            source_config,
            context.system.node_count,
            _source_voltage_constraint_time(source_config, context.t_s),
            context.source_function_runtime !== nothing &&
                _source_signal_stage_recording_active(
                    context.source_function_runtime,
                ),
            get(
                solver_over16_kwargs,
                :inactive_source_voltage_constraints,
                nothing,
            ),
        )
    presolve_line_voltage =
        _source_constrained_presolve_voltage(
            context.system.v,
            source_voltage_constraint_result,
            presolve_voltage,
            record_presolve_voltage_state,
        )
    electromagnetic_history_rhs =
        seed_current_source_values ?
        _electromagnetic_history_rhs_values(
            context,
            current_source_seed_element_count;
            history_voltage = presolve_line_voltage,
            advance_breqiv_history_currents = false,
            consume_breqiv_history_currents = false,
        ) :
        lean_empty_float_values
    current_injections =
        _electromagnetic_rhs_current_injections(
            over16_state,
            context.system.node_count,
            solver_over16_kwargs,
        )
    _apply_ideal_transformer_rhs_values!(
        current_injections,
        source_voltage_constraint_result,
    )
    stored_nonlinear_current_compensation_injections =
        nonlinear_current_compensation_enabled ?
        (
            dense_primary_nonlinear_compensation ?
            zeros(Float64, context.system.node_count) :
            _nonlinear_current_compensation_injection_values(
                over16_state,
                context.system.node_count,
            )
        ) :
        (
            haskey(
                solver_over16_kwargs,
                :inactive_nonlinear_current_compensation_values,
            ) ?
            solver_over16_kwargs.inactive_nonlinear_current_compensation_values :
            zeros(Float64, context.system.node_count)
        )
    use_initial_nonlinear_companion_current =
        nonlinear_current_compensation_enabled &&
        get(nonlinear_current_config, :seed_initial_nonlinear_state, false) &&
        over16_state.nonlinear_inverse.current_update_count == 0 &&
        over16_state.nonlinear_inverse.update_count == 0 &&
        over16_state.nonlinear_inverse.source_column_update_count == 0
    nonlinear_current_compensation_injections =
        use_initial_nonlinear_companion_current ?
        _pseudo_nonlinear_inductor_initial_companion_current_injections(
            nonlinear_current_config,
            context.system.node_count,
        ) :
        stored_nonlinear_current_compensation_injections
    current_source_seed_voltage_constraint_zero_count = 0
    if reset_current_source_values
        if seed_current_source_values
            _seed_current_source_values!(
                over16_state.source.f_values,
                electromagnetic_history_rhs,
            )
            current_source_seed_voltage_constraint_zero_count =
                _clear_current_source_nodes!(
                    over16_state.source.f_values,
                    current_source_seed_clear_nodes,
                )
        else
            _ensure_float_vector_length!(
                over16_state.source.f_values,
                context.system.node_count,
            )
            fill!(over16_state.source.f_values, 0.0)
        end
    end
    distributed_transposed_line_current_result =
        _distributed_transposed_line_deck_current_injection!(
            distributed_transposed_line_config,
            over16_state.source.f_values,
            collect_diagnostics = collect_step_diagnostics,
        )
    if distributed_transposed_line_current_result !== nothing
        _add_current_injection_delta!(
            current_injections,
            distributed_transposed_line_current_result.rhs_before_values,
            distributed_transposed_line_current_result.rhs_after_values,
        )
    end
    frequency_dependent_line_current_result =
        _frequency_dependent_line_deck_current_injection!(
            frequency_dependent_line_config,
            over16_state.source.f_values,
        )
    if frequency_dependent_line_current_result !== nothing
        _add_current_injection_delta!(
            current_injections,
            frequency_dependent_line_current_result.rhs_before_values,
            frequency_dependent_line_current_result.rhs_after_values,
        )
    end
    nonlinear_current_compensation_base_rhs =
        nonlinear_current_compensation_enabled ?
        (
            reset_current_source_values ?
            Float64.(over16_state.source.f_values) :
            _source_rhs_without_nonlinear_current_compensation(
                over16_state.source.f_values,
                dense_primary_nonlinear_compensation ?
                    over16_state.source.nonlinear_current_compensation_values :
                    nonlinear_current_compensation_injections,
            )
        ) :
        Float64[]
    if dense_primary_nonlinear_compensation
        nonlinear_current_compensation_base_rhs =
            _dense_primary_nonlinear_reference_rhs(
                nonlinear_current_compensation_base_rhs,
                nonlinear_current_config,
                context.system.node_count,
            )
    end
    nonlinear_pre_solve_update = nothing
    nonlinear_slope_sync_result = (matched_count = 0, mutation_count = 0)
    if nonlinear_current_compensation_enabled &&
       _has_pseudo_nonlinear_inductor_current(nonlinear_current_config) &&
       _has_live_saturated_transformer_nonlinear_slope_branch(context) &&
       get(nonlinear_current_config, :seed_initial_nonlinear_state, false) &&
       !dense_primary_nonlinear_compensation
        nonlinear_pre_solve_rhs = copy(nonlinear_current_compensation_base_rhs)
        _add_current_injection_prefix_values!(
            nonlinear_pre_solve_rhs,
            nonlinear_current_compensation_injections,
        )
        nonlinear_pre_solve_kwargs = _prepare_over16_step_kwargs!(
            context,
            over16_state,
            context.system.v,
            (
                # This boundary mirrors only the SUBTS1 nonlinear-current
                # mutation. Pending switch operations remain owned by the
                # normal accepted-timestep pass below and must run once.
                switch_operation_enabled = false,
                nonlinear_current_config = merge(
                    nonlinear_current_config,
                    (
                        nonlinear_current_compensation_base_rhs =
                            nonlinear_current_compensation_base_rhs,
                        rhs = nonlinear_pre_solve_rhs,
                    ),
                ),
            ),
            report_step_s = hybrid_report_step_s,
        )
        nonlinear_pre_solve_update = over16_accepted_timestep_update!(
            over16_state,
            context.system.v;
            collect_diagnostics = collect_step_diagnostics,
            nonlinear_pre_solve_kwargs...,
        )
        nonlinear_slope_sync_result =
            _sync_saturated_transformer_nonlinear_slope_branches!(
                context,
                nonlinear_pre_solve_kwargs.nonlinear_current_config,
                over16_state,
            )
        nonlinear_current_compensation_injections =
            _nonlinear_current_compensation_injection_values(
                over16_state,
                context.system.node_count,
            )
        solver_over16_kwargs =
            _without_nonlinear_current_config(solver_over16_kwargs)
    end
    if nonlinear_current_compensation_enabled
        _add_current_injection_values!(
            current_injections,
            nonlinear_current_compensation_injections,
        )
    end
    if nonlinear_current_compensation_enabled &&
       reset_current_source_values &&
       nonlinear_pre_solve_update === nothing
        _add_current_injection_prefix_values!(
            over16_state.source.f_values,
            nonlinear_current_compensation_injections,
        )
    end
    use_complete_source_rhs = reset_current_source_values
    complete_source_rhs =
        use_complete_source_rhs ? copy(over16_state.source.f_values) : Float64[]
    use_complete_source_rhs && _apply_ideal_transformer_rhs_values!(
        complete_source_rhs,
        source_voltage_constraint_result,
    )
    current_zero_switching =
        current_extinction_enabled(context.system.elements)
    voltage_max_before_solve = maximum(abs, context.system.v; init = 0.0)
    voltage =
        nonlinear_nodal_step_solver !== nothing ?
        _execute_nonlinear_nodal_step_solver(
            nonlinear_nodal_step_solver,
            context,
            current_injections,
            source_voltage_constraint_result,
        ) :
        sparse_node_group_enabled ?
        solve_step_with_sparse_node_groups!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
            node_group_successors = sparse_node_group_config.node_group_successors,
            grouped_switch_from_nodes = sparse_node_group_config.grouped_switch_from_nodes,
            grouped_switch_to_nodes = sparse_node_group_config.grouped_switch_to_nodes,
            complete_rhs_values =
                use_complete_source_rhs ? complete_source_rhs : nothing,
            workspace = context.sparse_node_group_workspace,
        ) :
        current_zero_switching ?
        solve_step_with_switch_state!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            switch_time_s = context.t_s,
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
        ) :
        source_voltage_constraint_result.applied ?
        solve_step_with_current_injections!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
        ) :
        solve_step_with_current_injections!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections,
        )
    if !all(isfinite, voltage)
        rhs_max = maximum(abs, context.system.rhs; init = 0.0)
        injection_max = maximum(abs, current_injections; init = 0.0)
        admittance_condition = all(isfinite, context.system.y) ?
            cond(context.system.y) : Inf
        throw(ArgumentError(
            "nodal solve produced non-finite voltage at t=$(t_before) s " *
            "(previous_max=$(voltage_max_before_solve), rhs_max=$(rhs_max), " *
            "injection_max=$(injection_max), admittance_condition=$(admittance_condition), " *
            "nonlinear_segments=$(over16_state.nonlinear_inverse.curr), " *
            "nonlinear_companions=$(over16_state.nonlinear_inverse.anonl), " *
            "nonlinear_flux_state=$(over16_state.nonlinear_inverse.vnonl))",
        ))
    end
    hysteretic_companion_admittance_stamp_count = 0
    if dense_primary_nonlinear_compensation
        hysteretic_companion_admittance_stamp_count =
            _apply_dense_primary_hysteretic_companion_admittance!(
                context,
                nonlinear_current_config,
                over16_state,
            )
        _dense_primary_nonlinear_inverse_columns!(
            over16_state,
            nonlinear_current_config,
            context,
        )
    end
    distributed_transposed_line_result =
        _distributed_transposed_line_deck_history_update!(
            distributed_transposed_line_config,
            voltage,
            distributed_transposed_line_current_result,
            collect_diagnostics = collect_step_diagnostics,
        )
    frequency_dependent_line_result =
        _frequency_dependent_line_deck_history_update!(
            frequency_dependent_line_config,
            voltage,
            frequency_dependent_line_current_result,
        )
    if sync_switch_base_admittance
        haskey(solver_over16_kwargs, :switch_admittance_config) ||
            throw(ArgumentError("step-context switch admittance sync requires switch_admittance_config"))
        _sync_step_context_switch_base_admittance!(
            over16_state,
            context,
            solver_over16_kwargs.switch_admittance_config,
        )
    end
    solver_over16_kwargs =
        _without_over16_step_context_switch_config(solver_over16_kwargs)
    if dense_primary_nonlinear_compensation &&
       haskey(solver_over16_kwargs, :nonlinear_current_config)
        nonlinear_current_config = _dense_primary_nonlinear_steady_state_seed(
            solver_over16_kwargs.nonlinear_current_config,
            over16_state,
            voltage,
        )
        solver_over16_kwargs = merge(
            solver_over16_kwargs,
            (nonlinear_current_config = nonlinear_current_config,),
        )
    end
    prepared_over16_kwargs = _prepare_over16_step_kwargs!(
            context,
            over16_state,
            voltage,
            solver_over16_kwargs,
            report_step_s = hybrid_report_step_s,
        )
    if nonlinear_current_compensation_enabled &&
       haskey(prepared_over16_kwargs, :nonlinear_current_config)
        prepared_over16_kwargs = merge(
            prepared_over16_kwargs,
            (
                nonlinear_current_config = merge(
                    prepared_over16_kwargs.nonlinear_current_config,
                    (
                        nonlinear_current_compensation_base_rhs =
                            nonlinear_current_compensation_base_rhs,
                        rhs =
                            dense_primary_nonlinear_compensation ?
                            nonlinear_current_compensation_base_rhs :
                            get(
                                prepared_over16_kwargs.nonlinear_current_config,
                                :rhs,
                                over16_state.source.f_values,
                            ),
                    ),
                ),
            ),
        )
    end
    prepared_over16_kwargs = _without_current_injection_values(prepared_over16_kwargs)
    sparse_switch_state_flow_result = nothing
    lean_sparse_switch_state_flow_result = nothing
    if _over16_step_sparse_switch_state_flow_enabled(prepared_over16_kwargs)
        sparse_step_config =
            _without_over16_sparse_switch_state_flow_config(prepared_over16_kwargs)
        sparse_max_passes = Int(get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_max_passes,
            4,
        ))
        sparse_initial_switch_operation_enabled = get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_initial_switch_operation_enabled,
            nothing,
        )
        sparse_repeat_after_post_current_queue = get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_repeat_after_post_current_queue,
            true,
        )
        if collect_step_diagnostics
            sparse_switch_state_flow_result = over16_sparse_switch_state_flow_update!(
                over16_state,
                voltage,
                sparse_step_config;
                collect_diagnostics = true,
                max_passes = sparse_max_passes,
                initial_switch_operation_enabled =
                    sparse_initial_switch_operation_enabled,
                repeat_after_post_current_queue =
                    sparse_repeat_after_post_current_queue,
            )
            over16_update = sparse_switch_state_flow_result.pass_updates[end]
        else
            lean_sparse_switch_state_flow_result =
                _over16_sparse_switch_state_flow_update_lean!(
                    over16_state,
                    voltage,
                    sparse_step_config;
                    max_passes = sparse_max_passes,
                    initial_switch_operation_enabled =
                        sparse_initial_switch_operation_enabled,
                    repeat_after_post_current_queue =
                        sparse_repeat_after_post_current_queue,
                )
            over16_update = lean_sparse_switch_state_flow_result.last_update
        end
    else
        step_config = _without_over16_sparse_switch_state_flow_config(prepared_over16_kwargs)
        over16_update = over16_accepted_timestep_update!(
            over16_state,
            voltage;
            collect_diagnostics = collect_step_diagnostics,
            step_config...,
        )
    end
    if dense_primary_nonlinear_compensation
        switched_topology_applied =
            sparse_switch_state_flow_result !== nothing ||
            lean_sparse_switch_state_flow_result !== nothing
        switched_topology_applied &&
            _apply_switched_topology_admittance!(context, over16_state)
        if any(
            ==(PSEUDO_NONLINEAR_INDUCTOR_TYPE),
            get(nonlinear_current_config, :nonlinear_types, Int[]),
        ) && _has_live_saturated_transformer_nonlinear_slope_branch(context)
            nonlinear_slope_sync_result =
                _sync_saturated_transformer_nonlinear_slope_branches!(
                    context,
                    nonlinear_current_config,
                    over16_state,
                )
        end
        voltage = _apply_dense_primary_nonlinear_solution!(
            context,
            nonlinear_current_config,
            over16_update.nonlinear_current_result,
            nonlinear_current_compensation_base_rhs,
            stored_nonlinear_current_compensation_injections,
            source_voltage_constraint_result,
        )
        if switched_topology_applied
            switch_current_config = get(
                solver_over16_kwargs,
                :switch_current_config,
                nothing,
            )
            voltage = _sync_switched_nonlinear_network_solution!(
                context,
                over16_state,
                switch_current_config,
            )
            over16_update = merge(
                over16_update,
                (
                    output_switch_network_solution =
                        copy(over16_state.switch_current.network_solution),
                    output_switch_currents =
                        copy(over16_state.switch_current.switch_currents),
                ),
            )
        end
    end
    over16_state_mutated = if sparse_switch_state_flow_result !== nothing
        sparse_switch_state_flow_result.accepted_timestep_state_mutation_count > 0
    elseif lean_sparse_switch_state_flow_result !== nothing
        lean_sparse_switch_state_flow_result.accepted_timestep_state_mutation_count > 0
    else
        over16_update.accepted_timestep_state_mutated
    end
    if nonlinear_pre_solve_update !== nothing
        over16_state_mutated =
            over16_state_mutated ||
            nonlinear_pre_solve_update.accepted_timestep_state_mutated
        over16_update = merge(
            over16_update,
            (
                nonlinear_current_result =
                    nonlinear_pre_solve_update.nonlinear_current_result,
                nonlinear_pre_solve_applied = true,
                accepted_timestep_state_mutated = over16_state_mutated,
            ),
        )
    end
    _accept_source_function_boundary_update!(
        context,
        over16_state,
        source_config,
        over16_update,
        idempotent_same_time = hybrid_substep,
    )
    trace_voltage = record_presolve_voltage_state ? presolve_voltage : voltage
    if hybrid_substep
        _update_deck_power_energy_state!(context, voltage)
    else
        record_step!(context, voltage)
    end
    collect_step_diagnostics || return nothing
    return (
        source = :emt_step_context_over16_boundary,
        outcome = :timestep_integration,
        step_index = step_before,
        t_s = t_before,
        dt_s = context.dt_s,
        voltage_pu = copy(voltage),
        trace_voltage_pu = copy(trace_voltage),
        presolve_voltage_state_recorded = record_presolve_voltage_state,
        over16_update = over16_update,
        over16_sparse_switch_state_flow_result = sparse_switch_state_flow_result,
        over16_sparse_switch_state_flow_applied =
            sparse_switch_state_flow_result !== nothing,
        step_context_recorded = !hybrid_substep,
        step_context_step_after = context.step_index,
        step_context_time_after = context.t_s,
        over16_state_mutated = over16_state_mutated,
        accepted_timestep_state_mutated = over16_state_mutated,
        saturated_transformer_nonlinear_slope_branch_count =
            nonlinear_slope_sync_result.matched_count,
        saturated_transformer_nonlinear_slope_mutation_count =
            nonlinear_slope_sync_result.mutation_count,
        nonlinear_current_pre_solve_applied =
            nonlinear_pre_solve_update !== nothing,
        source_voltage_constraint_result = source_voltage_constraint_result,
        source_voltage_constraint_count = source_voltage_constraint_result.count,
        source_voltage_constraints_applied =
            source_voltage_constraint_result.applied,
        dynamic_source_row_update_count =
            source_function_row_result === nothing ?
            0 :
            length(context.source_function_runtime.dynamic_row_indices),
        current_injection_count = length(current_injections),
        current_injection_nonzero_count = count(!=(0.0), current_injections),
        nonlinear_current_compensation_carry_values =
            copy(nonlinear_current_compensation_injections),
        hysteretic_companion_admittance_stamp_count =
            hysteretic_companion_admittance_stamp_count,
        nonlinear_current_compensation_carry_nonzero_count =
            count(!=(0.0), nonlinear_current_compensation_injections),
        nonlinear_current_compensation_base_rhs_values =
            copy(nonlinear_current_compensation_base_rhs),
        current_source_seeded = reset_current_source_values && seed_current_source_values,
        electromagnetic_history_rhs_values = copy(electromagnetic_history_rhs),
        electromagnetic_history_rhs_nonzero_count =
            count(!=(0.0), electromagnetic_history_rhs),
        current_source_seed_voltage_constraint_zero_count =
            current_source_seed_voltage_constraint_zero_count,
        distributed_transposed_line_result = distributed_transposed_line_result,
        distributed_transposed_line_update_count =
            distributed_transposed_line_result === nothing ? 0 :
            distributed_transposed_line_result.line_update_count,
        distributed_transposed_line_rhs_update_count =
            distributed_transposed_line_result === nothing ? 0 :
            distributed_transposed_line_result.rhs_update_count,
        distributed_transposed_line_state_mutated =
            distributed_transposed_line_result !== nothing &&
            distributed_transposed_line_result.state_mutated,
        frequency_dependent_line_result = frequency_dependent_line_result,
        frequency_dependent_line_update_count =
            frequency_dependent_line_result === nothing ? 0 :
            frequency_dependent_line_result.line_update_count,
        frequency_dependent_line_rhs_update_count =
            frequency_dependent_line_result === nothing ? 0 :
            frequency_dependent_line_result.rhs_update_count,
        frequency_dependent_line_state_mutated =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.state_mutated,
        frequency_dependent_line_runtime_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_dependent_line_runtime_executed,
        frequency_dependent_line_skin_effect_internal_impedance_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.skin_effect_internal_impedance_executed,
        frequency_dependent_line_earth_return_impedance_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.earth_return_impedance_executed,
        frequency_dependent_line_fitting_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_dependent_fitting_executed,
        frequency_dependent_line_frequency_loop_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_loop_executed,
        frequency_dependent_line_pipe_sheath_side_effects_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.pipe_sheath_side_effects_executed,
        legacy_fortran_in_loop = false,
        full_bpa_timestep_executed = false,
        full_deck_orchestration_executed = false,
        deferred_effects = (
            :full_init_step_calc_elec_orchestration,
            :deck_card_io,
            :solvum,
            :report_file_writers,
            :external_bpa_executable_waveform_comparison,
        ),
    )
end

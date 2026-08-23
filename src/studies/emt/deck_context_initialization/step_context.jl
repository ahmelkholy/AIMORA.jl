
function _saturated_transformer_nonlinear_slope_branch_batch(system::NodalSystem)
    branches = SaturatedTransformerNonlinearSlopeBranch[]
    for element in system.elements
        element isa SaturatedTransformerNonlinearSlopeBranch || continue
        push!(branches, element)
    end
    return branches
end

function _analytic_source_signal_batch(
    system::NodalSystem,
    element_names::AbstractVector{Symbol},
)
    signals = AnalyticSourceSignal[]
    names = Symbol[]
    for (element_index, element) in enumerate(system.elements)
        element isa Union{TheveninSource,CurrentInjection} || continue
        element.value isa AnalyticSourceSignal || continue
        push!(signals, element.value)
        push!(
            names,
            element_index <= length(element_names) ?
                element_names[element_index] :
                Symbol("analytic_source_", element_index),
        )
    end
    return signals, names
end

function _electromagnetic_history_execution_plan(system::NodalSystem)
    element_indices = Int[]
    kinds = ElectromagneticHistoryKind[]
    batch_indices = Int[]
    series_rl_branches = SeriesRLBranch[]
    series_rlc_branches = SeriesRLCBranch[]
    capacitor_branches = CapacitorBranch[]
    coupled_inductive_branches = CoupledInductiveBranch[]
    coupled_series_rl_branches = CoupledSeriesRLBranch[]
    breqiv_injections = BreqivHistoryInjection[]
    for (element_index, element) in enumerate(system.elements)
        kind = if element isa SeriesRLBranch
            push!(series_rl_branches, element)
            SERIES_RL_HISTORY
        elseif element isa SeriesRLCBranch
            push!(series_rlc_branches, element)
            SERIES_RLC_HISTORY
        elseif element isa CapacitorBranch
            push!(capacitor_branches, element)
            CAPACITOR_HISTORY
        elseif element isa CoupledInductiveBranch
            push!(coupled_inductive_branches, element)
            COUPLED_INDUCTIVE_HISTORY
        elseif element isa CoupledSeriesRLBranch
            push!(coupled_series_rl_branches, element)
            COUPLED_SERIES_RL_HISTORY
        elseif element isa BreqivHistoryInjection
            push!(breqiv_injections, element)
            BREQIV_HISTORY
        else
            continue
        end
        push!(element_indices, element_index)
        push!(kinds, kind)
        push!(batch_indices, if kind == SERIES_RL_HISTORY
            length(series_rl_branches)
        elseif kind == SERIES_RLC_HISTORY
            length(series_rlc_branches)
        elseif kind == CAPACITOR_HISTORY
            length(capacitor_branches)
        elseif kind == COUPLED_INDUCTIVE_HISTORY
            length(coupled_inductive_branches)
        elseif kind == COUPLED_SERIES_RL_HISTORY
            length(coupled_series_rl_branches)
        else
            length(breqiv_injections)
        end)
    end
    return ElectromagneticHistoryExecutionPlan(
        element_indices,
        kinds,
        batch_indices,
        series_rl_branches,
        series_rlc_branches,
        capacitor_branches,
        coupled_inductive_branches,
        coupled_series_rl_branches,
        breqiv_injections,
    )
end

function configure_series_rlc_alterations!(
    context::EMTStepContext,
    alterations::AbstractVector{<:SeriesRLCAlteration},
)
    isempty(context.series_rlc_alteration_records) ||
        throw(ArgumentError("series-RLC alterations cannot be reconfigured after execution"))
    context.next_series_rlc_alteration_index == 1 ||
        throw(ArgumentError("series-RLC alterations cannot be reconfigured after execution"))

    events = SeriesRLCAlteration[event for event in alterations]
    order = sortperm(
        eachindex(events);
        by = index -> (events[index].activation_time_s, index),
        alg = Base.Sort.MergeSort,
    )
    ordered_events = events[order]
    branch_indices = Vector{Int}(undef, length(ordered_events))
    horizon_tolerance = 8eps(max(abs(context.t_end_s), abs(context.dt_s)))
    for (event_index, event) in enumerate(ordered_events)
        event.activation_time_s <= context.t_end_s + horizon_tolerance ||
            throw(ArgumentError(
                "series-RLC alteration $(event.branch_name) at " *
                "$(event.activation_time_s) s is outside the EMT horizon",
            ))
        matching_indices = Int[]
        for (element_index, element) in enumerate(context.system.elements)
            element_index <= length(context.element_names) || continue
            context.element_names[element_index] == event.branch_name || continue
            element isa SeriesRLCBranch || continue
            push!(matching_indices, element_index)
        end
        length(matching_indices) == 1 ||
            throw(ArgumentError(
                "series-RLC alteration branch $(event.branch_name) must resolve " *
                "to exactly one SeriesRLCBranch",
            ))
        branch_indices[event_index] = only(matching_indices)
    end
    context.series_rlc_alterations = ordered_events
    context.series_rlc_alteration_branch_indices = branch_indices
    context.next_series_rlc_alteration_index = 1
    context.series_rlc_network_refactor_count = 0
    empty!(context.series_rlc_alteration_records)
    return context
end

function _apply_due_series_rlc_alterations!(context::EMTStepContext)
    next_index = context.next_series_rlc_alteration_index
    event_count = length(context.series_rlc_alterations)
    next_index > event_count && return 0
    tolerance = 8eps(max(abs(context.t_s), abs(context.dt_s)))
    context.series_rlc_alterations[next_index].activation_time_s <=
        context.t_s + tolerance || return 0

    refactor_count = context.series_rlc_network_refactor_count + 1
    applied_count = 0
    while next_index <= event_count
        event = context.series_rlc_alterations[next_index]
        event.activation_time_s <= context.t_s + tolerance || break
        branch_index = context.series_rlc_alteration_branch_indices[next_index]
        branch = context.system.elements[branch_index]
        branch isa SeriesRLCBranch ||
            throw(ArgumentError("resolved series-RLC alteration owner changed type"))

        previous_parameters = (branch.r, branch.l, branch.c)
        history_before = (
            branch.i_prev,
            branch.inductor_voltage_prev,
            branch.capacitor_voltage_prev,
            branch.v_prev,
            branch.i_last,
        )
        previous_conductance, _ = series_rlc_companion(
            branch.r,
            branch.l,
            branch.c,
            branch.i_prev,
            branch.inductor_voltage_prev,
            branch.capacitor_voltage_prev,
            context.dt_s,
        )
        branch.r = event.resistance_ohm
        branch.l = event.inductance_h
        branch.c = event.capacitance_f
        conductance, _ = series_rlc_companion(
            branch.r,
            branch.l,
            branch.c,
            branch.i_prev,
            branch.inductor_voltage_prev,
            branch.capacitor_voltage_prev,
            context.dt_s,
        )
        history_after = (
            branch.i_prev,
            branch.inductor_voltage_prev,
            branch.capacitor_voltage_prev,
            branch.v_prev,
            branch.i_last,
        )
        push!(
            context.series_rlc_alteration_records,
            SeriesRLCAlterationRecord(
                event.branch_name,
                branch_index,
                event.activation_time_s,
                context.t_s,
                context.step_index,
                previous_parameters[1],
                previous_parameters[2],
                previous_parameters[3],
                branch.r,
                branch.l,
                branch.c,
                previous_conductance,
                conductance,
                conductance - previous_conductance,
                history_before == history_after,
                refactor_count,
            ),
        )
        applied_count += 1
        next_index += 1
    end
    context.next_series_rlc_alteration_index = next_index
    context.series_rlc_network_refactor_count = refactor_count
    context.sparse_node_group_workspace.factor_valid = false
    return applied_count
end

function initialize_step_context(
    system::S;
    node_map::AbstractDict{Symbol,<:Integer},
    element_names::AbstractVector{Symbol} = Symbol[],
    deck_output_channel_names::AbstractVector{Symbol} = Symbol[],
    deck_output_node_indices::AbstractVector{<:Integer} = Int[],
    deck_branch_voltage_output_names::AbstractVector{Symbol} = Symbol[],
    deck_branch_voltage_output_branch_indices::AbstractVector{<:Integer} = Int[],
    deck_branch_current_output_names::AbstractVector{Symbol} = Symbol[],
    deck_branch_current_output_branch_indices::AbstractVector{<:Integer} = Int[],
    deck_branch_power_output_branch_indices::AbstractVector{<:Integer} = Int[],
    deck_time_switch_names::AbstractVector{Symbol} = Symbol[],
    deck_time_switch_from_node_indices::AbstractVector{<:Integer} = Int[],
    deck_time_switch_to_node_indices::AbstractVector{<:Integer} = Int[],
    deck_time_switch_close_time_s_values::AbstractVector{<:Real} = Float64[],
    deck_time_switch_open_time_s_values::AbstractVector{<:Real} = Float64[],
    deck_time_switch_initially_closed_flags::AbstractVector{Bool} = Bool[],
    deck_time_switch_on_conductance_values::AbstractVector{<:Real} = Float64[],
    deck_time_switch_off_conductance_values::AbstractVector{<:Real} = Float64[],
    deck_over5_switch_output_codes::AbstractVector{<:Integer} = Int[],
    source_function_runtime::Union{Nothing,SourceFunctionNetworkRuntime} = nothing,
    control_system_runtime::Union{Nothing,ControlSystemNetworkRuntime} = nothing,
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    source::AbstractString = "system",
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
) where {S<:NodalSystem}
    steps = fixed_step_count(dt_s, t_end_s)
    normalized_map = normalized_node_map(node_map, system.node_count)
    normalized_element_names = collect(element_names)
    element_output_channel_names =
        _trace_output_channel_names(system.elements, normalized_element_names)
    deck_output_names = Symbol.(deck_output_channel_names)
    deck_output_indices = Int.(deck_output_node_indices)
    length(deck_output_names) == length(deck_output_indices) ||
        throw(ArgumentError("deck output channel names and node indices must have the same length"))
    for node_index in deck_output_indices
        1 <= node_index <= system.node_count ||
            throw(ArgumentError("deck output node index $node_index is outside the nodal system"))
    end
    branch_voltage_output_names = Symbol.(deck_branch_voltage_output_names)
    branch_voltage_output_indices = Int.(deck_branch_voltage_output_branch_indices)
    length(branch_voltage_output_names) == length(branch_voltage_output_indices) ||
        throw(ArgumentError("deck branch voltage output names and branch indices must have the same length"))
    branch_current_output_names = Symbol.(deck_branch_current_output_names)
    branch_current_output_indices = Int.(deck_branch_current_output_branch_indices)
    length(branch_current_output_names) == length(branch_current_output_indices) ||
        throw(ArgumentError("deck branch current output names and branch indices must have the same length"))
    branch_power_output_indices = Int.(deck_branch_power_output_branch_indices)
    for branch_index in vcat(
        branch_voltage_output_indices,
        branch_current_output_indices,
        branch_power_output_indices,
    )
        1 <= branch_index <= length(system.elements) ||
            throw(ArgumentError("deck branch output index $branch_index is outside the element table"))
    end
    time_switch_count = length(deck_time_switch_names)
    time_switch_from_indices = Int.(deck_time_switch_from_node_indices)
    time_switch_to_indices = Int.(deck_time_switch_to_node_indices)
    time_switch_close_values = Float64.(deck_time_switch_close_time_s_values)
    time_switch_open_values = Float64.(deck_time_switch_open_time_s_values)
    time_switch_initial_flags = Bool[flag for flag in deck_time_switch_initially_closed_flags]
    time_switch_on_values = Float64.(deck_time_switch_on_conductance_values)
    time_switch_off_values = Float64.(deck_time_switch_off_conductance_values)
    switch_output_codes = Int.(deck_over5_switch_output_codes)
    length(time_switch_from_indices) == time_switch_count &&
        length(time_switch_to_indices) == time_switch_count &&
        length(time_switch_close_values) == time_switch_count &&
        length(time_switch_open_values) == time_switch_count &&
        length(time_switch_initial_flags) == time_switch_count &&
        length(time_switch_on_values) == time_switch_count &&
        length(time_switch_off_values) == time_switch_count ||
        throw(ArgumentError("deck time-switch metadata vectors must have the same length"))
    for node_index in vcat(time_switch_from_indices, time_switch_to_indices)
        0 <= node_index <= system.node_count ||
            throw(ArgumentError("deck time-switch node index $node_index is outside the nodal system"))
    end
    switch_voltage_output_indices =
        _deck_switch_voltage_report_indices(switch_output_codes, time_switch_count)
    switch_current_output_indices =
        _deck_switch_current_report_indices(switch_output_codes, time_switch_count)
    switch_voltage_output_names = _deck_switch_voltage_report_names(
        Symbol.(deck_time_switch_names),
        switch_voltage_output_indices,
    )
    switch_current_output_names = _deck_switch_current_report_names(
        Symbol.(deck_time_switch_names),
        switch_current_output_indices,
    )
    output_channel_names = vcat(
        element_output_channel_names,
        deck_output_names,
        branch_voltage_output_names,
        switch_voltage_output_names,
        switch_current_output_names,
        branch_current_output_names,
        control_system_runtime === nothing ? Symbol[] :
        control_system_runtime.switch_output_names,
        control_system_runtime === nothing ? Symbol[] :
        control_system_runtime.control_output_names,
    )
    recorded_steps = _recorded_trace_step_indices(steps, recorded_step_indices)
    sample_count = length(recorded_steps)
    analytic_source_signals, analytic_source_names =
        _analytic_source_signal_batch(system, normalized_element_names)
    context = EMTStepContext(
        system,
        analytic_source_signals,
        analytic_source_names,
        _saturated_transformer_nonlinear_slope_branch_batch(system),
        _electromagnetic_history_execution_plan(system),
        zeros(Float64, system.node_count),
        zeros(Float64, system.node_count + 1, system.node_count + 1),
        zeros(Float64, system.node_count + 1),
        Int[],
        Float64[],
        zeros(Int, system.node_count + 1),
        SparseNodeGroupSolveWorkspace(
            system.node_count;
            allocate_dense_storage = time_switch_count > 0 ||
                current_extinction_enabled(system.elements),
        ),
        String(source),
        dt_s,
        t_end_s,
        0.0,
        0,
        steps,
        normalized_map,
        ordered_node_names(normalized_map),
        normalized_element_names,
        zeros(Float64, sample_count),
        zeros(Float64, system.node_count, sample_count),
        output_channel_names,
        source_function_runtime,
        control_system_runtime,
        deck_output_indices,
        branch_voltage_output_indices,
        branch_current_output_indices,
        branch_power_output_indices,
        switch_voltage_output_indices,
        switch_current_output_indices,
        time_switch_count,
        Symbol.(deck_time_switch_names),
        time_switch_from_indices,
        time_switch_to_indices,
        time_switch_close_values,
        time_switch_open_values,
        time_switch_initial_flags,
        time_switch_on_values,
        time_switch_off_values,
        switch_output_codes,
        zeros(Int, time_switch_count),
        zeros(Float64, time_switch_count),
        zeros(Float64, time_switch_count),
        zeros(Float64, time_switch_count),
        zeros(Float64, time_switch_count),
        zeros(Float64, length(system.elements)),
        falses(length(system.elements)),
        zeros(Float64, length(system.elements)),
        zeros(Float64, time_switch_count),
        falses(time_switch_count),
        zeros(Float64, time_switch_count),
        zeros(Float64, length(output_channel_names), sample_count),
        zeros(Float64, length(output_channel_names), 1),
        fill(-Inf, system.node_count),
        zeros(Float64, system.node_count),
        fill(Inf, system.node_count),
        zeros(Float64, system.node_count),
        fill(-Inf, length(output_channel_names)),
        zeros(Float64, length(output_channel_names)),
        fill(Inf, length(output_channel_names)),
        zeros(Float64, length(output_channel_names)),
        recorded_steps,
        1,
        SeriesRLCAlteration[],
        Int[],
        1,
        SeriesRLCAlterationRecord[],
        0,
    )
    return configure_series_rlc_alterations!(context, series_rlc_alterations)
end

function _recorded_trace_step_indices(step_count::Integer, recorded_step_indices)
    final_step = Int(step_count)
    if recorded_step_indices === nothing
        return collect(0:final_step)
    end
    steps = Int[Int(step) for step in recorded_step_indices]
    isempty(steps) && throw(ArgumentError("recorded step indices must not be empty"))
    last_step = -1
    for step in steps
        0 <= step <= final_step ||
            throw(ArgumentError("recorded step index $step is outside 0:$final_step"))
        step > last_step ||
            throw(ArgumentError("recorded step indices must be strictly increasing"))
        last_step = step
    end
    return steps
end

function _recorded_trace_column!(context::EMTStepContext, step_index::Integer)
    next_index = context.trace_write_index
    next_index <= length(context.recorded_step_indices) || return 0
    requested_step = context.recorded_step_indices[next_index]
    current_step = Int(step_index)
    if requested_step == current_step
        context.trace_write_index += 1
        return next_index
    end
    requested_step > current_step ||
        throw(ArgumentError("recorded step index $requested_step was skipped"))
    return 0
end

function _update_context_extrema!(
    context::EMTStepContext,
    voltage::AbstractVector{<:Real},
    output_values::AbstractVector{<:Real},
)
    sample_time_s = context.t_s
    for node_index in eachindex(voltage)
        value = Float64(voltage[node_index])
        if value > context.node_maximum_values[node_index]
            context.node_maximum_values[node_index] = value
            context.node_maximum_times_s[node_index] = sample_time_s
        end
        if value < context.node_minimum_values[node_index]
            context.node_minimum_values[node_index] = value
            context.node_minimum_times_s[node_index] = sample_time_s
        end
    end
    for output_index in eachindex(output_values)
        value = Float64(output_values[output_index])
        if value > context.output_maximum_values[output_index]
            context.output_maximum_values[output_index] = value
            context.output_maximum_times_s[output_index] = sample_time_s
        end
        if value < context.output_minimum_values[output_index]
            context.output_minimum_values[output_index] = value
            context.output_minimum_times_s[output_index] = sample_time_s
        end
    end
    return context
end

function record_step!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64};
    update_power_energy::Bool = true,
)
    context.step_index <= context.step_count ||
        throw(ArgumentError("fixed-step EMT context is already complete"))
    length(voltage) == context.system.node_count ||
        throw(ArgumentError("voltage length must match context node count"))
    update_power_energy && _update_deck_power_energy_state!(context, voltage)
    output_values = @view context.output_step_values[:, 1]
    isempty(output_values) ||
        _record_context_outputs!(context.output_step_values, 1, context, voltage)
    _update_context_extrema!(context, voltage, output_values)
    column = _recorded_trace_column!(context, context.step_index)
    if column > 0
        context.time_s[column] = context.t_s
        for node in 1:context.system.node_count
            context.voltage_pu[node, column] = voltage[node]
        end
        context.output_pu[:, column] .= output_values
    end
    context.step_index += 1
    context.t_s = min(context.step_index, context.step_count) * context.dt_s
    return voltage
end

function step!(context::EMTStepContext)
    _apply_due_series_rlc_alterations!(context)
    _advance_source_function_network!(context)
    voltage = current_extinction_enabled(context.system.elements) ?
        solve_step_with_switch_state!(
            context.system,
            context.t_s,
            context.dt_s;
            switch_time_s = context.t_s,
        ) :
        solve_step!(context.system, context.t_s, context.dt_s)
    return record_step!(context, voltage)
end

function step!(
    context::EMTStepContext,
    current_injections::AbstractVector{<:Real},
)
    _apply_due_series_rlc_alterations!(context)
    _advance_source_function_network!(context)
    voltage = current_extinction_enabled(context.system.elements) ?
        solve_step_with_switch_state!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            switch_time_s = context.t_s,
        ) :
        solve_step_with_current_injections!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections,
        )
    return record_step!(context, voltage)
end

const OVER16_SPARSE_SWITCH_STATE_FLOW_CONFIG_KEYS = (
    :sparse_switch_state_flow_enabled,
    :sparse_switch_state_flow_max_passes,
    :sparse_switch_state_flow_initial_switch_operation_enabled,
    :sparse_switch_state_flow_repeat_after_post_current_queue,
)

const OVER16_STEP_CONTEXT_SWITCH_CONFIG_KEYS = (
    :switch_base_admittance_from_step_context,
)

function _without_over16_sparse_switch_state_flow_config(config::NamedTuple)
    return Base.structdiff(
        config,
        (
            sparse_switch_state_flow_enabled = nothing,
            sparse_switch_state_flow_max_passes = nothing,
            sparse_switch_state_flow_initial_switch_operation_enabled = nothing,
            sparse_switch_state_flow_repeat_after_post_current_queue = nothing,
        ),
    )
end

function _over16_step_sparse_switch_state_flow_enabled(config::NamedTuple)
    return haskey(config, :sparse_switch_state_flow_enabled) &&
           config.sparse_switch_state_flow_enabled
end

function _without_over16_step_context_switch_config(config::NamedTuple)
    return Base.structdiff(
        config,
        (switch_base_admittance_from_step_context = nothing,),
    )
end

function _without_distributed_transposed_line_config(config::NamedTuple)
    return Base.structdiff(config, (distributed_transposed_line_config = nothing,))
end

function _without_frequency_dependent_line_config(config::NamedTuple)
    return Base.structdiff(config, (frequency_dependent_line_config = nothing,))
end

function _without_current_injection_values(config::NamedTuple)
    return Base.structdiff(
        config,
        (
            current_injection_values = nothing,
            current_injection_workspace = nothing,
            lean_empty_float_values = nothing,
            inactive_nonlinear_current_compensation_values = nothing,
            inactive_source_voltage_constraints = nothing,
        ),
    )
end

function _without_current_source_reset(config::NamedTuple)
    return Base.structdiff(config, (reset_current_source_values = nothing,))
end

function _without_current_source_seed(config::NamedTuple)
    return Base.structdiff(
        config,
        (
            seed_current_source_values = nothing,
            current_source_seed_element_count = nothing,
            current_source_seed_clear_nodes = nothing,
        ),
    )
end

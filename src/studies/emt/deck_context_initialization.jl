
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
    return Base.structdiff(config, (current_injection_values = nothing,))
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

function _stamp_complex_branch_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    from_node::Integer,
    to_node::Integer,
    admittance::ComplexF64,
)
    a = Int(from_node)
    b = Int(to_node)
    if a != 0
        matrix[a, a] += admittance
    end
    if b != 0
        matrix[b, b] += admittance
    end
    if a != 0 && b != 0
        matrix[a, b] -= admittance
        matrix[b, a] -= admittance
    end
    return matrix
end

function _stamp_complex_phase_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    admittance::AbstractMatrix{ComplexF64},
)
    phase_count = length(from_nodes)
    length(to_nodes) == phase_count ||
        throw(ArgumentError("phase admittance terminal counts must match"))
    size(admittance) == (phase_count, phase_count) ||
        throw(ArgumentError("phase admittance matrix size must match terminal count"))
    for row in 1:phase_count, column in 1:phase_count
        value = admittance[row, column]
        from_row = Int(from_nodes[row])
        to_row = Int(to_nodes[row])
        from_column = Int(from_nodes[column])
        to_column = Int(to_nodes[column])
        if from_row != 0 && from_column != 0
            matrix[from_row, from_column] += value
        end
        if from_row != 0 && to_column != 0
            matrix[from_row, to_column] -= value
        end
        if to_row != 0 && from_column != 0
            matrix[to_row, from_column] -= value
        end
        if to_row != 0 && to_column != 0
            matrix[to_row, to_column] += value
        end
    end
    return matrix
end

function _deck_steady_state_frequency_hz(parsed::DeckParser.DeckParseResult)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    options.x_frequency_hz > 0.0 && return options.x_frequency_hz
    for row in DeckParser.deck_over5a_source_rows(parsed)
        row.sfreq > 0.0 && return Float64(row.sfreq) / (2.0 * pi)
    end
    return 60.0
end

function _deck_current_source_seed_clear_nodes(parsed::DeckParser.DeckParseResult)
    return unique(
        abs(Int(row.node_value))
        for row in DeckParser.deck_over5a_source_rows(parsed)
    )
end

function _steady_state_terminal_frequency_hz(
    partition::DeckParser.DeckSteadyStateFrequencyPartition,
    node_indices,
    default_frequency_hz::Float64,
)
    frequency_hz = nothing
    for node in node_indices
        node_index = abs(Int(node))
        node_index == 0 && continue
        1 <= node_index <= length(partition.node_frequencies_hz) ||
            throw(ArgumentError("steady-state terminal is outside the frequency partition"))
        node_frequency_hz = partition.node_frequencies_hz[node_index]
        node_frequency_hz == 0.0 && continue
        if frequency_hz === nothing
            frequency_hz = node_frequency_hz
        elseif !isapprox(
            frequency_hz,
            node_frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        )
            throw(ArgumentError("coupled steady-state terminals have different frequencies"))
        end
    end
    return frequency_hz === nothing ? default_frequency_hz : frequency_hz
end

function _deck_source_voltage_phasor(row)
    if row.iform == 11
        return complex(Float64(row.crest), 0.0)
    elseif row.iform == 14
        return Float64(row.crest) * cis(Float64(row.time1))
    end
    throw(ArgumentError("steady-state source phasors require constant or sinusoidal source rows"))
end

function _is_basic_harmonic_element(element)
    return element isa Union{
        ConductanceBranch,
        SeriesRLBranch,
        SeriesRLCBranch,
        CapacitorBranch,
        TheveninSource,
        CurrentInjection,
    }
end

function _claim_basic_harmonic_element_owner!(
    owners::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    owner::Symbol,
    name::Symbol,
    line_no::Int;
    required::Bool,
)
    candidates = Int[
        index for index in eachindex(parsed.elements)
        if owners[index] === :named_basic_element &&
           parsed.element_names[index] === name &&
           parsed.element_line_numbers[index] == line_no
    ]
    if isempty(candidates)
        required && throw(ArgumentError(
            "harmonic owner $owner has no matching basic element for $name on line $line_no",
        ))
        return owners
    end
    length(candidates) == 1 || throw(ArgumentError(
        "harmonic owner $owner is ambiguous for $name on line $line_no",
    ))
    owners[only(candidates)] = owner
    return owners
end

function _basic_harmonic_element_owners(parsed::DeckParser.DeckParseResult)
    owners = Symbol[
        _is_basic_harmonic_element(element) ?
        :named_basic_element : :not_a_basic_harmonic_element
        for element in parsed.elements
    ]
    for row in DeckParser.deck_over2_branch_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :fixed_card_scalar_branch,
            row.name,
            row.line_no;
            required=true,
        )
    end
    for row in DeckParser.deck_over5a_source_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :fixed_card_source,
            row.name,
            row.line_no;
            required=false,
        )
    end
    for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        row.reactance === missing && continue
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :universal_machine_terminal_reactance,
            Symbol(
                "machine_terminal_reactance_",
                row.machine_index,
                "_",
                row.branch_index,
            ),
            row.line_no;
            required=true,
        )
    end
    for row in DeckParser.deck_universal_machine_speed_capacitor_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :universal_machine_speed_branch,
            Symbol("machine_speed_capacitor_", row.machine_index),
            row.line_no;
            required=true,
        )
    end
    return owners
end

function _stamp_named_basic_harmonic_elements!(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    physical_frequency_hz::Float64,
    formulation::AbstractEMTHarmonicFormulation,
)
    owners = _basic_harmonic_element_owners(parsed)
    physical_angular_frequency = 2.0 * pi * physical_frequency_hz
    reactive_angular_frequency = _emt_reactive_angular_frequency(
        formulation,
        physical_angular_frequency,
    )
    for (element_index, element) in pairs(parsed.elements)
        owners[element_index] === :named_basic_element || continue
        if element isa ConductanceBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(element.g, 0.0),
            )
        elseif element isa SeriesRLBranch
            impedance = complex(
                element.r,
                reactive_angular_frequency * element.l,
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                inv(impedance),
            )
        elseif element isa SeriesRLCBranch
            physical_angular_frequency > 0.0 || throw(ArgumentError(
                "series R-L-C harmonic initialization requires positive frequency",
            ))
            impedance = complex(
                element.r,
                reactive_angular_frequency * element.l -
                    inv(reactive_angular_frequency * element.c),
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                inv(impedance),
            )
        elseif element isa CapacitorBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(0.0, reactive_angular_frequency * element.c),
            )
        elseif element isa TheveninSource
            element.value isa SinusoidalSourceSignal || throw(ArgumentError(
                "named Thevenin harmonic initialization requires a typed sinusoidal source",
            ))
            source_phasor = sinusoidal_source_peak_phasor(
                element.value,
                physical_frequency_hz,
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.node,
                0,
                complex(element.g, 0.0),
            )
            rhs[element.node] += element.g * source_phasor
        elseif element isa CurrentInjection
            element.value isa SinusoidalSourceSignal || throw(ArgumentError(
                "named harmonic current injection requires a typed sinusoidal source",
            ))
            rhs[element.node] += sinusoidal_source_peak_phasor(
                element.value,
                physical_frequency_hz,
            )
        end
    end
    return admittance
end

function _stamp_deck_branch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_over2_branch_rows(parsed)
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_value,
            row.to_node_value,
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        if row.branch_kind == :conductance
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                complex(Float64(row.conductance), 0.0),
            )
        elseif row.branch_kind == :series_rl
            impedance = complex(Float64(row.resistance), angular_frequency * Float64(row.inductance))
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                inv(impedance),
            )
        elseif row.branch_kind == :capacitor &&
               Float64(row.raw_capacitance) != Float64(row.capacitance)
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                complex(0.0, angular_frequency * Float64(row.capacitance)),
            )
        end
    end
    for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        row.machine_index in excluded_universal_machine_indices && continue
        row.reactance === missing && continue
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_value,
            row.to_node_value,
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        inductance = DeckParser.fixed_card_branch_timestep_inductance(
            parsed,
            Float64(row.reactance),
        )
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            inv(complex(0.0, angular_frequency * inductance)),
        )
    end
    for row in DeckParser.deck_universal_machine_speed_capacitor_rows(parsed)
        admittance = if row.resistance > 0.0 && row.capacitance == 0.0
            complex(inv(Float64(row.resistance)), 0.0)
        elseif row.resistance == 0.0 && row.capacitance > 0.0
            frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
                frequency_partition,
                row.capacitor_node_value,
                row.mass_node_value,
                default_frequency_hz,
            )
            angular_frequency = _emt_reactive_angular_frequency(
                formulation,
                2.0 * pi * frequency_hz,
            )
            complex(0.0, angular_frequency * Float64(row.capacitance))
        else
            throw(ArgumentError("unsupported universal-machine speed-capacitor impedance"))
        end
        _stamp_complex_branch_admittance!(
            matrix,
            row.capacitor_node_value,
            row.mass_node_value,
            admittance,
        )
    end
    for shunt in deck_switching_nonlinear_resistor_safety_shunts(parsed)
        _stamp_complex_branch_admittance!(
            matrix,
            shunt.from_node_index,
            shunt.to_node_index,
            complex(shunt.conductance_s, 0.0),
        )
    end
    return matrix
end

function _stamp_hysteretic_inductor_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_hysteretic_inductor_rows(parsed)
        current_a = Float64(row.steady_state_current_A)
        flux_wb = Float64(row.steady_state_flux_Wb)
        current_a == 0.0 && continue
        current_a > 0.0 && flux_wb > 0.0 || throw(ArgumentError(
            "a nonzero hysteretic-inductor steady-state current requires a positive steady-state flux",
        ))
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_index,
            row.to_node_index,
            default_frequency_hz,
        )
        reactive_angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        incremental_inductance_h = flux_wb / current_a
        impedance = complex(0.0, reactive_angular_frequency * incremental_inductance_h)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_index,
            row.to_node_index,
            inv(impedance),
        )
    end
    return matrix
end

function _stamp_piecewise_nonlinear_inductor_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)
        row.nonlinear_type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE || continue
        current_a = Float64(row.steady_state_current_a)
        flux_wb = Float64(row.steady_state_flux_wb)
        if current_a == 0.0 && flux_wb == 0.0
            continue
        end
        current_a != 0.0 && flux_wb != 0.0 && flux_wb / current_a > 0.0 ||
            throw(ArgumentError(
                "a nonzero piecewise nonlinear-inductor steady state requires finite current and flux with a positive secant inductance",
            ))
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_index,
            row.to_node_index,
            default_frequency_hz,
        )
        reactive_angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        reactive_angular_frequency > 0.0 || throw(ArgumentError(
            "piecewise nonlinear-inductor harmonic initialization requires positive frequency",
        ))
        secant_inductance_h = flux_wb / current_a
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_index,
            row.to_node_index,
            inv(complex(0.0, reactive_angular_frequency * secant_inductance_h)),
        )
    end
    return matrix
end

function _stamp_deck_induction_machine_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    definitions = DeckParser.deck_universal_machine_definition_rows(parsed)
    terminals = DeckParser.deck_universal_machine_terminal_rows(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for machine_index in sort(unique(row.machine_index for row in definitions))
        card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
        card1.machine_type in (3, 4, 5, 6, 7) || continue
        _deck_universal_machine_initialization_mode(parsed) == :manual &&
            continue
        machine_terminals = sort!(
            [row for row in terminals if row.machine_index == machine_index];
            by = row -> row.terminal_index,
        )
        active_terminals = card1.machine_type in (6, 7) ?
            machine_terminals[3:3] :
            card1.machine_type == 5 ? machine_terminals[2:3] : machine_terminals[1:3]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(
                [row.terminal_node_value for row in active_terminals],
                [row.reference_node_value for row in active_terminals],
            ),
            default_frequency_hz,
        )
        equivalent = _deck_induction_machine_steady_state_equivalent(
            parsed,
            machine_index;
            frequency_hz = frequency_hz,
        )
        for row in active_terminals
            _stamp_complex_branch_admittance!(
                matrix,
                row.terminal_node_value,
                row.reference_node_value,
                equivalent.terminal_admittance,
            )
        end
    end
    return matrix
end

function _deck_induction_machine_steady_state_equivalent(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int;
    frequency_hz::Real=_deck_steady_state_frequency_hz(parsed),
)
    card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    card2.value2 === missing &&
        throw(ArgumentError("induction-machine main inductance is missing"))
    machine_coils = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    coil_count = _deck_coupled_dq_coil_count(card1)
    length(machine_coils) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) induction-machine steady state requires $coil_count coils"))
    power_resistances = filter(>(0.0), Float64[row.resistance for row in machine_coils[1:3]])
    rotor_resistances = filter(>(0.0), Float64[row.resistance for row in machine_coils[4:coil_count]])
    rotor_inductances = filter(>(0.0), Float64[row.inductance for row in machine_coils[4:coil_count]])
    isempty(power_resistances) && throw(ArgumentError("power-coil resistance is missing"))
    isempty(rotor_resistances) && throw(ArgumentError("rotor resistance is missing"))
    isempty(rotor_inductances) && throw(ArgumentError("rotor leakage inductance is missing"))
    section = _deck_universal_machine_section(parsed)
    angular_frequency_rad_s = 2.0 * pi * Float64(frequency_hz)
    inductance_scale =
        section.parameter_basis == :power_frequency_normalized ?
        inv(angular_frequency_rad_s) : 1.0
    slip = (card4.value1 === missing ? 0.0 : Float64(card4.value1)) / 100.0
    return induction_machine_steady_state_equivalent(
        power_coil_resistance = first(power_resistances),
        rotor_resistance = first(rotor_resistances),
        rotor_leakage_inductance = first(rotor_inductances) * inductance_scale,
        main_inductance = Float64(card2.value2) * inductance_scale,
        slip = slip,
        frequency_hz = frequency_hz,
    )
end

function _stamp_deck_switch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        conductance = row.initially_closed ?
            Float64(row.on_conductance) :
            Float64(row.off_conductance)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _stamp_deck_open_switch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        if row.initially_closed
            if row.from_node_value != 0 && row.to_node_value != 0
                continue
            end
            conductance = Float64(row.on_conductance)
        else
            conductance = Float64(row.off_conductance)
        end
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _stamp_deck_switch_at_time_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    time_s::Float64,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        closed = _deck_time_switch_closed_at(
            row.initially_closed,
            Float64(row.close_time_s),
            Float64(row.open_time_s),
            time_s,
        )
        closed && row.from_node_value != 0 && row.to_node_value != 0 && continue
        conductance = closed ? Float64(row.on_conductance) :
            Float64(row.off_conductance)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _steady_state_closed_switch_representatives(
    parsed::DeckParser.DeckParseResult,
    node_count::Int,
    time_s::Union{Nothing,Float64}=nothing,
    ;
    additional_closed_node_pairs::AbstractVector{<:Tuple}=Tuple{Int,Int}[],
    include_grounded_constraints::Bool=false,
)
    representatives = collect(1:node_count)
    grounded_nodes = Int[]

    function representative(node::Int)
        root = node
        while representatives[root] != root
            root = representatives[root]
        end
        while representatives[node] != node
            parent = representatives[node]
            representatives[node] = root
            node = parent
        end
        return root
    end

    function union_nodes(left::Int, right::Int)
        left_root = representative(left)
        right_root = representative(right)
        left_root == right_root && return nothing
        retained = max(left_root, right_root)
        replaced = min(left_root, right_root)
        representatives[replaced] = retained
        return nothing
    end

    function connect_nodes(left::Int, right::Int)
        0 <= left <= node_count && 0 <= right <= node_count ||
            throw(ArgumentError("steady-state closed-switch node is outside the network"))
        left != right || return nothing
        if left == 0 || right == 0
            include_grounded_constraints && push!(grounded_nodes, max(left, right))
            return nothing
        end
        return union_nodes(left, right)
    end

    for row in DeckParser.deck_over5_switch_rows(parsed)
        closed = time_s === nothing ? row.initially_closed :
            _deck_time_switch_closed_at(
                row.initially_closed,
                Float64(row.close_time_s),
                Float64(row.open_time_s),
                time_s,
            )
        closed || continue
        from_node = Int(row.from_node_value)
        to_node = Int(row.to_node_value)
        connect_nodes(from_node, to_node)
    end
    for row in DeckParser.deck_control_system_switch_coupling_rows(parsed)
        _control_system_switch_initially_closed(row.initial_state) || continue
        row.from_node_index !== missing && row.to_node_index !== missing ||
            throw(ArgumentError(
                "steady-state controlled switch on line $(row.line_no) has unresolved nodes",
            ))
        from_node = Int(row.from_node_index)
        to_node = Int(row.to_node_index)
        connect_nodes(from_node, to_node)
    end
    for pair in additional_closed_node_pairs
        length(pair) == 2 ||
            throw(ArgumentError("additional steady-state closed-node pair must contain two nodes"))
        from_node = Int(pair[1])
        to_node = Int(pair[2])
        connect_nodes(from_node, to_node)
    end

    grounded_representatives = Set(
        representative(node) for node in grounded_nodes
    )
    for node in 1:node_count
        root = representative(node)
        representatives[node] = root in grounded_representatives ? 0 : root
    end
    return representatives
end

function _group_harmonic_network(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    node_count = length(rhs)
    size(admittance, 1) == node_count && size(admittance, 2) == node_count ||
        throw(ArgumentError("steady-state admittance and RHS dimensions must match"))
    length(representatives) == node_count ||
        throw(ArgumentError("steady-state switch representatives must cover every node"))
    all(representative -> 0 <= representative <= node_count, representatives) ||
        throw(ArgumentError("steady-state switch representative is outside the network"))
    identity_representatives = all(
        node -> Int(representatives[node]) == node,
        1:node_count,
    )
    if identity_representatives
        active_nodes = collect(1:node_count)
        reduced_index = Dict(node => node for node in active_nodes)
        active_node_groups = Vector{Int}[Int[node] for node in active_nodes]
        return (;
            active_nodes,
            reduced_index,
            active_node_groups,
            grounded_node_indices=Int[],
            switch_node_groups=active_node_groups,
            reduced_admittance=admittance,
            reduced_rhs=rhs,
        )
    end
    active_nodes = sort(filter(!iszero, unique(Int.(representatives))))
    reduced_index = Dict(node => index for (index, node) in enumerate(active_nodes))
    active_node_groups = Vector{Int}[
        findall(representative -> Int(representative) == active_node, representatives)
        for active_node in active_nodes
    ]
    grounded_node_indices = findall(iszero, representatives)
    switch_node_groups = isempty(grounded_node_indices) ?
        active_node_groups : vcat(active_node_groups, [grounded_node_indices])
    reduced_admittance = zeros(ComplexF64, length(active_nodes), length(active_nodes))
    reduced_rhs = zeros(ComplexF64, length(active_nodes))
    for row in 1:node_count
        representatives[row] == 0 && continue
        reduced_row = reduced_index[Int(representatives[row])]
        reduced_rhs[reduced_row] += rhs[row]
        for column in 1:node_count
            representatives[column] == 0 && continue
            reduced_column = reduced_index[Int(representatives[column])]
            reduced_admittance[reduced_row, reduced_column] +=
                admittance[row, column]
        end
    end
    return (;
        active_nodes,
        reduced_index,
        active_node_groups,
        grounded_node_indices,
        switch_node_groups,
        reduced_admittance,
        reduced_rhs,
    )
end

function _expand_harmonic_component(
    component::AbstractVector{<:Integer},
    active_node_groups::Vector{Vector{Int}},
)
    expanded_count = sum(
        length(active_node_groups[Int(group_index)])
        for group_index in component
    )
    expanded = Vector{Int}(undef, expanded_count)
    destination_index = 1
    for group_index in component
        group = active_node_groups[Int(group_index)]
        copyto!(expanded, destination_index, group, 1, length(group))
        destination_index += length(group)
    end
    return sort!(expanded)
end

function _solve_grouped_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
    ;
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
    passive_conductance_network::Bool=false,
)
    node_count = length(rhs)
    grouped = _group_harmonic_network(admittance, rhs, representatives)
    active_nodes = grouped.active_nodes
    reduced_index = grouped.reduced_index
    active_node_groups = grouped.active_node_groups

    reduced = isempty(active_nodes) ?
        (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            connected_components=Vector{Vector{Int}}(),
            referenced_components=BitVector(),
            unreferenced_components=Vector{Vector{Int}}(),
            classification=:unique,
        ) :
        passive_conductance_network ?
        _solve_passive_conductance_harmonic_linear_system(
            grouped.reduced_admittance,
            grouped.reduced_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        ) :
        _solve_harmonic_linear_system(
            grouped.reduced_admittance,
            grouped.reduced_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    solution = if reduced.solution === nothing
        nothing
    else
        expanded = Vector{ComplexF64}(undef, node_count)
        for node in 1:node_count
            representative = Int(representatives[node])
            expanded[node] = representative == 0 ? 0.0 + 0.0im :
                reduced.solution[reduced_index[representative]]
        end
        expanded
    end
    return (
        solution,
        numerical_rank=reduced.numerical_rank,
        condition_estimate=reduced.condition_estimate,
        maximum_residual_a=reduced.maximum_residual_a,
        relative_residual=reduced.relative_residual,
        connected_components=Vector{Int}[
            _expand_harmonic_component(component, active_node_groups)
            for component in reduced.connected_components
        ],
        referenced_components=copy(reduced.referenced_components),
        unreferenced_components=Vector{Int}[
            _expand_harmonic_component(component, active_node_groups)
            for component in reduced.unreferenced_components
        ],
        classification=reduced.classification,
        reduced_node_count=length(active_nodes),
        switch_node_groups=grouped.switch_node_groups,
    )
end

function _solve_grouped_constrained_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
    fixed_node_phasors::AbstractDict{<:Integer,<:Complex};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
    passive_conductance_network::Bool=false,
)
    isempty(fixed_node_phasors) && return _solve_grouped_harmonic_linear_system(
        admittance,
        rhs,
        representatives;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
        passive_conductance_network,
    )
    grouped = _group_harmonic_network(admittance, rhs, representatives)
    node_count = length(rhs)
    reduced_node_count = length(grouped.active_nodes)
    fixed_values = Dict{Int,ComplexF64}()
    for (node_value, phasor_value) in fixed_node_phasors
        node = Int(node_value)
        1 <= node <= node_count || throw(ArgumentError(
            "fixed harmonic node is outside the network",
        ))
        phasor = ComplexF64(phasor_value)
        all(isfinite, (real(phasor), imag(phasor))) || throw(ArgumentError(
            "fixed harmonic phasor must be finite",
        ))
        representative = Int(representatives[node])
        if representative == 0
            iszero(phasor) || throw(ArgumentError(
                "grounded closed-switch node group requires a zero fixed harmonic phasor",
            ))
            continue
        end
        reduced_node = grouped.reduced_index[representative]
        haskey(fixed_values, reduced_node) && fixed_values[reduced_node] != phasor &&
            throw(ArgumentError(
                "closed-switch node group has inconsistent fixed harmonic phasors",
            ))
        fixed_values[reduced_node] = phasor
    end
    fixed_indices = sort!(collect(keys(fixed_values)))
    unknown_indices = setdiff(collect(1:reduced_node_count), fixed_indices)
    reduced_solution = zeros(ComplexF64, reduced_node_count)
    reduced_solution[fixed_indices] .= ComplexF64[
        fixed_values[index] for index in fixed_indices
    ]
    solved = if isempty(unknown_indices)
        (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            classification=:unique,
        )
    else
        constrained_rhs = grouped.reduced_rhs[unknown_indices] -
            grouped.reduced_admittance[unknown_indices, fixed_indices] *
            reduced_solution[fixed_indices]
        passive_conductance_network ?
        _solve_passive_conductance_harmonic_linear_system(
            grouped.reduced_admittance[unknown_indices, unknown_indices],
            constrained_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        ) :
        _solve_harmonic_linear_system(
            grouped.reduced_admittance[unknown_indices, unknown_indices],
            constrained_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    end
    solved.solution === nothing ||
        (reduced_solution[unknown_indices] .= solved.solution)

    matrix_scale = maximum(abs, grouped.reduced_admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * reduced_node_count * eps(Float64) *
        Float64(rank_relative_threshold_multiplier),
        floatmin(Float64),
    )
    reduced_components = _harmonic_connected_components(
        grouped.reduced_admittance,
        structural_threshold,
    )
    natural_references = _harmonic_component_references(
        grouped.reduced_admittance,
        reduced_components,
        structural_threshold,
    )
    fixed_index_set = Set(fixed_indices)
    referenced_components = BitVector(
        natural_reference || any(node -> node in fixed_index_set, component)
        for (component, natural_reference) in
            zip(reduced_components, natural_references)
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in
            zip(reduced_components, referenced_components)
        if !referenced
    ]
    classification = if solved.classification === :infeasible
        :infeasible
    elseif solved.numerical_rank < length(unknown_indices)
        isempty(unreferenced_components) ? :nonunique : :islanded
    else
        solved.classification
    end
    expanded_solution = classification === :unique ? ComplexF64[
        Int(representatives[node]) == 0 ? 0.0 + 0.0im :
        reduced_solution[grouped.reduced_index[Int(representatives[node])]]
        for node in 1:node_count
    ] : nothing
    reaction_currents = expanded_solution === nothing ? ComplexF64[] :
        admittance * expanded_solution - rhs
    return (
        solution=expanded_solution,
        numerical_rank=solved.numerical_rank + length(fixed_indices),
        condition_estimate=solved.condition_estimate,
        maximum_residual_a=solved.maximum_residual_a,
        relative_residual=solved.relative_residual,
        connected_components=Vector{Int}[
            _expand_harmonic_component(
                component,
                grouped.active_node_groups,
            ) for component in reduced_components
        ],
        referenced_components,
        unreferenced_components=Vector{Int}[
            _expand_harmonic_component(
                component,
                grouped.active_node_groups,
            ) for component in unreferenced_components
        ],
        classification,
        reduced_node_count,
        switch_node_groups=grouped.switch_node_groups,
        constraint_reaction_current_phasors=reaction_currents,
    )
end

function _grouped_steady_state_diagnostics(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    result = _solve_grouped_harmonic_linear_system(
        admittance,
        rhs,
        representatives;
        current_absolute_a=1.0e-12,
        current_relative=1.0e-10,
        rank_relative_threshold_multiplier=10.0,
        maximum_condition_estimate=Inf,
    )
    result.classification === :unique || throw(ArgumentError(
        "grouped steady-state network classification $(result.classification): " *
        "rank $(result.numerical_rank)/$(result.reduced_node_count), " *
        "condition $(result.condition_estimate), residual $(result.maximum_residual_a) A",
    ))
    return result
end

function _solve_grouped_steady_state_admittance(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    return something(_grouped_steady_state_diagnostics(
        admittance,
        rhs,
        representatives,
    ).solution)
end

function _solve_steady_state_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
)
    result = _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a=1.0e-12,
        current_relative=1.0e-10,
        rank_relative_threshold_multiplier=10.0,
        maximum_condition_estimate=Inf,
    )
    result.classification === :unique || throw(ArgumentError(
        "steady-state network classification $(result.classification): " *
        "rank $(result.numerical_rank)/$(size(admittance, 1)), " *
        "condition $(result.condition_estimate), residual $(result.maximum_residual_a) A",
    ))
    return something(result.solution)
end

function _harmonic_connected_components(
    admittance::AbstractMatrix{ComplexF64},
    threshold::Float64,
)
    node_count = size(admittance, 1)
    visited = falses(node_count)
    components = Vector{Vector{Int}}()
    for initial_node in 1:node_count
        visited[initial_node] && continue
        component = Int[]
        stack = Int[initial_node]
        visited[initial_node] = true
        while !isempty(stack)
            node = pop!(stack)
            push!(component, node)
            for neighbor in 1:node_count
                neighbor == node && continue
                visited[neighbor] && continue
                if abs(admittance[node, neighbor]) > threshold ||
                   abs(admittance[neighbor, node]) > threshold
                    visited[neighbor] = true
                    push!(stack, neighbor)
                end
            end
        end
        sort!(component)
        push!(components, component)
    end
    return components
end

function _harmonic_component_references(
    admittance::AbstractMatrix{ComplexF64},
    components::Vector{Vector{Int}},
    threshold::Float64,
)
    references = falses(length(components))
    for (component_index, component) in enumerate(components)
        for node in component
            row_sum = sum(admittance[node, column] for column in axes(admittance, 2))
            if abs(row_sum) > threshold
                references[component_index] = true
                break
            end
        end
    end
    return references
end

function _solve_passive_conductance_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    size(admittance) == (node_count, node_count) || throw(DimensionMismatch(
        "passive harmonic admittance and right-hand side dimensions must match",
    ))
    node_count == 0 && return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    ishermitian(admittance) || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    if all(value -> iszero(imag(value)), admittance)
        return _solve_real_passive_conductance_harmonic_linear_system(
            admittance,
            rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    end
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    matrix_scale = maximum(abs, admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(admittance, node, :))
        column_magnitude = sum(abs, view(admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_admittance = Matrix{ComplexF64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        scaled_admittance[row, column] = admittance[row, column] /
            (symmetric_scales[row] * symmetric_scales[column])
    end
    scaled_eigenvalues = eigvals!(Hermitian(scaled_admittance))
    leading_eigenvalue = maximum(scaled_eigenvalues; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_eigenvalue * rank_multiplier
    minimum_eigenvalue = minimum(scaled_eigenvalues; init=0.0)
    minimum_eigenvalue > rank_threshold || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    condition_estimate = leading_eigenvalue / minimum_eigenvalue
    physical_admittance = Matrix{ComplexF64}(admittance)
    physical_factor = cholesky!(Hermitian(physical_admittance))
    candidate_solution = physical_factor \ rhs
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    components = _harmonic_connected_components(
        admittance,
        structural_threshold,
    )
    referenced_components = _harmonic_component_references(
        admittance,
        components,
        structural_threshold,
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    return (
        solution=classification === :unique ? candidate_solution : nothing,
        numerical_rank=node_count,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

function _solve_real_passive_conductance_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    real_admittance = Matrix{Float64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        real_admittance[row, column] = real(admittance[row, column])
    end
    tridiagonal_entries =
        _real_symmetric_tridiagonal_entries(real_admittance)
    matrix_scale = maximum(abs, real_admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(real_admittance, node, :))
        column_magnitude = sum(abs, view(real_admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_eigenvalues = if tridiagonal_entries === nothing
        scaled_admittance = Matrix{Float64}(undef, node_count, node_count)
        for column in 1:node_count, row in 1:node_count
            scaled_admittance[row, column] = real_admittance[row, column] /
                (symmetric_scales[row] * symmetric_scales[column])
        end
        eigvals!(Symmetric(scaled_admittance))
    else
        scaled_diagonal = Vector{Float64}(undef, node_count)
        for node in 1:node_count
            scaled_diagonal[node] = tridiagonal_entries.diagonal[node] /
                (symmetric_scales[node] * symmetric_scales[node])
        end
        scaled_off_diagonal = Vector{Float64}(undef, node_count - 1)
        for node in 1:(node_count - 1)
            scaled_off_diagonal[node] =
                tridiagonal_entries.off_diagonal[node] /
                (symmetric_scales[node] * symmetric_scales[node + 1])
        end
        eigvals!(SymTridiagonal(scaled_diagonal, scaled_off_diagonal))
    end
    leading_eigenvalue = maximum(scaled_eigenvalues; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_eigenvalue * rank_multiplier
    minimum_eigenvalue = minimum(scaled_eigenvalues; init=0.0)
    minimum_eigenvalue > rank_threshold || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    condition_estimate = leading_eigenvalue / minimum_eigenvalue
    physical_factor = cholesky!(Symmetric(copy(real_admittance)))
    real_solution = physical_factor \ real.(rhs)
    imaginary_solution = physical_factor \ imag.(rhs)
    candidate_solution = complex.(real_solution, imaginary_solution)
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    components = tridiagonal_entries === nothing ?
        _harmonic_connected_components(admittance, structural_threshold) :
        _tridiagonal_harmonic_components(
            tridiagonal_entries.off_diagonal,
            structural_threshold,
        )
    referenced_components = tridiagonal_entries === nothing ?
        _harmonic_component_references(
            admittance,
            components,
            structural_threshold,
        ) : _tridiagonal_harmonic_component_references(
            tridiagonal_entries,
            components,
            structural_threshold,
        )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    return (
        solution=classification === :unique ? candidate_solution : nothing,
        numerical_rank=node_count,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

function _real_symmetric_tridiagonal_entries(
    admittance::Matrix{Float64},
)
    node_count = size(admittance, 1)
    size(admittance, 2) == node_count || return nothing
    for column in 3:node_count, row in 1:(column - 2)
        iszero(admittance[row, column]) || return nothing
    end
    diagonal = Vector{Float64}(undef, node_count)
    off_diagonal = Vector{Float64}(undef, max(node_count - 1, 0))
    for node in 1:node_count
        diagonal[node] = admittance[node, node]
    end
    for node in 1:(node_count - 1)
        off_diagonal[node] = admittance[node, node + 1]
    end
    return (; diagonal, off_diagonal)
end

function _tridiagonal_harmonic_components(
    off_diagonal::Vector{Float64},
    threshold::Float64,
)
    node_count = length(off_diagonal) + 1
    components = Vector{Vector{Int}}()
    first_node = 1
    for node in eachindex(off_diagonal)
        abs(off_diagonal[node]) > threshold && continue
        push!(components, collect(first_node:node))
        first_node = node + 1
    end
    push!(components, collect(first_node:node_count))
    return components
end

function _tridiagonal_harmonic_component_references(
    entries::NamedTuple,
    components::Vector{Vector{Int}},
    threshold::Float64,
)
    node_count = length(entries.diagonal)
    references = falses(length(components))
    for (component_index, component) in enumerate(components)
        for node in component
            row_sum = entries.diagonal[node]
            node > 1 && (row_sum += entries.off_diagonal[node - 1])
            node < node_count && (row_sum += entries.off_diagonal[node])
            if abs(row_sum) > threshold
                references[component_index] = true
                break
            end
        end
    end
    return references
end

function _solve_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    size(admittance) == (node_count, node_count) || throw(DimensionMismatch(
        "harmonic admittance and right-hand side dimensions must match",
    ))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), admittance) ||
        throw(ArgumentError("harmonic admittance must be finite"))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), rhs) ||
        throw(ArgumentError("harmonic right-hand side must be finite"))
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    isfinite(absolute_tolerance) && absolute_tolerance > 0.0 ||
        throw(ArgumentError("harmonic absolute-current tolerance must be finite and positive"))
    isfinite(relative_tolerance) && relative_tolerance > 0.0 ||
        throw(ArgumentError("harmonic relative-current tolerance must be finite and positive"))
    isfinite(rank_multiplier) && rank_multiplier > 0.0 ||
        throw(ArgumentError("harmonic rank multiplier must be finite and positive"))
    (isfinite(condition_limit) || condition_limit == Inf) && condition_limit > 1.0 ||
        throw(ArgumentError("harmonic condition limit must exceed one"))
    if node_count == 0
        return (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            connected_components=Vector{Vector{Int}}(),
            referenced_components=BitVector(),
            unreferenced_components=Vector{Vector{Int}}(),
            classification=:unique,
        )
    end

    matrix_scale = maximum(abs, admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    components = _harmonic_connected_components(admittance, structural_threshold)
    referenced_components = _harmonic_component_references(
        admittance,
        components,
        structural_threshold,
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]

    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(admittance, node, :))
        column_magnitude = sum(abs, view(admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_admittance = Matrix{ComplexF64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        scaled_admittance[row, column] =
            admittance[row, column] /
            (symmetric_scales[row] * symmetric_scales[column])
    end
    scaled_rhs = rhs ./ symmetric_scales
    solve_components = _harmonic_connected_components(admittance, 0.0)
    decompositions = [
        svd(scaled_admittance[component, component])
        for component in solve_components
    ]
    singular_values = reduce(
        vcat,
        (decomposition.S for decomposition in decompositions);
        init=Float64[],
    )
    leading_singular_value = maximum(singular_values; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_singular_value * rank_multiplier
    numerical_rank = count(>(rank_threshold), singular_values)
    condition_estimate = numerical_rank == node_count ?
        leading_singular_value / minimum(singular_values) : Inf
    candidate_solution = zeros(ComplexF64, node_count)
    if numerical_rank == node_count
        physical_admittance = Matrix{ComplexF64}(admittance)
        physical_factor = lu(physical_admittance)
        candidate_solution .= physical_factor \ rhs
    else
        for (component, decomposition) in zip(solve_components, decompositions)
            component_rhs = scaled_rhs[component]
            all(iszero, component_rhs) && continue
            candidate_solution[component] .=
                (decomposition \ component_rhs) ./
                symmetric_scales[component]
        end
    end
    for component in solve_components
        if all(iszero, @view rhs[component])
            fill!(@view(candidate_solution[component]), 0.0 + 0.0im)
        end
    end
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif numerical_rank < node_count
        isempty(unreferenced_components) ? :nonunique : :islanded
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    solution = classification === :unique ? candidate_solution : nothing
    return (
        solution,
        numerical_rank,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

function _stamp_deck_source_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    source_conductance_by_node = Dict{Int,Float64}()
    for element in parsed.elements
        element isa TheveninSource || continue
        source_conductance_by_node[Int(element.node)] = Float64(element.g)
    end
    for row in DeckParser.deck_over5a_source_rows(parsed)
        node_value = Int(row.node_value)
        target_node = abs(node_value)
        if node_value < 0
            abs(Int(row.iform)) in (11, 14) || continue
            rhs[target_node] += _deck_source_voltage_phasor(row)
            continue
        end
        conductance = get(source_conductance_by_node, target_node, 1.0e12)
        admittance = complex(conductance, 0.0)
        _stamp_complex_branch_admittance!(matrix, target_node, 0, admittance)
        rhs[target_node] += admittance * _deck_source_voltage_phasor(row)
    end
    return matrix
end

function _stamp_coupled_lumped_sequence_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    input_frequency_hz = DeckParser.deck_fixed_time_horizon_options(parsed).x_frequency_hz
    input_frequency_hz > 0.0 || (input_frequency_hz = default_frequency_hz)
    for impedance in DeckParser.deck_coupled_lumped_sequence_impedances(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(impedance.from_node_indices, impedance.to_node_indices),
            default_frequency_hz,
        )
        reactive_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        input_angular_frequency = 2.0 * pi * input_frequency_hz
        phase_impedance = ComplexF64.(
            impedance.phase_resistance_matrix,
            (reactive_frequency / input_angular_frequency) .*
            impedance.phase_inductance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            impedance.from_node_indices,
            impedance.to_node_indices,
            inv(phase_impedance),
        )
    end
    return matrix
end

function _generator_equivalent_modal_admittance(branch, angular_frequency::Float64)
    series_inductance = im * angular_frequency * branch.inductance_h
    damping_resistance = branch.damping_resistance_ohm
    damped_inductance =
        damping_resistance > 0.0 && branch.inductance_h > 0.0 ?
        (series_inductance * damping_resistance) /
        (series_inductance + damping_resistance) :
        series_inductance
    series_capacitance =
        branch.capacitance_f > 0.0 ?
        inv(im * angular_frequency * branch.capacitance_f) :
        complex(0.0, 0.0)
    return inv(branch.resistance_ohm + damped_inductance + series_capacitance)
end

function _stamp_generator_equivalent_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_generator_equivalent_rows(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(row.from_node_indices, row.to_node_indices),
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        nph = length(row.from_node_indices)
        zero_admittance = sum(
            _generator_equivalent_modal_admittance(branch_row.branch, angular_frequency)
            for branch_row in row.zero_mode_branches
        )
        positive_admittance = sum(
            _generator_equivalent_modal_admittance(branch_row.branch, angular_frequency)
            for branch_row in row.positive_mode_branches
        )
        diagonal =
            (zero_admittance + (nph - 1) * positive_admittance) / nph
        off_diagonal = (zero_admittance - positive_admittance) / nph
        phase_admittance = fill(off_diagonal, nph, nph)
        for phase in 1:nph
            phase_admittance[phase, phase] = diagonal
        end
        _stamp_complex_phase_admittance!(
            matrix,
            row.from_node_indices,
            row.to_node_indices,
            phase_admittance,
        )
    end
    return matrix
end

function _stamp_coupled_lumped_phase_pi_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for section in DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(section.from_node_indices, section.to_node_indices),
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        phase_impedance = ComplexF64.(
            section.phase_resistance_matrix,
            angular_frequency .* section.phase_inductance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            section.from_node_indices,
            section.to_node_indices,
            inv(phase_impedance),
        )
        phase_shunt = (0.5im * angular_frequency) .* section.phase_capacitance_matrix
        ground_nodes = zeros(Int, section.phase_count)
        _stamp_complex_phase_admittance!(
            matrix,
            section.from_node_indices,
            ground_nodes,
            phase_shunt,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            section.to_node_indices,
            ground_nodes,
            phase_shunt,
        )
    end
    return matrix
end

function _stamp_cascaded_phase_pi_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for equivalent in DeckParser.deck_cascaded_phase_pi_equivalents(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(equivalent.from_node_indices, equivalent.to_node_indices),
            default_frequency_hz,
        )
        equivalent_frequency_hz = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        ) / (2.0 * pi)
        frequency_equivalent =
            isapprox(
                equivalent.frequency_hz,
                equivalent_frequency_hz;
                atol = 1.0e-12,
                rtol = 1.0e-12,
            ) ?
            equivalent :
            cascaded_phase_pi_equivalent(
                equivalent.name,
                equivalent.blocks,
                equivalent_frequency_hz,
            )
        node_indices = vcat(
            frequency_equivalent.from_node_indices,
            frequency_equivalent.to_node_indices,
        )
        for local_column in eachindex(node_indices)
            column = node_indices[local_column]
            column == 0 && continue
            for local_row in eachindex(node_indices)
                row = node_indices[local_row]
                row == 0 && continue
                matrix[row, column] +=
                    frequency_equivalent.terminal_admittance_s[
                        local_row,
                        local_column,
                    ]
            end
        end
    end
    return matrix
end

function _stamp_distributed_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    constants_rows = DeckParser.deck_distributed_transposed_line_constants(parsed)
    modal_states = DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)
    ground_nodes = fill(0, 3)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for (index, constants) in enumerate(constants_rows)
        modal_state = modal_states[index]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(modal_state.from_node_indices, modal_state.to_node_indices),
            default_frequency_hz,
        )
        equivalent = distributed_transposed_line_steady_state_pi_equivalent(
            constants;
            steady_state_frequency_hz = frequency_hz,
            storage_start_index = 1,
            name = Symbol("mixed_frequency_distributed_line_pi_", index),
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.from_node_indices,
            modal_state.to_node_indices,
            inv(equivalent.phase_series_impedance_matrix),
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.from_node_indices,
            ground_nodes,
            equivalent.phase_shunt_admittance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.to_node_indices,
            ground_nodes,
            equivalent.phase_shunt_admittance_matrix,
        )
    end
    for element in parsed.elements
        element isa ComplexModalBergeronLine || continue
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(element.from_nodes, element.to_nodes),
            default_frequency_hz,
        )
        terminal_admittance =
            complex_modal_bergeron_steady_state_terminal_admittance(
                element,
                frequency_hz,
            )
        nodes = vcat(element.from_nodes, element.to_nodes)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] +=
                terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_semlyen_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    lines = vcat(
        DeckParser.deck_semlyen_line_elements(parsed, options.dt_s),
        DeckParser.deck_rational_frequency_line_elements(parsed, options.dt_s),
    )
    for line in lines
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(line.from_nodes, line.to_nodes),
            default_frequency_hz,
        )
        terminal_admittance =
            semlyen_line_steady_state_terminal_admittance(line, frequency_hz)
        nodes = vcat(line.from_nodes, line.to_nodes)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] += terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_sampled_frequency_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for line in DeckParser.deck_sampled_frequency_line_elements(parsed, options.dt_s)
        nodes = line isa SampledFrequencyDependentLineGroup ?
            vcat(line.from_nodes, line.to_nodes) :
            [line.a, line.b]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            nodes,
            default_frequency_hz,
        )
        terminal_admittance =
            sampled_line_steady_state_terminal_admittance(line, frequency_hz)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] += terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_saturated_transformer_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    transformer_admittance,
)
    for index in eachindex(transformer_admittance.scalar_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.scalar_branch_from_node_indices[index],
            transformer_admittance.scalar_branch_to_node_indices[index],
            transformer_admittance.scalar_branch_admittances[index],
        )
    end
    for index in eachindex(transformer_admittance.ideal_low_branch_from_node_indices)
        _stamp_complex_phase_admittance!(
            matrix,
            [
                transformer_admittance.ideal_low_branch_from_node_indices[index],
                transformer_admittance.ideal_internal_branch_from_node_indices[index],
            ],
            [
                transformer_admittance.ideal_low_branch_to_node_indices[index],
                transformer_admittance.ideal_internal_branch_to_node_indices[index],
            ],
            ComplexF64[
                transformer_admittance.ideal_low_branch_admittances[index] transformer_admittance.ideal_mutual_branch_admittances[index]
                transformer_admittance.ideal_mutual_branch_admittances[index] transformer_admittance.ideal_internal_branch_admittances[index]
            ],
        )
    end
    for index in eachindex(transformer_admittance.linearized_nonlinear_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.linearized_nonlinear_branch_from_node_indices[index],
            transformer_admittance.linearized_nonlinear_branch_to_node_indices[index],
            transformer_admittance.linearized_nonlinear_branch_admittances[index],
        )
    end
    for index in eachindex(transformer_admittance.magnetizing_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.magnetizing_branch_from_node_indices[index],
            transformer_admittance.magnetizing_branch_to_node_indices[index],
            transformer_admittance.magnetizing_branch_admittances[index],
        )
    end
    return matrix
end

function _deck_steady_state_nodal_equations(
    parsed::DeckParser.DeckParseResult,
    node_count::Int;
    include_induction_machines::Bool=true,
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
    transformer_admittance=nothing,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
    default_frequency_hz::Float64=_deck_steady_state_frequency_hz(parsed),
)
    admittance = zeros(ComplexF64, node_count, node_count)
    rhs = zeros(ComplexF64, node_count)
    _stamp_deck_branch_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        excluded_universal_machine_indices =
            excluded_universal_machine_indices,
        formulation=formulation,
    )
    _stamp_hysteretic_inductor_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition;
        formulation=formulation,
    )
    _stamp_piecewise_nonlinear_inductor_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition;
        formulation=formulation,
    )
    include_induction_machines &&
        _stamp_deck_induction_machine_steady_state_admittance!(
            admittance,
            parsed,
            frequency_partition,
        )
    _stamp_deck_open_switch_steady_state_admittance!(admittance, parsed)
    _stamp_deck_source_steady_state_admittance!(admittance, rhs, parsed)
    _stamp_named_basic_harmonic_elements!(
        admittance,
        rhs,
        parsed,
        default_frequency_hz,
        formulation,
    )
    _stamp_coupled_lumped_sequence_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_generator_equivalent_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_coupled_lumped_phase_pi_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_cascaded_phase_pi_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_distributed_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    _stamp_sampled_frequency_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    _stamp_semlyen_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    transformer_admittance === nothing ||
        _stamp_saturated_transformer_steady_state_admittance!(admittance, transformer_admittance)
    zeroed_universal_machine_branch_pairs = Tuple{Int,Int}[
        (row.from_node_value, row.to_node_value)
        for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        if row.machine_index in excluded_universal_machine_indices &&
           row.from_node_value != row.to_node_value
    ]
    switch_representatives = _steady_state_closed_switch_representatives(
        parsed,
        node_count;
        additional_closed_node_pairs =
            zeroed_universal_machine_branch_pairs,
        include_grounded_constraints=true,
    )
    return admittance, rhs, switch_representatives
end

function _deck_external_steady_state_thevenin(
    parsed::DeckParser.DeckParseResult,
    node_index::Integer,
)
    node_count = maximum(values(parsed.node_map); init = 0)
    node = Int(node_index)
    1 <= node <= node_count ||
        throw(ArgumentError("steady-state Thevenin node is outside the deck network"))
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        include_induction_machines = false,
    )
    open_circuit_phasors = _solve_grouped_steady_state_admittance(
        admittance,
        rhs,
        switch_representatives,
    )
    unit_injection = zeros(ComplexF64, node_count)
    unit_injection[node] = 1.0
    driving_point_response = _solve_grouped_steady_state_admittance(
        admittance,
        unit_injection,
        switch_representatives,
    )
    return (
        source = :deck_external_steady_state_thevenin,
        node_index = node,
        voltage_phasor = open_circuit_phasors[node],
        impedance = driving_point_response[node],
        node_voltage_phasors = open_circuit_phasors,
        admittance = admittance,
        rhs = rhs,
        switch_representatives = switch_representatives,
    )
end

function _deck_external_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult,
    current_injections::AbstractDict{<:Integer,<:Complex},
    ;
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
)
    node_map = Dict{Symbol,Int}(parsed.node_map)
    node_count = maximum(values(node_map); init = 0)
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        include_induction_machines = false,
        excluded_universal_machine_indices,
    )
    for (node_value, injection) in current_injections
        node = Int(node_value)
        1 <= node <= node_count ||
            throw(ArgumentError("steady-state current injection node is outside the deck network"))
        rhs[node] += ComplexF64(injection)
    end
    topology_diagnostics = _grouped_steady_state_diagnostics(
        admittance,
        rhs,
        switch_representatives,
    )
    node_voltage_phasors = something(topology_diagnostics.solution)
    return (
        source = :deck_external_steady_state_voltage_phasors,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_map = node_map,
        node_names = ordered_node_names(node_map),
        node_voltage_phasors = node_voltage_phasors,
        node_voltage_values = real.(node_voltage_phasors),
        topology_diagnostics,
        steady_state_admittance = admittance,
        steady_state_source_injection_phasors = rhs,
        steady_state_switch_representatives = switch_representatives,
    )
end

function deck_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult;
    saturated_transformer_intake = nothing,
    winding_number::Int = 1,
)
    node_map = Dict{Symbol,Int}(parsed.node_map)
    frequency_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    branch_assembly = nothing
    transformer_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    if saturated_transformer_intake !== nothing
        arrays = saturated_transformer_nonlinear_arrays(saturated_transformer_intake)
        physical_node_map = _saturated_transformer_physical_node_map(parsed, arrays)
        transformer_frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            _saturated_transformer_frequency_node_indices(
                physical_node_map,
                arrays,
                length(frequency_partition.node_frequencies_hz),
            ),
            transformer_frequency_hz,
        )
        branch_assembly = saturated_transformer_winding_branch_assembly(
            arrays,
            physical_node_map;
            nonlinear_winding_number = winding_number,
            reference_node_index = 0,
        )
        node_map = _saturated_transformer_augmented_node_map(
            physical_node_map,
            branch_assembly,
        )
    end
    node_count = maximum(values(node_map); init = 0)
    transformer_admittance = saturated_transformer_intake === nothing ? nothing :
        saturated_transformer_steady_state_branch_admittance(
            saturated_transformer_intake,
            physical_node_map;
            nonlinear_winding_number = winding_number,
            reference_node_index = 0,
            reactance_units = 2.0 * pi * transformer_frequency_hz,
        )
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        transformer_admittance = transformer_admittance,
        frequency_partition = frequency_partition,
    )
    topology_diagnostics = _grouped_steady_state_diagnostics(
        admittance,
        rhs,
        switch_representatives,
    )
    node_voltage_phasors = something(topology_diagnostics.solution)
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    output_voltage_values = Float64[
        index == 0 ? 0.0 : real(node_voltage_phasors[index])
        for index in output_node_indices
    ]
    return (
        source = :deck_steady_state_voltage_phasors,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_steady_state_frequencies_hz = copy(frequency_partition.node_frequencies_hz),
        node_frequency_source_row_indices =
            copy(frequency_partition.node_source_row_indices),
        source_frequency_successor_indices =
            copy(frequency_partition.source_successor_indices),
        steady_state_frequency_subnetwork_count =
            length(frequency_partition.subnetwork_node_indices),
        timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_map = node_map,
        node_names = ordered_node_names(node_map),
        node_voltage_phasors = node_voltage_phasors,
        node_voltage_values = real.(node_voltage_phasors),
        topology_diagnostics,
        steady_state_admittance = admittance,
        steady_state_source_injection_phasors = rhs,
        steady_state_switch_representatives = switch_representatives,
        output_node_indices = output_node_indices,
        output_voltage_values = output_voltage_values,
        source_row_count = length(DeckParser.deck_over5a_source_rows(parsed)),
        coupled_lumped_sequence_count =
            length(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)),
        generator_equivalent_count =
            length(DeckParser.deck_generator_equivalent_rows(parsed)),
        coupled_lumped_phase_pi_section_count =
            length(DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)),
        cascaded_phase_pi_equivalent_count =
            length(DeckParser.deck_cascaded_phase_pi_equivalents(parsed)),
        distributed_transposed_line_count =
            length(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)),
        saturated_transformer_branch_count =
            transformer_admittance === nothing ? 0 : transformer_admittance.branch_count,
        saturated_transformer_ideal_branch_count =
            transformer_admittance === nothing ? 0 : transformer_admittance.ideal_branch_count,
        saturated_transformer_linearized_nonlinear_branch_count =
            transformer_admittance === nothing ? 0 :
            transformer_admittance.linearized_nonlinear_branch_count,
        saturated_transformer_ideal_primary_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_primary_storage_coefficients),
        saturated_transformer_ideal_mutual_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_mutual_storage_coefficients),
        saturated_transformer_ideal_secondary_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_secondary_storage_coefficients),
    )
end

function _deck_node_voltage_initial_sample(parsed::DeckParser.DeckParseResult)
    rows = [
        row for row in DeckParser.deck_node_initial_condition_rows(parsed)
        if row.condition_kind == :node_voltage_initial_condition && row.node_index > 0
    ]
    isempty(rows) && return nothing
    values = zeros(Float64, length(parsed.node_map))
    for row in rows
        reference_value =
            row.reference_node_index isa Missing ? 0.0 : values[Int(row.reference_node_index)]
        values[Int(row.node_index)] = reference_value + Float64(row.real_value)
    end
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    output_voltage_values = Float64[
        index == 0 ? 0.0 : values[index]
        for index in output_node_indices
    ]
    return (
        source = :deck_node_voltage_initial_conditions,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        output_node_indices = output_node_indices,
        output_voltage_values = output_voltage_values,
        node_voltage_condition_count = length(rows),
    )
end

function _deck_zero_initial_voltage_sample(parsed::DeckParser.DeckParseResult)
    values = zeros(Float64, length(parsed.node_map))
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    return (
        source = :zero_initial_network_state,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        output_node_indices = output_node_indices,
        output_voltage_values = zeros(Float64, length(output_node_indices)),
        node_voltage_condition_count = 0,
    )
end

function _synchronous_machine_terminal_phase_angle_deg(row, base_angle_deg::Float64)
    row.angle_deg isa Missing || return Float64(row.angle_deg)
    row.phase_index == 1 && return base_angle_deg
    row.phase_index == 2 && return base_angle_deg + 240.0
    row.phase_index == 3 && return base_angle_deg + 120.0
    return base_angle_deg
end

function _deck_synchronous_machine_terminal_voltage_initial_sample(
    parsed::DeckParser.DeckParseResult,
)
    rows = DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    isempty(rows) && return nothing
    values = zeros(Float64, length(parsed.node_map))
    phasors = zeros(ComplexF64, length(parsed.node_map))
    source_count = 0
    for machine_index in unique(row.machine_index for row in rows)
        machine_rows = [row for row in rows if row.machine_index == machine_index]
        reference_row = findfirst(row -> !(row.peak_terminal_voltage isa Missing), machine_rows)
        reference_row === nothing && continue
        reference = machine_rows[reference_row]
        peak_voltage = Float64(reference.peak_terminal_voltage)
        base_angle_deg =
            reference.angle_deg isa Missing ? 0.0 : Float64(reference.angle_deg)
        for row in machine_rows
            node = row.terminal_node_value
            1 <= node <= length(values) || continue
            local_peak =
                row.peak_terminal_voltage isa Missing ?
                peak_voltage :
                Float64(row.peak_terminal_voltage)
            angle_deg = _synchronous_machine_terminal_phase_angle_deg(
                row,
                base_angle_deg,
            )
            angle_rad = deg2rad(angle_deg)
            phasor = local_peak * ComplexF64(cos(angle_rad), sin(angle_rad))
            phasors[node] = phasor
            values[node] = real(phasor)
            source_count += 1
        end
    end
    source_count > 0 || return nothing
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    return (
        source = :synchronous_machine_terminal_voltage_initial_sample,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        node_voltage_phasors = phasors,
        output_node_indices = output_node_indices,
        output_voltage_values = Float64[
            index == 0 ? 0.0 : values[index]
            for index in output_node_indices
        ],
        synchronous_machine_terminal_voltage_count = source_count,
    )
end

function _stamp_synchronous_machine_transformer_steady_state!(
    admittance::Matrix{ComplexF64},
    context,
    angular_frequency::Float64,
)
    Base.@nospecialize context
    for (name, element) in zip(context.element_names, context.system.elements)
        name_text = String(name)
        if startswith(name_text, "saturated_transformer_")
            if element isa SeriesRLBranch
                branch_admittance = inv(complex(
                    element.r,
                    angular_frequency * element.l,
                ))
                _stamp_complex_branch_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    branch_admittance,
                )
            elseif element isa CoupledInductiveBranch
                _stamp_complex_phase_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    _coupled_inductive_steady_state_admittance(element),
                )
            elseif element isa ConductanceBranch
                _stamp_complex_branch_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    complex(element.g, 0.0),
                )
            end
        elseif startswith(name_text, "transformer_branch_shunt_capacitance_") &&
               element isa CapacitorBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(0.0, angular_frequency * element.c),
            )
        end
    end
    return admittance
end

function _deck_synchronous_machine_network_initial_sample(
    parsed::DeckParser.DeckParseResult,
    context,
    ;
    strict_topology_classification::Bool=false,
)
    Base.@nospecialize context
    terminal_sample = _deck_synchronous_machine_terminal_voltage_initial_sample(parsed)
    terminal_sample === nothing && return nothing
    node_count = context.system.node_count
    admittance = zeros(ComplexF64, node_count, node_count)
    rhs = zeros(ComplexF64, node_count)
    Base.inferencebarrier(_stamp_deck_branch_steady_state_admittance!)(
        admittance,
        parsed,
    )
    time_zero_ground_fault = any(
        row -> _deck_time_switch_closed_at(
                   row.initially_closed,
                   Float64(row.close_time_s),
                   Float64(row.open_time_s),
                   0.0,
               ) && (row.from_node_value == 0 || row.to_node_value == 0),
        DeckParser.deck_over5_switch_rows(parsed),
    )
    external_excitation_port_initialization = any(
        row -> row.coupling_kind in
               (:exciter_voltage_output, :exciter_current_output),
        DeckParser.deck_synchronous_machine_control_interface_rows(parsed),
    )
    use_time_zero_topology =
        time_zero_ground_fault && external_excitation_port_initialization
    if use_time_zero_topology
        Base.inferencebarrier(_stamp_deck_switch_at_time_steady_state_admittance!)(
            admittance,
            parsed,
            0.0,
        )
    else
        Base.inferencebarrier(_stamp_deck_open_switch_steady_state_admittance!)(
            admittance,
            parsed,
        )
    end
    Base.inferencebarrier(_stamp_deck_source_steady_state_admittance!)(
        admittance,
        rhs,
        parsed,
    )
    Base.inferencebarrier(_stamp_coupled_lumped_sequence_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_generator_equivalent_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_coupled_lumped_phase_pi_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_cascaded_phase_pi_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_distributed_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_sampled_frequency_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_semlyen_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    frequency_hz = Float64(terminal_sample.steady_state_frequency_hz)
    angular_frequency = 2.0 * pi * frequency_hz
    Base.inferencebarrier(_stamp_synchronous_machine_transformer_steady_state!)(
        admittance,
        context,
        angular_frequency,
    )

    terminal_rows = DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    fixed_nodes = use_time_zero_topology ? Int[] :
        sort(unique(Int(row.terminal_node_value) for row in terminal_rows))
    fixed_phasors = ComplexF64[terminal_sample.node_voltage_phasors[node] for node in fixed_nodes]
    topology_diagnostics = nothing
    switch_representatives = nothing
    phasors = zeros(ComplexF64, node_count)
    if strict_topology_classification
        switch_representatives = _steady_state_closed_switch_representatives(
            parsed,
            node_count,
            0.0,
            include_grounded_constraints=true,
        )
        topology_diagnostics = if use_time_zero_topology
            _grouped_steady_state_diagnostics(
                admittance,
                rhs,
                switch_representatives,
            )
        else
            fixed_node_phasors = Dict{Int,ComplexF64}(
                node => phasor for (node, phasor) in zip(fixed_nodes, fixed_phasors)
            )
            result = _solve_grouped_constrained_harmonic_linear_system(
                admittance,
                rhs,
                switch_representatives,
                fixed_node_phasors;
                current_absolute_a=1.0e-12,
                current_relative=1.0e-10,
                rank_relative_threshold_multiplier=10.0,
                maximum_condition_estimate=Inf,
            )
            result.classification === :unique || throw(ArgumentError(
                "constrained synchronous-machine steady-state network classification " *
                "$(result.classification): rank $(result.numerical_rank)/" *
                "$(result.reduced_node_count), condition $(result.condition_estimate), " *
                "residual $(result.maximum_residual_a) A",
            ))
            result
        end
        phasors .= something(topology_diagnostics.solution)
    else
        unknown_nodes = setdiff(collect(1:node_count), fixed_nodes)
        phasors[fixed_nodes] .= fixed_phasors
        if use_time_zero_topology
            switch_representatives = _steady_state_closed_switch_representatives(
                parsed,
                node_count,
                0.0,
            )
            phasors .= _solve_grouped_steady_state_admittance(
                admittance,
                rhs,
                switch_representatives,
            )
        elseif !isempty(unknown_nodes)
            reduced_rhs = rhs[unknown_nodes] -
                          admittance[unknown_nodes, fixed_nodes] * fixed_phasors
            phasors[unknown_nodes] .= _solve_steady_state_linear_system(
                admittance[unknown_nodes, unknown_nodes],
                reduced_rhs,
            )
        end
    end
    node_current_phasors = admittance * phasors - rhs
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    sample = (
        source = use_time_zero_topology ?
            :synchronous_machine_time_zero_topology_steady_state :
            :synchronous_machine_constrained_steady_state,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = frequency_hz,
        node_voltage_values = real.(phasors),
        node_voltage_phasors = phasors,
        node_current_phasors = node_current_phasors,
        output_node_indices = output_node_indices,
        output_voltage_values = Float64[
            index == 0 ? 0.0 : real(phasors[index])
            for index in output_node_indices
        ],
        synchronous_machine_terminal_voltage_count = length(fixed_nodes),
        time_zero_ground_fault,
        external_excitation_port_initialization,
    )
    strict_topology_classification || return sample
    return merge(
        sample,
        (;
            topology_diagnostics,
            steady_state_admittance=admittance,
            steady_state_source_injection_phasors=rhs,
            steady_state_switch_representatives=switch_representatives,
        ),
    )
end

function _initial_voltage_sample_for_context(sample, node_count::Int)
    sample === nothing && return nothing
    values = Float64.(sample.node_voltage_values)
    length(values) >= node_count && return sample
    padded_values = zeros(Float64, node_count)
    copyto!(padded_values, 1, values, 1, length(values))
    if hasproperty(sample, :node_voltage_phasors)
        phasors = ComplexF64.(sample.node_voltage_phasors)
        padded_phasors = zeros(ComplexF64, node_count)
        copyto!(padded_phasors, 1, phasors, 1, length(phasors))
        return merge(
            sample,
            (
                node_voltage_values = padded_values,
                node_voltage_phasors = padded_phasors,
            ),
        )
    end
    return merge(sample, (node_voltage_values = padded_values,))
end

function _current_injection_samples_for_context(
    current_injection_samples,
    node_count::Int,
    sample_count::Int,
)
    current_injection_samples === nothing && return nothing
    values = Matrix{Float64}(current_injection_samples)
    size(values, 1) >= node_count ||
        throw(ArgumentError("current injection samples must cover every node"))
    size(values, 2) >= sample_count ||
        throw(ArgumentError("current injection samples must cover every trace sample"))
    return values
end

function _deck_runtime_saturated_transformer_intake(
    parsed::DeckParser.DeckParseResult,
)
    get(parsed.card_counts, :fixed_card_saturated_transformer_intake, 0) > 0 ||
        return nothing
    source_path = DeckParser.deck_source_path(parsed)
    source_path === nothing && return nothing
    isfile(source_path) || return nothing
    intake = DeckParser.parse_saturated_transformer_branch_section_intake_file(source_path)
    assert_valid!(intake.validation)
    return intake
end

function _deck_runtime_initial_voltage_sample(
    parsed::DeckParser.DeckParseResult,
    initial_voltage_source::Symbol,
)
    if initial_voltage_source == :none
        return nothing
    elseif initial_voltage_source == :zero
        return _deck_zero_initial_voltage_sample(parsed)
    elseif initial_voltage_source == :node_conditions
        return _deck_node_voltage_initial_sample(parsed)
    elseif initial_voltage_source == :steady_state
        return _deck_runtime_steady_state_voltage_phasors(parsed)
    elseif initial_voltage_source == :synchronous_machine_terminals
        return _deck_synchronous_machine_terminal_voltage_initial_sample(parsed)
    end
    throw(ArgumentError("unsupported initial voltage source $(initial_voltage_source)"))
end

function _deck_runtime_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult,
)
    saturated_transformer_intake = _deck_runtime_saturated_transformer_intake(parsed)
    saturated_transformer_intake === nothing &&
        return deck_steady_state_voltage_phasors(parsed)
    try
        return deck_steady_state_voltage_phasors(
            parsed;
            saturated_transformer_intake = saturated_transformer_intake,
        )
    catch err
        message = sprint(showerror, err)
        if err isa ArgumentError &&
           occursin("missing saturated transformer winding", message)
            return deck_steady_state_voltage_phasors(parsed)
        end
        rethrow()
    end
end

function _distributed_transposed_line_modal_transform()
    return [
        inv(sqrt(3.0)) inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) -inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) 0.0 -2.0 * inv(sqrt(6.0))
    ]
end

function _steady_state_terminal_phasors(sample, node_indices::AbstractVector{<:Integer})
    return ComplexF64[
        node == 0 ? complex(0.0, 0.0) : sample.node_voltage_phasors[Int(node)]
        for node in node_indices
    ]
end

function _distributed_transposed_line_steady_state_history_waves(
    modal_state,
    sample,
)
    phase_count = modal_state.phase_count
    from_voltage =
        _steady_state_terminal_phasors(sample, modal_state.from_node_indices)
    to_voltage =
        _steady_state_terminal_phasors(sample, modal_state.to_node_indices)
    modal_transform = _distributed_transposed_line_modal_transform()
    modal_voltage_from =
        ComplexF64.(phase_count .* (transpose(modal_transform) * from_voltage))
    modal_voltage_to =
        ComplexF64.(phase_count .* (transpose(modal_transform) * to_voltage))
    coefficients = distributed_transposed_line_norton_coefficients(modal_state)
    angular_frequency = 2.0 * pi * Float64(sample.steady_state_frequency_hz)
    outgoing_wave_from = ComplexF64[]
    outgoing_wave_to = ComplexF64[]
    for mode_index in 1:modal_state.group_length
        admittance = coefficients.modal_companion_admittance_values[mode_index]
        damping = coefficients.modal_history_damping_values[mode_index]
        phase_shift = cis(
            -angular_frequency * (
                modal_state.modal_propagation_times_s[mode_index] -
                Float64(sample.timestep_s)
            ),
        )
        delayed_damping = damping * phase_shift
        denominator = 1.0 - delayed_damping^2
        push!(
            outgoing_wave_from,
            admittance * (
                modal_voltage_to[mode_index] -
                delayed_damping * modal_voltage_from[mode_index]
            ) / denominator,
        )
        push!(
            outgoing_wave_to,
            admittance * (
                modal_voltage_from[mode_index] -
                delayed_damping * modal_voltage_to[mode_index]
            ) / denominator,
        )
    end
    return (
        outgoing_wave_from = outgoing_wave_from,
        outgoing_wave_to = outgoing_wave_to,
    )
end

function _distributed_transposed_line_modal_component_terms(phase_terms, phase_count::Int)
    modal_transform = _distributed_transposed_line_modal_transform()
    term_modal_indices = Int[]
    term_weight_values = Float64[]
    term_from_real_values = Float64[]
    term_from_imag_values = Float64[]
    term_to_real_values = Float64[]
    term_to_imag_values = Float64[]
    for mode_index in 1:phase_count
        for phase_index in 1:phase_count
            push!(term_modal_indices, mode_index)
            push!(term_weight_values, modal_transform[phase_index, mode_index])
            push!(term_from_real_values, phase_terms.terminal_from_real_values[phase_index])
            push!(term_from_imag_values, phase_terms.terminal_from_imag_values[phase_index])
            push!(term_to_real_values, phase_terms.terminal_to_real_values[phase_index])
            push!(term_to_imag_values, phase_terms.terminal_to_imag_values[phase_index])
        end
    end
    return (
        term_modal_indices = term_modal_indices,
        term_weight_values = term_weight_values,
        term_from_real_values = term_from_real_values,
        term_from_imag_values = term_from_imag_values,
        term_to_real_values = term_to_real_values,
        term_to_imag_values = term_to_imag_values,
    )
end

function _distributed_transposed_line_steady_state_initial_history(
    modal_state,
    pi_equivalent,
    sample,
    timestep_s::Real,
    history_storage_start_index::Integer,
    index::Integer,
)
    phase_count = modal_state.phase_count
    from_voltage =
        _steady_state_terminal_phasors(sample, modal_state.from_node_indices)
    to_voltage =
        _steady_state_terminal_phasors(sample, modal_state.to_node_indices)
    series_admittance = inv(pi_equivalent.phase_series_impedance_matrix)
    shunt_admittance = pi_equivalent.phase_shunt_admittance_matrix
    terminal_current_from =
        series_admittance * (from_voltage - to_voltage) +
        shunt_admittance * from_voltage
    terminal_current_to =
        series_admittance * (to_voltage - from_voltage) +
        shunt_admittance * to_voltage
    phase_matrix = distributed_transposed_line_initial_history_phase_matrix(
        modal_state;
        name = Symbol(
            "distributed_transposed_line_initial_history_phase_matrix_",
            index,
        ),
    )
    phase_terms = distributed_transposed_line_initial_history_phase_terms(
        reduced_phase_matrix_upper_values = phase_matrix.phase_matrix_upper_values,
        matrix_terminal_from_real_values = real.(terminal_current_to),
        matrix_terminal_from_imag_values = imag.(terminal_current_to),
        matrix_terminal_to_real_values = real.(terminal_current_from),
        matrix_terminal_to_imag_values = imag.(terminal_current_from),
        initial_terminal_from_real_values = real.(to_voltage),
        initial_terminal_from_imag_values = imag.(to_voltage),
        initial_terminal_to_real_values = real.(from_voltage),
        initial_terminal_to_imag_values = imag.(from_voltage),
        phase_indices = modal_state.phase_indices,
        name = Symbol(
            "distributed_transposed_line_initial_history_phase_terms_",
            index,
        ),
    )
    component_terms =
        _distributed_transposed_line_modal_component_terms(phase_terms, phase_count)
    modal_indices = collect(1:modal_state.group_length)
    components = distributed_transposed_line_initial_history_modal_components(
        modal_indices = modal_indices,
        term_modal_indices = component_terms.term_modal_indices,
        term_weight_values = component_terms.term_weight_values,
        term_from_real_values = component_terms.term_from_real_values,
        term_from_imag_values = component_terms.term_from_imag_values,
        term_to_real_values = component_terms.term_to_real_values,
        term_to_imag_values = component_terms.term_to_imag_values,
        name = Symbol(
            "distributed_transposed_line_initial_history_modal_components_",
            index,
        ),
    )
    phasors = distributed_transposed_line_initial_history_phasors(
        modal_scale_values = phase_matrix.modal_history_scale_values,
        modal_component_from_real_values =
            components.modal_component_from_real_values,
        modal_component_from_imag_values =
            components.modal_component_from_imag_values,
        modal_component_to_real_values =
            components.modal_component_to_real_values,
        modal_component_to_imag_values =
            components.modal_component_to_imag_values,
        modal_indices = modal_indices,
        name = Symbol("distributed_transposed_line_initial_history_phasors_", index),
    )
    norton_coefficients = distributed_transposed_line_norton_coefficients(modal_state)
    delay_counts =
        floor.(Int, modal_state.modal_propagation_times_s ./ Float64(timestep_s))
    seed = distributed_transposed_line_initial_history_seed(
        timestep_s = timestep_s,
        steady_state_frequency_hz = sample.steady_state_frequency_hz,
        modal_delay_counts = delay_counts,
        outgoing_wave_from_amplitudes = phasors.modal_history_from_amplitudes,
        outgoing_wave_from_phases = phasors.modal_history_from_phases,
        outgoing_wave_to_amplitudes = phasors.modal_history_to_amplitudes,
        outgoing_wave_to_phases = phasors.modal_history_to_phases,
        modal_history_interpolation_factors =
            modal_state.modal_propagation_times_s ./ Float64(timestep_s) .-
            delay_counts,
        modal_history_damping_values =
            norton_coefficients.modal_history_damping_values,
        modal_signed_characteristic_values =
            modal_state.modal_signed_characteristic_impedances,
        history_storage_start_index = history_storage_start_index,
        modal_indices = modal_indices,
        name = Symbol("distributed_transposed_line_initial_history_seed_", index),
    )
    history = distributed_transposed_line_history_state(
        modal_state,
        seed;
        timestep_s = timestep_s,
        steady_state_frequency_hz = sample.steady_state_frequency_hz,
        initialized_from_steady_state = true,
        name = Symbol("distributed_transposed_line_history_state_", index),
    )
    return (
        phase_matrix = phase_matrix,
        phase_terms = phase_terms,
        modal_components = components,
        phasors = phasors,
        seed = seed,
        history = history,
        terminal_current_from = terminal_current_from,
        terminal_current_to = terminal_current_to,
    )
end

function _shared_distributed_transposed_line_history_states(
    histories::AbstractVector{DistributedTransposedLineHistoryState},
)
    isempty(histories) && return DistributedTransposedLineHistoryState[]
    global_start = minimum(history.history_storage_start_index for history in histories)
    global_end = maximum(maximum(history.storage_end_indices) for history in histories)
    shared_from = zeros(Float64, global_end - global_start + 1)
    shared_to = zeros(Float64, global_end - global_start + 1)

    for history in histories
        for mode_index in eachindex(history.storage_start_indices)
            start_slot = history.storage_start_indices[mode_index]
            first_slot = start_slot == global_start ? start_slot : start_slot + 1
            last_slot = history.storage_end_indices[mode_index]
            first_slot <= last_slot || continue
            for slot in first_slot:last_slot
                source_offset = slot - history.history_storage_start_index + 1
                1 <= source_offset <= length(history.packed_history_from_values) ||
                    throw(ArgumentError("distributed-line packed history source slot is outside the local storage vector"))
                shared_offset = slot - global_start + 1
                shared_from[shared_offset] = history.packed_history_from_values[source_offset]
                shared_to[shared_offset] = history.packed_history_to_values[source_offset]
            end
        end
    end

    shared_histories = DistributedTransposedLineHistoryState[]
    for history in histories
        push!(
            shared_histories,
            DistributedTransposedLineHistoryState(
                history.name,
                copy(history.phase_indices),
                copy(history.line_numbers),
                copy(history.modal_sequence_indices),
                history.timestep_s,
                history.steady_state_frequency_hz,
                history.angular_step_rad,
                global_start,
                global_end,
                copy(history.storage_start_indices),
                copy(history.storage_end_indices),
                copy(history.storage_lengths),
                copy(history.history_sample_counts),
                copy(history.modal_history_interpolation_factors),
                copy(history.history_read_indices),
                copy(history.outgoing_wave_from_real_values),
                copy(history.outgoing_wave_from_imag_values),
                copy(history.outgoing_wave_to_real_values),
                copy(history.outgoing_wave_to_imag_values),
                [copy(values) for values in history.modal_history_from_values],
                [copy(values) for values in history.modal_history_to_values],
                shared_from,
                shared_to,
                history.initialized_from_steady_state,
                history.phase_count,
            ),
        )
    end
    return shared_histories
end

function _deck_distributed_transposed_line_config(
    parsed::DeckParser.DeckParseResult;
    steady_state_initial_sample = nothing,
)
    modal_states = DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)
    isempty(modal_states) && return nothing
    pi_equivalents =
        DeckParser.deck_distributed_transposed_line_steady_state_pi_equivalents(parsed)
    if steady_state_initial_sample !== nothing &&
       length(pi_equivalents) != length(modal_states)
        throw(ArgumentError("steady-state distributed-line initialization requires a PI equivalent for each modal line state"))
    end
    history_states =
        steady_state_initial_sample === nothing ?
        DeckParser.deck_distributed_transposed_line_history_states(
            parsed;
            first_history_storage_index = 1,
        ) :
        begin
            options = DeckParser.deck_fixed_time_horizon_options(parsed)
            storage_index = 1
            histories = DistributedTransposedLineHistoryState[]
            initial_history_results = Any[]
            for index in eachindex(modal_states)
                initial_history =
                    _distributed_transposed_line_steady_state_initial_history(
                    modal_states[index],
                    pi_equivalents[index],
                    steady_state_initial_sample,
                    options.dt_s,
                    storage_index,
                    index,
                )
                history = initial_history.history
                push!(histories, history)
                push!(initial_history_results, initial_history)
                storage_index = history.next_history_storage_index
            end
            _shared_distributed_transposed_line_history_states(histories)
        end
    length(history_states) == length(modal_states) ||
        throw(ArgumentError("distributed transposed line modal/history state counts must match"))
    return (
        enabled = true,
        modal_states = modal_states,
        history_states = history_states,
        steady_state_initial_history_results =
            steady_state_initial_sample === nothing ? Any[] : initial_history_results,
        current_injection_values = Float64[],
    )
end

_deck_has_dynamic_distributed_line(parsed::DeckParser.DeckParseResult) =
    !isempty(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed))

_deck_has_frequency_dependent_line_runtime(parsed::DeckParser.DeckParseResult) =
    !isempty(DeckParser.deck_sampled_frequency_line_rows(parsed)) ||
    !isempty(DeckParser.deck_semlyen_line_groups(parsed)) ||
    !isempty(DeckParser.deck_rational_frequency_line_groups(parsed))

function _primary_dynamic_timestep_state(
    node_count::Integer,
    dt_s::Real,
    switch_count::Integer,
    ;
    source_card_read_requested::Bool=true,
)
    node_total = Int(node_count)
    node_total > 0 || throw(ArgumentError("deck dynamic timestep state requires at least one node"))
    switches = Int(switch_count)
    switches >= 0 || throw(ArgumentError("deck switch count must be nonnegative"))
    shifted_node_total = node_total + 1
    return OVER16AcceptedTimestepState(
        OVER16SourceUpdateState(
            fill(0.0, 10),
            zeros(node_total),
            zeros(node_total);
            iread = source_card_read_requested ? 1 : 0,
        ),
        OVER16TACSUtilityState(zeros(1)),
        OVER16OutputReportState(peaknd_values = [0.0, 0.0, 0.0]),
        OVER16PostExtremaControlState(0.0, Float64(dt_s), 0),
        OVER16SwitchScanState(fill(-1, switches), zeros(switches)),
        OVER16SwitchOperationState(Int[], 0),
        OVER16SwitchTopologyState(
            fill(false, switches);
            nextsw = zeros(Int, switches),
            kode = zeros(Int, shifted_node_total),
        ),
        OVER16SwitchCurrentState(zeros(shifted_node_total), zeros(switches)),
        OVER16SwitchPostCurrentState(fill(-1, switches), zeros(switches), zeros(switches)),
        OVER16SwitchBValueExportState(Float64[]),
        OVER16SwitchAlterationState(0, 0, 0, 0),
    )
end

function _deck_dynamic_timestep_state(
    parsed::DeckParser.DeckParseResult,
    node_count::Integer,
    dt_s::Real,
)
    return _primary_dynamic_timestep_state(
        node_count,
        dt_s,
        length(DeckParser.deck_time_switch_names(parsed));
        source_card_read_requested =
            !isempty(DeckParser.deck_over16_source_card_rows(parsed)),
    )
end

function run_nested_cable_primary_timestep!(
    context::EMTStepContext,
    line_state::NestedCableTransientLineState;
    from_node_indices,
    to_node_indices,
    target_frequency_hz::Real,
    current_projection::Symbol = :harmonic_phasor,
    current_projection_angle_rad::Real = 0.0,
    over16_state::Union{Nothing,OVER16AcceptedTimestepState} = nothing,
)
    abs(context.dt_s - line_state.dt_s) <=
        64.0 * eps(Float64) * max(context.dt_s, line_state.dt_s) ||
        throw(ArgumentError("nested cable state dt_s must match the primary timestep context"))
    line_state.physical_checks_passed ||
        throw(ArgumentError("nested cable transient state failed its physical checks"))
    state = over16_state === nothing ?
        _primary_dynamic_timestep_state(context.system.node_count, context.dt_s, 0) :
        over16_state
    current_values = Float64[]
    return step_with_over16_boundary!(
        context,
        state;
        frequency_dependent_line_config = (
            enabled = true,
            nested_cable_states = [line_state],
            from_node_indices = from_node_indices,
            to_node_indices = to_node_indices,
            target_frequency_hz = Float64(target_frequency_hz),
            current_projection = current_projection,
            current_projection_angle_rad = Float64(current_projection_angle_rad),
            current_injection_values = current_values,
        ),
    )
end

function _terminal_voltage_values(
    voltage::AbstractVector{Float64},
    node_indices::AbstractVector{Int},
)
    values = Float64[]
    for node in node_indices
        if node == 0
            push!(values, 0.0)
        else
            1 <= node <= length(voltage) ||
                throw(ArgumentError("distributed transposed line node index exceeds solved voltage length"))
            push!(values, voltage[node])
        end
    end
    return values
end

function _frequency_dependent_line_recursive_states(config::NamedTuple)
    if haskey(config, :nested_cable_states)
        cable_states = getproperty(config, :nested_cable_states)
        all(state -> state isa NestedCableTransientLineState, cable_states) ||
            throw(ArgumentError("nested_cable_states must contain NestedCableTransientLineState values"))
        return [state.recursive_state for state in cable_states]
    elseif haskey(config, :recursive_convolution_states)
        return getproperty(config, :recursive_convolution_states)
    elseif haskey(config, :states)
        return getproperty(config, :states)
    end
    throw(ArgumentError(
        "frequency_dependent_line_config requires nested_cable_states or recursive_convolution_states",
    ))
end

function _frequency_dependent_line_projection_tolerance(config::NamedTuple)
    tolerance = haskey(config, :real_current_projection_tolerance) ?
        Float64(getproperty(config, :real_current_projection_tolerance)) :
        1.0e-9
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("real_current_projection_tolerance must be finite and nonnegative"))
    return tolerance
end

function _frequency_dependent_line_current_projection(config::NamedTuple)
    kind = Symbol(get(config, :current_projection, :strict_real))
    kind in (:strict_real, :harmonic_phasor) ||
        throw(ArgumentError("current_projection must be :strict_real or :harmonic_phasor"))
    angle = Float64(get(config, :current_projection_angle_rad, 0.0))
    isfinite(angle) || throw(ArgumentError("current_projection_angle_rad must be finite"))
    return kind, angle
end

function _frequency_dependent_line_node_group(
    values,
    index::Int,
    line_count::Int,
    phase_count::Int,
    node_count::Int,
    label::AbstractString,
)
    raw_nodes =
        values isa AbstractVector{<:Integer} ?
        begin
            line_count == 1 ||
                throw(ArgumentError("$label must contain one node-index vector per frequency-dependent line"))
            values
        end :
        begin
            length(values) == line_count ||
                throw(ArgumentError("$label count must match frequency-dependent line count"))
            values[index]
        end
    nodes = Int.(collect(raw_nodes))
    length(nodes) == phase_count ||
        throw(ArgumentError("$label phase count must match frequency-dependent line state"))
    for node in nodes
        0 <= node <= node_count ||
            throw(ArgumentError("$label contains a node index outside the solved network"))
    end
    return nodes
end

function _frequency_dependent_line_target_frequency(
    config::NamedTuple,
    index::Int,
    line_count::Int,
)
    if haskey(config, :target_frequency_hz_values)
        values = getproperty(config, :target_frequency_hz_values)
        length(values) == line_count ||
            throw(ArgumentError("target_frequency_hz_values count must match frequency-dependent line count"))
        frequency = Float64(values[index])
    elseif haskey(config, :frequency_hz_values)
        values = getproperty(config, :frequency_hz_values)
        length(values) == line_count ||
            throw(ArgumentError("frequency_hz_values count must match frequency-dependent line count"))
        frequency = Float64(values[index])
    elseif haskey(config, :target_frequency_hz)
        frequency = Float64(getproperty(config, :target_frequency_hz))
    elseif haskey(config, :frequency_hz)
        frequency = Float64(getproperty(config, :frequency_hz))
    else
        throw(ArgumentError("frequency_dependent_line_config requires target_frequency_hz or frequency_hz"))
    end
    isfinite(frequency) && frequency >= 0.0 ||
        throw(ArgumentError("frequency-dependent line target frequency must be finite and nonnegative"))
    return frequency
end

function _frequency_dependent_line_real_currents(
    values::AbstractVector,
    tolerance::Float64,
    label::AbstractString,
    projection_kind::Symbol,
    projection_angle_rad::Float64,
)
    currents = Float64[]
    sizehint!(currents, length(values))
    max_imag = 0.0
    rotation = cis(projection_angle_rad)
    for (index, value) in pairs(values)
        current = ComplexF64(value) * rotation
        imag_abs = abs(imag(current))
        max_imag = max(max_imag, imag_abs)
        (projection_kind == :harmonic_phasor || imag_abs <= tolerance) ||
            throw(ArgumentError("$label current $index has non-real residual $imag_abs"))
        push!(currents, real(current))
    end
    return currents, max_imag
end

function _frequency_dependent_line_phase_voltages(
    voltage::AbstractVector{Float64},
    nodes::AbstractVector{Int},
)
    return ComplexF64.(_terminal_voltage_values(voltage, nodes))
end

function _frequency_dependent_line_deck_current_injection!(
    config,
    rhs::AbstractVector{Float64},
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("frequency_dependent_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    states = _frequency_dependent_line_recursive_states(config)
    line_count = length(states)
    line_count > 0 ||
        throw(ArgumentError("frequency_dependent_line_config requires at least one recursive state"))
    haskey(config, :from_node_indices) ||
        throw(ArgumentError("frequency_dependent_line_config requires from_node_indices"))
    haskey(config, :to_node_indices) ||
        throw(ArgumentError("frequency_dependent_line_config requires to_node_indices"))
    tolerance = _frequency_dependent_line_projection_tolerance(config)
    projection_kind, projection_angle =
        _frequency_dependent_line_current_projection(config)

    rhs_before = copy(rhs)
    rhs_update_count = 0
    sending_currents = Vector{Vector{Float64}}()
    receiving_currents = Vector{Vector{Float64}}()
    from_node_groups = Vector{Vector{Int}}()
    to_node_groups = Vector{Vector{Int}}()
    max_imag = 0.0
    for index in 1:line_count
        state = states[index]
        phase_count = length(state.sending_phase_current)
        length(state.receiving_phase_current) == phase_count ||
            throw(ArgumentError("frequency-dependent line sending/receiving current counts must match"))
        from_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :from_node_indices),
            index,
            line_count,
            phase_count,
            length(rhs),
            "frequency-dependent from_node_indices",
        )
        to_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :to_node_indices),
            index,
            line_count,
            phase_count,
            length(rhs),
            "frequency-dependent to_node_indices",
        )
        sending, sending_imag = _frequency_dependent_line_real_currents(
            state.sending_phase_current,
            tolerance,
            "frequency-dependent sending",
            projection_kind,
            projection_angle,
        )
        receiving, receiving_imag = _frequency_dependent_line_real_currents(
            state.receiving_phase_current,
            tolerance,
            "frequency-dependent receiving",
            projection_kind,
            projection_angle,
        )
        max_imag = max(max_imag, sending_imag, receiving_imag)
        for phase_index in 1:phase_count
            from_node = from_nodes[phase_index]
            if from_node != 0
                rhs[from_node] += sending[phase_index]
                rhs_update_count += 1
            end
            to_node = to_nodes[phase_index]
            if to_node != 0
                rhs[to_node] += receiving[phase_index]
                rhs_update_count += 1
            end
        end
        push!(sending_currents, sending)
        push!(receiving_currents, receiving)
        push!(from_node_groups, from_nodes)
        push!(to_node_groups, to_nodes)
    end
    if haskey(config, :current_injection_values)
        values = getproperty(config, :current_injection_values)
        values isa Vector{Float64} ||
            throw(ArgumentError("frequency-dependent current_injection_values must be Vector{Float64}"))
        resize!(values, length(rhs))
        for index in eachindex(rhs)
            values[index] = rhs[index] - rhs_before[index]
        end
    end
    return (
        source = :frequency_dependent_line_deck_current_injection,
        outcome = :history_current_injection,
        line_update_count = 0,
        rhs_update_count = rhs_update_count,
        rhs_before_values = rhs_before,
        rhs_after_values = copy(rhs),
        sending_phase_current_values = sending_currents,
        receiving_phase_current_values = receiving_currents,
        from_node_indices = from_node_groups,
        to_node_indices = to_node_groups,
        real_current_projection_max_imag_abs = max_imag,
        current_projection = projection_kind,
        current_projection_angle_rad = projection_angle,
        nested_cable_frequency_state_consumed = haskey(config, :nested_cable_states),
        state_mutated = false,
    )
end

function _frequency_dependent_line_deck_history_update!(
    config,
    voltage::AbstractVector{Float64},
    current_result,
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("frequency_dependent_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    states = _frequency_dependent_line_recursive_states(config)
    line_count = length(states)
    line_count > 0 ||
        throw(ArgumentError("frequency_dependent_line_config requires at least one recursive state"))
    updates = Any[]
    target_frequencies = Float64[]
    for index in 1:line_count
        state = states[index]
        phase_count = length(state.sending_phase_current)
        from_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :from_node_indices),
            index,
            line_count,
            phase_count,
            length(voltage),
            "frequency-dependent from_node_indices",
        )
        to_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :to_node_indices),
            index,
            line_count,
            phase_count,
            length(voltage),
            "frequency-dependent to_node_indices",
        )
        target_frequency = _frequency_dependent_line_target_frequency(config, index, line_count)
        update = frequency_dependent_line_recursive_convolution_update!(
            state,
            target_frequency,
            _frequency_dependent_line_phase_voltages(voltage, from_nodes),
            _frequency_dependent_line_phase_voltages(voltage, to_nodes),
        )
        push!(updates, update)
        push!(target_frequencies, target_frequency)
    end
    recursive_runtime_executed =
        any(update -> get(update, :bounded_recursive_convolution_runtime_executed, false), updates)
    skin_effect_internal_impedance_executed =
        any(update -> get(update, :skin_effect_internal_impedance_executed, false), updates)
    earth_return_impedance_executed =
        any(update -> get(update, :earth_return_impedance_executed, false), updates)
    frequency_dependent_fitting_executed =
        any(update -> get(update, :frequency_dependent_fitting_executed, false), updates)
    frequency_loop_executed =
        any(update -> get(update, :frequency_loop_executed, false), updates)
    pipe_sheath_side_effects_executed =
        any(update -> get(update, :pipe_sheath_side_effects_executed, false), updates)
    deferred_effects = Symbol[]
    frequency_dependent_fitting_executed ||
        push!(deferred_effects, :full_bpa_frequency_dependent_fitting)
    pipe_sheath_side_effects_executed ||
        push!(deferred_effects, :full_over47_pipe_or_sheath_coupling)
    frequency_loop_executed ||
        push!(deferred_effects, :full_over47_frequency_loop_side_effects)
    haskey(config, :nested_cable_states) ||
        push!(deferred_effects, :frequency_dependent_line_deck_grammar)
    return (
        source = :frequency_dependent_line_deck_timestep_update,
        outcome = :timestep_integration,
        line_update_count = length(updates),
        rhs_update_count = current_result === nothing ? 0 : current_result.rhs_update_count,
        recursive_convolution_updates = updates,
        phase_current_injections =
            current_result === nothing ? Any[] : [current_result],
        rhs_before_values =
            current_result === nothing ? Float64[] : current_result.rhs_before_values,
        rhs_after_values =
            current_result === nothing ? Float64[] : current_result.rhs_after_values,
        target_frequency_hz_values = target_frequencies,
        frequency_dependent_line_runtime_executed = !isempty(updates),
        recursive_convolution_runtime_executed =
            recursive_runtime_executed,
        skin_effect_internal_impedance_executed =
            skin_effect_internal_impedance_executed,
        earth_return_impedance_executed =
            earth_return_impedance_executed,
        frequency_dependent_fitting_executed =
            frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed =
            pipe_sheath_side_effects_executed,
        nested_cable_frequency_state_consumed = haskey(config, :nested_cable_states),
        state_mutated = !isempty(updates),
        deferred_effects = Tuple(deferred_effects),
    )
end

function _distributed_transposed_line_deck_current_injection!(
    config,
    rhs::AbstractVector{Float64},
    ;
    collect_diagnostics::Bool = true,
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("distributed_transposed_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    haskey(config, :modal_states) ||
        throw(ArgumentError("distributed_transposed_line_config requires modal_states"))
    haskey(config, :history_states) ||
        throw(ArgumentError("distributed_transposed_line_config requires history_states"))
    modal_states = config.modal_states
    history_states = config.history_states
    length(modal_states) == length(history_states) ||
        throw(ArgumentError("distributed transposed line modal/history state counts must match"))
    injections = collect_diagnostics ? Any[] : nothing
    rhs_before = copy(rhs)
    rhs_update_count = 0
    for index in eachindex(modal_states)
        modal_state = modal_states[index]
        history_state = history_states[index]
        injection = distributed_transposed_line_history_current_injection!(
            rhs,
            modal_state,
            history_state;
            name = collect_diagnostics ?
                Symbol("distributed_transposed_line_history_current_injection_", index) :
                history_state.name,
            collect_diagnostics = collect_diagnostics,
        )
        if collect_diagnostics
            rhs_update_count += injection.rhs_update_count
            push!(injections, injection)
        end
    end
    if haskey(config, :current_injection_values)
        current_injection_values = getproperty(config, :current_injection_values)
        current_injection_values isa Vector{Float64} ||
            throw(ArgumentError("distributed line current_injection_values must be Vector{Float64}"))
        resize!(current_injection_values, length(rhs))
        for index in eachindex(rhs)
            current_injection_values[index] = rhs[index] - rhs_before[index]
        end
    end
    if !collect_diagnostics
        return (
            rhs_before_values = rhs_before,
            rhs_after_values = rhs,
        )
    end
    return (
        source = :distributed_transposed_line_deck_current_injection,
        outcome = :history_current_injection,
        line_update_count = 0,
        rhs_update_count = rhs_update_count,
        modal_timestep_updates = Any[],
        phase_current_injections = injections,
        rhs_before_values = rhs_before,
        rhs_after_values = copy(rhs),
        state_mutated = false,
    )
end

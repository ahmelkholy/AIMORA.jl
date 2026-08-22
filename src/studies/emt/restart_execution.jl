function _restart_symbol_equal(left::Symbol, right::Symbol)
    normalize(name::Symbol) = let value = uppercase(String(name))
        value in ("0", "GROUND", "REF", "REFERENCE") ? "GROUND" : value
    end
    return normalize(left) == normalize(right)
end

function _restart_trace_matrix(
    values::Matrix{Float64},
    appended_column_count::Int,
)
    appended_column_count >= 0 || throw(ArgumentError(
        "appended trace column count must be nonnegative",
    ))
    result = Matrix{Float64}(
        undef,
        size(values, 1),
        size(values, 2) + appended_column_count,
    )
    isempty(values) || (result[:, 1:size(values, 2)] .= values)
    appended_column_count == 0 || fill!(
        @view(result[:, size(values, 2) + 1:end]),
        0.0,
    )
    return result
end

function _extend_emt_restart_horizon!(
    context::EMTStepContext,
    additional_time_s::Float64;
    recorded_step_indices = nothing,
)
    additional_time_s > 0.0 && isfinite(additional_time_s) || throw(ArgumentError(
        "restart additional_time_s must be finite and positive",
    ))
    context.step_index == context.step_count + 1 || throw(ArgumentError(
        "restart requires a completed EMT checkpoint",
    ))
    context.trace_write_index == length(context.recorded_step_indices) + 1 ||
        throw(ArgumentError("restart checkpoint has incomplete trace storage"))
    !isempty(context.recorded_step_indices) &&
        last(context.recorded_step_indices) == context.step_count ||
        throw(ArgumentError("restart checkpoint must record its final timestep"))

    appended_step_count = fixed_step_count(context.dt_s, additional_time_s)
    appended_step_count > 0 || throw(ArgumentError(
        "restart continuation must append at least one timestep",
    ))
    checkpoint_step = context.step_count
    final_step = checkpoint_step + appended_step_count
    future_steps = recorded_step_indices === nothing ?
        collect(checkpoint_step + 1:final_step) :
        Int[Int(step) for step in recorded_step_indices]
    isempty(future_steps) && throw(ArgumentError(
        "restart recorded_step_indices must not be empty",
    ))
    previous_step = checkpoint_step
    for step in future_steps
        checkpoint_step < step <= final_step || throw(ArgumentError(
            "restart recorded step $step is outside $(checkpoint_step + 1):$final_step",
        ))
        step > previous_step || throw(ArgumentError(
            "restart recorded_step_indices must be strictly increasing",
        ))
        previous_step = step
    end
    last(future_steps) == final_step || throw(ArgumentError(
        "restart continuation must record its final timestep",
    ))

    appended_column_count = length(future_steps)
    append!(context.recorded_step_indices, future_steps)
    append!(context.time_s, zeros(Float64, appended_column_count))
    context.voltage_pu = _restart_trace_matrix(
        context.voltage_pu,
        appended_column_count,
    )
    context.output_pu = _restart_trace_matrix(
        context.output_pu,
        appended_column_count,
    )
    context.step_count = final_step
    context.t_end_s = final_step * context.dt_s
    context.t_s = context.step_index * context.dt_s
    return appended_step_count
end

function _restart_switch_index(
    plan::DeckOVER16BoundaryPlan,
    mutation::DeckParser.DeckRestartSwitchMutation,
)
    if mutation.selector_kind == :index
        index = Int(mutation.switch_index)
        1 <= index <= plan.switch_count || throw(ArgumentError(
            "restart switch index $index is outside 1:$(plan.switch_count)",
        ))
        return index
    elseif mutation.selector_kind == :endpoints
        matches = Int[]
        for index in 1:plan.switch_count
            direct = _restart_symbol_equal(
                plan.switch_from_node_names[index],
                mutation.from_node,
            ) && _restart_symbol_equal(
                plan.switch_to_node_names[index],
                mutation.to_node,
            )
            reverse = _restart_symbol_equal(
                plan.switch_from_node_names[index],
                mutation.to_node,
            ) && _restart_symbol_equal(
                plan.switch_to_node_names[index],
                mutation.from_node,
            )
            (direct || reverse) && push!(matches, index)
        end
        length(matches) == 1 || throw(ArgumentError(
            "restart switch endpoints $(mutation.from_node)-$(mutation.to_node) matched $(length(matches)) switches",
        ))
        return only(matches)
    end
    throw(ArgumentError("unsupported restart switch selector $(mutation.selector_kind)"))
end

function _restart_replace_time_switch!(
    context::EMTStepContext,
    switch_name::Symbol,
    close_time_s::Float64,
    open_time_s::Float64,
    open_delay_time_s::Float64,
    critical_current_a::Float64,
)
    element_index = findfirst(
        name -> _restart_symbol_equal(name, switch_name),
        context.element_names,
    )
    element_index === nothing && throw(ArgumentError(
        "restart switch $switch_name is missing from the timestep system",
    ))
    element = context.system.elements[element_index]
    if element isa TimeSwitch
        previously_closed = switch_closed(element, context.t_s)
        element.close_time_s = close_time_s
        element.open_time_s = open_time_s
        configure_current_extinction!(
            element,
            open_delay_time_s,
            critical_current_a,
            context.t_s,
            currently_closed = previously_closed,
        )
    elseif element isa CurrentZeroSwitch
        element.close_time_s = close_time_s
        element.open_request_time_s = open_time_s
        element.open_delay_time_s = open_delay_time_s
        element.critical_current_a = critical_current_a
        scheduled_closed =
            context.t_s >= close_time_s && context.t_s < open_time_s
        if scheduled_closed && !element.closed
            element.closed = true
            element.opened = false
            element.current_initialized = false
            element.open_reason = :restart_topology_change
            element.operation_count += 1
        end
    else
        throw(ArgumentError(
            "restart switch $switch_name does not own a mutable timed-switch model",
        ))
    end
    return nothing
end

function _restart_over5_switch_index(
    plan::DeckOVER16BoundaryPlan,
    switch_name::Symbol,
)
    return findfirst(
        name -> _restart_symbol_equal(name, switch_name),
        plan.over5_switch_names,
    )
end

function _apply_restart_switch_mutation!(
    runtime::PreparedDynamicDeckRuntime,
    mutation::DeckParser.DeckRestartSwitchMutation,
)
    plan = runtime.plan
    context = runtime.context
    index = _restart_switch_index(plan, mutation)
    previous_close = plan.switch_close_time_s_values[index]
    previous_open = plan.switch_open_time_s_values[index]
    element_index = findfirst(
        name -> _restart_symbol_equal(name, plan.switch_names[index]),
        context.element_names,
    )
    element_index === nothing && throw(ArgumentError(
        "restart switch $(plan.switch_names[index]) is missing from the timestep context",
    ))
    element = context.system.elements[element_index]
    previous_delay =
        element isa CurrentZeroSwitch ? element.open_delay_time_s :
        element isa TimeSwitch && element.current_extinction !== nothing ?
            element.current_extinction.not_before_time_s : 0.0
    previous_critical =
        element isa CurrentZeroSwitch ? element.critical_current_a :
        element isa TimeSwitch && element.current_extinction !== nothing ?
            element.current_extinction.critical_current_a : 0.0
    applied_close = mutation.close_time_s === missing ?
        previous_close : Float64(mutation.close_time_s)
    applied_open = mutation.open_time_s === missing ?
        previous_open : Float64(mutation.open_time_s)
    applied_delay = mutation.open_delay_s === missing ?
        previous_delay : Float64(mutation.open_delay_s)
    applied_critical = mutation.critical_current === missing ?
        previous_critical : Float64(mutation.critical_current)
    applied_open > applied_close || throw(ArgumentError(
        "restart switch $(plan.switch_names[index]) open time must follow its close time",
    ))
    applied_delay >= 0.0 && (isfinite(applied_delay) || applied_delay == Inf) ||
        throw(ArgumentError("restart switch delay time must be nonnegative"))
    applied_critical >= 0.0 && isfinite(applied_critical) ||
        throw(ArgumentError("restart switch critical current must be finite and nonnegative"))

    plan.switch_close_time_s_values[index] = applied_close
    plan.switch_open_time_s_values[index] = applied_open
    context.deck_time_switch_close_time_s_values[index] = applied_close
    context.deck_time_switch_open_time_s_values[index] = applied_open
    switch_name = plan.switch_names[index]
    over5_index = _restart_over5_switch_index(plan, switch_name)
    if over5_index !== nothing
        plan.over5_switch_close_time_s_values[over5_index] = applied_close
        plan.over5_switch_open_time_s_values[over5_index] = applied_open
    end
    _restart_replace_time_switch!(
        context,
        switch_name,
        applied_close,
        applied_open,
        applied_delay,
        applied_critical,
    )
    extinction_parameters =
        applied_delay != 0.0 ||
        applied_critical != 0.0 ||
        previous_delay != 0.0 ||
        previous_critical != 0.0
    return EMTRestartMutationRecord(
        extinction_parameters ? :current_extinction_switch : :time_switch,
        switch_name,
        mutation.line_no,
        extinction_parameters ?
            [:close_time_s, :open_time_s, :open_delay_time_s, :critical_current_a] :
            [:close_time_s, :open_time_s],
        extinction_parameters ?
            [previous_close, previous_open, previous_delay, previous_critical] :
            [previous_close, previous_open],
        extinction_parameters ?
            [applied_close, applied_open, applied_delay, applied_critical] :
            [applied_close, applied_open],
    )
end

function _restart_control_source_row(
    parsed::DeckParser.DeckParseResult,
    name::Symbol,
)
    rows = [
        row for row in DeckParser.deck_control_system_source_rows(parsed)
        if _restart_symbol_equal(row.name, name)
    ]
    length(rows) == 1 || throw(ArgumentError(
        "restart control source $name matched $(length(rows)) accepted source rows",
    ))
    return only(rows)
end

function _restart_source_configuration(
    runtime::ControlSystemNetworkRuntime,
    row::DeckParser.DeckControlSystemSourceRow,
)
    window_index = findfirst(
        source -> _restart_symbol_equal(source.name, row.name),
        runtime.windowed_constant_sources,
    )
    sinusoidal_index = findfirst(
        source -> _restart_symbol_equal(source.name, row.name),
        runtime.sinusoidal_sources,
    )
    if window_index !== nothing
        source = runtime.windowed_constant_sources[window_index]
        return (
            source_type = 11,
            amplitude = source.value,
            frequency_or_delay = 0.0,
            phase_or_width = 0.0,
            start_time_s = source.start_time_s,
            stop_time_s = source.stop_time_s,
        )
    elseif sinusoidal_index !== nothing
        source = runtime.sinusoidal_sources[sinusoidal_index]
        return (
            source_type = 14,
            amplitude = source.amplitude,
            frequency_or_delay = source.frequency_hz,
            phase_or_width = rad2deg(source.phase_rad),
            start_time_s = source.start_time_s,
            stop_time_s = source.stop_time_s,
        )
    end
    return (
        source_type = row.source_type,
        amplitude = get(runtime.state.values, row.name, 0.0),
        frequency_or_delay = length(row.numeric_values) >= 2 ? row.numeric_values[2] : 0.0,
        phase_or_width = length(row.numeric_values) >= 3 ? row.numeric_values[3] : 0.0,
        start_time_s = row.activation_start_time_s,
        stop_time_s = row.activation_stop_time_s,
    )
end

function _restart_delete_source_configuration!(
    runtime::ControlSystemNetworkRuntime,
    name::Symbol,
)
    filter!(
        source -> !_restart_symbol_equal(source.name, name),
        runtime.windowed_constant_sources,
    )
    filter!(
        source -> !_restart_symbol_equal(source.name, name),
        runtime.sinusoidal_sources,
    )
    return nothing
end

function _restart_update_network_source_binding!(
    runtime::ControlSystemNetworkRuntime,
    parsed::DeckParser.DeckParseResult,
    name::Symbol,
    target_type::Int,
)
    indices = findall(
        source_name -> _restart_symbol_equal(source_name, name),
        runtime.network_voltage_source_names,
    )
    for index in Iterators.reverse(indices)
        deleteat!(runtime.network_voltage_source_names, index)
        deleteat!(runtime.network_voltage_source_node_indices, index)
    end
    if target_type == 90
        node_index = get(parsed.node_map, name, 0)
        node_index > 0 || throw(ArgumentError(
            "restart type-90 control source $name has no matching network node",
        ))
        push!(runtime.network_voltage_source_names, name)
        push!(runtime.network_voltage_source_node_indices, node_index)
    end
    return nothing
end

_restart_applied_value(value, previous::Float64) =
    value === missing ? previous : Float64(value)

function _apply_restart_control_source_mutation!(
    runtime::PreparedDynamicDeckRuntime,
    parsed::DeckParser.DeckParseResult,
    mutation::DeckParser.DeckRestartControlSourceMutation,
)
    control_runtime = runtime.context.control_system_runtime
    control_runtime === nothing && throw(ArgumentError(
        "restart control-source mutation requires an executable control system",
    ))
    row = _restart_control_source_row(parsed, mutation.name)
    previous = _restart_source_configuration(control_runtime, row)
    source_type = mutation.source_type === missing ?
        previous.source_type : Int(mutation.source_type)
    source_type in (11, 14, 90) || source_type == previous.source_type ||
        throw(ArgumentError(
            "restart control-source type change to $source_type has no executable Julia signal owner",
        ))
    amplitude = _restart_applied_value(mutation.amplitude, previous.amplitude)
    frequency_or_delay = _restart_applied_value(
        mutation.frequency_or_delay,
        previous.frequency_or_delay,
    )
    phase_or_width = _restart_applied_value(
        mutation.phase_or_width,
        previous.phase_or_width,
    )
    start_time_s = _restart_applied_value(
        mutation.activation_start_time_s,
        previous.start_time_s,
    )
    stop_time_s = _restart_applied_value(
        mutation.activation_stop_time_s,
        previous.stop_time_s,
    )
    stop_time_s == 0.0 && (stop_time_s = Inf)
    isfinite(amplitude) || throw(ArgumentError(
        "restart control-source amplitude must be finite",
    ))
    isfinite(start_time_s) || throw(ArgumentError(
        "restart control-source start time must be finite",
    ))
    (isfinite(stop_time_s) || stop_time_s == Inf) && stop_time_s > start_time_s ||
        throw(ArgumentError("restart control-source stop time must follow its start time"))

    _restart_delete_source_configuration!(control_runtime, row.name)
    _restart_update_network_source_binding!(
        control_runtime,
        parsed,
        row.name,
        source_type,
    )
    if source_type == 11
        push!(
            control_runtime.windowed_constant_sources,
            WindowedConstantControlSignal(
                row.name,
                amplitude,
                start_time_s,
                stop_time_s,
            ),
        )
        control_runtime.state.values[row.name] =
            _windowed_constant_control_signal_value(
                last(control_runtime.windowed_constant_sources),
                runtime.context.t_s,
            )
    elseif source_type == 14
        frequency_or_delay > 0.0 && isfinite(frequency_or_delay) ||
            throw(ArgumentError("restart sinusoidal source frequency must be positive"))
        push!(
            control_runtime.sinusoidal_sources,
            SinusoidalControlSignal(
                row.name,
                amplitude,
                frequency_or_delay,
                deg2rad(phase_or_width);
                start_time_s = start_time_s,
                stop_time_s = stop_time_s,
            ),
        )
        control_runtime.state.values[row.name] =
            sinusoidal_control_signal_value(
                last(control_runtime.sinusoidal_sources),
                runtime.context.t_s,
            )
    elseif source_type != 90
        control_runtime.state.values[row.name] = amplitude
    end

    return EMTRestartMutationRecord(
        :control_source,
        row.name,
        mutation.line_no,
        [
            :source_type,
            :amplitude,
            :frequency_or_delay,
            :phase_or_width,
            :activation_start_time_s,
            :activation_stop_time_s,
        ],
        [
            previous.source_type,
            previous.amplitude,
            previous.frequency_or_delay,
            previous.phase_or_width,
            previous.start_time_s,
            previous.stop_time_s,
        ],
        [
            source_type,
            amplitude,
            frequency_or_delay,
            phase_or_width,
            start_time_s,
            stop_time_s,
        ],
    )
end

function _restart_final_kcl_error(runtime::PreparedDynamicDeckRuntime)
    state = runtime.timestep_state
    admittance = state.switch_admittance.admittance
    solution = state.switch_current.network_solution
    rhs = state.switch_current.rhs
    if size(admittance, 1) == length(solution) == length(rhs)
        residual = admittance * solution - rhs
        return maximum(abs, @view(residual[2:end]); init = 0.0)
    end
    context = runtime.context
    isempty(context.system.v) && return 0.0
    residual = context.system.y * context.system.v - context.system.rhs
    return maximum(abs, residual; init = 0.0)
end

function _append_restart_control_outputs(
    trace::DeckEMTTrace,
    full_trace::DeckEMTTrace,
    parsed::DeckParser.DeckParseResult,
    checkpoint_sample_count::Int,
)
    names = Symbol[]
    switch_spec = _deck_control_system_switch_current_spec(parsed)
    switch_spec === nothing || append!(names, switch_spec.output_names)
    append!(names, _deck_control_system_trace_output_names(parsed))
    isempty(names) && return trace
    indices = Int[]
    for name in names
        index = findlast(==(name), full_trace.output_channel_names)
        index === nothing && throw(ArgumentError(
            "restart control output $name is missing from the runtime trace",
        ))
        push!(indices, index)
    end
    values = copy(full_trace.output_pu[indices, :])
    if checkpoint_sample_count > 0
        original_outputs = _deck_control_system_trace_outputs(parsed, trace)
        original_outputs === nothing && throw(ArgumentError(
            "restart checkpoint control outputs could not be reconstructed",
        ))
        original_outputs.output_names == names || throw(ArgumentError(
            "restart checkpoint control output order changed",
        ))
        values[:, 1:checkpoint_sample_count] .=
            original_outputs.values[:, 1:checkpoint_sample_count]
    end
    extrema = _sampled_trace_extrema(values, trace.time_s)
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        vcat(copy(trace.output_channel_names), names),
        copy(trace.output_node_indices),
        vcat(copy(trace.output_pu), values),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        vcat(copy(trace.output_maximum_values), extrema.maximum_values),
        vcat(copy(trace.output_maximum_times_s), extrema.maximum_times_s),
        vcat(copy(trace.output_minimum_values), extrema.minimum_values),
        vcat(copy(trace.output_minimum_times_s), extrema.minimum_times_s),
    )
end

"""
Continue a completed prepared EMT workspace from its typed Julia state.

Mutation cards are applied in deck order before the first appended timestep.
The original trace, extrema, electromagnetic histories, sparse switch state,
source state, and control-device histories remain owned by the workspace.
"""
function restart_emt_study!(
    workspace::EMTStudyWorkspace,
    request::DeckParser.DeckRestartRequest;
    additional_time_s::Real,
    recorded_step_indices = nothing,
)
    workspace.ready && throw(ArgumentError(
        "restart requires an evaluated EMT workspace checkpoint",
    ))
    workspace.execution_mode === :monolithic || throw(ArgumentError(
        "monolithic restart requires a monolithic EMT workspace checkpoint",
    ))
    context = workspace.runtime.context
    checkpoint_step = context.step_count
    checkpoint_time_s = checkpoint_step * context.dt_s
    checkpoint_column = findfirst(==(checkpoint_step), context.recorded_step_indices)
    checkpoint_column === nothing && throw(ArgumentError(
        "restart checkpoint timestep is absent from the stored trace",
    ))
    checkpoint_state_error = maximum(
        abs,
        context.system.v .- context.voltage_pu[:, checkpoint_column];
        init = 0.0,
    )
    checkpoint_sample_count = length(context.recorded_step_indices)

    mutation_records = EMTRestartMutationRecord[]
    ordered_mutations = Any[
        request.switch_mutations...,
        request.control_source_mutations...,
    ]
    sort!(ordered_mutations; by = mutation -> mutation.line_no)
    for mutation in ordered_mutations
        record = if mutation isa DeckParser.DeckRestartSwitchMutation
            _apply_restart_switch_mutation!(workspace.runtime, mutation)
        else
            _apply_restart_control_source_mutation!(
                workspace.runtime,
                workspace.parsed,
                mutation,
            )
        end
        push!(mutation_records, record)
    end
    appended_step_count = _extend_emt_restart_horizon!(
        context,
        Float64(additional_time_s);
        recorded_step_indices = recorded_step_indices,
    )
    boundary_run = _run_prepared_dynamic_deck!(workspace.runtime)
    requested_trace = _deck_requested_electrical_trace(
        workspace.parsed,
        boundary_run.trace,
    )
    trace = _append_restart_control_outputs(
        requested_trace,
        boundary_run.trace,
        workspace.parsed,
        checkpoint_sample_count,
    )
    workspace.evaluation_count += 1
    return EMTRestartRun(
        trace,
        request,
        checkpoint_time_s,
        context.t_end_s,
        appended_step_count,
        mutation_records,
        checkpoint_state_error,
        _restart_final_kcl_error(workspace.runtime),
    )
end

function write_emt_restart_report(path::AbstractString, run::EMTRestartRun)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"schema\": \"aimora.emt_restart_run.v1\",")
        @printf(io, "  \"checkpoint_time_s\": %.17g,\n", run.checkpoint_time_s)
        @printf(io, "  \"final_time_s\": %.17g,\n", run.final_time_s)
        @printf(io, "  \"appended_step_count\": %d,\n", run.appended_step_count)
        @printf(io, "  \"mutation_count\": %d,\n", length(run.mutation_records))
        @printf(io, "  \"checkpoint_state_error\": %.17g,\n", run.checkpoint_state_error)
        @printf(io, "  \"final_kcl_error\": %.17g,\n", run.final_kcl_error)
        println(io, "  \"mutations\": [")
        for (index, record) in enumerate(run.mutation_records)
            suffix = index == length(run.mutation_records) ? "" : ","
            println(io, "    {")
            @printf(io, "      \"kind\": \"%s\",\n", String(record.kind))
            @printf(io, "      \"target\": \"%s\",\n", String(record.target))
            @printf(io, "      \"line_no\": %d,\n", record.line_no)
            println(io, "      \"parameters\": [", join(("\"$(name)\"" for name in record.parameter_names), ", "), "],")
            previous_values = join(
                (isfinite(value) ? string(value) : "null" for value in record.previous_values),
                ", ",
            )
            applied_values = join(
                (isfinite(value) ? string(value) : "null" for value in record.applied_values),
                ", ",
            )
            println(io, "      \"previous_values\": [$previous_values],")
            println(io, "      \"applied_values\": [$applied_values]")
            println(io, "    }$suffix")
        end
        println(io, "  ]")
        println(io, "}")
    end
    return path
end

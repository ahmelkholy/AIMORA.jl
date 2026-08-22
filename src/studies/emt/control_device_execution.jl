
function _control_system_device_weighted_input(
    state::OVER16CSUPState,
    row::OVER16CSUPDeviceRow,
)
    row.device_type in (60, 61, 63, 67) && return 0.0
    return sum(
        term.xtcs_index == 0 ? 0.0 : state.xtcs_values[term.xtcs_index] * term.scale
        for term in row.input_terms;
        init = 0.0,
    )
end

function _control_system_device_term_value(
    state::OVER16CSUPState,
    term::OVER16CSUPDeviceInputTerm,
)
    term.xtcs_index == 0 && return 0.0
    return state.xtcs_values[term.xtcs_index] * term.scale
end

function _control_system_digitizer_value(
    state::OVER16CSUPState,
    row::OVER16CSUPDeviceRow,
    weighted_input::Float64,
)
    gain = state.parsup_values[row.parsup_index]
    value = gain == 0.0 ? weighted_input : gain * weighted_input
    output = state.parsup_values[row.table_start_index]
    for index in row.table_end_index:-1:row.table_start_index
        if value >= state.parsup_values[index]
            output = state.parsup_values[index]
            break
        end
    end
    return output
end

function _control_system_point_value(
    state::OVER16CSUPState,
    row::OVER16CSUPDeviceRow,
    weighted_input::Float64,
)
    gain = state.parsup_values[row.parsup_index]
    value = gain == 0.0 ? weighted_input : gain * weighted_input
    first_x = state.parsup_values[row.table_start_index]
    first_y = state.parsup_values[row.table_start_index + 1]
    value <= first_x && return first_y
    for x_index in (row.table_start_index + 2):2:row.table_end_index
        next_x = state.parsup_values[x_index]
        value <= next_x || continue
        prior_x = state.parsup_values[x_index - 2]
        prior_y = state.parsup_values[x_index - 1]
        next_y = state.parsup_values[x_index + 1]
        next_x != prior_x ||
            throw(ArgumentError("control-system point table abscissas must differ"))
        return prior_y + (next_y - prior_y) / (next_x - prior_x) * (value - prior_x)
    end
    return state.parsup_values[row.table_end_index + 1]
end

function _control_system_device_steady_state_value(
    state::OVER16CSUPState,
    row::OVER16CSUPDeviceRow,
)
    weighted_input = _control_system_device_weighted_input(state, row)
    p1 = state.parsup_values[row.parsup_index]
    p2 = row.device_type == 66 ? 0.0 : state.parsup_values[row.parsup_index + 1]
    p3 = row.device_type in (60, 66) ? 0.0 :
        state.parsup_values[row.parsup_index + 2]
    if row.device_type == 50
        return p1
    elseif row.device_type in (51, 52)
        return p3 <= -1.0 ? (p1 == 0.0 ? weighted_input : p1 * weighted_input) : 0.0
    elseif row.device_type == 53
        return 0.0
    elseif row.device_type == 55
        return p2 == 1.0 ? _control_system_digitizer_value(state, row, weighted_input) : 0.0
    elseif row.device_type == 56
        return p2 == 1.0 ? _control_system_point_value(state, row, weighted_input) : 0.0
    elseif row.device_type == 60
        selector = trunc(Int, p2)
        1 <= selector <= 3 || return 0.0
        return _control_system_device_term_value(state, row.input_terms[4 - selector])
    elseif row.device_type == 61
        selector = trunc(Int, p1)
        1 <= selector <= 8 || return 0.0
        if selector <= 5
            index = length(row.input_terms) - selector + 1
            return 1 <= index <= length(row.input_terms) ?
                _control_system_device_term_value(state, row.input_terms[index]) : 0.0
        elseif selector == 6
            return row.control_index == 0 ? 0.0 : state.xtcs_values[row.control_index]
        elseif selector == 7
            return p3
        end
        return p2
    elseif row.device_type == 62
        return p2 < 0.5 ? 0.0 : p2 > 1.5 ? weighted_input : p3
    elseif row.device_type == 64
        p1 == 0.0 && return 0.0
        if p1 == 99999.0
            state.parsup_values[row.parsup_index] = 0.0
            p1 = 0.0
        end
        return weighted_input + p1 * p2
    elseif row.device_type == 66
        sample_count = row.control_index
        sum_squares = sum(
            state.parsup_values[(row.parsup_index + 1):(row.parsup_index + sample_count)];
            init = 0.0,
        )
        return sqrt(sum_squares * p1 + weighted_input * weighted_input)
    elseif row.device_type == 67
        return 1.0
    end
    return 0.0
end

function _initialize_control_system_device_histories!(
    runtime::ControlSystemSupplementalDeviceRuntime,
)
    runtime.initialized && return runtime
    for (index, row) in enumerate(runtime.rows)
        output = _control_system_device_steady_state_value(runtime.state, row)
        runtime.state.xtcs_values[runtime.device_output_slots[index]] = output
    end
    for (index, row) in enumerate(runtime.rows)
        weighted_input = _control_system_device_weighted_input(runtime.state, row)
        if runtime.initialize_transport_delay_from_input[index]
            history_start = trunc(Int, runtime.state.parsup_values[row.parsup_index])
            history_count = trunc(Int, runtime.state.parsup_values[row.parsup_index + 2])
            runtime.state.parsup_values[history_start:(history_start + history_count - 1)] .=
                weighted_input
        elseif runtime.initialize_rms_from_input[index]
            sample_count = row.control_index
            runtime.state.parsup_values[(row.parsup_index + 1):(row.parsup_index + sample_count)] .=
                weighted_input * weighted_input
        end
        if row.device_type == 50
            runtime.state.device_integer_values[index] = weighted_input < 0.0 ? -1 : 1
            runtime.state.parsup_values[row.parsup_index + 1] = -1.0
            runtime.state.parsup_values[row.parsup_index + 2] = weighted_input
        elseif row.device_type == 58
            output = runtime.state.xtcs_values[runtime.device_output_slots[index]]
            denominator = runtime.state.parsup_values[row.parsup_index + 1]
            coefficient = runtime.state.parsup_values[row.parsup_index + 2]
            runtime.state.parsup_values[row.parsup_index] =
                (denominator - coefficient) / 2.0 * output
        elseif row.device_type == 59
            runtime.state.parsup_values[row.parsup_index + 1] = weighted_input
        elseif row.device_type == 62
            runtime.state.parsup_values[row.parsup_index + 2] =
                runtime.state.xtcs_values[runtime.device_output_slots[index]]
        elseif row.device_type in (64, 65)
            runtime.state.parsup_values[row.parsup_index] =
                runtime.state.xtcs_values[runtime.device_output_slots[index]]
        end
    end
    runtime.initialized = true
    return runtime
end

function _advance_control_system_supplemental_devices!(
    runtime::ControlSystemSupplementalDeviceRuntime,
    state::ControlSystemExecutionState,
    step::Int,
    time_s::Float64,
)
    for (index, name) in enumerate(runtime.ordinary_signal_names)
        runtime.state.xtcs_values[index] = get(state.values, name, 0.0)
    end
    _initialize_control_system_device_histories!(runtime)
    update = over16_csup_device_update!(
        runtime.state,
        runtime.rows;
        next_indices = runtime.next_indices,
        start_index = 1,
        t = time_s,
        deltat = state.deltat_s,
        onehalf = 0.5,
        nuk = length(runtime.ordinary_signal_names),
    )
    for (index, row) in enumerate(runtime.rows)
        row.device_type == 53 || continue
        pointer = runtime.transport_delay_pointers[index]
        runtime.state.parsup_values[pointer] = update.weighted_input_values[index]
        history_start = trunc(Int, runtime.state.parsup_values[row.parsup_index])
        history_count = trunc(Int, runtime.state.parsup_values[row.parsup_index + 2])
        pointer += 1
        runtime.transport_delay_pointers[index] =
            pointer == history_start + history_count ? history_start : pointer
        runtime.rows[index] = over16_csup_device_row(
            row.supplemental_index,
            row.device_type,
            row.parsup_index;
            input_terms = row.input_terms,
            control_index = row.control_index,
            reference_index = runtime.transport_delay_pointers[index],
        )
    end
    for index in eachindex(runtime.device_output_names)
        state.values[runtime.device_output_names[index]] =
            runtime.state.xtcs_values[runtime.device_output_slots[index]]
    end
    runtime.executed_step_count += 1
    return update
end

function _control_system_supplemental_step!(
    runtime::ControlSystemSupplementalDeviceRuntime,
    state::ControlSystemExecutionState,
    control_expressions::AbstractVector{ControlExpressionRuntime},
    step::Int,
    time_s::Real,
)
    time = Float64(time_s)
    advance_control_system_state!(state, step, time; execute_devices = false)
    _advance_control_system_expressions!(control_expressions, state)
    update = _advance_control_system_supplemental_devices!(runtime, state, step, time)
    assignment_outputs = Symbol[assignment.output_name for assignment in state.assignments]
    function_outputs = Symbol[function_row.output_name for function_row in state.functions]
    output_values = control_system_output_values(state)
    full_device_execution =
        !isempty(runtime.rows) &&
        update.device_executed &&
        length(update.processed_supplemental_indices) == length(runtime.rows)
    deferred_effects = Symbol[
        :full_use1_stack_replay,
        :full_xref1_interpreter,
        :full_errstp_equivalence,
        :full_tacs_solvum_coupling,
    ]
    return (
        step = step,
        time_s = time,
        deltat_s = state.deltat_s,
        frequency_hz = state.frequency_hz,
        assignment_outputs = assignment_outputs,
        assignment_values = Float64[state.values[name] for name in assignment_outputs],
        function_outputs = function_outputs,
        function_values = Float64[state.values[name] for name in function_outputs],
        device_outputs = copy(runtime.device_output_names),
        device_values = Float64[state.values[name] for name in runtime.device_output_names],
        unsupported_device_types = Int[],
        output_names = copy(state.output_names),
        output_values = output_values,
        assignment_count = length(state.assignments) + length(state.functions),
        function_count = length(state.functions),
        device_count = length(runtime.rows),
        executed_device_count = length(update.processed_supplemental_indices),
        output_count = length(output_values),
        control_system_executed = true,
        algebraic_assignments_executed =
            !isempty(state.assignments) || !isempty(state.functions),
        supported_assignment_complete = true,
        supported_function_complete = true,
        supported_device_complete = full_device_execution,
        complete_control_system_executed = full_device_execution,
        machine_coupling_executed = false,
        full_tacs3_device_execution = full_device_execution,
        full_tacs_solvum_coupling = false,
        deferred_effects = Tuple(deferred_effects),
    )
end

function _control_system_supplemental_step!(
    runtime::ControlSystemSupplementalDeviceRuntime,
    state::ControlSystemExecutionState,
    step::Int,
    time_s::Real,
)
    return _control_system_supplemental_step!(
        runtime,
        state,
        ControlExpressionRuntime[],
        step,
        time_s,
    )
end

function _source_function_network_runtime(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    source_signal_provider::AbstractSourceSignalProvider,
)
    plan = deck_over16_boundary_plan(parsed)
    source_rows = DeckParser.deck_over5a_source_rows(parsed)
    source_slots = Base.RefValue{Float64}[Ref(0.0) for _ in 1:10]
    row_slots = Base.RefValue{Float64}[Ref(0.0) for _ in source_rows]
    dynamic_row_indices = Int[]
    source_element_count = 0
    previous_iform = 0
    for (row_index, row) in enumerate(source_rows)
        source_type = abs(row.iform)
        type16_successor = previous_iform == 16
        type17_successor = previous_iform == 17
        previous_iform = row.iform
        if 1 <= source_type <= 10 && !type16_successor && !type17_successor
            node_index = abs(row.node_value)
            node_index > 0 ||
                throw(ArgumentError("source-function rows require a nonzero node"))
            slot = source_slots[source_type]
            element = row.node_value < 0 ?
                CurrentInjection(node_index, _time_s -> slot[]) :
                TheveninSource(
                    node_index,
                    DeckParser.BPA_FIXED_SOURCE_CONDUCTANCE,
                    _time_s -> slot[],
                )
            push!(elements, element)
            push!(element_names, Symbol("source_function_", String(row.name)))
            source_element_count += 1
            continue
        end
        dynamic_source = source_type >= 60 || type17_successor
        dynamic_source || continue
        node_index = abs(row.node_value)
        node_index > 0 || throw(ArgumentError("source-function rows require a nonzero node"))
        if source_type in 11:15
            element_index = findfirst(==(row.name), element_names)
            element_index === nothing || begin
                deleteat!(elements, element_index)
                deleteat!(element_names, element_index)
            end
        end
        slot = row_slots[row_index]
        element = row.node_value < 0 ?
            CurrentInjection(node_index, _time_s -> slot[]) :
            TheveninSource(
                node_index,
                DeckParser.BPA_FIXED_SOURCE_CONDUCTANCE,
                _time_s -> slot[],
        )
        push!(elements, element)
        push!(element_names, Symbol("source_row_", String(row.name)))
        push!(dynamic_row_indices, row_index)
        source_element_count += 1
    end
    source_element_count > 0 || return nothing
    validate_source_signal_program_units(
        source_signal_provider,
        plan.source_node_values,
        plan.source_iform_values,
    )

    internal_analytic_requested = any(
        row -> row.request_kind == :analytic_source_usage,
        DeckParser.deck_study_option_request_rows(parsed),
    )
    return SourceFunctionNetworkRuntime(
        OVER16SourceCardState(
            zeros(10);
            iread = plan.source_card_row_count > 0 ? 1 : 0,
        ),
        plan,
        source_slots,
        row_slots,
        dynamic_row_indices,
        source_signal_provider,
        internal_analytic_requested,
        0,
        0,
        0,
        0,
        0,
        0,
        SourceSignalStageSample[],
        -Inf,
        1,
    )
end

function _synchronize_source_function_row_slots!(
    runtime::SourceFunctionNetworkRuntime,
    time_s::Float64,
    xtcs_values::AbstractVector{<:Real},
)
    isempty(runtime.dynamic_row_indices) && return nothing
    plan = runtime.plan
    node_count = maximum(abs, plan.source_node_values; init = 0)
    row_result = over16_source_row_update!(
        zeros(Float64, node_count),
        zeros(Float64, node_count),
        copy(plan.source_node_values),
        copy(plan.source_iform_values),
        copy(plan.source_tstart_values),
        copy(plan.source_tstop_values),
        runtime.state.voltbc_values,
        time_s;
        kconst = plan.source_row_count,
        crest_values = copy(plan.source_crest_values),
        time1_values = copy(plan.source_time1_values),
        time2_values = copy(plan.source_time2_values),
        sfreq_values = copy(plan.source_sfreq_values),
        xtcs_values = xtcs_values,
    )
    unresolved_dynamic_row_indices =
        intersect(runtime.dynamic_row_indices, row_result.deferred_row_indices)
    isempty(unresolved_dynamic_row_indices) || throw(ArgumentError(
        "dynamic source rows remain unresolved: " *
        join(unresolved_dynamic_row_indices, ","),
    ))
    for row_index in runtime.dynamic_row_indices
        runtime.row_slot_values[row_index][] = row_result.source_values[row_index]
    end
    return row_result
end

function _control_system_signal_slot_values(
    context::EMTStepContext;
    maximum_index::Int=0,
)
    control_runtime = context.control_system_runtime
    control_runtime === nothing &&
        throw(ArgumentError("control-system feedback requires an executable control-system runtime"))
    control_runtime.halted && throw(ArgumentError(
        "control-system feedback is halted: $(control_runtime.last_diagnostic)",
    ))
    maximum_index <= length(control_runtime.signal_slot_names) || throw(ArgumentError(
        "control-system feedback index $maximum_index exceeds " *
        "$(length(control_runtime.signal_slot_names)) ordered signal slots",
    ))
    values = Vector{Float64}(undef, length(control_runtime.signal_slot_names))
    for index in eachindex(values)
        name = control_runtime.signal_slot_names[index]
        haskey(control_runtime.state.values, name) || throw(ArgumentError(
            "ordered control-system signal slot $index references unresolved signal $name",
        ))
        value = control_runtime.state.values[name]
        isfinite(value) || throw(ArgumentError(
            "ordered control-system signal slot $index ($name) is nonfinite",
        ))
        values[index] = value
    end
    control_runtime.feedback_application_count += 1
    return values
end

function _source_function_control_signal_values(context::EMTStepContext)
    runtime = context.source_function_runtime
    runtime === nothing && return Float64[]
    plan = runtime.plan
    row_indices = Int[
        round(Int, plan.source_sfreq_values[index])
        for index in eachindex(plan.source_iform_values)
        if (
            plan.source_iform_values[index] == 17 ||
            abs(plan.source_iform_values[index]) >= 60
        ) &&
           isinteger(plan.source_sfreq_values[index]) &&
           plan.source_sfreq_values[index] > 0.0
    ]
    override_indices =
        plan.source_tacs_override_count == 0 ?
        Int[] :
        copy(plan.source_tacs_override_xtcs_indices)
    referenced_indices = vcat(row_indices, override_indices)
    isempty(referenced_indices) && return Float64[]
    return _control_system_signal_slot_values(
        context;
        maximum_index = maximum(referenced_indices),
    )
end

function _advance_control_system_owner!(
    state::ControlSystemExecutionState,
    supplemental::Union{Nothing,ControlSystemSupplementalDeviceRuntime},
    control_expressions::AbstractVector{ControlExpressionRuntime},
    step::Int,
    time_s::Real,
)
    if supplemental === nothing
        isempty(state.devices) || throw(ArgumentError(
            "parsed control-system devices require the complete supplemental-device owner",
        ))
        advance_control_system_state!(state, step, time_s; execute_devices = false)
        _advance_control_system_expressions!(control_expressions, state)
        return nothing
    end
    return _control_system_supplemental_step!(
        supplemental,
        state,
        control_expressions,
        step,
        time_s,
    )
end

function _set_controlled_switch_inputs!(runtime::ControlSystemNetworkRuntime)
    for index in eachindex(runtime.switch_elements)
        runtime.switch_elements[index].control[] =
            get(runtime.state.values, runtime.switch_signal_names[index], 0.0)
        clamp_name = runtime.switch_clamp_signal_names[index]
        clamp_name === missing && continue
        clamp_control = runtime.switch_elements[index].clamp_control
        clamp_control === nothing && throw(ArgumentError(
            "controlled valve clamp signal $clamp_name has no physical switch input",
        ))
        clamp_control[] = get(runtime.state.values, clamp_name, 0.0)
    end
    return runtime
end

function _sample_control_switch_observations!(
    runtime::ControlSystemNetworkRuntime,
    network_voltage::AbstractVector{<:Real},
    time_s::Float64,
    ordinary_switch_currents::Union{Nothing,AbstractVector{<:Real}},
    ordinary_switch_closed_flags::Union{Nothing,AbstractVector{<:Integer}},
)
    for observation in runtime.switch_observations
        switch = observation.switch
        ordinary_switch =
            switch isa Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch}
        value = 0.0
        if ordinary_switch
            if observation.source_type == 91
                if ordinary_switch_currents !== nothing &&
                   observation.ordinary_switch_index > 0
                    value = Float64(
                        ordinary_switch_currents[observation.ordinary_switch_index],
                    )
                else
                    from_voltage = switch.a == 0 ? 0.0 :
                        Float64(network_voltage[switch.a])
                    to_voltage = switch.b == 0 ? 0.0 :
                        Float64(network_voltage[switch.b])
                    value =
                        switch_conductance(switch, time_s) *
                        (from_voltage - to_voltage)
                end
            elseif observation.source_type == 93
                if ordinary_switch_closed_flags !== nothing &&
                   observation.ordinary_switch_index > 0
                    value =
                        ordinary_switch_closed_flags[
                            observation.ordinary_switch_index
                        ] == 0 ? 0.0 : 1.0
                else
                    value = switch_closed(switch, time_s) ? 1.0 : 0.0
                end
            else
                throw(ArgumentError(
                    "unsupported control-system switch observation type " *
                    "$(observation.source_type)",
                ))
            end
        elseif observation.source_type == 91
            value = switch.last_current
        elseif observation.source_type == 93
            value = controlled_switch_closed(switch) ? 1.0 : 0.0
        else
            throw(ArgumentError(
                "unsupported control-system switch observation type " *
                "$(observation.source_type)",
            ))
        end
        runtime.state.values[observation.name] = value
    end
    return runtime
end

function _advance_control_system_network_runtime!(
    runtime::ControlSystemNetworkRuntime,
    step::Int,
    time_s::Real;
    network_voltage::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ordinary_switch_currents::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ordinary_switch_closed_flags::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)
    runtime.halted && throw(ArgumentError(
        "control-system execution is halted: $(runtime.last_diagnostic)",
    ))
    time = Float64(time_s)
    try
        for source in runtime.windowed_constant_sources
            runtime.state.values[source.name] =
                _windowed_constant_control_signal_value(source, time)
        end
        for source in runtime.sinusoidal_sources
            runtime.state.values[source.name] =
                sinusoidal_control_signal_value(source, time)
        end
        for source in runtime.waveform_sources
            runtime.state.values[source.name] =
                _control_waveform_signal_value(source, time)
        end
        if network_voltage !== nothing
            for index in eachindex(runtime.network_voltage_source_names)
                runtime.state.values[runtime.network_voltage_source_names[index]] =
                    Float64(network_voltage[runtime.network_voltage_source_node_indices[index]])
            end
            _project_closed_controlled_switch_voltages!(runtime, network_voltage)
            _sample_control_switch_observations!(
                runtime,
                network_voltage,
                time,
                ordinary_switch_currents,
                ordinary_switch_closed_flags,
            )
        end
        result = _advance_control_system_owner!(
            runtime.state,
            runtime.supplemental_devices,
            runtime.control_expressions,
            step,
            time,
        )
        for (index, name) in enumerate(runtime.signal_slot_names)
            haskey(runtime.state.values, name) || throw(ArgumentError(
                "ordered control-system signal slot $index references unresolved signal $name",
            ))
            isfinite(runtime.state.values[name]) || throw(ArgumentError(
                "ordered control-system signal slot $index ($name) became nonfinite",
            ))
        end
        _set_controlled_switch_inputs!(runtime)
        runtime.executed_step_count += 1
        return result
    catch error
        runtime.halted = true
        runtime.last_diagnostic =
            "step $step at $(time) s: $(sprint(showerror, error))"
        rethrow()
    end
end

function _initialize_control_system_network_steady_state!(
    runtime::ControlSystemNetworkRuntime,
    network_voltage::AbstractVector{<:Real},
)
    for index in eachindex(runtime.network_voltage_source_names)
        runtime.state.values[runtime.network_voltage_source_names[index]] =
            Float64(network_voltage[runtime.network_voltage_source_node_indices[index]])
    end
    _project_closed_controlled_switch_voltages!(runtime, network_voltage)
    initialize_control_function_steady_state!(runtime.state)
    supplemental = runtime.supplemental_devices
    if supplemental !== nothing
        for (index, name) in enumerate(supplemental.ordinary_signal_names)
            supplemental.state.xtcs_values[index] =
                get(runtime.state.values, name, 0.0)
        end
        _initialize_control_system_device_histories!(supplemental)
        for index in eachindex(supplemental.device_output_names)
            runtime.state.values[supplemental.device_output_names[index]] =
                supplemental.state.xtcs_values[
                    supplemental.device_output_slots[index]
                ]
        end
    end
    _apply_control_system_frequency_initializations!(runtime)
    _set_controlled_switch_inputs!(runtime)
    for switch in runtime.switch_elements
        sync_controlled_switch!(switch, 0.0)
    end
    return runtime
end

function _initialize_control_system_network_steady_state!(
    context::EMTStepContext,
    network_voltage::AbstractVector{<:Real},
)
    runtime = context.control_system_runtime
    runtime === nothing && return nothing
    return _initialize_control_system_network_steady_state!(
        runtime,
        network_voltage,
    )
end

function _source_signal_stage_input_values(
    state::OVER16SourceCardState,
    card_kwargs::NamedTuple,
)
    if haskey(card_kwargs, :interpolated_values)
        return Float64.(card_kwargs.interpolated_values)
    elseif state.iread != 0 && haskey(card_kwargs, :card_values)
        values = Float64.(card_kwargs.card_values)
        return !isempty(values) && values[1] == 9999.0 ?
            zeros(Float64, 10) :
            values
    elseif state.iread == 0
        return zeros(Float64, 10)
    end
    return copy(state.voltbc_values)
end

function _source_signal_post_tacs_values(
    input_values::AbstractVector{<:Real},
    tacs_kwargs::NamedTuple,
)
    values = Float64.(input_values)
    nstacs = Int(get(tacs_kwargs, :nstacs, 0))
    indices = get(tacs_kwargs, :vstacs_indices, Int[])
    control_values = get(tacs_kwargs, :xtcs_values, Float64[])
    for position in 1:nstacs
        source_index = indices[position]
        source_index == 0 && continue
        values[position] = Float64(control_values[source_index])
    end
    return values
end

function _source_signal_program_analytic_kwargs(
    provider::AbstractSourceSignalProvider,
    state::OVER16SourceCardState,
    card_kwargs::NamedTuple,
    tacs_kwargs::NamedTuple,
    time_s::Real,
)
    source_signal_analytic_active(provider) ||
        return (kwargs = NamedTuple(), assignment_indices = Int[])
    post_tacs = _source_signal_post_tacs_values(
        _source_signal_stage_input_values(state, card_kwargs),
        tacs_kwargs,
    )
    analytic = source_signal_analytic_values(provider, post_tacs, time_s)
    return (
        kwargs = (kanal = 1, analytic_values = analytic.values),
        assignment_indices = copy(analytic.assignment_indices),
    )
end

function _record_source_function_update!(
    runtime::SourceFunctionNetworkRuntime,
    update,
    step_index::Int,
    time_s::Real;
    interpolation_applied::Bool=false,
    analytic_assignment_indices::AbstractVector{Int}=Int[],
    accepted_values::AbstractVector{<:Real}=update.output_voltbc_values,
    record_stage_sample::Bool=true,
)
    time = Float64(time_s)
    time >= runtime.last_accepted_time_s || throw(ArgumentError(
        "source-signal accepted time moved backward from " *
        "$(runtime.last_accepted_time_s) s to $time s",
    ))
    for index in eachindex(runtime.slot_values)
        runtime.slot_values[index][] = runtime.state.voltbc_values[index]
    end
    runtime.executed_step_count += 1
    runtime.card_read_count += Int(update.card_decoded)
    runtime.signal_synchronization_count += Int(update.source_signal_synchronized)
    runtime.external_signal_count += Int(interpolation_applied)
    runtime.tacs_override_count += update.tacs_override_count
    runtime.analytic_execution_count += Int(update.analyt_executed)
    runtime.last_accepted_time_s = time
    runtime.next_input_row_index += 1
    if record_stage_sample
        push!(
            runtime.stage_samples,
            SourceSignalStageSample(
                step_index,
                time,
                copy(update.interp_input_values),
                copy(update.pre_tacs_voltbc_values),
                copy(update.analytic_input_voltbc_values),
                copy(update.analytic_output_voltbc_values),
                Float64.(accepted_values),
                update.card_decoded,
                interpolation_applied,
                update.tacs_override_count,
                Int.(analytic_assignment_indices),
            ),
        )
    end
    return update
end

function _source_signal_stage_recording_active(
    runtime::SourceFunctionNetworkRuntime,
)
    return source_signal_interpolation_active(runtime.signal_provider) ||
           source_signal_analytic_active(runtime.signal_provider)
end

function _source_signal_card_preview_kwargs(source_kwargs::NamedTuple)
    vstacs_indices = get(source_kwargs, :vstacs_indices, Int[])
    return (
        nchain = get(source_kwargs, :nchain, 16),
        kolbeg = get(source_kwargs, :kolbeg, 0),
        card_values = get(source_kwargs, :card_values, nothing),
        interpolated_values = get(
            source_kwargs,
            :interpolated_values,
            nothing,
        ),
        source_signal_provider_applied = get(
            source_kwargs,
            :source_signal_provider_applied,
            nothing,
        ),
        xtcs_values = get(source_kwargs, :xtcs_values, Float64[]),
        vstacs_indices = vstacs_indices,
        nstacs = get(source_kwargs, :nstacs, length(vstacs_indices)),
        kxtcs = get(source_kwargs, :kxtcs, 0),
        kanal = get(source_kwargs, :kanal, 0),
        analytic_values = get(source_kwargs, :analytic_values, nothing),
        analytic_source_count = get(
            source_kwargs,
            :analytic_source_count,
            0,
        ),
        analytic_branch_count = get(
            source_kwargs,
            :analytic_branch_count,
            0,
        ),
        analytic_time_s = get(source_kwargs, :analytic_time_s, 0.0),
        analytic_step_s = get(source_kwargs, :analytic_step_s, 0.0),
    )
end

function _accept_source_function_boundary_update!(
    context::EMTStepContext,
    accepted_state::OVER16AcceptedTimestepState,
    source_config,
    over16_update,
    ;
    idempotent_same_time::Bool = false,
)
    runtime = context.source_function_runtime
    (runtime === nothing || source_config === nothing) && return nothing
    source_kwargs = get(source_config, :kwargs, NamedTuple())
    source_result =
        hasproperty(over16_update, :source_result) ?
        over16_update.source_result :
        nothing
    record_stage_sample = _source_signal_stage_recording_active(runtime)
    source_time = Float64(get(source_config, :t, context.t_s))
    program_assignment_indices =
        source_signal_analytic_active(runtime.signal_provider) ?
        _source_signal_program_analytic_kwargs(
            runtime.signal_provider,
            runtime.state,
            source_kwargs,
            source_kwargs,
            source_time,
        ).assignment_indices :
        Int[]
    lean_preview =
        record_stage_sample &&
        (
            source_result === nothing ||
            !hasproperty(source_result, :source_card_result)
        ) ?
        over16_source_card_update_intent_preview(
            runtime.state.voltbc_values;
            iread = runtime.state.iread,
            _source_signal_card_preview_kwargs(source_kwargs)...,
        ) :
        nothing
    copyto!(
        runtime.state.voltbc_values,
        accepted_state.source.source_card.voltbc_values,
    )
    runtime.state.iread = accepted_state.source.source_card.iread
    runtime.state.nfrfld = accepted_state.source.source_card.nfrfld
    for index in eachindex(runtime.slot_values)
        runtime.slot_values[index][] = runtime.state.voltbc_values[index]
    end
    if idempotent_same_time &&
       abs(source_time - runtime.last_accepted_time_s) <=
       16.0 * eps(Float64) * max(1.0, abs(source_time))
        return source_result
    end
    if source_result === nothing ||
       !hasproperty(source_result, :source_card_result)
        if lean_preview !== nothing
            return _record_source_function_update!(
                runtime,
                lean_preview,
                context.step_index,
                source_time;
                interpolation_applied = get(
                    source_kwargs,
                    :source_signal_provider_applied,
                    false,
                ),
                analytic_assignment_indices =
                    source_signal_analytic_active(runtime.signal_provider) ?
                    program_assignment_indices :
                    lean_preview.analytic_assignment_indices,
                accepted_values =
                    accepted_state.source.source_card.voltbc_values,
            )
        end
        source_time >= runtime.last_accepted_time_s || throw(ArgumentError(
            "source-signal accepted time moved backward from " *
            "$(runtime.last_accepted_time_s) s to $source_time s",
        ))
        runtime.executed_step_count += 1
        runtime.signal_synchronization_count += 1
        runtime.external_signal_count += Int(
            get(source_kwargs, :source_signal_provider_applied, false),
        )
        runtime.tacs_override_count += count(
            index -> index != 0,
            get(source_kwargs, :vstacs_indices, Int[]),
        )
        runtime.analytic_execution_count += Int(
            get(source_kwargs, :kanal, 0) > 0,
        )
        runtime.last_accepted_time_s = source_time
        runtime.next_input_row_index += 1
        return source_result
    end
    card_result = source_result.source_card_result
    return _record_source_function_update!(
        runtime,
        card_result,
        context.step_index,
        source_result.source_update_time_s;
        interpolation_applied = card_result.source_signal_provider_applied,
        analytic_assignment_indices =
            source_signal_analytic_active(runtime.signal_provider) ?
            program_assignment_indices :
            card_result.analytic_assignment_indices,
        accepted_values = accepted_state.source.source_card.voltbc_values,
        record_stage_sample = record_stage_sample,
    )
end

function _initialize_source_function_state!(
    context::EMTStepContext,
    accepted_state::OVER16AcceptedTimestepState,
    plan::DeckOVER16BoundaryPlan,
)
    runtime = context.source_function_runtime
    runtime === nothing && return nothing
    _source_signal_stage_recording_active(runtime) || return nothing
    provider = runtime.signal_provider
    card_kwargs = _deck_over16_source_card_step_kwargs(plan, context)
    interpolation_applied =
        get(card_kwargs, :source_signal_provider_applied, false)
    if !interpolation_applied &&
       source_signal_interpolation_active(provider) &&
       !(
           accepted_state.source.source_card.iread != 0 &&
           get(card_kwargs, :nchain, 16) == 16
       )
        input_values = accepted_state.source.source_card.iread == 0 ?
            zeros(Float64, 10) :
            Float64.(get(
                card_kwargs,
                :card_values,
                accepted_state.source.source_card.voltbc_values,
            ))
        if !isempty(input_values) && input_values[1] == 9999.0
            fill!(input_values, 0.0)
        end
        card_kwargs = merge(
            card_kwargs,
            (
                interpolated_values =
                    source_signal_values(provider, input_values, 0.0),
                source_signal_provider_applied = true,
            ),
        )
        interpolation_applied = true
    end
    control_signal_values = _source_function_control_signal_values(context)
    tacs_kwargs =
        plan.source_tacs_override_count == 0 ?
        (xtcs_values = control_signal_values,) :
        merge(
            _deck_over16_source_tacs_override_kwargs(plan),
            (xtcs_values = control_signal_values,),
        )
    analytic_kwargs = _deck_over16_source_analytic_kwargs(plan, context)
    analytic_assignment_indices = isempty(analytic_kwargs) ? Int[] : collect(1:10)
    if isempty(analytic_kwargs) && source_signal_analytic_active(provider)
        analytic_program = _source_signal_program_analytic_kwargs(
            provider,
            accepted_state.source.source_card,
            card_kwargs,
            tacs_kwargs,
            0.0,
        )
        analytic_kwargs = analytic_program.kwargs
        analytic_assignment_indices = analytic_program.assignment_indices
    elseif isempty(analytic_kwargs) && runtime.internal_analytic_requested
        throw(ArgumentError(
            "analytic source usage requires a typed SourceSignalProgram",
        ))
    end
    update = over16_source_update!(
        accepted_state.source,
        plan.source_node_values,
        plan.source_iform_values,
        plan.source_tstart_values,
        plan.source_tstop_values,
        0.0;
        kconst = plan.source_row_count,
        crest_values = plan.source_crest_values,
        time1_values = plan.source_time1_values,
        time2_values = plan.source_time2_values,
        sfreq_values = plan.source_sfreq_values,
        delta2 = context.dt_s / 2.0,
        merge(card_kwargs, tacs_kwargs, analytic_kwargs)...,
    )
    copyto!(
        runtime.state.voltbc_values,
        accepted_state.source.source_card.voltbc_values,
    )
    runtime.state.iread = accepted_state.source.source_card.iread
    runtime.state.nfrfld = accepted_state.source.source_card.nfrfld
    return _record_source_function_update!(
        runtime,
        update.source_card_result,
        0,
        0.0;
        interpolation_applied = interpolation_applied,
        analytic_assignment_indices = analytic_assignment_indices,
        record_stage_sample = _source_signal_stage_recording_active(runtime),
    )
end

function _advance_source_function_network!(context::EMTStepContext)
    runtime = context.source_function_runtime
    runtime === nothing && return nothing
    plan = runtime.plan
    card_kwargs = _deck_over16_source_card_step_kwargs(plan, context)
    provider_applied = false
    if !get(card_kwargs, :source_signal_provider_applied, false) &&
       source_signal_interpolation_active(runtime.signal_provider) &&
       !(runtime.state.iread != 0 && get(card_kwargs, :nchain, 16) == 16)
        input_values = runtime.state.iread == 0 ?
            zeros(Float64, 10) :
            Float64.(get(card_kwargs, :card_values, runtime.state.voltbc_values))
        if !isempty(input_values) && input_values[1] == 9999.0
            fill!(input_values, 0.0)
        end
        card_kwargs = merge(
            card_kwargs,
            (
                interpolated_values = source_signal_values(
                    runtime.signal_provider,
                    input_values,
                    context.t_s,
                ),
                source_signal_provider_applied = true,
            ),
        )
        provider_applied = true
    end
    control_signal_values = _source_function_control_signal_values(context)
    tacs_kwargs =
        plan.source_tacs_override_count == 0 ?
        (xtcs_values = control_signal_values,) :
        merge(
            _deck_over16_source_tacs_override_kwargs(plan),
            (xtcs_values = control_signal_values,),
        )
    analytic_kwargs = _deck_over16_source_analytic_kwargs(plan, context)
    analytic_assignment_indices = Int[]
    if isempty(analytic_kwargs) &&
       source_signal_analytic_active(runtime.signal_provider)
        analytic_program = _source_signal_program_analytic_kwargs(
            runtime.signal_provider,
            runtime.state,
            card_kwargs,
            tacs_kwargs,
            context.t_s,
        )
        analytic_kwargs = analytic_program.kwargs
        analytic_assignment_indices = analytic_program.assignment_indices
    elseif isempty(analytic_kwargs) && runtime.internal_analytic_requested
        throw(ArgumentError(
            "analytic source usage requires a typed SourceSignalProgram",
        ))
    elseif !isempty(analytic_kwargs)
        analytic_assignment_indices = collect(1:10)
    end
    update = over16_source_card_update!(
        runtime.state;
        merge(card_kwargs, tacs_kwargs, analytic_kwargs)...,
    )
    recorded = _record_source_function_update!(
        runtime,
        update,
        context.step_index,
        context.t_s;
        interpolation_applied =
            provider_applied ||
            get(card_kwargs, :source_signal_provider_applied, false),
        analytic_assignment_indices = analytic_assignment_indices,
        record_stage_sample = _source_signal_stage_recording_active(runtime),
    )
    _synchronize_source_function_row_slots!(
        runtime,
        context.t_s,
        get(tacs_kwargs, :xtcs_values, Float64[]),
    )
    return recorded
end

function _append_dynamic_source_and_control_elements!(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    source_runtime = _source_function_network_runtime(
        elements,
        element_names,
        parsed,
        source_signal_provider,
    )
    control_runtime = _append_control_system_network_elements!(
        elements,
        element_names,
        parsed,
        dt_s,
    )
    return source_runtime, control_runtime
end

function _deck_control_system_switch_current_spec(parsed::DeckParser.DeckParseResult)
    rows = [
        row for row in DeckParser.deck_control_system_switch_coupling_rows(parsed)
        if _switch_current_report_selected(row.output_code) &&
           row.from_node_index !== missing &&
           row.to_node_index !== missing
    ]
    isempty(rows) && return nothing
    return (
        rows = rows,
        output_names = Symbol[
            Symbol(String(row.from_node), ".", String(row.to_node))
            for row in rows
        ],
    )
end

function _control_system_switch_element_names(parsed::DeckParser.DeckParseResult)
    return Symbol[
        _control_system_switch_element_name(row, index)
        for (index, row) in
            enumerate(DeckParser.deck_control_system_switch_coupling_rows(parsed))
    ]
end

function _deck_control_system_signal_slot_names(
    parsed::DeckParser.DeckParseResult,
)
    names = Symbol[
        row.name for row in DeckParser.deck_control_system_function_rows(parsed)
    ]
    append!(names, _CONTROL_SYSTEM_UTILITY_VALUES)
    append!(
        names,
        (row.name for row in DeckParser.deck_control_system_source_rows(parsed)),
    )
    supplemental_rows = Any[
        DeckParser.deck_control_system_expression_rows(parsed)...,
        DeckParser.deck_control_system_device_rows(parsed)...,
    ]
    sort!(supplemental_rows; by = row -> row.line_no)
    append!(names, (row.name for row in supplemental_rows))
    duplicate_names = sort!(
        unique(name for name in names if count(==(name), names) > 1);
        by = String,
    )
    isempty(duplicate_names) || throw(ArgumentError(
        "control-system ordered signal slots contain duplicate names: " *
        join(duplicate_names, ", "),
    ))
    return names
end

function _deck_control_switch_observations(
    parsed::DeckParser.DeckParseResult,
    network_elements::AbstractVector,
    controlled_switches::Vector{TACSControlledSwitch},
)
    candidates = ControlObservedSwitch[
        element
        for element in Iterators.flatten((network_elements, controlled_switches))
        if element isa ControlObservedSwitch
    ]
    observations = ControlSwitchObservation[]
    ordinary_candidates = ControlObservedSwitch[
        element
        for element in network_elements
        if element isa Union{TimeSwitch,CurrentZeroSwitch}
    ]
    for row in DeckParser.deck_control_system_source_rows(parsed)
        row.source_type in (91, 93) || continue
        node_index = get(parsed.node_map, row.name, 0)
        node_index > 0 || throw(ArgumentError(
            "type-$(row.source_type) control source $(row.name) has no switch node",
        ))
        switch_index = findfirst(candidates) do switch
            switch.a == node_index || switch.b == node_index
        end
        switch_index === nothing && throw(ArgumentError(
            "type-$(row.source_type) control source $(row.name) does not identify a switch endpoint",
        ))
        push!(
            observations,
            ControlSwitchObservation(
                row.name,
                row.source_type,
                candidates[switch_index],
                something(
                    findfirst(
                        switch -> switch === candidates[switch_index],
                        ordinary_candidates,
                    ),
                    0,
                ),
            ),
        )
    end
    return observations
end

function _deck_control_system_network_runtime(
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
    network_elements::AbstractVector = Any[],
)
    source_rows = DeckParser.deck_control_system_source_rows(parsed)
    function_rows = DeckParser.deck_control_system_function_rows(parsed)
    expression_rows = DeckParser.deck_control_system_expression_rows(parsed)
    device_rows = DeckParser.deck_control_system_device_rows(parsed)
    output_names = _deck_control_system_trace_output_names(parsed)
    switch_rows = DeckParser.deck_control_system_switch_coupling_rows(parsed)
    switch_config_rows = DeckParser.deck_control_system_switch_rows(parsed)
    length(switch_rows) == length(switch_config_rows) ||
        throw(ArgumentError("control-system switch coupling count does not match switch-card count"))
    isempty(source_rows) && isempty(function_rows) && isempty(expression_rows) && isempty(device_rows) &&
        isempty(output_names) && isempty(switch_rows) && return nothing

    initial_values = _deck_control_system_source_values(
        parsed,
        Dict{Symbol,Float64}(),
        nothing,
        0.0,
    )
    windowed_constant_sources = WindowedConstantControlSignal[
        _deck_windowed_constant_control_signal(row)
        for row in source_rows
        if row.source_type == 11
    ]
    sinusoidal_sources = _deck_control_system_sinusoidal_sources(parsed)
    for source in sinusoidal_sources
        initial_values[source.name] = 0.0
    end
    waveform_sources = ControlWaveformSignal[
        _deck_control_waveform_signal(row, dt_s)
        for row in source_rows
        if 12 <= row.source_type <= 24 && row.source_type != 14
    ]
    for source in waveform_sources
        initial_values[source.name] = 0.0
    end
    functions = _deck_control_system_functions(parsed, nothing)
    devices = _deck_control_system_devices(parsed)
    for function_row in functions
        get!(initial_values, function_row.output_name, 0.0)
    end
    for row in expression_rows
        get!(initial_values, row.name, 0.0)
    end
    for device in devices
        get!(initial_values, device.name, 0.0)
    end
    for name in output_names
        get!(initial_values, name, 0.0)
    end
    initial_values[:ALWAYS_ENABLED] = 1.0
    initial_values[:ALWAYS_DISABLED] = 0.0
    signals = ConstantControlSignal[
        ConstantControlSignal(name, initial_values[name])
        for name in sort!(collect(keys(initial_values)); by = String)
    ]
    state = ControlSystemExecutionState(
        signals,
        AlgebraicControlAssignment[];
        functions = functions,
        devices = devices,
        output_names = output_names,
        deltat_s = dt_s,
        frequency_hz = _deck_control_system_frequency_hz(parsed),
    )
    control_expressions = if isempty(expression_rows)
        initialize_control_system_utilities!(state, 0, 0.0)
        ControlExpressionRuntime[]
    else
        advance_control_system_state!(state, 0, 0.0; execute_devices = false)
        _deck_control_system_expression_runtimes(parsed, state)
    end
    supplemental_devices = _control_system_supplemental_device_runtime(
        parsed,
        state,
        dt_s,
    )
    frequency_initializations = _control_system_frequency_initializations(
        state,
        supplemental_devices,
        sinusoidal_sources,
    )
    signal_slot_names = _deck_control_system_signal_slot_names(parsed)
    slot_count = length(signal_slot_names)
    for row in DeckParser.deck_over5a_source_rows(parsed)
        source_type = abs(row.iform)
        source_type == 17 || source_type >= 60 || continue
        slot = round(Int, row.sfreq)
        1 <= slot <= slot_count || throw(ArgumentError(
            "source row on line $(row.line_no) references control-system signal " *
            "slot $slot, but only $slot_count slots are defined",
        ))
    end
    plan = deck_over16_boundary_plan(parsed)
    for (slot, line_no) in zip(
        plan.source_tacs_override_xtcs_indices,
        plan.source_tacs_override_line_numbers,
    )
        0 <= slot <= slot_count || throw(ArgumentError(
            "source override on line $line_no references control-system signal " *
            "slot $slot, but only $slot_count slots are defined",
        ))
    end

    network_voltage_source_names = Symbol[]
    network_voltage_source_node_indices = Int[]
    for row in source_rows
        row.source_type == 90 || continue
        node_index = get(parsed.node_map, row.name, 0)
        node_index > 0 ||
            throw(ArgumentError("control-system voltage source $(row.name) has no network node"))
        push!(network_voltage_source_names, row.name)
        push!(network_voltage_source_node_indices, node_index)
    end

    switch_elements = TACSControlledSwitch[]
    switch_signal_names = Symbol[]
    switch_clamp_signal_names = Union{Missing,Symbol}[]
    switch_output_indices = Int[]
    switch_output_names = Symbol[]
    for (index, row) in enumerate(switch_rows)
        switch_config = switch_config_rows[index]
        row.from_node_index !== missing && row.to_node_index !== missing ||
            throw(ArgumentError("control-system switch on line $(row.line_no) has unresolved nodes"))
        for signal_name in (switch_config.gate_signal, switch_config.clamp_signal)
            signal_name === missing && continue
            haskey(state.values, signal_name) || throw(ArgumentError(
                "control-system switch on line $(row.line_no) references unresolved signal $signal_name",
            ))
        end
        initial_control =
            _control_system_switch_initially_closed(row.initial_state) ? 1.0 :
            get(state.values, row.control_signal, 0.0)
        clamp_control =
            row.switch_type in (11, 12) && switch_config.clamp_signal !== missing ?
            Ref(get(state.values, switch_config.clamp_signal, 0.0)) : nothing
        push!(
            switch_elements,
            TACSControlledSwitch(
                Int(row.from_node_index),
                Int(row.to_node_index),
                Ref(initial_control);
                threshold = nextfloat(0.0),
                unidirectional_latching = row.switch_type == 11,
                bidirectional_latching = row.switch_type == 12,
                clamp_control = clamp_control,
                initially_closed =
                    _control_system_switch_initially_closed(row.initial_state),
                ignition_voltage = switch_config.ignition_voltage === missing ?
                    0.0 : Float64(switch_config.ignition_voltage),
                holding_current = switch_config.holding_current === missing ?
                    0.0 : Float64(switch_config.holding_current),
                deionization_time_s = switch_config.deionization_time_s === missing ?
                    0.0 : Float64(switch_config.deionization_time_s),
                delayed_arc = switch_config.delayed_arc === nothing ?
                    nothing :
                    ControlledSwitchDelayedArcState(
                        switch_config.delayed_arc.current_coefficient,
                        switch_config.delayed_arc.current_exponent,
                        switch_config.delayed_arc.time_scale_s,
                        switch_config.delayed_arc.cutoff_current_a,
                    ),
            ),
        )
        push!(switch_signal_names, row.control_signal)
        push!(
            switch_clamp_signal_names,
            row.switch_type in (11, 12) ? switch_config.clamp_signal : missing,
        )
        if _switch_current_report_selected(row.output_code)
            push!(switch_output_indices, index)
            push!(
                switch_output_names,
                Symbol(String(row.from_node), ".", String(row.to_node)),
            )
        end
    end
    switch_observations = _deck_control_switch_observations(
        parsed,
        network_elements,
        switch_elements,
    )
    runtime = ControlSystemNetworkRuntime(
        state,
        signal_slot_names,
        network_voltage_source_names,
        network_voltage_source_node_indices,
        switch_observations,
        switch_elements,
        switch_signal_names,
        switch_clamp_signal_names,
        switch_output_indices,
        switch_output_names,
        output_names,
        windowed_constant_sources,
        sinusoidal_sources,
        waveform_sources,
        control_expressions,
        supplemental_devices,
        frequency_initializations,
        0,
        0,
        false,
        nothing,
    )
    return runtime
end

function _append_control_system_network_elements!(
    elements::Vector,
    element_names::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    dt_s::Float64,
)
    runtime = _deck_control_system_network_runtime(parsed, dt_s, elements)
    runtime === nothing && return nothing
    append!(elements, runtime.switch_elements)
    append!(element_names, _control_system_switch_element_names(parsed))
    return runtime
end

function _record_control_system_network_outputs!(
    output::AbstractMatrix{Float64},
    sample::Int,
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
    first_channel::Int,
)
    runtime = context.control_system_runtime
    runtime === nothing && return first_channel
    if context.step_index > 0
        ordinary_switch_currents =
            isempty(runtime.switch_observations) ?
            context.switch_current_step_values :
            _deck_time_switch_report_current_values!(
                context,
                voltage,
                context.t_s,
            )
        _advance_control_system_network_runtime!(
            runtime,
            context.step_index,
            context.t_s;
            network_voltage = voltage,
            ordinary_switch_currents,
            ordinary_switch_closed_flags = context.switch_closed_step_flags,
        )
    end

    next_channel = first_channel
    for switch_index in runtime.switch_output_indices
        switch = runtime.switch_elements[switch_index]
        from_voltage = switch.a == 0 ? 0.0 : voltage[switch.a]
        to_voltage = switch.b == 0 ? 0.0 : voltage[switch.b]
        output[next_channel, sample] =
            switch.last_conductance * (from_voltage - to_voltage)
        next_channel += 1
    end
    for name in runtime.control_output_names
        haskey(runtime.state.values, name) ||
            throw(ArgumentError("missing control-system runtime output $name"))
        output[next_channel, sample] = runtime.state.values[name]
        next_channel += 1
    end
    return next_channel
end

function _project_closed_controlled_switch_voltages!(
    runtime::ControlSystemNetworkRuntime,
    voltage::AbstractVector{<:Real},
)
    for switch in runtime.switch_elements
        switch.last_conductance == switch.on_conductance || continue
        common_voltage = if switch.a == 0 || switch.b == 0
            0.0
        else
            0.5 * (voltage[switch.a] + voltage[switch.b])
        end
        for index in eachindex(runtime.network_voltage_source_names)
            node_index = runtime.network_voltage_source_node_indices[index]
            if node_index == switch.a || node_index == switch.b
                runtime.state.values[runtime.network_voltage_source_names[index]] =
                    common_voltage
            end
        end
    end
    return runtime
end

function _deck_control_system_trace_output_names(parsed::DeckParser.DeckParseResult)
    names = Symbol[]
    for row in parsed.control_system_output_request_rows
        if row.all_signals
            append!(names, (function_row.name for function_row in parsed.control_system_function_rows))
            append!(names, (expression_row.name for expression_row in parsed.control_system_expression_rows))
            append!(names, (device_row.name for device_row in parsed.control_system_device_rows))
            append!(names, _CONTROL_SYSTEM_UTILITY_VALUES)
            append!(names, (source_row.name for source_row in parsed.control_system_source_rows))
        end
        append!(names, row.signal_names)
    end
    return unique(names)
end

function _deck_control_system_output_names(
    parsed::DeckParser.DeckParseResult,
    _,
)
    names = Symbol[]
    for row in parsed.control_system_output_request_rows
        append!(names, row.signal_names)
    end
    return names
end

function _deck_control_system_source_values(
    parsed::DeckParser.DeckParseResult,
    steady_values::Dict{Symbol,Float64},
    _,
    time_s::Float64,
)
    values = copy(steady_values)
    network_source_values = Dict{Symbol,Float64}()
    for row in parsed.over5a_source_rows
        source_type = abs(row.iform)
        source_type in 11:15 || continue
        network_source_values[row.node] = analytic_source_value(
            source_type,
            row.crest,
            row.time1,
            row.sfreq,
            row.tstart,
            row.tstop,
            time_s,
        )
    end
    for row in parsed.control_system_source_rows
        if row.source_type == 90
            values[row.name] = get(network_source_values, row.name, 0.0)
        elseif row.source_type == 14
            values[row.name] = sinusoidal_control_signal_value(
                _deck_control_system_sinusoidal_signal(row),
                time_s,
            )
        elseif row.source_type == 11 && row.amplitude !== missing
            values[row.name] = _windowed_constant_control_signal_value(
                _deck_windowed_constant_control_signal(row),
                time_s,
            )
        elseif row.amplitude !== missing
            values[row.name] = Float64(row.amplitude)
        else
            values[row.name] = get(steady_values, row.name, 0.0)
        end
    end
    return values
end

function _deck_control_system_functions(
    parsed::DeckParser.DeckParseResult,
    _,
)
    return ControlTransferFunction[
        ControlTransferFunction(
            row.name,
            SignedControlSignalTerm[
                SignedControlSignalTerm(term.name, term.polarity)
                for term in row.input_terms
            ];
            gain = row.gain,
            order = row.order,
            numerator_coefficients = row.numerator_coefficients,
            denominator_coefficients = row.denominator_coefficients,
            lower_limit = row.lower_limit,
            upper_limit = row.upper_limit,
            lower_limit_signal = row.lower_limit_signal,
            upper_limit_signal = row.upper_limit_signal,
        )
        for row in parsed.control_system_function_rows
    ]
end

_deck_control_system_assignments(::DeckParser.DeckParseResult, _) =
    AlgebraicControlAssignment[]

function _deck_control_system_devices(parsed::DeckParser.DeckParseResult)
    return ControlSystemDevice[
        ControlSystemDevice(
            row.name,
            row.device_type,
            row.first_input === missing ? :ZERO : row.first_input;
            tail_signal_names = row.tail_signal_names,
            parameter_values = row.parameter_values,
        )
        for row in parsed.control_system_device_rows
    ]
end

function _control_system_switch_initially_closed(state::AbstractString)
    token = uppercase(strip(String(state)))
    return token in ("1", "CLOSE", "CLOSED", "ON")
end

function _control_system_switch_element_name(row, index::Int)
    return Symbol(
        row.switch_type == 11 ? "thyristor_" : "controlled_switch_",
        String(row.from_node),
        "_",
        String(row.to_node),
        "_",
        string(index),
    )
end

function _deck_windowed_constant_control_signal(
    row::DeckParser.DeckControlSystemSourceRow,
)
    row.source_type == 11 || throw(ArgumentError(
        "control source $(row.name) is not windowed constant type 11",
    ))
    row.amplitude === missing && throw(ArgumentError(
        "type-11 control source $(row.name) requires an amplitude",
    ))
    return WindowedConstantControlSignal(
        row.name,
        Float64(row.amplitude),
        row.activation_start_time_s,
        row.activation_stop_time_s,
    )
end

function _windowed_constant_control_signal_value(
    source::WindowedConstantControlSignal,
    time_s::Real,
)
    time = Float64(time_s)
    active = source.start_time_s <= time < source.stop_time_s
    time == 0.0 && source.start_time_s >= 0.0 && (active = false)
    return active ? source.value : 0.0
end

function _deck_control_waveform_signal(
    row::DeckParser.DeckControlSystemSourceRow,
    dt_s::Real,
)
    12 <= row.source_type <= 24 && row.source_type != 14 ||
        throw(ArgumentError(
            "control source $(row.name) is not a supported waveform source",
        ))
    amplitude = row.amplitude === missing ? 0.0 : row.amplitude
    cycle =
        row.delay_or_time_constant === missing ? 0.0 :
        row.delay_or_time_constant
    pulse_width =
        row.phase_or_width === missing ? 0.0 : row.phase_or_width
    if row.source_type in (23, 24)
        cycle > 0.0 || throw(ArgumentError(
            "type-$(row.source_type) control source $(row.name) requires a positive cycle",
        ))
    end
    if row.source_type == 23
        pulse_width = max(Float64(pulse_width), Float64(dt_s))
        pulse_width <= cycle || throw(ArgumentError(
            "type-23 control source $(row.name) requires pulse width no greater than its cycle",
        ))
    end
    return ControlWaveformSignal(
        row.name,
        row.source_type,
        Float64(amplitude),
        Float64(cycle),
        Float64(pulse_width),
        row.activation_start_time_s,
        row.activation_stop_time_s,
    )
end

function _control_waveform_signal_value(
    source::ControlWaveformSignal,
    time_s::Real,
)
    time = Float64(time_s)
    source.start_time_s <= time < source.stop_time_s || return 0.0
    elapsed = time - source.start_time_s
    if source.waveform_type == 23
        phase = mod(elapsed, source.cycle_s)
        return phase < source.pulse_width_s ? source.amplitude : 0.0
    elseif source.waveform_type == 24
        phase = mod(elapsed, source.cycle_s)
        return source.amplitude * phase / source.cycle_s
    end
    return source.amplitude
end

function _deck_trace_control_system_source_values(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
    sample_index::Int,
)
    values = _deck_control_system_source_values(
        parsed,
        Dict{Symbol,Float64}(),
        nothing,
        trace.time_s[sample_index],
    )
    for row in parsed.control_system_source_rows
        row.source_type == 90 || continue
        node_index = get(trace.node_map, row.name, 0)
        node_index > 0 || continue
        values[row.name] = trace.voltage_pu[node_index, sample_index]
    end
    return values
end

function _deck_control_system_initial_signals(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    isempty(trace.time_s) && return ConstantControlSignal[]
    source_values = _deck_trace_control_system_source_values(parsed, trace, 1)
    return ConstantControlSignal[
        ConstantControlSignal(name, source_values[name])
        for name in sort!(collect(keys(source_values)); by = String)
    ]
end

function _deck_control_system_switch_current_values(
    switch_spec,
    control_values::Dict{Symbol,Float64},
    voltage::AbstractVector{<:Real},
)
    switch_spec === nothing && return zeros(Float64, 0)
    values = Vector{Float64}(undef, length(switch_spec.rows))
    voltage_values = Float64.(voltage)
    for (index, row) in enumerate(switch_spec.rows)
        closed =
            _control_system_switch_initially_closed(row.initial_state) ||
            get(control_values, row.control_signal, 0.0) > 0.0
        if closed
            voltage_difference =
                _deck_node_voltage(voltage_values, row.from_node_index) -
                _deck_node_voltage(voltage_values, row.to_node_index)
            values[index] = 1.0e9 * voltage_difference
        else
            values[index] = 0.0
        end
    end
    return values
end

function _deck_control_system_trace_outputs(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
    initial_network_voltage::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    runtime = _deck_control_system_network_runtime(parsed, trace.dt_s, parsed.elements)
    runtime === nothing && return nothing
    switch_count = length(runtime.switch_output_names)
    output_count = length(runtime.control_output_names)
    switch_count == 0 && output_count == 0 && return nothing
    initial_network_voltage === nothing ||
        _initialize_control_system_network_steady_state!(
            runtime,
            initial_network_voltage,
        )
    values = Matrix{Float64}(undef, switch_count + output_count, length(trace.time_s))
    for sample_index in eachindex(trace.time_s)
        voltage = @view trace.voltage_pu[:, sample_index]
        if sample_index > 1
            ordinary_switch_currents, ordinary_switch_closed_flags =
                _deck_control_system_trace_switch_samples(
                    parsed,
                    voltage,
                    trace.time_s[sample_index],
                    trace.dt_s,
                )
            _advance_control_system_network_runtime!(
                runtime,
                sample_index - 1,
                trace.time_s[sample_index];
                network_voltage = voltage,
                ordinary_switch_currents,
                ordinary_switch_closed_flags,
            )
        end
        switch_values = zeros(Float64, switch_count)
        for (output_index, switch_index) in enumerate(runtime.switch_output_indices)
            switch = runtime.switch_elements[switch_index]
            sync_controlled_switch!(switch)
            from_voltage = switch.a == 0 ? 0.0 : voltage[switch.a]
            to_voltage = switch.b == 0 ? 0.0 : voltage[switch.b]
            switch_values[output_index] =
                switch.last_conductance * (from_voltage - to_voltage)
        end
        values[1:switch_count, sample_index] .= switch_values
        if output_count > 0
            values[(switch_count + 1):(switch_count + output_count), sample_index] .=
                control_system_output_values(runtime.state)
        end
        for switch in runtime.switch_elements
            update!(switch, voltage, trace.dt_s)
        end
    end
    names = vcat(runtime.switch_output_names, runtime.control_output_names)
    return (output_names = names, values = values)
end

function _deck_control_system_trace_switch_samples(
    parsed::DeckParser.DeckParseResult,
    voltage::AbstractVector{<:Real},
    time_s::Float64,
    dt_s::Float64,
)
    rows = DeckParser.deck_over5_switch_rows(parsed)
    currents = zeros(Float64, length(rows))
    closed_flags = zeros(Int, length(rows))
    voltage_values = Float64.(voltage)
    elements = Tuple(parsed.elements)
    for (index, row) in enumerate(rows)
        closed = _deck_time_switch_closed_at(
            row.initially_closed,
            Float64(row.close_time_s),
            Float64(row.open_time_s),
            time_s,
        )
        closed_flags[index] = closed ? 1 : 0
        closed || continue
        voltage_difference =
            _deck_node_voltage(voltage_values, row.from_node_value) -
            _deck_node_voltage(voltage_values, row.to_node_value)
        currents[index] = Float64(row.on_conductance) * voltage_difference
        to_ground_current = _deck_switch_grounded_conductance_current(
            elements,
            row.to_node_value,
            voltage_values,
            dt_s,
        )
        if to_ground_current !== nothing
            currents[index] = to_ground_current
            continue
        end
        from_ground_current = _deck_switch_grounded_conductance_current(
            elements,
            row.from_node_value,
            voltage_values,
            dt_s,
        )
        from_ground_current === nothing ||
            (currents[index] = -from_ground_current)
    end
    return currents, closed_flags
end

function _append_deck_control_system_outputs(
    trace::DeckEMTTrace,
    parsed::DeckParser.DeckParseResult,
    initial_network_voltage::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    switch_spec = _deck_control_system_switch_current_spec(parsed)
    embedded_names = Symbol[]
    switch_spec === nothing || append!(embedded_names, switch_spec.output_names)
    append!(embedded_names, _deck_control_system_trace_output_names(parsed))
    if !isempty(embedded_names) && length(trace.output_channel_names) >= length(embedded_names)
        first_embedded = length(trace.output_channel_names) - length(embedded_names) + 1
        trace.output_channel_names[first_embedded:end] == embedded_names && return trace
    end
    outputs = _deck_control_system_trace_outputs(
        parsed,
        trace,
        initial_network_voltage,
    )
    outputs === nothing && return trace
    appended_extrema = _sampled_trace_extrema(outputs.values, trace.time_s)
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        vcat(copy(trace.output_channel_names), outputs.output_names),
        copy(trace.output_node_indices),
        vcat(copy(trace.output_pu), outputs.values),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        vcat(copy(trace.output_maximum_values), appended_extrema.maximum_values),
        vcat(copy(trace.output_maximum_times_s), appended_extrema.maximum_times_s),
        vcat(copy(trace.output_minimum_values), appended_extrema.minimum_values),
        vcat(copy(trace.output_minimum_times_s), appended_extrema.minimum_times_s),
    )
end

function _deck_requested_electrical_output_names(
    parsed::DeckParser.DeckParseResult,
)
    switch_names = DeckParser.deck_over5_switch_names(parsed)
    switch_codes = DeckParser.deck_over5_switch_output_codes(parsed)
    switch_count = length(switch_names)
    switch_voltage_indices =
        _deck_switch_voltage_report_indices(switch_codes, switch_count)
    switch_current_indices =
        _deck_switch_current_report_indices(switch_codes, switch_count)
    return vcat(
        DeckParser.deck_over16_output_channel_names(parsed),
        DeckParser.deck_over16_branch_voltage_output_names(parsed),
        _deck_switch_voltage_report_names(switch_names, switch_voltage_indices),
        _deck_switch_current_report_names(switch_names, switch_current_indices),
        DeckParser.deck_over16_branch_current_output_names(parsed),
    )
end

function _deck_requested_electrical_trace(
    parsed::DeckParser.DeckParseResult,
    trace::DeckEMTTrace,
)
    output_names = Symbol[]
    for (element_name, element) in zip(parsed.element_names, parsed.elements)
        trace_output_is_public(element) || continue
        trace_output_channel_names!(output_names, element_name, element)
    end
    frequency_dependent_line_names = vcat(
        DeckParser.deck_sampled_frequency_line_element_names(parsed),
        DeckParser.deck_semlyen_line_element_names(parsed),
        DeckParser.deck_rational_frequency_line_element_names(parsed),
    )
    for element_name in frequency_dependent_line_names
        channel_prefix = string(element_name, '_')
        append!(
            output_names,
            (
                channel for channel in trace.output_channel_names
                if startswith(String(channel), channel_prefix)
            ),
        )
    end
    append!(output_names, _deck_requested_electrical_output_names(parsed))
    unique!(output_names)
    output_indices = Int[]
    for name in output_names
        index = findfirst(==(name), trace.output_channel_names)
        index === nothing && throw(ArgumentError(
            "deck-requested output channel $name is missing from the runtime trace",
        ))
        push!(output_indices, index)
    end
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        output_names,
        copy(trace.output_node_indices),
        copy(trace.output_pu[output_indices, :]),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        copy(trace.output_maximum_values[output_indices]),
        copy(trace.output_maximum_times_s[output_indices]),
        copy(trace.output_minimum_values[output_indices]),
        copy(trace.output_minimum_times_s[output_indices]),
    )
end

function scheduled_trace(trace::DeckEMTTrace, step_indices::AbstractVector{<:Integer})
    sample_indices = Int[]
    sample_count = length(trace.time_s)
    for step in step_indices
        step_value = Int(step)
        step_value >= 0 || throw(ArgumentError("scheduled trace step must be nonnegative"))
        sample_index = step_value + 1
        sample_index <= sample_count ||
            throw(ArgumentError("scheduled trace step $step_value exceeds trace range"))
        push!(sample_indices, sample_index)
    end
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        isempty(sample_indices) ? 0.0 : trace.time_s[sample_indices[end]],
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s[sample_indices]),
        copy(trace.voltage_pu[:, sample_indices]),
        copy(trace.output_channel_names),
        copy(trace.output_node_indices),
        copy(trace.output_pu[:, sample_indices]),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        copy(trace.output_maximum_values),
        copy(trace.output_maximum_times_s),
        copy(trace.output_minimum_values),
        copy(trace.output_minimum_times_s),
    )
end

function run_deck_emt(
    lines;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    source::AbstractString = "deck",
    initial_voltage_source::Symbol = :none,
    saturated_transformer_branch_runtime_enabled::Bool = false,
    coupled_lumped_sequence_history_enabled::Bool = false,
    distributed_transposed_line_runtime_enabled::Bool = true,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    time_horizon::Symbol = :arguments,
    output_schedule::Symbol = :all_steps,
    synchronous_machine_output_runtime_enabled::Bool = false,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    parsed = DeckParser.parse_deck_lines(lines; source = source)
    return run_deck_emt(
        parsed;
        dt_s = dt_s,
        t_end_s = t_end_s,
        initial_voltage_source = initial_voltage_source,
        saturated_transformer_branch_runtime_enabled =
            saturated_transformer_branch_runtime_enabled,
        coupled_lumped_sequence_history_enabled =
            coupled_lumped_sequence_history_enabled,
        distributed_transposed_line_runtime_enabled =
            distributed_transposed_line_runtime_enabled,
        recorded_step_indices = recorded_step_indices,
        time_horizon = time_horizon,
        output_schedule = output_schedule,
        synchronous_machine_output_runtime_enabled =
            synchronous_machine_output_runtime_enabled,
        source_signal_provider = source_signal_provider,
    )
end

function _deck_uses_type11_zero_state_initialization(
    parsed::DeckParser.DeckParseResult,
)
    source_rows = DeckParser.deck_over5a_source_rows(parsed)
    return !isempty(source_rows) && all(row -> abs(row.iform) == 11, source_rows)
end

function _deck_uses_node_condition_initialization(
    parsed::DeckParser.DeckParseResult,
)
    return any(
        row -> row.condition_kind == :node_voltage_initial_condition &&
               row.node_index > 0,
        DeckParser.deck_node_initial_condition_rows(parsed),
    )
end

function _deck_uses_current_zero_switching(parsed::DeckParser.DeckParseResult)
    return !isempty(_deck_current_zero_switch_names(parsed))
end

function _deck_has_nonlinear_runtime_owners(
    parsed::DeckParser.DeckParseResult,
)
    return !isempty(DeckParser.deck_zinc_oxide_nonlinear_rows(parsed)) ||
           !isempty(DeckParser.deck_nonlinear_resistance_rows(parsed)) ||
           !isempty(DeckParser.deck_triggered_timed_resistance_rows(parsed)) ||
           !isempty(DeckParser.deck_switching_nonlinear_resistor_rows(parsed)) ||
           !isempty(DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)) ||
           !isempty(DeckParser.deck_hysteretic_inductor_rows(parsed)) ||
           !isempty(DeckParser.deck_arrester_nonlinear_rows(parsed))
end

function _deck_uses_switched_nonlinear_runtime(
    parsed::DeckParser.DeckParseResult,
)
    return _deck_has_nonlinear_runtime_owners(parsed) &&
           deck_over16_boundary_plan(parsed).switch_count > 0
end

function _deck_uses_controlled_switch_nonlinear_runtime(
    parsed::DeckParser.DeckParseResult,
)
    return _deck_has_nonlinear_runtime_owners(parsed) &&
           !isempty(DeckParser.deck_control_system_switch_coupling_rows(parsed))
end

function _deck_uses_dynamic_nonlinear_runtime(
    parsed::DeckParser.DeckParseResult,
)
    return _deck_uses_switched_nonlinear_runtime(parsed) ||
           _deck_uses_controlled_switch_nonlinear_runtime(parsed)
end

function _deck_uses_control_system_feedback_runtime(
    parsed::DeckParser.DeckParseResult,
)
    isempty(DeckParser.deck_control_system_source_rows(parsed)) &&
    isempty(DeckParser.deck_control_system_function_rows(parsed)) &&
    isempty(DeckParser.deck_control_system_expression_rows(parsed)) &&
    isempty(DeckParser.deck_control_system_device_rows(parsed)) && return false
    return any(
        row -> abs(row.iform) == 17 || abs(row.iform) >= 60,
        DeckParser.deck_over5a_source_rows(parsed),
    ) || deck_over16_boundary_plan(parsed).source_tacs_override_count > 0
end

function _deck_current_zero_steady_state_supported(
    parsed::DeckParser.DeckParseResult,
)
    source_rows = DeckParser.deck_over5a_source_rows(parsed)
    !isempty(source_rows) && all(row -> abs(row.iform) in (11, 14), source_rows) ||
        return false
    current_zero_names = _deck_current_zero_switch_names(parsed)
    return any(
        row -> row.name in current_zero_names &&
               uppercase(strip(row.closed_marker)) == "CLOSED",
        DeckParser.deck_over5_switch_rows(parsed),
    )
end

function _deck_uses_primary_mixed_nonlinear_runtime(
    parsed::DeckParser.DeckParseResult,
)
    _deck_has_nonlinear_runtime_owners(parsed) || return false
    deck_over16_boundary_plan(parsed).switch_count == 0 || return false
    isempty(DeckParser.deck_universal_machine_definition_rows(parsed)) || return false
    isempty(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)) ||
        return false
    isempty(parsed.control_system_function_rows) || return false
    isempty(parsed.control_system_expression_rows) || return false
    isempty(parsed.control_system_device_rows) || return false
    source_rows = DeckParser.deck_over5a_source_rows(parsed)
    return isempty(source_rows) || all(row -> abs(row.iform) in (11, 14), source_rows)
end

function _deck_emt_execution_result(
    trace::DeckEMTTrace,
    parsed::DeckParser.DeckParseResult,
    ::Val{false};
    context=nothing,
    nonlinear_run=nothing,
)
    return trace
end

function _deck_emt_execution_result(
    trace::DeckEMTTrace,
    parsed::DeckParser.DeckParseResult,
    ::Val{true};
    context=nothing,
    nonlinear_run=nothing,
)
    return DeckEMTExecution(
        trace,
        electromagnetic_terminal_state(
            parsed,
            trace;
            context,
            nonlinear_run,
        ),
        context === nothing ?
            SeriesRLCAlterationRecord[] :
            copy(context.series_rlc_alteration_records),
    )
end

function _run_deck_emt(
    parsed::DeckParser.DeckParseResult;
    capture_terminal_state::Val{C}=Val(false),
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    initial_voltage_source::Symbol = :none,
    saturated_transformer_branch_runtime_enabled::Bool = false,
    coupled_lumped_sequence_history_enabled::Bool = false,
    distributed_transposed_line_runtime_enabled::Bool = true,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    time_horizon::Symbol = :arguments,
    output_schedule::Symbol = :all_steps,
    synchronous_machine_output_runtime_enabled::Bool = false,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
) where {C}
    time_horizon in (:arguments, :deck) ||
        throw(ArgumentError("unsupported EMT time horizon $time_horizon"))
    output_schedule in (:all_steps, :print, :plot, :print_and_plot) ||
        throw(ArgumentError("unsupported EMT output schedule $output_schedule"))
    recorded_step_indices !== nothing && output_schedule != :all_steps &&
        throw(ArgumentError("recorded_step_indices and a deck output schedule cannot both be supplied"))
    if !isempty(DeckParser.deck_fixed_source_constraint_rows(parsed))
        parsed = apply_deck_fixed_source_load_flow(parsed).deck
    end
    timing =
        time_horizon == :deck ?
        deck_fixed_step_horizon(parsed) :
        (dt_s = dt_s, t_end_s = t_end_s)
    runtime_dt_s = Float64(timing.dt_s)
    runtime_t_end_s = Float64(timing.t_end_s)
    runtime_recorded_step_indices =
        recorded_step_indices === nothing && output_schedule != :all_steps ?
        deck_output_step_indices(
            parsed,
            runtime_dt_s,
            runtime_t_end_s;
            schedule = output_schedule,
        ) :
        recorded_step_indices
    universal_machine_definitions = [
        row for row in DeckParser.deck_universal_machine_definition_rows(parsed)
        if row.card_index == 1 && row.machine_type in 1:12
    ]
    universal_machine_outputs =
        DeckParser.deck_universal_machine_output_summary_rows(parsed)
    if time_horizon == :deck && length(universal_machine_definitions) == 1 &&
       length(universal_machine_outputs) == 1 &&
       only(universal_machine_outputs).machine_count == 1 &&
       initial_voltage_source == :none &&
       !synchronous_machine_output_runtime_enabled &&
       isempty(series_rlc_alterations) &&
       source_signal_provider isa IdentitySourceSignalProvider
        machine_horizon = run_deck_universal_machine_horizon(
            parsed;
            machine_index = universal_machine_definitions[1].machine_index,
            time_step_s = runtime_dt_s,
        )
        machine_trace = _deck_universal_machine_trace(parsed, machine_horizon)
        trace = runtime_recorded_step_indices === nothing ? machine_trace :
            scheduled_trace(machine_trace, runtime_recorded_step_indices)
        return _deck_emt_execution_result(
            trace,
            parsed,
            capture_terminal_state,
        )
    end
    synchronous_machine_indices = unique(
        row.machine_index
        for row in DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    )
    if synchronous_machine_output_runtime_enabled &&
       length(synchronous_machine_indices) == 1 &&
       initial_voltage_source in (:none, :synchronous_machine_terminals) &&
       isempty(series_rlc_alterations) &&
       source_signal_provider isa IdentitySourceSignalProvider
        synchronous_step_count = fixed_step_count(runtime_dt_s, runtime_t_end_s)
        synchronous_step_count > 0 || throw(ArgumentError(
            "synchronous-machine runtime requires at least one dynamic step",
        ))
        synchronous_horizon = run_deck_synchronous_machine_horizon(
            parsed;
            time_step_s = runtime_dt_s,
            dynamic_step_count = synchronous_step_count,
            saturated_transformer_branch_runtime_enabled,
            coupled_lumped_sequence_history_enabled,
            recorded_step_indices = runtime_recorded_step_indices,
        )
        isempty(synchronous_horizon.deferred_effects) || throw(ArgumentError(
            "synchronous-machine runtime has unresolved effects: " *
            join(string.(synchronous_horizon.deferred_effects), ','),
        ))
        return _deck_emt_execution_result(
            synchronous_horizon.trace,
            parsed,
            capture_terminal_state,
        )
    elseif synchronous_machine_output_runtime_enabled &&
           !isempty(synchronous_machine_indices)
        throw(ArgumentError(
            "synchronous-machine runtime configuration is not owned by the native Julia horizon",
        ))
    end
    saturated_transformer_dynamic_runtime =
        saturated_transformer_branch_runtime_enabled &&
        let intake = _deck_runtime_saturated_transformer_intake(parsed)
            intake !== nothing && !isempty(intake.breakpoints)
        end
    ideal_transformer_source_runtime = any(
        element -> element isa IdealTransformerVoltageConstraint,
        parsed.elements,
    )
    switched_nonlinear_dynamic_runtime =
        _deck_uses_dynamic_nonlinear_runtime(parsed)
    control_system_feedback_runtime =
        _deck_uses_control_system_feedback_runtime(parsed)
    if (
        distributed_transposed_line_runtime_enabled &&
        (
            _deck_has_dynamic_distributed_line(parsed) ||
            !isempty(DeckParser.deck_bergeron_line_rows(parsed)) ||
            _deck_has_frequency_dependent_line_runtime(parsed)
        ) ||
        saturated_transformer_dynamic_runtime ||
        ideal_transformer_source_runtime ||
        switched_nonlinear_dynamic_runtime ||
        control_system_feedback_runtime
    ) &&
       initial_voltage_source in (:none, :steady_state) &&
       !synchronous_machine_output_runtime_enabled
        saturated_transformer_intake =
            _deck_runtime_saturated_transformer_intake(parsed)
        augmented_node_count =
            saturated_transformer_intake === nothing ?
            length(parsed.node_map) :
            _saturated_transformer_runtime_node_count(
                parsed,
                saturated_transformer_intake,
            )
        boundary_state = _deck_dynamic_timestep_state(
            parsed,
            augmented_node_count,
            runtime_dt_s,
        )
        boundary_run = run_deck_emt_with_over16_boundary(
            parsed,
            boundary_state;
            dt_s = runtime_dt_s,
            t_end_s = runtime_t_end_s,
            steady_state_initial_sample_enabled =
                initial_voltage_source == :steady_state,
            saturated_transformer_intake = saturated_transformer_intake,
            saturated_transformer_nonlinear_current_enabled =
                saturated_transformer_intake !== nothing,
            coupled_lumped_sequence_history_enabled =
                coupled_lumped_sequence_history_enabled ||
                !isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)),
            time_switch_event_delay_s =
                time_horizon == :deck ? runtime_dt_s : 0.0,
            current_zero_switching =
                time_horizon == :deck &&
                _deck_uses_current_zero_switching(parsed),
            recorded_step_indices = runtime_recorded_step_indices,
            series_rlc_alterations = series_rlc_alterations,
            store_step_updates = switched_nonlinear_dynamic_runtime,
            source_signal_provider = source_signal_provider,
        )
        requested_trace = _deck_requested_electrical_trace(parsed, boundary_run.trace)
        if switched_nonlinear_dynamic_runtime
            requested_trace = _append_deck_nonlinear_outputs(
                requested_trace,
                boundary_run,
            )
        end
        initial_control_voltage =
            boundary_run.steady_state_initial_sample_applied ?
            boundary_run.steady_state_initial_voltage_values :
            nothing
        trace = _append_deck_control_system_outputs(
            requested_trace,
            parsed,
            initial_control_voltage,
        )
        return _deck_emt_execution_result(
            trace,
            parsed,
            capture_terminal_state;
            context = boundary_run.context,
            nonlinear_run = boundary_run,
        )
    end

    primary_mixed_nonlinear_runtime =
        (
            time_horizon == :deck ||
            !isempty(DeckParser.deck_triggered_timed_resistance_rows(parsed)) ||
            !isempty(DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed))
        ) &&
        initial_voltage_source == :none &&
        !saturated_transformer_branch_runtime_enabled &&
        !coupled_lumped_sequence_history_enabled &&
        !synchronous_machine_output_runtime_enabled &&
        source_signal_provider isa IdentitySourceSignalProvider &&
        _deck_uses_primary_mixed_nonlinear_runtime(parsed)
    if primary_mixed_nonlinear_runtime
        nonlinear_run = _run_primary_nonlinear_deck(
            parsed;
            dt_s = runtime_dt_s,
            t_end_s = runtime_t_end_s,
            recorded_step_indices = runtime_recorded_step_indices,
            source_signal_provider = source_signal_provider,
        )
        trace = _append_deck_nonlinear_outputs(nonlinear_run.trace, nonlinear_run)
        return _deck_emt_execution_result(
            trace,
            parsed,
            capture_terminal_state;
            context = nonlinear_run.context,
            nonlinear_run,
        )
    end
    isempty(DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)) ||
        throw(ArgumentError(
            "piecewise nonlinear inductors require the primary memoryless-network runtime",
        ))

    current_zero_switching =
        time_horizon == :deck &&
        initial_voltage_source in (:none, :zero, :steady_state) &&
        !saturated_transformer_branch_runtime_enabled &&
        !coupled_lumped_sequence_history_enabled &&
        _deck_uses_current_zero_switching(parsed)
    type11_zero_state_initialization =
        time_horizon == :deck &&
        initial_voltage_source == :none &&
        !current_zero_switching &&
        !saturated_transformer_branch_runtime_enabled &&
        !coupled_lumped_sequence_history_enabled &&
        _deck_uses_type11_zero_state_initialization(parsed)
    node_condition_initialization =
        time_horizon == :deck &&
        initial_voltage_source == :none &&
        !current_zero_switching &&
        !saturated_transformer_branch_runtime_enabled &&
        !coupled_lumped_sequence_history_enabled &&
        _deck_uses_node_condition_initialization(parsed)
    runtime_initial_voltage_source =
        current_zero_switching && initial_voltage_source == :none ?
        (_deck_current_zero_steady_state_supported(parsed) ? :steady_state : :zero) :
        node_condition_initialization ?
        :node_conditions :
        type11_zero_state_initialization ?
        :zero :
        initial_voltage_source
    initial_sample =
        _deck_runtime_initial_voltage_sample(parsed, runtime_initial_voltage_source)
    saturated_transformer_intake =
        saturated_transformer_branch_runtime_enabled ?
        _deck_runtime_saturated_transformer_intake(parsed) :
        nothing
    context =
        saturated_transformer_intake !== nothing ?
        saturated_transformer_branch_augmented_step_context(
            parsed,
            saturated_transformer_intake;
            dt_s = runtime_dt_s,
            t_end_s = runtime_t_end_s,
            transformer_branch_shunt_capacitance_rows =
                _deck_transformer_branch_shunt_capacitance_rows(parsed, nothing),
            include_coupled_lumped_sequence_history =
                coupled_lumped_sequence_history_enabled,
            recorded_step_indices = runtime_recorded_step_indices,
            source_signal_provider = source_signal_provider,
        ) :
        coupled_lumped_sequence_history_enabled &&
        !isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)) ?
        coupled_lumped_sequence_augmented_step_context(
            parsed;
            dt_s = runtime_dt_s,
            t_end_s = runtime_t_end_s,
            recorded_step_indices = runtime_recorded_step_indices,
            source_signal_provider = source_signal_provider,
        ) :
        initialize_step_context(
            parsed;
            dt_s = runtime_dt_s,
            t_end_s = runtime_t_end_s,
            recorded_step_indices = runtime_recorded_step_indices,
            time_switch_event_delay_s =
                time_horizon == :deck ? runtime_dt_s : 0.0,
            current_zero_switching = current_zero_switching,
            source_signal_provider = source_signal_provider,
        )
    configure_series_rlc_alterations!(context, series_rlc_alterations)
    initial_sample = _initial_voltage_sample_for_context(
        initial_sample,
        context.system.node_count,
    )
    initial_control_voltage =
        initial_sample !== nothing &&
        hasproperty(initial_sample, :node_voltage_values) ?
        initial_sample.node_voltage_values :
        initial_sample
    initial_control_voltage === nothing ||
        _initialize_control_system_network_steady_state!(
            context,
            initial_control_voltage,
        )
    trace = run_deck_emt(context; initial_voltage_sample = initial_sample)
    trace = _append_deck_control_system_outputs(
        trace,
        parsed,
        initial_sample === nothing ? nothing : initial_control_voltage,
    )
    return _deck_emt_execution_result(
        trace,
        parsed,
        capture_terminal_state;
        context,
    )
end

function run_deck_emt(parsed::DeckParser.DeckParseResult; kwargs...)
    return _run_deck_emt(parsed; capture_terminal_state = Val(false), kwargs...)
end

function run_deck_emt_execution(parsed::DeckParser.DeckParseResult; kwargs...)
    return _run_deck_emt(parsed; capture_terminal_state = Val(true), kwargs...)
end

function prepare_emt_study(
    lines;
    source::AbstractString = "deck",
    kwargs...,
)
    parsed = DeckParser.parse_deck_lines(lines; source = source)
    return prepare_emt_study(parsed; kwargs...)
end

function prepare_emt_study(
    parsed::DeckParser.DeckParseResult;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    initial_voltage_source::Symbol = :none,
    initial_voltage_sample = nothing,
    saturated_transformer_branch_runtime_enabled::Bool = false,
    coupled_lumped_sequence_history_enabled::Bool = false,
    distributed_transposed_line_runtime_enabled::Bool = true,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    time_horizon::Symbol = :arguments,
    output_schedule::Symbol = :all_steps,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
    external_current_injection_provider = nothing,
)
    parsed = deepcopy(parsed)
    time_horizon in (:arguments, :deck) ||
        throw(ArgumentError("unsupported EMT time horizon $time_horizon"))
    output_schedule in (:all_steps, :print, :plot, :print_and_plot) ||
        throw(ArgumentError("unsupported EMT output schedule $output_schedule"))
    recorded_step_indices !== nothing && output_schedule != :all_steps &&
        throw(ArgumentError("recorded_step_indices and a deck output schedule cannot both be supplied"))
    initial_voltage_source in (:none, :steady_state) ||
        throw(ArgumentError(
            "prepared dynamic EMT execution supports :none or :steady_state initialization",
        ))
    initial_voltage_sample === nothing || initial_voltage_source == :none ||
        throw(ArgumentError(
            "an explicit initial voltage sample cannot be combined with initial_voltage_source",
        ))
    timing =
        time_horizon == :deck ?
        deck_fixed_step_horizon(parsed) :
        (dt_s = dt_s, t_end_s = t_end_s)
    runtime_dt_s = Float64(timing.dt_s)
    runtime_t_end_s = Float64(timing.t_end_s)
    runtime_recorded_step_indices =
        recorded_step_indices === nothing && output_schedule != :all_steps ?
        deck_output_step_indices(
            parsed,
            runtime_dt_s,
            runtime_t_end_s;
            schedule = output_schedule,
        ) :
        recorded_step_indices
    saturated_transformer_intake =
        saturated_transformer_branch_runtime_enabled ?
        _deck_runtime_saturated_transformer_intake(parsed) :
        nothing
    saturated_transformer_dynamic_runtime =
        saturated_transformer_intake !== nothing &&
        !isempty(saturated_transformer_intake.breakpoints)
    ideal_transformer_source_runtime = any(
        element -> element isa IdealTransformerVoltageConstraint,
        parsed.elements,
    )
    dynamic_network_runtime =
        distributed_transposed_line_runtime_enabled &&
        (
            _deck_has_dynamic_distributed_line(parsed) ||
            !isempty(DeckParser.deck_bergeron_line_rows(parsed)) ||
            _deck_has_frequency_dependent_line_runtime(parsed)
        ) ||
        saturated_transformer_dynamic_runtime ||
        ideal_transformer_source_runtime ||
        _deck_uses_dynamic_nonlinear_runtime(parsed) ||
        _deck_uses_control_system_feedback_runtime(parsed) ||
        any(
            element -> element isa PowerSemiconductorSwitch ||
                element isa PowerSemiconductorBridgeLeg,
            parsed.elements,
        ) ||
        !isempty(series_rlc_alterations) ||
        initial_voltage_sample !== nothing ||
        external_current_injection_provider !== nothing
    dynamic_network_runtime || throw(ArgumentError(
        "prepared EMT execution currently requires the production dynamic network runtime",
    ))
    augmented_node_count =
        saturated_transformer_intake === nothing ?
        length(parsed.node_map) :
        _saturated_transformer_runtime_node_count(
            parsed,
            saturated_transformer_intake,
        )
    timestep_state = _deck_dynamic_timestep_state(
        parsed,
        augmented_node_count,
        runtime_dt_s,
    )
    runtime = _prepare_dynamic_deck_runtime(
        parsed,
        timestep_state;
        dt_s = runtime_dt_s,
        t_end_s = runtime_t_end_s,
        steady_state_initial_sample_enabled =
            initial_voltage_source == :steady_state,
        supplied_initial_sample = initial_voltage_sample,
        saturated_transformer_intake = saturated_transformer_intake,
        saturated_transformer_nonlinear_current_enabled =
            saturated_transformer_intake !== nothing,
        coupled_lumped_sequence_history_enabled =
            coupled_lumped_sequence_history_enabled ||
            !isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)),
        time_switch_event_delay_s =
            time_horizon == :deck ? runtime_dt_s : 0.0,
        current_zero_switching =
            time_horizon == :deck &&
            _deck_uses_current_zero_switching(parsed),
        recorded_step_indices = runtime_recorded_step_indices,
        series_rlc_alterations = series_rlc_alterations,
        store_step_updates = _deck_uses_dynamic_nonlinear_runtime(parsed),
        source_signal_provider = source_signal_provider,
        over16_step_configs = external_current_injection_provider,
    )
    return PreparedEMTStudy(runtime, parsed)
end

@generated function _tuple_element_alias(
    elements::T,
    element_index::Int,
    expected,
) where {T<:Tuple}
    comparisons = [
        :(element_index == $index && elements[$index] === expected)
        for index in 1:fieldcount(T)
    ]
    return foldl(
        (left, right) -> :($left || $right),
        comparisons;
        init = :(false),
    )
end

function _tuple_element_alias(
    elements::NodalElementSequence,
    element_index::Int,
    expected,
)
    1 <= element_index <= length(elements) || return false
    return _nodal_element_index_alias(
        elements.contiguous_type_batches,
        element_index,
        0,
        expected,
    )
end

_nodal_element_index_alias(::Tuple{}, ::Int, ::Int, _expected) = false

function _nodal_element_index_alias(
    batches::Tuple,
    element_index::Int,
    preceding_count::Int,
    expected,
)
    batch = first(batches)
    final_index = preceding_count + length(batch)
    if element_index <= final_index
        return @inbounds batch[element_index - preceding_count] === expected
    end
    return _nodal_element_index_alias(
        Base.tail(batches),
        element_index,
        final_index,
        expected,
    )
end

_nodal_element_identity_alias(::Tuple{}, _expected) = false

function _nodal_element_identity_alias(batches::Tuple, expected)
    for element in first(batches)
        element === expected && return true
    end
    return _nodal_element_identity_alias(Base.tail(batches), expected)
end

_nodal_source_signal_alias(::Tuple{}, _signal) = false

function _nodal_source_signal_alias(batches::Tuple, signal)
    for element in first(batches)
        hasproperty(element, :value) && getproperty(element, :value) === signal &&
            return true
    end
    return _nodal_source_signal_alias(Base.tail(batches), signal)
end

function _check_prepared_runtime_aliases(runtime::PreparedDynamicDeckRuntime)
    step_configs = runtime.step_configs
    if step_configs isa DynamicDeckStepConfigProvider
        step_configs.timestep_state === runtime.timestep_state ||
            throw(ArgumentError("prepared timestep config must own the workspace timestep state"))
        step_configs.plan === runtime.plan ||
            throw(ArgumentError("prepared timestep config must own the workspace boundary plan"))
    end
    context = runtime.context
    elements = context.system.elements
    history_plan = context.electromagnetic_history_plan
    for index in eachindex(history_plan.element_indices)
        element_index = history_plan.element_indices[index]
        kind = history_plan.kinds[index]
        batch_index = history_plan.batch_indices[index]
        batch = if kind == SERIES_RL_HISTORY
            history_plan.series_rl_branches
        elseif kind == SERIES_RLC_HISTORY
            history_plan.series_rlc_branches
        elseif kind == CAPACITOR_HISTORY
            history_plan.capacitor_branches
        elseif kind == COUPLED_INDUCTIVE_HISTORY
            history_plan.coupled_inductive_branches
        elseif kind == COUPLED_SERIES_RL_HISTORY
            history_plan.coupled_series_rl_branches
        elseif kind == BREQIV_HISTORY
            history_plan.breqiv_injections
        else
            throw(ArgumentError(
                "prepared electromagnetic history kind has no alias owner",
            ))
        end
        _tuple_element_alias(elements, element_index, batch[batch_index]) ||
            throw(ArgumentError(
                "prepared electromagnetic history plan must alias the nodal element",
            ))
    end
    for branch in context.saturated_transformer_nonlinear_slope_branch_batch
        _nodal_element_identity_alias(
            elements.contiguous_type_batches,
            branch,
        ) || throw(ArgumentError(
            "prepared nonlinear slope batch must alias the nodal element",
        ))
    end
    for signal in context.analytic_source_signals
        _nodal_source_signal_alias(
            elements.contiguous_type_batches,
            signal,
        ) || throw(ArgumentError(
            "prepared analytic source list must alias the nodal source signal",
        ))
    end
    source_runtime = context.source_function_runtime
    if source_runtime !== nothing
        source_runtime.plan === runtime.plan || throw(ArgumentError(
            "prepared source-function runtime must alias the boundary plan",
        ))
    end
    control_runtime = context.control_system_runtime
    if control_runtime !== nothing
        for observation in control_runtime.switch_observations
            _nodal_element_identity_alias(
                elements.contiguous_type_batches,
                observation.switch,
            ) ||
                throw(ArgumentError(
                    "prepared control observation must alias the nodal switch",
                ))
        end
        for switch in control_runtime.switch_elements
            _nodal_element_identity_alias(
                elements.contiguous_type_batches,
                switch,
            ) || throw(ArgumentError(
                "prepared control switch runtime must alias the nodal switch",
            ))
        end
    end
    return runtime
end

function _restore_prepared_dynamic_runtime!(
    runtime::PreparedDynamicDeckRuntime,
    template::PreparedDynamicDeckRuntime,
    restorer::TimestepStateRestorer,
)
    typeof(runtime) === typeof(template) || throw(ArgumentError(
        "prepared runtime and workspace types must match for in-place reset",
    ))
    restore_timestep_state!(runtime, template, restorer)
    return _check_prepared_runtime_aliases(runtime)
end

function EMTStudyWorkspace(prepared::PreparedEMTStudy{R,P}) where {R,P}
    runtime = _check_prepared_runtime_aliases(deepcopy(prepared.runtime_template))
    reduced_output_names = _deck_requested_electrical_output_names(prepared.parsed)
    switch_spec = _deck_control_system_switch_current_spec(prepared.parsed)
    switch_spec === nothing || append!(
        reduced_output_names,
        switch_spec.output_names,
    )
    append!(
        reduced_output_names,
        _deck_control_system_trace_output_names(prepared.parsed),
    )
    reduced_output_indices = Vector{Int}(undef, length(reduced_output_names))
    for index in eachindex(reduced_output_names)
        context_index = findfirst(
            ==(reduced_output_names[index]),
            runtime.context.output_channel_names,
        )
        context_index === nothing && throw(ArgumentError(
            "prepared reduced output channel $(reduced_output_names[index]) is missing from the context",
        ))
        reduced_output_indices[index] = context_index
    end
    source_names = runtime.plan.source_names
    source_signal_plan_indices = Vector{Int}(
        undef,
        length(runtime.context.analytic_source_names),
    )
    for signal_index in eachindex(source_signal_plan_indices)
        source_index = findfirst(
            ==(runtime.context.analytic_source_names[signal_index]),
            source_names,
        )
        source_signal_plan_indices[signal_index] =
            source_index === nothing ? 0 : source_index
    end
    return EMTStudyWorkspace{R,P}(
        runtime,
        prepared.parsed,
        reduced_output_indices,
        source_signal_plan_indices,
        0,
        0,
        true,
        :unselected,
        TimestepStateRestorer(),
    )
end

function reset_emt_study!(
    workspace::EMTStudyWorkspace{R,P},
    prepared::PreparedEMTStudy{R,P},
) where {R,P}
    _restore_prepared_dynamic_runtime!(
        workspace.runtime,
        prepared.runtime_template,
        workspace.reset_restorer,
    )
    workspace.parsed = prepared.parsed
    workspace.reset_count += 1
    workspace.ready = true
    workspace.execution_mode = :unselected
    return workspace
end

function EMTStudyBatch(
    prepared::PreparedEMTStudy;
    backend::AbstractEMTExecutionBackend=EMTCPUBackend(),
    workspace_count::Int=backend isa EMTCPUBackend && backend.threaded ?
        Threads.nthreads() : 1,
)
    workspace_count > 0 ||
        throw(ArgumentError("workspace_count must be positive"))
    workspaces = [EMTStudyWorkspace(prepared) for _ in 1:workspace_count]
    return EMTStudyBatch(workspaces, backend)
end

function emt_candidate_parameter_names(prepared::PreparedEMTStudy)
    return copy(prepared.runtime_template.plan.source_names)
end

function emt_candidate_parameter_names(workspace::EMTStudyWorkspace)
    return copy(workspace.runtime.plan.source_names)
end

function _append_emt_parameter_names!(
    names::Vector{Symbol},
    source_names::AbstractVector{Symbol},
    ::EMTSourceCrestCandidate,
)
    for source_name in source_names
        push!(names, Symbol(source_name, :_crest))
    end
    return names
end

function _append_emt_parameter_names!(
    names::Vector{Symbol},
    source_names::AbstractVector{Symbol},
    ::EMTSourceFrequencyRateParameter,
)
    for source_name in source_names
        push!(names, Symbol(source_name, :_frequency_or_rate))
    end
    return names
end

_append_emt_parameter_tuple_names!(
    names::Vector{Symbol},
    source_names::AbstractVector{Symbol},
    ::Tuple{},
) = names

function _append_emt_parameter_tuple_names!(
    names::Vector{Symbol},
    source_names::AbstractVector{Symbol},
    parameters::Tuple,
)
    _append_emt_parameter_names!(names, source_names, first(parameters))
    return _append_emt_parameter_tuple_names!(
        names,
        source_names,
        Base.tail(parameters),
    )
end

function emt_candidate_parameter_names(
    prepared::PreparedEMTStudy,
    candidate::EMTModelParameterCandidate,
)
    names = Symbol[]
    _append_emt_parameter_tuple_names!(
        names,
        prepared.runtime_template.plan.source_names,
        candidate.parameters,
    )
    return names
end

function _apply_source_crest_parameter!(
    workspace::EMTStudyWorkspace,
    candidate::EMTSourceCrestCandidate,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be reset before candidate mutation",
    ))
    runtime = workspace.runtime
    source_crests = runtime.plan.source_crest_values
    length(candidate.crest_values) == length(source_crests) || throw(ArgumentError(
        "candidate crest count must match prepared source count",
    ))
    for index in eachindex(source_crests, candidate.crest_values)
        crest = candidate.crest_values[index]
        isfinite(crest) ||
            throw(ArgumentError("source crest candidates must be finite"))
        source_crests[index] = crest
    end
    signals = runtime.context.analytic_source_signals
    signal_plan_indices = workspace.source_signal_plan_indices
    @inbounds for signal_index in eachindex(signals, signal_plan_indices)
        source_index = signal_plan_indices[signal_index]
        source_index == 0 && continue
        signals[signal_index].crest = source_crests[source_index]
    end
    return workspace
end

function apply_emt_parameter!(
    workspace::EMTStudyWorkspace,
    candidate::EMTSourceCrestCandidate,
)
    return _apply_source_crest_parameter!(workspace, candidate)
end

function apply_emt_parameter!(
    workspace::EMTStudyWorkspace,
    parameter::EMTSourceFrequencyRateParameter,
)
    source_rates = workspace.runtime.plan.source_sfreq_values
    values = parameter.frequency_or_rate_values
    length(values) == length(source_rates) || throw(ArgumentError(
        "source frequency-or-rate count must match prepared source count",
    ))
    @inbounds for index in eachindex(source_rates, values)
        rate = values[index]
        isfinite(rate) || throw(ArgumentError(
            "source frequency-or-rate parameters must be finite",
        ))
        source_rates[index] = rate
    end
    signals = workspace.runtime.context.analytic_source_signals
    signal_plan_indices = workspace.source_signal_plan_indices
    @inbounds for signal_index in eachindex(signals, signal_plan_indices)
        source_index = signal_plan_indices[signal_index]
        source_index == 0 && continue
        signals[signal_index].angular_frequency_or_rate = source_rates[source_index]
    end
    return workspace
end

function apply_emt_candidate!(
    workspace::EMTStudyWorkspace,
    parameter::AbstractEMTModelParameter,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be reset before candidate mutation",
    ))
    return apply_emt_parameter!(workspace, parameter)
end

_apply_emt_parameter_tuple!(workspace::EMTStudyWorkspace, ::Tuple{}) = workspace

function _apply_emt_parameter_tuple!(
    workspace::EMTStudyWorkspace,
    parameters::Tuple,
)
    apply_emt_parameter!(workspace, first(parameters))
    return _apply_emt_parameter_tuple!(workspace, Base.tail(parameters))
end

function apply_emt_candidate!(
    workspace::EMTStudyWorkspace,
    candidate::EMTModelParameterCandidate,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be reset before candidate mutation",
    ))
    return _apply_emt_parameter_tuple!(workspace, candidate.parameters)
end

function evaluate_emt_candidate!(
    workspace::EMTStudyWorkspace,
    prepared::PreparedEMTStudy,
    candidate::AbstractEMTStudyCandidate,
)
    workspace.ready || reset_emt_study!(workspace, prepared)
    apply_emt_candidate!(workspace, candidate)
    return evaluate_emt_study!(workspace)
end

function _check_reducer_channel(
    channel_index::Int,
    output_values::AbstractMatrix,
)
    1 <= channel_index <= size(output_values, 1) || throw(ArgumentError(
        "output reducer channel index is outside the trace",
    ))
    size(output_values, 2) > 0 ||
        throw(ArgumentError("output reducer requires at least one sample"))
    return nothing
end

function reduce_emt_trace(reducer::EMTOutputRMSReducer, trace)
    output_values = trace.output_pu
    channel = reducer.channel_index
    _check_reducer_channel(channel, output_values)
    squared_sum = 0.0
    @inbounds for sample in axes(output_values, 2)
        squared_sum += abs2(output_values[channel, sample])
    end
    return sqrt(squared_sum / size(output_values, 2))
end

function reduce_emt_trace(reducer::EMTOutputRMSReducer, context::EMTStepContext)
    output_values = context.output_pu
    channel = reducer.channel_index
    _check_reducer_channel(channel, output_values)
    squared_sum = 0.0
    @inbounds for sample in axes(output_values, 2)
        squared_sum += abs2(output_values[channel, sample])
    end
    return sqrt(squared_sum / size(output_values, 2))
end

function reduce_emt_trace(reducer::EMTOutputPeakReducer, trace)
    output_values = trace.output_pu
    channel = reducer.channel_index
    _check_reducer_channel(channel, output_values)
    peak = 0.0
    @inbounds for sample in axes(output_values, 2)
        peak = max(peak, abs(output_values[channel, sample]))
    end
    return peak
end

function reduce_emt_trace(reducer::EMTOutputPeakReducer, context::EMTStepContext)
    output_values = context.output_pu
    channel = reducer.channel_index
    _check_reducer_channel(channel, output_values)
    peak = 0.0
    @inbounds for sample in axes(output_values, 2)
        peak = max(peak, abs(output_values[channel, sample]))
    end
    return peak
end

function _mapped_context_output_channel(
    channel_index::Int,
    reduced_output_indices::AbstractVector{Int},
    output_values::AbstractMatrix,
)
    1 <= channel_index <= length(reduced_output_indices) || throw(ArgumentError(
        "output reducer channel index is outside the prepared public trace",
    ))
    context_channel = reduced_output_indices[channel_index]
    _check_reducer_channel(context_channel, output_values)
    return context_channel
end

function _reduce_emt_context(
    reducer::EMTOutputRMSReducer,
    context::EMTStepContext,
    reduced_output_indices::AbstractVector{Int},
)
    output_values = context.output_pu
    channel = _mapped_context_output_channel(
        reducer.channel_index,
        reduced_output_indices,
        output_values,
    )
    squared_sum = 0.0
    @inbounds for sample in axes(output_values, 2)
        squared_sum += abs2(output_values[channel, sample])
    end
    return sqrt(squared_sum / size(output_values, 2))
end

function _reduce_emt_context(
    reducer::EMTOutputPeakReducer,
    context::EMTStepContext,
    reduced_output_indices::AbstractVector{Int},
)
    output_values = context.output_pu
    channel = _mapped_context_output_channel(
        reducer.channel_index,
        reduced_output_indices,
        output_values,
    )
    peak = 0.0
    @inbounds for sample in axes(output_values, 2)
        peak = max(peak, abs(output_values[channel, sample]))
    end
    return peak
end

function evaluate_emt_reduced!(
    workspace::EMTStudyWorkspace,
    reducer::AbstractEMTTraceReducer,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be reset before another evaluation",
    ))
    workspace.ready = false
    _run_prepared_dynamic_deck!(
        workspace.runtime;
        collect_run_diagnostics = false,
    )
    workspace.evaluation_count += 1
    return _reduce_emt_context(
        reducer,
        workspace.runtime.context,
        workspace.reduced_output_indices,
    )
end

function _evaluate_emt_candidate_reduction!(
    results::AbstractVector,
    result_index::Int,
    workspace::EMTStudyWorkspace,
    prepared::PreparedEMTStudy,
    candidate::AbstractEMTStudyCandidate,
    reducer::AbstractEMTTraceReducer,
)
    workspace.ready || reset_emt_study!(workspace, prepared)
    apply_emt_candidate!(workspace, candidate)
    results[result_index] = evaluate_emt_reduced!(workspace, reducer)
    return nothing
end

function _prepare_emt_candidate!(
    workspace::EMTStudyWorkspace,
    prepared::PreparedEMTStudy,
    candidate::AbstractEMTStudyCandidate,
)
    workspace.ready || reset_emt_study!(workspace, prepared)
    apply_emt_candidate!(workspace, candidate)
    return workspace
end

function evaluate_emt_batch!(
    results::AbstractVector,
    batch::EMTStudyBatch{W,EMTCPUBackend},
    prepared::PreparedEMTStudy,
    candidates::AbstractVector{C},
    reducer::AbstractEMTTraceReducer,
) where {W,C<:AbstractEMTStudyCandidate}
    length(results) == length(candidates) || throw(ArgumentError(
        "batch result and candidate counts must match",
    ))
    isconcretetype(C) || throw(ArgumentError(
        "EMT batch candidates must have one concrete candidate type",
    ))
    isempty(candidates) && return results
    workspaces = batch.workspaces
    if batch.backend.threaded && Threads.nthreads() > 1 && length(workspaces) > 1
        workspace_count = length(workspaces)
        for wave_start in firstindex(candidates):workspace_count:lastindex(candidates)
            active_workspace_count = min(
                workspace_count,
                lastindex(candidates) - wave_start + 1,
            )
            # Reset and candidate mutation are setup operations. Keep them on
            # the caller thread and put a barrier before numerical execution;
            # this prevents mutable setup from overlapping another workspace's
            # timestep loop while retaining parallel solver work.
            for workspace_index in 1:active_workspace_count
                candidate_index = wave_start + workspace_index - 1
                _prepare_emt_candidate!(
                    workspaces[workspace_index],
                    prepared,
                    candidates[candidate_index],
                )
            end
            Threads.@threads :static for workspace_index in 1:active_workspace_count
                candidate_index = wave_start + workspace_index - 1
                results[candidate_index] = evaluate_emt_reduced!(
                    workspaces[workspace_index],
                    reducer,
                )
            end
        end
    else
        workspace = first(workspaces)
        for candidate_index in eachindex(candidates)
            _evaluate_emt_candidate_reduction!(
                results,
                candidate_index,
                workspace,
                prepared,
                candidates[candidate_index],
                reducer,
            )
        end
    end
    return results
end

function evaluate_emt_study!(workspace::EMTStudyWorkspace)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be reset before another evaluation",
    ))
    workspace.execution_mode === :unselected || throw(ArgumentError(
        "EMT study workspace is already owned by $(workspace.execution_mode) execution",
    ))
    workspace.execution_mode = :monolithic
    workspace.ready = false
    boundary_run = _run_prepared_dynamic_deck!(workspace.runtime)
    requested_trace =
        _deck_requested_electrical_trace(workspace.parsed, boundary_run.trace)
    if _deck_uses_dynamic_nonlinear_runtime(workspace.parsed)
        requested_trace = _append_deck_nonlinear_outputs(
            requested_trace,
            boundary_run,
        )
    end
    initial_control_voltage =
        workspace.runtime.steady_state_initial_sample === nothing ?
        nothing :
        workspace.runtime.steady_state_initial_sample.node_voltage_values
    trace = _append_deck_control_system_outputs(
        requested_trace,
        workspace.parsed,
        initial_control_voltage,
    )
    workspace.evaluation_count += 1
    return trace
end

function run_deck_emt(
    context::EMTStepContext;
    initial_voltage_sample = nothing,
    current_injection_samples = nothing,
)
    current_samples = _current_injection_samples_for_context(
        current_injection_samples,
        context.system.node_count,
        context.step_count + 1,
    )
    _apply_due_series_rlc_alterations!(context)
    if initial_voltage_sample !== nothing
        _seed_steady_state_network_state!(context, initial_voltage_sample)
        initial_output_values = _steady_state_initial_output_values(context)
        _apply_steady_state_initial_sample!(
            context,
            initial_voltage_sample,
            initial_output_values,
        )
        context.step_index = 1
        context.t_s = min(context.step_index, context.step_count) * context.dt_s
    end
    while context.step_index <= context.step_count
        if current_samples === nothing
            step!(context)
        else
            step!(context, @view current_samples[:, context.step_index + 1])
        end
    end
    return deck_trace(context)
end

function final_voltage_pu(trace::DeckEMTTrace, node::Symbol)
    index = get(trace.node_map, node, 0)
    index == 0 && throw(ArgumentError("unknown node $node"))
    return trace.voltage_pu[index, end]
end

function final_output_pu(trace::DeckEMTTrace, channel::Symbol)
    index = findfirst(==(channel), trace.output_channel_names)
    index === nothing && throw(ArgumentError("unknown output channel $channel"))
    return trace.output_pu[index, end]
end

function _last_over16_accepted_update(boundary_run)
    isempty(boundary_run.over16_updates) && return nothing
    return boundary_run.over16_updates[end].over16_update
end

_nested_value_count(values::AbstractVector) = mapreduce(length, +, values; init = 0)

function _deck_over16_output_request_line_number_count(run)
    return length(run.deck_over16_output_channel_line_numbers) +
           length(run.deck_over16_branch_voltage_output_line_numbers) +
           length(run.deck_over16_branch_current_output_line_numbers) +
           length(run.deck_over16_branch_power_output_line_numbers)
end

function _last_over16_tacs_or_csup_update(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        if update.tacs_result !== nothing ||
           update.ntacs3_result !== nothing ||
           update.tacs_linear_result !== nothing ||
           update.tacs_post_solve_result !== nothing ||
           update.csup_pre_solve_result !== nothing ||
           update.csup_termination_result !== nothing
            return update
        end
    end
    return nothing
end

function _last_over16_ntacs3_result(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        result = update.ntacs3_result
        result !== nothing && return result
    end
    return nothing
end

function _last_over16_output_report_result(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        update.output_report_result === nothing || return update.output_report_result
    end
    return nothing
end

function _last_over16_branch_power_result(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return nothing
    return output_result.power_result
end

function _over16_output_report_values(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        update.output_report_result === nothing || return Float64.(update.output_volti_values)
    end
    return Float64[]
end

function _over16_voltage_output_values(boundary_run)
    values = _over16_output_report_values(boundary_run)
    count = boundary_run.deck_over16_output_channel_count
    count > 0 || return Float64[]
    length(values) >= count + 1 || return Float64[]
    return Float64[values[index] for index in 2:(count + 1)]
end

function _boundary_voltage_output_values(boundary_run)
    if hasproperty(boundary_run, :steady_state_initial_sample_applied) &&
       boundary_run.steady_state_initial_sample_applied
        return Float64.(boundary_run.steady_state_initial_output_voltage_values)
    end
    return _over16_voltage_output_values(boundary_run)
end

function _over16_output_energy_state_values(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return Float64[]
    return Float64.(output_result.updated_bnrg_values)
end

function _over16_branch_output_bvalue_values(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return Float64[]
    return Float64.(output_result.updated_bvalue_values)
end

function _over16_branch_voltage_output_values(boundary_run)
    bvalues = _over16_branch_output_bvalue_values(boundary_run)
    count = boundary_run.deck_over16_branch_voltage_output_count
    length(bvalues) >= count || return Float64[]
    return Float64[bvalues[index] for index in 1:count]
end

function _over16_branch_current_output_values(boundary_run)
    bvalues = _over16_branch_output_bvalue_values(boundary_run)
    voltage_count = boundary_run.deck_over16_branch_voltage_output_count
    current_count = boundary_run.deck_over16_branch_current_output_count
    end_index = voltage_count + current_count
    length(bvalues) >= end_index || return Float64[]
    return Float64[bvalues[index] for index in (voltage_count + 1):end_index]
end

function _over16_branch_power_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.power_values)
end

function _over16_branch_energy_output_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.energy_output_values)
end

function _over16_branch_power_voltage_factor_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.voltage_factor_values)
end

function _over16_branch_power_current_factor_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.current_factor_values)
end

function _over16_branch_power_voltage_selectors(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Int[]
    return Int.(power_result.voltage_selectors)
end

function _over16_branch_power_current_selectors(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Int[]
    return Int.(power_result.current_selectors)
end

function _over16_tacs_linear_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_linear_xtcs_values)
end

function _over16_tacs_linear_rsblk_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_rsblk_values)
end

function _over16_tacs_post_solve_parsup_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_parsup_values)
end

function _over16_tacs_utility_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_xtcs_values)
end

function _over16_ntacs3_elec_output_result(boundary_run)
    result = _last_over16_ntacs3_result(boundary_run)
    result === nothing && return nothing
    return result.elec_output_copy_result
end

function _over16_ntacs3_elec_output_values(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Float64[]
    return Float64.(result.values)
end

function _over16_ntacs3_elec_output_update_indices(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    return Int.(result.xtcs_update_indices)
end

function _over16_elec_tacsto_report_kind_symbols(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Symbol[]
    hasproperty(result, :elec_report_records) || return Symbol[]
    return Symbol[record.kind for record in result.elec_report_records]
end

function _over16_elec_tacsto_report_has_value(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_report_records) || return Int[]
    return Int[record.has_value ? 1 : 0 for record in result.elec_report_records]
end

function _over16_elec_tacsto_report_values(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Float64[]
    hasproperty(result, :elec_report_records) || return Float64[]
    return Float64[record.has_value ? record.value : 0.0 for record in result.elec_report_records]
end

function _over16_elec_tacsto_error_stpflg(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_error_records) || return Int[]
    return Int[record.stpflg for record in result.elec_error_records]
end

function _over16_elec_tacsto_error_labels(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_error_records) || return Int[]
    return Int[record.fortran_label for record in result.elec_error_records]
end

function _over16_csup_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_csup_xtcs_values)
end

function _over16_csup_parsup_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_csup_parsup_values)
end

function _last_over16_switch_update(boundary_run)
    for step_update in Iterators.reverse(boundary_run.over16_updates)
        update = step_update.over16_update
        if update.switch_state_mutated ||
           update.switch_scan_result !== nothing ||
           update.switch_operation_result !== nothing ||
           update.switch_status_result !== nothing ||
           update.switch_order_result !== nothing ||
           update.switch_admittance_result !== nothing ||
           update.switch_current_result !== nothing ||
           update.switch_post_current_result !== nothing ||
           update.switch_bvalue_result !== nothing ||
           update.switch_alteration_result !== nothing ||
           update.switch_retriangularization_result !== nothing ||
           update.switch_sparse_factor_result !== nothing ||
           update.switch_fortran_sparse_factor_result !== nothing ||
           update.switch_network_solution_result !== nothing
            return update
        end
    end
    return nothing
end

function _over16_switch_closed_mask_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int[value ? 1 : 0 for value in update.output_switch_closed_mask]
end

function _over16_switch_nextsw_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_nextsw)
end

function _over16_switch_kode_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_kode)
end

function _over16_switch_current_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_currents)
end

function _over16_switch_current_product_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_current_products)
end

function _over16_switch_position_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_positions)
end

function _over16_switch_operation_modswt_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_operation_modswt)
end

function _over16_switch_post_current_modswt_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_post_current_modswt)
end

function _over16_switch_bvalue_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_bvalue)
end

function _over16_switch_admittance_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(vec(update.output_switch_admittance))
end

function _over16_switch_conductance_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_conductances)
end

function _over16_switch_factor_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(vec(update.output_switch_factor))
end

function _over16_switch_factor_pivot_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_factor_pivots)
end

function _over16_switch_sparse_factor_km_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_sparse_factor_km)
end

function _over16_switch_sparse_factor_ykm_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_sparse_factor_ykm)
end

function _over16_switch_sparse_factor_kk_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_sparse_factor_kk)
end

function _over16_switch_fortran_sparse_factor_km_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_fortran_sparse_factor_km)
end

function _over16_switch_fortran_sparse_factor_ykm_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_fortran_sparse_factor_ykm)
end

function _over16_switch_fortran_sparse_factor_kk_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_fortran_sparse_factor_kk)
end

function _over16_switch_rhs_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_rhs)
end

function _over16_switch_network_solution_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_network_solution)
end

function _over16_sparse_switch_state_flow_results(boundary_run)
    flows = Any[]
    for update in boundary_run.over16_updates
        flow = _over16_step_sparse_switch_state_flow_result(update)
        flow === nothing || push!(flows, flow)
    end
    return flows
end

function _over16_sparse_switch_state_flow_pass_count_values(boundary_run)
    return Int[flow.pass_count for flow in _over16_sparse_switch_state_flow_results(boundary_run)]
end

function _over16_sparse_switch_state_flow_pass_total(boundary_run)
    return sum(_over16_sparse_switch_state_flow_pass_count_values(boundary_run); init = 0)
end

function _over16_sparse_switch_state_flow_repeated_after_post_current_count(boundary_run)
    return sum(
        Int[flow.repeated_after_post_current_count for
            flow in _over16_sparse_switch_state_flow_results(boundary_run)];
        init = 0,
    )
end

function _over16_sparse_switch_state_flow_pass_lengths(boundary_run, field::Symbol)
    lengths = Int[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for values in getproperty(flow, field)
            push!(lengths, length(values))
        end
    end
    return lengths
end

function _over16_sparse_switch_state_flow_pass_int_values(boundary_run, field::Symbol)
    values = Int[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for pass_values in getproperty(flow, field)
            append!(values, Int.(pass_values))
        end
    end
    return values
end

function _over16_sparse_switch_state_flow_pass_float_values(boundary_run, field::Symbol)
    values = Float64[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for pass_values in getproperty(flow, field)
            append!(values, Float64.(pass_values))
        end
    end
    return values
end

function _is_over16_nonlinear_update(update)
    return update.nonlinear_inverse_state_mutated ||
           update.nonlinear_source_column_state_mutated ||
           update.nonlinear_current_state_mutated ||
           update.nonlinear_source_column_result !== nothing ||
           update.nonlinear_inverse_result !== nothing ||
           update.nonlinear_current_result !== nothing
end

function _over16_step_owner_updates(step_update)
    flow = _over16_step_sparse_switch_state_flow_result(step_update)
    return flow === nothing ? (step_update.over16_update,) : flow.pass_updates
end

function _last_over16_nonlinear_update(boundary_run)
    for step_update in Iterators.reverse(boundary_run.over16_updates)
        for update in Iterators.reverse(_over16_step_owner_updates(step_update))
            _is_over16_nonlinear_update(update) && return update
        end
    end
    return nothing
end

function _over16_nonlinear_state_mutation_count(boundary_run)
    count = 0
    for step_update in boundary_run.over16_updates
        for update in _over16_step_owner_updates(step_update)
            count += (
                update.nonlinear_inverse_state_mutated ||
                update.nonlinear_source_column_state_mutated ||
                update.nonlinear_current_state_mutated
            ) ? 1 : 0
        end
    end
    return count
end

function _over16_nonlinear_inverse_znonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_znonl)
end

function _over16_nonlinear_inverse_anonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_anonl)
end

function _over16_nonlinear_inverse_voltbc_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_voltbc)
end

function _over16_nonlinear_inverse_vzero_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_vzero)
end

function _over16_nonlinear_inverse_vnonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_vnonl)
end

function _over16_nonlinear_inverse_ilast_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_nonlinear_inverse_ilast)
end

function _over16_nonlinear_current_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_current_values)
end

function _over16_nonlinear_cursub_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_cursub_values)
end

function _over16_nonlinear_vchar_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_vchar_values)
end

function _last_over16_nonlinear_current_result(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return nothing
    update.nonlinear_current_result === nothing && return nothing
    return update.nonlinear_current_result
end

function _over16_nonlinear_current_results(boundary_run)
    results = Any[]
    for step_update in boundary_run.over16_updates
        for update in _over16_step_owner_updates(step_update)
            update.nonlinear_current_result === nothing && continue
            push!(results, update.nonlinear_current_result)
        end
    end
    return results
end

function _over16_nonlinear_current_result_counts(boundary_run, field::Symbol)
    return Int[Int(getproperty(result, field))
               for result in _over16_nonlinear_current_results(boundary_run)]
end

function _over16_nonlinear_current_result_count(boundary_run, field::Symbol)
    return sum(_over16_nonlinear_current_result_counts(boundary_run, field); init = 0)
end

function _over16_nonlinear_current_result_lengths(boundary_run, field::Symbol)
    return Int[length(getproperty(result, field))
               for result in _over16_nonlinear_current_results(boundary_run)]
end

function _over16_nonlinear_current_result_float_values(boundary_run, field::Symbol)
    values = Float64[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(values, Float64.(getproperty(result, field)))
    end
    return values
end

function _over16_nonlinear_current_mutation_order(boundary_run)
    mutation_order = Symbol[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(mutation_order, result.mutation_order)
    end
    return mutation_order
end

function _over16_nonlinear_current_deferred_calls(boundary_run)
    deferred_calls = Symbol[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(deferred_calls, result.deferred_calls)
    end
    return deferred_calls
end

function _deck_time_switch_closed_at(initially_closed::Bool,
                                     close_time_s::Float64,
                                     open_time_s::Float64,
                                     time_s::Float64)
    closed = initially_closed
    if isfinite(close_time_s) && time_s >= close_time_s
        closed = true
    end
    if isfinite(open_time_s) && time_s >= open_time_s
        closed = false
    end
    return closed
end

function _deck_time_switch_final_closed_flags(boundary_run)
    final_time_s = boundary_run.trace.time_s[end]
    return Int[
        _deck_time_switch_closed_at(
            boundary_run.deck_time_switch_initially_closed_flags[index],
            boundary_run.deck_time_switch_close_time_s_values[index],
            boundary_run.deck_time_switch_open_time_s_values[index],
            final_time_s,
        ) ? 1 : 0
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_final_conductance_values(boundary_run,
                                                    final_closed_flags::Vector{Int})
    return Float64[
        final_closed_flags[index] == 1 ?
        boundary_run.deck_time_switch_on_conductance_values[index] :
        boundary_run.deck_time_switch_off_conductance_values[index]
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_switch_grounded_conductance_current(
    boundary_run,
    endpoint::Int,
    voltage::AbstractVector{<:Real},
)
    hasproperty(boundary_run, :context) ||
        throw(ArgumentError("deck run does not expose its typed timestep context"))
    return _deck_switch_grounded_conductance_current(
        boundary_run.context,
        endpoint,
        voltage,
    )
end

function _deck_switch_grounded_conductance_current(
    context::EMTStepContext,
    endpoint::Int,
    voltage::AbstractVector{<:Real},
)
    endpoint == 0 && return nothing
    voltage_values = _deck_float64_voltage_values(voltage)
    return _deck_switch_grounded_conductance_current(
        context.system.elements,
        endpoint,
        voltage_values,
        context.dt_s,
    )
end

_deck_switch_grounded_conductance_current(
    ::Tuple{},
    ::Int,
    ::AbstractVector{Float64},
    ::Float64,
) = nothing

function _deck_switch_grounded_conductance_current(
    elements::Tuple,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    element = first(elements)
    if applicable(branch_companion_snapshot, element, voltage, dt_s)
        snapshot = branch_companion_snapshot(element, voltage, dt_s)
        if snapshot !== nothing && snapshot.kind == :conductance &&
           (
               (snapshot.a == endpoint && snapshot.b == 0) ||
               (snapshot.b == endpoint && snapshot.a == 0)
           )
            return snapshot.conductance * _deck_node_voltage(voltage, endpoint)
        end
    end
    return _deck_switch_grounded_conductance_current(
        Base.tail(elements),
        endpoint,
        voltage,
        dt_s,
    )
end


function _deck_switch_grounded_conductance_current(
    elements::NodalElementSequence,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    return _deck_switch_grounded_conductance_current_batches(
        elements.contiguous_type_batches,
        endpoint,
        voltage,
        dt_s,
    )
end

_deck_switch_grounded_conductance_current_batches(
    ::Tuple{},
    ::Int,
    ::AbstractVector{Float64},
    ::Float64,
) = nothing

function _deck_switch_grounded_conductance_current_batches(
    batches::Tuple,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    for element in first(batches)
        if applicable(branch_companion_snapshot, element, voltage, dt_s)
            snapshot = branch_companion_snapshot(element, voltage, dt_s)
            if snapshot !== nothing && snapshot.kind == :conductance &&
               (
                   (snapshot.a == endpoint && snapshot.b == 0) ||
                   (snapshot.b == endpoint && snapshot.a == 0)
               )
                return snapshot.conductance *
                    _deck_node_voltage(voltage, endpoint)
            end
        end
    end
    return _deck_switch_grounded_conductance_current_batches(
        Base.tail(batches),
        endpoint,
        voltage,
        dt_s,
    )
end

function _deck_closed_switch_grounded_path_current(
    boundary_run,
    switch_index::Int,
    voltage::AbstractVector{<:Real},
)
    from_node = boundary_run.deck_time_switch_from_node_indices[switch_index]
    to_node = boundary_run.deck_time_switch_to_node_indices[switch_index]
    to_ground_current =
        _deck_switch_grounded_conductance_current(boundary_run, to_node, voltage)
    to_ground_current !== nothing && return to_ground_current
    from_ground_current =
        _deck_switch_grounded_conductance_current(boundary_run, from_node, voltage)
    from_ground_current !== nothing && return -from_ground_current
    return nothing
end

_deck_node_voltage(voltage::AbstractVector{<:Real}, index::Int) =
    index == 0 ? 0.0 : voltage[index]

_deck_float64_voltage_values(voltage::AbstractVector{Float64}) = voltage
_deck_float64_voltage_values(voltage::AbstractVector{<:Real}) = Float64.(voltage)

function _deck_final_voltage_values(boundary_run)
    required_node_count = max(
        maximum(boundary_run.deck_time_switch_from_node_indices; init = 0),
        maximum(boundary_run.deck_time_switch_to_node_indices; init = 0),
    )
    update = _last_over16_accepted_update(boundary_run)
    if update !== nothing && hasproperty(update, :output_e_values) &&
       length(update.output_e_values) >= required_node_count
        return Float64.(update.output_e_values)
    end
    return Float64.(boundary_run.trace.voltage_pu[:, end])
end

function _deck_time_switch_final_voltage_values(boundary_run)
    final_voltage = _deck_final_voltage_values(boundary_run)
    return Float64[
        _deck_node_voltage(final_voltage, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(final_voltage, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_report_voltage_values(boundary_run, voltage::AbstractVector{<:Real})
    voltage_values = _deck_float64_voltage_values(voltage)
    return Float64[
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_final_current_values(final_voltage_values::Vector{Float64},
                                                final_conductance_values::Vector{Float64})
    return Float64[
        final_conductance_values[index] * final_voltage_values[index]
        for index in 1:length(final_voltage_values)
    ]
end

function _deck_current_extinction_report_state(element)
    if element isa CurrentZeroSwitch
        return (
            closed = element.closed,
            current_initialized = element.current_initialized,
            previous_current = element.previous_current,
        )
    elseif element isa TimeSwitch && element.current_extinction !== nothing
        state = element.current_extinction
        return (
            closed = state.closed,
            current_initialized = state.current_initialized,
            previous_current = state.previous_current,
        )
    end
    return nothing
end

function _deck_runtime_switch_element(boundary_run, switch_index::Int)
    owner = hasproperty(boundary_run, :context) ? boundary_run.context : boundary_run
    hasproperty(owner, :system) || return nothing
    switch_index >= 1 || return nothing
    return _deck_runtime_switch_element(owner.system.elements, switch_index, 0)
end

_deck_runtime_switch_element(::Tuple{}, ::Int, ::Int) = nothing

function _deck_runtime_switch_element(
    elements::Tuple,
    switch_index::Int,
    found::Int,
)
    element = first(elements)
    next_found =
        element isa Union{TimeSwitch,CurrentZeroSwitch} ? found + 1 : found
    if next_found == switch_index && next_found != found
        return element
    end
    return _deck_runtime_switch_element(
        Base.tail(elements),
        switch_index,
        next_found,
    )
end


function _deck_runtime_switch_element(
    elements::NodalElementSequence,
    switch_index::Int,
    found::Int,
)
    return _deck_runtime_switch_element_batches(
        elements.contiguous_type_batches,
        switch_index,
        found,
    )
end

_deck_runtime_switch_element_batches(::Tuple{}, ::Int, ::Int) = nothing

function _deck_runtime_switch_element_batches(
    batches::Tuple,
    switch_index::Int,
    found::Int,
)
    batch = first(batches)
    next_found = found
    if eltype(batch) <: Union{TimeSwitch,CurrentZeroSwitch}
        next_found += length(batch)
        if found < switch_index <= next_found
            return @inbounds batch[switch_index - found]
        end
    end
    return _deck_runtime_switch_element_batches(
        Base.tail(batches),
        switch_index,
        next_found,
    )
end

function _deck_time_switch_report_current_values!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
    time_s::Real,
)
    closed_flags = context.switch_closed_step_flags
    conductance_values = context.switch_conductance_step_values
    voltage_differences = context.switch_voltage_step_values
    current_values = context.switch_current_step_values
    time = Float64(time_s)
    for index in 1:context.deck_time_switch_count
        closed = _deck_time_switch_closed_at(
            context.deck_time_switch_initially_closed_flags[index],
            context.deck_time_switch_close_time_s_values[index],
            context.deck_time_switch_open_time_s_values[index],
            time,
        )
        element = _deck_runtime_switch_element(context, index)
        extinction = _deck_current_extinction_report_state(element)
        extinction === nothing || (closed = extinction.closed)
        closed_flags[index] = closed ? 1 : 0
        conductance = closed ?
            context.deck_time_switch_on_conductance_values[index] :
            context.deck_time_switch_off_conductance_values[index]
        conductance_values[index] = conductance
        voltage_difference =
            _deck_node_voltage(voltage, context.deck_time_switch_from_node_indices[index]) -
            _deck_node_voltage(voltage, context.deck_time_switch_to_node_indices[index])
        voltage_differences[index] = voltage_difference
        current_values[index] = conductance * voltage_difference
    end
    for index in 1:context.deck_time_switch_count
        element = _deck_runtime_switch_element(context, index)
        extinction = _deck_current_extinction_report_state(element)
        if extinction !== nothing
            if !extinction.closed
                current_values[index] = 0.0
            elseif extinction.current_initialized
                current_values[index] = extinction.previous_current
            else
                grounded_current = _deck_closed_switch_grounded_path_current(
                    context,
                    index,
                    voltage,
                )
                grounded_current === nothing ||
                    (current_values[index] = grounded_current)
            end
            continue
        end
        closed_flags[index] == 1 || continue
        grounded_current = _deck_closed_switch_grounded_path_current(
            context,
            index,
            voltage,
        )
        grounded_current === nothing ||
            (current_values[index] = grounded_current)
    end
    return current_values
end

function _deck_time_switch_report_current_values(
    boundary_run,
    voltage::AbstractVector{<:Real},
    time_s::Real,
)
    voltage_values = _deck_float64_voltage_values(voltage)
    closed_flags = Int[
        _deck_time_switch_closed_at(
            boundary_run.deck_time_switch_initially_closed_flags[index],
            boundary_run.deck_time_switch_close_time_s_values[index],
            boundary_run.deck_time_switch_open_time_s_values[index],
            Float64(time_s),
        ) ? 1 : 0
        for index in 1:boundary_run.deck_time_switch_count
    ]
    for index in eachindex(closed_flags)
        element = _deck_runtime_switch_element(boundary_run, index)
        extinction = _deck_current_extinction_report_state(element)
        extinction === nothing && continue
        closed_flags[index] = extinction.closed ? 1 : 0
    end
    conductance_values =
        _deck_time_switch_final_conductance_values(boundary_run, closed_flags)
    voltage_differences = Float64[
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
    current_values = _deck_time_switch_final_current_values(
        voltage_differences,
        conductance_values,
    )
    for index in 1:boundary_run.deck_time_switch_count
        element = _deck_runtime_switch_element(boundary_run, index)
        extinction = _deck_current_extinction_report_state(element)
        if extinction !== nothing
            if !extinction.closed
                current_values[index] = 0.0
            elseif extinction.current_initialized
                current_values[index] = extinction.previous_current
            else
                grounded_current = _deck_closed_switch_grounded_path_current(
                    boundary_run,
                    index,
                    voltage_values,
                )
                grounded_current === nothing ||
                    (current_values[index] = grounded_current)
            end
            continue
        end
        closed_flags[index] == 1 || continue
        grounded_current =
            _deck_closed_switch_grounded_path_current(boundary_run, index, voltage_values)
        grounded_current === nothing && continue
        current_values[index] = grounded_current
    end
    return current_values
end

function _deck_time_switch_final_power_values(final_voltage_values::Vector{Float64},
                                              final_current_values::Vector{Float64})
    return Float64[
        final_voltage_values[index] * final_current_values[index]
        for index in 1:length(final_voltage_values)
    ]
end

function _deck_time_switch_report_power_values(
    boundary_run,
    voltage::AbstractVector{<:Real},
    current_values::AbstractVector{<:Real},
)
    return Float64[
        _deck_node_voltage(voltage, boundary_run.deck_time_switch_from_node_indices[index]) *
        Float64(current_values[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_report_power_values!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
    current_values::AbstractVector{Float64},
)
    power_values = context.switch_power_step_values
    for index in 1:context.deck_time_switch_count
        power_values[index] =
            _deck_node_voltage(voltage, context.deck_time_switch_from_node_indices[index]) *
            current_values[index]
    end
    return power_values
end

_switch_signal_output_code(output_code::Integer) =
    output_code > 3 ? 3 : output_code

_switch_current_report_selected(output_code::Integer) =
    _switch_signal_output_code(output_code) in (1, 3)

_switch_voltage_report_selected(output_code::Integer) =
    _switch_signal_output_code(output_code) >= 2

function _deck_switch_voltage_report_indices(
    output_codes::AbstractVector{<:Integer},
    switch_count::Int,
)
    row_count = min(length(output_codes), switch_count)
    return Int[
        index for index in 1:row_count
        if _switch_voltage_report_selected(output_codes[index])
    ]
end

function _deck_switch_current_report_indices(
    output_codes::AbstractVector{<:Integer},
    switch_count::Int,
)
    row_count = min(length(output_codes), switch_count)
    return Int[
        index for index in 1:row_count
        if _switch_current_report_selected(output_codes[index])
    ]
end

function _deck_switch_current_report_indices(boundary_run)
    return _deck_switch_current_report_indices(
        boundary_run.deck_over5_switch_output_codes,
        boundary_run.deck_time_switch_count,
    )
end

function _deck_switch_voltage_report_indices(boundary_run)
    return _deck_switch_voltage_report_indices(
        boundary_run.deck_over5_switch_output_codes,
        boundary_run.deck_time_switch_count,
    )
end

_selected_symbol_values(values::AbstractVector{Symbol}, indices::AbstractVector{Int}) =
    Symbol[values[index] for index in indices]

_selected_int_values(values::AbstractVector{<:Integer}, indices::AbstractVector{Int}) =
    Int[values[index] for index in indices]

_selected_float_values(values::AbstractVector{<:Real}, indices::AbstractVector{Int}) =
    Float64[values[index] for index in indices]

function _deck_switch_current_report_names(boundary_run, indices::AbstractVector{Int})
    return _deck_switch_current_report_names(boundary_run.deck_over5_switch_names, indices)
end

function _deck_switch_voltage_report_names(boundary_run, indices::AbstractVector{Int})
    return _deck_switch_voltage_report_names(boundary_run.deck_over5_switch_names, indices)
end

function _deck_switch_voltage_report_names(
    switch_names::AbstractVector{Symbol},
    indices::AbstractVector{Int},
)
    return Symbol[
        Symbol("switch_voltage_", String(switch_names[index]))
        for index in indices
    ]
end

function _deck_switch_current_report_names(
    switch_names::AbstractVector{Symbol},
    indices::AbstractVector{Int},
)
    return Symbol[
        Symbol("switch_current_", String(switch_names[index]))
        for index in indices
    ]
end

function _deck_electrical_report_channel_values(voltage_output_values::AbstractVector{<:Real},
                                                branch_voltage_output_values::AbstractVector{<:Real},
                                                switch_voltage_report_values::AbstractVector{<:Real},
                                                branch_current_output_values::AbstractVector{<:Real},
                                                switch_current_report_values::AbstractVector{<:Real})
    return Float64[
        Float64.(voltage_output_values)...,
        Float64.(branch_voltage_output_values)...,
        Float64.(switch_voltage_report_values)...,
        Float64.(switch_current_report_values)...,
        Float64.(branch_current_output_values)...,
    ]
end

function _over16_post_extrema_results(boundary_run)
    return [
        update.post_extrema_result
        for update in _over16_boundary_pass_updates(boundary_run.over16_updates)
        if update.post_extrema_result !== nothing
    ]
end

function _over16_post_extrema_time_after_values(boundary_run)
    return Float64.([result.time_after for result in _over16_post_extrema_results(boundary_run)])
end

function _over16_post_extrema_step_after_values(boundary_run)
    return Int[result.istep_after for result in _over16_post_extrema_results(boundary_run)]
end

function _over16_post_extrema_final_time_s(boundary_run)
    values = _over16_post_extrema_time_after_values(boundary_run)
    return isempty(values) ? 0.0 : values[end]
end

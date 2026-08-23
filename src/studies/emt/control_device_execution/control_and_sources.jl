
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

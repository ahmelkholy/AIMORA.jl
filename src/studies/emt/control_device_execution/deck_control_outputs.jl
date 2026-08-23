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

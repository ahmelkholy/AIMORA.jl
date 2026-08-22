function _portable_emt_control_assignment_record(
    assignment::AlgebraicControlAssignment,
)
    return PortableSnapshotRecord(
        "aimora.emt.control_assignment.v1",
        Pair{String,Any}[
            "gain" => assignment.gain,
            "input_name" => String(assignment.input_name),
            "order" => assignment.order,
            "output_name" => String(assignment.output_name),
        ],
    )
end

function _portable_emt_control_function_record(
    function_row::ControlTransferFunction,
    runtime,
)
    return PortableSnapshotRecord(
        "aimora.emt.control_transfer_function.v1",
        Pair{String,Any}[
            "denominator_coefficients" => function_row.denominator_coefficients,
            "feedback_coefficients" => runtime.feedback_coefficients,
            "feedforward_coefficients" => runtime.feedforward_coefficients,
            "gain" => function_row.gain,
            "input_names" => String.(getfield.(function_row.input_terms, :name)),
            "input_polarities" => getfield.(function_row.input_terms, :polarity),
            "lower_limit" => (
                function_row.lower_limit === missing ?
                nothing : function_row.lower_limit
            ),
            "lower_limit_signal" => (
                function_row.lower_limit_signal === missing ?
                nothing : String(function_row.lower_limit_signal)
            ),
            "numerator_coefficients" => function_row.numerator_coefficients,
            "order" => function_row.order,
            "output_name" => String(function_row.output_name),
            "upper_limit" => (
                function_row.upper_limit === missing ?
                nothing : function_row.upper_limit
            ),
            "upper_limit_signal" => (
                function_row.upper_limit_signal === missing ?
                nothing : String(function_row.upper_limit_signal)
            ),
        ],
    )
end

function _portable_emt_control_device_record(device::ControlSystemDevice)
    return PortableSnapshotRecord(
        "aimora.emt.control_device.v1",
        Pair{String,Any}[
            "device_type" => device.device_type,
            "input_name" => String(device.input_name),
            "name" => String(device.name),
            "parameter_values" => device.parameter_values,
            "tail_signal_names" => String.(device.tail_signal_names),
        ],
    )
end

function _portable_emt_control_expression_record(
    runtime::ControlExpressionRuntime,
)
    program = runtime.program
    instructions = PortableSnapshotRecord[
        PortableSnapshotRecord(
            "aimora.emt.control_expression_instruction.v1",
            Pair{String,Any}[
                "name" => String(instruction.name),
                "operation" => String(instruction.operation),
                "value" => instruction.value,
            ],
        ) for instruction in program.instructions
    ]
    return PortableSnapshotRecord(
        "aimora.emt.control_expression.v1",
        Pair{String,Any}[
            "input_names" => String.(program.input_names),
            "instructions" => instructions,
            "max_stack_depth" => program.max_stack_depth,
            "output_name" => String(program.output_name),
            "source_text" => program.source_text,
        ],
    )
end

function _portable_emt_control_supplemental_row_record(
    row::OVER16CSUPDeviceRow,
)
    terms = PortableSnapshotRecord[
        PortableSnapshotRecord(
            "aimora.emt.control_supplemental_input.v1",
            Pair{String,Any}[
                "scale" => term.scale,
                "xtcs_index" => term.xtcs_index,
            ],
        ) for term in row.input_terms
    ]
    return PortableSnapshotRecord(
        "aimora.emt.control_supplemental_row.v1",
        Pair{String,Any}[
            "control_index" => row.control_index,
            "device_type" => row.device_type,
            "input_terms" => terms,
            "parsup_index" => row.parsup_index,
            "supplemental_index" => row.supplemental_index,
            "table_end_index" => row.table_end_index,
            "table_start_index" => row.table_start_index,
        ],
    )
end

function _portable_emt_control_supplemental_identity(
    runtime::Union{Nothing,ControlSystemSupplementalDeviceRuntime},
)
    runtime === nothing && return nothing
    return PortableSnapshotRecord(
        "aimora.emt.control_supplemental_identity.v1",
        Pair{String,Any}[
            "device_output_names" => String.(runtime.device_output_names),
            "device_output_slots" => runtime.device_output_slots,
            "initialize_rms_from_input" => collect(runtime.initialize_rms_from_input),
            "initialize_transport_delay_from_input" =>
                collect(runtime.initialize_transport_delay_from_input),
            "next_indices" => runtime.next_indices,
            "ordinary_signal_names" => String.(runtime.ordinary_signal_names),
            "rows" => _portable_emt_control_supplemental_row_record.(runtime.rows),
        ],
    )
end

function _portable_emt_control_frequency_initialization_record(
    initialization::ControlSystemFrequencyInitialization,
)
    return PortableSnapshotRecord(
        "aimora.emt.control_frequency_initialization.v1",
        Pair{String,Any}[
            "device_output_names" => String.(initialization.device_output_names),
            "frequency_hz" => initialization.frequency_hz,
            "history_mutation_count" => initialization.history_mutation_count,
            "signal_names" => String.(initialization.signal_names),
            "signal_phasor_imaginary" => imag.(initialization.signal_phasors),
            "signal_phasor_real" => real.(initialization.signal_phasors),
            "source_names" => String.(initialization.source_names),
        ],
    )
end

function _portable_emt_control_system_identity_record(
    runtime::ControlSystemNetworkRuntime,
)
    state = runtime.state
    length(state.functions) == length(state.function_states) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "control-system function definitions and runtime states differ in count",
        )
    clamp_names = Any[
        name === missing ? nothing : String(name)
        for name in runtime.switch_clamp_signal_names
    ]
    observations = PortableSnapshotRecord[
        PortableSnapshotRecord(
            "aimora.emt.control_switch_observation.v1",
            Pair{String,Any}[
                "name" => String(observation.name),
                "ordinary_switch_index" => observation.ordinary_switch_index,
                "source_type" => observation.source_type,
                "switch_kind" => String(nameof(typeof(observation.switch))),
            ],
        ) for observation in runtime.switch_observations
    ]
    return PortableSnapshotRecord(
        "aimora.emt.control_system_identity.v1",
        Pair{String,Any}[
            "assignments" => _portable_emt_control_assignment_record.(state.assignments),
            "control_output_names" => String.(runtime.control_output_names),
            "devices" => _portable_emt_control_device_record.(state.devices),
            "expressions" => _portable_emt_control_expression_record.(runtime.control_expressions),
            "frequency_hz" => state.frequency_hz,
            "frequency_initializations" =>
                _portable_emt_control_frequency_initialization_record.(
                    runtime.frequency_initializations,
                ),
            "functions" => PortableSnapshotRecord[
                _portable_emt_control_function_record(
                    state.functions[index],
                    state.function_states[index],
                ) for index in eachindex(state.functions)
            ],
            "output_names" => String.(state.output_names),
            "signal_slot_names" => String.(runtime.signal_slot_names),
            "supplemental" =>
                _portable_emt_control_supplemental_identity(
                    runtime.supplemental_devices,
                ),
            "switch_clamp_signal_names" => clamp_names,
            "switch_observations" => observations,
            "switch_output_indices" => runtime.switch_output_indices,
            "switch_output_names" => String.(runtime.switch_output_names),
            "switch_signal_names" => String.(runtime.switch_signal_names),
            "time_step_s" => state.deltat_s,
        ],
    )
end

function _portable_emt_control_named_values(
    values::Dict{Symbol,Float64},
)
    names = sort!(collect(keys(values)); by = String)
    return names, Float64[values[name] for name in names]
end

function _portable_emt_control_nested_float_fields(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    rows::AbstractVector{<:AbstractVector{<:Real}},
)
    row_lengths, flattened_values =
        _portable_emt_nested_source_values(rows, identity)
    return PortableSnapshotStateField[
        _portable_emt_array_field(
            "$identity.row_lengths",
            owner,
            :checkpoint,
            "1",
            ["row"],
            row_lengths,
        ),
        _portable_emt_array_field(
            "$identity.values",
            owner,
            family,
            unit,
            ["value"],
            flattened_values,
        ),
    ]
end

function _portable_emt_control_time_sequence_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    values::AbstractVector{<:Real},
    axis::AbstractString,
)
    return _portable_emt_state_field(
        identity,
        owner,
        family,
        "s",
        _portable_emt_time_sequence(values, identity);
        axes = [axis],
    )
end

function _portable_emt_control_source_state_fields(
    runtime::ControlSystemNetworkRuntime,
)
    windowed = runtime.windowed_constant_sources
    sinusoidal = runtime.sinusoidal_sources
    waveforms = runtime.waveform_sources
    return PortableSnapshotStateField[
        _portable_emt_state_field("control.source.windowed.names", "emt.control.source", :discrete, "1", String.(getfield.(windowed, :name)); axes = ["source"]),
        _portable_emt_array_field("control.source.windowed.values", "emt.control.source", :discrete, "pu", ["source"], Float64[getfield(source, :value) for source in windowed]),
        _portable_emt_control_time_sequence_field("control.source.windowed.start_time", "emt.control.source", :scheduler, Float64[getfield(source, :start_time_s) for source in windowed], "source"),
        _portable_emt_control_time_sequence_field("control.source.windowed.stop_time", "emt.control.source", :scheduler, Float64[getfield(source, :stop_time_s) for source in windowed], "source"),
        _portable_emt_state_field("control.source.sinusoidal.names", "emt.control.source", :discrete, "1", String.(getfield.(sinusoidal, :name)); axes = ["source"]),
        _portable_emt_array_field("control.source.sinusoidal.amplitude", "emt.control.source", :discrete, "pu", ["source"], Float64[getfield(source, :amplitude) for source in sinusoidal]),
        _portable_emt_array_field("control.source.sinusoidal.frequency", "emt.control.source", :discrete, "Hz", ["source"], Float64[getfield(source, :frequency_hz) for source in sinusoidal]),
        _portable_emt_array_field("control.source.sinusoidal.phase", "emt.control.source", :discrete, "rad", ["source"], Float64[getfield(source, :phase_rad) for source in sinusoidal]),
        _portable_emt_control_time_sequence_field("control.source.sinusoidal.start_time", "emt.control.source", :scheduler, Float64[getfield(source, :start_time_s) for source in sinusoidal], "source"),
        _portable_emt_control_time_sequence_field("control.source.sinusoidal.stop_time", "emt.control.source", :scheduler, Float64[getfield(source, :stop_time_s) for source in sinusoidal], "source"),
        _portable_emt_state_field("control.source.waveform.names", "emt.control.source", :discrete, "1", String.(getfield.(waveforms, :name)); axes = ["source"]),
        _portable_emt_array_field("control.source.waveform.types", "emt.control.source", :discrete, "1", ["source"], Int[getfield(source, :waveform_type) for source in waveforms]),
        _portable_emt_array_field("control.source.waveform.amplitude", "emt.control.source", :discrete, "pu", ["source"], Float64[getfield(source, :amplitude) for source in waveforms]),
        _portable_emt_array_field("control.source.waveform.cycle", "emt.control.source", :discrete, "s", ["source"], Float64[getfield(source, :cycle_s) for source in waveforms]),
        _portable_emt_array_field("control.source.waveform.pulse_width", "emt.control.source", :discrete, "s", ["source"], Float64[getfield(source, :pulse_width_s) for source in waveforms]),
        _portable_emt_control_time_sequence_field("control.source.waveform.start_time", "emt.control.source", :scheduler, Float64[getfield(source, :start_time_s) for source in waveforms], "source"),
        _portable_emt_control_time_sequence_field("control.source.waveform.stop_time", "emt.control.source", :scheduler, Float64[getfield(source, :stop_time_s) for source in waveforms], "source"),
    ]
end

function _portable_emt_control_supplemental_state_fields(
    runtime::Union{Nothing,ControlSystemSupplementalDeviceRuntime},
)
    fields = PortableSnapshotStateField[
        _portable_emt_state_field(
            "control.supplemental.present",
            "emt.control.supplemental",
            :checkpoint,
            "1",
            runtime !== nothing,
        ),
    ]
    runtime === nothing && return fields
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("control.supplemental.xtcs_values", "emt.control.supplemental", :algebraic, "pu", ["slot"], runtime.state.xtcs_values),
        _portable_emt_array_field("control.supplemental.parsup_values", "emt.control.supplemental", :history, "pu", ["slot"], runtime.state.parsup_values),
        _portable_emt_array_field("control.supplemental.device_integer_values", "emt.control.supplemental", :discrete, "1", ["device"], runtime.state.device_integer_values),
        _portable_emt_array_field("control.supplemental.row_reference_indices", "emt.control.supplemental", :scheduler, "1", ["device"], Int[row.reference_index for row in runtime.rows]),
        _portable_emt_array_field("control.supplemental.transport_delay_pointers", "emt.control.supplemental", :scheduler, "1", ["device"], runtime.transport_delay_pointers),
        _portable_emt_state_field("control.supplemental.initialized", "emt.control.supplemental", :discrete, "1", runtime.initialized),
        _portable_emt_state_field("control.supplemental.executed_step_count", "emt.control.supplemental", :output, "1", runtime.executed_step_count),
    ])
    return fields
end

function _portable_emt_control_system_state_fields(
    runtime::Union{Nothing,ControlSystemNetworkRuntime},
)
    fields = PortableSnapshotStateField[
        _portable_emt_state_field(
            "control.present",
            "emt.control",
            :checkpoint,
            "1",
            runtime !== nothing,
        ),
    ]
    runtime === nothing && return fields
    state_names, state_values = _portable_emt_control_named_values(runtime.state.values)
    history_names, history_values =
        _portable_emt_control_named_values(runtime.state.device_input_history)
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field(
            "control.identity",
            "emt.control",
            :checkpoint,
            "1",
            _portable_emt_control_system_identity_record(runtime),
        ),
        _portable_emt_state_field("control.values.names", "emt.control", :checkpoint, "1", String.(state_names); axes = ["signal"]),
        _portable_emt_array_field("control.values.accepted", "emt.control", :algebraic, "pu", ["signal"], state_values),
        _portable_emt_state_field("control.device_input_history.names", "emt.control", :checkpoint, "1", String.(history_names); axes = ["device"]),
        _portable_emt_array_field("control.device_input_history.values", "emt.control", :history, "pu", ["device"], history_values),
        _portable_emt_state_field("control.network_voltage_source_names", "emt.control", :discrete, "1", String.(runtime.network_voltage_source_names); axes = ["source"]),
        _portable_emt_array_field("control.network_voltage_source_node_indices", "emt.control", :discrete, "1", ["source"], runtime.network_voltage_source_node_indices),
        _portable_emt_state_field("control.function.limiter_active", "emt.control.function", :discrete, "1", Bool[state.limiter_active for state in runtime.state.function_states]; axes = ["function"]),
        _portable_emt_array_field("control.function.crossed_limit_count", "emt.control.function", :output, "1", ["function"], Int[state.crossed_limit_count for state in runtime.state.function_states]),
        _portable_emt_state_field("control.halted", "emt.control", :discrete, "1", runtime.halted),
        _portable_emt_state_field("control.last_diagnostic", "emt.control", :output, "1", runtime.last_diagnostic),
        _portable_emt_state_field("control.executed_step_count", "emt.control", :output, "1", runtime.executed_step_count),
        _portable_emt_state_field("control.feedback_application_count", "emt.control", :output, "1", runtime.feedback_application_count),
    ])
    append!(fields, _portable_emt_control_nested_float_fields(
        "control.function.history",
        "emt.control.function",
        :history,
        "pu",
        Vector{Float64}[
            state.history_terms for state in runtime.state.function_states
        ],
    ))
    append!(fields, _portable_emt_control_nested_float_fields(
        "control.expression.stack",
        "emt.control.expression",
        :algebraic,
        "pu",
        Vector{Float64}[
            expression.stack for expression in runtime.control_expressions
        ],
    ))
    append!(fields, _portable_emt_control_source_state_fields(runtime))
    append!(fields, _portable_emt_control_supplemental_state_fields(
        runtime.supplemental_devices,
    ))
    return fields
end

function _portable_emt_element_state_fields(
    element::TACSControlledSwitch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(
        index,
        name,
        "tacs_controlled_switch",
    )
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.node_a", "emt.controlled_switch", :checkpoint, "1", element.a),
        _portable_emt_state_field("$prefix.node_b", "emt.controlled_switch", :checkpoint, "1", element.b),
        _portable_emt_state_field("$prefix.clamp_present", "emt.controlled_switch", :checkpoint, "1", element.clamp_control !== nothing),
        _portable_emt_state_field("$prefix.threshold", "emt.controlled_switch", :checkpoint, "pu", element.threshold),
        _portable_emt_state_field("$prefix.on_conductance", "emt.controlled_switch", :checkpoint, "S", element.on_conductance),
        _portable_emt_state_field("$prefix.off_conductance", "emt.controlled_switch", :checkpoint, "S", element.off_conductance),
        _portable_emt_state_field("$prefix.unidirectional_latching", "emt.controlled_switch", :checkpoint, "1", element.unidirectional_latching),
        _portable_emt_state_field("$prefix.bidirectional_latching", "emt.controlled_switch", :checkpoint, "1", element.bidirectional_latching),
        _portable_emt_state_field("$prefix.ignition_voltage", "emt.controlled_switch", :checkpoint, "V", element.ignition_voltage),
        _portable_emt_state_field("$prefix.holding_current", "emt.controlled_switch", :checkpoint, "A", element.holding_current),
        _portable_emt_state_field("$prefix.deionization_time", "emt.controlled_switch", :checkpoint, "s", element.deionization_time_s),
        _portable_emt_state_field("$prefix.control", "emt.controlled_switch", :algebraic, "pu", element.control[]),
        _portable_emt_state_field("$prefix.closed", "emt.controlled_switch", :discrete, "1", element.closed),
        _portable_emt_state_field("$prefix.last_control", "emt.controlled_switch", :history, "pu", element.last_control),
        _portable_emt_state_field("$prefix.last_conductance", "emt.controlled_switch", :history, "S", element.last_conductance),
        _portable_emt_time_field("$prefix.open_duration", "emt.controlled_switch", :history, element.open_duration_s),
        _portable_emt_state_field("$prefix.last_voltage", "emt.controlled_switch", :algebraic, "V", element.last_voltage),
        _portable_emt_state_field("$prefix.last_current", "emt.controlled_switch", :algebraic, "A", element.last_current),
        _portable_emt_state_field("$prefix.delayed_arc_present", "emt.controlled_switch", :checkpoint, "1", element.delayed_arc !== nothing),
    ])
    element.clamp_control === nothing || push!(
        fields,
        _portable_emt_state_field(
            "$prefix.clamp_control",
            "emt.controlled_switch",
            :algebraic,
            "pu",
            element.clamp_control[],
        ),
    )
    arc = element.delayed_arc
    arc === nothing && return fields
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.arc.current_coefficient", "emt.controlled_switch.arc", :checkpoint, "1", arc.current_coefficient),
        _portable_emt_state_field("$prefix.arc.current_exponent", "emt.controlled_switch.arc", :checkpoint, "1", arc.current_exponent),
        _portable_emt_state_field("$prefix.arc.time_scale", "emt.controlled_switch.arc", :checkpoint, "s", arc.time_scale_s),
        _portable_emt_state_field("$prefix.arc.cutoff_current", "emt.controlled_switch.arc", :checkpoint, "A", arc.cutoff_current_a),
        _portable_emt_state_field("$prefix.arc.previous_current", "emt.controlled_switch.arc", :history, "A", arc.previous_current_a),
        _portable_emt_state_field("$prefix.arc.tail_amplitude", "emt.controlled_switch.arc", :history, "A", arc.tail_amplitude_a),
        _portable_emt_state_field("$prefix.arc.tail_time_constant", "emt.controlled_switch.arc", :history, "s", arc.tail_time_constant_s),
        _portable_emt_time_field("$prefix.arc.scheduled_open_time", "emt.controlled_switch.arc", :scheduler, arc.scheduled_open_time_s),
        _portable_emt_state_field("$prefix.arc.opening_requested", "emt.controlled_switch.arc", :discrete, "1", arc.opening_requested),
        _portable_emt_state_field("$prefix.arc.transition_deferred", "emt.controlled_switch.arc", :discrete, "1", arc.transition_deferred),
        _portable_emt_state_field("$prefix.arc.tail_active", "emt.controlled_switch.arc", :discrete, "1", arc.tail_active),
        _portable_emt_state_field("$prefix.arc.tail_current", "emt.controlled_switch.arc", :algebraic, "A", arc.tail_current_a),
        _portable_emt_state_field("$prefix.arc.decay_factor", "emt.controlled_switch.arc", :algebraic, "1", arc.decay_factor),
        _portable_emt_state_field("$prefix.arc.transition_count", "emt.controlled_switch.arc", :output, "1", arc.transition_count),
        _portable_emt_state_field("$prefix.arc.absolute_tail_energy", "emt.controlled_switch.arc", :output, "J", arc.absolute_tail_energy_j),
    ])
    return fields
end

function _portable_emt_control_decode_time_value(
    value,
    identity::AbstractString,
)
    value isa Float64 && return value
    value == _PORTABLE_EMT_POSITIVE_INFINITY && return Inf
    value == _PORTABLE_EMT_NEGATIVE_INFINITY && return -Inf
    _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT time sequence $identity has an invalid value",
    )
end

function _portable_emt_control_time_sequence_values(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    axis::AbstractString,
)
    field = _portable_emt_inventory_field(fields, identity, family, "s")
    field.axes == [String(axis)] || _portable_emt_fail(
        :state_axis_mismatch,
        "portable EMT time sequence $identity has the wrong axis",
    )
    field.value isa AbstractVector || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT time sequence $identity is not a sequence",
    )
    return Float64[
        _portable_emt_control_decode_time_value(
            value,
            "$identity[$index]",
        ) for (index, value) in enumerate(field.value)
    ]
end

function _portable_emt_control_string_sequence(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    axis::AbstractString,
)
    field = _portable_emt_inventory_field(fields, identity, family, "1")
    field.axes == [String(axis)] || _portable_emt_fail(
        :state_axis_mismatch,
        "portable EMT string sequence $identity has the wrong axis",
    )
    value = field.value
    value isa AbstractVector && all(item -> item isa AbstractString, value) ||
        _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT state field $identity is not a string sequence",
        )
    return String[item for item in value]
end

function _portable_emt_restore_control_nested_float_rows!(
    destinations::AbstractVector{<:AbstractVector{Float64}},
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
)
    row_lengths = _portable_emt_array_values(
        fields,
        "$identity.row_lengths",
        :checkpoint,
        "1",
        ["row"],
        Int64,
    )
    values = _portable_emt_array_values(
        fields,
        "$identity.values",
        family,
        unit,
        ["value"],
        Float64,
    )
    length(row_lengths) == length(destinations) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT nested state $identity has the wrong row count",
    )
    all(row_lengths[index] == length(destinations[index]) for index in eachindex(destinations)) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT nested state $identity has the wrong row lengths",
        )
    sum(row_lengths; init = Int64(0)) == length(values) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT nested state $identity has an inconsistent value count",
        )
    offset = 0
    for index in eachindex(destinations)
        count = _portable_emt_checked_integer(
            row_lengths[index],
            "$identity.row_lengths",
        )
        copyto!(destinations[index], @view(values[(offset + 1):(offset + count)]))
        offset += count
    end
    return destinations
end

function _portable_emt_validate_control_system_identity(
    runtime::ControlSystemNetworkRuntime,
    fields::Dict{String,PortableSnapshotStateField},
)
    reconstructed = PortableSnapshotStateInventory(PortableSnapshotStateField[
        _portable_emt_state_field(
            "control.identity",
            "emt.control",
            :checkpoint,
            "1",
            _portable_emt_control_system_identity_record(runtime),
        ),
    ])
    captured = PortableSnapshotStateInventory(PortableSnapshotStateField[
        _portable_emt_inventory_field(
            fields,
            "control.identity",
            :checkpoint,
            "1",
        ),
    ])
    captured.signature_sha256 == reconstructed.signature_sha256 ||
        _portable_emt_fail(
            :settings_mismatch,
            "portable EMT control-system equations or owner mappings changed in the receiving study",
        )
    return runtime
end

function _restore_portable_emt_control_sources!(
    runtime::ControlSystemNetworkRuntime,
    fields::Dict{String,PortableSnapshotStateField},
)
    windowed_names = _portable_emt_control_string_sequence(
        fields,
        "control.source.windowed.names",
        :discrete,
        "source",
    )
    windowed_values = _portable_emt_array_values(fields, "control.source.windowed.values", :discrete, "pu", ["source"], Float64)
    windowed_start = _portable_emt_control_time_sequence_values(fields, "control.source.windowed.start_time", :scheduler, "source")
    windowed_stop = _portable_emt_control_time_sequence_values(fields, "control.source.windowed.stop_time", :scheduler, "source")
    windowed_count = length(windowed_names)
    all(length(values) == windowed_count for values in (
        windowed_values,
        windowed_start,
        windowed_stop,
    )) || _portable_emt_fail(:state_shape_mismatch, "portable EMT windowed control sources have inconsistent counts")
    empty!(runtime.windowed_constant_sources)
    for index in 1:windowed_count
        push!(runtime.windowed_constant_sources, WindowedConstantControlSignal(
            Symbol(windowed_names[index]),
            windowed_values[index],
            windowed_start[index],
            windowed_stop[index],
        ))
    end

    sinusoidal_names = _portable_emt_control_string_sequence(fields, "control.source.sinusoidal.names", :discrete, "source")
    sinusoidal_amplitudes = _portable_emt_array_values(fields, "control.source.sinusoidal.amplitude", :discrete, "pu", ["source"], Float64)
    sinusoidal_frequencies = _portable_emt_array_values(fields, "control.source.sinusoidal.frequency", :discrete, "Hz", ["source"], Float64)
    sinusoidal_phases = _portable_emt_array_values(fields, "control.source.sinusoidal.phase", :discrete, "rad", ["source"], Float64)
    sinusoidal_start = _portable_emt_control_time_sequence_values(fields, "control.source.sinusoidal.start_time", :scheduler, "source")
    sinusoidal_stop = _portable_emt_control_time_sequence_values(fields, "control.source.sinusoidal.stop_time", :scheduler, "source")
    sinusoidal_count = length(sinusoidal_names)
    all(length(values) == sinusoidal_count for values in (
        sinusoidal_amplitudes,
        sinusoidal_frequencies,
        sinusoidal_phases,
        sinusoidal_start,
        sinusoidal_stop,
    )) || _portable_emt_fail(:state_shape_mismatch, "portable EMT sinusoidal control sources have inconsistent counts")
    empty!(runtime.sinusoidal_sources)
    for index in 1:sinusoidal_count
        push!(runtime.sinusoidal_sources, SinusoidalControlSignal(
            Symbol(sinusoidal_names[index]),
            sinusoidal_amplitudes[index],
            sinusoidal_frequencies[index],
            sinusoidal_phases[index];
            start_time_s = sinusoidal_start[index],
            stop_time_s = sinusoidal_stop[index],
        ))
    end

    waveform_names = _portable_emt_control_string_sequence(fields, "control.source.waveform.names", :discrete, "source")
    waveform_types = _portable_emt_array_values(fields, "control.source.waveform.types", :discrete, "1", ["source"], Int64)
    waveform_amplitudes = _portable_emt_array_values(fields, "control.source.waveform.amplitude", :discrete, "pu", ["source"], Float64)
    waveform_cycles = _portable_emt_array_values(fields, "control.source.waveform.cycle", :discrete, "s", ["source"], Float64)
    waveform_widths = _portable_emt_array_values(fields, "control.source.waveform.pulse_width", :discrete, "s", ["source"], Float64)
    waveform_start = _portable_emt_control_time_sequence_values(fields, "control.source.waveform.start_time", :scheduler, "source")
    waveform_stop = _portable_emt_control_time_sequence_values(fields, "control.source.waveform.stop_time", :scheduler, "source")
    waveform_count = length(waveform_names)
    all(length(values) == waveform_count for values in (
        waveform_types,
        waveform_amplitudes,
        waveform_cycles,
        waveform_widths,
        waveform_start,
        waveform_stop,
    )) || _portable_emt_fail(:state_shape_mismatch, "portable EMT waveform control sources have inconsistent counts")
    empty!(runtime.waveform_sources)
    for index in 1:waveform_count
        push!(runtime.waveform_sources, ControlWaveformSignal(
            Symbol(waveform_names[index]),
            _portable_emt_checked_integer(waveform_types[index], "control.source.waveform.types"),
            waveform_amplitudes[index],
            waveform_cycles[index],
            waveform_widths[index],
            waveform_start[index],
            waveform_stop[index],
        ))
    end
    return runtime
end

function _restore_portable_emt_control_supplemental_state!(
    runtime::ControlSystemNetworkRuntime,
    fields::Dict{String,PortableSnapshotStateField},
)
    present = _portable_emt_scalar(
        fields,
        "control.supplemental.present",
        :checkpoint,
        "1",
        Bool,
    )
    supplemental = runtime.supplemental_devices
    (supplemental !== nothing) == present || _portable_emt_fail(
        :settings_mismatch,
        "portable EMT supplemental-device ownership changed in the receiving study",
    )
    supplemental === nothing && return runtime
    _restore_portable_emt_array!(supplemental.state.xtcs_values, fields, "control.supplemental.xtcs_values", :algebraic, "pu", ["slot"])
    _restore_portable_emt_array!(supplemental.state.parsup_values, fields, "control.supplemental.parsup_values", :history, "pu", ["slot"])
    _restore_portable_emt_array!(supplemental.state.device_integer_values, fields, "control.supplemental.device_integer_values", :discrete, "1", ["device"])
    reference_indices = Int[0 for _ in supplemental.rows]
    _restore_portable_emt_array!(reference_indices, fields, "control.supplemental.row_reference_indices", :scheduler, "1", ["device"])
    _restore_portable_emt_array!(supplemental.transport_delay_pointers, fields, "control.supplemental.transport_delay_pointers", :scheduler, "1", ["device"])
    for index in eachindex(supplemental.rows)
        row = supplemental.rows[index]
        supplemental.rows[index] = over16_csup_device_row(
            row.supplemental_index,
            row.device_type,
            row.parsup_index;
            input_terms = row.input_terms,
            control_index = row.control_index,
            reference_index = reference_indices[index],
            table_start_index = row.table_start_index,
            table_end_index = row.table_end_index,
        )
    end
    supplemental.initialized = _portable_emt_scalar(fields, "control.supplemental.initialized", :discrete, "1", Bool)
    supplemental.executed_step_count = _portable_emt_nonnegative_integer(fields, "control.supplemental.executed_step_count", :output)
    return runtime
end

function _restore_portable_emt_control_system_state!(
    context::EMTStepContext,
    fields::Dict{String,PortableSnapshotStateField},
)
    present = _portable_emt_scalar(
        fields,
        "control.present",
        :checkpoint,
        "1",
        Bool,
    )
    runtime = context.control_system_runtime
    (runtime !== nothing) == present || _portable_emt_fail(
        :settings_mismatch,
        "portable EMT control-system ownership changed in the receiving study",
    )
    runtime === nothing && return context
    _portable_emt_validate_control_system_identity(runtime, fields)

    value_names = _portable_emt_control_string_sequence(
        fields,
        "control.values.names",
        :checkpoint,
        "signal",
    )
    reconstructed_names = sort!(String.(collect(keys(runtime.state.values))))
    value_names == reconstructed_names || _portable_emt_fail(
        :settings_mismatch,
        "portable EMT control signal ownership changed in the receiving study",
    )
    value_data = _portable_emt_array_values(
        fields,
        "control.values.accepted",
        :algebraic,
        "pu",
        ["signal"],
        Float64,
    )
    length(value_data) == length(value_names) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT control values have the wrong signal count",
    )
    for index in eachindex(value_names)
        runtime.state.values[Symbol(value_names[index])] = value_data[index]
    end

    history_names = _portable_emt_control_string_sequence(
        fields,
        "control.device_input_history.names",
        :checkpoint,
        "device",
    )
    device_names = Set(getfield.(runtime.state.devices, :name))
    all(name -> Symbol(name) in device_names, history_names) ||
        _portable_emt_fail(
            :settings_mismatch,
            "portable EMT device-input history references an unknown device",
        )
    history_values = _portable_emt_array_values(
        fields,
        "control.device_input_history.values",
        :history,
        "pu",
        ["device"],
        Float64,
    )
    length(history_values) == length(history_names) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT device-input histories have inconsistent counts",
    )
    empty!(runtime.state.device_input_history)
    for index in eachindex(history_names)
        runtime.state.device_input_history[Symbol(history_names[index])] =
            history_values[index]
    end

    _portable_emt_restore_control_nested_float_rows!(
        Vector{Float64}[
            state.history_terms for state in runtime.state.function_states
        ],
        fields,
        "control.function.history",
        :history,
        "pu",
    )
    limiter_active = _portable_emt_boolean_sequence(
        fields,
        "control.function.limiter_active",
        :discrete;
        axes = ["function"],
    )
    length(limiter_active) == length(runtime.state.function_states) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT control limiter state has the wrong function count",
        )
    crossed_limit_count = Int[0 for _ in runtime.state.function_states]
    _restore_portable_emt_array!(crossed_limit_count, fields, "control.function.crossed_limit_count", :output, "1", ["function"])
    for index in eachindex(runtime.state.function_states)
        runtime.state.function_states[index].limiter_active = limiter_active[index]
        runtime.state.function_states[index].crossed_limit_count =
            crossed_limit_count[index]
    end
    _portable_emt_restore_control_nested_float_rows!(
        Vector{Float64}[
            expression.stack for expression in runtime.control_expressions
        ],
        fields,
        "control.expression.stack",
        :algebraic,
        "pu",
    )

    network_source_names = _portable_emt_control_string_sequence(
        fields,
        "control.network_voltage_source_names",
        :discrete,
        "source",
    )
    network_source_indices = _portable_emt_array_values(
        fields,
        "control.network_voltage_source_node_indices",
        :discrete,
        "1",
        ["source"],
        Int64,
    )
    length(network_source_names) == length(network_source_indices) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT network-control source mappings have inconsistent counts",
        )
    all(index -> 1 <= index <= context.system.node_count, network_source_indices) ||
        _portable_emt_fail(
            :state_value_mismatch,
            "portable EMT network-control source references an invalid node",
        )
    empty!(runtime.network_voltage_source_names)
    append!(runtime.network_voltage_source_names, Symbol.(network_source_names))
    empty!(runtime.network_voltage_source_node_indices)
    append!(runtime.network_voltage_source_node_indices, Int.(network_source_indices))

    _restore_portable_emt_control_sources!(runtime, fields)
    _restore_portable_emt_control_supplemental_state!(runtime, fields)
    runtime.halted = _portable_emt_scalar(fields, "control.halted", :discrete, "1", Bool)
    diagnostic = _portable_emt_inventory_field(
        fields,
        "control.last_diagnostic",
        :output,
        "1",
    ).value
    diagnostic === nothing || diagnostic isa AbstractString ||
        _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT control diagnostic is neither absent nor text",
        )
    runtime.last_diagnostic = diagnostic === nothing ? nothing : String(diagnostic)
    runtime.executed_step_count = _portable_emt_nonnegative_integer(
        fields,
        "control.executed_step_count",
        :output,
    )
    runtime.feedback_application_count = _portable_emt_nonnegative_integer(
        fields,
        "control.feedback_application_count",
        :output,
    )
    return context
end

function _portable_emt_validate_controlled_switch_setting(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
    expected,
)
    captured = _portable_emt_scalar(fields, identity, family, unit, typeof(expected))
    captured == expected || _portable_emt_fail(
        :model_mismatch,
        "portable EMT controlled-switch setting $identity changed in the receiving study",
    )
    return expected
end

function _restore_portable_emt_element_state!(
    element::TACSControlledSwitch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.node_a", :checkpoint, "1", element.a)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.node_b", :checkpoint, "1", element.b)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.clamp_present", :checkpoint, "1", element.clamp_control !== nothing)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.threshold", :checkpoint, "pu", element.threshold)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.on_conductance", :checkpoint, "S", element.on_conductance)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.off_conductance", :checkpoint, "S", element.off_conductance)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.unidirectional_latching", :checkpoint, "1", element.unidirectional_latching)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.bidirectional_latching", :checkpoint, "1", element.bidirectional_latching)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.ignition_voltage", :checkpoint, "V", element.ignition_voltage)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.holding_current", :checkpoint, "A", element.holding_current)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.deionization_time", :checkpoint, "s", element.deionization_time_s)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.delayed_arc_present", :checkpoint, "1", element.delayed_arc !== nothing)
    element.control[] = _portable_emt_scalar(fields, "$prefix.control", :algebraic, "pu", Float64)
    element.clamp_control === nothing || (
        element.clamp_control[] = _portable_emt_scalar(
            fields,
            "$prefix.clamp_control",
            :algebraic,
            "pu",
            Float64,
        )
    )
    element.closed = _portable_emt_scalar(fields, "$prefix.closed", :discrete, "1", Bool)
    element.last_control = _portable_emt_scalar(fields, "$prefix.last_control", :history, "pu", Float64)
    element.last_conductance = _portable_emt_scalar(fields, "$prefix.last_conductance", :history, "S", Float64)
    element.open_duration_s = _portable_emt_time_scalar(fields, "$prefix.open_duration", :history)
    element.last_voltage = _portable_emt_scalar(fields, "$prefix.last_voltage", :algebraic, "V", Float64)
    element.last_current = _portable_emt_scalar(fields, "$prefix.last_current", :algebraic, "A", Float64)
    arc = element.delayed_arc
    arc === nothing && return element
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.arc.current_coefficient", :checkpoint, "1", arc.current_coefficient)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.arc.current_exponent", :checkpoint, "1", arc.current_exponent)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.arc.time_scale", :checkpoint, "s", arc.time_scale_s)
    _portable_emt_validate_controlled_switch_setting(fields, "$prefix.arc.cutoff_current", :checkpoint, "A", arc.cutoff_current_a)
    arc.previous_current_a = _portable_emt_scalar(fields, "$prefix.arc.previous_current", :history, "A", Float64)
    arc.tail_amplitude_a = _portable_emt_scalar(fields, "$prefix.arc.tail_amplitude", :history, "A", Float64)
    arc.tail_time_constant_s = _portable_emt_scalar(fields, "$prefix.arc.tail_time_constant", :history, "s", Float64)
    arc.scheduled_open_time_s = _portable_emt_time_scalar(fields, "$prefix.arc.scheduled_open_time", :scheduler)
    arc.opening_requested = _portable_emt_scalar(fields, "$prefix.arc.opening_requested", :discrete, "1", Bool)
    arc.transition_deferred = _portable_emt_scalar(fields, "$prefix.arc.transition_deferred", :discrete, "1", Bool)
    arc.tail_active = _portable_emt_scalar(fields, "$prefix.arc.tail_active", :discrete, "1", Bool)
    arc.tail_current_a = _portable_emt_scalar(fields, "$prefix.arc.tail_current", :algebraic, "A", Float64)
    arc.decay_factor = _portable_emt_scalar(fields, "$prefix.arc.decay_factor", :algebraic, "1", Float64)
    arc.transition_count = _portable_emt_nonnegative_integer(fields, "$prefix.arc.transition_count", :output)
    arc.absolute_tail_energy_j = _portable_emt_scalar(fields, "$prefix.arc.absolute_tail_energy", :output, "J", Float64)
    return element
end

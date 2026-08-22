function _portable_emt_fail(code::Symbol, message::AbstractString)
    throw(PortableSnapshotFailure(code, String(message)))
end

function _portable_emt_state_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    value;
    axes::AbstractVector{<:AbstractString}=String[],
)
    return PortableSnapshotStateField(identity, owner, family, unit, axes, value)
end

const _PORTABLE_EMT_POSITIVE_INFINITY = "positive_infinity"
const _PORTABLE_EMT_NEGATIVE_INFINITY = "negative_infinity"

function _portable_emt_time_value(value::Real, identity::AbstractString)
    time = Float64(value)
    isnan(time) && _portable_emt_fail(
        :nonfinite_value,
        "portable EMT time field $identity is NaN",
    )
    isfinite(time) && return time
    time == Inf && return _PORTABLE_EMT_POSITIVE_INFINITY
    time == -Inf && return _PORTABLE_EMT_NEGATIVE_INFINITY
    _portable_emt_fail(
        :nonfinite_value,
        "portable EMT time field $identity has an unsupported value",
    )
end

function _portable_emt_time_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    value::Real,
)
    return _portable_emt_state_field(
        identity,
        owner,
        family,
        "s",
        _portable_emt_time_value(value, identity),
    )
end

function _portable_emt_time_sequence(
    values::AbstractVector{<:Real},
    identity::AbstractString,
)
    return Any[
        _portable_emt_time_value(value, "$identity[$index]")
        for (index, value) in enumerate(values)
    ]
end

function _portable_emt_float_array(
    values::AbstractArray{Float64},
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
)
    return portable_snapshot_array(values; unit, axes)
end

function _portable_emt_integer_array(
    values::AbstractArray{<:Integer},
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
)
    converted = try
        Int64.(values)
    catch error
        _portable_emt_fail(
            :integer_overflow,
            "EMT state array cannot be represented as portable Int64: $(sprint(showerror, error))",
        )
    end
    return portable_snapshot_array(converted; unit, axes)
end

function _portable_emt_array_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
    values::AbstractArray{Float64},
)
    return _portable_emt_state_field(
        identity,
        owner,
        family,
        unit,
        _portable_emt_float_array(values, unit, axes);
        axes,
    )
end

function _portable_emt_array_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
    values::AbstractArray{<:Integer},
)
    return _portable_emt_state_field(
        identity,
        owner,
        family,
        unit,
        _portable_emt_integer_array(values, unit, axes);
        axes,
    )
end

function _portable_emt_element_prefix(index::Integer)
    index > 0 || _portable_emt_fail(
        :invalid_state_owner,
        "EMT element state index must be positive",
    )
    return "model.i" * lpad(string(index), 8, '0')
end

function _portable_emt_element_identity_field(
    index::Integer,
    name::Symbol,
    kind::AbstractString,
)
    prefix = _portable_emt_element_prefix(index)
    return PortableSnapshotStateField[
        _portable_emt_state_field(
            "$prefix.owner_name",
            "emt.model",
            :checkpoint,
            "1",
            String(name),
        ),
        _portable_emt_state_field(
            "$prefix.model_kind",
            "emt.model",
            :checkpoint,
            "1",
            String(kind),
        ),
    ]
end

function _portable_emt_element_state_fields(
    element::BergeronLine,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "bergeron_line")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field(
            "$prefix.from_wave_history",
            "emt.line_history",
            :delayed,
            "pu",
            ["delay"],
            element.from_wave_history,
        ),
        _portable_emt_array_field(
            "$prefix.to_wave_history",
            "emt.line_history",
            :delayed,
            "pu",
            ["delay"],
            element.to_wave_history,
        ),
        _portable_emt_state_field("$prefix.write_index", "emt.line_history", :delayed, "1", element.write_index),
        _portable_emt_state_field("$prefix.history_from", "emt.line_history", :history, "pu", element.h_from),
        _portable_emt_state_field("$prefix.history_to", "emt.line_history", :history, "pu", element.h_to),
        _portable_emt_state_field("$prefix.terminal_voltage_from", "emt.line_state", :algebraic, "pu", element.v_from),
        _portable_emt_state_field("$prefix.terminal_voltage_to", "emt.line_state", :algebraic, "pu", element.v_to),
        _portable_emt_state_field("$prefix.terminal_current_from", "emt.line_state", :algebraic, "pu", element.i_from),
        _portable_emt_state_field("$prefix.terminal_current_to", "emt.line_state", :algebraic, "pu", element.i_to),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::SeriesRLBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "series_rl")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.previous_current", "emt.branch_history", :history, "A", element.i_prev),
        _portable_emt_state_field("$prefix.previous_voltage", "emt.branch_history", :history, "V", element.v_prev),
        _portable_emt_state_field("$prefix.last_current", "emt.branch_state", :algebraic, "A", element.i_last),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::SeriesRLCBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "series_rlc")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.resistance", "emt.branch_parameter", :discrete, "ohm", element.r),
        _portable_emt_state_field("$prefix.inductance", "emt.branch_parameter", :discrete, "H", element.l),
        _portable_emt_state_field("$prefix.capacitance", "emt.branch_parameter", :discrete, "F", element.c),
        _portable_emt_state_field("$prefix.previous_current", "emt.branch_history", :history, "A", element.i_prev),
        _portable_emt_state_field("$prefix.previous_inductor_voltage", "emt.branch_history", :history, "V", element.inductor_voltage_prev),
        _portable_emt_state_field("$prefix.previous_capacitor_voltage", "emt.branch_history", :history, "V", element.capacitor_voltage_prev),
        _portable_emt_state_field("$prefix.previous_terminal_voltage", "emt.branch_history", :history, "V", element.v_prev),
        _portable_emt_state_field("$prefix.last_current", "emt.branch_state", :algebraic, "A", element.i_last),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::CapacitorBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "capacitor")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.capacitance", "emt.branch_parameter", :checkpoint, "F", element.c),
        _portable_emt_state_field("$prefix.previous_current", "emt.branch_history", :history, "A", element.i_prev),
        _portable_emt_state_field("$prefix.previous_voltage", "emt.branch_history", :history, "V", element.v_prev),
        _portable_emt_state_field("$prefix.last_current", "emt.branch_state", :algebraic, "A", element.i_last),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::CoupledInductiveBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "coupled_inductive")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.previous_current", "emt.branch_history", :history, "A", ["port"], element.previous_current),
        _portable_emt_array_field("$prefix.previous_voltage", "emt.branch_history", :history, "V", ["port"], element.previous_voltage),
        _portable_emt_array_field("$prefix.last_current", "emt.branch_state", :algebraic, "A", ["port"], element.last_current),
        _portable_emt_array_field("$prefix.conductance_workspace", "emt.branch_state", :algebraic, "S", ["row", "column"], element.conductance_workspace),
        _portable_emt_array_field("$prefix.history_current_workspace", "emt.branch_state", :algebraic, "A", ["port"], element.history_current_workspace),
        _portable_emt_array_field("$prefix.port_voltage_workspace", "emt.branch_state", :algebraic, "V", ["port"], element.port_voltage_workspace),
        _portable_emt_array_field("$prefix.current_workspace", "emt.branch_state", :algebraic, "A", ["port"], element.current_workspace),
    ])
    return fields
end

function _portable_emt_cached_timestep_value(value::Float64)
    isnan(value) && return "not_initialized"
    isfinite(value) && value > 0.0 || _portable_emt_fail(
        :nonfinite_value,
        "coupled series R-L cached timestep is invalid",
    )
    return value
end

function _portable_emt_element_state_fields(
    element::CoupledSeriesRLBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "coupled_series_rl")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.previous_current", "emt.branch_history", :history, "A", ["port"], element.previous_current),
        _portable_emt_array_field("$prefix.previous_voltage", "emt.branch_history", :history, "V", ["port"], element.previous_voltage),
        _portable_emt_array_field("$prefix.last_current", "emt.branch_state", :algebraic, "A", ["port"], element.last_current),
        _portable_emt_array_field("$prefix.conductance_workspace", "emt.branch_state", :algebraic, "S", ["row", "column"], element.conductance_workspace),
        _portable_emt_array_field("$prefix.history_current_workspace", "emt.branch_state", :algebraic, "A", ["port"], element.history_current_workspace),
        _portable_emt_array_field("$prefix.port_voltage_workspace", "emt.branch_state", :algebraic, "V", ["port"], element.port_voltage_workspace),
        _portable_emt_array_field("$prefix.current_workspace", "emt.branch_state", :algebraic, "A", ["port"], element.current_workspace),
        _portable_emt_state_field("$prefix.cached_timestep", "emt.branch_state", :reconstruction, "s", _portable_emt_cached_timestep_value(element.cached_dt_s)),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::SaturableInductorBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "saturable_inductor")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.previous_current", "emt.nonlinear_history", :history, "A", element.i_prev),
        _portable_emt_state_field("$prefix.previous_voltage", "emt.nonlinear_history", :history, "V", element.v_prev),
        _portable_emt_state_field("$prefix.last_current", "emt.nonlinear_state", :algebraic, "A", element.i_last),
        _portable_emt_state_field("$prefix.last_inductance", "emt.nonlinear_state", :discrete, "H", element.last_inductance),
        _portable_emt_state_field("$prefix.saturated", "emt.nonlinear_state", :discrete, "1", element.saturated),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element::SaturatedTransformerNonlinearSlopeBranch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(
        index,
        name,
        "saturated_transformer_nonlinear_slope",
    )
    push!(
        fields,
        _portable_emt_state_field(
            "$prefix.conductance",
            "emt.nonlinear_state",
            :algebraic,
            "S",
            element.conductance,
        ),
    )
    return fields
end

function _portable_emt_alteration_state_fields(context::EMTStepContext)
    isempty(context.series_rlc_alterations) &&
        isempty(context.series_rlc_alteration_records) && return PortableSnapshotStateField[]
    schedule = context.series_rlc_alterations
    records = context.series_rlc_alteration_records
    return PortableSnapshotStateField[
        _portable_emt_state_field("alteration.schedule_branch_names", "emt.alteration", :checkpoint, "1", String.(getfield.(schedule, :branch_name))),
        _portable_emt_array_field("alteration.schedule_activation_time", "emt.alteration", :checkpoint, "s", ["alteration"], getfield.(schedule, :activation_time_s)),
        _portable_emt_array_field("alteration.schedule_resistance", "emt.alteration", :checkpoint, "ohm", ["alteration"], getfield.(schedule, :resistance_ohm)),
        _portable_emt_array_field("alteration.schedule_inductance", "emt.alteration", :checkpoint, "H", ["alteration"], getfield.(schedule, :inductance_h)),
        _portable_emt_array_field("alteration.schedule_capacitance", "emt.alteration", :checkpoint, "F", ["alteration"], getfield.(schedule, :capacitance_f)),
        _portable_emt_array_field("alteration.schedule_branch_indices", "emt.alteration", :checkpoint, "1", ["alteration"], context.series_rlc_alteration_branch_indices),
        _portable_emt_state_field("alteration.record_branch_names", "emt.alteration", :output, "1", String.(getfield.(records, :branch_name))),
        _portable_emt_array_field("alteration.record_branch_indices", "emt.alteration", :output, "1", ["record"], getfield.(records, :branch_index)),
        _portable_emt_array_field("alteration.record_requested_time", "emt.alteration", :output, "s", ["record"], getfield.(records, :requested_time_s)),
        _portable_emt_array_field("alteration.record_applied_time", "emt.alteration", :output, "s", ["record"], getfield.(records, :applied_time_s)),
        _portable_emt_array_field("alteration.record_applied_step", "emt.alteration", :output, "1", ["record"], getfield.(records, :applied_step_index)),
        _portable_emt_array_field("alteration.record_previous_resistance", "emt.alteration", :output, "ohm", ["record"], getfield.(records, :previous_resistance_ohm)),
        _portable_emt_array_field("alteration.record_previous_inductance", "emt.alteration", :output, "H", ["record"], getfield.(records, :previous_inductance_h)),
        _portable_emt_array_field("alteration.record_previous_capacitance", "emt.alteration", :output, "F", ["record"], getfield.(records, :previous_capacitance_f)),
        _portable_emt_array_field("alteration.record_resistance", "emt.alteration", :output, "ohm", ["record"], getfield.(records, :resistance_ohm)),
        _portable_emt_array_field("alteration.record_inductance", "emt.alteration", :output, "H", ["record"], getfield.(records, :inductance_h)),
        _portable_emt_array_field("alteration.record_capacitance", "emt.alteration", :output, "F", ["record"], getfield.(records, :capacitance_f)),
        _portable_emt_array_field("alteration.record_previous_conductance", "emt.alteration", :output, "S", ["record"], getfield.(records, :previous_conductance_s)),
        _portable_emt_array_field("alteration.record_conductance", "emt.alteration", :output, "S", ["record"], getfield.(records, :conductance_s)),
        _portable_emt_array_field("alteration.record_conductance_delta", "emt.alteration", :output, "S", ["record"], getfield.(records, :conductance_delta_s)),
        _portable_emt_state_field("alteration.record_history_preserved", "emt.alteration", :output, "1", getfield.(records, :history_state_preserved)),
        _portable_emt_array_field("alteration.record_refactor_count", "emt.alteration", :output, "1", ["record"], getfield.(records, :network_refactor_count)),
    ]
end

function _portable_emt_element_state_fields(
    element::TimeSwitch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "time_switch")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.close_time", "emt.switch", :discrete, "s", element.close_time_s),
        _portable_emt_state_field("$prefix.open_time", "emt.switch", :discrete, "s", element.open_time_s),
        _portable_emt_state_field("$prefix.initially_closed", "emt.switch", :discrete, "1", element.initially_closed),
    ])
    extinction = element.current_extinction
    if extinction !== nothing
        append!(fields, PortableSnapshotStateField[
            _portable_emt_state_field("$prefix.extinction_not_before_time", "emt.switch", :discrete, "s", extinction.not_before_time_s),
            _portable_emt_state_field("$prefix.extinction_critical_current", "emt.switch", :discrete, "A", extinction.critical_current_a),
            _portable_emt_state_field("$prefix.extinction_closed", "emt.switch", :discrete, "1", extinction.closed),
            _portable_emt_state_field("$prefix.extinction_opened", "emt.switch", :discrete, "1", extinction.opened),
            _portable_emt_state_field("$prefix.extinction_current_initialized", "emt.switch", :discrete, "1", extinction.current_initialized),
            _portable_emt_state_field("$prefix.extinction_previous_current", "emt.switch", :history, "A", extinction.previous_current),
            _portable_emt_state_field("$prefix.extinction_operation_count", "emt.switch", :discrete, "1", extinction.operation_count),
            _portable_emt_state_field("$prefix.extinction_open_reason", "emt.switch", :discrete, "1", String(extinction.open_reason)),
            _portable_emt_state_field("$prefix.extinction_opened_time", "emt.switch", :discrete, "s", extinction.opened_time_s),
        ])
    end
    return fields
end

function _portable_emt_element_state_fields(
    element::CurrentZeroSwitch,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "current_zero_switch")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.closed", "emt.switch", :discrete, "1", element.closed),
        _portable_emt_state_field("$prefix.opened", "emt.switch", :discrete, "1", element.opened),
        _portable_emt_state_field("$prefix.current_initialized", "emt.switch", :discrete, "1", element.current_initialized),
        _portable_emt_state_field("$prefix.previous_current", "emt.switch", :history, "A", element.previous_current),
        _portable_emt_state_field("$prefix.operation_count", "emt.switch", :discrete, "1", element.operation_count),
        _portable_emt_state_field("$prefix.open_reason", "emt.switch", :discrete, "1", String(element.open_reason)),
    ])
    return fields
end

function _portable_emt_element_state_fields(
    element,
    index::Integer,
    name::Symbol,
)
    ismutabletype(typeof(element)) && _portable_emt_fail(
        :unregistered_state_owner,
        "EMT model owner $(String(name)) has mutable state without a portable registration",
    )
    return _portable_emt_element_identity_field(index, name, "immutable_reconstructed")
end

function _portable_emt_source_state_fields(context::EMTStepContext)
    fields = PortableSnapshotStateField[]
    for (index, signal) in enumerate(context.analytic_source_signals)
        prefix = "source.i" * lpad(string(index), 8, '0')
        append!(fields, PortableSnapshotStateField[
            _portable_emt_state_field("$prefix.owner_name", "emt.source", :checkpoint, "1", String(context.analytic_source_names[index])),
            _portable_emt_state_field("$prefix.source_type", "emt.source", :discrete, "1", signal.source_type),
            _portable_emt_state_field("$prefix.crest", "emt.source", :continuous, "pu", signal.crest),
            _portable_emt_state_field("$prefix.time_parameter", "emt.source", :continuous, "s", signal.time1_s),
            _portable_emt_state_field("$prefix.angular_frequency_or_rate", "emt.source", :continuous, "rad_s_or_pu_s", signal.angular_frequency_or_rate),
            _portable_emt_time_field("$prefix.start_time", "emt.source", :scheduler, signal.start_time_s),
            _portable_emt_time_field("$prefix.stop_time", "emt.source", :scheduler, signal.stop_time_s),
        ])
    end
    return fields
end

function _portable_emt_nested_source_values(
    rows::AbstractVector{<:AbstractVector{<:Real}},
    identity::AbstractString,
)
    row_lengths = Int[length(row) for row in rows]
    flattened_values = Float64[]
    for (row_index, row) in enumerate(rows)
        converted = Float64.(row)
        all(isfinite, converted) || _portable_emt_fail(
            :nonfinite_value,
            "portable EMT source row $identity[$row_index] is nonfinite",
        )
        append!(flattened_values, converted)
    end
    return row_lengths, flattened_values
end

function _portable_emt_source_rows_identity_fields(
    identity::AbstractString,
    rows::AbstractVector{<:AbstractVector{<:Real}},
)
    row_lengths, flattened_values =
        _portable_emt_nested_source_values(rows, identity)
    return PortableSnapshotStateField[
        _portable_emt_array_field(
            "$identity.row_lengths",
            "emt.source_function.plan",
            :checkpoint,
            "1",
            ["row"],
            row_lengths,
        ),
        _portable_emt_array_field(
            "$identity.values",
            "emt.source_function.plan",
            :checkpoint,
            "pu",
            ["value"],
            flattened_values,
        ),
    ]
end

function _portable_emt_source_provider_identity_fields(
    ::IdentitySourceSignalProvider,
    prefix::AbstractString,
)
    return PortableSnapshotStateField[
        _portable_emt_state_field(
            "$prefix.kind",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            "identity",
        ),
    ]
end

function _portable_emt_source_provider_identity_fields(
    provider::TabulatedSourceSignalProvider,
    prefix::AbstractString,
)
    return PortableSnapshotStateField[
        _portable_emt_state_field(
            "$prefix.kind",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            "tabulated",
        ),
        _portable_emt_array_field(
            "$prefix.sample_time",
            "emt.source_function.provider",
            :checkpoint,
            "s",
            ["sample"],
            provider.times_s,
        ),
        _portable_emt_array_field(
            "$prefix.sample_value",
            "emt.source_function.provider",
            :checkpoint,
            "pu",
            ["slot", "sample"],
            provider.values,
        ),
        _portable_emt_state_field(
            "$prefix.extrapolation",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            String(provider.extrapolation),
        ),
    ]
end

function _portable_emt_source_provider_identity_fields(
    provider::SourceSignalProgram,
    prefix::AbstractString,
)
    fields = PortableSnapshotStateField[
        _portable_emt_state_field(
            "$prefix.kind",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            "program",
        ),
    ]
    append!(
        fields,
        _portable_emt_source_provider_identity_fields(
            provider.interpolation,
            "$prefix.interpolation",
        ),
    )
    slots = provider.analytic_slots
    signals = getfield.(slots, :signal)
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field(
            "$prefix.analytic_slot",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            ["assignment"],
            Int[slot.slot for slot in slots],
        ),
        _portable_emt_array_field(
            "$prefix.analytic_source_type",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            ["assignment"],
            Int[signal.source_type for signal in signals],
        ),
        _portable_emt_array_field(
            "$prefix.analytic_crest",
            "emt.source_function.provider",
            :checkpoint,
            "pu",
            ["assignment"],
            Float64[signal.crest for signal in signals],
        ),
        _portable_emt_array_field(
            "$prefix.analytic_time_parameter",
            "emt.source_function.provider",
            :checkpoint,
            "s",
            ["assignment"],
            Float64[signal.time1_s for signal in signals],
        ),
        _portable_emt_array_field(
            "$prefix.analytic_angular_frequency_or_rate",
            "emt.source_function.provider",
            :checkpoint,
            "rad_s_or_pu_s",
            ["assignment"],
            Float64[signal.angular_frequency_or_rate for signal in signals],
        ),
        _portable_emt_state_field(
            "$prefix.analytic_start_time",
            "emt.source_function.provider",
            :checkpoint,
            "s",
            _portable_emt_time_sequence(
                Float64[signal.start_time_s for signal in signals],
                "$prefix.analytic_start_time",
            );
            axes = ["assignment"],
        ),
        _portable_emt_state_field(
            "$prefix.analytic_stop_time",
            "emt.source_function.provider",
            :checkpoint,
            "s",
            _portable_emt_time_sequence(
                Float64[signal.stop_time_s for signal in signals],
                "$prefix.analytic_stop_time",
            );
            axes = ["assignment"],
        ),
        _portable_emt_state_field(
            "$prefix.analytic_assignment",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            String.(getfield.(slots, :assignment));
            axes = ["assignment"],
        ),
        _portable_emt_state_field(
            "$prefix.analytic_unit",
            "emt.source_function.provider",
            :checkpoint,
            "1",
            String.(getfield.(slots, :unit));
            axes = ["assignment"],
        ),
    ])
    return fields
end

function _portable_emt_source_provider_identity_fields(
    provider::AbstractSourceSignalProvider,
    prefix::AbstractString,
)
    _portable_emt_fail(
        :unregistered_state_owner,
        "source signal provider $(nameof(typeof(provider))) at $prefix has no portable registration",
    )
end

function _portable_emt_source_function_identity_fields(
    runtime::SourceFunctionNetworkRuntime,
)
    plan = runtime.plan
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("source_function.plan.source", "emt.source_function.plan", :checkpoint, "1", String(plan.source)),
        _portable_emt_state_field("source_function.plan.source_names", "emt.source_function.plan", :checkpoint, "1", String.(plan.source_names); axes = ["source"]),
        _portable_emt_state_field("source_function.plan.source_node_names", "emt.source_function.plan", :checkpoint, "1", String.(plan.source_node_names); axes = ["source"]),
        _portable_emt_array_field("source_function.plan.source_node_values", "emt.source_function.plan", :checkpoint, "1", ["source"], plan.source_node_values),
        _portable_emt_array_field("source_function.plan.source_iform_values", "emt.source_function.plan", :checkpoint, "1", ["source"], plan.source_iform_values),
        _portable_emt_array_field("source_function.plan.source_line_numbers", "emt.source_function.plan", :checkpoint, "1", ["source"], plan.source_line_numbers),
        _portable_emt_state_field("source_function.plan.source_layout_kinds", "emt.source_function.plan", :checkpoint, "1", String.(plan.source_layout_kinds); axes = ["source"]),
        _portable_emt_state_field("source_function.plan.source_row_count", "emt.source_function.plan", :checkpoint, "1", plan.source_row_count),
        _portable_emt_state_field("source_function.plan.source_card_kinds", "emt.source_function.plan", :checkpoint, "1", String.(plan.source_card_kinds); axes = ["row"]),
        _portable_emt_array_field("source_function.plan.source_card_provided_value_counts", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_card_provided_value_counts),
        _portable_emt_array_field("source_function.plan.source_card_line_numbers", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_card_line_numbers),
        _portable_emt_state_field("source_function.plan.source_card_row_count", "emt.source_function.plan", :checkpoint, "1", plan.source_card_row_count),
        _portable_emt_array_field("source_function.plan.source_interpolation_provided_value_counts", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_interpolation_provided_value_counts),
        _portable_emt_array_field("source_function.plan.source_interpolation_line_numbers", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_interpolation_line_numbers),
        _portable_emt_state_field("source_function.plan.source_interpolation_row_count", "emt.source_function.plan", :checkpoint, "1", plan.source_interpolation_row_count),
        _portable_emt_array_field("source_function.plan.source_tacs_override_positions", "emt.source_function.plan", :checkpoint, "1", ["override"], plan.source_tacs_override_positions),
        _portable_emt_array_field("source_function.plan.source_tacs_override_xtcs_indices", "emt.source_function.plan", :checkpoint, "1", ["override"], plan.source_tacs_override_xtcs_indices),
        _portable_emt_array_field("source_function.plan.source_tacs_override_line_numbers", "emt.source_function.plan", :checkpoint, "1", ["override"], plan.source_tacs_override_line_numbers),
        _portable_emt_state_field("source_function.plan.source_tacs_override_count", "emt.source_function.plan", :checkpoint, "1", plan.source_tacs_override_count),
        _portable_emt_array_field("source_function.plan.source_analytic_provided_value_counts", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_analytic_provided_value_counts),
        _portable_emt_array_field("source_function.plan.source_analytic_line_numbers", "emt.source_function.plan", :checkpoint, "1", ["row"], plan.source_analytic_line_numbers),
        _portable_emt_state_field("source_function.plan.source_analytic_row_count", "emt.source_function.plan", :checkpoint, "1", plan.source_analytic_row_count),
        _portable_emt_state_field("source_function.plan.fortran_scope", "emt.source_function.plan", :checkpoint, "1", String.(collect(plan.fortran_scope)); axes = ["owner"]),
        _portable_emt_state_field("source_function.plan.deferred_effects", "emt.source_function.plan", :checkpoint, "1", String.(collect(plan.deferred_effects)); axes = ["effect"]),
        _portable_emt_array_field("source_function.plan.dynamic_row_indices", "emt.source_function.plan", :checkpoint, "1", ["row"], runtime.dynamic_row_indices),
        _portable_emt_state_field("source_function.plan.internal_analytic_requested", "emt.source_function.plan", :checkpoint, "1", runtime.internal_analytic_requested),
    ]
    append!(fields, _portable_emt_source_rows_identity_fields(
        "source_function.plan.source_card",
        plan.source_card_values,
    ))
    append!(fields, _portable_emt_source_rows_identity_fields(
        "source_function.plan.source_interpolation",
        plan.source_interpolation_values,
    ))
    append!(fields, _portable_emt_source_rows_identity_fields(
        "source_function.plan.source_analytic",
        plan.source_analytic_values,
    ))
    append!(fields, _portable_emt_source_provider_identity_fields(
        runtime.signal_provider,
        "source_function.provider",
    ))
    return fields
end

function _portable_emt_source_stage_matrix(
    samples::AbstractVector{SourceSignalStageSample},
    field::Symbol,
)
    matrix = Matrix{Float64}(undef, 10, length(samples))
    for (sample_index, sample) in enumerate(samples)
        values = getfield(sample, field)
        length(values) == 10 || _portable_emt_fail(
            :state_shape_mismatch,
            "source-function stage $field does not contain ten slots",
        )
        all(isfinite, values) || _portable_emt_fail(
            :nonfinite_value,
            "source-function stage $field is nonfinite",
        )
        matrix[:, sample_index] .= values
    end
    return matrix
end

function _portable_emt_source_function_state_fields(
    runtime::Union{Nothing,SourceFunctionNetworkRuntime},
)
    fields = PortableSnapshotStateField[
        _portable_emt_state_field(
            "source_function.present",
            "emt.source_function",
            :checkpoint,
            "1",
            runtime !== nothing,
        ),
    ]
    runtime === nothing && return fields
    append!(fields, _portable_emt_source_function_identity_fields(runtime))
    samples = runtime.stage_samples
    assignment_counts = Int[length(sample.analytic_assignment_indices) for sample in samples]
    assignment_indices = Int[]
    for sample in samples
        all(index -> 1 <= index <= 10, sample.analytic_assignment_indices) ||
            _portable_emt_fail(
                :invalid_state_owner,
                "source-function analytic assignment index is outside the ten source slots",
            )
        append!(assignment_indices, sample.analytic_assignment_indices)
    end
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("source_function.plan.source_tstart_values", "emt.source_function.plan", :history, "s", ["source"], runtime.plan.source_tstart_values),
        _portable_emt_state_field(
            "source_function.plan.source_tstop_values",
            "emt.source_function.plan",
            :history,
            "s",
            _portable_emt_time_sequence(
                runtime.plan.source_tstop_values,
                "source_function.plan.source_tstop_values",
            );
            axes = ["source"],
        ),
        _portable_emt_array_field("source_function.plan.source_crest_values", "emt.source_function.plan", :history, "pu", ["source"], runtime.plan.source_crest_values),
        _portable_emt_array_field("source_function.plan.source_time1_values", "emt.source_function.plan", :history, "s", ["source"], runtime.plan.source_time1_values),
        _portable_emt_array_field("source_function.plan.source_time2_values", "emt.source_function.plan", :history, "s", ["source"], runtime.plan.source_time2_values),
        _portable_emt_array_field("source_function.plan.source_sfreq_values", "emt.source_function.plan", :history, "rad_s_or_pu_s", ["source"], runtime.plan.source_sfreq_values),
        _portable_emt_array_field("source_function.card.voltbc_values", "emt.source_function.card", :continuous, "pu", ["slot"], runtime.state.voltbc_values),
        _portable_emt_state_field("source_function.card.iread", "emt.source_function.card", :discrete, "1", runtime.state.iread),
        _portable_emt_state_field("source_function.card.nfrfld", "emt.source_function.card", :discrete, "1", runtime.state.nfrfld),
        _portable_emt_array_field("source_function.accepted_slot_values", "emt.source_function", :algebraic, "pu", ["slot"], Float64[value[] for value in runtime.slot_values]),
        _portable_emt_array_field("source_function.accepted_row_slot_values", "emt.source_function", :algebraic, "pu", ["row"], Float64[value[] for value in runtime.row_slot_values]),
        _portable_emt_state_field("source_function.executed_step_count", "emt.source_function", :output, "1", runtime.executed_step_count),
        _portable_emt_state_field("source_function.card_read_count", "emt.source_function", :output, "1", runtime.card_read_count),
        _portable_emt_state_field("source_function.signal_synchronization_count", "emt.source_function", :output, "1", runtime.signal_synchronization_count),
        _portable_emt_state_field("source_function.external_signal_count", "emt.source_function", :output, "1", runtime.external_signal_count),
        _portable_emt_state_field("source_function.tacs_override_count", "emt.source_function", :output, "1", runtime.tacs_override_count),
        _portable_emt_state_field("source_function.analytic_execution_count", "emt.source_function", :output, "1", runtime.analytic_execution_count),
        _portable_emt_time_field("source_function.last_accepted_time", "emt.source_function", :scheduler, runtime.last_accepted_time_s),
        _portable_emt_state_field("source_function.next_input_row_index", "emt.source_function", :scheduler, "1", runtime.next_input_row_index),
        _portable_emt_array_field("source_function.stage.step_index", "emt.source_function.stage", :output, "1", ["sample"], Int[sample.step_index for sample in samples]),
        _portable_emt_array_field("source_function.stage.time", "emt.source_function.stage", :output, "s", ["sample"], Float64[sample.time_s for sample in samples]),
        _portable_emt_array_field("source_function.stage.card_input", "emt.source_function.stage", :output, "pu", ["slot", "sample"], _portable_emt_source_stage_matrix(samples, :card_input_values)),
        _portable_emt_array_field("source_function.stage.interpolation_output", "emt.source_function.stage", :output, "pu", ["slot", "sample"], _portable_emt_source_stage_matrix(samples, :interpolation_output_values)),
        _portable_emt_array_field("source_function.stage.tacs_output", "emt.source_function.stage", :output, "pu", ["slot", "sample"], _portable_emt_source_stage_matrix(samples, :tacs_output_values)),
        _portable_emt_array_field("source_function.stage.analytic_output", "emt.source_function.stage", :output, "pu", ["slot", "sample"], _portable_emt_source_stage_matrix(samples, :analytic_output_values)),
        _portable_emt_array_field("source_function.stage.accepted", "emt.source_function.stage", :output, "pu", ["slot", "sample"], _portable_emt_source_stage_matrix(samples, :accepted_values)),
        _portable_emt_state_field("source_function.stage.card_read", "emt.source_function.stage", :output, "1", Bool[sample.card_read for sample in samples]; axes = ["sample"]),
        _portable_emt_state_field("source_function.stage.interpolation_applied", "emt.source_function.stage", :output, "1", Bool[sample.interpolation_applied for sample in samples]; axes = ["sample"]),
        _portable_emt_array_field("source_function.stage.tacs_override_count", "emt.source_function.stage", :output, "1", ["sample"], Int[sample.tacs_override_count for sample in samples]),
        _portable_emt_array_field("source_function.stage.analytic_assignment_count", "emt.source_function.stage", :output, "1", ["sample"], assignment_counts),
        _portable_emt_array_field("source_function.stage.analytic_assignment_index", "emt.source_function.stage", :output, "1", ["assignment"], assignment_indices),
    ])
    return fields
end

function _portable_emt_context_state_fields(context::EMTStepContext)
    accepted_step = min(max(context.step_index - 1, 0), context.step_count)
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("execution.timestep", "emt.execution", :checkpoint, "s", context.dt_s),
        _portable_emt_state_field("execution.horizon", "emt.execution", :checkpoint, "s", context.t_end_s),
        _portable_emt_state_field("execution.accepted_time", "emt.execution", :scheduler, "s", accepted_step * context.dt_s),
        _portable_emt_state_field("execution.next_time", "emt.execution", :scheduler, "s", context.t_s),
        _portable_emt_state_field("execution.accepted_step", "emt.execution", :scheduler, "1", accepted_step),
        _portable_emt_state_field("execution.next_step_index", "emt.execution", :scheduler, "1", context.step_index),
        _portable_emt_state_field("execution.horizon_step_count", "emt.execution", :checkpoint, "1", context.step_count),
        _portable_emt_state_field("output.trace_write_index", "emt.output", :output, "1", context.trace_write_index),
        _portable_emt_array_field("output.recorded_step_indices", "emt.output", :output, "1", ["sample"], context.recorded_step_indices),
        _portable_emt_array_field("network.accepted_node_voltage", "emt.network", :algebraic, "pu", ["node"], context.system.v),
        _portable_emt_array_field("network.accepted_rhs", "emt.network", :algebraic, "pu", ["node"], context.system.rhs),
        _portable_emt_array_field("history.electromagnetic_rhs", "emt.network", :history, "pu", ["node"], context.electromagnetic_history_rhs),
        _portable_emt_array_field("output.time", "emt.output", :output, "s", ["sample"], context.time_s),
        _portable_emt_array_field("output.node_voltage", "emt.output", :output, "pu", ["node", "sample"], context.voltage_pu),
        _portable_emt_array_field("output.channel_value", "emt.output", :output, "pu", ["channel", "sample"], context.output_pu),
        _portable_emt_array_field("output.accepted_step_value", "emt.output", :output, "pu", ["channel", "step"], context.output_step_values),
        _portable_emt_array_field("output.node_maximum", "emt.output", :output, "pu", ["node"], context.node_maximum_values),
        _portable_emt_array_field("output.node_maximum_time", "emt.output", :output, "s", ["node"], context.node_maximum_times_s),
        _portable_emt_array_field("output.node_minimum", "emt.output", :output, "pu", ["node"], context.node_minimum_values),
        _portable_emt_array_field("output.node_minimum_time", "emt.output", :output, "s", ["node"], context.node_minimum_times_s),
        _portable_emt_array_field("output.channel_maximum", "emt.output", :output, "pu", ["channel"], context.output_maximum_values),
        _portable_emt_array_field("output.channel_maximum_time", "emt.output", :output, "s", ["channel"], context.output_maximum_times_s),
        _portable_emt_array_field("output.channel_minimum", "emt.output", :output, "pu", ["channel"], context.output_minimum_values),
        _portable_emt_array_field("output.channel_minimum_time", "emt.output", :output, "s", ["channel"], context.output_minimum_times_s),
        _portable_emt_array_field("history.branch_previous_power", "emt.energy", :history, "W", ["branch"], context.branch_previous_power_values),
        _portable_emt_array_field("history.branch_energy", "emt.energy", :history, "J", ["branch"], context.branch_energy_values),
        _portable_emt_state_field("history.branch_power_valid", "emt.energy", :history, "1", collect(context.branch_power_history_valid)),
        _portable_emt_array_field("history.switch_previous_power", "emt.energy", :history, "W", ["switch"], context.switch_previous_power_values),
        _portable_emt_array_field("history.switch_energy", "emt.energy", :history, "J", ["switch"], context.switch_energy_values),
        _portable_emt_state_field("history.switch_power_valid", "emt.energy", :history, "1", collect(context.switch_power_history_valid)),
        _portable_emt_array_field("switch.closed_step_flags", "emt.switch", :discrete, "1", ["switch_step"], context.switch_closed_step_flags),
        _portable_emt_array_field("switch.conductance_step", "emt.switch", :history, "S", ["switch_step"], context.switch_conductance_step_values),
        _portable_emt_array_field("switch.voltage_step", "emt.switch", :history, "V", ["switch_step"], context.switch_voltage_step_values),
        _portable_emt_array_field("switch.current_step", "emt.switch", :history, "A", ["switch_step"], context.switch_current_step_values),
        _portable_emt_array_field("switch.power_step", "emt.switch", :history, "W", ["switch_step"], context.switch_power_step_values),
        _portable_emt_state_field("alteration.next_index", "emt.alteration", :scheduler, "1", context.next_series_rlc_alteration_index),
        _portable_emt_state_field("alteration.refactor_count", "emt.alteration", :reconstruction, "1", context.series_rlc_network_refactor_count),
        _portable_emt_state_field("identity.node_order", "emt.network", :checkpoint, "1", String.(context.node_names)),
        _portable_emt_state_field("identity.element_order", "emt.network", :checkpoint, "1", String.(context.element_names)),
        _portable_emt_state_field("identity.output_channel_order", "emt.output", :checkpoint, "1", String.(context.output_channel_names)),
    ]
    append!(fields, _portable_emt_source_state_fields(context))
    append!(fields, _portable_emt_source_function_state_fields(
        context.source_function_runtime,
    ))
    append!(fields, _portable_emt_control_system_state_fields(
        context.control_system_runtime,
    ))
    append!(fields, _portable_emt_alteration_state_fields(context))
    for (index, element) in enumerate(context.system.elements)
        append!(fields, _portable_emt_element_state_fields(
            element,
            index,
            context.element_names[index],
        ))
    end
    return fields
end

function _portable_emt_state_inventory(workspace::EMTStudyWorkspace)
    context = workspace.runtime.context
    fields = _portable_emt_context_state_fields(context)
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("workspace.evaluation_count", "emt.workspace", :output, "1", workspace.evaluation_count),
        _portable_emt_state_field("workspace.reset_count", "emt.workspace", :discrete, "1", workspace.reset_count),
        _portable_emt_state_field("workspace.ready", "emt.workspace", :discrete, "1", workspace.ready),
        _portable_emt_state_field("workspace.execution_mode", "emt.workspace", :discrete, "1", String(workspace.execution_mode)),
        _portable_emt_state_field("workspace.random_state_policy", "emt.workspace", :random, "1", "not_applicable"),
        _portable_emt_state_field("workspace.reduced_output_indices", "emt.workspace", :checkpoint, "1", _portable_emt_integer_array(workspace.reduced_output_indices, "1", ["channel"]); axes = ["channel"]),
        _portable_emt_state_field("workspace.source_signal_plan_indices", "emt.workspace", :checkpoint, "1", _portable_emt_integer_array(workspace.source_signal_plan_indices, "1", ["source"]); axes = ["source"]),
    ])
    return PortableSnapshotStateInventory(fields)
end

function _validate_portable_emt_workspace_boundary(workspace::EMTStudyWorkspace)
    workspace.ready && _portable_emt_fail(
        :unsynchronized_capture,
        "portable EMT capture requires an initialized accepted workspace boundary",
    )
    context = workspace.runtime.context
    1 <= context.step_index <= context.step_count + 1 || _portable_emt_fail(
        :unsynchronized_capture,
        "portable EMT next-step cursor is outside its fixed-step horizon",
    )
    accepted_step = min(max(context.step_index - 1, 0), context.step_count)
    expected_next_time = min(context.step_index, context.step_count) * context.dt_s
    context.t_s == expected_next_time || _portable_emt_fail(
        :unsynchronized_capture,
        "portable EMT next-time cursor is outside its accepted fixed-step boundary",
    )
    1 <= context.trace_write_index <= length(context.recorded_step_indices) + 1 ||
        _portable_emt_fail(
            :unsynchronized_capture,
            "portable EMT output cursor is outside its declared trace schedule",
        )
    previous_trace_steps = @view context.recorded_step_indices[1:(context.trace_write_index - 1)]
    future_trace_steps = @view context.recorded_step_indices[context.trace_write_index:end]
    all(step -> step <= accepted_step, previous_trace_steps) &&
        all(step -> step > accepted_step, future_trace_steps) || _portable_emt_fail(
            :unsynchronized_capture,
            "portable EMT output cursor disagrees with the accepted step",
        )
    all(isfinite, context.system.v) || _portable_emt_fail(
        :nonfinite_value,
        "portable EMT accepted node state is nonfinite",
    )
    _check_prepared_runtime_aliases(workspace.runtime)
    return accepted_step
end

"""Return the canonical public scientific-state inventory at one accepted EMT boundary."""
function portable_emt_state_inventory(workspace::EMTStudyWorkspace)
    workspace.execution_mode === :hybrid && _portable_emt_fail(
        :hybrid_owner_required,
        "portable hybrid execution must be captured through its coordinating integrator",
    )
    _validate_portable_emt_workspace_boundary(workspace)
    return _portable_emt_state_inventory(workspace)
end

function _portable_emt_inventory_fields(inventory::PortableSnapshotStateInventory)
    return Dict(field.identity => field for field in inventory.fields)
end

function _portable_emt_inventory_field(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
)
    key = String(identity)
    haskey(fields, key) || _portable_emt_fail(
        :missing_state_field,
        "portable EMT state field $key is missing",
    )
    field = fields[key]
    field.family == family || _portable_emt_fail(
        :state_family_mismatch,
        "portable EMT state field $key has the wrong state family",
    )
    field.unit == unit || _portable_emt_fail(
        :state_unit_mismatch,
        "portable EMT state field $key has the wrong unit",
    )
    return field
end

function _portable_emt_scalar(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
    ::Type{T},
) where {T}
    value = _portable_emt_inventory_field(fields, identity, family, unit).value
    value isa T || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT state field $identity has the wrong scalar type",
    )
    return value
end

function _portable_emt_time_scalar(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
)
    value = _portable_emt_inventory_field(fields, identity, family, "s").value
    value isa Float64 && return value
    value == _PORTABLE_EMT_POSITIVE_INFINITY && return Inf
    value == _PORTABLE_EMT_NEGATIVE_INFINITY && return -Inf
    _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT time field $identity has the wrong scalar type",
    )
end

function _portable_emt_time_sequence_values(
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
        value isa Float64 ? value :
        value == _PORTABLE_EMT_POSITIVE_INFINITY ? Inf :
        value == _PORTABLE_EMT_NEGATIVE_INFINITY ? -Inf :
        _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT time sequence $identity[$index] has an invalid value",
        ) for (index, value) in enumerate(field.value)
    ]
end

function _portable_emt_integer(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
)
    value = _portable_emt_inventory_field(fields, identity, family, "1").value
    value isa Integer || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT state field $identity is not an integer",
    )
    typemin(Int) <= value <= typemax(Int) || _portable_emt_fail(
        :integer_overflow,
        "portable EMT state field $identity does not fit this process",
    )
    return Int(value)
end

function _portable_emt_nonnegative_integer(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
)
    value = _portable_emt_integer(fields, identity, family)
    value >= 0 || _portable_emt_fail(
        :invalid_state_value,
        "portable EMT state field $identity must be nonnegative",
    )
    return value
end

function _restore_portable_emt_array!(
    destination::AbstractArray{T},
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
) where {T<:Union{Float64,Int}}
    field = _portable_emt_inventory_field(fields, identity, family, unit)
    field.axes == String.(axes) || _portable_emt_fail(
        :state_axis_mismatch,
        "portable EMT state field $identity has the wrong scientific axes",
    )
    field.value isa PortableSnapshotArray || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT state field $identity is not a portable array",
    )
    values = portable_snapshot_array_values(field.value)
    size(values) == size(destination) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT state field $identity does not fit its receiving owner",
    )
    if T === Float64
        values isa AbstractArray{Float64} || _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT state field $identity is not Float64",
        )
        copyto!(destination, values)
    else
        values isa AbstractArray{Int64} || _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT state field $identity is not Int64",
        )
        for index in eachindex(destination, values)
            typemin(Int) <= values[index] <= typemax(Int) || _portable_emt_fail(
                :integer_overflow,
                "portable EMT state field $identity contains an unaddressable integer",
            )
            destination[index] = Int(values[index])
        end
    end
    return destination
end

function _portable_emt_array_values(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
    ::Type{T},
) where {T<:Union{Float64,Int64}}
    field = _portable_emt_inventory_field(fields, identity, family, unit)
    field.axes == String.(axes) || _portable_emt_fail(
        :state_axis_mismatch,
        "portable EMT state field $identity has the wrong scientific axes",
    )
    field.value isa PortableSnapshotArray || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT state field $identity is not a portable array",
    )
    values = portable_snapshot_array_values(field.value)
    values isa AbstractArray{T} || _portable_emt_fail(
        :state_type_mismatch,
        "portable EMT state field $identity has the wrong array element type",
    )
    return vec(values)
end

function _portable_emt_string_sequence(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
)
    value = _portable_emt_inventory_field(fields, identity, family, "1").value
    value isa AbstractVector && all(item -> item isa AbstractString, value) ||
        _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT state field $identity is not a string sequence",
        )
    return String[item for item in value]
end

function _portable_emt_boolean_sequence(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol;
    axes::AbstractVector{<:AbstractString}=String[],
)
    field = _portable_emt_inventory_field(fields, identity, family, "1")
    field.axes == String.(axes) || _portable_emt_fail(
        :state_axis_mismatch,
        "portable EMT state field $identity has the wrong scientific axes",
    )
    value = field.value
    value isa AbstractVector && all(item -> item isa Bool, value) ||
        _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT state field $identity is not a Boolean sequence",
        )
    return Bool[item for item in value]
end

function _restore_portable_emt_boolean_vector!(
    destination::AbstractVector{Bool},
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
)
    value = _portable_emt_boolean_sequence(fields, identity, :history)
    length(value) == length(destination) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT state field $identity has the wrong Boolean count",
    )
    copyto!(destination, value)
    return destination
end

function _portable_emt_identity_sequence(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
)
    return _portable_emt_string_sequence(fields, identity, :checkpoint)
end

function _restore_portable_emt_element_state!(
    element::BergeronLine,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _restore_portable_emt_array!(element.from_wave_history, fields, "$prefix.from_wave_history", :delayed, "pu", ["delay"])
    _restore_portable_emt_array!(element.to_wave_history, fields, "$prefix.to_wave_history", :delayed, "pu", ["delay"])
    element.write_index = _portable_emt_integer(fields, "$prefix.write_index", :delayed)
    element.h_from = _portable_emt_scalar(fields, "$prefix.history_from", :history, "pu", Float64)
    element.h_to = _portable_emt_scalar(fields, "$prefix.history_to", :history, "pu", Float64)
    element.v_from = _portable_emt_scalar(fields, "$prefix.terminal_voltage_from", :algebraic, "pu", Float64)
    element.v_to = _portable_emt_scalar(fields, "$prefix.terminal_voltage_to", :algebraic, "pu", Float64)
    element.i_from = _portable_emt_scalar(fields, "$prefix.terminal_current_from", :algebraic, "pu", Float64)
    element.i_to = _portable_emt_scalar(fields, "$prefix.terminal_current_to", :algebraic, "pu", Float64)
    return element
end

function _restore_portable_emt_element_state!(
    element::SeriesRLBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.i_prev = _portable_emt_scalar(fields, "$prefix.previous_current", :history, "A", Float64)
    element.v_prev = _portable_emt_scalar(fields, "$prefix.previous_voltage", :history, "V", Float64)
    element.i_last = _portable_emt_scalar(fields, "$prefix.last_current", :algebraic, "A", Float64)
    return element
end

function _restore_portable_emt_element_state!(
    element::SeriesRLCBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.r = _portable_emt_scalar(fields, "$prefix.resistance", :discrete, "ohm", Float64)
    element.l = _portable_emt_scalar(fields, "$prefix.inductance", :discrete, "H", Float64)
    element.c = _portable_emt_scalar(fields, "$prefix.capacitance", :discrete, "F", Float64)
    element.i_prev = _portable_emt_scalar(fields, "$prefix.previous_current", :history, "A", Float64)
    element.inductor_voltage_prev = _portable_emt_scalar(fields, "$prefix.previous_inductor_voltage", :history, "V", Float64)
    element.capacitor_voltage_prev = _portable_emt_scalar(fields, "$prefix.previous_capacitor_voltage", :history, "V", Float64)
    element.v_prev = _portable_emt_scalar(fields, "$prefix.previous_terminal_voltage", :history, "V", Float64)
    element.i_last = _portable_emt_scalar(fields, "$prefix.last_current", :algebraic, "A", Float64)
    return element
end

function _restore_portable_emt_element_state!(
    element::CapacitorBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.c = _portable_emt_scalar(fields, "$prefix.capacitance", :checkpoint, "F", Float64)
    element.i_prev = _portable_emt_scalar(fields, "$prefix.previous_current", :history, "A", Float64)
    element.v_prev = _portable_emt_scalar(fields, "$prefix.previous_voltage", :history, "V", Float64)
    element.i_last = _portable_emt_scalar(fields, "$prefix.last_current", :algebraic, "A", Float64)
    return element
end

function _restore_portable_emt_coupled_branch_workspaces!(
    element,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _restore_portable_emt_array!(element.previous_current, fields, "$prefix.previous_current", :history, "A", ["port"])
    _restore_portable_emt_array!(element.previous_voltage, fields, "$prefix.previous_voltage", :history, "V", ["port"])
    _restore_portable_emt_array!(element.last_current, fields, "$prefix.last_current", :algebraic, "A", ["port"])
    _restore_portable_emt_array!(element.conductance_workspace, fields, "$prefix.conductance_workspace", :algebraic, "S", ["row", "column"])
    _restore_portable_emt_array!(element.history_current_workspace, fields, "$prefix.history_current_workspace", :algebraic, "A", ["port"])
    _restore_portable_emt_array!(element.port_voltage_workspace, fields, "$prefix.port_voltage_workspace", :algebraic, "V", ["port"])
    _restore_portable_emt_array!(element.current_workspace, fields, "$prefix.current_workspace", :algebraic, "A", ["port"])
    return element
end

function _restore_portable_emt_element_state!(
    element::CoupledInductiveBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    return _restore_portable_emt_coupled_branch_workspaces!(element, fields, index)
end

function _restore_portable_emt_element_state!(
    element::CoupledSeriesRLBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    _restore_portable_emt_coupled_branch_workspaces!(element, fields, index)
    prefix = _portable_emt_element_prefix(index)
    cached = _portable_emt_inventory_field(
        fields,
        "$prefix.cached_timestep",
        :reconstruction,
        "s",
    ).value
    element.cached_dt_s = cached == "not_initialized" ? NaN :
        cached isa Float64 && isfinite(cached) && cached > 0.0 ? cached :
        _portable_emt_fail(
            :state_type_mismatch,
            "portable coupled series R-L cached timestep is invalid",
        )
    return element
end

function _restore_portable_emt_element_state!(
    element::SaturableInductorBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.i_prev = _portable_emt_scalar(fields, "$prefix.previous_current", :history, "A", Float64)
    element.v_prev = _portable_emt_scalar(fields, "$prefix.previous_voltage", :history, "V", Float64)
    element.i_last = _portable_emt_scalar(fields, "$prefix.last_current", :algebraic, "A", Float64)
    element.last_inductance = _portable_emt_scalar(fields, "$prefix.last_inductance", :discrete, "H", Float64)
    element.saturated = _portable_emt_scalar(fields, "$prefix.saturated", :discrete, "1", Bool)
    return element
end

function _restore_portable_emt_element_state!(
    element::SaturatedTransformerNonlinearSlopeBranch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.conductance = _portable_emt_scalar(
        fields,
        "$prefix.conductance",
        :algebraic,
        "S",
        Float64,
    )
    return element
end

function _restore_portable_emt_element_state!(
    element::TimeSwitch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.close_time_s = _portable_emt_scalar(fields, "$prefix.close_time", :discrete, "s", Float64)
    element.open_time_s = _portable_emt_scalar(fields, "$prefix.open_time", :discrete, "s", Float64)
    element.initially_closed = _portable_emt_scalar(fields, "$prefix.initially_closed", :discrete, "1", Bool)
    extinction = element.current_extinction
    if extinction !== nothing
        extinction.not_before_time_s = _portable_emt_scalar(fields, "$prefix.extinction_not_before_time", :discrete, "s", Float64)
        extinction.critical_current_a = _portable_emt_scalar(fields, "$prefix.extinction_critical_current", :discrete, "A", Float64)
        extinction.closed = _portable_emt_scalar(fields, "$prefix.extinction_closed", :discrete, "1", Bool)
        extinction.opened = _portable_emt_scalar(fields, "$prefix.extinction_opened", :discrete, "1", Bool)
        extinction.current_initialized = _portable_emt_scalar(fields, "$prefix.extinction_current_initialized", :discrete, "1", Bool)
        extinction.previous_current = _portable_emt_scalar(fields, "$prefix.extinction_previous_current", :history, "A", Float64)
        extinction.operation_count = _portable_emt_integer(fields, "$prefix.extinction_operation_count", :discrete)
        reason = _portable_emt_scalar(fields, "$prefix.extinction_open_reason", :discrete, "1", String)
        extinction.open_reason = Symbol(reason)
        extinction.opened_time_s = _portable_emt_scalar(fields, "$prefix.extinction_opened_time", :discrete, "s", Float64)
    end
    return element
end

function _restore_portable_emt_element_state!(
    element::CurrentZeroSwitch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    element.closed = _portable_emt_scalar(fields, "$prefix.closed", :discrete, "1", Bool)
    element.opened = _portable_emt_scalar(fields, "$prefix.opened", :discrete, "1", Bool)
    element.current_initialized = _portable_emt_scalar(fields, "$prefix.current_initialized", :discrete, "1", Bool)
    element.previous_current = _portable_emt_scalar(fields, "$prefix.previous_current", :history, "A", Float64)
    element.operation_count = _portable_emt_integer(fields, "$prefix.operation_count", :discrete)
    element.open_reason = Symbol(_portable_emt_scalar(fields, "$prefix.open_reason", :discrete, "1", String))
    return element
end

_restore_portable_emt_element_state!(element, fields, index::Integer) = element

function _portable_emt_validate_source_function_identity(
    runtime::SourceFunctionNetworkRuntime,
    fields::Dict{String,PortableSnapshotStateField},
)
    reconstructed = PortableSnapshotStateInventory(
        _portable_emt_source_function_identity_fields(runtime),
    )
    captured = PortableSnapshotStateInventory(PortableSnapshotStateField[
        _portable_emt_inventory_field(
            fields,
            field.identity,
            field.family,
            field.unit,
        ) for field in reconstructed.fields
    ])
    captured.signature_sha256 == reconstructed.signature_sha256 ||
        _portable_emt_fail(
            :settings_mismatch,
            "portable EMT source-function plan or provider changed in the receiving study",
        )
    return runtime
end

function _portable_emt_source_stage_matrix(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    sample_count::Integer,
)
    values = zeros(Float64, 10, sample_count)
    _restore_portable_emt_array!(
        values,
        fields,
        identity,
        :output,
        "pu",
        ["slot", "sample"],
    )
    return values
end

function _restore_portable_emt_source_function_state!(
    context::EMTStepContext,
    fields::Dict{String,PortableSnapshotStateField},
)
    present = _portable_emt_scalar(
        fields,
        "source_function.present",
        :checkpoint,
        "1",
        Bool,
    )
    runtime = context.source_function_runtime
    (runtime !== nothing) == present || _portable_emt_fail(
        :settings_mismatch,
        "portable EMT source-function ownership changed in the receiving study",
    )
    runtime === nothing && return context
    _portable_emt_validate_source_function_identity(runtime, fields)
    _restore_portable_emt_array!(
        runtime.plan.source_tstart_values,
        fields,
        "source_function.plan.source_tstart_values",
        :history,
        "s",
        ["source"],
    )
    source_tstop_values = _portable_emt_time_sequence_values(
        fields,
        "source_function.plan.source_tstop_values",
        :history,
        "source",
    )
    length(source_tstop_values) == length(runtime.plan.source_tstop_values) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT source stop-time history has the wrong source count",
        )
    copyto!(runtime.plan.source_tstop_values, source_tstop_values)
    _restore_portable_emt_array!(
        runtime.plan.source_crest_values,
        fields,
        "source_function.plan.source_crest_values",
        :history,
        "pu",
        ["source"],
    )
    _restore_portable_emt_array!(
        runtime.plan.source_time1_values,
        fields,
        "source_function.plan.source_time1_values",
        :history,
        "s",
        ["source"],
    )
    _restore_portable_emt_array!(
        runtime.plan.source_time2_values,
        fields,
        "source_function.plan.source_time2_values",
        :history,
        "s",
        ["source"],
    )
    _restore_portable_emt_array!(
        runtime.plan.source_sfreq_values,
        fields,
        "source_function.plan.source_sfreq_values",
        :history,
        "rad_s_or_pu_s",
        ["source"],
    )
    _restore_portable_emt_array!(
        runtime.state.voltbc_values,
        fields,
        "source_function.card.voltbc_values",
        :continuous,
        "pu",
        ["slot"],
    )
    runtime.state.iread = _portable_emt_nonnegative_integer(
        fields,
        "source_function.card.iread",
        :discrete,
    )
    runtime.state.nfrfld = _portable_emt_nonnegative_integer(
        fields,
        "source_function.card.nfrfld",
        :discrete,
    )
    accepted_slot_values = Float64[value[] for value in runtime.slot_values]
    _restore_portable_emt_array!(
        accepted_slot_values,
        fields,
        "source_function.accepted_slot_values",
        :algebraic,
        "pu",
        ["slot"],
    )
    for index in eachindex(runtime.slot_values, accepted_slot_values)
        runtime.slot_values[index][] = accepted_slot_values[index]
    end
    accepted_row_slot_values = Float64[value[] for value in runtime.row_slot_values]
    _restore_portable_emt_array!(
        accepted_row_slot_values,
        fields,
        "source_function.accepted_row_slot_values",
        :algebraic,
        "pu",
        ["row"],
    )
    for index in eachindex(runtime.row_slot_values, accepted_row_slot_values)
        runtime.row_slot_values[index][] = accepted_row_slot_values[index]
    end
    runtime.executed_step_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.executed_step_count",
        :output,
    )
    runtime.card_read_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.card_read_count",
        :output,
    )
    runtime.signal_synchronization_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.signal_synchronization_count",
        :output,
    )
    runtime.external_signal_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.external_signal_count",
        :output,
    )
    runtime.tacs_override_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.tacs_override_count",
        :output,
    )
    runtime.analytic_execution_count = _portable_emt_nonnegative_integer(
        fields,
        "source_function.analytic_execution_count",
        :output,
    )
    runtime.last_accepted_time_s = _portable_emt_time_scalar(
        fields,
        "source_function.last_accepted_time",
        :scheduler,
    )
    runtime.next_input_row_index = _portable_emt_nonnegative_integer(
        fields,
        "source_function.next_input_row_index",
        :scheduler,
    )
    runtime.next_input_row_index >= 1 || _portable_emt_fail(
        :invalid_state_value,
        "portable EMT source-function input cursor must be positive",
    )

    step_indices = _portable_emt_array_values(
        fields,
        "source_function.stage.step_index",
        :output,
        "1",
        ["sample"],
        Int64,
    )
    sample_count = length(step_indices)
    times = _portable_emt_array_values(
        fields,
        "source_function.stage.time",
        :output,
        "s",
        ["sample"],
        Float64,
    )
    card_inputs = _portable_emt_source_stage_matrix(
        fields,
        "source_function.stage.card_input",
        sample_count,
    )
    interpolation_outputs = _portable_emt_source_stage_matrix(
        fields,
        "source_function.stage.interpolation_output",
        sample_count,
    )
    tacs_outputs = _portable_emt_source_stage_matrix(
        fields,
        "source_function.stage.tacs_output",
        sample_count,
    )
    analytic_outputs = _portable_emt_source_stage_matrix(
        fields,
        "source_function.stage.analytic_output",
        sample_count,
    )
    accepted_values = _portable_emt_source_stage_matrix(
        fields,
        "source_function.stage.accepted",
        sample_count,
    )
    card_reads = _portable_emt_boolean_sequence(
        fields,
        "source_function.stage.card_read",
        :output;
        axes = ["sample"],
    )
    interpolation_applied = _portable_emt_boolean_sequence(
        fields,
        "source_function.stage.interpolation_applied",
        :output;
        axes = ["sample"],
    )
    tacs_override_counts = _portable_emt_array_values(
        fields,
        "source_function.stage.tacs_override_count",
        :output,
        "1",
        ["sample"],
        Int64,
    )
    assignment_counts = _portable_emt_array_values(
        fields,
        "source_function.stage.analytic_assignment_count",
        :output,
        "1",
        ["sample"],
        Int64,
    )
    assignment_indices = _portable_emt_array_values(
        fields,
        "source_function.stage.analytic_assignment_index",
        :output,
        "1",
        ["assignment"],
        Int64,
    )
    all(length(values) == sample_count for values in (
        times,
        card_reads,
        interpolation_applied,
        tacs_override_counts,
        assignment_counts,
    )) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT source-function stage fields have inconsistent sample counts",
    )
    all(count -> count >= 0, assignment_counts) || _portable_emt_fail(
        :invalid_state_value,
        "portable EMT source-function assignment counts must be nonnegative",
    )
    sum(assignment_counts; init = Int64(0)) == length(assignment_indices) ||
        _portable_emt_fail(
            :state_shape_mismatch,
            "portable EMT source-function assignments have inconsistent counts",
        )
    all(index -> 1 <= index <= 10, assignment_indices) || _portable_emt_fail(
        :invalid_state_value,
        "portable EMT source-function assignment index is outside the ten source slots",
    )
    empty!(runtime.stage_samples)
    assignment_offset = 0
    for sample_index in 1:sample_count
        assignment_count = _portable_emt_checked_integer(
            assignment_counts[sample_index],
            "source_function.stage.analytic_assignment_count",
        )
        sample_assignments = Int[
            _portable_emt_checked_integer(
                assignment_indices[index],
                "source_function.stage.analytic_assignment_index",
            ) for index in (assignment_offset + 1):(assignment_offset + assignment_count)
        ]
        assignment_offset += assignment_count
        push!(
            runtime.stage_samples,
            SourceSignalStageSample(
                _portable_emt_checked_integer(
                    step_indices[sample_index],
                    "source_function.stage.step_index",
                ),
                times[sample_index],
                copy(card_inputs[:, sample_index]),
                copy(interpolation_outputs[:, sample_index]),
                copy(tacs_outputs[:, sample_index]),
                copy(analytic_outputs[:, sample_index]),
                copy(accepted_values[:, sample_index]),
                card_reads[sample_index],
                interpolation_applied[sample_index],
                _portable_emt_checked_integer(
                    tacs_override_counts[sample_index],
                    "source_function.stage.tacs_override_count",
                ),
                sample_assignments,
            ),
        )
    end
    return context
end

function _restore_portable_emt_source_state!(
    context::EMTStepContext,
    fields::Dict{String,PortableSnapshotStateField},
)
    for (index, signal) in enumerate(context.analytic_source_signals)
        prefix = "source.i" * lpad(string(index), 8, '0')
        signal.source_type = _portable_emt_integer(fields, "$prefix.source_type", :discrete)
        signal.crest = _portable_emt_scalar(fields, "$prefix.crest", :continuous, "pu", Float64)
        signal.time1_s = _portable_emt_scalar(fields, "$prefix.time_parameter", :continuous, "s", Float64)
        signal.angular_frequency_or_rate = _portable_emt_scalar(fields, "$prefix.angular_frequency_or_rate", :continuous, "rad_s_or_pu_s", Float64)
        signal.start_time_s = _portable_emt_time_scalar(fields, "$prefix.start_time", :scheduler)
        signal.stop_time_s = _portable_emt_time_scalar(fields, "$prefix.stop_time", :scheduler)
    end
    return context
end


function _portable_emt_checked_integer(value::Int64, identity::AbstractString)
    typemin(Int) <= value <= typemax(Int) || _portable_emt_fail(
        :integer_overflow,
        "portable EMT state field $identity does not fit this process",
    )
    return Int(value)
end

function _restore_portable_emt_alteration_state!(
    context::EMTStepContext,
    fields::Dict{String,PortableSnapshotStateField},
)
    isempty(context.series_rlc_alterations) && return context
    schedule_names = _portable_emt_string_sequence(
        fields,
        "alteration.schedule_branch_names",
        :checkpoint,
    )
    schedule_times = _portable_emt_array_values(fields, "alteration.schedule_activation_time", :checkpoint, "s", ["alteration"], Float64)
    schedule_resistance = _portable_emt_array_values(fields, "alteration.schedule_resistance", :checkpoint, "ohm", ["alteration"], Float64)
    schedule_inductance = _portable_emt_array_values(fields, "alteration.schedule_inductance", :checkpoint, "H", ["alteration"], Float64)
    schedule_capacitance = _portable_emt_array_values(fields, "alteration.schedule_capacitance", :checkpoint, "F", ["alteration"], Float64)
    schedule_indices = _portable_emt_array_values(fields, "alteration.schedule_branch_indices", :checkpoint, "1", ["alteration"], Int64)
    schedule_count = length(context.series_rlc_alterations)
    all(length(values) == schedule_count for values in (
        schedule_names,
        schedule_times,
        schedule_resistance,
        schedule_inductance,
        schedule_capacitance,
        schedule_indices,
    )) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT alteration schedule has inconsistent field counts",
    )
    restored_schedule = SeriesRLCAlteration[
        SeriesRLCAlteration(
            schedule_names[index],
            schedule_times[index],
            schedule_resistance[index],
            schedule_inductance[index],
            schedule_capacitance[index],
        ) for index in 1:schedule_count
    ]
    restored_indices = Int[
        _portable_emt_checked_integer(value, "alteration.schedule_branch_indices")
        for value in schedule_indices
    ]
    restored_schedule == context.series_rlc_alterations &&
        restored_indices == context.series_rlc_alteration_branch_indices ||
        _portable_emt_fail(
            :settings_mismatch,
            "portable EMT alteration schedule changed in the receiving study",
        )

    record_names = _portable_emt_string_sequence(
        fields,
        "alteration.record_branch_names",
        :output,
    )
    record_indices = _portable_emt_array_values(fields, "alteration.record_branch_indices", :output, "1", ["record"], Int64)
    requested_times = _portable_emt_array_values(fields, "alteration.record_requested_time", :output, "s", ["record"], Float64)
    applied_times = _portable_emt_array_values(fields, "alteration.record_applied_time", :output, "s", ["record"], Float64)
    applied_steps = _portable_emt_array_values(fields, "alteration.record_applied_step", :output, "1", ["record"], Int64)
    previous_resistance = _portable_emt_array_values(fields, "alteration.record_previous_resistance", :output, "ohm", ["record"], Float64)
    previous_inductance = _portable_emt_array_values(fields, "alteration.record_previous_inductance", :output, "H", ["record"], Float64)
    previous_capacitance = _portable_emt_array_values(fields, "alteration.record_previous_capacitance", :output, "F", ["record"], Float64)
    resistance = _portable_emt_array_values(fields, "alteration.record_resistance", :output, "ohm", ["record"], Float64)
    inductance = _portable_emt_array_values(fields, "alteration.record_inductance", :output, "H", ["record"], Float64)
    capacitance = _portable_emt_array_values(fields, "alteration.record_capacitance", :output, "F", ["record"], Float64)
    previous_conductance = _portable_emt_array_values(fields, "alteration.record_previous_conductance", :output, "S", ["record"], Float64)
    conductance = _portable_emt_array_values(fields, "alteration.record_conductance", :output, "S", ["record"], Float64)
    conductance_delta = _portable_emt_array_values(fields, "alteration.record_conductance_delta", :output, "S", ["record"], Float64)
    history_preserved = _portable_emt_inventory_field(
        fields,
        "alteration.record_history_preserved",
        :output,
        "1",
    ).value
    history_preserved isa AbstractVector &&
        all(value -> value isa Bool, history_preserved) || _portable_emt_fail(
            :state_type_mismatch,
            "portable EMT alteration history flags are not Boolean",
        )
    refactor_counts = _portable_emt_array_values(fields, "alteration.record_refactor_count", :output, "1", ["record"], Int64)
    record_count = length(record_names)
    all(length(values) == record_count for values in (
        record_indices,
        requested_times,
        applied_times,
        applied_steps,
        previous_resistance,
        previous_inductance,
        previous_capacitance,
        resistance,
        inductance,
        capacitance,
        previous_conductance,
        conductance,
        conductance_delta,
        history_preserved,
        refactor_counts,
    )) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable EMT alteration records have inconsistent field counts",
    )
    empty!(context.series_rlc_alteration_records)
    for index in 1:record_count
        push!(
            context.series_rlc_alteration_records,
            SeriesRLCAlterationRecord(
                Symbol(record_names[index]),
                _portable_emt_checked_integer(record_indices[index], "alteration.record_branch_indices"),
                requested_times[index],
                applied_times[index],
                _portable_emt_checked_integer(applied_steps[index], "alteration.record_applied_step"),
                previous_resistance[index],
                previous_inductance[index],
                previous_capacitance[index],
                resistance[index],
                inductance[index],
                capacitance[index],
                previous_conductance[index],
                conductance[index],
                conductance_delta[index],
                history_preserved[index],
                _portable_emt_checked_integer(refactor_counts[index], "alteration.record_refactor_count"),
            ),
        )
    end
    return context
end

"""Restore only registered public scientific state into an isolated, exactly prepared candidate."""
function restore_portable_emt_state_inventory!(
    candidate::EMTStudyWorkspace,
    inventory::PortableSnapshotStateInventory,
)
    fields = _portable_emt_inventory_fields(inventory)
    context = candidate.runtime.context
    _portable_emt_scalar(
        fields,
        "workspace.random_state_policy",
        :random,
        "1",
        String,
    ) == "not_applicable" || _portable_emt_fail(
        :random_state_mismatch,
        "portable deterministic EMT workspace has an unsupported random-state policy",
    )
    _portable_emt_identity_sequence(fields, "identity.node_order") == String.(context.node_names) ||
        _portable_emt_fail(:topology_mismatch, "portable EMT node order changed")
    _portable_emt_identity_sequence(fields, "identity.element_order") == String.(context.element_names) ||
        _portable_emt_fail(:topology_mismatch, "portable EMT element order changed")
    _portable_emt_identity_sequence(fields, "identity.output_channel_order") == String.(context.output_channel_names) ||
        _portable_emt_fail(:output_mismatch, "portable EMT output channel order changed")
    execution_mode = Symbol(_portable_emt_scalar(
        fields,
        "workspace.execution_mode",
        :discrete,
        "1",
        String,
    ))
    partition_execution_modes = (
        :local_subcycling,
        :partitioned_lagged,
        :partitioned_waveform,
        :combined_local_partitioned,
    )
    execution_mode in (:monolithic, :hybrid, partition_execution_modes...) ||
        _portable_emt_fail(
            :execution_mode_mismatch,
            "portable EMT execution mode is unsupported",
        )
    if execution_mode === :hybrid && candidate.execution_mode !== :hybrid
        code = candidate.execution_mode === :unselected ?
            :hybrid_owner_required : :execution_mode_mismatch
        _portable_emt_fail(
            code,
            "portable hybrid workspace state requires its coordinating integrator",
        )
    elseif execution_mode === :monolithic &&
           !(candidate.execution_mode in (:unselected, :monolithic))
        _portable_emt_fail(
            :execution_mode_mismatch,
            "portable monolithic state cannot replace a hybrid execution owner",
        )
    elseif execution_mode in partition_execution_modes &&
           !(
               candidate.execution_mode === execution_mode &&
               !candidate.ready
           )
        _portable_emt_fail(
            :partition_owner_required,
            "portable partition-region state requires its coordinating partition owner",
        )
    end
    context.dt_s = _portable_emt_scalar(fields, "execution.timestep", :checkpoint, "s", Float64)
    context.t_end_s = _portable_emt_scalar(fields, "execution.horizon", :checkpoint, "s", Float64)
    accepted_time = _portable_emt_scalar(fields, "execution.accepted_time", :scheduler, "s", Float64)
    context.t_s = _portable_emt_scalar(fields, "execution.next_time", :scheduler, "s", Float64)
    accepted_step = _portable_emt_integer(fields, "execution.accepted_step", :scheduler)
    context.step_index = _portable_emt_integer(fields, "execution.next_step_index", :scheduler)
    context.step_count = _portable_emt_integer(fields, "execution.horizon_step_count", :checkpoint)
    accepted_step == min(max(context.step_index - 1, 0), context.step_count) &&
        accepted_time == accepted_step * context.dt_s || _portable_emt_fail(
            :accepted_step_mismatch,
            "portable EMT accepted step, time, and next cursor disagree",
        )
    context.trace_write_index = _portable_emt_integer(fields, "output.trace_write_index", :output)
    _restore_portable_emt_array!(context.recorded_step_indices, fields, "output.recorded_step_indices", :output, "1", ["sample"])
    _restore_portable_emt_array!(context.system.v, fields, "network.accepted_node_voltage", :algebraic, "pu", ["node"])
    _restore_portable_emt_array!(context.system.rhs, fields, "network.accepted_rhs", :algebraic, "pu", ["node"])
    _restore_portable_emt_array!(context.electromagnetic_history_rhs, fields, "history.electromagnetic_rhs", :history, "pu", ["node"])
    _restore_portable_emt_array!(context.time_s, fields, "output.time", :output, "s", ["sample"])
    _restore_portable_emt_array!(context.voltage_pu, fields, "output.node_voltage", :output, "pu", ["node", "sample"])
    _restore_portable_emt_array!(context.output_pu, fields, "output.channel_value", :output, "pu", ["channel", "sample"])
    _restore_portable_emt_array!(context.output_step_values, fields, "output.accepted_step_value", :output, "pu", ["channel", "step"])
    _restore_portable_emt_array!(context.node_maximum_values, fields, "output.node_maximum", :output, "pu", ["node"])
    _restore_portable_emt_array!(context.node_maximum_times_s, fields, "output.node_maximum_time", :output, "s", ["node"])
    _restore_portable_emt_array!(context.node_minimum_values, fields, "output.node_minimum", :output, "pu", ["node"])
    _restore_portable_emt_array!(context.node_minimum_times_s, fields, "output.node_minimum_time", :output, "s", ["node"])
    _restore_portable_emt_array!(context.output_maximum_values, fields, "output.channel_maximum", :output, "pu", ["channel"])
    _restore_portable_emt_array!(context.output_maximum_times_s, fields, "output.channel_maximum_time", :output, "s", ["channel"])
    _restore_portable_emt_array!(context.output_minimum_values, fields, "output.channel_minimum", :output, "pu", ["channel"])
    _restore_portable_emt_array!(context.output_minimum_times_s, fields, "output.channel_minimum_time", :output, "s", ["channel"])
    _restore_portable_emt_array!(context.branch_previous_power_values, fields, "history.branch_previous_power", :history, "W", ["branch"])
    _restore_portable_emt_array!(context.branch_energy_values, fields, "history.branch_energy", :history, "J", ["branch"])
    _restore_portable_emt_boolean_vector!(context.branch_power_history_valid, fields, "history.branch_power_valid")
    _restore_portable_emt_array!(context.switch_previous_power_values, fields, "history.switch_previous_power", :history, "W", ["switch"])
    _restore_portable_emt_array!(context.switch_energy_values, fields, "history.switch_energy", :history, "J", ["switch"])
    _restore_portable_emt_boolean_vector!(context.switch_power_history_valid, fields, "history.switch_power_valid")
    _restore_portable_emt_array!(context.switch_closed_step_flags, fields, "switch.closed_step_flags", :discrete, "1", ["switch_step"])
    _restore_portable_emt_array!(context.switch_conductance_step_values, fields, "switch.conductance_step", :history, "S", ["switch_step"])
    _restore_portable_emt_array!(context.switch_voltage_step_values, fields, "switch.voltage_step", :history, "V", ["switch_step"])
    _restore_portable_emt_array!(context.switch_current_step_values, fields, "switch.current_step", :history, "A", ["switch_step"])
    _restore_portable_emt_array!(context.switch_power_step_values, fields, "switch.power_step", :history, "W", ["switch_step"])
    context.next_series_rlc_alteration_index = _portable_emt_integer(fields, "alteration.next_index", :scheduler)
    context.series_rlc_network_refactor_count = _portable_emt_integer(fields, "alteration.refactor_count", :reconstruction)
    _restore_portable_emt_alteration_state!(context, fields)
    _restore_portable_emt_source_state!(context, fields)
    _restore_portable_emt_source_function_state!(context, fields)
    _restore_portable_emt_control_system_state!(context, fields)
    for (index, element) in enumerate(context.system.elements)
        prefix = _portable_emt_element_prefix(index)
        _portable_emt_scalar(fields, "$prefix.owner_name", :checkpoint, "1", String) == String(context.element_names[index]) ||
            _portable_emt_fail(:model_mismatch, "portable EMT model owner order changed")
        _restore_portable_emt_element_state!(element, fields, index)
    end
    _restore_portable_emt_array!(candidate.reduced_output_indices, fields, "workspace.reduced_output_indices", :checkpoint, "1", ["channel"])
    _restore_portable_emt_array!(candidate.source_signal_plan_indices, fields, "workspace.source_signal_plan_indices", :checkpoint, "1", ["source"])
    candidate.evaluation_count = _portable_emt_integer(fields, "workspace.evaluation_count", :output)
    candidate.reset_count = _portable_emt_integer(fields, "workspace.reset_count", :discrete)
    candidate.ready = _portable_emt_scalar(fields, "workspace.ready", :discrete, "1", Bool)
    candidate.execution_mode = execution_mode
    restored = _portable_emt_state_inventory(candidate)
    getfield.(restored.fields, :identity) == getfield.(inventory.fields, :identity) &&
        restored.signature_sha256 == inventory.signature_sha256 ||
        _portable_emt_fail(
            :state_inventory_mismatch,
            "portable EMT state inventory does not match the isolated candidate",
        )
    return candidate
end

struct PortableEMTRestoreResult{W}
    workspace::W
    descriptor::PortableSnapshotDescriptor
    public_state_signature_sha256::String
    backend_state_signature_sha256::String
    reconstructed::Bool
end

function _portable_emt_parent_function(name::Symbol)
    parent = parentmodule(@__MODULE__)
    isdefined(parent, name) || _portable_emt_fail(
        :backend_unavailable,
        "portable EMT backend operation $name is unavailable",
    )
    return getfield(parent, name)
end

function _portable_emt_backend_snapshot(runtime)
    snapshot = Base.invokelatest(
        _portable_emt_parent_function(:snapshot_backend_state),
        runtime,
    )
    snapshot isa PortableSnapshotStateInventory || _portable_emt_fail(
        :backend_unavailable,
        "the active backend does not provide portable accepted-state reconstruction",
    )
    return snapshot
end

function _restore_portable_emt_backend_state!(runtime, state::PortableSnapshotStateInventory)
    restored = Base.invokelatest(
        _portable_emt_parent_function(:restore_backend_state!),
        runtime,
        state,
    )
    restored === runtime || _portable_emt_fail(
        :backend_reconstruction,
        "the active backend did not restore the isolated runtime in place",
    )
    return runtime
end

function _portable_emt_accepted_step(context::EMTStepContext)
    return min(max(context.step_index - 1, 0), context.step_count)
end

function _portable_emt_represented_time(context::EMTStepContext)
    value = _portable_emt_accepted_step(context) * context.dt_s
    isfinite(value) && value >= 0.0 || _portable_emt_fail(
        :invalid_represented_time,
        "portable EMT represented time is nonfinite or negative",
    )
    return try
        rationalize(Int128, value; tol = 0)
    catch error
        _portable_emt_fail(
            :time_overflow,
            "portable EMT represented time cannot fit Int128: $(sprint(showerror, error))",
        )
    end
end

function _portable_emt_section(
    snapshot::PortableEMTSnapshot,
    identity::AbstractString,
    visibility::Symbol,
)
    matches = [section for section in snapshot.sections if section.identity == identity]
    length(matches) == 1 || _portable_emt_fail(
        :missing_section,
        "portable EMT snapshot requires exactly one $identity section",
    )
    section = only(matches)
    section.version_major == 1 && section.version_minor == 0 || _portable_emt_fail(
        :unsupported_section_version,
        "portable EMT section $identity has an unsupported version",
    )
    section.visibility == visibility || _portable_emt_fail(
        :section_visibility,
        "portable EMT section $identity has the wrong visibility",
    )
    section.value isa PortableSnapshotRecord || _portable_emt_fail(
        :section_type,
        "portable EMT section $identity is not a registered inventory",
    )
    return section
end

function _portable_emt_validate_metadata(
    metadata::PortableSnapshotMetadata,
    candidate::EMTStudyWorkspace;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    source_descriptor::Union{Nothing,PortableSnapshotDescriptor}=nothing,
)
    metadata.profile == :portable_full || _portable_emt_fail(
        :unsupported_profile,
        "production EMT reconstruction requires the portable_full profile",
    )
    metadata.project_signature_sha256 == project_signature_sha256 || _portable_emt_fail(
        :project_mismatch,
        "portable EMT project signature does not match the receiving project",
    )
    metadata.model_signature_sha256 == model_signature_sha256 || _portable_emt_fail(
        :model_mismatch,
        "portable EMT model signature does not match the receiving model",
    )
    metadata.settings_signature_sha256 == settings_signature_sha256 || _portable_emt_fail(
        :settings_mismatch,
        "portable EMT settings signature does not match the receiving study",
    )
    metadata.topology_signature_sha256 == _emt_checkpoint_topology_fingerprint(candidate) ||
        _portable_emt_fail(
            :topology_mismatch,
            "portable EMT topology signature does not match the receiving study",
        )
    all(capability -> capability in metadata.capabilities, (
        "emt.fixed_step",
        "emt.portable_snapshot",
    )) || _portable_emt_fail(
        :capability_mismatch,
        "portable EMT snapshot lacks a required fixed-step reconstruction capability",
    )
    return metadata
end

"""Capture one completed accepted EMT synchronization point into canonical portable sections."""
function capture_portable_emt_snapshot(
    workspace::EMTStudyWorkspace;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    provenance::AbstractString,
    capabilities::AbstractVector{<:AbstractString}=String[
        "emt.fixed_step",
        "emt.portable_snapshot",
    ],
)
    public_state = portable_emt_state_inventory(workspace)
    backend_state = _portable_emt_backend_snapshot(workspace.runtime.timestep_state)
    context = workspace.runtime.context
    metadata = PortableSnapshotMetadata(
        :portable_full,
        project_signature_sha256,
        model_signature_sha256,
        _emt_checkpoint_topology_fingerprint(workspace),
        settings_signature_sha256,
        _portable_emt_represented_time(context),
        _portable_emt_accepted_step(context),
        capabilities,
        provenance,
    )
    snapshot = PortableEMTSnapshot(
        metadata,
        PortableSnapshotSection[
            PortableSnapshotSection(
                "backend.reconstruction_state",
                1,
                0,
                :private_reconstructible,
                portable_state_inventory_record(backend_state),
            ),
            PortableSnapshotSection(
                "emt.public_state",
                1,
                0,
                :public,
                portable_state_inventory_record(public_state),
            ),
        ],
    )
    portable_snapshot_descriptor(snapshot)
    return snapshot
end

function _restore_portable_emt_snapshot(
    candidate::EMTStudyWorkspace,
    snapshot::PortableEMTSnapshot;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    source_descriptor::Union{Nothing,PortableSnapshotDescriptor}=nothing,
)
    _portable_emt_validate_metadata(
        snapshot.metadata,
        candidate;
        project_signature_sha256,
        model_signature_sha256,
        settings_signature_sha256,
        source_descriptor,
    )
    public_section = _portable_emt_section(snapshot, "emt.public_state", :public)
    backend_section = _portable_emt_section(
        snapshot,
        "backend.reconstruction_state",
        :private_reconstructible,
    )
    public_state = portable_state_inventory(public_section.value)
    backend_state = portable_state_inventory(backend_section.value)
    restore_portable_emt_state_inventory!(candidate, public_state)
    _restore_portable_emt_backend_state!(candidate.runtime.timestep_state, backend_state)
    _portable_emt_accepted_step(candidate.runtime.context) == snapshot.metadata.accepted_step ||
        _portable_emt_fail(
            :accepted_step_mismatch,
            "portable EMT restored accepted step disagrees with its envelope",
        )
    _portable_emt_represented_time(candidate.runtime.context) ==
        snapshot.metadata.represented_time_s || _portable_emt_fail(
            :represented_time_mismatch,
            "portable EMT restored time disagrees with its envelope",
        )
    restored_public = portable_emt_state_inventory(candidate)
    restored_public.signature_sha256 == public_state.signature_sha256 ||
        _portable_emt_fail(
            :public_state_reconstruction,
            "portable EMT public state changed during isolated reconstruction",
        )
    restored_backend = _portable_emt_backend_snapshot(candidate.runtime.timestep_state)
    restored_backend.signature_sha256 == backend_state.signature_sha256 ||
        _portable_emt_fail(
            :backend_state_reconstruction,
            "portable EMT backend state changed during isolated reconstruction",
        )
    return PortableEMTRestoreResult(
        candidate,
        source_descriptor === nothing ? portable_snapshot_descriptor(snapshot) : source_descriptor,
        restored_public.signature_sha256,
        restored_backend.signature_sha256,
        true,
    )
end

"""Restore a portable snapshot into a newly prepared isolated EMT workspace."""
function restore_portable_emt_snapshot(
    prepared::PreparedEMTStudy,
    snapshot::PortableEMTSnapshot;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    source_descriptor::Union{Nothing,PortableSnapshotDescriptor}=nothing,
)
    return _restore_portable_emt_snapshot(
        EMTStudyWorkspace(prepared),
        snapshot;
        project_signature_sha256,
        model_signature_sha256,
        settings_signature_sha256,
        source_descriptor,
    )
end

function write_portable_emt_workspace_snapshot(
    path::AbstractString,
    workspace::EMTStudyWorkspace;
    kwargs...,
)
    snapshot = capture_portable_emt_snapshot(workspace; kwargs...)
    return write_portable_emt_snapshot(path, snapshot)
end

function read_portable_emt_workspace_snapshot(
    path::AbstractString,
    prepared::PreparedEMTStudy;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    maximum_file_bytes::Integer=2_000_000_000,
    maximum_payload_bytes::Integer=maximum_file_bytes,
    maximum_sections::Integer=4096,
    maximum_depth::Integer=64,
)
    snapshot, descriptor = read_portable_emt_snapshot_with_descriptor(
        path;
        allow_private = true,
        maximum_file_bytes,
        maximum_payload_bytes,
        maximum_sections,
        maximum_depth,
    )
    return restore_portable_emt_snapshot(
        prepared,
        snapshot;
        project_signature_sha256,
        model_signature_sha256,
        settings_signature_sha256,
        source_descriptor = descriptor,
    )
end

function _portable_emt_complex_array_fields(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
    values::AbstractArray{<:Complex},
)
    return PortableSnapshotStateField[
        _portable_emt_array_field(
            "$identity.real",
            owner,
            family,
            unit,
            axes,
            Float64.(real.(values)),
        ),
        _portable_emt_array_field(
            "$identity.imaginary",
            owner,
            family,
            unit,
            axes,
            Float64.(imag.(values)),
        ),
    ]
end

function _restore_portable_emt_complex_array!(
    destination::AbstractArray{ComplexF64},
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
)
    real_values = zeros(Float64, size(destination))
    imaginary_values = zeros(Float64, size(destination))
    _restore_portable_emt_array!(
        real_values,
        fields,
        "$identity.real",
        family,
        unit,
        axes,
    )
    _restore_portable_emt_array!(
        imaginary_values,
        fields,
        "$identity.imaginary",
        family,
        unit,
        axes,
    )
    destination .= complex.(real_values, imaginary_values)
    return destination
end

const _PORTABLE_BREQIV_RECORD_FIELDS = (
    (:r, "ohm", :checkpoint),
    (:l, "H", :checkpoint),
    (:c, "F", :checkpoint),
    (:rl, "ohm", :checkpoint),
    (:azr, "S", :reconstruction),
    (:ac1, "ohm", :reconstruction),
    (:azi, "1", :reconstruction),
    (:a2, "1", :reconstruction),
)

function _portable_breqiv_record_fields(
    records,
    identity::AbstractString,
)
    return PortableSnapshotStateField[
        _portable_emt_array_field(
            "$identity.$(String(field_name))",
            "emt.breqiv_branch",
            family,
            unit,
            ["branch"],
            Float64[getfield(record, field_name) for record in records],
        ) for (field_name, unit, family) in _PORTABLE_BREQIV_RECORD_FIELDS
    ]
end

function _restore_portable_breqiv_records!(
    records,
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
)
    for (field_name, unit, family) in _PORTABLE_BREQIV_RECORD_FIELDS
        values = _portable_emt_array_values(
            fields,
            "$identity.$(String(field_name))",
            family,
            unit,
            ["branch"],
            Float64,
        )
        length(values) == length(records) || _portable_emt_fail(
            :state_shape_mismatch,
            "portable BREQIV record field $identity has the wrong branch count",
        )
        for index in eachindex(records, values)
            setfield!(records[index], field_name, values[index])
        end
    end
    return records
end

function _portable_emt_element_state_fields(
    element::BreqivHistoryInjection,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    state = element.state
    fields = _portable_emt_element_identity_field(
        index,
        name,
        "breqiv_history_injection",
    )
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.from_nodes", "emt.breqiv", :checkpoint, "1", ["phase"], element.a),
        _portable_emt_array_field("$prefix.to_nodes", "emt.breqiv", :checkpoint, "1", ["phase"], element.b),
        _portable_emt_array_field("$prefix.initial_phase_voltage", "emt.breqiv", :checkpoint, "pu", ["phase"], element.initial_phase_voltage),
        _portable_emt_array_field("$prefix.phase_voltage", "emt.breqiv", :algebraic, "pu", ["phase"], element.phase_voltage),
        _portable_emt_array_field("$prefix.phase_current", "emt.breqiv", :history, "pu", ["phase"], element.phase_current),
        _portable_emt_array_field("$prefix.phase_admittance", "emt.breqiv", :algebraic, "pu", ["row", "column"], element.phase_admittance),
        _portable_emt_state_field("$prefix.history_current_scale", "emt.breqiv", :checkpoint, "1", element.history_current_scale),
        _portable_emt_state_field("$prefix.history_voltage_scale", "emt.breqiv", :checkpoint, "1", element.history_voltage_scale),
        _portable_emt_state_field("$prefix.frequency_history_voltage_scale", "emt.breqiv", :checkpoint, "1", element.frequency_history_voltage_scale),
        _portable_emt_state_field("$prefix.history_current_consumed_for_step", "emt.breqiv", :history, "1", element.history_current_consumed_for_step),
        _portable_emt_state_field("$prefix.initialized", "emt.breqiv", :discrete, "1", element.initialized),
        _portable_emt_array_field("$prefix.zero_history", "emt.breqiv_history", :history, "pu", ["slot", "branch"], state.zero_history),
        _portable_emt_array_field("$prefix.positive_history", "emt.breqiv_history", :history, "pu", ["slot", "branch", "mode"], state.positive_history),
        _portable_emt_array_field("$prefix.modal_voltage", "emt.breqiv", :algebraic, "pu", ["mode"], state.modal),
        _portable_emt_array_field("$prefix.modal_current", "emt.breqiv", :history, "pu", ["mode"], state.modal_current),
        _portable_emt_state_field("$prefix.start_isfd", "emt.breqiv", :checkpoint, "1", state.start_isfd),
        _portable_emt_state_field("$prefix.start_ibf", "emt.breqiv", :checkpoint, "1", state.start_ibf),
        _portable_emt_state_field("$prefix.ikf", "emt.breqiv", :scheduler, "1", state.ikf),
        _portable_emt_state_field("$prefix.isfd", "emt.breqiv", :scheduler, "1", state.isfd),
        _portable_emt_state_field("$prefix.ibf", "emt.breqiv", :scheduler, "1", state.ibf),
    ])
    append!(fields, _portable_breqiv_record_fields(state.zero, "$prefix.zero_branch"))
    append!(fields, _portable_breqiv_record_fields(state.positive, "$prefix.positive_branch"))
    return fields
end

function _restore_portable_emt_element_state!(
    element::BreqivHistoryInjection,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    state = element.state
    _restore_portable_emt_array!(element.initial_phase_voltage, fields, "$prefix.initial_phase_voltage", :checkpoint, "pu", ["phase"])
    _restore_portable_emt_array!(element.phase_voltage, fields, "$prefix.phase_voltage", :algebraic, "pu", ["phase"])
    _restore_portable_emt_array!(element.phase_current, fields, "$prefix.phase_current", :history, "pu", ["phase"])
    _restore_portable_emt_array!(element.phase_admittance, fields, "$prefix.phase_admittance", :algebraic, "pu", ["row", "column"])
    element.history_current_scale = _portable_emt_scalar(fields, "$prefix.history_current_scale", :checkpoint, "1", Float64)
    element.history_voltage_scale = _portable_emt_scalar(fields, "$prefix.history_voltage_scale", :checkpoint, "1", Float64)
    element.frequency_history_voltage_scale = _portable_emt_scalar(fields, "$prefix.frequency_history_voltage_scale", :checkpoint, "1", Float64)
    element.history_current_consumed_for_step = _portable_emt_scalar(fields, "$prefix.history_current_consumed_for_step", :history, "1", Bool)
    element.initialized = _portable_emt_scalar(fields, "$prefix.initialized", :discrete, "1", Bool)
    _restore_portable_emt_array!(state.zero_history, fields, "$prefix.zero_history", :history, "pu", ["slot", "branch"])
    _restore_portable_emt_array!(state.positive_history, fields, "$prefix.positive_history", :history, "pu", ["slot", "branch", "mode"])
    _restore_portable_emt_array!(state.modal, fields, "$prefix.modal_voltage", :algebraic, "pu", ["mode"])
    _restore_portable_emt_array!(state.modal_current, fields, "$prefix.modal_current", :history, "pu", ["mode"])
    state.start_isfd = _portable_emt_integer(fields, "$prefix.start_isfd", :checkpoint)
    state.start_ibf = _portable_emt_integer(fields, "$prefix.start_ibf", :checkpoint)
    state.ikf = _portable_emt_integer(fields, "$prefix.ikf", :scheduler)
    state.isfd = _portable_emt_integer(fields, "$prefix.isfd", :scheduler)
    state.ibf = _portable_emt_integer(fields, "$prefix.ibf", :scheduler)
    _restore_portable_breqiv_records!(state.zero, fields, "$prefix.zero_branch")
    _restore_portable_breqiv_records!(state.positive, fields, "$prefix.positive_branch")
    return element
end

function _portable_emt_element_state_fields(
    line::ComplexModalBergeronLine,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "complex_modal_bergeron_line")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.from_nodes", "emt.complex_modal_line", :checkpoint, "1", ["phase"], line.from_nodes),
        _portable_emt_array_field("$prefix.to_nodes", "emt.complex_modal_line", :checkpoint, "1", ["phase"], line.to_nodes),
        _portable_emt_array_field("$prefix.travel_time", "emt.complex_modal_line", :checkpoint, "s", ["mode"], line.travel_time_s),
        _portable_emt_array_field("$prefix.attenuation", "emt.complex_modal_line", :checkpoint, "1", ["mode"], line.attenuation),
        _portable_emt_array_field("$prefix.delay_steps", "emt.complex_modal_line", :checkpoint, "1", ["mode"], line.delay_steps),
        _portable_emt_array_field("$prefix.delay_interpolation_factor", "emt.complex_modal_line", :checkpoint, "1", ["mode"], line.delay_interpolation_factors),
        _portable_emt_array_field("$prefix.write_indices", "emt.complex_modal_line", :delayed, "1", ["mode"], line.write_indices),
        _portable_emt_array_field("$prefix.phase_admittance", "emt.complex_modal_line", :checkpoint, "S", ["row", "column"], line.phase_admittance),
        _portable_emt_array_field("$prefix.phase_current_from", "emt.complex_modal_line", :algebraic, "A", ["phase"], line.phase_current_from),
        _portable_emt_array_field("$prefix.phase_current_to", "emt.complex_modal_line", :algebraic, "A", ["phase"], line.phase_current_to),
        _portable_emt_array_field("$prefix.phase_history_from", "emt.complex_modal_line", :history, "A", ["phase"], line.phase_history_from),
        _portable_emt_array_field("$prefix.phase_history_to", "emt.complex_modal_line", :history, "A", ["phase"], line.phase_history_to),
        _portable_emt_state_field("$prefix.power_invariance_error", "emt.complex_modal_line", :checkpoint, "1", line.power_invariance_error),
        _portable_emt_state_field("$prefix.imaginary_stamp_residual", "emt.complex_modal_line", :checkpoint, "S", line.imaginary_stamp_residual),
        _portable_emt_state_field("$prefix.timestep", "emt.complex_modal_line", :checkpoint, "s", line.dt_s),
    ])
    for (identity, values, unit, axes, family) in (
        ("$prefix.phase_to_modal", line.transform.phase_to_modal, "1", ["row", "column"], :checkpoint),
        ("$prefix.modal_to_phase", line.transform.modal_to_phase, "1", ["row", "column"], :checkpoint),
        ("$prefix.characteristic_admittance", line.characteristic_admittance, "S", ["mode"], :checkpoint),
        ("$prefix.phase_admittance_complex", line.phase_admittance_complex, "S", ["row", "column"], :checkpoint),
        ("$prefix.modal_voltage_from", line.modal_voltage_from, "V", ["mode"], :algebraic),
        ("$prefix.modal_voltage_to", line.modal_voltage_to, "V", ["mode"], :algebraic),
        ("$prefix.modal_current_from", line.modal_current_from, "A", ["mode"], :algebraic),
        ("$prefix.modal_current_to", line.modal_current_to, "A", ["mode"], :algebraic),
        ("$prefix.modal_history_from", line.modal_history_from, "A", ["mode"], :history),
        ("$prefix.modal_history_to", line.modal_history_to, "A", ["mode"], :history),
    )
        append!(fields, _portable_emt_complex_array_fields(identity, "emt.complex_modal_line", family, unit, axes, values))
    end
    for mode in eachindex(line.from_wave_history)
        mode_prefix = "$prefix.mode.i" * lpad(string(mode), 8, '0')
        append!(fields, _portable_emt_complex_array_fields("$mode_prefix.from_wave_history", "emt.complex_modal_line", :delayed, "V", ["delay"], line.from_wave_history[mode]))
        append!(fields, _portable_emt_complex_array_fields("$mode_prefix.to_wave_history", "emt.complex_modal_line", :delayed, "V", ["delay"], line.to_wave_history[mode]))
    end
    return fields
end

function _restore_portable_emt_element_state!(
    line::ComplexModalBergeronLine,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _restore_portable_emt_array!(line.write_indices, fields, "$prefix.write_indices", :delayed, "1", ["mode"])
    for (identity, destination, unit, axes, family) in (
        ("$prefix.modal_voltage_from", line.modal_voltage_from, "V", ["mode"], :algebraic),
        ("$prefix.modal_voltage_to", line.modal_voltage_to, "V", ["mode"], :algebraic),
        ("$prefix.modal_current_from", line.modal_current_from, "A", ["mode"], :algebraic),
        ("$prefix.modal_current_to", line.modal_current_to, "A", ["mode"], :algebraic),
        ("$prefix.modal_history_from", line.modal_history_from, "A", ["mode"], :history),
        ("$prefix.modal_history_to", line.modal_history_to, "A", ["mode"], :history),
    )
        _restore_portable_emt_complex_array!(destination, fields, identity, family, unit, axes)
    end
    for mode in eachindex(line.from_wave_history)
        mode_prefix = "$prefix.mode.i" * lpad(string(mode), 8, '0')
        _restore_portable_emt_complex_array!(line.from_wave_history[mode], fields, "$mode_prefix.from_wave_history", :delayed, "V", ["delay"])
        _restore_portable_emt_complex_array!(line.to_wave_history[mode], fields, "$mode_prefix.to_wave_history", :delayed, "V", ["delay"])
    end
    _restore_portable_emt_array!(line.phase_current_from, fields, "$prefix.phase_current_from", :algebraic, "A", ["phase"])
    _restore_portable_emt_array!(line.phase_current_to, fields, "$prefix.phase_current_to", :algebraic, "A", ["phase"])
    _restore_portable_emt_array!(line.phase_history_from, fields, "$prefix.phase_history_from", :history, "A", ["phase"])
    _restore_portable_emt_array!(line.phase_history_to, fields, "$prefix.phase_history_to", :history, "A", ["phase"])
    return line
end

function _portable_emt_element_state_fields(
    line::CoupledFrequencyDependentLine,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    state = line.runtime_state
    preparation = state.preparation
    fields = _portable_emt_element_identity_field(index, name, "coupled_frequency_dependent_line")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.from_nodes", "emt.coupled_line", :checkpoint, "1", ["phase"], line.from_nodes),
        _portable_emt_array_field("$prefix.to_nodes", "emt.coupled_line", :checkpoint, "1", ["phase"], line.to_nodes),
        _portable_emt_state_field("$prefix.preparation_signature", "emt.coupled_line", :checkpoint, "1", preparation.deterministic_signature_sha256),
        _portable_emt_array_field("$prefix.rational_state", "emt.coupled_line", :continuous, "1", ["state"], state.rational_state),
        _portable_emt_array_field("$prefix.previous_incident_wave", "emt.coupled_line", :history, "V", ["port"], state.previous_incident_wave),
        _portable_emt_array_field("$prefix.incident_wave", "emt.coupled_line", :algebraic, "V", ["port"], state.incident_wave),
        _portable_emt_array_field("$prefix.outgoing_wave", "emt.coupled_line", :algebraic, "V", ["port"], state.outgoing_wave),
        _portable_emt_array_field("$prefix.terminal_voltage", "emt.coupled_line", :algebraic, "V", ["port"], state.terminal_voltage_v),
        _portable_emt_array_field("$prefix.terminal_current", "emt.coupled_line", :algebraic, "A", ["port"], state.terminal_current_a),
        _portable_emt_array_field("$prefix.history_current", "emt.coupled_line", :history, "A", ["port"], state.history_current_a),
        _portable_emt_array_field("$prefix.voltage_workspace", "emt.coupled_line", :reconstruction, "V", ["port"], line.voltage_workspace_v),
        _portable_emt_state_field("$prefix.previous_terminal_power", "emt.coupled_line", :history, "W", state.previous_terminal_power_w),
        _portable_emt_state_field("$prefix.terminal_power", "emt.coupled_line", :output, "W", state.terminal_power_w),
        _portable_emt_state_field("$prefix.cumulative_supplied_energy", "emt.coupled_line", :history, "J", state.cumulative_supplied_energy_j),
        _portable_emt_state_field("$prefix.minimum_cumulative_supplied_energy", "emt.coupled_line", :output, "J", state.minimum_cumulative_supplied_energy_j),
        _portable_emt_state_field("$prefix.maximum_kcl_residual", "emt.coupled_line", :output, "A", state.maximum_kcl_residual_a),
        _portable_emt_state_field("$prefix.maximum_state_magnitude", "emt.coupled_line", :output, "1", state.maximum_state_magnitude),
        _portable_emt_state_field("$prefix.accepted_time", "emt.coupled_line", :scheduler, "s", state.accepted_time_s),
        _portable_emt_state_field("$prefix.accepted_step_count", "emt.coupled_line", :scheduler, "1", state.accepted_step_count),
        _portable_emt_state_field("$prefix.initialization_kind", "emt.coupled_line", :discrete, "1", String(state.initialization_kind)),
        _portable_emt_state_field("$prefix.initialization_frequency", "emt.coupled_line", :checkpoint, "Hz", state.initialization_frequency_hz),
    ])
    return fields
end

function _restore_portable_emt_element_state!(
    line::CoupledFrequencyDependentLine,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    state = line.runtime_state
    _portable_emt_scalar(fields, "$prefix.preparation_signature", :checkpoint, "1", String) == state.preparation.deterministic_signature_sha256 || _portable_emt_fail(
        :model_mismatch,
        "portable coupled-line preparation changed",
    )
    _restore_portable_emt_array!(state.rational_state, fields, "$prefix.rational_state", :continuous, "1", ["state"])
    _restore_portable_emt_array!(state.previous_incident_wave, fields, "$prefix.previous_incident_wave", :history, "V", ["port"])
    _restore_portable_emt_array!(state.incident_wave, fields, "$prefix.incident_wave", :algebraic, "V", ["port"])
    _restore_portable_emt_array!(state.outgoing_wave, fields, "$prefix.outgoing_wave", :algebraic, "V", ["port"])
    _restore_portable_emt_array!(state.terminal_voltage_v, fields, "$prefix.terminal_voltage", :algebraic, "V", ["port"])
    _restore_portable_emt_array!(state.terminal_current_a, fields, "$prefix.terminal_current", :algebraic, "A", ["port"])
    _restore_portable_emt_array!(state.history_current_a, fields, "$prefix.history_current", :history, "A", ["port"])
    _restore_portable_emt_array!(line.voltage_workspace_v, fields, "$prefix.voltage_workspace", :reconstruction, "V", ["port"])
    state.previous_terminal_power_w = _portable_emt_scalar(fields, "$prefix.previous_terminal_power", :history, "W", Float64)
    state.terminal_power_w = _portable_emt_scalar(fields, "$prefix.terminal_power", :output, "W", Float64)
    state.cumulative_supplied_energy_j = _portable_emt_scalar(fields, "$prefix.cumulative_supplied_energy", :history, "J", Float64)
    state.minimum_cumulative_supplied_energy_j = _portable_emt_scalar(fields, "$prefix.minimum_cumulative_supplied_energy", :output, "J", Float64)
    state.maximum_kcl_residual_a = _portable_emt_scalar(fields, "$prefix.maximum_kcl_residual", :output, "A", Float64)
    state.maximum_state_magnitude = _portable_emt_scalar(fields, "$prefix.maximum_state_magnitude", :output, "1", Float64)
    state.accepted_time_s = _portable_emt_scalar(fields, "$prefix.accepted_time", :scheduler, "s", Float64)
    state.accepted_step_count = _portable_emt_nonnegative_integer(fields, "$prefix.accepted_step_count", :scheduler)
    state.initialization_kind = Symbol(_portable_emt_scalar(fields, "$prefix.initialization_kind", :discrete, "1", String))
    state.initialization_frequency_hz = _portable_emt_scalar(fields, "$prefix.initialization_frequency", :checkpoint, "Hz", Float64)
    fill!(state.state_workspace, 0.0)
    fill!(state.wave_history_workspace, 0.0)
    fill!(state.current_workspace, 0.0)
    fill!(state.incident_workspace, 0.0)
    fill!(state.outgoing_workspace, 0.0)
    return line
end

function _portable_sampled_line_state_fields(
    line::SampledFrequencyDependentLine,
    prefix::AbstractString,
)
    coefficients = line.coefficients
    return PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.timestep", "emt.sampled_line", :checkpoint, "s", coefficients.timestep_s),
        _portable_emt_state_field("$prefix.characteristic_impedance", "emt.sampled_line", :checkpoint, "ohm", coefficients.characteristic_impedance_ohm),
        _portable_emt_state_field("$prefix.eta", "emt.sampled_line", :checkpoint, "S", coefficients.eta),
        _portable_emt_state_field("$prefix.propagation_delay_steps", "emt.sampled_line", :checkpoint, "1", coefficients.propagation_delay_steps),
        _portable_emt_state_field("$prefix.history_sample_count", "emt.sampled_line", :checkpoint, "1", coefficients.history_sample_count),
        _portable_emt_array_field("$prefix.admittance_weights", "emt.sampled_line", :checkpoint, "1", ["weight"], coefficients.admittance_weights),
        _portable_emt_array_field("$prefix.propagation_weights", "emt.sampled_line", :checkpoint, "1", ["weight"], coefficients.propagation_weights),
        _portable_emt_array_field("$prefix.admittance_tail", "emt.sampled_line", :checkpoint, "1", ["coefficient"], collect(coefficients.admittance_tail)),
        _portable_emt_array_field("$prefix.propagation_tail", "emt.sampled_line", :checkpoint, "1", ["coefficient"], collect(coefficients.propagation_tail)),
        _portable_emt_state_field("$prefix.loss_factor", "emt.sampled_line", :checkpoint, "1", line.loss_factor),
        _portable_emt_array_field("$prefix.from_wave_history", "emt.sampled_line", :delayed, "V", ["delay"], line.from_wave_history),
        _portable_emt_array_field("$prefix.to_wave_history", "emt.sampled_line", :delayed, "V", ["delay"], line.to_wave_history),
        _portable_emt_state_field("$prefix.write_index", "emt.sampled_line", :delayed, "1", line.write_index),
        _portable_emt_state_field("$prefix.admittance_tail_from", "emt.sampled_line", :history, "V", line.admittance_tail_from),
        _portable_emt_state_field("$prefix.admittance_tail_to", "emt.sampled_line", :history, "V", line.admittance_tail_to),
        _portable_emt_state_field("$prefix.propagation_tail_from", "emt.sampled_line", :history, "V", line.propagation_tail_from),
        _portable_emt_state_field("$prefix.propagation_tail_to", "emt.sampled_line", :history, "V", line.propagation_tail_to),
        _portable_emt_state_field("$prefix.history_current_from", "emt.sampled_line", :history, "A", line.history_current_from),
        _portable_emt_state_field("$prefix.history_current_to", "emt.sampled_line", :history, "A", line.history_current_to),
        _portable_emt_state_field("$prefix.terminal_voltage_from", "emt.sampled_line", :algebraic, "V", line.terminal_voltage_from),
        _portable_emt_state_field("$prefix.terminal_voltage_to", "emt.sampled_line", :algebraic, "V", line.terminal_voltage_to),
        _portable_emt_state_field("$prefix.terminal_current_from", "emt.sampled_line", :algebraic, "A", line.terminal_current_from),
        _portable_emt_state_field("$prefix.terminal_current_to", "emt.sampled_line", :algebraic, "A", line.terminal_current_to),
        _portable_emt_state_field("$prefix.convolution_from", "emt.sampled_line", :history, "V", line.convolution_from),
        _portable_emt_state_field("$prefix.convolution_to", "emt.sampled_line", :history, "V", line.convolution_to),
        _portable_emt_state_field("$prefix.update_count", "emt.sampled_line", :scheduler, "1", line.update_count),
    ]
end

function _restore_portable_sampled_line_state!(
    line::SampledFrequencyDependentLine,
    fields::Dict{String,PortableSnapshotStateField},
    prefix::AbstractString,
)
    _restore_portable_emt_array!(line.from_wave_history, fields, "$prefix.from_wave_history", :delayed, "V", ["delay"])
    _restore_portable_emt_array!(line.to_wave_history, fields, "$prefix.to_wave_history", :delayed, "V", ["delay"])
    line.write_index = _portable_emt_integer(fields, "$prefix.write_index", :delayed)
    for (field_name, family, unit) in (
        (:admittance_tail_from, :history, "V"),
        (:admittance_tail_to, :history, "V"),
        (:propagation_tail_from, :history, "V"),
        (:propagation_tail_to, :history, "V"),
        (:history_current_from, :history, "A"),
        (:history_current_to, :history, "A"),
        (:terminal_voltage_from, :algebraic, "V"),
        (:terminal_voltage_to, :algebraic, "V"),
        (:terminal_current_from, :algebraic, "A"),
        (:terminal_current_to, :algebraic, "A"),
        (:convolution_from, :history, "V"),
        (:convolution_to, :history, "V"),
    )
        setfield!(line, field_name, _portable_emt_scalar(fields, "$prefix.$(String(field_name))", family, unit, Float64))
    end
    line.update_count = _portable_emt_nonnegative_integer(fields, "$prefix.update_count", :scheduler)
    return line
end

function _portable_emt_element_state_fields(
    line::SampledFrequencyDependentLine,
    index::Integer,
    name::Symbol,
)
    fields = _portable_emt_element_identity_field(index, name, "sampled_frequency_dependent_line")
    append!(fields, _portable_sampled_line_state_fields(line, _portable_emt_element_prefix(index)))
    return fields
end

function _restore_portable_emt_element_state!(
    line::SampledFrequencyDependentLine,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    return _restore_portable_sampled_line_state!(line, fields, _portable_emt_element_prefix(index))
end

function _portable_emt_element_state_fields(
    line::SampledFrequencyDependentLineGroup,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "sampled_frequency_dependent_line_group")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.modal_line_count", "emt.sampled_line_group", :checkpoint, "1", length(line.modal_lines)),
        _portable_emt_array_field("$prefix.phase_to_modal", "emt.sampled_line_group", :checkpoint, "1", ["row", "column"], line.phase_to_modal),
        _portable_emt_array_field("$prefix.modal_to_phase", "emt.sampled_line_group", :checkpoint, "1", ["row", "column"], line.modal_to_phase),
        _portable_emt_array_field("$prefix.phase_admittance", "emt.sampled_line_group", :checkpoint, "S", ["row", "column"], line.phase_admittance),
        _portable_emt_array_field("$prefix.history_current_from", "emt.sampled_line_group", :history, "A", ["phase"], line.history_current_from),
        _portable_emt_array_field("$prefix.history_current_to", "emt.sampled_line_group", :history, "A", ["phase"], line.history_current_to),
        _portable_emt_array_field("$prefix.terminal_voltage_from", "emt.sampled_line_group", :algebraic, "V", ["phase"], line.terminal_voltage_from),
        _portable_emt_array_field("$prefix.terminal_voltage_to", "emt.sampled_line_group", :algebraic, "V", ["phase"], line.terminal_voltage_to),
        _portable_emt_array_field("$prefix.terminal_current_from", "emt.sampled_line_group", :algebraic, "A", ["phase"], line.terminal_current_from),
        _portable_emt_array_field("$prefix.terminal_current_to", "emt.sampled_line_group", :algebraic, "A", ["phase"], line.terminal_current_to),
        _portable_emt_array_field("$prefix.modal_voltage_from", "emt.sampled_line_group", :algebraic, "V", ["mode"], line.modal_voltage_from),
        _portable_emt_array_field("$prefix.modal_voltage_to", "emt.sampled_line_group", :algebraic, "V", ["mode"], line.modal_voltage_to),
        _portable_emt_array_field("$prefix.modal_current_from", "emt.sampled_line_group", :algebraic, "A", ["mode"], line.modal_current_from),
        _portable_emt_array_field("$prefix.modal_current_to", "emt.sampled_line_group", :algebraic, "A", ["mode"], line.modal_current_to),
        _portable_emt_state_field("$prefix.update_count", "emt.sampled_line_group", :scheduler, "1", line.update_count),
    ])
    for mode in eachindex(line.modal_lines)
        mode_prefix = "$prefix.modal.i" * lpad(string(mode), 8, '0')
        append!(fields, _portable_sampled_line_state_fields(line.modal_lines[mode], mode_prefix))
    end
    return fields
end

function _restore_portable_emt_element_state!(
    line::SampledFrequencyDependentLineGroup,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _portable_emt_integer(fields, "$prefix.modal_line_count", :checkpoint) == length(line.modal_lines) || _portable_emt_fail(
        :model_mismatch,
        "portable sampled-line modal count changed",
    )
    for mode in eachindex(line.modal_lines)
        mode_prefix = "$prefix.modal.i" * lpad(string(mode), 8, '0')
        _restore_portable_sampled_line_state!(line.modal_lines[mode], fields, mode_prefix)
    end
    for (field_name, family, unit, axis) in (
        (:history_current_from, :history, "A", "phase"),
        (:history_current_to, :history, "A", "phase"),
        (:terminal_voltage_from, :algebraic, "V", "phase"),
        (:terminal_voltage_to, :algebraic, "V", "phase"),
        (:terminal_current_from, :algebraic, "A", "phase"),
        (:terminal_current_to, :algebraic, "A", "phase"),
        (:modal_voltage_from, :algebraic, "V", "mode"),
        (:modal_voltage_to, :algebraic, "V", "mode"),
        (:modal_current_from, :algebraic, "A", "mode"),
        (:modal_current_to, :algebraic, "A", "mode"),
    )
        _restore_portable_emt_array!(getfield(line, field_name), fields, "$prefix.$(String(field_name))", family, unit, [axis])
    end
    line.update_count = _portable_emt_nonnegative_integer(fields, "$prefix.update_count", :scheduler)
    return line
end

function _portable_semlyen_mode_state_fields(
    state,
    prefix::AbstractString,
)
    parameters = state.parameters
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.characteristic_admittance", "emt.semlyen_line", :checkpoint, "S", parameters.characteristic_admittance_s),
        _portable_emt_state_field("$prefix.travel_time", "emt.semlyen_line", :checkpoint, "s", parameters.travel_time_s),
        _portable_emt_state_field("$prefix.phasor_frequency", "emt.semlyen_line", :checkpoint, "Hz", parameters.phasor_frequency_hz),
        _portable_emt_state_field("$prefix.propagation_term_count", "emt.semlyen_line", :checkpoint, "1", length(parameters.propagation_terms)),
        _portable_emt_state_field("$prefix.admittance_term_count", "emt.semlyen_line", :checkpoint, "1", length(parameters.admittance_terms)),
        _portable_emt_array_field("$prefix.outgoing_from_history", "emt.semlyen_line", :delayed, "V", ["delay"], state.outgoing_from_history),
        _portable_emt_array_field("$prefix.outgoing_to_history", "emt.semlyen_line", :delayed, "V", ["delay"], state.outgoing_to_history),
        _portable_emt_state_field("$prefix.write_index", "emt.semlyen_line", :delayed, "1", state.write_index),
        _portable_emt_state_field("$prefix.delay_steps", "emt.semlyen_line", :checkpoint, "1", state.delay_steps),
        _portable_emt_state_field("$prefix.delay_fraction", "emt.semlyen_line", :checkpoint, "1", state.delay_fraction),
        _portable_emt_state_field("$prefix.history_current_from", "emt.semlyen_line", :history, "A", state.history_current_from),
        _portable_emt_state_field("$prefix.history_current_to", "emt.semlyen_line", :history, "A", state.history_current_to),
    ]
    for (identity, values, unit) in (
        ("$prefix.propagation_from_state", state.propagation_from_state, "V"),
        ("$prefix.propagation_to_state", state.propagation_to_state, "V"),
        ("$prefix.admittance_from_state", state.admittance_from_state, "A"),
        ("$prefix.admittance_to_state", state.admittance_to_state, "A"),
    )
        append!(fields, _portable_emt_complex_array_fields(identity, "emt.semlyen_line", :history, unit, ["term"], values))
    end
    return fields
end

function _restore_portable_semlyen_mode_state!(
    state,
    fields::Dict{String,PortableSnapshotStateField},
    prefix::AbstractString,
)
    for (identity, destination, unit) in (
        ("$prefix.propagation_from_state", state.propagation_from_state, "V"),
        ("$prefix.propagation_to_state", state.propagation_to_state, "V"),
        ("$prefix.admittance_from_state", state.admittance_from_state, "A"),
        ("$prefix.admittance_to_state", state.admittance_to_state, "A"),
    )
        _restore_portable_emt_complex_array!(destination, fields, identity, :history, unit, ["term"])
    end
    _restore_portable_emt_array!(state.outgoing_from_history, fields, "$prefix.outgoing_from_history", :delayed, "V", ["delay"])
    _restore_portable_emt_array!(state.outgoing_to_history, fields, "$prefix.outgoing_to_history", :delayed, "V", ["delay"])
    state.write_index = _portable_emt_integer(fields, "$prefix.write_index", :delayed)
    state.history_current_from = _portable_emt_scalar(fields, "$prefix.history_current_from", :history, "A", Float64)
    state.history_current_to = _portable_emt_scalar(fields, "$prefix.history_current_to", :history, "A", Float64)
    return state
end

function _portable_emt_element_state_fields(
    line::SemlyenFrequencyDependentLine,
    index::Integer,
    name::Symbol,
)
    prefix = _portable_emt_element_prefix(index)
    fields = _portable_emt_element_identity_field(index, name, "semlyen_frequency_dependent_line")
    append!(fields, PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.mode_count", "emt.semlyen_line", :checkpoint, "1", length(line.modes)),
        _portable_emt_state_field("$prefix.timestep", "emt.semlyen_line", :checkpoint, "s", line.timestep_s),
        _portable_emt_array_field("$prefix.runtime_current_modal_to_phase", "emt.semlyen_line", :checkpoint, "1", ["row", "column"], line.runtime_current_modal_to_phase),
        _portable_emt_array_field("$prefix.phase_to_modal_voltage", "emt.semlyen_line", :checkpoint, "1", ["row", "column"], line.phase_to_modal_voltage),
        _portable_emt_array_field("$prefix.phase_admittance", "emt.semlyen_line", :checkpoint, "S", ["row", "column"], line.phase_admittance),
        _portable_emt_array_field("$prefix.history_current_from", "emt.semlyen_line", :history, "A", ["phase"], line.history_current_from),
        _portable_emt_array_field("$prefix.history_current_to", "emt.semlyen_line", :history, "A", ["phase"], line.history_current_to),
        _portable_emt_array_field("$prefix.terminal_voltage_from", "emt.semlyen_line", :algebraic, "V", ["phase"], line.terminal_voltage_from),
        _portable_emt_array_field("$prefix.terminal_voltage_to", "emt.semlyen_line", :algebraic, "V", ["phase"], line.terminal_voltage_to),
        _portable_emt_array_field("$prefix.terminal_current_from", "emt.semlyen_line", :algebraic, "A", ["phase"], line.terminal_current_from),
        _portable_emt_array_field("$prefix.terminal_current_to", "emt.semlyen_line", :algebraic, "A", ["phase"], line.terminal_current_to),
        _portable_emt_array_field("$prefix.modal_voltage_from", "emt.semlyen_line", :algebraic, "V", ["mode"], line.modal_voltage_from),
        _portable_emt_array_field("$prefix.modal_voltage_to", "emt.semlyen_line", :algebraic, "V", ["mode"], line.modal_voltage_to),
        _portable_emt_array_field("$prefix.modal_current_from", "emt.semlyen_line", :algebraic, "A", ["mode"], line.modal_current_from),
        _portable_emt_array_field("$prefix.modal_current_to", "emt.semlyen_line", :algebraic, "A", ["mode"], line.modal_current_to),
        _portable_emt_array_field("$prefix.modal_history_current_from", "emt.semlyen_line", :history, "A", ["mode"], line.modal_history_current_from),
        _portable_emt_array_field("$prefix.modal_history_current_to", "emt.semlyen_line", :history, "A", ["mode"], line.modal_history_current_to),
        _portable_emt_state_field("$prefix.update_count", "emt.semlyen_line", :scheduler, "1", line.update_count),
    ])
    append!(fields, _portable_emt_complex_array_fields("$prefix.voltage_modal_to_phase", "emt.semlyen_line", :checkpoint, "1", ["row", "column"], line.voltage_modal_to_phase))
    append!(fields, _portable_emt_complex_array_fields("$prefix.current_modal_to_phase", "emt.semlyen_line", :checkpoint, "1", ["row", "column"], line.current_modal_to_phase))
    for mode in eachindex(line.modes)
        mode_prefix = "$prefix.mode.i" * lpad(string(mode), 8, '0')
        append!(fields, _portable_semlyen_mode_state_fields(line.modes[mode], mode_prefix))
    end
    return fields
end

function _restore_portable_emt_element_state!(
    line::SemlyenFrequencyDependentLine,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    prefix = _portable_emt_element_prefix(index)
    _portable_emt_integer(fields, "$prefix.mode_count", :checkpoint) == length(line.modes) || _portable_emt_fail(
        :model_mismatch,
        "portable Semlyen-line mode count changed",
    )
    for mode in eachindex(line.modes)
        mode_prefix = "$prefix.mode.i" * lpad(string(mode), 8, '0')
        _restore_portable_semlyen_mode_state!(line.modes[mode], fields, mode_prefix)
    end
    for (field_name, family, unit, axis) in (
        (:history_current_from, :history, "A", "phase"),
        (:history_current_to, :history, "A", "phase"),
        (:terminal_voltage_from, :algebraic, "V", "phase"),
        (:terminal_voltage_to, :algebraic, "V", "phase"),
        (:terminal_current_from, :algebraic, "A", "phase"),
        (:terminal_current_to, :algebraic, "A", "phase"),
        (:modal_voltage_from, :algebraic, "V", "mode"),
        (:modal_voltage_to, :algebraic, "V", "mode"),
        (:modal_current_from, :algebraic, "A", "mode"),
        (:modal_current_to, :algebraic, "A", "mode"),
        (:modal_history_current_from, :history, "A", "mode"),
        (:modal_history_current_to, :history, "A", "mode"),
    )
        _restore_portable_emt_array!(getfield(line, field_name), fields, "$prefix.$(String(field_name))", family, unit, [axis])
    end
    line.update_count = _portable_emt_nonnegative_integer(fields, "$prefix.update_count", :scheduler)
    return line
end

function _portable_semiconductor_provenance_fields(
    provenance,
    prefix::AbstractString,
    owner::AbstractString,
)
    return PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.source", owner, :checkpoint, "1", provenance.source),
        _portable_emt_state_field("$prefix.units", owner, :checkpoint, "1", provenance.units),
        _portable_emt_state_field("$prefix.transformation", owner, :checkpoint, "1", provenance.transformation),
        _portable_emt_state_field("$prefix.uncertainty", owner, :checkpoint, "1", provenance.uncertainty),
        _portable_emt_state_field("$prefix.validity_domain", owner, :checkpoint, "1", provenance.validity_domain),
        _portable_emt_state_field("$prefix.nature", owner, :checkpoint, "1", string(provenance.nature)),
    ]
end

function _portable_semiconductor_optional_presence_field(
    prefix::AbstractString,
    component::AbstractString,
    value,
)
    return _portable_emt_state_field(
        "$prefix.$component.present",
        "emt.power_semiconductor",
        :checkpoint,
        "1",
        value !== nothing,
    )
end

function _portable_semiconductor_require_presence(
    fields::Dict{String,PortableSnapshotStateField},
    prefix::AbstractString,
    component::AbstractString,
    value,
)
    captured = _portable_emt_scalar(
        fields,
        "$prefix.$component.present",
        :checkpoint,
        "1",
        Bool,
    )
    captured == (value !== nothing) || _portable_emt_fail(
        :model_mismatch,
        "portable power-semiconductor component $prefix.$component changed",
    )
    return value
end

function _portable_semiconductor_quantity_field(
    identity::AbstractString,
    owner::AbstractString,
    family::Symbol,
    unit::AbstractString,
    value::Real;
    allow_nan::Bool=false,
)
    quantity = Float64(value)
    if isnan(quantity)
        allow_nan || _portable_emt_fail(
            :nonfinite_value,
            "portable power-semiconductor field $identity is NaN",
        )
        return _portable_emt_state_field(
            identity,
            owner,
            family,
            unit,
            "not_initialized",
        )
    end
    encoded = isfinite(quantity) ? quantity :
        quantity == Inf ? _PORTABLE_EMT_POSITIVE_INFINITY :
        quantity == -Inf ? _PORTABLE_EMT_NEGATIVE_INFINITY :
        _portable_emt_fail(
            :nonfinite_value,
            "portable power-semiconductor field $identity is unsupported",
        )
    return _portable_emt_state_field(identity, owner, family, unit, encoded)
end

function _portable_semiconductor_quantity_scalar(
    fields::Dict{String,PortableSnapshotStateField},
    identity::AbstractString,
    family::Symbol,
    unit::AbstractString;
    allow_nan::Bool=false,
)
    value = _portable_emt_inventory_field(fields, identity, family, unit).value
    value isa Float64 && return value
    value == _PORTABLE_EMT_POSITIVE_INFINITY && return Inf
    value == _PORTABLE_EMT_NEGATIVE_INFINITY && return -Inf
    allow_nan && value == "not_initialized" && return NaN
    _portable_emt_fail(
        :state_type_mismatch,
        "portable power-semiconductor field $identity has the wrong scalar type",
    )
end

function _portable_semiconductor_gate_driver_fields(driver, prefix::AbstractString)
    owner = "emt.power_semiconductor_gate"
    return PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.turn_on_delay", owner, :checkpoint, "s", driver.turn_on_delay_s),
        _portable_emt_state_field("$prefix.turn_off_delay", owner, :checkpoint, "s", driver.turn_off_delay_s),
        _portable_emt_state_field("$prefix.dead_time", owner, :checkpoint, "s", driver.dead_time_s),
        _portable_emt_state_field("$prefix.minimum_pulse_width", owner, :checkpoint, "s", driver.minimum_pulse_width_s),
        _portable_emt_state_field("$prefix.commanded_on", owner, :discrete, "1", driver.commanded_on),
        _portable_emt_state_field("$prefix.applied_on", owner, :discrete, "1", driver.applied_on),
        _portable_emt_state_field(
            "$prefix.pending_state",
            owner,
            :scheduler,
            "1",
            driver.pending_state === nothing ? "none" : driver.pending_state ? "on" : "off",
        ),
        _portable_emt_time_field("$prefix.pending_transition_time", owner, :scheduler, driver.pending_transition_time_s),
        _portable_emt_time_field("$prefix.last_command_time", owner, :scheduler, driver.last_command_time_s),
        _portable_emt_time_field("$prefix.last_turn_on_time", owner, :scheduler, driver.last_turn_on_time_s),
        _portable_emt_time_field("$prefix.last_turn_off_time", owner, :scheduler, driver.last_turn_off_time_s),
        _portable_emt_state_field("$prefix.command_count", owner, :discrete, "1", driver.command_count),
        _portable_emt_state_field("$prefix.transition_count", owner, :discrete, "1", driver.transition_count),
        _portable_emt_state_field("$prefix.filtered_pulse_count", owner, :discrete, "1", driver.filtered_pulse_count),
    ]
end

function _restore_portable_semiconductor_gate_driver!(
    driver,
    fields::Dict{String,PortableSnapshotStateField},
    prefix::AbstractString,
)
    driver.commanded_on = _portable_emt_scalar(fields, "$prefix.commanded_on", :discrete, "1", Bool)
    driver.applied_on = _portable_emt_scalar(fields, "$prefix.applied_on", :discrete, "1", Bool)
    pending = _portable_emt_scalar(fields, "$prefix.pending_state", :scheduler, "1", String)
    driver.pending_state = pending == "none" ? nothing :
        pending == "on" ? true : pending == "off" ? false :
        _portable_emt_fail(:invalid_state_value, "portable gate pending state is invalid")
    driver.pending_transition_time_s = _portable_emt_time_scalar(fields, "$prefix.pending_transition_time", :scheduler)
    driver.last_command_time_s = _portable_emt_time_scalar(fields, "$prefix.last_command_time", :scheduler)
    driver.last_turn_on_time_s = _portable_emt_time_scalar(fields, "$prefix.last_turn_on_time", :scheduler)
    driver.last_turn_off_time_s = _portable_emt_time_scalar(fields, "$prefix.last_turn_off_time", :scheduler)
    driver.command_count = _portable_emt_nonnegative_integer(fields, "$prefix.command_count", :discrete)
    driver.transition_count = _portable_emt_nonnegative_integer(fields, "$prefix.transition_count", :discrete)
    driver.filtered_pulse_count = _portable_emt_nonnegative_integer(fields, "$prefix.filtered_pulse_count", :discrete)
    return driver
end

function _portable_semiconductor_snubber_fields(snubber, prefix::AbstractString)
    owner = "emt.power_semiconductor_snubber"
    return PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.resistance", owner, :checkpoint, "ohm", snubber.resistance_ohm),
        _portable_emt_state_field("$prefix.capacitance", owner, :checkpoint, "F", snubber.capacitance_f),
        _portable_emt_state_field("$prefix.previous_current", owner, :history, "A", snubber.previous_current_a),
        _portable_emt_state_field("$prefix.capacitor_voltage", owner, :continuous, "V", snubber.capacitor_voltage_v),
        _portable_emt_state_field("$prefix.last_branch_voltage", owner, :algebraic, "V", snubber.last_branch_voltage_v),
        _portable_emt_state_field("$prefix.last_current", owner, :algebraic, "A", snubber.last_current_a),
        _portable_emt_state_field("$prefix.last_resistor_loss", owner, :output, "W", snubber.last_resistor_loss_w),
        _portable_emt_state_field("$prefix.dissipated_energy", owner, :history, "J", snubber.dissipated_energy_j),
    ]
end

function _restore_portable_semiconductor_snubber!(snubber, fields, prefix)
    for (field_name, identity, family, unit) in (
        (:previous_current_a, :previous_current, :history, "A"),
        (:capacitor_voltage_v, :capacitor_voltage, :continuous, "V"),
        (:last_branch_voltage_v, :last_branch_voltage, :algebraic, "V"),
        (:last_current_a, :last_current, :algebraic, "A"),
        (:last_resistor_loss_w, :last_resistor_loss, :output, "W"),
        (:dissipated_energy_j, :dissipated_energy, :history, "J"),
    )
        setfield!(snubber, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    return snubber
end

function _portable_semiconductor_recovery_fields(recovery, prefix::AbstractString)
    owner = "emt.power_semiconductor_recovery"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.lifetime", owner, :checkpoint, "s", recovery.lifetime_s),
        _portable_emt_state_field("$prefix.stored_charge", owner, :continuous, "C", recovery.stored_charge_c),
        _portable_emt_state_field("$prefix.previous_stored_charge", owner, :history, "C", recovery.previous_stored_charge_c),
        _portable_emt_state_field("$prefix.active", owner, :discrete, "1", recovery.recovery_active),
        _portable_emt_state_field("$prefix.last_current", owner, :algebraic, "A", recovery.last_recovery_current_a),
        _portable_emt_state_field("$prefix.peak_reverse_current", owner, :output, "A", recovery.peak_reverse_current_a),
        _portable_emt_state_field("$prefix.cumulative_recovered_charge", owner, :output, "C", recovery.cumulative_recovered_charge_c),
        _portable_emt_time_field("$prefix.start_time", owner, :scheduler, recovery.recovery_start_time_s),
        _portable_emt_state_field("$prefix.last_duration", owner, :output, "s", recovery.last_recovery_duration_s),
        _portable_emt_state_field("$prefix.zero_event_count", owner, :discrete, "1", recovery.recovery_zero_event_count),
        _portable_emt_time_field("$prefix.last_zero_time", owner, :scheduler, recovery.last_recovery_zero_time_s),
    ]
    append!(fields, _portable_semiconductor_provenance_fields(recovery.provenance, "$prefix.provenance", owner))
    return fields
end

function _restore_portable_semiconductor_recovery!(recovery, fields, prefix)
    for (field_name, identity, family, unit) in (
        (:stored_charge_c, :stored_charge, :continuous, "C"),
        (:previous_stored_charge_c, :previous_stored_charge, :history, "C"),
        (:last_recovery_current_a, :last_current, :algebraic, "A"),
        (:peak_reverse_current_a, :peak_reverse_current, :output, "A"),
        (:cumulative_recovered_charge_c, :cumulative_recovered_charge, :output, "C"),
        (:last_recovery_duration_s, :last_duration, :output, "s"),
    )
        setfield!(recovery, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    recovery.recovery_active = _portable_emt_scalar(fields, "$prefix.active", :discrete, "1", Bool)
    recovery.recovery_start_time_s = _portable_emt_time_scalar(fields, "$prefix.start_time", :scheduler)
    recovery.recovery_zero_event_count = _portable_emt_nonnegative_integer(fields, "$prefix.zero_event_count", :discrete)
    recovery.last_recovery_zero_time_s = _portable_emt_time_scalar(fields, "$prefix.last_zero_time", :scheduler)
    return recovery
end

function _portable_semiconductor_junction_charge_fields(charge, prefix::AbstractString)
    owner = "emt.power_semiconductor_junction_charge"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.zero_bias_capacitance", owner, :checkpoint, "F", charge.zero_bias_capacitance_f),
        _portable_emt_state_field("$prefix.junction_voltage", owner, :checkpoint, "V", charge.junction_voltage_v),
        _portable_emt_state_field("$prefix.grading_exponent", owner, :checkpoint, "1", charge.grading_exponent),
        _portable_emt_state_field("$prefix.minimum_voltage", owner, :checkpoint, "V", charge.minimum_voltage_v),
        _portable_emt_state_field("$prefix.maximum_voltage", owner, :checkpoint, "V", charge.maximum_voltage_v),
        _portable_emt_state_field("$prefix.previous_voltage", owner, :history, "V", charge.previous_voltage_v),
        _portable_emt_state_field("$prefix.previous_charge", owner, :history, "C", charge.previous_charge_c),
        _portable_emt_state_field("$prefix.last_capacitance", owner, :algebraic, "F", charge.last_capacitance_f),
        _portable_emt_state_field("$prefix.last_charge", owner, :algebraic, "C", charge.last_charge_c),
        _portable_emt_state_field("$prefix.last_displacement_current", owner, :algebraic, "A", charge.last_displacement_current_a),
    ]
    append!(fields, _portable_semiconductor_provenance_fields(charge.provenance, "$prefix.provenance", owner))
    return fields
end

function _restore_portable_semiconductor_junction_charge!(charge, fields, prefix)
    for (field_name, identity, family, unit) in (
        (:previous_voltage_v, :previous_voltage, :history, "V"),
        (:previous_charge_c, :previous_charge, :history, "C"),
        (:last_capacitance_f, :last_capacitance, :algebraic, "F"),
        (:last_charge_c, :last_charge, :algebraic, "C"),
        (:last_displacement_current_a, :last_displacement_current, :algebraic, "A"),
    )
        setfield!(charge, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    return charge
end

function _portable_semiconductor_tail_fields(tail, prefix::AbstractString)
    owner = "emt.power_semiconductor_tail"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.decay_time", owner, :checkpoint, "s", tail.decay_time_s),
        _portable_emt_state_field("$prefix.cutoff_current", owner, :checkpoint, "A", tail.cutoff_current_a),
        _portable_emt_state_field("$prefix.active", owner, :discrete, "1", tail.active),
        _portable_emt_state_field("$prefix.current", owner, :continuous, "A", tail.current_a),
        _portable_emt_state_field("$prefix.initial_current", owner, :history, "A", tail.initial_current_a),
        _portable_emt_time_field("$prefix.turn_off_time", owner, :scheduler, tail.turn_off_time_s),
        _portable_emt_state_field("$prefix.last_duration", owner, :output, "s", tail.last_duration_s),
        _portable_emt_state_field("$prefix.cutoff_event_count", owner, :discrete, "1", tail.cutoff_event_count),
        _portable_emt_time_field("$prefix.last_cutoff_time", owner, :scheduler, tail.last_cutoff_time_s),
    ]
    append!(fields, _portable_semiconductor_provenance_fields(tail.provenance, "$prefix.provenance", owner))
    return fields
end

function _restore_portable_semiconductor_tail!(tail, fields, prefix)
    tail.active = _portable_emt_scalar(fields, "$prefix.active", :discrete, "1", Bool)
    tail.current_a = _portable_emt_scalar(fields, "$prefix.current", :continuous, "A", Float64)
    tail.initial_current_a = _portable_emt_scalar(fields, "$prefix.initial_current", :history, "A", Float64)
    tail.turn_off_time_s = _portable_emt_time_scalar(fields, "$prefix.turn_off_time", :scheduler)
    tail.last_duration_s = _portable_emt_scalar(fields, "$prefix.last_duration", :output, "s", Float64)
    tail.cutoff_event_count = _portable_emt_nonnegative_integer(fields, "$prefix.cutoff_event_count", :discrete)
    tail.last_cutoff_time_s = _portable_emt_time_scalar(fields, "$prefix.last_cutoff_time", :scheduler)
    return tail
end

function _portable_semiconductor_energy_fields(energy, prefix::AbstractString)
    owner = "emt.power_semiconductor_switching_energy"
    fields = PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.current_axis", owner, :checkpoint, "A", ["current"], energy.current_axis_a),
        _portable_emt_array_field("$prefix.blocking_voltage_axis", owner, :checkpoint, "V", ["voltage"], energy.blocking_voltage_axis_v),
        _portable_emt_array_field("$prefix.junction_temperature_axis", owner, :checkpoint, "K", ["temperature"], energy.junction_temperature_axis_k),
        _portable_emt_array_field("$prefix.turn_on_table", owner, :checkpoint, "J", ["current", "voltage", "temperature"], energy.turn_on_energy_j),
        _portable_emt_array_field("$prefix.turn_off_table", owner, :checkpoint, "J", ["current", "voltage", "temperature"], energy.turn_off_energy_j),
        _portable_emt_array_field("$prefix.reverse_recovery_table", owner, :checkpoint, "J", ["current", "voltage", "temperature"], energy.reverse_recovery_energy_j),
        _portable_emt_state_field("$prefix.cumulative_turn_on_energy", owner, :output, "J", energy.cumulative_turn_on_energy_j),
        _portable_emt_state_field("$prefix.cumulative_turn_off_energy", owner, :output, "J", energy.cumulative_turn_off_energy_j),
        _portable_emt_state_field("$prefix.cumulative_reverse_recovery_energy", owner, :output, "J", energy.cumulative_reverse_recovery_energy_j),
        _portable_emt_state_field("$prefix.last_event_kind", owner, :discrete, "1", String(energy.last_event_kind)),
        _portable_emt_state_field("$prefix.last_event_energy", owner, :output, "J", energy.last_event_energy_j),
        _portable_emt_state_field("$prefix.last_event_transition_count", owner, :discrete, "1", energy.last_event_transition_count),
        _portable_emt_time_field("$prefix.last_reverse_recovery_start_time", owner, :scheduler, energy.last_reverse_recovery_start_time_s),
    ]
    append!(fields, _portable_semiconductor_provenance_fields(energy.provenance, "$prefix.provenance", owner))
    return fields
end

function _restore_portable_semiconductor_energy!(energy, fields, prefix)
    energy.cumulative_turn_on_energy_j = _portable_emt_scalar(fields, "$prefix.cumulative_turn_on_energy", :output, "J", Float64)
    energy.cumulative_turn_off_energy_j = _portable_emt_scalar(fields, "$prefix.cumulative_turn_off_energy", :output, "J", Float64)
    energy.cumulative_reverse_recovery_energy_j = _portable_emt_scalar(fields, "$prefix.cumulative_reverse_recovery_energy", :output, "J", Float64)
    energy.last_event_kind = Symbol(_portable_emt_scalar(fields, "$prefix.last_event_kind", :discrete, "1", String))
    energy.last_event_energy_j = _portable_emt_scalar(fields, "$prefix.last_event_energy", :output, "J", Float64)
    energy.last_event_transition_count = _portable_emt_nonnegative_integer(fields, "$prefix.last_event_transition_count", :discrete)
    energy.last_reverse_recovery_start_time_s = _portable_emt_time_scalar(fields, "$prefix.last_reverse_recovery_start_time", :scheduler)
    return energy
end

function _portable_semiconductor_thermal_fields(thermal, prefix::AbstractString)
    owner = "emt.power_semiconductor_thermal"
    fields = PortableSnapshotStateField[
        _portable_emt_array_field("$prefix.capacitance", owner, :checkpoint, "J_K", ["stage"], thermal.capacitance_j_per_k),
        _portable_emt_array_field("$prefix.resistance", owner, :checkpoint, "K_W", ["stage"], thermal.resistance_k_per_w),
        _portable_emt_state_field("$prefix.ambient_temperature", owner, :checkpoint, "K", thermal.ambient_temperature_k),
        _portable_emt_array_field("$prefix.node_temperature", owner, :continuous, "K", ["stage"], thermal.node_temperature_k),
        _portable_emt_state_field("$prefix.minimum_temperature", owner, :checkpoint, "K", thermal.minimum_temperature_k),
        _portable_emt_state_field("$prefix.maximum_temperature", owner, :checkpoint, "K", thermal.maximum_temperature_k),
        _portable_emt_state_field("$prefix.last_loss_power", owner, :algebraic, "W", thermal.last_loss_power_w),
        _portable_emt_state_field("$prefix.last_ambient_heat_flow", owner, :algebraic, "W", thermal.last_ambient_heat_flow_w),
        _portable_emt_state_field("$prefix.last_stored_energy", owner, :history, "J", thermal.last_stored_energy_j),
        _portable_emt_state_field("$prefix.cumulative_input_energy", owner, :history, "J", thermal.cumulative_input_energy_j),
        _portable_emt_state_field("$prefix.cumulative_ambient_energy", owner, :history, "J", thermal.cumulative_ambient_energy_j),
        _portable_emt_array_field("$prefix.trial_lower_conductance", owner, :reconstruction, "W_K", ["coupling"], thermal.trial_lower_conductance_w_per_k),
        _portable_emt_array_field("$prefix.trial_diagonal_conductance", owner, :reconstruction, "W_K", ["stage"], thermal.trial_diagonal_conductance_w_per_k),
        _portable_emt_array_field("$prefix.trial_upper_conductance", owner, :reconstruction, "W_K", ["coupling"], thermal.trial_upper_conductance_w_per_k),
        _portable_emt_array_field("$prefix.trial_right_hand_side", owner, :reconstruction, "W", ["stage"], thermal.trial_right_hand_side_w),
        _portable_emt_array_field("$prefix.trial_temperature_rise", owner, :reconstruction, "K", ["stage"], thermal.trial_temperature_rise_k),
        _portable_emt_array_field("$prefix.trial_temperature", owner, :reconstruction, "K", ["stage"], thermal.trial_temperature_k),
        _portable_emt_array_field("$prefix.trial_storage_conductance", owner, :reconstruction, "W_K", ["stage"], thermal.trial_storage_conductance_w_per_k),
        _portable_semiconductor_quantity_field("$prefix.factorized_step", owner, :reconstruction, "s", thermal.factorized_step_s; allow_nan=true),
    ]
    append!(fields, _portable_semiconductor_provenance_fields(thermal.provenance, "$prefix.provenance", owner))
    return fields
end

function _restore_portable_semiconductor_thermal!(thermal, fields, prefix)
    for (field_name, identity, unit, axis) in (
        (:node_temperature_k, :node_temperature, "K", "stage"),
        (:trial_lower_conductance_w_per_k, :trial_lower_conductance, "W_K", "coupling"),
        (:trial_diagonal_conductance_w_per_k, :trial_diagonal_conductance, "W_K", "stage"),
        (:trial_upper_conductance_w_per_k, :trial_upper_conductance, "W_K", "coupling"),
        (:trial_right_hand_side_w, :trial_right_hand_side, "W", "stage"),
        (:trial_temperature_rise_k, :trial_temperature_rise, "K", "stage"),
        (:trial_temperature_k, :trial_temperature, "K", "stage"),
        (:trial_storage_conductance_w_per_k, :trial_storage_conductance, "W_K", "stage"),
    )
        _restore_portable_emt_array!(
            getfield(thermal, field_name),
            fields,
            "$prefix.$(String(identity))",
            field_name === :node_temperature_k ? :continuous : :reconstruction,
            unit,
            [axis],
        )
    end
    for (field_name, identity, family, unit) in (
        (:last_loss_power_w, :last_loss_power, :algebraic, "W"),
        (:last_ambient_heat_flow_w, :last_ambient_heat_flow, :algebraic, "W"),
        (:last_stored_energy_j, :last_stored_energy, :history, "J"),
        (:cumulative_input_energy_j, :cumulative_input_energy, :history, "J"),
        (:cumulative_ambient_energy_j, :cumulative_ambient_energy, :history, "J"),
    )
        setfield!(thermal, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    thermal.factorized_step_s = _portable_semiconductor_quantity_scalar(
        fields, "$prefix.factorized_step", :reconstruction, "s"; allow_nan=true,
    )
    return thermal
end

function _portable_semiconductor_extended_fields(fidelity, prefix::AbstractString)
    owner = "emt.power_semiconductor_extended"
    fields = PortableSnapshotStateField[]
    for (component, value) in (
        ("recovered_charge", fidelity.recovered_charge),
        ("junction_charge", fidelity.junction_charge),
        ("turn_off_tail", fidelity.turn_off_tail),
        ("switching_energy", fidelity.switching_energy),
        ("thermal", fidelity.thermal),
        ("declared_model", fidelity.declared_model),
    )
        push!(fields, _portable_semiconductor_optional_presence_field(prefix, component, value))
    end
    append!(fields, _portable_semiconductor_provenance_fields(fidelity.provenance, "$prefix.provenance", owner))
    fidelity.recovered_charge === nothing || append!(fields,
        _portable_semiconductor_recovery_fields(fidelity.recovered_charge, "$prefix.recovered_charge"))
    fidelity.junction_charge === nothing || append!(fields,
        _portable_semiconductor_junction_charge_fields(fidelity.junction_charge, "$prefix.junction_charge"))
    fidelity.turn_off_tail === nothing || append!(fields,
        _portable_semiconductor_tail_fields(fidelity.turn_off_tail, "$prefix.turn_off_tail"))
    fidelity.switching_energy === nothing || append!(fields,
        _portable_semiconductor_energy_fields(fidelity.switching_energy, "$prefix.switching_energy"))
    fidelity.thermal === nothing || append!(fields,
        _portable_semiconductor_thermal_fields(fidelity.thermal, "$prefix.thermal"))
    if fidelity.declared_model !== nothing
        model = fidelity.declared_model
        append!(fields, PortableSnapshotStateField[
            _portable_emt_state_field("$prefix.declared_model.schema", owner, :checkpoint, "1", String(model.schema)),
            _portable_emt_state_field("$prefix.declared_model.version", owner, :checkpoint, "1", string(model.version)),
            _portable_emt_state_field("$prefix.declared_model.source_identity", owner, :checkpoint, "1", model.source_identity),
            _portable_emt_state_field("$prefix.declared_model.content_sha256", owner, :checkpoint, "1", model.content_sha256),
            _portable_emt_state_field("$prefix.declared_model.licence", owner, :checkpoint, "1", model.licence),
            _portable_emt_state_field("$prefix.declared_model.redistribution", owner, :checkpoint, "1", String(model.redistribution)),
        ])
        append!(fields, _portable_semiconductor_provenance_fields(
            model.provenance, "$prefix.declared_model.provenance", owner,
        ))
    end
    append!(fields, PortableSnapshotStateField[
        _portable_emt_time_field("$prefix.candidate_time", owner, :reconstruction, fidelity.candidate_time_s),
        _portable_emt_state_field("$prefix.candidate_step", owner, :reconstruction, "s", fidelity.candidate_step_s),
        _portable_emt_state_field("$prefix.candidate_method", owner, :reconstruction, "1", String(fidelity.candidate_method)),
        _portable_emt_state_field("$prefix.candidate_prepared", owner, :reconstruction, "1", fidelity.candidate_prepared),
        _portable_emt_state_field("$prefix.previous_terminal_voltage", owner, :history, "V", fidelity.previous_terminal_voltage_v),
        _portable_emt_state_field("$prefix.previous_terminal_current", owner, :history, "A", fidelity.previous_terminal_current_a),
        _portable_emt_state_field("$prefix.companion_energy_residual", owner, :output, "J", fidelity.companion_energy_residual_j),
        _portable_emt_state_field("$prefix.accepted_topology_transition_count", owner, :discrete, "1", fidelity.accepted_topology_transition_count),
        _portable_emt_state_field("$prefix.pending_event_current", owner, :scheduler, "A", fidelity.pending_event_current_a),
        _portable_emt_state_field("$prefix.pending_event_blocking_voltage", owner, :scheduler, "V", fidelity.pending_event_blocking_voltage_v),
        _portable_emt_state_field("$prefix.candidate_recovery_charge", owner, :reconstruction, "C", fidelity.candidate_recovery_charge_c),
    ])
    return fields
end

function _restore_portable_semiconductor_extended!(fidelity, fields, prefix)
    for (component, value) in (
        ("recovered_charge", fidelity.recovered_charge),
        ("junction_charge", fidelity.junction_charge),
        ("turn_off_tail", fidelity.turn_off_tail),
        ("switching_energy", fidelity.switching_energy),
        ("thermal", fidelity.thermal),
        ("declared_model", fidelity.declared_model),
    )
        _portable_semiconductor_require_presence(fields, prefix, component, value)
    end
    fidelity.recovered_charge === nothing || _restore_portable_semiconductor_recovery!(
        fidelity.recovered_charge, fields, "$prefix.recovered_charge",
    )
    fidelity.junction_charge === nothing || _restore_portable_semiconductor_junction_charge!(
        fidelity.junction_charge, fields, "$prefix.junction_charge",
    )
    fidelity.turn_off_tail === nothing || _restore_portable_semiconductor_tail!(
        fidelity.turn_off_tail, fields, "$prefix.turn_off_tail",
    )
    fidelity.switching_energy === nothing || _restore_portable_semiconductor_energy!(
        fidelity.switching_energy, fields, "$prefix.switching_energy",
    )
    fidelity.thermal === nothing || _restore_portable_semiconductor_thermal!(
        fidelity.thermal, fields, "$prefix.thermal",
    )
    fidelity.candidate_time_s = _portable_emt_time_scalar(fields, "$prefix.candidate_time", :reconstruction)
    fidelity.candidate_step_s = _portable_emt_scalar(fields, "$prefix.candidate_step", :reconstruction, "s", Float64)
    fidelity.candidate_method = Symbol(_portable_emt_scalar(fields, "$prefix.candidate_method", :reconstruction, "1", String))
    fidelity.candidate_prepared = _portable_emt_scalar(fields, "$prefix.candidate_prepared", :reconstruction, "1", Bool)
    for (field_name, identity, family, unit) in (
        (:previous_terminal_voltage_v, :previous_terminal_voltage, :history, "V"),
        (:previous_terminal_current_a, :previous_terminal_current, :history, "A"),
        (:companion_energy_residual_j, :companion_energy_residual, :output, "J"),
        (:pending_event_current_a, :pending_event_current, :scheduler, "A"),
        (:pending_event_blocking_voltage_v, :pending_event_blocking_voltage, :scheduler, "V"),
        (:candidate_recovery_charge_c, :candidate_recovery_charge, :reconstruction, "C"),
    )
        setfield!(fidelity, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    fidelity.accepted_topology_transition_count = _portable_emt_nonnegative_integer(
        fields, "$prefix.accepted_topology_transition_count", :discrete,
    )
    return fidelity
end

function _portable_power_semiconductor_switch_fields(switch, prefix::AbstractString)
    owner = "emt.power_semiconductor_switch"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.device_kind", owner, :checkpoint, "1", string(typeof(switch).parameters[1])),
        _portable_emt_state_field("$prefix.from_node", owner, :checkpoint, "1", switch.a),
        _portable_emt_state_field("$prefix.to_node", owner, :checkpoint, "1", switch.b),
        _portable_emt_state_field("$prefix.threshold_voltage", owner, :checkpoint, "V", switch.threshold_v),
        _portable_emt_state_field("$prefix.holding_current", owner, :checkpoint, "A", switch.holding_current),
        _portable_emt_state_field("$prefix.on_conductance", owner, :checkpoint, "S", switch.on_conductance),
        _portable_emt_state_field("$prefix.off_conductance", owner, :checkpoint, "S", switch.off_conductance),
        _portable_emt_state_field("$prefix.forward_voltage_drop", owner, :checkpoint, "V", switch.forward_voltage_drop_v),
        _portable_semiconductor_optional_presence_field(prefix, "gate_driver", switch.gate_driver),
        _portable_semiconductor_optional_presence_field(prefix, "antiparallel_diode", switch.antiparallel_diode),
        _portable_semiconductor_optional_presence_field(prefix, "snubber", switch.snubber),
        _portable_semiconductor_optional_presence_field(prefix, "extended_fidelity", switch.extended_fidelity),
        _portable_emt_state_field("$prefix.closed", owner, :discrete, "1", switch.closed),
        _portable_emt_state_field("$prefix.last_voltage", owner, :algebraic, "V", switch.last_voltage),
        _portable_emt_state_field("$prefix.last_current", owner, :algebraic, "A", switch.last_current),
        _portable_emt_state_field("$prefix.last_conductance", owner, :algebraic, "S", switch.last_conductance),
        _portable_emt_state_field("$prefix.reverse_diode_conducting", owner, :discrete, "1", switch.reverse_diode_conducting),
        _portable_emt_state_field("$prefix.event_localization_enabled", owner, :discrete, "1", switch.event_localization_enabled),
        _portable_emt_time_field("$prefix.last_evaluation_time", owner, :scheduler, switch.last_evaluation_time_s),
        _portable_emt_state_field("$prefix.last_history_current", owner, :history, "A", switch.last_history_current_a),
        _portable_emt_state_field("$prefix.last_forward_current", owner, :algebraic, "A", switch.last_forward_current_a),
        _portable_emt_state_field("$prefix.last_reverse_diode_current", owner, :algebraic, "A", switch.last_reverse_diode_current_a),
        _portable_emt_state_field("$prefix.last_snubber_current", owner, :algebraic, "A", switch.last_snubber_current_a),
        _portable_emt_state_field("$prefix.last_semiconductor_loss", owner, :output, "W", switch.last_semiconductor_loss_w),
        _portable_emt_state_field("$prefix.previous_semiconductor_loss", owner, :history, "W", switch.previous_semiconductor_loss_w),
        _portable_emt_state_field("$prefix.semiconductor_dissipated_energy", owner, :history, "J", switch.semiconductor_dissipated_energy_j),
        _portable_emt_state_field("$prefix.topology_transition_count", owner, :discrete, "1", switch.topology_transition_count),
        _portable_emt_time_field("$prefix.last_transition_time", owner, :scheduler, switch.last_transition_time_s),
        _portable_emt_state_field("$prefix.conduction_direction", owner, :discrete, "1", Int(switch.conduction_direction)),
        _portable_semiconductor_quantity_field("$prefix.gate_turn_off_current_limit", owner, :checkpoint, "A", switch.gate_turn_off_current_limit_a),
        _portable_emt_state_field("$prefix.gate_turn_off_policy", owner, :checkpoint, "1", String(switch.gate_turn_off_policy)),
        _portable_emt_state_field("$prefix.gate_turn_off_disposition", owner, :discrete, "1", String(switch.gate_turn_off_disposition)),
    ]
    if switch.antiparallel_diode !== nothing
        diode = switch.antiparallel_diode
        append!(fields, PortableSnapshotStateField[
            _portable_emt_state_field("$prefix.antiparallel_diode.forward_voltage", owner, :checkpoint, "V", diode.forward_voltage_v),
            _portable_emt_state_field("$prefix.antiparallel_diode.holding_current", owner, :checkpoint, "A", diode.holding_current_a),
            _portable_emt_state_field("$prefix.antiparallel_diode.on_conductance", owner, :checkpoint, "S", diode.on_conductance_s),
        ])
    end
    switch.gate_driver === nothing || append!(fields,
        _portable_semiconductor_gate_driver_fields(switch.gate_driver, "$prefix.gate_driver"))
    switch.snubber === nothing || append!(fields,
        _portable_semiconductor_snubber_fields(switch.snubber, "$prefix.snubber"))
    switch.extended_fidelity === nothing || append!(fields,
        _portable_semiconductor_extended_fields(switch.extended_fidelity, "$prefix.extended_fidelity"))
    return fields
end

function _restore_portable_power_semiconductor_switch!(switch, fields, prefix)
    _portable_semiconductor_require_presence(fields, prefix, "gate_driver", switch.gate_driver)
    _portable_semiconductor_require_presence(fields, prefix, "antiparallel_diode", switch.antiparallel_diode)
    _portable_semiconductor_require_presence(fields, prefix, "snubber", switch.snubber)
    _portable_semiconductor_require_presence(fields, prefix, "extended_fidelity", switch.extended_fidelity)
    switch.gate_driver === nothing || _restore_portable_semiconductor_gate_driver!(
        switch.gate_driver, fields, "$prefix.gate_driver",
    )
    switch.snubber === nothing || _restore_portable_semiconductor_snubber!(
        switch.snubber, fields, "$prefix.snubber",
    )
    switch.extended_fidelity === nothing || _restore_portable_semiconductor_extended!(
        switch.extended_fidelity, fields, "$prefix.extended_fidelity",
    )
    for (field_name, identity, family, unit) in (
        (:last_voltage, :last_voltage, :algebraic, "V"),
        (:last_current, :last_current, :algebraic, "A"),
        (:last_conductance, :last_conductance, :algebraic, "S"),
        (:last_history_current_a, :last_history_current, :history, "A"),
        (:last_forward_current_a, :last_forward_current, :algebraic, "A"),
        (:last_reverse_diode_current_a, :last_reverse_diode_current, :algebraic, "A"),
        (:last_snubber_current_a, :last_snubber_current, :algebraic, "A"),
        (:last_semiconductor_loss_w, :last_semiconductor_loss, :output, "W"),
        (:previous_semiconductor_loss_w, :previous_semiconductor_loss, :history, "W"),
        (:semiconductor_dissipated_energy_j, :semiconductor_dissipated_energy, :history, "J"),
    )
        setfield!(switch, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", family, unit, Float64,
        ))
    end
    switch.closed = _portable_emt_scalar(fields, "$prefix.closed", :discrete, "1", Bool)
    switch.reverse_diode_conducting = _portable_emt_scalar(fields, "$prefix.reverse_diode_conducting", :discrete, "1", Bool)
    switch.event_localization_enabled = _portable_emt_scalar(fields, "$prefix.event_localization_enabled", :discrete, "1", Bool)
    switch.last_evaluation_time_s = _portable_emt_time_scalar(fields, "$prefix.last_evaluation_time", :scheduler)
    switch.topology_transition_count = _portable_emt_nonnegative_integer(fields, "$prefix.topology_transition_count", :discrete)
    switch.last_transition_time_s = _portable_emt_time_scalar(fields, "$prefix.last_transition_time", :scheduler)
    direction = _portable_emt_integer(fields, "$prefix.conduction_direction", :discrete)
    typemin(Int8) <= direction <= typemax(Int8) || _portable_emt_fail(
        :invalid_state_value,
        "portable semiconductor conduction direction is outside Int8",
    )
    switch.conduction_direction = Int8(direction)
    switch.gate_turn_off_disposition = Symbol(_portable_emt_scalar(
        fields, "$prefix.gate_turn_off_disposition", :discrete, "1", String,
    ))
    return switch
end

function _portable_emt_element_state_fields(
    switch::PowerSemiconductorSwitch,
    index::Integer,
    name::Symbol,
)
    fields = _portable_emt_element_identity_field(index, name, "power_semiconductor_switch")
    append!(fields, _portable_power_semiconductor_switch_fields(
        switch, _portable_emt_element_prefix(index),
    ))
    return fields
end

function _restore_portable_emt_element_state!(
    switch::PowerSemiconductorSwitch,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    return _restore_portable_power_semiconductor_switch!(
        switch, fields, _portable_emt_element_prefix(index),
    )
end

function _portable_power_semiconductor_bridge_leg_fields(
    bridge,
    prefix::AbstractString;
    include_switches::Bool,
)
    owner = "emt.power_semiconductor_bridge_leg"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.commutation_dead_time", owner, :checkpoint, "s", bridge.commutation_dead_time_s),
        _portable_emt_state_field("$prefix.blocked", owner, :discrete, "1", bridge.blocked),
        _portable_emt_state_field("$prefix.requested_upper_on", owner, :discrete, "1", bridge.requested_upper_on),
        _portable_emt_state_field("$prefix.requested_lower_on", owner, :discrete, "1", bridge.requested_lower_on),
        _portable_emt_time_field("$prefix.last_command_time", owner, :scheduler, bridge.last_command_time_s),
        _portable_emt_state_field("$prefix.command_count", owner, :discrete, "1", bridge.command_count),
        _portable_emt_state_field("$prefix.shoot_through_rejection_count", owner, :discrete, "1", bridge.shoot_through_rejection_count),
        _portable_emt_state_field("$prefix.block_count", owner, :discrete, "1", bridge.block_count),
        _portable_emt_state_field("$prefix.restart_count", owner, :discrete, "1", bridge.restart_count),
        _portable_emt_state_field("$prefix.last_dc_positive_voltage", owner, :algebraic, "V", bridge.last_dc_positive_voltage_v),
        _portable_emt_state_field("$prefix.last_ac_terminal_voltage", owner, :algebraic, "V", bridge.last_ac_terminal_voltage_v),
        _portable_emt_state_field("$prefix.last_dc_negative_voltage", owner, :algebraic, "V", bridge.last_dc_negative_voltage_v),
    ]
    if include_switches
        append!(fields, _portable_power_semiconductor_switch_fields(
            bridge.upper_switch, "$prefix.upper_switch",
        ))
        append!(fields, _portable_power_semiconductor_switch_fields(
            bridge.lower_switch, "$prefix.lower_switch",
        ))
    end
    return fields
end

function _restore_portable_power_semiconductor_bridge_leg!(
    bridge,
    fields,
    prefix;
    include_switches::Bool,
)
    include_switches && begin
        _restore_portable_power_semiconductor_switch!(
            bridge.upper_switch, fields, "$prefix.upper_switch",
        )
        _restore_portable_power_semiconductor_switch!(
            bridge.lower_switch, fields, "$prefix.lower_switch",
        )
    end
    bridge.blocked = _portable_emt_scalar(fields, "$prefix.blocked", :discrete, "1", Bool)
    bridge.requested_upper_on = _portable_emt_scalar(fields, "$prefix.requested_upper_on", :discrete, "1", Bool)
    bridge.requested_lower_on = _portable_emt_scalar(fields, "$prefix.requested_lower_on", :discrete, "1", Bool)
    bridge.last_command_time_s = _portable_emt_time_scalar(fields, "$prefix.last_command_time", :scheduler)
    for (field_name, identity) in (
        (:command_count, :command_count),
        (:shoot_through_rejection_count, :shoot_through_rejection_count),
        (:block_count, :block_count),
        (:restart_count, :restart_count),
    )
        setfield!(bridge, field_name, _portable_emt_nonnegative_integer(
            fields, "$prefix.$(String(identity))", :discrete,
        ))
    end
    for (field_name, identity) in (
        (:last_dc_positive_voltage_v, :last_dc_positive_voltage),
        (:last_ac_terminal_voltage_v, :last_ac_terminal_voltage),
        (:last_dc_negative_voltage_v, :last_dc_negative_voltage),
    )
        setfield!(bridge, field_name, _portable_emt_scalar(
            fields, "$prefix.$(String(identity))", :algebraic, "V", Float64,
        ))
    end
    return bridge
end

function _portable_emt_element_state_fields(
    bridge::PowerSemiconductorBridgeLeg,
    index::Integer,
    name::Symbol,
)
    fields = _portable_emt_element_identity_field(index, name, "power_semiconductor_bridge_leg")
    append!(fields, _portable_power_semiconductor_bridge_leg_fields(
        bridge, _portable_emt_element_prefix(index); include_switches=true,
    ))
    return fields
end

function _restore_portable_emt_element_state!(
    bridge::PowerSemiconductorBridgeLeg,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    return _restore_portable_power_semiconductor_bridge_leg!(
        bridge, fields, _portable_emt_element_prefix(index); include_switches=true,
    )
end

function _portable_semiconductor_rekey_field(
    field::PortableSnapshotStateField,
    source_prefix::AbstractString,
    destination_prefix::AbstractString,
)
    startswith(field.identity, source_prefix) || _portable_emt_fail(
        :invalid_state_owner,
        "portable nested element field does not use its expected source prefix",
    )
    identity = destination_prefix * field.identity[(lastindex(source_prefix) + 1):end]
    return PortableSnapshotStateField(
        identity,
        field.owner,
        field.family,
        field.unit,
        field.axes,
        field.value,
    )
end

function _portable_nested_element_fields(element, name::Symbol, prefix::AbstractString)
    source_prefix = _portable_emt_element_prefix(1)
    return PortableSnapshotStateField[
        _portable_semiconductor_rekey_field(field, source_prefix, prefix)
        for field in _portable_emt_element_state_fields(element, 1, name)
    ]
end

function _restore_portable_nested_element!(element, fields, prefix::AbstractString)
    source_prefix = _portable_emt_element_prefix(1)
    local_fields = Dict{String,PortableSnapshotStateField}()
    for field in values(fields)
        startswith(field.identity, prefix) || continue
        rekeyed = _portable_semiconductor_rekey_field(field, prefix, source_prefix)
        local_fields[rekeyed.identity] = rekeyed
    end
    _restore_portable_emt_element_state!(element, local_fields, 1)
    return element
end

function _portable_power_semiconductor_topology_fields(bridge, prefix::AbstractString)
    owner = "emt.power_semiconductor_bridge_topology"
    fields = PortableSnapshotStateField[
        _portable_emt_state_field("$prefix.topology_signature", owner, :checkpoint, "1", bridge.topology.topology_signature),
        _portable_emt_state_field("$prefix.valve_count", owner, :checkpoint, "1", length(bridge.valves)),
        _portable_emt_state_field("$prefix.passive_count", owner, :checkpoint, "1", length(bridge.passives)),
        _portable_emt_state_field("$prefix.bridge_leg_count", owner, :checkpoint, "1", length(bridge.bridge_legs)),
        _portable_emt_state_field("$prefix.blocked", owner, :discrete, "1", bridge.blocked),
        _portable_emt_array_field("$prefix.position_fault", owner, :discrete, "1", ["valve"], Int.(bridge.position_faults)),
        _portable_emt_state_field("$prefix.transition_count", owner, :discrete, "1", bridge.transition_count),
        _portable_emt_state_field("$prefix.refusal_count", owner, :discrete, "1", bridge.refusal_count),
        _portable_emt_array_field("$prefix.last_terminal_voltage", owner, :algebraic, "V", ["terminal"], bridge.last_terminal_voltage_v),
        _portable_emt_array_field("$prefix.last_terminal_current", owner, :algebraic, "A", ["terminal"], bridge.last_terminal_current_a),
        _portable_emt_state_field("$prefix.last_step", owner, :history, "s", bridge.last_step_s),
        _portable_emt_state_field("$prefix.dissipated_energy", owner, :history, "J", bridge.dissipated_energy_j),
    ]
    for (valve_index, valve) in enumerate(bridge.valves)
        valve_prefix = "$prefix.valve.i" * lpad(string(valve_index), 8, '0')
        append!(fields, _portable_power_semiconductor_switch_fields(valve, valve_prefix))
    end
    for (passive_index, passive) in enumerate(bridge.passives)
        passive_prefix = "$prefix.passive.i" * lpad(string(passive_index), 8, '0')
        passive_name = bridge.topology.passive_positions[passive_index].name
        append!(fields, _portable_nested_element_fields(passive, passive_name, passive_prefix))
    end
    for (leg_index, leg) in enumerate(bridge.bridge_legs)
        leg_prefix = "$prefix.bridge_leg.i" * lpad(string(leg_index), 8, '0')
        upper_index = findfirst(valve -> valve === leg.upper_switch, bridge.valves)
        lower_index = findfirst(valve -> valve === leg.lower_switch, bridge.valves)
        (upper_index === nothing || lower_index === nothing) && _portable_emt_fail(
            :invalid_state_owner,
            "portable bridge leg does not alias topology valves",
        )
        upper_index = something(upper_index)
        lower_index = something(lower_index)
        append!(fields, PortableSnapshotStateField[
            _portable_emt_state_field("$leg_prefix.upper_valve_index", owner, :checkpoint, "1", upper_index),
            _portable_emt_state_field("$leg_prefix.lower_valve_index", owner, :checkpoint, "1", lower_index),
        ])
        append!(fields, _portable_power_semiconductor_bridge_leg_fields(
            leg, leg_prefix; include_switches=false,
        ))
    end
    return fields
end

function _restore_portable_power_semiconductor_topology!(bridge, fields, prefix)
    _portable_emt_scalar(fields, "$prefix.topology_signature", :checkpoint, "1", String) ==
        bridge.topology.topology_signature || _portable_emt_fail(
            :model_mismatch,
            "portable semiconductor bridge topology signature changed",
        )
    for (identity, count) in (
        (:valve_count, length(bridge.valves)),
        (:passive_count, length(bridge.passives)),
        (:bridge_leg_count, length(bridge.bridge_legs)),
    )
        _portable_emt_integer(fields, "$prefix.$(String(identity))", :checkpoint) == count ||
            _portable_emt_fail(
                :model_mismatch,
                "portable semiconductor topology $identity changed",
            )
    end
    for (valve_index, valve) in enumerate(bridge.valves)
        valve_prefix = "$prefix.valve.i" * lpad(string(valve_index), 8, '0')
        _restore_portable_power_semiconductor_switch!(valve, fields, valve_prefix)
    end
    for (passive_index, passive) in enumerate(bridge.passives)
        passive_prefix = "$prefix.passive.i" * lpad(string(passive_index), 8, '0')
        _restore_portable_nested_element!(passive, fields, passive_prefix)
    end
    for (leg_index, leg) in enumerate(bridge.bridge_legs)
        leg_prefix = "$prefix.bridge_leg.i" * lpad(string(leg_index), 8, '0')
        upper_index = _portable_emt_integer(fields, "$leg_prefix.upper_valve_index", :checkpoint)
        lower_index = _portable_emt_integer(fields, "$leg_prefix.lower_valve_index", :checkpoint)
        1 <= upper_index <= length(bridge.valves) && leg.upper_switch === bridge.valves[upper_index] ||
            _portable_emt_fail(:model_mismatch, "portable upper bridge-leg valve alias changed")
        1 <= lower_index <= length(bridge.valves) && leg.lower_switch === bridge.valves[lower_index] ||
            _portable_emt_fail(:model_mismatch, "portable lower bridge-leg valve alias changed")
        _restore_portable_power_semiconductor_bridge_leg!(
            leg, fields, leg_prefix; include_switches=false,
        )
    end
    bridge.blocked = _portable_emt_scalar(fields, "$prefix.blocked", :discrete, "1", Bool)
    faults = _portable_emt_array_values(
        fields, "$prefix.position_fault", :discrete, "1", ["valve"], Int64,
    )
    length(faults) == length(bridge.position_faults) || _portable_emt_fail(
        :state_shape_mismatch,
        "portable semiconductor position-fault count changed",
    )
    admitted_faults = instances(eltype(bridge.position_faults))
    for index in eachindex(faults, bridge.position_faults)
        match = findfirst(fault -> Int(fault) == faults[index], admitted_faults)
        match === nothing && _portable_emt_fail(
            :invalid_state_value,
            "portable semiconductor position fault is invalid",
        )
        bridge.position_faults[index] = admitted_faults[match]
    end
    bridge.transition_count = _portable_emt_nonnegative_integer(fields, "$prefix.transition_count", :discrete)
    bridge.refusal_count = _portable_emt_nonnegative_integer(fields, "$prefix.refusal_count", :discrete)
    _restore_portable_emt_array!(bridge.last_terminal_voltage_v, fields, "$prefix.last_terminal_voltage", :algebraic, "V", ["terminal"])
    _restore_portable_emt_array!(bridge.last_terminal_current_a, fields, "$prefix.last_terminal_current", :algebraic, "A", ["terminal"])
    bridge.last_step_s = _portable_emt_scalar(fields, "$prefix.last_step", :history, "s", Float64)
    bridge.dissipated_energy_j = _portable_emt_scalar(fields, "$prefix.dissipated_energy", :history, "J", Float64)
    return bridge
end

function _portable_emt_element_state_fields(
    bridge::PowerSemiconductorBridgeTopology,
    index::Integer,
    name::Symbol,
)
    fields = _portable_emt_element_identity_field(index, name, "power_semiconductor_bridge_topology")
    append!(fields, _portable_power_semiconductor_topology_fields(
        bridge, _portable_emt_element_prefix(index),
    ))
    return fields
end

function _restore_portable_emt_element_state!(
    bridge::PowerSemiconductorBridgeTopology,
    fields::Dict{String,PortableSnapshotStateField},
    index::Integer,
)
    return _restore_portable_power_semiconductor_topology!(
        bridge, fields, _portable_emt_element_prefix(index),
    )
end

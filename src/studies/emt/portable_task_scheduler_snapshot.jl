const _PORTABLE_EMT_TASK_SCHEDULER_OWNER = "emt.task_scheduler"

portable_emt_task_state(::Nothing) = nothing
portable_emt_task_state(value::Bool) = value
portable_emt_task_state(value::Signed) = value
portable_emt_task_state(value::Unsigned) = value
portable_emt_task_state(value::AbstractFloat) = value
portable_emt_task_state(value::AbstractString) = String(value)
portable_emt_task_state(value::Rational) = value
portable_emt_task_state(value::Symbol) = String(value)
portable_emt_task_state(value::PortableSnapshotArray) = value
portable_emt_task_state(value::PortableSnapshotRecord) = value
portable_emt_task_state(value::Tuple) = tuple((portable_emt_task_state(item) for item in value)...)
portable_emt_task_state(value::AbstractVector) =
    Any[portable_emt_task_state(item) for item in value]

function portable_emt_task_state(
    value::SwitchDetailedVSC.ExtendedVSCControlCommand,
)
    return PortableSnapshotRecord(
        "aimora.emt.extended_vsc_control_command.v1",
        Pair{String,Any}[
            "active_power_w" => value.active_power_w,
            "angle_rad" => value.angle_rad,
            "controller_family" => Int64(Int(value.controller_family)),
            "duties" => portable_emt_task_state(value.duties),
            "frequency_hz" => value.frequency_hz,
            "limited" => value.limited,
            "mode" => Int64(Int(value.mode)),
            "negative_sequence_voltage_v" => value.negative_sequence_voltage_v,
            "phase_current_reference_a" => portable_emt_task_state(
                value.phase_current_reference_a,
            ),
            "pole_voltage_reference_v" => portable_emt_task_state(
                value.pole_voltage_reference_v,
            ),
            "positive_sequence_voltage_v" => value.positive_sequence_voltage_v,
            "reactive_power_var" => value.reactive_power_var,
            "request_disposition" => Int64(Int(value.request_disposition)),
            "sequence_extractor_settled" => value.sequence_extractor_settled,
            "wire_form" => Int64(Int(value.wire_form)),
            "zero_sequence_voltage_v" => value.zero_sequence_voltage_v,
        ],
    )
end

function portable_emt_task_state(value::Base.RefValue)
    return PortableSnapshotRecord(
        "aimora.emt.task_ref.v1",
        Pair{String,Any}["value" => portable_emt_task_state(value[])],
    )
end

function portable_emt_task_state(value)
    _portable_emt_fail(
        :unsupported_task_state,
        "task state $(typeof(value)) requires an explicit portable_emt_task_state method",
    )
end

restore_portable_emt_task_state(::Nothing, ::Nothing) = nothing
restore_portable_emt_task_state(::Bool, value::Bool) = value
restore_portable_emt_task_state(template::T, value::Signed) where {T<:Signed} =
    convert(T, value)
restore_portable_emt_task_state(template::T, value::Unsigned) where {T<:Unsigned} =
    convert(T, value)
restore_portable_emt_task_state(template::T, value::Real) where {T<:AbstractFloat} =
    convert(T, value)
restore_portable_emt_task_state(::T, value::AbstractString) where {T<:AbstractString} =
    convert(T, value)
restore_portable_emt_task_state(::Symbol, value::AbstractString) = Symbol(value)
restore_portable_emt_task_state(::T, value::Rational) where {T<:Rational} =
    convert(T, value)
restore_portable_emt_task_state(::PortableSnapshotArray, value::PortableSnapshotArray) = value
restore_portable_emt_task_state(::PortableSnapshotRecord, value::PortableSnapshotRecord) = value

function restore_portable_emt_task_state(
    template::Tuple,
    value::Union{Tuple,AbstractVector},
)
    length(template) == length(value) || _portable_emt_fail(
        :task_state_shape_mismatch,
        "portable task tuple state changed length",
    )
    return tuple((
        restore_portable_emt_task_state(template[index], value[index])
        for index in eachindex(template)
    )...)
end

function restore_portable_emt_task_state(template::AbstractVector{T}, value::AbstractVector) where {T}
    restored = Vector{T}(undef, length(value))
    if isempty(template) && !isempty(value)
        _portable_emt_fail(
            :task_state_shape_mismatch,
            "portable task vector needs a typed template element",
        )
    end
    for index in eachindex(value)
        exemplar = index <= length(template) ? template[index] : first(template)
        restored[index] = restore_portable_emt_task_state(exemplar, value[index])
    end
    return restored
end

function _restore_portable_emt_task_enum(template::T, value) where {T<:Enum}
    value isa Signed || _portable_emt_fail(
        :task_state_type_mismatch,
        "portable task enumeration has the wrong type",
    )
    return try
        T(value)
    catch error
        _portable_emt_fail(
            :task_state_type_mismatch,
            "portable task enumeration is invalid: $(sprint(showerror, error))",
        )
    end
end

function restore_portable_emt_task_state(
    template::SwitchDetailedVSC.ExtendedVSCControlCommand,
    value::PortableSnapshotRecord,
)
    required = [
        "active_power_w",
        "angle_rad",
        "controller_family",
        "duties",
        "frequency_hz",
        "limited",
        "mode",
        "negative_sequence_voltage_v",
        "phase_current_reference_a",
        "pole_voltage_reference_v",
        "positive_sequence_voltage_v",
        "reactive_power_var",
        "request_disposition",
        "sequence_extractor_settled",
        "wire_form",
        "zero_sequence_voltage_v",
    ]
    fields = _portable_emt_scheduler_record_fields(
        value,
        "aimora.emt.extended_vsc_control_command.v1",
        required,
    )
    return SwitchDetailedVSC.ExtendedVSCControlCommand(
        _restore_portable_emt_task_enum(
            template.controller_family,
            fields["controller_family"],
        ),
        _restore_portable_emt_task_enum(template.wire_form, fields["wire_form"]),
        restore_portable_emt_task_state(template.duties, fields["duties"]),
        restore_portable_emt_task_state(
            template.pole_voltage_reference_v,
            fields["pole_voltage_reference_v"],
        ),
        restore_portable_emt_task_state(
            template.phase_current_reference_a,
            fields["phase_current_reference_a"],
        ),
        restore_portable_emt_task_state(template.angle_rad, fields["angle_rad"]),
        restore_portable_emt_task_state(
            template.frequency_hz,
            fields["frequency_hz"],
        ),
        restore_portable_emt_task_state(
            template.active_power_w,
            fields["active_power_w"],
        ),
        restore_portable_emt_task_state(
            template.reactive_power_var,
            fields["reactive_power_var"],
        ),
        restore_portable_emt_task_state(
            template.positive_sequence_voltage_v,
            fields["positive_sequence_voltage_v"],
        ),
        restore_portable_emt_task_state(
            template.negative_sequence_voltage_v,
            fields["negative_sequence_voltage_v"],
        ),
        restore_portable_emt_task_state(
            template.zero_sequence_voltage_v,
            fields["zero_sequence_voltage_v"],
        ),
        restore_portable_emt_task_state(
            template.sequence_extractor_settled,
            fields["sequence_extractor_settled"],
        ),
        restore_portable_emt_task_state(template.limited, fields["limited"]),
        _restore_portable_emt_task_enum(template.mode, fields["mode"]),
        _restore_portable_emt_task_enum(
            template.request_disposition,
            fields["request_disposition"],
        ),
    )
end

function _portable_emt_scheduler_record_fields(
    record::PortableSnapshotRecord,
    schema::AbstractString,
    required::AbstractVector{<:AbstractString},
)
    record.schema_id == schema || _portable_emt_fail(
        :task_state_schema_mismatch,
        "portable task state expected schema $schema but found $(record.schema_id)",
    )
    fields = Dict{String,Any}(record.fields)
    Set(keys(fields)) == Set(String.(required)) || _portable_emt_fail(
        :task_state_field_mismatch,
        "portable task state $schema has missing or unknown fields",
    )
    return fields
end

function restore_portable_emt_task_state(
    template::Base.RefValue{T},
    value::PortableSnapshotRecord,
) where {T}
    fields = _portable_emt_scheduler_record_fields(
        value,
        "aimora.emt.task_ref.v1",
        ["value"],
    )
    restored = Ref{T}()
    restored[] = restore_portable_emt_task_state(template[], fields["value"])
    return restored
end

function restore_portable_emt_task_state(template, value)
    _portable_emt_fail(
        :unsupported_task_state,
        "task state $(typeof(template)) requires an explicit restore_portable_emt_task_state method",
    )
end

function _portable_emt_task_program_signature(value::AbstractString)
    signature = String(value)
    occursin(r"^[0-9a-f]{64}$", signature) || _portable_emt_fail(
        :invalid_task_program_signature,
        "portable task-program identity must be a lowercase SHA-256 digest",
    )
    return signature
end

function _portable_emt_task_time_value(value, identity::AbstractString)
    value isa Float64 && return value
    value == _PORTABLE_EMT_POSITIVE_INFINITY && return Inf
    value == _PORTABLE_EMT_NEGATIVE_INFINITY && return -Inf
    _portable_emt_fail(
        :task_state_type_mismatch,
        "portable task time field $identity has the wrong type",
    )
end

function _portable_emt_exact_occurrence_record(occurrence::SampledTaskOccurrence)
    return PortableSnapshotRecord(
        "aimora.emt.exact_task_occurrence.v1",
        Pair{String,Any}[
            "execution_index" => occurrence.execution_index,
            "name" => String(occurrence.name),
            "priority" => occurrence.priority,
            "tick" => occurrence.tick,
            "time_s" => occurrence.time_s,
        ],
    )
end

function _portable_emt_exact_occurrence(record::PortableSnapshotRecord)
    fields = _portable_emt_scheduler_record_fields(
        record,
        "aimora.emt.exact_task_occurrence.v1",
        ["execution_index", "name", "priority", "tick", "time_s"],
    )
    return SampledTaskOccurrence(
        Symbol(fields["name"]),
        Float64(fields["time_s"]),
        Int(fields["tick"]),
        Int(fields["priority"]),
        Int(fields["execution_index"]),
    )
end

function _portable_emt_pending_control_record(value)
    return PortableSnapshotRecord(
        "aimora.emt.pending_control_value.v1",
        Pair{String,Any}[
            "release_tick" => value.release_tick,
            "sample_index" => value.sample_index,
            "sample_tick" => value.sample_tick,
            "value" => portable_emt_task_state(value.value),
        ],
    )
end

function _portable_emt_control_sample_record(sample)
    return PortableSnapshotRecord(
        "aimora.emt.control_sample.v1",
        Pair{String,Any}[
            "computed_output" => portable_emt_task_state(sample.computed_output),
            "name" => String(sample.name),
            "release_tick" => sample.release_tick,
            "sample_index" => sample.sample_index,
            "sample_tick" => sample.sample_tick,
            "sample_time_s" => sample.sample_time_s,
        ],
    )
end

function _portable_emt_control_write_record(write)
    return PortableSnapshotRecord(
        "aimora.emt.control_write.v1",
        Pair{String,Any}[
            "name" => String(write.name),
            "output" => portable_emt_task_state(write.output),
            "sample_index" => write.sample_index,
            "sample_tick" => write.sample_tick,
            "write_tick" => write.write_tick,
            "write_time_s" => write.write_time_s,
        ],
    )
end

function _portable_emt_pwm_occurrence_record(occurrence)
    return PortableSnapshotRecord(
        "aimora.emt.pwm_edge_occurrence.v1",
        Pair{String,Any}[
            "cycle_index" => occurrence.cycle_index,
            "high" => occurrence.high,
            "name" => String(occurrence.name),
            "quantized_duty" => occurrence.quantized_duty,
            "requested_duty" => occurrence.requested_duty,
            "tick" => occurrence.tick,
            "time_s" => occurrence.time_s,
        ],
    )
end

function _portable_emt_exact_task_record(task::ExactSampledTask)
    return PortableSnapshotRecord(
        "aimora.emt.exact_action_task.v1",
        Pair{String,Any}[
            "execution_count" => task.execution_count,
            "last_execution_tick" => task.last_execution_tick,
            "name" => String(task.name),
            "next_tick" => task.next_tick,
            "period_ticks" => task.period_ticks,
            "power_history_invalidating" => task.power_history_invalidating,
            "priority" => task.priority,
            "tick_s" => task.tick_s,
        ],
    )
end

function _portable_emt_exact_task_record(task::ExactSampledControlTask)
    return PortableSnapshotRecord(
        "aimora.emt.exact_control_task.v1",
        Pair{String,Any}[
            "boundary_count" => task.boundary_count,
            "computational_delay_ticks" => task.computational_delay_ticks,
            "held_output" => portable_emt_task_state(task.held_output),
            "last_sample_tick" => task.last_sample_tick,
            "last_write_tick" => task.last_write_tick,
            "name" => String(task.name),
            "next_tick" => task.next_tick,
            "pending" => _portable_emt_pending_control_record.(task.pending),
            "period_ticks" => task.period_ticks,
            "power_history_invalidating" => task.power_history_invalidating,
            "priority" => task.priority,
            "sample_count" => task.sample_count,
            "samples" => _portable_emt_control_sample_record.(task.samples),
            "tick_s" => task.tick_s,
            "write_count" => task.write_count,
            "writes" => _portable_emt_control_write_record.(task.writes),
        ],
    )
end

function _portable_emt_exact_task_record(task::ExactPWMTask)
    return PortableSnapshotRecord(
        "aimora.emt.exact_pwm_task.v1",
        Pair{String,Any}[
            "boundary_count" => task.boundary_count,
            "carrier_period_ticks" => task.carrier_period_ticks,
            "cycle_count" => task.cycle_count,
            "edge_count" => task.edge_count,
            "gate_high" => task.gate_high,
            "last_edge_tick" => task.last_edge_tick,
            "last_quantized_duty" => _portable_emt_time_value(
                task.last_quantized_duty,
                "scheduler PWM quantized duty",
            ),
            "last_requested_duty" => _portable_emt_time_value(
                task.last_requested_duty,
                "scheduler PWM requested duty",
            ),
            "name" => String(task.name),
            "next_carrier_tick" => task.next_carrier_tick,
            "next_edge_tick" => task.next_edge_tick,
            "occurrences" => _portable_emt_pwm_occurrence_record.(task.occurrences),
            "power_history_invalidating" => task.power_history_invalidating,
            "priority" => task.priority,
            "tick_s" => task.tick_s,
        ],
    )
end

function _portable_emt_logical_time_value(value::EMTTaskPlatform.EMTLogicalTime)
    return value.numerator // value.denominator
end

function _portable_emt_general_occurrence_record(occurrence::EMTTaskPlatform.EMTTaskOccurrence)
    return PortableSnapshotRecord(
        "aimora.emt.general_task_occurrence.v1",
        Pair{String,Any}[
            "activation_index" => occurrence.activation_index,
            "exact_instant" => _portable_emt_logical_time_value(occurrence.exact_instant),
            "execution_index" => occurrence.execution_index,
            "family" => UInt8(occurrence.family),
            "instant_s" => occurrence.instant_s,
            "priority" => occurrence.priority,
            "release_index" => occurrence.release_index,
            "sample_index" => occurrence.sample_index,
            "stage" => UInt8(occurrence.stage),
            "task" => occurrence.task,
        ],
    )
end

function _portable_emt_general_occurrence(record::PortableSnapshotRecord)
    fields = _portable_emt_scheduler_record_fields(
        record,
        "aimora.emt.general_task_occurrence.v1",
        [
            "activation_index", "exact_instant", "execution_index", "family",
            "instant_s", "priority", "release_index", "sample_index", "stage", "task",
        ],
    )
    instant = fields["exact_instant"]
    instant isa Rational || _portable_emt_fail(
        :task_state_type_mismatch,
        "portable general-task occurrence has a non-rational instant",
    )
    return EMTTaskPlatform.EMTTaskOccurrence(
        String(fields["task"]),
        EMTTaskPlatform.EMTTaskFamily(UInt8(fields["family"])),
        EMTTaskPlatform.EMTLogicalTime(numerator(instant), denominator(instant)),
        Float64(fields["instant_s"]),
        EMTTaskPlatform.EMTTaskStage(UInt8(fields["stage"])),
        Int(fields["priority"]),
        Int(fields["activation_index"]),
        Int(fields["sample_index"]),
        Int(fields["release_index"]),
        Int(fields["execution_index"]),
    )
end

function _portable_emt_general_pending_record(value)
    return PortableSnapshotRecord(
        "aimora.emt.general_pending_value.v1",
        Pair{String,Any}[
            "release_tick" => value.release_tick,
            "sample_index" => value.sample_index,
            "sample_tick" => value.sample_tick,
            "value" => portable_emt_task_state(value.value),
        ],
    )
end

function _portable_emt_general_task_record(task::GeneralEMTTask)
    entry = task.plan_entry
    return PortableSnapshotRecord(
        "aimora.emt.general_task.v1",
        Pair{String,Any}[
            "activation_count" => task.activation_count,
            "boundary_count" => task.boundary_count,
            "delay_ticks" => entry.delay_ticks,
            "execution_rank" => entry.execution_rank,
            "family" => UInt8(entry.spec.family),
            "first_activation_tick" => entry.first_activation_tick,
            "held_output" => portable_emt_task_state(task.held_output),
            "last_sample_tick" => task.last_sample_tick,
            "last_write_tick" => task.last_write_tick,
            "maximum_pending_depth" => task.maximum_pending_depth,
            "name" => entry.spec.name,
            "next_activation_tick" => task.next_activation_tick,
            "pending" => _portable_emt_general_pending_record.(task.pending),
            "pending_head" => task.pending_head,
            "period_ticks" => entry.period_ticks,
            "registration_index" => entry.registration_index,
            "release_count" => task.release_count,
            "sample_count" => task.sample_count,
            "state" => portable_emt_task_state(task.state),
            "write_count" => task.write_count,
        ],
    )
end

function _portable_emt_scheduler_fields(
    scheduler::ExactSampledTaskScheduler,
    program_signature::String,
    adapter::Bool,
)
    return PortableSnapshotStateField[
        _portable_emt_state_field("scheduler.adapter", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", adapter),
        _portable_emt_state_field("scheduler.kind", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", "exact_sampled"),
        _portable_emt_state_field("scheduler.last_run_power_history_invalidating", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :discrete, "1", scheduler.last_run_power_history_invalidating),
        _portable_emt_state_field("scheduler.occurrences", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :output, "1", _portable_emt_exact_occurrence_record.(scheduler.occurrences); axes = ["occurrence"]),
        _portable_emt_state_field("scheduler.origin", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "s", scheduler.origin_s),
        _portable_emt_state_field("scheduler.program_signature", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", program_signature),
        _portable_emt_state_field("scheduler.retain_occurrences", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", scheduler.retain_occurrences),
        _portable_emt_state_field("scheduler.tasks", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", _portable_emt_exact_task_record.(scheduler.tasks); axes = ["task"]),
        _portable_emt_state_field("scheduler.tick", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "s", scheduler.tick_s),
    ]
end

function _portable_emt_scheduler_fields(
    scheduler::GeneralEMTTaskScheduler,
    program_signature::String,
    adapter::Bool,
)
    scheduler.terminal_failure === nothing || _portable_emt_fail(
        :terminal_task_scheduler,
        "terminally failed general task scheduler cannot be captured",
    )
    return PortableSnapshotStateField[
        _portable_emt_state_field("scheduler.adapter", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", adapter),
        _portable_emt_state_field("scheduler.current_tick", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", scheduler.current_tick),
        _portable_emt_state_field("scheduler.due_indices", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", _portable_emt_integer_array(scheduler.due_indices, "1", ["task"]); axes = ["task"]),
        _portable_emt_state_field("scheduler.execution_count", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :output, "1", scheduler.execution_count),
        _portable_emt_state_field("scheduler.kind", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", "general"),
        _portable_emt_state_field("scheduler.last_accepted_tick", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", scheduler.last_accepted_tick),
        _portable_emt_state_field("scheduler.last_run_effects", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :discrete, "1", copy(scheduler.last_run_effects); axes = ["effect"]),
        _portable_emt_state_field("scheduler.occurrences", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :output, "1", _portable_emt_general_occurrence_record.(scheduler.occurrences); axes = ["occurrence"]),
        _portable_emt_state_field("scheduler.plan_signature", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", scheduler.plan.signature_sha256),
        _portable_emt_state_field("scheduler.program_signature", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", program_signature),
        _portable_emt_state_field("scheduler.retain_occurrences", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :checkpoint, "1", scheduler.retain_occurrences),
        _portable_emt_state_field("scheduler.staged_occurrences", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", _portable_emt_general_occurrence_record.(scheduler.staged_occurrences); axes = ["occurrence"]),
        _portable_emt_state_field("scheduler.tasks", _PORTABLE_EMT_TASK_SCHEDULER_OWNER, :scheduler, "1", _portable_emt_general_task_record.(scheduler.tasks); axes = ["task"]),
    ]
end

"""Capture complete accepted exact/general task scheduler state without serializing callbacks."""
function portable_emt_task_scheduler_state_inventory(
    scheduler::Union{ExactSampledTaskScheduler,GeneralEMTTaskScheduler};
    program_signature_sha256::AbstractString,
)
    program_signature = _portable_emt_task_program_signature(program_signature_sha256)
    return PortableSnapshotStateInventory(
        _portable_emt_scheduler_fields(scheduler, program_signature, false),
    )
end

function portable_emt_task_scheduler_state_inventory(
    adapter::ExactSampledTaskCompatibilityAdapter;
    program_signature_sha256::AbstractString,
)
    program_signature = _portable_emt_task_program_signature(program_signature_sha256)
    return PortableSnapshotStateInventory(
        _portable_emt_scheduler_fields(adapter.scheduler, program_signature, true),
    )
end

function _portable_emt_scheduler_inventory_fields(inventory::PortableSnapshotStateInventory)
    return Dict(field.identity => field for field in inventory.fields)
end

function _portable_emt_scheduler_value(fields, identity, family, unit)
    return _portable_emt_inventory_field(fields, identity, family, unit).value
end

function _portable_emt_scheduler_records(fields, identity, family)
    value = _portable_emt_scheduler_value(fields, identity, family, "1")
    value isa AbstractVector && all(item -> item isa PortableSnapshotRecord, value) ||
        _portable_emt_fail(
            :task_state_type_mismatch,
            "portable scheduler field $identity is not a record sequence",
        )
    return PortableSnapshotRecord[item for item in value]
end

function _portable_emt_validate_exact_task_identity(task, values)
    String(task.name) == values["name"] || _portable_emt_fail(
        :task_plan_mismatch,
        "portable exact task identity changed",
    )
    task.tick_s == Float64(values["tick_s"]) || _portable_emt_fail(:task_plan_mismatch, "portable exact task tick changed")
    task.priority == Int(values["priority"]) || _portable_emt_fail(:task_plan_mismatch, "portable exact task priority changed")
    task.power_history_invalidating == values["power_history_invalidating"] ||
        _portable_emt_fail(:task_plan_mismatch, "portable exact task effects changed")
    return task
end

function _restore_portable_emt_exact_task!(task::ExactSampledTask, record)
    required = [
        "execution_count", "last_execution_tick", "name", "next_tick", "period_ticks",
        "power_history_invalidating", "priority", "tick_s",
    ]
    values = _portable_emt_scheduler_record_fields(record, "aimora.emt.exact_action_task.v1", required)
    _portable_emt_validate_exact_task_identity(task, values)
    task.period_ticks == Int(values["period_ticks"]) || _portable_emt_fail(:task_plan_mismatch, "portable exact task period changed")
    task.next_tick = Int(values["next_tick"])
    task.execution_count = Int(values["execution_count"])
    task.last_execution_tick = Int(values["last_execution_tick"])
    return task
end

function _restore_portable_emt_exact_task!(task::ExactSampledControlTask, record)
    required = [
        "boundary_count", "computational_delay_ticks", "held_output", "last_sample_tick",
        "last_write_tick", "name", "next_tick", "pending", "period_ticks",
        "power_history_invalidating", "priority", "sample_count", "samples", "tick_s",
        "write_count", "writes",
    ]
    values = _portable_emt_scheduler_record_fields(record, "aimora.emt.exact_control_task.v1", required)
    _portable_emt_validate_exact_task_identity(task, values)
    task.period_ticks == Int(values["period_ticks"]) || _portable_emt_fail(:task_plan_mismatch, "portable control period changed")
    task.computational_delay_ticks == Int(values["computational_delay_ticks"]) || _portable_emt_fail(:task_plan_mismatch, "portable control delay changed")
    held = restore_portable_emt_task_state(task.held_output, values["held_output"])
    pending_records = values["pending"]
    sample_records = values["samples"]
    write_records = values["writes"]
    all(records -> records isa AbstractVector && all(record -> record isa PortableSnapshotRecord, records), (pending_records, sample_records, write_records)) ||
        _portable_emt_fail(:task_state_type_mismatch, "portable control histories are not record sequences")
    pending = eltype(task.pending)[]
    for record in pending_records
        fields = _portable_emt_scheduler_record_fields(record, "aimora.emt.pending_control_value.v1", ["release_tick", "sample_index", "sample_tick", "value"])
        push!(pending, eltype(task.pending)(
            Int(fields["release_tick"]),
            Int(fields["sample_tick"]),
            Int(fields["sample_index"]),
            restore_portable_emt_task_state(task.held_output, fields["value"]),
        ))
    end
    samples = eltype(task.samples)[]
    for record in sample_records
        fields = _portable_emt_scheduler_record_fields(record, "aimora.emt.control_sample.v1", ["computed_output", "name", "release_tick", "sample_index", "sample_tick", "sample_time_s"])
        Symbol(fields["name"]) == task.name || _portable_emt_fail(:task_plan_mismatch, "portable control sample owner changed")
        push!(samples, eltype(task.samples)(
            task.name,
            Float64(fields["sample_time_s"]),
            Int(fields["sample_tick"]),
            Int(fields["release_tick"]),
            Int(fields["sample_index"]),
            restore_portable_emt_task_state(task.held_output, fields["computed_output"]),
        ))
    end
    writes = eltype(task.writes)[]
    for record in write_records
        fields = _portable_emt_scheduler_record_fields(record, "aimora.emt.control_write.v1", ["name", "output", "sample_index", "sample_tick", "write_tick", "write_time_s"])
        Symbol(fields["name"]) == task.name || _portable_emt_fail(:task_plan_mismatch, "portable control write owner changed")
        push!(writes, eltype(task.writes)(
            task.name,
            Float64(fields["write_time_s"]),
            Int(fields["write_tick"]),
            Int(fields["sample_tick"]),
            Int(fields["sample_index"]),
            restore_portable_emt_task_state(task.held_output, fields["output"]),
        ))
    end
    task.next_tick = Int(values["next_tick"])
    task.sample_count = Int(values["sample_count"])
    task.write_count = Int(values["write_count"])
    task.boundary_count = Int(values["boundary_count"])
    task.last_sample_tick = Int(values["last_sample_tick"])
    task.last_write_tick = Int(values["last_write_tick"])
    task.held_output = held
    empty!(task.pending); append!(task.pending, pending)
    empty!(task.samples); append!(task.samples, samples)
    empty!(task.writes); append!(task.writes, writes)
    return task
end

function _restore_portable_emt_exact_task!(task::ExactPWMTask, record)
    required = [
        "boundary_count", "carrier_period_ticks", "cycle_count", "edge_count", "gate_high",
        "last_edge_tick", "last_quantized_duty", "last_requested_duty", "name",
        "next_carrier_tick", "next_edge_tick", "occurrences", "power_history_invalidating",
        "priority", "tick_s",
    ]
    values = _portable_emt_scheduler_record_fields(record, "aimora.emt.exact_pwm_task.v1", required)
    _portable_emt_validate_exact_task_identity(task, values)
    task.carrier_period_ticks == Int(values["carrier_period_ticks"]) || _portable_emt_fail(:task_plan_mismatch, "portable PWM period changed")
    occurrence_records = values["occurrences"]
    occurrence_records isa AbstractVector && all(item -> item isa PortableSnapshotRecord, occurrence_records) ||
        _portable_emt_fail(:task_state_type_mismatch, "portable PWM occurrences are not records")
    occurrences = eltype(task.occurrences)[]
    for record in occurrence_records
        fields = _portable_emt_scheduler_record_fields(record, "aimora.emt.pwm_edge_occurrence.v1", ["cycle_index", "high", "name", "quantized_duty", "requested_duty", "tick", "time_s"])
        Symbol(fields["name"]) == task.name || _portable_emt_fail(:task_plan_mismatch, "portable PWM occurrence owner changed")
        push!(occurrences, eltype(task.occurrences)(
            task.name,
            Float64(fields["time_s"]),
            Int(fields["tick"]),
            Bool(fields["high"]),
            Int(fields["cycle_index"]),
            Float64(fields["requested_duty"]),
            Float64(fields["quantized_duty"]),
        ))
    end
    task.next_carrier_tick = Int(values["next_carrier_tick"])
    task.next_edge_tick = Int(values["next_edge_tick"])
    task.gate_high = Bool(values["gate_high"])
    task.cycle_count = Int(values["cycle_count"])
    task.edge_count = Int(values["edge_count"])
    task.boundary_count = Int(values["boundary_count"])
    task.last_edge_tick = Int(values["last_edge_tick"])
    task.last_requested_duty = _portable_emt_task_time_value(
        values["last_requested_duty"],
        "scheduler PWM requested duty",
    )
    task.last_quantized_duty = _portable_emt_task_time_value(
        values["last_quantized_duty"],
        "scheduler PWM quantized duty",
    )
    empty!(task.occurrences); append!(task.occurrences, occurrences)
    return task
end

function _restore_portable_emt_exact_scheduler!(scheduler, fields)
    _portable_emt_scheduler_value(fields, "scheduler.kind", :checkpoint, "1") == "exact_sampled" ||
        _portable_emt_fail(:task_scheduler_kind_mismatch, "portable task scheduler kind changed")
    scheduler.tick_s == _portable_emt_scalar(fields, "scheduler.tick", :checkpoint, "s", Float64) ||
        _portable_emt_fail(:task_plan_mismatch, "portable scheduler tick changed")
    scheduler.origin_s == _portable_emt_scalar(fields, "scheduler.origin", :checkpoint, "s", Float64) ||
        _portable_emt_fail(:task_plan_mismatch, "portable scheduler origin changed")
    scheduler.retain_occurrences == _portable_emt_scalar(fields, "scheduler.retain_occurrences", :checkpoint, "1", Bool) ||
        _portable_emt_fail(:task_plan_mismatch, "portable scheduler occurrence policy changed")
    task_records = _portable_emt_scheduler_records(fields, "scheduler.tasks", :scheduler)
    length(task_records) == length(scheduler.tasks) || _portable_emt_fail(:task_plan_mismatch, "portable scheduler task count changed")
    for index in eachindex(scheduler.tasks)
        _restore_portable_emt_exact_task!(scheduler.tasks[index], task_records[index])
    end
    occurrence_records = _portable_emt_scheduler_records(fields, "scheduler.occurrences", :output)
    occurrences = _portable_emt_exact_occurrence.(occurrence_records)
    empty!(scheduler.occurrences); append!(scheduler.occurrences, occurrences)
    scheduler.last_run_power_history_invalidating = _portable_emt_scalar(
        fields,
        "scheduler.last_run_power_history_invalidating",
        :discrete,
        "1",
        Bool,
    )
    return scheduler
end

function _restore_portable_emt_general_task!(task, record)
    required = [
        "activation_count", "boundary_count", "delay_ticks", "execution_rank", "family",
        "first_activation_tick", "held_output", "last_sample_tick", "last_write_tick",
        "maximum_pending_depth", "name", "next_activation_tick", "pending", "pending_head",
        "period_ticks", "registration_index", "release_count", "sample_count", "state",
        "write_count",
    ]
    values = _portable_emt_scheduler_record_fields(record, "aimora.emt.general_task.v1", required)
    entry = task.plan_entry
    entry.spec.name == values["name"] && UInt8(entry.spec.family) == UInt8(values["family"]) &&
        entry.registration_index == Int(values["registration_index"]) &&
        entry.first_activation_tick == Int64(values["first_activation_tick"]) &&
        entry.period_ticks == Int64(values["period_ticks"]) &&
        entry.delay_ticks == Int64(values["delay_ticks"]) &&
        entry.execution_rank == Int(values["execution_rank"]) ||
        _portable_emt_fail(:task_plan_mismatch, "portable general task plan entry changed")
    state = restore_portable_emt_task_state(task.state, values["state"])
    held = restore_portable_emt_task_state(task.held_output, values["held_output"])
    pending_records = values["pending"]
    pending_records isa AbstractVector && all(item -> item isa PortableSnapshotRecord, pending_records) ||
        _portable_emt_fail(:task_state_type_mismatch, "portable general pending values are not records")
    pending = eltype(task.pending)[]
    for record in pending_records
        fields = _portable_emt_scheduler_record_fields(record, "aimora.emt.general_pending_value.v1", ["release_tick", "sample_index", "sample_tick", "value"])
        push!(pending, eltype(task.pending)(
            Int64(fields["release_tick"]),
            Int64(fields["sample_tick"]),
            Int(fields["sample_index"]),
            restore_portable_emt_task_state(task.held_output, fields["value"]),
        ))
    end
    task.state = state
    task.next_activation_tick = Int64(values["next_activation_tick"])
    task.activation_count = Int(values["activation_count"])
    task.sample_count = Int(values["sample_count"])
    task.release_count = Int(values["release_count"])
    task.write_count = Int(values["write_count"])
    task.boundary_count = Int(values["boundary_count"])
    task.last_sample_tick = Int64(values["last_sample_tick"])
    task.last_write_tick = Int64(values["last_write_tick"])
    task.held_output = held
    empty!(task.pending); append!(task.pending, pending)
    task.pending_head = Int(values["pending_head"])
    task.maximum_pending_depth = Int(values["maximum_pending_depth"])
    1 <= task.pending_head <= length(task.pending) + 1 || _portable_emt_fail(
        :task_state_shape_mismatch,
        "portable general pending cursor is outside its queue",
    )
    return task
end

function _restore_portable_emt_general_scheduler!(scheduler, fields)
    _portable_emt_scheduler_value(fields, "scheduler.kind", :checkpoint, "1") == "general" ||
        _portable_emt_fail(:task_scheduler_kind_mismatch, "portable task scheduler kind changed")
    _portable_emt_scheduler_value(fields, "scheduler.plan_signature", :checkpoint, "1") == scheduler.plan.signature_sha256 ||
        _portable_emt_fail(:task_plan_mismatch, "portable general task plan changed")
    scheduler.retain_occurrences == _portable_emt_scalar(fields, "scheduler.retain_occurrences", :checkpoint, "1", Bool) ||
        _portable_emt_fail(:task_plan_mismatch, "portable scheduler occurrence policy changed")
    task_records = _portable_emt_scheduler_records(fields, "scheduler.tasks", :scheduler)
    length(task_records) == length(scheduler.tasks) || _portable_emt_fail(:task_plan_mismatch, "portable scheduler task count changed")
    for index in eachindex(scheduler.tasks)
        _restore_portable_emt_general_task!(scheduler.tasks[index], task_records[index])
    end
    occurrences = _portable_emt_general_occurrence.(_portable_emt_scheduler_records(fields, "scheduler.occurrences", :output))
    staged = _portable_emt_general_occurrence.(_portable_emt_scheduler_records(fields, "scheduler.staged_occurrences", :scheduler))
    due = _portable_emt_array_values(fields, "scheduler.due_indices", :scheduler, "1", ["task"], Int64)
    all(index -> 1 <= index <= length(scheduler.tasks), due) || _portable_emt_fail(:task_state_shape_mismatch, "portable due-task index is outside the task plan")
    effects = _portable_emt_scheduler_value(fields, "scheduler.last_run_effects", :discrete, "1")
    effects isa AbstractVector && all(value -> value isa Bool, effects) && length(effects) == length(scheduler.last_run_effects) ||
        _portable_emt_fail(:task_state_shape_mismatch, "portable task effect state changed shape")
    scheduler.current_tick = Int64(_portable_emt_integer(fields, "scheduler.current_tick", :scheduler))
    scheduler.last_accepted_tick = Int64(_portable_emt_integer(fields, "scheduler.last_accepted_tick", :scheduler))
    scheduler.execution_count = _portable_emt_integer(fields, "scheduler.execution_count", :output)
    empty!(scheduler.occurrences); append!(scheduler.occurrences, occurrences)
    empty!(scheduler.due_indices); append!(scheduler.due_indices, Int.(due))
    empty!(scheduler.staged_occurrences); append!(scheduler.staged_occurrences, staged)
    copyto!(scheduler.last_run_effects, Bool[value for value in effects])
    scheduler.terminal_failure = nothing
    return scheduler
end

function _restore_portable_emt_scheduler_inventory!(
    scheduler,
    inventory,
    program_signature,
    adapter,
)
    fields = _portable_emt_scheduler_inventory_fields(inventory)
    _portable_emt_scheduler_value(fields, "scheduler.program_signature", :checkpoint, "1") == program_signature ||
        _portable_emt_fail(:task_program_mismatch, "portable task callback program identity changed")
    _portable_emt_scalar(fields, "scheduler.adapter", :checkpoint, "1", Bool) == adapter ||
        _portable_emt_fail(:task_scheduler_kind_mismatch, "portable task compatibility-adapter boundary changed")
    backup = deepcopy(scheduler)
    try
        if scheduler isa ExactSampledTaskScheduler
            _restore_portable_emt_exact_scheduler!(scheduler, fields)
        else
            _restore_portable_emt_general_scheduler!(scheduler, fields)
        end
        restored = PortableSnapshotStateInventory(
            _portable_emt_scheduler_fields(scheduler, program_signature, adapter),
        )
        restored.signature_sha256 == inventory.signature_sha256 || _portable_emt_fail(
            :task_state_reconstruction,
            "portable task scheduler state changed during reconstruction",
        )
        return scheduler
    catch
        restore_timestep_state!(scheduler, backup)
        rethrow()
    end
end

"""Restore scheduler state only after its immutable plan and callback program are independently reconstructed."""
function restore_portable_emt_task_scheduler_state_inventory!(
    scheduler::Union{ExactSampledTaskScheduler,GeneralEMTTaskScheduler},
    inventory::PortableSnapshotStateInventory;
    program_signature_sha256::AbstractString,
)
    signature = _portable_emt_task_program_signature(program_signature_sha256)
    return _restore_portable_emt_scheduler_inventory!(scheduler, inventory, signature, false)
end

function restore_portable_emt_task_scheduler_state_inventory!(
    adapter::ExactSampledTaskCompatibilityAdapter,
    inventory::PortableSnapshotStateInventory;
    program_signature_sha256::AbstractString,
)
    signature = _portable_emt_task_program_signature(program_signature_sha256)
    _restore_portable_emt_scheduler_inventory!(adapter.scheduler, inventory, signature, true)
    return adapter
end

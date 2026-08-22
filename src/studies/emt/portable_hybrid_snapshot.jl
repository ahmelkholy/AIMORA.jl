const _PORTABLE_EMT_HYBRID_OWNER = "emt.hybrid_execution"
const _PORTABLE_EMT_HYBRID_NOT_INITIALIZED = "not_initialized"

function _portable_emt_hybrid_program_signature(
    value::AbstractString,
    label::AbstractString,
)
    signature = String(value)
    occursin(r"^[0-9a-f]{64}$", signature) || _portable_emt_fail(
        :invalid_hybrid_program_signature,
        "portable $label identity must be a lowercase SHA-256 digest",
    )
    return signature
end

function _portable_emt_hybrid_surface_record(surface::HybridEventSurface)
    return PortableSnapshotRecord(
        "aimora.emt.hybrid_surface.v1",
        Pair{String,Any}[
            "candidate_is_event" => surface.candidate_is_event,
            "direction" => Int8(surface.direction),
            "name" => String(surface.name),
            "priority" => surface.priority,
            "repeatable" => surface.repeatable,
            "topology_invalidating" => surface.topology_invalidating,
        ],
    )
end

function _portable_emt_hybrid_occurrence_record(occurrence::HybridEventOccurrence)
    return PortableSnapshotRecord(
        "aimora.emt.hybrid_occurrence.v1",
        Pair{String,Any}[
            "name" => String(occurrence.name),
            "priority" => occurrence.priority,
            "root_bracket_width_s" => occurrence.root_bracket_width_s,
            "root_iteration_count" => occurrence.root_iteration_count,
            "time_s" => occurrence.time_s,
            "topology_invalidating" => occurrence.topology_invalidating,
            "value" => occurrence.value,
        ],
    )
end

function _portable_emt_hybrid_occurrence(record::PortableSnapshotRecord)
    fields = _portable_emt_scheduler_record_fields(
        record,
        "aimora.emt.hybrid_occurrence.v1",
        [
            "name", "priority", "root_bracket_width_s", "root_iteration_count",
            "time_s", "topology_invalidating", "value",
        ],
    )
    return HybridEventOccurrence(
        Symbol(fields["name"]),
        Float64(fields["time_s"]),
        Float64(fields["value"]),
        Int(fields["priority"]),
        Bool(fields["topology_invalidating"]),
        Int(fields["root_iteration_count"]),
        Float64(fields["root_bracket_width_s"]),
    )
end

function _portable_emt_hybrid_endpoint_value(integrator::EMTHybridStepIntegrator)
    integrator.initialized && return _portable_emt_time_value(
        integrator.active_global_endpoint_s,
        "hybrid active global endpoint",
    )
    isnan(integrator.active_global_endpoint_s) || _portable_emt_fail(
        :hybrid_state_mismatch,
        "uninitialized hybrid execution has an active endpoint",
    )
    return _PORTABLE_EMT_HYBRID_NOT_INITIALIZED
end

function _portable_emt_hybrid_endpoint_from_value(value)
    value == _PORTABLE_EMT_HYBRID_NOT_INITIALIZED && return NaN
    return _portable_emt_task_time_value(value, "hybrid active global endpoint")
end

function _portable_emt_hybrid_accepted_step(context)
    return min(max(context.step_index - 1, 0), context.step_count)
end

function _portable_emt_validate_hybrid_boundary(
    integrator::EMTHybridStepIntegrator;
    allow_reconstruction_candidate::Bool=false,
)
    transaction = emt_step_transaction_status(integrator.transaction)
    transaction.active && _portable_emt_fail(
        :active_transaction,
        "portable hybrid capture requires an inactive EMT transaction",
    )
    integrator.last_failure === nothing || _portable_emt_fail(
        :terminal_hybrid_execution,
        "terminally failed hybrid execution cannot be captured",
    )
    runtime = _check_prepared_runtime_aliases(integrator.workspace.runtime)
    integrator.callback_owner.runtime === runtime || _portable_emt_fail(
        :hybrid_owner_alias,
        "hybrid callback owner is detached from its workspace runtime",
    )
    integrator.transaction.runtime === runtime || _portable_emt_fail(
        :hybrid_owner_alias,
        "hybrid transaction owner is detached from its workspace runtime",
    )
    context = runtime.context
    if integrator.initialized
        context.dt_s == integrator.nominal_dt_s || _portable_emt_fail(
            :unsynchronized_capture,
            "portable hybrid capture is inside a provisional local interval",
        )
        accepted_step = _portable_emt_hybrid_accepted_step(context)
        accepted_time = accepted_step * integrator.nominal_dt_s
        next_time = min(context.step_index, context.step_count) * integrator.nominal_dt_s
        context.t_s == next_time || _portable_emt_fail(
            :unsynchronized_capture,
            "portable hybrid workspace cursor is outside its accepted fixed-step boundary",
        )
        integrator.active_global_endpoint_s == accepted_time || _portable_emt_fail(
            :unsynchronized_capture,
            "portable hybrid capture is not at its accepted global endpoint",
        )
    else
        valid_uninitialized_owner = integrator.workspace.ready ||
            allow_reconstruction_candidate &&
                integrator.workspace.execution_mode === :hybrid
        valid_uninitialized_owner || _portable_emt_fail(
            :hybrid_state_mismatch,
            "uninitialized hybrid execution has already mutated its workspace",
        )
    end
    return transaction
end

function _portable_emt_hybrid_fields(
    integrator::EMTHybridStepIntegrator,
    event_program_signature::String,
)
    transaction = _portable_emt_validate_hybrid_boundary(integrator)
    context = integrator.workspace.runtime.context
    return PortableSnapshotStateField[
        _portable_emt_state_field("hybrid.accepted_step", _PORTABLE_EMT_HYBRID_OWNER, :scheduler, "1", _portable_emt_hybrid_accepted_step(context)),
        _portable_emt_state_field("hybrid.accepted_interval_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.accepted_interval_count),
        _portable_emt_state_field("hybrid.active_global_endpoint", _PORTABLE_EMT_HYBRID_OWNER, :scheduler, "s", _portable_emt_hybrid_endpoint_value(integrator)),
        _portable_emt_state_field("hybrid.completed", _PORTABLE_EMT_HYBRID_OWNER, :discrete, "1", integrator.completed),
        _portable_emt_state_field("hybrid.completed_global_step_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.completed_global_step_count),
        _portable_emt_state_field("hybrid.evaluation_recorded", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.evaluation_recorded),
        _portable_emt_state_field("hybrid.event_program_signature", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", event_program_signature),
        _portable_emt_state_field("hybrid.initialized", _PORTABLE_EMT_HYBRID_OWNER, :discrete, "1", integrator.initialized),
        _portable_emt_state_field("hybrid.last_global_step_event_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.last_global_step_event_count),
        _portable_emt_state_field("hybrid.localized_root_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.localized_root_count),
        _portable_emt_state_field("hybrid.nominal_timestep", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "s", integrator.nominal_dt_s),
        _portable_emt_state_field("hybrid.occurrences", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", _portable_emt_hybrid_occurrence_record.(integrator.occurrences); axes = ["occurrence"]),
        _portable_emt_state_field("hybrid.policy.maximum_events_per_step", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", integrator.policy.maximum_events_per_step),
        _portable_emt_state_field("hybrid.policy.maximum_root_iterations", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", integrator.policy.maximum_root_iterations),
        _portable_emt_state_field("hybrid.policy.minimum_progress", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "s", integrator.policy.minimum_progress_s),
        _portable_emt_state_field("hybrid.policy.root_time_tolerance", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "s", integrator.policy.root_time_tolerance_s),
        _portable_emt_state_field("hybrid.policy.root_value_tolerance", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", integrator.policy.root_value_tolerance),
        _portable_emt_state_field("hybrid.policy.simultaneity_tolerance", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "s", integrator.policy.simultaneity_tolerance_s),
        _portable_emt_state_field("hybrid.provisional_interval_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.provisional_interval_count),
        _portable_emt_state_field("hybrid.surface_fired", _PORTABLE_EMT_HYBRID_OWNER, :discrete, "1", collect(integrator.surface_fired); axes = ["surface"]),
        _portable_emt_state_field("hybrid.surfaces", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", _portable_emt_hybrid_surface_record.(integrator.surfaces); axes = ["surface"]),
        _portable_emt_state_field("hybrid.topology_invalidation_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", integrator.topology_invalidation_count),
        _portable_emt_state_field("hybrid.transaction.capture_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", transaction.capture_count),
        _portable_emt_state_field("hybrid.transaction.commit_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", transaction.commit_count),
        _portable_emt_state_field("hybrid.transaction.restore_count", _PORTABLE_EMT_HYBRID_OWNER, :output, "1", transaction.restore_count),
        _portable_emt_state_field("hybrid.workspace_step_index", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "1", context.step_index),
        _portable_emt_state_field("hybrid.workspace_time", _PORTABLE_EMT_HYBRID_OWNER, :checkpoint, "s", context.t_s),
    ]
end

function _portable_emt_hybrid_scheduler(integrator::EMTHybridStepIntegrator)
    scheduler = integrator.scheduler
    return scheduler isa ExactSampledTaskCompatibilityAdapter ?
        (scheduler.scheduler, true) : (scheduler, false)
end

"""Capture the accepted hybrid coordinator plus its complete exact/general task scheduler."""
function portable_emt_hybrid_state_inventory(
    integrator::EMTHybridStepIntegrator;
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
)
    task_signature = _portable_emt_task_program_signature(
        task_program_signature_sha256,
    )
    event_signature = _portable_emt_hybrid_program_signature(
        event_program_signature_sha256,
        "hybrid-event program",
    )
    scheduler, adapter = _portable_emt_hybrid_scheduler(integrator)
    fields = _portable_emt_scheduler_fields(scheduler, task_signature, adapter)
    append!(fields, _portable_emt_hybrid_fields(integrator, event_signature))
    return PortableSnapshotStateInventory(fields)
end

function _portable_emt_hybrid_records(fields, identity, family)
    value = _portable_emt_scheduler_value(fields, identity, family, "1")
    value isa AbstractVector && all(item -> item isa PortableSnapshotRecord, value) ||
        _portable_emt_fail(
            :hybrid_state_type_mismatch,
            "portable hybrid field $identity is not a record sequence",
        )
    return PortableSnapshotRecord[item for item in value]
end

function _portable_emt_validate_hybrid_surface(
    surface::HybridEventSurface,
    record::PortableSnapshotRecord,
)
    values = _portable_emt_scheduler_record_fields(
        record,
        "aimora.emt.hybrid_surface.v1",
        [
            "candidate_is_event", "direction", "name", "priority", "repeatable",
            "topology_invalidating",
        ],
    )
    String(surface.name) == values["name"] &&
        Int8(surface.direction) == Int8(values["direction"]) &&
        surface.priority == Int(values["priority"]) &&
        surface.topology_invalidating == values["topology_invalidating"] &&
        surface.repeatable == values["repeatable"] &&
        surface.candidate_is_event == values["candidate_is_event"] ||
        _portable_emt_fail(
            :hybrid_event_program_mismatch,
            "portable hybrid event-surface identity changed",
        )
    return surface
end

function _portable_emt_validate_hybrid_policy(integrator, fields)
    policy = integrator.policy
    policy.root_time_tolerance_s == _portable_emt_scalar(fields, "hybrid.policy.root_time_tolerance", :checkpoint, "s", Float64) &&
        policy.root_value_tolerance == _portable_emt_scalar(fields, "hybrid.policy.root_value_tolerance", :checkpoint, "1", Float64) &&
        policy.simultaneity_tolerance_s == _portable_emt_scalar(fields, "hybrid.policy.simultaneity_tolerance", :checkpoint, "s", Float64) &&
        policy.minimum_progress_s == _portable_emt_scalar(fields, "hybrid.policy.minimum_progress", :checkpoint, "s", Float64) &&
        policy.maximum_root_iterations == _portable_emt_integer(fields, "hybrid.policy.maximum_root_iterations", :checkpoint) &&
        policy.maximum_events_per_step == _portable_emt_integer(fields, "hybrid.policy.maximum_events_per_step", :checkpoint) ||
        _portable_emt_fail(
            :hybrid_policy_mismatch,
            "portable hybrid numerical policy changed",
        )
    return policy
end

function _portable_emt_restore_hybrid_fields!(
    integrator,
    fields;
    allow_reconstruction_candidate::Bool=false,
)
    _portable_emt_validate_hybrid_boundary(
        integrator;
        allow_reconstruction_candidate,
    )
    integrator.nominal_dt_s == _portable_emt_scalar(fields, "hybrid.nominal_timestep", :checkpoint, "s", Float64) ||
        _portable_emt_fail(:hybrid_policy_mismatch, "portable hybrid nominal timestep changed")
    _portable_emt_validate_hybrid_policy(integrator, fields)
    context = integrator.workspace.runtime.context
    context.t_s == _portable_emt_scalar(fields, "hybrid.workspace_time", :checkpoint, "s", Float64) &&
        context.step_index == _portable_emt_integer(fields, "hybrid.workspace_step_index", :checkpoint) ||
        _portable_emt_fail(
            :hybrid_workspace_mismatch,
            "portable hybrid workspace is not at the represented boundary",
        )
    surface_records = _portable_emt_hybrid_records(fields, "hybrid.surfaces", :checkpoint)
    length(surface_records) == length(integrator.surfaces) || _portable_emt_fail(
        :hybrid_event_program_mismatch,
        "portable hybrid event-surface count changed",
    )
    foreach(_portable_emt_validate_hybrid_surface, integrator.surfaces, surface_records)
    fired = _portable_emt_scheduler_value(fields, "hybrid.surface_fired", :discrete, "1")
    fired isa AbstractVector && all(value -> value isa Bool, fired) &&
        length(fired) == length(integrator.surface_fired) || _portable_emt_fail(
            :hybrid_state_shape_mismatch,
            "portable hybrid fired-surface state changed shape",
        )
    occurrences = _portable_emt_hybrid_occurrence.(
        _portable_emt_hybrid_records(fields, "hybrid.occurrences", :output),
    )
    initialized = _portable_emt_scalar(fields, "hybrid.initialized", :discrete, "1", Bool)
    endpoint = _portable_emt_hybrid_endpoint_from_value(
        _portable_emt_scheduler_value(fields, "hybrid.active_global_endpoint", :scheduler, "s"),
    )
    initialized == !isnan(endpoint) || _portable_emt_fail(
        :hybrid_state_mismatch,
        "portable hybrid initialization and endpoint state disagree",
    )
    if initialized
        accepted_step = _portable_emt_integer(fields, "hybrid.accepted_step", :scheduler)
        accepted_step == _portable_emt_hybrid_accepted_step(context) ||
            _portable_emt_fail(
                :hybrid_workspace_mismatch,
                "portable hybrid accepted-step identity differs from its restored workspace",
            )
        endpoint == accepted_step * integrator.nominal_dt_s || _portable_emt_fail(
            :hybrid_workspace_mismatch,
            "portable hybrid endpoint differs from its restored accepted-step time",
        )
    else
        integrator.workspace.ready || _portable_emt_fail(
            :hybrid_workspace_mismatch,
            "uninitialized portable hybrid state requires a ready workspace",
        )
    end
    transaction = integrator.transaction.transaction
    copyto!(integrator.surface_fired, Bool[value for value in fired])
    empty!(integrator.occurrences)
    append!(integrator.occurrences, occurrences)
    integrator.active_global_endpoint_s = endpoint
    integrator.initialized = initialized
    integrator.completed = _portable_emt_scalar(fields, "hybrid.completed", :discrete, "1", Bool)
    integrator.evaluation_recorded = _portable_emt_scalar(fields, "hybrid.evaluation_recorded", :output, "1", Bool)
    integrator.accepted_interval_count = _portable_emt_integer(fields, "hybrid.accepted_interval_count", :output)
    integrator.provisional_interval_count = _portable_emt_integer(fields, "hybrid.provisional_interval_count", :output)
    integrator.localized_root_count = _portable_emt_integer(fields, "hybrid.localized_root_count", :output)
    integrator.topology_invalidation_count = _portable_emt_integer(fields, "hybrid.topology_invalidation_count", :output)
    integrator.completed_global_step_count = _portable_emt_integer(fields, "hybrid.completed_global_step_count", :output)
    integrator.last_global_step_event_count = _portable_emt_integer(fields, "hybrid.last_global_step_event_count", :output)
    integrator.last_failure = nothing
    transaction.active = false
    transaction.capture_count = _portable_emt_integer(fields, "hybrid.transaction.capture_count", :output)
    transaction.restore_count = _portable_emt_integer(fields, "hybrid.transaction.restore_count", :output)
    transaction.commit_count = _portable_emt_integer(fields, "hybrid.transaction.commit_count", :output)
    return integrator
end

"""Restore a hybrid coordinator only after its workspace, callbacks, plan, and surfaces are reconstructed."""
function _restore_portable_emt_hybrid_state_inventory!(
    integrator::EMTHybridStepIntegrator,
    inventory::PortableSnapshotStateInventory;
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
    isolated_reconstruction::Bool,
)
    task_signature = _portable_emt_task_program_signature(
        task_program_signature_sha256,
    )
    event_signature = _portable_emt_hybrid_program_signature(
        event_program_signature_sha256,
        "hybrid-event program",
    )
    fields = _portable_emt_scheduler_inventory_fields(inventory)
    _portable_emt_scheduler_value(fields, "hybrid.event_program_signature", :checkpoint, "1") == event_signature ||
        _portable_emt_fail(
            :hybrid_event_program_mismatch,
            "portable hybrid event callback program identity changed",
        )
    scheduler, adapter = _portable_emt_hybrid_scheduler(integrator)
    scheduler_inventory = PortableSnapshotStateInventory(PortableSnapshotStateField[
        field for field in inventory.fields if startswith(field.identity, "scheduler.")
    ])
    if isolated_reconstruction
        _restore_portable_emt_scheduler_inventory!(
            scheduler,
            scheduler_inventory,
            task_signature,
            adapter,
        )
        _portable_emt_restore_hybrid_fields!(
            integrator,
            fields;
            allow_reconstruction_candidate = true,
        )
        restored = portable_emt_hybrid_state_inventory(
            integrator;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        )
        restored.signature_sha256 == inventory.signature_sha256 || _portable_emt_fail(
            :hybrid_state_reconstruction,
            "portable hybrid state changed during reconstruction",
        )
        return integrator
    end
    backup_scheduler = PortableSnapshotStateInventory(
        _portable_emt_scheduler_fields(scheduler, task_signature, adapter),
    )
    backup_fields = _portable_emt_hybrid_fields(integrator, event_signature)
    backup_inventory = PortableSnapshotStateInventory(vcat(
        backup_scheduler.fields,
        backup_fields,
    ))
    try
        _restore_portable_emt_scheduler_inventory!(
            scheduler,
            scheduler_inventory,
            task_signature,
            adapter,
        )
        _portable_emt_restore_hybrid_fields!(integrator, fields)
        restored = portable_emt_hybrid_state_inventory(
            integrator;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        )
        restored.signature_sha256 == inventory.signature_sha256 || _portable_emt_fail(
            :hybrid_state_reconstruction,
            "portable hybrid state changed during reconstruction",
        )
        return integrator
    catch
        backup_scheduler_only = PortableSnapshotStateInventory(PortableSnapshotStateField[
            field for field in backup_inventory.fields if startswith(field.identity, "scheduler.")
        ])
        _restore_portable_emt_scheduler_inventory!(
            scheduler,
            backup_scheduler_only,
            task_signature,
            adapter,
        )
        backup_dictionary = _portable_emt_scheduler_inventory_fields(backup_inventory)
        _portable_emt_restore_hybrid_fields!(integrator, backup_dictionary)
        rethrow()
    end
end

function restore_portable_emt_hybrid_state_inventory!(
    integrator::EMTHybridStepIntegrator,
    inventory::PortableSnapshotStateInventory;
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
)
    return _restore_portable_emt_hybrid_state_inventory!(
        integrator,
        inventory;
        task_program_signature_sha256,
        event_program_signature_sha256,
        isolated_reconstruction = false,
    )
end

struct PortableEMTHybridRestoreResult{I}
    integrator::I
    descriptor::PortableSnapshotDescriptor
    public_state_signature_sha256::String
    backend_state_signature_sha256::String
    hybrid_state_signature_sha256::String
    reconstructed::Bool
end

"""Capture a coordinated accepted hybrid EMT boundary into canonical portable sections."""
function capture_portable_emt_snapshot(
    integrator::EMTHybridStepIntegrator;
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
    provenance::AbstractString,
    capabilities::AbstractVector{<:AbstractString}=String[
        "emt.events",
        "emt.fixed_step",
        "emt.portable_snapshot",
        "emt.tasks",
    ],
)
    hybrid_state = portable_emt_hybrid_state_inventory(
        integrator;
        task_program_signature_sha256,
        event_program_signature_sha256,
    )
    workspace = integrator.workspace
    public_state = _portable_emt_state_inventory(workspace)
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
    all(capability -> capability in metadata.capabilities, ("emt.events", "emt.tasks")) ||
        _portable_emt_fail(
            :capability_mismatch,
            "portable hybrid snapshot requires event and task capabilities",
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
                "emt.hybrid_state",
                1,
                0,
                :public,
                portable_state_inventory_record(hybrid_state),
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

"""Reconstruct a hybrid integrator in isolation through a caller-owned fresh callback factory."""
function restore_portable_emt_hybrid_snapshot(
    prepared::PreparedEMTStudy,
    snapshot::PortableEMTSnapshot;
    hybrid_factory,
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
    source_descriptor::Union{Nothing,PortableSnapshotDescriptor}=nothing,
)
    candidate = EMTStudyWorkspace(prepared)
    _portable_emt_validate_metadata(
        snapshot.metadata,
        candidate;
        project_signature_sha256,
        model_signature_sha256,
        settings_signature_sha256,
        source_descriptor,
    )
    all(capability -> capability in snapshot.metadata.capabilities, ("emt.events", "emt.tasks")) ||
        _portable_emt_fail(
            :capability_mismatch,
            "portable hybrid restore requires event and task capabilities",
        )
    public_section = _portable_emt_section(snapshot, "emt.public_state", :public)
    backend_section = _portable_emt_section(
        snapshot,
        "backend.reconstruction_state",
        :private_reconstructible,
    )
    hybrid_section = _portable_emt_section(snapshot, "emt.hybrid_state", :public)
    public_state = portable_state_inventory(public_section.value)
    backend_state = portable_state_inventory(backend_section.value)
    hybrid_state = portable_state_inventory(hybrid_section.value)
    integrator = Base.invokelatest(hybrid_factory, candidate)
    integrator isa EMTHybridStepIntegrator || _portable_emt_fail(
        :hybrid_factory_mismatch,
        "portable hybrid factory did not return an EMT hybrid integrator",
    )
    integrator.workspace === candidate || _portable_emt_fail(
        :hybrid_factory_mismatch,
        "portable hybrid factory returned an integrator for another workspace",
    )
    !integrator.initialized && candidate.ready &&
        !emt_step_transaction_status(integrator.transaction).active ||
        _portable_emt_fail(
            :hybrid_factory_mismatch,
            "portable hybrid factory must return a fresh inactive integrator",
        )
    restore_portable_emt_state_inventory!(candidate, public_state)
    candidate.execution_mode === :hybrid || _portable_emt_fail(
        :execution_mode_mismatch,
        "portable hybrid public state has the wrong execution owner",
    )
    _restore_portable_emt_backend_state!(candidate.runtime.timestep_state, backend_state)
    _restore_portable_emt_hybrid_state_inventory!(
        integrator,
        hybrid_state;
        task_program_signature_sha256,
        event_program_signature_sha256,
        isolated_reconstruction = true,
    )
    context = candidate.runtime.context
    _portable_emt_accepted_step(context) == snapshot.metadata.accepted_step ||
        _portable_emt_fail(
            :accepted_step_mismatch,
            "portable hybrid restored accepted step disagrees with its envelope",
        )
    _portable_emt_represented_time(context) == snapshot.metadata.represented_time_s ||
        _portable_emt_fail(
            :represented_time_mismatch,
            "portable hybrid restored time disagrees with its envelope",
        )
    restored_public = _portable_emt_state_inventory(candidate)
    restored_public.signature_sha256 == public_state.signature_sha256 ||
        _portable_emt_fail(
            :public_state_reconstruction,
            "portable hybrid public state changed during isolated reconstruction",
        )
    restored_backend = _portable_emt_backend_snapshot(candidate.runtime.timestep_state)
    restored_backend.signature_sha256 == backend_state.signature_sha256 ||
        _portable_emt_fail(
            :backend_state_reconstruction,
            "portable hybrid backend state changed during isolated reconstruction",
        )
    restored_hybrid = portable_emt_hybrid_state_inventory(
        integrator;
        task_program_signature_sha256,
        event_program_signature_sha256,
    )
    restored_hybrid.signature_sha256 == hybrid_state.signature_sha256 ||
        _portable_emt_fail(
            :hybrid_state_reconstruction,
            "portable hybrid coordinator state changed during isolated reconstruction",
        )
    return PortableEMTHybridRestoreResult(
        integrator,
        source_descriptor === nothing ? portable_snapshot_descriptor(snapshot) : source_descriptor,
        restored_public.signature_sha256,
        restored_backend.signature_sha256,
        restored_hybrid.signature_sha256,
        true,
    )
end

function write_portable_emt_hybrid_snapshot(
    path::AbstractString,
    integrator::EMTHybridStepIntegrator;
    kwargs...,
)
    snapshot = capture_portable_emt_snapshot(integrator; kwargs...)
    return write_portable_emt_snapshot(path, snapshot)
end

function read_portable_emt_hybrid_snapshot(
    path::AbstractString,
    prepared::PreparedEMTStudy;
    hybrid_factory,
    project_signature_sha256::AbstractString,
    model_signature_sha256::AbstractString,
    settings_signature_sha256::AbstractString,
    task_program_signature_sha256::AbstractString,
    event_program_signature_sha256::AbstractString,
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
    return restore_portable_emt_hybrid_snapshot(
        prepared,
        snapshot;
        hybrid_factory,
        project_signature_sha256,
        model_signature_sha256,
        settings_signature_sha256,
        task_program_signature_sha256,
        event_program_signature_sha256,
        source_descriptor = descriptor,
    )
end

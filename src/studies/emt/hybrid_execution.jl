const EMTHybridEventSurface = HybridEventSurface
const EMTHybridEventPolicy = HybridEventPolicy
const EMTHybridEventOccurrence = HybridEventOccurrence
const EMTSampledTaskOccurrence = SampledTaskOccurrence
const EMTExactSampledTask = ExactSampledTask
const EMTExactSampledControlTask = ExactSampledControlTask
const EMTExactPWMTask = ExactPWMTask
const EMTExactSampledTaskScheduler = ExactSampledTaskScheduler
const next_emt_sampled_task_time = next_sampled_task_time
const run_due_emt_sampled_tasks! = run_due_sampled_tasks!

_emt_task_scheduler_empty(scheduler::ExactSampledTaskScheduler) =
    isempty(scheduler.tasks)
_emt_task_scheduler_empty(scheduler::GeneralEMTTaskScheduler) =
    isempty(scheduler.tasks)
_emt_task_scheduler_empty(adapter::ExactSampledTaskCompatibilityAdapter) =
    _emt_task_scheduler_empty(adapter.scheduler)

function _emt_next_task_time(
    scheduler::ExactSampledTaskScheduler,
    current_time_s::Float64,
    endpoint_time_s::Float64,
    tolerance_s::Float64,
)
    return next_sampled_task_time(
        scheduler,
        current_time_s,
        endpoint_time_s;
        tolerance_s,
    )
end

function _emt_next_task_time(
    adapter::ExactSampledTaskCompatibilityAdapter,
    current_time_s::Float64,
    endpoint_time_s::Float64,
    tolerance_s::Float64,
)
    return _emt_next_task_time(
        adapter.scheduler,
        current_time_s,
        endpoint_time_s,
        tolerance_s,
    )
end

function _emt_next_task_time(
    scheduler::GeneralEMTTaskScheduler,
    current_time_s::Float64,
    endpoint_time_s::Float64,
    tolerance_s::Float64,
)
    instant = next_general_task_instant(scheduler)
    instant === nothing && return Inf
    time_s = Float64(instant)
    time_s < current_time_s - tolerance_s && throw(EMTTaskPlatform.EMTTaskPlatformFailure(
        :missed_task_activation,
        "general EMT task boundary precedes the active EMT interval";
        instant,
    ))
    return time_s <= endpoint_time_s + tolerance_s ? time_s : Inf
end

_emt_task_scheduler_checkpoint(scheduler::ExactSampledTaskScheduler) = (
    task_states = sampled_task_checkpoint.(scheduler.tasks),
    occurrence_count = length(scheduler.occurrences),
    invalidated_power = scheduler.last_run_power_history_invalidating,
)
_emt_task_scheduler_checkpoint(scheduler::GeneralEMTTaskScheduler) =
    general_task_scheduler_checkpoint(scheduler)
_emt_task_scheduler_checkpoint(adapter::ExactSampledTaskCompatibilityAdapter) =
    _emt_task_scheduler_checkpoint(adapter.scheduler)

function _restore_emt_task_scheduler!(scheduler::ExactSampledTaskScheduler, checkpoint)
    length(checkpoint.task_states) == length(scheduler.tasks) || throw(ArgumentError(
        "hybrid sampled-task state count changed during a boundary action",
    ))
    for (task, state) in zip(scheduler.tasks, checkpoint.task_states)
        restore_sampled_task_checkpoint!(task, state)
    end
    resize!(scheduler.occurrences, checkpoint.occurrence_count)
    scheduler.last_run_power_history_invalidating = checkpoint.invalidated_power
    return scheduler
end
_restore_emt_task_scheduler!(scheduler::GeneralEMTTaskScheduler, checkpoint) =
    restore_general_task_scheduler_checkpoint!(scheduler, checkpoint)
_restore_emt_task_scheduler!(adapter::ExactSampledTaskCompatibilityAdapter, checkpoint) =
    _restore_emt_task_scheduler!(adapter.scheduler, checkpoint)

function _run_due_emt_tasks!(
    scheduler::ExactSampledTaskScheduler,
    owner,
    time_s::Float64,
    tolerance_s::Float64,
)
    return run_due_sampled_tasks!(scheduler, owner, time_s; tolerance_s)
end

function _run_due_emt_tasks!(
    adapter::ExactSampledTaskCompatibilityAdapter,
    owner,
    time_s::Float64,
    tolerance_s::Float64,
)
    return _run_due_emt_tasks!(
        adapter.scheduler,
        owner,
        time_s,
        tolerance_s,
    )
end

function _run_due_emt_tasks!(
    scheduler::GeneralEMTTaskScheduler,
    owner,
    time_s::Float64,
    tolerance_s::Float64,
)
    instant = next_general_task_instant(scheduler)
    instant === nothing && return 0
    abs(Float64(instant) - time_s) <= tolerance_s || return 0
    return run_due_general_tasks!(scheduler, owner, instant; checkpoint_owner = false)
end

_emt_task_scheduler_invalidated_power(scheduler::ExactSampledTaskScheduler) =
    sampled_task_scheduler_last_run_invalidated_power(scheduler)
_emt_task_scheduler_invalidated_power(scheduler::GeneralEMTTaskScheduler) =
    general_task_scheduler_last_run_invalidated_power(scheduler)
_emt_task_scheduler_invalidated_power(adapter::ExactSampledTaskCompatibilityAdapter) =
    _emt_task_scheduler_invalidated_power(adapter.scheduler)

_emt_task_scheduler_occurrence_count(scheduler::ExactSampledTaskScheduler) =
    length(scheduler.occurrences)
_emt_task_scheduler_occurrence_count(scheduler::GeneralEMTTaskScheduler) =
    general_task_scheduler_occurrence_count(scheduler)
_emt_task_scheduler_occurrence_count(adapter::ExactSampledTaskCompatibilityAdapter) =
    _emt_task_scheduler_occurrence_count(adapter.scheduler)

struct EMTHybridCallbackOwner{R}
    runtime::R
end

"""A scheduler callback that routes an exact gate edge to one named power-semiconductor owner in the active EMT runtime."""
struct PowerSemiconductorGateCommand
    element_name::Symbol
end

"""A scheduler callback that routes an exact pole selection through one named complementary bridge-leg interlock."""
struct PowerSemiconductorBridgePoleCommand
    element_name::Symbol
end

"""A scheduler callback that routes an exact ordered gate-state request to one bridge topology aggregate."""
struct PowerSemiconductorTopologyGateCommand
    element_name::Symbol
end

function (command::PowerSemiconductorGateCommand)(
    owner::EMTHybridCallbackOwner,
    commanded_on::Bool,
    time_s::Real,
    _edge_index::Int,
)
    context = owner.runtime.context
    element_index = findfirst(==(command.element_name), context.element_names)
    element_index === nothing && throw(ArgumentError(
        "power-semiconductor gate target $(command.element_name) is absent from the EMT runtime",
    ))
    device = context.system.elements[element_index]
    device isa PowerSemiconductorSwitch || throw(ArgumentError(
        "gate target $(command.element_name) is not a power-semiconductor device",
    ))
    request_power_semiconductor_gate!(device, commanded_on, time_s)
    return nothing
end

function (command::PowerSemiconductorTopologyGateCommand)(
    owner::EMTHybridCallbackOwner,
    requested_state::AbstractVector{Bool},
    time_s::Real,
    _activation_index::Int,
)
    context = owner.runtime.context
    element_index = findfirst(==(command.element_name), context.element_names)
    element_index === nothing && throw(ArgumentError(
        "power-semiconductor topology target $(command.element_name) is absent from the EMT runtime",
    ))
    bridge = context.system.elements[element_index]
    bridge isa PowerSemiconductorBridgeTopology || throw(ArgumentError(
        "bridge target $(command.element_name) is not a power-semiconductor topology aggregate",
    ))
    request_power_semiconductor_topology_gates!(bridge, requested_state, time_s)
    return nothing
end

function (command::PowerSemiconductorBridgePoleCommand)(
    owner::EMTHybridCallbackOwner,
    upper_on::Bool,
    time_s::Real,
    _edge_index::Int,
)
    context = owner.runtime.context
    element_index = findfirst(==(command.element_name), context.element_names)
    element_index === nothing && throw(ArgumentError(
        "power-semiconductor bridge target $(command.element_name) is absent from the EMT runtime",
    ))
    bridge = context.system.elements[element_index]
    bridge isa PowerSemiconductorBridgeLeg || throw(ArgumentError(
        "bridge target $(command.element_name) is not a complementary power-semiconductor bridge leg",
    ))
    request_power_semiconductor_bridge_pole!(bridge, upper_on, time_s)
    return nothing
end

mutable struct EMTHybridStepIntegrator{W,O,T,S}
    workspace::W
    callback_owner::O
    transaction::T
    surfaces::Vector{AbstractHybridEventSurface}
    scheduler::S
    surface_fired::BitVector
    policy::HybridEventPolicy
    nominal_dt_s::Float64
    active_global_endpoint_s::Float64
    occurrences::Vector{HybridEventOccurrence}
    initialized::Bool
    completed::Bool
    evaluation_recorded::Bool
    accepted_interval_count::Int
    provisional_interval_count::Int
    localized_root_count::Int
    topology_invalidation_count::Int
    completed_global_step_count::Int
    last_global_step_event_count::Int
    last_failure::Union{Nothing,String}
end

function _emt_hybrid_element(owner::EMTHybridCallbackOwner, index::Int)
    return owner.runtime.context.system.elements[index]
end

function _emt_hybrid_current_extinction_value(element)
    if element isa CurrentZeroSwitch
        return element.current_initialized ? element.previous_current : nothing
    elseif element isa TimeSwitch && element.current_extinction !== nothing
        state = element.current_extinction
        return state.current_initialized ? state.previous_current : nothing
    end
    return nothing
end

function _emt_hybrid_current_extinction_threshold(element)
    if element isa CurrentZeroSwitch
        return element.critical_current_a
    elseif element isa TimeSwitch && element.current_extinction !== nothing
        return element.current_extinction.critical_current_a
    end
    return 0.0
end

function _emt_hybrid_current_extinction_reason(element)
    return _emt_hybrid_current_extinction_threshold(element) > 0.0 ?
        :critical_current : :current_reversal
end

function _emt_hybrid_current_extinction_surface_value(element)
    value = _emt_hybrid_current_extinction_value(element)
    value === nothing && return nothing
    threshold = _emt_hybrid_current_extinction_threshold(element)
    return threshold > 0.0 ? abs(value) - threshold : value
end

function _emt_hybrid_current_extinction_direction(element)
    return _emt_hybrid_current_extinction_threshold(element) > 0.0 ?
        HYBRID_EVENT_FALLING : HYBRID_EVENT_ANY
end

function _emt_hybrid_positive_surface_candidate_time(owner, value)
    value === nothing && return nothing
    Float64(value) > 0.0 || return nothing
    return owner.runtime.context.t_s
end

function _emt_hybrid_negative_surface_candidate_time(owner, value)
    value === nothing && return nothing
    Float64(value) < 0.0 || return nothing
    return owner.runtime.context.t_s
end

function _emt_hybrid_time_surface(
    name::Symbol,
    element_index::Int,
    event_time_s::Float64,
    priority::Int,
)
    return HybridEventSurface(
        name,
        _workspace -> nothing,
        (_workspace, _time_s) -> nothing;
        direction = HYBRID_EVENT_RISING,
        priority = priority,
        topology_invalidating = true,
        repeatable = false,
        candidate_time = workspace -> begin
            current_time = workspace.runtime.context.t_s
            current_time <= event_time_s ? event_time_s : nothing
        end,
    )
end

function _append_emt_power_semiconductor_conduction_surfaces!(
    surfaces::Vector{AbstractHybridEventSurface},
    owner_name::Symbol,
    configured_device::PowerSemiconductorSwitch,
    runtime_device,
)
    push!(
        surfaces,
        HybridEventSurface(
            Symbol(owner_name, :_forward_turn_on),
            owner -> power_semiconductor_forward_turn_on_residual(runtime_device(owner)),
            (owner, time_s) -> apply_power_semiconductor_forward_turn_on!(
                runtime_device(owner),
                time_s,
            );
            direction = HYBRID_EVENT_RISING,
            priority = -8,
            topology_invalidating = true,
            repeatable = true,
            candidate_time = owner -> _emt_hybrid_positive_surface_candidate_time(
                owner,
                power_semiconductor_forward_turn_on_residual(runtime_device(owner)),
            ),
        ),
        HybridEventSurface(
            Symbol(owner_name, :_forward_extinction),
            owner -> power_semiconductor_forward_extinction_residual(runtime_device(owner)),
            (owner, time_s) -> apply_power_semiconductor_forward_extinction!(
                runtime_device(owner),
                time_s,
            );
            direction = HYBRID_EVENT_FALLING,
            priority = -7,
            topology_invalidating = true,
            repeatable = true,
            candidate_time = owner -> _emt_hybrid_negative_surface_candidate_time(
                owner,
                power_semiconductor_forward_extinction_residual(runtime_device(owner)),
            ),
        ),
    )
    configured_device.antiparallel_diode === nothing && return surfaces
    push!(
        surfaces,
        HybridEventSurface(
            Symbol(owner_name, :_reverse_diode_turn_on),
            owner -> power_semiconductor_reverse_turn_on_residual(runtime_device(owner)),
            (owner, time_s) -> apply_power_semiconductor_reverse_turn_on!(
                runtime_device(owner),
                time_s,
            );
            direction = HYBRID_EVENT_RISING,
            priority = -8,
            topology_invalidating = true,
            repeatable = true,
            candidate_time = owner -> _emt_hybrid_positive_surface_candidate_time(
                owner,
                power_semiconductor_reverse_turn_on_residual(runtime_device(owner)),
            ),
        ),
        HybridEventSurface(
            Symbol(owner_name, :_reverse_diode_extinction),
            owner -> power_semiconductor_reverse_extinction_residual(runtime_device(owner)),
            (owner, time_s) -> apply_power_semiconductor_reverse_extinction!(
                runtime_device(owner),
                time_s,
            );
            direction = HYBRID_EVENT_FALLING,
            priority = -7,
            topology_invalidating = true,
            repeatable = true,
            candidate_time = owner -> _emt_hybrid_negative_surface_candidate_time(
                owner,
                power_semiconductor_reverse_extinction_residual(runtime_device(owner)),
            ),
        ),
    )
    return surfaces
end

function _emt_hybrid_device_surfaces(workspace::EMTStudyWorkspace)
    context = workspace.runtime.context
    surfaces = AbstractHybridEventSurface[]
    for index in eachindex(context.system.elements)
        element = context.system.elements[index]
        owner_name = context.element_names[index]
        if element isa PowerSemiconductorSwitch
            power_semiconductor_event_localization!(element)
            if element.gate_driver !== nothing
                push!(
                    surfaces,
                    HybridEventSurface(
                        Symbol(owner_name, :_gate_transition),
                        _owner -> nothing,
                        (owner, time_s) -> apply_power_semiconductor_gate_transition!(
                            _emt_hybrid_element(owner, index),
                            time_s,
                        );
                        direction = HYBRID_EVENT_RISING,
                        priority = -15,
                        topology_invalidating = true,
                        repeatable = true,
                        candidate_time = owner -> power_semiconductor_gate_transition_time(
                            _emt_hybrid_element(owner, index),
                        ),
                    ),
                )
            end
            _append_emt_power_semiconductor_conduction_surfaces!(
                surfaces,
                owner_name,
                element,
                owner -> _emt_hybrid_element(owner, index),
            )
        elseif element isa PowerSemiconductorBridgeLeg
            power_semiconductor_event_localization!(element)
            push!(
                surfaces,
                HybridEventSurface(
                    Symbol(owner_name, :_gate_transition),
                    _owner -> nothing,
                    (owner, time_s) -> apply_power_semiconductor_bridge_gate_transitions!(
                        _emt_hybrid_element(owner, index),
                        time_s,
                    );
                    direction = HYBRID_EVENT_RISING,
                    priority = -15,
                    topology_invalidating = true,
                    repeatable = true,
                    candidate_time = owner -> power_semiconductor_bridge_gate_transition_time(
                        _emt_hybrid_element(owner, index),
                    ),
                ),
            )
            for position in (:upper, :lower)
                switch = power_semiconductor_bridge_switch(element, position)
                _append_emt_power_semiconductor_conduction_surfaces!(
                    surfaces,
                    Symbol(owner_name, :_, position),
                    switch,
                    owner -> power_semiconductor_bridge_switch(
                        _emt_hybrid_element(owner, index),
                        position,
                    ),
                )
            end
        elseif element isa PowerSemiconductorBridgeTopology
            power_semiconductor_event_localization!(element)
            push!(
                surfaces,
                HybridEventSurface(
                    Symbol(owner_name, :_gate_transition),
                    _owner -> nothing,
                    (owner, time_s) -> apply_power_semiconductor_bridge_gate_transitions!(
                        _emt_hybrid_element(owner, index),
                        time_s,
                    );
                    direction = HYBRID_EVENT_RISING,
                    priority = -15,
                    topology_invalidating = true,
                    repeatable = true,
                    candidate_time = owner -> power_semiconductor_bridge_gate_transition_time(
                        _emt_hybrid_element(owner, index),
                    ),
                ),
            )
            for (position_index, switch) in
                enumerate(power_semiconductor_bridge_topology_valves(element))
                position_name = element.topology.valve_positions[position_index].name
                _append_emt_power_semiconductor_conduction_surfaces!(
                    surfaces,
                    Symbol(owner_name, :_, position_name),
                    switch,
                    owner -> power_semiconductor_bridge_topology_valves(
                        _emt_hybrid_element(owner, index),
                    )[position_index],
                )
            end
        end
        if element isa TimeSwitch
            if isfinite(element.close_time_s)
                push!(
                    surfaces,
                    _emt_hybrid_time_surface(
                        Symbol(owner_name, :_close),
                        index,
                        element.close_time_s,
                        -20,
                    ),
                )
            end
            if element.current_extinction === nothing && isfinite(element.open_time_s)
                push!(
                    surfaces,
                    _emt_hybrid_time_surface(
                        Symbol(owner_name, :_open),
                        index,
                        element.open_time_s,
                        -10,
                    ),
                )
            end
        end
        if element isa CurrentZeroSwitch ||
           element isa TimeSwitch && element.current_extinction !== nothing
            direction = _emt_hybrid_current_extinction_direction(element)
            push!(
                surfaces,
                HybridEventSurface(
                    Symbol(owner_name, :_current_extinction),
                    owner -> _emt_hybrid_current_extinction_surface_value(
                        _emt_hybrid_element(owner, index),
                    ),
                    (owner, time_s) -> begin
                        target = _emt_hybrid_element(owner, index)
                        apply_current_zero_transition!(
                            target,
                            _emt_hybrid_current_extinction_reason(target),
                            time_s,
                        )
                    end;
                    direction = direction,
                    priority = -5,
                    topology_invalidating = true,
                    repeatable = false,
                ),
            )
        elseif element isa TACSControlledSwitch && element.delayed_arc !== nothing
            element.delayed_arc.transition_deferred = true
            push!(
                surfaces,
                HybridEventSurface(
                    Symbol(owner_name, :_delayed_arc_open),
                    _owner -> nothing,
                    (owner, time_s) -> begin
                        target = _emt_hybrid_element(owner, index)
                        apply_controlled_switch_delayed_arc_transition!(target, time_s)
                        target.delayed_arc.transition_deferred = true
                    end;
                    direction = HYBRID_EVENT_FALLING,
                    priority = -5,
                    topology_invalidating = true,
                    repeatable = true,
                    candidate_time = owner -> begin
                        target = _emt_hybrid_element(owner, index)
                        arc = target.delayed_arc
                        arc !== nothing && target.closed && arc.opening_requested &&
                        isfinite(arc.scheduled_open_time_s) ?
                            arc.scheduled_open_time_s : nothing
                    end,
                ),
            )
        end
    end
    return surfaces
end

function _check_emt_hybrid_surfaces!(surfaces)
    all(surface -> surface isa HybridEventSurface, surfaces) || throw(ArgumentError(
        "EMT hybrid execution currently requires typed HybridEventSurface owners",
    ))
    names = Symbol[surface.name for surface in surfaces]
    length(unique(names)) == length(names) || throw(ArgumentError(
        "EMT hybrid event surface names must be unique",
    ))
    return surfaces
end

"""Configure hybrid EMT execution; callbacks receive a runtime-only owner, surface reads must be pure, and successful transitions/tasks may mutate only that owner before returning."""
function configure_emt_hybrid_execution(
    workspace::EMTStudyWorkspace;
    event_surfaces::AbstractVector{<:AbstractHybridEventSurface} =
        AbstractHybridEventSurface[],
    sampled_tasks::AbstractVector{<:AbstractExactSampledTask} =
        AbstractExactSampledTask[],
    scheduler::Union{
        Nothing,
        ExactSampledTaskScheduler,
        GeneralEMTTaskScheduler,
        ExactSampledTaskCompatibilityAdapter,
    } = nothing,
    scheduler_tick_s::Union{Nothing,Real} = nothing,
    scheduler_origin_s::Real = 0.0,
    policy::HybridEventPolicy = HybridEventPolicy(),
    include_device_events::Bool = true,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be ready before configuring hybrid execution",
    ))
    workspace.execution_mode === :unselected || throw(ArgumentError(
        "EMT study workspace is already owned by $(workspace.execution_mode) execution",
    ))
    runtime = _check_prepared_runtime_aliases(workspace.runtime)
    nominal_dt = runtime.context.dt_s
    nominal_dt > 0.0 && isfinite(nominal_dt) || throw(ArgumentError(
        "EMT hybrid execution requires a finite positive nominal timestep",
    ))
    surfaces = AbstractHybridEventSurface[]
    include_device_events && append!(surfaces, _emt_hybrid_device_surfaces(workspace))
    append!(surfaces, event_surfaces)
    _check_emt_hybrid_surfaces!(surfaces)
    resolved_scheduler = if scheduler === nothing
        tick = scheduler_tick_s === nothing ? nominal_dt / 1_000_000.0 :
            Float64(scheduler_tick_s)
        ExactSampledTaskScheduler(
            tick;
            origin_s = scheduler_origin_s,
            tasks = sampled_tasks,
        )
    else
        isempty(sampled_tasks) || throw(ArgumentError(
            "sampled_tasks cannot be supplied together with an existing scheduler",
        ))
        scheduler_tick_s === nothing || throw(ArgumentError(
            "scheduler_tick_s cannot be supplied together with an existing scheduler",
        ))
        scheduler
    end
    integrator = EMTHybridStepIntegrator(
        workspace,
        EMTHybridCallbackOwner(runtime),
        EMTStepTransaction(workspace),
        surfaces,
        resolved_scheduler,
        falses(length(surfaces)),
        policy,
        nominal_dt,
        NaN,
        HybridEventOccurrence[],
        false,
        false,
        false,
        0,
        0,
        0,
        0,
        0,
        0,
        nothing,
    )
    workspace.execution_mode = :hybrid
    return integrator
end

function _emt_hybrid_initialize!(integrator::EMTHybridStepIntegrator)
    integrator.initialized && return integrator
    runtime = integrator.workspace.runtime
    context = runtime.context
    if runtime.steady_state_initial_sample !== nothing
        _apply_due_series_rlc_alterations!(context)
        _apply_steady_state_initial_sample!(
            context,
            runtime.steady_state_initial_sample,
            runtime.steady_state_initial_output_values,
        )
        context.step_index = 1
        context.t_s = min(context.step_index, context.step_count) * integrator.nominal_dt_s
    elseif context.step_index <= context.step_count
        _advance_prepared_emt_step!(runtime; collect_step_diagnostics = false)
    end
    integrator.completed = context.step_index > context.step_count
    accepted_time = context.step_count == 0 ? 0.0 :
        max(0.0, context.t_s - integrator.nominal_dt_s)
    next_target = context.t_s
    context.t_s = accepted_time
    try
        _emt_hybrid_apply_boundary_actions!(
            integrator,
            Int[],
            accepted_time,
            Dict{Int,Any}(),
        )
    finally
        context.t_s = next_target
    end
    integrator.initialized = true
    return integrator
end

function _emt_hybrid_advance_interval!(
    integrator::EMTHybridStepIntegrator,
    left_time_s::Float64,
    target_time_s::Float64,
)
    delta = target_time_s - left_time_s
    delta > 0.0 || throw(ArgumentError(
        "hybrid interval must make positive forward progress",
    ))
    runtime = integrator.workspace.runtime
    context = runtime.context
    context.dt_s = delta
    context.t_s = target_time_s
    _apply_due_series_rlc_alterations!(context)
    final_global_interval = abs(
        target_time_s - integrator.active_global_endpoint_s,
    ) <= integrator.policy.root_time_tolerance_s
    over16_kwargs = merge(
        _over16_step_kwargs(runtime.step_configs, context),
        (
            hybrid_substep = true,
            hybrid_report_step_s = final_global_interval ?
                integrator.nominal_dt_s : delta,
        ),
    )
    if !final_global_interval
        over16_kwargs = Base.structdiff(
            over16_kwargs,
            (output_report_config = nothing, post_extrema_config = nothing),
        )
    end
    if haskey(over16_kwargs, :source_config)
        source_config = over16_kwargs.source_config
        source_runtime = context.source_function_runtime
        staged_time = integrator.active_global_endpoint_s +
            integrator.nominal_dt_s
        source_time = source_runtime === nothing ? staged_time :
            max(staged_time, source_runtime.last_accepted_time_s)
        over16_kwargs = merge(
            over16_kwargs,
            (
                source_config = merge(
                    source_config,
                    (t = source_time, constraint_t = target_time_s),
                ),
            ),
        )
    end
    return _step_with_over16_config!(
        context,
        runtime.timestep_state,
        over16_kwargs,
        false,
    )
end

function _emt_hybrid_restore_if_active!(integrator::EMTHybridStepIntegrator)
    timestep_transaction_active(integrator.transaction.transaction) &&
        restore_emt_step_transaction!(integrator.transaction)
    return integrator
end

function _emt_hybrid_provisional_values!(
    integrator::EMTHybridStepIntegrator,
    left_time_s::Float64,
    target_time_s::Float64,
)
    workspace_ready = integrator.workspace.ready
    evaluation_count = integrator.workspace.evaluation_count
    begin_emt_step_transaction!(integrator.transaction)
    try
        _emt_hybrid_advance_interval!(integrator, left_time_s, target_time_s)
        values = Any[
            hybrid_event_value(surface, integrator.callback_owner)
            for surface in integrator.surfaces
        ]
        restore_emt_step_transaction!(integrator.transaction)
        integrator.provisional_interval_count += 1
        return values
    catch
        _emt_hybrid_restore_if_active!(integrator)
        rethrow()
    finally
        integrator.workspace.ready = workspace_ready
        integrator.workspace.evaluation_count = evaluation_count
    end
end

function _emt_hybrid_provisional_surface_value!(
    integrator::EMTHybridStepIntegrator,
    surface_index::Int,
    left_time_s::Float64,
    target_time_s::Float64,
)
    workspace_ready = integrator.workspace.ready
    evaluation_count = integrator.workspace.evaluation_count
    begin_emt_step_transaction!(integrator.transaction)
    try
        _emt_hybrid_advance_interval!(integrator, left_time_s, target_time_s)
        value = hybrid_event_value(
            integrator.surfaces[surface_index],
            integrator.callback_owner,
        )
        value === nothing && throw(ArgumentError(
            "hybrid event surface became unavailable inside its root bracket",
        ))
        restore_emt_step_transaction!(integrator.transaction)
        integrator.provisional_interval_count += 1
        return value
    catch
        _emt_hybrid_restore_if_active!(integrator)
        rethrow()
    finally
        integrator.workspace.ready = workspace_ready
        integrator.workspace.evaluation_count = evaluation_count
    end
end

function _emt_hybrid_source_snapshot(integrator::EMTHybridStepIntegrator)
    runtime = integrator.workspace.runtime
    source_runtime = runtime.context.source_function_runtime
    network_runtime = source_runtime === nothing ? nothing : (
        state = deepcopy(source_runtime.state),
        slot_values = Float64[slot[] for slot in source_runtime.slot_values],
        executed_step_count = source_runtime.executed_step_count,
        card_read_count = source_runtime.card_read_count,
        signal_synchronization_count = source_runtime.signal_synchronization_count,
        external_signal_count = source_runtime.external_signal_count,
        tacs_override_count = source_runtime.tacs_override_count,
        analytic_execution_count = source_runtime.analytic_execution_count,
        stage_samples = deepcopy(source_runtime.stage_samples),
        last_accepted_time_s = source_runtime.last_accepted_time_s,
        next_input_row_index = source_runtime.next_input_row_index,
    )
    return (
        accepted_source_card = deepcopy(runtime.timestep_state.source.source_card),
        network_runtime = network_runtime,
    )
end

function _restore_emt_hybrid_source_snapshot!(
    integrator::EMTHybridStepIntegrator,
    snapshot,
)
    runtime = integrator.workspace.runtime
    restore_timestep_state!(
        runtime.timestep_state.source.source_card,
        snapshot.accepted_source_card,
    )
    source_runtime = runtime.context.source_function_runtime
    if source_runtime !== nothing
        snapshot.network_runtime === nothing && throw(ArgumentError(
            "hybrid source snapshot is missing the network runtime",
        ))
        network = snapshot.network_runtime
        restore_timestep_state!(source_runtime.state, network.state)
        length(source_runtime.slot_values) == length(network.slot_values) ||
            throw(ArgumentError("hybrid source slot count changed inside a substep"))
        for (slot, value) in zip(source_runtime.slot_values, network.slot_values)
            slot[] = value
        end
        source_runtime.executed_step_count = network.executed_step_count
        source_runtime.card_read_count = network.card_read_count
        source_runtime.signal_synchronization_count =
            network.signal_synchronization_count
        source_runtime.external_signal_count = network.external_signal_count
        source_runtime.tacs_override_count = network.tacs_override_count
        source_runtime.analytic_execution_count = network.analytic_execution_count
        restore_timestep_state!(source_runtime.stage_samples, network.stage_samples)
        source_runtime.last_accepted_time_s = network.last_accepted_time_s
        source_runtime.next_input_row_index = network.next_input_row_index
    end
    _check_prepared_runtime_aliases(runtime)
    return integrator
end

function _emt_hybrid_accept_interval!(
    integrator::EMTHybridStepIntegrator,
    left_time_s::Float64,
    target_time_s::Float64,
)
    final_global_interval = abs(
        target_time_s - integrator.active_global_endpoint_s,
    ) <= integrator.policy.root_time_tolerance_s
    source_snapshot =
        final_global_interval ? nothing : _emt_hybrid_source_snapshot(integrator)
    begin_emt_step_transaction!(integrator.transaction)
    try
        _emt_hybrid_advance_interval!(integrator, left_time_s, target_time_s)
        commit_emt_step_transaction!(integrator.transaction)
        final_global_interval ||
            _restore_emt_hybrid_source_snapshot!(integrator, source_snapshot)
        integrator.accepted_interval_count += 1
        return integrator
    catch
        _emt_hybrid_restore_if_active!(integrator)
        rethrow()
    end
end

function _emt_hybrid_candidate_time(
    integrator::EMTHybridStepIntegrator,
    index::Int,
    left_time_s::Float64,
    endpoint_time_s::Float64,
)
    integrator.surface_fired[index] && return nothing
    candidate = hybrid_event_candidate_time(
        integrator.surfaces[index],
        integrator.callback_owner,
    )
    candidate === nothing && return nothing
    tolerance = integrator.policy.simultaneity_tolerance_s
    candidate < left_time_s - tolerance && throw(ArgumentError(
        "hybrid event $(integrator.surfaces[index].name) supplied a past candidate time",
    ))
    candidate <= endpoint_time_s + tolerance || return nothing
    return clamp(candidate, left_time_s, endpoint_time_s)
end

function _emt_hybrid_sorted_indices(
    integrator::EMTHybridStepIntegrator,
    indices::Vector{Int},
)
    sort!(indices; by = index -> begin
        surface = integrator.surfaces[index]
        (surface.priority, String(surface.name), index)
    end)
    return indices
end

function _invalidate_emt_hybrid_power_history!(context::EMTStepContext)
    fill!(context.branch_power_history_valid, false)
    fill!(context.switch_power_history_valid, false)
    return context
end

function _emt_hybrid_apply_events!(
    integrator::EMTHybridStepIntegrator,
    indices::Vector{Int},
    time_s::Float64,
    roots::Dict{Int,Any},
)
    isempty(indices) && return 0
    ordered_indices = _emt_hybrid_sorted_indices(integrator, unique(indices))
    pending_indices = Int[
        index for index in ordered_indices if !integrator.surface_fired[index]
    ]
    projected_event_count =
        integrator.last_global_step_event_count + length(pending_indices)
    projected_event_count <= integrator.policy.maximum_events_per_step ||
        throw(ArgumentError(
            "hybrid execution would exceed maximum_events_per_step=$(integrator.policy.maximum_events_per_step) at t=$time_s s",
        ))
    applied = 0
    for index in pending_indices
        surface = integrator.surfaces[index]
        apply_hybrid_event!(surface, integrator.callback_owner, time_s)
        root = get(roots, index, nothing)
        value = root === nothing ? something(
            hybrid_event_value(surface, integrator.callback_owner),
            0.0,
        ) : root.value
        push!(
            integrator.occurrences,
            HybridEventOccurrence(
                surface.name,
                time_s,
                Float64(value),
                surface.priority,
                surface.topology_invalidating,
                root === nothing ? 0 : root.iteration_count,
                root === nothing ? 0.0 : root.bracket_width_s,
            ),
        )
        surface.topology_invalidating &&
            (integrator.topology_invalidation_count += 1)
        surface.repeatable || (integrator.surface_fired[index] = true)
        applied += 1
    end
    integrator.last_global_step_event_count += applied
    any(
        index -> integrator.surfaces[index].topology_invalidating,
        pending_indices,
    ) && _invalidate_emt_hybrid_power_history!(integrator.workspace.runtime.context)
    return applied
end

function _emt_hybrid_has_due_tasks(
    integrator::EMTHybridStepIntegrator,
    time_s::Float64,
)
    _emt_task_scheduler_empty(integrator.scheduler) && return false
    task_time = _emt_next_task_time(
        integrator.scheduler,
        time_s,
        time_s,
        integrator.policy.simultaneity_tolerance_s,
    )
    return isfinite(task_time) &&
        abs(task_time - time_s) <= integrator.policy.simultaneity_tolerance_s
end

function _emt_hybrid_apply_boundary_actions!(
    integrator::EMTHybridStepIntegrator,
    indices::Vector{Int},
    time_s::Float64,
    roots::Dict{Int,Any},
)
    tasks_due = _emt_hybrid_has_due_tasks(integrator, time_s)
    isempty(indices) && !tasks_due && return (event_count = 0, task_count = 0)
    fired_before = copy(integrator.surface_fired)
    occurrence_count_before = length(integrator.occurrences)
    scheduler_before = _emt_task_scheduler_checkpoint(integrator.scheduler)
    topology_count_before = integrator.topology_invalidation_count
    step_event_count_before = integrator.last_global_step_event_count
    workspace_ready_before = integrator.workspace.ready
    evaluation_count_before = integrator.workspace.evaluation_count
    begin_emt_step_transaction!(integrator.transaction)
    try
        event_count = _emt_hybrid_apply_events!(integrator, indices, time_s, roots)
        task_count = _run_due_emt_tasks!(
            integrator.scheduler,
            integrator.callback_owner,
            time_s,
            integrator.policy.simultaneity_tolerance_s,
        )
        task_count > 0 &&
            _emt_task_scheduler_invalidated_power(integrator.scheduler) &&
            _invalidate_emt_hybrid_power_history!(integrator.workspace.runtime.context)
        commit_emt_step_transaction!(integrator.transaction)
        return (event_count = event_count, task_count = task_count)
    catch
        _emt_hybrid_restore_if_active!(integrator)
        copyto!(integrator.surface_fired, fired_before)
        resize!(integrator.occurrences, occurrence_count_before)
        _restore_emt_task_scheduler!(integrator.scheduler, scheduler_before)
        integrator.topology_invalidation_count = topology_count_before
        integrator.last_global_step_event_count = step_event_count_before
        integrator.workspace.ready = workspace_ready_before
        integrator.workspace.evaluation_count = evaluation_count_before
        rethrow()
    end
end

function _emt_hybrid_candidate_indices_at!(
    destination::Vector{Int},
    integrator::EMTHybridStepIntegrator,
    time_s::Float64,
    endpoint_time_s::Float64,
)
    empty!(destination)
    tolerance = integrator.policy.simultaneity_tolerance_s
    for index in eachindex(integrator.surfaces)
        candidate = _emt_hybrid_candidate_time(
            integrator,
            index,
            time_s,
            endpoint_time_s,
        )
        candidate === nothing && continue
        abs(candidate - time_s) <= tolerance || continue
        hybrid_event_candidate_is_event(integrator.surfaces[index]) &&
            push!(destination, index)
    end
    return destination
end

function _emt_hybrid_next_boundary(
    integrator::EMTHybridStepIntegrator,
    left_time_s::Float64,
    endpoint_time_s::Float64,
)
    tolerance = integrator.policy.simultaneity_tolerance_s
    boundary = endpoint_time_s
    if !_emt_task_scheduler_empty(integrator.scheduler) &&
       endpoint_time_s > left_time_s + tolerance
        task_time = _emt_next_task_time(
            integrator.scheduler,
            left_time_s + tolerance,
            endpoint_time_s,
            tolerance,
        )
        isfinite(task_time) && task_time > left_time_s + tolerance &&
            (boundary = min(boundary, task_time))
    end
    for index in eachindex(integrator.surfaces)
        candidate = _emt_hybrid_candidate_time(
            integrator,
            index,
            left_time_s,
            endpoint_time_s,
        )
        candidate === nothing && continue
        candidate > left_time_s + tolerance && (boundary = min(boundary, candidate))
    end
    return boundary
end

function _emt_hybrid_localized_roots!(
    integrator::EMTHybridStepIntegrator,
    left_time_s::Float64,
    right_time_s::Float64,
    left_values,
    right_values,
)
    roots = Dict{Int,Any}()
    earliest = Inf
    for index in eachindex(integrator.surfaces)
        integrator.surface_fired[index] && continue
        surface = integrator.surfaces[index]
        bracket = hybrid_event_bracket(
            surface.direction,
            left_time_s,
            left_values[index],
            right_time_s,
            right_values[index],
            integrator.policy,
        )
        bracket === nothing && continue
        root = localize_hybrid_event!(
            time_s -> _emt_hybrid_provisional_surface_value!(
                integrator,
                index,
                left_time_s,
                time_s,
            ),
            surface.direction,
            bracket,
            integrator.policy,
        )
        root.converged || throw(ArgumentError(
            "hybrid event $(surface.name) failed to converge after $(root.iteration_count) iterations with bracket width $(root.bracket_width_s) s",
        ))
        roots[index] = root
        earliest = min(earliest, root.time_s)
        integrator.localized_root_count += 1
    end
    return roots, earliest
end

function _emt_hybrid_root_indices(
    integrator::EMTHybridStepIntegrator,
    roots::Dict{Int,Any},
    event_time_s::Float64,
)
    tolerance = integrator.policy.simultaneity_tolerance_s
    return Int[
        index for (index, root) in roots
        if abs(root.time_s - event_time_s) <= tolerance
    ]
end

function _emt_hybrid_finish_global_step!(
    integrator::EMTHybridStepIntegrator,
    endpoint_time_s::Float64,
)
    context = integrator.workspace.runtime.context
    context.dt_s = integrator.nominal_dt_s
    context.t_s = endpoint_time_s
    record_step!(context, context.system.v; update_power_energy = false)
    integrator.completed_global_step_count += 1
    integrator.completed = context.step_index > context.step_count
    return integrator
end

function _advance_emt_hybrid_step_impl!(integrator::EMTHybridStepIntegrator)
    integrator.workspace.ready = false
    _emt_hybrid_initialize!(integrator)
    integrator.completed && return integrator
    workspace = integrator.workspace
    context = workspace.runtime.context
    if isempty(integrator.surfaces) && _emt_task_scheduler_empty(integrator.scheduler)
        context.dt_s = integrator.nominal_dt_s
        _advance_prepared_emt_step!(
            workspace.runtime;
            collect_step_diagnostics = false,
        )
        integrator.accepted_interval_count += 1
        integrator.completed_global_step_count += 1
        integrator.completed = context.step_index > context.step_count
        return integrator
    end
    endpoint_time = context.t_s
    integrator.active_global_endpoint_s = endpoint_time
    left_time = max(0.0, endpoint_time - integrator.nominal_dt_s)
    context.t_s = left_time
    context.dt_s = integrator.nominal_dt_s
    integrator.last_global_step_event_count = 0
    candidate_indices = Int[]
    try
        while left_time < endpoint_time - integrator.policy.root_time_tolerance_s
            _emt_hybrid_candidate_indices_at!(
                candidate_indices,
                integrator,
                left_time,
                endpoint_time,
            )
            _emt_hybrid_apply_boundary_actions!(
                integrator,
                candidate_indices,
                left_time,
                Dict{Int,Any}(),
            )
            left_values = Any[
                hybrid_event_value(surface, integrator.callback_owner)
                for surface in integrator.surfaces
            ]
            boundary = _emt_hybrid_next_boundary(
                integrator,
                left_time,
                endpoint_time,
            )
            boundary_candidate_indices = Int[]
            for index in eachindex(integrator.surfaces)
                candidate = _emt_hybrid_candidate_time(
                    integrator,
                    index,
                    left_time,
                    endpoint_time,
                )
                candidate === nothing && continue
                abs(candidate - boundary) <=
                    integrator.policy.simultaneity_tolerance_s || continue
                hybrid_event_candidate_is_event(integrator.surfaces[index]) &&
                    push!(boundary_candidate_indices, index)
            end
            boundary - left_time >= integrator.policy.minimum_progress_s ||
                throw(ArgumentError(
                    "hybrid execution cannot make minimum progress at t=$left_time s toward boundary $boundary s",
                ))
            if all(integrator.surface_fired)
                _emt_hybrid_accept_interval!(
                    integrator,
                    left_time,
                    boundary,
                )
                _emt_hybrid_apply_boundary_actions!(
                    integrator,
                    Int[],
                    boundary,
                    Dict{Int,Any}(),
                )
                left_time = boundary
                context.t_s = left_time
                continue
            end
            right_values = _emt_hybrid_provisional_values!(
                integrator,
                left_time,
                boundary,
            )
            roots, earliest_root = _emt_hybrid_localized_roots!(
                integrator,
                left_time,
                boundary,
                left_values,
                right_values,
            )
            target = isfinite(earliest_root) ? min(boundary, earliest_root) : boundary
            target - left_time >= integrator.policy.minimum_progress_s ||
                throw(ArgumentError(
                    "hybrid event localization produced chatter or zero progress at t=$left_time s",
                ))
            _emt_hybrid_accept_interval!(
                integrator,
                left_time,
                target,
            )
            root_indices = isfinite(earliest_root) ?
                _emt_hybrid_root_indices(integrator, roots, target) : Int[]
            if abs(target - boundary) <=
               integrator.policy.simultaneity_tolerance_s
                append!(root_indices, boundary_candidate_indices)
            end
            _emt_hybrid_apply_boundary_actions!(
                integrator,
                root_indices,
                target,
                roots,
            )
            left_time = target
            context.t_s = left_time
        end
        return _emt_hybrid_finish_global_step!(integrator, endpoint_time)
    catch error
        integrator.last_failure = sprint(showerror, error)
        _emt_hybrid_restore_if_active!(integrator)
        rethrow()
    end
end

function advance_emt_hybrid_step!(integrator::EMTHybridStepIntegrator)
    integrator.last_failure === nothing || throw(ArgumentError(
        "EMT hybrid integrator is terminally failed and cannot be reused: $(integrator.last_failure)",
    ))
    try
        return _advance_emt_hybrid_step_impl!(integrator)
    catch error
        integrator.last_failure = sprint(showerror, error)
        _emt_hybrid_restore_if_active!(integrator)
        rethrow()
    end
end

function _emt_hybrid_trace(integrator::EMTHybridStepIntegrator)
    workspace = integrator.workspace
    trace = _deck_requested_electrical_trace(
        workspace.parsed,
        deck_trace(workspace.runtime.context),
    )
    initial_control_voltage =
        workspace.runtime.steady_state_initial_sample === nothing ? nothing :
        workspace.runtime.steady_state_initial_sample.node_voltage_values
    return _append_deck_control_system_outputs(
        trace,
        workspace.parsed,
        initial_control_voltage,
    )
end

function evaluate_emt_hybrid_study!(integrator::EMTHybridStepIntegrator)
    integrator.evaluation_recorded && throw(ArgumentError(
        "EMT hybrid integrator has already returned its completed evaluation",
    ))
    while !integrator.completed
        advance_emt_hybrid_step!(integrator)
    end
    integrator.evaluation_recorded = true
    integrator.workspace.evaluation_count += 1
    return _emt_hybrid_trace(integrator)
end

function evaluate_emt_hybrid_study!(workspace::EMTStudyWorkspace; kwargs...)
    integrator = configure_emt_hybrid_execution(workspace; kwargs...)
    return evaluate_emt_hybrid_study!(integrator)
end

function emt_hybrid_execution_status(integrator::EMTHybridStepIntegrator)
    transaction = emt_step_transaction_status(integrator.transaction)
    status_scheduler = integrator.scheduler isa ExactSampledTaskCompatibilityAdapter ?
        integrator.scheduler.scheduler : integrator.scheduler
    sampled_controls = [
        task for task in status_scheduler.tasks
        if task isa ExactSampledControlTask
    ]
    pwm_tasks = [
        task for task in status_scheduler.tasks
        if task isa ExactPWMTask
    ]
    return (
        initialized = integrator.initialized,
        completed = integrator.completed,
        evaluation_recorded = integrator.evaluation_recorded,
        accepted_interval_count = integrator.accepted_interval_count,
        provisional_interval_count = integrator.provisional_interval_count,
        localized_root_count = integrator.localized_root_count,
        event_count = length(integrator.occurrences),
        sampled_task_execution_count = _emt_task_scheduler_occurrence_count(
            integrator.scheduler,
        ),
        sampled_control_sample_count = sum(
            task -> task.sample_count,
            sampled_controls;
            init = 0,
        ),
        sampled_control_write_count = sum(
            task -> task.write_count,
            sampled_controls;
            init = 0,
        ),
        pwm_cycle_count = sum(task -> task.cycle_count, pwm_tasks; init = 0),
        pwm_edge_count = sum(task -> task.edge_count, pwm_tasks; init = 0),
        boundary_action_order = (:events, :sampled_tasks),
        topology_invalidation_count = integrator.topology_invalidation_count,
        completed_global_step_count = integrator.completed_global_step_count,
        last_global_step_event_count = integrator.last_global_step_event_count,
        transaction_capture_count = transaction.capture_count,
        transaction_restore_count = transaction.restore_count,
        transaction_commit_count = transaction.commit_count,
        last_failure = integrator.last_failure,
    )
end

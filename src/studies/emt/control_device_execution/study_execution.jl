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

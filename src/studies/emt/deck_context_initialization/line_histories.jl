function _distributed_transposed_line_modal_transform()
    return [
        inv(sqrt(3.0)) inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) -inv(sqrt(2.0)) inv(sqrt(6.0))
        inv(sqrt(3.0)) 0.0 -2.0 * inv(sqrt(6.0))
    ]
end

function _steady_state_terminal_phasors(sample, node_indices::AbstractVector{<:Integer})
    return ComplexF64[
        node == 0 ? complex(0.0, 0.0) : sample.node_voltage_phasors[Int(node)]
        for node in node_indices
    ]
end

function _distributed_transposed_line_steady_state_history_waves(
    modal_state,
    sample,
)
    phase_count = modal_state.phase_count
    from_voltage =
        _steady_state_terminal_phasors(sample, modal_state.from_node_indices)
    to_voltage =
        _steady_state_terminal_phasors(sample, modal_state.to_node_indices)
    modal_transform = _distributed_transposed_line_modal_transform()
    modal_voltage_from =
        ComplexF64.(phase_count .* (transpose(modal_transform) * from_voltage))
    modal_voltage_to =
        ComplexF64.(phase_count .* (transpose(modal_transform) * to_voltage))
    coefficients = distributed_transposed_line_norton_coefficients(modal_state)
    angular_frequency = 2.0 * pi * Float64(sample.steady_state_frequency_hz)
    outgoing_wave_from = ComplexF64[]
    outgoing_wave_to = ComplexF64[]
    for mode_index in 1:modal_state.group_length
        admittance = coefficients.modal_companion_admittance_values[mode_index]
        damping = coefficients.modal_history_damping_values[mode_index]
        phase_shift = cis(
            -angular_frequency * (
                modal_state.modal_propagation_times_s[mode_index] -
                Float64(sample.timestep_s)
            ),
        )
        delayed_damping = damping * phase_shift
        denominator = 1.0 - delayed_damping^2
        push!(
            outgoing_wave_from,
            admittance * (
                modal_voltage_to[mode_index] -
                delayed_damping * modal_voltage_from[mode_index]
            ) / denominator,
        )
        push!(
            outgoing_wave_to,
            admittance * (
                modal_voltage_from[mode_index] -
                delayed_damping * modal_voltage_to[mode_index]
            ) / denominator,
        )
    end
    return (
        outgoing_wave_from = outgoing_wave_from,
        outgoing_wave_to = outgoing_wave_to,
    )
end

function _distributed_transposed_line_modal_component_terms(phase_terms, phase_count::Int)
    modal_transform = _distributed_transposed_line_modal_transform()
    term_modal_indices = Int[]
    term_weight_values = Float64[]
    term_from_real_values = Float64[]
    term_from_imag_values = Float64[]
    term_to_real_values = Float64[]
    term_to_imag_values = Float64[]
    for mode_index in 1:phase_count
        for phase_index in 1:phase_count
            push!(term_modal_indices, mode_index)
            push!(term_weight_values, modal_transform[phase_index, mode_index])
            push!(term_from_real_values, phase_terms.terminal_from_real_values[phase_index])
            push!(term_from_imag_values, phase_terms.terminal_from_imag_values[phase_index])
            push!(term_to_real_values, phase_terms.terminal_to_real_values[phase_index])
            push!(term_to_imag_values, phase_terms.terminal_to_imag_values[phase_index])
        end
    end
    return (
        term_modal_indices = term_modal_indices,
        term_weight_values = term_weight_values,
        term_from_real_values = term_from_real_values,
        term_from_imag_values = term_from_imag_values,
        term_to_real_values = term_to_real_values,
        term_to_imag_values = term_to_imag_values,
    )
end

function _distributed_transposed_line_steady_state_initial_history(
    modal_state,
    pi_equivalent,
    sample,
    timestep_s::Real,
    history_storage_start_index::Integer,
    index::Integer,
)
    phase_count = modal_state.phase_count
    from_voltage =
        _steady_state_terminal_phasors(sample, modal_state.from_node_indices)
    to_voltage =
        _steady_state_terminal_phasors(sample, modal_state.to_node_indices)
    series_admittance = inv(pi_equivalent.phase_series_impedance_matrix)
    shunt_admittance = pi_equivalent.phase_shunt_admittance_matrix
    terminal_current_from =
        series_admittance * (from_voltage - to_voltage) +
        shunt_admittance * from_voltage
    terminal_current_to =
        series_admittance * (to_voltage - from_voltage) +
        shunt_admittance * to_voltage
    phase_matrix = distributed_transposed_line_initial_history_phase_matrix(
        modal_state;
        name = Symbol(
            "distributed_transposed_line_initial_history_phase_matrix_",
            index,
        ),
    )
    phase_terms = distributed_transposed_line_initial_history_phase_terms(
        reduced_phase_matrix_upper_values = phase_matrix.phase_matrix_upper_values,
        matrix_terminal_from_real_values = real.(terminal_current_to),
        matrix_terminal_from_imag_values = imag.(terminal_current_to),
        matrix_terminal_to_real_values = real.(terminal_current_from),
        matrix_terminal_to_imag_values = imag.(terminal_current_from),
        initial_terminal_from_real_values = real.(to_voltage),
        initial_terminal_from_imag_values = imag.(to_voltage),
        initial_terminal_to_real_values = real.(from_voltage),
        initial_terminal_to_imag_values = imag.(from_voltage),
        phase_indices = modal_state.phase_indices,
        name = Symbol(
            "distributed_transposed_line_initial_history_phase_terms_",
            index,
        ),
    )
    component_terms =
        _distributed_transposed_line_modal_component_terms(phase_terms, phase_count)
    modal_indices = collect(1:modal_state.group_length)
    components = distributed_transposed_line_initial_history_modal_components(
        modal_indices = modal_indices,
        term_modal_indices = component_terms.term_modal_indices,
        term_weight_values = component_terms.term_weight_values,
        term_from_real_values = component_terms.term_from_real_values,
        term_from_imag_values = component_terms.term_from_imag_values,
        term_to_real_values = component_terms.term_to_real_values,
        term_to_imag_values = component_terms.term_to_imag_values,
        name = Symbol(
            "distributed_transposed_line_initial_history_modal_components_",
            index,
        ),
    )
    phasors = distributed_transposed_line_initial_history_phasors(
        modal_scale_values = phase_matrix.modal_history_scale_values,
        modal_component_from_real_values =
            components.modal_component_from_real_values,
        modal_component_from_imag_values =
            components.modal_component_from_imag_values,
        modal_component_to_real_values =
            components.modal_component_to_real_values,
        modal_component_to_imag_values =
            components.modal_component_to_imag_values,
        modal_indices = modal_indices,
        name = Symbol("distributed_transposed_line_initial_history_phasors_", index),
    )
    norton_coefficients = distributed_transposed_line_norton_coefficients(modal_state)
    delay_counts =
        floor.(Int, modal_state.modal_propagation_times_s ./ Float64(timestep_s))
    seed = distributed_transposed_line_initial_history_seed(
        timestep_s = timestep_s,
        steady_state_frequency_hz = sample.steady_state_frequency_hz,
        modal_delay_counts = delay_counts,
        outgoing_wave_from_amplitudes = phasors.modal_history_from_amplitudes,
        outgoing_wave_from_phases = phasors.modal_history_from_phases,
        outgoing_wave_to_amplitudes = phasors.modal_history_to_amplitudes,
        outgoing_wave_to_phases = phasors.modal_history_to_phases,
        modal_history_interpolation_factors =
            modal_state.modal_propagation_times_s ./ Float64(timestep_s) .-
            delay_counts,
        modal_history_damping_values =
            norton_coefficients.modal_history_damping_values,
        modal_signed_characteristic_values =
            modal_state.modal_signed_characteristic_impedances,
        history_storage_start_index = history_storage_start_index,
        modal_indices = modal_indices,
        name = Symbol("distributed_transposed_line_initial_history_seed_", index),
    )
    history = distributed_transposed_line_history_state(
        modal_state,
        seed;
        timestep_s = timestep_s,
        steady_state_frequency_hz = sample.steady_state_frequency_hz,
        initialized_from_steady_state = true,
        name = Symbol("distributed_transposed_line_history_state_", index),
    )
    return (
        phase_matrix = phase_matrix,
        phase_terms = phase_terms,
        modal_components = components,
        phasors = phasors,
        seed = seed,
        history = history,
        terminal_current_from = terminal_current_from,
        terminal_current_to = terminal_current_to,
    )
end

function _shared_distributed_transposed_line_history_states(
    histories::AbstractVector{DistributedTransposedLineHistoryState},
)
    isempty(histories) && return DistributedTransposedLineHistoryState[]
    global_start = minimum(history.history_storage_start_index for history in histories)
    global_end = maximum(maximum(history.storage_end_indices) for history in histories)
    shared_from = zeros(Float64, global_end - global_start + 1)
    shared_to = zeros(Float64, global_end - global_start + 1)

    for history in histories
        for mode_index in eachindex(history.storage_start_indices)
            start_slot = history.storage_start_indices[mode_index]
            first_slot = start_slot == global_start ? start_slot : start_slot + 1
            last_slot = history.storage_end_indices[mode_index]
            first_slot <= last_slot || continue
            for slot in first_slot:last_slot
                source_offset = slot - history.history_storage_start_index + 1
                1 <= source_offset <= length(history.packed_history_from_values) ||
                    throw(ArgumentError("distributed-line packed history source slot is outside the local storage vector"))
                shared_offset = slot - global_start + 1
                shared_from[shared_offset] = history.packed_history_from_values[source_offset]
                shared_to[shared_offset] = history.packed_history_to_values[source_offset]
            end
        end
    end

    shared_histories = DistributedTransposedLineHistoryState[]
    for history in histories
        push!(
            shared_histories,
            DistributedTransposedLineHistoryState(
                history.name,
                copy(history.phase_indices),
                copy(history.line_numbers),
                copy(history.modal_sequence_indices),
                history.timestep_s,
                history.steady_state_frequency_hz,
                history.angular_step_rad,
                global_start,
                global_end,
                copy(history.storage_start_indices),
                copy(history.storage_end_indices),
                copy(history.storage_lengths),
                copy(history.history_sample_counts),
                copy(history.modal_history_interpolation_factors),
                copy(history.history_read_indices),
                copy(history.outgoing_wave_from_real_values),
                copy(history.outgoing_wave_from_imag_values),
                copy(history.outgoing_wave_to_real_values),
                copy(history.outgoing_wave_to_imag_values),
                [copy(values) for values in history.modal_history_from_values],
                [copy(values) for values in history.modal_history_to_values],
                shared_from,
                shared_to,
                history.initialized_from_steady_state,
                history.phase_count,
            ),
        )
    end
    return shared_histories
end

function _deck_distributed_transposed_line_config(
    parsed::DeckParser.DeckParseResult;
    steady_state_initial_sample = nothing,
)
    modal_states = DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)
    isempty(modal_states) && return nothing
    pi_equivalents =
        DeckParser.deck_distributed_transposed_line_steady_state_pi_equivalents(parsed)
    if steady_state_initial_sample !== nothing &&
       length(pi_equivalents) != length(modal_states)
        throw(ArgumentError("steady-state distributed-line initialization requires a PI equivalent for each modal line state"))
    end
    history_states =
        steady_state_initial_sample === nothing ?
        DeckParser.deck_distributed_transposed_line_history_states(
            parsed;
            first_history_storage_index = 1,
        ) :
        begin
            options = DeckParser.deck_fixed_time_horizon_options(parsed)
            storage_index = 1
            histories = DistributedTransposedLineHistoryState[]
            initial_history_results = Any[]
            for index in eachindex(modal_states)
                initial_history =
                    _distributed_transposed_line_steady_state_initial_history(
                    modal_states[index],
                    pi_equivalents[index],
                    steady_state_initial_sample,
                    options.dt_s,
                    storage_index,
                    index,
                )
                history = initial_history.history
                push!(histories, history)
                push!(initial_history_results, initial_history)
                storage_index = history.next_history_storage_index
            end
            _shared_distributed_transposed_line_history_states(histories)
        end
    length(history_states) == length(modal_states) ||
        throw(ArgumentError("distributed transposed line modal/history state counts must match"))
    return (
        enabled = true,
        modal_states = modal_states,
        history_states = history_states,
        steady_state_initial_history_results =
            steady_state_initial_sample === nothing ? Any[] : initial_history_results,
        current_injection_values = Float64[],
    )
end

_deck_has_dynamic_distributed_line(parsed::DeckParser.DeckParseResult) =
    !isempty(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed))

_deck_has_frequency_dependent_line_runtime(parsed::DeckParser.DeckParseResult) =
    !isempty(DeckParser.deck_sampled_frequency_line_rows(parsed)) ||
    !isempty(DeckParser.deck_semlyen_line_groups(parsed)) ||
    !isempty(DeckParser.deck_rational_frequency_line_groups(parsed))

function _primary_dynamic_timestep_state(
    node_count::Integer,
    dt_s::Real,
    switch_count::Integer,
    ;
    source_card_read_requested::Bool=true,
)
    node_total = Int(node_count)
    node_total > 0 || throw(ArgumentError("deck dynamic timestep state requires at least one node"))
    switches = Int(switch_count)
    switches >= 0 || throw(ArgumentError("deck switch count must be nonnegative"))
    shifted_node_total = node_total + 1
    return OVER16AcceptedTimestepState(
        OVER16SourceUpdateState(
            fill(0.0, 10),
            zeros(node_total),
            zeros(node_total);
            iread = source_card_read_requested ? 1 : 0,
        ),
        OVER16TACSUtilityState(zeros(1)),
        OVER16OutputReportState(peaknd_values = [0.0, 0.0, 0.0]),
        OVER16PostExtremaControlState(0.0, Float64(dt_s), 0),
        OVER16SwitchScanState(fill(-1, switches), zeros(switches)),
        OVER16SwitchOperationState(Int[], 0),
        OVER16SwitchTopologyState(
            fill(false, switches);
            nextsw = zeros(Int, switches),
            kode = zeros(Int, shifted_node_total),
        ),
        OVER16SwitchCurrentState(zeros(shifted_node_total), zeros(switches)),
        OVER16SwitchPostCurrentState(fill(-1, switches), zeros(switches), zeros(switches)),
        OVER16SwitchBValueExportState(Float64[]),
        OVER16SwitchAlterationState(0, 0, 0, 0),
    )
end

function _deck_dynamic_timestep_state(
    parsed::DeckParser.DeckParseResult,
    node_count::Integer,
    dt_s::Real,
)
    return _primary_dynamic_timestep_state(
        node_count,
        dt_s,
        length(DeckParser.deck_time_switch_names(parsed));
        source_card_read_requested =
            !isempty(DeckParser.deck_over16_source_card_rows(parsed)),
    )
end

function run_nested_cable_primary_timestep!(
    context::EMTStepContext,
    line_state::NestedCableTransientLineState;
    from_node_indices,
    to_node_indices,
    target_frequency_hz::Real,
    current_projection::Symbol = :harmonic_phasor,
    current_projection_angle_rad::Real = 0.0,
    over16_state::Union{Nothing,OVER16AcceptedTimestepState} = nothing,
)
    abs(context.dt_s - line_state.dt_s) <=
        64.0 * eps(Float64) * max(context.dt_s, line_state.dt_s) ||
        throw(ArgumentError("nested cable state dt_s must match the primary timestep context"))
    line_state.physical_checks_passed ||
        throw(ArgumentError("nested cable transient state failed its physical checks"))
    state = over16_state === nothing ?
        _primary_dynamic_timestep_state(context.system.node_count, context.dt_s, 0) :
        over16_state
    current_values = Float64[]
    return step_with_over16_boundary!(
        context,
        state;
        frequency_dependent_line_config = (
            enabled = true,
            nested_cable_states = [line_state],
            from_node_indices = from_node_indices,
            to_node_indices = to_node_indices,
            target_frequency_hz = Float64(target_frequency_hz),
            current_projection = current_projection,
            current_projection_angle_rad = Float64(current_projection_angle_rad),
            current_injection_values = current_values,
        ),
    )
end

function _terminal_voltage_values(
    voltage::AbstractVector{Float64},
    node_indices::AbstractVector{Int},
)
    values = Float64[]
    for node in node_indices
        if node == 0
            push!(values, 0.0)
        else
            1 <= node <= length(voltage) ||
                throw(ArgumentError("distributed transposed line node index exceeds solved voltage length"))
            push!(values, voltage[node])
        end
    end
    return values
end

function _frequency_dependent_line_recursive_states(config::NamedTuple)
    if haskey(config, :nested_cable_states)
        cable_states = getproperty(config, :nested_cable_states)
        all(state -> state isa NestedCableTransientLineState, cable_states) ||
            throw(ArgumentError("nested_cable_states must contain NestedCableTransientLineState values"))
        return [state.recursive_state for state in cable_states]
    elseif haskey(config, :recursive_convolution_states)
        return getproperty(config, :recursive_convolution_states)
    elseif haskey(config, :states)
        return getproperty(config, :states)
    end
    throw(ArgumentError(
        "frequency_dependent_line_config requires nested_cable_states or recursive_convolution_states",
    ))
end

function _frequency_dependent_line_projection_tolerance(config::NamedTuple)
    tolerance = haskey(config, :real_current_projection_tolerance) ?
        Float64(getproperty(config, :real_current_projection_tolerance)) :
        1.0e-9
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("real_current_projection_tolerance must be finite and nonnegative"))
    return tolerance
end

function _frequency_dependent_line_current_projection(config::NamedTuple)
    kind = Symbol(get(config, :current_projection, :strict_real))
    kind in (:strict_real, :harmonic_phasor) ||
        throw(ArgumentError("current_projection must be :strict_real or :harmonic_phasor"))
    angle = Float64(get(config, :current_projection_angle_rad, 0.0))
    isfinite(angle) || throw(ArgumentError("current_projection_angle_rad must be finite"))
    return kind, angle
end

function _frequency_dependent_line_node_group(
    values,
    index::Int,
    line_count::Int,
    phase_count::Int,
    node_count::Int,
    label::AbstractString,
)
    raw_nodes =
        values isa AbstractVector{<:Integer} ?
        begin
            line_count == 1 ||
                throw(ArgumentError("$label must contain one node-index vector per frequency-dependent line"))
            values
        end :
        begin
            length(values) == line_count ||
                throw(ArgumentError("$label count must match frequency-dependent line count"))
            values[index]
        end
    nodes = Int.(collect(raw_nodes))
    length(nodes) == phase_count ||
        throw(ArgumentError("$label phase count must match frequency-dependent line state"))
    for node in nodes
        0 <= node <= node_count ||
            throw(ArgumentError("$label contains a node index outside the solved network"))
    end
    return nodes
end

function _frequency_dependent_line_target_frequency(
    config::NamedTuple,
    index::Int,
    line_count::Int,
)
    if haskey(config, :target_frequency_hz_values)
        values = getproperty(config, :target_frequency_hz_values)
        length(values) == line_count ||
            throw(ArgumentError("target_frequency_hz_values count must match frequency-dependent line count"))
        frequency = Float64(values[index])
    elseif haskey(config, :frequency_hz_values)
        values = getproperty(config, :frequency_hz_values)
        length(values) == line_count ||
            throw(ArgumentError("frequency_hz_values count must match frequency-dependent line count"))
        frequency = Float64(values[index])
    elseif haskey(config, :target_frequency_hz)
        frequency = Float64(getproperty(config, :target_frequency_hz))
    elseif haskey(config, :frequency_hz)
        frequency = Float64(getproperty(config, :frequency_hz))
    else
        throw(ArgumentError("frequency_dependent_line_config requires target_frequency_hz or frequency_hz"))
    end
    isfinite(frequency) && frequency >= 0.0 ||
        throw(ArgumentError("frequency-dependent line target frequency must be finite and nonnegative"))
    return frequency
end

function _frequency_dependent_line_real_currents(
    values::AbstractVector,
    tolerance::Float64,
    label::AbstractString,
    projection_kind::Symbol,
    projection_angle_rad::Float64,
)
    currents = Float64[]
    sizehint!(currents, length(values))
    max_imag = 0.0
    rotation = cis(projection_angle_rad)
    for (index, value) in pairs(values)
        current = ComplexF64(value) * rotation
        imag_abs = abs(imag(current))
        max_imag = max(max_imag, imag_abs)
        (projection_kind == :harmonic_phasor || imag_abs <= tolerance) ||
            throw(ArgumentError("$label current $index has non-real residual $imag_abs"))
        push!(currents, real(current))
    end
    return currents, max_imag
end

function _frequency_dependent_line_phase_voltages(
    voltage::AbstractVector{Float64},
    nodes::AbstractVector{Int},
)
    return ComplexF64.(_terminal_voltage_values(voltage, nodes))
end

function _frequency_dependent_line_deck_current_injection!(
    config,
    rhs::AbstractVector{Float64},
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("frequency_dependent_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    states = _frequency_dependent_line_recursive_states(config)
    line_count = length(states)
    line_count > 0 ||
        throw(ArgumentError("frequency_dependent_line_config requires at least one recursive state"))
    haskey(config, :from_node_indices) ||
        throw(ArgumentError("frequency_dependent_line_config requires from_node_indices"))
    haskey(config, :to_node_indices) ||
        throw(ArgumentError("frequency_dependent_line_config requires to_node_indices"))
    tolerance = _frequency_dependent_line_projection_tolerance(config)
    projection_kind, projection_angle =
        _frequency_dependent_line_current_projection(config)

    rhs_before = copy(rhs)
    rhs_update_count = 0
    sending_currents = Vector{Vector{Float64}}()
    receiving_currents = Vector{Vector{Float64}}()
    from_node_groups = Vector{Vector{Int}}()
    to_node_groups = Vector{Vector{Int}}()
    max_imag = 0.0
    for index in 1:line_count
        state = states[index]
        phase_count = length(state.sending_phase_current)
        length(state.receiving_phase_current) == phase_count ||
            throw(ArgumentError("frequency-dependent line sending/receiving current counts must match"))
        from_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :from_node_indices),
            index,
            line_count,
            phase_count,
            length(rhs),
            "frequency-dependent from_node_indices",
        )
        to_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :to_node_indices),
            index,
            line_count,
            phase_count,
            length(rhs),
            "frequency-dependent to_node_indices",
        )
        sending, sending_imag = _frequency_dependent_line_real_currents(
            state.sending_phase_current,
            tolerance,
            "frequency-dependent sending",
            projection_kind,
            projection_angle,
        )
        receiving, receiving_imag = _frequency_dependent_line_real_currents(
            state.receiving_phase_current,
            tolerance,
            "frequency-dependent receiving",
            projection_kind,
            projection_angle,
        )
        max_imag = max(max_imag, sending_imag, receiving_imag)
        for phase_index in 1:phase_count
            from_node = from_nodes[phase_index]
            if from_node != 0
                rhs[from_node] += sending[phase_index]
                rhs_update_count += 1
            end
            to_node = to_nodes[phase_index]
            if to_node != 0
                rhs[to_node] += receiving[phase_index]
                rhs_update_count += 1
            end
        end
        push!(sending_currents, sending)
        push!(receiving_currents, receiving)
        push!(from_node_groups, from_nodes)
        push!(to_node_groups, to_nodes)
    end
    if haskey(config, :current_injection_values)
        values = getproperty(config, :current_injection_values)
        values isa Vector{Float64} ||
            throw(ArgumentError("frequency-dependent current_injection_values must be Vector{Float64}"))
        resize!(values, length(rhs))
        for index in eachindex(rhs)
            values[index] = rhs[index] - rhs_before[index]
        end
    end
    return (
        source = :frequency_dependent_line_deck_current_injection,
        outcome = :history_current_injection,
        line_update_count = 0,
        rhs_update_count = rhs_update_count,
        rhs_before_values = rhs_before,
        rhs_after_values = copy(rhs),
        sending_phase_current_values = sending_currents,
        receiving_phase_current_values = receiving_currents,
        from_node_indices = from_node_groups,
        to_node_indices = to_node_groups,
        real_current_projection_max_imag_abs = max_imag,
        current_projection = projection_kind,
        current_projection_angle_rad = projection_angle,
        nested_cable_frequency_state_consumed = haskey(config, :nested_cable_states),
        state_mutated = false,
    )
end

function _frequency_dependent_line_deck_history_update!(
    config,
    voltage::AbstractVector{Float64},
    current_result,
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("frequency_dependent_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    states = _frequency_dependent_line_recursive_states(config)
    line_count = length(states)
    line_count > 0 ||
        throw(ArgumentError("frequency_dependent_line_config requires at least one recursive state"))
    updates = Any[]
    target_frequencies = Float64[]
    for index in 1:line_count
        state = states[index]
        phase_count = length(state.sending_phase_current)
        from_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :from_node_indices),
            index,
            line_count,
            phase_count,
            length(voltage),
            "frequency-dependent from_node_indices",
        )
        to_nodes = _frequency_dependent_line_node_group(
            getproperty(config, :to_node_indices),
            index,
            line_count,
            phase_count,
            length(voltage),
            "frequency-dependent to_node_indices",
        )
        target_frequency = _frequency_dependent_line_target_frequency(config, index, line_count)
        update = frequency_dependent_line_recursive_convolution_update!(
            state,
            target_frequency,
            _frequency_dependent_line_phase_voltages(voltage, from_nodes),
            _frequency_dependent_line_phase_voltages(voltage, to_nodes),
        )
        push!(updates, update)
        push!(target_frequencies, target_frequency)
    end
    recursive_runtime_executed =
        any(update -> get(update, :bounded_recursive_convolution_runtime_executed, false), updates)
    skin_effect_internal_impedance_executed =
        any(update -> get(update, :skin_effect_internal_impedance_executed, false), updates)
    earth_return_impedance_executed =
        any(update -> get(update, :earth_return_impedance_executed, false), updates)
    frequency_dependent_fitting_executed =
        any(update -> get(update, :frequency_dependent_fitting_executed, false), updates)
    frequency_loop_executed =
        any(update -> get(update, :frequency_loop_executed, false), updates)
    pipe_sheath_side_effects_executed =
        any(update -> get(update, :pipe_sheath_side_effects_executed, false), updates)
    deferred_effects = Symbol[]
    frequency_dependent_fitting_executed ||
        push!(deferred_effects, :full_bpa_frequency_dependent_fitting)
    pipe_sheath_side_effects_executed ||
        push!(deferred_effects, :full_over47_pipe_or_sheath_coupling)
    frequency_loop_executed ||
        push!(deferred_effects, :full_over47_frequency_loop_side_effects)
    haskey(config, :nested_cable_states) ||
        push!(deferred_effects, :frequency_dependent_line_deck_grammar)
    return (
        source = :frequency_dependent_line_deck_timestep_update,
        outcome = :timestep_integration,
        line_update_count = length(updates),
        rhs_update_count = current_result === nothing ? 0 : current_result.rhs_update_count,
        recursive_convolution_updates = updates,
        phase_current_injections =
            current_result === nothing ? Any[] : [current_result],
        rhs_before_values =
            current_result === nothing ? Float64[] : current_result.rhs_before_values,
        rhs_after_values =
            current_result === nothing ? Float64[] : current_result.rhs_after_values,
        target_frequency_hz_values = target_frequencies,
        frequency_dependent_line_runtime_executed = !isempty(updates),
        recursive_convolution_runtime_executed =
            recursive_runtime_executed,
        skin_effect_internal_impedance_executed =
            skin_effect_internal_impedance_executed,
        earth_return_impedance_executed =
            earth_return_impedance_executed,
        frequency_dependent_fitting_executed =
            frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed =
            pipe_sheath_side_effects_executed,
        nested_cable_frequency_state_consumed = haskey(config, :nested_cable_states),
        state_mutated = !isempty(updates),
        deferred_effects = Tuple(deferred_effects),
    )
end

function _distributed_transposed_line_deck_current_injection!(
    config,
    rhs::AbstractVector{Float64},
    ;
    collect_diagnostics::Bool = true,
)
    config === nothing && return nothing
    config isa NamedTuple ||
        throw(ArgumentError("distributed_transposed_line_config must be a NamedTuple"))
    get(config, :enabled, true) || return nothing
    haskey(config, :modal_states) ||
        throw(ArgumentError("distributed_transposed_line_config requires modal_states"))
    haskey(config, :history_states) ||
        throw(ArgumentError("distributed_transposed_line_config requires history_states"))
    modal_states = config.modal_states
    history_states = config.history_states
    length(modal_states) == length(history_states) ||
        throw(ArgumentError("distributed transposed line modal/history state counts must match"))
    injections = collect_diagnostics ? Any[] : nothing
    rhs_before = copy(rhs)
    rhs_update_count = 0
    for index in eachindex(modal_states)
        modal_state = modal_states[index]
        history_state = history_states[index]
        injection = distributed_transposed_line_history_current_injection!(
            rhs,
            modal_state,
            history_state;
            name = collect_diagnostics ?
                Symbol("distributed_transposed_line_history_current_injection_", index) :
                history_state.name,
            collect_diagnostics = collect_diagnostics,
        )
        if collect_diagnostics
            rhs_update_count += injection.rhs_update_count
            push!(injections, injection)
        end
    end
    if haskey(config, :current_injection_values)
        current_injection_values = getproperty(config, :current_injection_values)
        current_injection_values isa Vector{Float64} ||
            throw(ArgumentError("distributed line current_injection_values must be Vector{Float64}"))
        resize!(current_injection_values, length(rhs))
        for index in eachindex(rhs)
            current_injection_values[index] = rhs[index] - rhs_before[index]
        end
    end
    if !collect_diagnostics
        return (
            rhs_before_values = rhs_before,
            rhs_after_values = rhs,
        )
    end
    return (
        source = :distributed_transposed_line_deck_current_injection,
        outcome = :history_current_injection,
        line_update_count = 0,
        rhs_update_count = rhs_update_count,
        modal_timestep_updates = Any[],
        phase_current_injections = injections,
        rhs_before_values = rhs_before,
        rhs_after_values = copy(rhs),
        state_mutated = false,
    )
end

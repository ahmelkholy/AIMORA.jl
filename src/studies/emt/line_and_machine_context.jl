
function _distributed_transposed_line_deck_history_update!(
    config,
    voltage::AbstractVector{Float64},
    current_result,
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
    updates = collect_diagnostics ? Any[] : nothing
    for index in eachindex(modal_states)
        modal_state = modal_states[index]
        history_state = history_states[index]
        from_voltage = _terminal_voltage_values(voltage, modal_state.from_node_indices)
        to_voltage = _terminal_voltage_values(voltage, modal_state.to_node_indices)
        update = distributed_transposed_line_modal_timestep_update!(
            modal_state,
            history_state,
            from_voltage,
            to_voltage,
            name = collect_diagnostics ?
                Symbol(string(history_state.name), "_modal_timestep_update") :
                history_state.name,
            collect_diagnostics = collect_diagnostics,
        )
        collect_diagnostics && push!(updates, update)
    end
    collect_diagnostics || return nothing
    return (
        source = :distributed_transposed_line_deck_timestep_update,
        outcome = :timestep_integration,
        line_update_count = length(updates),
        rhs_update_count =
            current_result === nothing ? 0 : current_result.rhs_update_count,
        modal_timestep_updates = updates,
        phase_current_injections =
            current_result === nothing ? Any[] : current_result.phase_current_injections,
        rhs_before_values =
            current_result === nothing ? Float64[] : current_result.rhs_before_values,
        rhs_after_values =
            current_result === nothing ? Float64[] : current_result.rhs_after_values,
        state_mutated = !isempty(updates),
    )
end

function _source_config_field(config::NamedTuple, key::Symbol)
    haskey(config, key) ||
        throw(ArgumentError("source_config requires $key"))
    return getproperty(config, key)
end

function _electromagnetic_source_voltage_constraints(
    over16_state::OVER16AcceptedTimestepState,
    source_config,
    node_count::Int,
    solve_time_s::Float64,
    use_accepted_source_signal_state::Bool,
)
    empty_result = (
        source = :electromagnetic_source_voltage_constraints,
        outcome = :integration_boundary,
        nodes = Int[],
        values = Float64[],
        row_indices = Int[],
        count = 0,
        ideal_transformer_rhs_nodes = Int[],
        ideal_transformer_rhs_values = Float64[],
        ideal_transformer_rhs_row_indices = Int[],
        ideal_transformer_rhs_count = 0,
        preview = nothing,
        deferred_row_count = 0,
        applied = false,
    )
    source_config === nothing && return empty_result
    source_config isa NamedTuple ||
        throw(ArgumentError("source_config must be a NamedTuple"))
    source_kwargs = get(source_config, :kwargs, NamedTuple())
    source_kwargs isa NamedTuple ||
        throw(ArgumentError("source_config kwargs must be a NamedTuple"))
    # A typed source-signal program owns its initial accepted state. Each
    # STEP/CALC/ELEC pass solves with that state, then SUBTS3 composes the
    # next reader, interpolation, TACS, and analytic values after the solve.
    # Existing deck-card boundary rows retain their accepted same-pass timing.
    voltage_values =
        use_accepted_source_signal_state ?
        over16_state.source.source_card.voltbc_values :
        get(
            source_kwargs,
            :interpolated_values,
            over16_state.source.source_card.voltbc_values,
        )
    preview = over16_source_row_update_preview(
        over16_state.source.f_values,
        over16_state.source.e_values,
        _source_config_field(source_config, :node_values),
        _source_config_field(source_config, :iform_values),
        _source_config_field(source_config, :tstart_values),
        _source_config_field(source_config, :tstop_values),
        voltage_values,
        solve_time_s;
        kconst = get(
            source_kwargs,
            :kconst,
            length(_source_config_field(source_config, :node_values)),
        ),
        crest_values = get(source_kwargs, :crest_values, nothing),
        time1_values = get(source_kwargs, :time1_values, nothing),
        time2_values = get(source_kwargs, :time2_values, nothing),
        sfreq_values = get(source_kwargs, :sfreq_values, nothing),
        xtcs_values = get(source_kwargs, :xtcs_values, Float64[]),
        kxtcs = get(source_kwargs, :kxtcs, 0),
        delta2 = get(source_kwargs, :delta2, nothing),
        initial_source_value = get(source_kwargs, :initial_source_value, nothing),
    )
    constraint_nodes = Int[]
    constraint_values = Float64[]
    constraint_rows = Int[]
    ideal_transformer_predecessor_rows = Set(
        row - 1 for row in preview.ideal_transformer_source_rows if row > 1
    )
    for index in eachindex(preview.e_update_indices, preview.e_update_values)
        preview.e_update_rows[index] in ideal_transformer_predecessor_rows && continue
        node = preview.e_update_indices[index]
        1 <= node <= node_count || continue
        existing = findfirst(==(node), constraint_nodes)
        if existing === nothing
            push!(constraint_nodes, node)
            push!(constraint_values, Float64(preview.e_update_values[index]))
            push!(constraint_rows, preview.e_update_rows[index])
        else
            constraint_values[existing] = Float64(preview.e_update_values[index])
            constraint_rows[existing] = preview.e_update_rows[index]
        end
    end
    ideal_transformer_rhs_nodes = Int[]
    ideal_transformer_rhs_values = Float64[]
    ideal_transformer_rhs_rows = Int[]
    ideal_transformer_rows = Set(preview.ideal_transformer_source_rows)
    for index in eachindex(preview.f_update_indices, preview.f_update_values)
        row = preview.f_update_rows[index]
        row in ideal_transformer_rows || continue
        node = preview.f_update_indices[index]
        1 <= node <= node_count || continue
        existing = findfirst(==(node), ideal_transformer_rhs_nodes)
        if existing === nothing
            push!(ideal_transformer_rhs_nodes, node)
            push!(ideal_transformer_rhs_values, Float64(preview.f_update_values[index]))
            push!(ideal_transformer_rhs_rows, row)
        else
            ideal_transformer_rhs_values[existing] =
                Float64(preview.f_update_values[index])
            ideal_transformer_rhs_rows[existing] = row
        end
    end
    return (
        source = :electromagnetic_source_voltage_constraints,
        outcome = :integration_boundary,
        nodes = constraint_nodes,
        values = constraint_values,
        row_indices = constraint_rows,
        count = length(constraint_nodes),
        ideal_transformer_rhs_nodes = ideal_transformer_rhs_nodes,
        ideal_transformer_rhs_values = ideal_transformer_rhs_values,
        ideal_transformer_rhs_row_indices = ideal_transformer_rhs_rows,
        ideal_transformer_rhs_count = length(ideal_transformer_rhs_nodes),
        preview = preview,
        deferred_row_count = preview.deferred_row_count,
        applied = !isempty(constraint_nodes),
    )
end

function _apply_ideal_transformer_rhs_values!(
    rhs::Vector{Float64},
    source_voltage_constraint_result,
)
    for (node, value) in zip(
        source_voltage_constraint_result.ideal_transformer_rhs_nodes,
        source_voltage_constraint_result.ideal_transformer_rhs_values,
    )
        rhs[node] = Float64(value)
    end
    return rhs
end

function _source_voltage_constraint_time(source_config, fallback_time_s::Float64)
    source_config isa NamedTuple && haskey(source_config, :t) ||
        return fallback_time_s
    if haskey(source_config, :constraint_t)
        return Float64(getproperty(source_config, :constraint_t))
    end
    return Float64(getproperty(source_config, :t))
end

function _apply_voltage_constraint_values!(
    voltage::Vector{Float64},
    source_voltage_constraint_result,
)
    get(source_voltage_constraint_result, :applied, false) || return voltage
    for (node, value) in zip(
        source_voltage_constraint_result.nodes,
        source_voltage_constraint_result.values,
    )
        1 <= node <= length(voltage) || continue
        voltage[node] = Float64(value)
    end
    return voltage
end

function _source_constrained_presolve_voltage(
    voltage::Vector{Float64},
    source_voltage_constraint_result,
    presolve_voltage::Vector{Float64},
    record_presolve_voltage_state::Bool,
)
    if record_presolve_voltage_state
        _apply_voltage_constraint_values!(
            presolve_voltage,
            source_voltage_constraint_result,
        )
        return presolve_voltage
    end
    get(source_voltage_constraint_result, :applied, false) || return voltage
    constrained_voltage = copy(voltage)
    return _apply_voltage_constraint_values!(
        constrained_voltage,
        source_voltage_constraint_result,
    )
end

function _electromagnetic_rhs_current_injections(
    over16_state::OVER16AcceptedTimestepState,
    node_count::Int,
    raw_over16_kwargs::NamedTuple,
)
    if haskey(raw_over16_kwargs, :current_injection_values)
        configured = Float64.(getproperty(raw_over16_kwargs, :current_injection_values))
        length(configured) >= node_count ||
            throw(ArgumentError("current_injection_values must cover every node"))
        return configured[1:node_count]
    end
    return zeros(Float64, node_count)
end

function _nonlinear_current_compensation_injection_values(
    over16_state::OVER16AcceptedTimestepState,
    node_count::Int,
)
    values = zeros(Float64, node_count)
    source_values = over16_state.source.nonlinear_current_compensation_values
    shared_count = min(node_count, length(source_values))
    for index in 1:shared_count
        values[index] = source_values[index]
    end
    return values
end

function _add_current_injection_values!(
    current_injections::Vector{Float64},
    values::AbstractVector{<:Real},
)
    length(values) >= length(current_injections) ||
        throw(ArgumentError("current injection values must cover current injections"))
    for index in eachindex(current_injections)
        current_injections[index] += Float64(values[index])
    end
    return current_injections
end

function _add_current_injection_prefix_values!(
    current_injections::Vector{Float64},
    values::AbstractVector{<:Real},
)
    shared_count = min(length(current_injections), length(values))
    for index in 1:shared_count
        current_injections[index] += Float64(values[index])
    end
    return current_injections
end

function _source_rhs_without_nonlinear_current_compensation(
    rhs_values::AbstractVector{<:Real},
    compensation_values::AbstractVector{<:Real},
)
    base = Float64.(rhs_values)
    shared_count = min(length(base), length(compensation_values))
    for index in 1:shared_count
        base[index] -= Float64(compensation_values[index])
    end
    return base
end

function _node_voltage_phasor(
    node_voltage_phasors::AbstractVector{ComplexF64},
    node::Integer,
)
    index = Int(node)
    index == 0 && return 0.0 + 0.0im
    1 <= index <= length(node_voltage_phasors) ||
        throw(ArgumentError("nonlinear initial-state node index must address steady-state phasors"))
    return node_voltage_phasors[index]
end

function _nonlinear_initial_frequency_hz(
    steady_state_initial_sample,
    from_node::Int,
    to_node::Int,
)
    if hasproperty(
        steady_state_initial_sample,
        :node_steady_state_frequencies_hz,
    )
        frequencies = Float64.(
            steady_state_initial_sample.node_steady_state_frequencies_hz
        )
        nodes = filter(!=(0), (from_node, to_node))
        isempty(nodes) &&
            throw(ArgumentError("nonlinear initial state requires a non-reference endpoint"))
        all(node -> 1 <= node <= length(frequencies), nodes) ||
            throw(ArgumentError("nonlinear initial-state node frequency is unavailable"))
        branch_frequencies = frequencies[collect(nodes)]
        frequency = first(branch_frequencies)
        all(value -> value == frequency, branch_frequencies) ||
            throw(ArgumentError("nonlinear initial-state endpoints must share one frequency"))
        return frequency
    end
    return Float64(steady_state_initial_sample.steady_state_frequency_hz)
end

function _pseudo_nonlinear_inductor_initial_state_config(
    nonlinear_current_config::Union{Nothing,NamedTuple},
    steady_state_initial_sample,
)
    nonlinear_current_config === nothing && return nothing
    steady_state_initial_sample === nothing && return nonlinear_current_config
    nonlinear_types = Int.(get(nonlinear_current_config, :nonlinear_types, Int[]))
    count = length(nonlinear_types)
    count == 0 && return nonlinear_current_config
    pseudo_indices = findall(_is_pseudo_nonlinear_inductor_type, nonlinear_types)
    isempty(pseudo_indices) && return nonlinear_current_config
    get(
        nonlinear_current_config,
        :saturated_transformer_residual_flux_initialized,
        false,
    ) && return nonlinear_current_config
    from_nodes = Int.(get(nonlinear_current_config, :nonlinear_from_nodes, Int[]))
    to_nodes = Int.(get(nonlinear_current_config, :nonlinear_to_nodes, Int[]))
    deck_from_nodes =
        Int.(get(nonlinear_current_config, :nonlinear_deck_from_nodes, Int[]))
    deck_to_nodes =
        Int.(get(nonlinear_current_config, :nonlinear_deck_to_nodes, Int[]))
    table_start_indices =
        Int.(get(nonlinear_current_config, :nonlinear_admittance_nodes, Int[]))
    table_end_indices =
        Int.(get(nonlinear_current_config, :nonlinear_table_end_indices, Int[]))
    current_segments =
        Float64.(get(nonlinear_current_config, :nonlinear_current_segments, ones(count)))
    steady_currents = Float64.(
        get(
            nonlinear_current_config,
            :nonlinear_steady_state_current_values,
            zeros(count),
        ),
    )
    gslope = Float64.(get(nonlinear_current_config, :gslope, Float64[]))
    length(from_nodes) == count && length(to_nodes) == count ||
        throw(ArgumentError("nonlinear initial-state node vectors must match nonlinear_types"))
    if any(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
        length(deck_from_nodes) == count && length(deck_to_nodes) == count ||
            throw(ArgumentError("public pseudo-nonlinear inductor endpoints must match nonlinear_types"))
    end
    length(table_start_indices) == count && length(table_end_indices) == count ||
        throw(ArgumentError("nonlinear initial-state table vectors must match nonlinear_types"))
    length(current_segments) == count ||
        throw(ArgumentError("nonlinear current segment vector must match nonlinear_types"))
    length(steady_currents) == count ||
        throw(ArgumentError("nonlinear steady-state current vector must match nonlinear_types"))
    delta2 = Float64(get(nonlinear_current_config, :delta2, 0.0))
    delta2 > 0.0 ||
        throw(ArgumentError("delta2 must be positive for nonlinear initial state"))
    phasors = steady_state_initial_sample.node_voltage_phasors

    companion_currents = zeros(Float64, count)
    signed_segments = copy(current_segments)
    stored_voltages = zeros(Float64, count)
    table_indices = zeros(Int, count)
    for index in eachindex(nonlinear_types)
        segment_magnitude = abs(current_segments[index])
        segment_magnitude == 0.0 && (segment_magnitude = 1.0)
        table_index = table_start_indices[index] + Int(segment_magnitude) - 1
        table_start_indices[index] <= table_index <= table_end_indices[index] ||
            throw(ArgumentError("nonlinear initial current segment must address gslope table"))
        table_index <= length(gslope) ||
            throw(ArgumentError("gslope must cover nonlinear initial current segment"))
        table_indices[index] = table_index
        _is_pseudo_nonlinear_inductor_type(nonlinear_types[index]) || continue
        public_type_98 =
            nonlinear_types[index] == PSEUDO_NONLINEAR_INDUCTOR_TYPE
        from_node, to_node =
            public_type_98 ?
            (deck_from_nodes[index], deck_to_nodes[index]) :
            (from_nodes[index], abs(to_nodes[index]))
        branch_phasor =
            _node_voltage_phasor(phasors, from_node) -
            _node_voltage_phasor(phasors, to_node)
        frequency_hz =
            public_type_98 ?
            _nonlinear_initial_frequency_hz(
                steady_state_initial_sample,
                from_node,
                to_node,
            ) :
            Float64(steady_state_initial_sample.steady_state_frequency_hz)
        omega = 2.0 * pi * frequency_hz
        omega > 0.0 ||
            throw(ArgumentError("steady-state frequency must be positive for nonlinear initial state"))
        # OVER13 seeds VNONL as the steady-state flux minus the trapezoidal
        # half-step voltage. OVER14 then derives ANONL = GSLOPE * VNONL /
        # DELTA2 before SUBTS1 advances the nonlinear element ahead of solve.
        stored_voltage = imag(branch_phasor) / omega - delta2 * real(branch_phasor)
        stored_voltages[index] = stored_voltage
        signed_segments[index] =
            stored_voltage < 0.0 ? -segment_magnitude : segment_magnitude
        companion_currents[index] =
            stored_voltage * Float64(gslope[table_index]) / delta2
    end
    return merge(
        nonlinear_current_config,
        (
            seed_initial_nonlinear_state = true,
            initial_companion_current_values = companion_currents,
            initial_current_segment_values = signed_segments,
            initial_characteristic_current_values = steady_currents,
            initial_stored_voltage_values = stored_voltages,
            initial_table_index_values = table_indices,
        ),
    )
end

function _piecewise_nonlinear_inductor_characteristic_state(
    current_a::Float64,
    flux_wb::Float64,
    table_start::Int,
    table_end::Int,
    currents_a::AbstractVector{<:Real},
    fluxes_wb::AbstractVector{<:Real};
    flux_tolerance_wb::Float64,
)
    1 <= table_start < table_end <= length(currents_a) || throw(ArgumentError(
        "piecewise nonlinear-inductor characteristic range is invalid",
    ))
    table_end <= length(fluxes_wb) || throw(ArgumentError(
        "piecewise nonlinear-inductor flux characteristic is incomplete",
    ))
    lower_current_a = Float64(currents_a[table_start])
    upper_current_a = Float64(currents_a[table_end])
    current_scale = max(abs(lower_current_a), abs(upper_current_a), abs(current_a), 1.0)
    current_tolerance_a = 64.0 * eps(current_scale)
    lower_current_a - current_tolerance_a <= current_a <=
        upper_current_a + current_tolerance_a || throw(ArgumentError(
        "piecewise nonlinear-inductor harmonic current lies outside its characteristic",
    ))
    bounded_current_a = clamp(current_a, lower_current_a, upper_current_a)
    segment = if bounded_current_a <= lower_current_a
        table_start
    elseif bounded_current_a >= upper_current_a
        table_end - 1
    else
        clamp(
            searchsortedlast(
                @view(currents_a[table_start:table_end]),
                bounded_current_a;
                by=Float64,
            ) + table_start - 1,
            table_start,
            table_end - 1,
        )
    end
    current_delta_a = Float64(currents_a[segment + 1]) -
        Float64(currents_a[segment])
    flux_delta_wb = Float64(fluxes_wb[segment + 1]) -
        Float64(fluxes_wb[segment])
    current_delta_a > 0.0 && flux_delta_wb > 0.0 || throw(ArgumentError(
        "piecewise nonlinear-inductor characteristic must increase strictly",
    ))
    fraction = (bounded_current_a - Float64(currents_a[segment])) /
        current_delta_a
    characteristic_flux_wb = Float64(fluxes_wb[segment]) +
        fraction * flux_delta_wb
    flux_scale = max(abs(characteristic_flux_wb), abs(flux_wb), 1.0)
    allowance_wb = max(flux_tolerance_wb, 64.0 * eps(flux_scale))
    abs(flux_wb - characteristic_flux_wb) <= allowance_wb || throw(ArgumentError(
        "piecewise nonlinear-inductor harmonic flux-current point is not on its characteristic",
    ))
    return (
        current_a=bounded_current_a,
        flux_wb=characteristic_flux_wb,
        segment=segment,
        flux_residual_wb=abs(flux_wb - characteristic_flux_wb),
        flux_allowance_wb=allowance_wb,
    )
end

function _saturated_transformer_initial_state_config(
    nonlinear_current_config::Union{Nothing,NamedTuple},
    steady_state_initial_sample,
)
    nonlinear_current_config === nothing && return nothing
    steady_state_initial_sample === nothing && return nonlinear_current_config
    get(
        nonlinear_current_config,
        :saturated_transformer_residual_flux_initialized,
        false,
    ) && return nonlinear_current_config
    haskey(nonlinear_current_config, :saturated_transformer_branch_assembly) ||
        return nonlinear_current_config
    nonlinear_types = Int.(nonlinear_current_config.nonlinear_types)
    transformer_indices = findall(
        ==(SATURATED_TRANSFORMER_NONLINEAR_TYPE),
        nonlinear_types,
    )
    isempty(transformer_indices) && return nonlinear_current_config
    count = length(nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = Int.(nonlinear_current_config.nonlinear_to_nodes)
    table_starts = Int.(nonlinear_current_config.nonlinear_admittance_nodes)
    table_ends = Int.(nonlinear_current_config.nonlinear_table_end_indices)
    declared_currents = Float64.(
        nonlinear_current_config.nonlinear_steady_state_current_values,
    )
    declared_fluxes = Float64.(
        nonlinear_current_config.nonlinear_steady_state_flux_values,
    )
    characteristic_currents = Float64.(
        nonlinear_current_config.nonlinear_characteristic_current_values,
    )
    characteristic_fluxes = Float64.(
        nonlinear_current_config.nonlinear_characteristic_flux_values,
    )
    cchar = Float64.(nonlinear_current_config.cchar)
    gslope = Float64.(nonlinear_current_config.gslope)
    delta2 = Float64(nonlinear_current_config.delta2)
    length(from_nodes) == length(to_nodes) == length(table_starts) ==
        length(table_ends) == length(declared_currents) ==
        length(declared_fluxes) == count || throw(ArgumentError(
        "saturated-transformer initialization arrays must match nonlinear_types",
    ))
    length(characteristic_currents) == length(characteristic_fluxes) ==
        length(cchar) == length(gslope) || throw(ArgumentError(
        "saturated-transformer characteristic arrays must have equal lengths",
    ))
    delta2 > 0.0 || throw(ArgumentError(
        "saturated-transformer initialization requires a positive half timestep",
    ))
    companion_currents = Float64.(get(
        nonlinear_current_config,
        :initial_companion_current_values,
        zeros(count),
    ))
    characteristic_state_currents = Float64.(get(
        nonlinear_current_config,
        :initial_characteristic_current_values,
        declared_currents,
    ))
    stored_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_stored_voltage_values,
        zeros(count),
    ))
    runtime_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_runtime_voltage_values,
        stored_fluxes,
    ))
    current_segments = Float64.(get(
        nonlinear_current_config,
        :initial_current_segment_values,
        ones(count),
    ))
    table_indices = Int.(get(
        nonlinear_current_config,
        :initial_table_index_values,
        table_starts,
    ))
    phasors = steady_state_initial_sample.node_voltage_phasors
    for index in transformer_indices
        state = _piecewise_nonlinear_inductor_characteristic_state(
            declared_currents[index],
            declared_fluxes[index],
            table_starts[index],
            table_ends[index],
            characteristic_currents,
            characteristic_fluxes;
            flux_tolerance_wb=
                64.0 * eps(max(abs(declared_fluxes[index]), 1.0)),
        )
        branch_voltage = real(
            _node_voltage_phasor(phasors, from_nodes[index]) -
            _node_voltage_phasor(phasors, to_nodes[index]),
        )
        stored_flux = state.flux_wb - delta2 * branch_voltage
        table_index = state.segment
        stored_fluxes[index] = stored_flux
        runtime_fluxes[index] = stored_flux
        characteristic_state_currents[index] = state.current_a
        companion_currents[index] =
            (stored_flux - cchar[table_index]) * gslope[table_index] / delta2
        current_segments[index] =
            stored_flux < 0.0 ?
            -(table_index - table_starts[index] + 1) :
            table_index - table_starts[index] + 1
        table_indices[index] = table_index
    end
    return merge(
        nonlinear_current_config,
        (
            seed_initial_nonlinear_state=true,
            initial_companion_current_values=companion_currents,
            initial_characteristic_current_values=characteristic_state_currents,
            initial_stored_voltage_values=stored_fluxes,
            initial_runtime_voltage_values=runtime_fluxes,
            initial_current_segment_values=current_segments,
            initial_table_index_values=table_indices,
            saturated_transformer_residual_flux_initialized=true,
        ),
    )
end

function _piecewise_nonlinear_inductor_initial_state_config(
    nonlinear_current_config::Union{Nothing,NamedTuple},
    steady_state_initial_sample,
)
    nonlinear_current_config === nothing && return nothing
    steady_state_initial_sample === nothing && return nonlinear_current_config
    nonlinear_types = Int.(get(nonlinear_current_config, :nonlinear_types, Int[]))
    indices = findall(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
    isempty(indices) && return nonlinear_current_config
    count = length(nonlinear_types)
    deck_from_nodes = Int.(get(
        nonlinear_current_config,
        :nonlinear_deck_from_nodes,
        Int[],
    ))
    deck_to_nodes = Int.(get(
        nonlinear_current_config,
        :nonlinear_deck_to_nodes,
        Int[],
    ))
    table_starts = Int.(get(
        nonlinear_current_config,
        :nonlinear_admittance_nodes,
        Int[],
    ))
    table_ends = Int.(get(
        nonlinear_current_config,
        :nonlinear_table_end_indices,
        Int[],
    ))
    declared_currents = Float64.(get(
        nonlinear_current_config,
        :nonlinear_steady_state_current_values,
        zeros(Float64, count),
    ))
    declared_fluxes = Float64.(get(
        nonlinear_current_config,
        :nonlinear_steady_state_flux_values,
        zeros(Float64, count),
    ))
    owner_names = Symbol.(get(
        nonlinear_current_config,
        :nonlinear_owner_names,
        fill(Symbol(""), count),
    ))
    owner_line_numbers = Int.(get(
        nonlinear_current_config,
        :nonlinear_owner_line_numbers,
        zeros(Int, count),
    ))
    all(length(values) == count for values in (
        deck_from_nodes,
        deck_to_nodes,
        table_starts,
        table_ends,
        declared_currents,
        declared_fluxes,
        owner_names,
        owner_line_numbers,
    )) || throw(ArgumentError(
        "piecewise nonlinear-inductor initialization arrays must match nonlinear_types",
    ))
    currents_a = Float64.(nonlinear_current_config.cchar)
    fluxes_wb = Float64.(nonlinear_current_config.vchar)
    delta2 = Float64(nonlinear_current_config.delta2)
    isfinite(delta2) && delta2 > 0.0 || throw(ArgumentError(
        "piecewise nonlinear-inductor initialization half timestep must be finite and positive",
    ))
    flux_tolerance_wb = Float64(get(nonlinear_current_config, :flzero, 0.0))
    phasors = ComplexF64.(steady_state_initial_sample.node_voltage_phasors)
    initial_companion_currents = Float64.(get(
        nonlinear_current_config,
        :initial_companion_current_values,
        zeros(Float64, count),
    ))
    initial_characteristic_values = Float64.(get(
        nonlinear_current_config,
        :initial_characteristic_current_values,
        zeros(Float64, count),
    ))
    initial_stored_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_stored_voltage_values,
        zeros(Float64, count),
    ))
    initial_runtime_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_runtime_voltage_values,
        initial_stored_fluxes,
    ))
    initial_currents = Float64.(get(
        nonlinear_current_config,
        :initial_current_segment_values,
        zeros(Float64, count),
    ))
    initial_segments = Int.(get(
        nonlinear_current_config,
        :initial_table_index_values,
        table_starts,
    ))
    all(length(values) == count for values in (
        initial_companion_currents,
        initial_characteristic_values,
        initial_stored_fluxes,
        initial_runtime_fluxes,
        initial_currents,
        initial_segments,
    )) || throw(ArgumentError(
        "piecewise nonlinear-inductor runtime seeds must match nonlinear_types",
    ))
    initialized_fluxes = fill(NaN, count)
    initialized_predictor_fluxes = fill(NaN, count)
    initialized_currents = fill(NaN, count)
    initialized_segments = zeros(Int, count)
    characteristic_flux_residuals = fill(NaN, count)
    for index in indices
        from_node = deck_from_nodes[index]
        to_node = deck_to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(phasors, from_node) -
            _node_voltage_phasor(phasors, to_node)
        frequency_hz = _nonlinear_initial_frequency_hz(
            steady_state_initial_sample,
            from_node,
            to_node,
        )
        physical_angular_frequency = 2.0 * pi * frequency_hz
        physical_angular_frequency > 0.0 || throw(_EMTInitializationRefusal(
            :unsupported_dc_state,
            :piecewise_nonlinear_inductor_state,
            :flux,
            "piecewise nonlinear-inductor initialization requires a positive harmonic frequency",
            (owner=owner_names[index], line_no=owner_line_numbers[index]),
        ))
        reactive_angular_frequency = if get(
            steady_state_initial_sample,
            :exact_discrete_histories,
            false,
        )
            _steady_state_reactive_angular_frequency(
                steady_state_initial_sample,
                frequency_hz,
            )
        else
            physical_angular_frequency
        end
        declared_current_a = declared_currents[index]
        declared_flux_wb = declared_fluxes[index]
        harmonic_current_phasor = if declared_current_a == 0.0 && declared_flux_wb == 0.0
            0.0 + 0.0im
        else
            declared_current_a != 0.0 && declared_flux_wb != 0.0 &&
                declared_flux_wb / declared_current_a > 0.0 ||
                throw(_EMTInitializationRefusal(
                    :invalid_declared_state,
                    :piecewise_nonlinear_inductor_state,
                    :steady_state_flux_current,
                    "a nonzero piecewise nonlinear-inductor steady state requires current and flux with a positive secant inductance",
                    (
                        owner=owner_names[index],
                        line_no=owner_line_numbers[index],
                        steady_state_current_a=declared_current_a,
                        steady_state_flux_wb=declared_flux_wb,
                    ),
                ))
            secant_inductance_h = declared_flux_wb / declared_current_a
            branch_phasor /
                complex(0.0, reactive_angular_frequency * secant_inductance_h)
        end
        current_a = real(harmonic_current_phasor)
        flux_wb = imag(branch_phasor) / reactive_angular_frequency
        state = try
            _piecewise_nonlinear_inductor_characteristic_state(
                current_a,
                flux_wb,
                table_starts[index],
                table_ends[index],
                currents_a,
                fluxes_wb;
                flux_tolerance_wb,
            )
        catch error
            error isa _EMTInitializationRefusal && rethrow()
            throw(_EMTInitializationRefusal(
                :infeasible_model_state,
                :piecewise_nonlinear_inductor_state,
                :flux_current_point,
                sprint(showerror, error),
                (
                    owner=owner_names[index],
                    line_no=owner_line_numbers[index],
                    branch_current_a=current_a,
                    harmonic_flux_wb=flux_wb,
                ),
            ))
        end
        branch_voltage_v = real(branch_phasor)
        predictor_flux_wb = state.flux_wb + delta2 * branch_voltage_v
        initial_characteristic_values[index] = state.flux_wb
        initial_stored_fluxes[index] = state.flux_wb
        initial_runtime_fluxes[index] = predictor_flux_wb
        initial_currents[index] = state.current_a
        initial_segments[index] = state.segment
        initialized_fluxes[index] = state.flux_wb
        initialized_predictor_fluxes[index] = predictor_flux_wb
        initialized_currents[index] = state.current_a
        initialized_segments[index] = state.segment
        characteristic_flux_residuals[index] = state.flux_residual_wb
    end
    return merge(
        nonlinear_current_config,
        (
            seed_initial_nonlinear_state=true,
            initial_companion_current_values=initial_companion_currents,
            initial_characteristic_current_values=initial_characteristic_values,
            initial_stored_voltage_values=initial_stored_fluxes,
            initial_runtime_voltage_values=initial_runtime_fluxes,
            initial_current_segment_values=initial_currents,
            initial_table_index_values=initial_segments,
            piecewise_nonlinear_inductor_initial_flux_values=initialized_fluxes,
            piecewise_nonlinear_inductor_initial_predictor_flux_values=
                initialized_predictor_fluxes,
            piecewise_nonlinear_inductor_initial_current_values=initialized_currents,
            piecewise_nonlinear_inductor_initial_segment_values=initialized_segments,
            piecewise_nonlinear_inductor_characteristic_flux_residual_values=
                characteristic_flux_residuals,
        ),
    )
end

function _hysteretic_inductor_initial_state_config(
    nonlinear_current_config::Union{Nothing,NamedTuple},
    steady_state_initial_sample,
)
    nonlinear_current_config === nothing && return nothing
    steady_state_initial_sample === nothing && return nonlinear_current_config
    nonlinear_types = Int.(get(nonlinear_current_config, :nonlinear_types, Int[]))
    hysteretic_indices = findall(==(HYSTERETIC_INDUCTOR_NONLINEAR_TYPE), nonlinear_types)
    isempty(hysteretic_indices) && return nonlinear_current_config
    count = length(nonlinear_types)
    deck_from_nodes = Int.(get(
        nonlinear_current_config,
        :nonlinear_deck_from_nodes,
        Int[],
    ))
    deck_to_nodes = Int.(get(
        nonlinear_current_config,
        :nonlinear_deck_to_nodes,
        Int[],
    ))
    state_start_indices = Int.(get(
        nonlinear_current_config,
        :nonlinear_admittance_nodes,
        Int[],
    ))
    major_loop_start_indices = Int.(get(
        nonlinear_current_config,
        :initial_table_index_values,
        Int[],
    ))
    initial_companion_currents = Float64.(get(
        nonlinear_current_config,
        :initial_companion_current_values,
        zeros(Float64, count),
    ))
    initial_characteristic_values = Float64.(get(
        nonlinear_current_config,
        :initial_characteristic_current_values,
        zeros(Float64, count),
    ))
    initial_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_stored_voltage_values,
        zeros(Float64, count),
    ))
    runtime_fluxes = Float64.(get(
        nonlinear_current_config,
        :initial_runtime_voltage_values,
        initial_fluxes,
    ))
    initial_currents = Float64.(get(
        nonlinear_current_config,
        :initial_current_segment_values,
        zeros(Float64, count),
    ))
    residual_fluxes = Float64.(get(
        nonlinear_current_config,
        :fortran_gap_status_values,
        zeros(Float64, count),
    ))
    owner_names = Symbol.(get(
        nonlinear_current_config,
        :nonlinear_owner_names,
        fill(Symbol(""), count),
    ))
    owner_line_numbers = Int.(get(
        nonlinear_current_config,
        :nonlinear_owner_line_numbers,
        zeros(Int, count),
    ))
    all(length(values) == count for values in (
        deck_from_nodes,
        deck_to_nodes,
        state_start_indices,
        major_loop_start_indices,
        initial_companion_currents,
        initial_characteristic_values,
        initial_fluxes,
        runtime_fluxes,
        initial_currents,
        residual_fluxes,
        owner_names,
        owner_line_numbers,
    )) || throw(ArgumentError(
        "hysteretic-inductor initialization arrays must match nonlinear_types",
    ))
    cchar = Float64.(nonlinear_current_config.cchar)
    vchar = Float64.(nonlinear_current_config.vchar)
    gslope = Float64.(nonlinear_current_config.gslope)
    phasors = ComplexF64.(steady_state_initial_sample.node_voltage_phasors)
    delta2 = Float64(nonlinear_current_config.delta2)
    flux_tolerance = Float64(get(nonlinear_current_config, :flzero, 0.0))
    initialized_fluxes = fill(NaN, count)
    initialized_currents = fill(NaN, count)
    initialized_directions = zeros(Int, count)
    initialized_trace_indices = zeros(Int, count)
    lower_flux_bounds = fill(NaN, count)
    upper_flux_bounds = fill(NaN, count)
    for index in hysteretic_indices
        from_node = deck_from_nodes[index]
        to_node = deck_to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(phasors, from_node) -
            _node_voltage_phasor(phasors, to_node)
        physical_frequency_hz = _nonlinear_initial_frequency_hz(
            steady_state_initial_sample,
            from_node,
            to_node,
        )
        physical_angular_frequency = 2.0 * pi * physical_frequency_hz
        physical_angular_frequency > 0.0 || throw(_EMTInitializationRefusal(
            :unsupported_dc_state,
            :hysteretic_magnetic_state,
            :flux,
            "hysteretic-inductor initialization requires a positive harmonic frequency",
            (owner=owner_names[index], line_no=owner_line_numbers[index]),
        ))
        declared_current_a = initial_companion_currents[index]
        declared_flux_wb = initial_fluxes[index]
        reactive_angular_frequency = if get(
            steady_state_initial_sample,
            :exact_discrete_histories,
            false,
        )
            timestep_s = Float64(steady_state_initial_sample.timestep_s)
            angle = 0.5 * physical_angular_frequency * timestep_s
            abs(angle) < 0.5 * pi || throw(_EMTInitializationRefusal(
                :frequency_above_nyquist,
                :hysteretic_magnetic_state,
                :frequency_hz,
                "hysteretic-inductor initialization frequency must remain below Nyquist",
                (owner=owner_names[index], line_no=owner_line_numbers[index]),
            ))
            (2.0 / timestep_s) * tan(angle)
        else
            physical_angular_frequency
        end
        branch_current_phasor = if declared_current_a == 0.0
            0.0 + 0.0im
        else
            declared_current_a > 0.0 && declared_flux_wb > 0.0 ||
                throw(_EMTInitializationRefusal(
                    :invalid_declared_state,
                    :hysteretic_magnetic_state,
                    :steady_state_flux_current,
                    "a nonzero hysteretic-inductor current requires positive declared current and flux magnitudes",
                    (
                        owner=owner_names[index],
                        line_no=owner_line_numbers[index],
                        steady_state_current_a=declared_current_a,
                        steady_state_flux_wb=declared_flux_wb,
                    ),
                ))
            incremental_inductance_h = declared_flux_wb / declared_current_a
            branch_phasor /
                complex(0.0, reactive_angular_frequency * incremental_inductance_h)
        end
        harmonic_flux_wb = imag(branch_phasor) / physical_angular_frequency
        state = try
            _hysteretic_inductor_steady_state_table(
                cchar,
                vchar,
                gslope;
                state_start_index=state_start_indices[index],
                major_loop_start_index=major_loop_start_indices[index],
                branch_current_a=real(branch_current_phasor),
                branch_flux_wb=harmonic_flux_wb,
                branch_voltage_v=real(branch_phasor),
                residual_flux_wb=residual_fluxes[index],
                half_timestep_s=delta2,
                flux_tolerance_wb=flux_tolerance,
            )
        catch error
            error isa _EMTInitializationRefusal && rethrow()
            throw(_EMTInitializationRefusal(
                :infeasible_model_state,
                :hysteretic_magnetic_state,
                :flux_current_point,
                sprint(showerror, error),
                (
                    owner=owner_names[index],
                    line_no=owner_line_numbers[index],
                    branch_current_a=real(branch_current_phasor),
                    harmonic_flux_wb=harmonic_flux_wb,
                    residual_flux_wb=residual_fluxes[index],
                ),
            ))
        end
        cchar = state.cchar
        vchar = state.vchar
        gslope = state.gslope
        initial_companion_currents[index] = state.companion_current_a
        initial_characteristic_values[index] = state.flux_wb
        initial_fluxes[index] = state.flux_wb
        runtime_fluxes[index] = state.runtime_flux_wb
        initial_currents[index] = state.current_a
        initialized_fluxes[index] = state.flux_wb
        initialized_currents[index] = state.current_a
        initialized_directions[index] = state.direction
        initialized_trace_indices[index] = state.trace_index
        lower_flux_bounds[index] = state.lower_major_loop_flux_wb
        upper_flux_bounds[index] = state.upper_major_loop_flux_wb
    end
    return merge(
        nonlinear_current_config,
        (
            seed_initial_nonlinear_state=true,
            cchar=cchar,
            vchar=vchar,
            gslope=gslope,
            initial_companion_current_values=initial_companion_currents,
            initial_characteristic_current_values=initial_characteristic_values,
            initial_stored_voltage_values=initial_fluxes,
            initial_runtime_voltage_values=runtime_fluxes,
            initial_current_segment_values=initial_currents,
            hysteretic_initial_flux_values=initialized_fluxes,
            hysteretic_initial_current_values=initialized_currents,
            hysteretic_initial_direction_values=initialized_directions,
            hysteretic_initial_trace_indices=initialized_trace_indices,
            hysteretic_major_loop_lower_flux_values=lower_flux_bounds,
            hysteretic_major_loop_upper_flux_values=upper_flux_bounds,
        ),
    )
end

function _pseudo_nonlinear_inductor_initial_companion_current_injections(
    nonlinear_current_config::NamedTuple,
    node_count::Int,
)
    companion_currents = Float64.(
        get(nonlinear_current_config, :initial_companion_current_values, Float64[]),
    )
    isempty(companion_currents) && return zeros(Float64, node_count)
    nonlinear_types = Int.(get(nonlinear_current_config, :nonlinear_types, Int[]))
    from_nodes = Int.(get(nonlinear_current_config, :nonlinear_from_nodes, Int[]))
    to_nodes = Int.(get(nonlinear_current_config, :nonlinear_to_nodes, Int[]))
    deck_from_nodes =
        Int.(get(nonlinear_current_config, :nonlinear_deck_from_nodes, Int[]))
    deck_to_nodes =
        Int.(get(nonlinear_current_config, :nonlinear_deck_to_nodes, Int[]))
    count = length(nonlinear_types)
    length(companion_currents) == count &&
        length(from_nodes) == count &&
        length(to_nodes) == count ||
        throw(ArgumentError("initial nonlinear companion current vectors must match nonlinear_types"))
    if any(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
        length(deck_from_nodes) == count && length(deck_to_nodes) == count ||
            throw(ArgumentError("public pseudo-nonlinear inductor endpoints must match nonlinear_types"))
    end
    injections = zeros(Float64, node_count)
    for index in eachindex(nonlinear_types)
        _is_pseudo_nonlinear_inductor_type(nonlinear_types[index]) || continue
        from_node, to_node =
            nonlinear_types[index] == PSEUDO_NONLINEAR_INDUCTOR_TYPE ?
            (deck_from_nodes[index], deck_to_nodes[index]) :
            (from_nodes[index], abs(to_nodes[index]))
        1 <= from_node <= node_count ||
            throw(ArgumentError("initial nonlinear from-node must address RHS"))
        injections[from_node] -= companion_currents[index]
        if to_node != 0
            1 <= to_node <= node_count ||
                throw(ArgumentError("initial nonlinear to-node must address RHS"))
            injections[to_node] += companion_currents[index]
        end
    end
    return injections
end

function _add_current_injection_delta!(
    current_injections::Vector{Float64},
    rhs_before::AbstractVector{<:Real},
    rhs_after::AbstractVector{<:Real},
)
    length(rhs_before) >= length(current_injections) ||
        throw(ArgumentError("line RHS before-vector must cover current injections"))
    length(rhs_after) >= length(current_injections) ||
        throw(ArgumentError("line RHS after-vector must cover current injections"))
    for index in eachindex(current_injections)
        current_injections[index] += Float64(rhs_after[index]) - Float64(rhs_before[index])
    end
    return current_injections
end

function _electromagnetic_history_rhs_values(
    context::EMTStepContext,
    element_count::Union{Nothing,Integer}=nothing;
    history_voltage::Union{Nothing,AbstractVector{<:Real}}=nothing,
    advance_breqiv_history_currents::Bool=false,
    consume_breqiv_history_currents::Bool=false,
)
    rhs = context.electromagnetic_history_rhs
    fill!(rhs, 0.0)
    if history_voltage !== nothing
        length(history_voltage) >= context.system.node_count ||
            throw(ArgumentError("history voltage vector must cover every node"))
    end
    last_element =
        element_count === nothing ?
        length(context.system.elements) :
        min(Int(element_count), length(context.system.elements))
    plan = context.electromagnetic_history_plan
    for plan_index in eachindex(plan.kinds)
        plan.element_indices[plan_index] <= last_element || break
        kind = plan.kinds[plan_index]
        batch_index = plan.batch_indices[plan_index]
        if kind == SERIES_RL_HISTORY
            element = plan.series_rl_branches[batch_index]
            _, history_current = companion(element, context.dt_s)
            stamp_history_current!(rhs, element.a, element.b, history_current)
        elseif kind == SERIES_RLC_HISTORY
            element = plan.series_rlc_branches[batch_index]
            _, history_current = companion(element, context.dt_s)
            stamp_history_current!(rhs, element.a, element.b, history_current)
        elseif kind == CAPACITOR_HISTORY
            element = plan.capacitor_branches[batch_index]
            _, history_current = companion(element, context.dt_s)
            stamp_history_current!(rhs, element.a, element.b, history_current)
        elseif kind == COUPLED_INDUCTIVE_HISTORY
            element = plan.coupled_inductive_branches[batch_index]
            stamp_history_current!(rhs, element, context.dt_s)
        elseif kind == COUPLED_SERIES_RL_HISTORY
            element = plan.coupled_series_rl_branches[batch_index]
            stamp_history_current!(rhs, element, context.dt_s)
        else
            element = plan.breqiv_injections[batch_index]
            initialize_breqiv_history_injection!(element, context.dt_s)
            for phase in eachindex(element.phase_current)
                stamp_history_current!(
                    rhs,
                    element.a[phase],
                    element.b[phase],
                    element.history_current_scale * element.phase_current[phase],
                )
            end
            if advance_breqiv_history_currents
                history_voltage === nothing &&
                    throw(ArgumentError("advancing BREQIV history currents requires history_voltage"))
                advance_breqiv_history_current!(
                    element,
                    history_voltage,
                    context.dt_s;
                    consumed_for_step = consume_breqiv_history_currents,
                )
            else
                if consume_breqiv_history_currents
                    element.history_current_consumed_for_step = true
                end
            end
        end
    end
    return rhs
end

function _seed_current_source_values!(
    current_source_values::Vector{Float64},
    seed_values::AbstractVector{<:Real},
)
    length(current_source_values) >= length(seed_values) ||
        resize!(current_source_values, length(seed_values))
    fill!(current_source_values, 0.0)
    for index in eachindex(seed_values)
        current_source_values[index] = Float64(seed_values[index])
    end
    return current_source_values
end

function _clear_current_source_nodes!(
    current_source_values::Vector{Float64},
    nodes::AbstractVector{<:Integer},
)
    cleared_count = 0
    for node in nodes
        index = Int(node)
        1 <= index <= length(current_source_values) || continue
        current_source_values[index] = 0.0
        cleared_count += 1
    end
    return cleared_count
end


function _without_presolve_trace_config(config::NamedTuple)
    return Base.structdiff(config, (record_presolve_voltage_state = nothing,))
end

function _without_sparse_node_group_config(config::NamedTuple)
    return Base.structdiff(config, (sparse_node_group_config = nothing,))
end

function _without_dense_primary_nonlinear_config(config::NamedTuple)
    return Base.structdiff(config, (dense_primary_nonlinear_compensation = nothing,))
end

function _without_nonlinear_current_config(config::NamedTuple)
    return Base.structdiff(config, (nonlinear_current_config = nothing,))
end

function _has_pseudo_nonlinear_inductor_current(
    nonlinear_current_config::NamedTuple,
)
    return any(
        _is_pseudo_nonlinear_inductor_type,
        get(nonlinear_current_config, :nonlinear_types, Int[]),
    )
end

function _has_live_saturated_transformer_nonlinear_slope_branch(
    context::EMTStepContext,
)
    return !isempty(context.saturated_transformer_nonlinear_slope_branch_batch)
end

function _nonlinear_reference_node_mapping(
    nonlinear_current_config::NamedTuple,
    node_count::Int,
)
    nodal_to_reference = Int.(get(
        nonlinear_current_config,
        :deck_to_runtime_node_indices,
        Int[],
    ))
    length(nodal_to_reference) == node_count ||
        throw(ArgumentError("primary nonlinear node mapping must cover every nodal unknown"))
    reference_count = Int(get(
        nonlinear_current_config,
        :nonlinear_required_node_count,
        node_count + 1,
    ))
    reference_to_nodal = zeros(Int, reference_count)
    for nodal_index in 1:node_count
        reference_index = nodal_to_reference[nodal_index]
        2 <= reference_index <= reference_count ||
            throw(ArgumentError("primary nonlinear node mapping must reserve index one for ground"))
        reference_to_nodal[reference_index] == 0 ||
            throw(ArgumentError("primary nonlinear node mapping must be one-to-one"))
        reference_to_nodal[reference_index] = nodal_index
    end
    return nodal_to_reference, reference_to_nodal
end

function _dense_primary_nonlinear_inverse_columns!(
    over16_state::OVER16AcceptedTimestepState,
    nonlinear_current_config::NamedTuple,
    context::EMTStepContext,
)
    node_count = context.system.node_count
    nodal_to_reference, reference_to_nodal =
        _nonlinear_reference_node_mapping(nonlinear_current_config, node_count)
    reference_count = length(reference_to_nodal)
    nonlinear_types = Int.(nonlinear_current_config.nonlinear_types)
    component_count = length(nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = Int.(nonlinear_current_config.nonlinear_to_nodes)
    length(from_nodes) == component_count && length(to_nodes) == component_count ||
        throw(ArgumentError("primary nonlinear endpoint vectors must match nonlinear owners"))

    columns = zeros(Float64, reference_count * component_count)
    injection = zeros(Float64, node_count)
    solution = zeros(Float64, node_count)
    factor = similar(context.system.y)
    for component in 1:component_count
        fill!(injection, 0.0)
        from_node = from_nodes[component]
        to_node = abs(to_nodes[component])
        if from_node != 1
            nodal_index = reference_to_nodal[from_node]
            nodal_index != 0 ||
                throw(ArgumentError("primary nonlinear from-node is not a nodal unknown"))
            injection[nodal_index] -= 1.0
        end
        if to_node != 1
            nodal_index = reference_to_nodal[to_node]
            nodal_index != 0 ||
                throw(ArgumentError("primary nonlinear to-node is not a nodal unknown"))
            injection[nodal_index] += 1.0
        end
        copyto!(factor, context.system.y)
        Nodal.solve_dense!(solution, factor, injection)
        offset = (component - 1) * reference_count
        for nodal_index in 1:node_count
            columns[offset + nodal_to_reference[nodal_index]] = solution[nodal_index]
        end
    end

    inverse = over16_state.nonlinear_inverse
    resize!(inverse.znonl, length(columns))
    inverse.znonl .= columns
    inverse.ntot = reference_count
    inverse.ncomp = component_count
    inverse.update_count += 1
    return columns
end

function _dense_primary_nonlinear_injections(
    reference_values::AbstractVector{<:Real},
    nonlinear_current_config::NamedTuple,
    node_count::Int,
)
    nodal_to_reference, _ =
        _nonlinear_reference_node_mapping(nonlinear_current_config, node_count)
    injections = zeros(Float64, node_count)
    for nodal_index in 1:node_count
        reference_index = nodal_to_reference[nodal_index]
        reference_index <= length(reference_values) || continue
        injections[nodal_index] = Float64(reference_values[reference_index])
    end
    return injections
end

function _dense_primary_nonlinear_reference_rhs(
    values::AbstractVector{<:Real},
    nonlinear_current_config::NamedTuple,
    node_count::Int,
)
    nodal_to_reference, reference_to_nodal =
        _nonlinear_reference_node_mapping(nonlinear_current_config, node_count)
    reference_count = length(reference_to_nodal)
    length(values) == reference_count && return Float64.(values)
    length(values) == node_count || throw(ArgumentError(
        "primary nonlinear base RHS must use nodal or reference coordinates",
    ))
    reference_values = zeros(Float64, reference_count)
    for nodal_index in 1:node_count
        reference_values[nodal_to_reference[nodal_index]] = Float64(values[nodal_index])
    end
    return reference_values
end

function _apply_dense_primary_hysteretic_admittance_deltas!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
)
    hasproperty(nonlinear_current_result, :hysteretic_inductor_admittance_deltas) ||
        return 0
    deltas = Float64.(nonlinear_current_result.hysteretic_inductor_admittance_deltas)
    types = Int.(nonlinear_current_config.nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    length(deltas) == length(types) == length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("primary nonlinear restamp arrays must match nonlinear owners"))
    _, reference_to_nodal = _nonlinear_reference_node_mapping(
        nonlinear_current_config,
        context.system.node_count,
    )
    update_count = 0
    for index in eachindex(types)
        types[index] == -96 || continue
        delta = deltas[index]
        delta == 0.0 && continue
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        from_nodal = from_node == 1 ? 0 : reference_to_nodal[from_node]
        to_nodal = to_node == 1 ? 0 : reference_to_nodal[to_node]
        from_node == 1 || from_nodal != 0 ||
            throw(ArgumentError("hysteretic-inductor from-node is not a nodal unknown"))
        to_node == 1 || to_nodal != 0 ||
            throw(ArgumentError("hysteretic-inductor to-node is not a nodal unknown"))
        from_nodal != 0 && (context.system.y[from_nodal, from_nodal] += delta)
        to_nodal != 0 && (context.system.y[to_nodal, to_nodal] += delta)
        if from_nodal != 0 && to_nodal != 0
            context.system.y[from_nodal, to_nodal] -= delta
            context.system.y[to_nodal, from_nodal] -= delta
        end
        update_count += 1
    end
    return update_count
end

function _apply_dense_primary_switching_resistor_admittance_deltas!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
)
    hasproperty(nonlinear_current_result, :switching_resistor_admittance_deltas) ||
        return 0
    deltas = Float64.(nonlinear_current_result.switching_resistor_admittance_deltas)
    types = Int.(nonlinear_current_config.nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    length(deltas) == length(types) == length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("switching resistor restamp arrays must match nonlinear owners"))
    _, reference_to_nodal = _nonlinear_reference_node_mapping(
        nonlinear_current_config,
        context.system.node_count,
    )
    update_count = 0
    for index in eachindex(types)
        types[index] == SWITCHING_NONLINEAR_RESISTOR_TYPE || continue
        delta = deltas[index]
        delta == 0.0 && continue
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        from_nodal = from_node == 1 ? 0 : reference_to_nodal[from_node]
        to_nodal = to_node == 1 ? 0 : reference_to_nodal[to_node]
        from_node == 1 || from_nodal != 0 ||
            throw(ArgumentError("switching resistor from-node is not a nodal unknown"))
        to_node == 1 || to_nodal != 0 ||
            throw(ArgumentError("switching resistor to-node is not a nodal unknown"))
        from_nodal != 0 && (context.system.y[from_nodal, from_nodal] += delta)
        to_nodal != 0 && (context.system.y[to_nodal, to_nodal] += delta)
        if from_nodal != 0 && to_nodal != 0
            context.system.y[from_nodal, to_nodal] -= delta
            context.system.y[to_nodal, from_nodal] -= delta
        end
        update_count += 1
    end
    return update_count
end

function _apply_dense_primary_timed_resistance_admittance_deltas!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
)
    hasproperty(nonlinear_current_result, :timed_resistance_admittance_deltas) ||
        return 0
    deltas = Float64.(nonlinear_current_result.timed_resistance_admittance_deltas)
    types = Int.(nonlinear_current_config.nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    length(deltas) == length(types) == length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("timed-resistance restamp arrays must match nonlinear owners"))
    _, reference_to_nodal = _nonlinear_reference_node_mapping(
        nonlinear_current_config,
        context.system.node_count,
    )
    update_count = 0
    for index in eachindex(types)
        types[index] == TRIGGERED_TIMED_RESISTANCE_TYPE || continue
        delta = deltas[index]
        delta == 0.0 && continue
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        from_nodal = from_node == 1 ? 0 : reference_to_nodal[from_node]
        to_nodal = to_node == 1 ? 0 : reference_to_nodal[to_node]
        from_node == 1 || from_nodal != 0 ||
            throw(ArgumentError("timed-resistance from-node is not a nodal unknown"))
        to_node == 1 || to_nodal != 0 ||
            throw(ArgumentError("timed-resistance to-node is not a nodal unknown"))
        from_nodal != 0 && (context.system.y[from_nodal, from_nodal] += delta)
        to_nodal != 0 && (context.system.y[to_nodal, to_nodal] += delta)
        if from_nodal != 0 && to_nodal != 0
            context.system.y[from_nodal, to_nodal] -= delta
            context.system.y[to_nodal, from_nodal] -= delta
        end
        update_count += 1
    end
    return update_count
end

function _apply_dense_primary_pseudo_nonlinear_inductor_admittance_deltas!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
)
    hasproperty(nonlinear_current_result, :saturated_transformer_admittance_deltas) ||
        return 0
    deltas = Float64.(nonlinear_current_result.saturated_transformer_admittance_deltas)
    types = Int.(nonlinear_current_config.nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    length(deltas) == length(types) == length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("pseudo-nonlinear inductor restamp arrays must match nonlinear owners"))
    _, reference_to_nodal = _nonlinear_reference_node_mapping(
        nonlinear_current_config,
        context.system.node_count,
    )
    update_count = 0
    for index in eachindex(types)
        types[index] == PSEUDO_NONLINEAR_INDUCTOR_TYPE || continue
        delta = deltas[index]
        delta == 0.0 && continue
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        from_nodal = from_node == 1 ? 0 : reference_to_nodal[from_node]
        to_nodal = to_node == 1 ? 0 : reference_to_nodal[to_node]
        from_node == 1 || from_nodal != 0 ||
            throw(ArgumentError("pseudo-nonlinear inductor from-node is not a nodal unknown"))
        to_node == 1 || to_nodal != 0 ||
            throw(ArgumentError("pseudo-nonlinear inductor to-node is not a nodal unknown"))
        from_nodal != 0 && (context.system.y[from_nodal, from_nodal] += delta)
        to_nodal != 0 && (context.system.y[to_nodal, to_nodal] += delta)
        if from_nodal != 0 && to_nodal != 0
            context.system.y[from_nodal, to_nodal] -= delta
            context.system.y[to_nodal, from_nodal] -= delta
        end
        update_count += 1
    end
    return update_count
end

function _apply_dense_primary_hysteretic_companion_admittance!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    over16_state::OVER16AcceptedTimestepState,
)
    types = Int.(nonlinear_current_config.nonlinear_types)
    state_starts = Int.(nonlinear_current_config.nonlinear_admittance_nodes)
    state_gslope = over16_state.nonlinear_inverse.gslope
    gslope = isempty(state_gslope) ?
        Float64.(nonlinear_current_config.gslope) :
        state_gslope
    current_segments = isempty(over16_state.nonlinear_inverse.curr) ?
        Float64.(get(
            nonlinear_current_config,
            :initial_current_segment_values,
            zeros(Float64, length(types)),
        )) :
        over16_state.nonlinear_inverse.curr
    synthetic_result = (
        hysteretic_inductor_admittance_deltas = Float64[
            types[index] == -96 ? gslope[state_starts[index] + 1] : 0.0
            for index in eachindex(types)
        ],
    )
    update_count = _apply_dense_primary_hysteretic_admittance_deltas!(
        context,
        nonlinear_current_config,
        synthetic_result,
    )
    timed_initial_result = (
        timed_resistance_admittance_deltas = Float64[
            types[index] == TRIGGERED_TIMED_RESISTANCE_TYPE &&
            round(Int, current_segments[index]) > 0 ?
                gslope[
                    Int(nonlinear_current_config.nonlinear_admittance_nodes[index]) +
                    round(Int, current_segments[index]) - 1
                ] :
                0.0
            for index in eachindex(types)
        ],
    )
    timed_update_count =
        _apply_dense_primary_timed_resistance_admittance_deltas!(
            context,
            nonlinear_current_config,
            timed_initial_result,
        )
    switching_initial_result = (
        switching_resistor_admittance_deltas = Float64[
            types[index] == SWITCHING_NONLINEAR_RESISTOR_TYPE &&
            round(Int, current_segments[index]) != 0 ?
                gslope[
                    Int(nonlinear_current_config.nonlinear_admittance_nodes[index]) +
                    abs(round(Int, current_segments[index])) - 1
                ] :
                0.0
            for index in eachindex(types)
        ],
    )
    switching_update_count =
        _apply_dense_primary_switching_resistor_admittance_deltas!(
            context,
            nonlinear_current_config,
            switching_initial_result,
        )
    if update_count + timed_update_count + switching_update_count > 0
        copyto!(context.system.y_factor, context.system.y)
        Nodal.solve_dense!(
            context.system.v,
            context.system.y_factor,
            context.system.rhs,
        )
    end
    return update_count + timed_update_count + switching_update_count
end

function _dense_primary_hysteretic_companion_rhs!(
    reference_rhs::Vector{Float64},
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
)
    types = Int.(nonlinear_current_config.nonlinear_types)
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    source_deltas = Float64.(
        nonlinear_current_result.hysteretic_inductor_source_current_deltas,
    )
    companion_currents = Float64.(nonlinear_current_result.anonl)
    length(types) == length(from_nodes) == length(to_nodes) ==
        length(source_deltas) == length(companion_currents) || throw(ArgumentError(
        "primary hysteretic companion arrays must match nonlinear owners",
    ))
    for index in eachindex(types)
        types[index] == HYSTERETIC_INDUCTOR_NONLINEAR_TYPE || continue
        correction = companion_currents[index] - source_deltas[index]
        reference_rhs[from_nodes[index]] -= correction
        reference_rhs[to_nodes[index]] += correction
    end
    return reference_rhs
end

function _dense_primary_nonlinear_steady_state_seed(
    nonlinear_current_config::NamedTuple,
    over16_state::OVER16AcceptedTimestepState,
    voltage::AbstractVector{Float64},
)
    over16_state.nonlinear_inverse.current_update_count == 0 ||
        return nonlinear_current_config
    types = Int.(nonlinear_current_config.nonlinear_types)
    any(type -> type == -96 || type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE, types) ||
        return nonlinear_current_config
    voltage_context = _deck_nonlinear_voltage_context(
        nonlinear_current_config,
        voltage,
    )
    from_nodes = Int.(nonlinear_current_config.nonlinear_from_nodes)
    to_nodes = abs.(Int.(nonlinear_current_config.nonlinear_to_nodes))
    steady_flux = Float64.(nonlinear_current_config.initial_stored_voltage_values)
    delta2 = Float64(nonlinear_current_config.delta2)
    runtime_flux = copy(steady_flux)
    initial_currents = Float64.(get(
        nonlinear_current_config,
        :initial_current_segment_values,
        zeros(Float64, length(types)),
    ))
    length(initial_currents) == length(types) ||
        throw(ArgumentError("nonlinear initial-current values must match nonlinear owners"))
    inverse = over16_state.nonlinear_inverse
    response_node_count = inverse.ntot
    length(inverse.znonl) == response_node_count * inverse.ncomp ||
        throw(ArgumentError("nonlinear inverse columns must match their recorded dimensions"))
    for index in eachindex(types)
        (types[index] == -96 || types[index] == PIECEWISE_NONLINEAR_INDUCTOR_TYPE) ||
            continue
        branch_voltage = voltage_context[from_nodes[index]] - voltage_context[to_nodes[index]]
        if types[index] == PIECEWISE_NONLINEAR_INDUCTOR_TYPE
            for response_index in eachindex(types)
                types[response_index] == PIECEWISE_NONLINEAR_INDUCTOR_TYPE || continue
                offset = (response_index - 1) * response_node_count
                branch_voltage += initial_currents[response_index] * (
                    inverse.znonl[offset + from_nodes[index]] -
                    inverse.znonl[offset + to_nodes[index]]
                )
            end
        end
        runtime_flux[index] -= delta2 * branch_voltage
    end
    return merge(
        nonlinear_current_config,
        (initial_runtime_voltage_values = runtime_flux,),
    )
end

function _apply_dense_primary_nonlinear_solution!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    nonlinear_current_result,
    nonlinear_base_rhs::AbstractVector{<:Real},
    stored_injections::AbstractVector{<:Real},
    source_voltage_constraint_result,
)
    nonlinear_current_result === nothing &&
        throw(ArgumentError("primary nonlinear timestep did not produce a current update"))
    _apply_dense_primary_pseudo_nonlinear_inductor_admittance_deltas!(
        context,
        nonlinear_current_config,
        nonlinear_current_result,
    )
    _apply_dense_primary_hysteretic_admittance_deltas!(
        context,
        nonlinear_current_config,
        nonlinear_current_result,
    )
    _apply_dense_primary_switching_resistor_admittance_deltas!(
        context,
        nonlinear_current_config,
        nonlinear_current_result,
    )
    _apply_dense_primary_timed_resistance_admittance_deltas!(
        context,
        nonlinear_current_config,
        nonlinear_current_result,
    )
    reference_rhs = Float64.(nonlinear_current_result.rhs)
    _dense_primary_hysteretic_companion_rhs!(
        reference_rhs,
        nonlinear_current_config,
        nonlinear_current_result,
    )
    compensation_rhs = copy(reference_rhs)
    shared_count = min(length(compensation_rhs), length(nonlinear_base_rhs))
    for index in 1:shared_count
        compensation_rhs[index] -= Float64(nonlinear_base_rhs[index])
    end
    desired_injections = _dense_primary_nonlinear_injections(
        compensation_rhs,
        nonlinear_current_config,
        context.system.node_count,
    )
    length(stored_injections) == context.system.node_count ||
        throw(ArgumentError("stored primary nonlinear injections must cover every node"))
    for node in 1:context.system.node_count
        context.system.rhs[node] +=
            desired_injections[node] - Float64(stored_injections[node])
    end
    if get(source_voltage_constraint_result, :applied, false)
        Nodal._apply_voltage_constraints!(
            context.system.y,
            context.system.rhs,
            source_voltage_constraint_result.nodes,
            source_voltage_constraint_result.values,
        )
    end
    copyto!(context.system.y_factor, context.system.y)
    Nodal.solve_dense!(
        context.system.v,
        context.system.y_factor,
        context.system.rhs,
    )
    accepted_voltage_values = get(
        nonlinear_current_config,
        :accepted_voltage_values,
        nothing,
    )
    if accepted_voltage_values !== nothing
        length(accepted_voltage_values) == context.system.node_count ||
            throw(ArgumentError("accepted primary nonlinear voltage state must cover every node"))
        tolerance = Float64(get(
            nonlinear_current_config,
            :accepted_voltage_tolerance,
            1.0e-12,
        ))
        tolerance >= 0.0 && isfinite(tolerance) ||
            throw(ArgumentError("accepted primary nonlinear voltage tolerance must be finite and nonnegative"))
        for node in 1:context.system.node_count
            previous = Float64(accepted_voltage_values[node])
            current = context.system.v[node]
            if isfinite(previous) &&
               abs(current - previous) <= tolerance * max(1.0, abs(previous))
                context.system.v[node] = previous
            end
            accepted_voltage_values[node] = context.system.v[node]
        end
    end
    return context.system.v
end

function _apply_switched_topology_admittance!(
    context::EMTStepContext,
    state::OVER16AcceptedTimestepState,
)
    node_count = context.system.node_count
    admittance = state.switch_admittance.admittance
    size(admittance) == (node_count + 1, node_count + 1) ||
        throw(ArgumentError(
            "switched nonlinear admittance must include one reference node",
        ))
    @views context.system.y .= admittance[2:end, 2:end]
    return context.system.y
end

function _sync_switched_nonlinear_network_solution!(
    context::EMTStepContext,
    state::OVER16AcceptedTimestepState,
    switch_current_config,
)
    node_count = context.system.node_count

    solution = state.switch_current.network_solution
    resize!(solution, node_count + 1)
    solution[1] = 0.0
    @views solution[2:end] .= context.system.v

    if switch_current_config !== nothing
        from_nodes = Int.(switch_current_config.from_nodes)
        to_nodes = Int.(switch_current_config.to_nodes)
        conductances = state.switch_admittance.switch_conductances
        length(from_nodes) == length(to_nodes) == length(conductances) ||
            throw(ArgumentError(
                "switched nonlinear current metadata must match switch count",
            ))
        resize!(state.switch_current.switch_currents, length(from_nodes))
        resize!(state.switch_post_current.switch_currents, length(from_nodes))
        for index in eachindex(from_nodes)
            current = conductances[index] *
                      (solution[from_nodes[index]] - solution[to_nodes[index]])
            state.switch_current.switch_currents[index] = current
            state.switch_post_current.switch_currents[index] = current
        end
    end
    return context.system.v
end

function _reference_augmented_step_admittance(context::EMTStepContext)
    augmented = context.reference_augmented_admittance
    fill!(augmented, 0.0)
    @views augmented[2:end, 2:end] .= context.system.y
    return augmented
end

function _stamp_reference_conductance_delta!(
    admittance::AbstractMatrix{Float64},
    from_node::Integer,
    to_node::Integer,
    conductance::Real,
)
    a = Int(from_node)
    b = Int(to_node)
    g = Float64(conductance)
    1 <= a <= size(admittance, 1) ||
        throw(ArgumentError("switch from-node is outside augmented admittance"))
    1 <= b <= size(admittance, 1) ||
        throw(ArgumentError("switch to-node is outside augmented admittance"))
    a != b || throw(ArgumentError("switch endpoints must be distinct"))
    a != 1 && (admittance[a, a] += g)
    b != 1 && (admittance[b, b] += g)
    if a != 1 && b != 1
        admittance[a, b] -= g
        admittance[b, a] -= g
    end
    return admittance
end

function _switch_closed_mask_after_modswt(
    closed_mask::AbstractVector{Bool},
    modswt::AbstractVector{Int},
)
    updated = collect(closed_mask)
    for entry in modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        row = abs(entry)
        1 <= row <= length(updated) ||
            throw(ArgumentError("MODSWT row index is outside switch table"))
        updated[row] = entry > 0
    end
    return updated
end

function _sync_step_context_switch_base_admittance!(
    over16_state::OVER16AcceptedTimestepState,
    context::EMTStepContext,
    switch_admittance_config::NamedTuple,
)
    base = _reference_augmented_step_admittance(context)
    from_nodes = Int.(switch_admittance_config.from_nodes)
    to_nodes = Int.(switch_admittance_config.to_nodes)
    switch_count = length(from_nodes)
    length(to_nodes) == switch_count ||
        throw(ArgumentError("switch endpoint vector lengths must match"))
    open_conductances = Float64.(get(
        switch_admittance_config,
        :open_conductances,
        fill(Float64(get(switch_admittance_config, :open_conductance, 0.0)), switch_count),
    ))
    length(open_conductances) == switch_count ||
        throw(ArgumentError("open_conductances length must match switch endpoints"))
    requested_closed_mask = _switch_closed_mask_after_modswt(
        over16_state.switch_topology.closed_mask,
        Int.(get(switch_admittance_config, :modswt, Int[])),
    )
    length(requested_closed_mask) == switch_count ||
        throw(ArgumentError("switch closed-mask length must match switch endpoints"))
    for row in 1:switch_count
        if !requested_closed_mask[row] && open_conductances[row] != 0.0
            _stamp_reference_conductance_delta!(
                base,
                from_nodes[row],
                to_nodes[row],
                -open_conductances[row],
            )
        end
    end
    size(over16_state.switch_admittance.base_admittance) == size(base) ||
        throw(ArgumentError("switch admittance state size must match context admittance"))
    over16_state.switch_admittance.base_admittance .= base
    augmented_rhs = context.reference_augmented_rhs
    fill!(augmented_rhs, 0.0)
    @views augmented_rhs[2:end] .= context.system.rhs
    length(over16_state.switch_current.rhs) == length(augmented_rhs) ||
        resize!(over16_state.switch_current.rhs, length(augmented_rhs))
    over16_state.switch_current.rhs .= augmented_rhs
    return base
end

function _sync_saturated_transformer_nonlinear_slope_branches!(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
    over16_state::OVER16AcceptedTimestepState,
)
    nonlinear_types = get(nonlinear_current_config, :nonlinear_types, Int[])
    count = length(nonlinear_types)
    count == 0 && return (matched_count = 0, mutation_count = 0)
    from_nodes = get(nonlinear_current_config, :nonlinear_from_nodes, Int[])
    to_nodes = get(nonlinear_current_config, :nonlinear_to_nodes, Int[])
    table_start_indices =
        get(nonlinear_current_config, :nonlinear_admittance_nodes, Int[])
    table_end_indices =
        get(nonlinear_current_config, :nonlinear_table_end_indices, Int[])
    all(length(values) == count for values in (
        from_nodes,
        to_nodes,
        table_start_indices,
        table_end_indices,
    )) || throw(ArgumentError("nonlinear slope branch config vectors must match nonlinear_types"))
    current_segments = over16_state.nonlinear_inverse.curr
    length(current_segments) >= count ||
        throw(ArgumentError("nonlinear inverse current segments must cover nonlinear_types"))
    state_slopes = over16_state.nonlinear_inverse.gslope
    slopes = isempty(state_slopes) ?
        get(nonlinear_current_config, :gslope, Float64[]) :
        state_slopes

    matched_count = 0
    mutation_count = 0
    for index in eachindex(nonlinear_types)
        _is_pseudo_nonlinear_inductor_type(nonlinear_types[index]) || continue
        segment = Int(abs(current_segments[index]))
        segment > 0 ||
            throw(ArgumentError("pseudo-nonlinear inductor current segment must be nonzero"))
        table_start_index = Int(table_start_indices[index])
        table_end_index = Int(table_end_indices[index])
        table_index = table_start_index + segment - 1
        table_start_index <= table_index <= table_end_index ||
            throw(ArgumentError("pseudo-nonlinear inductor current segment must address gslope table"))
        table_index <= length(slopes) ||
            throw(ArgumentError("gslope must cover pseudo-nonlinear inductor slope branches"))
        if Int(nonlinear_types[index]) == PSEUDO_NONLINEAR_INDUCTOR_TYPE
            deck_from_nodes = get(
                nonlinear_current_config,
                :nonlinear_deck_from_nodes,
                Int[],
            )
            deck_to_nodes = get(
                nonlinear_current_config,
                :nonlinear_deck_to_nodes,
                Int[],
            )
            length(deck_from_nodes) == count == length(deck_to_nodes) ||
                throw(ArgumentError("public pseudo-nonlinear inductor endpoints must cover nonlinear owners"))
            expected_from_node = Int(deck_from_nodes[index])
            expected_to_node = Int(deck_to_nodes[index])
        else
            expected_from_node = Int(from_nodes[index])
            expected_to_node = abs(Int(to_nodes[index]))
        end
        matching_branch = nothing
        match_count = 0
        for element in context.saturated_transformer_nonlinear_slope_branch_batch
            if element.from_node == expected_from_node &&
               element.to_node == expected_to_node
                matching_branch = element
                match_count += 1
            end
        end
        match_count == 1 ||
            throw(ArgumentError("pseudo-nonlinear inductor slope branch must have one live nodal owner"))
        matched_count += 1
        mutation_count += set_saturated_transformer_nonlinear_slope!(
            matching_branch,
            Float64(slopes[table_index]),
        )
    end
    return (matched_count = matched_count, mutation_count = mutation_count)
end

function _execute_nonlinear_nodal_step_solver(
    nonlinear_nodal_step_solver,
    context,
    current_injections,
    source_voltage_constraint_result,
)
    applicable(
        nonlinear_nodal_step_solver,
        context,
        current_injections,
        source_voltage_constraint_result,
    ) || throw(ArgumentError(
        "nonlinear nodal step solver must be callable with the EMT context, current injections, and source-voltage constraints",
    ))
    return nonlinear_nodal_step_solver(
        context,
        current_injections,
        source_voltage_constraint_result,
    )
end

function step_with_over16_boundary!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState;
    collect_step_diagnostics::Bool = true,
    over16_kwargs...,
)
    context.step_index <= context.step_count ||
        throw(ArgumentError("fixed-step EMT context is already complete"))
    step_before = context.step_index
    t_before = context.t_s
    raw_over16_kwargs = (; over16_kwargs...)
    nonlinear_nodal_step_solver = get(
        raw_over16_kwargs,
        :nonlinear_nodal_step_solver,
        nothing,
    )
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (nonlinear_nodal_step_solver = nothing,),
    )
    hybrid_substep = get(raw_over16_kwargs, :hybrid_substep, false)
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (hybrid_substep = nothing,),
    )
    hybrid_report_step_s = Float64(get(
        raw_over16_kwargs,
        :hybrid_report_step_s,
        context.dt_s,
    ))
    hybrid_report_step_s > 0.0 && isfinite(hybrid_report_step_s) ||
        throw(ArgumentError("hybrid report step must be finite and positive"))
    raw_over16_kwargs = Base.structdiff(
        raw_over16_kwargs,
        (hybrid_report_step_s = nothing,),
    )
    record_presolve_voltage_state =
        get(raw_over16_kwargs, :record_presolve_voltage_state, false)
    presolve_voltage =
        record_presolve_voltage_state ? copy(context.system.v) : Float64[]
    solver_over16_kwargs = _without_presolve_trace_config(raw_over16_kwargs)
    sparse_node_group_config =
        get(solver_over16_kwargs, :sparse_node_group_config, nothing)
    sparse_node_group_enabled =
        sparse_node_group_config !== nothing &&
        get(sparse_node_group_config, :enabled, true)
    solver_over16_kwargs = _without_sparse_node_group_config(solver_over16_kwargs)
    dense_primary_nonlinear_compensation = get(
        solver_over16_kwargs,
        :dense_primary_nonlinear_compensation,
        false,
    )
    solver_over16_kwargs =
        _without_dense_primary_nonlinear_config(solver_over16_kwargs)
    seed_current_source_values =
        get(solver_over16_kwargs, :seed_current_source_values, false)
    current_source_seed_element_count =
        get(solver_over16_kwargs, :current_source_seed_element_count, nothing)
    current_source_seed_clear_nodes =
        get(solver_over16_kwargs, :current_source_seed_clear_nodes, Int[])
    solver_over16_kwargs = _without_current_source_seed(solver_over16_kwargs)
    source_config = get(solver_over16_kwargs, :source_config, nothing)
    distributed_transposed_line_config =
        get(solver_over16_kwargs, :distributed_transposed_line_config, nothing)
    frequency_dependent_line_config =
        get(solver_over16_kwargs, :frequency_dependent_line_config, nothing)
    nonlinear_current_config =
        get(solver_over16_kwargs, :nonlinear_current_config, nothing)
    nonlinear_current_compensation_enabled =
        nonlinear_current_config !== nothing &&
        get(nonlinear_current_config, :enabled, true)
    reset_current_source_values =
        get(solver_over16_kwargs, :reset_current_source_values, false)
    solver_over16_kwargs =
        _without_distributed_transposed_line_config(solver_over16_kwargs)
    solver_over16_kwargs =
        _without_frequency_dependent_line_config(solver_over16_kwargs)
    solver_over16_kwargs = _without_current_source_reset(solver_over16_kwargs)
    sync_switch_base_admittance =
        get(solver_over16_kwargs, :switch_base_admittance_from_step_context, false)
    source_function_control_values =
        source_config === nothing ?
        Float64[] :
        get(get(source_config, :kwargs, NamedTuple()), :xtcs_values, Float64[])
    source_function_row_result =
        context.source_function_runtime === nothing ?
        nothing :
        _synchronize_source_function_row_slots!(
            context.source_function_runtime,
            context.t_s,
            source_function_control_values,
        )
    source_voltage_constraint_result =
        _electromagnetic_source_voltage_constraints(
            over16_state,
            source_config,
            context.system.node_count,
            _source_voltage_constraint_time(source_config, context.t_s),
            context.source_function_runtime !== nothing &&
                _source_signal_stage_recording_active(
                    context.source_function_runtime,
                ),
        )
    presolve_line_voltage =
        _source_constrained_presolve_voltage(
            context.system.v,
            source_voltage_constraint_result,
            presolve_voltage,
            record_presolve_voltage_state,
        )
    electromagnetic_history_rhs =
        seed_current_source_values ?
        _electromagnetic_history_rhs_values(
            context,
            current_source_seed_element_count;
            history_voltage = presolve_line_voltage,
            advance_breqiv_history_currents = false,
            consume_breqiv_history_currents = false,
        ) :
        Float64[]
    current_injections =
        _electromagnetic_rhs_current_injections(
            over16_state,
            context.system.node_count,
            solver_over16_kwargs,
        )
    _apply_ideal_transformer_rhs_values!(
        current_injections,
        source_voltage_constraint_result,
    )
    stored_nonlinear_current_compensation_injections =
        nonlinear_current_compensation_enabled ?
        (
            dense_primary_nonlinear_compensation ?
            zeros(Float64, context.system.node_count) :
            _nonlinear_current_compensation_injection_values(
                over16_state,
                context.system.node_count,
            )
        ) :
        zeros(Float64, context.system.node_count)
    use_initial_nonlinear_companion_current =
        nonlinear_current_compensation_enabled &&
        get(nonlinear_current_config, :seed_initial_nonlinear_state, false) &&
        over16_state.nonlinear_inverse.current_update_count == 0 &&
        over16_state.nonlinear_inverse.update_count == 0 &&
        over16_state.nonlinear_inverse.source_column_update_count == 0
    nonlinear_current_compensation_injections =
        use_initial_nonlinear_companion_current ?
        _pseudo_nonlinear_inductor_initial_companion_current_injections(
            nonlinear_current_config,
            context.system.node_count,
        ) :
        stored_nonlinear_current_compensation_injections
    current_source_seed_voltage_constraint_zero_count = 0
    if reset_current_source_values
        if seed_current_source_values
            _seed_current_source_values!(
                over16_state.source.f_values,
                electromagnetic_history_rhs,
            )
            current_source_seed_voltage_constraint_zero_count =
                _clear_current_source_nodes!(
                    over16_state.source.f_values,
                    current_source_seed_clear_nodes,
                )
        else
            _ensure_float_vector_length!(
                over16_state.source.f_values,
                context.system.node_count,
            )
            fill!(over16_state.source.f_values, 0.0)
        end
    end
    distributed_transposed_line_current_result =
        _distributed_transposed_line_deck_current_injection!(
            distributed_transposed_line_config,
            over16_state.source.f_values,
            collect_diagnostics = collect_step_diagnostics,
        )
    if distributed_transposed_line_current_result !== nothing
        _add_current_injection_delta!(
            current_injections,
            distributed_transposed_line_current_result.rhs_before_values,
            distributed_transposed_line_current_result.rhs_after_values,
        )
    end
    frequency_dependent_line_current_result =
        _frequency_dependent_line_deck_current_injection!(
            frequency_dependent_line_config,
            over16_state.source.f_values,
        )
    if frequency_dependent_line_current_result !== nothing
        _add_current_injection_delta!(
            current_injections,
            frequency_dependent_line_current_result.rhs_before_values,
            frequency_dependent_line_current_result.rhs_after_values,
        )
    end
    nonlinear_current_compensation_base_rhs =
        nonlinear_current_compensation_enabled ?
        (
            reset_current_source_values ?
            Float64.(over16_state.source.f_values) :
            _source_rhs_without_nonlinear_current_compensation(
                over16_state.source.f_values,
                dense_primary_nonlinear_compensation ?
                    over16_state.source.nonlinear_current_compensation_values :
                    nonlinear_current_compensation_injections,
            )
        ) :
        Float64[]
    if dense_primary_nonlinear_compensation
        nonlinear_current_compensation_base_rhs =
            _dense_primary_nonlinear_reference_rhs(
                nonlinear_current_compensation_base_rhs,
                nonlinear_current_config,
                context.system.node_count,
            )
    end
    nonlinear_pre_solve_update = nothing
    nonlinear_slope_sync_result = (matched_count = 0, mutation_count = 0)
    if nonlinear_current_compensation_enabled &&
       _has_pseudo_nonlinear_inductor_current(nonlinear_current_config) &&
       _has_live_saturated_transformer_nonlinear_slope_branch(context) &&
       get(nonlinear_current_config, :seed_initial_nonlinear_state, false) &&
       !dense_primary_nonlinear_compensation
        nonlinear_pre_solve_rhs = copy(nonlinear_current_compensation_base_rhs)
        _add_current_injection_prefix_values!(
            nonlinear_pre_solve_rhs,
            nonlinear_current_compensation_injections,
        )
        nonlinear_pre_solve_kwargs = _prepare_over16_step_kwargs!(
            context,
            over16_state,
            context.system.v,
            (
                # This boundary mirrors only the SUBTS1 nonlinear-current
                # mutation. Pending switch operations remain owned by the
                # normal accepted-timestep pass below and must run once.
                switch_operation_enabled = false,
                nonlinear_current_config = merge(
                    nonlinear_current_config,
                    (
                        nonlinear_current_compensation_base_rhs =
                            nonlinear_current_compensation_base_rhs,
                        rhs = nonlinear_pre_solve_rhs,
                    ),
                ),
            ),
            report_step_s = hybrid_report_step_s,
        )
        nonlinear_pre_solve_update = over16_accepted_timestep_update!(
            over16_state,
            context.system.v;
            collect_diagnostics = collect_step_diagnostics,
            nonlinear_pre_solve_kwargs...,
        )
        nonlinear_slope_sync_result =
            _sync_saturated_transformer_nonlinear_slope_branches!(
                context,
                nonlinear_pre_solve_kwargs.nonlinear_current_config,
                over16_state,
            )
        nonlinear_current_compensation_injections =
            _nonlinear_current_compensation_injection_values(
                over16_state,
                context.system.node_count,
            )
        solver_over16_kwargs =
            _without_nonlinear_current_config(solver_over16_kwargs)
    end
    if nonlinear_current_compensation_enabled
        _add_current_injection_values!(
            current_injections,
            nonlinear_current_compensation_injections,
        )
    end
    if nonlinear_current_compensation_enabled &&
       reset_current_source_values &&
       nonlinear_pre_solve_update === nothing
        _add_current_injection_prefix_values!(
            over16_state.source.f_values,
            nonlinear_current_compensation_injections,
        )
    end
    use_complete_source_rhs = reset_current_source_values
    complete_source_rhs =
        use_complete_source_rhs ? copy(over16_state.source.f_values) : Float64[]
    use_complete_source_rhs && _apply_ideal_transformer_rhs_values!(
        complete_source_rhs,
        source_voltage_constraint_result,
    )
    current_zero_switching =
        current_extinction_enabled(context.system.elements)
    voltage_max_before_solve = maximum(abs, context.system.v; init = 0.0)
    voltage =
        nonlinear_nodal_step_solver !== nothing ?
        _execute_nonlinear_nodal_step_solver(
            nonlinear_nodal_step_solver,
            context,
            current_injections,
            source_voltage_constraint_result,
        ) :
        sparse_node_group_enabled ?
        solve_step_with_sparse_node_groups!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
            node_group_successors = sparse_node_group_config.node_group_successors,
            grouped_switch_from_nodes = sparse_node_group_config.grouped_switch_from_nodes,
            grouped_switch_to_nodes = sparse_node_group_config.grouped_switch_to_nodes,
            complete_rhs_values =
                use_complete_source_rhs ? complete_source_rhs : nothing,
            workspace = context.sparse_node_group_workspace,
        ) :
        current_zero_switching ?
        solve_step_with_switch_state!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            switch_time_s = context.t_s,
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
        ) :
        source_voltage_constraint_result.applied ?
        solve_step_with_current_injections!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections;
            voltage_constraint_nodes = source_voltage_constraint_result.nodes,
            voltage_constraint_values = source_voltage_constraint_result.values,
        ) :
        solve_step_with_current_injections!(
            context.system,
            context.t_s,
            context.dt_s,
            current_injections,
        )
    if !all(isfinite, voltage)
        rhs_max = maximum(abs, context.system.rhs; init = 0.0)
        injection_max = maximum(abs, current_injections; init = 0.0)
        admittance_condition = all(isfinite, context.system.y) ?
            cond(context.system.y) : Inf
        throw(ArgumentError(
            "nodal solve produced non-finite voltage at t=$(t_before) s " *
            "(previous_max=$(voltage_max_before_solve), rhs_max=$(rhs_max), " *
            "injection_max=$(injection_max), admittance_condition=$(admittance_condition), " *
            "nonlinear_segments=$(over16_state.nonlinear_inverse.curr), " *
            "nonlinear_companions=$(over16_state.nonlinear_inverse.anonl), " *
            "nonlinear_flux_state=$(over16_state.nonlinear_inverse.vnonl))",
        ))
    end
    hysteretic_companion_admittance_stamp_count = 0
    if dense_primary_nonlinear_compensation
        hysteretic_companion_admittance_stamp_count =
            _apply_dense_primary_hysteretic_companion_admittance!(
                context,
                nonlinear_current_config,
                over16_state,
            )
        _dense_primary_nonlinear_inverse_columns!(
            over16_state,
            nonlinear_current_config,
            context,
        )
    end
    distributed_transposed_line_result =
        _distributed_transposed_line_deck_history_update!(
            distributed_transposed_line_config,
            voltage,
            distributed_transposed_line_current_result,
            collect_diagnostics = collect_step_diagnostics,
        )
    frequency_dependent_line_result =
        _frequency_dependent_line_deck_history_update!(
            frequency_dependent_line_config,
            voltage,
            frequency_dependent_line_current_result,
        )
    if sync_switch_base_admittance
        haskey(solver_over16_kwargs, :switch_admittance_config) ||
            throw(ArgumentError("step-context switch admittance sync requires switch_admittance_config"))
        _sync_step_context_switch_base_admittance!(
            over16_state,
            context,
            solver_over16_kwargs.switch_admittance_config,
        )
    end
    solver_over16_kwargs =
        _without_over16_step_context_switch_config(solver_over16_kwargs)
    if dense_primary_nonlinear_compensation &&
       haskey(solver_over16_kwargs, :nonlinear_current_config)
        nonlinear_current_config = _dense_primary_nonlinear_steady_state_seed(
            solver_over16_kwargs.nonlinear_current_config,
            over16_state,
            voltage,
        )
        solver_over16_kwargs = merge(
            solver_over16_kwargs,
            (nonlinear_current_config = nonlinear_current_config,),
        )
    end
    prepared_over16_kwargs = _prepare_over16_step_kwargs!(
            context,
            over16_state,
            voltage,
            solver_over16_kwargs,
            report_step_s = hybrid_report_step_s,
        )
    if nonlinear_current_compensation_enabled &&
       haskey(prepared_over16_kwargs, :nonlinear_current_config)
        prepared_over16_kwargs = merge(
            prepared_over16_kwargs,
            (
                nonlinear_current_config = merge(
                    prepared_over16_kwargs.nonlinear_current_config,
                    (
                        nonlinear_current_compensation_base_rhs =
                            nonlinear_current_compensation_base_rhs,
                        rhs =
                            dense_primary_nonlinear_compensation ?
                            nonlinear_current_compensation_base_rhs :
                            get(
                                prepared_over16_kwargs.nonlinear_current_config,
                                :rhs,
                                over16_state.source.f_values,
                            ),
                    ),
                ),
            ),
        )
    end
    prepared_over16_kwargs = _without_current_injection_values(prepared_over16_kwargs)
    sparse_switch_state_flow_result = nothing
    lean_sparse_switch_state_flow_result = nothing
    if _over16_step_sparse_switch_state_flow_enabled(prepared_over16_kwargs)
        sparse_step_config =
            _without_over16_sparse_switch_state_flow_config(prepared_over16_kwargs)
        sparse_max_passes = Int(get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_max_passes,
            4,
        ))
        sparse_initial_switch_operation_enabled = get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_initial_switch_operation_enabled,
            nothing,
        )
        sparse_repeat_after_post_current_queue = get(
            prepared_over16_kwargs,
            :sparse_switch_state_flow_repeat_after_post_current_queue,
            true,
        )
        if collect_step_diagnostics
            sparse_switch_state_flow_result = over16_sparse_switch_state_flow_update!(
                over16_state,
                voltage,
                sparse_step_config;
                collect_diagnostics = true,
                max_passes = sparse_max_passes,
                initial_switch_operation_enabled =
                    sparse_initial_switch_operation_enabled,
                repeat_after_post_current_queue =
                    sparse_repeat_after_post_current_queue,
            )
            over16_update = sparse_switch_state_flow_result.pass_updates[end]
        else
            lean_sparse_switch_state_flow_result =
                _over16_sparse_switch_state_flow_update_lean!(
                    over16_state,
                    voltage,
                    sparse_step_config;
                    max_passes = sparse_max_passes,
                    initial_switch_operation_enabled =
                        sparse_initial_switch_operation_enabled,
                    repeat_after_post_current_queue =
                        sparse_repeat_after_post_current_queue,
                )
            over16_update = lean_sparse_switch_state_flow_result.last_update
        end
    else
        step_config = _without_over16_sparse_switch_state_flow_config(prepared_over16_kwargs)
        over16_update = over16_accepted_timestep_update!(
            over16_state,
            voltage;
            collect_diagnostics = collect_step_diagnostics,
            step_config...,
        )
    end
    if dense_primary_nonlinear_compensation
        switched_topology_applied =
            sparse_switch_state_flow_result !== nothing ||
            lean_sparse_switch_state_flow_result !== nothing
        switched_topology_applied &&
            _apply_switched_topology_admittance!(context, over16_state)
        if any(
            ==(PSEUDO_NONLINEAR_INDUCTOR_TYPE),
            get(nonlinear_current_config, :nonlinear_types, Int[]),
        ) && _has_live_saturated_transformer_nonlinear_slope_branch(context)
            nonlinear_slope_sync_result =
                _sync_saturated_transformer_nonlinear_slope_branches!(
                    context,
                    nonlinear_current_config,
                    over16_state,
                )
        end
        voltage = _apply_dense_primary_nonlinear_solution!(
            context,
            nonlinear_current_config,
            over16_update.nonlinear_current_result,
            nonlinear_current_compensation_base_rhs,
            stored_nonlinear_current_compensation_injections,
            source_voltage_constraint_result,
        )
        if switched_topology_applied
            switch_current_config = get(
                solver_over16_kwargs,
                :switch_current_config,
                nothing,
            )
            voltage = _sync_switched_nonlinear_network_solution!(
                context,
                over16_state,
                switch_current_config,
            )
            over16_update = merge(
                over16_update,
                (
                    output_switch_network_solution =
                        copy(over16_state.switch_current.network_solution),
                    output_switch_currents =
                        copy(over16_state.switch_current.switch_currents),
                ),
            )
        end
    end
    over16_state_mutated = if sparse_switch_state_flow_result !== nothing
        sparse_switch_state_flow_result.accepted_timestep_state_mutation_count > 0
    elseif lean_sparse_switch_state_flow_result !== nothing
        lean_sparse_switch_state_flow_result.accepted_timestep_state_mutation_count > 0
    else
        over16_update.accepted_timestep_state_mutated
    end
    if nonlinear_pre_solve_update !== nothing
        over16_state_mutated =
            over16_state_mutated ||
            nonlinear_pre_solve_update.accepted_timestep_state_mutated
        over16_update = merge(
            over16_update,
            (
                nonlinear_current_result =
                    nonlinear_pre_solve_update.nonlinear_current_result,
                nonlinear_pre_solve_applied = true,
                accepted_timestep_state_mutated = over16_state_mutated,
            ),
        )
    end
    _accept_source_function_boundary_update!(
        context,
        over16_state,
        source_config,
        over16_update,
        idempotent_same_time = hybrid_substep,
    )
    trace_voltage = record_presolve_voltage_state ? presolve_voltage : voltage
    if hybrid_substep
        _update_deck_power_energy_state!(context, voltage)
    else
        record_step!(context, voltage)
    end
    collect_step_diagnostics || return nothing
    return (
        source = :emt_step_context_over16_boundary,
        outcome = :timestep_integration,
        step_index = step_before,
        t_s = t_before,
        dt_s = context.dt_s,
        voltage_pu = copy(voltage),
        trace_voltage_pu = copy(trace_voltage),
        presolve_voltage_state_recorded = record_presolve_voltage_state,
        over16_update = over16_update,
        over16_sparse_switch_state_flow_result = sparse_switch_state_flow_result,
        over16_sparse_switch_state_flow_applied =
            sparse_switch_state_flow_result !== nothing,
        step_context_recorded = !hybrid_substep,
        step_context_step_after = context.step_index,
        step_context_time_after = context.t_s,
        over16_state_mutated = over16_state_mutated,
        accepted_timestep_state_mutated = over16_state_mutated,
        saturated_transformer_nonlinear_slope_branch_count =
            nonlinear_slope_sync_result.matched_count,
        saturated_transformer_nonlinear_slope_mutation_count =
            nonlinear_slope_sync_result.mutation_count,
        nonlinear_current_pre_solve_applied =
            nonlinear_pre_solve_update !== nothing,
        source_voltage_constraint_result = source_voltage_constraint_result,
        source_voltage_constraint_count = source_voltage_constraint_result.count,
        source_voltage_constraints_applied =
            source_voltage_constraint_result.applied,
        dynamic_source_row_update_count =
            source_function_row_result === nothing ?
            0 :
            length(context.source_function_runtime.dynamic_row_indices),
        current_injection_count = length(current_injections),
        current_injection_nonzero_count = count(!=(0.0), current_injections),
        nonlinear_current_compensation_carry_values =
            copy(nonlinear_current_compensation_injections),
        hysteretic_companion_admittance_stamp_count =
            hysteretic_companion_admittance_stamp_count,
        nonlinear_current_compensation_carry_nonzero_count =
            count(!=(0.0), nonlinear_current_compensation_injections),
        nonlinear_current_compensation_base_rhs_values =
            copy(nonlinear_current_compensation_base_rhs),
        current_source_seeded = reset_current_source_values && seed_current_source_values,
        electromagnetic_history_rhs_values = copy(electromagnetic_history_rhs),
        electromagnetic_history_rhs_nonzero_count =
            count(!=(0.0), electromagnetic_history_rhs),
        current_source_seed_voltage_constraint_zero_count =
            current_source_seed_voltage_constraint_zero_count,
        distributed_transposed_line_result = distributed_transposed_line_result,
        distributed_transposed_line_update_count =
            distributed_transposed_line_result === nothing ? 0 :
            distributed_transposed_line_result.line_update_count,
        distributed_transposed_line_rhs_update_count =
            distributed_transposed_line_result === nothing ? 0 :
            distributed_transposed_line_result.rhs_update_count,
        distributed_transposed_line_state_mutated =
            distributed_transposed_line_result !== nothing &&
            distributed_transposed_line_result.state_mutated,
        frequency_dependent_line_result = frequency_dependent_line_result,
        frequency_dependent_line_update_count =
            frequency_dependent_line_result === nothing ? 0 :
            frequency_dependent_line_result.line_update_count,
        frequency_dependent_line_rhs_update_count =
            frequency_dependent_line_result === nothing ? 0 :
            frequency_dependent_line_result.rhs_update_count,
        frequency_dependent_line_state_mutated =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.state_mutated,
        frequency_dependent_line_runtime_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_dependent_line_runtime_executed,
        frequency_dependent_line_skin_effect_internal_impedance_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.skin_effect_internal_impedance_executed,
        frequency_dependent_line_earth_return_impedance_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.earth_return_impedance_executed,
        frequency_dependent_line_fitting_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_dependent_fitting_executed,
        frequency_dependent_line_frequency_loop_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.frequency_loop_executed,
        frequency_dependent_line_pipe_sheath_side_effects_executed =
            frequency_dependent_line_result !== nothing &&
            frequency_dependent_line_result.pipe_sheath_side_effects_executed,
        legacy_fortran_in_loop = false,
        full_bpa_timestep_executed = false,
        full_deck_orchestration_executed = false,
        deferred_effects = (
            :full_init_step_calc_elec_orchestration,
            :deck_card_io,
            :solvum,
            :report_file_writers,
            :external_bpa_executable_waveform_comparison,
        ),
    )
end

function _apply_steady_state_initial_sample!(context::EMTStepContext, sample)
    return _apply_steady_state_initial_sample!(context, sample, nothing)
end

function _apply_steady_state_initial_sample!(
    context::EMTStepContext,
    sample,
    output_values::Union{Nothing,AbstractVector{<:Real}},
)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state initial sample must cover context nodes"))
    voltage = Float64[values[node] for node in 1:context.system.node_count]
    _update_deck_power_energy_state!(context, voltage)
    output_sample = @view context.output_step_values[:, 1]
    if !isempty(output_sample)
        if output_values === nothing
            _record_context_outputs!(context.output_step_values, 1, context, voltage)
        else
            length(output_values) == length(context.output_channel_names) ||
                throw(ArgumentError("steady-state initial output sample length mismatch"))
            for index in eachindex(output_values)
                output_sample[index] = Float64(output_values[index])
            end
        end
    end
    _update_context_extrema!(context, voltage, output_sample)
    column = _recorded_trace_column!(context, 0)
    if column > 0
        context.time_s[column] = 0.0
        for node in 1:context.system.node_count
            context.voltage_pu[node, column] = voltage[node]
        end
        context.output_pu[:, column] .= output_sample
    end
    return context
end

function _steady_state_initial_output_values(context::EMTStepContext)
    isempty(context.output_channel_names) && return Float64[]
    values = zeros(Float64, length(context.output_channel_names), 1)
    _record_context_outputs!(values, 1, context, context.system.v)
    return vec(values)
end

function _steady_state_node_voltage_phasor(sample, node::Int)
    node == 0 && return complex(0.0, 0.0)
    hasproperty(sample, :node_voltage_phasors) ||
        throw(ArgumentError("steady-state sample must include node voltage phasors"))
    1 <= node <= length(sample.node_voltage_phasors) ||
        throw(ArgumentError("steady-state phasor sample must cover branch nodes"))
    return sample.node_voltage_phasors[node]
end

function _steady_state_branch_frequency_hz(
    sample,
    from_node::Int,
    to_node::Int,
    default_frequency_hz::Float64,
)
    hasproperty(sample, :node_steady_state_frequencies_hz) ||
        return default_frequency_hz
    node_frequencies_hz = sample.node_steady_state_frequencies_hz
    selected_frequency_hz = nothing
    for node in (from_node, to_node)
        node == 0 && continue
        # Transformer and machine assembly can append internal nodes after the
        # parsed deck partition is formed. Their mixed-frequency paths are
        # rejected before assembly, so those nodes inherit the deck default.
        1 <= node <= length(node_frequencies_hz) || continue
        node_frequency_hz = Float64(node_frequencies_hz[node])
        node_frequency_hz > 0.0 || continue
        if selected_frequency_hz === nothing
            selected_frequency_hz = node_frequency_hz
        elseif selected_frequency_hz != node_frequency_hz
            throw(ArgumentError("steady-state branch endpoints have different frequencies"))
        end
    end
    return selected_frequency_hz === nothing ?
        default_frequency_hz : selected_frequency_hz
end

function _steady_state_reactive_angular_frequency(sample, frequency_hz::Float64)
    physical_angular_frequency = 2.0 * pi * frequency_hz
    formulation = hasproperty(sample, :harmonic_formulation) ?
        getproperty(sample, :harmonic_formulation) : :physical_frequency
    if formulation === :physical_frequency
        return physical_angular_frequency
    elseif formulation === :timestep_matched
        hasproperty(sample, :timestep_s) || throw(ArgumentError(
            "timestep-matched steady-state sample must declare timestep_s",
        ))
        return _emt_reactive_angular_frequency(
            TimestepMatchedFormulation(Float64(sample.timestep_s)),
            physical_angular_frequency,
        )
    end
    throw(ArgumentError("unsupported steady-state harmonic formulation $formulation"))
end

function _seed_steady_state_series_rl_branch!(
    branch::SeriesRLBranch,
    sample,
    frequency_hz::Float64,
)
    isfinite(frequency_hz) && frequency_hz >= 0.0 ||
        throw(ArgumentError("series R-L steady-state frequency must be finite and nonnegative"))
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    all(
        isfinite,
        (
            real(branch_voltage_phasor),
            imag(branch_voltage_phasor),
        ),
    ) || throw(ArgumentError("series R-L steady-state voltage must be finite"))
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    impedance =
        branch.l <= 0.0 ?
        complex(branch.r, 0.0) :
        complex(branch.r, reactive_angular_frequency * branch.l)
    current_phasor =
        abs(impedance) == 0.0 ? complex(0.0, 0.0) : branch_voltage_phasor / impedance
    all(
        isfinite,
        (
            real(current_phasor),
            imag(current_phasor),
        ),
    ) || throw(ArgumentError("series R-L steady-state current must be finite"))
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    return branch
end

function _seed_steady_state_series_rlc_branch!(
    branch::SeriesRLCBranch,
    sample,
    frequency_hz::Float64,
)
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    omega = _steady_state_reactive_angular_frequency(sample, frequency_hz)
    if omega <= 0.0
        branch.v_prev = real(branch_voltage_phasor)
        branch.i_prev = 0.0
        branch.i_last = 0.0
        branch.inductor_voltage_prev = 0.0
        branch.capacitor_voltage_prev = branch.v_prev
        return branch
    end
    impedance = complex(branch.r, omega * branch.l - inv(omega * branch.c))
    current_phasor =
        abs(impedance) == 0.0 ? complex(0.0, 0.0) : branch_voltage_phasor / impedance
    inductor_voltage_phasor = im * omega * branch.l * current_phasor
    capacitor_voltage_phasor = current_phasor / (im * omega * branch.c)
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    branch.inductor_voltage_prev = real(inductor_voltage_phasor)
    branch.capacitor_voltage_prev = real(capacitor_voltage_phasor)
    return branch
end

function _seed_steady_state_capacitor_branch!(
    branch::CapacitorBranch,
    sample,
    frequency_hz::Float64,
)
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    current_phasor = im * reactive_angular_frequency * branch.c * branch_voltage_phasor
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    return branch
end

function _coupled_inductive_steady_state_admittance(
    branch::CoupledInductiveBranch,
    angular_frequency::Float64=branch.angular_frequency,
)
    isfinite(angular_frequency) && angular_frequency > 0.0 || throw(ArgumentError(
        "coupled-inductive steady-state angular frequency must be finite and positive",
    ))
    frequency_scale = branch.angular_frequency / angular_frequency
    scaled_susceptance = frequency_scale .* branch.susceptance
    admittance = im .* scaled_susceptance
    if branch.series_resistance > 0.0
        reference_susceptance = scaled_susceptance[
            branch.resistance_reference_port,
            branch.resistance_reference_port,
        ]
        reference_reactance = -inv(reference_susceptance)
        admittance .*= (im * reference_reactance) /
                       (branch.series_resistance + im * reference_reactance)
    end
    return admittance
end

function _seed_steady_state_coupled_inductive_branch!(
    branch::CoupledInductiveBranch,
    sample,
)
    port_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, branch.a[index]) -
        _steady_state_node_voltage_phasor(sample, branch.b[index])
        for index in eachindex(branch.a)
    ]
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    current_phasors = _coupled_inductive_steady_state_admittance(
        branch,
        reactive_angular_frequency,
    ) *
                      port_voltage_phasors
    branch.previous_voltage .= real.(port_voltage_phasors)
    branch.previous_current .= real.(current_phasors)
    branch.last_current .= branch.previous_current
    return branch
end

function _seed_steady_state_coupled_series_rl_branch!(
    branch::CoupledSeriesRLBranch,
    sample,
    frequency_hz::Float64,
)
    port_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, branch.a[index]) -
        _steady_state_node_voltage_phasor(sample, branch.b[index])
        for index in eachindex(branch.a)
    ]
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    impedance = complex.(
        branch.resistance_matrix,
        reactive_angular_frequency .* branch.inductance_matrix,
    )
    current_phasors = impedance \ port_voltage_phasors
    branch.previous_voltage .= real.(port_voltage_phasors)
    branch.previous_current .= real.(current_phasors)
    branch.last_current .= branch.previous_current
    return branch
end

function _sequence_modal_phasors(phase_phasors::AbstractVector{<:Complex})
    nph = length(phase_phasors)
    nph > 0 || throw(ArgumentError("phase_phasors must not be empty"))
    modal = Vector{ComplexF64}(undef, nph)
    sum_phasor = complex(0.0, 0.0)
    for phase in 1:nph
        sum_phasor += phase_phasors[phase]
    end
    modal[1] = sum_phasor / nph
    first_phase = phase_phasors[1]
    for phase in 2:nph
        modal[phase] = (first_phase - phase_phasors[phase]) / nph
    end
    return modal
end

function _frequency_domain_sequence_admittance(record, angular_frequency::Float64)
    resistance = Float64(record.r)
    inductive_reactance = Float64(record.l) * angular_frequency
    damping_resistance = Float64(record.rl)
    capacitive_reactance =
        record.c > 0.0 ? inv(Float64(record.c) * angular_frequency) : 0.0
    real_impedance = resistance
    imaginary_impedance = inductive_reactance
    if damping_resistance != 0.0 && inductive_reactance != 0.0
        denominator = inv(damping_resistance^2 + inductive_reactance^2)
        real_impedance += damping_resistance * inductive_reactance^2 * denominator
        imaginary_impedance = inductive_reactance * damping_resistance^2 * denominator
    end
    imaginary_impedance -= capacitive_reactance
    denominator = inv(real_impedance^2 + imaginary_impedance^2)
    admittance = complex(
        real_impedance * denominator,
        -imaginary_impedance * denominator,
    )
    return admittance, capacitive_reactance
end

function _seed_lumped_sequence_frequency_histories!(
    element::BreqivHistoryInjection,
    sample,
    dt_s::Float64,
    frequency_hz::Float64,
)
    angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    branch_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, element.a[phase]) -
        _steady_state_node_voltage_phasor(sample, element.b[phase])
        for phase in eachindex(element.a)
    ]
    seed_breqiv_frequency_histories!(
        element,
        branch_voltage_phasors,
        angular_frequency,
    )
    initialize_breqiv_history_injection!(element, dt_s)
    return element
end

function _seed_semlyen_line_steady_state!(
    element::SemlyenFrequencyDependentLine,
    sample,
    default_frequency_hz::Float64,
)
    frequencies = Float64[
        _steady_state_branch_frequency_hz(
            sample,
            element.from_nodes[phase],
            element.to_nodes[phase],
            default_frequency_hz,
        )
        for phase in eachindex(element.from_nodes)
        if element.from_nodes[phase] != 0 || element.to_nodes[phase] != 0
    ]
    frequency_hz = isempty(frequencies) ? default_frequency_hz : first(frequencies)
    all(value -> isapprox(value, frequency_hz; atol = 1.0e-9, rtol = 1.0e-9), frequencies) ||
        throw(ArgumentError("Semlyen line phases have different steady-state frequencies"))
    from_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.from_nodes
    ]
    to_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.to_nodes
    ]
    initialize_semlyen_line_steady_state!(
        element,
        from_phasors,
        to_phasors,
        frequency_hz,
    )
    return element
end

function _seed_complex_modal_line_steady_state!(
    element::ComplexModalBergeronLine,
    sample,
    default_frequency_hz::Float64,
)
    frequencies = Float64[
        _steady_state_branch_frequency_hz(
            sample,
            element.from_nodes[phase],
            element.to_nodes[phase],
            default_frequency_hz,
        )
        for phase in eachindex(element.from_nodes)
        if element.from_nodes[phase] != 0 || element.to_nodes[phase] != 0
    ]
    frequency_hz = isempty(frequencies) ? default_frequency_hz : first(frequencies)
    all(value -> isapprox(value, frequency_hz; atol = 1.0e-9, rtol = 1.0e-9), frequencies) ||
        throw(ArgumentError("complex modal line phases have different steady-state frequencies"))
    from_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.from_nodes
    ]
    to_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.to_nodes
    ]
    initialize_complex_modal_bergeron_steady_state!(
        element,
        from_phasors,
        to_phasors,
        frequency_hz,
    )
    return element
end

function _steady_state_current_injections(context::EMTStepContext, sample)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state current injection seed must cover context nodes"))
    return nodal_current_injections_for_voltage!(
        context.system,
        0.0,
        context.dt_s,
        Float64.(values[1:context.system.node_count]),
    )
end

function _seed_steady_state_network_state!(context::EMTStepContext, sample)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state network seed must cover context nodes"))
    for node in 1:context.system.node_count
        context.system.v[node] = Float64(values[node])
    end
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    for element in context.system.elements
        if element isa ComplexModalBergeronLine &&
           hasproperty(sample, :node_voltage_phasors)
            _seed_complex_modal_line_steady_state!(element, sample, frequency_hz)
            continue
        elseif element isa SemlyenFrequencyDependentLine &&
           hasproperty(sample, :node_voltage_phasors)
            _seed_semlyen_line_steady_state!(element, sample, frequency_hz)
            continue
        elseif element isa BreqivHistoryInjection
            if hasproperty(sample, :node_voltage_phasors)
                element_frequency_hz = _steady_state_branch_frequency_hz(
                    sample,
                    first(element.a),
                    first(element.b),
                    frequency_hz,
                )
                _seed_lumped_sequence_frequency_histories!(
                    element,
                    sample,
                    context.dt_s,
                    element_frequency_hz,
                )
                # The recorded steady-state sample replaces the t=0 solve. Advance
                # the BREQIV history once so its Norton source is staged for the
                # first dynamic solve at t=dt, like the other companion branches.
                advance_breqiv_history_current!(
                    element,
                    context.system.v,
                    context.dt_s;
                    consumed_for_step = false,
                )
            else
                for phase in eachindex(element.initial_phase_voltage)
                    from_node = element.a[phase]
                    to_node = element.b[phase]
                    from_voltage = from_node == 0 ? 0.0 : context.system.v[from_node]
                    to_voltage = to_node == 0 ? 0.0 : context.system.v[to_node]
                    element.initial_phase_voltage[phase] = from_voltage - to_voltage
                end
            end
            continue
        elseif element isa SeriesRLBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_series_rl_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa SeriesRLCBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_series_rlc_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa CapacitorBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_capacitor_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa CoupledInductiveBranch && hasproperty(sample, :node_voltage_phasors)
            _seed_steady_state_coupled_inductive_branch!(element, sample)
            continue
        elseif element isa CoupledSeriesRLBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                first(element.a),
                first(element.b),
                frequency_hz,
            )
            _seed_steady_state_coupled_series_rl_branch!(
                element,
                sample,
                element_frequency_hz,
            )
            continue
        end
        if !hasproperty(sample, :node_voltage_phasors)
            if element isa SeriesRLBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                continue
            elseif element isa SeriesRLCBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                element.inductor_voltage_prev = 0.0
                element.capacitor_voltage_prev = element.v_prev
                continue
            elseif element isa CapacitorBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                continue
            elseif element isa CoupledInductiveBranch
                for port in eachindex(element.a)
                    element.previous_voltage[port] =
                        _deck_node_voltage(context.system.v, element.a[port]) -
                        _deck_node_voltage(context.system.v, element.b[port])
                end
                fill!(element.previous_current, 0.0)
                fill!(element.last_current, 0.0)
                continue
            elseif element isa CoupledSeriesRLBranch
                for port in eachindex(element.a)
                    element.previous_voltage[port] =
                        _deck_node_voltage(context.system.v, element.a[port]) -
                        _deck_node_voltage(context.system.v, element.b[port])
                end
                fill!(element.previous_current, 0.0)
                fill!(element.last_current, 0.0)
                continue
            end
        end
        update!(element, context.system.v, context.dt_s)
    end
    return context
end

function _seed_direct_machine_power_leakage_currents!(
    context::EMTStepContext,
    parsed::DeckParser.DeckParseResult,
    machine_indices::AbstractVector{<:Integer},
    power_terminal_currents::AbstractVector{<:Real},
)
    indices = Int.(machine_indices)
    currents = Float64.(power_terminal_currents)
    length(indices) == length(currents) ||
        throw(ArgumentError("direct-machine leakage-current seeds must align with machine indices"))
    for (machine_index, current) in zip(indices, currents)
        rows = [
            row
            for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
            if row.machine_index == machine_index &&
               row.from_node_value != row.to_node_value
        ]
        length(rows) == 1 ||
            throw(ArgumentError("automatic direct machine $machine_index requires one generated power-leakage branch"))
        row = only(rows)
        matches = SeriesRLBranch[
            element
            for element in context.system.elements
            if element isa SeriesRLBranch &&
               element.a == row.from_node_value &&
               element.b == row.to_node_value
        ]
        length(matches) == 1 ||
            throw(ArgumentError("automatic direct machine $machine_index power-leakage runtime branch is missing or ambiguous"))
        branch = only(matches)
        branch.i_prev = current
        branch.i_last = current
        branch.v_prev = branch.r * current
    end
    return context
end

function _deck_synchronous_machine_delta_connected(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    return any(
        row -> row.machine_index == machine_index &&
               row.parameter_kind == :delta_connection,
        DeckParser.deck_synchronous_machine_model_parameter_rows(parsed),
    )
end

function _synchronous_machine_winding_voltages(
    terminal_node_voltages::AbstractVector{<:Real},
    delta_connected::Bool,
)
    length(terminal_node_voltages) == 3 || throw(ArgumentError(
        "synchronous-machine terminal voltage must contain three phases",
    ))
    values = Float64.(terminal_node_voltages)
    delta_connected || return values
    return Float64[
        values[1] - values[2],
        values[2] - values[3],
        values[3] - values[1],
    ]
end

function _synchronous_machine_terminal_currents(
    winding_currents::AbstractVector{<:Real},
    delta_connected::Bool,
)
    length(winding_currents) == 3 || throw(ArgumentError(
        "synchronous-machine winding current must contain three phases",
    ))
    values = Float64.(winding_currents)
    delta_connected || return values
    return Float64[
        values[1] - values[3],
        values[2] - values[1],
        values[3] - values[2],
    ]
end

function _deck_synchronous_machine_terminal_admittance(
    parsed::DeckParser.DeckParseResult,
    state::SynchronousMachineDynamicState,
    machine_index::Int=1,
    delta_reference_phase_admittance::Union{Nothing,Real}=nothing,
)
    direct = Float64(state.electrical_coefficients[27])
    mutual = Float64(state.electrical_coefficients[28])
    delta_connected = _deck_synchronous_machine_delta_connected(parsed, machine_index)
    if delta_connected
        phase_admittance = direct - mutual
        if delta_reference_phase_admittance !== nothing
            # PAST retains the initial delta companion in the network base;
            # UPDATE maps only later winding-admittance changes through the delta.
            reference = Float64(delta_reference_phase_admittance)
            phase_admittance = reference + (phase_admittance - reference) / 3.0
        end
        terminal_admittance = fill(-phase_admittance, 3, 3)
        for phase in 1:3
            terminal_admittance[phase, phase] = 2.0 * phase_admittance
        end
        return terminal_admittance
    end
    terminal_admittance = fill(mutual, 3, 3)
    for phase in 1:3
        terminal_admittance[phase, phase] = direct
    end
    return terminal_admittance
end

function _over16_step_kwargs(over16_step_configs, context::EMTStepContext)
    config =
        over16_step_configs === nothing ? NamedTuple() :
        over16_step_configs isa Function ? over16_step_configs(context) :
        over16_step_configs isa AbstractVector ?
            (context.step_index + 1 <= length(over16_step_configs) ?
                over16_step_configs[context.step_index + 1] : NamedTuple()) :
        over16_step_configs
    config === nothing && return NamedTuple()
    config isa NamedTuple ||
        throw(ArgumentError("over16_step_configs entries must be NamedTuple or nothing"))
    return config
end

function _deck_reference_node_name(name::Symbol)
    normalized = lowercase(strip(String(name)))
    return normalized in ("", "0", "gnd", "ground", "ref")
end

function _deck_saturated_transformer_node_index(
    node_map::AbstractDict{Symbol,<:Integer},
    node::Symbol,
)
    _deck_reference_node_name(node) && return 0
    index = get(node_map, node, nothing)
    index !== nothing ||
        throw(ArgumentError("saturated transformer winding node $(String(node)) is not present in the runtime node map"))
    return Int(index)
end

function _deck_saturated_transformer_node_index(
    parsed::DeckParser.DeckParseResult,
    node::Symbol,
)
    return _deck_saturated_transformer_node_index(parsed.node_map, node)
end

function _deck_saturated_transformer_winding_node_names(arrays, winding_number::Int)
    winding_number > 0 ||
        throw(ArgumentError("saturated transformer winding_number must be positive"))
    from_nodes = Symbol[]
    to_nodes = Symbol[]
    for transformer_name in arrays.transformer_names
        matches = findall(eachindex(arrays.winding_transformer_names)) do index
            arrays.winding_transformer_names[index] == transformer_name &&
                arrays.winding_numbers[index] == winding_number
        end
        length(matches) == 1 ||
            throw(ArgumentError("saturated transformer $(String(transformer_name)) requires exactly one winding $winding_number row"))
        winding_index = only(matches)
        push!(from_nodes, arrays.winding_from_nodes[winding_index])
        push!(to_nodes, arrays.winding_to_nodes[winding_index])
    end
    return from_nodes, to_nodes
end

function _deck_saturated_transformer_sparse_node_index(
    node::Int;
    reference_node_index::Int = 0,
    sparse_reference_node_index::Int = 1,
)
    node >= 0 || throw(ArgumentError("saturated transformer sparse node source index must be nonnegative"))
    reference_node_index >= 0 ||
        throw(ArgumentError("saturated transformer reference_node_index must be nonnegative"))
    sparse_reference_node_index >= 1 ||
        throw(ArgumentError("saturated transformer sparse_reference_node_index must be positive"))
    node == reference_node_index && return sparse_reference_node_index
    return node + sparse_reference_node_index - reference_node_index
end

function _deck_saturated_transformer_sparse_config(
    sparse_config::Union{Nothing,NamedTuple},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
)
    length(from_nodes) == length(to_nodes) ||
        throw(ArgumentError("saturated transformer sparse endpoint lengths must match"))
    base_config =
        sparse_config === nothing ?
        (
            derive_from_step_admittance = true,
            commit_to_sparse_network_state = true,
            factor_after_restamp = true,
        ) :
        sparse_config
    has_from_nodes = haskey(base_config, :from_nodes)
    has_to_nodes = haskey(base_config, :to_nodes)
    has_from_nodes == has_to_nodes ||
        throw(ArgumentError("saturated transformer sparse config must provide both from_nodes and to_nodes"))
    has_from_nodes && return base_config
    return merge(
        base_config,
        (
            from_nodes = [
                _deck_saturated_transformer_sparse_node_index(node)
                for node in from_nodes
            ],
            to_nodes = [
                _deck_saturated_transformer_sparse_node_index(node)
                for node in to_nodes
            ],
            current_reference_node_index = 0,
            sparse_reference_node_index = 1,
        ),
    )
end

function _reference_admittance_sparse_rows(
    admittance::AbstractMatrix{<:Real};
    zero_tolerance::Real = 0.0,
)
    return _reference_admittance_sparse_rows!(
        Int[],
        Float64[],
        Int[],
        admittance;
        zero_tolerance = zero_tolerance,
    )
end

function _reference_admittance_sparse_rows!(
    km::Vector{Int},
    ykm::Vector{Float64},
    kks::Vector{Int},
    admittance::AbstractMatrix{<:Real};
    zero_tolerance::Real = 0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("saturated transformer step admittance must be square"))
    node_count = size(admittance, 1)
    tolerance = Float64(zero_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("saturated transformer sparse zero tolerance must be finite and nonnegative"))
    empty!(km)
    empty!(ykm)
    resize!(kks, node_count)
    for row in 1:node_count
        for column in node_count:-1:1
            value = Float64(admittance[row, column])
            isfinite(value) ||
                throw(ArgumentError("saturated transformer step admittance entries must be finite"))
            abs(value) > tolerance || continue
            push!(km, column == row ? -row : column)
            push!(ykm, value)
        end
        kks[row] = length(km) + 1
    end
    return (
        km = km,
        ykm = ykm,
        kks = kks,
        row_count = node_count,
        entry_count = length(km),
    )
end

function _deck_saturated_transformer_intake(
    saturated_transformer_intake,
    saturated_transformer_deck_lines,
    source::AbstractString,
)
    if saturated_transformer_intake !== nothing &&
       saturated_transformer_deck_lines !== nothing
        throw(ArgumentError("provide saturated_transformer_intake or saturated_transformer_deck_lines, not both"))
    end
    saturated_transformer_deck_lines === nothing && return saturated_transformer_intake
    return DeckParser.parse_saturated_transformer_branch_section_intake_lines(
        saturated_transformer_deck_lines;
        source = source,
    )
end

function _deck_transformer_branch_shunt_capacitance_rows(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_deck_lines,
)
    if saturated_transformer_deck_lines !== nothing
        return DeckParser.parse_saturated_transformer_branch_section_shunt_capacitance_rows(
            saturated_transformer_deck_lines;
            source = parsed.source,
        )
    end
    source_path = DeckParser.deck_source_path(parsed)
    source_path === nothing &&
        return DeckParser.DeckTransformerBranchShuntCapacitanceRow[]
    isfile(source_path) || return DeckParser.DeckTransformerBranchShuntCapacitanceRow[]
    return DeckParser.parse_saturated_transformer_branch_section_shunt_capacitance_rows(
        readlines(source_path);
        source = parsed.source,
    )
end

function _append_transformer_branch_shunt_capacitance_elements!(
    elements::Vector{Any},
    element_names::Vector{Symbol},
    node_map::Dict{Symbol,Int},
    parsed::DeckParser.DeckParseResult,
    rows,
)
    isempty(rows) && return nothing
    assigned_indices = Set{Int}(values(node_map))
    for (index, row) in enumerate(rows)
        node_index = _deck_add_runtime_node!(
            node_map,
            assigned_indices,
            row.from_node,
        )
        node_index == 0 && continue
        push!(
            elements,
            CapacitorBranch(
                node_index,
                0,
                DeckParser.fixed_card_branch_timestep_capacitance(
                    parsed,
                    Float64(row.capacitance),
                ),
            ),
        )
        push!(
            element_names,
            Symbol("transformer_branch_shunt_capacitance_", string(index)),
        )
    end
    return nothing
end

function _saturated_transformer_internal_top_node_name(transformer_name::Symbol)
    return Symbol(String(transformer_name), "_saturated_transformer_internal_top")
end

function _saturated_transformer_unique_node_name(base_name::Symbol,
                                                 node_map::AbstractDict{Symbol,<:Integer})
    !haskey(node_map, base_name) && return base_name
    suffix = 2
    while true
        candidate = Symbol(String(base_name), "_", suffix)
        !haskey(node_map, candidate) && return candidate
        suffix += 1
    end
end

function _saturated_transformer_augmented_node_map(
    parsed::DeckParser.DeckParseResult,
    branch_assembly,
)
    return _saturated_transformer_augmented_node_map(parsed.node_map, branch_assembly)
end

function _deck_add_runtime_node!(
    node_map::Dict{Symbol,Int},
    assigned_indices::Set{Int},
    node_name::Symbol,
)
    _deck_reference_node_name(node_name) && return 0
    existing = get(node_map, node_name, 0)
    existing != 0 && return existing
    next_index = maximum(assigned_indices; init = 0) + 1
    while next_index in assigned_indices
        next_index += 1
    end
    node_map[node_name] = next_index
    push!(assigned_indices, next_index)
    return next_index
end

function _saturated_transformer_physical_node_map(
    parsed::DeckParser.DeckParseResult,
    arrays::SaturatedTransformerNonlinearArrays,
)
    node_map = Dict{Symbol,Int}(name => Int(index) for (name, index) in parsed.node_map)
    assigned_indices = Set{Int}(values(node_map))
    for index in eachindex(arrays.winding_transformer_names)
        _deck_add_runtime_node!(
            node_map,
            assigned_indices,
            arrays.winding_from_nodes[index],
        )
        _deck_add_runtime_node!(
            node_map,
            assigned_indices,
            arrays.winding_to_nodes[index],
        )
    end
    return node_map
end

function _saturated_transformer_frequency_node_indices(
    physical_node_map::AbstractDict{Symbol,<:Integer},
    arrays::SaturatedTransformerNonlinearArrays,
    maximum_partition_node_index::Integer,
)
    maximum_index = Int(maximum_partition_node_index)
    maximum_index >= 0 ||
        throw(ArgumentError("maximum partition node index must be nonnegative"))
    node_indices = Int[]
    for node_name in Iterators.flatten((
        arrays.winding_from_nodes,
        arrays.winding_to_nodes,
    ))
        _deck_reference_node_name(node_name) && continue
        node_index = get(physical_node_map, node_name, 0)
        1 <= node_index <= maximum_index || continue
        push!(node_indices, Int(node_index))
    end
    return sort!(unique!(node_indices))
end

function _saturated_transformer_augmented_node_map(
    physical_node_map::AbstractDict{Symbol,<:Integer},
    branch_assembly,
)
    node_map = Dict{Symbol,Int}(name => Int(index) for (name, index) in physical_node_map)
    assigned_indices = Set(values(node_map))
    for (transformer_name, node_index) in zip(
        branch_assembly.internal_top_node_names,
        branch_assembly.internal_top_node_indices,
    )
        node_index = Int(node_index)
        node_index > 0 ||
            throw(ArgumentError("saturated transformer internal top node must be positive"))
        node_index in assigned_indices && continue
        base_name = _saturated_transformer_internal_top_node_name(transformer_name)
        node_name = _saturated_transformer_unique_node_name(base_name, node_map)
        node_map[node_name] = node_index
        push!(assigned_indices, node_index)
    end
    node_count = maximum(values(node_map); init = 0)
    for node_index in 1:node_count
        node_index in assigned_indices && continue
        base_name = Symbol("saturated_transformer_internal_node_", node_index)
        node_name = _saturated_transformer_unique_node_name(base_name, node_map)
        node_map[node_name] = node_index
        push!(assigned_indices, node_index)
    end
    return node_map
end

function _saturated_transformer_timestep_inductance(
    reactance::Real,
    reactance_units::Real,
)
    value = Float64(reactance)
    value == 0.0 && return 0.0
    units = Float64(reactance_units)
    isfinite(units) && units > 0.0 ||
        throw(ArgumentError("saturated transformer reactance units must be positive"))
    return value / units
end

function _saturated_transformer_winding_index(
    branch_assembly,
    transformer_name::Symbol,
    winding_number::Int,
)
    for index in eachindex(branch_assembly.winding_transformer_names)
        branch_assembly.winding_transformer_names[index] == transformer_name || continue
        Int(branch_assembly.winding_numbers[index]) == winding_number || continue
        return index
    end
    return 0
end

function _saturated_transformer_ideal_winding_branch(
    branch_assembly,
    primary_index::Int,
    ideal_index::Int,
    reactance_units::Real,
)
    low_reactance = Float64(branch_assembly.series_branch_inductances[ideal_index])
    low_resistance = Float64(branch_assembly.series_branch_resistances[ideal_index])
    primary_turns = Float64(branch_assembly.series_branch_turns[primary_index])
    ideal_turns = Float64(branch_assembly.series_branch_turns[ideal_index])
    low_reactance > 0.0 ||
        throw(ArgumentError("saturated transformer ideal winding reactance must be positive"))
    low_resistance >= 0.0 ||
        throw(ArgumentError("saturated transformer ideal winding resistance must be nonnegative"))
    primary_turns > 0.0 ||
        throw(ArgumentError("saturated transformer primary turns must be positive"))
    ideal_turns > 0.0 ||
        throw(ArgumentError("saturated transformer ideal winding turns must be positive"))
    turns_ratio = primary_turns / ideal_turns
    internal_reactance = low_reactance * turns_ratio^2
    susceptance = [
        -inv(low_reactance) inv(turns_ratio * low_reactance)
        inv(turns_ratio * low_reactance) -inv(internal_reactance)
    ]
    return CoupledInductiveBranch(
        [
            Int(branch_assembly.winding_from_node_indices[ideal_index]),
            Int(branch_assembly.winding_internal_top_node_indices[ideal_index]),
        ],
        [
            Int(branch_assembly.winding_terminal_node_indices[ideal_index]),
            Int(branch_assembly.winding_terminal_node_indices[primary_index]),
        ],
        susceptance,
        Float64(reactance_units),
        series_resistance = low_resistance,
    )
end

function saturated_transformer_branch_elements(
    branch_assembly;
    reactance_units::Real = 2.0 * pi * 60.0,
    primary_winding_number::Int = 1,
    ideal_winding_number::Int = 2,
)
    elements = Any[]
    element_names = Symbol[]
    element_winding_indices = Int[]
    element_magnetizing_indices = Int[]
    handled_windings = falses(length(branch_assembly.winding_transformer_names))

    function push_series!(index::Int)
        push!(
            elements,
            SeriesRLBranch(
                Int(branch_assembly.series_branch_from_node_indices[index]),
                Int(branch_assembly.series_branch_to_node_indices[index]),
                Float64(branch_assembly.series_branch_resistances[index]),
                _saturated_transformer_timestep_inductance(
                    branch_assembly.series_branch_inductances[index],
                    reactance_units,
                ),
            ),
        )
        transformer_name = branch_assembly.winding_transformer_names[index]
        winding_number = Int(branch_assembly.winding_numbers[index])
        push!(
            element_names,
            Symbol(
                "saturated_transformer_",
                String(transformer_name),
                "_winding_",
                string(winding_number),
                "_series",
            ),
        )
        push!(element_winding_indices, index)
        push!(element_magnetizing_indices, 0)
        handled_windings[index] = true
        return nothing
    end

    for transformer_name in branch_assembly.transformer_names
        primary_index = _saturated_transformer_winding_index(
            branch_assembly,
            transformer_name,
            primary_winding_number,
        )
        primary_index == 0 && continue
        push_series!(primary_index)
        ideal_index = _saturated_transformer_winding_index(
            branch_assembly,
            transformer_name,
            ideal_winding_number,
        )
        ideal_index == 0 && continue
        push!(
            elements,
            _saturated_transformer_ideal_winding_branch(
                branch_assembly,
                primary_index,
                ideal_index,
                reactance_units,
            ),
        )
        push!(
            element_names,
            Symbol("saturated_transformer_", String(transformer_name), "_ideal_winding"),
        )
        push!(element_winding_indices, ideal_index)
        push!(element_magnetizing_indices, 0)
        handled_windings[ideal_index] = true
    end

    for index in eachindex(branch_assembly.series_branch_from_node_indices)
        handled_windings[index] && continue
        push_series!(index)
    end
    for index in eachindex(branch_assembly.magnetizing_branch_from_node_indices)
        resistance = Float64(branch_assembly.magnetizing_branch_resistances[index])
        resistance != 0.0 ||
            throw(ArgumentError("saturated transformer magnetizing resistance must be nonzero"))
        push!(
            elements,
            ConductanceBranch(
                Int(branch_assembly.magnetizing_branch_from_node_indices[index]),
                Int(branch_assembly.magnetizing_branch_to_node_indices[index]),
                inv(resistance),
            ),
        )
        push!(
            element_names,
            Symbol("saturated_transformer_magnetizing_branch_", string(index)),
        )
        push!(element_winding_indices, 0)
        push!(element_magnetizing_indices, index)
    end
    return (
        elements = elements,
        element_names = element_names,
        element_count = length(elements),
        element_winding_indices = element_winding_indices,
        element_magnetizing_indices = element_magnetizing_indices,
        series_branch_count = count(!=(0), element_winding_indices),
        magnetizing_branch_count = length(branch_assembly.magnetizing_branch_from_node_indices),
    )
end

function saturated_transformer_augmented_step_context(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_current_config::NamedTuple;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    include_coupled_lumped_sequence_history::Bool = false,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    recorded_step_indices = nothing,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    branch_assembly = saturated_transformer_current_config.saturated_transformer_branch_assembly
    branch_elements = saturated_transformer_branch_elements(
        branch_assembly;
        reactance_units = 2.0 * pi * _deck_steady_state_frequency_hz(parsed),
    )
    nonlinear_slope_branches =
        saturated_transformer_nonlinear_slope_branches(saturated_transformer_current_config)
    node_map = _saturated_transformer_augmented_node_map(parsed, branch_assembly)
    node_count = maximum(values(node_map); init = 0)
    elements = Any[parsed.elements...; branch_elements.elements...]
    element_names = Symbol[parsed.element_names...; branch_elements.element_names...]
    _append_switching_nonlinear_resistor_safety_shunts!(
        elements,
        element_names,
        parsed,
    )
    _append_saturated_transformer_safety_shunts!(
        elements,
        element_names,
        parsed,
        saturated_transformer_current_config,
    )
    append!(elements, nonlinear_slope_branches.elements)
    append!(element_names, nonlinear_slope_branches.element_names)
    _delay_deck_time_switch_events!(elements, time_switch_event_delay_s, t_end_s)
    current_zero_switching &&
        _convert_deck_current_zero_switches!(elements, element_names, parsed)
    if include_coupled_lumped_sequence_history
        source_equivalent = coupled_lumped_sequence_history_injection_elements(parsed)
        append!(elements, source_equivalent.elements)
        append!(element_names, source_equivalent.element_names)
    end
    source_function_runtime, control_system_runtime =
        _append_dynamic_source_and_control_elements!(
        elements,
        element_names,
        parsed,
        dt_s,
        source_signal_provider,
    )
    system = NodalSystem(node_count, elements)
    return initialize_step_context(
        system;
        node_map = node_map,
        element_names = element_names,
        source_function_runtime = source_function_runtime,
        control_system_runtime = control_system_runtime,
        _deck_runtime_output_context_kwargs(
            parsed;
            time_switch_event_delay_s = time_switch_event_delay_s,
            event_horizon_s = t_end_s,
        )...,
        dt_s = dt_s,
        t_end_s = t_end_s,
        source = parsed.source,
        recorded_step_indices = recorded_step_indices,
    )
end

function _deck_source_voltage_guess(
    parsed::DeckParser.DeckParseResult,
    node_count::Int,
    time_s::Float64,
)
    voltages = zeros(Float64, node_count)
    for element in parsed.elements
        element isa TheveninSource || continue
        1 <= element.node <= node_count || continue
        voltages[element.node] = Float64(element.value(time_s))
    end
    return voltages
end

function _branch_phase_voltage(
    voltage::AbstractVector{Float64},
    from_node::Int,
    to_node::Int,
)
    from_voltage = from_node == 0 ? 0.0 : voltage[from_node]
    to_voltage = to_node == 0 ? 0.0 : voltage[to_node]
    return from_voltage - to_voltage
end

function _coupled_lumped_sequence_timestep_inductance(
    parsed::DeckParser.DeckParseResult,
    raw_inductance::Real,
)
    return DeckParser.fixed_card_branch_timestep_inductance(
        parsed,
        Float64(raw_inductance),
    )
end

function coupled_lumped_sequence_history_injection_elements(
    parsed::DeckParser.DeckParseResult;
    initial_time_s::Float64 = 0.0,
)
    DeckParser.assert_deck_valid!(parsed)
    node_count = maximum(values(parsed.node_map); init = 0)
    voltage_guess = _deck_source_voltage_guess(parsed, node_count, initial_time_s)
    elements = Any[]
    element_names = Symbol[]
    element_line_numbers = Int[]
    for impedance in DeckParser.deck_coupled_lumped_sequence_impedances(parsed)
        if impedance.input_kind == :triangular_matrix
            physical_inductance = map(
                value -> _coupled_lumped_sequence_timestep_inductance(parsed, value),
                impedance.phase_inductance_matrix,
            )
            push!(
                elements,
                CoupledSeriesRLBranch(
                    impedance.from_node_indices,
                    impedance.to_node_indices,
                    impedance.phase_resistance_matrix,
                    physical_inductance,
                ),
            )
            push!(element_names, impedance.name)
            push!(
                element_line_numbers,
                isempty(impedance.line_numbers) ? 0 : first(impedance.line_numbers),
            )
            continue
        end
        initial_voltage = [
            _branch_phase_voltage(
                voltage_guess,
                impedance.from_node_indices[index],
                impedance.to_node_indices[index],
            )
            for index in 1:impedance.phase_count
        ]
        push!(
            elements,
            three_phase_breqiv_history_injection(
                impedance.from_node_indices[1],
                impedance.to_node_indices[1],
                impedance.from_node_indices[2],
                impedance.to_node_indices[2],
                impedance.from_node_indices[3],
                impedance.to_node_indices[3],
                impedance.zero_sequence_resistance,
                _coupled_lumped_sequence_timestep_inductance(
                    parsed,
                    impedance.zero_sequence_inductance,
                ),
                0.0,
                0.0,
                impedance.positive_sequence_resistance,
                _coupled_lumped_sequence_timestep_inductance(
                    parsed,
                    impedance.positive_sequence_inductance,
                ),
                0.0,
                0.0,
                initial_voltage[1],
                initial_voltage[2],
                initial_voltage[3];
                history_current_scale = -1.0,
                history_voltage_scale = -1.0,
            ),
        )
        push!(
            element_names,
            Symbol(String(impedance.name), "_source_equivalent_history"),
        )
        push!(
            element_line_numbers,
            isempty(impedance.line_numbers) ? 0 : first(impedance.line_numbers),
        )
    end
    return (
        elements = elements,
        element_names = element_names,
        element_line_numbers = element_line_numbers,
        element_count = length(elements),
    )
end

function saturated_transformer_winding_node_map(
    arrays::SaturatedTransformerNonlinearArrays,
)
    node_map = Dict{Symbol,Int}()
    for index in eachindex(arrays.winding_transformer_names)
        for node_name in (
            arrays.winding_from_nodes[index],
            arrays.winding_to_nodes[index],
        )
            node_name == Symbol("") && continue
            haskey(node_map, node_name) && continue
            node_map[node_name] = length(node_map) + 1
        end
    end
    return node_map
end

function _saturated_transformer_linear_branch_arrays(saturated_transformer_intake)
    transformers = collect(getproperty(saturated_transformer_intake, :transformers))
    windings = collect(getproperty(saturated_transformer_intake, :windings))
    return SaturatedTransformerNonlinearArrays(
        String(getproperty(saturated_transformer_intake, :source)),
        Symbol[getproperty(row, :name) for row in transformers],
        Symbol[getproperty(row, :reference_name) for row in transformers],
        Int[],
        Int[],
        Int[],
        Float64[
            ismissing(getproperty(row, :initial_current)) ?
            0.0 :
            Float64(getproperty(row, :initial_current))
            for row in transformers
        ],
        Float64[
            ismissing(getproperty(row, :initial_flux)) ?
            0.0 :
            Float64(getproperty(row, :initial_flux))
            for row in transformers
        ],
        Union{Missing,Float64}[
            getproperty(row, :magnetizing_resistance)
            for row in transformers
        ],
        Float64[],
        Bool[],
        Symbol[],
        Int[],
        Float64[],
        Float64[],
        Symbol[getproperty(row, :transformer_name) for row in windings],
        Int[Int(getproperty(row, :winding_number)) for row in windings],
        Symbol[getproperty(row, :from_node) for row in windings],
        Symbol[getproperty(row, :to_node) for row in windings],
        Union{Missing,Float64}[getproperty(row, :resistance) for row in windings],
        Union{Missing,Float64}[getproperty(row, :inductance) for row in windings],
        Union{Missing,Float64}[getproperty(row, :turns) for row in windings],
        Bool[Bool(getproperty(row, :inherited_parameters)) for row in windings],
    )
end

function _assert_saturated_transformer_intake_valid!(saturated_transformer_intake)
    if hasproperty(saturated_transformer_intake, :validation)
        assert_valid!(getproperty(saturated_transformer_intake, :validation))
    end
    return saturated_transformer_intake
end

function _saturated_transformer_branch_arrays(saturated_transformer_intake)
    _assert_saturated_transformer_intake_valid!(saturated_transformer_intake)
    if hasproperty(saturated_transformer_intake, :breakpoints) &&
       isempty(getproperty(saturated_transformer_intake, :breakpoints))
        return _saturated_transformer_linear_branch_arrays(saturated_transformer_intake)
    end
    return saturated_transformer_nonlinear_arrays(saturated_transformer_intake)
end

function saturated_transformer_winding_node_map(saturated_transformer_intake)
    return saturated_transformer_winding_node_map(
        _saturated_transformer_branch_arrays(saturated_transformer_intake),
    )
end

function saturated_transformer_branch_augmented_step_context(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_intake;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    winding_number::Int = 1,
    include_coupled_lumped_sequence_history::Bool = false,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    transformer_branch_shunt_capacitance_rows = nothing,
    recorded_step_indices = nothing,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    DeckParser.assert_deck_valid!(parsed)
    arrays = _saturated_transformer_branch_arrays(saturated_transformer_intake)
    physical_node_map = _saturated_transformer_physical_node_map(parsed, arrays)
    frequency_partition =
        DeckParser.deck_steady_state_frequency_partition(parsed)
    transformer_frequency_hz = _steady_state_terminal_frequency_hz(
        frequency_partition,
        _saturated_transformer_frequency_node_indices(
            physical_node_map,
            arrays,
            length(frequency_partition.node_frequencies_hz),
        ),
        _deck_steady_state_frequency_hz(parsed),
    )
    branch_assembly = saturated_transformer_winding_branch_assembly(
        arrays,
        physical_node_map;
        nonlinear_winding_number = winding_number,
    )
    branch_elements = saturated_transformer_branch_elements(
        branch_assembly;
        reactance_units = 2.0 * pi * transformer_frequency_hz,
    )
    node_map = _saturated_transformer_augmented_node_map(physical_node_map, branch_assembly)
    deck_elements = Any[parsed.elements...]
    deck_element_names = copy(parsed.element_names)
    _append_switching_nonlinear_resistor_safety_shunts!(
        deck_elements,
        deck_element_names,
        parsed,
    )
    _delay_deck_time_switch_events!(
        deck_elements,
        time_switch_event_delay_s,
        t_end_s,
    )
    current_zero_switching &&
        _convert_deck_current_zero_switches!(deck_elements, deck_element_names, parsed)
    elements = Any[deck_elements...; branch_elements.elements...]
    element_names = Symbol[deck_element_names...; branch_elements.element_names...]
    _append_transformer_branch_shunt_capacitance_elements!(
        elements,
        element_names,
        node_map,
        parsed,
        transformer_branch_shunt_capacitance_rows === nothing ?
        DeckParser.DeckTransformerBranchShuntCapacitanceRow[] :
        transformer_branch_shunt_capacitance_rows,
    )
    if include_coupled_lumped_sequence_history
        source_equivalent = coupled_lumped_sequence_history_injection_elements(parsed)
        append!(elements, source_equivalent.elements)
        append!(element_names, source_equivalent.element_names)
    end
    source_function_runtime, control_system_runtime =
        _append_dynamic_source_and_control_elements!(
        elements,
        element_names,
        parsed,
        dt_s,
        source_signal_provider,
    )
    system = NodalSystem(maximum(values(node_map); init = 0), elements)
    return initialize_step_context(
        system;
        node_map = node_map,
        element_names = element_names,
        source_function_runtime = source_function_runtime,
        control_system_runtime = control_system_runtime,
        _deck_runtime_output_context_kwargs(
            parsed;
            time_switch_event_delay_s = time_switch_event_delay_s,
            event_horizon_s = t_end_s,
        )...,
        dt_s = dt_s,
        t_end_s = t_end_s,
        source = parsed.source,
        recorded_step_indices = recorded_step_indices,
    )
end

function _saturated_transformer_runtime_node_count(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_intake,
)
    arrays = _saturated_transformer_branch_arrays(saturated_transformer_intake)
    physical_node_map = _saturated_transformer_physical_node_map(parsed, arrays)
    branch_assembly = saturated_transformer_winding_branch_assembly(
        arrays,
        physical_node_map;
        nonlinear_winding_number = 1,
    )
    node_map = _saturated_transformer_augmented_node_map(
        physical_node_map,
        branch_assembly,
    )
    return maximum(values(node_map); init = 0)
end

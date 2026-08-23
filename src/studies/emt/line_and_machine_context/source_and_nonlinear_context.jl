
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
    inactive_result,
)
    source_config === nothing && inactive_result !== nothing &&
        return inactive_result
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
        configured = getproperty(raw_over16_kwargs, :current_injection_values)
        length(configured) >= node_count ||
            throw(ArgumentError("current_injection_values must cover every node"))
        if haskey(raw_over16_kwargs, :current_injection_workspace)
            workspace = getproperty(
                raw_over16_kwargs,
                :current_injection_workspace,
            )
            workspace isa Vector{Float64} || throw(ArgumentError(
                "current_injection_workspace must be a Vector{Float64}",
            ))
            length(workspace) == node_count || throw(ArgumentError(
                "current_injection_workspace must match the EMT node count",
            ))
            Base.mightalias(workspace, configured) && throw(ArgumentError(
                "current_injection_workspace must not alias current_injection_values",
            ))
            copyto!(workspace, 1, configured, 1, node_count)
            return workspace
        end
        return Float64[configured[index] for index in 1:node_count]
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

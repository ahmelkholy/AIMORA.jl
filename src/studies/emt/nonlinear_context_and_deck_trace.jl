
function _deck_offset_nonlinear_table_index(index::Integer, offset::Int)
    value = Int(index)
    value == 0 && return 0
    return value > 0 ? value + offset : -(abs(value) + offset)
end

function _deck_nonlinear_current_row_is_simultaneous(
    config::NamedTuple,
    row_index::Int,
)
    subsystem_indices = Int.(get(config, :nonlinear_subsystem_indices, Int[]))
    row_index <= length(subsystem_indices) || return false
    subsystem_index = subsystem_indices[row_index]
    begin_indices = Int.(get(config, :subsystem_begin_indices, Int[]))
    1 <= subsystem_index <= length(begin_indices) || return false
    head = begin_indices[subsystem_index]
    flags = Int.(get(config, :subsystem_simultaneous_flags, Int[]))
    return 1 <= head <= length(flags) && flags[head] != 0
end

function _deck_nonlinear_current_subnetwork_element_type(
    config::NamedTuple,
    row_index::Int,
    nonlinear_type::Int,
)
    row_indices = Int.(get(config, :subnetwork_nonlinear_indices, Int[]))
    element_types = Int.(get(config, :subnetwork_element_types, Int[]))
    slot = findfirst(==(row_index), row_indices)
    if slot !== nothing && slot <= length(element_types) && element_types[slot] != 0
        return element_types[slot]
    end
    nonlinear_type == 921 && return 1
    nonlinear_type == 922 && return 2
    nonlinear_type == 923 && return 3
    return 0
end

function _deck_append_unique_symbols!(dest::Vector{Symbol}, values)
    for value in values
        symbol = Symbol(value)
        symbol in dest || push!(dest, symbol)
    end
    return dest
end

function _deck_append_unique_ints!(dest::Vector{Int}, values)
    for value in values
        int_value = Int(value)
        int_value in dest || push!(dest, int_value)
    end
    return dest
end

function _deck_merge_nonlinear_current_configs(configs::Vector{NamedTuple})
    isempty(configs) && return nothing
    complete_source_loop = all(
        config -> get(config, :complete_nonlinear_source_loop, false),
        configs,
    )

    delta = Float64(get(first(configs), :delta2, 1.0))
    isfinite(delta) && delta > 0.0 ||
        throw(ArgumentError("delta2 must be finite and positive"))
    epszno = Float64(get(first(configs), :epszno, 1.0e-10))
    znolim = get(first(configs), :znolim, (0.2, 1.5))
    max_iterations = Int(get(first(configs), :max_iterations, 20))
    fltinf = Float64(get(first(configs), :fltinf, Inf))
    flzero = Float64(get(first(configs), :flzero, 0.0))
    time = Float64(get(first(configs), :t, 0.0))
    deltat = Float64(get(first(configs), :deltat, 2.0 * delta))

    runtime_node_indices = Int[]
    reference_node_index = nothing
    required_node_count = 1
    source_activity_flags = Int[]
    nonlinear_types = Int[]
    nonlinear_from_nodes = Int[]
    nonlinear_to_nodes = Int[]
    nonlinear_admittance_nodes = Int[]
    nonlinear_table_end_indices = Int[]
    nonlinear_subsystem_flags = Bool[]
    nonlinear_element_types = Int[]
    nonlinear_deck_from_nodes = Int[]
    nonlinear_deck_to_nodes = Int[]
    nonlinear_deck_from_node_names = Symbol[]
    nonlinear_deck_to_node_names = Symbol[]
    initial_companion_current_values = Float64[]
    initial_characteristic_current_values = Float64[]
    initial_stored_voltage_values = Float64[]
    initial_runtime_voltage_values = Float64[]
    initial_current_segment_values = Float64[]
    initial_table_index_values = Int[]
    minimum_on_time_values = Float64[]
    timed_resistance_arm_time_values = Float64[]
    rearm_time_state_indices = Int[]
    single_flash_flags = Bool[]
    nonlinear_output_codes = Int[]
    nonlinear_owner_names = Symbol[]
    nonlinear_owner_line_numbers = Int[]
    safety_shunt_resistor_names = Symbol[]
    safety_shunt_resistance_ohm_values = Float64[]
    safety_shunt_conductance_s_values = Float64[]
    nonlinear_current_segments = Float64[]
    nonlinear_steady_state_current_values = Float64[]
    nonlinear_steady_state_flux_values = Float64[]
    fortran_nonlinear_admittance_nodes = Int[]
    fortran_nonlinear_state_start_indices = Int[]
    fortran_gap_status_values = Float64[]
    cchar = Float64[]
    vchar = Float64[]
    gslope = Float64[]
    fortran_files = Symbol[]
    mutation_order = Symbol[]
    deferred_calls = Symbol[]
    fortran_labels = Int[]
    config_sources = Symbol[]
    saturated_transformer_sparse_config = nothing
    saturated_transformer_sparse_from_nodes = Int[]
    saturated_transformer_sparse_to_nodes = Int[]
    saturated_transformer_internal_top_node_indices = Int[]
    saturated_transformer_internal_top_voltage_source_nodes = Int[]
    saturated_transformer_terminal_node_indices = Int[]
    saw_arrester_state_rows = false

    for config in configs
        get(config, :deck_owned_voltage_context, false) ||
            throw(ArgumentError("mixed nonlinear current configs must be deck-owned"))
        get(config, :reference_shifted_voltage_context, false) ||
            throw(ArgumentError("mixed nonlinear current configs must use the reference-shifted voltage context"))
        isapprox(Float64(get(config, :delta2, delta)), delta; rtol = 0.0, atol = 0.0) ||
            throw(ArgumentError("mixed nonlinear current configs must use one delta2"))
        haskey(config, :simultaneous_zno_config) &&
            throw(ArgumentError("explicit simultaneous nonlinear inverse columns must be supplied after mixed config composition"))

        mapping = Int.(get(config, :deck_to_runtime_node_indices, Int[]))
        if !isempty(mapping)
            if isempty(runtime_node_indices)
                runtime_node_indices = copy(mapping)
            elseif mapping != runtime_node_indices
                throw(ArgumentError("mixed nonlinear current configs must share one deck node mapping"))
            end
        end
        local_reference = Int(get(config, :reference_node_index, 1))
        if reference_node_index === nothing
            reference_node_index = local_reference
        elseif reference_node_index != local_reference
            throw(ArgumentError("mixed nonlinear current configs must share one reference node"))
        end

        local_types = Int.(get(config, :nonlinear_types, Int[]))
        local_count = length(local_types)
        local_from_nodes = Int.(get(config, :nonlinear_from_nodes, Int[]))
        local_to_nodes = Int.(get(config, :nonlinear_to_nodes, Int[]))
        local_admittance_nodes = Int.(get(config, :nonlinear_admittance_nodes, Int[]))
        local_table_end_indices = Int.(get(config, :nonlinear_table_end_indices, Int[]))
        length(local_from_nodes) == local_count &&
            length(local_to_nodes) == local_count &&
            length(local_admittance_nodes) == local_count &&
            length(local_table_end_indices) == local_count ||
            throw(ArgumentError("mixed nonlinear current config row arrays must match nonlinear_types"))

        has_arrester_rows = any(==(94), local_types)
        has_table_rows = any(!=(94), local_types)
        saw_arrester_state_rows && has_table_rows &&
            throw(ArgumentError("simultaneous nonlinear table rows must precede arrester state rows"))

        local_cchar = Float64.(get(config, :cchar, Float64[]))
        local_vchar = Float64.(get(config, :vchar, Float64[]))
        local_gslope = Float64.(get(config, :gslope, Float64[]))
        cchar_offset = length(cchar)
        vchar_offset = length(vchar)

        local_initial_anonl =
            Float64.(get(config, :initial_companion_current_values, zeros(Float64, local_count)))
        local_initial_vzero =
            Float64.(get(config, :initial_characteristic_current_values, zeros(Float64, local_count)))
        local_initial_vnonl =
            Float64.(get(config, :initial_stored_voltage_values, zeros(Float64, local_count)))
        local_initial_runtime_vnonl =
            Float64.(get(config, :initial_runtime_voltage_values, local_initial_vnonl))
        local_initial_curr =
            Float64.(get(config, :initial_current_segment_values, zeros(Float64, local_count)))
        local_initial_ilast =
            Int.(get(config, :initial_table_index_values, zeros(Int, local_count)))
        local_segments =
            Float64.(get(config, :nonlinear_current_segments, zeros(Float64, local_count)))
        local_steady_currents = Float64.(get(
            config,
            :nonlinear_steady_state_current_values,
            zeros(Float64, local_count),
        ))
        local_steady_fluxes = Float64.(get(
            config,
            :nonlinear_steady_state_flux_values,
            zeros(Float64, local_count),
        ))
        local_minimum_on_time =
            Float64.(get(config, :minimum_on_time_values, zeros(Float64, local_count)))
        local_timed_resistance_arm_times = Float64.(get(
            config,
            :timed_resistance_arm_time_values,
            zeros(Float64, local_count),
        ))
        local_rearm_indices =
            Int.(get(config, :rearm_time_state_indices, zeros(Int, local_count)))
        local_single_flash =
            Bool.(get(config, :single_flash_flags, falses(local_count)))
        local_output_codes =
            Int.(get(config, :nonlinear_output_codes, zeros(Int, local_count)))
        local_owner_names =
            Symbol.(get(config, :nonlinear_owner_names, fill(Symbol(""), local_count)))
        local_owner_line_numbers =
            Int.(get(config, :nonlinear_owner_line_numbers, zeros(Int, local_count)))
        length(local_initial_anonl) == local_count &&
            length(local_initial_vzero) == local_count &&
            length(local_initial_vnonl) == local_count &&
            length(local_initial_runtime_vnonl) == local_count &&
            length(local_initial_curr) == local_count &&
            length(local_initial_ilast) == local_count &&
            length(local_segments) == local_count &&
            length(local_steady_currents) == local_count &&
            length(local_steady_fluxes) == local_count &&
            length(local_minimum_on_time) == local_count &&
            length(local_timed_resistance_arm_times) == local_count &&
            length(local_rearm_indices) == local_count &&
            length(local_single_flash) == local_count &&
            length(local_output_codes) == local_count &&
            length(local_owner_names) == local_count &&
            length(local_owner_line_numbers) == local_count ||
            throw(ArgumentError("mixed nonlinear current initial arrays must match nonlinear_types"))

        local_fortran_nonlad =
            Int.(get(config, :fortran_nonlinear_admittance_nodes, local_admittance_nodes))
        local_fortran_state =
            Int.(get(config, :fortran_nonlinear_state_start_indices, zeros(Int, local_count)))
        local_gap_status =
            Float64.(get(config, :fortran_gap_status_values, zeros(Float64, local_count)))
        length(local_fortran_nonlad) == local_count &&
            length(local_fortran_state) == local_count &&
            length(local_gap_status) == local_count ||
            throw(ArgumentError("mixed nonlinear current Fortran-side arrays must match nonlinear_types"))

        append!(cchar, local_cchar)
        append!(vchar, local_vchar)
        if length(local_gslope) >= length(local_cchar)
            append!(gslope, local_gslope[1:length(local_cchar)])
        else
            append!(gslope, local_gslope)
            append!(gslope, zeros(Float64, length(local_cchar) - length(local_gslope)))
        end

        required_node_count = max(
            required_node_count,
            Int(get(config, :nonlinear_required_node_count, 1)),
            maximum(vcat(local_from_nodes, local_to_nodes); init = 1),
        )
        local_activity = Int.(get(config, :nonlinear_source_activity_flags, Int[]))
        if length(source_activity_flags) < required_node_count
            old_count = length(source_activity_flags)
            resize!(source_activity_flags, required_node_count)
            for index in (old_count + 1):required_node_count
                source_activity_flags[index] = 0
            end
        end
        for index in eachindex(local_activity)
            index <= required_node_count || break
            local_activity[index] != 0 && (source_activity_flags[index] = 1)
        end
        for node in vcat(local_from_nodes, local_to_nodes)
            node > 1 && (source_activity_flags[node] = 1)
        end

        local_deck_from_nodes =
            Int.(get(config, :nonlinear_deck_from_nodes, zeros(Int, local_count)))
        local_deck_to_nodes =
            Int.(get(config, :nonlinear_deck_to_nodes, zeros(Int, local_count)))
        local_deck_from_names =
            Symbol.(get(config, :nonlinear_deck_from_node_names, fill(Symbol(""), local_count)))
        local_deck_to_names =
            Symbol.(get(config, :nonlinear_deck_to_node_names, fill(Symbol(""), local_count)))
        length(local_deck_from_nodes) == local_count &&
            length(local_deck_to_nodes) == local_count &&
            length(local_deck_from_names) == local_count &&
            length(local_deck_to_names) == local_count ||
            throw(ArgumentError("mixed nonlinear current deck-node arrays must match nonlinear_types"))
        local_sparse_config = get(config, :saturated_transformer_sparse_config, nothing)
        local_sparse_from_nodes = Int[]
        local_sparse_to_nodes = Int[]
        if local_sparse_config !== nothing
            local_sparse_config isa NamedTuple ||
                throw(ArgumentError("saturated transformer sparse config must be a NamedTuple"))
            saturated_transformer_sparse_config === nothing ||
                throw(ArgumentError("mixed nonlinear current configs must contain one saturated transformer sparse config"))
            saturated_transformer_sparse_config = local_sparse_config
            local_sparse_from_nodes = Int.(get(local_sparse_config, :from_nodes, Int[]))
            local_sparse_to_nodes = Int.(get(local_sparse_config, :to_nodes, Int[]))
            length(local_sparse_from_nodes) == local_count &&
                length(local_sparse_to_nodes) == local_count ||
                throw(ArgumentError("saturated transformer sparse endpoint arrays must match local nonlinear rows"))
        end
        append!(
            saturated_transformer_internal_top_node_indices,
            Int.(get(config, :saturated_transformer_internal_top_node_indices, Int[])),
        )
        append!(
            saturated_transformer_internal_top_voltage_source_nodes,
            Int.(get(
                config,
                :saturated_transformer_internal_top_voltage_source_nodes,
                Int[],
            )),
        )
        append!(
            saturated_transformer_terminal_node_indices,
            Int.(get(config, :saturated_transformer_terminal_node_indices, Int[])),
        )

        for local_index in 1:local_count
            nonlinear_type = local_types[local_index]
            table_offset = nonlinear_type == 94 ? vchar_offset : cchar_offset
            push!(nonlinear_types, nonlinear_type)
            push!(nonlinear_from_nodes, local_from_nodes[local_index])
            push!(nonlinear_to_nodes, local_to_nodes[local_index])
            push!(
                nonlinear_admittance_nodes,
                _deck_offset_nonlinear_table_index(
                    local_admittance_nodes[local_index],
                    cchar_offset,
                ),
            )
            push!(
                nonlinear_table_end_indices,
                _deck_offset_nonlinear_table_index(
                    local_table_end_indices[local_index],
                    table_offset,
                ),
            )
            push!(
                nonlinear_subsystem_flags,
                _deck_nonlinear_current_row_is_simultaneous(config, local_index),
            )
            push!(
                nonlinear_element_types,
                _deck_nonlinear_current_subnetwork_element_type(
                    config,
                    local_index,
                    nonlinear_type,
                ),
            )
            push!(initial_companion_current_values, local_initial_anonl[local_index])
            push!(initial_characteristic_current_values, local_initial_vzero[local_index])
            push!(initial_stored_voltage_values, local_initial_vnonl[local_index])
            push!(initial_runtime_voltage_values, local_initial_runtime_vnonl[local_index])
            push!(initial_current_segment_values, local_initial_curr[local_index])
            push!(
                initial_table_index_values,
                _deck_offset_nonlinear_table_index(
                    local_initial_ilast[local_index],
                    table_offset,
                ),
            )
            push!(nonlinear_current_segments, local_segments[local_index])
            push!(
                nonlinear_steady_state_current_values,
                local_steady_currents[local_index],
            )
            push!(
                nonlinear_steady_state_flux_values,
                local_steady_fluxes[local_index],
            )
            push!(minimum_on_time_values, local_minimum_on_time[local_index])
            push!(
                timed_resistance_arm_time_values,
                local_timed_resistance_arm_times[local_index],
            )
            push!(
                rearm_time_state_indices,
                _deck_offset_nonlinear_table_index(
                    local_rearm_indices[local_index],
                    cchar_offset,
                ),
            )
            push!(single_flash_flags, local_single_flash[local_index])
            push!(nonlinear_output_codes, local_output_codes[local_index])
            push!(nonlinear_owner_names, local_owner_names[local_index])
            push!(nonlinear_owner_line_numbers, local_owner_line_numbers[local_index])
            push!(
                fortran_nonlinear_admittance_nodes,
                _deck_offset_nonlinear_table_index(
                    local_fortran_nonlad[local_index],
                    cchar_offset,
                ),
            )
            push!(
                fortran_nonlinear_state_start_indices,
                _deck_offset_nonlinear_table_index(
                    local_fortran_state[local_index],
                    vchar_offset,
                ),
            )
            push!(fortran_gap_status_values, local_gap_status[local_index])
            push!(nonlinear_deck_from_nodes, local_deck_from_nodes[local_index])
            push!(nonlinear_deck_to_nodes, local_deck_to_nodes[local_index])
            push!(nonlinear_deck_from_node_names, local_deck_from_names[local_index])
            push!(nonlinear_deck_to_node_names, local_deck_to_names[local_index])
            if nonlinear_type == SATURATED_TRANSFORMER_NONLINEAR_TYPE &&
               local_sparse_config !== nothing
                push!(saturated_transformer_sparse_from_nodes, local_sparse_from_nodes[local_index])
                push!(saturated_transformer_sparse_to_nodes, local_sparse_to_nodes[local_index])
            else
                push!(saturated_transformer_sparse_from_nodes, max(local_from_nodes[local_index], 1))
                push!(saturated_transformer_sparse_to_nodes, max(local_to_nodes[local_index], 1))
            end
        end

        _deck_append_unique_symbols!(fortran_files, get(config, :fortran_files, ()))
        _deck_append_unique_symbols!(mutation_order, get(config, :mutation_order, ()))
        _deck_append_unique_symbols!(deferred_calls, get(config, :deferred_calls, ()))
        _deck_append_unique_ints!(fortran_labels, get(config, :fortran_labels, ()))
        push!(config_sources, Symbol(get(config, :source, :deck_nonlinear_current_config)))
        append!(
            safety_shunt_resistor_names,
            Symbol.(get(config, :safety_shunt_resistor_names, Symbol[])),
        )
        append!(
            safety_shunt_resistance_ohm_values,
            Float64.(get(config, :safety_shunt_resistance_ohm_values, Float64[])),
        )
        append!(
            safety_shunt_conductance_s_values,
            Float64.(get(config, :safety_shunt_conductance_s_values, Float64[])),
        )
        has_arrester_rows && (saw_arrester_state_rows = true)
    end

    total_count = length(nonlinear_types)
    total_count > 0 || throw(ArgumentError("mixed nonlinear current config is empty"))
    runtime_node_count = max(required_node_count - 1, 0)
    if length(runtime_node_indices) < runtime_node_count
        first_missing = length(runtime_node_indices) + 1
        append!(
            runtime_node_indices,
            (index + 1 for index in first_missing:runtime_node_count),
        )
    end
    if length(source_activity_flags) < required_node_count
        old_count = length(source_activity_flags)
        resize!(source_activity_flags, required_node_count)
        for index in (old_count + 1):required_node_count
            source_activity_flags[index] = 0
        end
    end

    source_last_slot = 1 + 5 * (total_count - 1)
    source_next_indices = zeros(Int, source_last_slot)
    source_from_nodes = zeros(Int, source_last_slot)
    source_to_nodes = zeros(Int, source_last_slot)
    for index in 1:total_count
        slot = 1 + 5 * (index - 1)
        source_next_indices[slot] = index == total_count ? 1 : slot + 5
        source_from_nodes[slot] = nonlinear_from_nodes[index]
        source_to_nodes[slot] = nonlinear_to_nodes[index]
    end

    simultaneous_rows = findall(identity, nonlinear_subsystem_flags)
    subnetwork_last_slot = isempty(simultaneous_rows) ? 0 : 1 + 5 * (length(simultaneous_rows) - 1)
    subnetwork_next_indices = zeros(Int, subnetwork_last_slot)
    subnetwork_from_nodes = zeros(Int, subnetwork_last_slot)
    subnetwork_to_nodes = zeros(Int, subnetwork_last_slot)
    subnetwork_nonlinear_indices = zeros(Int, subnetwork_last_slot)
    subnetwork_element_types = zeros(Int, subnetwork_last_slot)
    for (position, nonlinear_index) in enumerate(simultaneous_rows)
        slot = 1 + 5 * (position - 1)
        subnetwork_next_indices[slot] =
            position == length(simultaneous_rows) ? 1 : slot + 5
        subnetwork_from_nodes[slot] = nonlinear_from_nodes[nonlinear_index]
        subnetwork_to_nodes[slot] = nonlinear_to_nodes[nonlinear_index]
        subnetwork_nonlinear_indices[slot] = nonlinear_index
        subnetwork_element_types[slot] = nonlinear_element_types[nonlinear_index]
    end

    nonlinear_subsystem_indices = zeros(Int, total_count)
    subsystem_begin_indices = Int[]
    subsystem_owner_rows = zeros(Int, source_last_slot)
    subsystem_simultaneous_flags = zeros(Int, source_last_slot)
    if !isempty(simultaneous_rows)
        push!(subsystem_begin_indices, 1)
        subsystem_owner_rows[1] = first(simultaneous_rows)
        subsystem_simultaneous_flags[1] = 1
        for nonlinear_index in simultaneous_rows
            nonlinear_subsystem_indices[nonlinear_index] = 1
        end
    end
    for nonlinear_index in 1:total_count
        nonlinear_subsystem_indices[nonlinear_index] != 0 && continue
        slot = 1 + 5 * (nonlinear_index - 1)
        push!(subsystem_begin_indices, slot)
        subsystem_index = length(subsystem_begin_indices)
        nonlinear_subsystem_indices[nonlinear_index] = subsystem_index
        subsystem_owner_rows[slot] = nonlinear_index
    end

    initial_cursub_count = max(required_node_count, div(source_last_slot, 5) + 1)
    partition_boundary = max(required_node_count - 1, 1)
    source_begin_indices = [1]
    reference_node = reference_node_index === nothing ? 1 : Int(reference_node_index)
    merged_config = (
        source = :deck_nonlinear_current_config,
        outcome = :nonlinear_current_config,
        component_sources = Tuple(config_sources),
        nonlinear_types = nonlinear_types,
        nonlinear_from_nodes = nonlinear_from_nodes,
        nonlinear_to_nodes = nonlinear_to_nodes,
        nonlinear_admittance_nodes = nonlinear_admittance_nodes,
        nonlinear_table_end_indices = nonlinear_table_end_indices,
        nonlinear_subsystem_indices = nonlinear_subsystem_indices,
        subsystem_begin_indices = subsystem_begin_indices,
        subsystem_owner_rows = subsystem_owner_rows,
        subsystem_simultaneous_flags = subsystem_simultaneous_flags,
        cchar = cchar,
        vchar = vchar,
        gslope = gslope,
        delta2 = delta,
        epszno = epszno,
        znolim = znolim,
        max_iterations = max_iterations,
        ncomp = total_count,
        fltinf = fltinf,
        flzero = flzero,
        t = time,
        deltat = deltat,
        deck_owned_voltage_context = true,
        reference_shifted_voltage_context = true,
        deck_to_runtime_node_indices = runtime_node_indices,
        reference_node_index = reference_node,
        nonlinear_required_node_count = required_node_count,
        initialize_nonlinear_state = true,
        initial_companion_current_values = initial_companion_current_values,
        initial_characteristic_current_values = initial_characteristic_current_values,
        initial_stored_voltage_values = initial_stored_voltage_values,
        initial_runtime_voltage_values = initial_runtime_voltage_values,
        initial_current_segment_values = initial_current_segment_values,
        initial_table_index_values = initial_table_index_values,
        initial_cursub_values = zeros(Float64, initial_cursub_count),
        nonlinear_current_segments = nonlinear_current_segments,
        nonlinear_steady_state_current_values =
            nonlinear_steady_state_current_values,
        nonlinear_steady_state_flux_values = nonlinear_steady_state_flux_values,
        saturated_transformer_residual_flux_initialized = any(
            config -> get(
                config,
                :saturated_transformer_residual_flux_initialized,
                false,
            ),
            configs,
        ),
        minimum_on_time_values = minimum_on_time_values,
        timed_resistance_arm_time_values = timed_resistance_arm_time_values,
        rearm_time_state_indices = rearm_time_state_indices,
        single_flash_flags = single_flash_flags,
        nonlinear_output_codes = nonlinear_output_codes,
        nonlinear_owner_names = nonlinear_owner_names,
        nonlinear_owner_line_numbers = nonlinear_owner_line_numbers,
        safety_shunt_resistor_names = safety_shunt_resistor_names,
        safety_shunt_resistance_ohm_values = safety_shunt_resistance_ohm_values,
        safety_shunt_conductance_s_values = safety_shunt_conductance_s_values,
        dense_primary_nonlinear_compensation = any(
            config -> get(config, :dense_primary_nonlinear_compensation, false),
            configs,
        ),
        use_state_nonlinear_inverse_columns = any(
            config -> get(config, :use_state_nonlinear_inverse_columns, false),
            configs,
        ),
        nonlinear_inverse_config = (
            require_fortran_sparse_factor_result = true,
            ntot = required_node_count,
            ncomp = total_count,
            source_begin_indices = source_begin_indices,
            source_next_indices = source_next_indices,
            source_from_nodes = source_from_nodes,
            source_to_nodes = source_to_nodes,
            source_activity_flags = source_activity_flags,
            nonlinear_types = nonlinear_types,
            nonlinear_admittance_nodes = nonlinear_admittance_nodes,
            nonlinear_from_nodes = nonlinear_from_nodes,
            nonlinear_to_nodes = nonlinear_to_nodes,
            nonlinear_source_flags = zeros(Int, total_count),
            partition_boundary = partition_boundary,
            delta2 = delta,
            fltinf = fltinf,
        ),
        nonlinear_source_begin_indices = source_begin_indices,
        nonlinear_source_next_indices = source_next_indices,
        nonlinear_source_from_nodes = source_from_nodes,
        nonlinear_source_to_nodes = source_to_nodes,
        nonlinear_source_activity_flags = source_activity_flags,
        nonlinear_source_flags = zeros(Int, total_count),
        nonlinear_inverse_partition_boundary = partition_boundary,
        nonlinear_deck_from_nodes = nonlinear_deck_from_nodes,
        nonlinear_deck_to_nodes = nonlinear_deck_to_nodes,
        nonlinear_deck_from_node_names = nonlinear_deck_from_node_names,
        nonlinear_deck_to_node_names = nonlinear_deck_to_node_names,
        fortran_nonlinear_admittance_nodes = fortran_nonlinear_admittance_nodes,
        fortran_nonlinear_state_start_indices = fortran_nonlinear_state_start_indices,
        fortran_gap_status_values = fortran_gap_status_values,
        subnetwork_next_indices = subnetwork_next_indices,
        subnetwork_from_nodes = subnetwork_from_nodes,
        subnetwork_to_nodes = subnetwork_to_nodes,
        subnetwork_nonlinear_indices = subnetwork_nonlinear_indices,
        subnetwork_element_types = subnetwork_element_types,
        table_entry_count = length(cchar),
        state_entry_count = length(vchar),
        nonlinear_row_count = total_count,
        fortran_files = Tuple(fortran_files),
        fortran_labels = Tuple(fortran_labels),
        mutation_order = Tuple(vcat([:mixed_nonlinear_current_config], mutation_order)),
        deferred_calls = Tuple(deferred_calls),
        complete_nonlinear_source_loop = complete_source_loop,
        replacement_ready = complete_source_loop && isempty(deferred_calls),
    )
    if !isempty(saturated_transformer_internal_top_node_indices)
        merged_config = merge(
            merged_config,
            (
                saturated_transformer_internal_top_node_indices =
                    saturated_transformer_internal_top_node_indices,
                saturated_transformer_internal_top_voltage_source_nodes =
                    saturated_transformer_internal_top_voltage_source_nodes,
                saturated_transformer_terminal_node_indices =
                    saturated_transformer_terminal_node_indices,
                saturated_transformer_required_node_count = required_node_count,
            ),
        )
    end
    saturated_transformer_sparse_config !== nothing && return merge(
        merged_config,
        (
            saturated_transformer_sparse_config = merge(
                saturated_transformer_sparse_config,
                (
                    from_nodes = saturated_transformer_sparse_from_nodes,
                    to_nodes = saturated_transformer_sparse_to_nodes,
                ),
            ),
        ),
    )
    return merged_config
end

function _deck_saturated_transformer_merge_current_config(
    config::NamedTuple,
    deck_to_runtime_node_indices::AbstractVector{<:Integer},
)
    nonlinear_types = Int.(get(config, :nonlinear_types, Int[]))
    count = length(nonlinear_types)
    source_from_nodes = Int.(get(config, :nonlinear_from_nodes, Int[]))
    source_to_nodes = Int.(get(config, :nonlinear_to_nodes, Int[]))
    length(source_from_nodes) == count && length(source_to_nodes) == count ||
        throw(ArgumentError("saturated transformer nonlinear endpoint arrays must match nonlinear_types"))
    shifted_from_nodes = [
        _deck_saturated_transformer_runtime_node(node, deck_to_runtime_node_indices)
        for node in source_from_nodes
    ]
    shifted_to_nodes = [
        _deck_saturated_transformer_runtime_node(node, deck_to_runtime_node_indices)
        for node in source_to_nodes
    ]
    shifted_internal_nodes = Int[
        _deck_saturated_transformer_runtime_node(node, deck_to_runtime_node_indices)
        for node in Int.(get(config, :saturated_transformer_internal_top_node_indices, Int[]))
    ]
    shifted_source_nodes = Int[
        _deck_saturated_transformer_runtime_node(node, deck_to_runtime_node_indices)
        for node in Int.(get(
            config,
            :saturated_transformer_internal_top_voltage_source_nodes,
            Int[],
        ))
    ]
    shifted_terminal_nodes = Int[
        _deck_saturated_transformer_runtime_node(node, deck_to_runtime_node_indices)
        for node in Int.(get(config, :saturated_transformer_terminal_node_indices, Int[]))
    ]
    shifted_required_node_count = maximum(
        vcat(
            shifted_from_nodes,
            abs.(shifted_to_nodes),
            shifted_internal_nodes,
            shifted_source_nodes,
            shifted_terminal_nodes,
        );
        init = 1,
    )
    required_node_count = max(
        Int(get(config, :nonlinear_required_node_count, 1)),
        Int(get(config, :saturated_transformer_required_node_count, 1)),
        shifted_required_node_count,
    )
    current_segments =
        Float64.(get(config, :nonlinear_current_segments, ones(Float64, count)))
    table_start_indices = Int.(get(config, :nonlinear_admittance_nodes, ones(Int, count)))
    table_end_indices = Int.(get(config, :nonlinear_table_end_indices, table_start_indices))
    length(current_segments) == count &&
        length(table_start_indices) == count &&
        length(table_end_indices) == count ||
        throw(ArgumentError("saturated transformer nonlinear table arrays must match nonlinear_types"))
    initial_table_indices = Int[]
    for index in 1:count
        segment = max(Int(abs(current_segments[index])), 1)
        table_index = table_start_indices[index] + segment - 1
        table_start_indices[index] <= table_index <= table_end_indices[index] ||
            throw(ArgumentError("saturated transformer current seed must address its table"))
        push!(initial_table_indices, table_index)
    end
    initial_stored_voltages =
        Float64.(get(config, :initial_stored_voltage_values, zeros(Float64, count)))
    return merge(
        config,
        (
            nonlinear_from_nodes = shifted_from_nodes,
            nonlinear_to_nodes = shifted_to_nodes,
            nonlinear_top_nodes = shifted_internal_nodes,
            nonlinear_terminal_node_indices = shifted_terminal_nodes,
            saturated_transformer_internal_top_node_indices = shifted_internal_nodes,
            saturated_transformer_internal_top_voltage_source_nodes = shifted_source_nodes,
            saturated_transformer_terminal_node_indices = shifted_terminal_nodes,
            saturated_transformer_required_node_count = required_node_count,
            deck_owned_voltage_context = true,
            reference_shifted_voltage_context = true,
            reference_node_index = 1,
            nonlinear_required_node_count = required_node_count,
            initialize_nonlinear_state = true,
            initial_companion_current_values = Float64.(
                get(config, :initial_companion_current_values, zeros(Float64, count)),
            ),
            initial_characteristic_current_values = Float64.(
                get(
                    config,
                    :initial_characteristic_current_values,
                    get(config, :nonlinear_steady_state_current_values, zeros(Float64, count)),
                ),
            ),
            initial_stored_voltage_values = initial_stored_voltages,
            initial_runtime_voltage_values = Float64.(
                get(config, :initial_runtime_voltage_values, initial_stored_voltages),
            ),
            initial_current_segment_values = Float64.(
                get(config, :initial_current_segment_values, current_segments),
            ),
            initial_table_index_values = Int.(
                get(config, :initial_table_index_values, initial_table_indices),
            ),
            initial_cursub_values = Float64.(
                get(
                    config,
                    :initial_cursub_values,
                    zeros(Float64, max(required_node_count, count + 1)),
                ),
            ),
        ),
    )
end

function _deck_nonlinear_current_config(
    parsed::DeckParser.DeckParseResult,
    saturated_transformer_current_config;
    delta2::Real,
)
    deck_configs = NamedTuple[]
    saturated_transformer_config_index = 0
    if !isempty(DeckParser.deck_zinc_oxide_nonlinear_rows(parsed))
        push!(
            deck_configs,
            deck_zinc_oxide_nonlinear_current_config(parsed; delta2 = delta2),
        )
    end
    if !isempty(DeckParser.deck_nonlinear_resistance_rows(parsed))
        push!(
            deck_configs,
            deck_nonlinear_resistance_current_config(parsed; delta2 = delta2),
        )
    end
    if !isempty(DeckParser.deck_triggered_timed_resistance_rows(parsed))
        push!(
            deck_configs,
            deck_triggered_timed_resistance_current_config(parsed; delta2 = delta2),
        )
    end
    if !isempty(DeckParser.deck_switching_nonlinear_resistor_rows(parsed))
        push!(
            deck_configs,
            deck_switching_nonlinear_resistor_current_config(parsed; delta2 = delta2),
        )
    end
    if !isempty(DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed))
        push!(
            deck_configs,
            deck_piecewise_nonlinear_inductor_current_config(parsed; delta2 = delta2),
        )
    end
    if !isempty(DeckParser.deck_hysteretic_inductor_rows(parsed))
        push!(
            deck_configs,
            deck_hysteretic_inductor_nonlinear_current_config(parsed; delta2 = delta2),
        )
    end
    if saturated_transformer_current_config !== nothing
        push!(deck_configs, saturated_transformer_current_config)
        saturated_transformer_config_index = length(deck_configs)
    end
    if !isempty(DeckParser.deck_arrester_nonlinear_rows(parsed))
        push!(
            deck_configs,
            deck_arrester_nonlinear_current_config(parsed; delta2 = delta2),
        )
    end
    isempty(deck_configs) && return nothing
    length(deck_configs) == 1 && return first(deck_configs)
    if saturated_transformer_config_index != 0
        deck_to_runtime_node_indices = Int[]
        for config in deck_configs
            mapping = Int.(get(config, :deck_to_runtime_node_indices, Int[]))
            isempty(mapping) || (deck_to_runtime_node_indices = mapping; break)
        end
        deck_configs[saturated_transformer_config_index] =
            _deck_saturated_transformer_merge_current_config(
                deck_configs[saturated_transformer_config_index],
                deck_to_runtime_node_indices,
            )
    end
    return _deck_merge_nonlinear_current_configs(deck_configs)
end

const _REPORTED_NONLINEAR_BRANCH_TYPES = (921, 922, 923, 94, -96, -97)

function _nonlinear_branch_report(run, config::NamedTuple)
    types = Int.(get(config, :nonlinear_types, Int[]))
    indices = findall(type -> type in _REPORTED_NONLINEAR_BRANCH_TYPES, types)
    count = length(types)
    names = Symbol.(get(config, :nonlinear_owner_names, fill(Symbol(""), count)))
    line_numbers = Int.(get(config, :nonlinear_owner_line_numbers, zeros(Int, count)))
    output_codes = Int.(get(config, :nonlinear_output_codes, zeros(Int, count)))
    from_nodes = Int.(get(config, :nonlinear_from_nodes, Int[]))
    to_nodes = abs.(Int.(get(config, :nonlinear_to_nodes, Int[])))
    node_mapping = Int.(get(config, :deck_to_runtime_node_indices, Int[]))
    steps = length(run.over16_updates)
    times = zeros(Float64, steps)
    voltages = zeros(Float64, length(indices), steps)
    currents = zeros(Float64, length(indices), steps)
    fluxes = fill(NaN, length(indices), steps)
    segments = zeros(Int, length(indices), steps)
    for (step, update) in enumerate(run.over16_updates)
        times[step] = Float64(update.t_s)
        result = update.over16_update.nonlinear_current_result
        result === nothing && throw(ArgumentError(
            "nonlinear branch report requires a nonlinear timestep result",
        ))
        reference_voltage = zeros(Float64, Int(config.nonlinear_required_node_count))
        for deck_index in eachindex(node_mapping)
            reference_voltage[node_mapping[deck_index]] = update.voltage_pu[deck_index]
        end
        for (position, nonlinear_index) in enumerate(indices)
            voltage = reference_voltage[from_nodes[nonlinear_index]] -
                      reference_voltage[to_nodes[nonlinear_index]]
            nonlinear_type = types[nonlinear_index]
            voltages[position, step] = voltage
            currents[position, step] = if nonlinear_type == TRIGGERED_TIMED_RESISTANCE_TYPE
                segment = round(Int, result.curr[nonlinear_index])
                segment == 0 ? 0.0 :
                    result.gslope[
                        config.nonlinear_admittance_nodes[nonlinear_index] + segment - 1
                    ] * voltage
            elseif nonlinear_type == HYSTERETIC_INDUCTOR_NONLINEAR_TYPE
                state_start = config.nonlinear_admittance_nodes[nonlinear_index]
                result.anonl[nonlinear_index] + result.gslope[state_start + 1] * voltage
            else
                result.curr[nonlinear_index]
            end
            segments[position, step] = result.ilast[nonlinear_index]
            nonlinear_type == -96 &&
                (fluxes[position, step] = result.vnonl[nonlinear_index])
        end
    end
    return (
        source = :nonlinear_branch_report,
        outcome = :report_output,
        names = names[indices],
        line_numbers = line_numbers[indices],
        nonlinear_types = types[indices],
        output_codes = output_codes[indices],
        output_requested_flags = output_codes[indices] .> 0,
        time_s = times,
        voltage_v = voltages,
        current_a = currents,
        flux_wb = fluxes,
        active_segments = segments,
        deferred_effects = (),
    )
end

function _nonlinear_current_output_names(report, indices::AbstractVector{Int})
    names = Symbol[]
    used = Set{Symbol}()
    for index in indices
        owner = strip(String(report.names[index]))
        suffix = isempty(owner) ? "line_$(report.line_numbers[index])" : owner
        name = Symbol("nonlinear_current_", suffix)
        if name in used
            name = Symbol(name, "_line_", report.line_numbers[index])
        end
        name in used && (name = Symbol(name, "_row_", index))
        push!(names, name)
        push!(used, name)
    end
    return names
end

function _nonlinear_report_trace_columns(report, trace::DeckEMTTrace)
    columns = Vector{Int}(undef, length(trace.time_s))
    tolerance = max(abs(trace.dt_s) * 1.0e-12, 8.0 * eps(Float64))
    for (trace_index, time_s) in enumerate(trace.time_s)
        column = findfirst(
            report_time_s -> abs(report_time_s - time_s) <= tolerance,
            report.time_s,
        )
        column === nothing && throw(ArgumentError(
            "nonlinear report does not own the stored trace sample at t=$time_s s",
        ))
        columns[trace_index] = column
    end
    return columns
end

function _append_deck_nonlinear_outputs(trace::DeckEMTTrace, run)
    for report_name in (
        :nonlinear_branch_report,
        :switching_nonlinear_resistor_report,
        :piecewise_nonlinear_inductor_report,
    )
        hasproperty(run, report_name) || continue
        trace = _append_nonlinear_report_outputs(
            trace,
            getproperty(run, report_name),
        )
    end
    return trace
end

function _append_nonlinear_report_outputs(trace::DeckEMTTrace, report)
    requested_indices = findall(report.output_requested_flags)
    isempty(requested_indices) && return trace
    names = _nonlinear_current_output_names(report, requested_indices)
    existing = intersect(Set(trace.output_channel_names), Set(names))
    isempty(existing) || throw(ArgumentError(
        "nonlinear report output names duplicate existing trace channels: " *
        join(string.(sort!(collect(existing); by = String)), ','),
    ))
    report_columns = _nonlinear_report_trace_columns(report, trace)
    values = copy(report.current_a[requested_indices, report_columns])
    extrema = _sampled_trace_extrema(values, trace.time_s)
    return DeckEMTTrace(
        trace.source,
        trace.dt_s,
        trace.t_end_s,
        copy(trace.node_map),
        copy(trace.node_names),
        copy(trace.element_names),
        copy(trace.time_s),
        copy(trace.voltage_pu),
        vcat(copy(trace.output_channel_names), names),
        copy(trace.output_node_indices),
        vcat(copy(trace.output_pu), values),
        copy(trace.node_maximum_values),
        copy(trace.node_maximum_times_s),
        copy(trace.node_minimum_values),
        copy(trace.node_minimum_times_s),
        vcat(copy(trace.output_maximum_values), extrema.maximum_values),
        vcat(copy(trace.output_maximum_times_s), extrema.maximum_times_s),
        vcat(copy(trace.output_minimum_values), extrema.minimum_values),
        vcat(copy(trace.output_minimum_times_s), extrema.minimum_times_s),
    )
end

function _with_nonlinear_reports(run, nonlinear_config::NamedTuple)
    reports = NamedTuple()
    types = Int.(get(nonlinear_config, :nonlinear_types, Int[]))
    if any(==(SWITCHING_NONLINEAR_RESISTOR_TYPE), types)
        reports = merge(
            reports,
            (
                switching_nonlinear_resistor_report =
                    _switching_nonlinear_resistor_report(run, nonlinear_config),
            ),
        )
    end
    if any(==(TRIGGERED_TIMED_RESISTANCE_TYPE), types)
        reports = merge(
            reports,
            (
                triggered_timed_resistance_report =
                    _triggered_timed_resistance_report(run, nonlinear_config),
            ),
        )
    end
    if any(
        type_code ->
            type_code in (
                PIECEWISE_NONLINEAR_INDUCTOR_TYPE,
                PSEUDO_NONLINEAR_INDUCTOR_TYPE,
            ),
        types,
    )
        reports = merge(
            reports,
            (
                piecewise_nonlinear_inductor_report =
                    _piecewise_nonlinear_inductor_report(run, nonlinear_config),
            ),
        )
    end
    if any(type -> type in _REPORTED_NONLINEAR_BRANCH_TYPES, types)
        reports = merge(
            reports,
            (nonlinear_branch_report = _nonlinear_branch_report(run, nonlinear_config),),
        )
    end
    if !isempty(get(nonlinear_config, :safety_shunt_resistor_names, Symbol[]))
        reports = merge(
            reports,
            (
                nonlinear_isolation_shunt_report =
                    _nonlinear_isolation_shunt_report(nonlinear_config),
            ),
        )
    end
    return merge(run, reports)
end

function _run_primary_nonlinear_deck(
    parsed::DeckParser.DeckParseResult;
    dt_s::Float64,
    t_end_s::Float64,
    recorded_step_indices = nothing,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    context = initialize_step_context(
        parsed;
        dt_s = dt_s,
        t_end_s = t_end_s,
        recorded_step_indices = recorded_step_indices,
        source_signal_provider = source_signal_provider,
    )
    all(
        element -> element isa Union{
            ConductanceBranch,
            TheveninSource,
            CurrentInjection,
            TACSControlledSwitch,
            SaturatedTransformerNonlinearSlopeBranch,
        },
        context.system.elements,
    ) || throw(ArgumentError(
        "primary dense nonlinear runtime requires memoryless linear or controlled-switch owners",
    ))
    nonlinear_current_config = _deck_nonlinear_current_config(
        parsed,
        nothing;
        delta2 = dt_s / 2.0,
    )
    nonlinear_current_config === nothing &&
        throw(ArgumentError("primary nonlinear runtime requires parsed nonlinear owners"))
    if any(
        ==(PSEUDO_NONLINEAR_INDUCTOR_TYPE),
        nonlinear_current_config.nonlinear_types,
    )
        nonlinear_current_config =
            _pseudo_nonlinear_inductor_initial_state_config(
                nonlinear_current_config,
                deck_steady_state_voltage_phasors(parsed),
            )
    end
    includes_inductive_state = any(
        type ->
            type in (
                HYSTERETIC_INDUCTOR_NONLINEAR_TYPE,
                PIECEWISE_NONLINEAR_INDUCTOR_TYPE,
                PSEUDO_NONLINEAR_INDUCTOR_TYPE,
            ),
        nonlinear_current_config.nonlinear_types,
    )
    state = _deck_dynamic_timestep_state(
        parsed,
        context.system.node_count,
        dt_s,
    )
    if !includes_inductive_state
        initial_sample = _deck_runtime_initial_voltage_sample(parsed, :zero)
        _apply_steady_state_initial_sample!(context, initial_sample, nothing)
        context.step_index = 1
        context.t_s = min(context.step_index, context.step_count) * context.dt_s
    end

    updates = Any[]
    accepted_voltage_values = fill(NaN, context.system.node_count)
    while context.step_index <= context.step_count
        nonlinear_step_config = _deck_nonlinear_current_step_config(
            nonlinear_current_config,
            state,
            context,
        )
        if !includes_inductive_state
            nonlinear_step_config = merge(
                nonlinear_step_config,
                (
                    nonlinear_current_config = merge(
                        nonlinear_step_config.nonlinear_current_config,
                        (
                            accepted_voltage_values = accepted_voltage_values,
                            accepted_voltage_tolerance = 1.0e-12,
                        ),
                    ),
                ),
            )
        end
        update = step_with_over16_boundary!(
            context,
            state;
            nonlinear_step_config...,
            dense_primary_nonlinear_compensation = true,
        )
        push!(updates, update)
    end
    run = (
        context = context,
        trace = deck_trace(context),
        nonlinear_state = state.nonlinear_inverse,
        updates = updates,
        over16_updates = updates,
        nonlinear_owner_count = length(nonlinear_current_config.nonlinear_types),
        nonlinear_inverse_build_count = state.nonlinear_inverse.update_count,
        hysteretic_companion_admittance_stamp_count = sum(
            update -> update.hysteretic_companion_admittance_stamp_count,
            updates;
            init = 0,
        ),
    )
    return _with_nonlinear_reports(run, nonlinear_current_config)
end

function _ensure_float_vector_length!(values::Vector{Float64}, count::Int)
    old_count = length(values)
    resize!(values, count)
    for index in (old_count + 1):count
        values[index] = 0.0
    end
    return values
end

function _deck_branch_output_snapshot(
    context::EMTStepContext,
    branch_index::Int,
    voltage::AbstractVector{Float64},
)
    1 <= branch_index <= length(context.system.elements) ||
        throw(ArgumentError("deck branch output branch index $branch_index is outside deck elements"))
    element = context.system.elements[branch_index]
    applicable(branch_companion_snapshot, element, voltage, context.dt_s) ||
        throw(ArgumentError("deck branch output requires an accepted branch companion owner"))
    snapshot = branch_companion_snapshot(element, voltage, context.dt_s)
    snapshot !== nothing ||
        throw(ArgumentError("deck branch output requires an accepted branch companion snapshot"))
    return snapshot
end

function _deck_branch_power_snapshot(context::EMTStepContext,
                                     branch_index::Int,
                                     voltage::AbstractVector{Float64})
    1 <= branch_index <= length(context.system.elements) ||
        throw(ArgumentError("deck branch power output branch index $branch_index is outside deck elements"))
    snapshot = branch_companion_snapshot(context.system.elements[branch_index], voltage, context.dt_s)
    snapshot !== nothing ||
        throw(ArgumentError("deck branch power output requires an accepted scalar branch companion owner"))
    return snapshot
end

function _deck_branch_current_value(
    context::EMTStepContext,
    branch_index::Int,
    voltage::AbstractVector{Float64},
)
    1 <= branch_index <= length(context.system.elements) ||
        throw(ArgumentError("deck branch output branch index $branch_index is outside deck elements"))
    return branch_current_value(
        context.system.elements[branch_index],
        voltage,
        context.dt_s,
    )
end

function _prepare_output_report_branch_output_config!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    config::NamedTuple,
    voltage::AbstractVector{Float64},
    report_step_s::Float64 = context.dt_s,
)
    has_voltage_config = haskey(config, :deck_branch_voltage_output_config)
    has_current_config = haskey(config, :deck_branch_current_output_config)
    has_power_config = haskey(config, :deck_branch_power_output_config)
    (has_voltage_config || has_current_config || has_power_config) || return config

    voltage_branch_indices =
        has_voltage_config ? collect(config.deck_branch_voltage_output_config.branch_indices) : Int[]
    current_branch_indices =
        has_current_config ? collect(config.deck_branch_current_output_config.branch_indices) : Int[]
    power_branch_indices =
        has_power_config ? collect(config.deck_branch_power_output_config.branch_indices) : Int[]
    bvalue_branch_indices = Int[]
    voltage_slot_by_branch = Dict{Int,Int}()
    current_slots = Set{Int}()
    for branch_index in voltage_branch_indices
        push!(bvalue_branch_indices, branch_index)
        get!(voltage_slot_by_branch, branch_index, length(bvalue_branch_indices))
    end
    for branch_index in current_branch_indices
        push!(bvalue_branch_indices, branch_index)
        push!(current_slots, length(bvalue_branch_indices))
    end

    power_voltage_selectors = Int[]
    for branch_index in power_branch_indices
        slot = get(voltage_slot_by_branch, branch_index, 0)
        if slot == 0
            push!(bvalue_branch_indices, branch_index)
            slot = length(bvalue_branch_indices)
            voltage_slot_by_branch[branch_index] = slot
        end
        push!(power_voltage_selectors, slot)
    end

    power_current_selectors = Int[]
    for branch_index in power_branch_indices
        push!(bvalue_branch_indices, branch_index)
        slot = length(bvalue_branch_indices)
        push!(current_slots, slot)
        push!(power_current_selectors, slot)
    end

    _ensure_float_vector_length!(over16_state.output_report.bvalue_values, length(bvalue_branch_indices))
    _ensure_float_vector_length!(over16_state.output_report.bnrg_values, length(power_branch_indices))
    for (slot, branch_index) in enumerate(bvalue_branch_indices)
        snapshot = _deck_branch_power_snapshot(context, branch_index, voltage)
        over16_state.output_report.bvalue_values[slot] =
            slot in current_slots ?
            _deck_branch_current_value(context, branch_index, voltage) :
            snapshot.branch_voltage
    end
    base_kwargs = haskey(config, :kwargs) ? config.kwargs : NamedTuple()
    output_kwargs = merge(
        base_kwargs,
        (
            branch_output_count = length(bvalue_branch_indices),
            voltage_selectors = power_voltage_selectors,
            current_selectors = power_current_selectors,
            delta2 = report_step_s / 2.0,
        ),
    )
    prepared = (
        t = config.t,
        deltat = report_step_s,
        istep = config.istep,
        kwargs = output_kwargs,
    )
    return haskey(config, :voltages) ? merge(prepared, (voltages = config.voltages,)) : prepared
end

function _prepare_over16_post_extrema_config!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    config::NamedTuple,
    report_step_s::Float64 = context.dt_s,
)
    post_extrema = over16_state.post_extrema
    post_extrema.t = context.t_s
    post_extrema.deltat = report_step_s
    post_extrema.istep = context.step_index
    base_kwargs = haskey(config, :kwargs) ? config.kwargs : NamedTuple()
    base_kwargs isa NamedTuple ||
        throw(ArgumentError("post_extrema_config.kwargs must be a NamedTuple"))
    haskey(base_kwargs, :begmax_values) && return config
    return merge(
        config,
        (
            kwargs = merge(
                base_kwargs,
                (
                    begmax_values = _deck_over16_post_extrema_begmax_values(over16_state),
                ),
            ),
        ),
    )
end

function _prepare_saturated_transformer_sparse_config(
    context::EMTStepContext,
    sparse_config::NamedTuple,
)
    has_workspace =
        haskey(sparse_config, :km) &&
        haskey(sparse_config, :ykm) &&
        haskey(sparse_config, :kks)
    derive_from_step_admittance =
        get(sparse_config, :derive_from_step_admittance, !has_workspace)
    derive_from_step_admittance || return sparse_config
    admittance = _reference_augmented_step_admittance(context)
    row_data = _reference_admittance_sparse_rows!(
        context.reference_sparse_columns,
        context.reference_sparse_values,
        context.reference_sparse_row_ends,
        admittance;
        zero_tolerance = get(
            sparse_config,
            :sparse_admittance_zero_tolerance,
            get(sparse_config, :zero_tolerance, 0.0),
        ),
    )
    partition_boundary = Int(
        get(sparse_config, :partition_boundary, size(admittance, 1)),
    )
    return merge(
        sparse_config,
        (
            km = row_data.km,
            ykm = row_data.ykm,
            kks = row_data.kks,
            partition_boundary = partition_boundary,
            sparse_factor_partition_boundary = Int(
                get(
                    sparse_config,
                    :sparse_factor_partition_boundary,
                    partition_boundary,
                ),
            ),
            sparse_factor_first_factor_row = Int(
                get(sparse_config, :sparse_factor_first_factor_row, 2),
            ),
            step_admittance_sparse_row_count = row_data.row_count,
            step_admittance_sparse_entry_count = row_data.entry_count,
        ),
    )
end

function _prepare_saturated_transformer_nonlinear_sparse_config(
    context::EMTStepContext,
    nonlinear_current_config::NamedTuple,
)
    haskey(nonlinear_current_config, :saturated_transformer_sparse_config) ||
        return nonlinear_current_config
    sparse_config = nonlinear_current_config.saturated_transformer_sparse_config
    sparse_config isa NamedTuple ||
        throw(ArgumentError("saturated transformer sparse config must be a NamedTuple"))
    return merge(
        nonlinear_current_config,
        (
            saturated_transformer_sparse_config =
                _prepare_saturated_transformer_sparse_config(context, sparse_config),
        ),
    )
end

function _prepare_over16_step_kwargs!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    voltage::AbstractVector{Float64},
    kwargs::NamedTuple,
    ;
    report_step_s::Float64 = context.dt_s,
)
    prepared = kwargs
    if haskey(prepared, :nonlinear_current_config)
        nonlinear_current_config = prepared.nonlinear_current_config
        if get(nonlinear_current_config, :deck_owned_voltage_context, false) &&
           !haskey(nonlinear_current_config, :voltages)
            voltage_context = _deck_nonlinear_voltage_context(
                nonlinear_current_config,
                voltage,
            )
            _deck_saturated_transformer_source_context!(
                over16_state,
                length(voltage_context),
            )
            prepared = merge(
                prepared,
                (
                    nonlinear_current_config = merge(
                        nonlinear_current_config,
                        (voltages = voltage_context,),
                    ),
                ),
            )
        end
    end
    if haskey(prepared, :nonlinear_current_config)
        prepared = merge(
            prepared,
            (
                nonlinear_current_config =
                    _prepare_saturated_transformer_nonlinear_sparse_config(
                        context,
                        prepared.nonlinear_current_config,
                    ),
            ),
        )
    end
    if haskey(prepared, :output_report_config)
        output_report_config = _prepare_output_report_branch_output_config!(
            context,
            over16_state,
            prepared.output_report_config,
            voltage,
            report_step_s,
        )
        prepared = merge(prepared, (output_report_config = output_report_config,))
    end
    if haskey(prepared, :post_extrema_config)
        post_extrema_config = _prepare_over16_post_extrema_config!(
            context,
            over16_state,
            prepared.post_extrema_config,
            report_step_s,
        )
        prepared = merge(prepared, (post_extrema_config = post_extrema_config,))
    end
    return prepared
end

function run_deck_emt_with_over16_boundary(
    lines,
    over16_state::OVER16AcceptedTimestepState;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    source::AbstractString = "deck",
    over16_step_configs = nothing,
    merge_deck_owned_step_configs::Bool = true,
    saturated_transformer_intake = nothing,
    saturated_transformer_deck_lines = nothing,
    saturated_transformer_winding_number::Int = 1,
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple} = nothing,
    saturated_transformer_nonlinear_current_enabled::Bool = true,
    coupled_lumped_sequence_history_enabled::Bool = false,
    steady_state_initial_sample_enabled::Bool = false,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    parsed = DeckParser.parse_deck_lines(lines; source = source)
    return run_deck_emt_with_over16_boundary(
        parsed,
        over16_state;
        dt_s = dt_s,
        t_end_s = t_end_s,
        over16_step_configs = over16_step_configs,
        merge_deck_owned_step_configs = merge_deck_owned_step_configs,
        saturated_transformer_intake = saturated_transformer_intake,
        saturated_transformer_deck_lines = saturated_transformer_deck_lines,
        saturated_transformer_winding_number = saturated_transformer_winding_number,
        saturated_transformer_sparse_config = saturated_transformer_sparse_config,
        saturated_transformer_nonlinear_current_enabled =
            saturated_transformer_nonlinear_current_enabled,
        coupled_lumped_sequence_history_enabled =
            coupled_lumped_sequence_history_enabled,
        steady_state_initial_sample_enabled =
            steady_state_initial_sample_enabled,
        time_switch_event_delay_s = time_switch_event_delay_s,
        current_zero_switching = current_zero_switching,
        recorded_step_indices = recorded_step_indices,
        series_rlc_alterations = series_rlc_alterations,
        source_signal_provider = source_signal_provider,
    )
end

mutable struct DynamicDeckStepConfigProvider{N,D,I,E,F,P,T,O,L,R} <: Function
    plan::DeckOVER16BoundaryPlan
    timestep_state::OVER16AcceptedTimestepState
    nonlinear_current_config::N
    distributed_line_config::D
    steady_state_initial_current_injections::I
    user_step_configs::E
    features::F
    time_switch_workspace::DeckTimeSwitchStepWorkspace
    current_source_seed_clear_nodes::Vector{Int}
    first_dynamic_step_index::Int
    current_source_seed_element_count::Int
    include_post_extrema::Val{P}
    include_time_switch_sparse::Val{T}
    include_output::Val{O}
    include_distributed_line::Val{L}
    reset_current_source_values::Val{R}
end

_steady_state_step_current_config!(::Nothing, ::EMTStepContext, ::Int) =
    NamedTuple()

function _steady_state_step_current_config!(
    values::Vector{Float64},
    context::EMTStepContext,
    first_dynamic_step_index::Int,
)
    context.step_index == first_dynamic_step_index || fill!(values, 0.0)
    return (current_injection_values = values,)
end

function (provider::DynamicDeckStepConfigProvider{N,D,I,E,F,P,T,O,L,R})(
    context::EMTStepContext,
) where {N,D,I,E,F,P,T,O,L,R}
    return merge(
        _deck_over16_step_config(
            provider.plan,
            context,
            provider.timestep_state,
            provider.nonlinear_current_config,
            provider.features,
            provider.time_switch_workspace,
            provider.include_post_extrema,
            provider.include_time_switch_sparse,
            context.dt_s,
        ),
        provider.nonlinear_current_config !== nothing &&
        get(
            provider.nonlinear_current_config,
            :dense_primary_nonlinear_compensation,
            false,
        ) ?
            (dense_primary_nonlinear_compensation = true,) :
            NamedTuple(),
        O ?
            (record_presolve_voltage_state = true,) :
            NamedTuple(),
        L ?
            (distributed_transposed_line_config = provider.distributed_line_config,) :
            NamedTuple(),
        R ?
            (
                reset_current_source_values = true,
                seed_current_source_values = true,
                current_source_seed_element_count =
                    provider.current_source_seed_element_count,
                current_source_seed_clear_nodes =
                    provider.current_source_seed_clear_nodes,
            ) :
            NamedTuple(),
        _steady_state_step_current_config!(
            provider.steady_state_initial_current_injections,
            context,
            provider.first_dynamic_step_index,
        ),
        _over16_step_kwargs(provider.user_step_configs, context),
    )
end

struct PreparedDynamicDeckRuntime{C,F,SS,SO}
    context::C
    timestep_state::OVER16AcceptedTimestepState
    step_configs::F
    plan::DeckOVER16BoundaryPlan
    steady_state_initial_sample::SS
    steady_state_initial_output_values::SO
    deck_owned_output::Bool
    deck_owned_source::Bool
    deck_owned_source_card::Bool
    deck_owned_post_extrema::Bool
    deck_owned_time_switch_sparse::Bool
    store_step_updates::Bool
end

function _step_with_over16_config!(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState,
    config::C,
    collect_step_diagnostics::Bool,
) where {C<:NamedTuple}
    return step_with_over16_boundary!(
        context,
        over16_state;
        collect_step_diagnostics = collect_step_diagnostics,
        config...,
    )
end

function _advance_prepared_emt_step!(
    context::EMTStepContext,
    timestep_state::OVER16AcceptedTimestepState,
    step_configs,
    collect_step_diagnostics::Bool,
)
    context.step_index <= context.step_count || throw(ArgumentError(
        "prepared EMT runtime has no remaining timestep to advance",
    ))
    _apply_due_series_rlc_alterations!(context)
    over16_kwargs = _over16_step_kwargs(step_configs, context)
    return _step_with_over16_config!(
        context,
        timestep_state,
        over16_kwargs,
        collect_step_diagnostics,
    )
end

function _advance_prepared_emt_step!(
    runtime::PreparedDynamicDeckRuntime;
    collect_step_diagnostics::Bool=false,
)
    return _advance_prepared_emt_step!(
        runtime.context,
        runtime.timestep_state,
        runtime.step_configs,
        collect_step_diagnostics,
    )
end

struct EMTStepTransaction{R,T}
    runtime::R
    transaction::T
end

function _restore_prepared_dynamic_deck_runtime!(
    destination::PreparedDynamicDeckRuntime,
    source::PreparedDynamicDeckRuntime,
    restorer::TimestepStateRestorer,
)
    restore_timestep_state!(destination.context, source.context, restorer)
    restore_timestep_state!(
        destination.timestep_state,
        source.timestep_state,
        restorer,
    )
    restore_timestep_state!(
        destination.step_configs,
        source.step_configs,
        restorer,
    )
    return destination
end

function EMTStepTransaction(
    workspace::EMTStudyWorkspace;
    stable_structure::Bool=false,
)
    workspace.ready || throw(ArgumentError(
        "EMT study workspace must be ready before creating a step transaction",
    ))
    runtime = _check_prepared_runtime_aliases(workspace.runtime)
    transaction = stable_structure ?
        StableStructureTimestepTransaction(
            runtime,
            _restore_prepared_dynamic_deck_runtime!,
        ) :
        TimestepTransaction(runtime)
    return EMTStepTransaction(runtime, transaction)
end

function initialize_partitioned_emt_workspace!(
    workspace::EMTStudyWorkspace;
    execution_mode::Symbol=:partitioned_waveform,
)
    workspace.ready || throw(ArgumentError(
        "EMT partition region must begin from an unclaimed prepared workspace",
    ))
    workspace.execution_mode === :unselected || throw(ArgumentError(
        "EMT partition region workspace is already owned by $(workspace.execution_mode) execution",
    ))
    runtime = workspace.runtime
    context = runtime.context
    context.step_index == 0 && context.t_s == 0.0 || throw(ArgumentError(
        "EMT partition region must begin at its prepared time-zero boundary",
    ))
    context.step_count > 0 || throw(ArgumentError(
        "EMT partition region requires at least one physical timestep",
    ))
    _apply_due_series_rlc_alterations!(context)
    if runtime.steady_state_initial_sample === nothing
        record_step!(context, copy(context.system.v))
    else
        _apply_steady_state_initial_sample!(
            context,
            runtime.steady_state_initial_sample,
            runtime.steady_state_initial_output_values,
        )
        context.step_index = 1
        context.t_s = context.dt_s
    end
    execution_mode in (
        :monolithic,
        :local_subcycling,
        :partitioned_lagged,
        :partitioned_waveform,
        :combined_local_partitioned,
    ) || throw(ArgumentError(
        "EMT partition region execution mode is unsupported",
    ))
    workspace.execution_mode = execution_mode
    workspace.ready = false
    return workspace
end

function begin_emt_step_transaction!(transaction::EMTStepTransaction)
    begin_timestep_transaction!(transaction.transaction)
    return transaction
end

function provisional_emt_step!(
    transaction::EMTStepTransaction;
    collect_step_diagnostics::Bool=false,
)
    timestep_transaction_active(transaction.transaction) ||
        throw(ArgumentError("EMT step transaction must be active before advance"))
    return _advance_prepared_emt_step!(
        transaction.runtime;
        collect_step_diagnostics = collect_step_diagnostics,
    )
end

function restore_emt_step_transaction!(transaction::EMTStepTransaction)
    restore_timestep_transaction!(transaction.transaction)
    _check_prepared_runtime_aliases(transaction.runtime)
    return transaction
end

function commit_emt_step_transaction!(transaction::EMTStepTransaction)
    commit_timestep_transaction!(transaction.transaction)
    _check_prepared_runtime_aliases(transaction.runtime)
    return transaction
end

emt_step_transaction_status(transaction::EMTStepTransaction) =
    timestep_transaction_status(transaction.transaction)

function _run_prepared_dynamic_deck!(
    runtime::PreparedDynamicDeckRuntime;
    collect_run_diagnostics::Bool=true,
)
    run = run_deck_emt_with_over16_boundary(
        runtime.context,
        runtime.timestep_state;
        over16_step_configs = runtime.step_configs,
        deck_over16_boundary_plan = runtime.plan,
        deck_owned_over16_output_config = runtime.deck_owned_output,
        deck_owned_over5a_source_config = runtime.deck_owned_source,
        deck_owned_over16_source_card_config = runtime.deck_owned_source_card,
        deck_owned_over16_post_extrema_config = runtime.deck_owned_post_extrema,
        deck_owned_time_switch_sparse_config =
            runtime.deck_owned_time_switch_sparse,
        steady_state_initial_sample = runtime.steady_state_initial_sample,
        steady_state_initial_output_values =
            runtime.steady_state_initial_output_values,
        store_step_updates = runtime.store_step_updates,
        collect_run_diagnostics = collect_run_diagnostics,
    )
    collect_run_diagnostics || return run
    nonlinear_config =
        runtime.step_configs isa DynamicDeckStepConfigProvider ?
        runtime.step_configs.nonlinear_current_config : nothing
    nonlinear_config === nothing && return run
    return _with_nonlinear_reports(run, nonlinear_config)
end

function _prepare_dynamic_deck_runtime(
    parsed::DeckParser.DeckParseResult,
    over16_state::OVER16AcceptedTimestepState;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    over16_step_configs = nothing,
    merge_deck_owned_step_configs::Bool = true,
    saturated_transformer_intake = nothing,
    saturated_transformer_deck_lines = nothing,
    saturated_transformer_winding_number::Int = 1,
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple} = nothing,
    saturated_transformer_nonlinear_current_enabled::Bool = true,
    coupled_lumped_sequence_history_enabled::Bool = false,
    steady_state_initial_sample_enabled::Bool = false,
    supplied_initial_sample = nothing,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    store_step_updates::Bool = true,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    plan = deck_over16_boundary_plan(parsed)
    resolved_saturated_transformer_intake = _deck_saturated_transformer_intake(
        saturated_transformer_intake,
        saturated_transformer_deck_lines,
        parsed.source,
    )
    saturated_transformer_current_config =
        resolved_saturated_transformer_intake === nothing ||
        !saturated_transformer_nonlinear_current_enabled ? nothing :
        deck_saturated_transformer_nonlinear_current_config(
            parsed,
            resolved_saturated_transformer_intake;
            delta2 = dt_s / 2.0,
            winding_number = saturated_transformer_winding_number,
            saturated_transformer_sparse_config = saturated_transformer_sparse_config,
        )
    steady_state_initial_sample =
        supplied_initial_sample !== nothing ?
        deepcopy(supplied_initial_sample) :
        steady_state_initial_sample_enabled ?
        deck_steady_state_voltage_phasors(
            parsed;
            saturated_transformer_intake = resolved_saturated_transformer_intake,
            winding_number = saturated_transformer_winding_number,
        ) :
        nothing
    consistent_initial_state_enabled = supplied_initial_sample !== nothing
    if consistent_initial_state_enabled
        saturated_transformer_current_config =
            _saturated_transformer_initial_state_config(
                saturated_transformer_current_config,
                steady_state_initial_sample,
            )
    end
    saturated_transformer_current_config =
        _pseudo_nonlinear_inductor_initial_state_config(
            saturated_transformer_current_config,
            steady_state_initial_sample,
        )
    deck_nonlinear_current_config =
        _deck_nonlinear_current_config(
            parsed,
            saturated_transformer_current_config;
            delta2 = dt_s / 2.0,
        )
    if consistent_initial_state_enabled
        deck_nonlinear_current_config =
            _saturated_transformer_initial_state_config(
                deck_nonlinear_current_config,
                steady_state_initial_sample,
            )
        deck_nonlinear_current_config =
            _pseudo_nonlinear_inductor_initial_state_config(
                deck_nonlinear_current_config,
                steady_state_initial_sample,
            )
        deck_nonlinear_current_config =
            _piecewise_nonlinear_inductor_initial_state_config(
                deck_nonlinear_current_config,
                steady_state_initial_sample,
            )
        deck_nonlinear_current_config =
            _hysteretic_inductor_initial_state_config(
                deck_nonlinear_current_config,
                steady_state_initial_sample,
            )
    end
    if deck_nonlinear_current_config !== nothing &&
       saturated_transformer_current_config === nothing &&
       _deck_uses_dynamic_nonlinear_runtime(parsed)
        deck_nonlinear_current_config = merge(
            deck_nonlinear_current_config,
            (dense_primary_nonlinear_compensation = true,),
        )
    end
    distributed_transposed_line_config =
        merge_deck_owned_step_configs ?
        _deck_distributed_transposed_line_config(
            parsed;
            steady_state_initial_sample = steady_state_initial_sample,
        ) :
        nothing
    deck_owned_over16_output_config =
        merge_deck_owned_step_configs &&
        (
            !isempty(plan.output_node_indices) ||
            !isempty(plan.branch_voltage_branch_indices) ||
            !isempty(plan.branch_current_branch_indices) ||
            !isempty(plan.branch_power_branch_indices)
        )
    deck_owned_over5a_source_config =
        merge_deck_owned_step_configs && plan.source_row_count > 0
    deck_owned_over16_source_card_config =
        deck_owned_over5a_source_config && plan.source_card_row_count > 0
    deck_owned_time_switch_sparse_config =
        merge_deck_owned_step_configs && plan.switch_count > 0 &&
        !current_zero_switching
    deck_owned_saturated_transformer_current_config =
        merge_deck_owned_step_configs && saturated_transformer_current_config !== nothing
    deck_owned_nonlinear_current_config =
        merge_deck_owned_step_configs && deck_nonlinear_current_config !== nothing
    deck_owned_saturated_transformer_branch_context =
        merge_deck_owned_step_configs &&
        resolved_saturated_transformer_intake !== nothing
    deck_owned_coupled_lumped_sequence_history =
        merge_deck_owned_step_configs &&
        coupled_lumped_sequence_history_enabled &&
        !isempty(DeckParser.deck_coupled_lumped_sequence_impedances(parsed))
    deck_owned_distributed_transposed_line_config =
        distributed_transposed_line_config !== nothing
    deck_owned_over16_post_extrema_config =
        deck_owned_over16_output_config ||
        deck_owned_over5a_source_config ||
        deck_owned_over16_source_card_config
    context =
        deck_owned_saturated_transformer_current_config ?
        saturated_transformer_augmented_step_context(
            parsed,
            saturated_transformer_current_config;
            dt_s = dt_s,
            t_end_s = t_end_s,
            include_coupled_lumped_sequence_history =
                deck_owned_coupled_lumped_sequence_history,
            time_switch_event_delay_s = time_switch_event_delay_s,
            current_zero_switching = current_zero_switching,
            recorded_step_indices = recorded_step_indices,
            source_signal_provider = source_signal_provider,
        ) :
        deck_owned_saturated_transformer_branch_context ?
        saturated_transformer_branch_augmented_step_context(
            parsed,
            resolved_saturated_transformer_intake;
            dt_s = dt_s,
            t_end_s = t_end_s,
            winding_number = saturated_transformer_winding_number,
            transformer_branch_shunt_capacitance_rows =
                _deck_transformer_branch_shunt_capacitance_rows(
                    parsed,
                    saturated_transformer_deck_lines,
                ),
            include_coupled_lumped_sequence_history =
                deck_owned_coupled_lumped_sequence_history,
            time_switch_event_delay_s = time_switch_event_delay_s,
            current_zero_switching = current_zero_switching,
            recorded_step_indices = recorded_step_indices,
            source_signal_provider = source_signal_provider,
        ) :
        deck_owned_coupled_lumped_sequence_history ?
        coupled_lumped_sequence_augmented_step_context(
            parsed;
            dt_s = dt_s,
            t_end_s = t_end_s,
            time_switch_event_delay_s = time_switch_event_delay_s,
            current_zero_switching = current_zero_switching,
            recorded_step_indices = recorded_step_indices,
            source_signal_provider = source_signal_provider,
        ) :
        initialize_step_context(
            parsed;
            dt_s = dt_s,
            t_end_s = t_end_s,
            recorded_step_indices = recorded_step_indices,
            time_switch_event_delay_s = time_switch_event_delay_s,
            current_zero_switching = current_zero_switching,
            source_signal_provider = source_signal_provider,
        )
    if context.source_function_runtime !== nothing
        context.source_function_runtime.plan = plan
    end
    configure_series_rlc_alterations!(context, series_rlc_alterations)
    if distributed_transposed_line_config !== nothing &&
       haskey(distributed_transposed_line_config, :current_injection_values)
        resize!(
            distributed_transposed_line_config.current_injection_values,
            context.system.node_count,
        )
        fill!(distributed_transposed_line_config.current_injection_values, 0.0)
    end
    if steady_state_initial_sample !== nothing
        _seed_steady_state_network_state!(context, steady_state_initial_sample)
        length(over16_state.source.e_values) >= context.system.node_count ||
            throw(ArgumentError("source voltage state is smaller than the initialized network"))
        over16_state.source.e_values[1:context.system.node_count] .= context.system.v
    end
    steady_state_initial_sample === nothing ||
        _initialize_control_system_network_steady_state!(
            context,
            context.system.v,
        )
    initialize_dc_simulator_sources!(plan, over16_state.source.e_values)
    _initialize_source_function_state!(context, over16_state, plan)
    steady_state_initial_output_values =
        steady_state_initial_sample === nothing ?
        nothing :
        _steady_state_initial_output_values(context)
    steady_state_initial_current_injections =
        steady_state_initial_sample === nothing ||
        (
            hasproperty(steady_state_initial_sample, :exact_discrete_histories) &&
            steady_state_initial_sample.exact_discrete_histories
        ) ? nothing :
        _steady_state_current_injections(context, steady_state_initial_sample)
    first_dynamic_step_index = steady_state_initial_sample === nothing ? 0 : 1
    reset_current_source_values =
        deck_owned_distributed_transposed_line_config ||
        deck_owned_saturated_transformer_current_config
    current_source_seed_clear_nodes =
        reset_current_source_values ?
        _deck_current_source_seed_clear_nodes(parsed) :
        Int[]
    effective_step_configs =
        (deck_owned_over16_post_extrema_config ||
         deck_owned_time_switch_sparse_config ||
         deck_owned_nonlinear_current_config ||
         deck_owned_distributed_transposed_line_config) ?
        DynamicDeckStepConfigProvider(
            plan,
            over16_state,
            deck_owned_nonlinear_current_config ?
                deck_nonlinear_current_config : nothing,
            distributed_transposed_line_config,
            steady_state_initial_current_injections,
            over16_step_configs,
            DeckStepConfigFeatures{
                deck_owned_over16_output_config,
                !isempty(plan.branch_voltage_branch_indices),
                !isempty(plan.branch_current_branch_indices),
                !isempty(plan.branch_power_branch_indices),
                plan.source_row_count > 0,
                plan.source_card_row_count > 0 ||
                    plan.source_interpolation_row_count > 0,
                plan.source_tacs_override_count > 0 ||
                    any(
                        source_type ->
                            abs(source_type) == 17 || abs(source_type) >= 60,
                        plan.source_iform_values,
                    ),
                plan.source_analytic_row_count > 0 ||
                    (
                        context.source_function_runtime !== nothing &&
                        (
                            source_signal_analytic_active(
                                context.source_function_runtime.signal_provider,
                            ) ||
                            context.source_function_runtime.internal_analytic_requested
                        )
                    ),
                context.source_function_runtime !== nothing &&
                    _source_signal_stage_recording_active(
                        context.source_function_runtime,
                    ),
            }(),
            DeckTimeSwitchStepWorkspace(
                plan,
                context.system.node_count + 1,
            ),
            current_source_seed_clear_nodes,
            first_dynamic_step_index,
            length(context.system.elements),
            Val(deck_owned_over16_post_extrema_config),
            Val(deck_owned_time_switch_sparse_config),
            Val(deck_owned_over16_output_config),
            Val(deck_owned_distributed_transposed_line_config),
            Val(reset_current_source_values),
        ) :
        over16_step_configs
    return PreparedDynamicDeckRuntime(
        context,
        over16_state,
        effective_step_configs,
        plan,
        steady_state_initial_sample,
        steady_state_initial_output_values,
        deck_owned_over16_output_config,
        deck_owned_over5a_source_config,
        deck_owned_over16_source_card_config,
        deck_owned_over16_post_extrema_config,
        deck_owned_time_switch_sparse_config,
        store_step_updates,
    )
end

function run_deck_emt_with_over16_boundary(
    parsed::DeckParser.DeckParseResult,
    over16_state::OVER16AcceptedTimestepState;
    dt_s::Float64 = 20e-6,
    t_end_s::Float64 = 0.0,
    over16_step_configs = nothing,
    merge_deck_owned_step_configs::Bool = true,
    saturated_transformer_intake = nothing,
    saturated_transformer_deck_lines = nothing,
    saturated_transformer_winding_number::Int = 1,
    saturated_transformer_sparse_config::Union{Nothing,NamedTuple} = nothing,
    saturated_transformer_nonlinear_current_enabled::Bool = true,
    coupled_lumped_sequence_history_enabled::Bool = false,
    steady_state_initial_sample_enabled::Bool = false,
    time_switch_event_delay_s::Float64 = 0.0,
    current_zero_switching::Bool = false,
    recorded_step_indices = nothing,
    series_rlc_alterations::AbstractVector{<:SeriesRLCAlteration} =
        SeriesRLCAlteration[],
    store_step_updates::Bool = true,
    source_signal_provider::AbstractSourceSignalProvider = IdentitySourceSignalProvider(),
)
    runtime = _prepare_dynamic_deck_runtime(
        parsed,
        over16_state;
        dt_s = dt_s,
        t_end_s = t_end_s,
        over16_step_configs = over16_step_configs,
        merge_deck_owned_step_configs = merge_deck_owned_step_configs,
        saturated_transformer_intake = saturated_transformer_intake,
        saturated_transformer_deck_lines = saturated_transformer_deck_lines,
        saturated_transformer_winding_number =
            saturated_transformer_winding_number,
        saturated_transformer_sparse_config =
            saturated_transformer_sparse_config,
        saturated_transformer_nonlinear_current_enabled =
            saturated_transformer_nonlinear_current_enabled,
        coupled_lumped_sequence_history_enabled =
            coupled_lumped_sequence_history_enabled,
        steady_state_initial_sample_enabled = steady_state_initial_sample_enabled,
        time_switch_event_delay_s = time_switch_event_delay_s,
        current_zero_switching = current_zero_switching,
        recorded_step_indices = recorded_step_indices,
        series_rlc_alterations = series_rlc_alterations,
        store_step_updates = store_step_updates,
        source_signal_provider = source_signal_provider,
    )
    return _run_prepared_dynamic_deck!(runtime)
end

function run_deck_emt_with_over16_boundary(
    context::EMTStepContext,
    over16_state::OVER16AcceptedTimestepState;
    over16_step_configs = nothing,
    deck_over16_boundary_plan::Union{Nothing,DeckOVER16BoundaryPlan} = nothing,
    deck_owned_over16_output_config::Bool = false,
    deck_owned_over5a_source_config::Bool = false,
    deck_owned_over16_source_card_config::Bool = false,
    deck_owned_over16_post_extrema_config::Bool = false,
    deck_owned_time_switch_sparse_config::Bool = false,
    steady_state_initial_sample = nothing,
    steady_state_initial_output_values = nothing,
    store_step_updates::Bool = true,
    collect_run_diagnostics::Bool = true,
)
    updates = collect_run_diagnostics ? Any[] : nothing
    _apply_due_series_rlc_alterations!(context)
    if steady_state_initial_sample !== nothing && context.step_index == 0
        _apply_steady_state_initial_sample!(
            context,
            steady_state_initial_sample,
            steady_state_initial_output_values,
        )
        context.step_index = 1
        context.t_s = min(context.step_index, context.step_count) * context.dt_s
    end
    while context.step_index <= context.step_count
        update = _advance_prepared_emt_step!(
            context,
            over16_state,
            over16_step_configs,
            store_step_updates,
        )
        store_step_updates && collect_run_diagnostics && push!(updates, update)
    end
    collect_run_diagnostics || return context
    trace = deck_trace(context)
    over16_update_count = sum(update -> _over16_step_pass_count(update), updates; init = 0)
    mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :accepted_timestep_state_mutated,
            :accepted_timestep_state_mutation_count,
        ),
        updates;
        init = 0,
    )
    source_update_mutation_count =
        _over16_boundary_pass_mutation_count(updates, :source_update_mutated)
    post_extrema_mutation_count =
        _over16_boundary_pass_mutation_count(updates, :post_extrema_state_mutated)
    tacs_linear_mutation_count =
        _over16_boundary_pass_mutation_count(updates, :tacs_linear_solve_state_mutated)
    tacs_post_solve_mutation_count =
        _over16_boundary_pass_mutation_count(updates, :tacs_post_solve_state_mutated)
    csup_mutation_count =
        _over16_boundary_pass_mutation_count(updates, :csup_state_mutated)
    switch_state_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_state_mutated,
            :switch_state_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_topology_admittance_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_topology_admittance_state_mutated,
            :switch_topology_admittance_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_admittance_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_admittance_state_mutated,
            :switch_admittance_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_retriangularization_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_retriangularization_state_mutated,
            :switch_retriangularization_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_sparse_factor_workspace_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_sparse_factor_workspace_mutated,
            :switch_sparse_factor_workspace_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_fortran_sparse_factor_workspace_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_fortran_sparse_factor_workspace_mutated,
            :switch_fortran_sparse_factor_workspace_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_network_solution_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_network_solution_state_mutated,
            :switch_network_solution_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_current_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_current_state_mutated,
            :switch_current_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_post_current_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_post_current_state_mutated,
            :switch_post_current_mutation_count,
        ),
        updates;
        init = 0,
    )
    switch_bvalue_mutation_count = sum(
        update -> _over16_step_mutation_count(
            update,
            :switch_bvalue_state_mutated,
            :switch_bvalue_mutation_count,
        ),
        updates;
        init = 0,
    )
    pass_updates = _over16_boundary_pass_updates(updates)
    tacs_execution_count = count(update -> update.tacs_executed, pass_updates)
    ntacs_execution_count = count(update -> update.ntacs_executed, pass_updates)
    elec_call_intent_count = count(update -> update.elec_call_intent, pass_updates)
    elec_output_copy_count =
        count(update -> update.elec_output_copy_executed, pass_updates)
    elec_tacsto_interpreter_count = count(pass_updates) do update
        result = update.ntacs3_result
        result !== nothing && result.elec_tacsto_interpreter_executed
    end
    elec_report_count = count(pass_updates) do update
        result = update.ntacs3_result
        result !== nothing && result.elec_report_executed
    end
    elec_error_count = count(pass_updates) do update
        result = update.ntacs3_result
        result !== nothing && result.elec_error_executed
    end
    csup_execution_count = count(update -> update.csup_executed, pass_updates)
    machine_solution_pass_count =
        _over16_boundary_pass_bool_count(pass_updates, :machine_solution_pass_executed)
    machine_equation_pass_count =
        _over16_boundary_pass_bool_count(pass_updates, :machine_equation_pass_executed)
    solvum_execution_count =
        _over16_boundary_pass_bool_count(pass_updates, :solvum_executed)
    full_solvum_execution_count =
        _over16_boundary_pass_bool_count(pass_updates, :full_solvum_execution)
    complete_machine_network_solution_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :complete_machine_network_solution,
        )
    machine_terminal_network_coupling_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_terminal_network_coupling,
        )
    machine_solution_pass_network_context_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_solution_pass_network_context_passed,
        )
    machine_solution_pass_output_network_coupling_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_solution_pass_output_network_coupled,
        )
    machine_solution_pass_network_coupling_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_solution_pass_network_coupled,
        )
    machine_runtime_projection_mutation_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_runtime_projection_state_mutated,
        )
    machine_runtime_state_mutation_count =
        _over16_boundary_pass_bool_count(pass_updates, :machine_runtime_state_mutated)
    machine_runtime_storage_count =
        _over16_boundary_pass_bool_count(pass_updates, :machine_runtime_storage_populated)
    machine_terminal_output_mutation_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_terminal_output_state_mutated,
        )
    machine_terminal_network_current_mutation_count =
        _over16_boundary_pass_bool_count(
            pass_updates,
            :machine_terminal_network_current_mutated,
        )
    machine_terminal_output_values =
        _over16_boundary_pass_float_values(pass_updates, :machine_terminal_output_values)
    machine_terminal_current_substitution_values =
        _over16_boundary_pass_float_values(
            pass_updates,
            :machine_terminal_current_substitution_values,
        )
    machine_terminal_network_current_values =
        _over16_boundary_pass_float_values(
            pass_updates,
            :machine_terminal_network_current_values,
        )
    machine_runtime_dynamic_values =
        _over16_boundary_pass_float_values(pass_updates, :machine_runtime_dynamic_values)
    output_channel_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.output_channel_names)
    output_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.output_node_names)
    output_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.output_node_indices)
    output_channel_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.output_channel_line_numbers)
    branch_power_output_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_power_output_names)
    branch_voltage_output_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_voltage_output_names)
    branch_voltage_branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_voltage_branch_names)
    branch_voltage_branch_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_voltage_branch_indices)
    branch_voltage_output_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_voltage_output_line_numbers)
    branch_current_output_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_current_output_names)
    branch_current_branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_current_branch_names)
    branch_current_branch_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_current_branch_indices)
    branch_current_output_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_current_output_line_numbers)
    branch_power_branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_power_branch_names)
    branch_power_branch_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_power_branch_indices)
    branch_power_output_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_power_output_line_numbers)
    over15_output_request_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_names)
    over15_output_request_output_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_output_kinds)
    over15_output_request_request_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_request_kinds)
    over15_output_request_layout_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_layout_kinds)
    over15_output_request_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over15_output_request_line_numbers)
    over15_output_request_output_codes =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over15_output_request_output_codes)
    over15_output_request_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_node_names)
    over15_output_request_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over15_output_request_node_indices)
    over15_output_request_branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over15_output_request_branch_names)
    over15_output_request_branch_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over15_output_request_branch_indices)
    branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_names)
    branch_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_kinds)
    branch_from_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_from_node_names)
    branch_to_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.branch_to_node_names)
    branch_from_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_from_node_indices)
    branch_to_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_to_node_indices)
    branch_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_conductance_values)
    branch_resistance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_resistance_values)
    branch_inductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_inductance_values)
    branch_capacitance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_capacitance_values)
    branch_previous_current_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_previous_current_values)
    branch_previous_voltage_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.branch_previous_voltage_values)
    branch_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.branch_line_numbers)
    bergeron_line_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.bergeron_line_names)
    bergeron_line_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.bergeron_line_line_numbers)
    bergeron_line_from_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.bergeron_line_from_node_names)
    bergeron_line_to_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.bergeron_line_to_node_names)
    bergeron_line_from_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.bergeron_line_from_node_indices)
    bergeron_line_to_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.bergeron_line_to_node_indices)
    bergeron_line_surge_impedance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_surge_impedance_values)
    bergeron_line_surge_admittance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_surge_admittance_values)
    bergeron_line_travel_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_travel_time_s_values)
    bergeron_line_dt_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_dt_s_values)
    bergeron_line_attenuation_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_attenuation_values)
    bergeron_line_delay_step_counts =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.bergeron_line_delay_step_counts)
    bergeron_line_write_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.bergeron_line_write_indices)
    bergeron_line_history_current_from_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_history_current_from_values)
    bergeron_line_history_current_to_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_history_current_to_values)
    bergeron_line_terminal_voltage_from_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_terminal_voltage_from_values)
    bergeron_line_terminal_voltage_to_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_terminal_voltage_to_values)
    bergeron_line_terminal_current_from_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_terminal_current_from_values)
    bergeron_line_terminal_current_to_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_terminal_current_to_values)
    bergeron_line_traveling_wave_from_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_traveling_wave_from_values)
    bergeron_line_traveling_wave_to_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.bergeron_line_traveling_wave_to_values)
    over2_branch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_names)
    over2_branch_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over2_branch_line_numbers)
    over2_branch_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_kinds)
    over2_branch_layout_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_layout_kinds)
    over2_branch_source_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_source_kinds)
    over2_branch_reference_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_reference_kinds)
    over2_branch_reference_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_reference_names)
    over2_branch_reference_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over2_branch_reference_line_numbers)
    over2_branch_from_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_from_node_names)
    over2_branch_to_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over2_branch_to_node_names)
    over2_branch_from_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over2_branch_from_node_indices)
    over2_branch_to_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over2_branch_to_node_indices)
    over2_branch_raw_resistance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_raw_resistance_values)
    over2_branch_raw_inductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_raw_inductance_values)
    over2_branch_raw_capacitance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_raw_capacitance_values)
    over2_branch_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_conductance_values)
    over2_branch_resistance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_resistance_values)
    over2_branch_inductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_inductance_values)
    over2_branch_capacitance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over2_branch_capacitance_values)
    over2_branch_output_codes =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over2_branch_output_codes)
    switch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.switch_names)
    switch_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.switch_line_numbers)
    switch_from_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.switch_from_node_names)
    switch_to_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.switch_to_node_names)
    switch_from_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.switch_from_node_indices)
    switch_to_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.switch_to_node_indices)
    switch_close_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.switch_close_time_s_values)
    switch_open_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.switch_open_time_s_values)
    switch_initially_closed_flags =
        deck_over16_boundary_plan === nothing ? Bool[] :
        copy(deck_over16_boundary_plan.switch_initially_closed_flags)
    switch_on_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.switch_on_conductance_values)
    switch_off_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.switch_off_conductance_values)
    over5_switch_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over5_switch_names)
    over5_switch_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over5_switch_line_numbers)
    over5_switch_from_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over5_switch_from_node_names)
    over5_switch_to_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over5_switch_to_node_names)
    over5_switch_from_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over5_switch_from_node_indices)
    over5_switch_to_node_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over5_switch_to_node_indices)
    over5_switch_layout_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.over5_switch_layout_kinds)
    over5_switch_raw_close_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_raw_close_time_s_values)
    over5_switch_raw_open_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_raw_open_time_s_values)
    over5_switch_close_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_close_time_s_values)
    over5_switch_open_time_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_open_time_s_values)
    over5_switch_initially_closed_flags =
        deck_over16_boundary_plan === nothing ? Bool[] :
        copy(deck_over16_boundary_plan.over5_switch_initially_closed_flags)
    over5_switch_measuring_flags =
        deck_over16_boundary_plan === nothing ? Bool[] :
        copy(deck_over16_boundary_plan.over5_switch_measuring_flags)
    over5_switch_closed_markers =
        deck_over16_boundary_plan === nothing ? String[] :
        copy(deck_over16_boundary_plan.over5_switch_closed_markers)
    over5_switch_marker_texts =
        deck_over16_boundary_plan === nothing ? String[] :
        copy(deck_over16_boundary_plan.over5_switch_marker_texts)
    over5_switch_type_values =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over5_switch_type_values)
    over5_switch_critical_current_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_critical_current_values)
    over5_switch_random_opening_standard_deviation_s_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_random_opening_standard_deviation_s_values)
    over5_switch_on_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_on_conductance_values)
    over5_switch_off_conductance_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.over5_switch_off_conductance_values)
    over5_switch_output_codes =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.over5_switch_output_codes)
    source_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.source_names)
    source_node_names =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.source_node_names)
    source_node_values =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_node_values)
    source_iform_values =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_iform_values)
    source_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_line_numbers)
    source_layout_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.source_layout_kinds)
    source_tstart_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_tstart_values)
    source_tstop_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_tstop_values)
    source_crest_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_crest_values)
    source_time1_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_time1_values)
    source_time2_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_time2_values)
    source_sfreq_values =
        deck_over16_boundary_plan === nothing ? Float64[] :
        copy(deck_over16_boundary_plan.source_sfreq_values)
    source_card_kinds =
        deck_over16_boundary_plan === nothing ? Symbol[] :
        copy(deck_over16_boundary_plan.source_card_kinds)
    source_card_values =
        deck_over16_boundary_plan === nothing ? Vector{Float64}[] :
        [copy(values) for values in deck_over16_boundary_plan.source_card_values]
    source_card_provided_value_counts =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_card_provided_value_counts)
    source_card_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_card_line_numbers)
    source_interpolation_values =
        deck_over16_boundary_plan === nothing ? Vector{Float64}[] :
        [copy(values) for values in deck_over16_boundary_plan.source_interpolation_values]
    source_interpolation_provided_value_counts =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_interpolation_provided_value_counts)
    source_interpolation_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_interpolation_line_numbers)
    source_tacs_override_positions =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_tacs_override_positions)
    source_tacs_override_xtcs_indices =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_tacs_override_xtcs_indices)
    source_tacs_override_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_tacs_override_line_numbers)
    source_analytic_values =
        deck_over16_boundary_plan === nothing ? Vector{Float64}[] :
        [copy(values) for values in deck_over16_boundary_plan.source_analytic_values]
    source_analytic_provided_value_counts =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_analytic_provided_value_counts)
    source_analytic_line_numbers =
        deck_over16_boundary_plan === nothing ? Int[] :
        copy(deck_over16_boundary_plan.source_analytic_line_numbers)
    return (
        source = :emt_deck_over16_boundary_run,
        outcome = :timestep_integration,
        context = context,
        trace = trace,
        over16_updates = updates,
        sample_count = length(trace.time_s),
        over16_update_count = over16_update_count,
        accepted_timestep_state_mutation_count = mutation_count,
        over16_state_mutated = mutation_count > 0,
        legacy_fortran_in_loop = false,
        full_bpa_timestep_executed = false,
        full_deck_orchestration_executed = false,
        replacement_ready = false,
        deck_over16_boundary_plan = deck_over16_boundary_plan,
        deck_owned_over16_output_config = deck_owned_over16_output_config,
        deck_owned_over5a_source_config = deck_owned_over5a_source_config,
        deck_owned_over16_source_card_config = deck_owned_over16_source_card_config,
        deck_owned_over16_post_extrema_config = deck_owned_over16_post_extrema_config,
        deck_owned_time_switch_sparse_config = deck_owned_time_switch_sparse_config,
        steady_state_initial_sample_applied = steady_state_initial_sample !== nothing,
        steady_state_network_state_seeded = steady_state_initial_sample !== nothing,
        steady_state_initial_frequency_hz =
            steady_state_initial_sample === nothing ? 0.0 :
            Float64(steady_state_initial_sample.steady_state_frequency_hz),
        steady_state_initial_voltage_values =
            steady_state_initial_sample === nothing ? Float64[] :
            Float64.(steady_state_initial_sample.node_voltage_values),
        steady_state_initial_output_voltage_values =
            steady_state_initial_sample === nothing ? Float64[] :
            Float64.(steady_state_initial_sample.output_voltage_values),
        steady_state_initial_source_row_count =
            steady_state_initial_sample === nothing ? 0 :
            Int(steady_state_initial_sample.source_row_count),
        steady_state_initial_distributed_transposed_line_count =
            steady_state_initial_sample === nothing ? 0 :
            Int(steady_state_initial_sample.distributed_transposed_line_count),
        steady_state_initial_saturated_transformer_branch_count =
            steady_state_initial_sample === nothing ? 0 :
            Int(steady_state_initial_sample.saturated_transformer_branch_count),
        steady_state_initial_saturated_transformer_ideal_branch_count =
            steady_state_initial_sample === nothing ? 0 :
            Int(steady_state_initial_sample.saturated_transformer_ideal_branch_count),
        steady_state_initial_saturated_transformer_linearized_nonlinear_branch_count =
            steady_state_initial_sample === nothing ? 0 :
            Int(steady_state_initial_sample.saturated_transformer_linearized_nonlinear_branch_count),
        over16_source_update_mutation_count = source_update_mutation_count,
        over16_post_extrema_state_mutation_count = post_extrema_mutation_count,
        over16_tacs_linear_solve_state_mutation_count = tacs_linear_mutation_count,
        over16_tacs_post_solve_state_mutation_count = tacs_post_solve_mutation_count,
        over16_csup_state_mutation_count = csup_mutation_count,
        over16_switch_state_mutation_count = switch_state_mutation_count,
        over16_switch_topology_admittance_mutation_count =
            switch_topology_admittance_mutation_count,
        over16_switch_admittance_mutation_count = switch_admittance_mutation_count,
        over16_switch_retriangularization_mutation_count =
            switch_retriangularization_mutation_count,
        over16_switch_sparse_factor_workspace_mutation_count =
            switch_sparse_factor_workspace_mutation_count,
        over16_switch_fortran_sparse_factor_workspace_mutation_count =
            switch_fortran_sparse_factor_workspace_mutation_count,
        over16_switch_network_solution_mutation_count =
            switch_network_solution_mutation_count,
        over16_switch_current_mutation_count = switch_current_mutation_count,
        over16_switch_post_current_mutation_count = switch_post_current_mutation_count,
        over16_switch_bvalue_mutation_count = switch_bvalue_mutation_count,
        over16_tacs_execution_step_count = tacs_execution_count,
        over16_ntacs_execution_step_count = ntacs_execution_count,
        over16_elec_call_intent_step_count = elec_call_intent_count,
        over16_elec_output_copy_step_count = elec_output_copy_count,
        over16_elec_tacsto_interpreter_step_count = elec_tacsto_interpreter_count,
        over16_elec_report_step_count = elec_report_count,
        over16_elec_error_step_count = elec_error_count,
        over16_csup_execution_step_count = csup_execution_count,
        over16_machine_solution_pass_step_count = machine_solution_pass_count,
        over16_machine_equation_pass_step_count = machine_equation_pass_count,
        over16_solvum_execution_step_count = solvum_execution_count,
        over16_full_solvum_execution_step_count = full_solvum_execution_count,
        over16_complete_machine_network_solution_step_count =
            complete_machine_network_solution_count,
        over16_machine_terminal_network_coupling_step_count =
            machine_terminal_network_coupling_count,
        over16_machine_solution_pass_network_context_step_count =
            machine_solution_pass_network_context_count,
        over16_machine_solution_pass_output_network_coupling_step_count =
            machine_solution_pass_output_network_coupling_count,
        over16_machine_solution_pass_network_coupling_step_count =
            machine_solution_pass_network_coupling_count,
        over16_machine_runtime_projection_state_mutation_count =
            machine_runtime_projection_mutation_count,
        over16_machine_runtime_state_mutation_count = machine_runtime_state_mutation_count,
        over16_machine_runtime_storage_step_count = machine_runtime_storage_count,
        over16_machine_terminal_output_state_mutation_count =
            machine_terminal_output_mutation_count,
        over16_machine_terminal_network_current_mutation_count =
            machine_terminal_network_current_mutation_count,
        over16_machine_terminal_output_values = machine_terminal_output_values,
        over16_machine_terminal_current_substitution_values =
            machine_terminal_current_substitution_values,
        over16_machine_terminal_network_current_values =
            machine_terminal_network_current_values,
        over16_machine_runtime_dynamic_values = machine_runtime_dynamic_values,
        deck_over16_output_channel_names = output_channel_names,
        deck_over16_output_node_names = output_node_names,
        deck_over16_output_node_indices = output_node_indices,
        deck_over16_output_channel_line_numbers = output_channel_line_numbers,
        deck_over16_output_channel_count = length(output_channel_names),
        deck_over16_branch_voltage_output_names = branch_voltage_output_names,
        deck_over16_branch_voltage_branch_names = branch_voltage_branch_names,
        deck_over16_branch_voltage_branch_indices = branch_voltage_branch_indices,
        deck_over16_branch_voltage_output_line_numbers = branch_voltage_output_line_numbers,
        deck_over16_branch_voltage_output_count = length(branch_voltage_output_names),
        deck_over16_branch_current_output_names = branch_current_output_names,
        deck_over16_branch_current_branch_names = branch_current_branch_names,
        deck_over16_branch_current_branch_indices = branch_current_branch_indices,
        deck_over16_branch_current_output_line_numbers = branch_current_output_line_numbers,
        deck_over16_branch_current_output_count = length(branch_current_output_names),
        deck_over16_branch_power_output_names = branch_power_output_names,
        deck_over16_branch_power_branch_names = branch_power_branch_names,
        deck_over16_branch_power_branch_indices = branch_power_branch_indices,
        deck_over16_branch_power_output_line_numbers = branch_power_output_line_numbers,
        deck_over16_branch_power_output_count = length(branch_power_output_names),
        deck_over15_output_request_names = over15_output_request_names,
        deck_over15_output_request_output_kinds = over15_output_request_output_kinds,
        deck_over15_output_request_request_kinds = over15_output_request_request_kinds,
        deck_over15_output_request_layout_kinds = over15_output_request_layout_kinds,
        deck_over15_output_request_line_numbers = over15_output_request_line_numbers,
        deck_over15_output_request_output_codes = over15_output_request_output_codes,
        deck_over15_output_request_node_names = over15_output_request_node_names,
        deck_over15_output_request_node_indices = over15_output_request_node_indices,
        deck_over15_output_request_branch_names = over15_output_request_branch_names,
        deck_over15_output_request_branch_indices = over15_output_request_branch_indices,
        deck_over15_output_request_count = length(over15_output_request_names),
        deck_branch_names = branch_names,
        deck_branch_kinds = branch_kinds,
        deck_branch_from_node_names = branch_from_node_names,
        deck_branch_to_node_names = branch_to_node_names,
        deck_branch_from_node_indices = branch_from_node_indices,
        deck_branch_to_node_indices = branch_to_node_indices,
        deck_branch_conductance_values = branch_conductance_values,
        deck_branch_resistance_values = branch_resistance_values,
        deck_branch_inductance_values = branch_inductance_values,
        deck_branch_capacitance_values = branch_capacitance_values,
        deck_branch_previous_current_values = branch_previous_current_values,
        deck_branch_previous_voltage_values = branch_previous_voltage_values,
        deck_branch_line_numbers = branch_line_numbers,
        deck_branch_count = length(branch_names),
        deck_bergeron_line_names = bergeron_line_names,
        deck_bergeron_line_line_numbers = bergeron_line_line_numbers,
        deck_bergeron_line_from_node_names = bergeron_line_from_node_names,
        deck_bergeron_line_to_node_names = bergeron_line_to_node_names,
        deck_bergeron_line_from_node_indices = bergeron_line_from_node_indices,
        deck_bergeron_line_to_node_indices = bergeron_line_to_node_indices,
        deck_bergeron_line_surge_impedance_values = bergeron_line_surge_impedance_values,
        deck_bergeron_line_surge_admittance_values = bergeron_line_surge_admittance_values,
        deck_bergeron_line_travel_time_s_values = bergeron_line_travel_time_s_values,
        deck_bergeron_line_dt_s_values = bergeron_line_dt_s_values,
        deck_bergeron_line_attenuation_values = bergeron_line_attenuation_values,
        deck_bergeron_line_delay_step_counts = bergeron_line_delay_step_counts,
        deck_bergeron_line_write_indices = bergeron_line_write_indices,
        deck_bergeron_line_history_current_from_values = bergeron_line_history_current_from_values,
        deck_bergeron_line_history_current_to_values = bergeron_line_history_current_to_values,
        deck_bergeron_line_terminal_voltage_from_values =
            bergeron_line_terminal_voltage_from_values,
        deck_bergeron_line_terminal_voltage_to_values =
            bergeron_line_terminal_voltage_to_values,
        deck_bergeron_line_terminal_current_from_values =
            bergeron_line_terminal_current_from_values,
        deck_bergeron_line_terminal_current_to_values =
            bergeron_line_terminal_current_to_values,
        deck_bergeron_line_traveling_wave_from_values =
            bergeron_line_traveling_wave_from_values,
        deck_bergeron_line_traveling_wave_to_values =
            bergeron_line_traveling_wave_to_values,
        deck_bergeron_line_count = length(bergeron_line_names),
        deck_over2_branch_names = over2_branch_names,
        deck_over2_branch_line_numbers = over2_branch_line_numbers,
        deck_over2_branch_kinds = over2_branch_kinds,
        deck_over2_branch_layout_kinds = over2_branch_layout_kinds,
        deck_over2_branch_source_kinds = over2_branch_source_kinds,
        deck_over2_branch_reference_kinds = over2_branch_reference_kinds,
        deck_over2_branch_reference_names = over2_branch_reference_names,
        deck_over2_branch_reference_line_numbers = over2_branch_reference_line_numbers,
        deck_over2_branch_from_node_names = over2_branch_from_node_names,
        deck_over2_branch_to_node_names = over2_branch_to_node_names,
        deck_over2_branch_from_node_indices = over2_branch_from_node_indices,
        deck_over2_branch_to_node_indices = over2_branch_to_node_indices,
        deck_over2_branch_raw_resistance_values = over2_branch_raw_resistance_values,
        deck_over2_branch_raw_inductance_values = over2_branch_raw_inductance_values,
        deck_over2_branch_raw_capacitance_values = over2_branch_raw_capacitance_values,
        deck_over2_branch_conductance_values = over2_branch_conductance_values,
        deck_over2_branch_resistance_values = over2_branch_resistance_values,
        deck_over2_branch_inductance_values = over2_branch_inductance_values,
        deck_over2_branch_capacitance_values = over2_branch_capacitance_values,
        deck_over2_branch_output_codes = over2_branch_output_codes,
        deck_over2_branch_count = length(over2_branch_names),
        deck_time_switch_names = switch_names,
        deck_time_switch_line_numbers = switch_line_numbers,
        deck_time_switch_from_node_names = switch_from_node_names,
        deck_time_switch_to_node_names = switch_to_node_names,
        deck_time_switch_from_node_indices = switch_from_node_indices,
        deck_time_switch_to_node_indices = switch_to_node_indices,
        deck_time_switch_close_time_s_values = switch_close_time_s_values,
        deck_time_switch_open_time_s_values = switch_open_time_s_values,
        deck_time_switch_initially_closed_flags = switch_initially_closed_flags,
        deck_time_switch_on_conductance_values = switch_on_conductance_values,
        deck_time_switch_off_conductance_values = switch_off_conductance_values,
        deck_time_switch_count = length(switch_names),
        deck_over5_switch_names = over5_switch_names,
        deck_over5_switch_line_numbers = over5_switch_line_numbers,
        deck_over5_switch_from_node_names = over5_switch_from_node_names,
        deck_over5_switch_to_node_names = over5_switch_to_node_names,
        deck_over5_switch_from_node_indices = over5_switch_from_node_indices,
        deck_over5_switch_to_node_indices = over5_switch_to_node_indices,
        deck_over5_switch_layout_kinds = over5_switch_layout_kinds,
        deck_over5_switch_raw_close_time_s_values = over5_switch_raw_close_time_s_values,
        deck_over5_switch_raw_open_time_s_values = over5_switch_raw_open_time_s_values,
        deck_over5_switch_close_time_s_values = over5_switch_close_time_s_values,
        deck_over5_switch_open_time_s_values = over5_switch_open_time_s_values,
        deck_over5_switch_initially_closed_flags = over5_switch_initially_closed_flags,
        deck_over5_switch_measuring_flags = over5_switch_measuring_flags,
        deck_over5_switch_closed_markers = over5_switch_closed_markers,
        deck_over5_switch_marker_texts = over5_switch_marker_texts,
        deck_over5_switch_type_values = over5_switch_type_values,
        deck_over5_switch_critical_current_values =
            over5_switch_critical_current_values,
        deck_over5_switch_random_opening_standard_deviation_s_values =
            over5_switch_random_opening_standard_deviation_s_values,
        deck_over5_switch_on_conductance_values = over5_switch_on_conductance_values,
        deck_over5_switch_off_conductance_values = over5_switch_off_conductance_values,
        deck_over5_switch_output_codes = over5_switch_output_codes,
        deck_over5_switch_count = length(over5_switch_names),
        deck_over5a_source_names = source_names,
        deck_over5a_source_node_names = source_node_names,
        deck_over5a_source_node_values = source_node_values,
        deck_over5a_source_iform_values = source_iform_values,
        deck_over5a_source_line_numbers = source_line_numbers,
        deck_over5a_source_layout_kinds = source_layout_kinds,
        deck_over5a_source_tstart_values = source_tstart_values,
        deck_over5a_source_tstop_values = source_tstop_values,
        deck_over5a_source_crest_values = source_crest_values,
        deck_over5a_source_time1_values = source_time1_values,
        deck_over5a_source_time2_values = source_time2_values,
        deck_over5a_source_sfreq_values = source_sfreq_values,
        deck_over5a_source_count = length(source_names),
        deck_over16_source_card_kinds = source_card_kinds,
        deck_over16_source_card_values = source_card_values,
        deck_over16_source_card_provided_value_counts = source_card_provided_value_counts,
        deck_over16_source_card_line_numbers = source_card_line_numbers,
        deck_over16_source_card_count = length(source_card_kinds),
        deck_over16_source_interpolation_values = source_interpolation_values,
        deck_over16_source_interpolation_provided_value_counts =
            source_interpolation_provided_value_counts,
        deck_over16_source_interpolation_line_numbers = source_interpolation_line_numbers,
        deck_over16_source_interpolation_count = length(source_interpolation_values),
        deck_over16_source_tacs_override_positions = source_tacs_override_positions,
        deck_over16_source_tacs_override_xtcs_indices = source_tacs_override_xtcs_indices,
        deck_over16_source_tacs_override_line_numbers = source_tacs_override_line_numbers,
        deck_over16_source_tacs_override_count = length(source_tacs_override_positions),
        deck_over16_source_analytic_values = source_analytic_values,
        deck_over16_source_analytic_provided_value_counts =
            source_analytic_provided_value_counts,
        deck_over16_source_analytic_line_numbers = source_analytic_line_numbers,
        deck_over16_source_analytic_count = length(source_analytic_values),
        deferred_effects = (
            :full_init_step_calc_elec_orchestration,
            :full_bpa_deck_grammar,
            :solvum,
            :report_file_writers,
            :external_bpa_executable_waveform_comparison,
        ),
    )
end

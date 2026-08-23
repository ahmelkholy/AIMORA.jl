function _over16_check_length(name::String, values::AbstractVector, expected::Int)
    length(values) == expected || throw(ArgumentError("$name length must be $expected"))
    return nothing
end

function _over16_check_optional_length(name::String, values::AbstractVector, expected::Int)
    (isempty(values) || length(values) == expected) ||
        throw(ArgumentError("$name length must be 0 or $expected"))
    return nothing
end

function _over16_optional_float(values::AbstractVector{Float64}, index::Int, default::Float64)
    return isempty(values) ? default : values[index]
end

function _over16_switch_conductance_vector(
    name::String,
    values::AbstractVector{<:Real},
    switch_count::Int,
    default::Real,
)
    default_value = Float64(default)
    isfinite(default_value) && default_value >= 0.0 ||
        throw(ArgumentError("$(name) default must be finite and nonnegative"))
    if isempty(values)
        return fill(default_value, switch_count)
    end
    _over16_check_length(name, values, switch_count)
    result = Float64.(values)
    for value in result
        isfinite(value) && value >= 0.0 ||
            throw(ArgumentError("$(name) entries must be finite and nonnegative"))
    end
    return result
end

function _over16_optional_int(values::AbstractVector{Int}, index::Int, default::Int)
    return isempty(values) ? default : values[index]
end

function over16_controlled_switch_table_scan(
    switch_types::AbstractVector{Int},
    positions::AbstractVector{Int},
    elapsed_open_times::AbstractVector{Float64},
    maximum_open_times::AbstractVector{Float64},
    voltage_differences::AbstractVector{Float64},
    switch_currents::AbstractVector{Float64},
    close_thresholds::AbstractVector{Float64},
    open_thresholds::AbstractVector{Float64},
    dt::Float64;
    controlled_nodes::AbstractVector{Int}=Int[],
    control_indices::AbstractVector{Int}=Int[],
    clamp_signals::AbstractVector{Float64}=Float64[],
    close_control_signals::AbstractVector{Float64}=Float64[],
    gap_control_signals::AbstractVector{Float64}=Float64[],
    previous_gap_currents::AbstractVector{Float64}=Float64[],
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
    gap_storage_limit::Int=100,
)
    row_count = length(switch_types)
    _over16_check_length("positions", positions, row_count)
    _over16_check_length("elapsed_open_times", elapsed_open_times, row_count)
    _over16_check_length("maximum_open_times", maximum_open_times, row_count)
    _over16_check_length("voltage_differences", voltage_differences, row_count)
    _over16_check_length("switch_currents", switch_currents, row_count)
    _over16_check_length("close_thresholds", close_thresholds, row_count)
    _over16_check_length("open_thresholds", open_thresholds, row_count)
    _over16_check_optional_length("controlled_nodes", controlled_nodes, row_count)
    _over16_check_optional_length("control_indices", control_indices, row_count)
    _over16_check_optional_length("clamp_signals", clamp_signals, row_count)
    _over16_check_optional_length("close_control_signals", close_control_signals, row_count)
    _over16_check_optional_length("gap_control_signals", gap_control_signals, row_count)
    gap_storage_limit >= 0 || throw(ArgumentError("gap_storage_limit must be nonnegative"))

    updated_positions = collect(positions)
    updated_elapsed = collect(elapsed_open_times)
    updated_gap_currents = collect(previous_gap_currents)
    modswt = Int[]
    scan_actions = Vector{Symbol}(undef, row_count)
    transitions = Vector{Symbol}(undef, row_count)
    gap_count = 0

    for row in 1:row_count
        switch_type = switch_types[row]
        if switch_type == 8890
            gap_count += 1
            gap_count <= gap_storage_limit ||
                throw(ArgumentError("OVER16 controlled-switch gap storage overflow"))
            if length(updated_gap_currents) < gap_count
                push!(updated_gap_currents, 0.0)
            end
        end

        if !(switch_type == 8888 || switch_type == 8890 || switch_type == 8891)
            scan_actions[row] = :skip
            transitions[row] = :none
            continue
        end

        previous_gap_current = switch_type == 8890 ? updated_gap_currents[gap_count] : 0.0
        step = over16_controlled_switch_scan_step(
            switch_type,
            updated_positions[row],
            updated_elapsed[row],
            maximum_open_times[row],
            voltage_differences[row],
            switch_currents[row],
            close_thresholds[row],
            open_thresholds[row],
            dt;
            controlled_node = _over16_optional_int(controlled_nodes, row, 0),
            control_index = _over16_optional_int(control_indices, row, 0),
            clamp_signal = _over16_optional_float(clamp_signals, row, 0.0),
            close_control_signal = _over16_optional_float(close_control_signals, row, 0.0),
            gap_control_signal = _over16_optional_float(gap_control_signals, row, 0.0),
            previous_gap_current = previous_gap_current,
            flzero = flzero,
            fltinf = fltinf,
        )

        scan_actions[row] = step.scan_action
        transitions[row] = step.transition
        updated_positions[row] = step.position
        updated_elapsed[row] = step.elapsed_open_time
        if switch_type == 8890
            updated_gap_currents[gap_count] = step.gap_current
        end
        if step.altered
            push!(modswt, step.modswt_sign * row)
        end
    end

    return (
        positions = updated_positions,
        elapsed_open_times = updated_elapsed,
        gap_currents = updated_gap_currents,
        modswt = modswt,
        ktrlsw_count = length(modswt),
        altered = !isempty(modswt),
        gap_count = gap_count,
        scan_actions = scan_actions,
        transitions = transitions,
    )
end

function over16_controlled_switch_table_scan!(
    state::OVER16SwitchScanState,
    switch_types::AbstractVector{Int},
    maximum_open_times::AbstractVector{Float64},
    voltage_differences::AbstractVector{Float64},
    switch_currents::AbstractVector{Float64},
    close_thresholds::AbstractVector{Float64},
    open_thresholds::AbstractVector{Float64},
    dt::Float64;
    controlled_nodes::AbstractVector{Int}=Int[],
    control_indices::AbstractVector{Int}=Int[],
    clamp_signals::AbstractVector{Float64}=Float64[],
    close_control_signals::AbstractVector{Float64}=Float64[],
    gap_control_signals::AbstractVector{Float64}=Float64[],
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
    gap_storage_limit::Int=100,
)
    positions_before = copy(state.positions)
    elapsed_before = copy(state.elapsed_open_times)
    gap_currents_before = copy(state.gap_currents)
    modswt_before = copy(state.modswt)
    preview = over16_controlled_switch_table_scan(
        switch_types,
        state.positions,
        state.elapsed_open_times,
        maximum_open_times,
        voltage_differences,
        switch_currents,
        close_thresholds,
        open_thresholds,
        dt;
        controlled_nodes = controlled_nodes,
        control_indices = control_indices,
        clamp_signals = clamp_signals,
        close_control_signals = close_control_signals,
        gap_control_signals = gap_control_signals,
        previous_gap_currents = state.gap_currents,
        flzero = flzero,
        fltinf = fltinf,
        gap_storage_limit = gap_storage_limit,
    )

    state.positions .= preview.positions
    state.elapsed_open_times .= preview.elapsed_open_times
    resize!(state.gap_currents, length(preview.gap_currents))
    state.gap_currents .= preview.gap_currents
    empty!(state.modswt)
    append!(state.modswt, preview.modswt)

    positions_mutated = state.positions != positions_before
    elapsed_open_times_mutated = state.elapsed_open_times != elapsed_before
    gap_currents_mutated = state.gap_currents != gap_currents_before
    modswt_mutated = state.modswt != modswt_before
    switch_scan_state_mutated =
        positions_mutated || elapsed_open_times_mutated ||
        gap_currents_mutated || modswt_mutated
    return merge(
        preview,
        (
            positions = copy(state.positions),
            elapsed_open_times = copy(state.elapsed_open_times),
            gap_currents = copy(state.gap_currents),
            modswt = copy(state.modswt),
            positions_mutated = positions_mutated,
            elapsed_open_times_mutated = elapsed_open_times_mutated,
            gap_currents_mutated = gap_currents_mutated,
            modswt_mutated = modswt_mutated,
            switch_scan_state_mutated = switch_scan_state_mutated,
            switch_graph_state_mutated = switch_scan_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_operation_schedule(
    modswt::AbstractVector{Int},
    closed_switch_count::Int,
    accumulated_operation_count::Int=0,
)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))

    opening_rows = Int[]
    closing_rows = Int[]
    for entry in modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        if entry < 0
            push!(opening_rows, -entry)
        else
            push!(closing_rows, entry)
        end
    end

    operation_count = length(modswt)
    updated_closed_count = closed_switch_count - length(opening_rows) + length(closing_rows)
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    return (
        opening_rows = opening_rows,
        closing_rows = closing_rows,
        processed_modswt = vcat(-opening_rows, closing_rows),
        ktrlsw_count = operation_count,
        closed_switch_count = updated_closed_count,
        accumulated_operation_count = accumulated_operation_count + operation_count,
        should_call_switch = operation_count > 0,
    )
end

function over16_switch_operation_schedule!(state::OVER16SwitchOperationState)
    modswt_before = copy(state.modswt)
    closed_count_before = state.closed_switch_count
    accumulated_count_before = state.accumulated_operation_count
    preview = over16_switch_operation_schedule(
        state.modswt,
        state.closed_switch_count,
        state.accumulated_operation_count,
    )

    empty!(state.modswt)
    append!(state.modswt, preview.processed_modswt)
    state.closed_switch_count = preview.closed_switch_count
    state.accumulated_operation_count = preview.accumulated_operation_count

    modswt_mutated = state.modswt != modswt_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    accumulated_operation_count_mutated =
        state.accumulated_operation_count != accumulated_count_before
    switch_operation_state_mutated =
        modswt_mutated || closed_switch_count_mutated ||
        accumulated_operation_count_mutated
    return merge(
        preview,
        (
            processed_modswt = copy(state.modswt),
            closed_switch_count = state.closed_switch_count,
            accumulated_operation_count = state.accumulated_operation_count,
            modswt_mutated = modswt_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            accumulated_operation_count_mutated = accumulated_operation_count_mutated,
            switch_operation_state_mutated = switch_operation_state_mutated,
            switch_graph_state_mutated = switch_operation_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_operation_schedule_lean!(state::OVER16SwitchOperationState)
    state.closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    state.accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))

    opening_count = 0
    closing_count = 0
    @inbounds for entry in state.modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        if entry < 0
            opening_count += 1
        else
            closing_count += 1
        end
    end
    operation_count = opening_count + closing_count
    updated_closed_count =
        state.closed_switch_count - opening_count + closing_count
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))

    # Stable in-place partition: the reference schedule processes openings before
    # closings while retaining the original order within each group.
    next_opening = firstindex(state.modswt)
    @inbounds for source in eachindex(state.modswt)
        if state.modswt[source] < 0
            entry = state.modswt[source]
            for destination in source:-1:(next_opening + 1)
                state.modswt[destination] = state.modswt[destination - 1]
            end
            state.modswt[next_opening] = entry
            next_opening += 1
        end
    end

    state.closed_switch_count = updated_closed_count
    state.accumulated_operation_count += operation_count
    return SwitchOperationStepResult(
        state.modswt,
        operation_count,
        operation_count > 0,
    )
end

function over16_switch_status_update(
    modswt::AbstractVector{Int},
    closed_mask::AbstractVector{Bool},
    closed_switch_count::Int,
    first_group_head::Int=0;
    strict_consistency::Bool=true,
)
    switch_count = length(closed_mask)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    closed_switch_count <= switch_count ||
        throw(ArgumentError("closed_switch_count cannot exceed switch count"))
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    if strict_consistency && count(identity, closed_mask) != closed_switch_count
        throw(ArgumentError("closed_switch_count must match closed_mask in strict mode"))
    end

    updated_closed_mask = collect(closed_mask)
    updated_closed_count = closed_switch_count
    for entry in modswt
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        row = abs(entry)
        1 <= row <= switch_count ||
            throw(ArgumentError("MODSWT row index out of switch range"))
        if entry < 0
            if strict_consistency && !updated_closed_mask[row]
                throw(ArgumentError("cannot open an already-open switch in strict mode"))
            end
            updated_closed_count -= 1
            updated_closed_mask[row] = false
        else
            if strict_consistency && updated_closed_mask[row]
                throw(ArgumentError("cannot close an already-closed switch in strict mode"))
            end
            updated_closed_count += 1
            updated_closed_mask[row] = true
        end
    end
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    updated_closed_count <= switch_count ||
        throw(ArgumentError("switch operations produce a closed-switch count above switch count"))

    return (
        closed_mask = updated_closed_mask,
        closed_switch_count = updated_closed_count,
        first_group_head = updated_closed_count > 0 ? first_group_head : 0,
        clear_nextsw = updated_closed_count == 0,
        clear_kode = updated_closed_count == 0,
        requires_order_rebuild = !isempty(modswt) && updated_closed_count > 0,
        topology_mutated = false,
    )
end

function over16_switch_status_update!(
    state::OVER16SwitchTopologyState,
    modswt::AbstractVector{Int};
    strict_consistency::Bool=true,
)
    closed_mask_before = copy(state.closed_mask)
    closed_count_before = state.closed_switch_count
    first_group_head_before = state.first_group_head
    nextsw_before = copy(state.nextsw)
    kode_before = copy(state.kode)
    preview = over16_switch_status_update(
        modswt,
        state.closed_mask,
        state.closed_switch_count,
        state.first_group_head;
        strict_consistency = strict_consistency,
    )

    state.closed_mask .= preview.closed_mask
    state.closed_switch_count = preview.closed_switch_count
    state.first_group_head = preview.first_group_head
    if preview.clear_nextsw
        fill!(state.nextsw, 0)
    end
    if preview.clear_kode
        fill!(state.kode, 0)
    end

    closed_mask_mutated = state.closed_mask != closed_mask_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    first_group_head_mutated = state.first_group_head != first_group_head_before
    nextsw_cleared = state.nextsw != nextsw_before
    kode_cleared = state.kode != kode_before
    switch_status_state_mutated =
        closed_mask_mutated || closed_switch_count_mutated ||
        first_group_head_mutated || nextsw_cleared || kode_cleared
    return merge(
        preview,
        (
            closed_mask = copy(state.closed_mask),
            closed_switch_count = state.closed_switch_count,
            first_group_head = state.first_group_head,
            nextsw = copy(state.nextsw),
            kode = copy(state.kode),
            closed_mask_mutated = closed_mask_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            first_group_head_mutated = first_group_head_mutated,
            nextsw_cleared = nextsw_cleared,
            kode_cleared = kode_cleared,
            switch_status_state_mutated = switch_status_state_mutated,
            switch_graph_state_mutated = switch_status_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_status_update_lean!(
    state::OVER16SwitchTopologyState,
    modswt::AbstractVector{Int};
    strict_consistency::Bool=true,
)
    switch_count = length(state.closed_mask)
    state.closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    state.closed_switch_count <= switch_count ||
        throw(ArgumentError("closed_switch_count cannot exceed switch count"))
    0 <= state.first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    if strict_consistency && count(identity, state.closed_mask) != state.closed_switch_count
        throw(ArgumentError("closed_switch_count must match closed_mask in strict mode"))
    end

    updated_closed_count = state.closed_switch_count
    @inbounds for operation_index in eachindex(modswt)
        entry = modswt[operation_index]
        entry != 0 || throw(ArgumentError("MODSWT entries must be nonzero"))
        row = abs(entry)
        1 <= row <= switch_count ||
            throw(ArgumentError("MODSWT row index out of switch range"))
        row_closed = state.closed_mask[row]
        for previous_index in firstindex(modswt):(operation_index - 1)
            previous_entry = modswt[previous_index]
            abs(previous_entry) == row && (row_closed = previous_entry > 0)
        end
        if entry < 0
            strict_consistency && !row_closed &&
                throw(ArgumentError("cannot open an already-open switch in strict mode"))
            updated_closed_count -= 1
        else
            strict_consistency && row_closed &&
                throw(ArgumentError("cannot close an already-closed switch in strict mode"))
            updated_closed_count += 1
        end
    end
    updated_closed_count >= 0 ||
        throw(ArgumentError("switch operations produce a negative closed-switch count"))
    updated_closed_count <= switch_count ||
        throw(ArgumentError("switch operations produce a closed-switch count above switch count"))

    closed_mask_mutated = false
    @inbounds for row in eachindex(state.closed_mask)
        updated = state.closed_mask[row]
        for entry in modswt
            abs(entry) == row && (updated = entry > 0)
        end
        closed_mask_mutated |= updated != state.closed_mask[row]
        state.closed_mask[row] = updated
    end
    closed_switch_count_mutated = updated_closed_count != state.closed_switch_count
    state.closed_switch_count = updated_closed_count

    first_group_head_mutated = false
    nextsw_cleared = false
    kode_cleared = false
    if updated_closed_count == 0
        first_group_head_mutated = state.first_group_head != 0
        state.first_group_head = 0
        @inbounds for index in eachindex(state.nextsw)
            nextsw_cleared |= state.nextsw[index] != 0
            state.nextsw[index] = 0
        end
        @inbounds for index in eachindex(state.kode)
            kode_cleared |= state.kode[index] != 0
            state.kode[index] = 0
        end
    end
    mutated =
        closed_mask_mutated || closed_switch_count_mutated ||
        first_group_head_mutated || nextsw_cleared || kode_cleared
    return SwitchStatusStepResult(
        !isempty(modswt) && updated_closed_count > 0,
        mutated,
    )
end

function _over16_endpoint_is_known(node::Int, partition_boundary::Int, reference_node::Int)
    return node == reference_node || node > partition_boundary
end

function _over16_shares_endpoint(
    a_from::Int,
    a_to::Int,
    b_from::Int,
    b_to::Int,
)
    return a_from == b_from || a_from == b_to || a_to == b_from || a_to == b_to
end

function _over16_insert_kode_chain!(
    kode::Vector{Int},
    from_node::Int,
    to_node::Int,
    node_count::Int,
)
    if kode[from_node] == 0 && kode[to_node] == 0
        kode[from_node] = to_node
        kode[to_node] = from_node
        return nothing
    end

    anchor = from_node
    if kode[from_node] == 0
        anchor = to_node
    end
    insert_node = anchor == from_node ? to_node : from_node
    first_insert = insert_node

    steps = 0
    while kode[anchor] >= anchor
        anchor = kode[anchor]
        steps += 1
        steps <= node_count ||
            throw(ArgumentError("KODE chain search did not terminate"))
    end

    while true
        next_insert = kode[insert_node]
        if !(insert_node <= anchor && insert_node >= kode[anchor])
            old_next = kode[anchor]
            kode[anchor] = insert_node
            kode[insert_node] = old_next
            if insert_node > anchor
                anchor = insert_node
            end
        else
            predecessor = kode[anchor]
            steps = 0
            while !(insert_node < kode[predecessor])
                predecessor = kode[predecessor]
                steps += 1
                steps <= node_count ||
                    throw(ArgumentError("KODE insertion search did not terminate"))
            end
            old_next = kode[predecessor]
            kode[predecessor] = insert_node
            kode[insert_node] = old_next
        end

        if next_insert == 0 || next_insert == first_insert
            return nothing
        end
        insert_node = next_insert
    end
end

function over16_switch_simple_ordering(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool},
    partition_boundary::Int;
    expected_closed_switch_count::Union{Nothing,Int}=nothing,
    node_count::Int=max(partition_boundary, maximum(vcat([1], collect(from_nodes), collect(to_nodes)))),
    reference_node::Int=1,
)
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    partition_boundary >= 1 ||
        throw(ArgumentError("partition_boundary must be positive"))
    reference_node >= 1 ||
        throw(ArgumentError("reference_node must be positive"))
    node_count >= partition_boundary ||
        throw(ArgumentError("node_count must be at least partition_boundary"))

    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        1 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be within node_count"))
        1 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be within node_count"))
        from_node != to_node ||
            throw(ArgumentError("closed switch ordering requires distinct endpoints"))
    end

    closed_switch_count = count(identity, closed_mask)
    if expected_closed_switch_count !== nothing &&
       closed_switch_count != expected_closed_switch_count
        throw(ArgumentError("expected_closed_switch_count must match closed_mask"))
    end
    nextsw = zeros(Int, switch_count)
    kode = zeros(Int, node_count)
    kcl_endpoint_indices = zeros(Int, switch_count)
    if closed_switch_count == 0
        return (
            nextsw = nextsw,
            kode = kode,
            first_group_head = 0,
            closed_switch_count = 0,
            ordered_switch_count = 0,
            ordered_rows = Int[],
            kcl_endpoint_indices = kcl_endpoint_indices,
            pass_count = 0,
            topology_mutated = false,
        )
    end

    ordered_switch_count = 0
    previous_count = 0
    previous_switch = 0
    first_switch = 0
    pass_count = 0
    ordered_rows = Int[]

    while ordered_switch_count < closed_switch_count
        for row in 1:switch_count
            if !closed_mask[row] || nextsw[row] != 0
                continue
            end

            from_node = from_nodes[row]
            to_node = to_nodes[row]
            from_blocked =
                _over16_endpoint_is_known(from_node, partition_boundary, reference_node)
            to_blocked =
                _over16_endpoint_is_known(to_node, partition_boundary, reference_node)
            blocked_count = (from_blocked ? 1 : 0) + (to_blocked ? 1 : 0)

            if blocked_count < 2
                for neighbor in 1:switch_count
                    if neighbor == row || !closed_mask[neighbor] || nextsw[neighbor] != 0
                        continue
                    end
                    if !_over16_shares_endpoint(
                        from_node,
                        to_node,
                        from_nodes[neighbor],
                        to_nodes[neighbor],
                    )
                        continue
                    end

                    if !from_blocked &&
                       (from_node == from_nodes[neighbor] || from_node == to_nodes[neighbor])
                        from_blocked = true
                        blocked_count += 1
                    elseif !to_blocked &&
                           (to_node == from_nodes[neighbor] || to_node == to_nodes[neighbor])
                        to_blocked = true
                        blocked_count += 1
                    end

                    if blocked_count == 2
                        break
                    end
                end
            end
            if blocked_count == 2
                continue
            end

            if ordered_switch_count > 0
                nextsw[previous_switch] *= row
            end
            previous_switch = row
            nextsw[row] = from_blocked ? -1 : 1
            kcl_endpoint_indices[row] = from_blocked ? 2 : 1
            ordered_switch_count += 1
            push!(ordered_rows, row)
            if ordered_switch_count == 1
                first_switch = row
            end
            _over16_insert_kode_chain!(kode, from_node, to_node, node_count)
        end

        pass_count += 1
        if ordered_switch_count <= previous_count
            throw(ArgumentError("closed switch ordering did not progress"))
        end
        previous_count = ordered_switch_count
    end

    nextsw[previous_switch] *= first_switch
    return (
        nextsw = nextsw,
        kode = kode,
        first_group_head = first_switch,
        closed_switch_count = closed_switch_count,
        ordered_switch_count = ordered_switch_count,
        ordered_rows = ordered_rows,
        kcl_endpoint_indices = kcl_endpoint_indices,
        pass_count = pass_count,
        topology_mutated = false,
    )
end

function over16_switch_simple_ordering!(
    state::OVER16SwitchTopologyState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    node_count::Int=length(state.kode),
    reference_node::Int=1,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    nextsw_before = copy(state.nextsw)
    kode_before = copy(state.kode)
    first_group_head_before = state.first_group_head
    closed_count_before = state.closed_switch_count
    preview = over16_switch_simple_ordering(
        from_nodes,
        to_nodes,
        state.closed_mask,
        partition_boundary;
        expected_closed_switch_count = state.closed_switch_count,
        node_count = node_count,
        reference_node = reference_node,
    )

    resize!(state.nextsw, length(preview.nextsw))
    state.nextsw .= preview.nextsw
    resize!(state.kode, length(preview.kode))
    state.kode .= preview.kode
    state.first_group_head = preview.first_group_head
    state.closed_switch_count = preview.closed_switch_count

    nextsw_mutated = state.nextsw != nextsw_before
    kode_mutated = state.kode != kode_before
    first_group_head_mutated = state.first_group_head != first_group_head_before
    closed_switch_count_mutated = state.closed_switch_count != closed_count_before
    switch_order_state_mutated =
        nextsw_mutated || kode_mutated || first_group_head_mutated ||
        closed_switch_count_mutated
    return merge(
        preview,
        (
            nextsw = copy(state.nextsw),
            kode = copy(state.kode),
            first_group_head = state.first_group_head,
            closed_switch_count = state.closed_switch_count,
            nextsw_mutated = nextsw_mutated,
            kode_mutated = kode_mutated,
            first_group_head_mutated = first_group_head_mutated,
            closed_switch_count_mutated = closed_switch_count_mutated,
            switch_order_state_mutated = switch_order_state_mutated,
            switch_graph_state_mutated = switch_order_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_simple_ordering_lean!(
    state::OVER16SwitchTopologyState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    node_count::Int=length(state.kode),
    reference_node::Int=1,
)
    result = over16_switch_simple_ordering!(
        state,
        from_nodes,
        to_nodes,
        partition_boundary;
        node_count = node_count,
        reference_node = reference_node,
    )
    return SwitchOrderStepResult(result.switch_order_state_mutated)
end

function over16_switch_admittance_update(
    base_admittance::AbstractMatrix{<:Real},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    node_count = size(base_admittance, 1)
    node_count > 0 || throw(ArgumentError("base_admittance must not be empty"))
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    closed_values = _over16_switch_conductance_vector(
        "closed_conductances",
        closed_conductances,
        switch_count,
        closed_conductance,
    )
    open_values = _over16_switch_conductance_vector(
        "open_conductances",
        open_conductances,
        switch_count,
        open_conductance,
    )

    admittance = Float64.(base_admittance)
    for value in admittance
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    switch_conductances = zeros(Float64, switch_count)
    closed_rows = Int[]
    open_rows = Int[]
    stamped_rows = Int[]
    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        0 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be between 0 and node count"))
        0 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be between 0 and node count"))
        from_node != to_node ||
            throw(ArgumentError("switch admittance endpoints must be distinct"))
        from_node != 0 || to_node != 0 ||
            throw(ArgumentError("switch admittance cannot connect ground to ground"))

        is_closed = closed_mask[row]
        conductance = is_closed ? closed_values[row] : open_values[row]
        switch_conductances[row] = conductance
        push!(is_closed ? closed_rows : open_rows, row)
        if conductance != 0.0
            stamp_conductance!(admittance, from_node, to_node, conductance)
            push!(stamped_rows, row)
        end
    end

    deferred_calls = Symbol[]
    if request_retriangularization
        push!(deferred_calls, :last14)
        push!(deferred_calls, :sparse_factor_update)
        push!(deferred_calls, :retriangularization_execution)
    end

    return (
        source = :over16_switch_admittance_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SWITCH, :LAST14),
        fortran_labels = (312, 317, 601, 1010, 3111),
        node_count = node_count,
        switch_count = switch_count,
        closed_switch_count = length(closed_rows),
        closed_rows = closed_rows,
        open_rows = open_rows,
        stamped_rows = stamped_rows,
        base_admittance = Float64.(base_admittance),
        admittance = admittance,
        switch_conductances = switch_conductances,
        should_retriangularize = request_retriangularization,
        should_clear_factor_workspace = request_retriangularization,
        deferred_calls = deferred_calls,
        topology_mutated = false,
        admittance_mutated = false,
        sparse_factor_mutated = false,
        retriangularized = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_admittance_update(
    node_count::Int,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    kwargs...,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return over16_switch_admittance_update(
        zeros(Float64, node_count, node_count),
        from_nodes,
        to_nodes,
        closed_mask;
        kwargs...,
    )
end

function over16_switch_admittance_update!(
    state::OVER16SwitchAdmittanceState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    admittance_before = copy(state.admittance)
    switch_conductances_before = copy(state.switch_conductances)
    retriangularization_count_before = state.retriangularization_count
    preview = over16_switch_admittance_update(
        state.base_admittance,
        from_nodes,
        to_nodes,
        closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )

    state.admittance .= preview.admittance
    resize!(state.switch_conductances, length(preview.switch_conductances))
    state.switch_conductances .= preview.switch_conductances
    if preview.should_retriangularize
        state.retriangularization_count += 1
    end

    admittance_mutated = state.admittance != admittance_before
    switch_conductances_mutated =
        state.switch_conductances != switch_conductances_before
    retriangularization_count_mutated =
        state.retriangularization_count != retriangularization_count_before
    switch_admittance_state_mutated =
        admittance_mutated || switch_conductances_mutated ||
        retriangularization_count_mutated
    return merge(
        preview,
        (
            base_admittance = copy(state.base_admittance),
            admittance = copy(state.admittance),
            switch_conductances = copy(state.switch_conductances),
            retriangularization_count = state.retriangularization_count,
            admittance_mutated = admittance_mutated,
            switch_conductances_mutated = switch_conductances_mutated,
            retriangularization_count_mutated = retriangularization_count_mutated,
            switch_admittance_state_mutated = switch_admittance_state_mutated,
            topology_mutated = false,
            sparse_factor_mutated = false,
            retriangularized = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function switch_admittance_update_lean!(
    state::OVER16SwitchAdmittanceState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    closed_mask::AbstractVector{Bool};
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=true,
)
    base_admittance = state.base_admittance
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    node_count = size(base_admittance, 1)
    node_count > 0 || throw(ArgumentError("base_admittance must not be empty"))
    size(state.admittance) == size(base_admittance) ||
        throw(ArgumentError("admittance and base_admittance sizes must match"))
    switch_count = length(closed_mask)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_optional_length(
        "closed_conductances",
        closed_conductances,
        switch_count,
    )
    _over16_check_optional_length(
        "open_conductances",
        open_conductances,
        switch_count,
    )
    closed_default = Float64(closed_conductance)
    open_default = Float64(open_conductance)
    isfinite(closed_default) && closed_default >= 0.0 ||
        throw(ArgumentError("closed_conductances default must be finite and nonnegative"))
    isfinite(open_default) && open_default >= 0.0 ||
        throw(ArgumentError("open_conductances default must be finite and nonnegative"))

    for value in base_admittance
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    @inbounds for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        0 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be between 0 and node count"))
        0 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be between 0 and node count"))
        from_node != to_node ||
            throw(ArgumentError("switch admittance endpoints must be distinct"))
        from_node != 0 || to_node != 0 ||
            throw(ArgumentError("switch admittance cannot connect ground to ground"))
        closed_value = isempty(closed_conductances) ?
            closed_default : Float64(closed_conductances[row])
        open_value = isempty(open_conductances) ?
            open_default : Float64(open_conductances[row])
        isfinite(closed_value) && closed_value >= 0.0 ||
            throw(ArgumentError("closed_conductances entries must be finite and nonnegative"))
        isfinite(open_value) && open_value >= 0.0 ||
            throw(ArgumentError("open_conductances entries must be finite and nonnegative"))
    end

    old_conductance_count = length(state.switch_conductances)
    switch_conductances_mutated = old_conductance_count != switch_count
    resize!(state.switch_conductances, switch_count)
    @inbounds for row in 1:switch_count
        conductance = closed_mask[row] ?
            (isempty(closed_conductances) ? closed_default : Float64(closed_conductances[row])) :
            (isempty(open_conductances) ? open_default : Float64(open_conductances[row]))
        if row <= old_conductance_count
            switch_conductances_mutated |=
                state.switch_conductances[row] != conductance
        end
        state.switch_conductances[row] = conductance
    end

    if size(state.admittance_workspace) != size(base_admittance)
        state.admittance_workspace = similar(base_admittance)
    end
    copyto!(state.admittance_workspace, base_admittance)
    @inbounds for row in 1:switch_count
        conductance = state.switch_conductances[row]
        conductance == 0.0 && continue
        stamp_conductance!(
            state.admittance_workspace,
            from_nodes[row],
            to_nodes[row],
            conductance,
        )
    end
    admittance_mutated = state.admittance != state.admittance_workspace
    copyto!(state.admittance, state.admittance_workspace)

    if request_retriangularization
        state.retriangularization_count += 1
    end
    state_mutated =
        admittance_mutated || switch_conductances_mutated ||
        request_retriangularization
    return SwitchAdmittanceStepResult(
        request_retriangularization,
        admittance_mutated,
        state_mutated,
    )
end

function over16_switch_topology_admittance_update!(
    topology_state::OVER16SwitchTopologyState,
    admittance_state::OVER16SwitchAdmittanceState,
    modswt::AbstractVector{Int},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    strict_consistency::Bool=true,
    node_count::Int=size(admittance_state.base_admittance, 1),
    reference_node::Int=1,
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=!isempty(modswt),
)
    status_result = over16_switch_status_update!(
        topology_state,
        modswt;
        strict_consistency = strict_consistency,
    )
    should_rebuild_order =
        status_result.requires_order_rebuild ||
        (topology_state.closed_switch_count > 0 && all(entry -> entry == 0, topology_state.nextsw))
    order_result = nothing
    if should_rebuild_order
        order_result = over16_switch_simple_ordering!(
            topology_state,
            from_nodes,
            to_nodes,
            partition_boundary;
            node_count = node_count,
            reference_node = reference_node,
        )
    end
    admittance_result = over16_switch_admittance_update!(
        admittance_state,
        from_nodes,
        to_nodes,
        topology_state.closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )

    switch_topology_state_mutated =
        status_result.switch_status_state_mutated ||
        (order_result !== nothing && order_result.switch_order_state_mutated)
    switch_admittance_state_mutated =
        admittance_result.switch_admittance_state_mutated
    return (
        source = :over16_switch_topology_admittance_update,
        outcome = :integration_boundary,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SWITCH, :LAST14),
        fortran_labels = (312, 317, 601, 1010, 3111),
        status_result = status_result,
        order_result = order_result,
        admittance_result = admittance_result,
        closed_mask = copy(topology_state.closed_mask),
        closed_switch_count = topology_state.closed_switch_count,
        first_group_head = topology_state.first_group_head,
        nextsw = copy(topology_state.nextsw),
        kode = copy(topology_state.kode),
        admittance = copy(admittance_state.admittance),
        switch_conductances = copy(admittance_state.switch_conductances),
        retriangularization_count = admittance_state.retriangularization_count,
        switch_topology_state_mutated = switch_topology_state_mutated,
        switch_admittance_state_mutated = switch_admittance_state_mutated,
        switch_topology_admittance_state_mutated =
            switch_topology_state_mutated || switch_admittance_state_mutated,
        topology_mutated = switch_topology_state_mutated,
        admittance_mutated = admittance_result.admittance_mutated,
        sparse_factor_mutated = false,
        should_retriangularize = admittance_result.should_retriangularize,
        retriangularized = false,
        deferred_calls = admittance_result.deferred_calls,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function switch_topology_admittance_update_lean!(
    topology_state::OVER16SwitchTopologyState,
    admittance_state::OVER16SwitchAdmittanceState,
    modswt::AbstractVector{Int},
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    partition_boundary::Int;
    strict_consistency::Bool=true,
    node_count::Int=size(admittance_state.base_admittance, 1),
    reference_node::Int=1,
    closed_conductances::AbstractVector{<:Real}=Float64[],
    open_conductances::AbstractVector{<:Real}=Float64[],
    closed_conductance::Real=1.0e9,
    open_conductance::Real=0.0,
    request_retriangularization::Bool=!isempty(modswt),
)
    status_result = switch_status_update_lean!(
        topology_state,
        modswt;
        strict_consistency = strict_consistency,
    )
    should_rebuild_order =
        status_result.requires_order_rebuild ||
        (topology_state.closed_switch_count > 0 &&
         all(iszero, topology_state.nextsw))
    order_result = should_rebuild_order ?
        switch_simple_ordering_lean!(
            topology_state,
            from_nodes,
            to_nodes,
            partition_boundary;
            node_count = node_count,
            reference_node = reference_node,
        ) : nothing
    admittance_result = switch_admittance_update_lean!(
        admittance_state,
        from_nodes,
        to_nodes,
        topology_state.closed_mask;
        closed_conductances = closed_conductances,
        open_conductances = open_conductances,
        closed_conductance = closed_conductance,
        open_conductance = open_conductance,
        request_retriangularization = request_retriangularization,
    )
    topology_mutated =
        status_result.switch_status_state_mutated ||
        (order_result !== nothing && order_result.switch_order_state_mutated)
    admittance_mutated = admittance_result.switch_admittance_state_mutated
    return SwitchTopologyAdmittanceStepResult(
        status_result,
        order_result,
        admittance_result,
        topology_mutated,
        admittance_mutated,
        topology_mutated || admittance_mutated,
        topology_mutated,
        admittance_result.admittance_mutated,
        admittance_result.should_retriangularize,
    )
end

function over16_switch_retriangularization_update(
    admittance::AbstractMatrix{<:Real};
    pivot_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 0 || throw(ArgumentError("admittance must not be empty"))
    tolerance = Float64(pivot_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))

    factor = Float64.(admittance)
    for value in factor
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end
    pivot_values = zeros(Float64, node_count)

    for k in 1:node_count
        pivot = factor[k, k]
        pivot_values[k] = pivot
        abs(pivot) > tolerance ||
            throw(ArgumentError("switch retriangularization pivot is zero or below tolerance"))
        if k < node_count
            for i in (k + 1):node_count
                multiplier = factor[i, k] / pivot
                factor[i, k] = multiplier
                for j in (k + 1):node_count
                    factor[i, j] -= multiplier * factor[k, j]
                end
            end
        end
    end

    return (
        source = :over16_switch_retriangularization_update,
        outcome = :integration_boundary,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:LAST14, :SUBTS1),
        fortran_labels = (929, 930, 934, 935, 936, 940, 946, 963, 996, 1069, 1071, 1072, 1076, 1077, 1082, 1096),
        node_count = node_count,
        factor = factor,
        pivot_values = pivot_values,
        pivot_tolerance = tolerance,
        should_retriangularize = true,
        dense_factor_mutated = false,
        sparse_factor_mutated = false,
        fortran_sparse_factor_mutated = false,
        retriangularized = true,
        fortran_sparse_retriangularized = false,
        deferred_calls = [:last14, :fortran_sparse_factor_workspace, :nonlinear_inverse_columns],
        topology_mutated = false,
        admittance_mutated = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_retriangularization_update!(
    state::OVER16SwitchRetriangularizationState,
    admittance::AbstractMatrix{<:Real};
    pivot_tolerance::Real=0.0,
)
    factor_before = copy(state.factor)
    pivots_before = copy(state.pivot_values)
    factorization_count_before = state.factorization_count
    preview = over16_switch_retriangularization_update(
        admittance;
        pivot_tolerance = pivot_tolerance,
    )

    if size(state.factor) == size(preview.factor)
        state.factor .= preview.factor
    else
        state.factor = copy(preview.factor)
    end
    resize!(state.pivot_values, length(preview.pivot_values))
    state.pivot_values .= preview.pivot_values
    state.factorization_count += 1

    dense_factor_mutated = state.factor != factor_before
    pivot_values_mutated = state.pivot_values != pivots_before
    factorization_count_mutated =
        state.factorization_count != factorization_count_before
    switch_retriangularization_state_mutated =
        dense_factor_mutated || pivot_values_mutated || factorization_count_mutated
    return merge(
        preview,
        (
            factor = copy(state.factor),
            pivot_values = copy(state.pivot_values),
            factorization_count = state.factorization_count,
            dense_factor_mutated = dense_factor_mutated,
            pivot_values_mutated = pivot_values_mutated,
            factorization_count_mutated = factorization_count_mutated,
            switch_retriangularization_state_mutated =
                switch_retriangularization_state_mutated,
            sparse_factor_mutated = false,
            fortran_sparse_factor_mutated = false,
            retriangularized = true,
            fortran_sparse_retriangularized = false,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_retriangularization_update!(
    state::OVER16SwitchRetriangularizationState,
    admittance_state::OVER16SwitchAdmittanceState;
    kwargs...,
)
    return over16_switch_retriangularization_update!(
        state,
        admittance_state.admittance;
        kwargs...,
    )
end

function over16_switch_retriangularization_solve(
    factor::AbstractMatrix{<:Real},
    rhs::AbstractVector{<:Real};
    pivot_tolerance::Real=0.0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    node_count = size(factor, 1)
    length(rhs) == node_count ||
        throw(ArgumentError("rhs length must match factor size"))
    tolerance = Float64(pivot_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))

    dense_factor = Float64.(factor)
    for value in dense_factor
        isfinite(value) ||
            throw(ArgumentError("factor entries must be finite"))
    end
    solution = Float64.(rhs)

    for i in 1:node_count
        for j in 1:(i - 1)
            solution[i] -= dense_factor[i, j] * solution[j]
        end
    end
    for i in node_count:-1:1
        total = solution[i]
        for j in (i + 1):node_count
            total -= dense_factor[i, j] * solution[j]
        end
        pivot = dense_factor[i, i]
        abs(pivot) > tolerance ||
            throw(ArgumentError("switch retriangularization factor pivot is zero or below tolerance"))
        solution[i] = total / pivot
    end
    return solution
end

function over16_switch_retriangularization_solve(
    state::OVER16SwitchRetriangularizationState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_switch_retriangularization_solve(state.factor, rhs; kwargs...)
end

function over16_switch_sparse_factor_matrix(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(kk)
    node_count > 0 || throw(ArgumentError("kk must not be empty"))
    dense = zeros(Float64, node_count, node_count)
    previous_end = 0
    for row in 1:node_count
        row_end = kk[row]
        previous_end <= row_end <= length(km) ||
            throw(ArgumentError("kk entries must be nondecreasing row ends within km"))
        row_start = previous_end + 1
        row_start <= row_end ||
            throw(ArgumentError("each sparse factor row must include a diagonal marker"))
        diagonal_seen = false
        for index in row_start:row_end
            marker = km[index]
            marker != 0 ||
                throw(ArgumentError("km entries must not be zero"))
            column = marker < 0 ? -marker : marker
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            value = Float64(ykm[index])
            isfinite(value) ||
                throw(ArgumentError("ykm entries must be finite"))
            if marker < 0
                column == row ||
                    throw(ArgumentError("negative km entries must mark the row diagonal"))
                !diagonal_seen ||
                    throw(ArgumentError("sparse factor row has duplicate diagonal markers"))
                diagonal_seen = true
            end
            dense[row, column] = value
        end
        diagonal_seen ||
            throw(ArgumentError("sparse factor row is missing a diagonal marker"))
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("kk must consume every km/ykm entry"))
    return dense
end

function over16_switch_sparse_factor_matrix(
    state::OVER16SwitchSparseFactorWorkspaceState,
)
    return over16_switch_sparse_factor_matrix(state.km, state.ykm, state.kk)
end

function over16_switch_sparse_factor_update(
    factor::AbstractMatrix{<:Real};
    zero_tolerance::Real=0.0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    node_count = size(factor, 1)
    node_count > 0 || throw(ArgumentError("factor must not be empty"))
    tolerance = Float64(zero_tolerance)
    isfinite(tolerance) && tolerance >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    for row in 1:node_count
        diagonal = Float64(factor[row, row])
        isfinite(diagonal) ||
            throw(ArgumentError("factor entries must be finite"))
        push!(km, -row)
        push!(ykm, diagonal)
        for column in 1:node_count
            column == row && continue
            value = Float64(factor[row, column])
            isfinite(value) ||
                throw(ArgumentError("factor entries must be finite"))
            if abs(value) > tolerance
                push!(km, column)
                push!(ykm, value)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :over16_switch_sparse_factor_update,
        outcome = :state_mutation,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:LAST14, :SUBTS1),
        fortran_labels = (929, 930, 934, 935, 936, 940, 946, 963, 996, 1069, 1071, 1072, 1076, 1077, 1082, 1096),
        node_count = node_count,
        km = km,
        ykm = ykm,
        kk = kk,
        zero_tolerance = tolerance,
        sparse_factor_workspace_built = true,
        sparse_factor_workspace_mutated = false,
        switch_sparse_factor_workspace_state_mutated = false,
        sparse_factor_mutated = false,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:fortran_sparse_factor_ordering, :nonlinear_inverse_columns],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_sparse_factor_update(
    state::OVER16SwitchRetriangularizationState;
    kwargs...,
)
    return over16_switch_sparse_factor_update(state.factor; kwargs...)
end

function over16_switch_sparse_factor_update!(
    state::OVER16SwitchSparseFactorWorkspaceState,
    factor::AbstractMatrix{<:Real};
    zero_tolerance::Real=0.0,
)
    km_before = copy(state.km)
    ykm_before = copy(state.ykm)
    kk_before = copy(state.kk)
    count_before = state.workspace_update_count
    preview = over16_switch_sparse_factor_update(
        factor;
        zero_tolerance = zero_tolerance,
    )

    resize!(state.km, length(preview.km))
    state.km .= preview.km
    resize!(state.ykm, length(preview.ykm))
    state.ykm .= preview.ykm
    resize!(state.kk, length(preview.kk))
    state.kk .= preview.kk
    state.workspace_update_count += 1

    km_mutated = state.km != km_before
    ykm_mutated = state.ykm != ykm_before
    kk_mutated = state.kk != kk_before
    workspace_update_count_mutated =
        state.workspace_update_count != count_before
    sparse_factor_workspace_mutated =
        km_mutated || ykm_mutated || kk_mutated || workspace_update_count_mutated
    return merge(
        preview,
        (
            km = copy(state.km),
            ykm = copy(state.ykm),
            kk = copy(state.kk),
            workspace_update_count = state.workspace_update_count,
            km_mutated = km_mutated,
            ykm_mutated = ykm_mutated,
            kk_mutated = kk_mutated,
            workspace_update_count_mutated = workspace_update_count_mutated,
            sparse_factor_workspace_mutated = sparse_factor_workspace_mutated,
            switch_sparse_factor_workspace_state_mutated =
                sparse_factor_workspace_mutated,
            sparse_factor_mutated = sparse_factor_workspace_mutated,
            fortran_sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_sparse_factor_update!(
    state::OVER16SwitchSparseFactorWorkspaceState,
    retriangularization_state::OVER16SwitchRetriangularizationState;
    kwargs...,
)
    return over16_switch_sparse_factor_update!(
        state,
        retriangularization_state.factor;
        kwargs...,
    )
end

function over16_switch_sparse_factor_solve(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    factor = over16_switch_sparse_factor_matrix(km, ykm, kk)
    return over16_switch_retriangularization_solve(factor, rhs; kwargs...)
end

function over16_switch_sparse_factor_solve(
    state::OVER16SwitchSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_switch_sparse_factor_solve(state.km, state.ykm, state.kk, rhs; kwargs...)
end

function _over16_sparse_row_balance(
    node::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
)
    node_count = length(rhs)
    length(kks) >= node_count ||
        throw(ArgumentError("kks length must cover rhs nodes"))
    length(voltages) >= node_count ||
        throw(ArgumentError("voltages length must cover rhs nodes"))
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))

    pointer = kks[node]
    2 <= pointer <= length(km) + 1 ||
        throw(ArgumentError("kks row pointer must reference one past a stored row"))

    balance_current = -Float64(rhs[node])
    while pointer > 1
        pointer -= 1
        row_node = km[pointer]
        voltage_node = abs(row_node)
        1 <= voltage_node <= length(voltages) ||
            throw(ArgumentError("km row node is outside the voltage vector"))
        balance_current += Float64(ykm[pointer]) * Float64(voltages[voltage_node])
        if row_node < 0
            return balance_current
        end
    end
    throw(ArgumentError("km row is missing the negative row-start marker"))
end

function over16_switch_current_reconstruction(
    from_node::Int,
    to_node::Int,
    nextsw_entry::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
    previous_switch_current::Real=0.0,
)
    node_count = length(rhs)
    1 <= from_node <= node_count ||
        throw(ArgumentError("from_node must be within rhs"))
    1 <= to_node <= node_count ||
        throw(ArgumentError("to_node must be within rhs"))
    from_node != to_node ||
        throw(ArgumentError("switch endpoints must be distinct"))
    nextsw_entry != 0 ||
        throw(ArgumentError("nextsw_entry must be nonzero for a closed switch"))

    kcl_node = nextsw_entry < 0 ? to_node : from_node
    opposite_node = kcl_node == from_node ? to_node : from_node
    row_balance_current =
        _over16_sparse_row_balance(kcl_node, rhs, kks, km, ykm, voltages)

    updated_rhs = Float64.(rhs)
    updated_rhs[opposite_node] -= row_balance_current
    switch_current = nextsw_entry > 0 ? -row_balance_current : row_balance_current
    current_product = switch_current * Float64(previous_switch_current)
    if current_product == 0.0 && previous_switch_current != 0.0
        current_product = -1.0
    end

    return (
        kcl_node = kcl_node,
        opposite_node = opposite_node,
        selected_endpoint_index = kcl_node == from_node ? 1 : 2,
        row_balance_current = row_balance_current,
        switch_current = switch_current,
        previous_current_product = current_product,
        updated_rhs = updated_rhs,
        topology_mutated = false,
        admittance_mutated = false,
    )
end

function over16_switch_current_reconstruction_table(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    first_group_head::Int,
    rhs::AbstractVector{<:Real},
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real};
    previous_switch_currents::AbstractVector{<:Real}=Float64[],
)
    switch_count = length(nextsw)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_optional_length("previous_switch_currents", previous_switch_currents, switch_count)
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))

    updated_rhs = Float64.(rhs)
    switch_currents = zeros(Float64, switch_count)
    row_balance_currents = zeros(Float64, switch_count)
    current_products = zeros(Float64, switch_count)
    kcl_endpoint_indices = zeros(Int, switch_count)
    kcl_nodes = zeros(Int, switch_count)
    opposite_nodes = zeros(Int, switch_count)
    ordered_rows = Int[]
    visited = falses(switch_count)

    if first_group_head == 0
        any(!=(0), nextsw) &&
            throw(ArgumentError("first_group_head is zero but nextsw contains closed rows"))
        return (
            switch_currents = switch_currents,
            row_balance_currents = row_balance_currents,
            current_products = current_products,
            updated_rhs = updated_rhs,
            ordered_rows = ordered_rows,
            kcl_endpoint_indices = kcl_endpoint_indices,
            kcl_nodes = kcl_nodes,
            opposite_nodes = opposite_nodes,
            topology_mutated = false,
            admittance_mutated = false,
        )
    end
    nextsw[first_group_head] != 0 ||
        throw(ArgumentError("first_group_head must reference a closed nextsw row"))

    row = first_group_head
    for _ in 1:switch_count
        nextsw[row] != 0 ||
            throw(ArgumentError("NEXTSW ring reached an open row"))
        !visited[row] ||
            throw(ArgumentError("NEXTSW ring repeats before returning to first_group_head"))
        previous_current =
            isempty(previous_switch_currents) ? 0.0 : previous_switch_currents[row]
        reconstruction = over16_switch_current_reconstruction(
            from_nodes[row],
            to_nodes[row],
            nextsw[row],
            updated_rhs,
            kks,
            km,
            ykm,
            voltages,
            previous_current,
        )

        visited[row] = true
        push!(ordered_rows, row)
        updated_rhs = reconstruction.updated_rhs
        switch_currents[row] = reconstruction.switch_current
        row_balance_currents[row] = reconstruction.row_balance_current
        current_products[row] = reconstruction.previous_current_product
        kcl_endpoint_indices[row] = reconstruction.selected_endpoint_index
        kcl_nodes[row] = reconstruction.kcl_node
        opposite_nodes[row] = reconstruction.opposite_node

        next_row = abs(nextsw[row])
        1 <= next_row <= switch_count ||
            throw(ArgumentError("NEXTSW successor row is out of range"))
        if next_row == first_group_head
            for i in 1:switch_count
                if nextsw[i] != 0 && !visited[i]
                    throw(ArgumentError("NEXTSW contains a closed row outside the circular order"))
                end
            end
            return (
                switch_currents = switch_currents,
                row_balance_currents = row_balance_currents,
                current_products = current_products,
                updated_rhs = updated_rhs,
                ordered_rows = ordered_rows,
                kcl_endpoint_indices = kcl_endpoint_indices,
                kcl_nodes = kcl_nodes,
                opposite_nodes = opposite_nodes,
                topology_mutated = false,
                admittance_mutated = false,
            )
        end
        row = next_row
    end
    throw(ArgumentError("NEXTSW ring did not return to first_group_head"))
end

function over16_switch_current_reconstruction_table!(
    state::OVER16SwitchCurrentState,
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    first_group_head::Int,
    kks::AbstractVector{Int},
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    voltages::AbstractVector{<:Real},
)
    rhs_before = copy(state.rhs)
    currents_before = copy(state.switch_currents)
    products_before = copy(state.current_products)
    preview = over16_switch_current_reconstruction_table(
        from_nodes,
        to_nodes,
        nextsw,
        first_group_head,
        state.rhs,
        kks,
        km,
        ykm,
        voltages;
        previous_switch_currents = state.switch_currents,
    )

    resize!(state.rhs, length(preview.updated_rhs))
    state.rhs .= preview.updated_rhs
    resize!(state.switch_currents, length(preview.switch_currents))
    state.switch_currents .= preview.switch_currents
    resize!(state.current_products, length(preview.current_products))
    state.current_products .= preview.current_products

    rhs_mutated = state.rhs != rhs_before
    switch_currents_mutated = state.switch_currents != currents_before
    current_products_mutated = state.current_products != products_before
    switch_current_state_mutated =
        rhs_mutated || switch_currents_mutated || current_products_mutated
    return merge(
        preview,
        (
            updated_rhs = copy(state.rhs),
            switch_currents = copy(state.switch_currents),
            current_products = copy(state.current_products),
            rhs_mutated = rhs_mutated,
            switch_currents_mutated = switch_currents_mutated,
            current_products_mutated = current_products_mutated,
            switch_current_state_mutated = switch_current_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_post_current_transition(
    position::Int,
    switch_current::Real,
    current_product::Real,
    energy::Real,
    source_voltage_difference::Real,
    source_index::Int,
    open_time::Real,
    critical_current::Real,
    delay_time::Real,
    t::Real,
    dt::Real,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    state = abs(position)
    current = Float64(switch_current)
    product_signal = state == 3 ? -current : Float64(current_product)
    vsl_value = position <= 0 ? current : 0.0
    updated_energy = Float64(energy)
    source_power_increment = 0.0
    blocked_by_topen = false
    blocked_by_delay = false
    open_reason = :none

    if !(1 <= state <= 3)
        return (
            position = position,
            switch_current = current,
            energy = updated_energy,
            opened = false,
            open_reason = :skipped_state,
            current_product_signal = product_signal,
            source_power_increment = source_power_increment,
            vsl_value = vsl_value,
            tstop_reset = false,
            energy_report_value = 0.0,
            blocked_by_topen = false,
            blocked_by_delay = false,
            modswt_sign = 0,
            altered = false,
        )
    end

    if source_index > 0
        source_power_increment = Float64(source_voltage_difference) * current
        updated_energy += source_power_increment
        if source_power_increment < 0.0
            open_reason = :source_energy_reversal
        end
    end

    if open_reason == :none
        if state > 1 && Float64(t) < Float64(open_time)
            blocked_by_topen = true
        else
            critical_requested = abs(current) < Float64(critical_current)
            if critical_requested
                product_signal = -1.0
            end
            if Float64(t) < Float64(delay_time)
                product_signal = 1.0
                blocked_by_delay = true
            end

            if product_signal < 0.0
                if critical_requested
                    open_reason = :critical_current
                elseif state == 3
                    open_reason = :state3_current
                else
                    open_reason = :current_reversal
                end
            end
        end
    end

    if open_reason == :none
        return (
            position = position,
            switch_current = current,
            energy = updated_energy,
            opened = false,
            open_reason = blocked_by_topen ? :blocked_by_topen :
                          blocked_by_delay ? :blocked_by_delay : :none,
            current_product_signal = product_signal,
            source_power_increment = source_power_increment,
            vsl_value = vsl_value,
            tstop_reset = false,
            energy_report_value = 0.0,
            blocked_by_topen = blocked_by_topen,
            blocked_by_delay = blocked_by_delay,
            modswt_sign = 0,
            altered = false,
        )
    end

    new_state = state + 1
    if new_state == 2
        new_state = 10
    elseif new_state == 3
        new_state = 5
    end
    updated_position = position < 0 ? -new_state : new_state
    return (
        position = updated_position,
        switch_current = 0.0,
        energy = 0.0,
        opened = true,
        open_reason = open_reason,
        current_product_signal = product_signal,
        source_power_increment = source_power_increment,
        vsl_value = vsl_value,
        tstop_reset = source_index != 0 && state == 1,
        energy_report_value = state == 1 ? updated_energy * Float64(dt) : 0.0,
        blocked_by_topen = blocked_by_topen,
        blocked_by_delay = blocked_by_delay,
        modswt_sign = -1,
        altered = true,
    )
end

function over16_switch_post_current_transition_table(
    positions::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    current_products::AbstractVector{<:Real},
    energies::AbstractVector{<:Real},
    source_voltage_differences::AbstractVector{<:Real},
    source_indices::AbstractVector{Int},
    open_times::AbstractVector{<:Real},
    critical_currents::AbstractVector{<:Real},
    delay_times::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    scan_rows::AbstractVector{Int}=Int[],
)
    switch_count = length(positions)
    _over16_check_length("nextsw", nextsw, switch_count)
    _over16_check_length("switch_currents", switch_currents, switch_count)
    _over16_check_length("current_products", current_products, switch_count)
    _over16_check_length("energies", energies, switch_count)
    _over16_check_length(
        "source_voltage_differences",
        source_voltage_differences,
        switch_count,
    )
    _over16_check_length("source_indices", source_indices, switch_count)
    _over16_check_length("open_times", open_times, switch_count)
    _over16_check_length("critical_currents", critical_currents, switch_count)
    _over16_check_length("delay_times", delay_times, switch_count)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    ordered_rows = isempty(scan_rows) ? findall(!=(0), nextsw) : collect(scan_rows)
    seen_rows = falses(switch_count)
    for row in ordered_rows
        1 <= row <= switch_count ||
            throw(ArgumentError("scan_rows entries must be within switch table"))
        !seen_rows[row] ||
            throw(ArgumentError("scan_rows entries must not repeat"))
        seen_rows[row] = true
        nextsw[row] != 0 ||
            throw(ArgumentError("scan_rows entries must reference closed NEXTSW rows"))
    end

    updated_positions = collect(positions)
    updated_switch_currents = Float64.(switch_currents)
    updated_energies = Float64.(energies)
    source_power_increments = zeros(Float64, switch_count)
    vsl_values = zeros(Float64, switch_count)
    energy_report_values = zeros(Float64, switch_count)
    open_reasons = fill(:none, switch_count)
    blocked_by_topen = falses(switch_count)
    blocked_by_delay = falses(switch_count)
    opened_rows = Int[]
    tstop_reset_rows = Int[]
    modswt = Int[]

    for row in ordered_rows
        transition = over16_switch_post_current_transition(
            updated_positions[row],
            updated_switch_currents[row],
            current_products[row],
            updated_energies[row],
            source_voltage_differences[row],
            source_indices[row],
            open_times[row],
            critical_currents[row],
            delay_times[row],
            t,
            dt,
        )
        updated_positions[row] = transition.position
        updated_switch_currents[row] = transition.switch_current
        updated_energies[row] = transition.energy
        source_power_increments[row] = transition.source_power_increment
        vsl_values[row] = transition.vsl_value
        energy_report_values[row] = transition.energy_report_value
        open_reasons[row] = transition.open_reason
        blocked_by_topen[row] = transition.blocked_by_topen
        blocked_by_delay[row] = transition.blocked_by_delay
        if transition.opened
            push!(opened_rows, row)
            push!(modswt, -row)
        end
        if transition.tstop_reset
            push!(tstop_reset_rows, row)
        end
    end

    return (
        positions = updated_positions,
        switch_currents = updated_switch_currents,
        energies = updated_energies,
        source_power_increments = source_power_increments,
        vsl_values = vsl_values,
        energy_report_values = energy_report_values,
        open_reasons = open_reasons,
        blocked_by_topen = blocked_by_topen,
        blocked_by_delay = blocked_by_delay,
        opened_rows = opened_rows,
        tstop_reset_rows = tstop_reset_rows,
        modswt = modswt,
        scan_rows = ordered_rows,
        altered = !isempty(modswt),
        topology_mutated = false,
        admittance_mutated = false,
        output_mutated = false,
        tacs_executed = false,
    )
end

function over16_switch_post_current_transition_table!(
    state::OVER16SwitchPostCurrentState,
    nextsw::AbstractVector{Int},
    current_products::AbstractVector{<:Real},
    source_voltage_differences::AbstractVector{<:Real},
    source_indices::AbstractVector{Int},
    open_times::AbstractVector{<:Real},
    critical_currents::AbstractVector{<:Real},
    delay_times::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    scan_rows::AbstractVector{Int}=Int[],
)
    positions_before = copy(state.positions)
    currents_before = copy(state.switch_currents)
    energies_before = copy(state.energies)
    modswt_before = copy(state.modswt)
    preview = over16_switch_post_current_transition_table(
        state.positions,
        nextsw,
        state.switch_currents,
        current_products,
        state.energies,
        source_voltage_differences,
        source_indices,
        open_times,
        critical_currents,
        delay_times,
        t,
        dt;
        scan_rows = scan_rows,
    )

    state.positions .= preview.positions
    state.switch_currents .= preview.switch_currents
    state.energies .= preview.energies
    empty!(state.modswt)
    append!(state.modswt, preview.modswt)

    positions_mutated = state.positions != positions_before
    switch_currents_mutated = state.switch_currents != currents_before
    energies_mutated = state.energies != energies_before
    modswt_mutated = state.modswt != modswt_before
    switch_post_current_state_mutated =
        positions_mutated || switch_currents_mutated ||
        energies_mutated || modswt_mutated
    return merge(
        preview,
        (
            positions = copy(state.positions),
            switch_currents = copy(state.switch_currents),
            energies = copy(state.energies),
            modswt = copy(state.modswt),
            positions_mutated = positions_mutated,
            switch_currents_mutated = switch_currents_mutated,
            energies_mutated = energies_mutated,
            modswt_mutated = modswt_mutated,
            switch_post_current_state_mutated = switch_post_current_state_mutated,
            switch_graph_state_mutated = positions_mutated || modswt_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            output_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_alteration_rebuild_intent(
    ialter::Int,
    ktrlsw1::Int,
    ktrlsw2::Int,
    ktrlsw3::Int,
    ktrlsw4::Int=0,
    ktrlsw5::Int=0,
    ktrlsw6::Int=1;
    m4plot::Bool=false,
    yserlc_altered::Bool=false,
    kanal::Int=0,
)
    ialter in (0, 1) ||
        throw(ArgumentError("ialter must be 0 or 1"))
    ktrlsw1 >= 0 ||
        throw(ArgumentError("ktrlsw1 operation count must be nonnegative"))
    ktrlsw2 >= 0 ||
        throw(ArgumentError("ktrlsw2 closed-switch count must be nonnegative"))
    ktrlsw3 >= 0 ||
        throw(ArgumentError("ktrlsw3 triangularization count must be nonnegative"))
    ktrlsw4 >= 0 ||
        throw(ArgumentError("ktrlsw4 first-group head must be nonnegative"))
    ktrlsw5 >= 0 ||
        throw(ArgumentError("ktrlsw5 total operation count must be nonnegative"))
    ktrlsw6 in (0, 1) ||
        throw(ArgumentError("ktrlsw6 switch-logic mode must be 0 or 1"))
    (!yserlc_altered || m4plot) ||
        throw(ArgumentError("yserlc_altered requires m4plot=true"))

    effective_ialter = (ialter == 1 || (m4plot && yserlc_altered)) ? 1 : 0
    should_call_switch = effective_ialter != 0 && ktrlsw1 > 0
    switch_path = should_call_switch ? (ktrlsw6 == 0 ? :sophisticated : :simple) : :none
    updated_triangularization_count =
        effective_ialter != 0 ? ktrlsw3 + 1 : ktrlsw3
    updated_total_operation_count =
        should_call_switch ? ktrlsw5 + ktrlsw1 : ktrlsw5

    deferred_calls = Symbol[]
    if should_call_switch
        push!(deferred_calls, :switch)
    end
    if effective_ialter != 0 && kanal == 2
        push!(deferred_calls, :last14)
    end
    if effective_ialter != 0
        push!(deferred_calls, :retriangularization)
    end

    return (
        should_call_yserlc = m4plot,
        series_rlc_parameter_mutated = m4plot && yserlc_altered,
        ialter = ialter,
        effective_ialter = effective_ialter,
        ialter_after_last14 = effective_ialter != 0 && kanal == 2 ? 1 : effective_ialter,
        altered = effective_ialter != 0,
        operation_count = ktrlsw1,
        closed_switch_count = ktrlsw2,
        triangularization_count = updated_triangularization_count,
        first_group_head = ktrlsw4,
        total_operation_count = updated_total_operation_count,
        ktrlsw6 = ktrlsw6,
        should_call_switch = should_call_switch,
        switch_path = switch_path,
        lastsw_required = should_call_switch && ktrlsw6 == 0,
        should_call_last14 = effective_ialter != 0 && kanal == 2,
        should_retriangularize = effective_ialter != 0,
        should_clear_factor_workspace = effective_ialter != 0,
        should_continue_to_subts1_exit = true,
        normal_next_nchain = 17,
        kill_next_nchain = 51,
        deferred_calls = deferred_calls,
        topology_mutated = false,
        admittance_mutated = m4plot && yserlc_altered,
        sparse_factor_mutated = false,
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_switch_alteration_rebuild_update!(
    state::OVER16SwitchAlterationState;
    m4plot::Bool=false,
    yserlc_altered::Bool=false,
    kanal::Int=0,
)
    ialter_before = state.ialter
    triangularization_count_before = state.triangularization_count
    total_operation_count_before = state.total_operation_count
    preview = over16_switch_alteration_rebuild_intent(
        state.ialter,
        state.operation_count,
        state.closed_switch_count,
        state.triangularization_count,
        state.first_group_head,
        state.total_operation_count,
        state.ktrlsw6;
        m4plot = m4plot,
        yserlc_altered = yserlc_altered,
        kanal = kanal,
    )

    state.ialter = preview.ialter_after_last14
    state.triangularization_count = preview.triangularization_count
    state.total_operation_count = preview.total_operation_count

    ialter_mutated = state.ialter != ialter_before
    triangularization_count_mutated =
        state.triangularization_count != triangularization_count_before
    total_operation_count_mutated =
        state.total_operation_count != total_operation_count_before
    switch_alteration_state_mutated =
        ialter_mutated || triangularization_count_mutated ||
        total_operation_count_mutated
    return merge(
        preview,
        (
            ialter = state.ialter,
            triangularization_count = state.triangularization_count,
            total_operation_count = state.total_operation_count,
            ialter_mutated = ialter_mutated,
            triangularization_count_mutated = triangularization_count_mutated,
            total_operation_count_mutated = total_operation_count_mutated,
            switch_alteration_state_mutated = switch_alteration_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            sparse_factor_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function over16_switch_bvalue_export(
    kpos::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    base_count::Int=0,
)
    switch_count = length(kpos)
    _over16_check_length("nextsw", nextsw, switch_count)
    _over16_check_length("switch_currents", switch_currents, switch_count)
    base_count >= 0 || throw(ArgumentError("base_count must be nonnegative"))

    export_rows = Int[]
    bvalue_indices = Int[]
    values = Float64[]
    count = base_count
    for row in 1:switch_count
        if kpos[row] >= 0
            continue
        end
        count += 1
        push!(export_rows, row)
        push!(bvalue_indices, count)
        push!(values, nextsw[row] == 0 ? 0.0 : Float64(switch_currents[row]))
    end

    return (
        export_rows = export_rows,
        bvalue_indices = bvalue_indices,
        values = values,
        final_count = count,
        output_mutated = false,
        tacs_executed = false,
    )
end

function over16_switch_bvalue_export!(
    state::OVER16SwitchBValueExportState,
    kpos::AbstractVector{Int},
    nextsw::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
)
    bvalue_before = copy(state.bvalue)
    output_count_before = state.output_count
    preview = over16_switch_bvalue_export(
        kpos,
        nextsw,
        switch_currents,
        state.output_count,
    )

    if length(state.bvalue) < preview.final_count
        resize!(state.bvalue, preview.final_count)
    end
    for (index, value) in zip(preview.bvalue_indices, preview.values)
        state.bvalue[index] = value
    end
    state.output_count = preview.final_count

    bvalue_mutated = state.bvalue != bvalue_before
    output_count_mutated = state.output_count != output_count_before
    switch_bvalue_state_mutated = bvalue_mutated || output_count_mutated
    return merge(
        preview,
        (
            bvalue = copy(state.bvalue),
            output_count = state.output_count,
            bvalue_mutated = bvalue_mutated,
            output_count_mutated = output_count_mutated,
            switch_bvalue_state_mutated = switch_bvalue_state_mutated,
            output_mutated = switch_bvalue_state_mutated,
            topology_mutated = false,
            admittance_mutated = false,
            tacs_executed = false,
            solvum_executed = false,
        ),
    )
end

function stamp!(y::AbstractMatrix{Float64}, _rhs::AbstractVector{Float64},
                s::Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch}, t::Float64, _dt::Float64)
    stamp_conductance!(y, s.a, s.b, switch_conductance(s, t))
    return nothing
end

backward_euler_companion_supported(
    ::Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch},
) = true

update!(::Union{IdealSwitch,TimeSwitch}, _voltages::AbstractVector{Float64},
        _dt::Float64) = nothing

update!(::CurrentZeroSwitch, _voltages::AbstractVector{Float64}, _dt::Float64) = nothing

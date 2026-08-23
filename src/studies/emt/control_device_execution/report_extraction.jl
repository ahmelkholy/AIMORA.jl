function final_voltage_pu(trace::DeckEMTTrace, node::Symbol)
    index = get(trace.node_map, node, 0)
    index == 0 && throw(ArgumentError("unknown node $node"))
    return trace.voltage_pu[index, end]
end

function final_output_pu(trace::DeckEMTTrace, channel::Symbol)
    index = findfirst(==(channel), trace.output_channel_names)
    index === nothing && throw(ArgumentError("unknown output channel $channel"))
    return trace.output_pu[index, end]
end

function _last_over16_accepted_update(boundary_run)
    isempty(boundary_run.over16_updates) && return nothing
    return boundary_run.over16_updates[end].over16_update
end

_nested_value_count(values::AbstractVector) = mapreduce(length, +, values; init = 0)

function _deck_over16_output_request_line_number_count(run)
    return length(run.deck_over16_output_channel_line_numbers) +
           length(run.deck_over16_branch_voltage_output_line_numbers) +
           length(run.deck_over16_branch_current_output_line_numbers) +
           length(run.deck_over16_branch_power_output_line_numbers)
end

function _last_over16_tacs_or_csup_update(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        if update.tacs_result !== nothing ||
           update.ntacs3_result !== nothing ||
           update.tacs_linear_result !== nothing ||
           update.tacs_post_solve_result !== nothing ||
           update.csup_pre_solve_result !== nothing ||
           update.csup_termination_result !== nothing
            return update
        end
    end
    return nothing
end

function _last_over16_ntacs3_result(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        result = update.ntacs3_result
        result !== nothing && return result
    end
    return nothing
end

function _last_over16_output_report_result(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        update.output_report_result === nothing || return update.output_report_result
    end
    return nothing
end

function _last_over16_branch_power_result(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return nothing
    return output_result.power_result
end

function _over16_output_report_values(boundary_run)
    for update in Iterators.reverse(_over16_boundary_pass_updates(boundary_run.over16_updates))
        update.output_report_result === nothing || return Float64.(update.output_volti_values)
    end
    return Float64[]
end

function _over16_voltage_output_values(boundary_run)
    values = _over16_output_report_values(boundary_run)
    count = boundary_run.deck_over16_output_channel_count
    count > 0 || return Float64[]
    length(values) >= count + 1 || return Float64[]
    return Float64[values[index] for index in 2:(count + 1)]
end

function _boundary_voltage_output_values(boundary_run)
    if hasproperty(boundary_run, :steady_state_initial_sample_applied) &&
       boundary_run.steady_state_initial_sample_applied
        return Float64.(boundary_run.steady_state_initial_output_voltage_values)
    end
    return _over16_voltage_output_values(boundary_run)
end

function _over16_output_energy_state_values(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return Float64[]
    return Float64.(output_result.updated_bnrg_values)
end

function _over16_branch_output_bvalue_values(boundary_run)
    output_result = _last_over16_output_report_result(boundary_run)
    output_result === nothing && return Float64[]
    return Float64.(output_result.updated_bvalue_values)
end

function _over16_branch_voltage_output_values(boundary_run)
    bvalues = _over16_branch_output_bvalue_values(boundary_run)
    count = boundary_run.deck_over16_branch_voltage_output_count
    length(bvalues) >= count || return Float64[]
    return Float64[bvalues[index] for index in 1:count]
end

function _over16_branch_current_output_values(boundary_run)
    bvalues = _over16_branch_output_bvalue_values(boundary_run)
    voltage_count = boundary_run.deck_over16_branch_voltage_output_count
    current_count = boundary_run.deck_over16_branch_current_output_count
    end_index = voltage_count + current_count
    length(bvalues) >= end_index || return Float64[]
    return Float64[bvalues[index] for index in (voltage_count + 1):end_index]
end

function _over16_branch_power_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.power_values)
end

function _over16_branch_energy_output_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.energy_output_values)
end

function _over16_branch_power_voltage_factor_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.voltage_factor_values)
end

function _over16_branch_power_current_factor_values(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Float64[]
    return Float64.(power_result.current_factor_values)
end

function _over16_branch_power_voltage_selectors(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Int[]
    return Int.(power_result.voltage_selectors)
end

function _over16_branch_power_current_selectors(boundary_run)
    power_result = _last_over16_branch_power_result(boundary_run)
    power_result === nothing && return Int[]
    return Int.(power_result.current_selectors)
end

function _over16_tacs_linear_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_linear_xtcs_values)
end

function _over16_tacs_linear_rsblk_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_rsblk_values)
end

function _over16_tacs_post_solve_parsup_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_tacs_parsup_values)
end

function _over16_tacs_utility_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_xtcs_values)
end

function _over16_ntacs3_elec_output_result(boundary_run)
    result = _last_over16_ntacs3_result(boundary_run)
    result === nothing && return nothing
    return result.elec_output_copy_result
end

function _over16_ntacs3_elec_output_values(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Float64[]
    return Float64.(result.values)
end

function _over16_ntacs3_elec_output_update_indices(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    return Int.(result.xtcs_update_indices)
end

function _over16_elec_tacsto_report_kind_symbols(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Symbol[]
    hasproperty(result, :elec_report_records) || return Symbol[]
    return Symbol[record.kind for record in result.elec_report_records]
end

function _over16_elec_tacsto_report_has_value(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_report_records) || return Int[]
    return Int[record.has_value ? 1 : 0 for record in result.elec_report_records]
end

function _over16_elec_tacsto_report_values(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Float64[]
    hasproperty(result, :elec_report_records) || return Float64[]
    return Float64[record.has_value ? record.value : 0.0 for record in result.elec_report_records]
end

function _over16_elec_tacsto_error_stpflg(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_error_records) || return Int[]
    return Int[record.stpflg for record in result.elec_error_records]
end

function _over16_elec_tacsto_error_labels(boundary_run)
    result = _over16_ntacs3_elec_output_result(boundary_run)
    result === nothing && return Int[]
    hasproperty(result, :elec_error_records) || return Int[]
    return Int[record.fortran_label for record in result.elec_error_records]
end

function _over16_csup_xtcs_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_csup_xtcs_values)
end

function _over16_csup_parsup_values(boundary_run)
    update = _last_over16_tacs_or_csup_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_csup_parsup_values)
end

function _last_over16_switch_update(boundary_run)
    for step_update in Iterators.reverse(boundary_run.over16_updates)
        update = step_update.over16_update
        if update.switch_state_mutated ||
           update.switch_scan_result !== nothing ||
           update.switch_operation_result !== nothing ||
           update.switch_status_result !== nothing ||
           update.switch_order_result !== nothing ||
           update.switch_admittance_result !== nothing ||
           update.switch_current_result !== nothing ||
           update.switch_post_current_result !== nothing ||
           update.switch_bvalue_result !== nothing ||
           update.switch_alteration_result !== nothing ||
           update.switch_retriangularization_result !== nothing ||
           update.switch_sparse_factor_result !== nothing ||
           update.switch_fortran_sparse_factor_result !== nothing ||
           update.switch_network_solution_result !== nothing
            return update
        end
    end
    return nothing
end

function _over16_switch_closed_mask_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int[value ? 1 : 0 for value in update.output_switch_closed_mask]
end

function _over16_switch_nextsw_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_nextsw)
end

function _over16_switch_kode_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_kode)
end

function _over16_switch_current_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_currents)
end

function _over16_switch_current_product_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_current_products)
end

function _over16_switch_position_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_positions)
end

function _over16_switch_operation_modswt_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_operation_modswt)
end

function _over16_switch_post_current_modswt_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_post_current_modswt)
end

function _over16_switch_bvalue_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_bvalue)
end

function _over16_switch_admittance_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(vec(update.output_switch_admittance))
end

function _over16_switch_conductance_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_conductances)
end

function _over16_switch_factor_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(vec(update.output_switch_factor))
end

function _over16_switch_factor_pivot_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_factor_pivots)
end

function _over16_switch_sparse_factor_km_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_sparse_factor_km)
end

function _over16_switch_sparse_factor_ykm_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_sparse_factor_ykm)
end

function _over16_switch_sparse_factor_kk_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_sparse_factor_kk)
end

function _over16_switch_fortran_sparse_factor_km_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_fortran_sparse_factor_km)
end

function _over16_switch_fortran_sparse_factor_ykm_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_fortran_sparse_factor_ykm)
end

function _over16_switch_fortran_sparse_factor_kk_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_switch_fortran_sparse_factor_kk)
end

function _over16_switch_rhs_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_rhs)
end

function _over16_switch_network_solution_values(boundary_run)
    update = _last_over16_switch_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_switch_network_solution)
end

function _over16_sparse_switch_state_flow_results(boundary_run)
    flows = Any[]
    for update in boundary_run.over16_updates
        flow = _over16_step_sparse_switch_state_flow_result(update)
        flow === nothing || push!(flows, flow)
    end
    return flows
end

function _over16_sparse_switch_state_flow_pass_count_values(boundary_run)
    return Int[flow.pass_count for flow in _over16_sparse_switch_state_flow_results(boundary_run)]
end

function _over16_sparse_switch_state_flow_pass_total(boundary_run)
    return sum(_over16_sparse_switch_state_flow_pass_count_values(boundary_run); init = 0)
end

function _over16_sparse_switch_state_flow_repeated_after_post_current_count(boundary_run)
    return sum(
        Int[flow.repeated_after_post_current_count for
            flow in _over16_sparse_switch_state_flow_results(boundary_run)];
        init = 0,
    )
end

function _over16_sparse_switch_state_flow_pass_lengths(boundary_run, field::Symbol)
    lengths = Int[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for values in getproperty(flow, field)
            push!(lengths, length(values))
        end
    end
    return lengths
end

function _over16_sparse_switch_state_flow_pass_int_values(boundary_run, field::Symbol)
    values = Int[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for pass_values in getproperty(flow, field)
            append!(values, Int.(pass_values))
        end
    end
    return values
end

function _over16_sparse_switch_state_flow_pass_float_values(boundary_run, field::Symbol)
    values = Float64[]
    for flow in _over16_sparse_switch_state_flow_results(boundary_run)
        for pass_values in getproperty(flow, field)
            append!(values, Float64.(pass_values))
        end
    end
    return values
end

function _is_over16_nonlinear_update(update)
    return update.nonlinear_inverse_state_mutated ||
           update.nonlinear_source_column_state_mutated ||
           update.nonlinear_current_state_mutated ||
           update.nonlinear_source_column_result !== nothing ||
           update.nonlinear_inverse_result !== nothing ||
           update.nonlinear_current_result !== nothing
end

function _over16_step_owner_updates(step_update)
    flow = _over16_step_sparse_switch_state_flow_result(step_update)
    return flow === nothing ? (step_update.over16_update,) : flow.pass_updates
end

function _last_over16_nonlinear_update(boundary_run)
    for step_update in Iterators.reverse(boundary_run.over16_updates)
        for update in Iterators.reverse(_over16_step_owner_updates(step_update))
            _is_over16_nonlinear_update(update) && return update
        end
    end
    return nothing
end

function _over16_nonlinear_state_mutation_count(boundary_run)
    count = 0
    for step_update in boundary_run.over16_updates
        for update in _over16_step_owner_updates(step_update)
            count += (
                update.nonlinear_inverse_state_mutated ||
                update.nonlinear_source_column_state_mutated ||
                update.nonlinear_current_state_mutated
            ) ? 1 : 0
        end
    end
    return count
end

function _over16_nonlinear_inverse_znonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_znonl)
end

function _over16_nonlinear_inverse_anonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_anonl)
end

function _over16_nonlinear_inverse_voltbc_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_voltbc)
end

function _over16_nonlinear_inverse_vzero_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_vzero)
end

function _over16_nonlinear_inverse_vnonl_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_inverse_vnonl)
end

function _over16_nonlinear_inverse_ilast_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Int[]
    return Int.(update.output_nonlinear_inverse_ilast)
end

function _over16_nonlinear_current_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_current_values)
end

function _over16_nonlinear_cursub_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_cursub_values)
end

function _over16_nonlinear_vchar_values(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return Float64[]
    return Float64.(update.output_nonlinear_vchar_values)
end

function _last_over16_nonlinear_current_result(boundary_run)
    update = _last_over16_nonlinear_update(boundary_run)
    update === nothing && return nothing
    update.nonlinear_current_result === nothing && return nothing
    return update.nonlinear_current_result
end

function _over16_nonlinear_current_results(boundary_run)
    results = Any[]
    for step_update in boundary_run.over16_updates
        for update in _over16_step_owner_updates(step_update)
            update.nonlinear_current_result === nothing && continue
            push!(results, update.nonlinear_current_result)
        end
    end
    return results
end

function _over16_nonlinear_current_result_counts(boundary_run, field::Symbol)
    return Int[Int(getproperty(result, field))
               for result in _over16_nonlinear_current_results(boundary_run)]
end

function _over16_nonlinear_current_result_count(boundary_run, field::Symbol)
    return sum(_over16_nonlinear_current_result_counts(boundary_run, field); init = 0)
end

function _over16_nonlinear_current_result_lengths(boundary_run, field::Symbol)
    return Int[length(getproperty(result, field))
               for result in _over16_nonlinear_current_results(boundary_run)]
end

function _over16_nonlinear_current_result_float_values(boundary_run, field::Symbol)
    values = Float64[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(values, Float64.(getproperty(result, field)))
    end
    return values
end

function _over16_nonlinear_current_mutation_order(boundary_run)
    mutation_order = Symbol[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(mutation_order, result.mutation_order)
    end
    return mutation_order
end

function _over16_nonlinear_current_deferred_calls(boundary_run)
    deferred_calls = Symbol[]
    for result in _over16_nonlinear_current_results(boundary_run)
        append!(deferred_calls, result.deferred_calls)
    end
    return deferred_calls
end

function _deck_time_switch_closed_at(initially_closed::Bool,
                                     close_time_s::Float64,
                                     open_time_s::Float64,
                                     time_s::Float64)
    closed = initially_closed
    if isfinite(close_time_s) && time_s >= close_time_s
        closed = true
    end
    if isfinite(open_time_s) && time_s >= open_time_s
        closed = false
    end
    return closed
end

function _deck_time_switch_final_closed_flags(boundary_run)
    final_time_s = boundary_run.trace.time_s[end]
    return Int[
        _deck_time_switch_closed_at(
            boundary_run.deck_time_switch_initially_closed_flags[index],
            boundary_run.deck_time_switch_close_time_s_values[index],
            boundary_run.deck_time_switch_open_time_s_values[index],
            final_time_s,
        ) ? 1 : 0
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_final_conductance_values(boundary_run,
                                                    final_closed_flags::Vector{Int})
    return Float64[
        final_closed_flags[index] == 1 ?
        boundary_run.deck_time_switch_on_conductance_values[index] :
        boundary_run.deck_time_switch_off_conductance_values[index]
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_switch_grounded_conductance_current(
    boundary_run,
    endpoint::Int,
    voltage::AbstractVector{<:Real},
)
    hasproperty(boundary_run, :context) ||
        throw(ArgumentError("deck run does not expose its typed timestep context"))
    return _deck_switch_grounded_conductance_current(
        boundary_run.context,
        endpoint,
        voltage,
    )
end

function _deck_switch_grounded_conductance_current(
    context::EMTStepContext,
    endpoint::Int,
    voltage::AbstractVector{<:Real},
)
    endpoint == 0 && return nothing
    voltage_values = _deck_float64_voltage_values(voltage)
    return _deck_switch_grounded_conductance_current(
        context.system.elements,
        endpoint,
        voltage_values,
        context.dt_s,
    )
end

_deck_switch_grounded_conductance_current(
    ::Tuple{},
    ::Int,
    ::AbstractVector{Float64},
    ::Float64,
) = nothing

function _deck_switch_grounded_conductance_current(
    elements::Tuple,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    element = first(elements)
    if applicable(branch_companion_snapshot, element, voltage, dt_s)
        snapshot = branch_companion_snapshot(element, voltage, dt_s)
        if snapshot !== nothing && snapshot.kind == :conductance &&
           (
               (snapshot.a == endpoint && snapshot.b == 0) ||
               (snapshot.b == endpoint && snapshot.a == 0)
           )
            return snapshot.conductance * _deck_node_voltage(voltage, endpoint)
        end
    end
    return _deck_switch_grounded_conductance_current(
        Base.tail(elements),
        endpoint,
        voltage,
        dt_s,
    )
end


function _deck_switch_grounded_conductance_current(
    elements::NodalElementSequence,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    return _deck_switch_grounded_conductance_current_batches(
        elements.contiguous_type_batches,
        endpoint,
        voltage,
        dt_s,
    )
end

_deck_switch_grounded_conductance_current_batches(
    ::Tuple{},
    ::Int,
    ::AbstractVector{Float64},
    ::Float64,
) = nothing

function _deck_switch_grounded_conductance_current_batches(
    batches::Tuple,
    endpoint::Int,
    voltage::AbstractVector{Float64},
    dt_s::Float64,
)
    for element in first(batches)
        if applicable(branch_companion_snapshot, element, voltage, dt_s)
            snapshot = branch_companion_snapshot(element, voltage, dt_s)
            if snapshot !== nothing && snapshot.kind == :conductance &&
               (
                   (snapshot.a == endpoint && snapshot.b == 0) ||
                   (snapshot.b == endpoint && snapshot.a == 0)
               )
                return snapshot.conductance *
                    _deck_node_voltage(voltage, endpoint)
            end
        end
    end
    return _deck_switch_grounded_conductance_current_batches(
        Base.tail(batches),
        endpoint,
        voltage,
        dt_s,
    )
end

function _deck_closed_switch_grounded_path_current(
    boundary_run,
    switch_index::Int,
    voltage::AbstractVector{<:Real},
)
    from_node = boundary_run.deck_time_switch_from_node_indices[switch_index]
    to_node = boundary_run.deck_time_switch_to_node_indices[switch_index]
    to_ground_current =
        _deck_switch_grounded_conductance_current(boundary_run, to_node, voltage)
    to_ground_current !== nothing && return to_ground_current
    from_ground_current =
        _deck_switch_grounded_conductance_current(boundary_run, from_node, voltage)
    from_ground_current !== nothing && return -from_ground_current
    return nothing
end

_deck_node_voltage(voltage::AbstractVector{<:Real}, index::Int) =
    index == 0 ? 0.0 : voltage[index]

_deck_float64_voltage_values(voltage::AbstractVector{Float64}) = voltage
_deck_float64_voltage_values(voltage::AbstractVector{<:Real}) = Float64.(voltage)

function _deck_final_voltage_values(boundary_run)
    required_node_count = max(
        maximum(boundary_run.deck_time_switch_from_node_indices; init = 0),
        maximum(boundary_run.deck_time_switch_to_node_indices; init = 0),
    )
    update = _last_over16_accepted_update(boundary_run)
    if update !== nothing && hasproperty(update, :output_e_values) &&
       length(update.output_e_values) >= required_node_count
        return Float64.(update.output_e_values)
    end
    return Float64.(boundary_run.trace.voltage_pu[:, end])
end

function _deck_time_switch_final_voltage_values(boundary_run)
    final_voltage = _deck_final_voltage_values(boundary_run)
    return Float64[
        _deck_node_voltage(final_voltage, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(final_voltage, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_report_voltage_values(boundary_run, voltage::AbstractVector{<:Real})
    voltage_values = _deck_float64_voltage_values(voltage)
    return Float64[
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_final_current_values(final_voltage_values::Vector{Float64},
                                                final_conductance_values::Vector{Float64})
    return Float64[
        final_conductance_values[index] * final_voltage_values[index]
        for index in 1:length(final_voltage_values)
    ]
end

function _deck_current_extinction_report_state(element)
    if element isa CurrentZeroSwitch
        return (
            closed = element.closed,
            current_initialized = element.current_initialized,
            previous_current = element.previous_current,
        )
    elseif element isa TimeSwitch && element.current_extinction !== nothing
        state = element.current_extinction
        return (
            closed = state.closed,
            current_initialized = state.current_initialized,
            previous_current = state.previous_current,
        )
    end
    return nothing
end

function _deck_runtime_switch_element(boundary_run, switch_index::Int)
    owner = hasproperty(boundary_run, :context) ? boundary_run.context : boundary_run
    hasproperty(owner, :system) || return nothing
    switch_index >= 1 || return nothing
    return _deck_runtime_switch_element(owner.system.elements, switch_index, 0)
end

_deck_runtime_switch_element(::Tuple{}, ::Int, ::Int) = nothing

function _deck_runtime_switch_element(
    elements::Tuple,
    switch_index::Int,
    found::Int,
)
    element = first(elements)
    next_found =
        element isa Union{TimeSwitch,CurrentZeroSwitch} ? found + 1 : found
    if next_found == switch_index && next_found != found
        return element
    end
    return _deck_runtime_switch_element(
        Base.tail(elements),
        switch_index,
        next_found,
    )
end


function _deck_runtime_switch_element(
    elements::NodalElementSequence,
    switch_index::Int,
    found::Int,
)
    return _deck_runtime_switch_element_batches(
        elements.contiguous_type_batches,
        switch_index,
        found,
    )
end

_deck_runtime_switch_element_batches(::Tuple{}, ::Int, ::Int) = nothing

function _deck_runtime_switch_element_batches(
    batches::Tuple,
    switch_index::Int,
    found::Int,
)
    batch = first(batches)
    next_found = found
    if eltype(batch) <: Union{TimeSwitch,CurrentZeroSwitch}
        next_found += length(batch)
        if found < switch_index <= next_found
            return @inbounds batch[switch_index - found]
        end
    end
    return _deck_runtime_switch_element_batches(
        Base.tail(batches),
        switch_index,
        next_found,
    )
end

function _deck_time_switch_report_current_values!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
    time_s::Real,
)
    closed_flags = context.switch_closed_step_flags
    conductance_values = context.switch_conductance_step_values
    voltage_differences = context.switch_voltage_step_values
    current_values = context.switch_current_step_values
    time = Float64(time_s)
    for index in 1:context.deck_time_switch_count
        closed = _deck_time_switch_closed_at(
            context.deck_time_switch_initially_closed_flags[index],
            context.deck_time_switch_close_time_s_values[index],
            context.deck_time_switch_open_time_s_values[index],
            time,
        )
        element = _deck_runtime_switch_element(context, index)
        extinction = _deck_current_extinction_report_state(element)
        extinction === nothing || (closed = extinction.closed)
        closed_flags[index] = closed ? 1 : 0
        conductance = closed ?
            context.deck_time_switch_on_conductance_values[index] :
            context.deck_time_switch_off_conductance_values[index]
        conductance_values[index] = conductance
        voltage_difference =
            _deck_node_voltage(voltage, context.deck_time_switch_from_node_indices[index]) -
            _deck_node_voltage(voltage, context.deck_time_switch_to_node_indices[index])
        voltage_differences[index] = voltage_difference
        current_values[index] = conductance * voltage_difference
    end
    for index in 1:context.deck_time_switch_count
        element = _deck_runtime_switch_element(context, index)
        extinction = _deck_current_extinction_report_state(element)
        if extinction !== nothing
            if !extinction.closed
                current_values[index] = 0.0
            elseif extinction.current_initialized
                current_values[index] = extinction.previous_current
            else
                grounded_current = _deck_closed_switch_grounded_path_current(
                    context,
                    index,
                    voltage,
                )
                grounded_current === nothing ||
                    (current_values[index] = grounded_current)
            end
            continue
        end
        closed_flags[index] == 1 || continue
        grounded_current = _deck_closed_switch_grounded_path_current(
            context,
            index,
            voltage,
        )
        grounded_current === nothing ||
            (current_values[index] = grounded_current)
    end
    return current_values
end

function _deck_time_switch_report_current_values(
    boundary_run,
    voltage::AbstractVector{<:Real},
    time_s::Real,
)
    voltage_values = _deck_float64_voltage_values(voltage)
    closed_flags = Int[
        _deck_time_switch_closed_at(
            boundary_run.deck_time_switch_initially_closed_flags[index],
            boundary_run.deck_time_switch_close_time_s_values[index],
            boundary_run.deck_time_switch_open_time_s_values[index],
            Float64(time_s),
        ) ? 1 : 0
        for index in 1:boundary_run.deck_time_switch_count
    ]
    for index in eachindex(closed_flags)
        element = _deck_runtime_switch_element(boundary_run, index)
        extinction = _deck_current_extinction_report_state(element)
        extinction === nothing && continue
        closed_flags[index] = extinction.closed ? 1 : 0
    end
    conductance_values =
        _deck_time_switch_final_conductance_values(boundary_run, closed_flags)
    voltage_differences = Float64[
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_from_node_indices[index]) -
        _deck_node_voltage(voltage_values, boundary_run.deck_time_switch_to_node_indices[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
    current_values = _deck_time_switch_final_current_values(
        voltage_differences,
        conductance_values,
    )
    for index in 1:boundary_run.deck_time_switch_count
        element = _deck_runtime_switch_element(boundary_run, index)
        extinction = _deck_current_extinction_report_state(element)
        if extinction !== nothing
            if !extinction.closed
                current_values[index] = 0.0
            elseif extinction.current_initialized
                current_values[index] = extinction.previous_current
            else
                grounded_current = _deck_closed_switch_grounded_path_current(
                    boundary_run,
                    index,
                    voltage_values,
                )
                grounded_current === nothing ||
                    (current_values[index] = grounded_current)
            end
            continue
        end
        closed_flags[index] == 1 || continue
        grounded_current =
            _deck_closed_switch_grounded_path_current(boundary_run, index, voltage_values)
        grounded_current === nothing && continue
        current_values[index] = grounded_current
    end
    return current_values
end

function _deck_time_switch_final_power_values(final_voltage_values::Vector{Float64},
                                              final_current_values::Vector{Float64})
    return Float64[
        final_voltage_values[index] * final_current_values[index]
        for index in 1:length(final_voltage_values)
    ]
end

function _deck_time_switch_report_power_values(
    boundary_run,
    voltage::AbstractVector{<:Real},
    current_values::AbstractVector{<:Real},
)
    return Float64[
        _deck_node_voltage(voltage, boundary_run.deck_time_switch_from_node_indices[index]) *
        Float64(current_values[index])
        for index in 1:boundary_run.deck_time_switch_count
    ]
end

function _deck_time_switch_report_power_values!(
    context::EMTStepContext,
    voltage::AbstractVector{Float64},
    current_values::AbstractVector{Float64},
)
    power_values = context.switch_power_step_values
    for index in 1:context.deck_time_switch_count
        power_values[index] =
            _deck_node_voltage(voltage, context.deck_time_switch_from_node_indices[index]) *
            current_values[index]
    end
    return power_values
end

_switch_signal_output_code(output_code::Integer) =
    output_code > 3 ? 3 : output_code

_switch_current_report_selected(output_code::Integer) =
    _switch_signal_output_code(output_code) in (1, 3)

_switch_voltage_report_selected(output_code::Integer) =
    _switch_signal_output_code(output_code) >= 2

function _deck_switch_voltage_report_indices(
    output_codes::AbstractVector{<:Integer},
    switch_count::Int,
)
    row_count = min(length(output_codes), switch_count)
    return Int[
        index for index in 1:row_count
        if _switch_voltage_report_selected(output_codes[index])
    ]
end

function _deck_switch_current_report_indices(
    output_codes::AbstractVector{<:Integer},
    switch_count::Int,
)
    row_count = min(length(output_codes), switch_count)
    return Int[
        index for index in 1:row_count
        if _switch_current_report_selected(output_codes[index])
    ]
end

function _deck_switch_current_report_indices(boundary_run)
    return _deck_switch_current_report_indices(
        boundary_run.deck_over5_switch_output_codes,
        boundary_run.deck_time_switch_count,
    )
end

function _deck_switch_voltage_report_indices(boundary_run)
    return _deck_switch_voltage_report_indices(
        boundary_run.deck_over5_switch_output_codes,
        boundary_run.deck_time_switch_count,
    )
end

_selected_symbol_values(values::AbstractVector{Symbol}, indices::AbstractVector{Int}) =
    Symbol[values[index] for index in indices]

_selected_int_values(values::AbstractVector{<:Integer}, indices::AbstractVector{Int}) =
    Int[values[index] for index in indices]

_selected_float_values(values::AbstractVector{<:Real}, indices::AbstractVector{Int}) =
    Float64[values[index] for index in indices]

function _deck_switch_current_report_names(boundary_run, indices::AbstractVector{Int})
    return _deck_switch_current_report_names(boundary_run.deck_over5_switch_names, indices)
end

function _deck_switch_voltage_report_names(boundary_run, indices::AbstractVector{Int})
    return _deck_switch_voltage_report_names(boundary_run.deck_over5_switch_names, indices)
end

function _deck_switch_voltage_report_names(
    switch_names::AbstractVector{Symbol},
    indices::AbstractVector{Int},
)
    return Symbol[
        Symbol("switch_voltage_", String(switch_names[index]))
        for index in indices
    ]
end

function _deck_switch_current_report_names(
    switch_names::AbstractVector{Symbol},
    indices::AbstractVector{Int},
)
    return Symbol[
        Symbol("switch_current_", String(switch_names[index]))
        for index in indices
    ]
end

function _deck_electrical_report_channel_values(voltage_output_values::AbstractVector{<:Real},
                                                branch_voltage_output_values::AbstractVector{<:Real},
                                                switch_voltage_report_values::AbstractVector{<:Real},
                                                branch_current_output_values::AbstractVector{<:Real},
                                                switch_current_report_values::AbstractVector{<:Real})
    return Float64[
        Float64.(voltage_output_values)...,
        Float64.(branch_voltage_output_values)...,
        Float64.(switch_voltage_report_values)...,
        Float64.(switch_current_report_values)...,
        Float64.(branch_current_output_values)...,
    ]
end

function _over16_post_extrema_results(boundary_run)
    return [
        update.post_extrema_result
        for update in _over16_boundary_pass_updates(boundary_run.over16_updates)
        if update.post_extrema_result !== nothing
    ]
end

function _over16_post_extrema_time_after_values(boundary_run)
    return Float64.([result.time_after for result in _over16_post_extrema_results(boundary_run)])
end

function _over16_post_extrema_step_after_values(boundary_run)
    return Int[result.istep_after for result in _over16_post_extrema_results(boundary_run)]
end

function _over16_post_extrema_final_time_s(boundary_run)
    values = _over16_post_extrema_time_after_values(boundary_run)
    return isempty(values) ? 0.0 : values[end]
end

import ..Branches: EMTElement,
                   backward_euler_companion_supported,
                   stamp!,
                   stamp_conductance!,
                   update!

export IdealSwitch,
       TimeSwitch,
       CurrentZeroSwitch,
       OVER16SwitchAdmittanceState,
       OVER16SwitchAlterationState,
       OVER16SwitchBValueExportState,
       OVER16SwitchCurrentState,
       OVER16SwitchOperationState,
       OVER16SwitchPostCurrentState,
       OVER16SwitchRetriangularizationState,
       OVER16SwitchScanState,
       OVER16SwitchSparseFactorWorkspaceState,
       OVER16FortranSparseFactorWorkspaceState,
       GroupedSparseNetworkWorkspace,
       OVER16SwitchTopologyState,
       switch_closed,
       switch_conductance,
       configure_current_extinction!,
       current_extinction_enabled,
       prepare_current_zero_switch!,
       apply_current_zero_transition!,
       update_current_zero_switch!,
       over16_switch_margin_history,
       over16_switch_close_critical_current,
       over16_gap_energy_check,
       over16_voltage_open_delay,
       over16_switch_tail_current_injection,
       over16_open_switch_close_decision,
       over16_closed_switch_open_decision,
       over16_switch_clamp_action,
       over16_switch_position_update,
       over16_a8sw_delayed_open_step,
       over16_a8sw_tail_current_table,
       over16_controlled_switch_scan_step,
       over16_controlled_switch_table_scan,
       over16_controlled_switch_table_scan!,
       over16_switch_operation_schedule,
       over16_switch_operation_schedule!,
       over16_switch_status_update,
       over16_switch_status_update!,
       over16_switch_simple_ordering,
       over16_switch_simple_ordering!,
       over16_switch_admittance_update,
       over16_switch_admittance_update!,
       over16_switch_topology_admittance_update!,
       over16_switch_retriangularization_update,
       over16_switch_retriangularization_update!,
       over16_switch_retriangularization_solve,
       over16_switch_sparse_factor_matrix,
       over16_switch_sparse_factor_update,
       over16_switch_sparse_factor_update!,
       over16_switch_sparse_factor_solve,
       sparse_network_admittance_matrix,
       sparse_network_factorization_update,
       sparse_network_factorization_update!,
       grouped_sparse_network_factorization_update!,
       sparse_network_solution,
       grouped_sparse_network_solution!,
       sparse_network_solution_update!,
       over16_fortran_sparse_factor_update,
       over16_fortran_sparse_factor_update!,
       over16_fortran_sparse_factor_solve,
       over16_network_solution_update!,
       over16_switch_current_reconstruction,
       over16_switch_current_reconstruction_table,
       over16_switch_current_reconstruction_table!,
       over16_switch_post_current_transition,
       over16_switch_post_current_transition_table,
       over16_switch_post_current_transition_table!,
       over16_switch_alteration_rebuild_intent,
       over16_switch_alteration_rebuild_update!,
       over16_switch_bvalue_export,
       over16_switch_bvalue_export!

struct IdealSwitch <: EMTElement
    a::Int
    b::Int
    closed::Bool
    on_conductance::Float64
    off_conductance::Float64
end

function IdealSwitch(a::Int, b::Int, closed::Bool; on_conductance::Real=1.0e9,
                     off_conductance::Real=0.0)
    return IdealSwitch(a, b, closed, Float64(on_conductance), Float64(off_conductance))
end

mutable struct CurrentExtinctionState
    not_before_time_s::Float64
    critical_current_a::Float64
    closed::Bool
    opened::Bool
    current_initialized::Bool
    previous_current::Float64
    operation_count::Int
    open_reason::Symbol
    opened_time_s::Float64
end

mutable struct TimeSwitch <: EMTElement
    a::Int
    b::Int
    close_time_s::Float64
    open_time_s::Float64
    initially_closed::Bool
    on_conductance::Float64
    off_conductance::Float64
    current_extinction::Union{Nothing,CurrentExtinctionState}
end

mutable struct CurrentZeroSwitch <: EMTElement
    a::Int
    b::Int
    close_time_s::Float64
    open_request_time_s::Float64
    open_delay_time_s::Float64
    critical_current_a::Float64
    initially_closed::Bool
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    opened::Bool
    current_initialized::Bool
    previous_current::Float64
    operation_count::Int
    open_reason::Symbol
end

function CurrentZeroSwitch(
    switch::TimeSwitch;
    open_delay_time_s::Real = 0.0,
    critical_current_a::Real = 0.0,
)
    delay_time = Float64(open_delay_time_s)
    critical_current = Float64(critical_current_a)
    delay_time >= 0.0 && (isfinite(delay_time) || delay_time == Inf) ||
        throw(ArgumentError("current-extinction delay time must be nonnegative"))
    critical_current >= 0.0 && isfinite(critical_current) ||
        throw(ArgumentError("critical current must be finite and nonnegative"))
    initially_active = switch_closed(switch, 0.0)
    return CurrentZeroSwitch(
        switch.a,
        switch.b,
        switch.close_time_s,
        switch.open_time_s,
        delay_time,
        critical_current,
        switch.initially_closed,
        switch.on_conductance,
        switch.off_conductance,
        initially_active,
        false,
        false,
        0.0,
        0,
        :none,
    )
end

mutable struct OVER16SwitchScanState
    positions::Vector{Int}
    elapsed_open_times::Vector{Float64}
    gap_currents::Vector{Float64}
    modswt::Vector{Int}
end

function OVER16SwitchScanState(
    positions::AbstractVector{Int},
    elapsed_open_times::AbstractVector{<:Real};
    gap_currents::AbstractVector{<:Real}=Float64[],
    modswt::AbstractVector{Int}=Int[],
)
    length(positions) == length(elapsed_open_times) ||
        throw(ArgumentError("positions and elapsed_open_times lengths must match"))
    return OVER16SwitchScanState(
        collect(positions),
        Float64.(elapsed_open_times),
        Float64.(gap_currents),
        collect(modswt),
    )
end

mutable struct OVER16SwitchOperationState
    modswt::Vector{Int}
    closed_switch_count::Int
    accumulated_operation_count::Int
end

function OVER16SwitchOperationState(
    modswt::AbstractVector{Int},
    closed_switch_count::Int;
    accumulated_operation_count::Int=0,
)
    closed_switch_count >= 0 ||
        throw(ArgumentError("closed_switch_count must be nonnegative"))
    accumulated_operation_count >= 0 ||
        throw(ArgumentError("accumulated_operation_count must be nonnegative"))
    return OVER16SwitchOperationState(
        collect(modswt),
        closed_switch_count,
        accumulated_operation_count,
    )
end

mutable struct OVER16SwitchTopologyState
    closed_mask::Vector{Bool}
    closed_switch_count::Int
    first_group_head::Int
    nextsw::Vector{Int}
    kode::Vector{Int}
end

function OVER16SwitchTopologyState(
    closed_mask::AbstractVector{Bool};
    closed_switch_count::Int=count(identity, closed_mask),
    first_group_head::Int=0,
    nextsw::AbstractVector{Int}=zeros(Int, length(closed_mask)),
    kode::AbstractVector{Int}=Int[],
)
    switch_count = length(closed_mask)
    closed_switch_count == count(identity, closed_mask) ||
        throw(ArgumentError("closed_switch_count must match closed_mask"))
    0 <= first_group_head <= switch_count ||
        throw(ArgumentError("first_group_head must be between 0 and switch count"))
    length(nextsw) == switch_count ||
        throw(ArgumentError("nextsw length must match closed_mask"))
    return OVER16SwitchTopologyState(
        collect(closed_mask),
        closed_switch_count,
        first_group_head,
        collect(nextsw),
        collect(kode),
    )
end

mutable struct OVER16SwitchAdmittanceState
    base_admittance::Matrix{Float64}
    admittance::Matrix{Float64}
    admittance_workspace::Matrix{Float64}
    switch_conductances::Vector{Float64}
    retriangularization_count::Int
end

struct SwitchOperationStepResult
    processed_modswt::Vector{Int}
    ktrlsw_count::Int
    switch_operation_state_mutated::Bool
end

struct SwitchStatusStepResult
    requires_order_rebuild::Bool
    switch_status_state_mutated::Bool
end

struct SwitchOrderStepResult
    switch_order_state_mutated::Bool
end

struct SwitchAdmittanceStepResult
    should_retriangularize::Bool
    admittance_mutated::Bool
    switch_admittance_state_mutated::Bool
end

struct SwitchTopologyAdmittanceStepResult
    status_result::SwitchStatusStepResult
    order_result::Union{Nothing,SwitchOrderStepResult}
    admittance_result::SwitchAdmittanceStepResult
    switch_topology_state_mutated::Bool
    switch_admittance_state_mutated::Bool
    switch_topology_admittance_state_mutated::Bool
    topology_mutated::Bool
    admittance_mutated::Bool
    should_retriangularize::Bool
end

function OVER16SwitchAdmittanceState(
    base_admittance::AbstractMatrix{<:Real};
    switch_conductances::AbstractVector{<:Real}=Float64[],
    retriangularization_count::Int=0,
)
    size(base_admittance, 1) == size(base_admittance, 2) ||
        throw(ArgumentError("base_admittance must be square"))
    retriangularization_count >= 0 ||
        throw(ArgumentError("retriangularization_count must be nonnegative"))
    base = Float64.(base_admittance)
    for value in base
        isfinite(value) ||
            throw(ArgumentError("base_admittance entries must be finite"))
    end
    conductances = Float64.(switch_conductances)
    for value in conductances
        isfinite(value) && value >= 0.0 ||
            throw(ArgumentError("switch_conductances entries must be finite and nonnegative"))
    end
    return OVER16SwitchAdmittanceState(
        base,
        copy(base),
        copy(base),
        conductances,
        retriangularization_count,
    )
end

function OVER16SwitchAdmittanceState(
    node_count::Int;
    switch_count::Int=0,
    retriangularization_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    switch_count >= 0 || throw(ArgumentError("switch_count must be nonnegative"))
    return OVER16SwitchAdmittanceState(
        zeros(Float64, node_count, node_count);
        switch_conductances = zeros(Float64, switch_count),
        retriangularization_count = retriangularization_count,
    )
end

mutable struct OVER16SwitchRetriangularizationState
    factor::Matrix{Float64}
    pivot_values::Vector{Float64}
    factorization_count::Int
end

function OVER16SwitchRetriangularizationState(
    factor::AbstractMatrix{<:Real};
    pivot_values::AbstractVector{<:Real}=Float64[],
    factorization_count::Int=0,
)
    size(factor, 1) == size(factor, 2) ||
        throw(ArgumentError("factor must be square"))
    factorization_count >= 0 ||
        throw(ArgumentError("factorization_count must be nonnegative"))
    dense_factor = Float64.(factor)
    for value in dense_factor
        isfinite(value) ||
            throw(ArgumentError("factor entries must be finite"))
    end
    pivots = isempty(pivot_values) ? zeros(Float64, size(factor, 1)) : Float64.(pivot_values)
    length(pivots) == size(factor, 1) ||
        throw(ArgumentError("pivot_values length must match factor size"))
    for value in pivots
        isfinite(value) ||
            throw(ArgumentError("pivot_values entries must be finite"))
    end
    return OVER16SwitchRetriangularizationState(dense_factor, pivots, factorization_count)
end

function OVER16SwitchRetriangularizationState(
    node_count::Int;
    factorization_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return OVER16SwitchRetriangularizationState(
        zeros(Float64, node_count, node_count);
        pivot_values = zeros(Float64, node_count),
        factorization_count = factorization_count,
    )
end

mutable struct OVER16SwitchSparseFactorWorkspaceState
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    workspace_update_count::Int
end

function OVER16SwitchSparseFactorWorkspaceState(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int};
    workspace_update_count::Int=0,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    workspace_update_count >= 0 ||
        throw(ArgumentError("workspace_update_count must be nonnegative"))
    dense_ykm = Float64.(ykm)
    for value in dense_ykm
        isfinite(value) ||
            throw(ArgumentError("ykm entries must be finite"))
    end
    return OVER16SwitchSparseFactorWorkspaceState(
        collect(km),
        dense_ykm,
        collect(kk),
        workspace_update_count,
    )
end

function OVER16SwitchSparseFactorWorkspaceState(
    node_count::Int;
    workspace_update_count::Int=0,
)
    node_count > 0 || throw(ArgumentError("node_count must be positive"))
    return OVER16SwitchSparseFactorWorkspaceState(
        Int[],
        Float64[],
        zeros(Int, node_count);
        workspace_update_count = workspace_update_count,
    )
end

mutable struct OVER16FortranSparseFactorWorkspaceState
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    workspace_update_count::Int
end

"""
Reusable storage for the grouped sparse factorization and solve kernels.

The workspace owns every temporary whose size is determined by the augmented
network. A caller may refactor it whenever admittance or node grouping changes,
then solve any number of right-hand sides without allocating.
"""
mutable struct GroupedSparseNetworkWorkspace
    km::Vector{Int}
    ykm::Vector{Float64}
    kk::Vector{Int}
    pivot_values::Vector{Float64}
    inverse_diagonal_values::Vector{Float64}
    contracted_admittance::Matrix{Float64}
    node_group_representatives::Vector{Int}
    active_rows::Vector{Int}
    row_starts::Vector{Int}
    factor_row::Vector{Float64}
    visit_marks::Vector{Int}
    visit_positions::Vector{Int}
    visit_path::Vector{Int}
    alias_successors::Vector{Int}
    solution::Vector{Float64}
    factorization_count::Int
end

function GroupedSparseNetworkWorkspace(node_count::Int)
    node_count > 1 || throw(ArgumentError(
        "grouped sparse workspace must include a reference and at least one node",
    ))
    km = Int[]
    ykm = Float64[]
    active_rows = Int[]
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]
    sizehint!(km, node_count * node_count)
    sizehint!(ykm, node_count * node_count)
    sizehint!(active_rows, node_count)
    sizehint!(pivot_values, node_count)
    sizehint!(inverse_diagonal_values, node_count)
    return GroupedSparseNetworkWorkspace(
        km,
        ykm,
        zeros(Int, node_count),
        pivot_values,
        inverse_diagonal_values,
        zeros(Float64, node_count, node_count),
        collect(1:node_count),
        active_rows,
        zeros(Int, node_count),
        zeros(Float64, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Int, node_count),
        zeros(Float64, node_count),
        0,
    )
end

function OVER16FortranSparseFactorWorkspaceState(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int};
    workspace_update_count::Int=0,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    workspace_update_count >= 0 ||
        throw(ArgumentError("workspace_update_count must be nonnegative"))
    values = Float64.(ykm)
    for value in values
        isfinite(value) ||
            throw(ArgumentError("ykm entries must be finite"))
    end
    return OVER16FortranSparseFactorWorkspaceState(
        collect(km),
        values,
        collect(kk),
        workspace_update_count,
    )
end

function OVER16FortranSparseFactorWorkspaceState(
    node_count::Int;
    workspace_update_count::Int=0,
)
    node_count > 1 || throw(ArgumentError("node_count must be greater than one"))
    return OVER16FortranSparseFactorWorkspaceState(
        Int[],
        Float64[],
        zeros(Int, node_count);
        workspace_update_count = workspace_update_count,
    )
end

mutable struct OVER16SwitchCurrentState
    rhs::Vector{Float64}
    switch_currents::Vector{Float64}
    current_products::Vector{Float64}
    network_solution::Vector{Float64}
end

function OVER16SwitchCurrentState(
    rhs::AbstractVector{<:Real},
    switch_currents::AbstractVector{<:Real};
    current_products::AbstractVector{<:Real}=zeros(Float64, length(switch_currents)),
    network_solution::AbstractVector{<:Real}=Float64[],
)
    length(current_products) == length(switch_currents) ||
        throw(ArgumentError("current_products length must match switch_currents"))
    isempty(network_solution) || length(network_solution) == length(rhs) ||
        throw(ArgumentError("network_solution length must match rhs when provided"))
    return OVER16SwitchCurrentState(
        Float64.(rhs),
        Float64.(switch_currents),
        Float64.(current_products),
        Float64.(network_solution),
    )
end

mutable struct OVER16SwitchPostCurrentState
    positions::Vector{Int}
    switch_currents::Vector{Float64}
    energies::Vector{Float64}
    modswt::Vector{Int}
end

function OVER16SwitchPostCurrentState(
    positions::AbstractVector{Int},
    switch_currents::AbstractVector{<:Real},
    energies::AbstractVector{<:Real};
    modswt::AbstractVector{Int}=Int[],
)
    switch_count = length(positions)
    length(switch_currents) == switch_count ||
        throw(ArgumentError("switch_currents length must match positions"))
    length(energies) == switch_count ||
        throw(ArgumentError("energies length must match positions"))
    return OVER16SwitchPostCurrentState(
        collect(positions),
        Float64.(switch_currents),
        Float64.(energies),
        collect(modswt),
    )
end

mutable struct OVER16SwitchBValueExportState
    bvalue::Vector{Float64}
    output_count::Int
end

function OVER16SwitchBValueExportState(
    bvalue::AbstractVector{<:Real};
    output_count::Int=length(bvalue),
)
    output_count >= 0 || throw(ArgumentError("output_count must be nonnegative"))
    output_count <= length(bvalue) ||
        throw(ArgumentError("output_count cannot exceed bvalue length"))
    return OVER16SwitchBValueExportState(Float64.(bvalue), output_count)
end

mutable struct OVER16SwitchAlterationState
    ialter::Int
    operation_count::Int
    closed_switch_count::Int
    triangularization_count::Int
    first_group_head::Int
    total_operation_count::Int
    ktrlsw6::Int
end

function OVER16SwitchAlterationState(
    ialter::Int,
    operation_count::Int,
    closed_switch_count::Int,
    triangularization_count::Int;
    first_group_head::Int=0,
    total_operation_count::Int=0,
    ktrlsw6::Int=1,
)
    over16_switch_alteration_rebuild_intent(
        ialter,
        operation_count,
        closed_switch_count,
        triangularization_count,
        first_group_head,
        total_operation_count,
        ktrlsw6,
    )
    return OVER16SwitchAlterationState(
        ialter,
        operation_count,
        closed_switch_count,
        triangularization_count,
        first_group_head,
        total_operation_count,
        ktrlsw6,
    )
end

function TimeSwitch(a::Int, b::Int; close_time_s::Real=Inf, open_time_s::Real=Inf,
                    initially_closed::Bool=false, on_conductance::Real=1.0e9,
                    off_conductance::Real=0.0)
    return TimeSwitch(
        a,
        b,
        Float64(close_time_s),
        Float64(open_time_s),
        initially_closed,
        Float64(on_conductance),
        Float64(off_conductance),
        nothing,
    )
end

switch_closed(s::IdealSwitch, _t::Real)::Bool = s.closed

function switch_closed(s::TimeSwitch, t::Real)::Bool
    s.current_extinction === nothing ||
        return s.current_extinction.closed
    time = Float64(t)
    closed = s.initially_closed
    if isfinite(s.close_time_s) && time >= s.close_time_s
        closed = true
    end
    if isfinite(s.open_time_s) && time >= s.open_time_s
        closed = false
    end
    return closed
end

function switch_conductance(s::Union{IdealSwitch,TimeSwitch}, t::Real)::Float64
    return switch_closed(s, t) ? s.on_conductance : s.off_conductance
end

current_extinction_enabled(::Any) = false
current_extinction_enabled(::CurrentZeroSwitch) = true
current_extinction_enabled(s::TimeSwitch) = s.current_extinction !== nothing
current_extinction_enabled(::Tuple{}) = false
current_extinction_enabled(elements::Tuple) =
    any(current_extinction_enabled, elements)

function configure_current_extinction!(
    switch::TimeSwitch,
    not_before_time_s::Real,
    critical_current_a::Real,
    time_s::Real,
    ;
    currently_closed::Bool = switch_closed(switch, time_s),
)
    not_before = Float64(not_before_time_s)
    critical_current = Float64(critical_current_a)
    time = Float64(time_s)
    not_before >= 0.0 && (isfinite(not_before) || not_before == Inf) ||
        throw(ArgumentError("current-extinction not-before time must be nonnegative"))
    critical_current >= 0.0 && isfinite(critical_current) ||
        throw(ArgumentError("critical current must be finite and nonnegative"))
    if not_before == 0.0 && critical_current == 0.0
        switch.current_extinction = nothing
        return switch
    end
    scheduled_closed =
        currently_closed ||
        (
            isfinite(switch.close_time_s) &&
            time >= switch.close_time_s &&
            time < switch.open_time_s
        )
    previous = switch.current_extinction
    reclosed = !currently_closed && scheduled_closed
    switch.current_extinction = CurrentExtinctionState(
        not_before,
        critical_current,
        scheduled_closed,
        reclosed || previous === nothing ? false : previous.opened,
        reclosed || previous === nothing ? false : previous.current_initialized,
        reclosed || previous === nothing ? 0.0 : previous.previous_current,
        (previous === nothing ? 0 : previous.operation_count) + (reclosed ? 1 : 0),
        reclosed ? :restart_topology_change :
            previous === nothing ? :none : previous.open_reason,
        reclosed || previous === nothing ? Inf : previous.opened_time_s,
    )
    return switch
end

switch_closed(s::CurrentZeroSwitch, _t::Real)::Bool = s.closed

function switch_conductance(s::CurrentZeroSwitch, _t::Real)::Float64
    return s.closed ? s.on_conductance : s.off_conductance
end

function prepare_current_zero_switch!(switch::CurrentZeroSwitch, time_s::Real)
    time = Float64(time_s)
    if !switch.opened && !switch.closed &&
       isfinite(switch.close_time_s) && time >= switch.close_time_s
        switch.closed = true
        switch.operation_count += 1
    end
    return switch
end

function prepare_current_zero_switch!(switch::TimeSwitch, time_s::Real)
    state = switch.current_extinction
    state === nothing && return switch
    time = Float64(time_s)
    if !state.opened && !state.closed &&
       isfinite(switch.close_time_s) && time >= switch.close_time_s
        state.closed = true
        state.operation_count += 1
    end
    return switch
end

function _check_current_zero_transition_reason(reason::Symbol)
    reason in (:current_reversal, :critical_current) || throw(ArgumentError(
        "current-zero transition reason must be :current_reversal or :critical_current",
    ))
    return reason
end

function apply_current_zero_transition!(
    switch::CurrentZeroSwitch,
    reason::Symbol,
    time_s::Real,
)
    _check_current_zero_transition_reason(reason)
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("current-zero transition time must be finite"))
    if switch.closed
        switch.closed = false
        switch.opened = true
        switch.operation_count += 1
        switch.open_reason = reason
    end
    return switch
end

function apply_current_zero_transition!(
    switch::TimeSwitch,
    reason::Symbol,
    time_s::Real,
)
    state = switch.current_extinction
    state === nothing && throw(ArgumentError(
        "time switch has no current-extinction transition owner",
    ))
    _check_current_zero_transition_reason(reason)
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("current-zero transition time must be finite"))
    if state.closed
        state.closed = false
        state.opened = true
        state.operation_count += 1
        state.open_reason = reason
        state.opened_time_s = time
    end
    return switch
end

function update_current_zero_switch!(
    switch::CurrentZeroSwitch,
    current::Real,
    time_s::Real,
)
    switch.closed || return switch
    current_value = Float64(current)
    time = Float64(time_s)
    current_reversed =
        switch.current_initialized &&
        current_value * switch.previous_current < 0.0
    critical_reached =
        switch.critical_current_a > 0.0 &&
        abs(current_value) < switch.critical_current_a
    if switch.current_initialized &&
       time >= switch.open_request_time_s &&
       time >= switch.open_delay_time_s &&
       (current_reversed || critical_reached)
        apply_current_zero_transition!(
            switch,
            critical_reached ? :critical_current : :current_reversal,
            time,
        )
    end
    switch.previous_current = current_value
    switch.current_initialized = true
    return switch
end

function update_current_zero_switch!(
    switch::TimeSwitch,
    current::Real,
    time_s::Real,
)
    state = switch.current_extinction
    state === nothing && return switch
    state.closed || return switch
    current_value = Float64(current)
    time = Float64(time_s)
    current_reversed =
        state.current_initialized &&
        current_value * state.previous_current < 0.0
    critical_reached =
        state.critical_current_a > 0.0 &&
        abs(current_value) < state.critical_current_a
    if state.current_initialized &&
       time >= switch.open_time_s &&
       time >= state.not_before_time_s &&
       (current_reversed || critical_reached)
        apply_current_zero_transition!(
            switch,
            critical_reached ? :critical_current : :current_reversal,
            time,
        )
    end
    state.previous_current = current_value
    state.current_initialized = true
    return switch
end

function over16_switch_margin_history(
    previous_ck::Float64,
    voltage_difference::Float64,
    delta2::Float64,
    critical_current::Float64,
    branch_multiplier::Float64,
    cik_current::Float64,
    cik_next::Float64,
)
    area = voltage_difference * delta2
    midpoint_ck = previous_ck + area
    updated_ck = midpoint_ck + area
    should_update = !(abs(midpoint_ck) >= critical_current || abs(midpoint_ck) > abs(previous_ck))
    if !should_update
        return updated_ck, cik_current, cik_next, false
    end

    history_current = branch_multiplier * voltage_difference
    return updated_ck, -history_current, history_current + cik_current + cik_next, true
end

function over16_switch_close_critical_current(
    previous_ck::Float64,
    voltage_difference::Float64,
    delta2::Float64,
    open_threshold::Float64,
    branch_multiplier::Float64,
    branch_delay_multiplier::Float64,
    cik_next::Float64,
    critical_current::Float64,
)
    delta2 > 0.0 || throw(ArgumentError("delta2 must be positive"))
    branch_delay_multiplier != 0.0 || throw(ArgumentError("branch_delay_multiplier must be nonzero"))

    history_current = branch_multiplier * voltage_difference
    area = voltage_difference * delta2
    midpoint_ck = previous_ck + area
    updated_ck = midpoint_ck + area
    should_update = !(abs(midpoint_ck) < open_threshold || abs(midpoint_ck) < abs(previous_ck))
    if !should_update
        return updated_ck, -history_current, critical_current, false
    end

    d2 = branch_delay_multiplier / delta2
    multiplier_ratio = branch_multiplier / (delta2 * d2)
    multiplier_ratio != 0.0 || throw(ArgumentError("branch_multiplier must be nonzero"))
    ci1 = cik_next + d2 * area
    d3 = midpoint_ck * (multiplier_ratio + 1.0) - ci1 / d2
    return updated_ck, -history_current, abs(d3 / multiplier_ratio), true
end

function over16_gap_energy_check(previous_energy::Float64, voltage_difference::Float64)
    average_voltage = (voltage_difference + previous_energy) / 2.0
    return average_voltage, voltage_difference, average_voltage >= 0.0
end

function over16_voltage_open_delay(
    voltage_difference::Float64,
    open_threshold::Float64,
    t::Float64,
    delay_offset::Float64,
    current_delay::Float64,
)
    if abs(voltage_difference) < open_threshold
        return current_delay, false
    end
    return t + delay_offset, true
end

function over16_switch_tail_current_injection(
    decay_end_time::Float64,
    t::Float64,
    dt::Float64,
    amplitude::Float64,
    time_constant::Float64,
    cutoff_current::Float64,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    time_constant != 0.0 || throw(ArgumentError("time_constant must be nonzero"))

    decay_factor = exp((decay_end_time - t - dt) / time_constant)
    current = amplitude * decay_factor
    clear_marker = current <= cutoff_current ||
                   (decay_factor <= 1.0e-4 && cutoff_current == 0.0)
    return current, -current, decay_factor, clear_marker
end

function over16_open_switch_close_decision(
    elapsed_open_time::Float64,
    dt::Float64,
    maximum_open_time::Float64,
    voltage_difference::Float64,
    close_threshold::Float64;
    is_gap::Bool=false,
    control_enabled::Bool=false,
    control_signal::Float64=0.0,
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))

    updated_elapsed = elapsed_open_time + dt
    if updated_elapsed > maximum_open_time
        updated_elapsed = -fltinf
    end

    threshold_voltage = is_gap ? abs(voltage_difference) : voltage_difference
    if threshold_voltage < close_threshold
        return updated_elapsed, false
    end

    if !is_gap && updated_elapsed >= 0.0
        return -fltinf, true
    end

    if !control_enabled
        return updated_elapsed, true
    end
    return updated_elapsed, control_signal > 10.0 * flzero
end

function over16_closed_switch_open_decision(
    switch_current::Float64,
    previous_gap_current::Float64,
    open_threshold::Float64;
    is_gap::Bool=false,
    gap_control_signal::Float64=0.0,
    flzero::Float64=0.0,
)
    if is_gap
        if gap_control_signal > 10.0 * flzero
            return previous_gap_current, false
        end
        updated_gap_current = switch_current
        if switch_current * previous_gap_current < 0.0
            return updated_gap_current, true
        end
        return updated_gap_current, abs(switch_current) < open_threshold
    end
    return previous_gap_current, switch_current < open_threshold
end

function over16_switch_clamp_action(
    switch_type::Int,
    position::Int,
    controlled_node::Int,
    control_index::Int,
    control_signal::Float64,
    flzero::Float64,
)
    m1 = abs(position)
    tolerance = 10.0 * flzero
    if control_index != 0
        if control_signal < -tolerance
            return m1 == 2 ? :open : :skip
        elseif control_signal > tolerance
            return m1 == 2 ? :skip : :close
        end
    elseif switch_type == 8891 || (switch_type == 8890 && controlled_node == 0)
        return :skip
    end

    if switch_type == 8891
        return m1 == 2 ? :open : :skip
    elseif switch_type == 8890 && controlled_node == 0
        return :skip
    elseif m1 == 2
        return :evaluate_closed
    end
    return :evaluate_open
end

function over16_switch_position_update(
    position::Int,
    action::Symbol,
    elapsed_open_time::Float64;
    reset_elapsed::Bool=false,
)
    if action == :open
        magnitude = 5
        modifier = -1
    elseif action == :close
        magnitude = 2
        modifier = 1
    else
        throw(ArgumentError("action must be :open or :close"))
    end
    signed_position = position < 0 ? -magnitude : magnitude
    updated_elapsed = action == :open && reset_elapsed ? 0.0 : elapsed_open_time
    return signed_position, modifier, updated_elapsed
end

function over16_a8sw_delayed_open_step(
    switch_current::Float64,
    previous_current::Float64,
    shape_current::Float64,
    shape_delay::Float64,
    scheduled_open_time::Float64,
    t::Float64,
    dt::Float64,
    coefficient::Float64,
    exponent::Float64,
    time_scale::Float64;
    threshold_marker::Float64=9999.0,
)
    dt > 0.0 || throw(ArgumentError("dt must be positive"))
    time_scale != 0.0 || throw(ArgumentError("time_scale must be nonzero"))

    updated_shape_current = shape_current
    updated_shape_delay = shape_delay
    updated_scheduled_open_time = scheduled_open_time

    if switch_current > 0.0
        delta_current = previous_current - switch_current
        if delta_current < switch_current
            return (
                previous_current = switch_current,
                shape_current = updated_shape_current,
                shape_delay = updated_shape_delay,
                scheduled_open_time = updated_scheduled_open_time,
                open_threshold = nothing,
                opens = false,
            )
        end
        didt = delta_current / dt
        updated_shape_current = didt^exponent * coefficient
        updated_shape_delay = updated_shape_current / (didt * time_scale)
        updated_scheduled_open_time =
            t + switch_current * dt / delta_current + updated_shape_delay * time_scale
    elseif previous_current != 0.0 && previous_current * switch_current <= 0.0
        delta_current = previous_current - switch_current
        didt = delta_current / dt
        updated_shape_current = didt^exponent * coefficient
        updated_shape_delay = updated_shape_current / (didt * time_scale)
        updated_scheduled_open_time =
            t + switch_current * dt / delta_current + updated_shape_delay * time_scale
    end

    if t + dt < updated_scheduled_open_time
        return (
            previous_current = switch_current,
            shape_current = updated_shape_current,
            shape_delay = updated_shape_delay,
            scheduled_open_time = updated_scheduled_open_time,
            open_threshold = nothing,
            opens = false,
        )
    end

    return (
        previous_current = previous_current,
        shape_current = updated_shape_current,
        shape_delay = updated_shape_delay,
        scheduled_open_time = updated_scheduled_open_time,
        open_threshold = threshold_marker,
        opens = true,
    )
end

function over16_a8sw_tail_current_table(
    from_nodes::AbstractVector{Int},
    to_nodes::AbstractVector{Int},
    switch_types::AbstractVector{Int},
    positions::AbstractVector{Int},
    tail_indices::AbstractVector{Int},
    tail_markers::AbstractVector{<:Real},
    scheduled_open_times::AbstractVector{<:Real},
    amplitudes::AbstractVector{<:Real},
    time_constants::AbstractVector{<:Real},
    cutoff_currents::AbstractVector{<:Real},
    t::Real,
    dt::Real;
    node_count::Int=max(maximum(vcat([1], collect(from_nodes), collect(to_nodes))), 1),
    tail_marker_value::Real=9999.0,
)
    switch_count = length(switch_types)
    _over16_check_length("from_nodes", from_nodes, switch_count)
    _over16_check_length("to_nodes", to_nodes, switch_count)
    _over16_check_length("positions", positions, switch_count)
    _over16_check_length("tail_indices", tail_indices, switch_count)
    _over16_check_length("tail_markers", tail_markers, switch_count)
    _over16_check_length("scheduled_open_times", scheduled_open_times, switch_count)
    _over16_check_length("amplitudes", amplitudes, switch_count)
    _over16_check_length("time_constants", time_constants, switch_count)
    _over16_check_length("cutoff_currents", cutoff_currents, switch_count)
    node_count >= 1 || throw(ArgumentError("node_count must be positive"))

    rhs_assignments = zeros(Float64, node_count)
    tail_currents = zeros(Float64, switch_count)
    decay_factors = zeros(Float64, switch_count)
    updated_tail_markers = Float64.(tail_markers)
    active_rows = Int[]
    cleared_rows = Int[]

    for row in 1:switch_count
        from_node = from_nodes[row]
        to_node = to_nodes[row]
        1 <= from_node <= node_count ||
            throw(ArgumentError("from_nodes entries must be within node_count"))
        1 <= to_node <= node_count ||
            throw(ArgumentError("to_nodes entries must be within node_count"))
        from_node != to_node ||
            throw(ArgumentError("tail-current switch endpoints must be distinct"))

        switch_type = switch_types[row]
        if !(switch_type == 8888 || switch_type == 8890 || switch_type == 8891)
            continue
        end
        if tail_indices[row] <= 0 || abs(positions[row]) != 5 ||
           Float64(tail_markers[row]) != Float64(tail_marker_value)
            continue
        end

        current_from, current_to, decay_factor, clear_marker =
            over16_switch_tail_current_injection(
                Float64(scheduled_open_times[row]),
                Float64(t),
                Float64(dt),
                Float64(amplitudes[row]),
                Float64(time_constants[row]),
                Float64(cutoff_currents[row]),
            )
        rhs_assignments[from_node] = current_from
        rhs_assignments[to_node] = current_to
        tail_currents[row] = current_from
        decay_factors[row] = decay_factor
        push!(active_rows, row)
        if clear_marker
            updated_tail_markers[row] = 0.0
            push!(cleared_rows, row)
        end
    end

    return (
        rhs_assignments = rhs_assignments,
        tail_currents = tail_currents,
        decay_factors = decay_factors,
        updated_tail_markers = updated_tail_markers,
        active_rows = active_rows,
        cleared_rows = cleared_rows,
        rhs_mutated = false,
        topology_mutated = false,
        admittance_mutated = false,
    )
end

function over16_controlled_switch_scan_step(
    switch_type::Int,
    position::Int,
    elapsed_open_time::Float64,
    maximum_open_time::Float64,
    voltage_difference::Float64,
    switch_current::Float64,
    close_threshold::Float64,
    open_threshold::Float64,
    dt::Float64;
    controlled_node::Int=0,
    control_index::Int=0,
    clamp_signal::Float64=0.0,
    close_control_signal::Float64=0.0,
    gap_control_signal::Float64=0.0,
    previous_gap_current::Float64=0.0,
    flzero::Float64=0.0,
    fltinf::Float64=Inf,
)
    scan_action = over16_switch_clamp_action(
        switch_type,
        position,
        controlled_node,
        control_index,
        clamp_signal,
        flzero,
    )

    transition = :none
    updated_elapsed = elapsed_open_time
    updated_gap_current = previous_gap_current
    if scan_action == :skip
        transition = :none
    elseif scan_action == :open || scan_action == :close
        transition = scan_action
    elseif scan_action == :evaluate_closed
        updated_gap_current, opens = over16_closed_switch_open_decision(
            switch_current,
            previous_gap_current,
            open_threshold;
            is_gap = switch_type == 8890,
            gap_control_signal = gap_control_signal,
            flzero = flzero,
        )
        transition = opens ? :open : :none
    elseif scan_action == :evaluate_open
        updated_elapsed, closes = over16_open_switch_close_decision(
            elapsed_open_time,
            dt,
            maximum_open_time,
            voltage_difference,
            close_threshold;
            is_gap = switch_type == 8890,
            control_enabled = controlled_node != 0,
            control_signal = close_control_signal,
            flzero = flzero,
            fltinf = fltinf,
        )
        transition = closes ? :close : :none
    else
        throw(ArgumentError("unsupported switch scan action"))
    end

    if transition == :none
        return (
            scan_action = scan_action,
            transition = transition,
            position = position,
            elapsed_open_time = updated_elapsed,
            gap_current = updated_gap_current,
            modswt_sign = 0,
            altered = false,
        )
    end

    updated_position, modswt_sign, final_elapsed = over16_switch_position_update(
        position,
        transition,
        updated_elapsed;
        reset_elapsed = transition == :open && maximum_open_time != 0.0,
    )
    return (
        scan_action = scan_action,
        transition = transition,
        position = updated_position,
        elapsed_open_time = final_elapsed,
        gap_current = updated_gap_current,
        modswt_sign = modswt_sign,
        altered = true,
    )
end

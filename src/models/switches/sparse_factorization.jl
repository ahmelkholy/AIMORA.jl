const OVER16_SPARSE_FACTOR_LABELS = (
    2205, 2220, 2227, 2237, 2229, 2230, 2231, 2232, 2233, 2240,
    2260, 2265, 2270, 2278, 2272, 2273, 2275, 2276, 2277, 2279,
    2280, 4312, 2283, 2288, 2285, 2290,
)

function sparse_network_admittance_matrix(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    row_end_pointers::AbstractVector{Int};
    node_count::Int=length(row_end_pointers),
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count == length(row_end_pointers) ||
        throw(ArgumentError("node_count must match row_end_pointers length"))
    node_count > 1 ||
        throw(ArgumentError("row_end_pointers must include reference and at least one sparse row"))
    previous_end = 0
    admittance = zeros(Float64, node_count, node_count)
    for row in 1:node_count
        one_past_row_end = row_end_pointers[row]
        1 <= one_past_row_end <= length(km) + 1 ||
            throw(ArgumentError("row_end_pointers entries must point one past a sparse row"))
        row_end = one_past_row_end - 1
        previous_end <= row_end ||
            throw(ArgumentError("row_end_pointers entries must be nondecreasing"))
        for index in (previous_end + 1):row_end
            column = abs(km[index])
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            value = Float64(ykm[index])
            isfinite(value) || throw(ArgumentError("ykm entries must be finite"))
            admittance[row, column] += value
        end
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("row_end_pointers must consume every km/ykm entry"))
    return admittance
end

function _check_fortran_sparse_workspace(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    first_factor_row::Int,
    partition_boundary::Int,
)
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(kk)
    node_count > 1 || throw(ArgumentError("kk must cover at least the reference and one factor row"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    all(index -> index == 0, kk[1:(first_factor_row - 1)]) ||
        throw(ArgumentError("kk entries before first_factor_row must be zero"))
    previous_end = 0
    for row in first_factor_row:partition_boundary
        row_end = kk[row]
        previous_end < row_end <= length(km) ||
            throw(ArgumentError("kk entries must provide nonempty increasing factor rows"))
        row_start = previous_end + 1
        km[row_start] == -row ||
            throw(ArgumentError("each Fortran sparse factor row must start with its negative diagonal marker"))
        for index in row_start:row_end
            marker = km[index]
            marker != 0 || throw(ArgumentError("km entries must not be zero"))
            column = abs(marker)
            1 <= column <= node_count ||
                throw(ArgumentError("km column entries must be within node count"))
            isfinite(Float64(ykm[index])) ||
                throw(ArgumentError("ykm entries must be finite"))
            if index > row_start
                marker > row ||
                    throw(ArgumentError("off-diagonal Fortran sparse factor entries must point to later columns"))
            end
        end
        previous_end = row_end
    end
    previous_end == length(km) ||
        throw(ArgumentError("kk must consume every km/ykm entry"))
    return nothing
end

function _node_group_representatives(
    node_group_successors::AbstractVector{<:Integer},
    node_count::Int,
)
    length(node_group_successors) >= node_count ||
        throw(ArgumentError("node_group_successors must cover every network node"))
    representatives = collect(1:node_count)
    for start in 1:node_count
        successor = Int(node_group_successors[start])
        successor == 0 && continue
        1 <= successor <= node_count ||
            throw(ArgumentError("node_group_successors entries must be valid node indices"))
        path = Int[]
        seen = Dict{Int,Int}()
        node = start
        while node != 0 && !haskey(seen, node)
            1 <= node <= node_count ||
                throw(ArgumentError("node group chain left the valid node range"))
            seen[node] = length(path) + 1
            push!(path, node)
            node = Int(node_group_successors[node])
        end
        node == 0 && continue
        group = path[seen[node]:end]
        representative = maximum(group)
        for grouped_node in group
            representatives[grouped_node] = representative
        end
    end
    return representatives
end

function _node_group_contracted_admittance(
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer},
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    dense = Float64.(admittance)
    for value in dense
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end
    representatives = _node_group_representatives(node_group_successors, node_count)
    contracted = zeros(Float64, node_count, node_count)
    for column in 1:node_count
        representative_column = representatives[column]
        for row in 1:node_count
            contracted[representatives[row], representative_column] += dense[row, column]
        end
    end
    return contracted, representatives
end

function _grouped_sparse_network_factorization_update(
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 1 || throw(ArgumentError("admittance must include a reference node and factor rows"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    dense, representatives =
        _node_group_contracted_admittance(admittance, node_group_successors)
    active_rows = Int[
        row for row in first_factor_row:partition_boundary
        if representatives[row] == row
    ]
    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    row_starts = zeros(Int, node_count)
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]

    for row in active_rows
        f = zeros(Float64, node_count)
        f[1] = dense[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            representatives[column] == column || continue
            value = dense[row, column]
            if abs(value) > zero_tol
                f[column] = value
            end
        end

        for previous_row in active_rows
            previous_row < row || break
            a = f[previous_row]
            abs(a) > zero_tol || continue
            previous_start = row_starts[previous_row]
            previous_end = kk[previous_row]
            previous_start <= previous_end ||
                throw(ArgumentError("invalid previously built grouped sparse factor row"))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    f[1] -= a * ykm[index]
                else
                    f[column] -= a * ykm[index]
                end
            end
        end

        pivot = f[1]
        abs(pivot) > pivot_tol ||
            throw(ArgumentError("grouped sparse factor pivot is zero or below tolerance"))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        row_starts[row] = length(km) + 1
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            representatives[column] == column || continue
            value = f[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :grouped_sparse_network_factorization_update,
        outcome = :equation_oracle,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1, :SUBTS3),
        fortran_labels = OVER16_SPARSE_FACTOR_LABELS,
        node_count = node_count,
        partition_boundary = partition_boundary,
        first_factor_row = first_factor_row,
        km = km,
        ykm = ykm,
        kk = kk,
        pivot_values = pivot_values,
        inverse_diagonal_values = inverse_diagonal_values,
        pivot_tolerance = pivot_tol,
        zero_tolerance = zero_tol,
        node_group_successors = Int.(node_group_successors[1:node_count]),
        node_group_representatives = representatives,
        sparse_factor_workspace_built = true,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:nonlinear_inverse_columns, :full_timestep_oracle],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_fortran_sparse_factor_update(
    admittance::AbstractMatrix{<:Real};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
    node_group_successors::Union{Nothing,AbstractVector{<:Integer}}=nothing,
)
    if node_group_successors !== nothing
        return _grouped_sparse_network_factorization_update(
            admittance,
            node_group_successors;
            partition_boundary = partition_boundary,
            first_factor_row = first_factor_row,
            pivot_tolerance = pivot_tolerance,
            zero_tolerance = zero_tolerance,
        )
    end
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count > 1 || throw(ArgumentError("admittance must include a reference node and factor rows"))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    dense = Float64.(admittance)
    for value in dense
        isfinite(value) ||
            throw(ArgumentError("admittance entries must be finite"))
    end

    km = Int[]
    ykm = Float64[]
    kk = zeros(Int, node_count)
    pivot_values = Float64[]
    inverse_diagonal_values = Float64[]

    for row in first_factor_row:partition_boundary
        f = zeros(Float64, node_count)
        f[1] = dense[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            value = dense[row, column]
            if abs(value) > zero_tol
                f[column] = value
            end
        end

        for previous_row in first_factor_row:(row - 1)
            a = f[previous_row]
            abs(a) > zero_tol || continue
            previous_start = previous_row == first_factor_row ? 1 : kk[previous_row - 1] + 1
            previous_end = kk[previous_row]
            previous_start <= previous_end ||
                throw(ArgumentError("invalid previously built Fortran sparse factor row"))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    f[1] -= a * ykm[index]
                else
                    f[column] -= a * ykm[index]
                end
            end
        end

        pivot = f[1]
        abs(pivot) > pivot_tol ||
            throw(ArgumentError("Fortran sparse factor pivot is zero or below tolerance"))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            value = f[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end

    return (
        source = :over16_fortran_sparse_factor_update,
        outcome = :equation_oracle,
        fortran_files = (:OVER16_FOR,),
        fortran_routines = (:SUBTS1,),
        fortran_labels = OVER16_SPARSE_FACTOR_LABELS,
        node_count = node_count,
        partition_boundary = partition_boundary,
        first_factor_row = first_factor_row,
        km = km,
        ykm = ykm,
        kk = kk,
        pivot_values = pivot_values,
        inverse_diagonal_values = inverse_diagonal_values,
        pivot_tolerance = pivot_tol,
        zero_tolerance = zero_tol,
        sparse_factor_workspace_built = true,
        fortran_sparse_factor_mutated = false,
        deferred_calls = [:kode_group_ordering, :nonlinear_inverse_columns, :full_timestep_oracle],
        tacs_executed = false,
        solvum_executed = false,
    )
end

function over16_fortran_sparse_factor_update!(
    state::OVER16FortranSparseFactorWorkspaceState,
    admittance::AbstractMatrix{<:Real};
    kwargs...,
)
    km_before = copy(state.km)
    ykm_before = copy(state.ykm)
    kk_before = copy(state.kk)
    count_before = state.workspace_update_count
    preview = over16_fortran_sparse_factor_update(admittance; kwargs...)

    resize!(state.km, length(preview.km))
    state.km .= preview.km
    resize!(state.ykm, length(preview.ykm))
    state.ykm .= preview.ykm
    resize!(state.kk, length(preview.kk))
    state.kk .= preview.kk
    state.workspace_update_count += 1

    workspace_mutated =
        state.km != km_before ||
        state.ykm != ykm_before ||
        state.kk != kk_before ||
        state.workspace_update_count != count_before
    return merge(
        preview,
        (
            km = copy(state.km),
            ykm = copy(state.ykm),
            kk = copy(state.kk),
            workspace_update_count = state.workspace_update_count,
            fortran_sparse_factor_workspace_state_mutated = workspace_mutated,
            fortran_sparse_factor_mutated = workspace_mutated,
        ),
    )
end

function over16_fortran_sparse_factor_solve(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(kk),
)
    _check_fortran_sparse_workspace(km, ykm, kk, first_factor_row, partition_boundary)
    node_count = length(kk)
    length(rhs) == node_count ||
        throw(ArgumentError("rhs length must match kk node count"))
    solution = Float64.(rhs)
    for value in solution
        isfinite(value) || throw(ArgumentError("rhs entries must be finite"))
    end
    solution[1] = 0.0

    for row in first_factor_row:partition_boundary
        row_start = row == first_factor_row ? 1 : kk[row - 1] + 1
        row_end = kk[row]
        a = solution[row]
        solution[row] = a * Float64(ykm[row_start])
        for index in (row_start + 1):row_end
            column = km[index]
            if column <= partition_boundary
                solution[column] -= a * Float64(ykm[index])
            end
        end
    end

    for row in partition_boundary:-1:first_factor_row
        row_start = row == first_factor_row ? 1 : kk[row - 1] + 1
        row_end = kk[row]
        correction = 0.0
        for index in row_end:-1:(row_start + 1)
            column = km[index]
            correction -= solution[column] * Float64(ykm[index])
        end
        solution[row] += correction
    end

    return solution
end

function over16_fortran_sparse_factor_solve(
    state::OVER16FortranSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return over16_fortran_sparse_factor_solve(state.km, state.ykm, state.kk, rhs; kwargs...)
end

function _check_node_group_successors(
    node_group_successors::AbstractVector{<:Integer},
    node_count::Int,
)
    length(node_group_successors) >= node_count ||
        throw(ArgumentError("node_group_successors must cover every network node"))
    for index in 1:node_count
        successor = node_group_successors[index]
        0 <= successor <= node_count ||
            throw(ArgumentError("node_group_successors entries must be zero or a valid node index"))
    end
    return nothing
end

function _node_group_representatives!(
    workspace::GroupedSparseNetworkWorkspace,
    node_group_successors::AbstractVector{<:Integer},
)
    node_count = length(workspace.node_group_representatives)
    _check_node_group_successors(node_group_successors, node_count)
    representatives = workspace.node_group_representatives
    for node in 1:node_count
        representatives[node] = node
    end

    marks = workspace.visit_marks
    positions = workspace.visit_positions
    path = workspace.visit_path
    for start in 1:node_count
        generation = start
        path_length = 0
        node = start
        while node != 0 && marks[node] != generation
            path_length += 1
            path[path_length] = node
            marks[node] = generation
            positions[node] = path_length
            node = Int(node_group_successors[node])
        end
        node == 0 && continue
        cycle_start = positions[node]
        representative = path[cycle_start]
        for path_index in (cycle_start + 1):path_length
            representative = max(representative, path[path_index])
        end
        for path_index in cycle_start:path_length
            representatives[path[path_index]] = representative
        end
    end
    return representatives
end

function grouped_sparse_network_factorization_update!(
    workspace::GroupedSparseNetworkWorkspace,
    admittance::AbstractMatrix{<:Real},
    node_group_successors::AbstractVector{<:Integer};
    partition_boundary::Int=size(admittance, 1),
    first_factor_row::Int=2,
    pivot_tolerance::Real=0.0,
    zero_tolerance::Real=0.0,
)
    size(admittance, 1) == size(admittance, 2) ||
        throw(ArgumentError("admittance must be square"))
    node_count = size(admittance, 1)
    node_count == length(workspace.kk) || throw(ArgumentError(
        "grouped sparse workspace size must match admittance",
    ))
    2 <= first_factor_row <= partition_boundary <= node_count ||
        throw(ArgumentError("invalid first_factor_row or partition_boundary"))
    pivot_tol = Float64(pivot_tolerance)
    zero_tol = Float64(zero_tolerance)
    isfinite(pivot_tol) && pivot_tol >= 0.0 ||
        throw(ArgumentError("pivot_tolerance must be finite and nonnegative"))
    isfinite(zero_tol) && zero_tol >= 0.0 ||
        throw(ArgumentError("zero_tolerance must be finite and nonnegative"))

    representatives =
        _node_group_representatives!(workspace, node_group_successors)
    contracted = workspace.contracted_admittance
    fill!(contracted, 0.0)
    for column in 1:node_count
        representative_column = representatives[column]
        for row in 1:node_count
            value = Float64(admittance[row, column])
            isfinite(value) ||
                throw(ArgumentError("admittance entries must be finite"))
            contracted[representatives[row], representative_column] += value
        end
    end

    active_rows = workspace.active_rows
    resize!(active_rows, 0)
    for row in first_factor_row:partition_boundary
        representatives[row] == row && push!(active_rows, row)
    end
    km = workspace.km
    ykm = workspace.ykm
    kk = workspace.kk
    row_starts = workspace.row_starts
    pivot_values = workspace.pivot_values
    inverse_diagonal_values = workspace.inverse_diagonal_values
    factor_row = workspace.factor_row
    resize!(km, 0)
    resize!(ykm, 0)
    resize!(pivot_values, 0)
    resize!(inverse_diagonal_values, 0)
    fill!(kk, 0)
    fill!(row_starts, 0)

    for row in active_rows
        fill!(factor_row, 0.0)
        factor_row[1] = contracted[row, row]
        for column in first_factor_row:node_count
            column == row && continue
            representatives[column] == column || continue
            value = contracted[row, column]
            abs(value) > zero_tol && (factor_row[column] = value)
        end

        for previous_row in active_rows
            previous_row < row || break
            a = factor_row[previous_row]
            abs(a) > zero_tol || continue
            previous_start = row_starts[previous_row]
            previous_end = kk[previous_row]
            previous_start <= previous_end || throw(ArgumentError(
                "invalid previously built grouped sparse factor row",
            ))
            for index in (previous_start + 1):previous_end
                column = km[index]
                if column == row
                    factor_row[1] -= a * ykm[index]
                else
                    factor_row[column] -= a * ykm[index]
                end
            end
        end

        pivot = factor_row[1]
        abs(pivot) > pivot_tol || throw(ArgumentError(
            "grouped sparse factor pivot is zero or below tolerance",
        ))
        inverse_pivot = 1.0 / pivot
        push!(pivot_values, pivot)
        push!(inverse_diagonal_values, inverse_pivot)
        row_starts[row] = length(km) + 1
        push!(km, -row)
        push!(ykm, inverse_pivot)
        for column in (row + 1):node_count
            representatives[column] == column || continue
            value = factor_row[column]
            if abs(value) > zero_tol
                push!(km, column)
                push!(ykm, value * inverse_pivot)
            end
        end
        kk[row] = length(km)
    end
    workspace.factorization_count += 1
    return workspace
end

function grouped_sparse_network_solution!(
    workspace::GroupedSparseNetworkWorkspace,
    rhs::AbstractVector{<:Real},
    node_group_successors::AbstractVector{<:Integer},
    initial_solution::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(workspace.kk),
)
    first_factor_row >= 1 ||
        throw(ArgumentError("first_factor_row must be positive"))
    node_count = length(workspace.solution)
    length(initial_solution) == node_count || throw(ArgumentError(
        "initial solution length must match grouped sparse workspace",
    ))
    length(rhs) >= partition_boundary ||
        throw(ArgumentError("rhs must cover the active solved partition"))
    1 <= partition_boundary <= node_count ||
        throw(ArgumentError("partition_boundary must be within network nodes"))
    _check_node_group_successors(node_group_successors, node_count)

    solution = workspace.solution
    for index in 1:node_count
        value = Float64(initial_solution[index])
        isfinite(value) || throw(ArgumentError("solution entries must be finite"))
        solution[index] = value
    end
    for index in 1:partition_boundary
        value = Float64(rhs[index])
        isfinite(value) || throw(ArgumentError("rhs entries must be finite"))
        solution[index] = value
    end
    alias_successors = workspace.alias_successors
    for index in 1:node_count
        alias_successors[index] = Int(node_group_successors[index])
    end

    active_node = 1
    visited = 0
    while true
        solution[active_node] = 0.0
        alias_successors[active_node] <= active_node && break
        active_node = alias_successors[active_node]
        visited += 1
        visited <= node_count || throw(ArgumentError(
            "node group reference-chain search did not terminate",
        ))
    end

    for node in 2:partition_boundary
        successor = alias_successors[node]
        successor == 0 && continue
        successor > partition_boundary && continue
        successor > node && (solution[successor] += solution[node])
    end

    km = workspace.km
    ykm = workspace.ykm
    kk = workspace.kk
    entry_index = 1
    entry_count = length(km)
    while entry_index <= entry_count
        row = abs(km[entry_index])
        scale = solution[row]
        solution[row] = scale * ykm[entry_index]
        row_end = kk[row]
        while true
            entry_index += 1
            entry_index > row_end && break
            column = km[entry_index]
            column > partition_boundary && continue
            solution[column] -= scale * ykm[entry_index]
        end
    end

    node = node_count
    while true
        successor = alias_successors[node]
        if successor != 0 && successor <= node
            copied_from = node
            copied_to = successor
            visited = 0
            while true
                solution[copied_to] = solution[copied_from]
                copied_from = copied_to
                copied_to = alias_successors[copied_from]
                copied_to == node && break
                visited += 1
                visited <= node_count || throw(ArgumentError(
                    "node group copy cycle did not terminate",
                ))
            end
        end
        if node <= partition_boundary
            entry_index == 1 && break
            correction = 0.0
            while true
                entry_index -= 1
                column = km[entry_index]
                if column < 0
                    node = abs(column)
                    solution[node] += correction
                    break
                end
                correction -= solution[column] * ykm[entry_index]
            end
        else
            node -= 1
        end
    end
    return solution
end

function _grouped_sparse_network_solution(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real},
    node_group_successors::AbstractVector{<:Integer},
    initial_solution::AbstractVector{<:Real};
    first_factor_row::Int=2,
    partition_boundary::Int=length(kk),
)
    first_factor_row >= 1 ||
        throw(ArgumentError("first_factor_row must be positive"))
    length(km) == length(ykm) ||
        throw(ArgumentError("km and ykm lengths must match"))
    node_count = length(initial_solution)
    length(kk) >= node_count ||
        throw(ArgumentError("kk must cover every network node"))
    length(rhs) >= partition_boundary ||
        throw(ArgumentError("rhs must cover the active solved partition"))
    1 <= partition_boundary <= node_count ||
        throw(ArgumentError("partition_boundary must be within network nodes"))
    _check_node_group_successors(node_group_successors, node_count)
    for entry in eachindex(km, ykm)
        node = abs(km[entry])
        1 <= node <= node_count ||
            throw(ArgumentError("km entries must point to a valid network node"))
        isfinite(Float64(ykm[entry])) ||
            throw(ArgumentError("ykm entries must be finite"))
    end

    solution = Float64.(initial_solution)
    solution[1:partition_boundary] .= Float64.(rhs[1:partition_boundary])
    for value in solution
        isfinite(value) || throw(ArgumentError("solution entries must be finite"))
    end

    alias_successors = Int.(node_group_successors)
    active_node = 1
    visited = 0
    while true
        solution[active_node] = 0.0
        alias_successors[active_node] <= active_node && break
        active_node = alias_successors[active_node]
        visited += 1
        visited <= node_count ||
            throw(ArgumentError("node group reference-chain search did not terminate"))
    end

    for node in 2:partition_boundary
        successor = alias_successors[node]
        successor == 0 && continue
        successor > partition_boundary && continue
        successor > node && (solution[successor] += solution[node])
    end

    entry_index = 1
    entry_count = length(km)
    while entry_index <= entry_count
        row = abs(km[entry_index])
        scale = solution[row]
        solution[row] = scale * Float64(ykm[entry_index])
        row_end = kk[row]
        while true
            entry_index += 1
            entry_index > row_end && break
            column = km[entry_index]
            column > partition_boundary && continue
            solution[column] -= scale * Float64(ykm[entry_index])
        end
    end

    node = node_count
    while true
        successor = alias_successors[node]
        if successor != 0 && successor <= node
            copied_from = node
            copied_to = successor
            visited = 0
            while true
                solution[copied_to] = solution[copied_from]
                copied_from = copied_to
                copied_to = alias_successors[copied_from]
                copied_to == node && break
                visited += 1
                visited <= node_count ||
                    throw(ArgumentError("node group copy cycle did not terminate"))
            end
        end
        if node <= partition_boundary
            entry_index == 1 && break
            correction = 0.0
            while true
                entry_index -= 1
                column = km[entry_index]
                if column < 0
                    node = abs(column)
                    solution[node] += correction
                    break
                end
                correction -= solution[column] * Float64(ykm[entry_index])
            end
        else
            node -= 1
        end
    end

    return solution
end

function over16_network_solution_update!(
    state::OVER16SwitchCurrentState,
    factor::OVER16FortranSparseFactorWorkspaceState;
    kwargs...,
)
    before = copy(state.network_solution)
    solution = sparse_network_solution(factor, state.rhs; kwargs...)
    resize!(state.network_solution, length(solution))
    state.network_solution .= solution
    network_solution_mutated = state.network_solution != before
    return (
        rhs = copy(state.rhs),
        network_solution = copy(state.network_solution),
        network_solution_mutated = network_solution_mutated,
        switch_current_state_mutated = network_solution_mutated,
    )
end

function sparse_network_factorization_update(admittance::AbstractMatrix{<:Real}; kwargs...)
    return over16_fortran_sparse_factor_update(admittance; kwargs...)
end

function sparse_network_factorization_update!(
    state::OVER16FortranSparseFactorWorkspaceState,
    admittance::AbstractMatrix{<:Real};
    kwargs...,
)
    return over16_fortran_sparse_factor_update!(state, admittance; kwargs...)
end

function sparse_network_solution(
    km::AbstractVector{Int},
    ykm::AbstractVector{<:Real},
    kk::AbstractVector{Int},
    rhs::AbstractVector{<:Real};
    node_group_successors::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    initial_solution::Union{Nothing,AbstractVector{<:Real}}=nothing,
    kwargs...,
)
    if node_group_successors !== nothing
        seed = initial_solution === nothing ? zeros(Float64, length(kk)) : initial_solution
        return _grouped_sparse_network_solution(
            km,
            ykm,
            kk,
            rhs,
            node_group_successors,
            seed;
            kwargs...,
        )
    end
    return over16_fortran_sparse_factor_solve(km, ykm, kk, rhs; kwargs...)
end

function sparse_network_solution(
    state::OVER16FortranSparseFactorWorkspaceState,
    rhs::AbstractVector{<:Real};
    kwargs...,
)
    return sparse_network_solution(state.km, state.ykm, state.kk, rhs; kwargs...)
end

function sparse_network_solution_update!(
    state::OVER16SwitchCurrentState,
    factor::OVER16FortranSparseFactorWorkspaceState;
    kwargs...,
)
    return over16_network_solution_update!(state, factor; kwargs...)
end

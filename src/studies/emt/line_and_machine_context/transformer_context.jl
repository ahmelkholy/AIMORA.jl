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

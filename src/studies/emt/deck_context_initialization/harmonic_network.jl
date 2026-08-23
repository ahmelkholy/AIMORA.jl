function _stamp_complex_branch_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    from_node::Integer,
    to_node::Integer,
    admittance::ComplexF64,
)
    a = Int(from_node)
    b = Int(to_node)
    if a != 0
        matrix[a, a] += admittance
    end
    if b != 0
        matrix[b, b] += admittance
    end
    if a != 0 && b != 0
        matrix[a, b] -= admittance
        matrix[b, a] -= admittance
    end
    return matrix
end

function _stamp_complex_phase_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    admittance::AbstractMatrix{ComplexF64},
)
    phase_count = length(from_nodes)
    length(to_nodes) == phase_count ||
        throw(ArgumentError("phase admittance terminal counts must match"))
    size(admittance) == (phase_count, phase_count) ||
        throw(ArgumentError("phase admittance matrix size must match terminal count"))
    for row in 1:phase_count, column in 1:phase_count
        value = admittance[row, column]
        from_row = Int(from_nodes[row])
        to_row = Int(to_nodes[row])
        from_column = Int(from_nodes[column])
        to_column = Int(to_nodes[column])
        if from_row != 0 && from_column != 0
            matrix[from_row, from_column] += value
        end
        if from_row != 0 && to_column != 0
            matrix[from_row, to_column] -= value
        end
        if to_row != 0 && from_column != 0
            matrix[to_row, from_column] -= value
        end
        if to_row != 0 && to_column != 0
            matrix[to_row, to_column] += value
        end
    end
    return matrix
end

function _deck_steady_state_frequency_hz(parsed::DeckParser.DeckParseResult)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    options.x_frequency_hz > 0.0 && return options.x_frequency_hz
    for row in DeckParser.deck_over5a_source_rows(parsed)
        row.sfreq > 0.0 && return Float64(row.sfreq) / (2.0 * pi)
    end
    return 60.0
end

function _deck_current_source_seed_clear_nodes(parsed::DeckParser.DeckParseResult)
    return unique(
        abs(Int(row.node_value))
        for row in DeckParser.deck_over5a_source_rows(parsed)
    )
end
function _steady_state_terminal_frequency_hz(
    partition::DeckParser.DeckSteadyStateFrequencyPartition,
    node_indices,
    default_frequency_hz::Float64,
)
    frequency_hz = nothing
    for node in node_indices
        node_index = abs(Int(node))
        node_index == 0 && continue
        1 <= node_index <= length(partition.node_frequencies_hz) ||
            throw(ArgumentError("steady-state terminal is outside the frequency partition"))
        node_frequency_hz = partition.node_frequencies_hz[node_index]
        node_frequency_hz == 0.0 && continue
        if frequency_hz === nothing
            frequency_hz = node_frequency_hz
        elseif !isapprox(
            frequency_hz,
            node_frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        )
            throw(ArgumentError("coupled steady-state terminals have different frequencies"))
        end
    end
    return frequency_hz === nothing ? default_frequency_hz : frequency_hz
end

function _deck_source_voltage_phasor(row)
    if row.iform == 11
        return complex(Float64(row.crest), 0.0)
    elseif row.iform == 14
        return Float64(row.crest) * cis(Float64(row.time1))
    end
    throw(ArgumentError("steady-state source phasors require constant or sinusoidal source rows"))
end

function _is_basic_harmonic_element(element)
    return element isa Union{
        ConductanceBranch,
        SeriesRLBranch,
        SeriesRLCBranch,
        CapacitorBranch,
        TheveninSource,
        CurrentInjection,
    }
end

function _claim_basic_harmonic_element_owner!(
    owners::Vector{Symbol},
    parsed::DeckParser.DeckParseResult,
    owner::Symbol,
    name::Symbol,
    line_no::Int;
    required::Bool,
)
    candidates = Int[
        index for index in eachindex(parsed.elements)
        if owners[index] === :named_basic_element &&
           parsed.element_names[index] === name &&
           parsed.element_line_numbers[index] == line_no
    ]
    if isempty(candidates)
        required && throw(ArgumentError(
            "harmonic owner $owner has no matching basic element for $name on line $line_no",
        ))
        return owners
    end
    length(candidates) == 1 || throw(ArgumentError(
        "harmonic owner $owner is ambiguous for $name on line $line_no",
    ))
    owners[only(candidates)] = owner
    return owners
end

function _basic_harmonic_element_owners(parsed::DeckParser.DeckParseResult)
    owners = Symbol[
        _is_basic_harmonic_element(element) ?
        :named_basic_element : :not_a_basic_harmonic_element
        for element in parsed.elements
    ]
    for row in DeckParser.deck_over2_branch_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :fixed_card_scalar_branch,
            row.name,
            row.line_no;
            required=true,
        )
    end
    for row in DeckParser.deck_over5a_source_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :fixed_card_source,
            row.name,
            row.line_no;
            required=false,
        )
    end
    for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        row.reactance === missing && continue
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :universal_machine_terminal_reactance,
            Symbol(
                "machine_terminal_reactance_",
                row.machine_index,
                "_",
                row.branch_index,
            ),
            row.line_no;
            required=true,
        )
    end
    for row in DeckParser.deck_universal_machine_speed_capacitor_rows(parsed)
        _claim_basic_harmonic_element_owner!(
            owners,
            parsed,
            :universal_machine_speed_branch,
            Symbol("machine_speed_capacitor_", row.machine_index),
            row.line_no;
            required=true,
        )
    end
    return owners
end

function _stamp_named_basic_harmonic_elements!(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    physical_frequency_hz::Float64,
    formulation::AbstractEMTHarmonicFormulation,
)
    owners = _basic_harmonic_element_owners(parsed)
    physical_angular_frequency = 2.0 * pi * physical_frequency_hz
    reactive_angular_frequency = _emt_reactive_angular_frequency(
        formulation,
        physical_angular_frequency,
    )
    for (element_index, element) in pairs(parsed.elements)
        owners[element_index] === :named_basic_element || continue
        if element isa ConductanceBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(element.g, 0.0),
            )
        elseif element isa SeriesRLBranch
            impedance = complex(
                element.r,
                reactive_angular_frequency * element.l,
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                inv(impedance),
            )
        elseif element isa SeriesRLCBranch
            physical_angular_frequency > 0.0 || throw(ArgumentError(
                "series R-L-C harmonic initialization requires positive frequency",
            ))
            impedance = complex(
                element.r,
                reactive_angular_frequency * element.l -
                    inv(reactive_angular_frequency * element.c),
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                inv(impedance),
            )
        elseif element isa CapacitorBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(0.0, reactive_angular_frequency * element.c),
            )
        elseif element isa TheveninSource
            element.value isa SinusoidalSourceSignal || throw(ArgumentError(
                "named Thevenin harmonic initialization requires a typed sinusoidal source",
            ))
            source_phasor = sinusoidal_source_peak_phasor(
                element.value,
                physical_frequency_hz,
            )
            _stamp_complex_branch_admittance!(
                admittance,
                element.node,
                0,
                complex(element.g, 0.0),
            )
            rhs[element.node] += element.g * source_phasor
        elseif element isa CurrentInjection
            element.value isa SinusoidalSourceSignal || throw(ArgumentError(
                "named harmonic current injection requires a typed sinusoidal source",
            ))
            rhs[element.node] += sinusoidal_source_peak_phasor(
                element.value,
                physical_frequency_hz,
            )
        end
    end
    return admittance
end

function _stamp_deck_branch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_over2_branch_rows(parsed)
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_value,
            row.to_node_value,
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        if row.branch_kind == :conductance
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                complex(Float64(row.conductance), 0.0),
            )
        elseif row.branch_kind == :series_rl
            impedance = complex(Float64(row.resistance), angular_frequency * Float64(row.inductance))
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                inv(impedance),
            )
        elseif row.branch_kind == :capacitor &&
               Float64(row.raw_capacitance) != Float64(row.capacitance)
            _stamp_complex_branch_admittance!(
                matrix,
                row.from_node_value,
                row.to_node_value,
                complex(0.0, angular_frequency * Float64(row.capacitance)),
            )
        end
    end
    for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        row.machine_index in excluded_universal_machine_indices && continue
        row.reactance === missing && continue
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_value,
            row.to_node_value,
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        inductance = DeckParser.fixed_card_branch_timestep_inductance(
            parsed,
            Float64(row.reactance),
        )
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            inv(complex(0.0, angular_frequency * inductance)),
        )
    end
    for row in DeckParser.deck_universal_machine_speed_capacitor_rows(parsed)
        admittance = if row.resistance > 0.0 && row.capacitance == 0.0
            complex(inv(Float64(row.resistance)), 0.0)
        elseif row.resistance == 0.0 && row.capacitance > 0.0
            frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
                frequency_partition,
                row.capacitor_node_value,
                row.mass_node_value,
                default_frequency_hz,
            )
            angular_frequency = _emt_reactive_angular_frequency(
                formulation,
                2.0 * pi * frequency_hz,
            )
            complex(0.0, angular_frequency * Float64(row.capacitance))
        else
            throw(ArgumentError("unsupported universal-machine speed-capacitor impedance"))
        end
        _stamp_complex_branch_admittance!(
            matrix,
            row.capacitor_node_value,
            row.mass_node_value,
            admittance,
        )
    end
    for shunt in deck_switching_nonlinear_resistor_safety_shunts(parsed)
        _stamp_complex_branch_admittance!(
            matrix,
            shunt.from_node_index,
            shunt.to_node_index,
            complex(shunt.conductance_s, 0.0),
        )
    end
    return matrix
end

function _stamp_hysteretic_inductor_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_hysteretic_inductor_rows(parsed)
        current_a = Float64(row.steady_state_current_A)
        flux_wb = Float64(row.steady_state_flux_Wb)
        current_a == 0.0 && continue
        current_a > 0.0 && flux_wb > 0.0 || throw(ArgumentError(
            "a nonzero hysteretic-inductor steady-state current requires a positive steady-state flux",
        ))
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_index,
            row.to_node_index,
            default_frequency_hz,
        )
        reactive_angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        incremental_inductance_h = flux_wb / current_a
        impedance = complex(0.0, reactive_angular_frequency * incremental_inductance_h)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_index,
            row.to_node_index,
            inv(impedance),
        )
    end
    return matrix
end

function _stamp_piecewise_nonlinear_inductor_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_piecewise_nonlinear_inductor_rows(parsed)
        row.nonlinear_type == PIECEWISE_NONLINEAR_INDUCTOR_TYPE || continue
        current_a = Float64(row.steady_state_current_a)
        flux_wb = Float64(row.steady_state_flux_wb)
        if current_a == 0.0 && flux_wb == 0.0
            continue
        end
        current_a != 0.0 && flux_wb != 0.0 && flux_wb / current_a > 0.0 ||
            throw(ArgumentError(
                "a nonzero piecewise nonlinear-inductor steady state requires finite current and flux with a positive secant inductance",
            ))
        frequency_hz = DeckParser.deck_node_steady_state_frequency_hz(
            frequency_partition,
            row.from_node_index,
            row.to_node_index,
            default_frequency_hz,
        )
        reactive_angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        reactive_angular_frequency > 0.0 || throw(ArgumentError(
            "piecewise nonlinear-inductor harmonic initialization requires positive frequency",
        ))
        secant_inductance_h = flux_wb / current_a
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_index,
            row.to_node_index,
            inv(complex(0.0, reactive_angular_frequency * secant_inductance_h)),
        )
    end
    return matrix
end

function _stamp_deck_induction_machine_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    definitions = DeckParser.deck_universal_machine_definition_rows(parsed)
    terminals = DeckParser.deck_universal_machine_terminal_rows(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for machine_index in sort(unique(row.machine_index for row in definitions))
        card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
        card1.machine_type in (3, 4, 5, 6, 7) || continue
        _deck_universal_machine_initialization_mode(parsed) == :manual &&
            continue
        machine_terminals = sort!(
            [row for row in terminals if row.machine_index == machine_index];
            by = row -> row.terminal_index,
        )
        active_terminals = card1.machine_type in (6, 7) ?
            machine_terminals[3:3] :
            card1.machine_type == 5 ? machine_terminals[2:3] : machine_terminals[1:3]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(
                [row.terminal_node_value for row in active_terminals],
                [row.reference_node_value for row in active_terminals],
            ),
            default_frequency_hz,
        )
        equivalent = _deck_induction_machine_steady_state_equivalent(
            parsed,
            machine_index;
            frequency_hz = frequency_hz,
        )
        for row in active_terminals
            _stamp_complex_branch_admittance!(
                matrix,
                row.terminal_node_value,
                row.reference_node_value,
                equivalent.terminal_admittance,
            )
        end
    end
    return matrix
end

function _deck_induction_machine_steady_state_equivalent(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int;
    frequency_hz::Real=_deck_steady_state_frequency_hz(parsed),
)
    card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    card2.value2 === missing &&
        throw(ArgumentError("induction-machine main inductance is missing"))
    machine_coils = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    coil_count = _deck_coupled_dq_coil_count(card1)
    length(machine_coils) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) induction-machine steady state requires $coil_count coils"))
    power_resistances = filter(>(0.0), Float64[row.resistance for row in machine_coils[1:3]])
    rotor_resistances = filter(>(0.0), Float64[row.resistance for row in machine_coils[4:coil_count]])
    rotor_inductances = filter(>(0.0), Float64[row.inductance for row in machine_coils[4:coil_count]])
    isempty(power_resistances) && throw(ArgumentError("power-coil resistance is missing"))
    isempty(rotor_resistances) && throw(ArgumentError("rotor resistance is missing"))
    isempty(rotor_inductances) && throw(ArgumentError("rotor leakage inductance is missing"))
    section = _deck_universal_machine_section(parsed)
    angular_frequency_rad_s = 2.0 * pi * Float64(frequency_hz)
    inductance_scale =
        section.parameter_basis == :power_frequency_normalized ?
        inv(angular_frequency_rad_s) : 1.0
    slip = (card4.value1 === missing ? 0.0 : Float64(card4.value1)) / 100.0
    return induction_machine_steady_state_equivalent(
        power_coil_resistance = first(power_resistances),
        rotor_resistance = first(rotor_resistances),
        rotor_leakage_inductance = first(rotor_inductances) * inductance_scale,
        main_inductance = Float64(card2.value2) * inductance_scale,
        slip = slip,
        frequency_hz = frequency_hz,
    )
end

function _stamp_deck_switch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        conductance = row.initially_closed ?
            Float64(row.on_conductance) :
            Float64(row.off_conductance)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _stamp_deck_open_switch_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        if row.initially_closed
            if row.from_node_value != 0 && row.to_node_value != 0
                continue
            end
            conductance = Float64(row.on_conductance)
        else
            conductance = Float64(row.off_conductance)
        end
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _stamp_deck_switch_at_time_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    time_s::Float64,
)
    for row in DeckParser.deck_over5_switch_rows(parsed)
        closed = _deck_time_switch_closed_at(
            row.initially_closed,
            Float64(row.close_time_s),
            Float64(row.open_time_s),
            time_s,
        )
        closed && row.from_node_value != 0 && row.to_node_value != 0 && continue
        conductance = closed ? Float64(row.on_conductance) :
            Float64(row.off_conductance)
        _stamp_complex_branch_admittance!(
            matrix,
            row.from_node_value,
            row.to_node_value,
            complex(conductance, 0.0),
        )
    end
    return matrix
end

function _steady_state_closed_switch_representatives(
    parsed::DeckParser.DeckParseResult,
    node_count::Int,
    time_s::Union{Nothing,Float64}=nothing,
    ;
    additional_closed_node_pairs::AbstractVector{<:Tuple}=Tuple{Int,Int}[],
    include_grounded_constraints::Bool=false,
)
    representatives = collect(1:node_count)
    grounded_nodes = Int[]

    function representative(node::Int)
        root = node
        while representatives[root] != root
            root = representatives[root]
        end
        while representatives[node] != node
            parent = representatives[node]
            representatives[node] = root
            node = parent
        end
        return root
    end

    function union_nodes(left::Int, right::Int)
        left_root = representative(left)
        right_root = representative(right)
        left_root == right_root && return nothing
        retained = max(left_root, right_root)
        replaced = min(left_root, right_root)
        representatives[replaced] = retained
        return nothing
    end

    function connect_nodes(left::Int, right::Int)
        0 <= left <= node_count && 0 <= right <= node_count ||
            throw(ArgumentError("steady-state closed-switch node is outside the network"))
        left != right || return nothing
        if left == 0 || right == 0
            include_grounded_constraints && push!(grounded_nodes, max(left, right))
            return nothing
        end
        return union_nodes(left, right)
    end

    for row in DeckParser.deck_over5_switch_rows(parsed)
        closed = time_s === nothing ? row.initially_closed :
            _deck_time_switch_closed_at(
                row.initially_closed,
                Float64(row.close_time_s),
                Float64(row.open_time_s),
                time_s,
            )
        closed || continue
        from_node = Int(row.from_node_value)
        to_node = Int(row.to_node_value)
        connect_nodes(from_node, to_node)
    end
    for row in DeckParser.deck_control_system_switch_coupling_rows(parsed)
        _control_system_switch_initially_closed(row.initial_state) || continue
        row.from_node_index !== missing && row.to_node_index !== missing ||
            throw(ArgumentError(
                "steady-state controlled switch on line $(row.line_no) has unresolved nodes",
            ))
        from_node = Int(row.from_node_index)
        to_node = Int(row.to_node_index)
        connect_nodes(from_node, to_node)
    end
    for pair in additional_closed_node_pairs
        length(pair) == 2 ||
            throw(ArgumentError("additional steady-state closed-node pair must contain two nodes"))
        from_node = Int(pair[1])
        to_node = Int(pair[2])
        connect_nodes(from_node, to_node)
    end

    grounded_representatives = Set(
        representative(node) for node in grounded_nodes
    )
    for node in 1:node_count
        root = representative(node)
        representatives[node] = root in grounded_representatives ? 0 : root
    end
    return representatives
end

function _group_harmonic_network(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    node_count = length(rhs)
    size(admittance, 1) == node_count && size(admittance, 2) == node_count ||
        throw(ArgumentError("steady-state admittance and RHS dimensions must match"))
    length(representatives) == node_count ||
        throw(ArgumentError("steady-state switch representatives must cover every node"))
    all(representative -> 0 <= representative <= node_count, representatives) ||
        throw(ArgumentError("steady-state switch representative is outside the network"))
    identity_representatives = all(
        node -> Int(representatives[node]) == node,
        1:node_count,
    )
    if identity_representatives
        active_nodes = collect(1:node_count)
        reduced_index = Dict(node => node for node in active_nodes)
        active_node_groups = Vector{Int}[Int[node] for node in active_nodes]
        return (;
            active_nodes,
            reduced_index,
            active_node_groups,
            grounded_node_indices=Int[],
            switch_node_groups=active_node_groups,
            reduced_admittance=admittance,
            reduced_rhs=rhs,
        )
    end
    active_nodes = sort(filter(!iszero, unique(Int.(representatives))))
    reduced_index = Dict(node => index for (index, node) in enumerate(active_nodes))
    active_node_groups = Vector{Int}[
        findall(representative -> Int(representative) == active_node, representatives)
        for active_node in active_nodes
    ]
    grounded_node_indices = findall(iszero, representatives)
    switch_node_groups = isempty(grounded_node_indices) ?
        active_node_groups : vcat(active_node_groups, [grounded_node_indices])
    reduced_admittance = zeros(ComplexF64, length(active_nodes), length(active_nodes))
    reduced_rhs = zeros(ComplexF64, length(active_nodes))
    for row in 1:node_count
        representatives[row] == 0 && continue
        reduced_row = reduced_index[Int(representatives[row])]
        reduced_rhs[reduced_row] += rhs[row]
        for column in 1:node_count
            representatives[column] == 0 && continue
            reduced_column = reduced_index[Int(representatives[column])]
            reduced_admittance[reduced_row, reduced_column] +=
                admittance[row, column]
        end
    end
    return (;
        active_nodes,
        reduced_index,
        active_node_groups,
        grounded_node_indices,
        switch_node_groups,
        reduced_admittance,
        reduced_rhs,
    )
end

function _expand_harmonic_component(
    component::AbstractVector{<:Integer},
    active_node_groups::Vector{Vector{Int}},
)
    expanded_count = sum(
        length(active_node_groups[Int(group_index)])
        for group_index in component
    )
    expanded = Vector{Int}(undef, expanded_count)
    destination_index = 1
    for group_index in component
        group = active_node_groups[Int(group_index)]
        copyto!(expanded, destination_index, group, 1, length(group))
        destination_index += length(group)
    end
    return sort!(expanded)
end

function _solve_grouped_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
    ;
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
    passive_conductance_network::Bool=false,
)
    node_count = length(rhs)
    grouped = _group_harmonic_network(admittance, rhs, representatives)
    active_nodes = grouped.active_nodes
    reduced_index = grouped.reduced_index
    active_node_groups = grouped.active_node_groups

    reduced = isempty(active_nodes) ?
        (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            connected_components=Vector{Vector{Int}}(),
            referenced_components=BitVector(),
            unreferenced_components=Vector{Vector{Int}}(),
            classification=:unique,
        ) :
        passive_conductance_network ?
        _solve_passive_conductance_harmonic_linear_system(
            grouped.reduced_admittance,
            grouped.reduced_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        ) :
        _solve_harmonic_linear_system(
            grouped.reduced_admittance,
            grouped.reduced_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    solution = if reduced.solution === nothing
        nothing
    else
        expanded = Vector{ComplexF64}(undef, node_count)
        for node in 1:node_count
            representative = Int(representatives[node])
            expanded[node] = representative == 0 ? 0.0 + 0.0im :
                reduced.solution[reduced_index[representative]]
        end
        expanded
    end
    return (
        solution,
        numerical_rank=reduced.numerical_rank,
        condition_estimate=reduced.condition_estimate,
        maximum_residual_a=reduced.maximum_residual_a,
        relative_residual=reduced.relative_residual,
        connected_components=Vector{Int}[
            _expand_harmonic_component(component, active_node_groups)
            for component in reduced.connected_components
        ],
        referenced_components=copy(reduced.referenced_components),
        unreferenced_components=Vector{Int}[
            _expand_harmonic_component(component, active_node_groups)
            for component in reduced.unreferenced_components
        ],
        classification=reduced.classification,
        reduced_node_count=length(active_nodes),
        switch_node_groups=grouped.switch_node_groups,
    )
end

function _solve_grouped_constrained_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
    fixed_node_phasors::AbstractDict{<:Integer,<:Complex};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
    passive_conductance_network::Bool=false,
)
    isempty(fixed_node_phasors) && return _solve_grouped_harmonic_linear_system(
        admittance,
        rhs,
        representatives;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
        passive_conductance_network,
    )
    grouped = _group_harmonic_network(admittance, rhs, representatives)
    node_count = length(rhs)
    reduced_node_count = length(grouped.active_nodes)
    fixed_values = Dict{Int,ComplexF64}()
    for (node_value, phasor_value) in fixed_node_phasors
        node = Int(node_value)
        1 <= node <= node_count || throw(ArgumentError(
            "fixed harmonic node is outside the network",
        ))
        phasor = ComplexF64(phasor_value)
        all(isfinite, (real(phasor), imag(phasor))) || throw(ArgumentError(
            "fixed harmonic phasor must be finite",
        ))
        representative = Int(representatives[node])
        if representative == 0
            iszero(phasor) || throw(ArgumentError(
                "grounded closed-switch node group requires a zero fixed harmonic phasor",
            ))
            continue
        end
        reduced_node = grouped.reduced_index[representative]
        haskey(fixed_values, reduced_node) && fixed_values[reduced_node] != phasor &&
            throw(ArgumentError(
                "closed-switch node group has inconsistent fixed harmonic phasors",
            ))
        fixed_values[reduced_node] = phasor
    end
    fixed_indices = sort!(collect(keys(fixed_values)))
    unknown_indices = setdiff(collect(1:reduced_node_count), fixed_indices)
    reduced_solution = zeros(ComplexF64, reduced_node_count)
    reduced_solution[fixed_indices] .= ComplexF64[
        fixed_values[index] for index in fixed_indices
    ]
    solved = if isempty(unknown_indices)
        (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            classification=:unique,
        )
    else
        constrained_rhs = grouped.reduced_rhs[unknown_indices] -
            grouped.reduced_admittance[unknown_indices, fixed_indices] *
            reduced_solution[fixed_indices]
        passive_conductance_network ?
        _solve_passive_conductance_harmonic_linear_system(
            grouped.reduced_admittance[unknown_indices, unknown_indices],
            constrained_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        ) :
        _solve_harmonic_linear_system(
            grouped.reduced_admittance[unknown_indices, unknown_indices],
            constrained_rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    end
    solved.solution === nothing ||
        (reduced_solution[unknown_indices] .= solved.solution)

    matrix_scale = maximum(abs, grouped.reduced_admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * reduced_node_count * eps(Float64) *
        Float64(rank_relative_threshold_multiplier),
        floatmin(Float64),
    )
    reduced_components = _harmonic_connected_components(
        grouped.reduced_admittance,
        structural_threshold,
    )
    natural_references = _harmonic_component_references(
        grouped.reduced_admittance,
        reduced_components,
        structural_threshold,
    )
    fixed_index_set = Set(fixed_indices)
    referenced_components = BitVector(
        natural_reference || any(node -> node in fixed_index_set, component)
        for (component, natural_reference) in
            zip(reduced_components, natural_references)
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in
            zip(reduced_components, referenced_components)
        if !referenced
    ]
    classification = if solved.classification === :infeasible
        :infeasible
    elseif solved.numerical_rank < length(unknown_indices)
        isempty(unreferenced_components) ? :nonunique : :islanded
    else
        solved.classification
    end
    expanded_solution = classification === :unique ? ComplexF64[
        Int(representatives[node]) == 0 ? 0.0 + 0.0im :
        reduced_solution[grouped.reduced_index[Int(representatives[node])]]
        for node in 1:node_count
    ] : nothing
    reaction_currents = expanded_solution === nothing ? ComplexF64[] :
        admittance * expanded_solution - rhs
    return (
        solution=expanded_solution,
        numerical_rank=solved.numerical_rank + length(fixed_indices),
        condition_estimate=solved.condition_estimate,
        maximum_residual_a=solved.maximum_residual_a,
        relative_residual=solved.relative_residual,
        connected_components=Vector{Int}[
            _expand_harmonic_component(
                component,
                grouped.active_node_groups,
            ) for component in reduced_components
        ],
        referenced_components,
        unreferenced_components=Vector{Int}[
            _expand_harmonic_component(
                component,
                grouped.active_node_groups,
            ) for component in unreferenced_components
        ],
        classification,
        reduced_node_count,
        switch_node_groups=grouped.switch_node_groups,
        constraint_reaction_current_phasors=reaction_currents,
    )
end

function _grouped_steady_state_diagnostics(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    result = _solve_grouped_harmonic_linear_system(
        admittance,
        rhs,
        representatives;
        current_absolute_a=1.0e-12,
        current_relative=1.0e-10,
        rank_relative_threshold_multiplier=10.0,
        maximum_condition_estimate=Inf,
    )
    result.classification === :unique || throw(ArgumentError(
        "grouped steady-state network classification $(result.classification): " *
        "rank $(result.numerical_rank)/$(result.reduced_node_count), " *
        "condition $(result.condition_estimate), residual $(result.maximum_residual_a) A",
    ))
    return result
end

function _solve_grouped_steady_state_admittance(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    representatives::AbstractVector{<:Integer},
)
    return something(_grouped_steady_state_diagnostics(
        admittance,
        rhs,
        representatives,
    ).solution)
end

function _solve_steady_state_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
)
    result = _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a=1.0e-12,
        current_relative=1.0e-10,
        rank_relative_threshold_multiplier=10.0,
        maximum_condition_estimate=Inf,
    )
    result.classification === :unique || throw(ArgumentError(
        "steady-state network classification $(result.classification): " *
        "rank $(result.numerical_rank)/$(size(admittance, 1)), " *
        "condition $(result.condition_estimate), residual $(result.maximum_residual_a) A",
    ))
    return something(result.solution)
end

function _harmonic_connected_components(
    admittance::AbstractMatrix{ComplexF64},
    threshold::Float64,
)
    node_count = size(admittance, 1)
    visited = falses(node_count)
    components = Vector{Vector{Int}}()
    for initial_node in 1:node_count
        visited[initial_node] && continue
        component = Int[]
        stack = Int[initial_node]
        visited[initial_node] = true
        while !isempty(stack)
            node = pop!(stack)
            push!(component, node)
            for neighbor in 1:node_count
                neighbor == node && continue
                visited[neighbor] && continue
                if abs(admittance[node, neighbor]) > threshold ||
                   abs(admittance[neighbor, node]) > threshold
                    visited[neighbor] = true
                    push!(stack, neighbor)
                end
            end
        end
        sort!(component)
        push!(components, component)
    end
    return components
end

function _harmonic_component_references(
    admittance::AbstractMatrix{ComplexF64},
    components::Vector{Vector{Int}},
    threshold::Float64,
)
    references = falses(length(components))
    for (component_index, component) in enumerate(components)
        for node in component
            row_sum = sum(admittance[node, column] for column in axes(admittance, 2))
            if abs(row_sum) > threshold
                references[component_index] = true
                break
            end
        end
    end
    return references
end

function _solve_passive_conductance_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    size(admittance) == (node_count, node_count) || throw(DimensionMismatch(
        "passive harmonic admittance and right-hand side dimensions must match",
    ))
    node_count == 0 && return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    ishermitian(admittance) || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    if all(value -> iszero(imag(value)), admittance)
        return _solve_real_passive_conductance_harmonic_linear_system(
            admittance,
            rhs;
            current_absolute_a,
            current_relative,
            rank_relative_threshold_multiplier,
            maximum_condition_estimate,
        )
    end
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    matrix_scale = maximum(abs, admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(admittance, node, :))
        column_magnitude = sum(abs, view(admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_admittance = Matrix{ComplexF64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        scaled_admittance[row, column] = admittance[row, column] /
            (symmetric_scales[row] * symmetric_scales[column])
    end
    scaled_eigenvalues = eigvals!(Hermitian(scaled_admittance))
    leading_eigenvalue = maximum(scaled_eigenvalues; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_eigenvalue * rank_multiplier
    minimum_eigenvalue = minimum(scaled_eigenvalues; init=0.0)
    minimum_eigenvalue > rank_threshold || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    condition_estimate = leading_eigenvalue / minimum_eigenvalue
    physical_admittance = Matrix{ComplexF64}(admittance)
    physical_factor = cholesky!(Hermitian(physical_admittance))
    candidate_solution = physical_factor \ rhs
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    components = _harmonic_connected_components(
        admittance,
        structural_threshold,
    )
    referenced_components = _harmonic_component_references(
        admittance,
        components,
        structural_threshold,
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    return (
        solution=classification === :unique ? candidate_solution : nothing,
        numerical_rank=node_count,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

function _solve_real_passive_conductance_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    real_admittance = Matrix{Float64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        real_admittance[row, column] = real(admittance[row, column])
    end
    tridiagonal_entries =
        _real_symmetric_tridiagonal_entries(real_admittance)
    matrix_scale = maximum(abs, real_admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(real_admittance, node, :))
        column_magnitude = sum(abs, view(real_admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_eigenvalues = if tridiagonal_entries === nothing
        scaled_admittance = Matrix{Float64}(undef, node_count, node_count)
        for column in 1:node_count, row in 1:node_count
            scaled_admittance[row, column] = real_admittance[row, column] /
                (symmetric_scales[row] * symmetric_scales[column])
        end
        eigvals!(Symmetric(scaled_admittance))
    else
        scaled_diagonal = Vector{Float64}(undef, node_count)
        for node in 1:node_count
            scaled_diagonal[node] = tridiagonal_entries.diagonal[node] /
                (symmetric_scales[node] * symmetric_scales[node])
        end
        scaled_off_diagonal = Vector{Float64}(undef, node_count - 1)
        for node in 1:(node_count - 1)
            scaled_off_diagonal[node] =
                tridiagonal_entries.off_diagonal[node] /
                (symmetric_scales[node] * symmetric_scales[node + 1])
        end
        eigvals!(SymTridiagonal(scaled_diagonal, scaled_off_diagonal))
    end
    leading_eigenvalue = maximum(scaled_eigenvalues; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_eigenvalue * rank_multiplier
    minimum_eigenvalue = minimum(scaled_eigenvalues; init=0.0)
    minimum_eigenvalue > rank_threshold || return _solve_harmonic_linear_system(
        admittance,
        rhs;
        current_absolute_a,
        current_relative,
        rank_relative_threshold_multiplier,
        maximum_condition_estimate,
    )
    condition_estimate = leading_eigenvalue / minimum_eigenvalue
    physical_factor = cholesky!(Symmetric(copy(real_admittance)))
    real_solution = physical_factor \ real.(rhs)
    imaginary_solution = physical_factor \ imag.(rhs)
    candidate_solution = complex.(real_solution, imaginary_solution)
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    components = tridiagonal_entries === nothing ?
        _harmonic_connected_components(admittance, structural_threshold) :
        _tridiagonal_harmonic_components(
            tridiagonal_entries.off_diagonal,
            structural_threshold,
        )
    referenced_components = tridiagonal_entries === nothing ?
        _harmonic_component_references(
            admittance,
            components,
            structural_threshold,
        ) : _tridiagonal_harmonic_component_references(
            tridiagonal_entries,
            components,
            structural_threshold,
        )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    return (
        solution=classification === :unique ? candidate_solution : nothing,
        numerical_rank=node_count,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

function _real_symmetric_tridiagonal_entries(
    admittance::Matrix{Float64},
)
    node_count = size(admittance, 1)
    size(admittance, 2) == node_count || return nothing
    for column in 3:node_count, row in 1:(column - 2)
        iszero(admittance[row, column]) || return nothing
    end
    diagonal = Vector{Float64}(undef, node_count)
    off_diagonal = Vector{Float64}(undef, max(node_count - 1, 0))
    for node in 1:node_count
        diagonal[node] = admittance[node, node]
    end
    for node in 1:(node_count - 1)
        off_diagonal[node] = admittance[node, node + 1]
    end
    return (; diagonal, off_diagonal)
end

function _tridiagonal_harmonic_components(
    off_diagonal::Vector{Float64},
    threshold::Float64,
)
    node_count = length(off_diagonal) + 1
    components = Vector{Vector{Int}}()
    first_node = 1
    for node in eachindex(off_diagonal)
        abs(off_diagonal[node]) > threshold && continue
        push!(components, collect(first_node:node))
        first_node = node + 1
    end
    push!(components, collect(first_node:node_count))
    return components
end

function _tridiagonal_harmonic_component_references(
    entries::NamedTuple,
    components::Vector{Vector{Int}},
    threshold::Float64,
)
    node_count = length(entries.diagonal)
    references = falses(length(components))
    for (component_index, component) in enumerate(components)
        for node in component
            row_sum = entries.diagonal[node]
            node > 1 && (row_sum += entries.off_diagonal[node - 1])
            node < node_count && (row_sum += entries.off_diagonal[node])
            if abs(row_sum) > threshold
                references[component_index] = true
                break
            end
        end
    end
    return references
end

function _solve_harmonic_linear_system(
    admittance::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64};
    current_absolute_a::Real,
    current_relative::Real,
    rank_relative_threshold_multiplier::Real,
    maximum_condition_estimate::Real,
)
    node_count = length(rhs)
    size(admittance) == (node_count, node_count) || throw(DimensionMismatch(
        "harmonic admittance and right-hand side dimensions must match",
    ))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), admittance) ||
        throw(ArgumentError("harmonic admittance must be finite"))
    all(value -> isfinite(real(value)) && isfinite(imag(value)), rhs) ||
        throw(ArgumentError("harmonic right-hand side must be finite"))
    absolute_tolerance = Float64(current_absolute_a)
    relative_tolerance = Float64(current_relative)
    rank_multiplier = Float64(rank_relative_threshold_multiplier)
    condition_limit = Float64(maximum_condition_estimate)
    isfinite(absolute_tolerance) && absolute_tolerance > 0.0 ||
        throw(ArgumentError("harmonic absolute-current tolerance must be finite and positive"))
    isfinite(relative_tolerance) && relative_tolerance > 0.0 ||
        throw(ArgumentError("harmonic relative-current tolerance must be finite and positive"))
    isfinite(rank_multiplier) && rank_multiplier > 0.0 ||
        throw(ArgumentError("harmonic rank multiplier must be finite and positive"))
    (isfinite(condition_limit) || condition_limit == Inf) && condition_limit > 1.0 ||
        throw(ArgumentError("harmonic condition limit must exceed one"))
    if node_count == 0
        return (
            solution=ComplexF64[],
            numerical_rank=0,
            condition_estimate=1.0,
            maximum_residual_a=0.0,
            relative_residual=0.0,
            connected_components=Vector{Vector{Int}}(),
            referenced_components=BitVector(),
            unreferenced_components=Vector{Vector{Int}}(),
            classification=:unique,
        )
    end

    matrix_scale = maximum(abs, admittance; init=0.0)
    structural_threshold = max(
        matrix_scale * node_count * eps(Float64) * rank_multiplier,
        floatmin(Float64),
    )
    components = _harmonic_connected_components(admittance, structural_threshold)
    referenced_components = _harmonic_component_references(
        admittance,
        components,
        structural_threshold,
    )
    unreferenced_components = Vector{Int}[
        copy(component)
        for (component, referenced) in zip(components, referenced_components)
        if !referenced
    ]

    symmetric_scales = Vector{Float64}(undef, node_count)
    for node in 1:node_count
        row_magnitude = sum(abs, view(admittance, node, :))
        column_magnitude = sum(abs, view(admittance, :, node))
        symmetric_scales[node] = sqrt(max(
            row_magnitude,
            column_magnitude,
            structural_threshold,
        ))
    end
    scaled_admittance = Matrix{ComplexF64}(undef, node_count, node_count)
    for column in 1:node_count, row in 1:node_count
        scaled_admittance[row, column] =
            admittance[row, column] /
            (symmetric_scales[row] * symmetric_scales[column])
    end
    scaled_rhs = rhs ./ symmetric_scales
    solve_components = _harmonic_connected_components(admittance, 0.0)
    decompositions = [
        svd(scaled_admittance[component, component])
        for component in solve_components
    ]
    singular_values = reduce(
        vcat,
        (decomposition.S for decomposition in decompositions);
        init=Float64[],
    )
    leading_singular_value = maximum(singular_values; init=0.0)
    rank_threshold = node_count * eps(Float64) *
        leading_singular_value * rank_multiplier
    numerical_rank = count(>(rank_threshold), singular_values)
    condition_estimate = numerical_rank == node_count ?
        leading_singular_value / minimum(singular_values) : Inf
    candidate_solution = zeros(ComplexF64, node_count)
    if numerical_rank == node_count
        physical_admittance = Matrix{ComplexF64}(admittance)
        physical_factor = lu(physical_admittance)
        candidate_solution .= physical_factor \ rhs
    else
        for (component, decomposition) in zip(solve_components, decompositions)
            component_rhs = scaled_rhs[component]
            all(iszero, component_rhs) && continue
            candidate_solution[component] .=
                (decomposition \ component_rhs) ./
                symmetric_scales[component]
        end
    end
    for component in solve_components
        if all(iszero, @view rhs[component])
            fill!(@view(candidate_solution[component]), 0.0 + 0.0im)
        end
    end
    residual = admittance * candidate_solution - rhs
    maximum_residual = norm(residual, Inf)
    rhs_scale = max(norm(rhs, Inf), absolute_tolerance)
    relative_residual = maximum_residual / rhs_scale
    roundoff_allowance = rank_multiplier * node_count * eps(Float64) * (
        norm(admittance, Inf) * norm(candidate_solution, Inf) + norm(rhs, Inf)
    )
    residual_tolerance = absolute_tolerance +
        relative_tolerance * norm(rhs, Inf) + roundoff_allowance
    classification = if maximum_residual > residual_tolerance
        :infeasible
    elseif numerical_rank < node_count
        isempty(unreferenced_components) ? :nonunique : :islanded
    elseif condition_estimate > condition_limit
        :ill_conditioned
    else
        :unique
    end
    solution = classification === :unique ? candidate_solution : nothing
    return (
        solution,
        numerical_rank,
        condition_estimate,
        maximum_residual_a=maximum_residual,
        relative_residual,
        connected_components=components,
        referenced_components,
        unreferenced_components,
        classification,
    )
end

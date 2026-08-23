function _stamp_deck_source_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    rhs::AbstractVector{ComplexF64},
    parsed::DeckParser.DeckParseResult,
)
    source_conductance_by_node = Dict{Int,Float64}()
    for element in parsed.elements
        element isa TheveninSource || continue
        source_conductance_by_node[Int(element.node)] = Float64(element.g)
    end
    for row in DeckParser.deck_over5a_source_rows(parsed)
        node_value = Int(row.node_value)
        target_node = abs(node_value)
        if node_value < 0
            abs(Int(row.iform)) in (11, 14) || continue
            rhs[target_node] += _deck_source_voltage_phasor(row)
            continue
        end
        conductance = get(source_conductance_by_node, target_node, 1.0e12)
        admittance = complex(conductance, 0.0)
        _stamp_complex_branch_admittance!(matrix, target_node, 0, admittance)
        rhs[target_node] += admittance * _deck_source_voltage_phasor(row)
    end
    return matrix
end

function _stamp_coupled_lumped_sequence_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    input_frequency_hz = DeckParser.deck_fixed_time_horizon_options(parsed).x_frequency_hz
    input_frequency_hz > 0.0 || (input_frequency_hz = default_frequency_hz)
    for impedance in DeckParser.deck_coupled_lumped_sequence_impedances(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(impedance.from_node_indices, impedance.to_node_indices),
            default_frequency_hz,
        )
        reactive_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        input_angular_frequency = 2.0 * pi * input_frequency_hz
        phase_impedance = ComplexF64.(
            impedance.phase_resistance_matrix,
            (reactive_frequency / input_angular_frequency) .*
            impedance.phase_inductance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            impedance.from_node_indices,
            impedance.to_node_indices,
            inv(phase_impedance),
        )
    end
    return matrix
end

function _generator_equivalent_modal_admittance(branch, angular_frequency::Float64)
    series_inductance = im * angular_frequency * branch.inductance_h
    damping_resistance = branch.damping_resistance_ohm
    damped_inductance =
        damping_resistance > 0.0 && branch.inductance_h > 0.0 ?
        (series_inductance * damping_resistance) /
        (series_inductance + damping_resistance) :
        series_inductance
    series_capacitance =
        branch.capacitance_f > 0.0 ?
        inv(im * angular_frequency * branch.capacitance_f) :
        complex(0.0, 0.0)
    return inv(branch.resistance_ohm + damped_inductance + series_capacitance)
end

function _stamp_generator_equivalent_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for row in DeckParser.deck_generator_equivalent_rows(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(row.from_node_indices, row.to_node_indices),
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        nph = length(row.from_node_indices)
        zero_admittance = sum(
            _generator_equivalent_modal_admittance(branch_row.branch, angular_frequency)
            for branch_row in row.zero_mode_branches
        )
        positive_admittance = sum(
            _generator_equivalent_modal_admittance(branch_row.branch, angular_frequency)
            for branch_row in row.positive_mode_branches
        )
        diagonal =
            (zero_admittance + (nph - 1) * positive_admittance) / nph
        off_diagonal = (zero_admittance - positive_admittance) / nph
        phase_admittance = fill(off_diagonal, nph, nph)
        for phase in 1:nph
            phase_admittance[phase, phase] = diagonal
        end
        _stamp_complex_phase_admittance!(
            matrix,
            row.from_node_indices,
            row.to_node_indices,
            phase_admittance,
        )
    end
    return matrix
end

function _stamp_coupled_lumped_phase_pi_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for section in DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(section.from_node_indices, section.to_node_indices),
            default_frequency_hz,
        )
        angular_frequency = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        )
        phase_impedance = ComplexF64.(
            section.phase_resistance_matrix,
            angular_frequency .* section.phase_inductance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            section.from_node_indices,
            section.to_node_indices,
            inv(phase_impedance),
        )
        phase_shunt = (0.5im * angular_frequency) .* section.phase_capacitance_matrix
        ground_nodes = zeros(Int, section.phase_count)
        _stamp_complex_phase_admittance!(
            matrix,
            section.from_node_indices,
            ground_nodes,
            phase_shunt,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            section.to_node_indices,
            ground_nodes,
            phase_shunt,
        )
    end
    return matrix
end

function _stamp_cascaded_phase_pi_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed);
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for equivalent in DeckParser.deck_cascaded_phase_pi_equivalents(parsed)
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(equivalent.from_node_indices, equivalent.to_node_indices),
            default_frequency_hz,
        )
        equivalent_frequency_hz = _emt_reactive_angular_frequency(
            formulation,
            2.0 * pi * frequency_hz,
        ) / (2.0 * pi)
        frequency_equivalent =
            isapprox(
                equivalent.frequency_hz,
                equivalent_frequency_hz;
                atol = 1.0e-12,
                rtol = 1.0e-12,
            ) ?
            equivalent :
            cascaded_phase_pi_equivalent(
                equivalent.name,
                equivalent.blocks,
                equivalent_frequency_hz,
            )
        node_indices = vcat(
            frequency_equivalent.from_node_indices,
            frequency_equivalent.to_node_indices,
        )
        for local_column in eachindex(node_indices)
            column = node_indices[local_column]
            column == 0 && continue
            for local_row in eachindex(node_indices)
                row = node_indices[local_row]
                row == 0 && continue
                matrix[row, column] +=
                    frequency_equivalent.terminal_admittance_s[
                        local_row,
                        local_column,
                    ]
            end
        end
    end
    return matrix
end

function _stamp_distributed_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    constants_rows = DeckParser.deck_distributed_transposed_line_constants(parsed)
    modal_states = DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)
    ground_nodes = fill(0, 3)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for (index, constants) in enumerate(constants_rows)
        modal_state = modal_states[index]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(modal_state.from_node_indices, modal_state.to_node_indices),
            default_frequency_hz,
        )
        equivalent = distributed_transposed_line_steady_state_pi_equivalent(
            constants;
            steady_state_frequency_hz = frequency_hz,
            storage_start_index = 1,
            name = Symbol("mixed_frequency_distributed_line_pi_", index),
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.from_node_indices,
            modal_state.to_node_indices,
            inv(equivalent.phase_series_impedance_matrix),
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.from_node_indices,
            ground_nodes,
            equivalent.phase_shunt_admittance_matrix,
        )
        _stamp_complex_phase_admittance!(
            matrix,
            modal_state.to_node_indices,
            ground_nodes,
            equivalent.phase_shunt_admittance_matrix,
        )
    end
    for element in parsed.elements
        element isa ComplexModalBergeronLine || continue
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(element.from_nodes, element.to_nodes),
            default_frequency_hz,
        )
        terminal_admittance =
            complex_modal_bergeron_steady_state_terminal_admittance(
                element,
                frequency_hz,
            )
        nodes = vcat(element.from_nodes, element.to_nodes)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] +=
                terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_semlyen_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    lines = vcat(
        DeckParser.deck_semlyen_line_elements(parsed, options.dt_s),
        DeckParser.deck_rational_frequency_line_elements(parsed, options.dt_s),
    )
    for line in lines
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            vcat(line.from_nodes, line.to_nodes),
            default_frequency_hz,
        )
        terminal_admittance =
            semlyen_line_steady_state_terminal_admittance(line, frequency_hz)
        nodes = vcat(line.from_nodes, line.to_nodes)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] += terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_sampled_frequency_line_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    parsed::DeckParser.DeckParseResult,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    default_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    for line in DeckParser.deck_sampled_frequency_line_elements(parsed, options.dt_s)
        nodes = line isa SampledFrequencyDependentLineGroup ?
            vcat(line.from_nodes, line.to_nodes) :
            [line.a, line.b]
        frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            nodes,
            default_frequency_hz,
        )
        terminal_admittance =
            sampled_line_steady_state_terminal_admittance(line, frequency_hz)
        for column in eachindex(nodes), row in eachindex(nodes)
            row_node = nodes[row]
            column_node = nodes[column]
            row_node == 0 && continue
            column_node == 0 && continue
            matrix[row_node, column_node] += terminal_admittance[row, column]
        end
    end
    return matrix
end

function _stamp_saturated_transformer_steady_state_admittance!(
    matrix::AbstractMatrix{ComplexF64},
    transformer_admittance,
)
    for index in eachindex(transformer_admittance.scalar_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.scalar_branch_from_node_indices[index],
            transformer_admittance.scalar_branch_to_node_indices[index],
            transformer_admittance.scalar_branch_admittances[index],
        )
    end
    for index in eachindex(transformer_admittance.ideal_low_branch_from_node_indices)
        _stamp_complex_phase_admittance!(
            matrix,
            [
                transformer_admittance.ideal_low_branch_from_node_indices[index],
                transformer_admittance.ideal_internal_branch_from_node_indices[index],
            ],
            [
                transformer_admittance.ideal_low_branch_to_node_indices[index],
                transformer_admittance.ideal_internal_branch_to_node_indices[index],
            ],
            ComplexF64[
                transformer_admittance.ideal_low_branch_admittances[index] transformer_admittance.ideal_mutual_branch_admittances[index]
                transformer_admittance.ideal_mutual_branch_admittances[index] transformer_admittance.ideal_internal_branch_admittances[index]
            ],
        )
    end
    for index in eachindex(transformer_admittance.linearized_nonlinear_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.linearized_nonlinear_branch_from_node_indices[index],
            transformer_admittance.linearized_nonlinear_branch_to_node_indices[index],
            transformer_admittance.linearized_nonlinear_branch_admittances[index],
        )
    end
    for index in eachindex(transformer_admittance.magnetizing_branch_from_node_indices)
        _stamp_complex_branch_admittance!(
            matrix,
            transformer_admittance.magnetizing_branch_from_node_indices[index],
            transformer_admittance.magnetizing_branch_to_node_indices[index],
            transformer_admittance.magnetizing_branch_admittances[index],
        )
    end
    return matrix
end

function _deck_steady_state_nodal_equations(
    parsed::DeckParser.DeckParseResult,
    node_count::Int;
    include_induction_machines::Bool=true,
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
    transformer_admittance=nothing,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition=
        DeckParser.deck_steady_state_frequency_partition(parsed),
    formulation::AbstractEMTHarmonicFormulation=PhysicalFrequencyFormulation(),
    default_frequency_hz::Float64=_deck_steady_state_frequency_hz(parsed),
)
    admittance = zeros(ComplexF64, node_count, node_count)
    rhs = zeros(ComplexF64, node_count)
    _stamp_deck_branch_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        excluded_universal_machine_indices =
            excluded_universal_machine_indices,
        formulation=formulation,
    )
    _stamp_hysteretic_inductor_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition;
        formulation=formulation,
    )
    _stamp_piecewise_nonlinear_inductor_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition;
        formulation=formulation,
    )
    include_induction_machines &&
        _stamp_deck_induction_machine_steady_state_admittance!(
            admittance,
            parsed,
            frequency_partition,
        )
    _stamp_deck_open_switch_steady_state_admittance!(admittance, parsed)
    _stamp_deck_source_steady_state_admittance!(admittance, rhs, parsed)
    _stamp_named_basic_harmonic_elements!(
        admittance,
        rhs,
        parsed,
        default_frequency_hz,
        formulation,
    )
    _stamp_coupled_lumped_sequence_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_generator_equivalent_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_coupled_lumped_phase_pi_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_cascaded_phase_pi_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
        formulation=formulation,
    )
    _stamp_distributed_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    _stamp_sampled_frequency_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    _stamp_semlyen_line_steady_state_admittance!(
        admittance,
        parsed,
        frequency_partition,
    )
    transformer_admittance === nothing ||
        _stamp_saturated_transformer_steady_state_admittance!(admittance, transformer_admittance)
    zeroed_universal_machine_branch_pairs = Tuple{Int,Int}[
        (row.from_node_value, row.to_node_value)
        for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
        if row.machine_index in excluded_universal_machine_indices &&
           row.from_node_value != row.to_node_value
    ]
    switch_representatives = _steady_state_closed_switch_representatives(
        parsed,
        node_count;
        additional_closed_node_pairs =
            zeroed_universal_machine_branch_pairs,
        include_grounded_constraints=true,
    )
    return admittance, rhs, switch_representatives
end

function _deck_external_steady_state_thevenin(
    parsed::DeckParser.DeckParseResult,
    node_index::Integer,
)
    node_count = maximum(values(parsed.node_map); init = 0)
    node = Int(node_index)
    1 <= node <= node_count ||
        throw(ArgumentError("steady-state Thevenin node is outside the deck network"))
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        include_induction_machines = false,
    )
    open_circuit_phasors = _solve_grouped_steady_state_admittance(
        admittance,
        rhs,
        switch_representatives,
    )
    unit_injection = zeros(ComplexF64, node_count)
    unit_injection[node] = 1.0
    driving_point_response = _solve_grouped_steady_state_admittance(
        admittance,
        unit_injection,
        switch_representatives,
    )
    return (
        source = :deck_external_steady_state_thevenin,
        node_index = node,
        voltage_phasor = open_circuit_phasors[node],
        impedance = driving_point_response[node],
        node_voltage_phasors = open_circuit_phasors,
        admittance = admittance,
        rhs = rhs,
        switch_representatives = switch_representatives,
    )
end

function _deck_external_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult,
    current_injections::AbstractDict{<:Integer,<:Complex},
    ;
    excluded_universal_machine_indices::AbstractVector{<:Integer}=Int[],
)
    node_map = Dict{Symbol,Int}(parsed.node_map)
    node_count = maximum(values(node_map); init = 0)
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        include_induction_machines = false,
        excluded_universal_machine_indices,
    )
    for (node_value, injection) in current_injections
        node = Int(node_value)
        1 <= node <= node_count ||
            throw(ArgumentError("steady-state current injection node is outside the deck network"))
        rhs[node] += ComplexF64(injection)
    end
    topology_diagnostics = _grouped_steady_state_diagnostics(
        admittance,
        rhs,
        switch_representatives,
    )
    node_voltage_phasors = something(topology_diagnostics.solution)
    return (
        source = :deck_external_steady_state_voltage_phasors,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_map = node_map,
        node_names = ordered_node_names(node_map),
        node_voltage_phasors = node_voltage_phasors,
        node_voltage_values = real.(node_voltage_phasors),
        topology_diagnostics,
        steady_state_admittance = admittance,
        steady_state_source_injection_phasors = rhs,
        steady_state_switch_representatives = switch_representatives,
    )
end

function deck_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult;
    saturated_transformer_intake = nothing,
    winding_number::Int = 1,
)
    node_map = Dict{Symbol,Int}(parsed.node_map)
    frequency_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    branch_assembly = nothing
    transformer_frequency_hz = _deck_steady_state_frequency_hz(parsed)
    if saturated_transformer_intake !== nothing
        arrays = saturated_transformer_nonlinear_arrays(saturated_transformer_intake)
        physical_node_map = _saturated_transformer_physical_node_map(parsed, arrays)
        transformer_frequency_hz = _steady_state_terminal_frequency_hz(
            frequency_partition,
            _saturated_transformer_frequency_node_indices(
                physical_node_map,
                arrays,
                length(frequency_partition.node_frequencies_hz),
            ),
            transformer_frequency_hz,
        )
        branch_assembly = saturated_transformer_winding_branch_assembly(
            arrays,
            physical_node_map;
            nonlinear_winding_number = winding_number,
            reference_node_index = 0,
        )
        node_map = _saturated_transformer_augmented_node_map(
            physical_node_map,
            branch_assembly,
        )
    end
    node_count = maximum(values(node_map); init = 0)
    transformer_admittance = saturated_transformer_intake === nothing ? nothing :
        saturated_transformer_steady_state_branch_admittance(
            saturated_transformer_intake,
            physical_node_map;
            nonlinear_winding_number = winding_number,
            reference_node_index = 0,
            reactance_units = 2.0 * pi * transformer_frequency_hz,
        )
    admittance, rhs, switch_representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        transformer_admittance = transformer_admittance,
        frequency_partition = frequency_partition,
    )
    topology_diagnostics = _grouped_steady_state_diagnostics(
        admittance,
        rhs,
        switch_representatives,
    )
    node_voltage_phasors = something(topology_diagnostics.solution)
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    output_voltage_values = Float64[
        index == 0 ? 0.0 : real(node_voltage_phasors[index])
        for index in output_node_indices
    ]
    return (
        source = :deck_steady_state_voltage_phasors,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_steady_state_frequencies_hz = copy(frequency_partition.node_frequencies_hz),
        node_frequency_source_row_indices =
            copy(frequency_partition.node_source_row_indices),
        source_frequency_successor_indices =
            copy(frequency_partition.source_successor_indices),
        steady_state_frequency_subnetwork_count =
            length(frequency_partition.subnetwork_node_indices),
        timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_map = node_map,
        node_names = ordered_node_names(node_map),
        node_voltage_phasors = node_voltage_phasors,
        node_voltage_values = real.(node_voltage_phasors),
        topology_diagnostics,
        steady_state_admittance = admittance,
        steady_state_source_injection_phasors = rhs,
        steady_state_switch_representatives = switch_representatives,
        output_node_indices = output_node_indices,
        output_voltage_values = output_voltage_values,
        source_row_count = length(DeckParser.deck_over5a_source_rows(parsed)),
        coupled_lumped_sequence_count =
            length(DeckParser.deck_coupled_lumped_sequence_impedances(parsed)),
        generator_equivalent_count =
            length(DeckParser.deck_generator_equivalent_rows(parsed)),
        coupled_lumped_phase_pi_section_count =
            length(DeckParser.deck_coupled_lumped_phase_pi_sections(parsed)),
        cascaded_phase_pi_equivalent_count =
            length(DeckParser.deck_cascaded_phase_pi_equivalents(parsed)),
        distributed_transposed_line_count =
            length(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)),
        saturated_transformer_branch_count =
            transformer_admittance === nothing ? 0 : transformer_admittance.branch_count,
        saturated_transformer_ideal_branch_count =
            transformer_admittance === nothing ? 0 : transformer_admittance.ideal_branch_count,
        saturated_transformer_linearized_nonlinear_branch_count =
            transformer_admittance === nothing ? 0 :
            transformer_admittance.linearized_nonlinear_branch_count,
        saturated_transformer_ideal_primary_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_primary_storage_coefficients),
        saturated_transformer_ideal_mutual_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_mutual_storage_coefficients),
        saturated_transformer_ideal_secondary_storage_coefficients =
            transformer_admittance === nothing ? Float64[] :
            copy(transformer_admittance.ideal_secondary_storage_coefficients),
    )
end

function _deck_node_voltage_initial_sample(parsed::DeckParser.DeckParseResult)
    rows = [
        row for row in DeckParser.deck_node_initial_condition_rows(parsed)
        if row.condition_kind == :node_voltage_initial_condition && row.node_index > 0
    ]
    isempty(rows) && return nothing
    values = zeros(Float64, length(parsed.node_map))
    for row in rows
        reference_value =
            row.reference_node_index isa Missing ? 0.0 : values[Int(row.reference_node_index)]
        values[Int(row.node_index)] = reference_value + Float64(row.real_value)
    end
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    output_voltage_values = Float64[
        index == 0 ? 0.0 : values[index]
        for index in output_node_indices
    ]
    return (
        source = :deck_node_voltage_initial_conditions,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        output_node_indices = output_node_indices,
        output_voltage_values = output_voltage_values,
        node_voltage_condition_count = length(rows),
    )
end

function _deck_zero_initial_voltage_sample(parsed::DeckParser.DeckParseResult)
    values = zeros(Float64, length(parsed.node_map))
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    return (
        source = :zero_initial_network_state,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        output_node_indices = output_node_indices,
        output_voltage_values = zeros(Float64, length(output_node_indices)),
        node_voltage_condition_count = 0,
    )
end

function _synchronous_machine_terminal_phase_angle_deg(row, base_angle_deg::Float64)
    row.angle_deg isa Missing || return Float64(row.angle_deg)
    row.phase_index == 1 && return base_angle_deg
    row.phase_index == 2 && return base_angle_deg + 240.0
    row.phase_index == 3 && return base_angle_deg + 120.0
    return base_angle_deg
end

function _deck_synchronous_machine_terminal_voltage_initial_sample(
    parsed::DeckParser.DeckParseResult,
)
    rows = DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    isempty(rows) && return nothing
    values = zeros(Float64, length(parsed.node_map))
    phasors = zeros(ComplexF64, length(parsed.node_map))
    source_count = 0
    for machine_index in unique(row.machine_index for row in rows)
        machine_rows = [row for row in rows if row.machine_index == machine_index]
        reference_row = findfirst(row -> !(row.peak_terminal_voltage isa Missing), machine_rows)
        reference_row === nothing && continue
        reference = machine_rows[reference_row]
        peak_voltage = Float64(reference.peak_terminal_voltage)
        base_angle_deg =
            reference.angle_deg isa Missing ? 0.0 : Float64(reference.angle_deg)
        for row in machine_rows
            node = row.terminal_node_value
            1 <= node <= length(values) || continue
            local_peak =
                row.peak_terminal_voltage isa Missing ?
                peak_voltage :
                Float64(row.peak_terminal_voltage)
            angle_deg = _synchronous_machine_terminal_phase_angle_deg(
                row,
                base_angle_deg,
            )
            angle_rad = deg2rad(angle_deg)
            phasor = local_peak * ComplexF64(cos(angle_rad), sin(angle_rad))
            phasors[node] = phasor
            values[node] = real(phasor)
            source_count += 1
        end
    end
    source_count > 0 || return nothing
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    return (
        source = :synchronous_machine_terminal_voltage_initial_sample,
        outcome = :initial_voltage_sample,
        steady_state_frequency_hz = _deck_steady_state_frequency_hz(parsed),
        node_voltage_values = values,
        node_voltage_phasors = phasors,
        output_node_indices = output_node_indices,
        output_voltage_values = Float64[
            index == 0 ? 0.0 : values[index]
            for index in output_node_indices
        ],
        synchronous_machine_terminal_voltage_count = source_count,
    )
end

function _stamp_synchronous_machine_transformer_steady_state!(
    admittance::Matrix{ComplexF64},
    context,
    angular_frequency::Float64,
)
    Base.@nospecialize context
    for (name, element) in zip(context.element_names, context.system.elements)
        name_text = String(name)
        if startswith(name_text, "saturated_transformer_")
            if element isa SeriesRLBranch
                branch_admittance = inv(complex(
                    element.r,
                    angular_frequency * element.l,
                ))
                _stamp_complex_branch_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    branch_admittance,
                )
            elseif element isa CoupledInductiveBranch
                _stamp_complex_phase_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    _coupled_inductive_steady_state_admittance(element),
                )
            elseif element isa ConductanceBranch
                _stamp_complex_branch_admittance!(
                    admittance,
                    element.a,
                    element.b,
                    complex(element.g, 0.0),
                )
            end
        elseif startswith(name_text, "transformer_branch_shunt_capacitance_") &&
               element isa CapacitorBranch
            _stamp_complex_branch_admittance!(
                admittance,
                element.a,
                element.b,
                complex(0.0, angular_frequency * element.c),
            )
        end
    end
    return admittance
end

function _deck_synchronous_machine_network_initial_sample(
    parsed::DeckParser.DeckParseResult,
    context,
    ;
    strict_topology_classification::Bool=false,
)
    Base.@nospecialize context
    terminal_sample = _deck_synchronous_machine_terminal_voltage_initial_sample(parsed)
    terminal_sample === nothing && return nothing
    node_count = context.system.node_count
    admittance = zeros(ComplexF64, node_count, node_count)
    rhs = zeros(ComplexF64, node_count)
    Base.inferencebarrier(_stamp_deck_branch_steady_state_admittance!)(
        admittance,
        parsed,
    )
    time_zero_ground_fault = any(
        row -> _deck_time_switch_closed_at(
                   row.initially_closed,
                   Float64(row.close_time_s),
                   Float64(row.open_time_s),
                   0.0,
               ) && (row.from_node_value == 0 || row.to_node_value == 0),
        DeckParser.deck_over5_switch_rows(parsed),
    )
    external_excitation_port_initialization = any(
        row -> row.coupling_kind in
               (:exciter_voltage_output, :exciter_current_output),
        DeckParser.deck_synchronous_machine_control_interface_rows(parsed),
    )
    use_time_zero_topology =
        time_zero_ground_fault && external_excitation_port_initialization
    if use_time_zero_topology
        Base.inferencebarrier(_stamp_deck_switch_at_time_steady_state_admittance!)(
            admittance,
            parsed,
            0.0,
        )
    else
        Base.inferencebarrier(_stamp_deck_open_switch_steady_state_admittance!)(
            admittance,
            parsed,
        )
    end
    Base.inferencebarrier(_stamp_deck_source_steady_state_admittance!)(
        admittance,
        rhs,
        parsed,
    )
    Base.inferencebarrier(_stamp_coupled_lumped_sequence_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_generator_equivalent_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_coupled_lumped_phase_pi_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_cascaded_phase_pi_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_distributed_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_sampled_frequency_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    Base.inferencebarrier(_stamp_semlyen_line_steady_state_admittance!)(
        admittance,
        parsed,
    )
    frequency_hz = Float64(terminal_sample.steady_state_frequency_hz)
    angular_frequency = 2.0 * pi * frequency_hz
    Base.inferencebarrier(_stamp_synchronous_machine_transformer_steady_state!)(
        admittance,
        context,
        angular_frequency,
    )

    terminal_rows = DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    fixed_nodes = use_time_zero_topology ? Int[] :
        sort(unique(Int(row.terminal_node_value) for row in terminal_rows))
    fixed_phasors = ComplexF64[terminal_sample.node_voltage_phasors[node] for node in fixed_nodes]
    topology_diagnostics = nothing
    switch_representatives = nothing
    phasors = zeros(ComplexF64, node_count)
    if strict_topology_classification
        switch_representatives = _steady_state_closed_switch_representatives(
            parsed,
            node_count,
            0.0,
            include_grounded_constraints=true,
        )
        topology_diagnostics = if use_time_zero_topology
            _grouped_steady_state_diagnostics(
                admittance,
                rhs,
                switch_representatives,
            )
        else
            fixed_node_phasors = Dict{Int,ComplexF64}(
                node => phasor for (node, phasor) in zip(fixed_nodes, fixed_phasors)
            )
            result = _solve_grouped_constrained_harmonic_linear_system(
                admittance,
                rhs,
                switch_representatives,
                fixed_node_phasors;
                current_absolute_a=1.0e-12,
                current_relative=1.0e-10,
                rank_relative_threshold_multiplier=10.0,
                maximum_condition_estimate=Inf,
            )
            result.classification === :unique || throw(ArgumentError(
                "constrained synchronous-machine steady-state network classification " *
                "$(result.classification): rank $(result.numerical_rank)/" *
                "$(result.reduced_node_count), condition $(result.condition_estimate), " *
                "residual $(result.maximum_residual_a) A",
            ))
            result
        end
        phasors .= something(topology_diagnostics.solution)
    else
        unknown_nodes = setdiff(collect(1:node_count), fixed_nodes)
        phasors[fixed_nodes] .= fixed_phasors
        if use_time_zero_topology
            switch_representatives = _steady_state_closed_switch_representatives(
                parsed,
                node_count,
                0.0,
            )
            phasors .= _solve_grouped_steady_state_admittance(
                admittance,
                rhs,
                switch_representatives,
            )
        elseif !isempty(unknown_nodes)
            reduced_rhs = rhs[unknown_nodes] -
                          admittance[unknown_nodes, fixed_nodes] * fixed_phasors
            phasors[unknown_nodes] .= _solve_steady_state_linear_system(
                admittance[unknown_nodes, unknown_nodes],
                reduced_rhs,
            )
        end
    end
    node_current_phasors = admittance * phasors - rhs
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    sample = (
        source = use_time_zero_topology ?
            :synchronous_machine_time_zero_topology_steady_state :
            :synchronous_machine_constrained_steady_state,
        outcome = :steady_state_initial_voltage_sample,
        steady_state_frequency_hz = frequency_hz,
        node_voltage_values = real.(phasors),
        node_voltage_phasors = phasors,
        node_current_phasors = node_current_phasors,
        output_node_indices = output_node_indices,
        output_voltage_values = Float64[
            index == 0 ? 0.0 : real(phasors[index])
            for index in output_node_indices
        ],
        synchronous_machine_terminal_voltage_count = length(fixed_nodes),
        time_zero_ground_fault,
        external_excitation_port_initialization,
    )
    strict_topology_classification || return sample
    return merge(
        sample,
        (;
            topology_diagnostics,
            steady_state_admittance=admittance,
            steady_state_source_injection_phasors=rhs,
            steady_state_switch_representatives=switch_representatives,
        ),
    )
end

function _initial_voltage_sample_for_context(sample, node_count::Int)
    sample === nothing && return nothing
    values = Float64.(sample.node_voltage_values)
    length(values) >= node_count && return sample
    padded_values = zeros(Float64, node_count)
    copyto!(padded_values, 1, values, 1, length(values))
    if hasproperty(sample, :node_voltage_phasors)
        phasors = ComplexF64.(sample.node_voltage_phasors)
        padded_phasors = zeros(ComplexF64, node_count)
        copyto!(padded_phasors, 1, phasors, 1, length(phasors))
        return merge(
            sample,
            (
                node_voltage_values = padded_values,
                node_voltage_phasors = padded_phasors,
            ),
        )
    end
    return merge(sample, (node_voltage_values = padded_values,))
end

function _current_injection_samples_for_context(
    current_injection_samples,
    node_count::Int,
    sample_count::Int,
)
    current_injection_samples === nothing && return nothing
    values = Matrix{Float64}(current_injection_samples)
    size(values, 1) >= node_count ||
        throw(ArgumentError("current injection samples must cover every node"))
    size(values, 2) >= sample_count ||
        throw(ArgumentError("current injection samples must cover every trace sample"))
    return values
end

function _deck_runtime_saturated_transformer_intake(
    parsed::DeckParser.DeckParseResult,
)
    get(parsed.card_counts, :fixed_card_saturated_transformer_intake, 0) > 0 ||
        return nothing
    source_path = DeckParser.deck_source_path(parsed)
    source_path === nothing && return nothing
    isfile(source_path) || return nothing
    intake = DeckParser.parse_saturated_transformer_branch_section_intake_file(source_path)
    assert_valid!(intake.validation)
    return intake
end

function _deck_runtime_initial_voltage_sample(
    parsed::DeckParser.DeckParseResult,
    initial_voltage_source::Symbol,
)
    if initial_voltage_source == :none
        return nothing
    elseif initial_voltage_source == :zero
        return _deck_zero_initial_voltage_sample(parsed)
    elseif initial_voltage_source == :node_conditions
        return _deck_node_voltage_initial_sample(parsed)
    elseif initial_voltage_source == :steady_state
        return _deck_runtime_steady_state_voltage_phasors(parsed)
    elseif initial_voltage_source == :synchronous_machine_terminals
        return _deck_synchronous_machine_terminal_voltage_initial_sample(parsed)
    end
    throw(ArgumentError("unsupported initial voltage source $(initial_voltage_source)"))
end

function _deck_runtime_steady_state_voltage_phasors(
    parsed::DeckParser.DeckParseResult,
)
    saturated_transformer_intake = _deck_runtime_saturated_transformer_intake(parsed)
    saturated_transformer_intake === nothing &&
        return deck_steady_state_voltage_phasors(parsed)
    try
        return deck_steady_state_voltage_phasors(
            parsed;
            saturated_transformer_intake = saturated_transformer_intake,
        )
    catch err
        message = sprint(showerror, err)
        if err isa ArgumentError &&
           occursin("missing saturated transformer winding", message)
            return deck_steady_state_voltage_phasors(parsed)
        end
        rethrow()
    end
end

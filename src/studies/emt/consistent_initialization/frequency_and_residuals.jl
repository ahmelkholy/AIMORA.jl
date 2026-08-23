struct _EMTInitializationRefusal <: Exception
    code::Symbol
    owner::Symbol
    quantity::Symbol
    message::String
    context::NamedTuple
end

function Base.showerror(io::IO, refusal::_EMTInitializationRefusal)
    print(io, refusal.message)
end

function _empty_emt_initialization_topology_report()
    return EMTInitializationTopologyReport(
        0,
        0,
        Vector{Vector{Int}}(),
        Vector{Vector{Int}}(),
        BitVector(),
        Vector{Vector{Int}}(),
        0,
        Inf,
        Inf,
        Inf,
        :not_evaluated,
    )
end

function _emt_initialization_topology_report(diagnostics, node_count::Int)
    return EMTInitializationTopologyReport(
        node_count,
        diagnostics.reduced_node_count,
        copy.(diagnostics.switch_node_groups),
        copy.(diagnostics.connected_components),
        copy(diagnostics.referenced_components),
        copy.(diagnostics.unreferenced_components),
        diagnostics.numerical_rank,
        diagnostics.condition_estimate,
        diagnostics.maximum_residual_a,
        diagnostics.relative_residual,
        diagnostics.classification,
    )
end

function _emt_initialization_classification_failure(topology::EMTInitializationTopologyReport)
    code = topology.classification in (
        :islanded,
        :nonunique,
        :infeasible,
        :ill_conditioned,
    ) ? topology.classification : :network_equilibrium
    return _EMTInitializationRefusal(
        code,
        :network_topology,
        :nodal_voltage,
        "EMT initialization network classification $(topology.classification): " *
        "rank $(topology.numerical_rank)/$(topology.reduced_node_count), " *
        "condition $(topology.condition_estimate), residual " *
        "$(topology.maximum_residual_a) A",
        (
            node_count=topology.node_count,
            reduced_node_count=topology.reduced_node_count,
            unreferenced_components=copy.(topology.unreferenced_components),
        ),
    )
end

function _emt_admittance_symmetry_error(
    admittance::AbstractMatrix{ComplexF64},
)
    error = 0.0
    for column in axes(admittance, 2), row in axes(admittance, 1)
        error = max(
            error,
            abs(admittance[row, column] - admittance[column, row]),
        )
    end
    return error
end

function _emt_minimum_dissipative_eigenvalue(
    admittance::AbstractMatrix{ComplexF64},
)
    row_count, column_count = size(admittance)
    row_count == column_count || throw(DimensionMismatch(
        "dissipative admittance diagnostic requires a square matrix",
    ))
    if all(value -> iszero(imag(value)), admittance)
        dissipative = Matrix{Float64}(undef, row_count, column_count)
        for column in 1:column_count, row in 1:row_count
            dissipative[row, column] = 0.5 * (
                real(admittance[row, column]) +
                real(admittance[column, row])
            )
        end
        return minimum(eigvals!(Symmetric(dissipative)); init=0.0)
    end
    dissipative = Matrix{ComplexF64}(undef, row_count, column_count)
    for column in 1:column_count, row in 1:row_count
        dissipative[row, column] = 0.5 * (
            admittance[row, column] + conj(admittance[column, row])
        )
    end
    return minimum(eigvals!(Hermitian(dissipative)); init=0.0)
end

function _emt_initialization_frequency_point_evidence(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition;
    frequency_assignment::Symbol,
    apply_operating_constraints::Bool,
    passive_conductance_network::Bool=false,
)
    node_count = maximum(values(parsed.node_map); init=0)
    node_count > 0 || throw(_EMTInitializationRefusal(
        :missing_network,
        :network_topology,
        :node,
        "EMT initialization requires at least one network node",
        (source=parsed.source,),
    ))
    admittance, rhs, representatives = _deck_steady_state_nodal_equations(
        parsed,
        node_count;
        frequency_partition,
        formulation=request.formulation,
        default_frequency_hz=frequency_hz,
    )
    fixed_node_phasors = apply_operating_constraints ?
        _emt_operating_point_voltage_constraints(
            parsed,
            request,
            frequency_hz,
            representatives,
            frequency_partition,
        ) : Dict{Int,ComplexF64}()
    tolerances = request.tolerances
    natural_diagnostics = if passive_conductance_network &&
                             !isempty(fixed_node_phasors)
        _solve_grouped_harmonic_linear_system(
            admittance,
            rhs,
            representatives;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    else
        nothing
    end
    redundant_operating_constraints =
        natural_diagnostics !== nothing &&
        natural_diagnostics.solution !== nothing &&
        all(fixed_node_phasors) do (representative, target)
            natural = natural_diagnostics.solution[representative]
            allowance = tolerances.voltage_absolute_v +
                tolerances.voltage_relative * max(abs(natural), abs(target))
            abs(natural - target) <= allowance
        end
    diagnostics = if redundant_operating_constraints
        natural_diagnostics
    elseif isempty(fixed_node_phasors)
        _solve_grouped_harmonic_linear_system(
            admittance,
            rhs,
            representatives;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    else
        _solve_grouped_constrained_harmonic_linear_system(
            admittance,
            rhs,
            representatives,
            fixed_node_phasors;
            current_absolute_a=tolerances.current_absolute_a,
            current_relative=tolerances.current_relative,
            rank_relative_threshold_multiplier=
                tolerances.rank_relative_threshold_multiplier,
            maximum_condition_estimate=tolerances.maximum_condition_estimate,
            passive_conductance_network,
        )
    end
    topology = _emt_initialization_topology_report(diagnostics, node_count)
    symmetry_error = _emt_admittance_symmetry_error(admittance)
    minimum_dissipative_eigenvalue =
        passive_conductance_network && topology.classification === :unique ?
        0.0 : _emt_minimum_dissipative_eigenvalue(admittance)
    physical_angular_frequency = 2.0 * pi * frequency_hz
    reactive_angular_frequency = _emt_reactive_angular_frequency(
        request.formulation,
        physical_angular_frequency,
    )
    node_voltage_phasors = diagnostics.solution === nothing ? ComplexF64[] :
        ComplexF64.(diagnostics.solution)
    source_injection_phasors = ComplexF64.(rhs)
    operating_constraint_current_phasors =
        hasproperty(diagnostics, :constraint_reaction_current_phasors) ?
        ComplexF64.(diagnostics.constraint_reaction_current_phasors) :
        zeros(ComplexF64, node_count)
    if request.time_origin_s != 0.0
        for node in 1:node_count
            rotation = cis(
                2.0 * pi *
                frequency_partition.node_frequencies_hz[node] *
                request.time_origin_s,
            )
            node_voltage_phasors[node] *= rotation
            source_injection_phasors[node] *= rotation
            operating_constraint_current_phasors[node] *= rotation
        end
    end
    matrix_scale = max(norm(admittance, Inf), 1.0)
    symmetry_passed = symmetry_error <= 1.0e-11 * matrix_scale
    passivity_floor = -1.0e-11 * matrix_scale
    passed = topology.classification === :unique &&
        symmetry_passed &&
        minimum_dissipative_eigenvalue >= passivity_floor &&
        !isempty(node_voltage_phasors)
    point = EMTInitializationFrequencyPoint(
        _emt_harmonic_formulation_symbol(request.formulation),
        frequency_assignment,
        frequency_hz,
        reactive_angular_frequency,
        copy(frequency_partition.node_frequencies_hz),
        copy(frequency_partition.node_source_row_indices),
        copy(frequency_partition.source_successor_indices),
        length(frequency_partition.subnetwork_node_indices),
        node_voltage_phasors,
        source_injection_phasors,
        operating_constraint_current_phasors,
        topology,
        symmetry_error,
        minimum_dissipative_eigenvalue,
        passed,
    )
    return (;
        point,
        natural_equilibrium_confirmed=redundant_operating_constraints ||
            isempty(fixed_node_phasors),
    )
end

function _emt_initialization_frequency_point(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition;
    frequency_assignment::Symbol,
    apply_operating_constraints::Bool,
    passive_conductance_network::Bool=false,
)
    return _emt_initialization_frequency_point_evidence(
        parsed,
        request,
        frequency_hz,
        frequency_partition;
        frequency_assignment,
        apply_operating_constraints,
        passive_conductance_network,
    ).point
end

function _emt_initialization_frequency_point(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    return _emt_initialization_frequency_point(
        parsed,
        request,
        frequency_hz,
        _network_frequency_scan_partition(parsed, frequency_hz);
        frequency_assignment=:uniform_scan,
        apply_operating_constraints=false,
        passive_conductance_network=
            _emt_frequency_invariant_passive_conductance_network(parsed),
    )
end

function _emt_operating_point_quantity_target(
    quantity::OperatingPointQuantity,
    operating_point::EMTOperatingPoint,
)
    quantity.quantity === :node_voltage_peak_phasor || throw(
        _EMTInitializationRefusal(
            :unsupported_quantity,
            :operating_point_mapping,
            quantity.quantity,
            "only explicit peak node-voltage phasors are admitted for EMT operating-point constraints",
            (asset=quantity.asset, phase=quantity.phase),
        ),
    )
    quantity.unit in ("V", "kV", "MV", "pu") || throw(
        _EMTInitializationRefusal(
            :unknown_unit,
            :operating_point_mapping,
            quantity.quantity,
            "unsupported operating-point voltage unit $(quantity.unit)",
            (asset=quantity.asset, unit=quantity.unit),
        ),
    )
    quantity.basis in (
        "peak_node_to_ground",
        "peak_phase_to_ground",
        "absolute_si_peak",
        "per_unit_peak_node_to_ground",
    ) || throw(_EMTInitializationRefusal(
        :unknown_basis,
        :operating_point_mapping,
        quantity.quantity,
        "unsupported operating-point voltage basis $(quantity.basis)",
        (asset=quantity.asset, basis=quantity.basis),
    ))
    quantity.orientation == "node_to_ground" || throw(
        _EMTInitializationRefusal(
            :unknown_orientation,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point node voltage must use node_to_ground orientation",
            (asset=quantity.asset, orientation=quantity.orientation),
        ),
    )
    quantity.phase === :not_applicable ||
        quantity.phase in operating_point.phase_order || throw(
        _EMTInitializationRefusal(
            :invalid_phase,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point phase is absent from the declared phase order",
            (asset=quantity.asset, phase=quantity.phase),
        ),
    )
    quantity.provenance.units == quantity.unit || throw(
        _EMTInitializationRefusal(
            :provenance_mismatch,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity and provenance units do not match",
            (
                asset=quantity.asset,
                quantity_unit=quantity.unit,
                provenance_units=quantity.provenance.units,
            ),
        ),
    )
    quantity.provenance.nature in (
        PhysicalModelParameter,
        ScalingBasisParameter,
    ) || throw(_EMTInitializationRefusal(
        :wrong_parameter_nature,
        :operating_point_mapping,
        quantity.quantity,
        "an imported operating point cannot be numerical-policy data",
        (asset=quantity.asset, nature=quantity.provenance.nature),
    ))
    return quantity.orientation_sign * quantity.scale_to_si * quantity.value
end

function _emt_operating_point_voltage_constraints(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
    representatives::AbstractVector{<:Integer},
    frequency_partition::DeckParser.DeckSteadyStateFrequencyPartition,
)
    operating_point = request.operating_point
    operating_point === nothing && return Dict{Int,ComplexF64}()
    frequency_hz == request.frequency_hz || return Dict{Int,ComplexF64}()
    operating_point isa EMTOperatingPoint || throw(_EMTInitializationRefusal(
        :unsupported_schema,
        :operating_point_mapping,
        :schema,
        "operating_point must be an EMTOperatingPoint",
        (runtime_type=string(typeof(operating_point)),),
    ))
    _validate_emt_operating_point_signatures(operating_point, request)
    fixed_by_representative = Dict{Int,Tuple{ComplexF64,Float64,Symbol}}()
    for quantity in operating_point.quantities
        node = get(parsed.node_map, quantity.asset, 0)
        node > 0 || throw(_EMTInitializationRefusal(
            :missing_asset,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point asset $(quantity.asset) is not a network node",
            (asset=quantity.asset,),
        ))
        node_frequency_hz = frequency_partition.node_frequencies_hz[node]
        isapprox(
            node_frequency_hz,
            operating_point.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity $(quantity.asset) belongs to a different frequency subnetwork",
            (
                asset=quantity.asset,
                operating_point_frequency_hz=operating_point.frequency_hz,
                node_frequency_hz,
            ),
        ))
        target_at_origin = _emt_operating_point_quantity_target(
            quantity,
            operating_point,
        )
        target_reference = target_at_origin /
            cis(2.0 * pi * node_frequency_hz * request.time_origin_s)
        representative = Int(representatives[node])
        uncertainty = quantity.absolute_uncertainty * quantity.scale_to_si
        if haskey(fixed_by_representative, representative)
            previous, previous_uncertainty, previous_asset =
                fixed_by_representative[representative]
            allowance = request.tolerances.voltage_absolute_v +
                request.tolerances.voltage_relative *
                    max(abs(previous), abs(target_reference)) +
                previous_uncertainty + uncertainty
            abs(previous - target_reference) <= allowance || throw(
                _EMTInitializationRefusal(
                    :mode_inconsistent,
                    :operating_point_mapping,
                    quantity.quantity,
                    "closed-switch node group has inconsistent operating-point voltages",
                    (
                        first_asset=previous_asset,
                        second_asset=quantity.asset,
                        residual_v=abs(previous - target_reference),
                        allowance_v=allowance,
                    ),
                ),
            )
        else
            fixed_by_representative[representative] = (
                target_reference,
                uncertainty,
                quantity.asset,
            )
        end
    end
    return Dict{Int,ComplexF64}(
        representative => value[1]
        for (representative, value) in fixed_by_representative
    )
end

function _emt_initialization_scan(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    frequencies = sort!(unique(vcat(
        request.frequency_hz,
        request.frequency_grid_hz,
    )))
    declared_partition = DeckParser.deck_steady_state_frequency_partition(parsed)
    declared_source_frequencies = sort!(unique(Float64[
        declared_partition.source_frequencies_hz[source_index]
        for source_index in declared_partition.active_source_row_indices
    ]))
    mixed_declared_subnetworks = length(declared_source_frequencies) > 1
    primary_partition = mixed_declared_subnetworks ? declared_partition :
        _network_frequency_scan_partition(parsed, request.frequency_hz)
    invariant_named_network = !mixed_declared_subnetworks &&
        _emt_frequency_invariant_named_network(parsed)
    passive_conductance_network = invariant_named_network &&
        _emt_frequency_invariant_passive_conductance_network(parsed)
    primary_evidence = _emt_initialization_frequency_point_evidence(
        parsed,
        request,
        request.frequency_hz,
        primary_partition;
        frequency_assignment=:initial_operating_point,
        apply_operating_constraints=true,
        passive_conductance_network,
    )
    primary_point = primary_evidence.point
    invariant_source_frequencies_hz = invariant_named_network ? Float64[
        element.value.frequency
        for element in parsed.elements
        if element isa Union{TheveninSource,CurrentInjection}
    ] : Float64[]
    points = EMTInitializationFrequencyPoint[primary_point]
    invariant_unforced_point = nothing
    for frequency_hz in frequencies
        !mixed_declared_subnetworks && frequency_hz == request.frequency_hz &&
            continue
        source_active = any(invariant_source_frequencies_hz) do source_frequency_hz
            isapprox(
                frequency_hz,
                source_frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            )
        end
        if invariant_named_network &&
           primary_evidence.natural_equilibrium_confirmed &&
           !source_active
            if invariant_unforced_point === nothing
                invariant_unforced_point =
                    _emt_frequency_invariant_unforced_scan_point(
                        primary_point,
                        request,
                        frequency_hz,
                    )
                push!(points, invariant_unforced_point)
            else
                push!(
                    points,
                    _emt_frequency_invariant_scan_point(
                        invariant_unforced_point,
                        request,
                        frequency_hz,
                    ),
                )
            end
            continue
        end
        frequency_partition = _network_frequency_scan_partition(
            parsed,
            frequency_hz,
        )
        point = _emt_initialization_frequency_point(
            parsed,
            request,
            frequency_hz,
            frequency_partition;
            frequency_assignment=:uniform_scan,
            apply_operating_constraints=false,
            passive_conductance_network,
        )
        push!(
            points,
            point,
        )
    end
    return points
end

function _emt_frequency_invariant_named_network(
    parsed::DeckParser.DeckParseResult,
)
    isempty(parsed.elements) && return false
    all(
        element -> element isa Union{
            ConductanceBranch,
            TheveninSource,
            CurrentInjection,
        },
        parsed.elements,
    ) || return false
    owners = _basic_harmonic_element_owners(parsed)
    return length(owners) == length(parsed.elements) &&
        all(owner -> owner === :named_basic_element, owners)
end

function _emt_frequency_invariant_passive_conductance_network(
    parsed::DeckParser.DeckParseResult,
)
    _emt_frequency_invariant_named_network(parsed) || return false
    return all(parsed.elements) do element
        if element isa ConductanceBranch
            isfinite(element.g) && element.g >= 0.0
        elseif element isa TheveninSource
            isfinite(element.g) && element.g >= 0.0
        else
            element isa CurrentInjection
        end
    end
end

function _emt_frequency_invariant_scan_point(
    template::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    node_count = length(template.node_voltage_phasors)
    node_physical_frequencies_hz = fill(frequency_hz, node_count)
    node_voltage_phasors = template.node_voltage_phasors
    source_injection_phasors = template.source_injection_phasors
    operating_constraint_current_phasors =
        template.operating_constraint_current_phasors
    if request.time_origin_s != 0.0
        node_voltage_phasors = copy(node_voltage_phasors)
        source_injection_phasors = copy(source_injection_phasors)
        operating_constraint_current_phasors =
            copy(operating_constraint_current_phasors)
        for node in 1:node_count
            rotation = cis(
                2.0 * pi *
                (node_physical_frequencies_hz[node] -
                 template.node_physical_frequencies_hz[node]) *
                request.time_origin_s,
            )
            node_voltage_phasors[node] *= rotation
            source_injection_phasors[node] *= rotation
            operating_constraint_current_phasors[node] *= rotation
        end
    end
    return EMTInitializationFrequencyPoint(
        template.formulation,
        :uniform_scan,
        frequency_hz,
        _emt_reactive_angular_frequency(
            request.formulation,
            2.0 * pi * frequency_hz,
        ),
        node_physical_frequencies_hz,
        template.node_frequency_source_row_indices,
        template.source_frequency_successor_indices,
        template.frequency_subnetwork_count,
        node_voltage_phasors,
        source_injection_phasors,
        operating_constraint_current_phasors,
        template.topology,
        template.admittance_symmetry_max_abs_error,
        template.minimum_dissipative_eigenvalue_s,
        template.passed,
    )
end

function _emt_frequency_invariant_unforced_scan_point(
    template::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
    frequency_hz::Float64,
)
    node_count = length(template.node_voltage_phasors)
    return EMTInitializationFrequencyPoint(
        template.formulation,
        :uniform_scan,
        frequency_hz,
        _emt_reactive_angular_frequency(
            request.formulation,
            2.0 * pi * frequency_hz,
        ),
        fill(frequency_hz, node_count),
        template.node_frequency_source_row_indices,
        template.source_frequency_successor_indices,
        template.frequency_subnetwork_count,
        zeros(ComplexF64, node_count),
        zeros(ComplexF64, node_count),
        zeros(ComplexF64, node_count),
        template.topology,
        template.admittance_symmetry_max_abs_error,
        template.minimum_dissipative_eigenvalue_s,
        template.passed,
    )
end

function _validate_emt_initialization_source_frequency(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    deck_source_frequencies = Float64[
        abs(Float64(row.sfreq)) / (2.0 * pi)
        for row in DeckParser.deck_over5a_source_rows(parsed)
        if abs(Int(row.iform)) == 14 &&
           (Float64(row.tstart) == 5432.0 || Float64(row.tstart) < 0.0)
    ]
    for element in parsed.elements
        element isa Union{TheveninSource,CurrentInjection} || continue
        element.value isa SinusoidalSourceSignal || continue
        signal = element.value
        if signal.frequency > 0.0 && abs(signal.offset) >
           64.0 * eps(Float64) * max(abs(signal.amplitude), 1.0)
            throw(_EMTInitializationRefusal(
                :unsupported_source_spectrum,
                :source_state,
                :frequency_hz,
                "a named source with simultaneous DC offset and AC amplitude requires a multifrequency initializer",
                (
                    source_frequency_hz=signal.frequency,
                    offset=signal.offset,
                    amplitude=signal.amplitude,
                ),
            ))
        end
        signal.frequency > 0.0 && abs(signal.amplitude) > 0.0 &&
            isapprox(
                signal.frequency,
                request.frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            ) || throw(_EMTInitializationRefusal(
                :source_frequency_mismatch,
                :source_state,
                :frequency_hz,
                "a named sinusoidal source does not match the primary initialization frequency",
                (
                    requested_frequency_hz=request.frequency_hz,
                    source_frequency_hz=signal.frequency,
                ),
            ))
    end
    isempty(deck_source_frequencies) && return request
    unique_frequencies = sort!(unique(deck_source_frequencies))
    any(
        frequency -> isapprox(
            frequency,
            request.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ),
        unique_frequencies,
    ) || throw(_EMTInitializationRefusal(
        :source_frequency_mismatch,
        :source_state,
        :frequency_hz,
        "the primary initialization frequency must match one declared source subnetwork",
        (
            requested_frequency_hz=request.frequency_hz,
            source_frequencies_hz=unique_frequencies,
        ),
    ))
    return request
end

function _emt_initialization_primary_point(
    points::Vector{EMTInitializationFrequencyPoint},
    frequency_hz::Float64,
)
    index = findfirst(
        point -> point.frequency_assignment === :initial_operating_point &&
            point.physical_frequency_hz == frequency_hz,
        points,
    )
    index === nothing && error("initialization scan omitted its primary frequency")
    return points[index]
end

function _validate_emt_operating_point_signatures(
    operating_point::EMTOperatingPoint,
    request::EMTInitializationRequest,
)
    fields = (
        (:project_signature, operating_point.project_signature,
            request.project_signature),
        (:settings_signature, operating_point.settings_signature,
            request.settings_signature),
        (:model_signature, operating_point.model_signature,
            request.model_signature),
    )
    for (field, source, target) in fields
        source == target || throw(_EMTInitializationRefusal(
            :stale_signature,
            :operating_point_mapping,
            field,
            "operating-point $field does not match the initialization request",
            (source_signature=source, target_signature=target),
        ))
    end
    isapprox(
        operating_point.frequency_hz,
        request.frequency_hz;
        atol=0.0,
        rtol=64.0 * eps(Float64),
    ) || throw(_EMTInitializationRefusal(
        :frequency_mismatch,
        :operating_point_mapping,
        :frequency_hz,
        "operating-point frequency does not match the initialization frequency",
        (
            source_frequency_hz=operating_point.frequency_hz,
            target_frequency_hz=request.frequency_hz,
        ),
    ))
    operating_point.time_origin_s == request.time_origin_s || throw(
        _EMTInitializationRefusal(
            :time_origin_mismatch,
            :operating_point_mapping,
            :time_origin_s,
            "operating-point time origin does not match the initialization request",
            (
                source_time_origin_s=operating_point.time_origin_s,
                target_time_origin_s=request.time_origin_s,
            ),
        ),
    )
    return operating_point
end

function _emt_operating_point_mappings(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    point::EMTInitializationFrequencyPoint,
)
    operating_point = request.operating_point
    operating_point === nothing && return OperatingPointMappingRecord[]
    operating_point isa EMTOperatingPoint || throw(_EMTInitializationRefusal(
        :unsupported_schema,
        :operating_point_mapping,
        :schema,
        "operating_point must be an EMTOperatingPoint",
        (runtime_type=string(typeof(operating_point)),),
    ))
    _validate_emt_operating_point_signatures(operating_point, request)
    mappings = OperatingPointMappingRecord[]
    for quantity in operating_point.quantities
        target_value = _emt_operating_point_quantity_target(
            quantity,
            operating_point,
        )
        node_index = get(parsed.node_map, quantity.asset, 0)
        node_index > 0 || throw(_EMTInitializationRefusal(
            :missing_asset,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point asset $(quantity.asset) is not a network node",
            (asset=quantity.asset,),
        ))
        required_value = point.node_voltage_phasors[node_index]
        uncertainty = quantity.absolute_uncertainty * quantity.scale_to_si
        residual = abs(target_value - required_value)
        allowance = request.tolerances.voltage_absolute_v +
            request.tolerances.voltage_relative * abs(required_value) +
            uncertainty
        passed = isfinite(residual) && residual <= allowance
        switch_group_index = findfirst(
            group -> node_index in group,
            point.topology.switch_node_groups,
        )
        reaction_nodes = switch_group_index === nothing ? Int[node_index] :
            point.topology.switch_node_groups[switch_group_index]
        constraint_current = sum(
            point.operating_constraint_current_phasors[reaction_nodes];
            init=0.0 + 0.0im,
        )
        push!(
            mappings,
            OperatingPointMappingRecord(
                quantity.asset,
                quantity.quantity,
                quantity.phase,
                quantity.value,
                target_value,
                quantity.unit,
                "V",
                quantity.basis,
                quantity.orientation,
                quantity.scale_to_si,
                quantity.orientation_sign,
                uncertainty,
                residual,
                constraint_current,
                passed,
            ),
        )
        passed || throw(_EMTInitializationRefusal(
            :infeasible_operating_target,
            :operating_point_mapping,
            quantity.quantity,
            "operating-point quantity $(quantity.asset)/$(quantity.quantity) " *
            "does not close the EMT equilibrium",
            (asset=quantity.asset, residual=residual, allowance=allowance),
        ))
    end
    return mappings
end

function _emt_initial_voltage_sample(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    point::EMTInitializationFrequencyPoint,
)
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    values = real.(point.node_voltage_phasors)
    return (
        source=:consistent_emt_initialization,
        outcome=:steady_state_initial_voltage_sample,
        harmonic_formulation=point.formulation,
        exact_discrete_histories=
            request.formulation isa TimestepMatchedFormulation,
        steady_state_frequency_hz=point.physical_frequency_hz,
        reactive_angular_frequency_rad_s=
            point.reactive_angular_frequency_rad_s,
        time_origin_s=request.time_origin_s,
        timestep_s=request.formulation isa TimestepMatchedFormulation ?
            request.formulation.timestep_s :
            DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        node_steady_state_frequencies_hz=
            copy(point.node_physical_frequencies_hz),
        node_frequency_source_row_indices=
            copy(point.node_frequency_source_row_indices),
        source_frequency_successor_indices=
            copy(point.source_frequency_successor_indices),
        steady_state_frequency_subnetwork_count=
            point.frequency_subnetwork_count,
        node_map=Dict{Symbol,Int}(parsed.node_map),
        node_names=ordered_node_names(parsed.node_map),
        node_voltage_phasors=copy(point.node_voltage_phasors),
        node_voltage_values=values,
        output_node_indices=output_node_indices,
        output_voltage_values=Float64[
            node == 0 ? 0.0 : values[node]
            for node in output_node_indices
        ],
        source_row_count=length(DeckParser.deck_over5a_source_rows(parsed)),
        distributed_transposed_line_count=length(
            DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed),
        ),
        saturated_transformer_branch_count=0,
        saturated_transformer_ideal_branch_count=0,
        saturated_transformer_linearized_nonlinear_branch_count=0,
    )
end

function _shift_emt_initialization_time_origin!(
    prepared::PreparedEMTStudy,
    time_origin_s::Float64,
)
    time_origin_s == 0.0 && return prepared
    runtime = prepared.runtime_template
    plan = runtime.plan
    for source_index in eachindex(plan.source_iform_values)
        abs(plan.source_iform_values[source_index]) == 14 || continue
        plan.source_time1_values[source_index] +=
            plan.source_sfreq_values[source_index] * time_origin_s
    end
    for signal in runtime.context.analytic_source_signals
        signal.source_type == 14 || continue
        signal.time1_s += signal.angular_frequency_or_rate * time_origin_s
    end
    for element in runtime.context.system.elements
        element isa Union{TheveninSource,CurrentInjection} || continue
        element.value isa SinusoidalSourceSignal || continue
        signal = element.value
        signal.phase += 2.0 * pi * signal.frequency * time_origin_s
    end
    return prepared
end

function _append_emt_initialization_state!(
    records::Vector{EMTInitializationStateRecord},
    owner::Symbol,
    state_family::Symbol,
    instance_count::Integer,
    initialization_basis::Symbol,
)
    instance_count == 0 && return records
    push!(
        records,
        EMTInitializationStateRecord(
            owner,
            state_family,
            instance_count,
            initialization_basis,
        ),
    )
    return records
end

function _emt_initialization_state_inventory(prepared::PreparedEMTStudy)
    context = prepared.runtime_template.context
    parsed = prepared.parsed
    elements = context.system.elements
    records = EMTInitializationStateRecord[]
    _append_emt_initialization_state!(
        records,
        :network_voltage,
        :algebraic,
        context.system.node_count,
        :harmonic_network_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :network_topology,
        :discrete,
        1,
        :ranked_connected_component_classification,
    )
    source_element_count = count(
        element -> element isa Union{TheveninSource,CurrentInjection},
        elements,
    )
    source_program_count = context.source_function_runtime === nothing ? 0 :
        max(context.source_function_runtime.plan.source_row_count, 1)
    source_count = max(
        source_element_count,
        length(context.analytic_source_signals),
        source_program_count,
    )
    _append_emt_initialization_state!(
        records,
        :source_state,
        :algebraic,
        source_count,
        :phasor_time_mapping,
    )
    history_kinds = context.electromagnetic_history_plan.kinds
    step_configs = prepared.runtime_template.step_configs
    algebraic_named_network =
        _emt_frequency_invariant_named_network(parsed) &&
        isempty(history_kinds) &&
        context.control_system_runtime === nothing &&
        !(step_configs isa DynamicDeckStepConfigProvider &&
          step_configs.nonlinear_current_config !== nothing) &&
        context.deck_time_switch_count == 0 &&
        isempty(context.series_rlc_alterations)
    if algebraic_named_network
        _append_emt_initialization_state!(
            records,
            :energy_accumulator,
            :algebraic,
            length(context.branch_energy_values) +
                length(context.switch_energy_values),
            :time_zero_energy_reference,
        )
        _append_emt_initialization_state!(
            records,
            :output_cursor,
            :discrete,
            1,
            :time_zero_output_epoch,
        )
        _append_emt_initialization_state!(
            records,
            :checkpoint_continuation_state,
            :checkpoint,
            1,
            :complete_prepared_runtime_before_first_advance,
        )
        sort!(records; by=record -> (String(record.state_family), String(record.owner)))
        return records
    end
    for (kind, owner, basis) in (
        (SERIES_RL_HISTORY, :series_rl_history, :exact_companion_recurrence),
        (SERIES_RLC_HISTORY, :series_rlc_history, :exact_companion_recurrence),
        (CAPACITOR_HISTORY, :capacitor_history, :exact_companion_recurrence),
        (COUPLED_INDUCTIVE_HISTORY, :coupled_inductive_history, :exact_coupled_recurrence),
        (COUPLED_SERIES_RL_HISTORY, :coupled_series_rl_history, :exact_coupled_recurrence),
        (BREQIV_HISTORY, :breqiv_history, :frequency_domain_history_mapping),
    )
        _append_emt_initialization_state!(
            records,
            owner,
            :delayed_history,
            count(==(kind), history_kinds),
            basis,
        )
    end
    _append_emt_initialization_state!(
        records,
        :complex_modal_line_history,
        :delayed_history,
        count(element -> element isa ComplexModalBergeronLine, elements),
        :traveling_wave_harmonic_prehistory,
    )
    _append_emt_initialization_state!(
        records,
        :frequency_dependent_line_history,
        :delayed_history,
        count(element -> element isa Union{
            SemlyenFrequencyDependentLine,
            SampledFrequencyDependentLine,
            SampledFrequencyDependentLineGroup,
        }, elements),
        :recursive_convolution_harmonic_prehistory,
    )
    _append_emt_initialization_state!(
        records,
        :distributed_line_history,
        :delayed_history,
        length(DeckParser.deck_distributed_transposed_line_modal_branch_states(parsed)),
        :modal_traveling_wave_harmonic_prehistory,
    )
    control_runtime = context.control_system_runtime
    if control_runtime !== nothing
        control_state_count = length(control_runtime.state.values) +
            length(control_runtime.state.function_states)
        _append_emt_initialization_state!(
            records,
            :control_state,
            :discrete,
            max(control_state_count, 1),
            :control_network_steady_state,
        )
        _append_emt_initialization_state!(
            records,
            :control_frequency_history,
            :delayed_history,
            length(control_runtime.frequency_initializations),
            :control_frequency_response,
        )
    end
    if step_configs isa DynamicDeckStepConfigProvider &&
       step_configs.nonlinear_current_config !== nothing
        nonlinear_types = Int.(
            step_configs.nonlinear_current_config.nonlinear_types,
        )
        saturated_transformer_count = get(
            step_configs.nonlinear_current_config,
            :saturated_transformer_residual_flux_initialized,
            false,
        ) ? length(get(
            step_configs.nonlinear_current_config,
            :saturated_transformer_internal_top_node_indices,
            Int[],
        )) : 0
        _append_emt_initialization_state!(
            records,
            :pseudo_nonlinear_inductor_state,
            :continuous,
            count(==(PSEUDO_NONLINEAR_INDUCTOR_TYPE), nonlinear_types),
            :harmonic_flux_companion_prehistory,
        )
        _append_emt_initialization_state!(
            records,
            :saturated_transformer_magnetic_state,
            :continuous,
            saturated_transformer_count,
            :declared_flux_characteristic_equilibrium,
        )
        _append_emt_initialization_state!(
            records,
            :piecewise_nonlinear_inductor_state,
            :continuous,
            count(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types),
            :harmonic_characteristic_equilibrium,
        )
        _append_emt_initialization_state!(
            records,
            :hysteretic_magnetic_state,
            :continuous,
            count(==(HYSTERETIC_INDUCTOR_NONLINEAR_TYPE), nonlinear_types),
            :model_owned_hysteresis_equilibrium,
        )
    end
    switch_count = count(element -> element isa Union{
        IdealSwitch,
        TimeSwitch,
        CurrentZeroSwitch,
        TACSControlledSwitch,
    }, elements)
    _append_emt_initialization_state!(
        records,
        :switch_mode,
        :discrete,
        switch_count,
        :declared_initial_topology,
    )
    _append_emt_initialization_state!(
        records,
        :switch_event_state,
        :scheduler,
        context.deck_time_switch_count,
        :initial_event_surface_classification,
    )
    _append_emt_initialization_state!(
        records,
        :energy_accumulator,
        :algebraic,
        length(context.branch_energy_values) + length(context.switch_energy_values),
        :time_zero_energy_reference,
    )
    _append_emt_initialization_state!(
        records,
        :series_rlc_alteration_schedule,
        :scheduler,
        length(context.series_rlc_alterations),
        :ordered_future_event_cursor,
    )
    _append_emt_initialization_state!(
        records,
        :output_cursor,
        :discrete,
        1,
        :time_zero_output_epoch,
    )
    _append_emt_initialization_state!(
        records,
        :checkpoint_continuation_state,
        :checkpoint,
        1,
        :complete_prepared_runtime_before_first_advance,
    )
    sort!(records; by=record -> (String(record.state_family), String(record.owner)))
    return records
end

function _emt_initialized_state_owners(
    records::Vector{EMTInitializationStateRecord},
)
    return sort!(unique(record.owner for record in records); by=String)
end

function _emt_unsupported_initialization_owners(
    parsed::DeckParser.DeckParseResult,
)
    elements = parsed.elements
    unsupported = Symbol[]
    any(element -> element isa PowerSemiconductorSwitch ||
        element isa PowerSemiconductorBridgeLeg ||
        element isa PowerSemiconductorBridgeTopology, elements) &&
        push!(unsupported, :switch_detailed_periodic_prehistory)
    isempty(DeckParser.deck_universal_machine_definition_rows(parsed)) ||
        push!(unsupported, :universal_machine_state)
    isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) ||
        push!(unsupported, :synchronous_machine_state)
    nonlinear_periodic_count = sum(length, (
        DeckParser.deck_zinc_oxide_nonlinear_rows(parsed),
        DeckParser.deck_nonlinear_resistance_rows(parsed),
        DeckParser.deck_triggered_timed_resistance_rows(parsed),
        DeckParser.deck_switching_nonlinear_resistor_rows(parsed),
        DeckParser.deck_arrester_nonlinear_rows(parsed),
    ))
    nonlinear_periodic_count == 0 ||
        push!(unsupported, :nonlinear_periodic_state)
    any(
        row -> row.condition_kind == :node_voltage_initial_condition,
        DeckParser.deck_node_initial_condition_rows(parsed),
    ) && push!(unsupported, :declared_node_initial_condition)
    any(
        row -> row.random_opening_standard_deviation_s > 0.0,
        DeckParser.deck_over5_switch_rows(parsed),
    ) && push!(unsupported, :stochastic_switch_state)
    return sort!(unique(unsupported); by=String)
end

function _emt_initialization_probe_runtime(prepared::PreparedEMTStudy)
    runtime = deepcopy(prepared.runtime_template)
    context = runtime.context
    context.step_count >= 1 || throw(_EMTInitializationRefusal(
        :missing_probe_horizon,
        :no_artificial_transient,
        :time_window,
        "initialization requires at least one timestep for its coupled transient probe",
        (step_count=context.step_count,),
    ))
    context.t_end_s = context.dt_s
    context.step_count = 1
    context.time_s = Vector{Float64}(undef, 2)
    context.voltage_pu = Matrix{Float64}(undef, context.system.node_count, 2)
    context.output_pu = Matrix{Float64}(
        undef,
        length(context.output_channel_names),
        2,
    )
    context.recorded_step_indices = Int[0, 1]
    context.trace_write_index = 1
    return runtime
end

function _emt_initialization_probe_metric(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    runtime = _emt_initialization_probe_runtime(prepared)
    run = _run_prepared_dynamic_deck!(runtime; collect_run_diagnostics=true)
    trace = run.trace
    size(trace.voltage_pu, 2) >= 2 || throw(_EMTInitializationRefusal(
        :incomplete_probe,
        :no_artificial_transient,
        :node_voltage,
        "initialization transient probe did not produce two voltage samples",
        (sample_count=size(trace.voltage_pu, 2),),
    ))
    timestep_s = runtime.context.dt_s
    rotations = ComplexF64[
        cis(2.0 * pi * frequency_hz * timestep_s)
        for frequency_hz in point.node_physical_frequencies_hz
    ]
    expected = real.(point.node_voltage_phasors .* rotations)
    actual = Float64.(trace.voltage_pu[:, 2])
    error = actual - expected
    reference_scale = max(norm(expected), request.tolerances.voltage_absolute_v)
    raw_normalized_rms = norm(error) / reference_scale
    raw_scaled_discontinuity = maximum(abs, error; init=0.0) /
        max(maximum(abs, expected; init=0.0), request.tolerances.voltage_absolute_v)
    physical_to_trapezoidal_warping = if request.formulation isa
                                         PhysicalFrequencyFormulation
        maximum(point.node_physical_frequencies_hz; init=0.0) do frequency_hz
            frequency_step = pi * frequency_hz * timestep_s
            frequency_step == 0.0 ? 0.0 :
                abs(tan(frequency_step) / frequency_step - 1.0)
        end
    else
        0.0
    end
    physical_frequency_probe_bound =
        2.0 * physical_to_trapezoidal_warping /
        max(1.0 - physical_to_trapezoidal_warping, eps(Float64))
    maximum_physical_step_angle = 2.0 * pi *
        maximum(point.node_physical_frequencies_hz; init=0.0) * timestep_s
    physical_frequency_probe_bound = max(
        physical_frequency_probe_bound,
        0.5 * maximum_physical_step_angle^2,
    )
    normalized_rms = max(
        0.0,
        raw_normalized_rms - physical_frequency_probe_bound,
    )
    scaled_discontinuity = max(
        0.0,
        raw_scaled_discontinuity - physical_frequency_probe_bound,
    )
    threshold = request.tolerances.no_artificial_transient_normalized_rms
    passed = isfinite(normalized_rms) && isfinite(scaled_discontinuity) &&
        normalized_rms <= threshold &&
        scaled_discontinuity <=
            request.tolerances.first_step_scaled_discontinuity
    return NoArtificialTransientMetric(
        :node_voltage,
        "V",
        request.time_origin_s,
        request.time_origin_s + timestep_s,
        normalized_rms,
        scaled_discontinuity,
        max(raw_normalized_rms, raw_scaled_discontinuity),
        0.0,
        threshold,
        passed,
    )
end

function _emt_initialization_residuals(
    point::EMTInitializationFrequencyPoint,
    mappings::Vector{OperatingPointMappingRecord},
    request::EMTInitializationRequest,
    fixed_source_load_flow::Union{Nothing,FixedSourceLoadFlowResult}=nothing,
)
    tolerances = request.tolerances
    current_scale = max(
        norm(point.source_injection_phasors, Inf),
        tolerances.current_absolute_a,
    )
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    residuals = EMTInitializationResidual[
        EMTInitializationResidual(
            :network_equilibrium,
            :nodal_kcl,
            "A",
            point.topology.maximum_residual_a,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            point.topology.maximum_residual_a / current_allowance,
            point.topology.maximum_residual_a <= current_allowance,
        ),
    ]
    for mapping in mappings
        allowance = tolerances.voltage_absolute_v +
            tolerances.voltage_relative * abs(mapping.target_value_si) +
            mapping.absolute_uncertainty_si
        push!(
            residuals,
            EMTInitializationResidual(
                :operating_point_mapping,
                mapping.quantity,
                mapping.target_unit,
                mapping.residual,
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(mapping.target_value_si),
                mapping.absolute_uncertainty_si,
                mapping.residual / allowance,
                mapping.passed,
            ),
        )
    end
    if fixed_source_load_flow !== nothing
        for constraint_index in eachindex(
            fixed_source_load_flow.constraint_kinds,
        )
            for (quantity, unit, target, actual) in (
                (
                    :active_power,
                    "W",
                    fixed_source_load_flow.constraint_target_active_powers[
                        constraint_index
                    ],
                    fixed_source_load_flow.constraint_active_powers[
                        constraint_index
                    ],
                ),
                (
                    :reactive_power,
                    "var",
                    fixed_source_load_flow.constraint_target_reactive_powers[
                        constraint_index
                    ],
                    fixed_source_load_flow.constraint_reactive_powers[
                        constraint_index
                    ],
                ),
            )
                target === missing && continue
                reference_scale = max(abs(Float64(target)), 1.0)
                allowance = tolerances.power_absolute_w +
                    tolerances.power_relative * reference_scale
                residual = abs(Float64(actual) - Float64(target))
                push!(
                    residuals,
                    EMTInitializationResidual(
                        :fixed_source_operating_point,
                        quantity,
                        unit,
                        residual,
                        tolerances.power_absolute_w,
                        tolerances.power_relative,
                        reference_scale,
                        0.0,
                        residual / allowance,
                        residual <= allowance,
                    ),
                )
            end
        end
    end
    return residuals
end

function _emt_pseudo_nonlinear_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    nonlinear_indices = findall(
        ==(PSEUDO_NONLINEAR_INDUCTOR_TYPE),
        nonlinear_types,
    )
    isempty(nonlinear_indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    stored_fluxes = Float64.(config.initial_stored_voltage_values)
    companion_currents = Float64.(config.initial_companion_current_values)
    table_indices = Int.(config.initial_table_index_values)
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_residuals = Float64[]
    flux_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    for index in nonlinear_indices
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_node) -
            _node_voltage_phasor(point.node_voltage_phasors, to_node)
        endpoint_nodes = filter(!=(0), (from_node, to_node))
        frequency_hz = point.node_physical_frequencies_hz[first(endpoint_nodes)]
        all(
            node -> point.node_physical_frequencies_hz[node] == frequency_hz,
            endpoint_nodes,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :pseudo_nonlinear_inductor_state,
            :frequency_hz,
            "pseudo-nonlinear inductor endpoints belong to different frequency subnetworks",
            (owner_index=index, from_node, to_node),
        ))
        angular_frequency = 2.0 * pi * frequency_hz
        angular_frequency > 0.0 || throw(_EMTInitializationRefusal(
            :unsupported_dc_state,
            :pseudo_nonlinear_inductor_state,
            :flux,
            "pseudo-nonlinear inductor prehistory requires a positive harmonic frequency",
            (owner_index=index, frequency_hz),
        ))
        expected_flux = imag(branch_phasor) / angular_frequency -
            delta2 * real(branch_phasor)
        actual_flux = stored_fluxes[index]
        push!(flux_residuals, abs(actual_flux - expected_flux))
        push!(flux_scales, abs(expected_flux))
        table_index = table_indices[index]
        1 <= table_index <= length(slopes) || throw(
            _EMTInitializationRefusal(
                :invalid_characteristic_state,
                :pseudo_nonlinear_inductor_state,
                :table_index,
                "pseudo-nonlinear inductor initialization selected an invalid characteristic segment",
                (owner_index=index, table_index),
            ),
        )
        expected_current = expected_flux * slopes[table_index] / delta2
        actual_current = companion_currents[index]
        push!(current_residuals, abs(actual_current - expected_current))
        push!(current_scales, abs(expected_current))
    end
    tolerances = request.tolerances
    flux_residual = maximum(flux_residuals)
    flux_scale = maximum(flux_scales; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_residual = maximum(current_residuals)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :pseudo_nonlinear_inductor_state,
            :flux_history_recurrence,
            "Wb",
            flux_residual,
            tolerances.flux_absolute_wb,
            0.0,
            flux_scale,
            0.0,
            flux_residual / flux_allowance,
            flux_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :pseudo_nonlinear_inductor_state,
            :companion_current_recurrence,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
    ]
end
function _emt_saturated_transformer_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    get(
        config,
        :saturated_transformer_residual_flux_initialized,
        false,
    ) || return EMTInitializationResidual[]
    count = length(get(
        config,
        :saturated_transformer_internal_top_node_indices,
        Int[],
    ))
    count > 0 || return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_from_nodes[1:count])
    to_nodes = Int.(config.nonlinear_to_nodes[1:count])
    declared_currents = Float64.(
        config.nonlinear_steady_state_current_values[1:count],
    )
    declared_fluxes = Float64.(
        config.nonlinear_steady_state_flux_values[1:count],
    )
    stored_fluxes = Float64.(config.initial_stored_voltage_values[1:count])
    companion_currents = Float64.(
        config.initial_companion_current_values[1:count],
    )
    table_indices = Int.(config.initial_table_index_values[1:count])
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_errors = Float64[]
    current_errors = Float64[]
    current_scales = Float64[]
    residual_fluxes = Float64[]
    for index in 1:count
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_nodes[index]) -
            _node_voltage_phasor(point.node_voltage_phasors, to_nodes[index])
        branch_voltage = real(branch_phasor)
        reconstructed_flux = stored_fluxes[index] + delta2 * branch_voltage
        push!(flux_errors, abs(reconstructed_flux - declared_fluxes[index]))
        table_index = table_indices[index]
        1 <= table_index <= length(slopes) || throw(
            _EMTInitializationRefusal(
                :invalid_characteristic_state,
                :saturated_transformer_magnetic_state,
                :table_index,
                "saturated-transformer initialization selected an invalid characteristic segment",
                (transformer_index=index, table_index),
            ),
        )
        reconstructed_current = companion_currents[index] +
            slopes[table_index] * branch_voltage
        push!(current_errors, abs(reconstructed_current - declared_currents[index]))
        push!(current_scales, abs(declared_currents[index]))
        harmonic_flux = point.physical_frequency_hz == 0.0 ? 0.0 :
            imag(branch_phasor) / (2.0 * pi * point.physical_frequency_hz)
        push!(residual_fluxes, declared_fluxes[index] - harmonic_flux)
    end
    tolerances = request.tolerances
    flux_error = maximum(flux_errors; init=0.0)
    current_error = maximum(current_errors; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    maximum_residual_flux = maximum(abs, residual_fluxes; init=0.0)
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :declared_time_zero_flux,
            "Wb-turn",
            flux_error,
            tolerances.flux_absolute_wb,
            0.0,
            maximum(abs, declared_fluxes; init=0.0),
            0.0,
            flux_error / tolerances.flux_absolute_wb,
            flux_error <= tolerances.flux_absolute_wb,
        ),
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :characteristic_current_equilibrium,
            "A",
            current_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_error / current_allowance,
            current_error <= current_allowance,
        ),
        EMTInitializationResidual(
            :saturated_transformer_magnetic_state,
            :residual_flux_magnitude,
            "Wb-turn",
            maximum_residual_flux,
            max(maximum_residual_flux, tolerances.flux_absolute_wb),
            0.0,
            maximum_residual_flux,
            0.0,
            maximum_residual_flux == 0.0 ? 0.0 : 1.0,
            true,
        ),
    ]
end

function _emt_piecewise_nonlinear_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    indices = findall(==(PIECEWISE_NONLINEAR_INDUCTOR_TYPE), nonlinear_types)
    isempty(indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    table_starts = Int.(config.nonlinear_admittance_nodes)
    table_ends = Int.(config.nonlinear_table_end_indices)
    currents_a = Float64.(config.cchar)
    fluxes_wb = Float64.(config.vchar)
    initialized_fluxes = Float64.(
        config.piecewise_nonlinear_inductor_initial_flux_values,
    )
    predictor_fluxes = Float64.(
        config.piecewise_nonlinear_inductor_initial_predictor_flux_values,
    )
    initialized_currents = Float64.(
        config.piecewise_nonlinear_inductor_initial_current_values,
    )
    initialized_segments = Int.(
        config.piecewise_nonlinear_inductor_initial_segment_values,
    )
    delta2 = Float64(config.delta2)
    characteristic_residuals = Float64[]
    characteristic_scales = Float64[]
    predictor_residuals = Float64[]
    predictor_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    segment_residuals = Float64[]
    for index in indices
        from_node = from_nodes[index]
        to_node = to_nodes[index]
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_node) -
            _node_voltage_phasor(point.node_voltage_phasors, to_node)
        frequency_hz = _nonlinear_initial_frequency_hz(
            prepared.runtime_template.steady_state_initial_sample,
            from_node,
            to_node,
        )
        reactive_angular_frequency = _steady_state_reactive_angular_frequency(
            prepared.runtime_template.steady_state_initial_sample,
            frequency_hz,
        )
        expected_flux_wb = imag(branch_phasor) / reactive_angular_frequency
        expected_predictor_flux_wb = expected_flux_wb + delta2 * real(branch_phasor)
        state = _piecewise_nonlinear_inductor_characteristic_state(
            initialized_currents[index],
            initialized_fluxes[index],
            table_starts[index],
            table_ends[index],
            currents_a,
            fluxes_wb;
            flux_tolerance_wb=Float64(config.flzero),
        )
        push!(
            characteristic_residuals,
            abs(initialized_fluxes[index] - state.flux_wb),
        )
        push!(characteristic_scales, abs(state.flux_wb))
        push!(
            predictor_residuals,
            abs(predictor_fluxes[index] - expected_predictor_flux_wb),
        )
        push!(predictor_scales, abs(expected_predictor_flux_wb))
        expected_current_a = begin
            declared_current_a = Float64(
                config.nonlinear_steady_state_current_values[index],
            )
            declared_flux_wb = Float64(
                config.nonlinear_steady_state_flux_values[index],
            )
            if declared_current_a == 0.0 && declared_flux_wb == 0.0
                0.0
            else
                secant_inductance_h = declared_flux_wb / declared_current_a
                real(
                    branch_phasor /
                    complex(0.0, reactive_angular_frequency * secant_inductance_h),
                )
            end
        end
        push!(current_residuals, abs(initialized_currents[index] - expected_current_a))
        push!(current_scales, abs(expected_current_a))
        push!(segment_residuals, abs(initialized_segments[index] - state.segment))
    end
    tolerances = request.tolerances
    characteristic_residual = maximum(characteristic_residuals; init=0.0)
    characteristic_scale = maximum(characteristic_scales; init=0.0)
    predictor_residual = maximum(predictor_residuals; init=0.0)
    predictor_scale = maximum(predictor_scales; init=0.0)
    current_residual = maximum(current_residuals; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    segment_residual = maximum(segment_residuals; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :characteristic_flux_equilibrium,
            "Wb",
            characteristic_residual,
            tolerances.flux_absolute_wb,
            0.0,
            characteristic_scale,
            0.0,
            characteristic_residual / flux_allowance,
            characteristic_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :flux_half_step_recurrence,
            "Wb",
            predictor_residual,
            tolerances.flux_absolute_wb,
            0.0,
            predictor_scale,
            0.0,
            predictor_residual / flux_allowance,
            predictor_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :harmonic_current_equilibrium,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
        EMTInitializationResidual(
            :piecewise_nonlinear_inductor_state,
            :active_characteristic_segment,
            "index",
            segment_residual,
            0.0,
            0.0,
            1.0,
            0.0,
            segment_residual == 0.0 ? 0.0 : Inf,
            segment_residual == 0.0,
        ),
    ]
end

function _emt_hysteretic_initialization_residuals(
    prepared::PreparedEMTStudy,
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    step_configs = prepared.runtime_template.step_configs
    step_configs isa DynamicDeckStepConfigProvider ||
        return EMTInitializationResidual[]
    config = step_configs.nonlinear_current_config
    config === nothing && return EMTInitializationResidual[]
    nonlinear_types = Int.(config.nonlinear_types)
    indices = findall(==(HYSTERETIC_INDUCTOR_NONLINEAR_TYPE), nonlinear_types)
    isempty(indices) && return EMTInitializationResidual[]
    from_nodes = Int.(config.nonlinear_deck_from_nodes)
    to_nodes = Int.(config.nonlinear_deck_to_nodes)
    state_starts = Int.(config.nonlinear_admittance_nodes)
    time_zero_fluxes = Float64.(config.hysteretic_initial_flux_values)
    runtime_fluxes = Float64.(config.initial_runtime_voltage_values)
    currents = Float64.(config.hysteretic_initial_current_values)
    companion_currents = Float64.(config.initial_companion_current_values)
    slopes = Float64.(config.gslope)
    delta2 = Float64(config.delta2)
    flux_residuals = Float64[]
    flux_scales = Float64[]
    current_residuals = Float64[]
    current_scales = Float64[]
    for index in indices
        branch_phasor =
            _node_voltage_phasor(point.node_voltage_phasors, from_nodes[index]) -
            _node_voltage_phasor(point.node_voltage_phasors, to_nodes[index])
        branch_voltage_v = real(branch_phasor)
        reconstructed_flux_wb =
            runtime_fluxes[index] + delta2 * branch_voltage_v
        push!(
            flux_residuals,
            abs(reconstructed_flux_wb - time_zero_fluxes[index]),
        )
        push!(flux_scales, abs(time_zero_fluxes[index]))
        companion_admittance_s = slopes[state_starts[index] + 1]
        reconstructed_current_a =
            companion_currents[index] + companion_admittance_s * branch_voltage_v
        push!(
            current_residuals,
            abs(reconstructed_current_a - currents[index]),
        )
        push!(current_scales, abs(currents[index]))
    end
    tolerances = request.tolerances
    flux_residual = maximum(flux_residuals; init=0.0)
    flux_scale = maximum(flux_scales; init=0.0)
    flux_allowance = tolerances.flux_absolute_wb
    current_residual = maximum(current_residuals; init=0.0)
    current_scale = maximum(current_scales; init=0.0)
    current_allowance = tolerances.current_absolute_a +
        tolerances.current_relative * current_scale
    return EMTInitializationResidual[
        EMTInitializationResidual(
            :hysteretic_magnetic_state,
            :flux_half_step_recurrence,
            "Wb",
            flux_residual,
            tolerances.flux_absolute_wb,
            0.0,
            flux_scale,
            0.0,
            flux_residual / flux_allowance,
            flux_residual <= flux_allowance,
        ),
        EMTInitializationResidual(
            :hysteretic_magnetic_state,
            :companion_current_equilibrium,
            "A",
            current_residual,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
            0.0,
            current_residual / current_allowance,
            current_residual <= current_allowance,
        ),
    ]
end

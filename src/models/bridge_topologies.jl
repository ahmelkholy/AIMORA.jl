module BridgeTopologies

using SHA
using ..StudyCore: ParameterProvenance, PhysicalModelParameter

export BridgeNode,
       BridgeValvePosition,
       BridgePassivePosition,
       BridgeStateGroup,
       BridgeTopologyDescriptor,
       generic_bridge_topology_provenance,
       bridge_topology_signature,
       bridge_topology_incidence,
       bridge_topology_state_is_allowed,
       single_phase_graetz_topology,
       polyphase_bridge_topology,
       multigroup_bridge_topology,
       two_level_bridge_topology,
       full_bridge_topology,
       step_down_chopper_topology,
       step_up_chopper_topology,
       inverting_buck_boost_topology,
       bidirectional_chopper_topology,
       neutral_point_clamped_leg_topology,
       t_type_leg_topology,
       flying_capacitor_leg_topology,
       cascaded_h_bridge_phase_topology,
       matrix_converter_topology,
       cycloconverter_topology

const BRIDGE_TOPOLOGY_SCHEMA = :aimora_bridge_topology_v1
const BRIDGE_FAMILIES = (
    :single_phase_graetz,
    :polyphase_bridge,
    :multi_group_bridge,
    :two_level_bridge,
    :full_bridge,
    :step_down_chopper,
    :step_up_chopper,
    :inverting_buck_boost,
    :bidirectional_chopper,
    :neutral_point_clamped_leg,
    :t_type_leg,
    :flying_capacitor_leg,
    :cascaded_h_bridge_phase,
    :matrix_converter,
    :cycloconverter,
)
const BRIDGE_VALVE_CLASSES = (:diode, :thyristor, :self_commutated)
const BRIDGE_PASSIVE_KINDS = (:conductance, :series_rl, :series_rlc, :capacitor)
const BRIDGE_STATE_GROUP_KINDS = (
    :natural_conduction,
    :complementary_leg,
    :admissible_gate_state,
)

"""One named physical node in a bridge topology contract."""
struct BridgeNode
    name::Symbol
    node::Int
    role::Symbol

    function BridgeNode(name::Symbol, node::Integer, role::Symbol)
        isempty(String(name)) && throw(ArgumentError("bridge node name must not be empty"))
        Int(node) >= 0 || throw(ArgumentError("bridge node index must be nonnegative"))
        role in (:external, :internal) || throw(ArgumentError(
            "bridge node role must be :external or :internal",
        ))
        return new(name, Int(node), role)
    end
end

"""One ordered semiconductor position with positive branch orientation from `from_node` to `to_node`."""
struct BridgeValvePosition
    name::Symbol
    from_node::Int
    to_node::Int
    valve_class::Symbol
    group::Int
    cell::Int

    function BridgeValvePosition(
        name::Symbol,
        from_node::Integer,
        to_node::Integer,
        valve_class::Symbol;
        group::Integer=1,
        cell::Integer=0,
    )
        isempty(String(name)) && throw(ArgumentError("bridge valve name must not be empty"))
        from = Int(from_node)
        to = Int(to_node)
        from >= 0 && to >= 0 && from != to || throw(ArgumentError(
            "bridge valve terminals must be distinct nonnegative nodes",
        ))
        valve_class in BRIDGE_VALVE_CLASSES || throw(ArgumentError(
            "bridge valve class must be :diode, :thyristor, or :self_commutated",
        ))
        Int(group) > 0 || throw(ArgumentError("bridge valve group must be positive"))
        Int(cell) >= 0 || throw(ArgumentError("bridge valve cell must be nonnegative"))
        return new(name, from, to, valve_class, Int(group), Int(cell))
    end
end

"""One ordered finite passive position composed through an existing branch owner."""
struct BridgePassivePosition
    name::Symbol
    from_node::Int
    to_node::Int
    kind::Symbol
    group::Int
    cell::Int

    function BridgePassivePosition(
        name::Symbol,
        from_node::Integer,
        to_node::Integer,
        kind::Symbol;
        group::Integer=1,
        cell::Integer=0,
    )
        isempty(String(name)) && throw(ArgumentError("bridge passive name must not be empty"))
        from = Int(from_node)
        to = Int(to_node)
        from >= 0 && to >= 0 && from != to || throw(ArgumentError(
            "bridge passive terminals must be distinct nonnegative nodes",
        ))
        kind in BRIDGE_PASSIVE_KINDS || throw(ArgumentError(
            "unsupported bridge passive kind $kind",
        ))
        Int(group) > 0 || throw(ArgumentError("bridge passive group must be positive"))
        Int(cell) >= 0 || throw(ArgumentError("bridge passive cell must be nonnegative"))
        return new(name, from, to, kind, Int(group), Int(cell))
    end
end

"""A local finite-state constraint over ordered valve positions; columns are admitted states."""
struct BridgeStateGroup
    name::Symbol
    kind::Symbol
    position_indices::Vector{Int}
    state_names::Vector{Symbol}
    admitted_states::BitMatrix

    function BridgeStateGroup(
        name::Symbol,
        kind::Symbol,
        position_indices::AbstractVector{<:Integer},
        state_names::AbstractVector{Symbol},
        admitted_states::AbstractMatrix{Bool},
    )
        isempty(String(name)) && throw(ArgumentError("bridge state-group name must not be empty"))
        kind in BRIDGE_STATE_GROUP_KINDS || throw(ArgumentError(
            "unsupported bridge state-group kind $kind",
        ))
        positions = Int.(position_indices)
        !isempty(positions) && all(>(0), positions) &&
            length(unique(positions)) == length(positions) || throw(ArgumentError(
            "bridge state-group positions must be unique positive indices",
        ))
        names = collect(state_names)
        !isempty(names) && length(unique(names)) == length(names) || throw(ArgumentError(
            "bridge state names must be nonempty and unique",
        ))
        size(admitted_states) == (length(positions), length(names)) ||
            throw(DimensionMismatch(
                "bridge admitted-state rows must match positions and columns must match names",
            ))
        return new(name, kind, positions, names, BitMatrix(admitted_states))
    end
end

"""Dependency-light, terminal-oriented switching-detailed topology contract."""
struct BridgeTopologyDescriptor
    schema::Symbol
    family::Symbol
    nodes::Vector{BridgeNode}
    valve_positions::Vector{BridgeValvePosition}
    passive_positions::Vector{BridgePassivePosition}
    state_groups::Vector{BridgeStateGroup}
    provenance::ParameterProvenance
    licence::String
    redistribution::Symbol
    topology_signature::String

    function BridgeTopologyDescriptor(
        family::Symbol,
        nodes::AbstractVector{BridgeNode},
        valve_positions::AbstractVector{BridgeValvePosition},
        passive_positions::AbstractVector{BridgePassivePosition},
        state_groups::AbstractVector{BridgeStateGroup};
        provenance::ParameterProvenance=generic_bridge_topology_provenance(family),
        licence::AbstractString="PolyForm-Noncommercial-1.0.0",
        redistribution::Symbol=:public,
    )
        family in BRIDGE_FAMILIES || throw(ArgumentError(
            "unsupported bridge topology family $family",
        ))
        node_values = collect(nodes)
        valves = collect(valve_positions)
        passives = collect(passive_positions)
        groups = collect(state_groups)
        !isempty(node_values) || throw(ArgumentError("bridge topology requires named nodes"))
        node_names = getfield.(node_values, :name)
        length(unique(node_names)) == length(node_names) || throw(ArgumentError(
            "bridge topology node names must be unique",
        ))
        external_nodes = [node.node for node in node_values if node.role === :external]
        !isempty(external_nodes) && length(unique(external_nodes)) == length(external_nodes) ||
            throw(ArgumentError("bridge external nodes must be nonempty and unique"))
        declared_nodes = Set(getfield.(node_values, :node))
        all(position -> position.from_node in declared_nodes &&
            position.to_node in declared_nodes, valves) || throw(ArgumentError(
            "every bridge valve terminal must be a declared topology node",
        ))
        all(position -> position.from_node in declared_nodes &&
            position.to_node in declared_nodes, passives) || throw(ArgumentError(
            "every bridge passive terminal must be a declared topology node",
        ))
        1 <= length(valves) <= 256 || throw(ArgumentError(
            "bridge topology must contain between one and 256 valve positions",
        ))
        length(passives) <= 128 || throw(ArgumentError(
            "bridge topology cannot exceed 128 passive positions",
        ))
        position_names = getfield.(valves, :name)
        length(unique(position_names)) == length(position_names) || throw(ArgumentError(
            "bridge valve-position names must be unique",
        ))
        passive_names = getfield.(passives, :name)
        length(unique(passive_names)) == length(passive_names) || throw(ArgumentError(
            "bridge passive-position names must be unique",
        ))
        group_names = getfield.(groups, :name)
        length(unique(group_names)) == length(group_names) || throw(ArgumentError(
            "bridge state-group names must be unique",
        ))
        all(group -> all(index -> index <= length(valves), group.position_indices), groups) ||
            throw(ArgumentError("bridge state group references an absent valve position"))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "bridge topology provenance must describe physical model parameters",
        ))
        isempty(strip(licence)) && throw(ArgumentError("bridge topology licence must not be empty"))
        redistribution in (:public, :private, :prohibited) || throw(ArgumentError(
            "bridge topology redistribution must be :public, :private, or :prohibited",
        ))
        signature = _bridge_topology_signature(
            family,
            node_values,
            valves,
            passives,
            groups,
            provenance,
            String(licence),
            redistribution,
        )
        return new(
            BRIDGE_TOPOLOGY_SCHEMA,
            family,
            node_values,
            valves,
            passives,
            groups,
            provenance,
            String(licence),
            redistribution,
            signature,
        )
    end
end

function generic_bridge_topology_provenance(family::Symbol)
    return ParameterProvenance(
        "AIMORA generic $(replace(String(family), '_' => ' ')) topology equations and caller-supplied node identities",
        "SI terminal/branch orientation; incidence dimensionless",
        "none; topology is assembled directly from the published oriented connection matrix",
        "synthetic generic topology; parameter uncertainty is caller-declared",
        "frozen AIMORA bridge family/count/state domain",
        PhysicalModelParameter,
    )
end

function _bridge_topology_signature(family, nodes, valves, passives, groups, provenance, licence, redistribution)
    io = IOBuffer()
    print(io, BRIDGE_TOPOLOGY_SCHEMA, '|', family, '|')
    for node in nodes
        print(io, node.name, ':', node.node, ':', node.role, ';')
    end
    print(io, '|')
    for valve in valves
        print(io, valve.name, ':', valve.from_node, ':', valve.to_node, ':',
            valve.valve_class, ':', valve.group, ':', valve.cell, ';')
    end
    print(io, '|')
    for passive in passives
        print(io, passive.name, ':', passive.from_node, ':', passive.to_node, ':',
            passive.kind, ':', passive.group, ':', passive.cell, ';')
    end
    print(io, '|')
    for group in groups
        print(io, group.name, ':', group.kind, ':', join(group.position_indices, ','), ':')
        for column in axes(group.admitted_states, 2)
            print(io, group.state_names[column], '=')
            for row in axes(group.admitted_states, 1)
                print(io, group.admitted_states[row, column] ? '1' : '0')
            end
            print(io, ';')
        end
    end
    print(io, '|', provenance.source, '|', provenance.units, '|', provenance.transformation,
        '|', provenance.uncertainty, '|', provenance.validity_domain, '|', provenance.nature,
        '|', licence, '|', redistribution)
    return bytes2hex(sha256(take!(io)))
end

function bridge_topology_signature(topology::BridgeTopologyDescriptor)
    signature = _bridge_topology_signature(
        topology.family,
        topology.nodes,
        topology.valve_positions,
        topology.passive_positions,
        topology.state_groups,
        topology.provenance,
        topology.licence,
        topology.redistribution,
    )
    signature == topology.topology_signature || throw(ArgumentError(
        "bridge topology contract was mutated after construction",
    ))
    return signature
end

"""Return the declared-node by ordered-branch incidence matrix with +1 at each from node and -1 at each to node."""
function bridge_topology_incidence(topology::BridgeTopologyDescriptor; include_passives::Bool=true)
    bridge_topology_signature(topology)
    node_order = unique(getfield.(topology.nodes, :node))
    node_lookup = Dict(node => index for (index, node) in enumerate(node_order))
    branches = include_passives ? Any[topology.valve_positions...; topology.passive_positions...] :
        Any[topology.valve_positions...]
    incidence = zeros(Int8, length(node_order), length(branches))
    for (column, branch) in enumerate(branches)
        incidence[node_lookup[branch.from_node], column] = 1
        incidence[node_lookup[branch.to_node], column] = -1
    end
    return (node_order=node_order, incidence=incidence)
end

"""Check one complete requested gate vector against every local finite-state constraint."""
function bridge_topology_state_is_allowed(
    topology::BridgeTopologyDescriptor,
    requested_state::AbstractVector{Bool},
)
    length(requested_state) == length(topology.valve_positions) || throw(DimensionMismatch(
        "bridge requested-state length must match valve-position count",
    ))
    for group in topology.state_groups
        group.kind === :natural_conduction && continue
        local_state = requested_state[group.position_indices]
        any(column -> all(local_state .== view(group.admitted_states, :, column)),
            axes(group.admitted_states, 2)) || return false
    end
    return true
end

_node(name::Symbol, node::Integer, role::Symbol=:external) = BridgeNode(name, node, role)
_valve(name, from, to, class; group=1, cell=0) =
    BridgeValvePosition(name, from, to, class; group, cell)
_passive(name, from, to, kind; group=1, cell=0) =
    BridgePassivePosition(name, from, to, kind; group, cell)

function _rectifier_valves(ac_nodes, dc_positive, dc_negative, mode; group=1)
    mode in (:diode, :thyristor, :half_controlled, :self_commutated) ||
        throw(ArgumentError(
            "bridge mode must be :diode, :thyristor, :half_controlled, or :self_commutated",
        ))
    valves = BridgeValvePosition[]
    for (phase, ac) in enumerate(ac_nodes)
        upper_class = mode === :diode ? :diode :
            mode === :self_commutated ? :self_commutated : :thyristor
        lower_class = mode in (:diode, :half_controlled) ? :diode :
            mode === :self_commutated ? :self_commutated : :thyristor
        if mode === :self_commutated
            push!(valves,
                _valve(Symbol(:phase_, phase, :_upper), dc_positive, ac, upper_class; group),
                _valve(Symbol(:phase_, phase, :_lower), ac, dc_negative, lower_class; group),
            )
        else
            push!(valves,
                _valve(Symbol(:phase_, phase, :_upper), ac, dc_positive, upper_class; group),
                _valve(Symbol(:phase_, phase, :_lower), dc_negative, ac, lower_class; group),
            )
        end
    end
    return valves
end

function _natural_rectifier_state_group(phase_count::Int, first_position::Int=1; name=:natural_paths)
    positions = collect(first_position:(first_position + 2 * phase_count - 1))
    states = falses(length(positions), phase_count * max(1, phase_count - 1))
    names = Symbol[]
    column = 0
    for upper_phase in 1:phase_count, lower_phase in 1:phase_count
        upper_phase == lower_phase && phase_count > 1 && continue
        column += 1
        states[2 * upper_phase - 1, column] = true
        states[2 * lower_phase, column] = true
        push!(names, Symbol(:path_, upper_phase, :_, lower_phase))
    end
    return BridgeStateGroup(name, :natural_conduction, positions, names, states)
end

function _complementary_leg_group(name::Symbol, upper_index::Int, lower_index::Int)
    return BridgeStateGroup(
        name,
        :complementary_leg,
        [upper_index, lower_index],
        [:upper, :lower, :blocked],
        Bool[1 0 0; 0 1 0],
    )
end

function single_phase_graetz_topology(
    ac_nodes::NTuple{2,<:Integer},
    dc_positive::Integer,
    dc_negative::Integer;
    mode::Symbol=:diode,
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:single_phase_graetz),
)
    topology = polyphase_bridge_topology(
        collect(ac_nodes),
        dc_positive,
        dc_negative;
        mode,
        provenance,
    )
    return BridgeTopologyDescriptor(
        :single_phase_graetz,
        topology.nodes,
        topology.valve_positions,
        topology.passive_positions,
        topology.state_groups;
        provenance,
    )
end

function polyphase_bridge_topology(
    ac_nodes::AbstractVector{<:Integer},
    dc_positive::Integer,
    dc_negative::Integer;
    mode::Symbol=:diode,
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:polyphase_bridge),
)
    phase_nodes = Int.(ac_nodes)
    1 <= length(phase_nodes) <= 12 || throw(ArgumentError(
        "polyphase bridge phase count must be between one and twelve",
    ))
    terminals = Int[dc_positive, phase_nodes..., dc_negative]
    length(unique(terminals)) == length(terminals) || throw(ArgumentError(
        "polyphase bridge DC and AC terminals must be distinct",
    ))
    nodes = BridgeNode[_node(:dc_positive, dc_positive)]
    append!(nodes, [_node(Symbol(:ac_, phase), node) for (phase, node) in enumerate(phase_nodes)])
    push!(nodes, _node(:dc_negative, dc_negative))
    valves = _rectifier_valves(phase_nodes, dc_positive, dc_negative, mode)
    groups = if mode === :self_commutated
        [_complementary_leg_group(Symbol(:phase_, phase), 2 * phase - 1, 2 * phase)
            for phase in eachindex(phase_nodes)]
    else
        [_natural_rectifier_state_group(length(phase_nodes))]
    end
    return BridgeTopologyDescriptor(
        :polyphase_bridge,
        nodes,
        valves,
        BridgePassivePosition[],
        groups;
        provenance,
    )
end

function two_level_bridge_topology(
    ac_nodes::AbstractVector{<:Integer},
    dc_positive::Integer,
    dc_negative::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:two_level_bridge),
)
    base = polyphase_bridge_topology(
        ac_nodes,
        dc_positive,
        dc_negative;
        mode=:self_commutated,
        provenance,
    )
    return BridgeTopologyDescriptor(
        :two_level_bridge,
        base.nodes,
        base.valve_positions,
        base.passive_positions,
        base.state_groups;
        provenance,
    )
end

function full_bridge_topology(
    ac_positive::Integer,
    ac_negative::Integer,
    dc_positive::Integer,
    dc_negative::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:full_bridge),
)
    base = two_level_bridge_topology(
        [ac_positive, ac_negative], dc_positive, dc_negative; provenance,
    )
    return BridgeTopologyDescriptor(
        :full_bridge,
        base.nodes,
        base.valve_positions,
        base.passive_positions,
        base.state_groups;
        provenance,
    )
end

function step_down_chopper_topology(
    dc_positive::Integer,
    output::Integer,
    dc_negative::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:step_down_chopper),
)
    nodes = [_node(:dc_positive, dc_positive), _node(:output, output), _node(:dc_negative, dc_negative)]
    valves = [
        _valve(:controlled, dc_positive, output, :self_commutated),
        _valve(:freewheel, dc_negative, output, :diode),
    ]
    group = BridgeStateGroup(:chopper_paths, :admissible_gate_state, [1],
        [:energize, :freewheel], Bool[1 0])
    return BridgeTopologyDescriptor(:step_down_chopper, nodes, valves,
        BridgePassivePosition[], [group]; provenance)
end

function step_up_chopper_topology(
    input_positive::Integer,
    output_positive::Integer,
    dc_negative::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:step_up_chopper),
)
    nodes = [_node(:input_positive, input_positive), _node(:output_positive, output_positive),
        _node(:dc_negative, dc_negative)]
    valves = [
        _valve(:controlled, input_positive, dc_negative, :self_commutated),
        _valve(:boost_diode, input_positive, output_positive, :diode),
    ]
    group = BridgeStateGroup(:chopper_paths, :admissible_gate_state, [1],
        [:charge, :transfer], Bool[1 0])
    return BridgeTopologyDescriptor(:step_up_chopper, nodes, valves,
        BridgePassivePosition[], [group]; provenance)
end

function inverting_buck_boost_topology(
    input_positive::Integer,
    switching::Integer,
    output_negative::Integer,
    reference::Integer;
    provenance::ParameterProvenance=
        generic_bridge_topology_provenance(:inverting_buck_boost),
)
    nodes = [
        _node(:input_positive, input_positive),
        _node(:switching, switching),
        _node(:output_negative, output_negative),
        _node(:reference, reference),
    ]
    valves = [
        _valve(:controlled, input_positive, switching, :self_commutated),
        _valve(:inverting_diode, output_negative, switching, :diode),
    ]
    group = BridgeStateGroup(
        :chopper_paths,
        :admissible_gate_state,
        [1],
        [:charge, :transfer],
        Bool[1 0],
    )
    return BridgeTopologyDescriptor(
        :inverting_buck_boost,
        nodes,
        valves,
        BridgePassivePosition[],
        [group];
        provenance,
    )
end

function bidirectional_chopper_topology(
    dc_positive::Integer,
    output::Integer,
    dc_negative::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:bidirectional_chopper),
)
    base = two_level_bridge_topology([output], dc_positive, dc_negative; provenance)
    return BridgeTopologyDescriptor(:bidirectional_chopper, base.nodes, base.valve_positions,
        base.passive_positions, base.state_groups; provenance)
end

function multigroup_bridge_topology(
    ac_groups::AbstractVector{<:AbstractVector{<:Integer}},
    group_dc_nodes::AbstractVector{<:Tuple};
    mode::Symbol=:diode,
    composition::Symbol=:series,
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:multi_group_bridge),
)
    group_count = length(ac_groups)
    1 <= group_count <= 4 && length(group_dc_nodes) == group_count || throw(ArgumentError(
        "multi-group bridge requires one through four equally described AC/DC groups",
    ))
    composition in (:series, :parallel) || throw(ArgumentError(
        "multi-group DC composition must be :series or :parallel",
    ))
    phases = length(first(ac_groups))
    1 <= phases <= 12 && all(group -> length(group) == phases, ac_groups) ||
        throw(ArgumentError("multi-group AC groups must have the same one-to-twelve phase count"))
    dc_pairs = [(Int(pair[1]), Int(pair[2])) for pair in group_dc_nodes]
    if composition === :series
        all(index -> dc_pairs[index][2] == dc_pairs[index + 1][1], 1:(group_count - 1)) ||
            throw(ArgumentError("series bridge groups require an explicit negative-to-positive DC chain"))
    else
        all(pair -> pair == first(dc_pairs), dc_pairs) || throw(ArgumentError(
            "parallel bridge groups must share exact positive and negative DC nodes",
        ))
    end
    nodes = BridgeNode[
        _node(:dc_positive, first(dc_pairs)[1]),
        _node(:dc_negative, last(dc_pairs)[2]),
    ]
    if composition === :series && group_count > 1
        append!(nodes, [_node(Symbol(:dc_series_, index), dc_pairs[index][2], :internal)
            for index in 1:(group_count - 1)])
    end
    valves = BridgeValvePosition[]
    groups = BridgeStateGroup[]
    for group_index in 1:group_count
        phase_nodes = Int.(ac_groups[group_index])
        append!(nodes, [_node(Symbol(:group_, group_index, :_ac_, phase), node)
            for (phase, node) in enumerate(phase_nodes)])
        first_position = length(valves) + 1
        group_valves = _rectifier_valves(
            phase_nodes,
            dc_pairs[group_index][1],
            dc_pairs[group_index][2],
            mode;
            group=group_index,
        )
        for valve in group_valves
            push!(valves, BridgeValvePosition(
                Symbol(:group_, group_index, :_, valve.name),
                valve.from_node,
                valve.to_node,
                valve.valve_class;
                group=group_index,
            ))
        end
        if mode === :self_commutated
            for phase in 1:phases
                push!(groups, _complementary_leg_group(
                    Symbol(:group_, group_index, :_phase_, phase),
                    first_position + 2 * phase - 2,
                    first_position + 2 * phase - 1,
                ))
            end
        else
            push!(groups, _natural_rectifier_state_group(
                phases,
                first_position;
                name=Symbol(:group_, group_index, :_natural_paths),
            ))
        end
    end
    unique_by_name = Dict(node.name => node for node in nodes)
    return BridgeTopologyDescriptor(
        :multi_group_bridge,
        collect(values(unique_by_name)),
        valves,
        BridgePassivePosition[],
        groups;
        provenance,
    )
end

function neutral_point_clamped_leg_topology(
    dc_positive::Integer,
    midpoint::Integer,
    ac_terminal::Integer,
    dc_negative::Integer,
    upper_clamp_node::Integer,
    lower_clamp_node::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:neutral_point_clamped_leg),
)
    nodes = [
        _node(:dc_positive, dc_positive), _node(:midpoint, midpoint),
        _node(:ac_terminal, ac_terminal), _node(:dc_negative, dc_negative),
        _node(:upper_clamp_node, upper_clamp_node, :internal),
        _node(:lower_clamp_node, lower_clamp_node, :internal),
    ]
    valves = [
        _valve(:outer_upper, dc_positive, upper_clamp_node, :self_commutated),
        _valve(:inner_upper, upper_clamp_node, ac_terminal, :self_commutated),
        _valve(:inner_lower, ac_terminal, lower_clamp_node, :self_commutated),
        _valve(:outer_lower, lower_clamp_node, dc_negative, :self_commutated),
        _valve(:upper_clamp, midpoint, upper_clamp_node, :diode),
        _valve(:lower_clamp, lower_clamp_node, midpoint, :diode),
    ]
    states = Bool[
        1 0 0 0;
        1 1 0 0;
        0 1 1 0;
        0 0 1 0;
    ]
    group = BridgeStateGroup(:three_level_gate_states, :admissible_gate_state,
        [1, 2, 3, 4], [:positive, :neutral, :negative, :blocked], states)
    return BridgeTopologyDescriptor(:neutral_point_clamped_leg, nodes, valves,
        BridgePassivePosition[], [group]; provenance)
end

function t_type_leg_topology(
    dc_positive::Integer,
    midpoint::Integer,
    ac_terminal::Integer,
    dc_negative::Integer,
    midpoint_path_node::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:t_type_leg),
)
    nodes = [
        _node(:dc_positive, dc_positive), _node(:midpoint, midpoint),
        _node(:ac_terminal, ac_terminal), _node(:dc_negative, dc_negative),
        _node(:midpoint_path_node, midpoint_path_node, :internal),
    ]
    valves = [
        _valve(:outer_upper, dc_positive, ac_terminal, :self_commutated),
        _valve(:midpoint_ac_side, ac_terminal, midpoint_path_node, :self_commutated),
        _valve(:midpoint_dc_side, midpoint, midpoint_path_node, :self_commutated),
        _valve(:outer_lower, ac_terminal, dc_negative, :self_commutated),
    ]
    states = Bool[
        1 0 0 0;
        0 1 0 0;
        0 1 0 0;
        0 0 1 0;
    ]
    group = BridgeStateGroup(:three_level_gate_states, :admissible_gate_state,
        collect(1:4), [:positive, :neutral, :negative, :blocked], states)
    return BridgeTopologyDescriptor(:t_type_leg, nodes, valves,
        BridgePassivePosition[], [group]; provenance)
end

function flying_capacitor_leg_topology(
    dc_positive::Integer,
    ac_terminal::Integer,
    dc_negative::Integer,
    upper_flying_node::Integer,
    lower_flying_node::Integer;
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:flying_capacitor_leg),
)
    nodes = [
        _node(:dc_positive, dc_positive), _node(:ac_terminal, ac_terminal),
        _node(:dc_negative, dc_negative),
        _node(:upper_flying_node, upper_flying_node, :internal),
        _node(:lower_flying_node, lower_flying_node, :internal),
    ]
    valves = [
        _valve(:outer_upper, dc_positive, upper_flying_node, :self_commutated),
        _valve(:inner_upper, upper_flying_node, ac_terminal, :self_commutated),
        _valve(:inner_lower, ac_terminal, lower_flying_node, :self_commutated),
        _valve(:outer_lower, lower_flying_node, dc_negative, :self_commutated),
    ]
    passives = [_passive(:flying_capacitor, upper_flying_node, lower_flying_node, :capacitor)]
    states = Bool[
        1 1 0 0 0;
        1 0 1 0 0;
        0 1 0 1 0;
        0 0 1 1 0;
    ]
    group = BridgeStateGroup(:three_level_gate_states, :admissible_gate_state,
        collect(1:4), [:positive, :zero_charge, :zero_discharge, :negative, :blocked], states)
    return BridgeTopologyDescriptor(:flying_capacitor_leg, nodes, valves,
        passives, [group]; provenance)
end

function cascaded_h_bridge_phase_topology(
    phase_positive::Integer,
    phase_negative::Integer,
    cell_dc_nodes::AbstractVector{<:Tuple},
    series_junction_nodes::AbstractVector{<:Integer};
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:cascaded_h_bridge_phase),
)
    cell_count = length(cell_dc_nodes)
    1 <= cell_count <= 16 || throw(ArgumentError(
        "cascaded H-bridge phase requires one through sixteen cells",
    ))
    length(series_junction_nodes) == max(0, cell_count - 1) || throw(ArgumentError(
        "cascaded H-bridge requires exactly one fewer series junction than cells",
    ))
    junctions = Int[phase_positive, Int.(series_junction_nodes)..., phase_negative]
    nodes = BridgeNode[_node(:phase_positive, phase_positive), _node(:phase_negative, phase_negative)]
    append!(nodes, [_node(Symbol(:series_junction_, index), node, :internal)
        for (index, node) in enumerate(series_junction_nodes)])
    valves = BridgeValvePosition[]
    groups = BridgeStateGroup[]
    for cell in 1:cell_count
        dc_positive, dc_negative = Int.(cell_dc_nodes[cell])
        push!(nodes,
            _node(Symbol(:cell_, cell, :_dc_positive), dc_positive),
            _node(Symbol(:cell_, cell, :_dc_negative), dc_negative),
        )
        left = junctions[cell]
        right = junctions[cell + 1]
        first_position = length(valves) + 1
        append!(valves, [
            _valve(Symbol(:cell_, cell, :_left_upper), dc_positive, left,
                :self_commutated; group=cell, cell),
            _valve(Symbol(:cell_, cell, :_left_lower), left, dc_negative,
                :self_commutated; group=cell, cell),
            _valve(Symbol(:cell_, cell, :_right_upper), dc_positive, right,
                :self_commutated; group=cell, cell),
            _valve(Symbol(:cell_, cell, :_right_lower), right, dc_negative,
                :self_commutated; group=cell, cell),
        ])
        states = Bool[
            1 0 1 0 0;
            0 1 0 1 0;
            1 0 0 1 0;
            0 1 1 0 0;
        ]
        push!(groups, BridgeStateGroup(
            Symbol(:cell_, cell, :_states),
            :admissible_gate_state,
            collect(first_position:(first_position + 3)),
            [:zero_upper, :zero_lower, :positive, :negative, :blocked],
            states,
        ))
    end
    return BridgeTopologyDescriptor(:cascaded_h_bridge_phase, nodes, valves,
        BridgePassivePosition[], groups; provenance)
end

function matrix_converter_topology(
    input_nodes::NTuple{3,<:Integer},
    output_nodes::NTuple{3,<:Integer};
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:matrix_converter),
)
    inputs = Int.(input_nodes)
    outputs = Int.(output_nodes)
    all(>=(0), (inputs..., outputs...)) &&
        length(unique((inputs..., outputs...))) == 6 || throw(ArgumentError(
        "matrix-converter input and output terminals must be six distinct nonnegative nodes",
    ))
    nodes = BridgeNode[
        [_node(Symbol(:input_, phase), inputs[phase]) for phase in 1:3]...,
        [_node(Symbol(:output_, phase), outputs[phase]) for phase in 1:3]...,
    ]
    valves = BridgeValvePosition[]
    for output_phase in 1:3, input_phase in 1:3
        connection = 3 * (output_phase - 1) + input_phase
        push!(
            valves,
            _valve(
                Symbol(:output_, output_phase, :_input_, input_phase, :_forward),
                inputs[input_phase],
                outputs[output_phase],
                :self_commutated;
                group=output_phase,
                cell=connection,
            ),
            _valve(
                Symbol(:output_, output_phase, :_input_, input_phase, :_reverse),
                outputs[output_phase],
                inputs[input_phase],
                :self_commutated;
                group=output_phase,
                cell=connection,
            ),
        )
    end
    states = falses(length(valves), 27)
    state_names = Symbol[]
    column = 0
    for input_for_output_1 in 1:3, input_for_output_2 in 1:3, input_for_output_3 in 1:3
        column += 1
        selected_inputs = (input_for_output_1, input_for_output_2, input_for_output_3)
        for output_phase in 1:3
            connection = 3 * (output_phase - 1) + selected_inputs[output_phase]
            states[2 * connection - 1, column] = true
            states[2 * connection, column] = true
        end
        push!(state_names, Symbol(:connection_, join(selected_inputs, :_)))
    end
    group = BridgeStateGroup(
        :matrix_connection_states,
        :admissible_gate_state,
        collect(eachindex(valves)),
        state_names,
        states,
    )
    return BridgeTopologyDescriptor(
        :matrix_converter,
        nodes,
        valves,
        BridgePassivePosition[],
        [group];
        provenance,
    )
end

function _cycloconverter_bridge_paths(first_position::Int)
    paths = NTuple{2,Int}[]
    for upper_phase in 1:3, lower_phase in 1:3
        upper_phase == lower_phase && continue
        push!(paths, (
            first_position + 2 * upper_phase - 2,
            first_position + 2 * lower_phase - 1,
        ))
    end
    return paths
end

function _cycloconverter_state_group(
    output_phase::Int,
    positive_first::Int,
    negative_first::Int,
    circulating_current::Bool,
)
    positive_paths = _cycloconverter_bridge_paths(positive_first)
    negative_paths = _cycloconverter_bridge_paths(negative_first)
    positions = collect(positive_first:(negative_first + 5))
    state_count = circulating_current ? 49 : 13
    states = falses(length(positions), state_count)
    names = Symbol[:blocked]
    column = 1
    function add_path!(prefix, positive_path, negative_path)
        column += 1
        positive_path === nothing || foreach(index -> states[index - positive_first + 1, column] = true, positive_path)
        negative_path === nothing || foreach(index -> states[index - positive_first + 1, column] = true, negative_path)
        push!(names, Symbol(prefix, :_, column - 1))
    end
    for path in positive_paths
        add_path!(:positive, path, nothing)
    end
    for path in negative_paths
        add_path!(:negative, nothing, path)
    end
    if circulating_current
        for positive_path in positive_paths, negative_path in negative_paths
            add_path!(:circulating, positive_path, negative_path)
        end
    end
    return BridgeStateGroup(
        Symbol(:output_, output_phase, :_bridge_group_states),
        :admissible_gate_state,
        positions,
        names,
        states,
    )
end

function cycloconverter_topology(
    input_nodes::NTuple{3,<:Integer},
    output_nodes,
    neutral_node::Integer;
    circulating_current::Bool=false,
    provenance::ParameterProvenance=generic_bridge_topology_provenance(:cycloconverter),
)
    inputs = Int.(input_nodes)
    outputs = Tuple(Int.(output_nodes))
    length(outputs) in (1, 3) || throw(ArgumentError(
        "cycloconverter output must contain one or three phases",
    ))
    neutral = Int(neutral_node)
    terminals = (inputs..., outputs..., neutral)
    all(>=(0), terminals) && length(unique(terminals)) == length(terminals) ||
        throw(ArgumentError("cycloconverter terminals must be distinct nonnegative nodes"))
    nodes = BridgeNode[
        [_node(Symbol(:input_, phase), inputs[phase]) for phase in 1:3]...,
        [_node(Symbol(:output_, phase), outputs[phase]) for phase in eachindex(outputs)]...,
        _node(:output_neutral, neutral),
    ]
    valves = BridgeValvePosition[]
    groups = BridgeStateGroup[]
    for output_phase in eachindex(outputs)
        positive_first = length(valves) + 1
        positive = _rectifier_valves(
            inputs,
            outputs[output_phase],
            neutral,
            :thyristor;
            group=2 * output_phase - 1,
        )
        for valve in positive
            push!(valves, BridgeValvePosition(
                Symbol(:output_, output_phase, :_positive_, valve.name),
                valve.from_node,
                valve.to_node,
                valve.valve_class;
                group=2 * output_phase - 1,
                cell=output_phase,
            ))
        end
        negative_first = length(valves) + 1
        negative = _rectifier_valves(
            inputs,
            neutral,
            outputs[output_phase],
            :thyristor;
            group=2 * output_phase,
        )
        for valve in negative
            push!(valves, BridgeValvePosition(
                Symbol(:output_, output_phase, :_negative_, valve.name),
                valve.from_node,
                valve.to_node,
                valve.valve_class;
                group=2 * output_phase,
                cell=output_phase,
            ))
        end
        push!(groups, _cycloconverter_state_group(
            output_phase,
            positive_first,
            negative_first,
            circulating_current,
        ))
    end
    return BridgeTopologyDescriptor(
        :cycloconverter,
        nodes,
        valves,
        BridgePassivePosition[],
        groups;
        provenance,
    )
end

end

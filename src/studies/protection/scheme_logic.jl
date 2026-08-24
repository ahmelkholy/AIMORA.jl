abstract type AbstractProtectionLogicNode end

struct ProtectionLogicInputNode <: AbstractProtectionLogicNode
    id::Symbol
    input::Symbol

    function ProtectionLogicInputNode(id::Symbol, input::Symbol)
        all(name -> !isempty(String(name)), (id, input)) || throw(ArgumentError(
            "protection logic input identities must not be empty",
        ))
        return new(id, input)
    end
end

struct ProtectionLogicNotNode <: AbstractProtectionLogicNode
    id::Symbol
    source::Symbol

    function ProtectionLogicNotNode(id::Symbol, source::Symbol)
        all(name -> !isempty(String(name)), (id, source)) || throw(ArgumentError(
            "protection logic NOT identities must not be empty",
        ))
        id != source || throw(ArgumentError("protection logic NOT node cannot depend on itself"))
        return new(id, source)
    end
end

struct ProtectionLogicVoteNode <: AbstractProtectionLogicNode
    id::Symbol
    sources::Vector{Symbol}
    minimum_true::Int

    function ProtectionLogicVoteNode(
        id::Symbol,
        sources::AbstractVector{Symbol},
        minimum_true::Integer,
    )
        inputs = Symbol.(sources)
        !isempty(String(id)) && !isempty(inputs) &&
            all(name -> !isempty(String(name)), inputs) || throw(ArgumentError(
            "protection vote node identities must not be empty",
        ))
        allunique(inputs) || throw(ArgumentError(
            "protection vote node repeats a source",
        ))
        threshold = Int(minimum_true)
        1 <= threshold <= length(inputs) || throw(ArgumentError(
            "protection vote threshold must be within its source count",
        ))
        id in inputs && throw(ArgumentError("protection vote node cannot depend on itself"))
        return new(id, inputs, threshold)
    end
end

struct ProtectionLogicLatchNode <: AbstractProtectionLogicNode
    id::Symbol
    set_source::Symbol
    reset_source::Symbol
    reset_dominant::Bool

    function ProtectionLogicLatchNode(
        id::Symbol,
        set_source::Symbol,
        reset_source::Symbol;
        reset_dominant::Bool=true,
    )
        all(name -> !isempty(String(name)), (id, set_source, reset_source)) ||
            throw(ArgumentError("protection latch identities must not be empty"))
        set_source != reset_source || throw(ArgumentError(
            "protection latch set and reset sources must differ",
        ))
        id in (set_source, reset_source) && throw(ArgumentError(
            "protection latch cannot depend on itself",
        ))
        return new(id, set_source, reset_source, reset_dominant)
    end
end

_protection_logic_node_id(node::AbstractProtectionLogicNode) = node.id
_protection_logic_dependencies(::ProtectionLogicInputNode) = Symbol[]
_protection_logic_dependencies(node::ProtectionLogicNotNode) = Symbol[node.source]
_protection_logic_dependencies(node::ProtectionLogicVoteNode) = node.sources
_protection_logic_dependencies(node::ProtectionLogicLatchNode) =
    Symbol[node.set_source, node.reset_source]

struct ProtectionLogicDefinition
    id::Symbol
    nodes::Vector{AbstractProtectionLogicNode}
    output_node::Symbol
    deterministic_signature_sha256::String

    function ProtectionLogicDefinition(
        id::Symbol,
        nodes::AbstractVector{<:AbstractProtectionLogicNode},
        output_node::Symbol,
    )
        !isempty(String(id)) && !isempty(String(output_node)) || throw(ArgumentError(
            "protection logic definition identities must not be empty",
        ))
        ordered_nodes = AbstractProtectionLogicNode[nodes...]
        isempty(ordered_nodes) && throw(ArgumentError(
            "protection logic definition must contain at least one node",
        ))
        node_ids = _protection_logic_node_id.(ordered_nodes)
        allunique(node_ids) || throw(ArgumentError(
            "protection logic definition repeats a node identity",
        ))
        available = Set{Symbol}()
        for node in ordered_nodes
            dependencies = _protection_logic_dependencies(node)
            all(in(available), dependencies) || throw(ArgumentError(
                "protection logic nodes must be supplied in acyclic dependency order",
            ))
            push!(available, _protection_logic_node_id(node))
        end
        output_node in available || throw(ArgumentError(
            "protection logic output node is absent",
        ))
        io = IOBuffer()
        print(io, id, '|', output_node)
        for node in ordered_nodes
            print(io, '|', repr(node))
        end
        return new(id, ordered_nodes, output_node, bytes2hex(sha256(take!(io))))
    end
end

mutable struct ProtectionLogicRuntime
    latches::Dict{Symbol,Bool}
    evaluation_count::Int
    last_values::Dict{Symbol,Bool}

    function ProtectionLogicRuntime(definition::ProtectionLogicDefinition)
        latches = Dict{Symbol,Bool}(
            node.id => false
            for node in definition.nodes if node isa ProtectionLogicLatchNode
        )
        return new(latches, 0, Dict{Symbol,Bool}())
    end
end

struct ProtectionLogicSnapshot
    definition_signature_sha256::String
    latches::Dict{Symbol,Bool}
    evaluation_count::Int
    last_values::Dict{Symbol,Bool}
end

struct ProtectionLogicResult
    definition_id::Symbol
    values::Dict{Symbol,Bool}
    output::Bool
    evaluation_index::Int
    deterministic_signature_sha256::String
end

function protection_logic_snapshot(
    runtime::ProtectionLogicRuntime,
    definition::ProtectionLogicDefinition,
)
    return ProtectionLogicSnapshot(
        definition.deterministic_signature_sha256,
        copy(runtime.latches),
        runtime.evaluation_count,
        copy(runtime.last_values),
    )
end

function restore_protection_logic_snapshot!(
    runtime::ProtectionLogicRuntime,
    definition::ProtectionLogicDefinition,
    snapshot::ProtectionLogicSnapshot,
)
    snapshot.definition_signature_sha256 == definition.deterministic_signature_sha256 ||
        throw(ArgumentError("protection logic snapshot definition identity is stale"))
    latch_ids = Set(
        node.id for node in definition.nodes if node isa ProtectionLogicLatchNode
    )
    Set(keys(snapshot.latches)) == latch_ids || throw(ArgumentError(
        "protection logic snapshot latch identities are invalid",
    ))
    node_ids = Set(_protection_logic_node_id.(definition.nodes))
    isempty(snapshot.last_values) || Set(keys(snapshot.last_values)) == node_ids ||
        throw(ArgumentError("protection logic snapshot value identities are invalid"))
    snapshot.evaluation_count >= 0 &&
        (snapshot.evaluation_count == 0) == isempty(snapshot.last_values) ||
        throw(ArgumentError("protection logic snapshot evaluation state is invalid"))
    runtime.latches = copy(snapshot.latches)
    runtime.evaluation_count = snapshot.evaluation_count
    runtime.last_values = copy(snapshot.last_values)
    return runtime
end

function evaluate_protection_logic!(
    runtime::ProtectionLogicRuntime,
    definition::ProtectionLogicDefinition,
    inputs::AbstractDict{Symbol,Bool},
)
    values = Dict{Symbol,Bool}()
    for node in definition.nodes
        if node isa ProtectionLogicInputNode
            haskey(inputs, node.input) || throw(ArgumentError(
                "protection logic input $(node.input) is unavailable",
            ))
            values[node.id] = inputs[node.input]
        elseif node isa ProtectionLogicNotNode
            values[node.id] = !values[node.source]
        elseif node isa ProtectionLogicVoteNode
            values[node.id] = count(source -> values[source], node.sources) >=
                node.minimum_true
        else
            latch = node::ProtectionLogicLatchNode
            set_value = values[latch.set_source]
            reset_value = values[latch.reset_source]
            prior = get(runtime.latches, latch.id, false)
            updated = if set_value && reset_value
                !latch.reset_dominant
            elseif reset_value
                false
            elseif set_value
                true
            else
                prior
            end
            runtime.latches[latch.id] = updated
            values[latch.id] = updated
        end
    end
    runtime.evaluation_count += 1
    runtime.last_values = copy(values)
    io = IOBuffer()
    print(io, definition.deterministic_signature_sha256, '|', runtime.evaluation_count)
    for node in definition.nodes
        id = _protection_logic_node_id(node)
        print(io, '|', id, '=', values[id])
    end
    return ProtectionLogicResult(
        definition.id,
        values,
        values[definition.output_node],
        runtime.evaluation_count,
        bytes2hex(sha256(take!(io))),
    )
end

@enum ProtectionMessageDisposition::UInt8 begin
    PassProtectionMessage = 0x01
    DropProtectionMessage = 0x02
end

struct ProtectionMessageRule
    sequence_number::Int
    disposition::ProtectionMessageDisposition
    copy_count::Int
    additional_delay_ticks::Int

    function ProtectionMessageRule(
        sequence_number::Integer;
        disposition::ProtectionMessageDisposition=PassProtectionMessage,
        copy_count::Integer=1,
        additional_delay_ticks::Integer=0,
    )
        sequence = Int(sequence_number)
        copies = Int(copy_count)
        additional_delay = Int(additional_delay_ticks)
        sequence > 0 || throw(ArgumentError(
            "protection message rule sequence must be positive",
        ))
        1 <= copies <= 16 || throw(ArgumentError(
            "protection message copy count must be in 1:16",
        ))
        additional_delay >= 0 || throw(ArgumentError(
            "protection message additional delay must be nonnegative",
        ))
        disposition === DropProtectionMessage && copies != 1 && throw(ArgumentError(
            "a dropped protection message cannot declare duplicates",
        ))
        return new(sequence, disposition, copies, additional_delay)
    end
end

struct ProtectionCommunicationLink
    id::Symbol
    sender::Symbol
    receiver::Symbol
    allowed_payloads::Vector{Symbol}
    fixed_delay_ticks::Int
    rules::Dict{Int,ProtectionMessageRule}
    provenance::ParameterProvenance
    deterministic_signature_sha256::String

    function ProtectionCommunicationLink(
        id::Symbol,
        sender::Symbol,
        receiver::Symbol;
        allowed_payloads::AbstractVector{Symbol},
        fixed_delay_ticks::Integer,
        rules::AbstractVector{ProtectionMessageRule}=ProtectionMessageRule[],
        provenance::ParameterProvenance,
    )
        all(name -> !isempty(String(name)), (id, sender, receiver)) ||
            throw(ArgumentError("protection communication identities must not be empty"))
        sender != receiver || throw(ArgumentError(
            "protection communication endpoints must differ",
        ))
        payloads = Symbol.(allowed_payloads)
        !isempty(payloads) && all(name -> !isempty(String(name)), payloads) &&
            allunique(payloads) || throw(ArgumentError(
            "protection communication payloads must be nonempty and unique",
        ))
        delay = Int(fixed_delay_ticks)
        delay >= 0 || throw(ArgumentError(
            "protection communication delay must be nonnegative",
        ))
        rule_map = Dict(rule.sequence_number => rule for rule in rules)
        length(rule_map) == length(rules) || throw(ArgumentError(
            "protection communication repeats a sequence rule",
        ))
        provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
            "protection communication policy provenance must be numerical",
        ))
        io = IOBuffer()
        print(io, id, '|', sender, '|', receiver, '|', join(string.(payloads), ','), '|', delay)
        for sequence in sort(collect(keys(rule_map)))
            print(io, '|', repr(rule_map[sequence]))
        end
        return new(
            id,
            sender,
            receiver,
            payloads,
            delay,
            rule_map,
            provenance,
            bytes2hex(sha256(take!(io))),
        )
    end
end

struct QueuedProtectionMessage
    link_id::Symbol
    sender::Symbol
    receiver::Symbol
    payload::Symbol
    sequence_number::Int
    copy_index::Int
    send_tick::Int
    delivery_tick::Int
end

mutable struct ProtectionCommunicationRuntime
    next_sequence_number::Int
    queue::Vector{QueuedProtectionMessage}
    last_send_tick::Int
    last_delivery_tick::Int
    sent_count::Int
    delivered_count::Int
    dropped_count::Int

    ProtectionCommunicationRuntime() = new(1, QueuedProtectionMessage[], -1, -1, 0, 0, 0)
end

struct ProtectionMessageSendResult
    sequence_number::Int
    disposition::ProtectionMessageDisposition
    queued_copy_count::Int
    delivery_tick::Union{Nothing,Int}
end

struct ProtectionCommunicationSnapshot
    link_signature_sha256::String
    next_sequence_number::Int
    queue::Vector{QueuedProtectionMessage}
    last_send_tick::Int
    last_delivery_tick::Int
    sent_count::Int
    delivered_count::Int
    dropped_count::Int
end

function send_protection_message!(
    runtime::ProtectionCommunicationRuntime,
    link::ProtectionCommunicationLink,
    payload::Symbol,
    send_tick::Integer,
)
    payload in link.allowed_payloads || throw(ArgumentError(
        "protection message payload is not admitted by its link",
    ))
    tick = Int(send_tick)
    tick >= 0 && tick >= runtime.last_send_tick || throw(ArgumentError(
        "protection message send ticks must be nonnegative and monotone",
    ))
    sequence = runtime.next_sequence_number
    runtime.next_sequence_number += 1
    runtime.last_send_tick = tick
    runtime.sent_count += 1
    rule = get(
        link.rules,
        sequence,
        ProtectionMessageRule(sequence),
    )
    if rule.disposition === DropProtectionMessage
        runtime.dropped_count += 1
        return ProtectionMessageSendResult(sequence, rule.disposition, 0, nothing)
    end
    delivery_tick = tick + link.fixed_delay_ticks + rule.additional_delay_ticks
    for copy_index in 1:rule.copy_count
        push!(runtime.queue, QueuedProtectionMessage(
            link.id,
            link.sender,
            link.receiver,
            payload,
            sequence,
            copy_index,
            tick,
            delivery_tick,
        ))
    end
    sort!(runtime.queue; by=message -> (
        message.delivery_tick,
        message.sequence_number,
        message.copy_index,
    ))
    return ProtectionMessageSendResult(
        sequence,
        rule.disposition,
        rule.copy_count,
        delivery_tick,
    )
end

function deliver_due_protection_messages!(
    runtime::ProtectionCommunicationRuntime,
    link::ProtectionCommunicationLink,
    delivery_tick::Integer,
)
    tick = Int(delivery_tick)
    tick >= 0 && tick >= runtime.last_delivery_tick || throw(ArgumentError(
        "protection message delivery ticks must be nonnegative and monotone",
    ))
    !isempty(runtime.queue) && first(runtime.queue).delivery_tick < tick &&
        throw(ArgumentError("a queued protection-message delivery tick was skipped"))
    due_count = count(message -> message.delivery_tick == tick, runtime.queue)
    due = due_count == 0 ? QueuedProtectionMessage[] : splice!(runtime.queue, 1:due_count)
    all(message -> message.link_id == link.id, due) || throw(ArgumentError(
        "protection communication runtime contains a foreign link message",
    ))
    runtime.last_delivery_tick = tick
    runtime.delivered_count += length(due)
    return due
end

function protection_communication_snapshot(
    runtime::ProtectionCommunicationRuntime,
    link::ProtectionCommunicationLink,
)
    return ProtectionCommunicationSnapshot(
        link.deterministic_signature_sha256,
        runtime.next_sequence_number,
        deepcopy(runtime.queue),
        runtime.last_send_tick,
        runtime.last_delivery_tick,
        runtime.sent_count,
        runtime.delivered_count,
        runtime.dropped_count,
    )
end

function restore_protection_communication_snapshot!(
    runtime::ProtectionCommunicationRuntime,
    link::ProtectionCommunicationLink,
    snapshot::ProtectionCommunicationSnapshot,
)
    snapshot.link_signature_sha256 == link.deterministic_signature_sha256 ||
        throw(ArgumentError("protection communication snapshot link identity is stale"))
    snapshot.next_sequence_number > 0 &&
        snapshot.last_send_tick >= -1 && snapshot.last_delivery_tick >= -1 &&
        all(>=(0), (
            snapshot.sent_count,
            snapshot.delivered_count,
            snapshot.dropped_count,
        )) || throw(ArgumentError("protection communication snapshot counters are invalid"))
    issorted(snapshot.queue; by=message -> (
        message.delivery_tick,
        message.sequence_number,
        message.copy_index,
    )) || throw(ArgumentError("protection communication snapshot queue is not ordered"))
    all(message -> message.link_id == link.id &&
        message.delivery_tick >= message.send_tick,
        snapshot.queue,
    ) || throw(ArgumentError("protection communication snapshot queue is invalid"))
    runtime.next_sequence_number = snapshot.next_sequence_number
    runtime.queue = deepcopy(snapshot.queue)
    runtime.last_send_tick = snapshot.last_send_tick
    runtime.last_delivery_tick = snapshot.last_delivery_tick
    runtime.sent_count = snapshot.sent_count
    runtime.delivered_count = snapshot.delivered_count
    runtime.dropped_count = snapshot.dropped_count
    return runtime
end

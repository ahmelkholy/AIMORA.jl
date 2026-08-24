@enum ProtectionTaskStage::UInt8 begin
    ProtectionMeasurementReleaseStage = 0x01
    ProtectionElementEvaluationStage = 0x02
    ProtectionLocalLogicStage = 0x03
    ProtectionMessageSendStage = 0x04
    ProtectionMessageDeliveryStage = 0x05
    ProtectionRemoteLogicStage = 0x06
    ProtectionTripCommandStage = 0x07
end

const _PROTECTION_TASK_STAGES = (
    ProtectionMeasurementReleaseStage,
    ProtectionElementEvaluationStage,
    ProtectionLocalLogicStage,
    ProtectionMessageSendStage,
    ProtectionMessageDeliveryStage,
    ProtectionRemoteLogicStage,
    ProtectionTripCommandStage,
)

struct ProtectionTaskOperations{M,E,L,S,D,R,T}
    measurement_release::M
    element_evaluation::E
    local_logic::L
    message_send::S
    message_delivery::D
    remote_logic::R
    trip_command::T
end

function protection_task_operation(
    operations::ProtectionTaskOperations,
    stage::ProtectionTaskStage,
)
    stage === ProtectionMeasurementReleaseStage && return operations.measurement_release
    stage === ProtectionElementEvaluationStage && return operations.element_evaluation
    stage === ProtectionLocalLogicStage && return operations.local_logic
    stage === ProtectionMessageSendStage && return operations.message_send
    stage === ProtectionMessageDeliveryStage && return operations.message_delivery
    stage === ProtectionRemoteLogicStage && return operations.remote_logic
    return operations.trip_command
end

struct ProtectionTaskPipeline{O,S}
    plan::EMTTaskPlan
    operations::O
    initial_state::S
    deterministic_signature_sha256::String
end

struct ProtectionTaskAdvanceResult{S}
    instant::EMTLogicalTime
    task_count::Int
    stages::Vector{ProtectionTaskStage}
    state::S
    deterministic_signature_sha256::String

    function ProtectionTaskAdvanceResult(
        instant::EMTLogicalTime,
        task_count::Integer,
        stages::AbstractVector{ProtectionTaskStage},
        state::S,
        deterministic_signature_sha256::AbstractString,
    ) where {S}
        count = Int(task_count)
        count >= 0 || throw(ArgumentError(
            "protection task advance count must be nonnegative",
        ))
        ordered_stages = collect(stages)
        count == length(ordered_stages) || throw(ArgumentError(
            "protection task advance count and stage trace differ",
        ))
        signature = lowercase(String(deterministic_signature_sha256))
        occursin(r"^[0-9a-f]{64}$", signature) || throw(ArgumentError(
            "protection task advance signature must be a 64-hex SHA-256",
        ))
        return new{S}(instant, count, ordered_stages, state, signature)
    end
end

function _protection_task_specification(
    name::String,
    family::EMTTaskFamily,
    epoch::EMTLogicalTime,
    period::EMTLogicalTime,
    priority::Int,
    read_resource::String,
    write_resource::String;
    predecessor::Union{Nothing,String}=nothing,
    effects::Vector{EMTTaskEffect}=EMTTaskEffect[],
)
    return EMTTaskSpec(
        name,
        family,
        epoch,
        period,
        emt_logical_time(0),
        emt_logical_time(0);
        priority,
        read_resources=[read_resource],
        write_resources=[write_resource],
        predecessors=predecessor === nothing ? String[] : [predecessor],
        effects,
    )
end

function ProtectionTaskPipeline(
    operations::ProtectionTaskOperations,
    initial_state;
    epoch::EMTLogicalTime,
    period::EMTLogicalTime,
    start::EMTLogicalTime,
    stop::EMTLogicalTime,
)
    names = (
        "protection_measurement_release",
        "protection_element_evaluation",
        "protection_local_logic",
        "protection_message_send",
        "protection_message_delivery",
        "protection_remote_logic",
        "protection_trip_command",
    )
    families = (
        ProtectionEMTTask,
        ProtectionEMTTask,
        ProtectionEMTTask,
        CarrierEMTTask,
        CarrierEMTTask,
        ProtectionEMTTask,
        ProtectionEMTTask,
    )
    resources = (
        "protection.accepted_measurement",
        "protection.element_decision",
        "protection.local_logic",
        "protection.message_queue",
        "protection.delivered_message",
        "protection.remote_logic",
        "protection.trip_command",
        "emt.topology",
    )
    specifications = EMTTaskSpec[]
    for index in eachindex(names)
        push!(
            specifications,
            _protection_task_specification(
                names[index],
                families[index],
                epoch,
                period,
                index * 10,
                resources[index],
                resources[index + 1];
                predecessor=index == 1 ? nothing : names[index - 1],
                effects=index == length(names) ? [InvalidateEMTTopology] : EMTTaskEffect[],
            ),
        )
    end
    plan = emt_task_plan(specifications; start, stop)
    stage_names = join(UInt8.(_PROTECTION_TASK_STAGES), ',')
    signature = bytes2hex(sha256(join((
        "aimora-protection-task-pipeline-v1",
        plan.signature_sha256,
        stage_names,
        string(typeof(initial_state)),
    ), '\n')))
    return ProtectionTaskPipeline(plan, operations, initial_state, signature)
end

protection_task_plan(pipeline::ProtectionTaskPipeline) = pipeline.plan

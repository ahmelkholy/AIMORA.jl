export RealtimeTargetKind,
       NativeLinuxSoftwareTarget,
       LocalUDPControllerTarget,
       SharedLibraryControllerTarget,
       PhysicalHILTarget,
       RealtimeOverrunPolicy,
       FailSafeOnOverrun,
       MeasureOnlyOverrun,
       RealtimeTarget,
       RealtimeChannel,
       RealtimeUnavailableResult,
       RealtimePreparation,
       RealtimeExecutionMetadata,
       RealtimeSafetyState,
       prepare_realtime_target,
       physical_channel_value,
       raw_channel_value,
       latch_realtime_safe_shutdown!

@enum RealtimeTargetKind begin
    NativeLinuxSoftwareTarget
    LocalUDPControllerTarget
    SharedLibraryControllerTarget
    PhysicalHILTarget
end

@enum RealtimeOverrunPolicy begin
    FailSafeOnOverrun
    MeasureOnlyOverrun
end

struct RealtimeTarget
    id::Symbol
    kind::RealtimeTargetKind
    period_ns::Int64
    step_count::Int
    overrun_policy::RealtimeOverrunPolicy
    require_hard_realtime::Bool
    scheduler_policy::Symbol
    scheduler_priority::Int
    cpu_index::Union{Nothing,Int}
    lock_memory::Bool
    maximum_channels::Int
    maximum_payload_bytes::Int

    function RealtimeTarget(
        id::Symbol,
        kind::RealtimeTargetKind;
        period_ns::Integer,
        step_count::Integer,
        overrun_policy::RealtimeOverrunPolicy=FailSafeOnOverrun,
        require_hard_realtime::Bool=false,
        scheduler_policy::Symbol=:other,
        scheduler_priority::Integer=0,
        cpu_index::Union{Nothing,Integer}=nothing,
        lock_memory::Bool=false,
        maximum_channels::Integer=4096,
        maximum_payload_bytes::Integer=64 * 1024,
    )
        id == Symbol("") && throw(ArgumentError("real-time target id must not be empty"))
        period = Int64(period_ns)
        steps = Int(step_count)
        priority = Int(scheduler_priority)
        cpu = cpu_index === nothing ? nothing : Int(cpu_index)
        channels = Int(maximum_channels)
        payload = Int(maximum_payload_bytes)
        10_000 <= period <= 100_000_000 || throw(ArgumentError(
            "real-time period_ns must be between 10 microseconds and 100 milliseconds",
        ))
        0 < steps <= 1_000_000 || throw(ArgumentError(
            "real-time step_count must be between one and one million",
        ))
        Base.Checked.checked_mul(period, Int64(steps))
        scheduler_policy in (:other, :fifo, :round_robin) || throw(ArgumentError(
            "real-time scheduler_policy must be :other, :fifo, or :round_robin",
        ))
        if scheduler_policy == :other
            priority == 0 || throw(ArgumentError(
                "the ordinary scheduler requires priority zero",
            ))
        else
            1 <= priority <= 99 || throw(ArgumentError(
                "real-time scheduler priority must be between one and 99",
            ))
        end
        cpu === nothing || cpu >= 0 || throw(ArgumentError(
            "real-time cpu_index must be nonnegative",
        ))
        0 < channels <= 4096 || throw(ArgumentError(
            "real-time maximum_channels must be between one and 4096",
        ))
        0 < payload <= 64 * 1024 || throw(ArgumentError(
            "real-time maximum_payload_bytes must be between one and 65536",
        ))
        kind == PhysicalHILTarget && overrun_policy == MeasureOnlyOverrun &&
            throw(ArgumentError("physical HIL cannot use measure-only overruns"))
        return new(
            id,
            kind,
            period,
            steps,
            overrun_policy,
            require_hard_realtime,
            scheduler_policy,
            priority,
            cpu,
            lock_memory,
            channels,
            payload,
        )
    end
end

struct RealtimeChannel
    id::Symbol
    direction::Symbol
    physical_unit::String
    scale::Float64
    offset::Float64
    raw_minimum::Float64
    raw_maximum::Float64
    safe_physical_value::Float64

    function RealtimeChannel(
        id::Symbol,
        direction::Symbol,
        physical_unit::AbstractString,
        scale::Real,
        offset::Real,
        raw_minimum::Real,
        raw_maximum::Real,
        safe_physical_value::Real,
    )
        id == Symbol("") && throw(ArgumentError("real-time channel id must not be empty"))
        direction in (:input, :output) || throw(ArgumentError(
            "real-time channel direction must be :input or :output",
        ))
        isempty(strip(physical_unit)) && throw(ArgumentError(
            "real-time channel physical_unit must not be empty",
        ))
        numeric = Float64.((scale, offset, raw_minimum, raw_maximum, safe_physical_value))
        all(isfinite, numeric) || throw(ArgumentError(
            "real-time channel numeric fields must be finite",
        ))
        numeric[1] != 0.0 || throw(ArgumentError(
            "real-time channel scale must be nonzero",
        ))
        numeric[3] < numeric[4] || throw(ArgumentError(
            "real-time channel raw range must be increasing",
        ))
        safe_raw = (numeric[5] - numeric[2]) / numeric[1]
        numeric[3] <= safe_raw <= numeric[4] || throw(ArgumentError(
            "real-time safe physical value is outside the declared raw range",
        ))
        return new(id, direction, String(physical_unit), numeric...)
    end
end

function physical_channel_value(channel::RealtimeChannel, raw_value::Real)
    raw = Float64(raw_value)
    isfinite(raw) || throw(ArgumentError("real-time raw channel value must be finite"))
    channel.raw_minimum <= raw <= channel.raw_maximum || throw(DomainError(
        raw,
        "real-time raw channel value is outside the declared range",
    ))
    return muladd(channel.scale, raw, channel.offset)
end

function raw_channel_value(channel::RealtimeChannel, physical_value::Real)
    physical = Float64(physical_value)
    isfinite(physical) || throw(ArgumentError(
        "real-time physical channel value must be finite",
    ))
    raw = (physical - channel.offset) / channel.scale
    channel.raw_minimum <= raw <= channel.raw_maximum || throw(DomainError(
        physical,
        "real-time physical channel value is outside the declared raw range",
    ))
    return raw
end

struct RealtimeUnavailableResult
    operation::Symbol
    target_id::Symbol
    target_kind::RealtimeTargetKind
    code::Symbol
    message::String
end

struct RealtimePreparation
    target::RealtimeTarget
    channels::Vector{RealtimeChannel}
    profile::Symbol
    hard_realtime::Bool
    physical_hardware::Bool
end

struct RealtimeExecutionMetadata
    target_id::Symbol
    target_kind::RealtimeTargetKind
    profile::Symbol
    configuration_sha256::String
    interface_kind::Symbol
    clock::Symbol
    scheduler_requested::Symbol
    scheduler_granted::Symbol
    cpu_requested::Union{Nothing,Int}
    cpu_granted::Union{Nothing,Int}
    memory_lock_requested::Bool
    memory_lock_granted::Bool
    channel_count::Int
    hard_realtime::Bool
    physical_hardware::Bool
end

function realtime_configuration_sha256(
    target::RealtimeTarget,
    channels::AbstractVector{RealtimeChannel},
)
    context = SHA.SHA256_CTX()
    SHA.update!(context, codeunits("aimora-realtime-configuration-v1\n"))
    for value in (
        target.id,
        target.kind,
        target.period_ns,
        target.step_count,
        target.overrun_policy,
        target.require_hard_realtime,
        target.scheduler_policy,
        target.scheduler_priority,
        target.cpu_index,
        target.lock_memory,
    )
        SHA.update!(context, codeunits(repr(value)))
        SHA.update!(context, UInt8[0x00])
    end
    for channel in channels
        for value in (
            channel.id,
            channel.direction,
            channel.physical_unit,
            channel.scale,
            channel.offset,
            channel.raw_minimum,
            channel.raw_maximum,
            channel.safe_physical_value,
        )
            SHA.update!(context, codeunits(repr(value)))
            SHA.update!(context, UInt8[0x00])
        end
    end
    return bytes2hex(SHA.digest!(context))
end

mutable struct RealtimeSafetyState
    latched::Bool
    failure_code::Symbol
    safe_output_values::Vector{Float64}
end

function RealtimeSafetyState(channels::AbstractVector{RealtimeChannel})
    outputs = filter(channel -> channel.direction == :output, channels)
    return RealtimeSafetyState(
        false,
        :none,
        getfield.(outputs, :safe_physical_value),
    )
end

function latch_realtime_safe_shutdown!(
    state::RealtimeSafetyState,
    failure_code::Symbol,
)
    failure_code == :none && throw(ArgumentError(
        "safe shutdown requires a nonempty failure code",
    ))
    state.latched = true
    state.failure_code = failure_code
    return state
end

function prepare_realtime_target(
    target::RealtimeTarget,
    channels::AbstractVector{RealtimeChannel},
)
    length(channels) <= target.maximum_channels || throw(ArgumentError(
        "real-time channel count exceeds the target limit",
    ))
    length(unique(getfield.(channels, :id))) == length(channels) ||
        throw(ArgumentError("real-time channel ids must be unique"))
    target.kind == PhysicalHILTarget && return RealtimeUnavailableResult(
        :prepare_realtime_target,
        target.id,
        target.kind,
        :physical_hil_unavailable,
        "Physical HIL requires exact hardware, driver, calibration, interlock, safe-state, and external evidence metadata.",
    )
    target.require_hard_realtime && return RealtimeUnavailableResult(
        :prepare_realtime_target,
        target.id,
        target.kind,
        :hard_realtime_unavailable,
        "This public checkout cannot establish an admitted PREEMPT_RT hard-real-time target.",
    )
    if target.scheduler_policy != :other || target.cpu_index !== nothing ||
            target.lock_memory
        return RealtimeUnavailableResult(
            :prepare_realtime_target,
            target.id,
            target.kind,
            :os_facility_unavailable,
            "Requested scheduler, affinity, or memory-lock facilities are unavailable until their granted state is applied and verified.",
        )
    end
    Sys.islinux() || return RealtimeUnavailableResult(
        :prepare_realtime_target,
        target.id,
        target.kind,
        :unsupported_operating_system,
        "The initial real-time software profile is available only on native Linux.",
    )
    target.kind == SharedLibraryControllerTarget && return RealtimeUnavailableResult(
        :prepare_realtime_target,
        target.id,
        target.kind,
        :controller_library_not_bound,
        "A shared-library target requires an exact ABI, content hash, and resolved AIMORA controller fixture.",
    )
    target.kind == LocalUDPControllerTarget && return RealtimeUnavailableResult(
        :prepare_realtime_target,
        target.id,
        target.kind,
        :local_udp_interface_not_bound,
        "A local UDP target requires an explicit bounded loopback interface.",
    )
    return RealtimePreparation(
        target,
        collect(channels),
        :linux_in_process_soft,
        false,
        false,
    )
end

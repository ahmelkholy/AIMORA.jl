export RealtimeExecutionFailure,
       RealtimeExecutionResult,
       run_paced_realtime!

struct RealtimeExecutionFailure
    code::Symbol
    logical_step::Int
    message::String
end

struct RealtimeExecutionResult
    preparation::RealtimePreparation
    metadata::RealtimeExecutionMetadata
    accepted_steps::Int
    timing_samples::Vector{RealtimeTimingSample}
    timing_summary::Union{Nothing,RealtimeTimingSummary}
    failure::Union{Nothing,RealtimeExecutionFailure}
    safety::RealtimeSafetyState
    final_input_values::Vector{Float64}
    final_output_values::Vector{Float64}
end

function _realtime_execution_failure(
    preparation::RealtimePreparation,
    metadata::RealtimeExecutionMetadata,
    accepted_steps::Int,
    samples::Vector{RealtimeTimingSample},
    populated_samples::Int,
    safety::RealtimeSafetyState,
    input_values::Vector{Float64},
    output_values::Vector{Float64},
    code::Symbol,
    logical_step::Int,
    message::AbstractString,
)
    populated = samples[1:populated_samples]
    summary = isempty(populated) ? nothing : summarize_realtime_timing(populated)
    return RealtimeExecutionResult(
        preparation,
        metadata,
        accepted_steps,
        populated,
        summary,
        RealtimeExecutionFailure(code, logical_step, String(message)),
        safety,
        copy(input_values),
        copy(output_values),
    )
end


function _realtime_execution_metadata(
    preparation::RealtimePreparation,
    interface::AbstractRealtimeInterface,
)
    target = preparation.target
    return RealtimeExecutionMetadata(
        target.id,
        target.kind,
        preparation.profile,
        realtime_configuration_sha256(target, preparation.channels),
        realtime_interface_kind(interface),
        :CLOCK_MONOTONIC,
        target.scheduler_policy,
        :other,
        target.cpu_index,
        nothing,
        target.lock_memory,
        false,
        length(preparation.channels),
        preparation.hard_realtime,
        preparation.physical_hardware,
    )
end

function run_paced_realtime!(
    model_step!::F,
    model_state,
    interface::I,
    preparation::RealtimePreparation;
    checkpoint,
    capture_state!::C,
    restore_state!::R,
    input_values::Vector{Float64},
    output_values::Vector{Float64},
    disable_gc::Bool=true,
) where {F,I<:AbstractRealtimeInterface,C,R}
    target = preparation.target
    target.kind in (
        NativeLinuxSoftwareTarget,
        LocalUDPControllerTarget,
        SharedLibraryControllerTarget,
    ) ||
        throw(ArgumentError(
            "paced public execution requires an admitted software interface target",
        ))
    length(input_values) == count(channel -> channel.direction == :input, preparation.channels) ||
        throw(DimensionMismatch("paced real-time input channel count mismatch"))
    length(output_values) == count(channel -> channel.direction == :output, preparation.channels) ||
        throw(DimensionMismatch("paced real-time output channel count mismatch"))
    all(isfinite, input_values) && all(isfinite, output_values) || throw(ArgumentError(
        "paced real-time initial frames must be finite",
    ))
    samples = Vector{RealtimeTimingSample}(undef, target.step_count)
    safety = RealtimeSafetyState(preparation.channels)
    metadata = _realtime_execution_metadata(preparation, interface)
    prepare_realtime_interface!(interface)
    interface_checkpoint = realtime_interface_checkpoint(interface)
    period_ns = target.period_ns
    period_s = period_ns / 1.0e9
    epoch_ns = Base.Checked.checked_add(monotonic_time_ns(), period_ns)
    accepted_steps = 0
    populated_samples = 0
    old_gc = GC.enable(!disable_gc)
    try
        for logical_step in 0:(target.step_count - 1)
            release_ns = realtime_release_ns(epoch_ns, period_ns, logical_step)
            wait_until_monotonic_ns(release_ns)
            capture_state!(checkpoint, model_state)
            capture_realtime_interface_state!(interface_checkpoint, interface)
            start_ns = monotonic_time_ns()
            try
                model_step!(
                    model_state,
                    input_values,
                    output_values,
                    logical_step,
                    logical_step * period_s,
                    period_s,
                )
                all(isfinite, output_values) || throw(ArgumentError(
                    "paced real-time model produced nonfinite output",
                ))
                exchange_realtime_interface!(
                    interface,
                    input_values,
                    output_values,
                    logical_step,
                    release_ns,
                )
            catch error
                restore_state!(model_state, checkpoint)
                restore_realtime_interface_state!(interface, interface_checkpoint)
                latch_realtime_safe_shutdown!(safety, :execution_failure)
                safe_shutdown_realtime_interface!(interface)
                return _realtime_execution_failure(
                    preparation,
                    metadata,
                    accepted_steps,
                    samples,
                    populated_samples,
                    safety,
                    input_values,
                    output_values,
                    :execution_failure,
                    logical_step,
                    sprint(showerror, error),
                )
            end
            completion_ns = monotonic_time_ns()
            sample = realtime_timing_sample(
                logical_step,
                release_ns,
                start_ns,
                completion_ns,
                period_ns,
            )
            populated_samples += 1
            samples[populated_samples] = sample
            if sample.overrun && target.overrun_policy == FailSafeOnOverrun
                restore_state!(model_state, checkpoint)
                restore_realtime_interface_state!(interface, interface_checkpoint)
                latch_realtime_safe_shutdown!(safety, :deadline_overrun)
                safe_shutdown_realtime_interface!(interface)
                return _realtime_execution_failure(
                    preparation,
                    metadata,
                    accepted_steps,
                    samples,
                    populated_samples,
                    safety,
                    input_values,
                    output_values,
                    :deadline_overrun,
                    logical_step,
                    "paced real-time response exceeded the declared period",
                )
            end
            accepted_steps += 1
        end
    finally
        GC.enable(old_gc)
    end
    populated = samples[1:populated_samples]
    return RealtimeExecutionResult(
        preparation,
        metadata,
        accepted_steps,
        populated,
        summarize_realtime_timing(populated),
        nothing,
        safety,
        copy(input_values),
        copy(output_values),
    )
end

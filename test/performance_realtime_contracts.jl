@testset "performance execution public contracts" begin
    performance = AIMORA.PerformanceExecution
    request = performance.PerformanceExecutionRequest(
        performance.ThreadedBatchExecution;
        job_count=8,
        worker_count=4,
        memory_limit_bytes=14 * 1024^3,
    )
    @test request.deterministic
    @test request.worker_count == 4
    @test !request.allow_pre_execution_cpu_fallback
    preparation = performance.PerformancePreparation(
        request,
        performance.ThreadedBatchExecution,
        :test_backend,
        nothing,
    )
    report = performance.PerformanceExecutionReport(
        preparation;
        complete_results=8,
        result_signature_sha256=repeat("a", 64),
        exact_discrete_identity=true,
        maximum_scaled_residual=1.0e-14,
        maximum_scaled_error=2.0e-14,
        setup_seconds=0.1,
        execution_seconds=0.2,
        output_seconds=0.01,
        retained_bytes=1024,
        speedup=2.0,
        parallel_efficiency=0.5,
        crossover_jobs=8,
    )
    @test report.complete_results == request.job_count
    @test report.exact_discrete_identity
    @test_throws ArgumentError performance.PerformanceExecutionReport(
        preparation;
        complete_results=7,
        result_signature_sha256=repeat("a", 64),
        exact_discrete_identity=true,
        maximum_scaled_residual=0.0,
        maximum_scaled_error=0.0,
        setup_seconds=0.0,
        execution_seconds=0.0,
        output_seconds=0.0,
        retained_bytes=0,
    )
    @test_throws ArgumentError performance.PerformanceExecutionRequest(
        performance.SerialCPUExecution;
        job_count=2,
        worker_count=2,
    )
    @test_throws ArgumentError performance.PerformanceExecutionRequest(
        performance.ThreadedBatchExecution;
        job_count=2,
        worker_count=2,
        deterministic=false,
    )

    if !AIMORA.solver_available()
        capabilities = performance.performance_capabilities()
        @test length(capabilities) == length(instances(performance.PerformanceExecutionMode))
        @test all(!capability.available for capability in capabilities)
        @test all(capability.backend == :open_core for capability in capabilities)
        unavailable = performance.prepare_performance_execution(request)
        @test unavailable isa performance.PerformanceUnavailableResult
        @test unavailable.code == :solver_backend_unavailable
        @test unavailable.requested_mode == performance.ThreadedBatchExecution
        @test !unavailable.cpu_fallback_permitted
    end
end

@testset "real-time target, channel, timing, and refusal contracts" begin
    realtime = AIMORA.RealtimeLoop
    input_channel = realtime.RealtimeChannel(
        :measured_voltage,
        :input,
        "V",
        2.0,
        -1.0,
        -10.0,
        10.0,
        0.0,
    )
    output_channel = realtime.RealtimeChannel(
        :control_command,
        :output,
        "1",
        0.5,
        0.0,
        -2.0,
        2.0,
        0.0,
    )
    @test realtime.physical_channel_value(input_channel, 3.0) == 5.0
    @test realtime.raw_channel_value(input_channel, 5.0) == 3.0
    @test_throws DomainError realtime.physical_channel_value(input_channel, 11.0)
    @test_throws ArgumentError realtime.RealtimeChannel(
        :invalid_scale,
        :input,
        "V",
        0.0,
        0.0,
        -1.0,
        1.0,
        0.0,
    )

    software_target = realtime.RealtimeTarget(
        :linux_software_loop,
        realtime.NativeLinuxSoftwareTarget;
        period_ns=20_000,
        step_count=2_500,
        overrun_policy=realtime.MeasureOnlyOverrun,
    )
    preparation = realtime.prepare_realtime_target(
        software_target,
        [input_channel, output_channel],
    )
    @test preparation isa realtime.RealtimePreparation
    @test preparation.profile == :linux_in_process_soft
    @test !preparation.hard_realtime
    @test !preparation.physical_hardware
    @test_throws ArgumentError realtime.prepare_realtime_target(
        software_target,
        [input_channel, input_channel],
    )

    hard_target = realtime.RealtimeTarget(
        :requested_hard_target,
        realtime.NativeLinuxSoftwareTarget;
        period_ns=100_000,
        step_count=10,
        require_hard_realtime=true,
        scheduler_policy=:fifo,
        scheduler_priority=50,
        lock_memory=true,
    )
    hard_unavailable = realtime.prepare_realtime_target(hard_target, [output_channel])
    @test hard_unavailable isa realtime.RealtimeUnavailableResult
    @test hard_unavailable.code == :hard_realtime_unavailable
    facility_target = realtime.RealtimeTarget(
        :requested_unverified_facilities,
        realtime.NativeLinuxSoftwareTarget;
        period_ns=100_000,
        step_count=10,
        scheduler_policy=:fifo,
        scheduler_priority=20,
        cpu_index=1,
        lock_memory=true,
    )
    facility_unavailable = realtime.prepare_realtime_target(
        facility_target,
        [output_channel],
    )
    @test facility_unavailable isa realtime.RealtimeUnavailableResult
    @test facility_unavailable.code == :os_facility_unavailable

    physical_target = realtime.RealtimeTarget(
        :physical_controller_hil,
        realtime.PhysicalHILTarget;
        period_ns=50_000,
        step_count=100,
    )
    physical_unavailable = realtime.prepare_realtime_target(
        physical_target,
        [output_channel],
    )
    @test physical_unavailable isa realtime.RealtimeUnavailableResult
    @test physical_unavailable.code == :physical_hil_unavailable
    @test_throws ArgumentError realtime.RealtimeTarget(
        :unsafe_physical_measurement,
        realtime.PhysicalHILTarget;
        period_ns=50_000,
        step_count=100,
        overrun_policy=realtime.MeasureOnlyOverrun,
    )

    @test realtime.realtime_release_ns(1_000, 1_000, 3) == 4_000
    if Sys.islinux()
        clock_start = realtime.monotonic_time_ns()
        clock_target = clock_start + 1_000_000
        clock_wake = realtime.wait_until_monotonic_ns(clock_target)
        @test clock_wake >= clock_target
    end
    timing_samples = [
        realtime.realtime_timing_sample(0, 1_000, 1_000, 1_500, 1_000),
        realtime.realtime_timing_sample(1, 2_000, 2_100, 2_900, 1_000),
        realtime.realtime_timing_sample(2, 3_000, 3_200, 4_300, 1_000),
    ]
    @test getfield.(timing_samples, :response_ns) == [500, 900, 1_300]
    @test getfield.(timing_samples, :slack_ns) == [500, 100, -300]
    @test getfield.(timing_samples, :overrun) == [false, false, true]
    timing_summary = realtime.summarize_realtime_timing(timing_samples)
    @test timing_summary.samples == 3
    @test timing_summary.maximum_response_ns == 1_300
    @test timing_summary.maximum_jitter_ns == 200
    @test timing_summary.overruns == 1
    @test_throws ArgumentError realtime.realtime_timing_sample(
        0,
        1_000,
        999,
        1_500,
        1_000,
    )

    safety = realtime.RealtimeSafetyState([input_channel, output_channel])
    @test !safety.latched
    @test safety.safe_output_values == [0.0]
    @test realtime.latch_realtime_safe_shutdown!(safety, :deadline_overrun) === safety
    @test safety.latched
    @test safety.failure_code == :deadline_overrun

    controller! = function (
        destination_inputs,
        source_outputs,
        _controller_state,
        _logical_step,
        _logical_time_ns,
    )
        destination_inputs[1] = 1.0 - 0.1 * source_outputs[1]
        return nothing
    end
    model_step! = function (
        state,
        inputs,
        outputs,
        _logical_step,
        _logical_time_s,
        step_s,
    )
        state[1] += inputs[1] * step_s
        outputs[1] = state[1]
        return nothing
    end
    capture_state! = (checkpoint, state) -> copyto!(checkpoint, state)
    restore_state! = (state, checkpoint) -> copyto!(state, checkpoint)
    paced_target = realtime.RealtimeTarget(
        :paced_in_process,
        realtime.NativeLinuxSoftwareTarget;
        period_ns=10_000_000,
        step_count=5,
        overrun_policy=realtime.FailSafeOnOverrun,
    )
    paced_preparation = realtime.prepare_realtime_target(
        paced_target,
        [input_channel, output_channel],
    )
    paced_interface = realtime.InProcessControllerInterface(
        controller!,
        nothing,
        paced_preparation.channels,
    )
    paced_state = [0.0]
    paced_result = realtime.run_paced_realtime!(
        model_step!,
        paced_state,
        paced_interface,
        paced_preparation;
        checkpoint=[0.0],
        capture_state!,
        restore_state!,
        input_values=[1.0],
        output_values=[0.0],
    )
    @test paced_result.failure === nothing
    @test paced_result.accepted_steps == 5
    @test length(paced_result.timing_samples) == 5
    @test paced_result.timing_summary.overruns == 0
    @test paced_interface.sequence == 5
    @test paced_result.metadata.interface_kind == :in_process
    @test paced_result.metadata.clock == :CLOCK_MONOTONIC
    @test occursin(r"^[0-9a-f]{64}$", paced_result.metadata.configuration_sha256)
    @test !paced_result.metadata.hard_realtime
    @test !paced_result.metadata.physical_hardware

    udp_target = realtime.RealtimeTarget(
        :paced_local_udp,
        realtime.LocalUDPControllerTarget;
        period_ns=10_000_000,
        step_count=5,
        overrun_policy=realtime.FailSafeOnOverrun,
    )
    udp_unbound = realtime.prepare_realtime_target(
        udp_target,
        [input_channel, output_channel],
    )
    @test udp_unbound isa realtime.RealtimeUnavailableResult
    @test udp_unbound.code == :local_udp_interface_not_bound
    udp_controller! = function (
        destination_inputs,
        source_outputs,
        controller_state,
        _logical_step,
        _logical_time_ns,
    )
        destination_inputs[1] = 1.0 - 0.1 * source_outputs[1]
        return controller_state + 1
    end
    udp_interface = realtime.open_local_udp_controller(
        udp_controller!,
        0,
        [input_channel, output_channel],
    )
    udp_preparation = realtime.prepare_realtime_target(
        udp_target,
        [input_channel, output_channel],
        udp_interface,
    )
    @test udp_preparation isa realtime.RealtimePreparation
    @test udp_preparation.profile == :linux_local_udp_soft
    udp_state = [0.0]
    udp_result = realtime.run_paced_realtime!(
        model_step!,
        udp_state,
        udp_interface,
        udp_preparation;
        checkpoint=[0.0],
        capture_state!,
        restore_state!,
        input_values=[1.0],
        output_values=[0.0],
    )
    @test udp_result.failure === nothing
    @test udp_result.accepted_steps == 5
    @test udp_result.timing_summary.overruns == 0
    @test udp_interface.sequence == 5
    @test udp_interface.controller_state == 5
    @test realtime.close_realtime_interface!(udp_interface) === udp_interface
    @test udp_interface.closed

    slow_model_step! = function (
        state,
        inputs,
        outputs,
        logical_step,
        logical_time_s,
        step_s,
    )
        model_step!(state, inputs, outputs, logical_step, logical_time_s, step_s)
        realtime.wait_until_monotonic_ns(realtime.monotonic_time_ns() + 2_000_000)
        return nothing
    end
    overrun_target = realtime.RealtimeTarget(
        :fail_safe_overrun,
        realtime.NativeLinuxSoftwareTarget;
        period_ns=1_000_000,
        step_count=2,
    )
    overrun_preparation = realtime.prepare_realtime_target(
        overrun_target,
        [input_channel, output_channel],
    )
    overrun_interface = realtime.InProcessControllerInterface(
        controller!,
        nothing,
        overrun_preparation.channels,
    )
    overrun_state = [0.0]
    overrun_result = realtime.run_paced_realtime!(
        slow_model_step!,
        overrun_state,
        overrun_interface,
        overrun_preparation;
        checkpoint=[0.0],
        capture_state!,
        restore_state!,
        input_values=[1.0],
        output_values=[0.0],
    )
    @test overrun_result.failure.code == :deadline_overrun
    @test overrun_result.accepted_steps == 0
    @test overrun_state == [0.0]
    @test overrun_result.safety.latched
    @test overrun_interface.output_values == [0.0]
end

@testset "public isolated study worker contracts" begin
    limits = AIMORA.StudyWorkerLimits(
        maximum_jobs_per_process=8,
        maximum_recoveries=2,
        cancellation_poll_interval_steps=32,
        maximum_binary_window_bytes=4096,
        maximum_result_windows=4,
        heartbeat_interval_ms=100,
        idle_timeout_ms=1000,
    )
    manifest = AIMORA.StudyWorkerManifest(
        AIMORA.STUDY_WORKER_PROTOCOL_VERSION,
        v"0.1.0",
        :aimora_production,
        (
            AIMORA.StudyWorkerCapability(:result_window_descriptor, v"1.0.0"),
            AIMORA.StudyWorkerCapability(:job_lifecycle, v"1.0.0"),
            AIMORA.StudyWorkerCapability(:request_cancel, v"1.0.0"),
        ),
        limits,
    )
    record = AIMORA.study_worker_manifest_record(manifest)
    @test record.protocol_version == "1.0.0"
    @test record.backend_id == "aimora_production"
    @test getfield.(record.capabilities, :id) == (
        "job_lifecycle",
        "request_cancel",
        "result_window_descriptor",
    )
    @test record.limits.maximum_binary_window_bytes == 4096
    public_text = sprint(show, record)
    @test !occursin("AIMORASolvers", public_text)
    @test !occursin("ProductionBackend", public_text)
    @test !occursin(homedir(), public_text)
    @test_throws ArgumentError AIMORA.StudyWorkerManifest(
        v"2.0.0",
        v"0.1.0",
        :aimora_production,
        (),
        limits,
    )
    @test_throws ArgumentError AIMORA.StudyWorkerCapability(
        "../private",
        v"1.0.0",
    )

    token = AIMORA.StudyWorkerCancellation()
    @test !AIMORA.study_worker_cancellation_requested(token)
    @test AIMORA.request_study_worker_cancellation!(token, "operator request")
    @test !AIMORA.request_study_worker_cancellation!(token, "duplicate")
    @test AIMORA.study_worker_cancellation_requested(token)
    @test AIMORA.study_worker_cancellation_reason(token) == "operator request"
    cancelled = try
        AIMORA.check_study_worker_cancellation(token)
        nothing
    catch failure
        failure
    end
    @test cancelled isa AIMORA.StudyWorkerCancelled
    @test cancelled.reason == "operator request"
    @test_throws ArgumentError AIMORA.request_study_worker_cancellation!(
        AIMORA.StudyWorkerCancellation(),
        "line one\nline two",
    )

    failure = AIMORA.StudyWorkerFailure(
        :resource_limit,
        :result_window,
        false,
        "bounded result window limit reached",
    )
    @test failure.code == :resource_limit
    @test failure.phase == :result_window
    @test !failure.retryable
    @test_throws ArgumentError AIMORA.StudyWorkerFailure(
        :private_solver_failure,
        :execute,
        false,
        "not public",
    )
    @test_throws ArgumentError AIMORA.StudyWorkerFailure(
        :worker_failed,
        :execute,
        false,
        "failure at /private/solver/state",
    )

    window = AIMORA.StudyWorkerResultWindow(
        "artifact-7",
        1024,
        2048,
        8192,
        "application/octet-stream",
        "float64",
        (256,),
        repeat("a", 64);
        maximum_window_bytes=4096,
    )
    window_record = AIMORA.study_worker_result_window_record(window)
    @test window_record.offset_bytes == 1024
    @test window_record.window_bytes == 2048
    @test window_record.shape == (256,)
    @test !hasproperty(window_record, :path)
    @test_throws ArgumentError AIMORA.StudyWorkerResultWindow(
        "../artifact",
        0,
        1,
        1,
        "application/octet-stream",
        "uint8",
        (1,),
        repeat("0", 64),
    )
    @test_throws ArgumentError AIMORA.StudyWorkerResultWindow(
        "artifact-8",
        0,
        4097,
        8192,
        "application/octet-stream",
        "uint8",
        (4097,),
        repeat("0", 64);
        maximum_window_bytes=4096,
    )

    if !AIMORA.solver_available()
        unavailable = AIMORA.study_worker_manifest()
        @test unavailable isa AIMORA.SolverUnavailableResult
        @test unavailable.required_capability == :study_worker_lifecycle
    end
end

module PerformanceExecution

export PerformanceExecutionMode,
       SerialCPUExecution,
       SparseCPUExecution,
       ThreadedBatchExecution,
       DistributedBatchExecution,
       NativeCUDABatchExecution,
       PerformanceExecutionRequest,
       PerformanceCapability,
       PerformanceUnavailableResult,
       PerformancePreparation,
       PerformanceExecutionReport,
       performance_capabilities,
       performance_mode_available,
       prepare_performance_execution

using ..AIMORA: AbstractAIMORASolverBackend, active_solver_backend

@enum PerformanceExecutionMode begin
    SerialCPUExecution
    SparseCPUExecution
    ThreadedBatchExecution
    DistributedBatchExecution
    NativeCUDABatchExecution
end

struct PerformanceExecutionRequest
    mode::PerformanceExecutionMode
    job_count::Int
    worker_count::Int
    memory_limit_bytes::Int
    deterministic::Bool
    allow_pre_execution_cpu_fallback::Bool

    function PerformanceExecutionRequest(
        mode::PerformanceExecutionMode;
        job_count::Integer=1,
        worker_count::Integer=1,
        memory_limit_bytes::Integer=14 * 1024^3,
        deterministic::Bool=true,
        allow_pre_execution_cpu_fallback::Bool=false,
    )
        jobs = Int(job_count)
        workers = Int(worker_count)
        memory = Int(memory_limit_bytes)
        jobs > 0 || throw(ArgumentError("performance job_count must be positive"))
        workers > 0 || throw(ArgumentError("performance worker_count must be positive"))
        workers <= jobs || throw(ArgumentError(
            "performance worker_count cannot exceed job_count",
        ))
        memory > 0 || throw(ArgumentError(
            "performance memory_limit_bytes must be positive",
        ))
        deterministic || throw(ArgumentError(
            "AIMORA performance execution must preserve deterministic collection",
        ))
        mode in (SerialCPUExecution, SparseCPUExecution) && workers != 1 &&
            throw(ArgumentError("serial and sparse CPU requests require one worker"))
        return new(
            mode,
            jobs,
            workers,
            memory,
            deterministic,
            allow_pre_execution_cpu_fallback,
        )
    end
end

struct PerformanceCapability
    mode::PerformanceExecutionMode
    available::Bool
    deterministic::Bool
    maximum_workers::Int
    minimum_jobs::Int
    backend::Symbol
    reason::String

    function PerformanceCapability(
        mode::PerformanceExecutionMode,
        available::Bool,
        deterministic::Bool,
        maximum_workers::Integer,
        backend::Symbol,
        reason::AbstractString,
        minimum_jobs::Integer=1,
    )
        workers = Int(maximum_workers)
        jobs = Int(minimum_jobs)
        workers >= 0 || throw(ArgumentError(
            "performance capability maximum_workers must be nonnegative",
        ))
        available && workers == 0 && throw(ArgumentError(
            "available performance capability requires at least one worker",
        ))
        jobs > 0 || throw(ArgumentError(
            "performance capability minimum_jobs must be positive",
        ))
        backend == Symbol("") && throw(ArgumentError(
            "performance capability backend must not be empty",
        ))
        isempty(strip(reason)) && throw(ArgumentError(
            "performance capability reason must not be empty",
        ))
        return new(
            mode,
            available,
            deterministic,
            workers,
            jobs,
            backend,
            String(reason),
        )
    end
end

struct PerformanceUnavailableResult
    operation::Symbol
    requested_mode::PerformanceExecutionMode
    code::Symbol
    message::String
    cpu_fallback_permitted::Bool
end

struct PerformancePreparation
    request::PerformanceExecutionRequest
    selected_mode::PerformanceExecutionMode
    backend::Symbol
    fallback_reason::Union{Nothing,String}
end

struct PerformanceExecutionReport
    preparation::PerformancePreparation
    complete_results::Int
    result_signature_sha256::String
    exact_discrete_identity::Bool
    maximum_scaled_residual::Float64
    maximum_scaled_error::Float64
    setup_seconds::Float64
    execution_seconds::Float64
    output_seconds::Float64
    retained_bytes::Int
    speedup::Union{Nothing,Float64}
    parallel_efficiency::Union{Nothing,Float64}
    crossover_jobs::Union{Nothing,Int}

    function PerformanceExecutionReport(
        preparation::PerformancePreparation;
        complete_results::Integer,
        result_signature_sha256::AbstractString,
        exact_discrete_identity::Bool,
        maximum_scaled_residual::Real,
        maximum_scaled_error::Real,
        setup_seconds::Real,
        execution_seconds::Real,
        output_seconds::Real,
        retained_bytes::Integer,
        speedup::Union{Nothing,Real}=nothing,
        parallel_efficiency::Union{Nothing,Real}=nothing,
        crossover_jobs::Union{Nothing,Integer}=nothing,
    )
        results = Int(complete_results)
        signature = lowercase(String(result_signature_sha256))
        residual = Float64(maximum_scaled_residual)
        error = Float64(maximum_scaled_error)
        timings = Float64.((setup_seconds, execution_seconds, output_seconds))
        bytes = Int(retained_bytes)
        measured_speedup = speedup === nothing ? nothing : Float64(speedup)
        efficiency = parallel_efficiency === nothing ? nothing :
            Float64(parallel_efficiency)
        crossover = crossover_jobs === nothing ? nothing : Int(crossover_jobs)
        results == preparation.request.job_count || throw(ArgumentError(
            "performance report must contain every requested result",
        ))
        occursin(r"^[0-9a-f]{64}$", signature) || throw(ArgumentError(
            "performance result signature must be 64 lowercase hex characters",
        ))
        exact_discrete_identity || throw(ArgumentError(
            "performance report requires exact discrete identity",
        ))
        isfinite(residual) && residual >= 0.0 || throw(ArgumentError(
            "performance maximum residual must be finite and nonnegative",
        ))
        isfinite(error) && error >= 0.0 || throw(ArgumentError(
            "performance maximum error must be finite and nonnegative",
        ))
        all(value -> isfinite(value) && value >= 0.0, timings) ||
            throw(ArgumentError(
                "performance phase timings must be finite and nonnegative",
            ))
        bytes >= 0 || throw(ArgumentError(
            "performance retained bytes must be nonnegative",
        ))
        measured_speedup === nothing ||
            isfinite(measured_speedup) && measured_speedup > 0.0 ||
            throw(ArgumentError("performance speedup must be finite and positive"))
        efficiency === nothing || isfinite(efficiency) && efficiency > 0.0 ||
            throw(ArgumentError(
                "performance parallel efficiency must be finite and positive",
            ))
        crossover === nothing || crossover > 0 || throw(ArgumentError(
            "performance crossover job count must be positive",
        ))
        return new(
            preparation,
            results,
            signature,
            exact_discrete_identity,
            residual,
            error,
            timings...,
            bytes,
            measured_speedup,
            efficiency,
            crossover,
        )
    end
end

performance_capabilities(::AbstractAIMORASolverBackend) = PerformanceCapability[]

function performance_capabilities()
    backend = active_solver_backend()
    if backend === nothing
        reason = "No AIMORA production solver backend is active."
        return [
            PerformanceCapability(mode, false, true, 0, :open_core, reason)
            for mode in instances(PerformanceExecutionMode)
        ]
    end
    capabilities = performance_capabilities(backend)
    all(capability -> capability isa PerformanceCapability, capabilities) ||
        throw(ArgumentError(
            "performance backend capabilities must use PerformanceCapability",
        ))
    length(unique(getfield.(capabilities, :mode))) == length(capabilities) ||
        throw(ArgumentError("performance backend repeats an execution mode"))
    return capabilities
end

function performance_mode_available(mode::PerformanceExecutionMode)
    capabilities = performance_capabilities()
    index = findfirst(capability -> capability.mode == mode, capabilities)
    return index !== nothing && capabilities[index].available
end

function prepare_performance_execution(
    ::AbstractAIMORASolverBackend,
    request::PerformanceExecutionRequest,
)
    return PerformanceUnavailableResult(
        :prepare_performance_execution,
        request.mode,
        :unsupported_performance_mode,
        "The active solver backend does not implement the requested performance mode.",
        request.allow_pre_execution_cpu_fallback,
    )
end


function prepare_performance_execution(request::PerformanceExecutionRequest)
    backend = active_solver_backend()
    backend === nothing && return PerformanceUnavailableResult(
        :prepare_performance_execution,
        request.mode,
        :solver_backend_unavailable,
        "No AIMORA production solver backend is active.",
        request.allow_pre_execution_cpu_fallback,
    )
    return prepare_performance_execution(backend, request)
end

end

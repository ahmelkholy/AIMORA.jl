export AbstractAIMORASolverBackend,
       SolverCapability,
       SolverUnavailableResult,
       activate_solver!,
       active_solver_backend,
       advance_converter_system!,
       backend_metadata,
       advance_partitioned_emt!,
       execute_study!,
       execute_partitioned_emt!,
       execute_converter_system!,
       materialize_emt_breaker_poles,
       materialize_measurement_branches,
       prepare_protection_task_pipeline,
       advance_protection_task_pipeline!,
       partitioned_emt_checkpoint,
       partitioned_emt_status,
       prepare_line_fit,
       prepare_partitioned_emt,
       prepare_converter_system,
       prepare_study,
       require_solver,
       restore_partitioned_emt_checkpoint!,
       restore_backend_state!,
       snapshot_backend_state,
       solver_available,
       solver_capabilities,
       solver_status

"""Public contract implemented by an explicitly loaded AIMORA numerical backend."""
abstract type AbstractAIMORASolverBackend end

"""One versioned study capability advertised by a registered solver backend."""
struct SolverCapability
    id::Symbol
    representation::Symbol
    fidelity::Symbol
    version::VersionNumber

    function SolverCapability(
        id::Symbol,
        representation::Symbol,
        fidelity::Symbol,
        version::VersionNumber,
    )
        id == Symbol("") && throw(ArgumentError("solver capability id must not be empty"))
        representation == Symbol("") && throw(ArgumentError(
            "solver capability representation must not be empty",
        ))
        fidelity == Symbol("") && throw(ArgumentError(
            "solver capability fidelity must not be empty",
        ))
        return new(id, representation, fidelity, version)
    end
end

"""Typed result returned when a requested numerical backend operation is unavailable."""
struct SolverUnavailableResult
    operation::Symbol
    required_capability::Symbol
    mode::Symbol
    message::String
end

const _ACTIVE_SOLVER_BACKEND = Ref{Union{Nothing,AbstractAIMORASolverBackend}}(nothing)
const _SOLVER_ACTIVATION_LOCK = ReentrantLock()

function _solver_unavailable_result(
    operation::Symbol,
    required_capability::Symbol = :production_solve;
    message::AbstractString = "No AIMORA production solver backend is active.",
)
    return SolverUnavailableResult(
        operation,
        required_capability,
        :open_core,
        String(message),
    )
end

"""Return the explicitly activated solver backend, or `nothing` in the open core."""
active_solver_backend() = _ACTIVE_SOLVER_BACKEND[]

"""Return `true` only after a backend has completed explicit runtime activation."""
solver_available() = active_solver_backend() !== nothing

"""Capabilities implemented by a backend. Backend packages extend this method."""
solver_capabilities(::AbstractAIMORASolverBackend) = SolverCapability[]

"""Nonsecret backend identity and version metadata. Backend packages extend this method."""
backend_metadata(::AbstractAIMORASolverBackend) = (
    name = :unregistered,
    version = v"0.0.0",
)

function _install_solver_runtime!(::AbstractAIMORASolverBackend)
    throw(ArgumentError(
        "the requested solver backend does not implement AIMORA runtime activation",
    ))
end

"""Install and activate one backend explicitly for the lifetime of this Julia process."""
function activate_solver!(backend::AbstractAIMORASolverBackend)
    lock(_SOLVER_ACTIVATION_LOCK) do
        active = active_solver_backend()
        if active !== nothing
            typeof(active) === typeof(backend) || throw(ArgumentError(
                "a different AIMORA solver backend is already active in this process",
            ))
            return active
        end
        capabilities = solver_capabilities(backend)
        all(capability -> capability isa SolverCapability, capabilities) ||
            throw(ArgumentError("solver backend capabilities must use SolverCapability"))
        _install_solver_runtime!(backend)
        _ACTIVE_SOLVER_BACKEND[] = backend
        return backend
    end
end

"""Report the backward-compatible open-core or full-engine solver state."""
function solver_status()
    backend = active_solver_backend()
    available = backend !== nothing
    return (
        available = available,
        mode = available ? :full_engine : :open_core,
        source = available ? :licensed_component : :not_installed,
        backend = available ? backend_metadata(backend) : nothing,
    )
end

"""Return `nothing` for an active backend or throw the established installation error."""
function require_solver()
    solver_available() && return nothing
    error(
        "The licensed AIMORA numerical component is not active. " *
        "Use the installation and activation instructions supplied with your distribution.",
    )
end

"""Return the active backend capabilities without exposing implementation details."""
function solver_capabilities()
    backend = active_solver_backend()
    return backend === nothing ? SolverCapability[] : solver_capabilities(backend)
end

"""Return active backend metadata or a typed unavailable result."""
function backend_metadata()
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:backend_metadata) :
           backend_metadata(backend)
end

prepare_study(::AbstractAIMORASolverBackend, project, study) =
    _solver_unavailable_result(
        :prepare_study,
        :study_preparation;
        message = "The active backend does not implement this study preparation contract.",
    )

execute_study!(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :execute_study,
        :study_execution;
        message = "The active backend does not implement this prepared study contract.",
    )

materialize_measurement_branches(::AbstractAIMORASolverBackend, definitions) =
    _solver_unavailable_result(
        :materialize_measurement_branches,
        :measurement_network_materialization;
        message = "The active backend does not implement measurement-network branches.",
    )

materialize_emt_breaker_poles(
    ::AbstractAIMORASolverBackend,
    runtime,
    specification,
    terminal_nodes,
) = _solver_unavailable_result(
    :materialize_emt_breaker_poles,
    :emt_protection_breaker;
    message = "The active backend does not implement physical EMT breaker poles.",
)

prepare_protection_task_pipeline(::AbstractAIMORASolverBackend, pipeline) =
    _solver_unavailable_result(
        :prepare_protection_task_pipeline,
        :emt_protection_breaker;
        message = "The active backend does not implement exact protection-task dispatch.",
    )

advance_protection_task_pipeline!(::AbstractAIMORASolverBackend, prepared, instant) =
    _solver_unavailable_result(
        :advance_protection_task_pipeline,
        :emt_protection_breaker;
        message = "The active backend does not implement exact protection-task dispatch.",
    )

prepare_line_fit(::AbstractAIMORASolverBackend, request) =
    _solver_unavailable_result(
        :prepare_line_fit,
        :coupled_line_fitting;
        message = "The active backend does not implement coupled line fitting.",
    )

snapshot_backend_state(::AbstractAIMORASolverBackend, runtime) =
    _solver_unavailable_result(
        :snapshot_backend_state,
        :backend_snapshot;
        message = "The active backend does not implement runtime snapshots.",
    )

restore_backend_state!(::AbstractAIMORASolverBackend, runtime, snapshot) =
    _solver_unavailable_result(
        :restore_backend_state,
        :backend_snapshot;
        message = "The active backend does not implement runtime restoration.",
    )

prepare_partitioned_emt(::AbstractAIMORASolverBackend, study) =
    _solver_unavailable_result(
        :prepare_partitioned_emt,
        :local_multirate_partitioned_emt;
        message = "The active backend does not implement partitioned EMT preparation.",
    )

advance_partitioned_emt!(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :advance_partitioned_emt,
        :local_multirate_partitioned_emt;
        message = "The active backend does not implement partitioned EMT advancement.",
    )

execute_partitioned_emt!(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :execute_partitioned_emt,
        :local_multirate_partitioned_emt;
        message = "The active backend does not implement partitioned EMT execution.",
    )

partitioned_emt_checkpoint(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :partitioned_emt_checkpoint,
        :local_multirate_partitioned_emt;
        message = "The active backend does not implement partitioned EMT checkpoints.",
    )

restore_partitioned_emt_checkpoint!(
    ::AbstractAIMORASolverBackend,
    prepared,
    checkpoint,
) = _solver_unavailable_result(
    :restore_partitioned_emt_checkpoint,
    :local_multirate_partitioned_emt;
    message = "The active backend does not implement partitioned EMT restoration.",
)

partitioned_emt_status(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :partitioned_emt_status,
        :local_multirate_partitioned_emt;
        message = "The active backend does not implement partitioned EMT status.",
    )

prepare_converter_system(::AbstractAIMORASolverBackend, study) =
    _solver_unavailable_result(
        :prepare_converter_system,
        :extended_converter_systems;
        message = "The active backend does not implement converter-system preparation.",
    )

advance_converter_system!(
    ::AbstractAIMORASolverBackend,
    prepared,
    accepted_step_count::Integer=1,
) = _solver_unavailable_result(
    :advance_converter_system,
    :extended_converter_systems;
    message = "The active backend does not implement bounded converter-system advancement.",
)

execute_converter_system!(::AbstractAIMORASolverBackend, prepared) =
    _solver_unavailable_result(
        :execute_converter_system,
        :extended_converter_systems;
        message = "The active backend does not implement converter-system execution.",
    )

function prepare_study(project, study)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:prepare_study, :study_preparation) :
           prepare_study(backend, project, study)
end

function execute_study!(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:execute_study, :study_execution) :
           execute_study!(backend, prepared)
end

function materialize_measurement_branches(definitions)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :materialize_measurement_branches,
        :measurement_network_materialization,
    ) : materialize_measurement_branches(backend, definitions)
end

function materialize_emt_breaker_poles(runtime, specification, terminal_nodes)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :materialize_emt_breaker_poles,
        :emt_protection_breaker,
    ) : materialize_emt_breaker_poles(
        backend,
        runtime,
        specification,
        terminal_nodes,
    )
end

function prepare_protection_task_pipeline(pipeline)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :prepare_protection_task_pipeline,
        :emt_protection_breaker,
    ) : prepare_protection_task_pipeline(backend, pipeline)
end

function advance_protection_task_pipeline!(prepared, instant)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :advance_protection_task_pipeline,
        :emt_protection_breaker,
    ) : advance_protection_task_pipeline!(backend, prepared, instant)
end

function prepare_line_fit(request)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:prepare_line_fit, :coupled_line_fitting) :
           prepare_line_fit(backend, request)
end

function snapshot_backend_state(runtime)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:snapshot_backend_state, :backend_snapshot) :
           snapshot_backend_state(backend, runtime)
end

function restore_backend_state!(runtime, snapshot)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(:restore_backend_state, :backend_snapshot) :
           restore_backend_state!(backend, runtime, snapshot)
end

function prepare_partitioned_emt(study)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :prepare_partitioned_emt,
        :local_multirate_partitioned_emt,
    ) : prepare_partitioned_emt(backend, study)
end

function advance_partitioned_emt!(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :advance_partitioned_emt,
        :local_multirate_partitioned_emt,
    ) : advance_partitioned_emt!(backend, prepared)
end

function execute_partitioned_emt!(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :execute_partitioned_emt,
        :local_multirate_partitioned_emt,
    ) : execute_partitioned_emt!(backend, prepared)
end

function partitioned_emt_checkpoint(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :partitioned_emt_checkpoint,
        :local_multirate_partitioned_emt,
    ) : partitioned_emt_checkpoint(backend, prepared)
end

function restore_partitioned_emt_checkpoint!(prepared, checkpoint)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :restore_partitioned_emt_checkpoint,
        :local_multirate_partitioned_emt,
    ) : restore_partitioned_emt_checkpoint!(backend, prepared, checkpoint)
end

function partitioned_emt_status(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :partitioned_emt_status,
        :local_multirate_partitioned_emt,
    ) : partitioned_emt_status(backend, prepared)
end

function prepare_converter_system(study)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :prepare_converter_system,
        :extended_converter_systems,
    ) : prepare_converter_system(backend, study)
end

function advance_converter_system!(prepared, accepted_step_count::Integer=1)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :advance_converter_system,
        :extended_converter_systems,
    ) : advance_converter_system!(backend, prepared, accepted_step_count)
end

function execute_converter_system!(prepared)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :execute_converter_system,
        :extended_converter_systems,
    ) : execute_converter_system!(backend, prepared)
end

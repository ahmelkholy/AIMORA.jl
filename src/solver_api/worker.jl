export STUDY_WORKER_PROTOCOL_VERSION,
       StudyWorkerCapability,
       StudyWorkerLimits,
       StudyWorkerManifest,
       StudyWorkerLifecycle,
       StudyWorkerStarting,
       StudyWorkerReady,
       StudyWorkerRunning,
       StudyWorkerCancelling,
       StudyWorkerStopping,
       StudyWorkerStopped,
       StudyWorkerFailed,
       StudyWorkerCancellation,
       StudyWorkerCancelled,
       StudyWorkerFailure,
       StudyWorkerResultWindow,
       study_worker_manifest,
       study_worker_manifest_record,
       request_study_worker_cancellation!,
       study_worker_cancellation_requested,
       study_worker_cancellation_reason,
       check_study_worker_cancellation,
       study_worker_result_window_record

const STUDY_WORKER_PROTOCOL_VERSION = v"1.0.0"
const _STUDY_WORKER_IDENTIFIER = r"^[a-z][a-z0-9_.-]{0,63}$"
const _STUDY_WORKER_OPAQUE_ID = r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$"
const _STUDY_WORKER_CHECKSUM = r"^[0-9a-f]{64}$"
const _STUDY_WORKER_FAILURE_CODES = (
    :cancelled,
    :invalid_request,
    :protocol_mismatch,
    :resource_limit,
    :worker_crashed,
    :worker_failed,
    :internal_error,
)

function _study_worker_identifier(value, label::AbstractString)
    text = string(value)
    occursin(_STUDY_WORKER_IDENTIFIER, text) || throw(ArgumentError(
        "$label must use a lower-case dotted identifier without path separators",
    ))
    return Symbol(text)
end

function _study_worker_bounded_integer(
    value::Integer,
    minimum::Int,
    maximum::Int,
    label::AbstractString,
)
    converted = Int(value)
    minimum <= converted <= maximum || throw(ArgumentError(
        "$label must be between $minimum and $maximum",
    ))
    return converted
end

function _study_worker_public_text(
    value::AbstractString,
    label::AbstractString;
    maximum_bytes::Int=512,
    allow_empty::Bool=false,
    path_free::Bool=true,
)
    text = String(value)
    (!allow_empty && isempty(text)) && throw(ArgumentError("$label must not be empty"))
    ncodeunits(text) <= maximum_bytes || throw(ArgumentError(
        "$label exceeds the public worker boundary",
    ))
    any(character -> character in ('\0', '\n', '\r'), text) && throw(ArgumentError(
        "$label must not contain control-line characters",
    ))
    if path_free
        contains_private_detail =
            occursin("AIMORASolvers", text) ||
            occursin("ProductionBackend", text) ||
            occursin(
                r"(^|[[:space:](])(?:/|~[/\\]|[A-Za-z]:[\\/])",
                text,
            )
        contains_private_detail && throw(ArgumentError(
            "$label must not contain private type names or filesystem paths",
        ))
    end
    return text
end

_study_worker_version_string(version::VersionNumber) =
    "$(version.major).$(version.minor).$(version.patch)"

"""One public, versioned capability advertised by an isolated study worker."""
struct StudyWorkerCapability
    id::Symbol
    version::VersionNumber

    function StudyWorkerCapability(id, version::VersionNumber)
        version >= v"0.1.0" || throw(ArgumentError(
            "study worker capability versions must be positive",
        ))
        return new(_study_worker_identifier(id, "study worker capability"), version)
    end
end

StudyWorkerCapability(id, version::AbstractString) =
    StudyWorkerCapability(id, VersionNumber(version))

"""Public resource ceilings enforced by the service and private worker runtime."""
struct StudyWorkerLimits
    maximum_jobs_per_process::Int
    maximum_recoveries::Int
    cancellation_poll_interval_steps::Int
    maximum_binary_window_bytes::Int
    maximum_result_windows::Int
    heartbeat_interval_ms::Int
    idle_timeout_ms::Int

    function StudyWorkerLimits(;
        maximum_jobs_per_process::Integer=64,
        maximum_recoveries::Integer=3,
        cancellation_poll_interval_steps::Integer=256,
        maximum_binary_window_bytes::Integer=16 * 1024 * 1024,
        maximum_result_windows::Integer=64,
        heartbeat_interval_ms::Integer=1000,
        idle_timeout_ms::Integer=300_000,
    )
        jobs = _study_worker_bounded_integer(
            maximum_jobs_per_process,
            1,
            10_000,
            "maximum jobs per process",
        )
        recoveries = _study_worker_bounded_integer(
            maximum_recoveries,
            0,
            16,
            "maximum worker recoveries",
        )
        cancellation_steps = _study_worker_bounded_integer(
            cancellation_poll_interval_steps,
            1,
            1_000_000,
            "cancellation poll interval",
        )
        window_bytes = _study_worker_bounded_integer(
            maximum_binary_window_bytes,
            1,
            64 * 1024 * 1024,
            "maximum binary window bytes",
        )
        windows = _study_worker_bounded_integer(
            maximum_result_windows,
            1,
            4096,
            "maximum result windows",
        )
        heartbeat = _study_worker_bounded_integer(
            heartbeat_interval_ms,
            10,
            60_000,
            "worker heartbeat interval",
        )
        idle_timeout = _study_worker_bounded_integer(
            idle_timeout_ms,
            heartbeat,
            24 * 60 * 60 * 1000,
            "worker idle timeout",
        )
        return new(
            jobs,
            recoveries,
            cancellation_steps,
            window_bytes,
            windows,
            heartbeat,
            idle_timeout,
        )
    end
end

"""Public worker identity and limits; it contains no private solver types or paths."""
struct StudyWorkerManifest
    protocol_version::VersionNumber
    worker_version::VersionNumber
    backend_id::Symbol
    capabilities::Tuple{Vararg{StudyWorkerCapability}}
    limits::StudyWorkerLimits

    function StudyWorkerManifest(
        protocol_version::VersionNumber,
        worker_version::VersionNumber,
        backend_id,
        capabilities,
        limits::StudyWorkerLimits,
    )
        protocol_version.major == STUDY_WORKER_PROTOCOL_VERSION.major ||
            throw(ArgumentError("unsupported study worker protocol major version"))
        worker_version >= v"0.1.0" || throw(ArgumentError(
            "study worker version must be positive",
        ))
        public_backend = _study_worker_identifier(backend_id, "study worker backend")
        converted = StudyWorkerCapability[]
        for capability in capabilities
            capability isa StudyWorkerCapability || throw(ArgumentError(
                "worker capabilities must use StudyWorkerCapability",
            ))
            push!(converted, capability)
        end
        length(converted) <= 64 || throw(ArgumentError(
            "study worker capability count exceeds the public boundary",
        ))
        sort!(converted; by=capability -> string(capability.id))
        identifiers = getfield.(converted, :id)
        length(unique(identifiers)) == length(identifiers) || throw(ArgumentError(
            "study worker capabilities must be unique",
        ))
        return new(
            protocol_version,
            worker_version,
            public_backend,
            Tuple(converted),
            limits,
        )
    end
end

@enum StudyWorkerLifecycle::UInt8 begin
    StudyWorkerStarting = 0x01
    StudyWorkerReady = 0x02
    StudyWorkerRunning = 0x03
    StudyWorkerCancelling = 0x04
    StudyWorkerStopping = 0x05
    StudyWorkerStopped = 0x06
    StudyWorkerFailed = 0x07
end

"""Thread-safe, cooperative cancellation token for bounded worker execution."""
mutable struct StudyWorkerCancellation
    requested::Threads.Atomic{Bool}
    reason::Base.RefValue{String}
    guard::ReentrantLock
end

StudyWorkerCancellation() = StudyWorkerCancellation(
    Threads.Atomic{Bool}(false),
    Ref(""),
    ReentrantLock(),
)

"""Public cancellation exception with a sanitized, path-free reason."""
struct StudyWorkerCancelled <: Exception
    reason::String
end

Base.showerror(io::IO, failure::StudyWorkerCancelled) =
    print(io, "study worker request cancelled: ", failure.reason)

"""Typed failure safe to cross the service boundary."""
struct StudyWorkerFailure <: Exception
    code::Symbol
    phase::Symbol
    retryable::Bool
    message::String

    function StudyWorkerFailure(
        code,
        phase,
        retryable::Bool,
        message::AbstractString,
    )
        public_code = code isa Symbol ? code : Symbol(code)
        public_code in _STUDY_WORKER_FAILURE_CODES || throw(ArgumentError(
            "unsupported public study worker failure code",
        ))
        public_phase = _study_worker_identifier(phase, "study worker failure phase")
        public_message = _study_worker_public_text(
            message,
            "study worker failure message";
            maximum_bytes=512,
        )
        return new(public_code, public_phase, retryable, public_message)
    end
end

Base.showerror(io::IO, failure::StudyWorkerFailure) =
    print(io, "study worker ", failure.code, " during ", failure.phase, ": ", failure.message)

"""Path-free descriptor for one bounded binary result window."""
struct StudyWorkerResultWindow
    artifact_id::String
    offset_bytes::Int
    window_bytes::Int
    total_bytes::Int
    content_type::String
    element_type::String
    shape::Tuple{Vararg{Int}}
    checksum_sha256::String

    function StudyWorkerResultWindow(
        artifact_id::AbstractString,
        offset_bytes::Integer,
        window_bytes::Integer,
        total_bytes::Integer,
        content_type::AbstractString,
        element_type::AbstractString,
        shape,
        checksum_sha256::AbstractString;
        maximum_window_bytes::Integer=16 * 1024 * 1024,
    )
        public_artifact_id = String(artifact_id)
        occursin(_STUDY_WORKER_OPAQUE_ID, public_artifact_id) || throw(ArgumentError(
            "artifact identifiers must be opaque and path-free",
        ))
        offset = _study_worker_bounded_integer(
            offset_bytes,
            0,
            typemax(Int),
            "result window offset",
        )
        total = _study_worker_bounded_integer(
            total_bytes,
            0,
            typemax(Int),
            "result artifact size",
        )
        maximum = _study_worker_bounded_integer(
            maximum_window_bytes,
            1,
            64 * 1024 * 1024,
            "result window ceiling",
        )
        window = _study_worker_bounded_integer(
            window_bytes,
            0,
            maximum,
            "result window size",
        )
        offset <= total || throw(ArgumentError("result window offset exceeds artifact size"))
        window <= total - offset || throw(ArgumentError(
            "result window exceeds artifact size",
        ))
        public_content_type = _study_worker_public_text(
            content_type,
            "result content type";
            maximum_bytes=128,
            path_free=false,
        )
        occursin(r"^[A-Za-z0-9.+-]+/[A-Za-z0-9.+-]+$", public_content_type) ||
            throw(ArgumentError("result content type is invalid"))
        public_element_type = string(
            _study_worker_identifier(element_type, "result element type"),
        )
        public_shape = tuple((Int(dimension) for dimension in shape)...)
        length(public_shape) <= 8 || throw(ArgumentError(
            "result window rank exceeds the public boundary",
        ))
        all(dimension -> dimension >= 0, public_shape) || throw(ArgumentError(
            "result window shape must be nonnegative",
        ))
        checksum = String(checksum_sha256)
        occursin(_STUDY_WORKER_CHECKSUM, checksum) || throw(ArgumentError(
            "result window checksum must be lower-case SHA-256",
        ))
        return new(
            public_artifact_id,
            offset,
            window,
            total,
            public_content_type,
            public_element_type,
            public_shape,
            checksum,
        )
    end
end

function study_worker_manifest(::Val{backend_id}) where {backend_id}
    return _solver_unavailable_result(
        :study_worker_manifest,
        :study_worker_lifecycle;
        message="The active backend does not publish an isolated study worker contract.",
    )
end

"""Return the active backend's public worker manifest without exposing its private type."""
function study_worker_manifest()
    backend = active_solver_backend()
    backend === nothing && return _solver_unavailable_result(
        :study_worker_manifest,
        :study_worker_lifecycle;
        message="No AIMORA production solver backend is active.",
    )
    metadata = backend_metadata(backend)
    hasproperty(metadata, :name) || return _solver_unavailable_result(
        :study_worker_manifest,
        :study_worker_lifecycle;
        message="The active backend does not publish a stable public identity.",
    )
    identity = getproperty(metadata, :name)
    identity isa Symbol || return _solver_unavailable_result(
        :study_worker_manifest,
        :study_worker_lifecycle;
        message="The active backend identity is not a public symbol.",
    )
    return study_worker_manifest(Val(identity))
end

"""Convert a manifest to schema-friendly primitives only."""
function study_worker_manifest_record(manifest::StudyWorkerManifest)
    return (
        protocol_version=_study_worker_version_string(manifest.protocol_version),
        worker_version=_study_worker_version_string(manifest.worker_version),
        backend_id=string(manifest.backend_id),
        capabilities=map(manifest.capabilities) do capability
            (
                id=string(capability.id),
                version=_study_worker_version_string(capability.version),
            )
        end,
        limits=(
            maximum_jobs_per_process=manifest.limits.maximum_jobs_per_process,
            maximum_recoveries=manifest.limits.maximum_recoveries,
            cancellation_poll_interval_steps=
                manifest.limits.cancellation_poll_interval_steps,
            maximum_binary_window_bytes=manifest.limits.maximum_binary_window_bytes,
            maximum_result_windows=manifest.limits.maximum_result_windows,
            heartbeat_interval_ms=manifest.limits.heartbeat_interval_ms,
            idle_timeout_ms=manifest.limits.idle_timeout_ms,
        ),
    )
end

"""Request cancellation once; repeated requests are idempotent."""
function request_study_worker_cancellation!(
    token::StudyWorkerCancellation,
    reason::AbstractString="requested by the authenticated client",
)
    public_reason = _study_worker_public_text(
        reason,
        "study worker cancellation reason";
        maximum_bytes=256,
    )
    return lock(token.guard) do
        token.requested[] && return false
        token.reason[] = public_reason
        token.requested[] = true
        return true
    end
end

study_worker_cancellation_requested(token::StudyWorkerCancellation) = token.requested[]

function study_worker_cancellation_reason(token::StudyWorkerCancellation)
    return lock(token.guard) do
        String(token.reason[])
    end
end

"""Throw at a solver-owned safe point when cancellation has been requested."""
function check_study_worker_cancellation(token::StudyWorkerCancellation)
    study_worker_cancellation_requested(token) || return nothing
    throw(StudyWorkerCancelled(study_worker_cancellation_reason(token)))
end

"""Convert a binary result descriptor to schema-friendly primitives only."""
function study_worker_result_window_record(window::StudyWorkerResultWindow)
    return (
        artifact_id=window.artifact_id,
        offset_bytes=window.offset_bytes,
        window_bytes=window.window_bytes,
        total_bytes=window.total_bytes,
        content_type=window.content_type,
        element_type=window.element_type,
        shape=window.shape,
        checksum_sha256=window.checksum_sha256,
    )
end

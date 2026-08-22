const _EMT_CHECKPOINT_MAGIC = collect(codeunits("AIMORA-EMT-CHECKPOINT"))
const _EMT_CHECKPOINT_SCHEMA = UInt32(1)
const _EMT_CHECKPOINT_DIGEST_BYTES = 32
const _EMT_CHECKPOINT_HEADER_BYTES =
    length(_EMT_CHECKPOINT_MAGIC) + sizeof(UInt32) + sizeof(UInt64) +
    _EMT_CHECKPOINT_DIGEST_BYTES

struct EMTBinaryCheckpoint{W}
    schema::UInt32
    julia_major::UInt32
    julia_minor::UInt32
    julia_patch::UInt32
    word_size::Int
    little_endian::Bool
    topology_fingerprint::String
    workspace::W
end

function _write_checkpoint_integer(io::IO, value::T) where {T<:Unsigned}
    for shift in 0:8:(8 * sizeof(T) - 8)
        write(io, UInt8((value >> shift) & T(0xff)))
    end
    return nothing
end

function _read_checkpoint_integer(io::IO, ::Type{T}) where {T<:Unsigned}
    value = zero(T)
    for shift in 0:8:(8 * sizeof(T) - 8)
        eof(io) && throw(ArgumentError("EMT checkpoint header is truncated"))
        value |= T(read(io, UInt8)) << shift
    end
    return value
end

function _emt_checkpoint_topology_fingerprint(workspace::EMTStudyWorkspace)
    context = workspace.runtime.context
    element_owner_names = String[
        String(nameof(typeof(element)))
        for element in context.system.elements
    ]
    fields = String[
        "aimora.emt.topology.v1",
        string(context.dt_s),
        string(context.system.node_count),
        join(string.(context.node_names), "\u001f"),
        join(string.(context.element_names), "\u001f"),
        join(element_owner_names, "\u001f"),
        join(string.(workspace.runtime.plan.switch_names), "\u001f"),
        join(string.(workspace.runtime.plan.source_names), "\u001f"),
    ]
    return bytes2hex(sha256(codeunits(join(fields, "\u001e"))))
end

function _validate_emt_checkpoint_workspace(workspace::EMTStudyWorkspace)
    workspace.execution_mode === :monolithic || throw(ArgumentError(
        workspace.execution_mode === :hybrid ?
        "local trusted EMT checkpoints cannot capture hybrid execution without its coordinating integrator; use a portable hybrid snapshot" :
        "local trusted EMT checkpoints require a completed monolithic execution",
    ))
    workspace.ready && throw(ArgumentError(
        "EMT checkpoint requires an evaluated workspace",
    ))
    context = workspace.runtime.context
    context.step_index == context.step_count + 1 || throw(ArgumentError(
        "EMT checkpoint requires a completed timestep horizon",
    ))
    context.trace_write_index == length(context.recorded_step_indices) + 1 ||
        throw(ArgumentError("EMT checkpoint trace storage is incomplete"))
    isempty(context.recorded_step_indices) && throw(ArgumentError(
        "EMT checkpoint requires at least one accepted trace sample",
    ))
    last(context.recorded_step_indices) == context.step_count || throw(ArgumentError(
        "EMT checkpoint does not contain its final accepted timestep",
    ))
    all(isfinite, context.system.v) ||
        throw(ArgumentError("EMT checkpoint node state is nonfinite"))
    _check_prepared_runtime_aliases(workspace.runtime)
    return workspace
end

"""
    write_emt_checkpoint(path, workspace)

Write a completed typed Julia EMT workspace to an integrity-checked binary
checkpoint. The format owns Julia state and deliberately does not reproduce a
Fortran COMMON-memory dump. Checkpoint files are trusted simulation artifacts,
not a safe interchange format for untrusted input.
"""
function write_emt_checkpoint(
    path::AbstractString,
    workspace::EMTStudyWorkspace,
)
    _validate_emt_checkpoint_workspace(workspace)
    payload = EMTBinaryCheckpoint(
        _EMT_CHECKPOINT_SCHEMA,
        VERSION.major,
        VERSION.minor,
        VERSION.patch,
        Sys.WORD_SIZE,
        ENDIAN_BOM == 0x04030201,
        _emt_checkpoint_topology_fingerprint(workspace),
        workspace,
    )
    payload_io = IOBuffer()
    serialize(payload_io, payload)
    payload_bytes = take!(payload_io)
    digest = sha256(payload_bytes)
    output_path = abspath(String(path))
    mkpath(dirname(output_path))
    temporary_path, io = mktemp(dirname(output_path))
    try
        write(io, _EMT_CHECKPOINT_MAGIC)
        _write_checkpoint_integer(io, _EMT_CHECKPOINT_SCHEMA)
        _write_checkpoint_integer(io, UInt64(length(payload_bytes)))
        write(io, digest)
        write(io, payload_bytes)
        close(io)
        mv(temporary_path, output_path; force = true)
    catch
        isopen(io) && close(io)
        isfile(temporary_path) && rm(temporary_path; force = true)
        rethrow()
    end
    return output_path
end

"""
    read_emt_checkpoint(path; max_payload_bytes=8_000_000_000)

Read and validate an AIMORA EMT checkpoint in a fresh Julia process.
"""
function read_emt_checkpoint(
    path::AbstractString;
    max_payload_bytes::Integer = 8_000_000_000,
)
    max_payload_bytes > 0 ||
        throw(ArgumentError("checkpoint payload limit must be positive"))
    input_path = abspath(String(path))
    isfile(input_path) ||
        throw(ArgumentError("EMT checkpoint file does not exist: $input_path"))
    filesize(input_path) >= _EMT_CHECKPOINT_HEADER_BYTES ||
        throw(ArgumentError("EMT checkpoint file is shorter than its header"))
    payload_bytes = open(input_path, "r") do io
        magic = read(io, length(_EMT_CHECKPOINT_MAGIC))
        magic == _EMT_CHECKPOINT_MAGIC ||
            throw(ArgumentError("file is not an AIMORA EMT checkpoint"))
        schema = _read_checkpoint_integer(io, UInt32)
        schema == _EMT_CHECKPOINT_SCHEMA || throw(ArgumentError(
            "unsupported EMT checkpoint schema $schema",
        ))
        payload_length = _read_checkpoint_integer(io, UInt64)
        payload_length <= UInt64(max_payload_bytes) || throw(ArgumentError(
            "EMT checkpoint payload exceeds the configured byte limit",
        ))
        payload_length <= UInt64(typemax(Int)) || throw(ArgumentError(
            "EMT checkpoint payload cannot be addressed by this Julia process",
        ))
        expected_size = UInt64(_EMT_CHECKPOINT_HEADER_BYTES) + payload_length
        UInt64(filesize(input_path)) == expected_size || throw(ArgumentError(
            "EMT checkpoint payload length does not match the file size",
        ))
        expected_digest = read(io, _EMT_CHECKPOINT_DIGEST_BYTES)
        length(expected_digest) == _EMT_CHECKPOINT_DIGEST_BYTES ||
            throw(ArgumentError("EMT checkpoint digest is truncated"))
        bytes = read(io, Int(payload_length))
        sha256(bytes) == expected_digest ||
            throw(ArgumentError("EMT checkpoint integrity digest does not match"))
        return bytes
    end
    checkpoint = try
        deserialize(IOBuffer(payload_bytes))
    catch error
        throw(ArgumentError(
            "EMT checkpoint payload could not be decoded: " *
            sprint(showerror, error),
        ))
    end
    checkpoint isa EMTBinaryCheckpoint || throw(ArgumentError(
        "EMT checkpoint payload has an unexpected root type",
    ))
    checkpoint.schema == _EMT_CHECKPOINT_SCHEMA ||
        throw(ArgumentError("EMT checkpoint payload schema is inconsistent"))
    (checkpoint.julia_major, checkpoint.julia_minor, checkpoint.julia_patch) ==
    (VERSION.major, VERSION.minor, VERSION.patch) || throw(ArgumentError(
        "EMT checkpoint was written by Julia " *
        "$(checkpoint.julia_major).$(checkpoint.julia_minor).$(checkpoint.julia_patch), " *
        "but this process uses $VERSION",
    ))
    checkpoint.word_size == Sys.WORD_SIZE || throw(ArgumentError(
        "EMT checkpoint word size $(checkpoint.word_size) does not match " *
        "this $(Sys.WORD_SIZE)-bit process",
    ))
    checkpoint.little_endian == (ENDIAN_BOM == 0x04030201) ||
        throw(ArgumentError("EMT checkpoint byte order does not match this process"))
    workspace = checkpoint.workspace
    workspace isa EMTStudyWorkspace || throw(ArgumentError(
        "EMT checkpoint payload does not contain an EMT workspace",
    ))
    _validate_emt_checkpoint_workspace(workspace)
    fingerprint = _emt_checkpoint_topology_fingerprint(workspace)
    fingerprint == checkpoint.topology_fingerprint || throw(ArgumentError(
        "EMT checkpoint topology fingerprint does not match its workspace",
    ))
    return workspace
end

"""
    restart_emt_checkpoint(path, request; additional_time_s, ...)

Import a completed checkpoint, apply restart cards, resume the solver, and
optionally write the restart report and the continued checkpoint.
"""
function restart_emt_checkpoint(
    path::AbstractString,
    request::DeckParser.DeckRestartRequest;
    additional_time_s::Real,
    recorded_step_indices = nothing,
    report_path::Union{Nothing,AbstractString} = nothing,
    continued_checkpoint_path::Union{Nothing,AbstractString} = nothing,
)
    workspace = read_emt_checkpoint(path)
    run = restart_emt_study!(
        workspace,
        request;
        additional_time_s = additional_time_s,
        recorded_step_indices = recorded_step_indices,
    )
    written_report = report_path === nothing ?
        nothing : write_emt_restart_report(report_path, run)
    written_checkpoint = continued_checkpoint_path === nothing ?
        nothing : write_emt_checkpoint(continued_checkpoint_path, workspace)
    return (
        workspace = workspace,
        run = run,
        report_path = written_report,
        continued_checkpoint_path = written_checkpoint,
    )
end

function restart_emt_checkpoint(
    request::DeckParser.DeckRestartRequest;
    kwargs...,
)
    request.checkpoint_path === nothing && throw(ArgumentError(
        "START AGAIN request does not name an AIMORA checkpoint",
    ))
    return restart_emt_checkpoint(request.checkpoint_path, request; kwargs...)
end

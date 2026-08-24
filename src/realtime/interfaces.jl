export AbstractRealtimeInterface,
       InProcessControllerInterface,
       LocalUDPControllerInterface,
       SharedLibraryControllerInterface,
       open_local_udp_controller,
       open_shared_library_controller,
       shared_library_controller_state,
       restore_shared_library_controller_state!,
       prepare_realtime_interface!,
       exchange_realtime_interface!,
       realtime_interface_checkpoint,
       capture_realtime_interface_state!,
       restore_realtime_interface_state!,
       safe_shutdown_realtime_interface!,
       close_realtime_interface!

abstract type AbstractRealtimeInterface end

realtime_interface_kind(::AbstractRealtimeInterface) = :unknown

mutable struct InProcessControllerInterface{F,S} <: AbstractRealtimeInterface
    controller!::F
    controller_state::S
    input_values::Vector{Float64}
    output_values::Vector{Float64}
    safe_output_values::Vector{Float64}
    sequence::UInt64
    prepared::Bool
    closed::Bool
end

realtime_interface_kind(::InProcessControllerInterface) = :in_process

mutable struct LocalUDPControllerCheckpoint{S}
    controller_state::S
    sequence::UInt64
end

mutable struct LocalUDPControllerInterface{F,S} <: AbstractRealtimeInterface
    controller!::F
    controller_state::S
    model_socket::Cint
    controller_socket::Cint
    request_frame::Vector{UInt8}
    response_frame::Vector{UInt8}
    receive_frame::Vector{UInt8}
    semantic_sha256::NTuple{32,UInt8}
    input_values::Vector{Float64}
    output_values::Vector{Float64}
    safe_output_values::Vector{Float64}
    sequence::UInt64
    prepared::Bool
    closed::Bool
end

realtime_interface_kind(::LocalUDPControllerInterface) = :local_udp

const _REALTIME_UDP_MAGIC = UInt32(0x41494d52)
const _REALTIME_UDP_ABI = UInt32(1)
const _REALTIME_UDP_HEADER_BYTES = 64

function _local_udp_address(port::UInt16)
    address = zeros(UInt8, 16)
    address[1] = 0x02
    address[2] = 0x00
    address[3] = UInt8(port >> 8)
    address[4] = UInt8(port & 0xff)
    address[5:8] .= (0x7f, 0x00, 0x00, 0x01)
    return address
end

function _bound_local_udp_socket()
    descriptor = ccall(:socket, Cint, (Cint, Cint, Cint), 2, 2, 0)
    descriptor >= 0 || error("failed to create local UDP socket")
    address = _local_udp_address(0x0000)
    status = GC.@preserve address ccall(
        :bind,
        Cint,
        (Cint, Ptr{UInt8}, Cuint),
        descriptor,
        pointer(address),
        Cuint(length(address)),
    )
    if status != 0
        ccall(:close, Cint, (Cint,), descriptor)
        error("failed to bind local UDP socket")
    end
    length_ref = Ref{Cuint}(Cuint(length(address)))
    status = GC.@preserve address ccall(
        :getsockname,
        Cint,
        (Cint, Ptr{UInt8}, Ref{Cuint}),
        descriptor,
        pointer(address),
        length_ref,
    )
    if status != 0 || length_ref[] < 8
        ccall(:close, Cint, (Cint,), descriptor)
        error("failed to identify local UDP socket")
    end
    port = (UInt16(address[3]) << 8) | UInt16(address[4])
    return descriptor, port
end

function _connect_local_udp_socket(descriptor::Cint, port::UInt16)
    address = _local_udp_address(port)
    status = GC.@preserve address ccall(
        :connect,
        Cint,
        (Cint, Ptr{UInt8}, Cuint),
        descriptor,
        pointer(address),
        Cuint(length(address)),
    )
    status == 0 || error("failed to connect local UDP socket")
    return descriptor
end

function _local_udp_semantic_sha256(channels::AbstractVector{RealtimeChannel})
    context = SHA.SHA256_CTX()
    SHA.update!(context, codeunits("aimora-local-udp-controller-v1\n"))
    for channel in channels
        SHA.update!(context, codeunits(join((
            channel.id,
            channel.direction,
            channel.physical_unit,
            repr(channel.scale),
            repr(channel.offset),
            repr(channel.raw_minimum),
            repr(channel.raw_maximum),
            repr(channel.safe_physical_value),
        ), '\t')))
        SHA.update!(context, UInt8['\n'])
    end
    return Tuple(SHA.digest!(context))
end

function _local_udp_frame_size(value_count::Integer)
    return Base.Checked.checked_add(
        _REALTIME_UDP_HEADER_BYTES,
        Base.Checked.checked_mul(Int(value_count), sizeof(Float64)),
    )
end

function open_local_udp_controller(
    controller!::F,
    controller_state::S,
    channels::AbstractVector{RealtimeChannel},
) where {F,S}
    Sys.islinux() || throw(ArgumentError(
        "the admitted local UDP controller is native-Linux only",
    ))
    isbitstype(S) || throw(ArgumentError(
        "local UDP controller state must be an immutable bits value for checkpointing",
    ))
    inputs = filter(channel -> channel.direction == :input, channels)
    outputs = filter(channel -> channel.direction == :output, channels)
    model_socket, model_port = _bound_local_udp_socket()
    controller_socket = Cint(-1)
    try
        controller_socket, controller_port = _bound_local_udp_socket()
        _connect_local_udp_socket(model_socket, controller_port)
        _connect_local_udp_socket(controller_socket, model_port)
        request_size = _local_udp_frame_size(length(outputs))
        response_size = _local_udp_frame_size(length(inputs))
        return LocalUDPControllerInterface{F,S}(
            controller!,
            controller_state,
            model_socket,
            controller_socket,
            zeros(UInt8, request_size),
            zeros(UInt8, response_size),
            zeros(UInt8, max(request_size, response_size)),
            _local_udp_semantic_sha256(channels),
            zeros(length(inputs)),
            zeros(length(outputs)),
            getfield.(outputs, :safe_physical_value),
            0,
            false,
            false,
        )
    catch
        ccall(:close, Cint, (Cint,), model_socket)
        controller_socket >= 0 && ccall(:close, Cint, (Cint,), controller_socket)
        rethrow()
    end
end

realtime_interface_checkpoint(interface::LocalUDPControllerInterface) =
    LocalUDPControllerCheckpoint(interface.controller_state, interface.sequence)

function capture_realtime_interface_state!(
    checkpoint::LocalUDPControllerCheckpoint,
    interface::LocalUDPControllerInterface,
)
    checkpoint.controller_state = interface.controller_state
    checkpoint.sequence = interface.sequence
    return checkpoint
end

function restore_realtime_interface_state!(
    interface::LocalUDPControllerInterface,
    checkpoint::LocalUDPControllerCheckpoint,
)
    interface.controller_state = checkpoint.controller_state
    interface.sequence = checkpoint.sequence
    return interface
end

function realtime_interface_checkpoint(interface::InProcessControllerInterface)
    interface.controller_state === nothing || isbitstype(typeof(interface.controller_state)) ||
        throw(ArgumentError(
            "mutable in-process controller state requires a separately registered snapshot owner",
        ))
    return nothing
end

capture_realtime_interface_state!(::Nothing, ::InProcessControllerInterface) = nothing
restore_realtime_interface_state!(::InProcessControllerInterface, ::Nothing) = nothing

mutable struct SharedLibraryControllerInterface <: AbstractRealtimeInterface
    library::Ptr{Cvoid}
    step_function::Ptr{Cvoid}
    reset_function::Ptr{Cvoid}
    content_sha256::String
    state_words::Vector{UInt64}
    input_values::Vector{Float64}
    output_values::Vector{Float64}
    safe_output_values::Vector{Float64}
    sequence::UInt64
    prepared::Bool
    closed::Bool
end

realtime_interface_kind(::SharedLibraryControllerInterface) = :shared_library_controller

realtime_interface_checkpoint(interface::SharedLibraryControllerInterface) =
    similar(interface.state_words)

function capture_realtime_interface_state!(
    checkpoint::Vector{UInt64},
    interface::SharedLibraryControllerInterface,
)
    copyto!(checkpoint, interface.state_words)
    return checkpoint
end

function restore_realtime_interface_state!(
    interface::SharedLibraryControllerInterface,
    checkpoint::Vector{UInt64},
)
    copyto!(interface.state_words, checkpoint)
    return interface
end

function _shared_library_symbol(library::Ptr{Cvoid}, name::Symbol)
    symbol = Libdl.dlsym_e(library, name)
    symbol == C_NULL && throw(ArgumentError(
        "shared-library real-time controller omits required symbol $name",
    ))
    return symbol
end

function open_shared_library_controller(
    path::AbstractString,
    expected_sha256::AbstractString,
    channels::AbstractVector{RealtimeChannel},
)
    isfile(path) || throw(ArgumentError(
        "shared-library real-time controller file is unavailable",
    ))
    expected = lowercase(String(expected_sha256))
    occursin(r"^[0-9a-f]{64}$", expected) || throw(ArgumentError(
        "shared-library real-time controller SHA-256 must be 64 lowercase hex characters",
    ))
    actual = bytes2hex(sha256(read(path)))
    actual == expected || throw(ArgumentError(
        "shared-library real-time controller content hash mismatch",
    ))
    inputs = filter(channel -> channel.direction == :input, channels)
    outputs = filter(channel -> channel.direction == :output, channels)
    library = Libdl.dlopen(path)
    try
        abi_function = _shared_library_symbol(
            library,
            :aimora_realtime_controller_abi_version,
        )
        state_size_function = _shared_library_symbol(
            library,
            :aimora_realtime_controller_state_bytes,
        )
        reset_function = _shared_library_symbol(
            library,
            :aimora_realtime_controller_reset,
        )
        step_function = _shared_library_symbol(
            library,
            :aimora_realtime_controller_step,
        )
        abi_version = ccall(abi_function, Cuint, ())
        abi_version == 1 || throw(ArgumentError(
            "shared-library real-time controller ABI version is unsupported",
        ))
        state_bytes = Int(ccall(state_size_function, Csize_t, ()))
        0 < state_bytes <= 1024 * 1024 || throw(ArgumentError(
            "shared-library real-time controller state size is outside the admitted range",
        ))
        state_words = zeros(UInt64, cld(state_bytes, sizeof(UInt64)))
        status = ccall(
            reset_function,
            Cint,
            (Ptr{Cvoid}, Csize_t),
            pointer(state_words),
            state_bytes,
        )
        status == 0 || throw(ArgumentError(
            "shared-library real-time controller reset failed with status $status",
        ))
        return SharedLibraryControllerInterface(
            library,
            step_function,
            reset_function,
            actual,
            state_words,
            zeros(length(inputs)),
            zeros(length(outputs)),
            getfield.(outputs, :safe_physical_value),
            0,
            false,
            false,
        )
    catch
        Libdl.dlclose(library)
        rethrow()
    end
end

function prepare_realtime_target(
    target::RealtimeTarget,
    channels::AbstractVector{RealtimeChannel},
    interface::AbstractRealtimeInterface,
)
    preparation = prepare_realtime_target(target, channels)
    if target.kind == LocalUDPControllerTarget &&
            preparation isa RealtimeUnavailableResult &&
            preparation.code == :local_udp_interface_not_bound
        interface isa LocalUDPControllerInterface || return preparation
        length(interface.input_values) == count(
            channel -> channel.direction == :input,
            channels,
        ) || throw(DimensionMismatch(
            "local UDP controller input channel count mismatch",
        ))
        length(interface.output_values) == count(
            channel -> channel.direction == :output,
            channels,
        ) || throw(DimensionMismatch(
            "local UDP controller output channel count mismatch",
        ))
        interface.semantic_sha256 == _local_udp_semantic_sha256(channels) ||
            throw(ArgumentError("local UDP controller channel identity mismatch"))
        maximum(length.((interface.request_frame, interface.response_frame))) <=
            target.maximum_payload_bytes || throw(ArgumentError(
                "local UDP controller frame exceeds the target payload limit",
            ))
        return RealtimePreparation(
            target,
            collect(channels),
            :linux_local_udp_soft,
            false,
            false,
        )
    end
    if target.kind == SharedLibraryControllerTarget &&
            preparation isa RealtimeUnavailableResult &&
            preparation.code == :controller_library_not_bound
        interface isa SharedLibraryControllerInterface || return preparation
        length(interface.input_values) == count(
            channel -> channel.direction == :input,
            channels,
        ) || throw(DimensionMismatch(
            "shared-library controller input channel count mismatch",
        ))
        length(interface.output_values) == count(
            channel -> channel.direction == :output,
            channels,
        ) || throw(DimensionMismatch(
            "shared-library controller output channel count mismatch",
        ))
        return RealtimePreparation(
            target,
            collect(channels),
            :linux_shared_library_soft,
            false,
            false,
        )
    end
    return preparation
end

function shared_library_controller_state(
    interface::SharedLibraryControllerInterface,
)
    interface.closed && throw(ArgumentError(
        "closed shared-library controller has no restorable state",
    ))
    return copy(interface.state_words)
end

function restore_shared_library_controller_state!(
    interface::SharedLibraryControllerInterface,
    state_words::AbstractVector{UInt64},
)
    interface.closed && throw(ArgumentError(
        "closed shared-library controller cannot restore state",
    ))
    length(state_words) == length(interface.state_words) || throw(DimensionMismatch(
        "shared-library controller state size mismatch",
    ))
    copyto!(interface.state_words, state_words)
    return interface
end

function InProcessControllerInterface(
    controller!::F,
    controller_state::S,
    channels::AbstractVector{RealtimeChannel},
) where {F,S}
    inputs = filter(channel -> channel.direction == :input, channels)
    outputs = filter(channel -> channel.direction == :output, channels)
    return InProcessControllerInterface{F,S}(
        controller!,
        controller_state,
        zeros(length(inputs)),
        zeros(length(outputs)),
        getfield.(outputs, :safe_physical_value),
        0,
        false,
        false,
    )
end

function prepare_realtime_interface!(interface::InProcessControllerInterface)
    interface.closed && throw(ArgumentError(
        "closed in-process real-time interface cannot be prepared",
    ))
    interface.sequence = 0
    fill!(interface.input_values, 0.0)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = true
    return interface
end

function prepare_realtime_interface!(interface::LocalUDPControllerInterface)
    interface.closed && throw(ArgumentError(
        "closed local UDP real-time interface cannot be prepared",
    ))
    interface.sequence = 0
    fill!(interface.input_values, 0.0)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = true
    return interface
end

function prepare_realtime_interface!(
    interface::SharedLibraryControllerInterface,
)
    interface.closed && throw(ArgumentError(
        "closed shared-library real-time interface cannot be prepared",
    ))
    interface.sequence = 0
    fill!(interface.input_values, 0.0)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = true
    return interface
end

function exchange_realtime_interface!(
    interface::InProcessControllerInterface,
    destination_inputs::AbstractVector{Float64},
    source_outputs::AbstractVector{Float64},
    logical_step::Int,
    logical_time_ns::Int64,
)
    interface.prepared || throw(ArgumentError(
        "in-process real-time interface is not prepared",
    ))
    interface.closed && throw(ArgumentError(
        "in-process real-time interface is closed",
    ))
    length(destination_inputs) == length(interface.input_values) ||
        throw(DimensionMismatch("real-time input frame size mismatch"))
    length(source_outputs) == length(interface.output_values) ||
        throw(DimensionMismatch("real-time output frame size mismatch"))
    all(isfinite, source_outputs) || throw(ArgumentError(
        "real-time model outputs must be finite",
    ))
    copyto!(interface.output_values, source_outputs)
    interface.controller!(
        interface.input_values,
        interface.output_values,
        interface.controller_state,
        logical_step,
        logical_time_ns,
    )
    all(isfinite, interface.input_values) || throw(ArgumentError(
        "real-time controller inputs must be finite",
    ))
    copyto!(destination_inputs, interface.input_values)
    interface.sequence = Base.Checked.checked_add(interface.sequence, UInt64(1))
    return interface
end

function _store_local_udp_integer!(
    frame::Vector{UInt8},
    byte_offset::Int,
    value::T,
) where {T<:Union{UInt32,UInt64,Int64}}
    GC.@preserve frame unsafe_store!(
        Ptr{T}(pointer(frame, byte_offset + 1)),
        htol(value),
    )
    return frame
end

function _load_local_udp_integer(
    frame::Vector{UInt8},
    byte_offset::Int,
    ::Type{T},
) where {T<:Union{UInt32,UInt64,Int64}}
    value = GC.@preserve frame unsafe_load(Ptr{T}(pointer(frame, byte_offset + 1)))
    return ltoh(value)
end

function _encode_local_udp_frame!(
    frame::Vector{UInt8},
    semantic_sha256::NTuple{32,UInt8},
    sequence::UInt64,
    logical_time_ns::Int64,
    values::Vector{Float64},
)
    length(frame) == _local_udp_frame_size(length(values)) || throw(DimensionMismatch(
        "local UDP transmit frame size mismatch",
    ))
    _store_local_udp_integer!(frame, 0, _REALTIME_UDP_MAGIC)
    _store_local_udp_integer!(frame, 4, _REALTIME_UDP_ABI)
    _store_local_udp_integer!(frame, 8, sequence)
    _store_local_udp_integer!(frame, 16, logical_time_ns)
    _store_local_udp_integer!(frame, 24, UInt64(length(values)))
    for index in eachindex(semantic_sha256)
        frame[32 + index] = semantic_sha256[index]
    end
    for index in eachindex(values)
        isfinite(values[index]) || throw(ArgumentError(
            "local UDP frame contains a nonfinite value",
        ))
        _store_local_udp_integer!(
            frame,
            _REALTIME_UDP_HEADER_BYTES + (index - 1) * sizeof(Float64),
            reinterpret(UInt64, values[index]),
        )
    end
    return frame
end

function _decode_local_udp_frame!(
    destination::Vector{Float64},
    frame::Vector{UInt8},
    semantic_sha256::NTuple{32,UInt8},
    sequence::UInt64,
    logical_time_ns::Int64,
)
    _load_local_udp_integer(frame, 0, UInt32) == _REALTIME_UDP_MAGIC ||
        throw(ArgumentError("local UDP frame magic mismatch"))
    _load_local_udp_integer(frame, 4, UInt32) == _REALTIME_UDP_ABI ||
        throw(ArgumentError("local UDP frame ABI mismatch"))
    _load_local_udp_integer(frame, 8, UInt64) == sequence ||
        throw(ArgumentError("local UDP frame sequence mismatch"))
    _load_local_udp_integer(frame, 16, Int64) == logical_time_ns ||
        throw(ArgumentError("local UDP frame timestamp mismatch"))
    _load_local_udp_integer(frame, 24, UInt64) == UInt64(length(destination)) ||
        throw(DimensionMismatch("local UDP frame channel count mismatch"))
    for index in eachindex(semantic_sha256)
        frame[32 + index] == semantic_sha256[index] || throw(ArgumentError(
            "local UDP frame configuration identity mismatch",
        ))
    end
    for index in eachindex(destination)
        bits = _load_local_udp_integer(
            frame,
            _REALTIME_UDP_HEADER_BYTES + (index - 1) * sizeof(Float64),
            UInt64,
        )
        destination[index] = reinterpret(Float64, bits)
        isfinite(destination[index]) || throw(ArgumentError(
            "local UDP frame contains a nonfinite value",
        ))
    end
    return destination
end

function _send_local_udp_frame(descriptor::Cint, frame::Vector{UInt8})
    transferred = GC.@preserve frame ccall(
        :send,
        Cssize_t,
        (Cint, Ptr{UInt8}, Csize_t, Cint),
        descriptor,
        pointer(frame),
        length(frame),
        0,
    )
    transferred == length(frame) || error("local UDP frame send failed")
    return nothing
end

function _receive_local_udp_frame!(descriptor::Cint, frame::Vector{UInt8}, size::Int)
    transferred = GC.@preserve frame ccall(
        :recv,
        Cssize_t,
        (Cint, Ptr{UInt8}, Csize_t, Cint),
        descriptor,
        pointer(frame),
        size,
        0,
    )
    transferred == size || error("local UDP frame receive failed")
    return frame
end

function exchange_realtime_interface!(
    interface::LocalUDPControllerInterface,
    destination_inputs::AbstractVector{Float64},
    source_outputs::AbstractVector{Float64},
    logical_step::Int,
    logical_time_ns::Int64,
)
    interface.prepared || throw(ArgumentError(
        "local UDP real-time interface is not prepared",
    ))
    interface.closed && throw(ArgumentError(
        "local UDP real-time interface is closed",
    ))
    length(destination_inputs) == length(interface.input_values) ||
        throw(DimensionMismatch("local UDP input frame size mismatch"))
    length(source_outputs) == length(interface.output_values) ||
        throw(DimensionMismatch("local UDP output frame size mismatch"))
    copyto!(interface.output_values, source_outputs)
    sequence = Base.Checked.checked_add(interface.sequence, UInt64(1))
    _encode_local_udp_frame!(
        interface.request_frame,
        interface.semantic_sha256,
        sequence,
        logical_time_ns,
        interface.output_values,
    )
    _send_local_udp_frame(interface.model_socket, interface.request_frame)
    _receive_local_udp_frame!(
        interface.controller_socket,
        interface.receive_frame,
        length(interface.request_frame),
    )
    _decode_local_udp_frame!(
        interface.output_values,
        interface.receive_frame,
        interface.semantic_sha256,
        sequence,
        logical_time_ns,
    )
    updated_state = interface.controller!(
        interface.input_values,
        interface.output_values,
        interface.controller_state,
        logical_step,
        logical_time_ns,
    )
    updated_state === nothing || (interface.controller_state = updated_state)
    _encode_local_udp_frame!(
        interface.response_frame,
        interface.semantic_sha256,
        sequence,
        logical_time_ns,
        interface.input_values,
    )
    _send_local_udp_frame(interface.controller_socket, interface.response_frame)
    _receive_local_udp_frame!(
        interface.model_socket,
        interface.receive_frame,
        length(interface.response_frame),
    )
    _decode_local_udp_frame!(
        interface.input_values,
        interface.receive_frame,
        interface.semantic_sha256,
        sequence,
        logical_time_ns,
    )
    copyto!(destination_inputs, interface.input_values)
    interface.sequence = sequence
    return interface
end

function exchange_realtime_interface!(
    interface::SharedLibraryControllerInterface,
    destination_inputs::AbstractVector{Float64},
    source_outputs::AbstractVector{Float64},
    logical_step::Int,
    logical_time_ns::Int64,
)
    interface.prepared || throw(ArgumentError(
        "shared-library real-time interface is not prepared",
    ))
    interface.closed && throw(ArgumentError(
        "shared-library real-time interface is closed",
    ))
    length(destination_inputs) == length(interface.input_values) ||
        throw(DimensionMismatch("real-time input frame size mismatch"))
    length(source_outputs) == length(interface.output_values) ||
        throw(DimensionMismatch("real-time output frame size mismatch"))
    all(isfinite, source_outputs) || throw(ArgumentError(
        "real-time model outputs must be finite",
    ))
    copyto!(interface.output_values, source_outputs)
    status = ccall(
        interface.step_function,
        Cint,
        (
            Ptr{Cvoid},
            Ptr{Cdouble},
            Csize_t,
            Ptr{Cdouble},
            Csize_t,
            UInt64,
            Int64,
        ),
        pointer(interface.state_words),
        pointer(interface.output_values),
        length(interface.output_values),
        pointer(interface.input_values),
        length(interface.input_values),
        UInt64(logical_step),
        logical_time_ns,
    )
    status == 0 || throw(ArgumentError(
        "shared-library real-time controller step failed with status $status",
    ))
    all(isfinite, interface.input_values) || throw(ArgumentError(
        "shared-library real-time controller returned nonfinite input",
    ))
    copyto!(destination_inputs, interface.input_values)
    interface.sequence = Base.Checked.checked_add(interface.sequence, UInt64(1))
    return interface
end

function safe_shutdown_realtime_interface!(
    interface::InProcessControllerInterface,
)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = false
    return interface
end

function safe_shutdown_realtime_interface!(interface::LocalUDPControllerInterface)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = false
    return interface
end

function safe_shutdown_realtime_interface!(
    interface::SharedLibraryControllerInterface,
)
    copyto!(interface.output_values, interface.safe_output_values)
    interface.prepared = false
    return interface
end

function close_realtime_interface!(interface::InProcessControllerInterface)
    safe_shutdown_realtime_interface!(interface)
    interface.closed = true
    return interface
end

function close_realtime_interface!(interface::LocalUDPControllerInterface)
    interface.closed && return interface
    safe_shutdown_realtime_interface!(interface)
    ccall(:close, Cint, (Cint,), interface.model_socket)
    ccall(:close, Cint, (Cint,), interface.controller_socket)
    interface.model_socket = Cint(-1)
    interface.controller_socket = Cint(-1)
    interface.closed = true
    return interface
end

function close_realtime_interface!(
    interface::SharedLibraryControllerInterface,
)
    interface.closed && return interface
    safe_shutdown_realtime_interface!(interface)
    Libdl.dlclose(interface.library)
    interface.library = C_NULL
    interface.step_function = C_NULL
    interface.reset_function = C_NULL
    interface.closed = true
    return interface
end

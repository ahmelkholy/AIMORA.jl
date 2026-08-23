struct _EMTInitializationDigestWriter
    context::SHA.SHA2_256_CTX
    pending_bytes::Vector{UInt8}
    symbol_identifiers::Dict{Symbol,UInt64}
    type_identifiers::Dict{Tuple{DataType,Bool},UInt64}
end

const _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES = 16 * 1024

_EMTInitializationDigestWriter() =
    _EMTInitializationDigestWriter(
        SHA.SHA2_256_CTX(),
        UInt8[],
        Dict{Symbol,UInt64}(),
        Dict{Tuple{DataType,Bool},UInt64}(),
    )

function _flush_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
)
    isempty(writer.pending_bytes) && return writer
    SHA.update!(writer.context, writer.pending_bytes)
    empty!(writer.pending_bytes)
    return writer
end

function _update_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
    bytes::NTuple{N,UInt8},
) where {N}
    pending_bytes = writer.pending_bytes
    for byte in bytes
        push!(pending_bytes, byte)
    end
    length(pending_bytes) >=
        _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES &&
        _flush_emt_initialization_digest!(writer)
    return writer
end

function _update_emt_initialization_digest!(
    writer::_EMTInitializationDigestWriter,
    bytes::AbstractVector{UInt8},
)
    if length(bytes) >= _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES
        _flush_emt_initialization_digest!(writer)
        SHA.update!(writer.context, bytes)
    else
        append!(writer.pending_bytes, bytes)
        length(writer.pending_bytes) >=
            _EMT_INITIALIZATION_DIGEST_BUFFER_BYTES &&
            _flush_emt_initialization_digest!(writer)
    end
    return writer
end

function _write_emt_initialization_digest_tag!(
    writer::_EMTInitializationDigestWriter,
    tag::UInt8,
)
    return _update_emt_initialization_digest!(writer, (tag,))
end

function _write_emt_initialization_digest_uint64!(
    writer::_EMTInitializationDigestWriter,
    value::UInt64,
)
    bytes = ntuple(
        index -> UInt8((value >> (8 * (index - 1))) & 0xff),
        Val(8),
    )
    return _update_emt_initialization_digest!(writer, bytes)
end

function _write_emt_initialization_digest_text!(
    writer::_EMTInitializationDigestWriter,
    value::AbstractString,
)
    bytes = codeunits(value)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(bytes)))
    return _update_emt_initialization_digest!(writer, bytes)
end

function _write_emt_initialization_digest_symbol!(
    writer::_EMTInitializationDigestWriter,
    value::Symbol,
)
    identifier = get(writer.symbol_identifiers, value, UInt64(0))
    if identifier != 0
        _write_emt_initialization_digest_tag!(writer, 0x00)
        return _write_emt_initialization_digest_uint64!(writer, identifier)
    end
    identifier = UInt64(length(writer.symbol_identifiers) + 1)
    writer.symbol_identifiers[value] = identifier
    _write_emt_initialization_digest_tag!(writer, 0x01)
    _write_emt_initialization_digest_uint64!(writer, identifier)
    return _write_emt_initialization_digest_text!(writer, String(value))
end

function _write_emt_initialization_digest_type!(
    writer::_EMTInitializationDigestWriter,
    value::DataType;
    include_parameters::Bool=true,
)
    key = (value, include_parameters)
    identifier = get(writer.type_identifiers, key, UInt64(0))
    if identifier != 0
        _write_emt_initialization_digest_tag!(writer, 0x00)
        return _write_emt_initialization_digest_uint64!(writer, identifier)
    end
    identifier = UInt64(length(writer.type_identifiers) + 1)
    writer.type_identifiers[key] = identifier
    _write_emt_initialization_digest_tag!(writer, 0x01)
    _write_emt_initialization_digest_uint64!(writer, identifier)
    _write_emt_initialization_digest_text!(writer, string(parentmodule(value)))
    _write_emt_initialization_digest_symbol!(writer, nameof(value))
    include_parameters || return writer
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(value.parameters)))
    for parameter in value.parameters
        if parameter isa DataType
            _write_emt_initialization_digest_type!(
                writer,
                parameter;
                include_parameters=false,
            )
        elseif parameter isa UnionAll
            _write_emt_initialization_digest_text!(writer, string(parameter))
        else
            _write_canonical_emt_initialization_state!(
                writer,
                parameter,
                IdDict{Any,Nothing}(),
            )
        end
    end
    return writer
end

function _emt_initialization_bulk_digest_eltype(::Type{T}) where {T}
    return T <: Union{
        Bool,
        Int8,
        Int16,
        Int32,
        Int64,
        Int128,
        UInt8,
        UInt16,
        UInt32,
        UInt64,
        UInt128,
        Float16,
        Float32,
        Float64,
        Complex{Float16},
        Complex{Float32},
        Complex{Float64},
    }
end

function _write_emt_initialization_digest_array_values!(
    writer::_EMTInitializationDigestWriter,
    value::AbstractArray,
    active_objects::IdDict{Any,Nothing},
)
    if value isa Array &&
       _emt_initialization_bulk_digest_eltype(eltype(value)) &&
       ENDIAN_BOM == 0x04030201
        _write_emt_initialization_digest_tag!(writer, 0x01)
        _update_emt_initialization_digest!(
            writer,
            reinterpret(UInt8, vec(value)),
        )
        return writer
    end
    _write_emt_initialization_digest_tag!(writer, 0x00)
    for index in eachindex(value)
        isassigned(value, index) || throw(ArgumentError(
            "deterministic EMT initialization state cannot contain an unassigned array entry",
        ))
        _write_canonical_emt_initialization_state!(
            writer,
            value[index],
            active_objects,
        )
    end
    return writer
end

function _write_emt_initialization_tuple_values_loop!(
    writer::_EMTInitializationDigestWriter,
    value::Tuple,
    active_objects::IdDict{Any,Nothing},
)
    for entry in value
        _write_canonical_emt_initialization_state!(
            writer,
            entry,
            active_objects,
        )
    end
    return writer
end

function _write_emt_initialization_homogeneous_tuple_values!(
    writer::_EMTInitializationDigestWriter,
    value::Tuple{Vararg{T}},
    active_objects::IdDict{Any,Nothing},
) where {T}
    for index in eachindex(value)
        _write_canonical_emt_initialization_state!(
            writer,
            value[index],
            active_objects,
        )
    end
    return writer
end

@generated function _write_emt_initialization_tuple_values!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:Tuple}
    field_count = fieldcount(T)
    if field_count > 32
        field_types = fieldtypes(T)
        if all(==(first(field_types)), field_types)
            return :(_write_emt_initialization_homogeneous_tuple_values!(
                writer,
                value,
                active_objects,
            ))
        elseif all(==(field_types[2]), field_types[2:end])
            return quote
                _write_canonical_emt_initialization_state!(
                    writer,
                    getfield(value, 1),
                    active_objects,
                )
                _write_emt_initialization_homogeneous_tuple_values!(
                    writer,
                    Base.tail(value),
                    active_objects,
                )
            end
        end
        return :(_write_emt_initialization_tuple_values_loop!(
            writer,
            value,
            active_objects,
        ))
    end
    expressions = [
        :(_write_canonical_emt_initialization_state!(
            writer,
            getfield(value, $index),
            active_objects,
        ))
        for index in 1:field_count
    ]
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_named_tuple_fields!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:NamedTuple}
    names = fieldnames(T)
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
            :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(value, $index),
                active_objects,
            )),
        )
    end
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_struct_fields!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T}
    names = fieldnames(T)
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
            :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(value, $index),
                active_objects,
            )),
        )
    end
    return quote
        $(expressions...)
        writer
    end
end

@generated function _write_emt_initialization_record_values!(
    writer::_EMTInitializationDigestWriter,
    value::T,
    active_objects::IdDict{Any,Nothing},
) where {T}
    expressions = [
        :(_write_canonical_emt_initialization_state!(
            writer,
            getfield(value, $index),
            active_objects,
        ))
        for index in 1:fieldcount(T)
    ]
    return quote
        $(expressions...)
        writer
    end
end

function _write_emt_initialization_record_schema!(
    writer::_EMTInitializationDigestWriter,
    ::Type{T},
) where {T}
    names = fieldnames(T)
    types = fieldtypes(T)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
    for (name, field_type) in zip(names, types)
        _write_emt_initialization_digest_symbol!(writer, name)
        if field_type isa DataType
            _write_emt_initialization_digest_type!(writer, field_type)
        else
            _write_emt_initialization_digest_text!(writer, string(field_type))
        end
    end
    return writer
end

function _write_emt_initialization_nodal_element_batch!(
    writer::_EMTInitializationDigestWriter,
    batch::AbstractVector{T},
    active_objects::IdDict{Any,Nothing},
) where {T}
    _write_emt_initialization_digest_type!(writer, T)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(batch)))
    compact_records = isbitstype(T) && fieldcount(T) > 0
    _write_emt_initialization_digest_tag!(
        writer,
        compact_records ? 0x01 : 0x00,
    )
    if compact_records
        _write_emt_initialization_record_schema!(writer, T)
        for element in batch
            _write_emt_initialization_record_values!(
                writer,
                element,
                active_objects,
            )
        end
    else
        for element in batch
            _write_canonical_emt_initialization_state!(
                writer,
                element,
                active_objects,
            )
        end
    end
    return writer
end

_write_emt_initialization_nodal_element_batches!(
    writer::_EMTInitializationDigestWriter,
    ::Tuple{},
    _active_objects::IdDict{Any,Nothing},
) = writer

function _write_emt_initialization_nodal_element_batches!(
    writer::_EMTInitializationDigestWriter,
    batches::Tuple,
    active_objects::IdDict{Any,Nothing},
)
    _write_emt_initialization_nodal_element_batch!(
        writer,
        first(batches),
        active_objects,
    )
    return _write_emt_initialization_nodal_element_batches!(
        writer,
        Base.tail(batches),
        active_objects,
    )
end

function _write_emt_initialization_nodal_elements!(
    writer::_EMTInitializationDigestWriter,
    elements::NodalElementSequence,
    active_objects::IdDict{Any,Nothing},
)
    batches = elements.contiguous_type_batches
    sum(length, batches; init=0) == length(elements) || throw(ArgumentError(
        "deterministic EMT initialization state found an inconsistent nodal element sequence",
    ))
    _write_emt_initialization_digest_tag!(writer, 0x19)
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(elements)))
    _write_emt_initialization_digest_uint64!(writer, UInt64(length(batches)))
    return _write_emt_initialization_nodal_element_batches!(
        writer,
        batches,
        active_objects,
    )
end

function _emt_initialization_arrays_bitwise_equal(
    first::Array{T},
    second::Array{T},
) where {T}
    size(first) == size(second) || return false
    isbitstype(T) || return false
    return reinterpret(UInt8, vec(first)) == reinterpret(UInt8, vec(second))
end

@generated function _write_emt_initialization_nodal_system_fields!(
    writer::_EMTInitializationDigestWriter,
    system::T,
    active_objects::IdDict{Any,Nothing},
) where {T<:NodalSystem}
    names = fieldnames(T)
    admittance_index = findfirst(==(:y), names)
    factor_index = findfirst(==(:y_factor), names)
    if isnothing(admittance_index) || isnothing(factor_index)
        return :(throw(ArgumentError(
            "deterministic EMT initialization requires nodal admittance and factor workspaces",
        )))
    end
    expressions = Expr[]
    for (index, name) in enumerate(names)
        push!(
            expressions,
            :(_write_emt_initialization_digest_symbol!(
                writer,
                $(QuoteNode(name)),
            )),
        )
        if index == factor_index
            push!(expressions, quote
                factor_matches_admittance =
                    _emt_initialization_arrays_bitwise_equal(
                        getfield(system, $admittance_index),
                        getfield(system, $factor_index),
                    )
                _write_emt_initialization_digest_tag!(
                    writer,
                    factor_matches_admittance ? 0x01 : 0x00,
                )
                factor_matches_admittance ||
                    _write_canonical_emt_initialization_state!(
                        writer,
                        getfield(system, $factor_index),
                        active_objects,
                    )
            end)
        else
            push!(expressions, :(_write_canonical_emt_initialization_state!(
                writer,
                getfield(system, $index),
                active_objects,
            )))
        end
    end
    return quote
        $(expressions...)
        writer
    end
end

function _write_emt_initialization_nodal_system!(
    writer::_EMTInitializationDigestWriter,
    system::NodalSystem,
    active_objects::IdDict{Any,Nothing},
)
    _write_emt_initialization_digest_tag!(writer, 0x18)
    _write_emt_initialization_digest_type!(writer, typeof(system))
    _write_emt_initialization_digest_uint64!(
        writer,
        UInt64(fieldcount(typeof(system))),
    )
    return _write_emt_initialization_nodal_system_fields!(
        writer,
        system,
        active_objects,
    )
end

function _emt_initialization_state_sort_key(value)
    return Tuple(_emt_initialization_state_digest(value))
end

function _write_canonical_emt_initialization_state!(
    writer::_EMTInitializationDigestWriter,
    value,
    active_objects::IdDict{Any,Nothing},
)
    if value === nothing
        return _write_emt_initialization_digest_tag!(writer, 0x00)
    elseif value === missing
        return _write_emt_initialization_digest_tag!(writer, 0x01)
    elseif value isa Bool
        _write_emt_initialization_digest_tag!(writer, 0x02)
        return _write_emt_initialization_digest_tag!(writer, value ? 0x01 : 0x00)
    elseif value isa Signed && sizeof(value) <= 8
        _write_emt_initialization_digest_tag!(writer, 0x03)
        return _write_emt_initialization_digest_uint64!(
            writer,
            reinterpret(UInt64, Int64(value)),
        )
    elseif value isa Unsigned && sizeof(value) <= 8
        _write_emt_initialization_digest_tag!(writer, 0x04)
        return _write_emt_initialization_digest_uint64!(writer, UInt64(value))
    elseif value isa Float64
        _write_emt_initialization_digest_tag!(writer, 0x05)
        return _write_emt_initialization_digest_uint64!(
            writer,
            reinterpret(UInt64, value),
        )
    elseif value isa Float32
        _write_emt_initialization_digest_tag!(writer, 0x06)
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(reinterpret(UInt32, value)),
        )
    elseif value isa Float16
        _write_emt_initialization_digest_tag!(writer, 0x07)
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(reinterpret(UInt16, value)),
        )
    elseif value isa AbstractFloat
        _write_emt_initialization_digest_tag!(writer, 0x08)
        return _write_emt_initialization_digest_text!(writer, repr(value))
    elseif value isa Char
        _write_emt_initialization_digest_tag!(writer, 0x09)
        return _write_emt_initialization_digest_uint64!(writer, UInt64(value))
    elseif value isa Symbol
        _write_emt_initialization_digest_tag!(writer, 0x0a)
        return _write_emt_initialization_digest_symbol!(writer, value)
    elseif value isa AbstractString
        _write_emt_initialization_digest_tag!(writer, 0x0b)
        return _write_emt_initialization_digest_text!(writer, value)
    elseif value isa Enum
        _write_emt_initialization_digest_tag!(writer, 0x0c)
        _write_emt_initialization_digest_type!(writer, typeof(value))
        return _write_emt_initialization_digest_uint64!(
            writer,
            UInt64(Integer(value)),
        )
    elseif value isa DataType
        _write_emt_initialization_digest_tag!(writer, 0x0d)
        return _write_emt_initialization_digest_type!(writer, value)
    elseif value isa Type
        _write_emt_initialization_digest_tag!(writer, 0x17)
        return _write_emt_initialization_digest_text!(writer, string(value))
    elseif value isa Module
        _write_emt_initialization_digest_tag!(writer, 0x0e)
        return _write_emt_initialization_digest_text!(writer, string(value))
    elseif value isa Complex
        _write_emt_initialization_digest_tag!(writer, 0x0f)
        _write_canonical_emt_initialization_state!(writer, real(value), active_objects)
        return _write_canonical_emt_initialization_state!(
            writer,
            imag(value),
            active_objects,
        )
    elseif value isa Ptr
        value == C_NULL || throw(ArgumentError(
            "deterministic EMT initialization state cannot contain a live pointer",
        ))
        return _write_emt_initialization_digest_tag!(writer, 0x10)
    end

    tracked = ismutabletype(typeof(value))
    if tracked
        haskey(active_objects, value) && throw(ArgumentError(
            "deterministic EMT initialization state cannot contain a reference cycle",
        ))
        active_objects[value] = nothing
    end
    try
        if value isa NodalSystem
            _write_emt_initialization_nodal_system!(
                writer,
                value,
                active_objects,
            )
        elseif value isa NodalElementSequence
            _write_emt_initialization_nodal_elements!(
                writer,
                value,
                active_objects,
            )
        elseif value isa NamedTuple
            _write_emt_initialization_digest_tag!(writer, 0x11)
            names = keys(value)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
            _write_emt_initialization_named_tuple_fields!(
                writer,
                value,
                active_objects,
            )
        elseif value isa Tuple
            _write_emt_initialization_digest_tag!(writer, 0x12)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(value)))
            _write_emt_initialization_tuple_values!(
                writer,
                value,
                active_objects,
            )
        elseif value isa AbstractArray
            _write_emt_initialization_digest_tag!(writer, 0x13)
            element_type = eltype(value)
            if element_type isa DataType
                _write_emt_initialization_digest_type!(
                    writer,
                    element_type;
                    include_parameters=isempty(value),
                )
            else
                _write_emt_initialization_digest_text!(writer, string(element_type))
            end
            _write_emt_initialization_digest_uint64!(writer, UInt64(ndims(value)))
            for dimension in size(value)
                _write_emt_initialization_digest_uint64!(writer, UInt64(dimension))
            end
            _write_emt_initialization_digest_array_values!(
                writer,
                value,
                active_objects,
            )
        elseif value isa AbstractDict
            _write_emt_initialization_digest_tag!(writer, 0x14)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(value)))
            if keytype(value) === Symbol
                ordered_keys = sort!(collect(keys(value)))
                for key in ordered_keys
                    _write_canonical_emt_initialization_state!(
                        writer,
                        key,
                        active_objects,
                    )
                    _write_canonical_emt_initialization_state!(
                        writer,
                        value[key],
                        active_objects,
                    )
                end
            else
                entries = [
                    (_emt_initialization_state_sort_key(key), key, entry)
                    for (key, entry) in value
                ]
                sort!(entries; by=first)
                for (_, key, entry) in entries
                    _write_canonical_emt_initialization_state!(
                        writer,
                        key,
                        active_objects,
                    )
                    _write_canonical_emt_initialization_state!(
                        writer,
                        entry,
                        active_objects,
                    )
                end
            end
        elseif value isa AbstractSet
            _write_emt_initialization_digest_tag!(writer, 0x15)
            entries = [
                (_emt_initialization_state_sort_key(entry), entry)
                for entry in value
            ]
            sort!(entries; by=first)
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(entries)))
            for (_, entry) in entries
                _write_canonical_emt_initialization_state!(writer, entry, active_objects)
            end
        else
            _write_emt_initialization_digest_tag!(writer, 0x16)
            _write_emt_initialization_digest_type!(writer, typeof(value))
            names = fieldnames(typeof(value))
            _write_emt_initialization_digest_uint64!(writer, UInt64(length(names)))
            if isempty(names)
                _write_emt_initialization_digest_text!(writer, repr(value))
            else
                _write_emt_initialization_struct_fields!(
                    writer,
                    value,
                    active_objects,
                )
            end
        end
    finally
        tracked && delete!(active_objects, value)
    end
    return writer
end

function _emt_initialization_state_digest(value)
    writer = _EMTInitializationDigestWriter()
    _write_canonical_emt_initialization_state!(
        writer,
        value,
        IdDict{Any,Nothing}(),
    )
    _flush_emt_initialization_digest!(writer)
    return SHA.digest!(writer.context)
end

function _emt_initialization_state_signature(
    prepared::Union{
        PreparedEMTStudy,
        PreparedMachineEMTStudy,
        PreparedAverageValueGridFollowingEMTStudy,
    },
    request::EMTInitializationRequest,
    mappings::Vector{OperatingPointMappingRecord},
)
    io = IOBuffer()
    write(io, "aimora.emt.initialization.state.v4\n")
    for value in (
        request.project_signature,
        request.settings_signature,
        request.model_signature,
        String(_emt_harmonic_formulation_symbol(request.formulation)),
        repr(request.frequency_hz),
        repr(request.time_origin_s),
    )
        write(io, value, '\n')
    end
    if request.operating_point isa EMTOperatingPoint
        write(
            io,
            String(request.operating_point.source_representation),
            '\n',
            request.operating_point.source_state_signature,
            '\n',
        )
    end
    accepted_state = if prepared isa PreparedEMTStudy
        prepared.runtime_template
    elseif prepared isa PreparedMachineEMTStudy
        (
            machine_family=prepared.machine_family,
            initialization_state=prepared.initialization_state,
        )
    else
        (
            network=prepared.network.runtime_template,
            converter_equilibrium=prepared.equilibrium,
            converter_ownership=prepared.converter,
        )
    end
    write(
        io,
        "complete_accepted_state_sha256=",
        bytes2hex(_emt_initialization_state_digest(accepted_state)),
        '\n',
    )
    for mapping in mappings
        write(
            io,
            String(mapping.asset),
            '/',
            String(mapping.quantity),
            '/',
            String(mapping.phase),
            '=',
            bitstring(real(mapping.target_value_si)),
            ',',
            bitstring(imag(mapping.target_value_si)),
            '\n',
        )
    end
    return bytes2hex(sha256(take!(io)))
end

function _protection_portable_type_name(type)
    name = String(nameof(type))
    snake = lowercase(replace(name, r"([a-z0-9])([A-Z])" => s"\1_\2"))
    occursin(r"^[a-z][a-z0-9_]*$", snake) || throw(ArgumentError(
        "protection portable type name is not a stable identity",
    ))
    return snake
end

_protection_portable_schema(type) =
    "aimora.protection.$(_protection_portable_type_name(type)).v1"

function _protection_portable_value(value)
    value === nothing && return nothing
    value isa Bool && return value
    value isa Signed && return Int64(value)
    value isa Unsigned && return UInt64(value)
    value isa AbstractFloat && return Float64(value)
    value isa AbstractString && return String(value)
    value isa Symbol && return PortableSnapshots.PortableSnapshotRecord(
        "aimora.protection.symbol.v1",
        Pair{String,Any}["value" => String(value)],
    )
    value isa Enum && return PortableSnapshots.PortableSnapshotRecord(
        _protection_portable_schema(typeof(value)),
        Pair{String,Any}["value" => Int64(UInt8(value))],
    )
    value isa Complex && return PortableSnapshots.PortableSnapshotRecord(
        "aimora.protection.complex.v1",
        Pair{String,Any}[
            "imaginary" => Float64(imag(value)),
            "real" => Float64(real(value)),
        ],
    )
    if value isa Tuple || value isa AbstractVector
        return Any[_protection_portable_value(item) for item in value]
    end
    if value isa AbstractDict
        records = PortableSnapshots.PortableSnapshotRecord[]
        for key in sort!(collect(keys(value)); by=repr)
            push!(records, PortableSnapshots.PortableSnapshotRecord(
                "aimora.protection.dictionary_entry.v1",
                Pair{String,Any}[
                    "key" => _protection_portable_value(key),
                    "value" => _protection_portable_value(value[key]),
                ],
            ))
        end
        return records
    end
    type = typeof(value)
    isstructtype(type) || throw(ArgumentError(
        "protection portable snapshot cannot encode $(type)",
    ))
    return PortableSnapshots.PortableSnapshotRecord(
        _protection_portable_schema(type),
        Pair{String,Any}[
            String(field) => _protection_portable_value(getfield(value, field)) for
            field in fieldnames(type)
        ],
    )
end

function _protection_portable_record_values(record, schema)
    record isa PortableSnapshots.PortableSnapshotRecord && record.schema_id == schema ||
        throw(ArgumentError("protection portable snapshot record schema changed"))
    return Dict(record.fields)
end

function _restore_protection_portable_value(encoded, type)
    if type isa Union
        encoded === nothing && Nothing <: type && return nothing
        candidates = filter(candidate -> candidate !== Nothing, Base.uniontypes(type))
        length(candidates) == 1 || throw(ArgumentError(
            "protection portable snapshot union is ambiguous",
        ))
        return _restore_protection_portable_value(encoded, only(candidates))
    end
    type === Nothing && return nothing
    type === Bool && return Bool(encoded)
    type <: Signed && return type(encoded)
    type <: Unsigned && return type(encoded)
    type <: AbstractFloat && return type(encoded)
    type <: AbstractString && return type(encoded)
    if type === Symbol
        values = _protection_portable_record_values(
            encoded,
            "aimora.protection.symbol.v1",
        )
        Set(keys(values)) == Set(["value"]) || throw(ArgumentError(
            "protection portable symbol record changed",
        ))
        return Symbol(values["value"])
    end
    if type <: Enum
        values = _protection_portable_record_values(
            encoded,
            _protection_portable_schema(type),
        )
        Set(keys(values)) == Set(["value"]) || throw(ArgumentError(
            "protection portable enum record changed",
        ))
        return type(values["value"])
    end
    if type <: Complex
        values = _protection_portable_record_values(
            encoded,
            "aimora.protection.complex.v1",
        )
        Set(keys(values)) == Set(["imaginary", "real"]) || throw(ArgumentError(
            "protection portable complex record changed",
        ))
        return type(complex(values["real"], values["imaginary"]))
    end
    if type <: Tuple
        encoded isa AbstractVector || throw(ArgumentError(
            "protection portable tuple has the wrong representation",
        ))
        types = fieldtypes(type)
        length(encoded) == length(types) || throw(ArgumentError(
            "protection portable tuple length changed",
        ))
        return tuple((
            _restore_protection_portable_value(encoded[index], types[index]) for
            index in eachindex(types)
        )...)
    end
    if type <: AbstractVector
        encoded isa AbstractVector || throw(ArgumentError(
            "protection portable vector has the wrong representation",
        ))
        element_type = eltype(type)
        return element_type[
            _restore_protection_portable_value(item, element_type) for item in encoded
        ]
    end
    if type <: AbstractDict
        encoded isa AbstractVector || throw(ArgumentError(
            "protection portable dictionary has the wrong representation",
        ))
        key_type, value_type = type.parameters
        restored = type()
        for record in encoded
            values = _protection_portable_record_values(
                record,
                "aimora.protection.dictionary_entry.v1",
            )
            Set(keys(values)) == Set(["key", "value"]) || throw(ArgumentError(
                "protection portable dictionary entry changed",
            ))
            key = _restore_protection_portable_value(values["key"], key_type)
            haskey(restored, key) && throw(ArgumentError(
                "protection portable dictionary repeats a key",
            ))
            restored[key] = _restore_protection_portable_value(
                values["value"],
                value_type,
            )
        end
        return restored
    end
    values = _protection_portable_record_values(
        encoded,
        _protection_portable_schema(type),
    )
    fields = fieldnames(type)
    Set(keys(values)) == Set(String.(fields)) || throw(ArgumentError(
        "protection portable $(type) fields changed",
    ))
    restored_fields = Any[
        _restore_protection_portable_value(values[String(field)], fieldtype(type, index)) for
        (index, field) in enumerate(fields)
    ]
    return type(restored_fields...)
end

function protection_product_portable_snapshot(
    runtime::ProtectionProductRuntime;
    project_signature_sha256::AbstractString,
    topology_signature_sha256::AbstractString,
)
    typed_snapshot = protection_product_runtime_snapshot(runtime)
    specification = runtime.preparation.specification
    represented_time = runtime.accepted_tick *
        (specification.network_timestep_logical.numerator //
         specification.network_timestep_logical.denominator)
    metadata = PortableSnapshots.PortableSnapshotMetadata(
        :portable_public_reference,
        project_signature_sha256,
        specification.deterministic_signature_sha256,
        topology_signature_sha256,
        runtime.preparation.preparation_signature_sha256,
        represented_time,
        runtime.accepted_tick,
        ["emt.protection", "emt.tasks", "emt.breaker", "emt.snapshot"],
        "AIMORA-authored generic public protection product runtime state",
    )
    section = PortableSnapshots.PortableSnapshotSection(
        "protection.product_runtime",
        1,
        0,
        :public,
        _protection_portable_value(typed_snapshot),
    )
    return PortableSnapshots.PortableEMTSnapshot(metadata, [section])
end

function write_protection_product_portable_snapshot(
    path::AbstractString,
    runtime::ProtectionProductRuntime;
    project_signature_sha256::AbstractString,
    topology_signature_sha256::AbstractString,
)
    snapshot = protection_product_portable_snapshot(
        runtime;
        project_signature_sha256,
        topology_signature_sha256,
    )
    return PortableSnapshots.write_portable_emt_snapshot(path, snapshot)
end

function restore_protection_product_portable_snapshot!(
    runtime::ProtectionProductRuntime,
    path::AbstractString;
    project_signature_sha256::AbstractString,
    topology_signature_sha256::AbstractString,
)
    snapshot, descriptor = PortableSnapshots.read_portable_emt_snapshot_with_descriptor(path)
    specification = runtime.preparation.specification
    metadata = snapshot.metadata
    metadata.profile === :portable_public_reference || throw(ArgumentError(
        "protection portable snapshot profile changed",
    ))
    metadata.project_signature_sha256 == lowercase(String(project_signature_sha256)) ||
        throw(ArgumentError("protection portable snapshot project identity is stale"))
    metadata.model_signature_sha256 == specification.deterministic_signature_sha256 ||
        throw(ArgumentError("protection portable snapshot product identity is stale"))
    metadata.topology_signature_sha256 == lowercase(String(topology_signature_sha256)) ||
        throw(ArgumentError("protection portable snapshot topology identity is stale"))
    metadata.settings_signature_sha256 ==
        runtime.preparation.preparation_signature_sha256 || throw(ArgumentError(
        "protection portable snapshot preparation identity is stale",
    ))
    length(snapshot.sections) == 1 &&
        only(snapshot.sections).identity == "protection.product_runtime" ||
        throw(ArgumentError("protection portable snapshot section inventory changed"))
    template = protection_product_runtime_snapshot(runtime)
    restored_snapshot = _restore_protection_portable_value(
        only(snapshot.sections).value,
        typeof(template),
    )
    restore_protection_product_runtime_snapshot!(runtime, restored_snapshot)
    runtime.accepted_tick == metadata.accepted_step || throw(ArgumentError(
        "protection portable snapshot accepted tick changed",
    ))
    return runtime, descriptor
end

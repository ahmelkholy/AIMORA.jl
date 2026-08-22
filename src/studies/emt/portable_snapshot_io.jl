module PortableSnapshots

using SHA
using Unicode

export PortableSnapshotFailure,
       PortableSnapshotArray,
       PortableSnapshotRecord,
       PortableSnapshotSection,
       PortableSnapshotMetadata,
       PortableEMTSnapshot,
       PortableSnapshotSectionDescriptor,
       PortableSnapshotDescriptor,
       PortableSnapshotStateField,
       PortableSnapshotStateInventory,
       portable_snapshot_array,
       portable_snapshot_array_values,
       portable_state_inventory_record,
       portable_state_inventory,
       portable_snapshot_bytes,
       portable_snapshot_descriptor,
       inspect_portable_emt_snapshot,
       read_portable_emt_snapshot,
       read_portable_emt_snapshot_with_descriptor,
       write_portable_emt_snapshot

const _PORTABLE_SNAPSHOT_MAGIC = collect(codeunits("AIMORA-PORTABLE-EMT"))
const _PORTABLE_SNAPSHOT_SCHEMA_MAJOR = UInt16(1)
const _PORTABLE_SNAPSHOT_SCHEMA_MINOR = UInt16(1)
const _PORTABLE_SNAPSHOT_DIGEST_BYTES = 32
const _PORTABLE_SNAPSHOT_ID = r"^[a-z][a-z0-9]*(?:[._-][a-z0-9]+)*$"

const _SNAPSHOT_NULL = UInt8(0)
const _SNAPSHOT_FALSE = UInt8(1)
const _SNAPSHOT_TRUE = UInt8(2)
const _SNAPSHOT_SIGNED = UInt8(3)
const _SNAPSHOT_UNSIGNED = UInt8(4)
const _SNAPSHOT_FLOAT64 = UInt8(5)
const _SNAPSHOT_TEXT = UInt8(6)
const _SNAPSHOT_BYTES = UInt8(7)
const _SNAPSHOT_SEQUENCE = UInt8(8)
const _SNAPSHOT_RECORD = UInt8(9)
const _SNAPSHOT_ARRAY = UInt8(10)
const _SNAPSHOT_RATIONAL = UInt8(11)

function _write_checkpoint_integer(io::IO, value::T) where {T<:Unsigned}
    for shift in 0:8:(8 * sizeof(T) - 8)
        write(io, UInt8((value >> shift) & T(0xff)))
    end
    return nothing
end

function _read_checkpoint_integer(io::IO, ::Type{T}) where {T<:Unsigned}
    value = zero(T)
    for shift in 0:8:(8 * sizeof(T) - 8)
        eof(io) && _portable_snapshot_fail(
            :truncated,
            "portable snapshot integer is truncated",
        )
        value |= T(read(io, UInt8)) << shift
    end
    return value
end

struct PortableSnapshotFailure <: Exception
    code::Symbol
    message::String
end

Base.showerror(io::IO, failure::PortableSnapshotFailure) =
    print(io, "portable EMT snapshot ", failure.code, ": ", failure.message)

function _portable_snapshot_fail(code::Symbol, message::AbstractString)
    throw(PortableSnapshotFailure(code, String(message)))
end

function _portable_snapshot_identity(value::AbstractString, label::AbstractString)
    identity = String(value)
    occursin(_PORTABLE_SNAPSHOT_ID, identity) || _portable_snapshot_fail(
        :invalid_identity,
        "$label must be a lowercase stable semantic identity",
    )
    return identity
end

function _portable_snapshot_text(value::AbstractString, label::AbstractString)
    text = String(value)
    Unicode.normalize(text, :NFC) == text || _portable_snapshot_fail(
        :noncanonical_text,
        "$label must use normalized UTF-8 text",
    )
    return text
end

function _portable_snapshot_sha256(value::AbstractString, label::AbstractString)
    digest = String(value)
    occursin(r"^[0-9a-f]{64}$", digest) || _portable_snapshot_fail(
        :invalid_digest,
        "$label must be a lowercase SHA-256 identity",
    )
    return digest
end

struct PortableSnapshotArray
    element_kind::Symbol
    shape::Vector{Int}
    unit::String
    axes::Vector{String}
    bytes::Vector{UInt8}

    function PortableSnapshotArray(
        element_kind::Symbol,
        shape::AbstractVector{<:Integer},
        unit::AbstractString,
        axes::AbstractVector{<:AbstractString},
        bytes::AbstractVector{UInt8},
    )
        element_kind in (:float64, :int64, :uint64) || _portable_snapshot_fail(
            :unsupported_array_element,
            "portable array element kind $element_kind is unsupported",
        )
        dimensions = Int[]
        for dimension in shape
            0 <= dimension <= typemax(Int) || _portable_snapshot_fail(
                :invalid_array_shape,
                "portable array dimensions must be nonnegative and addressable",
            )
            push!(dimensions, Int(dimension))
        end
        axis_names = String[
            _portable_snapshot_identity(axis, "portable array axis") for axis in axes
        ]
        length(axis_names) == length(dimensions) || _portable_snapshot_fail(
            :invalid_array_axes,
            "portable array axis count must equal its rank",
        )
        length(unique(axis_names)) == length(axis_names) || _portable_snapshot_fail(
            :duplicate_array_axis,
            "portable array axes must be unique",
        )
        element_count = 1
        for dimension in dimensions
            element_count = try
                Base.checked_mul(element_count, dimension)
            catch error
                error isa OverflowError || rethrow()
                _portable_snapshot_fail(
                    :invalid_array_shape,
                    "portable array element count exceeds the addressable range",
                )
            end
        end
        expected_bytes = try
            Base.checked_mul(element_count, 8)
        catch error
            error isa OverflowError || rethrow()
            _portable_snapshot_fail(
                :invalid_array_shape,
                "portable array payload length exceeds the addressable range",
            )
        end
        length(bytes) == expected_bytes || _portable_snapshot_fail(
            :invalid_array_payload,
            "portable array payload length does not match its shape",
        )
        if element_kind == :float64
            payload = IOBuffer(bytes)
            for _ in 1:element_count
                value = reinterpret(Float64, _read_checkpoint_integer(payload, UInt64))
                isfinite(value) || _portable_snapshot_fail(
                    :nonfinite_value,
                    "portable physical floating array contains a nonfinite value",
                )
            end
        end
        return new(
            element_kind,
            dimensions,
            _portable_snapshot_text(unit, "portable array unit"),
            axis_names,
            Vector{UInt8}(bytes),
        )
    end
end

struct PortableSnapshotRecord
    schema_id::String
    fields::Vector{Pair{String,Any}}

    function PortableSnapshotRecord(
        schema_id::AbstractString,
        fields::AbstractVector{<:Pair},
    )
        identity = _portable_snapshot_identity(schema_id, "portable record schema")
        normalized = Pair{String,Any}[
            _portable_snapshot_identity(String(first(field)), "portable record field") =>
                last(field)
            for field in fields
        ]
        sort!(normalized; by = first)
        names = first.(normalized)
        length(unique(names)) == length(names) || _portable_snapshot_fail(
            :duplicate_record_field,
            "portable record $identity repeats a field identity",
        )
        return new(identity, normalized)
    end
end

struct PortableSnapshotSection
    identity::String
    version_major::UInt16
    version_minor::UInt16
    visibility::Symbol
    value::Any

    function PortableSnapshotSection(
        identity::AbstractString,
        version_major::Integer,
        version_minor::Integer,
        visibility::Symbol,
        value,
    )
        0 <= version_major <= typemax(UInt16) || _portable_snapshot_fail(
            :invalid_section_version,
            "portable section major version is outside UInt16",
        )
        0 <= version_minor <= typemax(UInt16) || _portable_snapshot_fail(
            :invalid_section_version,
            "portable section minor version is outside UInt16",
        )
        visibility in (:public, :private_reconstructible) || _portable_snapshot_fail(
            :invalid_section_visibility,
            "portable section visibility must be public or private_reconstructible",
        )
        return new(
            _portable_snapshot_identity(identity, "portable section"),
            UInt16(version_major),
            UInt16(version_minor),
            visibility,
            value,
        )
    end
end

struct PortableSnapshotMetadata
    profile::Symbol
    project_signature_sha256::String
    model_signature_sha256::String
    topology_signature_sha256::String
    settings_signature_sha256::String
    represented_time_s::Rational{Int128}
    accepted_step::Int64
    capabilities::Vector{String}
    provenance::String
    writer_version::String
    creator_platform::String
    numeric_profile::String
    compression::String
    minimum_reader_version::String

    function PortableSnapshotMetadata(
        profile::Symbol,
        project_signature_sha256::AbstractString,
        model_signature_sha256::AbstractString,
        topology_signature_sha256::AbstractString,
        settings_signature_sha256::AbstractString,
        represented_time_s::Rational,
        accepted_step::Integer,
        capabilities::AbstractVector{<:AbstractString},
        provenance::AbstractString;
        writer_version::AbstractString = "AIMORA.jl/0.1.0",
        creator_platform::AbstractString = "$(Sys.KERNEL)/$(Sys.ARCH)/Julia-$(VERSION)",
        numeric_profile::AbstractString = "ieee754_binary64_finite_preserve_signed_zero",
        compression::AbstractString = "none",
        minimum_reader_version::AbstractString = "1.0",
    )
        profile in (:portable_full, :portable_public_reference) ||
            _portable_snapshot_fail(:unsupported_profile, "portable snapshot profile is unsupported")
        0 <= accepted_step <= typemax(Int64) || _portable_snapshot_fail(
            :invalid_accepted_step,
            "portable snapshot accepted step must be nonnegative and fit Int64",
        )
        numerator_big = BigInt(numerator(represented_time_s))
        denominator_big = BigInt(denominator(represented_time_s))
        typemin(Int128) <= numerator_big <= typemax(Int128) ||
            _portable_snapshot_fail(:time_overflow, "portable snapshot time numerator exceeds Int128")
        0 < denominator_big <= typemax(Int128) ||
            _portable_snapshot_fail(:time_overflow, "portable snapshot time denominator exceeds Int128")
        capability_ids = sort!(unique(String[
            _portable_snapshot_identity(value, "portable snapshot capability")
            for value in capabilities
        ]))
        numeric_profile == "ieee754_binary64_finite_preserve_signed_zero" ||
            _portable_snapshot_fail(
                :unsupported_numeric_profile,
                "portable snapshot numeric profile is unsupported",
            )
        compression == "none" || _portable_snapshot_fail(
            :unsupported_compression,
            "portable snapshot compression profile is unsupported",
        )
        minimum_reader_version in ("1.0", "1.1") || _portable_snapshot_fail(
            :invalid_reader_version,
            "portable snapshot minimum reader version is unsupported",
        )
        return new(
            profile,
            _portable_snapshot_sha256(project_signature_sha256, "project signature"),
            _portable_snapshot_sha256(model_signature_sha256, "model signature"),
            _portable_snapshot_sha256(topology_signature_sha256, "topology signature"),
            _portable_snapshot_sha256(settings_signature_sha256, "settings signature"),
            Int128(numerator_big) // Int128(denominator_big),
            Int64(accepted_step),
            capability_ids,
            _portable_snapshot_text(provenance, "portable snapshot provenance"),
            _portable_snapshot_text(writer_version, "portable snapshot writer version"),
            _portable_snapshot_text(creator_platform, "portable snapshot creator platform"),
            String(numeric_profile),
            String(compression),
            String(minimum_reader_version),
        )
    end
end

struct PortableEMTSnapshot
    metadata::PortableSnapshotMetadata
    sections::Vector{PortableSnapshotSection}

    function PortableEMTSnapshot(
        metadata::PortableSnapshotMetadata,
        sections::AbstractVector{PortableSnapshotSection},
    )
        normalized = sort!(collect(sections); by = section -> section.identity)
        identities = getfield.(normalized, :identity)
        length(unique(identities)) == length(identities) || _portable_snapshot_fail(
            :duplicate_section,
            "portable snapshot repeats a section identity",
        )
        isempty(normalized) && _portable_snapshot_fail(
            :missing_section,
            "portable snapshot requires at least one scientific state section",
        )
        metadata.profile == :portable_public_reference &&
            any(section -> section.visibility != :public, normalized) &&
            _portable_snapshot_fail(
                :private_section_in_public_profile,
                "public-reference snapshot cannot contain private reconstructible sections",
            )
        return new(metadata, normalized)
    end
end

struct PortableSnapshotSectionDescriptor
    identity::String
    version_major::UInt16
    version_minor::UInt16
    visibility::Symbol
    payload_bytes::UInt64
    payload_sha256::String
end

function _portable_section_digest(
    identity::AbstractString,
    version_major::UInt16,
    version_minor::UInt16,
    visibility::Symbol,
    payload::AbstractVector{UInt8},
)
    framed = IOBuffer()
    _write_portable_string(framed, identity)
    _write_checkpoint_integer(framed, version_major)
    _write_checkpoint_integer(framed, version_minor)
    write(framed, visibility == :public ? UInt8(0) : UInt8(1))
    _write_checkpoint_integer(framed, UInt64(length(payload)))
    write(framed, payload)
    return sha256(take!(framed))
end

struct PortableSnapshotDescriptor
    schema_major::UInt16
    schema_minor::UInt16
    metadata::PortableSnapshotMetadata
    sections::Vector{PortableSnapshotSectionDescriptor}
    canonical_bytes::UInt64
    content_sha256::String
    migration_path::Vector{String}
end

const _PORTABLE_STATE_FAMILIES = Set((
    :continuous,
    :algebraic,
    :discrete,
    :delayed,
    :scheduler,
    :random,
    :history,
    :output,
    :checkpoint,
    :reconstruction,
))

struct PortableSnapshotStateField
    identity::String
    owner::String
    family::Symbol
    unit::String
    axes::Vector{String}
    value::Any

    function PortableSnapshotStateField(
        identity::AbstractString,
        owner::AbstractString,
        family::Symbol,
        unit::AbstractString,
        axes::AbstractVector{<:AbstractString},
        value,
    )
        family in _PORTABLE_STATE_FAMILIES || _portable_snapshot_fail(
            :invalid_state_family,
            "portable state family $family is unsupported",
        )
        axis_names = String[
            _portable_snapshot_identity(axis, "portable state axis") for axis in axes
        ]
        length(unique(axis_names)) == length(axis_names) || _portable_snapshot_fail(
            :duplicate_state_axis,
            "portable state field axes must be unique",
        )
        normalized_unit = _portable_snapshot_text(unit, "portable state unit")
        if value isa PortableSnapshotArray
            value.unit == normalized_unit || _portable_snapshot_fail(
                :state_array_unit,
                "portable state field and array units disagree",
            )
            value.axes == axis_names || _portable_snapshot_fail(
                :state_array_axes,
                "portable state field and array axes disagree",
            )
        end
        _portable_value_bytes(value)
        return new(
            _portable_snapshot_identity(identity, "portable state field"),
            _portable_snapshot_identity(owner, "portable state owner"),
            family,
            normalized_unit,
            axis_names,
            value,
        )
    end
end

function _portable_state_scalar_count(value)
    if value === nothing || value isa AbstractString || value isa AbstractVector{UInt8}
        return 0
    elseif value isa Number || value isa Bool
        return 1
    elseif value isa PortableSnapshotArray
        return isempty(value.shape) ? 1 : prod(value.shape; init = 1)
    elseif value isa PortableSnapshotRecord
        return sum(_portable_state_scalar_count(last(field)) for field in value.fields; init = 0)
    elseif value isa AbstractVector || value isa Tuple
        return sum(_portable_state_scalar_count(item) for item in value; init = 0)
    end
    _portable_snapshot_fail(
        :unsupported_state_value,
        "portable state value type $(typeof(value)) cannot be counted",
    )
end

struct PortableSnapshotStateInventory
    fields::Vector{PortableSnapshotStateField}
    scalar_count::UInt64
    signature_sha256::String

    function PortableSnapshotStateInventory(
        fields::AbstractVector{PortableSnapshotStateField},
    )
        normalized = sort!(collect(fields); by = field -> field.identity)
        identities = getfield.(normalized, :identity)
        length(unique(identities)) == length(identities) || _portable_snapshot_fail(
            :duplicate_state_field,
            "portable state inventory repeats a field identity",
        )
        isempty(normalized) && _portable_snapshot_fail(
            :missing_state_field,
            "portable state inventory requires at least one registered field",
        )
        scalar_count = sum(
            _portable_state_scalar_count(field.value) for field in normalized;
            init = 0,
        )
        scalar_count <= typemax(UInt64) || _portable_snapshot_fail(
            :state_count_overflow,
            "portable state scalar count exceeds UInt64",
        )
        unsigned_count = UInt64(scalar_count)
        provisional = new(normalized, unsigned_count, "")
        signature = bytes2hex(sha256(_portable_value_bytes(
            _portable_state_inventory_record(provisional; include_signature = false),
        )))
        return new(normalized, unsigned_count, signature)
    end
end

function _portable_state_field_record(field::PortableSnapshotStateField)
    return PortableSnapshotRecord(
        "aimora.snapshot.state_field.v1",
        Pair{String,Any}[
            "axes" => field.axes,
            "family" => String(field.family),
            "identity" => field.identity,
            "owner" => field.owner,
            "unit" => field.unit,
            "value" => field.value,
        ],
    )
end

function _portable_state_inventory_record(
    inventory::PortableSnapshotStateInventory;
    include_signature::Bool,
)
    fields = Pair{String,Any}[
        "fields" => _portable_state_field_record.(inventory.fields),
        "scalar_count" => inventory.scalar_count,
    ]
    include_signature && push!(fields, "signature_sha256" => inventory.signature_sha256)
    return PortableSnapshotRecord("aimora.snapshot.state_inventory.v1", fields)
end

portable_state_inventory_record(inventory::PortableSnapshotStateInventory) =
    _portable_state_inventory_record(inventory; include_signature = true)

function portable_state_inventory(record::PortableSnapshotRecord)
    fields = _record_dictionary(record, "aimora.snapshot.state_inventory.v1")
    Set(keys(fields)) == Set(("fields", "scalar_count", "signature_sha256")) ||
        _portable_snapshot_fail(
            :state_inventory_fields,
            "portable state inventory fields are incomplete or unknown",
        )
    encoded_fields = fields["fields"]
    encoded_fields isa AbstractVector || _portable_snapshot_fail(
        :state_inventory_type,
        "portable state inventory field list has the wrong type",
    )
    state_fields = PortableSnapshotStateField[]
    for encoded in encoded_fields
        encoded isa PortableSnapshotRecord || _portable_snapshot_fail(
            :state_inventory_type,
            "portable state inventory entry is not a record",
        )
        values = _record_dictionary(encoded, "aimora.snapshot.state_field.v1")
        Set(keys(values)) == Set(("axes", "family", "identity", "owner", "unit", "value")) ||
            _portable_snapshot_fail(
                :state_inventory_fields,
                "portable state field is incomplete or unknown",
            )
        axes = values["axes"]
        axes isa AbstractVector && all(axis -> axis isa AbstractString, axes) ||
            _portable_snapshot_fail(:state_inventory_type, "portable state axes have the wrong type")
        push!(state_fields, PortableSnapshotStateField(
            values["identity"],
            values["owner"],
            Symbol(values["family"]),
            values["unit"],
            String[axis for axis in axes],
            values["value"],
        ))
    end
    inventory = PortableSnapshotStateInventory(state_fields)
    fields["scalar_count"] == inventory.scalar_count || _portable_snapshot_fail(
        :state_inventory_count,
        "portable state inventory scalar count is inconsistent",
    )
    fields["signature_sha256"] == inventory.signature_sha256 || _portable_snapshot_fail(
        :state_inventory_integrity,
        "portable state inventory signature is inconsistent",
    )
    return inventory
end

function _portable_array_kind(::Type{Float64})
    return :float64
end
function _portable_array_kind(::Type{Int64})
    return :int64
end
function _portable_array_kind(::Type{UInt64})
    return :uint64
end

function portable_snapshot_array(
    values::AbstractArray{T};
    unit::AbstractString,
    axes::AbstractVector{<:AbstractString},
) where {T<:Union{Float64,Int64,UInt64}}
    io = IOBuffer(; sizehint = Base.checked_mul(length(values), 8))
    for value in values
        encoded = if T === Float64
            reinterpret(UInt64, value)
        elseif T === Int64
            reinterpret(UInt64, value)
        else
            value
        end
        _write_checkpoint_integer(io, encoded)
    end
    return PortableSnapshotArray(_portable_array_kind(T), collect(size(values)), unit, axes, take!(io))
end

function portable_snapshot_array_values(array::PortableSnapshotArray)
    io = IOBuffer(array.bytes)
    count = isempty(array.shape) ? 1 : prod(array.shape; init = 1)
    if array.element_kind == :float64
        values = Float64[reinterpret(Float64, _read_checkpoint_integer(io, UInt64)) for _ in 1:count]
    elseif array.element_kind == :int64
        values = Int64[reinterpret(Int64, _read_checkpoint_integer(io, UInt64)) for _ in 1:count]
    else
        values = UInt64[_read_checkpoint_integer(io, UInt64) for _ in 1:count]
    end
    return reshape(values, Tuple(array.shape))
end

function _write_portable_string(io::IO, value::AbstractString)
    bytes = codeunits(value)
    length(bytes) <= typemax(UInt32) || _portable_snapshot_fail(
        :text_too_large,
        "portable snapshot text exceeds UInt32 framing",
    )
    _write_checkpoint_integer(io, UInt32(length(bytes)))
    write(io, bytes)
    return nothing
end

function _write_portable_value(io::IO, value)
    if value === nothing
        write(io, _SNAPSHOT_NULL)
    elseif value === false
        write(io, _SNAPSHOT_FALSE)
    elseif value === true
        write(io, _SNAPSHOT_TRUE)
    elseif value isa Signed
        typemin(Int64) <= value <= typemax(Int64) ||
            _portable_snapshot_fail(:integer_overflow, "portable signed value exceeds Int64")
        write(io, _SNAPSHOT_SIGNED)
        _write_checkpoint_integer(io, reinterpret(UInt64, Int64(value)))
    elseif value isa Unsigned
        value <= typemax(UInt64) ||
            _portable_snapshot_fail(:integer_overflow, "portable unsigned value exceeds UInt64")
        write(io, _SNAPSHOT_UNSIGNED)
        _write_checkpoint_integer(io, UInt64(value))
    elseif value isa AbstractFloat
        isfinite(value) || _portable_snapshot_fail(
            :nonfinite_value,
            "portable physical floating value must be finite",
        )
        write(io, _SNAPSHOT_FLOAT64)
        _write_checkpoint_integer(io, reinterpret(UInt64, Float64(value)))
    elseif value isa AbstractString
        write(io, _SNAPSHOT_TEXT)
        _write_portable_string(io, _portable_snapshot_text(value, "portable value"))
    elseif value isa AbstractVector{UInt8}
        write(io, _SNAPSHOT_BYTES)
        _write_checkpoint_integer(io, UInt64(length(value)))
        write(io, value)
    elseif value isa Rational
        numerator_big = BigInt(numerator(value))
        denominator_big = BigInt(denominator(value))
        typemin(Int128) <= numerator_big <= typemax(Int128) ||
            _portable_snapshot_fail(:integer_overflow, "portable rational numerator exceeds Int128")
        0 < denominator_big <= typemax(Int128) ||
            _portable_snapshot_fail(:integer_overflow, "portable rational denominator exceeds Int128")
        write(io, _SNAPSHOT_RATIONAL)
        _write_checkpoint_integer(io, reinterpret(UInt128, Int128(numerator_big)))
        _write_checkpoint_integer(io, UInt128(denominator_big))
    elseif value isa PortableSnapshotArray
        write(io, _SNAPSHOT_ARRAY)
        _write_portable_string(io, String(value.element_kind))
        _write_checkpoint_integer(io, UInt32(length(value.shape)))
        for dimension in value.shape
            _write_checkpoint_integer(io, UInt64(dimension))
        end
        _write_portable_string(io, value.unit)
        for axis in value.axes
            _write_portable_string(io, axis)
        end
        _write_checkpoint_integer(io, UInt64(length(value.bytes)))
        write(io, value.bytes)
    elseif value isa PortableSnapshotRecord
        write(io, _SNAPSHOT_RECORD)
        _write_portable_string(io, value.schema_id)
        _write_checkpoint_integer(io, UInt32(length(value.fields)))
        for field in value.fields
            _write_portable_string(io, first(field))
            _write_portable_value(io, last(field))
        end
    elseif value isa AbstractVector || value isa Tuple
        write(io, _SNAPSHOT_SEQUENCE)
        _write_checkpoint_integer(io, UInt64(length(value)))
        for item in value
            _write_portable_value(io, item)
        end
    else
        _portable_snapshot_fail(
            :unsupported_value,
            "portable snapshot value type $(typeof(value)) has no registered stable encoding",
        )
    end
    return nothing
end

function _portable_value_bytes(value)
    io = IOBuffer()
    _write_portable_value(io, value)
    return take!(io)
end

function _read_portable_exact(io::IO, count::Integer, label::AbstractString)
    count >= 0 || _portable_snapshot_fail(:invalid_length, "$label length is negative")
    count <= typemax(Int) || _portable_snapshot_fail(
        :resource_limit,
        "$label length exceeds the addressable range",
    )
    requested = Int(count)
    bytes = read(io, requested)
    length(bytes) == requested || _portable_snapshot_fail(:truncated, "$label is truncated")
    return bytes
end

function _read_portable_string(io::IO, label::AbstractString; maximum_bytes::Integer)
    count = Int(_read_checkpoint_integer(io, UInt32))
    count <= maximum_bytes || _portable_snapshot_fail(:resource_limit, "$label exceeds its byte limit")
    bytes = _read_portable_exact(io, count, label)
    isvalid(String, bytes) || _portable_snapshot_fail(:invalid_text, "$label is not valid UTF-8")
    text = String(bytes)
    return _portable_snapshot_text(text, label)
end

function _read_portable_value(
    io::IO;
    maximum_payload_bytes::Integer,
    maximum_depth::Integer,
    depth::Integer = 0,
)
    depth <= maximum_depth || _portable_snapshot_fail(
        :resource_limit,
        "portable value nesting exceeds its limit",
    )
    eof(io) && _portable_snapshot_fail(:truncated, "portable value tag is missing")
    tag = read(io, UInt8)
    tag == _SNAPSHOT_NULL && return nothing
    tag == _SNAPSHOT_FALSE && return false
    tag == _SNAPSHOT_TRUE && return true
    tag == _SNAPSHOT_SIGNED && return reinterpret(Int64, _read_checkpoint_integer(io, UInt64))
    tag == _SNAPSHOT_UNSIGNED && return _read_checkpoint_integer(io, UInt64)
    if tag == _SNAPSHOT_FLOAT64
        value = reinterpret(Float64, _read_checkpoint_integer(io, UInt64))
        isfinite(value) || _portable_snapshot_fail(
            :nonfinite_value,
            "portable physical floating value is nonfinite",
        )
        return value
    end
    if tag == _SNAPSHOT_TEXT
        return _read_portable_string(io, "portable value text"; maximum_bytes = maximum_payload_bytes)
    elseif tag == _SNAPSHOT_BYTES
        count = _read_checkpoint_integer(io, UInt64)
        count <= UInt64(maximum_payload_bytes) ||
            _portable_snapshot_fail(:resource_limit, "portable byte value exceeds its limit")
        return _read_portable_exact(io, Int(count), "portable byte value")
    elseif tag == _SNAPSHOT_RATIONAL
        numerator_value = reinterpret(Int128, _read_checkpoint_integer(io, UInt128))
        denominator_value = _read_checkpoint_integer(io, UInt128)
        denominator_value <= UInt128(typemax(Int128)) ||
            _portable_snapshot_fail(:integer_overflow, "portable rational denominator exceeds Int128")
        denominator_value > 0 || _portable_snapshot_fail(:invalid_rational, "portable rational denominator is zero")
        return numerator_value // Int128(denominator_value)
    elseif tag == _SNAPSHOT_SEQUENCE
        count = _read_checkpoint_integer(io, UInt64)
        count <= UInt64(maximum_payload_bytes) ||
            _portable_snapshot_fail(:resource_limit, "portable sequence count exceeds its limit")
        return Any[
            _read_portable_value(
                io;
                maximum_payload_bytes,
                maximum_depth,
                depth = depth + 1,
            ) for _ in 1:Int(count)
        ]
    elseif tag == _SNAPSHOT_RECORD
        schema_id = _read_portable_string(io, "portable record schema"; maximum_bytes = maximum_payload_bytes)
        count = Int(_read_checkpoint_integer(io, UInt32))
        count <= maximum_payload_bytes ||
            _portable_snapshot_fail(:resource_limit, "portable record field count exceeds its limit")
        fields = Pair{String,Any}[]
        for _ in 1:count
            name = _read_portable_string(io, "portable record field"; maximum_bytes = maximum_payload_bytes)
            value = _read_portable_value(
                io;
                maximum_payload_bytes,
                maximum_depth,
                depth = depth + 1,
            )
            push!(fields, name => value)
        end
        record = PortableSnapshotRecord(schema_id, fields)
        first.(record.fields) == first.(fields) || _portable_snapshot_fail(
            :noncanonical_order,
            "portable record fields are not in canonical order",
        )
        return record
    elseif tag == _SNAPSHOT_ARRAY
        kind_text = _read_portable_string(io, "portable array element kind"; maximum_bytes = 32)
        element_kind = if kind_text == "float64"
            :float64
        elseif kind_text == "int64"
            :int64
        elseif kind_text == "uint64"
            :uint64
        else
            _portable_snapshot_fail(
                :unsupported_array_element,
                "portable array element kind is unsupported",
            )
        end
        rank = Int(_read_checkpoint_integer(io, UInt32))
        rank <= 32 || _portable_snapshot_fail(:resource_limit, "portable array rank exceeds 32")
        shape = Int[]
        for _ in 1:rank
            dimension = _read_checkpoint_integer(io, UInt64)
            dimension <= UInt64(typemax(Int)) ||
                _portable_snapshot_fail(:resource_limit, "portable array dimension exceeds Int")
            push!(shape, Int(dimension))
        end
        unit = _read_portable_string(io, "portable array unit"; maximum_bytes = 1024)
        axes = String[
            _read_portable_string(io, "portable array axis"; maximum_bytes = 1024)
            for _ in 1:rank
        ]
        byte_count = _read_checkpoint_integer(io, UInt64)
        byte_count <= UInt64(maximum_payload_bytes) ||
            _portable_snapshot_fail(:resource_limit, "portable array payload exceeds its limit")
        bytes = _read_portable_exact(io, Int(byte_count), "portable array payload")
        return PortableSnapshotArray(element_kind, shape, unit, axes, bytes)
    end
    _portable_snapshot_fail(:unknown_value_tag, "portable value tag $tag is unsupported")
end

function _metadata_record(
    metadata::PortableSnapshotMetadata;
    schema_minor::UInt16=_PORTABLE_SNAPSHOT_SCHEMA_MINOR,
)
    schema_minor <= _PORTABLE_SNAPSHOT_SCHEMA_MINOR || _portable_snapshot_fail(
        :unsupported_minor_version,
        "portable snapshot metadata minor version is unsupported",
    )
    fields = Pair{String,Any}[
        "accepted_step" => metadata.accepted_step,
        "capabilities" => metadata.capabilities,
        "creator_platform" => metadata.creator_platform,
        "model_signature_sha256" => metadata.model_signature_sha256,
        "profile" => String(metadata.profile),
        "project_signature_sha256" => metadata.project_signature_sha256,
        "provenance" => metadata.provenance,
        "represented_time_s" => metadata.represented_time_s,
        "settings_signature_sha256" => metadata.settings_signature_sha256,
        "topology_signature_sha256" => metadata.topology_signature_sha256,
        "writer_version" => metadata.writer_version,
    ]
    if schema_minor >= 1
        append!(fields, Pair{String,Any}[
            "compression" => metadata.compression,
            "minimum_reader_version" => metadata.minimum_reader_version,
            "numeric_profile" => metadata.numeric_profile,
        ])
    end
    return PortableSnapshotRecord(
        "aimora.snapshot.metadata.v1",
        fields,
    )
end

function _record_dictionary(record::PortableSnapshotRecord, expected_schema::AbstractString)
    record.schema_id == expected_schema || _portable_snapshot_fail(
        :unexpected_schema,
        "portable record schema $(record.schema_id) is not $expected_schema",
    )
    return Dict(record.fields)
end

function _metadata_from_record(record::PortableSnapshotRecord, schema_minor::UInt16)
    fields = _record_dictionary(record, "aimora.snapshot.metadata.v1")
    expected = Set(String[
        "accepted_step",
        "capabilities",
        "creator_platform",
        "model_signature_sha256",
        "profile",
        "project_signature_sha256",
        "provenance",
        "represented_time_s",
        "settings_signature_sha256",
        "topology_signature_sha256",
        "writer_version",
    ])
    schema_minor >= 1 && union!(expected, (
        "compression",
        "minimum_reader_version",
        "numeric_profile",
    ))
    Set(keys(fields)) == expected || _portable_snapshot_fail(
        :metadata_fields,
        "portable snapshot metadata fields are incomplete or unknown",
    )
    capabilities = fields["capabilities"]
    capabilities isa AbstractVector && all(value -> value isa AbstractString, capabilities) ||
        _portable_snapshot_fail(:metadata_type, "portable snapshot capabilities have the wrong type")
    represented_time = fields["represented_time_s"]
    represented_time isa Rational ||
        _portable_snapshot_fail(:metadata_type, "portable snapshot time has the wrong type")
    return PortableSnapshotMetadata(
        Symbol(fields["profile"]),
        fields["project_signature_sha256"],
        fields["model_signature_sha256"],
        fields["topology_signature_sha256"],
        fields["settings_signature_sha256"],
        represented_time,
        fields["accepted_step"],
        String[value for value in capabilities],
        fields["provenance"];
        writer_version = fields["writer_version"],
        creator_platform = fields["creator_platform"],
        numeric_profile = get(
            fields,
            "numeric_profile",
            "ieee754_binary64_finite_preserve_signed_zero",
        ),
        compression = get(fields, "compression", "none"),
        minimum_reader_version = get(fields, "minimum_reader_version", "1.0"),
    )
end

function _portable_snapshot_bytes(
    snapshot::PortableEMTSnapshot;
    schema_minor::UInt16=_PORTABLE_SNAPSHOT_SCHEMA_MINOR,
)
    schema_minor <= _PORTABLE_SNAPSHOT_SCHEMA_MINOR || _portable_snapshot_fail(
        :unsupported_minor_version,
        "portable snapshot writer minor version is unsupported",
    )
    body = IOBuffer()
    write(body, _PORTABLE_SNAPSHOT_MAGIC)
    _write_checkpoint_integer(body, _PORTABLE_SNAPSHOT_SCHEMA_MAJOR)
    _write_checkpoint_integer(body, schema_minor)
    metadata_bytes = _portable_value_bytes(_metadata_record(snapshot.metadata; schema_minor))
    _write_checkpoint_integer(body, UInt64(length(metadata_bytes)))
    write(body, sha256(metadata_bytes))
    write(body, metadata_bytes)
    _write_checkpoint_integer(body, UInt32(length(snapshot.sections)))
    for section in snapshot.sections
        _write_portable_string(body, section.identity)
        _write_checkpoint_integer(body, section.version_major)
        _write_checkpoint_integer(body, section.version_minor)
        write(body, section.visibility == :public ? UInt8(0) : UInt8(1))
        payload = _portable_value_bytes(section.value)
        _write_checkpoint_integer(body, UInt64(length(payload)))
        write(body, _portable_section_digest(
            section.identity,
            section.version_major,
            section.version_minor,
            section.visibility,
            payload,
        ))
        write(body, payload)
    end
    body_bytes = take!(body)
    return vcat(body_bytes, sha256(body_bytes))
end


portable_snapshot_bytes(snapshot::PortableEMTSnapshot) = _portable_snapshot_bytes(snapshot)

function _decode_portable_snapshot(
    bytes::AbstractVector{UInt8};
    decode_sections::Bool,
    allow_private::Bool,
    maximum_payload_bytes::Integer,
    maximum_sections::Integer,
    maximum_depth::Integer,
)
    maximum_payload_bytes >= 0 || _portable_snapshot_fail(
        :resource_limit,
        "portable payload limit must be nonnegative",
    )
    maximum_sections >= 0 || _portable_snapshot_fail(
        :resource_limit,
        "portable section limit must be nonnegative",
    )
    maximum_depth >= 0 || _portable_snapshot_fail(
        :resource_limit,
        "portable nesting limit must be nonnegative",
    )
    minimum_bytes = length(_PORTABLE_SNAPSHOT_MAGIC) + 2 * sizeof(UInt16) +
        sizeof(UInt64) + _PORTABLE_SNAPSHOT_DIGEST_BYTES + sizeof(UInt32) +
        _PORTABLE_SNAPSHOT_DIGEST_BYTES
    length(bytes) >= minimum_bytes || _portable_snapshot_fail(
        :truncated,
        "portable snapshot is shorter than its minimum envelope",
    )
    body_bytes = @view bytes[1:(end - _PORTABLE_SNAPSHOT_DIGEST_BYTES)]
    expected_digest = @view bytes[(end - _PORTABLE_SNAPSHOT_DIGEST_BYTES + 1):end]
    sha256(body_bytes) == expected_digest || _portable_snapshot_fail(
        :integrity,
        "portable snapshot envelope digest does not match",
    )
    io = IOBuffer(body_bytes)
    _read_portable_exact(io, length(_PORTABLE_SNAPSHOT_MAGIC), "portable snapshot magic") ==
        _PORTABLE_SNAPSHOT_MAGIC || _portable_snapshot_fail(
        :not_portable_snapshot,
        "file does not contain the portable EMT snapshot magic",
    )
    schema_major = _read_checkpoint_integer(io, UInt16)
    schema_minor = _read_checkpoint_integer(io, UInt16)
    schema_major == _PORTABLE_SNAPSHOT_SCHEMA_MAJOR || _portable_snapshot_fail(
        :unsupported_major_version,
        "portable snapshot major version $schema_major is unsupported",
    )
    schema_minor <= _PORTABLE_SNAPSHOT_SCHEMA_MINOR || _portable_snapshot_fail(
        :unsupported_minor_version,
        "portable snapshot minor version $schema_minor is unsupported",
    )
    metadata_length = _read_checkpoint_integer(io, UInt64)
    metadata_length <= UInt64(maximum_payload_bytes) ||
        _portable_snapshot_fail(:resource_limit, "portable metadata exceeds its byte limit")
    metadata_digest = _read_portable_exact(io, _PORTABLE_SNAPSHOT_DIGEST_BYTES, "metadata digest")
    metadata_bytes = _read_portable_exact(io, Int(metadata_length), "portable metadata")
    sha256(metadata_bytes) == metadata_digest ||
        _portable_snapshot_fail(:integrity, "portable metadata digest does not match")
    metadata_io = IOBuffer(metadata_bytes)
    metadata_value = _read_portable_value(
        metadata_io;
        maximum_payload_bytes,
        maximum_depth,
    )
    eof(metadata_io) || _portable_snapshot_fail(:trailing_bytes, "portable metadata has trailing bytes")
    metadata_value isa PortableSnapshotRecord ||
        _portable_snapshot_fail(:metadata_type, "portable metadata root is not a record")
    metadata = _metadata_from_record(metadata_value, schema_minor)
    section_count = Int(_read_checkpoint_integer(io, UInt32))
    section_count > 0 || _portable_snapshot_fail(
        :missing_section,
        "portable snapshot requires at least one section",
    )
    section_count <= maximum_sections ||
        _portable_snapshot_fail(:resource_limit, "portable section count exceeds its limit")
    descriptors = PortableSnapshotSectionDescriptor[]
    sections = PortableSnapshotSection[]
    previous_identity = ""
    for _ in 1:section_count
        identity = _portable_snapshot_identity(
            _read_portable_string(io, "portable section identity"; maximum_bytes = 1024),
            "portable section",
        )
        isempty(previous_identity) || previous_identity < identity ||
            _portable_snapshot_fail(:noncanonical_order, "portable sections are not strictly ordered")
        previous_identity = identity
        version_major = _read_checkpoint_integer(io, UInt16)
        version_minor = _read_checkpoint_integer(io, UInt16)
        eof(io) && _portable_snapshot_fail(:truncated, "portable section visibility is missing")
        visibility_code = read(io, UInt8)
        visibility = visibility_code == 0 ? :public : visibility_code == 1 ?
            :private_reconstructible : _portable_snapshot_fail(
                :invalid_section_visibility,
                "portable section visibility code is unsupported",
            )
        payload_length = _read_checkpoint_integer(io, UInt64)
        payload_length <= UInt64(maximum_payload_bytes) ||
            _portable_snapshot_fail(:resource_limit, "portable section payload exceeds its limit")
        payload_digest = _read_portable_exact(io, _PORTABLE_SNAPSHOT_DIGEST_BYTES, "section digest")
        payload = _read_portable_exact(io, Int(payload_length), "portable section payload")
        _portable_section_digest(
            identity,
            version_major,
            version_minor,
            visibility,
            payload,
        ) == payload_digest ||
            _portable_snapshot_fail(:integrity, "portable section $identity digest does not match")
        push!(descriptors, PortableSnapshotSectionDescriptor(
            identity,
            version_major,
            version_minor,
            visibility,
            payload_length,
            bytes2hex(payload_digest),
        ))
        if decode_sections
            visibility == :private_reconstructible && !allow_private &&
                _portable_snapshot_fail(
                    :private_section_unavailable,
                    "portable section $identity requires private reconstruction authority",
                )
            payload_io = IOBuffer(payload)
            value = _read_portable_value(
                payload_io;
                maximum_payload_bytes,
                maximum_depth,
            )
            eof(payload_io) || _portable_snapshot_fail(
                :trailing_bytes,
                "portable section $identity has trailing bytes",
            )
            push!(sections, PortableSnapshotSection(
                identity,
                version_major,
                version_minor,
                visibility,
                value,
            ))
        end
    end
    eof(io) || _portable_snapshot_fail(:trailing_bytes, "portable snapshot body has trailing bytes")
    descriptor = PortableSnapshotDescriptor(
        schema_major,
        schema_minor,
        metadata,
        descriptors,
        UInt64(length(bytes)),
        bytes2hex(expected_digest),
        schema_minor == 0 ? ["aimora.portable_emt.1.0_to_1.1"] : String[],
    )
    return decode_sections ? (PortableEMTSnapshot(metadata, sections), descriptor) : descriptor
end

function portable_snapshot_descriptor(snapshot::PortableEMTSnapshot)
    bytes = portable_snapshot_bytes(snapshot)
    return _decode_portable_snapshot(
        bytes;
        decode_sections = false,
        allow_private = true,
        maximum_payload_bytes = length(bytes),
        maximum_sections = max(length(snapshot.sections), 1),
        maximum_depth = 64,
    )
end

function _read_portable_snapshot_file(
    path::AbstractString;
    maximum_file_bytes::Integer,
)
    maximum_file_bytes > 0 || _portable_snapshot_fail(
        :resource_limit,
        "portable snapshot file limit must be positive",
    )
    input_path = abspath(String(path))
    isfile(input_path) || _portable_snapshot_fail(
        :missing_file,
        "portable snapshot file does not exist: $input_path",
    )
    file_bytes = filesize(input_path)
    file_bytes <= maximum_file_bytes || _portable_snapshot_fail(
        :resource_limit,
        "portable snapshot file exceeds its configured byte limit",
    )
    return read(input_path)
end

function inspect_portable_emt_snapshot(
    path::AbstractString;
    maximum_file_bytes::Integer = 2_000_000_000,
    maximum_payload_bytes::Integer = maximum_file_bytes,
    maximum_sections::Integer = 4096,
    maximum_depth::Integer = 64,
)
    bytes = _read_portable_snapshot_file(path; maximum_file_bytes)
    return _decode_portable_snapshot(
        bytes;
        decode_sections = false,
        allow_private = true,
        maximum_payload_bytes,
        maximum_sections,
        maximum_depth,
    )
end

function read_portable_emt_snapshot(
    path::AbstractString;
    allow_private::Bool = false,
    maximum_file_bytes::Integer = 2_000_000_000,
    maximum_payload_bytes::Integer = maximum_file_bytes,
    maximum_sections::Integer = 4096,
    maximum_depth::Integer = 64,
)
    snapshot, _ = read_portable_emt_snapshot_with_descriptor(
        path;
        allow_private,
        maximum_file_bytes,
        maximum_payload_bytes,
        maximum_sections,
        maximum_depth,
    )
    return snapshot
end

function read_portable_emt_snapshot_with_descriptor(
    path::AbstractString;
    allow_private::Bool = false,
    maximum_file_bytes::Integer = 2_000_000_000,
    maximum_payload_bytes::Integer = maximum_file_bytes,
    maximum_sections::Integer = 4096,
    maximum_depth::Integer = 64,
)
    bytes = _read_portable_snapshot_file(path; maximum_file_bytes)
    return _decode_portable_snapshot(
        bytes;
        decode_sections = true,
        allow_private,
        maximum_payload_bytes,
        maximum_sections,
        maximum_depth,
    )
end

function write_portable_emt_snapshot(
    path::AbstractString,
    snapshot::PortableEMTSnapshot,
)
    bytes = portable_snapshot_bytes(snapshot)
    output_path = abspath(String(path))
    mkpath(dirname(output_path))
    temporary_path, io = mktemp(dirname(output_path))
    try
        write(io, bytes)
        close(io)
        decoded, descriptor = _decode_portable_snapshot(
            read(temporary_path);
            decode_sections = true,
            allow_private = true,
            maximum_payload_bytes = length(bytes),
            maximum_sections = max(length(snapshot.sections), 1),
            maximum_depth = 64,
        )
        portable_snapshot_bytes(decoded) == bytes || _portable_snapshot_fail(
            :noncanonical_roundtrip,
            "portable snapshot did not reproduce canonical bytes before publication",
        )
        mv(temporary_path, output_path; force = true)
        return descriptor
    catch
        isopen(io) && close(io)
        isfile(temporary_path) && rm(temporary_path; force = true)
        rethrow()
    end
end

end

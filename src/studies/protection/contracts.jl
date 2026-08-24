@enum ProtectionMeasurementStage::UInt8 begin
    ProtectionInstantaneousStage = 0x01
    ProtectionSlidingRMSStage = 0x02
    ProtectionFundamentalPhasorStage = 0x03
    ProtectionSequencePhasorStage = 0x04
    ProtectionFrequencyStage = 0x05
    ProtectionROCOFStage = 0x06
end

@enum MagnitudeRelayDirection::UInt8 begin
    OverMagnitude = 0x01
    UnderMagnitude = 0x02
end

@enum MagnitudeRelayValueMode::UInt8 begin
    AbsoluteMagnitudeValue = 0x01
    SignedScalarValue = 0x02
end

@enum ProtectionTimerMode::UInt8 begin
    ProtectionTimerInstantaneous = 0x01
    ProtectionDefiniteTimer = 0x02
    ProtectionInverseTimer = 0x03
end

"""One typed, timestamped A200 value selected for causal protection use."""
struct ProtectionMeasurement
    measurement_id::Symbol
    channel::Symbol
    quantity::Symbol
    unit::String
    orientation::String
    stage::ProtectionMeasurementStage
    source_tick::Int
    release_tick::Int
    tick_s::Float64
    value::Union{Nothing,ComplexF64}
    quality::Symbol
    unavailable_reason::Union{Nothing,Symbol}
    measurement_signature_sha256::String
    deterministic_signature_sha256::String

    function ProtectionMeasurement(
        measurement_id::Symbol,
        channel::Symbol,
        quantity::Symbol,
        unit::AbstractString,
        orientation::AbstractString,
        stage::ProtectionMeasurementStage,
        source_tick::Integer,
        release_tick::Integer,
        tick_s::Real,
        value::Union{Nothing,Number},
        quality::Symbol,
        unavailable_reason::Union{Nothing,Symbol},
        measurement_signature_sha256::AbstractString,
    )
        isempty(String(measurement_id)) && throw(ArgumentError(
            "protection measurement id must not be empty",
        ))
        isempty(String(channel)) && throw(ArgumentError(
            "protection measurement channel must not be empty",
        ))
        isempty(String(quantity)) && throw(ArgumentError(
            "protection measurement quantity must not be empty",
        ))
        unit_string = String(unit)
        orientation_string = String(orientation)
        isempty(strip(unit_string)) && throw(ArgumentError(
            "protection measurement unit must not be empty",
        ))
        isempty(strip(orientation_string)) && throw(ArgumentError(
            "protection measurement orientation must not be empty",
        ))
        source = Int(source_tick)
        release = Int(release_tick)
        source >= 0 && release >= source || throw(ArgumentError(
            "protection measurement ticks must be nonnegative and causally ordered",
        ))
        tick = Float64(tick_s)
        isfinite(tick) && tick > 0.0 || throw(ArgumentError(
            "protection measurement tick must be finite and positive",
        ))
        isempty(String(quality)) && throw(ArgumentError(
            "protection measurement quality must not be empty",
        ))
        available = value !== nothing
        available == (unavailable_reason === nothing) || throw(ArgumentError(
            "protection measurement availability and refusal reason are inconsistent",
        ))
        complex_value = value === nothing ? nothing : ComplexF64(value)
        complex_value === nothing ||
            (isfinite(real(complex_value)) && isfinite(imag(complex_value))) ||
            throw(ArgumentError("protection measurement value must be finite"))
        measurement_signature = lowercase(String(measurement_signature_sha256))
        occursin(r"^[0-9a-f]{64}$", measurement_signature) || throw(ArgumentError(
            "protection measurement owner signature must be a 64-hex SHA-256",
        ))
        io = IOBuffer()
        print(
            io,
            measurement_id,
            '|',
            channel,
            '|',
            quantity,
            '|',
            unit_string,
            '|',
            orientation_string,
            '|',
            UInt8(stage),
            '|',
            source,
            '|',
            release,
            '|',
            repr(tick),
            '|',
            complex_value === nothing ? "unavailable" : repr(complex_value),
            '|',
            quality,
            '|',
            something(unavailable_reason, :available),
            '|',
            measurement_signature,
        )
        signature = bytes2hex(sha256(take!(io)))
        return new(
            measurement_id,
            channel,
            quantity,
            unit_string,
            orientation_string,
            stage,
            source,
            release,
            tick,
            complex_value,
            quality,
            unavailable_reason,
            measurement_signature,
            signature,
        )
    end
end

function _selected_protection_value(
    sample::MeasurementSample,
    channel_index::Int,
    stage::ProtectionMeasurementStage,
    sequence::Symbol,
)
    clipped = stage in (ProtectionSequencePhasorStage, ProtectionFrequencyStage) ?
        any(sample.clipped) : sample.clipped[channel_index]
    clipped && return nothing, :measurement_clipped
    if stage === ProtectionInstantaneousStage
        return ComplexF64(sample.instantaneous[channel_index]), nothing
    elseif stage === ProtectionSlidingRMSStage
        sample.sliding_rms === nothing && return nothing, :window_incomplete
        return ComplexF64(sample.sliding_rms[channel_index]), nothing
    elseif stage === ProtectionFundamentalPhasorStage
        sample.fundamental_rms_phasors === nothing &&
            return nothing, :phasor_unavailable
        return sample.fundamental_rms_phasors[channel_index], nothing
    elseif stage === ProtectionSequencePhasorStage
        sample.sequence_phasors === nothing && return nothing, :sequence_unavailable
        sequence in (:zero, :positive, :negative) || throw(ArgumentError(
            "protection sequence selection must be zero, positive, or negative",
        ))
        return getproperty(sample.sequence_phasors, sequence), nothing
    end
    stage === ProtectionROCOFStage && throw(ArgumentError(
        "ROCOF is protection-owned state and cannot be selected directly from A200",
    ))
    sample.frequency_hz === nothing && return nothing, :frequency_unavailable
    return ComplexF64(sample.frequency_hz), nothing
end

"""Select one declared channel/stage without rewriting the accepted measurement owner."""
function protection_measurement(
    sample::MeasurementSample,
    specification::MeasurementChainSpecification,
    channel::Symbol;
    stage::ProtectionMeasurementStage=ProtectionInstantaneousStage,
    sequence::Symbol=:positive,
)
    channel_index = findfirst(==(channel), specification.channel_names)
    channel_index === nothing && throw(ArgumentError(
        "protection channel $channel is absent from measurement $(specification.id)",
    ))
    value, unavailable_reason = _selected_protection_value(
        sample,
        channel_index,
        stage,
        sequence,
    )
    selected_channel = stage === ProtectionSequencePhasorStage ? sequence :
        (stage === ProtectionFrequencyStage ? :frequency : channel)
    return ProtectionMeasurement(
        specification.id,
        selected_channel,
        specification.quantity,
        stage === ProtectionFrequencyStage ? "Hz" : specification.unit,
        specification.orientation,
        stage,
        sample.source_tick,
        sample.release_tick,
        specification.acquisition.tick_s,
        value,
        sample.quality,
        unavailable_reason,
        sample.deterministic_signature_sha256,
    )
end

struct MagnitudeRelaySettings
    id::Symbol
    measurement_id::Symbol
    channel::Symbol
    stage::ProtectionMeasurementStage
    value_mode::MagnitudeRelayValueMode
    orientation_polarity::Int
    direction::MagnitudeRelayDirection
    timer_mode::ProtectionTimerMode
    pickup::Float64
    dropout_ratio::Float64
    definite_time_s::Float64
    inverse_a::Float64
    inverse_b::Float64
    inverse_p::Float64
    time_dial_s::Float64
    reset_time_s::Float64
    unit::String
    setting_provenance::ParameterProvenance
    timer_provenance::ParameterProvenance

    function MagnitudeRelaySettings(
        id::Symbol,
        measurement_id::Symbol,
        channel::Symbol;
        stage::ProtectionMeasurementStage=ProtectionSlidingRMSStage,
        value_mode::MagnitudeRelayValueMode=AbsoluteMagnitudeValue,
        orientation_polarity::Integer=1,
        direction::MagnitudeRelayDirection=OverMagnitude,
        timer_mode::ProtectionTimerMode=ProtectionTimerInstantaneous,
        pickup::Real,
        dropout_ratio::Real=0.95,
        definite_time_s::Real=0.0,
        inverse_a::Real=1.0,
        inverse_b::Real=0.0,
        inverse_p::Real=1.0,
        time_dial_s::Real=1.0,
        reset_time_s::Real=0.0,
        unit::AbstractString,
        setting_provenance::ParameterProvenance,
        timer_provenance::ParameterProvenance,
    )
        all(name -> !isempty(String(name)), (id, measurement_id, channel)) ||
            throw(ArgumentError("magnitude relay identities must not be empty"))
        values = Float64.((
            pickup,
            dropout_ratio,
            definite_time_s,
            inverse_a,
            inverse_b,
            inverse_p,
            time_dial_s,
            reset_time_s,
        ))
        all(isfinite, values) || throw(ArgumentError(
            "magnitude relay settings must be finite",
        ))
        pickup_value, dropout, definite_time, coefficient_a, coefficient_b,
            exponent, time_dial, reset_time = values
        polarity = Int(orientation_polarity)
        polarity in (-1, 1) || throw(ArgumentError(
            "magnitude relay orientation polarity must be -1 or 1",
        ))
        pickup_value > 0.0 || throw(ArgumentError(
            "magnitude relay pickup must be positive",
        ))
        0.0 < dropout <= 1.0 || throw(ArgumentError(
            "magnitude relay dropout ratio must be in (0, 1]",
        ))
        definite_time >= 0.0 && coefficient_a > 0.0 && coefficient_b >= 0.0 &&
            exponent > 0.0 && time_dial > 0.0 && reset_time >= 0.0 ||
            throw(ArgumentError("magnitude relay timer settings are outside their domain"))
        timer_mode === ProtectionDefiniteTimer && definite_time <= 0.0 &&
            throw(ArgumentError("definite-time magnitude relay requires positive delay"))
        unit_string = String(unit)
        isempty(strip(unit_string)) && throw(ArgumentError(
            "magnitude relay unit must not be empty",
        ))
        setting_provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "magnitude relay pickup provenance must be physical",
        ))
        timer_provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
            "magnitude relay timer provenance must be numerical policy",
        ))
        return new(
            id,
            measurement_id,
            channel,
            stage,
            value_mode,
            polarity,
            direction,
            timer_mode,
            pickup_value,
            dropout,
            definite_time,
            coefficient_a,
            coefficient_b,
            exponent,
            time_dial,
            reset_time,
            unit_string,
            setting_provenance,
            timer_provenance,
        )
    end
end

function magnitude_relay_signature(settings::MagnitudeRelaySettings)
    io = IOBuffer()
    for field in fieldnames(MagnitudeRelaySettings)
        value = getfield(settings, field)
        print(io, field, '=', repr(value), '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

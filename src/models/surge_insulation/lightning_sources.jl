export AbstractLightningImpulse,
       DoubleExponentialLightningImpulse,
       HeidlerLightningImpulse,
       LightningStroke,
       LightningStrokeSequence,
       LightningImpulseMetrics,
       lightning_current_a,
       lightning_impulse_metrics,
       lightning_sequence_current_a

abstract type AbstractLightningImpulse end

"""Normalized difference-of-exponentials current impulse with explicit polarity."""
struct DoubleExponentialLightningImpulse <: AbstractLightningImpulse
    peak_current_a::Float64
    slow_decay_rate_per_s::Float64
    fast_decay_rate_per_s::Float64
    peak_time_s::Float64
    normalization::Float64
    provenance::SurgeParameterProvenance

    function DoubleExponentialLightningImpulse(
        peak_current_a::Real,
        slow_decay_rate_per_s::Real,
        fast_decay_rate_per_s::Real;
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "ampere and inverse second",
            "analytic double-exponential lightning current with explicit peak and polarity",
        ),
    )
        peak = _finite_value(peak_current_a, "lightning peak current")
        peak != 0.0 || throw(ArgumentError("lightning peak current must be nonzero"))
        slow = _positive_finite(slow_decay_rate_per_s, "slow lightning decay rate")
        fast = _positive_finite(fast_decay_rate_per_s, "fast lightning decay rate")
        fast > slow || throw(ArgumentError(
            "fast lightning decay rate must exceed the slow decay rate",
        ))
        peak_time = log(fast / slow) / (fast - slow)
        normalization = exp(-slow * peak_time) - exp(-fast * peak_time)
        normalization > 0.0 || throw(ArgumentError(
            "double-exponential lightning normalization must be positive",
        ))
        _require_physical_provenance(provenance, "lightning impulse")
        return new(peak, slow, fast, peak_time, normalization, provenance)
    end
end

"""Numerically peak-normalized Heidler current impulse."""
struct HeidlerLightningImpulse <: AbstractLightningImpulse
    peak_current_a::Float64
    front_time_s::Float64
    tail_time_s::Float64
    exponent::Float64
    peak_time_s::Float64
    normalization::Float64
    provenance::SurgeParameterProvenance

    function HeidlerLightningImpulse(
        peak_current_a::Real,
        front_time_s::Real,
        tail_time_s::Real,
        exponent::Real;
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "ampere, second, and dimensionless exponent",
            "analytic Heidler lightning current with numerically verified peak",
        ),
    )
        peak = _finite_value(peak_current_a, "Heidler peak current")
        peak != 0.0 || throw(ArgumentError("Heidler peak current must be nonzero"))
        front = _positive_finite(front_time_s, "Heidler front time")
        tail = _positive_finite(tail_time_s, "Heidler tail time")
        order = _positive_finite(exponent, "Heidler exponent")
        tail > front || throw(ArgumentError("Heidler tail time must exceed front time"))
        raw(time_s) = begin
            ratio = time_s / front
            powered = ratio^order
            powered / (1.0 + powered) * exp(-time_s / tail)
        end
        lower = 0.0
        upper = order * tail
        for _ in 1:180
            middle = 0.5 * (lower + upper)
            stationarity = middle * (1.0 + (middle / front)^order) - order * tail
            if stationarity < 0.0
                lower = middle
            else
                upper = middle
            end
        end
        peak_time = 0.5 * (lower + upper)
        normalization = raw(peak_time)
        isfinite(normalization) && normalization > 0.0 || throw(ArgumentError(
            "Heidler lightning normalization must be finite and positive",
        ))
        _require_physical_provenance(provenance, "Heidler impulse")
        return new(peak, front, tail, order, peak_time, normalization, provenance)
    end
end

function lightning_current_a(impulse::DoubleExponentialLightningImpulse, time_s::Real)
    time = _finite_value(time_s, "lightning evaluation time")
    time < 0.0 && return 0.0
    return impulse.peak_current_a * (
        exp(-impulse.slow_decay_rate_per_s * time) -
        exp(-impulse.fast_decay_rate_per_s * time)
    ) / impulse.normalization
end

function lightning_current_a(impulse::HeidlerLightningImpulse, time_s::Real)
    time = _finite_value(time_s, "Heidler evaluation time")
    time < 0.0 && return 0.0
    ratio = time / impulse.front_time_s
    powered = ratio^impulse.exponent
    return impulse.peak_current_a * powered /
        (1.0 + powered) * exp(-time / impulse.tail_time_s) /
        impulse.normalization
end

struct LightningImpulseMetrics
    duration_s::Float64
    peak_current_a::Float64
    peak_time_s::Float64
    charge_c::Float64
    specific_energy_a2s::Float64
    sample_count::Int
end

function lightning_impulse_metrics(
    impulse::AbstractLightningImpulse;
    duration_s::Real,
    intervals::Integer=8192,
)
    duration = _positive_finite(duration_s, "lightning integration duration")
    interval_count = Int(intervals)
    interval_count >= 2 && iseven(interval_count) || throw(ArgumentError(
        "lightning integration intervals must be an even integer of at least two",
    ))
    step = duration / interval_count
    charge_sum = 0.0
    energy_sum = 0.0
    peak_current = 0.0
    peak_time = 0.0
    for index in 0:interval_count
        time = index * step
        current = lightning_current_a(impulse, time)
        weight = index == 0 || index == interval_count ? 1.0 : isodd(index) ? 4.0 : 2.0
        charge_sum += weight * current
        energy_sum += weight * current * current
        if abs(current) > abs(peak_current)
            peak_current = current
            peak_time = time
        end
    end
    return LightningImpulseMetrics(
        duration,
        peak_current,
        peak_time,
        charge_sum * step / 3.0,
        energy_sum * step / 3.0,
        interval_count + 1,
    )
end

struct LightningStroke{I<:AbstractLightningImpulse}
    identity::Symbol
    release_time_s::Float64
    impulse::I

    function LightningStroke(
        identity::Symbol,
        release_time_s::Real,
        impulse::I,
    ) where {I<:AbstractLightningImpulse}
        identity == Symbol("") && throw(ArgumentError("lightning stroke identity must not be empty"))
        release = _nonnegative_finite(release_time_s, "lightning stroke release time")
        return new{I}(identity, release, impulse)
    end
end

struct LightningStrokeSequence
    identity::Symbol
    strokes::Vector{LightningStroke}
    signature::String

    function LightningStrokeSequence(
        identity::Symbol,
        strokes::AbstractVector{<:LightningStroke},
    )
        identity == Symbol("") && throw(ArgumentError("lightning sequence identity must not be empty"))
        checked = LightningStroke[strokes...]
        isempty(checked) && throw(ArgumentError("lightning sequence must contain a stroke"))
        issorted(getfield.(checked, :release_time_s)) || throw(ArgumentError(
            "lightning stroke releases must be sorted",
        ))
        length(unique(getfield.(checked, :identity))) == length(checked) ||
            throw(ArgumentError("lightning stroke identities must be unique"))
        signature = bytes2hex(sha256(codeunits(join(
            ["aimora-lightning-sequence-v1", String(identity), repr(checked)],
            '\n',
        ))))
        return new(identity, checked, signature)
    end
end

function lightning_sequence_current_a(sequence::LightningStrokeSequence, time_s::Real)
    time = _finite_value(time_s, "lightning sequence evaluation time")
    current = 0.0
    for stroke in sequence.strokes
        stroke.release_time_s > time && break
        current += lightning_current_a(stroke.impulse, time - stroke.release_time_s)
    end
    return current
end

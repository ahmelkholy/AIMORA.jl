export InsulationStudyPlan,
       InsulationSampleResult,
       InsulationStudySummary,
       deterministic_insulation_margin,
       run_insulation_study

"""Preregistered Gaussian stress/strength sampling plan with stable indexed draws."""
struct InsulationStudyPlan
    identity::Symbol
    sample_count::Int
    seed::UInt64
    stress_mean_v::Float64
    stress_standard_deviation_v::Float64
    strength_mean_v::Float64
    strength_standard_deviation_v::Float64
    stress_strength_correlation::Float64
    confidence_level::Float64
    provenance::SurgeParameterProvenance

    function InsulationStudyPlan(
        identity::Symbol;
        sample_count::Integer,
        seed::Integer,
        stress_mean_v::Real,
        stress_standard_deviation_v::Real,
        strength_mean_v::Real,
        strength_standard_deviation_v::Real,
        stress_strength_correlation::Real=0.0,
        confidence_level::Real=0.95,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "volt and dimensionless probability",
            "synthetic preregistered correlated Gaussian insulation stress/strength study",
        ),
    )
        identity == Symbol("") && throw(ArgumentError("insulation study identity must not be empty"))
        count = Int(sample_count)
        count > 0 || throw(ArgumentError("insulation study sample count must be positive"))
        seed >= 0 || throw(ArgumentError("insulation study seed must be nonnegative"))
        checked_seed = UInt64(seed)
        stress_mean = _positive_finite(stress_mean_v, "mean insulation stress")
        stress_deviation = _nonnegative_finite(
            stress_standard_deviation_v,
            "insulation stress standard deviation",
        )
        strength_mean = _positive_finite(strength_mean_v, "mean insulation strength")
        strength_deviation = _nonnegative_finite(
            strength_standard_deviation_v,
            "insulation strength standard deviation",
        )
        correlation = _finite_value(
            stress_strength_correlation,
            "stress-strength correlation",
        )
        -1.0 <= correlation <= 1.0 || throw(ArgumentError(
            "stress-strength correlation must lie in [-1, 1]",
        ))
        confidence = _finite_value(confidence_level, "insulation confidence level")
        0.0 < confidence < 1.0 || throw(ArgumentError(
            "insulation confidence level must lie strictly between zero and one",
        ))
        _require_physical_provenance(provenance, "insulation study")
        return new(
            identity,
            count,
            checked_seed,
            stress_mean,
            stress_deviation,
            strength_mean,
            strength_deviation,
            correlation,
            confidence,
            provenance,
        )
    end
end

struct InsulationSampleResult
    index::Int
    first_standard_draw::Float64
    second_standard_draw::Float64
    stress_v::Float64
    strength_v::Float64
    margin_v::Float64
    failed::Bool
end

struct InsulationStudySummary
    identity::Symbol
    sample_count::Int
    failure_count::Int
    empirical_failure_probability::Float64
    confidence_lower::Float64
    confidence_upper::Float64
    minimum_margin_v::Float64
    maximum_stress_v::Float64
    minimum_strength_v::Float64
    samples::Vector{InsulationSampleResult}
    signature::String
end

deterministic_insulation_margin(stress_v::Real, strength_v::Real) =
    _positive_finite(strength_v, "deterministic insulation strength") -
    _nonnegative_finite(stress_v, "deterministic insulation stress")

function _standard_normal_pair(rng::AbstractRNG)
    first_uniform = max(rand(rng), nextfloat(0.0))
    second_uniform = rand(rng)
    radius = sqrt(-2.0 * log(first_uniform))
    angle = 2.0 * pi * second_uniform
    return radius * cos(angle), radius * sin(angle)
end

function _normal_quantile(probability::Float64)
    # Acklam's rational inverse-normal approximation; deterministic and adequate for confidence reporting.
    a = (-39.69683028665376, 220.9460984245205, -275.9285104469687,
         138.3577518672690, -30.66479806614716, 2.506628277459239)
    b = (-54.47609879822406, 161.5858368580409, -155.6989798598866,
         66.80131188771972, -13.28068155288572)
    c = (-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
         -2.549732539343734, 4.374664141464968, 2.938163982698783)
    d = (0.007784695709041462, 0.3224671290700398, 2.445134137142996,
         3.754408661907416)
    lower = 0.02425
    upper = 1.0 - lower
    if probability < lower
        q = sqrt(-2.0 * log(probability))
        return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
            ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0)
    elseif probability <= upper
        q = probability - 0.5
        r = q * q
        return (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
            (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1.0)
    end
    q = sqrt(-2.0 * log(1.0 - probability))
    return -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
        ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0)
end

function _wilson_interval(failures::Int, samples::Int, confidence_level::Float64)
    estimate = failures / samples
    z = _normal_quantile(0.5 + 0.5 * confidence_level)
    denominator = 1.0 + z^2 / samples
    centre = (estimate + z^2 / (2.0 * samples)) / denominator
    half_width = z / denominator * sqrt(
        estimate * (1.0 - estimate) / samples + z^2 / (4.0 * samples^2),
    )
    return max(centre - half_width, 0.0), min(centre + half_width, 1.0)
end

function run_insulation_study(plan::InsulationStudyPlan)
    rng = Random.Xoshiro(plan.seed)
    samples = Vector{InsulationSampleResult}(undef, plan.sample_count)
    correlation_complement = sqrt(max(1.0 - plan.stress_strength_correlation^2, 0.0))
    failure_count = 0
    minimum_margin = Inf
    maximum_stress = -Inf
    minimum_strength = Inf
    for index in 1:plan.sample_count
        first_draw, second_draw = _standard_normal_pair(rng)
        strength_draw = plan.stress_strength_correlation * first_draw +
            correlation_complement * second_draw
        stress = plan.stress_mean_v + plan.stress_standard_deviation_v * first_draw
        strength = plan.strength_mean_v +
            plan.strength_standard_deviation_v * strength_draw
        isfinite(stress) && isfinite(strength) || throw(ArgumentError(
            "insulation sampling produced a nonfinite value",
        ))
        margin = strength - stress
        failed = margin <= 0.0
        failure_count += failed
        minimum_margin = min(minimum_margin, margin)
        maximum_stress = max(maximum_stress, stress)
        minimum_strength = min(minimum_strength, strength)
        samples[index] = InsulationSampleResult(
            index,
            first_draw,
            second_draw,
            stress,
            strength,
            margin,
            failed,
        )
    end
    probability = failure_count / plan.sample_count
    lower, upper = _wilson_interval(
        failure_count,
        plan.sample_count,
        plan.confidence_level,
    )
    signature = bytes2hex(sha256(codeunits(join(
        [
            "aimora-insulation-study-v1",
            String(plan.identity),
            string(plan.sample_count),
            string(plan.seed),
            repr(samples),
        ],
        '\n',
    ))))
    return InsulationStudySummary(
        plan.identity,
        plan.sample_count,
        failure_count,
        probability,
        lower,
        upper,
        minimum_margin,
        maximum_stress,
        minimum_strength,
        samples,
        signature,
    )
end

export SurgeParameterProvenance,
       SurgeInsulationRefusal,
       SurgeInsulationDiagnostic,
       SurgeInsulationResult,
       surge_result_accepted,
       surge_physical_parameter_provenance

const SurgeParameterProvenance = ParameterProvenance

"""Typed refusal for an invalid, unavailable, or unsupported surge/insulation request."""
struct SurgeInsulationRefusal
    operation::Symbol
    reason::Symbol
    message::String
    accepted_state_unchanged::Bool

    function SurgeInsulationRefusal(
        operation::Symbol,
        reason::Symbol,
        message::AbstractString;
        accepted_state_unchanged::Bool=true,
    )
        operation == Symbol("") && throw(ArgumentError("surge refusal operation must not be empty"))
        reason == Symbol("") && throw(ArgumentError("surge refusal reason must not be empty"))
        text = strip(String(message))
        isempty(text) && throw(ArgumentError("surge refusal message must not be empty"))
        return new(operation, reason, text, accepted_state_unchanged)
    end
end

"""One quantity-specific physical or numerical diagnostic."""
struct SurgeInsulationDiagnostic
    quantity::Symbol
    value::Float64
    scale::Float64
    tolerance::Float64
    passed::Bool

    function SurgeInsulationDiagnostic(
        quantity::Symbol,
        value::Real,
        scale::Real,
        tolerance::Real,
    )
        quantity == Symbol("") && throw(ArgumentError("diagnostic quantity must not be empty"))
        checked_value = Float64(value)
        checked_scale = Float64(scale)
        checked_tolerance = Float64(tolerance)
        all(isfinite, (checked_value, checked_scale, checked_tolerance)) ||
            throw(ArgumentError("surge diagnostic values must be finite"))
        checked_scale > 0.0 || throw(ArgumentError("surge diagnostic scale must be positive"))
        checked_tolerance >= 0.0 || throw(ArgumentError("surge diagnostic tolerance must be nonnegative"))
        return new(
            quantity,
            checked_value,
            checked_scale,
            checked_tolerance,
            abs(checked_value) <= checked_tolerance * checked_scale,
        )
    end
end

"""Public aggregate result boundary for one deterministic or statistical surge study."""
struct SurgeInsulationResult{T}
    identity::Symbol
    status::Symbol
    payload::T
    diagnostics::Vector{SurgeInsulationDiagnostic}
    warnings::Vector{String}
    signature::String

    function SurgeInsulationResult(
        identity::Symbol,
        status::Symbol,
        payload::T,
        diagnostics::AbstractVector{SurgeInsulationDiagnostic}=SurgeInsulationDiagnostic[],
        warnings::AbstractVector{<:AbstractString}=String[],
    ) where {T}
        identity == Symbol("") && throw(ArgumentError("surge result identity must not be empty"))
        status in (:accepted, :failed, :unavailable, :unsupported) ||
            throw(ArgumentError("surge result status is not registered"))
        checked_diagnostics = collect(diagnostics)
        checked_warnings = String.(warnings)
        context = SHA.SHA2_256_CTX()
        SHA.update!(context, codeunits("aimora-surge-insulation-result-v1\n"))
        SHA.update!(context, codeunits(String(identity) * "\n" * String(status) * "\n"))
        SHA.update!(context, codeunits(repr(payload) * "\n"))
        for diagnostic in checked_diagnostics
            SHA.update!(context, codeunits(repr(diagnostic) * "\n"))
        end
        for warning in checked_warnings
            SHA.update!(context, codeunits(warning * "\n"))
        end
        signature = bytes2hex(SHA.digest!(context))
        return new{T}(
            identity,
            status,
            payload,
            checked_diagnostics,
            checked_warnings,
            signature,
        )
    end
end

surge_result_accepted(result::SurgeInsulationResult) =
    result.status === :accepted && all(diagnostic -> diagnostic.passed, result.diagnostics)

function surge_physical_parameter_provenance(
    source::AbstractString,
    units::AbstractString,
    transformation::AbstractString,
    uncertainty::AbstractString,
    validity_domain::AbstractString,
)
    return SurgeParameterProvenance(
        source,
        units,
        transformation,
        uncertainty,
        validity_domain,
        PhysicalModelParameter,
    )
end

function _surge_default_provenance(units::AbstractString, domain::AbstractString)
    return surge_physical_parameter_provenance(
        "caller-supplied generic surge/insulation parameters",
        units,
        "converted once to Float64 SI values without inferred bases or signs",
        "not supplied; caller retains parameter and model-form uncertainty ownership",
        domain,
    )
end

function _surge_nodes(positive_node::Integer, negative_node::Integer)
    positive = Int(positive_node)
    negative = Int(negative_node)
    positive >= 0 || throw(ArgumentError("positive surge terminal must be nonnegative"))
    negative >= 0 || throw(ArgumentError("negative surge terminal must be nonnegative"))
    positive != negative || throw(ArgumentError("surge terminals must be distinct"))
    return positive, negative
end

function _positive_finite(value::Real, name::AbstractString)
    checked = Float64(value)
    isfinite(checked) && checked > 0.0 ||
        throw(ArgumentError("$name must be finite and positive"))
    return checked
end

function _nonnegative_finite(value::Real, name::AbstractString)
    checked = Float64(value)
    isfinite(checked) && checked >= 0.0 ||
        throw(ArgumentError("$name must be finite and nonnegative"))
    return checked
end

function _finite_value(value::Real, name::AbstractString)
    checked = Float64(value)
    isfinite(checked) || throw(ArgumentError("$name must be finite"))
    return checked
end

function _require_physical_provenance(provenance::SurgeParameterProvenance, owner::AbstractString)
    provenance.nature === PhysicalModelParameter ||
        throw(ArgumentError("$owner provenance must describe physical model parameters"))
    return provenance
end

function _surge_companion_method(companion_method::Symbol, owner::AbstractString)
    companion_method in (:backward_euler, :BackwardEulerCompanion) && return :backward_euler
    companion_method in (:trapezoidal, :TrapezoidalCompanion) && return :trapezoidal
    throw(ArgumentError(
        "$owner supports trapezoidal execution and event-localization backward Euler only",
    ))
end

function _fill_two_terminal_stamp!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    branch_current_a::Float64,
    differential_conductance_s::Float64,
)
    length(terminal_current_a) >= 2 ||
        throw(DimensionMismatch("surge current workspace must contain two terminals"))
    size(terminal_jacobian_s, 1) >= 2 && size(terminal_jacobian_s, 2) >= 2 ||
        throw(DimensionMismatch("surge Jacobian workspace must be at least 2x2"))
    terminal_current_a[1] = branch_current_a
    terminal_current_a[2] = -branch_current_a
    terminal_jacobian_s[1, 1] = differential_conductance_s
    terminal_jacobian_s[1, 2] = -differential_conductance_s
    terminal_jacobian_s[2, 1] = -differential_conductance_s
    terminal_jacobian_s[2, 2] = differential_conductance_s
    return nothing
end

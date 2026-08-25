export ConverterSystemEventRecord,
       ConverterSystemState,
       ConverterSystemResult,
       converter_system_result

struct ConverterSystemEventRecord
    time_s::Float64
    kind::Symbol
    owner::Symbol
    accepted::Bool
    message::String

    function ConverterSystemEventRecord(time_s, kind::Symbol, owner::Symbol, accepted::Bool; message="")
        time = Float64(time_s)
        isfinite(time) && time >= 0.0 || throw(ArgumentError("converter event time must be finite and nonnegative"))
        isempty(String(kind)) && throw(ArgumentError("converter event kind must not be empty"))
        isempty(String(owner)) && throw(ArgumentError("converter event owner must not be empty"))
        return new(time, kind, owner, accepted, String(message))
    end
end

struct ConverterSystemState
    time_s::Float64
    port_voltage_v::Vector{Float64}
    port_current_a::Vector{Float64}
    requested_valve_state::BitVector
    applied_valve_state::BitVector
    conducting_valve_state::BitVector
    capacitor_charge_c::Vector{Float64}
    inductor_flux_wb::Vector{Float64}
    thermal_energy_j::Vector{Float64}
    control_state::Vector{Float64}
    stored_energy_j::Float64
    dissipated_energy_j::Float64
    accepted_step_count::Int
    event_count::Int
    deterministic_signature_sha256::String
end

struct ConverterSystemResult
    schema::Symbol
    specification_signature_sha256::String
    family::ConverterSystemFamily
    fidelity::ModelFidelity
    application::ConverterApplication
    accepted::Bool
    status::Symbol
    state::ConverterSystemState
    events::Vector{ConverterSystemEventRecord}
    maximum_kcl_residual_a::Float64
    relative_charge_residual::Float64
    relative_energy_residual::Float64
    harmonic_metrics::Dict{Symbol,Float64}
    warnings::Vector{String}
end

function converter_system_result(
    specification::ConverterSystemSpecification,
    state::ConverterSystemState;
    accepted::Bool,
    status::Symbol,
    events=ConverterSystemEventRecord[],
    maximum_kcl_residual_a::Real,
    relative_charge_residual::Real,
    relative_energy_residual::Real,
    harmonic_metrics=Dict{Symbol,Float64}(),
    warnings=String[],
)
    converter_system_is_ready(converter_system_readiness(specification)) || throw(ArgumentError(
        "a converter result cannot be built for a specification that is not ready",
    ))
    residuals = Float64.((maximum_kcl_residual_a, relative_charge_residual, relative_energy_residual))
    all(isfinite, residuals) && all(>=(0.0), residuals) || throw(ArgumentError(
        "converter residual metrics must be finite and nonnegative",
    ))
    return ConverterSystemResult(
        :aimora_converter_system_result_v1,
        specification.signature_sha256,
        specification.selection.family,
        specification.selection.fidelity,
        specification.selection.application,
        accepted,
        status,
        state,
        collect(events),
        residuals...,
        Dict{Symbol,Float64}(harmonic_metrics),
        String.(warnings),
    )
end

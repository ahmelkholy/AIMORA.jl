export SurgeInsulationProductFamily,
       InterruptionRestrikeProduct,
       ArresterProtectedTerminalProduct,
       TowerBackflashProduct,
       GISCoronaTerminalProduct,
       StatisticalInsulationProduct,
       SurgeInsulationProductSpecification,
       SurgeInsulationProductPreparation,
       SurgeInsulationProductReadiness,
       prepare_surge_insulation_product,
       surge_insulation_product_readiness

@enum SurgeInsulationProductFamily::UInt8 begin
    InterruptionRestrikeProduct = 0x01
    ArresterProtectedTerminalProduct = 0x02
    TowerBackflashProduct = 0x03
    GISCoronaTerminalProduct = 0x04
    StatisticalInsulationProduct = 0x05
end

const _SURGE_PRODUCT_UNSUPPORTED = (
    :field_resolved_electromagnetics,
    :grounding_safety_design,
    :insulation_design_or_lifetime,
    :protected_standard_conformance,
    :manufacturer_or_utility_equivalence,
    :field_or_laboratory_validation,
    :atp_or_pscad_equivalence,
    :safety_integrity,
    :certification,
)

function _validate_surge_product_components(
    family::SurgeInsulationProductFamily,
    components::NamedTuple,
)
    if family === InterruptionRestrikeProduct
        all(haskey(components, key) for key in (
            :breaker_specification,
            :arc_branches,
            :vacuum_gaps,
            :vacuum_branches,
        )) ||
            throw(ArgumentError("interruption product omits breaker, arc, or vacuum-gap owners"))
        components.breaker_specification isa EMTBreakerSpecification ||
            throw(ArgumentError("interruption breaker owner is invalid"))
        components.arc_branches isa Tuple && length(components.arc_branches) == 3 &&
            all(branch -> branch isa CombinedArcBranch, components.arc_branches) ||
            throw(ArgumentError("interruption product requires three combined arc branches"))
        components.vacuum_gaps isa Tuple && length(components.vacuum_gaps) == 3 &&
            all(gap -> gap isa VacuumInterruptionState, components.vacuum_gaps) ||
            throw(ArgumentError("interruption product requires three vacuum-gap states"))
        components.vacuum_branches isa Tuple && length(components.vacuum_branches) == 3 &&
            all(branch -> branch isa VacuumInterruptionBranch, components.vacuum_branches) &&
            all(index ->
                components.vacuum_branches[index].state === components.vacuum_gaps[index] &&
                components.vacuum_branches[index].arc_owner === components.arc_branches[index],
                1:3,
            ) || throw(ArgumentError(
            "interruption vacuum branches must bind the declared gaps and arcs",
        ))
    elseif family === ArresterProtectedTerminalProduct
        all(haskey(components, key) for key in (
            :apparatus_identity,
            :apparatus_preparation,
            :arrester,
            :dynamic_equivalent,
            :terminal_section,
        )) || throw(ArgumentError("arrester-terminal product omits an apparatus owner"))
        components.apparatus_identity isa Symbol &&
            components.apparatus_preparation isa TransformerApparatusPreparation &&
            components.arrester isa MetalOxideArrester &&
            components.dynamic_equivalent isa DynamicArresterEquivalent &&
            components.terminal_section isa SurgePropagationSection ||
            throw(ArgumentError("arrester-terminal product components have invalid types"))
    elseif family === TowerBackflashProduct
        all(haskey(components, key) for key in (
            :line_terminal_identity,
            :line_runtime_preparation,
            :lightning_sequence,
            :tower,
            :ground,
            :insulator,
        )) || throw(ArgumentError("tower/backflash product omits a surge-path owner"))
        components.line_terminal_identity isa Symbol &&
            components.line_runtime_preparation isa CoupledLineRuntimePreparation &&
            components.lightning_sequence isa LightningStrokeSequence &&
            components.tower isa TransmissionTowerModel &&
            components.ground isa IonizingGroundBranch &&
            components.insulator isa Union{DisruptiveEffectInsulator,LeaderProgressionInsulator} ||
            throw(ArgumentError("tower/backflash product components have invalid types"))
    elseif family === GISCoronaTerminalProduct
        all(haskey(components, key) for key in (
            :terminal_identity,
            :gis_section,
            :corona,
            :insulator,
        )) || throw(ArgumentError("GIS/corona product omits a terminal owner"))
        components.terminal_identity isa Symbol &&
            components.gis_section isa GISGILSection &&
            components.corona isa DynamicCoronaBranch &&
            components.insulator isa Union{DisruptiveEffectInsulator,LeaderProgressionInsulator} ||
            throw(ArgumentError("GIS/corona product components have invalid types"))
    else
        haskey(components, :study_plan) && components.study_plan isa InsulationStudyPlan ||
            throw(ArgumentError("statistical insulation product requires a preregistered study plan"))
    end
    return components
end

"""One original generic public product that binds complete surge owners without embedding a solver."""
struct SurgeInsulationProductSpecification{C<:NamedTuple}
    id::Symbol
    family::SurgeInsulationProductFamily
    components::C
    timestep_s::Float64
    stop_time_s::Float64
    provenance::SurgeParameterProvenance
    uncertainty::String
    validity_domain::String
    unsupported_phenomena::Tuple{Vararg{Symbol}}
    deterministic_signature_sha256::String

    function SurgeInsulationProductSpecification(
        id::Symbol,
        family::SurgeInsulationProductFamily,
        components::C;
        timestep_s::Real,
        stop_time_s::Real,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "product-specific SI surge quantities",
            "generic public surge/insulation product",
        ),
        uncertainty::AbstractString,
        validity_domain::AbstractString,
        unsupported_phenomena::AbstractVector{Symbol}=collect(_SURGE_PRODUCT_UNSUPPORTED),
    ) where {C<:NamedTuple}
        isempty(String(id)) && throw(ArgumentError("surge product identity must not be empty"))
        _validate_surge_product_components(family, components)
        timestep = _positive_finite(timestep_s, "surge product timestep")
        stop_time = _positive_finite(stop_time_s, "surge product stop time")
        stop_time >= timestep || throw(ArgumentError(
            "surge product stop time must contain at least one timestep",
        ))
        _require_physical_provenance(provenance, "surge product")
        uncertainty_text = strip(String(uncertainty))
        validity_text = strip(String(validity_domain))
        isempty(uncertainty_text) && throw(ArgumentError("surge product uncertainty must be explicit"))
        isempty(validity_text) && throw(ArgumentError("surge product validity must be explicit"))
        unsupported = Tuple(Symbol.(unsupported_phenomena))
        allunique(unsupported) && all(in(unsupported), _SURGE_PRODUCT_UNSUPPORTED) ||
            throw(ArgumentError("surge product must retain every global unsupported boundary"))
        io = IOBuffer()
        println(io, "aimora-surge-insulation-product-v1")
        println(io, id, '|', UInt8(family), '|', bitstring(timestep), '|', bitstring(stop_time))
        println(io, repr(components))
        println(io, repr(provenance))
        println(io, uncertainty_text)
        println(io, validity_text)
        println(io, join(String.(unsupported), ','))
        return new{C}(
            id,
            family,
            components,
            timestep,
            stop_time,
            provenance,
            uncertainty_text,
            validity_text,
            unsupported,
            bytes2hex(sha256(take!(io))),
        )
    end
end

struct SurgeInsulationProductPreparation{S,C}
    specification::S
    initial_time_s::Float64
    components::C
    preparation_signature_sha256::String
end

struct SurgeInsulationProductReadiness
    ready::Bool
    code::Symbol
    product::Symbol
    family::SurgeInsulationProductFamily
    production_backend_available::Bool
    unsupported_phenomena::Tuple{Vararg{Symbol}}
    deterministic_signature_sha256::String
end

function prepare_surge_insulation_product(
    specification::SurgeInsulationProductSpecification;
    initial_time_s::Real=0.0,
)
    initial_time = _nonnegative_finite(initial_time_s, "surge product initial time")
    components = deepcopy(specification.components)
    io = IOBuffer()
    println(io, specification.deterministic_signature_sha256)
    println(io, bitstring(initial_time))
    println(io, repr(components))
    return SurgeInsulationProductPreparation(
        specification,
        initial_time,
        components,
        bytes2hex(sha256(take!(io))),
    )
end

function surge_insulation_product_readiness(
    preparation::SurgeInsulationProductPreparation;
    production_backend_available::Bool=false,
)
    specification = preparation.specification
    return SurgeInsulationProductReadiness(
        true,
        production_backend_available ? :ready_for_coupled_execution :
            :ready_for_solver_free_inspection,
        specification.id,
        specification.family,
        production_backend_available,
        specification.unsupported_phenomena,
        preparation.preparation_signature_sha256,
    )
end

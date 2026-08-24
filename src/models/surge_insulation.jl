module SurgeInsulation

using LinearAlgebra
using Random
using SHA
using Statistics

using ..NonlinearNetwork
using ..CoupledLineRuntime: CoupledLineRuntimePreparation
using ..ProtectionStudy: EMTBreakerSpecification
using ..StudyCore: ParameterProvenance,
                   PhysicalModelParameter,
                   ScalingBasisParameter,
                   NumericalPolicyParameter
using ..TransformerApparatus: TransformerApparatusPreparation

include(joinpath(@__DIR__, "surge_insulation", "types.jl"))
include(joinpath(@__DIR__, "surge_insulation", "lightning_sources.jl"))
include(joinpath(@__DIR__, "surge_insulation", "arc_models.jl"))
include(joinpath(@__DIR__, "surge_insulation", "arrester_models.jl"))
include(joinpath(@__DIR__, "surge_insulation", "grounding_and_corona.jl"))
include(joinpath(@__DIR__, "surge_insulation", "insulation_models.jl"))
include(joinpath(@__DIR__, "surge_insulation", "traveling_wave_models.jl"))
include(joinpath(@__DIR__, "surge_insulation", "insulation_studies.jl"))
include(joinpath(@__DIR__, "surge_insulation", "public_products.jl"))
include(joinpath(@__DIR__, "surge_insulation", "portable_state.jl"))

end

module ConverterSystems

using SHA
using ..BridgeTopologies: BridgeTopologyDescriptor, bridge_topology_signature
using ..SwitchDetailedVSC
using ..StudyCore: AverageValue,
                   ModelFidelity,
                   ModelValidityDomain,
                   ParameterProvenance,
                   PhysicalModelParameter,
                   SwitchingDetailed,
                   SwitchingStateEquivalent

include(joinpath(@__DIR__, "converter_systems", "contracts.jl"))
include(joinpath(@__DIR__, "converter_systems", "selection.jl"))
include(joinpath(@__DIR__, "converter_systems", "equations.jl"))
include(joinpath(@__DIR__, "converter_systems", "modulation.jl"))
include(joinpath(@__DIR__, "converter_systems", "results.jl"))
include(joinpath(@__DIR__, "converter_systems", "average_buck_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "average_boost_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "average_inverting_buck_boost_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "switching_chopper_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "four_quadrant_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "three_phase_two_level_inverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "three_level_neutral_point_clamped_inverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "three_level_t_type_inverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "flying_capacitor_inverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "cascaded_h_bridge_inverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "matrix_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "cycloconverter.jl"))
include(joinpath(@__DIR__, "converter_systems", "interleaved_chopper_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "dual_active_bridge_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "average_dual_active_bridge_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "line_commutated_rectifier_converter.jl"))
include(joinpath(@__DIR__, "converter_systems", "converter_applications.jl"))

end

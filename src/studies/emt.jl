module EMTStudy
using Printf
using SHA
using Serialization
using Unicode
using ..PortableSnapshots
using ..Sources: AnalyticSourceSignal

include(joinpath(@__DIR__, "emt", "source_signal_program.jl"))
include(joinpath(@__DIR__, "emt", "types_and_horizons.jl"))
include(joinpath(@__DIR__, "emt", "nonlinear_network_execution.jl"))
include(joinpath(@__DIR__, "emt", "control_and_machine_boundaries.jl"))
include(joinpath(@__DIR__, "emt", "direct_machine_fleet.jl"))
include(joinpath(@__DIR__, "emt", "deck_context_initialization.jl"))
include(joinpath(@__DIR__, "emt", "fixed_source_load_flow.jl"))
include(joinpath(@__DIR__, "emt", "fault_equivalent_fitting.jl"))
include(joinpath(@__DIR__, "emt", "line_and_machine_context.jl"))
include(joinpath(@__DIR__, "emt", "dc_simulator_sources.jl"))
include(joinpath(@__DIR__, "emt", "coupled_branch_context.jl"))
include(joinpath(@__DIR__, "emt", "switching_resistor_context.jl"))
include(joinpath(@__DIR__, "emt", "triggered_timed_resistance_context.jl"))
include(joinpath(@__DIR__, "emt", "piecewise_nonlinear_inductor_context.jl"))
include(joinpath(@__DIR__, "emt", "hysteretic_inductor_context.jl"))
include(joinpath(@__DIR__, "emt", "nonlinear_context_and_deck_trace.jl"))
include(joinpath(@__DIR__, "emt", "hybrid_execution.jl"))
include(joinpath(@__DIR__, "emt", "three_phase_vsc.jl"))
include(joinpath(@__DIR__, "emt", "extended_vsc_platform.jl"))
include(joinpath(@__DIR__, "emt", "terminal_state.jl"))
include(joinpath(@__DIR__, "emt", "trace_and_control_inputs.jl"))
include(joinpath(@__DIR__, "emt", "control_expression_execution.jl"))
include(joinpath(@__DIR__, "emt", "control_frequency_initialization.jl"))
include(joinpath(@__DIR__, "emt", "control_device_execution.jl"))
include(joinpath(@__DIR__, "emt", "ensemble_execution.jl"))
include(joinpath(@__DIR__, "emt", "checkpoint_io.jl"))
include(joinpath(@__DIR__, "emt", "portable_workspace_snapshot.jl"))
include(joinpath(@__DIR__, "emt", "portable_nodal_element_snapshot.jl"))
include(joinpath(@__DIR__, "emt", "portable_control_snapshot.jl"))
include(joinpath(@__DIR__, "emt", "portable_task_scheduler_snapshot.jl"))
include(joinpath(@__DIR__, "emt", "portable_hybrid_snapshot.jl"))
include(joinpath(@__DIR__, "emt", "restart_execution.jl"))
include(joinpath(@__DIR__, "emt", "deck_timestep_orchestration.jl"))
include(joinpath(@__DIR__, "emt", "case_sequence_execution.jl"))
include(joinpath(@__DIR__, "emt", "consistent_initialization.jl"))
include(joinpath(@__DIR__, "emt", "snapshot_results_and_json.jl"))
include(joinpath(@__DIR__, "emt", "result_serialization.jl"))

end

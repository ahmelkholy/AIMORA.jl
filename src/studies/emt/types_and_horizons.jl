
using LinearAlgebra
using Printf
using TOML

using ..Branches: BranchCompanionSnapshot,
                  BreqivHistoryInjection,
                  GeneratorEquivalentModalBranch,
                  IdealTransformerVoltageConstraint,
                  TheveninSource,
                  CurrentInjection,
                  ConductanceBranch,
                  SeriesRLBranch,
                  SeriesRLCBranch,
                  CoupledInductiveBranch,
                  CoupledSeriesRLBranch,
                  CapacitorBranch,
                  advance_breqiv_history_current!,
                  branch_companion_snapshot,
                  branch_current_value,
                  companion,
                  generator_equivalent_history_injection,
                  initialize_breqiv_history_injection!,
                  seed_breqiv_frequency_histories!,
                  stamp_history_current!,
                  trace_output_channel_count,
                  trace_output_channel_names!,
                  trace_output_is_public,
                  trace_output_values!,
                  three_phase_breqiv_history_injection,
                  update!
using ..DeckParser
using ..Inverter
using ..Lines: distributed_transposed_line_modal_timestep_update!,
               BergeronLine,
               ComplexModalBergeronLine,
               CoupledFrequencyDependentLine,
               DistributedTransposedLineHistoryState,
               NestedCableTransientLineState,
               SampledFrequencyDependentLine,
               SampledFrequencyDependentLineGroup,
               SemlyenFrequencyDependentLine,
               distributed_transposed_line_steady_state_pi_equivalent,
               distributed_transposed_line_history_state,
               distributed_transposed_line_initial_history_phase_matrix,
               distributed_transposed_line_initial_history_phase_terms,
               distributed_transposed_line_initial_history_modal_components,
               distributed_transposed_line_initial_history_phasors,
               distributed_transposed_line_initial_history_seed,
               distributed_transposed_line_history_current_injection!,
               distributed_transposed_line_norton_coefficients,
               distributed_transposed_line_phase_current_injection!,
               frequency_dependent_line_recursive_convolution_update!,
               initialize_complex_modal_bergeron_steady_state!,
               initialize_semlyen_line_steady_state!,
               cascaded_phase_pi_equivalent,
               complex_modal_bergeron_steady_state_terminal_admittance,
               sampled_line_steady_state_terminal_admittance,
               semlyen_line_steady_state_terminal_admittance,
               LineStepResponseExponentialFitResult,
               line_step_response_exponential_fit,
               pole_residue_transfer_value
using ..Nodal
using ..NonlinearNetwork: NonlinearChatterDecision,
                          NonlinearNetworkScales,
                          NonlinearSolveDiagnostics
using ..NonlinearNodal: NonlinearNodalSystem,
                        advance_nonlinear_step!,
                        nonlinear_nodal_checkpoint,
                        nonlinear_linear_system,
                        restore_nonlinear_nodal_checkpoint!,
                        solve_nonlinear_algebraic_state!
using ..Sources: AnalyticSourceSignal,
                 SinusoidalSourceSignal,
                 sinusoidal_source_peak_phasor
import ..EMTTaskPlatform
using ..Machines: MachineNetworkCouplingState,
                  MachineTerminalCurrentState,
                  SynchronousMachineDynamicState,
                  SynchronousMachineTACSInterfaceState,
                  SynchronousMachineRotorMass,
                  SynchronousMachineParameters,
                  SynchronousMachineInitialization,
                  InductionMachineParameters,
                  InductionMachineState,
                  MachineTerminalPredictionState,
                  CoupledDQMachineParameters,
                  CoupledDQMachineState,
                  DirectCurrentMachineParameters,
                  DirectCurrentMachineState,
                  UniversalMachineType4Parameters,
                  UniversalMachineType4State,
                  induction_machine_steady_state_equivalent,
                  induction_machine_initial_state,
                  induction_machine_axis_fluxes,
                  machine_network_coupling_state_update!,
                  machine_terminal_current_state_update!,
                  synchronous_machine_initial_state,
                  synchronous_machine_dynamic_step!,
                  synchronous_machine_tacs_transfer_update!,
                  induction_machine_step!,
                  coupled_dq_machine_step!,
                  coupled_dq_power_transform,
                  predict_machine_terminal_currents!,
                  universal_machine_type4_step!,
                  single_phase_induction_steady_state_initialization
using ..Nonlinear: SaturatedTransformerNonlinearArrays,
                   SaturatedTransformerNonlinearSlopeBranch,
                   SaturableInductorBranch,
                   PowerSemiconductorSwitch,
                   PowerSemiconductorBridgeLeg,
                   PowerSemiconductorBridgeTopology,
                   power_semiconductor_event_localization!,
                   power_semiconductor_gate_transition_time,
                   power_semiconductor_bridge_gate_transition_time,
                   power_semiconductor_forward_turn_on_residual,
                   power_semiconductor_forward_extinction_residual,
                   power_semiconductor_reverse_turn_on_residual,
                   power_semiconductor_reverse_extinction_residual,
                   apply_power_semiconductor_gate_transition!,
                   apply_power_semiconductor_forward_turn_on!,
                   apply_power_semiconductor_forward_extinction!,
                   apply_power_semiconductor_reverse_turn_on!,
                   apply_power_semiconductor_reverse_extinction!,
                   apply_power_semiconductor_bridge_gate_transitions!,
                   request_power_semiconductor_gate!,
                   request_power_semiconductor_bridge_pole!,
                   request_power_semiconductor_topology_gates!,
                   power_semiconductor_bridge_switch,
                   power_semiconductor_bridge_topology_valves,
                   HysteresisLoopPreprocessResult,
                   SATURATED_TRANSFORMER_NONLINEAR_TYPE,
                   HYSTERETIC_INDUCTOR_NONLINEAR_TYPE,
                   SWITCHING_NONLINEAR_RESISTOR_TYPE,
                   TRIGGERED_TIMED_RESISTANCE_TYPE,
                   PIECEWISE_NONLINEAR_INDUCTOR_TYPE,
                   PSEUDO_NONLINEAR_INDUCTOR_TYPE,
                   _is_pseudo_nonlinear_inductor_type,
                   set_saturated_transformer_nonlinear_slope!,
                   saturated_transformer_nonlinear_arrays,
                   saturated_transformer_nonlinear_slope_branches,
                   saturated_transformer_steady_state_branch_admittance,
                   saturated_transformer_winding_branch_assembly,
                   saturated_transformer_winding_current_config
using ..OVER16TimestepIntegration: OVER16AcceptedTimestepState,
                                  HybridEventDirection,
                                  HYBRID_EVENT_FALLING,
                                  HYBRID_EVENT_ANY,
                                  HYBRID_EVENT_RISING,
                                  HybridEventPolicy,
                                  AbstractHybridEventSurface,
                                  HybridEventSurface,
                                  HybridEventOccurrence,
                                  SampledTaskOccurrence,
                                  hybrid_event_value,
                                  hybrid_event_candidate_time,
                                  hybrid_event_candidate_is_event,
                                  hybrid_event_bracket,
                                  localize_hybrid_event!,
                                  apply_hybrid_event!,
                                  AbstractExactSampledTask,
                                  ExactSampledTask,
                                  ExactSampledControlTask,
                                  ExactPWMTask,
                                  ExactSampledTaskScheduler,
                                  sampled_task_next_tick,
                                  sampled_task_power_history_invalidating,
                                  sampled_task_checkpoint,
                                  restore_sampled_task_checkpoint!,
                                  sampled_task_scheduler_last_run_invalidated_power,
                                  next_sampled_task_time,
                                  run_due_sampled_tasks!,
                                  GeneralEMTTask,
                                  GeneralEMTTaskScheduler,
                                  GeneralEMTTaskSchedulerCheckpoint,
                                  ExactSampledTaskCompatibilityAdapter,
                                  adapt_exact_sampled_task_scheduler,
                                  general_task_scheduler_checkpoint,
                                  general_task_scheduler_result,
                                  general_task_scheduler_last_run_invalidated_power,
                                  general_task_scheduler_occurrence_count,
                                  next_general_task_instant,
                                  restore_general_task_scheduler_checkpoint!,
                                  run_due_general_tasks!,
                                  TimestepTransaction,
                                  StableStructureTimestepTransaction,
                                  TimestepStateRestorer,
                                  begin_timestep_transaction!,
                                  restore_timestep_transaction!,
                                  commit_timestep_transaction!,
                                  timestep_transaction_active,
                                  timestep_transaction_status,
                                  restore_timestep_state!,
                                  _over16_sparse_switch_state_flow_update_lean!,
                                  over16_accepted_timestep_update!,
                                  over16_sparse_switch_state_flow_update!
using ..Companion: past_machine_history_update!,
                   series_rlc_companion
using ..Switches: OVER16FortranSparseFactorWorkspaceState,
                 CurrentZeroSwitch,
                 configure_current_extinction!,
                 current_extinction_enabled,
                 IdealSwitch,
                 TimeSwitch,
                 OVER16SwitchAdmittanceState,
                 OVER16SwitchAlterationState,
                 OVER16SwitchBValueExportState,
                 OVER16SwitchCurrentState,
                 OVER16SwitchOperationState,
                 OVER16SwitchPostCurrentState,
                 OVER16SwitchRetriangularizationState,
                 OVER16SwitchScanState,
                 OVER16SwitchSparseFactorWorkspaceState,
                 OVER16SwitchTopologyState,
                 prepare_current_zero_switch!,
                 apply_current_zero_transition!,
                 switch_closed,
                 switch_conductance,
                 over16_network_solution_update!,
                 over16_switch_simple_ordering
using ..TACS: AlgebraicControlAssignment,
              ControlledSwitchDelayedArcState,
              apply_controlled_switch_delayed_arc_transition!,
              ControlTransferFunction,
              ConstantControlSignal,
              SinusoidalControlSignal,
              ControlExpressionRuntime,
              ControlSystemExecutionState,
              ControlSystemExecutionTraceState,
              ControlSystemDevice,
              SignedControlSignalTerm,
              TACSControlledSwitch,
              advance_control_system_state!,
              compile_control_expression,
              control_system_output_values,
              control_system_execution_trace_update!,
              control_system_step!,
              controlled_switch_closed,
              sync_controlled_switch!,
              evaluate_control_expression!,
              initialize_control_function_steady_state!,
              initialize_control_system_utilities!,
              _CONTROL_SYSTEM_UTILITY_VALUES
using ..TACS: sinusoidal_control_signal_phasor,
              sinusoidal_control_signal_value
using ..Sources: analytic_source_value
using ..Timestep: OVER16CSUPDeviceInputTerm,
                  OVER16CSUPDeviceRow,
                  OVER16CSUPState,
                  ControlOutputReferenceResolutionState,
                  OVER16OutputReportState,
                  OVER16OutputWriterControlState,
                  OVER16PostExtremaControlState,
                  OVER16ELECTACSTOOutputNode,
                  OVER16SourceCardState,
                  OVER16SourceUpdateState,
                  OVER16TACSUtilityState,
                  UniversalMachineRuntimeState,
                  UniversalMachineRuntimeStorageLayout,
                  control_output_reference_resolution_state,
                  control_output_reference_resolution_update!,
                  over16_csup_device_input_term,
                  over16_csup_device_row,
                  over16_csup_device_update!,
                  over16_csup_device_update_from_legacy_descriptors!,
                  over16_csup_update_from_legacy_descriptors!,
                  over16_elec_tacsto_xref1_propagation_request,
                  over16_elec_tacsto_xref1_reference_entry,
                  over16_elec_tacsto_xref1_reference_request,
                  over16_elec_tacsto_xref1_table_entry,
                  over16_source_row_update!,
                  over16_source_row_update_preview,
                  over16_elec_tacsto_error_state,
                  over16_elec_tacsto_output_node,
                  over16_elec_tacsto_output_update!,
                  over16_elec_tacsto_report_state,
                  over16_output_report_update!,
                  over16_synchronous_machine_output_preview,
                  over16_source_card_update_intent_preview,
                  over16_source_card_update!,
                  over16_source_update!,
                  universal_machine_runtime_state_update!,
                  universal_machine_runtime_storage_layout
using ..StudyCore
using ..LineConstantsStudy:
    LineConstantsStudyResult,
    run_line_constants_study
using ..ValidationCore: assert_valid!, validation_result

export UnifiedEMTConfig,
       NonlinearEMTDiscontinuity,
       NonlinearEMTStudySchedule,
       NonlinearEMTStudyTrace,
       evaluate_nonlinear_emt_network!,
       EMTStepContext,
       EMTTerminalNodeState,
       EMTTerminalBranchState,
       EMTTerminalNonlinearState,
       EMTTerminalSwitchState,
       EMTTerminalState,
       SeriesRLCAlteration,
       SeriesRLCAlterationRecord,
       DeckEMTExecution,
       PreparedEMTStudy,
       PreparedMachineEMTStudy,
       AverageValueGridFollowingConverterInitialization,
       PreparedAverageValueGridFollowingEMTStudy,
       AbstractEMTHarmonicFormulation,
       PhysicalFrequencyFormulation,
       TimestepMatchedFormulation,
       EMTInitializationTolerances,
       OperatingPointQuantity,
       EMTOperatingPoint,
       EMTInitializationRequest,
       EMTInitializationResidual,
       OperatingPointMappingRecord,
       EMTInitializationTopologyReport,
       EMTInitializationFrequencyPoint,
       NoArtificialTransientMetric,
       EMTInitializationStateRecord,
       EMTInitializationFailure,
       EMTInitializationReport,
       EMTInitializationResult,
       DeckCaseInitialization,
       initialize_emt_study,
       initialization_accepted,
       EMTStudyWorkspace,
       EMTStepTransaction,
       EMTHybridEventSurface,
       EMTHybridEventPolicy,
       EMTHybridEventOccurrence,
       EMTSampledTaskOccurrence,
       EMTExactSampledTask,
       EMTExactSampledControlTask,
       EMTExactPWMTask,
       EMTExactSampledTaskScheduler,
       next_emt_sampled_task_time,
       run_due_emt_sampled_tasks!,
       GeneralEMTTask,
       GeneralEMTTaskScheduler,
       GeneralEMTTaskSchedulerCheckpoint,
       ExactSampledTaskCompatibilityAdapter,
       adapt_exact_sampled_task_scheduler,
       general_task_scheduler_checkpoint,
       general_task_scheduler_result,
       next_general_task_instant,
       restore_general_task_scheduler_checkpoint!,
       run_due_general_tasks!,
       EMTHybridStepIntegrator,
       PowerSemiconductorGateCommand,
       PowerSemiconductorBridgePoleCommand,
       PowerSemiconductorTopologyGateCommand,
       configure_emt_hybrid_execution,
       advance_emt_hybrid_step!,
       evaluate_emt_hybrid_study!,
       emt_hybrid_execution_status,
       EMTRestartMutationRecord,
       EMTRestartRun,
       EMTStudyBatch,
       AbstractEMTStudyCandidate,
       AbstractEMTModelParameter,
       EMTSourceCrestCandidate,
       EMTSourceFrequencyRateParameter,
       EMTModelParameterCandidate,
       AbstractEMTExecutionBackend,
       EMTCPUBackend,
       AbstractEMTTraceReducer,
       EMTOutputRMSReducer,
       EMTOutputPeakReducer,
       prepare_emt_study,
       reset_emt_study!,
       begin_emt_step_transaction!,
       provisional_emt_step!,
       restore_emt_step_transaction!,
       commit_emt_step_transaction!,
       emt_step_transaction_status,
       emt_candidate_parameter_names,
       apply_emt_parameter!,
       apply_emt_candidate!,
       evaluate_emt_candidate!,
       evaluate_emt_reduced!,
       evaluate_emt_batch!,
       reduce_emt_trace,
       evaluate_emt_study!,
       run_deck_emt_execution,
       electromagnetic_terminal_state,
       restart_emt_study!,
       restart_emt_checkpoint,
       read_emt_checkpoint,
       write_emt_checkpoint,
       PortableSnapshotFailure,
       PortableSnapshotArray,
       PortableSnapshotRecord,
       PortableSnapshotSection,
       PortableSnapshotMetadata,
       PortableEMTSnapshot,
       PortableSnapshotSectionDescriptor,
       PortableSnapshotDescriptor,
       portable_snapshot_array,
       portable_snapshot_array_values,
       portable_snapshot_bytes,
       portable_snapshot_descriptor,
       portable_emt_state_inventory,
       restore_portable_emt_state_inventory!,
       portable_emt_task_state,
       restore_portable_emt_task_state,
       portable_emt_task_scheduler_state_inventory,
       restore_portable_emt_task_scheduler_state_inventory!,
       portable_emt_hybrid_state_inventory,
       restore_portable_emt_hybrid_state_inventory!,
       PortableEMTHybridRestoreResult,
       restore_portable_emt_hybrid_snapshot,
       write_portable_emt_hybrid_snapshot,
       read_portable_emt_hybrid_snapshot,
       PortableEMTRestoreResult,
       capture_portable_emt_snapshot,
       restore_portable_emt_snapshot,
       write_portable_emt_workspace_snapshot,
       read_portable_emt_workspace_snapshot,
       inspect_portable_emt_snapshot,
       read_portable_emt_snapshot,
       read_portable_emt_snapshot_with_descriptor,
       write_portable_emt_snapshot,
       write_emt_restart_report,
       DeckEMTTrace,
       DeckUniversalMachineHorizon,
       DeckDirectMachineFleetHorizon,
       DeckSynchronousMachineHorizon,
       DeckSynchronousMachineFleetHorizon,
       DeckOVER16BoundaryPlan,
       FixedSourceLoadFlowResult,
       deck_steady_state_voltage_phasors,
       deck_fixed_source_load_flow,
       apply_deck_fixed_source_load_flow,
       write_fixed_source_load_flow_report,
       deck_over16_boundary_plan,
       deck_arrester_nonlinear_current_config,
       deck_piecewise_nonlinear_inductor_current_config,
       deck_hysteretic_inductor_nonlinear_current_config,
       deck_nonlinear_resistance_current_config,
       deck_saturated_transformer_nonlinear_current_config,
       deck_zinc_oxide_nonlinear_current_config,
       deck_fixed_step_horizon,
       deck_induction_machine_parameters,
       deck_induction_machine_state,
       deck_induction_machine_initial_state,
       deck_coupled_dq_machine_parameters,
       deck_coupled_dq_machine_state,
       deck_coupled_dq_machine_initial_state,
       deck_direct_current_machine_parameters,
       deck_direct_current_machine_state,
       deck_direct_current_machine_initial_state,
       deck_direct_current_machine_automatic_initialization,
       deck_separately_excited_dc_initialization,
       deck_universal_machine_type4_parameters,
       deck_universal_machine_type4_state,
       deck_universal_machine_type4_initial_state,
       deck_synchronous_machine_parameters,
       deck_synchronous_machine_initial_state,
       deck_output_step_indices,
       DeckOutputChannelMetadata,
       DeckFrequencyScanSchedule,
       DeckNetworkFrequencyScanPoint,
       DeckNetworkFrequencyScanStudy,
       DeckImpulseResponseFitControl,
       DeckLineFrequencyScanStudy,
       DeckLineImpulseResponseStudy,
       DeckAuxiliaryStudyRun,
       DeckCaseRun,
       DeckCaseSequenceRun,
       deck_frequency_scan_schedule,
       run_deck_auxiliary_studies,
       run_deck_case_sequence_emt,
       deck_case_sequence_result,
       write_deck_case_sequence_summary,
       deck_report_output_trace,
       deck_output_channel_metadata,
       AbstractSourceSignalProvider,
       IdentitySourceSignalProvider,
       TabulatedSourceSignalProvider,
       AnalyticSourceSlot,
       SourceSignalProgram,
       SourceSignalStageSample,
       source_signal_values,
       source_signal_analytic_values,
       source_signal_interpolation_active,
       source_signal_analytic_active,
       source_signal_provenance,
       write_source_signal_stage_report,
       initialize_step_context,
       run_reduced_feeder_inverter,
       run_reduced_feeder_inverter_with_trace,
       run_deck_emt,
       run_deck_universal_machine_horizon,
       run_deck_direct_machine_fleet_horizon,
       run_deck_synchronous_machine_horizon,
       run_deck_synchronous_machine_fleet_horizon,
       scheduled_trace,
       final_voltage_pu,
       final_output_pu,
       coupled_lumped_sequence_history_injection_elements,
       saturated_transformer_branch_augmented_step_context,
       saturated_transformer_winding_node_map,
       deck_trace_result,
       step!,
       step_with_over16_boundary!,
       run_nested_cable_primary_timestep!,
       run_deck_emt_with_over16_boundary,
       reduced_feeder_result,
       write_deck_trace_summary,
       write_universal_machine_horizon_csv,
       write_direct_machine_fleet_csv,
       write_unified_csv,
       write_unified_summary

const DECK_TRACE_SUMMARY_SCHEMA = "aimora.deck_trace_summary.v1"

Base.@kwdef struct UnifiedEMTConfig
    t_end_s::Float64 = 0.150
    dt_s::Float64 = 20e-6
    v_ll_base_v::Float64 = 4160.0
    source_pu::Float64 = 1.05
    initial_bus_pu::Float64 = 0.972
    feeder_r_pu::Float64 = 0.045
    feeder_l_pu_s::Float64 = 1.0e-7
    source_stiffness_pu::Float64 = 1000.0
    load_conductance_pu::Float64 = 0.843
    inverter_p0_pu::Float64 = 0.35
    inverter_p1_pu::Float64 = 0.80
    inverter_q_pu::Float64 = 0.05
end

struct DeckEMTTrace
    source::String
    dt_s::Float64
    t_end_s::Float64
    node_map::Dict{Symbol,Int}
    node_names::Vector{Symbol}
    element_names::Vector{Symbol}
    time_s::Vector{Float64}
    voltage_pu::Matrix{Float64}
    output_channel_names::Vector{Symbol}
    output_node_indices::Vector{Int}
    output_pu::Matrix{Float64}
    node_maximum_values::Vector{Float64}
    node_maximum_times_s::Vector{Float64}
    node_minimum_values::Vector{Float64}
    node_minimum_times_s::Vector{Float64}
    output_maximum_values::Vector{Float64}
    output_maximum_times_s::Vector{Float64}
    output_minimum_values::Vector{Float64}
    output_minimum_times_s::Vector{Float64}
end

struct DeckOutputChannelMetadata
    name::Symbol
    quantity::Symbol
    unit::String
    upper_name::String
    lower_name::String
end

struct EMTTerminalNodeState
    name::Symbol
    voltage_v::Float64
end

struct EMTTerminalBranchState
    name::Symbol
    kind::Symbol
    from_node::Symbol
    to_node::Symbol
    conductance_s::Float64
    history_current_a::Float64
    voltage_v::Float64
    current_a::Float64
    previous_current_a::Float64
    previous_voltage_v::Float64
    energy_j::Float64
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
end

struct EMTTerminalNonlinearState
    name::Symbol
    kind::Symbol
    from_node::Symbol
    to_node::Symbol
    voltage_v::Float64
    current_a::Float64
    flux_wb::Float64
    active_segment::Int
    energy_j::Float64
end

struct EMTTerminalSwitchState
    name::Symbol
    from_node::Symbol
    to_node::Symbol
    closed::Bool
    conductance_s::Float64
    voltage_v::Float64
    current_a::Float64
    power_w::Float64
    energy_j::Float64
    close_time_s::Float64
    open_time_s::Float64
end

struct EMTTerminalState
    source::String
    time_s::Float64
    nodes::Vector{EMTTerminalNodeState}
    branches::Vector{EMTTerminalBranchState}
    nonlinear_elements::Vector{EMTTerminalNonlinearState}
    switches::Vector{EMTTerminalSwitchState}
    physical_checks_passed::Bool
end

struct SeriesRLCAlteration
    branch_name::Symbol
    activation_time_s::Float64
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64

    function SeriesRLCAlteration(
        branch_name,
        activation_time_s::Real,
        resistance_ohm::Real,
        inductance_h::Real,
        capacitance_f::Real,
    )
        name = Symbol(branch_name)
        isempty(String(name)) &&
            throw(ArgumentError("series-RLC alteration branch name must not be empty"))
        activation_time = Float64(activation_time_s)
        resistance = Float64(resistance_ohm)
        inductance = Float64(inductance_h)
        capacitance = Float64(capacitance_f)
        isfinite(activation_time) && activation_time >= 0.0 ||
            throw(ArgumentError("series-RLC alteration time must be finite and nonnegative"))
        isfinite(resistance) && resistance >= 0.0 ||
            throw(ArgumentError("series-RLC alteration resistance must be finite and nonnegative"))
        isfinite(inductance) && inductance >= 0.0 ||
            throw(ArgumentError("series-RLC alteration inductance must be finite and nonnegative"))
        isfinite(capacitance) && capacitance > 0.0 ||
            throw(ArgumentError("series-RLC alteration capacitance must be finite and positive"))
        return new(name, activation_time, resistance, inductance, capacitance)
    end
end

struct SeriesRLCAlterationRecord
    branch_name::Symbol
    branch_index::Int
    requested_time_s::Float64
    applied_time_s::Float64
    applied_step_index::Int
    previous_resistance_ohm::Float64
    previous_inductance_h::Float64
    previous_capacitance_f::Float64
    resistance_ohm::Float64
    inductance_h::Float64
    capacitance_f::Float64
    previous_conductance_s::Float64
    conductance_s::Float64
    conductance_delta_s::Float64
    history_state_preserved::Bool
    network_refactor_count::Int
end

struct DeckEMTExecution{T}
    trace::T
    terminal_state::EMTTerminalState
    series_rlc_alteration_records::Vector{SeriesRLCAlterationRecord}
end

struct DeckUniversalMachineHorizon
    source::String
    machine_type::Int
    input_layout::Symbol
    parameter_basis::Symbol
    remanent_flux_enabled::Bool
    initialization_mode::Symbol
    maximum_shaft_mass_count::Int
    requested_terminal_coupling::Symbol
    time_s::Vector{Float64}
    output_values::Matrix{Float64}
    current_values::Matrix{Float64}
    history_currents::Matrix{Float64}
    current_substitution_values::Matrix{Float64}
    predicted_terminal_current_injections_A::Matrix{Float64}
    power_terminal_network_currents_A::Matrix{Float64}
    terminal_prediction_mutation_count::Int
    d_axis_flux::Vector{Float64}
    q_axis_flux::Vector{Float64}
    generated_torque::Vector{Float64}
    mechanical_speed_rad_s::Vector{Float64}
    mechanical_angle_rad::Vector{Float64}
    iteration_counts::Vector{Int}
    call_count::Int
    terminal_rhs_mutation_count::Int
    output_mutation_count::Int
    node_names::Vector{Symbol}
    power_terminal_voltages::Matrix{Float64}
    rotor_thevenin_matrices::Array{Float64,3}
    mechanical_speed_thevenin_rad_s::Vector{Float64}
    generated_torque_impedance::Vector{Float64}
    series_path_leakage_inductance_h::Float64
    effective_armature_leakage_inductance_h::Float64
    effective_compound_field_leakage_inductance_h::Float64
    compensated_voltage_values::Matrix{Float64}
    network_solve_count::Int
    network_correction_count::Int
    drive_source_values::Vector{Float64}
    excitation_source_values::Vector{Float64}
    excitation_source_frequency_hz::Float64
    q_axis_excitation_source_values::Vector{Float64}
    q_axis_excitation_source_frequency_hz::Float64
    network_control_mutation_count::Int
    report_mutation_count::Int
    complete_induction_machine_path::Bool
    deferred_effects::Vector{Symbol}
    report_output_names::Vector{Symbol}
    report_output_values::Matrix{Float64}
    control_system_step_count::Int
end

struct DeckSynchronousMachineHorizon
    source::String
    trace::DeckEMTTrace
    time_s::Vector{Float64}
    source_type::Int
    terminal_source_node_values::Vector{Int}
    terminal_voltage_values::Matrix{Float64}
    terminal_current_values::Matrix{Float64}
    terminal_open_circuit_voltage_values::Matrix{Float64}
    terminal_impedance_matrices::Array{Float64,3}
    terminal_admittance_matrices::Array{Float64,3}
    machine_output_values::Matrix{Float64}
    iteration_counts::Vector{Int}
    state::SynchronousMachineDynamicState
    network_solve_count::Int
    machine_call_count::Int
    terminal_rhs_mutation_count::Int
    output_mutation_count::Int
    control_output_names::Vector{Symbol}
    control_output_values::Matrix{Float64}
    machine_control_transfer_values::Matrix{Float64}
    field_voltage_multiplier_values::Vector{Float64}
    external_field_voltage_input_values_pu::Vector{Float64}
    mechanical_torque_multiplier_values::Matrix{Float64}
    total_applied_torque_values::Vector{Float64}
    mechanical_history_values::Matrix{Float64}
    control_system_step_count::Int
    machine_control_input_count::Int
    machine_control_transfer_count::Int
    complete_machine_control_coupling::Bool
    deferred_effects::Vector{Symbol}
end

struct DeckSynchronousMachineFleetHorizon
    source::String
    trace::DeckEMTTrace
    time_s::Vector{Float64}
    source_types::Vector{Int}
    terminal_node_indices::Matrix{Int}
    terminal_source_node_values::Matrix{Int}
    terminal_voltage_values::Array{Float64,3}
    terminal_current_values::Array{Float64,3}
    terminal_open_circuit_voltage_values::Array{Float64,3}
    terminal_impedance_matrices::Array{Float64,3}
    terminal_admittance_matrices::Array{Float64,3}
    machine_output_names::Vector{Symbol}
    machine_output_values::Matrix{Float64}
    mechanical_history_values::Vector{Matrix{Float64}}
    electrical_coefficient_values::Vector{Matrix{Float64}}
    saturation_enabled::Vector{Bool}
    d_axis_saturation_region_values::Matrix{Int}
    q_axis_saturation_region_values::Matrix{Int}
    saturation_refactor_count_values::Matrix{Int}
    states::Vector{SynchronousMachineDynamicState}
    network_solve_count::Int
    machine_call_count::Int
    terminal_rhs_mutation_count::Int
    output_mutation_count::Int
    control_output_names::Vector{Symbol}
    control_output_values::Matrix{Float64}
    machine_control_transfer_values::Matrix{Float64}
    field_voltage_multiplier_values::Matrix{Float64}
    external_field_voltage_input_values_pu::Matrix{Float64}
    mechanical_torque_multiplier_values::Vector{Matrix{Float64}}
    total_applied_torque_values::Matrix{Float64}
    control_system_step_count::Int
    machine_control_input_count::Int
    machine_control_transfer_count::Int
    complete_machine_network_coupling::Bool
    complete_machine_control_coupling::Bool
    deferred_effects::Vector{Symbol}
end

function _sampled_trace_extrema(
    values::AbstractMatrix{<:Real},
    time_s::AbstractVector{<:Real},
)
    size(values, 2) == length(time_s) ||
        throw(ArgumentError("trace values and times must have equal sample counts"))
    variable_count = size(values, 1)
    maxima = fill(-Inf, variable_count)
    maxima_times_s = zeros(Float64, variable_count)
    minima = fill(Inf, variable_count)
    minima_times_s = zeros(Float64, variable_count)
    for sample_index in eachindex(time_s)
        sample_time_s = Float64(time_s[sample_index])
        for variable_index in 1:variable_count
            value = Float64(values[variable_index, sample_index])
            if value > maxima[variable_index]
                maxima[variable_index] = value
                maxima_times_s[variable_index] = sample_time_s
            end
            if value < minima[variable_index]
                minima[variable_index] = value
                minima_times_s[variable_index] = sample_time_s
            end
        end
    end
    return (
        maximum_values = maxima,
        maximum_times_s = maxima_times_s,
        minimum_values = minima,
        minimum_times_s = minima_times_s,
    )
end

function DeckEMTTrace(
    source::String,
    dt_s::Float64,
    t_end_s::Float64,
    node_map::Dict{Symbol,Int},
    node_names::Vector{Symbol},
    element_names::Vector{Symbol},
    time_s::Vector{Float64},
    voltage_pu::Matrix{Float64},
    output_channel_names::Vector{Symbol},
    output_node_indices::Vector{Int},
    output_pu::Matrix{Float64},
)
    node_extrema = _sampled_trace_extrema(voltage_pu, time_s)
    output_extrema = _sampled_trace_extrema(output_pu, time_s)
    return DeckEMTTrace(
        source,
        dt_s,
        t_end_s,
        node_map,
        node_names,
        element_names,
        time_s,
        voltage_pu,
        output_channel_names,
        output_node_indices,
        output_pu,
        node_extrema.maximum_values,
        node_extrema.maximum_times_s,
        node_extrema.minimum_values,
        node_extrema.minimum_times_s,
        output_extrema.maximum_values,
        output_extrema.maximum_times_s,
        output_extrema.minimum_values,
        output_extrema.minimum_times_s,
    )
end

mutable struct DeckOVER16BoundaryPlan
    source::Symbol
    branch_names::Vector{Symbol}
    branch_kinds::Vector{Symbol}
    branch_from_node_names::Vector{Symbol}
    branch_to_node_names::Vector{Symbol}
    branch_from_node_indices::Vector{Int}
    branch_to_node_indices::Vector{Int}
    branch_conductance_values::Vector{Float64}
    branch_resistance_values::Vector{Float64}
    branch_inductance_values::Vector{Float64}
    branch_capacitance_values::Vector{Float64}
    branch_previous_current_values::Vector{Float64}
    branch_previous_voltage_values::Vector{Float64}
    branch_line_numbers::Vector{Int}
    branch_count::Int
    bergeron_line_names::Vector{Symbol}
    bergeron_line_line_numbers::Vector{Int}
    bergeron_line_from_node_names::Vector{Symbol}
    bergeron_line_to_node_names::Vector{Symbol}
    bergeron_line_from_node_indices::Vector{Int}
    bergeron_line_to_node_indices::Vector{Int}
    bergeron_line_surge_impedance_values::Vector{Float64}
    bergeron_line_surge_admittance_values::Vector{Float64}
    bergeron_line_travel_time_s_values::Vector{Float64}
    bergeron_line_dt_s_values::Vector{Float64}
    bergeron_line_attenuation_values::Vector{Float64}
    bergeron_line_delay_step_counts::Vector{Int}
    bergeron_line_write_indices::Vector{Int}
    bergeron_line_history_current_from_values::Vector{Float64}
    bergeron_line_history_current_to_values::Vector{Float64}
    bergeron_line_terminal_voltage_from_values::Vector{Float64}
    bergeron_line_terminal_voltage_to_values::Vector{Float64}
    bergeron_line_terminal_current_from_values::Vector{Float64}
    bergeron_line_terminal_current_to_values::Vector{Float64}
    bergeron_line_traveling_wave_from_values::Vector{Float64}
    bergeron_line_traveling_wave_to_values::Vector{Float64}
    bergeron_line_count::Int
    over2_branch_names::Vector{Symbol}
    over2_branch_line_numbers::Vector{Int}
    over2_branch_kinds::Vector{Symbol}
    over2_branch_layout_kinds::Vector{Symbol}
    over2_branch_source_kinds::Vector{Symbol}
    over2_branch_reference_kinds::Vector{Symbol}
    over2_branch_reference_names::Vector{Symbol}
    over2_branch_reference_line_numbers::Vector{Int}
    over2_branch_from_node_names::Vector{Symbol}
    over2_branch_to_node_names::Vector{Symbol}
    over2_branch_from_node_indices::Vector{Int}
    over2_branch_to_node_indices::Vector{Int}
    over2_branch_raw_resistance_values::Vector{Float64}
    over2_branch_raw_inductance_values::Vector{Float64}
    over2_branch_raw_capacitance_values::Vector{Float64}
    over2_branch_conductance_values::Vector{Float64}
    over2_branch_resistance_values::Vector{Float64}
    over2_branch_inductance_values::Vector{Float64}
    over2_branch_capacitance_values::Vector{Float64}
    over2_branch_output_codes::Vector{Int}
    over2_branch_count::Int
    output_channel_names::Vector{Symbol}
    output_node_names::Vector{Symbol}
    output_node_indices::Vector{Int}
    output_channel_line_numbers::Vector{Int}
    branch_voltage_output_names::Vector{Symbol}
    branch_voltage_branch_names::Vector{Symbol}
    branch_voltage_branch_indices::Vector{Int}
    branch_voltage_output_line_numbers::Vector{Int}
    branch_current_output_names::Vector{Symbol}
    branch_current_branch_names::Vector{Symbol}
    branch_current_branch_indices::Vector{Int}
    branch_current_output_line_numbers::Vector{Int}
    branch_power_output_names::Vector{Symbol}
    branch_power_branch_names::Vector{Symbol}
    branch_power_branch_indices::Vector{Int}
    branch_power_output_line_numbers::Vector{Int}
    over15_output_request_names::Vector{Symbol}
    over15_output_request_output_kinds::Vector{Symbol}
    over15_output_request_request_kinds::Vector{Symbol}
    over15_output_request_layout_kinds::Vector{Symbol}
    over15_output_request_line_numbers::Vector{Int}
    over15_output_request_output_codes::Vector{Int}
    over15_output_request_node_names::Vector{Symbol}
    over15_output_request_node_indices::Vector{Int}
    over15_output_request_branch_names::Vector{Symbol}
    over15_output_request_branch_indices::Vector{Int}
    switch_names::Vector{Symbol}
    switch_line_numbers::Vector{Int}
    switch_from_node_names::Vector{Symbol}
    switch_to_node_names::Vector{Symbol}
    switch_from_node_indices::Vector{Int}
    switch_to_node_indices::Vector{Int}
    switch_close_time_s_values::Vector{Float64}
    switch_open_time_s_values::Vector{Float64}
    switch_initially_closed_flags::Vector{Bool}
    switch_on_conductance_values::Vector{Float64}
    switch_off_conductance_values::Vector{Float64}
    switch_count::Int
    over5_switch_names::Vector{Symbol}
    over5_switch_line_numbers::Vector{Int}
    over5_switch_from_node_names::Vector{Symbol}
    over5_switch_to_node_names::Vector{Symbol}
    over5_switch_from_node_indices::Vector{Int}
    over5_switch_to_node_indices::Vector{Int}
    over5_switch_layout_kinds::Vector{Symbol}
    over5_switch_raw_close_time_s_values::Vector{Float64}
    over5_switch_raw_open_time_s_values::Vector{Float64}
    over5_switch_close_time_s_values::Vector{Float64}
    over5_switch_open_time_s_values::Vector{Float64}
    over5_switch_initially_closed_flags::Vector{Bool}
    over5_switch_measuring_flags::Vector{Bool}
    over5_switch_closed_markers::Vector{String}
    over5_switch_marker_texts::Vector{String}
    over5_switch_type_values::Vector{Int}
    over5_switch_critical_current_values::Vector{Float64}
    over5_switch_random_opening_standard_deviation_s_values::Vector{Float64}
    over5_switch_on_conductance_values::Vector{Float64}
    over5_switch_off_conductance_values::Vector{Float64}
    over5_switch_output_codes::Vector{Int}
    over5_switch_count::Int
    source_names::Vector{Symbol}
    source_node_names::Vector{Symbol}
    source_node_values::Vector{Int}
    source_iform_values::Vector{Int}
    source_line_numbers::Vector{Int}
    source_layout_kinds::Vector{Symbol}
    source_tstart_values::Vector{Float64}
    source_tstop_values::Vector{Float64}
    source_crest_values::Vector{Float64}
    source_time1_values::Vector{Float64}
    source_time2_values::Vector{Float64}
    source_sfreq_values::Vector{Float64}
    source_row_count::Int
    source_card_kinds::Vector{Symbol}
    source_card_values::Vector{Vector{Float64}}
    source_card_provided_value_counts::Vector{Int}
    source_card_line_numbers::Vector{Int}
    source_card_row_count::Int
    source_interpolation_values::Vector{Vector{Float64}}
    source_interpolation_provided_value_counts::Vector{Int}
    source_interpolation_line_numbers::Vector{Int}
    source_interpolation_row_count::Int
    source_tacs_override_positions::Vector{Int}
    source_tacs_override_xtcs_indices::Vector{Int}
    source_tacs_override_line_numbers::Vector{Int}
    source_tacs_override_count::Int
    source_analytic_values::Vector{Vector{Float64}}
    source_analytic_provided_value_counts::Vector{Int}
    source_analytic_line_numbers::Vector{Int}
    source_analytic_row_count::Int
    fortran_scope::Tuple{Vararg{Symbol}}
    deferred_effects::Tuple{Vararg{Symbol}}
end

mutable struct ControlSystemSupplementalDeviceRuntime
    state::OVER16CSUPState
    rows::Vector{OVER16CSUPDeviceRow}
    next_indices::Vector{Int}
    ordinary_signal_names::Vector{Symbol}
    device_output_names::Vector{Symbol}
    device_output_slots::Vector{Int}
    transport_delay_pointers::Vector{Int}
    initialize_transport_delay_from_input::BitVector
    initialize_rms_from_input::BitVector
    initialized::Bool
    executed_step_count::Int
end

struct ControlSystemFrequencyInitialization
    frequency_hz::Float64
    signal_names::Vector{Symbol}
    signal_phasors::Vector{ComplexF64}
    source_names::Vector{Symbol}
    device_output_names::Vector{Symbol}
    history_mutation_count::Int
end

struct WindowedConstantControlSignal
    name::Symbol
    value::Float64
    start_time_s::Float64
    stop_time_s::Float64
end

struct ControlWaveformSignal
    name::Symbol
    waveform_type::Int
    amplitude::Float64
    cycle_s::Float64
    pulse_width_s::Float64
    start_time_s::Float64
    stop_time_s::Float64
end

const ControlObservedSwitch =
    Union{IdealSwitch,TimeSwitch,CurrentZeroSwitch,TACSControlledSwitch}

struct ControlSwitchObservation
    name::Symbol
    source_type::Int
    switch::ControlObservedSwitch
    ordinary_switch_index::Int
end

mutable struct ControlSystemNetworkRuntime
    state::ControlSystemExecutionState
    signal_slot_names::Vector{Symbol}
    network_voltage_source_names::Vector{Symbol}
    network_voltage_source_node_indices::Vector{Int}
    switch_observations::Vector{ControlSwitchObservation}
    switch_elements::Vector{TACSControlledSwitch}
    switch_signal_names::Vector{Symbol}
    switch_clamp_signal_names::Vector{Union{Missing,Symbol}}
    switch_output_indices::Vector{Int}
    switch_output_names::Vector{Symbol}
    control_output_names::Vector{Symbol}
    windowed_constant_sources::Vector{WindowedConstantControlSignal}
    sinusoidal_sources::Vector{SinusoidalControlSignal}
    waveform_sources::Vector{ControlWaveformSignal}
    control_expressions::Vector{ControlExpressionRuntime}
    supplemental_devices::Union{Nothing,ControlSystemSupplementalDeviceRuntime}
    frequency_initializations::Vector{ControlSystemFrequencyInitialization}
    executed_step_count::Int
    feedback_application_count::Int
    halted::Bool
    last_diagnostic::Union{Nothing,String}
end

mutable struct SourceFunctionNetworkRuntime
    state::OVER16SourceCardState
    plan::DeckOVER16BoundaryPlan
    slot_values::Vector{Base.RefValue{Float64}}
    row_slot_values::Vector{Base.RefValue{Float64}}
    dynamic_row_indices::Vector{Int}
    signal_provider::AbstractSourceSignalProvider
    internal_analytic_requested::Bool
    executed_step_count::Int
    card_read_count::Int
    signal_synchronization_count::Int
    external_signal_count::Int
    tacs_override_count::Int
    analytic_execution_count::Int
    stage_samples::Vector{SourceSignalStageSample}
    last_accepted_time_s::Float64
    next_input_row_index::Int
end

@enum ElectromagneticHistoryKind::UInt8 begin
    SERIES_RL_HISTORY = 1
    SERIES_RLC_HISTORY = 2
    CAPACITOR_HISTORY = 3
    COUPLED_INDUCTIVE_HISTORY = 4
    COUPLED_SERIES_RL_HISTORY = 5
    BREQIV_HISTORY = 6
end

struct ElectromagneticHistoryExecutionPlan
    element_indices::Vector{Int}
    kinds::Vector{ElectromagneticHistoryKind}
    batch_indices::Vector{Int}
    series_rl_branches::Vector{SeriesRLBranch}
    series_rlc_branches::Vector{SeriesRLCBranch}
    capacitor_branches::Vector{CapacitorBranch}
    coupled_inductive_branches::Vector{CoupledInductiveBranch}
    coupled_series_rl_branches::Vector{CoupledSeriesRLBranch}
    breqiv_injections::Vector{BreqivHistoryInjection}
end

mutable struct EMTStepContext{S<:NodalSystem}
    system::S
    analytic_source_signals::Vector{AnalyticSourceSignal}
    analytic_source_names::Vector{Symbol}
    saturated_transformer_nonlinear_slope_branch_batch::Vector{SaturatedTransformerNonlinearSlopeBranch}
    electromagnetic_history_plan::ElectromagneticHistoryExecutionPlan
    electromagnetic_history_rhs::Vector{Float64}
    reference_augmented_admittance::Matrix{Float64}
    reference_augmented_rhs::Vector{Float64}
    reference_sparse_columns::Vector{Int}
    reference_sparse_values::Vector{Float64}
    reference_sparse_row_ends::Vector{Int}
    sparse_node_group_workspace::SparseNodeGroupSolveWorkspace
    source::String
    dt_s::Float64
    t_end_s::Float64
    t_s::Float64
    step_index::Int
    step_count::Int
    node_map::Dict{Symbol,Int}
    node_names::Vector{Symbol}
    element_names::Vector{Symbol}
    time_s::Vector{Float64}
    voltage_pu::Matrix{Float64}
    output_channel_names::Vector{Symbol}
    source_function_runtime::Union{Nothing,SourceFunctionNetworkRuntime}
    control_system_runtime::Union{Nothing,ControlSystemNetworkRuntime}
    output_node_indices::Vector{Int}
    branch_voltage_output_branch_indices::Vector{Int}
    branch_current_output_branch_indices::Vector{Int}
    branch_power_output_branch_indices::Vector{Int}
    switch_voltage_output_switch_indices::Vector{Int}
    switch_current_output_switch_indices::Vector{Int}
    deck_time_switch_count::Int
    deck_time_switch_names::Vector{Symbol}
    deck_time_switch_from_node_indices::Vector{Int}
    deck_time_switch_to_node_indices::Vector{Int}
    deck_time_switch_close_time_s_values::Vector{Float64}
    deck_time_switch_open_time_s_values::Vector{Float64}
    deck_time_switch_initially_closed_flags::Vector{Bool}
    deck_time_switch_on_conductance_values::Vector{Float64}
    deck_time_switch_off_conductance_values::Vector{Float64}
    deck_over5_switch_output_codes::Vector{Int}
    switch_closed_step_flags::Vector{Int}
    switch_conductance_step_values::Vector{Float64}
    switch_voltage_step_values::Vector{Float64}
    switch_current_step_values::Vector{Float64}
    switch_power_step_values::Vector{Float64}
    branch_previous_power_values::Vector{Float64}
    branch_power_history_valid::BitVector
    branch_energy_values::Vector{Float64}
    switch_previous_power_values::Vector{Float64}
    switch_power_history_valid::BitVector
    switch_energy_values::Vector{Float64}
    output_pu::Matrix{Float64}
    output_step_values::Matrix{Float64}
    node_maximum_values::Vector{Float64}
    node_maximum_times_s::Vector{Float64}
    node_minimum_values::Vector{Float64}
    node_minimum_times_s::Vector{Float64}
    output_maximum_values::Vector{Float64}
    output_maximum_times_s::Vector{Float64}
    output_minimum_values::Vector{Float64}
    output_minimum_times_s::Vector{Float64}
    recorded_step_indices::Vector{Int}
    trace_write_index::Int
    series_rlc_alterations::Vector{SeriesRLCAlteration}
    series_rlc_alteration_branch_indices::Vector{Int}
    next_series_rlc_alteration_index::Int
    series_rlc_alteration_records::Vector{SeriesRLCAlterationRecord}
    series_rlc_network_refactor_count::Int
end

struct PreparedEMTStudy{R,P}
    runtime_template::R
    parsed::P
end

"""Atomic model-owned machine state together with a one-step coupled EMT probe."""
struct PreparedMachineEMTStudy{I,H,P}
    machine_family::Symbol
    initialization_state::I
    one_step_horizon::H
    parsed::P
end

"""Explicit ownership binding between one admitted average-value converter and its three network terminal/current-source phases."""
struct AverageValueGridFollowingConverterInitialization{P<:InverterParams}
    parameters::P
    terminal_nodes::NTuple{3,Symbol}
    current_source_names::NTuple{3,Symbol}

    function AverageValueGridFollowingConverterInitialization(
        parameters::P;
        terminal_nodes::NTuple{3,Symbol},
        current_source_names::NTuple{3,Symbol},
    ) where {P<:InverterParams}
        length(unique(terminal_nodes)) == 3 || throw(ArgumentError(
            "average-value converter initialization requires three distinct phase terminals",
        ))
        length(unique(current_source_names)) == 3 || throw(ArgumentError(
            "average-value converter initialization requires three distinct phase-current owners",
        ))
        return new{P}(parameters, terminal_nodes, current_source_names)
    end
end

"""Accepted harmonic network plus the model-owned average-value converter/controller state and its one-step equilibrium probe."""
struct PreparedAverageValueGridFollowingEMTStudy{N,E,C}
    network::N
    equilibrium::E
    converter::C
end

"""Frequency-domain formulation used to construct an instantaneous-EMT state."""
abstract type AbstractEMTHarmonicFormulation end

"""Evaluate every admitted model at the requested physical frequency."""
struct PhysicalFrequencyFormulation <: AbstractEMTHarmonicFormulation end

"""Evaluate each admitted discrete recurrence at the requested physical waveform frequency and timestep."""
struct TimestepMatchedFormulation <: AbstractEMTHarmonicFormulation
    timestep_s::Float64

    function TimestepMatchedFormulation(timestep_s::Real)
        timestep = Float64(timestep_s)
        isfinite(timestep) && timestep > 0.0 || throw(ArgumentError(
            "timestep-matched initialization requires a finite positive timestep",
        ))
        return new(timestep)
    end
end

"""Quantity-specific numerical limits for accepting one EMT initialization."""
Base.@kwdef struct EMTInitializationTolerances
    voltage_absolute_v::Float64 = 1.0e-9
    voltage_relative::Float64 = 1.0e-10
    current_absolute_a::Float64 = 1.0e-10
    current_relative::Float64 = 1.0e-10
    power_absolute_w::Float64 = 1.0e-7
    power_relative::Float64 = 1.0e-9
    energy_absolute_j::Float64 = 1.0e-9
    energy_relative::Float64 = 1.0e-8
    flux_absolute_wb::Float64 = 1.0e-12
    operating_point_maximum_iterations::Int = 2000
    scaled_residual_maximum::Float64 = 1.0
    rank_relative_threshold_multiplier::Float64 = 10.0
    maximum_condition_estimate::Float64 = 1.0e12
    no_artificial_transient_normalized_rms::Float64 = 1.0e-8
    first_step_scaled_discontinuity::Float64 = 1.0e-8
end

"""One versioned operating-point quantity with an explicit conversion into its target SI orientation."""
struct OperatingPointQuantity
    asset::Symbol
    quantity::Symbol
    phase::Symbol
    value::ComplexF64
    unit::String
    basis::String
    scale_to_si::Float64
    orientation::String
    orientation_sign::Float64
    absolute_uncertainty::Float64
    provenance::ParameterProvenance

    function OperatingPointQuantity(
        asset::Symbol,
        quantity::Symbol,
        value::Number;
        phase::Symbol=:not_applicable,
        unit::AbstractString,
        basis::AbstractString="absolute_si",
        scale_to_si::Real=1.0,
        orientation::AbstractString,
        orientation_sign::Real=1.0,
        absolute_uncertainty::Real=0.0,
        provenance::ParameterProvenance,
    )
        complex_value = ComplexF64(value)
        scale = Float64(scale_to_si)
        sign = Float64(orientation_sign)
        uncertainty = Float64(absolute_uncertainty)
        isempty(String(asset)) && throw(ArgumentError("operating-point asset cannot be empty"))
        isempty(String(quantity)) && throw(ArgumentError("operating-point quantity cannot be empty"))
        all(isfinite, (real(complex_value), imag(complex_value))) || throw(ArgumentError(
            "operating-point quantity value must be finite",
        ))
        isfinite(scale) && scale > 0.0 || throw(ArgumentError(
            "operating-point SI scale must be finite and positive",
        ))
        sign in (-1.0, 1.0) || throw(ArgumentError(
            "operating-point orientation sign must be -1 or 1",
        ))
        isfinite(uncertainty) && uncertainty >= 0.0 || throw(ArgumentError(
            "operating-point uncertainty must be finite and nonnegative",
        ))
        isempty(strip(unit)) && throw(ArgumentError("operating-point unit cannot be empty"))
        isempty(strip(basis)) && throw(ArgumentError("operating-point basis cannot be empty"))
        isempty(strip(orientation)) && throw(ArgumentError(
            "operating-point orientation cannot be empty",
        ))
        return new(
            asset,
            quantity,
            phase,
            complex_value,
            String(unit),
            String(basis),
            scale,
            String(orientation),
            sign,
            uncertainty,
            provenance,
        )
    end
end

"""Immutable versioned operating point that may be mapped into an EMT realization."""
struct EMTOperatingPoint
    schema::UInt32
    source_representation::Symbol
    project_signature::String
    settings_signature::String
    model_signature::String
    source_state_signature::String
    frequency_hz::Float64
    time_origin_s::Float64
    phase_order::NTuple{3,Symbol}
    quantities::Vector{OperatingPointQuantity}

    function EMTOperatingPoint(
        source_representation::Symbol,
        project_signature::AbstractString,
        settings_signature::AbstractString,
        model_signature::AbstractString,
        frequency_hz::Real,
        quantities::AbstractVector{<:OperatingPointQuantity};
        schema::Integer=1,
        time_origin_s::Real=0.0,
        phase_order::NTuple{3,Symbol}=(:a, :b, :c),
        source_state_signature::AbstractString="",
    )
        schema_value = UInt32(schema)
        schema_value == UInt32(1) || throw(ArgumentError(
            "unsupported EMT operating-point schema $schema_value",
        ))
        isempty(String(source_representation)) && throw(ArgumentError(
            "operating-point source representation cannot be empty",
        ))
        frequency = Float64(frequency_hz)
        time_origin = Float64(time_origin_s)
        isfinite(frequency) && frequency >= 0.0 || throw(ArgumentError(
            "operating-point frequency must be finite and nonnegative",
        ))
        isfinite(time_origin) || throw(ArgumentError(
            "operating-point time origin must be finite",
        ))
        length(unique(phase_order)) == 3 || throw(ArgumentError(
            "operating-point phase order must contain three distinct phases",
        ))
        signatures = String.((project_signature, settings_signature, model_signature))
        all(value -> !isempty(strip(value)), signatures) || throw(ArgumentError(
            "operating-point signatures cannot be empty",
        ))
        owned_quantities = OperatingPointQuantity[quantity for quantity in quantities]
        isempty(owned_quantities) && throw(ArgumentError(
            "operating point must contain at least one mapped quantity",
        ))
        keys = [(quantity.asset, quantity.quantity, quantity.phase) for quantity in owned_quantities]
        length(unique(keys)) == length(keys) || throw(ArgumentError(
            "operating-point asset/quantity/phase keys must be unique",
        ))
        source_state = String(source_state_signature)
        source_representation === :accepted_emt_initialization &&
            isempty(strip(source_state)) && throw(ArgumentError(
            "an accepted-EMT operating point requires its source-state signature",
        ))
        return new(
            schema_value,
            source_representation,
            signatures...,
            source_state,
            frequency,
            time_origin,
            phase_order,
            owned_quantities,
        )
    end
end

"""Complete request for physical or timestep-matched EMT state construction."""
struct EMTInitializationRequest{F<:AbstractEMTHarmonicFormulation,O}
    formulation::F
    frequency_hz::Float64
    frequency_grid_hz::Vector{Float64}
    time_origin_s::Float64
    operating_point::O
    project_signature::String
    settings_signature::String
    model_signature::String
    tolerances::EMTInitializationTolerances
end

function EMTInitializationRequest(
    formulation::F;
    frequency_hz::Real,
    frequency_grid_hz::AbstractVector{<:Real}=Float64[frequency_hz],
    time_origin_s::Real=0.0,
    operating_point::O=nothing,
    project_signature::AbstractString,
    settings_signature::AbstractString,
    model_signature::AbstractString,
    tolerances::EMTInitializationTolerances=EMTInitializationTolerances(),
) where {F<:AbstractEMTHarmonicFormulation,O}
    frequency = Float64(frequency_hz)
    time_origin = Float64(time_origin_s)
    isfinite(frequency) && frequency >= 0.0 || throw(ArgumentError(
        "initialization frequency must be finite and nonnegative",
    ))
    isfinite(time_origin) || throw(ArgumentError(
        "initialization time origin must be finite",
    ))
    frequency_grid = sort!(unique(Float64.(frequency_grid_hz)))
    isempty(frequency_grid) && throw(ArgumentError(
        "initialization frequency grid must contain at least one point",
    ))
    all(value -> isfinite(value) && value >= 0.0, frequency_grid) ||
        throw(ArgumentError(
            "initialization frequency grid must be finite and nonnegative",
        ))
    if formulation isa TimestepMatchedFormulation
        all(value -> 2.0 * pi * value * formulation.timestep_s < pi,
            frequency_grid) || throw(ArgumentError(
            "timestep-matched initialization frequency grid must remain below Nyquist",
        ))
    end
    signatures = String.((project_signature, settings_signature, model_signature))
    all(value -> !isempty(strip(value)), signatures) || throw(ArgumentError(
        "initialization signatures cannot be empty",
    ))
    _validate_emt_initialization_tolerances(tolerances)
    return EMTInitializationRequest(
        formulation,
        frequency,
        frequency_grid,
        time_origin,
        operating_point,
        signatures...,
        tolerances,
    )
end

function _validate_emt_initialization_tolerances(tolerances::EMTInitializationTolerances)
    tolerances.operating_point_maximum_iterations > 0 || throw(ArgumentError(
        "EMT operating-point maximum iterations must be positive",
    ))
    values = Float64[
        tolerances.voltage_absolute_v,
        tolerances.voltage_relative,
        tolerances.current_absolute_a,
        tolerances.current_relative,
        tolerances.power_absolute_w,
        tolerances.power_relative,
        tolerances.energy_absolute_j,
        tolerances.energy_relative,
        tolerances.flux_absolute_wb,
        tolerances.scaled_residual_maximum,
        tolerances.rank_relative_threshold_multiplier,
        tolerances.maximum_condition_estimate,
        tolerances.no_artificial_transient_normalized_rms,
        tolerances.first_step_scaled_discontinuity,
    ]
    all(value -> isfinite(value) && value > 0.0, values) || throw(ArgumentError(
        "EMT initialization tolerances must be finite and positive",
    ))
    return tolerances
end

"""One unscaled and scaled initialization residual owned by a physical equation or mapping."""
struct EMTInitializationResidual
    owner::Symbol
    quantity::Symbol
    unit::String
    value::Float64
    absolute_tolerance::Float64
    relative_tolerance::Float64
    reference_scale::Float64
    uncertainty_allowance::Float64
    scaled_value::Float64
    passed::Bool
end

"""One explicit source-to-target operating-point quantity conversion."""
struct OperatingPointMappingRecord
    asset::Symbol
    quantity::Symbol
    phase::Symbol
    source_value::ComplexF64
    target_value_si::ComplexF64
    source_unit::String
    target_unit::String
    basis::String
    orientation::String
    scale_to_si::Float64
    orientation_sign::Float64
    absolute_uncertainty_si::Float64
    residual::Float64
    constraint_current_phasor_a::ComplexF64
    passed::Bool
end

"""Connected-component, rank, conditioning, and residual diagnostics for the harmonic network."""
struct EMTInitializationTopologyReport
    node_count::Int
    reduced_node_count::Int
    switch_node_groups::Vector{Vector{Int}}
    connected_components::Vector{Vector{Int}}
    referenced_components::BitVector
    unreferenced_components::Vector{Vector{Int}}
    numerical_rank::Int
    condition_estimate::Float64
    maximum_residual_a::Float64
    relative_residual::Float64
    classification::Symbol
end

"""One physical-frequency or timestep-matched network equilibrium with strict topology diagnostics."""
struct EMTInitializationFrequencyPoint
    formulation::Symbol
    frequency_assignment::Symbol
    physical_frequency_hz::Float64
    reactive_angular_frequency_rad_s::Float64
    node_physical_frequencies_hz::Vector{Float64}
    node_frequency_source_row_indices::Vector{Int}
    source_frequency_successor_indices::Vector{Int}
    frequency_subnetwork_count::Int
    node_voltage_phasors::Vector{ComplexF64}
    source_injection_phasors::Vector{ComplexF64}
    operating_constraint_current_phasors::Vector{ComplexF64}
    topology::EMTInitializationTopologyReport
    admittance_symmetry_max_abs_error::Float64
    minimum_dissipative_eigenvalue_s::Float64
    passed::Bool
end

"""One independently defined unwanted-transient metric over an undisturbed window."""
struct NoArtificialTransientMetric
    quantity::Symbol
    unit::String
    window_start_s::Float64
    window_end_s::Float64
    normalized_rms::Float64
    maximum_scaled_discontinuity::Float64
    low_frequency_envelope_drift::Float64
    energy_residual_j::Float64
    threshold::Float64
    passed::Bool
end

"""One state family proven present and initialized in the accepted EMT candidate."""
struct EMTInitializationStateRecord
    owner::Symbol
    state_family::Symbol
    instance_count::Int
    initialization_basis::Symbol

    function EMTInitializationStateRecord(
        owner::Symbol,
        state_family::Symbol,
        instance_count::Integer,
        initialization_basis::Symbol,
    )
        count = Int(instance_count)
        count > 0 || throw(ArgumentError(
            "an EMT initialization state record must describe at least one instance",
        ))
        return new(owner, state_family, count, initialization_basis)
    end
end

"""Typed initialization refusal with the exact scientific owner and context."""
struct EMTInitializationFailure
    code::Symbol
    owner::Symbol
    quantity::Symbol
    message::String
    context::NamedTuple
end

"""Complete accepted-or-refused initialization evidence without private factor objects."""
struct EMTInitializationReport
    status::Symbol
    formulation::Symbol
    frequency_hz::Float64
    time_origin_s::Float64
    topology::EMTInitializationTopologyReport
    frequency_scan::Vector{EMTInitializationFrequencyPoint}
    residuals::Vector{EMTInitializationResidual}
    mappings::Vector{OperatingPointMappingRecord}
    state_inventory::Vector{EMTInitializationStateRecord}
    initialized_state_owners::Vector{Symbol}
    unsupported_state_owners::Vector{Symbol}
    transient_metrics::Vector{NoArtificialTransientMetric}
    warnings::Vector{String}
    project_signature::String
    settings_signature::String
    model_signature::String
    deterministic_state_signature::String
end

"""Atomic initialization result; a refused result never contains a mutable prepared study."""
struct EMTInitializationResult{P}
    prepared::P
    report::EMTInitializationReport
    failure::Union{Nothing,EMTInitializationFailure}
end

"""One direct or prior-case-mapped initialization instruction for a deck sequence case."""
struct DeckCaseInitialization{R<:EMTInitializationRequest}
    request::R
    mapped_from_case_index::Union{Nothing,Int}

    function DeckCaseInitialization(
        request::R;
        mapped_from_case_index::Union{Nothing,Integer}=nothing,
    ) where {R<:EMTInitializationRequest}
        source_case = mapped_from_case_index === nothing ? nothing :
            Int(mapped_from_case_index)
        source_case === nothing || source_case > 0 || throw(ArgumentError(
            "mapped source case index must be positive",
        ))
        if source_case !== nothing
            request.operating_point isa EMTOperatingPoint || throw(ArgumentError(
                "mapped case initialization requires a typed operating point",
            ))
            request.operating_point.source_representation ===
                :accepted_emt_initialization || throw(ArgumentError(
                "mapped case initialization must identify accepted EMT state as its source",
            ))
        end
        return new{R}(request, source_case)
    end
end

initialization_accepted(result::EMTInitializationResult) =
    result.failure === nothing && result.report.status === :accepted

_emt_harmonic_formulation_symbol(::PhysicalFrequencyFormulation) =
    :physical_frequency
_emt_harmonic_formulation_symbol(::TimestepMatchedFormulation) =
    :timestep_matched

function _emt_reactive_angular_frequency(
    ::PhysicalFrequencyFormulation,
    physical_angular_frequency::Float64,
)
    return physical_angular_frequency
end

function _emt_reactive_angular_frequency(
    formulation::TimestepMatchedFormulation,
    physical_angular_frequency::Float64,
)
    angle = 0.5 * physical_angular_frequency * formulation.timestep_s
    abs(angle) < 0.5 * pi || throw(ArgumentError(
        "timestep-matched initialization frequency must remain below Nyquist",
    ))
    value = (2.0 / formulation.timestep_s) * tan(angle)
    isfinite(value) || throw(ArgumentError(
        "timestep-matched reactive frequency is nonfinite",
    ))
    return value
end

struct EMTRestartMutationRecord
    kind::Symbol
    target::Symbol
    line_no::Int
    parameter_names::Vector{Symbol}
    previous_values::Vector{Float64}
    applied_values::Vector{Float64}
end

struct EMTRestartRun{T,R}
    trace::T
    request::R
    checkpoint_time_s::Float64
    final_time_s::Float64
    appended_step_count::Int
    mutation_records::Vector{EMTRestartMutationRecord}
    checkpoint_state_error::Float64
    final_kcl_error::Float64
end

mutable struct EMTStudyWorkspace{R,P}
    runtime::R
    parsed::P
    reduced_output_indices::Vector{Int}
    source_signal_plan_indices::Vector{Int}
    evaluation_count::Int
    reset_count::Int
    ready::Bool
    execution_mode::Symbol
    reset_restorer::TimestepStateRestorer
end

abstract type AbstractEMTStudyCandidate end

"""
Typed setup-time extension point for repeated-study model parameters.

A model package adds a concrete subtype and an `apply_emt_parameter!` method
that mutates only the reset workspace before numerical execution. Candidate
vectors retain one concrete composite type, so optimizer and ensemble loops do
not introduce dynamic dispatch into the timestep solver.
"""
abstract type AbstractEMTModelParameter <: AbstractEMTStudyCandidate end

struct EMTSourceCrestCandidate{V<:AbstractVector{Float64}} <:
       AbstractEMTModelParameter
    crest_values::V
end

struct EMTSourceFrequencyRateParameter{V<:AbstractVector{Float64}} <:
       AbstractEMTModelParameter
    frequency_or_rate_values::V
end

struct EMTModelParameterCandidate{P<:Tuple} <: AbstractEMTStudyCandidate
    parameters::P

    function EMTModelParameterCandidate(parameters::P) where {P<:Tuple}
        all(parameter -> parameter isa AbstractEMTModelParameter, parameters) ||
            throw(ArgumentError(
                "EMT model-parameter candidates require typed EMT parameters",
            ))
        return new{P}(parameters)
    end
end

mutable struct DeckTimeSwitchStepWorkspace
    requested_closed_mask::Vector{Bool}
    previous_closed_mask::Vector{Bool}
    shifted_from_nodes::Vector{Int}
    shifted_to_nodes::Vector{Int}
    operation_queue::Vector{Int}
    source_voltage_differences::Vector{Float64}
    source_indices::Vector{Int}
    open_times::Vector{Float64}
    critical_currents::Vector{Float64}
    delay_times::Vector{Float64}
    node_group_successors::Vector{Int}
    grouped_from_nodes::Vector{Int}
    grouped_to_nodes::Vector{Int}
end

struct DeckStepConfigFeatures{O,BV,BC,BP,S,SC,ST,SA,SX} end

function EMTSourceCrestCandidate(values::AbstractVector{<:Real})
    crests = Float64.(values)
    all(isfinite, crests) ||
        throw(ArgumentError("source crest candidates must be finite"))
    return EMTSourceCrestCandidate(crests)
end

function EMTSourceFrequencyRateParameter(values::AbstractVector{<:Real})
    rates = Float64.(values)
    all(isfinite, rates) || throw(ArgumentError(
        "source frequency-or-rate parameters must be finite",
    ))
    return EMTSourceFrequencyRateParameter(rates)
end

EMTModelParameterCandidate(parameters::AbstractEMTModelParameter...) =
    EMTModelParameterCandidate(parameters)

abstract type AbstractEMTExecutionBackend end

struct EMTCPUBackend <: AbstractEMTExecutionBackend
    threaded::Bool
end

EMTCPUBackend(; threaded::Bool=false) = EMTCPUBackend(threaded)

abstract type AbstractEMTTraceReducer end

struct EMTOutputRMSReducer <: AbstractEMTTraceReducer
    channel_index::Int
end

struct EMTOutputPeakReducer <: AbstractEMTTraceReducer
    channel_index::Int
end

struct EMTStudyBatch{W,B<:AbstractEMTExecutionBackend}
    workspaces::Vector{W}
    backend::B
end

function ordered_node_names(node_map::Dict{Symbol,Int})
    names = Vector{Symbol}(undef, length(node_map))
    seen = falses(length(node_map))
    for (name, index) in node_map
        index > 0 || throw(ArgumentError("node indices must be positive"))
        index <= length(names) || throw(ArgumentError("node index $index exceeds node count"))
        !seen[index] || throw(ArgumentError("node index $index is assigned more than once"))
        seen[index] = true
        names[index] = name
    end
    return names
end

function normalized_node_map(node_map::AbstractDict{Symbol,<:Integer}, node_count::Int)
    length(node_map) == node_count ||
        throw(ArgumentError("node_map count must match system node_count"))
    normalized = Dict{Symbol,Int}()
    for (name, index) in node_map
        normalized[name] = Int(index)
    end
    ordered_node_names(normalized)
    return normalized
end

function fixed_step_count(dt_s::Float64, t_end_s::Float64)::Int
    dt_s > 0.0 || throw(ArgumentError("dt_s must be positive"))
    t_end_s >= 0.0 || throw(ArgumentError("t_end_s must be non-negative"))
    steps = Int(round(t_end_s / dt_s))
    abs(steps * dt_s - t_end_s) <= max(1.0e-15, 16.0 * eps(Float64) * max(abs(t_end_s), dt_s)) ||
        throw(ArgumentError("t_end_s must be an integer multiple of dt_s for this fixed-step deck runner"))
    return steps
end

function deck_fixed_step_horizon(parsed::DeckParser.DeckParseResult)
    options = DeckParser.deck_fixed_time_horizon_options(parsed)
    dt_s = Float64(options.dt_s)
    requested_t_end_s = Float64(options.tmax_s)
    dt_s > 0.0 || throw(ArgumentError("deck time step must be positive"))
    requested_t_end_s >= 0.0 ||
        throw(ArgumentError("deck stop time must be non-negative"))
    step_count = Int(round(requested_t_end_s / dt_s))
    step_count >= 0 || throw(ArgumentError("deck fixed-step count must be non-negative"))
    return (
        dt_s = dt_s,
        requested_t_end_s = requested_t_end_s,
        step_count = step_count,
        t_end_s = step_count * dt_s,
    )
end

function _deck_universal_machine_definition(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
    card_index::Int,
)
    row = findfirst(
        candidate -> candidate.machine_index == machine_index &&
                     candidate.card_index == card_index,
        DeckParser.deck_universal_machine_definition_rows(parsed),
    )
    row === nothing &&
        throw(ArgumentError("universal-machine definition card $card_index is missing"))
    return DeckParser.deck_universal_machine_definition_rows(parsed)[row]
end

function _deck_universal_machine_initialization_mode(parsed::DeckParser.DeckParseResult)
    return _deck_universal_machine_section(parsed).initialization_mode
end

function _deck_universal_machine_section(parsed::DeckParser.DeckParseResult)
    sections = DeckParser.deck_universal_machine_section_rows(parsed)
    length(sections) == 1 ||
        throw(ArgumentError("coupled d/q machine deck requires one universal-machine class-1 card"))
    return only(sections)
end

function deck_coupled_dq_machine_parameters(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
)
    DeckParser.assert_deck_valid!(parsed)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
    card3 = _deck_universal_machine_definition(parsed, machine_index, 3)
    card1.machine_type in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck coupled d/q machine runtime requires a type between 1 and 12"))
    coils = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    d_axis_coil_count, q_axis_coil_count =
        _deck_coupled_dq_axis_coil_counts(card1)
    coil_count = _deck_coupled_dq_coil_count(card1)
    length(coils) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) coupled d/q machine deck must define $coil_count coils"))
    card1.pole_pair_count > 0 ||
        throw(ArgumentError("induction-machine pole-pair count is missing"))
    card1.speed_convergence_margin === missing &&
        throw(ArgumentError("induction-machine speed convergence margin is missing"))
    card2.value2 === missing &&
        throw(ArgumentError("induction-machine d-axis inductance is missing"))
    card3.value2 === missing &&
        throw(ArgumentError("induction-machine q-axis inductance is missing"))
    section = _deck_universal_machine_section(parsed)
    angular_frequency_rad_s = 2.0 * pi * _deck_steady_state_frequency_hz(parsed)
    inductance_scale =
        section.parameter_basis == :power_frequency_normalized ?
        inv(angular_frequency_rad_s) : 1.0
    remanent_flux_scale = section.remanent_flux_enabled ? inductance_scale : 0.0
    raw_speed_tolerance = Float64(card1.speed_convergence_margin)
    speed_tolerance_rad_s = if raw_speed_tolerance == 0.0
        section.parameter_basis == :power_frequency_normalized ?
        0.001 * angular_frequency_rad_s :
        0.001 * angular_frequency_rad_s / card1.pole_pair_count
    elseif section.parameter_basis == :power_frequency_normalized
        raw_speed_tolerance * angular_frequency_rad_s
    else
        raw_speed_tolerance
    end
    return InductionMachineParameters(
        Float64[row.resistance for row in coils],
        Float64[row.inductance * inductance_scale for row in coils];
        time_step_s = time_step_s,
        d_axis_unsaturated_inductance = Float64(card2.value2) * inductance_scale,
        q_axis_unsaturated_inductance = Float64(card3.value2) * inductance_scale,
        d_axis_saturated_inductance =
            card2.saturated_inductance === missing ? 0.0 :
            Float64(card2.saturated_inductance) * inductance_scale,
        q_axis_saturated_inductance =
            card3.saturated_inductance === missing ? 0.0 :
            Float64(card3.saturated_inductance) * inductance_scale,
        d_axis_saturation_flux =
            card2.saturation_flux === missing ? 0.0 :
            Float64(card2.saturation_flux) * inductance_scale,
        q_axis_saturation_flux =
            card3.saturation_flux === missing ? 0.0 :
            Float64(card3.saturation_flux) * inductance_scale,
        d_axis_remanent_flux =
            card2.remanent_flux === missing ? 0.0 :
            Float64(card2.remanent_flux) * remanent_flux_scale,
        q_axis_remanent_flux =
            card3.remanent_flux === missing ? 0.0 :
            Float64(card3.remanent_flux) * remanent_flux_scale,
        d_axis_saturation_mode = card2.saturation_mode,
        q_axis_saturation_mode = card3.saturation_mode,
        pole_pair_count = card1.pole_pair_count,
        d_axis_coil_count = d_axis_coil_count,
        q_axis_coil_count = q_axis_coil_count,
        synchronous_electrical_speed_rad_s =
            angular_frequency_rad_s,
        speed_tolerance_rad_s = speed_tolerance_rad_s,
        machine_type = card1.machine_type,
        power_leakage_owned_by_network =
            section.initialization_mode == :automatic,
    )
end

function deck_coupled_dq_machine_state(
    parsed::DeckParser.DeckParseResult,
    initial_history_currents::AbstractVector{<:Real};
    machine_index::Int=1,
    power_frequency_hz::Real=_deck_steady_state_frequency_hz(parsed),
    initial_current_values::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck coupled d/q machine runtime requires a type between 1 and 12"))
    coil_count = _deck_coupled_dq_coil_count(card1)
    length(initial_history_currents) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) coupled d/q initial history requires $coil_count currents"))
    card1.pole_pair_count > 0 ||
        throw(ArgumentError("coupled d/q machine pole-pair count is missing"))
    histories = Float64.(initial_history_currents)
    if card1.machine_type in (1, 2, 3, 4, 5, 6, 7) &&
       _deck_universal_machine_initialization_mode(parsed) == :manual
        card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
        card3 = _deck_universal_machine_definition(parsed, machine_index, 3)
        card2.value1 === missing &&
            throw(ArgumentError("type-$(card1.machine_type) manual initialization requires mechanical speed on machine card 2"))
        card3.value1 === missing &&
            throw(ArgumentError("type-$(card1.machine_type) manual initialization requires rotor angle on machine card 3"))
        section = _deck_universal_machine_section(parsed)
        mechanical_speed_rad_s = Float64(card2.value1)
        if section.parameter_basis == :power_frequency_normalized
            mechanical_speed_rad_s *= 2.0 * pi * Float64(power_frequency_hz)
        end
        mechanical_angle_rad = card1.machine_type in (1, 2) ?
            (Float64(card3.value1) + pi / 2.0) / card1.pole_pair_count :
            Float64(card3.value1)
        transform = coupled_dq_power_transform(
            card1.pole_pair_count * mechanical_angle_rad;
            machine_type = card1.machine_type,
        )
        internal_histories = vcat(transform * histories[1:3], histories[4:end])
        currents = initial_current_values === nothing ?
            vcat(transform * histories[1:3], -histories[4:end]) :
            Float64.(initial_current_values)
        return InductionMachineState(
            currents,
            internal_histories;
            mechanical_speed_rad_s = mechanical_speed_rad_s,
            previous_mechanical_speed_rad_s = mechanical_speed_rad_s,
            mechanical_angle_rad = mechanical_angle_rad,
        )
    end
    if card1.machine_type in (8, 9, 10, 11, 12)
        _deck_universal_machine_initialization_mode(parsed) == :manual ||
            throw(ArgumentError("type-$(card1.machine_type) direct-current initialization requires explicit manual state ownership"))
        card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
        card3 = _deck_universal_machine_definition(parsed, machine_index, 3)
        card2.value1 === missing &&
            throw(ArgumentError("type-$(card1.machine_type) manual initialization requires mechanical speed on machine card 2"))
        card3.value1 === missing &&
            throw(ArgumentError("type-$(card1.machine_type) manual initialization requires mechanical angle on machine card 3"))
        transform = coupled_dq_power_transform(
            card1.pole_pair_count * Float64(card3.value1);
            machine_type = card1.machine_type,
        )
        internal_histories = vcat(transform * histories[1:3], histories[4:end])
        currents = initial_current_values === nothing ?
            vcat(transform * histories[1:3], -histories[4:end]) :
            Float64.(initial_current_values)
        return InductionMachineState(
            currents,
            internal_histories;
            mechanical_speed_rad_s = Float64(card2.value1),
            previous_mechanical_speed_rad_s = Float64(card2.value1),
            mechanical_angle_rad = Float64(card3.value1),
        )
    end
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    slip_percent = card1.machine_type in (1, 2) ? 0.0 :
        (card4.value1 === missing ? 0.0 : Float64(card4.value1))
    initial_angle_deg = card4.value2 === missing ? 0.0 : Float64(card4.value2)
    synchronous_electrical_speed = 2.0 * pi * Float64(power_frequency_hz)
    mechanical_speed =
        (1.0 - slip_percent / 100.0) * synchronous_electrical_speed /
        card1.pole_pair_count
    mechanical_angle =
        (deg2rad(initial_angle_deg) + pi / 2.0) / card1.pole_pair_count
    currents = initial_current_values === nothing ?
        vcat(histories[1:3], -histories[4:coil_count]) :
        Float64.(initial_current_values)
    return InductionMachineState(
        currents,
        histories;
        mechanical_speed_rad_s = mechanical_speed,
        previous_mechanical_speed_rad_s = mechanical_speed,
        mechanical_angle_rad = mechanical_angle,
    )
end

function _synchronous_field_current_for_d_axis_flux(
    parameters::InductionMachineParameters,
    power_currents::Vector{Float64},
    target_d_axis_flux::Float64,
    initial_field_current::Float64,
)
    length(power_currents) == 3 || throw(ArgumentError(
        "synchronous field initialization requires three power currents",
    ))
    all(isfinite, (target_d_axis_flux, initial_field_current)) || throw(ArgumentError(
        "synchronous field initialization inputs must be finite",
    ))
    function flux_residual(field_current::Float64)
        currents = vcat(power_currents, field_current, 0.0)
        d_axis_flux, _ = induction_machine_axis_fluxes(parameters, currents)
        return d_axis_flux - target_d_axis_flux
    end

    initial_residual = flux_residual(initial_field_current)
    tolerance = 1.0e-12 * max(1.0, abs(target_d_axis_flux))
    abs(initial_residual) <= tolerance && return initial_field_current
    interval = max(1.0, abs(initial_field_current))
    lower = initial_field_current - interval
    upper = initial_field_current + interval
    lower_residual = flux_residual(lower)
    upper_residual = flux_residual(upper)
    for _ in 1:64
        signbit(lower_residual) != signbit(upper_residual) && break
        interval *= 2.0
        lower = initial_field_current - interval
        upper = initial_field_current + interval
        lower_residual = flux_residual(lower)
        upper_residual = flux_residual(upper)
    end
    signbit(lower_residual) != signbit(upper_residual) || throw(ArgumentError(
        "synchronous field initialization could not bracket the requested d-axis flux",
    ))
    for _ in 1:80
        midpoint = (lower + upper) / 2.0
        midpoint_residual = flux_residual(midpoint)
        abs(midpoint_residual) <= tolerance && return midpoint
        if signbit(midpoint_residual) == signbit(lower_residual)
            lower = midpoint
            lower_residual = midpoint_residual
        else
            upper = midpoint
            upper_residual = midpoint_residual
        end
    end
    return (lower + upper) / 2.0
end

function _state_with_induction_machine_axis_fluxes!(
    state::InductionMachineState,
    parameters::InductionMachineParameters,
)
    state.d_axis_flux, state.q_axis_flux =
        induction_machine_axis_fluxes(parameters, state.current_values)
    return state
end

function deck_coupled_dq_machine_initial_state(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
    steady_state=deck_steady_state_voltage_phasors(parsed),
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck coupled d/q machine initialization requires a type between 1 and 12"))
    d_axis_coil_count, q_axis_coil_count =
        _deck_coupled_dq_axis_coil_counts(card1)
    coil_count = _deck_coupled_dq_coil_count(card1)
    if card1.machine_type in (8, 9, 10, 11, 12) &&
       _deck_universal_machine_initialization_mode(parsed) == :automatic
        return _deck_direct_machine_automatic_initialization(
            parsed,
            machine_index,
        ).state
    end
    if card1.machine_type in (8, 9, 10, 11, 12)
        coils = sort!(
            [
                row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
                if row.machine_index == machine_index
            ];
            by = row -> row.coil_index,
        )
        length(coils) == coil_count ||
            throw(ArgumentError("type-$(card1.machine_type) manual initialization requires $coil_count coil rows"))
        state = deck_coupled_dq_machine_state(
            parsed,
            Float64[row.initial_history_current for row in coils];
            machine_index = machine_index,
        )
        parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
        return _state_with_induction_machine_axis_fluxes!(state, parameters)
    end
    coils = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    length(coils) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) initialization requires $coil_count coil rows"))
    kinematics = deck_coupled_dq_machine_state(
        parsed,
        _deck_universal_machine_initialization_mode(parsed) == :manual ?
            Float64[row.initial_history_current for row in coils] :
            zeros(coil_count);
        machine_index = machine_index,
    )
    if _deck_universal_machine_initialization_mode(parsed) == :manual
        parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
        return _state_with_induction_machine_axis_fluxes!(kinematics, parameters)
    end
    terminal_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_terminal_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.terminal_index,
    )
    length(terminal_rows) == coil_count ||
        throw(ArgumentError("type-$(card1.machine_type) induction-machine initialization requires $coil_count terminals"))
    terminal_phasors = ComplexF64[
        row.terminal_node_value == 0 ? 0.0 :
        steady_state.node_voltage_phasors[row.terminal_node_value]
        for row in terminal_rows[1:3]
    ]
    if card1.machine_type in (1, 2)
        parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
        electrical_speed = 2.0 * pi * _deck_steady_state_frequency_hz(parsed)
        if !hasproperty(steady_state, :node_current_phasors)
            power_transform = coupled_dq_power_transform(
                card1.pole_pair_count * kinematics.mechanical_angle_rad,
                machine_type = card1.machine_type,
            )
            transformed_voltages = power_transform * real.(terminal_phasors)
            d_axis_flux = transformed_voltages[3] / electrical_speed
            power_currents = zeros(Float64, 3)
            unsaturated_field_current =
                d_axis_flux / parameters.d_axis_unsaturated_inductance
            field_current = _synchronous_field_current_for_d_axis_flux(
                parameters,
                power_currents,
                d_axis_flux,
                unsaturated_field_current,
            )
            currents = vcat(power_currents, field_current, 0.0)
            histories = vcat(power_currents, -field_current, -0.0)
            state = InductionMachineState(
                currents,
                histories;
                mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
                previous_mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
                mechanical_angle_rad = kinematics.mechanical_angle_rad,
            )
            return _state_with_induction_machine_axis_fluxes!(state, parameters)
        end
        terminal_current_phasors = ComplexF64[
            steady_state.node_current_phasors[row.terminal_node_value]
            for row in terminal_rows[1:3]
        ]
        sequence_voltage, sequence_current, axis_voltage_scale =
            if card1.machine_type == 1
                phase_rotation = cis(2.0 * pi / 3.0)
                (
                    (
                        terminal_phasors[1] +
                        phase_rotation * terminal_phasors[2] +
                        phase_rotation^2 * terminal_phasors[3]
                    ) / 3.0,
                    (
                        terminal_current_phasors[1] +
                        phase_rotation * terminal_current_phasors[2] +
                        phase_rotation^2 * terminal_current_phasors[3]
                    ) / 3.0,
                    sqrt(3.0 / 2.0),
                )
            else
                (
                    (terminal_phasors[2] - im * terminal_phasors[3]) / 2.0,
                    (terminal_current_phasors[2] - im * terminal_current_phasors[3]) / 2.0,
                    1.0,
                )
            end
        power_resistance = iszero(parameters.coil_conductances[3]) ?
            0.0 : inv(parameters.coil_conductances[3])
        quadrature_series_inductance =
            parameters.q_axis_unsaturated_inductance +
            parameters.coil_reactances[3]
        quadrature_internal_voltage = sequence_voltage +
            complex(power_resistance, electrical_speed * quadrature_series_inductance) *
            sequence_current
        rotor_electrical_angle = angle(quadrature_internal_voltage)
        mechanical_angle =
            (rotor_electrical_angle + pi / 2.0) / card1.pole_pair_count
        power_transform = coupled_dq_power_transform(
            card1.pole_pair_count * mechanical_angle,
            machine_type = card1.machine_type,
        )
        power_currents = power_transform * real.(terminal_current_phasors)
        direct_axis_leakage = parameters.coil_reactances[2]
        quadrature_axis_leakage = parameters.coil_reactances[3]
        saliency_current_flux = (
            parameters.d_axis_unsaturated_inductance -
            parameters.q_axis_unsaturated_inductance +
            direct_axis_leakage - quadrature_axis_leakage
        ) * power_currents[2]
        internal_voltage_axis_magnitude =
            axis_voltage_scale * abs(quadrature_internal_voltage)
        unsaturated_field_current = (
            internal_voltage_axis_magnitude / electrical_speed -
            saliency_current_flux
        ) / parameters.d_axis_unsaturated_inductance
        target_d_axis_flux = parameters.d_axis_unsaturated_inductance *
                             (power_currents[2] + unsaturated_field_current)
        field_current = _synchronous_field_current_for_d_axis_flux(
            parameters,
            power_currents,
            target_d_axis_flux,
            unsaturated_field_current,
        )
        q_axis_current = 0.0
        currents = vcat(power_currents, field_current, q_axis_current)
        histories = vcat(power_currents, -field_current, -q_axis_current)
        state = InductionMachineState(
            currents,
            histories;
            mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
            previous_mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
            mechanical_angle_rad = mechanical_angle,
        )
        return _state_with_induction_machine_axis_fluxes!(state, parameters)
    end
    if card1.machine_type in (6, 7)
        return _deck_single_phase_induction_initialization(
            parsed,
            machine_index,
        ).state
    end
    equivalent = _deck_induction_machine_steady_state_equivalent(
        parsed,
        machine_index,
    )
    return induction_machine_initial_state(
        terminal_phasors,
        equivalent;
        mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
        mechanical_angle_rad = kinematics.mechanical_angle_rad,
        pole_pair_count = card1.pole_pair_count,
        machine_type = card1.machine_type,
        d_axis_coil_count = d_axis_coil_count,
        q_axis_coil_count = q_axis_coil_count,
    )
end

function _deck_type8_automatic_initialization(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type == 8 ||
        throw(ArgumentError("automatic DC initialization requires machine type 8"))
    _deck_universal_machine_initialization_mode(parsed) == :automatic ||
        throw(ArgumentError("automatic type-8 initialization requires INITUM=1"))
    card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    card2.value1 === missing &&
        throw(ArgumentError("automatic type-8 initialization requires mechanical speed on machine card 2"))
    card4.value1 === missing &&
        throw(ArgumentError("automatic type-8 initialization requires requested armature voltage on machine card 4"))
    mechanical_speed = Float64(card2.value1)
    requested_armature_voltage = Float64(card4.value1)
    mechanical_speed != 0.0 ||
        throw(ArgumentError("automatic type-8 initialization requires nonzero mechanical speed"))

    parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
    coil_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    length(coil_rows) == 4 ||
        throw(ArgumentError("automatic type-8 initialization requires four coil rows"))
    armature_resistance = Float64(coil_rows[3].resistance)
    field_resistance = Float64(coil_rows[4].resistance)
    armature_resistance > 0.0 && field_resistance > 0.0 ||
        throw(ArgumentError("automatic type-8 initialization requires positive armature and field resistance"))

    terminal_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_terminal_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.terminal_index,
    )
    length(terminal_rows) == 4 ||
        throw(ArgumentError("automatic type-8 initialization requires four machine terminals"))
    power_node = Int(terminal_rows[3].terminal_node_value)
    field_node = Int(terminal_rows[4].terminal_node_value)
    power_thevenin = _deck_external_steady_state_thevenin(parsed, power_node)
    power_resistance = real(power_thevenin.impedance)
    power_resistance > 0.0 ||
        throw(ArgumentError("automatic type-8 armature requires a positive-real external Thevenin impedance"))
    physical_armature_current =
        (real(power_thevenin.voltage_phasor) - requested_armature_voltage) /
        power_resistance

    network_nodes = _deck_universal_machine_network_nodes(parsed, machine_index)
    field_source = network_nodes.field_source
    field_source === nothing &&
        throw(ArgumentError("automatic type-8 initialization requires a field source"))
    field_voltage = real(_deck_source_voltage_phasor(field_source))
    field_current = field_voltage / field_resistance
    d_axis_flux = parameters.d_axis_unsaturated_inductance * field_current
    q_axis_current = -physical_armature_current
    q_axis_flux = parameters.q_axis_unsaturated_inductance * q_axis_current
    electrical_speed = parameters.pole_pair_count * mechanical_speed
    back_emf = electrical_speed * d_axis_flux
    armature_kvl_residual =
        requested_armature_voltage -
        (back_emf + armature_resistance * physical_armature_current)
    kvl_tolerance = 1.0e-8 + 1.0e-8 * abs(requested_armature_voltage)
    abs(armature_kvl_residual) <= kvl_tolerance ||
        throw(ArgumentError(
            "automatic type-8 inputs violate steady armature KVL by $armature_kvl_residual V",
        ))

    currents = [0.0, 0.0, q_axis_current, field_current]
    histories = [0.0, 0.0, q_axis_current, -field_current]
    mechanical_angle = pi / (4.0 * parameters.pole_pair_count)
    state = InductionMachineState(
        currents,
        histories;
        mechanical_speed_rad_s = mechanical_speed,
        previous_mechanical_speed_rad_s = mechanical_speed,
        mechanical_angle_rad = mechanical_angle,
    )
    state.d_axis_flux = d_axis_flux
    state.q_axis_flux = q_axis_flux
    state.generated_torque =
        parameters.pole_pair_count * d_axis_flux * q_axis_current
    steady_state = _deck_external_steady_state_voltage_phasors(
        parsed,
        Dict(
            power_node => complex(q_axis_current, 0.0),
            field_node => complex(-field_current, 0.0),
        ),
    )
    return (
        source = :automatic_separately_excited_dc_machine_initialization,
        state = state,
        steady_state = steady_state,
        armature_node = power_node,
        field_node = field_node,
        requested_armature_voltage = requested_armature_voltage,
        armature_thevenin_voltage = real(power_thevenin.voltage_phasor),
        armature_thevenin_resistance = power_resistance,
        armature_resistance = armature_resistance,
        physical_armature_current = physical_armature_current,
        field_voltage = field_voltage,
        field_resistance = field_resistance,
        field_current = field_current,
        electrical_speed_rad_s = electrical_speed,
        back_emf = back_emf,
        armature_kvl_residual = armature_kvl_residual,
        field_kvl_residual = field_voltage - field_resistance * field_current,
        defect_isolation = :compiled_type8_label8800_zero_state_excluded,
    )
end

function _direct_machine_steady_state_step(
    parameters::InductionMachineParameters,
    current_values::Vector{Float64},
    power_terminal_open_circuit_voltage::Float64,
    mechanical_speed_rad_s::Float64,
    mechanical_angle_rad::Float64,
    power_terminal_thevenin_impedance::Float64,
)
    histories = vcat(current_values[1:3], -current_values[4:end])
    state = InductionMachineState(
        current_values,
        histories;
        mechanical_speed_rad_s,
        previous_mechanical_speed_rad_s = mechanical_speed_rad_s,
        mechanical_angle_rad,
    )
    state.d_axis_flux, state.q_axis_flux =
        induction_machine_axis_fluxes(parameters, current_values)
    power_terminal_voltages =
        [0.0, 0.0, power_terminal_open_circuit_voltage]
    rotor_thevenin_matrix = zeros(3, 3)
    rotor_thevenin_matrix[3, 3] =
        power_terminal_thevenin_impedance
    stator_count = length(current_values) - 3
    stator_terminal_voltages = zeros(stator_count)
    stator_thevenin_matrix = zeros(stator_count, stator_count)
    coupled_dq_machine_step!(
        state,
        parameters;
        power_terminal_voltages,
        rotor_thevenin_matrix,
        mechanical_speed_thevenin_rad_s = mechanical_speed_rad_s,
        generated_torque_impedance = 0.0,
        stator_terminal_voltages,
        stator_thevenin_matrix,
        initial_step = true,
    )
    result = coupled_dq_machine_step!(
        state,
        parameters;
        power_terminal_voltages,
        rotor_thevenin_matrix,
        mechanical_speed_thevenin_rad_s = mechanical_speed_rad_s,
        generated_torque_impedance = 0.0,
        stator_terminal_voltages,
        stator_thevenin_matrix,
        initial_step = false,
    )
    return result
end

function _direct_machine_steady_state_current_seed(
    parameters::InductionMachineParameters,
    armature_voltage::Float64,
)
    machine_type = parameters.machine_type
    machine_type in (9, 10, 11, 12) ||
        throw(ArgumentError("direct-machine fixed-point seed requires type 9 through 12"))
    currents = zeros(5)
    armature_resistance = parameters.coil_conductances[3] == 0.0 ?
        0.0 : inv(parameters.coil_conductances[3])
    series_resistance = parameters.coil_conductances[5] == 0.0 ?
        0.0 : inv(parameters.coil_conductances[5])
    resistance_scale = max(armature_resistance + series_resistance, eps(Float64))
    currents[3] = armature_voltage / resistance_scale
    currents[4] = parameters.coil_conductances[4] == 0.0 ?
        0.0 : armature_voltage * parameters.coil_conductances[4]
    currents[5] = machine_type == 12 ? 0.0 :
        machine_type == 11 ? currents[3] - currents[4] : currents[3]
    return currents
end

function _direct_machine_runtime_fixed_point(
    parameters::InductionMachineParameters,
    power_terminal_open_circuit_voltage::Float64,
    mechanical_speed_rad_s::Float64;
    power_terminal_thevenin_impedance::Real = 0.0,
    initial_current_values::Union{Nothing,AbstractVector{<:Real}} = nothing,
    maximum_iterations::Int = 24,
    maximum_history_iterations::Int = 200,
    absolute_tolerance::Float64 = 1.0e-10,
    relative_tolerance::Float64 = 1.0e-10,
)
    parameters.machine_type in (9, 10, 11, 12) ||
        throw(ArgumentError("direct-machine runtime fixed point requires type 9 through 12"))
    maximum_iterations > 0 ||
        throw(ArgumentError("direct-machine fixed-point iteration count must be positive"))
    maximum_history_iterations > 0 ||
        throw(ArgumentError("direct-machine history-settling iteration count must be positive"))
    all(
        value -> isfinite(value) && value >= 0.0,
        (absolute_tolerance, relative_tolerance),
    ) || throw(ArgumentError("direct-machine fixed-point tolerances must be finite and nonnegative"))
    thevenin_impedance = Float64(power_terminal_thevenin_impedance)
    isfinite(power_terminal_open_circuit_voltage) &&
        isfinite(thevenin_impedance) && thevenin_impedance >= 0.0 ||
        throw(ArgumentError("direct-machine fixed-point electrical boundary must be finite with nonnegative impedance"))
    mechanical_angle = pi / (2.0 * parameters.pole_pair_count)
    currents = initial_current_values === nothing ?
        _direct_machine_steady_state_current_seed(
            parameters,
            power_terminal_open_circuit_voltage,
        ) :
        Float64.(initial_current_values)
    length(currents) == length(parameters.coil_conductances) ||
        throw(ArgumentError("direct-machine fixed-point current seed has the wrong coil count"))
    all(isfinite, currents) ||
        throw(ArgumentError("direct-machine fixed-point current seed must be finite"))
    active = 3:5
    residual = zeros(3)
    completed_iterations = 0
    converged = false
    final_result = nothing
    for iteration in 1:maximum_iterations
        completed_iterations = iteration
        result = _direct_machine_steady_state_step(
            parameters,
            currents,
            power_terminal_open_circuit_voltage,
            mechanical_speed_rad_s,
            mechanical_angle,
            thevenin_impedance,
        )
        final_result = result
        residual .= result.current_values[active] .- currents[active]
        residual_norm = maximum(abs, residual; init = 0.0)
        current_scale = maximum(abs, result.current_values[active]; init = 0.0)
        allowance = absolute_tolerance + relative_tolerance * max(current_scale, 1.0)
        if residual_norm <= allowance
            currents .= result.current_values
            converged = true
            break
        end

        jacobian = zeros(3, 3)
        for column in eachindex(residual)
            current_index = first(active) + column - 1
            perturbation =
                1.0e-4 * max(abs(currents[current_index]), 1.0)
            upper = copy(currents)
            lower = copy(currents)
            upper[current_index] += perturbation
            lower[current_index] -= perturbation
            upper_result = _direct_machine_steady_state_step(
                parameters,
                upper,
                power_terminal_open_circuit_voltage,
                mechanical_speed_rad_s,
                mechanical_angle,
                thevenin_impedance,
            )
            lower_result = _direct_machine_steady_state_step(
                parameters,
                lower,
                power_terminal_open_circuit_voltage,
                mechanical_speed_rad_s,
                mechanical_angle,
                thevenin_impedance,
            )
            upper_residual =
                upper_result.current_values[active] .- upper[active]
            lower_residual =
                lower_result.current_values[active] .- lower[active]
            jacobian[:, column] .=
                (upper_residual .- lower_residual) ./ (2.0 * perturbation)
        end
        correction = try
            -(jacobian \ residual)
        catch error
            error isa SingularException || rethrow()
            -(pinv(jacobian) * residual)
        end
        all(isfinite, correction) ||
            throw(ArgumentError("direct-machine fixed-point correction is nonfinite"))
        accepted = false
        best_currents = currents
        best_residual_norm = residual_norm
        for direction in (correction, copy(residual))
            for damping in (
                1.0,
                0.5,
                0.25,
                0.125,
                0.0625,
                0.03125,
                0.015625,
                0.0078125,
                0.00390625,
            )
                candidate = copy(currents)
                candidate[active] .+= damping .* direction
                candidate_result = _direct_machine_steady_state_step(
                    parameters,
                    candidate,
                    power_terminal_open_circuit_voltage,
                    mechanical_speed_rad_s,
                    mechanical_angle,
                    thevenin_impedance,
                )
                candidate_residual =
                    candidate_result.current_values[active] .- candidate[active]
                candidate_residual_norm =
                    maximum(abs, candidate_residual; init = 0.0)
                if candidate_residual_norm < best_residual_norm
                    best_currents = candidate
                    best_residual_norm = candidate_residual_norm
                    accepted = true
                end
            end
        end
        if accepted
            currents = best_currents
        elseif residual_norm <= 10.0 * allowance
            currents .= result.current_values
            converged = true
            break
        else
            throw(ArgumentError(
                "direct-machine fixed-point correction did not reduce runtime residual $residual_norm",
            ))
        end
    end
    converged || throw(ArgumentError(
        "direct-machine fixed-point initialization did not converge after $maximum_iterations iterations",
    ))
    state = InductionMachineState(
        currents,
        vcat(currents[1:3], -currents[4:end]);
        mechanical_speed_rad_s,
        previous_mechanical_speed_rad_s = mechanical_speed_rad_s,
        mechanical_angle_rad = mechanical_angle,
    )
    power_terminal_voltages =
        [0.0, 0.0, power_terminal_open_circuit_voltage]
    rotor_thevenin_matrix = zeros(Float64, 3, 3)
    rotor_thevenin_matrix[3, 3] = thevenin_impedance
    stator_count = length(currents) - 3
    stator_terminal_voltages = zeros(Float64, stator_count)
    stator_thevenin_matrix = zeros(Float64, stator_count, stator_count)
    coupled_dq_machine_step!(
        state,
        parameters;
        power_terminal_voltages,
        rotor_thevenin_matrix,
        mechanical_speed_thevenin_rad_s = mechanical_speed_rad_s,
        generated_torque_impedance = 0.0,
        stator_terminal_voltages,
        stator_thevenin_matrix,
        initial_step = true,
    )
    fixed_current_prefix = copy(state.current_values[1:2])
    fixed_history_prefix = copy(state.history_currents[1:2])
    variables = vcat(
        state.current_values[active],
        state.history_currents[active],
    )

    function evaluate_runtime_state(candidate_variables::Vector{Float64})
        candidate_currents = vcat(
            fixed_current_prefix,
            candidate_variables[1:3],
        )
        candidate_histories = vcat(
            fixed_history_prefix,
            candidate_variables[4:6],
        )
        candidate_state = InductionMachineState(
            candidate_currents,
            candidate_histories;
            mechanical_speed_rad_s,
            previous_mechanical_speed_rad_s =
                mechanical_speed_rad_s,
            mechanical_angle_rad = mechanical_angle,
        )
        candidate_state.d_axis_flux, candidate_state.q_axis_flux =
            induction_machine_axis_fluxes(
                parameters,
                candidate_currents,
            )
        candidate_result = coupled_dq_machine_step!(
            candidate_state,
            parameters;
            power_terminal_voltages,
            rotor_thevenin_matrix,
            mechanical_speed_thevenin_rad_s =
                mechanical_speed_rad_s,
            generated_torque_impedance = 0.0,
            stator_terminal_voltages,
            stator_thevenin_matrix,
            initial_step = false,
        )
        mapped_variables = vcat(
            candidate_result.current_values[active],
            candidate_result.history_currents[active],
        )
        return (
            residual = mapped_variables - candidate_variables,
            state = candidate_state,
            result = candidate_result,
        )
    end

    history_iterations = 0
    history_converged = false
    final_result = nothing
    final_current_residual = zeros(Float64, length(currents))
    final_history_residual = zeros(Float64, length(currents))
    for iteration in 1:maximum_history_iterations
        history_iterations = iteration
        evaluated = evaluate_runtime_state(variables)
        residual = evaluated.residual
        residual_norm = maximum(abs, residual; init = 0.0)
        current_residual_norm =
            maximum(abs, residual[1:3]; init = 0.0)
        history_residual_norm =
            maximum(abs, residual[4:6]; init = 0.0)
        current_allowance = absolute_tolerance + relative_tolerance *
            max(maximum(abs, variables[1:3]; init = 0.0), 1.0)
        history_allowance = absolute_tolerance + relative_tolerance *
            max(maximum(abs, variables[4:6]; init = 0.0), 1.0)
        if current_residual_norm <= current_allowance &&
           history_residual_norm <= history_allowance
            state = evaluated.state
            final_result = evaluated.result
            final_current_residual[active] .= residual[1:3]
            final_history_residual[active] .= residual[4:6]
            history_converged = true
            break
        end
        jacobian = zeros(Float64, 6, 6)
        for column in 1:6
            perturbation =
                1.0e-5 * max(abs(variables[column]), 1.0)
            upper = copy(variables)
            lower = copy(variables)
            upper[column] += perturbation
            lower[column] -= perturbation
            upper_residual =
                evaluate_runtime_state(upper).residual
            lower_residual =
                evaluate_runtime_state(lower).residual
            jacobian[:, column] .= (
                upper_residual - lower_residual
            ) ./ (2.0 * perturbation)
        end
        correction = try
            -(jacobian \ residual)
        catch error
            error isa SingularException || rethrow()
            -(pinv(jacobian) * residual)
        end
        all(isfinite, correction) ||
            throw(ArgumentError("direct-machine history fixed-point correction is nonfinite"))
        accepted = false
        best_variables = variables
        best_residual_norm = residual_norm
        for direction in (correction, copy(residual))
            for damping in (
                1.0,
                0.5,
                0.25,
                0.125,
                0.0625,
                0.03125,
                0.015625,
                0.0078125,
            )
                candidate = variables + damping * direction
                candidate_residual =
                    evaluate_runtime_state(candidate).residual
                candidate_norm =
                    maximum(abs, candidate_residual; init = 0.0)
                if candidate_norm < best_residual_norm
                    best_variables = candidate
                    best_residual_norm = candidate_norm
                    accepted = true
                end
            end
        end
        if accepted
            variables = best_variables
        elseif current_residual_norm <= 20.0 * current_allowance &&
               history_residual_norm <= 20.0 * history_allowance
            state = evaluated.state
            final_result = evaluated.result
            final_current_residual[active] .= residual[1:3]
            final_history_residual[active] .= residual[4:6]
            history_converged = true
            break
        else
            throw(ArgumentError(
                "direct-machine current/history fixed-point correction did not reduce residual $residual_norm",
            ))
        end
    end
    history_converged || throw(ArgumentError(
        "direct-machine runtime current/history fixed point did not converge after $maximum_history_iterations iterations",
    ))
    state.mechanical_speed_rad_s = mechanical_speed_rad_s
    state.previous_mechanical_speed_rad_s = mechanical_speed_rad_s
    state.mechanical_angle_rad = mechanical_angle
    state.call_count = 0
    return (
        state,
        runtime_result = final_result,
        current_residual = collect(final_current_residual),
        maximum_current_residual =
            maximum(abs, final_current_residual; init = 0.0),
        history_residual = collect(final_history_residual),
        maximum_history_residual =
            maximum(abs, final_history_residual; init = 0.0),
        iteration_count = completed_iterations,
        history_iteration_count = history_iterations,
        converged = converged && history_converged,
    )
end

function _deck_type9_through_12_automatic_electrical_state(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type in (9, 10, 11, 12) ||
        throw(ArgumentError("automatic coupled-field DC initialization requires type 9 through 12"))
    _deck_universal_machine_initialization_mode(parsed) == :automatic ||
        throw(ArgumentError("automatic coupled-field DC initialization requires INITUM=1"))
    card2 = _deck_universal_machine_definition(parsed, machine_index, 2)
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    card2.value1 === missing &&
        throw(ArgumentError("automatic type-$(card1.machine_type) initialization requires mechanical speed on machine card 2"))
    card4.value1 === missing &&
        throw(ArgumentError("automatic type-$(card1.machine_type) initialization requires requested armature voltage on machine card 4"))
    mechanical_speed = Float64(card2.value1)
    requested_armature_voltage = Float64(card4.value1)
    mechanical_speed != 0.0 ||
        throw(ArgumentError("automatic type-$(card1.machine_type) initialization requires nonzero mechanical speed"))
    requested_armature_voltage > 0.0 ||
        throw(ArgumentError("automatic type-$(card1.machine_type) initialization requires positive armature voltage"))
    parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
    fixed_point = _direct_machine_runtime_fixed_point(
        parameters,
        requested_armature_voltage,
        mechanical_speed,
    )
    network_nodes = _deck_universal_machine_network_nodes(parsed, machine_index)
    power_node = Int(network_nodes.power[only(network_nodes.active_power_positions)])
    power_current = Float64(fixed_point.runtime_result.power_terminal_currents[3])
    return (
        source = :automatic_coupled_field_dc_machine_electrical_state,
        machine_index,
        machine_type = card1.machine_type,
        state = fixed_point.state,
        armature_node = power_node,
        requested_armature_voltage,
        armature_current_injection = power_current,
        current_injections = Dict(power_node => complex(power_current, 0.0)),
        fixed_point_iteration_count = fixed_point.iteration_count,
        fixed_point_current_residual = fixed_point.current_residual,
        fixed_point_maximum_current_residual =
            fixed_point.maximum_current_residual,
        fixed_point_history_iteration_count =
            fixed_point.history_iteration_count,
        fixed_point_history_residual =
            fixed_point.history_residual,
        fixed_point_maximum_history_residual =
            fixed_point.maximum_history_residual,
        fixed_point_converged = fixed_point.converged,
        compiled_defect_isolation =
            :compiled_direct_machine_label8800_uninitialized_state_excluded,
    )
end

function _deck_type9_through_12_automatic_initialization(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    group = _deck_type9_through_12_automatic_initialization_group(
        parsed,
        [machine_index],
    )
    electrical = only(group.electrical_states)
    solved_voltage = only(group.solved_armature_voltages)
    voltage_residual = only(group.armature_voltage_residuals)
    return merge(
        electrical,
        (
            source = :automatic_coupled_field_dc_machine_initialization,
            steady_state = group.steady_state,
            solved_armature_voltage = solved_voltage,
            armature_voltage_residual = voltage_residual,
            runtime_group_fixed_point_iteration_count =
                group.runtime_group_fixed_point_iteration_count,
            runtime_group_current_residual =
                group.runtime_group_current_residual,
            runtime_group_maximum_current_residual =
                group.runtime_group_maximum_current_residual,
        ),
    )
end

function _deck_direct_machine_group_runtime_fixed_points(
    parsed::DeckParser.DeckParseResult,
    electrical_states::AbstractVector;
    maximum_iterations::Int = 12,
    absolute_tolerance::Float64 = 1.0e-10,
    relative_tolerance::Float64 = 1.0e-10,
)
    machine_count = length(electrical_states)
    machine_count > 0 ||
        throw(ArgumentError("automatic direct-machine runtime group must not be empty"))
    dt_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s
    power_nodes = Int[
        electrical.armature_node for electrical in electrical_states
    ]
    length(unique(power_nodes)) == machine_count ||
        throw(ArgumentError("automatic direct machines must own distinct armature nodes"))
    probe = initialize_step_context(
        parsed;
        dt_s,
        t_end_s = dt_s,
        recorded_step_indices = [0, 1],
    )
    calibration = solve_step_with_compensated_current_injections!(
        probe.system,
        dt_s,
        dt_s,
        zeros(Float64, probe.system.node_count),
        power_nodes,
        (_voltage, _impedance) -> zeros(Float64, machine_count),
    )
    thevenin_impedance =
        Matrix{Float64}(calibration.compensation_impedance)
    maximum(abs, thevenin_impedance - transpose(thevenin_impedance); init = 0.0) <=
        1.0e-10 * max(maximum(abs, thevenin_impedance; init = 0.0), 1.0) ||
        throw(ArgumentError("automatic direct-machine runtime Thévenin matrix is not reciprocal"))
    all(diag(thevenin_impedance) .>= 0.0) ||
        throw(ArgumentError("automatic direct-machine runtime Thévenin matrix has negative driving impedance"))
    parameters = [
        deck_coupled_dq_machine_parameters(
            parsed;
            machine_index = electrical.machine_index,
            time_step_s = dt_s,
        ) for electrical in electrical_states
    ]
    requested_voltages = Float64[
        electrical.requested_armature_voltage
        for electrical in electrical_states
    ]
    mechanical_speeds = Float64[
        electrical.state.mechanical_speed_rad_s
        for electrical in electrical_states
    ]
    initial_current_values = [
        copy(electrical.state.current_values)
        for electrical in electrical_states
    ]

    function evaluate(guesses::Vector{Float64})
        open_circuit_voltages =
            requested_voltages - thevenin_impedance * guesses
        fixed_points = [
            _direct_machine_runtime_fixed_point(
                parameters[index],
                open_circuit_voltages[index],
                mechanical_speeds[index];
                power_terminal_thevenin_impedance =
                    thevenin_impedance[index, index],
                initial_current_values = initial_current_values[index],
                absolute_tolerance,
                relative_tolerance,
            ) for index in 1:machine_count
        ]
        mapped_currents = Float64[
            fixed_point.runtime_result.power_terminal_currents[3]
            for fixed_point in fixed_points
        ]
        return mapped_currents, fixed_points, open_circuit_voltages
    end

    currents = Float64[
        electrical.armature_current_injection
        for electrical in electrical_states
    ]
    converged = false
    completed_iterations = 0
    fixed_points = Any[]
    open_circuit_voltages = zeros(Float64, machine_count)
    residual = zeros(Float64, machine_count)
    for iteration in 1:maximum_iterations
        completed_iterations = iteration
        mapped, trial_fixed_points, trial_open_circuit =
            evaluate(currents)
        residual .= mapped .- currents
        residual_norm = maximum(abs, residual; init = 0.0)
        allowance = absolute_tolerance + relative_tolerance *
            max(maximum(abs, mapped; init = 0.0), 1.0)
        if residual_norm <= allowance
            currents .= mapped
            fixed_points = trial_fixed_points
            open_circuit_voltages .= trial_open_circuit
            converged = true
            break
        end
        jacobian = zeros(Float64, machine_count, machine_count)
        for column in 1:machine_count
            perturbation =
                1.0e-6 * max(abs(currents[column]), 1.0)
            perturbed = copy(currents)
            perturbed[column] += perturbation
            perturbed_mapped, _, _ = evaluate(perturbed)
            jacobian[:, column] .= (
                perturbed_mapped .-
                perturbed .-
                residual
            ) ./ perturbation
        end
        correction = try
            -(jacobian \ residual)
        catch error
            error isa SingularException || rethrow()
            -(pinv(jacobian) * residual)
        end
        all(isfinite, correction) ||
            throw(ArgumentError("automatic direct-machine group correction is nonfinite"))
        accepted = false
        for damping in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            candidate = currents + damping * correction
            candidate_mapped, _, _ = evaluate(candidate)
            candidate_residual = candidate_mapped - candidate
            if maximum(abs, candidate_residual; init = 0.0) < residual_norm
                currents = candidate
                accepted = true
                break
            end
        end
        accepted || throw(ArgumentError(
            "automatic direct-machine group correction did not reduce the runtime residual",
        ))
    end
    converged || throw(ArgumentError(
        "automatic direct-machine runtime group did not converge after $maximum_iterations iterations",
    ))
    updated_electrical_states = [
        merge(
            electrical_states[index],
            (
                state = fixed_points[index].state,
                armature_current_injection = currents[index],
                current_injections = Dict(
                    power_nodes[index] => complex(currents[index], 0.0),
                ),
                fixed_point_iteration_count =
                    fixed_points[index].iteration_count,
                fixed_point_current_residual =
                    fixed_points[index].current_residual,
                fixed_point_maximum_current_residual =
                    fixed_points[index].maximum_current_residual,
                fixed_point_history_iteration_count =
                    fixed_points[index].history_iteration_count,
                fixed_point_history_residual =
                    fixed_points[index].history_residual,
                fixed_point_maximum_history_residual =
                    fixed_points[index].maximum_history_residual,
                runtime_power_terminal_open_circuit_voltage =
                    open_circuit_voltages[index],
                runtime_power_terminal_thevenin_impedance =
                    thevenin_impedance[index, index],
            ),
        ) for index in 1:machine_count
    ]
    return (
        source = :automatic_direct_machine_group_runtime_fixed_points,
        electrical_states = updated_electrical_states,
        power_terminal_thevenin_impedance = thevenin_impedance,
        power_terminal_open_circuit_voltages = open_circuit_voltages,
        power_terminal_currents = currents,
        current_residual = residual,
        maximum_current_residual =
            maximum(abs, residual; init = 0.0),
        iteration_count = completed_iterations,
        converged,
    )
end

function _deck_type9_through_12_automatic_initialization_group(
    parsed::DeckParser.DeckParseResult,
    machine_indices::AbstractVector{<:Integer},
)
    indices = Int.(machine_indices)
    !isempty(indices) && length(unique(indices)) == length(indices) ||
        throw(ArgumentError("automatic direct-machine group indices must be nonempty and unique"))
    preliminary_electrical_states = [
        _deck_type9_through_12_automatic_electrical_state(
            parsed,
            machine_index,
        ) for machine_index in indices
    ]
    runtime_fixed_points =
        _deck_direct_machine_group_runtime_fixed_points(
            parsed,
            preliminary_electrical_states,
        )
    electrical_states = runtime_fixed_points.electrical_states
    current_injections = Dict{Int,ComplexF64}()
    for electrical in electrical_states
        for (node, injection) in electrical.current_injections
            current_injections[Int(node)] =
                get(current_injections, Int(node), 0.0 + 0.0im) +
                ComplexF64(injection)
        end
    end
    steady_state = _deck_external_steady_state_voltage_phasors(
        parsed,
        current_injections,
        excluded_universal_machine_indices = indices,
    )
    solved_armature_voltages = ComplexF64[
        steady_state.node_voltage_phasors[electrical.armature_node]
        for electrical in electrical_states
    ]
    armature_voltage_residuals = ComplexF64[
        solved_armature_voltages[index] -
        electrical_states[index].requested_armature_voltage
        for index in eachindex(electrical_states)
    ]
    for index in eachindex(electrical_states)
        electrical = electrical_states[index]
        residual = armature_voltage_residuals[index]
        tolerance =
            2.0e-7 + 1.0e-8 * abs(electrical.requested_armature_voltage)
        abs(residual) <= tolerance || throw(ArgumentError(
            "automatic type-$(electrical.machine_type) machine $(electrical.machine_index) armature/network steady state is inconsistent by $residual V",
        ))
    end
    return (
        source = :automatic_coupled_field_dc_machine_group_initialization,
        machine_indices = indices,
        electrical_states,
        current_injections,
        steady_state,
        solved_armature_voltages,
        armature_voltage_residuals,
        maximum_armature_voltage_residual =
            maximum(abs, armature_voltage_residuals; init = 0.0),
        runtime_power_terminal_thevenin_impedance =
            runtime_fixed_points.power_terminal_thevenin_impedance,
        runtime_power_terminal_open_circuit_voltages =
            runtime_fixed_points.power_terminal_open_circuit_voltages,
        runtime_group_fixed_point_iteration_count =
            runtime_fixed_points.iteration_count,
        runtime_group_current_residual =
            runtime_fixed_points.current_residual,
        runtime_group_maximum_current_residual =
            runtime_fixed_points.maximum_current_residual,
    )
end

function _deck_direct_machine_automatic_initialization(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type == 8 &&
        return _deck_type8_automatic_initialization(parsed, machine_index)
    card1.machine_type in (9, 10, 11, 12) &&
        return _deck_type9_through_12_automatic_initialization(
            parsed,
            machine_index,
        )
    throw(ArgumentError("automatic direct-machine initialization requires type 8 through 12"))
end

function deck_separately_excited_dc_initialization(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int=1,
)
    return _deck_type8_automatic_initialization(parsed, machine_index)
end

function deck_direct_current_machine_automatic_initialization(
    parsed::DeckParser.DeckParseResult;
    machine_index::Int = 1,
)
    return _deck_direct_machine_automatic_initialization(
        parsed,
        machine_index,
    )
end

function _deck_coupled_dq_axis_coil_counts(card1)
    if card1.machine_type in (3, 4, 5)
        return (1, 1)
    end
    return (card1.d_axis_coil_count, card1.q_axis_coil_count)
end

function _deck_coupled_dq_coil_count(card1)
    d_axis_coil_count, q_axis_coil_count =
        _deck_coupled_dq_axis_coil_counts(card1)
    return 3 + d_axis_coil_count + q_axis_coil_count +
           (card1.machine_type == 4 ? 1 : 0)
end

function _deck_single_phase_excitation_resistance_ohm(
    parsed::DeckParser.DeckParseResult,
    field_node::Int,
)
    shunt_conductance = 0.0
    for row in DeckParser.deck_over2_branch_rows(parsed)
        row.branch_kind == :conductance || continue
        from_node = Int(row.from_node_value)
        to_node = Int(row.to_node_value)
        if (from_node == field_node && to_node == 0) ||
           (to_node == field_node && from_node == 0)
            shunt_conductance += Float64(row.conductance)
        end
    end
    shunt_conductance > 0.0 ||
        throw(ArgumentError("single-phase excitation source requires a grounded resistance"))
    return inv(shunt_conductance)
end

function _deck_single_phase_excitation_steady_state_impedance(
    parsed::DeckParser.DeckParseResult,
    field_node::Int,
)
    physical_resistance_ohm =
        _deck_single_phase_excitation_resistance_ohm(parsed, field_node)
    excitation_resistance_base_ohm = 1.0e6
    return complex(physical_resistance_ohm / excitation_resistance_base_ohm, 0.0)
end

function _deck_single_phase_induction_initialization(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    card1.machine_type in (6, 7) ||
        throw(ArgumentError("single-phase induction initialization requires machine type 6 or 7"))
    parameters = deck_coupled_dq_machine_parameters(parsed; machine_index)
    expected_terminal_count = card1.machine_type == 6 ? 4 : 5
    kinematics = deck_coupled_dq_machine_state(
        parsed,
        zeros(expected_terminal_count);
        machine_index = machine_index,
    )
    terminal_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_terminal_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.terminal_index,
    )
    length(terminal_rows) == expected_terminal_count ||
        throw(ArgumentError("type-$(card1.machine_type) single-phase induction initialization requires $expected_terminal_count terminals"))
    power_terminal_node = Int(terminal_rows[3].terminal_node_value)
    power_thevenin = _deck_external_steady_state_thevenin(parsed, power_terminal_node)
    network_nodes = _deck_universal_machine_network_nodes(parsed, machine_index)
    field_source = network_nodes.field_source
    field_source === nothing &&
        throw(ArgumentError("single-phase induction initialization requires a field source"))
    coil_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    length(coil_rows) == expected_terminal_count ||
        throw(ArgumentError("type-$(card1.machine_type) single-phase induction initialization requires $expected_terminal_count coils"))
    card4 = _deck_universal_machine_definition(parsed, machine_index, 4)
    slip = card4.value1 === missing ? 0.0 : Float64(card4.value1) / 100.0
    initialization = single_phase_induction_steady_state_initialization(
        power_thevenin_voltage_phasor = power_thevenin.voltage_phasor,
        power_thevenin_impedance = power_thevenin.impedance,
        excitation_voltage_phasor = -_deck_source_voltage_phasor(field_source),
        excitation_thevenin_impedance =
            _deck_single_phase_excitation_steady_state_impedance(
                parsed,
                Int(network_nodes.field),
            ),
        q_axis_excitation_voltage_phasor = card1.machine_type == 7 ?
            -_deck_source_voltage_phasor(network_nodes.q_axis_field_source) : nothing,
        q_axis_excitation_thevenin_impedance = card1.machine_type == 7 ?
            _deck_single_phase_excitation_steady_state_impedance(
                parsed,
                Int(network_nodes.q_axis_field),
            ) : 0.0 + 0.0im,
        coil_resistances = Float64[row.resistance for row in coil_rows],
        coil_inductances = Float64[row.inductance for row in coil_rows],
        d_axis_main_inductance = parameters.d_axis_unsaturated_inductance,
        q_axis_main_inductance = parameters.q_axis_unsaturated_inductance,
        slip = slip,
        frequency_hz = _deck_steady_state_frequency_hz(parsed),
        mechanical_speed_rad_s = kinematics.mechanical_speed_rad_s,
        mechanical_angle_rad = kinematics.mechanical_angle_rad,
        pole_pair_count = card1.pole_pair_count,
    )
    steady_state = _deck_external_steady_state_voltage_phasors(
        parsed,
        Dict(power_terminal_node => initialization.power_terminal_current_phasor),
    )
    return merge(initialization, (; steady_state))
end

function deck_induction_machine_parameters(parsed::DeckParser.DeckParseResult; kwargs...)
    parameters = deck_coupled_dq_machine_parameters(parsed; kwargs...)
    parameters.machine_type in (3, 4, 5, 6, 7) ||
        throw(ArgumentError("deck_induction_machine_parameters requires machine type 3 through 7"))
    return parameters
end

function deck_induction_machine_state(parsed::DeckParser.DeckParseResult, histories; kwargs...)
    state = deck_coupled_dq_machine_state(parsed, histories; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type in (3, 4, 5, 6, 7) ||
        throw(ArgumentError("deck_induction_machine_state requires machine type 3 through 7"))
    return state
end

function deck_induction_machine_initial_state(parsed::DeckParser.DeckParseResult; kwargs...)
    state = deck_coupled_dq_machine_initial_state(parsed; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type in (3, 4, 5, 6, 7) ||
        throw(ArgumentError("deck_induction_machine_initial_state requires machine type 3 through 7"))
    return state
end

function deck_direct_current_machine_parameters(
    parsed::DeckParser.DeckParseResult;
    kwargs...,
)
    parameters = deck_coupled_dq_machine_parameters(parsed; kwargs...)
    parameters.machine_type in (8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck_direct_current_machine_parameters requires machine type 8, 9, 10, 11, or 12"))
    return parameters
end

function deck_direct_current_machine_state(
    parsed::DeckParser.DeckParseResult,
    histories;
    kwargs...,
)
    state = deck_coupled_dq_machine_state(parsed, histories; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type in (8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck_direct_current_machine_state requires machine type 8, 9, 10, 11, or 12"))
    return state
end

function deck_direct_current_machine_initial_state(
    parsed::DeckParser.DeckParseResult;
    kwargs...,
)
    state = deck_coupled_dq_machine_initial_state(parsed; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type in (8, 9, 10, 11, 12) ||
        throw(ArgumentError("deck_direct_current_machine_initial_state requires machine type 8, 9, 10, 11, or 12"))
    return state
end

function deck_universal_machine_type4_parameters(parsed::DeckParser.DeckParseResult; kwargs...)
    parameters = deck_induction_machine_parameters(parsed; kwargs...)
    parameters.machine_type == 4 ||
        throw(ArgumentError("deck_universal_machine_type4_parameters requires machine type 4"))
    return parameters
end

function deck_universal_machine_type4_state(parsed::DeckParser.DeckParseResult, histories; kwargs...)
    state = deck_induction_machine_state(parsed, histories; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type == 4 ||
        throw(ArgumentError("deck_universal_machine_type4_state requires machine type 4"))
    return state
end

function deck_universal_machine_type4_initial_state(parsed::DeckParser.DeckParseResult; kwargs...)
    state = deck_induction_machine_initial_state(parsed; kwargs...)
    card1 = _deck_universal_machine_definition(parsed, get(kwargs, :machine_index, 1), 1)
    card1.machine_type == 4 ||
        throw(ArgumentError("deck_universal_machine_type4_initial_state requires machine type 4"))
    return state
end

function run_deck_universal_machine_horizon(
    parsed::DeckParser.DeckParseResult,
    state::InductionMachineState,
    network_steps::AbstractVector;
    machine_index::Int=1,
    time_step_s::Real=deck_fixed_step_horizon(parsed).dt_s,
)
    parameters = deck_coupled_dq_machine_parameters(
        parsed;
        machine_index = machine_index,
        time_step_s = time_step_s,
    )
    call_count = length(network_steps)
    call_count > 0 || throw(ArgumentError("network_steps must not be empty"))
    coil_count = length(parameters.coil_conductances)
    outputs = zeros(Float64, coil_count + 3, call_count)
    currents = zeros(Float64, coil_count, call_count)
    histories = zeros(Float64, coil_count, call_count)
    substitution_count = if parameters.machine_type in (1, 2)
        5
    elseif parameters.machine_type == 6
        3
    elseif parameters.machine_type == 7
        4
    elseif parameters.machine_type == 8
        3
    elseif parameters.machine_type == 9
        2
    elseif parameters.machine_type == 10
        2
    elseif parameters.machine_type == 11
        2
    elseif parameters.machine_type == 12
        2
    elseif parameters.machine_type in (3, 5)
        4
    else
        7
    end
    substitutions = zeros(Float64, substitution_count, call_count)
    d_flux = zeros(Float64, call_count)
    q_flux = zeros(Float64, call_count)
    torque = zeros(Float64, call_count)
    speed = zeros(Float64, call_count)
    angle = zeros(Float64, call_count)
    iterations = zeros(Int, call_count)
    times = zeros(Float64, call_count)
    terminal_voltages = zeros(Float64, 3, call_count)
    rotor_thevenin = zeros(Float64, 3, 3, call_count)
    speed_thevenin = zeros(Float64, call_count)
    torque_impedance = zeros(Float64, call_count)
    for index in eachindex(network_steps)
        step = network_steps[index]
        result = coupled_dq_machine_step!(
            state,
            parameters;
            power_terminal_voltages = step.power_terminal_voltages,
            rotor_thevenin_matrix = step.rotor_thevenin_matrix,
            mechanical_speed_thevenin_rad_s =
                step.mechanical_speed_thevenin_rad_s,
            generated_torque_impedance = step.generated_torque_impedance,
            stator_terminal_voltages =
                get(step, :stator_terminal_voltages, zeros(coil_count - 3)),
            stator_thevenin_matrix =
                get(
                    step,
                    :stator_thevenin_matrix,
                    zeros(coil_count - 3, coil_count - 3),
                ),
            initial_step = index == firstindex(network_steps),
        )
        times[index] = Float64(get(step, :time_s, (index - 1) * Float64(time_step_s)))
        outputs[:, index] .= result.output_values
        currents[:, index] .= result.current_values
        histories[:, index] .= result.history_currents
        substitutions[:, index] .= result.current_substitution_values
        d_flux[index] = result.d_axis_flux
        q_flux[index] = result.q_axis_flux
        torque[index] = result.generated_torque
        speed[index] = result.mechanical_speed_rad_s
        angle[index] = result.mechanical_angle_rad
        iterations[index] = result.iteration_count
        terminal_voltages[:, index] .= step.power_terminal_voltages
        rotor_thevenin[:, :, index] .= step.rotor_thevenin_matrix
        speed_thevenin[index] = step.mechanical_speed_thevenin_rad_s
        torque_impedance[index] = step.generated_torque_impedance
    end
    machine_section = _deck_universal_machine_section(parsed)
    return DeckUniversalMachineHorizon(
        parsed.source,
        parameters.machine_type,
        machine_section.input_layout,
        machine_section.parameter_basis,
        machine_section.remanent_flux_enabled,
        machine_section.initialization_mode,
        machine_section.maximum_shaft_mass_count,
        machine_section.terminal_coupling,
        times,
        outputs,
        currents,
        histories,
        substitutions,
        zeros(Float64, 3, call_count),
        zeros(Float64, 3, call_count),
        0,
        d_flux,
        q_flux,
        torque,
        speed,
        angle,
        iterations,
        state.call_count,
        call_count,
        call_count,
        Symbol[],
        terminal_voltages,
        rotor_thevenin,
        speed_thevenin,
        torque_impedance,
        parameters.series_path_leakage_inductance_h,
        parameters.effective_armature_leakage_inductance_h,
        parameters.effective_compound_field_leakage_inductance_h,
        zeros(Float64, 0, call_count),
        0,
        0,
        zeros(Float64, call_count),
        zeros(Float64, call_count),
        0.0,
        zeros(Float64, call_count),
        0.0,
        0,
        call_count,
        false,
        [:deck_network_control_orchestration],
        Symbol[],
        zeros(Float64, 0, call_count),
        0,
    )
end

function _deck_universal_machine_network_nodes(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    all_terminal_rows = sort!(
        [
            row for row in DeckParser.deck_universal_machine_terminal_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.terminal_index,
    )
    terminal_rows = all_terminal_rows[1:min(3, length(all_terminal_rows))]
    card1 = _deck_universal_machine_definition(parsed, machine_index, 1)
    length(terminal_rows) == 3 ||
        throw(ArgumentError("compensated coupled d/q machine requires three stored power-terminal slots"))
    active_power_positions = card1.machine_type in (6, 7, 8, 9, 10, 11, 12) ? [3] :
        card1.machine_type in (2, 5) ? [2, 3] : [1, 2, 3]
    all(position -> terminal_rows[position].terminal_node_value > 0, active_power_positions) ||
        throw(ArgumentError("type-$(card1.machine_type) active power terminals are missing from the deck network"))
    summaries = [
        row for row in DeckParser.deck_universal_machine_node_summary_rows(parsed)
        if row.machine_index == machine_index
    ]
    length(summaries) == 1 ||
        throw(ArgumentError("coupled d/q machine requires one node/source summary"))
    summary = only(summaries)
    haskey(parsed.node_map, card1.node) ||
        throw(ArgumentError("machine mechanical node is missing from the deck network"))
    ismissing(summary.mechanical_slack_source) &&
        throw(ArgumentError("machine torque-source node is missing from the deck"))
    haskey(parsed.node_map, summary.mechanical_slack_source) ||
        throw(ArgumentError("machine torque-source node is missing from the deck network"))
    source_node = parsed.node_map[summary.mechanical_slack_source]
    source_rows = [
        row for row in DeckParser.deck_over5a_source_rows(parsed)
        if row.node == summary.mechanical_slack_source
    ]
    length(source_rows) == 1 ||
        throw(ArgumentError("coupled d/q machine requires one source at its torque-source node"))
    allowed_drive_source_types = card1.machine_type in (8, 9, 10, 11, 12) ? (11, 14) : (14,)
    only(source_rows).iform in allowed_drive_source_types ||
        throw(ArgumentError("machine torque initialization requires an accepted current source"))
    field_node = ismissing(summary.field_slack_source) ? nothing :
        get(parsed.node_map, summary.field_slack_source, nothing)
    field_source_rows = field_node === nothing ? eltype(source_rows)[] : [
        row for row in DeckParser.deck_over5a_source_rows(parsed)
        if row.node == summary.field_slack_source
    ]
    card1.machine_type in (1, 2, 6, 7, 8) && length(field_source_rows) != 1 &&
        throw(ArgumentError("type-$(card1.machine_type) externally excited machine requires one field source"))
    machine_coils = sort!(
        [
            row for row in DeckParser.deck_universal_machine_coil_rows(parsed)
            if row.machine_index == machine_index
        ];
        by = row -> row.coil_index,
    )
    card1.machine_type == 7 && length(machine_coils) < 5 &&
        throw(ArgumentError("type-7 machine requires five coil rows"))
    q_axis_field_node = card1.machine_type == 7 ? machine_coils[5].terminal_node : nothing
    card1.machine_type == 7 && ismissing(q_axis_field_node) &&
        throw(ArgumentError("type-7 machine q-axis rotor terminal is missing"))
    q_axis_field_source_rows = q_axis_field_node === nothing ? eltype(source_rows)[] : [
        row for row in DeckParser.deck_over5a_source_rows(parsed)
        if row.node == q_axis_field_node
    ]
    card1.machine_type == 7 && length(q_axis_field_source_rows) != 1 &&
        throw(ArgumentError("type-7 machine requires one q-axis rotor source"))
    return (
        power = Int[row.terminal_node_value for row in terminal_rows],
        active_power_positions = active_power_positions,
        mechanical = parsed.node_map[card1.node],
        drive = source_node,
        drive_source = only(source_rows),
        field = field_node,
        field_source = isempty(field_source_rows) ? nothing : only(field_source_rows),
        q_axis_field = q_axis_field_node === nothing ? nothing :
            parsed.node_map[q_axis_field_node],
        q_axis_field_source = isempty(q_axis_field_source_rows) ? nothing :
            only(q_axis_field_source_rows),
    )
end

function _machine_drive_source_value(
    source,
    crest::Real,
    time_s::Real;
    angular_frequency_rad_s::Real=source.sfreq,
)
    return analytic_source_value(
        source.iform,
        Float64(crest),
        source.time1,
        Float64(angular_frequency_rad_s),
        source.tstart,
        source.tstop,
        Float64(time_s),
    )
end

function _seed_machine_mechanical_inertia!(
    system::NodalSystem,
    mechanical_node::Int,
    speed_rad_s::Float64,
)
    seeded_count = 0
    for element in system.elements
        element isa CapacitorBranch || continue
        if element.a == mechanical_node && element.b == 0
            element.v_prev = speed_rad_s
        elseif element.b == mechanical_node && element.a == 0
            element.v_prev = -speed_rad_s
        else
            continue
        end
        element.i_prev = 0.0
        element.i_last = 0.0
        seeded_count += 1
    end
    seeded_count == 1 ||
        throw(ArgumentError("coupled d/q machine mechanical network requires one grounded inertia capacitor"))
    system.v[mechanical_node] = speed_rad_s
    return system
end

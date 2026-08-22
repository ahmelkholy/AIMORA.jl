export ExtendedVSCBoundaryOccurrence,
       ExtendedVSCMetrics,
       ExtendedVSCTrace,
       ExtendedVSCIntegrator,
       configure_extended_vsc,
       advance_extended_vsc!,
       run_extended_vsc!,
       extended_vsc_trace,
       simulate_extended_vsc,
       extended_vsc_task_program_signature_sha256,
       write_extended_vsc_checkpoint,
       read_extended_vsc_checkpoint

const _EXTENDED_VSC_CHECKPOINT_MAGIC = collect(codeunits("AIMORA-EXTENDED-VSC"))
const _EXTENDED_VSC_CHECKPOINT_SCHEMA = UInt32(1)
const _EXTENDED_VSC_TASK_PROGRAM_SIGNATURE_SHA256 = bytes2hex(sha256(codeunits(join((
    "aimora.extended_vsc.task_program.v1",
    "protection:ExtendedVSCProtectionTask",
    "dispatch:ExtendedVSCDispatchTask",
    "measurement:ExtendedVSCMeasurementTask",
    "control:ExtendedVSCControlReader+ExtendedVSCControlComputer+ExtendedVSCControlWriter",
    "pwm:ExtendedVSCDutyReader+ExtendedVSCGateWriter",
), '\n'))))

"""Return the stable callback-program identity required by portable scheduler state."""
extended_vsc_task_program_signature_sha256() =
    _EXTENDED_VSC_TASK_PROGRAM_SIGNATURE_SHA256

struct ExtendedVSCNodeLayout
    dc_positive::Int
    dc_negative::Int
    pole::NTuple{3,Int}
    neutral_pole::Int
    filter::NTuple{3,Int}
    grid::NTuple{3,Int}
    grid_neutral::Int
    node_count::Int
end

function ExtendedVSCNodeLayout(parameters::SwitchDetailedVSC.ExtendedVSCParameters)
    next_node = 6
    neutral_pole = if parameters.wire_form === SwitchDetailedVSC.FourWireForm
        node = next_node
        next_node += 1
        node
    else
        0
    end
    filter = ntuple(phase -> next_node + phase - 1, 3)
    next_node += 3
    grid = if parameters.filter.family === SwitchDetailedVSC.LCLFilter
        nodes = ntuple(phase -> next_node + phase - 1, 3)
        next_node += 3
        nodes
    else
        filter
    end
    grid_neutral = if parameters.wire_form === SwitchDetailedVSC.FourWireForm
        node = next_node
        next_node += 1
        node
    else
        0
    end
    return ExtendedVSCNodeLayout(
        1,
        2,
        (3, 4, 5),
        neutral_pole,
        filter,
        grid,
        grid_neutral,
        next_node - 1,
    )
end

struct ExtendedVSCGridSource <: Branches.EMTElement
    phase_node::Int
    neutral_node::Int
    phase::Int
    resistance_ohm::Float64
    parameters::SwitchDetailedVSC.ExtendedVSCParameters
end

function _extended_vsc_grid_source_voltage(source::ExtendedVSCGridSource, time_s::Float64)
    p = source.parameters
    scenario = p.scenario
    angle = 2.0 * pi * p.frequency_hz * time_s
    phase_offset = (source.phase - 1) * 2.0 * pi / 3.0
    crest = sqrt(2.0 / 3.0) * p.grid_line_line_rms_v
    voltage = crest * cos(angle - phase_offset)
    voltage += crest * scenario.negative_sequence_voltage_ratio * cos(angle + phase_offset)
    voltage += crest * scenario.zero_sequence_voltage_ratio * cos(angle)
    for (order, ratio) in zip(
        scenario.harmonic_orders,
        scenario.harmonic_voltage_ratios,
    )
        sequence_sign = mod(order, 6) == 5 ? -1.0 : 1.0
        voltage += crest * ratio * cos(order * angle - sequence_sign * phase_offset)
    end
    faulted = scenario.fault_kind === :three_phase || source.phase in scenario.faulted_phases
    scenario.fault_start_s <= time_s < scenario.fault_end_s && faulted &&
        (voltage *= scenario.fault_voltage_factor)
    return voltage
end

function _extended_vsc_grid_source_connected(source::ExtendedVSCGridSource, time_s::Float64)
    scenario = source.parameters.scenario
    return !(scenario.island_start_s <= time_s < scenario.reconnect_time_s)
end

function Branches.stamp!(
    admittance,
    right_hand_side,
    source::ExtendedVSCGridSource,
    time_s::Float64,
    _step_s::Float64,
)
    _extended_vsc_grid_source_connected(source, time_s) || return nothing
    conductance = inv(source.resistance_ohm)
    voltage = _extended_vsc_grid_source_voltage(source, time_s)
    Branches.stamp_conductance!(
        admittance,
        source.phase_node,
        source.neutral_node,
        conductance,
    )
    Branches.stamp_history_current!(
        right_hand_side,
        source.phase_node,
        source.neutral_node,
        -conductance * voltage,
    )
    return nothing
end

Branches.stamp!(admittance, right_hand_side, source::ExtendedVSCGridSource,
    time_s::Float64, step_s::Float64, ::Val{Branches.TrapezoidalCompanion}) =
    Branches.stamp!(admittance, right_hand_side, source, time_s, step_s)

Branches.stamp!(admittance, right_hand_side, source::ExtendedVSCGridSource,
    time_s::Float64, step_s::Float64, ::Val{Branches.BackwardEulerCompanion}) =
    Branches.stamp!(admittance, right_hand_side, source, time_s, step_s)

Branches.backward_euler_companion_supported(::ExtendedVSCGridSource) = true
Branches.update!(::ExtendedVSCGridSource, _voltage, _step_s::Float64) = nothing
Branches.update!(::ExtendedVSCGridSource, _voltage, _step_s::Float64,
    ::Val{Branches.TrapezoidalCompanion}) = nothing
Branches.update!(::ExtendedVSCGridSource, _voltage, _step_s::Float64,
    ::Val{Branches.BackwardEulerCompanion}) = nothing

struct ExtendedVSCDCSourceSignal
    parameters::SwitchDetailedVSC.ExtendedVSCParameters
end

function (signal::ExtendedVSCDCSourceSignal)(time_s::Real)
    p = signal.parameters
    time = Float64(time_s)
    factor = p.scenario.dc_sag_start_s <= time < p.scenario.dc_sag_end_s ?
        p.scenario.dc_sag_factor : 1.0
    return factor * p.dc_source_voltage_v
end

"""One exact disturbance, protection, block, or reconnection occurrence."""
struct ExtendedVSCBoundaryOccurrence
    name::Symbol
    time_s::Float64
    tick::Int
    topology_invalidating::Bool
end

"""Bounded physical, numerical, task, event, and energy summary."""
struct ExtendedVSCMetrics
    finite_output::Bool
    exact_boundary_alignment::Bool
    maximum_nodal_kcl_residual_a::Float64
    maximum_nonlinear_kcl_residual_a::Float64
    maximum_terminal_kcl_residual_a::Float64
    relative_energy_residual::Float64
    minimum_dc_link_voltage_v::Float64
    maximum_phase_current_a::Float64
    maximum_neutral_kcl_residual_a::Float64
    controller_sample_count::Int
    controller_write_count::Int
    task_occurrence_count::Int
    boundary_count::Int
    protection_trip_count::Int
    protection_restart_count::Int
    bridge_transition_count::Int
    bridge_refusal_count::Int
    sequence_extractor_settled::Bool
    deterministic_signature::String
end

"""Typed complete retained result for one extended switch-detailed VSC horizon."""
struct ExtendedVSCTrace
    parameters::SwitchDetailedVSC.ExtendedVSCParameters
    time_s::Vector{Float64}
    pole_voltage_v::Matrix{Float64}
    filter_voltage_v::Matrix{Float64}
    grid_voltage_v::Matrix{Float64}
    converter_current_a::Matrix{Float64}
    grid_current_a::Matrix{Float64}
    neutral_current_a::Vector{Float64}
    dc_link_voltage_v::Vector{Float64}
    active_power_w::Vector{Float64}
    reactive_power_var::Vector{Float64}
    positive_sequence_voltage_v::Vector{Float64}
    negative_sequence_voltage_v::Vector{Float64}
    zero_sequence_voltage_v::Vector{Float64}
    controller_frequency_hz::Vector{Float64}
    duty::Matrix{Float64}
    operating_mode::Vector{SwitchDetailedVSC.ExtendedVSCOperatingMode}
    request_disposition::Vector{SwitchDetailedVSC.VSCPlantRequestDisposition}
    stored_energy_j::Vector{Float64}
    dissipated_power_w::Vector{Float64}
    external_power_w::Vector{Float64}
    companion_energy_residual_j::Vector{Float64}
    nodal_kcl_residual_a::Vector{Float64}
    nonlinear_kcl_residual_a::Vector{Float64}
    boundary_occurrences::Vector{ExtendedVSCBoundaryOccurrence}
    sampled_task_occurrences::Vector{OVER16TimestepIntegration.SampledTaskOccurrence}
    metrics::ExtendedVSCMetrics
end

mutable struct ExtendedVSCProtectionState
    violation_start_s::Float64
    healthy_start_s::Float64
    tripped::Bool
    trip_count::Int
    restart_count::Int
end

mutable struct ExtendedVSCRuntime{S,T,L,C,G,D,N}
    parameters::SwitchDetailedVSC.ExtendedVSCParameters
    plant_request::SwitchDetailedVSC.ExtendedVSCPlantRequest
    dispatched_request::SwitchDetailedVSC.ExtendedVSCPlantRequest
    control_state::SwitchDetailedVSC.ExtendedVSCControlState
    held_measurement::SwitchDetailedVSC.ExtendedVSCMeasurement
    held_command::SwitchDetailedVSC.ExtendedVSCControlCommand
    nonlinear_system::S
    bridge_topology::T
    bridge_legs::L
    dc_link_capacitor::C
    converter_filters::G
    grid_filters::D
    shunt_capacitors::N
    shunt_damping::Tuple
    neutral_filter::Union{Nothing,Branches.SeriesRLBranch}
    grid_sources::Tuple
    grid_loads::Tuple
    neutral_grounding::Union{Nothing,Branches.ConductanceBranch}
    dc_source::Branches.TwoTerminalTheveninSource
    layout::ExtendedVSCNodeLayout
    protection::ExtendedVSCProtectionState
    boundary_plan::Vector{Tuple{Symbol,Float64,Bool}}
    boundary_occurrences::Vector{ExtendedVSCBoundaryOccurrence}
    next_boundary_index::Int
    dispatch_count::Int
    measurement_count::Int
    maximum_nonlinear_kcl_residual_a::Float64
    last_nonlinear_kcl_residual_a::Float64
    cumulative_energy_residual_j::Float64
    cumulative_linear_companion_energy_residual_j::Float64
end

struct ExtendedVSCRuntimeCheckpoint{N,C,M,H,P}
    nonlinear_network::N
    control_state::C
    dispatched_request::SwitchDetailedVSC.ExtendedVSCPlantRequest
    held_measurement::M
    held_command::H
    protection::P
    boundary_occurrences::Vector{ExtendedVSCBoundaryOccurrence}
    next_boundary_index::Int
    dispatch_count::Int
    measurement_count::Int
    maximum_nonlinear_kcl_residual_a::Float64
    last_nonlinear_kcl_residual_a::Float64
    cumulative_energy_residual_j::Float64
    cumulative_linear_companion_energy_residual_j::Float64
end

function _extended_vsc_runtime_checkpoint(runtime::ExtendedVSCRuntime)
    return ExtendedVSCRuntimeCheckpoint(
        nonlinear_nodal_checkpoint(runtime.nonlinear_system),
        deepcopy(runtime.control_state),
        runtime.dispatched_request,
        runtime.held_measurement,
        runtime.held_command,
        deepcopy(runtime.protection),
        copy(runtime.boundary_occurrences),
        runtime.next_boundary_index,
        runtime.dispatch_count,
        runtime.measurement_count,
        runtime.maximum_nonlinear_kcl_residual_a,
        runtime.last_nonlinear_kcl_residual_a,
        runtime.cumulative_energy_residual_j,
        runtime.cumulative_linear_companion_energy_residual_j,
    )
end

function _restore_extended_vsc_runtime!(runtime, checkpoint)
    restore_nonlinear_nodal_checkpoint!(
        runtime.nonlinear_system,
        checkpoint.nonlinear_network,
    )
    restore_timestep_state!(runtime.control_state, checkpoint.control_state)
    runtime.dispatched_request = checkpoint.dispatched_request
    runtime.held_measurement = checkpoint.held_measurement
    runtime.held_command = checkpoint.held_command
    restore_timestep_state!(runtime.protection, checkpoint.protection)
    resize!(runtime.boundary_occurrences, length(checkpoint.boundary_occurrences))
    copyto!(runtime.boundary_occurrences, checkpoint.boundary_occurrences)
    runtime.next_boundary_index = checkpoint.next_boundary_index
    runtime.dispatch_count = checkpoint.dispatch_count
    runtime.measurement_count = checkpoint.measurement_count
    runtime.maximum_nonlinear_kcl_residual_a =
        checkpoint.maximum_nonlinear_kcl_residual_a
    runtime.last_nonlinear_kcl_residual_a = checkpoint.last_nonlinear_kcl_residual_a
    runtime.cumulative_energy_residual_j = checkpoint.cumulative_energy_residual_j
    runtime.cumulative_linear_companion_energy_residual_j =
        checkpoint.cumulative_linear_companion_energy_residual_j
    return runtime
end

mutable struct ExtendedVSCRuntimeTransaction{R,C}
    runtime::R
    checkpoint::C
    active::Bool
end

function ExtendedVSCRuntimeTransaction(runtime::ExtendedVSCRuntime)
    return ExtendedVSCRuntimeTransaction(
        runtime,
        _extended_vsc_runtime_checkpoint(runtime),
        false,
    )
end

function _begin_extended_vsc_transaction!(transaction::ExtendedVSCRuntimeTransaction)
    transaction.active && throw(ArgumentError(
        "extended VSC transaction is already active",
    ))
    transaction.checkpoint = _extended_vsc_runtime_checkpoint(transaction.runtime)
    transaction.active = true
    return transaction
end

function _restore_extended_vsc_transaction!(transaction::ExtendedVSCRuntimeTransaction)
    transaction.active || throw(ArgumentError(
        "extended VSC transaction has no active capture",
    ))
    _restore_extended_vsc_runtime!(transaction.runtime, transaction.checkpoint)
    transaction.active = false
    return transaction
end

function _commit_extended_vsc_transaction!(transaction::ExtendedVSCRuntimeTransaction)
    transaction.active || throw(ArgumentError(
        "extended VSC transaction has no active capture",
    ))
    transaction.active = false
    return transaction
end

mutable struct ExtendedVSCRecorder
    time_s::Vector{Float64}
    pole_voltage_v::Matrix{Float64}
    filter_voltage_v::Matrix{Float64}
    grid_voltage_v::Matrix{Float64}
    converter_current_a::Matrix{Float64}
    grid_current_a::Matrix{Float64}
    neutral_current_a::Vector{Float64}
    dc_link_voltage_v::Vector{Float64}
    active_power_w::Vector{Float64}
    reactive_power_var::Vector{Float64}
    positive_sequence_voltage_v::Vector{Float64}
    negative_sequence_voltage_v::Vector{Float64}
    zero_sequence_voltage_v::Vector{Float64}
    controller_frequency_hz::Vector{Float64}
    duty::Matrix{Float64}
    operating_mode::Vector{SwitchDetailedVSC.ExtendedVSCOperatingMode}
    request_disposition::Vector{SwitchDetailedVSC.VSCPlantRequestDisposition}
    stored_energy_j::Vector{Float64}
    dissipated_power_w::Vector{Float64}
    external_power_w::Vector{Float64}
    companion_energy_residual_j::Vector{Float64}
    nodal_kcl_residual_a::Vector{Float64}
    nonlinear_kcl_residual_a::Vector{Float64}
    write_index::Int
end

function ExtendedVSCRecorder(sample_count::Int)
    vector() = Vector{Float64}(undef, sample_count)
    matrix(rows) = Matrix{Float64}(undef, rows, sample_count)
    return ExtendedVSCRecorder(
        vector(), matrix(4), matrix(3), matrix(3), matrix(3), matrix(3),
        vector(), vector(), vector(), vector(), vector(), vector(), vector(),
        vector(), matrix(4),
        Vector{SwitchDetailedVSC.ExtendedVSCOperatingMode}(undef, sample_count),
        Vector{SwitchDetailedVSC.VSCPlantRequestDisposition}(undef, sample_count),
        vector(), vector(), vector(), vector(), vector(), vector(), 0,
    )
end

mutable struct ExtendedVSCIntegrator{R,S,T,C,P}
    runtime::R
    scheduler::S
    transaction::T
    controller_task::C
    pwm_tasks::P
    recorder::ExtendedVSCRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

function _extended_vsc_bridge(parameters, layout)
    ac_nodes = parameters.wire_form === SwitchDetailedVSC.FourWireForm ?
        [layout.pole..., layout.neutral_pole] : collect(layout.pole)
    descriptor = BridgeTopologies.two_level_bridge_topology(
        ac_nodes,
        layout.dc_positive,
        layout.dc_negative,
    )
    diode = Nonlinear.AntiparallelDiodeParameters(
        forward_voltage_v=parameters.diode_forward_voltage_v,
        on_conductance_s=parameters.diode_on_conductance_s,
    )
    valves = map(descriptor.valve_positions) do position
        fidelity = Nonlinear.PowerSemiconductorExtendedFidelity(
            junction_charge=Nonlinear.NonlinearJunctionChargeFidelity(
                1.0e-9,
                100.0,
                0.4;
                voltage_domain_v=(-20.0e3, 20.0e3),
            ),
        )
        Nonlinear.IGBTSwitch(
            position.from_node,
            position.to_node;
            gate_driver=Nonlinear.PowerSemiconductorGateDriver(
                minimum_pulse_width_s=parameters.minimum_gate_pulse_width_s,
            ),
            forward_voltage_drop_v=parameters.semiconductor_forward_voltage_v,
            on_conductance=parameters.semiconductor_on_conductance_s,
            off_conductance=parameters.semiconductor_off_conductance_s,
            antiparallel_diode=diode,
            snubber=Nonlinear.SeriesRCSnubber(
                parameters.snubber_resistance_ohm,
                parameters.snubber_capacitance_f,
            ),
            extended_fidelity=fidelity,
        )
    end
    legs = ntuple(length(ac_nodes)) do leg_index
        Nonlinear.PowerSemiconductorBridgeLeg(
            valves[2 * leg_index - 1],
            valves[2 * leg_index];
            commutation_dead_time_s=parameters.commutation_dead_time_s,
        )
    end
    topology = Nonlinear.PowerSemiconductorBridgeTopology(
        descriptor,
        valves;
        bridge_legs=collect(legs),
    )
    Nonlinear.power_semiconductor_event_localization!(topology)
    return topology, legs, Tuple(valves)
end

function _extended_vsc_initial_command(parameters)
    state = SwitchDetailedVSC.ExtendedVSCControlState(parameters)
    return SwitchDetailedVSC.ExtendedVSCControlCommand(
        parameters.controller.family,
        parameters.wire_form,
        (0.5, 0.5, 0.5, 0.5),
        (0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        state.angle_rad,
        parameters.frequency_hz,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        false,
        false,
        SwitchDetailedVSC.VSCNormalOperation,
        SwitchDetailedVSC.VSCPlantRequestApplied,
    )
end

function _extended_vsc_boundary_plan(parameters)
    scenario = parameters.scenario
    return sort!(
        [
            (:fault_begins, scenario.fault_start_s, false),
            (:fault_clears, scenario.fault_end_s, false),
            (:dc_sag_begins, scenario.dc_sag_start_s, false),
            (:dc_sag_clears, scenario.dc_sag_end_s, false),
            (:island_begins, scenario.island_start_s, true),
            (:grid_reconnects, scenario.reconnect_time_s, true),
            (:commanded_block, scenario.block_time_s, true),
            (:commanded_restart, scenario.restart_time_s, true),
        ];
        by=row -> (row[2], String(row[1])),
    )
end

function _extended_vsc_runtime(parameters, plant_request)
    p = SwitchDetailedVSC.validate_extended_vsc_parameters(parameters)
    layout = ExtendedVSCNodeLayout(p)
    topology, legs, nonlinear_valves = _extended_vsc_bridge(p, layout)
    dc_link_capacitor = Branches.CapacitorBranch(
        layout.dc_positive,
        layout.dc_negative,
        p.dc_link_capacitance_f,
        0.0,
        p.dc_source_voltage_v,
        0.0,
    )
    converter_filters = ntuple(phase -> Branches.SeriesRLBranch(
        layout.pole[phase],
        layout.filter[phase],
        p.filter.converter_resistance_ohm,
        p.filter.converter_inductance_h,
    ), 3)
    grid_filters = p.filter.family === SwitchDetailedVSC.LCLFilter ?
        ntuple(phase -> Branches.SeriesRLBranch(
            layout.filter[phase],
            layout.grid[phase],
            p.filter.grid_resistance_ohm,
            p.filter.grid_inductance_h,
        ), 3) : nothing
    phase_crest = sqrt(2.0 / 3.0) * p.grid_line_line_rms_v
    shunt_capacitors = p.filter.family === SwitchDetailedVSC.SeriesLFilter ? nothing :
        ntuple(phase -> Branches.CapacitorBranch(
            layout.filter[phase],
            layout.grid_neutral,
            p.filter.shunt_capacitance_f,
            0.0,
            phase_crest * cos(-(phase - 1) * 2.0 * pi / 3.0),
            0.0,
        ), 3)
    shunt_damping = p.filter.family === SwitchDetailedVSC.SeriesLFilter ? () :
        ntuple(phase -> Branches.ConductanceBranch(
            layout.filter[phase],
            layout.grid_neutral,
            p.filter.shunt_damping_conductance_s,
        ), 3)
    neutral_filter = p.wire_form === SwitchDetailedVSC.FourWireForm ?
        Branches.SeriesRLBranch(
            layout.neutral_pole,
            layout.grid_neutral,
            p.filter.neutral_resistance_ohm,
            p.filter.neutral_inductance_h,
        ) : nothing
    grid_sources = ntuple(phase -> ExtendedVSCGridSource(
        layout.grid[phase],
        layout.grid_neutral,
        phase,
        p.grid_source_resistance_ohm,
        p,
    ), 3)
    grid_loads = ntuple(phase -> Branches.ConductanceBranch(
        layout.grid[phase],
        layout.grid_neutral,
        inv(p.grid_load_resistance_ohm),
    ), 3)
    neutral_grounding = layout.grid_neutral == 0 ? nothing :
        Branches.ConductanceBranch(
            layout.grid_neutral,
            0,
            inv(p.neutral_grounding_resistance_ohm),
        )
    dc_source = Branches.TwoTerminalTheveninSource(
        layout.dc_positive,
        layout.dc_negative,
        inv(p.dc_source_resistance_ohm),
        ExtendedVSCDCSourceSignal(p),
    )
    elements = Any[topology, dc_link_capacitor, converter_filters...]
    grid_filters === nothing || append!(elements, grid_filters)
    shunt_capacitors === nothing || append!(elements, shunt_capacitors)
    append!(elements, shunt_damping)
    neutral_filter === nothing || push!(elements, neutral_filter)
    push!(elements, dc_source)
    append!(elements, grid_sources)
    append!(elements, grid_loads)
    neutral_grounding === nothing || push!(elements, neutral_grounding)
    linear_system = Nodal.NodalSystem(layout.node_count, elements)
    linear_system.v[layout.dc_positive] = 0.5 * p.dc_source_voltage_v
    linear_system.v[layout.dc_negative] = -0.5 * p.dc_source_voltage_v
    dc_link_capacitor.v_prev = p.dc_source_voltage_v
    for phase in 1:3
        voltage = _extended_vsc_grid_source_voltage(grid_sources[phase], 0.0)
        linear_system.v[layout.pole[phase]] = voltage
        linear_system.v[layout.filter[phase]] = voltage
        linear_system.v[layout.grid[phase]] = voltage
    end
    layout.neutral_pole == 0 || (linear_system.v[layout.neutral_pole] = 0.0)
    layout.grid_neutral == 0 || (linear_system.v[layout.grid_neutral] = 0.0)
    for valve in nonlinear_valves
        Nonlinear.initialize_power_semiconductor_junction_state!(
            valve,
            Branches.branch_voltage(linear_system.v, valve.a, valve.b),
        )
    end
    scales = NonlinearNetworkScales(
        layout.node_count,
        0;
        nominal_voltage_v=p.dc_source_voltage_v,
        nominal_current_a=p.current_limit_a,
    )
    nonlinear_system = NonlinearNodalSystem(
        linear_system,
        nonlinear_valves;
        scales,
    )
    phase_voltage = ntuple(phase -> linear_system.v[layout.grid[phase]] -
        (layout.grid_neutral == 0 ? 0.0 : linear_system.v[layout.grid_neutral]), 3)
    measurement = SwitchDetailedVSC.ExtendedVSCMeasurement(
        phase_voltage,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        phase_voltage,
        0.0,
        p.dc_source_voltage_v,
        0.0,
    )
    return ExtendedVSCRuntime(
        p,
        plant_request,
        plant_request,
        SwitchDetailedVSC.ExtendedVSCControlState(p),
        measurement,
        _extended_vsc_initial_command(p),
        nonlinear_system,
        topology,
        legs,
        dc_link_capacitor,
        converter_filters,
        grid_filters,
        shunt_capacitors,
        shunt_damping,
        neutral_filter,
        grid_sources,
        grid_loads,
        neutral_grounding,
        dc_source,
        layout,
        ExtendedVSCProtectionState(Inf, Inf, false, 0, 0),
        _extended_vsc_boundary_plan(p),
        ExtendedVSCBoundaryOccurrence[],
        1,
        0,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

function _extended_vsc_linear_system(runtime::ExtendedVSCRuntime)
    return nonlinear_linear_system(runtime.nonlinear_system)
end

function _extended_vsc_grid_current(runtime, phase, time_s)
    runtime.grid_filters === nothing || return runtime.grid_filters[phase].i_last
    p = runtime.parameters
    p.filter.family === SwitchDetailedVSC.SeriesLFilter &&
        return runtime.converter_filters[phase].i_last
    source = runtime.grid_sources[phase]
    _extended_vsc_grid_source_connected(source, time_s) || return 0.0
    system = _extended_vsc_linear_system(runtime)
    neutral_voltage = source.neutral_node == 0 ? 0.0 : system.v[source.neutral_node]
    terminal_voltage = system.v[source.phase_node] - neutral_voltage
    return (_extended_vsc_grid_source_voltage(source, time_s) - terminal_voltage) /
        source.resistance_ohm
end

function _extended_vsc_grid_source_current(runtime, phase, time_s)
    source = runtime.grid_sources[phase]
    _extended_vsc_grid_source_connected(source, time_s) || return 0.0
    system = _extended_vsc_linear_system(runtime)
    neutral_voltage = source.neutral_node == 0 ? 0.0 : system.v[source.neutral_node]
    terminal_voltage = system.v[source.phase_node] - neutral_voltage
    return (_extended_vsc_grid_source_voltage(source, time_s) - terminal_voltage) /
        source.resistance_ohm
end

function _extended_vsc_measurement(runtime, time_s)
    system = _extended_vsc_linear_system(runtime)
    layout = runtime.layout
    neutral_voltage = layout.grid_neutral == 0 ? 0.0 : system.v[layout.grid_neutral]
    phase_voltage = ntuple(phase -> system.v[layout.grid[phase]] - neutral_voltage, 3)
    filter_voltage = ntuple(phase -> system.v[layout.filter[phase]] - neutral_voltage, 3)
    converter_current = ntuple(phase -> runtime.converter_filters[phase].i_last, 3)
    grid_current = ntuple(phase -> _extended_vsc_grid_current(runtime, phase, time_s), 3)
    neutral_current = runtime.neutral_filter === nothing ? 0.0 :
        runtime.neutral_filter.i_last
    dc_voltage = system.v[layout.dc_positive] - system.v[layout.dc_negative]
    return SwitchDetailedVSC.ExtendedVSCMeasurement(
        phase_voltage,
        converter_current,
        grid_current,
        filter_voltage,
        neutral_current,
        dc_voltage,
        time_s,
    )
end

struct ExtendedVSCDispatchTask end
function (::ExtendedVSCDispatchTask)(runtime, _time_s, _execution_index)
    runtime.dispatched_request = runtime.plant_request
    runtime.dispatch_count += 1
    return nothing
end

struct ExtendedVSCMeasurementTask end
function (::ExtendedVSCMeasurementTask)(runtime, time_s, _execution_index)
    runtime.held_measurement = _extended_vsc_measurement(runtime, Float64(time_s))
    runtime.measurement_count += 1
    return nothing
end

function _extended_vsc_record_boundary!(runtime, name, time_s, topology_invalidating)
    push!(
        runtime.boundary_occurrences,
        ExtendedVSCBoundaryOccurrence(
            name,
            time_s,
            round(Int, time_s / runtime.parameters.scheduler_tick_s),
            topology_invalidating,
        ),
    )
    return runtime
end

function _extended_vsc_block!(runtime, time_s, name)
    runtime.bridge_topology.blocked && return false
    foreach(leg -> Nonlinear.block_power_semiconductor_bridge!(leg, time_s), runtime.bridge_legs)
    Nonlinear.block_power_semiconductor_topology!(runtime.bridge_topology, time_s)
    runtime.control_state.mode = SwitchDetailedVSC.VSCBlockedOperation
    _extended_vsc_record_boundary!(runtime, name, time_s, true)
    return true
end

function _extended_vsc_restart!(runtime, time_s, name)
    runtime.bridge_topology.blocked || return false
    name === :protection_restart || !runtime.protection.tripped || return false
    Nonlinear.restart_power_semiconductor_topology!(runtime.bridge_topology, time_s)
    foreach(leg -> Nonlinear.restart_power_semiconductor_bridge!(leg, time_s), runtime.bridge_legs)
    runtime.control_state.mode = SwitchDetailedVSC.VSCNormalOperation
    if name === :protection_restart
        runtime.protection.tripped = false
        runtime.protection.violation_start_s = Inf
    end
    _extended_vsc_record_boundary!(runtime, name, time_s, true)
    return true
end

struct ExtendedVSCProtectionTask end
function (::ExtendedVSCProtectionTask)(runtime, time_s, _execution_index)
    time = Float64(time_s)
    measurement = _extended_vsc_measurement(runtime, time)
    p = runtime.parameters
    voltage_magnitude = hypot(
        SwitchDetailedVSC.clarke_transform(measurement.phase_voltage_v).alpha,
        SwitchDetailedVSC.clarke_transform(measurement.phase_voltage_v).beta,
    ) / sqrt(2.0)
    frequency_hz = runtime.control_state.frequency_rad_per_s / (2.0 * pi)
    violated = maximum(abs, measurement.converter_current_a) > p.protection.ac_overcurrent_a ||
        !(p.protection.dc_undervoltage_v <= measurement.dc_link_voltage_v <=
          p.protection.dc_overvoltage_v) ||
        !(p.protection.minimum_frequency_hz <= frequency_hz <=
          p.protection.maximum_frequency_hz) ||
        !(p.protection.minimum_phase_rms_voltage_v <= voltage_magnitude <=
          p.protection.maximum_phase_rms_voltage_v)
    if violated
        runtime.protection.healthy_start_s = Inf
        isfinite(runtime.protection.violation_start_s) ||
            (runtime.protection.violation_start_s = time)
        if !runtime.protection.tripped &&
           time - runtime.protection.violation_start_s >= p.protection.trip_delay_s
            if _extended_vsc_block!(runtime, time, :protection_trip)
                runtime.protection.tripped = true
                runtime.protection.trip_count += 1
            end
        end
    else
        runtime.protection.violation_start_s = Inf
        isfinite(runtime.protection.healthy_start_s) ||
            (runtime.protection.healthy_start_s = time)
        if runtime.protection.tripped &&
           time - runtime.protection.healthy_start_s >= p.protection.restart_delay_s
            _extended_vsc_restart!(runtime, time, :protection_restart) &&
                (runtime.protection.restart_count += 1)
        end
    end
    return nothing
end

struct ExtendedVSCControlReader end
(::ExtendedVSCControlReader)(runtime, _time_s, _sample_index) = runtime.held_measurement

struct ExtendedVSCControlComputer end
function (::ExtendedVSCControlComputer)(runtime, measurement, _time_s, _sample_index)
    if runtime.bridge_topology.blocked
        command = _extended_vsc_initial_command(runtime.parameters)
        runtime.control_state.mode = SwitchDetailedVSC.VSCBlockedOperation
        return SwitchDetailedVSC.ExtendedVSCControlCommand(
            command.controller_family,
            command.wire_form,
            command.duties,
            command.pole_voltage_reference_v,
            command.phase_current_reference_a,
            runtime.control_state.angle_rad,
            runtime.control_state.frequency_rad_per_s / (2.0 * pi),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            false,
            false,
            SwitchDetailedVSC.VSCBlockedOperation,
            runtime.control_state.request_disposition,
        )
    end
    return SwitchDetailedVSC.compute_extended_vsc_control!(
        runtime.control_state,
        measurement,
        runtime.dispatched_request,
        runtime.parameters,
    )
end

struct ExtendedVSCControlWriter end
function (::ExtendedVSCControlWriter)(runtime, command, _time_s, _sample_index)
    runtime.held_command = command
    return nothing
end

struct ExtendedVSCDutyReader
    leg::Int
end
(reader::ExtendedVSCDutyReader)(runtime, _time_s, _cycle_index) =
    runtime.held_command.duties[reader.leg]

struct ExtendedVSCGateWriter
    leg::Int
end
function (writer::ExtendedVSCGateWriter)(runtime, high, time_s, _edge_index)
    Nonlinear.request_power_semiconductor_bridge_pole!(
        runtime.bridge_legs[writer.leg],
        high,
        time_s,
    )
    return nothing
end

function _extended_vsc_scheduler(runtime)
    p = runtime.parameters
    protection = ExactSampledTask(
        :extended_vsc_protection,
        p.controller.control_period_s,
        ExtendedVSCProtectionTask();
        tick_s=p.scheduler_tick_s,
        priority=-40,
        power_history_invalidating=false,
    )
    dispatch = ExactSampledTask(
        :extended_vsc_plant_dispatch,
        p.controller.control_period_s,
        ExtendedVSCDispatchTask();
        tick_s=p.scheduler_tick_s,
        priority=-30,
    )
    measurement = ExactSampledTask(
        :extended_vsc_measurement,
        p.controller.control_period_s,
        ExtendedVSCMeasurementTask();
        tick_s=p.scheduler_tick_s,
        priority=-20,
    )
    controller = ExactSampledControlTask(
        :extended_vsc_control,
        p.controller.control_period_s,
        ExtendedVSCControlReader(),
        ExtendedVSCControlComputer(),
        ExtendedVSCControlWriter();
        tick_s=p.scheduler_tick_s,
        computational_delay_s=p.controller.control_delay_s,
        initial_output=runtime.held_command,
        priority=-10,
        power_history_invalidating=false,
    )
    leg_count = length(runtime.bridge_legs)
    pwm = ntuple(leg_count) do leg
        ExactPWMTask(
            Symbol(:extended_vsc_leg_, leg, :_pwm),
            inv(p.carrier_frequency_hz),
            ExtendedVSCDutyReader(leg),
            ExtendedVSCGateWriter(leg);
            tick_s=p.scheduler_tick_s,
            priority=0,
            power_history_invalidating=false,
        )
    end
    scheduler = ExactSampledTaskScheduler(
        p.scheduler_tick_s;
        tasks=[protection, dispatch, measurement, controller, pwm...],
    )
    return scheduler, controller, pwm
end

function _extended_vsc_task_states(integrator)
    return sampled_task_checkpoint.(integrator.scheduler.tasks)
end

function _restore_extended_vsc_task_states!(integrator, states, occurrence_count)
    foreach(restore_sampled_task_checkpoint!, integrator.scheduler.tasks, states)
    resize!(integrator.scheduler.occurrences, occurrence_count)
    return integrator
end

function _apply_extended_vsc_boundaries!(integrator, time_s)
    runtime = integrator.runtime
    tolerance = 64.0 * eps(Float64) * max(1.0, abs(time_s))
    while runtime.next_boundary_index <= length(runtime.boundary_plan)
        name, boundary_time, topology_invalidating =
            runtime.boundary_plan[runtime.next_boundary_index]
        boundary_time > time_s + tolerance && break
        boundary_time < time_s - tolerance && throw(ArgumentError(
            "extended VSC boundary $name was missed before $time_s s",
        ))
        if name === :commanded_block
            _extended_vsc_block!(runtime, time_s, name)
        elseif name === :commanded_restart
            _extended_vsc_restart!(runtime, time_s, name)
        else
            _extended_vsc_record_boundary!(runtime, name, time_s, topology_invalidating)
        end
        runtime.next_boundary_index += 1
    end
    return runtime
end

function _next_extended_vsc_boundary(runtime, current_time_s, endpoint_time_s)
    runtime.next_boundary_index <= length(runtime.boundary_plan) || return Inf
    boundary_time = runtime.boundary_plan[runtime.next_boundary_index][2]
    boundary_time < current_time_s - runtime.parameters.scheduler_tick_s * 1.0e-9 &&
        throw(ArgumentError("extended VSC boundary precedes the active interval"))
    return boundary_time <= endpoint_time_s ? boundary_time : Inf
end

function _extended_vsc_network_window!(integrator, start_time_s, endpoint_time_s)
    runtime = integrator.runtime
    current_time = start_time_s
    tolerance = runtime.parameters.scheduler_tick_s * 1.0e-9
    while current_time < endpoint_time_s - tolerance
        task_time = next_sampled_task_time(
            integrator.scheduler,
            current_time,
            endpoint_time_s;
            tolerance_s=tolerance,
        )
        boundary_time = _next_extended_vsc_boundary(runtime, current_time, endpoint_time_s)
        segment_endpoint = min(endpoint_time_s, task_time, boundary_time)
        segment_step = segment_endpoint - current_time
        segment_step > tolerance || throw(ArgumentError(
            "extended VSC scheduler produced a nonadvancing interval",
        ))
        initial_measurement = _extended_vsc_measurement(runtime, current_time)
        initial_accounting = _extended_vsc_power(
            runtime,
            initial_measurement,
            current_time,
        )
        initial_stored_energy_j = _extended_vsc_stored_energy(runtime)
        initial_companion_energy_residual_j =
            initial_accounting.bridge.companion_energy_residual_j
        initial_linear_companion = _extended_vsc_linear_companion_accounting(runtime)
        _apply_extended_vsc_boundaries!(integrator, segment_endpoint)
        run_due_sampled_tasks!(integrator.scheduler, runtime, segment_endpoint)
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
            runtime.bridge_topology,
            segment_endpoint,
        )
        topology_invalidating = sampled_task_scheduler_last_run_invalidated_power(
            integrator.scheduler,
        )
        result = advance_nonlinear_step!(
            runtime.nonlinear_system,
            segment_endpoint,
            segment_step;
            discontinuity_treatment=topology_invalidating ?
                :two_backward_euler_half_steps : :none,
            discontinuity_reason=topology_invalidating ? :topology_change : :none,
        )
        result.accepted || throw(something(result.failure))
        runtime.maximum_nonlinear_kcl_residual_a = max(
            runtime.maximum_nonlinear_kcl_residual_a,
            result.diagnostics.maximum_kcl_residual_a,
        )
        runtime.last_nonlinear_kcl_residual_a =
            result.diagnostics.maximum_kcl_residual_a
        final_measurement = _extended_vsc_measurement(runtime, segment_endpoint)
        final_accounting = _extended_vsc_power(
            runtime,
            final_measurement,
            segment_endpoint,
        )
        final_stored_energy_j = _extended_vsc_stored_energy(runtime)
        final_linear_companion = _extended_vsc_linear_companion_accounting(runtime)
        linear_companion_energy_residual_j = 0.5 * segment_step * (
            initial_linear_companion.capacitor_power_w +
            final_linear_companion.capacitor_power_w +
            initial_linear_companion.inductor_power_w +
            final_linear_companion.inductor_power_w -
            initial_linear_companion.inductor_dissipated_power_w -
            final_linear_companion.inductor_dissipated_power_w
        ) - (
            final_linear_companion.capacitor_stored_energy_j -
            initial_linear_companion.capacitor_stored_energy_j +
            final_linear_companion.inductor_stored_energy_j -
            initial_linear_companion.inductor_stored_energy_j
        )
        runtime.cumulative_linear_companion_energy_residual_j +=
            linear_companion_energy_residual_j
        runtime.cumulative_energy_residual_j += 0.5 * segment_step * (
            initial_accounting.external_power_w -
            initial_accounting.dissipated_power_w +
            final_accounting.external_power_w -
            final_accounting.dissipated_power_w
        ) - (final_stored_energy_j - initial_stored_energy_j) -
            (final_accounting.bridge.companion_energy_residual_j -
             initial_companion_energy_residual_j) -
            linear_companion_energy_residual_j
        current_time = segment_endpoint
    end
    return runtime
end

function _extended_vsc_bridge_state(runtime)
    return Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge_topology,
        _extended_vsc_linear_system(runtime).v,
        runtime.parameters.timestep_s,
    )
end

function _extended_vsc_stored_energy(runtime)
    p = runtime.parameters
    system = _extended_vsc_linear_system(runtime)
    layout = runtime.layout
    dc_voltage = system.v[layout.dc_positive] - system.v[layout.dc_negative]
    energy = 0.5 * p.dc_link_capacitance_f * dc_voltage^2
    energy += sum(0.5 * branch.l * branch.i_last^2 for branch in runtime.converter_filters)
    runtime.grid_filters === nothing ||
        (energy += sum(0.5 * branch.l * branch.i_last^2 for branch in runtime.grid_filters))
    runtime.shunt_capacitors === nothing || (energy += sum(runtime.shunt_capacitors) do capacitor
        voltage = system.v[capacitor.a] -
            (capacitor.b == 0 ? 0.0 : system.v[capacitor.b])
        0.5 * capacitor.c * voltage^2
    end)
    runtime.neutral_filter === nothing ||
        (energy += 0.5 * runtime.neutral_filter.l * runtime.neutral_filter.i_last^2)
    return energy + _extended_vsc_bridge_state(runtime).stored_energy_j
end

function _extended_vsc_linear_companion_accounting(runtime)
    system = _extended_vsc_linear_system(runtime)
    capacitor_power_w = 0.0
    capacitor_stored_energy_j = 0.0
    for capacitor in (runtime.dc_link_capacitor,)
        voltage = Branches.branch_voltage(system.v, capacitor.a, capacitor.b)
        capacitor_power_w += voltage * capacitor.i_last
        capacitor_stored_energy_j += 0.5 * capacitor.c * voltage^2
    end
    if runtime.shunt_capacitors !== nothing
        for capacitor in runtime.shunt_capacitors
            voltage = Branches.branch_voltage(system.v, capacitor.a, capacitor.b)
            capacitor_power_w += voltage * capacitor.i_last
            capacitor_stored_energy_j += 0.5 * capacitor.c * voltage^2
        end
    end
    inductor_power_w = 0.0
    inductor_dissipated_power_w = 0.0
    inductor_stored_energy_j = 0.0
    for branch in runtime.converter_filters
        voltage = Branches.branch_voltage(system.v, branch.a, branch.b)
        inductor_power_w += voltage * branch.i_last
        inductor_dissipated_power_w += branch.r * branch.i_last^2
        inductor_stored_energy_j += 0.5 * branch.l * branch.i_last^2
    end
    if runtime.grid_filters !== nothing
        for branch in runtime.grid_filters
            voltage = Branches.branch_voltage(system.v, branch.a, branch.b)
            inductor_power_w += voltage * branch.i_last
            inductor_dissipated_power_w += branch.r * branch.i_last^2
            inductor_stored_energy_j += 0.5 * branch.l * branch.i_last^2
        end
    end
    if runtime.neutral_filter !== nothing
        branch = runtime.neutral_filter
        voltage = Branches.branch_voltage(system.v, branch.a, branch.b)
        inductor_power_w += voltage * branch.i_last
        inductor_dissipated_power_w += branch.r * branch.i_last^2
        inductor_stored_energy_j += 0.5 * branch.l * branch.i_last^2
    end
    return (
        capacitor_power_w=capacitor_power_w,
        capacitor_stored_energy_j=capacitor_stored_energy_j,
        inductor_power_w=inductor_power_w,
        inductor_dissipated_power_w=inductor_dissipated_power_w,
        inductor_stored_energy_j=inductor_stored_energy_j,
    )
end

function _extended_vsc_power(runtime, measurement, time_s)
    p = runtime.parameters
    system = _extended_vsc_linear_system(runtime)
    source_signal = ExtendedVSCDCSourceSignal(p)
    dc_source_voltage = source_signal(time_s)
    dc_link_voltage = measurement.dc_link_voltage_v
    dc_current = (dc_source_voltage - dc_link_voltage) / p.dc_source_resistance_ohm
    external_power = dc_source_voltage * dc_current
    source_loss = p.dc_source_resistance_ohm * dc_current^2
    for phase in 1:3
        source = runtime.grid_sources[phase]
        source_current = _extended_vsc_grid_source_current(runtime, phase, time_s)
        source_voltage = _extended_vsc_grid_source_voltage(source, time_s)
        external_power += source_voltage * source_current
        source_loss += source.resistance_ohm * source_current^2
    end
    filter_loss = sum(branch.r * branch.i_last^2 for branch in runtime.converter_filters)
    runtime.grid_filters === nothing ||
        (filter_loss += sum(branch.r * branch.i_last^2 for branch in runtime.grid_filters))
    runtime.neutral_filter === nothing ||
        (filter_loss += runtime.neutral_filter.r * runtime.neutral_filter.i_last^2)
    isempty(runtime.shunt_damping) || (filter_loss += sum(runtime.shunt_damping) do damping
        voltage = system.v[damping.a] - (damping.b == 0 ? 0.0 : system.v[damping.b])
        damping.g * voltage^2
    end)
    filter_loss += sum(runtime.grid_loads) do load
        voltage = system.v[load.a] - (load.b == 0 ? 0.0 : system.v[load.b])
        load.g * voltage^2
    end
    runtime.neutral_grounding === nothing || begin
        grounding = runtime.neutral_grounding
        voltage = system.v[grounding.a]
        filter_loss += grounding.g * voltage^2
    end
    bridge = _extended_vsc_bridge_state(runtime)
    return (
        external_power_w=external_power,
        dissipated_power_w=source_loss + filter_loss + bridge.semiconductor_loss_w +
            bridge.passive_dissipated_power_w,
        bridge,
    )
end

function _record_extended_vsc_sample!(recorder, runtime, time_s)
    sample = recorder.write_index + 1
    system = _extended_vsc_linear_system(runtime)
    layout = runtime.layout
    measurement = _extended_vsc_measurement(runtime, time_s)
    command = runtime.held_command
    power = SwitchDetailedVSC.instantaneous_three_phase_power(
        measurement.phase_voltage_v,
        measurement.grid_current_a,
        command.angle_rad,
    )
    accounting = _extended_vsc_power(runtime, measurement, time_s)
    neutral_pole_voltage = layout.neutral_pole == 0 ? 0.0 : system.v[layout.neutral_pole]
    for phase in 1:3
        recorder.pole_voltage_v[phase, sample] = system.v[layout.pole[phase]]
        recorder.filter_voltage_v[phase, sample] = measurement.filter_voltage_v[phase]
        recorder.grid_voltage_v[phase, sample] = measurement.phase_voltage_v[phase]
        recorder.converter_current_a[phase, sample] = measurement.converter_current_a[phase]
        recorder.grid_current_a[phase, sample] = measurement.grid_current_a[phase]
        recorder.duty[phase, sample] = command.duties[phase]
    end
    recorder.pole_voltage_v[4, sample] = neutral_pole_voltage
    recorder.duty[4, sample] = command.duties[4]
    recorder.time_s[sample] = time_s
    recorder.neutral_current_a[sample] = measurement.neutral_current_a
    recorder.dc_link_voltage_v[sample] = measurement.dc_link_voltage_v
    recorder.active_power_w[sample] = power.active_w
    recorder.reactive_power_var[sample] = power.reactive_var
    recorder.positive_sequence_voltage_v[sample] = command.positive_sequence_voltage_v
    recorder.negative_sequence_voltage_v[sample] = command.negative_sequence_voltage_v
    recorder.zero_sequence_voltage_v[sample] = command.zero_sequence_voltage_v
    recorder.controller_frequency_hz[sample] = command.frequency_hz
    recorder.operating_mode[sample] = command.mode
    recorder.request_disposition[sample] = command.request_disposition
    recorder.stored_energy_j[sample] = _extended_vsc_stored_energy(runtime)
    recorder.dissipated_power_w[sample] = accounting.dissipated_power_w
    recorder.external_power_w[sample] = accounting.external_power_w
    recorder.companion_energy_residual_j[sample] =
        accounting.bridge.companion_energy_residual_j +
        runtime.cumulative_linear_companion_energy_residual_j
    recorder.nodal_kcl_residual_a[sample] = runtime.last_nonlinear_kcl_residual_a
    recorder.nonlinear_kcl_residual_a[sample] = runtime.last_nonlinear_kcl_residual_a
    recorder.write_index = sample
    return recorder
end

"""Construct the complete coupled extended-VSC private runtime."""
function configure_extended_vsc(
    parameters::SwitchDetailedVSC.ExtendedVSCParameters=
        SwitchDetailedVSC.ExtendedVSCParameters();
    plant_request::SwitchDetailedVSC.ExtendedVSCPlantRequest=
        SwitchDetailedVSC.ExtendedVSCPlantRequest(
            available_active_power_w=parameters.rated_power_va,
        ),
)
    runtime = _extended_vsc_runtime(parameters, plant_request)
    scheduler, controller, pwm = _extended_vsc_scheduler(runtime)
    sample_count = round(Int, parameters.scenario.end_time_s / parameters.timestep_s) + 1
    integrator = ExtendedVSCIntegrator(
        runtime,
        scheduler,
        ExtendedVSCRuntimeTransaction(runtime),
        controller,
        pwm,
        ExtendedVSCRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    task_states = _extended_vsc_task_states(integrator)
    occurrence_count = length(scheduler.occurrences)
    _begin_extended_vsc_transaction!(integrator.transaction)
    try
        _apply_extended_vsc_boundaries!(integrator, 0.0)
        run_due_sampled_tasks!(scheduler, runtime, 0.0)
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
            runtime.bridge_topology,
            0.0,
        )
        result = solve_nonlinear_algebraic_state!(
            runtime.nonlinear_system,
            0.0,
            parameters.timestep_s,
        )
        result.accepted || throw(something(result.failure))
        accepted_voltage = _extended_vsc_linear_system(runtime).v
        Branches.update!(
            runtime.dc_link_capacitor,
            accepted_voltage,
            parameters.timestep_s,
        )
        foreach(runtime.converter_filters) do branch
            Branches.update!(branch, accepted_voltage, parameters.timestep_s)
        end
        runtime.grid_filters === nothing || foreach(runtime.grid_filters) do branch
            Branches.update!(branch, accepted_voltage, parameters.timestep_s)
        end
        runtime.shunt_capacitors === nothing || foreach(runtime.shunt_capacitors) do branch
            Branches.update!(branch, accepted_voltage, parameters.timestep_s)
        end
        runtime.neutral_filter === nothing || Branches.update!(
            runtime.neutral_filter,
            accepted_voltage,
            parameters.timestep_s,
        )
        runtime.maximum_nonlinear_kcl_residual_a = max(
            runtime.maximum_nonlinear_kcl_residual_a,
            result.diagnostics.maximum_kcl_residual_a,
        )
        runtime.last_nonlinear_kcl_residual_a =
            result.diagnostics.maximum_kcl_residual_a
        _commit_extended_vsc_transaction!(integrator.transaction)
    catch
        _restore_extended_vsc_transaction!(integrator.transaction)
        _restore_extended_vsc_task_states!(integrator, task_states, occurrence_count)
        rethrow()
    end
    _record_extended_vsc_sample!(integrator.recorder, runtime, 0.0)
    return integrator
end

"""Advance one exact accepted electrical step atomically."""
function advance_extended_vsc!(integrator::ExtendedVSCIntegrator)
    integrator.failed && throw(ArgumentError(
        "extended VSC integrator is terminally failed: $(integrator.last_failure)",
    ))
    integrator.completed && return false
    p = integrator.runtime.parameters
    final_step = round(Int, p.scenario.end_time_s / p.timestep_s)
    integrator.accepted_step_index >= final_step && begin
        integrator.completed = true
        return false
    end
    next_step = integrator.accepted_step_index + 1
    time_s = next_step * p.timestep_s
    task_states = _extended_vsc_task_states(integrator)
    occurrence_count = length(integrator.scheduler.occurrences)
    recorder_index = integrator.recorder.write_index
    _begin_extended_vsc_transaction!(integrator.transaction)
    try
        _extended_vsc_network_window!(
            integrator,
            (next_step - 1) * p.timestep_s,
            time_s,
        )
        _record_extended_vsc_sample!(integrator.recorder, integrator.runtime, time_s)
        _commit_extended_vsc_transaction!(integrator.transaction)
    catch error
        integrator.transaction.active &&
            _restore_extended_vsc_transaction!(integrator.transaction)
        _restore_extended_vsc_task_states!(integrator, task_states, occurrence_count)
        integrator.recorder.write_index = recorder_index
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow()
    end
    integrator.accepted_step_index = next_step
    integrator.completed = next_step == final_step
    return true
end

function run_extended_vsc!(integrator::ExtendedVSCIntegrator; stop_time_s=nothing)
    stop_step = stop_time_s === nothing ? typemax(Int) :
        round(Int, Float64(stop_time_s) / integrator.runtime.parameters.timestep_s)
    while !integrator.completed && integrator.accepted_step_index < stop_step
        advance_extended_vsc!(integrator)
    end
    return integrator
end

function _extended_vsc_integral(values, time_s)
    total = 0.0
    for index in 2:length(values)
        total += 0.5 * (time_s[index] - time_s[index - 1]) *
            (values[index] + values[index - 1])
    end
    return total
end

function _extended_vsc_signature(recorder, runtime)
    io = IOBuffer()
    print(io, "aimora-extended-vsc-trace-v1\n")
    for values in (
        recorder.time_s,
        vec(recorder.grid_voltage_v),
        vec(recorder.grid_current_a),
        recorder.dc_link_voltage_v,
        vec(recorder.duty),
    )
        for value in values
            write(io, reinterpret(UInt8, [value]))
        end
    end
    print(io, '\n', runtime.parameters.controller.family, '\n',
        runtime.parameters.filter.family, '\n', runtime.parameters.wire_form)
    return bytes2hex(sha256(take!(io)))
end

function _extended_vsc_metrics(integrator)
    recorder = integrator.recorder
    runtime = integrator.runtime
    energy_residual = runtime.cumulative_energy_residual_j
    energy_scale = max(
        abs(_extended_vsc_integral(recorder.external_power_w, recorder.time_s)),
        abs(_extended_vsc_integral(recorder.dissipated_power_w, recorder.time_s)) +
            abs(recorder.stored_energy_j[end] - recorder.stored_energy_j[1]),
        eps(Float64),
    )
    bridge = _extended_vsc_bridge_state(runtime)
    neutral_kcl = runtime.parameters.wire_form === SwitchDetailedVSC.FourWireForm ?
        maximum(abs, vec(sum(recorder.converter_current_a; dims=1)) .+
            recorder.neutral_current_a) :
        maximum(abs, vec(sum(recorder.converter_current_a; dims=1)))
    finite_output = all(isfinite, recorder.grid_voltage_v) &&
        all(isfinite, recorder.grid_current_a) &&
        all(isfinite, recorder.dc_link_voltage_v) &&
        all(isfinite, recorder.duty)
    exact_boundaries = all(runtime.boundary_occurrences) do occurrence
        occurrence.time_s == occurrence.tick * runtime.parameters.scheduler_tick_s
    end
    return ExtendedVSCMetrics(
        finite_output,
        exact_boundaries,
        maximum(recorder.nodal_kcl_residual_a),
        maximum(recorder.nonlinear_kcl_residual_a),
        abs(bridge.terminal_kcl_residual_a),
        abs(energy_residual) / energy_scale,
        minimum(recorder.dc_link_voltage_v),
        maximum(abs, recorder.converter_current_a),
        neutral_kcl,
        integrator.controller_task.sample_count,
        integrator.controller_task.write_count,
        length(integrator.scheduler.occurrences),
        length(runtime.boundary_occurrences),
        runtime.protection.trip_count,
        runtime.protection.restart_count,
        bridge.transition_count,
        bridge.refusal_count,
        runtime.held_command.sequence_extractor_settled,
        _extended_vsc_signature(recorder, runtime),
    )
end

function extended_vsc_trace(integrator::ExtendedVSCIntegrator)
    integrator.completed || throw(ArgumentError(
        "extended VSC trace requires a completed horizon",
    ))
    recorder = integrator.recorder
    return ExtendedVSCTrace(
        integrator.runtime.parameters,
        copy(recorder.time_s),
        copy(recorder.pole_voltage_v),
        copy(recorder.filter_voltage_v),
        copy(recorder.grid_voltage_v),
        copy(recorder.converter_current_a),
        copy(recorder.grid_current_a),
        copy(recorder.neutral_current_a),
        copy(recorder.dc_link_voltage_v),
        copy(recorder.active_power_w),
        copy(recorder.reactive_power_var),
        copy(recorder.positive_sequence_voltage_v),
        copy(recorder.negative_sequence_voltage_v),
        copy(recorder.zero_sequence_voltage_v),
        copy(recorder.controller_frequency_hz),
        copy(recorder.duty),
        copy(recorder.operating_mode),
        copy(recorder.request_disposition),
        copy(recorder.stored_energy_j),
        copy(recorder.dissipated_power_w),
        copy(recorder.external_power_w),
        copy(recorder.companion_energy_residual_j),
        copy(recorder.nodal_kcl_residual_a),
        copy(recorder.nonlinear_kcl_residual_a),
        copy(integrator.runtime.boundary_occurrences),
        copy(integrator.scheduler.occurrences),
        _extended_vsc_metrics(integrator),
    )
end

function simulate_extended_vsc(
    parameters::SwitchDetailedVSC.ExtendedVSCParameters=
        SwitchDetailedVSC.ExtendedVSCParameters();
    plant_request::SwitchDetailedVSC.ExtendedVSCPlantRequest=
        SwitchDetailedVSC.ExtendedVSCPlantRequest(
            available_active_power_w=parameters.rated_power_va,
        ),
)
    integrator = configure_extended_vsc(parameters; plant_request)
    run_extended_vsc!(integrator)
    return extended_vsc_trace(integrator)
end

struct ExtendedVSCCheckpointPayload{I}
    schema::UInt32
    julia_version::VersionNumber
    integrator::I
end

function write_extended_vsc_checkpoint(path::AbstractString, integrator::ExtendedVSCIntegrator)
    integrator.transaction.active && throw(ArgumentError(
        "extended VSC checkpoint cannot be written during a transaction",
    ))
    payload_io = IOBuffer()
    serialize(payload_io, ExtendedVSCCheckpointPayload(
        _EXTENDED_VSC_CHECKPOINT_SCHEMA,
        VERSION,
        integrator,
    ))
    payload = take!(payload_io)
    digest = sha256(payload)
    output_path = abspath(String(path))
    mkpath(dirname(output_path))
    temporary_path, io = mktemp(dirname(output_path))
    try
        write(io, _EXTENDED_VSC_CHECKPOINT_MAGIC)
        _write_checkpoint_integer(io, UInt64(length(payload)))
        write(io, digest)
        write(io, payload)
        close(io)
        mv(temporary_path, output_path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
    return output_path
end

function read_extended_vsc_checkpoint(path::AbstractString)
    input_path = abspath(String(path))
    open(input_path, "r") do io
        read(io, length(_EXTENDED_VSC_CHECKPOINT_MAGIC)) ==
            _EXTENDED_VSC_CHECKPOINT_MAGIC || throw(ArgumentError(
            "file is not an extended VSC checkpoint",
        ))
        payload_length = _read_checkpoint_integer(io, UInt64)
        digest = read(io, 32)
        payload = read(io, Int(payload_length))
        sha256(payload) == digest || throw(ArgumentError(
            "extended VSC checkpoint digest mismatch",
        ))
        object = deserialize(IOBuffer(payload))
        object isa ExtendedVSCCheckpointPayload || throw(ArgumentError(
            "extended VSC checkpoint payload type mismatch",
        ))
        object.schema == _EXTENDED_VSC_CHECKPOINT_SCHEMA || throw(ArgumentError(
            "extended VSC checkpoint schema is unsupported",
        ))
        return object.integrator
    end
end

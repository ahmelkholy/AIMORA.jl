using ..ConverterSystems

mutable struct DualActiveBridgeRuntime{S,N,B,L,C}
    study::S
    network::N
    primary_bridge::B
    secondary_bridge::B
    primary_legs::NTuple{2,Any}
    secondary_legs::NTuple{2,Any}
    primary_controlled::NTuple{4,Any}
    secondary_controlled::NTuple{4,Any}
    primary_diodes::NTuple{4,Any}
    secondary_diodes::NTuple{4,Any}
    leakage_branch::L
    transformer_constraint::C
    primary_dc_positive_node::Int
    primary_dc_negative_node::Int
    secondary_dc_positive_node::Int
    secondary_dc_negative_node::Int
    primary_output_positive_node::Int
    primary_output_negative_node::Int
    secondary_output_positive_node::Int
    secondary_output_negative_node::Int
    time_s::Float64
    primary_source_energy_j::Float64
    secondary_source_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct DualActiveBridgeRecorder
    time_s::Vector{Float64}
    primary_dc_voltage_v::Vector{Float64}
    secondary_dc_voltage_v::Vector{Float64}
    primary_dc_current_a::Vector{Float64}
    secondary_dc_current_a::Vector{Float64}
    primary_bridge_voltage_v::Vector{Float64}
    secondary_bridge_voltage_v::Vector{Float64}
    leakage_current_a::Vector{Float64}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    controlled_conducting_state::BitMatrix
    diode_conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Matrix{Float64}
    diode_junction_temperature_k::Matrix{Float64}
    diode_recovered_charge_c::Matrix{Float64}
    controlled_turn_off_tail_current_a::Matrix{Float64}
    write_index::Int
end

mutable struct DualActiveBridgeIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::DualActiveBridgeRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

mutable struct DualActiveBridgeStepEnergyAccounting{P,L}
    previous_power::P
    previous_linear_companion::L
    primary_source_energy_j::Float64
    secondary_source_energy_j::Float64
    source_and_leakage_dissipated_energy_j::Float64
    semiconductor_dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
end

DualActiveBridgeStepEnergyAccounting(previous_power, previous_linear_companion) =
    DualActiveBridgeStepEnergyAccounting(
        previous_power,
        previous_linear_companion,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )

_dab_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

function _dab_topology_nodes(study::SwitchingDualActiveBridgeStudy)
    primary = Dict(node.name => node.node for node in study.primary_topology.nodes)
    secondary = Dict(node.name => node.node for node in study.secondary_topology.nodes)
    return (
        primary_dc_positive=primary[:dc_positive],
        primary_dc_negative=primary[:dc_negative],
        primary_output_positive=primary[:ac_1],
        primary_output_negative=primary[:ac_2],
        secondary_dc_positive=secondary[:dc_positive],
        secondary_dc_negative=secondary[:dc_negative],
        secondary_output_positive=secondary[:ac_1],
        secondary_output_negative=secondary[:ac_2],
    )
end

function _dab_gate_state(study, time_s)
    modulation = study.specification.modulation
    return ConverterSystems.dual_active_bridge_gate_state(
        time_s,
        study.specification.timing.carrier_frequency_hz,
        modulation.phase_shift_rad,
        modulation.kind;
        primary_inner_phase_shift_rad=modulation.primary_inner_phase_shift_rad,
        secondary_inner_phase_shift_rad=modulation.secondary_inner_phase_shift_rad,
    )
end

function _dab_controlled_valve(study, bridge_index, position, initially_closed)
    detailed = _dab_is_detailed(study)
    parameters = detailed ? study.detailed_semiconductor[bridge_index] : nothing
    timing = study.specification.timing
    return Nonlinear.IGBTSwitch(
        position.from_node,
        position.to_node;
        gate_driver=Nonlinear.PowerSemiconductorGateDriver(
            minimum_pulse_width_s=timing.minimum_pulse_s,
            initially_on=initially_closed,
        ),
        on_conductance=detailed ? parameters.controlled_on_conductance_s : 1.0e6,
        off_conductance=detailed ? parameters.controlled_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? parameters.controlled_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            parameters.snubber_resistance_ohm,
            parameters.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ?
            _detailed_chopper_controlled_fidelity(parameters) : nothing,
        initially_closed,
    )
end

function _dab_diode(study, bridge_index, position, initially_closed)
    detailed = _dab_is_detailed(study)
    parameters = detailed ? study.detailed_semiconductor[bridge_index] : nothing
    return Nonlinear.DiodeValveSwitch(
        position.to_node,
        position.from_node;
        on_conductance=detailed ? parameters.freewheel_on_conductance_s : 1.0e6,
        off_conductance=detailed ? parameters.freewheel_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? parameters.freewheel_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            parameters.snubber_resistance_ohm,
            parameters.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ?
            _detailed_chopper_freewheel_fidelity(parameters) : nothing,
        initially_closed,
    )
end

function _dab_initial_conduction(gates, current, tolerance=0.0)
    controlled = falses(4)
    diodes = falses(4)
    if current > tolerance
        gates[1] ? (controlled[1] = true) : (diodes[2] = true)
        gates[4] ? (controlled[4] = true) : (diodes[3] = true)
    elseif current < -tolerance
        gates[2] ? (controlled[2] = true) : (diodes[1] = true)
        gates[3] ? (controlled[3] = true) : (diodes[4] = true)
    else
        controlled .= gates
    end
    return controlled, diodes
end

function _dab_bridge_devices(study, bridge_index, topology, gates, current)
    controlled_paths, diode_paths = _dab_initial_conduction(gates, current)
    controlled = ntuple(4) do index
        _dab_controlled_valve(
            study,
            bridge_index,
            topology.valve_positions[index],
            controlled_paths[index],
        )
    end
    diodes = ntuple(4) do index
        _dab_diode(
            study,
            bridge_index,
            topology.valve_positions[index],
            diode_paths[index],
        )
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, controlled)
    foreach(Nonlinear.power_semiconductor_event_localization!, diodes)
    legs = (
        Nonlinear.PowerSemiconductorBridgeLeg(
            controlled[1],
            controlled[2];
            commutation_dead_time_s=study.specification.timing.dead_time_s,
        ),
        Nonlinear.PowerSemiconductorBridgeLeg(
            controlled[3],
            controlled[4];
            commutation_dead_time_s=study.specification.timing.dead_time_s,
        ),
    )
    bridge = Nonlinear.PowerSemiconductorBridgeTopology(
        topology,
        collect(controlled);
        bridge_legs=collect(legs),
    )
    return bridge, legs, controlled, diodes
end

_dab_bridge_factor(gates) = Float64(gates[1]) - Float64(gates[3])

function _seed_dual_active_bridge(study::SwitchingDualActiveBridgeStudy)
    nodes = _dab_topology_nodes(study)
    gate = _dab_gate_state(study, study.start_time_s)
    ratio = study.transformer_link.turns_ratio
    leakage_current = study.initial_state.leakage_current_a
    secondary_current = -leakage_current / ratio
    primary_bridge, primary_legs, primary_controlled, primary_diodes =
        _dab_bridge_devices(
            study,
            1,
            study.primary_topology,
            gate.primary_state,
            leakage_current,
        )
    secondary_bridge, secondary_legs, secondary_controlled, secondary_diodes =
        _dab_bridge_devices(
            study,
            2,
            study.secondary_topology,
            gate.secondary_state,
            secondary_current,
        )
    primary_dc_current = _dab_bridge_factor(gate.primary_state) * leakage_current
    secondary_dc_current = _dab_bridge_factor(gate.secondary_state) * secondary_current
    primary_dc_voltage = study.primary_dc_voltage_v -
        study.primary_source_resistance_ohm * primary_dc_current
    secondary_dc_voltage = study.secondary_dc_voltage_v -
        study.secondary_source_resistance_ohm * secondary_dc_current
    primary_leg_voltage(upper_on) = upper_on ? primary_dc_voltage : 0.0
    secondary_leg_voltage(upper_on) = upper_on ? secondary_dc_voltage : 0.0
    primary_output_positive = primary_leg_voltage(gate.primary_state[1])
    primary_output_negative = primary_leg_voltage(gate.primary_state[3])
    secondary_output_positive = secondary_leg_voltage(gate.secondary_state[1])
    secondary_output_negative = secondary_leg_voltage(gate.secondary_state[3])
    transformer_primary_node = maximum((
        nodes.primary_dc_positive,
        nodes.primary_output_positive,
        nodes.primary_output_negative,
        nodes.secondary_dc_positive,
        nodes.secondary_output_positive,
        nodes.secondary_output_negative,
    )) + 1
    constraint_node = transformer_primary_node + 1
    transformer_primary_voltage = primary_output_negative +
        (secondary_output_positive - secondary_output_negative) / ratio
    leakage_voltage = primary_output_positive - transformer_primary_voltage
    primary_source = Branches.TwoTerminalTheveninSource(
        nodes.primary_dc_positive,
        nodes.primary_dc_negative,
        inv(study.primary_source_resistance_ohm),
        _time_s -> study.primary_dc_voltage_v,
    )
    secondary_source = Branches.TwoTerminalTheveninSource(
        nodes.secondary_dc_positive,
        nodes.secondary_dc_negative,
        inv(study.secondary_source_resistance_ohm),
        _time_s -> study.secondary_dc_voltage_v,
    )
    leakage = Branches.SeriesRLBranch(
        nodes.primary_output_positive,
        transformer_primary_node,
        study.transformer_link.leakage_resistance_ohm,
        study.transformer_link.leakage_inductance_h,
        leakage_current,
        leakage_voltage,
        leakage_current,
    )
    constraint = Branches.IdealTransformerVoltageConstraint(
        transformer_primary_node,
        nodes.secondary_output_positive,
        nodes.secondary_output_negative,
        nodes.primary_output_negative,
        constraint_node,
        -ratio,
    )
    devices = (
        primary_controlled...,
        primary_diodes...,
        secondary_controlled...,
        secondary_diodes...,
    )
    node_count = constraint_node
    branches = Any[
            primary_source,
            secondary_source,
            primary_bridge,
            secondary_bridge,
            leakage,
            constraint,
        ]
    _dab_is_detailed(study) || append!(
        branches,
        (primary_diodes..., secondary_diodes...),
    )
    linear_network = Nodal.NodalSystem(node_count, branches)
    linear_network.v[nodes.primary_dc_positive] = primary_dc_voltage
    nodes.primary_dc_negative == 0 ||
        (linear_network.v[nodes.primary_dc_negative] = 0.0)
    linear_network.v[nodes.secondary_dc_positive] = secondary_dc_voltage
    nodes.secondary_dc_negative == 0 ||
        (linear_network.v[nodes.secondary_dc_negative] = 0.0)
    linear_network.v[nodes.primary_output_positive] = primary_output_positive
    linear_network.v[nodes.primary_output_negative] = primary_output_negative
    linear_network.v[nodes.secondary_output_positive] = secondary_output_positive
    linear_network.v[nodes.secondary_output_negative] = secondary_output_negative
    linear_network.v[transformer_primary_node] = transformer_primary_voltage
    linear_network.v[constraint_node] = -leakage_current
    if _dab_is_detailed(study)
        for device in devices
            Nonlinear.initialize_power_semiconductor_junction_state!(
                device,
                Branches.branch_voltage(linear_network.v, device.a, device.b),
            )
        end
    end
    network = _dab_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        devices;
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=max(study.primary_dc_voltage_v,
                study.secondary_dc_voltage_v),
            nominal_current_a=max(abs(leakage_current), 1.0),
        ),
    ) : linear_network
    return DualActiveBridgeRuntime(
        study,
        network,
        primary_bridge,
        secondary_bridge,
        primary_legs,
        secondary_legs,
        primary_controlled,
        secondary_controlled,
        primary_diodes,
        secondary_diodes,
        leakage,
        constraint,
        nodes.primary_dc_positive,
        nodes.primary_dc_negative,
        nodes.secondary_dc_positive,
        nodes.secondary_dc_negative,
        nodes.primary_output_positive,
        nodes.primary_output_negative,
        nodes.secondary_output_positive,
        nodes.secondary_output_negative,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function DualActiveBridgeRecorder(sample_count::Int)
    sample_count > 1 || throw(ArgumentError("DAB recorder requires at least two samples"))
    vectors = ntuple(_ -> zeros(Float64, sample_count), 12)
    matrices = ntuple(_ -> zeros(Float64, 8, sample_count), 4)
    boolean_matrices = ntuple(_ -> falses(8, sample_count), 4)
    return DualActiveBridgeRecorder(
        vectors[1:8]...,
        boolean_matrices...,
        vectors[9:12]...,
        matrices...,
        0,
    )
end

function _dab_linear_network(runtime)
    return _dab_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network
end

function _dab_devices(runtime)
    return (
        runtime.primary_controlled...,
        runtime.primary_diodes...,
        runtime.secondary_controlled...,
        runtime.secondary_diodes...,
    )
end

function _dab_bridge_states(runtime)
    network = _dab_linear_network(runtime)
    step = runtime.study.specification.timing.fixed_step_s
    return (
        Nonlinear.power_semiconductor_bridge_topology_state(
            runtime.primary_bridge,
            network.v,
            step,
        ),
        Nonlinear.power_semiconductor_bridge_topology_state(
            runtime.secondary_bridge,
            network.v,
            step,
        ),
    )
end

function _dab_stored_energy(runtime)
    leakage_energy = 0.5 * runtime.study.transformer_link.leakage_inductance_h *
        runtime.leakage_branch.i_last^2
    _dab_is_detailed(runtime.study) || return leakage_energy
    device_energy = sum(_dab_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
    return leakage_energy + device_energy
end

function _dab_linear_companion(runtime)
    network = _dab_linear_network(runtime)
    voltage = Branches.branch_voltage(
        network.v,
        runtime.leakage_branch.a,
        runtime.leakage_branch.b,
    )
    current = runtime.leakage_branch.i_last
    return (
        terminal_power_w=voltage * current,
        dissipated_power_w=
            runtime.study.transformer_link.leakage_resistance_ohm * current^2,
        stored_energy_j=
            0.5 * runtime.study.transformer_link.leakage_inductance_h * current^2,
    )
end

function _dab_device_companion_residual(runtime)
    _dab_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _dab_devices(runtime)
    )
end

function _dab_power(runtime)
    study = runtime.study
    network = _dab_linear_network(runtime)
    primary_dc_voltage = network.v[runtime.primary_dc_positive_node] -
        (runtime.primary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.primary_dc_negative_node])
    secondary_dc_voltage = network.v[runtime.secondary_dc_positive_node] -
        (runtime.secondary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.secondary_dc_negative_node])
    primary_current = (study.primary_dc_voltage_v - primary_dc_voltage) /
        study.primary_source_resistance_ohm
    secondary_current = (study.secondary_dc_voltage_v - secondary_dc_voltage) /
        study.secondary_source_resistance_ohm
    source_loss = study.primary_source_resistance_ohm * primary_current^2 +
        study.secondary_source_resistance_ohm * secondary_current^2
    leakage_loss = study.transformer_link.leakage_resistance_ohm *
        runtime.leakage_branch.i_last^2
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _dab_devices(runtime)
    )
    return (
        primary_source_w=study.primary_dc_voltage_v * primary_current,
        secondary_source_w=study.secondary_dc_voltage_v * secondary_current,
        source_and_leakage_loss_w=source_loss + leakage_loss,
        semiconductor_loss_w=semiconductor_loss,
    )
end

function _observe_dab_substep!(
    accounting::DualActiveBridgeStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _dab_power(runtime)
    final_linear = _dab_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 : throw(ArgumentError(
            "DAB energy accounting requires an accepted companion method",
        ))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    accounting.primary_source_energy_j += step_s *
        blend(previous_power.primary_source_w, final_power.primary_source_w)
    accounting.secondary_source_energy_j += step_s *
        blend(previous_power.secondary_source_w, final_power.secondary_source_w)
    accounting.source_and_leakage_dissipated_energy_j += step_s * blend(
        previous_power.source_and_leakage_loss_w,
        final_power.source_and_leakage_loss_w,
    )
    accounting.semiconductor_dissipated_energy_j += step_s *
        blend(previous_power.semiconductor_loss_w, final_power.semiconductor_loss_w)
    accounting.linear_companion_energy_residual_j += step_s *
        blend(previous_linear.terminal_power_w, final_linear.terminal_power_w) -
        step_s * blend(
            previous_linear.dissipated_power_w,
            final_linear.dissipated_power_w,
        ) - (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _dab_set_device_state!(device, desired_closed, time_s)
    if desired_closed && !device.closed
        Nonlinear.apply_power_semiconductor_forward_turn_on!(device, time_s)
        return true
    elseif !desired_closed && device.closed
        Nonlinear.apply_power_semiconductor_forward_extinction!(device, time_s)
        return true
    end
    return false
end

function _synchronize_dab_bridge_conduction!(controlled, diodes, current, time_s)
    gates = BitVector(valve.gate_driver.applied_on for valve in controlled)
    tolerance = maximum(valve.holding_current for valve in (controlled..., diodes...))
    controlled_target, diode_target = _dab_initial_conduction(gates, current, tolerance)
    changed = false
    for index in 1:4
        controlled_target[index] || (changed |=
            _dab_set_device_state!(controlled[index], false, time_s))
        diode_target[index] || (changed |=
            _dab_set_device_state!(diodes[index], false, time_s))
    end
    for index in 1:4
        controlled_target[index] && (changed |=
            _dab_set_device_state!(controlled[index], true, time_s))
        diode_target[index] && (changed |=
            _dab_set_device_state!(diodes[index], true, time_s))
    end
    return changed
end

function _synchronize_dab_conduction!(runtime, time_s)
    _dab_is_detailed(runtime.study) || return false
    current = runtime.leakage_branch.i_last
    ratio = runtime.study.transformer_link.turns_ratio
    primary_changed = _synchronize_dab_bridge_conduction!(
        runtime.primary_controlled,
        runtime.primary_diodes,
        current,
        time_s,
    )
    secondary_changed = _synchronize_dab_bridge_conduction!(
        runtime.secondary_controlled,
        runtime.secondary_diodes,
        -current / ratio,
        time_s,
    )
    return primary_changed || secondary_changed
end

function _stabilize_dab_topology!(runtime, time_s, step_s; maximum_iterations=32)
    _dab_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed DAB execution must use the nonlinear D200 device solver",
    ))
    network = _dab_linear_network(runtime)
    for bridge in (runtime.primary_bridge, runtime.secondary_bridge)
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s)
    end
    devices = _dab_devices(runtime)
    for iteration in 1:maximum_iterations
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        selected = nothing
        for (device_index, device) in enumerate(devices)
            action = _three_phase_vsc_topology_action(device, network.v)
            action === nothing && continue
            selected = (device_index, device, action)
            break
        end
        selected === nothing && return iteration
        device_index, device, action = selected
        _apply_three_phase_vsc_topology_action!(device, action.transition, time_s)
        bridge_name = device_index <= 8 ? :primary : :secondary
        local_index = mod1(device_index, 8)
        controlled = local_index <= 4
        position_index = controlled ? local_index : local_index - 4
        topology = bridge_name === :primary ?
            runtime.study.primary_topology : runtime.study.secondary_topology
        position_name = topology.valve_positions[position_index].name
        owner = Symbol(bridge_name, "_", position_name,
            controlled ? "" : "_freewheel")
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            owner,
            true;
            message="accepted DAB switching-state commutation",
        ))
    end
    throw(ArgumentError("DAB switching-state topology failed to stabilize"))
end

function _request_dab_bridge!(bridge, legs, gates, time_s)
    previous = BitVector((legs[1].requested_upper_on, legs[2].requested_upper_on))
    disposition_a = Nonlinear.request_power_semiconductor_bridge_pole!(
        legs[1], gates[1], time_s,
    )
    disposition_b = Nonlinear.request_power_semiconductor_bridge_pole!(
        legs[2], gates[3], time_s,
    )
    changed = previous != BitVector((gates[1], gates[3]))
    return changed, (disposition_a, disposition_b)
end

function _request_dab_gates!(runtime, time_s)
    gate = _dab_gate_state(runtime.study, time_s)
    primary_changed, _ = _request_dab_bridge!(
        runtime.primary_bridge,
        runtime.primary_legs,
        gate.primary_state,
        time_s,
    )
    secondary_changed, _ = _request_dab_bridge!(
        runtime.secondary_bridge,
        runtime.secondary_legs,
        gate.secondary_state,
        time_s,
    )
    changed = primary_changed || secondary_changed
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :dab_phase_shift_edge,
        :dual_active_bridge,
        true;
        message="accepted DAB bridge-leg command with complementary interlock",
    ))
    return changed
end

function _record_dab_sample!(recorder, runtime, sample_index, energy_residual_w)
    study = runtime.study
    network = _dab_linear_network(runtime)
    primary_state, secondary_state = _dab_bridge_states(runtime)
    primary_dc_voltage = network.v[runtime.primary_dc_positive_node] -
        (runtime.primary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.primary_dc_negative_node])
    secondary_dc_voltage = network.v[runtime.secondary_dc_positive_node] -
        (runtime.secondary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.secondary_dc_negative_node])
    primary_dc_current = (study.primary_dc_voltage_v - primary_dc_voltage) /
        study.primary_source_resistance_ohm
    secondary_dc_current = (study.secondary_dc_voltage_v - secondary_dc_voltage) /
        study.secondary_source_resistance_ohm
    recorder.time_s[sample_index] = runtime.time_s
    recorder.primary_dc_voltage_v[sample_index] = primary_dc_voltage
    recorder.secondary_dc_voltage_v[sample_index] = secondary_dc_voltage
    recorder.primary_dc_current_a[sample_index] = primary_dc_current
    recorder.secondary_dc_current_a[sample_index] = secondary_dc_current
    recorder.primary_bridge_voltage_v[sample_index] =
        network.v[runtime.primary_output_positive_node] -
        network.v[runtime.primary_output_negative_node]
    recorder.secondary_bridge_voltage_v[sample_index] =
        network.v[runtime.secondary_output_positive_node] -
        network.v[runtime.secondary_output_negative_node]
    recorder.leakage_current_a[sample_index] = runtime.leakage_branch.i_last
    controlled = (runtime.primary_controlled..., runtime.secondary_controlled...)
    diodes = (runtime.primary_diodes..., runtime.secondary_diodes...)
    requested = vcat(primary_state.requested_gate_state,
        secondary_state.requested_gate_state)
    applied = vcat(primary_state.applied_gate_state, secondary_state.applied_gate_state)
    for index in 1:8
        recorder.requested_gate_state[index, sample_index] = requested[index]
        recorder.applied_gate_state[index, sample_index] = applied[index]
        recorder.controlled_conducting_state[index, sample_index] = controlled[index].closed
        recorder.diode_conducting_state[index, sample_index] = diodes[index].closed
        if _dab_is_detailed(study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled[index])
            diode_state = Nonlinear.power_semiconductor_extended_state(diodes[index])
            recorder.controlled_junction_temperature_k[index, sample_index] =
                controlled_state.junction_temperature_k
            recorder.diode_junction_temperature_k[index, sample_index] =
                diode_state.junction_temperature_k
            recorder.diode_recovered_charge_c[index, sample_index] =
                diode_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[index, sample_index] =
                controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[index, sample_index] = 0.0
            recorder.diode_junction_temperature_k[index, sample_index] = 0.0
            recorder.diode_recovered_charge_c[index, sample_index] = 0.0
            recorder.controlled_turn_off_tail_current_a[index, sample_index] = 0.0
        end
    end
    recorder.stored_energy_j[sample_index] = _dab_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample_index] =
        primary_state.semiconductor_loss_w + secondary_state.semiconductor_loss_w
    recorder.kcl_residual_a[sample_index] = if _dab_is_detailed(study)
        runtime.last_kcl_residual_a
    else
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    end
    recorder.energy_residual_w[sample_index] = energy_residual_w
    recorder.write_index = sample_index
    return recorder
end

function prepare_switching_dual_active_bridge(study::SwitchingDualActiveBridgeStudy)
    runtime = _seed_dual_active_bridge(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    ) + 1
    integrator = DualActiveBridgeIntegrator(
        runtime,
        TimestepTransaction(runtime),
        DualActiveBridgeRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_dab_sample!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_dual_active_bridge!(integrator)
    integrator.failed && throw(ArgumentError(
        "DAB integrator is terminally failed: $(integrator.last_failure)",
    ))
    integrator.completed && return false
    runtime = integrator.runtime
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    final_step = length(integrator.recorder.time_s) - 1
    next_step = integrator.accepted_step_index + 1
    next_step > final_step && begin
        integrator.completed = true
        return false
    end
    endpoint = study.start_time_s + next_step * step
    previous_energy = integrator.recorder.stored_energy_j[next_step]
    previous_power = _dab_power(runtime)
    previous_linear = _dab_linear_companion(runtime)
    previous_device_residual = _dab_device_companion_residual(runtime)
    previous_semiconductor_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _dab_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_dab_gates!(runtime, endpoint)
        accounting = DualActiveBridgeStepEnergyAccounting(previous_power, previous_linear)
        semiconductor_increment = 0.0
        if _dab_is_detailed(study)
            Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                runtime.primary_bridge,
                endpoint,
            )
            Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                runtime.secondary_bridge,
                endpoint,
            )
            conduction_changed = _synchronize_dab_conduction!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) -> _observe_dab_substep!(
                        accounting,
                        runtime,
                        system,
                        time_s,
                        substep_s,
                        companion_method,
                        substep_index,
                    ),
            )
            nonlinear_result.accepted || throw(something(nonlinear_result.failure))
            runtime.last_kcl_residual_a =
                nonlinear_result.diagnostics.maximum_kcl_residual_a
            final_semiconductor_dissipation = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in _dab_devices(runtime)
            )
            semiconductor_increment =
                final_semiconductor_dissipation - previous_semiconductor_dissipation
        else
            event_count_before = length(runtime.events)
            _stabilize_dab_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count_before
            network = _dab_linear_network(runtime)
            if topology_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint - half_step,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_dab_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_dab_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_dab_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    2,
                )
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_dab_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint,
                    step,
                    Branches.TrapezoidalCompanion,
                    1,
                )
            end
            runtime.last_kcl_residual_a =
                maximum(abs, network.y * network.v - network.rhs; init=0.0)
            semiconductor_increment = accounting.semiconductor_dissipated_energy_j
        end
        dissipated_increment =
            accounting.source_and_leakage_dissipated_energy_j + semiconductor_increment
        runtime.primary_source_energy_j += accounting.primary_source_energy_j
        runtime.secondary_source_energy_j += accounting.secondary_source_energy_j
        runtime.dissipated_energy_j += dissipated_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        stored_energy = _dab_stored_energy(runtime)
        device_residual_increment =
            _dab_device_companion_residual(runtime) - previous_device_residual
        energy_residual = (
            accounting.primary_source_energy_j + accounting.secondary_source_energy_j -
            dissipated_increment - (stored_energy - previous_energy) -
            accounting.linear_companion_energy_residual_j - device_residual_increment
        ) / step
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_dab_sample!(integrator.recorder, runtime, next_step + 1, energy_residual)
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _dual_active_bridge_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _dab_linear_network(runtime)
    primary_state, secondary_state = _dab_bridge_states(runtime)
    primary_dc_voltage = network.v[runtime.primary_dc_positive_node] -
        (runtime.primary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.primary_dc_negative_node])
    secondary_dc_voltage = network.v[runtime.secondary_dc_positive_node] -
        (runtime.secondary_dc_negative_node == 0 ? 0.0 :
         network.v[runtime.secondary_dc_negative_node])
    primary_current = (study.primary_dc_voltage_v - primary_dc_voltage) /
        study.primary_source_resistance_ohm
    secondary_current = (study.secondary_dc_voltage_v - secondary_dc_voltage) /
        study.secondary_source_resistance_ohm
    stored_energy = _dab_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        primary_state.deterministic_signature,
        secondary_state.deterministic_signature,
        study.transformer_link.deterministic_signature_sha256,
        repr(runtime.time_s),
        repr(runtime.leakage_branch.i_last),
        repr(stored_energy),
        repr(runtime.dissipated_energy_j),
        string(integrator.accepted_step_index),
        string(length(runtime.events)),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [primary_dc_voltage, secondary_dc_voltage],
        [primary_current, secondary_current],
        BitVector((
            primary_state.requested_gate_state...,
            secondary_state.requested_gate_state...,
        )),
        BitVector((
            primary_state.applied_gate_state...,
            secondary_state.applied_gate_state...,
        )),
        BitVector((
            getfield.(runtime.primary_controlled, :closed)...,
            getfield.(runtime.secondary_controlled, :closed)...,
            getfield.(runtime.primary_diodes, :closed)...,
            getfield.(runtime.secondary_diodes, :closed)...,
        )),
        Float64[],
        [study.transformer_link.leakage_inductance_h * runtime.leakage_branch.i_last],
        Float64[],
        [
            study.specification.modulation.phase_shift_rad,
            study.specification.modulation.primary_inner_phase_shift_rad,
            study.specification.modulation.secondary_inner_phase_shift_rad,
        ],
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    stored_energy_change = stored_energy - recorder.stored_energy_j[1]
    integrated_residual = runtime.primary_source_energy_j +
        runtime.secondary_source_energy_j - runtime.dissipated_energy_j -
        stored_energy_change - runtime.linear_companion_energy_residual_j -
        _dab_device_companion_residual(runtime)
    scale = max(
        abs(runtime.primary_source_energy_j),
        abs(runtime.secondary_source_energy_j),
        abs(runtime.dissipated_energy_j),
        abs(stored_energy_change),
        eps(Float64),
    )
    result = ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=integrator.completed && !integrator.failed,
        status=integrator.completed && !integrator.failed ? :ok : :incomplete,
        events=runtime.events,
        maximum_kcl_residual_a=maximum(recorder.kcl_residual_a),
        relative_charge_residual=maximum(recorder.kcl_residual_a),
        relative_energy_residual=abs(integrated_residual) / scale,
    )
    return ConverterSystems.SwitchingDualActiveBridgeTrace(
        recorder.time_s,
        recorder.primary_dc_voltage_v,
        recorder.secondary_dc_voltage_v,
        recorder.primary_dc_current_a,
        recorder.secondary_dc_current_a,
        recorder.primary_bridge_voltage_v,
        recorder.secondary_bridge_voltage_v,
        recorder.leakage_current_a,
        recorder.requested_gate_state,
        recorder.applied_gate_state,
        recorder.controlled_conducting_state,
        recorder.diode_conducting_state,
        recorder.stored_energy_j,
        recorder.semiconductor_loss_w,
        recorder.kcl_residual_a,
        recorder.energy_residual_w,
        recorder.controlled_junction_temperature_k,
        recorder.diode_junction_temperature_k,
        recorder.diode_recovered_charge_c,
        recorder.controlled_turn_off_tail_current_a,
        result,
    )
end

function execute_switching_dual_active_bridge!(integrator::DualActiveBridgeIntegrator)
    while _advance_dual_active_bridge!(integrator)
    end
    return _dual_active_bridge_result(integrator)
end

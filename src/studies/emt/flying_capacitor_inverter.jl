mutable struct FlyingCapacitorInverterRuntime{S,N,B,C,A,P,L}
    study::S
    network::N
    bridges::B
    controlled_valves::C
    antiparallel_diodes::A
    flying_capacitors::P
    load_branches::L
    dc_positive_node::Int
    dc_negative_node::Int
    phase_nodes::NTuple{3,Int}
    load_neutral_node::Int
    time_s::Float64
    input_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct FlyingCapacitorInverterRecorder
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    flying_capacitor_voltage_v::Matrix{Float64}
    flying_capacitor_current_a::Matrix{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    requested_level::Matrix{Int8}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    controlled_conducting_state::BitMatrix
    antiparallel_diode_conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Matrix{Float64}
    antiparallel_junction_temperature_k::Matrix{Float64}
    antiparallel_recovered_charge_c::Matrix{Float64}
    controlled_turn_off_tail_current_a::Matrix{Float64}
    write_index::Int
end

mutable struct FlyingCapacitorInverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::FlyingCapacitorInverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_flying_capacitor_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_flying_capacitor_linear_network(runtime) =
    _flying_capacitor_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network

function _flying_capacitor_nodes(study)
    leg = ntuple(phase -> Dict(node.name => node.node for node in study.topologies[phase].nodes), 3)
    required = (:dc_positive, :ac_terminal, :dc_negative, :upper_flying_node, :lower_flying_node)
    all(nodes -> all(haskey(nodes, name) for name in required), leg) ||
        throw(ArgumentError("flying-capacitor topology is missing a canonical node"))
    all(phase -> leg[phase][:dc_positive] == leg[1][:dc_positive] &&
        leg[phase][:dc_negative] == leg[1][:dc_negative], 2:3) ||
        throw(ArgumentError("flying-capacitor legs must share exact DC rails"))
    phase_nodes = ntuple(phase -> leg[phase][:ac_terminal], 3)
    length(unique(phase_nodes)) == 3 || throw(ArgumentError(
        "flying-capacitor phase terminals must be distinct",
    ))
    input = only(port for port in study.specification.ports if port.identity === :input_dc)
    output = only(port for port in study.specification.ports if port.identity === :output_ac)
    input.ordered_nodes == (leg[1][:dc_positive], leg[1][:dc_negative]) ||
        throw(ArgumentError("flying-capacitor DC port does not match its rails"))
    output.ordered_nodes == phase_nodes || throw(ArgumentError(
        "flying-capacitor AC port does not match its phase terminals",
    ))
    maximum_node = maximum(node.node for topology in study.topologies for node in topology.nodes)
    return (
        dc_positive=leg[1][:dc_positive],
        dc_negative=leg[1][:dc_negative],
        phase=phase_nodes,
        neutral=maximum_node + 1,
        leg,
    )
end

function _flying_capacitor_reference(study, time_s)
    modulation = study.specification.modulation
    angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
        modulation.phase_shift_rad
    return ntuple(phase -> modulation.modulation_index *
        sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
end

function _flying_capacitor_modulation_state(study, time_s, capacitor_voltage, phase_current)
    return ConverterSystems.converter_flying_capacitor_gate_state(
        _flying_capacitor_reference(study, time_s),
        capacitor_voltage .- 0.5 * study.input_voltage_v,
        phase_current,
        time_s,
        study.specification.timing.carrier_frequency_hz,
        balance_voltage_tolerance_v=study.balance_voltage_tolerance_v,
    )
end

function _flying_capacitor_path_targets(gates, current)
    controlled = falses(4)
    diodes = falses(4)
    active = findall(gates)
    length(active) == 2 || return controlled, diodes
    state = Tuple(active)
    if current >= 0.0
        if state == (1, 2)
            controlled[1:2] .= true
        elseif state == (1, 3)
            controlled[1] = true
            diodes[3] = true
        elseif state == (2, 4)
            controlled[2] = true
            diodes[4] = true
        else
            diodes[3:4] .= true
        end
    else
        if state == (1, 2)
            diodes[1:2] .= true
        elseif state == (1, 3)
            diodes[1] = true
            controlled[3] = true
        elseif state == (2, 4)
            diodes[2] = true
            controlled[4] = true
        else
            controlled[3:4] .= true
        end
    end
    return controlled, diodes
end

function _seed_flying_capacitor_inverter(study)
    nodes = _flying_capacitor_nodes(study)
    initial_voltage = study.initial_state.flying_capacitor_voltage_v
    initial_current = study.initial_state.phase_current_a
    modulation = _flying_capacitor_modulation_state(
        study,
        study.start_time_s,
        collect(initial_voltage),
        collect(initial_current),
    )
    controlled_paths = falses(12)
    diode_paths = falses(12)
    for phase in 1:3
        first = 4 * (phase - 1) + 1
        controlled, diodes = _flying_capacitor_path_targets(
            modulation.requested_valve_state[first:(first + 3)],
            initial_current[phase],
        )
        controlled_paths[first:(first + 3)] .= controlled
        diode_paths[first:(first + 3)] .= diodes
    end
    controlled = ntuple(12) do index
        phase = div(index - 1, 4) + 1
        position = study.topologies[phase].valve_positions[mod(index - 1, 4) + 1]
        _three_level_split_link_controlled_valve(
            study,
            position,
            modulation.requested_valve_state[index],
            controlled_paths[index],
        )
    end
    antiparallel = ntuple(12) do index
        phase = div(index - 1, 4) + 1
        position = study.topologies[phase].valve_positions[mod(index - 1, 4) + 1]
        _three_level_split_link_diode(
            study,
            position.to_node,
            position.from_node,
            diode_paths[index],
        )
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, controlled)
    foreach(Nonlinear.power_semiconductor_event_localization!, antiparallel)
    capacitors = ntuple(3) do phase
        passive = only(study.topologies[phase].passive_positions)
        Branches.CapacitorBranch(
            passive.from_node,
            passive.to_node,
            study.flying_capacitance_f,
            0.0,
            initial_voltage[phase],
            0.0,
        )
    end
    bridges = ntuple(3) do phase
        first = 4 * (phase - 1) + 1
        Nonlinear.PowerSemiconductorBridgeTopology(
            study.topologies[phase],
            Nonlinear.PowerSemiconductorSwitch[controlled[first:(first + 3)]...];
            passives=Branches.CapacitorBranch[capacitors[phase]],
        )
    end
    source = Branches.TwoTerminalTheveninSource(
        nodes.dc_positive,
        nodes.dc_negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    pole_voltage = zeros(3)
    upper_node_voltage = zeros(3)
    lower_node_voltage = zeros(3)
    for phase in 1:3
        first = 4 * (phase - 1) + 1
        gates = modulation.requested_valve_state[first:(first + 3)]
        if gates == BitVector((true, true, false, false))
            pole_voltage[phase] = study.input_voltage_v
            upper_node_voltage[phase] = study.input_voltage_v
            lower_node_voltage[phase] = study.input_voltage_v - initial_voltage[phase]
        elseif gates == BitVector((true, false, true, false))
            upper_node_voltage[phase] = study.input_voltage_v
            lower_node_voltage[phase] = study.input_voltage_v - initial_voltage[phase]
            pole_voltage[phase] = lower_node_voltage[phase]
        elseif gates == BitVector((false, true, false, true))
            lower_node_voltage[phase] = 0.0
            upper_node_voltage[phase] = initial_voltage[phase]
            pole_voltage[phase] = upper_node_voltage[phase]
        else
            lower_node_voltage[phase] = 0.0
            upper_node_voltage[phase] = initial_voltage[phase]
            pole_voltage[phase] = 0.0
        end
    end
    neutral_voltage = sum(pole_voltage) / 3.0
    loads = ntuple(3) do phase
        Branches.SeriesRLBranch(
            nodes.phase[phase],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            initial_current[phase],
            pole_voltage[phase] - neutral_voltage,
            initial_current[phase],
        )
    end
    node_count = nodes.neutral
    elements = _flying_capacitor_is_detailed(study) ?
        Any[source, bridges..., loads...] : Any[source, bridges..., antiparallel..., loads...]
    linear_network = Nodal.NodalSystem(node_count, elements)
    linear_network.v[nodes.dc_positive] = study.input_voltage_v
    nodes.dc_negative == 0 || (linear_network.v[nodes.dc_negative] = 0.0)
    for phase in 1:3
        linear_network.v[nodes.phase[phase]] = pole_voltage[phase]
        linear_network.v[nodes.leg[phase][:upper_flying_node]] = upper_node_voltage[phase]
        linear_network.v[nodes.leg[phase][:lower_flying_node]] = lower_node_voltage[phase]
    end
    linear_network.v[nodes.neutral] = neutral_voltage
    devices = (controlled..., antiparallel...)
    if _flying_capacitor_is_detailed(study)
        for valve in devices
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _flying_capacitor_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        devices;
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=study.specification.rated_bases.current_a,
        ),
        options=NonlinearSolveOptions(
            maximum_iterations=50,
            current_absolute_tolerance_a=1.0e-9,
            current_relative_tolerance=2.5e-9,
            voltage_absolute_tolerance_v=1.0e-9,
            scaled_step_tolerance=1.0e-12,
            maximum_condition_estimate=2.0e13,
        ),
    ) : linear_network
    return FlyingCapacitorInverterRuntime(
        study,
        network,
        bridges,
        controlled,
        antiparallel,
        capacitors,
        loads,
        nodes.dc_positive,
        nodes.dc_negative,
        nodes.phase,
        nodes.neutral,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function FlyingCapacitorInverterRecorder(sample_count)
    return FlyingCapacitorInverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Int8, 3, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, 12, sample_count),
        0,
    )
end

_flying_capacitor_devices(runtime) =
    (runtime.controlled_valves..., runtime.antiparallel_diodes...)

function _flying_capacitor_stored_energy(runtime)
    network = _flying_capacitor_linear_network(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    energy += sum(runtime.flying_capacitors; init=0.0) do capacitor
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        0.5 * capacitor.c * voltage^2
    end
    _flying_capacitor_is_detailed(runtime.study) || return energy
    return energy + sum(_flying_capacitor_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _flying_capacitor_power(runtime)
    study = runtime.study
    network = _flying_capacitor_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    source_current = (study.input_voltage_v - dc_voltage) / study.source_resistance_ohm
    load_loss = study.load_resistance_ohm *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _flying_capacitor_devices(runtime)
    )
    return (
        input_w=study.input_voltage_v * source_current,
        load_w=load_loss,
        source_loss_w=study.source_resistance_ohm * source_current^2,
        semiconductor_loss_w=semiconductor_loss,
        source_current_a=source_current,
    )
end

function _flying_capacitor_linear_companion(runtime)
    network = _flying_capacitor_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for branch in runtime.load_branches
        voltage = Branches.branch_voltage(network.v, branch.a, branch.b)
        terminal_power += voltage * branch.i_last
        dissipated_power += runtime.study.load_resistance_ohm * branch.i_last^2
        stored_energy += 0.5 * runtime.study.load_inductance_h * branch.i_last^2
    end
    for capacitor in runtime.flying_capacitors
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        terminal_power += voltage * capacitor.i_last
        stored_energy += 0.5 * capacitor.c * voltage^2
    end
    return (; terminal_power_w=terminal_power, dissipated_power_w=dissipated_power, stored_energy_j=stored_energy)
end

function _flying_capacitor_device_residual(runtime)
    _flying_capacitor_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _flying_capacitor_devices(runtime)
    )
end

function _observe_flying_capacitor_substep!(accounting, runtime, _system, _time, step, method, _index)
    final_power = _flying_capacitor_power(runtime)
    final_linear = _flying_capacitor_linear_companion(runtime)
    weight = method === Branches.TrapezoidalCompanion ? 0.5 :
        method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("flying-capacitor inverter received an unknown companion method"))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    accounting.input_energy_j += step * blend(previous_power.input_w, final_power.input_w)
    accounting.load_energy_j += step * blend(previous_power.load_w, final_power.load_w)
    accounting.source_dissipated_energy_j += step *
        blend(previous_power.source_loss_w, final_power.source_loss_w)
    accounting.semiconductor_dissipated_energy_j += step *
        blend(previous_power.semiconductor_loss_w, final_power.semiconductor_loss_w)
    accounting.linear_companion_energy_residual_j +=
        step * blend(previous_linear.terminal_power_w, final_linear.terminal_power_w) -
        step * blend(previous_linear.dissipated_power_w, final_linear.dissipated_power_w) -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _current_flying_capacitor_modulation(runtime, time_s)
    network = _flying_capacitor_linear_network(runtime)
    capacitor_voltage = [
        Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        for capacitor in runtime.flying_capacitors
    ]
    phase_current = [branch.i_last for branch in runtime.load_branches]
    return _flying_capacitor_modulation_state(
        runtime.study,
        time_s,
        capacitor_voltage,
        phase_current,
    )
end

function _request_flying_capacitor_gates!(runtime, time_s)
    modulation = _current_flying_capacitor_modulation(runtime, time_s)
    timing = runtime.study.specification.timing
    changed = false
    for phase in 1:3
        first = 4 * (phase - 1) + 1
        requested = modulation.requested_valve_state[first:(first + 3)]
        BridgeTopologies.bridge_topology_state_is_allowed(
            runtime.study.topologies[phase],
            requested,
        ) || throw(ArgumentError("flying-capacitor modulation requested a prohibited state"))
        for position in 1:4
            valve = runtime.controlled_valves[first + position - 1]
            desired = requested[position]
            desired == valve.gate_driver.commanded_on && continue
            Nonlinear.request_power_semiconductor_gate!(
                valve,
                desired,
                time_s;
                earliest_transition_time_s=desired ? time_s + timing.dead_time_s : time_s,
            )
            runtime.bridges[phase].transition_count += 1
            changed = true
        end
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :flying_capacitor_pwm_command,
        :flying_capacitor_bridge,
        true;
        message="accepted flying-capacitor carrier-PWM gate command",
    ))
    return changed
end

function _synchronize_flying_capacitor_conduction!(runtime, time_s)
    controlled_target = falses(12)
    diode_target = falses(12)
    for phase in 1:3
        first = 4 * (phase - 1) + 1
        applied = BitVector(
            runtime.controlled_valves[first + position - 1].gate_driver.applied_on
            for position in 1:4
        )
        current = runtime.load_branches[phase].i_last
        controlled, diodes = if count(applied) == 2 &&
            BridgeTopologies.bridge_topology_state_is_allowed(runtime.study.topologies[phase], applied)
            _flying_capacitor_path_targets(applied, current)
        elseif current >= 0.0
            targets = falses(4)
            targets[3:4] .= true
            falses(4), targets
        else
            targets = falses(4)
            targets[1:2] .= true
            falses(4), targets
        end
        controlled_target[first:(first + 3)] .= controlled
        diode_target[first:(first + 3)] .= diodes
    end
    changed = false
    for index in 1:12
        changed |= _set_full_bridge_device_state!(
            runtime.controlled_valves[index],
            controlled_target[index] && runtime.controlled_valves[index].gate_driver.applied_on,
            time_s,
        )
        changed |= _set_full_bridge_device_state!(
            runtime.antiparallel_diodes[index],
            diode_target[index],
            time_s,
        )
    end
    return changed
end

function _stabilize_flying_capacitor!(runtime, time_s, step_s)
    _flying_capacitor_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed flying-capacitor execution requires the nonlinear solver",
    ))
    network = _flying_capacitor_linear_network(runtime)
    foreach(bridge -> Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s), runtime.bridges)
    for _iteration in 1:16
        changed = _synchronize_flying_capacitor_conduction!(runtime, time_s)
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        changed || return nothing
    end
    throw(ArgumentError("flying-capacitor topology failed to stabilize"))
end

function _record_flying_capacitor!(recorder, runtime, sample, energy_residual)
    network = _flying_capacitor_linear_network(runtime)
    power = _flying_capacitor_power(runtime)
    modulation = _current_flying_capacitor_modulation(runtime, runtime.time_s)
    neutral_voltage = network.v[runtime.load_neutral_node]
    recorder.time_s[sample] = runtime.time_s
    recorder.source_current_a[sample] = power.source_current_a
    recorder.requested_level[:, sample] .= modulation.requested_level
    for phase in 1:3
        capacitor = runtime.flying_capacitors[phase]
        recorder.flying_capacitor_voltage_v[phase, sample] =
            Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        recorder.flying_capacitor_current_a[phase, sample] = capacitor.i_last
        recorder.phase_voltage_v[phase, sample] = network.v[runtime.phase_nodes[phase]] - neutral_voltage
        recorder.phase_current_a[phase, sample] = runtime.load_branches[phase].i_last
    end
    for index in 1:12
        controlled = runtime.controlled_valves[index]
        diode = runtime.antiparallel_diodes[index]
        recorder.requested_gate_state[index, sample] = controlled.gate_driver.commanded_on
        recorder.applied_gate_state[index, sample] = controlled.gate_driver.applied_on
        recorder.controlled_conducting_state[index, sample] = controlled.closed
        recorder.antiparallel_diode_conducting_state[index, sample] = diode.closed
        if _flying_capacitor_is_detailed(runtime.study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            diode_state = Nonlinear.power_semiconductor_extended_state(diode)
            recorder.controlled_junction_temperature_k[index, sample] = controlled_state.junction_temperature_k
            recorder.antiparallel_junction_temperature_k[index, sample] = diode_state.junction_temperature_k
            recorder.antiparallel_recovered_charge_c[index, sample] = diode_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[index, sample] = controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_recovered_charge_c[index, sample] = 0.0
            recorder.controlled_turn_off_tail_current_a[index, sample] = 0.0
        end
    end
    recorder.stored_energy_j[sample] = _flying_capacitor_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _flying_capacitor_is_detailed(runtime.study) ?
        runtime.last_kcl_residual_a : maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual
    recorder.write_index = sample
    return recorder
end

function prepare_switching_flying_capacitor(study::ConverterSystems.SwitchingFlyingCapacitorStudy)
    runtime = _seed_flying_capacitor_inverter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    ) + 1
    integrator = FlyingCapacitorInverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        FlyingCapacitorInverterRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_flying_capacitor!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_flying_capacitor_inverter!(integrator)
    integrator.failed && throw(ArgumentError(
        "flying-capacitor integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _flying_capacitor_power(runtime)
    previous_linear = _flying_capacitor_linear_companion(runtime)
    previous_stored = _flying_capacitor_stored_energy(runtime)
    previous_device_residual = _flying_capacitor_device_residual(runtime)
    previous_device_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _flying_capacitor_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_flying_capacitor_gates!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        device_increment = 0.0
        if _flying_capacitor_is_detailed(study)
            foreach(bridge -> Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, endpoint), runtime.bridges)
            conduction_changed = _synchronize_flying_capacitor_conduction!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s, method, index) ->
                    _observe_flying_capacitor_substep!(
                        accounting,
                        runtime,
                        system,
                        time_s,
                        substep_s,
                        method,
                        index,
                    ),
            )
            nonlinear_result.accepted || throw(something(nonlinear_result.failure))
            runtime.last_kcl_residual_a = nonlinear_result.diagnostics.maximum_kcl_residual_a
            final_device_dissipation = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in _flying_capacitor_devices(runtime)
            )
            device_increment = final_device_dissipation - previous_device_dissipation
        else
            _stabilize_flying_capacitor!(runtime, endpoint, step)
            network = _flying_capacitor_linear_network(runtime)
            if gate_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                Nodal.solve_algebraic_state!(network, endpoint - half_step, half_step, backward_euler)
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_flying_capacitor_substep!(accounting, runtime, network,
                    endpoint - half_step, half_step, Branches.BackwardEulerCompanion, 1)
                _stabilize_flying_capacitor!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(network, endpoint, half_step, backward_euler)
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_flying_capacitor_substep!(accounting, runtime, network,
                    endpoint, half_step, Branches.BackwardEulerCompanion, 2)
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_flying_capacitor_substep!(accounting, runtime, network,
                    endpoint, step, Branches.TrapezoidalCompanion, 1)
            end
            runtime.last_kcl_residual_a = maximum(abs, network.y * network.v - network.rhs; init=0.0)
            device_increment = accounting.semiconductor_dissipated_energy_j
        end
        final_stored = _flying_capacitor_stored_energy(runtime)
        device_residual_increment = _flying_capacitor_device_residual(runtime) - previous_device_residual
        energy_residual = (
            accounting.input_energy_j - accounting.load_energy_j -
            accounting.source_dissipated_energy_j - device_increment -
            (final_stored - previous_stored) -
            accounting.linear_companion_energy_residual_j - device_residual_increment
        ) / step
        runtime.input_energy_j += accounting.input_energy_j
        runtime.dissipated_energy_j += accounting.load_energy_j +
            accounting.source_dissipated_energy_j + device_increment
        runtime.linear_companion_energy_residual_j += accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_flying_capacitor!(integrator.recorder, runtime, next_step + 1, energy_residual)
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _flying_capacitor_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _flying_capacitor_linear_network(runtime)
    power = _flying_capacitor_power(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    phase_voltage = [network.v[node] - neutral_voltage for node in runtime.phase_nodes]
    phase_current = [branch.i_last for branch in runtime.load_branches]
    capacitor_voltage = [
        Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        for capacitor in runtime.flying_capacitors
    ]
    bridge_states = ntuple(phase -> Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridges[phase], network.v, study.specification.timing.fixed_step_s,
    ), 3)
    stored_energy = _flying_capacitor_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        (state.deterministic_signature for state in bridge_states)...,
        repr(runtime.time_s),
        repr(phase_voltage),
        repr(phase_current),
        repr(capacitor_voltage),
        string(integrator.accepted_step_index),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [capacitor_voltage..., phase_voltage...],
        [power.source_current_a, phase_current...],
        BitVector(valve.gate_driver.commanded_on for valve in runtime.controlled_valves),
        BitVector(valve.gate_driver.applied_on for valve in runtime.controlled_valves),
        BitVector((
            getfield.(runtime.controlled_valves, :closed)...,
            getfield.(runtime.antiparallel_diodes, :closed)...,
        )),
        study.flying_capacitance_f .* capacitor_voltage,
        study.load_inductance_h .* phase_current,
        Float64[],
        Float64.(recorder.requested_level[:, end]),
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    energy_scale = max(abs(runtime.input_energy_j), abs(runtime.dissipated_energy_j),
        abs(stored_energy), eps(Float64))
    integrated_residual = runtime.input_energy_j - runtime.dissipated_energy_j -
        (stored_energy - recorder.stored_energy_j[1]) -
        runtime.linear_companion_energy_residual_j - _flying_capacitor_device_residual(runtime)
    result = ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=integrator.completed && !integrator.failed,
        status=integrator.completed && !integrator.failed ? :ok : :incomplete,
        events=runtime.events,
        maximum_kcl_residual_a=maximum(recorder.kcl_residual_a),
        relative_charge_residual=maximum(abs, vec(sum(recorder.phase_current_a; dims=1))),
        relative_energy_residual=abs(integrated_residual) / energy_scale,
    )
    return ConverterSystems.SwitchingFlyingCapacitorTrace(
        recorder.time_s,
        recorder.source_current_a,
        recorder.flying_capacitor_voltage_v,
        recorder.flying_capacitor_current_a,
        recorder.phase_voltage_v,
        recorder.phase_current_a,
        recorder.requested_level,
        recorder.requested_gate_state,
        recorder.applied_gate_state,
        recorder.controlled_conducting_state,
        recorder.antiparallel_diode_conducting_state,
        recorder.stored_energy_j,
        recorder.semiconductor_loss_w,
        recorder.kcl_residual_a,
        recorder.energy_residual_w,
        recorder.controlled_junction_temperature_k,
        recorder.antiparallel_junction_temperature_k,
        recorder.antiparallel_recovered_charge_c,
        recorder.controlled_turn_off_tail_current_a,
        result,
    )
end

function execute_switching_flying_capacitor!(integrator::FlyingCapacitorInverterIntegrator)
    while _advance_flying_capacitor_inverter!(integrator)
    end
    return _flying_capacitor_result(integrator)
end

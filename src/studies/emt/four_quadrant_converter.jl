mutable struct FullBridgeConverterRuntime{S,N,T,L,V,D}
    study::S
    network::N
    bridge::T
    legs::L
    controlled_valves::V
    freewheel_diodes::D
    load_branch::Branches.SeriesRLBranch
    dc_positive_node::Int
    dc_negative_node::Int
    output_positive_node::Int
    output_negative_node::Int
    time_s::Float64
    input_energy_j::Float64
    load_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct FullBridgeConverterRecorder
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    output_voltage_v::Vector{Float64}
    load_current_a::Vector{Float64}
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

mutable struct FullBridgeConverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::FullBridgeConverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

mutable struct FullBridgeStepEnergyAccounting{P,L}
    previous_power::P
    previous_linear_companion::L
    input_energy_j::Float64
    load_energy_j::Float64
    source_dissipated_energy_j::Float64
    semiconductor_dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
end

FullBridgeStepEnergyAccounting(previous_power, previous_linear_companion) =
    FullBridgeStepEnergyAccounting(
        previous_power,
        previous_linear_companion,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )

function _full_bridge_nodes(study::SwitchingFullBridgeConverterStudy)
    topology_nodes = Dict(node.name => node.node for node in study.topology.nodes)
    all(haskey(topology_nodes, name) for name in
        (:dc_positive, :dc_negative, :ac_1, :ac_2)) || throw(ArgumentError(
            "full-bridge full bridge is missing a canonical terminal",
        ))
    input_port = only(port for port in study.specification.ports if port.identity === :input_dc)
    output_identity = study.specification.selection.family ===
        ConverterSystems.SinglePhaseTwoLevelBridge ? :output_ac : :output_dc
    output_port = only(port for port in study.specification.ports if
        port.identity === output_identity)
    input_port.ordered_nodes ==
        (topology_nodes[:dc_positive], topology_nodes[:dc_negative]) ||
        throw(ArgumentError("full-bridge input port does not match its DC rails"))
    output_port.ordered_nodes == (topology_nodes[:ac_1], topology_nodes[:ac_2]) ||
        throw(ArgumentError("full-bridge output port does not match its bridge poles"))
    return (
        dc_positive=topology_nodes[:dc_positive],
        dc_negative=topology_nodes[:dc_negative],
        output_positive=topology_nodes[:ac_1],
        output_negative=topology_nodes[:ac_2],
    )
end

function _full_bridge_positive_command(study, time_s)
    modulation = study.specification.modulation
    duty = if study.specification.selection.family ===
            ConverterSystems.SinglePhaseTwoLevelBridge
        angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
            modulation.phase_shift_rad
        0.5 * (1.0 + modulation.modulation_index * sin(angle))
    else
        modulation.duty
    end
    pwm = ConverterSystems.converter_pwm_gate_state(
        [duty],
        time_s,
        study.specification.timing.carrier_frequency_hz,
        modulation.kind,
    )
    return pwm.requested_valve_state[1]
end

_full_bridge_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

function _full_bridge_linear_network(runtime)
    return _full_bridge_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network
end

_full_bridge_gate_pattern(positive::Bool) =
    positive ? Bool[true, false, false, true] : Bool[false, true, true, false]

function _full_bridge_controlled_valve(
    study,
    position,
    initially_commanded,
    initially_closed,
)
    parameters = study.detailed_semiconductor
    timing = study.specification.timing
    detailed = _full_bridge_is_detailed(study)
    device = Nonlinear.IGBTSwitch(
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
    if initially_commanded != initially_closed
        Nonlinear.request_power_semiconductor_gate!(
            device,
            initially_commanded,
            study.start_time_s;
            earliest_transition_time_s=study.start_time_s + timing.dead_time_s,
        )
    end
    return device
end

function _full_bridge_freewheel_diode(study, position, initially_closed)
    parameters = study.detailed_semiconductor
    detailed = _full_bridge_is_detailed(study)
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

function _full_bridge_initial_paths(positive_command, current_a)
    aligned = iszero(current_a) || (positive_command ? current_a > 0.0 : current_a < 0.0)
    controlled = aligned ? _full_bridge_gate_pattern(positive_command) : falses(4)
    diode = aligned ? falses(4) :
        current_a > 0.0 ? _full_bridge_gate_pattern(false) :
        _full_bridge_gate_pattern(true)
    return BitVector(controlled), BitVector(diode)
end

function _seed_full_bridge_converter(study::SwitchingFullBridgeConverterStudy)
    nodes = _full_bridge_nodes(study)
    positive_command = _full_bridge_positive_command(study, study.start_time_s)
    gate_pattern = _full_bridge_gate_pattern(positive_command)
    initial_current = study.initial_state.load_current_a
    controlled_paths, diode_paths =
        _full_bridge_initial_paths(positive_command, initial_current)
    controlled = ntuple(4) do index
        _full_bridge_controlled_valve(
            study,
            study.topology.valve_positions[index],
            gate_pattern[index],
            controlled_paths[index],
        )
    end
    diodes = ntuple(4) do index
        _full_bridge_freewheel_diode(
            study,
            study.topology.valve_positions[index],
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
        study.topology,
        collect(controlled);
        bridge_legs=collect(legs),
    )
    source_current = positive_command ? initial_current : -initial_current
    source_drop = study.source_resistance_ohm * source_current
    dc_voltage = study.input_voltage_v - source_drop
    output_voltage = positive_command ? dc_voltage : -dc_voltage
    source = Branches.TwoTerminalTheveninSource(
        nodes.dc_positive,
        nodes.dc_negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    load = Branches.SeriesRLBranch(
        nodes.output_positive,
        nodes.output_negative,
        study.load_resistance_ohm,
        study.load_inductance_h,
        initial_current,
        output_voltage,
        initial_current,
    )
    node_count = maximum((
        nodes.dc_positive,
        nodes.dc_negative,
        nodes.output_positive,
        nodes.output_negative,
    ))
    elements = _full_bridge_is_detailed(study) ?
        Any[source, bridge, load] : Any[source, bridge, diodes..., load]
    linear_network = Nodal.NodalSystem(node_count, elements)
    linear_network.v[nodes.dc_positive] = dc_voltage
    nodes.dc_negative == 0 || (linear_network.v[nodes.dc_negative] = 0.0)
    linear_network.v[nodes.output_positive] = 0.5 * output_voltage
    linear_network.v[nodes.output_negative] = -0.5 * output_voltage
    if _full_bridge_is_detailed(study)
        for valve in (controlled..., diodes...)
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _full_bridge_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        (controlled..., diodes...);
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=max(abs(initial_current), 1.0),
        ),
    ) : linear_network
    return FullBridgeConverterRuntime(
        study,
        network,
        bridge,
        legs,
        controlled,
        diodes,
        load,
        nodes.dc_positive,
        nodes.dc_negative,
        nodes.output_positive,
        nodes.output_negative,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function FullBridgeConverterRecorder(sample_count::Int)
    sample_count > 1 || throw(ArgumentError(
        "full-bridge recorder requires at least two samples",
    ))
    vectors = ntuple(_ -> zeros(Float64, sample_count), 9)
    matrices = ntuple(_ -> zeros(Float64, 4, sample_count), 4)
    return FullBridgeConverterRecorder(
        vectors[1:5]...,
        falses(4, sample_count),
        falses(4, sample_count),
        falses(4, sample_count),
        falses(4, sample_count),
        vectors[6:9]...,
        matrices...,
        0,
    )
end

function _full_bridge_devices(runtime)
    return (runtime.controlled_valves..., runtime.freewheel_diodes...)
end

function _full_bridge_stored_energy(runtime)
    load_energy = 0.5 * runtime.study.load_inductance_h * runtime.load_branch.i_last^2
    _full_bridge_is_detailed(runtime.study) || return load_energy
    device_energy = sum(_full_bridge_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
    return load_energy + device_energy
end

function _full_bridge_linear_companion(runtime)
    network = _full_bridge_linear_network(runtime)
    voltage = Branches.branch_voltage(
        network.v,
        runtime.load_branch.a,
        runtime.load_branch.b,
    )
    current = runtime.load_branch.i_last
    return (
        terminal_power_w=voltage * current,
        dissipated_power_w=runtime.study.load_resistance_ohm * current^2,
        stored_energy_j=0.5 * runtime.study.load_inductance_h * current^2,
    )
end

function _full_bridge_device_companion_residual(runtime)
    _full_bridge_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _full_bridge_devices(runtime)
    )
end

function _full_bridge_power(runtime)
    study = runtime.study
    network = _full_bridge_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    input_current = (study.input_voltage_v - dc_voltage) / study.source_resistance_ohm
    load_current = runtime.load_branch.i_last
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _full_bridge_devices(runtime)
    )
    return (
        input_w=study.input_voltage_v * input_current,
        load_w=study.load_resistance_ohm * load_current^2,
        source_loss_w=study.source_resistance_ohm * input_current^2,
        semiconductor_loss_w=semiconductor_loss,
    )
end

function _observe_full_bridge_substep!(
    accounting::FullBridgeStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _full_bridge_power(runtime)
    final_linear = _full_bridge_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 : throw(ArgumentError(
            "full-bridge energy accounting requires an accepted companion method",
        ))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    accounting.input_energy_j += step_s *
        blend(previous_power.input_w, final_power.input_w)
    accounting.load_energy_j += step_s *
        blend(previous_power.load_w, final_power.load_w)
    accounting.source_dissipated_energy_j += step_s *
        blend(previous_power.source_loss_w, final_power.source_loss_w)
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

function _stabilize_full_bridge_topology!(
    runtime,
    time_s,
    step_s;
    maximum_iterations=32,
)
    _full_bridge_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed full-bridge execution must use the nonlinear D200 device solver",
    ))
    network = _full_bridge_linear_network(runtime)
    Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
        runtime.bridge,
        time_s,
    )
    devices = _full_bridge_devices(runtime)
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
        controlled = device_index <= length(runtime.controlled_valves)
        position_index = controlled ? device_index :
            device_index - length(runtime.controlled_valves)
        position_name = runtime.study.topology.valve_positions[position_index].name
        owner = controlled ? position_name : Symbol(position_name, :_freewheel)
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            owner,
            true;
            message="accepted full-bridge switching-state commutation",
        ))
    end
    throw(ArgumentError("full-bridge switching-state topology failed to stabilize"))
end

function _set_full_bridge_device_state!(device, desired_closed, time_s)
    if desired_closed && !device.closed
        Nonlinear.apply_power_semiconductor_forward_turn_on!(device, time_s)
        return true
    elseif !desired_closed && device.closed
        Nonlinear.apply_power_semiconductor_forward_extinction!(device, time_s)
        return true
    end
    return false
end

function _synchronize_full_bridge_conduction!(runtime, time_s)
    controlled = runtime.controlled_valves
    applied = BitVector(valve.gate_driver.applied_on for valve in controlled)
    current = runtime.load_branch.i_last
    tolerance = maximum(valve.holding_current for valve in _full_bridge_devices(runtime))
    positive_applied = applied[1] && applied[4]
    negative_applied = applied[2] && applied[3]
    controlled_target = falses(4)
    diode_target = falses(4)
    if positive_applied && current >= -tolerance
        controlled_target[[1, 4]] .= true
    elseif negative_applied && current <= tolerance
        controlled_target[[2, 3]] .= true
    elseif current > tolerance
        diode_target[[2, 3]] .= true
    elseif current < -tolerance
        diode_target[[1, 4]] .= true
    elseif positive_applied
        controlled_target[[1, 4]] .= true
    elseif negative_applied
        controlled_target[[2, 3]] .= true
    end
    changed = false
    for index in 1:4
        changed |= _set_full_bridge_device_state!(
            controlled[index],
            controlled_target[index],
            time_s,
        )
        changed |= _set_full_bridge_device_state!(
            runtime.freewheel_diodes[index],
            diode_target[index],
            time_s,
        )
    end
    return changed
end

function _request_full_bridge_gates!(runtime, time_s)
    positive = _full_bridge_positive_command(runtime.study, time_s)
    previous = runtime.legs[1].requested_upper_on
    disposition_a = Nonlinear.request_power_semiconductor_bridge_pole!(
        runtime.legs[1],
        positive,
        time_s,
    )
    disposition_b = Nonlinear.request_power_semiconductor_bridge_pole!(
        runtime.legs[2],
        !positive,
        time_s,
    )
    changed = positive != previous
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        positive ? :positive_bridge_command : :negative_bridge_command,
        :full_bridge,
        true;
        message="accepted bipolar PWM command with complementary-leg interlock",
    ))
    return changed, (disposition_a, disposition_b)
end

function _record_full_bridge_sample!(recorder, runtime, sample_index, energy_residual_w)
    network = _full_bridge_linear_network(runtime)
    study = runtime.study
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    input_current = (study.input_voltage_v - dc_voltage) / study.source_resistance_ohm
    output_voltage = network.v[runtime.output_positive_node] -
        network.v[runtime.output_negative_node]
    power = _full_bridge_power(runtime)
    recorder.time_s[sample_index] = runtime.time_s
    recorder.input_voltage_v[sample_index] = study.input_voltage_v
    recorder.input_current_a[sample_index] = input_current
    recorder.output_voltage_v[sample_index] = output_voltage
    recorder.load_current_a[sample_index] = runtime.load_branch.i_last
    for index in 1:4
        controlled = runtime.controlled_valves[index]
        diode = runtime.freewheel_diodes[index]
        recorder.requested_gate_state[index, sample_index] =
            controlled.gate_driver.commanded_on
        recorder.applied_gate_state[index, sample_index] =
            controlled.gate_driver.applied_on
        recorder.controlled_conducting_state[index, sample_index] = controlled.closed
        recorder.diode_conducting_state[index, sample_index] = diode.closed
        if _full_bridge_is_detailed(study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            diode_state = Nonlinear.power_semiconductor_extended_state(diode)
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
    recorder.stored_energy_j[sample_index] = _full_bridge_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample_index] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample_index] = if _full_bridge_is_detailed(study)
        runtime.last_kcl_residual_a
    else
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    end
    recorder.energy_residual_w[sample_index] = energy_residual_w
    recorder.write_index = sample_index
    return recorder
end

function prepare_switching_full_bridge_converter(
    study::SwitchingFullBridgeConverterStudy,
)
    runtime = _seed_full_bridge_converter(study)
    sample_count = Int(round(
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    )) + 1
    integrator = FullBridgeConverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        FullBridgeConverterRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_full_bridge_sample!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_full_bridge_converter!(integrator)
    integrator.failed && throw(ArgumentError(
        "full-bridge integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _full_bridge_power(runtime)
    previous_linear = _full_bridge_linear_companion(runtime)
    previous_device_residual = _full_bridge_device_companion_residual(runtime)
    previous_semiconductor_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _full_bridge_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed, _ = _request_full_bridge_gates!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        semiconductor_increment = 0.0
        if _full_bridge_is_detailed(study)
            Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                runtime.bridge,
                endpoint,
            )
            conduction_changed = _synchronize_full_bridge_conduction!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) -> _observe_full_bridge_substep!(
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
                for valve in _full_bridge_devices(runtime)
            )
            semiconductor_increment =
                final_semiconductor_dissipation - previous_semiconductor_dissipation
        else
            event_count_before = length(runtime.events)
            _stabilize_full_bridge_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count_before
            network = _full_bridge_linear_network(runtime)
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
                _observe_full_bridge_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_full_bridge_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_full_bridge_substep!(
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
                _observe_full_bridge_substep!(
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
        source_and_semiconductor_increment =
            accounting.source_dissipated_energy_j + semiconductor_increment
        runtime.input_energy_j += accounting.input_energy_j
        runtime.load_energy_j += accounting.load_energy_j
        runtime.dissipated_energy_j += source_and_semiconductor_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        stored_energy = _full_bridge_stored_energy(runtime)
        device_residual_increment =
            _full_bridge_device_companion_residual(runtime) - previous_device_residual
        energy_residual = (
            accounting.input_energy_j - accounting.load_energy_j -
            source_and_semiconductor_increment - (stored_energy - previous_energy) -
            accounting.linear_companion_energy_residual_j - device_residual_increment
        ) / step
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_full_bridge_sample!(
            integrator.recorder,
            runtime,
            next_step + 1,
            energy_residual,
        )
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _full_bridge_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _full_bridge_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    input_current = (study.input_voltage_v - dc_voltage) / study.source_resistance_ohm
    output_voltage = network.v[runtime.output_positive_node] -
        network.v[runtime.output_negative_node]
    stored_energy = _full_bridge_stored_energy(runtime)
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        study.specification.timing.fixed_step_s,
    )
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        bridge_state.deterministic_signature,
        repr(runtime.time_s),
        repr(runtime.load_branch.i_last),
        repr(output_voltage),
        repr(stored_energy),
        string(integrator.accepted_step_index),
        string(length(runtime.events)),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, output_voltage],
        [input_current, runtime.load_branch.i_last],
        bridge_state.requested_gate_state,
        bridge_state.applied_gate_state,
        BitVector((
            getfield.(runtime.controlled_valves, :closed)...,
            getfield.(runtime.freewheel_diodes, :closed)...,
        )),
        Float64[],
        [study.load_inductance_h * runtime.load_branch.i_last],
        Float64[],
        [_average_full_bridge_polarity(study, runtime.time_s)],
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    stored_energy_change = stored_energy - recorder.stored_energy_j[1]
    integrated_residual = runtime.input_energy_j - runtime.load_energy_j -
        runtime.dissipated_energy_j - stored_energy_change -
        runtime.linear_companion_energy_residual_j -
        _full_bridge_device_companion_residual(runtime)
    scale = max(
        abs(runtime.input_energy_j),
        abs(runtime.load_energy_j),
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
    return ConverterSystems.SwitchingFullBridgeConverterTrace(
        recorder.time_s,
        recorder.input_voltage_v,
        recorder.input_current_a,
        recorder.output_voltage_v,
        recorder.load_current_a,
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

function execute_switching_full_bridge_converter!(
    integrator::FullBridgeConverterIntegrator,
)
    while _advance_full_bridge_converter!(integrator)
    end
    return _full_bridge_result(integrator)
end

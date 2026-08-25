using ..ConverterSystems

mutable struct InterleavedChopperRuntime{S,N,C}
    study::S
    network::N
    bridges::Vector{Any}
    inductors::Vector{Any}
    capacitor::C
    input_positive_node::Int
    switching_nodes::Vector{Int}
    output_node::Int
    time_s::Float64
    input_energy_j::Float64
    load_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct InterleavedChopperRecorder
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    output_voltage_v::Vector{Float64}
    load_current_a::Vector{Float64}
    channel_switch_voltage_v::Matrix{Float64}
    channel_inductor_current_a::Matrix{Float64}
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

mutable struct InterleavedChopperIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::InterleavedChopperRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

mutable struct InterleavedChopperStepEnergyAccounting{P,L}
    previous_power::P
    previous_linear_companion::L
    input_energy_j::Float64
    load_energy_j::Float64
    series_dissipated_energy_j::Float64
    semiconductor_dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
end

InterleavedChopperStepEnergyAccounting(previous_power, previous_linear_companion) =
    InterleavedChopperStepEnergyAccounting(
        previous_power,
        previous_linear_companion,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )

_interleaved_chopper_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

function _interleaved_chopper_nodes(study::SwitchingInterleavedChopperStudy)
    channel_nodes = [Dict(node.name => node.node for node in topology.nodes)
        for topology in study.topologies]
    input_port = only(port for port in study.specification.ports if port.identity === :input_dc)
    output_port = only(port for port in study.specification.ports if port.identity === :output_dc)
    return (
        input_positive=input_port.ordered_nodes[1],
        negative=input_port.ordered_nodes[2],
        switching=[nodes[:output] for nodes in channel_nodes],
        output=output_port.ordered_nodes[1],
    )
end

function _interleaved_chopper_gate_state(study, time_s)
    count = study.specification.selection.channel_count
    phases = ConverterSystems.interleaved_carrier_phases_rad(count)
    pwm = ConverterSystems.converter_pwm_gate_state(
        fill(study.specification.modulation.duty, count),
        time_s,
        study.specification.timing.carrier_frequency_hz,
        study.specification.modulation.kind;
        carrier_phases_rad=phases,
    )
    return BitVector(pwm.requested_valve_state[1:2:end])
end

function _interleaved_controlled_valve(study, channel, position, initially_closed)
    detailed = _interleaved_chopper_is_detailed(study)
    parameters = detailed ? study.detailed_semiconductor[channel] : nothing
    return Nonlinear.IGBTSwitch(
        position.from_node,
        position.to_node;
        gate_driver=Nonlinear.PowerSemiconductorGateDriver(initially_on=initially_closed),
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

function _interleaved_freewheel_diode(study, channel, position, initially_closed)
    detailed = _interleaved_chopper_is_detailed(study)
    parameters = detailed ? study.detailed_semiconductor[channel] : nothing
    return Nonlinear.DiodeValveSwitch(
        position.from_node,
        position.to_node;
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

function _seed_interleaved_chopper(study::SwitchingInterleavedChopperStudy)
    nodes = _interleaved_chopper_nodes(study)
    initial = study.initial_state
    requested = _interleaved_chopper_gate_state(study, study.start_time_s)
    input_current = sum(initial.channel_inductor_current_a[requested])
    input_node_voltage = study.input_voltage_v -
        study.source_resistance_ohm * input_current
    bridges = Vector{Any}(undef, length(study.topologies))
    devices = Any[]
    for channel in eachindex(study.topologies)
        controlled_position, diode_position = study.topologies[channel].valve_positions
        controlled = _interleaved_controlled_valve(
            study,
            channel,
            controlled_position,
            requested[channel],
        )
        diode = _interleaved_freewheel_diode(
            study,
            channel,
            diode_position,
            !requested[channel],
        )
        Nonlinear.power_semiconductor_event_localization!(controlled)
        Nonlinear.power_semiconductor_event_localization!(diode)
        bridges[channel] = Nonlinear.PowerSemiconductorBridgeTopology(
            study.topologies[channel],
            [controlled, diode],
        )
        append!(devices, (controlled, diode))
    end
    source = Branches.TwoTerminalTheveninSource(
        nodes.input_positive,
        nodes.negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    inductors = Vector{Any}(undef, length(study.topologies))
    for channel in eachindex(inductors)
        switch_voltage = requested[channel] ? input_node_voltage : 0.0
        inductor_voltage = switch_voltage - initial.output_voltage_v
        current = initial.channel_inductor_current_a[channel]
        inductors[channel] = Branches.SeriesRLBranch(
            nodes.switching[channel],
            nodes.output,
            study.inductor_resistance_ohm[channel],
            study.inductance_h[channel],
            current,
            inductor_voltage,
            current,
        )
    end
    load_current = initial.output_voltage_v / study.load_resistance_ohm
    capacitor_current = sum(initial.channel_inductor_current_a) - load_current
    capacitor = Branches.CapacitorBranch(
        nodes.output,
        nodes.negative,
        study.capacitance_f,
        capacitor_current,
        initial.output_voltage_v,
        capacitor_current,
    )
    load = Branches.ConductanceBranch(
        nodes.output,
        nodes.negative,
        inv(study.load_resistance_ohm),
    )
    node_count = maximum((
        nodes.input_positive,
        nodes.switching...,
        nodes.output,
        (node.node for topology in study.topologies for node in topology.nodes)...,
    ))
    branches = Any[source]
    append!(branches, bridges)
    append!(branches, inductors)
    append!(branches, (capacitor, load))
    linear_network = Nodal.NodalSystem(node_count, branches)
    nodes.negative == 0 || (linear_network.v[nodes.negative] = 0.0)
    linear_network.v[nodes.input_positive] = input_node_voltage
    linear_network.v[nodes.output] = initial.output_voltage_v
    for channel in eachindex(nodes.switching)
        linear_network.v[nodes.switching[channel]] =
            requested[channel] ? input_node_voltage : 0.0
    end
    if _interleaved_chopper_is_detailed(study)
        for device in devices
            Nonlinear.initialize_power_semiconductor_junction_state!(
                device,
                Branches.branch_voltage(linear_network.v, device.a, device.b),
            )
        end
    end
    network = _interleaved_chopper_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        Tuple(devices);
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=max(sum(abs, initial.channel_inductor_current_a), 1.0),
        ),
    ) : linear_network
    return InterleavedChopperRuntime(
        study,
        network,
        bridges,
        inductors,
        capacitor,
        nodes.input_positive,
        nodes.switching,
        nodes.output,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function InterleavedChopperRecorder(channel_count::Int, sample_count::Int)
    2 <= channel_count <= 8 && sample_count > 1 || throw(ArgumentError(
        "interleaved recorder requires two through eight channels and at least two samples",
    ))
    vector_traces = ntuple(_ -> zeros(Float64, sample_count), 9)
    matrix_traces = ntuple(_ -> zeros(Float64, channel_count, sample_count), 6)
    boolean_traces = ntuple(_ -> falses(channel_count, sample_count), 4)
    return InterleavedChopperRecorder(
        vector_traces[1:5]...,
        matrix_traces[1:2]...,
        boolean_traces...,
        vector_traces[6:9]...,
        matrix_traces[3:6]...,
        0,
    )
end

function _interleaved_linear_network(runtime)
    return _interleaved_chopper_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network
end

function _interleaved_bridge_states(runtime)
    network = _interleaved_linear_network(runtime)
    step = runtime.study.specification.timing.fixed_step_s
    return [Nonlinear.power_semiconductor_bridge_topology_state(bridge, network.v, step)
        for bridge in runtime.bridges]
end

function _interleaved_stored_energy(runtime)
    network = _interleaved_linear_network(runtime)
    output_voltage = network.v[runtime.output_node]
    passive_energy = 0.5 * runtime.study.capacitance_f * output_voltage^2 +
        sum(0.5 * runtime.study.inductance_h[channel] *
            runtime.inductors[channel].i_last^2 for channel in eachindex(runtime.inductors))
    return passive_energy + sum(state.stored_energy_j for state in
        _interleaved_bridge_states(runtime))
end

function _interleaved_linear_companion(runtime)
    network = _interleaved_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for channel in eachindex(runtime.inductors)
        branch = runtime.inductors[channel]
        voltage = Branches.branch_voltage(network.v, branch.a, branch.b)
        current = branch.i_last
        terminal_power += voltage * current
        dissipated_power += runtime.study.inductor_resistance_ohm[channel] * current^2
        stored_energy += 0.5 * runtime.study.inductance_h[channel] * current^2
    end
    capacitor_voltage = Branches.branch_voltage(
        network.v,
        runtime.capacitor.a,
        runtime.capacitor.b,
    )
    terminal_power += capacitor_voltage * runtime.capacitor.i_last
    stored_energy += 0.5 * runtime.study.capacitance_f * capacitor_voltage^2
    return (
        terminal_power_w=terminal_power,
        dissipated_power_w=dissipated_power,
        stored_energy_j=stored_energy,
    )
end

function _interleaved_device_companion_residual(runtime)
    _interleaved_chopper_is_detailed(runtime.study) || return 0.0
    return sum(state.companion_energy_residual_j for state in
        _interleaved_bridge_states(runtime))
end

function _interleaved_power(runtime)
    study = runtime.study
    network = _interleaved_linear_network(runtime)
    input_current = (study.input_voltage_v - network.v[runtime.input_positive_node]) /
        study.source_resistance_ohm
    output_voltage = network.v[runtime.output_node]
    inductor_loss = sum(study.inductor_resistance_ohm[channel] *
        runtime.inductors[channel].i_last^2 for channel in eachindex(runtime.inductors))
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w for bridge in runtime.bridges for valve in bridge.valves
    )
    series_loss = study.source_resistance_ohm * input_current^2 + inductor_loss
    return (
        input_w=study.input_voltage_v * input_current,
        load_w=output_voltage^2 / study.load_resistance_ohm,
        series_loss_w=series_loss,
        semiconductor_loss_w=semiconductor_loss,
        loss_w=series_loss + semiconductor_loss,
    )
end

function _observe_interleaved_substep!(
    accounting::InterleavedChopperStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _interleaved_power(runtime)
    final_linear = _interleaved_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    if companion_method === Branches.TrapezoidalCompanion
        input_energy = 0.5 * step_s * (previous_power.input_w + final_power.input_w)
        load_energy = 0.5 * step_s * (previous_power.load_w + final_power.load_w)
        series_energy = 0.5 * step_s *
            (previous_power.series_loss_w + final_power.series_loss_w)
        terminal_work = 0.5 * step_s *
            (previous_linear.terminal_power_w + final_linear.terminal_power_w)
        linear_loss = 0.5 * step_s *
            (previous_linear.dissipated_power_w + final_linear.dissipated_power_w)
    elseif companion_method === Branches.BackwardEulerCompanion
        input_energy = step_s * final_power.input_w
        load_energy = step_s * final_power.load_w
        series_energy = step_s * final_power.series_loss_w
        terminal_work = step_s * final_linear.terminal_power_w
        linear_loss = step_s * final_linear.dissipated_power_w
    else
        throw(ArgumentError(
            "interleaved energy accounting requires an accepted companion method",
        ))
    end
    accounting.input_energy_j += input_energy
    accounting.load_energy_j += load_energy
    accounting.series_dissipated_energy_j += series_energy
    semiconductor_energy = if companion_method === Branches.TrapezoidalCompanion
        0.5 * step_s * (
            previous_power.semiconductor_loss_w + final_power.semiconductor_loss_w
        )
    else
        step_s * final_power.semiconductor_loss_w
    end
    accounting.semiconductor_dissipated_energy_j += semiconductor_energy
    accounting.linear_companion_energy_residual_j += terminal_work - linear_loss -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _request_interleaved_gates!(runtime, time_s)
    requested = _interleaved_chopper_gate_state(runtime.study, time_s)
    changed = false
    for channel in eachindex(runtime.bridges)
        bridge = runtime.bridges[channel]
        previous = bridge.valves[1].gate_driver.commanded_on
        Nonlinear.request_power_semiconductor_topology_gates!(
            bridge,
            Bool[requested[channel], false],
            time_s,
        )
        if requested[channel] != previous
            changed = true
            push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                time_s,
                requested[channel] ? :channel_gate_turn_on : :channel_gate_turn_off,
                Symbol("channel_", channel),
                true;
                message="accepted phase-shifted carrier edge on the fixed-step calendar",
            ))
        end
    end
    return changed
end

function _synchronize_interleaved_commutation!(runtime, time_s)
    _interleaved_chopper_is_detailed(runtime.study) || return false
    changed = false
    for channel in eachindex(runtime.bridges)
        controlled, diode = runtime.bridges[channel].valves
        requested_on = controlled.gate_driver.applied_on
        if requested_on
            if diode.closed
                Nonlinear.apply_power_semiconductor_forward_extinction!(diode, time_s)
                changed = true
            end
            if !controlled.closed
                Nonlinear.apply_power_semiconductor_forward_turn_on!(controlled, time_s)
                changed = true
            end
        else
            if controlled.closed
                Nonlinear.apply_power_semiconductor_forward_extinction!(controlled, time_s)
                changed = true
            end
            should_freewheel = runtime.inductors[channel].i_last > diode.holding_current
            if should_freewheel && !diode.closed
                Nonlinear.apply_power_semiconductor_forward_turn_on!(diode, time_s)
                changed = true
            elseif !should_freewheel && diode.closed
                Nonlinear.apply_power_semiconductor_forward_extinction!(diode, time_s)
                changed = true
            end
        end
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :interleaved_channel_commutation,
        :interleaved_bridge,
        true;
        message="accepted detailed per-channel commutation with individual current continuity",
    ))
    return changed
end

function _stabilize_interleaved_topology!(
    runtime,
    time_s,
    step_s;
    maximum_iterations=16,
)
    _interleaved_chopper_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed interleaved execution must use the nonlinear D200 device solver",
    ))
    network = _interleaved_linear_network(runtime)
    for bridge in runtime.bridges
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s)
    end
    for iteration in 1:maximum_iterations
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        selected = nothing
        for (channel, bridge) in enumerate(runtime.bridges)
            for (position_index, valve) in enumerate(bridge.valves)
                action = _three_phase_vsc_topology_action(valve, network.v)
                action === nothing && continue
                selected = (channel, position_index, valve, action)
                break
            end
            selected === nothing || break
        end
        selected === nothing && return iteration
        channel, position_index, valve, action = selected
        _apply_three_phase_vsc_topology_action!(valve, action.transition, time_s)
        position_name = runtime.study.topologies[channel].valve_positions[position_index].name
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            Symbol("channel_", channel, "_", position_name),
            true;
            message="accepted interleaved switching-state commutation",
        ))
    end
    throw(ArgumentError("interleaved switching-state topology failed to stabilize"))
end

function _record_interleaved_sample!(recorder, runtime, sample_index, energy_residual_w)
    study = runtime.study
    network = _interleaved_linear_network(runtime)
    states = _interleaved_bridge_states(runtime)
    input_current = (study.input_voltage_v - network.v[runtime.input_positive_node]) /
        study.source_resistance_ohm
    output_voltage = network.v[runtime.output_node]
    recorder.time_s[sample_index] = runtime.time_s
    recorder.input_voltage_v[sample_index] = study.input_voltage_v
    recorder.input_current_a[sample_index] = input_current
    recorder.output_voltage_v[sample_index] = output_voltage
    recorder.load_current_a[sample_index] = output_voltage / study.load_resistance_ohm
    for channel in eachindex(runtime.bridges)
        state = states[channel]
        controlled, diode = runtime.bridges[channel].valves
        recorder.channel_switch_voltage_v[channel, sample_index] =
            network.v[runtime.switching_nodes[channel]]
        recorder.channel_inductor_current_a[channel, sample_index] =
            runtime.inductors[channel].i_last
        recorder.requested_gate_state[channel, sample_index] =
            state.requested_gate_state[1]
        recorder.applied_gate_state[channel, sample_index] = state.applied_gate_state[1]
        recorder.controlled_conducting_state[channel, sample_index] =
            state.conducting_state[1]
        recorder.diode_conducting_state[channel, sample_index] = state.conducting_state[2]
        if _interleaved_chopper_is_detailed(study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            diode_state = Nonlinear.power_semiconductor_extended_state(diode)
            recorder.controlled_junction_temperature_k[channel, sample_index] =
                controlled_state.junction_temperature_k
            recorder.diode_junction_temperature_k[channel, sample_index] =
                diode_state.junction_temperature_k
            recorder.diode_recovered_charge_c[channel, sample_index] =
                diode_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[channel, sample_index] =
                controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[channel, sample_index] = 0.0
            recorder.diode_junction_temperature_k[channel, sample_index] = 0.0
            recorder.diode_recovered_charge_c[channel, sample_index] = 0.0
            recorder.controlled_turn_off_tail_current_a[channel, sample_index] = 0.0
        end
    end
    recorder.stored_energy_j[sample_index] = _interleaved_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample_index] =
        sum(state.semiconductor_loss_w for state in states)
    recorder.kcl_residual_a[sample_index] = if _interleaved_chopper_is_detailed(study)
        runtime.last_kcl_residual_a
    else
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    end
    recorder.energy_residual_w[sample_index] = energy_residual_w
    recorder.write_index = sample_index
    return recorder
end

function prepare_switching_interleaved_chopper(study::SwitchingInterleavedChopperStudy)
    runtime = _seed_interleaved_chopper(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    ) + 1
    channel_count = study.specification.selection.channel_count
    integrator = InterleavedChopperIntegrator(
        runtime,
        TimestepTransaction(runtime),
        InterleavedChopperRecorder(channel_count, sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_interleaved_sample!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_interleaved_chopper!(integrator)
    integrator.failed && throw(ArgumentError(
        "interleaved chopper integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _interleaved_power(runtime)
    previous_linear = _interleaved_linear_companion(runtime)
    previous_device_residual = _interleaved_device_companion_residual(runtime)
    previous_semiconductor_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for bridge in runtime.bridges for valve in bridge.valves
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_interleaved_gates!(runtime, endpoint)
        accounting = InterleavedChopperStepEnergyAccounting(previous_power, previous_linear)
        semiconductor_increment = 0.0
        if _interleaved_chopper_is_detailed(study)
            commutation_changed = _synchronize_interleaved_commutation!(runtime, endpoint)
            for bridge in runtime.bridges
                Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, endpoint)
            end
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || commutation_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || commutation_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) -> _observe_interleaved_substep!(
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
                for bridge in runtime.bridges for valve in bridge.valves
            )
            semiconductor_increment =
                final_semiconductor_dissipation - previous_semiconductor_dissipation
        else
            event_count_before = length(runtime.events)
            _stabilize_interleaved_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count_before
            network = _interleaved_linear_network(runtime)
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
                _observe_interleaved_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_interleaved_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_interleaved_substep!(
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
                _observe_interleaved_substep!(
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
        dissipated_increment = accounting.series_dissipated_energy_j +
            semiconductor_increment
        runtime.input_energy_j += accounting.input_energy_j
        runtime.load_energy_j += accounting.load_energy_j
        runtime.dissipated_energy_j += dissipated_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        stored_energy = _interleaved_stored_energy(runtime)
        device_residual_increment =
            _interleaved_device_companion_residual(runtime) - previous_device_residual
        energy_residual = (
            accounting.input_energy_j - accounting.load_energy_j -
            dissipated_increment - (stored_energy - previous_energy) -
            accounting.linear_companion_energy_residual_j - device_residual_increment
        ) / step
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_interleaved_sample!(
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

function _interleaved_chopper_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _interleaved_linear_network(runtime)
    states = _interleaved_bridge_states(runtime)
    output_voltage = network.v[runtime.output_node]
    input_current = (study.input_voltage_v - network.v[runtime.input_positive_node]) /
        study.source_resistance_ohm
    stored_energy = _interleaved_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        (state.deterministic_signature for state in states)...,
        repr(runtime.time_s),
        repr(getfield.(runtime.inductors, :i_last)),
        repr(output_voltage),
        repr(stored_energy),
        repr(runtime.dissipated_energy_j),
        string(integrator.accepted_step_index),
        string(length(runtime.events)),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, output_voltage],
        [input_current, output_voltage / study.load_resistance_ohm],
        BitVector(Iterators.flatten(getfield.(states, :requested_gate_state))),
        BitVector(Iterators.flatten(getfield.(states, :applied_gate_state))),
        BitVector(Iterators.flatten(getfield.(states, :conducting_state))),
        [study.capacitance_f * output_voltage],
        study.inductance_h .* getfield.(runtime.inductors, :i_last),
        Float64[],
        vcat(
            study.specification.modulation.duty,
            ConverterSystems.interleaved_carrier_phases_rad(
                study.specification.selection.channel_count,
            ),
        ),
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
        _interleaved_device_companion_residual(runtime)
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
    return ConverterSystems.SwitchingInterleavedChopperTrace(
        recorder.time_s,
        recorder.input_voltage_v,
        recorder.input_current_a,
        recorder.output_voltage_v,
        recorder.load_current_a,
        recorder.channel_switch_voltage_v,
        recorder.channel_inductor_current_a,
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

function execute_switching_interleaved_chopper!(integrator::InterleavedChopperIntegrator)
    while _advance_interleaved_chopper!(integrator)
    end
    return _interleaved_chopper_result(integrator)
end

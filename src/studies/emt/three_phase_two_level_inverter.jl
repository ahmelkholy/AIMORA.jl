mutable struct AverageThreePhaseTwoLevelInverterRuntime{S}
    study::S
    phase_current_a::Vector{Float64}
    time_s::Float64
    accepted_step_index::Int
    input_energy_j::Float64
    load_energy_j::Float64
    transition_factor::Float64
    source_gain_a_per_v::Float64
    time_trace_s::Vector{Float64}
    input_voltage_trace_v::Vector{Float64}
    input_current_trace_a::Vector{Float64}
    phase_voltage_trace_v::Matrix{Float64}
    phase_current_trace_a::Matrix{Float64}
    stored_energy_trace_j::Vector{Float64}
    circuit_residual_trace_v::Matrix{Float64}
    energy_residual_trace_w::Vector{Float64}
end

function _three_phase_inverter_phase_reference(study, time_s)
    modulation = study.specification.modulation
    angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
        modulation.phase_shift_rad
    amplitude = 0.5 * modulation.modulation_index * study.input_voltage_v
    return ntuple(phase ->
        amplitude * sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
end

function prepare_average_three_phase_two_level_inverter(
    study::ConverterSystems.AverageThreePhaseTwoLevelInverterStudy,
)
    step = study.specification.timing.fixed_step_s
    rate = -study.load_resistance_ohm / study.load_inductance_h
    left = 1.0 - 0.5 * step * rate
    transition = (1.0 + 0.5 * step * rate) / left
    source_gain = 0.5 * step / study.load_inductance_h / left
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / step,
    ) + 1
    runtime = AverageThreePhaseTwoLevelInverterRuntime(
        study,
        collect(study.initial_state.phase_current_a),
        study.start_time_s,
        0,
        0.0,
        0.0,
        transition,
        source_gain,
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, sample_count),
    )
    _record_average_three_phase_inverter!(runtime, 1, zeros(3), 0.0)
    return runtime
end

function _record_average_three_phase_inverter!(
    runtime,
    sample,
    circuit_residual_v,
    energy_residual_w,
)
    study = runtime.study
    phase_voltage = _three_phase_inverter_phase_reference(study, runtime.time_s)
    current = runtime.phase_current_a
    phase_power = sum(phase_voltage[index] * current[index] for index in 1:3)
    runtime.time_trace_s[sample] = runtime.time_s
    runtime.input_voltage_trace_v[sample] = study.input_voltage_v
    runtime.input_current_trace_a[sample] = phase_power / study.input_voltage_v
    runtime.phase_voltage_trace_v[:, sample] .= phase_voltage
    runtime.phase_current_trace_a[:, sample] .= current
    runtime.stored_energy_trace_j[sample] =
        0.5 * study.load_inductance_h * sum(abs2, current)
    runtime.circuit_residual_trace_v[:, sample] .= circuit_residual_v
    runtime.energy_residual_trace_w[sample] = energy_residual_w
    return runtime
end

function _advance_average_three_phase_inverter!(runtime, sample)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    previous_time = runtime.time_s
    endpoint = study.start_time_s + (sample - 1) * step
    previous_voltage = _three_phase_inverter_phase_reference(study, previous_time)
    endpoint_voltage = _three_phase_inverter_phase_reference(study, endpoint)
    previous_current = copy(runtime.phase_current_a)
    previous_stored_energy =
        0.5 * study.load_inductance_h * sum(abs2, previous_current)
    for phase in 1:3
        runtime.phase_current_a[phase] =
            runtime.transition_factor * previous_current[phase] +
            runtime.source_gain_a_per_v *
            (previous_voltage[phase] + endpoint_voltage[phase])
    end
    current = runtime.phase_current_a
    average_current = 0.5 .* (previous_current .+ current)
    average_voltage = 0.5 .* (collect(previous_voltage) .+ collect(endpoint_voltage))
    derivative = (current .- previous_current) ./ step
    circuit_residual = average_voltage .-
        study.load_resistance_ohm .* average_current .-
        study.load_inductance_h .* derivative
    input_energy = step * dot(average_voltage, average_current)
    load_energy = step * study.load_resistance_ohm * sum(abs2, average_current)
    stored_energy = 0.5 * study.load_inductance_h * sum(abs2, current)
    energy_residual = (input_energy - load_energy -
        (stored_energy - previous_stored_energy)) / step
    runtime.input_energy_j += input_energy
    runtime.load_energy_j += load_energy
    runtime.time_s = endpoint
    return _record_average_three_phase_inverter!(
        runtime,
        sample,
        circuit_residual,
        energy_residual,
    )
end

function _average_three_phase_inverter_result(runtime)
    study = runtime.study
    phase_voltage = _three_phase_inverter_phase_reference(study, runtime.time_s)
    phase_current = runtime.phase_current_a
    input_current = dot(phase_voltage, phase_current) / study.input_voltage_v
    stored_energy = 0.5 * study.load_inductance_h * sum(abs2, phase_current)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(phase_voltage),
        repr(phase_current),
        repr(stored_energy),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [study.input_voltage_v, phase_voltage...],
        [input_current, phase_current...],
        BitVector(),
        BitVector(),
        BitVector(),
        Float64[],
        study.load_inductance_h .* phase_current,
        Float64[],
        [study.specification.modulation.modulation_index],
        stored_energy,
        runtime.load_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
    voltage_scale = max(maximum(abs, runtime.phase_voltage_trace_v), 1.0)
    power_scale = max(maximum(abs,
        runtime.input_voltage_trace_v .* runtime.input_current_trace_a), 1.0)
    return ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=maximum(abs,
            vec(sum(runtime.phase_current_trace_a; dims=1))),
        relative_charge_residual=maximum(abs, runtime.circuit_residual_trace_v) /
            voltage_scale,
        relative_energy_residual=maximum(abs, runtime.energy_residual_trace_w) /
            power_scale,
        harmonic_metrics=Dict(:phase_voltage_thd => 0.0),
    )
end

function execute_average_three_phase_two_level_inverter!(
    runtime::AverageThreePhaseTwoLevelInverterRuntime,
)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    return ConverterSystems.AverageThreePhaseTwoLevelInverterTrace(
        runtime.time_trace_s,
        runtime.input_voltage_trace_v,
        runtime.input_current_trace_a,
        runtime.phase_voltage_trace_v,
        runtime.phase_current_trace_a,
        runtime.stored_energy_trace_j,
        runtime.circuit_residual_trace_v,
        runtime.energy_residual_trace_w,
        _average_three_phase_inverter_result(runtime),
    )
end

mutable struct ThreePhaseTwoLevelInverterRuntime{S,N,T,L,V,D,B}
    study::S
    network::N
    bridge::T
    legs::L
    controlled_valves::V
    freewheel_diodes::D
    load_branches::B
    dc_positive_node::Int
    dc_negative_node::Int
    phase_nodes::NTuple{3,Int}
    load_neutral_node::Int
    time_s::Float64
    input_energy_j::Float64
    load_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct ThreePhaseTwoLevelInverterRecorder
    time_s::Vector{Float64}
    dc_link_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
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

mutable struct ThreePhaseTwoLevelInverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::ThreePhaseTwoLevelInverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_three_phase_inverter_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_three_phase_inverter_linear_network(runtime) =
    _three_phase_inverter_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network

function _three_phase_inverter_nodes(study)
    nodes = Dict(node.name => node.node for node in study.topology.nodes)
    required = (:dc_positive, :dc_negative, :ac_1, :ac_2, :ac_3)
    all(haskey(nodes, name) for name in required) || throw(ArgumentError(
        "three-phase two-level bridge is missing a canonical terminal",
    ))
    input = only(port for port in study.specification.ports if port.identity === :input_dc)
    output = only(port for port in study.specification.ports if port.identity === :output_ac)
    input.ordered_nodes == (nodes[:dc_positive], nodes[:dc_negative]) ||
        throw(ArgumentError("three-phase inverter DC port does not match its rails"))
    output.ordered_nodes == (nodes[:ac_1], nodes[:ac_2], nodes[:ac_3]) ||
        throw(ArgumentError("three-phase inverter AC port does not match its poles"))
    phase_nodes = (nodes[:ac_1], nodes[:ac_2], nodes[:ac_3])
    neutral = maximum((nodes[:dc_positive], nodes[:dc_negative], phase_nodes...)) + 1
    return (
        dc_positive=nodes[:dc_positive],
        dc_negative=nodes[:dc_negative],
        phase=phase_nodes,
        neutral,
    )
end

function _three_phase_inverter_duties(study, time_s)
    modulation = study.specification.modulation
    angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
        modulation.phase_shift_rad
    phase_reference = ntuple(phase ->
        0.5 * study.input_voltage_v * modulation.modulation_index *
        sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
    return ConverterSystems.converter_two_level_duties(
        phase_reference,
        study.input_voltage_v,
        modulation.kind,
    )
end

function _three_phase_inverter_gate_state(study, time_s)
    return ConverterSystems.converter_pwm_gate_state(
        _three_phase_inverter_duties(study, time_s),
        time_s,
        study.specification.timing.carrier_frequency_hz,
        study.specification.modulation.kind,
    ).requested_valve_state
end

function _three_phase_inverter_initial_paths(gates, phase_current)
    controlled = falses(6)
    diodes = falses(6)
    for phase in 1:3
        upper = 2 * phase - 1
        lower = 2 * phase
        current = phase_current[phase]
        if current >= 0.0
            gates[upper] ? (controlled[upper] = true) : (diodes[lower] = true)
        else
            gates[lower] ? (controlled[lower] = true) : (diodes[upper] = true)
        end
    end
    return controlled, diodes
end

function _seed_three_phase_two_level_inverter(study)
    nodes = _three_phase_inverter_nodes(study)
    gates = _three_phase_inverter_gate_state(study, study.start_time_s)
    initial_current = study.initial_state.phase_current_a
    controlled_paths, diode_paths =
        _three_phase_inverter_initial_paths(gates, initial_current)
    controlled = ntuple(6) do index
        _full_bridge_controlled_valve(
            study,
            study.topology.valve_positions[index],
            gates[index],
            controlled_paths[index],
        )
    end
    diodes = ntuple(6) do index
        _full_bridge_freewheel_diode(
            study,
            study.topology.valve_positions[index],
            diode_paths[index],
        )
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, controlled)
    foreach(Nonlinear.power_semiconductor_event_localization!, diodes)
    legs = ntuple(3) do phase
        Nonlinear.PowerSemiconductorBridgeLeg(
            controlled[2 * phase - 1],
            controlled[2 * phase];
            commutation_dead_time_s=study.specification.timing.dead_time_s,
        )
    end
    bridge = Nonlinear.PowerSemiconductorBridgeTopology(
        study.topology,
        collect(controlled);
        bridge_legs=collect(legs),
    )
    dc_voltage = study.input_voltage_v
    pole_voltage = ntuple(phase -> gates[2 * phase - 1] ? dc_voltage : 0.0, 3)
    neutral_voltage = sum(pole_voltage) / 3.0
    phase_voltage = ntuple(phase -> pole_voltage[phase] - neutral_voltage, 3)
    phase_power = sum(phase_voltage[phase] * initial_current[phase] for phase in 1:3)
    source_current = phase_power / dc_voltage
    dc_link_voltage = study.input_voltage_v -
        study.source_resistance_ohm * source_current
    source = Branches.TwoTerminalTheveninSource(
        nodes.dc_positive,
        nodes.dc_negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    loads = ntuple(3) do phase
        Branches.SeriesRLBranch(
            nodes.phase[phase],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            initial_current[phase],
            phase_voltage[phase],
            initial_current[phase],
        )
    end
    node_count = maximum((nodes.dc_positive, nodes.phase..., nodes.neutral))
    elements = _three_phase_inverter_is_detailed(study) ?
        Any[source, bridge, loads...] : Any[source, bridge, diodes..., loads...]
    linear_network = Nodal.NodalSystem(node_count, elements)
    linear_network.v[nodes.dc_positive] = dc_link_voltage
    nodes.dc_negative == 0 || (linear_network.v[nodes.dc_negative] = 0.0)
    for phase in 1:3
        linear_network.v[nodes.phase[phase]] =
            (gates[2 * phase - 1] ? dc_link_voltage : 0.0)
    end
    linear_network.v[nodes.neutral] =
        sum(linear_network.v[node] for node in nodes.phase) / 3.0
    if _three_phase_inverter_is_detailed(study)
        for valve in (controlled..., diodes...)
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _three_phase_inverter_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        (controlled..., diodes...);
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=max(maximum(abs, initial_current), 1.0),
        ),
    ) : linear_network
    return ThreePhaseTwoLevelInverterRuntime(
        study,
        network,
        bridge,
        legs,
        controlled,
        diodes,
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
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function ThreePhaseTwoLevelInverterRecorder(sample_count)
    sample_count > 1 || throw(ArgumentError(
        "three-phase inverter recorder requires at least two samples",
    ))
    return ThreePhaseTwoLevelInverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        falses(6, sample_count),
        falses(6, sample_count),
        falses(6, sample_count),
        falses(6, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 6, sample_count),
        zeros(Float64, 6, sample_count),
        zeros(Float64, 6, sample_count),
        zeros(Float64, 6, sample_count),
        0,
    )
end

_three_phase_inverter_devices(runtime) =
    (runtime.controlled_valves..., runtime.freewheel_diodes...)

function _three_phase_inverter_stored_energy(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    _three_phase_inverter_is_detailed(runtime.study) || return energy
    return energy + sum(_three_phase_inverter_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _three_phase_inverter_power(runtime)
    study = runtime.study
    network = _three_phase_inverter_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    input_current = (study.input_voltage_v - dc_voltage) /
        study.source_resistance_ohm
    load_loss = study.load_resistance_ohm *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _three_phase_inverter_devices(runtime)
    )
    return (
        input_w=study.input_voltage_v * input_current,
        load_w=load_loss,
        source_loss_w=study.source_resistance_ohm * input_current^2,
        semiconductor_loss_w=semiconductor_loss,
        input_current_a=input_current,
        dc_voltage_v=dc_voltage,
    )
end

function _three_phase_inverter_linear_companion(runtime)
    network = _three_phase_inverter_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for branch in runtime.load_branches
        voltage = Branches.branch_voltage(network.v, branch.a, branch.b)
        current = branch.i_last
        terminal_power += voltage * current
        dissipated_power += runtime.study.load_resistance_ohm * current^2
        stored_energy += 0.5 * runtime.study.load_inductance_h * current^2
    end
    return (
        terminal_power_w=terminal_power,
        dissipated_power_w=dissipated_power,
        stored_energy_j=stored_energy,
    )
end

function _three_phase_inverter_device_companion_residual(runtime)
    _three_phase_inverter_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _three_phase_inverter_devices(runtime)
    )
end

function _observe_three_phase_inverter_substep!(
    accounting::FullBridgeStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _three_phase_inverter_power(runtime)
    final_linear = _three_phase_inverter_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("three-phase inverter received an unknown companion method"))
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

function _stabilize_three_phase_inverter_topology!(runtime, time_s, step_s)
    _three_phase_inverter_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed three-phase inverter must use the nonlinear D200 solver",
    ))
    network = _three_phase_inverter_linear_network(runtime)
    Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(runtime.bridge, time_s)
    devices = _three_phase_inverter_devices(runtime)
    for _iteration in 1:48
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        selected = nothing
        for (index, device) in enumerate(devices)
            action = _three_phase_vsc_topology_action(device, network.v)
            action === nothing && continue
            selected = (index, device, action)
            break
        end
        selected === nothing && return nothing
        index, device, action = selected
        _apply_three_phase_vsc_topology_action!(device, action.transition, time_s)
        controlled = index <= 6
        position = controlled ? index : index - 6
        owner = runtime.study.topology.valve_positions[position].name
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            controlled ? owner : Symbol(owner, :_freewheel),
            true;
            message="accepted three-phase inverter switching-state commutation",
        ))
    end
    throw(ArgumentError("three-phase inverter switching topology failed to stabilize"))
end

function _synchronize_three_phase_inverter_conduction!(runtime, time_s)
    controlled_target = falses(6)
    diode_target = falses(6)
    tolerance = maximum(
        valve.holding_current for valve in _three_phase_inverter_devices(runtime)
    )
    for phase in 1:3
        upper = 2 * phase - 1
        lower = 2 * phase
        current = runtime.load_branches[phase].i_last
        upper_applied = runtime.controlled_valves[upper].gate_driver.applied_on
        lower_applied = runtime.controlled_valves[lower].gate_driver.applied_on
        if upper_applied && current >= -tolerance
            controlled_target[upper] = true
        elseif lower_applied && current <= tolerance
            controlled_target[lower] = true
        elseif current > tolerance
            diode_target[lower] = true
        elseif current < -tolerance
            diode_target[upper] = true
        elseif upper_applied
            controlled_target[upper] = true
        elseif lower_applied
            controlled_target[lower] = true
        end
    end
    changed = false
    for index in 1:6
        changed |= _set_full_bridge_device_state!(
            runtime.controlled_valves[index],
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

function _request_three_phase_inverter_gates!(runtime, time_s)
    requested = _three_phase_inverter_gate_state(runtime.study, time_s)
    changed = false
    for phase in 1:3
        previous = runtime.legs[phase].requested_upper_on
        upper = requested[2 * phase - 1]
        Nonlinear.request_power_semiconductor_bridge_pole!(
            runtime.legs[phase],
            upper,
            time_s,
        )
        changed |= upper != previous
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :three_phase_pwm_command,
        :three_phase_two_level_bridge,
        true;
        message="accepted balanced carrier-PWM gate command",
    ))
    return changed
end

function _record_three_phase_inverter!(recorder, runtime, sample, energy_residual_w)
    network = _three_phase_inverter_linear_network(runtime)
    power = _three_phase_inverter_power(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    recorder.time_s[sample] = runtime.time_s
    recorder.dc_link_voltage_v[sample] = power.dc_voltage_v
    recorder.input_current_a[sample] = power.input_current_a
    for phase in 1:3
        recorder.phase_voltage_v[phase, sample] =
            network.v[runtime.phase_nodes[phase]] - neutral_voltage
        recorder.phase_current_a[phase, sample] = runtime.load_branches[phase].i_last
    end
    for index in 1:6
        controlled = runtime.controlled_valves[index]
        diode = runtime.freewheel_diodes[index]
        recorder.requested_gate_state[index, sample] =
            controlled.gate_driver.commanded_on
        recorder.applied_gate_state[index, sample] = controlled.gate_driver.applied_on
        recorder.controlled_conducting_state[index, sample] = controlled.closed
        recorder.diode_conducting_state[index, sample] = diode.closed
        if _three_phase_inverter_is_detailed(runtime.study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            diode_state = Nonlinear.power_semiconductor_extended_state(diode)
            recorder.controlled_junction_temperature_k[index, sample] =
                controlled_state.junction_temperature_k
            recorder.diode_junction_temperature_k[index, sample] =
                diode_state.junction_temperature_k
            recorder.diode_recovered_charge_c[index, sample] =
                diode_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[index, sample] =
                controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[index, sample] = 0.0
            recorder.diode_junction_temperature_k[index, sample] = 0.0
            recorder.diode_recovered_charge_c[index, sample] = 0.0
            recorder.controlled_turn_off_tail_current_a[index, sample] = 0.0
        end
    end
    recorder.stored_energy_j[sample] = _three_phase_inverter_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _three_phase_inverter_is_detailed(runtime.study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual_w
    recorder.write_index = sample
    return recorder
end

function prepare_switching_three_phase_two_level_inverter(
    study::ConverterSystems.SwitchingThreePhaseTwoLevelInverterStudy,
)
    runtime = _seed_three_phase_two_level_inverter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) /
        study.specification.timing.fixed_step_s,
    ) + 1
    integrator = ThreePhaseTwoLevelInverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        ThreePhaseTwoLevelInverterRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_three_phase_inverter!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_three_phase_two_level_inverter!(integrator)
    integrator.failed && throw(ArgumentError(
        "three-phase inverter integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _three_phase_inverter_power(runtime)
    previous_linear = _three_phase_inverter_linear_companion(runtime)
    previous_stored = _three_phase_inverter_stored_energy(runtime)
    previous_device_residual = _three_phase_inverter_device_companion_residual(runtime)
    previous_device_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _three_phase_inverter_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_three_phase_inverter_gates!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        device_increment = 0.0
        if _three_phase_inverter_is_detailed(study)
            Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                runtime.bridge,
                endpoint,
            )
            conduction_changed =
                _synchronize_three_phase_inverter_conduction!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) ->
                    _observe_three_phase_inverter_substep!(
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
            final_device_dissipation = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in _three_phase_inverter_devices(runtime)
            )
            device_increment = final_device_dissipation - previous_device_dissipation
        else
            event_count = length(runtime.events)
            _stabilize_three_phase_inverter_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count
            network = _three_phase_inverter_linear_network(runtime)
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
                _observe_three_phase_inverter_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_three_phase_inverter_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_three_phase_inverter_substep!(
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
                _observe_three_phase_inverter_substep!(
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
            device_increment = accounting.semiconductor_dissipated_energy_j
        end
        final_stored = _three_phase_inverter_stored_energy(runtime)
        device_residual_increment =
            _three_phase_inverter_device_companion_residual(runtime) -
            previous_device_residual
        energy_residual = (
            accounting.input_energy_j - accounting.load_energy_j -
            accounting.source_dissipated_energy_j - device_increment -
            (final_stored - previous_stored) -
            accounting.linear_companion_energy_residual_j -
            device_residual_increment
        ) / step
        runtime.input_energy_j += accounting.input_energy_j
        runtime.load_energy_j += accounting.load_energy_j
        runtime.dissipated_energy_j +=
            accounting.load_energy_j + accounting.source_dissipated_energy_j +
            device_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_three_phase_inverter!(
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

function _three_phase_two_level_inverter_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _three_phase_inverter_linear_network(runtime)
    power = _three_phase_inverter_power(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    phase_voltage = [
        network.v[node] - neutral_voltage for node in runtime.phase_nodes
    ]
    phase_current = [branch.i_last for branch in runtime.load_branches]
    stored_energy = _three_phase_inverter_stored_energy(runtime)
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        study.specification.timing.fixed_step_s,
    )
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        bridge_state.deterministic_signature,
        repr(runtime.time_s),
        repr(phase_voltage),
        repr(phase_current),
        repr(stored_energy),
        string(integrator.accepted_step_index),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [power.dc_voltage_v, phase_voltage...],
        [power.input_current_a, phase_current...],
        bridge_state.requested_gate_state,
        bridge_state.applied_gate_state,
        BitVector((
            getfield.(runtime.controlled_valves, :closed)...,
            getfield.(runtime.freewheel_diodes, :closed)...,
        )),
        Float64[],
        study.load_inductance_h .* phase_current,
        Float64[],
        collect(_three_phase_inverter_duties(study, runtime.time_s)),
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    energy_scale = max(abs(runtime.input_energy_j), abs(runtime.load_energy_j),
        abs(runtime.dissipated_energy_j), abs(stored_energy), eps(Float64))
    integrated_residual = runtime.input_energy_j - runtime.dissipated_energy_j -
        (stored_energy - recorder.stored_energy_j[1]) -
        runtime.linear_companion_energy_residual_j -
        _three_phase_inverter_device_companion_residual(runtime)
    result = ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=integrator.completed && !integrator.failed,
        status=integrator.completed && !integrator.failed ? :ok : :incomplete,
        events=runtime.events,
        maximum_kcl_residual_a=maximum(recorder.kcl_residual_a),
        relative_charge_residual=maximum(abs,
            vec(sum(recorder.phase_current_a; dims=1))),
        relative_energy_residual=abs(integrated_residual) / energy_scale,
    )
    return ConverterSystems.SwitchingThreePhaseTwoLevelInverterTrace(
        recorder.time_s,
        recorder.dc_link_voltage_v,
        recorder.input_current_a,
        recorder.phase_voltage_v,
        recorder.phase_current_a,
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

function execute_switching_three_phase_two_level_inverter!(
    integrator::ThreePhaseTwoLevelInverterIntegrator,
)
    while _advance_three_phase_two_level_inverter!(integrator)
    end
    return _three_phase_two_level_inverter_result(integrator)
end

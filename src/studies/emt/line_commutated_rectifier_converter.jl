mutable struct LineCommutatedRectifierRuntime{S,N,B,V,L}
    study::S
    network::N
    bridge::B
    valves::V
    load_branch::L
    ac_nodes::Vector{Int}
    source_neutral_nodes::Vector{Int}
    dc_positive_node::Int
    dc_negative_node::Int
    time_s::Float64
    source_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct LineCommutatedRectifierRecorder
    time_s::Vector{Float64}
    source_voltage_v::Matrix{Float64}
    source_current_a::Matrix{Float64}
    dc_voltage_v::Vector{Float64}
    dc_current_a::Vector{Float64}
    requested_firing_state::BitMatrix
    applied_firing_state::BitMatrix
    conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    junction_temperature_k::Matrix{Float64}
    recovered_charge_c::Matrix{Float64}
    turn_off_tail_current_a::Matrix{Float64}
    write_index::Int
end

mutable struct LineCommutatedRectifierIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::LineCommutatedRectifierRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

mutable struct LineCommutatedRectifierStepEnergyAccounting{P,L}
    previous_power::P
    previous_linear_companion::L
    source_energy_j::Float64
    source_dissipated_energy_j::Float64
    load_dissipated_energy_j::Float64
    semiconductor_dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
end

LineCommutatedRectifierStepEnergyAccounting(previous_power, previous_linear_companion) =
    LineCommutatedRectifierStepEnergyAccounting(
        previous_power,
        previous_linear_companion,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )

_rectifier_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_rectifier_is_natural_diode(study) = study.specification.selection.family in (
    ConverterSystems.SinglePhaseDiodeBridge,
    ConverterSystems.ThreePhaseDiodeBridge,
    ConverterSystems.MultipulseDiodeBridge,
)

_rectifier_linear_network(runtime) = _rectifier_is_detailed(runtime.study) ?
    nonlinear_linear_system(runtime.network) : runtime.network

function _line_commutated_rectifier_nodes(study)
    names = Dict(node.name => node.node for node in study.topology.nodes)
    phase_count = study.specification.selection.phase_count
    multipulse = study.specification.selection.family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    )
    ac_nodes = if multipulse
        group_count = study.specification.selection.pulse_count ÷ 6
        [names[Symbol(:group_, group, :_ac_, phase)]
            for group in 1:group_count for phase in 1:3]
    else
        [names[Symbol(:ac_, phase)] for phase in 1:(phase_count == 1 ? 2 : 3)]
    end
    dc_positive = names[:dc_positive]
    dc_negative = names[:dc_negative]
    ac_port = only(port for port in study.specification.ports if
        port.kind === ConverterSystems.AlternatingCurrentPort)
    dc_port = only(port for port in study.specification.ports if
        port.kind === ConverterSystems.DirectCurrentPort)
    ac_port.ordered_nodes == Tuple(ac_nodes) || throw(ArgumentError(
        "natural rectifier AC port does not match its bridge phase terminals",
    ))
    dc_port.ordered_nodes == (dc_positive, dc_negative) || throw(ArgumentError(
        "natural rectifier DC port does not match its bridge rails",
    ))
    return ac_nodes, dc_positive, dc_negative
end

function _rectifier_source_voltage(study, channel::Int, time_s::Real)
    group = study.specification.selection.family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    ) ? (channel - 1) ÷ 3 + 1 : 1
    phase = study.specification.selection.phase_count == 1 ? channel :
        (channel - 1) % 3 + 1
    angle = 2pi * study.source_frequency_hz * Float64(time_s) +
        study.group_phase_shift_rad[group]
    if study.specification.selection.phase_count == 1
        phase == 1 || throw(BoundsError())
        return study.source_peak_voltage_v * sin(angle)
    end
    1 <= phase <= 3 || throw(BoundsError())
    return study.source_peak_voltage_v * sin(angle - 2pi * (phase - 1) / 3)
end

function _rectifier_source_values(study, time_s)
    if study.specification.selection.phase_count == 1
        line_voltage = _rectifier_source_voltage(study, 1, time_s)
        return [0.5 * line_voltage, -0.5 * line_voltage]
    end
    channel_count = study.specification.selection.family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    ) ? 3 * length(study.group_phase_shift_rad) : 3
    return [_rectifier_source_voltage(study, channel, time_s)
        for channel in 1:channel_count]
end

function _rectifier_valve(study, position, initially_closed)
    detailed = _rectifier_is_detailed(study)
    parameters = study.detailed_semiconductor
    snubber = detailed ? Nonlinear.SeriesRCSnubber(
        parameters.snubber_resistance_ohm,
        parameters.snubber_capacitance_f,
    ) : nothing
    if position.valve_class === :diode
        return Nonlinear.DiodeValveSwitch(
            position.from_node,
            position.to_node;
            on_conductance=detailed ? parameters.freewheel_on_conductance_s : 1.0e6,
            off_conductance=detailed ? parameters.freewheel_off_conductance_s : 1.0e-9,
            forward_voltage_drop_v=detailed ? parameters.freewheel_forward_voltage_v : 0.0,
            snubber,
            extended_fidelity=detailed ?
                _detailed_chopper_freewheel_fidelity(parameters) : nothing,
            initially_closed,
        )
    end
    thyristor_fidelity = detailed ? Nonlinear.PowerSemiconductorExtendedFidelity(
        junction_charge=Nonlinear.NonlinearJunctionChargeFidelity(
            parameters.controlled_junction_capacitance_f,
            parameters.junction_potential_v,
            parameters.junction_grading_exponent;
            voltage_domain_v=(-parameters.junction_voltage_limit_v,
                parameters.junction_voltage_limit_v),
            provenance=parameters.provenance,
        ),
        switching_energy=_detailed_chopper_switching_energy_table(
            parameters;
            controlled=true,
        ),
        thermal=_detailed_chopper_thermal_fidelity(parameters),
        provenance=parameters.provenance,
    ) : nothing
    return Nonlinear.ThyristorValveSwitch(
        position.from_node,
        position.to_node;
        gate_driver=Nonlinear.PowerSemiconductorGateDriver(
            minimum_pulse_width_s=study.specification.timing.minimum_pulse_s,
            initially_on=initially_closed,
        ),
        on_conductance=detailed ? parameters.controlled_on_conductance_s : 1.0e6,
        off_conductance=detailed ? parameters.controlled_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? parameters.controlled_forward_voltage_v : 0.0,
        snubber,
        extended_fidelity=thyristor_fidelity,
        initially_closed,
    )
end

function _rectifier_conduction_path(study, time_s)
    source_values = _rectifier_source_values(study, time_s)
    family = study.specification.selection.family
    diode_bridge = family in (
        ConverterSystems.SinglePhaseDiodeBridge,
        ConverterSystems.ThreePhaseDiodeBridge,
        ConverterSystems.MultipulseDiodeBridge,
    )
    half_controlled = family in (
        ConverterSystems.SinglePhaseHalfControlledBridge,
        ConverterSystems.ThreePhaseHalfControlledBridge,
    )
    delayed_values = if diode_bridge
        source_values
    else
        delay_s = study.specification.modulation.firing_angle_rad /
            (2pi * study.source_frequency_hz)
        _rectifier_source_values(study, time_s - delay_s)
    end
    multipulse = family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    )
    if multipulse
        group_count = length(study.group_phase_shift_rad)
        group_values = half_controlled ? source_values : delayed_values
        group_index = argmax([
            maximum(view(group_values, (3group - 2):(3group))) -
                minimum(view(group_values, (3group - 2):(3group)))
            for group in 1:group_count
        ])
        offset = 3 * (group_index - 1)
        upper_phase = offset + argmax(view(delayed_values, (offset + 1):(offset + 3)))
        lower_phase = offset + argmin(view(group_values, (offset + 1):(offset + 3)))
    else
        upper_phase = argmax(delayed_values)
        lower_phase = argmin(half_controlled ? source_values : delayed_values)
    end
    if length(source_values) == 2 && upper_phase == lower_phase && !half_controlled
        lower_phase = upper_phase == 1 ? 2 : 1
    end
    state = falses(length(study.topology.valve_positions))
    state[2 * upper_phase - 1] = true
    state[2 * lower_phase] = true
    return state
end

function _rectifier_open_circuit_dc_voltage(study, source_values)
    if study.specification.selection.family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    )
        return maximum(
            maximum(view(source_values, (3group - 2):(3group))) -
                minimum(view(source_values, (3group - 2):(3group)))
            for group in eachindex(study.group_phase_shift_rad)
        )
    end
    return maximum(source_values) - minimum(source_values)
end

function _seed_line_commutated_rectifier(study::SwitchingLineCommutatedRectifierStudy)
    ac_nodes, dc_positive, dc_negative = _line_commutated_rectifier_nodes(study)
    natural_path = _rectifier_conduction_path(study, study.start_time_s)
    valves = ntuple(length(study.topology.valve_positions)) do index
        _rectifier_valve(
            study,
            study.topology.valve_positions[index],
            natural_path[index] && study.initial_state.dc_current_a > 0.0,
        )
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, valves)
    bridge = Nonlinear.PowerSemiconductorBridgeTopology(study.topology, collect(valves))
    multipulse = study.specification.selection.family in (
        ConverterSystems.MultipulseDiodeBridge,
        ConverterSystems.MultipulseThyristorBridge,
    )
    source_neutral_nodes = if multipulse
        first_neutral = maximum((ac_nodes..., dc_positive, dc_negative)) + 1
        collect(first_neutral:(first_neutral + length(study.group_phase_shift_rad) - 1))
    else
        Int[0]
    end
    sources = if study.specification.selection.phase_count == 1
        [Branches.TwoTerminalTheveninSource(
            ac_nodes[1],
            ac_nodes[2],
            inv(study.source_resistance_ohm),
            time_s -> _rectifier_source_voltage(study, 1, time_s),
        )]
    else
        [Branches.TwoTerminalTheveninSource(
            ac_nodes[channel],
            multipulse ? source_neutral_nodes[(channel - 1) ÷ 3 + 1] : 0,
            inv(study.source_resistance_ohm),
            time_s -> _rectifier_source_voltage(study, channel, time_s),
        ) for channel in eachindex(ac_nodes)]
    end
    initial_current = study.initial_state.dc_current_a
    initial_sources = _rectifier_source_values(study, study.start_time_s)
    source_by_node = Dict(zip(ac_nodes, initial_sources))
    selected_positions = findall(natural_path)
    upper_position = only(index for index in selected_positions if isodd(index))
    lower_position = only(index for index in selected_positions if iseven(index))
    initial_applied_voltage =
        source_by_node[study.topology.valve_positions[upper_position].from_node] -
        source_by_node[study.topology.valve_positions[lower_position].to_node]
    initial_dc_voltage = initial_applied_voltage -
        2.0 * study.source_resistance_ohm * initial_current
    load = Branches.SeriesRLBranch(
        dc_positive,
        dc_negative,
        study.load_resistance_ohm,
        study.load_inductance_h,
        initial_current,
        initial_dc_voltage,
        initial_current,
    )
    node_count = maximum((ac_nodes..., source_neutral_nodes..., dc_positive, dc_negative))
    elements = Any[sources..., bridge, load]
    linear_network = Nodal.NodalSystem(node_count, elements)
    for (node, voltage) in zip(ac_nodes, initial_sources)
        linear_network.v[node] = voltage
    end
    linear_network.v[dc_positive] = max(initial_dc_voltage, 0.0)
    dc_negative == 0 || (linear_network.v[dc_negative] = 0.0)
    if _rectifier_is_detailed(study)
        for valve in valves
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _rectifier_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        valves;
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.source_peak_voltage_v,
            nominal_current_a=max(initial_current, study.source_peak_voltage_v /
                study.load_resistance_ohm, 1.0),
        ),
    ) : linear_network
    runtime = LineCommutatedRectifierRuntime(
        study,
        network,
        bridge,
        valves,
        load,
        ac_nodes,
        source_neutral_nodes,
        dc_positive,
        dc_negative,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
    if _rectifier_is_detailed(study) || !_rectifier_is_natural_diode(study)
        _synchronize_rectifier_path!(runtime, study.start_time_s; record=false)
    else
        _stabilize_rectifier_topology!(runtime, study.start_time_s,
            study.specification.timing.fixed_step_s)
    end
    return runtime
end

function LineCommutatedRectifierRecorder(phase_count, valve_count, sample_count)
    sample_count > 1 || throw(ArgumentError(
        "line-commutated rectifier recorder requires at least two samples",
    ))
    return LineCommutatedRectifierRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        falses(valve_count, sample_count),
        falses(valve_count, sample_count),
        falses(valve_count, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, valve_count, sample_count),
        zeros(Float64, valve_count, sample_count),
        zeros(Float64, valve_count, sample_count),
        0,
    )
end

function _rectifier_source_currents(runtime, time_s=runtime.time_s)
    study = runtime.study
    network = _rectifier_linear_network(runtime)
    values = _rectifier_source_values(study, time_s)
    if study.specification.selection.phase_count == 1
        terminal_voltage = network.v[runtime.ac_nodes[1]] -
            network.v[runtime.ac_nodes[2]]
        return [(values[1] - values[2] - terminal_voltage) /
            study.source_resistance_ohm]
    end
    return [(values[channel] - (
        network.v[runtime.ac_nodes[channel]] -
        (runtime.source_neutral_nodes == Int[0] ? 0.0 :
            network.v[runtime.source_neutral_nodes[(channel - 1) ÷ 3 + 1]])
    )) / study.source_resistance_ohm for channel in eachindex(runtime.ac_nodes)]
end

function _rectifier_stored_energy(runtime)
    energy = 0.5 * runtime.study.load_inductance_h * runtime.load_branch.i_last^2
    _rectifier_is_detailed(runtime.study) || return energy
    return energy + sum(runtime.valves; init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _rectifier_linear_companion(runtime)
    network = _rectifier_linear_network(runtime)
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

function _rectifier_device_companion_residual(runtime)
    _rectifier_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in runtime.valves
    )
end

function _rectifier_power(runtime, time_s=runtime.time_s)
    study = runtime.study
    source_voltage = study.specification.selection.phase_count == 1 ?
        [_rectifier_source_voltage(study, 1, time_s)] :
        _rectifier_source_values(study, time_s)
    source_current = _rectifier_source_currents(runtime, time_s)
    load_current = runtime.load_branch.i_last
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in runtime.valves
    )
    return (
        source_w=dot(source_voltage, source_current),
        source_loss_w=study.source_resistance_ohm * sum(abs2, source_current),
        load_loss_w=study.load_resistance_ohm * load_current^2,
        semiconductor_loss_w=semiconductor_loss,
    )
end

function _observe_rectifier_substep!(
    accounting::LineCommutatedRectifierStepEnergyAccounting,
    runtime,
    _system,
    time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _rectifier_power(runtime, time_s)
    final_linear = _rectifier_linear_companion(runtime)
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 : throw(ArgumentError(
            "rectifier energy accounting requires an accepted companion method",
        ))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    accounting.source_energy_j += step_s *
        blend(accounting.previous_power.source_w, final_power.source_w)
    accounting.source_dissipated_energy_j += step_s *
        blend(accounting.previous_power.source_loss_w, final_power.source_loss_w)
    accounting.load_dissipated_energy_j += step_s *
        blend(accounting.previous_power.load_loss_w, final_power.load_loss_w)
    accounting.semiconductor_dissipated_energy_j += step_s * blend(
        accounting.previous_power.semiconductor_loss_w,
        final_power.semiconductor_loss_w,
    )
    accounting.linear_companion_energy_residual_j += step_s *
        blend(accounting.previous_linear_companion.terminal_power_w,
            final_linear.terminal_power_w) -
        step_s * blend(accounting.previous_linear_companion.dissipated_power_w,
            final_linear.dissipated_power_w) -
        (final_linear.stored_energy_j -
            accounting.previous_linear_companion.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _stabilize_rectifier_topology!(runtime, time_s, step_s; maximum_iterations=32)
    _rectifier_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed rectifier execution must use the nonlinear device solver",
    ))
    network = _rectifier_linear_network(runtime)
    for iteration in 1:maximum_iterations
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        selected = nothing
        for (index, valve) in enumerate(runtime.valves)
            action = _three_phase_vsc_topology_action(valve, network.v)
            action === nothing && continue
            selected = (index, valve, action)
            break
        end
        selected === nothing && return iteration
        index, valve, action = selected
        _apply_three_phase_vsc_topology_action!(valve, action.transition, time_s)
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            runtime.study.topology.valve_positions[index].name,
            true;
            message="accepted natural diode commutation",
        ))
    end
    throw(ArgumentError("line-commutated rectifier topology failed to stabilize"))
end

function _synchronize_rectifier_path!(runtime, time_s; record=true)
    target = _rectifier_conduction_path(runtime.study, time_s)
    changed = false
    for (index, valve) in enumerate(runtime.valves)
        if valve.gate_driver !== nothing
            Nonlinear.request_power_semiconductor_gate!(valve, target[index], time_s)
        end
        previous = valve.closed
        if target[index] && !previous
            Nonlinear.apply_power_semiconductor_forward_turn_on!(valve, time_s)
        elseif !target[index] && previous
            Nonlinear.apply_power_semiconductor_forward_extinction!(valve, time_s)
        end
        if previous != valve.closed
            changed = true
            record && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                time_s,
                valve.closed ? :forward_turn_on : :forward_extinction,
                runtime.study.topology.valve_positions[index].name,
                true;
                message=_rectifier_is_natural_diode(runtime.study) ?
                    "accepted source-driven natural diode commutation" :
                    "accepted phase-controlled line commutation",
            ))
        end
    end
    return changed
end

function _record_rectifier_sample!(recorder, runtime, index, energy_residual_w)
    study = runtime.study
    network = _rectifier_linear_network(runtime)
    source_values = _rectifier_source_values(study, runtime.time_s)
    source_currents = _rectifier_source_currents(runtime)
    if study.specification.selection.phase_count == 1
        recorder.source_voltage_v[1, index] = source_values[1] - source_values[2]
        recorder.source_current_a[1, index] = source_currents[1]
    else
        recorder.source_voltage_v[:, index] .= source_values
        recorder.source_current_a[:, index] .= source_currents
    end
    recorder.time_s[index] = runtime.time_s
    recorder.dc_voltage_v[index] = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    recorder.dc_current_a[index] = runtime.load_branch.i_last
    for (valve_index, valve) in enumerate(runtime.valves)
        if valve.gate_driver !== nothing
            recorder.requested_firing_state[valve_index, index] =
                valve.gate_driver.commanded_on
            recorder.applied_firing_state[valve_index, index] =
                valve.gate_driver.applied_on
        end
        recorder.conducting_state[valve_index, index] = valve.closed
        if _rectifier_is_detailed(study)
            state = Nonlinear.power_semiconductor_extended_state(valve)
            recorder.junction_temperature_k[valve_index, index] =
                state.junction_temperature_k
            recorder.recovered_charge_c[valve_index, index] =
                state.stored_recovery_charge_c
            recorder.turn_off_tail_current_a[valve_index, index] =
                state.tail_current_a
        else
            recorder.junction_temperature_k[valve_index, index] = 0.0
            recorder.recovered_charge_c[valve_index, index] = 0.0
            recorder.turn_off_tail_current_a[valve_index, index] = 0.0
        end
    end
    power = _rectifier_power(runtime)
    recorder.stored_energy_j[index] = _rectifier_stored_energy(runtime)
    recorder.semiconductor_loss_w[index] = power.semiconductor_loss_w
    recorder.kcl_residual_a[index] = _rectifier_is_detailed(study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[index] = energy_residual_w
    recorder.write_index = index
    return recorder
end

function prepare_switching_line_commutated_rectifier(
    study::SwitchingLineCommutatedRectifierStudy,
)
    runtime = _seed_line_commutated_rectifier(study)
    sample_count = Int(round(
        (study.stop_time_s - study.start_time_s) /
            study.specification.timing.fixed_step_s,
    )) + 1
    integrator = LineCommutatedRectifierIntegrator(
        runtime,
        TimestepTransaction(runtime),
        LineCommutatedRectifierRecorder(
            study.specification.selection.phase_count == 1 ? 1 :
                length(runtime.ac_nodes),
            length(runtime.valves),
            sample_count,
        ),
        0,
        false,
        false,
        nothing,
    )
    _record_rectifier_sample!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_line_commutated_rectifier!(integrator)
    integrator.failed && throw(ArgumentError(
        "line-commutated rectifier integrator is terminally failed: $(integrator.last_failure)",
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
    previous_energy = _rectifier_stored_energy(runtime)
    previous_power = _rectifier_power(runtime)
    previous_linear = _rectifier_linear_companion(runtime)
    previous_device_residual = _rectifier_device_companion_residual(runtime)
    previous_semiconductor_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in runtime.valves
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        accounting = LineCommutatedRectifierStepEnergyAccounting(
            previous_power,
            previous_linear,
        )
        semiconductor_increment = 0.0
        if _rectifier_is_detailed(study)
            topology_changed = _synchronize_rectifier_path!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=topology_changed ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=topology_changed ? :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) -> _observe_rectifier_substep!(
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
                for valve in runtime.valves
            )
            semiconductor_increment =
                final_semiconductor_dissipation - previous_semiconductor_dissipation
        else
            event_count_before = length(runtime.events)
            controlled_path_changed = !_rectifier_is_natural_diode(study) &&
                _synchronize_rectifier_path!(runtime, endpoint)
            _rectifier_is_natural_diode(study) &&
                _stabilize_rectifier_topology!(runtime, endpoint, step)
            topology_changed = controlled_path_changed ||
                length(runtime.events) != event_count_before
            network = _rectifier_linear_network(runtime)
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
                _observe_rectifier_substep!(accounting, runtime, network,
                    endpoint - half_step, half_step,
                    Branches.BackwardEulerCompanion, 1)
                _rectifier_is_natural_diode(study) &&
                    _stabilize_rectifier_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(network, endpoint, half_step, backward_euler)
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_rectifier_substep!(accounting, runtime, network,
                    endpoint, half_step, Branches.BackwardEulerCompanion, 2)
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_rectifier_substep!(accounting, runtime, network,
                    endpoint, step, Branches.TrapezoidalCompanion, 1)
            end
            runtime.last_kcl_residual_a =
                maximum(abs, network.y * network.v - network.rhs; init=0.0)
            semiconductor_increment = accounting.semiconductor_dissipated_energy_j
        end
        dissipation_increment = accounting.source_dissipated_energy_j +
            accounting.load_dissipated_energy_j + semiconductor_increment
        runtime.source_energy_j += accounting.source_energy_j
        runtime.dissipated_energy_j += dissipation_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        stored_energy = _rectifier_stored_energy(runtime)
        device_residual_increment = _rectifier_device_companion_residual(runtime) -
            previous_device_residual
        energy_residual = (
            accounting.source_energy_j - dissipation_increment -
            (stored_energy - previous_energy) -
            accounting.linear_companion_energy_residual_j -
            device_residual_increment
        ) / step
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_rectifier_sample!(
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

function _line_commutated_rectifier_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _rectifier_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    source_currents = _rectifier_source_currents(runtime)
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        study.specification.timing.fixed_step_s,
    )
    stored_energy = _rectifier_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        bridge_state.deterministic_signature,
        repr(runtime.time_s),
        repr(runtime.load_branch.i_last),
        repr(dc_voltage),
        repr(stored_energy),
        string(integrator.accepted_step_index),
        string(length(runtime.events)),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        vcat(_rectifier_source_values(study, runtime.time_s), dc_voltage),
        vcat(source_currents, runtime.load_branch.i_last),
        BitVector(map(runtime.valves) do valve
            valve.gate_driver === nothing ? false : valve.gate_driver.commanded_on
        end),
        BitVector(map(runtime.valves) do valve
            valve.gate_driver === nothing ? false : valve.gate_driver.applied_on
        end),
        BitVector(getfield.(runtime.valves, :closed)),
        Float64[],
        [study.load_inductance_h * runtime.load_branch.i_last],
        Float64[],
        [2pi * study.source_frequency_hz * runtime.time_s],
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    stored_energy_change = stored_energy - recorder.stored_energy_j[1]
    integrated_residual = runtime.source_energy_j - runtime.dissipated_energy_j -
        stored_energy_change - runtime.linear_companion_energy_residual_j -
        _rectifier_device_companion_residual(runtime)
    scale = max(
        abs(runtime.source_energy_j),
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
    return ConverterSystems.SwitchingLineCommutatedRectifierTrace(
        recorder.time_s,
        recorder.source_voltage_v,
        recorder.source_current_a,
        recorder.dc_voltage_v,
        recorder.dc_current_a,
        recorder.requested_firing_state,
        recorder.applied_firing_state,
        recorder.conducting_state,
        recorder.stored_energy_j,
        recorder.semiconductor_loss_w,
        recorder.kcl_residual_a,
        recorder.energy_residual_w,
        recorder.junction_temperature_k,
        recorder.recovered_charge_c,
        recorder.turn_off_tail_current_a,
        result,
    )
end

function execute_switching_line_commutated_rectifier!(
    integrator::LineCommutatedRectifierIntegrator,
)
    while _advance_line_commutated_rectifier!(integrator)
    end
    return _line_commutated_rectifier_result(integrator)
end

using ..ConverterSystems

mutable struct SwitchingChopperPhysicalRuntime{S,N,T,B,C}
    study::S
    network::N
    bridge::T
    inductor::B
    capacitor::C
    input_positive_node::Int
    switching_node::Int
    output_node::Int
    time_s::Float64
    input_energy_j::Float64
    load_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

const SwitchingBuckPhysicalRuntime = SwitchingChopperPhysicalRuntime

mutable struct SwitchingChopperRecorder
    time_s::Vector{Float64}
    input_voltage_v::Vector{Float64}
    input_current_a::Vector{Float64}
    switch_node_voltage_v::Vector{Float64}
    output_voltage_v::Vector{Float64}
    inductor_current_a::Vector{Float64}
    load_current_a::Vector{Float64}
    requested_gate_state::BitVector
    applied_gate_state::BitVector
    controlled_valve_conducting_state::BitVector
    freewheel_diode_conducting_state::BitVector
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Vector{Float64}
    freewheel_junction_temperature_k::Vector{Float64}
    freewheel_recovered_charge_c::Vector{Float64}
    controlled_turn_off_tail_current_a::Vector{Float64}
    controlled_junction_stored_energy_j::Vector{Float64}
    freewheel_junction_stored_energy_j::Vector{Float64}
    write_index::Int
end


const SwitchingBuckRecorder = SwitchingChopperRecorder

mutable struct SwitchingChopperConverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::SwitchingChopperRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end


const SwitchingBuckConverterIntegrator = SwitchingChopperConverterIntegrator

mutable struct SwitchingChopperStepEnergyAccounting{P,L}
    previous_power::P
    previous_linear_companion::L
    input_energy_j::Float64
    load_energy_j::Float64
    series_dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
end


const SwitchingBuckStepEnergyAccounting = SwitchingChopperStepEnergyAccounting

SwitchingChopperStepEnergyAccounting(previous_power, previous_linear_companion) =
    SwitchingChopperStepEnergyAccounting(
        previous_power,
        previous_linear_companion,
        0.0,
        0.0,
        0.0,
        0.0,
    )

function _switching_chopper_topology_nodes(study::SwitchingChopperConverterStudy)
    topology_nodes = Dict(node.name => node.node for node in study.topology.nodes)
    input_port = only(port for port in study.specification.ports if port.identity === :input_dc)
    output_port = only(port for port in study.specification.ports if port.identity === :output_dc)
    length(input_port.ordered_nodes) == 2 && length(output_port.ordered_nodes) == 2 ||
        throw(ArgumentError("switching chopper input and output ports must each have two terminals"))
    family = study.specification.selection.family
    if family === ConverterSystems.BuckChopper
        all(haskey(topology_nodes, name) for name in (:dc_positive, :output, :dc_negative)) ||
            throw(ArgumentError("switching chopper topology is missing a canonical named terminal"))
        input_port.ordered_nodes == (topology_nodes[:dc_positive], topology_nodes[:dc_negative]) ||
            throw(ArgumentError("switching chopper input port must match the exact topology rails"))
        output_port.ordered_nodes[2] == topology_nodes[:dc_negative] || throw(ArgumentError(
            "switching chopper output negative must match the topology DC-negative rail",
        ))
        output_port.ordered_nodes[1] != topology_nodes[:output] || throw(ArgumentError(
            "switching chopper output port must follow rather than bypass its series inductor",
        ))
        return (
            input_positive=topology_nodes[:dc_positive],
            switching=topology_nodes[:output],
            negative=topology_nodes[:dc_negative],
            output=output_port.ordered_nodes[1],
        )
    elseif family === ConverterSystems.BoostChopper
        all(haskey(topology_nodes, name) for name in
            (:input_positive, :output_positive, :dc_negative)) || throw(ArgumentError(
                "switching boost topology is missing a canonical named terminal",
            ))
        input_port.ordered_nodes[2] == topology_nodes[:dc_negative] &&
            output_port.ordered_nodes ==
                (topology_nodes[:output_positive], topology_nodes[:dc_negative]) ||
            throw(ArgumentError("switching boost ports do not match the topology rails"))
        input_port.ordered_nodes[1] != topology_nodes[:input_positive] || throw(ArgumentError(
            "switching boost input port must precede rather than bypass its series inductor",
        ))
        return (
            input_positive=input_port.ordered_nodes[1],
            switching=topology_nodes[:input_positive],
            negative=topology_nodes[:dc_negative],
            output=topology_nodes[:output_positive],
        )
    elseif family === ConverterSystems.InvertingBuckBoostChopper
        all(haskey(topology_nodes, name) for name in
            (:input_positive, :switching, :output_negative, :reference)) ||
            throw(ArgumentError(
                "switching inverting buck-boost topology is missing a canonical named terminal",
            ))
        input_port.ordered_nodes ==
            (topology_nodes[:input_positive], topology_nodes[:reference]) ||
            throw(ArgumentError(
                "switching inverting buck-boost input port does not match its source rails",
            ))
        output_port.ordered_nodes ==
            (topology_nodes[:output_negative], topology_nodes[:reference]) ||
            throw(ArgumentError(
                "switching inverting buck-boost output port must retain its signed negative reference",
            ))
        return (
            input_positive=topology_nodes[:input_positive],
            switching=topology_nodes[:switching],
            negative=topology_nodes[:reference],
            output=topology_nodes[:output_negative],
        )
    end
    throw(ArgumentError("switching chopper topology nodes are unavailable for family $family"))
end

function _switching_chopper_gate_state(study, time_s)
    pwm = ConverterSystems.converter_pwm_gate_state(
        [study.specification.modulation.duty],
        time_s,
        study.specification.timing.carrier_frequency_hz,
        study.specification.modulation.kind,
    )
    return !ConverterSystems.converter_system_is_blocked(
        study.specification,
        time_s,
        1,
    ) && pwm.requested_valve_state[1]
end

function _detailed_chopper_switching_energy_table(parameters; controlled::Bool)
    current_axis = [0.0, parameters.switching_energy_current_domain_maximum_a]
    voltage_axis = [0.0, parameters.switching_energy_blocking_voltage_domain_maximum_v]
    temperature_axis = [parameters.ambient_temperature_k, 600.0]
    dimensions = (2, 2, 2)
    normalized_energy(rated_energy_j) = [
        rated_energy_j * (current_index - 1) * (voltage_index - 1) *
            (0.8 + 0.2 * (temperature_index - 1))
        for current_index in 1:dimensions[1],
            voltage_index in 1:dimensions[2],
            temperature_index in 1:dimensions[3]
    ]
    return Nonlinear.SwitchingEnergyTable(
        current_axis,
        voltage_axis,
        temperature_axis;
        turn_on_energy_j=normalized_energy(
            controlled ? parameters.turn_on_energy_at_domain_maximum_j : 0.0,
        ),
        turn_off_energy_j=normalized_energy(
            controlled ? parameters.turn_off_energy_at_domain_maximum_j : 0.0,
        ),
        reverse_recovery_energy_j=normalized_energy(
            controlled ? 0.0 : parameters.reverse_recovery_energy_at_domain_maximum_j,
        ),
        provenance=parameters.provenance,
    )
end

function _detailed_chopper_thermal_fidelity(parameters)
    return Nonlinear.CauerThermalFidelity(
        collect(parameters.thermal_capacitance_j_per_k),
        collect(parameters.thermal_resistance_k_per_w);
        ambient_temperature_k=parameters.ambient_temperature_k,
        provenance=parameters.provenance,
    )
end

function _detailed_chopper_controlled_fidelity(parameters)
    return Nonlinear.PowerSemiconductorExtendedFidelity(
        junction_charge=Nonlinear.NonlinearJunctionChargeFidelity(
            parameters.controlled_junction_capacitance_f,
            parameters.junction_potential_v,
            parameters.junction_grading_exponent;
            voltage_domain_v=(-parameters.junction_voltage_limit_v,
                parameters.junction_voltage_limit_v),
            provenance=parameters.provenance,
        ),
        turn_off_tail=Nonlinear.TurnOffTailFidelity(
            parameters.turn_off_tail_time_s;
            cutoff_current_a=parameters.turn_off_tail_cutoff_current_a,
            provenance=parameters.provenance,
        ),
        switching_energy=_detailed_chopper_switching_energy_table(
            parameters;
            controlled=true,
        ),
        thermal=_detailed_chopper_thermal_fidelity(parameters),
        provenance=parameters.provenance,
    )
end

function _detailed_chopper_freewheel_fidelity(parameters)
    return Nonlinear.PowerSemiconductorExtendedFidelity(
        recovered_charge=Nonlinear.RecoveredChargeFidelity(
            parameters.recovered_charge_lifetime_s;
            provenance=parameters.provenance,
        ),
        junction_charge=Nonlinear.NonlinearJunctionChargeFidelity(
            parameters.freewheel_junction_capacitance_f,
            parameters.junction_potential_v,
            parameters.junction_grading_exponent;
            voltage_domain_v=(-parameters.junction_voltage_limit_v,
                parameters.junction_voltage_limit_v),
            provenance=parameters.provenance,
        ),
        switching_energy=_detailed_chopper_switching_energy_table(
            parameters;
            controlled=false,
        ),
        thermal=_detailed_chopper_thermal_fidelity(parameters),
        provenance=parameters.provenance,
    )
end

_switching_chopper_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

function _switching_chopper_linear_network(runtime)
    return _switching_chopper_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network
end

function _seed_switching_chopper_network(study::SwitchingChopperConverterStudy)
    nodes = _switching_chopper_topology_nodes(study)
    initial = study.initial_state
    controlled_on = _switching_chopper_gate_state(study, study.start_time_s)
    detailed = _switching_chopper_is_detailed(study)
    semiconductor = study.detailed_semiconductor
    family = study.specification.selection.family
    controlled_position, freewheel_position = study.topology.valve_positions
    controlled_terminals = (controlled_position.from_node, controlled_position.to_node)
    freewheel_terminals = (freewheel_position.from_node, freewheel_position.to_node)
    controlled = Nonlinear.IGBTSwitch(
        controlled_terminals...;
        gate_driver=Nonlinear.PowerSemiconductorGateDriver(initially_on=controlled_on),
        on_conductance=detailed ? semiconductor.controlled_on_conductance_s : 1.0e6,
        off_conductance=detailed ? semiconductor.controlled_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? semiconductor.controlled_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            semiconductor.snubber_resistance_ohm,
            semiconductor.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ? _detailed_chopper_controlled_fidelity(semiconductor) : nothing,
        initially_closed=controlled_on,
    )
    freewheel = Nonlinear.DiodeValveSwitch(
        freewheel_terminals...;
        on_conductance=detailed ? semiconductor.freewheel_on_conductance_s : 1.0e6,
        off_conductance=detailed ? semiconductor.freewheel_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? semiconductor.freewheel_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            semiconductor.snubber_resistance_ohm,
            semiconductor.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ? _detailed_chopper_freewheel_fidelity(semiconductor) : nothing,
    )
    Nonlinear.power_semiconductor_event_localization!(controlled)
    Nonlinear.power_semiconductor_event_localization!(freewheel)
    bridge = Nonlinear.PowerSemiconductorBridgeTopology(
        study.topology,
        [controlled, freewheel],
    )
    source_carries_inductor_current = family === ConverterSystems.BoostChopper ||
        controlled_on
    initial_input_current = source_carries_inductor_current ?
        initial.inductor_current_a : 0.0
    initial_source_drop = study.source_resistance_ohm * initial_input_current
    initial_input_voltage = study.input_voltage_v - initial_source_drop
    initial_switch_voltage = if family === ConverterSystems.BuckChopper
        controlled_on ? initial_input_voltage : 0.0
    elseif family === ConverterSystems.BoostChopper
        controlled_on ? 0.0 : initial.output_voltage_v
    else
        controlled_on ? initial_input_voltage : initial.output_voltage_v
    end
    initial_inductor_voltage = if family === ConverterSystems.BuckChopper
        initial_switch_voltage - initial.output_voltage_v
    elseif family === ConverterSystems.BoostChopper
        initial_input_voltage - initial_switch_voltage
    else
        initial_switch_voltage
    end
    load_current = initial.output_voltage_v / study.load_resistance_ohm
    initial_capacitor_current = if family === ConverterSystems.BuckChopper
        initial.inductor_current_a - load_current
    elseif family === ConverterSystems.BoostChopper
        controlled_on ? -load_current : initial.inductor_current_a - load_current
    else
        controlled_on ? -load_current : -initial.inductor_current_a - load_current
    end
    source = Branches.TwoTerminalTheveninSource(
        nodes.input_positive,
        nodes.negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    inductor = Branches.SeriesRLBranch(
        family === ConverterSystems.BuckChopper ? nodes.switching :
            family === ConverterSystems.BoostChopper ? nodes.input_positive : nodes.switching,
        family === ConverterSystems.BuckChopper ? nodes.output :
            family === ConverterSystems.BoostChopper ? nodes.switching : nodes.negative,
        study.inductor_resistance_ohm,
        study.inductance_h,
        initial.inductor_current_a,
        initial_inductor_voltage,
        initial.inductor_current_a,
    )
    capacitor = Branches.CapacitorBranch(
        nodes.output,
        nodes.negative,
        study.capacitance_f,
        initial_capacitor_current,
        initial.output_voltage_v,
        initial_capacitor_current,
    )
    load = Branches.ConductanceBranch(
        nodes.output,
        nodes.negative,
        inv(study.load_resistance_ohm),
    )
    node_count = maximum((
        getfield.(study.topology.nodes, :node)...,
        nodes.input_positive,
        nodes.switching,
        nodes.output,
    ))
    linear_network = Nodal.NodalSystem(
        node_count,
        Any[source, bridge, inductor, capacitor, load],
    )
    nodes.negative == 0 || (linear_network.v[nodes.negative] = 0.0)
    linear_network.v[nodes.input_positive] = initial_input_voltage
    linear_network.v[nodes.switching] = initial_switch_voltage
    linear_network.v[nodes.output] = initial.output_voltage_v
    if detailed
        for valve in (controlled, freewheel)
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = detailed ? NonlinearNodalSystem(
        linear_network,
        (controlled, freewheel);
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=max(abs(initial.inductor_current_a), 1.0),
        ),
    ) : linear_network
    return SwitchingChopperPhysicalRuntime(
        study,
        network,
        bridge,
        inductor,
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

function SwitchingChopperRecorder(sample_count::Int)
    sample_count > 1 || throw(ArgumentError("switching chopper recorder requires at least two samples"))
    real_traces = ntuple(_ -> zeros(Float64, sample_count), 17)
    boolean_traces = ntuple(_ -> falses(sample_count), 4)
    return SwitchingChopperRecorder(
        real_traces[1:7]...,
        boolean_traces...,
        real_traces[8:17]...,
        0,
    )
end

function _switching_chopper_stored_energy(runtime)
    network = _switching_chopper_linear_network(runtime)
    output_voltage = network.v[runtime.output_node]
    stored_energy = 0.5 * runtime.study.inductance_h * runtime.inductor.i_last^2 +
        0.5 * runtime.study.capacitance_f * output_voltage^2
    _switching_chopper_is_detailed(runtime.study) || return stored_energy
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        runtime.study.specification.timing.fixed_step_s,
    )
    return stored_energy + bridge_state.stored_energy_j
end

function _switching_chopper_linear_companion_accounting(runtime)
    network = _switching_chopper_linear_network(runtime)
    inductor_voltage = Branches.branch_voltage(
        network.v,
        runtime.inductor.a,
        runtime.inductor.b,
    )
    capacitor_voltage = Branches.branch_voltage(
        network.v,
        runtime.capacitor.a,
        runtime.capacitor.b,
    )
    inductor_current = runtime.inductor.i_last
    capacitor_current = runtime.capacitor.i_last
    return (
        terminal_power_w=inductor_voltage * inductor_current +
            capacitor_voltage * capacitor_current,
        dissipated_power_w=runtime.study.inductor_resistance_ohm * inductor_current^2,
        stored_energy_j=0.5 * runtime.study.inductance_h * inductor_current^2 +
            0.5 * runtime.study.capacitance_f * capacitor_voltage^2,
    )
end

function _switching_chopper_bridge_companion_energy_residual(runtime)
    _switching_chopper_is_detailed(runtime.study) || return 0.0
    network = _switching_chopper_linear_network(runtime)
    return Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        runtime.study.specification.timing.fixed_step_s,
    ).companion_energy_residual_j
end

function _switching_chopper_power(runtime)
    study = runtime.study
    network = _switching_chopper_linear_network(runtime)
    input_node_voltage = network.v[runtime.input_positive_node]
    input_current = (study.input_voltage_v - input_node_voltage) /
        study.source_resistance_ohm
    inductor_current = runtime.inductor.i_last
    output_voltage = network.v[runtime.output_node]
    semiconductor_loss_w = sum(
        valve.last_semiconductor_loss_w for valve in runtime.bridge.valves
    )
    series_loss_w = study.source_resistance_ohm * input_current^2 +
        study.inductor_resistance_ohm * inductor_current^2
    return (
        input_w=study.input_voltage_v * input_current,
        load_w=output_voltage^2 / study.load_resistance_ohm,
        series_loss_w,
        semiconductor_loss_w,
        loss_w=series_loss_w + semiconductor_loss_w,
    )
end

function _observe_switching_chopper_substep!(
    accounting::SwitchingChopperStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _switching_chopper_power(runtime)
    final_linear = _switching_chopper_linear_companion_accounting(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    if companion_method === Branches.TrapezoidalCompanion
        input_energy_j = 0.5 * step_s * (previous_power.input_w + final_power.input_w)
        load_energy_j = 0.5 * step_s * (previous_power.load_w + final_power.load_w)
        series_dissipated_energy_j = 0.5 * step_s * (
            previous_power.series_loss_w + final_power.series_loss_w
        )
        linear_terminal_work_j = 0.5 * step_s * (
            previous_linear.terminal_power_w + final_linear.terminal_power_w
        )
        linear_dissipated_energy_j = 0.5 * step_s * (
            previous_linear.dissipated_power_w + final_linear.dissipated_power_w
        )
    elseif companion_method === Branches.BackwardEulerCompanion
        input_energy_j = step_s * final_power.input_w
        load_energy_j = step_s * final_power.load_w
        series_dissipated_energy_j = step_s * final_power.series_loss_w
        linear_terminal_work_j = step_s * final_linear.terminal_power_w
        linear_dissipated_energy_j = step_s * final_linear.dissipated_power_w
    else
        throw(ArgumentError(
            "switching chopper energy accounting requires an accepted trapezoidal or backward-Euler companion substep",
        ))
    end
    accounting.input_energy_j += input_energy_j
    accounting.load_energy_j += load_energy_j
    accounting.series_dissipated_energy_j += series_dissipated_energy_j
    accounting.linear_companion_energy_residual_j += linear_terminal_work_j -
        linear_dissipated_energy_j -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _stabilize_switching_chopper_topology!(runtime, time_s, step_s; maximum_iterations=16)
    _switching_chopper_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed chopper execution must use the nonlinear D200 device solver",
    ))
    network = _switching_chopper_linear_network(runtime)
    Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(runtime.bridge, time_s)
    valves = Nonlinear.power_semiconductor_bridge_topology_valves(runtime.bridge)
    for iteration in 1:maximum_iterations
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        selected = nothing
        for (position_index, valve) in enumerate(valves)
            action = _three_phase_vsc_topology_action(valve, network.v)
            action === nothing && continue
            selected = (position_index, valve, action)
            break
        end
        selected === nothing && return iteration
        position_index, valve, action = selected
        _apply_three_phase_vsc_topology_action!(valve, action.transition, time_s)
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            action.transition,
            runtime.study.topology.valve_positions[position_index].name,
            true;
            message="accepted coupled switching-state commutation",
        ))
    end
    throw(ArgumentError("switching chopper topology failed to stabilize"))
end

function _apply_switching_chopper_event_commands!(runtime, time_s)
    commands = ConverterSystems.converter_system_commands_at_boundary(
        runtime.study.specification,
        time_s,
    )
    gate_changed = false
    for command in commands
        controlled_driver = runtime.bridge.valves[1].gate_driver
        previous_gate_command = controlled_driver.commanded_on
        kind = command.kind === ConverterSystems.ConverterBlockEvent ?
            :converter_block : :converter_restart
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            kind,
            command.id,
            true;
            message="accepted converter command before the same-boundary PWM transition",
        ))
        changed = if command.kind === ConverterSystems.ConverterBlockEvent
            Nonlinear.block_power_semiconductor_topology!(runtime.bridge, time_s)
        else
            Nonlinear.restart_power_semiconductor_topology!(runtime.bridge, time_s)
        end
        changed || throw(ArgumentError(
            "converter-system command did not change its declared block state",
        ))
        current_gate_command = controlled_driver.commanded_on
        if current_gate_command != previous_gate_command
            gate_changed = true
            push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                time_s,
                current_gate_command ? :gate_turn_on : :gate_turn_off,
                :controlled,
                true;
                message="accepted block-command gate transition on the fixed-step calendar",
            ))
        end
    end
    return !isempty(commands), gate_changed
end

function _request_switching_chopper_gate!(runtime, time_s)
    requested_on = _switching_chopper_gate_state(runtime.study, time_s)
    previous = runtime.bridge.valves[1].gate_driver.commanded_on
    Nonlinear.request_power_semiconductor_topology_gates!(
        runtime.bridge,
        Bool[requested_on, false],
        time_s,
    )
    if requested_on != previous
        push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
            time_s,
            requested_on ? :gate_turn_on : :gate_turn_off,
            :controlled,
            true;
            message="accepted carrier-state change on the fixed-step calendar",
        ))
    end
    return requested_on, requested_on != previous
end

function _synchronize_detailed_chopper_commutation!(runtime, time_s)
    _switching_chopper_is_detailed(runtime.study) || return false
    controlled, freewheel = runtime.bridge.valves
    requested_on = controlled.gate_driver.applied_on
    changed = false
    if requested_on
        if freewheel.closed
            Nonlinear.apply_power_semiconductor_forward_extinction!(freewheel, time_s)
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
        should_freewheel = runtime.inductor.i_last > freewheel.holding_current
        if should_freewheel && !freewheel.closed
            Nonlinear.apply_power_semiconductor_forward_turn_on!(freewheel, time_s)
            changed = true
        elseif !should_freewheel && freewheel.closed
            Nonlinear.apply_power_semiconductor_forward_extinction!(freewheel, time_s)
            changed = true
        end
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        requested_on ? :controlled_conduction : :freewheel_commutation,
        requested_on ? :controlled : :freewheel,
        true;
        message="accepted detailed chopper commutation from gate state and inductor-current continuity",
    ))
    return changed
end

function _record_switching_chopper_sample!(recorder, runtime, sample_index, energy_residual_w)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    network = _switching_chopper_linear_network(runtime)
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        step,
    )
    input_node_voltage = network.v[runtime.input_positive_node]
    input_current = (study.input_voltage_v - input_node_voltage) / study.source_resistance_ohm
    output_voltage = network.v[runtime.output_node]
    load_current = output_voltage / study.load_resistance_ohm
    kcl_residual = if _switching_chopper_is_detailed(study)
        runtime.last_kcl_residual_a
    else
        kcl = network.y * network.v - network.rhs
        maximum(abs, kcl; init=0.0)
    end
    recorder.time_s[sample_index] = runtime.time_s
    recorder.input_voltage_v[sample_index] = study.input_voltage_v
    recorder.input_current_a[sample_index] = input_current
    recorder.switch_node_voltage_v[sample_index] = network.v[runtime.switching_node]
    recorder.output_voltage_v[sample_index] = output_voltage
    recorder.inductor_current_a[sample_index] = runtime.inductor.i_last
    recorder.load_current_a[sample_index] = load_current
    recorder.requested_gate_state[sample_index] = bridge_state.requested_gate_state[1]
    recorder.applied_gate_state[sample_index] = bridge_state.applied_gate_state[1]
    recorder.controlled_valve_conducting_state[sample_index] = bridge_state.conducting_state[1]
    recorder.freewheel_diode_conducting_state[sample_index] = bridge_state.conducting_state[2]
    recorder.stored_energy_j[sample_index] = _switching_chopper_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample_index] = bridge_state.semiconductor_loss_w
    recorder.kcl_residual_a[sample_index] = kcl_residual
    recorder.energy_residual_w[sample_index] = energy_residual_w
    if _switching_chopper_is_detailed(study)
        controlled_state = Nonlinear.power_semiconductor_extended_state(runtime.bridge.valves[1])
        freewheel_state = Nonlinear.power_semiconductor_extended_state(runtime.bridge.valves[2])
        recorder.controlled_junction_temperature_k[sample_index] =
            controlled_state.junction_temperature_k
        recorder.freewheel_junction_temperature_k[sample_index] =
            freewheel_state.junction_temperature_k
        recorder.freewheel_recovered_charge_c[sample_index] =
            freewheel_state.stored_recovery_charge_c
        recorder.controlled_turn_off_tail_current_a[sample_index] =
            controlled_state.tail_current_a
        recorder.controlled_junction_stored_energy_j[sample_index] =
            controlled_state.junction_stored_energy_j
        recorder.freewheel_junction_stored_energy_j[sample_index] =
            freewheel_state.junction_stored_energy_j
    else
        recorder.controlled_junction_temperature_k[sample_index] = 0.0
        recorder.freewheel_junction_temperature_k[sample_index] = 0.0
        recorder.freewheel_recovered_charge_c[sample_index] = 0.0
        recorder.controlled_turn_off_tail_current_a[sample_index] = 0.0
        recorder.controlled_junction_stored_energy_j[sample_index] = 0.0
        recorder.freewheel_junction_stored_energy_j[sample_index] = 0.0
    end
    recorder.write_index = sample_index
    return recorder
end

function prepare_switching_chopper_converter(study::SwitchingChopperConverterStudy)
    runtime = _seed_switching_chopper_network(study)
    sample_count = Int(round(
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    )) + 1
    integrator = SwitchingChopperConverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        SwitchingChopperRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_switching_chopper_sample!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

prepare_switching_buck_converter(study::SwitchingChopperConverterStudy) =
    prepare_switching_chopper_converter(study)

function _advance_switching_chopper_converter!(integrator)
    integrator.failed && throw(ArgumentError(
        "switching chopper integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _switching_chopper_power(runtime)
    previous_semiconductor_dissipated_energy_j = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in runtime.bridge.valves
    )
    previous_linear_companion = _switching_chopper_linear_companion_accounting(runtime)
    previous_linear_companion_energy_residual_j =
        runtime.linear_companion_energy_residual_j
    previous_bridge_companion_energy_residual_j =
        _switching_chopper_bridge_companion_energy_residual(runtime)
    detailed_input_energy_increment_j = 0.0
    detailed_load_energy_increment_j = 0.0
    detailed_dissipated_energy_increment_j = 0.0
    detailed_linear_companion_energy_residual_increment_j = 0.0
    begin_timestep_transaction!(integrator.transaction)
    try
        _, command_gate_changed = _apply_switching_chopper_event_commands!(runtime, endpoint)
        _, pwm_gate_changed = _request_switching_chopper_gate!(runtime, endpoint)
        gate_changed = command_gate_changed || pwm_gate_changed
        detailed = _switching_chopper_is_detailed(study)
        if detailed
            commutation_changed = _synchronize_detailed_chopper_commutation!(runtime, endpoint)
            Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                runtime.bridge,
                endpoint,
            )
            energy_accounting = SwitchingChopperStepEnergyAccounting(
                previous_power,
                previous_linear_companion,
            )
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || commutation_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || commutation_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) ->
                    _observe_switching_chopper_substep!(
                        energy_accounting,
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
            final_semiconductor_dissipated_energy_j = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in runtime.bridge.valves
            )
            detailed_input_energy_increment_j = energy_accounting.input_energy_j
            detailed_load_energy_increment_j = energy_accounting.load_energy_j
            detailed_dissipated_energy_increment_j =
                energy_accounting.series_dissipated_energy_j +
                (final_semiconductor_dissipated_energy_j -
                 previous_semiconductor_dissipated_energy_j)
            detailed_linear_companion_energy_residual_increment_j =
                energy_accounting.linear_companion_energy_residual_j
            runtime.input_energy_j += detailed_input_energy_increment_j
            runtime.load_energy_j += detailed_load_energy_increment_j
            runtime.dissipated_energy_j += detailed_dissipated_energy_increment_j
            runtime.linear_companion_energy_residual_j +=
                detailed_linear_companion_energy_residual_increment_j
        else
            energy_accounting = SwitchingChopperStepEnergyAccounting(
                previous_power,
                previous_linear_companion,
            )
            event_count_before = length(runtime.events)
            _stabilize_switching_chopper_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count_before
            if topology_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                network = _switching_chopper_linear_network(runtime)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint - half_step,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_switching_chopper_substep!(
                    energy_accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_switching_chopper_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_switching_chopper_substep!(
                    energy_accounting,
                    runtime,
                    network,
                    endpoint,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    2,
                )
            else
                network = _switching_chopper_linear_network(runtime)
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_switching_chopper_substep!(
                    energy_accounting,
                    runtime,
                    network,
                    endpoint,
                    step,
                    Branches.TrapezoidalCompanion,
                    1,
                )
            end
            runtime.input_energy_j += energy_accounting.input_energy_j
            runtime.load_energy_j += energy_accounting.load_energy_j
            runtime.dissipated_energy_j += energy_accounting.series_dissipated_energy_j
            runtime.linear_companion_energy_residual_j +=
                energy_accounting.linear_companion_energy_residual_j
        end
        runtime.time_s = endpoint
        network = _switching_chopper_linear_network(runtime)
        current = runtime.inductor.i_last
        output_voltage = network.v[runtime.output_node]
        input_current = (study.input_voltage_v -
            network.v[runtime.input_positive_node]) / study.source_resistance_ohm
        load_power = output_voltage^2 / study.load_resistance_ohm
        series_loss = study.source_resistance_ohm * input_current^2 +
            study.inductor_resistance_ohm * current^2
        semiconductor_loss = sum(
            valve.last_semiconductor_loss_w for valve in runtime.bridge.valves
        )
        stored_energy = _switching_chopper_stored_energy(runtime)
        bridge_companion_energy_residual_increment_j =
            _switching_chopper_bridge_companion_energy_residual(runtime) -
            previous_bridge_companion_energy_residual_j
        linear_companion_energy_residual_increment_j =
            runtime.linear_companion_energy_residual_j -
                previous_linear_companion_energy_residual_j
        energy_residual = if detailed
            (
                detailed_input_energy_increment_j -
                detailed_load_energy_increment_j -
                detailed_dissipated_energy_increment_j -
                (stored_energy - previous_energy) -
                bridge_companion_energy_residual_increment_j -
                detailed_linear_companion_energy_residual_increment_j
            ) / step
        else
            study.input_voltage_v * input_current - load_power -
                series_loss - semiconductor_loss -
                (stored_energy - previous_energy) / step -
                (bridge_companion_energy_residual_increment_j +
                 linear_companion_energy_residual_increment_j) / step
        end
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_switching_chopper_sample!(
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

function _switching_chopper_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _switching_chopper_linear_network(runtime)
    bridge_state = Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridge,
        network.v,
        study.specification.timing.fixed_step_s,
    )
    output_voltage = network.v[runtime.output_node]
    input_current = (study.input_voltage_v -
        network.v[runtime.input_positive_node]) / study.source_resistance_ohm
    stored_energy = _switching_chopper_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        bridge_state.deterministic_signature,
        repr(runtime.time_s),
        repr(runtime.inductor.i_last),
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
        bridge_state.requested_gate_state,
        bridge_state.applied_gate_state,
        bridge_state.conducting_state,
        [study.capacitance_f * output_voltage],
        [study.inductance_h * runtime.inductor.i_last],
        Float64[],
        [study.specification.modulation.duty],
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    stored_energy_change_j = stored_energy - recorder.stored_energy_j[1]
    integrated_energy_residual_j = runtime.input_energy_j - runtime.load_energy_j -
        runtime.dissipated_energy_j - stored_energy_change_j -
        _switching_chopper_bridge_companion_energy_residual(runtime) -
        runtime.linear_companion_energy_residual_j
    energy_scale_j = max(
        abs(runtime.input_energy_j),
        abs(runtime.load_energy_j),
        abs(runtime.dissipated_energy_j),
        abs(stored_energy_change_j),
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
        relative_energy_residual=abs(integrated_energy_residual_j) / energy_scale_j,
    )
    return ConverterSystems.SwitchingChopperConverterTrace(
        study.specification.selection.fidelity,
        recorder.time_s,
        recorder.input_voltage_v,
        recorder.input_current_a,
        recorder.switch_node_voltage_v,
        recorder.output_voltage_v,
        recorder.inductor_current_a,
        recorder.load_current_a,
        recorder.requested_gate_state,
        recorder.applied_gate_state,
        recorder.controlled_valve_conducting_state,
        recorder.freewheel_diode_conducting_state,
        recorder.stored_energy_j,
        recorder.semiconductor_loss_w,
        recorder.kcl_residual_a,
        recorder.energy_residual_w,
        recorder.controlled_junction_temperature_k,
        recorder.freewheel_junction_temperature_k,
        recorder.freewheel_recovered_charge_c,
        recorder.controlled_turn_off_tail_current_a,
        recorder.controlled_junction_stored_energy_j,
        recorder.freewheel_junction_stored_energy_j,
        result,
    )
end

function execute_switching_chopper_converter!(integrator::SwitchingChopperConverterIntegrator)
    while _advance_switching_chopper_converter!(integrator)
    end
    return _switching_chopper_result(integrator)
end


execute_switching_buck_converter!(integrator::SwitchingChopperConverterIntegrator) =
    execute_switching_chopper_converter!(integrator)

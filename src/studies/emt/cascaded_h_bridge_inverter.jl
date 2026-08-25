mutable struct CascadedHBridgeInverterRuntime{S,N,B,C,A,P,L}
    study::S
    network::N
    bridges::B
    controlled_valves::C
    antiparallel_diodes::A
    cell_capacitors::P
    load_branches::L
    phase_nodes::NTuple{3,Int}
    phase_return_node::Int
    load_neutral_node::Int
    cell_nodes::Vector{NamedTuple{(:dc_positive, :dc_negative, :left, :right),NTuple{4,Int}}}
    time_s::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct CascadedHBridgeInverterRecorder
    time_s::Vector{Float64}
    cell_dc_voltage_v::Array{Float64,3}
    cell_dc_current_a::Array{Float64,3}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    requested_level::Matrix{Int8}
    requested_cell_state::Array{Int8,3}
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

mutable struct CascadedHBridgeInverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::CascadedHBridgeInverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_cascaded_h_bridge_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_cascaded_h_bridge_cell_count(study) = study.specification.selection.cell_count

function _cascaded_h_bridge_linear_network(runtime)
    return _cascaded_h_bridge_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network
end

function _cascaded_h_bridge_nodes(study)
    cell_count = _cascaded_h_bridge_cell_count(study)
    phase_nodes = ntuple(phase -> begin
        nodes = Dict(node.name => node.node for node in study.topologies[phase].nodes)
        nodes[:phase_positive]
    end, 3)
    phase_returns = ntuple(phase -> begin
        nodes = Dict(node.name => node.node for node in study.topologies[phase].nodes)
        nodes[:phase_negative]
    end, 3)
    length(unique(phase_nodes)) == 3 || throw(ArgumentError(
        "cascaded H-bridge phase outputs must be distinct",
    ))
    length(unique(phase_returns)) == 1 || throw(ArgumentError(
        "three-phase cascaded H-bridge phases must share one exact series return",
    ))
    output_port = only(port for port in study.specification.ports if
        port.identity === :output_ac)
    output_port.ordered_nodes == phase_nodes || throw(ArgumentError(
        "cascaded H-bridge AC port does not match its phase outputs",
    ))
    cell_nodes = NamedTuple{(:dc_positive, :dc_negative, :left, :right),NTuple{4,Int}}[]
    expected_dc_ports = Tuple{Int,Int}[]
    maximum_node = maximum((phase_nodes..., phase_returns...))
    for phase in 1:3
        nodes = Dict(node.name => node.node for node in study.topologies[phase].nodes)
        for cell in 1:cell_count
            dc_positive = nodes[Symbol(:cell_, cell, :_dc_positive)]
            dc_negative = nodes[Symbol(:cell_, cell, :_dc_negative)]
            left = cell == 1 ? phase_nodes[phase] :
                nodes[Symbol(:series_junction_, cell - 1)]
            right = cell == cell_count ? phase_returns[phase] :
                nodes[Symbol(:series_junction_, cell)]
            push!(cell_nodes, (; dc_positive, dc_negative, left, right))
            push!(expected_dc_ports, (dc_positive, dc_negative))
            maximum_node = max(maximum_node, dc_positive, dc_negative, left, right)
        end
    end
    actual_dc_ports = [port.ordered_nodes for port in study.specification.ports if
        port.kind === ConverterSystems.IsolatedDirectCurrentPort]
    sort(actual_dc_ports) == sort(expected_dc_ports) || throw(ArgumentError(
        "cascaded H-bridge isolated DC ports do not match every cell rail",
    ))
    load_neutral = maximum_node + 1
    return (
        phase=phase_nodes,
        phase_return=only(unique(phase_returns)),
        neutral=load_neutral,
        cell=cell_nodes,
        node_count=load_neutral,
    )
end

function _cascaded_h_bridge_reference(study, time_s)
    modulation = study.specification.modulation
    frequency = study.specification.rated_bases.frequency_hz
    return [modulation.modulation_index * sin(
        2.0 * pi * frequency * time_s + modulation.phase_shift_rad -
            (phase - 1) * 2.0 * pi / 3.0,
    ) for phase in 1:3]
end

function _cascaded_h_bridge_modulation_state(
    study,
    time_s,
    cell_voltage,
    phase_current,
)
    phase_mean = sum(cell_voltage; dims=2) ./ size(cell_voltage, 2)
    voltage_error = cell_voltage .- phase_mean
    modulation = study.specification.modulation
    return ConverterSystems.converter_cascaded_h_bridge_gate_state(
        _cascaded_h_bridge_reference(study, time_s),
        voltage_error,
        phase_current,
        time_s,
        study.specification.timing.carrier_frequency_hz,
        modulation.kind;
        modulation_index=modulation.modulation_index,
        selective_harmonic_angles_rad=modulation.selective_harmonic_angles_rad,
    )
end

function _cascaded_h_bridge_path_targets(gates, current, tolerance=0.0)
    controlled = falses(4)
    diodes = falses(4)
    active = findall(gates)
    if length(active) == 2
        state = Tuple(active)
        if state == (1, 4) || state == (2, 3)
            target = current >= -tolerance ? controlled : diodes
            target[collect(state)] .= true
        elseif state == (1, 3)
            if current >= 0.0
                diodes[1] = true
                controlled[3] = true
            else
                controlled[1] = true
                diodes[3] = true
            end
        elseif state == (2, 4)
            if current >= 0.0
                controlled[2] = true
                diodes[4] = true
            else
                diodes[2] = true
                controlled[4] = true
            end
        end
    elseif current > tolerance
        diodes[[2, 3]] .= true
    elseif current < -tolerance
        diodes[[1, 4]] .= true
    end
    return controlled, diodes
end

function _cascaded_h_bridge_initial_node_voltage!(
    voltage,
    study,
    nodes,
    modulation,
)
    cell_count = _cascaded_h_bridge_cell_count(study)
    phase_output = zeros(3)
    for phase in 1:3
        right_voltage = nodes.phase_return == 0 ? 0.0 : voltage[nodes.phase_return]
        for cell in cell_count:-1:1
            index = (phase - 1) * cell_count + cell
            cell_nodes = nodes.cell[index]
            cell_voltage = study.initial_state.cell_dc_voltage_v[phase, cell]
            state = modulation.requested_cell_state[phase, cell]
            left_voltage = right_voltage + state * cell_voltage
            gates = modulation.requested_valve_state[
                (4 * (index - 1) + 1):(4 * index),
            ]
            if gates == BitVector((true, false, false, true))
                voltage[cell_nodes.dc_positive] = left_voltage
                voltage[cell_nodes.dc_negative] = right_voltage
            elseif gates == BitVector((false, true, true, false))
                voltage[cell_nodes.dc_negative] = left_voltage
                voltage[cell_nodes.dc_positive] = right_voltage
            elseif gates == BitVector((true, false, true, false))
                voltage[cell_nodes.dc_positive] = left_voltage
                voltage[cell_nodes.dc_negative] = left_voltage - cell_voltage
            else
                voltage[cell_nodes.dc_negative] = left_voltage
                voltage[cell_nodes.dc_positive] = left_voltage + cell_voltage
            end
            cell_nodes.left == 0 || (voltage[cell_nodes.left] = left_voltage)
            cell_nodes.right == 0 || (voltage[cell_nodes.right] = right_voltage)
            right_voltage = left_voltage
        end
        phase_output[phase] = right_voltage
        voltage[nodes.phase[phase]] = right_voltage
    end
    return phase_output
end

function _seed_cascaded_h_bridge_inverter(study)
    nodes = _cascaded_h_bridge_nodes(study)
    cell_count = _cascaded_h_bridge_cell_count(study)
    initial_current = collect(study.initial_state.phase_current_a)
    initial_voltage = study.initial_state.cell_dc_voltage_v
    modulation = _cascaded_h_bridge_modulation_state(
        study,
        study.start_time_s,
        initial_voltage,
        initial_current,
    )
    valve_count = 12 * cell_count
    controlled_paths = falses(valve_count)
    diode_paths = falses(valve_count)
    for phase in 1:3, cell in 1:cell_count
        cell_index = (phase - 1) * cell_count + cell
        first = 4 * (cell_index - 1) + 1
        controlled, diodes = _cascaded_h_bridge_path_targets(
            modulation.requested_valve_state[first:(first + 3)],
            initial_current[phase],
        )
        controlled_paths[first:(first + 3)] .= controlled
        diode_paths[first:(first + 3)] .= diodes
    end
    controlled = [begin
        phase = div(index - 1, 4 * cell_count) + 1
        position_index = mod(index - 1, 4 * cell_count) + 1
        position = study.topologies[phase].valve_positions[position_index]
        _three_level_split_link_controlled_valve(
            study,
            position,
            modulation.requested_valve_state[index],
            controlled_paths[index],
        )
    end for index in 1:valve_count]
    antiparallel = [begin
        phase = div(index - 1, 4 * cell_count) + 1
        position_index = mod(index - 1, 4 * cell_count) + 1
        position = study.topologies[phase].valve_positions[position_index]
        _three_level_split_link_diode(
            study,
            position.to_node,
            position.from_node,
            diode_paths[index],
        )
    end for index in 1:valve_count]
    foreach(Nonlinear.power_semiconductor_event_localization!, controlled)
    foreach(Nonlinear.power_semiconductor_event_localization!, antiparallel)
    capacitors = [begin
        cell_nodes = nodes.cell[index]
        phase = div(index - 1, cell_count) + 1
        cell = mod(index - 1, cell_count) + 1
        Branches.CapacitorBranch(
            cell_nodes.dc_positive,
            cell_nodes.dc_negative,
            study.cell_dc_capacitance_f[phase, cell],
            0.0,
            initial_voltage[phase, cell],
            0.0,
        )
    end for index in 1:(3 * cell_count)]
    bridges = ntuple(3) do phase
        first = (phase - 1) * 4 * cell_count + 1
        last = phase * 4 * cell_count
        Nonlinear.PowerSemiconductorBridgeTopology(
            study.topologies[phase],
            Nonlinear.PowerSemiconductorSwitch[controlled[first:last]...],
        )
    end
    initial_node_voltage = zeros(nodes.node_count)
    phase_output = _cascaded_h_bridge_initial_node_voltage!(
        initial_node_voltage,
        study,
        nodes,
        modulation,
    )
    neutral_voltage = sum(phase_output) / 3.0
    initial_node_voltage[nodes.neutral] = neutral_voltage
    loads = ntuple(3) do phase
        Branches.SeriesRLBranch(
            nodes.phase[phase],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            initial_current[phase],
            phase_output[phase] - neutral_voltage,
            initial_current[phase],
        )
    end
    linear_elements = _cascaded_h_bridge_is_detailed(study) ?
        Any[bridges..., capacitors..., loads...] :
        Any[bridges..., antiparallel..., capacitors..., loads...]
    linear_network = Nodal.NodalSystem(nodes.node_count, linear_elements)
    linear_network.v .= initial_node_voltage
    devices = (controlled..., antiparallel...)
    if _cascaded_h_bridge_is_detailed(study)
        for valve in devices
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = if _cascaded_h_bridge_is_detailed(study)
        NonlinearNodalSystem(
            linear_network,
            devices;
            scales=NonlinearNetworkScales(
                nodes.node_count,
                0;
                nominal_voltage_v=study.specification.rated_bases.voltage_v,
                nominal_current_a=study.specification.rated_bases.current_a,
            ),
            options=NonlinearSolveOptions(
                maximum_iterations=60,
                current_absolute_tolerance_a=1.0e-9,
                current_relative_tolerance=2.5e-9,
                voltage_absolute_tolerance_v=1.0e-9,
                scaled_step_tolerance=1.0e-12,
                maximum_condition_estimate=2.0e13,
            ),
        )
    else
        linear_network
    end
    return CascadedHBridgeInverterRuntime(
        study,
        network,
        bridges,
        controlled,
        antiparallel,
        capacitors,
        loads,
        nodes.phase,
        nodes.phase_return,
        nodes.neutral,
        nodes.cell,
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function CascadedHBridgeInverterRecorder(cell_count, sample_count)
    valve_count = 12 * cell_count
    return CascadedHBridgeInverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, 3, cell_count, sample_count),
        zeros(Float64, 3, cell_count, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Int8, 3, sample_count),
        zeros(Int8, 3, cell_count, sample_count),
        falses(valve_count, sample_count),
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
        zeros(Float64, valve_count, sample_count),
        0,
    )
end

_cascaded_h_bridge_devices(runtime) =
    (runtime.controlled_valves..., runtime.antiparallel_diodes...)

function _cascaded_h_bridge_stored_energy(runtime)
    network = _cascaded_h_bridge_linear_network(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    energy += sum(runtime.cell_capacitors; init=0.0) do capacitor
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        0.5 * capacitor.c * voltage^2
    end
    _cascaded_h_bridge_is_detailed(runtime.study) || return energy
    return energy + sum(_cascaded_h_bridge_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _cascaded_h_bridge_power(runtime)
    load_loss = runtime.study.load_resistance_ohm *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _cascaded_h_bridge_devices(runtime)
    )
    return (; load_w=load_loss, semiconductor_loss_w=semiconductor_loss)
end

function _cascaded_h_bridge_linear_companion(runtime)
    network = _cascaded_h_bridge_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for branch in runtime.load_branches
        voltage = Branches.branch_voltage(network.v, branch.a, branch.b)
        terminal_power += voltage * branch.i_last
        dissipated_power += runtime.study.load_resistance_ohm * branch.i_last^2
        stored_energy += 0.5 * runtime.study.load_inductance_h * branch.i_last^2
    end
    for capacitor in runtime.cell_capacitors
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        terminal_power += voltage * capacitor.i_last
        stored_energy += 0.5 * capacitor.c * voltage^2
    end
    return (; terminal_power_w=terminal_power,
        dissipated_power_w=dissipated_power, stored_energy_j=stored_energy)
end

function _cascaded_h_bridge_device_residual(runtime)
    _cascaded_h_bridge_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _cascaded_h_bridge_devices(runtime)
    )
end

function _observe_cascaded_h_bridge_substep!(
    accounting,
    runtime,
    _system,
    _time,
    step,
    method,
    _index,
)
    final_power = _cascaded_h_bridge_power(runtime)
    final_linear = _cascaded_h_bridge_linear_companion(runtime)
    weight = method === Branches.TrapezoidalCompanion ? 0.5 :
        method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("cascaded H-bridge received an unknown companion method"))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    accounting.load_energy_j += step * blend(previous_power.load_w, final_power.load_w)
    accounting.semiconductor_dissipated_energy_j += step * blend(
        previous_power.semiconductor_loss_w,
        final_power.semiconductor_loss_w,
    )
    accounting.linear_companion_energy_residual_j +=
        step * blend(previous_linear.terminal_power_w, final_linear.terminal_power_w) -
        step * blend(previous_linear.dissipated_power_w, final_linear.dissipated_power_w) -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _current_cascaded_h_bridge_modulation(runtime, time_s)
    network = _cascaded_h_bridge_linear_network(runtime)
    cell_count = _cascaded_h_bridge_cell_count(runtime.study)
    cell_voltage = Matrix{Float64}(undef, 3, cell_count)
    for phase in 1:3, cell in 1:cell_count
        index = (phase - 1) * cell_count + cell
        capacitor = runtime.cell_capacitors[index]
        cell_voltage[phase, cell] = Branches.branch_voltage(
            network.v,
            capacitor.a,
            capacitor.b,
        )
    end
    phase_current = [branch.i_last for branch in runtime.load_branches]
    return _cascaded_h_bridge_modulation_state(
        runtime.study,
        time_s,
        cell_voltage,
        phase_current,
    )
end

function _request_cascaded_h_bridge_gates!(runtime, time_s)
    modulation = _current_cascaded_h_bridge_modulation(runtime, time_s)
    timing = runtime.study.specification.timing
    cell_count = _cascaded_h_bridge_cell_count(runtime.study)
    changed = false
    for phase in 1:3
        first = (phase - 1) * 4 * cell_count + 1
        last = phase * 4 * cell_count
        requested = modulation.requested_valve_state[first:last]
        BridgeTopologies.bridge_topology_state_is_allowed(
            runtime.study.topologies[phase],
            requested,
        ) || throw(ArgumentError(
            "cascaded H-bridge modulation requested a prohibited phase state",
        ))
        for index in first:last
            valve = runtime.controlled_valves[index]
            desired = modulation.requested_valve_state[index]
            desired == valve.gate_driver.commanded_on && continue
            Nonlinear.request_power_semiconductor_gate!(
                valve,
                desired,
                time_s;
                earliest_transition_time_s=desired ?
                    time_s + timing.dead_time_s : time_s,
            )
            runtime.bridges[phase].transition_count += 1
            changed = true
        end
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :cascaded_h_bridge_modulation_command,
        :cascaded_h_bridge,
        true;
        message="accepted cascaded H-bridge cell-state command",
    ))
    return changed
end

function _synchronize_cascaded_h_bridge_conduction!(runtime, time_s)
    cell_count = _cascaded_h_bridge_cell_count(runtime.study)
    controlled_target = falses(length(runtime.controlled_valves))
    diode_target = falses(length(runtime.antiparallel_diodes))
    tolerance = maximum(
        valve.holding_current for valve in _cascaded_h_bridge_devices(runtime)
    )
    for phase in 1:3, cell in 1:cell_count
        index = (phase - 1) * cell_count + cell
        first = 4 * (index - 1) + 1
        applied = BitVector(
            runtime.controlled_valves[position].gate_driver.applied_on
            for position in first:(first + 3)
        )
        controlled, diodes = _cascaded_h_bridge_path_targets(
            applied,
            runtime.load_branches[phase].i_last,
            tolerance,
        )
        controlled_target[first:(first + 3)] .= controlled
        diode_target[first:(first + 3)] .= diodes
    end
    changed = false
    for index in eachindex(runtime.controlled_valves)
        changed |= _set_full_bridge_device_state!(
            runtime.controlled_valves[index],
            controlled_target[index] &&
                runtime.controlled_valves[index].gate_driver.applied_on,
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

function _stabilize_cascaded_h_bridge!(runtime, time_s, step_s)
    _cascaded_h_bridge_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed cascaded H-bridge requires the nonlinear solver",
    ))
    network = _cascaded_h_bridge_linear_network(runtime)
    foreach(bridge ->
        Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(bridge, time_s),
        runtime.bridges,
    )
    for _iteration in 1:24
        changed = _synchronize_cascaded_h_bridge_conduction!(runtime, time_s)
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        changed || return nothing
    end
    throw(ArgumentError("cascaded H-bridge topology failed to stabilize"))
end

function _record_cascaded_h_bridge!(recorder, runtime, sample, energy_residual)
    study = runtime.study
    network = _cascaded_h_bridge_linear_network(runtime)
    cell_count = _cascaded_h_bridge_cell_count(study)
    neutral_voltage = network.v[runtime.load_neutral_node]
    modulation = _current_cascaded_h_bridge_modulation(runtime, runtime.time_s)
    recorder.time_s[sample] = runtime.time_s
    recorder.requested_level[:, sample] .= modulation.requested_level
    recorder.requested_cell_state[:, :, sample] .= modulation.requested_cell_state
    for phase in 1:3
        recorder.phase_voltage_v[phase, sample] =
            network.v[runtime.phase_nodes[phase]] - neutral_voltage
        recorder.phase_current_a[phase, sample] = runtime.load_branches[phase].i_last
        for cell in 1:cell_count
            index = (phase - 1) * cell_count + cell
            capacitor = runtime.cell_capacitors[index]
            recorder.cell_dc_voltage_v[phase, cell, sample] =
                Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
            recorder.cell_dc_current_a[phase, cell, sample] = capacitor.i_last
        end
    end
    for index in eachindex(runtime.controlled_valves)
        controlled = runtime.controlled_valves[index]
        diode = runtime.antiparallel_diodes[index]
        recorder.requested_gate_state[index, sample] =
            controlled.gate_driver.commanded_on
        recorder.applied_gate_state[index, sample] = controlled.gate_driver.applied_on
        recorder.controlled_conducting_state[index, sample] = controlled.closed
        recorder.antiparallel_diode_conducting_state[index, sample] = diode.closed
        if _cascaded_h_bridge_is_detailed(study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            diode_state = Nonlinear.power_semiconductor_extended_state(diode)
            recorder.controlled_junction_temperature_k[index, sample] =
                controlled_state.junction_temperature_k
            recorder.antiparallel_junction_temperature_k[index, sample] =
                diode_state.junction_temperature_k
            recorder.antiparallel_recovered_charge_c[index, sample] =
                diode_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[index, sample] =
                controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_recovered_charge_c[index, sample] = 0.0
            recorder.controlled_turn_off_tail_current_a[index, sample] = 0.0
        end
    end
    power = _cascaded_h_bridge_power(runtime)
    recorder.stored_energy_j[sample] = _cascaded_h_bridge_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _cascaded_h_bridge_is_detailed(study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual
    recorder.write_index = sample
    return recorder
end

function prepare_switching_cascaded_h_bridge(
    study::ConverterSystems.SwitchingCascadedHBridgeStudy,
)
    runtime = _seed_cascaded_h_bridge_inverter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) /
            study.specification.timing.fixed_step_s,
    ) + 1
    integrator = CascadedHBridgeInverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        CascadedHBridgeInverterRecorder(
            _cascaded_h_bridge_cell_count(study),
            sample_count,
        ),
        0,
        false,
        false,
        nothing,
    )
    _record_cascaded_h_bridge!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_cascaded_h_bridge_inverter!(integrator)
    integrator.failed && throw(ArgumentError(
        "cascaded H-bridge integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _cascaded_h_bridge_power(runtime)
    previous_linear = _cascaded_h_bridge_linear_companion(runtime)
    previous_stored = _cascaded_h_bridge_stored_energy(runtime)
    previous_device_residual = _cascaded_h_bridge_device_residual(runtime)
    previous_device_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _cascaded_h_bridge_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_cascaded_h_bridge_gates!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        device_increment = 0.0
        if _cascaded_h_bridge_is_detailed(study)
            foreach(bridge ->
                Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                    bridge,
                    endpoint,
                ), runtime.bridges)
            conduction_changed = _synchronize_cascaded_h_bridge_conduction!(
                runtime,
                endpoint,
            )
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s, method, index) ->
                    _observe_cascaded_h_bridge_substep!(
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
            runtime.last_kcl_residual_a =
                nonlinear_result.diagnostics.maximum_kcl_residual_a
            final_device_dissipation = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in _cascaded_h_bridge_devices(runtime)
            )
            device_increment = final_device_dissipation - previous_device_dissipation
        else
            _stabilize_cascaded_h_bridge!(runtime, endpoint, step)
            network = _cascaded_h_bridge_linear_network(runtime)
            if gate_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint - half_step,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_cascaded_h_bridge_substep!(accounting, runtime, network,
                    endpoint - half_step, half_step,
                    Branches.BackwardEulerCompanion, 1)
                _stabilize_cascaded_h_bridge!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(
                    network,
                    endpoint,
                    half_step,
                    backward_euler,
                )
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_cascaded_h_bridge_substep!(accounting, runtime, network,
                    endpoint, half_step, Branches.BackwardEulerCompanion, 2)
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_cascaded_h_bridge_substep!(accounting, runtime, network,
                    endpoint, step, Branches.TrapezoidalCompanion, 1)
            end
            runtime.last_kcl_residual_a = maximum(
                abs,
                network.y * network.v - network.rhs;
                init=0.0,
            )
            device_increment = accounting.semiconductor_dissipated_energy_j
        end
        final_stored = _cascaded_h_bridge_stored_energy(runtime)
        device_residual_increment = _cascaded_h_bridge_device_residual(runtime) -
            previous_device_residual
        energy_residual = (
            -accounting.load_energy_j - device_increment -
            (final_stored - previous_stored) -
            accounting.linear_companion_energy_residual_j -
            device_residual_increment
        ) / step
        runtime.dissipated_energy_j += accounting.load_energy_j + device_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_cascaded_h_bridge!(
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

function _cascaded_h_bridge_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _cascaded_h_bridge_linear_network(runtime)
    cell_count = _cascaded_h_bridge_cell_count(study)
    cell_voltage = [
        Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        for capacitor in runtime.cell_capacitors
    ]
    phase_voltage = [
        network.v[node] - network.v[runtime.load_neutral_node]
        for node in runtime.phase_nodes
    ]
    phase_current = [branch.i_last for branch in runtime.load_branches]
    bridge_states = [
        Nonlinear.power_semiconductor_bridge_topology_state(
            bridge,
            network.v,
            study.specification.timing.fixed_step_s,
        ) for bridge in runtime.bridges
    ]
    stored_energy = _cascaded_h_bridge_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        (state.deterministic_signature for state in bridge_states)...,
        repr(runtime.time_s),
        repr(cell_voltage),
        repr(phase_voltage),
        repr(phase_current),
        string(integrator.accepted_step_index),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [cell_voltage..., phase_voltage...],
        [
            (capacitor.i_last for capacitor in runtime.cell_capacitors)...,
            phase_current...,
        ],
        BitVector(valve.gate_driver.commanded_on for valve in runtime.controlled_valves),
        BitVector(valve.gate_driver.applied_on for valve in runtime.controlled_valves),
        BitVector((
            getfield.(runtime.controlled_valves, :closed)...,
            getfield.(runtime.antiparallel_diodes, :closed)...,
        )),
        vec(permutedims(study.cell_dc_capacitance_f)) .* cell_voltage,
        study.load_inductance_h .* phase_current,
        Float64[],
        Float64.(recorder.requested_level[:, end]),
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    initial_stored = recorder.stored_energy_j[1]
    integrated_residual = -runtime.dissipated_energy_j -
        (stored_energy - initial_stored) -
        runtime.linear_companion_energy_residual_j -
        _cascaded_h_bridge_device_residual(runtime)
    energy_scale = max(
        abs(initial_stored),
        abs(stored_energy),
        abs(runtime.dissipated_energy_j),
        eps(Float64),
    )
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
    return ConverterSystems.SwitchingCascadedHBridgeTrace(
        recorder.time_s,
        recorder.cell_dc_voltage_v,
        recorder.cell_dc_current_a,
        recorder.phase_voltage_v,
        recorder.phase_current_a,
        recorder.requested_level,
        recorder.requested_cell_state,
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

function execute_switching_cascaded_h_bridge!(
    integrator::CascadedHBridgeInverterIntegrator,
)
    while _advance_cascaded_h_bridge_inverter!(integrator)
    end
    return _cascaded_h_bridge_result(integrator)
end

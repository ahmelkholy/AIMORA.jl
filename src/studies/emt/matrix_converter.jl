using ..ConverterSystems

mutable struct MatrixConverterRuntime{S,N,T,V,R,L}
    study::S
    network::N
    bridge::T
    valves::V
    sources::R
    load_branches::L
    input_nodes::NTuple{3,Int}
    output_nodes::NTuple{3,Int}
    load_neutral_node::Int
    requested_input_for_output::Vector{Int}
    applied_input_for_output::Vector{Int}
    commutation_outgoing_input::Vector{Int}
    commutation_incoming_input::Vector{Int}
    commutation_stage::Vector{UInt8}
    commutation_direction::Vector{Int8}
    time_s::Float64
    input_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct MatrixConverterRecorder
    time_s::Vector{Float64}
    input_source_voltage_v::Matrix{Float64}
    input_terminal_voltage_v::Matrix{Float64}
    input_current_a::Matrix{Float64}
    output_reference_voltage_v::Matrix{Float64}
    output_phase_voltage_v::Matrix{Float64}
    output_phase_current_a::Matrix{Float64}
    requested_connection::Array{Bool,3}
    applied_connection::Array{Bool,3}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    conducting_state::BitMatrix
    commutation_stage::Matrix{UInt8}
    commutation_direction::Matrix{Int8}
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    junction_temperature_k::Matrix{Float64}
    recovered_charge_c::Matrix{Float64}
    turn_off_tail_current_a::Matrix{Float64}
    write_index::Int
end

mutable struct MatrixConverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::MatrixConverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_matrix_converter_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_matrix_converter_linear_network(runtime) =
    _matrix_converter_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network

_matrix_converter_position(output::Int, input::Int, direction::Int) =
    2 * (3 * (output - 1) + input - 1) + direction

function _matrix_converter_nodes(study)
    nodes = Dict(node.name => node.node for node in study.topology.nodes)
    input = ntuple(phase -> nodes[Symbol(:input_, phase)], 3)
    output = ntuple(phase -> nodes[Symbol(:output_, phase)], 3)
    all(>(0), (input..., output...)) || throw(ArgumentError(
        "matrix-converter phase terminals must use positive node identities",
    ))
    input_port = only(port for port in study.specification.ports if
        port.identity === :input_ac)
    output_port = only(port for port in study.specification.ports if
        port.identity === :output_ac)
    input_port.ordered_nodes == input || throw(ArgumentError(
        "matrix-converter input port does not match its topology terminals",
    ))
    output_port.ordered_nodes == output || throw(ArgumentError(
        "matrix-converter output port does not match its topology terminals",
    ))
    return (; input, output, neutral=maximum((input..., output...)) + 1)
end

function _matrix_converter_source_voltage(study, time_s)
    angle = 2.0 * pi * study.input_frequency_hz * time_s
    return ntuple(phase -> study.input_phase_voltage_peak_v *
        sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
end

function _matrix_converter_modulation(study, time_s, input_voltage)
    modulation = study.specification.modulation
    return ConverterSystems.converter_matrix_space_vector_state(
        input_voltage,
        time_s,
        study.specification.rated_bases.frequency_hz,
        study.specification.timing.carrier_frequency_hz,
        modulation.modulation_index;
        output_phase_rad=modulation.phase_shift_rad,
    )
end

function _matrix_converter_connection(selected)
    connection = falses(3, 3)
    for output in 1:3
        connection[output, selected[output]] = true
    end
    return connection
end

function _matrix_converter_steady_gates(selected)
    gates = falses(18)
    for output in 1:3
        first = _matrix_converter_position(output, selected[output], 1)
        gates[first:(first + 1)] .= true
    end
    return gates
end

function _matrix_converter_electrical_input(runtime, output)
    stage = runtime.commutation_stage[output]
    return stage == 0x00 ? runtime.applied_input_for_output[output] :
        stage >= 0x04 ? runtime.commutation_incoming_input[output] :
        runtime.commutation_outgoing_input[output]
end

function _matrix_converter_safe_gates(runtime)
    gates = falses(18)
    for output in 1:3
        stage = Int(runtime.commutation_stage[output])
        if stage == 0
            first = _matrix_converter_position(
                output,
                runtime.applied_input_for_output[output],
                1,
            )
            gates[first:(first + 1)] .= true
            continue
        end
        direction_current = runtime.commutation_direction[output] > 0 ? 1.0 : -1.0
        sequence = ConverterSystems.matrix_safe_commutation_sequence(
            runtime.commutation_outgoing_input[output],
            runtime.commutation_incoming_input[output],
            direction_current,
        )
        local_state = sequence.stages[stage]
        for input in 1:3, direction in 1:2
            gates[_matrix_converter_position(output, input, direction)] =
                local_state[2 * (input - 1) + direction]
        end
    end
    return gates
end

function _matrix_converter_conduction_target(runtime)
    target = falses(18)
    tolerance = maximum(valve.holding_current for valve in runtime.valves)
    for output in 1:3
        current = runtime.load_branches[output].i_last
        direction = current < -tolerance ? 2 : 1
        input = _matrix_converter_electrical_input(runtime, output)
        target[_matrix_converter_position(output, input, direction)] = true
    end
    return target
end

function _seed_matrix_converter(study)
    nodes = _matrix_converter_nodes(study)
    source_voltage = _matrix_converter_source_voltage(study, study.start_time_s)
    initial_current = collect(study.initial_state.output_phase_current_a)
    selected = collect(study.initial_state.input_for_output)
    input_current = zeros(3)
    for output in 1:3
        input_current[selected[output]] += initial_current[output]
    end
    input_terminal_voltage = collect(source_voltage) .-
        study.source_resistance_ohm .* input_current
    raw_output_voltage = [input_terminal_voltage[selected[output]] for output in 1:3]
    output_common_mode = sum(raw_output_voltage) / 3.0
    output_phase_voltage = raw_output_voltage .- output_common_mode
    steady_gates = _matrix_converter_steady_gates(selected)
    conducting = falses(18)
    for output in 1:3
        direction = initial_current[output] < 0.0 ? 2 : 1
        conducting[_matrix_converter_position(output, selected[output], direction)] = true
    end
    valves = ntuple(18) do index
        _full_bridge_controlled_valve(
            study,
            study.topology.valve_positions[index],
            steady_gates[index],
            conducting[index],
        )
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, valves)
    bridge = Nonlinear.PowerSemiconductorBridgeTopology(
        study.topology,
        collect(valves),
    )
    sources = ntuple(3) do phase
        phase_shift = (phase - 1) * 2.0 * pi / 3.0
        Branches.TwoTerminalTheveninSource(
            nodes.input[phase],
            0,
            inv(study.source_resistance_ohm),
            time_s -> study.input_phase_voltage_peak_v *
                sin(2.0 * pi * study.input_frequency_hz * time_s - phase_shift),
        )
    end
    loads = ntuple(3) do phase
        Branches.SeriesRLBranch(
            nodes.output[phase],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            initial_current[phase],
            output_phase_voltage[phase],
            initial_current[phase],
        )
    end
    # Safe matrix commutation intentionally passes through directional gate
    # subsets that are not steady 3x3 connection states.  Stamp the canonical
    # directional valve positions directly while the runtime enforces the
    # five-stage interlock; the bridge descriptor remains the immutable
    # topology identity and steady-state admission owner.
    elements = _matrix_converter_is_detailed(study) ?
        Any[sources..., loads...] : Any[sources..., valves..., loads...]
    linear_network = Nodal.NodalSystem(nodes.neutral, elements)
    for phase in 1:3
        linear_network.v[nodes.input[phase]] = input_terminal_voltage[phase]
        linear_network.v[nodes.output[phase]] = raw_output_voltage[phase]
    end
    linear_network.v[nodes.neutral] = output_common_mode
    if _matrix_converter_is_detailed(study)
        for valve in valves
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _matrix_converter_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        valves;
        scales=NonlinearNetworkScales(
            nodes.neutral,
            0;
            nominal_voltage_v=study.input_phase_voltage_peak_v,
            nominal_current_a=max(maximum(abs, initial_current), 1.0),
        ),
    ) : linear_network
    return MatrixConverterRuntime(
        study,
        network,
        bridge,
        valves,
        sources,
        loads,
        nodes.input,
        nodes.output,
        nodes.neutral,
        copy(selected),
        copy(selected),
        copy(selected),
        copy(selected),
        zeros(UInt8, 3),
        ones(Int8, 3),
        study.start_time_s,
        0.0,
        0.0,
        0.0,
        0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function MatrixConverterRecorder(sample_count)
    sample_count > 1 || throw(ArgumentError(
        "matrix-converter recorder requires at least two samples",
    ))
    return MatrixConverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        falses(3, 3, sample_count),
        falses(3, 3, sample_count),
        falses(18, sample_count),
        falses(18, sample_count),
        falses(18, sample_count),
        zeros(UInt8, 3, sample_count),
        zeros(Int8, 3, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 18, sample_count),
        zeros(Float64, 18, sample_count),
        zeros(Float64, 18, sample_count),
        0,
    )
end

function _matrix_converter_stored_energy(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    _matrix_converter_is_detailed(runtime.study) || return energy
    return energy + sum(runtime.valves; init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _matrix_converter_power(runtime, time_s=runtime.time_s)
    study = runtime.study
    network = _matrix_converter_linear_network(runtime)
    source_voltage = _matrix_converter_source_voltage(study, time_s)
    input_current = ntuple(phase ->
        (source_voltage[phase] - network.v[runtime.input_nodes[phase]]) /
            study.source_resistance_ohm,
        3,
    )
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in runtime.valves
    )
    return (
        input_w=dot(source_voltage, input_current),
        load_w=study.load_resistance_ohm *
            sum(branch.i_last^2 for branch in runtime.load_branches),
        source_loss_w=study.source_resistance_ohm * sum(abs2, input_current),
        semiconductor_loss_w=semiconductor_loss,
        source_voltage_v=source_voltage,
        input_current_a=input_current,
    )
end

function _matrix_converter_linear_companion(runtime)
    network = _matrix_converter_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for branch in runtime.load_branches
        current = branch.i_last
        terminal_power += Branches.branch_voltage(network.v, branch.a, branch.b) * current
        dissipated_power += runtime.study.load_resistance_ohm * current^2
        stored_energy += 0.5 * runtime.study.load_inductance_h * current^2
    end
    return (; terminal_power_w=terminal_power,
        dissipated_power_w=dissipated_power, stored_energy_j=stored_energy)
end

function _matrix_converter_device_companion_residual(runtime)
    _matrix_converter_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in runtime.valves
    )
end

function _observe_matrix_converter_substep!(
    accounting::FullBridgeStepEnergyAccounting,
    runtime,
    _system,
    time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _matrix_converter_power(runtime, time_s)
    final_linear = _matrix_converter_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("matrix converter received an unknown companion method"))
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
        step_s * blend(previous_linear.dissipated_power_w,
            final_linear.dissipated_power_w) -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _update_matrix_converter_commutation!(runtime, endpoint)
    network = _matrix_converter_linear_network(runtime)
    input_voltage = ntuple(phase -> network.v[runtime.input_nodes[phase]], 3)
    modulation = _matrix_converter_modulation(runtime.study, endpoint, input_voltage)
    runtime.requested_input_for_output .= modulation.requested_input_for_output
    changed = false
    for output in 1:3
        stage = runtime.commutation_stage[output]
        if stage == 0x05
            runtime.applied_input_for_output[output] =
                runtime.commutation_incoming_input[output]
            runtime.commutation_stage[output] = 0x00
            push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                endpoint,
                :matrix_safe_commutation_complete,
                Symbol(:matrix_output_, output),
                true;
                message="accepted current-direction-dependent five-stage commutation",
            ))
            changed = true
        elseif stage != 0x00
            runtime.commutation_stage[output] += 0x01
            changed = true
        end
        if runtime.commutation_stage[output] == 0x00 &&
           runtime.requested_input_for_output[output] !=
                runtime.applied_input_for_output[output]
            runtime.commutation_outgoing_input[output] =
                runtime.applied_input_for_output[output]
            runtime.commutation_incoming_input[output] =
                runtime.requested_input_for_output[output]
            current = runtime.load_branches[output].i_last
            runtime.commutation_direction[output] = current < 0.0 ? Int8(-1) : Int8(1)
            runtime.commutation_stage[output] = 0x01
            push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                endpoint,
                :matrix_safe_commutation_start,
                Symbol(:matrix_output_, output),
                true;
                message="accepted requested matrix switching-vector transition",
            ))
            changed = true
        end
    end
    requested_gates = _matrix_converter_safe_gates(runtime)
    for index in 1:18
        changed |= Nonlinear.request_power_semiconductor_gate!(
            runtime.valves[index],
            requested_gates[index],
            endpoint,
        )
        transition_time = Nonlinear.power_semiconductor_gate_transition_time(
            runtime.valves[index],
        )
        if transition_time !== nothing && transition_time <= endpoint +
                16.0 * eps(max(1.0, endpoint))
            changed |= Nonlinear.apply_power_semiconductor_gate_transition!(
                runtime.valves[index],
                endpoint,
            )
        end
    end
    conduction = _matrix_converter_conduction_target(runtime)
    for index in 1:18
        changed |= _set_full_bridge_device_state!(
            runtime.valves[index],
            conduction[index],
            endpoint,
        )
    end
    return modulation, changed
end

function _record_matrix_converter!(recorder, runtime, sample, energy_residual_w)
    network = _matrix_converter_linear_network(runtime)
    power = _matrix_converter_power(runtime)
    modulation = _matrix_converter_modulation(
        runtime.study,
        runtime.time_s,
        ntuple(phase -> network.v[runtime.input_nodes[phase]], 3),
    )
    neutral_voltage = network.v[runtime.load_neutral_node]
    requested = _matrix_converter_connection(runtime.requested_input_for_output)
    applied = _matrix_converter_connection([
        _matrix_converter_electrical_input(runtime, output) for output in 1:3
    ])
    recorder.time_s[sample] = runtime.time_s
    recorder.input_source_voltage_v[:, sample] .= power.source_voltage_v
    recorder.input_terminal_voltage_v[:, sample] .=
        network.v[collect(runtime.input_nodes)]
    recorder.input_current_a[:, sample] .= power.input_current_a
    recorder.output_reference_voltage_v[:, sample] .=
        modulation.output_reference_voltage_v
    for output in 1:3
        recorder.output_phase_voltage_v[output, sample] =
            network.v[runtime.output_nodes[output]] - neutral_voltage
        recorder.output_phase_current_a[output, sample] =
            runtime.load_branches[output].i_last
    end
    recorder.requested_connection[:, :, sample] .= requested
    recorder.applied_connection[:, :, sample] .= applied
    safe_gates = _matrix_converter_safe_gates(runtime)
    for index in 1:18
        valve = runtime.valves[index]
        recorder.requested_gate_state[index, sample] = safe_gates[index]
        recorder.applied_gate_state[index, sample] = valve.gate_driver.applied_on
        recorder.conducting_state[index, sample] = valve.closed
        if _matrix_converter_is_detailed(runtime.study)
            extended = Nonlinear.power_semiconductor_extended_state(valve)
            recorder.junction_temperature_k[index, sample] =
                extended.junction_temperature_k
            recorder.recovered_charge_c[index, sample] =
                extended.stored_recovery_charge_c
            recorder.turn_off_tail_current_a[index, sample] = extended.tail_current_a
        else
            recorder.junction_temperature_k[index, sample] = 0.0
            recorder.recovered_charge_c[index, sample] = 0.0
            recorder.turn_off_tail_current_a[index, sample] = 0.0
        end
    end
    recorder.commutation_stage[:, sample] .= runtime.commutation_stage
    recorder.commutation_direction[:, sample] .= runtime.commutation_direction
    recorder.stored_energy_j[sample] = _matrix_converter_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _matrix_converter_is_detailed(runtime.study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual_w
    recorder.write_index = sample
    return recorder
end

function prepare_switching_matrix_converter(
    study::ConverterSystems.SwitchingMatrixConverterStudy,
)
    runtime = _seed_matrix_converter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) /
            study.specification.timing.fixed_step_s,
    ) + 1
    integrator = MatrixConverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        MatrixConverterRecorder(sample_count),
        0,
        false,
        false,
        nothing,
    )
    _record_matrix_converter!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_matrix_converter!(integrator)
    integrator.failed && throw(ArgumentError(
        "matrix-converter integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _matrix_converter_power(runtime)
    previous_linear = _matrix_converter_linear_companion(runtime)
    previous_stored = _matrix_converter_stored_energy(runtime)
    previous_device_residual = _matrix_converter_device_companion_residual(runtime)
    previous_device_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in runtime.valves
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        _, topology_changed = _update_matrix_converter_commutation!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        device_increment = 0.0
        if _matrix_converter_is_detailed(study)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=topology_changed ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=topology_changed ? :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) ->
                    _observe_matrix_converter_substep!(
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
                for valve in runtime.valves
            )
            device_increment = final_device_dissipation - previous_device_dissipation
        else
            network = _matrix_converter_linear_network(runtime)
            if topology_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                for (time_s, substep_index) in (
                    (endpoint - half_step, 1),
                    (endpoint, 2),
                )
                    Nodal.solve_algebraic_state!(network, time_s, half_step, backward_euler)
                    Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                    _observe_matrix_converter_substep!(accounting, runtime, network,
                        time_s, half_step, Branches.BackwardEulerCompanion,
                        substep_index)
                end
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_matrix_converter_substep!(accounting, runtime, network,
                    endpoint, step, Branches.TrapezoidalCompanion, 1)
            end
            runtime.last_kcl_residual_a =
                maximum(abs, network.y * network.v - network.rhs; init=0.0)
            device_increment = accounting.semiconductor_dissipated_energy_j
        end
        final_stored = _matrix_converter_stored_energy(runtime)
        device_residual_increment =
            _matrix_converter_device_companion_residual(runtime) -
            previous_device_residual
        energy_residual = (
            accounting.input_energy_j - accounting.load_energy_j -
            accounting.source_dissipated_energy_j - device_increment -
            (final_stored - previous_stored) -
            accounting.linear_companion_energy_residual_j -
            device_residual_increment
        ) / step
        runtime.input_energy_j += accounting.input_energy_j
        runtime.dissipated_energy_j += accounting.load_energy_j +
            accounting.source_dissipated_energy_j + device_increment
        runtime.linear_companion_energy_residual_j +=
            accounting.linear_companion_energy_residual_j
        runtime.time_s = endpoint
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_matrix_converter!(integrator.recorder, runtime, next_step + 1,
            energy_residual)
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _matrix_converter_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _matrix_converter_linear_network(runtime)
    power = _matrix_converter_power(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    output_voltage = [network.v[node] - neutral_voltage for node in runtime.output_nodes]
    output_current = [branch.i_last for branch in runtime.load_branches]
    stored_energy = _matrix_converter_stored_energy(runtime)
    requested_gates = _matrix_converter_safe_gates(runtime)
    applied_gates = BitVector(valve.gate_driver.applied_on for valve in runtime.valves)
    conduction = BitVector(valve.closed for valve in runtime.valves)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(runtime.requested_input_for_output),
        repr(runtime.applied_input_for_output),
        repr(runtime.commutation_stage),
        repr(output_voltage),
        repr(output_current),
        repr(stored_energy),
        string(integrator.accepted_step_index),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [network.v[collect(runtime.input_nodes)]..., output_voltage...],
        [power.input_current_a..., output_current...],
        requested_gates,
        applied_gates,
        conduction,
        Float64[],
        study.load_inductance_h .* output_current,
        Float64[],
        Float64.(runtime.requested_input_for_output),
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
        runtime.linear_companion_energy_residual_j -
        _matrix_converter_device_companion_residual(runtime)
    result = ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=integrator.completed && !integrator.failed,
        status=integrator.completed && !integrator.failed ? :ok : :incomplete,
        events=runtime.events,
        maximum_kcl_residual_a=maximum(recorder.kcl_residual_a),
        relative_charge_residual=maximum(abs,
            vec(sum(recorder.output_phase_current_a; dims=1))),
        relative_energy_residual=abs(integrated_residual) / energy_scale,
    )
    return ConverterSystems.SwitchingMatrixConverterTrace(
        recorder.time_s,
        recorder.input_source_voltage_v,
        recorder.input_terminal_voltage_v,
        recorder.input_current_a,
        recorder.output_reference_voltage_v,
        recorder.output_phase_voltage_v,
        recorder.output_phase_current_a,
        recorder.requested_connection,
        recorder.applied_connection,
        recorder.requested_gate_state,
        recorder.applied_gate_state,
        recorder.conducting_state,
        recorder.commutation_stage,
        recorder.commutation_direction,
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

function execute_switching_matrix_converter!(integrator::MatrixConverterIntegrator)
    while _advance_matrix_converter!(integrator)
    end
    return _matrix_converter_result(integrator)
end

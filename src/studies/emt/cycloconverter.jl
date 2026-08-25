using ..ConverterSystems

mutable struct CycloconverterRuntime{S,N,T,V,R,L}
    study::S
    network::N
    bridges::T
    valves::V
    sources::R
    load_branches::L
    input_nodes::Vector{NTuple{3,Int}}
    source_neutral_nodes::Vector{Int}
    output_nodes::Vector{Int}
    output_neutral_node::Int
    requested_bridge_group::Vector{Int8}
    active_bridge_group::Vector{Int8}
    previous_firing_state::BitVector
    failed_transfer_latched::BitVector
    time_s::Float64
    source_energy_j::Float64
    dissipated_energy_j::Float64
    linear_companion_energy_residual_j::Float64
    last_kcl_residual_a::Float64
    events::Vector{ConverterSystems.ConverterSystemEventRecord}
end

mutable struct CycloconverterRecorder
    time_s::Vector{Float64}
    source_voltage_v::Matrix{Float64}
    source_current_a::Matrix{Float64}
    output_reference_voltage_v::Matrix{Float64}
    output_voltage_v::Matrix{Float64}
    output_current_a::Matrix{Float64}
    firing_angle_rad::Matrix{Float64}
    requested_bridge_group::Matrix{Int8}
    active_bridge_group::Matrix{Int8}
    requested_firing_state::BitMatrix
    applied_firing_state::BitMatrix
    conducting_state::BitMatrix
    commutation_overlap_valve_count::Matrix{UInt8}
    failed_commutation::BitMatrix
    reactive_power_var::Vector{Float64}
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    junction_temperature_k::Matrix{Float64}
    recovered_charge_c::Matrix{Float64}
    write_index::Int
end

mutable struct CycloconverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::CycloconverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_cycloconverter_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_cycloconverter_linear_network(runtime) =
    _cycloconverter_is_detailed(runtime.study) ?
        nonlinear_linear_system(runtime.network) : runtime.network

function _cycloconverter_nodes(study)
    phase_count = study.specification.selection.phase_count
    input = NTuple{3,Int}[]
    output = Int[]
    neutrals = Int[]
    for topology in study.topologies
        names = Dict(node.name => node.node for node in topology.nodes)
        push!(input, ntuple(phase -> names[Symbol(:input_, phase)], 3))
        push!(output, names[:output_1])
        push!(neutrals, names[:output_neutral])
    end
    length(unique(neutrals)) == 1 || throw(ArgumentError(
        "cycloconverter output topologies must share one exact load neutral",
    ))
    neutral = only(unique(neutrals))
    all(>(0), (Iterators.flatten(input)..., output...)) || throw(ArgumentError(
        "cycloconverter input and output phase terminals must be positive nodes",
    ))
    length(unique((Iterators.flatten(input)..., output...))) == 3 * phase_count +
        phase_count || throw(ArgumentError(
            "cycloconverter transformer-secondary and output phase nodes must be distinct",
        ))
    input_port = only(port for port in study.specification.ports if
        port.identity === :input_ac)
    output_port = only(port for port in study.specification.ports if
        port.identity === :output_ac)
    input_port.ordered_nodes == Tuple(Iterators.flatten(input)) || throw(ArgumentError(
        "cycloconverter input port does not match its ordered transformer-secondary terminals",
    ))
    output_port.ordered_nodes == Tuple((output..., neutral)) || throw(ArgumentError(
        "cycloconverter output port must bind every phase followed by its neutral",
    ))
    return (; input, output, neutral)
end

function _cycloconverter_source_voltage(study, time_s)
    angle = 2.0 * pi * study.input_frequency_hz * time_s
    return ntuple(phase -> study.input_phase_voltage_peak_v *
        sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
end

function _cycloconverter_reference(study, time_s)
    phase_count = study.specification.selection.phase_count
    modulation = study.specification.modulation
    output_frequency = study.specification.rated_bases.frequency_hz
    maximum_average_voltage = 3.0 * sqrt(3.0) / pi *
        study.input_phase_voltage_peak_v
    return [modulation.modulation_index * maximum_average_voltage *
        sin(2.0 * pi * output_frequency * time_s + modulation.phase_shift_rad -
            (output - 1) * 2.0 * pi / 3.0)
        for output in 1:phase_count]
end

function _cycloconverter_requested_group(reference, current, previous, tolerance)
    abs(current) > tolerance && return Int8(sign(current))
    abs(reference) <= 64.0 * eps(max(abs(reference), 1.0)) && return
        (previous == 0 ? Int8(1) : previous)
    return Int8(sign(reference))
end

function _cycloconverter_firing_state(study, time_s, active_group)
    phase_count = study.specification.selection.phase_count
    reference = _cycloconverter_reference(study, time_s)
    maximum_average_voltage = 3.0 * sqrt(3.0) / pi *
        study.input_phase_voltage_peak_v
    angles = [begin
        group = active_group[output] == 0 ?
            (reference[output] < 0.0 ? Int8(-1) : Int8(1)) :
            active_group[output]
        acos(clamp(group * reference[output] / maximum_average_voltage, -1.0, 1.0))
    end for output in eachindex(reference)]
    state = falses(12 * phase_count)
    for output in 1:phase_count
        group = active_group[output]
        group == 0 && continue
        delayed_time = time_s - angles[output] /
            (2.0 * pi * study.input_frequency_hz)
        delayed_voltage = _cycloconverter_source_voltage(study, delayed_time)
        upper = argmax(delayed_voltage)
        lower = argmin(delayed_voltage)
        base = 12 * (output - 1) + (group > 0 ? 0 : 6)
        state[base + 2 * upper - 1] = true
        state[base + 2 * lower] = true
    end
    return (; state, reference, firing_angle_rad=angles)
end

function _seed_cycloconverter(study)
    nodes = _cycloconverter_nodes(study)
    phase_count = study.specification.selection.phase_count
    current = study.initial_state.output_current_a
    active_group = copy(study.initial_state.active_bridge_group)
    reference = _cycloconverter_reference(study, study.start_time_s)
    for output in 1:phase_count
        active_group[output] == 0 && !iszero(current[output]) &&
            (active_group[output] = Int8(sign(current[output])))
    end
    firing = _cycloconverter_firing_state(study, study.start_time_s, active_group)
    positions = tuple((position for topology in study.topologies for position in
        topology.valve_positions)...)
    valves = ntuple(length(positions)) do index
        _rectifier_valve(study, positions[index], firing.state[index])
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, valves)
    bridges = ntuple(phase_count) do output
        first = 12 * (output - 1) + 1
        Nonlinear.PowerSemiconductorBridgeTopology(
            study.topologies[output],
            collect(valves[first:(first + 11)]),
        )
    end
    maximum_external_node = maximum((Iterators.flatten(nodes.input)...,
        nodes.output..., nodes.neutral))
    source_neutral_nodes = vcat(0, collect((maximum_external_node + 1):
        (maximum_external_node + phase_count - 1)))
    sources = ntuple(3 * phase_count) do channel
        output = (channel - 1) ÷ 3 + 1
        phase = (channel - 1) % 3 + 1
        phase_shift = (phase - 1) * 2.0 * pi / 3.0
        Branches.TwoTerminalTheveninSource(
            nodes.input[output][phase],
            source_neutral_nodes[output],
            inv(study.source_resistance_ohm),
            time_s -> study.input_phase_voltage_peak_v *
                sin(2.0 * pi * study.input_frequency_hz * time_s - phase_shift),
        )
    end
    source_voltage = _cycloconverter_source_voltage(study, study.start_time_s)
    output_voltage = zeros(phase_count)
    for output in 1:phase_count
        local_positions = findall(firing.state[
            (12 * (output - 1) + 1):(12 * output)])
        isempty(local_positions) && continue
        upper_position = only(index for index in local_positions if isodd(index))
        lower_position = only(index for index in local_positions if iseven(index))
        upper_node = study.topologies[output].valve_positions[upper_position].from_node
        lower_node = study.topologies[output].valve_positions[lower_position].to_node
        source_by_node = Dict(zip(nodes.input[output], source_voltage))
        output_voltage[output] = get(source_by_node, upper_node, 0.0) -
            get(source_by_node, lower_node, 0.0)
    end
    loads = ntuple(phase_count) do output
        Branches.SeriesRLBranch(
            nodes.output[output],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            current[output],
            output_voltage[output],
            current[output],
        )
    end
    node_count = maximum((Iterators.flatten(nodes.input)..., source_neutral_nodes...,
        nodes.output..., nodes.neutral))
    linear_network = Nodal.NodalSystem(node_count, Any[sources..., bridges..., loads...])
    for output in 1:phase_count, phase in 1:3
        linear_network.v[nodes.input[output][phase]] = source_voltage[phase]
    end
    for node in source_neutral_nodes
        node == 0 || (linear_network.v[node] = 0.0)
    end
    nodes.neutral == 0 || (linear_network.v[nodes.neutral] = 0.0)
    for output in 1:phase_count
        linear_network.v[nodes.output[output]] = output_voltage[output]
    end
    if _cycloconverter_is_detailed(study)
        for valve in valves
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _cycloconverter_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        valves;
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_phase_voltage_peak_v,
            nominal_current_a=max(maximum(abs, current),
                study.input_phase_voltage_peak_v / study.load_resistance_ohm, 1.0),
        ),
    ) : linear_network
    requested_group = [_cycloconverter_requested_group(reference[output],
        current[output], active_group[output], study.current_zero_tolerance_a)
        for output in 1:phase_count]
    return CycloconverterRuntime(
        study, network, bridges, valves, sources, loads, nodes.input,
        source_neutral_nodes, nodes.output, nodes.neutral, requested_group,
        active_group, firing.state,
        falses(phase_count), study.start_time_s, 0.0, 0.0, 0.0, 0.0,
        ConverterSystems.ConverterSystemEventRecord[],
    )
end

function CycloconverterRecorder(phase_count, valve_count, sample_count)
    return CycloconverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, 3 * phase_count, sample_count),
        zeros(Float64, 3 * phase_count, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Float64, phase_count, sample_count),
        zeros(Int8, phase_count, sample_count),
        zeros(Int8, phase_count, sample_count),
        falses(valve_count, sample_count),
        falses(valve_count, sample_count),
        falses(valve_count, sample_count),
        zeros(UInt8, phase_count, sample_count),
        falses(phase_count, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, valve_count, sample_count),
        zeros(Float64, valve_count, sample_count),
        0,
    )
end

function _cycloconverter_source_currents(runtime, time_s=runtime.time_s)
    study = runtime.study
    network = _cycloconverter_linear_network(runtime)
    source_voltage = _cycloconverter_source_voltage(study, time_s)
    return tuple((
        (source_voltage[phase] -
            (network.v[runtime.input_nodes[output][phase]] -
                (runtime.source_neutral_nodes[output] == 0 ? 0.0 :
                    network.v[runtime.source_neutral_nodes[output]]))) /
            study.source_resistance_ohm
        for output in eachindex(runtime.input_nodes) for phase in 1:3
    )...)
end

function _cycloconverter_stored_energy(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    _cycloconverter_is_detailed(runtime.study) || return energy
    return energy + sum(runtime.valves; init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _cycloconverter_linear_companion(runtime)
    network = _cycloconverter_linear_network(runtime)
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

function _cycloconverter_device_companion_residual(runtime)
    _cycloconverter_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in runtime.valves
    )
end

function _cycloconverter_power(runtime, time_s=runtime.time_s)
    study = runtime.study
    source_voltage = _cycloconverter_source_voltage(study, time_s)
    source_current = _cycloconverter_source_currents(runtime, time_s)
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in runtime.valves
    )
    return (
        source_w=sum(source_voltage[phase] *
            source_current[3 * (output - 1) + phase]
            for output in eachindex(runtime.input_nodes) for phase in 1:3),
        source_loss_w=study.source_resistance_ohm * sum(abs2, source_current),
        load_loss_w=study.load_resistance_ohm *
            sum(branch.i_last^2 for branch in runtime.load_branches),
        semiconductor_loss_w=semiconductor_loss,
        source_voltage_v=tuple((source_voltage[phase]
            for _output in eachindex(runtime.input_nodes) for phase in 1:3)...),
        source_current_a=source_current,
    )
end

function _observe_cycloconverter_substep!(
    accounting::LineCommutatedRectifierStepEnergyAccounting,
    runtime,
    _system,
    time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _cycloconverter_power(runtime, time_s)
    final_linear = _cycloconverter_linear_companion(runtime)
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("cycloconverter received an unknown companion method"))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    accounting.source_energy_j += step_s *
        blend(accounting.previous_power.source_w, final_power.source_w)
    accounting.source_dissipated_energy_j += step_s *
        blend(accounting.previous_power.source_loss_w, final_power.source_loss_w)
    accounting.load_dissipated_energy_j += step_s *
        blend(accounting.previous_power.load_loss_w, final_power.load_loss_w)
    accounting.semiconductor_dissipated_energy_j += step_s *
        blend(accounting.previous_power.semiconductor_loss_w,
            final_power.semiconductor_loss_w)
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

function _update_cycloconverter_firing!(runtime, endpoint)
    study = runtime.study
    reference = _cycloconverter_reference(study, endpoint)
    failed = falses(length(runtime.load_branches))
    for output in eachindex(runtime.load_branches)
        current = runtime.load_branches[output].i_last
        requested = _cycloconverter_requested_group(reference[output], current,
            runtime.requested_bridge_group[output], study.current_zero_tolerance_a)
        runtime.requested_bridge_group[output] = requested
        active = runtime.active_bridge_group[output]
        active == requested && begin
            runtime.failed_transfer_latched[output] = false
            continue
        end
        if abs(current) <= study.current_zero_tolerance_a ||
           (active != 0 && active * current < 0.0)
            runtime.active_bridge_group[output] = requested
            runtime.failed_transfer_latched[output] = false
            push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                endpoint,
                :cycloconverter_group_transfer,
                Symbol(:cycloconverter_output_, output),
                true;
                message="accepted noncirculating bridge-group transfer at a resolved or bracketed current zero",
            ))
        else
            failed[output] = true
            if !runtime.failed_transfer_latched[output]
                push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
                    endpoint,
                    :cycloconverter_commutation_deferred,
                    Symbol(:cycloconverter_output_, output),
                    false;
                    message="opposing bridge group remained blocked away from current zero",
                ))
            end
            runtime.failed_transfer_latched[output] = true
        end
    end
    firing = _cycloconverter_firing_state(study, endpoint,
        runtime.active_bridge_group)
    changed = firing.state != runtime.previous_firing_state
    for (index, valve) in enumerate(runtime.valves)
        Nonlinear.request_power_semiconductor_gate!(valve, firing.state[index], endpoint)
        previous = valve.closed
        if firing.state[index] && !previous
            Nonlinear.apply_power_semiconductor_forward_turn_on!(valve, endpoint)
        elseif !firing.state[index] && previous
            Nonlinear.apply_power_semiconductor_forward_extinction!(valve, endpoint)
        end
        previous == valve.closed || push!(runtime.events,
            ConverterSystems.ConverterSystemEventRecord(
                endpoint,
                valve.closed ? :thyristor_firing : :natural_current_extinction,
                study.topologies[(index - 1) ÷ 12 + 1].valve_positions[
                    (index - 1) % 12 + 1].name,
                true;
                message="accepted cycloconverter line-commutated valve transition",
            ))
    end
    runtime.previous_firing_state .= firing.state
    return firing, failed, changed
end

function _cycloconverter_reactive_power(source_voltage, source_current)
    v_alpha = 2.0 / 3.0 * (source_voltage[1] -
        0.5 * source_voltage[2] - 0.5 * source_voltage[3])
    v_beta = 2.0 / 3.0 * (sqrt(3.0) / 2.0 *
        (source_voltage[2] - source_voltage[3]))
    i_alpha = 2.0 / 3.0 * (source_current[1] -
        0.5 * source_current[2] - 0.5 * source_current[3])
    i_beta = 2.0 / 3.0 * (sqrt(3.0) / 2.0 *
        (source_current[2] - source_current[3]))
    return 1.5 * (v_beta * i_alpha - v_alpha * i_beta)
end

function _record_cycloconverter!(recorder, runtime, sample, firing, failed,
    energy_residual_w)
    network = _cycloconverter_linear_network(runtime)
    power = _cycloconverter_power(runtime)
    neutral_voltage = runtime.output_neutral_node == 0 ? 0.0 :
        network.v[runtime.output_neutral_node]
    recorder.time_s[sample] = runtime.time_s
    recorder.source_voltage_v[:, sample] .= power.source_voltage_v
    recorder.source_current_a[:, sample] .= power.source_current_a
    recorder.output_reference_voltage_v[:, sample] .= firing.reference
    recorder.firing_angle_rad[:, sample] .= firing.firing_angle_rad
    for output in eachindex(runtime.load_branches)
        recorder.output_voltage_v[output, sample] =
            network.v[runtime.output_nodes[output]] - neutral_voltage
        recorder.output_current_a[output, sample] =
            runtime.load_branches[output].i_last
        recorder.requested_bridge_group[output, sample] =
            runtime.requested_bridge_group[output]
        recorder.active_bridge_group[output, sample] =
            runtime.active_bridge_group[output]
        local_range = (12 * (output - 1) + 1):(12 * output)
        conducting_count = count(valve.closed for valve in runtime.valves[local_range])
        recorder.commutation_overlap_valve_count[output, sample] =
            UInt8(max(conducting_count - 2, 0))
        recorder.failed_commutation[output, sample] = failed[output]
    end
    for (index, valve) in enumerate(runtime.valves)
        recorder.requested_firing_state[index, sample] = valve.gate_driver.commanded_on
        recorder.applied_firing_state[index, sample] = valve.gate_driver.applied_on
        recorder.conducting_state[index, sample] = valve.closed
        if _cycloconverter_is_detailed(runtime.study)
            extended = Nonlinear.power_semiconductor_extended_state(valve)
            recorder.junction_temperature_k[index, sample] =
                extended.junction_temperature_k
            recorder.recovered_charge_c[index, sample] =
                extended.stored_recovery_charge_c
        else
            recorder.junction_temperature_k[index, sample] = 0.0
            recorder.recovered_charge_c[index, sample] = 0.0
        end
    end
    recorder.reactive_power_var[sample] = sum(
        _cycloconverter_reactive_power(
            power.source_voltage_v[(3 * (output - 1) + 1):(3 * output)],
            power.source_current_a[(3 * (output - 1) + 1):(3 * output)],
        ) for output in eachindex(runtime.input_nodes)
    )
    recorder.stored_energy_j[sample] = _cycloconverter_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _cycloconverter_is_detailed(runtime.study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual_w
    recorder.write_index = sample
    return recorder
end

function prepare_switching_cycloconverter(
    study::ConverterSystems.SwitchingCycloconverterStudy,
)
    runtime = _seed_cycloconverter(study)
    sample_count = round(Int, (study.stop_time_s - study.start_time_s) /
        study.specification.timing.fixed_step_s) + 1
    integrator = CycloconverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        CycloconverterRecorder(length(runtime.load_branches),
            length(runtime.valves), sample_count),
        0, false, false, nothing,
    )
    firing = _cycloconverter_firing_state(study, study.start_time_s,
        runtime.active_bridge_group)
    _record_cycloconverter!(integrator.recorder, runtime, 1, firing,
        falses(length(runtime.load_branches)), 0.0)
    return integrator
end

function _advance_cycloconverter!(integrator)
    integrator.failed && throw(ArgumentError(
        "cycloconverter integrator is terminally failed: $(integrator.last_failure)",
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
    previous_energy = _cycloconverter_stored_energy(runtime)
    previous_power = _cycloconverter_power(runtime)
    previous_linear = _cycloconverter_linear_companion(runtime)
    previous_device_residual = _cycloconverter_device_companion_residual(runtime)
    previous_semiconductor_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in runtime.valves
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        firing, failed, topology_changed =
            _update_cycloconverter_firing!(runtime, endpoint)
        accounting = LineCommutatedRectifierStepEnergyAccounting(previous_power,
            previous_linear)
        semiconductor_increment = 0.0
        if _cycloconverter_is_detailed(study)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network, endpoint, step;
                discontinuity_treatment=topology_changed ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=topology_changed ? :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) ->
                    _observe_cycloconverter_substep!(accounting, runtime, system,
                        time_s, substep_s, companion_method, substep_index),
            )
            nonlinear_result.accepted || throw(something(nonlinear_result.failure))
            runtime.last_kcl_residual_a =
                nonlinear_result.diagnostics.maximum_kcl_residual_a
            final_semiconductor_dissipation = sum(
                valve.semiconductor_dissipated_energy_j +
                    (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
                for valve in runtime.valves
            )
            semiconductor_increment = final_semiconductor_dissipation -
                previous_semiconductor_dissipation
        else
            network = _cycloconverter_linear_network(runtime)
            if topology_changed
                half_step = step / 2.0
                backward_euler = Val(Branches.BackwardEulerCompanion)
                for (time_s, substep_index) in ((endpoint - half_step, 1), (endpoint, 2))
                    Nodal.solve_algebraic_state!(network, time_s, half_step,
                        backward_euler)
                    Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                    _observe_cycloconverter_substep!(accounting, runtime, network,
                        time_s, half_step, Branches.BackwardEulerCompanion,
                        substep_index)
                end
            else
                Nodal.solve_algebraic_state!(network, endpoint, step)
                Nodal.accept_algebraic_state!(network, step)
                _observe_cycloconverter_substep!(accounting, runtime, network,
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
        stored_energy = _cycloconverter_stored_energy(runtime)
        energy_residual = (accounting.source_energy_j - dissipation_increment -
            (stored_energy - previous_energy) -
            accounting.linear_companion_energy_residual_j -
            (_cycloconverter_device_companion_residual(runtime) -
                previous_device_residual)) / step
        commit_timestep_transaction!(integrator.transaction)
        integrator.accepted_step_index = next_step
        integrator.completed = next_step == final_step
        _record_cycloconverter!(integrator.recorder, runtime, next_step + 1,
            firing, failed, energy_residual)
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _cycloconverter_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _cycloconverter_linear_network(runtime)
    power = _cycloconverter_power(runtime)
    neutral_voltage = runtime.output_neutral_node == 0 ? 0.0 :
        network.v[runtime.output_neutral_node]
    output_voltage = [network.v[node] - neutral_voltage for node in runtime.output_nodes]
    output_current = [branch.i_last for branch in runtime.load_branches]
    stored_energy = _cycloconverter_stored_energy(runtime)
    signature = bytes2hex(sha256(join((study.specification.signature_sha256,
        repr(runtime.time_s), repr(runtime.requested_bridge_group),
        repr(runtime.active_bridge_group), repr(output_voltage),
        repr(output_current), repr(stored_energy),
        string(integrator.accepted_step_index)), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [power.source_voltage_v..., output_voltage...],
        [power.source_current_a..., output_current...],
        BitVector(valve.gate_driver.commanded_on for valve in runtime.valves),
        BitVector(valve.gate_driver.applied_on for valve in runtime.valves),
        BitVector(valve.closed for valve in runtime.valves),
        Float64[], study.load_inductance_h .* output_current, Float64[],
        Float64.(runtime.active_bridge_group), stored_energy,
        runtime.dissipated_energy_j, integrator.accepted_step_index,
        length(runtime.events), signature,
    )
    stored_change = stored_energy - recorder.stored_energy_j[1]
    integrated_residual = runtime.source_energy_j - runtime.dissipated_energy_j -
        stored_change - runtime.linear_companion_energy_residual_j -
        _cycloconverter_device_companion_residual(runtime)
    scale = max(abs(runtime.source_energy_j), abs(runtime.dissipated_energy_j),
        abs(stored_change), eps(Float64))
    result = ConverterSystems.converter_system_result(
        study.specification, state;
        accepted=integrator.completed && !integrator.failed,
        status=integrator.completed && !integrator.failed ? :ok : :incomplete,
        events=runtime.events,
        maximum_kcl_residual_a=maximum(recorder.kcl_residual_a),
        relative_charge_residual=maximum(recorder.kcl_residual_a),
        relative_energy_residual=abs(integrated_residual) / scale,
        harmonic_metrics=Dict(
            :output_input_frequency_ratio =>
                study.specification.rated_bases.frequency_hz / study.input_frequency_hz,
            :failed_commutation_count => count(recorder.failed_commutation),
        ),
    )
    return ConverterSystems.SwitchingCycloconverterTrace(
        recorder.time_s, recorder.source_voltage_v, recorder.source_current_a,
        recorder.output_reference_voltage_v, recorder.output_voltage_v,
        recorder.output_current_a, recorder.firing_angle_rad,
        recorder.requested_bridge_group, recorder.active_bridge_group,
        recorder.requested_firing_state, recorder.applied_firing_state,
        recorder.conducting_state, recorder.commutation_overlap_valve_count,
        recorder.failed_commutation, recorder.reactive_power_var,
        recorder.stored_energy_j, recorder.semiconductor_loss_w,
        recorder.kcl_residual_a, recorder.energy_residual_w,
        recorder.junction_temperature_k, recorder.recovered_charge_c, result,
    )
end

function execute_switching_cycloconverter!(integrator::CycloconverterIntegrator)
    while _advance_cycloconverter!(integrator)
    end
    return _cycloconverter_result(integrator)
end

function _three_level_split_link_normalized_reference(study, time_s)
    modulation = study.specification.modulation
    angle = 2.0 * pi * study.specification.rated_bases.frequency_hz * time_s +
        modulation.phase_shift_rad
    return ntuple(phase ->
        modulation.modulation_index * sin(angle - (phase - 1) * 2.0 * pi / 3.0), 3)
end

function _npc_level_duties(study, time_s)
    reference = _three_level_split_link_normalized_reference(study, time_s)
    positive = ntuple(phase -> max(reference[phase], 0.0), 3)
    negative = ntuple(phase -> max(-reference[phase], 0.0), 3)
    neutral = ntuple(phase -> 1.0 - abs(reference[phase]), 3)
    return (positive=positive, neutral=neutral, negative=negative)
end

function _three_level_split_link_modulation_state(study, time_s)
    leg_topology = study.specification.selection.family ===
        ConverterSystems.ThreeLevelTTypeBridge ? :t_type : :neutral_point_clamped
    return ConverterSystems.converter_three_level_gate_state(
        _three_level_split_link_normalized_reference(study, time_s),
        time_s,
        study.specification.timing.carrier_frequency_hz,
        ; leg_topology,
    )
end

function _npc_average_system(study, time_s)
    duties = _npc_level_duties(study, time_s)
    resistance = study.load_resistance_ohm
    inductance = study.load_inductance_h
    capacitance = study.dc_link_capacitance_f
    source_resistance = study.source_resistance_ohm
    system = zeros(5, 5)
    forcing = zeros(5)
    mean_positive = sum(duties.positive) / 3.0
    mean_negative = sum(duties.negative) / 3.0
    for phase in 1:3
        system[phase, phase] = -resistance / inductance
        system[phase, 4] = (duties.positive[phase] - mean_positive) / inductance
        system[phase, 5] = (-duties.negative[phase] + mean_negative) / inductance
        system[4, phase] = -duties.positive[phase] / capacitance
        system[5, phase] = duties.negative[phase] / capacitance
    end
    link_coefficient = -inv(source_resistance * capacitance)
    system[4, 4] = link_coefficient
    system[4, 5] = link_coefficient
    system[5, 4] = link_coefficient
    system[5, 5] = link_coefficient
    forcing[4] = study.input_voltage_v / (source_resistance * capacitance)
    forcing[5] = forcing[4]
    return system, forcing, duties
end

function _npc_average_observables(study, state, time_s)
    current = view(state, 1:3)
    upper_voltage = state[4]
    lower_voltage = state[5]
    duties = _npc_level_duties(study, time_s)
    pole_voltage = [
        duties.positive[phase] * upper_voltage -
            duties.negative[phase] * lower_voltage
        for phase in 1:3
    ]
    neutral_voltage = sum(pole_voltage) / 3.0
    phase_voltage = pole_voltage .- neutral_voltage
    source_current = (study.input_voltage_v - upper_voltage - lower_voltage) /
        study.source_resistance_ohm
    midpoint_current = sum(duties.neutral[phase] * current[phase] for phase in 1:3)
    stored_energy = 0.5 * study.load_inductance_h * sum(abs2, current) +
        0.5 * study.dc_link_capacitance_f * (upper_voltage^2 + lower_voltage^2)
    return (
        source_current=source_current,
        upper_voltage=upper_voltage,
        lower_voltage=lower_voltage,
        midpoint_current=midpoint_current,
        phase_voltage=phase_voltage,
        phase_current=current,
        stored_energy=stored_energy,
    )
end

mutable struct AverageThreeLevelSplitLinkInverterRuntime{S}
    study::S
    state::Vector{Float64}
    time_s::Float64
    accepted_step_index::Int
    input_energy_j::Float64
    dissipated_energy_j::Float64
    time_trace_s::Vector{Float64}
    source_current_trace_a::Vector{Float64}
    upper_dc_link_voltage_trace_v::Vector{Float64}
    lower_dc_link_voltage_trace_v::Vector{Float64}
    midpoint_current_trace_a::Vector{Float64}
    phase_voltage_trace_v::Matrix{Float64}
    phase_current_trace_a::Matrix{Float64}
    stored_energy_trace_j::Vector{Float64}
    circuit_residual_trace_v::Matrix{Float64}
    energy_residual_trace_w::Vector{Float64}
end

function _record_average_npc!(runtime, sample, circuit_residual, energy_residual)
    observed = _npc_average_observables(runtime.study, runtime.state, runtime.time_s)
    runtime.time_trace_s[sample] = runtime.time_s
    runtime.source_current_trace_a[sample] = observed.source_current
    runtime.upper_dc_link_voltage_trace_v[sample] = observed.upper_voltage
    runtime.lower_dc_link_voltage_trace_v[sample] = observed.lower_voltage
    runtime.midpoint_current_trace_a[sample] = observed.midpoint_current
    runtime.phase_voltage_trace_v[:, sample] .= observed.phase_voltage
    runtime.phase_current_trace_a[:, sample] .= observed.phase_current
    runtime.stored_energy_trace_j[sample] = observed.stored_energy
    runtime.circuit_residual_trace_v[:, sample] .= circuit_residual
    runtime.energy_residual_trace_w[sample] = energy_residual
    return runtime
end

function prepare_average_three_level_neutral_point_clamped(
    study::ConverterSystems.AverageThreeLevelNeutralPointClampedStudy,
)
    step = study.specification.timing.fixed_step_s
    sample_count = round(Int, (study.stop_time_s - study.start_time_s) / step) + 1
    state = [
        study.initial_state.phase_current_a...,
        study.initial_state.upper_dc_link_voltage_v,
        study.initial_state.lower_dc_link_voltage_v,
    ]
    runtime = AverageThreeLevelSplitLinkInverterRuntime(
        study,
        state,
        study.start_time_s,
        0,
        0.0,
        0.0,
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, sample_count),
    )
    _record_average_npc!(runtime, 1, zeros(3), 0.0)
    return runtime
end

function _advance_average_npc!(runtime, sample)
    study = runtime.study
    step = study.specification.timing.fixed_step_s
    previous_time = runtime.time_s
    endpoint = study.start_time_s + (sample - 1) * step
    previous_state = copy(runtime.state)
    previous_observed = _npc_average_observables(study, previous_state, previous_time)
    previous_system, previous_forcing, _ = _npc_average_system(study, previous_time)
    endpoint_system, endpoint_forcing, _ = _npc_average_system(study, endpoint)
    identity_matrix = Matrix{Float64}(I, 5, 5)
    runtime.state .= (identity_matrix - 0.5 * step * endpoint_system) \
        ((identity_matrix + 0.5 * step * previous_system) * previous_state +
         0.5 * step * (previous_forcing + endpoint_forcing))
    endpoint_observed = _npc_average_observables(study, runtime.state, endpoint)
    average_current = 0.5 .* (previous_state[1:3] .+ runtime.state[1:3])
    average_voltage = 0.5 .* (
        previous_observed.phase_voltage .+ endpoint_observed.phase_voltage
    )
    derivative = (runtime.state[1:3] .- previous_state[1:3]) ./ step
    circuit_residual = average_voltage .-
        study.load_resistance_ohm .* average_current .-
        study.load_inductance_h .* derivative
    input_energy = 0.5 * step * study.input_voltage_v *
        (previous_observed.source_current + endpoint_observed.source_current)
    load_energy = 0.5 * step * study.load_resistance_ohm * (
        sum(abs2, previous_state[1:3]) + sum(abs2, runtime.state[1:3])
    )
    source_energy = 0.5 * step * study.source_resistance_ohm * (
        previous_observed.source_current^2 + endpoint_observed.source_current^2
    )
    energy_residual = (input_energy - load_energy - source_energy -
        (endpoint_observed.stored_energy - previous_observed.stored_energy)) / step
    runtime.input_energy_j += input_energy
    runtime.dissipated_energy_j += load_energy + source_energy
    runtime.time_s = endpoint
    return _record_average_npc!(runtime, sample, circuit_residual, energy_residual)
end

function _average_three_level_split_link_result(runtime)
    study = runtime.study
    observed = _npc_average_observables(study, runtime.state, runtime.time_s)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        repr(runtime.time_s),
        repr(runtime.state),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [observed.upper_voltage, observed.lower_voltage, observed.phase_voltage...],
        [observed.source_current, observed.midpoint_current, observed.phase_current...],
        BitVector(),
        BitVector(),
        BitVector(),
        study.dc_link_capacitance_f .* [observed.upper_voltage, observed.lower_voltage],
        study.load_inductance_h .* collect(observed.phase_current),
        Float64[],
        collect(_three_level_split_link_normalized_reference(study, runtime.time_s)),
        observed.stored_energy,
        runtime.dissipated_energy_j,
        length(runtime.time_trace_s) - 1,
        0,
        signature,
    )
    power_scale = max(maximum(abs,
        study.input_voltage_v .* runtime.source_current_trace_a), 1.0)
    result = ConverterSystems.converter_system_result(
        study.specification,
        state;
        accepted=true,
        status=:ok,
        maximum_kcl_residual_a=maximum(abs,
            vec(sum(runtime.phase_current_trace_a; dims=1))),
        relative_charge_residual=maximum(abs,
            runtime.upper_dc_link_voltage_trace_v .+
            runtime.lower_dc_link_voltage_trace_v .- study.input_voltage_v) /
            study.input_voltage_v,
        relative_energy_residual=abs(
            runtime.input_energy_j - runtime.dissipated_energy_j -
            (observed.stored_energy - runtime.stored_energy_trace_j[1]),
        ) / max(abs(runtime.input_energy_j), observed.stored_energy, eps(Float64)),
        harmonic_metrics=Dict(:phase_voltage_thd => 0.0),
    )
    return ConverterSystems.AverageThreeLevelNeutralPointClampedTrace(
        runtime.time_trace_s,
        runtime.source_current_trace_a,
        runtime.upper_dc_link_voltage_trace_v,
        runtime.lower_dc_link_voltage_trace_v,
        runtime.midpoint_current_trace_a,
        runtime.phase_voltage_trace_v,
        runtime.phase_current_trace_a,
        runtime.stored_energy_trace_j,
        runtime.circuit_residual_trace_v,
        runtime.energy_residual_trace_w,
        result,
    )
end

function execute_average_three_level_neutral_point_clamped!(
    runtime::AverageThreeLevelSplitLinkInverterRuntime,
)
    advance_prepared_converter_system!(
        runtime,
        length(runtime.time_trace_s) - 1 - runtime.accepted_step_index,
    )
    return _average_three_level_split_link_result(runtime)
end

mutable struct ThreeLevelSplitLinkInverterRuntime{S,N,B,C,A,D,L,P}
    study::S
    network::N
    bridges::B
    controlled_valves::C
    antiparallel_diodes::A
    clamp_diodes::D
    load_branches::L
    dc_link_capacitors::P
    dc_positive_node::Int
    midpoint_node::Int
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

mutable struct ThreeLevelSplitLinkInverterRecorder
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    upper_dc_link_voltage_v::Vector{Float64}
    lower_dc_link_voltage_v::Vector{Float64}
    midpoint_current_a::Vector{Float64}
    phase_voltage_v::Matrix{Float64}
    phase_current_a::Matrix{Float64}
    requested_level::Matrix{Int8}
    requested_gate_state::BitMatrix
    applied_gate_state::BitMatrix
    controlled_conducting_state::BitMatrix
    antiparallel_diode_conducting_state::BitMatrix
    clamp_diode_conducting_state::BitMatrix
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    controlled_junction_temperature_k::Matrix{Float64}
    antiparallel_junction_temperature_k::Matrix{Float64}
    clamp_junction_temperature_k::Matrix{Float64}
    antiparallel_recovered_charge_c::Matrix{Float64}
    clamp_recovered_charge_c::Matrix{Float64}
    controlled_turn_off_tail_current_a::Matrix{Float64}
    write_index::Int
end

mutable struct ThreeLevelSplitLinkInverterIntegrator{R,T}
    runtime::R
    transaction::T
    recorder::ThreeLevelSplitLinkInverterRecorder
    accepted_step_index::Int
    completed::Bool
    failed::Bool
    last_failure::Union{Nothing,String}
end

_three_level_split_link_is_detailed(study) =
    study.specification.selection.fidelity === StudyCore.SwitchingDetailed

_three_level_split_link_linear_network(runtime) =
    _three_level_split_link_is_detailed(runtime.study) ? nonlinear_linear_system(runtime.network) : runtime.network

function _three_level_split_link_nodes(study)
    leg_nodes = ntuple(3) do phase
        Dict(node.name => node.node for node in study.topologies[phase].nodes)
    end
    family = study.specification.selection.family
    required = family === ConverterSystems.ThreeLevelNeutralPointClampedBridge ? (
        :dc_positive,
        :midpoint,
        :ac_terminal,
        :dc_negative,
        :upper_clamp_node,
        :lower_clamp_node,
    ) : (
        :dc_positive,
        :midpoint,
        :ac_terminal,
        :dc_negative,
        :midpoint_path_node,
    )
    all(nodes -> all(haskey(nodes, name) for name in required), leg_nodes) ||
        throw(ArgumentError("three-level split-link topology is missing a canonical node"))
    all(phase -> (
        leg_nodes[phase][:dc_positive] == leg_nodes[1][:dc_positive] &&
        leg_nodes[phase][:midpoint] == leg_nodes[1][:midpoint] &&
        leg_nodes[phase][:dc_negative] == leg_nodes[1][:dc_negative]
    ), 2:3) || throw(ArgumentError("three-level legs must share exact DC rails and midpoint"))
    phase_nodes = ntuple(phase -> leg_nodes[phase][:ac_terminal], 3)
    length(unique(phase_nodes)) == 3 || throw(ArgumentError("three-level phase terminals must be distinct"))
    input = only(port for port in study.specification.ports if port.identity === :input_dc)
    output = only(port for port in study.specification.ports if port.identity === :output_ac)
    input.ordered_nodes == (leg_nodes[1][:dc_positive], leg_nodes[1][:dc_negative]) ||
        throw(ArgumentError("three-level DC port does not match its rails"))
    output.ordered_nodes == phase_nodes ||
        throw(ArgumentError("three-level AC port does not match its phase terminals"))
    maximum_node = maximum(node.node for topology in study.topologies for node in topology.nodes)
    return (
        dc_positive=leg_nodes[1][:dc_positive],
        midpoint=leg_nodes[1][:midpoint],
        dc_negative=leg_nodes[1][:dc_negative],
        phase=phase_nodes,
        neutral=maximum_node + 1,
        leg=leg_nodes,
    )
end

function _three_level_split_link_controlled_valve(study, position, commanded, closed)
    parameters = study.detailed_semiconductor
    detailed = _three_level_split_link_is_detailed(study)
    timing = study.specification.timing
    device = Nonlinear.IGBTSwitch(
        position.from_node,
        position.to_node;
        gate_driver=Nonlinear.PowerSemiconductorGateDriver(
            dead_time_s=timing.dead_time_s,
            minimum_pulse_width_s=timing.minimum_pulse_s,
            initially_on=closed,
        ),
        on_conductance=detailed ? parameters.controlled_on_conductance_s : 1.0e6,
        off_conductance=detailed ? parameters.controlled_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? parameters.controlled_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            parameters.snubber_resistance_ohm,
            parameters.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ? _detailed_chopper_controlled_fidelity(parameters) : nothing,
        initially_closed=closed,
    )
    commanded == closed || Nonlinear.request_power_semiconductor_gate!(
        device,
        commanded,
        study.start_time_s;
        earliest_transition_time_s=commanded ?
            study.start_time_s + timing.dead_time_s : study.start_time_s,
    )
    return device
end

function _three_level_split_link_diode(study, from_node, to_node, closed)
    parameters = study.detailed_semiconductor
    detailed = _three_level_split_link_is_detailed(study)
    return Nonlinear.DiodeValveSwitch(
        from_node,
        to_node;
        on_conductance=detailed ? parameters.freewheel_on_conductance_s : 1.0e6,
        off_conductance=detailed ? parameters.freewheel_off_conductance_s : 1.0e-9,
        forward_voltage_drop_v=detailed ? parameters.freewheel_forward_voltage_v : 0.0,
        snubber=detailed ? Nonlinear.SeriesRCSnubber(
            parameters.snubber_resistance_ohm,
            parameters.snubber_capacitance_f,
        ) : nothing,
        extended_fidelity=detailed ? _detailed_chopper_freewheel_fidelity(parameters) : nothing,
        initially_closed=closed,
    )
end

function _npc_initial_paths(levels, phase_current)
    controlled = falses(12)
    antiparallel = falses(12)
    clamp = falses(6)
    for phase in 1:3
        offset = 4 * (phase - 1)
        clamp_offset = 2 * (phase - 1)
        current = phase_current[phase]
        if levels[phase] == 1
            target = current >= 0.0 ? controlled : antiparallel
            target[(offset + 1):(offset + 2)] .= true
        elseif levels[phase] == 0
            if current >= 0.0
                controlled[offset + 2] = true
                clamp[clamp_offset + 1] = true
            else
                controlled[offset + 3] = true
                clamp[clamp_offset + 2] = true
            end
        else
            target = current < 0.0 ? controlled : antiparallel
            target[(offset + 3):(offset + 4)] .= true
        end
    end
    return controlled, antiparallel, clamp
end

function _t_type_initial_paths(levels, phase_current)
    controlled = falses(12)
    antiparallel = falses(12)
    for phase in 1:3
        offset = 4 * (phase - 1)
        current = phase_current[phase]
        if levels[phase] == 1
            (current >= 0.0 ? controlled : antiparallel)[offset + 1] = true
        elseif levels[phase] == 0
            if current >= 0.0
                antiparallel[offset + 2] = true
                controlled[offset + 3] = true
            else
                controlled[offset + 2] = true
                antiparallel[offset + 3] = true
            end
        else
            (current < 0.0 ? controlled : antiparallel)[offset + 4] = true
        end
    end
    return controlled, antiparallel, falses(0)
end

function _seed_three_level_split_link_inverter(study)
    nodes = _three_level_split_link_nodes(study)
    modulation = _three_level_split_link_modulation_state(study, study.start_time_s)
    is_npc = study.specification.selection.family ===
        ConverterSystems.ThreeLevelNeutralPointClampedBridge
    controlled_paths, antiparallel_paths, clamp_paths =
        (is_npc ? _npc_initial_paths : _t_type_initial_paths)(
            modulation.requested_level,
            study.initial_state.phase_current_a,
        )
    controlled = ntuple(12) do index
        phase = div(index - 1, 4) + 1
        position_index = mod(index - 1, 4) + 1
        position = study.topologies[phase].valve_positions[position_index]
        _three_level_split_link_controlled_valve(
            study,
            position,
            modulation.requested_valve_state[index],
            controlled_paths[index],
        )
    end
    antiparallel = ntuple(12) do index
        phase = div(index - 1, 4) + 1
        position_index = mod(index - 1, 4) + 1
        position = study.topologies[phase].valve_positions[position_index]
        _three_level_split_link_diode(study, position.to_node, position.from_node, antiparallel_paths[index])
    end
    clamps = is_npc ? ntuple(6) do index
            phase = div(index - 1, 2) + 1
            position_index = mod(index - 1, 2) + 5
            position = study.topologies[phase].valve_positions[position_index]
            _three_level_split_link_diode(
                study,
                position.from_node,
                position.to_node,
                clamp_paths[index],
            )
        end : ()
    # The T-type midpoint devices form one coupled back-to-back commutation path;
    # their individual floating-node residuals are not physical localization events.
    for index in eachindex(controlled)
        position = mod(index - 1, 4) + 1
        if is_npc || position in (1, 4)
            Nonlinear.power_semiconductor_event_localization!(controlled[index])
            Nonlinear.power_semiconductor_event_localization!(antiparallel[index])
        end
    end
    foreach(Nonlinear.power_semiconductor_event_localization!, clamps)
    bridges = ntuple(3) do phase
        first_controlled = 4 * (phase - 1) + 1
        devices = if is_npc
            first_clamp = 2 * (phase - 1) + 1
            Nonlinear.PowerSemiconductorSwitch[
                controlled[first_controlled:(first_controlled + 3)]...,
                clamps[first_clamp:(first_clamp + 1)]...,
            ]
        else
            Nonlinear.PowerSemiconductorSwitch[
                controlled[first_controlled:(first_controlled + 3)]...,
            ]
        end
        Nonlinear.PowerSemiconductorBridgeTopology(study.topologies[phase], devices)
    end
    upper_voltage = study.initial_state.upper_dc_link_voltage_v
    lower_voltage = study.initial_state.lower_dc_link_voltage_v
    source = Branches.TwoTerminalTheveninSource(
        nodes.dc_positive,
        nodes.dc_negative,
        inv(study.source_resistance_ohm),
        _time_s -> study.input_voltage_v,
    )
    capacitors = (
        Branches.CapacitorBranch(
            nodes.dc_positive,
            nodes.midpoint,
            study.dc_link_capacitance_f,
            0.0,
            upper_voltage,
            0.0,
        ),
        Branches.CapacitorBranch(
            nodes.midpoint,
            nodes.dc_negative,
            study.dc_link_capacitance_f,
            0.0,
            lower_voltage,
            0.0,
        ),
    )
    pole_voltage = ntuple(phase -> modulation.requested_level[phase] == 1 ?
        upper_voltage + lower_voltage :
        modulation.requested_level[phase] == 0 ? lower_voltage : 0.0, 3)
    neutral_voltage = sum(pole_voltage) / 3.0
    loads = ntuple(3) do phase
        phase_voltage = pole_voltage[phase] - neutral_voltage
        current = study.initial_state.phase_current_a[phase]
        Branches.SeriesRLBranch(
            nodes.phase[phase],
            nodes.neutral,
            study.load_resistance_ohm,
            study.load_inductance_h,
            current,
            phase_voltage,
            current,
        )
    end
    node_count = nodes.neutral
    elements = _three_level_split_link_is_detailed(study) ?
        Any[source, bridges..., capacitors..., loads...] :
        Any[source, bridges..., antiparallel..., capacitors..., loads...]
    linear_network = Nodal.NodalSystem(node_count, elements)
    linear_network.v[nodes.dc_positive] = upper_voltage + lower_voltage
    linear_network.v[nodes.midpoint] = lower_voltage
    nodes.dc_negative == 0 || (linear_network.v[nodes.dc_negative] = 0.0)
    for phase in 1:3
        linear_network.v[nodes.phase[phase]] = pole_voltage[phase]
        if is_npc
            linear_network.v[nodes.leg[phase][:upper_clamp_node]] =
                modulation.requested_level[phase] == 1 ?
                    upper_voltage + lower_voltage : lower_voltage
            linear_network.v[nodes.leg[phase][:lower_clamp_node]] =
                modulation.requested_level[phase] == -1 ? 0.0 : lower_voltage
        else
            linear_network.v[nodes.leg[phase][:midpoint_path_node]] = lower_voltage
        end
    end
    linear_network.v[nodes.neutral] = neutral_voltage
    devices = (controlled..., antiparallel..., clamps...)
    if _three_level_split_link_is_detailed(study)
        for valve in devices
            Nonlinear.initialize_power_semiconductor_junction_state!(
                valve,
                Branches.branch_voltage(linear_network.v, valve.a, valve.b),
            )
        end
    end
    network = _three_level_split_link_is_detailed(study) ? NonlinearNodalSystem(
        linear_network,
        devices;
        scales=NonlinearNetworkScales(
            node_count,
            0;
            nominal_voltage_v=study.input_voltage_v,
            nominal_current_a=max(maximum(abs, study.initial_state.phase_current_a), 1.0),
        ),
        options=NonlinearSolveOptions(
            maximum_iterations=50,
            current_absolute_tolerance_a=1.0e-9,
            voltage_absolute_tolerance_v=1.0e-9,
            scaled_step_tolerance=1.0e-12,
        ),
    ) : linear_network
    return ThreeLevelSplitLinkInverterRuntime(
        study,
        network,
        bridges,
        controlled,
        antiparallel,
        clamps,
        loads,
        capacitors,
        nodes.dc_positive,
        nodes.midpoint,
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

function ThreeLevelSplitLinkInverterRecorder(sample_count, clamp_count)
    sample_count > 1 || throw(ArgumentError(
        "three-level split-link recorder requires at least two samples",
    ))
    clamp_count in (0, 6) || throw(ArgumentError(
        "three-level split-link recorder received an invalid clamp-device count",
    ))
    return ThreeLevelSplitLinkInverterRecorder(
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Float64, 3, sample_count),
        zeros(Int8, 3, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        falses(12, sample_count),
        falses(clamp_count, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, clamp_count, sample_count),
        zeros(Float64, 12, sample_count),
        zeros(Float64, clamp_count, sample_count),
        zeros(Float64, 12, sample_count),
        0,
    )
end

_three_level_split_link_devices(runtime) = (
    runtime.controlled_valves...,
    runtime.antiparallel_diodes...,
    runtime.clamp_diodes...,
)

function _three_level_split_link_stored_energy(runtime)
    network = _three_level_split_link_linear_network(runtime)
    energy = 0.5 * runtime.study.load_inductance_h *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    for capacitor in runtime.dc_link_capacitors
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        energy += 0.5 * capacitor.c * voltage^2
    end
    _three_level_split_link_is_detailed(runtime.study) || return energy
    return energy + sum(_three_level_split_link_devices(runtime); init=0.0) do valve
        terminal = Nonlinear.power_semiconductor_terminal_state(valve)
        extended = Nonlinear.power_semiconductor_extended_state(valve)
        terminal.snubber_capacitor_energy_j + extended.junction_stored_energy_j
    end
end

function _three_level_split_link_power(runtime)
    study = runtime.study
    network = _three_level_split_link_linear_network(runtime)
    dc_voltage = network.v[runtime.dc_positive_node] -
        (runtime.dc_negative_node == 0 ? 0.0 : network.v[runtime.dc_negative_node])
    source_current = (study.input_voltage_v - dc_voltage) / study.source_resistance_ohm
    load_loss = study.load_resistance_ohm *
        sum(branch.i_last^2 for branch in runtime.load_branches)
    semiconductor_loss = sum(
        valve.last_semiconductor_loss_w +
            (valve.snubber === nothing ? 0.0 : valve.snubber.last_resistor_loss_w)
        for valve in _three_level_split_link_devices(runtime)
    )
    return (
        input_w=study.input_voltage_v * source_current,
        load_w=load_loss,
        source_loss_w=study.source_resistance_ohm * source_current^2,
        semiconductor_loss_w=semiconductor_loss,
        source_current_a=source_current,
        dc_voltage_v=dc_voltage,
    )
end

function _three_level_split_link_linear_companion(runtime)
    network = _three_level_split_link_linear_network(runtime)
    terminal_power = 0.0
    dissipated_power = 0.0
    stored_energy = 0.0
    for branch in runtime.load_branches
        voltage = Branches.branch_voltage(network.v, branch.a, branch.b)
        terminal_power += voltage * branch.i_last
        dissipated_power += runtime.study.load_resistance_ohm * branch.i_last^2
        stored_energy += 0.5 * runtime.study.load_inductance_h * branch.i_last^2
    end
    for capacitor in runtime.dc_link_capacitors
        voltage = Branches.branch_voltage(network.v, capacitor.a, capacitor.b)
        terminal_power += voltage * capacitor.i_last
        stored_energy += 0.5 * capacitor.c * voltage^2
    end
    return (
        terminal_power_w=terminal_power,
        dissipated_power_w=dissipated_power,
        stored_energy_j=stored_energy,
    )
end

function _three_level_split_link_device_companion_residual(runtime)
    _three_level_split_link_is_detailed(runtime.study) || return 0.0
    return sum(
        Nonlinear.power_semiconductor_extended_state(valve).companion_energy_residual_j
        for valve in _three_level_split_link_devices(runtime)
    )
end

function _observe_three_level_split_link_substep!(
    accounting::FullBridgeStepEnergyAccounting,
    runtime,
    _system,
    _time_s,
    step_s,
    companion_method,
    _substep_index,
)
    final_power = _three_level_split_link_power(runtime)
    final_linear = _three_level_split_link_linear_companion(runtime)
    previous_power = accounting.previous_power
    previous_linear = accounting.previous_linear_companion
    weight = companion_method === Branches.TrapezoidalCompanion ? 0.5 :
        companion_method === Branches.BackwardEulerCompanion ? 1.0 :
        throw(ArgumentError("three-level split-link inverter received an unknown companion method"))
    blend(previous, final) = weight == 0.5 ? 0.5 * (previous + final) : final
    accounting.input_energy_j += step_s * blend(previous_power.input_w, final_power.input_w)
    accounting.load_energy_j += step_s * blend(previous_power.load_w, final_power.load_w)
    accounting.source_dissipated_energy_j += step_s *
        blend(previous_power.source_loss_w, final_power.source_loss_w)
    accounting.semiconductor_dissipated_energy_j += step_s *
        blend(previous_power.semiconductor_loss_w, final_power.semiconductor_loss_w)
    accounting.linear_companion_energy_residual_j += step_s *
        blend(previous_linear.terminal_power_w, final_linear.terminal_power_w) -
        step_s * blend(previous_linear.dissipated_power_w, final_linear.dissipated_power_w) -
        (final_linear.stored_energy_j - previous_linear.stored_energy_j)
    accounting.previous_power = final_power
    accounting.previous_linear_companion = final_linear
    return nothing
end

function _request_three_level_split_link_gates!(runtime, time_s)
    modulation = _three_level_split_link_modulation_state(runtime.study, time_s)
    changed = false
    timing = runtime.study.specification.timing
    for phase in 1:3
        first = 4 * (phase - 1) + 1
        requested = modulation.requested_valve_state[first:(first + 3)]
        is_npc = runtime.study.specification.selection.family ===
            ConverterSystems.ThreeLevelNeutralPointClampedBridge
        full_state = is_npc ? BitVector((requested..., false, false)) :
            BitVector(requested)
        BridgeTopologies.bridge_topology_state_is_allowed(
            runtime.study.topologies[phase],
            full_state,
        ) || throw(ArgumentError("three-level modulation requested a prohibited leg state"))
        transition_count = 0
        for local_index in 1:4
            valve = runtime.controlled_valves[first + local_index - 1]
            desired = requested[local_index]
            desired == valve.gate_driver.commanded_on && continue
            Nonlinear.request_power_semiconductor_gate!(
                valve,
                desired,
                time_s;
                earliest_transition_time_s=desired ? time_s + timing.dead_time_s : time_s,
            )
            transition_count += 1
        end
        runtime.bridges[phase].transition_count += transition_count
        changed |= transition_count > 0
    end
    changed && push!(runtime.events, ConverterSystems.ConverterSystemEventRecord(
        time_s,
        :three_level_pwm_command,
        runtime.study.specification.selection.family ===
            ConverterSystems.ThreeLevelNeutralPointClampedBridge ?
                :three_level_neutral_point_clamped_bridge : :three_level_t_type_bridge,
        true;
        message="accepted three-level carrier-PWM gate command",
    ))
    return changed
end

function _synchronize_three_level_split_link_conduction!(runtime, time_s)
    network = _three_level_split_link_linear_network(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    currents = ntuple(phase -> runtime.load_branches[phase].i_last, 3)
    modulation = _three_level_split_link_modulation_state(runtime.study, time_s)
    is_npc = runtime.study.specification.selection.family ===
        ConverterSystems.ThreeLevelNeutralPointClampedBridge
    controlled_target = falses(12)
    antiparallel_target = falses(12)
    clamp_target = falses(length(runtime.clamp_diodes))
    for phase in 1:3
        offset = 4 * (phase - 1)
        applied = ntuple(index ->
            runtime.controlled_valves[offset + index].gate_driver.applied_on, 4)
        level = modulation.requested_level[phase]
        current = currents[phase]
        if is_npc
            clamp_offset = 2 * (phase - 1)
            if current >= 0.0
                if level == 1 && applied[1] && applied[2]
                    controlled_target[(offset + 1):(offset + 2)] .= true
                elseif level >= 0 && applied[2]
                    controlled_target[offset + 2] = true
                    clamp_target[clamp_offset + 1] = true
                else
                    antiparallel_target[(offset + 3):(offset + 4)] .= true
                end
            else
                if level == -1 && applied[3] && applied[4]
                    controlled_target[(offset + 3):(offset + 4)] .= true
                elseif level <= 0 && applied[3]
                    controlled_target[offset + 3] = true
                    clamp_target[clamp_offset + 2] = true
                else
                    antiparallel_target[(offset + 1):(offset + 2)] .= true
                end
            end
        elseif current >= 0.0
            if level == 1 && applied[1]
                controlled_target[offset + 1] = true
            elseif level >= 0 && applied[2] && applied[3]
                antiparallel_target[offset + 2] = true
                controlled_target[offset + 3] = true
            else
                antiparallel_target[offset + 4] = true
            end
        else
            if level == -1 && applied[4]
                controlled_target[offset + 4] = true
            elseif level <= 0 && applied[2] && applied[3]
                controlled_target[offset + 2] = true
                antiparallel_target[offset + 3] = true
            else
                antiparallel_target[offset + 1] = true
            end
        end
    end
    changed = false
    for index in 1:12
        changed |= _set_full_bridge_device_state!(
            runtime.controlled_valves[index],
            controlled_target[index] &&
                runtime.controlled_valves[index].gate_driver.applied_on,
            time_s,
        )
        changed |= _set_full_bridge_device_state!(
            runtime.antiparallel_diodes[index],
            antiparallel_target[index],
            time_s,
        )
    end
    for index in eachindex(runtime.clamp_diodes)
        changed |= _set_full_bridge_device_state!(
            runtime.clamp_diodes[index],
            clamp_target[index],
            time_s,
        )
    end
    isfinite(neutral_voltage) || throw(ArgumentError(
        "three-level split-link neutral voltage became nonfinite",
    ))
    return changed
end

function _stabilize_three_level_split_link_topology!(runtime, time_s, step_s)
    _three_level_split_link_is_detailed(runtime.study) && throw(ArgumentError(
        "switching-detailed three-level split-link execution must use the nonlinear D200 solver",
    ))
    network = _three_level_split_link_linear_network(runtime)
    foreach(bridge -> Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
        bridge,
        time_s,
    ), runtime.bridges)
    for _iteration in 1:12
        changed = _synchronize_three_level_split_link_conduction!(runtime, time_s)
        Nodal.solve_algebraic_state!(network, time_s, step_s)
        changed || return nothing
    end
    throw(ArgumentError("three-level split-link switching topology failed to stabilize"))
end

function _record_three_level_split_link!(recorder, runtime, sample, energy_residual)
    network = _three_level_split_link_linear_network(runtime)
    power = _three_level_split_link_power(runtime)
    modulation = _three_level_split_link_modulation_state(runtime.study, runtime.time_s)
    neutral_voltage = network.v[runtime.load_neutral_node]
    upper_voltage = Branches.branch_voltage(
        network.v,
        runtime.dc_link_capacitors[1].a,
        runtime.dc_link_capacitors[1].b,
    )
    lower_voltage = Branches.branch_voltage(
        network.v,
        runtime.dc_link_capacitors[2].a,
        runtime.dc_link_capacitors[2].b,
    )
    recorder.time_s[sample] = runtime.time_s
    recorder.source_current_a[sample] = power.source_current_a
    recorder.upper_dc_link_voltage_v[sample] = upper_voltage
    recorder.lower_dc_link_voltage_v[sample] = lower_voltage
    recorder.midpoint_current_a[sample] =
        runtime.dc_link_capacitors[1].i_last - runtime.dc_link_capacitors[2].i_last
    recorder.requested_level[:, sample] .= modulation.requested_level
    for phase in 1:3
        recorder.phase_voltage_v[phase, sample] =
            network.v[runtime.phase_nodes[phase]] - neutral_voltage
        recorder.phase_current_a[phase, sample] = runtime.load_branches[phase].i_last
    end
    for index in 1:12
        controlled = runtime.controlled_valves[index]
        antiparallel = runtime.antiparallel_diodes[index]
        recorder.requested_gate_state[index, sample] = controlled.gate_driver.commanded_on
        recorder.applied_gate_state[index, sample] = controlled.gate_driver.applied_on
        recorder.controlled_conducting_state[index, sample] = controlled.closed
        recorder.antiparallel_diode_conducting_state[index, sample] = antiparallel.closed
        if _three_level_split_link_is_detailed(runtime.study)
            controlled_state = Nonlinear.power_semiconductor_extended_state(controlled)
            antiparallel_state = Nonlinear.power_semiconductor_extended_state(antiparallel)
            recorder.controlled_junction_temperature_k[index, sample] =
                controlled_state.junction_temperature_k
            recorder.antiparallel_junction_temperature_k[index, sample] =
                antiparallel_state.junction_temperature_k
            recorder.antiparallel_recovered_charge_c[index, sample] =
                antiparallel_state.stored_recovery_charge_c
            recorder.controlled_turn_off_tail_current_a[index, sample] =
                controlled_state.tail_current_a
        else
            recorder.controlled_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_junction_temperature_k[index, sample] = 0.0
            recorder.antiparallel_recovered_charge_c[index, sample] = 0.0
            recorder.controlled_turn_off_tail_current_a[index, sample] = 0.0
        end
    end
    for index in eachindex(runtime.clamp_diodes)
        clamp = runtime.clamp_diodes[index]
        recorder.clamp_diode_conducting_state[index, sample] = clamp.closed
        if _three_level_split_link_is_detailed(runtime.study)
            clamp_state = Nonlinear.power_semiconductor_extended_state(clamp)
            recorder.clamp_junction_temperature_k[index, sample] =
                clamp_state.junction_temperature_k
            recorder.clamp_recovered_charge_c[index, sample] =
                clamp_state.stored_recovery_charge_c
        else
            recorder.clamp_junction_temperature_k[index, sample] = 0.0
            recorder.clamp_recovered_charge_c[index, sample] = 0.0
        end
    end
    recorder.stored_energy_j[sample] = _three_level_split_link_stored_energy(runtime)
    recorder.semiconductor_loss_w[sample] = power.semiconductor_loss_w
    recorder.kcl_residual_a[sample] = _three_level_split_link_is_detailed(runtime.study) ?
        runtime.last_kcl_residual_a :
        maximum(abs, network.y * network.v - network.rhs; init=0.0)
    recorder.energy_residual_w[sample] = energy_residual
    recorder.write_index = sample
    return recorder
end

function prepare_switching_three_level_neutral_point_clamped(
    study::ConverterSystems.SwitchingThreeLevelNeutralPointClampedStudy,
)
    runtime = _seed_three_level_split_link_inverter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    ) + 1
    integrator = ThreeLevelSplitLinkInverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        ThreeLevelSplitLinkInverterRecorder(sample_count, 6),
        0,
        false,
        false,
        nothing,
    )
    _record_three_level_split_link!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end


function prepare_switching_three_level_t_type(
    study::ConverterSystems.SwitchingThreeLevelTTypeStudy,
)
    runtime = _seed_three_level_split_link_inverter(study)
    sample_count = round(Int,
        (study.stop_time_s - study.start_time_s) / study.specification.timing.fixed_step_s,
    ) + 1
    integrator = ThreeLevelSplitLinkInverterIntegrator(
        runtime,
        TimestepTransaction(runtime),
        ThreeLevelSplitLinkInverterRecorder(sample_count, 0),
        0,
        false,
        false,
        nothing,
    )
    _record_three_level_split_link!(integrator.recorder, runtime, 1, 0.0)
    return integrator
end

function _advance_three_level_split_link_inverter!(integrator)
    integrator.failed && throw(ArgumentError(
        "NPC integrator is terminally failed: $(integrator.last_failure)",
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
    previous_power = _three_level_split_link_power(runtime)
    previous_linear = _three_level_split_link_linear_companion(runtime)
    previous_stored = _three_level_split_link_stored_energy(runtime)
    previous_device_residual = _three_level_split_link_device_companion_residual(runtime)
    previous_device_dissipation = sum(
        valve.semiconductor_dissipated_energy_j +
            (valve.snubber === nothing ? 0.0 : valve.snubber.dissipated_energy_j)
        for valve in _three_level_split_link_devices(runtime)
    )
    begin_timestep_transaction!(integrator.transaction)
    try
        gate_changed = _request_three_level_split_link_gates!(runtime, endpoint)
        accounting = FullBridgeStepEnergyAccounting(previous_power, previous_linear)
        device_increment = 0.0
        if _three_level_split_link_is_detailed(study)
            foreach(bridge -> Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
                bridge,
                endpoint,
            ), runtime.bridges)
            conduction_changed = _synchronize_three_level_split_link_conduction!(runtime, endpoint)
            nonlinear_result = advance_nonlinear_step!(
                runtime.network,
                endpoint,
                step;
                discontinuity_treatment=(gate_changed || conduction_changed) ?
                    :two_backward_euler_half_steps : :none,
                discontinuity_reason=(gate_changed || conduction_changed) ?
                    :topology_change : :none,
                accepted_substep_observer=(system, time_s, substep_s,
                    companion_method, substep_index) -> _observe_three_level_split_link_substep!(
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
                for valve in _three_level_split_link_devices(runtime)
            )
            device_increment = final_device_dissipation - previous_device_dissipation
        else
            event_count = length(runtime.events)
            _stabilize_three_level_split_link_topology!(runtime, endpoint, step)
            topology_changed = gate_changed || length(runtime.events) != event_count
            network = _three_level_split_link_linear_network(runtime)
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
                _observe_three_level_split_link_substep!(
                    accounting,
                    runtime,
                    network,
                    endpoint - half_step,
                    half_step,
                    Branches.BackwardEulerCompanion,
                    1,
                )
                _stabilize_three_level_split_link_topology!(runtime, endpoint, half_step)
                Nodal.solve_algebraic_state!(network, endpoint, half_step, backward_euler)
                Nodal.accept_algebraic_state!(network, half_step, backward_euler)
                _observe_three_level_split_link_substep!(
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
                _observe_three_level_split_link_substep!(
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
        final_stored = _three_level_split_link_stored_energy(runtime)
        device_residual_increment =
            _three_level_split_link_device_companion_residual(runtime) - previous_device_residual
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
        _record_three_level_split_link!(integrator.recorder, runtime, next_step + 1, energy_residual)
        return true
    catch error
        timestep_transaction_active(integrator.transaction) &&
            restore_timestep_transaction!(integrator.transaction)
        integrator.failed = true
        integrator.last_failure = sprint(showerror, error)
        rethrow(error)
    end
end

function _three_level_split_link_result(integrator)
    runtime = integrator.runtime
    study = runtime.study
    recorder = integrator.recorder
    network = _three_level_split_link_linear_network(runtime)
    power = _three_level_split_link_power(runtime)
    neutral_voltage = network.v[runtime.load_neutral_node]
    phase_voltage = [network.v[node] - neutral_voltage for node in runtime.phase_nodes]
    phase_current = [branch.i_last for branch in runtime.load_branches]
    upper_voltage = Branches.branch_voltage(
        network.v,
        runtime.dc_link_capacitors[1].a,
        runtime.dc_link_capacitors[1].b,
    )
    lower_voltage = Branches.branch_voltage(
        network.v,
        runtime.dc_link_capacitors[2].a,
        runtime.dc_link_capacitors[2].b,
    )
    bridge_states = ntuple(phase -> Nonlinear.power_semiconductor_bridge_topology_state(
        runtime.bridges[phase],
        network.v,
        study.specification.timing.fixed_step_s,
    ), 3)
    stored_energy = _three_level_split_link_stored_energy(runtime)
    signature = bytes2hex(sha256(join((
        study.specification.signature_sha256,
        (state.deterministic_signature for state in bridge_states)...,
        repr(runtime.time_s),
        repr(phase_voltage),
        repr(phase_current),
        repr(upper_voltage),
        repr(lower_voltage),
        string(integrator.accepted_step_index),
    ), '\n')))
    state = ConverterSystems.ConverterSystemState(
        runtime.time_s,
        [upper_voltage, lower_voltage, phase_voltage...],
        [power.source_current_a, recorder.midpoint_current_a[end], phase_current...],
        BitVector(getfield.(runtime.controlled_valves, :gate_driver) .|> driver -> driver.commanded_on),
        BitVector(getfield.(runtime.controlled_valves, :gate_driver) .|> driver -> driver.applied_on),
        BitVector((
            getfield.(runtime.controlled_valves, :closed)...,
            getfield.(runtime.antiparallel_diodes, :closed)...,
            getfield.(runtime.clamp_diodes, :closed)...,
        )),
        study.dc_link_capacitance_f .* [upper_voltage, lower_voltage],
        study.load_inductance_h .* phase_current,
        Float64[],
        Float64.(recorder.requested_level[:, end]),
        stored_energy,
        runtime.dissipated_energy_j,
        integrator.accepted_step_index,
        length(runtime.events),
        signature,
    )
    energy_scale = max(
        abs(runtime.input_energy_j),
        abs(runtime.dissipated_energy_j),
        abs(stored_energy),
        eps(Float64),
    )
    integrated_residual = runtime.input_energy_j - runtime.dissipated_energy_j -
        (stored_energy - recorder.stored_energy_j[1]) -
        runtime.linear_companion_energy_residual_j - _three_level_split_link_device_companion_residual(runtime)
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
    common = (
        recorder.time_s,
        recorder.source_current_a,
        recorder.upper_dc_link_voltage_v,
        recorder.lower_dc_link_voltage_v,
        recorder.midpoint_current_a,
        recorder.phase_voltage_v,
        recorder.phase_current_a,
        recorder.requested_level,
        recorder.requested_gate_state,
        recorder.applied_gate_state,
        recorder.controlled_conducting_state,
        recorder.antiparallel_diode_conducting_state,
    )
    tail = (
        recorder.stored_energy_j,
        recorder.semiconductor_loss_w,
        recorder.kcl_residual_a,
        recorder.energy_residual_w,
        recorder.controlled_junction_temperature_k,
        recorder.antiparallel_junction_temperature_k,
    )
    detailed_tail = (
        recorder.antiparallel_recovered_charge_c,
        recorder.controlled_turn_off_tail_current_a,
        result,
    )
    if study.specification.selection.family ===
       ConverterSystems.ThreeLevelNeutralPointClampedBridge
        return ConverterSystems.SwitchingThreeLevelNeutralPointClampedTrace(
            common...,
            recorder.clamp_diode_conducting_state,
            tail...,
            recorder.clamp_junction_temperature_k,
            recorder.antiparallel_recovered_charge_c,
            recorder.clamp_recovered_charge_c,
            recorder.controlled_turn_off_tail_current_a,
            result,
        )
    end
    return ConverterSystems.SwitchingThreeLevelTTypeTrace(
        common...,
        tail...,
        detailed_tail...,
    )
end

function execute_switching_three_level_split_link_inverter!(
    integrator::ThreeLevelSplitLinkInverterIntegrator,
)
    while _advance_three_level_split_link_inverter!(integrator)
    end
    return _three_level_split_link_result(integrator)
end

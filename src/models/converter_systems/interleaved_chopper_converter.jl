export InterleavedChopperInitialState,
       SwitchingInterleavedChopperStudy,
       SwitchingInterleavedChopperTrace,
       interleaved_chopper_semiconductor_signatures

struct InterleavedChopperInitialState
    channel_inductor_current_a::Vector{Float64}
    output_voltage_v::Float64

    function InterleavedChopperInitialState(channel_inductor_current_a, output_voltage_v::Real)
        currents = Float64.(channel_inductor_current_a)
        voltage = Float64(output_voltage_v)
        2 <= length(currents) <= 8 && all(isfinite, currents) && all(>=(0.0), currents) ||
            throw(ArgumentError(
                "interleaved chopper requires two through eight finite nonnegative initial channel currents",
            ))
        isfinite(voltage) && voltage >= 0.0 || throw(ArgumentError(
            "interleaved chopper initial output voltage must be finite and nonnegative",
        ))
        return new(currents, voltage)
    end
end

function interleaved_chopper_semiconductor_signatures(parameters)
    devices = DetailedChopperSemiconductorParameters[parameters...]
    2 <= length(devices) <= 8 || throw(ArgumentError(
        "interleaved chopper requires two through eight channel semiconductor owners",
    ))
    return Tuple(Iterators.flatten(detailed_chopper_semiconductor_signatures(device)
        for device in devices))
end

struct SwitchingInterleavedChopperStudy
    specification::ConverterSystemSpecification
    topologies::Vector{BridgeTopologyDescriptor}
    input_voltage_v::Float64
    source_resistance_ohm::Float64
    inductor_resistance_ohm::Vector{Float64}
    inductance_h::Vector{Float64}
    capacitance_f::Float64
    load_resistance_ohm::Float64
    initial_state::InterleavedChopperInitialState
    detailed_semiconductor::Union{Nothing,Vector{DetailedChopperSemiconductorParameters}}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingInterleavedChopperStudy(
        specification::ConverterSystemSpecification;
        topologies,
        input_voltage_v::Real,
        source_resistance_ohm::Real,
        inductor_resistance_ohm,
        inductance_h,
        capacitance_f::Real,
        load_resistance_ohm::Real,
        initial_state::InterleavedChopperInitialState,
        detailed_semiconductor=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === InterleavedChopper || throw(ArgumentError(
            "interleaved chopper study requires the canonical interleaved family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "interleaved chopper execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "interleaved chopper study is a standalone conversion owner",
        ))
        topology_owners = BridgeTopologyDescriptor[topologies...]
        count = selection.channel_count
        length(topology_owners) == count || throw(DimensionMismatch(
            "interleaved chopper requires one B200 topology per selected channel",
        ))
        all(topology -> topology.family === :step_down_chopper, topology_owners) ||
            throw(ArgumentError(
                "interleaved chopper channels require canonical step-down B200 topologies",
            ))
        specification.topology_signatures ==
            Tuple(bridge_topology_signature.(topology_owners)) || throw(ArgumentError(
                "interleaved chopper specification does not bind every exact channel topology",
            ))
        semiconductor_owners = if selection.fidelity === SwitchingDetailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed interleaved execution requires channel semiconductor owners",
            ))
            owners = DetailedChopperSemiconductorParameters[detailed_semiconductor...]
            length(owners) == count || throw(DimensionMismatch(
                "interleaved chopper requires one detailed semiconductor owner per channel",
            ))
            specification.device_fidelity_signatures ==
                interleaved_chopper_semiconductor_signatures(owners) ||
                throw(ArgumentError(
                    "interleaved chopper specification does not bind every exact D200 channel device",
                ))
            owners
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state interleaved execution cannot accept detailed semiconductor owners",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state interleaved execution cannot claim detailed device signatures",
            ))
            nothing
        end
        specification.modulation.kind === PhaseShiftedCarrierPulseWidthModulation ||
            throw(ArgumentError(
                "interleaved chopper requires phase-shifted carrier modulation",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("interleaved chopper specification is not ready"))
        resistances = Float64.(inductor_resistance_ohm)
        inductances = Float64.(inductance_h)
        length(resistances) == count && length(inductances) == count ||
            throw(DimensionMismatch(
                "interleaved chopper requires one resistance and inductance per channel",
            ))
        length(initial_state.channel_inductor_current_a) == count ||
            throw(DimensionMismatch(
                "interleaved chopper initial state does not match its channel count",
            ))
        scalar_values = Float64.((
            input_voltage_v,
            source_resistance_ohm,
            capacitance_f,
            load_resistance_ohm,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, (scalar_values..., resistances..., inductances...)) ||
            throw(ArgumentError("interleaved chopper parameters must be finite"))
        all(>(0.0), scalar_values[1:4]) || throw(ArgumentError(
            "interleaved chopper voltage, source resistance, capacitance, and load resistance must be positive",
        ))
        all(>=(0.0), resistances) && all(>(0.0), inductances) || throw(ArgumentError(
            "interleaved chopper channel resistances must be nonnegative and inductances positive",
        ))
        scalar_values[5] >= 0.0 && scalar_values[6] > scalar_values[5] ||
            throw(ArgumentError(
                "interleaved chopper stop time must follow its nonnegative start time",
            ))
        input_port = only(port for port in specification.ports if port.identity === :input_dc)
        output_port = only(port for port in specification.ports if port.identity === :output_dc)
        length(input_port.ordered_nodes) == 2 && length(output_port.ordered_nodes) == 2 ||
            throw(ArgumentError("interleaved chopper ports must each have two terminals"))
        channel_nodes = [Dict(node.name => node.node for node in topology.nodes)
            for topology in topology_owners]
        all(nodes -> all(haskey(nodes, name) for name in
            (:dc_positive, :output, :dc_negative)), channel_nodes) ||
            throw(ArgumentError("interleaved channel topology is missing a canonical terminal"))
        all(nodes -> input_port.ordered_nodes ==
            (nodes[:dc_positive], nodes[:dc_negative]), channel_nodes) ||
            throw(ArgumentError("interleaved channel input rails do not match the common input port"))
        all(nodes -> output_port.ordered_nodes[2] == nodes[:dc_negative], channel_nodes) ||
            throw(ArgumentError("interleaved channel and output negative rails must match"))
        switching_nodes = [nodes[:output] for nodes in channel_nodes]
        length(unique(switching_nodes)) == count &&
            output_port.ordered_nodes[1] ∉ switching_nodes || throw(ArgumentError(
                "interleaved channels require distinct switching nodes ahead of the shared output",
            ))
        timing = specification.timing
        timing.dead_time_s == 0.0 && timing.minimum_pulse_s == 0.0 ||
            throw(ArgumentError(
                "interleaved single-ended channels require zero bridge dead time and minimum pulse",
            ))
        if selection.fidelity === SwitchingDetailed
            timing.fixed_step_s <= minimum(min(
                device.recovered_charge_lifetime_s,
                device.turn_off_tail_time_s,
            ) for device in semiconductor_owners) / 10.0 || throw(ArgumentError(
                "interleaved timestep must resolve every recovered-charge and tail-current state",
            ))
        end
        carrier_steps = inv(timing.carrier_frequency_hz * timing.fixed_step_s)
        isapprox(carrier_steps, round(carrier_steps); atol=1.0e-10, rtol=1.0e-10) &&
            isapprox(carrier_steps / count, round(carrier_steps / count);
                atol=1.0e-10, rtol=1.0e-10) || throw(ArgumentError(
                    "interleaved carrier period and channel phase offsets must lie on the fixed-step calendar",
                ))
        for phase_step in (0:(count - 1)) .* (carrier_steps / count)
            for edge_step in (
                phase_step + carrier_steps * specification.modulation.duty / 2.0,
                phase_step + carrier_steps *
                    (1.0 - specification.modulation.duty / 2.0),
            )
                isapprox(edge_step, round(edge_step); atol=1.0e-10, rtol=1.0e-10) ||
                    throw(ArgumentError(
                        "interleaved PWM edges must lie on the fixed-step calendar",
                    ))
            end
        end
        horizon_steps = (scalar_values[6] - scalar_values[5]) / timing.fixed_step_s
        isapprox(horizon_steps, round(horizon_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "interleaved chopper horizon must contain integer fixed steps",
            ))
        if selection.fidelity === SwitchingDetailed
            all(device ->
                device.switching_energy_current_domain_maximum_a >=
                    specification.rated_bases.current_a / count &&
                device.switching_energy_blocking_voltage_domain_maximum_v >= scalar_values[1],
                semiconductor_owners) || throw(ArgumentError(
                "interleaved D200 switching-energy domains must include per-channel rated current and input voltage",
            ))
        end
        return new(
            specification,
            topology_owners,
            scalar_values[1],
            scalar_values[2],
            resistances,
            inductances,
            scalar_values[3],
            scalar_values[4],
            initial_state,
            semiconductor_owners,
            scalar_values[5],
            scalar_values[6],
        )
    end
end

struct SwitchingInterleavedChopperTrace
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
    result::ConverterSystemResult
end

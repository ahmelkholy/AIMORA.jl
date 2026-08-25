export DualActiveBridgeTransformerLink,
       DualActiveBridgeInitialState,
       SwitchingDualActiveBridgeStudy,
       SwitchingDualActiveBridgeTrace,
       dual_active_bridge_semiconductor_signatures

struct DualActiveBridgeTransformerLink
    turns_ratio::Float64
    leakage_resistance_ohm::Float64
    leakage_inductance_h::Float64
    provenance::ParameterProvenance
    deterministic_signature_sha256::String

    function DualActiveBridgeTransformerLink(
        turns_ratio::Real,
        leakage_resistance_ohm::Real,
        leakage_inductance_h::Real;
        provenance::ParameterProvenance,
    )
        ratio, resistance, inductance = Float64.((
            turns_ratio,
            leakage_resistance_ohm,
            leakage_inductance_h,
        ))
        isfinite(ratio) && ratio > 0.0 && isfinite(resistance) && resistance >= 0.0 &&
            isfinite(inductance) && inductance > 0.0 || throw(ArgumentError(
                "DAB transformer ratio and leakage inductance must be positive and leakage resistance nonnegative",
            ))
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "DAB transformer-link provenance must describe physical model parameters",
        ))
        signature = bytes2hex(sha256(join((
            "aimora_dab_transformer_link_v1",
            bitstring(ratio),
            bitstring(resistance),
            bitstring(inductance),
            repr(provenance),
        ), '\n')))
        return new(ratio, resistance, inductance, provenance, signature)
    end
end

struct DualActiveBridgeInitialState
    leakage_current_a::Float64

    function DualActiveBridgeInitialState(leakage_current_a::Real=0.0)
        current = Float64(leakage_current_a)
        isfinite(current) || throw(ArgumentError(
            "DAB initial leakage current must be finite",
        ))
        return new(current)
    end
end

function dual_active_bridge_semiconductor_signatures(parameters)
    devices = DetailedChopperSemiconductorParameters[parameters...]
    length(devices) == 2 || throw(DimensionMismatch(
        "DAB requires one detailed semiconductor owner per active bridge",
    ))
    return Tuple(Iterators.flatten(detailed_chopper_semiconductor_signatures(device)
        for device in devices))
end

struct SwitchingDualActiveBridgeStudy
    specification::ConverterSystemSpecification
    primary_topology::BridgeTopologyDescriptor
    secondary_topology::BridgeTopologyDescriptor
    primary_dc_voltage_v::Float64
    secondary_dc_voltage_v::Float64
    primary_source_resistance_ohm::Float64
    secondary_source_resistance_ohm::Float64
    transformer_link::DualActiveBridgeTransformerLink
    initial_state::DualActiveBridgeInitialState
    detailed_semiconductor::Union{Nothing,Vector{DetailedChopperSemiconductorParameters}}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingDualActiveBridgeStudy(
        specification::ConverterSystemSpecification;
        primary_topology::BridgeTopologyDescriptor,
        secondary_topology::BridgeTopologyDescriptor,
        primary_dc_voltage_v::Real,
        secondary_dc_voltage_v::Real,
        primary_source_resistance_ohm::Real,
        secondary_source_resistance_ohm::Real,
        transformer_link::DualActiveBridgeTransformerLink,
        initial_state::DualActiveBridgeInitialState,
        detailed_semiconductor=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === DualActiveBridge || throw(ArgumentError(
            "DAB study requires the canonical dual-active-bridge family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "DAB execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "DAB study is a standalone isolated conversion owner",
        ))
        primary_topology.family === :full_bridge &&
            secondary_topology.family === :full_bridge || throw(ArgumentError(
                "DAB requires two canonical B200 full-bridge topologies",
            ))
        specification.topology_signatures == (
            bridge_topology_signature(primary_topology),
            bridge_topology_signature(secondary_topology),
        ) || throw(ArgumentError(
            "DAB specification does not bind its exact primary and secondary topologies",
        ))
        semiconductor = if selection.fidelity === SwitchingDetailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed DAB execution requires bridge semiconductor owners",
            ))
            owners = DetailedChopperSemiconductorParameters[detailed_semiconductor...]
            specification.device_fidelity_signatures ==
                dual_active_bridge_semiconductor_signatures(owners) ||
                throw(ArgumentError(
                    "DAB specification does not bind both exact D200 bridge device owners",
                ))
            owners
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state DAB execution cannot accept detailed semiconductor owners",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state DAB execution cannot claim detailed device signatures",
            ))
            nothing
        end
        specification.passive_and_transformer_signatures ==
            (transformer_link.deterministic_signature_sha256,) || throw(ArgumentError(
                "DAB specification does not bind its exact ideal-ratio and leakage transformer link",
            ))
        specification.modulation.kind in (
            SinglePhaseShiftModulation,
            DualPhaseShiftModulation,
            TriplePhaseShiftModulation,
        ) || throw(ArgumentError(
            "DAB study requires a declared single-, dual-, or triple-phase-shift modulation",
        ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("DAB specification is not ready"))
        values = Float64.((
            primary_dc_voltage_v,
            secondary_dc_voltage_v,
            primary_source_resistance_ohm,
            secondary_source_resistance_ohm,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:4]) || throw(ArgumentError(
            "DAB DC voltages and source resistances must be finite and positive",
        ))
        values[5] >= 0.0 && values[6] > values[5] || throw(ArgumentError(
            "DAB stop time must follow its nonnegative start time",
        ))
        primary_nodes = Dict(node.name => node.node for node in primary_topology.nodes)
        secondary_nodes = Dict(node.name => node.node for node in secondary_topology.nodes)
        all(haskey(primary_nodes, name) for name in
            (:dc_positive, :dc_negative, :ac_1, :ac_2)) &&
            all(haskey(secondary_nodes, name) for name in
                (:dc_positive, :dc_negative, :ac_1, :ac_2)) || throw(ArgumentError(
                    "DAB bridge topology is missing a canonical terminal",
                ))
        input_port = only(port for port in specification.ports if port.identity === :primary_dc)
        output_port = only(port for port in specification.ports if port.identity === :secondary_dc)
        input_port.kind === IsolatedDirectCurrentPort &&
            output_port.kind === IsolatedDirectCurrentPort || throw(ArgumentError(
                "DAB ports must declare isolated direct-current ownership",
            ))
        input_port.ordered_nodes ==
            (primary_nodes[:dc_positive], primary_nodes[:dc_negative]) &&
            output_port.ordered_nodes ==
                (secondary_nodes[:dc_positive], secondary_nodes[:dc_negative]) ||
            throw(ArgumentError("DAB DC ports do not match their exact bridge rails"))
        isempty(intersect(
            Set(filter(!=(0), getfield.(primary_topology.nodes, :node))),
            Set(filter(!=(0), getfield.(secondary_topology.nodes, :node))),
        )) || throw(ArgumentError(
            "DAB primary and secondary bridge terminals must be galvanically distinct",
        ))
        modulation = specification.modulation
        primary_inner = modulation.primary_inner_phase_shift_rad
        secondary_inner = modulation.secondary_inner_phase_shift_rad
        all(angle -> 0.0 <= angle < pi, (primary_inner, secondary_inner)) &&
            abs(modulation.phase_shift_rad) <= pi || throw(ArgumentError(
                "DAB outer and inner phase shifts exceed their admitted domains",
            ))
        modulation.kind === SinglePhaseShiftModulation &&
            (!iszero(primary_inner) || !iszero(secondary_inner)) &&
            throw(ArgumentError("single-phase-shift DAB requires zero inner shifts"))
        modulation.kind === DualPhaseShiftModulation &&
            (iszero(primary_inner) == iszero(secondary_inner)) &&
            throw(ArgumentError("dual-phase-shift DAB requires one nonzero inner shift"))
        modulation.kind === TriplePhaseShiftModulation &&
            (iszero(primary_inner) || iszero(secondary_inner)) &&
            throw(ArgumentError("triple-phase-shift DAB requires two nonzero inner shifts"))
        timing = specification.timing
        timing.dead_time_s >= timing.fixed_step_s || throw(ArgumentError(
            "DAB commutation dead time must span at least one fixed step",
        ))
        if selection.fidelity === SwitchingDetailed
            timing.fixed_step_s <= minimum(min(
                device.recovered_charge_lifetime_s,
                device.turn_off_tail_time_s,
            ) for device in semiconductor) / 10.0 || throw(ArgumentError(
                "DAB timestep must resolve both bridges' recovery and tail-current state",
            ))
        end
        carrier_steps = inv(timing.carrier_frequency_hz * timing.fixed_step_s)
        isapprox(carrier_steps, round(carrier_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("DAB carrier period must contain integer fixed steps"))
        for angle in (
            modulation.phase_shift_rad,
            modulation.primary_inner_phase_shift_rad,
            modulation.secondary_inner_phase_shift_rad,
        )
            phase_steps = abs(angle) / (2.0 * pi) * carrier_steps
            isapprox(phase_steps, round(phase_steps); atol=1.0e-10, rtol=1.0e-10) ||
                throw(ArgumentError("DAB phase-shift edges must lie on the fixed-step calendar"))
        end
        dead_time_steps = timing.dead_time_s / timing.fixed_step_s
        isapprox(dead_time_steps, round(dead_time_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("DAB dead time must lie on the fixed-step calendar"))
        horizon_steps = (values[6] - values[5]) / timing.fixed_step_s
        isapprox(horizon_steps, round(horizon_steps); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError("DAB horizon must contain integer fixed steps"))
        if selection.fidelity === SwitchingDetailed
            all(device ->
                device.switching_energy_current_domain_maximum_a >=
                    specification.rated_bases.current_a &&
                device.switching_energy_blocking_voltage_domain_maximum_v >=
                    max(values[1], values[2]), semiconductor) || throw(ArgumentError(
                        "DAB switching-energy domains must include rated current and both DC voltages",
                    ))
        end
        return new(
            specification,
            primary_topology,
            secondary_topology,
            values[1:4]...,
            transformer_link,
            initial_state,
            semiconductor,
            values[5:6]...,
        )
    end
end

struct SwitchingDualActiveBridgeTrace
    time_s::Vector{Float64}
    primary_dc_voltage_v::Vector{Float64}
    secondary_dc_voltage_v::Vector{Float64}
    primary_dc_current_a::Vector{Float64}
    secondary_dc_current_a::Vector{Float64}
    primary_bridge_voltage_v::Vector{Float64}
    secondary_bridge_voltage_v::Vector{Float64}
    leakage_current_a::Vector{Float64}
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

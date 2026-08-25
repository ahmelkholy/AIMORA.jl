export AverageDualActiveBridgeStudy,
       AverageDualActiveBridgeTrace

struct AverageDualActiveBridgeStudy
    specification::ConverterSystemSpecification
    primary_topology::BridgeTopologyDescriptor
    secondary_topology::BridgeTopologyDescriptor
    primary_dc_voltage_v::Float64
    secondary_dc_voltage_v::Float64
    primary_source_resistance_ohm::Float64
    secondary_source_resistance_ohm::Float64
    transformer_link::DualActiveBridgeTransformerLink
    start_time_s::Float64
    stop_time_s::Float64

    function AverageDualActiveBridgeStudy(
        specification::ConverterSystemSpecification;
        primary_topology::BridgeTopologyDescriptor,
        secondary_topology::BridgeTopologyDescriptor,
        primary_dc_voltage_v::Real,
        secondary_dc_voltage_v::Real,
        primary_source_resistance_ohm::Real,
        secondary_source_resistance_ohm::Real,
        transformer_link::DualActiveBridgeTransformerLink,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === DualActiveBridge || throw(ArgumentError(
            "average DAB study requires the canonical dual-active-bridge family",
        ))
        selection.fidelity === AverageValue || throw(ArgumentError(
            "average DAB study requires explicit AverageValue fidelity",
        ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "average DAB study is a standalone conversion owner",
        ))
        primary_topology.family === :full_bridge &&
            secondary_topology.family === :full_bridge || throw(ArgumentError(
                "average DAB requires two canonical B200 full bridges",
            ))
        specification.topology_signatures == (
            bridge_topology_signature(primary_topology),
            bridge_topology_signature(secondary_topology),
        ) || throw(ArgumentError(
            "average DAB specification does not bind its exact bridge topologies",
        ))
        isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
            "average DAB execution cannot claim switch-detailed device identities",
        ))
        specification.passive_and_transformer_signatures ==
            (transformer_link.deterministic_signature_sha256,) || throw(ArgumentError(
                "average DAB specification does not bind its transformer link",
            ))
        specification.modulation.kind === SinglePhaseShiftModulation ||
            throw(ArgumentError(
                "average DAB fundamental-transfer execution requires single-phase-shift modulation",
            ))
        modulation = specification.modulation
        iszero(modulation.primary_inner_phase_shift_rad) &&
            iszero(modulation.secondary_inner_phase_shift_rad) || throw(ArgumentError(
                "average DAB fundamental-transfer execution requires zero inner phase shifts",
            ))
        specification.timing.dead_time_s == 0.0 &&
            specification.timing.minimum_pulse_s == 0.0 || throw(ArgumentError(
                "average DAB execution does not admit switching-edge timing",
            ))
        iszero(transformer_link.leakage_resistance_ohm) || throw(ArgumentError(
            "average DAB fundamental-transfer execution requires the admitted lossless leakage limit",
        ))
        values = Float64.((
            primary_dc_voltage_v,
            secondary_dc_voltage_v,
            primary_source_resistance_ohm,
            secondary_source_resistance_ohm,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:4]) || throw(ArgumentError(
            "average DAB voltages and source resistances must be finite and positive",
        ))
        values[5] >= 0.0 && values[6] > values[5] || throw(ArgumentError(
            "average DAB stop time must follow its nonnegative start time",
        ))
        step_count = (values[6] - values[5]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "average DAB horizon must contain an integer number of fixed steps",
            ))
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("average DAB specification is not ready"))
        angular_frequency = 2.0 * pi * specification.timing.carrier_frequency_hz
        power = dual_active_bridge_power_w(
            values[1],
            values[2],
            transformer_link.turns_ratio,
            modulation.phase_shift_rad,
            angular_frequency,
            transformer_link.leakage_inductance_h,
        )
        values[1]^2 - 4.0 * values[3] * power >= 0.0 &&
            values[2]^2 + 4.0 * values[4] * power >= 0.0 || throw(ArgumentError(
                "average DAB requested transfer exceeds a source-resistance power limit",
            ))
        return new(
            specification,
            primary_topology,
            secondary_topology,
            values[1:4]...,
            transformer_link,
            values[5:6]...,
        )
    end
end

struct AverageDualActiveBridgeTrace
    time_s::Vector{Float64}
    primary_dc_voltage_v::Vector{Float64}
    secondary_dc_voltage_v::Vector{Float64}
    primary_dc_current_a::Vector{Float64}
    secondary_dc_current_a::Vector{Float64}
    transferred_power_w::Vector{Float64}
    source_dissipated_energy_j::Vector{Float64}
    energy_residual_w::Vector{Float64}
    result::ConverterSystemResult
end

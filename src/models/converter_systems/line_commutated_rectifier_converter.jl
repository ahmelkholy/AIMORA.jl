export LineCommutatedRectifierInitialState,
       SwitchingLineCommutatedRectifierStudy,
       SwitchingLineCommutatedRectifierTrace,
       detailed_rectifier_semiconductor_signatures

struct LineCommutatedRectifierInitialState
    dc_current_a::Float64

    function LineCommutatedRectifierInitialState(dc_current_a::Real=0.0)
        current = Float64(dc_current_a)
        isfinite(current) && current >= 0.0 || throw(ArgumentError(
            "line-commutated rectifier initial DC current must be finite and nonnegative",
        ))
        return new(current)
    end
end

function detailed_rectifier_semiconductor_signatures(
    parameters::DetailedConverterSemiconductorParameters,
    valve_count::Integer,
)
    count = Int(valve_count)
    count in (4, 6, 12, 18, 24) || throw(ArgumentError(
        "line-commutated rectifier detailed fidelity requires a supported valve count",
    ))
    body = repr(parameters)
    return ntuple(count) do position
        bytes2hex(sha256("line_commutated_rectifier_valve_$(position)\n" * body))
    end
end

struct SwitchingLineCommutatedRectifierStudy
    specification::ConverterSystemSpecification
    topology::BridgeTopologyDescriptor
    source_peak_voltage_v::Float64
    source_frequency_hz::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::LineCommutatedRectifierInitialState
    detailed_semiconductor::Union{Nothing,DetailedConverterSemiconductorParameters}
    group_phase_shift_rad::Vector{Float64}
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingLineCommutatedRectifierStudy(
        specification::ConverterSystemSpecification;
        topology::BridgeTopologyDescriptor,
        source_peak_voltage_v::Real,
        source_frequency_hz::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::LineCommutatedRectifierInitialState,
        detailed_semiconductor::Union{
            Nothing,
            DetailedConverterSemiconductorParameters,
        }=nothing,
        group_phase_shift_rad=nothing,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family in (
            SinglePhaseDiodeBridge,
            ThreePhaseDiodeBridge,
            SinglePhaseThyristorBridge,
            ThreePhaseThyristorBridge,
            SinglePhaseHalfControlledBridge,
            ThreePhaseHalfControlledBridge,
            MultipulseDiodeBridge,
            MultipulseThyristorBridge,
        ) ||
            throw(ArgumentError(
                "line-commutated rectifier study requires a supported single- or three-phase bridge family",
            ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "natural rectifier execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "natural rectifier study is a standalone conversion owner",
        ))
        single_phase = selection.family in (
            SinglePhaseDiodeBridge,
            SinglePhaseThyristorBridge,
            SinglePhaseHalfControlledBridge,
        )
        multipulse = selection.family in (
            MultipulseDiodeBridge,
            MultipulseThyristorBridge,
        )
        expected_topology_family = multipulse ? :multi_group_bridge :
            single_phase ? :single_phase_graetz : :polyphase_bridge
        topology.family === expected_topology_family || throw(ArgumentError(
            "natural rectifier family does not match its B200 topology",
        ))
        expected_valve_count = multipulse ? selection.pulse_count : single_phase ? 4 : 6
        expected_classes = if selection.family in (
            SinglePhaseDiodeBridge,
            ThreePhaseDiodeBridge,
            MultipulseDiodeBridge,
        )
            fill(:diode, expected_valve_count)
        elseif selection.family in (
            SinglePhaseThyristorBridge,
            ThreePhaseThyristorBridge,
            MultipulseThyristorBridge,
        )
            fill(:thyristor, expected_valve_count)
        else
            [isodd(index) ? :thyristor : :diode for index in 1:expected_valve_count]
        end
        length(topology.valve_positions) == expected_valve_count &&
            getfield.(topology.valve_positions, :valve_class) == expected_classes ||
            throw(ArgumentError(
                "line-commutated rectifier topology does not match its exact valve classes",
            ))
        specification.topology_signatures == (bridge_topology_signature(topology),) ||
            throw(ArgumentError(
                "natural rectifier specification does not bind its exact topology",
            ))
        semiconductor = if selection.fidelity === SwitchingDetailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed natural rectifier requires semiconductor fidelity",
            ))
            specification.device_fidelity_signatures ==
                detailed_rectifier_semiconductor_signatures(
                    detailed_semiconductor,
                    expected_valve_count,
                ) || throw(ArgumentError(
                    "natural rectifier specification does not bind every detailed valve",
                ))
            detailed_semiconductor
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state natural rectifier cannot accept detailed semiconductor data",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state natural rectifier cannot claim detailed device identities",
            ))
            nothing
        end
        expected_modulation = selection.family in (
            SinglePhaseDiodeBridge,
            ThreePhaseDiodeBridge,
            MultipulseDiodeBridge,
        ) ? NaturalDiodeCommutation : PhaseControlledFiring
        specification.modulation.kind === expected_modulation || throw(ArgumentError(
            "line-commutated rectifier modulation does not match its valve family",
        ))
        if expected_modulation === PhaseControlledFiring
            0.0 <= specification.modulation.firing_angle_rad <= pi ||
                throw(ArgumentError(
                    "line-commutated rectifier firing angle must lie from zero through pi",
                ))
        end
        group_count = multipulse ? selection.pulse_count ÷ 6 : 1
        phase_shifts = group_phase_shift_rad === nothing ?
            (multipulse ? throw(ArgumentError(
                "multipulse rectifier requires one explicit transformer phase shift per group",
            )) : [0.0]) : Float64.(group_phase_shift_rad)
        length(phase_shifts) == group_count && all(isfinite, phase_shifts) ||
            throw(ArgumentError(
                "line-commutated rectifier group phase shifts are incomplete or nonfinite",
            ))
        iszero(first(phase_shifts)) || throw(ArgumentError(
            "line-commutated rectifier first group must define the zero phase reference",
        ))
        length(unique(phase_shifts)) == length(phase_shifts) || throw(ArgumentError(
            "line-commutated rectifier group phase shifts must be unique",
        ))
        values = Float64.((
            source_peak_voltage_v,
            source_frequency_hz,
            source_resistance_ohm,
            load_resistance_ohm,
            load_inductance_h,
            start_time_s,
            stop_time_s,
        ))
        all(isfinite, values) && all(>(0.0), values[1:5]) || throw(ArgumentError(
            "natural rectifier source and load parameters must be finite and positive",
        ))
        values[6] >= 0.0 && values[7] > values[6] || throw(ArgumentError(
            "natural rectifier stop time must follow its nonnegative start time",
        ))
        isapprox(specification.timing.firing_frequency_hz, values[2];
            atol=0.0, rtol=16eps(Float64)) || throw(ArgumentError(
                "natural rectifier timing must bind the AC source frequency",
            ))
        step_count = (values[7] - values[6]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "natural rectifier horizon must contain integer fixed steps",
            ))
        if selection.fidelity === SwitchingDetailed
            specification.timing.fixed_step_s <= min(
                semiconductor.recovered_charge_lifetime_s,
                semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "natural rectifier timestep must resolve detailed stored-charge state",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("natural rectifier specification is not ready"))
        return new(
            specification,
            topology,
            values[1:5]...,
            initial_state,
            semiconductor,
            phase_shifts,
            values[6:7]...,
        )
    end
end

struct SwitchingLineCommutatedRectifierTrace
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
    result::ConverterSystemResult
end

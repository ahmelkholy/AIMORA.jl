export ConverterSystemReadiness,
       converter_family_category,
       converter_executable_fidelities,
       converter_supported_fidelities,
       converter_supported_modulations,
       converter_system_readiness,
       converter_system_is_ready

struct ConverterSystemReadiness
    accepted::Bool
    reason::Symbol
    message::String
    family::ConverterSystemFamily
    fidelity::ModelFidelity
    application::ConverterApplication
    specification_signature_sha256::String
end

converter_system_is_ready(readiness::ConverterSystemReadiness) = readiness.accepted

function converter_family_category(family::ConverterSystemFamily)
    family in (
        SinglePhaseDiodeBridge,
        ThreePhaseDiodeBridge,
        SinglePhaseThyristorBridge,
        ThreePhaseThyristorBridge,
        SinglePhaseHalfControlledBridge,
        ThreePhaseHalfControlledBridge,
        MultipulseDiodeBridge,
        MultipulseThyristorBridge,
    ) && return :ac_dc
    family in (
        BuckChopper,
        BoostChopper,
        InvertingBuckBoostChopper,
        FourQuadrantChopper,
        InterleavedChopper,
        DualActiveBridge,
    ) && return :dc_dc
    family in (
        SinglePhaseTwoLevelBridge,
        ThreePhaseTwoLevelBridge,
        ThreeLevelNeutralPointClampedBridge,
        ThreeLevelTTypeBridge,
        FlyingCapacitorBridge,
        CascadedHBridge,
    ) && return :dc_ac
    return :direct_ac_ac
end

const _AVERAGE_VALUE_FAMILIES = (
    BuckChopper,
    BoostChopper,
    InvertingBuckBoostChopper,
    FourQuadrantChopper,
    DualActiveBridge,
    SinglePhaseTwoLevelBridge,
    ThreePhaseTwoLevelBridge,
    ThreeLevelNeutralPointClampedBridge,
)

function converter_supported_fidelities(
    family::ConverterSystemFamily,
    application::ConverterApplication=StandaloneConversion,
)
    application !== StandaloneConversion && return family === ThreePhaseTwoLevelBridge ?
        (AverageValue,) : ()
    if family in _AVERAGE_VALUE_FAMILIES
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    end
    return (SwitchingStateEquivalent, SwitchingDetailed)
end

function converter_executable_fidelities(
    family::ConverterSystemFamily,
    application::ConverterApplication=StandaloneConversion,
)
    application === StandaloneConversion || return family === ThreePhaseTwoLevelBridge ?
        (AverageValue,) : ()
    family in (
        SinglePhaseDiodeBridge,
        ThreePhaseDiodeBridge,
        SinglePhaseThyristorBridge,
        ThreePhaseThyristorBridge,
        SinglePhaseHalfControlledBridge,
        ThreePhaseHalfControlledBridge,
        MultipulseDiodeBridge,
        MultipulseThyristorBridge,
    ) &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    family === BuckChopper && return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === BoostChopper && return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === InvertingBuckBoostChopper &&
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === FourQuadrantChopper &&
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === InterleavedChopper && return (SwitchingStateEquivalent, SwitchingDetailed)
    family === DualActiveBridge &&
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family in (SinglePhaseTwoLevelBridge, ThreePhaseTwoLevelBridge) &&
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === ThreeLevelNeutralPointClampedBridge &&
        return (AverageValue, SwitchingStateEquivalent, SwitchingDetailed)
    family === ThreeLevelTTypeBridge &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    family === FlyingCapacitorBridge &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    family === CascadedHBridge &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    family === ThreePhaseMatrixConverter &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    family === LineCommutatedCycloconverter &&
        return (SwitchingStateEquivalent, SwitchingDetailed)
    return ()
end

function converter_supported_modulations(family::ConverterSystemFamily)
    family in (SinglePhaseDiodeBridge, ThreePhaseDiodeBridge, MultipulseDiodeBridge) &&
        return (NaturalDiodeCommutation,)
    family in (
        SinglePhaseThyristorBridge,
        ThreePhaseThyristorBridge,
        SinglePhaseHalfControlledBridge,
        ThreePhaseHalfControlledBridge,
        MultipulseThyristorBridge,
    ) && return (PhaseControlledFiring,)
    family in (BuckChopper, BoostChopper, InvertingBuckBoostChopper, FourQuadrantChopper) &&
        return (CarrierSinusoidalPulseWidthModulation,)
    family === InterleavedChopper && return (PhaseShiftedCarrierPulseWidthModulation,)
    family === DualActiveBridge && return (
        SinglePhaseShiftModulation,
        DualPhaseShiftModulation,
        TriplePhaseShiftModulation,
    )
    family in (SinglePhaseTwoLevelBridge, ThreePhaseTwoLevelBridge) && return (
        CarrierSinusoidalPulseWidthModulation,
        SpaceVectorPulseWidthModulation,
        SelectiveHarmonicElimination,
    )
    family in (
        ThreeLevelNeutralPointClampedBridge,
        ThreeLevelTTypeBridge,
        FlyingCapacitorBridge,
        CascadedHBridge,
    ) && return (
        CarrierSinusoidalPulseWidthModulation,
        PhaseShiftedCarrierPulseWidthModulation,
        SelectiveHarmonicElimination,
        NearestLevelModulation,
    )
    family === ThreePhaseMatrixConverter && return (MatrixSpaceVectorModulation,)
    return (CycloconverterFiringSynthesis,)
end

function _converter_refusal(specification, reason, message)
    selection = specification.selection
    return ConverterSystemReadiness(
        false,
        reason,
        message,
        selection.family,
        selection.fidelity,
        selection.application,
        specification.signature_sha256,
    )
end

function _validate_converter_timing(specification)
    timing = specification.timing
    values = (
        timing.fixed_step_s,
        timing.control_period_s,
        timing.carrier_frequency_hz,
        timing.firing_frequency_hz,
        timing.dead_time_s,
        timing.minimum_pulse_s,
    )
    all(isfinite, values) || return "converter timing values must be finite"
    timing.fixed_step_s > 0.0 || return "fixed EMT step must be positive"
    timing.control_period_s > 0.0 || return "control period must be positive"
    all(>=(0.0), values[3:end]) || return "carrier, firing, dead-time, and pulse values must be nonnegative"
    if specification.selection.fidelity !== AverageValue
        active_frequency = max(timing.carrier_frequency_hz, timing.firing_frequency_hz)
        active_frequency > 0.0 || return "switching execution requires a positive carrier or firing frequency"
        timing.fixed_step_s <= inv(40.0 * active_frequency) || return "fixed step must provide at least forty points per carrier or firing period"
    end
    timing.dead_time_s + timing.minimum_pulse_s <= timing.control_period_s ||
        return "dead time and minimum pulse exceed the control period"
    return nothing
end

function _validate_converter_counts(selection::ConverterSystemSelection)
    family = selection.family
    expected_phase_count = family in (
        SinglePhaseDiodeBridge,
        SinglePhaseThyristorBridge,
        SinglePhaseHalfControlledBridge,
        SinglePhaseTwoLevelBridge,
    ) ? (1,) : family in (
        BuckChopper,
        BoostChopper,
        InvertingBuckBoostChopper,
        FourQuadrantChopper,
        InterleavedChopper,
        DualActiveBridge,
    ) ? (1,) : family === LineCommutatedCycloconverter ? (1, 3) : (3,)
    selection.phase_count in expected_phase_count || return "phase count is incompatible with the selected family"
    if family in (MultipulseDiodeBridge, MultipulseThyristorBridge)
        selection.pulse_count in (12, 18, 24) || return "multipulse systems require 12, 18, or 24 pulses"
    elseif family in (ThreePhaseDiodeBridge, ThreePhaseThyristorBridge, ThreePhaseHalfControlledBridge)
        selection.pulse_count == 6 || return "three-phase bridge systems require six pulses"
    elseif family in (SinglePhaseDiodeBridge, SinglePhaseThyristorBridge, SinglePhaseHalfControlledBridge)
        selection.pulse_count == 2 || return "single-phase bridge systems require two pulses"
    end
    family === InterleavedChopper && !(2 <= selection.channel_count <= 8) &&
        return "interleaved choppers require two through eight channels"
    family !== InterleavedChopper && selection.channel_count != 1 &&
        return "channel count is only selectable for an interleaved chopper"
    family === CascadedHBridge && !(2 <= selection.cell_count <= 8) &&
        return "cascaded H-bridges require two through eight cells"
    family !== CascadedHBridge && selection.cell_count != 1 &&
        return "cell count is only selectable for a cascaded H-bridge"
    1 <= selection.thermal_stage_count <= 8 || return "device thermal-stage count must be one through eight"
    return nothing
end

function _validate_application_composition(specification)
    application = specification.selection.application
    application === StandaloneConversion ||
        specification.selection.family === ThreePhaseTwoLevelBridge ||
        return "released converter applications require the three-phase two-level grid-interface family"
    topology_count = length(specification.topology_signatures)
    minimum_count = application === StandaloneConversion ? 1 :
        application in (ShuntActiveHarmonicFilter, DynamicVoltageRestorer) ? 1 :
        application in (DoubleConversionUninterruptiblePowerSupply, ConductiveElectricVehicleCharger) ? 2 : 3
    topology_count >= minimum_count || return "application composition omits a required converter stage"
    application === DynamicVoltageRestorer &&
        isempty(specification.passive_and_transformer_signatures) &&
        return "dynamic-voltage-restorer composition requires an explicit series transformer owner"
    application === ThreeStageSolidStateTransformer &&
        isempty(specification.passive_and_transformer_signatures) &&
        return "solid-state-transformer composition requires an explicit transformer owner"
    return nothing
end

function converter_system_readiness(specification::ConverterSystemSpecification)
    selection = specification.selection
    if !isempty(specification.event_commands)
        selection.family in (BuckChopper, BoostChopper, InvertingBuckBoostChopper) &&
            selection.fidelity !== AverageValue || return _converter_refusal(
                specification,
                :unsupported_event_domain,
                "converter block/restart events are released only for physical switching choppers",
            )
        all(command -> command.kind in (ConverterBlockEvent, ConverterRestartEvent),
            specification.event_commands) || return _converter_refusal(
                specification,
                :unsupported_event_kind,
                "converter event calendar contains an unreleased event kind",
            )
    end
    selection.fidelity in converter_supported_fidelities(selection.family, selection.application) ||
        return _converter_refusal(
            specification,
            :unsupported_fidelity,
            "selected converter family and application do not release the requested fidelity",
        )
    selection.fidelity in converter_executable_fidelities(
        selection.family,
        selection.application,
    ) || return _converter_refusal(
        specification,
        :execution_unavailable,
        "selected converter family, application, and fidelity do not yet have a canonical executable converter-system runtime",
    )
    specification.modulation.kind in converter_supported_modulations(selection.family) ||
        return _converter_refusal(
            specification,
            :unsupported_modulation,
            "selected converter family does not release the requested modulation or commutation method",
        )
    timing_failure = _validate_converter_timing(specification)
    timing_failure === nothing || return _converter_refusal(
        specification,
        :invalid_timing,
        timing_failure,
    )
    count_failure = _validate_converter_counts(selection)
    count_failure === nothing || return _converter_refusal(
        specification,
        :invalid_family_count,
        count_failure,
    )
    application_failure = _validate_application_composition(specification)
    application_failure === nothing || return _converter_refusal(
        specification,
        :incomplete_application_composition,
        application_failure,
    )
    selection.fidelity === SwitchingDetailed && isempty(specification.device_fidelity_signatures) &&
        return _converter_refusal(
            specification,
            :missing_device_fidelity,
            "switching-detailed execution requires selected D200 device-fidelity owners",
        )
    modulation = specification.modulation
    0.0 < modulation.duty < 1.0 || return _converter_refusal(
        specification,
        :invalid_duty,
        "converter duty must lie strictly between zero and one",
    )
    isfinite(modulation.modulation_index) && modulation.modulation_index >= 0.0 ||
        return _converter_refusal(specification, :invalid_modulation_index, "modulation index must be finite and nonnegative")
    all(isfinite, (
        modulation.firing_angle_rad,
        modulation.phase_shift_rad,
        modulation.primary_inner_phase_shift_rad,
        modulation.secondary_inner_phase_shift_rad,
    )) || return _converter_refusal(
        specification,
        :invalid_angle,
        "firing, outer phase-shift, and inner phase-shift angles must be finite",
    )
    return ConverterSystemReadiness(
        true,
        :ready,
        "converter system is ready for its declared backend capability",
        selection.family,
        selection.fidelity,
        selection.application,
        specification.signature_sha256,
    )
end

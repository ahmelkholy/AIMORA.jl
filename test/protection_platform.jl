using Test
using AIMORA.ProtectionStudy
using AIMORA.MeasurementChains: MeasurementAcquisitionSettings,
                                MeasurementChainSpecification,
                                MeasurementSample,
                                ThreePhaseSampledMeasurement
using AIMORA.StudyCore: NumericalPolicyParameter,
                        ParameterProvenance,
                        PhysicalModelParameter
using AIMORA.EMTTaskPlatform: CarrierEMTTask,
                              InvalidateEMTTopology,
                              ProtectionEMTTask,
                              emt_logical_time

const _PROTECTION_SETTING_PROVENANCE = ParameterProvenance(
    "AIMORA-authored generic protection setting",
    "declared relay engineering unit",
    "supplied directly in SI engineering units",
    "synthetic setting with no field uncertainty claim",
    "generic fixed-step magnitude relay tests",
    PhysicalModelParameter,
)

const _PROTECTION_TIMER_PROVENANCE = ParameterProvenance(
    "AIMORA-authored exact relay timer policy",
    "seconds and dimensionless fraction",
    "finite Float64 setting on exact accepted sample ticks",
    "deterministic numerical policy",
    "instantaneous, definite, or generic inverse timer",
    NumericalPolicyParameter,
)

const _PROTECTION_MEASUREMENT_PROVENANCE = ParameterProvenance(
    "AIMORA-authored generic sampled-current input",
    "ampere",
    "supplied directly as accepted primary RMS and instantaneous SI values",
    "synthetic deterministic input",
    "three-phase protection measurement adapter test",
    PhysicalModelParameter,
)

function _magnitude_settings(;
    direction=OverMagnitude,
    value_mode=AbsoluteMagnitudeValue,
    orientation_polarity=1,
    timer_mode=ProtectionTimerInstantaneous,
    definite_time_s=0.0,
    reset_time_s=0.0,
)
    return MagnitudeRelaySettings(
        :phase_overcurrent,
        :feeder_current,
        :phase_a;
        stage=ProtectionSlidingRMSStage,
        value_mode,
        orientation_polarity,
        direction,
        timer_mode,
        pickup=10.0,
        dropout_ratio=0.9,
        definite_time_s,
        inverse_a=1.0,
        inverse_b=0.0,
        inverse_p=1.0,
        time_dial_s=1.0,
        reset_time_s,
        unit="A",
        setting_provenance=_PROTECTION_SETTING_PROVENANCE,
        timer_provenance=_PROTECTION_TIMER_PROVENANCE,
    )
end

function _protection_measurement(tick, value; unavailable_reason=nothing)
    return ProtectionMeasurement(
        :feeder_current,
        :phase_a,
        :current,
        "A",
        "positive into protected feeder",
        ProtectionSlidingRMSStage,
        tick,
        tick,
        0.01,
        value,
        unavailable_reason === nothing ? :valid : :window_incomplete,
        unavailable_reason,
        repeat("a", 64),
    )
end

function _phasor_measurement(id, channel, unit, tick, value; reason=nothing)
    return ProtectionMeasurement(
        id,
        channel,
        unit == "V" ? :voltage : :current,
        unit,
        "positive into protected zone",
        ProtectionFundamentalPhasorStage,
        tick,
        tick,
        0.001,
        value,
        reason === nothing ? :valid : :phasor_unavailable,
        reason,
        repeat("d", 64),
    )
end

function _instantaneous_measurement(id, channel, unit, tick, value; reason=nothing)
    return ProtectionMeasurement(
        id,
        channel,
        unit == "V" ? :voltage : :current,
        unit,
        "positive into protected line terminal",
        ProtectionInstantaneousStage,
        tick,
        tick,
        1.0e-6,
        value,
        reason === nothing ? :valid : :measurement_unavailable,
        reason,
        repeat("e", 64),
    )
end

@testset "causal magnitude protection element" begin
    instantaneous = _magnitude_settings()
    state = MagnitudeRelayState()
    below = evaluate_magnitude_relay!(state, instantaneous, _protection_measurement(0, 9.0))
    @test !below.pickup_active
    @test !below.operated
    pickup = evaluate_magnitude_relay!(state, instantaneous, _protection_measurement(1, 10.0))
    @test pickup.pickup_active
    @test pickup.operated
    @test pickup.pickup_multiple == 1.0
    @test state.operation_count == 1
    hysteresis = evaluate_magnitude_relay!(state, instantaneous, _protection_measurement(2, 9.5))
    @test hysteresis.pickup_active
    @test hysteresis.operated
    dropout = evaluate_magnitude_relay!(state, instantaneous, _protection_measurement(3, 8.9))
    @test !dropout.pickup_active
    @test !dropout.operated
    @test_throws ArgumentError evaluate_magnitude_relay!(
        state,
        instantaneous,
        _protection_measurement(3, 12.0),
    )

    definite = _magnitude_settings(
        timer_mode=ProtectionDefiniteTimer,
        definite_time_s=0.02,
    )
    definite_state = MagnitudeRelayState()
    first = evaluate_magnitude_relay!(definite_state, definite, _protection_measurement(0, 20.0))
    @test first.timer_fraction == 0.0
    @test !first.operated
    halfway = evaluate_magnitude_relay!(definite_state, definite, _protection_measurement(1, 20.0))
    @test halfway.timer_fraction ≈ 0.5
    @test !halfway.operated
    complete = evaluate_magnitude_relay!(definite_state, definite, _protection_measurement(2, 20.0))
    @test complete.timer_fraction == 1.0
    @test complete.operated

    inverse = _magnitude_settings(timer_mode=ProtectionInverseTimer)
    inverse_state = MagnitudeRelayState()
    evaluate_magnitude_relay!(inverse_state, inverse, _protection_measurement(0, 20.0))
    for tick in 1:99
        decision = evaluate_magnitude_relay!(
            inverse_state,
            inverse,
            _protection_measurement(tick, 20.0),
        )
        @test !decision.operated
    end
    inverse_operation = evaluate_magnitude_relay!(
        inverse_state,
        inverse,
        _protection_measurement(100, 20.0),
    )
    @test inverse_operation.operated
    @test inverse_operation.timer_fraction == 1.0

    under = _magnitude_settings(direction=UnderMagnitude)
    under_state = MagnitudeRelayState()
    @test evaluate_magnitude_relay!(
        under_state,
        under,
        _protection_measurement(0, 10.0),
    ).operated
    @test evaluate_magnitude_relay!(
        under_state,
        under,
        _protection_measurement(1, 10.5),
    ).operated
    @test !evaluate_magnitude_relay!(
        under_state,
        under,
        _protection_measurement(2, 11.2),
    ).operated

    falling_rate = _magnitude_settings(
        value_mode=SignedScalarValue,
        orientation_polarity=-1,
    )
    falling_state = MagnitudeRelayState()
    @test evaluate_magnitude_relay!(
        falling_state,
        falling_rate,
        _protection_measurement(0, -12.0),
    ).operated
    @test_throws ArgumentError evaluate_magnitude_relay!(
        MagnitudeRelayState(),
        falling_rate,
        ProtectionMeasurement(
            :feeder_current,
            :phase_a,
            :current,
            "A",
            "positive into protected feeder",
            ProtectionSlidingRMSStage,
            0,
            0,
            0.01,
            12.0 + 1.0im,
            :valid,
            nothing,
            repeat("a", 64),
        ),
    )

    blocked_state = MagnitudeRelayState()
    blocked = evaluate_magnitude_relay!(
        blocked_state,
        definite,
        _protection_measurement(0, nothing; unavailable_reason=:window_incomplete),
    )
    @test blocked.blocked_reason == :window_incomplete
    @test blocked.measured_value === nothing
    @test blocked_state.blocked_count == 1
    @test blocked_state.timer_fraction == 0.0

    snapshot = magnitude_relay_snapshot(inverse_state, inverse)
    restored = MagnitudeRelayState()
    restore_magnitude_relay_snapshot!(restored, inverse, snapshot)
    @test restored.pickup_active == inverse_state.pickup_active
    @test restored.operated == inverse_state.operated
    @test restored.timer_fraction == inverse_state.timer_fraction
    @test restored.last_release_tick == inverse_state.last_release_tick
    @test restored.operation_count == inverse_state.operation_count
    @test length(inverse_operation.deterministic_signature_sha256) == 64
end

@testset "accepted measurement to protection adapter" begin
    acquisition = MeasurementAcquisitionSettings(
        tick_s=1.0e-3,
        sample_period_ticks=1,
        window_weights_newest_first=ones(4),
        nominal_frequency_hz=50.0,
        maximum_retained_samples=4,
    )
    specification = MeasurementChainSpecification(
        :feeder_current,
        ThreePhaseSampledMeasurement;
        channel_names=[:phase_a, :phase_b, :phase_c],
        quantity=:current,
        unit="A",
        orientation="positive into protected feeder",
        phase_order=[:a, :b, :c],
        acquisition,
        minimum_input=-100.0,
        maximum_input=100.0,
        maximum_spectral_frequency_hz=100.0,
        maximum_timestep_s=1.0e-3,
        provenance=_PROTECTION_MEASUREMENT_PROVENANCE,
    )
    sample = MeasurementSample(
        10,
        12,
        0.010,
        0.012,
        [1.0, -2.0, 3.0],
        nothing,
        falses(3),
        [4.0, 5.0, 6.0],
        ComplexF64[1.0 + 2.0im, -2.0 + 3.0im, 4.0 - 1.0im],
        (zero=0.1 + 0.2im, positive=4.0 + 3.0im, negative=0.3 - 0.1im),
        49.9,
        :valid,
        repeat("b", 64),
    )
    instantaneous = protection_measurement(
        sample,
        specification,
        :phase_b;
        stage=ProtectionInstantaneousStage,
    )
    @test instantaneous.value == -2.0 + 0.0im
    @test instantaneous.source_tick == 10
    @test instantaneous.release_tick == 12
    rms = protection_measurement(
        sample,
        specification,
        :phase_c;
        stage=ProtectionSlidingRMSStage,
    )
    @test rms.value == 6.0 + 0.0im
    phasor = protection_measurement(
        sample,
        specification,
        :phase_a;
        stage=ProtectionFundamentalPhasorStage,
    )
    @test phasor.value == 1.0 + 2.0im
    sequence = protection_measurement(
        sample,
        specification,
        :phase_a;
        stage=ProtectionSequencePhasorStage,
        sequence=:positive,
    )
    @test sequence.channel == :positive
    @test sequence.value == 4.0 + 3.0im
    frequency = protection_measurement(
        sample,
        specification,
        :phase_a;
        stage=ProtectionFrequencyStage,
    )
    @test frequency.channel == :frequency
    @test frequency.unit == "Hz"
    @test frequency.value == 49.9 + 0.0im
    @test length(frequency.deterministic_signature_sha256) == 64

    clipped_sample = MeasurementSample(
        sample.source_tick,
        sample.release_tick,
        sample.source_time_s,
        sample.release_time_s,
        sample.instantaneous,
        sample.codes,
        BitVector([false, true, false]),
        sample.sliding_rms,
        sample.fundamental_rms_phasors,
        sample.sequence_phasors,
        sample.frequency_hz,
        sample.quality,
        sample.deterministic_signature_sha256,
    )
    clipped_frequency = protection_measurement(
        clipped_sample,
        specification,
        :phase_a;
        stage=ProtectionFrequencyStage,
    )
    @test clipped_frequency.value === nothing
    @test clipped_frequency.unavailable_reason == :measurement_clipped
end

@testset "causal frequency and ROCOF protection" begin
    settings = ROCOFEstimatorSettings(
        :feeder_rocof,
        :feeder_frequency;
        sample_period_ticks=1,
        window_intervals=2,
        minimum_frequency_hz=40.0,
        maximum_frequency_hz=70.0,
        orientation="positive for rising frequency",
        numerical_provenance=_PROTECTION_TIMER_PROVENANCE,
    )
    state = ROCOFEstimatorState()
    frequency(tick, value; reason=nothing) = ProtectionMeasurement(
        :feeder_frequency,
        :frequency,
        :frequency,
        "Hz",
        "positive frequency",
        ProtectionFrequencyStage,
        tick,
        tick,
        0.01,
        value,
        reason === nothing ? :valid : :frequency_unavailable,
        reason,
        repeat("c", 64),
    )
    first = estimate_rocof!(state, settings, frequency(0, 50.0))
    @test first.value === nothing
    @test first.unavailable_reason == :rocof_window_incomplete
    second = estimate_rocof!(state, settings, frequency(1, 49.9))
    @test second.value === nothing
    estimate = estimate_rocof!(state, settings, frequency(2, 49.8))
    @test real(estimate.value) ≈ -10.0 atol=1.0e-12
    @test estimate.stage === ProtectionROCOFStage
    @test state.estimate_count == 1

    falling = MagnitudeRelaySettings(
        :underfrequency_rate,
        :feeder_rocof,
        :rocof;
        stage=ProtectionROCOFStage,
        value_mode=SignedScalarValue,
        orientation_polarity=-1,
        direction=OverMagnitude,
        pickup=5.0,
        unit="Hz/s",
        setting_provenance=_PROTECTION_SETTING_PROVENANCE,
        timer_provenance=_PROTECTION_TIMER_PROVENANCE,
    )
    @test evaluate_magnitude_relay!(
        MagnitudeRelayState(),
        falling,
        estimate,
    ).operated

    snapshot = rocof_estimator_snapshot(state, settings)
    restored = ROCOFEstimatorState()
    restore_rocof_estimator_snapshot!(restored, settings, snapshot)
    @test restored.source_ticks == state.source_ticks
    @test restored.frequencies_hz == state.frequencies_hz
    @test restored.estimate_count == state.estimate_count

    missing = estimate_rocof!(
        state,
        settings,
        frequency(3, nothing; reason=:frequency_unavailable),
    )
    @test missing.value === nothing
    @test isempty(state.source_ticks)
    after_gap = estimate_rocof!(state, settings, frequency(5, 49.7))
    @test after_gap.unavailable_reason == :rocof_window_incomplete
end

@testset "directional distance differential and incremental-wave elements" begin
    directional_settings = DirectionalRelaySettings(
        :forward_power_direction,
        :line_voltage,
        :phase_a,
        :line_current,
        :phase_a;
        characteristic_angle_rad=0.0,
        minimum_polarizing_voltage_v=10.0,
        minimum_operating_torque_w=1.0,
        memory_decay_per_tick=0.99,
        maximum_memory_age_ticks=2,
        provenance=_PROTECTION_SETTING_PROVENANCE,
    )
    directional_state = DirectionalRelayState()
    forward = evaluate_directional_relay!(
        directional_state,
        directional_settings,
        _phasor_measurement(:line_voltage, :phase_a, "V", 0, 100.0 + 0.0im),
        _phasor_measurement(:line_current, :phase_a, "A", 0, 10.0 + 0.0im),
    )
    @test forward.forward
    @test forward.operating_torque_w == 1_000.0
    reverse = evaluate_directional_relay!(
        directional_state,
        directional_settings,
        _phasor_measurement(:line_voltage, :phase_a, "V", 1, 100.0 + 0.0im),
        _phasor_measurement(:line_current, :phase_a, "A", 1, -10.0 + 0.0im),
    )
    @test !reverse.forward
    @test reverse.operating_torque_w == -1_000.0
    memory = evaluate_directional_relay!(
        directional_state,
        directional_settings,
        _phasor_measurement(:line_voltage, :phase_a, "V", 2, 1.0 + 0.0im),
        _phasor_measurement(:line_current, :phase_a, "A", 2, 10.0 + 0.0im),
    )
    @test memory.used_memory_polarization
    @test memory.polarizing_voltage_v ≈ 99.0 + 0.0im
    directional_snapshot = directional_relay_snapshot(
        directional_state,
        directional_settings,
    )
    restored_directional = DirectionalRelayState()
    restore_directional_relay_snapshot!(
        restored_directional,
        directional_settings,
        directional_snapshot,
    )
    uninterrupted_directional = evaluate_directional_relay!(
        directional_state,
        directional_settings,
        _phasor_measurement(:line_voltage, :phase_a, "V", 3, 1.0 + 0.0im),
        _phasor_measurement(:line_current, :phase_a, "A", 3, 10.0 + 0.0im),
    )
    restarted_directional = evaluate_directional_relay!(
        restored_directional,
        directional_settings,
        _phasor_measurement(:line_voltage, :phase_a, "V", 3, 1.0 + 0.0im),
        _phasor_measurement(:line_current, :phase_a, "A", 3, 10.0 + 0.0im),
    )
    @test restarted_directional == uninterrupted_directional

    mho_settings = DistanceRelaySettings(
        :line_zone_one,
        :loop_voltage,
        :ab,
        :loop_current,
        :ab,
        MhoDistanceZone(5.0 + 0.0im, 5.0);
        minimum_loop_current_a=0.1,
        provenance=_PROTECTION_SETTING_PROVENANCE,
    )
    mho_state = DistanceRelayState()
    inside = evaluate_distance_relay!(
        mho_state,
        mho_settings,
        _phasor_measurement(:loop_voltage, :ab, "V", 0, 50.0 + 0.0im),
        _phasor_measurement(:loop_current, :ab, "A", 0, 10.0 + 0.0im),
    )
    @test inside.asserted
    @test inside.apparent_impedance_ohm == 5.0 + 0.0im
    @test inside.zone_margin_ohm == 5.0
    outside = evaluate_distance_relay!(
        mho_state,
        mho_settings,
        _phasor_measurement(:loop_voltage, :ab, "V", 1, 120.0 + 0.0im),
        _phasor_measurement(:loop_current, :ab, "A", 1, 10.0 + 0.0im),
    )
    @test !outside.asserted
    @test outside.zone_margin_ohm == -2.0
    distance_snapshot = distance_relay_snapshot(mho_state, mho_settings)
    restored_distance = DistanceRelayState()
    restore_distance_relay_snapshot!(restored_distance, mho_settings, distance_snapshot)
    uninterrupted_distance = evaluate_distance_relay!(
        mho_state,
        mho_settings,
        _phasor_measurement(:loop_voltage, :ab, "V", 2, 60.0 + 0.0im),
        _phasor_measurement(:loop_current, :ab, "A", 2, 10.0 + 0.0im),
    )
    restarted_distance = evaluate_distance_relay!(
        restored_distance,
        mho_settings,
        _phasor_measurement(:loop_voltage, :ab, "V", 2, 60.0 + 0.0im),
        _phasor_measurement(:loop_current, :ab, "A", 2, 10.0 + 0.0im),
    )
    @test restarted_distance == uninterrupted_distance

    polygon = PolygonDistanceZone(
        ComplexF64[1.0, -1.0, 1.0im, -1.0im],
        [10.0, 0.0, 5.0, 5.0],
    )
    polygon_settings = DistanceRelaySettings(
        :line_polygon_zone,
        :loop_voltage,
        :ag,
        :loop_current,
        :ag,
        polygon;
        minimum_loop_current_a=0.1,
        provenance=_PROTECTION_SETTING_PROVENANCE,
    )
    polygon_decision = evaluate_distance_relay!(
        DistanceRelayState(),
        polygon_settings,
        _phasor_measurement(:loop_voltage, :ag, "V", 0, 40.0 + 20.0im),
        _phasor_measurement(:loop_current, :ag, "A", 0, 10.0 + 0.0im),
    )
    @test polygon_decision.asserted
    @test polygon_decision.apparent_impedance_ohm == 4.0 + 2.0im

    differential_settings = DifferentialRelaySettings(
        :bus_differential,
        [:terminal_one_current, :terminal_two_current],
        [:phase_a, :phase_a];
        compensation=ComplexF64[1.0, 1.0],
        restraint_mode=HalfSumDifferentialRestraint,
        minimum_operate_a=1.0,
        initial_bias_a=0.0,
        restraint_breakpoints_a=[10.0],
        region_slopes=[0.2, 0.5],
        provenance=_PROTECTION_SETTING_PROVENANCE,
    )
    differential_state = DifferentialRelayState()
    external_fault = evaluate_differential_relay!(
        differential_state,
        differential_settings,
        ProtectionMeasurement[
            _phasor_measurement(:terminal_one_current, :phase_a, "A", 0, 10.0),
            _phasor_measurement(:terminal_two_current, :phase_a, "A", 0, -10.0),
        ],
    )
    @test !external_fault.asserted
    @test external_fault.operate_current_a == 0.0
    @test external_fault.restraint_current_a == 10.0
    internal_fault = evaluate_differential_relay!(
        differential_state,
        differential_settings,
        ProtectionMeasurement[
            _phasor_measurement(:terminal_one_current, :phase_a, "A", 1, 10.0),
            _phasor_measurement(:terminal_two_current, :phase_a, "A", 1, 5.0),
        ],
    )
    @test internal_fault.asserted
    @test internal_fault.operate_current_a == 15.0
    @test internal_fault.restraint_current_a == 7.5
    @test internal_fault.threshold_current_a == 1.5
    differential_snapshot = differential_relay_snapshot(
        differential_state,
        differential_settings,
    )
    restored_differential = DifferentialRelayState()
    restore_differential_relay_snapshot!(
        restored_differential,
        differential_settings,
        differential_snapshot,
    )
    next_terminal_currents = ProtectionMeasurement[
        _phasor_measurement(:terminal_one_current, :phase_a, "A", 2, 8.0),
        _phasor_measurement(:terminal_two_current, :phase_a, "A", 2, -8.0),
    ]
    uninterrupted_differential = evaluate_differential_relay!(
        differential_state,
        differential_settings,
        next_terminal_currents,
    )
    restarted_differential = evaluate_differential_relay!(
        restored_differential,
        differential_settings,
        next_terminal_currents,
    )
    @test all(
        field -> getfield(restarted_differential, field) ==
            getfield(uninterrupted_differential, field),
        fieldnames(DifferentialRelayDecision),
    )

    wave_settings = IncrementalWaveSettings(
        :line_forward_wave,
        :terminal_voltage,
        :phase_a,
        :terminal_current,
        :phase_a;
        sample_period_ticks=1,
        reference_impedance_ohm=10.0,
        direction=ForwardIncrementalWave,
        polarity=1,
        threshold_v=15.0,
        provenance=_PROTECTION_SETTING_PROVENANCE,
    )
    wave_state = IncrementalWaveState()
    warmup = evaluate_incremental_wave!(
        wave_state,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 0, 100.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 0, 5.0),
    )
    @test warmup.blocked_reason == :incremental_history_unavailable
    wave = evaluate_incremental_wave!(
        wave_state,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 1, 110.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 1, 6.0),
    )
    @test wave.forward_wave_v == 20.0
    @test wave.reverse_wave_v == 0.0
    @test wave.asserted
    wave_snapshot = incremental_wave_snapshot(wave_state, wave_settings)
    restored_wave = IncrementalWaveState()
    restore_incremental_wave_snapshot!(restored_wave, wave_settings, wave_snapshot)
    uninterrupted_wave = evaluate_incremental_wave!(
        wave_state,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 2, 120.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 2, 7.0),
    )
    restarted_wave = evaluate_incremental_wave!(
        restored_wave,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 2, 120.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 2, 7.0),
    )
    @test restarted_wave == uninterrupted_wave
    calendar_gap = evaluate_incremental_wave!(
        wave_state,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 4, 130.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 4, 8.0),
    )
    @test calendar_gap.blocked_reason == :incremental_calendar_gap
    @test evaluate_incremental_wave!(
        wave_state,
        wave_settings,
        _instantaneous_measurement(:terminal_voltage, :phase_a, "V", 5, 140.0),
        _instantaneous_measurement(:terminal_current, :phase_a, "A", 5, 9.0),
    ).asserted
end

@testset "deterministic protection logic and communication" begin
    definition = ProtectionLogicDefinition(
        :permissive_trip_logic,
        AbstractProtectionLogicNode[
            ProtectionLogicInputNode(:local_pickup, :local_pickup),
            ProtectionLogicInputNode(:remote_permissive, :remote_permissive),
            ProtectionLogicVoteNode(
                :trip_condition,
                [:local_pickup, :remote_permissive],
                2,
            ),
            ProtectionLogicInputNode(:manual_reset, :manual_reset),
            ProtectionLogicLatchNode(
                :trip_latch,
                :trip_condition,
                :manual_reset;
                reset_dominant=true,
            ),
        ],
        :trip_latch,
    )
    logic_runtime = ProtectionLogicRuntime(definition)
    blocked = evaluate_protection_logic!(logic_runtime, definition, Dict(
        :local_pickup => true,
        :remote_permissive => false,
        :manual_reset => false,
    ))
    @test !blocked.output
    trip = evaluate_protection_logic!(logic_runtime, definition, Dict(
        :local_pickup => true,
        :remote_permissive => true,
        :manual_reset => false,
    ))
    @test trip.output
    @test trip.values[:trip_condition]
    retained = evaluate_protection_logic!(logic_runtime, definition, Dict(
        :local_pickup => false,
        :remote_permissive => false,
        :manual_reset => false,
    ))
    @test retained.output
    reset = evaluate_protection_logic!(logic_runtime, definition, Dict(
        :local_pickup => true,
        :remote_permissive => true,
        :manual_reset => true,
    ))
    @test !reset.output
    @test length(reset.deterministic_signature_sha256) == 64
    logic_snapshot = protection_logic_snapshot(logic_runtime, definition)
    restored_logic = ProtectionLogicRuntime(definition)
    restore_protection_logic_snapshot!(restored_logic, definition, logic_snapshot)
    next_inputs = Dict(
        :local_pickup => true,
        :remote_permissive => true,
        :manual_reset => false,
    )
    restarted_logic = evaluate_protection_logic!(restored_logic, definition, next_inputs)
    uninterrupted_logic = evaluate_protection_logic!(logic_runtime, definition, next_inputs)
    @test restarted_logic.values == uninterrupted_logic.values
    @test restarted_logic.output == uninterrupted_logic.output
    @test restarted_logic.evaluation_index == uninterrupted_logic.evaluation_index
    @test restarted_logic.deterministic_signature_sha256 ==
        uninterrupted_logic.deterministic_signature_sha256
    @test_throws ArgumentError ProtectionLogicDefinition(
        :cyclic_logic,
        AbstractProtectionLogicNode[
            ProtectionLogicNotNode(:first, :second),
            ProtectionLogicNotNode(:second, :first),
        ],
        :second,
    )

    link = ProtectionCommunicationLink(
        :line_permissive_link,
        :relay_local,
        :relay_remote;
        allowed_payloads=[:permissive, :blocking, :direct_trip],
        fixed_delay_ticks=2,
        rules=ProtectionMessageRule[
            ProtectionMessageRule(2; disposition=DropProtectionMessage),
            ProtectionMessageRule(3; copy_count=2, additional_delay_ticks=2),
        ],
        provenance=_PROTECTION_TIMER_PROVENANCE,
    )
    communication = ProtectionCommunicationRuntime()
    first_send = send_protection_message!(communication, link, :permissive, 0)
    @test first_send.sequence_number == 1
    @test first_send.delivery_tick == 2
    dropped_send = send_protection_message!(communication, link, :blocking, 0)
    @test dropped_send.disposition === DropProtectionMessage
    @test dropped_send.delivery_tick === nothing
    duplicated_send = send_protection_message!(communication, link, :direct_trip, 1)
    @test duplicated_send.queued_copy_count == 2
    @test duplicated_send.delivery_tick == 5
    snapshot = protection_communication_snapshot(communication, link)
    restored_communication = ProtectionCommunicationRuntime()
    restore_protection_communication_snapshot!(restored_communication, link, snapshot)
    @test restored_communication.queue == communication.queue
    due_first = deliver_due_protection_messages!(communication, link, 2)
    @test length(due_first) == 1
    @test only(due_first).payload == :permissive
    due_duplicates = deliver_due_protection_messages!(communication, link, 5)
    @test length(due_duplicates) == 2
    @test getfield.(due_duplicates, :copy_index) == [1, 2]
    @test communication.sent_count == 3
    @test communication.dropped_count == 1
    @test communication.delivered_count == 3
    @test isempty(communication.queue)
    @test_throws ArgumentError send_protection_message!(
        communication,
        link,
        :unknown_payload,
        6,
    )
end

@testset "exact protection task order" begin
    no_op = (_state, _instant, _activation_index) -> nothing
    operations = ProtectionTaskOperations(
        no_op,
        no_op,
        no_op,
        no_op,
        no_op,
        no_op,
        no_op,
    )
    pipeline = ProtectionTaskPipeline(
        operations,
        (accepted=false,);
        epoch=emt_logical_time(0),
        period=emt_logical_time(1 // 1_000),
        start=emt_logical_time(0),
        stop=emt_logical_time(2 // 1_000),
    )
    plan = protection_task_plan(pipeline)
    ordered_entries = plan.entries[plan.execution_order]
    @test getfield.(getfield.(ordered_entries, :spec), :name) == [
        "protection_measurement_release",
        "protection_element_evaluation",
        "protection_local_logic",
        "protection_message_send",
        "protection_message_delivery",
        "protection_remote_logic",
        "protection_trip_command",
    ]
    @test getfield.(getfield.(ordered_entries, :spec), :family) == [
        ProtectionEMTTask,
        ProtectionEMTTask,
        ProtectionEMTTask,
        CarrierEMTTask,
        CarrierEMTTask,
        ProtectionEMTTask,
        ProtectionEMTTask,
    ]
    @test only(last(ordered_entries).spec.effects) === InvalidateEMTTopology
    @test length(pipeline.deterministic_signature_sha256) == 64
    if !AIMORA.solver_available()
        unavailable = AIMORA.prepare_protection_task_pipeline(pipeline)
        @test unavailable isa AIMORA.SolverUnavailableResult
        @test unavailable.required_capability == :emt_protection_breaker
    end
end

@testset "three-pole EMT breaker failure reclose and restart" begin
    specification = EMTBreakerSpecification(
        :feeder_breaker;
        closed_conductance_s=1.0e6,
        open_conductance_s=1.0e-6,
        opening_travel_ticks=2,
        closing_travel_ticks=1,
        current_zero_required=true,
        current_zero_threshold_a=0.1,
        failure_delay_ticks=3,
        failure_current_threshold_a=1.0,
        reclose_dead_ticks=2,
        reclaim_ticks=3,
        maximum_reclose_shots=1,
        physical_provenance=_PROTECTION_SETTING_PROVENANCE,
        timing_provenance=_PROTECTION_TIMER_PROVENANCE,
    )
    failed_breaker = EMTBreakerRuntime(specification; tick_s=1.0e-3)
    @test request_breaker_trip!(failed_breaker, specification, 0)
    @test breaker_pole_conductances(failed_breaker, specification) ==
        (1.0e6, 1.0e6, 1.0e6)
    advance_emt_breaker!(
        failed_breaker,
        specification,
        0,
        (100.0, 100.0, 100.0),
        (10.0, 10.0, 10.0),
    )
    awaiting = advance_emt_breaker!(
        failed_breaker,
        specification,
        2,
        (100.0, 100.0, 100.0),
        (10.0, 10.0, 10.0),
    )
    @test all(==(BreakerPoleAwaitingCurrentZero), awaiting.pole_positions)
    failure = advance_emt_breaker!(
        failed_breaker,
        specification,
        3,
        (100.0, 100.0, 100.0),
        (10.0, 10.0, 10.0),
    )
    @test failure.failure_asserted
    @test failure.backup_trip_required
    opened_after_failure = advance_emt_breaker!(
        failed_breaker,
        specification,
        4,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0);
        current_zero_reached=(true, true, true),
    )
    @test all(==(BreakerPoleOpen), opened_after_failure.pole_positions)
    @test sum(opened_after_failure.contact_energy_j) > 0.0
    @test_throws ArgumentError request_breaker_reclose!(
        failed_breaker,
        specification,
        4,
    )

    successful_breaker = EMTBreakerRuntime(specification; tick_s=1.0e-3)
    request_breaker_trip!(successful_breaker, specification, 0)
    advance_emt_breaker!(
        successful_breaker,
        specification,
        0,
        (100.0, 100.0, 100.0),
        (10.0, 10.0, 10.0),
    )
    opened = advance_emt_breaker!(
        successful_breaker,
        specification,
        2,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0);
        current_zero_reached=(true, true, true),
    )
    @test opened.reclose_available
    @test request_breaker_reclose!(successful_breaker, specification, 2) == 4
    trip_success = advance_emt_breaker!(
        successful_breaker,
        specification,
        3,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    )
    @test !trip_success.failure_asserted
    closing = advance_emt_breaker!(
        successful_breaker,
        specification,
        4,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    )
    @test all(==(BreakerPoleClosing), closing.pole_positions)
    closed = advance_emt_breaker!(
        successful_breaker,
        specification,
        5,
        (100.0, 100.0, 100.0),
        (1.0, 1.0, 1.0),
    )
    @test all(==(BreakerPoleClosed), closed.pole_positions)
    @test closed.shot_count == 1
    snapshot = emt_breaker_snapshot(successful_breaker, specification)
    restored_breaker = EMTBreakerRuntime(specification; tick_s=1.0e-3)
    restore_emt_breaker_snapshot!(restored_breaker, specification, snapshot)
    @test getfield.(restored_breaker.poles, :position) ==
        getfield.(successful_breaker.poles, :position)
    @test getfield.(restored_breaker.poles, :contact_energy_j) ==
        getfield.(successful_breaker.poles, :contact_energy_j)
    @test restored_breaker.shot_count == successful_breaker.shot_count

    request_breaker_trip!(successful_breaker, specification, 6)
    @test successful_breaker.lockout
    advance_emt_breaker!(
        successful_breaker,
        specification,
        6,
        (100.0, 100.0, 100.0),
        (10.0, 10.0, 10.0),
    )
    locked_open = advance_emt_breaker!(
        successful_breaker,
        specification,
        8,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0);
        current_zero_reached=(true, true, true),
    )
    @test locked_open.lockout
    @test !locked_open.reclose_available
    @test_throws ArgumentError request_breaker_reclose!(
        successful_breaker,
        specification,
        8,
    )
end

@testset "generic public protection product contracts" begin
    terminal = ProtectionTerminalDefinition(
        :source_terminal,
        :radial_feeder_source,
        [:phase_a_current, :phase_b_current, :phase_c_current];
        inward_current_orientation="positive from source into the protected feeder",
    )
    product_logic = ProtectionLogicDefinition(
        :generic_radial_product_logic,
        AbstractProtectionLogicNode[
            ProtectionLogicInputNode(:trip_input, :trip_input),
        ],
        :trip_input,
    )
    product_breaker = EMTBreakerSpecification(
        :feeder_breaker;
        closed_conductance_s=1.0e6,
        open_conductance_s=1.0e-9,
        opening_travel_ticks=1,
        closing_travel_ticks=1,
        current_zero_required=false,
        current_zero_threshold_a=0.0,
        failure_delay_ticks=2,
        failure_current_threshold_a=0.1,
        reclose_dead_ticks=2,
        reclaim_ticks=3,
        maximum_reclose_shots=1,
        physical_provenance=_PROTECTION_SETTING_PROVENANCE,
        timing_provenance=_PROTECTION_TIMER_PROVENANCE,
    )
    specification = ProtectionProductSpecification(
        :generic_radial_feeder_protection,
        RadialFeederProtectionProduct,
        :radial_feeder,
        :feeder_zone;
        terminals=[terminal],
        measurement_products=[:feeder_phase_currents, :feeder_residual_current],
        element_families=[
            :phase_overcurrent,
            :residual_earth_fault,
            :breaker_failure,
            :autoreclose,
        ],
        trip_breakers=[:feeder_breaker],
        network_timestep_s=10.0e-6,
        network_timestep_logical=emt_logical_time(1 // 100_000),
        protection_sample_period_ticks=10,
        setting_provenance=_PROTECTION_SETTING_PROVENANCE,
        timing_provenance=_PROTECTION_TIMER_PROVENANCE,
        uncertainty="synthetic settings; field and manufacturer uncertainty unknown",
        validity_domain="generic radial AC feeder from 45 through 65 Hz",
        configuration=(
            element_settings=(_magnitude_settings(),),
            logic_definitions=(product_logic,),
            communication_links=(),
            breaker_specifications=(product_breaker,),
        ),
    )
    preparation = prepare_protection_product(specification)
    readiness = protection_product_readiness(preparation)
    @test readiness.ready
    @test readiness.code == :ready_for_solver_free_inspection
    @test !readiness.production_backend_available
    @test readiness.product == specification.id
    @test readiness.deterministic_signature_sha256 ==
        preparation.preparation_signature_sha256
    @test :certification in readiness.unsupported_phenomena
    @test length(specification.deterministic_signature_sha256) == 64

    runtime = ProtectionProductRuntime(
        preparation;
        task_plan_signature_sha256=repeat("f", 64),
    )
    evaluate_magnitude_relay!(
        only(runtime.element_states),
        only(specification.configuration.element_settings),
        _protection_measurement(0, 12.0),
    )
    evaluate_protection_logic!(
        only(runtime.logic_states),
        only(specification.configuration.logic_definitions),
        Dict(:trip_input => true),
    )
    request_breaker_trip!(
        only(runtime.breaker_states),
        only(specification.configuration.breaker_specifications),
        0,
    )
    runtime.output_cursor = 1
    push!(runtime.event_trace, ProtectionProductEvent(0, :feeder_breaker, :trip, :phase_a))
    runtime_snapshot = protection_product_runtime_snapshot(runtime)
    @test length(runtime_snapshot.deterministic_signature_sha256) == 64
    runtime.output_cursor = 9
    empty!(runtime.event_trace)
    restore_protection_product_runtime_snapshot!(runtime, runtime_snapshot)
    @test runtime.output_cursor == 1
    @test runtime.event_trace == runtime_snapshot.event_trace
    @test protection_product_runtime_snapshot(runtime).deterministic_signature_sha256 ==
        runtime_snapshot.deterministic_signature_sha256
    mktempdir() do directory
        snapshot_path = joinpath(directory, "generic-radial-product.aimora")
        project_signature = repeat("a", 64)
        topology_signature = repeat("b", 64)
        descriptor = write_protection_product_portable_snapshot(
            snapshot_path,
            runtime;
            project_signature_sha256=project_signature,
            topology_signature_sha256=topology_signature,
        )
        inspected = AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(snapshot_path)
        @test inspected.content_sha256 == descriptor.content_sha256
        @test inspected.metadata.profile == :portable_public_reference
        @test inspected.metadata.accepted_step == runtime.accepted_tick
        @test only(inspected.sections).identity == "protection.product_runtime"

        restored_runtime = ProtectionProductRuntime(
            preparation;
            task_plan_signature_sha256=repeat("f", 64),
        )
        restored_runtime, restored_descriptor =
            restore_protection_product_portable_snapshot!(
                restored_runtime,
                snapshot_path;
                project_signature_sha256=project_signature,
                topology_signature_sha256=topology_signature,
            )
        @test restored_descriptor.content_sha256 == descriptor.content_sha256
        restored_snapshot = protection_product_runtime_snapshot(restored_runtime)
        @test restored_snapshot.deterministic_signature_sha256 ==
            runtime_snapshot.deterministic_signature_sha256
        @test repr(restored_snapshot) == repr(runtime_snapshot)

        original_decision = evaluate_magnitude_relay!(
            only(runtime.element_states),
            only(specification.configuration.element_settings),
            _protection_measurement(10, 13.0),
        )
        restored_decision = evaluate_magnitude_relay!(
            only(restored_runtime.element_states),
            only(specification.configuration.element_settings),
            _protection_measurement(10, 13.0),
        )
        runtime.accepted_tick = 10
        restored_runtime.accepted_tick = 10
        runtime.output_cursor += 1
        restored_runtime.output_cursor += 1
        @test repr(original_decision) == repr(restored_decision)
        continued_original = protection_product_runtime_snapshot(runtime)
        continued_restored = protection_product_runtime_snapshot(restored_runtime)
        @test continued_restored.deterministic_signature_sha256 ==
            continued_original.deterministic_signature_sha256
        @test repr(continued_restored) == repr(continued_original)

        @test_throws ArgumentError restore_protection_product_portable_snapshot!(
            deepcopy(restored_runtime),
            snapshot_path;
            project_signature_sha256=repeat("c", 64),
            topology_signature_sha256=topology_signature,
        )
        @test_throws ArgumentError restore_protection_product_portable_snapshot!(
            deepcopy(restored_runtime),
            snapshot_path;
            project_signature_sha256=project_signature,
            topology_signature_sha256=repeat("d", 64),
        )
        changed_preparation = prepare_protection_product(
            specification;
            initial_tick=1,
        )
        @test_throws ArgumentError restore_protection_product_portable_snapshot!(
            ProtectionProductRuntime(
                changed_preparation;
                task_plan_signature_sha256=repeat("f", 64),
            ),
            snapshot_path;
            project_signature_sha256=project_signature,
            topology_signature_sha256=topology_signature,
        )
        corrupted = read(snapshot_path)
        corrupted[end - 40] = xor(corrupted[end - 40], 0x01)
        corrupted_path = joinpath(directory, "corrupted.aimora")
        write(corrupted_path, corrupted)
        @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure begin
            restore_protection_product_portable_snapshot!(
                deepcopy(restored_runtime),
                corrupted_path;
                project_signature_sha256=project_signature,
                topology_signature_sha256=topology_signature,
            )
        end
    end
    stale_runtime = ProtectionProductRuntime(
        preparation;
        task_plan_signature_sha256=repeat("e", 64),
    )
    @test_throws ArgumentError restore_protection_product_runtime_snapshot!(
        stale_runtime,
        runtime_snapshot,
    )

    declared_preparation = prepare_protection_product(
        specification;
        initialization_mode=DeclaredProtectionStateInitialization,
        initially_closed_breakers=Symbol[],
        initially_locked_out_breakers=[:feeder_breaker],
    )
    declared_runtime = ProtectionProductRuntime(
        declared_preparation;
        task_plan_signature_sha256=repeat("f", 64),
    )
    @test only(declared_runtime.breaker_states).lockout
    @test all(
        ==(BreakerPoleOpen),
        getfield.(only(declared_runtime.breaker_states).poles, :position),
    )
    operating_point_preparation = prepare_protection_product(
        specification;
        initialization_mode=OperatingPointProtectionInitialization,
        initial_tick=5,
    )
    operating_point_runtime = ProtectionProductRuntime(
        operating_point_preparation;
        task_plan_signature_sha256=repeat("f", 64),
    )
    @test operating_point_runtime.accepted_tick == 5
    @test all(
        ==(BreakerPoleClosed),
        getfield.(only(operating_point_runtime.breaker_states).poles, :position),
    )

    no_op = (_state, _instant, _activation_index) -> nothing
    study_pipeline = ProtectionTaskPipeline(
        ProtectionTaskOperations(no_op, no_op, no_op, no_op, no_op, no_op, no_op),
        (product_runtime=runtime,);
        epoch=emt_logical_time(0),
        period=emt_logical_time(1 // 10_000),
        start=emt_logical_time(0),
        stop=emt_logical_time(0),
    )
    invalid_study = prepare_protection_study(
        preparation,
        study_pipeline;
        execution_instant=emt_logical_time(0),
        topology_signature_sha256="not-a-sha256",
        output_ids=[:trip_command],
    )
    @test invalid_study isa ProtectionStudyRefusal
    @test invalid_study.code == :invalid_preparation
    valid_study = prepare_protection_study(
        preparation,
        study_pipeline;
        execution_instant=emt_logical_time(0),
        topology_signature_sha256=repeat("0", 64),
        output_ids=[:trip_command],
    )
    @test valid_study isa ProtectionStudyPreparation
    backend_result = run_protection(valid_study)
    if AIMORA.solver_available()
        @test backend_result isa ProtectionStudyResult
    else
        @test backend_result isa ProtectionStudyRefusal
        @test backend_result.code == :production_backend_unavailable
        @test backend_result.last_accepted_tick == preparation.initial_tick
    end

    @test_throws UndefKeywordError ProtectionProductSpecification(
        :incomplete_line_protection,
        DirectionalDistanceLineProtectionProduct,
        :line,
        :line_zone;
        terminals=[terminal],
        measurement_products=[:line_measurements],
        element_families=[:directional, :distance, :incremental_wave],
        trip_breakers=[:line_breaker],
        network_timestep_s=10.0e-6,
        network_timestep_logical=emt_logical_time(1 // 100_000),
        protection_sample_period_ticks=1,
        setting_provenance=_PROTECTION_SETTING_PROVENANCE,
        timing_provenance=_PROTECTION_TIMER_PROVENANCE,
        uncertainty="synthetic",
        validity_domain="generic line",
    )
    @test_throws ArgumentError prepare_protection_product(
        specification;
        initially_locked_out_breakers=[:feeder_breaker],
    )
end

@testset "converter-system family and fidelity contracts" begin
    converters = AIMORA.ConverterSystems
    core = AIMORA.StudyCore
    provenance = core.ParameterProvenance(
        "AIMORA generic converter contract",
        "SI",
        "identity",
        "synthetic bounded parameters",
        "declared family-specific domain",
        core.PhysicalModelParameter,
    )
    validity = core.ModelValidityDomain(
        :generic_dc_chopper;
        description="finite fixed-step generic DC chopper domain",
    )
    ports = (
        converters.ConverterPortDefinition(
            :input_dc,
            converters.DirectCurrentPort,
            (1, 0);
            voltage_orientation="positive from input positive to input negative",
        ),
        converters.ConverterPortDefinition(
            :output_dc,
            converters.DirectCurrentPort,
            (2, 0);
            voltage_orientation="positive from output positive to output negative",
        ),
    )
    selection = converters.ConverterSystemSelection(
        converters.BuckChopper,
        core.AverageValue;
        phase_count=1,
    )
    specification = converters.ConverterSystemSpecification(
        :generic_buck,
        selection,
        ports,
        ("step_down_chopper_topology",),
        (),
        ("input_inductor_output_capacitor",),
        converters.ConverterRatedBases(800.0, 100.0, 80_000.0, 10_000.0),
        converters.ConverterTimingParameters(
            fixed_step_s=1.0e-6,
            control_period_s=100.0e-6,
            carrier_frequency_hz=10_000.0,
        ),
        converters.ConverterModulationParameters(
            kind=converters.CarrierSinusoidalPulseWidthModulation,
            duty=0.4,
        ),
        (provenance,),
        validity,
    )
    readiness = converters.converter_system_readiness(specification)
    @test converters.converter_system_is_ready(readiness)
    @test readiness.reason == :ready
    @test converters.converter_executable_fidelities(converters.BuckChopper) ==
        (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(converters.BoostChopper) ==
        (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.InvertingBuckBoostChopper,
    ) == (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(converters.FourQuadrantChopper) ==
        (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(converters.InterleavedChopper) ==
        (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(converters.DualActiveBridge) ==
        (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.SinglePhaseTwoLevelBridge,
    ) == (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.ThreePhaseTwoLevelBridge,
    ) == (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.ThreeLevelNeutralPointClampedBridge,
    ) == (core.AverageValue, core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.ThreePhaseMatrixConverter,
    ) == (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.LineCommutatedCycloconverter,
    ) == (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.ThreeLevelTTypeBridge,
    ) == (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.FlyingCapacitorBridge,
    ) == (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    @test converters.converter_executable_fidelities(
        converters.CascadedHBridge,
    ) == (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    for family in (
        converters.SinglePhaseDiodeBridge,
        converters.ThreePhaseDiodeBridge,
        converters.SinglePhaseThyristorBridge,
        converters.ThreePhaseThyristorBridge,
        converters.SinglePhaseHalfControlledBridge,
        converters.ThreePhaseHalfControlledBridge,
        converters.MultipulseDiodeBridge,
        converters.MultipulseThyristorBridge,
    )
        @test converters.converter_executable_fidelities(family) ==
            (core.SwitchingStateEquivalent, core.SwitchingDetailed)
    end
    @test length(converters.converter_system_signature(specification)) == 64
    @test converters.converter_family_category(converters.BuckChopper) == :dc_dc
    @test converters.ideal_dcdc_conversion_ratio(converters.BuckChopper, 0.4) == 0.4
    @test converters.ideal_dcdc_conversion_ratio(converters.BoostChopper, 0.4) ≈ 5 / 3
    @test converters.ideal_dcdc_conversion_ratio(
        converters.InvertingBuckBoostChopper,
        0.4,
    ) ≈ -2 / 3
    boost_state = converters.boost_converter_operating_point(100.0, 0.4, 0.1, 0.2, 20.0)
    transfer_fraction = 0.6
    @test 100.0 - 0.3 * boost_state.inductor_current_a -
        transfer_fraction * boost_state.output_voltage_v ≈ 0.0 atol=1.0e-12
    @test transfer_fraction * boost_state.inductor_current_a -
        boost_state.output_voltage_v / 20.0 ≈ 0.0 atol=1.0e-12
    lossless_boost = converters.boost_converter_operating_point(
        100.0,
        0.4,
        eps(Float64),
        0.0,
        20.0,
    )
    @test lossless_boost.output_voltage_v / 100.0 ≈ 5 / 3 rtol=1.0e-14
    inverting_state = converters.inverting_buck_boost_operating_point(
        100.0,
        0.4,
        0.1,
        0.2,
        20.0,
    )
    @test 0.4 * (100.0 - 0.3 * inverting_state.inductor_current_a) +
        0.6 * (inverting_state.output_voltage_v -
               0.2 * inverting_state.inductor_current_a) ≈ 0.0 atol=1.0e-12
    @test 0.4 * (-inverting_state.output_voltage_v / 20.0) +
        0.6 * (-inverting_state.inductor_current_a -
               inverting_state.output_voltage_v / 20.0) ≈ 0.0 atol=1.0e-12
    lossless_inverting = converters.inverting_buck_boost_operating_point(
        100.0,
        0.4,
        eps(Float64),
        0.0,
        20.0,
    )
    @test lossless_inverting.output_voltage_v / 100.0 ≈ -2 / 3 rtol=1.0e-14
    positive_four_quadrant = converters.average_four_quadrant_operating_point(
        100.0,
        0.6,
        0.1,
        1.0,
    )
    negative_four_quadrant = converters.average_four_quadrant_operating_point(
        100.0,
        0.4,
        0.1,
        1.0,
    )
    @test positive_four_quadrant.load_current_a ≈ -negative_four_quadrant.load_current_a
    @test positive_four_quadrant.load_current_a > 0.0

    unsupported_selection = converters.ConverterSystemSelection(
        converters.ThreePhaseDiodeBridge,
        core.AverageValue;
        phase_count=3,
        pulse_count=6,
    )
    unsupported = converters.ConverterSystemSpecification(
        :unsupported_average_diode_bridge,
        unsupported_selection,
        ports,
        ("polyphase_diode_bridge",),
        (),
        (),
        converters.ConverterRatedBases(690.0, 100.0, 100_000.0, 50.0),
        converters.ConverterTimingParameters(
            fixed_step_s=1.0e-6,
            control_period_s=100.0e-6,
            firing_frequency_hz=300.0,
        ),
        converters.ConverterModulationParameters(kind=converters.NaturalDiodeCommutation),
        (provenance,),
        validity,
    )
    refusal = converters.converter_system_readiness(unsupported)
    @test !converters.converter_system_is_ready(refusal)
    @test refusal.reason == :unsupported_fidelity

    average_boost = converters.ConverterSystemSpecification(
        :average_boost_contract,
        converters.ConverterSystemSelection(
            converters.BoostChopper,
            core.AverageValue;
            phase_count=1,
        ),
        ports,
        ("step_up_chopper_topology",),
        (),
        ("input_inductor_output_capacitor",),
        converters.ConverterRatedBases(800.0, 100.0, 80_000.0, 10_000.0),
        specification.timing,
        specification.modulation,
        (provenance,),
        validity,
    )
    average_boost_readiness = converters.converter_system_readiness(average_boost)
    @test converters.converter_system_is_ready(average_boost_readiness)
    @test average_boost_readiness.reason == :ready
end

@testset "converter-system block and restart event contracts" begin
    converters = AIMORA.ConverterSystems
    bridges = AIMORA.BridgeTopologies
    core = AIMORA.StudyCore
    topology = bridges.step_down_chopper_topology(1, 2, 0)
    provenance = core.ParameterProvenance(
        "AIMORA converter event contract",
        "SI",
        "identity",
        "synthetic deterministic event calendar",
        "fixed-step switching chopper block and restart domain",
        core.NumericalPolicyParameter,
    )
    function event_commands(; block_time_s=0.2e-3, restart_time_s=0.8e-3,
        target_valve_indices=(1,))
        return (
            converters.ConverterSystemEventCommand(
                :controlled_path_block,
                converters.ConverterBlockEvent,
                block_time_s;
                target_valve_indices,
            ),
            converters.ConverterSystemEventCommand(
                :controlled_path_restart,
                converters.ConverterRestartEvent,
                restart_time_s;
                target_valve_indices,
                reference_id=:controlled_path_block,
            ),
        )
    end
    function event_specification(commands)
        return converters.ConverterSystemSpecification(
            :switching_chopper_event_contract,
            converters.ConverterSystemSelection(
                converters.BuckChopper,
                core.SwitchingStateEquivalent;
                phase_count=1,
            ),
            (
                converters.ConverterPortDefinition(
                    :input_dc,
                    converters.DirectCurrentPort,
                    (1, 0);
                    voltage_orientation="positive from input positive to input negative",
                ),
                converters.ConverterPortDefinition(
                    :output_dc,
                    converters.DirectCurrentPort,
                    (3, 0);
                    voltage_orientation="positive from output positive to output negative",
                ),
            ),
            (bridges.bridge_topology_signature(topology),),
            (),
            ("physical_series_inductor_and_output_capacitor",),
            converters.ConverterRatedBases(100.0, 10.0, 1_000.0, 1_000.0),
            converters.ConverterTimingParameters(
                fixed_step_s=10.0e-6,
                control_period_s=100.0e-6,
                carrier_frequency_hz=1_000.0,
            ),
            converters.ConverterModulationParameters(
                kind=converters.CarrierSinusoidalPulseWidthModulation,
                duty=0.5,
            ),
            (provenance,),
            core.ModelValidityDomain(
                :switching_chopper_event_contract;
                description="fixed-step switching chopper block and restart contract",
            );
            event_commands=commands,
        )
    end

    commands = event_commands()
    specification = event_specification(commands)
    @test converters.converter_system_is_ready(
        converters.converter_system_readiness(specification),
    )
    @test converters.validate_converter_system_event_calendar(
        specification;
        start_time_s=0.0,
        stop_time_s=2.0e-3,
        allowed_target_valve_indices=(1,),
    ) === specification
    @test isempty(converters.converter_system_commands_at_boundary(specification, 0.1e-3))
    @test converters.converter_system_commands_at_boundary(specification, 0.2e-3) ==
        (commands[1],)
    @test !converters.converter_system_is_blocked(specification, 0.19e-3, 1)
    @test converters.converter_system_is_blocked(specification, 0.2e-3, 1)
    @test converters.converter_system_is_blocked(specification, 0.79e-3)
    @test !converters.converter_system_is_blocked(specification, 0.8e-3, 1)
    @test_throws ArgumentError event_specification((commands[1],))
    duplicate_restart = converters.ConverterSystemEventCommand(
        :duplicate_controlled_path_restart,
        converters.ConverterRestartEvent,
        0.9e-3;
        target_valve_indices=(1,),
        reference_id=:controlled_path_block,
    )
    @test_throws ArgumentError event_specification((commands..., duplicate_restart))
    off_calendar = event_specification(event_commands(
        block_time_s=0.205e-3,
        restart_time_s=0.805e-3,
    ))
    @test_throws ArgumentError converters.validate_converter_system_event_calendar(
        off_calendar;
        start_time_s=0.0,
        stop_time_s=2.0e-3,
        allowed_target_valve_indices=(1,),
    )
    wrong_target = event_specification(event_commands(target_valve_indices=(2,)))
    @test_throws ArgumentError converters.validate_converter_system_event_calendar(
        wrong_target;
        start_time_s=0.0,
        stop_time_s=2.0e-3,
        allowed_target_valve_indices=(1,),
    )
    outside_horizon = event_specification(event_commands(
        block_time_s=2.1e-3,
        restart_time_s=2.2e-3,
    ))
    @test_throws ArgumentError converters.validate_converter_system_event_calendar(
        outside_horizon;
        start_time_s=0.0,
        stop_time_s=2.0e-3,
        allowed_target_valve_indices=(1,),
    )
    overlapping_commands = (
        commands[1],
        converters.ConverterSystemEventCommand(
            :overlapping_controlled_path_block,
            converters.ConverterBlockEvent,
            0.3e-3;
            target_valve_indices=(1,),
        ),
        converters.ConverterSystemEventCommand(
            :overlapping_controlled_path_restart,
            converters.ConverterRestartEvent,
            0.7e-3;
            target_valve_indices=(1,),
            reference_id=:overlapping_controlled_path_block,
        ),
        commands[2],
    )
    overlapping = event_specification(overlapping_commands)
    @test_throws ArgumentError converters.validate_converter_system_event_calendar(
        overlapping;
        start_time_s=0.0,
        stop_time_s=2.0e-3,
        allowed_target_valve_indices=(1,),
    )
end

@testset "line-commutated cycloconverter executable contracts" begin
    converters = AIMORA.ConverterSystems
    bridges = AIMORA.BridgeTopologies
    core = AIMORA.StudyCore
    topology = bridges.cycloconverter_topology((1, 2, 3), (4,), 5)
    provenance = core.ParameterProvenance(
        "AIMORA cycloconverter executable contract",
        "SI",
        "identity",
        "synthetic deterministic data",
        "noncirculating one-output fixed-step cycloconverter domain",
        core.PhysicalModelParameter,
    )
    specification = converters.ConverterSystemSpecification(
        :line_commutated_cycloconverter_contract,
        converters.ConverterSystemSelection(
            converters.LineCommutatedCycloconverter,
            core.SwitchingStateEquivalent;
            phase_count=1,
        ),
        (
            converters.ConverterPortDefinition(
                :input_ac,
                converters.AlternatingCurrentPort,
                (1, 2, 3);
                voltage_orientation="input phase order positive into converter",
            ),
            converters.ConverterPortDefinition(
                :output_ac,
                converters.AlternatingCurrentPort,
                (4, 5);
                voltage_orientation="output phase followed by neutral",
            ),
        ),
        (bridges.bridge_topology_signature(topology),),
        (),
        ("finite_three_phase_source_impedance", "series_rl_output_load"),
        converters.ConverterRatedBases(100.0, 20.0, 2_000.0, 20.0),
        converters.ConverterTimingParameters(
            fixed_step_s=10.0e-6,
            control_period_s=100.0e-6,
            firing_frequency_hz=50.0,
        ),
        converters.ConverterModulationParameters(
            kind=converters.CycloconverterFiringSynthesis,
            modulation_index=0.8,
            phase_shift_rad=0.02,
        ),
        (provenance,),
        core.ModelValidityDomain(
            :line_commutated_cycloconverter_contract;
            description="noncirculating one-output fixed-step cycloconverter contract",
        ),
    )
    study = converters.SwitchingCycloconverterStudy(
        specification;
        topology,
        input_phase_voltage_peak_v=100.0,
        input_frequency_hz=50.0,
        source_resistance_ohm=0.1,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        stop_time_s=1.0e-3,
    )
    @test converters.converter_system_is_ready(
        converters.converter_system_readiness(specification),
    )
    @test study.specification.rated_bases.frequency_hz /
        study.input_frequency_hz == 0.4
    prepared = AIMORA.prepare_converter_system(study)
    if AIMORA.solver_available()
        @test !(prepared isa AIMORA.SolverUnavailableResult)
    else
        @test prepared isa AIMORA.SolverUnavailableResult
    end
    @test_throws ArgumentError converters.SwitchingCycloconverterStudy(
        specification;
        topology,
        input_phase_voltage_peak_v=100.0,
        input_frequency_hz=50.0,
        source_resistance_ohm=0.1,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        circulating_current=true,
        stop_time_s=1.0e-3,
    )
end

@testset "average converter-application public contracts" begin
    converters = AIMORA.ConverterSystems
    bridges = AIMORA.BridgeTopologies
    core = AIMORA.StudyCore
    provenance = core.ParameterProvenance(
        "AIMORA generic converter-application contract",
        "SI",
        "identity",
        "synthetic deterministic data",
        "average-value staged application contract",
        core.PhysicalModelParameter,
    )
    grid = bridges.two_level_bridge_topology([2, 3, 4], 1, 0)
    applications = (
        converters.ShuntActiveHarmonicFilter,
        converters.DynamicVoltageRestorer,
        converters.DoubleConversionUninterruptiblePowerSupply,
        converters.ConductiveElectricVehicleCharger,
        converters.ThreeStageSolidStateTransformer,
    )
    for application in applications
        topologies, parameters, initial_state = if application ===
            converters.ShuntActiveHarmonicFilter
            (
                (grid,),
                converters.ShuntActiveFilterParameters(
                    grid_phase_voltage_peak_v=100.0,
                    grid_frequency_hz=50.0,
                    load_fundamental_current_peak_a=5.0,
                    load_fifth_harmonic_current_peak_a=1.0,
                    load_seventh_harmonic_current_peak_a=0.5,
                    filter_resistance_ohm=0.1,
                    filter_inductance_h=1.0e-3,
                    dc_link_capacitance_f=1.0e-3,
                    dc_link_reference_voltage_v=300.0,
                    current_control_resistance_ohm=10.0,
                    current_control_integral_ohm_per_s=1_000.0,
                    dc_voltage_control_a_per_v=0.01,
                    reference_loss_time_s=0.4e-3,
                    reference_restore_time_s=0.7e-3,
                ),
                converters.ShuntActiveFilterInitialState((0.0, 0.0, 0.0), 300.0),
            )
        elseif application === converters.DynamicVoltageRestorer
            (
                (grid,),
                converters.DynamicVoltageRestorerParameters(
                    source_phase_voltage_peak_v=100.0,
                    source_frequency_hz=50.0,
                    sag_start_time_s=0.2e-3,
                    sag_stop_time_s=0.8e-3,
                    sag_retained_voltage_pu=0.0,
                    bypass_start_time_s=0.4e-3,
                    bypass_stop_time_s=0.6e-3,
                    target_phase_voltage_peak_v=100.0,
                    load_resistance_ohm=10.0,
                    load_inductance_h=1.0e-3,
                    dc_link_capacitance_f=1.0e-3,
                    dc_link_reference_voltage_v=300.0,
                    injection_control_gain=1.0,
                    series_transformer_ratio=1.0,
                ),
                converters.DynamicVoltageRestorerInitialState((0.0, 0.0, 0.0), 300.0),
            )
        elseif application === converters.DoubleConversionUninterruptiblePowerSupply
            (
                (grid, bridges.two_level_bridge_topology([5, 6, 7], 1, 0)),
                converters.UninterruptiblePowerSupplyParameters(
                    source_phase_voltage_peak_v=100.0,
                    source_frequency_hz=50.0,
                    input_filter_resistance_ohm=0.1,
                    input_filter_inductance_h=1.0e-3,
                    dc_link_capacitance_f=1.0e-3,
                    dc_link_reference_voltage_v=300.0,
                    output_phase_voltage_peak_v=80.0,
                    output_frequency_hz=50.0,
                    load_resistance_ohm=10.0,
                    load_inductance_h=1.0e-3,
                    input_current_control_bandwidth_rad_s=1_000.0,
                    dc_voltage_control_w_per_v=10.0,
                    source_loss_time_s=0.2e-3,
                    bypass_transfer_time_s=0.3e-3,
                    source_recovery_time_s=0.7e-3,
                    double_conversion_restore_time_s=0.8e-3,
                    bypass_source_phase_voltage_peak_v=80.0,
                    bypass_source_frequency_hz=50.0,
                ),
                converters.UninterruptiblePowerSupplyInitialState(
                    (0.0, 0.0, 0.0), 300.0, (0.0, 0.0, 0.0),
                ),
            )
        elseif application === converters.ConductiveElectricVehicleCharger
            (
                (grid, bridges.step_down_chopper_topology(1, 5, 0)),
                converters.ConductiveElectricVehicleChargerParameters(
                    source_phase_voltage_peak_v=100.0,
                    source_frequency_hz=50.0,
                    input_filter_resistance_ohm=0.1,
                    input_filter_inductance_h=1.0e-3,
                    dc_link_capacitance_f=1.0e-3,
                    dc_link_reference_voltage_v=300.0,
                    output_inductance_h=1.0e-3,
                    output_capacitance_f=1.0e-3,
                    output_load_resistance_ohm=20.0,
                    output_voltage_reference_v=100.0,
                    input_current_control_bandwidth_rad_s=1_000.0,
                    dc_voltage_control_w_per_v=10.0,
                    output_voltage_control_a_per_v=0.1,
                    load_step_time_s=0.2e-3,
                    stepped_output_load_resistance_ohm=15.0,
                    output_short_start_time_s=0.5e-3,
                    output_short_clear_time_s=0.6e-3,
                    output_short_resistance_ohm=1.0,
                ),
                converters.ConductiveElectricVehicleChargerInitialState(
                    (0.0, 0.0, 0.0), 300.0, 5.0, 100.0,
                ),
            )
        else
            (
                (
                    grid,
                    bridges.full_bridge_topology(5, 6, 1, 0),
                    bridges.full_bridge_topology(8, 9, 7, 13),
                    bridges.two_level_bridge_topology([10, 11, 12], 7, 13),
                ),
                converters.SolidStateTransformerParameters(
                    source_phase_voltage_peak_v=100.0,
                    source_frequency_hz=50.0,
                    input_filter_resistance_ohm=0.1,
                    input_filter_inductance_h=1.0e-3,
                    primary_dc_link_capacitance_f=1.0e-3,
                    primary_dc_link_reference_voltage_v=400.0,
                    transformer_ratio=1.0,
                    transformer_resistance_ohm=0.1,
                    transformer_leakage_inductance_h=1.0e-3,
                    secondary_dc_link_capacitance_f=1.0e-3,
                    secondary_dc_link_reference_voltage_v=300.0,
                    output_phase_voltage_peak_v=80.0,
                    output_frequency_hz=50.0,
                    load_resistance_ohm=10.0,
                    load_inductance_h=1.0e-3,
                    input_current_control_bandwidth_rad_s=1_000.0,
                    primary_dc_voltage_control_w_per_v=10.0,
                    secondary_dc_voltage_control_v_per_v=1.0,
                    transformer_side_fault_start_time_s=0.4e-3,
                    transformer_side_fault_clear_time_s=0.7e-3,
                    transformer_side_fault_resistance_ohm=2.0,
                ),
                converters.SolidStateTransformerInitialState(
                    (0.0, 0.0, 0.0), 400.0, 1.0, 300.0,
                    (0.0, 0.0, 0.0),
                ),
            )
        end
        specification = converters.ConverterSystemSpecification(
            Symbol(:average_application_, lowercase(string(application))),
            converters.ConverterSystemSelection(
                converters.ThreePhaseTwoLevelBridge,
                core.AverageValue;
                application,
                phase_count=3,
            ),
            (
                converters.ConverterPortDefinition(
                    :input,
                    converters.AlternatingCurrentPort,
                    (2, 3, 4, 0);
                    voltage_orientation="phase-to-neutral input voltage",
                ),
                converters.ConverterPortDefinition(
                    :output,
                    converters.AlternatingCurrentPort,
                    (5, 6, 7, 0);
                    voltage_orientation="phase-to-neutral output voltage",
                ),
            ),
            Tuple(bridges.bridge_topology_signature.(topologies)),
            (),
            ("physical_passive_graph", "declared_transformer_boundary"),
            converters.ConverterRatedBases(300.0, 20.0, 6_000.0, 50.0),
            converters.ConverterTimingParameters(
                fixed_step_s=10.0e-6,
                control_period_s=50.0e-6,
                carrier_frequency_hz=10_000.0,
            ),
            converters.ConverterModulationParameters(
                kind=converters.SpaceVectorPulseWidthModulation,
                modulation_index=0.9,
            ),
            (provenance,),
            core.ModelValidityDomain(
                Symbol(:average_application_, lowercase(string(application)));
                description="generic average converter-application contract",
            ),
        )
        study = converters.AverageConverterApplicationStudy(
            specification;
            topologies,
            parameters,
            initial_state,
            stop_time_s=1.0e-3,
        )
        @test converters.converter_system_is_ready(
            converters.converter_system_readiness(specification),
        )
        @test study.specification.selection.application === application
        prepared = AIMORA.prepare_converter_system(study)
        if AIMORA.solver_available()
            @test !(prepared isa AIMORA.SolverUnavailableResult)
        else
            @test prepared isa AIMORA.SolverUnavailableResult
        end
        @test converters.converter_supported_fidelities(
            converters.ThreePhaseTwoLevelBridge,
            application,
        ) == (core.AverageValue,)
        if application === converters.DynamicVoltageRestorer
            @test_throws ArgumentError converters.AverageConverterApplicationStudy(
                specification;
                topologies,
                parameters,
                initial_state,
                start_time_s=5.0e-6,
                stop_time_s=1.005e-3,
            )
        elseif application === converters.DoubleConversionUninterruptiblePowerSupply
            @test parameters.source_loss_time_s < parameters.bypass_transfer_time_s <
                parameters.source_recovery_time_s <
                parameters.double_conversion_restore_time_s
            @test parameters.bypass_source_phase_voltage_peak_v ==
                parameters.output_phase_voltage_peak_v
        elseif application === converters.ShuntActiveHarmonicFilter
            @test parameters.reference_loss_time_s <
                parameters.reference_restore_time_s
        end
    end
end

@testset "matrix-converter executable contracts" begin
    converters = AIMORA.ConverterSystems
    bridges = AIMORA.BridgeTopologies
    core = AIMORA.StudyCore
    topology = bridges.matrix_converter_topology((1, 2, 3), (4, 5, 6))
    provenance = core.ParameterProvenance(
        "AIMORA matrix-converter executable contract",
        "SI",
        "identity",
        "synthetic deterministic data",
        "three-phase 3x3 fixed-step matrix-converter domain",
        core.PhysicalModelParameter,
    )
    specification = converters.ConverterSystemSpecification(
        :three_phase_matrix_converter_contract,
        converters.ConverterSystemSelection(
            converters.ThreePhaseMatrixConverter,
            core.SwitchingStateEquivalent;
            phase_count=3,
        ),
        (
            converters.ConverterPortDefinition(
                :input_ac,
                converters.AlternatingCurrentPort,
                (1, 2, 3);
                voltage_orientation="input phase order positive into converter",
            ),
            converters.ConverterPortDefinition(
                :output_ac,
                converters.AlternatingCurrentPort,
                (4, 5, 6);
                voltage_orientation="output phase order positive from converter",
            ),
        ),
        (bridges.bridge_topology_signature(topology),),
        (),
        ("finite_three_phase_source_impedance", "balanced_series_rl_load"),
        converters.ConverterRatedBases(100.0, 20.0, 2_000.0, 30.0),
        converters.ConverterTimingParameters(
            fixed_step_s=0.25e-6,
            control_period_s=100.0e-6,
            carrier_frequency_hz=10_000.0,
            dead_time_s=0.25e-6,
        ),
        converters.ConverterModulationParameters(
            kind=converters.MatrixSpaceVectorModulation,
            modulation_index=0.75,
        ),
        (provenance,),
        core.ModelValidityDomain(
            :three_phase_matrix_converter_contract;
            description="three-phase 3x3 fixed-step matrix-converter contract",
        ),
    )
    study = converters.SwitchingMatrixConverterStudy(
        specification;
        topology,
        input_phase_voltage_peak_v=100.0,
        input_frequency_hz=50.0,
        source_resistance_ohm=0.1,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        stop_time_s=20.0e-6,
    )
    @test converters.converter_system_is_ready(
        converters.converter_system_readiness(specification),
    )
    modulation = converters.converter_matrix_space_vector_state(
        (0.0, -sqrt(3.0) * 50.0, sqrt(3.0) * 50.0),
        0.0,
        30.0,
        10_000.0,
        0.75,
    )
    @test all(sum(modulation.connection; dims=2) .== 1)
    @test count(modulation.requested_valve_state) == 6
    @test length(modulation.deterministic_signature_sha256) == 64
    prepared = AIMORA.prepare_converter_system(study)
    if AIMORA.solver_available()
        @test !(prepared isa AIMORA.SolverUnavailableResult)
    else
        @test prepared isa AIMORA.SolverUnavailableResult
    end
end

@testset "direct AC/AC topology state contracts" begin
    bridges = AIMORA.BridgeTopologies
    matrix_topology = bridges.matrix_converter_topology((1, 2, 3), (4, 5, 6))
    @test matrix_topology.family == :matrix_converter
    @test length(matrix_topology.valve_positions) == 18
    @test size(only(matrix_topology.state_groups).admitted_states) == (18, 27)
    first_matrix_state = only(matrix_topology.state_groups).admitted_states[:, 1]
    @test bridges.bridge_topology_state_is_allowed(matrix_topology, first_matrix_state)
    invalid_matrix_state = copy(first_matrix_state)
    invalid_matrix_state[1] = false
    @test !bridges.bridge_topology_state_is_allowed(matrix_topology, invalid_matrix_state)

    noncirculating = bridges.cycloconverter_topology((1, 2, 3), (4,), 5)
    @test noncirculating.family == :cycloconverter
    @test length(noncirculating.valve_positions) == 12
    @test size(only(noncirculating.state_groups).admitted_states) == (12, 13)
    circulating = bridges.cycloconverter_topology(
        (1, 2, 3),
        (4, 5, 6),
        7;
        circulating_current=true,
    )
    @test length(circulating.valve_positions) == 36
    @test all(group -> size(group.admitted_states) == (12, 49), circulating.state_groups)
end

@testset "inverting buck-boost topology contract" begin
    bridges = AIMORA.BridgeTopologies
    topology = bridges.inverting_buck_boost_topology(1, 2, 3, 0)
    @test topology.family === :inverting_buck_boost
    @test getfield.(topology.nodes, :name) ==
        [:input_positive, :switching, :output_negative, :reference]
    @test [(valve.from_node, valve.to_node) for valve in topology.valve_positions] ==
        [(1, 2), (3, 2)]
    @test bridges.bridge_topology_state_is_allowed(topology, Bool[true, false])
    @test bridges.bridge_topology_state_is_allowed(topology, Bool[false, true])
end

@testset "converter-system governing equation checks" begin
    converters = AIMORA.ConverterSystems
    dab_power = converters.dual_active_bridge_power_w(
        800.0,
        400.0,
        1.0,
        pi / 6,
        2pi * 10_000.0,
        50.0e-6,
    )
    @test dab_power > 0.0
    @test converters.dual_active_bridge_power_w(
        800.0,
        400.0,
        1.0,
        -pi / 6,
        2pi * 10_000.0,
        50.0e-6,
    ) ≈ -dab_power
    @test converters.controlled_rectifier_average_voltage_v(
        converters.SinglePhaseThyristorBridge,
        100.0,
        0.0,
    ) ≈ 200.0 / pi
    @test converters.controlled_rectifier_average_voltage_v(
        converters.ThreePhaseThyristorBridge,
        400.0,
        0.0,
    ) ≈ 1200.0 * sqrt(2.0) / pi
    @test converters.interleaved_carrier_phases_rad(4) ≈ [0.0, pi / 2, pi, 3pi / 2]

    connection = Bool[
        1 0 0
        0 0 1
        0 1 0
    ]
    mapped = converters.matrix_converter_terminal_map(
        connection,
        [100.0, -40.0, -60.0],
        [2.0, -1.0, -1.0],
    )
    @test mapped.output_voltage_v == [100.0, -60.0, -40.0]
    @test mapped.input_current_a == [2.0, -1.0, -1.0]
    @test dot(mapped.output_voltage_v, [2.0, -1.0, -1.0]) ≈
        dot([100.0, -40.0, -60.0], mapped.input_current_a)
    @test converters.converter_stage_energy_residual_w(100.0, 80.0, 5.0, 15.0) == 0.0
    @test_throws ArgumentError converters.matrix_converter_terminal_map(
        trues(3, 3),
        zeros(3),
        zeros(3),
    )
end

@testset "converter modulation and safe commutation contracts" begin
    converters = AIMORA.ConverterSystems
    sinusoidal = converters.converter_two_level_duties(
        (100.0, -50.0, -50.0),
        400.0,
        converters.CarrierSinusoidalPulseWidthModulation,
    )
    space_vector = converters.converter_two_level_duties(
        (100.0, -50.0, -50.0),
        400.0,
        converters.SpaceVectorPulseWidthModulation,
    )
    @test sinusoidal == (0.75, 0.375, 0.375)
    @test space_vector == (0.6875, 0.3125, 0.3125)
    @test converters.converter_triangular_carrier(0.0, 1_000.0) == 0.0
    @test converters.converter_triangular_carrier(0.5e-3, 1_000.0) == 1.0
    pwm = converters.converter_pwm_gate_state(
        [0.25, 0.5, 0.75],
        0.25e-3,
        1_000.0,
        converters.CarrierSinusoidalPulseWidthModulation,
    )
    @test pwm.carrier_values == fill(0.5, 3)
    @test pwm.requested_valve_state == BitVector((false, true, true, false, true, false))
    @test length(pwm.deterministic_signature_sha256) == 64
    tuple_pwm = converters.converter_pwm_gate_state(
        (0.25, 0.5, 0.75),
        0.25e-3,
        1_000.0,
        converters.CarrierSinusoidalPulseWidthModulation,
    )
    @test tuple_pwm.duties == [0.25, 0.5, 0.75]
    three_level = converters.converter_three_level_gate_state(
        (0.75, 0.0, -0.75),
        0.25e-3,
        1_000.0,
    )
    @test three_level.carrier_values == fill(0.5, 3)
    @test three_level.requested_level == Int8[1, 0, -1]
    @test three_level.requested_valve_state == BitVector((
        true, true, false, false,
        false, true, true, false,
        false, false, true, true,
    ))
    @test length(three_level.deterministic_signature_sha256) == 64
    flying_capacitor = converters.converter_flying_capacitor_gate_state(
        (0.75, 0.0, -0.75),
        (0.0, 0.0, 0.0),
        (1.0, -1.0, 0.0),
        0.25e-3,
        1_000.0,
    )
    repeated_flying_capacitor = converters.converter_flying_capacitor_gate_state(
        (0.75, 0.0, -0.75),
        (0.0, 0.0, 0.0),
        (1.0, -1.0, 0.0),
        0.25e-3,
        1_000.0,
    )
    @test flying_capacitor.requested_level == Int8[1, 0, -1]
    @test flying_capacitor.requested_valve_state == BitVector((
        true, true, false, false,
        false, true, false, true,
        false, false, true, true,
    ))
    @test repeated_flying_capacitor.requested_valve_state ==
        flying_capacitor.requested_valve_state
    @test repeated_flying_capacitor.deterministic_signature_sha256 ==
        flying_capacitor.deterministic_signature_sha256
    @test length(flying_capacitor.deterministic_signature_sha256) == 64
    t_type_level = converters.converter_three_level_gate_state(
        (0.75, 0.0, -0.75),
        0.25e-3,
        1_000.0;
        leg_topology=:t_type,
    )
    @test t_type_level.requested_level == Int8[1, 0, -1]
    @test t_type_level.requested_valve_state == BitVector((
        true, false, false, false,
        false, true, true, false,
        false, false, false, true,
    ))

    she = converters.selective_harmonic_elimination_residual(
        [0.2, 0.5, 0.8],
        0.9,
        [5, 7],
    )
    @test length(she) == 3
    @test all(isfinite, she)
    nearest = converters.nearest_level_cell_state(0.62, 5)
    @test nearest.requested_level == 3
    @test nearest.cell_state == Int8[1, 1, 1, 0, 0]
    @test nearest.normalized_voltage == 0.6
    dab = converters.dual_active_bridge_gate_state(0.0, 10_000.0, pi / 2)
    @test count(dab.primary_state) == 2
    @test count(dab.secondary_state) == 2
    @test length(dab.requested_valve_state) == 8
    half_cycle = converters.dual_active_bridge_gate_state(
        0.5 / 10_000.0,
        10_000.0,
        pi / 2,
    )
    arithmetically_equivalent_half_cycle = converters.dual_active_bridge_gate_state(
        100 * 0.5e-6,
        10_000.0,
        pi / 2,
    )
    @test half_cycle.primary_state == BitVector((false, true, true, false))
    @test arithmetically_equivalent_half_cycle.primary_state == half_cycle.primary_state

    for current in (-10.0, 10.0)
        commutation = converters.matrix_safe_commutation_sequence(1, 3, current)
        @test length(commutation.stages) == 5
        @test commutation.stages[1][1:2] == trues(2)
        @test commutation.stages[end][5:6] == trues(2)
        @test all(stage -> any(stage), commutation.stages)
        @test all(stage -> begin
            active_inputs = [any(stage[(2 * input - 1):(2 * input)]) for input in 1:3]
            count(active_inputs) <= 2
        end, commutation.stages)
    end
    @test_throws ArgumentError converters.matrix_safe_commutation_sequence(1, 1, 10.0)
    @test_throws ArgumentError converters.matrix_safe_commutation_sequence(1, 2, 0.0)
end

@testset "cascaded H-bridge modulation and public solver-free contracts" begin
    converters = AIMORA.ConverterSystems
    core = AIMORA.StudyCore
    bridges = AIMORA.BridgeTopologies
    cell_count = 2

    voltage_error_v = [0.2 -0.2; -0.1 0.1; 0.05 -0.05]
    phase_current_a = [2.0, -1.0, -1.0]
    normalized_reference = [0.7, -0.4, 0.0]
    for kind in (
        converters.CarrierSinusoidalPulseWidthModulation,
        converters.PhaseShiftedCarrierPulseWidthModulation,
        converters.SelectiveHarmonicElimination,
        converters.NearestLevelModulation,
    )
        angles = kind === converters.SelectiveHarmonicElimination ?
            (0.2, 0.5) : ()
        first = converters.converter_cascaded_h_bridge_gate_state(
            normalized_reference,
            voltage_error_v,
            phase_current_a,
            7.25e-6,
            50_000.0,
            kind;
            modulation_index=0.8,
            selective_harmonic_angles_rad=angles,
        )
        repeated = converters.converter_cascaded_h_bridge_gate_state(
            normalized_reference,
            voltage_error_v,
            phase_current_a,
            7.25e-6,
            50_000.0,
            kind;
            modulation_index=0.8,
            selective_harmonic_angles_rad=angles,
        )
        @test first.requested_level == vec(sum(first.requested_cell_state; dims=2))
        @test all(value -> value in Int8[-1, 0, 1], first.requested_cell_state)
        @test all(cell -> count(identity,
            first.requested_valve_state[(4 * (cell - 1) + 1):(4 * cell)],
        ) == 2, 1:(3 * cell_count))
        @test repeated.requested_cell_state == first.requested_cell_state
        @test repeated.requested_valve_state == first.requested_valve_state
        @test repeated.deterministic_signature_sha256 ==
            first.deterministic_signature_sha256
        @test length(first.deterministic_signature_sha256) == 64
    end
    @test_throws ArgumentError converters.converter_cascaded_h_bridge_gate_state(
        [0.9, 0.0, 0.0],
        voltage_error_v,
        phase_current_a,
        0.0,
        50_000.0,
        converters.NearestLevelModulation;
        modulation_index=0.8,
    )

    next_node = Ref(4)
    topologies = ntuple(3) do phase
        series_nodes = [next_node[]]
        next_node[] += 1
        dc_nodes = Tuple{Int,Int}[]
        for _cell in 1:cell_count
            push!(dc_nodes, (next_node[], next_node[] + 1))
            next_node[] += 2
        end
        bridges.cascaded_h_bridge_phase_topology(
            phase,
            0,
            dc_nodes,
            series_nodes,
        )
    end
    ports = converters.ConverterPortDefinition[
        converters.ConverterPortDefinition(
            :output_ac,
            converters.AlternatingCurrentPort,
            (1, 2, 3);
            voltage_orientation="phase-to-floating-neutral positive sequence",
        ),
    ]
    for phase in 1:3
        nodes = Dict(node.name => node.node for node in topologies[phase].nodes)
        for cell in 1:cell_count
            push!(ports, converters.ConverterPortDefinition(
                Symbol(:phase_, phase, :_cell_, cell, :_dc),
                converters.IsolatedDirectCurrentPort,
                (
                    nodes[Symbol(:cell_, cell, :_dc_positive)],
                    nodes[Symbol(:cell_, cell, :_dc_negative)],
                );
                voltage_orientation="positive from isolated cell positive to negative rail",
            ))
        end
    end
    provenance = core.ParameterProvenance(
        "AIMORA generic cascaded H-bridge inverter",
        "SI",
        "identity",
        "synthetic bounded parameters",
        "balanced three-wire isolated-cell fixed-step domain",
        core.PhysicalModelParameter,
    )
    specification = converters.ConverterSystemSpecification(
        :cascaded_h_bridge_public_contract,
        converters.ConverterSystemSelection(
            converters.CascadedHBridge,
            core.SwitchingStateEquivalent;
            phase_count=3,
            cell_count,
        ),
        Tuple(ports),
        Tuple(bridges.bridge_topology_signature.(topologies)),
        (),
        ("isolated_cell_dc_capacitors", "physical_balanced_three_wire_series_rl_load"),
        converters.ConverterRatedBases(100.0, 20.0, 2_000.0, 50.0),
        converters.ConverterTimingParameters(
            fixed_step_s=0.25e-6,
            control_period_s=20.0e-6,
            carrier_frequency_hz=50_000.0,
            dead_time_s=0.25e-6,
        ),
        converters.ConverterModulationParameters(
            kind=converters.PhaseShiftedCarrierPulseWidthModulation,
            modulation_index=0.8,
        ),
        (provenance,),
        core.ModelValidityDomain(
            :cascaded_h_bridge_public_contract;
            description="balanced three-wire isolated-cell cascaded H-bridge contract",
        ),
    )
    readiness = converters.converter_system_readiness(specification)
    @test converters.converter_system_is_ready(readiness)
    @test readiness.reason === :ready
    initial_state = converters.CascadedHBridgeInitialState(
        (0.0, 0.0, 0.0);
        cell_dc_voltage_v=fill(50.0, 3, cell_count),
    )
    study = converters.SwitchingCascadedHBridgeStudy(
        specification;
        topologies,
        cell_dc_capacitance_f=2.0e-3,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        initial_state,
        stop_time_s=0.1e-3,
    )
    @test study.specification.selection.family === converters.CascadedHBridge
    @test study.detailed_semiconductor === nothing
    prepared = AIMORA.prepare_converter_system(study)
    if AIMORA.solver_available()
        @test !(prepared isa AIMORA.SolverUnavailableResult)
    else
        @test prepared isa AIMORA.SolverUnavailableResult
        @test prepared.required_capability === :extended_converter_systems
    end
    @test_throws DimensionMismatch converters.CascadedHBridgeInitialState(
        (0.0, 0.0, 0.0);
        cell_dc_voltage_v=fill(50.0, 3, 1),
    )
    @test_throws ArgumentError converters.SwitchingCascadedHBridgeStudy(
        specification;
        topologies,
        cell_dc_capacitance_f=2.0e-3,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        initial_state,
        detailed_semiconductor=converters.DetailedChopperSemiconductorParameters(
            ; provenance,
        ),
        stop_time_s=0.1e-3,
    )
end

@testset "flying-capacitor public study and solver-free contracts" begin
    converters = AIMORA.ConverterSystems
    core = AIMORA.StudyCore
    bridges = AIMORA.BridgeTopologies
    topologies = (
        bridges.flying_capacitor_leg_topology(1, 2, 0, 5, 6),
        bridges.flying_capacitor_leg_topology(1, 3, 0, 7, 8),
        bridges.flying_capacitor_leg_topology(1, 4, 0, 9, 10),
    )
    provenance = core.ParameterProvenance(
        "AIMORA generic flying-capacitor inverter",
        "SI",
        "identity",
        "synthetic bounded parameters",
        "balanced three-wire fixed-step carrier-PWM domain",
        core.PhysicalModelParameter,
    )
    function specification(fidelity)
        return converters.ConverterSystemSpecification(
            :flying_capacitor_public_contract,
            converters.ConverterSystemSelection(
                converters.FlyingCapacitorBridge,
                fidelity;
                phase_count=3,
            ),
            (
                converters.ConverterPortDefinition(
                    :input_dc,
                    converters.DirectCurrentPort,
                    (1, 0);
                    voltage_orientation="positive from DC positive to DC negative",
                ),
                converters.ConverterPortDefinition(
                    :output_ac,
                    converters.AlternatingCurrentPort,
                    (2, 3, 4);
                    voltage_orientation="phase-to-floating-neutral positive sequence",
                ),
            ),
            Tuple(bridges.bridge_topology_signature.(topologies)),
            (),
            ("physical_flying_capacitors", "physical_balanced_three_wire_series_rl_load"),
            converters.ConverterRatedBases(100.0, 20.0, 2_000.0, 50.0),
            converters.ConverterTimingParameters(
                fixed_step_s=0.25e-6,
                control_period_s=20.0e-6,
                carrier_frequency_hz=50_000.0,
                dead_time_s=0.25e-6,
            ),
            converters.ConverterModulationParameters(
                kind=converters.CarrierSinusoidalPulseWidthModulation,
                duty=0.5,
                modulation_index=0.8,
            ),
            (provenance,),
            core.ModelValidityDomain(
                :flying_capacitor_public_contract;
                description="balanced three-wire flying-capacitor contract",
            ),
        )
    end
    switching_state_specification = specification(core.SwitchingStateEquivalent)
    @test converters.converter_system_is_ready(
        converters.converter_system_readiness(switching_state_specification),
    )
    initial_state = converters.FlyingCapacitorInitialState(
        (0.0, 0.0, 0.0);
        flying_capacitor_voltage_v=(50.0, 50.0, 50.0),
    )
    study = converters.SwitchingFlyingCapacitorStudy(
        switching_state_specification;
        topologies,
        input_voltage_v=100.0,
        source_resistance_ohm=0.05,
        flying_capacitance_f=2.0e-3,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        initial_state,
        stop_time_s=0.1e-3,
    )
    @test study.specification.selection.family === converters.FlyingCapacitorBridge
    @test study.detailed_semiconductor === nothing
    @test study.balance_voltage_tolerance_v == 1.0e-3
    prepared = AIMORA.prepare_converter_system(study)
    if AIMORA.solver_available()
        @test !(prepared isa AIMORA.SolverUnavailableResult)
    else
        @test prepared isa AIMORA.SolverUnavailableResult
    end
    @test_throws ArgumentError converters.SwitchingFlyingCapacitorStudy(
        specification(core.SwitchingDetailed);
        topologies,
        input_voltage_v=100.0,
        source_resistance_ohm=0.05,
        flying_capacitance_f=2.0e-3,
        load_resistance_ohm=5.0,
        load_inductance_h=2.0e-3,
        initial_state,
        stop_time_s=0.1e-3,
    )
    @test_throws ArgumentError converters.FlyingCapacitorInitialState(
        (1.0, 0.0, 0.0);
        flying_capacitor_voltage_v=(50.0, 50.0, 50.0),
    )
    @test_throws ArgumentError converters.FlyingCapacitorInitialState(
        (0.0, 0.0, 0.0);
        flying_capacitor_voltage_v=(50.0, 0.0, 50.0),
    )
end

@testset "line-commutated diode rectifier contracts" begin
    converters = AIMORA.ConverterSystems
    core = AIMORA.StudyCore
    bridges = AIMORA.BridgeTopologies
    topology = bridges.single_phase_graetz_topology((1, 2), 3, 0)
    provenance = core.ParameterProvenance(
        "AIMORA generic line-commutated diode rectifier",
        "SI",
        "identity",
        "synthetic bounded parameters",
        "fixed-step ideal natural-commutation domain",
        core.PhysicalModelParameter,
    )
    specification = converters.ConverterSystemSpecification(
        :single_phase_diode_rectifier_contract,
        converters.ConverterSystemSelection(
            converters.SinglePhaseDiodeBridge,
            core.SwitchingStateEquivalent;
            phase_count=1,
            pulse_count=2,
        ),
        (
            converters.ConverterPortDefinition(
                :input_ac,
                converters.AlternatingCurrentPort,
                (1, 2);
                voltage_orientation="positive from AC terminal one to terminal two",
            ),
            converters.ConverterPortDefinition(
                :output_dc,
                converters.DirectCurrentPort,
                (3, 0);
                voltage_orientation="positive from DC positive to DC negative",
            ),
        ),
        (bridges.bridge_topology_signature(topology),),
        (),
        ("physical_series_rl_dc_load",),
        converters.ConverterRatedBases(100.0, 20.0, 2_000.0, 1_000.0),
        converters.ConverterTimingParameters(
            fixed_step_s=0.5e-6,
            control_period_s=10.0e-6,
            firing_frequency_hz=1_000.0,
        ),
        converters.ConverterModulationParameters(
            kind=converters.NaturalDiodeCommutation,
        ),
        (provenance,),
        core.ModelValidityDomain(
            :single_phase_diode_rectifier_contract;
            description="fixed-step ideal natural-commutation contract",
        ),
    )
    study = converters.SwitchingLineCommutatedRectifierStudy(
        specification;
        topology,
        source_peak_voltage_v=100.0,
        source_frequency_hz=1_000.0,
        source_resistance_ohm=0.1,
        load_resistance_ohm=10.0,
        load_inductance_h=2.0e-3,
        initial_state=converters.LineCommutatedRectifierInitialState(),
        stop_time_s=0.1e-3,
    )
    @test study.specification.selection.family === converters.SinglePhaseDiodeBridge
    @test study.detailed_semiconductor === nothing
    @test converters.converter_system_is_ready(
        converters.converter_system_readiness(study.specification),
    )
    @test length(converters.detailed_rectifier_semiconductor_signatures(
        converters.DetailedConverterSemiconductorParameters(; provenance),
        4,
    )) == 4
    @test length(converters.detailed_rectifier_semiconductor_signatures(
        converters.DetailedConverterSemiconductorParameters(; provenance),
        24,
    )) == 24
    @test_throws ArgumentError converters.LineCommutatedRectifierInitialState(-1.0)
    prepared = AIMORA.prepare_converter_system(study)
    if AIMORA.solver_available()
        @test !(prepared isa AIMORA.SolverUnavailableResult)
    else
        @test prepared isa AIMORA.SolverUnavailableResult
    end
end

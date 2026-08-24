using AIMORA.SurgeInsulation

function surge_test_line_runtime_preparation(timestep_s::Float64)
    fitting = AIMORA.CoupledLineFitting
    runtime = AIMORA.CoupledLineRuntime
    rate_per_s = 2.0 * pi * 500.0
    direct = 0.03 .* Matrix{Float64}(I, 4, 4) .+ 0.001 .* ones(4, 4)
    normalized_residue = 0.15 .* Matrix{Float64}(I, 4, 4) .+ 0.004 .* ones(4, 4)
    model = fitting.coupled_line_rational_model(
        [rate_per_s],
        direct,
        [rate_per_s .* normalized_residue];
        port_order=[:sending_a, :sending_b, :receiving_a, :receiving_b],
        reference_impedance_ohm=[40.0, 45.0, 40.0, 45.0],
    )
    frequencies_hz = [1.0, 100.0, 1.0e3]
    return runtime._coupled_line_runtime_preparation(
        model,
        runtime.CoupledLineRuntimeSettings(; timestep_s);
        source_signature_sha256=repeat("a", 64),
        response_signature_sha256=repeat("b", 64),
        fit_signature_sha256=repeat("c", 64),
        phase_order=[:a, :b],
        frequencies_hz,
        continuous_passivity_passed=true,
        uncertainty_alternative_fit_signatures_sha256=[repeat("d", 64)],
        uncertainty_set_complete=false,
    )
end

function surge_test_transformer_preparation(timestep_s::Float64, provenance)
    transformers = AIMORA.TransformerApparatus
    source = transformers.TransformerSourceRecord(
        :surge_test_transformer,
        repeat("e", 64),
        provenance,
    )
    connection = transformers.TransformerConnectionTopology(
        node_order=[:primary_terminal, :secondary_terminal],
        coil_order=[:primary_coil, :secondary_coil],
        winding_order=[:primary_winding, :secondary_winding],
        phase_order=[:phase_a],
        coil_winding=[:primary_winding, :secondary_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=Matrix{Float64}(I, 2, 2),
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
    matrices = transformers.TransformerTerminalMatrices(
        [0.4 0.02; 0.02 0.5],
        [0.12 0.01; 0.01 0.10];
        capacitance_f=[2.0e-9 0.0; 0.0 1.0e-9],
        conductance_s=[1.0e-7 0.0; 0.0 2.0e-7],
    )
    specification = transformers.TransformerApparatusSpecification(
        :surge_test_transformer,
        transformers.LowFrequencyTerminalTier,
        connection,
        transformers.LowFrequencyTransformerModel(matrices),
        transformers.TransformerRuntimeSettings(
            timestep_s=timestep_s,
            initialization_frequency_hz=60.0,
        );
        phase_count=1,
        rated_power_va=1.0e6,
        rated_voltage_v=13.8e3,
        rated_frequency_hz=60.0,
        sources=[source],
        uncertainty="exact synthetic contract fixture",
        validity_domain="one-phase surge owner contract test",
    )
    return transformers.prepare_transformer_apparatus(specification)
end

@testset "surge and insulation public scientific contracts" begin
    double_exponential = DoubleExponentialLightningImpulse(30.0e3, 2.0e4, 2.0e6)
    @test lightning_current_a(double_exponential, -1.0e-6) == 0.0
    @test lightning_current_a(
        double_exponential,
        double_exponential.peak_time_s,
    ) ≈ 30.0e3 rtol=2.0e-15
    double_metrics = lightning_impulse_metrics(
        double_exponential;
        duration_s=500.0e-6,
        intervals=4096,
    )
    @test double_metrics.peak_current_a ≈ 30.0e3 rtol=3.0e-4
    @test double_metrics.charge_c > 0.0
    @test double_metrics.specific_energy_a2s > 0.0

    heidler = HeidlerLightningImpulse(-50.0e3, 1.0e-6, 50.0e-6, 4.0)
    @test lightning_current_a(heidler, heidler.peak_time_s) ≈ -50.0e3 rtol=2.0e-14
    sequence = LightningStrokeSequence(
        :two_stroke,
        [
            LightningStroke(:first, 0.0, double_exponential),
            LightningStroke(:second, 1.0e-3, heidler),
        ],
    )
    @test lightning_sequence_current_a(sequence, 0.0) == 0.0
    @test lightning_sequence_current_a(
        sequence,
        double_exponential.peak_time_s,
    ) ≈ 30.0e3 rtol=2.0e-15
    @test_throws ArgumentError LightningStrokeSequence(
        :unsorted,
        [
            LightningStroke(:late, 1.0, double_exponential),
            LightningStroke(:early, 0.0, heidler),
        ],
    )

    current_workspace = zeros(Float64, 2)
    jacobian_workspace = zeros(Float64, 2, 2)
    combined_arc = CombinedArcBranch(
        1,
        0;
        cassie_power_w=2.0e5,
        mayr_power_w=2.0e3,
        cassie_time_constant_s=20.0e-6,
        mayr_time_constant_s=5.0e-6,
        transition_power_w=2.0e4,
        initial_conductance_s=10.0,
        extinction_current_a=0.1,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        combined_arc,
        0.0,
        1.0e-6,
        :backward_euler,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        combined_arc,
        [100.0, 0.0],
        1.0e-6,
    )
    @test current_workspace[1] == -current_workspace[2]
    @test jacobian_workspace[1, 1] > 0.0
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        combined_arc,
        [100.0, 0.0],
        current_workspace,
        1.0e-6,
    )
    @test arc_conductance_s(combined_arc) > 0.0
    @test arc_dissipated_energy_j(combined_arc) > 0.0
    arc_snapshot = snapshot_surge_device(combined_arc)
    old_arc_conductance = arc_conductance_s(combined_arc)
    ignite_arc!(combined_arc, 20.0)
    @test arc_conductance_s(combined_arc) == 20.0
    restore_surge_device!(combined_arc, arc_snapshot)
    @test arc_conductance_s(combined_arc) == old_arc_conductance

    fault_arc = FaultArcBranch(
        2,
        0;
        cooling_power_w=1.0e4,
        time_constant_s=10.0e-6,
        arc_length_m=0.5,
        initial_conductance_s=1.0,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        fault_arc,
        0.0,
        1.0e-6,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        fault_arc,
        [50.0, 0.0],
        1.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        fault_arc,
        [50.0, 0.0],
        current_workspace,
        1.0e-6,
    )
    @test arc_dissipated_energy_j(fault_arc) > 0.0

    vacuum = VacuumInterruptionState(
        :pole_a;
        chopping_current_a=5.0,
        initial_dielectric_strength_v=1.0e3,
        dielectric_recovery_rate_v_per_s=1.0e9,
        maximum_dielectric_strength_v=100.0e3,
    )
    request_vacuum_open!(vacuum, 1.0e-3)
    values = vacuum_interruption_event_values(vacuum, 4.0, 0.0, 1.001e-3)
    @test values.chop_surface_a < 0.0
    apply_vacuum_chop!(vacuum, 1.001e-3)
    @test vacuum.state === :open
    restrike_values = vacuum_interruption_event_values(vacuum, 0.0, 200.0e3, 1.002e-3)
    @test restrike_values.restrike_surface_v > 0.0
    apply_vacuum_restrike!(vacuum, 1.002e-3)
    @test vacuum.state === :restruck

    characteristic = MetalOxideCharacteristic(
        [1.0, 10.0, 100.0],
        [1.0e3, 2.0e3, 3.0e3];
        extrapolation=:power_law,
    )
    @test metal_oxide_current_and_derivative(characteristic, 2.0e3).current_a ≈ 10.0
    @test metal_oxide_current_and_derivative(characteristic, -2.0e3).current_a ≈ -10.0
    arrester = MetalOxideArrester(
        1,
        0,
        characteristic;
        thermal_capacitance_j_per_k=100.0,
        thermal_resistance_k_per_w=10.0,
        maximum_temperature_k=500.0,
        maximum_energy_j=1.0e6,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        arrester,
        0.0,
        1.0e-6,
        :backward_euler,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        arrester,
        [2.0e3, 0.0],
        1.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        arrester,
        [2.0e3, 0.0],
        current_workspace,
        1.0e-6,
    )
    @test arrester_charge_c(arrester) > 0.0
    @test arrester_absorbed_energy_j(arrester) > 0.0
    @test arrester_temperature_k(arrester) >= 293.15
    @test !arrester_duty_margin(arrester).failed
    arrester_snapshot = snapshot_surge_device(arrester)
    arrester.charge_c = -1.0
    restore_surge_device!(arrester, arrester_snapshot)
    @test arrester_charge_c(arrester) > 0.0

    dynamic_arrester = DynamicArresterEquivalent(
        characteristic,
        characteristic;
        series_resistance_ohm=0.1,
        series_inductance_h=1.0e-6,
        filter_resistance_ohm=100.0,
        filter_inductance_h=10.0e-6,
        shunt_capacitance_f=100.0e-12,
    )
    @test dynamic_arrester.filter_inductance_h == 10.0e-6

    grounding = PositiveRealGroundingModel(
        0.01,
        [1.0e3, 1.0e5],
        [0.02, 0.03];
        maximum_frequency_hz=1.0e6,
    )
    @test real(grounding_admittance_s(grounding, 10.0e3)) > 0.0
    @test real(grounding_impedance_ohm(grounding, 10.0e3)) > 0.0

    ionizing_ground = IonizingGroundBranch(
        1;
        linear_resistance_ohm=20.0,
        electrode_radius_m=0.01,
        maximum_ionized_radius_m=1.0,
        critical_field_v_per_m=1.0e5,
        expansion_rate_m_per_v_s=1.0e-6,
        recovery_rate_per_s=1.0e3,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        ionizing_ground,
        0.0,
        1.0e-6,
        :backward_euler,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        ionizing_ground,
        [10.0e3, 0.0],
        1.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        ionizing_ground,
        [10.0e3, 0.0],
        current_workspace,
        1.0e-6,
    )
    @test ionized_radius_m(ionizing_ground) > 0.01
    @test ground_potential_rise_v(ionizing_ground) == 10.0e3
    @test grounding_dissipated_energy_j(ionizing_ground) > 0.0

    corona = DynamicCoronaBranch(
        1,
        0;
        base_capacitance_f=1.0e-9,
        incremental_capacitance_f_per_v=1.0e-15,
        onset_voltage_v=100.0e3,
        extinction_voltage_v=80.0e3,
        loss_conductance_s=1.0e-7,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        corona,
        0.0,
        1.0e-6,
        :backward_euler,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        corona,
        [120.0e3, 0.0],
        1.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        corona,
        [120.0e3, 0.0],
        current_workspace,
        1.0e-6,
    )
    @test corona.active
    @test corona_charge_c(corona) > 0.0
    @test corona_dissipated_energy_j(corona) > 0.0

    disruptive = DisruptiveEffectInsulator(
        1,
        0;
        positive_threshold_voltage_v=100.0,
        negative_threshold_voltage_v=100.0,
        positive_exponent=1.0,
        negative_exponent=1.0,
        positive_critical_effect=1.0,
        negative_critical_effect=1.0,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        disruptive,
        0.0,
        0.1,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        disruptive,
        [200.0, 0.0],
        0.1,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        disruptive,
        [200.0, 0.0],
        current_workspace,
        0.1,
    )
    @test insulation_flashover_margin(disruptive) < 0.0
    apply_insulation_flashover!(disruptive, 0.1)
    @test disruptive_effect_state(disruptive).flashed
    @test disruptive_effect_state(disruptive).polarity == 1

    leader = LeaderProgressionInsulator(
        1,
        0;
        gap_length_m=1.0,
        positive_inception_field_v_per_m=100.0,
        negative_inception_field_v_per_m=100.0,
        positive_velocity_coefficient=0.01,
        negative_velocity_coefficient=0.01,
        velocity_exponent=1.0,
    )
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        leader,
        0.0,
        0.25,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current_workspace,
        jacobian_workspace,
        leader,
        [1.0e3, 0.0],
        0.25,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        leader,
        [1.0e3, 0.0],
        current_workspace,
        0.25,
    )
    @test leader_progression_state(leader).length_m == 1.0
    apply_insulation_flashover!(leader, 0.25)
    @test leader_progression_state(leader).flashed

    section = SurgePropagationSection(
        :tower_top,
        [400.0 20.0; 20.0 350.0];
        length_m=30.0,
        propagation_speed_m_per_s=2.5e8,
    )
    waves = traveling_wave_components(section, [100.0, 50.0], [0.2, -0.1])
    reconstructed = reconstruct_terminal_state(
        section,
        waves.forward_voltage_v,
        waves.reverse_voltage_v,
    )
    @test reconstructed.voltage_v ≈ [100.0, 50.0]
    @test reconstructed.current_a ≈ [0.2, -0.1]
    matched = traveling_wave_reflection(
        section,
        section.characteristic_impedance_ohm,
        [100.0, 50.0],
    )
    @test matched.reflected_voltage_v ≈ zeros(2) atol=1.0e-14
    tower = TransmissionTowerModel(
        :two_section_tower,
        [section, section],
        [:top, :crossarm, :footing],
    )
    @test tower_total_travel_time_s(tower) == 2.0 * section.travel_time_s

    gis = GISGILSection(
        :gis_bus,
        [1.0e-4 0.0; 0.0 1.0e-4],
        [1.0e-6 0.1e-6; 0.1e-6 1.0e-6],
        [1.0e-9 0.0; 0.0 1.0e-9],
        [100.0e-12 -10.0e-12; -10.0e-12 100.0e-12];
        length_m=10.0,
        conductor_names=[:phase, :enclosure],
        enclosure_reference=:ground,
    )
    @test size(gis_gil_series_impedance_ohm(gis, 1.0e6)) == (2, 2)
    @test size(gis_gil_shunt_admittance_s(gis, 1.0e6)) == (2, 2)

    plan = InsulationStudyPlan(
        :synthetic_coordination;
        sample_count=1024,
        seed=20260824,
        stress_mean_v=800.0e3,
        stress_standard_deviation_v=80.0e3,
        strength_mean_v=1.0e6,
        strength_standard_deviation_v=100.0e3,
        stress_strength_correlation=0.2,
    )
    first_summary = run_insulation_study(plan)
    second_summary = run_insulation_study(plan)
    @test first_summary.signature == second_summary.signature
    @test first_summary.failure_count == second_summary.failure_count
    @test 0.0 <= first_summary.confidence_lower <=
        first_summary.empirical_failure_probability <=
        first_summary.confidence_upper <= 1.0
    @test deterministic_insulation_margin(800.0e3, 1.0e6) == 200.0e3

    diagnostic = SurgeInsulationDiagnostic(:energy_residual, 1.0e-10, 1.0, 1.0e-9)
    result = SurgeInsulationResult(
        :public_contract,
        :accepted,
        (sample_count=first_summary.sample_count,),
        [diagnostic],
    )
    @test surge_result_accepted(result)
    @test length(result.signature) == 64

    product_provenance = surge_physical_parameter_provenance(
        "AIMORA-authored generic public surge product",
        "SI peak surge quantities",
        "direct deterministic synthetic inputs",
        "field and device uncertainty unknown",
        "public scientific demonstration only",
    )
    product_arc() = CombinedArcBranch(
        1,
        0;
        cassie_power_w=2.0e5,
        mayr_power_w=2.0e3,
        cassie_time_constant_s=20.0e-6,
        mayr_time_constant_s=5.0e-6,
        transition_power_w=2.0e4,
        initial_conductance_s=10.0,
        initially_ignited=false,
    )
    product_gap(id) = VacuumInterruptionState(
        id;
        chopping_current_a=5.0,
        initial_dielectric_strength_v=1.0e3,
        dielectric_recovery_rate_v_per_s=1.0e9,
        maximum_dielectric_strength_v=100.0e3,
    )
    breaker_timing_provenance = AIMORA.StudyCore.ParameterProvenance(
        "AIMORA-authored generic public surge product",
        "integer ticks",
        "direct deterministic fixed-step calendar",
        "exact synthetic timing values",
        "three-pole public interruption product",
        AIMORA.StudyCore.NumericalPolicyParameter,
    )
    breaker_specification = AIMORA.ProtectionStudy.EMTBreakerSpecification(
        :generic_breaker;
        closed_conductance_s=1.0e3,
        open_conductance_s=1.0e-9,
        opening_travel_ticks=0,
        closing_travel_ticks=0,
        current_zero_required=false,
        current_zero_threshold_a=5.0,
        failure_delay_ticks=100,
        failure_current_threshold_a=5.0,
        reclose_dead_ticks=10,
        reclaim_ticks=10,
        maximum_reclose_shots=1,
        contact_tail_enabled=true,
        physical_provenance=product_provenance,
        timing_provenance=breaker_timing_provenance,
    )
    line_runtime_preparation = surge_test_line_runtime_preparation(0.1e-6)
    apparatus_preparation = surge_test_transformer_preparation(
        0.1e-6,
        product_provenance,
    )
    interruption_arcs = (product_arc(), product_arc(), product_arc())
    interruption_gaps = (product_gap(:pole_a), product_gap(:pole_b), product_gap(:pole_c))
    interruption_branches = ntuple(
        index -> VacuumInterruptionBranch(
            1,
            0,
            interruption_gaps[index],
            interruption_arcs[index],
        ),
        3,
    )
    product_specs = [
        SurgeInsulationProductSpecification(
            :generic_three_pole_interruption,
            InterruptionRestrikeProduct,
            (
                breaker_specification,
                arc_branches=interruption_arcs,
                vacuum_gaps=interruption_gaps,
                vacuum_branches=interruption_branches,
            );
            timestep_s=1.0e-6,
            stop_time_s=100.0e-6,
            provenance=product_provenance,
            uncertainty="generic arc and vacuum parameters with unknown apparatus uncertainty",
            validity_domain="three synthetic poles on fixed-step instantaneous EMT",
        ),
        SurgeInsulationProductSpecification(
            :generic_arrester_terminal,
            ArresterProtectedTerminalProduct,
            (
                apparatus_identity=:generic_transformer_terminal,
                apparatus_preparation,
                arrester=deepcopy(arrester),
                dynamic_equivalent=dynamic_arrester,
                terminal_section=section,
            );
            timestep_s=0.1e-6,
            stop_time_s=100.0e-6,
            provenance=product_provenance,
            uncertainty="generic arrester and terminal parameters with unknown apparatus uncertainty",
            validity_domain="one synthetic protected apparatus terminal",
        ),
        SurgeInsulationProductSpecification(
            :generic_tower_backflash,
            TowerBackflashProduct,
            (
                line_terminal_identity=:generic_line_terminal,
                line_runtime_preparation,
                lightning_sequence=sequence,
                tower=tower,
                ground=deepcopy(ionizing_ground),
                insulator=deepcopy(disruptive),
            );
            timestep_s=0.1e-6,
            stop_time_s=500.0e-6,
            provenance=product_provenance,
            uncertainty="generic lightning tower ground and insulation uncertainty",
            validity_domain="one synthetic direct-strike or backflash surge path",
        ),
        SurgeInsulationProductSpecification(
            :generic_gis_corona_terminal,
            GISCoronaTerminalProduct,
            (
                terminal_identity=:generic_gis_terminal,
                gis_section=gis,
                corona=deepcopy(corona),
                insulator=deepcopy(leader),
            );
            timestep_s=0.01e-6,
            stop_time_s=20.0e-6,
            provenance=product_provenance,
            uncertainty="generic GIS corona and spacer uncertainty",
            validity_domain="one synthetic GIS/GIL terminal section",
        ),
        SurgeInsulationProductSpecification(
            :generic_statistical_insulation,
            StatisticalInsulationProduct,
            (study_plan=plan,);
            timestep_s=1.0,
            stop_time_s=1.0,
            provenance=product_provenance,
            uncertainty="explicit synthetic Gaussian stress-strength distributions",
            validity_domain="preregistered deterministic seeded synthetic ensemble",
        ),
    ]
    @test length(unique(spec.deterministic_signature_sha256 for spec in product_specs)) == 5
    for specification in product_specs
        preparation = prepare_surge_insulation_product(specification)
        public_readiness = surge_insulation_product_readiness(preparation)
        private_readiness = surge_insulation_product_readiness(
            preparation;
            production_backend_available=true,
        )
        @test public_readiness.ready
        @test public_readiness.code === :ready_for_solver_free_inspection
        @test private_readiness.code === :ready_for_coupled_execution
        @test length(preparation.preparation_signature_sha256) == 64
    end
    @test_throws ArgumentError SurgeInsulationProductSpecification(
        :missing_tower_components,
        TowerBackflashProduct,
        (line_terminal_identity=:line,);
        timestep_s=1.0e-6,
        stop_time_s=2.0e-6,
        uncertainty="explicit",
        validity_domain="synthetic",
    )

    wrong_arrester = MetalOxideArrester(2, 0, characteristic)
    @test_throws ArgumentError restore_surge_device!(wrong_arrester, arrester_snapshot)
end

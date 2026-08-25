using Test
using AIMORA
using LinearAlgebra
using Printf

include("nonlinear_network.jl")
include("native_extensions.jl")
include("coupled_line_fitting.jl")
include("coupled_line_runtime.jl")
include("emt_task_platform.jl")
include("emt_partitioning.jl")
include("dassl_class_contracts.jl")
include("performance_realtime_contracts.jl")
include("extended_vsc_control_filter_platform.jl")
include("transformer_apparatus.jl")
include("modern_machine_families.jl")
include("measurement_chains.jl")
include("protection_platform.jl")
include("surge_insulation.jl")
include("portable_emt_snapshots.jl")
include("converter_systems.jl")

@testset "public package isolation" begin
    for module_name in (
        :LegacyEMTP,
        :ExecutableComparison,
        :TranslationProgress,
        :ValidationBenchmarks,
    )
        @test !isdefined(AIMORA, module_name)
    end
end

@testset "public solver backend contract" begin
    @test AIMORA.AbstractAIMORASolverBackend isa DataType
    if AIMORA.solver_available()
        @test :emt in getfield.(AIMORA.solver_capabilities(), :id)
        @test AIMORA.backend_metadata().name == :aimora_production
    else
        @test isempty(AIMORA.solver_capabilities())
        @test AIMORA.backend_metadata() isa AIMORA.SolverUnavailableResult
        unavailable = AIMORA.prepare_study(nothing, :emt)
        @test unavailable isa AIMORA.SolverUnavailableResult
        @test unavailable.operation == :prepare_study
        @test unavailable.required_capability == :study_preparation
        @test unavailable.mode == :open_core
        measurement_unavailable = AIMORA.materialize_measurement_branches(Any[])
        @test measurement_unavailable isa AIMORA.SolverUnavailableResult
        @test measurement_unavailable.operation == :materialize_measurement_branches
        @test measurement_unavailable.required_capability ==
            :measurement_network_materialization
        breaker_unavailable = AIMORA.materialize_emt_breaker_poles(
            nothing,
            nothing,
            ((1, 2), (3, 4), (5, 6)),
        )
        @test breaker_unavailable isa AIMORA.SolverUnavailableResult
        @test breaker_unavailable.operation == :materialize_emt_breaker_poles
        @test breaker_unavailable.required_capability == :emt_protection_breaker
        converter_unavailable = AIMORA.prepare_converter_system(:average_buck_converter)
        @test converter_unavailable isa AIMORA.SolverUnavailableResult
        @test converter_unavailable.operation == :prepare_converter_system
        @test converter_unavailable.required_capability == :extended_converter_systems
    end
end

if AIMORA.solver_available()
    @testset "nonlinear nodal step callback admission" begin
        callback = (context, injections, constraints) ->
            (context, injections, constraints)
        @test AIMORA.EMTStudy._execute_nonlinear_nodal_step_solver(
            callback,
            :context,
            [1.0],
            :constraints,
        ) == (:context, [1.0], :constraints)
        @test_throws ArgumentError AIMORA.EMTStudy._execute_nonlinear_nodal_step_solver(
            :not_callable,
            :context,
            [1.0],
            :constraints,
        )
    end
end

@testset "open engineering core" begin
    issue = AIMORA.ValidationCore.missing_data("line", "length is required")
    result = AIMORA.ValidationCore.validation_result(source = "public package test")
    AIMORA.ValidationCore.add_issue!(result, issue)
    @test !AIMORA.ValidationCore.is_valid(result)

    scenario = AIMORA.ProjectData.Scenario(id = :base, name = "Base")
    row = AIMORA.InverterAssets.inverter_row(
        id = :inv1,
        bus = :bus1,
        rated_kva = 1000.0,
        v_ll_rms_v = 4160.0,
    )
    AIMORA.InverterAssets.set_inverter_table!(scenario, [row])
    @test only(AIMORA.InverterAssets.inverter_table(scenario))[:id] == :inv1

    profile = AIMORA.StudyInputProfiles.input_profile(:power_flow)
    @test :buses in AIMORA.StudyInputs.required_keys(profile)
    @test :emt in (study.id for study in AIMORA.StudyCatalog.available_studies())
end

if AIMORA.solver_available()
    @testset "complementary power-semiconductor bridge public contract" begin
        upper = AIMORA.Nonlinear.MOSFETSwitch(
            1,
            2;
            gate_driver = AIMORA.Nonlinear.PowerSemiconductorGateDriver(),
            antiparallel_diode = AIMORA.Nonlinear.AntiparallelDiodeParameters(),
        )
        lower = AIMORA.Nonlinear.MOSFETSwitch(
            2,
            0;
            gate_driver = AIMORA.Nonlinear.PowerSemiconductorGateDriver(),
            antiparallel_diode = AIMORA.Nonlinear.AntiparallelDiodeParameters(),
        )
        bridge = AIMORA.Nonlinear.PowerSemiconductorBridgeLeg(
            upper,
            lower;
            commutation_dead_time_s = 1.0e-6,
        )
        @test AIMORA.Nonlinear.request_power_semiconductor_bridge_pole!(
            bridge,
            true,
            0.0,
        ) == AIMORA.Nonlinear.BRIDGE_GATE_ACCEPTED
        @test AIMORA.Nonlinear.power_semiconductor_bridge_gate_transition_time(bridge) ==
            1.0e-6
        @test AIMORA.Nonlinear.apply_power_semiconductor_bridge_gate_transitions!(
            bridge,
            1.0e-6,
        ) == 1
        @test upper.gate_driver.applied_on
        @test !lower.gate_driver.applied_on
        @test AIMORA.Nonlinear.request_power_semiconductor_bridge_gates!(
            bridge,
            true,
            true,
            2.0e-6,
        ) == AIMORA.Nonlinear.BRIDGE_GATE_SHOOT_THROUGH_REJECTED
        @test AIMORA.Branches.trace_output_channel_count(bridge) == 24
        @test AIMORA.Branches.trace_output_is_public(bridge)
        @test AIMORA.EMTStudy.PowerSemiconductorBridgePoleCommand(:main_bridge).element_name ==
            :main_bridge
    end

    @testset "floating Thevenin source and coupled VSC runtime" begin
        floating_source = AIMORA.Branches.TwoTerminalTheveninSource(
            1,
            2,
            2.0,
            _time_s -> 10.0,
        )
        source_system = AIMORA.Nodal.NodalSystem(
            2,
            Any[
                floating_source,
                AIMORA.Branches.ConductanceBranch(1, 2, 1.0),
                AIMORA.Branches.ConductanceBranch(2, 0, 1.0e3),
            ],
        )
        AIMORA.Nodal.solve_algebraic_state!(source_system, 0.0, 1.0e-6)
        @test source_system.v[1] - source_system.v[2] ≈ 20.0 / 3.0 atol = 2.0e-12
        @test source_system.v[2] ≈ 0.0 atol = 2.0e-12
        @test AIMORA.Nodal.accept_algebraic_state!(source_system, 1.0e-6) ===
            source_system.v

        parameters = AIMORA.SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters(
            end_time_s = 20.0e-3,
            sag_start_s = 1.0e-3,
            sag_end_s = 2.0e-3,
            fault_start_s = 8.0e-3,
            block_time_s = 9.0e-3,
            fault_end_s = 10.0e-3,
            restart_time_s = 11.0e-3,
            harmonic_window_start_s = 0.0,
            harmonic_window_end_s = 20.0e-3,
        )
        trace = AIMORA.EMTStudy.simulate_three_phase_vsc(parameters)
        metrics = trace.metrics
        @test length(trace.time_s) == 20_001
        @test metrics.finite_output
        @test metrics.exact_boundary_alignment
        @test metrics.boundary_count == 6
        @test metrics.block_count == 1
        @test metrics.restart_count == 1
        @test metrics.commutation_count > 0
        @test metrics.maximum_topology_iterations <= 12
        @test metrics.maximum_nodal_kcl_residual_a <= 1.0e-7
        @test metrics.relative_energy_residual <= 1.0e-5
        @test metrics.relative_dc_ac_energy_residual <= 1.0e-5
        @test minimum(trace.dc_link_voltage_v) >= 650.0
        @test maximum(abs, trace.filter_current_a) <= parameters.current_limit_a
        @test maximum(abs, vec(sum(trace.filter_current_a; dims = 1))) <= 1.0e-7
        @test any(trace.upper_antiparallel_diode_conducting .> 0.5) ||
            any(trace.lower_antiparallel_diode_conducting .> 0.5)
    end

    @testset "grounded scalar branch references" begin
        fixed_control = [
            "BEGIN NEW DATA CASE",
            "POWER FREQUENCY                     60.0",
            "ABSOLUTE U.M. DIMENSIONS              20       2      50      60",
            "C PRINTED NUMBER WIDTH  13  2",
            " .000200    .100",
            "       1       1       1       1       1      -1                               1",
            "       5       5      20      20     100     100",
            "C BRANCHES",
        ]
        fixed_terminators = [
            "BLANK card ending branch cards",
            "BLANK card ending nonexistent switch cards",
            "BLANK card ending all electric-network sources",
            "BLANK card terminating output variable requests",
            "BLANK card terminating plot cards",
            "BEGIN NEW DATA CASE",
        ]
        parsed = AIMORA.DeckParser.parse_deck_lines(
            vcat(
                fixed_control,
                [
                    "  BUSAS2                  1.0E+6",
                    "  BUSBS2      BUSAS2",
                    "  BUSCS2      BUSAS2",
                ],
                fixed_terminators,
            );
            source = "grounded-scalar-reference-contract",
        )
        @test AIMORA.ValidationCore.is_valid(parsed.validation)
        @test length(parsed.elements) == 3
        @test all(element -> element isa AIMORA.Branches.ConductanceBranch, parsed.elements)
        @test getfield.(parsed.elements, :g) == fill(1.0e-6, 3)
        @test get(
            parsed.card_counts,
            :fixed_grounded_scalar_branch_reference,
            0,
        ) == 2
        @test get(
            parsed.card_counts,
            :deferred_single_terminal_capacitance_group,
            0,
        ) == 0
        @test getfield.(parsed.over2_branch_rows, :layout_kind) == [
            :fixed_sparse_numeric,
            :fixed_grounded_scalar_branch_reference,
            :fixed_grounded_scalar_branch_reference,
        ]
        @test getfield.(parsed.over2_branch_rows, :reference_name) == [
            :none,
            :branch_fixed_1,
            :branch_fixed_1,
        ]

        missing_reference = AIMORA.DeckParser.parse_deck_lines(
            vcat(
                fixed_control,
                ["  BUSBS2      ABSENT"],
                fixed_terminators,
            );
            source = "missing-grounded-scalar-reference-contract",
        )
        @test !AIMORA.ValidationCore.is_valid(missing_reference.validation)
        @test any(
            issue -> occursin(
                "does not match a prior accepted scalar branch owner",
                issue.message,
            ),
            missing_reference.validation.issues,
        )
    end
end

@testset "public inverter model" begin
    contract = AIMORA.Inverter.inverter_contract()
    @test contract.id == :average_value_grid_following_inverter
    @test contract.maturity == :prototype
    @test contract.fidelity == AIMORA.StudyCore.AverageValue
    @test :switching_ripple in contract.validity_domain.unsupported_phenomena
    @test contract.state_inventory.differential == (:id_a, :iq_a, :xid, :xiq)
    @test contract.mutation_order == (
        :sample_outputs,
        :read_power_references,
        :evaluate_rk4_stages,
        :apply_voltage_limit,
        :commit_state,
    )
    rows = AIMORA.Inverter.simulate_inverter(t_end = 2.0e-4, dt = 20.0e-6)
    summary = AIMORA.Inverter.inverter_summary(rows)
    @test summary.samples == 11
    @test all(isfinite, rows[end])
    @test rows == AIMORA.Inverter.simulate_inverter(t_end = 2.0e-4, dt = 20.0e-6)
    @test_throws AIMORA.StudyCore.ValidityDomainError AIMORA.Inverter.simulate_inverter(
        t_end = 2.0e-4,
        dt = 20.0e-6,
        fidelity = AIMORA.StudyCore.SwitchingDetailed,
    )
    @test_throws AIMORA.StudyCore.ValidityDomainError AIMORA.Inverter.simulate_inverter(
        t_end = 2.0e-4,
        dt = 20.0e-6,
        p = AIMORA.Inverter.InverterParams(f_hz = 0.0),
    )
end

@testset "switch-detailed three-phase VSC public model" begin
    contract = AIMORA.SwitchDetailedVSC.switch_detailed_vsc_contract()
    @test contract.id == :three_phase_two_level_switch_detailed_vsc
    @test contract.maturity == :implemented
    @test contract.fidelity == AIMORA.StudyCore.SwitchingDetailed
    @test Set(quantity.key for quantity in contract.outputs) >= Set((
        :dc_source_energy_j,
        :ac_terminal_energy_j,
        :dc_ac_energy_residual_j,
    ))
    @test :phase_locked_loop_dynamics in contract.validity_domain.unsupported_phenomena
    @test :transformer_magnetizing_saturation in
        contract.validity_domain.unsupported_phenomena

    phase = (310.0, -122.0, -188.0)
    angle = 0.37
    stationary = AIMORA.SwitchDetailedVSC.clarke_transform(phase)
    synchronous = AIMORA.SwitchDetailedVSC.park_transform(stationary, angle)
    recovered = AIMORA.SwitchDetailedVSC.inverse_clarke_transform(
        AIMORA.SwitchDetailedVSC.inverse_park_transform(synchronous, angle),
    )
    @test collect(recovered) ≈ collect(phase) atol = 2.0e-13
    @test AIMORA.SwitchDetailedVSC.instantaneous_three_phase_power(
        (100.0, -50.0, -50.0),
        (10.0, -5.0, -5.0),
        0.0,
    ) == (active_w = 1500.0, reactive_var = 0.0)
    sinusoidal = AIMORA.SwitchDetailedVSC.modulation_duties(
        (240.0, -120.0, -120.0),
        800.0,
        AIMORA.SwitchDetailedVSC.SinusoidalPulseWidthModulation,
    )
    injected = AIMORA.SwitchDetailedVSC.modulation_duties(
        (240.0, -120.0, -120.0),
        800.0,
        AIMORA.SwitchDetailedVSC.ZeroSequenceInjectedPulseWidthModulation,
    )
    @test collect(sinusoidal) ≈ [0.8, 0.35, 0.35] atol = 2.0e-16
    @test collect(injected) ≈ [0.725, 0.275, 0.275] atol = 2.0e-16
    @test_throws ArgumentError AIMORA.SwitchDetailedVSC.validate_three_phase_vsc_parameters(
        AIMORA.SwitchDetailedVSC.ThreePhaseTwoLevelVSCParameters(
            scheduler_tick_s = 3.0e-6,
        ),
    )
end

@testset "transformer parameter studies" begin
    short_circuit_lines = [
        "XFORMER",
        "22.        700.",
        "139.4      13.6     2100.     12.       700.",
        "\$PUNCH",
        "BRANCH  NAME1 NAME2 NAME3 NAME4 NAME5 NAME6",
        "3.3        83.3",
        "132.8      250.     6.7       83.3",
        "66.4       56.8     5.1       18.96",
        "13.2       56.8     3.2       18.96",
        "BLANK card ending XFORMER cases",
    ]
    short_deck =
        AIMORA.TransformerParameterInput.parse_over41_transformer_parameter_lines(
            short_circuit_lines;
            source = "short-circuit-contract",
        )
    @test length(short_deck.cases) == 2
    short_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study(
            short_deck,
        )
    @test short_study.generated_branch_count == 5
    @test short_study.physical_checks_passed
    two_winding = short_study.case_results[1]
    three_winding = short_study.case_results[2]
    @test two_winding.admittance_matrix_s *
          two_winding.impedance_matrix_ohm ≈
          -I atol = 2.0e-12
    @test three_winding.admittance_matrix_s *
          three_winding.impedance_matrix_ohm ≈
          -I atol = 1.0e-10
    @test two_winding.impedance_matrix_ohm ≈
          transpose(two_winding.impedance_matrix_ohm) atol = 1.0e-12
    @test eigmin(Symmetric(real.(two_winding.impedance_matrix_ohm))) >=
          -1.0e-12
    @test two_winding.generated_branches[2].resistance_values_ohm[1] ≈
          -0.002028822586562 atol = 5.0e-13
    @test two_winding.generated_branches[2].
          inductance_or_reactance_values[1] ≈
          135.3359666262 atol = 5.0e-10
    @test getfield.(two_winding.generated_branches, :from_node) ==
          [Symbol(""), Symbol("")]
    @test getfield.(two_winding.generated_branches, :to_node) ==
          [Symbol(""), Symbol("")]
    @test getfield.(
        three_winding.generated_branches,
        :from_node,
    ) == [:NAME1, :NAME3, :NAME5]

    saturable_header =
        rpad("XFORMER", 32) * @sprintf("%8.1f%8.1f", 60.0, 0.0)
    saturable_case =
        @sprintf(
            "  %-6s%18s%6.2f%6.2f%6s%6.2f",
            "TRANSF",
            "",
            0.5,
            1.2,
            "",
            1000.0,
        )
    saturable_winding_1 =
        @sprintf(
            "%2d%-6s%-6s%12s%6.2f%6.2f%6.2f",
            1,
            "H1",
            "H0",
            "",
            0.2,
            5.0,
            132.0,
        )
    saturable_winding_2 =
        @sprintf(
            "%2d%-6s%-6s%12s%6.2f%6.2f%6.2f",
            2,
            "L1",
            "L0",
            "",
            0.1,
            2.0,
            33.0,
        )
    saturable_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
            [
                saturable_header,
                saturable_case,
                saturable_winding_1,
                saturable_winding_2,
                "BLANK",
            ];
            source = "saturable-contract",
        )
    saturable = only(saturable_study.case_results)
    @test saturable.magnetizing_parallel_inductance_h ≈ 2.4
    @test saturable.branch_resistance_matrix_ohm ≈
          transpose(saturable.branch_resistance_matrix_ohm)
    @test saturable.branch_inductance_or_reactance_matrix ≈
          1000.0 .* saturable.physical_inductance_matrix_h
    @test eigmin(Symmetric(saturable.physical_inductance_matrix_h)) >=
          -1.0e-12
    @test saturable.physical_checks_passed

    bctran_lines = [
        "ACCESS MODULE BCTRAN",
        " 360.       .428      300.      135.73    .428      300.      135.73       1 3 1",
        "  1132.79056 .2054666   H-1         H-2         H-3",
        "  263.393059 .0742333   L-1         L-2         L-3",
        "  350.       .0822      T-1   T-2   T-2               T-1",
        " 1 20.        8.74      300.      7.3431941 300.       3 1",
        " 1 30.        8.68      76.       26.258183 300.",
        " 2 30.        5.31      76.       18.552824 300.",
        "BLANK card that terminates short-circuit test data of BCTRAN",
        "BEGIN NEW DATA CASE",
    ]
    bctran_study =
        AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
            bctran_lines;
            source = "bctran-contract",
        )
    bctran = only(bctran_study.case_results)
    @test bctran.input.output_representation == :reactance
    @test size(bctran.reactance_matrix_ohm) == (9, 9)
    @test bctran.winding_resistances_ohm ≈
          [0.2054666, 0.0742333, 0.0822]
    @test bctran.inverse_inductance_matrix_per_h[1, 1] ≈
          26.512692374898 atol = 5.0e-11
    @test bctran.inverse_inductance_matrix_per_h[2, 1] ≈
          -59.57848438329 atol = 5.0e-11
    @test bctran.reactance_matrix_ohm[1, 1] ≈
          41432.097487193 atol = 2.0e-8
    @test bctran.reactance_matrix_ohm[4, 1] ≈
          -0.0533106395901 atol = 5.0e-8
    @test only(bctran.magnetizing_shunts).self_resistance_ohm ≈
          55098.277352343 atol = 5.0e-8
    @test bctran.positive_pair_reconstruction_residual <= 1.0e-14
    @test bctran.zero_pair_reconstruction_residual <= 1.0e-14
    @test bctran.inverse_inductance_matrix_per_h ≈
          transpose(bctran.inverse_inductance_matrix_per_h) atol = 1.0e-12
    @test eigmin(Symmetric(bctran.inverse_inductance_matrix_per_h)) >
          0.0
    @test bctran.physical_checks_passed

    legacy_bctran_lines = copy(bctran_lines)
    legacy_bctran_lines[1] = "XFORMER, 44."
    legacy_bctran =
        only(
            AIMORA.TransformerParameterStudy.run_transformer_parameter_study_lines(
                legacy_bctran_lines;
                source = "legacy-bctran-contract",
            ).case_results,
        )
    @test legacy_bctran.reactance_matrix_ohm ≈
          bctran.reactance_matrix_ohm atol = 0.0

    report =
        AIMORA.TransformerParameterReport.transformer_parameter_report_text(
            bctran_study,
        )
    @test report ==
          AIMORA.TransformerParameterReport.transformer_parameter_report_text(
        bctran_study,
    )
    @test occursin("CASE 1 MULTIPHASE_TRANSFORMER", report)
    @test occursin("GENERATED_BRANCH_COUNT 9", report)
    mktempdir() do directory
        path = joinpath(directory, "transformer_parameters.txt")
        @test AIMORA.TransformerParameterReport.write_transformer_parameter_report(
            path,
            bctran_study,
        ) == abspath(path)
        @test read(path, String) == report
    end

    @test_throws ArgumentError AIMORA.TransformerParameterInput.
                               parse_over41_transformer_parameter_lines(
        bctran_lines[1:6],
    )
    @test_throws ArgumentError AIMORA.TransformerParameters.
                               TransformerShortCircuitCase(
        1,
        [132.0, 33.0],
        [10000.0],
        [1.0],
        [100.0],
        1.0,
        100.0,
        Tuple[],
    ) |>
                               AIMORA.TransformerParameters.
                               transformer_short_circuit_parameters
end

if AIMORA.solver_available()
    @testset "power-semiconductor public model and gate surface" begin
        driver = AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            turn_on_delay_s = 2.0e-6,
            minimum_pulse_width_s = 3.0e-6,
        )
        device = AIMORA.Nonlinear.IGBTSwitch(
            1,
            0;
            gate_driver = driver,
            forward_voltage_drop_v = 0.7,
            on_conductance = 5.0,
            antiparallel_diode = AIMORA.Nonlinear.AntiparallelDiodeParameters(
                forward_voltage_v = 0.6,
                on_conductance_s = 5.0,
            ),
            snubber = AIMORA.Nonlinear.SeriesRCSnubber(10.0, 1.0e-6),
        )
        device.last_voltage = 2.0
        @test !AIMORA.Nonlinear.request_power_semiconductor_gate!(device, true, 0.0)
        @test AIMORA.Nonlinear.power_semiconductor_gate_transition_time(device) == 2.0e-6
        @test AIMORA.Nonlinear.apply_power_semiconductor_gate_transition!(
            device,
            2.0e-6,
        )
        @test device.closed
        terminal = AIMORA.Nonlinear.power_semiconductor_terminal_state(device)
        @test terminal.device_kind == :igbt
        @test terminal.gate_applied_on
        @test AIMORA.Branches.trace_output_channel_count(device) == 10
        @test AIMORA.Branches.trace_output_is_public(device)
        @test AIMORA.EMTStudy.PowerSemiconductorGateCommand(:main_valve).element_name ==
            :main_valve
    end

    @testset "sampled line steady-state terminal admittance" begin
        propagation = AIMORA.Lines.line_weighting_samples(
            (2:7) .* 1.0e-6,
            [0.01, 0.15, 0.50, 0.20, 0.04, 0.01],
        )
        admittance = AIMORA.Lines.line_weighting_samples(
            (0:5) .* 1.0e-6,
            [0.01, 0.10, 0.40, 0.20, 0.05, 0.01],
        )
        coefficients = AIMORA.Lines.sampled_line_weighting_coefficients(
            propagation,
            admittance,
            1.0e-6,
            300.0;
            propagation_peak_index = 3,
            admittance_rise_index = 2,
        )
        line = AIMORA.Lines.sampled_frequency_dependent_line(
            1,
            2,
            coefficients,
        )
        terminal_admittance =
            AIMORA.Lines.sampled_line_steady_state_terminal_admittance(
                line,
                60.0,
            )
        @test size(terminal_admittance) == (2, 2)
        @test all(isfinite, terminal_admittance)
        @test terminal_admittance ≈ transpose(terminal_admittance)
        @test_throws ArgumentError AIMORA.Lines.
                                   sampled_line_steady_state_terminal_admittance(
            line,
            -60.0,
        )
    end

    @testset "retired input-converter dispositions" begin
        parsed = AIMORA.DeckParser.parse_deck_lines(
            [
                "BEGIN NEW DATA CASE",
                "CS",
                "CZ",
                "BLANK CARD TERMINATING THE CASE",
            ];
            source = "retired-input-converter-contract",
        )
        @test AIMORA.ValidationCore.is_valid(parsed.validation)
        auxiliary = AIMORA.EMTStudy.run_deck_auxiliary_studies(parsed)
        @test auxiliary.compatibility_exclusions == [
            :pre_m37_zinc_oxide_card_converter,
            :pre_m37_switch_pseudononlinear_card_converter,
        ]
        @test isempty(auxiliary.deferred_requests)
    end

    @testset "private solver integration" begin
        system = AIMORA.Nodal.NodalSystem(2, [
            AIMORA.Branches.TheveninSource(1, 1.0e9, _ -> 1.0),
            AIMORA.Branches.ConductanceBranch(1, 2, 1.0),
            AIMORA.Branches.ConductanceBranch(2, 0, 1.0),
        ])
        voltage = AIMORA.Nodal.solve_step!(system, 0.0, 20.0e-6)
        @test voltage[2] ≈ 0.5 atol = 1.0e-6
        @test AIMORA.solver_status().mode == :full_engine
    end

    @testset "reversible timestep transaction foundation" begin
        shared_history = [1.0, 2.0]
        owner = (
            primary_history = shared_history,
            aliased_history = shared_history,
            cursor = Ref(3),
        )
        transaction = AIMORA.OVER16TimestepIntegration.TimestepTransaction(owner)
        AIMORA.OVER16TimestepIntegration.begin_timestep_transaction!(transaction)
        push!(owner.primary_history, 4.0)
        owner.cursor[] = 9
        @test_throws ArgumentError AIMORA.OVER16TimestepIntegration.
            begin_timestep_transaction!(transaction)
        AIMORA.OVER16TimestepIntegration.restore_timestep_transaction!(transaction)
        @test owner.primary_history == [1.0, 2.0]
        @test owner.primary_history === owner.aliased_history
        @test owner.cursor[] == 3

        AIMORA.OVER16TimestepIntegration.begin_timestep_transaction!(transaction)
        owner.primary_history[1] = 7.0
        AIMORA.OVER16TimestepIntegration.commit_timestep_transaction!(transaction)
        @test owner.aliased_history[1] == 7.0
        AIMORA.OVER16TimestepIntegration.begin_timestep_transaction!(transaction)
        owner.primary_history[1] = 11.0
        AIMORA.OVER16TimestepIntegration.restore_timestep_transaction!(transaction)
        @test owner.primary_history[1] == 7.0
        @test AIMORA.OVER16TimestepIntegration.timestep_transaction_status(
            transaction,
        ) == (
            active = false,
            capture_count = 3,
            restore_count = 2,
            commit_count = 1,
        )
    end

    @testset "hybrid event and sampled-task public surface" begin
        policy = AIMORA.EMTStudy.EMTHybridEventPolicy(
            root_time_tolerance_s = 1.0e-10,
            root_value_tolerance = 1.0e-12,
            simultaneity_tolerance_s = 1.0e-10,
        )
        surface = AIMORA.EMTStudy.EMTHybridEventSurface(
            :public_rising_root,
            owner -> owner.value,
            (owner, time_s) -> begin
                owner.transition_time_s[] = time_s
                owner.transition_count[] += 1
            end;
            direction = :rising,
            priority = -1,
            repeatable = false,
        )
        owner = (value = -1.0, transition_time_s = Ref(Inf), transition_count = Ref(0))
        mutable_owner = (
            value = owner.value,
            transition_time_s = owner.transition_time_s,
            transition_count = owner.transition_count,
        )
        AIMORA.OVER16TimestepIntegration.apply_hybrid_event!(
            surface,
            (
                value = mutable_owner.value,
                transition_time_s = mutable_owner.transition_time_s,
                transition_count = mutable_owner.transition_count,
            ),
            0.25,
        )
        @test policy.maximum_events_per_step == 64
        @test mutable_owner.transition_time_s[] == 0.25
        @test mutable_owner.transition_count[] == 1

        task_times = Float64[]
        task = AIMORA.EMTStudy.EMTExactSampledTask(
            :public_task,
            3.0e-6,
            (_owner, time_s, _execution_index) -> push!(task_times, time_s);
            tick_s = 1.0e-6,
            first_time_s = 0.0,
        )
        scheduler = AIMORA.EMTStudy.EMTExactSampledTaskScheduler(
            1.0e-6;
            tasks = [task],
        )
        for tick in 0:9
            AIMORA.EMTStudy.run_due_emt_sampled_tasks!(
                scheduler,
                nothing,
                tick * 1.0e-6,
            )
        end
        @test task_times == [0.0, 3.0e-6, 6.0e-6, 9.0e-6]

        sampled_owner = (
            input = Ref(0.4),
            output = Ref(0.0),
            gate = Ref(false),
            edges = Tuple{Bool,Float64,Int}[],
        )
        sampled_control = AIMORA.EMTStudy.EMTExactSampledControlTask(
            :public_sampled_control,
            10.0e-6,
            (owner, _time_s, _sample_index) -> owner.input[],
            (_owner, input, _time_s, _sample_index) -> input,
            (owner, output, _time_s, _sample_index) -> (owner.output[] = output);
            tick_s = 1.0e-6,
            computational_delay_s = 2.0e-6,
            initial_output = 0.0,
            priority = -1,
        )
        pwm = AIMORA.EMTStudy.EMTExactPWMTask(
            :public_pwm,
            10.0e-6,
            (owner, _time_s, _cycle_index) -> owner.output[],
            (owner, high, time_s, edge_index) -> begin
                owner.gate[] = high
                push!(owner.edges, (high, time_s, edge_index))
            end;
            tick_s = 1.0e-6,
            priority = 1,
            power_history_invalidating = false,
        )
        multirate_scheduler = AIMORA.EMTStudy.EMTExactSampledTaskScheduler(
            1.0e-6;
            tasks = [pwm, sampled_control],
        )
        for tick in 0:25
            AIMORA.EMTStudy.run_due_emt_sampled_tasks!(
                multirate_scheduler,
                sampled_owner,
                tick * 1.0e-6,
            )
        end
        @test sampled_control.sample_count == 3
        @test sampled_control.write_count == 3
        @test sampled_control.held_output == 0.4
        @test pwm.cycle_count == 3
        @test getindex.(sampled_owner.edges, 1) == [true, false, true, false]
        @test getfield.(pwm.occurrences, :tick) == [10, 14, 20, 24]
        @test !sampled_owner.gate[]

        arc = AIMORA.TACS.ControlledSwitchDelayedArcState(1.0, 1.0, 1.0, 0.0)
        arc.opening_requested = true
        controlled = AIMORA.TACS.TACSControlledSwitch(
            1,
            0,
            Ref(1.0);
            initially_closed = true,
            delayed_arc = arc,
        )
        AIMORA.TACS.apply_controlled_switch_delayed_arc_transition!(
            controlled,
            2.0e-6,
        )
        @test !controlled.closed
        @test controlled.delayed_arc.tail_active
        @test controlled.delayed_arc.scheduled_open_time_s == 2.0e-6
        @test controlled.delayed_arc.transition_count == 1
        AIMORA.TACS.apply_controlled_switch_delayed_arc_transition!(
            controlled,
            2.0e-6,
        )
        @test controlled.delayed_arc.transition_count == 1
    end
else
    @testset "public checkout has no solver source" begin
        @test AIMORA.solver_status().mode == :open_core
        @test AIMORA.solver_status().backend === nothing
        @test !isdefined(AIMORA, :Nodal)
        @test !isdefined(AIMORA, :Nonlinear)
        @test_throws ErrorException AIMORA.require_solver()
    end
end

if AIMORA.solver_available()
    include("extended_power_semiconductor_fidelity.jl")
    include("performance.jl")
end

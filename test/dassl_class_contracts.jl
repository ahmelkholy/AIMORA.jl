@testset "DASSL-class variable-step public contracts" begin
    @test AIMORA.EMTIntegrationSelection().mode == AIMORA.FixedStepEMT
    @test AIMORA.EMTIntegrationSelection().settings === nothing

    settings = AIMORA.DASSLClassEMTSettings(
        initial_step_s=1.0e-5,
        minimum_step_s=1.0e-9,
        maximum_step_s=1.0e-2,
        maximum_order=5,
        relative_tolerance=1.0e-7,
    )
    selection = AIMORA.EMTIntegrationSelection(settings)
    @test selection.mode == AIMORA.DASSLClassVariableStep
    @test selection.settings === settings
    @test settings.maximum_order == 5
    @test settings.deterministic
    checkpoint_state = [1.0, 2.0]
    checkpoint_derivative = [3.0, 4.0]
    @test AIMORA.dassl_class_checkpoint_boundary!(
        checkpoint_state,
        checkpoint_derivative,
        0.5,
    ) === nothing
    @test checkpoint_state == [1.0, 2.0]
    @test checkpoint_derivative == [3.0, 4.0]

    public_network_request = AIMORA.DASSLClassEMTNetworkRequest(
        "public_owner_inventory",
        settings,
        1,
        0.0,
        1.0e-2,
        ["opaque_public_owner"],
        ((; family=:public_fixture),);
        initial_node_voltage_v=[0.0],
        voltage_absolute_tolerance_v=1.0e-8,
        current_absolute_tolerance_a=1.0e-9,
        voltage_residual_scale_v=1.0e-8,
        current_residual_scale_a=1.0e-9,
        requested_output_times_s=[0.0, 5.0e-3, 1.0e-2],
    )
    @test public_network_request.node_count == 1
    @test public_network_request.owner_identities == ["opaque_public_owner"]
    @test public_network_request.owners isa Vector{Any}
    @test only(public_network_request.owners).family == :public_fixture
    @test public_network_request.requested_output_times_s ==
        [0.0, 5.0e-3, 1.0e-2]
    @test_throws DimensionMismatch AIMORA.DASSLClassEMTNetworkRequest(
        "invalid_owner_inventory",
        settings,
        1,
        0.0,
        1.0e-2,
        ["first", "second"],
        ((; family=:public_fixture),);
        initial_node_voltage_v=[0.0],
        voltage_absolute_tolerance_v=1.0e-8,
        current_absolute_tolerance_a=1.0e-9,
        voltage_residual_scale_v=1.0e-8,
        current_residual_scale_a=1.0e-9,
    )

    @test_throws ArgumentError AIMORA.DASSLClassEMTSettings(
        initial_step_s=1.0e-5,
        minimum_step_s=1.0e-4,
        maximum_step_s=1.0e-2,
    )
    @test_throws ArgumentError AIMORA.DASSLClassEMTSettings(
        initial_step_s=1.0e-5,
        minimum_step_s=1.0e-9,
        maximum_step_s=1.0e-2,
        maximum_order=6,
    )
    @test_throws ArgumentError AIMORA.DASSLClassEMTSettings(
        initial_step_s=1.0e-5,
        minimum_step_s=1.0e-9,
        maximum_step_s=1.0e-2,
        deterministic=false,
    )

    layout = AIMORA.DASSLClassEMTStateLayout(
        [:inductor_current, :capacitor_voltage, :node_constraint],
        [true, true, false],
        ["A", "V", "A"],
        [1.0e-9, 1.0e-8, 1.0e-9],
    )
    @test layout.differential_mask == BitVector([true, true, false])
    @test layout.absolute_tolerances == [1.0e-9, 1.0e-8, 1.0e-9]
    @test_throws ArgumentError AIMORA.DASSLClassEMTStateLayout(
        [:state, :state],
        [true, false],
        ["A", "V"],
        [1.0e-9, 1.0e-8],
    )

    root = AIMORA.DASSLClassEMTValidationRoot(
        :zero_crossing,
        1,
        10,
        (time_s, _state, _derivative) -> time_s - 5.0e-3,
        (_state, _derivative, _time_s) -> nothing,
    )
    task = AIMORA.DASSLClassEMTValidationTask(
        :sample_release,
        5.0e-3,
        20,
        (_state, _derivative, _time_s) -> nothing,
    )
    @test root.direction == 1
    @test task.time_s == 5.0e-3
    @test_throws ArgumentError AIMORA.DASSLClassEMTValidationRoot(
        :bad_direction,
        2,
        1,
        (_time_s, _state, _derivative) -> 0.0,
        (_state, _derivative, _time_s) -> nothing,
    )
    @test_throws ArgumentError AIMORA.DASSLClassEMTValidationTask(
        Symbol(""),
        5.0e-3,
        1,
        (_state, _derivative, _time_s) -> nothing,
    )

    admitted = AIMORA.DASSLClassEMTOwnerDisposition(
        "load_capacitor",
        "AIMORA.Branches.CapacitorBranch",
        true,
        :admitted,
        "complete residual, initialization, tolerance, rollback, and snapshot contract",
    )
    refused = AIMORA.DASSLClassEMTOwnerDisposition(
        "switch_detailed_bridge",
        "AIMORA.Nonlinear.PowerSemiconductorBridgeTopology",
        false,
        :fixed_step_trial_state,
        "switch-detailed bridge state is fixed-step-only in the initial admission matrix",
    )
    readiness = AIMORA.DASSLClassEMTReadiness(
        [refused, admitted],
        3,
        2,
        1,
        "readiness-signature",
    )
    @test !readiness.compatible
    @test getfield.(readiness.owners, :identity) == [
        "load_capacitor",
        "switch_detailed_bridge",
    ]
    @test_throws ArgumentError AIMORA.DASSLClassEMTReadiness(
        [admitted],
        3,
        1,
        1,
        "invalid-counts",
    )

    if !AIMORA.solver_available()
        preparation = AIMORA.prepare_dassl_class_emt((; settings, layout))
        @test preparation isa AIMORA.SolverUnavailableResult
        @test preparation.operation == :prepare_dassl_class_emt
        @test preparation.required_capability == :dassl_class_variable_step_emt
        readiness_unavailable = AIMORA.dassl_class_emt_readiness((; layout))
        @test readiness_unavailable isa AIMORA.SolverUnavailableResult
        @test readiness_unavailable.operation == :dassl_class_emt_readiness
    end
end

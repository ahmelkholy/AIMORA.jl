using AIMORA.NonlinearNetwork

function synthetic_semiconductor_energy_table()
    current_axis_a = [0.0, 10.0]
    voltage_axis_v = [0.0, 100.0]
    temperature_axis_k = [250.0, 350.0]
    values = Array{Float64}(undef, 2, 2, 2)
    for current_index in 1:2, voltage_index in 1:2, temperature_index in 1:2
        values[current_index, voltage_index, temperature_index] =
            1.0e-6 * (
                current_axis_a[current_index] +
                voltage_axis_v[voltage_index] / 10.0 +
                (temperature_axis_k[temperature_index] - 250.0) / 10.0
            )
    end
    return AIMORA.Nonlinear.SwitchingEnergyTable(
        current_axis_a,
        voltage_axis_v,
        temperature_axis_k;
        turn_on_energy_j=values,
        turn_off_energy_j=2.0 .* values,
        reverse_recovery_energy_j=3.0 .* values,
    )
end

@testset "extended power-semiconductor constitutive fidelity" begin
    nonlinear_charge = AIMORA.Nonlinear.NonlinearJunctionChargeFidelity(
        2.0e-6,
        10.0,
        0.4;
        voltage_domain_v=(-100.0, 100.0),
    )
    for voltage_v in (-50.0, -1.0, 0.0, 20.0)
        perturbation_v = 1.0e-6
        central_charge_derivative = (
            AIMORA.Nonlinear.power_semiconductor_junction_charge(
                nonlinear_charge,
                voltage_v + perturbation_v,
            ) -
            AIMORA.Nonlinear.power_semiconductor_junction_charge(
                nonlinear_charge,
                voltage_v - perturbation_v,
            )
        ) / (2.0 * perturbation_v)
        @test central_charge_derivative ≈
            AIMORA.Nonlinear.power_semiconductor_junction_capacitance(
                nonlinear_charge,
                voltage_v,
            ) rtol=1.0e-7
    end

    energy_table = synthetic_semiconductor_energy_table()
    @test AIMORA.Nonlinear.power_semiconductor_switching_energy(
        energy_table,
        :turn_on,
        10.0,
        100.0,
        350.0,
    ) == energy_table.turn_on_energy_j[2, 2, 2]
    @test_throws DomainError AIMORA.Nonlinear.power_semiconductor_switching_energy(
        energy_table,
        :turn_on,
        11.0,
        100.0,
        350.0,
    )

    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        recovered_charge=AIMORA.Nonlinear.RecoveredChargeFidelity(
            2.0e-6;
            initial_charge_c=2.0e-6,
        ),
        junction_charge=nonlinear_charge,
        thermal=AIMORA.Nonlinear.CauerThermalFidelity([1.0, 2.0], [0.5, 1.0]),
    )
    recovered_diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=1.0,
        off_conductance=0.0,
        extended_fidelity=fidelity,
    )
    prepare_nonlinear_device_step!(
        recovered_diode,
        1.0e-6,
        1.0e-6,
        :TrapezoidalCompanion,
    )
    terminal_current_a = zeros(Float64, 2)
    terminal_jacobian_s = zeros(Float64, 2, 2)
    nonlinear_current_jacobian!(
        terminal_current_a,
        terminal_jacobian_s,
        recovered_diode,
        [-1.0, 0.0],
        1.0e-6,
    )
    perturbation_v = 1.0e-7
    positive_current_a = zeros(Float64, 2)
    negative_current_a = zeros(Float64, 2)
    nonlinear_current_jacobian!(
        positive_current_a,
        zeros(Float64, 2, 2),
        recovered_diode,
        [-1.0 + perturbation_v, 0.0],
        1.0e-6,
    )
    nonlinear_current_jacobian!(
        negative_current_a,
        zeros(Float64, 2, 2),
        recovered_diode,
        [-1.0 - perturbation_v, 0.0],
        1.0e-6,
    )
    @test terminal_jacobian_s[1, 1] ≈
        (positive_current_a[1] - negative_current_a[1]) / (2.0 * perturbation_v) rtol=1.0e-7
    fidelity_snapshot = sprint(show, fidelity)
    nonlinear_current_jacobian!(
        terminal_current_a,
        terminal_jacobian_s,
        recovered_diode,
        [-1.0, 0.0],
        1.0e-6,
    )
    @test sprint(show, fidelity) == fidelity_snapshot
    fallback_acceptance_diode = deepcopy(recovered_diode)
    jacobian_acceptance_diode = deepcopy(recovered_diode)
    accept_nonlinear_device_state!(
        fallback_acceptance_diode,
        [-1.0, 0.0],
        terminal_current_a,
        1.0e-6,
    )
    accept_nonlinear_device_state!(
        jacobian_acceptance_diode,
        [-1.0, 0.0],
        terminal_current_a,
        terminal_jacobian_s,
        1.0e-6,
    )
    finish_nonlinear_device_step!(fallback_acceptance_diode)
    finish_nonlinear_device_step!(jacobian_acceptance_diode)
    fallback_acceptance_state = AIMORA.Nonlinear.power_semiconductor_extended_state(
        fallback_acceptance_diode,
    )
    jacobian_acceptance_state = AIMORA.Nonlinear.power_semiconductor_extended_state(
        jacobian_acceptance_diode,
    )
    for field in fieldnames(typeof(fallback_acceptance_state))
        fallback_value = getfield(fallback_acceptance_state, field)
        jacobian_value = getfield(jacobian_acceptance_state, field)
        if fallback_value isa AbstractArray
            @test fallback_value ≈ jacobian_value rtol=16 * eps(Float64) atol=16 * eps(Float64)
        elseif fallback_value isa AbstractFloat
            @test isequal(fallback_value, jacobian_value) ||
                isapprox(
                    fallback_value,
                    jacobian_value;
                    rtol=16 * eps(Float64),
                    atol=16 * eps(Float64),
                )
        else
            @test fallback_value == jacobian_value
        end
    end
    @test fallback_acceptance_diode.last_voltage ≈ jacobian_acceptance_diode.last_voltage
    @test fallback_acceptance_diode.last_current ≈ jacobian_acceptance_diode.last_current
    @test fallback_acceptance_diode.semiconductor_dissipated_energy_j ≈
        jacobian_acceptance_diode.semiconductor_dissipated_energy_j
    fallback_fidelity = something(fallback_acceptance_diode.extended_fidelity)
    @test fallback_fidelity.recovered_charge.stored_charge_c ≈ 2.0e-6 / 3.0
    @test fallback_fidelity.recovered_charge.recovery_active
    @test fallback_fidelity.thermal.node_temperature_k[1] >=
        fallback_fidelity.thermal.ambient_temperature_k

    triac = AIMORA.Nonlinear.TriacSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        threshold_v=0.5,
        holding_current=0.1,
    )
    triac.last_voltage = -1.0
    AIMORA.Nonlinear.apply_power_semiconductor_forward_turn_on!(triac, 0.0)
    @test triac.closed
    @test triac.conduction_direction == -1
    triac.last_forward_current_a = 0.05
    AIMORA.Nonlinear.apply_power_semiconductor_forward_extinction!(triac, 1.0e-6)
    @test !triac.closed
    @test triac.conduction_direction == 0

    refusing_gto = AIMORA.Nonlinear.GateTurnOffThyristorSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        initially_closed=true,
        interruptible_current_a=1.0,
        turn_off_policy=:refuse,
    )
    refusing_gto.last_forward_current_a = 2.0
    refusing_driver_state = sprint(show, refusing_gto.gate_driver)
    @test_throws DomainError AIMORA.Nonlinear.request_power_semiconductor_gate!(
        refusing_gto,
        false,
        0.0,
    )
    @test refusing_gto.closed
    @test sprint(show, refusing_gto.gate_driver) == refusing_driver_state

    remaining_gto = AIMORA.Nonlinear.GateTurnOffThyristorSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(initially_on=true),
        initially_closed=true,
        interruptible_current_a=1.0,
        turn_off_policy=:remain_latched,
    )
    remaining_gto.last_forward_current_a = 2.0
    @test AIMORA.Nonlinear.request_power_semiconductor_gate!(
        remaining_gto,
        false,
        0.0,
    )
    @test remaining_gto.closed
    @test remaining_gto.gate_turn_off_disposition == :remained_latched
end

@testset "private nonlinear ownership and accepted integration" begin
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        junction_charge=AIMORA.Nonlinear.NonlinearJunctionChargeFidelity(
            1.0e-6,
            10.0,
            0.5;
            voltage_domain_v=(-100.0, 100.0),
        ),
    )
    diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=10.0,
        off_conductance=0.01,
        extended_fidelity=fidelity,
    )
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
    )
    nonlinear_system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        1.0e-6,
        1.0e-6,
    )
    @test result.accepted
    @test diode.last_current ≈ 1.0 atol=1.0e-10
    @test !fidelity.candidate_prepared
    accepted_voltage_v = copy(linear_system.v)
    accepted_charge_c = fidelity.junction_charge.previous_charge_c
    # Force the next 0.249 V candidate outside the admitted domain.  A 0.5 V
    # ceiling admits that candidate and therefore cannot exercise rollback.
    fidelity.junction_charge.maximum_voltage_v = 0.2
    rejected_result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-6,
        1.0e-6,
    )
    @test !rejected_result.accepted
    @test linear_system.v == accepted_voltage_v
    @test fidelity.junction_charge.previous_charge_c == accepted_charge_c
    @test_throws ArgumentError AIMORA.NonlinearNodal.NonlinearNodalSystem(
        AIMORA.Nodal.NodalSystem(1, [diode]),
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
end

@testset "localized recovered-charge exhaustion and transactional retry" begin
    recovery = AIMORA.Nonlinear.RecoveredChargeFidelity(
        1.0e-3;
        initial_charge_c=0.5e-3,
    )
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        recovered_charge=recovery,
    )
    diode = AIMORA.Nonlinear.DiodeValveSwitch(
        1,
        0;
        on_conductance=1.0,
        off_conductance=0.1,
        extended_fidelity=fidelity,
    )
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> -1.0)],
    )
    nonlinear_system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [diode];
        scales=NonlinearNetworkScales([1.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        1.0e-3,
        1.0e-3;
        event_policy=AIMORA.OVER16TimestepIntegration.HybridEventPolicy(
            root_time_tolerance_s=1.0e-10,
            root_value_tolerance=1.0e-12,
            simultaneity_tolerance_s=1.0e-10,
        ),
    )
    @test result.accepted
    @test result.diagnostics.discontinuity_reason == :localized_event
    @test result.diagnostics.accepted_substep_count == 2
    @test recovery.recovery_zero_event_count == 1
    @test recovery.last_recovery_zero_time_s ≈ 0.5e-3 atol=1.0e-8
    @test recovery.stored_charge_c == 0.0
    @test fidelity.candidate_prepared == false

    checkpoint = AIMORA.NonlinearNodal.nonlinear_nodal_checkpoint(nonlinear_system)
    accepted_energy_j = diode.semiconductor_dissipated_energy_j
    accepted_event_count = recovery.recovery_zero_event_count
    retry = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-3,
        1.0e-3,
    )
    @test retry.accepted
    @test recovery.recovery_zero_event_count == accepted_event_count
    AIMORA.NonlinearNodal.restore_nonlinear_nodal_checkpoint!(
        nonlinear_system,
        checkpoint,
    )
    replay = AIMORA.NonlinearNodal.advance_nonlinear_step!(
        nonlinear_system,
        2.0e-3,
        1.0e-3,
    )
    @test replay.accepted
    @test recovery.recovery_zero_event_count == accepted_event_count
    @test diode.semiconductor_dissipated_energy_j >= accepted_energy_j
end

@testset "analytic tail cutoff localization" begin
    tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        1.0e-3;
        cutoff_current_a=exp(-1.0),
    )
    tail.active = true
    tail.current_a = 1.0
    tail.initial_current_a = 1.0
    tail.turn_off_time_s = 0.0
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        turn_off_tail=tail,
    )
    igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(),
        off_conductance=0.2,
        extended_fidelity=fidelity,
    )
    system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        AIMORA.Nodal.NodalSystem(
            1,
            [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
        ),
        [igbt];
        scales=NonlinearNetworkScales([10.0], [1.0], Float64[], Float64[]),
    )
    result = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 2.0e-3, 2.0e-3)
    @test result.accepted
    @test result.diagnostics.discontinuity_reason == :localized_event
    @test tail.cutoff_event_count == 1
    @test tail.last_cutoff_time_s ≈ 1.0e-3 atol=1.0e-12
    @test !tail.active
    @test tail.current_a == 0.0

    rounded_tail = AIMORA.Nonlinear.TurnOffTailFidelity(
        1.0e-3;
        cutoff_current_a=exp(-1.0),
    )
    rounded_tail.active = true
    rounded_tail.current_a = 1.0
    rounded_tail.initial_current_a = 1.0
    rounded_tail.turn_off_time_s = 0.0
    rounded_igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        extended_fidelity=AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
            turn_off_tail=rounded_tail,
        ),
    )
    AIMORA.Nonlinear.apply_power_semiconductor_tail_cutoff!(
        rounded_igbt,
        1.0e-3 + 0.5e-12,
    )
    @test rounded_tail.cutoff_event_count == 1
    @test !rounded_tail.active
end

@testset "physical event energy rollback and single retry deposition" begin
    energy_table = synthetic_semiconductor_energy_table()
    fidelity = AIMORA.Nonlinear.PowerSemiconductorExtendedFidelity(
        switching_energy=energy_table,
        thermal=AIMORA.Nonlinear.CauerThermalFidelity(
            [1.0],
            [1.0];
            initial_temperature_k=[300.0],
        ),
    )
    igbt = AIMORA.Nonlinear.IGBTSwitch(
        1,
        0;
        gate_driver=AIMORA.Nonlinear.PowerSemiconductorGateDriver(
            initially_on=true,
        ),
        initially_closed=true,
        on_conductance=1.0,
        off_conductance=0.1,
        extended_fidelity=fidelity,
    )
    AIMORA.Nonlinear.power_semiconductor_event_localization!(igbt)
    igbt.last_voltage = 1.0
    igbt.last_forward_current_a = 1.0
    AIMORA.Nonlinear.request_power_semiconductor_gate!(
        igbt,
        false,
        0.0;
        earliest_transition_time_s=0.4,
    )
    failure_enabled = Ref(true)
    failure = ToggleableFailureCurrentBranch(1, 0, failure_enabled, 0.75, 0)
    linear_system = AIMORA.Nodal.NodalSystem(
        1,
        [AIMORA.Branches.CurrentInjection(1, _time_s -> 1.0)],
    )
    system = AIMORA.NonlinearNodal.NonlinearNodalSystem(
        linear_system,
        [igbt, failure];
        scales=NonlinearNetworkScales([10.0], [1.0], Float64[], Float64[]),
    )
    baseline_voltage_v = copy(linear_system.v)
    baseline_temperature_k = copy(fidelity.thermal.node_temperature_k)
    baseline_energy_j = igbt.semiconductor_dissipated_energy_j
    baseline_transition_count = igbt.topology_transition_count
    baseline_driver_transition_count = igbt.gate_driver.transition_count

    rejected = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 1.0, 1.0)
    @test !rejected.accepted
    @test rejected.failure.code == :nonfinite_device_current
    @test linear_system.v == baseline_voltage_v
    @test igbt.closed
    @test igbt.gate_driver.pending_state === false
    @test igbt.gate_driver.transition_count == baseline_driver_transition_count
    @test igbt.topology_transition_count == baseline_transition_count
    @test energy_table.cumulative_turn_off_energy_j == 0.0
    @test energy_table.last_event_kind == :none
    @test fidelity.thermal.node_temperature_k == baseline_temperature_k
    @test igbt.semiconductor_dissipated_energy_j == baseline_energy_j
    @test failure.accepted_state_count == 0

    failure_enabled[] = false
    accepted = AIMORA.NonlinearNodal.advance_nonlinear_step!(system, 1.0, 1.0)
    @test accepted.accepted
    @test accepted.diagnostics.discontinuity_reason == :localized_event
    @test !igbt.closed
    @test igbt.gate_driver.pending_state === nothing
    @test igbt.gate_driver.transition_count == baseline_driver_transition_count + 1
    @test igbt.topology_transition_count == baseline_transition_count + 1
    expected_event_energy_j = energy_table.last_event_energy_j
    @test expected_event_energy_j > 0.0
    @test energy_table.cumulative_turn_off_energy_j == expected_event_energy_j
    @test energy_table.last_event_kind == :turn_off
    @test energy_table.last_event_energy_j == expected_event_energy_j
    @test igbt.semiconductor_dissipated_energy_j >= expected_event_energy_j
    @test isfinite(fidelity.thermal.node_temperature_k[1])
    @test fidelity.thermal.node_temperature_k[1] != baseline_temperature_k[1]
end

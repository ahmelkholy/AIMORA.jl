using AIMORA.EMTPartitioning
using AIMORA.EMTTaskPlatform

function public_partition_plan(;
    communication_step = 1 // 1_000_000,
    source_step = 1 // 4_000_000,
    load_step = 1 // 1_000_000,
    stop = 20 // 1_000_000,
    method = IteratedWaveformExchange,
)
    regions = (
        EMTPartitionRegion("source", ["source_rl"], emt_logical_time(source_step)),
        EMTPartitionRegion("load", ["load_rc"], emt_logical_time(load_step)),
    )
    port = EMTInterfacePort(
        "interface",
        VoltageCurrentInterfacePort,
        "source",
        "load",
        "source_terminal",
        "load_terminal",
    )
    return emt_partition_plan(
        regions,
        (port,);
        start = emt_logical_time(0),
        stop = emt_logical_time(stop),
        communication_step = emt_logical_time(communication_step),
        exchange = EMTPartitionExchangePolicy(method),
    )
end

@testset "public local multirate partition contracts" begin
    plan = public_partition_plan()
    @test plan.execution_mode == CombinedLocalPartitionedExecution
    @test plan.rate_ratios == (4, 1)
    @test plan.communication_window_count == 20
    @test partition_plan_signature_sha256(plan) == plan.signature_sha256
    @test occursin(r"^[0-9a-f]{64}$", plan.signature_sha256)
    study = PassiveTwoRegionRLCStudy(
        plan;
        source_voltage_v = 120.0,
        source_resistance_ohm = 0.8,
        source_inductance_h = 0.015,
        load_resistance_ohm = 12.0,
        load_capacitance_f = 2.0e-4,
    )
    @test partition_study_signature_sha256(study) == study.signature_sha256
    @test study.signature_sha256 == PassiveTwoRegionRLCStudy(
        plan;
        source_voltage_v = 120.0,
        source_resistance_ohm = 0.8,
        source_inductance_h = 0.015,
        load_resistance_ohm = 12.0,
        load_capacitance_f = 2.0e-4,
    ).signature_sha256

    source_region = EMTDeckRegion(
        "source",
        (
            "source grid SOURCE 1.0e9 0.0 0.0 0.0 120.0",
            "rl source_rl SOURCE SOURCE_TERMINAL 0.8 0.015",
        );
        source_identity = "source_region_deck",
    )
    load_region = EMTDeckRegion(
        "load",
        (
            "resistor load_r LOAD_TERMINAL 0 12.0",
            "capacitor load_c LOAD_TERMINAL 0 2.0e-4",
        );
        source_identity = "load_region_deck",
    )
    deck_plan = emt_partition_plan(
        (
            EMTPartitionRegion(
                "source",
                ["grid", "source_rl"],
                emt_logical_time(1 // 4_000_000),
            ),
            EMTPartitionRegion(
                "load",
                ["load_c", "load_r"],
                emt_logical_time(1 // 1_000_000),
            ),
        ),
        (
            EMTInterfacePort(
                "interface",
                VoltageCurrentInterfacePort,
                "source",
                "load",
                "SOURCE_TERMINAL",
                "LOAD_TERMINAL";
                reference_impedance_ohm = 12.8,
            ),
        );
        start = emt_logical_time(0),
        stop = emt_logical_time(20 // 1_000_000),
        communication_step = emt_logical_time(1 // 1_000_000),
    )
    deck_study = PartitionedDeckEMTStudy(
        deck_plan,
        (load_region, source_region),
    )
    @test getfield.(deck_study.regions, :identity) == ("source", "load")
    @test deck_study.initial_interface_current_a == (0.0,)
    @test partition_study_signature_sha256(deck_study) ==
        deck_study.signature_sha256
    @test occursin(r"^[0-9a-f]{64}$", source_region.signature_sha256)
    @test_throws ArgumentError PartitionedDeckEMTStudy(
        deck_plan,
        (source_region,),
    )

    duplicate_owner = (
        EMTPartitionRegion("one", ["shared_model"], emt_logical_time(1 // 1_000_000)),
        EMTPartitionRegion("two", ["shared_model"], emt_logical_time(1 // 1_000_000)),
    )
    duplicate_port = EMTInterfacePort(
        "duplicate_interface",
        VoltageCurrentInterfacePort,
        "one",
        "two",
        "one_terminal",
        "two_terminal",
    )
    @test_throws ArgumentError emt_partition_plan(
        duplicate_owner,
        (duplicate_port,);
        start = emt_logical_time(0),
        stop = emt_logical_time(1 // 1_000),
        communication_step = emt_logical_time(1 // 1_000_000),
    )
    @test_throws ArgumentError public_partition_plan(
        communication_step = 1 // 1_000_000,
        source_step = 1 // 3_000_001,
    )
    @test_throws ArgumentError public_partition_plan(
        source_step = 1 // 4_000_000,
        load_step = 1 // 1_000_000,
        method = DirectCoupledExchange,
    )

    equal_rate_waveform = public_partition_plan(
        source_step = 1 // 1_000_000,
        load_step = 1 // 1_000_000,
    )
    @test equal_rate_waveform.execution_mode == PartitionedWaveformExecution
    lagged = public_partition_plan(
        source_step = 1 // 1_000_000,
        load_step = 1 // 1_000_000,
        method = LaggedCausalExchange,
    )
    @test lagged.execution_mode == PartitionedLaggedExecution
    local_subcycling = emt_partition_plan(
        (
            EMTPartitionRegion(
                "connected_network",
                ["frequency_dependent_line"],
                emt_logical_time(1 // 4_000_000),
            ),
        ),
        ();
        start=emt_logical_time(0),
        stop=emt_logical_time(4 // 1_000_000),
        communication_step=emt_logical_time(1 // 1_000_000),
        exchange=EMTPartitionExchangePolicy(DirectCoupledExchange),
        execution_mode=LocalSubcyclingExecution,
    )
    @test local_subcycling.execution_mode == LocalSubcyclingExecution

    for kind in (
        VoltageCurrentInterfacePort,
        NortonInterfacePort,
        TheveninInterfacePort,
        ScatteringInterfacePort,
        TravelingWaveInterfacePort,
    )
        equivalent_port = EMTInterfacePort(
            "equivalent_port",
            kind,
            "positive",
            "negative",
            "positive_terminal",
            "negative_terminal";
            reference_impedance_ohm=12.8,
        )
        coordinates = emt_interface_coordinates(equivalent_port, 87.5, -3.25)
        physical = emt_interface_physical_values(
            equivalent_port,
            coordinates.first,
            coordinates.second,
        )
        @test physical.voltage_v ≈ 87.5 rtol=8eps(Float64)
        @test physical.outward_current_a ≈ -3.25 rtol=8eps(Float64)
    end

    if !AIMORA.solver_available()
        unavailable = AIMORA.prepare_partitioned_emt(study)
        @test unavailable isa AIMORA.SolverUnavailableResult
        @test unavailable.operation == :prepare_partitioned_emt
        @test unavailable.required_capability == :local_multirate_partitioned_emt
        unavailable_deck = AIMORA.prepare_partitioned_emt(deck_study)
        @test unavailable_deck isa AIMORA.SolverUnavailableResult
    end
end

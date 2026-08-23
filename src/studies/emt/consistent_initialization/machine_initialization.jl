function _emt_machine_frequency_point(
    sample,
    request::EMTInitializationRequest,
    ;
    frequency_assignment::Symbol=:model_owned_machine_equilibrium,
)
    node_voltage_phasors = ComplexF64.(sample.node_voltage_phasors)
    node_count = length(node_voltage_phasors)
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    isapprox(
        frequency_hz,
        request.frequency_hz;
        atol=1.0e-12,
        rtol=1.0e-12,
    ) || throw(_EMTInitializationRefusal(
        :source_frequency_mismatch,
        :machine_equilibrium,
        :frequency_hz,
        "machine-owned steady-state frequency does not match the initialization request",
        (requested_frequency_hz=request.frequency_hz, machine_frequency_hz=frequency_hz),
    ))
    diagnostics = sample.topology_diagnostics
    topology = _emt_initialization_topology_report(diagnostics, node_count)
    topology.classification === :unique || throw(
        _emt_initialization_classification_failure(topology),
    )
    topology.condition_estimate <= request.tolerances.maximum_condition_estimate ||
        throw(_EMTInitializationRefusal(
            :ill_conditioned,
            :network_topology,
            :nodal_voltage,
            "machine-owned harmonic network exceeds the requested condition limit",
            (
                condition_estimate=topology.condition_estimate,
                maximum_condition_estimate=
                    request.tolerances.maximum_condition_estimate,
            ),
        ))
    admittance = ComplexF64.(sample.steady_state_admittance)
    source_injections = ComplexF64.(
        sample.steady_state_source_injection_phasors,
    )
    size(admittance) == (node_count, node_count) || throw(
        _EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :admittance,
            "machine-owned steady-state sample omitted its complete nodal admittance",
            (; node_count, admittance_size=size(admittance)),
        ),
    )
    length(source_injections) == node_count || throw(
        _EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :source_injection,
            "machine-owned steady-state sample omitted its complete source injection",
            (; node_count, source_injection_count=length(source_injections)),
        ),
    )
    symmetry_error = maximum(abs, admittance - transpose(admittance); init=0.0)
    dissipative = Hermitian(0.5 .* (admittance + adjoint(admittance)))
    minimum_dissipative_eigenvalue = minimum(eigvals(dissipative); init=0.0)
    node_frequencies = hasproperty(sample, :node_steady_state_frequencies_hz) ?
        Float64.(sample.node_steady_state_frequencies_hz) :
        fill(frequency_hz, node_count)
    source_rows = hasproperty(sample, :node_frequency_source_row_indices) ?
        Int.(sample.node_frequency_source_row_indices) : zeros(Int, node_count)
    if length(node_frequencies) < node_count
        append!(
            node_frequencies,
            fill(frequency_hz, node_count - length(node_frequencies)),
        )
    end
    if length(source_rows) < node_count
        append!(source_rows, zeros(Int, node_count - length(source_rows)))
    end
    length(node_frequencies) == node_count && length(source_rows) == node_count ||
        throw(_EMTInitializationRefusal(
            :incomplete_model_state,
            :machine_network_equilibrium,
            :frequency_partition,
            "steady-state frequency ownership does not cover the augmented network",
            (;
                node_count,
                frequency_count=length(node_frequencies),
                source_owner_count=length(source_rows),
            ),
        ))
    successors = hasproperty(sample, :source_frequency_successor_indices) ?
        Int.(sample.source_frequency_successor_indices) : Int[]
    subnetwork_count = hasproperty(sample, :steady_state_frequency_subnetwork_count) ?
        Int(sample.steady_state_frequency_subnetwork_count) : 1
    return EMTInitializationFrequencyPoint(
        :physical_frequency,
        frequency_assignment,
        frequency_hz,
        2.0 * pi * frequency_hz,
        node_frequencies,
        source_rows,
        successors,
        subnetwork_count,
        node_voltage_phasors,
        source_injections,
        zeros(ComplexF64, node_count),
        topology,
        symmetry_error,
        minimum_dissipative_eigenvalue,
        true,
    )
end

function _emt_model_initialization_residual(
    owner::Symbol,
    quantity::Symbol,
    unit::AbstractString,
    value::Real,
    absolute_tolerance::Real,
    relative_tolerance::Real,
    reference_scale::Real,
)
    residual = Float64(value)
    absolute = Float64(absolute_tolerance)
    relative = Float64(relative_tolerance)
    scale = Float64(reference_scale)
    allowance = absolute + relative * scale
    scaled = allowance > 0.0 ? residual / allowance :
        (residual == 0.0 ? 0.0 : Inf)
    return EMTInitializationResidual(
        owner,
        quantity,
        String(unit),
        residual,
        absolute,
        relative,
        scale,
        0.0,
        scaled,
        isfinite(residual) && residual <= allowance,
    )
end

function _emt_machine_network_residual(
    point::EMTInitializationFrequencyPoint,
    request::EMTInitializationRequest,
)
    tolerance = request.tolerances.current_relative
    return _emt_model_initialization_residual(
        :machine_network_equilibrium,
        :scaled_nodal_backward_error,
        "1",
        point.topology.relative_residual,
        tolerance,
        0.0,
        1.0,
    )
end

function _emt_validate_machine_initialization_request(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    request.formulation isa PhysicalFrequencyFormulation || throw(
        _EMTInitializationRefusal(
            :unsupported_formulation,
            :machine_equilibrium,
            :harmonic_formulation,
            "the admitted machine initializer owns a physical-frequency equilibrium; timestep-matched machine recurrence has not been proven",
            (requested_formulation=_emt_harmonic_formulation_symbol(request.formulation),),
        ),
    )
    all(frequency -> isapprox(
        frequency,
        request.frequency_hz;
        atol=1.0e-12,
        rtol=1.0e-12,
    ), request.frequency_grid_hz) || throw(_EMTInitializationRefusal(
        :unsupported_frequency_scan,
        :machine_equilibrium,
        :frequency_grid_hz,
        "model-owned machine initialization accepts one operating frequency; whole-network scan points must be requested from the network scan owner",
        (frequency_grid_hz=copy(request.frequency_grid_hz),),
    ))
    request.operating_point === nothing || throw(_EMTInitializationRefusal(
        :unsupported_operating_point_mapping,
        :machine_equilibrium,
        :operating_point,
        "machine operating-point import requires machine-owned current, torque, flux, and control quantities rather than voltage-only constraints",
        (operating_point_type=string(typeof(request.operating_point)),),
    ))
    request.time_origin_s == 0.0 || throw(_EMTInitializationRefusal(
        :unsupported_time_origin,
        :machine_equilibrium,
        :time_origin_s,
        "machine-owned initialization currently requires the deck time origin",
        (time_origin_s=request.time_origin_s,),
    ))
    _validate_emt_initialization_source_frequency(parsed, request)
    return request
end

function _emt_machine_state_inventory(
    parsed::DeckParser.DeckParseResult,
    machine_family::Symbol,
    electrical_state_count::Int,
    mechanical_state_count::Int,
)
    records = EMTInitializationStateRecord[]
    _append_emt_initialization_state!(
        records,
        :network_voltage,
        :algebraic,
        maximum(values(parsed.node_map); init=0),
        :model_owned_harmonic_network_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :network_topology,
        :discrete,
        1,
        :ranked_connected_component_classification,
    )
    _append_emt_initialization_state!(
        records,
        Symbol(machine_family, :_electrical_state),
        :continuous,
        electrical_state_count,
        :model_owned_machine_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        Symbol(machine_family, :_mechanical_state),
        :continuous,
        mechanical_state_count,
        :model_owned_torque_speed_equilibrium,
    )
    switch_count = length(DeckParser.deck_over5_switch_rows(parsed)) +
        length(DeckParser.deck_control_system_switch_coupling_rows(parsed))
    _append_emt_initialization_state!(
        records,
        :switch_mode,
        :discrete,
        switch_count,
        :declared_initial_topology,
    )
    _append_emt_initialization_state!(
        records,
        :switch_event_state,
        :scheduler,
        length(DeckParser.deck_over5_switch_rows(parsed)),
        :initial_event_surface_classification,
    )
    control_count = length(
        DeckParser.deck_synchronous_machine_control_interface_rows(parsed),
    )
    _append_emt_initialization_state!(
        records,
        :machine_control_state,
        :discrete,
        control_count,
        :declared_control_equilibrium,
    )
    transformer_count = get(
        parsed.card_counts,
        :fixed_card_saturated_transformer_intake,
        0,
    )
    _append_emt_initialization_state!(
        records,
        :machine_terminal_transformer_state,
        :continuous,
        transformer_count,
        :coupled_transformer_branch_equilibrium,
    )
    _append_emt_initialization_state!(
        records,
        :output_cursor,
        :discrete,
        1,
        :time_zero_output_epoch,
    )
    sort!(records; by=record -> (String(record.state_family), String(record.owner)))
    return records
end

function _emt_synchronous_machine_preparation(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    machine_indices = sort!(unique(
        row.machine_index for row in
        DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
    ))
    length(machine_indices) == 1 || throw(_EMTInitializationRefusal(
        :unsupported_machine_fleet_initialization,
        :synchronous_machine_state,
        :machine_count,
        "the cohesive synchronous-machine initializer currently requires one admitted machine",
        (machine_indices=copy(machine_indices),),
    ))
    machine_index = only(machine_indices)
    timestep_s = DeckParser.deck_fixed_time_horizon_options(parsed).dt_s
    context = _deck_synchronous_machine_runtime_context(
        parsed,
        timestep_s,
        timestep_s;
        saturated_transformer_branch_runtime_enabled=true,
        coupled_lumped_sequence_history_enabled=true,
        recorded_step_indices=[0],
    )
    sample = _deck_synchronous_machine_network_initial_sample(
        parsed,
        context;
        strict_topology_classification=true,
    )
    sample === nothing && throw(_EMTInitializationRefusal(
        :missing_model_state,
        :synchronous_machine_state,
        :terminal_voltage,
        "synchronous-machine terminal equilibrium is missing",
        (machine_index=machine_index,),
    ))
    point = _emt_machine_frequency_point(sample, request)
    initialization = _deck_synchronous_machine_initial_state(
        parsed,
        context,
        sample;
        machine_index,
    )
    horizon = run_deck_synchronous_machine_horizon(
        parsed,
        deepcopy(initialization.state);
        numask=initialization.numask,
        nlocg=initialization.nlocg,
        nloce=initialization.nloce,
        time_step_s=timestep_s,
        dynamic_step_count=1,
        angle_half_step_inverse=initialization.angle_half_step_inverse,
        speed_tolerance=initialization.speed_tolerance,
        omega_tolerance=initialization.omega_tolerance,
        speed_floor=initialization.speed_floor,
        max_iterations=initialization.max_iterations,
        damping_ratio=initialization.damping_ratio,
        rotor_angle_extrapolation_interval=
            initialization.rotor_angle_extrapolation_interval,
        speed_voltage_factor=initialization.speed_voltage_factor,
        electrical_speed_rad_s=initialization.electrical_speed_rad_s,
        electrical_angle_increment=initialization.electrical_angle_increment,
        saturated_transformer_branch_runtime_enabled=true,
        coupled_lumped_sequence_history_enabled=true,
        recorded_step_indices=[0, 1],
    )
    terminal_rows = sort!(
        [
            row for row in DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)
            if row.machine_index == machine_index
        ];
        by=row -> row.phase_index,
    )
    terminal_nodes = Int[row.terminal_node_value for row in terminal_rows]
    tolerances = request.tolerances
    current_error = maximum(
        abs,
        initialization.phase_current_phasors .-
            sample.node_current_phasors[terminal_nodes];
        init=0.0,
    )
    current_scale = maximum(abs, initialization.phase_current_phasors; init=0.0)
    voltage_error = maximum(
        abs,
        horizon.terminal_voltage_values[:, 1] .-
            real.(sample.node_voltage_phasors[terminal_nodes]);
        init=0.0,
    )
    voltage_scale = maximum(
        abs,
        sample.node_voltage_phasors[terminal_nodes];
        init=0.0,
    )
    terminal_kcl = abs(sum(horizon.terminal_current_values[:, 1]))
    residuals = EMTInitializationResidual[
        _emt_machine_network_residual(point, request),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_current_equilibrium,
            "A",
            current_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
        ),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_voltage_equilibrium,
            "V",
            voltage_error,
            tolerances.voltage_absolute_v,
            tolerances.voltage_relative,
            voltage_scale,
        ),
        _emt_model_initialization_residual(
            :synchronous_machine_state,
            :terminal_zero_sequence_kcl,
            "A",
            terminal_kcl,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            maximum(abs, horizon.terminal_current_values[:, 1]; init=0.0),
        ),
    ]
    maximum_scaled = maximum(residual.scaled_value for residual in residuals)
    metric = NoArtificialTransientMetric(
        :synchronous_machine_time_zero_equilibrium,
        "1",
        request.time_origin_s,
        request.time_origin_s,
        maximum_scaled <= 1.0 ? 0.0 : maximum_scaled - 1.0,
        maximum_scaled <= 1.0 ? 0.0 : maximum_scaled - 1.0,
        0.0,
        0.0,
        request.tolerances.no_artificial_transient_normalized_rms,
        all(residual -> residual.passed, residuals),
    )
    accepted_state = (
        machine_initialization=deepcopy(initialization),
        network_voltage_phasors=copy(sample.node_voltage_phasors),
        terminal_voltage_values=copy(horizon.terminal_voltage_values[:, 1]),
        terminal_current_values=copy(horizon.terminal_current_values[:, 1]),
        machine_output_values=copy(horizon.machine_output_values[:, 1]),
        mechanical_history_values=copy(horizon.mechanical_history_values[:, 1]),
        control_output_names=copy(horizon.control_output_names),
        control_output_values=copy(horizon.control_output_values[:, 1]),
        switch_node_groups=copy.(point.topology.switch_node_groups),
        output_cursor=0,
    )
    prepared = PreparedMachineEMTStudy(
        :synchronous_machine,
        accepted_state,
        horizon,
        parsed,
    )
    inventory = _emt_machine_state_inventory(
        parsed,
        :synchronous_machine,
        length(initialization.state.current_history),
        length(initialization.state.equation_state.histq_values),
    )
    warnings = sample.time_zero_ground_fault ? String[
        "The first advance contains the deck-declared switch transition; no-artificial-transient acceptance is therefore evaluated at the complete time-zero machine/network equilibrium.",
    ] : String[]
    return prepared, point, residuals, inventory, metric, warnings
end

function _emt_universal_machine_preparation(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    machine_indices = sort!(unique(
        row.machine_index for row in
        DeckParser.deck_universal_machine_definition_rows(parsed)
        if row.card_index == 1
    ))
    length(machine_indices) == 1 || throw(_EMTInitializationRefusal(
        :unsupported_machine_fleet_initialization,
        :universal_machine_state,
        :machine_count,
        "the cohesive universal-machine initializer currently requires one admitted machine",
        (machine_indices=copy(machine_indices),),
    ))
    machine_index = only(machine_indices)
    card = _deck_universal_machine_definition(parsed, machine_index, 1)
    card.machine_type in 3:12 || throw(_EMTInitializationRefusal(
        :unsupported_machine_family,
        :universal_machine_state,
        :machine_type,
        "universal wound-field synchronous types use a separate admitted initialization owner",
        (machine_type=card.machine_type,),
    ))
    automatic_direct = card.machine_type in 8:12 &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        deck_direct_current_machine_automatic_initialization(
            parsed;
            machine_index,
        ) : nothing
    single_phase = card.machine_type in (6, 7) &&
        _deck_universal_machine_initialization_mode(parsed) == :automatic ?
        _deck_single_phase_induction_initialization(parsed, machine_index) : nothing
    steady_state = automatic_direct !== nothing ? automatic_direct.steady_state :
        single_phase !== nothing ? single_phase.steady_state :
        deck_steady_state_voltage_phasors(parsed)
    seed_state = automatic_direct !== nothing ? automatic_direct.state :
        single_phase !== nothing ? single_phase.state :
        card.machine_type in 3:7 ?
        deck_induction_machine_initial_state(
            parsed;
            machine_index,
            steady_state,
        ) :
        deck_direct_current_machine_initial_state(
            parsed;
            machine_index,
            steady_state,
        )
    point = _emt_machine_frequency_point(steady_state, request)
    horizon = run_deck_universal_machine_horizon(
        parsed;
        machine_index,
        dynamic_step_count=1,
    )
    time_zero_state = InductionMachineState(
        copy(horizon.current_values[:, 1]),
        copy(horizon.history_currents[:, 1]);
        mechanical_speed_rad_s=horizon.mechanical_speed_rad_s[1],
        previous_mechanical_speed_rad_s=horizon.mechanical_speed_rad_s[1],
        mechanical_angle_rad=horizon.mechanical_angle_rad[1],
    )
    time_zero_state.d_axis_flux = horizon.d_axis_flux[1]
    time_zero_state.q_axis_flux = horizon.q_axis_flux[1]
    time_zero_state.generated_torque = horizon.generated_torque[1]
    time_zero_state.output_values .= horizon.output_values[:, 1]
    time_zero_state.call_count = 1
    tolerances = request.tolerances
    current_seed_error = maximum(
        abs,
        seed_state.current_values .- horizon.current_values[:, 1];
        init=0.0,
    )
    current_scale = maximum(abs, horizon.current_values[:, 1]; init=0.0)
    residuals = EMTInitializationResidual[
        _emt_machine_network_residual(point, request),
        _emt_model_initialization_residual(
            card.machine_type in 3:7 ? :induction_machine_state :
                :direct_current_machine_state,
            :time_zero_current_equilibrium,
            "A",
            current_seed_error,
            tolerances.current_absolute_a,
            tolerances.current_relative,
            current_scale,
        ),
        _emt_model_initialization_residual(
            card.machine_type in 3:7 ? :induction_machine_state :
                :direct_current_machine_state,
            :coupled_runtime_completion,
            "count",
            horizon.complete_induction_machine_path ? 0.0 : 1.0,
            eps(Float64),
            0.0,
            1.0,
        ),
    ]
    if automatic_direct !== nothing &&
       hasproperty(automatic_direct, :armature_kvl_residual)
        push!(
            residuals,
            _emt_model_initialization_residual(
                :direct_current_machine_state,
                :armature_kvl,
                "V",
                abs(Float64(automatic_direct.armature_kvl_residual)),
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(Float64(automatic_direct.requested_armature_voltage)),
            ),
        )
        push!(
            residuals,
            _emt_model_initialization_residual(
                :direct_current_machine_state,
                :field_kvl,
                "V",
                abs(Float64(automatic_direct.field_kvl_residual)),
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                abs(Float64(automatic_direct.field_voltage)),
            ),
        )
    end
    current_envelope_drift = abs(
        norm(horizon.current_values[:, 2]) -
            norm(horizon.current_values[:, 1]),
    ) / max(norm(horizon.current_values[:, 1]), tolerances.current_absolute_a)
    initial_flux_envelope = hypot(horizon.d_axis_flux[1], horizon.q_axis_flux[1])
    flux_envelope_drift = abs(
        hypot(horizon.d_axis_flux[2], horizon.q_axis_flux[2]) -
            initial_flux_envelope,
    ) / max(initial_flux_envelope, tolerances.flux_absolute_wb)
    torque_envelope_drift = abs(
        horizon.generated_torque[2] - horizon.generated_torque[1],
    ) / max(abs(horizon.generated_torque[1]), tolerances.power_absolute_w)
    envelope_drift = maximum((
        current_envelope_drift,
        flux_envelope_drift,
        torque_envelope_drift,
    ))
    frequency_step = pi * request.frequency_hz *
        DeckParser.deck_fixed_time_horizon_options(parsed).dt_s
    physical_to_trapezoidal_warping = frequency_step == 0.0 ? 0.0 :
        abs(tan(frequency_step) / frequency_step - 1.0)
    excess_envelope_drift = max(
        0.0,
        envelope_drift - physical_to_trapezoidal_warping,
    )
    transient_threshold = request.tolerances.no_artificial_transient_normalized_rms
    metric = NoArtificialTransientMetric(
        :machine_periodic_envelope,
        "1",
        request.time_origin_s,
        request.time_origin_s +
            DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
        excess_envelope_drift,
        excess_envelope_drift,
        envelope_drift,
        0.0,
        transient_threshold,
        isfinite(excess_envelope_drift) &&
            excess_envelope_drift <= transient_threshold,
    )
    accepted_state = (
        machine_state=time_zero_state,
        network_voltage_phasors=copy(steady_state.node_voltage_phasors),
        compensated_voltage_values=copy(horizon.compensated_voltage_values[:, 1]),
        power_terminal_voltage_values=copy(horizon.power_terminal_voltages[:, 1]),
        current_substitution_values=copy(horizon.current_substitution_values[:, 1]),
        drive_source_value=horizon.drive_source_values[1],
        excitation_source_value=horizon.excitation_source_values[1],
        report_output_names=copy(horizon.report_output_names),
        report_output_values=copy(horizon.report_output_values[:, 1]),
        switch_node_groups=copy.(point.topology.switch_node_groups),
        output_cursor=0,
    )
    machine_family = card.machine_type in 3:7 ? :induction_machine :
        :direct_current_machine
    prepared = PreparedMachineEMTStudy(
        machine_family,
        accepted_state,
        horizon,
        parsed,
    )
    inventory = _emt_machine_state_inventory(
        parsed,
        machine_family,
        length(time_zero_state.current_values) +
            length(time_zero_state.history_currents) + 2,
        3,
    )
    return prepared, point, residuals, inventory, metric, String[]
end

function _initialize_model_owned_machine_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
)
    points = EMTInitializationFrequencyPoint[]
    residuals = EMTInitializationResidual[]
    state_inventory = EMTInitializationStateRecord[]
    transient_metrics = NoArtificialTransientMetric[]
    try
        working_parsed = deepcopy(parsed)
        _emt_validate_machine_initialization_request(working_parsed, request)
        detailed_synchronous = !isempty(
            DeckParser.deck_synchronous_machine_terminal_voltage_rows(working_parsed),
        )
        universal = !isempty(
            DeckParser.deck_universal_machine_definition_rows(working_parsed),
        )
        detailed_synchronous != universal || throw(_EMTInitializationRefusal(
            :ambiguous_machine_owner,
            :machine_equilibrium,
            :machine_family,
            "one cohesive initialization request must resolve to exactly one machine owner",
            (detailed_synchronous, universal),
        ))
        prepared, point, residuals, state_inventory, metric, warnings =
            detailed_synchronous ?
            _emt_synchronous_machine_preparation(working_parsed, request) :
            _emt_universal_machine_preparation(working_parsed, request)
        push!(points, point)
        push!(transient_metrics, metric)
        all(residual -> residual.passed, residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :machine_equilibrium,
                :scaled_residual,
                "one or more machine initialization residuals exceeded their quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, residuals),),
            ),
        )
        metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :machine_equilibrium,
            metric.quantity,
            "machine initialization exceeded its physical-frequency first-step envelope allowance",
            (
                normalized_rms=metric.normalized_rms,
                envelope_drift=metric.low_frequency_envelope_drift,
                threshold=metric.threshold,
            ),
        ))
        initialized_state_owners = _emt_initialized_state_owners(state_inventory)
        mappings = OperatingPointMappingRecord[]
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            :physical_frequency,
            request.frequency_hz,
            request.time_origin_s,
            point.topology,
            points,
            residuals,
            mappings,
            state_inventory,
            initialized_state_owners,
            Symbol[],
            transient_metrics,
            warnings,
            request.project_signature,
            request.settings_signature,
            request.model_signature,
            signature,
        )
        return EMTInitializationResult(prepared, report, nothing)
    catch error
        refusal = error isa _EMTInitializationRefusal ? error :
            _EMTInitializationRefusal(
                :initialization_error,
                :machine_equilibrium,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            residuals,
            state_inventory,
            transient_metrics,
        )
    end
end

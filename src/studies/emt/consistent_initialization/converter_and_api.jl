function _average_value_converter_current_source_phasors(
    parsed::DeckParser.DeckParseResult,
    converter::AverageValueGridFollowingConverterInitialization,
    frequency_hz::Float64,
)
    phasors = ComplexF64[]
    for owner in converter.current_source_names
        element_index = findfirst(==(owner), parsed.element_names)
        element_index === nothing && throw(_EMTInitializationRefusal(
            :missing_model_input,
            :average_value_grid_following_converter,
            :phase_current_source,
            "average-value converter phase-current owner $(owner) is missing from the EMT network",
            (current_source=owner,),
        ))
        element = parsed.elements[element_index]
        element isa CurrentInjection &&
            element.value isa SinusoidalSourceSignal || throw(
                _EMTInitializationRefusal(
                    :unsupported_model_input,
                    :average_value_grid_following_converter,
                    :phase_current_source,
                    "average-value converter phase-current owners must be typed sinusoidal current injections",
                    (current_source=owner, element_type=string(typeof(element))),
                ),
            )
        push!(
            phasors,
            sinusoidal_source_peak_phasor(element.value, frequency_hz),
        )
    end
    return Tuple(phasors)
end

function _average_value_converter_terminal_phasors(
    parsed::DeckParser.DeckParseResult,
    point::EMTInitializationFrequencyPoint,
    converter::AverageValueGridFollowingConverterInitialization,
)
    phasors = ComplexF64[]
    for terminal in converter.terminal_nodes
        node = get(parsed.node_map, terminal, 0)
        node > 0 || throw(_EMTInitializationRefusal(
            :missing_model_input,
            :average_value_grid_following_converter,
            :phase_terminal,
            "average-value converter phase terminal $(terminal) is missing from the EMT network",
            (terminal,),
        ))
        push!(phasors, point.node_voltage_phasors[node])
    end
    return Tuple(phasors)
end

function _average_value_converter_state_values(state::InverterState)
    return (state.id_a, state.iq_a, state.xid, state.xiq)
end

function initialize_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest,
    converter::AverageValueGridFollowingConverterInitialization;
    kwargs...,
)
    network_result = initialize_emt_study(parsed, request; kwargs...)
    initialization_accepted(network_result) || return network_result
    points = copy(network_result.report.frequency_scan)
    mappings = copy(network_result.report.mappings)
    network_residuals = copy(network_result.report.residuals)
    network_metrics = copy(network_result.report.transient_metrics)
    try
        parameters = converter.parameters
        isapprox(
            parameters.f_hz,
            request.frequency_hz;
            atol=1.0e-12,
            rtol=1.0e-12,
        ) || throw(_EMTInitializationRefusal(
            :frequency_mapping_mismatch,
            :average_value_grid_following_converter,
            :frequency_hz,
            "average-value converter frequency must equal the requested EMT operating frequency",
            (converter_frequency_hz=parameters.f_hz,
             request_frequency_hz=request.frequency_hz),
        ))
        primary = _emt_initialization_primary_point(points, request.frequency_hz)
        timestep_s = network_result.prepared.runtime_template.context.dt_s
        equilibrium = try
            grid_following_inverter_equilibrium(
                parameters;
                timestep_s,
            )
        catch error
            throw(_EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :dq_equilibrium,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            ))
        end
        terminal_phasors = _average_value_converter_terminal_phasors(
            parsed,
            primary,
            converter,
        )
        expected_voltage_magnitude_v =
            sqrt(2.0) * parameters.v_ll_rms_v / sqrt(3.0)
        phase_rotations = (
            1.0 + 0.0im,
            cis(-2.0 * pi / 3.0),
            cis(2.0 * pi / 3.0),
        )
        expected_terminal_phasors = ntuple(
            phase -> expected_voltage_magnitude_v * phase_rotations[phase],
            3,
        )
        terminal_voltage_error_v = maximum(
            abs(terminal_phasors[phase] - expected_terminal_phasors[phase])
            for phase in 1:3
        )
        voltage_allowance_v = request.tolerances.voltage_absolute_v +
            request.tolerances.voltage_relative * expected_voltage_magnitude_v
        terminal_voltage_error_v <= voltage_allowance_v || throw(
            _EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :terminal_voltage_phasors,
                "converter terminal phasors do not match the declared balanced grid-voltage basis",
                (error_v=terminal_voltage_error_v,
                 allowance_v=voltage_allowance_v),
            ),
        )
        source_phasors = _average_value_converter_current_source_phasors(
            parsed,
            converter,
            request.frequency_hz,
        )
        current_injection_error_a = maximum(
            abs(
                source_phasors[phase] -
                equilibrium.phase_current_phasors_a[phase],
            )
            for phase in 1:3
        )
        current_scale_a = max(
            maximum(abs, equilibrium.phase_current_phasors_a),
            request.tolerances.current_absolute_a,
        )
        current_allowance_a = request.tolerances.current_absolute_a +
            request.tolerances.current_relative * current_scale_a
        current_injection_error_a <= current_allowance_a || throw(
            _EMTInitializationRefusal(
                :infeasible_model_state,
                :average_value_grid_following_converter,
                :terminal_current_phasors,
                "network current-source phasors do not match the model-owned dq current equilibrium",
                (error_a=current_injection_error_a,
                 allowance_a=current_allowance_a),
            ),
        )
        derivative = equilibrium.derivative
        current_derivative_a_per_s = max(abs(derivative.id_a), abs(derivative.iq_a))
        controller_error_a = max(abs(derivative.xid), abs(derivative.xiq))
        initial_values = _average_value_converter_state_values(equilibrium.state)
        advanced_values = _average_value_converter_state_values(
            equilibrium.one_step_state,
        )
        current_recurrence_error_a = max(
            abs(advanced_values[1] - initial_values[1]),
            abs(advanced_values[2] - initial_values[2]),
        )
        integral_recurrence_error_as = max(
            abs(advanced_values[3] - initial_values[3]),
            abs(advanced_values[4] - initial_values[4]),
        )
        power_row = inverter_row(equilibrium.state, 0.0, parameters)
        active_power_error_w = abs(
            power_row[6] * parameters.s_base_va - equilibrium.active_power_w,
        )
        reactive_power_error_var = abs(
            power_row[7] * parameters.s_base_va - equilibrium.reactive_power_var,
        )
        tolerances = request.tolerances
        converter_residuals = EMTInitializationResidual[
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :balanced_terminal_voltage_mapping,
                "V",
                terminal_voltage_error_v,
                tolerances.voltage_absolute_v,
                tolerances.voltage_relative,
                expected_voltage_magnitude_v,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :terminal_current_injection,
                "A",
                current_injection_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :dq_current_derivative,
                "A/s",
                current_derivative_a_per_s,
                tolerances.current_absolute_a / timestep_s,
                tolerances.current_relative,
                current_scale_a / timestep_s,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :controller_current_error,
                "A",
                controller_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :one_step_current_recurrence,
                "A",
                current_recurrence_error_a,
                tolerances.current_absolute_a,
                tolerances.current_relative,
                current_scale_a,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :one_step_controller_recurrence,
                "A*s",
                integral_recurrence_error_as,
                tolerances.current_absolute_a * timestep_s,
                tolerances.current_relative,
                current_scale_a * timestep_s,
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :active_power_equilibrium,
                "W",
                active_power_error_w,
                tolerances.power_absolute_w,
                tolerances.power_relative,
                abs(equilibrium.active_power_w),
            ),
            _emt_model_initialization_residual(
                :average_value_grid_following_converter,
                :reactive_power_equilibrium,
                "var",
                reactive_power_error_var,
                tolerances.power_absolute_w,
                tolerances.power_relative,
                abs(equilibrium.reactive_power_var),
            ),
        ]
        all(residual -> residual.passed, converter_residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :average_value_grid_following_converter,
                :scaled_residual,
                "average-value converter equilibrium exceeded one or more quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, converter_residuals),),
            ),
        )
        state_inventory = copy(network_result.report.state_inventory)
        _append_emt_initialization_state!(
            state_inventory,
            :converter_filter_current_state,
            :continuous,
            2,
            :balanced_dq_zero_derivative_equilibrium,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_current_controller_state,
            :continuous,
            2,
            :zero_error_integral_bias,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_terminal_current_injection,
            :algebraic,
            3,
            :balanced_positive_sequence_current_phasors,
        )
        _append_emt_initialization_state!(
            state_inventory,
            :converter_dc_voltage_input,
            :algebraic,
            1,
            :declared_constant_dc_boundary,
        )
        sort!(state_inventory; by=record -> (
            String(record.state_family),
            String(record.owner),
        ))
        initial_state_scale = max(
            maximum(abs, initial_values),
            tolerances.current_absolute_a,
        )
        state_recurrence_error = maximum(
            abs(advanced_values[index] - initial_values[index])
            for index in eachindex(initial_values)
        )
        normalized_state_recurrence = state_recurrence_error / initial_state_scale
        converter_metric = NoArtificialTransientMetric(
            :average_value_converter_state,
            "scaled state",
            request.time_origin_s,
            request.time_origin_s + timestep_s,
            normalized_state_recurrence,
            normalized_state_recurrence,
            normalized_state_recurrence,
            0.0,
            tolerances.no_artificial_transient_normalized_rms,
            normalized_state_recurrence <=
                tolerances.no_artificial_transient_normalized_rms,
        )
        converter_metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :average_value_grid_following_converter,
            :state_recurrence,
            "average-value converter state changed during its undisturbed first-step probe",
            (normalized_state_recurrence,
             threshold=tolerances.no_artificial_transient_normalized_rms),
        ))
        prepared = PreparedAverageValueGridFollowingEMTStudy(
            network_result.prepared,
            equilibrium,
            converter,
        )
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            network_result.report.formulation,
            request.frequency_hz,
            request.time_origin_s,
            network_result.report.topology,
            points,
            [network_residuals; converter_residuals],
            mappings,
            state_inventory,
            _emt_initialized_state_owners(state_inventory),
            Symbol[],
            [network_metrics; converter_metric],
            vcat(
                network_result.report.warnings,
                [
                    "The admitted average-value converter has a constant DC-voltage boundary and no switching, PLL, grid-forming, or sampled-control state.",
                ],
            ),
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
                :average_value_grid_following_converter,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            mappings,
            residuals=network_residuals,
            transient_metrics=network_metrics,
        )
    end
end

function _emt_initialization_failure_result(
    request::EMTInitializationRequest,
    refusal::_EMTInitializationRefusal;
    points::Vector{EMTInitializationFrequencyPoint}=
        EMTInitializationFrequencyPoint[],
    mappings::Vector{OperatingPointMappingRecord}=
        OperatingPointMappingRecord[],
    residuals::Vector{EMTInitializationResidual}=
        EMTInitializationResidual[],
    initialized_state_owners::Vector{Symbol}=Symbol[],
    state_inventory::Vector{EMTInitializationStateRecord}=
        EMTInitializationStateRecord[],
    unsupported_state_owners::Vector{Symbol}=Symbol[],
    transient_metrics::Vector{NoArtificialTransientMetric}=
        NoArtificialTransientMetric[],
)
    topology = isempty(points) ?
        _empty_emt_initialization_topology_report() :
        last(points).topology
    report = EMTInitializationReport(
        :failed,
        _emt_harmonic_formulation_symbol(request.formulation),
        request.frequency_hz,
        request.time_origin_s,
        topology,
        points,
        residuals,
        mappings,
        state_inventory,
        initialized_state_owners,
        unsupported_state_owners,
        transient_metrics,
        String[],
        request.project_signature,
        request.settings_signature,
        request.model_signature,
        "",
    )
    failure = EMTInitializationFailure(
        refusal.code,
        refusal.owner,
        refusal.quantity,
        refusal.message,
        refusal.context,
    )
    return EMTInitializationResult(nothing, report, failure)
end

function initialize_emt_study(
    parsed::DeckParser.DeckParseResult,
    request::EMTInitializationRequest;
    timestep_s::Union{Nothing,Real}=nothing,
    t_end_s::Real=timestep_s !== nothing ? Float64(timestep_s) :
        request.formulation isa TimestepMatchedFormulation ?
        request.formulation.timestep_s :
        DeckParser.deck_fixed_time_horizon_options(parsed).dt_s,
    recorded_step_indices=nothing,
    output_schedule::Symbol=:all_steps,
    source_signal_provider::AbstractSourceSignalProvider=
        IdentitySourceSignalProvider(),
)
    points = EMTInitializationFrequencyPoint[]
    mappings = OperatingPointMappingRecord[]
    residuals = EMTInitializationResidual[]
    initialized_state_owners = Symbol[]
    state_inventory = EMTInitializationStateRecord[]
    unsupported_state_owners = Symbol[]
    transient_metrics = NoArtificialTransientMetric[]
    fixed_source_load_flow = nothing
    try
        DeckParser.assert_deck_valid!(parsed)
        if !isempty(DeckParser.deck_synchronous_machine_terminal_voltage_rows(parsed)) ||
           !isempty(DeckParser.deck_universal_machine_definition_rows(parsed))
            return _initialize_model_owned_machine_emt_study(parsed, request)
        end
        working_parsed = deepcopy(parsed)
        saturated_transformer_declared = get(
            working_parsed.card_counts,
            :fixed_card_saturated_transformer_intake,
            0,
        ) > 0
        saturated_transformer_intake =
            _deck_runtime_saturated_transformer_intake(working_parsed)
        saturated_transformer_declared &&
            saturated_transformer_intake === nothing && throw(
                _EMTInitializationRefusal(
                    :missing_source_backed_intake,
                    :saturated_transformer_magnetic_state,
                    :source_path,
                    "saturated-transformer initialization requires its exact source-backed characteristic and winding intake",
                    (source=working_parsed.source,),
                ),
            )
        unsupported_state_owners =
            _emt_unsupported_initialization_owners(working_parsed)
        isempty(unsupported_state_owners) || throw(_EMTInitializationRefusal(
            :unsupported_state_owner,
            :model_state,
            first(unsupported_state_owners),
            "the requested model requires unsupported initialization state: " *
            join(String.(unsupported_state_owners), ", "),
            (unsupported=copy(unsupported_state_owners),),
        ))
        if !isempty(DeckParser.deck_fixed_source_constraint_rows(working_parsed))
            applied = try
                apply_deck_fixed_source_load_flow(
                    working_parsed;
                    relative_power_tolerance=
                        _fixed_source_normalized_power_tolerance(
                            working_parsed,
                            request.tolerances.power_absolute_w,
                            request.tolerances.power_relative,
                        ),
                    maximum_iterations=
                        request.tolerances.operating_point_maximum_iterations,
                )
            catch error
                throw(_EMTInitializationRefusal(
                    :nonconvergent_operating_point,
                    :fixed_source_operating_point,
                    :active_reactive_power,
                    sprint(showerror, error),
                    (exception_type=string(typeof(error)),),
                ))
            end
            working_parsed = applied.deck
            fixed_source_load_flow = applied.load_flow
        end
        _validate_emt_initialization_source_frequency(working_parsed, request)
        transformer_initial_sample = nothing
        if saturated_transformer_intake !== nothing
            request.formulation isa PhysicalFrequencyFormulation || throw(
                _EMTInitializationRefusal(
                    :unsupported_formulation,
                    :saturated_transformer_magnetic_state,
                    :harmonic_formulation,
                    "saturated-transformer residual-flux initialization currently owns a physical-frequency equilibrium",
                    (requested_formulation=
                        _emt_harmonic_formulation_symbol(request.formulation),),
                ),
            )
            request.operating_point === nothing || throw(
                _EMTInitializationRefusal(
                    :unsupported_operating_point_mapping,
                    :saturated_transformer_magnetic_state,
                    :operating_point,
                    "saturated-transformer operating-point import requires explicit winding current and residual-flux quantities",
                    (operating_point_type=string(typeof(request.operating_point)),),
                ),
            )
            request.time_origin_s == 0.0 || throw(
                _EMTInitializationRefusal(
                    :unsupported_time_origin,
                    :saturated_transformer_magnetic_state,
                    :time_origin_s,
                    "declared saturated-transformer residual flux is referenced to the deck time origin",
                    (time_origin_s=request.time_origin_s,),
                ),
            )
            all(frequency -> isapprox(
                frequency,
                request.frequency_hz;
                atol=1.0e-12,
                rtol=1.0e-12,
            ), request.frequency_grid_hz) || throw(
                _EMTInitializationRefusal(
                    :unsupported_frequency_scan,
                    :saturated_transformer_magnetic_state,
                    :frequency_grid_hz,
                    "one residual-flux state cannot be reused across a frequency grid",
                    (frequency_grid_hz=copy(request.frequency_grid_hz),),
                ),
            )
            transformer_initial_sample = deck_steady_state_voltage_phasors(
                working_parsed;
                saturated_transformer_intake,
            )
            points = EMTInitializationFrequencyPoint[
                _emt_machine_frequency_point(
                    transformer_initial_sample,
                    request;
                    frequency_assignment=:initial_operating_point,
                ),
            ]
        else
            points = _emt_initialization_scan(working_parsed, request)
        end
        primary = _emt_initialization_primary_point(points, request.frequency_hz)
        primary.passed || throw(_emt_initialization_classification_failure(
            primary.topology,
        ))
        mappings = _emt_operating_point_mappings(
            working_parsed,
            request,
            primary,
        )
        sample = transformer_initial_sample === nothing ?
            _emt_initial_voltage_sample(working_parsed, request, primary) :
            transformer_initial_sample
        resolved_timestep_s = timestep_s === nothing ?
            request.formulation isa TimestepMatchedFormulation ?
                request.formulation.timestep_s :
                DeckParser.deck_fixed_time_horizon_options(working_parsed).dt_s :
            Float64(timestep_s)
        isfinite(resolved_timestep_s) && resolved_timestep_s > 0.0 || throw(
            _EMTInitializationRefusal(
                :missing_probe_timestep,
                :study_horizon,
                :timestep_s,
                "initialized EMT state requires a finite positive probe timestep",
                (timestep_s=resolved_timestep_s,),
            ),
        )
        horizon = Float64(t_end_s)
        isfinite(horizon) && horizon >= resolved_timestep_s || throw(
            _EMTInitializationRefusal(
                :missing_probe_horizon,
                :study_horizon,
                :t_end_s,
                "initialized EMT horizon must include at least one timestep",
                (timestep_s=resolved_timestep_s, t_end_s=horizon),
            ),
        )
        prepared = try
            prepare_emt_study(
                working_parsed;
                dt_s=resolved_timestep_s,
                t_end_s=horizon,
                initial_voltage_sample=sample,
                saturated_transformer_branch_runtime_enabled=
                    saturated_transformer_intake !== nothing,
                coupled_lumped_sequence_history_enabled=true,
                recorded_step_indices,
                output_schedule,
                source_signal_provider,
            )
        catch error
            message = sprint(showerror, error)
            if saturated_transformer_intake !== nothing &&
               error isa ArgumentError && occursin("characteristic", message)
                throw(_EMTInitializationRefusal(
                    :invalid_characteristic_state,
                    :saturated_transformer_magnetic_state,
                    :declared_current_flux,
                    message,
                    (exception_type=string(typeof(error)),),
                ))
            end
            rethrow()
        end
        _shift_emt_initialization_time_origin!(prepared, request.time_origin_s)
        transient_metric = _emt_initialization_probe_metric(
            prepared,
            primary,
            request,
        )
        push!(transient_metrics, transient_metric)
        transient_metric.passed || throw(_EMTInitializationRefusal(
            :excessive_artificial_transient,
            :no_artificial_transient,
            transient_metric.quantity,
            "initialized state exceeded the no-artificial-transient threshold",
            (
                normalized_rms=transient_metric.normalized_rms,
                maximum_scaled_discontinuity=
                    transient_metric.maximum_scaled_discontinuity,
                threshold=transient_metric.threshold,
            ),
        ))
        residuals = _emt_initialization_residuals(
            primary,
            mappings,
            request,
            fixed_source_load_flow,
        )
        append!(
            residuals,
            _emt_saturated_transformer_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_pseudo_nonlinear_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_piecewise_nonlinear_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        append!(
            residuals,
            _emt_hysteretic_initialization_residuals(
                prepared,
                primary,
                request,
            ),
        )
        all(residual -> residual.passed, residuals) || throw(
            _EMTInitializationRefusal(
                :residual_failure,
                :initialization_residual,
                :scaled_residual,
                "one or more initialization residuals exceeded their quantity-specific limits",
                (failed_count=count(residual -> !residual.passed, residuals),),
            ),
        )
        state_inventory = _emt_initialization_state_inventory(prepared)
        if fixed_source_load_flow !== nothing
            _append_emt_initialization_state!(
                state_inventory,
                :fixed_source_operating_point,
                :algebraic,
                length(fixed_source_load_flow.constraint_kinds),
                :coupled_active_reactive_power_solution,
            )
            sort!(
                state_inventory;
                by=record -> (
                    String(record.state_family),
                    String(record.owner),
                ),
            )
        end
        initialized_state_owners =
            _emt_initialized_state_owners(state_inventory)
        signature = _emt_initialization_state_signature(
            prepared,
            request,
            mappings,
        )
        report = EMTInitializationReport(
            :accepted,
            _emt_harmonic_formulation_symbol(request.formulation),
            request.frequency_hz,
            request.time_origin_s,
            primary.topology,
            points,
            residuals,
            mappings,
            state_inventory,
            initialized_state_owners,
            Symbol[],
            transient_metrics,
            String[],
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
                :initialization_orchestration,
                :state,
                sprint(showerror, error),
                (exception_type=string(typeof(error)),),
            )
        return _emt_initialization_failure_result(
            request,
            refusal;
            points,
            mappings,
            residuals,
            initialized_state_owners,
            state_inventory,
            unsupported_state_owners,
            transient_metrics,
        )
    end
end

function initialize_emt_study(
    lines,
    request::EMTInitializationRequest;
    source::AbstractString="deck",
    kwargs...,
)
    parsed = DeckParser.parse_deck_lines(lines; source)
    return initialize_emt_study(parsed, request; kwargs...)
end

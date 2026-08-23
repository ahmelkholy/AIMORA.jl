function _apply_steady_state_initial_sample!(context::EMTStepContext, sample)
    return _apply_steady_state_initial_sample!(context, sample, nothing)
end

function _apply_steady_state_initial_sample!(
    context::EMTStepContext,
    sample,
    output_values::Union{Nothing,AbstractVector{<:Real}},
)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state initial sample must cover context nodes"))
    voltage = Float64[values[node] for node in 1:context.system.node_count]
    _update_deck_power_energy_state!(context, voltage)
    output_sample = @view context.output_step_values[:, 1]
    if !isempty(output_sample)
        if output_values === nothing
            _record_context_outputs!(context.output_step_values, 1, context, voltage)
        else
            length(output_values) == length(context.output_channel_names) ||
                throw(ArgumentError("steady-state initial output sample length mismatch"))
            for index in eachindex(output_values)
                output_sample[index] = Float64(output_values[index])
            end
        end
    end
    _update_context_extrema!(context, voltage, output_sample)
    column = _recorded_trace_column!(context, 0)
    if column > 0
        context.time_s[column] = 0.0
        for node in 1:context.system.node_count
            context.voltage_pu[node, column] = voltage[node]
        end
        context.output_pu[:, column] .= output_sample
    end
    return context
end

function _steady_state_initial_output_values(context::EMTStepContext)
    isempty(context.output_channel_names) && return Float64[]
    values = zeros(Float64, length(context.output_channel_names), 1)
    _record_context_outputs!(values, 1, context, context.system.v)
    return vec(values)
end

function _steady_state_node_voltage_phasor(sample, node::Int)
    node == 0 && return complex(0.0, 0.0)
    hasproperty(sample, :node_voltage_phasors) ||
        throw(ArgumentError("steady-state sample must include node voltage phasors"))
    1 <= node <= length(sample.node_voltage_phasors) ||
        throw(ArgumentError("steady-state phasor sample must cover branch nodes"))
    return sample.node_voltage_phasors[node]
end

function _steady_state_branch_frequency_hz(
    sample,
    from_node::Int,
    to_node::Int,
    default_frequency_hz::Float64,
)
    hasproperty(sample, :node_steady_state_frequencies_hz) ||
        return default_frequency_hz
    node_frequencies_hz = sample.node_steady_state_frequencies_hz
    selected_frequency_hz = nothing
    for node in (from_node, to_node)
        node == 0 && continue
        # Transformer and machine assembly can append internal nodes after the
        # parsed deck partition is formed. Their mixed-frequency paths are
        # rejected before assembly, so those nodes inherit the deck default.
        1 <= node <= length(node_frequencies_hz) || continue
        node_frequency_hz = Float64(node_frequencies_hz[node])
        node_frequency_hz > 0.0 || continue
        if selected_frequency_hz === nothing
            selected_frequency_hz = node_frequency_hz
        elseif selected_frequency_hz != node_frequency_hz
            throw(ArgumentError("steady-state branch endpoints have different frequencies"))
        end
    end
    return selected_frequency_hz === nothing ?
        default_frequency_hz : selected_frequency_hz
end

function _steady_state_reactive_angular_frequency(sample, frequency_hz::Float64)
    physical_angular_frequency = 2.0 * pi * frequency_hz
    formulation = hasproperty(sample, :harmonic_formulation) ?
        getproperty(sample, :harmonic_formulation) : :physical_frequency
    if formulation === :physical_frequency
        return physical_angular_frequency
    elseif formulation === :timestep_matched
        hasproperty(sample, :timestep_s) || throw(ArgumentError(
            "timestep-matched steady-state sample must declare timestep_s",
        ))
        return _emt_reactive_angular_frequency(
            TimestepMatchedFormulation(Float64(sample.timestep_s)),
            physical_angular_frequency,
        )
    end
    throw(ArgumentError("unsupported steady-state harmonic formulation $formulation"))
end

function _seed_steady_state_series_rl_branch!(
    branch::SeriesRLBranch,
    sample,
    frequency_hz::Float64,
)
    isfinite(frequency_hz) && frequency_hz >= 0.0 ||
        throw(ArgumentError("series R-L steady-state frequency must be finite and nonnegative"))
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    all(
        isfinite,
        (
            real(branch_voltage_phasor),
            imag(branch_voltage_phasor),
        ),
    ) || throw(ArgumentError("series R-L steady-state voltage must be finite"))
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    impedance =
        branch.l <= 0.0 ?
        complex(branch.r, 0.0) :
        complex(branch.r, reactive_angular_frequency * branch.l)
    current_phasor =
        abs(impedance) == 0.0 ? complex(0.0, 0.0) : branch_voltage_phasor / impedance
    all(
        isfinite,
        (
            real(current_phasor),
            imag(current_phasor),
        ),
    ) || throw(ArgumentError("series R-L steady-state current must be finite"))
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    return branch
end

function _seed_steady_state_series_rlc_branch!(
    branch::SeriesRLCBranch,
    sample,
    frequency_hz::Float64,
)
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    omega = _steady_state_reactive_angular_frequency(sample, frequency_hz)
    if omega <= 0.0
        branch.v_prev = real(branch_voltage_phasor)
        branch.i_prev = 0.0
        branch.i_last = 0.0
        branch.inductor_voltage_prev = 0.0
        branch.capacitor_voltage_prev = branch.v_prev
        return branch
    end
    impedance = complex(branch.r, omega * branch.l - inv(omega * branch.c))
    current_phasor =
        abs(impedance) == 0.0 ? complex(0.0, 0.0) : branch_voltage_phasor / impedance
    inductor_voltage_phasor = im * omega * branch.l * current_phasor
    capacitor_voltage_phasor = current_phasor / (im * omega * branch.c)
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    branch.inductor_voltage_prev = real(inductor_voltage_phasor)
    branch.capacitor_voltage_prev = real(capacitor_voltage_phasor)
    return branch
end

function _seed_steady_state_capacitor_branch!(
    branch::CapacitorBranch,
    sample,
    frequency_hz::Float64,
)
    branch_voltage_phasor =
        _steady_state_node_voltage_phasor(sample, branch.a) -
        _steady_state_node_voltage_phasor(sample, branch.b)
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    current_phasor = im * reactive_angular_frequency * branch.c * branch_voltage_phasor
    branch.v_prev = real(branch_voltage_phasor)
    branch.i_prev = real(current_phasor)
    branch.i_last = branch.i_prev
    return branch
end

function _coupled_inductive_steady_state_admittance(
    branch::CoupledInductiveBranch,
    angular_frequency::Float64=branch.angular_frequency,
)
    isfinite(angular_frequency) && angular_frequency > 0.0 || throw(ArgumentError(
        "coupled-inductive steady-state angular frequency must be finite and positive",
    ))
    frequency_scale = branch.angular_frequency / angular_frequency
    scaled_susceptance = frequency_scale .* branch.susceptance
    admittance = im .* scaled_susceptance
    if branch.series_resistance > 0.0
        reference_susceptance = scaled_susceptance[
            branch.resistance_reference_port,
            branch.resistance_reference_port,
        ]
        reference_reactance = -inv(reference_susceptance)
        admittance .*= (im * reference_reactance) /
                       (branch.series_resistance + im * reference_reactance)
    end
    return admittance
end

function _seed_steady_state_coupled_inductive_branch!(
    branch::CoupledInductiveBranch,
    sample,
)
    port_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, branch.a[index]) -
        _steady_state_node_voltage_phasor(sample, branch.b[index])
        for index in eachindex(branch.a)
    ]
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    current_phasors = _coupled_inductive_steady_state_admittance(
        branch,
        reactive_angular_frequency,
    ) *
                      port_voltage_phasors
    branch.previous_voltage .= real.(port_voltage_phasors)
    branch.previous_current .= real.(current_phasors)
    branch.last_current .= branch.previous_current
    return branch
end

function _seed_steady_state_coupled_series_rl_branch!(
    branch::CoupledSeriesRLBranch,
    sample,
    frequency_hz::Float64,
)
    port_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, branch.a[index]) -
        _steady_state_node_voltage_phasor(sample, branch.b[index])
        for index in eachindex(branch.a)
    ]
    reactive_angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    impedance = complex.(
        branch.resistance_matrix,
        reactive_angular_frequency .* branch.inductance_matrix,
    )
    current_phasors = impedance \ port_voltage_phasors
    branch.previous_voltage .= real.(port_voltage_phasors)
    branch.previous_current .= real.(current_phasors)
    branch.last_current .= branch.previous_current
    return branch
end

function _sequence_modal_phasors(phase_phasors::AbstractVector{<:Complex})
    nph = length(phase_phasors)
    nph > 0 || throw(ArgumentError("phase_phasors must not be empty"))
    modal = Vector{ComplexF64}(undef, nph)
    sum_phasor = complex(0.0, 0.0)
    for phase in 1:nph
        sum_phasor += phase_phasors[phase]
    end
    modal[1] = sum_phasor / nph
    first_phase = phase_phasors[1]
    for phase in 2:nph
        modal[phase] = (first_phase - phase_phasors[phase]) / nph
    end
    return modal
end

function _frequency_domain_sequence_admittance(record, angular_frequency::Float64)
    resistance = Float64(record.r)
    inductive_reactance = Float64(record.l) * angular_frequency
    damping_resistance = Float64(record.rl)
    capacitive_reactance =
        record.c > 0.0 ? inv(Float64(record.c) * angular_frequency) : 0.0
    real_impedance = resistance
    imaginary_impedance = inductive_reactance
    if damping_resistance != 0.0 && inductive_reactance != 0.0
        denominator = inv(damping_resistance^2 + inductive_reactance^2)
        real_impedance += damping_resistance * inductive_reactance^2 * denominator
        imaginary_impedance = inductive_reactance * damping_resistance^2 * denominator
    end
    imaginary_impedance -= capacitive_reactance
    denominator = inv(real_impedance^2 + imaginary_impedance^2)
    admittance = complex(
        real_impedance * denominator,
        -imaginary_impedance * denominator,
    )
    return admittance, capacitive_reactance
end

function _seed_lumped_sequence_frequency_histories!(
    element::BreqivHistoryInjection,
    sample,
    dt_s::Float64,
    frequency_hz::Float64,
)
    angular_frequency =
        _steady_state_reactive_angular_frequency(sample, frequency_hz)
    branch_voltage_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, element.a[phase]) -
        _steady_state_node_voltage_phasor(sample, element.b[phase])
        for phase in eachindex(element.a)
    ]
    seed_breqiv_frequency_histories!(
        element,
        branch_voltage_phasors,
        angular_frequency,
    )
    initialize_breqiv_history_injection!(element, dt_s)
    return element
end

function _seed_semlyen_line_steady_state!(
    element::SemlyenFrequencyDependentLine,
    sample,
    default_frequency_hz::Float64,
)
    frequencies = Float64[
        _steady_state_branch_frequency_hz(
            sample,
            element.from_nodes[phase],
            element.to_nodes[phase],
            default_frequency_hz,
        )
        for phase in eachindex(element.from_nodes)
        if element.from_nodes[phase] != 0 || element.to_nodes[phase] != 0
    ]
    frequency_hz = isempty(frequencies) ? default_frequency_hz : first(frequencies)
    all(value -> isapprox(value, frequency_hz; atol = 1.0e-9, rtol = 1.0e-9), frequencies) ||
        throw(ArgumentError("Semlyen line phases have different steady-state frequencies"))
    from_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.from_nodes
    ]
    to_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.to_nodes
    ]
    initialize_semlyen_line_steady_state!(
        element,
        from_phasors,
        to_phasors,
        frequency_hz,
    )
    return element
end

function _seed_complex_modal_line_steady_state!(
    element::ComplexModalBergeronLine,
    sample,
    default_frequency_hz::Float64,
)
    frequencies = Float64[
        _steady_state_branch_frequency_hz(
            sample,
            element.from_nodes[phase],
            element.to_nodes[phase],
            default_frequency_hz,
        )
        for phase in eachindex(element.from_nodes)
        if element.from_nodes[phase] != 0 || element.to_nodes[phase] != 0
    ]
    frequency_hz = isempty(frequencies) ? default_frequency_hz : first(frequencies)
    all(value -> isapprox(value, frequency_hz; atol = 1.0e-9, rtol = 1.0e-9), frequencies) ||
        throw(ArgumentError("complex modal line phases have different steady-state frequencies"))
    from_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.from_nodes
    ]
    to_phasors = ComplexF64[
        _steady_state_node_voltage_phasor(sample, node) for node in element.to_nodes
    ]
    initialize_complex_modal_bergeron_steady_state!(
        element,
        from_phasors,
        to_phasors,
        frequency_hz,
    )
    return element
end

function _steady_state_current_injections(context::EMTStepContext, sample)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state current injection seed must cover context nodes"))
    return nodal_current_injections_for_voltage!(
        context.system,
        0.0,
        context.dt_s,
        Float64.(values[1:context.system.node_count]),
    )
end

function _seed_steady_state_network_state!(context::EMTStepContext, sample)
    values = sample.node_voltage_values
    length(values) >= context.system.node_count ||
        throw(ArgumentError("steady-state network seed must cover context nodes"))
    for node in 1:context.system.node_count
        context.system.v[node] = Float64(values[node])
    end
    frequency_hz = Float64(sample.steady_state_frequency_hz)
    for element in context.system.elements
        if element isa ComplexModalBergeronLine &&
           hasproperty(sample, :node_voltage_phasors)
            _seed_complex_modal_line_steady_state!(element, sample, frequency_hz)
            continue
        elseif element isa SemlyenFrequencyDependentLine &&
           hasproperty(sample, :node_voltage_phasors)
            _seed_semlyen_line_steady_state!(element, sample, frequency_hz)
            continue
        elseif element isa BreqivHistoryInjection
            if hasproperty(sample, :node_voltage_phasors)
                element_frequency_hz = _steady_state_branch_frequency_hz(
                    sample,
                    first(element.a),
                    first(element.b),
                    frequency_hz,
                )
                _seed_lumped_sequence_frequency_histories!(
                    element,
                    sample,
                    context.dt_s,
                    element_frequency_hz,
                )
                # The recorded steady-state sample replaces the t=0 solve. Advance
                # the BREQIV history once so its Norton source is staged for the
                # first dynamic solve at t=dt, like the other companion branches.
                advance_breqiv_history_current!(
                    element,
                    context.system.v,
                    context.dt_s;
                    consumed_for_step = false,
                )
            else
                for phase in eachindex(element.initial_phase_voltage)
                    from_node = element.a[phase]
                    to_node = element.b[phase]
                    from_voltage = from_node == 0 ? 0.0 : context.system.v[from_node]
                    to_voltage = to_node == 0 ? 0.0 : context.system.v[to_node]
                    element.initial_phase_voltage[phase] = from_voltage - to_voltage
                end
            end
            continue
        elseif element isa SeriesRLBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_series_rl_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa SeriesRLCBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_series_rlc_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa CapacitorBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                element.a,
                element.b,
                frequency_hz,
            )
            _seed_steady_state_capacitor_branch!(element, sample, element_frequency_hz)
            continue
        elseif element isa CoupledInductiveBranch && hasproperty(sample, :node_voltage_phasors)
            _seed_steady_state_coupled_inductive_branch!(element, sample)
            continue
        elseif element isa CoupledSeriesRLBranch && hasproperty(sample, :node_voltage_phasors)
            element_frequency_hz = _steady_state_branch_frequency_hz(
                sample,
                first(element.a),
                first(element.b),
                frequency_hz,
            )
            _seed_steady_state_coupled_series_rl_branch!(
                element,
                sample,
                element_frequency_hz,
            )
            continue
        end
        if !hasproperty(sample, :node_voltage_phasors)
            if element isa SeriesRLBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                continue
            elseif element isa SeriesRLCBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                element.inductor_voltage_prev = 0.0
                element.capacitor_voltage_prev = element.v_prev
                continue
            elseif element isa CapacitorBranch
                element.v_prev =
                    _deck_node_voltage(context.system.v, element.a) -
                    _deck_node_voltage(context.system.v, element.b)
                element.i_prev = 0.0
                element.i_last = 0.0
                continue
            elseif element isa CoupledInductiveBranch
                for port in eachindex(element.a)
                    element.previous_voltage[port] =
                        _deck_node_voltage(context.system.v, element.a[port]) -
                        _deck_node_voltage(context.system.v, element.b[port])
                end
                fill!(element.previous_current, 0.0)
                fill!(element.last_current, 0.0)
                continue
            elseif element isa CoupledSeriesRLBranch
                for port in eachindex(element.a)
                    element.previous_voltage[port] =
                        _deck_node_voltage(context.system.v, element.a[port]) -
                        _deck_node_voltage(context.system.v, element.b[port])
                end
                fill!(element.previous_current, 0.0)
                fill!(element.last_current, 0.0)
                continue
            end
        end
        update!(element, context.system.v, context.dt_s)
    end
    return context
end

function _seed_direct_machine_power_leakage_currents!(
    context::EMTStepContext,
    parsed::DeckParser.DeckParseResult,
    machine_indices::AbstractVector{<:Integer},
    power_terminal_currents::AbstractVector{<:Real},
)
    indices = Int.(machine_indices)
    currents = Float64.(power_terminal_currents)
    length(indices) == length(currents) ||
        throw(ArgumentError("direct-machine leakage-current seeds must align with machine indices"))
    for (machine_index, current) in zip(indices, currents)
        rows = [
            row
            for row in DeckParser.deck_universal_machine_generated_branch_rows(parsed)
            if row.machine_index == machine_index &&
               row.from_node_value != row.to_node_value
        ]
        length(rows) == 1 ||
            throw(ArgumentError("automatic direct machine $machine_index requires one generated power-leakage branch"))
        row = only(rows)
        matches = SeriesRLBranch[
            element
            for element in context.system.elements
            if element isa SeriesRLBranch &&
               element.a == row.from_node_value &&
               element.b == row.to_node_value
        ]
        length(matches) == 1 ||
            throw(ArgumentError("automatic direct machine $machine_index power-leakage runtime branch is missing or ambiguous"))
        branch = only(matches)
        branch.i_prev = current
        branch.i_last = current
        branch.v_prev = branch.r * current
    end
    return context
end

function _deck_synchronous_machine_delta_connected(
    parsed::DeckParser.DeckParseResult,
    machine_index::Int,
)
    return any(
        row -> row.machine_index == machine_index &&
               row.parameter_kind == :delta_connection,
        DeckParser.deck_synchronous_machine_model_parameter_rows(parsed),
    )
end

function _synchronous_machine_winding_voltages(
    terminal_node_voltages::AbstractVector{<:Real},
    delta_connected::Bool,
)
    length(terminal_node_voltages) == 3 || throw(ArgumentError(
        "synchronous-machine terminal voltage must contain three phases",
    ))
    values = Float64.(terminal_node_voltages)
    delta_connected || return values
    return Float64[
        values[1] - values[2],
        values[2] - values[3],
        values[3] - values[1],
    ]
end

function _synchronous_machine_terminal_currents(
    winding_currents::AbstractVector{<:Real},
    delta_connected::Bool,
)
    length(winding_currents) == 3 || throw(ArgumentError(
        "synchronous-machine winding current must contain three phases",
    ))
    values = Float64.(winding_currents)
    delta_connected || return values
    return Float64[
        values[1] - values[3],
        values[2] - values[1],
        values[3] - values[2],
    ]
end

function _deck_synchronous_machine_terminal_admittance(
    parsed::DeckParser.DeckParseResult,
    state::SynchronousMachineDynamicState,
    machine_index::Int=1,
    delta_reference_phase_admittance::Union{Nothing,Real}=nothing,
)
    direct = Float64(state.electrical_coefficients[27])
    mutual = Float64(state.electrical_coefficients[28])
    delta_connected = _deck_synchronous_machine_delta_connected(parsed, machine_index)
    if delta_connected
        phase_admittance = direct - mutual
        if delta_reference_phase_admittance !== nothing
            # PAST retains the initial delta companion in the network base;
            # UPDATE maps only later winding-admittance changes through the delta.
            reference = Float64(delta_reference_phase_admittance)
            phase_admittance = reference + (phase_admittance - reference) / 3.0
        end
        terminal_admittance = fill(-phase_admittance, 3, 3)
        for phase in 1:3
            terminal_admittance[phase, phase] = 2.0 * phase_admittance
        end
        return terminal_admittance
    end
    terminal_admittance = fill(mutual, 3, 3)
    for phase in 1:3
        terminal_admittance[phase, phase] = direct
    end
    return terminal_admittance
end

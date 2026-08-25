export ConverterPulseWidthModulationState,
       ConverterThreeLevelModulationState,
       converter_two_level_duties,
       converter_triangular_carrier,
       converter_pwm_gate_state,
       converter_three_level_gate_state,
       converter_flying_capacitor_gate_state,
       ConverterCascadedHBridgeModulationState,
       converter_cascaded_h_bridge_gate_state,
       selective_harmonic_elimination_residual,
       nearest_level_cell_state,
       dual_active_bridge_gate_state,
       ConverterMatrixSpaceVectorState,
       converter_matrix_space_vector_state,
       matrix_safe_commutation_sequence

struct ConverterPulseWidthModulationState
    time_s::Float64
    kind::ConverterModulationKind
    duties::Vector{Float64}
    carrier_values::Vector{Float64}
    requested_valve_state::BitVector
    deterministic_signature_sha256::String
end

function converter_flying_capacitor_gate_state(
    normalized_reference,
    flying_capacitor_voltage_error_v,
    phase_current_a,
    time_s::Real,
    carrier_frequency_hz::Real;
    carrier_phases_rad=nothing,
    balance_voltage_tolerance_v::Real=1.0e-3,
)
    references = Float64[value for value in normalized_reference]
    voltage_error = Float64[value for value in flying_capacitor_voltage_error_v]
    current = Float64[value for value in phase_current_a]
    count = length(references)
    count > 0 && length(voltage_error) == count && length(current) == count ||
        throw(DimensionMismatch(
            "flying-capacitor modulation requires one reference, voltage error, and current per leg",
        ))
    all(isfinite, (references..., voltage_error..., current...)) &&
        all(value -> -1.0 <= value <= 1.0, references) || throw(ArgumentError(
            "flying-capacitor modulation inputs must be finite with references in [-1, 1]",
        ))
    phases = carrier_phases_rad === nothing ? zeros(count) :
        Float64.(carrier_phases_rad)
    length(phases) == count && all(isfinite, phases) || throw(ArgumentError(
        "flying-capacitor modulation requires one finite carrier phase per leg",
    ))
    time = Float64(time_s)
    frequency = Float64(carrier_frequency_hz)
    balance_tolerance = Float64(balance_voltage_tolerance_v)
    isfinite(balance_tolerance) && balance_tolerance >= 0.0 || throw(ArgumentError(
        "flying-capacitor balance tolerance must be finite and nonnegative",
    ))
    carriers = [
        converter_triangular_carrier(time, frequency; phase_rad=phase)
        for phase in phases
    ]
    levels = Vector{Int8}(undef, count)
    gates = falses(4 * count)
    for index in eachindex(references)
        reference = references[index]
        level = reference >= carriers[index] ? Int8(1) :
            -reference >= carriers[index] ? Int8(-1) : Int8(0)
        levels[index] = level
        offset = 4 * (index - 1)
        if level == 1
            gates[(offset + 1):(offset + 2)] .= true
        elseif level == -1
            gates[(offset + 3):(offset + 4)] .= true
        else
            charge_state = abs(voltage_error[index]) > balance_tolerance ?
                voltage_error[index] * current[index] < 0.0 :
                iseven(floor(Int, 2.0 * frequency * time + index - 1))
            if charge_state
                gates[offset + 1] = true
                gates[offset + 3] = true
            else
                gates[offset + 2] = true
                gates[offset + 4] = true
            end
        end
    end
    io = IOBuffer()
    println(io, "aimora_flying_capacitor_modulation_state_v1")
    println(io, bitstring(time))
    foreach(value -> println(io, bitstring(value)), references)
    foreach(value -> println(io, bitstring(value)), voltage_error)
    foreach(value -> println(io, bitstring(value)), current)
    println(io, join(levels, ','))
    println(io, join(Int.(gates)))
    return ConverterThreeLevelModulationState(
        time,
        references,
        carriers,
        levels,
        gates,
        bytes2hex(sha256(take!(io))),
    )
end

struct ConverterThreeLevelModulationState
    time_s::Float64
    normalized_reference::Vector{Float64}
    carrier_values::Vector{Float64}
    requested_level::Vector{Int8}
    requested_valve_state::BitVector
    deterministic_signature_sha256::String
end

function converter_three_level_gate_state(
    normalized_reference,
    time_s::Real,
    carrier_frequency_hz::Real;
    carrier_phases_rad=nothing,
    leg_topology::Symbol=:neutral_point_clamped,
)
    references = Float64[value for value in normalized_reference]
    !isempty(references) && all(value -> isfinite(value) && -1.0 <= value <= 1.0, references) ||
        throw(ArgumentError("three-level references must be finite values in [-1, 1]"))
    phases = carrier_phases_rad === nothing ? zeros(length(references)) :
        Float64.(carrier_phases_rad)
    length(phases) == length(references) && all(isfinite, phases) ||
        throw(ArgumentError("three-level PWM requires one finite carrier phase per reference"))
    leg_topology in (:neutral_point_clamped, :t_type) || throw(ArgumentError(
        "three-level PWM requires neutral-point-clamped or T-type leg mapping",
    ))
    time = Float64(time_s)
    carriers = [
        converter_triangular_carrier(time, carrier_frequency_hz; phase_rad=phase)
        for phase in phases
    ]
    levels = Vector{Int8}(undef, length(references))
    gates = falses(4 * length(references))
    for index in eachindex(references)
        reference = references[index]
        level = reference >= carriers[index] ? Int8(1) :
            -reference >= carriers[index] ? Int8(-1) : Int8(0)
        levels[index] = level
        offset = 4 * (index - 1)
        if level == 1
            gates[offset + 1] = true
            leg_topology === :neutral_point_clamped && (gates[offset + 2] = true)
        elseif level == 0
            gates[(offset + 2):(offset + 3)] .= true
        else
            leg_topology === :neutral_point_clamped && (gates[offset + 3] = true)
            gates[offset + 4] = true
        end
    end
    io = IOBuffer()
    println(io, "aimora_three_level_modulation_state_v1")
    println(io, leg_topology)
    println(io, bitstring(time))
    foreach(value -> println(io, bitstring(value)), references)
    foreach(value -> println(io, bitstring(value)), carriers)
    println(io, join(levels, ','))
    println(io, join(Int.(gates)))
    return ConverterThreeLevelModulationState(
        time,
        references,
        carriers,
        levels,
        gates,
        bytes2hex(sha256(take!(io))),
    )
end

function converter_two_level_duties(
    phase_voltage_reference_v::NTuple{3,<:Real},
    dc_link_voltage_v::Real,
    kind::ConverterModulationKind;
    minimum_duty::Real=0.0,
    maximum_duty::Real=1.0,
)
    modulation = kind === CarrierSinusoidalPulseWidthModulation ?
        SwitchDetailedVSC.SinusoidalPulseWidthModulation :
        kind === SpaceVectorPulseWidthModulation ?
            SwitchDetailedVSC.ZeroSequenceInjectedPulseWidthModulation :
            throw(ArgumentError("two-level duties require sinusoidal or space-vector PWM"))
    return SwitchDetailedVSC.modulation_duties(
        phase_voltage_reference_v,
        dc_link_voltage_v,
        modulation;
        minimum_duty,
        maximum_duty,
    )
end

function converter_triangular_carrier(
    time_s::Real,
    frequency_hz::Real;
    phase_rad::Real=0.0,
)
    time = Float64(time_s)
    frequency = Float64(frequency_hz)
    phase = Float64(phase_rad)
    all(isfinite, (time, frequency, phase)) || throw(ArgumentError(
        "converter carrier time, frequency, and phase must be finite",
    ))
    time >= 0.0 && frequency > 0.0 || throw(ArgumentError(
        "converter carrier requires nonnegative time and positive frequency",
    ))
    cycle = mod(frequency * time + phase / (2.0 * pi), 1.0)
    return 1.0 - abs(2.0 * cycle - 1.0)
end

function _converter_pwm_signature(time, kind, duties, carrier_values, requested)
    io = IOBuffer()
    println(io, "aimora_converter_pwm_state_v1")
    println(io, bitstring(time))
    println(io, Int(kind))
    foreach(value -> println(io, bitstring(value)), duties)
    foreach(value -> println(io, bitstring(value)), carrier_values)
    println(io, join(Int.(requested)))
    return bytes2hex(sha256(take!(io)))
end

function converter_pwm_gate_state(
    duties,
    time_s::Real,
    carrier_frequency_hz::Real,
    kind::ConverterModulationKind;
    carrier_phases_rad=nothing,
)
    kind in (
        CarrierSinusoidalPulseWidthModulation,
        SpaceVectorPulseWidthModulation,
        PhaseShiftedCarrierPulseWidthModulation,
    ) || throw(ArgumentError("gate-state synthesis requires a carrier PWM method"))
    duty_values = Float64[value for value in duties]
    !isempty(duty_values) && all(value -> isfinite(value) && 0.0 <= value <= 1.0, duty_values) ||
        throw(ArgumentError("converter PWM duties must be finite values in [0, 1]"))
    phases = carrier_phases_rad === nothing ? zeros(length(duty_values)) :
        Float64.(carrier_phases_rad)
    length(phases) == length(duty_values) && all(isfinite, phases) || throw(ArgumentError(
        "converter PWM requires one finite carrier phase per duty",
    ))
    time = Float64(time_s)
    carrier_values = [
        converter_triangular_carrier(time, carrier_frequency_hz; phase_rad=phase)
        for phase in phases
    ]
    requested = BitVector(undef, 2 * length(duty_values))
    for index in eachindex(duty_values)
        upper_on = duty_values[index] >= carrier_values[index]
        requested[2 * index - 1] = upper_on
        requested[2 * index] = !upper_on
    end
    signature = _converter_pwm_signature(
        time,
        kind,
        duty_values,
        carrier_values,
        requested,
    )
    return ConverterPulseWidthModulationState(
        time,
        kind,
        duty_values,
        carrier_values,
        requested,
        signature,
    )
end

function selective_harmonic_elimination_residual(
    switching_angles_rad,
    fundamental_ratio::Real,
    eliminated_harmonic_orders,
)
    angles = Float64.(switching_angles_rad)
    !isempty(angles) && all(isfinite, angles) &&
        all(angle -> 0.0 < angle < pi / 2.0, angles) && issorted(angles) &&
        length(unique(angles)) == length(angles) || throw(ArgumentError(
        "SHE switching angles must be unique, increasing, finite, and inside (0, pi/2)",
    ))
    target = Float64(fundamental_ratio)
    isfinite(target) && 0.0 <= target <= 4.0 * length(angles) / pi || throw(ArgumentError(
        "SHE fundamental ratio is outside the admitted quarter-wave range",
    ))
    orders = Int.(eliminated_harmonic_orders)
    length(orders) == length(angles) - 1 && all(order -> isodd(order) && order > 1, orders) &&
        length(unique(orders)) == length(orders) || throw(ArgumentError(
        "SHE requires one unique odd harmonic order above one per nonfundamental angle",
    ))
    signs = [isodd(index) ? 1.0 : -1.0 for index in eachindex(angles)]
    residual = Vector{Float64}(undef, length(angles))
    residual[1] = 4.0 / pi * sum(signs .* cos.(angles)) - target
    for (index, order) in enumerate(orders)
        residual[index + 1] = sum(signs .* cos.(order .* angles))
    end
    return residual
end

function nearest_level_cell_state(reference_pu::Real, cell_count::Integer)
    reference = Float64(reference_pu)
    count = Int(cell_count)
    isfinite(reference) && -1.0 <= reference <= 1.0 || throw(ArgumentError(
        "nearest-level reference must be finite and inside [-1, 1]",
    ))
    2 <= count <= 8 || throw(ArgumentError(
        "nearest-level cascaded conversion requires two through eight cells",
    ))
    requested_level = clamp(round(Int, count * reference), -count, count)
    cells = zeros(Int8, count)
    if requested_level > 0
        cells[1:requested_level] .= Int8(1)
    elseif requested_level < 0
        cells[1:(-requested_level)] .= Int8(-1)
    end
    return (
        requested_level=requested_level,
        cell_state=cells,
        normalized_voltage=requested_level / count,
        quantization_error=reference - requested_level / count,
    )
end

struct ConverterCascadedHBridgeModulationState
    time_s::Float64
    normalized_reference::Vector{Float64}
    requested_level::Vector{Int8}
    requested_cell_state::Matrix{Int8}
    requested_valve_state::BitVector
    deterministic_signature_sha256::String
end

function _cascaded_h_bridge_selected_cells(
    cell_voltage_error_v,
    phase_current_a,
    requested_polarity,
    selected_count,
    rotation,
)
    cell_count = length(cell_voltage_error_v)
    tie_order = mod.((0:(cell_count - 1)) .- rotation, cell_count)
    discharge = requested_polarity * phase_current_a > 0.0
    order = sortperm(1:cell_count; by=index -> (
        discharge ? -cell_voltage_error_v[index] : cell_voltage_error_v[index],
        tie_order[index],
    ))
    return order[1:selected_count]
end

function converter_cascaded_h_bridge_gate_state(
    normalized_reference,
    cell_voltage_error_v,
    phase_current_a,
    time_s::Real,
    carrier_frequency_hz::Real,
    kind::ConverterModulationKind;
    modulation_index::Real=1.0,
    selective_harmonic_angles_rad=(),
)
    references = Float64[value for value in normalized_reference]
    voltage_error = Matrix{Float64}(cell_voltage_error_v)
    current = Float64[value for value in phase_current_a]
    phase_count, cell_count = size(voltage_error)
    phase_count > 0 && 2 <= cell_count <= 8 && length(references) == phase_count &&
        length(current) == phase_count || throw(DimensionMismatch(
        "cascaded H-bridge modulation requires one reference/current per phase and two through eight cell voltages",
    ))
    all(isfinite, (references..., voltage_error..., current...)) &&
        all(value -> -1.0 <= value <= 1.0, references) || throw(ArgumentError(
        "cascaded H-bridge modulation inputs must be finite with references in [-1, 1]",
    ))
    kind in (
        CarrierSinusoidalPulseWidthModulation,
        PhaseShiftedCarrierPulseWidthModulation,
        SelectiveHarmonicElimination,
        NearestLevelModulation,
    ) || throw(ArgumentError(
        "cascaded H-bridge modulation requires carrier, selective-harmonic, or nearest-level synthesis",
    ))
    time = Float64(time_s)
    frequency = Float64(carrier_frequency_hz)
    modulation = Float64(modulation_index)
    isfinite(time) && time >= 0.0 && isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError(
            "cascaded H-bridge modulation requires nonnegative time and positive carrier frequency",
        ))
    isfinite(modulation) && 0.0 < modulation <= 1.0 &&
        all(reference -> abs(reference) <= modulation + 64.0eps(1.0), references) ||
        throw(ArgumentError(
            "cascaded H-bridge references must lie inside their positive modulation index",
        ))
    angles = Float64.(selective_harmonic_angles_rad)
    if kind === SelectiveHarmonicElimination
        1 <= length(angles) <= cell_count && issorted(angles) &&
            all(angle -> 0.0 < angle < pi / 2.0, angles) || throw(ArgumentError(
            "cascaded H-bridge selective-harmonic angles must be ordered inside (0, pi/2)",
        ))
    elseif !isempty(angles)
        throw(ArgumentError(
            "selective-harmonic angles are valid only for selective-harmonic modulation",
        ))
    end
    cell_state = zeros(Int8, phase_count, cell_count)
    requested_level = zeros(Int8, phase_count)
    carrier_cycle = floor(Int, frequency * time + 64.0eps(max(1.0, frequency * time)))
    for phase in 1:phase_count
        reference = references[phase]
        rotation = mod(carrier_cycle + phase - 1, cell_count)
        if kind in (
            CarrierSinusoidalPulseWidthModulation,
            PhaseShiftedCarrierPulseWidthModulation,
        )
            for cell in 1:cell_count
                carrier_phase = kind === PhaseShiftedCarrierPulseWidthModulation ?
                    2.0 * pi * (cell - 1) / cell_count : 0.0
                carrier = converter_triangular_carrier(
                    time,
                    frequency;
                    phase_rad=carrier_phase,
                )
                cell_state[phase, cell] = abs(reference) >= carrier ?
                    Int8(sign(reference)) : Int8(0)
            end
        else
            selected_count = if kind === NearestLevelModulation
                abs(nearest_level_cell_state(reference, cell_count).requested_level)
            else
                quarter_angle = asin(clamp(abs(reference) / modulation, 0.0, 1.0))
                count(angle -> quarter_angle >= angle, angles)
            end
            polarity = Int8(sign(reference))
            if selected_count > 0 && polarity != 0
                selected = _cascaded_h_bridge_selected_cells(
                    view(voltage_error, phase, :),
                    current[phase],
                    polarity,
                    selected_count,
                    rotation,
                )
                cell_state[phase, selected] .= polarity
            end
        end
        requested_level[phase] = sum(view(cell_state, phase, :))
    end
    gates = falses(4 * phase_count * cell_count)
    for phase in 1:phase_count, cell in 1:cell_count
        first = 4 * ((phase - 1) * cell_count + cell - 1) + 1
        state = cell_state[phase, cell]
        if state == 1
            gates[first] = true
            gates[first + 3] = true
        elseif state == -1
            gates[first + 1] = true
            gates[first + 2] = true
        elseif iseven(carrier_cycle + phase + cell)
            gates[first] = true
            gates[first + 2] = true
        else
            gates[first + 1] = true
            gates[first + 3] = true
        end
    end
    io = IOBuffer()
    println(io, "aimora_cascaded_h_bridge_modulation_state_v1")
    println(io, bitstring(time))
    println(io, Int(kind))
    println(io, bitstring(modulation))
    foreach(value -> println(io, bitstring(value)), references)
    foreach(value -> println(io, bitstring(value)), voltage_error)
    foreach(value -> println(io, bitstring(value)), current)
    foreach(value -> println(io, bitstring(value)), angles)
    println(io, join(requested_level, ','))
    println(io, join(vec(cell_state), ','))
    println(io, join(Int.(gates)))
    return ConverterCascadedHBridgeModulationState(
        time,
        references,
        requested_level,
        cell_state,
        gates,
        bytes2hex(sha256(take!(io))),
    )
end

function _full_bridge_square_state(angle_rad::Float64)
    half_cycles = angle_rad / pi
    nearest_boundary = round(half_cycles)
    boundary_tolerance = 32.0 * eps(max(1.0, abs(half_cycles)))
    half_cycle_index = abs(half_cycles - nearest_boundary) <= boundary_tolerance ?
        Int(nearest_boundary) : floor(Int, half_cycles)
    positive = iseven(half_cycle_index)
    return positive ? BitVector((true, false, false, true)) :
        BitVector((false, true, true, false))
end

function dual_active_bridge_gate_state(
    time_s::Real,
    switching_frequency_hz::Real,
    phase_shift_rad::Real,
    kind::ConverterModulationKind=SinglePhaseShiftModulation;
    primary_inner_phase_shift_rad::Real=0.0,
    secondary_inner_phase_shift_rad::Real=0.0,
)
    time = Float64(time_s)
    frequency = Float64(switching_frequency_hz)
    shift = Float64(phase_shift_rad)
    primary_inner = Float64(primary_inner_phase_shift_rad)
    secondary_inner = Float64(secondary_inner_phase_shift_rad)
    kind in (SinglePhaseShiftModulation, DualPhaseShiftModulation,
        TriplePhaseShiftModulation) || throw(ArgumentError(
            "DAB gate synthesis requires single-, dual-, or triple-phase-shift modulation",
        ))
    all(isfinite, (time, frequency, shift, primary_inner, secondary_inner)) &&
        time >= 0.0 && frequency > 0.0 && abs(shift) <= pi &&
        0.0 <= primary_inner < pi && 0.0 <= secondary_inner < pi ||
        throw(ArgumentError(
            "DAB gate synthesis requires nonnegative time, positive frequency, |outer shift| <= pi, and inner shifts in [0, pi)",
        ))
    kind === SinglePhaseShiftModulation &&
        (!iszero(primary_inner) || !iszero(secondary_inner)) && throw(ArgumentError(
            "single-phase-shift DAB modulation requires zero inner bridge shifts",
        ))
    kind === DualPhaseShiftModulation &&
        (iszero(primary_inner) == iszero(secondary_inner)) && throw(ArgumentError(
            "dual-phase-shift DAB modulation requires exactly one nonzero inner bridge shift",
        ))
    kind === TriplePhaseShiftModulation &&
        (iszero(primary_inner) || iszero(secondary_inner)) && throw(ArgumentError(
            "triple-phase-shift DAB modulation requires two nonzero inner bridge shifts",
    ))
    primary_angle = 2.0 * pi * frequency * time
    function bridge_state(angle, inner_shift)
        leg_a_upper = _full_bridge_square_state(angle)[1]
        leg_b_upper = iszero(inner_shift) ? !leg_a_upper :
            _full_bridge_square_state(angle - (pi - inner_shift))[1]
        return BitVector((leg_a_upper, !leg_a_upper, leg_b_upper, !leg_b_upper))
    end
    primary = bridge_state(primary_angle, primary_inner)
    secondary = bridge_state(primary_angle - shift, secondary_inner)
    return (
        primary_state=primary,
        secondary_state=secondary,
        requested_valve_state=vcat(primary, secondary),
    )
end

function matrix_safe_commutation_sequence(
    outgoing_input::Integer,
    incoming_input::Integer,
    output_current_a::Real,
)
    outgoing = Int(outgoing_input)
    incoming = Int(incoming_input)
    outgoing in 1:3 && incoming in 1:3 && outgoing != incoming || throw(ArgumentError(
        "matrix commutation requires distinct outgoing and incoming inputs in 1:3",
    ))
    current = Float64(output_current_a)
    isfinite(current) && current != 0.0 || throw(ArgumentError(
        "matrix safe commutation requires a finite nonzero measured output current",
    ))
    conducting_direction = current > 0.0 ? 1 : 2
    nonconducting_direction = 3 - conducting_direction
    function state(outgoing_pair, incoming_pair)
        connection = falses(6)
        connection[2 * outgoing - 2 + 1] = outgoing_pair[1]
        connection[2 * outgoing - 2 + 2] = outgoing_pair[2]
        connection[2 * incoming - 2 + 1] = incoming_pair[1]
        connection[2 * incoming - 2 + 2] = incoming_pair[2]
        return connection
    end
    outgoing_conducting = ntuple(index -> index == conducting_direction, 2)
    incoming_conducting = ntuple(index -> index == conducting_direction, 2)
    all_on = (true, true)
    all_off = (false, false)
    stages = (
        state(all_on, all_off),
        state(outgoing_conducting, all_off),
        state(outgoing_conducting, incoming_conducting),
        state(all_off, incoming_conducting),
        state(all_off, all_on),
    )
    return (
        current_direction=current > 0.0 ? :positive : :negative,
        conducting_direction=conducting_direction,
        nonconducting_direction=nonconducting_direction,
        stages=stages,
    )
end

struct ConverterMatrixSpaceVectorState
    time_s::Float64
    input_voltage_v::NTuple{3,Float64}
    output_reference_voltage_v::NTuple{3,Float64}
    requested_input_for_output::NTuple{3,Int}
    connection::Matrix{Bool}
    requested_valve_state::BitVector
    switching_vector_rank::Int
    carrier_subslot::Int
    deterministic_signature_sha256::String
end

function converter_matrix_space_vector_state(
    input_voltage_v,
    time_s::Real,
    output_frequency_hz::Real,
    carrier_frequency_hz::Real,
    modulation_index::Real;
    output_phase_rad::Real=0.0,
)
    length(input_voltage_v) == 3 || throw(DimensionMismatch(
        "matrix space-vector synthesis requires three input voltages",
    ))
    input_voltage = ntuple(index -> Float64(input_voltage_v[index]), 3)
    time, output_frequency, carrier_frequency, modulation, phase = Float64.((
        time_s,
        output_frequency_hz,
        carrier_frequency_hz,
        modulation_index,
        output_phase_rad,
    ))
    all(isfinite, (input_voltage..., time, output_frequency, carrier_frequency,
        modulation, phase)) && time >= 0.0 && output_frequency > 0.0 &&
        carrier_frequency > 0.0 && 0.0 < modulation <= sqrt(3.0) / 2.0 ||
        throw(ArgumentError(
            "matrix space-vector inputs require finite voltages, nonnegative time, positive frequencies, and modulation in (0, sqrt(3)/2]",
        ))
    input_peak = sqrt(2.0 / 3.0 * sum(abs2, input_voltage))
    input_peak > sqrt(eps(Float64)) || throw(ArgumentError(
        "matrix space-vector synthesis requires a nonzero three-phase input vector",
    ))
    angle = 2.0 * pi * output_frequency * time + phase
    reference = ntuple(output_phase ->
        modulation * input_peak * sin(angle - (output_phase - 1) * 2.0 * pi / 3.0),
        3,
    )
    candidates = NamedTuple[]
    for input_for_output_1 in 1:3,
        input_for_output_2 in 1:3,
        input_for_output_3 in 1:3
        selected = (input_for_output_1, input_for_output_2, input_for_output_3)
        raw = ntuple(output -> input_voltage[selected[output]], 3)
        common_mode = sum(raw) / 3.0
        phase_voltage = ntuple(output -> raw[output] - common_mode, 3)
        error = sum((phase_voltage[index] - reference[index])^2 for index in 1:3) /
            max(input_peak^2, eps(Float64))
        push!(candidates, (; selected, error))
    end
    sort!(candidates; by=candidate -> (candidate.error, candidate.selected))
    carrier_fraction = mod(carrier_frequency * time, 1.0)
    subslot = min(floor(Int, 8.0 * carrier_fraction) + 1, 8)
    rank = (1, 2, 1, 3, 1, 4, 1, 2)[subslot]
    selected = candidates[rank].selected
    connection = falses(3, 3)
    gates = falses(18)
    for output in 1:3
        input = selected[output]
        connection[output, input] = true
        pair = 2 * (3 * (output - 1) + input) - 1
        gates[pair:(pair + 1)] .= true
    end
    io = IOBuffer()
    println(io, "aimora_matrix_space_vector_state_v1")
    println(io, bitstring(time))
    foreach(value -> println(io, bitstring(value)), input_voltage)
    foreach(value -> println(io, bitstring(value)), reference)
    println(io, join(selected, ','))
    println(io, rank)
    println(io, subslot)
    return ConverterMatrixSpaceVectorState(
        time,
        input_voltage,
        reference,
        selected,
        connection,
        gates,
        rank,
        subslot,
        bytes2hex(sha256(take!(io))),
    )
end

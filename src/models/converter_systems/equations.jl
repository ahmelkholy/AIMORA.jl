export ideal_dcdc_conversion_ratio,
       dual_active_bridge_power_w,
       controlled_rectifier_average_voltage_v,
       interleaved_carrier_phases_rad,
       matrix_converter_terminal_map,
       converter_stage_energy_residual_w

function ideal_dcdc_conversion_ratio(family::ConverterSystemFamily, duty::Real)
    d = Float64(duty)
    isfinite(d) && 0.0 < d < 1.0 || throw(ArgumentError("duty must lie strictly between zero and one"))
    family === BuckChopper && return d
    family === BoostChopper && return inv(1.0 - d)
    family === InvertingBuckBoostChopper && return -d / (1.0 - d)
    throw(ArgumentError("the selected converter family has no admitted canonical DC/DC ratio"))
end

function dual_active_bridge_power_w(
    primary_voltage_v::Real,
    secondary_voltage_v::Real,
    transformer_ratio::Real,
    phase_shift_rad::Real,
    angular_frequency_rad_per_s::Real,
    leakage_inductance_h::Real,
)
    values = Float64.((
        primary_voltage_v,
        secondary_voltage_v,
        transformer_ratio,
        phase_shift_rad,
        angular_frequency_rad_per_s,
        leakage_inductance_h,
    ))
    all(isfinite, values) || throw(ArgumentError("dual-active-bridge inputs must be finite"))
    v1, v2, ratio, phase_shift, angular_frequency, leakage = values
    v1 > 0.0 && v2 > 0.0 && ratio > 0.0 && angular_frequency > 0.0 && leakage > 0.0 ||
        throw(ArgumentError("dual-active-bridge voltages, ratio, frequency, and leakage must be positive"))
    abs(phase_shift) <= pi || throw(ArgumentError("dual-active-bridge phase shift must not exceed pi radians"))
    return ratio * v1 * v2 * phase_shift * (1.0 - abs(phase_shift) / pi) /
        (angular_frequency * leakage)
end

function controlled_rectifier_average_voltage_v(
    family::ConverterSystemFamily,
    source_voltage_v::Real,
    firing_angle_rad::Real,
)
    voltage = Float64(source_voltage_v)
    angle = Float64(firing_angle_rad)
    isfinite(voltage) && voltage > 0.0 || throw(ArgumentError("rectifier source voltage must be finite and positive"))
    isfinite(angle) && 0.0 <= angle <= pi || throw(ArgumentError("rectifier firing angle must be between zero and pi"))
    family === SinglePhaseThyristorBridge && return 2.0 * voltage * cos(angle) / pi
    family in (ThreePhaseThyristorBridge, MultipulseThyristorBridge) &&
        return 3.0 * sqrt(2.0) * voltage * cos(angle) / pi
    throw(ArgumentError("the selected converter family has no admitted controlled-rectifier average relation"))
end

function interleaved_carrier_phases_rad(channel_count::Integer; initial_phase_rad::Real=0.0)
    count = Int(channel_count)
    2 <= count <= 8 || throw(ArgumentError("interleaved carrier count must be two through eight"))
    initial = Float64(initial_phase_rad)
    isfinite(initial) || throw(ArgumentError("initial carrier phase must be finite"))
    return [mod(initial + 2.0 * pi * (channel - 1) / count, 2.0 * pi) for channel in 1:count]
end

function matrix_converter_terminal_map(
    connection::AbstractMatrix{Bool},
    input_voltage_v::AbstractVector{<:Real},
    output_current_a::AbstractVector{<:Real},
)
    size(connection) == (3, 3) || throw(DimensionMismatch("matrix converter connection must be 3x3"))
    length(input_voltage_v) == 3 || throw(DimensionMismatch("matrix converter requires three input voltages"))
    length(output_current_a) == 3 || throw(DimensionMismatch("matrix converter requires three output currents"))
    all(row -> count(connection[row, :]) == 1, axes(connection, 1)) || throw(ArgumentError(
        "every matrix-converter output must connect to exactly one input",
    ))
    voltages = Float64.(input_voltage_v)
    currents = Float64.(output_current_a)
    all(isfinite, voltages) && all(isfinite, currents) || throw(ArgumentError(
        "matrix-converter terminal values must be finite",
    ))
    incidence = Float64.(connection)
    return (
        output_voltage_v = incidence * voltages,
        input_current_a = transpose(incidence) * currents,
    )
end

function converter_stage_energy_residual_w(
    input_power_w::Real,
    output_power_w::Real,
    loss_power_w::Real,
    stored_energy_rate_w::Real,
)
    values = Float64.((input_power_w, output_power_w, loss_power_w, stored_energy_rate_w))
    all(isfinite, values) || throw(ArgumentError("converter energy-balance terms must be finite"))
    values[3] >= 0.0 || throw(ArgumentError("converter loss power must be nonnegative"))
    return values[1] - values[2] - values[3] - values[4]
end

module TACS

using ..Branches
using ..Switches: over16_a8sw_delayed_open_step,
                  over16_switch_tail_current_injection

import ..Branches: EMTElement,
                   backward_euler_companion_supported,
                   stamp!,
                   stamp_conductance!,
                   update!

export TACSControlledSwitch,
       ControlledSwitchDelayedArcState,
       apply_controlled_switch_delayed_arc_transition!,
       AlgebraicControlAssignment,
       ControlTransferFunction,
       ConstantControlSignal,
       SinusoidalControlSignal,
       ControlExpressionProgram,
       ControlExpressionRuntime,
       ControlSystemExecutionState,
       ControlSystemExecutionTraceState,
       ControlSystemDevice,
       SignedControlSignalTerm,
       control_system_output_values,
       control_system_execution_trace_preview,
       control_system_execution_trace_update!,
       initialize_control_system_utilities!,
       initialize_control_function_steady_state!,
       advance_control_system_state!,
       sinusoidal_control_signal_phasor,
       sinusoidal_control_signal_value,
       compile_control_expression,
       evaluate_control_expression!,
       control_system_step!,
       controlled_switch_closed,
       controlled_switch_conductance

# TACS treats magnitudes through ten times its 1e-12 machine-zero as zero.
const CONTROL_SYSTEM_ZERO_TOLERANCE = 1.0e-11

include("control_expressions.jl")
include(joinpath(@__DIR__, "tacs", "delayed_arc.jl"))

mutable struct TACSControlledSwitch <: EMTElement
    a::Int
    b::Int
    control::Base.RefValue{Float64}
    clamp_control::Union{Nothing,Base.RefValue{Float64}}
    threshold::Float64
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    last_control::Float64
    last_conductance::Float64
    unidirectional_latching::Bool
    bidirectional_latching::Bool
    ignition_voltage::Float64
    holding_current::Float64
    deionization_time_s::Float64
    delayed_arc::Union{Nothing,ControlledSwitchDelayedArcState}
    open_duration_s::Float64
    last_voltage::Float64
    last_current::Float64
end

function TACSControlledSwitch(
    a::Int,
    b::Int,
    control::Base.RefValue{Float64};
    threshold::Real=0.0,
    on_conductance::Real=1.0e9,
    off_conductance::Real=0.0,
    unidirectional_latching::Bool=false,
    bidirectional_latching::Bool=false,
    clamp_control::Union{Nothing,Base.RefValue{Float64}}=nothing,
    initially_closed::Union{Nothing,Bool}=nothing,
    ignition_voltage::Real=0.0,
    holding_current::Real=0.0,
    deionization_time_s::Real=0.0,
    delayed_arc::Union{Nothing,ControlledSwitchDelayedArcState}=nothing,
)
    a >= 0 && b >= 0 || throw(ArgumentError("switch nodes must be nonnegative"))
    limit = Float64(threshold)
    on = Float64(on_conductance)
    off = Float64(off_conductance)
    ignition = Float64(ignition_voltage)
    holding = Float64(holding_current)
    deionization = Float64(deionization_time_s)
    isfinite(limit) || throw(ArgumentError("threshold must be finite"))
    isfinite(on) && isfinite(off) && on >= 0.0 && off >= 0.0 ||
        throw(ArgumentError("switch conductances must be finite and nonnegative"))
    on >= off || throw(ArgumentError("on_conductance must be at least off_conductance"))
    isfinite(ignition) && ignition >= 0.0 ||
        throw(ArgumentError("ignition_voltage must be finite and nonnegative"))
    isfinite(holding) && holding >= 0.0 ||
        throw(ArgumentError("holding_current must be finite and nonnegative"))
    isfinite(deionization) && deionization >= 0.0 ||
        throw(ArgumentError("deionization_time_s must be finite and nonnegative"))
    !(unidirectional_latching && bidirectional_latching) ||
        throw(ArgumentError("a controlled switch cannot use both latching polarities"))
    clamp_control === nothing || isfinite(clamp_control[]) ||
        throw(ArgumentError("clamp control must be finite"))
    latching = unidirectional_latching || bidirectional_latching
    closed = initially_closed === nothing ?
        (latching ? control[] > limit : control[] >= limit) :
        initially_closed
    conductance = closed ? on : off
    return TACSControlledSwitch(
        a,
        b,
        control,
        clamp_control,
        limit,
        on,
        off,
        closed,
        control[],
        conductance,
        unidirectional_latching,
        bidirectional_latching,
        ignition,
        holding,
        deionization,
        delayed_arc,
        closed ? 0.0 : Inf,
        0.0,
        0.0,
    )
end

function TACSControlledSwitch(
    a::Int,
    b::Int,
    control::Base.RefValue{T};
    threshold::Real=0.0,
    on_conductance::Real=1.0e9,
    off_conductance::Real=0.0,
    unidirectional_latching::Bool=false,
    bidirectional_latching::Bool=false,
    clamp_control::Union{Nothing,Base.RefValue{Float64}}=nothing,
    initially_closed::Union{Nothing,Bool}=nothing,
    ignition_voltage::Real=0.0,
    holding_current::Real=0.0,
    deionization_time_s::Real=0.0,
    delayed_arc::Union{Nothing,ControlledSwitchDelayedArcState}=nothing,
) where {T<:Real}
    return TACSControlledSwitch(
        a,
        b,
        Ref(Float64(control[]));
        threshold = threshold,
        on_conductance = on_conductance,
        off_conductance = off_conductance,
        unidirectional_latching = unidirectional_latching,
        bidirectional_latching = bidirectional_latching,
        clamp_control = clamp_control,
        initially_closed = initially_closed,
        ignition_voltage = ignition_voltage,
        holding_current = holding_current,
        deionization_time_s = deionization_time_s,
        delayed_arc = delayed_arc,
    )
end

controlled_switch_closed(s::TACSControlledSwitch)::Bool = s.closed

controlled_switch_conductance(s::TACSControlledSwitch)::Float64 =
    controlled_switch_closed(s) ? s.on_conductance : s.off_conductance

# The controlled switch stamps an instantaneous conductance and an analytic
# delayed-arc source; neither history term changes its equation between the
# trapezoidal and backward-Euler event substeps. Transaction rollback owns the
# mutable calendar and tail state during localization probes.
backward_euler_companion_supported(::TACSControlledSwitch) = true

function sync_controlled_switch!(s::TACSControlledSwitch, dt_s::Real=0.0)
    dt = Float64(dt_s)
    isfinite(dt) && dt >= 0.0 ||
        throw(ArgumentError("controlled-switch synchronization step must be finite and nonnegative"))
    s.last_control = s.control[]
    if s.delayed_arc !== nothing && s.delayed_arc.tail_active
        s.closed = false
        s.last_conductance = s.off_conductance
        return s
    end
    previously_closed = s.closed
    if s.clamp_control !== nothing
        clamp_value = s.clamp_control[]
        if clamp_value > CONTROL_SYSTEM_ZERO_TOLERANCE
            s.closed = true
            s.delayed_arc === nothing ||
                (s.delayed_arc.opening_requested = false)
            s.open_duration_s = 0.0
            s.last_conductance = s.on_conductance
            return s
        elseif clamp_value < -CONTROL_SYSTEM_ZERO_TOLERANCE
            s.closed = false
            s.open_duration_s = 0.0
            _defer_controlled_switch_arc_opening!(s, previously_closed)
            s.last_conductance = controlled_switch_conductance(s)
            return s
        end
    end
    if s.unidirectional_latching || s.bidirectional_latching
        !s.closed && isfinite(s.open_duration_s) &&
            (s.open_duration_s += dt)
        ignition_voltage =
            s.bidirectional_latching ? abs(s.last_voltage) : s.last_voltage
        if !s.closed &&
           s.open_duration_s > s.deionization_time_s &&
           s.last_control > s.threshold &&
           ignition_voltage >= s.ignition_voltage
            s.closed = true
            s.open_duration_s = 0.0
        end
    else
        s.closed = s.last_control >= s.threshold
    end
    _defer_controlled_switch_arc_opening!(s, previously_closed)
    s.last_conductance = controlled_switch_conductance(s)
    return s
end

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    s::TACSControlledSwitch,
    t::Float64,
    dt::Float64,
)
    sync_controlled_switch!(s, dt)
    _advance_controlled_switch_delayed_arc!(s, t, dt)
    s.last_conductance = controlled_switch_conductance(s)
    stamp_conductance!(y, s.a, s.b, s.last_conductance)
    _stamp_controlled_switch_arc_tail!(rhs, s, t, dt)
    return nothing
end

function update!(
    s::TACSControlledSwitch,
    voltages::AbstractVector{Float64},
    dt::Float64,
)
    previously_closed = s.closed
    prior_current = s.last_current
    from_voltage = s.a == 0 ? 0.0 : voltages[s.a]
    to_voltage = s.b == 0 ? 0.0 : voltages[s.b]
    s.last_voltage = from_voltage - to_voltage
    conduction_current = s.last_conductance * s.last_voltage
    s.last_current = conduction_current
    if s.clamp_control !== nothing &&
       abs(s.clamp_control[]) > CONTROL_SYSTEM_ZERO_TOLERANCE
        return nothing
    end
    if s.unidirectional_latching
        if s.closed && s.last_current < s.holding_current
            s.closed = false
            s.open_duration_s = 0.0
        end
    elseif s.bidirectional_latching
        gate_holds_closed = s.control[] > CONTROL_SYSTEM_ZERO_TOLERANCE
        if s.closed && !gate_holds_closed &&
           (s.last_current * prior_current < 0.0 ||
            abs(s.last_current) < s.holding_current)
            s.closed = false
            s.open_duration_s = 0.0
        end
    end
    _defer_controlled_switch_arc_opening!(s, previously_closed)
    if s.delayed_arc !== nothing
        s.last_current -= s.delayed_arc.tail_current_a
        s.delayed_arc.absolute_tail_energy_j +=
            abs(s.last_voltage * s.delayed_arc.tail_current_a) * dt
    end
    return nothing
end

struct ConstantControlSignal
    name::Symbol
    value::Float64
end

function ConstantControlSignal(name::Symbol, value::Real)
    signal_value = Float64(value)
    isfinite(signal_value) || throw(ArgumentError("control signal value must be finite"))
    return ConstantControlSignal(name, signal_value)
end

struct SinusoidalControlSignal
    name::Symbol
    amplitude::Float64
    frequency_hz::Float64
    phase_rad::Float64
    start_time_s::Float64
    stop_time_s::Float64
end

function SinusoidalControlSignal(
    name::Symbol,
    amplitude::Real,
    frequency_hz::Real,
    phase_rad::Real;
    start_time_s::Real = 0.0,
    stop_time_s::Real = Inf,
)
    magnitude = Float64(amplitude)
    frequency = Float64(frequency_hz)
    phase = Float64(phase_rad)
    start_time = Float64(start_time_s)
    stop_time = Float64(stop_time_s)
    isfinite(magnitude) || throw(ArgumentError("sinusoidal control amplitude must be finite"))
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("sinusoidal control frequency must be finite and positive"))
    isfinite(phase) || throw(ArgumentError("sinusoidal control phase must be finite"))
    isfinite(start_time) || throw(ArgumentError("sinusoidal control start time must be finite"))
    (isfinite(stop_time) || stop_time == Inf) && stop_time > start_time ||
        throw(ArgumentError("sinusoidal control stop time must follow its start time"))
    return SinusoidalControlSignal(
        name,
        magnitude,
        frequency,
        phase,
        start_time,
        stop_time,
    )
end

sinusoidal_control_signal_phasor(signal::SinusoidalControlSignal) =
    signal.amplitude * cis(signal.phase_rad)

function sinusoidal_control_signal_value(
    signal::SinusoidalControlSignal,
    time_s::Real,
)
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("sinusoidal control time must be finite"))
    if time < signal.start_time_s || time >= signal.stop_time_s ||
       (time == 0.0 && signal.start_time_s >= 0.0)
        return 0.0
    end
    return signal.amplitude * cospi(
        2.0 * signal.frequency_hz * time + signal.phase_rad / pi,
    )
end

struct AlgebraicControlAssignment
    output_name::Symbol
    input_name::Symbol
    gain::Float64
    order::Int
end

function AlgebraicControlAssignment(
    output_name::Symbol,
    input_name::Symbol;
    gain::Real=1.0,
    order::Int=0,
)
    order == 0 ||
        throw(ArgumentError("only zero-order algebraic control assignments are supported"))
    assignment_gain = Float64(gain)
    isfinite(assignment_gain) || throw(ArgumentError("control assignment gain must be finite"))
    return AlgebraicControlAssignment(output_name, input_name, assignment_gain, order)
end

struct SignedControlSignalTerm
    name::Symbol
    polarity::Int
end

function SignedControlSignalTerm(name::Symbol, polarity::Integer)
    sign = Int(polarity)
    sign == 1 || sign == -1 ||
        throw(ArgumentError("control signal term polarity must be +1 or -1"))
    return SignedControlSignalTerm(name, sign)
end

struct ControlTransferFunction
    output_name::Symbol
    input_terms::Vector{SignedControlSignalTerm}
    gain::Float64
    order::Int
    numerator_coefficients::Vector{Float64}
    denominator_coefficients::Vector{Float64}
    lower_limit::Union{Missing,Float64}
    upper_limit::Union{Missing,Float64}
    lower_limit_signal::Union{Missing,Symbol}
    upper_limit_signal::Union{Missing,Symbol}
end

function ControlTransferFunction(
    output_name::Symbol,
    input_terms::AbstractVector{SignedControlSignalTerm};
    gain::Real=1.0,
    order::Int=0,
    numerator_coefficients::AbstractVector{<:Real}=order == 0 ? [1.0] : Float64[],
    denominator_coefficients::AbstractVector{<:Real}=order == 0 ? [1.0] : Float64[],
    lower_limit::Union{Missing,Real}=missing,
    upper_limit::Union{Missing,Real}=missing,
    lower_limit_signal::Union{Missing,Symbol}=missing,
    upper_limit_signal::Union{Missing,Symbol}=missing,
)
    0 <= order <= 7 ||
        throw(ArgumentError("control function order must be between zero and seven"))
    !isempty(input_terms) ||
        throw(ArgumentError("control function must have at least one input term"))
    function_gain = Float64(gain)
    isfinite(function_gain) || throw(ArgumentError("control function gain must be finite"))
    numerator = Float64.(numerator_coefficients)
    denominator = Float64.(denominator_coefficients)
    length(numerator) == order + 1 ||
        throw(ArgumentError("control function numerator length must equal order plus one"))
    length(denominator) == order + 1 ||
        throw(ArgumentError("control function denominator length must equal order plus one"))
    all(isfinite, numerator) && all(isfinite, denominator) ||
        throw(ArgumentError("control function coefficients must be finite"))
    any(!iszero, denominator) ||
        throw(ArgumentError("control function denominator must not be identically zero"))
    low = lower_limit === missing ? missing : Float64(lower_limit)
    high = upper_limit === missing ? missing : Float64(upper_limit)
    low === missing || isfinite(low) ||
        throw(ArgumentError("control function lower limit must be finite"))
    high === missing || isfinite(high) ||
        throw(ArgumentError("control function upper limit must be finite"))
    return ControlTransferFunction(
        output_name,
        collect(input_terms),
        function_gain,
        order,
        numerator,
        denominator,
        low,
        high,
        lower_limit_signal,
        upper_limit_signal,
    )
end

mutable struct ControlFunctionRuntimeState
    feedforward_coefficients::Vector{Float64}
    feedback_coefficients::Vector{Float64}
    history_terms::Vector{Float64}
    limiter_active::Bool
    crossed_limit_count::Int
end

struct ControlSystemDevice
    name::Symbol
    device_type::Int
    input_name::Symbol
    tail_signal_names::Vector{Symbol}
    parameter_values::Vector{Float64}
end

function ControlSystemDevice(
    name::Symbol,
    device_type::Integer,
    input_name::Symbol;
    tail_signal_names::AbstractVector{Symbol}=Symbol[],
    parameter_values::AbstractVector{<:Real}=Float64[],
)
    typed_parameters = Float64.(parameter_values)
    all(isfinite, typed_parameters) ||
        throw(ArgumentError("control-system device parameters must be finite"))
    return ControlSystemDevice(
        name,
        Int(device_type),
        input_name,
        collect(tail_signal_names),
        typed_parameters,
    )
end

mutable struct ControlSystemExecutionState
    values::Dict{Symbol,Float64}
    assignments::Vector{AlgebraicControlAssignment}
    functions::Vector{ControlTransferFunction}
    function_states::Vector{ControlFunctionRuntimeState}
    devices::Vector{ControlSystemDevice}
    device_input_history::Dict{Symbol,Float64}
    output_names::Vector{Symbol}
    deltat_s::Float64
    frequency_hz::Float64
end

function _control_polynomial_product_coefficient(
    derivative_order::Int,
    total_order::Int,
    output_order::Int,
)
    coefficient = 0.0
    first_power = max(0, output_order - (total_order - derivative_order))
    last_power = min(derivative_order, output_order)
    for negative_power in first_power:last_power
        positive_power = output_order - negative_power
        coefficient +=
            (-1.0)^negative_power *
            binomial(derivative_order, negative_power) *
            binomial(total_order - derivative_order, positive_power)
    end
    return coefficient
end

function _bilinear_control_coefficients(
    function_row::ControlTransferFunction,
    deltat_s::Float64,
)
    order = function_row.order
    numerator = zeros(Float64, order + 1)
    denominator = zeros(Float64, order + 1)
    derivative_scale = 2.0 / deltat_s
    for derivative_order in 0:order
        scale = derivative_scale^derivative_order
        numerator_value =
            function_row.gain *
            function_row.numerator_coefficients[derivative_order + 1] *
            scale
        denominator_value =
            function_row.denominator_coefficients[derivative_order + 1] *
            scale
        for output_order in 0:order
            basis = _control_polynomial_product_coefficient(
                derivative_order,
                order,
                output_order,
            )
            numerator[output_order + 1] += numerator_value * basis
            denominator[output_order + 1] += denominator_value * basis
        end
    end
    abs(denominator[1]) > eps(Float64) ||
        throw(ArgumentError(
            "control function $(function_row.output_name) has a singular bilinear denominator",
        ))
    numerator ./= denominator[1]
    denominator ./= denominator[1]
    return numerator, denominator
end

function ControlFunctionRuntimeState(
    function_row::ControlTransferFunction,
    deltat_s::Float64,
)
    feedforward, feedback =
        _bilinear_control_coefficients(function_row, deltat_s)
    return ControlFunctionRuntimeState(
        feedforward,
        feedback,
        zeros(Float64, function_row.order),
        false,
        0,
    )
end

mutable struct ControlSystemExecutionTraceState
    stage_names::Vector{Symbol}
    stage_codes::Vector{Int}
    step_indices::Vector{Int}
    times_s::Vector{Float64}
    time_steps_s::Vector{Float64}
    network_voltage_samples::Vector{Float64}
    active_variable_counts::Vector{Int}
    output_counts::Vector{Int}
    machine_counts::Vector{Int}
    machine_output_counts::Vector{Int}
end

ControlSystemExecutionTraceState() = ControlSystemExecutionTraceState(
    Symbol[],
    Int[],
    Int[],
    Float64[],
    Float64[],
    Float64[],
    Int[],
    Int[],
    Int[],
    Int[],
)

function _control_system_trace_stage_names(stage_names)
    return [Symbol(name) for name in stage_names]
end

function _control_system_trace_ints(name::AbstractString, values, count::Int)
    output = Int.(values)
    length(output) == count ||
        throw(ArgumentError("control system execution trace $name length must match stage count"))
    all(value -> value >= 0, output) ||
        throw(ArgumentError("control system execution trace $name values must be nonnegative"))
    return output
end

function _control_system_trace_floats(name::AbstractString, values, count::Int)
    output = Float64.(values)
    length(output) == count ||
        throw(ArgumentError("control system execution trace $name length must match stage count"))
    all(isfinite, output) ||
        throw(ArgumentError("control system execution trace $name values must be finite"))
    return output
end

function _control_system_trace_is_nondecreasing(values::AbstractVector{<:Real})
    length(values) <= 1 && return true
    return all(values[index + 1] >= values[index] for index in 1:(length(values) - 1))
end

function control_system_execution_trace_preview(
    stage_names;
    stage_codes,
    step_indices,
    times_s,
    time_steps_s,
    network_voltage_samples,
    active_variable_counts,
    output_counts,
    machine_counts,
    machine_output_counts,
)
    stages = _control_system_trace_stage_names(stage_names)
    count = length(stages)
    count > 0 || throw(ArgumentError("control system execution trace must contain at least one stage"))
    codes = _control_system_trace_ints("stage_codes", stage_codes, count)
    steps = _control_system_trace_ints("step_indices", step_indices, count)
    times = _control_system_trace_floats("times_s", times_s, count)
    time_steps = _control_system_trace_floats("time_steps_s", time_steps_s, count)
    voltages = _control_system_trace_floats("network_voltage_samples", network_voltage_samples, count)
    active_counts = _control_system_trace_ints("active_variable_counts", active_variable_counts, count)
    output_counts_vector = _control_system_trace_ints("output_counts", output_counts, count)
    machine_counts_vector = _control_system_trace_ints("machine_counts", machine_counts, count)
    machine_output_counts_vector =
        _control_system_trace_ints("machine_output_counts", machine_output_counts, count)
    all(value -> value > 0.0, time_steps) ||
        throw(ArgumentError("control system execution trace time steps must be positive"))

    return (
        source = :control_system_execution_trace,
        outcome = :state_mutation,
        stage_names = stages,
        stage_codes = codes,
        step_indices = steps,
        times_s = times,
        time_steps_s = time_steps,
        network_voltage_samples = voltages,
        active_variable_counts = active_counts,
        output_counts = output_counts_vector,
        machine_counts = machine_counts_vector,
        machine_output_counts = machine_output_counts_vector,
        stage_count = count,
        time_nondecreasing = _control_system_trace_is_nondecreasing(times),
        step_index_nondecreasing = _control_system_trace_is_nondecreasing(steps),
        control_system_present = any(>(0), active_counts) && any(>(0), output_counts_vector),
        machine_solution_present =
            any(>(0), machine_counts_vector) && any(>(0), machine_output_counts_vector),
        trace_state_mutated = false,
        execution_trace_replayed = false,
        full_control_system_equation_execution = false,
        full_machine_equation_solution = false,
        full_output_reference_interpreter = false,
        deferred_effects = (
            :control_system_equation_execution,
            :machine_equation_solution,
            :output_reference_interpreter,
            :report_file_generation,
        ),
    )
end

function control_system_execution_trace_update!(
    state::ControlSystemExecutionTraceState,
    stage_names;
    kwargs...,
)
    preview = control_system_execution_trace_preview(stage_names; kwargs...)
    empty!(state.stage_names)
    empty!(state.stage_codes)
    empty!(state.step_indices)
    empty!(state.times_s)
    empty!(state.time_steps_s)
    empty!(state.network_voltage_samples)
    empty!(state.active_variable_counts)
    empty!(state.output_counts)
    empty!(state.machine_counts)
    empty!(state.machine_output_counts)
    append!(state.stage_names, preview.stage_names)
    append!(state.stage_codes, preview.stage_codes)
    append!(state.step_indices, preview.step_indices)
    append!(state.times_s, preview.times_s)
    append!(state.time_steps_s, preview.time_steps_s)
    append!(state.network_voltage_samples, preview.network_voltage_samples)
    append!(state.active_variable_counts, preview.active_variable_counts)
    append!(state.output_counts, preview.output_counts)
    append!(state.machine_counts, preview.machine_counts)
    append!(state.machine_output_counts, preview.machine_output_counts)
    return merge(
        preview,
        (
            trace_state_mutated = true,
            execution_trace_replayed = true,
        ),
    )
end

function ControlSystemExecutionState(
    signals::AbstractVector{ConstantControlSignal},
    assignments::AbstractVector{AlgebraicControlAssignment};
    functions::AbstractVector{ControlTransferFunction}=ControlTransferFunction[],
    devices::AbstractVector{ControlSystemDevice}=ControlSystemDevice[],
    output_names::AbstractVector{Symbol}=Symbol[],
    deltat_s::Real,
    frequency_hz::Real=0.0,
)
    step_s = Float64(deltat_s)
    frequency = Float64(frequency_hz)
    step_s > 0.0 && isfinite(step_s) ||
        throw(ArgumentError("deltat_s must be finite and positive"))
    isfinite(frequency) || throw(ArgumentError("frequency_hz must be finite"))
    values = Dict{Symbol,Float64}()
    for signal in signals
        values[signal.name] = signal.value
    end
    return ControlSystemExecutionState(
        values,
        collect(assignments),
        collect(functions),
        ControlFunctionRuntimeState[
            ControlFunctionRuntimeState(function_row, step_s)
            for function_row in functions
        ],
        collect(devices),
        Dict{Symbol,Float64}(),
        collect(output_names),
        step_s,
        frequency,
    )
end

const _CONTROL_SYSTEM_UTILITY_VALUES = (
    :TIMEX,
    :ISTEP,
    :DELTAT,
    :FREQHZ,
    :OMEGAR,
    :ZERO,
    :PLUS1,
    :MINUS1,
    :UNITY,
    :INFNTY,
    :PI,
)

function _write_control_system_utilities!(
    values::Dict{Symbol,Float64},
    step::Int,
    time_s::Float64,
    deltat_s::Float64,
    frequency_hz::Float64,
)
    values[:TIMEX] = time_s
    values[:ISTEP] = Float64(step)
    values[:DELTAT] = deltat_s
    values[:FREQHZ] = frequency_hz
    values[:OMEGAR] = 2.0 * pi * frequency_hz
    values[:ZERO] = 0.0
    values[:PLUS1] = 1.0
    values[:MINUS1] = -1.0
    values[:UNITY] = 1.0
    values[:INFNTY] = 1.0e20
    values[:PI] = pi
    return values
end

function initialize_control_system_utilities!(
    state::ControlSystemExecutionState,
    step::Int,
    time_s::Real,
)
    step >= 0 || throw(ArgumentError("step must be nonnegative"))
    time_value = Float64(time_s)
    isfinite(time_value) || throw(ArgumentError("time_s must be finite"))
    _write_control_system_utilities!(
        state.values,
        step,
        time_value,
        state.deltat_s,
        state.frequency_hz,
    )
    return state
end

function _control_system_function_history_value(
    runtime::ControlFunctionRuntimeState,
)
    return isempty(runtime.history_terms) ? 0.0 : first(runtime.history_terms)
end

function _update_control_function_history!(
    runtime::ControlFunctionRuntimeState,
    input_value::Float64,
    output_value::Float64;
    reset_for_limit::Bool,
)
    order = length(runtime.history_terms)
    order == 0 && return runtime
    if reset_for_limit
        runtime.history_terms[end] =
            runtime.feedforward_coefficients[end] * input_value -
            runtime.feedback_coefficients[end] * output_value
        for index in (order - 1):-1:1
            runtime.history_terms[index] = (
                runtime.feedforward_coefficients[index + 1] * input_value -
                runtime.feedback_coefficients[index + 1] * output_value +
                runtime.history_terms[index + 1]
            ) / 2.0
        end
    else
        for index in 1:(order - 1)
            runtime.history_terms[index] =
                runtime.feedforward_coefficients[index + 1] * input_value -
                runtime.feedback_coefficients[index + 1] * output_value -
                runtime.history_terms[index] +
                runtime.history_terms[index + 1]
        end
        runtime.history_terms[end] =
            runtime.feedforward_coefficients[end] * input_value -
            runtime.feedback_coefficients[end] * output_value
    end
    return runtime
end

function _control_system_function_limit(
    values::Dict{Symbol,Float64},
    numeric_limit::Union{Missing,Float64},
    signal_name::Union{Missing,Symbol},
    default::Float64,
)
    signal_name === missing &&
        return numeric_limit === missing ? default : numeric_limit
    haskey(values, signal_name) ||
        throw(ArgumentError("missing control function limit signal $signal_name"))
    value = values[signal_name]
    isfinite(value) ||
        throw(ArgumentError("control function limit signal $signal_name is nonfinite"))
    return value
end

function _control_system_function_input(
    values::Dict{Symbol,Float64},
    function_row::ControlTransferFunction,
)
    value = 0.0
    for term in function_row.input_terms
        haskey(values, term.name) ||
            throw(ArgumentError("missing control function input $(term.name)"))
        value += term.polarity * values[term.name]
    end
    return value
end

function _control_system_function_solution!(
    state::ControlSystemExecutionState,
    direct_gains::AbstractVector{<:Number},
    history_values::AbstractVector{<:Number},
)
    count = length(state.functions)
    length(direct_gains) == count && length(history_values) == count ||
        throw(ArgumentError("control function solution vectors must match function count"))
    count == 0 && return Float64[], Dict{Int,Float64}()
    scalar_type = promote_type(eltype(direct_gains), eltype(history_values))
    output_indices = Dict(
        function_row.output_name => index
        for (index, function_row) in enumerate(state.functions)
    )
    fixed_values = Dict{Int,scalar_type}()
    solution = zeros(scalar_type, count)
    for _ in 1:(count + 1)
        matrix = zeros(scalar_type, count, count)
        right_hand_side = zeros(scalar_type, count)
        for index in 1:count
            if haskey(fixed_values, index)
                matrix[index, index] = 1.0
                right_hand_side[index] = fixed_values[index]
                continue
            end
            function_row = state.functions[index]
            matrix[index, index] = one(scalar_type)
            right_hand_side[index] = history_values[index]
            direct_gain = direct_gains[index]
            for term in function_row.input_terms
                coupled_index = get(output_indices, term.name, 0)
                if coupled_index == 0
                    haskey(state.values, term.name) ||
                        throw(ArgumentError("missing control function input $(term.name)"))
                    right_hand_side[index] +=
                        direct_gain * term.polarity * state.values[term.name]
                else
                    matrix[index, coupled_index] -=
                        direct_gain * term.polarity
                end
            end
        end
        solution .= matrix \ right_hand_side
        all(isfinite, solution) ||
            throw(ArgumentError("control function simultaneous solution is nonfinite"))
        scalar_type <: Real || break
        solved_values = copy(state.values)
        for index in 1:count
            solved_values[state.functions[index].output_name] = solution[index]
        end
        newly_fixed = false
        for index in 1:count
            haskey(fixed_values, index) && continue
            function_row = state.functions[index]
            lower = _control_system_function_limit(
                solved_values,
                function_row.lower_limit,
                function_row.lower_limit_signal,
                -Inf,
            )
            upper = _control_system_function_limit(
                solved_values,
                function_row.upper_limit,
                function_row.upper_limit_signal,
                Inf,
            )
            if lower > upper
                state.function_states[index].crossed_limit_count += 1
                fixed_values[index] = scalar_type(upper)
                newly_fixed = true
            elseif solution[index] < lower
                fixed_values[index] = scalar_type(lower)
                newly_fixed = true
            elseif solution[index] > upper
                fixed_values[index] = scalar_type(upper)
                newly_fixed = true
            end
        end
        newly_fixed || break
    end
    return solution, fixed_values
end

function _advance_control_system_functions!(
    state::ControlSystemExecutionState,
)
    count = length(state.functions)
    count == 0 && return state
    direct_gains = Float64[
        runtime.feedforward_coefficients[1] for runtime in state.function_states
    ]
    history_values = Float64[
        _control_system_function_history_value(runtime)
        for runtime in state.function_states
    ]
    solution, fixed_values = _control_system_function_solution!(
        state,
        direct_gains,
        history_values,
    )
    for index in 1:count
        function_row = state.functions[index]
        runtime = state.function_states[index]
        value = get(fixed_values, index, solution[index])
        state.values[function_row.output_name] = value
        runtime.limiter_active = haskey(fixed_values, index)
    end
    for index in 1:count
        function_row = state.functions[index]
        runtime = state.function_states[index]
        input_value = _control_system_function_input(state.values, function_row)
        _update_control_function_history!(
            runtime,
            input_value,
            state.values[function_row.output_name];
            reset_for_limit = runtime.limiter_active,
        )
    end
    return state
end

function initialize_control_function_steady_state!(
    state::ControlSystemExecutionState,
)
    initialize_control_system_utilities!(state, 0, 0.0)
    for assignment in state.assignments
        haskey(state.values, assignment.input_name) ||
            throw(ArgumentError("missing control assignment input $(assignment.input_name)"))
        state.values[assignment.output_name] =
            assignment.gain * state.values[assignment.input_name]
    end
    direct_gains = Vector{Float64}(undef, length(state.functions))
    for (index, function_row) in enumerate(state.functions)
        denominator = function_row.denominator_coefficients[1]
        abs(denominator) > eps(Float64) || throw(ArgumentError(
            "control function $(function_row.output_name) has no finite DC equilibrium",
        ))
        direct_gains[index] =
            function_row.gain *
            function_row.numerator_coefficients[1] /
            denominator
    end
    solution, fixed_values = _control_system_function_solution!(
        state,
        direct_gains,
        zeros(Float64, length(state.functions)),
    )
    for index in eachindex(state.functions)
        function_row = state.functions[index]
        runtime = state.function_states[index]
        output = get(fixed_values, index, solution[index])
        state.values[function_row.output_name] = output
        runtime.limiter_active = haskey(fixed_values, index)
    end
    for index in eachindex(state.functions)
        function_row = state.functions[index]
        runtime = state.function_states[index]
        input = _control_system_function_input(state.values, function_row)
        _update_control_function_history!(
            runtime,
            input,
            state.values[function_row.output_name];
            reset_for_limit = true,
        )
    end
    return state
end

_control_device_parameter(device::ControlSystemDevice, index::Int, default::Float64) =
    index <= length(device.parameter_values) ? device.parameter_values[index] : default

function _control_system_device_value!(
    values::Dict{Symbol,Float64},
    input_history::Dict{Symbol,Float64},
    device::ControlSystemDevice,
    deltat_s::Float64,
)
    haskey(values, device.input_name) ||
        throw(ArgumentError("missing control-system device input $(device.input_name)"))
    input_value = values[device.input_name]
    if device.device_type == 58
        gain = _control_device_parameter(device, 1, 1.0)
        reset_value = _control_device_parameter(device, 2, 0.0)
        if !isempty(device.tail_signal_names)
            control_name = device.tail_signal_names[1]
            haskey(values, control_name) ||
                throw(ArgumentError("missing type-58 control signal $control_name"))
            if values[control_name] <= CONTROL_SYSTEM_ZERO_TOLERANCE
                input_history[device.name] = 0.0
                return reset_value
            end
        end
        previous_input = get(input_history, device.name, 0.0)
        previous_output = get(values, device.name, reset_value)
        output_value =
            previous_output + 0.5 * deltat_s * gain * (input_value + previous_input)
        if length(device.parameter_values) >= 3
            lower = device.parameter_values[2]
            upper = device.parameter_values[3]
            if upper >= lower
                output_value = clamp(output_value, lower, upper)
            end
        end
        input_history[device.name] = input_value
        return output_value
    elseif device.device_type == 52
        length(device.tail_signal_names) >= 2 ||
            throw(ArgumentError("type-52 control-system device requires reference and observed signal names"))
        reference_name = device.tail_signal_names[1]
        observed_name = device.tail_signal_names[2]
        haskey(values, reference_name) ||
            throw(ArgumentError("missing type-52 reference signal $(reference_name)"))
        haskey(values, observed_name) ||
            throw(ArgumentError("missing type-52 observed signal $(observed_name)"))
        gain = _control_device_parameter(device, 1, 1.0)
        offset = _control_device_parameter(device, 2, 0.0)
        mode = abs(_control_device_parameter(device, 3, 0.0))
        comparison_value = values[reference_name]
        switching_level = offset + values[observed_name]
        active = mode > 1.5 ?
            comparison_value >= switching_level :
            comparison_value < switching_level
        return active ? gain * input_value : 0.0
    end
    throw(ArgumentError("unsupported control-system device type $(device.device_type)"))
end

function advance_control_system_state!(
    state::ControlSystemExecutionState,
    step::Int,
    time_s::Real,
    ;
    execute_devices::Bool=true,
)
    step >= 0 || throw(ArgumentError("step must be nonnegative"))
    time_value = Float64(time_s)
    isfinite(time_value) || throw(ArgumentError("time_s must be finite"))
    initialize_control_system_utilities!(state, step, time_value)
    for assignment in state.assignments
        haskey(state.values, assignment.input_name) ||
            throw(ArgumentError("missing control assignment input $(assignment.input_name)"))
        value = assignment.gain * state.values[assignment.input_name]
        state.values[assignment.output_name] = value
    end
    _advance_control_system_functions!(state)
    if execute_devices
        for device in state.devices
            if device.device_type in (52, 58)
                value = _control_system_device_value!(
                    state.values,
                    state.device_input_history,
                    device,
                    state.deltat_s,
                )
                state.values[device.name] = value
            end
        end
    end
    return state
end

function control_system_step!(
    state::ControlSystemExecutionState,
    step::Int,
    time_s::Real,
)
    advance_control_system_state!(state, step, time_s)
    time_value = Float64(time_s)
    assignment_outputs = Symbol[assignment.output_name for assignment in state.assignments]
    assignment_values = Float64[state.values[name] for name in assignment_outputs]
    function_outputs = Symbol[function_row.output_name for function_row in state.functions]
    function_values = Float64[state.values[name] for name in function_outputs]
    supported_devices = [device for device in state.devices if device.device_type in (52, 58)]
    device_outputs = Symbol[device.name for device in supported_devices]
    device_values = Float64[state.values[name] for name in device_outputs]
    unsupported_device_types = Int[
        device.device_type for device in state.devices if !(device.device_type in (52, 58))
    ]
    output_values = control_system_output_values(state)
    supported_assignment_complete = length(assignment_outputs) == length(state.assignments)
    supported_function_complete = length(function_outputs) == length(state.functions)
    supported_device_complete =
        length(device_outputs) == length(state.devices) && isempty(unsupported_device_types)
    full_device_execution = !isempty(state.devices) && supported_device_complete
    deferred_effects = Symbol[
        :full_use1_stack_replay,
        :full_xref1_interpreter,
        :full_errstp_equivalence,
        :full_tacs_solvum_coupling,
    ]
    full_device_execution || pushfirst!(deferred_effects, :full_tacs3_device_execution)
    return (
        step = step,
        time_s = time_value,
        deltat_s = state.deltat_s,
        frequency_hz = state.frequency_hz,
        utility_names = _CONTROL_SYSTEM_UTILITY_VALUES,
        assignment_outputs = assignment_outputs,
        assignment_values = assignment_values,
        function_outputs = function_outputs,
        function_values = function_values,
        device_outputs = device_outputs,
        device_values = device_values,
        unsupported_device_types = unsupported_device_types,
        output_names = copy(state.output_names),
        output_values = output_values,
        assignment_count = length(state.assignments) + length(state.functions),
        function_count = length(state.functions),
        device_count = length(state.devices),
        executed_device_count = length(device_outputs),
        output_count = length(output_values),
        control_system_executed = true,
        algebraic_assignments_executed =
            !isempty(state.assignments) || !isempty(state.functions),
        supported_assignment_complete = supported_assignment_complete,
        supported_function_complete = supported_function_complete,
        supported_device_complete = supported_device_complete,
        complete_control_system_executed =
            supported_assignment_complete &&
            supported_function_complete &&
            (isempty(state.devices) || supported_device_complete),
        machine_coupling_executed = false,
        full_tacs3_device_execution = full_device_execution,
        full_tacs_solvum_coupling = false,
        deferred_effects = tuple(deferred_effects...),
    )
end

function control_system_output_values(state::ControlSystemExecutionState)
    values = Vector{Float64}(undef, length(state.output_names))
    for i in eachindex(state.output_names)
        name = state.output_names[i]
        haskey(state.values, name) ||
            throw(ArgumentError("missing control output $(name)"))
        values[i] = state.values[name]
    end
    return values
end

end

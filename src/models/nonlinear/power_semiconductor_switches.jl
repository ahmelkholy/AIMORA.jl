abstract type PowerSemiconductorJunction end

struct DiodeJunction <: PowerSemiconductorJunction end
struct ThyristorJunction <: PowerSemiconductorJunction end
struct TriacJunction <: PowerSemiconductorJunction end
struct GateTurnOffThyristorJunction <: PowerSemiconductorJunction end
struct InsulatedGateBipolarTransistorJunction <: PowerSemiconductorJunction end
struct MetalOxideSemiconductorFieldEffectTransistorJunction <: PowerSemiconductorJunction end

"""Inertial gate-command state with separate turn-on/turn-off delay, minimum off-time dead time, and minimum applied pulse width."""
mutable struct PowerSemiconductorGateDriver
    turn_on_delay_s::Float64
    turn_off_delay_s::Float64
    dead_time_s::Float64
    minimum_pulse_width_s::Float64
    commanded_on::Bool
    applied_on::Bool
    pending_state::Union{Nothing,Bool}
    pending_transition_time_s::Float64
    last_command_time_s::Float64
    last_turn_on_time_s::Float64
    last_turn_off_time_s::Float64
    command_count::Int
    transition_count::Int
    filtered_pulse_count::Int
end

function PowerSemiconductorGateDriver(;
    turn_on_delay_s::Real=0.0,
    turn_off_delay_s::Real=0.0,
    dead_time_s::Real=0.0,
    minimum_pulse_width_s::Real=0.0,
    initially_on::Bool=false,
)
    turn_on_delay = Float64(turn_on_delay_s)
    turn_off_delay = Float64(turn_off_delay_s)
    dead_time = Float64(dead_time_s)
    minimum_pulse_width = Float64(minimum_pulse_width_s)
    all(isfinite, (turn_on_delay, turn_off_delay, dead_time, minimum_pulse_width)) ||
        throw(ArgumentError("gate-driver times must be finite"))
    all(>=(0.0), (turn_on_delay, turn_off_delay, dead_time, minimum_pulse_width)) ||
        throw(ArgumentError("gate-driver times must be nonnegative"))
    return PowerSemiconductorGateDriver(
        turn_on_delay,
        turn_off_delay,
        dead_time,
        minimum_pulse_width,
        initially_on,
        initially_on,
        nothing,
        Inf,
        0.0,
        initially_on ? 0.0 : -Inf,
        initially_on ? -Inf : 0.0,
        0,
        0,
        0,
    )
end

"""Generic piecewise-linear antiparallel-diode parameters; current is positive from the semiconductor cathode-side terminal back toward its anode-side terminal."""
struct AntiparallelDiodeParameters
    forward_voltage_v::Float64
    holding_current_a::Float64
    on_conductance_s::Float64
end

function AntiparallelDiodeParameters(;
    forward_voltage_v::Real=0.0,
    holding_current_a::Real=0.0,
    on_conductance_s::Real=1.0e3,
)
    forward_voltage = Float64(forward_voltage_v)
    holding_current = Float64(holding_current_a)
    on_conductance = Float64(on_conductance_s)
    isfinite(forward_voltage) && forward_voltage >= 0.0 || throw(ArgumentError(
        "antiparallel-diode forward voltage must be finite and nonnegative",
    ))
    isfinite(holding_current) && holding_current >= 0.0 || throw(ArgumentError(
        "antiparallel-diode holding current must be finite and nonnegative",
    ))
    isfinite(on_conductance) && on_conductance > 0.0 || throw(ArgumentError(
        "antiparallel-diode on conductance must be finite and positive",
    ))
    return AntiparallelDiodeParameters(forward_voltage, holding_current, on_conductance)
end

"""Trapezoidal companion state for a series resistance-capacitance snubber connected across a semiconductor device."""
mutable struct SeriesRCSnubber
    resistance_ohm::Float64
    capacitance_f::Float64
    previous_current_a::Float64
    capacitor_voltage_v::Float64
    last_branch_voltage_v::Float64
    last_current_a::Float64
    last_resistor_loss_w::Float64
    dissipated_energy_j::Float64
end

function SeriesRCSnubber(
    resistance_ohm::Real,
    capacitance_f::Real;
    initial_current_a::Real=0.0,
    initial_capacitor_voltage_v::Real=0.0,
)
    resistance = Float64(resistance_ohm)
    capacitance = Float64(capacitance_f)
    current = Float64(initial_current_a)
    capacitor_voltage = Float64(initial_capacitor_voltage_v)
    isfinite(resistance) && resistance > 0.0 ||
        throw(ArgumentError("snubber resistance must be finite and positive"))
    isfinite(capacitance) && capacitance > 0.0 ||
        throw(ArgumentError("snubber capacitance must be finite and positive"))
    isfinite(current) && isfinite(capacitor_voltage) || throw(ArgumentError(
        "snubber initial current and capacitor voltage must be finite",
    ))
    initial_loss = resistance * current^2
    return SeriesRCSnubber(
        resistance,
        capacitance,
        current,
        capacitor_voltage,
        capacitor_voltage + resistance * current,
        current,
        initial_loss,
        0.0,
    )
end

"""Typed terminal, path-state, loss, and stored-energy snapshot for one accepted semiconductor state."""
struct PowerSemiconductorTerminalState
    device_kind::Symbol
    forward_conducting::Bool
    reverse_diode_conducting::Bool
    gate_commanded_on::Bool
    gate_applied_on::Bool
    terminal_voltage_v::Float64
    terminal_current_a::Float64
    forward_current_a::Float64
    reverse_diode_current_a::Float64
    snubber_current_a::Float64
    semiconductor_loss_w::Float64
    snubber_resistor_loss_w::Float64
    semiconductor_dissipated_energy_j::Float64
    snubber_dissipated_energy_j::Float64
    snubber_capacitor_energy_j::Float64
end

"""
Evidence-bounded semiconductor switch whose unchanged baseline owns piecewise-linear diode, thyristor, IGBT, and MOSFET behavior. Explicit optional fidelity components add TRIAC/GTO semantics, recovered charge, nonlinear junction charge, tail current, event-energy maps, and passive lumped thermal state without changing the baseline constructor or trace.

Validity requires an EMT timestep fine enough to resolve the external circuit and every admitted dynamic component. Generic extended parameters remain AIMORA-authored physical contracts and never establish manufacturer prediction, arbitrary compact-model compatibility, standard conformance, ATP/PSCAD equivalence, or certification.
"""
mutable struct PowerSemiconductorSwitch{K<:PowerSemiconductorJunction} <: EMTElement
    a::Int
    b::Int
    threshold_v::Float64
    holding_current::Float64
    on_conductance::Float64
    off_conductance::Float64
    closed::Bool
    last_voltage::Float64
    last_current::Float64
    last_conductance::Float64
    forward_voltage_drop_v::Float64
    gate_driver::Union{Nothing,PowerSemiconductorGateDriver}
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}
    snubber::Union{Nothing,SeriesRCSnubber}
    reverse_diode_conducting::Bool
    event_localization_enabled::Bool
    last_evaluation_time_s::Float64
    last_history_current_a::Float64
    last_forward_current_a::Float64
    last_reverse_diode_current_a::Float64
    last_snubber_current_a::Float64
    last_semiconductor_loss_w::Float64
    previous_semiconductor_loss_w::Float64
    semiconductor_dissipated_energy_j::Float64
    topology_transition_count::Int
    last_transition_time_s::Float64
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}
    conduction_direction::Int8
    gate_turn_off_current_limit_a::Float64
    gate_turn_off_policy::Symbol
    gate_turn_off_disposition::Symbol
end

"""Piecewise-linear naturally commutated diode specialization."""
const DiodeValveSwitch = PowerSemiconductorSwitch{DiodeJunction}
"""Piecewise-linear gate-triggered, holding-current-commutated thyristor specialization."""
const ThyristorValveSwitch = PowerSemiconductorSwitch{ThyristorJunction}
"""Bidirectional latching TRIAC composed by the accepted semiconductor terminal owner."""
const TriacSwitch = PowerSemiconductorSwitch{TriacJunction}
"""Forward latching gate-turn-off thyristor with an explicit interruptible-current policy."""
const GateTurnOffThyristorSwitch = PowerSemiconductorSwitch{GateTurnOffThyristorJunction}
"""Piecewise-linear forward-channel IGBT specialization with an optional antiparallel diode."""
const IGBTSwitch = PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}
"""Piecewise-linear bidirectional-on-channel MOSFET specialization with an optional body-diode path."""
const MOSFETSwitch = PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}

power_semiconductor_kind(::PowerSemiconductorSwitch{DiodeJunction}) = :diode
power_semiconductor_kind(::PowerSemiconductorSwitch{ThyristorJunction}) = :thyristor
power_semiconductor_kind(::PowerSemiconductorSwitch{TriacJunction}) = :triac
power_semiconductor_kind(::PowerSemiconductorSwitch{GateTurnOffThyristorJunction}) = :gto
power_semiconductor_kind(::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}) = :igbt
power_semiconductor_kind(::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}) = :mosfet

function _power_semiconductor_switch(
    ::Type{K},
    a::Int,
    b::Int;
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    gate_driver::Union{Nothing,PowerSemiconductorGateDriver}=nothing,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
    gate_turn_off_current_limit_a::Real=Inf,
    gate_turn_off_policy::Symbol=:not_applicable,
) where {K<:PowerSemiconductorJunction}
    a >= 0 && b >= 0 ||
        throw(ArgumentError("power-semiconductor nodes must be nonnegative"))
    a != b || throw(ArgumentError("power-semiconductor terminals must be distinct"))
    threshold = Float64(threshold_v)
    forward_drop = Float64(forward_voltage_drop_v)
    holding = Float64(holding_current)
    on = Float64(on_conductance)
    off = Float64(off_conductance)
    isfinite(threshold) && threshold >= 0.0 || throw(ArgumentError(
        "power-semiconductor turn-on threshold must be finite and nonnegative",
    ))
    isfinite(forward_drop) && forward_drop >= 0.0 || throw(ArgumentError(
        "power-semiconductor forward voltage drop must be finite and nonnegative",
    ))
    isfinite(holding) && holding >= 0.0 || throw(ArgumentError(
        "power-semiconductor holding current must be finite and nonnegative",
    ))
    isfinite(on) && on > 0.0 || throw(ArgumentError(
        "power-semiconductor on conductance must be finite and positive",
    ))
    isfinite(off) && off >= 0.0 && off <= on || throw(ArgumentError(
        "power-semiconductor off conductance must be finite, nonnegative, and no greater than on conductance",
    ))
    K === DiodeJunction && gate_driver !== nothing && throw(ArgumentError(
        "a diode junction cannot own a gate driver",
    ))
    K !== DiodeJunction && gate_driver === nothing && throw(ArgumentError(
        "a controlled power semiconductor requires a gate driver",
    ))
    if K === TriacJunction
        antiparallel_diode === nothing || throw(ArgumentError(
            "a TRIAC owns its two oriented latching channels and cannot add a baseline antiparallel diode",
        ))
    end
    turn_off_current_limit = Float64(gate_turn_off_current_limit_a)
    if K === GateTurnOffThyristorJunction
        isfinite(turn_off_current_limit) && turn_off_current_limit >= 0.0 ||
            throw(ArgumentError(
                "GTO interruptible-current limit must be finite and nonnegative",
            ))
        gate_turn_off_policy in (:refuse, :remain_latched) || throw(ArgumentError(
            "GTO turn-off policy must be :refuse or :remain_latched",
        ))
    else
        turn_off_current_limit == Inf || throw(ArgumentError(
            "only a GTO may declare an interruptible-current limit",
        ))
        gate_turn_off_policy === :not_applicable || throw(ArgumentError(
            "only a GTO may declare a gate-turn-off policy",
        ))
    end
    if extended_fidelity !== nothing
        recovery = extended_fidelity.recovered_charge
        recovery === nothing || K === DiodeJunction || antiparallel_diode !== nothing ||
            throw(ArgumentError(
                "recovered charge requires a diode junction or an explicit antiparallel diode path",
            ))
        tail = extended_fidelity.turn_off_tail
        tail === nothing || K === InsulatedGateBipolarTransistorJunction ||
            K === GateTurnOffThyristorJunction || throw(ArgumentError(
            "turn-off tail fidelity is admitted only for IGBT or GTO junctions",
        ))
    end
    if K === InsulatedGateBipolarTransistorJunction ||
       K === MetalOxideSemiconductorFieldEffectTransistorJunction
        gate_driver.applied_on == initially_closed || throw(ArgumentError(
            "initial gated-channel state must match the applied gate-driver state",
        ))
    end
    conductance = initially_closed ? on : off
    return PowerSemiconductorSwitch{K}(
        a,
        b,
        threshold,
        holding,
        on,
        off,
        initially_closed,
        0.0,
        0.0,
        conductance,
        forward_drop,
        gate_driver,
        antiparallel_diode,
        snubber,
        false,
        false,
        0.0,
        initially_closed ? -on * forward_drop : 0.0,
        0.0,
        0.0,
        snubber === nothing ? 0.0 : snubber.last_current_a,
        0.0,
        0.0,
        0.0,
        0,
        0.0,
        extended_fidelity,
        initially_closed ? Int8(1) : Int8(0),
        turn_off_current_limit,
        gate_turn_off_policy,
        :not_requested,
    )
end

function PowerSemiconductorSwitch{DiodeJunction}(
    a::Int,
    b::Int;
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e9,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        DiodeJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{ThyristorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        ThyristorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{TriacJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    initial_direction::Integer=initially_closed ? 1 : 0,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    direction = Int(initial_direction)
    direction in (-1, 0, 1) || throw(ArgumentError(
        "TRIAC initial direction must be -1, 0, or 1",
    ))
    initially_closed == (direction != 0) || throw(ArgumentError(
        "TRIAC initial direction must be nonzero exactly when initially closed",
    ))
    device = _power_semiconductor_switch(
        TriacJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        snubber,
        extended_fidelity,
    )
    device.conduction_direction = Int8(direction)
    device.last_history_current_a = initially_closed ?
        -device.on_conductance * device.forward_voltage_drop_v * direction : 0.0
    return device
end

function PowerSemiconductorSwitch{GateTurnOffThyristorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    interruptible_current_a::Real,
    turn_off_policy::Symbol=:refuse,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        GateTurnOffThyristorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
        gate_turn_off_current_limit_a=interruptible_current_a,
        gate_turn_off_policy=turn_off_policy,
    )
end

function PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    threshold_v::Real=0.0,
    forward_voltage_drop_v::Real=0.0,
    holding_current::Real=0.0,
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        InsulatedGateBipolarTransistorJunction,
        a,
        b;
        threshold_v,
        forward_voltage_drop_v,
        holding_current,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction}(
    a::Int,
    b::Int;
    gate_driver::PowerSemiconductorGateDriver=PowerSemiconductorGateDriver(),
    on_conductance::Real=1.0e3,
    off_conductance::Real=0.0,
    initially_closed::Bool=false,
    antiparallel_diode::Union{Nothing,AntiparallelDiodeParameters}=nothing,
    snubber::Union{Nothing,SeriesRCSnubber}=nothing,
    extended_fidelity::Union{Nothing,PowerSemiconductorExtendedFidelity}=nothing,
)
    return _power_semiconductor_switch(
        MetalOxideSemiconductorFieldEffectTransistorJunction,
        a,
        b;
        threshold_v=0.0,
        forward_voltage_drop_v=0.0,
        holding_current=0.0,
        on_conductance,
        off_conductance,
        initially_closed,
        gate_driver,
        antiparallel_diode,
        snubber,
        extended_fidelity,
    )
end

function diode_next_closed(
    device::DiodeValveSwitch,
    voltage::Real,
    current::Real,
)::Bool
    return device.closed ? Float64(current) >= device.holding_current :
        Float64(voltage) >= power_semiconductor_forward_turn_on_voltage(device)
end

diode_conductance(device::DiodeValveSwitch)::Float64 =
    device.closed ? device.on_conductance : device.off_conductance

power_semiconductor_forward_turn_on_voltage(device::PowerSemiconductorSwitch) =
    max(
        device.threshold_v,
        device.forward_voltage_drop_v +
        device.holding_current / device.on_conductance,
    )

function power_semiconductor_reverse_turn_on_voltage(
    device::PowerSemiconductorSwitch,
)
    diode = something(device.antiparallel_diode)
    return diode.forward_voltage_v +
        diode.holding_current_a / diode.on_conductance_s
end

function power_semiconductor_event_localization!(device::PowerSemiconductorSwitch)
    device.event_localization_enabled = true
    return device
end

function _record_power_semiconductor_topology_transition!(
    device::PowerSemiconductorSwitch,
    previous_forward_state::Bool,
    previous_reverse_state::Bool,
    time_s::Real,
)
    if previous_forward_state != device.closed ||
       previous_reverse_state != device.reverse_diode_conducting
        device.topology_transition_count += 1
        device.last_transition_time_s = Float64(time_s)
    end
    return device
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{DiodeJunction},
    _time_s::Float64,
)
    return device
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{ThyristorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{TriacJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = false
    if device.gate_driver.applied_on && !device.closed &&
       abs(device.last_voltage) >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = device.last_voltage < 0.0 ? Int8(-1) : Int8(1)
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on && !device.closed &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
        device.gate_turn_off_disposition = :not_requested
    elseif !device.gate_driver.applied_on && device.closed
        if abs(device.last_forward_current_a) <= device.gate_turn_off_current_limit_a
            _activate_power_semiconductor_tail!(device, time_s)
            device.closed = false
            device.conduction_direction = Int8(0)
            device.gate_turn_off_disposition = :accepted
        elseif device.gate_turn_off_policy === :remain_latched
            device.gate_turn_off_disposition = :remained_latched
        else
            throw(DomainError(
                device.last_forward_current_a,
                "GTO gate turn-off current exceeds its declared interruptible domain",
            ))
        end
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    if device.gate_driver.applied_on &&
       !device.reverse_diode_conducting &&
       device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
    elseif !device.gate_driver.applied_on
        device.closed && _activate_power_semiconductor_tail!(device, time_s)
        device.closed = false
        device.conduction_direction = Int8(0)
    end
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function _apply_power_semiconductor_gate_state!(
    device::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
    time_s::Float64,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.closed = device.gate_driver.applied_on
    device.closed && (device.reverse_diode_conducting = false)
    return _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
end

function request_power_semiconductor_gate!(
    device::PowerSemiconductorSwitch,
    commanded_on::Bool,
    time_s::Real,
    ;
    earliest_transition_time_s::Union{Nothing,Real}=nothing,
)
    driver = device.gate_driver
    driver === nothing && throw(ArgumentError("a diode junction has no gate command"))
    command_time = Float64(time_s)
    isfinite(command_time) && command_time >= 0.0 || throw(ArgumentError(
        "gate command time must be finite and nonnegative",
    ))
    earliest_transition_time = if earliest_transition_time_s === nothing
        nothing
    else
        boundary = Float64(earliest_transition_time_s)
        isfinite(boundary) && boundary >= command_time || throw(ArgumentError(
            "earliest gate-transition time must be finite and no earlier than the command",
        ))
        boundary
    end
    if device isa PowerSemiconductorSwitch{GateTurnOffThyristorJunction} &&
       !commanded_on && device.closed &&
       abs(device.last_forward_current_a) > device.gate_turn_off_current_limit_a &&
       device.gate_turn_off_policy === :refuse
        device.gate_turn_off_disposition = :refused
        throw(DomainError(
            device.last_forward_current_a,
            "GTO gate turn-off current exceeds its declared interruptible domain",
        ))
    end
    commanded_on == driver.commanded_on && return false
    driver.commanded_on = commanded_on
    driver.last_command_time_s = command_time
    driver.command_count += 1
    if commanded_on == driver.applied_on
        driver.pending_state !== nothing && (driver.filtered_pulse_count += 1)
        driver.pending_state = nothing
        driver.pending_transition_time_s = Inf
        return false
    end
    transition_time = if commanded_on
        max(
            command_time + driver.turn_on_delay_s,
            driver.last_turn_off_time_s + driver.dead_time_s,
        )
    else
        max(
            command_time + driver.turn_off_delay_s,
            driver.last_turn_on_time_s + driver.minimum_pulse_width_s,
        )
    end
    if earliest_transition_time !== nothing
        transition_time = max(transition_time, earliest_transition_time)
    end
    driver.pending_state = commanded_on
    driver.pending_transition_time_s = transition_time
    if transition_time <= command_time + 16 * eps(max(1.0, command_time))
        apply_power_semiconductor_gate_transition!(device, command_time)
        return true
    end
    return false
end

function power_semiconductor_gate_transition_time(device::PowerSemiconductorSwitch)
    driver = device.gate_driver
    driver === nothing && return nothing
    driver.pending_state === nothing && return nothing
    return driver.pending_transition_time_s
end

function power_semiconductor_gate_transition_residual(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    transition_time_s = power_semiconductor_gate_transition_time(device)
    transition_time_s === nothing && return nothing
    return transition_time_s - Float64(time_s)
end

function _apply_power_semiconductor_gate_event!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = something(device.extended_fidelity)
    fidelity.pending_event_current_a = abs(device.last_forward_current_a)
    fidelity.pending_event_blocking_voltage_v = abs(device.last_voltage)
    transitioned = apply_power_semiconductor_gate_transition!(device, time_s)
    transitioned || return device
    event_energy_j = _deposit_power_semiconductor_switching_energy!(
        device,
        device.last_voltage,
        device.last_forward_current_a,
        Float64(time_s),
    )
    device.semiconductor_dissipated_energy_j += event_energy_j
    return device
end

function apply_power_semiconductor_gate_transition!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    driver = device.gate_driver
    driver === nothing && throw(ArgumentError("a diode junction has no gate transition"))
    driver.pending_state === nothing && return false
    transition_time = Float64(time_s)
    tolerance = 16 * eps(max(1.0, abs(transition_time)))
    transition_time + tolerance >= driver.pending_transition_time_s || throw(ArgumentError(
        "gate transition was requested before its delay/dead-time/minimum-pulse boundary",
    ))
    applied_on = something(driver.pending_state)
    driver.applied_on = applied_on
    driver.pending_state = nothing
    driver.pending_transition_time_s = Inf
    driver.transition_count += 1
    if applied_on
        driver.last_turn_on_time_s = transition_time
    else
        driver.last_turn_off_time_s = transition_time
    end
    _apply_power_semiconductor_gate_state!(device, transition_time)
    return true
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{DiodeJunction},
)
    device.closed || device.reverse_diode_conducting ? nothing :
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device)
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{ThyristorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    driver = device.gate_driver
    !device.closed && driver.applied_on ?
        abs(device.last_voltage) - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

function power_semiconductor_forward_turn_on_residual(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
)
    driver = device.gate_driver
    !device.closed && !device.reverse_diode_conducting && driver.applied_on ?
        device.last_voltage - power_semiconductor_forward_turn_on_voltage(device) : nothing
end

power_semiconductor_forward_turn_on_residual(
    ::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
) = nothing

function power_semiconductor_forward_extinction_residual(
    device::PowerSemiconductorSwitch,
)
    device.closed ? device.last_forward_current_a - device.holding_current : nothing
end

function power_semiconductor_forward_extinction_residual(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    device.closed ?
        device.conduction_direction * device.last_forward_current_a - device.holding_current :
        nothing
end

power_semiconductor_forward_extinction_residual(
    ::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
) = nothing

function power_semiconductor_reverse_turn_on_residual(
    device::PowerSemiconductorSwitch,
)
    diode = device.antiparallel_diode
    diode === nothing && return nothing
    device.closed || device.reverse_diode_conducting ? nothing :
        -device.last_voltage - power_semiconductor_reverse_turn_on_voltage(device)
end

function power_semiconductor_reverse_extinction_residual(
    device::PowerSemiconductorSwitch,
)
    diode = device.antiparallel_diode
    diode === nothing && return nothing
    device.reverse_diode_conducting ?
        device.last_reverse_diode_current_a - diode.holding_current_a : nothing
end

"""Return the accepted stored-charge exhaustion surface without mutating device state."""
function power_semiconductor_recovery_zero_residual(
    device::PowerSemiconductorSwitch,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    recovery = fidelity.recovered_charge
    recovery === nothing && return nothing
    if recovery.stored_charge_c > 0.0 ||
       recovery.previous_stored_charge_c > 0.0 &&
       recovery.last_recovery_zero_time_s < fidelity.candidate_time_s
        return fidelity.candidate_recovery_charge_c
    end
    return nothing
end

"""Return the within-companion exhaustion boundary after an endpoint recovery probe."""
function power_semiconductor_recovery_zero_time(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    recovery = fidelity.recovered_charge
    recovery === nothing && return nothing
    return nothing
end

"""Accept one localized stored-charge exhaustion transition exactly once."""
function apply_power_semiconductor_recovery_zero!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no recovered-charge transition",
    ))
    recovery = fidelity.recovered_charge
    recovery === nothing && throw(ArgumentError(
        "power semiconductor has no recovered-charge component",
    ))
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("recovery-zero event time must be finite"))
    recovery.last_recovery_zero_time_s < time || throw(ArgumentError(
        "recovery-zero event must advance beyond the last accepted occurrence",
    ))
    localized_charge_c = fidelity.candidate_recovery_charge_c
    localization_scale_c = max(
        abs(recovery.previous_stored_charge_c),
        abs(recovery.last_recovery_current_a) * fidelity.candidate_step_s,
        eps(Float64),
    )
    root_time_tolerance_s = 1.0e-12
    tolerance_c = max(
        8192.0 * eps(Float64) * max(1.0, localization_scale_c),
        abs(recovery.last_recovery_current_a) * root_time_tolerance_s,
    )
    abs(localized_charge_c) <= tolerance_c || throw(DomainError(
        localized_charge_c,
        "recovery-zero transition requires localized exhausted charge",
    ))
    recovery.stored_charge_c = 0.0
    recovery.previous_stored_charge_c = 0.0
    recovery.recovery_active = false
    recovery.last_recovery_current_a = 0.0
    recovery.last_recovery_duration_s = isfinite(recovery.recovery_start_time_s) ?
        max(0.0, time - recovery.recovery_start_time_s) : 0.0
    recovery.recovery_zero_event_count += 1
    recovery.last_recovery_zero_time_s = time
    return device
end

"""Return the analytic tail-current cutoff surface at one candidate time."""
function power_semiconductor_tail_cutoff_residual(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    tail = fidelity.turn_off_tail
    tail === nothing && return nothing
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("tail-cutoff evaluation time must be finite"))
    tail.initial_current_a > tail.cutoff_current_a || return nothing
    isfinite(tail.turn_off_time_s) || return nothing
    tail.last_cutoff_time_s < tail.turn_off_time_s || return nothing
    elapsed_time_s = max(0.0, time - tail.turn_off_time_s)
    return tail.initial_current_a * exp(-elapsed_time_s / tail.decay_time_s) -
        tail.cutoff_current_a
end

"""Return the exact positive-cutoff time for one active exponential tail."""
function power_semiconductor_tail_cutoff_time(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return nothing
    tail = fidelity.turn_off_tail
    tail === nothing && return nothing
    tail.active || return nothing
    tail.cutoff_current_a > 0.0 || return nothing
    tail.initial_current_a > tail.cutoff_current_a || return nothing
    tail.last_cutoff_time_s < tail.turn_off_time_s || return nothing
    return tail.turn_off_time_s + tail.decay_time_s *
        log(tail.initial_current_a / tail.cutoff_current_a)
end

"""Accept one localized exponential-tail cutoff exactly once."""
function apply_power_semiconductor_tail_cutoff!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no tail-cutoff transition",
    ))
    tail = fidelity.turn_off_tail
    tail === nothing && throw(ArgumentError(
        "power semiconductor has no turn-off-tail component",
    ))
    time = Float64(time_s)
    isfinite(time) || throw(ArgumentError("tail-cutoff event time must be finite"))
    tail.last_cutoff_time_s < time || throw(ArgumentError(
        "tail-cutoff event must advance beyond the last accepted occurrence",
    ))
    residual = power_semiconductor_tail_cutoff_residual(device, time)
    residual === nothing && throw(ArgumentError(
        "turn-off tail has no pending localized cutoff",
    ))
    root_time_tolerance_s = 1.0e-12
    tolerance = max(
        8192.0 * eps(Float64) * max(1.0, abs(tail.cutoff_current_a)),
        abs(tail.cutoff_current_a) *
            expm1(root_time_tolerance_s / tail.decay_time_s),
    )
    abs(residual) <= tolerance || throw(DomainError(
        residual,
        "tail-cutoff transition requires its localized current surface",
    ))
    tail.active = false
    tail.current_a = 0.0
    tail.initial_current_a = 0.0
    tail.last_duration_s = max(0.0, time - tail.turn_off_time_s)
    tail.cutoff_event_count += 1
    tail.last_cutoff_time_s = time
    return device
end

const _POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT = NonlinearDeviceEventSurface(
    :stored_charge_exhaustion,
    (device, _time_s) -> power_semiconductor_recovery_zero_residual(device),
    apply_power_semiconductor_recovery_zero!;
    direction=:falling,
    priority=10,
    candidate_time=power_semiconductor_recovery_zero_time,
)

const _POWER_SEMICONDUCTOR_GATE_TRANSITION_EVENT = NonlinearDeviceEventSurface(
    :gate_transition,
    power_semiconductor_gate_transition_residual,
    _apply_power_semiconductor_gate_event!;
    direction=:falling,
    priority=0,
    topology_invalidating=true,
    candidate_time=power_semiconductor_gate_transition_time,
)

const _POWER_SEMICONDUCTOR_TAIL_CUTOFF_EVENT = NonlinearDeviceEventSurface(
    :turn_off_tail_cutoff,
    power_semiconductor_tail_cutoff_residual,
    apply_power_semiconductor_tail_cutoff!;
    direction=:falling,
    priority=20,
    candidate_time=power_semiconductor_tail_cutoff_time,
)

function nonlinear_device_event_surfaces(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && return ()
    gate_event = power_semiconductor_gate_transition_time(device) === nothing ? () :
        (_POWER_SEMICONDUCTOR_GATE_TRANSITION_EVENT,)
    recovery_event = fidelity.recovered_charge === nothing ? () :
        (_POWER_SEMICONDUCTOR_RECOVERY_ZERO_EVENT,)
    tail_event = fidelity.turn_off_tail === nothing ? () :
        (_POWER_SEMICONDUCTOR_TAIL_CUTOFF_EVENT,)
    return (gate_event..., recovery_event..., tail_event...)
end

function apply_power_semiconductor_forward_turn_on!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.reverse_diode_conducting && throw(ArgumentError(
        "forward conduction cannot begin while the antiparallel diode conducts",
    ))
    device.closed = true
    device.conduction_direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
        (device.last_voltage < 0.0 ? Int8(-1) : Int8(1)) : Int8(1)
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_forward_extinction!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    _activate_power_semiconductor_tail!(device, Float64(time_s))
    device.closed = false
    device.conduction_direction = Int8(0)
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_reverse_turn_on!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    device.antiparallel_diode === nothing && throw(ArgumentError(
        "power semiconductor has no antiparallel diode",
    ))
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.closed && throw(ArgumentError(
        "antiparallel conduction cannot begin while the forward path conducts",
    ))
    device.reverse_diode_conducting = true
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function apply_power_semiconductor_reverse_extinction!(
    device::PowerSemiconductorSwitch,
    time_s::Real,
)
    previous_forward = device.closed
    previous_reverse = device.reverse_diode_conducting
    device.reverse_diode_conducting = false
    _record_power_semiconductor_topology_transition!(
        device,
        previous_forward,
        previous_reverse,
        time_s,
    )
    return device
end

function _series_rc_snubber_companion(snubber::SeriesRCSnubber, dt_s::Float64)
    dt_s > 0.0 && isfinite(dt_s) ||
        throw(ArgumentError("snubber timestep must be finite and positive"))
    capacitive_resistance = dt_s / (2.0 * snubber.capacitance_f)
    conductance = inv(snubber.resistance_ohm + capacitive_resistance)
    history_current = conductance * (
        -capacitive_resistance * snubber.previous_current_a -
        snubber.capacitor_voltage_v
    )
    return conductance, history_current
end

function _series_rc_snubber_backward_euler_companion(
    snubber::SeriesRCSnubber,
    dt_s::Float64,
)
    dt_s > 0.0 && isfinite(dt_s) || throw(ArgumentError(
        "snubber timestep must be finite and positive",
    ))
    capacitive_resistance = dt_s / snubber.capacitance_f
    conductance = inv(snubber.resistance_ohm + capacitive_resistance)
    history_current = -conductance * snubber.capacitor_voltage_v
    return conductance, history_current
end

function _stamp_power_semiconductor_gate_due!(
    device::PowerSemiconductorSwitch,
    time_s::Float64,
)
    device.event_localization_enabled && return device
    transition_time = power_semiconductor_gate_transition_time(device)
    transition_time === nothing && return device
    transition_time <= time_s + 16 * eps(max(1.0, abs(time_s))) &&
        apply_power_semiconductor_gate_transition!(device, time_s)
    return device
end

function stamp!(
    admittance::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
)
    _stamp_power_semiconductor_gate_due!(device, time_s)
    device.last_evaluation_time_s = time_s
    conductance = device.off_conductance
    history_current = 0.0
    if device.closed
        conductance = device.on_conductance
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        history_current = -conductance * device.forward_voltage_drop_v * direction
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        conductance = diode.on_conductance_s
        history_current = conductance * diode.forward_voltage_v
    end
    device.last_conductance = conductance
    device.last_history_current_a = history_current
    stamp_conductance!(admittance, device.a, device.b, conductance)
    stamp_history_current!(rhs, device.a, device.b, history_current)
    if device.snubber !== nothing
        snubber_conductance, snubber_history_current =
            _series_rc_snubber_companion(device.snubber, dt_s)
        stamp_conductance!(admittance, device.a, device.b, snubber_conductance)
        stamp_history_current!(rhs, device.a, device.b, snubber_history_current)
    end
    return nothing
end

Branches.backward_euler_companion_supported(::PowerSemiconductorSwitch) = true

function stamp!(
    admittance::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
    ::Val{TrapezoidalCompanion},
)
    return stamp!(admittance, rhs, device, time_s, dt_s)
end

function stamp!(
    admittance::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    dt_s::Float64,
    ::Val{BackwardEulerCompanion},
)
    _stamp_power_semiconductor_gate_due!(device, time_s)
    device.last_evaluation_time_s = time_s
    conductance = device.off_conductance
    history_current = 0.0
    if device.closed
        conductance = device.on_conductance
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        history_current = -conductance * device.forward_voltage_drop_v * direction
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        conductance = diode.on_conductance_s
        history_current = conductance * diode.forward_voltage_v
    end
    device.last_conductance = conductance
    device.last_history_current_a = history_current
    stamp_conductance!(admittance, device.a, device.b, conductance)
    stamp_history_current!(rhs, device.a, device.b, history_current)
    if device.snubber !== nothing
        snubber_conductance, snubber_history_current =
            _series_rc_snubber_backward_euler_companion(device.snubber, dt_s)
        stamp_conductance!(admittance, device.a, device.b, snubber_conductance)
        stamp_history_current!(rhs, device.a, device.b, snubber_history_current)
    end
    return nothing
end

function _update_series_rc_snubber!(
    snubber::SeriesRCSnubber,
    terminal_voltage_v::Float64,
    dt_s::Float64,
)
    previous_current = snubber.previous_current_a
    previous_loss = snubber.last_resistor_loss_w
    conductance, history_current = _series_rc_snubber_companion(snubber, dt_s)
    current = conductance * terminal_voltage_v + history_current
    capacitor_voltage = snubber.capacitor_voltage_v +
        dt_s / (2.0 * snubber.capacitance_f) * (current + previous_current)
    resistor_loss = snubber.resistance_ohm * current^2
    snubber.dissipated_energy_j += 0.5 * dt_s * (previous_loss + resistor_loss)
    snubber.previous_current_a = current
    snubber.capacitor_voltage_v = capacitor_voltage
    snubber.last_branch_voltage_v = terminal_voltage_v
    snubber.last_current_a = current
    snubber.last_resistor_loss_w = resistor_loss
    return snubber
end

function _update_series_rc_snubber_backward_euler!(
    snubber::SeriesRCSnubber,
    terminal_voltage_v::Float64,
    dt_s::Float64,
)
    previous_loss = snubber.last_resistor_loss_w
    conductance, history_current =
        _series_rc_snubber_backward_euler_companion(snubber, dt_s)
    current = conductance * terminal_voltage_v + history_current
    capacitor_voltage = snubber.capacitor_voltage_v +
        dt_s * current / snubber.capacitance_f
    resistor_loss = snubber.resistance_ohm * current^2
    snubber.dissipated_energy_j += 0.5 * dt_s * (previous_loss + resistor_loss)
    snubber.previous_current_a = current
    snubber.capacitor_voltage_v = capacitor_voltage
    snubber.last_branch_voltage_v = terminal_voltage_v
    snubber.last_current_a = current
    snubber.last_resistor_loss_w = resistor_loss
    return snubber
end

function _update_power_semiconductor_sampled_state!(device::DiodeValveSwitch)
    device.closed = diode_next_closed(
        device,
        device.last_voltage,
        device.last_forward_current_a,
    )
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{ThyristorJunction},
)
    if device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false)
    elseif device.gate_driver.applied_on &&
           !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{TriacJunction},
)
    if device.closed
        device.conduction_direction * device.last_forward_current_a < device.holding_current &&
            (device.closed = false; device.conduction_direction = Int8(0))
    elseif device.gate_driver.applied_on &&
           abs(device.last_voltage) >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = device.last_voltage < 0.0 ? Int8(-1) : Int8(1)
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{GateTurnOffThyristorJunction},
)
    if !device.gate_driver.applied_on && device.closed
        if abs(device.last_forward_current_a) <= device.gate_turn_off_current_limit_a
            _activate_power_semiconductor_tail!(device, device.last_evaluation_time_s)
            device.closed = false
            device.conduction_direction = Int8(0)
            device.gate_turn_off_disposition = :accepted
        elseif device.gate_turn_off_policy === :remain_latched
            device.gate_turn_off_disposition = :remained_latched
        else
            throw(DomainError(
                device.last_forward_current_a,
                "GTO gate turn-off current exceeds its declared interruptible domain",
            ))
        end
    elseif device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false; device.conduction_direction = Int8(0))
    elseif device.gate_driver.applied_on && !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
        device.conduction_direction = Int8(1)
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{InsulatedGateBipolarTransistorJunction},
)
    if !device.gate_driver.applied_on
        device.closed = false
    elseif device.closed
        device.last_forward_current_a < device.holding_current &&
            (device.closed = false)
    elseif !device.reverse_diode_conducting &&
           device.last_voltage >= power_semiconductor_forward_turn_on_voltage(device)
        device.closed = true
    end
    return device
end

function _update_power_semiconductor_sampled_state!(
    device::PowerSemiconductorSwitch{MetalOxideSemiconductorFieldEffectTransistorJunction},
)
    device.closed = device.gate_driver.applied_on
    device.closed && (device.reverse_diode_conducting = false)
    return device
end

function _update_power_semiconductor_reverse_state!(device::PowerSemiconductorSwitch)
    diode = device.antiparallel_diode
    diode === nothing && return device
    if device.reverse_diode_conducting
        device.last_reverse_diode_current_a < diode.holding_current_a &&
            (device.reverse_diode_conducting = false)
    elseif !device.closed &&
           -device.last_voltage >= power_semiconductor_reverse_turn_on_voltage(device)
        device.reverse_diode_conducting = true
    end
    return device
end

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
)
    terminal_voltage = Branches.branch_voltage(voltages, device.a, device.b)
    semiconductor_current = device.last_conductance * terminal_voltage +
        device.last_history_current_a
    forward_current = device.closed ? semiconductor_current : 0.0
    reverse_current = device.reverse_diode_conducting ? -semiconductor_current : 0.0
    snubber_current = 0.0
    if device.snubber !== nothing
        _update_series_rc_snubber!(device.snubber, terminal_voltage, dt_s)
        snubber_current = device.snubber.last_current_a
    end
    semiconductor_loss = max(0.0, terminal_voltage * semiconductor_current)
    device.semiconductor_dissipated_energy_j += 0.5 * dt_s * (
        device.previous_semiconductor_loss_w + semiconductor_loss
    )
    device.last_voltage = terminal_voltage
    device.last_forward_current_a = forward_current
    device.last_reverse_diode_current_a = reverse_current
    device.last_snubber_current_a = snubber_current
    device.last_current = semiconductor_current + snubber_current
    device.last_semiconductor_loss_w = semiconductor_loss
    device.previous_semiconductor_loss_w = semiconductor_loss
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        _update_power_semiconductor_sampled_state!(device)
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            device.last_evaluation_time_s,
        )
    end
    return nothing
end


nonlinear_device_formulation(::PowerSemiconductorSwitch) = PhysicalConstitutiveCurrent

function nonlinear_device_provenance(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "only a power semiconductor with explicit extended fidelity may enter the nonlinear-device owner list",
    ))
    return fidelity.provenance
end

nonlinear_terminal_nodes(device::PowerSemiconductorSwitch) = (device.a, device.b)

function _power_semiconductor_candidate_method(method::Symbol)
    if method === :TrapezoidalCompanion || method === :trapezoidal
        return :trapezoidal
    elseif method === :BackwardEulerCompanion || method === :backward_euler
        return :backward_euler
    end
    throw(ArgumentError(
        "extended semiconductor companion method must be trapezoidal or backward Euler",
    ))
end

function prepare_nonlinear_device_step!(
    device::PowerSemiconductorSwitch,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor is linear and cannot be prepared as a nonlinear device",
    ))
    isfinite(time_s) || throw(ArgumentError(
        "extended semiconductor candidate time must be finite",
    ))
    isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
        "extended semiconductor candidate step must be finite and positive",
    ))
    method = _power_semiconductor_candidate_method(companion_method)
    fidelity.candidate_time_s = time_s
    fidelity.candidate_step_s = step_s
    fidelity.candidate_method = method
    fidelity.candidate_prepared = true
    if device.topology_transition_count > fidelity.accepted_topology_transition_count
        fidelity.pending_event_current_a = abs(device.last_forward_current_a)
        fidelity.pending_event_blocking_voltage_v = abs(device.last_voltage)
    end
    return nothing
end

function _power_semiconductor_terminal_voltage(
    terminal_voltage_v::AbstractVector{Float64},
)
    length(terminal_voltage_v) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two terminal voltages",
    ))
    voltage_a_v = terminal_voltage_v[1]
    voltage_b_v = terminal_voltage_v[2]
    isfinite(voltage_a_v) && isfinite(voltage_b_v) || throw(ArgumentError(
        "power-semiconductor terminal voltages must be finite",
    ))
    return voltage_a_v - voltage_b_v
end

function _power_semiconductor_channel_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
)
    if device.closed
        direction = device isa PowerSemiconductorSwitch{TriacJunction} ?
            Float64(device.conduction_direction) : 1.0
        return device.on_conductance * (
            voltage_v - direction * device.forward_voltage_drop_v
        ), device.on_conductance
    elseif device.reverse_diode_conducting
        diode = something(device.antiparallel_diode)
        return diode.on_conductance_s * (
            voltage_v + diode.forward_voltage_v
        ), diode.on_conductance_s
    end
    fidelity = device.extended_fidelity
    recovery = fidelity === nothing ? nothing : fidelity.recovered_charge
    if recovery !== nothing && recovery.stored_charge_c > 0.0 && voltage_v < 0.0
        fidelity.candidate_prepared || throw(ArgumentError(
            "recovered-charge evaluation requires one prepared candidate step",
        ))
        return device.on_conductance * voltage_v, device.on_conductance
    end
    return device.off_conductance * voltage_v, device.off_conductance
end

function _power_semiconductor_dynamic_current_jacobian(
    fidelity::PowerSemiconductorExtendedFidelity,
    voltage_v::Float64,
    time_s::Float64,
)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor must be prepared for an exact candidate step before evaluation",
    ))
    tolerance = 16 * eps(max(1.0, abs(time_s), abs(fidelity.candidate_time_s)))
    abs(time_s - fidelity.candidate_time_s) <= tolerance || throw(ArgumentError(
        "extended semiconductor candidate time does not match the nonlinear evaluation time",
    ))
    step_s = fidelity.candidate_step_s
    current_a = 0.0
    jacobian_s = 0.0
    junction = fidelity.junction_charge
    if junction !== nothing
        charge, capacitance = _power_semiconductor_junction_charge_capacitance(
            junction,
            voltage_v,
        )
        if fidelity.candidate_method === :trapezoidal
            current_a += 2.0 * (charge - junction.previous_charge_c) / step_s -
                junction.last_displacement_current_a
            jacobian_s += 2.0 * capacitance / step_s
        else
            current_a += (charge - junction.previous_charge_c) / step_s
            jacobian_s += capacitance / step_s
        end
    end
    tail = fidelity.turn_off_tail
    if tail !== nothing && tail.active
        duration_s = max(0.0, time_s - tail.turn_off_time_s)
        tail_current = tail.initial_current_a * exp(-duration_s / tail.decay_time_s)
        tail_current > tail.cutoff_current_a && (current_a += tail_current)
    end
    return current_a, jacobian_s
end

function _power_semiconductor_snubber_current_jacobian(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
)
    snubber = device.snubber
    snubber === nothing && return 0.0, 0.0
    fidelity = something(device.extended_fidelity)
    conductance, history_current = if fidelity.candidate_method === :trapezoidal
        _series_rc_snubber_companion(snubber, fidelity.candidate_step_s)
    else
        _series_rc_snubber_backward_euler_companion(snubber, fidelity.candidate_step_s)
    end
    return conductance * voltage_v + history_current, conductance
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two terminal currents",
    ))
    size(terminal_jacobian_s) == (2, 2) || throw(DimensionMismatch(
        "power semiconductor requires a two-by-two terminal Jacobian",
    ))
    voltage_v = _power_semiconductor_terminal_voltage(terminal_voltage_v)
    channel_current_a, channel_jacobian_s =
        _power_semiconductor_channel_current_jacobian(device, voltage_v)
    fidelity = something(device.extended_fidelity)
    dynamic_current_a, dynamic_jacobian_s =
        _power_semiconductor_dynamic_current_jacobian(fidelity, voltage_v, time_s)
    snubber_current_a, snubber_jacobian_s =
        _power_semiconductor_snubber_current_jacobian(device, voltage_v)
    current_a = channel_current_a + dynamic_current_a + snubber_current_a
    jacobian_s = channel_jacobian_s + dynamic_jacobian_s + snubber_jacobian_s
    terminal_current_a[1] = current_a
    terminal_current_a[2] = -current_a
    terminal_jacobian_s[1, 1] = jacobian_s
    terminal_jacobian_s[1, 2] = -jacobian_s
    terminal_jacobian_s[2, 1] = -jacobian_s
    terminal_jacobian_s[2, 2] = jacobian_s
    return nothing
end

function _power_semiconductor_switching_temperature_k(
    fidelity::PowerSemiconductorExtendedFidelity,
)
    thermal = fidelity.thermal
    thermal === nothing && return 293.15
    return thermal.node_temperature_k[1]
end

function _deposit_power_semiconductor_switching_energy!(
    device::PowerSemiconductorSwitch,
    voltage_v::Float64,
    current_a::Float64,
    time_s::Float64,
)
    fidelity = something(device.extended_fidelity)
    energy_table = fidelity.switching_energy
    energy_table === nothing && return 0.0
    event_kind = :none
    recovery = fidelity.recovered_charge
    if recovery !== nothing && recovery.recovery_start_time_s == time_s &&
       energy_table.last_reverse_recovery_start_time_s != time_s
        event_kind = :reverse_recovery
    elseif device.topology_transition_count > energy_table.last_event_transition_count
        event_kind = device.closed ? :turn_on : :turn_off
    end
    event_kind === :none && return 0.0
    interpolation_current_a = event_kind in (:turn_on, :turn_off) ?
        fidelity.pending_event_current_a : abs(current_a)
    interpolation_voltage_v = event_kind in (:turn_on, :turn_off) ?
        fidelity.pending_event_blocking_voltage_v : abs(voltage_v)
    event_energy_j = power_semiconductor_switching_energy(
        energy_table,
        event_kind,
        interpolation_current_a,
        interpolation_voltage_v,
        _power_semiconductor_switching_temperature_k(fidelity),
    )
    if event_kind === :turn_on
        energy_table.cumulative_turn_on_energy_j += event_energy_j
        energy_table.last_event_transition_count = device.topology_transition_count
    elseif event_kind === :turn_off
        energy_table.cumulative_turn_off_energy_j += event_energy_j
        energy_table.last_event_transition_count = device.topology_transition_count
    else
        energy_table.cumulative_reverse_recovery_energy_j += event_energy_j
        energy_table.last_reverse_recovery_start_time_s = time_s
    end
    energy_table.last_event_kind = event_kind
    energy_table.last_event_energy_j = event_energy_j
    fidelity.thermal === nothing ||
        _deposit_cauer_event_energy!(fidelity.thermal, event_energy_j)
    return event_energy_j
end

function _accept_power_semiconductor_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
    terminal_jacobian_s::Union{Nothing,AbstractMatrix{Float64}},
)
    fidelity = something(device.extended_fidelity)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor acceptance requires one prepared candidate step",
    ))
    voltage_v = _power_semiconductor_terminal_voltage(terminal_voltage_v)
    length(terminal_current_a) == 2 || throw(DimensionMismatch(
        "power semiconductor requires two accepted terminal currents",
    ))
    current_a = terminal_current_a[1]
    isfinite(current_a) && isfinite(time_s) || throw(ArgumentError(
        "accepted power-semiconductor current and time must be finite",
    ))
    step_s = fidelity.candidate_step_s
    previous_terminal_voltage_v = fidelity.previous_terminal_voltage_v
    previous_terminal_current_a = fidelity.previous_terminal_current_a
    previous_snubber_current_a = device.snubber === nothing ? 0.0 :
        device.snubber.last_current_a
    previous_snubber_stored_energy_j = device.snubber === nothing ? 0.0 :
        0.5 * device.snubber.capacitance_f * device.snubber.capacitor_voltage_v^2
    previous_snubber_dissipated_energy_j = device.snubber === nothing ? 0.0 :
        device.snubber.dissipated_energy_j
    previous_semiconductor_loss_w = device.previous_semiconductor_loss_w
    channel_current_a, channel_jacobian_s =
        _power_semiconductor_channel_current_jacobian(device, voltage_v)
    tail = fidelity.turn_off_tail
    evaluated_tail_current_a = if tail !== nothing && tail.active
        duration_s = max(0.0, time_s - tail.turn_off_time_s)
        tail.initial_current_a * exp(-duration_s / tail.decay_time_s)
    else
        0.0
    end
    accepted_tail_current_a = tail !== nothing &&
        evaluated_tail_current_a > tail.cutoff_current_a ? evaluated_tail_current_a : 0.0
    snubber_current_a, snubber_jacobian_s =
        _power_semiconductor_snubber_current_jacobian(device, voltage_v)
    junction = fidelity.junction_charge
    if junction !== nothing
        previous_charge_c = junction.previous_charge_c
        previous_voltage_v = junction.previous_voltage_v
        previous_energy_j = power_semiconductor_junction_stored_energy(
            junction,
            previous_voltage_v,
        )
        previous_displacement_current_a = junction.last_displacement_current_a
        displacement_current_a = current_a - channel_current_a -
            accepted_tail_current_a - snubber_current_a
        if terminal_jacobian_s === nothing
            charge, capacitance = _power_semiconductor_junction_charge_capacitance(
                junction,
                voltage_v,
            )
        else
            size(terminal_jacobian_s) == (2, 2) || throw(DimensionMismatch(
                "power semiconductor requires a two-by-two accepted terminal Jacobian",
            ))
            dynamic_jacobian_s = (
                terminal_jacobian_s[1, 1] -
                channel_jacobian_s -
                snubber_jacobian_s
            )
            if fidelity.candidate_method === :trapezoidal
                capacitance = 0.5 * dynamic_jacobian_s * step_s
                charge = previous_charge_c + 0.5 * step_s * (
                    junction.last_displacement_current_a + displacement_current_a
                )
            else
                capacitance = dynamic_jacobian_s * step_s
                charge = previous_charge_c + step_s * displacement_current_a
            end
        end
        junction.previous_voltage_v = voltage_v
        junction.previous_charge_c = charge
        junction.last_capacitance_f = capacitance
        junction.last_charge_c = charge
        junction.last_displacement_current_a = displacement_current_a
        current_energy_j = power_semiconductor_junction_stored_energy(
            junction,
            voltage_v,
        )
        junction_work_j = if fidelity.candidate_method === :trapezoidal
            0.5 * step_s * (
                previous_voltage_v * previous_displacement_current_a +
                voltage_v * displacement_current_a
            )
        else
            step_s * voltage_v * displacement_current_a
        end
        fidelity.companion_energy_residual_j +=
            junction_work_j - (current_energy_j - previous_energy_j)
    end
    recovery = fidelity.recovered_charge
    if recovery !== nothing
        previous_charge_c = recovery.stored_charge_c
        candidate_charge_c =
            (previous_charge_c + step_s * channel_current_a) /
            (1.0 + step_s / recovery.lifetime_s)
        fidelity.candidate_recovery_charge_c = candidate_charge_c
        accepted_charge_c = max(0.0, candidate_charge_c)
        charge_change_c = accepted_charge_c - previous_charge_c
        if channel_current_a >= 0.0
            recovery.last_recovery_current_a = 0.0
            recovery.recovery_active = false
        else
            recovered_charge_c = max(0.0, -charge_change_c)
            recovery.last_recovery_current_a = min(0.0, channel_current_a)
            if recovered_charge_c > 0.0 && !recovery.recovery_active
                recovery.recovery_start_time_s = time_s
            end
            recovery.recovery_active = accepted_charge_c > 0.0 && recovered_charge_c > 0.0
            recovery.peak_reverse_current_a = min(
                recovery.peak_reverse_current_a,
                recovery.last_recovery_current_a,
            )
            recovery.cumulative_recovered_charge_c += recovered_charge_c
            recovery.last_recovery_duration_s = isfinite(recovery.recovery_start_time_s) ?
                max(0.0, time_s - recovery.recovery_start_time_s) : 0.0
        end
        recovery.previous_stored_charge_c = previous_charge_c
        recovery.stored_charge_c = max(0.0, accepted_charge_c)
    end
    if tail !== nothing && tail.active
        duration_s = max(0.0, time_s - tail.turn_off_time_s)
        if evaluated_tail_current_a <= tail.cutoff_current_a
            tail.active = false
            tail.current_a = 0.0
        else
            tail.current_a = evaluated_tail_current_a
        end
        tail.last_duration_s = duration_s
    end
    tail_current_a = tail === nothing ? 0.0 : tail.current_a
    instantaneous_loss_w = max(
        0.0,
        voltage_v * (channel_current_a + tail_current_a),
    )
    thermal = fidelity.thermal
    device.last_voltage = voltage_v
    device.last_current = current_a
    device.last_forward_current_a = device.closed ? channel_current_a : 0.0
    device.last_reverse_diode_current_a =
        device.reverse_diode_conducting ? -channel_current_a : 0.0
    if device.snubber !== nothing
        if fidelity.candidate_method === :trapezoidal
            _update_series_rc_snubber!(device.snubber, voltage_v, step_s)
        else
            _update_series_rc_snubber_backward_euler!(device.snubber, voltage_v, step_s)
        end
        device.last_snubber_current_a = device.snubber.last_current_a
    end
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        _update_power_semiconductor_sampled_state!(device)
        if previous_forward && !device.closed
            _activate_power_semiconductor_tail!(device, time_s)
        end
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            time_s,
        )
    end
    event_energy_j = _deposit_power_semiconductor_switching_energy!(
        device,
        voltage_v,
        channel_current_a,
        time_s,
    )
    thermal === nothing || _accept_cauer_thermal_step!(
        thermal,
        instantaneous_loss_w,
        step_s,
    )
    device.semiconductor_dissipated_energy_j += 0.5 * step_s * (
        previous_semiconductor_loss_w + instantaneous_loss_w
    ) + event_energy_j
    previous_junction_current_a = junction === nothing ? 0.0 :
        previous_displacement_current_a
    channel_work_j = 0.5 * step_s * (
        previous_terminal_voltage_v * (
            previous_terminal_current_a - previous_junction_current_a -
            previous_snubber_current_a
        ) +
        voltage_v * (channel_current_a + tail_current_a)
    )
    semiconductor_energy_increment_j = 0.5 * step_s * (
        previous_semiconductor_loss_w + instantaneous_loss_w
    ) + event_energy_j
    fidelity.companion_energy_residual_j +=
        channel_work_j - semiconductor_energy_increment_j
    if device.snubber !== nothing
        snubber = device.snubber
        snubber_work_j = 0.5 * step_s * (
            previous_terminal_voltage_v * previous_snubber_current_a +
            voltage_v * snubber.last_current_a
        )
        snubber_stored_energy_j =
            0.5 * snubber.capacitance_f * snubber.capacitor_voltage_v^2
        fidelity.companion_energy_residual_j += snubber_work_j -
            (snubber.dissipated_energy_j - previous_snubber_dissipated_energy_j) -
            (snubber_stored_energy_j - previous_snubber_stored_energy_j)
    end
    device.last_semiconductor_loss_w = instantaneous_loss_w
    device.previous_semiconductor_loss_w = instantaneous_loss_w
    fidelity.previous_terminal_voltage_v = voltage_v
    fidelity.previous_terminal_current_a = current_a
    fidelity.accepted_topology_transition_count = device.topology_transition_count
    return nothing
end


function accept_nonlinear_device_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    time_s::Float64,
)
    return _accept_power_semiconductor_state!(
        device,
        terminal_voltage_v,
        terminal_current_a,
        time_s,
        nothing,
    )
end

function accept_nonlinear_device_state!(
    device::PowerSemiconductorSwitch,
    terminal_voltage_v::AbstractVector{Float64},
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    return _accept_power_semiconductor_state!(
        device,
        terminal_voltage_v,
        terminal_current_a,
        time_s,
        terminal_jacobian_s,
    )
end

function finish_nonlinear_device_step!(device::PowerSemiconductorSwitch)
    fidelity = something(device.extended_fidelity)
    fidelity.candidate_prepared || throw(ArgumentError(
        "extended semiconductor candidate step was already finished or never prepared",
    ))
    fidelity.candidate_prepared = false
    return nothing
end

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
    ::Val{TrapezoidalCompanion},
)
    return update!(device, voltages, dt_s)
end

function update!(
    device::PowerSemiconductorSwitch,
    voltages::AbstractVector{Float64},
    dt_s::Float64,
    ::Val{BackwardEulerCompanion},
)
    terminal_voltage = Branches.branch_voltage(voltages, device.a, device.b)
    semiconductor_current = device.last_conductance * terminal_voltage +
        device.last_history_current_a
    forward_current = device.closed ? semiconductor_current : 0.0
    reverse_current = device.reverse_diode_conducting ? -semiconductor_current : 0.0
    snubber_current = 0.0
    if device.snubber !== nothing
        _update_series_rc_snubber_backward_euler!(device.snubber, terminal_voltage, dt_s)
        snubber_current = device.snubber.last_current_a
    end
    semiconductor_loss = max(0.0, terminal_voltage * semiconductor_current)
    device.semiconductor_dissipated_energy_j += 0.5 * dt_s * (
        device.previous_semiconductor_loss_w + semiconductor_loss
    )
    device.last_voltage = terminal_voltage
    device.last_forward_current_a = forward_current
    device.last_reverse_diode_current_a = reverse_current
    device.last_snubber_current_a = snubber_current
    device.last_current = semiconductor_current + snubber_current
    device.last_semiconductor_loss_w = semiconductor_loss
    device.previous_semiconductor_loss_w = semiconductor_loss
    if !device.event_localization_enabled
        previous_forward = device.closed
        previous_reverse = device.reverse_diode_conducting
        _update_power_semiconductor_sampled_state!(device)
        _update_power_semiconductor_reverse_state!(device)
        _record_power_semiconductor_topology_transition!(
            device,
            previous_forward,
            previous_reverse,
            device.last_evaluation_time_s,
        )
    end
    return nothing
end

function power_semiconductor_terminal_state(device::PowerSemiconductorSwitch)
    driver = device.gate_driver
    snubber = device.snubber
    return PowerSemiconductorTerminalState(
        power_semiconductor_kind(device),
        device.closed,
        device.reverse_diode_conducting,
        driver === nothing ? false : driver.commanded_on,
        driver === nothing ? false : driver.applied_on,
        device.last_voltage,
        device.last_current,
        device.last_forward_current_a,
        device.last_reverse_diode_current_a,
        device.last_snubber_current_a,
        device.last_semiconductor_loss_w,
        snubber === nothing ? 0.0 : snubber.last_resistor_loss_w,
        device.semiconductor_dissipated_energy_j,
        snubber === nothing ? 0.0 : snubber.dissipated_energy_j,
        snubber === nothing ? 0.0 :
            0.5 * snubber.capacitance_f * snubber.capacitor_voltage_v^2,
    )
end

power_semiconductor_has_extended_fidelity(device::PowerSemiconductorSwitch) =
    device.extended_fidelity !== nothing

function initialize_power_semiconductor_junction_state!(
    device::PowerSemiconductorSwitch,
    voltage_v::Real,
)
    fidelity = something(device.extended_fidelity)
    junction = something(fidelity.junction_charge)
    voltage = Float64(voltage_v)
    charge, capacitance = _power_semiconductor_junction_charge_capacitance(
        junction,
        voltage,
    )
    junction.previous_voltage_v = voltage
    junction.previous_charge_c = charge
    junction.last_capacitance_f = capacitance
    junction.last_charge_c = charge
    junction.last_displacement_current_a = 0.0
    fidelity.companion_energy_residual_j = 0.0
    fidelity.previous_terminal_voltage_v = voltage
    return device
end

function _activate_power_semiconductor_tail!(
    device::PowerSemiconductorSwitch,
    time_s::Float64,
)
    fidelity = device.extended_fidelity
    fidelity === nothing && return device
    tail = fidelity.turn_off_tail
    tail === nothing && return device
    initial_current = max(0.0, device.last_forward_current_a)
    if initial_current > tail.cutoff_current_a
        tail.active = true
        tail.current_a = initial_current
        tail.initial_current_a = initial_current
        tail.turn_off_time_s = time_s
        tail.last_duration_s = 0.0
    else
        tail.active = false
        tail.current_a = 0.0
        tail.initial_current_a = 0.0
        tail.turn_off_time_s = time_s
        tail.last_duration_s = 0.0
    end
    return device
end

function _power_semiconductor_fidelity_components(
    fidelity::PowerSemiconductorExtendedFidelity,
)
    components = Symbol[]
    fidelity.recovered_charge === nothing || push!(components, :recovered_charge)
    fidelity.junction_charge === nothing || push!(components, :nonlinear_junction_charge)
    fidelity.turn_off_tail === nothing || push!(components, :turn_off_tail)
    fidelity.switching_energy === nothing || push!(components, :event_energy_map)
    fidelity.thermal === nothing || push!(components, :lumped_electrothermal)
    fidelity.declared_model === nothing || push!(components, :declared_compact_model)
    return Tuple(components)
end

function power_semiconductor_extended_state(device::PowerSemiconductorSwitch)
    fidelity = device.extended_fidelity
    fidelity === nothing && throw(ArgumentError(
        "baseline power semiconductor has no extended state",
    ))
    recovery = fidelity.recovered_charge
    junction = fidelity.junction_charge
    tail = fidelity.turn_off_tail
    energy = fidelity.switching_energy
    thermal = fidelity.thermal
    return PowerSemiconductorExtendedState(
        1,
        _power_semiconductor_fidelity_components(fidelity),
        device.conduction_direction,
        device.gate_turn_off_disposition,
        recovery === nothing ? false : recovery.recovery_active,
        recovery === nothing ? 0.0 : recovery.stored_charge_c,
        recovery === nothing ? 0.0 : recovery.last_recovery_current_a,
        recovery === nothing ? 0.0 : recovery.peak_reverse_current_a,
        recovery === nothing ? 0.0 : recovery.cumulative_recovered_charge_c,
        recovery === nothing ? 0 : recovery.recovery_zero_event_count,
        recovery === nothing ? -Inf : recovery.last_recovery_zero_time_s,
        junction === nothing ? 0.0 : junction.last_capacitance_f,
        junction === nothing ? 0.0 : junction.last_charge_c,
        junction === nothing ? 0.0 : junction.last_displacement_current_a,
        junction === nothing ? 0.0 :
            power_semiconductor_junction_stored_energy(junction, junction.last_charge_c == junction.previous_charge_c ? junction.previous_voltage_v : device.last_voltage),
        fidelity.companion_energy_residual_j,
        tail === nothing ? false : tail.active,
        tail === nothing ? 0.0 : tail.current_a,
        tail === nothing ? 0 : tail.cutoff_event_count,
        tail === nothing ? -Inf : tail.last_cutoff_time_s,
        energy === nothing ? :none : energy.last_event_kind,
        energy === nothing ? 0.0 : energy.last_event_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_turn_on_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_turn_off_energy_j,
        energy === nothing ? 0.0 : energy.cumulative_reverse_recovery_energy_j,
        thermal === nothing ? 0.0 : thermal.node_temperature_k[1],
        thermal === nothing ? Float64[] : copy(thermal.node_temperature_k),
        thermal === nothing ? 0.0 : thermal.last_ambient_heat_flow_w,
        thermal === nothing ? 0.0 : power_semiconductor_thermal_stored_energy(thermal),
    )
end

trace_output_channel_count(::PowerSemiconductorSwitch) = 10
trace_output_is_public(::PowerSemiconductorSwitch) = true

function trace_output_channel_names!(
    names::Vector{Symbol},
    element_name::Symbol,
    ::PowerSemiconductorSwitch,
)
    append!(
        names,
        Symbol[
            Symbol(element_name, :_terminal_voltage_v),
            Symbol(element_name, :_terminal_current_a),
            Symbol(element_name, :_forward_current_a),
            Symbol(element_name, :_reverse_diode_current_a),
            Symbol(element_name, :_gate_commanded),
            Symbol(element_name, :_gate_applied),
            Symbol(element_name, :_semiconductor_loss_w),
            Symbol(element_name, :_snubber_resistor_loss_w),
            Symbol(element_name, :_semiconductor_dissipated_energy_j),
            Symbol(element_name, :_snubber_dissipated_energy_j),
        ],
    )
    return names
end

function trace_output_values!(
    output::AbstractMatrix{Float64},
    first_channel::Int,
    sample::Int,
    device::PowerSemiconductorSwitch,
    _voltage::AbstractVector{Float64},
)
    terminal = power_semiconductor_terminal_state(device)
    values = (
        terminal.terminal_voltage_v,
        terminal.terminal_current_a,
        terminal.forward_current_a,
        terminal.reverse_diode_current_a,
        terminal.gate_commanded_on ? 1.0 : 0.0,
        terminal.gate_applied_on ? 1.0 : 0.0,
        terminal.semiconductor_loss_w,
        terminal.snubber_resistor_loss_w,
        terminal.semiconductor_dissipated_energy_j,
        terminal.snubber_dissipated_energy_j,
    )
    for offset in eachindex(values)
        output[first_channel + offset - 1, sample] = values[offset]
    end
    return first_channel + length(values)
end

@enum BreakerPolePosition::UInt8 begin
    BreakerPoleOpen = 0x01
    BreakerPoleClosed = 0x02
    BreakerPoleOpening = 0x03
    BreakerPoleClosing = 0x04
    BreakerPoleAwaitingCurrentZero = 0x05
end

struct EMTBreakerSpecification
    id::Symbol
    pole_ids::NTuple{3,Symbol}
    closed_conductance_s::Float64
    open_conductance_s::Float64
    opening_travel_ticks::Int
    closing_travel_ticks::Int
    current_zero_required::Bool
    current_zero_threshold_a::Float64
    failure_delay_ticks::Int
    failure_current_threshold_a::Float64
    reclose_dead_ticks::Int
    reclaim_ticks::Int
    maximum_reclose_shots::Int
    contact_tail_enabled::Bool
    physical_provenance::ParameterProvenance
    timing_provenance::ParameterProvenance
    deterministic_signature_sha256::String

    function EMTBreakerSpecification(
        id::Symbol;
        pole_ids::NTuple{3,Symbol}=(:phase_a, :phase_b, :phase_c),
        closed_conductance_s::Real,
        open_conductance_s::Real,
        opening_travel_ticks::Integer,
        closing_travel_ticks::Integer,
        current_zero_required::Bool,
        current_zero_threshold_a::Real,
        failure_delay_ticks::Integer,
        failure_current_threshold_a::Real,
        reclose_dead_ticks::Integer,
        reclaim_ticks::Integer,
        maximum_reclose_shots::Integer,
        contact_tail_enabled::Bool=false,
        physical_provenance::ParameterProvenance,
        timing_provenance::ParameterProvenance,
    )
        !isempty(String(id)) && all(name -> !isempty(String(name)), pole_ids) &&
            allunique(pole_ids) || throw(ArgumentError(
            "EMT breaker and pole identities must be nonempty and unique",
        ))
        closed_conductance, open_conductance, current_zero_threshold,
            failure_threshold = Float64.((
            closed_conductance_s,
            open_conductance_s,
            current_zero_threshold_a,
            failure_current_threshold_a,
        ))
        all(isfinite, (
            closed_conductance,
            open_conductance,
            current_zero_threshold,
            failure_threshold,
        )) && closed_conductance > open_conductance >= 0.0 &&
            current_zero_threshold >= 0.0 && failure_threshold >= 0.0 ||
            throw(ArgumentError("EMT breaker conductance and current thresholds are invalid"))
        opening_travel = Int(opening_travel_ticks)
        closing_travel = Int(closing_travel_ticks)
        failure_delay = Int(failure_delay_ticks)
        dead_time = Int(reclose_dead_ticks)
        reclaim = Int(reclaim_ticks)
        shots = Int(maximum_reclose_shots)
        opening_travel >= 0 && closing_travel >= 0 && failure_delay >= 0 &&
            dead_time >= 0 && reclaim >= 0 && 0 <= shots <= 16 || throw(ArgumentError(
            "EMT breaker timing and shot settings are outside their domain",
        ))
        physical_provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "EMT breaker conductance provenance must be physical",
        ))
        timing_provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
            "EMT breaker timing provenance must be numerical policy",
        ))
        io = IOBuffer()
        print(
            io,
            id,
            '|',
            join(string.(pole_ids), ','),
            '|',
            repr((
                closed_conductance,
                open_conductance,
                opening_travel,
                closing_travel,
                current_zero_required,
                current_zero_threshold,
                failure_delay,
                failure_threshold,
                dead_time,
                reclaim,
                shots,
                contact_tail_enabled,
            )),
            '|',
            repr(physical_provenance),
            '|',
            repr(timing_provenance),
        )
        return new(
            id,
            pole_ids,
            closed_conductance,
            open_conductance,
            opening_travel,
            closing_travel,
            current_zero_required,
            current_zero_threshold,
            failure_delay,
            failure_threshold,
            dead_time,
            reclaim,
            shots,
            contact_tail_enabled,
            physical_provenance,
            timing_provenance,
            bytes2hex(sha256(take!(io))),
        )
    end
end

mutable struct EMTBreakerPoleRuntime
    position::BreakerPolePosition
    transition_due_tick::Int
    previous_contact_power_w::Float64
    contact_energy_j::Float64
    opening_count::Int
    closing_count::Int
end

struct EMTBreakerEvent
    tick::Int
    kind::Symbol
    pole::Union{Nothing,Symbol}
    shot_count::Int
end

mutable struct EMTBreakerRuntime
    poles::Vector{EMTBreakerPoleRuntime}
    tick_s::Float64
    last_advance_tick::Int
    last_command_tick::Int
    trip_initiation_tick::Int
    failure_deadline_tick::Int
    reclose_due_tick::Int
    reclaim_due_tick::Int
    shot_count::Int
    failure_asserted::Bool
    backup_trip_required::Bool
    lockout::Bool
    events::Vector{EMTBreakerEvent}

    function EMTBreakerRuntime(
        specification::EMTBreakerSpecification;
        tick_s::Real,
        initially_closed::Bool=true,
    )
        tick = Float64(tick_s)
        isfinite(tick) && tick > 0.0 || throw(ArgumentError(
            "EMT breaker tick must be finite and positive",
        ))
        position = initially_closed ? BreakerPoleClosed : BreakerPoleOpen
        poles = [EMTBreakerPoleRuntime(position, -1, 0.0, 0.0, 0, 0) for _ in 1:3]
        return new(poles, tick, -1, -1, -1, -1, -1, -1, 0, false, false, false, EMTBreakerEvent[])
    end
end

struct EMTBreakerAdvanceResult
    tick::Int
    pole_positions::NTuple{3,BreakerPolePosition}
    failure_asserted::Bool
    backup_trip_required::Bool
    reclose_available::Bool
    shot_count::Int
    lockout::Bool
    event_count::Int
    contact_energy_j::NTuple{3,Float64}
end

struct EMTBreakerSnapshot
    specification_signature_sha256::String
    tick_s::Float64
    poles::Vector{EMTBreakerPoleRuntime}
    last_advance_tick::Int
    last_command_tick::Int
    trip_initiation_tick::Int
    failure_deadline_tick::Int
    reclose_due_tick::Int
    reclaim_due_tick::Int
    shot_count::Int
    failure_asserted::Bool
    backup_trip_required::Bool
    lockout::Bool
    events::Vector{EMTBreakerEvent}
end

struct EMTBreakerNetworkBinding
    specification_signature_sha256::String
    terminal_nodes::NTuple{3,NTuple{2,Int}}
    control_refs::NTuple{3,Base.RefValue{Float64}}
    network_elements::Vector{Any}

    function EMTBreakerNetworkBinding(
        specification_signature_sha256::AbstractString,
        terminal_nodes::NTuple{3,NTuple{2,Int}},
        control_refs::NTuple{3,Base.RefValue{Float64}},
        network_elements::AbstractVector,
    )
        signature = lowercase(String(specification_signature_sha256))
        occursin(r"^[0-9a-f]{64}$", signature) || throw(ArgumentError(
            "EMT breaker binding signature must be a 64-hex SHA-256",
        ))
        all(nodes -> nodes[1] >= 0 && nodes[2] >= 0 && nodes[1] != nodes[2], terminal_nodes) ||
            throw(ArgumentError("EMT breaker pole terminals must be distinct nonnegative nodes"))
        length(network_elements) == 3 || throw(ArgumentError(
            "EMT breaker binding requires exactly three physical network elements",
        ))
        all(reference -> isfinite(reference[]), control_refs) || throw(ArgumentError(
            "EMT breaker binding controls must be finite",
        ))
        return new(signature, terminal_nodes, control_refs, Any[network_elements...])
    end
end

function synchronize_emt_breaker_poles!(
    binding::EMTBreakerNetworkBinding,
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
)
    binding.specification_signature_sha256 ==
        specification.deterministic_signature_sha256 || throw(ArgumentError(
        "EMT breaker network binding specification identity is stale",
    ))
    conductances = breaker_pole_conductances(runtime, specification)
    for index in 1:3
        binding.control_refs[index][] =
            conductances[index] == specification.closed_conductance_s ? 1.0 : -1.0
    end
    return binding
end

function _validate_breaker_command_tick(runtime::EMTBreakerRuntime, command_tick::Integer)
    tick = Int(command_tick)
    tick >= 0 && tick >= runtime.last_command_tick && tick >= runtime.last_advance_tick ||
        throw(ArgumentError("EMT breaker command tick is noncausal"))
    runtime.last_command_tick = tick
    return tick
end

function request_breaker_trip!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    command_tick::Integer,
)
    tick = _validate_breaker_command_tick(runtime, command_tick)
    runtime.trip_initiation_tick = tick
    runtime.failure_deadline_tick = tick + specification.failure_delay_ticks
    runtime.failure_asserted = false
    runtime.backup_trip_required = false
    runtime.reclose_due_tick = -1
    runtime.reclaim_due_tick = -1
    if runtime.shot_count >= specification.maximum_reclose_shots && runtime.shot_count > 0
        runtime.lockout = true
    end
    transitioned = false
    for (index, pole) in pairs(runtime.poles)
        if pole.position in (BreakerPoleClosed, BreakerPoleClosing)
            pole.position = BreakerPoleOpening
            pole.transition_due_tick = tick + specification.opening_travel_ticks
            pole.opening_count += 1
            push!(runtime.events, EMTBreakerEvent(
                tick,
                :trip_requested,
                specification.pole_ids[index],
                runtime.shot_count,
            ))
            transitioned = true
        end
    end
    return transitioned
end

function request_breaker_reclose!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    command_tick::Integer,
)
    tick = _validate_breaker_command_tick(runtime, command_tick)
    all(pole -> pole.position === BreakerPoleOpen, runtime.poles) || throw(ArgumentError(
        "EMT breaker reclose requires every pole open",
    ))
    !runtime.failure_asserted && !runtime.lockout || throw(ArgumentError(
        "EMT breaker reclose is unavailable during failure or lockout",
    ))
    runtime.shot_count < specification.maximum_reclose_shots || begin
        runtime.lockout = true
        throw(ArgumentError("EMT breaker reclose shot limit is exhausted"))
    end
    runtime.reclose_due_tick < 0 || throw(ArgumentError(
        "EMT breaker already has a pending reclose",
    ))
    runtime.shot_count += 1
    runtime.reclose_due_tick = tick + specification.reclose_dead_ticks
    push!(runtime.events, EMTBreakerEvent(tick, :reclose_scheduled, nothing, runtime.shot_count))
    return runtime.reclose_due_tick
end

function breaker_pole_conductances(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
)
    return ntuple(3) do index
        runtime.poles[index].position in (
            BreakerPoleClosed,
            BreakerPoleOpening,
            BreakerPoleAwaitingCurrentZero,
        ) ? specification.closed_conductance_s : specification.open_conductance_s
    end
end

function _open_breaker_pole!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    index::Int,
    tick::Int,
)
    pole = runtime.poles[index]
    pole.position = BreakerPoleOpen
    pole.transition_due_tick = -1
    push!(runtime.events, EMTBreakerEvent(
        tick,
        :pole_opened,
        specification.pole_ids[index],
        runtime.shot_count,
    ))
    return nothing
end

function _close_breaker_pole!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    index::Int,
    tick::Int,
)
    pole = runtime.poles[index]
    pole.position = BreakerPoleClosed
    pole.transition_due_tick = -1
    push!(runtime.events, EMTBreakerEvent(
        tick,
        :pole_closed,
        specification.pole_ids[index],
        runtime.shot_count,
    ))
    return nothing
end

function advance_emt_breaker!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    accepted_tick::Integer,
    pole_voltage_v::NTuple{3,<:Real},
    pole_current_a::NTuple{3,<:Real};
    current_zero_reached::NTuple{3,Bool}=(false, false, false),
)
    tick = Int(accepted_tick)
    tick >= 0 && tick > runtime.last_advance_tick || throw(ArgumentError(
        "EMT breaker accepted ticks must advance strictly",
    ))
    voltages = Float64.(pole_voltage_v)
    currents = Float64.(pole_current_a)
    all(isfinite, voltages) && all(isfinite, currents) || throw(ArgumentError(
        "EMT breaker accepted pole voltage and current must be finite",
    ))
    if runtime.last_advance_tick >= 0
        elapsed_s = (tick - runtime.last_advance_tick) * runtime.tick_s
        for index in 1:3
            power = voltages[index] * currents[index]
            pole = runtime.poles[index]
            pole.contact_energy_j +=
                0.5 * (pole.previous_contact_power_w + power) * elapsed_s
            pole.previous_contact_power_w = power
        end
    else
        for index in 1:3
            runtime.poles[index].previous_contact_power_w =
                voltages[index] * currents[index]
        end
    end

    if runtime.reclose_due_tick >= 0 && tick >= runtime.reclose_due_tick
        runtime.lockout && throw(ArgumentError(
            "EMT breaker cannot execute a pending reclose after lockout",
        ))
        for (index, pole) in pairs(runtime.poles)
            pole.position === BreakerPoleOpen || throw(ArgumentError(
                "EMT breaker reclose found a pole that is not open",
            ))
            pole.position = BreakerPoleClosing
            pole.transition_due_tick = tick + specification.closing_travel_ticks
            pole.closing_count += 1
            push!(runtime.events, EMTBreakerEvent(
                tick,
                :close_requested,
                specification.pole_ids[index],
                runtime.shot_count,
            ))
        end
        runtime.reclose_due_tick = -1
    end

    for (index, pole) in pairs(runtime.poles)
        if pole.position === BreakerPoleOpening && tick >= pole.transition_due_tick
            if specification.current_zero_required &&
               !current_zero_reached[index] &&
               abs(currents[index]) > specification.current_zero_threshold_a
                pole.position = BreakerPoleAwaitingCurrentZero
                pole.transition_due_tick = -1
            else
                _open_breaker_pole!(runtime, specification, index, tick)
            end
        elseif pole.position === BreakerPoleAwaitingCurrentZero &&
               (current_zero_reached[index] ||
                abs(currents[index]) <= specification.current_zero_threshold_a)
            _open_breaker_pole!(runtime, specification, index, tick)
        elseif pole.position === BreakerPoleClosing && tick >= pole.transition_due_tick
            _close_breaker_pole!(runtime, specification, index, tick)
        end
    end

    if runtime.failure_deadline_tick >= 0 && tick >= runtime.failure_deadline_tick
        runtime.failure_asserted = any(index ->
            runtime.poles[index].position !== BreakerPoleOpen ||
            abs(currents[index]) > specification.failure_current_threshold_a,
            1:3,
        )
        runtime.backup_trip_required = runtime.failure_asserted
        push!(runtime.events, EMTBreakerEvent(
            tick,
            runtime.failure_asserted ? :breaker_failure : :breaker_trip_success,
            nothing,
            runtime.shot_count,
        ))
        runtime.failure_deadline_tick = -1
    end

    if all(pole -> pole.position === BreakerPoleClosed, runtime.poles) &&
       runtime.shot_count > 0 && runtime.reclaim_due_tick < 0
        runtime.reclaim_due_tick = tick + specification.reclaim_ticks
    end
    if runtime.reclaim_due_tick >= 0 && tick >= runtime.reclaim_due_tick
        runtime.shot_count = 0
        runtime.reclaim_due_tick = -1
        push!(runtime.events, EMTBreakerEvent(tick, :reclaim_complete, nothing, 0))
    end
    runtime.last_advance_tick = tick
    return EMTBreakerAdvanceResult(
        tick,
        Tuple(pole.position for pole in runtime.poles),
        runtime.failure_asserted,
        runtime.backup_trip_required,
        all(pole -> pole.position === BreakerPoleOpen, runtime.poles) &&
            !runtime.failure_asserted && !runtime.lockout &&
            runtime.shot_count < specification.maximum_reclose_shots,
        runtime.shot_count,
        runtime.lockout,
        length(runtime.events),
        Tuple(pole.contact_energy_j for pole in runtime.poles),
    )
end

function emt_breaker_snapshot(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
)
    return EMTBreakerSnapshot(
        specification.deterministic_signature_sha256,
        runtime.tick_s,
        deepcopy(runtime.poles),
        runtime.last_advance_tick,
        runtime.last_command_tick,
        runtime.trip_initiation_tick,
        runtime.failure_deadline_tick,
        runtime.reclose_due_tick,
        runtime.reclaim_due_tick,
        runtime.shot_count,
        runtime.failure_asserted,
        runtime.backup_trip_required,
        runtime.lockout,
        deepcopy(runtime.events),
    )
end

function restore_emt_breaker_snapshot!(
    runtime::EMTBreakerRuntime,
    specification::EMTBreakerSpecification,
    snapshot::EMTBreakerSnapshot,
)
    snapshot.specification_signature_sha256 ==
        specification.deterministic_signature_sha256 || throw(ArgumentError(
        "EMT breaker snapshot specification identity is stale",
    ))
    snapshot.tick_s == runtime.tick_s || throw(ArgumentError(
        "EMT breaker snapshot tick differs from its runtime",
    ))
    length(snapshot.poles) == 3 && all(pole ->
        isfinite(pole.previous_contact_power_w) && isfinite(pole.contact_energy_j) &&
        pole.transition_due_tick >= -1 && pole.opening_count >= 0 && pole.closing_count >= 0,
        snapshot.poles,
    ) || throw(ArgumentError("EMT breaker snapshot pole state is invalid"))
    snapshot.last_advance_tick >= -1 && snapshot.last_command_tick >= -1 &&
        snapshot.shot_count >= 0 || throw(ArgumentError(
        "EMT breaker snapshot calendar or counter is invalid",
    ))
    runtime.poles = deepcopy(snapshot.poles)
    runtime.last_advance_tick = snapshot.last_advance_tick
    runtime.last_command_tick = snapshot.last_command_tick
    runtime.trip_initiation_tick = snapshot.trip_initiation_tick
    runtime.failure_deadline_tick = snapshot.failure_deadline_tick
    runtime.reclose_due_tick = snapshot.reclose_due_tick
    runtime.reclaim_due_tick = snapshot.reclaim_due_tick
    runtime.shot_count = snapshot.shot_count
    runtime.failure_asserted = snapshot.failure_asserted
    runtime.backup_trip_required = snapshot.backup_trip_required
    runtime.lockout = snapshot.lockout
    runtime.events = deepcopy(snapshot.events)
    return runtime
end

module EMTPartitioning

using SHA
using ..EMTTaskPlatform: EMTLogicalTime
import ..MeasurementChains
import ..ModernMachines
import ..TransformerApparatus

export CausalLinearReconstruction,
       CausalZeroOrderReconstruction,
       DirectCoupledExchange,
       EMTDeckPartitionCheckpoint,
       EMTDeckPartitionResult,
       EMTDeckRegion,
       EMTInterfacePort,
       EMTInterfacePortKind,
       EMTPartitionCheckpoint,
       EMTPartitionExecutionMode,
       EMTPartitionExchangeMethod,
       EMTPartitionExchangePolicy,
       EMTPartitionFailure,
       EMTPartitionInterpolation,
       EMTPartitionPlan,
       EMTPartitionRegion,
       EMTPartitionResult,
       EMTPartitionTolerancePolicy,
       EMTRegionalBridgeLegOwner,
       EMTRegionalInstrumentTransformerOwner,
       EMTRegionalModernMachineOwner,
       EMTRegionalSemlyenLineOwner,
       IteratedWaveformExchange,
       LaggedCausalExchange,
       LocalSubcyclingExecution,
       MonolithicReferenceExecution,
       NortonInterfacePort,
       PartitionedLaggedExecution,
       PartitionedWaveformExecution,
       CombinedLocalPartitionedExecution,
       PartitionedDeckEMTStudy,
       PassiveTwoRegionRLCStudy,
       ScatteringInterfacePort,
       TheveninInterfacePort,
       TravelingWaveInterfacePort,
       VoltageCurrentInterfacePort,
       emt_interface_coordinates,
       emt_interface_physical_values,
       emt_partition_plan,
       partition_checkpoint_signature_sha256,
       partition_plan_signature_sha256,
       partition_result,
       partition_study_signature_sha256

const _PARTITION_IDENTITY = r"^[A-Za-z][A-Za-z0-9_.:-]*$"
const _MAXIMUM_REGIONS = 64
const _MAXIMUM_LOCAL_RATES = 16
const _MAXIMUM_RATE_RATIO = 10_000

@enum EMTPartitionExchangeMethod::UInt8 begin
    DirectCoupledExchange = 0x01
    IteratedWaveformExchange = 0x02
    LaggedCausalExchange = 0x03
end

@enum EMTPartitionExecutionMode::UInt8 begin
    MonolithicReferenceExecution = 0x01
    LocalSubcyclingExecution = 0x02
    PartitionedLaggedExecution = 0x03
    PartitionedWaveformExecution = 0x04
    CombinedLocalPartitionedExecution = 0x05
end

"""One public Semlyen frequency-dependent line bound to named regional terminals."""
struct EMTRegionalSemlyenLineOwner{M<:Tuple}
    identity::String
    from_terminal_names::Tuple
    to_terminal_names::Tuple
    modes::M
    voltage_modal_to_phase::Matrix{ComplexF64}
    current_modal_to_phase::Matrix{ComplexF64}
    local_step::Union{Nothing,EMTLogicalTime}
    signature_sha256::String

    function EMTRegionalSemlyenLineOwner(
        identity::AbstractString,
        from_terminal_names,
        to_terminal_names,
        modes;
        voltage_modal_to_phase,
        current_modal_to_phase,
        local_step::Union{Nothing,EMTLogicalTime}=nothing,
    )
        lines = _partition_lines_module()
        name = _partition_identity(
            identity,
            "regional Semlyen line identity",
        )
        from_terminals = Tuple(
            _partition_terminal_identity(
                String(terminal),
                "regional Semlyen sending terminal",
            ) for terminal in from_terminal_names
        )
        to_terminals = Tuple(
            _partition_terminal_identity(
                String(terminal),
                "regional Semlyen receiving terminal",
            ) for terminal in to_terminal_names
        )
        mode_values = Tuple(modes)
        mode_type = getfield(lines, :SemlyenModeParameters)
        !isempty(mode_values) && all(mode -> mode isa mode_type, mode_values) ||
            throw(ArgumentError(
                "regional Semlyen line requires public mode parameters",
            ))
        length(from_terminals) == length(to_terminals) == length(mode_values) ||
            throw(DimensionMismatch(
                "regional Semlyen terminal and mode counts must match",
            ))
        all(zip(from_terminals, to_terminals)) do (from_terminal, to_terminal)
            from_terminal != to_terminal
        end || throw(ArgumentError(
            "regional Semlyen sending and receiving terminals must differ",
        ))
        voltage_transform = Matrix{ComplexF64}(voltage_modal_to_phase)
        current_transform = Matrix{ComplexF64}(current_modal_to_phase)
        Base.invokelatest(
            getfield(lines, :semlyen_line_physical_checks),
            collect(mode_values),
            voltage_transform,
            current_transform,
        ).physical_checks_passed || throw(ArgumentError(
            "regional Semlyen line parameters fail their physical checks",
        ))
        local_step !== nothing && local_step <= zero(local_step) && throw(
            ArgumentError("regional Semlyen local step must be positive"),
        )
        signature_lines = String[
            "schema=aimora.emt.regional_semlyen_line.v1",
            "identity=$name",
            "from=$(join(from_terminals, ','))",
            "to=$(join(to_terminals, ','))",
            "voltage_transform=$(join(bitstring.(reinterpret(Float64, vec(voltage_transform))), ','))",
            "current_transform=$(join(bitstring.(reinterpret(Float64, vec(current_transform))), ','))",
            "local_step=" * (local_step === nothing ? "network" :
                "$(local_step.numerator)/$(local_step.denominator)"),
        ]
        for (index, mode) in enumerate(mode_values)
            push!(
                signature_lines,
                "mode=$index|$(repr(mode.characteristic_admittance_s))|" *
                "$(repr(mode.travel_time_s))|$(repr(mode.phasor_series_impedance))|" *
                "$(repr(mode.phasor_characteristic_admittance))|" *
                "$(repr(mode.phasor_frequency_hz))",
            )
            for (kind, terms) in (
                (:propagation, mode.propagation_terms),
                (:admittance, mode.admittance_terms),
            )
                for (term_index, term) in enumerate(terms)
                    push!(
                        signature_lines,
                        "term=$index|$kind|$term_index|$(repr(term.pole))|" *
                        "$(repr(term.residue))|$(term.conjugate_pair)",
                    )
                end
            end
        end
        return new{typeof(mode_values)}(
            name,
            from_terminals,
            to_terminals,
            mode_values,
            voltage_transform,
            current_transform,
            local_step,
            _partition_runtime_signature(signature_lines),
        )
    end
end

@enum EMTPartitionInterpolation::UInt8 begin
    CausalZeroOrderReconstruction = 0x01
    CausalLinearReconstruction = 0x02
end

@enum EMTInterfacePortKind::UInt8 begin
    VoltageCurrentInterfacePort = 0x01
    NortonInterfacePort = 0x02
    TheveninInterfacePort = 0x03
    ScatteringInterfacePort = 0x04
    TravelingWaveInterfacePort = 0x05
end

function _partition_identity(value::AbstractString, owner::AbstractString)
    identity = String(value)
    occursin(_PARTITION_IDENTITY, identity) || throw(ArgumentError(
        "$owner must be portable nonempty text beginning with a letter",
    ))
    return identity
end

function _partition_terminal_identity(value::AbstractString, owner::AbstractString)
    identity = String(value)
    identity == "0" && return identity
    return _partition_identity(identity, owner)
end

function _partition_runtime_signature(lines)
    return bytes2hex(sha256(join(String.(lines), '\n')))
end

function _partition_lines_module()
    root = parentmodule(@__MODULE__)
    isdefined(root, :Lines) || throw(ArgumentError(
        "regional Semlyen line construction requires an activated line runtime",
    ))
    return getfield(root, :Lines)
end

"""One prepared modern-machine owner bound to named nodes inside a regional deck."""
struct EMTRegionalModernMachineOwner{P,E<:Tuple}
    identity::String
    terminal_names::NTuple{4,String}
    preparation::P
    events::E
    signature_sha256::String

    function EMTRegionalModernMachineOwner(
        identity::AbstractString,
        terminal_names,
        preparation::P;
        events=(),
    ) where {P}
        preparation isa ModernMachines.ModernMachinePreparation || throw(
            ArgumentError(
                "regional modern-machine owner requires a prepared modern machine",
            ),
        )
        terminals = Tuple(String.(terminal_names))
        length(terminals) == 4 || throw(DimensionMismatch(
            "regional modern-machine owner requires three phases and one neutral",
        ))
        normalized_terminals = ntuple(4) do index
            _partition_terminal_identity(
                terminals[index],
                "regional modern-machine terminal",
            )
        end
        all(terminal -> terminal != "0", normalized_terminals[1:3]) || throw(
            ArgumentError("regional modern-machine phase terminals cannot be ground"),
        )
        length(unique(normalized_terminals)) == 4 || throw(ArgumentError(
            "regional modern-machine terminals must be distinct",
        ))
        ordered_events = Tuple(sort!(
            ModernMachines.ModernMachineEvent[events...];
            by=event -> (event.time_s, event.priority, String(event.id)),
        ))
        length(unique(getfield.(ordered_events, :id))) == length(ordered_events) ||
            throw(ArgumentError("regional modern-machine event ids must be unique"))
        name = _partition_identity(identity, "regional modern-machine identity")
        signature_lines = String[
            "schema=aimora.emt.regional_modern_machine.v1",
            "identity=$name",
            "terminals=$(join(normalized_terminals, ','))",
            "preparation=$(preparation.deterministic_signature_sha256)",
        ]
        for event in ordered_events
            push!(
                signature_lines,
                "event=$(event.id)|$(repr(event.time_s))|$(UInt8(event.kind))|" *
                "$(join(repr.(event.values), ','))|$(event.enabled_value)|$(event.priority)",
            )
        end
        return new{P,typeof(ordered_events)}(
            name,
            normalized_terminals,
            preparation,
            ordered_events,
            _partition_runtime_signature(signature_lines),
        )
    end
end

"""One physical transformer and sampled measurement chain bound to named regional nodes."""
struct EMTRegionalInstrumentTransformerOwner{D,P,M}
    identity::String
    terminal_names::Tuple
    definition::D
    preparation::P
    measurement_specification::M
    signature_sha256::String

    function EMTRegionalInstrumentTransformerOwner(
        identity::AbstractString,
        terminal_names,
        definition::D,
        preparation::P,
        measurement_specification::M,
    ) where {D,P,M}
        definition isa MeasurementChains.InstrumentTransformerMeasurementDefinition ||
            throw(ArgumentError(
                "regional instrument owner requires a physical transformer measurement definition",
            ))
        preparation isa TransformerApparatus.TransformerApparatusPreparation ||
            throw(ArgumentError(
                "regional instrument owner requires a prepared transformer apparatus",
            ))
        measurement_specification isa MeasurementChains.MeasurementChainSpecification ||
            throw(ArgumentError(
                "regional instrument owner requires a sampled measurement specification",
            ))
        terminals = Tuple(
            _partition_terminal_identity(
                String(terminal),
                "regional instrument-transformer terminal",
            ) for terminal in terminal_names
        )
        isempty(terminals) && throw(ArgumentError(
            "regional instrument-transformer terminals must not be empty",
        ))
        length(terminals) == length(unique(terminals)) || throw(ArgumentError(
            "regional instrument-transformer terminals must be distinct",
        ))
        name = _partition_identity(identity, "regional instrument-transformer identity")
        preparation.specification.deterministic_signature_sha256 ==
            definition.apparatus.deterministic_signature_sha256 || throw(
            ArgumentError(
                "regional instrument preparation does not match its physical definition",
            ),
        )
        signature = _partition_runtime_signature((
            "schema=aimora.emt.regional_instrument_transformer.v1",
            "identity=$name",
            "terminals=$(join(terminals, ','))",
            "definition=$(definition.deterministic_signature_sha256)",
            "preparation=$(preparation.preparation_signature_sha256)",
            "measurement=$(MeasurementChains.measurement_chain_signature(measurement_specification))",
        ))
        return new{D,P,M}(
            name,
            terminals,
            definition,
            preparation,
            measurement_specification,
            signature,
        )
    end
end

"""One switch-detailed complementary bridge leg and its exact trailing-edge PWM task."""
struct EMTRegionalBridgeLegOwner
    identity::String
    dc_positive_terminal::String
    ac_terminal::String
    dc_negative_terminal::String
    initially_upper_on::Bool
    on_conductance_s::Float64
    off_conductance_s::Float64
    forward_voltage_drop_v::Float64
    diode_forward_voltage_v::Float64
    diode_conductance_s::Float64
    snubber_resistance_ohm::Float64
    snubber_capacitance_f::Float64
    gate_turn_on_delay_s::Float64
    gate_turn_off_delay_s::Float64
    commutation_dead_time_s::Float64
    minimum_pulse_width_s::Float64
    pwm_task_identity::String
    pwm_tick::EMTLogicalTime
    pwm_period::EMTLogicalTime
    pwm_first_time::EMTLogicalTime
    pwm_duty::Float64
    pwm_priority::Int
    signature_sha256::String

    function EMTRegionalBridgeLegOwner(
        identity::AbstractString,
        dc_positive_terminal::AbstractString,
        ac_terminal::AbstractString,
        dc_negative_terminal::AbstractString;
        initially_upper_on::Bool=false,
        on_conductance_s::Real=20.0,
        off_conductance_s::Real=0.0,
        forward_voltage_drop_v::Real=0.8,
        diode_forward_voltage_v::Real=0.6,
        diode_conductance_s::Real=20.0,
        snubber_resistance_ohm::Real=5.0,
        snubber_capacitance_f::Real=2.0e-6,
        gate_turn_on_delay_s::Real=0.0,
        gate_turn_off_delay_s::Real=0.0,
        commutation_dead_time_s::Real=0.0,
        minimum_pulse_width_s::Real=0.0,
        pwm_task_identity::AbstractString="converter_pwm",
        pwm_tick::EMTLogicalTime,
        pwm_period::EMTLogicalTime,
        pwm_first_time::EMTLogicalTime=zero(pwm_tick),
        pwm_duty::Real=0.5,
        pwm_priority::Integer=0,
    )
        name = _partition_identity(identity, "regional bridge-leg identity")
        terminals = (
            _partition_terminal_identity(
                dc_positive_terminal,
                "regional bridge-leg positive DC terminal",
            ),
            _partition_terminal_identity(
                ac_terminal,
                "regional bridge-leg AC terminal",
            ),
            _partition_terminal_identity(
                dc_negative_terminal,
                "regional bridge-leg negative DC terminal",
            ),
        )
        length(unique(terminals)) == 3 || throw(ArgumentError(
            "regional bridge-leg terminals must be distinct",
        ))
        task_identity = _partition_identity(
            pwm_task_identity,
            "regional bridge-leg PWM task identity",
        )
        pwm_tick > zero(pwm_tick) || throw(ArgumentError(
            "regional bridge-leg PWM tick must be positive",
        ))
        pwm_period > zero(pwm_period) || throw(ArgumentError(
            "regional bridge-leg PWM period must be positive",
        ))
        pwm_first_time >= zero(pwm_first_time) || throw(ArgumentError(
            "regional bridge-leg PWM first time must be nonnegative",
        ))
        _logical_ratio(pwm_period, pwm_tick) === nothing && throw(ArgumentError(
            "regional bridge-leg PWM period must be an integer number of ticks",
        ))
        _logical_ratio(pwm_first_time, pwm_tick) === nothing && throw(ArgumentError(
            "regional bridge-leg PWM first time must be an integer number of ticks",
        ))
        numeric_values = Float64.((
            on_conductance_s,
            off_conductance_s,
            forward_voltage_drop_v,
            diode_forward_voltage_v,
            diode_conductance_s,
            snubber_resistance_ohm,
            snubber_capacitance_f,
            gate_turn_on_delay_s,
            gate_turn_off_delay_s,
            commutation_dead_time_s,
            minimum_pulse_width_s,
            pwm_duty,
        ))
        all(isfinite, numeric_values) || throw(ArgumentError(
            "regional bridge-leg parameters must be finite",
        ))
        numeric_values[1] > 0.0 || throw(ArgumentError(
            "regional bridge-leg on conductance must be positive",
        ))
        numeric_values[2] >= 0.0 || throw(ArgumentError(
            "regional bridge-leg off conductance must be nonnegative",
        ))
        all(>=(0.0), numeric_values[3:11]) || throw(ArgumentError(
            "regional bridge-leg drops, diode, snubber, and timing parameters must be nonnegative",
        ))
        0.0 <= numeric_values[12] <= 1.0 || throw(ArgumentError(
            "regional bridge-leg PWM duty must be within [0, 1]",
        ))
        signature_lines = String[
            "schema=aimora.emt.regional_bridge_leg.v1",
            "identity=$name",
            "terminals=$(join(terminals, ','))",
            "initially_upper_on=$initially_upper_on",
            "pwm_task=$task_identity",
            "pwm_tick=$(pwm_tick.numerator)/$(pwm_tick.denominator)",
            "pwm_period=$(pwm_period.numerator)/$(pwm_period.denominator)",
            "pwm_first=$(pwm_first_time.numerator)/$(pwm_first_time.denominator)",
            "pwm_priority=$(Int(pwm_priority))",
        ]
        for (field, value) in zip(
            (
                :on_conductance_s,
                :off_conductance_s,
                :forward_voltage_drop_v,
                :diode_forward_voltage_v,
                :diode_conductance_s,
                :snubber_resistance_ohm,
                :snubber_capacitance_f,
                :gate_turn_on_delay_s,
                :gate_turn_off_delay_s,
                :commutation_dead_time_s,
                :minimum_pulse_width_s,
                :pwm_duty,
            ),
            numeric_values,
        )
            push!(signature_lines, "$(field)=$(repr(value))")
        end
        return new(
            name,
            terminals...,
            initially_upper_on,
            numeric_values[1:11]...,
            task_identity,
            pwm_tick,
            pwm_period,
            pwm_first_time,
            numeric_values[12],
            Int(pwm_priority),
            _partition_runtime_signature(signature_lines),
        )
    end
end

function _positive_finite(value::Real, owner::AbstractString)
    normalized = Float64(value)
    isfinite(normalized) && normalized > 0.0 || throw(ArgumentError(
        "$owner must be finite and positive",
    ))
    return normalized
end

function _nonnegative_finite(value::Real, owner::AbstractString)
    normalized = Float64(value)
    isfinite(normalized) && normalized >= 0.0 || throw(ArgumentError(
        "$owner must be finite and nonnegative",
    ))
    return normalized
end

function _logical_ratio(numerator::EMTLogicalTime, denominator::EMTLogicalTime)
    denominator > zero(denominator) || throw(ArgumentError(
        "partition logical-time denominator must be positive",
    ))
    wide_numerator = BigInt(numerator.numerator) * denominator.denominator
    wide_denominator = BigInt(numerator.denominator) * denominator.numerator
    wide_denominator > 0 || throw(ArgumentError(
        "partition logical-time ratio must be positive",
    ))
    iszero(rem(wide_numerator, wide_denominator)) || return nothing
    quotient = div(wide_numerator, wide_denominator)
    typemin(Int) <= quotient <= typemax(Int) || throw(OverflowError(
        "partition logical-time ratio exceeds Int",
    ))
    return Int(quotient)
end

"""One exact fixed-step regional owner; each physical model identity may occur in only one region."""
struct EMTPartitionRegion
    identity::String
    model_identities::Tuple
    local_step::EMTLogicalTime

    function EMTPartitionRegion(
        identity::AbstractString,
        model_identities,
        local_step::EMTLogicalTime,
    )
        name = _partition_identity(identity, "partition region identity")
        local_step > zero(local_step) || throw(ArgumentError(
            "partition region local step must be positive",
        ))
        models = Tuple(sort!(String[
            _partition_identity(String(model), "partition model identity")
            for model in model_identities
        ]))
        isempty(models) && throw(ArgumentError(
            "partition region must own at least one model",
        ))
        length(models) == length(unique(models)) || throw(ArgumentError(
            "partition region repeats a model identity",
        ))
        return new(name, models, local_step)
    end
end

"""One explicitly oriented interface; both regional currents are positive outward toward the interface."""
struct EMTInterfacePort
    identity::String
    kind::EMTInterfacePortKind
    positive_region::String
    negative_region::String
    positive_terminal::String
    negative_terminal::String
    voltage_unit::String
    current_unit::String
    voltage_base_v::Float64
    current_base_a::Float64
    reference_impedance_ohm::Float64

    function EMTInterfacePort(
        identity::AbstractString,
        kind::EMTInterfacePortKind,
        positive_region::AbstractString,
        negative_region::AbstractString,
        positive_terminal::AbstractString,
        negative_terminal::AbstractString;
        voltage_unit::AbstractString = "V",
        current_unit::AbstractString = "A",
        voltage_base_v::Real = 1.0,
        current_base_a::Real = 1.0,
        reference_impedance_ohm::Real = 1.0,
    )
        positive = _partition_identity(positive_region, "positive port region")
        negative = _partition_identity(negative_region, "negative port region")
        positive != negative || throw(ArgumentError(
            "partition interface must join two distinct regions",
        ))
        String(voltage_unit) == "V" || throw(ArgumentError(
            "partition voltage/current ports currently require SI volts",
        ))
        String(current_unit) == "A" || throw(ArgumentError(
            "partition voltage/current ports currently require SI amperes",
        ))
        return new(
            _partition_identity(identity, "partition interface identity"),
            kind,
            positive,
            negative,
            _partition_identity(positive_terminal, "positive port terminal"),
            _partition_identity(negative_terminal, "negative port terminal"),
            String(voltage_unit),
            String(current_unit),
            _positive_finite(voltage_base_v, "partition voltage base"),
            _positive_finite(current_base_a, "partition current base"),
            _positive_finite(
                reference_impedance_ohm,
                "partition reference impedance",
            ),
        )
    end
end

"""Map physical oriented voltage/current into the declared equivalent-port coordinates."""
function emt_interface_coordinates(
    port::EMTInterfacePort,
    voltage_v::Real,
    outward_current_a::Real,
)
    voltage = Float64(voltage_v)
    current = Float64(outward_current_a)
    all(isfinite, (voltage, current)) || throw(ArgumentError(
        "partition interface values must be finite",
    ))
    impedance = port.reference_impedance_ohm
    first_coordinate, second_coordinate = if port.kind == VoltageCurrentInterfacePort
        voltage, current
    elseif port.kind == NortonInterfacePort
        voltage, current + voltage / impedance
    elseif port.kind == TheveninInterfacePort
        voltage + impedance * current, current
    elseif port.kind == ScatteringInterfacePort
        scale = 2.0 * sqrt(impedance)
        (voltage + impedance * current) / scale,
            (voltage - impedance * current) / scale
    elseif port.kind == TravelingWaveInterfacePort
        0.5 * (voltage + impedance * current),
            0.5 * (voltage - impedance * current)
    else
        throw(ArgumentError("partition interface kind is unsupported"))
    end
    return (
        kind=port.kind,
        first=first_coordinate,
        second=second_coordinate,
        voltage_unit=port.voltage_unit,
        current_unit=port.current_unit,
        reference_impedance_ohm=impedance,
    )
end

"""Recover the common physical voltage and outward current from equivalent-port coordinates."""
function emt_interface_physical_values(
    port::EMTInterfacePort,
    first_coordinate::Real,
    second_coordinate::Real,
)
    first_value = Float64(first_coordinate)
    second_value = Float64(second_coordinate)
    all(isfinite, (first_value, second_value)) || throw(ArgumentError(
        "partition interface coordinates must be finite",
    ))
    impedance = port.reference_impedance_ohm
    voltage, current = if port.kind == VoltageCurrentInterfacePort
        first_value, second_value
    elseif port.kind == NortonInterfacePort
        first_value, second_value - first_value / impedance
    elseif port.kind == TheveninInterfacePort
        first_value - impedance * second_value, second_value
    elseif port.kind == ScatteringInterfacePort
        scale = sqrt(impedance)
        scale * (first_value + second_value),
            (first_value - second_value) / scale
    elseif port.kind == TravelingWaveInterfacePort
        first_value + second_value,
            (first_value - second_value) / impedance
    else
        throw(ArgumentError("partition interface kind is unsupported"))
    end
    return (voltage_v=voltage, outward_current_a=current)
end

"""Independent physical and scaled acceptance budgets for one interface window."""
struct EMTPartitionTolerancePolicy
    voltage_absolute_v::Float64
    voltage_relative::Float64
    current_absolute_a::Float64
    current_relative::Float64
    kcl_absolute_a::Float64
    kcl_relative::Float64
    interface_energy_absolute_j::Float64
    interface_energy_relative::Float64
    communication_error_absolute_v::Float64
    communication_error_relative::Float64

    function EMTPartitionTolerancePolicy(;
        voltage_absolute_v::Real = 1.0e-8,
        voltage_relative::Real = 1.0e-7,
        current_absolute_a::Real = 1.0e-9,
        current_relative::Real = 1.0e-7,
        kcl_absolute_a::Real = 1.0e-9,
        kcl_relative::Real = 1.0e-8,
        interface_energy_absolute_j::Real = 1.0e-9,
        interface_energy_relative::Real = 1.0e-6,
        communication_error_absolute_v::Real = 1.0e-8,
        communication_error_relative::Real = 2.0e-2,
    )
        values = Float64.(tuple(
            voltage_absolute_v,
            voltage_relative,
            current_absolute_a,
            current_relative,
            kcl_absolute_a,
            kcl_relative,
            interface_energy_absolute_j,
            interface_energy_relative,
            communication_error_absolute_v,
            communication_error_relative,
        ))
        all(isfinite, values) && all(>=(0.0), values) || throw(ArgumentError(
            "partition tolerances must be finite and nonnegative",
        ))
        values[1] > 0.0 && values[3] > 0.0 && values[5] > 0.0 &&
            values[7] > 0.0 && values[9] > 0.0 || throw(ArgumentError(
                "partition absolute tolerances must be positive",
            ))
        return new(values...)
    end
end

"""A declared deterministic exchange algorithm and its bounded convergence policy."""
struct EMTPartitionExchangePolicy
    method::EMTPartitionExchangeMethod
    interpolation::EMTPartitionInterpolation
    maximum_iterations::Int
    relaxation::Float64
    tolerances::EMTPartitionTolerancePolicy

    function EMTPartitionExchangePolicy(
        method::EMTPartitionExchangeMethod = IteratedWaveformExchange;
        interpolation::EMTPartitionInterpolation = CausalLinearReconstruction,
        maximum_iterations::Integer = 32,
        relaxation::Real = 1.0,
        tolerances::EMTPartitionTolerancePolicy = EMTPartitionTolerancePolicy(),
    )
        iteration_count = Int(maximum_iterations)
        1 <= iteration_count <= 10_000 || throw(ArgumentError(
            "partition maximum iterations must be from 1 through 10,000",
        ))
        damping = Float64(relaxation)
        isfinite(damping) && 0.0 < damping <= 1.0 || throw(ArgumentError(
            "partition relaxation must be in (0, 1]",
        ))
        return new(method, interpolation, iteration_count, damping, tolerances)
    end
end

struct EMTPartitionPlan{R<:Tuple,P<:Tuple}
    execution_mode::EMTPartitionExecutionMode
    start::EMTLogicalTime
    stop::EMTLogicalTime
    communication_step::EMTLogicalTime
    regions::R
    ports::P
    exchange::EMTPartitionExchangePolicy
    rate_ratios::Tuple
    communication_window_count::Int
    signature_sha256::String
end

function _partition_signature_lines(
    execution_mode,
    start,
    stop,
    communication_step,
    regions,
    ports,
    exchange,
    ratios,
)
    lines = String[
        "schema=aimora.emt.partition_plan.v2",
        "execution_mode=$(UInt8(execution_mode))",
        "start=$(start.numerator)/$(start.denominator)",
        "stop=$(stop.numerator)/$(stop.denominator)",
        "communication_step=$(communication_step.numerator)/$(communication_step.denominator)",
        "method=$(UInt8(exchange.method))",
        "interpolation=$(UInt8(exchange.interpolation))",
        "maximum_iterations=$(exchange.maximum_iterations)",
        "relaxation=$(repr(exchange.relaxation))",
    ]
    tolerance = exchange.tolerances
    for field in fieldnames(EMTPartitionTolerancePolicy)
        push!(lines, "tolerance.$field=$(repr(getfield(tolerance, field)))")
    end
    for (region, ratio) in zip(regions, ratios)
        push!(lines, "region=$(region.identity)|$(join(region.model_identities, ','))|$(region.local_step.numerator)/$(region.local_step.denominator)|$ratio")
    end
    for port in ports
        push!(lines, "port=$(port.identity)|$(UInt8(port.kind))|$(port.positive_region)|$(port.negative_region)|$(port.positive_terminal)|$(port.negative_terminal)|$(port.voltage_unit)|$(port.current_unit)|$(repr(port.voltage_base_v))|$(repr(port.current_base_a))|$(repr(port.reference_impedance_ohm))")
    end
    return lines
end

"""Solver-free declaration of one canonical EMT deck owned by one partition region."""
struct EMTDeckRegion
    identity::String
    deck_lines::Tuple
    source_identity::String
    initial_voltage_source::Symbol
    saturated_transformer_runtime::Bool
    coupled_lumped_history_runtime::Bool
    distributed_line_runtime::Bool
    runtime_owners::Tuple
    signature_sha256::String

    function EMTDeckRegion(
        identity::AbstractString,
        deck_lines;
        source_identity::AbstractString = identity,
        initial_voltage_source::Symbol = :none,
        saturated_transformer_runtime::Bool = false,
        coupled_lumped_history_runtime::Bool = false,
        distributed_line_runtime::Bool = true,
        runtime_owners=(),
    )
        region_identity = _partition_identity(identity, "deck region identity")
        source = _partition_identity(
            source_identity,
            "deck region source identity",
        )
        initial_voltage_source in (:none, :steady_state) || throw(ArgumentError(
            "deck region initial voltage source must be :none or :steady_state",
        ))
        lines = Tuple(String(line) for line in deck_lines)
        isempty(lines) && throw(ArgumentError(
            "deck region must contain at least one deck line",
        ))
        all(line -> !occursin('\0', line), lines) || throw(ArgumentError(
            "deck region lines must not contain NUL bytes",
        ))
        owners = Tuple(runtime_owners)
        all(owners) do owner
            owner isa EMTRegionalModernMachineOwner ||
                owner isa EMTRegionalInstrumentTransformerOwner ||
                owner isa EMTRegionalSemlyenLineOwner ||
                owner isa EMTRegionalBridgeLegOwner
        end || throw(ArgumentError(
            "deck region runtime owners must use a supported typed regional declaration",
        ))
        owner_identities = getfield.(owners, :identity)
        length(owner_identities) == length(unique(owner_identities)) || throw(
            ArgumentError("deck region repeats a runtime-owner identity"),
        )
        signature_lines = String[
            "schema=aimora.emt.deck_region.v1",
            "identity=$region_identity",
            "source_identity=$source",
            "initial_voltage_source=$initial_voltage_source",
            "saturated_transformer_runtime=$saturated_transformer_runtime",
            "coupled_lumped_history_runtime=$coupled_lumped_history_runtime",
            "distributed_line_runtime=$distributed_line_runtime",
        ]
        for (index, line) in enumerate(lines)
            push!(
                signature_lines,
                "line=$index|$(ncodeunits(line))|$(bytes2hex(sha256(line)))",
            )
        end
        for owner in owners
            push!(
                signature_lines,
                "runtime_owner=$(owner.identity)|$(owner.signature_sha256)",
            )
        end
        return new(
            region_identity,
            lines,
            source,
            initial_voltage_source,
            saturated_transformer_runtime,
            coupled_lumped_history_runtime,
            distributed_line_runtime,
            owners,
            bytes2hex(sha256(join(signature_lines, '\n'))),
        )
    end
end

"""A partition plan bound to complete canonical regional decks and initial port currents."""
struct PartitionedDeckEMTStudy{P<:EMTPartitionPlan,R<:Tuple}
    plan::P
    regions::R
    initial_interface_current_a::Tuple
    signature_sha256::String

    function PartitionedDeckEMTStudy(
        plan::P,
        regions;
        initial_interface_current_a = zeros(length(plan.ports)),
    ) where {P<:EMTPartitionPlan}
        declarations = Tuple(regions)
        all(region -> region isa EMTDeckRegion, declarations) || throw(
            ArgumentError("partitioned deck regions must use EMTDeckRegion"),
        )
        declared_names = getfield.(declarations, :identity)
        length(declared_names) == length(unique(declared_names)) || throw(
            ArgumentError("partitioned deck study repeats a region identity"),
        )
        plan_names = getfield.(plan.regions, :identity)
        Set(declared_names) == Set(plan_names) || throw(ArgumentError(
            "partitioned deck declarations must match every plan region exactly",
        ))
        ordered = Tuple(
            declarations[only(findall(==(name), declared_names))]
            for name in plan_names
        )
        currents = Tuple(Float64.(initial_interface_current_a))
        length(currents) == length(plan.ports) || throw(DimensionMismatch(
            "partitioned deck initial currents must contain one value per interface",
        ))
        all(isfinite, currents) || throw(ArgumentError(
            "partitioned deck initial interface currents must be finite",
        ))
        signature_lines = String[
            "schema=aimora.emt.partitioned_deck_study.v1",
            "plan=$(plan.signature_sha256)",
        ]
        for region in ordered
            push!(signature_lines, "region=$(region.identity)|$(region.signature_sha256)")
        end
        for (port, current) in zip(plan.ports, currents)
            push!(signature_lines, "initial_current=$(port.identity)|$(repr(current))")
        end
        return new{P,typeof(ordered)}(
            plan,
            ordered,
            currents,
            bytes2hex(sha256(join(signature_lines, '\n'))),
        )
    end
end

partition_study_signature_sha256(study::PartitionedDeckEMTStudy) =
    study.signature_sha256

"""Accepted synchronization samples and conservative interface diagnostics."""
struct EMTDeckPartitionResult
    execution_mode::EMTPartitionExecutionMode
    plan_signature_sha256::String
    study_signature_sha256::String
    region_identities::Tuple
    port_identities::Tuple
    time_s::Vector{Float64}
    positive_terminal_voltage_v::Matrix{Float64}
    negative_terminal_voltage_v::Matrix{Float64}
    initialization_voltage_residual_v::Vector{Float64}
    initialization_consistent::Bool
    interface_current_a::Matrix{Float64}
    voltage_residual_v::Matrix{Float64}
    kcl_residual_a::Matrix{Float64}
    interface_energy_defect_j::Matrix{Float64}
    interface_value_age_s::Matrix{Float64}
    communication_error_estimate_v::Matrix{Float64}
    topology_revalidated::Vector{Bool}
    fixed_point_iterations::Vector{Int}
    regional_local_step_counts::Tuple
    accepted_window_count::Int
    rejected_window_count::Int
    refinement_count::Int
    fallback_count::Int
    last_failure::Union{Nothing,String}
    accepted::Bool
    deterministic_signature_sha256::String
end

function emt_partition_plan(
    regions,
    ports;
    start::EMTLogicalTime,
    stop::EMTLogicalTime,
    communication_step::EMTLogicalTime,
    exchange::EMTPartitionExchangePolicy = EMTPartitionExchangePolicy(),
    execution_mode::Union{Nothing,EMTPartitionExecutionMode} = nothing,
)
    start < stop || throw(ArgumentError(
        "partition stop must be later than start",
    ))
    communication_step > zero(communication_step) || throw(ArgumentError(
        "partition communication step must be positive",
    ))
    typed_regions = Tuple(regions)
    typed_ports = Tuple(ports)
    all(region -> region isa EMTPartitionRegion, typed_regions) || throw(ArgumentError(
        "partition regions must use EMTPartitionRegion",
    ))
    all(port -> port isa EMTInterfacePort, typed_ports) || throw(ArgumentError(
        "partition ports must use EMTInterfacePort",
    ))
    1 <= length(typed_regions) <= _MAXIMUM_REGIONS || throw(ArgumentError(
        "partition plan requires from 1 through 64 regions",
    ))
    region_names = getfield.(typed_regions, :identity)
    length(region_names) == length(unique(region_names)) || throw(ArgumentError(
        "partition plan repeats a region identity",
    ))
    model_names = reduce(vcat, (collect(region.model_identities) for region in typed_regions))
    length(model_names) == length(unique(model_names)) || throw(ArgumentError(
        "partition plan assigns one physical model to multiple regions",
    ))
    port_names = getfield.(typed_ports, :identity)
    length(port_names) == length(unique(port_names)) || throw(ArgumentError(
        "partition plan repeats an interface identity",
    ))
    for port in typed_ports
        port.positive_region in region_names && port.negative_region in region_names ||
            throw(ArgumentError("partition interface names an absent region"))
        coordinates = emt_interface_coordinates(port, 0.375, -0.125)
        physical = emt_interface_physical_values(
            port,
            coordinates.first,
            coordinates.second,
        )
        isapprox(physical.voltage_v, 0.375; rtol=8eps(Float64), atol=0.0) &&
            isapprox(
                physical.outward_current_a,
                -0.125;
                rtol=8eps(Float64),
                atol=0.0,
            ) ||
            throw(ArgumentError(
                "partition equivalent-port mapping is not exactly reversible",
            ))
    end
    ratios = Tuple(map(typed_regions) do region
        ratio = _logical_ratio(communication_step, region.local_step)
        ratio === nothing && throw(ArgumentError(
            "partition local steps must divide the communication step exactly",
        ))
        1 <= ratio <= _MAXIMUM_RATE_RATIO || throw(ArgumentError(
            "partition local-to-communication rate ratio must be from 1 through 10,000",
        ))
        ratio
    end)
    length(unique(ratios)) <= _MAXIMUM_LOCAL_RATES || throw(ArgumentError(
        "partition plan exceeds 16 distinct local rates",
    ))
    window_count = _logical_ratio(stop - start, communication_step)
    window_count === nothing && throw(ArgumentError(
        "partition horizon must contain an exact integer number of communication windows",
    ))
    window_count > 0 || throw(ArgumentError(
        "partition horizon must contain at least one communication window",
    ))
    if exchange.method == DirectCoupledExchange
        length(unique(getfield.(typed_regions, :local_step))) == 1 || throw(ArgumentError(
            "direct coupled exchange requires one equal regional step",
        ))
    end
    selected_mode = if execution_mode === nothing
        if exchange.method == DirectCoupledExchange
            MonolithicReferenceExecution
        elseif exchange.method == LaggedCausalExchange
            any(>(1), ratios) ?
                CombinedLocalPartitionedExecution : PartitionedLaggedExecution
        else
            any(>(1), ratios) ?
                CombinedLocalPartitionedExecution : PartitionedWaveformExecution
        end
    else
        execution_mode
    end
    if selected_mode == MonolithicReferenceExecution
        exchange.method == DirectCoupledExchange || throw(ArgumentError(
            "monolithic reference execution requires direct coupled exchange",
        ))
    elseif selected_mode == LocalSubcyclingExecution
        length(typed_regions) == 1 && isempty(typed_ports) ||
            throw(ArgumentError(
                "local subcycling requires one connected region without a cut port",
            ))
        exchange.method == DirectCoupledExchange || throw(ArgumentError(
            "local subcycling uses the connected direct network owner",
        ))
    elseif selected_mode == PartitionedLaggedExecution
        length(typed_regions) > 1 && !isempty(typed_ports) &&
            exchange.method == LaggedCausalExchange && all(==(1), ratios) ||
            throw(ArgumentError(
                "partitioned lagged execution requires multiple equal-rate regions and causal lagged exchange",
            ))
    elseif selected_mode == PartitionedWaveformExecution
        length(typed_regions) > 1 && !isempty(typed_ports) &&
            exchange.method == IteratedWaveformExchange && all(==(1), ratios) ||
            throw(ArgumentError(
                "partitioned waveform execution requires multiple equal-rate regions and iterated exchange",
            ))
    elseif selected_mode == CombinedLocalPartitionedExecution
        length(typed_regions) > 1 && !isempty(typed_ports) &&
            exchange.method in (IteratedWaveformExchange, LaggedCausalExchange) &&
            any(>(1), ratios) || throw(ArgumentError(
                "combined local/partitioned execution requires multiple regions and at least one local rate above one",
            ))
    end
    signature = bytes2hex(sha256(join(
        _partition_signature_lines(
            selected_mode,
            start,
            stop,
            communication_step,
            typed_regions,
            typed_ports,
            exchange,
            ratios,
        ),
        '\n',
    )))
    return EMTPartitionPlan(
        selected_mode,
        start,
        stop,
        communication_step,
        typed_regions,
        typed_ports,
        exchange,
        ratios,
        window_count,
        signature,
    )
end

partition_plan_signature_sha256(plan::EMTPartitionPlan) = plan.signature_sha256

"""A synthetic passive source-RL and load-RC network split by one oriented voltage/current interface."""
struct PassiveTwoRegionRLCStudy{P<:EMTPartitionPlan}
    plan::P
    source_voltage_v::Float64
    source_resistance_ohm::Float64
    source_inductance_h::Float64
    load_resistance_ohm::Float64
    load_capacitance_f::Float64
    initial_source_current_a::Float64
    initial_interface_voltage_v::Float64
    signature_sha256::String

    function PassiveTwoRegionRLCStudy(
        plan::P;
        source_voltage_v::Real,
        source_resistance_ohm::Real,
        source_inductance_h::Real,
        load_resistance_ohm::Real,
        load_capacitance_f::Real,
        initial_source_current_a::Real = 0.0,
        initial_interface_voltage_v::Real = 0.0,
    ) where {P<:EMTPartitionPlan}
        length(plan.regions) == 2 || throw(ArgumentError(
            "passive two-region study requires exactly two regions",
        ))
        length(plan.ports) == 1 || throw(ArgumentError(
            "passive two-region study requires exactly one interface",
        ))
        source_voltage = Float64(source_voltage_v)
        source_resistance = _nonnegative_finite(
            source_resistance_ohm,
            "source resistance",
        )
        source_inductance = _positive_finite(source_inductance_h, "source inductance")
        load_resistance = _positive_finite(load_resistance_ohm, "load resistance")
        load_capacitance = _positive_finite(load_capacitance_f, "load capacitance")
        initial_current = Float64(initial_source_current_a)
        initial_voltage = Float64(initial_interface_voltage_v)
        all(isfinite, (source_voltage, initial_current, initial_voltage)) || throw(
            ArgumentError("passive two-region sources and initial state must be finite"),
        )
        lines = (
            "schema=aimora.emt.passive_two_region_rlc.v1",
            "plan=$(plan.signature_sha256)",
            "source_voltage_v=$(repr(source_voltage))",
            "source_resistance_ohm=$(repr(source_resistance))",
            "source_inductance_h=$(repr(source_inductance))",
            "load_resistance_ohm=$(repr(load_resistance))",
            "load_capacitance_f=$(repr(load_capacitance))",
            "initial_source_current_a=$(repr(initial_current))",
            "initial_interface_voltage_v=$(repr(initial_voltage))",
        )
        signature = bytes2hex(sha256(join(lines, '\n')))
        return new{P}(
            plan,
            source_voltage,
            source_resistance,
            source_inductance,
            load_resistance,
            load_capacitance,
            initial_current,
            initial_voltage,
            signature,
        )
    end
end

partition_study_signature_sha256(study::PassiveTwoRegionRLCStudy) =
    study.signature_sha256

struct EMTPartitionFailure <: Exception
    code::Symbol
    last_accepted_time_s::Float64
    message::String
end

function Base.showerror(io::IO, failure::EMTPartitionFailure)
    print(
        io,
        String(failure.code),
        ": ",
        failure.message,
        " [last_accepted_time_s=",
        failure.last_accepted_time_s,
        ']'
    )
end

struct EMTPartitionCheckpoint
    schema::Symbol
    execution_mode::EMTPartitionExecutionMode
    plan_signature_sha256::String
    study_signature_sha256::String
    accepted_window_count::Int
    rejected_window_count::Int
    source_local_step_count::Int
    load_local_step_count::Int
    time_s::Float64
    source_current_a::Float64
    interface_voltage_v::Float64
    load_capacitor_current_a::Float64
    previous_voltage_slope_v_per_s::Float64
    time_trace_s::Tuple
    current_trace_a::Tuple
    voltage_trace_v::Tuple
    voltage_residual_trace_v::Tuple
    kcl_residual_trace_a::Tuple
    interface_energy_defect_trace_j::Tuple
    interface_value_age_trace_s::Tuple
    communication_error_estimate_trace_v::Tuple
    topology_revalidation_trace::Tuple
    fixed_point_iteration_trace::Tuple
    refinement_count::Int
    fallback_count::Int
    last_failure::Union{Nothing,String}
    signature_sha256::String
end

function partition_checkpoint_signature_sha256(checkpoint::EMTPartitionCheckpoint)
    return checkpoint.signature_sha256
end

"""One accepted multi-region synchronization point with portable regional state."""
struct EMTDeckPartitionCheckpoint
    schema::Symbol
    execution_mode::EMTPartitionExecutionMode
    plan_signature_sha256::String
    study_signature_sha256::String
    region_identities::Tuple
    port_identities::Tuple
    accepted_window_count::Int
    rejected_window_count::Int
    regional_local_step_counts::Tuple
    time_s::Float64
    interface_current_a::Tuple
    positive_terminal_voltage_v::Tuple
    negative_terminal_voltage_v::Tuple
    initialization_voltage_residual_v::Tuple
    initialization_consistent::Bool
    time_trace_s::Tuple
    positive_terminal_voltage_trace_v::Tuple
    negative_terminal_voltage_trace_v::Tuple
    interface_current_trace_a::Tuple
    voltage_residual_trace_v::Tuple
    kcl_residual_trace_a::Tuple
    interface_energy_defect_trace_j::Tuple
    interface_value_age_trace_s::Tuple
    communication_error_estimate_trace_v::Tuple
    topology_revalidation_trace::Tuple
    fixed_point_iteration_trace::Tuple
    refinement_count::Int
    fallback_count::Int
    regional_snapshots::Tuple
    regional_runtime_snapshots::Tuple
    last_failure::Union{Nothing,String}
    signature_sha256::String
end

function partition_checkpoint_signature_sha256(
    checkpoint::EMTDeckPartitionCheckpoint,
)
    return checkpoint.signature_sha256
end

struct EMTPartitionResult
    execution_mode::EMTPartitionExecutionMode
    plan_signature_sha256::String
    study_signature_sha256::String
    time_s::Vector{Float64}
    source_current_a::Vector{Float64}
    interface_voltage_v::Vector{Float64}
    voltage_residual_v::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    interface_energy_defect_j::Vector{Float64}
    interface_value_age_s::Vector{Float64}
    communication_error_estimate_v::Vector{Float64}
    topology_revalidated::Vector{Bool}
    fixed_point_iterations::Vector{Int}
    accepted_window_count::Int
    rejected_window_count::Int
    source_local_step_count::Int
    load_local_step_count::Int
    refinement_count::Int
    fallback_count::Int
    last_failure::Union{Nothing,String}
    accepted::Bool
    deterministic_signature_sha256::String
end

function partition_result(
    execution_mode::EMTPartitionExecutionMode,
    plan_signature::AbstractString,
    study_signature::AbstractString,
    time_s,
    source_current_a,
    interface_voltage_v,
    voltage_residual_v,
    kcl_residual_a,
    interface_energy_defect_j,
    interface_value_age_s,
    communication_error_estimate_v,
    topology_revalidated,
    fixed_point_iterations;
    accepted_window_count::Integer,
    rejected_window_count::Integer,
    source_local_step_count::Integer,
    load_local_step_count::Integer,
    refinement_count::Integer=0,
    fallback_count::Integer=0,
    last_failure::Union{Nothing,AbstractString}=nothing,
    accepted::Bool,
)
    times = Float64[time_s...]
    currents = Float64[source_current_a...]
    voltages = Float64[interface_voltage_v...]
    voltage_residuals = Float64[voltage_residual_v...]
    kcl_residuals = Float64[kcl_residual_a...]
    energy_defects = Float64[interface_energy_defect_j...]
    value_ages = Float64[interface_value_age_s...]
    communication_errors = Float64[communication_error_estimate_v...]
    topology_revalidations = Bool[topology_revalidated...]
    iterations = Int[fixed_point_iterations...]
    sample_count = length(times)
    length(currents) == sample_count && length(voltages) == sample_count || throw(
        DimensionMismatch("partition state traces must have equal lengths"),
    )
    diagnostic_count = max(sample_count - 1, 0)
    all(length(values) == diagnostic_count for values in (
        voltage_residuals,
        kcl_residuals,
        energy_defects,
        value_ages,
        communication_errors,
        topology_revalidations,
        iterations,
    )) || throw(DimensionMismatch(
        "partition window diagnostics must contain one row per accepted interval",
    ))
    all(isfinite, vcat(
        times,
        currents,
        voltages,
        voltage_residuals,
        kcl_residuals,
        energy_defects,
        value_ages,
        communication_errors,
    )) || throw(ArgumentError("partition result traces must be finite"))
    lines = String[
        "schema=aimora.emt.partition_result.v2",
        "execution_mode=$(UInt8(execution_mode))",
        "plan=$(String(plan_signature))",
        "study=$(String(study_signature))",
        "accepted_window_count=$(Int(accepted_window_count))",
        "rejected_window_count=$(Int(rejected_window_count))",
        "source_local_step_count=$(Int(source_local_step_count))",
        "load_local_step_count=$(Int(load_local_step_count))",
        "refinement_count=$(Int(refinement_count))",
        "fallback_count=$(Int(fallback_count))",
        "last_failure=$(repr(last_failure))",
        "accepted=$accepted",
    ]
    for index in eachindex(times)
        push!(lines, "state=$index|$(repr(times[index]))|$(repr(currents[index]))|$(repr(voltages[index]))")
    end
    for index in eachindex(iterations)
        push!(lines, "window=$index|$(repr(voltage_residuals[index]))|$(repr(kcl_residuals[index]))|$(repr(energy_defects[index]))|$(repr(value_ages[index]))|$(repr(communication_errors[index]))|$(topology_revalidations[index])|$(iterations[index])")
    end
    signature = bytes2hex(sha256(join(lines, '\n')))
    return EMTPartitionResult(
        execution_mode,
        String(plan_signature),
        String(study_signature),
        times,
        currents,
        voltages,
        voltage_residuals,
        kcl_residuals,
        energy_defects,
        value_ages,
        communication_errors,
        topology_revalidations,
        iterations,
        Int(accepted_window_count),
        Int(rejected_window_count),
        Int(source_local_step_count),
        Int(load_local_step_count),
        Int(refinement_count),
        Int(fallback_count),
        last_failure === nothing ? nothing : String(last_failure),
        accepted,
        signature,
    )
end

end

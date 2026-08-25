export ConverterSystemFamily,
       SinglePhaseDiodeBridge,
       ThreePhaseDiodeBridge,
       SinglePhaseThyristorBridge,
       ThreePhaseThyristorBridge,
       SinglePhaseHalfControlledBridge,
       ThreePhaseHalfControlledBridge,
       MultipulseDiodeBridge,
       MultipulseThyristorBridge,
       BuckChopper,
       BoostChopper,
       InvertingBuckBoostChopper,
       FourQuadrantChopper,
       InterleavedChopper,
       DualActiveBridge,
       SinglePhaseTwoLevelBridge,
       ThreePhaseTwoLevelBridge,
       ThreeLevelNeutralPointClampedBridge,
       ThreeLevelTTypeBridge,
       FlyingCapacitorBridge,
       CascadedHBridge,
       ThreePhaseMatrixConverter,
       LineCommutatedCycloconverter,
       ConverterApplication,
       StandaloneConversion,
       ShuntActiveHarmonicFilter,
       DynamicVoltageRestorer,
       DoubleConversionUninterruptiblePowerSupply,
       ConductiveElectricVehicleCharger,
       ThreeStageSolidStateTransformer,
       ConverterModulationKind,
       NaturalDiodeCommutation,
       PhaseControlledFiring,
       CarrierSinusoidalPulseWidthModulation,
       SpaceVectorPulseWidthModulation,
       PhaseShiftedCarrierPulseWidthModulation,
       SelectiveHarmonicElimination,
       NearestLevelModulation,
       SinglePhaseShiftModulation,
       DualPhaseShiftModulation,
       TriplePhaseShiftModulation,
       MatrixSpaceVectorModulation,
       CycloconverterFiringSynthesis,
       ConverterPortKind,
       AlternatingCurrentPort,
       DirectCurrentPort,
       IsolatedDirectCurrentPort,
       ConverterPortDefinition,
       ConverterRatedBases,
       ConverterTimingParameters,
       ConverterModulationParameters,
       ConverterSystemEventKind,
       ConverterSystemEventCommand,
       ConverterBlockEvent,
       ConverterRestartEvent,
       validate_converter_system_event_calendar,
       converter_system_commands_at_boundary,
       converter_system_is_blocked,
       ConverterSystemSelection,
       ConverterSystemSpecification,
       converter_system_signature

@enum ConverterSystemFamily begin
    SinglePhaseDiodeBridge
    ThreePhaseDiodeBridge
    SinglePhaseThyristorBridge
    ThreePhaseThyristorBridge
    SinglePhaseHalfControlledBridge
    ThreePhaseHalfControlledBridge
    MultipulseDiodeBridge
    MultipulseThyristorBridge
    BuckChopper
    BoostChopper
    InvertingBuckBoostChopper
    FourQuadrantChopper
    InterleavedChopper
    DualActiveBridge
    SinglePhaseTwoLevelBridge
    ThreePhaseTwoLevelBridge
    ThreeLevelNeutralPointClampedBridge
    ThreeLevelTTypeBridge
    FlyingCapacitorBridge
    CascadedHBridge
    ThreePhaseMatrixConverter
    LineCommutatedCycloconverter
end

@enum ConverterApplication begin
    StandaloneConversion
    ShuntActiveHarmonicFilter
    DynamicVoltageRestorer
    DoubleConversionUninterruptiblePowerSupply
    ConductiveElectricVehicleCharger
    ThreeStageSolidStateTransformer
end

@enum ConverterModulationKind begin
    NaturalDiodeCommutation
    PhaseControlledFiring
    CarrierSinusoidalPulseWidthModulation
    SpaceVectorPulseWidthModulation
    PhaseShiftedCarrierPulseWidthModulation
    SelectiveHarmonicElimination
    NearestLevelModulation
    SinglePhaseShiftModulation
    DualPhaseShiftModulation
    TriplePhaseShiftModulation
    MatrixSpaceVectorModulation
    CycloconverterFiringSynthesis
end

@enum ConverterPortKind begin
    AlternatingCurrentPort
    DirectCurrentPort
    IsolatedDirectCurrentPort
end

struct ConverterPortDefinition
    identity::Symbol
    kind::ConverterPortKind
    ordered_nodes::Tuple{Vararg{Int}}
    voltage_orientation::String
    current_orientation::String

    function ConverterPortDefinition(
        identity::Symbol,
        kind::ConverterPortKind,
        ordered_nodes;
        voltage_orientation::AbstractString,
        current_orientation::AbstractString="positive into the converter system",
    )
        isempty(String(identity)) && throw(ArgumentError("converter port identity must not be empty"))
        nodes = Tuple(Int.(ordered_nodes))
        length(nodes) >= 2 || throw(ArgumentError("converter ports require at least two ordered nodes"))
        all(>=(0), nodes) || throw(ArgumentError("converter port nodes must be nonnegative"))
        length(unique(nodes)) == length(nodes) || throw(ArgumentError("converter port nodes must be unique"))
        isempty(strip(voltage_orientation)) && throw(ArgumentError("converter voltage orientation must not be empty"))
        isempty(strip(current_orientation)) && throw(ArgumentError("converter current orientation must not be empty"))
        return new(
            identity,
            kind,
            nodes,
            String(voltage_orientation),
            String(current_orientation),
        )
    end
end

struct ConverterRatedBases
    voltage_v::Float64
    current_a::Float64
    power_va::Float64
    frequency_hz::Float64

    function ConverterRatedBases(voltage_v, current_a, power_va, frequency_hz)
        values = Float64.((voltage_v, current_a, power_va, frequency_hz))
        all(isfinite, values) && all(>(0.0), values) || throw(ArgumentError(
            "converter rated voltage, current, power, and frequency must be finite and positive",
        ))
        return new(values...)
    end
end

Base.@kwdef struct ConverterTimingParameters
    fixed_step_s::Float64
    control_period_s::Float64
    carrier_frequency_hz::Float64 = 0.0
    firing_frequency_hz::Float64 = 0.0
    dead_time_s::Float64 = 0.0
    minimum_pulse_s::Float64 = 0.0
end

Base.@kwdef struct ConverterModulationParameters
    kind::ConverterModulationKind
    duty::Float64 = 0.5
    modulation_index::Float64 = 0.8
    firing_angle_rad::Float64 = 0.0
    phase_shift_rad::Float64 = 0.0
    primary_inner_phase_shift_rad::Float64 = 0.0
    secondary_inner_phase_shift_rad::Float64 = 0.0
    selective_harmonic_angles_rad::Tuple{Vararg{Float64}} = ()
end

@enum ConverterSystemEventKind begin
    ConverterBlockEvent
    ConverterRestartEvent
end

_converter_system_event_priority(kind::ConverterSystemEventKind) =
    kind === ConverterBlockEvent ? -20 : 10

struct ConverterSystemEventCommand
    id::Symbol
    kind::ConverterSystemEventKind
    time_s::Float64
    target_valve_indices::Tuple{Vararg{Int}}
    reference_id::Union{Nothing,Symbol}
    priority::Int
    topology_invalidating::Bool
end

function ConverterSystemEventCommand(
    id::Symbol,
    kind::ConverterSystemEventKind,
    time_s::Real;
    target_valve_indices=Int[],
    reference_id::Union{Nothing,Symbol}=nothing,
)
    isempty(String(id)) && throw(ArgumentError(
        "converter-system event identity must not be empty",
    ))
    time = Float64(time_s)
    isfinite(time) && time >= 0.0 || throw(ArgumentError(
        "converter-system event time must be finite and nonnegative",
    ))
    targets = Tuple(Int.(target_valve_indices))
    all(>(0), targets) && length(unique(targets)) == length(targets) ||
        isempty(targets) || throw(ArgumentError(
            "converter-system event target valve indices must be unique and positive",
        ))
    if kind === ConverterRestartEvent
        reference_id === nothing && throw(ArgumentError(
            "converter restart event must reference its block event",
        ))
    else
        reference_id === nothing || throw(ArgumentError(
            "converter block event cannot carry a clear-event reference",
        ))
    end
    return ConverterSystemEventCommand(
        id,
        kind,
        time,
        targets,
        reference_id,
        _converter_system_event_priority(kind),
        true,
    )
end

_converter_system_event_order_key(command::ConverterSystemEventCommand) =
    (command.time_s, command.priority, String(command.id))

function _validate_converter_system_event_commands(commands)
    ids = getfield.(commands, :id)
    length(unique(ids)) == length(ids) || throw(ArgumentError(
        "converter-system event identities must be unique",
    ))
    issorted(_converter_system_event_order_key.(commands)) || throw(ArgumentError(
        "converter-system events must retain time, physical-priority, and identity order",
    ))
    block_by_id = Dict{Symbol,ConverterSystemEventCommand}()
    restarted_block_ids = Set{Symbol}()
    for command in commands
        if command.kind === ConverterBlockEvent
            block_by_id[command.id] = command
        else
            block = get(block_by_id, command.reference_id, nothing)
            block === nothing && throw(ArgumentError(
                "converter restart event references an unavailable prior block event",
            ))
            block.time_s < command.time_s || throw(ArgumentError(
                "converter restart must occur after its referenced block event",
            ))
            block.target_valve_indices == command.target_valve_indices ||
                throw(ArgumentError(
                    "converter restart targets must match its referenced block event",
                ))
            command.reference_id in restarted_block_ids && throw(ArgumentError(
                "each converter block event may have only one restart event",
            ))
            push!(restarted_block_ids, something(command.reference_id))
        end
    end
    Set(keys(block_by_id)) == restarted_block_ids || throw(ArgumentError(
        "every converter block event must have exactly one later restart event",
    ))
    return commands
end

struct ConverterSystemSelection
    family::ConverterSystemFamily
    fidelity::ModelFidelity
    application::ConverterApplication
    phase_count::Int
    pulse_count::Int
    channel_count::Int
    cell_count::Int
    thermal_stage_count::Int

    function ConverterSystemSelection(
        family::ConverterSystemFamily,
        fidelity::ModelFidelity;
        application::ConverterApplication=StandaloneConversion,
        phase_count::Integer=3,
        pulse_count::Integer=6,
        channel_count::Integer=1,
        cell_count::Integer=1,
        thermal_stage_count::Integer=1,
    )
        counts = Int.((phase_count, pulse_count, channel_count, cell_count, thermal_stage_count))
        all(>(0), counts) || throw(ArgumentError("converter selection counts must be positive"))
        return new(family, fidelity, application, counts...)
    end
end

struct ConverterSystemSpecification{P<:Tuple,T<:Tuple,D<:Tuple,R<:Tuple,E<:Tuple}
    schema::Symbol
    identity::Symbol
    selection::ConverterSystemSelection
    ports::P
    topology_signatures::T
    device_fidelity_signatures::D
    passive_and_transformer_signatures::R
    rated_bases::ConverterRatedBases
    timing::ConverterTimingParameters
    modulation::ConverterModulationParameters
    event_commands::E
    provenance::Tuple{Vararg{ParameterProvenance}}
    validity_domain::ModelValidityDomain
    output_channels::Tuple{Vararg{Symbol}}
    signature_sha256::String
end

function _converter_system_signature(parts...)
    io = IOBuffer()
    for part in parts
        print(io, repr(part), '\n')
    end
    return bytes2hex(sha256(take!(io)))
end

function ConverterSystemSpecification(
    identity::Symbol,
    selection::ConverterSystemSelection,
    ports,
    topology_signatures,
    device_fidelity_signatures,
    passive_and_transformer_signatures,
    rated_bases::ConverterRatedBases,
    timing::ConverterTimingParameters,
    modulation::ConverterModulationParameters,
    provenance,
    validity_domain::ModelValidityDomain;
    event_commands=ConverterSystemEventCommand[],
    output_channels=(:terminal_voltage_v, :terminal_current_a, :power_w, :stored_energy_j),
)
    isempty(String(identity)) && throw(ArgumentError("converter system identity must not be empty"))
    port_values = Tuple(ports)
    !isempty(port_values) && all(port -> port isa ConverterPortDefinition, port_values) ||
        throw(ArgumentError("converter systems require typed electrical ports"))
    port_names = getfield.(port_values, :identity)
    length(unique(port_names)) == length(port_names) || throw(ArgumentError(
        "converter port identities must be unique",
    ))
    topologies = Tuple(String.(topology_signatures))
    devices = Tuple(String.(device_fidelity_signatures))
    physical_owners = Tuple(String.(passive_and_transformer_signatures))
    !isempty(topologies) || throw(ArgumentError("converter systems require canonical topology identities"))
    all(value -> !isempty(strip(value)), (topologies..., devices..., physical_owners...)) ||
        throw(ArgumentError("converter owner signatures must not be empty"))
    provenance_values = Tuple(provenance)
    !isempty(provenance_values) && all(value -> value isa ParameterProvenance, provenance_values) ||
        throw(ArgumentError("converter systems require typed parameter provenance"))
    outputs = Tuple(Symbol.(output_channels))
    !isempty(outputs) && length(unique(outputs)) == length(outputs) || throw(ArgumentError(
        "converter output channels must be nonempty and unique",
    ))
    commands = Tuple(event_commands)
    all(command -> command isa ConverterSystemEventCommand, commands) ||
        throw(ArgumentError("converter event calendar requires typed commands"))
    _validate_converter_system_event_commands(commands)
    signature = _converter_system_signature(
        :aimora_converter_system_v1,
        identity,
        selection,
        port_values,
        topologies,
        devices,
        physical_owners,
        rated_bases,
        timing,
        modulation,
        commands,
        provenance_values,
        validity_domain,
        outputs,
    )
    return ConverterSystemSpecification(
        :aimora_converter_system_v1,
        identity,
        selection,
        port_values,
        topologies,
        devices,
        physical_owners,
        rated_bases,
        timing,
        modulation,
        commands,
        provenance_values,
        validity_domain,
        outputs,
        signature,
    )
end

function _converter_system_same_event_boundary(specification, left, right)
    step = specification.timing.fixed_step_s
    return isapprox(
        left,
        right;
        atol=16.0 * eps(Float64) * max(abs(right), step),
        rtol=0.0,
    )
end

function converter_system_commands_at_boundary(
    specification::ConverterSystemSpecification,
    time_s::Real,
)
    time = Float64(time_s)
    isfinite(time) && time >= 0.0 || throw(ArgumentError(
        "converter-system event query time must be finite and nonnegative",
    ))
    return Tuple(command for command in specification.event_commands if
        _converter_system_same_event_boundary(specification, time, command.time_s))
end

function validate_converter_system_event_calendar(
    specification::ConverterSystemSpecification;
    start_time_s::Real,
    stop_time_s::Real,
    allowed_target_valve_indices,
)
    start_time = Float64(start_time_s)
    stop_time = Float64(stop_time_s)
    isfinite(start_time) && isfinite(stop_time) &&
        0.0 <= start_time < stop_time || throw(ArgumentError(
        "converter-system event horizon must be finite, nonnegative, and increasing",
    ))
    step = specification.timing.fixed_step_s
    isfinite(step) && step > 0.0 || throw(ArgumentError(
        "converter-system event calendar requires a finite positive fixed step",
    ))
    allowed_targets = Int.(collect(allowed_target_valve_indices))
    !isempty(allowed_targets) && all(>(0), allowed_targets) &&
        length(unique(allowed_targets)) == length(allowed_targets) ||
        throw(ArgumentError(
            "converter-system event calendar requires unique positive admitted valve indices",
    ))
    allowed_target_set = Set(allowed_targets)
    active_target_sets = Dict{Symbol,Set{Int}}()
    for command in specification.event_commands
        start_time < command.time_s <= stop_time || throw(ArgumentError(
            "converter-system events must occur after the study start and no later than its stop",
        ))
        step_index = (command.time_s - start_time) / step
        isapprox(step_index, round(step_index); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "converter-system events must lie exactly on the fixed-step calendar",
            ))
        targets = isempty(command.target_valve_indices) ?
            allowed_target_set : Set(command.target_valve_indices)
        !isempty(targets) && issubset(targets, allowed_target_set) ||
            throw(ArgumentError(
                "converter-system event targets are outside the study's admitted controlled valves",
            ))
        if command.kind === ConverterBlockEvent
            any(active_targets -> !isdisjoint(targets, active_targets),
                values(active_target_sets)) && throw(ArgumentError(
                "converter-system block intervals must not overlap on one controlled valve",
            ))
            active_target_sets[command.id] = targets
        else
            delete!(active_target_sets, something(command.reference_id))
        end
    end
    return specification
end

function converter_system_is_blocked(
    specification::ConverterSystemSpecification,
    time_s::Real,
    target_valve_index::Union{Nothing,Integer}=nothing,
)
    time = Float64(time_s)
    isfinite(time) && time >= 0.0 || throw(ArgumentError(
        "converter-system block-state query time must be finite and nonnegative",
    ))
    tolerance = 16.0 * eps(Float64) * max(abs(time), specification.timing.fixed_step_s)
    target = if target_valve_index === nothing
        nothing
    else
        index = Int(target_valve_index)
        index > 0 || throw(ArgumentError(
            "converter-system block-state target valve index must be positive",
        ))
        index
    end
    active_blocks = Dict{Symbol,ConverterSystemEventCommand}()
    for command in specification.event_commands
        command.time_s <= time + tolerance || break
        if command.kind === ConverterBlockEvent
            active_blocks[command.id] = command
        else
            delete!(active_blocks, something(command.reference_id))
        end
    end
    target === nothing && return !isempty(active_blocks)
    return any(values(active_blocks)) do block
        isempty(block.target_valve_indices) || target in block.target_valve_indices
    end
end

converter_system_signature(specification::ConverterSystemSpecification) =
    specification.signature_sha256

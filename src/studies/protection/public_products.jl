@enum ProtectionProductFamily::UInt8 begin
    RadialFeederProtectionProduct = 0x01
    DirectionalDistanceLineProtectionProduct = 0x02
    TransformerBusDifferentialProtectionProduct = 0x03
    MachineConverterFrequencyProtectionProduct = 0x04
    DCDifferentialOvercurrentProtectionProduct = 0x05
end

@enum ProtectionInitializationMode::UInt8 begin
    HealthyInactiveProtectionInitialization = 0x01
    DeclaredProtectionStateInitialization = 0x02
    OperatingPointProtectionInitialization = 0x03
end

const _PROTECTION_PRODUCT_REQUIRED_ELEMENTS = Dict(
    RadialFeederProtectionProduct => (
        :phase_overcurrent,
        :residual_earth_fault,
        :breaker_failure,
        :autoreclose,
    ),
    DirectionalDistanceLineProtectionProduct => (
        :directional,
        :distance,
        :incremental_wave,
        :communication,
    ),
    TransformerBusDifferentialProtectionProduct => (
        :transformer_differential,
        :bus_differential,
    ),
    MachineConverterFrequencyProtectionProduct => (
        :frequency,
        :rocof,
        :asset_logic,
    ),
    DCDifferentialOvercurrentProtectionProduct => (
        :dc_differential,
        :dc_overcurrent,
        :breaker_failure,
    ),
)

const _PROTECTION_PRODUCT_UNSUPPORTED = (
    :protected_standard_conformance,
    :vendor_relay_equivalence,
    :field_setting_coordination,
    :field_or_laboratory_timing,
    :telecommunication_protocol_conformance,
    :cybersecurity_assurance,
    :detailed_arc_or_insulation_physics,
    :hard_realtime_or_physical_hil,
    :atp_or_pscad_equivalence,
    :certification,
)

"""One explicit terminal whose current is positive into the protected zone."""
struct ProtectionTerminalDefinition
    id::Symbol
    asset_terminal::Symbol
    measurement_channels::Tuple{Vararg{Symbol}}
    inward_current_orientation::String

    function ProtectionTerminalDefinition(
        id::Symbol,
        asset_terminal::Symbol,
        measurement_channels::AbstractVector{Symbol};
        inward_current_orientation::AbstractString,
    )
        channels = Tuple(Symbol.(measurement_channels))
        all(name -> !isempty(String(name)), (id, asset_terminal)) &&
            !isempty(channels) && all(name -> !isempty(String(name)), channels) ||
            throw(ArgumentError("protection terminal identities must not be empty"))
        allunique(channels) || throw(ArgumentError(
            "protection terminal measurement channels must be unique",
        ))
        orientation = String(inward_current_orientation)
        isempty(strip(orientation)) && throw(ArgumentError(
            "protection terminal inward-current orientation must be explicit",
        ))
        return new(id, asset_terminal, channels, orientation)
    end
end

"""A generic, nonvendor public protection product with explicit owner identities."""
struct ProtectionProductSpecification
    id::Symbol
    family::ProtectionProductFamily
    protected_asset::Symbol
    protected_zone::Symbol
    terminals::Tuple{Vararg{ProtectionTerminalDefinition}}
    measurement_products::Tuple{Vararg{Symbol}}
    element_families::Tuple{Vararg{Symbol}}
    communication_links::Tuple{Vararg{Symbol}}
    trip_breakers::Tuple{Vararg{Symbol}}
    network_timestep_s::Float64
    network_timestep_logical::EMTLogicalTime
    protection_sample_period_ticks::Int
    minimum_frequency_hz::Float64
    maximum_frequency_hz::Float64
    setting_provenance::ParameterProvenance
    timing_provenance::ParameterProvenance
    uncertainty::String
    validity_domain::String
    unsupported_phenomena::Tuple{Vararg{Symbol}}
    configuration::NamedTuple
    deterministic_signature_sha256::String

    function ProtectionProductSpecification(
        id::Symbol,
        family::ProtectionProductFamily,
        protected_asset::Symbol,
        protected_zone::Symbol;
        terminals::AbstractVector{ProtectionTerminalDefinition},
        measurement_products::AbstractVector{Symbol},
        element_families::AbstractVector{Symbol},
        communication_links::AbstractVector{Symbol}=Symbol[],
        trip_breakers::AbstractVector{Symbol},
        network_timestep_s::Real,
        network_timestep_logical::EMTLogicalTime,
        protection_sample_period_ticks::Integer,
        minimum_frequency_hz::Real=45.0,
        maximum_frequency_hz::Real=65.0,
        setting_provenance::ParameterProvenance,
        timing_provenance::ParameterProvenance,
        uncertainty::AbstractString,
        validity_domain::AbstractString,
        unsupported_phenomena::AbstractVector{Symbol}=collect(_PROTECTION_PRODUCT_UNSUPPORTED),
        configuration::NamedTuple,
    )
        all(name -> !isempty(String(name)), (id, protected_asset, protected_zone)) ||
            throw(ArgumentError("protection product identities must not be empty"))
        terminal_tuple = Tuple(terminals)
        measurements = Tuple(Symbol.(measurement_products))
        elements = Tuple(Symbol.(element_families))
        links = Tuple(Symbol.(communication_links))
        breakers = Tuple(Symbol.(trip_breakers))
        unsupported = Tuple(Symbol.(unsupported_phenomena))
        !isempty(terminal_tuple) && !isempty(measurements) && !isempty(elements) &&
            !isempty(breakers) || throw(ArgumentError(
            "protection products require terminals, measurements, elements, and trip breakers",
        ))
        for identities in (
            getfield.(terminal_tuple, :id),
            measurements,
            elements,
            links,
            breakers,
            unsupported,
        )
            allunique(identities) || throw(ArgumentError(
                "protection product identity lists must not contain duplicates",
            ))
            all(name -> !isempty(String(name)), identities) || throw(ArgumentError(
                "protection product identity lists must not contain empty names",
            ))
        end
        required = _PROTECTION_PRODUCT_REQUIRED_ELEMENTS[family]
        all(in(elements), required) || throw(ArgumentError(
            "protection product omits a required family for $(family)",
        ))
        family === DirectionalDistanceLineProtectionProduct && isempty(links) &&
            throw(ArgumentError(
                "directional-distance line protection requires an explicit communication link",
            ))
        timestep = Float64(network_timestep_s)
        sample_ticks = Int(protection_sample_period_ticks)
        minimum_frequency = Float64(minimum_frequency_hz)
        maximum_frequency = Float64(maximum_frequency_hz)
        isfinite(timestep) && timestep > 0.0 && 1 <= sample_ticks <= 200_000 ||
            throw(ArgumentError("protection product clocks are outside their domain"))
        Float64(network_timestep_logical) == timestep || throw(ArgumentError(
            "protection product floating and exact network timesteps differ",
        ))
        isfinite(minimum_frequency) && isfinite(maximum_frequency) &&
            0.0 <= minimum_frequency <= maximum_frequency || throw(ArgumentError(
            "protection product frequency domain is invalid",
        ))
        setting_provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "protection product setting provenance must be physical",
        ))
        timing_provenance.nature === NumericalPolicyParameter || throw(ArgumentError(
            "protection product timing provenance must be numerical policy",
        ))
        uncertainty_string = String(uncertainty)
        validity_string = String(validity_domain)
        isempty(strip(uncertainty_string)) && throw(ArgumentError(
            "protection product uncertainty must be explicit",
        ))
        isempty(strip(validity_string)) && throw(ArgumentError(
            "protection product validity domain must be explicit",
        ))
        all(in(unsupported), _PROTECTION_PRODUCT_UNSUPPORTED) || throw(ArgumentError(
            "generic protection products must retain every global unsupported boundary",
        ))
        required_configuration = (
            :element_settings,
            :logic_definitions,
            :communication_links,
            :breaker_specifications,
        )
        all(key -> haskey(configuration, key), required_configuration) ||
            throw(ArgumentError(
                "protection product configuration omits an exact owner collection",
            ))
        element_settings = configuration.element_settings
        logic_definitions = configuration.logic_definitions
        communication_definitions = configuration.communication_links
        breaker_definitions = configuration.breaker_specifications
        element_settings isa Tuple && !isempty(element_settings) &&
            all(setting -> setting isa Union{
                MagnitudeRelaySettings,
                DirectionalRelaySettings,
                DistanceRelaySettings,
                DifferentialRelaySettings,
                ROCOFEstimatorSettings,
                IncrementalWaveSettings,
            }, element_settings) || throw(ArgumentError(
            "protection product element settings must use public typed relay owners",
        ))
        logic_definitions isa Tuple && !isempty(logic_definitions) &&
            all(definition -> definition isa ProtectionLogicDefinition, logic_definitions) ||
            throw(ArgumentError(
                "protection product logic must use public typed definitions",
            ))
        communication_definitions isa Tuple &&
            all(link -> link isa ProtectionCommunicationLink, communication_definitions) ||
            throw(ArgumentError(
                "protection product communication must use public typed links",
            ))
        breaker_definitions isa Tuple && !isempty(breaker_definitions) &&
            all(specification -> specification isa EMTBreakerSpecification, breaker_definitions) ||
            throw(ArgumentError(
                "protection product breakers must use public typed specifications",
            ))
        Tuple(getfield.(communication_definitions, :id)) == links || throw(ArgumentError(
            "protection product communication identities differ from their definitions",
        ))
        Tuple(getfield.(breaker_definitions, :id)) == breakers || throw(ArgumentError(
            "protection product breaker identities differ from their definitions",
        ))
        io = IOBuffer()
        println(io, "aimora-protection-product-v1")
        println(io, id, '|', UInt8(family), '|', protected_asset, '|', protected_zone)
        for terminal in terminal_tuple
            println(io, repr(terminal))
        end
        for identities in (measurements, elements, links, breakers, unsupported)
            println(io, join(String.(identities), ','))
        end
        println(io, bitstring(timestep), '|', sample_ticks)
        println(
            io,
            network_timestep_logical.numerator,
            '/',
            network_timestep_logical.denominator,
        )
        println(io, bitstring(minimum_frequency), '|', bitstring(maximum_frequency))
        println(io, repr(setting_provenance))
        println(io, repr(timing_provenance))
        println(io, uncertainty_string)
        println(io, validity_string)
        println(io, repr(configuration))
        signature = bytes2hex(sha256(take!(io)))
        return new(
            id,
            family,
            protected_asset,
            protected_zone,
            terminal_tuple,
            measurements,
            elements,
            links,
            breakers,
            timestep,
            network_timestep_logical,
            sample_ticks,
            minimum_frequency,
            maximum_frequency,
            setting_provenance,
            timing_provenance,
            uncertainty_string,
            validity_string,
            unsupported,
            configuration,
            signature,
        )
    end
end

struct ProtectionProductPreparation
    specification::ProtectionProductSpecification
    initialization_mode::ProtectionInitializationMode
    initial_tick::Int
    initially_closed_breakers::Tuple{Vararg{Symbol}}
    initially_locked_out_breakers::Tuple{Vararg{Symbol}}
    pending_messages::Tuple{Vararg{QueuedProtectionMessage}}
    preparation_signature_sha256::String
end

struct ProtectionProductReadiness
    ready::Bool
    code::Symbol
    product::Symbol
    family::ProtectionProductFamily
    production_backend_available::Bool
    represented_element_families::Tuple{Vararg{Symbol}}
    unsupported_phenomena::Tuple{Vararg{Symbol}}
    deterministic_signature_sha256::String
end

function prepare_protection_product(
    specification::ProtectionProductSpecification;
    initialization_mode::ProtectionInitializationMode=
        HealthyInactiveProtectionInitialization,
    initial_tick::Integer=0,
    initially_closed_breakers::AbstractVector{Symbol}=collect(specification.trip_breakers),
    initially_locked_out_breakers::AbstractVector{Symbol}=Symbol[],
    pending_messages::AbstractVector{QueuedProtectionMessage}=QueuedProtectionMessage[],
)
    tick = Int(initial_tick)
    tick >= 0 || throw(ArgumentError("protection product initial tick must be nonnegative"))
    closed = Tuple(Symbol.(initially_closed_breakers))
    locked = Tuple(Symbol.(initially_locked_out_breakers))
    messages = Tuple(pending_messages)
    sequences = getfield.(messages, :sequence_number)
    allunique(closed) && allunique(locked) && allunique(sequences) || throw(ArgumentError(
        "protection product initial state identities must be unique",
    ))
    all(in(specification.trip_breakers), closed) &&
        all(in(specification.trip_breakers), locked) || throw(ArgumentError(
        "protection product initial breaker state names an unknown breaker",
    ))
    isempty(intersect(Set(closed), Set(locked))) || throw(ArgumentError(
        "a protection product breaker cannot be both initially closed and locked out",
    ))
    all(>(0), sequences) || throw(ArgumentError(
        "pending protection message sequences must be positive",
    ))
    configured_links = Set(specification.communication_links)
    all(message -> message.link_id in configured_links, messages) || throw(ArgumentError(
        "pending protection message names an unknown communication link",
    ))
    initialization_mode === HealthyInactiveProtectionInitialization &&
        (!isempty(locked) || !isempty(messages)) && throw(ArgumentError(
            "healthy protection initialization cannot contain lockout or pending messages",
        ))
    io = IOBuffer()
    println(io, specification.deterministic_signature_sha256)
    println(io, UInt8(initialization_mode), '|', tick)
    println(io, join(String.(closed), ','))
    println(io, join(String.(locked), ','))
    for message in messages
        println(io, repr(message))
    end
    return ProtectionProductPreparation(
        specification,
        initialization_mode,
        tick,
        closed,
        locked,
        messages,
        bytes2hex(sha256(take!(io))),
    )
end

function protection_product_readiness(
    preparation::ProtectionProductPreparation;
    production_backend_available::Bool=false,
)
    specification = preparation.specification
    return ProtectionProductReadiness(
        true,
        production_backend_available ? :ready_for_coupled_execution :
            :ready_for_solver_free_inspection,
        specification.id,
        specification.family,
        production_backend_available,
        specification.element_families,
        specification.unsupported_phenomena,
        preparation.preparation_signature_sha256,
    )
end

struct ProtectionProductEvent
    tick::Int
    owner::Symbol
    kind::Symbol
    detail::Symbol

    function ProtectionProductEvent(
        tick::Integer,
        owner::Symbol,
        kind::Symbol,
        detail::Symbol,
    )
        exact_tick = Int(tick)
        exact_tick >= 0 && all(name -> !isempty(String(name)), (owner, kind, detail)) ||
            throw(ArgumentError("protection product event is invalid"))
        return new(exact_tick, owner, kind, detail)
    end
end

struct ProtectionStudyRefusal
    code::Symbol
    product::Symbol
    instant::Union{Nothing,EMTLogicalTime}
    last_accepted_tick::Int
    message::String

    function ProtectionStudyRefusal(
        code::Symbol,
        product::Symbol,
        message::AbstractString;
        instant::Union{Nothing,EMTLogicalTime}=nothing,
        last_accepted_tick::Integer=-1,
    )
        !isempty(String(code)) && !isempty(String(product)) &&
            !isempty(strip(String(message))) || throw(ArgumentError(
            "protection study refusal must retain code, product, and message",
        ))
        return new(code, product, instant, Int(last_accepted_tick), String(message))
    end
end

struct ProtectionStudyPreparation{P}
    product_preparation::ProtectionProductPreparation
    task_pipeline::P
    execution_instant::EMTLogicalTime
    topology_signature_sha256::String
    output_ids::Tuple{Vararg{Symbol}}
    deterministic_signature_sha256::String
end

struct ProtectionStudyResult{S}
    product::Symbol
    preparation_signature_sha256::String
    task_result::S
    represented_element_families::Tuple{Vararg{Symbol}}
    output_ids::Tuple{Vararg{Symbol}}
    event_trace::Vector{ProtectionProductEvent}
    warnings::Tuple{Vararg{Symbol}}
    deterministic_signature_sha256::String
end

function prepare_protection_study(
    product_preparation::ProtectionProductPreparation,
    task_pipeline::ProtectionTaskPipeline;
    execution_instant::EMTLogicalTime,
    topology_signature_sha256::AbstractString,
    output_ids::AbstractVector{Symbol},
)
    specification = product_preparation.specification
    try
        topology_signature = lowercase(String(topology_signature_sha256))
        occursin(r"^[0-9a-f]{64}$", topology_signature) || throw(ArgumentError(
            "protection study topology signature must be a 64-hex SHA-256",
        ))
        outputs = Tuple(Symbol.(output_ids))
        !isempty(outputs) && allunique(outputs) &&
            all(name -> !isempty(String(name)), outputs) || throw(ArgumentError(
            "protection study output identities must be nonempty and unique",
        ))
        task_pipeline.plan.start <= execution_instant <= task_pipeline.plan.stop ||
            throw(ArgumentError(
                "protection study execution instant is outside its exact task horizon",
            ))
        expected_period = specification.protection_sample_period_ticks *
            specification.network_timestep_logical
        all(entry -> entry.spec.period == expected_period, task_pipeline.plan.entries) ||
            throw(ArgumentError(
                "protection study task periods differ from the product sample calendar",
            ))
        io = IOBuffer()
        println(io, "aimora-protection-study-preparation-v1")
        println(io, product_preparation.preparation_signature_sha256)
        println(io, task_pipeline.deterministic_signature_sha256)
        println(io, task_pipeline.plan.signature_sha256)
        println(io, execution_instant.numerator, '/', execution_instant.denominator)
        println(io, topology_signature)
        println(io, join(String.(outputs), ','))
        return ProtectionStudyPreparation(
            product_preparation,
            task_pipeline,
            execution_instant,
            topology_signature,
            outputs,
            bytes2hex(sha256(take!(io))),
        )
    catch error
        error isa InterruptException && rethrow()
        return ProtectionStudyRefusal(
            :invalid_preparation,
            specification.id,
            sprint(showerror, error);
            instant=execution_instant,
            last_accepted_tick=product_preparation.initial_tick,
        )
    end
end

function run_protection(preparation::ProtectionStudyPreparation)
    product = preparation.product_preparation.specification
    aimora = parentmodule(@__MODULE__)
    getfield(aimora, :solver_available)() || return ProtectionStudyRefusal(
        :production_backend_unavailable,
        product.id,
        "Coupled protection execution requires an explicitly activated production backend.";
        instant=preparation.execution_instant,
        last_accepted_tick=preparation.product_preparation.initial_tick,
    )
    try
        prepared_pipeline = getfield(aimora, :prepare_protection_task_pipeline)(
            preparation.task_pipeline,
        )
        task_result = getfield(aimora, :advance_protection_task_pipeline!)(
            prepared_pipeline,
            preparation.execution_instant,
        )
        state = task_result.state
        events = hasproperty(state, :product_runtime) ?
            copy(state.product_runtime.event_trace) : ProtectionProductEvent[]
        io = IOBuffer()
        println(io, "aimora-protection-study-result-v1")
        println(io, preparation.deterministic_signature_sha256)
        println(io, task_result.deterministic_signature_sha256)
        println(io, join(String.(product.element_families), ','))
        println(io, join(String.(preparation.output_ids), ','))
        for event in events
            println(io, repr(event))
        end
        return ProtectionStudyResult(
            product.id,
            preparation.deterministic_signature_sha256,
            task_result,
            product.element_families,
            preparation.output_ids,
            events,
            (:generic_noncertifying_settings,),
            bytes2hex(sha256(take!(io))),
        )
    catch error
        error isa InterruptException && rethrow()
        return ProtectionStudyRefusal(
            :execution_failed,
            product.id,
            sprint(showerror, error);
            instant=preparation.execution_instant,
            last_accepted_tick=preparation.product_preparation.initial_tick,
        )
    end
end

mutable struct ProtectionProductRuntime{E,L,C,B}
    preparation::ProtectionProductPreparation
    element_states::E
    logic_states::L
    communication_states::C
    breaker_states::B
    task_plan_signature_sha256::String
    accepted_tick::Int
    output_cursor::Int
    event_trace::Vector{ProtectionProductEvent}
end

struct ProtectionProductRuntimeSnapshot{E,L,C,B}
    preparation_signature_sha256::String
    task_plan_signature_sha256::String
    accepted_tick::Int
    output_cursor::Int
    element_snapshots::E
    logic_snapshots::L
    communication_snapshots::C
    breaker_snapshots::B
    event_trace::Vector{ProtectionProductEvent}
    deterministic_signature_sha256::String
end

_protection_element_state(::MagnitudeRelaySettings) = MagnitudeRelayState()
_protection_element_state(::DirectionalRelaySettings) = DirectionalRelayState()
_protection_element_state(::DistanceRelaySettings) = DistanceRelayState()
_protection_element_state(::DifferentialRelaySettings) = DifferentialRelayState()
_protection_element_state(::ROCOFEstimatorSettings) = ROCOFEstimatorState()
_protection_element_state(::IncrementalWaveSettings) = IncrementalWaveState()

_protection_element_snapshot(state::MagnitudeRelayState, settings::MagnitudeRelaySettings) =
    magnitude_relay_snapshot(state, settings)
_protection_element_snapshot(state::DirectionalRelayState, settings::DirectionalRelaySettings) =
    directional_relay_snapshot(state, settings)
_protection_element_snapshot(state::DistanceRelayState, settings::DistanceRelaySettings) =
    distance_relay_snapshot(state, settings)
_protection_element_snapshot(state::DifferentialRelayState, settings::DifferentialRelaySettings) =
    differential_relay_snapshot(state, settings)
_protection_element_snapshot(state::ROCOFEstimatorState, settings::ROCOFEstimatorSettings) =
    rocof_estimator_snapshot(state, settings)
_protection_element_snapshot(state::IncrementalWaveState, settings::IncrementalWaveSettings) =
    incremental_wave_snapshot(state, settings)

_restore_protection_element_snapshot!(state::MagnitudeRelayState, settings::MagnitudeRelaySettings, snapshot) =
    restore_magnitude_relay_snapshot!(state, settings, snapshot)
_restore_protection_element_snapshot!(state::DirectionalRelayState, settings::DirectionalRelaySettings, snapshot) =
    restore_directional_relay_snapshot!(state, settings, snapshot)
_restore_protection_element_snapshot!(state::DistanceRelayState, settings::DistanceRelaySettings, snapshot) =
    restore_distance_relay_snapshot!(state, settings, snapshot)
_restore_protection_element_snapshot!(state::DifferentialRelayState, settings::DifferentialRelaySettings, snapshot) =
    restore_differential_relay_snapshot!(state, settings, snapshot)
_restore_protection_element_snapshot!(state::ROCOFEstimatorState, settings::ROCOFEstimatorSettings, snapshot) =
    restore_rocof_estimator_snapshot!(state, settings, snapshot)
_restore_protection_element_snapshot!(state::IncrementalWaveState, settings::IncrementalWaveSettings, snapshot) =
    restore_incremental_wave_snapshot!(state, settings, snapshot)

function ProtectionProductRuntime(
    preparation::ProtectionProductPreparation;
    task_plan_signature_sha256::AbstractString,
)
    task_signature = lowercase(String(task_plan_signature_sha256))
    occursin(r"^[0-9a-f]{64}$", task_signature) || throw(ArgumentError(
        "protection product task-plan signature must be a 64-hex SHA-256",
    ))
    configuration = preparation.specification.configuration
    element_states = Tuple(_protection_element_state(setting) for
                           setting in configuration.element_settings)
    logic_states = Tuple(ProtectionLogicRuntime(definition) for
                         definition in configuration.logic_definitions)
    communication_states = Tuple(
        ProtectionCommunicationRuntime() for _ in configuration.communication_links
    )
    for (link, state) in zip(configuration.communication_links, communication_states)
        queued = QueuedProtectionMessage[
            message for message in preparation.pending_messages if message.link_id == link.id
        ]
        sort!(queued; by=message -> (message.delivery_tick, message.sequence_number, message.copy_index))
        append!(state.queue, queued)
        if !isempty(queued)
            state.next_sequence_number = maximum(getfield.(queued, :sequence_number)) + 1
            state.last_send_tick = maximum(getfield.(queued, :send_tick))
            state.sent_count = length(unique(getfield.(queued, :sequence_number)))
        end
    end
    closed = Set(preparation.initially_closed_breakers)
    locked = Set(preparation.initially_locked_out_breakers)
    breaker_states = Tuple(
        EMTBreakerRuntime(
            specification;
            tick_s=preparation.specification.network_timestep_s,
            initially_closed=specification.id in closed,
        ) for specification in configuration.breaker_specifications
    )
    for (specification, state) in zip(configuration.breaker_specifications, breaker_states)
        state.lockout = specification.id in locked
    end
    return ProtectionProductRuntime(
        preparation,
        element_states,
        logic_states,
        communication_states,
        breaker_states,
        task_signature,
        preparation.initial_tick,
        0,
        ProtectionProductEvent[],
    )
end

function _protection_product_snapshot_signature(
    preparation_signature,
    task_plan_signature,
    accepted_tick,
    output_cursor,
    element_snapshots,
    logic_snapshots,
    communication_snapshots,
    breaker_snapshots,
    event_trace,
)
    io = IOBuffer()
    println(io, "aimora-protection-product-runtime-snapshot-v1")
    println(io, preparation_signature)
    println(io, task_plan_signature)
    println(io, accepted_tick, '|', output_cursor)
    for collection in (
        element_snapshots,
        logic_snapshots,
        communication_snapshots,
        breaker_snapshots,
        event_trace,
    )
        for value in collection
            println(io, repr(value))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function protection_product_runtime_snapshot(runtime::ProtectionProductRuntime)
    configuration = runtime.preparation.specification.configuration
    element_snapshots = Tuple(
        _protection_element_snapshot(state, settings) for
        (state, settings) in zip(runtime.element_states, configuration.element_settings)
    )
    logic_snapshots = Tuple(
        protection_logic_snapshot(state, definition) for
        (state, definition) in zip(runtime.logic_states, configuration.logic_definitions)
    )
    communication_snapshots = Tuple(
        protection_communication_snapshot(state, link) for
        (state, link) in zip(runtime.communication_states, configuration.communication_links)
    )
    breaker_snapshots = Tuple(
        emt_breaker_snapshot(state, specification) for
        (state, specification) in zip(runtime.breaker_states, configuration.breaker_specifications)
    )
    events = copy(runtime.event_trace)
    signature = _protection_product_snapshot_signature(
        runtime.preparation.preparation_signature_sha256,
        runtime.task_plan_signature_sha256,
        runtime.accepted_tick,
        runtime.output_cursor,
        element_snapshots,
        logic_snapshots,
        communication_snapshots,
        breaker_snapshots,
        events,
    )
    return ProtectionProductRuntimeSnapshot(
        runtime.preparation.preparation_signature_sha256,
        runtime.task_plan_signature_sha256,
        runtime.accepted_tick,
        runtime.output_cursor,
        element_snapshots,
        logic_snapshots,
        communication_snapshots,
        breaker_snapshots,
        events,
        signature,
    )
end

function restore_protection_product_runtime_snapshot!(
    runtime::ProtectionProductRuntime,
    snapshot::ProtectionProductRuntimeSnapshot,
)
    snapshot.preparation_signature_sha256 ==
        runtime.preparation.preparation_signature_sha256 || throw(ArgumentError(
        "protection product snapshot preparation identity is stale",
    ))
    snapshot.task_plan_signature_sha256 == runtime.task_plan_signature_sha256 ||
        throw(ArgumentError("protection product snapshot task-plan identity is stale"))
    snapshot.accepted_tick >= runtime.preparation.initial_tick && snapshot.output_cursor >= 0 ||
        throw(ArgumentError("protection product snapshot cursors are invalid"))
    expected_signature = _protection_product_snapshot_signature(
        snapshot.preparation_signature_sha256,
        snapshot.task_plan_signature_sha256,
        snapshot.accepted_tick,
        snapshot.output_cursor,
        snapshot.element_snapshots,
        snapshot.logic_snapshots,
        snapshot.communication_snapshots,
        snapshot.breaker_snapshots,
        snapshot.event_trace,
    )
    expected_signature == snapshot.deterministic_signature_sha256 || throw(ArgumentError(
        "protection product snapshot content identity is corrupt",
    ))
    configuration = runtime.preparation.specification.configuration
    counts = (
        length(snapshot.element_snapshots) == length(runtime.element_states),
        length(snapshot.logic_snapshots) == length(runtime.logic_states),
        length(snapshot.communication_snapshots) == length(runtime.communication_states),
        length(snapshot.breaker_snapshots) == length(runtime.breaker_states),
    )
    all(counts) || throw(ArgumentError("protection product snapshot owner count changed"))
    trial = deepcopy(runtime)
    for (state, settings, state_snapshot) in zip(
        trial.element_states,
        configuration.element_settings,
        snapshot.element_snapshots,
    )
        _restore_protection_element_snapshot!(state, settings, state_snapshot)
    end
    for (state, definition, state_snapshot) in zip(
        trial.logic_states,
        configuration.logic_definitions,
        snapshot.logic_snapshots,
    )
        restore_protection_logic_snapshot!(state, definition, state_snapshot)
    end
    for (state, link, state_snapshot) in zip(
        trial.communication_states,
        configuration.communication_links,
        snapshot.communication_snapshots,
    )
        restore_protection_communication_snapshot!(state, link, state_snapshot)
    end
    for (state, specification, state_snapshot) in zip(
        trial.breaker_states,
        configuration.breaker_specifications,
        snapshot.breaker_snapshots,
    )
        restore_emt_breaker_snapshot!(state, specification, state_snapshot)
    end
    trial.accepted_tick = snapshot.accepted_tick
    trial.output_cursor = snapshot.output_cursor
    trial.event_trace = copy(snapshot.event_trace)
    protection_product_runtime_snapshot(trial).deterministic_signature_sha256 ==
        snapshot.deterministic_signature_sha256 || throw(ArgumentError(
        "protection product snapshot reconstruction is not exact",
    ))
    runtime.element_states = trial.element_states
    runtime.logic_states = trial.logic_states
    runtime.communication_states = trial.communication_states
    runtime.breaker_states = trial.breaker_states
    runtime.accepted_tick = trial.accepted_tick
    runtime.output_cursor = trial.output_cursor
    runtime.event_trace = trial.event_trace
    return runtime
end

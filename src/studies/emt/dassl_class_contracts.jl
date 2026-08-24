export DASSLClassEMTDiagnostics,
       DASSLClassEMTFailure,
       DASSLClassEMTNetworkRequest,
       DASSLClassEMTOwnerDisposition,
       DASSLClassEMTReadiness,
       DASSLClassEMTResult,
       DASSLClassEMTSettings,
       DASSLClassEMTStateLayout,
       DASSLClassEMTValidationRoot,
       DASSLClassEMTValidationTask,
       DASSLClassEMTValidationProblem,
       DASSLClassVariableStep,
       EMTIntegrationMode,
       EMTIntegrationSelection,
       FixedStepEMT,
       dassl_class_checkpoint_boundary!,
       dassl_class_emt_readiness,
       execute_dassl_class_emt!,
       prepare_dassl_class_emt

@enum EMTIntegrationMode::UInt8 begin
    FixedStepEMT = 0x01
    DASSLClassVariableStep = 0x02
end

"""Explicit EMT integration selection; the zero-argument form preserves fixed-step EMT."""
struct EMTIntegrationSelection{S}
    mode::EMTIntegrationMode
    settings::S
end

EMTIntegrationSelection() = EMTIntegrationSelection(FixedStepEMT, nothing)

"""Public numerical policy for the optional DASSL-class variable-step mode."""
struct DASSLClassEMTSettings
    initial_step_s::Float64
    minimum_step_s::Float64
    maximum_step_s::Float64
    maximum_order::Int
    relative_tolerance::Float64
    maximum_newton_iterations::Int
    maximum_rejected_steps::Int
    root_time_tolerance_s::Float64
    residual_wrms_tolerance::Float64
    correction_wrms_tolerance::Float64
    safety_factor::Float64
    minimum_step_factor::Float64
    maximum_step_factor::Float64
    deterministic::Bool

    function DASSLClassEMTSettings(;
        initial_step_s::Real,
        minimum_step_s::Real,
        maximum_step_s::Real,
        maximum_order::Integer = 5,
        relative_tolerance::Real = 1.0e-7,
        maximum_newton_iterations::Integer = 12,
        maximum_rejected_steps::Integer = 32,
        root_time_tolerance_s::Real = 1.0e-10,
        residual_wrms_tolerance::Real = 1.0,
        correction_wrms_tolerance::Real = 1.0,
        safety_factor::Real = 0.9,
        minimum_step_factor::Real = 0.1,
        maximum_step_factor::Real = 5.0,
        deterministic::Bool = true,
    )
        initial_step = Float64(initial_step_s)
        minimum_step = Float64(minimum_step_s)
        maximum_step = Float64(maximum_step_s)
        all(isfinite, (initial_step, minimum_step, maximum_step)) ||
            throw(ArgumentError("DASSL-class steps must be finite"))
        1.0e-12 <= minimum_step <= initial_step <= maximum_step <= 1.0 ||
            throw(ArgumentError(
                "DASSL-class steps must satisfy 1e-12 <= minimum <= initial <= maximum <= 1 second",
            ))
        order = Int(maximum_order)
        1 <= order <= 5 || throw(ArgumentError(
            "DASSL-class maximum BDF order must be from one through five",
        ))
        relative = Float64(relative_tolerance)
        isfinite(relative) && 1.0e-12 <= relative <= 1.0e-2 ||
            throw(ArgumentError(
                "DASSL-class relative tolerance must be from 1e-12 through 1e-2",
            ))
        newton_limit = Int(maximum_newton_iterations)
        1 <= newton_limit <= 1_000 || throw(ArgumentError(
            "DASSL-class Newton iteration limit must be from one through 1000",
        ))
        rejection_limit = Int(maximum_rejected_steps)
        1 <= rejection_limit <= 100_000 || throw(ArgumentError(
            "DASSL-class rejection limit must be from one through 100000",
        ))
        root_tolerance = Float64(root_time_tolerance_s)
        isfinite(root_tolerance) && 1.0e-12 <= root_tolerance <= maximum_step ||
            throw(ArgumentError(
                "DASSL-class root-time tolerance must be finite and within the step domain",
            ))
        residual_tolerance = Float64(residual_wrms_tolerance)
        correction_tolerance = Float64(correction_wrms_tolerance)
        all(value -> isfinite(value) && 0.0 < value <= 1.0, (
            residual_tolerance,
            correction_tolerance,
        )) || throw(ArgumentError(
            "DASSL-class residual and correction WRMS tolerances must be in (0, 1]",
        ))
        safety = Float64(safety_factor)
        shrink = Float64(minimum_step_factor)
        growth = Float64(maximum_step_factor)
        isfinite(safety) && 0.0 < safety < 1.0 || throw(ArgumentError(
            "DASSL-class safety factor must be finite and in (0, 1)",
        ))
        isfinite(shrink) && 0.0 < shrink < 1.0 || throw(ArgumentError(
            "DASSL-class minimum step factor must be finite and in (0, 1)",
        ))
        isfinite(growth) && 1.0 < growth <= 100.0 || throw(ArgumentError(
            "DASSL-class maximum step factor must be finite and in (1, 100]",
        ))
        deterministic || throw(ArgumentError(
            "the admitted DASSL-class mode currently requires deterministic execution",
        ))
        return new(
            initial_step,
            minimum_step,
            maximum_step,
            order,
            relative,
            newton_limit,
            rejection_limit,
            root_tolerance,
            residual_tolerance,
            correction_tolerance,
            safety,
            shrink,
            growth,
            deterministic,
        )
    end
end

EMTIntegrationSelection(settings::DASSLClassEMTSettings) =
    EMTIntegrationSelection(DASSLClassVariableStep, settings)

"""Exact no-mutation task callback for a portable DASSL-class checkpoint boundary."""
function dassl_class_checkpoint_boundary!(_state, _derivative, _time_s)
    return nothing
end

"""Production-owner network request for the explicitly selected DASSL-class mode."""
struct DASSLClassEMTNetworkRequest
    identity::String
    settings::DASSLClassEMTSettings
    node_count::Int
    initial_time_s::Float64
    stop_time_s::Float64
    owner_identities::Vector{String}
    owners::Vector{Any}
    initial_node_voltage_v::Vector{Float64}
    voltage_absolute_tolerance_v::Float64
    current_absolute_tolerance_a::Float64
    flux_absolute_tolerance_wb::Float64
    angle_absolute_tolerance_rad::Float64
    speed_absolute_tolerance_rad_s::Float64
    voltage_residual_scale_v::Float64
    current_residual_scale_a::Float64
    flux_residual_scale_v::Float64
    angle_residual_scale_rad_s::Float64
    speed_residual_scale_rad_s2::Float64
    requested_output_times_s::Vector{Float64}
    directed_roots::Tuple
    exact_tasks::Tuple

    function DASSLClassEMTNetworkRequest(
        identity::AbstractString,
        settings::DASSLClassEMTSettings,
        node_count::Integer,
        initial_time_s::Real,
        stop_time_s::Real,
        owner_identities,
        owners;
        initial_node_voltage_v,
        voltage_absolute_tolerance_v::Real,
        current_absolute_tolerance_a::Real,
        flux_absolute_tolerance_wb::Real=1.0e-10,
        angle_absolute_tolerance_rad::Real=1.0e-10,
        speed_absolute_tolerance_rad_s::Real=1.0e-8,
        voltage_residual_scale_v::Real,
        current_residual_scale_a::Real,
        flux_residual_scale_v::Real=1.0e-8,
        angle_residual_scale_rad_s::Real=1.0e-8,
        speed_residual_scale_rad_s2::Real=1.0e-7,
        requested_output_times_s=Float64[],
        directed_roots=(),
        exact_tasks=(),
    )
        normalized_identity = String(identity)
        isempty(strip(normalized_identity)) && throw(ArgumentError(
            "DASSL-class network identity must not be empty",
        ))
        nodes = Int(node_count)
        1 <= nodes <= 100_000 || throw(ArgumentError(
            "DASSL-class network node count must be from one through 100000",
        ))
        initial_time = Float64(initial_time_s)
        stop_time = Float64(stop_time_s)
        isfinite(initial_time) && isfinite(stop_time) && stop_time > initial_time ||
            throw(ArgumentError(
                "DASSL-class network stop time must be finite and after initial time",
            ))
        names = String[String(owner) for owner in owner_identities]
        owner_vector = Any[owner for owner in owners]
        length(names) == length(owner_vector) || throw(DimensionMismatch(
            "DASSL-class network owner identities and owners differ in length",
        ))
        isempty(names) && throw(ArgumentError(
            "DASSL-class network must contain at least one physical owner",
        ))
        any(name -> isempty(strip(name)), names) && throw(ArgumentError(
            "DASSL-class network owner identities must not be empty",
        ))
        length(names) == length(unique(names)) || throw(ArgumentError(
            "DASSL-class network owner identities must be unique",
        ))
        node_voltage = Float64[Float64(value) for value in initial_node_voltage_v]
        length(node_voltage) == nodes || throw(DimensionMismatch(
            "DASSL-class initial node-voltage count differs from node count",
        ))
        all(isfinite, node_voltage) || throw(ArgumentError(
            "DASSL-class initial node voltages must be finite",
        ))
        numerical_scales = Float64[
            Float64(voltage_absolute_tolerance_v),
            Float64(current_absolute_tolerance_a),
            Float64(flux_absolute_tolerance_wb),
            Float64(angle_absolute_tolerance_rad),
            Float64(speed_absolute_tolerance_rad_s),
            Float64(voltage_residual_scale_v),
            Float64(current_residual_scale_a),
            Float64(flux_residual_scale_v),
            Float64(angle_residual_scale_rad_s),
            Float64(speed_residual_scale_rad_s2),
        ]
        all(value -> isfinite(value) && value > 0.0, numerical_scales) ||
            throw(ArgumentError(
                "DASSL-class physical tolerance and residual scales must be finite and positive",
            ))
        output_times = Float64[Float64(value) for value in requested_output_times_s]
        length(output_times) <= 10_000 || throw(ArgumentError(
            "DASSL-class network output grid exceeds 10000 points",
        ))
        all(isfinite, output_times) && issorted(output_times) &&
            length(output_times) == length(unique(output_times)) ||
            throw(ArgumentError(
                "DASSL-class network output times must be finite and strictly increasing",
            ))
        all(time -> initial_time <= time <= stop_time, output_times) ||
            throw(ArgumentError(
                "DASSL-class network output times must lie inside the integration interval",
            ))
        roots = Tuple(directed_roots)
        tasks = Tuple(sort!(collect(exact_tasks); by=task -> (
            task.time_s,
            task.priority,
            String(task.identity),
        )))
        all(root -> root isa DASSLClassEMTValidationRoot, roots) ||
            throw(ArgumentError(
                "DASSL-class initial network roots must use the public directed-root contract",
            ))
        length(roots) == length(unique(root.identity for root in roots)) ||
            throw(ArgumentError(
                "DASSL-class initial network root identities must be unique",
            ))
        all(task -> task isa DASSLClassEMTValidationTask, tasks) ||
            throw(ArgumentError(
                "DASSL-class initial network tasks must use the public exact-task contract",
            ))
        all(task -> initial_time < task.time_s <= stop_time, tasks) ||
            throw(ArgumentError(
                "DASSL-class initial network tasks must lie after start and within the interval",
            ))
        task_keys = [(task.time_s, task.identity) for task in tasks]
        length(task_keys) == length(unique(task_keys)) || throw(ArgumentError(
            "DASSL-class initial network repeats one task identity at one time",
        ))
        return new(
            normalized_identity,
            settings,
            nodes,
            initial_time,
            stop_time,
            names,
            owner_vector,
            node_voltage,
            numerical_scales...,
            output_times,
            roots,
            tasks,
        )
    end
end

"""Stable identity, differential mask, SI unit, and absolute tolerance for one DAE state vector."""
struct DASSLClassEMTStateLayout
    identities::Vector{Symbol}
    differential_mask::BitVector
    units::Vector{String}
    absolute_tolerances::Vector{Float64}

    function DASSLClassEMTStateLayout(
        identities,
        differential_mask,
        units,
        absolute_tolerances,
    )
        names = Symbol[Symbol(identity) for identity in identities]
        isempty(names) && throw(ArgumentError(
            "DASSL-class state layout must contain at least one state",
        ))
        any(==(Symbol("")), names) && throw(ArgumentError(
            "DASSL-class state identity must not be empty",
        ))
        length(names) == length(unique(names)) || throw(ArgumentError(
            "DASSL-class state identities must be unique",
        ))
        mask = BitVector(differential_mask)
        typed_units = String[String(unit) for unit in units]
        tolerances = Float64[Float64(value) for value in absolute_tolerances]
        length(mask) == length(names) || throw(DimensionMismatch(
            "DASSL-class differential mask length differs from state layout",
        ))
        length(typed_units) == length(names) || throw(DimensionMismatch(
            "DASSL-class unit count differs from state layout",
        ))
        length(tolerances) == length(names) || throw(DimensionMismatch(
            "DASSL-class absolute-tolerance count differs from state layout",
        ))
        any(unit -> isempty(strip(unit)), typed_units) && throw(ArgumentError(
            "DASSL-class state units must be explicit",
        ))
        all(value -> isfinite(value) && value > 0.0, tolerances) ||
            throw(ArgumentError(
                "DASSL-class absolute tolerances must be finite and positive",
            ))
        return new(names, mask, typed_units, tolerances)
    end
end

"""Validation-only directed root and atomic reset contract."""
struct DASSLClassEMTValidationRoot{F,A}
    identity::Symbol
    direction::Int8
    priority::Int
    root::F
    reset!::A

    function DASSLClassEMTValidationRoot(
        identity::Symbol,
        direction::Integer,
        priority::Integer,
        root,
        reset!,
    )
        identity == Symbol("") && throw(ArgumentError(
            "DASSL-class validation-root identity must not be empty",
        ))
        direction in (-1, 0, 1) || throw(ArgumentError(
            "DASSL-class validation-root direction must be -1, 0, or 1",
        ))
        root isa Function || throw(ArgumentError(
            "DASSL-class validation-root evaluator must be callable",
        ))
        reset! isa Function || throw(ArgumentError(
            "DASSL-class validation-root reset must be callable",
        ))
        return new{typeof(root),typeof(reset!)}(
            identity,
            Int8(direction),
            Int(priority),
            root,
            reset!,
        )
    end
end

"""Validation-only exact scheduled boundary and atomic task action."""
struct DASSLClassEMTValidationTask{A}
    identity::Symbol
    time_s::Float64
    priority::Int
    apply!::A

    function DASSLClassEMTValidationTask(
        identity::Symbol,
        time_s::Real,
        priority::Integer,
        apply!,
    )
        identity == Symbol("") && throw(ArgumentError(
            "DASSL-class validation-task identity must not be empty",
        ))
        typed_time = Float64(time_s)
        isfinite(typed_time) || throw(ArgumentError(
            "DASSL-class validation-task time must be finite",
        ))
        apply! isa Function || throw(ArgumentError(
            "DASSL-class validation-task action must be callable",
        ))
        return new{typeof(apply!)}(
            identity,
            typed_time,
            Int(priority),
            apply!,
        )
    end
end

"""Validation-only residual system; it is never a production physical-owner admission."""
struct DASSLClassEMTValidationProblem{R,J}
    identity::String
    settings::DASSLClassEMTSettings
    layout::DASSLClassEMTStateLayout
    initial_time_s::Float64
    stop_time_s::Float64
    initial_state::Vector{Float64}
    initial_derivative::Vector{Float64}
    residual_scales::Vector{Float64}
    requested_output_times_s::Vector{Float64}
    directed_roots::Tuple
    exact_tasks::Tuple
    residual!::R
    jacobian!::J

    function DASSLClassEMTValidationProblem(
        identity::AbstractString,
        settings::DASSLClassEMTSettings,
        layout::DASSLClassEMTStateLayout,
        initial_time_s::Real,
        stop_time_s::Real,
        initial_state,
        initial_derivative,
        residual_scales,
        residual!,
        jacobian!,
        ;
        requested_output_times_s=Float64[],
        directed_roots=(),
        exact_tasks=(),
    )
        normalized_identity = String(identity)
        isempty(strip(normalized_identity)) && throw(ArgumentError(
            "DASSL-class validation problem identity must not be empty",
        ))
        initial_time = Float64(initial_time_s)
        stop_time = Float64(stop_time_s)
        isfinite(initial_time) && isfinite(stop_time) && stop_time > initial_time ||
            throw(ArgumentError(
                "DASSL-class validation stop time must be finite and after initial time",
            ))
        state = Float64[Float64(value) for value in initial_state]
        derivative = Float64[Float64(value) for value in initial_derivative]
        scales = Float64[Float64(value) for value in residual_scales]
        output_times = Float64[Float64(value) for value in requested_output_times_s]
        roots = Tuple(directed_roots)
        tasks = Tuple(sort!(collect(exact_tasks); by=task -> (
            task.time_s,
            task.priority,
            String(task.identity),
        )))
        dimension = length(layout.identities)
        all(values -> length(values) == dimension, (state, derivative, scales)) ||
            throw(DimensionMismatch(
                "DASSL-class validation state, derivative, and residual scales must match the layout",
            ))
        all(isfinite, state) && all(isfinite, derivative) || throw(ArgumentError(
            "DASSL-class validation initial values must be finite",
        ))
        all(value -> isfinite(value) && value > 0.0, scales) ||
            throw(ArgumentError(
                "DASSL-class validation residual scales must be finite and positive",
            ))
        length(output_times) <= 10_000 || throw(ArgumentError(
            "DASSL-class validation output grid exceeds 10000 points",
        ))
        all(isfinite, output_times) || throw(ArgumentError(
            "DASSL-class validation output times must be finite",
        ))
        issorted(output_times) && length(output_times) == length(unique(output_times)) ||
            throw(ArgumentError(
                "DASSL-class validation output times must be strictly increasing",
            ))
        all(time -> initial_time <= time <= stop_time, output_times) ||
            throw(ArgumentError(
                "DASSL-class validation output times must lie inside the integration interval",
            ))
        all(root -> root isa DASSLClassEMTValidationRoot, roots) ||
            throw(ArgumentError(
                "DASSL-class validation roots must use validation-root contracts",
            ))
        length(roots) == length(unique(root.identity for root in roots)) ||
            throw(ArgumentError(
                "DASSL-class validation-root identities must be unique",
            ))
        all(task -> task isa DASSLClassEMTValidationTask, tasks) ||
            throw(ArgumentError(
                "DASSL-class validation tasks must use validation-task contracts",
            ))
        all(task -> initial_time < task.time_s <= stop_time, tasks) ||
            throw(ArgumentError(
                "DASSL-class validation tasks must lie after start and within the interval",
            ))
        task_keys = [(task.time_s, task.identity) for task in tasks]
        length(task_keys) == length(unique(task_keys)) || throw(ArgumentError(
            "DASSL-class validation repeats one task identity at one time",
        ))
        any(index -> !layout.differential_mask[index] && !iszero(derivative[index]),
            eachindex(derivative)) && throw(ArgumentError(
            "DASSL-class validation algebraic derivative slots must be zero",
        ))
        residual! isa Function || throw(ArgumentError(
            "DASSL-class validation residual owner must be callable",
        ))
        jacobian! isa Function || throw(ArgumentError(
            "DASSL-class validation Jacobian owner must be callable",
        ))
        return new{typeof(residual!),typeof(jacobian!)}(
            normalized_identity,
            settings,
            layout,
            initial_time,
            stop_time,
            state,
            derivative,
            scales,
            output_times,
            roots,
            tasks,
            residual!,
            jacobian!,
        )
    end
end

"""One deterministic admitted or refused physical-owner readiness record."""
struct DASSLClassEMTOwnerDisposition
    identity::String
    owner_type::String
    admitted::Bool
    code::Symbol
    message::String

    function DASSLClassEMTOwnerDisposition(
        identity::AbstractString,
        owner_type::AbstractString,
        admitted::Bool,
        code::Symbol,
        message::AbstractString,
    )
        normalized_identity = String(identity)
        normalized_type = String(owner_type)
        normalized_message = String(message)
        isempty(strip(normalized_identity)) && throw(ArgumentError(
            "DASSL-class owner identity must not be empty",
        ))
        isempty(strip(normalized_type)) && throw(ArgumentError(
            "DASSL-class owner type must not be empty",
        ))
        code == Symbol("") && throw(ArgumentError(
            "DASSL-class owner disposition code must not be empty",
        ))
        isempty(strip(normalized_message)) && throw(ArgumentError(
            "DASSL-class owner disposition message must not be empty",
        ))
        return new(
            normalized_identity,
            normalized_type,
            admitted,
            code,
            normalized_message,
        )
    end
end

"""Complete sorted owner inventory returned before private runtime allocation."""
struct DASSLClassEMTReadiness
    compatible::Bool
    owners::Vector{DASSLClassEMTOwnerDisposition}
    state_count::Int
    differential_state_count::Int
    algebraic_state_count::Int
    signature::String

    function DASSLClassEMTReadiness(
        owners,
        state_count::Integer,
        differential_state_count::Integer,
        algebraic_state_count::Integer,
        signature::AbstractString,
    )
        ordered = sort!(DASSLClassEMTOwnerDisposition[owners...]; by=owner -> (
            owner.identity,
            owner.owner_type,
            String(owner.code),
        ))
        keys = [(owner.identity, owner.owner_type) for owner in ordered]
        length(keys) == length(unique(keys)) || throw(ArgumentError(
            "DASSL-class readiness repeats a physical owner",
        ))
        states = Int(state_count)
        differential = Int(differential_state_count)
        algebraic = Int(algebraic_state_count)
        states >= 0 && differential >= 0 && algebraic >= 0 ||
            throw(ArgumentError("DASSL-class readiness counts must be nonnegative"))
        differential + algebraic == states || throw(ArgumentError(
            "DASSL-class differential and algebraic counts must equal total state count",
        ))
        normalized_signature = String(signature)
        isempty(strip(normalized_signature)) && throw(ArgumentError(
            "DASSL-class readiness signature must not be empty",
        ))
        return new(
            all(owner -> owner.admitted, ordered),
            ordered,
            states,
            differential,
            algebraic,
            normalized_signature,
        )
    end
end

"""Typed atomic failure from the last unchanged accepted DAE state."""
struct DASSLClassEMTFailure <: Exception
    code::Symbol
    owner::Union{Nothing,String}
    time_s::Float64
    attempted_step_s::Float64
    order::Int
    residual_wrms::Float64
    error_wrms::Float64
    last_accepted_time_s::Float64
    message::String
end

function Base.showerror(io::IO, failure::DASSLClassEMTFailure)
    print(io, String(failure.code), ": ", failure.message)
    isnothing(failure.owner) || print(io, " [owner=", failure.owner, ']')
    print(io, " [last_accepted_time_s=", failure.last_accepted_time_s, ']')
end

"""Deterministic numerical work summary; physical results remain separately typed."""
struct DASSLClassEMTDiagnostics
    accepted_steps::Int
    rejected_steps::Int
    newton_iterations::Int
    residual_evaluations::Int
    jacobian_evaluations::Int
    factorization_count::Int
    factorization_reuse_count::Int
    root_evaluations::Int
    localized_roots::Int
    consistency_restarts::Int
    minimum_accepted_step_s::Float64
    maximum_accepted_step_s::Float64
    order_counts::NTuple{5,Int}
end

"""Public variable-step result at requested and explicitly labelled transition times."""
struct DASSLClassEMTResult
    readiness::DASSLClassEMTReadiness
    layout::DASSLClassEMTStateLayout
    time_s::Vector{Float64}
    state::Matrix{Float64}
    derivative::Matrix{Float64}
    transition_side::Vector{Symbol}
    diagnostics::DASSLClassEMTDiagnostics
    numerical_snapshot_signature::String
    deterministic_signature::String

    function DASSLClassEMTResult(
        readiness::DASSLClassEMTReadiness,
        layout::DASSLClassEMTStateLayout,
        time_s,
        state,
        derivative,
        transition_side,
        diagnostics::DASSLClassEMTDiagnostics,
        numerical_snapshot_signature::AbstractString,
        deterministic_signature::AbstractString,
    )
        times = Float64[Float64(value) for value in time_s]
        states = Matrix{Float64}(state)
        derivatives = Matrix{Float64}(derivative)
        sides = Symbol[Symbol(side) for side in transition_side]
        size(states) == (length(layout.identities), length(times)) ||
            throw(DimensionMismatch("DASSL-class result state shape is invalid"))
        size(derivatives) == size(states) || throw(DimensionMismatch(
            "DASSL-class result derivative shape differs from state",
        ))
        length(sides) == length(times) || throw(DimensionMismatch(
            "DASSL-class transition labels differ from output times",
        ))
        issorted(times) || throw(ArgumentError(
            "DASSL-class result times must be nondecreasing",
        ))
        all(isfinite, times) && all(isfinite, states) && all(isfinite, derivatives) ||
            throw(ArgumentError("DASSL-class result contains nonfinite values"))
        all(side -> side in (:smooth, :left, :right), sides) || throw(ArgumentError(
            "DASSL-class transition labels must be smooth, left, or right",
        ))
        snapshot_signature = String(numerical_snapshot_signature)
        isempty(strip(snapshot_signature)) && throw(ArgumentError(
            "DASSL-class numerical snapshot signature must not be empty",
        ))
        signature = String(deterministic_signature)
        isempty(strip(signature)) && throw(ArgumentError(
            "DASSL-class deterministic signature must not be empty",
        ))
        return new(
            readiness,
            layout,
            times,
            states,
            derivatives,
            sides,
            diagnostics,
            snapshot_signature,
            signature,
        )
    end
end

prepare_dassl_class_emt(::AbstractAIMORASolverBackend, request) =
    _solver_unavailable_result(
        :prepare_dassl_class_emt,
        :dassl_class_variable_step_emt;
        message = "The active backend does not implement DASSL-class EMT preparation.",
    )

execute_dassl_class_emt!(::AbstractAIMORASolverBackend, prepared; kwargs...) =
    _solver_unavailable_result(
        :execute_dassl_class_emt,
        :dassl_class_variable_step_emt;
        message = "The active backend does not implement DASSL-class EMT execution.",
    )

dassl_class_emt_readiness(::AbstractAIMORASolverBackend, request) =
    _solver_unavailable_result(
        :dassl_class_emt_readiness,
        :dassl_class_variable_step_emt;
        message = "The active backend does not implement DASSL-class owner readiness.",
    )

function prepare_dassl_class_emt(request)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :prepare_dassl_class_emt,
        :dassl_class_variable_step_emt,
    ) : prepare_dassl_class_emt(backend, request)
end

function execute_dassl_class_emt!(prepared; kwargs...)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :execute_dassl_class_emt,
        :dassl_class_variable_step_emt,
    ) : execute_dassl_class_emt!(backend, prepared; kwargs...)
end

function dassl_class_emt_readiness(request)
    backend = active_solver_backend()
    return backend === nothing ?
           _solver_unavailable_result(
        :dassl_class_emt_readiness,
        :dassl_class_variable_step_emt,
    ) : dassl_class_emt_readiness(backend, request)
end

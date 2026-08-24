export SurgePropagationSection,
       TransmissionTowerModel,
       GISGILSection,
       traveling_wave_components,
       reconstruct_terminal_state,
       traveling_wave_reflection,
       tower_total_travel_time_s,
       gis_gil_series_impedance_ohm,
       gis_gil_shunt_admittance_s

function _symmetric_positive_definite_matrix(
    values::AbstractMatrix{<:Real},
    name::AbstractString,
)
    rows, columns = size(values)
    rows == columns && rows > 0 || throw(ArgumentError("$name must be nonempty and square"))
    matrix = Matrix{Float64}(values)
    all(isfinite, matrix) || throw(ArgumentError("$name must contain finite values"))
    scale = max(opnorm(matrix, Inf), 1.0)
    maximum(abs, matrix - transpose(matrix)) <= 128.0 * eps(Float64) * scale ||
        throw(ArgumentError("$name must be symmetric"))
    matrix = 0.5 .* (matrix .+ transpose(matrix))
    isposdef(Symmetric(matrix)) || throw(ArgumentError("$name must be positive definite"))
    return matrix
end

function _symmetric_positive_semidefinite_matrix(
    values::AbstractMatrix{<:Real},
    name::AbstractString,
)
    rows, columns = size(values)
    rows == columns && rows > 0 || throw(ArgumentError("$name must be nonempty and square"))
    matrix = Matrix{Float64}(values)
    all(isfinite, matrix) || throw(ArgumentError("$name must contain finite values"))
    scale = max(opnorm(matrix, Inf), 1.0)
    maximum(abs, matrix - transpose(matrix)) <= 128.0 * eps(Float64) * scale ||
        throw(ArgumentError("$name must be symmetric"))
    matrix = 0.5 .* (matrix .+ transpose(matrix))
    minimum(eigvals(Symmetric(matrix))) >= -256.0 * eps(Float64) * scale ||
        throw(ArgumentError("$name must be positive semidefinite"))
    return matrix
end

"""One oriented positive-real multiconductor surge propagation section."""
struct SurgePropagationSection
    identity::Symbol
    characteristic_impedance_ohm::Matrix{Float64}
    characteristic_admittance_s::Matrix{Float64}
    length_m::Float64
    propagation_speed_m_per_s::Float64
    travel_time_s::Float64
    attenuation_neper::Float64
    provenance::SurgeParameterProvenance

    function SurgePropagationSection(
        identity::Symbol,
        characteristic_impedance_ohm::AbstractMatrix{<:Real};
        length_m::Real,
        propagation_speed_m_per_s::Real,
        attenuation_neper::Real=0.0,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "ohm, metre, metre per second, and neper",
            "positive-real oriented surge propagation section",
        ),
    )
        identity == Symbol("") && throw(ArgumentError("surge section identity must not be empty"))
        impedance = _symmetric_positive_definite_matrix(
            characteristic_impedance_ohm,
            "surge characteristic impedance",
        )
        admittance = inv(impedance)
        length = _positive_finite(length_m, "surge section length")
        speed = _positive_finite(propagation_speed_m_per_s, "surge propagation speed")
        attenuation = _nonnegative_finite(attenuation_neper, "surge attenuation")
        _require_physical_provenance(provenance, "surge propagation section")
        return new(
            identity,
            impedance,
            admittance,
            length,
            speed,
            length / speed,
            attenuation,
            provenance,
        )
    end
end

function traveling_wave_components(
    section::SurgePropagationSection,
    voltage_v::AbstractVector{<:Real},
    current_a::AbstractVector{<:Real},
)
    phase_count = size(section.characteristic_impedance_ohm, 1)
    length(voltage_v) == length(current_a) == phase_count || throw(DimensionMismatch(
        "surge voltage/current vectors must match section conductor count",
    ))
    voltage = Float64.(voltage_v)
    current = Float64.(current_a)
    all(isfinite, voltage) && all(isfinite, current) || throw(ArgumentError(
        "surge terminal state must be finite",
    ))
    impedance_current = section.characteristic_impedance_ohm * current
    return (
        forward_voltage_v=0.5 .* (voltage .+ impedance_current),
        reverse_voltage_v=0.5 .* (voltage .- impedance_current),
    )
end

function reconstruct_terminal_state(
    section::SurgePropagationSection,
    forward_voltage_v::AbstractVector{<:Real},
    reverse_voltage_v::AbstractVector{<:Real},
)
    phase_count = size(section.characteristic_impedance_ohm, 1)
    length(forward_voltage_v) == length(reverse_voltage_v) == phase_count ||
        throw(DimensionMismatch("surge wave vectors must match section conductor count"))
    forward = Float64.(forward_voltage_v)
    reverse = Float64.(reverse_voltage_v)
    return (
        voltage_v=forward .+ reverse,
        current_a=section.characteristic_admittance_s * (forward .- reverse),
    )
end

function traveling_wave_reflection(
    incident_section::SurgePropagationSection,
    terminating_impedance_ohm::AbstractMatrix{<:Real},
    incident_voltage_v::AbstractVector{<:Real},
)
    termination = _symmetric_positive_definite_matrix(
        terminating_impedance_ohm,
        "surge terminating impedance",
    )
    size(termination) == size(incident_section.characteristic_impedance_ohm) ||
        throw(DimensionMismatch("terminating impedance size must match surge section"))
    incident = Float64.(incident_voltage_v)
    length(incident) == size(termination, 1) || throw(DimensionMismatch(
        "incident wave size must match surge section",
    ))
    characteristic = incident_section.characteristic_impedance_ohm
    reflection_operator = (termination - characteristic) /
        (termination + characteristic)
    reflected = reflection_operator * incident
    transmitted = incident + reflected
    return (
        reflected_voltage_v=reflected,
        transmitted_voltage_v=transmitted,
        reflection_operator=reflection_operator,
    )
end

"""Ordered multisection tower or lightning-channel surge path."""
struct TransmissionTowerModel
    identity::Symbol
    sections::Vector{SurgePropagationSection}
    attachment_nodes::Vector{Symbol}
    total_travel_time_s::Float64
    provenance::SurgeParameterProvenance

    function TransmissionTowerModel(
        identity::Symbol,
        sections::AbstractVector{SurgePropagationSection},
        attachment_nodes::AbstractVector{Symbol};
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "section-specific SI surge parameters",
            "oriented multisection tower or channel path",
        ),
    )
        identity == Symbol("") && throw(ArgumentError("tower identity must not be empty"))
        checked_sections = collect(sections)
        isempty(checked_sections) && throw(ArgumentError("tower must contain a surge section"))
        conductor_count = size(first(checked_sections).characteristic_impedance_ohm, 1)
        all(
            section -> size(section.characteristic_impedance_ohm, 1) == conductor_count,
            checked_sections,
        ) || throw(ArgumentError("tower sections must share one conductor count"))
        nodes = Symbol.(attachment_nodes)
        length(nodes) == length(checked_sections) + 1 || throw(ArgumentError(
            "tower attachment nodes must delimit every section",
        ))
        length(unique(nodes)) == length(nodes) || throw(ArgumentError(
            "tower attachment node identities must be unique",
        ))
        _require_physical_provenance(provenance, "transmission tower")
        return new(
            identity,
            checked_sections,
            nodes,
            sum(section.travel_time_s for section in checked_sections),
            provenance,
        )
    end
end

tower_total_travel_time_s(tower::TransmissionTowerModel) = tower.total_travel_time_s

"""Passive multiconductor gas-insulated line terminal section."""
struct GISGILSection
    identity::Symbol
    resistance_ohm_per_m::Matrix{Float64}
    inductance_h_per_m::Matrix{Float64}
    conductance_s_per_m::Matrix{Float64}
    capacitance_f_per_m::Matrix{Float64}
    length_m::Float64
    conductor_names::Vector{Symbol}
    enclosure_reference::Symbol
    provenance::SurgeParameterProvenance

    function GISGILSection(
        identity::Symbol,
        resistance_ohm_per_m::AbstractMatrix{<:Real},
        inductance_h_per_m::AbstractMatrix{<:Real},
        conductance_s_per_m::AbstractMatrix{<:Real},
        capacitance_f_per_m::AbstractMatrix{<:Real};
        length_m::Real,
        conductor_names::AbstractVector{Symbol},
        enclosure_reference::Symbol,
        provenance::SurgeParameterProvenance=_surge_default_provenance(
            "ohm/metre, henry/metre, siemens/metre, farad/metre, and metre",
            "passive multiconductor GIS/GIL terminal section",
        ),
    )
        identity == Symbol("") && throw(ArgumentError("GIS/GIL identity must not be empty"))
        resistance = _symmetric_positive_semidefinite_matrix(
            resistance_ohm_per_m,
            "GIS/GIL resistance",
        )
        inductance = _symmetric_positive_definite_matrix(
            inductance_h_per_m,
            "GIS/GIL inductance",
        )
        conductance = _symmetric_positive_semidefinite_matrix(
            conductance_s_per_m,
            "GIS/GIL conductance",
        )
        capacitance = _symmetric_positive_definite_matrix(
            capacitance_f_per_m,
            "GIS/GIL capacitance",
        )
        size(resistance) == size(inductance) == size(conductance) == size(capacitance) ||
            throw(DimensionMismatch("GIS/GIL RLCG matrices must share one size"))
        names = Symbol.(conductor_names)
        length(names) == size(resistance, 1) || throw(DimensionMismatch(
            "GIS/GIL conductor names must match matrix size",
        ))
        length(unique(names)) == length(names) || throw(ArgumentError(
            "GIS/GIL conductor names must be unique",
        ))
        enclosure_reference == Symbol("") && throw(ArgumentError(
            "GIS/GIL enclosure reference must not be empty",
        ))
        section_length = _positive_finite(length_m, "GIS/GIL length")
        _require_physical_provenance(provenance, "GIS/GIL section")
        return new(
            identity,
            resistance,
            inductance,
            conductance,
            capacitance,
            section_length,
            names,
            enclosure_reference,
            provenance,
        )
    end
end

function gis_gil_series_impedance_ohm(section::GISGILSection, frequency_hz::Real)
    frequency = _nonnegative_finite(frequency_hz, "GIS/GIL frequency")
    return section.length_m .* (
        section.resistance_ohm_per_m .+
        complex(0.0, 2.0 * pi * frequency) .* section.inductance_h_per_m
    )
end

function gis_gil_shunt_admittance_s(section::GISGILSection, frequency_hz::Real)
    frequency = _nonnegative_finite(frequency_hz, "GIS/GIL frequency")
    return section.length_m .* (
        section.conductance_s_per_m .+
        complex(0.0, 2.0 * pi * frequency) .* section.capacitance_f_per_m
    )
end

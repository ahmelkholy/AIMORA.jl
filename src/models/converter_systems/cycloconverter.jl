export CycloconverterInitialState,
       SwitchingCycloconverterStudy,
       SwitchingCycloconverterTrace,
       detailed_cycloconverter_semiconductor_signatures

struct CycloconverterInitialState
    output_current_a::Vector{Float64}
    active_bridge_group::Vector{Int8}

    function CycloconverterInitialState(
        output_current_a=(0.0,);
        active_bridge_group=nothing,
    )
        current = collect(Float64.(output_current_a))
        length(current) in (1, 3) && all(isfinite, current) || throw(ArgumentError(
            "cycloconverter initial state requires one or three finite output currents",
        ))
        length(current) == 3 &&
            abs(sum(current)) > 1.0e-10 * max(maximum(abs, current), 1.0) &&
            throw(ArgumentError(
                "three-wire cycloconverter initial output currents must sum to zero",
            ))
        groups = active_bridge_group === nothing ?
            fill(Int8(0), length(current)) :
            collect(Int8.(active_bridge_group))
        length(groups) == length(current) && all(value -> value in (-1, 0, 1), groups) ||
            throw(ArgumentError(
                "cycloconverter initial bridge groups must be negative, blocked, or positive",
            ))
        all(index -> iszero(current[index]) ||
            groups[index] == Int8(sign(current[index])), eachindex(current)) ||
            throw(ArgumentError(
                "energized cycloconverter output current must match its active bridge group",
            ))
        return new(current, groups)
    end
end

function detailed_cycloconverter_semiconductor_signatures(
    parameters::DetailedConverterSemiconductorParameters,
    valve_count::Integer,
)
    count = Int(valve_count)
    count in (12, 36) || throw(ArgumentError(
        "cycloconverter detailed fidelity requires 12 or 36 thyristor positions",
    ))
    body = repr(parameters)
    return ntuple(count) do position
        bytes2hex(sha256("line_commutated_cycloconverter_valve_$(position)\n" * body))
    end
end

struct SwitchingCycloconverterStudy{T<:Tuple}
    specification::ConverterSystemSpecification
    topologies::T
    input_phase_voltage_peak_v::Float64
    input_frequency_hz::Float64
    source_resistance_ohm::Float64
    load_resistance_ohm::Float64
    load_inductance_h::Float64
    initial_state::CycloconverterInitialState
    detailed_semiconductor::Union{Nothing,DetailedConverterSemiconductorParameters}
    circulating_current::Bool
    current_zero_tolerance_a::Float64
    start_time_s::Float64
    stop_time_s::Float64

    function SwitchingCycloconverterStudy(
        specification::ConverterSystemSpecification;
        topology,
        input_phase_voltage_peak_v::Real,
        input_frequency_hz::Real,
        source_resistance_ohm::Real,
        load_resistance_ohm::Real,
        load_inductance_h::Real,
        initial_state::CycloconverterInitialState=CycloconverterInitialState(
            zeros(specification.selection.phase_count),
        ),
        detailed_semiconductor::Union{
            Nothing,
            DetailedConverterSemiconductorParameters,
        }=nothing,
        circulating_current::Bool=false,
        current_zero_tolerance_a::Real=1.0e-6,
        start_time_s::Real=0.0,
        stop_time_s::Real,
    )
        selection = specification.selection
        selection.family === LineCommutatedCycloconverter || throw(ArgumentError(
            "cycloconverter study requires the canonical line-commutated family",
        ))
        selection.fidelity in (SwitchingStateEquivalent, SwitchingDetailed) ||
            throw(ArgumentError(
                "cycloconverter execution requires switching-state or switching-detailed fidelity",
            ))
        selection.application === StandaloneConversion || throw(ArgumentError(
            "cycloconverter study is a standalone conversion owner",
        ))
        selection.phase_count in (1, 3) &&
            length(initial_state.output_current_a) == selection.phase_count ||
            throw(ArgumentError(
                "cycloconverter output and initial-state phase counts do not agree",
            ))
        topologies = topology isa BridgeTopologyDescriptor ? (topology,) : Tuple(topology)
        length(topologies) == selection.phase_count && all(topology ->
            topology isa BridgeTopologyDescriptor &&
            topology.family === :cycloconverter &&
            length(topology.valve_positions) == 12 &&
            length(topology.state_groups) == 1 &&
            size(only(topology.state_groups).admitted_states) == (12, 13) &&
            all(==( :thyristor), getfield.(topology.valve_positions, :valve_class)),
            topologies) ||
            throw(ArgumentError(
                "cycloconverter study requires one canonical noncirculating one-output topology per output phase",
            ))
        circulating_current && throw(ArgumentError(
            "circulating-current cycloconversion requires an explicit intergroup reactor and is not released by this noncirculating study",
        ))
        specification.topology_signatures ==
            Tuple(bridge_topology_signature.(topologies)) ||
            throw(ArgumentError(
                "cycloconverter specification does not bind its exact topology",
            ))
        specification.modulation.kind === CycloconverterFiringSynthesis ||
            throw(ArgumentError(
                "cycloconverter execution requires cycloconverter firing synthesis",
            ))
        0.0 < specification.modulation.modulation_index <= 1.0 ||
            throw(ArgumentError(
                "cycloconverter modulation index must lie in (0, 1]",
            ))
        detailed = selection.fidelity === SwitchingDetailed
        valve_count = 12 * selection.phase_count
        if detailed
            detailed_semiconductor === nothing && throw(ArgumentError(
                "switching-detailed cycloconversion requires typed D200 parameters",
            ))
            specification.device_fidelity_signatures ==
                detailed_cycloconverter_semiconductor_signatures(
                    detailed_semiconductor,
                    valve_count,
                ) || throw(ArgumentError(
                    "cycloconverter specification does not bind every detailed valve",
                ))
        else
            detailed_semiconductor === nothing || throw(ArgumentError(
                "switching-state cycloconversion cannot accept detailed semiconductor data",
            ))
            isempty(specification.device_fidelity_signatures) || throw(ArgumentError(
                "switching-state cycloconversion cannot claim D200 identities",
            ))
        end
        values = Float64.((input_phase_voltage_peak_v, input_frequency_hz,
            source_resistance_ohm, load_resistance_ohm, load_inductance_h,
            current_zero_tolerance_a, start_time_s, stop_time_s))
        all(isfinite, values) && all(>(0.0), values[1:6]) &&
            values[7] >= 0.0 && values[8] > values[7] || throw(ArgumentError(
                "cycloconverter source, load, tolerance, and horizon domain is invalid",
            ))
        output_frequency = specification.rated_bases.frequency_hz
        0.0 < output_frequency < values[2] && output_frequency / values[2] <= 0.5 ||
            throw(ArgumentError(
                "cycloconverter output/input frequency ratio must lie in (0, 0.5]",
            ))
        isapprox(specification.timing.firing_frequency_hz, values[2];
            atol=0.0, rtol=16eps(Float64)) || throw(ArgumentError(
                "cycloconverter firing calendar must bind the input frequency",
            ))
        step_count = (values[8] - values[7]) / specification.timing.fixed_step_s
        isapprox(step_count, round(step_count); atol=1.0e-10, rtol=1.0e-10) ||
            throw(ArgumentError(
                "cycloconverter horizon must contain integer fixed steps",
            ))
        if detailed
            specification.timing.fixed_step_s <= min(
                detailed_semiconductor.recovered_charge_lifetime_s,
                detailed_semiconductor.turn_off_tail_time_s,
            ) / 10.0 || throw(ArgumentError(
                "cycloconverter timestep must resolve detailed stored-charge state",
            ))
        end
        converter_system_is_ready(converter_system_readiness(specification)) ||
            throw(ArgumentError("cycloconverter specification is not ready"))
        return new{typeof(topologies)}(
            specification,
            topologies,
            values[1:5]...,
            initial_state,
            detailed_semiconductor,
            false,
            values[6],
            values[7:8]...,
        )
    end
end

struct SwitchingCycloconverterTrace
    time_s::Vector{Float64}
    source_voltage_v::Matrix{Float64}
    source_current_a::Matrix{Float64}
    output_reference_voltage_v::Matrix{Float64}
    output_voltage_v::Matrix{Float64}
    output_current_a::Matrix{Float64}
    firing_angle_rad::Matrix{Float64}
    requested_bridge_group::Matrix{Int8}
    active_bridge_group::Matrix{Int8}
    requested_firing_state::BitMatrix
    applied_firing_state::BitMatrix
    conducting_state::BitMatrix
    commutation_overlap_valve_count::Matrix{UInt8}
    failed_commutation::BitMatrix
    reactive_power_var::Vector{Float64}
    stored_energy_j::Vector{Float64}
    semiconductor_loss_w::Vector{Float64}
    kcl_residual_a::Vector{Float64}
    energy_residual_w::Vector{Float64}
    junction_temperature_k::Matrix{Float64}
    recovered_charge_c::Matrix{Float64}
    result::ConverterSystemResult
end

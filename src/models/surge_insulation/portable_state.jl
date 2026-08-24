export SurgeDeviceSnapshot,
       snapshot_surge_device,
       restore_surge_device!

"""Portable, identity-bound scientific state for one surge/insulation owner."""
struct SurgeDeviceSnapshot{S<:NamedTuple}
    schema_version::Int
    model_type::Symbol
    identity_signature::String
    state::S
    checksum::String
end

function _surge_snapshot_checksum(
    schema_version::Int,
    model_type::Symbol,
    identity_signature::String,
    state::NamedTuple,
)
    return bytes2hex(sha256(codeunits(join(
        [
            "aimora-surge-device-snapshot-v1",
            string(schema_version),
            String(model_type),
            identity_signature,
            repr(state),
        ],
        '\n',
    ))))
end

function _surge_identity_signature(model_type::Symbol, identity_payload)
    return bytes2hex(sha256(codeunits(join(
        ["aimora-surge-device-identity-v1", String(model_type), repr(identity_payload)],
        '\n',
    ))))
end

_surge_identity_payload(device::CombinedArcBranch) = (
    nodes=(device.positive_node, device.negative_node),
    cassie_power_w=device.cassie_power_w,
    mayr_power_w=device.mayr_power_w,
    cassie_time_constant_s=device.cassie_time_constant_s,
    mayr_time_constant_s=device.mayr_time_constant_s,
    transition_power_w=device.transition_power_w,
    conductance_bounds=(device.minimum_conductance_s, device.maximum_conductance_s),
    extinction_current_a=device.extinction_current_a,
    provenance=device.provenance,
)

_surge_state(device::CombinedArcBranch) = (
    accepted_conductance_s=device.accepted_conductance_s,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    dissipated_energy_j=device.dissipated_energy_j,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
    mode=device.mode,
    ignition_count=device.ignition_count,
    extinction_count=device.extinction_count,
    last_transition_time_s=device.last_transition_time_s,
)

_surge_identity_payload(device::FaultArcBranch) = (
    nodes=(device.positive_node, device.negative_node),
    cooling_power_w=device.cooling_power_w,
    time_constant_s=device.time_constant_s,
    arc_length_m=device.arc_length_m,
    conductance_bounds=(device.minimum_conductance_s, device.maximum_conductance_s),
    provenance=device.provenance,
)

_surge_state(device::FaultArcBranch) = (
    accepted_conductance_s=device.accepted_conductance_s,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    dissipated_energy_j=device.dissipated_energy_j,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
)

_surge_identity_payload(device::MetalOxideArrester) = (
    nodes=(device.positive_node, device.negative_node),
    characteristic=device.characteristic,
    gap_sparkover_voltage_v=device.gap_sparkover_voltage_v,
    leakage_conductance_s=device.leakage_conductance_s,
    thermal=(
        device.thermal_capacitance_j_per_k,
        device.thermal_resistance_k_per_w,
        device.ambient_temperature_k,
        device.maximum_temperature_k,
        device.maximum_energy_j,
    ),
    provenance=device.provenance,
)

_surge_state(device::MetalOxideArrester) = (
    conducting=device.conducting,
    failed=device.failed,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    charge_c=device.charge_c,
    absorbed_energy_j=device.absorbed_energy_j,
    temperature_k=device.temperature_k,
    peak_current_a=device.peak_current_a,
    peak_residual_voltage_v=device.peak_residual_voltage_v,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
    sparkover_count=device.sparkover_count,
)

_surge_identity_payload(device::IonizingGroundBranch) = (
    node=device.node,
    linear_conductance_s=device.linear_conductance_s,
    electrode_radius_m=device.electrode_radius_m,
    maximum_ionized_radius_m=device.maximum_ionized_radius_m,
    critical_field_v_per_m=device.critical_field_v_per_m,
    expansion_rate_m_per_v_s=device.expansion_rate_m_per_v_s,
    recovery_rate_per_s=device.recovery_rate_per_s,
    provenance=device.provenance,
)

_surge_state(device::IonizingGroundBranch) = (
    accepted_ionized_radius_m=device.accepted_ionized_radius_m,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    dissipated_energy_j=device.dissipated_energy_j,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
)

_surge_identity_payload(device::DynamicCoronaBranch) = (
    nodes=(device.positive_node, device.negative_node),
    base_capacitance_f=device.base_capacitance_f,
    incremental_capacitance_f_per_v=device.incremental_capacitance_f_per_v,
    onset_voltage_v=device.onset_voltage_v,
    extinction_voltage_v=device.extinction_voltage_v,
    loss_conductance_s=device.loss_conductance_s,
    provenance=device.provenance,
)

_surge_state(device::DynamicCoronaBranch) = (
    active=device.active,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    accepted_charge_c=device.accepted_charge_c,
    dissipated_energy_j=device.dissipated_energy_j,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
)

_surge_identity_payload(device::DisruptiveEffectInsulator) = (
    nodes=(device.positive_node, device.negative_node),
    thresholds=(device.positive_threshold_voltage_v, device.negative_threshold_voltage_v),
    exponents=(device.positive_exponent, device.negative_exponent),
    critical_effects=(device.positive_critical_effect, device.negative_critical_effect),
    conductances=(device.leakage_conductance_s, device.flashover_conductance_s),
    provenance=device.provenance,
)

_surge_state(device::DisruptiveEffectInsulator) = (
    positive_effect=device.positive_effect,
    negative_effect=device.negative_effect,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
    flashed=device.flashed,
    flashover_time_s=device.flashover_time_s,
    flashover_polarity=device.flashover_polarity,
    flashover_count=device.flashover_count,
)

_surge_identity_payload(device::LeaderProgressionInsulator) = (
    nodes=(device.positive_node, device.negative_node),
    gap_length_m=device.gap_length_m,
    inception_fields=(
        device.positive_inception_field_v_per_m,
        device.negative_inception_field_v_per_m,
    ),
    velocity_coefficients=(
        device.positive_velocity_coefficient,
        device.negative_velocity_coefficient,
    ),
    velocity_exponent=device.velocity_exponent,
    conductances=(device.leakage_conductance_s, device.flashover_conductance_s),
    provenance=device.provenance,
)

_surge_state(device::LeaderProgressionInsulator) = (
    leader_length_m=device.leader_length_m,
    leader_polarity=device.leader_polarity,
    incepted=device.incepted,
    flashed=device.flashed,
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
    flashover_time_s=device.flashover_time_s,
)

_surge_identity_payload(device::VacuumInterruptionState) = (
    identity=device.identity,
    chopping_current_a=device.chopping_current_a,
    initial_dielectric_strength_v=device.initial_dielectric_strength_v,
    dielectric_recovery_rate_v_per_s=device.dielectric_recovery_rate_v_per_s,
    maximum_dielectric_strength_v=device.maximum_dielectric_strength_v,
    provenance=device.provenance,
)

_surge_state(device::VacuumInterruptionState) = (
    separation_time_s=device.separation_time_s,
    last_transition_time_s=device.last_transition_time_s,
    last_chop_time_s=device.last_chop_time_s,
    last_restrike_time_s=device.last_restrike_time_s,
    state=device.state,
    chop_count=device.chop_count,
    restrike_count=device.restrike_count,
)

_surge_identity_payload(device::VacuumInterruptionBranch) = (
    nodes=(device.positive_node, device.negative_node),
    state=_surge_identity_payload(device.state),
    arc=_surge_identity_payload(device.arc_owner),
    conductances=(device.closed_conductance_s, device.open_conductance_s),
    reignition_conductance_s=device.reignition_conductance_s,
    provenance=device.provenance,
)

_surge_state(device::VacuumInterruptionBranch) = (
    interruption_state=_surge_state(device.state),
    accepted_voltage_v=device.accepted_voltage_v,
    accepted_current_a=device.accepted_current_a,
    accepted_time_s=device.accepted_time_s,
    trial_step_s=device.trial_step_s,
)

function snapshot_surge_device(device)
    model_type = Symbol(nameof(typeof(device)))
    identity_signature = _surge_identity_signature(
        model_type,
        _surge_identity_payload(device),
    )
    state = _surge_state(device)
    checksum = _surge_snapshot_checksum(1, model_type, identity_signature, state)
    return SurgeDeviceSnapshot(1, model_type, identity_signature, state, checksum)
end

function _validate_surge_snapshot(device, snapshot::SurgeDeviceSnapshot)
    snapshot.schema_version == 1 || throw(ArgumentError(
        "surge device snapshot schema is unsupported",
    ))
    model_type = Symbol(nameof(typeof(device)))
    snapshot.model_type === model_type || throw(ArgumentError(
        "surge device snapshot model type does not match destination",
    ))
    identity_signature = _surge_identity_signature(
        model_type,
        _surge_identity_payload(device),
    )
    snapshot.identity_signature == identity_signature || throw(ArgumentError(
        "surge device snapshot identity does not match destination",
    ))
    expected_checksum = _surge_snapshot_checksum(
        snapshot.schema_version,
        snapshot.model_type,
        snapshot.identity_signature,
        snapshot.state,
    )
    snapshot.checksum == expected_checksum || throw(ArgumentError(
        "surge device snapshot checksum is invalid",
    ))
    return snapshot.state
end

function restore_surge_device!(device::CombinedArcBranch, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.accepted_conductance_s = state.accepted_conductance_s
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.dissipated_energy_j = state.dissipated_energy_j
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    device.mode = state.mode
    device.ignition_count = state.ignition_count
    device.extinction_count = state.extinction_count
    device.last_transition_time_s = state.last_transition_time_s
    return device
end

function restore_surge_device!(device::FaultArcBranch, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.accepted_conductance_s = state.accepted_conductance_s
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.dissipated_energy_j = state.dissipated_energy_j
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    return device
end

function restore_surge_device!(device::MetalOxideArrester, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.conducting = state.conducting
    device.failed = state.failed
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.charge_c = state.charge_c
    device.absorbed_energy_j = state.absorbed_energy_j
    device.temperature_k = state.temperature_k
    device.peak_current_a = state.peak_current_a
    device.peak_residual_voltage_v = state.peak_residual_voltage_v
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    device.sparkover_count = state.sparkover_count
    return device
end

function restore_surge_device!(device::IonizingGroundBranch, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.accepted_ionized_radius_m = state.accepted_ionized_radius_m
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.dissipated_energy_j = state.dissipated_energy_j
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    return device
end

function restore_surge_device!(device::DynamicCoronaBranch, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.active = state.active
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.accepted_charge_c = state.accepted_charge_c
    device.dissipated_energy_j = state.dissipated_energy_j
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    return device
end

function restore_surge_device!(device::DisruptiveEffectInsulator, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.positive_effect = state.positive_effect
    device.negative_effect = state.negative_effect
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    device.flashed = state.flashed
    device.flashover_time_s = state.flashover_time_s
    device.flashover_polarity = state.flashover_polarity
    device.flashover_count = state.flashover_count
    return device
end

function restore_surge_device!(device::LeaderProgressionInsulator, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.leader_length_m = state.leader_length_m
    device.leader_polarity = state.leader_polarity
    device.incepted = state.incepted
    device.flashed = state.flashed
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    device.flashover_time_s = state.flashover_time_s
    return device
end

function restore_surge_device!(device::VacuumInterruptionState, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    device.separation_time_s = state.separation_time_s
    device.last_transition_time_s = state.last_transition_time_s
    device.last_chop_time_s = state.last_chop_time_s
    device.last_restrike_time_s = state.last_restrike_time_s
    device.state = state.state
    device.chop_count = state.chop_count
    device.restrike_count = state.restrike_count
    return device
end

function restore_surge_device!(device::VacuumInterruptionBranch, snapshot::SurgeDeviceSnapshot)
    state = _validate_surge_snapshot(device, snapshot)
    interruption = state.interruption_state
    device.state.separation_time_s = interruption.separation_time_s
    device.state.last_transition_time_s = interruption.last_transition_time_s
    device.state.last_chop_time_s = interruption.last_chop_time_s
    device.state.last_restrike_time_s = interruption.last_restrike_time_s
    device.state.state = interruption.state
    device.state.chop_count = interruption.chop_count
    device.state.restrike_count = interruption.restrike_count
    device.accepted_voltage_v = state.accepted_voltage_v
    device.accepted_current_a = state.accepted_current_a
    device.accepted_time_s = state.accepted_time_s
    device.trial_step_s = state.trial_step_s
    return device
end

//! Finite native motion with complete corotated spatial Maxwell histories.
//!
//! Two exact simultaneous sources surround the original loaded body drift.
//! Every memory amplitude retains its original first-body material-frame basis;
//! the drift carries those coordinates unchanged and therefore changes no
//! material storage. Endpoint operators use actual attachment geometry. This
//! defines an objective corotated rate law, not a finite-strain rest spring.
//! Original Kelvin loss, memory heat and signed drive work remain separate.

use super::super::{
    BodyLoads, Change, ExternalLoadAccount, HeldExternalLoads, HeldSpatialLoads,
    InternalWrenchStepError, LoadChange, accept, add, finite, poll, rounded, state_change,
};
use crate::spatial_loads::internal_impulse::{
    network::{InternalWrenchNetworkError, admit_limit, reserve},
    viscous::{
        ViscousConnectorControls,
        network::{
            driven::DrivenViscousNetworkLimits,
            relaxation::{
                PreparedRelaxationNetwork, RelaxationNetworkError, RelaxationNetworkLink,
                RelaxationNetworkMemory, RelaxationNetworkReceipt, branch_count,
                memory_energy_change_retained, planned_costs,
            },
        },
    },
};
pub use crate::spatial_loads::{
    internal_impulse::viscous::network::relaxation::{
        HeldRelaxationKinetics, RelaxationNetworkLimits,
    },
    internal_motion::viscous::driven::HeldInternalEffort,
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, RigidKineticEnergyChange, Vec3};
use fs_math::dd::Dd;

/// Native physical admission and COMPLETE finite-step resource limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxationNetworkMotionControls {
    pub source: ViscousConnectorControls,
    /// Both endpoint preparations' dense work and complete peak storage,
    /// including the retained initial receipt and native drift state bank.
    pub limits: RelaxationNetworkLimits,
}

/// Accepted original trajectory and separately integrated material heat.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct RelaxationNetworkMotionStep {
    pub duration_s: f64,
    /// The endpoint source owns the accepted state bank.
    pub half_sources: [RelaxationNetworkReceipt; 2],
    pub connector_actuator_work_j: Vec<f64>,
    pub connector_dissipation_j: Vec<f64>,
    pub actuator_work_j: f64,
    pub dissipated_energy_j: f64,
    pub angular_origin_world_m: Vec3,
    pub kinetic_energy_change_j: f64,
    pub memory_energy_change_j: f64,
    pub external_loads: ExternalLoadAccount,
    pub gravitational_potential_energy_change_j: f64,
    pub mechanical_energy_change_j: f64,
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    /// Kinetic + gravity + memory change + heat - all applied source work [J].
    pub work_heat_residual_j: f64,
    pub planned_live_bytes: usize,
    pub planned_preparation_terms: u64,
}
impl RelaxationNetworkMotionStep {
    pub fn states_after(&self) -> &[RigidBodyState] {
        &self.half_sources[1].states_after
    }
    pub fn memory_after(&self) -> &RelaxationNetworkMemory {
        &self.half_sources[1].memory_after
    }
    pub fn into_state_and_memory(self) -> (Vec<RigidBodyState>, RelaxationNetworkMemory) {
        let [_, endpoint] = self.half_sources;
        (endpoint.states_after, endpoint.memory_after)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RelaxationNetworkMotionError {
    Motion(InternalWrenchStepError),
    Network(RelaxationNetworkError),
    Cancelled,
}
impl From<InternalWrenchStepError> for RelaxationNetworkMotionError {
    fn from(error: InternalWrenchStepError) -> Self {
        match error {
            InternalWrenchStepError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}
impl From<DynamicsError> for RelaxationNetworkMotionError {
    fn from(error: DynamicsError) -> Self {
        InternalWrenchStepError::from(error).into()
    }
}
impl From<RelaxationNetworkError> for RelaxationNetworkMotionError {
    fn from(error: RelaxationNetworkError) -> Self {
        match error {
            RelaxationNetworkError::Cancelled => Self::Cancelled,
            error => Self::Network(error),
        }
    }
}
impl From<InternalWrenchNetworkError> for RelaxationNetworkMotionError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        RelaxationNetworkError::from(error).into()
    }
}
impl core::fmt::Display for RelaxationNetworkMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial relaxation network motion refused: {self:?}")
    }
}
impl std::error::Error for RelaxationNetworkMotionError {}

/// Advance the complete graph and all original material histories under held drive.
///
/// Fixed-pose source evolution is simultaneous and exact up to the original
/// numerical admission. Moving geometry uses a symmetric second-order
/// composition, and still needs independent time refinement. The explicitly
/// supplied material factors remain fixed; changing them requires a separate
/// thermodynamic history transfer. No relaxed spring or receiver is inferred.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        None,
        BodyLoads::Unloaded,
        duration_s,
        controls,
        cancelled,
    )
}

/// The same source under native held COM forces, body couples and gravity.
/// Actuator recovery, applied load work and gravitational storage stay distinct
/// from each original connector's independently integrated positive heat.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network_with_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        None,
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// The same source under the original prepared spatial load owner, including
/// material-point force arms and world/body-held free couples.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network_with_spatial_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        None,
        BodyLoads::Spatial(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// The same complete native motion with explicitly held viscous kinetics.
/// Original reference laws and memory bases stay fixed; this bank changes
/// Kelvin viscosity and Maxwell times only. Both source halves use the same
/// admitted coefficients. The composing thermal/time owner resolves their
/// variation and supplies the next interval's bank without resetting history.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network_with_kinetics(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    kinetics: HeldRelaxationKinetics<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        Some(kinetics),
        BodyLoads::Unloaded,
        duration_s,
        controls,
        cancelled,
    )
}

/// The same complete native motion with explicitly held viscous kinetics.
/// Original reference laws and memory bases stay fixed; this bank changes
/// Kelvin viscosity and Maxwell times only. Both source halves use the same
/// admitted coefficients. The composing thermal/time owner resolves their
/// variation and supplies the next interval's bank without resetting history.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network_with_loads_and_kinetics(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    kinetics: HeldRelaxationKinetics<'_>,
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        Some(kinetics),
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// The same complete native motion with explicitly held viscous kinetics.
/// Original reference laws and memory bases stay fixed; this bank changes
/// Kelvin viscosity and Maxwell times only. Both source halves use the same
/// admitted coefficients. The composing thermal/time owner resolves their
/// variation and supplies the next interval's bank without resetting history.
#[allow(clippy::too_many_arguments)]
pub fn step_relaxation_network_with_spatial_loads_and_kinetics(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    kinetics: HeldRelaxationKinetics<'_>,
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        memory,
        efforts,
        Some(kinetics),
        BodyLoads::Spatial(loads),
        duration_s,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn step_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    efforts: &[HeldInternalEffort],
    kinetics: Option<HeldRelaxationKinetics<'_>>,
    loads: BodyLoads<'_>,
    duration_s: f64,
    controls: RelaxationNetworkMotionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkMotionError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse.into());
    }
    if states.len() != properties.len() || links.len() != efforts.len() {
        return Err(RelaxationNetworkError::Invalid(
            "matching original body/inertia banks and one effort per connector required",
        )
        .into());
    }
    loads.check_count(states.len())?;
    let (n, e) = (states.len(), links.len());
    admit_limit(
        "bodies",
        n as u128,
        controls.limits.source.maximum_bodies as u128,
    )?;
    admit_limit(
        "connectors",
        e as u128,
        controls.limits.source.maximum_connectors as u128,
    )?;
    for link in links {
        poll(&mut cancelled)?;
        if link.viscous.bodies[0] == link.viscous.bodies[1]
            || link.viscous.bodies.iter().any(|&i| i >= n)
        {
            return Err(RelaxationNetworkError::Invalid(
                "every connector requires two distinct existing original bodies",
            )
            .into());
        }
    }
    let branches = branch_count(links)?;
    admit_limit(
        "branches",
        branches as u128,
        controls.limits.maximum_branches as u128,
    )?;
    let (source_bytes, source_terms) = planned_costs(n, e, branches, half)?;
    let bytes = source_bytes
        + n as u128
            * (2 * core::mem::size_of::<RigidBodyState>() + core::mem::size_of::<Pose>()) as u128
        + e as u128
            * (4 * core::mem::size_of::<[f64; 6]>() + 8 * core::mem::size_of::<f64>()) as u128;
    let bytes = bytes
        + branches as u128 * (4 * core::mem::size_of::<[f64; 6]>()) as u128
        + e as u128 * (2 * core::mem::size_of::<Vec<[f64; 6]>>()) as u128;
    let terms = 2 * source_terms;
    admit_limit(
        "live bytes",
        bytes,
        controls.limits.source.maximum_live_bytes as u128,
    )?;
    admit_limit(
        "preparation terms",
        terms,
        controls.limits.source.maximum_preparation_terms as u128,
    )?;
    if let Some(kinetics) = kinetics {
        kinetics.validate(links)?;
    }
    let source_limits = RelaxationNetworkLimits {
        source: DrivenViscousNetworkLimits {
            maximum_bodies: n,
            maximum_connectors: e,
            maximum_live_bytes: source_bytes as usize,
            maximum_preparation_terms: source_terms as u64,
        },
        maximum_branches: branches,
    };
    let mut poses = reserve(n)?;
    poses.extend(states.iter().map(|s| s.pose()));
    let mut commands = reserve(e)?;
    commands.resize(e, [0.0; 6]);
    let source = |poses: &[Pose],
                  states: &[RigidBodyState],
                  memory: &[Vec<[f64; 6]>],
                  commands: &mut [[f64; 6]],
                  cancelled: &mut dyn FnMut() -> bool|
     -> Result<RelaxationNetworkReceipt, RelaxationNetworkMotionError> {
        for (i, link) in links.iter().enumerate() {
            poll(&mut || cancelled())?;
            commands[i] = efforts[i].in_first_frame(states[link.viscous.bodies[0]])?;
        }
        let mut prepared = match kinetics {
            Some(kinetics) => PreparedRelaxationNetwork::new_with_kinetics(
                poses,
                properties,
                links,
                kinetics,
                half,
                source_limits,
                &mut *cancelled,
            )?,
            None => PreparedRelaxationNetwork::new(
                poses,
                properties,
                links,
                half,
                source_limits,
                &mut *cancelled,
            )?,
        };
        Ok(prepared.advance(states, memory, commands, controls.source, cancelled)?)
    };
    let initial = source(&poses, states, memory, &mut commands, &mut cancelled)?;
    let mut drifted = reserve(n)?;
    let mut external = LoadChange::ZERO;
    let mut drift_energy = Dd::ZERO;
    let origin = controls.source.impulse.angular_origin_world_m;
    for i in 0..n {
        poll(&mut cancelled)?;
        let before = initial.states_after[i];
        let (step, load) =
            loads.step(i, before, properties[i], duration_s, origin, &mut cancelled)?;
        let after = step.state_after;
        let mut change = state_change(before, after, properties[i], origin)?;
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(properties[i], before, after)?;
        change.energy = kinetic.finish_retained()?;
        accept(
            load.residual(change),
            load.work(),
            controls.source.impulse,
            Some(i),
        )?;
        drift_energy = drift_energy + change.energy - load.work();
        external = external.add(load);
        poses[i] = after.pose();
        drifted.push(after);
    }
    let endpoint = source(
        &poses,
        &drifted,
        &initial.memory_after,
        &mut commands,
        &mut cancelled,
    )?;
    let mut total = Change {
        linear: [Dd::ZERO; 3],
        angular: [Dd::ZERO; 3],
        energy: Dd::ZERO,
    };
    let mut kinetic = RigidKineticEnergyChange::new();
    for i in 0..n {
        poll(&mut cancelled)?;
        let after = endpoint.states_after[i];
        let change = state_change(states[i], after, properties[i], origin)?;
        total.linear = add(total.linear, change.linear);
        total.angular = add(total.angular, change.angular);
        kinetic.add_body(properties[i], states[i], after)?;
    }
    total.energy = kinetic.finish_retained()?;
    let heat =
        Dd::from_f64(initial.dissipated_energy_j) + Dd::from_f64(endpoint.dissipated_energy_j);
    let work = Dd::from_f64(initial.actuator_work_j) + Dd::from_f64(endpoint.actuator_work_j);
    let memory_change = memory_energy_change_retained(memory, &endpoint.memory_after)?;
    let residual = external.residual(total);
    accept(
        residual,
        external.work() + work - heat - memory_change,
        controls.source.impulse,
        None,
    )?;
    let mut connector_dissipation_j = reserve(e)?;
    let mut connector_actuator_work_j = reserve(e)?;
    for i in 0..e {
        poll(&mut cancelled)?;
        connector_dissipation_j.push(finite(
            (Dd::from_f64(initial.connector_dissipation_j[i])
                + Dd::from_f64(endpoint.connector_dissipation_j[i]))
            .to_f64(),
        )?);
        connector_actuator_work_j.push(finite(
            (Dd::from_f64(initial.connector_actuator_work_j[i])
                + Dd::from_f64(endpoint.connector_actuator_work_j[i]))
            .to_f64(),
        )?);
    }
    let result = RelaxationNetworkMotionStep {
        duration_s,
        half_sources: [initial, endpoint],
        connector_actuator_work_j,
        connector_dissipation_j,
        actuator_work_j: finite(work.to_f64())?,
        dissipated_energy_j: finite(heat.to_f64())?,
        angular_origin_world_m: origin,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        memory_energy_change_j: finite(memory_change.to_f64())?,
        external_loads: external.receipt()?,
        gravitational_potential_energy_change_j: finite((-external.gravity_work).to_f64())?,
        mechanical_energy_change_j: finite(
            (total.energy - external.gravity_work + memory_change).to_f64(),
        )?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(residual.linear)?,
        angular_momentum_residual_nms: rounded(residual.angular)?,
        work_heat_residual_j: finite(
            (total.energy + memory_change + heat - work - external.work()).to_f64(),
        )?,
        planned_live_bytes: bytes as usize,
        planned_preparation_terms: terms as u64,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

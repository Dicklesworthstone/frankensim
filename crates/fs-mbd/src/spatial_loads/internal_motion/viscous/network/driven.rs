//! Finite motion with simultaneous driven spatial damping on the original graph.
//!
//! Two complete fixed-pose graph sources surround one original native loaded
//! drift per body. Shared bodies, cycles and parallel links are solved together
//! inside each source. Endpoint geometry and world-held actuator components are
//! recomputed from actual poses; the graph never becomes a sequence of pairs.

use super::super::{
    BodyLoads, Change, ExternalLoadAccount, HeldExternalLoads, HeldSpatialLoads,
    InternalWrenchStepError, LoadChange, accept, add, finite, poll, rounded, state_change,
};
use crate::spatial_loads::internal_impulse::{
    network::{InternalWrenchNetworkError, admit_limit, reserve},
    viscous::{
        ViscousConnectorControls,
        network::{
            ViscousNetworkLink,
            driven::{
                DrivenViscousNetworkError, DrivenViscousNetworkReceipt,
                PreparedDrivenViscousNetwork, planned_costs,
            },
        },
    },
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, RigidKineticEnergyChange, Vec3};
use fs_math::dd::Dd;

pub use crate::spatial_loads::{
    internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits,
    internal_motion::viscous::driven::HeldInternalEffort,
};

/// Native physical admission and COMPLETE finite-step resource limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenViscousNetworkMotionControls {
    pub source: ViscousConnectorControls,
    /// Both endpoint preparations' dense work and complete peak storage,
    /// including the retained initial receipt and native drift state bank.
    pub limits: DrivenViscousNetworkLimits,
}

/// Accepted original trajectory and separately integrated material heat.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct DrivenViscousNetworkMotionStep {
    pub duration_s: f64,
    /// The endpoint source owns the accepted state bank.
    pub half_sources: [DrivenViscousNetworkReceipt; 2],
    pub connector_actuator_work_j: Vec<f64>,
    pub connector_dissipation_j: Vec<f64>,
    pub actuator_work_j: f64,
    pub dissipated_energy_j: f64,
    pub angular_origin_world_m: Vec3,
    pub kinetic_energy_change_j: f64,
    pub external_loads: ExternalLoadAccount,
    pub gravitational_potential_energy_change_j: f64,
    pub mechanical_energy_change_j: f64,
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    /// Mechanical change + heat - internal drive - applied external work [J].
    pub work_heat_residual_j: f64,
    pub planned_live_bytes: usize,
    pub planned_preparation_terms: u64,
}
impl DrivenViscousNetworkMotionStep {
    pub fn states_after(&self) -> &[RigidBodyState] {
        &self.half_sources[1].states_after
    }
    pub fn into_states_after(self) -> Vec<RigidBodyState> {
        let [_, endpoint] = self.half_sources;
        endpoint.states_after
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum DrivenViscousNetworkMotionError {
    Motion(InternalWrenchStepError),
    Network(DrivenViscousNetworkError),
    Cancelled,
}
impl From<InternalWrenchStepError> for DrivenViscousNetworkMotionError {
    fn from(error: InternalWrenchStepError) -> Self {
        match error {
            InternalWrenchStepError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}
impl From<DynamicsError> for DrivenViscousNetworkMotionError {
    fn from(error: DynamicsError) -> Self {
        InternalWrenchStepError::from(error).into()
    }
}
impl From<DrivenViscousNetworkError> for DrivenViscousNetworkMotionError {
    fn from(error: DrivenViscousNetworkError) -> Self {
        match error {
            DrivenViscousNetworkError::Cancelled => Self::Cancelled,
            error => Self::Network(error),
        }
    }
}
impl From<InternalWrenchNetworkError> for DrivenViscousNetworkMotionError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl core::fmt::Display for DrivenViscousNetworkMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven spatial viscous network motion refused: {self:?}")
    }
}
impl std::error::Error for DrivenViscousNetworkMotionError {}

/// Advance the complete graph under held internal actuator effort.
///
/// Fixed-pose source evolution is simultaneous and exact up to the original
/// numerical admission. Moving geometry uses a symmetric second-order
/// composition, and still needs independent time refinement. Thermal receivers,
/// compliant springs and contact constraints are not invented by this entry.
#[allow(clippy::too_many_arguments)]
pub fn step_driven_viscous_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    efforts: &[HeldInternalEffort],
    duration_s: f64,
    controls: DrivenViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousNetworkMotionStep, DrivenViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        efforts,
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
pub fn step_driven_viscous_network_with_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    efforts: &[HeldInternalEffort],
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: DrivenViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousNetworkMotionStep, DrivenViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        efforts,
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// The same source under the original prepared spatial load owner, including
/// material-point force arms and world/body-held free couples.
#[allow(clippy::too_many_arguments)]
pub fn step_driven_viscous_network_with_spatial_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    efforts: &[HeldInternalEffort],
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: DrivenViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousNetworkMotionStep, DrivenViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        efforts,
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
    links: &[ViscousNetworkLink],
    efforts: &[HeldInternalEffort],
    loads: BodyLoads<'_>,
    duration_s: f64,
    controls: DrivenViscousNetworkMotionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousNetworkMotionStep, DrivenViscousNetworkMotionError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse.into());
    }
    if states.len() != properties.len() || links.len() != efforts.len() {
        return Err(DrivenViscousNetworkError::Invalid(
            "matching original body/inertia banks and one effort per connector required",
        )
        .into());
    }
    loads.check_count(states.len())?;
    let (n, e) = (states.len(), links.len());
    admit_limit("bodies", n as u128, controls.limits.maximum_bodies as u128)?;
    admit_limit(
        "connectors",
        e as u128,
        controls.limits.maximum_connectors as u128,
    )?;
    for link in links {
        poll(&mut cancelled)?;
        if link.bodies[0] == link.bodies[1] || link.bodies.iter().any(|&i| i >= n) {
            return Err(DrivenViscousNetworkError::Invalid(
                "every connector requires two distinct existing original bodies",
            )
            .into());
        }
    }
    let (source_bytes, source_terms) = planned_costs(n, e, half)?;
    let bytes = source_bytes
        + n as u128
            * (2 * core::mem::size_of::<RigidBodyState>() + core::mem::size_of::<Pose>()) as u128
        + e as u128
            * (4 * core::mem::size_of::<[f64; 6]>() + 8 * core::mem::size_of::<f64>()) as u128;
    let terms = 2 * source_terms;
    admit_limit(
        "live bytes",
        bytes,
        controls.limits.maximum_live_bytes as u128,
    )?;
    admit_limit(
        "preparation terms",
        terms,
        controls.limits.maximum_preparation_terms as u128,
    )?;
    let source_limits = DrivenViscousNetworkLimits {
        maximum_bodies: n,
        maximum_connectors: e,
        maximum_live_bytes: source_bytes as usize,
        maximum_preparation_terms: source_terms as u64,
    };
    let mut poses = reserve(n)?;
    poses.extend(states.iter().map(|s| s.pose()));
    let mut commands = reserve(e)?;
    commands.resize(e, [0.0; 6]);
    let source = |poses: &[Pose],
                  states: &[RigidBodyState],
                  commands: &mut [[f64; 6]],
                  cancelled: &mut dyn FnMut() -> bool|
     -> Result<DrivenViscousNetworkReceipt, DrivenViscousNetworkMotionError> {
        for (i, link) in links.iter().enumerate() {
            poll(&mut || cancelled())?;
            commands[i] = efforts[i].in_first_frame(states[link.bodies[0]])?;
        }
        let mut prepared = PreparedDrivenViscousNetwork::new(
            poses,
            properties,
            links,
            half,
            source_limits,
            &mut *cancelled,
        )?;
        Ok(prepared.advance(states, commands, controls.source, cancelled)?)
    };
    let initial = source(&poses, states, &mut commands, &mut cancelled)?;
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
    let endpoint = source(&poses, &drifted, &mut commands, &mut cancelled)?;
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
    let residual = external.residual(total);
    accept(
        residual,
        external.work() + work - heat,
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
    let result = DrivenViscousNetworkMotionStep {
        duration_s,
        half_sources: [initial, endpoint],
        connector_actuator_work_j,
        connector_dissipation_j,
        actuator_work_j: finite(work.to_f64())?,
        dissipated_energy_j: finite(heat.to_f64())?,
        angular_origin_world_m: origin,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        external_loads: external.receipt()?,
        gravitational_potential_energy_change_j: finite((-external.gravity_work).to_f64())?,
        mechanical_energy_change_j: finite((total.energy - external.gravity_work).to_f64())?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(residual.linear)?,
        angular_momentum_residual_nms: rounded(residual.angular)?,
        work_heat_residual_j: finite((total.energy + heat - work - external.work()).to_f64())?,
        planned_live_bytes: bytes as usize,
        planned_preparation_terms: terms as u64,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

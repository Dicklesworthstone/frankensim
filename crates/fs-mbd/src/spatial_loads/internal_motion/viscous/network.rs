//! Moving shared-body damping networks with original native work and heat.

/// Simultaneous held internal actuation and full spatial damping.
#[path = "network/driven.rs"]
pub mod driven;

/// Complete corotated Maxwell memory under native body motion.
#[path = "network/relaxation.rs"]
pub mod relaxation;

use super::{
    BodyLoads, Change, Dd, DynamicsError, ExternalLoadAccount, HeldExternalLoads, HeldSpatialLoads,
    InternalWrenchStepError, LoadChange, MassProperties, RigidBodyState, RigidKineticEnergyChange,
    Vec3, ViscousConnectorControls, accept, add, finite, poll, rounded, state_change,
};
use crate::Pose;
use crate::spatial_loads::internal_impulse::network::{
    InternalWrenchNetworkError, admit_limit, reserve,
};
use crate::spatial_loads::internal_impulse::viscous::network::{
    PreparedViscousNetwork, ViscousNetworkError, ViscousNetworkLimits, ViscousNetworkLink,
    ViscousNetworkReceipt, nonzero, planned_live_bytes,
};

/// Source accuracy, physical admission and complete finite-step resource limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousNetworkMotionControls {
    /// Original per-connector, per-body and whole-network physical balance gates.
    pub source: ViscousConnectorControls,
    /// Positive palindromic subdivisions in EACH fixed-geometry half source.
    /// Temporal refinement of moving geometry remains a separate obligation.
    pub source_substeps: usize,
    /// These caps apply to the COMPLETE step: both source stages' pair calls
    /// and peak storage, including retained receipts and native drift states.
    pub limits: ViscousNetworkLimits,
}

/// Accepted moving network, with original body/connector order retained.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct ViscousNetworkMotionStep {
    pub duration_s: f64,
    /// Actual initial and endpoint geometry sources. The second receipt owns
    /// the authoritative endpoint state bank, avoiding an extra body-state copy.
    pub half_sources: [ViscousNetworkReceipt; 2],
    /// Sum of independent source integrals for each original connector [J].
    /// Numerical drift defects and external work are never reclassified as heat.
    pub connector_dissipation_j: Vec<f64>,
    pub dissipated_energy_j: f64,
    pub angular_origin_world_m: Vec3,
    pub kinetic_energy_change_j: f64,
    pub external_loads: ExternalLoadAccount,
    pub gravitational_potential_energy_change_j: f64,
    pub mechanical_energy_change_j: f64,
    /// Sum of separately admitted original native drift defects [J].
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub work_heat_residual_j: f64,
    /// Native pair relaxations across both endpoint sources.
    pub pair_relaxations: usize,
    pub planned_live_bytes: usize,
}
impl ViscousNetworkMotionStep {
    /// All original poses and momenta after the complete accepted interval.
    pub fn states_after(&self) -> &[RigidBodyState] {
        &self.half_sources[1].states_after
    }

    /// Reuse the accepted state bank for the next interval without copying it.
    pub fn into_states_after(self) -> Vec<RigidBodyState> {
        let [_, endpoint] = self.half_sources;
        endpoint.states_after
    }
}

/// Failed motion never publishes a partially advanced body or heat bank.
#[derive(Clone, Debug, PartialEq)]
pub enum ViscousNetworkMotionError {
    Motion(InternalWrenchStepError),
    Network(ViscousNetworkError),
    Cancelled,
}
impl From<InternalWrenchStepError> for ViscousNetworkMotionError {
    fn from(error: InternalWrenchStepError) -> Self {
        match error {
            InternalWrenchStepError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}
impl From<ViscousNetworkError> for ViscousNetworkMotionError {
    fn from(error: ViscousNetworkError) -> Self {
        match error {
            ViscousNetworkError::Cancelled => Self::Cancelled,
            error => Self::Network(error),
        }
    }
}
impl From<InternalWrenchNetworkError> for ViscousNetworkMotionError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        ViscousNetworkError::from(error).into()
    }
}
impl From<DynamicsError> for ViscousNetworkMotionError {
    fn from(error: DynamicsError) -> Self {
        Self::Motion(error.into())
    }
}
impl core::fmt::Display for ViscousNetworkMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "viscous network motion refused: {self:?}")
    }
}
impl std::error::Error for ViscousNetworkMotionError {}

/// Advance arbitrary shared-body, cyclic or parallel passive connector graphs.
///
/// Two original fixed-geometry network half sources surround one native rigid
/// motion per body, including isolated bodies. Endpoint damping is re-prepared
/// from actual evolved poses, full inertias and transported material arms. No
/// global dense matrix, new damping law, constraint or hidden time step is used.
///
/// The composition is second order for smooth motion under time refinement;
/// the original source splitting is separately second order under its explicit
/// subdivisions. Source stability and closed heat accounts do not establish
/// accuracy at unresolved moving geometry or coarse shared-body coupling.
///
/// Storage is linear in bodies plus connectors. Both source stages' pair-call
/// counts and the complete peak payload are admitted before allocation. Each
/// source, EACH native body drift and the final original-to-endpoint account
/// must pass physical gates. Refusal/cancellation leaves every input unchanged.
pub fn step_viscous_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    duration_s: f64,
    controls: ViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousNetworkMotionStep, ViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        BodyLoads::Unloaded,
        duration_s,
        controls,
        cancelled,
    )
}

/// The same complete damping network under held native COM forces, body couples
/// and uniform gravity. Applied source work and gravity remain distinct from heat.
pub fn step_viscous_network_with_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: ViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousNetworkMotionStep, ViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// The same network under prepared material-point forces and world/body couples.
/// The original spatial midpoint owner retains rotating-arm work and its native
/// rotation/nonlinear-solve gates, including loads on isolated bodies.
pub fn step_viscous_network_with_spatial_loads(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: ViscousNetworkMotionControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousNetworkMotionStep, ViscousNetworkMotionError> {
    step_network(
        states,
        properties,
        links,
        BodyLoads::Spatial(loads),
        duration_s,
        controls,
        cancelled,
    )
}

fn step_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    loads: BodyLoads<'_>,
    duration_s: f64,
    controls: ViscousNetworkMotionControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ViscousNetworkMotionStep, ViscousNetworkMotionError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse.into());
    }
    if states.len() != properties.len() || controls.source_substeps == 0 {
        return Err(ViscousNetworkError::Invalid(
            "matching bodies/inertias and positive source substeps required",
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
    let mut active = 0_usize;
    for link in links {
        poll(&mut cancelled)?;
        if nonzero(&link.damping) {
            active += 1;
        }
    }
    let calls = if active == 0 {
        0
    } else {
        (2 * active as u128 - 1)
            .checked_mul(controls.source_substeps as u128)
            .and_then(|calls| calls.checked_mul(2))
            .ok_or(ViscousNetworkError::Invalid(
                "complete pair-call count overflows",
            ))?
    };
    admit_limit(
        "pair relaxations",
        calls,
        controls.limits.maximum_pair_relaxations as u128,
    )?;
    let source_bytes = planned_live_bytes(n, e, active);
    // One source preparation/relaxation coexists with the first receipt,
    // evolved drift bank and reusable pose bank. Include the combined edge
    // heat output conservatively even before that output is allocated.
    let bytes = source_bytes
        + n as u128
            * (2 * core::mem::size_of::<RigidBodyState>() + core::mem::size_of::<Pose>()) as u128
        + e as u128 * (2 * core::mem::size_of::<f64>()) as u128;
    admit_limit(
        "live bytes",
        bytes,
        controls.limits.maximum_live_bytes as u128,
    )?;
    let mut source_limits = controls.limits;
    source_limits.maximum_pair_relaxations = (calls / 2) as usize;
    source_limits.maximum_live_bytes = source_bytes as usize;
    let mut poses = reserve(n)?;
    for state in states {
        poll(&mut cancelled)?;
        poses.push(state.pose());
    }
    let source =
        |poses: &[Pose], states: &[RigidBodyState], cancelled: &mut dyn FnMut() -> bool| {
            // The prepared source is released before any next preparation, so the
            // complete peak never contains two cached connector banks.
            let prepared = PreparedViscousNetwork::new(
                poses,
                properties,
                links,
                half,
                controls.source_substeps,
                source_limits,
                &mut *cancelled,
            )?;
            prepared.relax(states, controls.source, cancelled)
        };
    let initial = source(&poses, states, &mut cancelled)?;
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
    let endpoint = source(&poses, &drifted, &mut cancelled)?;
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
    let residual = external.residual(total);
    accept(
        residual,
        external.work() - heat,
        controls.source.impulse,
        None,
    )?;
    let mut connector_dissipation_j = reserve(e)?;
    for (&first, &second) in initial
        .connector_dissipation_j
        .iter()
        .zip(&endpoint.connector_dissipation_j)
    {
        poll(&mut cancelled)?;
        connector_dissipation_j.push(finite(
            (Dd::from_f64(first) + Dd::from_f64(second)).to_f64(),
        )?);
    }
    let result = ViscousNetworkMotionStep {
        duration_s,
        half_sources: [initial, endpoint],
        connector_dissipation_j,
        dissipated_energy_j: finite(heat.to_f64())?,
        angular_origin_world_m: origin,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        external_loads: external.receipt()?,
        gravitational_potential_energy_change_j: finite((-external.gravity_work).to_f64())?,
        mechanical_energy_change_j: finite((total.energy - external.gravity_work).to_f64())?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(residual.linear)?,
        angular_momentum_residual_nms: rounded(residual.angular)?,
        work_heat_residual_j: finite((total.energy + heat - external.work()).to_f64())?,
        pair_relaxations: calls as usize,
        planned_live_bytes: bytes as usize,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

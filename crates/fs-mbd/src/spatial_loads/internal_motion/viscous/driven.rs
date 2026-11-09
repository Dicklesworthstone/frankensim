//! Moving finite bodies with simultaneous held actuation and spatial damping.
//!
//! Each frozen-pose source solves the held force/couple and full passive law
//! together. Two such half sources surround the original native loaded motion;
//! the endpoint source uses the actual new attachments, inertia orientation and
//! actuator frame. Actuator work and physical heat are independent integrals.
//! This second-order composition still requires temporal refinement for moving
//! geometry; an exact stiff local source is not an exact finite trajectory.

use super::{
    BodyLoads, Change, ExternalLoadAccount, HeldExternalLoads, HeldSpatialLoads,
    InternalWrenchStepError, LoadChange, ViscousMotionLaw, accept, add, finite, poll, rounded,
    state_change,
};
use crate::spatial_loads::internal_impulse::viscous::{
    ViscousConnectorControls,
    driven::{
        DrivenViscousConnectorError, DrivenViscousConnectorReceipt, PreparedDrivenViscousConnector,
    },
};
use crate::{DynamicsError, MassProperties, RigidBodyState, RigidKineticEnergyChange, Vec3};
use fs_math::dd::Dd;

/// A force and free couple on attachment B, held in the declared frame.
///
/// Both variants contain [force N; 3], followed by [free couple N m; 3].
/// Attachment A receives the transported opposite wrench. These are internal
/// actuator commands: their signed work may supply or recover mechanical energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HeldInternalEffort {
    /// Components follow body A's original material axes as that body rotates.
    MaterialFirst([f64; 6]),
    /// Components remain fixed in world axes throughout the step.
    World([f64; 6]),
}

impl HeldInternalEffort {
    pub(in crate::spatial_loads) fn in_first_frame(
        self,
        first: RigidBodyState,
    ) -> Result<[f64; 6], DynamicsError> {
        let authored = match self {
            Self::MaterialFirst(value) | Self::World(value) => value,
        };
        if authored.iter().any(|x| !x.is_finite()) {
            return Err(DynamicsError::NonFinite("driven_viscous_motion.effort"));
        }
        let value = match self {
            Self::MaterialFirst(value) => value,
            Self::World(value) => {
                let rotation = first.pose().orientation();
                let force = rotation.rotate_world_to_body(Vec3::new(value[0], value[1], value[2]));
                let couple = rotation.rotate_world_to_body(Vec3::new(value[3], value[4], value[5]));
                [force.x, force.y, force.z, couple.x, couple.y, couple.z]
            }
        };
        if value.iter().any(|x| !x.is_finite()) {
            return Err(DynamicsError::NonFinite(
                "driven_viscous_motion.rotated_effort",
            ));
        }
        Ok(value)
    }
}

/// Accepted finite motion with distinct actuation, passive heat and environment.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct DrivenViscousMotionStep {
    pub duration_s: f64,
    pub states_after: [RigidBodyState; 2],
    /// Original native geometry, simultaneous impulses and integrals at each end.
    pub half_sources: [DrivenViscousConnectorReceipt; 2],
    /// Sum of independently integrated held-actuator work in the two sources [J].
    /// Negative work is actuator recovery, without changing the sign of heat.
    pub actuator_work_j: f64,
    pub dissipated_energy_j: f64,
    pub angular_origin_world_m: Vec3,
    pub kinetic_energy_change_j: f64,
    pub external_loads: ExternalLoadAccount,
    pub gravitational_potential_energy_change_j: f64,
    pub mechanical_energy_change_j: f64,
    /// Original native motion's kinetic change minus its source work [J].
    /// This numerical defect is separately admitted and is never physical heat.
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    /// Mechanical change + heat - actuator work - external applied work [J].
    pub work_heat_residual_j: f64,
}

/// A refusal returns no partial endpoint, accepted work or heat.
#[derive(Clone, Debug, PartialEq)]
pub enum DrivenViscousMotionError {
    Motion(InternalWrenchStepError),
    Source(DrivenViscousConnectorError),
}
impl From<InternalWrenchStepError> for DrivenViscousMotionError {
    fn from(error: InternalWrenchStepError) -> Self {
        Self::Motion(error)
    }
}
impl From<DynamicsError> for DrivenViscousMotionError {
    fn from(error: DynamicsError) -> Self {
        Self::Motion(InternalWrenchStepError::Dynamics(error))
    }
}
impl From<DrivenViscousConnectorError> for DrivenViscousMotionError {
    fn from(error: DrivenViscousConnectorError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for DrivenViscousMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven viscous connector motion refused: {self:?}")
    }
}
impl std::error::Error for DrivenViscousMotionError {}

/// Advance a freely moving pair with simultaneous internal actuation and damping.
///
/// Full finite-body mobility and all six original passive force/couple channels
/// participate in both sources. A world-held command is re-expressed in the
/// actual endpoint material axes; a material-held command follows body A.
/// No force kick is separated from its damping response.
///
/// Two bounded six-coordinate preparations and two original native body steps
/// perform the work. Each preparation is released before the next. Cancellation
/// is observed within source preparation, between native steps and before
/// publication. All accepted state remains local until every source, native
/// drift and whole-trajectory momentum/work gate passes.
///
/// Refine time for finite-rotation, changing-geometry and stiff-trajectory
/// accuracy. No contact constraint, spring storage or thermal destination is
/// introduced. The initial state, material law and effort are immutable.
pub fn step_driven_viscous_connector(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    effort: HeldInternalEffort,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousMotionStep, DrivenViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        effort,
        BodyLoads::Unloaded,
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance the driven damped pair under native held COM forces, body couples
/// and uniform gravity. Environment work and gravity potential remain separate
/// from internal actuator work and viscous heat.
#[allow(clippy::too_many_arguments)]
pub fn step_driven_viscous_connector_with_loads(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    effort: HeldInternalEffort,
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousMotionStep, DrivenViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        effort,
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance under original prepared material-point forces and body/world free
/// couples. The native spatial midpoint owner retains its force-arm work,
/// source angular impulse and finite-rotation/nonlinear-solve limits.
#[allow(clippy::too_many_arguments)]
pub fn step_driven_viscous_connector_with_spatial_loads(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    effort: HeldInternalEffort,
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousMotionStep, DrivenViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        effort,
        BodyLoads::Spatial(loads),
        duration_s,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn step_pair(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    effort: HeldInternalEffort,
    loads: BodyLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DrivenViscousMotionStep, DrivenViscousMotionError> {
    poll(&mut cancelled)?;
    loads.check_count(2)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse.into());
    }
    let origin = controls.impulse.angular_origin_world_m;
    let source = |states: [RigidBodyState; 2],
                  cancelled: &mut dyn FnMut() -> bool|
     -> Result<DrivenViscousConnectorReceipt, DrivenViscousMotionError> {
        let components = effort.in_first_frame(states[0])?;
        let mut prepared = PreparedDrivenViscousConnector::new(
            states.map(|body| body.pose()),
            properties,
            law.arms_body_m,
            law.damping.clone(),
            half,
            &mut *cancelled,
        )?;
        Ok(prepared.advance(states, components, controls, cancelled)?)
    };
    let initial = source(states, &mut cancelled)?;
    let kicked = initial.impulse.bodies.map(|body| body.state_after);
    let mut drifted = kicked;
    let mut external = LoadChange::ZERO;
    let mut drift_energy = Dd::ZERO;
    for i in 0..2 {
        poll(&mut cancelled)?;
        let (step, load) = loads.step(
            i,
            kicked[i],
            properties[i],
            duration_s,
            origin,
            &mut cancelled,
        )?;
        drifted[i] = step.state_after;
        let mut change = state_change(kicked[i], drifted[i], properties[i], origin)?;
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(properties[i], kicked[i], drifted[i])?;
        change.energy = kinetic.finish_retained()?;
        accept(
            load.residual(change),
            load.work(),
            controls.impulse,
            Some(i),
        )?;
        drift_energy = drift_energy + change.energy - load.work();
        external = external.add(load);
    }
    poll(&mut cancelled)?;
    let endpoint = source(drifted, &mut cancelled)?;
    let states_after = endpoint.impulse.bodies.map(|body| body.state_after);
    let mut total = Change {
        linear: [Dd::ZERO; 3],
        angular: [Dd::ZERO; 3],
        energy: Dd::ZERO,
    };
    let mut kinetic = RigidKineticEnergyChange::new();
    for i in 0..2 {
        let change = state_change(states[i], states_after[i], properties[i], origin)?;
        total.linear = add(total.linear, change.linear);
        total.angular = add(total.angular, change.angular);
        kinetic.add_body(properties[i], states[i], states_after[i])?;
    }
    total.energy = kinetic.finish_retained()?;
    let heat =
        Dd::from_f64(initial.dissipated_energy_j) + Dd::from_f64(endpoint.dissipated_energy_j);
    let actuator = Dd::from_f64(initial.actuator_work_j) + Dd::from_f64(endpoint.actuator_work_j);
    let residual = external.residual(total);
    accept(
        residual,
        external.work() + actuator - heat,
        controls.impulse,
        None,
    )?;
    let result = DrivenViscousMotionStep {
        duration_s,
        states_after,
        half_sources: [initial, endpoint],
        actuator_work_j: finite(actuator.to_f64())?,
        dissipated_energy_j: finite(heat.to_f64())?,
        angular_origin_world_m: origin,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        external_loads: external.receipt()?,
        gravitational_potential_energy_change_j: finite((-external.gravity_work).to_f64())?,
        mechanical_energy_change_j: finite((total.energy - external.gravity_work).to_f64())?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(residual.linear)?,
        angular_momentum_residual_nms: rounded(residual.angular)?,
        work_heat_residual_j: finite((total.energy + heat - actuator - external.work()).to_f64())?,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

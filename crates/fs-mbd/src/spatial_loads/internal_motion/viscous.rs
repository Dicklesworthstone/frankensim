//! Finite motion under passive full spatial damping and original external loads.
//!
//! Exact fixed-geometry half relaxations surround the original loaded native
//! rigid-body step. The second source is prepared at the actual endpoint poses:
//! material axes, attachment transport and both full inertia mobilities evolve.
//! This symmetric composition is second order for the smooth frozen-coefficient
//! constitutive law under temporal refinement. A stiff source remains passive,
//! but a closed momentum/heat budget does not establish trajectory accuracy.

use super::{
    BodyLoads, Change, ExternalLoadAccount, HeldExternalLoads, HeldSpatialLoads,
    InternalWrenchStepError, LoadChange, accept, add, finite, poll, rounded, state_change,
};
use crate::spatial_loads::internal_impulse::viscous::{
    PassiveBushingMatrix, PreparedViscousConnector, ViscousConnectorControls,
    ViscousConnectorError, ViscousConnectorReceipt,
};
use crate::{DynamicsError, MassProperties, RigidBodyState, RigidKineticEnergyChange, Vec3};
use fs_math::dd::Dd;

/// Finite motion with shared-body networks of the same passive spatial laws.
#[path = "viscous/network.rs"]
pub mod network;

/// Finite motion with simultaneous held actuation and full spatial damping.
#[path = "viscous/driven.rs"]
pub mod driven;

/// One passive spatial law attached to the original material frames.
///
/// The full SI damping matrix acts on the transported relative twist at B,
/// expressed in A's material axes. The same authored matrix is used at both
/// endpoint poses; no stiffness, rest orientation or thermal law is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscousMotionLaw {
    /// Original COM-relative material attachment coordinates [m].
    pub arms_body_m: [Vec3; 2],
    pub damping: PassiveBushingMatrix,
}

/// Accepted finite motion, independent heat and original external-work accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct ViscousMotionStep {
    pub duration_s: f64,
    pub states_after: [RigidBodyState; 2],
    /// Actual initial and endpoint source geometry, native impulses and heat.
    pub half_sources: [ViscousConnectorReceipt; 2],
    /// Sum of the two independently integrated source losses [J].
    /// Numerical native-motion defects are never assigned to this heat.
    pub dissipated_energy_j: f64,
    pub angular_origin_world_m: Vec3,
    /// Original native fixed-inertia kinetic storage difference [J].
    pub kinetic_energy_change_j: f64,
    pub external_loads: ExternalLoadAccount,
    pub gravitational_potential_energy_change_j: f64,
    pub mechanical_energy_change_j: f64,
    /// Native loaded motion's kinetic change minus its physical source work [J].
    /// Each body's contribution is admitted separately before summation.
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    /// Actual mechanical energy change + viscous heat - external applied work.
    pub work_heat_residual_j: f64,
}

/// Neither input state, source law nor prepared external load is mutated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ViscousMotionError {
    Motion(InternalWrenchStepError),
    Source(ViscousConnectorError),
}
impl From<InternalWrenchStepError> for ViscousMotionError {
    fn from(error: InternalWrenchStepError) -> Self {
        Self::Motion(error)
    }
}
impl From<DynamicsError> for ViscousMotionError {
    fn from(error: DynamicsError) -> Self {
        Self::Motion(InternalWrenchStepError::Dynamics(error))
    }
}
impl From<ViscousConnectorError> for ViscousMotionError {
    fn from(error: ViscousConnectorError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for ViscousMotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "viscous connector motion refused: {self:?}")
    }
}
impl std::error::Error for ViscousMotionError {}

/// Advance a freely moving pair coupled by a full passive spatial damper.
///
/// Uses the original finite-body inertia, quaternion and implicit free-rotation
/// owners. Both full six-coordinate source maps use the native internal-wrench
/// impulse owner and independently integrated dissipation. Initial geometry is
/// never reused for the endpoint source after either body has moved.
///
/// The original impulse controls gate each source, EACH native free-body drift,
/// and the complete original-to-final linear/angular momentum and heat balance.
/// Actual quadratic kinetic increments use the original compensated storage
/// owner, including complete inertia tensors and cancellation between bodies.
///
/// Two bounded six-coordinate source preparations surround two native body
/// steps. Spectral scratch stays with the original small-eigenproblem owner;
/// there is no hidden subdivision or growing state history. Cancellation
/// is observed between these bounded operations and before return. Refusal
/// returns no candidate endpoint; both input values remain reusable.
///
/// Refine time for moving-geometry accuracy, including stiff modes. The source
/// duration is not a free-flight clock for either isolated half relaxation.
/// No collision, joint constraint or thermal destination is imposed.
pub fn step_viscous_connector(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousMotionStep, ViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        BodyLoads::Unloaded,
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance the damped pair under original held COM forces, body couples and
/// uniform gravity. The native loaded-motion owner supplies independent work
/// and spatial impulses. Gravity potential and viscous heat remain distinct.
pub fn step_viscous_connector_with_loads(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    loads: HeldExternalLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousMotionStep, ViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        BodyLoads::Held(loads),
        duration_s,
        controls,
        cancelled,
    )
}

/// Advance the damped pair under the original prepared material-point forces
/// and body/world free couples. Rotating lever-arm work and independent source
/// angular impulses remain with the native spatial midpoint solve; its original
/// finite-rotation and nonlinear-solve limits still apply.
pub fn step_viscous_connector_with_spatial_loads(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    loads: HeldSpatialLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ViscousMotionStep, ViscousMotionError> {
    step_pair(
        states,
        properties,
        law,
        BodyLoads::Spatial(loads),
        duration_s,
        controls,
        cancelled,
    )
}

fn step_pair(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    law: &ViscousMotionLaw,
    loads: BodyLoads<'_>,
    duration_s: f64,
    controls: ViscousConnectorControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ViscousMotionStep, ViscousMotionError> {
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
    let source = |states: [RigidBodyState; 2], cancelled: &mut dyn FnMut() -> bool| {
        PreparedViscousConnector::new(
            states.map(|body| body.pose()),
            properties,
            law.arms_body_m,
            law.damping.clone(),
            half,
            cancelled,
        )
    };
    let prepared = source(states, &mut cancelled)?;
    let initial = prepared.relax(states, controls, &mut cancelled)?;
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
    let prepared = source(drifted, &mut cancelled)?;
    let endpoint = prepared.relax(drifted, controls, &mut cancelled)?;
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
    let residual = external.residual(total);
    accept(residual, external.work() - heat, controls.impulse, None)?;
    let result = ViscousMotionStep {
        duration_s,
        states_after,
        half_sources: [initial, endpoint],
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
    };
    poll(&mut cancelled)?;
    Ok(result)
}

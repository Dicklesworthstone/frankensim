//! One native free drift per body between complete internal network kicks.

use super::{
    Change, Dd, DynamicsError, Gravity, InternalWrenchImpulse, InternalWrenchStepError,
    MassProperties, PrescribedInternalWrench, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
    accept, add, finite, poll, rounded, scaled_impulse, state_change,
};
use crate::StepReceipt;
use crate::spatial_loads::internal_impulse::network::{
    InternalWrenchLink, InternalWrenchNetworkControls, InternalWrenchNetworkError,
    InternalWrenchNetworkReceipt, admit_limit, apply_internal_wrench_network, planned_live_bytes,
    reserve, retained_bytes,
};

/// A world-held internal drive joining two original native bodies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedInternalLink {
    pub bodies: [usize; 2],
    pub drive: PrescribedInternalWrench,
}

/// Complete finite motion and its internal momentum and actuator-work account.
#[derive(Clone, Debug, PartialEq)]
pub struct InternalWrenchNetworkStep {
    pub duration_s: f64,
    pub states_after: Vec<RigidBodyState>,
    /// Each stage coalesces all connectors before any native body update.
    pub half_impulses: [InternalWrenchNetworkReceipt; 2],
    /// Exactly one original unconstrained free step per body, in body order.
    pub free_steps: Vec<StepReceipt>,
    pub angular_origin_world_m: Vec3,
    pub kinetic_energy_change_j: f64,
    pub actuator_work_j: f64,
    /// Actual native free-drift energy defect, not a heat or dissipation law.
    pub native_drift_energy_defect_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub work_energy_residual_j: f64,
    /// Complete peak including both stage receipts and the active native bank.
    pub planned_live_bytes: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum InternalWrenchNetworkStepError {
    Motion(InternalWrenchStepError),
    Network(InternalWrenchNetworkError),
    Cancelled,
}
impl From<InternalWrenchStepError> for InternalWrenchNetworkStepError {
    fn from(error: InternalWrenchStepError) -> Self {
        match error {
            InternalWrenchStepError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}
impl From<InternalWrenchNetworkError> for InternalWrenchNetworkStepError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        match error {
            InternalWrenchNetworkError::Cancelled => Self::Cancelled,
            error => Self::Network(error),
        }
    }
}
impl From<DynamicsError> for InternalWrenchNetworkStepError {
    fn from(error: DynamicsError) -> Self {
        Self::Motion(error.into())
    }
}
impl core::fmt::Display for InternalWrenchNetworkStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "internal wrench network motion refused: {self:?}")
    }
}
impl std::error::Error for InternalWrenchNetworkStepError {}

/// Advance a whole network of finite internally driven bodies on one clock.
///
/// Each endpoint stage uses the original simultaneous internal-impulse network.
/// The transported opposite couples follow the actual moving attachments.
/// Shared-body resultants are applied once, and each connector's signed work
/// uses the complete midpoint velocity, including every other active connector.
/// The intervening free motion is exactly the original pair trajectory's
/// `RigidBodyIntegrator::step` with zero external wrench and zero gravity,
/// called once per body, including isolated bodies.
///
/// The original pair's compensated finite-state account checks EACH native
/// free drift before the whole network's momentum/energy/work admission. No
/// individual drift defect may cancel against another body to evade its gate.
/// The symmetric composition is second order for smooth held world wrenches;
/// continuous source work and orientation still require step refinement.
///
/// All storage is reserved after admitting the complete peak payload. Work is
/// linear in body and connector counts, bounded by the caller's two count caps.
/// There is no per-connector drift, pairwise trajectory splitting, dense graph
/// solve or inferred constraint. Any refusal or cancellation discards all
/// candidate bodies, impulses and work. Springs, joints, contact, controller
/// feedback, external loads and thermal allocation retain their own laws.
pub fn step_internal_wrench_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[PrescribedInternalLink],
    duration_s: f64,
    controls: InternalWrenchNetworkControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InternalWrenchNetworkStep, InternalWrenchNetworkStepError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse.into());
    }
    if states.len() != properties.len() {
        return Err(InternalWrenchNetworkError::Invalid("body/property counts differ").into());
    }
    let (n, e) = (states.len(), links.len());
    admit_limit("bodies", n as u128, controls.maximum_bodies as u128)?;
    admit_limit("connectors", e as u128, controls.maximum_connectors as u128)?;
    // At the second kick, one retained receipt coexists with a complete active
    // impulse transaction, the only state vector, native steps and commands.
    let native_bytes = planned_live_bytes(n, e);
    let bytes = native_bytes
        + retained_bytes(n, e)
        + n as u128
            * (core::mem::size_of::<RigidBodyState>() + core::mem::size_of::<StepReceipt>())
                as u128
        + e as u128 * core::mem::size_of::<InternalWrenchLink>() as u128;
    admit_limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut half_commands = reserve(e)?;
    let mut candidate = reserve(n)?;
    let mut free_steps = reserve(n)?;
    for link in links {
        poll(&mut cancelled)?;
        half_commands.push(InternalWrenchLink {
            bodies: link.bodies,
            impulse: InternalWrenchImpulse {
                arms_body_m: link.drive.arms_body_m,
                impulse_on_second_world_ns: scaled_impulse(
                    link.drive.force_on_second_world_n,
                    half,
                )?,
                couple_on_second_world_nms: scaled_impulse(
                    link.drive.couple_on_second_world_nm,
                    half,
                )?,
            },
        });
    }
    let mut native_controls = controls;
    native_controls.maximum_live_bytes = native_bytes as usize;
    let first = apply_internal_wrench_network(
        states,
        properties,
        &half_commands,
        native_controls,
        &mut cancelled,
    )?;
    let free = RigidBodyIntegrator::new(Gravity::ZERO);
    let mut drift_energy = Dd::ZERO;
    for (i, event) in first.bodies.iter().enumerate() {
        poll(&mut cancelled)?;
        let before = event.state_after;
        let step = free.step(before, properties[i], Wrench::ZERO, duration_s)?;
        poll(&mut cancelled)?;
        let change = state_change(
            before,
            step.state_after,
            properties[i],
            controls.balance.angular_origin_world_m,
        )?;
        accept(change, Dd::ZERO, controls.balance, Some(i))?;
        drift_energy = drift_energy + change.energy;
        candidate.push(step.state_after);
        free_steps.push(step);
    }
    let second = apply_internal_wrench_network(
        &candidate,
        properties,
        &half_commands,
        native_controls,
        &mut cancelled,
    )?;
    let mut total = Change {
        linear: [Dd::ZERO; 3],
        angular: [Dd::ZERO; 3],
        energy: Dd::ZERO,
    };
    for i in 0..n {
        poll(&mut cancelled)?;
        candidate[i] = second.bodies[i].state_after;
        let change = state_change(
            states[i],
            candidate[i],
            properties[i],
            controls.balance.angular_origin_world_m,
        )?;
        total.linear = add(total.linear, change.linear);
        total.angular = add(total.angular, change.angular);
        total.energy = total.energy + change.energy;
    }
    let work = Dd::from_f64(first.connector_work_j) + Dd::from_f64(second.connector_work_j);
    accept(total, work, controls.balance, None)?;
    let result = InternalWrenchNetworkStep {
        duration_s,
        states_after: candidate,
        half_impulses: [first, second],
        free_steps,
        angular_origin_world_m: controls.balance.angular_origin_world_m,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        actuator_work_j: finite(work.to_f64())?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(total.linear)?,
        angular_momentum_residual_nms: rounded(total.angular)?,
        work_energy_residual_j: finite((total.energy - work).to_f64())?,
        planned_live_bytes: bytes as usize,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

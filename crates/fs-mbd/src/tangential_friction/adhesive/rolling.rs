//! Source-owned adhesive rolling on two finite rigid bodies at held geometry.
//!
//! The longitudinal tangent is one fixed bending axis. Its actual inverse-
//! inertia response couples the source's exact work-mean resisting moment to
//! a discrete-gradient midpoint impulse phase. Equal/opposite free couples
//! change the original momenta; compensated native kinetic change independently
//! checks recoverable rolling storage and irreversible rolling loss.
//!
//! This supplies no force reaction enforcing no-slip rolling. The normal source,
//! attached area, frame, material arms, poses and mass properties remain fixed.
//! Full inertia can change angular rates on other axes; those are not clamped.
//! Constitutive path integration is exact, while the inertial midpoint phase is
//! a time discretization. A reversing step is admitted only when a guarded
//! native mode-energy bound excludes yield anywhere on its elastic orbit.
//! Otherwise it refuses: an elastic net step can hide yielding at a turn.
//! Held applied COM forces and body torques can drive the same solve. Their
//! native work is kept separate from contact storage and loss; a driven
//! reversal requires a refined interval rather than the unforced orbit bound.
//! Pose evolution, normal-area evolution, history transport, arbitrary vector
//! rolling and heat allocation remain with their physical owners.

use super::super::hertz::RigidHertzError;
use super::super::{RigidContactBinding, RigidFrictionError, poll, vector};
use super::moment::{self, HeldGeometry, MomentError, MomentReport};
use crate::{DynamicsError, MassProperties, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::rolling::{
    SourceJkrRollingContact, SourceJkrRollingError, SourceJkrRollingState, SourceJkrRollingStep,
};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;

/// Physical-step and native work/momentum gates shared with bonded torsion.
pub use super::moment::RigidAdhesiveMomentControls as RigidAdhesiveRollingControls;

/// Physical pair and coincident reaction markers; longitudinal is the rolling axis.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveRollingContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceJkrRollingContact,
}

/// Sealed source memory bound to the original pair, frame and held geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAdhesiveRollingState {
    source: SourceJkrRollingState,
    geometry: HeldGeometry,
}
impl RigidAdhesiveRollingState {
    pub const fn source_state(&self) -> SourceJkrRollingState {
        self.source
    }
    pub fn accumulated_angle_rad(&self) -> f64 {
        self.source.accumulated_angle_rad()
    }
    pub fn elastic_shift_m(&self) -> f64 {
        self.source.elastic_shift_m()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingReport {
    pub contact: SourceJkrRollingStep,
    pub angle_increment_rad: f64,
    pub rolling_rate_before_rad_s: f64,
    pub rolling_rate_after_rad_s: f64,
    /// Axis dot (I_A^-1 + I_B^-1) axis, with each actual world inertia [1/(kg m2)].
    pub rotational_mobility_per_kg_m2: f64,
    /// Contact free-couple impulse only, excluding the external body torques.
    pub angular_impulse_on_a_world_nms: Vec3,
    /// Held world COM forces and body-frame torques, in the binding's A/B order.
    pub external_wrenches: [Wrench; 2],
    /// Actual midpoint COM travel, in A/B order; poses remain held.
    pub body_displacement_world_m: [Vec3; 2],
    /// Integrated held-frame angular velocities, not finite orientation updates.
    pub angular_displacement_body_rad: [Vec3; 2],
    /// Native actual-momentum change, without subtracting large kinetic totals.
    pub kinetic_energy_change_j: f64,
    /// Contact free-couple work on the actual native midpoint rates.
    pub midpoint_impulse_work_j: f64,
    /// Applied force/torque work on the actual COM/angular travel.
    pub external_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub angle_residual_rad: f64,
    pub angular_velocity_residual_rad_s: f64,
    pub body_angular_impulse_residual_nms: [Vec3; 2],
    pub body_linear_impulse_residual_ns: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    /// Actual pair momentum change minus commanded external impulse.
    pub linear_momentum_residual_ns: Vec3,
    /// Orbital-plus-spin change minus external moment about the world origin.
    pub angular_momentum_residual_nms: Vec3,
    pub evaluations: usize,
}

/// Complete staged proposal; the A/B ordering is the original binding ordering.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidAdhesiveRollingState,
    pub report: RigidAdhesiveRollingReport,
}

#[derive(Debug)]
pub enum RigidAdhesiveRollingError {
    Invalid(&'static str),
    Source(SourceJkrRollingError),
    Mechanics(RigidHertzError),
    GeometryMismatch,
    StepLimit,
    Unresolved {
        evaluations: usize,
        angle_residual_rad: f64,
    },
    /// The turning orbit could not be shown to remain below rolling yield.
    UnresolvedRollingReversal,
    Audit(Box<RigidAdhesiveRollingReport>),
}
impl From<SourceJkrRollingError> for RigidAdhesiveRollingError {
    fn from(x: SourceJkrRollingError) -> Self {
        Self::Source(x)
    }
}
impl From<RigidHertzError> for RigidAdhesiveRollingError {
    fn from(x: RigidHertzError) -> Self {
        Self::Mechanics(x)
    }
}
impl From<DynamicsError> for RigidAdhesiveRollingError {
    fn from(x: DynamicsError) -> Self {
        Self::Mechanics(x.into())
    }
}
impl From<RigidFrictionError> for RigidAdhesiveRollingError {
    fn from(x: RigidFrictionError) -> Self {
        Self::Mechanics(x.into())
    }
}
impl core::fmt::Display for RigidAdhesiveRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native adhesive rolling: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveRollingError {}

impl RigidAdhesiveRollingContact<'_> {
    /// Start unstrained rolling on this exact pair and accepted normal source.
    pub fn initial_state(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        controls: RigidAdhesiveRollingControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveRollingState, RigidAdhesiveRollingError> {
        let native = moment::capture::<SourceJkrRollingContact>(
            self.binding,
            self.frame,
            vector(self.frame.longitudinal()),
            bodies,
            properties,
            controls,
            gate,
        )?;
        Ok(RigidAdhesiveRollingState {
            source: self.source.initial_state(),
            geometry: native.geometry,
        })
    }
}

/// Prepare original-source rolling and actual native recoil without publication.
pub fn prepare_rigid_adhesive_rolling(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingCandidate, RigidAdhesiveRollingError> {
    prepare_rigid_adhesive_rolling_with_wrenches(
        bodies,
        properties,
        state,
        contact,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        gate,
    )
}

/// Stage rolling driven by held actual COM forces and body-frame torques.
/// Loads are in the original binding's A/B order. Full inertia enters both the
/// contact solve and external work; no intermediate free-body state is accepted.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_rolling_with_wrenches(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingCandidate, RigidAdhesiveRollingError> {
    let native = moment::capture::<SourceJkrRollingContact>(
        contact.binding,
        contact.frame,
        vector(contact.frame.longitudinal()),
        bodies,
        properties,
        controls,
        gate,
    )?;
    let candidate = moment::prepare(
        native,
        state.geometry,
        contact.source,
        &state.source,
        wrenches,
        dt_s,
        controls,
        gate,
    )?;
    Ok(RigidAdhesiveRollingCandidate {
        bodies_after: candidate.bodies_after,
        state_after: RigidAdhesiveRollingState {
            source: candidate.source_after,
            geometry: state.geometry,
        },
        report: candidate.report.into(),
    })
}

/// Publish both native momenta and original rolling memory only after all gates.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingReport, RigidAdhesiveRollingError> {
    apply_rigid_adhesive_rolling_with_wrenches(
        bodies,
        properties,
        state,
        contact,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        gate,
    )
}

/// Publish the loaded pair and its original rolling memory atomically.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling_with_wrenches(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingReport, RigidAdhesiveRollingError> {
    let candidate = prepare_rigid_adhesive_rolling_with_wrenches(
        bodies, properties, state, contact, wrenches, dt_s, controls, gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *state = candidate.state_after;
    Ok(candidate.report)
}

impl From<MomentReport<SourceJkrRollingContact>> for RigidAdhesiveRollingReport {
    fn from(r: MomentReport<SourceJkrRollingContact>) -> Self {
        Self {
            contact: r.contact,
            angle_increment_rad: r.angle_increment_rad,
            rolling_rate_before_rad_s: r.rate_before_rad_s,
            rolling_rate_after_rad_s: r.rate_after_rad_s,
            rotational_mobility_per_kg_m2: r.rotational_mobility_per_kg_m2,
            angular_impulse_on_a_world_nms: r.angular_impulse_on_a_world_nms,
            external_wrenches: r.external_wrenches,
            body_displacement_world_m: r.body_displacement_world_m,
            angular_displacement_body_rad: r.angular_displacement_body_rad,
            kinetic_energy_change_j: r.kinetic_energy_change_j,
            midpoint_impulse_work_j: r.midpoint_impulse_work_j,
            external_work_j: r.external_work_j,
            work_energy_residual_j: r.work_energy_residual_j,
            total_energy_residual_j: r.total_energy_residual_j,
            angle_residual_rad: r.angle_residual_rad,
            angular_velocity_residual_rad_s: r.angular_velocity_residual_rad_s,
            body_angular_impulse_residual_nms: r.body_angular_impulse_residual_nms,
            body_linear_impulse_residual_ns: r.body_linear_impulse_residual_ns,
            linear_momentum_change_ns: r.linear_momentum_change_ns,
            angular_momentum_change_nms: r.angular_momentum_change_nms,
            linear_momentum_residual_ns: r.linear_momentum_residual_ns,
            angular_momentum_residual_nms: r.angular_momentum_residual_nms,
            evaluations: r.evaluations,
        }
    }
}
impl From<MomentError<SourceJkrRollingContact>> for RigidAdhesiveRollingError {
    fn from(error: MomentError<SourceJkrRollingContact>) -> Self {
        match error {
            MomentError::Invalid(field) => Self::Invalid(field),
            MomentError::Source(e) => Self::Source(e),
            MomentError::Mechanics(e) => Self::Mechanics(e),
            MomentError::GeometryMismatch => Self::GeometryMismatch,
            MomentError::StepLimit => Self::StepLimit,
            MomentError::Unresolved {
                evaluations,
                angle_residual_rad,
            } => Self::Unresolved {
                evaluations,
                angle_residual_rad,
            },
            MomentError::Reversal => Self::UnresolvedRollingReversal,
            MomentError::Audit(report) => Self::Audit(Box::new((*report).into())),
        }
    }
}

//! Source-owned bonded torsion on a finite rigid pair at held geometry.
//!
//! Relative angular velocity projected onto the original contact normal drives
//! elastic twist and calibrated irreversible yield. The actual two-body inverse
//! inertia supplies the scalar recoil; equal/opposite free couples change both
//! original angular momenta. The shared rolling/torsion phase independently
//! audits native kinetic change against stored twist and irreversible loss.
//!
//! This is one fixed-normal angular impulse phase. It does not impose normal
//! or tangential force reactions, advance poses, transport history, evolve area,
//! solve local Coulomb partial slip, combine sliding capacities, or allocate
//! heat. The source's mode-III/adhesion ceiling screens the frozen-area model.
//! Work integration is exact on each monotone constitutive path; inertial
//! midpoint motion remains a time discretization. Unresolved sign reversals
//! refuse unless the guarded scalar-mode energy excludes hidden yield.
//! Held COM forces and body-frame torques drive the same native solve and
//! contribute their actual work. A driven reversal requires interval refinement.
//! [`normal`] continues a loaded spring through a prescribed normal-area change.

pub mod normal;

use super::super::hertz::RigidHertzError;
use super::super::{RigidContactBinding, RigidFrictionError, poll, vector};
use super::moment::{self, HeldGeometry, MomentError, MomentReport};
use crate::{MassProperties, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::torsion::{
    SourceJkrTorsionContact, SourceJkrTorsionError, SourceJkrTorsionState, SourceJkrTorsionStep,
};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;

/// Physical-step and native work/momentum gates shared with bonded rolling.
pub use super::moment::RigidAdhesiveMomentControls as RigidAdhesiveTorsionControls;

/// Original coincident material markers, frame and source; normal is the twist axis.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveTorsionContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceJkrTorsionContact,
}

/// Source memory bound to the original pair, frame, poses and mass properties.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionState {
    source: SourceJkrTorsionState,
    geometry: HeldGeometry,
}
impl RigidAdhesiveTorsionState {
    pub const fn source_state(&self) -> SourceJkrTorsionState {
        self.source
    }
    pub fn elastic_twist_rad(&self) -> f64 {
        self.source.elastic_twist_rad()
    }
    pub fn accumulated_angle_rad(&self) -> f64 {
        self.source.accumulated_angle_rad()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionReport {
    pub contact: SourceJkrTorsionStep,
    pub angle_increment_rad: f64,
    pub twisting_rate_before_rad_s: f64,
    pub twisting_rate_after_rad_s: f64,
    /// Normal dot (I_A^-1 + I_B^-1) normal, using actual world inertias [1/(kg m2)].
    pub rotational_mobility_per_kg_m2: f64,
    /// Contact free-couple impulse only, excluding the external body torques.
    pub angular_impulse_on_a_world_nms: Vec3,
    /// Held world COM forces and body-frame torques in the binding's A/B order.
    pub external_wrenches: [Wrench; 2],
    /// Actual midpoint COM travel, in A/B order; poses remain held.
    pub body_displacement_world_m: [Vec3; 2],
    /// Integrated held-frame angular velocities, not finite orientation updates.
    pub angular_displacement_body_rad: [Vec3; 2],
    /// Actual momentum-based kinetic change; large background totals are not subtracted.
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

/// Complete staged proposal in the original A/B binding order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidAdhesiveTorsionState,
    pub report: RigidAdhesiveTorsionReport,
}

#[derive(Debug)]
pub enum RigidAdhesiveTorsionError {
    Invalid(&'static str),
    Source(SourceJkrTorsionError),
    Mechanics(RigidHertzError),
    GeometryMismatch,
    StepLimit,
    Unresolved {
        evaluations: usize,
        angle_residual_rad: f64,
    },
    /// Refine an interval whose turning orbit might cross the twist-yield limit.
    UnresolvedTorsionReversal,
    Audit(Box<RigidAdhesiveTorsionReport>),
}
impl core::fmt::Display for RigidAdhesiveTorsionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native adhesive torsion: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveTorsionError {}
impl From<RigidFrictionError> for RigidAdhesiveTorsionError {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}

impl RigidAdhesiveTorsionContact<'_> {
    /// Start unstrained twist on this exact pair and accepted original normal area.
    pub fn initial_state(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        controls: RigidAdhesiveTorsionControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTorsionState, RigidAdhesiveTorsionError> {
        let native = moment::capture::<SourceJkrTorsionContact>(
            self.binding,
            self.frame,
            vector(self.frame.normal()),
            bodies,
            properties,
            controls,
            gate,
        )?;
        Ok(RigidAdhesiveTorsionState {
            source: self.source.initial_state(),
            geometry: native.geometry,
        })
    }
}

/// Prepare original-source twist and actual finite-body recoil without publication.
pub fn prepare_rigid_adhesive_torsion(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveTorsionState,
    contact: RigidAdhesiveTorsionContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionCandidate, RigidAdhesiveTorsionError> {
    prepare_rigid_adhesive_torsion_with_wrenches(
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

/// Stage bonded twist driven by held actual COM forces and body-frame torques.
/// Full native inertia maps applied loads into the original contact coordinate;
/// all impulses are summed before either native endpoint is rounded.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_torsion_with_wrenches(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveTorsionState,
    contact: RigidAdhesiveTorsionContact<'_>,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionCandidate, RigidAdhesiveTorsionError> {
    let native = moment::capture::<SourceJkrTorsionContact>(
        contact.binding,
        contact.frame,
        vector(contact.frame.normal()),
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
    Ok(RigidAdhesiveTorsionCandidate {
        bodies_after: candidate.bodies_after,
        state_after: RigidAdhesiveTorsionState {
            source: candidate.source_after,
            geometry: state.geometry,
        },
        report: candidate.report.into(),
    })
}

/// Atomically publish both native momenta and original torsional history.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_torsion(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveTorsionState,
    contact: RigidAdhesiveTorsionContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionReport, RigidAdhesiveTorsionError> {
    apply_rigid_adhesive_torsion_with_wrenches(
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

/// Publish the driven native momenta and source-owned twist history together.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_torsion_with_wrenches(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveTorsionState,
    contact: RigidAdhesiveTorsionContact<'_>,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionReport, RigidAdhesiveTorsionError> {
    let candidate = prepare_rigid_adhesive_torsion_with_wrenches(
        bodies, properties, state, contact, wrenches, dt_s, controls, gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *state = candidate.state_after;
    Ok(candidate.report)
}

impl From<MomentReport<SourceJkrTorsionContact>> for RigidAdhesiveTorsionReport {
    fn from(r: MomentReport<SourceJkrTorsionContact>) -> Self {
        Self {
            contact: r.contact,
            angle_increment_rad: r.angle_increment_rad,
            twisting_rate_before_rad_s: r.rate_before_rad_s,
            twisting_rate_after_rad_s: r.rate_after_rad_s,
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
impl From<MomentError<SourceJkrTorsionContact>> for RigidAdhesiveTorsionError {
    fn from(error: MomentError<SourceJkrTorsionContact>) -> Self {
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
            MomentError::Reversal => Self::UnresolvedTorsionReversal,
            MomentError::Audit(report) => Self::Audit(Box::new((*report).into())),
        }
    }
}

//! Native recoil under the original local-speed Hertz pressure integral.
//!
//! Each constitutive stage evaluates the original speed curve across the
//! actual translating/spinning patch. Preparation, impulse publication and
//! independent momentum/work audits share the existing native Hertz owners.
//! No centre-speed effective coefficient or substitute Coulomb dynamics is
//! introduced. The law reports its admitted speed domain and force bound.
//!
//! Geometry, load, radius and inertia stay fixed; there is no external drive.
//! This follows passive motion and does not select a driven static reaction.
//! Source applicability and temperature/history evolution remain with their
//! owners. The original spatial quadrature and physical-time errors require
//! separate refinement; neither a small work residual nor a stopping-tail
//! bound certifies the entire resolved trajectory.

/// Original local-speed friction curves under held external force and torque.
pub mod driven;

use super::trajectory::PreparedRigidHertzTrajectory;
use super::transient::{
    HertzTransientControls, PreparedRigidHertzTransient, RigidHertzTransientError,
    RigidHertzTransientReport,
};
use super::*;
use fs_tribo::sliding_patch::contact::relaxation::HertzRelaxationInput;
use fs_tribo::sliding_patch::rate_dependent::transient::HertzRateTransientError;
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

/// An original pressure footprint and explicit original sliding curve.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRateContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub law: &'a HertzFrictionPatch,
}

/// Native mechanics plus the actual law's admitted finite-speed bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRateTransientReport {
    pub mechanical: RigidHertzTransientReport,
    pub local_speed_bound_m_s: f64,
    pub kinetic_coefficient_bounds: [f64; 2],
    pub maximum_scaled_wrench_norm_n: f64,
}

/// Original A/B endpoints staged for joint publication by a coupled owner.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRateTransientCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzRateTransientReport,
}

/// Original geometry and full native response retained for repeated law trials.
#[derive(Clone, Debug)]
pub struct PreparedRigidHertzRateTransient {
    mechanical: PreparedRigidHertzTransient,
    quadrature: HertzFrictionQuadrature,
}

#[derive(Debug)]
pub enum RigidHertzRateTransientError {
    Mechanical(RigidHertzTransientError),
    Law(HertzRateTransientError),
}
impl From<RigidHertzTransientError> for RigidHertzRateTransientError {
    fn from(error: RigidHertzTransientError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<RigidHertzError> for RigidHertzRateTransientError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanical(error.into())
    }
}
impl From<DynamicsError> for RigidHertzRateTransientError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanical(error.into())
    }
}
impl From<HertzRateTransientError> for RigidHertzRateTransientError {
    fn from(error: HertzRateTransientError) -> Self {
        Self::Law(error)
    }
}
impl core::fmt::Display for RigidHertzRateTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rate-dependent Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for RigidHertzRateTransientError {}

impl PreparedRigidHertzRateTransient {
    /// Capture real bodies and geometry without evaluating a substitute law.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzRateContact<'_>,
        dt_s: f64,
        controls: RigidHertzControls,
        quadrature: HertzFrictionQuadrature,
        transient: HertzTransientControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRateTransientError> {
        let pair = PreparedRigidHertzTrajectory::from_geometry(
            bodies,
            properties,
            contact.binding,
            contact.frame,
            contact.law.normal_force_n(),
            contact.law.contact_radius_m(),
            dt_s,
            controls,
            gate,
        )?;
        Ok(Self {
            mechanical: PreparedRigidHertzTransient::from_pair(pair, transient, gate)?,
            quadrature,
        })
    }

    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.mechanical.contact_points_world()
    }

    /// Resolve the original rate law from the same accepted native momenta.
    /// Normal load and radius must match the captured contact. Every native
    /// body/energy audit is shared with the constant-coefficient trajectory.
    pub fn advance(
        &self,
        law: &HertzFrictionPatch,
        gate: &CancelGate,
    ) -> Result<RigidHertzRateTransientCandidate, RigidHertzRateTransientError> {
        poll(gate).map_err(RigidHertzError::from)?;
        let pair = &self.mechanical.pair;
        pair.admit_geometry(law.normal_force_n(), law.contact_radius_m())?;
        let velocity = relative_point(pair.points_before)?;
        let spin = relative_spin(pair.points_before);
        let resolved = law.advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: [velocity.dot(pair.axes[0]), velocity.dot(pair.axes[1])],
                spin_rate_rad_s: spin.dot(pair.axes[2]),
                scaled_mobility_per_kg: pair.mobility,
                dt_s: pair.dt_s,
            },
            pair.controls.contact,
            self.quadrature,
            self.mechanical.transient,
            || gate.is_requested(),
        )?;
        let native = self.mechanical.finish_path(
            resolved.trajectory,
            resolved.maximum_scaled_wrench_norm_n,
            gate,
        )?;
        Ok(RigidHertzRateTransientCandidate {
            bodies_after: native.bodies_after,
            report: RigidHertzRateTransientReport {
                mechanical: native.report,
                local_speed_bound_m_s: resolved.local_speed_bound_m_s,
                kinetic_coefficient_bounds: resolved.kinetic_coefficient_bounds,
                maximum_scaled_wrench_norm_n: resolved.maximum_scaled_wrench_norm_n,
            },
        })
    }
}

/// Prepare one complete original rate-law candidate without accepted mutation.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_hertz_rate_transient(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    quadrature: HertzFrictionQuadrature,
    transient: HertzTransientControls,
    gate: &CancelGate,
) -> Result<RigidHertzRateTransientCandidate, RigidHertzRateTransientError> {
    PreparedRigidHertzRateTransient::new(
        bodies, properties, contact, dt_s, controls, quadrature, transient, gate,
    )?
    .advance(contact.law, gate)
}

/// Publish both original body momenta after the final cooperative checkpoint.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_hertz_rate_transient(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateContact<'_>,
    dt_s: f64,
    controls: RigidHertzControls,
    quadrature: HertzFrictionQuadrature,
    transient: HertzTransientControls,
    gate: &CancelGate,
) -> Result<RigidHertzRateTransientReport, RigidHertzRateTransientError> {
    let candidate = prepare_rigid_hertz_rate_transient(
        bodies, properties, contact, dt_s, controls, quadrature, transient, gate,
    )?;
    poll(gate).map_err(RigidHertzError::from)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    Ok(candidate.report)
}

//! Original local-speed Hertz friction under actual external body loads.
//!
//! The physical initial-value trajectory follows the original pressure-resolved
//! curve. The existing driven native pair owns external acceleration, one summed
//! impulse per body, complementary source work and independent momentum/energy
//! admission. No centre-speed coefficient or alternative Coulomb dynamics is
//! inserted for a nonlinear curve.
//!
//! Source parameters, geometry, inertia, normal load, radius and external loads
//! stay fixed within the phase. Spatial pressure and adaptive physical-time
//! errors remain separate. An unresolved driven arrest is an explicit original
//! law refusal; the adapter never clips rates or selects an algebraic branch.
//! Admitted affine invariant paths retain their closed-form stopping interval
//! and event exactness through native impulse, work and atomic-mode publication.

use super::super::driven::{
    DrivenRigidPair, HertzDrivenMode, RigidHertzDrivenError, RigidHertzDrivenReport,
};
use super::*;
use crate::Wrench;
use fs_tribo::sliding_patch::rate_dependent::driven::HertzRateDrivenError;

pub use fs_tribo::sliding_patch::rate_dependent::driven::HertzRateDrivenControls;

/// Actual native motion and the original curve's admitted numerical domain.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRateDrivenReport {
    pub mechanical: RigidHertzDrivenReport,
    pub local_speed_bound_m_s: f64,
    pub kinetic_coefficient_bounds: [f64; 2],
    /// Includes resting static reactions and the complete kinetic trajectory.
    pub maximum_scaled_wrench_norm_n: f64,
    /// False for analytic invariant paths and the delegated Coulomb method.
    /// `exact_invariant` distinguishes the closed-form solution from Coulomb's
    /// curved time stepping, whose maximum step remains a refinement control.
    pub adaptive_error_estimates_available: bool,
    /// Original-law closed-form evolution on an admitted invariant direction.
    pub exact_invariant: bool,
    /// Original law's first arrest interval, relative to this phase [s].
    /// A point interval denotes an analytic event, not an adaptive estimate.
    /// It does not certify a floating-point enclosure of the physical clock.
    pub stopping_time_bracket_s: Option<[f64; 2]>,
    pub estimated_scaled_rate_error_m_s: f64,
    pub estimated_scaled_impulse_error_ns: f64,
    pub estimated_scaled_travel_error_m: f64,
    pub estimated_heat_error_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzRateDrivenCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub report: RigidHertzRateDrivenReport,
}

#[derive(Clone, Debug)]
pub struct PreparedRigidHertzRateDriven {
    mechanical: DrivenRigidPair,
    quadrature: HertzFrictionQuadrature,
    driven: HertzRateDrivenControls,
}

#[derive(Debug)]
pub enum RigidHertzRateDrivenError {
    Mechanical(RigidHertzDrivenError),
    Law(HertzRateDrivenError),
}
impl From<RigidHertzDrivenError> for RigidHertzRateDrivenError {
    fn from(error: RigidHertzDrivenError) -> Self {
        Self::Mechanical(error)
    }
}
impl From<RigidHertzError> for RigidHertzRateDrivenError {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanical(error.into())
    }
}
impl From<HertzRateDrivenError> for RigidHertzRateDrivenError {
    fn from(error: HertzRateDrivenError) -> Self {
        Self::Law(error)
    }
}
impl core::fmt::Display for RigidHertzRateDrivenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native driven rate-dependent Hertz trajectory: {self:?}")
    }
}
impl std::error::Error for RigidHertzRateDrivenError {}

impl PreparedRigidHertzRateDriven {
    /// Capture original geometry and held external loads once, without a
    /// replacement friction law. Every coefficient trial uses these same bodies.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contact: RigidHertzRateContact<'_>,
        external: [Wrench; 2],
        mode: HertzDrivenMode,
        dt_s: f64,
        controls: RigidHertzControls,
        quadrature: HertzFrictionQuadrature,
        driven: HertzRateDrivenControls,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRateDrivenError> {
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
            mechanical: DrivenRigidPair::new(pair, external, mode, gate)?,
            quadrature,
            driven,
        })
    }

    #[must_use]
    pub fn contact_points_world(&self) -> [Vec3; 2] {
        self.mechanical.contact_points_world()
    }

    /// Advance the original local-speed law from the captured driven state.
    pub fn advance(
        &self,
        law: &HertzFrictionPatch,
        gate: &CancelGate,
    ) -> Result<RigidHertzRateDrivenCandidate, RigidHertzRateDrivenError> {
        poll(gate).map_err(RigidHertzError::from)?;
        let pair = &self.mechanical.pair;
        pair.admit_geometry(law.normal_force_n(), law.contact_radius_m())?;
        let resolved = law.advance_driven(
            &self.mechanical.input()?,
            pair.controls.contact,
            self.quadrature,
            self.driven,
            || gate.is_requested(),
        )?;
        let native = self.mechanical.finish_path(
            resolved.trajectory,
            resolved.maximum_scaled_wrench_norm_n,
            gate,
        )?;
        Ok(RigidHertzRateDrivenCandidate {
            bodies_after: native.bodies_after,
            report: RigidHertzRateDrivenReport {
                mechanical: native.report,
                local_speed_bound_m_s: resolved.local_speed_bound_m_s,
                kinetic_coefficient_bounds: resolved.kinetic_coefficient_bounds,
                maximum_scaled_wrench_norm_n: resolved.maximum_scaled_wrench_norm_n,
                adaptive_error_estimates_available: resolved.adaptive_error_estimates_available,
                exact_invariant: resolved.exact_invariant,
                stopping_time_bracket_s: resolved.stopping_time_bracket_s,
                estimated_scaled_rate_error_m_s: resolved.estimated_scaled_rate_error_m_s,
                estimated_scaled_impulse_error_ns: resolved.estimated_scaled_impulse_error_ns,
                estimated_scaled_travel_error_m: resolved.estimated_scaled_travel_error_m,
                estimated_heat_error_j: resolved.estimated_heat_error_j,
            },
        })
    }
}

/// Stage the complete original pair, mode and physical accounts without mutation.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_hertz_rate_driven(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateContact<'_>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    dt_s: f64,
    controls: RigidHertzControls,
    quadrature: HertzFrictionQuadrature,
    driven: HertzRateDrivenControls,
    gate: &CancelGate,
) -> Result<RigidHertzRateDrivenCandidate, RigidHertzRateDrivenError> {
    PreparedRigidHertzRateDriven::new(
        bodies, properties, contact, external, mode, dt_s, controls, quadrature, driven, gate,
    )?
    .advance(contact.law, gate)
}

/// Publish both original momenta and the retained mode after the final gate.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_hertz_rate_driven(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateContact<'_>,
    external: [Wrench; 2],
    mode: &mut HertzDrivenMode,
    dt_s: f64,
    controls: RigidHertzControls,
    quadrature: HertzFrictionQuadrature,
    driven: HertzRateDrivenControls,
    gate: &CancelGate,
) -> Result<RigidHertzRateDrivenReport, RigidHertzRateDrivenError> {
    let candidate = prepare_rigid_hertz_rate_driven(
        bodies, properties, contact, external, *mode, dt_s, controls, quadrature, driven, gate,
    )?;
    poll(gate).map_err(RigidHertzError::from)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *mode = candidate.report.mechanical.driven.mode;
    Ok(candidate.report)
}

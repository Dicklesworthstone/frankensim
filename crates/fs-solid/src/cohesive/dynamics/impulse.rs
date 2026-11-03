//! Instantaneous generalized momentum transfer on the original reduced mass.
//!
//! The caller supplies an already resolved impulse, in the same free-coordinate
//! ordering as the retained solid. This operation does not infer restitution,
//! a collision partner, or a new contact topology. Displacement, damage, bond
//! viscosity, friction/contact bindings, physical time and step count survive
//! unchanged. The following finite-duration step owns constitutive evolution
//! and reconsiders affected Coulomb occupancy at its actual incoming velocity.

use fs_exec::Cx;
use fs_sparse::precond::Precond;

use super::CohesiveMotion;

mod velocity;
pub use velocity::{ImpulseControls, ImpulseError, ImpulseReport};

/// Accepted instantaneous candidate and its incremental mechanical work.
/// Add `report.applied_impulse_work_j` to the caller's trajectory work ledger.
pub struct ImpulseStep<'m, 'a> {
    pub motion: CohesiveMotion<'m, 'a>,
    pub report: ImpulseReport,
}

impl<'m, 'a> CohesiveMotion<'m, 'a> {
    /// Apply a finite generalized impulse without rebuilding the motion state.
    ///
    /// The original mass and supplied preconditioner must be unchanged SPD
    /// operators. The event changes only velocity; no finite-time fracture,
    /// viscous or frictional work is fabricated at zero duration. Signed impulse
    /// work may remove energy. A physical impact/contact owner must account for
    /// the other body and any restitution loss separately.
    pub fn impulse_candidate<P: Precond>(
        &self,
        cx: &Cx<'_>,
        generalized_impulse_ns: &[f64],
        preconditioner: &P,
        controls: ImpulseControls,
    ) -> Result<ImpulseStep<'m, 'a>, ImpulseError> {
        self.impulse_candidate_controlled(
            generalized_impulse_ns,
            preconditioner,
            controls,
            &mut || cx.checkpoint().is_err(),
        )
    }

    /// Identical event with caller-owned cancellation. A refusal leaves the
    /// accepted motion intact. Polls between vector batches, around mass
    /// applications and within checked CG; opaque operator/preconditioner
    /// application latency remains their owners' responsibility.
    pub fn impulse_candidate_controlled<P: Precond, C: FnMut() -> bool>(
        &self,
        generalized_impulse_ns: &[f64],
        preconditioner: &P,
        controls: ImpulseControls,
        cancelled: &mut C,
    ) -> Result<ImpulseStep<'m, 'a>, ImpulseError> {
        let (velocity_m_s, report) = velocity::velocity_candidate(
            self.mass,
            &self.velocity_m_s,
            generalized_impulse_ns,
            preconditioner,
            self.config.maximum_velocity_m_s,
            controls,
            cancelled,
        )?;
        // Cloning the accepted owner, rather than re-adopting its material
        // state, preserves every retained binding and its irreversible state.
        let mut motion = self.clone();
        if let Some(binding) = &self.coulomb {
            motion.coulomb = Some(binding.after_impulse(self, &velocity_m_s, cancelled)?);
        }
        motion.velocity_m_s = velocity_m_s;
        if cancelled() {
            return Err(ImpulseError::Cancelled);
        }
        Ok(ImpulseStep { motion, report })
    }
}

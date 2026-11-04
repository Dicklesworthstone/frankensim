//! Prescribed-temperature continuation of the original normal-motion owner.
//!
//! The source refresh owns the temperature queries, stable branch and signed
//! change of contact potential. At fixed approach this is parameter work, not
//! an impulse, viscous loss, snap release or heat. The transferred token retains
//! the entire mechanical chronology and numerical policy; it is never rebuilt
//! as a fresh initial condition. The original immutable refresh must outlive
//! the returned driver, which borrows its successor source directly.

use super::{
    CancelGate, ContentHash, JkrBranch, JkrNormalMotion, JkrNormalMotionError,
    JkrNormalMotionState, finite, poll,
};
use crate::normal_patch::adhesive::{
    JkrError, JkrResponse,
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};

/// Original potential changes once; inertial and support storage stay fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalTemperatureReport {
    pub original_model_identity: ContentHash,
    pub refreshed_model_identity: ContentHash,
    pub before: JkrResponse,
    pub after: JkrResponse,
    /// Signed supply from the explicitly prescribed temperature driver [J].
    /// This is neither a caloric transfer nor an entropy-production estimate.
    pub parameter_work_j: f64,
    /// Actual represented mechanical storage change at unchanged coordinates [J].
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

/// Complete source, state and work proposal. Accepted owners remain untouched.
pub struct JkrNormalTemperatureContinuation<'a> {
    previous: JkrNormalMotionState,
    motion: JkrNormalMotion<'a>,
    state: JkrNormalMotionState,
    report: JkrNormalTemperatureReport,
}
impl<'a> JkrNormalTemperatureContinuation<'a> {
    #[must_use]
    pub const fn previous_state(&self) -> JkrNormalMotionState {
        self.previous
    }
    #[must_use]
    pub const fn motion(&self) -> &JkrNormalMotion<'a> {
        &self.motion
    }
    #[must_use]
    pub const fn state(&self) -> JkrNormalMotionState {
        self.state
    }
    #[must_use]
    pub const fn report(&self) -> &JkrNormalTemperatureReport {
        &self.report
    }
    /// Publish only after the native coordinates and any retained heat history
    /// admit this same source successor and its explicit parameter-work port.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        JkrNormalMotion<'a>,
        JkrNormalMotionState,
        JkrNormalTemperatureReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl JkrNormalMotion<'_> {
    /// Retain normal inertia, branch and chronology through an original-source
    /// temperature change at exactly held indentation.
    ///
    /// `previous` must belong to this complete mass/source/support/policy owner,
    /// and `refresh` must start from that exact contact state. Its original
    /// immutable source controls the refreshed material claims and geometry.
    /// Velocity, time, step count, mass, spring, anchor, dashpot and every motion
    /// budget are preserved. The successor rechecks all source and motion limits.
    /// A newly marginal attached branch refuses: this split cannot create an
    /// attachment or detachment event, even at zero duration.
    ///
    /// The actual contact-potential change is the sole parameter-work supply.
    /// Kinetic and support storage do not change, and no branch release or
    /// damping work occurs. An enclosing caloric coupling needs its own contact
    /// entropy/internal-energy law; this operation never relabels work as heat.
    /// All source, history, force, range, energy and cancellation failures return
    /// without changing either original owner or token.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: JkrNormalMotionState,
        refresh: &'new SourceJkrThermalRefresh,
        gate: &CancelGate,
    ) -> Result<JkrNormalTemperatureContinuation<'new>, JkrNormalMotionError> {
        poll(gate)?;
        let before = self.response(previous)?;
        if refresh.original_model_identity() != self.normal.identity()
            || refresh.before() != &before
            || refresh.before().state() != previous.contact
        {
            return Err(JkrNormalMotionError::StateMismatch);
        }
        if refresh.control() != AdhesiveThermalControl::HeldIndentation
            || refresh.mechanical_work_j() != 0.0
            || refresh.response().indentation_m.to_bits() != before.indentation_m.to_bits()
            || refresh.response().branch != before.branch
            || refresh.source().effective_radius_m() != self.normal.effective_radius_m()
            || refresh.source().interaction_range_m() != self.normal.interaction_range_m()
            || refresh.source().envelope() != self.normal.envelope()
        {
            return Err(JkrNormalMotionError::Invalid(
                "original held-indentation temperature successor and unchanged geometry required",
            ));
        }
        let motion = JkrNormalMotion::new_with_damped_support(
            refresh.source(),
            self.mass_kg,
            self.budget,
            self.support.stiffness,
            self.support.anchor,
            self.damping_n_s_per_m,
        )?;
        let contact = refresh.response().state();
        if contact.branch() == JkrBranch::Attached
            && contact.indentation_m() <= motion.normal.minimum_attached_indentation_m()
        {
            return Err(JkrError::SingularTangent.into());
        }
        let state = JkrNormalMotionState {
            identity: motion.identity,
            contact,
            velocity_m_s: previous.velocity_m_s,
            time_s: previous.time_s,
            steps: previous.steps,
        };
        let after = motion.response(state)?;
        if refresh.response() != &after {
            return Err(JkrNormalMotionError::StateMismatch);
        }
        let kinetic = finite(
            0.5 * self.mass_kg * previous.velocity_m_s * previous.velocity_m_s,
            "temperature continuation kinetic energy",
        )?;
        let support = self.support.energy(before.indentation_m)?;
        // Check representability without subtracting common kinetic/support
        // datums: only the original contact potential changes in this split.
        finite(
            kinetic + before.reversible_energy_j + support,
            "original temperature continuation energy",
        )?;
        finite(
            kinetic + after.reversible_energy_j + support,
            "refreshed temperature continuation energy",
        )?;
        let change = finite(
            after.reversible_energy_j - before.reversible_energy_j,
            "temperature continuation storage change",
        )?;
        let parameter_work_j = refresh.parameter_work_j();
        let residual = finite(
            change - parameter_work_j,
            "temperature continuation energy balance",
        )?;
        let scale = finite(
            before.reversible_energy_j.abs()
                + after.reversible_energy_j.abs()
                + parameter_work_j.abs(),
            "temperature continuation contact scale",
        )?;
        let tolerance = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "temperature continuation tolerance",
        )?;
        for error in [
            residual,
            finite(
                change - refresh.potential_change_j(),
                "original temperature refresh storage balance",
            )?,
        ] {
            if error.abs() > tolerance {
                return Err(JkrNormalMotionError::EnergyBalance {
                    residual_j: error,
                    tolerance_j: tolerance,
                });
            }
        }
        poll(gate)?;
        Ok(JkrNormalTemperatureContinuation {
            previous,
            motion,
            state,
            report: JkrNormalTemperatureReport {
                original_model_identity: self.normal.identity(),
                refreshed_model_identity: refresh.source().identity(),
                before,
                after,
                parameter_work_j,
                storage_change_j: change,
                energy_residual_j: residual,
                energy_tolerance_j: tolerance,
            },
        })
    }
}

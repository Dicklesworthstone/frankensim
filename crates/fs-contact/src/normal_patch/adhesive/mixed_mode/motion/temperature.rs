//! Original thermal sources at held normal AND elastic shear coordinates.
//!
//! The sealed normal refresh authenticates material cards, selected claims,
//! temperature queries and geometry. Its pure-normal work is not the work of a
//! sheared contact. Both mixed-mode endpoints are therefore evaluated under
//! their own original source, and all three storage contributions enter signed
//! parameter work. No mechanical impulse, clock increment, dissipation, caloric
//! transfer or continuous thermal-path stability is inferred.

use super::{JkrMixedModeMotion, JkrMixedModeMotionError, JkrMixedModeMotionState, finite, poll};
use crate::interface_binding::adhesive_shear::SourceAdhesiveShear;
use crate::normal_patch::adhesive::{
    JkrBranch,
    mixed_mode::{JkrMixedModeResponse, SourceJkrMixedMode, path::JkrMixedModePathError},
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Original mixed-mode storage changes; inertia and held coordinates do not.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModeTemperatureReport {
    pub original_model_identity: ContentHash,
    pub refreshed_model_identity: ContentHash,
    pub before: JkrMixedModeResponse,
    pub after: JkrMixedModeResponse,
    pub normal_elastic_change_j: f64,
    pub tangential_elastic_change_j: f64,
    pub surface_energy_change_j: f64,
    /// Signed supply from the prescribed parameter change, not heat or loss.
    pub parameter_work_j: f64,
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub evaluations: usize,
}

/// Complete immutable proposal; retain its parts only after native acceptance.
pub struct JkrMixedModeTemperatureContinuation<'a> {
    previous: JkrMixedModeMotionState,
    motion: JkrMixedModeMotion<'a>,
    state: JkrMixedModeMotionState,
    report: JkrMixedModeTemperatureReport,
}
impl<'a> JkrMixedModeTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> JkrMixedModeMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrMixedModeMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrMixedModeMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrMixedModeTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrMixedModeMotion<'a>,
        JkrMixedModeMotionState,
        JkrMixedModeTemperatureReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl JkrMixedModeMotion<'_> {
    /// Transfer the exact accepted full-stick history to an authenticated
    /// original-source thermal successor at fixed [approach,shear1,shear2].
    ///
    /// Construct `successor` from `refresh.source().clone()` and the unchanged
    /// mixed-mode envelope; construct `shear` with the original shear owner's
    /// `rebind_jkr_temperature`. Both are borrowed by the returned driver.
    /// The sealed refresh must belong to the exact original normal source at
    /// this state's held approach, with the attached branch on both sides.
    /// Independent bulk/interface temperatures and pinned claims are retained.
    ///
    /// Velocity, clock, step count, full mobility and every numerical/physical
    /// budget stay unchanged. Refreshed radius stability, source applicability,
    /// static capacity and force limits are rechecked before a token is issued.
    /// A lost equilibrium or capacity refuses instead of inventing slip or snap.
    /// This endpoint parameter update makes no claim about an intervening
    /// thermal path, thermal expansion, entropy or a caloric heat law.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: JkrMixedModeMotionState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        shear: &'new SourceAdhesiveShear,
        gate: &CancelGate,
    ) -> Result<JkrMixedModeTemperatureContinuation<'new>, JkrMixedModeMotionError> {
        poll(gate)?;
        let normal = self.source.normal_source();
        if previous.identity != self.identity
            || refresh.original_model_identity() != normal.identity()
            || refresh.before().normal_state_identity() != normal.normal_state().identity()
        {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        if refresh.control() != AdhesiveThermalControl::HeldIndentation
            || refresh.mechanical_work_j() != 0.0
            || refresh.before().branch != JkrBranch::Attached
            || refresh.response().branch != JkrBranch::Attached
            || refresh.before().indentation_m.to_bits() != previous.position_m[0].to_bits()
            || refresh.response().indentation_m.to_bits() != previous.position_m[0].to_bits()
            || successor.normal_source() != refresh.source()
            || successor.envelope() != self.source.envelope()
            || refresh.source().effective_radius_m() != normal.effective_radius_m()
            || refresh.source().interaction_range_m() != normal.interaction_range_m()
            || refresh.source().envelope() != normal.envelope()
        {
            return Err(JkrMixedModeMotionError::Invalid(
                "original held-coordinate thermal successor and unchanged geometry required",
            ));
        }
        // Reserve the canonical normal query and both coincident mixed queries;
        // all remaining corner probes share the SAME original work budget.
        if self.budget.maximum_evaluations <= 3 {
            return Err(JkrMixedModeMotionError::Budget);
        }
        let canonical = normal
            .evaluate_attached(previous.position_m[0])
            .map_err(|error| JkrMixedModePathError::Source(error.into()))?;
        if refresh.before() != &canonical {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        poll(gate)?;
        let expected_shear = self
            .shear
            .rebind_jkr_temperature(normal, refresh)
            .map_err(JkrMixedModePathError::from)?;
        if shear != &expected_shear {
            return Err(JkrMixedModeMotionError::StateMismatch);
        }
        let before = self
            .source
            .discrete_path(previous.position_m, previous.position_m, gate)?
            .contact;
        let motion = JkrMixedModeMotion::new(successor, shear, self.mobility, self.budget)?;
        motion.admit_velocity(previous.velocity_m_s)?;
        let bounds = successor.admit_quadratic_path(
            shear,
            [previous.position_m; 3],
            self.budget.maximum_evaluations - 3,
            gate,
        )?;
        if bounds.maximum_force_n > self.budget.maximum_force_n {
            return Err(JkrMixedModeMotionError::Range("refreshed source force"));
        }
        let after = successor
            .discrete_path(previous.position_m, previous.position_m, gate)?
            .contact;
        let state = JkrMixedModeMotionState {
            identity: motion.identity,
            ..previous
        };
        let kinetic = finite(
            0.5 * super::dot(
                previous.velocity_m_s,
                super::mul(self.mass, previous.velocity_m_s),
            ),
            "thermal continuation kinetic energy",
        )?;
        finite(
            kinetic + before.reversible_energy_j,
            "original thermal mechanical storage",
        )?;
        finite(
            kinetic + after.reversible_energy_j,
            "refreshed thermal mechanical storage",
        )?;
        let normal_elastic_change_j =
            after.normal_elastic_energy_j - before.normal_elastic_energy_j;
        let tangential_elastic_change_j =
            after.tangential_elastic_energy_j - before.tangential_elastic_energy_j;
        let surface_energy_change_j = after.surface_energy_j - before.surface_energy_j;
        let parameter_work_j = finite(
            (Dd::from_f64(normal_elastic_change_j)
                + Dd::from_f64(tangential_elastic_change_j)
                + Dd::from_f64(surface_energy_change_j))
            .to_f64(),
            "mixed-mode thermal parameter work",
        )?;
        let storage_change_j = finite(
            after.reversible_energy_j - before.reversible_energy_j,
            "mixed-mode thermal storage change",
        )?;
        let energy_residual_j = finite(
            storage_change_j - parameter_work_j,
            "mixed-mode thermal work balance",
        )?;
        let scale = finite(
            before.normal_elastic_energy_j + before.tangential_elastic_energy_j
                - before.surface_energy_j
                + after.normal_elastic_energy_j
                + after.tangential_elastic_energy_j
                - after.surface_energy_j
                + parameter_work_j.abs(),
            "mixed-mode thermal storage scale",
        )?;
        let energy_tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "mixed-mode thermal work tolerance",
        )?;
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(JkrMixedModeMotionError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(JkrMixedModeTemperatureContinuation {
            previous,
            motion,
            state,
            report: JkrMixedModeTemperatureReport {
                original_model_identity: self.source.identity(),
                refreshed_model_identity: successor.identity(),
                before,
                after,
                normal_elastic_change_j,
                tangential_elastic_change_j,
                surface_energy_change_j,
                parameter_work_j,
                storage_change_j,
                energy_residual_j,
                energy_tolerance_j,
                evaluations: 3 + bounds.evaluations,
            },
        })
    }
}

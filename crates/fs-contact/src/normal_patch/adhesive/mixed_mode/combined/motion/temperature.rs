//! Original thermal sources at held approach, vector shear and twist.
//!
//! The original normal refresh authenticates the independent bulk/interface
//! queries. Its pure-normal work is not the work of this loaded contact. The
//! complete I/II/III potential is therefore evaluated with one relaxed radius
//! under the successor source, retaining all four signed storage changes.
//! This endpoint operation supplies no impulse, time advance, dissipation,
//! caloric heat or stability claim along an unspecified temperature path.

use super::{
    JkrCombinedMotion, JkrCombinedMotionError, JkrCombinedMotionState, JkrCombinedPathError,
    JkrCombinedResponse, SourceJkrMixedMode, finite, poll,
};
use crate::normal_patch::adhesive::{
    JkrBranch,
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Signed original normal, shear, twist and surface storage changes. Coordinates,
/// rates and native mobility remain fixed, so mechanical actuator work is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedTemperatureReport {
    pub original_model_identity: ContentHash,
    pub refreshed_model_identity: ContentHash,
    pub before: JkrCombinedResponse,
    pub after: JkrCombinedResponse,
    pub normal_elastic_change_j: f64,
    pub tangential_elastic_change_j: f64,
    pub torsional_elastic_change_j: f64,
    pub surface_energy_change_j: f64,
    /// Signed supply from the prescribed coefficient changes, not heat or loss.
    pub parameter_work_j: f64,
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Canonical normal witness and one original coincident I/II/III query.
    pub evaluations: usize,
}

/// Original history and authenticated successor. Adopting the proposal requires
/// retaining its source, state and signed parameter-work accounting together.
#[derive(Clone, Debug)]
pub struct JkrCombinedTemperatureContinuation<'a> {
    previous: JkrCombinedMotionState,
    motion: JkrCombinedMotion<'a>,
    state: JkrCombinedMotionState,
    report: JkrCombinedTemperatureReport,
}
impl<'a> JkrCombinedTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> JkrCombinedMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrCombinedMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrCombinedMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrCombinedTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrCombinedMotion<'a>,
        JkrCombinedMotionState,
        JkrCombinedTemperatureReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl JkrCombinedMotion<'_> {
    /// Continue accepted full-stick motion at held [approach,shear1,shear2,twist].
    ///
    /// Construct `successor` from `refresh.source().clone()` and the unchanged
    /// mixed-mode envelope. The sealed refresh must belong to the exact original
    /// normal owner at this state's held approach. Its selected claims and
    /// independent bulk/interface temperatures remain authoritative, while the
    /// complete successor I/II/III potential owns relaxed area and stored work.
    ///
    /// Rates, full mixed-unit mobility, all budgets, time and step count remain
    /// unchanged. One original coincident combined query rechecks strict stable
    /// attachment, geometry/material applicability and every conjugate effort.
    /// Lost equilibrium, changed source/envelope, exhausted work or cancellation
    /// publish nothing. This does not infer slip, a snap, thermal expansion,
    /// entropy, caloric feedback or stability along an unobserved thermal path.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: JkrCombinedMotionState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<JkrCombinedTemperatureContinuation<'new>, JkrCombinedMotionError> {
        poll(gate)?;
        self.check(previous)?;
        let normal = self.source.normal_source();
        if refresh.original_model_identity() != normal.identity()
            || refresh.before().normal_state_identity() != normal.normal_state().identity()
        {
            return Err(JkrCombinedMotionError::StateMismatch);
        }
        let coordinates = previous.coordinates();
        if refresh.control() != AdhesiveThermalControl::HeldIndentation
            || refresh.mechanical_work_j() != 0.0
            || refresh.before().branch != JkrBranch::Attached
            || refresh.response().branch != JkrBranch::Attached
            || refresh.before().indentation_m.to_bits() != coordinates[0].to_bits()
            || refresh.response().indentation_m.to_bits() != coordinates[0].to_bits()
            || successor.normal_source() != refresh.source()
            || successor.envelope() != self.source.envelope()
            || refresh.source().effective_radius_m() != normal.effective_radius_m()
            || refresh.source().interaction_range_m() != normal.interaction_range_m()
            || refresh.source().envelope() != normal.envelope()
        {
            return Err(JkrCombinedMotionError::Invalid(
                "original held-coordinate thermal successor and unchanged geometry required",
            ));
        }
        // The old combined endpoint is already sealed in accepted history.
        // Only the canonical normal witness and new combined endpoint are queried.
        if self.budget.maximum_evaluations < 2 {
            return Err(JkrCombinedMotionError::Budget);
        }
        let canonical = normal
            .evaluate_attached(coordinates[0])
            .map_err(|error| JkrCombinedPathError::Source(error.into()))?;
        if refresh.before() != &canonical {
            return Err(JkrCombinedMotionError::StateMismatch);
        }
        poll(gate)?;
        let motion = JkrCombinedMotion::new(successor, self.mobility, self.budget)?;
        let admitted = motion.initial_state(coordinates, previous.rates, gate)?;
        let before = previous.contact;
        let after = admitted.contact;
        let state = JkrCombinedMotionState {
            identity: motion.identity,
            contact: after,
            ..previous
        };
        let kinetic = motion.kinetic_energy(previous.rates)?;
        finite(
            kinetic + before.reversible_energy_j,
            "original thermal mechanical storage",
        )?;
        finite(
            kinetic + after.reversible_energy_j,
            "refreshed thermal mechanical storage",
        )?;
        let normal_elastic_change_j = finite(
            after.normal_elastic_energy_j - before.normal_elastic_energy_j,
            "thermal normal elastic change",
        )?;
        let tangential_elastic_change_j = finite(
            after.tangential_elastic_energy_j - before.tangential_elastic_energy_j,
            "thermal tangential elastic change",
        )?;
        let torsional_elastic_change_j = finite(
            after.torsional_elastic_energy_j - before.torsional_elastic_energy_j,
            "thermal torsional elastic change",
        )?;
        let surface_energy_change_j = finite(
            after.surface_energy_j - before.surface_energy_j,
            "thermal surface energy change",
        )?;
        let parameter_work_j = finite(
            (Dd::from_f64(normal_elastic_change_j)
                + Dd::from_f64(tangential_elastic_change_j)
                + Dd::from_f64(torsional_elastic_change_j)
                + Dd::from_f64(surface_energy_change_j))
            .to_f64(),
            "combined thermal parameter work",
        )?;
        let storage_change_j = finite(
            after.reversible_energy_j - before.reversible_energy_j,
            "combined thermal storage change",
        )?;
        let energy_residual_j = finite(
            storage_change_j - parameter_work_j,
            "combined thermal work balance",
        )?;
        let scale = finite(
            before.normal_elastic_energy_j
                + before.tangential_elastic_energy_j
                + before.torsional_elastic_energy_j
                - before.surface_energy_j
                + after.normal_elastic_energy_j
                + after.tangential_elastic_energy_j
                + after.torsional_elastic_energy_j
                - after.surface_energy_j
                + parameter_work_j.abs(),
            "combined thermal storage scale",
        )?;
        let energy_tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "combined thermal work tolerance",
        )?;
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(JkrCombinedMotionError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(JkrCombinedTemperatureContinuation {
            previous,
            motion,
            state,
            report: JkrCombinedTemperatureReport {
                original_model_identity: self.source.identity(),
                refreshed_model_identity: successor.identity(),
                before,
                after,
                normal_elastic_change_j,
                tangential_elastic_change_j,
                torsional_elastic_change_j,
                surface_energy_change_j,
                parameter_work_j,
                storage_change_j,
                energy_residual_j,
                energy_tolerance_j,
                evaluations: 2,
            },
        })
    }
}

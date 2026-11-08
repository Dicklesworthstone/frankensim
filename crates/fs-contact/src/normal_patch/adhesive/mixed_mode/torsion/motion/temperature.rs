//! Source-authenticated temperature changes at held approach AND twist.
//!
//! The original normal refresh owns bulk/interface queries and selected claims.
//! Its pure-normal work cannot stand in for the loaded I/III potential. The
//! cached original torsional endpoint and the successor's original torsional
//! response supply all three storage changes, including area relaxation.
//! This endpoint parameter update supplies no mechanical impulse, time advance,
//! caloric heat or entropy law and makes no claim about an intervening thermal
//! path. Losing the admitted attached branch refuses without changing history.

use super::{
    JkrTorsionalMotion, JkrTorsionalMotionError, JkrTorsionalMotionState, JkrTorsionalPathError,
    JkrTorsionalResponse, SourceJkrMixedMode, finite, poll,
};
use crate::normal_patch::adhesive::{
    JkrBranch,
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Signed changes of the original normal/torsional/surface stores. Mechanical
/// coordinates, rates and native mobility do not change, so their work is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalTemperatureReport {
    pub original_model_identity: ContentHash,
    pub refreshed_model_identity: ContentHash,
    pub before: JkrTorsionalResponse,
    pub after: JkrTorsionalResponse,
    pub normal_elastic_change_j: f64,
    pub torsional_elastic_change_j: f64,
    pub surface_energy_change_j: f64,
    /// Supply from prescribed coefficient changes; neither heat nor loss [J].
    pub parameter_work_j: f64,
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Canonical normal-witness query plus one original successor I/III query.
    pub evaluations: usize,
}

/// Immutable complete proposal; a native owner accepts all parts together.
#[derive(Clone, Debug)]
pub struct JkrTorsionalTemperatureContinuation<'a> {
    previous: JkrTorsionalMotionState,
    motion: JkrTorsionalMotion<'a>,
    state: JkrTorsionalMotionState,
    report: JkrTorsionalTemperatureReport,
}
impl<'a> JkrTorsionalTemperatureContinuation<'a> {
    pub const fn previous_state(&self) -> JkrTorsionalMotionState {
        self.previous
    }
    pub const fn motion(&self) -> &JkrTorsionalMotion<'a> {
        &self.motion
    }
    pub const fn state(&self) -> JkrTorsionalMotionState {
        self.state
    }
    pub const fn report(&self) -> &JkrTorsionalTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrTorsionalMotion<'a>,
        JkrTorsionalMotionState,
        JkrTorsionalTemperatureReport,
    ) {
        (self.motion, self.state, self.report)
    }
}

impl JkrTorsionalMotion<'_> {
    /// Continue the exact accepted history through an authenticated original
    /// thermal refresh at fixed [approach,twist].
    ///
    /// Construct the borrowed successor from refresh.source().clone() and the
    /// unchanged mixed-mode envelope. The refresh must be the original normal
    /// owner's held-indentation operation at this state's exact approach; it
    /// retains independent bulk/interface temperatures and original claim
    /// selection. A normal-only refresh does not establish torsional stability,
    /// so the complete successor I/III response is evaluated and admitted.
    ///
    /// Rates, complete mobility, budget, time and step count are unchanged.
    /// Only the source identity and original response change. Refreshed source
    /// or effort-limit failure, exhausted work allowance and cancellation issue
    /// no successor. This operation does not infer thermal expansion, a thermal
    /// path, slip, detachment, kinetic exchange, caloric heat or dissipation.
    pub fn continue_at_temperature<'new>(
        &self,
        previous: JkrTorsionalMotionState,
        refresh: &SourceJkrThermalRefresh,
        successor: &'new SourceJkrMixedMode,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalTemperatureContinuation<'new>, JkrTorsionalMotionError> {
        poll(gate)?;
        self.check(previous)?;
        let normal = self.source.normal_source();
        if refresh.original_model_identity() != normal.identity()
            || refresh.before().normal_state_identity() != normal.normal_state().identity()
        {
            return Err(JkrTorsionalMotionError::StateMismatch);
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
            return Err(JkrTorsionalMotionError::Invalid(
                "original held-coordinate thermal successor and unchanged geometry required",
            ));
        }
        // The accepted original torsional endpoint is already sealed. Only the
        // normal witness and the successor endpoint need constitutive queries.
        if self.budget.maximum_evaluations < 2 {
            return Err(JkrTorsionalMotionError::Budget);
        }
        let canonical = normal
            .evaluate_attached(coordinates[0])
            .map_err(|error| JkrTorsionalPathError::Source(error.into()))?;
        if refresh.before() != &canonical {
            return Err(JkrTorsionalMotionError::StateMismatch);
        }
        poll(gate)?;
        let motion = JkrTorsionalMotion::new(successor, self.mobility, self.budget)?;
        // Reuse original coincidence admission: strict stable attachment,
        // original geometry/material limits, both efforts and finite inertia.
        let admitted = motion.initial_state(coordinates, previous.rates, gate)?;
        let before = previous.contact;
        let after = admitted.contact;
        let state = JkrTorsionalMotionState {
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
                + Dd::from_f64(torsional_elastic_change_j)
                + Dd::from_f64(surface_energy_change_j))
            .to_f64(),
            "torsional thermal parameter work",
        )?;
        let storage_change_j = finite(
            after.reversible_energy_j - before.reversible_energy_j,
            "torsional thermal storage change",
        )?;
        let energy_residual_j = finite(
            storage_change_j - parameter_work_j,
            "torsional thermal work balance",
        )?;
        let scale = finite(
            before.normal_elastic_energy_j + before.torsional_elastic_energy_j
                - before.surface_energy_j
                + after.normal_elastic_energy_j
                + after.torsional_elastic_energy_j
                - after.surface_energy_j
                + parameter_work_j.abs(),
            "torsional thermal storage scale",
        )?;
        let energy_tolerance_j = finite(
            self.budget.energy_absolute_tolerance_j + self.budget.energy_relative_tolerance * scale,
            "torsional thermal work tolerance",
        )?;
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(JkrTorsionalMotionError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(JkrTorsionalTemperatureContinuation {
            previous,
            motion,
            state,
            report: JkrTorsionalTemperatureReport {
                original_model_identity: self.source.identity(),
                refreshed_model_identity: successor.identity(),
                before,
                after,
                normal_elastic_change_j,
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

//! Loaded twist follows original material claims at actual cell temperatures.
//!
//! This same-time attached continuation retains native geometry, full twist
//! memory and every heat ledger. Signed coefficient work stays with its named
//! external driver; caloric enthalpy is unchanged by the refresh.

use super::{JkrTorsionHeat, JkrTorsionHeatStamp};
pub use crate::friction_heat::adhesive::snap::angular_temperature::{
    JkrAngularTemperatureHeatControls,
    JkrAngularTemperatureHeatError as JkrTorsionTemperatureHeatError,
};
use crate::friction_heat::adhesive::snap::angular_temperature::{check, temperatures};
use fs_contact::normal_patch::adhesive::{
    temperature::AdhesiveTemperatures, torsion::temperature::SourceJkrTorsionTemperatureReport,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_mbd::tangential_friction::adhesive::torsion::events::RigidAdhesiveTorsionEventState;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
pub struct JkrTorsionTemperatureHeatReport {
    /// Complete source identities, normal endpoints and separate signed work.
    pub temperature: SourceJkrTorsionTemperatureReport,
    pub temperatures: AdhesiveTemperatures,
    pub parameter_work_source: String,
}

pub struct JkrTorsionTemperatureHeatCandidate<'s> {
    history: JkrTorsionHeat<'s>,
    native_state: RigidAdhesiveTorsionEventState,
    report: JkrTorsionTemperatureHeatReport,
}
impl<'s> JkrTorsionTemperatureHeatCandidate<'s> {
    pub const fn history(&self) -> &JkrTorsionHeat<'s> {
        &self.history
    }
    pub const fn native_state(&self) -> &RigidAdhesiveTorsionEventState {
        &self.native_state
    }
    pub const fn report(&self) -> &JkrTorsionTemperatureHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrTorsionHeat<'s>,
        RigidAdhesiveTorsionEventState,
        JkrTorsionTemperatureHeatReport,
    ) {
        (self.history, self.native_state, self.report)
    }
}

impl<'s> JkrTorsionHeat<'s> {
    /// Query the original bulk materials at the actual two receiving-cell
    /// temperatures and the interface at its explicit prescribed temperature.
    /// Normal approach, total/elastic twist, native geometry and calibration
    /// stay held. The caller owns validity of the unchanged critical angle;
    /// the original small-twist and mode-III applicability gates run again.
    ///
    /// Requires the accepted thermal time and exact native/source history.
    /// This endpoint operation supplies no thermal branch event, entropy or
    /// transport. Its signed normal-plus-torsion parameter work is bounded and
    /// assigned to an external driver; caloric enthalpy is never debited here.
    #[allow(clippy::too_many_arguments)]
    pub fn propose_native_temperature_continuation(
        &self,
        state: &RigidAdhesiveTorsionEventState,
        thermal: &DepositionSimulation<'_>,
        expected: JkrTorsionHeatStamp,
        material_cards: [&MaterialCard; 2],
        interface_temperature_k: f64,
        controls: JkrAngularTemperatureHeatControls<'_>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTorsionTemperatureHeatCandidate<'s>, JkrTorsionTemperatureHeatError> {
        check(gate, &mut checkpoint)?;
        controls.validate()?;
        if expected != self.stamp
            || self.state() != state.source_state()
            || thermal.time_s() != self.accepted_time_s
        {
            return Err(JkrTorsionTemperatureHeatError::HistoryMismatch);
        }
        self.admit_receivers(thermal)?;
        let revision =
            self.stamp
                .revision
                .checked_add(1)
                .ok_or(JkrTorsionTemperatureHeatError::Invalid(
                    "heat history version exhausted",
                ))?;
        let temperatures = temperatures(
            thermal,
            self.receivers.map(|receiver| receiver.cell),
            interface_temperature_k,
        )?;
        let normal = self.source().normal_source().refresh_at_held_indentation(
            self.state().normal_response().state(),
            material_cards,
            temperatures,
        )?;
        let continuation = self
            .source()
            .continue_at_temperature(self.state(), normal)?;
        controls.admit_work(continuation.report().parameter_work_j)?;
        let native_state = state.continue_at_temperature(&continuation, gate)?;
        let (source, source_state, temperature) = continuation.into_parts();
        let mut history = self.clone();
        history.stamp.source = source.normal_source().identity();
        history.source = Arc::new(source);
        history.stamp.state = source_state;
        history.stamp.revision = revision;
        let report = JkrTorsionTemperatureHeatReport {
            temperature,
            temperatures,
            parameter_work_source: controls.parameter_work_source.to_owned(),
        };
        check(gate, &mut checkpoint)?;
        Ok(JkrTorsionTemperatureHeatCandidate {
            history,
            native_state,
            report,
        })
    }

    /// Publish both original histories only after the final caller checkpoint
    /// and cancellation check. Bodies, thermal cells and time are read-only.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_native_temperature_continuation(
        &mut self,
        state: &mut RigidAdhesiveTorsionEventState,
        thermal: &DepositionSimulation<'_>,
        expected: JkrTorsionHeatStamp,
        material_cards: [&MaterialCard; 2],
        interface_temperature_k: f64,
        controls: JkrAngularTemperatureHeatControls<'_>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrTorsionTemperatureHeatReport, JkrTorsionTemperatureHeatError> {
        let candidate = self.propose_native_temperature_continuation(
            state,
            thermal,
            expected,
            material_cards,
            interface_temperature_k,
            controls,
            gate,
            &mut checkpoint,
        )?;
        check(gate, &mut checkpoint)?;
        let (history, native_state, report) = candidate.into_parts();
        *self = history;
        *state = native_state;
        Ok(report)
    }
}

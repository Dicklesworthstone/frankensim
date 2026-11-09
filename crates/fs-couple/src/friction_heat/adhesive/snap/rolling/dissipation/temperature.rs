//! Loaded rolling follows original materials at actual receiver temperatures.
//!
//! This attached, same-time endpoint continuation preserves complete native
//! history and all three heat ledgers. The temperature driver accounts for
//! signed normal and rolling coefficient work; no heat is deposited here.

use super::JkrRollingDissipationHeat;
pub use crate::friction_heat::adhesive::snap::angular_temperature::{
    JkrAngularTemperatureHeatControls,
    JkrAngularTemperatureHeatError as JkrRollingTemperatureHeatError,
};
use crate::friction_heat::adhesive::snap::{
    angular_temperature::{check, temperatures},
    rolling::JkrRollingSnapHeatStamp,
};
use fs_contact::normal_patch::adhesive::{
    rolling::temperature::SourceJkrRollingTemperatureReport, temperature::AdhesiveTemperatures,
};
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_mbd::tangential_friction::adhesive::rolling::events::RigidAdhesiveRollingEventState;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
pub struct JkrRollingTemperatureHeatReport {
    /// Complete source identities, normal endpoints and separate signed work.
    pub temperature: SourceJkrRollingTemperatureReport,
    pub temperatures: AdhesiveTemperatures,
    pub parameter_work_source: String,
}

/// Source, native history and heat binding are accepted together. Thermal
/// cells and rigid bodies are read-only throughout this coefficient update.
pub struct JkrRollingTemperatureHeatCandidate<'s> {
    history: JkrRollingDissipationHeat<'s>,
    native_state: RigidAdhesiveRollingEventState,
    report: JkrRollingTemperatureHeatReport,
}
impl<'s> JkrRollingTemperatureHeatCandidate<'s> {
    pub const fn history(&self) -> &JkrRollingDissipationHeat<'s> {
        &self.history
    }
    pub const fn native_state(&self) -> &RigidAdhesiveRollingEventState {
        &self.native_state
    }
    pub const fn report(&self) -> &JkrRollingTemperatureHeatReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        JkrRollingDissipationHeat<'s>,
        RigidAdhesiveRollingEventState,
        JkrRollingTemperatureHeatReport,
    ) {
        (self.history, self.native_state, self.report)
    }
}

impl<'s> JkrRollingDissipationHeat<'s> {
    /// Refresh the original ordered bulk sources at the actual receiving-cell
    /// temperatures and the interface at an explicitly supplied temperature.
    /// Normal indentation, native pair geometry, rolling angle and calibration
    /// stay held. The caller must establish the critical shift's continuing
    /// validity; all original approximation and small-shift gates run again.
    ///
    /// The accepted thermal time must be unchanged. No thermal snap, contact
    /// entropy, time integration or caloric debit is inferred. Signed parameter
    /// work belongs to the named external driver; prior plastic and release
    /// accounts, partitions, clocks and exact caloric bindings survive intact.
    #[allow(clippy::too_many_arguments)]
    pub fn propose_native_temperature_continuation(
        &self,
        state: &RigidAdhesiveRollingEventState,
        thermal: &DepositionSimulation<'_>,
        expected: JkrRollingSnapHeatStamp,
        material_cards: [&MaterialCard; 2],
        interface_temperature_k: f64,
        controls: JkrAngularTemperatureHeatControls<'_>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingTemperatureHeatCandidate<'s>, JkrRollingTemperatureHeatError> {
        check(gate, &mut checkpoint)?;
        controls.validate()?;
        if expected != self.snap.stamp
            || self.state() != state.source_state()
            || thermal.time_s() != self.snap.accepted_time_s
        {
            return Err(JkrRollingTemperatureHeatError::HistoryMismatch);
        }
        self.snap.admit_receivers(thermal)?;
        let revision = self.snap.stamp.revision.checked_add(1).ok_or(
            JkrRollingTemperatureHeatError::Invalid("heat history version exhausted"),
        )?;
        let temperatures = temperatures(
            thermal,
            self.snap.receivers.map(|receiver| receiver.cell),
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
        history.snap.stamp.source = source.normal_source().identity();
        history.snap.source = Arc::new(source);
        history.snap.stamp.state = source_state;
        history.snap.stamp.revision = revision;
        let report = JkrRollingTemperatureHeatReport {
            temperature,
            temperatures,
            parameter_work_source: controls.parameter_work_source.to_owned(),
        };
        check(gate, &mut checkpoint)?;
        Ok(JkrRollingTemperatureHeatCandidate {
            history,
            native_state,
            report,
        })
    }

    /// Publish the original heat binding and native history after the final
    /// caller checkpoint and cancellation check. Refusal changes neither owner.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_native_temperature_continuation(
        &mut self,
        state: &mut RigidAdhesiveRollingEventState,
        thermal: &DepositionSimulation<'_>,
        expected: JkrRollingSnapHeatStamp,
        material_cards: [&MaterialCard; 2],
        interface_temperature_k: f64,
        controls: JkrAngularTemperatureHeatControls<'_>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrRollingTemperatureHeatReport, JkrRollingTemperatureHeatError> {
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

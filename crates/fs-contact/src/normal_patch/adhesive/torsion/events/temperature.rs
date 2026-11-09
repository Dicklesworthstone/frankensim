//! Held-coordinate temperature continuation of a complete attached history.
//!
//! Original normal and torsional owners admit the new material state and
//! supply its signed parameter work, including the refreshed G_III/W screen.
//! Every accumulated release and plastic loss survives unchanged. This endpoint
//! operation supplies no twist, branch event, heat or thermal-path stability.

use super::{JkrTorsionEventState, SourceJkrTorsionError, SourceJkrTorsionEvents};
use crate::normal_patch::adhesive::{
    JkrError, temperature::SourceJkrThermalRefresh,
    torsion::temperature::SourceJkrTorsionTemperatureReport,
};

/// Sealed successor retaining the complete accepted junction history.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrTorsionEventTemperatureStep {
    previous: JkrTorsionEventState,
    state: JkrTorsionEventState,
    source: SourceJkrTorsionEvents,
    report: SourceJkrTorsionTemperatureReport,
}
impl SourceJkrTorsionEventTemperatureStep {
    pub const fn previous_state(&self) -> &JkrTorsionEventState {
        &self.previous
    }
    pub const fn state(&self) -> &JkrTorsionEventState {
        &self.state
    }
    pub const fn source(&self) -> &SourceJkrTorsionEvents {
        &self.source
    }
    pub const fn report(&self) -> &SourceJkrTorsionTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrTorsionEvents,
        JkrTorsionEventState,
        SourceJkrTorsionTemperatureReport,
    ) {
        (self.source, self.state, self.report)
    }
}

impl SourceJkrTorsionEvents {
    /// Transfer an attached, possibly yielded history through an original
    /// normal temperature refresh at held indentation and twist angle.
    ///
    /// The original torsional owner retains its exact elastic and total angular
    /// memory and independently calibrated yield twist. Its material, radius
    /// and G_III/W admission rules still apply. Normal snap release, released
    /// twist storage and plastic torsional loss retain their complete
    /// compensated sums. Signed parameter work stays separate.
    ///
    /// Detached histories refuse: this endpoint transfer cannot choose a
    /// thermal attachment or separation event. Subsequent prescribed mechanical
    /// travel uses the returned source, including stress-free recontact.
    pub fn continue_at_temperature(
        &self,
        previous: &JkrTorsionEventState,
        refresh: SourceJkrThermalRefresh,
    ) -> Result<SourceJkrTorsionEventTemperatureStep, SourceJkrTorsionError> {
        self.admit(previous)?;
        let (contact, memory) = previous.attached().ok_or(JkrError::InvalidInput {
            field: "attached branch for temperature continuation",
        })?;
        let (contact, memory, report) = contact
            .continue_at_temperature(&memory, refresh)?
            .into_parts();
        // A fresh junction uses this admitted contact's unstrained initial
        // state. Do not impose another admission at zero indentation on a
        // valid temperature endpoint at the current prescribed indentation.
        let source = Self {
            seed: contact.clone(),
        };
        // Preserve the original Dd accounts directly, never via f64 getters.
        let mut state = previous.clone();
        state.identity = source.seed.identity();
        state.normal = *contact.normal_response();
        state.torsion_storage_j = report.torsion_storage_after_j;
        state.attached = Some((contact, memory));
        Ok(SourceJkrTorsionEventTemperatureStep {
            previous: previous.clone(),
            state,
            source,
            report,
        })
    }
}

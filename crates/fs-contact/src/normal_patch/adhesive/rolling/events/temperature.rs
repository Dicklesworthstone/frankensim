//! Held-coordinate temperature continuation of a complete attached history.
//!
//! The original normal and rolling owners admit the new material state and
//! supply its signed parameter work. Every accumulated release and plastic
//! loss survives unchanged. This endpoint operation supplies no angular travel,
//! branch event, heat or assertion of stability along a thermal path.

use super::{JkrRollingEventState, SourceJkrRollingError, SourceJkrRollingEvents};
use crate::normal_patch::adhesive::{
    JkrError, rolling::temperature::SourceJkrRollingTemperatureReport,
    temperature::SourceJkrThermalRefresh,
};

/// Sealed successor retaining the complete accepted junction history.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingEventTemperatureStep {
    previous: JkrRollingEventState,
    state: JkrRollingEventState,
    source: SourceJkrRollingEvents,
    report: SourceJkrRollingTemperatureReport,
}
impl SourceJkrRollingEventTemperatureStep {
    pub const fn previous_state(&self) -> &JkrRollingEventState {
        &self.previous
    }
    pub const fn state(&self) -> &JkrRollingEventState {
        &self.state
    }
    pub const fn source(&self) -> &SourceJkrRollingEvents {
        &self.source
    }
    pub const fn report(&self) -> &SourceJkrRollingTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrRollingEvents,
        JkrRollingEventState,
        SourceJkrRollingTemperatureReport,
    ) {
        (self.source, self.state, self.report)
    }
}

impl SourceJkrRollingEvents {
    /// Transfer an attached, possibly yielded history through an original
    /// normal temperature refresh at held indentation and rolling angle.
    ///
    /// The original rolling owner retains its exact elastic and total angular
    /// memory and independently calibrated yield displacement. Its admission
    /// rules still apply to the refreshed material and contact radius. Normal
    /// snap release, released rolling storage and plastic rolling loss retain
    /// their complete compensated sums. Signed parameter work stays separate.
    ///
    /// Detached histories refuse: this endpoint transfer cannot choose a
    /// thermal attachment or separation event. Subsequent prescribed mechanical
    /// travel uses the returned source, including stress-free recontact.
    pub fn continue_at_temperature(
        &self,
        previous: &JkrRollingEventState,
        refresh: SourceJkrThermalRefresh,
    ) -> Result<SourceJkrRollingEventTemperatureStep, SourceJkrRollingError> {
        self.admit(previous)?;
        let (contact, memory) = previous.attached().ok_or(JkrError::InvalidInput {
            field: "attached branch for temperature continuation",
        })?;
        let (contact, memory, report) = contact
            .continue_at_temperature(&memory, refresh)?
            .into_parts();
        // The admitted loaded contact can seed fresh, unstrained junctions.
        // Requiring a second admission at zero indentation here would reject
        // valid endpoint transfers outside that unrelated operating point.
        let source = Self {
            seed: contact.clone(),
        };
        // Preserve the original Dd accounts directly, never via f64 getters.
        let mut state = previous.clone();
        state.identity = source.seed.identity();
        state.normal = *contact.normal_response();
        state.rolling_storage_j = report.rolling_storage_after_j;
        state.attached = Some((contact, memory));
        Ok(SourceJkrRollingEventTemperatureStep {
            previous: previous.clone(),
            state,
            source,
            report,
        })
    }
}

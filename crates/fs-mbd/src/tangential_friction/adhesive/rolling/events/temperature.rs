//! Carry the native pair through an authenticated held-coordinate refresh.

use super::{RigidAdhesiveRollingError, RigidAdhesiveRollingEventState, poll};
use fs_contact::normal_patch::adhesive::rolling::{
    SourceJkrRollingError, events::temperature::SourceJkrRollingEventTemperatureStep,
};
use fs_exec::CancelGate;

impl RigidAdhesiveRollingEventState {
    /// Retain the original pair geometry and complete angular/event history
    /// through a source-owned material-temperature continuation.
    ///
    /// The caller publishes the returned state with `continuation.source()`.
    /// No body impulse, pose motion, time advance or heat is introduced. The
    /// exact predecessor check prevents adopting another loaded pair's history.
    pub fn continue_at_temperature(
        &self,
        continuation: &SourceJkrRollingEventTemperatureStep,
        gate: &CancelGate,
    ) -> Result<Self, RigidAdhesiveRollingError> {
        poll(gate)?;
        if &self.source != continuation.previous_state() {
            return Err(SourceJkrRollingError::StateMismatch.into());
        }
        let next = Self {
            source: continuation.state().clone(),
            geometry: self.geometry,
        };
        poll(gate)?;
        Ok(next)
    }
}

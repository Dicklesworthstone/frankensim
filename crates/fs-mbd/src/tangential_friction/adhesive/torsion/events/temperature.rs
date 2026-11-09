//! Carry the native pair through an authenticated held-coordinate refresh.

use super::{RigidAdhesiveTorsionError, RigidAdhesiveTorsionEventState, poll};
use fs_contact::normal_patch::adhesive::torsion::{
    SourceJkrTorsionError, events::temperature::SourceJkrTorsionEventTemperatureStep,
};
use fs_exec::CancelGate;

impl RigidAdhesiveTorsionEventState {
    /// Retain the original pair geometry and complete angular/event history
    /// through a source-owned material-temperature continuation.
    ///
    /// The caller publishes the returned state with `continuation.source()`.
    /// No body impulse, pose motion, time advance or heat is introduced. The
    /// exact predecessor check prevents adopting another loaded pair's history.
    pub fn continue_at_temperature(
        &self,
        continuation: &SourceJkrTorsionEventTemperatureStep,
        gate: &CancelGate,
    ) -> Result<Self, RigidAdhesiveTorsionError> {
        poll(gate)?;
        if &self.source != continuation.previous_state() {
            return Err(RigidAdhesiveTorsionError::Source(
                SourceJkrTorsionError::StateMismatch,
            ));
        }
        let next = Self {
            source: continuation.state().clone(),
            geometry: self.geometry,
        };
        poll(gate)?;
        Ok(next)
    }
}

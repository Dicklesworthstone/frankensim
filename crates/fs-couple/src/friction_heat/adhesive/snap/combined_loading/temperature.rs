//! Retain combined failure and caloric history across an original source refresh.
//!
//! The sealed source continuation owns the signed change of all four contact
//! stores. That parameter work belongs to its prescribed temperature driver;
//! transferring the history here never deposits it as heat or releases the
//! contact a second time.

use super::{JkrCombinedHeat, JkrCombinedHeatError, JkrCombinedHeatStamp};
use crate::friction_heat::adhesive::snap::JkrSnapHeatError;
use fs_contact::normal_patch::adhesive::mixed_mode::combined::loading::temperature::JkrCombinedLoadingTemperatureContinuation;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

impl<'s> JkrCombinedHeat<'s> {
    /// Stage the exact source successor while retaining every caloric binding
    /// and previously routed release. Publish the returned history together
    /// with the original continuation and its signed parameter-work account.
    ///
    /// The thermal owner may have advanced through its own transport since the
    /// last heat transaction; its time cannot go backwards and the original
    /// cells, masses, geometry and registered caloric laws must still match.
    /// The source operation keeps held approach, shear/twist loading, release
    /// history and step count, and recomputes attached stability itself. A
    /// detached junction remains detached even if no attached branch survives.
    ///
    /// This is a prescribed source update. It does not solve thermal feedback,
    /// equate source temperatures with receiver temperatures, or create caloric
    /// heat from parameter work. No heat cell changes and no conducting link is
    /// introduced. The new revision invalidates the old stamp even for an exact
    /// no-op refresh. Every refusal leaves both accepted owners untouched.
    pub fn propose_temperature_continuation<'new>(
        &self,
        thermal: &DepositionSimulation<'_>,
        expected: JkrCombinedHeatStamp,
        continuation: &JkrCombinedLoadingTemperatureContinuation<'new>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<JkrCombinedHeat<'new>, JkrCombinedHeatError>
    where
        's: 'new,
    {
        checkpoint()?;
        if expected != self.stamp
            || continuation.previous_state() != self.state()
            || continuation.report().original_model_identity != self.loading.source().identity()
            || thermal.time_s() < self.accepted_time_s
        {
            return Err(JkrSnapHeatError::HistoryMismatch.into());
        }
        self.admit_receivers(thermal)?;
        let loading = *continuation.loading();
        let normal = loading.source().normal_source();
        let elastic = normal.normal_state().elastic();
        if self.caloric[0].material_card_identity()
            != elastic.surface_a_state().resolved().card_identity()
            || self.caloric[1].material_card_identity()
                != elastic.surface_b_state().resolved().card_identity()
        {
            return Err(JkrSnapHeatError::Invalid(
                "temperature continuation changed the original ordered material cards",
            )
            .into());
        }
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
        let history = JkrCombinedHeat {
            loading,
            caloric: self.caloric,
            receivers: self.receivers,
            partition: self.partition,
            outside_destination: self.outside_destination.clone(),
            stamp: JkrCombinedHeatStamp {
                source: loading.source().identity(),
                revision,
                state: continuation.state(),
            },
            accepted_time_s: thermal.time_s(),
            ledger: self.ledger,
        };
        checkpoint()?;
        Ok(history)
    }
}

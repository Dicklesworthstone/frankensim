//! Transfer snap history to the original temperature-refreshed normal source.
//!
//! This operation deposits no heat. The source's signed parameter work must be
//! handled by the prescribed temperature driver; existing snap and support
//! losses keep their own explicit thermal destinations.

use super::{JkrSnapHeat, JkrSnapHeatError, JkrSnapHeatStamp};
use fs_contact::normal_patch::adhesive::motion::temperature::JkrNormalTemperatureContinuation;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

impl<'s> JkrSnapHeat<'s> {
    /// Stage a source change without resetting accepted branch or heat history.
    ///
    /// Requires the original sealed normal-motion continuation, the current
    /// heat stamp and the unchanged accepted thermal time. The returned history
    /// retains the exact receivers, caloric objects, partitions and cumulative
    /// release account. Its new version invalidates the old heat stamp even for
    /// an otherwise identical source refresh. No caloric cell is cloned or heated.
    /// Publish this candidate together with its mechanical continuation.
    pub fn propose_temperature_continuation(
        &self,
        thermal: &DepositionSimulation<'_>,
        expected: JkrSnapHeatStamp,
        continuation: &JkrNormalTemperatureContinuation<'s>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, JkrSnapHeatError> {
        checkpoint()?;
        let report = continuation.report();
        if expected != self.stamp
            || report.original_model_identity != self.source.identity()
            || continuation.previous_state().contact_state() != self.state()
            || self
                .normal_motion_state
                .is_some_and(|state| state != continuation.previous_state())
            || thermal.time_s() != self.accepted_time_s
        {
            return Err(JkrSnapHeatError::HistoryMismatch);
        }
        let source = continuation.motion().normal();
        let elastic = source.normal_state().elastic();
        if self.caloric[0].material_card_identity()
            != elastic.surface_a_state().resolved().card_identity()
            || self.caloric[1].material_card_identity()
                != elastic.surface_b_state().resolved().card_identity()
        {
            return Err(JkrSnapHeatError::Invalid(
                "temperature continuation changed the original ordered material cards",
            ));
        }
        self.admit_receivers(thermal)?;
        let revision = self
            .stamp
            .revision
            .checked_add(1)
            .ok_or(JkrSnapHeatError::Invalid("history version exhausted"))?;
        let history = Self {
            source,
            caloric: self.caloric,
            receivers: self.receivers,
            partition: self.partition,
            outside_destination: self.outside_destination.clone(),
            stamp: JkrSnapHeatStamp {
                source: source.identity(),
                revision,
                state: continuation.state().contact_state(),
            },
            normal_motion_state: Some(continuation.state()),
            accepted_time_s: self.accepted_time_s,
            ledger: self.ledger,
        };
        checkpoint()?;
        Ok(history)
    }
}

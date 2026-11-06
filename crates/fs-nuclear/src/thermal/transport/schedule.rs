//! Authored first-order nuclide movement on the original absolute thermal clock.

use std::sync::Arc;

use crate::FissionError;
pub use crate::burnup::transport::schedule::NuclideTransportChange;
use crate::burnup::transport::schedule::{NuclideTransportTimeline, NuclideTransportTimelineError};
use crate::burnup::transport::{BurnupTransportLimits, BurnupTransportState, NuclideTransport};
use crate::burnup::{BurnupError, BurnupModel};

/// Bounded movement history bound to the same nuclides, volumes and path order.
///
/// Rates are first-order hazards [1/s]. A well-mixed volumetric flow Q requires
/// an explicitly authored Q/V_donor; no velocity or bulk enthalpy is inferred.
/// Flux and external nuclide production retain their existing separate schedules.
#[derive(Clone, Debug)]
pub struct NuclideTransportSchedule {
    model: Arc<BurnupModel>,
    volumes_m3: Vec<f64>,
    timeline: NuclideTransportTimeline,
}

impl NuclideTransportSchedule {
    pub fn new(
        state: &BurnupTransportState,
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        limits: BurnupTransportLimits,
        maximum_changes: usize,
    ) -> Result<Self, BurnupError> {
        Self::new_with_checkpoint(state, initial, changes, limits, maximum_changes, || Ok(()))
    }

    /// Admit the original physical shape and every retained rate vector.
    /// Cancellation leaves both the state and schedules unchanged.
    pub fn new_with_checkpoint(
        state: &BurnupTransportState,
        initial: NuclideTransport,
        changes: Vec<NuclideTransportChange>,
        limits: BurnupTransportLimits,
        maximum_changes: usize,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, BurnupError> {
        checkpoint()?;
        if changes.len() > maximum_changes {
            return Err(BurnupError::Budget("nuclide movement schedule changes"));
        }
        let zero_sources =
            vec![vec![0.0; state.model().nuclides().len()]; state.inventories().len()];
        state.validate_inputs(&initial, &zero_sources, limits, &mut checkpoint)?;
        // The initial physical admission binds every dimension and path. The
        // shared timeline then retains those shapes and validates each future
        // coefficient once, using this owner's original callback/error type.
        let timeline = NuclideTransportTimeline::new_with_checkpoint(
            initial,
            changes,
            maximum_changes,
            || checkpoint().map_err(BurnupError::from),
        )?;
        Ok(Self {
            model: Arc::clone(state.model()),
            volumes_m3: state.volumes_m3().to_vec(),
            timeline,
        })
    }

    pub(super) fn constant(state: &BurnupTransportState) -> Self {
        Self {
            model: Arc::clone(state.model()),
            volumes_m3: state.volumes_m3().to_vec(),
            timeline: NuclideTransportTimeline::admitted_constant(state.transport().clone()),
        }
    }

    /// Active movement rates at finite nonnegative absolute time [s].
    pub fn transport_at(&self, time_s: f64) -> Result<&NuclideTransport, BurnupError> {
        self.timeline.transport_at(time_s).map_err(Into::into)
    }

    #[must_use]
    pub fn changes(&self) -> &[NuclideTransportChange] {
        self.timeline.changes()
    }

    /// First authored discontinuity strictly after this absolute time.
    #[must_use]
    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        self.timeline.next_after(time_s)
    }

    pub(super) fn validate_binding(&self, state: &BurnupTransportState) -> Result<(), BurnupError> {
        if !Arc::ptr_eq(&self.model, state.model())
            || self.volumes_m3 != state.volumes_m3()
            || !self.timeline.initial().same_topology(state.transport())
        {
            return Err(BurnupError::Invalid(
                "movement schedule belongs to different nuclides, volumes or paths",
            ));
        }
        Ok(())
    }
}

impl From<NuclideTransportTimelineError> for BurnupError {
    fn from(error: NuclideTransportTimelineError) -> Self {
        match error {
            NuclideTransportTimelineError::Invalid(reason) => Self::Invalid(reason),
            NuclideTransportTimelineError::Budget(reason) => Self::Budget(reason),
        }
    }
}

//! Authored first-order nuclide movement on the original absolute thermal clock.

use std::sync::Arc;

use crate::burnup::transport::{BurnupTransportLimits, BurnupTransportState, NuclideTransport};
use crate::burnup::{BurnupError, BurnupModel};
use crate::{FissionError, nonnegative};

/// Right-continuous replacement of all movement rates in the original path order.
#[derive(Clone, Debug, PartialEq)]
pub struct NuclideTransportChange {
    pub time_s: f64,
    pub transport: NuclideTransport,
}

/// Bounded movement history bound to the same nuclides, volumes and path order.
///
/// Rates are first-order hazards [1/s]. A well-mixed volumetric flow Q requires
/// an explicitly authored Q/V_donor; no velocity or bulk enthalpy is inferred.
/// Flux and external nuclide production retain their existing separate schedules.
#[derive(Clone, Debug)]
pub struct NuclideTransportSchedule {
    model: Arc<BurnupModel>,
    volumes_m3: Vec<f64>,
    initial: NuclideTransport,
    changes: Vec<NuclideTransportChange>,
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

    /// Admit every rate vector through the original transport owner before
    /// retaining it. Cancellation leaves both the state and schedules unchanged.
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
        let mut previous = 0.0;
        for change in &changes {
            checkpoint()?;
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(BurnupError::Invalid(
                    "movement changes need strictly increasing positive times",
                ));
            }
            state.validate_inputs(&change.transport, &zero_sources, limits, &mut checkpoint)?;
            previous = change.time_s;
        }
        Ok(Self {
            model: Arc::clone(state.model()),
            volumes_m3: state.volumes_m3().to_vec(),
            initial,
            changes,
        })
    }

    pub(super) fn constant(state: &BurnupTransportState) -> Self {
        Self {
            model: Arc::clone(state.model()),
            volumes_m3: state.volumes_m3().to_vec(),
            initial: state.transport().clone(),
            changes: Vec::new(),
        }
    }

    /// Active movement rates at finite nonnegative absolute time [s].
    pub fn transport_at(&self, time_s: f64) -> Result<&NuclideTransport, BurnupError> {
        if !nonnegative(time_s) {
            return Err(BurnupError::Invalid(
                "movement schedule time must be finite and nonnegative",
            ));
        }
        let after = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if after == 0 {
            &self.initial
        } else {
            &self.changes[after - 1].transport
        })
    }

    #[must_use]
    pub fn changes(&self) -> &[NuclideTransportChange] {
        &self.changes
    }

    /// First authored discontinuity strictly after this absolute time.
    #[must_use]
    pub fn next_after(&self, time_s: f64) -> Option<f64> {
        if !nonnegative(time_s) {
            return None;
        }
        self.changes
            .get(
                self.changes
                    .partition_point(|change| change.time_s <= time_s),
            )
            .map(|change| change.time_s)
    }

    pub(super) fn validate_binding(&self, state: &BurnupTransportState) -> Result<(), BurnupError> {
        if !Arc::ptr_eq(&self.model, state.model())
            || self.volumes_m3 != state.volumes_m3()
            || !self.initial.same_topology(state.transport())
        {
            return Err(BurnupError::Invalid(
                "movement schedule belongs to different nuclides, volumes or paths",
            ));
        }
        Ok(())
    }
}

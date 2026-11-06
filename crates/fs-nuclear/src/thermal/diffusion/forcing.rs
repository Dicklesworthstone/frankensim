//! Explicit external particle and nuclide histories for a spatial cohort.
//!
//! Values are extensive expected particles/s in the original spatial order.
//! They describe incoming nuclear populations only: neither their kinetic
//! energy nor carrier mass, fluid flow or incoming enthalpy is inferred.

use std::sync::Arc;

use crate::burnup::BurnupModel;
use crate::diffusion::DiffusionModel;
use crate::diffusion::burnup::{
    DiffusionBurnupError, DiffusionBurnupSources, DiffusionBurnupState,
};
use crate::{FissionError, nonnegative};

/// Right-continuous replacement of all original-order external source rates.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionSourceChange {
    pub time_s: f64,
    pub sources: DiffusionBurnupSources,
}

/// An immutable bounded history shared by every cell in one spatial source.
/// The initial rates apply from time zero. Binding later uses the explicitly
/// supplied current inventories; it never replays an unobserved prior history.
#[derive(Clone, Debug)]
pub struct DiffusionSourceSchedule {
    neutron_model: Arc<DiffusionModel>,
    burnup_models: Vec<Arc<BurnupModel>>,
    initial: DiffusionBurnupSources,
    changes: Vec<DiffusionSourceChange>,
}

impl DiffusionSourceSchedule {
    pub fn new(
        state: &DiffusionBurnupState,
        initial: DiffusionBurnupSources,
        changes: Vec<DiffusionSourceChange>,
        maximum_changes: usize,
    ) -> Result<Self, DiffusionBurnupError> {
        Self::new_with_checkpoint(state, initial, changes, maximum_changes, || Ok(()))
    }

    /// Validate every vector before retaining it, with bounded checkpoints.
    pub fn new_with_checkpoint(
        state: &DiffusionBurnupState,
        initial: DiffusionBurnupSources,
        changes: Vec<DiffusionSourceChange>,
        maximum_changes: usize,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DiffusionBurnupError> {
        checkpoint()?;
        if changes.len() > maximum_changes {
            return Err(FissionError::Budget("spatial source schedule changes").into());
        }
        initial.validate_with_checkpoint(state, &mut checkpoint)?;
        let mut previous = 0.0;
        for change in &changes {
            checkpoint()?;
            if !change.time_s.is_finite() || change.time_s <= previous {
                return Err(DiffusionBurnupError::Invalid(
                    "spatial source changes need strictly increasing positive times",
                ));
            }
            change
                .sources
                .validate_with_checkpoint(state, &mut checkpoint)?;
            previous = change.time_s;
        }
        let mut burnup_models = Vec::with_capacity(state.inventories().len());
        for inventory in state.inventories() {
            checkpoint()?;
            burnup_models.push(Arc::clone(inventory.model()));
        }
        Ok(Self {
            neutron_model: Arc::clone(state.neutrons().model()),
            burnup_models,
            initial,
            changes,
        })
    }

    /// Active rates at a finite nonnegative absolute time; lookup is logarithmic.
    pub fn sources_at(&self, time_s: f64) -> Result<&DiffusionBurnupSources, DiffusionBurnupError> {
        if !nonnegative(time_s) {
            return Err(DiffusionBurnupError::Invalid(
                "spatial source schedule time must be finite and nonnegative",
            ));
        }
        let after = self
            .changes
            .partition_point(|change| change.time_s <= time_s);
        Ok(if after == 0 {
            &self.initial
        } else {
            &self.changes[after - 1].sources
        })
    }

    #[must_use]
    pub fn changes(&self) -> &[DiffusionSourceChange] {
        &self.changes
    }

    /// First discontinuity strictly after the admitted time.
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

    pub(super) fn validate_binding(
        &self,
        state: &DiffusionBurnupState,
    ) -> Result<(), DiffusionBurnupError> {
        if self.burnup_models.len() != state.inventories().len()
            || !self
                .burnup_models
                .iter()
                .zip(state.inventories())
                .all(|(model, inventory)| Arc::ptr_eq(model, inventory.model()))
        {
            return Err(DiffusionBurnupError::Invalid(
                "spatial source schedule belongs to different ordered burnup models",
            ));
        }
        // Reuse the original owner's geometry, neutron-grid and precursor
        // identity guard. Frozen coefficients may have changed after burnup.
        state
            .neutrons()
            .rebind_model(Arc::clone(&self.neutron_model))?;
        Ok(())
    }
}

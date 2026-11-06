//! One spatial neutron/burnup cohort on the original thermal transaction.
//!
//! All cell handles share an immutable complete physical checkpoint. The batch
//! source hook advances that checkpoint exactly once per native thermal step,
//! then the original NuclearHeating owner credits each cell's deposited heat.

use std::{collections::BTreeSet, sync::Arc};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating};
use crate::burnup::BurnupStep;
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState, DiffusionBurnupStep};
use crate::{EnergyChannels, FissionError};
use fs_thermal::enthalpy::{CellId, MaterialId};

/// A spatial cell handle. Use `bind` to create the complete shared collection.
#[derive(Clone, Debug)]
pub struct DiffusionCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    cell_index: usize,
    state: Arc<DiffusionBurnupState>,
}

/// Spatial neutron diffusion, finite burnup and heat on the native clock.
pub type DiffusionHeating = NuclearHeating<DiffusionCellSource>;

#[derive(Clone, Debug)]
pub struct DiffusionCellStep {
    cell_index: usize,
    spatial: Arc<DiffusionBurnupStep>,
}
impl DiffusionCellStep {
    #[must_use]
    pub fn burnup(&self) -> &BurnupStep {
        &self.spatial.burnup[self.cell_index]
    }
    #[must_use]
    pub fn spatial(&self) -> &DiffusionBurnupStep {
        &self.spatial
    }
    #[must_use]
    pub fn cell_index(&self) -> usize {
        self.cell_index
    }
}

impl DiffusionCellSource {
    /// Bind diffusion-cell order to distinct active thermal cells/materials.
    /// The original thermal owner also checks each reference geometry volume.
    pub fn bind(
        state: Arc<DiffusionBurnupState>,
        cells: Vec<(CellId, MaterialId)>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        if cells.len() != state.inventories().len() {
            return Err(FissionThermalError::Invalid(
                "complete spatial-to-thermal cell correspondence required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut sources = Vec::with_capacity(cells.len());
        for (cell_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !used.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct spatial thermal cells required",
                ));
            }
            sources.push(Self {
                cell_id,
                material_id,
                cell_index,
                state: Arc::clone(&state),
            });
        }
        Ok(sources)
    }
    #[must_use]
    pub fn state(&self) -> &Arc<DiffusionBurnupState> {
        &self.state
    }
    #[must_use]
    pub fn cell_index(&self) -> usize {
        self.cell_index
    }
}

impl sealed::Sealed for DiffusionCellSource {}
impl NuclearHeatSource for DiffusionCellSource {
    type Budget = DiffusionBurnupBudget;
    type Step = DiffusionCellStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn expected_volume_m3(&self) -> Option<f64> {
        Some(self.state.neutrons().model().cell_volumes_m3()[self.cell_index])
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.burnup.energy_tolerance_j
    }
    fn nuclei_tolerance(budget: Self::Budget) -> Option<f64> {
        Some(budget.burnup.nuclei_tolerance + budget.reaction_tolerance)
    }
    fn next_transition(&self, _: f64) -> Option<f64> {
        None
    }
    fn admit(
        &self,
        _: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        budget.validate()?;
        if self.cell_index == 0 {
            self.state
                .neutrons()
                .trial(0.0, budget.diffusion, || Ok(()))
                .map_err(crate::diffusion::burnup::DiffusionBurnupError::from)?;
        }
        self.state.inventories()[self.cell_index].trial(
            0.0,
            &self
                .state
                .neutrons()
                .group_flux(self.cell_index)
                .map_err(crate::diffusion::burnup::DiffusionBurnupError::from)?,
            temperature_k,
            budget.burnup,
            || Ok(()),
        )?;
        Ok(())
    }
    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.state.inventories()[self.cell_index]
            .model()
            .validate_temperature(temperature_k)?;
        Ok(())
    }
    fn trial_source(
        &self,
        _: f64,
        _: f64,
        _: f64,
        _: Self::Budget,
        _: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        Err(FissionThermalError::Invalid(
            "spatial neutron evolution requires the complete cell collection",
        ))
    }
    fn trial_sources(
        sources: &[&Self],
        _: f64,
        duration_s: f64,
        temperatures_k: &[f64],
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(Self, Self::Step)>, FissionThermalError> {
        checkpoint()?;
        let first = sources.first().ok_or(FissionThermalError::Invalid(
            "spatial collection must be nonempty",
        ))?;
        if sources.len() != first.state.inventories().len() || temperatures_k.len() != sources.len()
        {
            return Err(FissionThermalError::Invalid(
                "all spatial cells must advance together",
            ));
        }
        let mut seen = vec![false; sources.len()];
        let mut temperatures = vec![0.0; sources.len()];
        for (source, &temperature) in sources.iter().zip(temperatures_k) {
            checkpoint()?;
            if !Arc::ptr_eq(&first.state, &source.state)
                || source.cell_index >= seen.len()
                || seen[source.cell_index]
            {
                return Err(FissionThermalError::Invalid(
                    "spatial cells must share one complete accepted history",
                ));
            }
            seen[source.cell_index] = true;
            temperatures[source.cell_index] = temperature;
        }
        let (state, report) =
            first
                .state
                .trial(duration_s, &temperatures, budget, &mut *checkpoint)?;
        let state = Arc::new(state);
        let report = Arc::new(report);
        let mut candidates = Vec::with_capacity(sources.len());
        for source in sources {
            checkpoint()?;
            let mut candidate = (*source).clone();
            candidate.state = Arc::clone(&state);
            candidates.push((
                candidate,
                DiffusionCellStep {
                    cell_index: source.cell_index,
                    spatial: Arc::clone(&report),
                },
            ));
        }
        Ok(candidates)
    }
    fn released(step: &Self::Step) -> EnergyChannels {
        step.burnup().released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.burnup().energy_roundoff_j
    }
    fn population_roundoff_nuclei(step: &Self::Step) -> f64 {
        step.burnup().population_roundoff_nuclei
            + step.spatial.reaction_defect_by_cell[step.cell_index]
    }
    fn decays(step: &Self::Step) -> f64 {
        step.burnup().decays
    }
    fn fissions(step: &Self::Step) -> f64 {
        step.burnup().fissions
    }
    fn captures(step: &Self::Step) -> f64 {
        step.burnup().captures
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.inventories()[self.cell_index].cumulative_released())
    }
}

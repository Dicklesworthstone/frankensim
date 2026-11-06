//! Scheduled Poisson supply and radioactive heat on the original material clock.
//!
//! Births, source changes and decays clip the shared thermal timeline. Births
//! carry no inferred heat: only a subsequent original decay branch releases its
//! authored energy. The complete production history remains inside the existing
//! nuclear/material transaction, including pending clocks and random-stream state.

use std::collections::BTreeSet;

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating};
use crate::decay::DecayError;
use crate::decay::stochastic::production::{
    StochasticDecayProduction, StochasticDecayProductionBudget, StochasticDecayProductionStep,
};
use crate::{EnergyChannels, FissionError};

/// An explicitly Poisson-produced integer inventory in one original thermal cell.
///
/// The production owner's absolute clock must match the thermal clock. Distinct
/// cells require independent logical streams. Bulk carrier mass and enthalpy,
/// activation and creation energy remain with their separately authored owners.
#[derive(Clone, Debug)]
pub struct StochasticDecayProductionCellSource {
    pub cell_id: CellId,
    pub material_id: MaterialId,
    pub production: StochasticDecayProduction,
}

/// The existing atomic thermal transaction with sampled supply and decay events.
pub type StochasticDecayProductionHeating = NuclearHeating<StochasticDecayProductionCellSource>;

impl StochasticDecayProductionCellSource {
    fn check_clock(&self, thermal_time_s: f64) -> Result<(), FissionThermalError> {
        let source_time = self.production.inventory().time_s();
        if source_time != thermal_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s,
                fission_time_s: source_time,
            });
        }
        Ok(())
    }

    fn trial_to(
        &self,
        start_s: f64,
        end_s: f64,
        budget: StochasticDecayProductionBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, StochasticDecayProductionStep), FissionThermalError> {
        self.check_clock(start_s)?;
        let mut checkpoint_error = None;
        let result = self.production.trial_to(end_s, budget, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                DecayError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (production, step) = result?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                production,
            },
            step,
        ))
    }
}

impl sealed::Sealed for StochasticDecayProductionCellSource {}

impl NuclearHeatSource for StochasticDecayProductionCellSource {
    type Budget = StochasticDecayProductionBudget;
    type Step = StochasticDecayProductionStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.decay.energy_tolerance_j
    }

    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.production
            .next_transition_time_s()
            .filter(|&time| time > time_s)
    }

    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.check_clock(time_s)?;
        self.production.trial_to(time_s, budget, || Ok(()))?;
        Ok(())
    }

    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
        // Constant authored nuclear laws do not invent a temperature correction.
        Ok(())
    }

    fn validate_sources<'source>(
        sources: impl ExactSizeIterator<Item = &'source Self>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionThermalError>,
    ) -> Result<(), FissionThermalError>
    where
        Self: 'source,
    {
        let mut keys = BTreeSet::new();
        for source in sources {
            checkpoint()?;
            let key = source.production.inventory().rng_checkpoint().key;
            if !keys.insert((key.seed, key.kernel, key.tile)) {
                return Err(FissionThermalError::Invalid(
                    "independent stochastic production cells require distinct stream keys",
                ));
            }
        }
        Ok(())
    }

    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        _: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(FissionThermalError::Invalid(
                "invalid stochastic production interval",
            ));
        }
        self.trial_to(time_s, time_s + duration_s, budget, checkpoint)
    }

    fn trial_sources_to(
        sources: &[&Self],
        time_s: f64,
        _: f64,
        end_s: f64,
        temperatures_k: &[f64],
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Vec<(Self, Self::Step)>, FissionThermalError> {
        if sources.len() != temperatures_k.len() {
            return Err(FissionThermalError::Invalid(
                "one temperature per production source required",
            ));
        }
        sources
            .iter()
            .map(|source| {
                checkpoint()?;
                source.trial_to(time_s, end_s, budget, checkpoint)
            })
            .collect()
    }

    fn released(step: &Self::Step) -> EnergyChannels {
        step.decay.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.decay.energy_roundoff_j
    }
    fn decays(step: &Self::Step) -> f64 {
        step.decay.events as f64
    }
    fn produced_nuclei(step: &Self::Step) -> f64 {
        step.births as f64
    }
    fn fissions(_: &Self::Step) -> f64 {
        0.0
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.production.inventory().cumulative_released())
    }
}

#[cfg(test)]
mod tests;

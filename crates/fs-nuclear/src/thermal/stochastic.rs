//! Individual radioactive transformations on the original thermal timeline.
//!
//! The next sampled decay is an explicit thermal step boundary. Heat is added
//! only at the endpoint containing that event, while the pre-event material
//! state supplies the intervening heat fluxes. Daughters retain their actual
//! sampled birth times and original half-life laws. This is a finite-population
//! source for the existing uniform-cell thermal model, not particle or photon
//! transport. Event energies and their deposition fractions are supplied data.

use std::collections::BTreeSet;

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating};
use crate::decay::DecayError;
use crate::decay::stochastic::{
    StochasticDecayBudget, StochasticDecayInventory, StochasticDecayStep,
};
use crate::{EnergyChannels, FissionError};

/// An integer radioactive inventory bound to one original thermal cell.
///
/// The inventory's absolute clock must equal the thermal clock at binding.
/// Independent cells require distinct logical random-stream keys. Cloning an
/// entire heating history deliberately retains its streams for replay.
#[derive(Clone, Debug)]
pub struct StochasticDecayCellSource {
    pub cell_id: CellId,
    pub material_id: MaterialId,
    pub inventory: StochasticDecayInventory,
}

/// Event-driven decay heat through the existing nuclear/thermal transaction.
pub type StochasticDecayHeating = NuclearHeating<StochasticDecayCellSource>;

impl StochasticDecayCellSource {
    fn trial_to(
        &self,
        start_s: f64,
        end_s: f64,
        budget: StochasticDecayBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, StochasticDecayStep), FissionThermalError> {
        self.check_clock(start_s)?;
        let mut checkpoint_error = None;
        let trial = self.inventory.trial_to(end_s, budget, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                DecayError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (inventory, step) = trial?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                inventory,
            },
            step,
        ))
    }

    fn check_clock(&self, thermal_time_s: f64) -> Result<(), FissionThermalError> {
        if self.inventory.time_s() != thermal_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s,
                fission_time_s: self.inventory.time_s(),
            });
        }
        Ok(())
    }
}

impl sealed::Sealed for StochasticDecayCellSource {}
impl NuclearHeatSource for StochasticDecayCellSource {
    type Budget = StochasticDecayBudget;
    type Step = StochasticDecayStep;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.energy_tolerance_j
    }

    fn next_transition(&self, time_s: f64) -> Option<f64> {
        self.inventory
            .next_event_time_s()
            .filter(|&next| next > time_s)
    }

    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.check_clock(time_s)?;
        self.inventory.trial_to(time_s, budget, || Ok(()))?;
        Ok(())
    }

    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
        // These are the original constant nuclear laws; no temperature law
        // or environmental correction to a measured half-life is invented.
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
            let key = source.inventory.rng_checkpoint().key;
            if !keys.insert((key.seed, key.kernel, key.tile)) {
                return Err(FissionThermalError::Invalid(
                    "independent stochastic decay cells require distinct stream keys",
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
                "invalid stochastic decay interval",
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
                "one temperature per stochastic source required",
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
        step.released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.energy_roundoff_j
    }
    fn decays(step: &Self::Step) -> f64 {
        // The stochastic owner admits at most 2^53 events per call, keeping
        // this shared audit conversion exact. Integer branch accounts remain
        // available on the original report and inventory.
        step.events as f64
    }
    fn fissions(_: &Self::Step) -> f64 {
        0.0
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.inventory.cumulative_released())
    }
}

#[cfg(test)]
mod tests;

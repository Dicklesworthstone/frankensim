//! Integer capture and radioactive decay on the original material transaction.
//!
//! Capture consumes a represented target, forms its explicit daughter at the
//! sampled event time and releases only that channel's prompt energy. Daughter
//! half-lives continue after irradiation shuts down. Decay and capture retain
//! separate counts and energy; material heat receives their deposited sum once.
//!
//! The nuclear flux and temperature are the explicitly authored history in the
//! prepared activation schedule. They are fixed between its changes. Material
//! temperature does not silently redraw future nuclear clocks or supply an
//! implicit temperature-feedback law. Coupled temperature-dependent irradiation
//! requires an appropriately authored/refined nuclear history.

use std::{collections::BTreeSet, sync::Arc};

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating};
use crate::burnup::stochastic::{
    StochasticActivationBudget, StochasticActivationState, StochasticActivationStep,
};
use crate::endf::EndfError;
use crate::radiation::decay::{BurnupDecayRadiation, DecayRadiationStep, PreparedDecayRadiation};
use crate::{EnergyChannels, FissionError};

/// One finite target/product population with an explicit nuclear forcing history.
/// Independent cells need distinct stream keys; cloning the entire coupled
/// history deliberately preserves its randomness for checkpoint/restart.
#[derive(Clone, Debug)]
pub struct StochasticActivationCellSource {
    pub cell_id: CellId,
    pub material_id: MaterialId,
    pub state: StochasticActivationState,
}

/// Event-resolved activation and daughter heat through the shared thermal owner.
pub type StochasticActivationHeating = NuclearHeating<StochasticActivationCellSource>;

impl StochasticActivationCellSource {
    fn check_clock(&self, thermal_time_s: f64) -> Result<(), FissionThermalError> {
        if self.state.time_s() != thermal_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s,
                fission_time_s: self.state.time_s(),
            });
        }
        Ok(())
    }

    fn trial_to(
        &self,
        start_s: f64,
        end_s: f64,
        budget: StochasticActivationBudget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, StochasticActivationStep), FissionThermalError> {
        self.check_clock(start_s)?;
        let (state, step) = self.state.trial_to(end_s, budget, checkpoint)?;
        Ok((
            Self {
                cell_id: self.cell_id,
                material_id: self.material_id,
                state,
            },
            step,
        ))
    }
}

impl sealed::Sealed for StochasticActivationCellSource {}

impl NuclearHeatSource for StochasticActivationCellSource {
    type Budget = StochasticActivationBudget;
    type Step = StochasticActivationStep;

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
        self.state
            .next_transition_time_s()
            .filter(|&next| next > time_s)
    }

    fn admit(
        &self,
        time_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
    ) -> Result<(), FissionThermalError> {
        self.check_clock(time_s)?;
        self.validate_temperature(temperature_k)?;
        self.state.trial_to(time_s, budget, || Ok(()))?;
        Ok(())
    }

    fn validate_temperature(&self, temperature_k: f64) -> Result<(), FissionThermalError> {
        self.state.model().validate_temperature(temperature_k)?;
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
            let key = source.state.rng_checkpoint().key;
            if !keys.insert((key.seed, key.kernel, key.tile)) {
                return Err(FissionThermalError::Invalid(
                    "independent stochastic activation cells require distinct stream keys",
                ));
            }
        }
        Ok(())
    }

    fn trial_source(
        &self,
        time_s: f64,
        duration_s: f64,
        temperature_k: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, Self::Step), FissionThermalError> {
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(FissionThermalError::Invalid(
                "invalid stochastic activation interval",
            ));
        }
        self.validate_temperature(temperature_k)?;
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
                "one temperature per stochastic activation source required",
            ));
        }
        sources
            .iter()
            .zip(temperatures_k)
            .map(|(source, &temperature_k)| {
                checkpoint()?;
                source.validate_temperature(temperature_k)?;
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
        step.decays as f64
    }

    fn captures(step: &Self::Step) -> f64 {
        step.captures as f64
    }

    fn fissions(_: &Self::Step) -> f64 {
        0.0
    }

    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.cumulative_released())
    }
}

impl super::spectral::DecayRadiationHeatSource for StochasticActivationCellSource {
    type RadiationBinding = BurnupDecayRadiation;

    fn bind_decay_radiation(
        &self,
        partition: Arc<PreparedDecayRadiation>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<Self::RadiationBinding, EndfError> {
        partition.bind_burnup(Arc::clone(self.state.model()), checkpoint)
    }

    fn decay_radiation_step(
        &self,
        binding: &Self::RadiationBinding,
        step: &Self::Step,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionError>,
    ) -> Result<DecayRadiationStep, EndfError> {
        binding.from_stochastic_activation_step(self.state.model(), step, checkpoint)
    }
}

#[cfg(test)]
mod tests;

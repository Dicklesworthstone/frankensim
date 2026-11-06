//! Integer radioactive movement and event-local heat in a complete cell cohort.
//!
//! Every cell handle shares one accepted physical history. The grouped source
//! trial advances the cohort once, so a transferred nucleus can decay in its
//! recipient during the same interval. Earlier emissions stay at their original
//! cell. Movement itself adds no nuclear heat or inferred carrier enthalpy.

use std::{collections::BTreeSet, sync::Arc};

use fs_thermal::enthalpy::{CellId, MaterialId};

use super::source::{NuclearHeatSource, sealed};
use super::{FissionThermalError, NuclearHeating};
use crate::decay::DecayError;
use crate::decay::stochastic::StochasticDecayStep;
use crate::decay::stochastic::transport::{
    StochasticDecayTransportBudget, StochasticDecayTransportState, StochasticDecayTransportStep,
};
use crate::{EnergyChannels, FissionError};

/// One thermal cell's immutable view of the complete stochastic cohort.
/// Construct all handles together with [`Self::bind`], then advance them together.
#[derive(Clone, Debug)]
pub struct StochasticDecayTransportCellSource {
    cell_id: CellId,
    material_id: MaterialId,
    compartment_index: usize,
    state: Arc<StochasticDecayTransportState>,
}

/// Sampled movement, radioactive transformation and original thermal storage.
pub type StochasticDecayTransportHeating = NuclearHeating<StochasticDecayTransportCellSource>;

/// The original compartment emission account plus the shared movement report.
#[derive(Clone, Debug)]
pub struct StochasticDecayTransportCellStep {
    compartment_index: usize,
    transport: Arc<StochasticDecayTransportStep>,
}

impl StochasticDecayTransportCellStep {
    pub fn decay(&self) -> &StochasticDecayStep {
        &self.transport.decay[self.compartment_index]
    }
    pub fn transport(&self) -> &StochasticDecayTransportStep {
        &self.transport
    }
    pub fn compartment_index(&self) -> usize {
        self.compartment_index
    }
}

impl StochasticDecayTransportCellSource {
    /// Bind one distinct original thermal cell to every compartment. Rates are
    /// already authored per-nucleus hazards; binding does not infer carrier flow.
    pub fn bind(
        state: Arc<StochasticDecayTransportState>,
        cells: Vec<(CellId, MaterialId)>,
    ) -> Result<Vec<Self>, FissionThermalError> {
        if cells.len() != state.inventories().len() {
            return Err(FissionThermalError::Invalid(
                "one thermal cell per stochastic compartment required",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut sources = Vec::with_capacity(cells.len());
        for (compartment_index, (cell_id, material_id)) in cells.into_iter().enumerate() {
            if !seen.insert(cell_id) {
                return Err(FissionThermalError::Invalid(
                    "distinct stochastic movement thermal cells required",
                ));
            }
            sources.push(Self {
                cell_id,
                material_id,
                compartment_index,
                state: Arc::clone(&state),
            });
        }
        Ok(sources)
    }

    pub fn state(&self) -> &Arc<StochasticDecayTransportState> {
        &self.state
    }
    pub fn compartment_index(&self) -> usize {
        self.compartment_index
    }

    fn check_clock(&self, thermal_time_s: f64) -> Result<(), FissionThermalError> {
        if self.state.time_s() != thermal_time_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s,
                fission_time_s: self.state.time_s(),
            });
        }
        Ok(())
    }
}

impl sealed::Sealed for StochasticDecayTransportCellSource {}

impl NuclearHeatSource for StochasticDecayTransportCellSource {
    type Budget = StochasticDecayTransportBudget;
    type Step = StochasticDecayTransportCellStep;

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
        // Complete-cohort validation guarantees this representative exists.
        // The thermal owner takes the union of source boundaries: advertising
        // the same shared event from every cell would rescan all clocks C times.
        if self.compartment_index == 0 {
            self.state.next_event_time_s().filter(|&next| next > time_s)
        } else {
            None
        }
    }

    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.check_clock(time_s)?;
        if self.compartment_index == 0 {
            self.state.trial_to(time_s, budget, || Ok(()))?;
        }
        Ok(())
    }
    fn validate_temperature(&self, _: f64) -> Result<(), FissionThermalError> {
        Ok(())
    }

    fn validate_sources<'source>(
        mut sources: impl ExactSizeIterator<Item = &'source Self>,
        checkpoint: &mut dyn FnMut() -> Result<(), FissionThermalError>,
    ) -> Result<(), FissionThermalError>
    where
        Self: 'source,
    {
        checkpoint()?;
        let count = sources.len();
        let first = sources.next().ok_or(FissionThermalError::Invalid(
            "stochastic movement collection must be nonempty",
        ))?;
        if count != first.state.inventories().len() {
            return Err(FissionThermalError::Invalid(
                "all stochastic movement compartments must advance together",
            ));
        }
        let mut seen = vec![false; count];
        for source in std::iter::once(first).chain(sources) {
            checkpoint()?;
            if !Arc::ptr_eq(&first.state, &source.state)
                || source.compartment_index >= count
                || seen[source.compartment_index]
            {
                return Err(FissionThermalError::Invalid(
                    "stochastic movement cells must share one complete accepted history",
                ));
            }
            seen[source.compartment_index] = true;
        }
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
            "stochastic movement requires the complete cohort",
        ))
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
        Self::validate_sources(sources.iter().copied(), &mut || {
            checkpoint().map_err(Into::into)
        })?;
        if temperatures_k.len() != sources.len() {
            return Err(FissionThermalError::Invalid(
                "one temperature per stochastic compartment required",
            ));
        }
        let first = sources[0];
        first.check_clock(time_s)?;
        let mut checkpoint_error = None;
        let trial = first.state.trial_to(end_s, budget, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                DecayError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let (state, report) = trial?;
        let state = Arc::new(state);
        let report = Arc::new(report);
        let mut candidates = Vec::with_capacity(sources.len());
        for source in sources {
            checkpoint()?;
            candidates.push((
                Self {
                    cell_id: source.cell_id,
                    material_id: source.material_id,
                    compartment_index: source.compartment_index,
                    state: Arc::clone(&state),
                },
                StochasticDecayTransportCellStep {
                    compartment_index: source.compartment_index,
                    transport: Arc::clone(&report),
                },
            ));
        }
        Ok(candidates)
    }

    fn released(step: &Self::Step) -> EnergyChannels {
        step.decay().released
    }
    fn energy_roundoff_j(step: &Self::Step) -> f64 {
        step.decay().energy_roundoff_j
    }
    fn decays(step: &Self::Step) -> f64 {
        step.decay().events as f64
    }
    fn fissions(_: &Self::Step) -> f64 {
        0.0
    }
    fn captures(_: &Self::Step) -> f64 {
        0.0
    }
    fn cumulative_released(&self) -> Result<EnergyChannels, FissionThermalError> {
        Ok(self.state.inventories()[self.compartment_index].cumulative_released())
    }
}

#[cfg(test)]
mod tests;

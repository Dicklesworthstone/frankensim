//! Integer radioactive supply, movement and heat in a complete cell cohort.
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
use crate::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget, StochasticDecayProductionTransportState,
    StochasticDecayProductionTransportStep,
};
use crate::decay::stochastic::transport::{
    StochasticDecayTransportBudget, StochasticDecayTransportState, StochasticDecayTransportStep,
};
use crate::decay::stochastic::{StochasticDecayInventory, StochasticDecayStep};
use crate::{EnergyChannels, FissionError};

mod cohort_sealed {
    pub trait Sealed {}
}

/// Admitted stochastic compartment histories supported by one thermal adapter.
///
/// The implementations are sealed to the original transport owners. They retain
/// fixed compartment identities and return one local decay account per cell.
/// This lets supplied and unsupplied cohorts share binding, event clipping,
/// material admission, report mapping and all-or-nothing thermal publication.
pub trait StochasticDecayCohort: cohort_sealed::Sealed + Clone + std::fmt::Debug {
    type Budget: Copy + std::fmt::Debug;
    type Step: Clone + std::fmt::Debug;

    fn time_s(&self) -> f64;
    fn inventories(&self) -> &[StochasticDecayInventory];
    fn next_transition_time_s(&self) -> Option<f64>;
    fn energy_tolerance_j(budget: Self::Budget) -> f64;
    fn trial_cohort_to(
        &self,
        end_s: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, Self::Step), DecayError>;
    fn transport_report(step: &Self::Step) -> &StochasticDecayTransportStep;
    fn produced_nuclei(step: &Self::Step, compartment: usize) -> u64;
}

impl cohort_sealed::Sealed for StochasticDecayTransportState {}

impl StochasticDecayCohort for StochasticDecayTransportState {
    type Budget = StochasticDecayTransportBudget;
    type Step = StochasticDecayTransportStep;

    fn time_s(&self) -> f64 {
        StochasticDecayTransportState::time_s(self)
    }
    fn inventories(&self) -> &[StochasticDecayInventory] {
        StochasticDecayTransportState::inventories(self)
    }
    fn next_transition_time_s(&self) -> Option<f64> {
        StochasticDecayTransportState::next_event_time_s(self)
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.decay.energy_tolerance_j
    }
    fn trial_cohort_to(
        &self,
        end_s: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, Self::Step), DecayError> {
        StochasticDecayTransportState::trial_to(self, end_s, budget, checkpoint)
    }
    fn transport_report(step: &Self::Step) -> &StochasticDecayTransportStep {
        step
    }
    fn produced_nuclei(_: &Self::Step, _: usize) -> u64 {
        0
    }
}

impl cohort_sealed::Sealed for StochasticDecayProductionTransportState {}

impl StochasticDecayCohort for StochasticDecayProductionTransportState {
    type Budget = StochasticDecayProductionTransportBudget;
    type Step = StochasticDecayProductionTransportStep;

    fn time_s(&self) -> f64 {
        StochasticDecayProductionTransportState::time_s(self)
    }
    fn inventories(&self) -> &[StochasticDecayInventory] {
        StochasticDecayProductionTransportState::inventories(self)
    }
    fn next_transition_time_s(&self) -> Option<f64> {
        StochasticDecayProductionTransportState::next_transition_time_s(self)
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        budget.transport.decay.energy_tolerance_j
    }
    fn trial_cohort_to(
        &self,
        end_s: f64,
        budget: Self::Budget,
        checkpoint: &mut dyn FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, Self::Step), DecayError> {
        StochasticDecayProductionTransportState::trial_to(self, end_s, budget, checkpoint)
    }
    fn transport_report(step: &Self::Step) -> &StochasticDecayTransportStep {
        &step.transport
    }
    fn produced_nuclei(step: &Self::Step, compartment: usize) -> u64 {
        // The core admits the complete interval's birth count at or below 2^53.
        step.produced_nuclei[compartment].iter().copied().sum()
    }
}

/// One thermal cell's immutable view of the complete stochastic cohort.
/// Construct all handles together with [`Self::bind`], then advance them together.
#[derive(Clone, Debug)]
pub struct StochasticDecayCohortCellSource<C: StochasticDecayCohort> {
    cell_id: CellId,
    material_id: MaterialId,
    compartment_index: usize,
    state: Arc<C>,
}

/// An original transport cohort bound to a thermal cell.
pub type StochasticDecayTransportCellSource =
    StochasticDecayCohortCellSource<StochasticDecayTransportState>;

/// A supplied transport cohort bound to a thermal cell.
pub type StochasticDecayProductionTransportCellSource =
    StochasticDecayCohortCellSource<StochasticDecayProductionTransportState>;

/// Sampled movement, radioactive transformation and original thermal storage.
pub type StochasticDecayTransportHeating = NuclearHeating<StochasticDecayTransportCellSource>;

/// Sampled supply, movement, decay and original thermal storage in one history.
pub type StochasticDecayProductionTransportHeating =
    NuclearHeating<StochasticDecayProductionTransportCellSource>;

/// The original compartment emission account plus the shared movement report.
#[derive(Clone, Debug)]
pub struct StochasticDecayCohortCellStep<C: StochasticDecayCohort> {
    compartment_index: usize,
    report: Arc<C::Step>,
}

pub type StochasticDecayTransportCellStep =
    StochasticDecayCohortCellStep<StochasticDecayTransportState>;
pub type StochasticDecayProductionTransportCellStep =
    StochasticDecayCohortCellStep<StochasticDecayProductionTransportState>;

impl<C: StochasticDecayCohort> StochasticDecayCohortCellStep<C> {
    pub fn decay(&self) -> &StochasticDecayStep {
        &self.transport().decay[self.compartment_index]
    }
    pub fn transport(&self) -> &StochasticDecayTransportStep {
        C::transport_report(&self.report)
    }
    /// The entire shared interval report, including supply for a produced cohort.
    pub fn report(&self) -> &C::Step {
        &self.report
    }
    pub fn compartment_index(&self) -> usize {
        self.compartment_index
    }
}

impl<C: StochasticDecayCohort> StochasticDecayCohortCellSource<C> {
    /// Bind one distinct original thermal cell to every compartment. Rates are
    /// already authored per-nucleus hazards; binding does not infer carrier flow.
    pub fn bind(
        state: Arc<C>,
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

    pub fn state(&self) -> &Arc<C> {
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

impl<C: StochasticDecayCohort> sealed::Sealed for StochasticDecayCohortCellSource<C> {}

impl<C: StochasticDecayCohort> NuclearHeatSource for StochasticDecayCohortCellSource<C> {
    type Budget = C::Budget;
    type Step = StochasticDecayCohortCellStep<C>;

    fn cell_id(&self) -> CellId {
        self.cell_id
    }
    fn material_id(&self) -> MaterialId {
        self.material_id
    }
    fn energy_tolerance_j(budget: Self::Budget) -> f64 {
        C::energy_tolerance_j(budget)
    }
    fn next_transition(&self, time_s: f64) -> Option<f64> {
        // Complete-cohort validation guarantees this representative exists.
        // The thermal owner takes the union of source boundaries: advertising
        // the same shared event from every cell would rescan all clocks C times.
        if self.compartment_index == 0 {
            self.state
                .next_transition_time_s()
                .filter(|&next| next > time_s)
        } else {
            None
        }
    }

    fn admit(&self, time_s: f64, _: f64, budget: Self::Budget) -> Result<(), FissionThermalError> {
        self.check_clock(time_s)?;
        if self.compartment_index == 0 {
            self.state.trial_cohort_to(time_s, budget, &mut || Ok(()))?;
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
        let trial = first.state.trial_cohort_to(end_s, budget, &mut || {
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
                StochasticDecayCohortCellStep {
                    compartment_index: source.compartment_index,
                    report: Arc::clone(&report),
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
    fn produced_nuclei(step: &Self::Step) -> f64 {
        C::produced_nuclei(step.report(), step.compartment_index) as f64
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

#[cfg(test)]
mod production_tests;

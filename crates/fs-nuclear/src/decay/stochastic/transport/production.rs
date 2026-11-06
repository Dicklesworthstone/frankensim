//! Scheduled Poisson supply, movement and decay on one compensated event race.
//!
//! An arriving nucleus can move, decay, create a daughter and leave the tracked
//! cohort during the same observation interval. These operations share their
//! actual event times, rather than advancing separate owners to rounded event
//! endpoints. Original local decay histories keep their energy; neither supply
//! nor movement invents heat, carrier mass or enthalpy.
//!
//! Schedules supply independent Poisson intensities, with the same explicit
//! interpretation as the standalone production owner. Movement changes win exact
//! ties, then source changes, original decay/movement order and births in
//! compartment/node order. Unchanged rates and observations never redraw clocks. Clone
//! retains the entire cohort, supply history and all future random events.

use std::sync::Arc;

use super::super::production::PoissonSource;
use super::schedule::StochasticDecayTransportSchedule;
use super::{
    DecayError, MAX_STOCHASTIC_NUCLEI, StochasticDecayInventory, StochasticDecayTransportBudget,
    StochasticDecayTransportState, StochasticDecayTransportStep, ceiling,
};
use crate::decay::production::DecayProductionSchedule;

/// Work and energy allowances are shared across the complete cohort and trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StochasticDecayProductionTransportBudget {
    pub transport: StochasticDecayTransportBudget,
    /// Maximum total external arrivals, at most 2^53, across all compartments.
    pub maximum_births: u64,
    /// Maximum traversed source transitions, including unchanged-rate entries.
    pub maximum_schedule_changes: usize,
}

impl StochasticDecayProductionTransportBudget {
    fn validate(self) -> Result<(), DecayError> {
        self.transport.validate()?;
        if self.maximum_births > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget(
                "stochastic production birth-count ceiling",
            ));
        }
        Ok(())
    }
}

/// The original local decay/movement accounts plus externally supplied nuclei.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayProductionTransportStep {
    pub transport: StochasticDecayTransportStep,
    /// Compartment order, then original nuclide order; absent sources have zeros.
    pub produced_nuclei: Vec<Vec<u64>>,
    pub births: u64,
    pub schedule_changes: usize,
}

/// Complete produced and transported inventory, including pending source clocks.
///
/// Each bound source adds one arrival clock per nuclide and one schedule slot to
/// the original construction clock limit. Trials charge those same slots and
/// source rescheduling to the existing clock-work budget. Source draws use their
/// original compartment's decay stream; the distinct movement stream is retained.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayProductionTransportState {
    transport: StochasticDecayTransportState,
    sources: Vec<Option<PoissonSource>>,
    source_clock_count: usize,
}

/// Optional supply context for the original transport trial kernel. The caller
/// owns candidate clocks and interval counts; failure never publishes them.
pub(super) struct ProductionTrial<'a> {
    pub(super) sources: &'a mut [Option<PoissonSource>],
    pub(super) source_clock_count: usize,
    pub(super) maximum_births: u64,
    pub(super) maximum_schedule_changes: usize,
    pub(super) produced_nuclei: Vec<Vec<u64>>,
    pub(super) births: u64,
    pub(super) schedule_changes: usize,
}

impl StochasticDecayProductionTransportState {
    /// Bind one optional independent Poisson source per existing compartment.
    ///
    /// Schedules must share the exact original model Arc. Rates active at the
    /// cohort's current time prepare future arrivals; past production is not
    /// replayed. Existing decay/movement clocks survive. `None` adds no clocks,
    /// random draws or source work; an all-None collection retains transport's
    /// original state and reports exactly.
    pub fn poisson(
        mut transport: StochasticDecayTransportState,
        schedules: Vec<Option<Arc<DecayProductionSchedule>>>,
    ) -> Result<Self, DecayError> {
        if schedules.len() != transport.inventories.len() {
            return Err(DecayError::Invalid(
                "one optional Poisson schedule per stochastic compartment required",
            ));
        }
        let source_clock_count = transport
            .nodes
            .checked_add(1)
            .and_then(|slots| slots.checked_mul(schedules.iter().filter(|s| s.is_some()).count()))
            .ok_or(DecayError::Budget("stochastic production transport clocks"))?;
        let total_clocks = transport
            .clock_count
            .checked_add(source_clock_count)
            .ok_or(DecayError::Budget("stochastic production transport clocks"))?;
        if total_clocks > transport.limits.maximum_clocks {
            return Err(DecayError::Budget("stochastic production transport clocks"));
        }
        let mut sources = Vec::with_capacity(schedules.len());
        for (schedule, inventory) in schedules.into_iter().zip(&mut transport.inventories) {
            sources.push(match schedule {
                Some(schedule) => Some(PoissonSource::new(inventory, schedule)?),
                None => None,
            });
        }
        Ok(Self {
            transport,
            sources,
            source_clock_count,
        })
    }

    /// Bind or replace movement rates at the current time without touching any
    /// pending source clock or birth counter. The original combined storage bound
    /// includes one movement schedule slot, which is reused on later rebindings.
    pub fn with_transport_schedule(
        mut self,
        schedule: Arc<StochasticDecayTransportSchedule>,
    ) -> Result<Self, DecayError> {
        let total_clocks = self
            .transport
            .clock_count
            .checked_add(usize::from(self.transport.transport_schedule().is_none()))
            .and_then(|count| count.checked_add(self.source_clock_count))
            .ok_or(DecayError::Budget("stochastic production transport clocks"))?;
        if total_clocks > self.transport.limits.maximum_clocks {
            return Err(DecayError::Budget("stochastic production transport clocks"));
        }
        self.transport = self.transport.with_transport_schedule(schedule)?;
        Ok(self)
    }

    #[must_use]
    pub fn transport_state(&self) -> &StochasticDecayTransportState {
        &self.transport
    }

    #[must_use]
    pub fn inventories(&self) -> &[StochasticDecayInventory] {
        self.transport.inventories()
    }

    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.transport.time_s()
    }

    /// None means an absent source or an out-of-range compartment.
    #[must_use]
    pub fn schedule(&self, compartment: usize) -> Option<&Arc<DecayProductionSchedule>> {
        self.sources
            .get(compartment)?
            .as_ref()
            .map(PoissonSource::schedule)
    }

    /// Original-order external supply history; initial and transferred nuclei
    /// are excluded. None means an absent source or out-of-range compartment.
    #[must_use]
    pub fn cumulative_produced_nuclei(&self, compartment: usize) -> Option<&[u64]> {
        self.sources
            .get(compartment)?
            .as_ref()
            .map(PoissonSource::cumulative_produced_nuclei)
    }

    /// First representable endpoint containing any supply, movement or decay
    /// event, or an authored source transition. No observation consumes a draw.
    #[must_use]
    pub fn next_transition_time_s(&self) -> Option<f64> {
        if self.source_clock_count == 0 {
            self.transport.next_event_time_s()
        } else {
            self.transport
                .next_event_with_sources(&self.sources)
                .map(|(_, time)| ceiling(time))
        }
    }

    /// Evaluate a complete candidate with exact integer supply/outflow balance.
    /// Every event uses its compensated time before endpoint conversion. Work,
    /// energy, population, counter, random-stream and callback refusals retain
    /// the entire original cohort and all pending source clocks.
    pub fn trial_to(
        &self,
        end_time_s: f64,
        budget: StochasticDecayProductionTransportBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, StochasticDecayProductionTransportStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        let mut sources = self.sources.clone();
        let mut production = ProductionTrial {
            sources: &mut sources,
            source_clock_count: self.source_clock_count,
            maximum_births: budget.maximum_births,
            maximum_schedule_changes: budget.maximum_schedule_changes,
            produced_nuclei: vec![vec![0; self.transport.nodes]; self.sources.len()],
            births: 0,
            schedule_changes: 0,
        };
        let context = if self.source_clock_count == 0 {
            None
        } else {
            Some(&mut production)
        };
        let (transport, report) = self.transport.trial_with_sources(
            end_time_s,
            budget.transport,
            context,
            &mut checkpoint,
        )?;
        let ProductionTrial {
            produced_nuclei,
            births,
            schedule_changes,
            ..
        } = production;
        Ok((
            Self {
                transport,
                sources,
                source_clock_count: self.source_clock_count,
            },
            StochasticDecayProductionTransportStep {
                transport: report,
                produced_nuclei,
                births,
                schedule_changes,
            },
        ))
    }

    /// Publish only the entirely accepted supply/movement/decay candidate.
    pub fn advance_to(
        &mut self,
        end_time_s: f64,
        budget: StochasticDecayProductionTransportBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<StochasticDecayProductionTransportStep, DecayError> {
        let (next, report) = self.trial_to(end_time_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;

//! Prescribed Poisson supply on the original integer radioactive trajectory.
//!
//! The existing [`DecayProductionSchedule`] supplies independent per-nuclide
//! arrival intensities in nuclei/s. This is an explicit Poisson-source model,
//! not a claim that every prescribed feed or exact dose has Poisson fluctuations.
//! Each arrival retains its actual time, so it and its daughters may decay before
//! the next observation. Supply itself adds no heat or carrier enthalpy.
//!
//! Pending birth, decay and authored rate-change events share one chronological
//! race. Only a changed source rate redraws its clock at a schedule boundary, by
//! memorylessness. Observation endpoints and unchanged rates draw nothing.
//! Compensated event times are compared before rounding to an endpoint. Exact
//! ties process schedule changes first, then decay, then birth; node ties retain
//! original order. In particular, the previous supply applies on a half-open
//! interval ending at its authored change. Clone retains the entire history.

use std::sync::Arc;

use super::{
    Accumulator, DecayError, EventClock, MAX_STOCHASTIC_NUCLEI, StochasticDecayBudget,
    StochasticDecayInventory, StochasticDecayStep, before, ceiling, checked, defect_bound,
    export_channels, increment, nonnegative, sample_clock,
};
use crate::decay::production::DecayProductionSchedule;
use crate::{EnergyChannels, channel_values};

/// Separate bounds for sampled births, transformations and authored transitions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StochasticDecayProductionBudget {
    pub decay: StochasticDecayBudget,
    /// Maximum arrivals in this trial, at most 2^53; zero permits no births.
    pub maximum_births: u64,
    /// Maximum traversed changes, including changes with identical rates.
    pub maximum_schedule_changes: usize,
}

impl StochasticDecayProductionBudget {
    fn validate(self) -> Result<(), DecayError> {
        self.decay.validate()?;
        if self.maximum_births > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget(
                "stochastic production birth-count ceiling",
            ));
        }
        Ok(())
    }
}

/// Accepted integer arrivals and the original per-branch emission account.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayProductionStep {
    /// Arrivals in the original model's nuclide order; excludes initial nuclei.
    pub produced_nuclei: Vec<u64>,
    pub births: u64,
    pub schedule_changes: usize,
    /// Its RNG draw count includes source clocks on the shared inventory stream.
    pub decay: StochasticDecayStep,
}

/// A complete produced-inventory checkpoint, including all future event clocks.
///
/// The original closed inventory still owns decay branches, energy and randomness.
/// Counts remain integers, with total population bounded by 2^53 after each birth.
/// Construction and trial storage are O(nodes); event work scans at most twice
/// the original bounded node list and one parent's branches, without allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayProduction {
    inventory: StochasticDecayInventory,
    schedule: Arc<DecayProductionSchedule>,
    rates_nuclei_s: Vec<f64>,
    clocks: Vec<EventClock>,
    produced_nuclei: Vec<u64>,
    total_population: u64,
    next_change: usize,
}

#[derive(Clone, Copy)]
enum Event {
    Change,
    Decay(usize),
    Birth(usize),
}

impl StochasticDecayProduction {
    /// Bind independent Poisson arrival intensities to an existing trajectory.
    ///
    /// The schedule must share the exact original model Arc. Rates active at the
    /// inventory's current time initialize future clocks; earlier supply is not
    /// replayed. The original decay clocks survive. Zero supply consumes no draws.
    pub fn poisson(
        inventory: StochasticDecayInventory,
        schedule: Arc<DecayProductionSchedule>,
    ) -> Result<Self, DecayError> {
        if !Arc::ptr_eq(inventory.model(), schedule.model()) {
            return Err(DecayError::Invalid(
                "stochastic production schedule must share the original decay model",
            ));
        }
        let time_s = inventory.time_s();
        let rates_nuclei_s = schedule.rates_at(time_s)?.to_vec();
        let nodes = rates_nuclei_s.len();
        let total_population = inventory.populations.iter().copied().sum();
        let next_change = schedule
            .changes()
            .partition_point(|change| change.time_s <= time_s);
        let mut result = Self {
            inventory,
            schedule,
            rates_nuclei_s,
            clocks: vec![EventClock::Inactive; nodes],
            produced_nuclei: vec![0; nodes],
            total_population,
            next_change,
        };
        for node in 0..nodes {
            result.reschedule_birth(node, Accumulator::new(time_s))?;
        }
        Ok(result)
    }

    #[must_use]
    pub fn inventory(&self) -> &StochasticDecayInventory {
        &self.inventory
    }

    #[must_use]
    pub fn schedule(&self) -> &Arc<DecayProductionSchedule> {
        &self.schedule
    }

    #[must_use]
    pub fn cumulative_produced_nuclei(&self) -> &[u64] {
        &self.produced_nuclei
    }

    /// First representable endpoint containing a birth, decay or rate change.
    /// A thermal owner must clip at births even when they release no energy:
    /// they alter the future decay history inside the following interval.
    #[must_use]
    pub fn next_transition_time_s(&self) -> Option<f64> {
        self.next_event().map(|(_, time)| ceiling(time))
    }

    /// Evaluate an atomic candidate through an absolute, nondecreasing endpoint.
    ///
    /// Births and all descendant decays use their sampled ages within the same
    /// interval. Endpoint changes are installed once and prepare future clocks;
    /// zero-duration observations are identities. Every refusal retains all
    /// original populations, supply counters, emitted energy and random state.
    pub fn trial_to(
        &self,
        end_time_s: f64,
        budget: StochasticDecayProductionBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, StochasticDecayProductionStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        let start_time_s = self.inventory.time_s();
        if !nonnegative(end_time_s) || end_time_s < start_time_s {
            return Err(DecayError::Invalid(
                "stochastic production endpoint must be finite and not precede the current time",
            ));
        }
        let mut next = self.clone();
        let mut report = StochasticDecayProductionStep {
            produced_nuclei: vec![0; self.rates_nuclei_s.len()],
            births: 0,
            schedule_changes: 0,
            decay: StochasticDecayStep {
                start_time_s,
                end_time_s,
                duration_s: end_time_s - start_time_s,
                events: 0,
                decays: vec![0; self.rates_nuclei_s.len()],
                branch_decays: vec![0; self.inventory.model.branches.len()],
                released: EnergyChannels::default(),
                energy_roundoff_j: 0.0,
                rng_draws: 0,
            },
        };
        let mut released = [Accumulator::new(0.0); 3];
        let mut energy_error = Accumulator::new(0.0);
        while let Some((event, time)) = next.next_event() {
            if ceiling(time) > end_time_s || end_time_s == start_time_s {
                break;
            }
            checkpoint()?;
            match event {
                Event::Change => {
                    if report.schedule_changes == budget.maximum_schedule_changes {
                        return Err(DecayError::Budget("stochastic production schedule changes"));
                    }
                    next.change_rates(time)?;
                    report.schedule_changes += 1;
                }
                Event::Decay(parent) => {
                    if report.decay.events == budget.decay.maximum_events {
                        return Err(DecayError::Budget("stochastic decay events"));
                    }
                    next.inventory.apply_decay(
                        parent,
                        time,
                        &mut report.decay,
                        &mut released,
                        &mut energy_error,
                    )?;
                }
                Event::Birth(node) => {
                    if report.births == budget.maximum_births {
                        return Err(DecayError::Budget("stochastic production births"));
                    }
                    if next.total_population == MAX_STOCHASTIC_NUCLEI {
                        return Err(DecayError::Budget("stochastic population ceiling"));
                    }
                    increment(&mut next.total_population)?;
                    increment(&mut next.inventory.populations[node])?;
                    increment(&mut next.produced_nuclei[node])?;
                    increment(&mut report.produced_nuclei[node])?;
                    increment(&mut report.births)?;
                    next.inventory.reschedule(node, time)?;
                    next.reschedule_birth(node, time)?;
                }
            }
        }
        report.decay.released = export_channels(released, &mut energy_error)?;
        checked(channel_values(next.inventory.released).total_j())?;
        report.decay.energy_roundoff_j = defect_bound(energy_error)?;
        if report.decay.energy_roundoff_j > budget.decay.energy_tolerance_j {
            return Err(DecayError::Budget("stochastic emission energy roundoff"));
        }
        report.decay.rng_draws = next
            .inventory
            .rng
            .index
            .checked_sub(self.inventory.rng.index)
            .ok_or(DecayError::Arithmetic)?;
        next.inventory.time_s = end_time_s;
        checkpoint()?;
        Ok((next, report))
    }

    /// Publish only an entirely admitted produced-inventory trajectory.
    pub fn advance_to(
        &mut self,
        end_time_s: f64,
        budget: StochasticDecayProductionBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<StochasticDecayProductionStep, DecayError> {
        let (next, report) = self.trial_to(end_time_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn next_event(&self) -> Option<(Event, Accumulator)> {
        let mut next = self
            .schedule
            .changes()
            .get(self.next_change)
            .map(|change| (Event::Change, Accumulator::new(change.time_s)));
        if let Some((node, time)) = self.inventory.next_event()
            && next.is_none_or(|(_, earliest)| before(time, earliest))
        {
            next = Some((Event::Decay(node), time));
        }
        for (node, clock) in self.clocks.iter().enumerate() {
            if let EventClock::Finite(time) = clock
                && next.is_none_or(|(_, earliest)| before(*time, earliest))
            {
                next = Some((Event::Birth(node), *time));
            }
        }
        next
    }

    fn reschedule_birth(&mut self, node: usize, origin: Accumulator) -> Result<(), DecayError> {
        self.clocks[node] = sample_clock(
            &mut self.inventory.rng,
            self.rates_nuclei_s[node],
            1,
            origin,
        )?;
        Ok(())
    }

    fn change_rates(&mut self, time: Accumulator) -> Result<(), DecayError> {
        let change = &self.schedule.changes()[self.next_change];
        for (node, &rate) in change.rates_nuclei_s.iter().enumerate() {
            if rate != self.rates_nuclei_s[node] {
                self.clocks[node] = sample_clock(&mut self.inventory.rng, rate, 1, time)?;
                self.rates_nuclei_s[node] = rate;
            }
        }
        self.next_change += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests;

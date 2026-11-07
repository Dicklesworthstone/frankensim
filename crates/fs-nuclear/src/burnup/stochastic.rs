//! Integer capture activation and radioactive decay on the original nuclear data.
//!
//! Every admitted channel removes one nucleus and produces one explicit retained
//! daughter. Independent reaction clocks compete at their compensated actual
//! times, so newly activated nuclei can decay or capture again inside the same
//! observation interval. The original BurnupModel owns microscopic rates,
//! half-lives, product identities and immediate event energy. Capture releases
//! its prompt energy; subsequent radioactive events own their delayed energy.
//!
//! The existing radioactive clock sampler, counter-based stream checkpoint and
//! endpoint convention are reused. Observation boundaries consume no randomness.
//! Authored forcing changes win exact ties; strictly earlier sub-ULP reactions
//! remain earlier. Unchanged rates retain their pending clocks. Clone preserves
//! populations, clocks, counters and energy on the same ISA/toolchain profile.
//!
//! This is a finite-population, one-product activation lane. Fission and expected
//! fractional/multiple yields require joint event outcomes and are refused. It
//! does not infer neutron feedback, carrier movement, emitted-particle mass,
//! event spectra, or temperature changes beyond the supplied forcing history.

use std::sync::Arc;

use super::{
    BurnupError, BurnupModel, BurnupRate, add, channel_array, checked, nonnegative_bound,
    projection_error,
};
use crate::decay::DecayError;
use crate::decay::stochastic::{EventClock, MAX_STOCHASTIC_NUCLEI, before, ceiling, sample_clock};
pub use crate::decay::stochastic::{StreamCheckpoint, StreamKey};
use crate::{Accumulator, EnergyChannels, FissionError, channel_values, nonnegative, sum};

mod schedule;
pub use schedule::{StochasticActivationChange, StochasticActivationSchedule};

/// Whole-interval bounds on nuclear events, clock work and numerical energy error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StochasticActivationBudget {
    /// Capture plus decay events in this trial; at most 2^53.
    pub maximum_events: u64,
    /// Reaction/schedule slots scanned, coefficients visited and clocks resampled.
    pub maximum_clock_checks: u64,
    /// Noncancelling energy-account rounding allowance [J], not data uncertainty.
    pub energy_tolerance_j: f64,
}

impl StochasticActivationBudget {
    fn validate(self) -> Result<(), BurnupError> {
        if self.maximum_events > MAX_STOCHASTIC_NUCLEI {
            return Err(BurnupError::Budget(
                "stochastic activation event-count ceiling",
            ));
        }
        if !nonnegative(self.energy_tolerance_j) {
            return Err(BurnupError::Invalid(
                "stochastic activation energy tolerance must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

/// Accepted original-order integer reactions and their immediate energy channels.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticActivationStep {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub duration_s: f64,
    pub reaction_counts: Vec<u64>,
    pub captures: u64,
    pub decays: u64,
    pub events: u64,
    /// Installed forcing frames, including frames with unchanged rates.
    pub schedule_changes: usize,
    pub capture_released: EnergyChannels,
    pub decay_released: EnergyChannels,
    pub released: EnergyChannels,
    pub energy_roundoff_j: f64,
    pub rng_draws: u64,
    pub clock_checks: u64,
}

/// One finite physical trajectory and its pending capture/decay events.
///
/// All populations are explicit integers with total at most 2^53. The immutable
/// forcing also retains parent-to-channel adjacency. Trials clone mutable state
/// and allocate reports once; individual reactions and frame changes allocate
/// nothing and never recollapse cross sections or copy a rate vector.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticActivationState {
    schedule: Arc<StochasticActivationSchedule>,
    position: usize,
    populations: Vec<u64>,
    reaction_counts: Vec<u64>,
    clocks: Vec<EventClock>,
    rng: StreamCheckpoint,
    time_s: f64,
    released: [Accumulator; 3],
    capture_released: [Accumulator; 3],
    decay_released: [Accumulator; 3],
}

#[derive(Clone, Copy)]
enum Event {
    Change,
    Reaction(usize),
}

impl StochasticActivationState {
    pub fn new(
        schedule: Arc<StochasticActivationSchedule>,
        populations: Vec<u64>,
        key: StreamKey,
    ) -> Result<Self, BurnupError> {
        Self::new_at(schedule, populations, key, 0.0)
    }

    /// Start at the supplied absolute time; earlier forcing or reactions are not replayed.
    pub fn new_at(
        schedule: Arc<StochasticActivationSchedule>,
        populations: Vec<u64>,
        key: StreamKey,
        time_s: f64,
    ) -> Result<Self, BurnupError> {
        let position = schedule.position_at(time_s)?;
        if populations.len() != schedule.model().nuclides().len() {
            return Err(BurnupError::Invalid(
                "one integer activation population per original nuclide required",
            ));
        }
        population_total(&populations)?;
        let reactions = schedule.model().reactions().len();
        let mut state = Self {
            schedule,
            position,
            populations,
            reaction_counts: vec![0; reactions],
            clocks: vec![EventClock::Inactive; reactions],
            rng: StreamCheckpoint::current(key, 0),
            time_s,
            released: [Accumulator::new(0.0); 3],
            capture_released: [Accumulator::new(0.0); 3],
            decay_released: [Accumulator::new(0.0); 3],
        };
        for reaction in 0..reactions {
            state.reschedule(reaction, Accumulator::new(time_s))?;
        }
        Ok(state)
    }

    #[must_use]
    pub fn model(&self) -> &Arc<BurnupModel> {
        self.schedule.model()
    }

    #[must_use]
    pub fn schedule(&self) -> &Arc<StochasticActivationSchedule> {
        &self.schedule
    }

    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    #[must_use]
    pub fn populations(&self) -> &[u64] {
        &self.populations
    }

    #[must_use]
    pub fn cumulative_reaction_counts(&self) -> &[u64] {
        &self.reaction_counts
    }

    #[must_use]
    pub fn rng_checkpoint(&self) -> StreamCheckpoint {
        self.rng
    }

    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        channel_values(self.released)
    }

    #[must_use]
    pub fn cumulative_capture_released(&self) -> EnergyChannels {
        channel_values(self.capture_released)
    }

    #[must_use]
    pub fn cumulative_decay_released(&self) -> EnergyChannels {
        channel_values(self.decay_released)
    }

    /// First representable endpoint containing a reaction or authored rate change.
    #[must_use]
    pub fn next_transition_time_s(&self) -> Option<f64> {
        self.next_event().map(|(_, time)| ceiling(time))
    }

    /// Bind new forcing at the actual current time, retaining all physical history.
    /// Only changed rates redraw clocks. Earlier frames are not replayed.
    pub fn with_schedule(
        mut self,
        schedule: Arc<StochasticActivationSchedule>,
    ) -> Result<Self, BurnupError> {
        if !Arc::ptr_eq(self.model(), schedule.model())
            || self.schedule.limits() != schedule.limits()
        {
            return Err(BurnupError::Invalid(
                "activation forcing belongs to a different original model or limits",
            ));
        }
        let position = schedule.position_at(self.time_s)?;
        for reaction in 0..self.clocks.len() {
            let rate = schedule.rates(position)[reaction];
            if rate != self.schedule.rates(self.position)[reaction] {
                self.reschedule_at_rate(reaction, rate, Accumulator::new(self.time_s))?;
            }
        }
        self.schedule = schedule;
        self.position = position;
        Ok(self)
    }

    /// Evaluate one atomic candidate; every callback error keeps its original cause.
    pub fn trial_to(
        &self,
        end_time_s: f64,
        budget: StochasticActivationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, StochasticActivationStep), BurnupError> {
        checkpoint()?;
        budget.validate()?;
        if !nonnegative(end_time_s) || end_time_s < self.time_s {
            return Err(BurnupError::Invalid(
                "activation endpoint must be finite and nondecreasing",
            ));
        }
        let initial_total = population_total(&self.populations)?;
        let mut next = self.clone();
        let mut report = StochasticActivationStep {
            start_time_s: self.time_s,
            end_time_s,
            duration_s: end_time_s - self.time_s,
            reaction_counts: vec![0; self.clocks.len()],
            captures: 0,
            decays: 0,
            events: 0,
            schedule_changes: 0,
            capture_released: EnergyChannels::default(),
            decay_released: EnergyChannels::default(),
            released: EnergyChannels::default(),
            energy_roundoff_j: 0.0,
            rng_draws: 0,
            clock_checks: 0,
        };
        let mut released = [Accumulator::new(0.0); 3];
        let mut capture_released = [Accumulator::new(0.0); 3];
        let mut decay_released = [Accumulator::new(0.0); 3];
        let mut energy_error = Accumulator::new(0.0);
        let slots = self
            .clocks
            .len()
            .checked_add(1)
            .ok_or(BurnupError::Budget("stochastic activation clock storage"))?;
        if end_time_s > self.time_s {
            loop {
                checkpoint()?;
                charge_clocks(&mut report, slots, budget.maximum_clock_checks)?;
                let Some((event, time)) = next.next_event() else {
                    break;
                };
                if ceiling(time) > end_time_s {
                    break;
                }
                match event {
                    Event::Change => {
                        next.apply_change(time, &mut report, budget, &mut checkpoint)?;
                    }
                    Event::Reaction(reaction) => {
                        if report.events == budget.maximum_events {
                            return Err(BurnupError::Budget("stochastic activation events"));
                        }
                        next.apply_reaction(
                            reaction,
                            time,
                            &mut report,
                            &mut released,
                            &mut capture_released,
                            &mut decay_released,
                            &mut energy_error,
                            budget,
                            &mut checkpoint,
                        )?;
                    }
                }
            }
        }
        report.released = export_energy(released, &mut energy_error)?;
        report.capture_released = export_energy(capture_released, &mut energy_error)?;
        report.decay_released = export_energy(decay_released, &mut energy_error)?;
        for ((total, capture), decay) in channel_array(report.released)
            .into_iter()
            .zip(channel_array(report.capture_released))
            .zip(channel_array(report.decay_released))
        {
            add(
                &mut energy_error,
                checked(sum(&[total, -capture, -decay]))?.abs(),
            )?;
        }
        checked(next.cumulative_released().total_j())?;
        checked(next.cumulative_capture_released().total_j())?;
        checked(next.cumulative_decay_released().total_j())?;
        if population_total(&next.populations)? != initial_total
            || report.captures.checked_add(report.decays) != Some(report.events)
        {
            return Err(BurnupError::Arithmetic(
                "stochastic activation integer conservation",
            ));
        }
        report.energy_roundoff_j = nonnegative_bound(energy_error)?;
        if report.energy_roundoff_j > budget.energy_tolerance_j {
            return Err(BurnupError::Budget("stochastic activation energy roundoff"));
        }
        report.rng_draws =
            next.rng
                .index
                .checked_sub(self.rng.index)
                .ok_or(BurnupError::Arithmetic(
                    "stochastic activation random history",
                ))?;
        next.time_s = end_time_s;
        checkpoint()?;
        Ok((next, report))
    }

    /// Publish only a completely admitted candidate, including its pending clocks.
    pub fn advance_to(
        &mut self,
        end_time_s: f64,
        budget: StochasticActivationBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<StochasticActivationStep, BurnupError> {
        let (next, report) = self.trial_to(end_time_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn next_event(&self) -> Option<(Event, Accumulator)> {
        // A new forcing interval owns exact ties; original channel order breaks
        // reaction ties. Compensated low parts keep truly earlier events earlier.
        let mut next = self
            .schedule
            .changes()
            .get(self.position)
            .map(|change| (Event::Change, Accumulator::new(change.time_s)));
        for (reaction, clock) in self.clocks.iter().enumerate() {
            if let EventClock::Finite(time) = clock
                && next.is_none_or(|(_, earlier)| before(*time, earlier))
            {
                next = Some((Event::Reaction(reaction), *time));
            }
        }
        next
    }

    fn apply_change(
        &mut self,
        time: Accumulator,
        report: &mut StochasticActivationStep,
        budget: StochasticActivationBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        let position = self.position.checked_add(1).ok_or(BurnupError::Arithmetic(
            "stochastic activation forcing position",
        ))?;
        for reaction in 0..self.clocks.len() {
            if reaction % 256 == 0 {
                checkpoint()?;
            }
            let old = self.schedule.rates(self.position)[reaction];
            let rate = self.schedule.rates(position)[reaction];
            charge_clocks(
                report,
                1 + usize::from(rate != old),
                budget.maximum_clock_checks,
            )?;
            if rate != old {
                self.reschedule_at_rate(reaction, rate, time)?;
            }
        }
        self.position = position;
        report.schedule_changes =
            report
                .schedule_changes
                .checked_add(1)
                .ok_or(BurnupError::Arithmetic(
                    "stochastic activation forcing count",
                ))?;
        Ok(())
    }

    fn apply_reaction(
        &mut self,
        reaction: usize,
        time: Accumulator,
        report: &mut StochasticActivationStep,
        released: &mut [Accumulator; 3],
        capture_released: &mut [Accumulator; 3],
        decay_released: &mut [Accumulator; 3],
        energy_error: &mut Accumulator,
        budget: StochasticActivationBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        let data = &self.model().reactions()[reaction];
        let (parent, daughter, is_decay, energy) = (
            data.parent,
            data.products[0].nuclide,
            matches!(&data.rate, BurnupRate::Decay { .. }),
            data.energy_per_event,
        );
        self.populations[parent] =
            self.populations[parent]
                .checked_sub(1)
                .ok_or(BurnupError::Arithmetic(
                    "stochastic activation parent debit",
                ))?;
        increment(&mut self.populations[daughter])?;
        increment(&mut self.reaction_counts[reaction])?;
        increment(&mut report.reaction_counts[reaction])?;
        increment(&mut report.events)?;
        record_energy(&mut self.released, released, energy, energy_error)?;
        if is_decay {
            increment(&mut report.decays)?;
            record_energy(
                &mut self.decay_released,
                decay_released,
                energy,
                energy_error,
            )?;
        } else {
            increment(&mut report.captures)?;
            record_energy(
                &mut self.capture_released,
                capture_released,
                energy,
                energy_error,
            )?;
        }
        if parent == daughter {
            // A retained self-product leaves other propensities unchanged; only
            // the firing channel needs its next event sampled.
            charge_clocks(report, 1, budget.maximum_clock_checks)?;
            self.reschedule(reaction, time)?;
        } else {
            self.reschedule_outgoing(parent, time, report, budget, checkpoint)?;
            self.reschedule_outgoing(daughter, time, report, budget, checkpoint)?;
        }
        Ok(())
    }

    fn reschedule_outgoing(
        &mut self,
        node: usize,
        time: Accumulator,
        report: &mut StochasticActivationStep,
        budget: StochasticActivationBudget,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), BurnupError> {
        let count = self.schedule.outgoing(node).len();
        charge_clocks(report, count, budget.maximum_clock_checks)?;
        for offset in 0..count {
            if offset % 256 == 0 {
                checkpoint()?;
            }
            let reaction = self.schedule.outgoing(node)[offset];
            self.reschedule(reaction, time)?;
        }
        Ok(())
    }

    fn reschedule(&mut self, reaction: usize, time: Accumulator) -> Result<(), BurnupError> {
        let rate = self.schedule.rates(self.position)[reaction];
        self.reschedule_at_rate(reaction, rate, time)
    }

    fn reschedule_at_rate(
        &mut self,
        reaction: usize,
        rate: f64,
        time: Accumulator,
    ) -> Result<(), BurnupError> {
        let parent = self.model().reactions()[reaction].parent;
        self.clocks[reaction] = sample_clock(&mut self.rng, rate, self.populations[parent], time)
            .map_err(clock_error)?;
        Ok(())
    }
}

fn record_energy(
    lifetime: &mut [Accumulator; 3],
    interval: &mut [Accumulator; 3],
    energy: EnergyChannels,
    defects: &mut Accumulator,
) -> Result<(), BurnupError> {
    for ((lifetime, interval), value) in
        lifetime.iter_mut().zip(interval).zip(channel_array(energy))
    {
        add(defects, add(lifetime, value)?.abs())?;
        add(defects, add(interval, value)?.abs())?;
    }
    Ok(())
}

fn export_energy(
    values: [Accumulator; 3],
    defects: &mut Accumulator,
) -> Result<EnergyChannels, BurnupError> {
    let energy = channel_values(values);
    checked(energy.total_j())?;
    for value in values {
        add(defects, projection_error(value)?)?;
    }
    Ok(energy)
}

fn population_total(populations: &[u64]) -> Result<u64, BurnupError> {
    let mut total = 0_u64;
    for &count in populations {
        total = total.checked_add(count).ok_or(BurnupError::Budget(
            "stochastic activation population ceiling",
        ))?;
        if total > MAX_STOCHASTIC_NUCLEI {
            return Err(BurnupError::Budget(
                "stochastic activation population ceiling",
            ));
        }
    }
    Ok(total)
}

fn increment(value: &mut u64) -> Result<(), BurnupError> {
    *value = value.checked_add(1).ok_or(BurnupError::Arithmetic(
        "stochastic activation integer counter",
    ))?;
    Ok(())
}

fn charge_clocks(
    report: &mut StochasticActivationStep,
    count: usize,
    maximum: u64,
) -> Result<(), BurnupError> {
    let count = u64::try_from(count)
        .map_err(|_| BurnupError::Budget("stochastic activation clock work"))?;
    report.clock_checks = report
        .clock_checks
        .checked_add(count)
        .ok_or(BurnupError::Budget("stochastic activation clock work"))?;
    if report.clock_checks > maximum {
        return Err(BurnupError::Budget("stochastic activation clock work"));
    }
    Ok(())
}

fn clock_error(error: DecayError) -> BurnupError {
    match error {
        DecayError::Invalid(reason) => BurnupError::Invalid(reason),
        DecayError::Budget(reason) => BurnupError::Budget(reason),
        DecayError::Arithmetic => BurnupError::Arithmetic("stochastic activation event clock"),
        DecayError::Propagation(error) => BurnupError::Propagation(error),
        DecayError::Cancelled => BurnupError::Fission(FissionError::Cancelled),
    }
}

#[cfg(test)]
mod tests;

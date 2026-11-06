//! Finite, integer radioactive populations on the original decay model.
//!
//! Each occupied unstable nuclide has an exponential clock with hazard
//! `lambda * nuclei`. The earliest clock selects the next transformation;
//! the admitted branch probabilities select its daughter and energy. Only the
//! changed parent and daughter clocks are resampled, by memorylessness, at that
//! event's actual time. Other clocks survive unchanged. This is the continuous-
//! time first-reaction construction, not a fixed-step Bernoulli approximation.
//!
//! Clocks and the original `fs_rand::Stream` checkpoint belong to the inventory.
//! Observation endpoints never redraw a clock. Compensated event times retain
//! positive waits smaller than one ULP of the absolute observation clock. Clone
//! and accepted interval partitioning therefore preserve the same trajectory,
//! including cumulative energy bits, on the same ISA/toolchain profile.
//! Independent realizations need distinct caller-supplied stream keys; cloning
//! a trajectory or reusing its key deliberately reproduces its randomness.
//!
//! This lane models demographic fluctuations, not uncertainty in half-lives or
//! per-event radiation spectra. It shares the source model's mean branch energy
//! and one represented daughter per event. Probabilities inherit the RNG's
//! 53-bit uniform grid and floating-point clock arithmetic; no exact probability
//! law over arbitrary binary64 extremes or cross-ISA trajectory proof is claimed.

use std::sync::Arc;

use fs_rand::Stream;
pub use fs_rand::{StreamCheckpoint, StreamKey};

use super::{
    Accumulator, DecayError, DecayModel, add, channel_array, checked, defect_bound,
    export_channels, nonnegative,
};
use crate::{EnergyChannels, channel_values};

/// Largest admitted total integer population, exactly representable in binary64.
/// Larger macroscopic inventories belong to the deterministic expectation owner.
pub const MAX_STOCHASTIC_NUCLEI: u64 = 1_u64 << 53;

/// Explicit event-work and emitted-energy accounting allowances for one trial.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StochasticDecayBudget {
    /// Maximum transformations, including all newly born daughters' decays.
    /// Zero permits an interval containing no events. At most 2^53 is admitted.
    pub maximum_events: u64,
    /// Maximum accumulated noncancelling energy-account rounding defect [J].
    /// This is neither nuclear-data uncertainty nor an error in sampled counts.
    pub energy_tolerance_j: f64,
}

impl StochasticDecayBudget {
    fn validate(self) -> Result<(), DecayError> {
        if self.maximum_events > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget("stochastic event-count ceiling"));
        }
        if !nonnegative(self.energy_tolerance_j) {
            return Err(DecayError::Invalid(
                "stochastic energy tolerance must be finite and nonnegative",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum EventClock {
    Inactive,
    Finite(Accumulator),
    // A positive waiting time beyond every finite binary64 observation endpoint.
    // This is not stability: the original law and population remain unchanged.
    BeyondFiniteTime,
}

/// One accepted trajectory interval and its integer transformations.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayStep {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub duration_s: f64,
    /// Total transformations; exactly representable as f64 under the work limit.
    pub events: u64,
    /// Per-nuclide transformations in the model's original order.
    pub decays: Vec<u64>,
    /// Per-branch transformations in the model's original order.
    pub branch_decays: Vec<u64>,
    pub released: EnergyChannels,
    pub energy_roundoff_j: f64,
    /// Includes resampling clocks for future events, not only completed events.
    pub rng_draws: u64,
}

/// Closed finite-population trajectory, with pending events and RNG owned together.
///
/// Construction bounds total nuclei to 2^53, admits complete branch closure,
/// and prepares the first independent clocks. Empty/stable populations consume
/// no randomness. Clone retains the complete trajectory; an RNG checkpoint alone
/// is not a substitute for this physical state when resuming a simulation.
#[derive(Clone, Debug, PartialEq)]
pub struct StochasticDecayInventory {
    model: Arc<DecayModel>,
    outgoing: Arc<Vec<Vec<usize>>>,
    populations: Vec<u64>,
    decays: Vec<u64>,
    branch_decays: Vec<u64>,
    released: [Accumulator; 3],
    time_s: f64,
    clocks: Vec<EventClock>,
    rng: StreamCheckpoint,
}

impl StochasticDecayInventory {
    /// Bind an explicit integer population and logical random-stream identity at t=0.
    pub fn new(
        model: Arc<DecayModel>,
        populations: Vec<u64>,
        key: StreamKey,
    ) -> Result<Self, DecayError> {
        Self::new_at(model, populations, key, 0.0)
    }

    /// Bind supplied populations at an explicit current time. No preceding decay
    /// history, equilibrium or extra random-stream identity is inferred.
    pub fn new_at(
        model: Arc<DecayModel>,
        populations: Vec<u64>,
        key: StreamKey,
        time_s: f64,
    ) -> Result<Self, DecayError> {
        if !nonnegative(time_s) || populations.len() != model.nuclides.len() {
            return Err(DecayError::Invalid(
                "one integer population per nuclide and a finite nonnegative time required",
            ));
        }
        let total = populations
            .iter()
            .try_fold(0_u64, |sum, &count| sum.checked_add(count))
            .ok_or(DecayError::Budget("stochastic population ceiling"))?;
        if total > MAX_STOCHASTIC_NUCLEI {
            return Err(DecayError::Budget("stochastic population ceiling"));
        }
        let nodes = populations.len();
        let mut outgoing = vec![Vec::new(); nodes];
        let mut probabilities = vec![Accumulator::new(0.0); nodes];
        for (index, branch) in model.branches.iter().enumerate() {
            outgoing[branch.parent].push(index);
            add(&mut probabilities[branch.parent], branch.probability)?;
        }
        for (index, probability) in probabilities.iter().enumerate() {
            // A deterministic model admits a small closure residual against its
            // numerical accounting budget. An integer event must actually select
            // one daughter: do not normalize an incomplete authored distribution.
            if model.rates[index] > 0.0 && probability.value() != 1.0 {
                return Err(DecayError::Invalid(
                    "stochastic branch probabilities must sum to one without normalization",
                ));
            }
        }
        let branches = model.branches.len();
        let mut result = Self {
            model,
            outgoing: Arc::new(outgoing),
            populations,
            decays: vec![0; nodes],
            branch_decays: vec![0; branches],
            released: [Accumulator::new(0.0); 3],
            time_s,
            clocks: vec![EventClock::Inactive; nodes],
            rng: StreamCheckpoint::current(key, 0),
        };
        for node in 0..nodes {
            result.reschedule(node, Accumulator::new(time_s))?;
        }
        Ok(result)
    }

    #[must_use]
    pub fn model(&self) -> &Arc<DecayModel> {
        &self.model
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
    pub fn cumulative_decays(&self) -> &[u64] {
        &self.decays
    }
    #[must_use]
    pub fn cumulative_branch_decays(&self) -> &[u64] {
        &self.branch_decays
    }
    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        channel_values(self.released)
    }
    /// Original versioned stream position, including all pending-event draws.
    #[must_use]
    pub fn rng_checkpoint(&self) -> StreamCheckpoint {
        self.rng
    }

    /// Earliest representable observation endpoint containing the next event.
    /// A positive compensated low part rounds UP, so a thermal owner clips to
    /// an endpoint that actually contains the event. Several sub-ULP events may
    /// share this endpoint; `advance_to` processes them in actual event order.
    /// None means no event at any finite representable future endpoint.
    #[must_use]
    pub fn next_event_time_s(&self) -> Option<f64> {
        self.next_event().map(|(_, time)| ceiling(time))
    }

    /// Evaluate a complete candidate at an absolute, nondecreasing endpoint.
    /// Zero-duration observation consumes no draws. Events at the endpoint are
    /// committed once, together with their branch energy and future clocks.
    /// Work, energy and callback refusal retain the entire original trajectory.
    pub fn trial_to(
        &self,
        end_time_s: f64,
        budget: StochasticDecayBudget,
        mut checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<(Self, StochasticDecayStep), DecayError> {
        checkpoint()?;
        budget.validate()?;
        if !nonnegative(end_time_s) || end_time_s < self.time_s {
            return Err(DecayError::Invalid(
                "stochastic decay endpoint must be finite and not precede the current time",
            ));
        }
        let mut next = self.clone();
        let mut report = StochasticDecayStep {
            start_time_s: self.time_s,
            end_time_s,
            duration_s: end_time_s - self.time_s,
            events: 0,
            decays: vec![0; self.populations.len()],
            branch_decays: vec![0; self.model.branches.len()],
            released: EnergyChannels::default(),
            energy_roundoff_j: 0.0,
            rng_draws: 0,
        };
        let mut released = [Accumulator::new(0.0); 3];
        let mut energy_error = Accumulator::new(0.0);
        while let Some((parent, time)) = next.next_event() {
            if ceiling(time) > end_time_s || end_time_s == self.time_s {
                break;
            }
            checkpoint()?;
            if report.events == budget.maximum_events {
                return Err(DecayError::Budget("stochastic decay events"));
            }
            let branch_index = next.choose_branch(parent)?;
            let branch = next.model.branches[branch_index];
            next.populations[parent] = next.populations[parent]
                .checked_sub(1)
                .ok_or(DecayError::Arithmetic)?;
            increment(&mut next.populations[branch.daughter])?;
            increment(&mut next.decays[parent])?;
            increment(&mut next.branch_decays[branch_index])?;
            increment(&mut report.events)?;
            increment(&mut report.decays[parent])?;
            increment(&mut report.branch_decays[branch_index])?;
            for ((interval, lifetime), value) in released
                .iter_mut()
                .zip(&mut next.released)
                .zip(channel_array(branch.energy_per_decay))
            {
                add(&mut energy_error, add(interval, value)?.abs())?;
                // Update the shared compensated energy account per actual event,
                // never by a partition-dependent sum of interval reports.
                add(&mut energy_error, add(lifetime, value)?.abs())?;
            }
            next.reschedule(parent, time)?;
            next.reschedule(branch.daughter, time)?;
        }
        report.released = export_channels(released, &mut energy_error)?;
        checked(channel_values(next.released).total_j())?;
        report.energy_roundoff_j = defect_bound(energy_error)?;
        if report.energy_roundoff_j > budget.energy_tolerance_j {
            return Err(DecayError::Budget("stochastic emission energy roundoff"));
        }
        report.rng_draws = next
            .rng
            .index
            .checked_sub(self.rng.index)
            .ok_or(DecayError::Arithmetic)?;
        next.time_s = end_time_s;
        checkpoint()?;
        Ok((next, report))
    }

    /// Atomically accept a fully admitted trajectory interval.
    pub fn advance_to(
        &mut self,
        end_time_s: f64,
        budget: StochasticDecayBudget,
        checkpoint: impl FnMut() -> Result<(), DecayError>,
    ) -> Result<StochasticDecayStep, DecayError> {
        let (next, report) = self.trial_to(end_time_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }

    fn next_event(&self) -> Option<(usize, Accumulator)> {
        let mut next: Option<(usize, Accumulator)> = None;
        for (node, clock) in self.clocks.iter().enumerate() {
            if let EventClock::Finite(time) = clock
                && next.is_none_or(|(_, earliest)| before(*time, earliest))
            {
                // Equal clocks retain original node order as the tie-break.
                next = Some((node, *time));
            }
        }
        next
    }

    fn reschedule(&mut self, node: usize, origin: Accumulator) -> Result<(), DecayError> {
        let rate = self.model.rates[node];
        let nuclei = self.populations[node];
        if nuclei == 0 || rate == 0.0 {
            self.clocks[node] = EventClock::Inactive;
            return Ok(());
        }
        let mut stream = self.stream_for_draw()?;
        let unit_wait = stream.next_exponential();
        self.rng = stream.checkpoint();
        if !(unit_wait > 0.0) || !unit_wait.is_finite() {
            // The RNG's finite grid includes zero. Do not invent a positive
            // lifetime, or permit a nonprogressing event, for that boundary.
            return Err(DecayError::Arithmetic);
        }
        // Count is exact and >=1. This order avoids overflowing lambda*N or
        // unit_wait/lambda when the final physical wait is representable.
        let wait = (unit_wait / nuclei as f64) / rate;
        if wait.is_infinite() || (origin.hi + wait).is_infinite() {
            self.clocks[node] = EventClock::BeyondFiniteTime;
            return Ok(());
        }
        if wait == 0.0 {
            return Err(DecayError::Arithmetic);
        }
        let mut time = origin;
        add(&mut time, wait)?;
        if !before(origin, time) {
            return Err(DecayError::Arithmetic);
        }
        self.clocks[node] = if ceiling(time).is_finite() {
            EventClock::Finite(time)
        } else {
            EventClock::BeyondFiniteTime
        };
        Ok(())
    }

    fn choose_branch(&mut self, parent: usize) -> Result<usize, DecayError> {
        let choices = &self.outgoing[parent];
        if choices.len() == 1 {
            return Ok(choices[0]);
        }
        let mut stream = self.stream_for_draw()?;
        let draw = stream.next_f64();
        self.rng = stream.checkpoint();
        let mut cumulative = Accumulator::new(0.0);
        for &index in choices {
            add(&mut cumulative, self.model.branches[index].probability)?;
            if draw < cumulative.value() {
                return Ok(index);
            }
        }
        // Complete source probabilities were admitted at construction. Never
        // repair an uncovered draw by assigning it to an arbitrary last branch.
        Err(DecayError::Arithmetic)
    }

    fn stream_for_draw(&self) -> Result<Stream, DecayError> {
        // fs-rand intentionally wraps modulo 2^64. This physical owner refuses
        // before wrap, so no accepted trajectory silently reuses its first draw.
        if self.rng.index == u64::MAX {
            return Err(DecayError::Budget("stochastic random stream exhausted"));
        }
        Stream::resume(self.rng)
            .map_err(|_| DecayError::Invalid("unsupported stochastic random checkpoint"))
    }
}

fn increment(value: &mut u64) -> Result<(), DecayError> {
    *value = value.checked_add(1).ok_or(DecayError::Arithmetic)?;
    Ok(())
}

fn before(a: Accumulator, b: Accumulator) -> bool {
    a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo)
}

fn ceiling(time: Accumulator) -> f64 {
    if time.lo > 0.0 {
        time.hi.next_up()
    } else {
        time.hi
    }
}

#[cfg(test)]
mod tests;

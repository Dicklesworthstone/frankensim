//! Bounded DE/rand/1/bin with frozen generations and resumable physical callbacks.
//!
//! Storn and Price (1997), doi:10.1023/A:1008202821328: three distinct donors,
//! excluding the target; one weighted difference; binomial crossover with one
//! forced donor coordinate; greedy selection. Mutation uses unit-box coordinates.
//! Out-of-box coordinates are independently redrawn uniformly, not clipped or
//! passed to the objective. Equal scores accept the trial, so a flat region
//! does not freeze every population member in place.
//!
//! Every generation reads the previous complete population. Successful trials
//! publish together. Completed callbacks survive cancellation; failing callbacks
//! spend an attempt, retain the candidate for retry, and never become invented
//! penalties. Keep the objective and its non-design physical inputs fixed.
//!
//! Philox draws are keyed by seed, candidate and generation. Pausing, retrying or
//! extending budgets does not redraw candidates. No flat-objective convergence
//! rule is used: a flat phase fraction can hide unexplored heating regimes.
//! Stops mean target reached, work exhausted or paused, not KKT stationarity,
//! guaranteed coverage or global optimality. The serial path and optional
//! tile-backed parallel waves share candidate generation and selection.

/// Parallel objective waves on the existing cancellation-aware tile executor.
pub mod parallel;
/// Portable byte transport for quiescent search states, including cached waves.
pub mod snapshot;
/// Feasibility-first selection without a weighted penalty objective.
pub mod constrained;

use fs_exec::Cx;
use fs_rand::{Stream, StreamKey};
use std::collections::VecDeque;
use std::fmt;

/// Kernel identity for the candidate streams of this algorithm.
pub const DIFFERENTIAL_EVOLUTION_STREAM_KERNEL: u32 = 0xDE01;

/// Immutable numerical policy, apart from explicit lifetime-budget extension.
#[derive(Debug, Clone, PartialEq)]
pub struct DifferentialEvolutionConfig {
    /// At least four members; a count, not a multiplier of dimension.
    pub population_size: usize,
    /// Fixed F in (0,2].
    pub differential_weight: f64,
    /// CR in [0,1]; a donor coordinate is forced even at CR=0.
    pub crossover_probability: f64,
    pub seed: u64,
    /// Lifetime callback attempts, including failures. May be zero.
    pub max_evaluations: usize,
    /// Population-size times dimension. Two populations, at most one parallel
    /// wave of candidate coordinates, and O(dimension) candidate/best buffers
    /// coexist; this is a shape cap, not a total allocation or byte bound.
    pub max_population_entries: usize,
    /// Optional finite objective threshold, not an inferred global minimum.
    pub target: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DifferentialEvolutionError {
    Invalid(&'static str),
    NonFinite(&'static str),
    Allocation,
    /// A candidate's bounded random-draw segment was exhausted.
    RandomBudget,
    Cancelled,
}
impl fmt::Display for DifferentialEvolutionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) | Self::NonFinite(s) => f.write_str(s),
            Self::Allocation => f.write_str("differential-evolution allocation refused"),
            Self::RandomBudget => f.write_str("differential-evolution random draw allowance exhausted"),
            Self::Cancelled => f.write_str("differential evolution cancelled; completed work retained"),
        }
    }
}
impl std::error::Error for DifferentialEvolutionError {}

#[derive(Debug)]
pub enum DifferentialEvolutionRunError<E> {
    Solver(DifferentialEvolutionError),
    /// Original callback error; no artificial fitness replacement.
    Evaluation(E),
}
impl<E> From<DifferentialEvolutionError> for DifferentialEvolutionRunError<E> {
    fn from(error: DifferentialEvolutionError) -> Self { Self::Solver(error) }
}
impl<E: fmt::Display> fmt::Display for DifferentialEvolutionRunError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Solver(e) => write!(f, "{e}"),
            Self::Evaluation(e) => write!(f, "differential-evolution objective: {e}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DifferentialEvolutionRunError<E> {}

/// A finite observation in the coordinates actually supplied to the callback.
#[derive(Debug, Clone, PartialEq)]
pub struct DifferentialEvolutionPoint<S = f64> {
    pub coordinates: Vec<f64>,
    pub value: S,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DifferentialEvolutionStop { TargetReached, EvaluationLimit, Paused }
#[derive(Debug, Clone, PartialEq)]
pub struct DifferentialEvolutionReport<S = f64> {
    pub stop: DifferentialEvolutionStop,
    /// Best consumed observation, including unfinished generations.
    /// Exact equal scores preserve the earliest best observation.
    pub best: Option<DifferentialEvolutionPoint<S>>,
    pub evaluations: usize,
    /// Complete evolutionary generations, excluding initialization.
    pub generations: usize,
    pub initialized_members: usize,
    pub pending_members: usize,
}

#[derive(Debug)]
struct Pending<S = f64> { coordinates: Vec<f64>, value: Option<S> }

// Keep the public scalar API and its snapshot format source-compatible. Both
// score types run ONE candidate/cache/budget/publication implementation.
/// Original scalar DE search. Its numerical policy and checkpoint bytes are unchanged.
pub type DifferentialEvolutionState = RankedDifferentialEvolutionState<f64>;

mod score_sealed {
    pub trait Sealed {}
    impl Sealed for f64 {}
    impl Sealed for super::constrained::FeasibilityScore {}
}

/// Internal, sealed score algebra for the shared serial/parallel DE worker.
/// Use DifferentialEvolutionState or ConstrainedDifferentialEvolutionState;
/// user-defined orderings cannot bypass the finite-score or target rules.
#[doc(hidden)]
pub trait DifferentialEvolutionScore: score_sealed::Sealed + Copy + fmt::Debug + PartialEq + Send + Sync {
    fn checked(self) -> Result<Self, DifferentialEvolutionError>;
    fn rank(self, other: Self) -> std::cmp::Ordering;
    fn accepts(self, parent: Self) -> bool;
    fn reaches(self, target: f64) -> bool;
}
impl DifferentialEvolutionScore for f64 {
    fn checked(self) -> Result<Self, DifferentialEvolutionError> { finite(self) }
    fn rank(self, other: Self) -> std::cmp::Ordering { self.total_cmp(&other) }
    // Preserve the original scalar replacement semantics, including signed zero.
    fn accepts(self, parent: Self) -> bool { self <= parent }
    fn reaches(self, target: f64) -> bool { self <= target }
}

/// Reports are not search state. Snapshots retain the work spent at capture;
/// replaying an older snapshot cannot account for work done after that capture.
#[doc(hidden)]
#[derive(Debug)]
pub struct RankedDifferentialEvolutionState<S: DifferentialEvolutionScore> {
    bounds: Vec<[f64; 2]>,
    initial: Option<Vec<f64>>,
    config: DifferentialEvolutionConfig,
    population: Vec<DifferentialEvolutionPoint<S>>,
    next: Vec<DifferentialEvolutionPoint<S>>,
    pending: Option<Pending<S>>,
    // Future members of the SAME frozen generation. Successful parallel scores
    // survive cancellation; the serial path can also drain/retry this queue.
    queued: VecDeque<Pending<S>>,
    best: Option<DifferentialEvolutionPoint<S>>,
    evaluations: usize,
    generations: usize,
    random_stride: u64,
}

impl<S: DifferentialEvolutionScore> RankedDifferentialEvolutionState<S> {
    /// Admit strictly ordered finite bounds without invoking a callback. A
    /// supplied initial point is evaluated exactly as supplied as member zero;
    /// remaining members sample the whole box. Remove fixed variables upstream.
    /// Each interval must have a finite representable width.
    pub fn new(cx: &Cx<'_>, bounds: &[[f64; 2]], initial: Option<&[f64]>,
        config: DifferentialEvolutionConfig) -> Result<Self, DifferentialEvolutionError>
    { Self::new_with(bounds, initial, config, &mut || cx.checkpoint().is_ok()) }

    fn new_with(bounds: &[[f64; 2]], initial: Option<&[f64]>,
        config: DifferentialEvolutionConfig, poll: &mut impl FnMut() -> bool)
        -> Result<Self, DifferentialEvolutionError>
    {
        checkpoint(poll)?;
        let n = bounds.len();
        if n == 0 || config.population_size < 4 || config.population_size > u32::MAX as usize
            || n.checked_mul(config.population_size).is_none_or(|v| v > config.max_population_entries)
            || !config.differential_weight.is_finite() || config.differential_weight <= 0.0
            || config.differential_weight > 2.0 || !config.crossover_probability.is_finite()
            || !(0.0..=1.0).contains(&config.crossover_probability)
            || config.target.is_some_and(|v| !v.is_finite())
            || initial.is_some_and(|p| p.len() != n)
        { return Err(invalid("DE requires a bounded nonempty shape, population >=4, F in (0,2], CR in [0,1] and finite target")); }
        for (i, &[lo, hi]) in bounds.iter().enumerate() {
            if i % 256 == 0 { checkpoint(poll)?; }
            if !lo.is_finite() || !hi.is_finite() || lo >= hi || !(hi-lo).is_finite()
                || initial.is_some_and(|p| !p[i].is_finite() || p[i] < lo || p[i] > hi)
            { return Err(invalid("DE needs finite ordered bounds with representable widths containing its initial point")); }
        }
        // Four bounded integer draws, one crossover and at most one repair per coordinate.
        let random_stride = u64::try_from(n).ok().and_then(|n| n.checked_mul(2))
            .and_then(|n| n.checked_add(4*64)).ok_or(DifferentialEvolutionError::RandomBudget)?;
        let mut stored_bounds = reserved(n)?;
        for (i, &bound) in bounds.iter().enumerate() {
            if i % 256 == 0 { checkpoint(poll)?; }
            stored_bounds.push(bound);
        }
        let initial = initial.map(|p| copied(p, poll)).transpose()?;
        let population = reserved(config.population_size)?;
        let next = reserved(config.population_size)?;
        checkpoint(poll)?;
        Ok(Self { bounds: stored_bounds, initial, config, population, next,
            pending: None, queued: VecDeque::new(), best: None,
            evaluations: 0, generations: 0, random_stride })
    }

    pub const fn evaluations(&self) -> usize { self.evaluations }
    pub const fn evaluation_limit(&self) -> usize { self.config.max_evaluations }
    pub const fn generations(&self) -> usize { self.generations }
    pub fn best(&self) -> Option<&DifferentialEvolutionPoint<S>> { self.best.as_ref() }
    /// Completed initialization members or the last complete generation.
    /// Partial next-generation results never alter this slice.
    pub fn population(&self) -> &[DifferentialEvolutionPoint<S>] { &self.population }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), DifferentialEvolutionError> {
        if total <= self.config.max_evaluations { return Err(invalid("DE evaluation allowance must increase")); }
        self.config.max_evaluations = total;
        Ok(())
    }

    /// Respect per-call AND lifetime callback ceilings. A completed finite
    /// callback is cached before checking cancellation; retries never repeat
    /// it. A failing callback retains its candidate, but retries cost another
    /// attempt. Only complete generations replace parents. Callback side
    /// effects or panics cannot be rolled back.
    pub fn advance<E>(&mut self, cx: &Cx<'_>, additional_evaluations: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<S, E>)
        -> Result<DifferentialEvolutionReport<S>, DifferentialEvolutionRunError<E>>
    { self.advance_with(additional_evaluations, evaluate, &mut || cx.checkpoint().is_ok()) }

    fn advance_with<E>(&mut self, additional: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<S, E>, poll: &mut impl FnMut() -> bool)
        -> Result<DifferentialEvolutionReport<S>, DifferentialEvolutionRunError<E>>
    {
        let mut spent = 0;
        loop {
            checkpoint(poll)?;
            if self.pending.is_none() { self.pending = self.queued.pop_front(); }
            if self.pending.as_ref().is_some_and(|p| p.value.is_some()) {
                self.consume(poll)?;
                continue;
            }
            let stop = if self.best.as_ref().is_some_and(|p| self.config.target.is_some_and(|t| p.value.reaches(t))) {
                Some(DifferentialEvolutionStop::TargetReached)
            } else if self.evaluations >= self.config.max_evaluations {
                Some(DifferentialEvolutionStop::EvaluationLimit)
            } else if spent >= additional { Some(DifferentialEvolutionStop::Paused) } else { None };
            if let Some(stop) = stop { return Ok(self.report(stop, poll)?); }
            if self.pending.is_none() {
                let coordinates = self.candidate(poll)?;
                self.pending = Some(Pending { coordinates, value: None });
            }
            checkpoint(poll)?;
            self.evaluations += 1;
            spent += 1;
            let value = evaluate(&self.pending.as_ref().expect("prepared DE candidate").coordinates)
                .map_err(DifferentialEvolutionRunError::Evaluation)?;
            let value = value.checked()?;
            self.pending.as_mut().expect("retained DE candidate").value = Some(value);
        }
    }

    fn candidate(&self, poll: &mut impl FnMut() -> bool) -> Result<Vec<f64>, DifferentialEvolutionError> {
        let initializing = self.population.len() < self.config.population_size;
        let member = if initializing { self.population.len() } else { self.next.len() };
        self.candidate_at(member, poll)
    }

    // Random identity is the logical member, never its wave/worker position.
    fn candidate_at(&self, member: usize, poll: &mut impl FnMut() -> bool)
        -> Result<Vec<f64>, DifferentialEvolutionError>
    {
        let initializing = self.population.len() < self.config.population_size;
        if initializing && member == 0 {
            if let Some(initial) = &self.initial { return copied(initial, poll); }
        }
        let generation = if initializing { 0 } else {
            self.generations.checked_add(1).ok_or(DifferentialEvolutionError::RandomBudget)?
        };
        let base = u64::try_from(generation).ok().and_then(|g| g.checked_mul(self.random_stride))
            .ok_or(DifferentialEvolutionError::RandomBudget)?;
        let end = base.checked_add(self.random_stride).ok_or(DifferentialEvolutionError::RandomBudget)?;
        let mut rng = CandidateRandom { key: StreamKey { seed: self.config.seed,
            kernel: DIFFERENTIAL_EVOLUTION_STREAM_KERNEL, tile: member as u32 }, next: base, end };
        let mut coordinates = reserved(self.bounds.len())?;
        if initializing {
            for (i, &[lo, hi]) in self.bounds.iter().enumerate() {
                if i % 256 == 0 { checkpoint(poll)?; }
                coordinates.push(interpolate(lo, hi, rng.uniform()?)?);
            }
        } else {
            let donors = donors(&mut rng, self.config.population_size, member, poll)?;
            let forced = rng.below(self.bounds.len(), poll)?;
            for (i, &[lo, hi]) in self.bounds.iter().enumerate() {
                if i % 256 == 0 { checkpoint(poll)?; }
                let cross = rng.uniform()? < self.config.crossover_probability;
                if cross || i == forced {
                    let unit = |k: usize| (self.population[k].coordinates[i]-lo)/(hi-lo);
                    let mut value = self.config.differential_weight.mul_add(
                        unit(donors[1])-unit(donors[2]), unit(donors[0]));
                    if !(0.0..=1.0).contains(&value) { value = rng.uniform()?; }
                    coordinates.push(interpolate(lo, hi, value)?);
                } else { coordinates.push(self.population[member].coordinates[i]); }
            }
        }
        checkpoint(poll)?;
        Ok(coordinates)
    }

    fn consume(&mut self, poll: &mut impl FnMut() -> bool) -> Result<(), DifferentialEvolutionError> {
        let pending = self.pending.as_ref().expect("completed DE candidate");
        let value = pending.value.expect("finite DE score");
        let initializing = self.population.len() < self.config.population_size;
        let next_generation = if !initializing && self.next.len()+1 == self.config.population_size {
            Some(self.generations.checked_add(1).ok_or_else(|| invalid("DE generation count overflow"))?)
        } else { None };
        let replacement = if !initializing && !value.accepts(self.population[self.next.len()].value) {
            Some(copy_point(&self.population[self.next.len()], poll)?)
        } else { None };
        let best = if self.best.as_ref().is_none_or(|p| value.rank(p.value).is_lt()) {
            Some(DifferentialEvolutionPoint { coordinates: copied(&pending.coordinates, poll)?, value })
        } else { None };
        checkpoint(poll)?;
        // All fallible work precedes publication. Both buffers were reserved.
        let pending = self.pending.take().expect("completed DE candidate retained");
        let point = replacement.unwrap_or(DifferentialEvolutionPoint { coordinates: pending.coordinates, value });
        if let Some(best) = best { self.best = Some(best); }
        if initializing { self.population.push(point); } else { self.next.push(point); }
        if let Some(generation) = next_generation {
            std::mem::swap(&mut self.population, &mut self.next);
            self.next.clear();
            self.generations = generation;
        }
        Ok(())
    }

    fn report(&self, stop: DifferentialEvolutionStop, poll: &mut impl FnMut() -> bool)
        -> Result<DifferentialEvolutionReport<S>, DifferentialEvolutionError>
    {
        Ok(DifferentialEvolutionReport { stop,
            best: self.best.as_ref().map(|p| copy_point(p, poll)).transpose()?,
            evaluations: self.evaluations, generations: self.generations,
            initialized_members: self.population.len(), pending_members: self.next.len() })
    }
}

// Random access avoids RNG-state mutation on a cancelled candidate build.
// Bounded Lemire reduction uses original fs-rand Philox words. Exhausting the
// explicit 64-draw rejection allowance refuses instead of introducing bias.
struct CandidateRandom { key: StreamKey, next: u64, end: u64 }
impl CandidateRandom {
    fn bits(&mut self) -> Result<u64, DifferentialEvolutionError> {
        if self.next >= self.end { return Err(DifferentialEvolutionError::RandomBudget); }
        let words = Stream::at(self.key, self.next);
        self.next += 1;
        Ok(u64::from(words[0]) | (u64::from(words[1]) << 32))
    }
    fn uniform(&mut self) -> Result<f64, DifferentialEvolutionError> {
        Ok((self.bits()? >> 11) as f64 * (1.0/9_007_199_254_740_992.0))
    }
    fn below(&mut self, n: usize, poll: &mut impl FnMut() -> bool) -> Result<usize, DifferentialEvolutionError> {
        let n = u64::try_from(n).ok().filter(|&v| v > 0).ok_or(DifferentialEvolutionError::RandomBudget)?;
        let threshold = n.wrapping_neg() % n;
        for i in 0..64 {
            if i % 8 == 0 { checkpoint(poll)?; }
            let product = u128::from(self.bits()?) * u128::from(n);
            if product as u64 >= threshold { return Ok((product >> 64) as usize); }
            // rejected: next bounded draw retains the original candidate identity
        }
        Err(DifferentialEvolutionError::RandomBudget)
    }
}
fn donors(rng: &mut CandidateRandom, count: usize, target: usize, poll: &mut impl FnMut() -> bool)
    -> Result<[usize; 3], DifferentialEvolutionError>
{
    let mut excluded = [target, 0, 0, 0];
    let mut selected = [0; 3];
    for i in 0..3 {
        excluded[..i+1].sort_unstable();
        let mut rank = rng.below(count-i-1, poll)?;
        for &index in &excluded[..i+1] { if rank >= index { rank += 1; } }
        selected[i] = rank;
        excluded[i+1] = rank;
    }
    Ok(selected)
}
fn interpolate(lo: f64, hi: f64, unit: f64) -> Result<f64, DifferentialEvolutionError> {
    // Last-bit interpolation rounding may reach an endpoint but not leave the box.
    Ok(finite(unit.mul_add(hi-lo, lo))?.clamp(lo, hi))
}
fn copy_point<S: DifferentialEvolutionScore>(point: &DifferentialEvolutionPoint<S>, poll: &mut impl FnMut() -> bool)
    -> Result<DifferentialEvolutionPoint<S>, DifferentialEvolutionError>
{ Ok(DifferentialEvolutionPoint { coordinates: copied(&point.coordinates, poll)?, value: point.value }) }
fn copied(input: &[f64], poll: &mut impl FnMut() -> bool) -> Result<Vec<f64>, DifferentialEvolutionError> {
    let mut out = reserved(input.len())?;
    for (i, &value) in input.iter().enumerate() {
        if i % 256 == 0 { checkpoint(poll)?; }
        out.push(value);
    }
    Ok(out)
}
fn reserved<T>(n: usize) -> Result<Vec<T>, DifferentialEvolutionError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| DifferentialEvolutionError::Allocation)?;
    Ok(out)
}
fn invalid(what: &'static str) -> DifferentialEvolutionError { DifferentialEvolutionError::Invalid(what) }
fn finite(value: f64) -> Result<f64, DifferentialEvolutionError> {
    if value.is_finite() { Ok(value) } else {
        Err(DifferentialEvolutionError::NonFinite("DE requires finite arithmetic and objective values"))
    }
}
fn checkpoint(poll: &mut impl FnMut() -> bool) -> Result<(), DifferentialEvolutionError> {
    if poll() { Ok(()) } else { Err(DifferentialEvolutionError::Cancelled) }
}

#[cfg(test)]
mod tests;

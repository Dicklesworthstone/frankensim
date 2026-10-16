//! Fine solver clocks inside coarse reconstruction windows.
//!
//! A window's knot states and model-error penalties are unchanged when its
//! forecast is refined. Inner `record` receives the GLOBAL fine-step index, so
//! prescribed source/boundary schedules and their derivatives use the same
//! address in forward and reverse. This composes existing maps, not a new solver.
//!
//! Forward storage is two states plus one small endpoint witness per substep.
//! Binary replay retains O(log(substeps)) parked states and ONE inner tape at a
//! time. Endpoint hashes diagnose changed replay; they do not certify a model,
//! derivative, or physical history. All data/solvers must remain immutable.

use super::{IntervalScheme, IntervalTape};
use super::super::{WindowError, copy, poll, zeros};
use fs_blake3::{Blake3, ContentHash};
use fs_time::adaptive::adjoint::trajectory::TrajectoryGradient;

/// Explicit fine clock and the indices of its reconstruction knots. Knots must
/// include the first/last fine endpoints. No nearest-time lookup is performed.
#[derive(Debug, Clone)]
pub struct SubstepGrid {
    fine: Vec<f64>,
    indices: Vec<usize>,
    knots: Vec<f64>,
}
impl SubstepGrid {
    pub fn new(fine: &[f64], indices: &[usize], max_steps: usize) -> Result<Self, WindowError> {
        if fine.len() < 2 || fine.len()-1 > max_steps || indices.len() < 2 || indices.len() > fine.len()
            || indices.first() != Some(&0) || indices.last() != Some(&(fine.len()-1))
            || indices.windows(2).any(|p| p[0] >= p[1])
            || fine.iter().any(|v| !v.is_finite())
            || fine.windows(2).any(|p| p[0] >= p[1] || !(p[1]-p[0]).is_finite())
        { return Err(WindowError::Invalid("ordered finite substep clock with bounded, covering knot indices required")); }
        // Increasing indices ending at fine.len()-1 are all in bounds.
        let mut owned = Vec::new();
        owned.try_reserve_exact(indices.len()).map_err(|_| WindowError::Allocation)?;
        owned.extend_from_slice(indices);
        let mut knots = zeros(indices.len())?;
        for (value, &index) in knots.iter_mut().zip(indices) { *value = fine[index]; }
        if knots.windows(2).any(|p| !(p[1]-p[0]).is_finite()) {
            return Err(WindowError::Invalid("substep knot duration overflow"));
        }
        Ok(Self { fine: copy(fine)?, indices: owned, knots })
    }

    /// Equal fractions WITHIN each knot interval, preserving its exact endpoints.
    /// The actual adjacent binary64 differences are the inner step durations;
    /// equal representable spacing is not assumed. An unresolvable cut refuses.
    /// One subdivision returns the original clock bits unchanged.
    pub fn uniform(knots: &[f64], subdivisions: usize, max_steps: usize) -> Result<Self, WindowError> {
        if knots.len() < 2 || subdivisions == 0 || knots.iter().any(|v| !v.is_finite())
            || knots.windows(2).any(|p| p[0] >= p[1] || !(p[1]-p[0]).is_finite())
        { return Err(WindowError::Invalid("positive subdivision count and ordered finite knots required")); }
        let steps = (knots.len()-1).checked_mul(subdivisions)
            .filter(|v| *v <= max_steps && *v < usize::MAX)
            .ok_or(WindowError::Invalid("substep clock exceeds step allowance"))?;
        let mut fine = zeros(steps+1)?;
        let mut indices = Vec::new();
        indices.try_reserve_exact(knots.len()).map_err(|_| WindowError::Allocation)?;
        fine[0] = knots[0]; indices.push(0);
        for (k, interval) in knots.windows(2).enumerate() {
            let start = k*subdivisions;
            for j in 1..subdivisions {
                fine[start+j] = (interval[1]-interval[0]).mul_add(j as f64/subdivisions as f64, interval[0]);
            }
            fine[start+subdivisions] = interval[1];
            indices.push(start+subdivisions);
        }
        Self::new(&fine, &indices, max_steps)
    }
    pub fn fine_times(&self) -> &[f64] { &self.fine }
    pub fn knot_times(&self) -> &[f64] { &self.knots }
    pub fn knot_indices(&self) -> &[usize] { &self.indices }
}

/// Limits for one coarse interval. The inner policy separately bounds each
/// solve/recording and its replay. These are not whole-process memory limits.
#[derive(Debug, Clone, Copy)]
pub struct SubstepBudget {
    pub max_state_components: usize,
    pub max_parameters: usize,
    /// Parked full states, including the initial anchor. Binary replay needs
    /// floor(log2(substeps))+1. Current working vectors/inner tape are separate.
    pub checkpoints: usize,
    /// Inner `record` calls made by ONE reverse sweep, including its leaves.
    /// Inner-recording replay work is additional and bounded by its own policy.
    pub replayed_substeps: usize,
}

#[derive(Clone)]
pub struct CheckpointedIntervals<S> {
    inner: S,
    grid: SubstepGrid,
    budget: SubstepBudget,
}
impl<S> CheckpointedIntervals<S> {
    pub fn new(inner: S, grid: SubstepGrid, budget: SubstepBudget) -> Self { Self { inner, grid, budget } }
    pub fn inner(&self) -> &S { &self.inner }
    pub fn grid(&self) -> &SubstepGrid { &self.grid }
}

#[derive(Clone)]
struct Witness { endpoint: ContentHash, accepted_steps: usize }
/// Immutable completed coarse map; no inner linearizations are retained.
pub struct SubstepTape<'a, M, S> {
    policy: &'a CheckpointedIntervals<S>,
    model: &'a M,
    interval: usize,
    begin: usize,
    parameters: usize,
    initial: Vec<f64>,
    endpoint: Vec<f64>,
    witnesses: Vec<Witness>,
    accepted_steps: usize,
}
fn hash(end: f64, values: &[f64], cancelled: &mut dyn FnMut() -> bool) -> Result<ContentHash, WindowError> {
    let mut h = Blake3::new();
    h.update(b"fs-ascent/interval-substeps/endpoint/v1\0");
    h.update(&end.to_bits().to_le_bytes());
    for chunk in values.chunks(256) {
        poll(cancelled)?;
        for value in chunk { h.update(&value.to_bits().to_le_bytes()); }
    }
    poll(cancelled)?;
    Ok(h.finalize())
}
fn add_count(total: &mut usize, value: usize) -> Result<(), WindowError> {
    *total = total.checked_add(value).ok_or(WindowError::Invalid("substep work count overflow"))?;
    Ok(())
}
fn checked_endpoint(tape: &impl IntervalTape, n: usize, end: f64, interval: usize) -> Result<(), WindowError> {
    if tape.end_time() != end || tape.endpoint().len() != n || tape.endpoint().iter().any(|v| !v.is_finite()) {
        return Err(WindowError::IntervalOutput { interval, what: "substep endpoint time/state" });
    }
    Ok(())
}
impl<M, S: IntervalScheme<M>> IntervalScheme<M> for CheckpointedIntervals<S> {
    type Tape<'a> = SubstepTape<'a, M, S> where Self: 'a, M: 'a;
    fn dimension(&self, model: &M) -> usize { self.inner.dimension(model) }
    fn parameter_count(&self, model: &M) -> usize { self.inner.parameter_count(model) }
    fn validate(&self, times: &[f64]) -> Result<(), WindowError> {
        if times != self.grid.knots { return Err(WindowError::Invalid("window must use the substep grid's knots")); }
        self.inner.validate(&self.grid.fine)
    }
    fn record<'a>(&'a self, model: &'a M, interval: usize, start: f64, end: f64,
        initial: &[f64], cancelled: &mut dyn FnMut() -> bool) -> Result<Self::Tape<'a>, WindowError>
    {
        poll(cancelled)?;
        if interval >= self.grid.knots.len()-1 || start != self.grid.knots[interval] || end != self.grid.knots[interval+1] {
            return Err(WindowError::Invalid("substep interval clock mismatch"));
        }
        let n = self.inner.dimension(model); let p = self.inner.parameter_count(model);
        let begin = self.grid.indices[interval]; let count = self.grid.indices[interval+1]-begin;
        let required = (usize::BITS-count.leading_zeros()) as usize;
        if n == 0 || n > self.budget.max_state_components || p > self.budget.max_parameters
            || initial.len() != n || initial.iter().any(|v| !v.is_finite()) || required > self.budget.checkpoints
        { return Err(WindowError::Invalid("substep state/parameter/checkpoint allowance or initial field")); }
        self.inner.validate(&self.grid.fine)?;
        let mut witnesses = Vec::new();
        witnesses.try_reserve_exact(count).map_err(|_| WindowError::Allocation)?;
        let original = copy(initial)?; let mut state = copy(initial)?; let mut accepted_steps = 0;
        let mut stopped = false;
        let mut check = || { stopped |= cancelled(); stopped };
        for i in begin..begin+count {
            poll(&mut check)?;
            if self.inner.dimension(model) != n || self.inner.parameter_count(model) != p {
                return Err(WindowError::Invalid("substep model dimensions changed"));
            }
            let tape = self.inner.record(model, i, self.grid.fine[i], self.grid.fine[i+1], &state, &mut check)?;
            checked_endpoint(&tape, n, self.grid.fine[i+1], interval)?;
            let endpoint = hash(tape.end_time(), tape.endpoint(), &mut check)?;
            add_count(&mut accepted_steps, tape.accepted_steps())?;
            witnesses.push(Witness { endpoint, accepted_steps: tape.accepted_steps() });
            state = copy(tape.endpoint())?;
        }
        poll(&mut check)?;
        Ok(SubstepTape { policy: self, model, interval, begin, parameters: p, initial: original,
            endpoint: state, witnesses, accepted_steps })
    }
}
struct Progress { calls: usize, work: usize, peak: usize }
struct Bar { state: Vec<f64>, parameters: Vec<f64> }
impl<M, S: IntervalScheme<M>> SubstepTape<'_, M, S> {
    fn replay<'s>(&'s self, step: usize, state: &[f64], progress: &mut Progress,
        cancelled: &mut dyn FnMut() -> bool) -> Result<S::Tape<'s>, WindowError>
    {
        poll(cancelled)?;
        if progress.calls >= self.policy.budget.replayed_substeps {
            return Err(WindowError::Integrator { interval: self.interval, phase: "substep replay", diagnostic: "replayed-substep allowance exhausted".into() });
        }
        if self.policy.inner.dimension(self.model) != self.initial.len() || self.policy.inner.parameter_count(self.model) != self.parameters {
            return Err(WindowError::Invalid("substep replay model dimensions changed"));
        }
        progress.calls += 1;
        let i = self.begin+step;
        let tape = self.policy.inner.record(self.model, i, self.policy.grid.fine[i], self.policy.grid.fine[i+1], state, cancelled)?;
        checked_endpoint(&tape, self.initial.len(), self.policy.grid.fine[i+1], self.interval)?;
        if tape.accepted_steps() != self.witnesses[step].accepted_steps
            || hash(tape.end_time(), tape.endpoint(), cancelled)? != self.witnesses[step].endpoint
        {
            return Err(WindowError::Integrator { interval: self.interval, phase: "substep replay", diagnostic: format!("endpoint replay mismatch at fine step {i}") });
        }
        add_count(&mut progress.work, tape.accepted_steps())?;
        Ok(tape)
    }
    #[allow(clippy::too_many_arguments)]
    fn segment(&self, state: &[f64], begin: usize, end: usize, bar: Bar, depth: usize,
        progress: &mut Progress, cancelled: &mut dyn FnMut() -> bool) -> Result<Bar, WindowError>
    {
        poll(cancelled)?;
        progress.peak = progress.peak.max(depth);
        if end-begin == 1 {
            let tape = self.replay(begin, state, progress, cancelled)?;
            let gradient = tape.pullback(&bar.state, &bar.parameters, cancelled)?;
            poll(cancelled)?;
            if gradient.initial.len() != self.initial.len() || gradient.parameters.len() != self.parameters
                || gradient.initial.iter().chain(&gradient.parameters).any(|v| !v.is_finite())
            { return Err(WindowError::IntervalOutput { interval: self.interval, what: "substep pullback shape/values" }); }
            add_count(&mut progress.work, gradient.replayed_steps)?;
            let peak = depth.checked_add(gradient.peak_checkpoints).ok_or(WindowError::Invalid("checkpoint count overflow"))?;
            progress.peak = progress.peak.max(peak);
            return Ok(Bar { state: gradient.initial, parameters: gradient.parameters });
        }
        let span = end-begin; let mid = begin+span/2+span%2;
        let mut midpoint = copy(state)?;
        for i in begin..mid { midpoint = copy(self.replay(i, &midpoint, progress, cancelled)?.endpoint())?; }
        let bar = self.segment(&midpoint, mid, end, bar, depth+1, progress, cancelled)?;
        drop(midpoint);
        self.segment(state, begin, mid, bar, depth, progress, cancelled)
    }
}
impl<M, S: IntervalScheme<M>> IntervalTape for SubstepTape<'_, M, S> {
    fn endpoint(&self) -> &[f64] { &self.endpoint }
    fn end_time(&self) -> f64 { self.policy.grid.knots[self.interval+1] }
    fn accepted_steps(&self) -> usize { self.accepted_steps }
    fn pullback(&self, seed: &[f64], direct: &[f64], cancelled: &mut dyn FnMut() -> bool) -> Result<TrajectoryGradient, WindowError> {
        poll(cancelled)?;
        if seed.len() != self.initial.len() || direct.len() != self.parameters
            || seed.iter().chain(direct).any(|v| !v.is_finite())
        { return Err(WindowError::Invalid("invalid substep objective cotangents")); }
        let mut progress = Progress { calls: 0, work: 0, peak: 0 };
        let bar = Bar { state: copy(seed)?, parameters: copy(direct)? };
        let mut stopped = false; let mut check = || { stopped |= cancelled(); stopped };
        let bar = self.segment(&self.initial, 0, self.witnesses.len(), bar, 1, &mut progress, &mut check)?;
        poll(&mut check)?;
        // Report actual underlying step work, including inner reverse replay,
        // not only the number of coarse window intervals or fine record calls.
        Ok(TrajectoryGradient { initial: bar.state, parameters: bar.parameters,
            replayed_steps: progress.work, peak_checkpoints: progress.peak })
    }
}

#[cfg(test)]
mod tests;

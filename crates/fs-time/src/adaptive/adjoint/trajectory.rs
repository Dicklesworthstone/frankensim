//! Accepted-mesh recording and fallible binary-checkpoint RK45 pullbacks.
//!
//! Stores O(N) small step records, not O(N*n) states or stages. The backward
//! sweep uses the binary split/snapshot convention of `fs_ad::revolve`, with
//! immediate Result/cancellation propagation rather than its infallible
//! callbacks. Replays are bounded explicitly and checked against each recorded
//! endpoint's hash. Hash equality is a replay diagnostic, not authentication
//! of model derivatives or a proof of their correctness.
//!
//! The borrowed model must remain pure and unchanged, including parameters,
//! across recording, cloned checkpoints and reverse sweeps. All accepted times
//! and step sizes are frozen: this is NOT the derivative of the PI controller,
//! rejected-trial decisions, event locations or reset maps. Integral objectives
//! can be included as extra ODE state components and seeded at the endpoint.

use super::{AdjointError, OdeVjp, check, finite, poll, replay, reverse, workspace_size};
use super::super::{AdaptiveError, AdaptiveState, PiController, Workspace, commit, validate};
use fs_blake3::{Blake3, ContentHash};

pub mod samples;

#[derive(Debug, Clone)]
pub struct RecordingConfig {
    pub end: f64,
    pub rtol: f64,
    pub atol: f64,
    pub controller: PiController,
    /// Per-step scalar scratch ceiling, as in `step_vjp` (17n + 2p).
    /// Does not include tape records, parked states or the sweep cotangent.
    pub max_workspace_components: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordingStatus { ReachedEnd, AttemptLimit, RecordLimit, Cancelled }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordingReport {
    pub status: RecordingStatus,
    pub attempts: usize,
    pub accepted: usize,
    pub rejected: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct ReplayBudget {
    /// Parked state limit, including the borrowed initial state. The binary
    /// schedule requires `fs_ad::revolve::min_budget(accepted_steps)` slots.
    /// Transient RHS workspace and the live cotangent are separate.
    pub checkpoints: usize,
    /// All forward step replays, including the replay within every leaf VJP.
    pub replayed_steps: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrajectoryError {
    Step(AdjointError),
    Incomplete,
    ModelDimensionsChanged,
    Allocation,
    ReplayMismatch { step: usize },
    CheckpointLimit { required: usize, limit: usize },
    ReplayLimit,
    Observation(String),
}
impl From<AdjointError> for TrajectoryError {
    fn from(error: AdjointError) -> Self { Self::Step(error) }
}
impl From<AdaptiveError> for TrajectoryError {
    fn from(error: AdaptiveError) -> Self { Self::Step(error.into()) }
}
impl std::fmt::Display for TrajectoryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "recorded RK45 pullback failed: {self:?}")
    }
}
impl std::error::Error for TrajectoryError {}

#[derive(Debug, Clone, PartialEq)]
struct StepRecord { start: f64, end: f64, h: f64, endpoint: ContentHash }

/// A cloneable forward checkpoint. Mutation is restricted to `advance`;
/// callers cannot replace a stored state or time grid behind a pullback.
/// Reverse sweeps are atomic/retryable, not checkpointed mid-sweep: failure
/// leaves this recording intact and returns no partially accumulated gradient.
pub struct RecordedRk45<'a, M> {
    model: &'a M,
    initial: Vec<f64>,
    state: AdaptiveState,
    config: RecordingConfig,
    parameters: usize,
    records: Vec<StepRecord>,
    sample_times: Vec<f64>,
}
impl<M> Clone for RecordedRk45<'_, M> {
    fn clone(&self) -> Self {
        Self { model: self.model, initial: self.initial.clone(), state: self.state.clone(),
            config: self.config.clone(), parameters: self.parameters, records: self.records.clone(),
            sample_times: self.sample_times.clone() }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryGradient {
    pub initial: Vec<f64>,
    pub parameters: Vec<f64>,
    pub replayed_steps: usize,
    pub peak_checkpoints: usize,
}
struct Cotangent { initial: Vec<f64>, parameters: Vec<f64> }
struct Progress { replays: usize, peak: usize, budget: ReplayBudget }

fn fingerprint<Cancel: FnMut() -> bool>(values: &[f64], cancelled: &mut Cancel) -> Option<ContentHash> {
    let mut hasher = Blake3::new();
    hasher.update(b"fs-time/rk45/accepted-endpoint/v1\0");
    for chunk in values.chunks(256) {
        if cancelled() { return None; }
        for value in chunk { hasher.update(&value.to_bits().to_le_bytes()); }
    }
    let result = hasher.finalize();
    if cancelled() { None } else { Some(result) }
}

impl<'a, M: OdeVjp> RecordedRk45<'a, M> {
    pub fn new(model: &'a M, state: AdaptiveState, config: RecordingConfig) -> Result<Self, TrajectoryError> {
        validate(&state, config.end, config.rtol, config.atol, &config.controller)?;
        check(model, &state.u, config.max_workspace_components)?;
        Ok(Self { model, initial: state.u.clone(), state, config,
            parameters: model.parameter_count(), records: Vec::new(), sample_times: Vec::new() })
    }
    pub fn state(&self) -> &AdaptiveState { &self.state }
    pub fn accepted_steps(&self) -> usize { self.records.len() }
    /// (start, end, actual h), including endpoint clipping/rounding conventions.
    pub fn step_schedule(&self) -> impl ExactSizeIterator<Item = (f64, f64, f64)> + '_ {
        self.records.iter().map(|r| (r.start, r.end, r.h))
    }
    fn model_shape(&self) -> Result<(), TrajectoryError> {
        if self.model.dimension() != self.initial.len() || self.model.parameter_count() != self.parameters {
            Err(TrajectoryError::ModelDimensionsChanged)
        } else { Ok(()) }
    }

    /// Uses the production trial/commit path. Both limits are explicit:
    /// `max_attempts` is this call's work allowance; `max_records` caps the
    /// total number of live records (not allocator metadata/capacity). Increase
    /// either limit on a later call to resume without repeating accepted work.
    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self, max_attempts: usize, max_records: usize, cancelled: &mut Cancel,
    ) -> Result<RecordingReport, TrajectoryError> {
        self.model_shape()?;
        let mut report = RecordingReport { status: RecordingStatus::AttemptLimit,
            attempts: 0, accepted: 0, rejected: 0 };
        if self.state.t == self.config.end {
            report.status = RecordingStatus::ReachedEnd; return Ok(report);
        }
        let mut work = Workspace::new(self.initial.len());
        while self.state.t < self.config.end && report.attempts < max_attempts {
            if cancelled() { report.status = RecordingStatus::Cancelled; return Ok(report); }
            if self.records.len() >= max_records {
                report.status = RecordingStatus::RecordLimit; return Ok(report);
            }
            self.records.try_reserve_exact(1).map_err(|_| TrajectoryError::Allocation)?;
            let start = self.state.t;
            // Land on observation times with the production stepper, not an
            // interpolated surrogate. Duplicates share one accepted endpoint.
            let next = self.sample_times.partition_point(|time| *time <= start);
            let stop = self.sample_times.get(next).copied().unwrap_or(self.config.end);
            let h = self.state.h.min(stop - start);
            report.attempts += 1;
            let Some(trial) = work.trial(&self.state, &|t,u,out| self.model.rhs(t,u,out), stop,
                self.config.rtol, self.config.atol, &self.config.controller, cancelled)? else {
                report.status = RecordingStatus::Cancelled; return Ok(report);
            };
            let record = if trial.err <= 1.0 {
                let Some(endpoint) = fingerprint(&work.next, cancelled) else {
                    report.status = RecordingStatus::Cancelled; return Ok(report);
                };
                Some(StepRecord { start, end: trial.t, h, endpoint })
            } else { None };
            commit(&mut self.state, &mut work, trial)?;
            if let Some(record) = record {
                self.records.push(record); report.accepted += 1;
            } else { report.rejected += 1; }
        }
        if self.state.t == self.config.end { report.status = RecordingStatus::ReachedEnd; }
        Ok(report)
    }

    /// Pull back a terminal objective on a COMPLETE recording. `direct_parameters`
    /// is its explicit partial derivative with state held fixed (supply zeros
    /// when absent). Initial-condition parameter dependence is NOT implicit:
    /// contract the returned initial cotangent with that separate Jacobian.
    /// No derivatives with respect to time, tolerances or controller parameters
    /// are returned. The model and its derivatives must remain unchanged.
    pub fn pullback<Cancel: FnMut() -> bool>(
        &self, terminal: &[f64], direct_parameters: &[f64], budget: ReplayBudget, cancelled: &mut Cancel,
    ) -> Result<TrajectoryGradient, TrajectoryError> {
        self.model_shape()?;
        if self.state.t != self.config.end { return Err(TrajectoryError::Incomplete); }
        if terminal.len() != self.initial.len() || direct_parameters.len() != self.parameters
            || !finite(terminal) || !finite(direct_parameters)
        {
            return Err(AdjointError::InvalidInput("invalid objective cotangent dimensions/values").into());
        }
        workspace_size(self.initial.len(), self.parameters, self.config.max_workspace_components)?;
        let required = fs_ad::revolve::min_budget(self.records.len());
        if budget.checkpoints < required {
            return Err(TrajectoryError::CheckpointLimit { required, limit: budget.checkpoints });
        }
        poll(cancelled)?;
        let mut progress = Progress { replays: 0, peak: 0, budget };
        let bar = Cotangent { initial: terminal.to_vec(), parameters: direct_parameters.to_vec() };
        let bar = if self.records.is_empty() { bar } else {
            self.segment(&self.initial, 0, self.records.len(), bar, 1, &mut progress,
                &mut |_, _, _, _| Ok(()), cancelled)?
        };
        poll(cancelled)?;
        Ok(TrajectoryGradient { initial: bar.initial, parameters: bar.parameters,
            replayed_steps: progress.replays, peak_checkpoints: progress.peak })
    }

    fn checked_replay<Cancel: FnMut() -> bool>(
        &self, step: usize, state: &[f64], progress: &mut Progress, cancelled: &mut Cancel,
    ) -> Result<Workspace, TrajectoryError> {
        poll(cancelled)?;
        if progress.replays >= progress.budget.replayed_steps { return Err(TrajectoryError::ReplayLimit); }
        progress.replays += 1;
        let r = &self.records[step];
        let work = replay(self.model, state, r.start, r.end, r.h, cancelled)?;
        let endpoint = fingerprint(&work.next, cancelled).ok_or(AdjointError::Cancelled)?;
        if endpoint != r.endpoint { return Err(TrajectoryError::ReplayMismatch { step }); }
        Ok(work)
    }

    #[allow(clippy::too_many_arguments)]
    fn segment<Cancel, Observe>(
        &self, state: &[f64], begin: usize, end: usize, mut bar: Cotangent,
        depth: usize, progress: &mut Progress, observe: &mut Observe, cancelled: &mut Cancel,
    ) -> Result<Cotangent, TrajectoryError>
    where
        Cancel: FnMut() -> bool,
        Observe: FnMut(usize, &[f64], &mut Cotangent, &mut Cancel) -> Result<(), TrajectoryError>,
    {
        poll(cancelled)?;
        progress.peak = progress.peak.max(depth);
        if end - begin == 1 {
            let work = self.checked_replay(begin, state, progress, cancelled)?;
            observe(begin, &work.next, &mut bar, cancelled)?;
            let r = &self.records[begin];
            let result = reverse(self.model, state, r.start, r.end, r.h, &bar.initial, work, cancelled)?;
            for (value, update) in bar.parameters.iter_mut().zip(result.parameters) { *value += update; }
            if !finite(&bar.parameters) { return Err(AdjointError::NonFiniteAccumulation.into()); }
            return Ok(Cotangent { initial: result.initial, parameters: bar.parameters });
        }
        let span = end - begin;
        let mid = begin + span / 2 + span % 2;
        let mut midpoint = state.to_vec();
        for i in begin..mid { midpoint = self.checked_replay(i, &midpoint, progress, cancelled)?.next; }
        let bar = self.segment(&midpoint, mid, end, bar, depth + 1, progress, observe, cancelled)?;
        drop(midpoint);
        self.segment(state, begin, mid, bar, depth, progress, observe, cancelled)
    }
}

#[cfg(test)]
mod tests;

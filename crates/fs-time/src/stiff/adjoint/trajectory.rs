//! Memory-bounded reverse sweeps of fixed-step implicit-explicit trajectories.
//!
//! A recording retains its initial/current states and O(N) small endpoint
//! records, not N full states or Krylov histories. The reverse sweep follows
//! the binary schedule of `fs_ad::revolve`, with fallible replay and immediate
//! cancellation propagation. Each replay must reproduce the recorded endpoint
//! bits before its derivative is used. Models and preconditioners must remain
//! pure and unchanged; this replay diagnostic does not certify derivatives.

use super::{ImexAdjointError, ImexVjp, accumulate, finite};
use crate::stiff::{ImexSolveError, ImexStages, ImexState, OperatorImex2, imex_poll};
use fs_blake3::{Blake3, ContentHash};
use fs_solver::FlexiblePreconditioner;

pub mod samples;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImexRecordingConfig {
    /// Number of fixed steps in the complete trajectory, including zero.
    pub steps: usize,
    /// Per-step bound returned by `OperatorImex2::adjoint_workspace_components`.
    /// Parked checkpoints, tape records, the current/initial recorded state
    /// and the live trajectory cotangent are accounted for separately.
    pub max_workspace_components: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImexRecordingStatus {
    ReachedEnd,
    StepLimit,
    RecordLimit,
    Cancelled,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImexRecordingReport {
    pub status: ImexRecordingStatus,
    /// Complete steps appended by this call, excluding any cancelled attempt.
    pub advanced: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImexReplayBudget {
    /// Parked states, including the borrowed initial state. The binary schedule
    /// requires `fs_ad::revolve::min_budget(steps)`; a zero-step sweep needs none.
    pub checkpoints: usize,
    /// Maximum forward replays, including each leaf's primal stage calculation.
    /// Each replay performs the method's two bounded forward linear solves.
    pub forward_steps: usize,
}

#[derive(Debug, Clone)]
pub enum ImexTrajectoryError {
    Step(ImexAdjointError),
    InvalidInput(&'static str),
    Incomplete,
    ModelDimensionsChanged,
    Allocation,
    ReplayMismatch { step: usize },
    CheckpointLimit { required: usize, limit: usize },
    ReplayLimit,
    Observation(String),
}
impl From<ImexAdjointError> for ImexTrajectoryError {
    fn from(error: ImexAdjointError) -> Self {
        Self::Step(error)
    }
}
impl From<ImexSolveError> for ImexTrajectoryError {
    fn from(error: ImexSolveError) -> Self {
        Self::Step(error.into())
    }
}
impl std::fmt::Display for ImexTrajectoryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "IMEX trajectory failed: {self:?}")
    }
}
impl std::error::Error for ImexTrajectoryError {}

#[derive(Debug, Clone, PartialEq)]
struct Endpoint {
    time: f64,
    fingerprint: ContentHash,
}

/// A resumable forward recording with immutable model/method bindings.
/// Clone it to fork an accepted prefix. Reverse sweeps are atomic/retryable,
/// not resumable mid-sweep; failure never mutates this recording or returns
/// a partial gradient. Models and preconditioners may not change across forks.
pub struct RecordedImex2<'a, M, P> {
    method: OperatorImex2,
    model: &'a M,
    primal_preconditioner: &'a P,
    config: ImexRecordingConfig,
    parameters: usize,
    initial: Vec<f64>,
    initial_time: f64,
    current: Vec<f64>,
    time: f64,
    records: Vec<Endpoint>,
}
impl<M, P> Clone for RecordedImex2<'_, M, P> {
    fn clone(&self) -> Self {
        Self {
            method: self.method,
            model: self.model,
            primal_preconditioner: self.primal_preconditioner,
            config: self.config,
            parameters: self.parameters,
            initial: self.initial.clone(),
            initial_time: self.initial_time,
            current: self.current.clone(),
            time: self.time,
            records: self.records.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ImexTrajectoryGradient {
    pub initial: Vec<f64>,
    pub parameters: Vec<f64>,
    pub replayed_steps: usize,
    pub peak_checkpoints: usize,
}
struct Cotangent {
    initial: Vec<f64>,
    parameters: Vec<f64>,
}
struct Progress {
    replays: usize,
    peak: usize,
    budget: ImexReplayBudget,
}

fn fingerprint<Cancel: FnMut() -> bool>(
    time: f64,
    state: &[f64],
    cancelled: &mut Cancel,
) -> Result<ContentHash, ImexTrajectoryError> {
    let mut hash = Blake3::new();
    hash.update(b"fs-time/imex/accepted-endpoint/v1\0");
    hash.update(&time.to_bits().to_le_bytes());
    for chunk in state.chunks(256) {
        imex_poll(cancelled)?;
        for value in chunk {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
    imex_poll(cancelled)?;
    Ok(hash.finalize())
}

impl<'a, M: ImexVjp, P: FlexiblePreconditioner> RecordedImex2<'a, M, P> {
    pub fn new(
        method: OperatorImex2,
        model: &'a M,
        primal_preconditioner: &'a P,
        time: f64,
        initial: &[f64],
        config: ImexRecordingConfig,
    ) -> Result<Self, ImexTrajectoryError> {
        let parameters = model.parameter_count();
        let required = method.adjoint_workspace_components(parameters).ok_or(
            ImexTrajectoryError::InvalidInput("workspace dimension overflow"),
        )?;
        if required > config.max_workspace_components {
            return Err(ImexAdjointError::WorkspaceLimit {
                required,
                limit: config.max_workspace_components,
            }
            .into());
        }
        if model.n() != method.n
            || initial.len() != method.n
            || !finite(initial)
            || !time.is_finite()
        {
            return Err(ImexTrajectoryError::InvalidInput(
                "finite time and dimension-matched initial state required",
            ));
        }
        Ok(Self {
            method,
            model,
            primal_preconditioner,
            config,
            parameters,
            initial: initial.to_vec(),
            initial_time: time,
            current: initial.to_vec(),
            time,
            records: Vec::new(),
        })
    }

    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.current
    }
    #[must_use]
    pub fn time(&self) -> f64 {
        self.time
    }
    #[must_use]
    pub fn accepted_steps(&self) -> usize {
        self.records.len()
    }
    #[must_use]
    pub fn required_checkpoints(&self) -> usize {
        fs_ad::revolve::min_budget(self.config.steps)
    }

    fn model_shape(&self) -> Result<(), ImexTrajectoryError> {
        if self.model.n() != self.initial.len() || self.model.parameter_count() != self.parameters {
            Err(ImexTrajectoryError::ModelDimensionsChanged)
        } else {
            Ok(())
        }
    }

    /// Append at most `max_steps` complete steps, retaining at most `max_records`
    /// total endpoint records. Raise either allowance to continue. No accepted
    /// work repeats on forward continuation; a failed/cancelled attempt retains
    /// the last complete prefix. Only one step's Krylov reports is live here.
    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        max_steps: usize,
        max_records: usize,
        cancelled: &mut Cancel,
    ) -> Result<ImexRecordingReport, ImexTrajectoryError> {
        self.model_shape()?;
        let mut report = ImexRecordingReport {
            status: ImexRecordingStatus::StepLimit,
            advanced: 0,
        };
        while self.records.len() < self.config.steps && report.advanced < max_steps {
            if cancelled() {
                report.status = ImexRecordingStatus::Cancelled;
                return Ok(report);
            }
            if self.records.len() >= max_records {
                report.status = ImexRecordingStatus::RecordLimit;
                return Ok(report);
            }
            self.records
                .try_reserve(1)
                .map_err(|_| ImexTrajectoryError::Allocation)?;
            let next = self.advance_one(cancelled);
            match next {
                Ok((state, record)) => {
                    self.time = state.t;
                    self.current = state.u;
                    self.records.push(record);
                    report.advanced += 1;
                }
                Err(ImexTrajectoryError::Step(ImexAdjointError::Step(
                    ImexSolveError::Cancelled,
                ))) => {
                    report.status = ImexRecordingStatus::Cancelled;
                    return Ok(report);
                }
                Err(error) => return Err(error),
            }
        }
        if self.records.len() == self.config.steps {
            report.status = ImexRecordingStatus::ReachedEnd;
        }
        Ok(report)
    }

    fn advance_one<Cancel: FnMut() -> bool>(
        &self,
        cancelled: &mut Cancel,
    ) -> Result<(ImexState, Endpoint), ImexTrajectoryError> {
        let mut state = ImexState::new(self.time, &self.current);
        state.steps = self.records.len();
        self.method.step_controlled(
            &mut state,
            self.model,
            self.primal_preconditioner,
            &|u, out| self.model.nonlinear(u, out),
            cancelled,
        )?;
        let fingerprint = fingerprint(state.t, &state.u, cancelled)?;
        imex_poll(cancelled)?;
        let record = Endpoint {
            time: state.t,
            fingerprint,
        };
        Ok((state, record))
    }

    /// Reverse a complete recording. The terminal seed is dJ/du_N and
    /// `direct_parameters` is the objective's explicit partial with u_N fixed.
    /// Initial-condition parameter dependence is a separate caller chain rule.
    /// Stored state memory during replay is O(n log N), with O(N) compact
    /// fingerprints. Replayed stages are shared with the leaf reverse action.
    /// Budget/refusal/cancellation errors return no partially accumulated result.
    pub fn pullback<Q: FlexiblePreconditioner, Cancel: FnMut() -> bool>(
        &self,
        terminal: &[f64],
        direct_parameters: &[f64],
        adjoint_preconditioner: &Q,
        budget: ImexReplayBudget,
        cancelled: &mut Cancel,
    ) -> Result<ImexTrajectoryGradient, ImexTrajectoryError> {
        self.model_shape()?;
        if self.records.len() != self.config.steps {
            return Err(ImexTrajectoryError::Incomplete);
        }
        if terminal.len() != self.initial.len()
            || direct_parameters.len() != self.parameters
            || !finite(terminal)
            || !finite(direct_parameters)
        {
            return Err(ImexTrajectoryError::InvalidInput(
                "finite dimension-matched objective partials required",
            ));
        }
        let required = self.required_checkpoints();
        if budget.checkpoints < required {
            return Err(ImexTrajectoryError::CheckpointLimit {
                required,
                limit: budget.checkpoints,
            });
        }
        imex_poll(cancelled)?;
        let mut progress = Progress {
            replays: 0,
            peak: 0,
            budget,
        };
        let bar = Cotangent {
            initial: terminal.to_vec(),
            parameters: direct_parameters.to_vec(),
        };
        let bar = if self.records.is_empty() {
            bar
        } else {
            self.segment(
                &self.initial,
                0,
                self.records.len(),
                bar,
                1,
                &mut progress,
                adjoint_preconditioner,
                &mut |_, _, _, _| Ok(()),
                cancelled,
            )?
        };
        imex_poll(cancelled)?;
        Ok(ImexTrajectoryGradient {
            initial: bar.initial,
            parameters: bar.parameters,
            replayed_steps: progress.replays,
            peak_checkpoints: progress.peak,
        })
    }

    fn replay<Cancel: FnMut() -> bool>(
        &self,
        step: usize,
        state: &[f64],
        progress: &mut Progress,
        cancelled: &mut Cancel,
    ) -> Result<ImexStages, ImexTrajectoryError> {
        imex_poll(cancelled)?;
        if progress.replays >= progress.budget.forward_steps {
            return Err(ImexTrajectoryError::ReplayLimit);
        }
        progress.replays += 1;
        let stages = self.method.stages(
            state,
            self.model,
            self.primal_preconditioner,
            &|u, out| self.model.nonlinear(u, out),
            cancelled,
        )?;
        let record = &self.records[step];
        if fingerprint(record.time, &stages.next, cancelled)? != record.fingerprint {
            return Err(ImexTrajectoryError::ReplayMismatch { step });
        }
        Ok(stages)
    }

    #[allow(clippy::too_many_arguments)]
    fn segment<Q, Cancel, Observe>(
        &self,
        state: &[f64],
        begin: usize,
        end: usize,
        mut bar: Cotangent,
        depth: usize,
        progress: &mut Progress,
        adjoint_preconditioner: &Q,
        observe: &mut Observe,
        cancelled: &mut Cancel,
    ) -> Result<Cotangent, ImexTrajectoryError>
    where
        Q: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
        Observe:
            FnMut(usize, &[f64], &mut Cotangent, &mut Cancel) -> Result<(), ImexTrajectoryError>,
    {
        imex_poll(cancelled)?;
        progress.peak = progress.peak.max(depth);
        if end - begin == 1 {
            let stages = self.replay(begin, state, progress, cancelled)?;
            observe(end, &stages.next, &mut bar, cancelled)?;
            let gradient = self.method.reverse_stages(
                state,
                self.model,
                adjoint_preconditioner,
                &bar.initial,
                self.parameters,
                stages,
                cancelled,
            )?;
            accumulate(&mut bar.parameters, &gradient.parameters)?;
            return Ok(Cotangent {
                initial: gradient.initial,
                parameters: bar.parameters,
            });
        }
        let span = end - begin;
        let mid = begin + span / 2 + span % 2;
        let mut midpoint = state.to_vec();
        for i in begin..mid {
            midpoint = self.replay(i, &midpoint, progress, cancelled)?.next;
        }
        let bar = self.segment(
            &midpoint,
            mid,
            end,
            bar,
            depth + 1,
            progress,
            adjoint_preconditioner,
            observe,
            cancelled,
        )?;
        drop(midpoint);
        self.segment(
            state,
            begin,
            mid,
            bar,
            depth,
            progress,
            adjoint_preconditioner,
            observe,
            cancelled,
        )
    }
}

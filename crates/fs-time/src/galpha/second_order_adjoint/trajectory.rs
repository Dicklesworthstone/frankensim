//! Bounded structural trajectory recordings with parameterized transient loads.
//!
//! Forward continuation retains only live q/v/a and compact endpoint records.
//! Reverse sweeps use the binary schedule in `fs_ad::revolve`; every replayed
//! endpoint must match recorded time and q/v/a bits before its pullback is used.
//! These fingerprints diagnose replay consistency, not derivative correctness.

use super::{SecondOrderAdjointConfig, SecondOrderAdjointError, SecondOrderVjp, finite};
use crate::galpha::{
    ImplicitStepTelemetry, OperatorGeneralizedAlpha, SecondOrderState, TimeSolveError,
    structural_poll,
};
use fs_blake3::{Blake3, ContentHash};
use fs_solver::FlexiblePreconditioner;

#[path = "trajectory/samples.rs"]
pub mod samples;

/// One immutable structural model and its time-dependent applied load.
/// Both callbacks must overwrite all outputs, be pure for fixed time and
/// parameter point, and bound their own work. The load is independent of q/v/a;
/// state-dependent forces belong in the structural residual, not this callback.
pub trait StructuralTrajectoryModel: SecondOrderVjp {
    /// Applied load at the integrator's actual intermediate forcing time.
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String>;
    /// Overwrite `parameters` with `(df(time,p)/dp)^T seed` at fixed time.
    /// Explicitly write zeros for parameter-independent applied loads.
    fn forcing_vjp(&self, time: f64, seed: &[f64], parameters: &mut [f64]) -> Result<(), String>;
}

/// Fixed recording policy. Only accepted-step/record caps may change on resume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StructuralRecordingConfig {
    /// Number of fixed steps to record, including zero.
    pub steps: usize,
    /// Immutable policy for the transposed solve at each reverse leaf.
    pub adjoint: SecondOrderAdjointConfig,
    /// Per-step numerical workspace ceiling from `adjoint_workspace_components`.
    /// Initial/current states, parked checkpoints, compact endpoint records,
    /// and the live trajectory cotangent are accounted for separately.
    pub max_workspace_components: usize,
}

/// Reason a forward continuation returned its accepted prefix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructuralRecordingStatus {
    /// The configured complete trajectory is recorded.
    ReachedEnd,
    /// This call's additional-step allowance ended.
    StepLimit,
    /// The maximum retained record count was reached.
    RecordLimit,
    /// Cancellation preserved the most recent accepted prefix.
    Cancelled,
}

/// Complete steps appended by one continuation call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StructuralRecordingReport {
    /// Why forward continuation stopped.
    pub status: StructuralRecordingStatus,
    /// Complete steps appended by this call.
    pub advanced: usize,
}

/// Explicit storage and work limits for one atomic reverse sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StructuralReplayBudget {
    /// Parked q/v/a checkpoints, including the borrowed initial checkpoint.
    /// Requires `fs_ad::revolve::min_budget(steps)`; zero steps need none.
    pub checkpoints: usize,
    /// Maximum complete forward replays, including each leaf's primal solve.
    /// Leaf reverse reuses that solve instead of running Newton a second time.
    pub forward_steps: usize,
}

/// A refused sweep leaves the recording unchanged and returns no partial gradient.
#[derive(Debug, Clone)]
pub enum StructuralTrajectoryError {
    /// The production step or discrete adjoint refused, including cancellation.
    Step(SecondOrderAdjointError),
    /// Invalid input shape, controls or finite-data precondition.
    InvalidInput(&'static str),
    /// Reverse requires all configured steps to have been recorded.
    Incomplete,
    /// The borrowed model's state or parameter dimension changed.
    ModelDimensionsChanged,
    /// The endpoint record vector could not reserve capacity.
    Allocation,
    /// Recomputed endpoint bits differ from the recording.
    ReplayMismatch {
        /// Zero-based local step producing the mismatched endpoint.
        step: usize,
    },
    /// The binary replay schedule requires more parked states.
    CheckpointLimit {
        /// Minimum number of parked q/v/a checkpoints.
        required: usize,
        /// Caller-supplied checkpoint allowance.
        limit: usize,
    },
    /// No more forward replays were admitted.
    ReplayLimit,
    /// The time-dependent load provider refused.
    Forcing(String),
    /// The time-dependent load derivative provider refused.
    ForcingDerivative(String),
    /// A load provider returned nonfinite or unwritten components.
    NonFiniteForcing,
    /// A sampled objective refused or returned nonfinite/unwritten outputs.
    Observation(String),
}
impl From<SecondOrderAdjointError> for StructuralTrajectoryError {
    fn from(error: SecondOrderAdjointError) -> Self {
        Self::Step(error)
    }
}
impl From<TimeSolveError> for StructuralTrajectoryError {
    fn from(error: TimeSolveError) -> Self {
        Self::Step(error.into())
    }
}
impl std::fmt::Display for StructuralTrajectoryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "structural trajectory failed: {self:?}")
    }
}
impl std::error::Error for StructuralTrajectoryError {}

#[derive(Debug, Clone, PartialEq)]
struct Endpoint {
    time: f64,
    fingerprint: ContentHash,
}

/// Resumable forward recording tied to one immutable method/model point.
/// Forks copy live q/v/a plus O(N) small fingerprints, never cumulative solver
/// histories. Reverse sweeps are atomic/retryable, not resumable mid-sweep.
/// The caller's original `SecondOrderState` may contain history; it is borrowed
/// during construction and its history is never copied or retained here.
pub struct RecordedStructural<'a, M: ?Sized> {
    method: OperatorGeneralizedAlpha,
    model: &'a M,
    config: StructuralRecordingConfig,
    parameters: usize,
    initial: SecondOrderState,
    current: SecondOrderState,
    records: Vec<Endpoint>,
}
impl<M: ?Sized> Clone for RecordedStructural<'_, M> {
    fn clone(&self) -> Self {
        Self {
            method: self.method,
            model: self.model,
            config: self.config,
            parameters: self.parameters,
            initial: self.initial.clone(),
            current: self.current.clone(),
            records: self.records.clone(),
        }
    }
}

/// Terminal-objective gradient and consumed replay resources.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralTrajectoryGradient {
    /// Cotangent of initial displacement, holding velocity/acceleration fixed.
    pub initial_q: Vec<f64>,
    /// Cotangent of initial velocity, holding displacement/acceleration fixed.
    pub initial_v: Vec<f64>,
    /// Cotangent of initial acceleration, holding displacement/velocity fixed.
    pub initial_a: Vec<f64>,
    /// Sum of structural, load and supplied direct objective parameter partials.
    pub parameters: Vec<f64>,
    /// Forward steps replayed, including one production solve per leaf.
    pub replayed_steps: usize,
    /// Maximum number of simultaneously parked checkpoints.
    pub peak_checkpoints: usize,
}
struct Cotangent {
    q: Vec<f64>,
    v: Vec<f64>,
    a: Vec<f64>,
    parameters: Vec<f64>,
}
struct Progress {
    replays: usize,
    peak: usize,
    budget: StructuralReplayBudget,
}
struct Advanced {
    state: SecondOrderState,
    primal: ImplicitStepTelemetry,
    forcing_time: f64,
}

fn accumulate(total: &mut [f64], contribution: &[f64]) -> Result<(), StructuralTrajectoryError> {
    if !finite(contribution) {
        return Err(SecondOrderAdjointError::NonFiniteDerivative.into());
    }
    for (value, delta) in total.iter_mut().zip(contribution) {
        *value += delta;
    }
    if !finite(total) {
        return Err(SecondOrderAdjointError::NonFiniteAccumulation.into());
    }
    Ok(())
}

fn fingerprint<Cancel: FnMut() -> bool>(
    state: &SecondOrderState,
    cancelled: &mut Cancel,
) -> Result<ContentHash, StructuralTrajectoryError> {
    let mut hash = Blake3::new();
    hash.update(b"fs-time/structural/accepted-endpoint/v1\0");
    hash.update(&state.t.to_bits().to_le_bytes());
    hash.update(&(state.steps as u64).to_le_bytes());
    hash.update(&(state.q.len() as u64).to_le_bytes());
    for vector in [&state.q, &state.v, &state.a] {
        for chunk in vector.chunks(256) {
            structural_poll(cancelled)?;
            for value in chunk {
                hash.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    structural_poll(cancelled)?;
    Ok(hash.finalize())
}

impl<'a, M: StructuralTrajectoryModel + ?Sized> RecordedStructural<'a, M> {
    /// Bind an immutable model/method and copy only the live initial state.
    pub fn new(
        method: OperatorGeneralizedAlpha,
        model: &'a M,
        initial: &SecondOrderState,
        config: StructuralRecordingConfig,
    ) -> Result<Self, StructuralTrajectoryError> {
        let parameters = model.parameter_count();
        if model.dimension() != method.n
            || !initial.t.is_finite()
            || [&initial.q, &initial.v, &initial.a]
                .iter()
                .any(|state| state.len() != method.n || !finite(state))
        {
            return Err(StructuralTrajectoryError::InvalidInput(
                "finite dimension-matched initial q/v/a and time required",
            ));
        }
        if config.adjoint.restart == 0
            || config.adjoint.max_cycles == 0
            || !config.adjoint.tolerance.is_finite()
            || config.adjoint.tolerance <= 0.0
        {
            return Err(StructuralTrajectoryError::InvalidInput(
                "positive finite adjoint controls required",
            ));
        }
        initial
            .steps
            .checked_add(config.steps)
            .ok_or(StructuralTrajectoryError::InvalidInput(
                "recorded step counter overflow",
            ))?;
        method.adjoint_coefficients()?;
        let required = method
            .adjoint_workspace_components(parameters, config.adjoint)
            .ok_or(StructuralTrajectoryError::InvalidInput(
                "workspace dimension overflow",
            ))?;
        if required > config.max_workspace_components {
            return Err(SecondOrderAdjointError::WorkspaceLimit {
                required,
                limit: config.max_workspace_components,
            }
            .into());
        }
        let mut initial_copy = SecondOrderState::new(initial.t, &initial.q, &initial.v, &initial.a);
        initial_copy.steps = initial.steps;
        Ok(Self {
            method,
            model,
            config,
            parameters,
            current: initial_copy.clone(),
            initial: initial_copy,
            records: Vec::new(),
        })
    }

    /// Current accepted q/v/a, clock and absolute step counter. The history is
    /// deliberately empty; the recorder retains compact endpoint records.
    #[must_use]
    pub fn state(&self) -> &SecondOrderState {
        &self.current
    }
    #[must_use]
    /// Clock of the last completely accepted state.
    pub fn time(&self) -> f64 {
        self.current.t
    }
    #[must_use]
    /// Number of steps recorded since this recording's initial state.
    pub fn accepted_steps(&self) -> usize {
        self.records.len()
    }
    #[must_use]
    /// Minimum parked checkpoint count for reversing the complete trajectory.
    pub fn required_checkpoints(&self) -> usize {
        fs_ad::revolve::min_budget(self.config.steps)
    }

    fn model_shape(&self) -> Result<(), StructuralTrajectoryError> {
        if self.model.dimension() != self.initial.q.len()
            || self.model.parameter_count() != self.parameters
        {
            Err(StructuralTrajectoryError::ModelDimensionsChanged)
        } else {
            Ok(())
        }
    }

    /// Append at most `max_steps` complete steps and retain at most `max_records`
    /// total endpoint records. Raise those allowances to continue the accepted
    /// prefix. A failed or cancelled attempt publishes no partial step. One
    /// step's Newton/Krylov reports are live at a time and discarded afterward.
    pub fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        max_steps: usize,
        max_records: usize,
        cancelled: &mut Cancel,
    ) -> Result<StructuralRecordingReport, StructuralTrajectoryError> {
        self.model_shape()?;
        let mut report = StructuralRecordingReport {
            status: StructuralRecordingStatus::StepLimit,
            advanced: 0,
        };
        while self.records.len() < self.config.steps && report.advanced < max_steps {
            if cancelled() {
                report.status = StructuralRecordingStatus::Cancelled;
                return Ok(report);
            }
            if self.records.len() >= max_records {
                report.status = StructuralRecordingStatus::RecordLimit;
                return Ok(report);
            }
            self.records
                .try_reserve(1)
                .map_err(|_| StructuralTrajectoryError::Allocation)?;
            let attempt = self.forward(&self.current, cancelled).and_then(|advanced| {
                let record = Endpoint {
                    time: advanced.state.t,
                    fingerprint: fingerprint(&advanced.state, cancelled)?,
                };
                structural_poll(cancelled)?;
                Ok((advanced.state, record))
            });
            match attempt {
                Ok((state, record)) => {
                    self.current = state;
                    self.records.push(record);
                    report.advanced += 1;
                }
                Err(StructuralTrajectoryError::Step(SecondOrderAdjointError::Step(
                    TimeSolveError::Cancelled,
                ))) => {
                    report.status = StructuralRecordingStatus::Cancelled;
                    return Ok(report);
                }
                Err(error) => return Err(error),
            }
        }
        if self.records.len() == self.config.steps {
            report.status = StructuralRecordingStatus::ReachedEnd;
        }
        Ok(report)
    }

    fn forward<Cancel: FnMut() -> bool>(
        &self,
        initial: &SecondOrderState,
        cancelled: &mut Cancel,
    ) -> Result<Advanced, StructuralTrajectoryError> {
        structural_poll(cancelled)?;
        let forcing_time = self.method.forcing_time(initial.t)?;
        let mut forcing = vec![f64::NAN; self.initial.q.len()];
        let result = self.model.forcing(forcing_time, &mut forcing);
        structural_poll(cancelled)?;
        result.map_err(StructuralTrajectoryError::Forcing)?;
        if !finite(&forcing) {
            return Err(StructuralTrajectoryError::NonFiniteForcing);
        }
        let mut state = SecondOrderState::new(initial.t, &initial.q, &initial.v, &initial.a);
        state.steps = initial.steps;
        let primal = self
            .method
            .step_controlled(&mut state, self.model, &forcing, cancelled)?;
        if !finite(&state.q) || !finite(&state.v) || !finite(&state.a) {
            return Err(SecondOrderAdjointError::NonFiniteAccumulation.into());
        }
        state.history.clear();
        Ok(Advanced {
            state,
            primal,
            forcing_time,
        })
    }

    /// Pull back a complete recording with a terminal `(q,v,a)` cotangent and
    /// explicit direct objective-parameter partials. Parameterized forcing is
    /// differentiated at the actual intermediate time exactly once per leaf.
    /// Initial-condition/acceleration consistency derivatives belong to the
    /// caller. Stored numerical state is O(n log N); records use O(N) hashes.
    /// Checkpoints, records, initial/current state and the live cotangent are
    /// separate from the admitted per-step workspace. No mid-sweep resume,
    /// adaptive timestep/event derivative or gradient enclosure is claimed.
    pub fn pullback<P: FlexiblePreconditioner, Cancel: FnMut() -> bool>(
        &self,
        terminal: (&[f64], &[f64], &[f64]),
        direct_parameters: &[f64],
        adjoint_preconditioner: &P,
        budget: StructuralReplayBudget,
        cancelled: &mut Cancel,
    ) -> Result<StructuralTrajectoryGradient, StructuralTrajectoryError> {
        self.model_shape()?;
        if self.records.len() != self.config.steps {
            return Err(StructuralTrajectoryError::Incomplete);
        }
        if [terminal.0, terminal.1, terminal.2]
            .iter()
            .any(|seed| seed.len() != self.initial.q.len() || !finite(seed))
            || direct_parameters.len() != self.parameters
            || !finite(direct_parameters)
        {
            return Err(StructuralTrajectoryError::InvalidInput(
                "finite dimension-matched terminal and direct objective partials required",
            ));
        }
        let required = self.required_checkpoints();
        if budget.checkpoints < required {
            return Err(StructuralTrajectoryError::CheckpointLimit {
                required,
                limit: budget.checkpoints,
            });
        }
        structural_poll(cancelled)?;
        let mut progress = Progress {
            replays: 0,
            peak: 0,
            budget,
        };
        let bar = Cotangent {
            q: terminal.0.to_vec(),
            v: terminal.1.to_vec(),
            a: terminal.2.to_vec(),
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
        structural_poll(cancelled)?;
        Ok(StructuralTrajectoryGradient {
            initial_q: bar.q,
            initial_v: bar.v,
            initial_a: bar.a,
            parameters: bar.parameters,
            replayed_steps: progress.replays,
            peak_checkpoints: progress.peak,
        })
    }

    fn replay<Cancel: FnMut() -> bool>(
        &self,
        step: usize,
        state: &SecondOrderState,
        progress: &mut Progress,
        cancelled: &mut Cancel,
    ) -> Result<Advanced, StructuralTrajectoryError> {
        structural_poll(cancelled)?;
        if progress.replays >= progress.budget.forward_steps {
            return Err(StructuralTrajectoryError::ReplayLimit);
        }
        progress.replays += 1;
        let advanced = self.forward(state, cancelled)?;
        if fingerprint(&advanced.state, cancelled)? != self.records[step].fingerprint {
            return Err(StructuralTrajectoryError::ReplayMismatch { step });
        }
        Ok(advanced)
    }

    #[allow(clippy::too_many_arguments)]
    fn segment<P, Cancel, Observe>(
        &self,
        state: &SecondOrderState,
        begin: usize,
        end: usize,
        mut bar: Cotangent,
        depth: usize,
        progress: &mut Progress,
        adjoint_preconditioner: &P,
        observe: &mut Observe,
        cancelled: &mut Cancel,
    ) -> Result<Cotangent, StructuralTrajectoryError>
    where
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
        Observe: FnMut(
            usize,
            &SecondOrderState,
            &mut Cotangent,
            &mut Cancel,
        ) -> Result<(), StructuralTrajectoryError>,
    {
        structural_poll(cancelled)?;
        progress.peak = progress.peak.max(depth);
        if end - begin == 1 {
            let advanced = self.replay(begin, state, progress, cancelled)?;
            observe(end, &advanced.state, &mut bar, cancelled)?;
            let gradient = self.method.reverse_endpoint(
                state,
                self.model,
                (&bar.q, &bar.v, &bar.a),
                adjoint_preconditioner,
                self.config.adjoint,
                self.parameters,
                advanced.state,
                advanced.primal,
                cancelled,
            )?;
            accumulate(&mut bar.parameters, &gradient.parameters)?;
            let mut forcing_bar = vec![f64::NAN; self.parameters];
            structural_poll(cancelled)?;
            let result =
                self.model
                    .forcing_vjp(advanced.forcing_time, &gradient.forcing, &mut forcing_bar);
            structural_poll(cancelled)?;
            result.map_err(StructuralTrajectoryError::ForcingDerivative)?;
            accumulate(&mut bar.parameters, &forcing_bar)?;
            return Ok(Cotangent {
                q: gradient.initial_q,
                v: gradient.initial_v,
                a: gradient.initial_a,
                parameters: bar.parameters,
            });
        }
        let span = end - begin;
        let mid = begin + span / 2 + span % 2;
        let mut midpoint = state.clone();
        for i in begin..mid {
            midpoint = self.replay(i, &midpoint, progress, cancelled)?.state;
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

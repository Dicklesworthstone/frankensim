//! Box-constrained calibration of stiff, fixed-step IMEX trajectories.
//!
//! One trial records the production ARS(2,2,2) trajectory, then obtains all
//! observation derivatives in one checkpointed reverse sweep. The existing
//! SQP engine owns search, constraints, BFGS updates and the local KKT stop.
//! Initial conditions and direct sensor parameters contribute their full
//! chain rules. Step size, endpoint indices and solver policy remain fixed.

use super::box_sample;
use crate::sqp::{SqpError, SqpRunReport, SqpState, SqpStop};
use fs_solver::{FlexiblePreconditioner, LinearOp};
use fs_time::stiff::adjoint::trajectory::{
    ImexRecordingConfig, ImexRecordingReport, ImexRecordingStatus, ImexReplayBudget,
    ImexTrajectoryError, RecordedImex2, samples::SampleObjective,
};
use fs_time::stiff::adjoint::{ImexAdjointError, ImexVjp};
use fs_time::stiff::{ImexSolveConfig, ImexSolveError, OperatorImex2};

/// An immutable parameter point with matched physical and sensor derivatives.
/// All partials use the family's decision coordinates; unit transformations
/// belong in the callbacks. The actual transpose of the linear operator is
/// required, including for unequal-capacity thermal and other nonsymmetric
/// systems. Callback work and callback-owned allocations must be bounded.
pub trait ImexTransientModel: ImexVjp + SampleObjective {
    /// Finite initial state, matching `LinearOp::n`.
    fn initial_values(&self) -> &[f64];
    /// `(du_initial/dparameters)^T * initial_bar`. Overwrite every component;
    /// explicitly write zeros for parameter-independent initial conditions.
    fn initial_vjp(&self, initial_bar: &[f64], parameter_bar: &mut [f64]) -> Result<(), String>;
}

/// A pure model factory over fixed bounds and nondecreasing endpoint indices.
/// Index zero is the initial state; N is the endpoint after N steps. Repeated
/// indices permit simultaneous sensors. The observation callback's sample
/// number is its position in this array, so repeated endpoints remain distinct.
/// Bounds, observations, model meaning and preconditioners must remain unchanged
/// across study continuation. No interpolation or implicit time rounding occurs.
pub trait ImexTransientFamily {
    /// One immutable physical and observation model at a decision point.
    type Model: ImexTransientModel;
    /// Finite lower/upper faces in the same coordinates as every derivative.
    fn bounds(&self) -> &[[f64; 2]];
    /// Ordered accepted endpoints, with repetition allowed for several sensors.
    fn sample_indices(&self) -> &[usize];
    /// Construct one model; refusal is preserved as a model error.
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String>;
}

/// Fixed method and per-evaluation resource ceilings. The dense SQP cap is
/// independent of state dimension: it bounds decisions plus both box faces.
/// Tape records and checkpoint states are accounted for separately from the
/// per-step workspace. Model/factory/observation memory is caller-owned.
#[derive(Debug, Clone, Copy)]
pub struct ImexTransientConfig {
    /// Clock of the initial state.
    pub start: f64,
    /// Fixed positive step, held constant while differentiating and optimizing.
    pub step: f64,
    /// Bounds both primal and transposed stage solves.
    pub solve: ImexSolveConfig,
    /// Complete step count and per-step scalar workspace ceiling.
    pub recording: ImexRecordingConfig,
    /// Maximum physical state dimension admitted before trajectory allocation.
    pub max_state_components: usize,
    /// Maximum number of observation terms, including repeated endpoints.
    pub max_samples: usize,
    /// Forward steps permitted per objective evaluation.
    pub max_forward_steps: usize,
    /// Total accepted endpoint records permitted per evaluation.
    pub max_records: usize,
    /// Checkpoint-state and total replay-work ceilings.
    pub replay: ImexReplayBudget,
    /// Decision dimension plus both faces of every parameter bound.
    pub max_kkt_dimension: usize,
}

/// A refused trial never becomes a penalty value or a partial gradient.
#[derive(Debug, Clone)]
pub enum ImexTransientError {
    /// Invalid configuration, shape or nonfinite derivative accumulation.
    Invalid(&'static str),
    /// Model factory or initial-condition derivative refusal.
    Model(String),
    /// Original stage, derivative, observation, replay or workspace refusal.
    Trajectory(ImexTrajectoryError),
    /// A bounded forward recording stopped before the prescribed endpoint.
    ForwardStopped(ImexRecordingStatus),
    /// Cancellation was observed before a complete trial could be returned.
    Cancelled,
}
impl std::fmt::Display for ImexTransientError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "IMEX transient calibration failed: {self:?}")
    }
}
impl std::error::Error for ImexTransientError {}
impl From<ImexTrajectoryError> for ImexTransientError {
    fn from(error: ImexTrajectoryError) -> Self {
        match error {
            ImexTrajectoryError::Step(ImexAdjointError::Step(ImexSolveError::Cancelled)) => {
                Self::Cancelled
            }
            other => Self::Trajectory(other),
        }
    }
}
/// SQP failure preserving the original IMEX calibration error.
pub type ImexTransientStudyError = SqpError<ImexTransientError>;

/// One complete simulated evaluation at exactly `point`.
#[derive(Debug, Clone, PartialEq)]
pub struct ImexTransientEvaluation {
    /// Decision point shared by the value, gradient and physical endpoint.
    pub point: Vec<f64>,
    /// Sum of all requested observation objectives.
    pub value: f64,
    /// Dynamics, direct observation AND initial-condition derivatives.
    pub gradient: Vec<f64>,
    /// Physical state after the configured number of steps.
    pub final_state: Vec<f64>,
    /// Complete forward status and number of accepted steps.
    pub forward: ImexRecordingReport,
    /// Number of observation terms included in this evaluation.
    pub observations: usize,
    /// Forward step replays spent during the single reverse sweep.
    pub replayed_steps: usize,
    /// Maximum parked checkpoint states used during that sweep.
    pub peak_checkpoints: usize,
}

fn poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), ImexTransientError> {
    if cancelled() {
        Err(ImexTransientError::Cancelled)
    } else {
        Ok(())
    }
}

fn normalize(error: ImexTransientStudyError) -> ImexTransientStudyError {
    match error {
        SqpError::Evaluation(ImexTransientError::Cancelled) => SqpError::Cancelled,
        other => other,
    }
}

fn admit<F: ImexTransientFamily, Cancel: FnMut() -> bool>(
    family: &F,
    config: &ImexTransientConfig,
    cancelled: &mut Cancel,
) -> Result<(), ImexTransientError> {
    let n = family.bounds().len();
    let dim = n.checked_mul(3).filter(|d| *d <= config.max_kkt_dimension);
    if n == 0 || dim.and_then(|d| d.checked_mul(d)).is_none() {
        return Err(ImexTransientError::Invalid(
            "parameter box exceeds the dense KKT cap",
        ));
    }
    for chunk in family.bounds().chunks(256) {
        poll(cancelled)?;
        if chunk.iter().any(|b| {
            !b[0].is_finite() || !b[1].is_finite() || b[0] >= b[1] || !(b[1] - b[0]).is_finite()
        }) {
            return Err(ImexTransientError::Invalid(
                "parameter bounds must have finite positive width",
            ));
        }
    }
    let indices = family.sample_indices();
    if indices.is_empty() || indices.len() > config.max_samples {
        return Err(ImexTransientError::Invalid(
            "nonempty endpoint timetable must fit the sample cap",
        ));
    }
    let mut previous = 0;
    for chunk in indices.chunks(256) {
        poll(cancelled)?;
        for &index in chunk {
            if index < previous || index > config.recording.steps {
                return Err(ImexTransientError::Invalid(
                    "observation indices must be ordered recorded endpoints",
                ));
            }
            previous = index;
        }
    }
    if !config.start.is_finite()
        || !config.step.is_finite()
        || config.step <= 0.0
        || !config.solve.tolerance.is_finite()
        || config.solve.tolerance <= 0.0
        || config.solve.restart == 0
        || config.solve.max_cycles == 0
    {
        return Err(ImexTransientError::Invalid(
            "finite start, positive step and bounded IMEX solver policy required",
        ));
    }
    Ok(())
}

/// Evaluate one point without changing an optimizer. `Ok(None)` means a finite
/// point outside the declared box; all model, numerical, budget and cancellation
/// failures propagate. The two preconditioners are explicit: the reverse solve
/// uses `I - gamma*h*L^T` and must not silently reuse a nonsymmetric primal one.
///
/// Each usable point uses one forward recording and one sampled reverse sweep,
/// regardless of the number of sensors. The derivative is of the fixed discrete
/// trajectory, with accuracy limited by stage residuals and supplied VJPs.
pub fn evaluate_imex_transient<F, P, Q, Cancel>(
    family: &F,
    config: &ImexTransientConfig,
    point: &[f64],
    primal_preconditioner: &P,
    adjoint_preconditioner: &Q,
    cancelled: &mut Cancel,
) -> Result<Option<ImexTransientEvaluation>, ImexTransientError>
where
    F: ImexTransientFamily,
    P: FlexiblePreconditioner,
    Q: FlexiblePreconditioner,
    Cancel: FnMut() -> bool,
{
    poll(cancelled)?;
    admit(family, config, cancelled)?;
    if point.len() != family.bounds().len() || point.iter().any(|x| !x.is_finite()) {
        return Err(ImexTransientError::Invalid(
            "decision point must be finite and dimension-matched",
        ));
    }
    if point
        .iter()
        .zip(family.bounds())
        .any(|(x, b)| *x < b[0] || *x > b[1])
    {
        return Ok(None);
    }
    let model = family
        .instantiate(point)
        .map_err(ImexTransientError::Model)?;
    poll(cancelled)?;
    if model.n() == 0
        || model.n() > config.max_state_components
        || model.parameter_count() != point.len()
        || model.initial_values().len() != model.n()
    {
        return Err(ImexTransientError::Invalid(
            "model state/parameter dimensions do not match the admitted problem",
        ));
    }
    let method = OperatorImex2::new(model.n(), config.step, config.solve);
    let mut tape = RecordedImex2::new(
        method,
        &model,
        primal_preconditioner,
        config.start,
        model.initial_values(),
        config.recording,
    )?;
    let forward = tape.advance(config.max_forward_steps, config.max_records, cancelled)?;
    if forward.status == ImexRecordingStatus::Cancelled {
        return Err(ImexTransientError::Cancelled);
    }
    if forward.status != ImexRecordingStatus::ReachedEnd {
        return Err(ImexTransientError::ForwardStopped(forward.status));
    }
    let result = tape.pullback_samples(
        family.sample_indices(),
        config.max_samples,
        &model,
        adjoint_preconditioner,
        config.replay,
        cancelled,
    )?;
    let mut initial_partial = vec![f64::NAN; point.len()];
    poll(cancelled)?;
    model
        .initial_vjp(&result.gradient.initial, &mut initial_partial)
        .map_err(ImexTransientError::Model)?;
    poll(cancelled)?;
    if initial_partial.iter().any(|v| !v.is_finite()) {
        return Err(ImexTransientError::Invalid(
            "non-finite or unwritten initial-condition partial",
        ));
    }
    let mut gradient = result.gradient.parameters;
    for (total, initial) in gradient.iter_mut().zip(initial_partial) {
        *total += initial;
    }
    if gradient.iter().any(|v| !v.is_finite()) {
        return Err(ImexTransientError::Invalid(
            "non-finite total parameter gradient",
        ));
    }
    poll(cancelled)?;
    Ok(Some(ImexTransientEvaluation {
        point: point.to_vec(),
        value: result.value,
        gradient,
        final_state: tape.state().to_vec(),
        forward,
        observations: result.observations,
        replayed_steps: result.gradient.replayed_steps,
        peak_checkpoints: result.gradient.peak_checkpoints,
    }))
}

/// Accepted SQP checkpoint with exactly its matching physical evaluation.
/// A failed or interrupted trial leaves this accepted pair intact; callback
/// attempts remain charged. Continuation reuses the accepted derivatives/BFGS
/// matrix. It restarts unfinished searches and their transient recordings.
///
/// Dense SQP kernels are bounded by `max_kkt_dimension` but not internally
/// interruptible. Model and solver callbacks are bounded/polled separately.
/// Convergence is a local KKT statement using supplied derivatives, not a claim
/// of global optimality, parameter identifiability, or physical validation.
pub struct ImexTransientStudy<'a, F: ImexTransientFamily, P, Q> {
    family: &'a F,
    primal_preconditioner: &'a P,
    adjoint_preconditioner: &'a Q,
    config: ImexTransientConfig,
    state: SqpState,
    accepted: ImexTransientEvaluation,
}
impl<F: ImexTransientFamily, P, Q> Clone for ImexTransientStudy<'_, F, P, Q> {
    fn clone(&self) -> Self {
        Self {
            family: self.family,
            primal_preconditioner: self.primal_preconditioner,
            adjoint_preconditioner: self.adjoint_preconditioner,
            config: self.config,
            state: self.state.clone(),
            accepted: self.accepted.clone(),
        }
    }
}
impl<'a, F, P, Q> ImexTransientStudy<'a, F, P, Q>
where
    F: ImexTransientFamily,
    P: FlexiblePreconditioner,
    Q: FlexiblePreconditioner,
{
    /// Evaluate the initial point once and retain its matching SQP checkpoint.
    /// Both preconditioners are borrowed and remain fixed across continuation.
    pub fn new<Cancel: FnMut() -> bool>(
        family: &'a F,
        point: &[f64],
        config: ImexTransientConfig,
        primal_preconditioner: &'a P,
        adjoint_preconditioner: &'a Q,
        cancelled: &mut Cancel,
    ) -> Result<Self, ImexTransientStudyError> {
        poll(cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        admit(family, &config, cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        let mut accepted = None;
        let state = SqpState::try_new(
            point,
            config.max_kkt_dimension,
            &mut |x| {
                let Some(evaluation) = evaluate_imex_transient(
                    family,
                    &config,
                    x,
                    primal_preconditioner,
                    adjoint_preconditioner,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                let result = box_sample(
                    family.bounds(),
                    &evaluation.point,
                    evaluation.value,
                    &evaluation.gradient,
                );
                accepted = Some(evaluation);
                Ok(Some(result))
            },
            None,
        )
        .map_err(normalize)?;
        let accepted = accepted.ok_or(SqpError::Invalid(
            "missing complete initial IMEX trajectory",
        ))?;
        poll(cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        Ok(Self {
            family,
            primal_preconditioner,
            adjoint_preconditioner,
            config,
            state,
            accepted,
        })
    }
    /// Accepted optimizer checkpoint, including cumulative work counters.
    #[must_use]
    pub fn optimizer(&self) -> &SqpState {
        &self.state
    }
    /// Physical value and derivatives at exactly `optimizer().point()`.
    #[must_use]
    pub fn accepted(&self) -> &ImexTransientEvaluation {
        &self.accepted
    }

    /// Continue with a cumulative sample-attempt ceiling, including failed and
    /// rejected trials. Per-attempt forward, record, replay and workspace caps
    /// remain those admitted at construction; limits never refund spent work.
    pub fn run<Cancel: FnMut() -> bool>(
        &mut self,
        tolerance: f64,
        additional_iterations: usize,
        maximum_evaluations: usize,
        cancelled: &mut Cancel,
    ) -> Result<SqpRunReport, ImexTransientStudyError> {
        poll(cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        let mut report = self.advance(tolerance, 0, maximum_evaluations, cancelled)?;
        for _ in 0..additional_iterations {
            if report.stop != SqpStop::IterationLimit {
                break;
            }
            poll(cancelled)
                .map_err(SqpError::Evaluation)
                .map_err(normalize)?;
            report = self.advance(tolerance, 1, maximum_evaluations, cancelled)?;
        }
        poll(cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        Ok(report)
    }

    fn advance<Cancel: FnMut() -> bool>(
        &mut self,
        tolerance: f64,
        steps: usize,
        maximum_evaluations: usize,
        cancelled: &mut Cancel,
    ) -> Result<SqpRunReport, ImexTransientStudyError> {
        let mut candidate = None;
        let outcome = self.state.try_run(
            &mut |point| {
                let Some(evaluation) = evaluate_imex_transient(
                    self.family,
                    &self.config,
                    point,
                    self.primal_preconditioner,
                    self.adjoint_preconditioner,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                let result = box_sample(
                    self.family.bounds(),
                    &evaluation.point,
                    evaluation.value,
                    &evaluation.gradient,
                );
                candidate = Some(evaluation);
                Ok(Some(result))
            },
            tolerance,
            steps,
            maximum_evaluations,
            None,
        );
        // This call permits at most one accepted step. Publish its complete
        // physical result only when SQP retained that exact decision point.
        if let Some(evaluation) = candidate
            && evaluation.point.as_slice() == self.state.point()
        {
            self.accepted = evaluation;
        }
        outcome.map_err(normalize)
    }
}

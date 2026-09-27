//! Bounded structural-trajectory calibration through the existing SQP engine.
//!
//! Each trial records one production generalized-alpha trajectory and obtains
//! all displacement, velocity and acceleration observation derivatives in one
//! checkpointed reverse sweep. Model, time-dependent forcing, direct sensor
//! and initial-state derivatives share the same decision coordinates. Initial
//! acceleration consistency belongs to the model's initial-state chain rule.
//! This numerical binding supplies no plate/shell or apparatus validation.

use super::box_sample;
use crate::sqp::{SqpError, SqpRunReport, SqpState, SqpStop};
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::second_order_adjoint::trajectory::{
    RecordedStructural, StructuralRecordingConfig, StructuralRecordingReport,
    StructuralRecordingStatus, StructuralReplayBudget, StructuralTrajectoryError,
    StructuralTrajectoryModel, samples::StructuralSampleObjective,
};
use fs_time::galpha::second_order_adjoint::{SecondOrderAdjointError, SecondOrderVjp};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    TimeSolveError,
};

/// One immutable structural parameter point and its observation model.
/// All derivatives use the family's decision coordinates. q/v/a and the load
/// must describe a consistent physical initialization when that is required
/// by the model. Callback work and callback-owned allocations are bounded by
/// the caller; the trajectory workspace does not account for them.
pub trait StructuralTransientModel: StructuralTrajectoryModel + StructuralSampleObjective {
    /// Finite initial q/v/a at the configured start, with zero accepted steps
    /// and empty history. Construct parameter-dependent acceleration from the
    /// same mass, damping, internal-force and initial-load model used later.
    fn initial_state(&self) -> &SecondOrderState;
    /// Total initial-state pullback into the decision parameters. Overwrite
    /// every output. When `a0=M^-1(f0-C v0-r(q0))`, this includes the derivative
    /// of that solve in addition to explicit q0/v0 dependence; dropping `abar`
    /// silently loses part of the transient objective's parameter gradient.
    fn initial_vjp(
        &self,
        qbar: &[f64],
        vbar: &[f64],
        abar: &[f64],
        parameter_bar: &mut [f64],
    ) -> Result<(), String>;
}

/// Pure model factory over fixed finite parameter bounds and ordered endpoint
/// indices. Zero denotes the initial state, and N the endpoint after N steps.
/// Repeated indices permit several sensors or state components at one time.
/// Model meaning, observation data, initial clock, bounds and preconditioner
/// must remain unchanged across continuation. No interpolation is inferred.
pub trait StructuralTransientFamily {
    /// Immutable physical/observation model at one decision point.
    type Model: StructuralTransientModel;
    /// Lower and upper faces in the coordinates used by all partials.
    fn bounds(&self) -> &[[f64; 2]];
    /// Nondecreasing accepted endpoints, with repetitions allowed.
    fn sample_indices(&self) -> &[usize];
    /// Instantiate exactly one parameter point, preserving model refusals.
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String>;
}

/// Fixed integration policy and independent per-evaluation resource ceilings.
/// Initial/current/checkpoint q/v/a, endpoint records and step scratch have
/// separate bounds. Model/factory/observation storage is caller-owned. The
/// SQP dimension cap bounds decisions plus both faces of their parameter box.
#[derive(Debug, Clone, Copy)]
pub struct StructuralTransientConfig {
    /// Exact initial clock, independent of the decision point.
    pub start: f64,
    /// Positive step held fixed while differentiating and optimizing.
    pub step: f64,
    /// Generalized-alpha high-frequency spectral radius in [0,1].
    pub rho_inf: f64,
    /// Existing Newton/Krylov primal controls, including the outer work limit.
    pub solve: ImplicitSolveConfig,
    /// Trajectory length, adjoint solve policy and per-step workspace ceiling.
    pub recording: StructuralRecordingConfig,
    /// Maximum dimension of EACH of displacement, velocity and acceleration.
    pub max_state_components: usize,
    /// Maximum observation terms, including repetitions at one endpoint.
    pub max_samples: usize,
    /// Maximum complete forward steps spent per objective evaluation.
    pub max_forward_steps: usize,
    /// Maximum retained accepted endpoint records per evaluation.
    pub max_records: usize,
    /// Checkpoint-state and forward-replay work ceilings for the reverse sweep.
    pub replay: StructuralReplayBudget,
    /// Decision dimension plus both box faces, before dense SQP allocations.
    pub max_kkt_dimension: usize,
}

/// Refusal never becomes a penalty value or a partially usable derivative.
#[derive(Debug, Clone)]
pub enum StructuralTransientError {
    /// Invalid configuration, shape or nonfinite derivative accumulation.
    Invalid(&'static str),
    /// Model factory or initial-state derivative refusal.
    Model(String),
    /// Original forcing, integration, derivative, observation or replay error.
    Trajectory(StructuralTrajectoryError),
    /// A forward work/record limit prevented a complete trajectory.
    ForwardStopped(StructuralRecordingStatus),
    /// Cancellation was observed before a complete trial could be returned.
    Cancelled,
}
impl std::fmt::Display for StructuralTransientError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "structural transient calibration failed: {self:?}")
    }
}
impl std::error::Error for StructuralTransientError {}
impl From<StructuralTrajectoryError> for StructuralTransientError {
    fn from(error: StructuralTrajectoryError) -> Self {
        match error {
            StructuralTrajectoryError::Step(SecondOrderAdjointError::Step(
                TimeSolveError::Cancelled,
            )) => Self::Cancelled,
            other => Self::Trajectory(other),
        }
    }
}
/// SQP failure preserving the original structural calibration refusal.
pub type StructuralTransientStudyError = SqpError<StructuralTransientError>;

/// Complete objective and physical endpoint at exactly one decision point.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralTransientEvaluation {
    /// Decision point shared by all retained values and derivatives.
    pub point: Vec<f64>,
    /// Sum of the requested scalar observation objectives.
    pub value: f64,
    /// Model, forcing, direct sensor and consistent-initial-state derivatives.
    pub gradient: Vec<f64>,
    /// Complete endpoint q/v/a, clock and counter, with empty history.
    pub final_state: SecondOrderState,
    /// Forward completion status and accepted step count.
    pub forward: StructuralRecordingReport,
    /// Number of observation terms accumulated in the reverse sweep.
    pub observations: usize,
    /// Forward steps replayed during the single checkpointed reverse sweep.
    pub replayed_steps: usize,
    /// Maximum number of parked q/v/a checkpoint states.
    pub peak_checkpoints: usize,
}

fn poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), StructuralTransientError> {
    if cancelled() {
        Err(StructuralTransientError::Cancelled)
    } else {
        Ok(())
    }
}
fn normalize(error: StructuralTransientStudyError) -> StructuralTransientStudyError {
    match error {
        SqpError::Evaluation(StructuralTransientError::Cancelled) => SqpError::Cancelled,
        other => other,
    }
}

fn admit<F: StructuralTransientFamily, Cancel: FnMut() -> bool>(
    family: &F,
    config: &StructuralTransientConfig,
    cancelled: &mut Cancel,
) -> Result<(), StructuralTransientError> {
    let n = family.bounds().len();
    let dimension = n.checked_mul(3).filter(|d| *d <= config.max_kkt_dimension);
    if n == 0 || dimension.and_then(|d| d.checked_mul(d)).is_none() {
        return Err(StructuralTransientError::Invalid(
            "parameter box exceeds the dense KKT cap",
        ));
    }
    for chunk in family.bounds().chunks(256) {
        poll(cancelled)?;
        if chunk.iter().any(|b| {
            !b[0].is_finite() || !b[1].is_finite() || b[0] >= b[1] || !(b[1] - b[0]).is_finite()
        }) {
            return Err(StructuralTransientError::Invalid(
                "parameter bounds must have finite positive width",
            ));
        }
    }
    let indices = family.sample_indices();
    if indices.is_empty() || indices.len() > config.max_samples {
        return Err(StructuralTransientError::Invalid(
            "nonempty endpoint timetable must fit the sample cap",
        ));
    }
    let mut previous = 0;
    for chunk in indices.chunks(256) {
        poll(cancelled)?;
        for &index in chunk {
            if index < previous || index > config.recording.steps {
                return Err(StructuralTransientError::Invalid(
                    "observation indices must be ordered recorded endpoints",
                ));
            }
            previous = index;
        }
    }
    let adjoint = config.recording.adjoint;
    if !config.start.is_finite()
        || !config.step.is_finite()
        || config.step <= 0.0
        || !(0.0..=1.0).contains(&config.rho_inf)
        || config.solve.max_newton_iterations == 0
        || adjoint.restart == 0
        || adjoint.max_cycles == 0
        || !adjoint.tolerance.is_finite()
        || adjoint.tolerance <= 0.0
    {
        return Err(StructuralTransientError::Invalid(
            "finite clock, positive step/solve budgets and spectral radius in [0,1] required",
        ));
    }
    Ok(())
}

/// Evaluate one point without changing an optimizer. `Ok(None)` means only a
/// finite point outside its declared parameter box. Every other refusal
/// propagates. The explicit adjoint preconditioner must suit the transposed
/// effective structural operator; the primal uses the existing Newton solver's
/// identity inner preconditioner, without a symmetry assumption.
///
/// One forward recording and one sampled reverse supply all observations.
/// Derivatives hold step, spectral radius, initial clock and endpoint timetable
/// fixed and are limited by supplied partials and primal/adjoint residuals.
#[allow(clippy::too_many_lines)]
pub fn evaluate_structural_transient<F, P, Cancel>(
    family: &F,
    config: &StructuralTransientConfig,
    point: &[f64],
    adjoint_preconditioner: &P,
    cancelled: &mut Cancel,
) -> Result<Option<StructuralTransientEvaluation>, StructuralTransientError>
where
    F: StructuralTransientFamily,
    P: FlexiblePreconditioner,
    Cancel: FnMut() -> bool,
{
    poll(cancelled)?;
    admit(family, config, cancelled)?;
    if point.len() != family.bounds().len() || point.iter().any(|x| !x.is_finite()) {
        return Err(StructuralTransientError::Invalid(
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
        .map_err(StructuralTransientError::Model)?;
    poll(cancelled)?;
    let initial = model.initial_state();
    if model.dimension() == 0
        || model.dimension() > config.max_state_components
        || model.parameter_count() != point.len()
        || initial.t.to_bits() != config.start.to_bits()
        || initial.steps != 0
        || !initial.history.is_empty()
        || [&initial.q, &initial.v, &initial.a].iter().any(|values| {
            values.len() != model.dimension() || values.iter().any(|v| !v.is_finite())
        })
    {
        return Err(StructuralTransientError::Invalid(
            "model dimensions and fresh finite initial q/v/a must match the admitted clock/problem",
        ));
    }
    let method =
        OperatorGeneralizedAlpha::new(model.dimension(), config.step, config.rho_inf, config.solve);
    let mut tape = RecordedStructural::new(method, &model, initial, config.recording)?;
    let forward = tape.advance(config.max_forward_steps, config.max_records, cancelled)?;
    if forward.status == StructuralRecordingStatus::Cancelled {
        return Err(StructuralTransientError::Cancelled);
    }
    if forward.status != StructuralRecordingStatus::ReachedEnd {
        return Err(StructuralTransientError::ForwardStopped(forward.status));
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
        .initial_vjp(
            &result.gradient.initial_q,
            &result.gradient.initial_v,
            &result.gradient.initial_a,
            &mut initial_partial,
        )
        .map_err(StructuralTransientError::Model)?;
    poll(cancelled)?;
    if initial_partial.iter().any(|v| !v.is_finite()) {
        return Err(StructuralTransientError::Invalid(
            "non-finite or unwritten initial-state partial",
        ));
    }
    let mut gradient = result.gradient.parameters;
    for (total, initial) in gradient.iter_mut().zip(initial_partial) {
        *total += initial;
    }
    if gradient.iter().any(|v| !v.is_finite()) {
        return Err(StructuralTransientError::Invalid(
            "non-finite total parameter gradient",
        ));
    }
    poll(cancelled)?;
    Ok(Some(StructuralTransientEvaluation {
        point: point.to_vec(),
        value: result.value,
        gradient,
        final_state: tape.state().clone(),
        forward,
        observations: result.observations,
        replayed_steps: result.gradient.replayed_steps,
        peak_checkpoints: result.gradient.peak_checkpoints,
    }))
}

/// Accepted SQP checkpoint paired with exactly its matching physical result.
/// Failed/cancelled attempts preserve that pair and remain charged. Cloning or
/// continuation retains accepted derivatives and the existing BFGS state;
/// unfinished searches restart their trial recordings without refunding work.
///
/// Dense QP/BFGS kernels have dimension bounds but are not internally
/// interruptible. Callback/step/replay work has independent bounds. A local
/// KKT stop is not a global optimum, identifiability or physical-validation
/// statement; in particular this does not validate plate/shell dynamics.
pub struct StructuralTransientStudy<'a, F: StructuralTransientFamily, P> {
    family: &'a F,
    adjoint_preconditioner: &'a P,
    config: StructuralTransientConfig,
    state: SqpState,
    accepted: StructuralTransientEvaluation,
}
impl<F: StructuralTransientFamily, P> Clone for StructuralTransientStudy<'_, F, P> {
    fn clone(&self) -> Self {
        Self {
            family: self.family,
            adjoint_preconditioner: self.adjoint_preconditioner,
            config: self.config,
            state: self.state.clone(),
            accepted: self.accepted.clone(),
        }
    }
}
impl<'a, F: StructuralTransientFamily, P: FlexiblePreconditioner>
    StructuralTransientStudy<'a, F, P>
{
    /// Evaluate the initial point once and retain its matching SQP sample.
    /// The borrowed model family and adjoint preconditioner remain fixed.
    pub fn new<Cancel: FnMut() -> bool>(
        family: &'a F,
        point: &[f64],
        config: StructuralTransientConfig,
        adjoint_preconditioner: &'a P,
        cancelled: &mut Cancel,
    ) -> Result<Self, StructuralTransientStudyError> {
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
                let Some(evaluation) = evaluate_structural_transient(
                    family,
                    &config,
                    x,
                    adjoint_preconditioner,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                let sample = box_sample(
                    family.bounds(),
                    &evaluation.point,
                    evaluation.value,
                    &evaluation.gradient,
                );
                accepted = Some(evaluation);
                Ok(Some(sample))
            },
            None,
        )
        .map_err(normalize)?;
        let accepted = accepted.ok_or(SqpError::Invalid(
            "missing complete initial structural trajectory",
        ))?;
        poll(cancelled)
            .map_err(SqpError::Evaluation)
            .map_err(normalize)?;
        Ok(Self {
            family,
            adjoint_preconditioner,
            config,
            state,
            accepted,
        })
    }

    /// Accepted optimizer checkpoint and cumulative work counters.
    #[must_use]
    pub fn optimizer(&self) -> &SqpState {
        &self.state
    }
    /// Physical value/gradient and endpoint at exactly `optimizer().point()`.
    #[must_use]
    pub fn accepted(&self) -> &StructuralTransientEvaluation {
        &self.accepted
    }

    /// Continue with an additional accepted-iteration allowance and cumulative
    /// sample-attempt limit. Failed/rejected attempts count toward that limit;
    /// each attempt obeys the fixed trajectory/workspace/replay caps above.
    pub fn run<Cancel: FnMut() -> bool>(
        &mut self,
        tolerance: f64,
        additional_iterations: usize,
        maximum_evaluations: usize,
        cancelled: &mut Cancel,
    ) -> Result<SqpRunReport, StructuralTransientStudyError> {
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
    ) -> Result<SqpRunReport, StructuralTransientStudyError> {
        let mut candidate = None;
        let outcome = self.state.try_run(
            &mut |point| {
                let Some(evaluation) = evaluate_structural_transient(
                    self.family,
                    &self.config,
                    point,
                    self.adjoint_preconditioner,
                    cancelled,
                )?
                else {
                    return Ok(None);
                };
                let sample = box_sample(
                    self.family.bounds(),
                    &evaluation.point,
                    evaluation.value,
                    &evaluation.gradient,
                );
                candidate = Some(evaluation);
                Ok(Some(sample))
            },
            tolerance,
            steps,
            maximum_evaluations,
            None,
        );
        // At most one step can commit here; a rejected trial cannot replace
        // the physical result associated with the retained SQP point.
        if let Some(evaluation) = candidate
            && evaluation.point.as_slice() == self.state.point()
        {
            self.accepted = evaluation;
        }
        outcome.map_err(normalize)
    }
}

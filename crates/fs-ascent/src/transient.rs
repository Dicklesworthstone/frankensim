//! Box-constrained transient calibration through the existing SQP engine.
//!
//! Each usable trial records one production RK45 trajectory and evaluates all
//! sensor terms with one checkpointed discrete adjoint. SQP owns the search,
//! BFGS updates and local KKT stop; this module does not implement a second
//! optimizer, finite-difference trajectories, or infer parameter identifiability.

use crate::sqp::{SqpError, SqpRunReport, SqpSample, SqpState, SqpStop};
use fs_time::AdaptiveState;
use fs_time::adaptive::adjoint::{AdjointError, OdeVjp};
use fs_time::adaptive::adjoint::trajectory::{
    RecordedRk45, RecordingConfig, RecordingReport, RecordingStatus, ReplayBudget, TrajectoryError,
    samples::SampleObjective,
};

/// One immutable parameter point, including its observation model and data.
/// `OdeVjp` and `SampleObjective` partials use the SAME decision coordinates as
/// the family's bounds. Apply physical-unit/parameter-transform chain rules in
/// those callbacks; the study never guesses scales or units.
pub trait TransientModel: OdeVjp + SampleObjective {
    fn initial_values(&self) -> &[f64];
    /// (dx_initial/dparameters)^T * initial_bar. Overwrite every entry; supply
    /// explicit zeros for parameter-independent initial conditions.
    fn initial_vjp(&self, initial_bar: &[f64], parameter_bar: &mut [f64]) -> Result<(), String>;
}

/// A pure model factory over a fixed, nondecreasing observation timetable and
/// finite parameter box. All model meaning, data, times and bounds must remain
/// unchanged on continuation. No material/posterior authority is minted here.
pub trait TransientFamily {
    type Model: TransientModel;
    fn bounds(&self) -> &[[f64; 2]];
    fn sample_times(&self) -> &[f64];
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String>;
}

#[derive(Debug, Clone)]
pub struct TransientConfig {
    pub start: f64,
    pub initial_step: f64,
    pub recording: RecordingConfig,
    pub max_state_components: usize,
    pub max_samples: usize,
    /// Per objective evaluation; failed attempts do not produce a usable trial.
    pub max_attempts: usize,
    pub max_records: usize,
    pub replay: ReplayBudget,
    /// Decisions plus both box faces; SQP remains a small-dense optimizer.
    pub max_kkt_dimension: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TransientError {
    Invalid(&'static str),
    Model(String),
    Trajectory(TrajectoryError),
    ForwardStopped(RecordingStatus),
    Cancelled,
}
impl std::fmt::Display for TransientError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "transient calibration failed: {self:?}")
    }
}
impl std::error::Error for TransientError {}
impl From<TrajectoryError> for TransientError {
    fn from(error: TrajectoryError) -> Self {
        if error == TrajectoryError::Step(AdjointError::Cancelled) { Self::Cancelled }
        else { Self::Trajectory(error) }
    }
}
pub type TransientStudyError = SqpError<TransientError>;

#[derive(Debug, Clone, PartialEq)]
pub struct TransientEvaluation {
    pub point: Vec<f64>,
    pub value: f64,
    /// Includes RHS, direct observation, AND initial-condition derivatives.
    pub gradient: Vec<f64>,
    pub final_state: Vec<f64>,
    pub forward: RecordingReport,
    pub observations: usize,
    pub replayed_steps: usize,
    pub peak_checkpoints: usize,
}

fn poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), TransientError> {
    if cancelled() { Err(TransientError::Cancelled) } else { Ok(()) }
}
fn normalize(error: TransientStudyError) -> TransientStudyError {
    match error { SqpError::Evaluation(TransientError::Cancelled) => SqpError::Cancelled, other => other }
}

fn admit<F: TransientFamily>(family: &F, config: &TransientConfig) -> Result<(), TransientError> {
    let n = family.bounds().len();
    let dim = n.checked_mul(3).filter(|d| *d <= config.max_kkt_dimension);
    if n == 0 || dim.and_then(|d| d.checked_mul(d)).is_none() {
        return Err(TransientError::Invalid("parameter box exceeds the dense KKT cap"));
    }
    if family.bounds().iter().any(|b| !b[0].is_finite() || !b[1].is_finite()
        || b[0] >= b[1] || !(b[1] - b[0]).is_finite()) {
        return Err(TransientError::Invalid("parameter bounds must have finite positive width"));
    }
    let times = family.sample_times();
    if times.is_empty() || times.len() > config.max_samples {
        return Err(TransientError::Invalid("nonempty observation timetable must fit the sample cap"));
    }
    if !config.start.is_finite() || !config.recording.end.is_finite()
        || config.recording.end < config.start || !(config.recording.end - config.start).is_finite()
        || !config.initial_step.is_finite() || config.initial_step <= 0.0
        || times.iter().any(|t| !t.is_finite() || *t < config.start || *t > config.recording.end)
        || times.windows(2).any(|p| p[0] > p[1])
    {
        return Err(TransientError::Invalid("invalid forward interval, step or observation timetable"));
    }
    Ok(())
}

/// Evaluate a point without mutating an optimizer. `Ok(None)` means only that
/// a finite trial lies outside its declared box. Numerical/model/budget errors
/// propagate, never become a fake penalty. Caller-owned model/callback memory
/// is outside the RK45 scratch and small-dense SQP dimension caps.
pub fn evaluate_transient<F: TransientFamily, Cancel: FnMut() -> bool>(
    family: &F, config: &TransientConfig, point: &[f64], cancelled: &mut Cancel,
) -> Result<Option<TransientEvaluation>, TransientError> {
    poll(cancelled)?;
    admit(family, config)?;
    if point.len() != family.bounds().len() || point.iter().any(|x| !x.is_finite()) {
        return Err(TransientError::Invalid("decision point must be finite and dimension-matched"));
    }
    if point.iter().zip(family.bounds()).any(|(x, b)| *x < b[0] || *x > b[1]) { return Ok(None); }
    let model = family.instantiate(point).map_err(TransientError::Model)?;
    poll(cancelled)?;
    if model.dimension() == 0 || model.dimension() > config.max_state_components
        || model.parameter_count() != point.len() || model.initial_values().len() != model.dimension()
    {
        return Err(TransientError::Invalid("model state/parameter dimensions do not match the admitted problem"));
    }
    let initial = AdaptiveState::new(config.start, model.initial_values(), config.initial_step);
    let mut tape = RecordedRk45::new_sampled(&model, initial, config.recording.clone(),
        family.sample_times(), config.max_samples)?;
    let forward = tape.advance(config.max_attempts, config.max_records, cancelled)?;
    if forward.status == RecordingStatus::Cancelled { return Err(TransientError::Cancelled); }
    if forward.status != RecordingStatus::ReachedEnd { return Err(TransientError::ForwardStopped(forward.status)); }
    let result = tape.pullback_samples(&model, config.replay, cancelled)?;
    let mut initial_partial = vec![f64::NAN; point.len()];
    poll(cancelled)?;
    model.initial_vjp(&result.gradient.initial, &mut initial_partial).map_err(TransientError::Model)?;
    poll(cancelled)?;
    if initial_partial.iter().any(|v| !v.is_finite()) {
        return Err(TransientError::Invalid("non-finite or unwritten initial-condition partial"));
    }
    let mut gradient = result.gradient.parameters;
    for (total, initial) in gradient.iter_mut().zip(initial_partial) { *total += initial; }
    if gradient.iter().any(|v| !v.is_finite()) {
        return Err(TransientError::Invalid("non-finite total parameter gradient"));
    }
    poll(cancelled)?;
    Ok(Some(TransientEvaluation { point: point.to_vec(), value: result.value, gradient,
        final_state: tape.state().u.clone(), forward, observations: result.observations,
        replayed_steps: result.gradient.replayed_steps, peak_checkpoints: result.gradient.peak_checkpoints,
    }))
}

fn sample(bounds: &[[f64; 2]], evaluation: &TransientEvaluation) -> SqpSample {
    let n = bounds.len(); let mut ci = Vec::with_capacity(2*n); let mut ji = vec![0.0; 2*n*n];
    for (i, (b, x)) in bounds.iter().zip(&evaluation.point).enumerate() {
        ci.push(b[0]-x); ci.push(x-b[1]); ji[2*i*n+i]=-1.0; ji[(2*i+1)*n+i]=1.0;
    }
    SqpSample { f: evaluation.value, gradient: evaluation.gradient.clone(), ce: Vec::new(), ci, je: Vec::new(), ji }
}

/// Retains the accepted SQP state and exactly its matching simulated evaluation.
/// Bounds are real SQP inequalities, not post-step clamping. Trial evaluations
/// use fresh trajectories, while accepted-step continuation reuses the cached
/// value, gradient, BFGS matrix and cumulative evaluation counter.
///
/// `run`'s cumulative evaluation limit counts failed and rejected callbacks;
/// each callback has the fixed simulation/replay caps above. Interrupted
/// searches restart, not refund work. Dense QP/BFGS kernels are bounded but not
/// internally cancellable. Convergence is a local KKT statement using supplied
/// derivatives and frozen-mesh adjoints, not validation or identifiability.
pub struct TransientStudy<'a, F: TransientFamily> {
    family: &'a F,
    config: TransientConfig,
    state: SqpState,
    accepted: TransientEvaluation,
}
impl<'a, F: TransientFamily> TransientStudy<'a, F> {
    pub fn new<Cancel: FnMut() -> bool>(
        family: &'a F, point: &[f64], config: TransientConfig, cancelled: &mut Cancel,
    ) -> Result<Self, TransientStudyError> {
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        admit(family, &config).map_err(SqpError::Evaluation)?;
        let mut accepted = None;
        let state = SqpState::try_new(point, config.max_kkt_dimension, &mut |x| {
            let Some(evaluation) = evaluate_transient(family, &config, x, cancelled)? else { return Ok(None); };
            let result = sample(family.bounds(), &evaluation); accepted = Some(evaluation); Ok(Some(result))
        }, None).map_err(normalize)?;
        let accepted = accepted.ok_or(SqpError::Invalid("missing complete initial trajectory"))?;
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(Self { family, config, state, accepted })
    }
    pub fn optimizer(&self) -> &SqpState { &self.state }
    pub fn accepted(&self) -> &TransientEvaluation { &self.accepted }

    pub fn run<Cancel: FnMut() -> bool>(
        &mut self, tolerance: f64, additional_iterations: usize,
        maximum_evaluations: usize, cancelled: &mut Cancel,
    ) -> Result<SqpRunReport, TransientStudyError> {
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        let mut report = self.advance(tolerance, 0, maximum_evaluations, cancelled)?;
        for _ in 0..additional_iterations {
            if report.stop != SqpStop::IterationLimit { break; }
            poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
            report = self.advance(tolerance, 1, maximum_evaluations, cancelled)?;
        }
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(report)
    }

    fn advance<Cancel: FnMut() -> bool>(
        &mut self, tolerance: f64, steps: usize, maximum_evaluations: usize, cancelled: &mut Cancel,
    ) -> Result<SqpRunReport, TransientStudyError> {
        let family = self.family; let config = &self.config; let mut candidate = None;
        let outcome = self.state.try_run(&mut |point| {
            let Some(evaluation) = evaluate_transient(family, config, point, cancelled)? else { return Ok(None); };
            let result = sample(family.bounds(), &evaluation); candidate = Some(evaluation); Ok(Some(result))
        }, tolerance, steps, maximum_evaluations, None);
        // At most one accepted step: never associate a rejected trial with the
        // retained optimizer. Errors keep the last accepted physical result.
        if let Some(evaluation) = candidate {
            if evaluation.point.as_slice() == self.state.point() { self.accepted = evaluation; }
        }
        outcome.map_err(normalize)
    }
}

#[cfg(test)]
mod tests;

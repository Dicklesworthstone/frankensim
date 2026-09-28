//! State-form weak-constraint 4D-Var over checked discrete interval maps.
//!
//! Each knot state is a decision variable. The objective is the observation
//! loss plus Gaussian penalties on the initial departure and on
//! `x[k+1] - M[k](x[k])`. Model-error scales describe ENDPOINT state increments,
//! not white-noise spectral densities; no implicit dt scaling is invented.
//! Background and model-error covariances are fixed, diagonal and mutually
//! independent. A joint observation callback can describe cross-time noise.
//!
//! Each interval uses the actual recorded numerical map and its checkpointed
//! transpose action, not differentiation of solver/controller iterations or a
//! finite-difference trajectory. The accepted RK mesh is frozen in a gradient.
//! Only one interval tape is live; neither a state covariance nor a trajectory
//! Jacobian is assembled. The state-form control vector still costs O(knots*n).
//!
//! This is a numerical objective, not an authenticated/certified adjoint, an
//! exact nonlinear posterior, or an uncertainty-coverage/physical certificate.
//! Reference: Tremolet, QJRMS 132 (2006), 2483-2504, doi:10.1256/qj.05.224.

use fs_time::adaptive::adjoint::{AdjointError, OdeVjp};
#[cfg(test)]
use fs_time::{AdaptiveState, adaptive::adjoint::trajectory::RecordedRk45};
use fs_time::adaptive::adjoint::trajectory::{
    RecordingConfig, RecordingStatus, ReplayBudget, TrajectoryError,
};

/// Accepted-state composition with the existing fallible L-BFGS engine.
pub mod study;
/// Shared model-parameter and state estimation in one numerical study.
pub mod joint;
/// Shared interval-map boundary for explicit and implicit discrete adjoints.
pub mod intervals;
use intervals::{IntervalScheme, IntervalTape};

#[derive(Debug, Clone)]
pub struct IntervalPolicy {
    /// The end field is replaced by each declared knot endpoint. The other
    /// fields retain their production RK45 meanings.
    pub recording: RecordingConfig,
    pub initial_step: f64,
    pub max_attempts: usize,
    pub max_records: usize,
    pub replay: ReplayBudget,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WindowError {
    Invalid(&'static str),
    Allocation,
    NonFinite(&'static str),
    WorkspaceLimit { required: usize, limit: usize },
    EvaluationLimit,
    IntervalLimit,
    Observation(String),
    Model(String),
    ForwardStopped { interval: usize, status: RecordingStatus },
    Trajectory { interval: usize, source: TrajectoryError },
    /// A backend refused a particular interval; its original diagnostic is retained.
    Integrator { interval: usize, phase: &'static str, diagnostic: String },
    /// A completed map or pullback violated its declared shape/time contract.
    IntervalOutput { interval: usize, what: &'static str },
    Cancelled,
}
impl std::fmt::Display for WindowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "weak-constraint window: {self:?}")
    }
}
impl std::error::Error for WindowError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Trajectory { source, .. } => Some(source), _ => None }
    }
}
fn trajectory_error(interval: usize, source: TrajectoryError) -> WindowError {
    if source == TrajectoryError::Step(AdjointError::Cancelled) { WindowError::Cancelled }
    else { WindowError::Trajectory { interval, source } }
}
fn poll(cancelled: &mut dyn FnMut() -> bool) -> Result<(), WindowError> {
    if cancelled() { Err(WindowError::Cancelled) } else { Ok(()) }
}
fn finite(value: f64, stage: &'static str) -> Result<f64, WindowError> {
    if value.is_finite() { Ok(value) } else { Err(WindowError::NonFinite(stage)) }
}
fn zeros(n: usize) -> Result<Vec<f64>, WindowError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| WindowError::Allocation)?;
    out.resize(n, 0.0); Ok(out)
}
fn copy(values: &[f64]) -> Result<Vec<f64>, WindowError> {
    let mut out = zeros(values.len())?; out.copy_from_slice(values); Ok(out)
}
#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, x: f64) -> Result<(), WindowError> {
        let next = finite(self.value + x, "objective sum")?;
        self.correction += if self.value.abs() >= x.abs() {
            (self.value-next)+x
        } else { (x-next)+self.value };
        self.value = next; finite(self.correction, "objective compensation")?; Ok(())
    }
    fn finish(&self) -> Result<f64, WindowError> { finite(self.value+self.correction, "objective total") }
}

/// Cumulative allowances survive failed initialization and failed optimizer
/// trials. Not Clone: checkpointing a numerical study never refunds callbacks.
/// Each interval has the separate recording/replay caps in IntervalPolicy.
/// Memory counts f64 components, excludes fixed window data, the borrowed point,
/// optimizer storage, model-owned scratch, and production interval-tape memory.
#[derive(Debug)]
pub struct WindowControl {
    max_evaluations: usize,
    max_intervals: usize,
    workspace: usize,
    evaluations: usize,
    intervals: usize,
}
impl WindowControl {
    pub fn new(max_evaluations: usize, max_intervals: usize, workspace: usize) -> Self {
        Self { max_evaluations, max_intervals, workspace, evaluations: 0, intervals: 0 }
    }
    pub fn evaluations(&self) -> usize { self.evaluations }
    pub fn interval_attempts(&self) -> usize { self.intervals }
    pub fn extend(&mut self, max_evaluations: usize, max_intervals: usize, workspace: usize)
        -> Result<(), WindowError>
    {
        if max_evaluations < self.max_evaluations || max_intervals < self.max_intervals || workspace < self.workspace {
            return Err(WindowError::Invalid("allowances may only increase"));
        }
        self.max_evaluations = max_evaluations; self.max_intervals = max_intervals; self.workspace = workspace;
        Ok(())
    }
    fn admit(&mut self, intervals: usize, workspace: usize) -> Result<(), WindowError> {
        if workspace > self.workspace { return Err(WindowError::WorkspaceLimit { required: workspace, limit: self.workspace }); }
        if self.evaluations == self.max_evaluations { return Err(WindowError::EvaluationLimit); }
        if intervals > self.max_intervals-self.intervals { return Err(WindowError::IntervalLimit); }
        self.evaluations += 1; Ok(())
    }
}

/// Joint scalar observation loss on the entire knot-major physical state.
/// Overwrite every state partial, including zeros at unobserved knots. Times
/// and model parameters are fixed. Callbacks must remain pure and unchanged
/// across optimizer continuation; their own scratch/work needs caller limits.
pub trait WindowObjective {
    fn evaluate(&self, times: &[f64], dimension: usize, states: &[f64],
        state_bar: &mut [f64], cancelled: &mut dyn FnMut() -> bool) -> Result<f64, String>;

    /// Explicit observation partials in the OdeVjp model's parameter coordinates,
    /// holding every knot state fixed. Override for parameter-dependent sensors.
    /// The default declares parameter-independent observation loss. These do not
    /// include state-prior or process-noise derivatives: those scales are fixed.
    fn parameter_partials(&self, _times: &[f64], _dimension: usize, _states: &[f64],
        parameter_bar: &mut [f64], _cancelled: &mut dyn FnMut() -> bool) -> Result<(), String>
    {
        parameter_bar.fill(0.0); Ok(())
    }
}

/// Immutable knot layout and diagonal background/model-error declarations.
/// Decisions z are dimensionless: x[k,i] = reference[k,i] + scale[i]*z[k,i].
/// The initial reference is the background mean. Other reference states ONLY
/// center optimizer coordinates; they add no prior penalty or pseudo-data.
#[derive(Debug, Clone)]
pub struct WeakConstraintWindow {
    times: Vec<f64>,
    reference: Vec<f64>,
    scale: Vec<f64>,
    background_sigma: Vec<f64>,
    model_sigma: Vec<f64>,
}
impl WeakConstraintWindow {
    /// Model sigma is interval-major, one positive finite scale per state per
    /// interval. All scales use that component's physical units. Singular or
    /// hard model constraints are not approximated by silently added jitter.
    /// max_components bounds all owned f64 entries BEFORE allocation.
    pub fn new(times: &[f64], reference: &[f64], scale: &[f64], background_sigma: &[f64],
        model_sigma: &[f64], max_components: usize) -> Result<Self, WindowError>
    {
        let n = scale.len();
        let length = n.checked_mul(times.len()).ok_or(WindowError::Invalid("window extent overflow"))?;
        let required = length.checked_mul(2).and_then(|v| v.checked_add(n)).and_then(|v| v.checked_add(times.len()))
            .ok_or(WindowError::Invalid("window storage overflow"))?;
        if required > max_components { return Err(WindowError::WorkspaceLimit { required, limit: max_components }); }
        if n == 0 || times.len() < 2 || reference.len() != length || background_sigma.len() != n
            || model_sigma.len() != length-n || times.iter().any(|t| !t.is_finite())
            || times.windows(2).any(|t| t[1] <= t[0] || !(t[1]-t[0]).is_finite())
            || reference.iter().any(|x| !x.is_finite())
            || scale.iter().chain(background_sigma).chain(model_sigma).any(|x| !x.is_finite() || *x <= 0.0)
        { return Err(WindowError::Invalid("ordered times, finite states and positive shape-matched scales required")); }
        Ok(Self { times: copy(times)?, reference: copy(reference)?, scale: copy(scale)?,
            background_sigma: copy(background_sigma)?, model_sigma: copy(model_sigma)? })
    }
    pub fn times(&self) -> &[f64] { &self.times }
    pub fn dimension(&self) -> usize { self.scale.len() }
    pub fn control_dimension(&self) -> usize { self.reference.len() }
    /// Includes point/states/gradient/defect outputs and all non-tape scratch.
    pub fn workspace_components(&self, parameters: usize) -> Result<usize, WindowError> {
        self.reference.len().checked_mul(4)
            .and_then(|v| self.dimension().checked_mul(4).and_then(|n| v.checked_add(n)))
            .and_then(|v| parameters.checked_mul(3).and_then(|p| v.checked_add(p)))
            .ok_or(WindowError::Invalid("evaluation workspace overflow"))
    }

    /// One complete numerical objective and gradient. Every interval uses its
    /// own accepted RK45 map and transposed initial-state derivative. Return
    /// parameter partials too, holding knots and all covariance scales fixed;
    /// the state-only study does not optimize these parameters.
    /// The window is immutable and no partial objective/gradient escapes errors.
    /// Successful earlier INTERVAL work is charged but is recomputed on retry.
    pub fn evaluate<M: OdeVjp, O: WindowObjective, C: FnMut() -> bool>(
        &self, model: &M, objective: &O, point: &[f64], policy: &IntervalPolicy,
        control: &mut WindowControl, cancelled: &mut C,
    ) -> Result<WindowEvaluation, WindowError> {
        self.evaluate_using(model, objective, point, policy, control, cancelled)
    }

    /// Evaluate with a declared discrete interval scheme. The same observation,
    /// background, defect and coordinate-chain-rule code serves every scheme.
    /// Only a complete endpoint and its matching transpose action are usable.
    /// Numerical step choices and model/solver policy remain fixed for a gradient.
    pub fn evaluate_using<M, O: WindowObjective, S: IntervalScheme<M>, C: FnMut() -> bool>(
        &self, model: &M, objective: &O, point: &[f64], policy: &S,
        control: &mut WindowControl, cancelled: &mut C,
    ) -> Result<WindowEvaluation, WindowError> {
        poll(cancelled)?;
        if policy.dimension(model) != self.dimension() || point.len() != self.control_dimension()
            || point.iter().any(|x| !x.is_finite())
        { return Err(WindowError::Invalid("model/point shape or values")); }
        policy.validate(&self.times)?;
        control.admit(self.times.len()-1, self.workspace_components(policy.parameter_count(model))?)?;
        self.evaluate_admitted_using(model, objective, point, policy, control, cancelled)
    }

    // The joint parameter path admits and charges the whole trial BEFORE its
    // model factory runs, then enters here without charging a second evaluation.
    #[allow(clippy::too_many_arguments)]
    fn evaluate_admitted<M: OdeVjp, O: WindowObjective, C: FnMut() -> bool>(
        &self, model: &M, objective: &O, point: &[f64], policy: &IntervalPolicy,
        control: &mut WindowControl, cancelled: &mut C,
    ) -> Result<WindowEvaluation, WindowError> {
        self.evaluate_admitted_using(model, objective, point, policy, control, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_admitted_using<M, O: WindowObjective, S: IntervalScheme<M>, C: FnMut() -> bool>(
        &self, model: &M, objective: &O, point: &[f64], policy: &S,
        control: &mut WindowControl, cancelled: &mut C,
    ) -> Result<WindowEvaluation, WindowError> {
        let n = self.dimension(); let length = self.control_dimension();
        let steps = self.times.len()-1; let p = policy.parameter_count(model);
        let mut states = zeros(length)?;
        for (i, (x, z)) in states.iter_mut().zip(point).enumerate() {
            if i % 256 == 0 { poll(cancelled)?; }
            *x = finite(self.scale[i%n].mul_add(*z, self.reference[i]), "physical decision state")?;
        }
        let mut gradient = zeros(length)?; gradient.fill(f64::NAN);
        let mut stopped = false;
        let mut check = || { stopped |= cancelled(); stopped };
        let observed = objective.evaluate(&self.times, n, &states, &mut gradient, &mut check);
        if check() { return Err(WindowError::Cancelled); }
        let observation_value = finite(observed.map_err(WindowError::Observation)?, "observation objective")?;
        for chunk in gradient.chunks(256) {
            poll(cancelled)?;
            if chunk.iter().any(|g| !g.is_finite()) { return Err(WindowError::NonFinite("observation partials")); }
        }
        let mut parameter_gradient = zeros(p)?;
        if p != 0 {
            parameter_gradient.fill(f64::NAN);
            let mut stopped = false;
            let mut check = || { stopped |= cancelled(); stopped };
            let result = objective.parameter_partials(&self.times, n, &states, &mut parameter_gradient, &mut check);
            if check() { return Err(WindowError::Cancelled); }
            result.map_err(WindowError::Observation)?;
            for chunk in parameter_gradient.chunks(256) {
                poll(cancelled)?;
                if chunk.iter().any(|g| !g.is_finite()) { return Err(WindowError::NonFinite("observation parameter partials")); }
            }
        }
        let mut background = Sum::default(); let mut model_error = Sum::default();
        for i in 0..n {
            if i % 256 == 0 { poll(cancelled)?; }
            let r = finite((states[i]-self.reference[i])/self.background_sigma[i], "background residual")?;
            background.add((0.5*r)*r)?;
            gradient[i] = finite(gradient[i] + r/self.background_sigma[i], "background gradient")?;
        }
        let mut defects = zeros(steps*n)?;
        let mut seed = zeros(n)?; let direct = zeros(p)?;
        let (mut accepted_steps, mut replayed_steps) = (0usize, 0usize);
        for k in 0..steps {
            poll(cancelled)?;
            if policy.dimension(model) != n || policy.parameter_count(model) != p { return Err(WindowError::Invalid("model dimensions changed")); }
            // Admission above reserves enough remaining allowance for the whole
            // trial, but only intervals actually attempted are spent.
            control.intervals += 1;
            let tape = policy.record(model, k, self.times[k], self.times[k+1],
                &states[k*n..(k+1)*n], &mut || cancelled())?;
            poll(cancelled)?;
            if tape.end_time() != self.times[k+1] || tape.endpoint().len() != n {
                return Err(WindowError::IntervalOutput { interval: k, what: "endpoint time or dimension" });
            }
            for i in 0..n {
                if i % 256 == 0 { poll(cancelled)?; }
                let offset = k*n+i; let next = (k+1)*n+i;
                let defect = finite(states[next]-tape.endpoint()[i], "model defect")?;
                defects[offset] = defect;
                let r = finite(defect/self.model_sigma[offset], "scaled model defect")?;
                model_error.add((0.5*r)*r)?;
                let partial = finite(r/self.model_sigma[offset], "model defect partial")?;
                gradient[next] = finite(gradient[next]+partial, "right endpoint gradient")?;
                seed[i] = -partial;
            }
            let bar = tape.pullback(&seed, &direct, &mut || cancelled())?;
            poll(cancelled)?;
            if bar.initial.len() != n || bar.parameters.len() != p {
                return Err(WindowError::IntervalOutput { interval: k, what: "pullback dimension" });
            }
            for (i, update) in bar.initial.iter().enumerate() {
                if i % 256 == 0 { poll(cancelled)?; }
                gradient[k*n+i] = finite(gradient[k*n+i]+update, "left endpoint gradient")?;
            }
            for (i, update) in bar.parameters.iter().enumerate() {
                if i % 256 == 0 { poll(cancelled)?; }
                parameter_gradient[i] = finite(parameter_gradient[i]+update, "model parameter gradient")?;
            }
            accepted_steps = accepted_steps.checked_add(tape.accepted_steps()).ok_or(WindowError::Invalid("step count overflow"))?;
            replayed_steps = replayed_steps.checked_add(bar.replayed_steps).ok_or(WindowError::Invalid("replay count overflow"))?;
        }
        for (i, g) in gradient.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cancelled)?; }
            *g = finite(*g*self.scale[i%n], "decision-coordinate gradient")?;
        }
        let background_value = background.finish()?; let model_value = model_error.finish()?;
        let mut total = Sum::default(); total.add(observation_value)?; total.add(background_value)?; total.add(model_value)?;
        let controls = copy(point)?;
        poll(cancelled)?;
        Ok(WindowEvaluation { controls, states, gradient, parameter_gradient, defects, value: total.finish()?,
            observation_value, background_value, model_value, accepted_steps, replayed_steps })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct WindowEvaluation {
    pub controls: Vec<f64>,
    pub states: Vec<f64>,
    /// Derivative with respect to dimensionless controls, NOT physical states.
    pub gradient: Vec<f64>,
    /// Partial derivative in OdeVjp parameter coordinates, with knots fixed.
    /// Includes explicit observation partials and every interval's dynamics.
    pub parameter_gradient: Vec<f64>,
    /// Actual knot-major endpoint increments, `x[k+1]-M[k](x[k])`.
    pub defects: Vec<f64>,
    pub value: f64,
    pub observation_value: f64,
    pub background_value: f64,
    pub model_value: f64,
    pub accepted_steps: usize,
    pub replayed_steps: usize,
}

#[cfg(test)]
mod tests;

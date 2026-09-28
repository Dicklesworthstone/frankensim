//! Joint estimation of static model coordinates and weak-constraint knot states.
//!
//! One parameter vector is shared by every interval, not an independently
//! drifting parameter at every knot. A fixed diagonal parameter prior declares
//! the tradeoff with freely estimated model-error increments. Parameter/state
//! confounding is not solved by adding a prior and no identifiability is claimed.
//! Model coordinates can be physical parameters or explicit transforms (e.g.
//! log rate); OdeVjp and observation partials must use those SAME coordinates.

use super::{IntervalPolicy, WeakConstraintWindow, WindowControl, WindowError,
    WindowEvaluation, WindowObjective, Sum, copy, finite, poll, zeros};
use super::study::{StudySettings, WindowStudyError};
use fs_time::adaptive::adjoint::OdeVjp;
use crate::{LbfgsError, LbfgsReport, LbfgsState, StopReason, StopRule};

/// Pure factory at one parameter point. The returned model owns the observation
/// loss and its explicit parameter partials as well as the ODE and its VJP.
/// Parameters are unbounded finite coordinates; implement a smooth transform
/// inside the model for positive physical quantities. Domain/producer errors
/// propagate, never become artificial objective penalties. The factory receives
/// cancellation and is charged as part of its window evaluation before running.
pub trait ParameterFamily {
    type Model: OdeVjp + WindowObjective;
    fn instantiate(&self, coordinates: &[f64], cancelled: &mut dyn FnMut() -> bool)
        -> Result<Self::Model, String>;
}

/// Control layout: all existing knot-state controls, then p shared parameter
/// controls. theta[j] = mean[j] + scale[j]*z[j]; sigma[j] describes the fixed
/// prior in theta coordinates. No process/background covariance is inferred.
pub struct JointWindow<'a, F> {
    window: &'a WeakConstraintWindow,
    family: &'a F,
    mean: Vec<f64>,
    scale: Vec<f64>,
    sigma: Vec<f64>,
}
impl<'a, F: ParameterFamily> JointWindow<'a, F> {
    pub fn new(window: &'a WeakConstraintWindow, family: &'a F,
        mean: &[f64], scale: &[f64], sigma: &[f64], max_components: usize,
    ) -> Result<Self, WindowError> {
        let required = mean.len().checked_mul(3).ok_or(WindowError::Invalid("parameter prior extent"))?;
        if required > max_components {
            return Err(WindowError::WorkspaceLimit { required, limit: max_components });
        }
        if mean.is_empty() || scale.len() != mean.len() || sigma.len() != mean.len()
            || mean.iter().any(|v| !v.is_finite())
            || scale.iter().chain(sigma).any(|v| !v.is_finite() || *v <= 0.0)
        { return Err(WindowError::Invalid("finite parameter means and positive matching scales required")); }
        window.control_dimension().checked_add(mean.len()).ok_or(WindowError::Invalid("joint control extent"))?;
        Ok(Self { window, family, mean: copy(mean)?, scale: copy(scale)?, sigma: copy(sigma)? })
    }
    pub fn control_dimension(&self) -> usize { self.window.control_dimension()+self.mean.len() }
    pub fn parameter_count(&self) -> usize { self.mean.len() }
    /// Includes the complete state result plus joint controls/gradient/parameters;
    /// excludes fixed window/factory data, optimizer storage and interval tapes.
    pub fn workspace_components(&self) -> Result<usize, WindowError> {
        self.window.workspace_components(self.mean.len())?
            .checked_add(self.control_dimension().checked_mul(2).ok_or(WindowError::Invalid("joint workspace extent"))?)
            .and_then(|v| self.mean.len().checked_mul(3).and_then(|p| v.checked_add(p)))
            .ok_or(WindowError::Invalid("joint workspace extent"))
    }
    pub fn evaluate<C: FnMut() -> bool>(&self, point: &[f64], policy: &IntervalPolicy,
        control: &mut WindowControl, cancelled: &mut C) -> Result<JointEvaluation, WindowError>
    {
        poll(cancelled)?;
        if point.len() != self.control_dimension() || point.iter().any(|v| !v.is_finite())
            || !policy.initial_step.is_finite() || policy.initial_step <= 0.0
        { return Err(WindowError::Invalid("invalid joint controls or initial step")); }
        control.admit(self.window.times.len()-1, self.workspace_components()?)?;
        let n = self.window.control_dimension(); let p = self.mean.len();
        let mut parameters = zeros(p)?;
        for j in 0..p {
            if j % 256 == 0 { poll(cancelled)?; }
            parameters[j] = finite(self.scale[j].mul_add(point[n+j], self.mean[j]), "parameter coordinate")?;
        }
        let mut stopped = false;
        let mut check = || { stopped |= cancelled(); stopped };
        let model = self.family.instantiate(&parameters, &mut check);
        if check() { return Err(WindowError::Cancelled); }
        let model = model.map_err(WindowError::Model)?;
        if model.dimension() != self.window.dimension() || model.parameter_count() != p {
            return Err(WindowError::Invalid("instantiated joint model dimensions"));
        }
        let window = self.window.evaluate_admitted(&model, &model, &point[..n], policy, control, cancelled)?;
        let mut gradient = zeros(self.control_dimension())?;
        gradient[..n].copy_from_slice(&window.gradient);
        let mut prior = Sum::default();
        for j in 0..p {
            if j % 256 == 0 { poll(cancelled)?; }
            let r = finite((parameters[j]-self.mean[j])/self.sigma[j], "parameter prior residual")?;
            prior.add((0.5*r)*r)?;
            let partial = finite(window.parameter_gradient[j]+r/self.sigma[j], "joint parameter partial")?;
            gradient[n+j] = finite(partial*self.scale[j], "scaled joint parameter partial")?;
        }
        let parameter_penalty = prior.finish()?;
        let mut value = Sum::default(); value.add(window.value)?; value.add(parameter_penalty)?;
        let controls = copy(point)?;
        poll(cancelled)?;
        Ok(JointEvaluation { controls, parameters, gradient, value: value.finish()?, parameter_penalty, window })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct JointEvaluation {
    pub controls: Vec<f64>,
    /// Coordinates passed to the factory, not implicitly transformed quantities.
    pub parameters: Vec<f64>,
    pub gradient: Vec<f64>,
    pub value: f64,
    pub parameter_penalty: f64,
    /// State/observation/process decomposition excludes the parameter prior.
    pub window: WindowEvaluation,
}

/// The same fallible L-BFGS engine used by the state-only study. Accepted model
/// coordinates, states, gradients and penalty components remain paired. Clone
/// retains curvature and evaluation count but cannot clone the work allowance.
/// Factories/data/scales/policies must remain unchanged across continuation.
pub struct JointWindowStudy<'a, 'w, F> {
    window: &'a JointWindow<'w, F>,
    policy: IntervalPolicy,
    settings: StudySettings,
    optimizer: LbfgsState,
    accepted: JointEvaluation,
}
impl<F> Clone for JointWindowStudy<'_, '_, F> {
    fn clone(&self) -> Self {
        Self { window: self.window, policy: self.policy.clone(), settings: self.settings,
            optimizer: self.optimizer.clone(), accepted: self.accepted.clone() }
    }
}
impl<'a, 'w, F: ParameterFamily> JointWindowStudy<'a, 'w, F> {
    pub fn new(window: &'a JointWindow<'w, F>, point: &[f64], policy: IntervalPolicy,
        settings: StudySettings, control: &mut WindowControl, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindowStudyError> {
        poll(cancelled).map_err(LbfgsError::Evaluation)?;
        settings.validate(window.control_dimension())?;
        let mut accepted = None;
        let optimizer = LbfgsState::try_new(point, settings.memory, &mut |z| {
            let evaluation = window.evaluate(z, &policy, control, cancelled)?;
            let result = (evaluation.value, evaluation.gradient.clone()); accepted = Some(evaluation);
            Ok::<_, WindowError>(result)
        })?;
        let accepted = accepted.ok_or(LbfgsError::InvalidInput("missing joint initial evaluation"))?;
        Ok(Self { window, policy, settings, optimizer, accepted })
    }
    pub fn optimizer(&self) -> &LbfgsState { &self.optimizer }
    pub fn accepted(&self) -> &JointEvaluation { &self.accepted }
    pub fn run(&mut self, additional_iterations: usize, control: &mut WindowControl,
        cancelled: &mut impl FnMut() -> bool) -> Result<LbfgsReport, WindowStudyError>
    {
        poll(cancelled).map_err(LbfgsError::Evaluation)?;
        let mut result = self.advance(0, control, cancelled)?;
        for _ in 0..additional_iterations {
            if result.reason != StopReason::IterationCap { break; }
            poll(cancelled).map_err(LbfgsError::Evaluation)?;
            result = self.advance(1, control, cancelled)?;
        }
        Ok(result)
    }
    fn advance(&mut self, iterations: usize, control: &mut WindowControl,
        cancelled: &mut impl FnMut() -> bool) -> Result<LbfgsReport, WindowStudyError>
    {
        let window = self.window; let policy = &self.policy; let mut candidate = None;
        let result = self.optimizer.try_run(&mut |z| {
            let evaluation = window.evaluate(z, policy, control, cancelled)?;
            let output = (evaluation.value, evaluation.gradient.clone()); candidate = Some(evaluation);
            Ok::<_, WindowError>(output)
        }, &StopRule::GradNorm(self.settings.gradient_tolerance), iterations, self.settings.max_evaluations);
        if let Some(candidate) = candidate {
            if candidate.controls == self.optimizer.x { self.accepted = candidate; }
        }
        result
    }
}

#[cfg(test)]
mod tests;

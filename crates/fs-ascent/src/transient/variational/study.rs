//! The existing fallible L-BFGS engine driving a weak-constraint window.

use super::{IntervalPolicy, WeakConstraintWindow, WindowControl, WindowError,
    WindowEvaluation, WindowObjective, poll};
use crate::{LbfgsError, LbfgsReport, LbfgsState, StopReason, StopRule};
use super::intervals::IntervalScheme;
#[cfg(test)]
use fs_time::adaptive::adjoint::OdeVjp;

pub type WindowStudyError = LbfgsError<WindowError>;

#[derive(Debug, Clone, Copy)]
pub struct StudySettings {
    pub memory: usize,
    pub gradient_tolerance: f64,
    /// Absolute cumulative optimizer evaluation ceiling, including failures.
    pub max_evaluations: usize,
    /// Conservative L-BFGS/vector/accepted-result envelope:
    /// (2*memory+32)*controls + 4*memory + max_evaluations scalars.
    /// Fixed problem data, producer scratch and interval-tape memory are separate.
    pub max_optimizer_components: usize,
}
impl StudySettings {
    pub(super) fn validate(self, dimension: usize) -> Result<(), WindowStudyError> {
        if self.memory == 0 || self.max_evaluations == 0 || dimension == 0
            || !self.gradient_tolerance.is_finite() || self.gradient_tolerance <= 0.0
        { return Err(LbfgsError::InvalidInput("positive memory, tolerance, dimension and evaluation cap required")); }
        let required = self.memory.checked_mul(2).and_then(|v| v.checked_add(32))
            .and_then(|v| v.checked_mul(dimension))
            .and_then(|v| self.memory.checked_mul(4).and_then(|m| v.checked_add(m)))
            .and_then(|v| v.checked_add(self.max_evaluations))
            .ok_or(LbfgsError::InvalidInput("optimizer workspace extent overflow"))?;
        if required > self.max_optimizer_components {
            return Err(LbfgsError::Evaluation(WindowError::WorkspaceLimit { required, limit: self.max_optimizer_components }));
        }
        Ok(())
    }
}

/// Coherent accepted optimizer and trajectory. No mutable optimizer or candidate
/// result escapes. One accepted iteration is committed at a time; errors retain
/// earlier accepted iterations and their matching objective/defect decomposition.
/// A cancelled trial never replaces them. Bounded L-BFGS algebra itself has no
/// internal cancellation checkpoints; cancellation is polled around evaluations
/// and at accepted-iteration boundaries, not with a claimed latency bound.
///
/// Clone checkpoints numerical work, not the separately borrowed WindowControl.
/// The model, objective, times, scales and numerical policy must remain unchanged
/// on continuation. GradNorm means local numerical stationarity, not a posterior
/// covariance, exact continuous sensitivity, or certified/global optimality.
pub struct WeakConstraintStudy<'a, M, O, S = IntervalPolicy> {
    window: &'a WeakConstraintWindow,
    model: &'a M,
    objective: &'a O,
    policy: S,
    settings: StudySettings,
    optimizer: LbfgsState,
    accepted: WindowEvaluation,
}
impl<M, O, S: Clone> Clone for WeakConstraintStudy<'_, M, O, S> {
    fn clone(&self) -> Self {
        Self { window: self.window, model: self.model, objective: self.objective,
            policy: self.policy.clone(), settings: self.settings,
            optimizer: self.optimizer.clone(), accepted: self.accepted.clone() }
    }
}
impl<'a, M, O: WindowObjective, S: IntervalScheme<M>> WeakConstraintStudy<'a, M, O, S> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(window: &'a WeakConstraintWindow, model: &'a M, objective: &'a O,
        point: &[f64], policy: S, settings: StudySettings,
        control: &mut WindowControl, cancelled: &mut impl FnMut()->bool,
    ) -> Result<Self, WindowStudyError> {
        poll(cancelled).map_err(LbfgsError::Evaluation)?;
        settings.validate(window.control_dimension())?;
        let mut accepted = None;
        let optimizer = LbfgsState::try_new(point, settings.memory, &mut |z| {
            let result = window.evaluate_using(model, objective, z, &policy, control, cancelled)?;
            let response = (result.value, result.gradient.clone()); accepted = Some(result); Ok::<_,WindowError>(response)
        })?;
        let accepted = accepted.ok_or(LbfgsError::InvalidInput("missing initial trajectory"))?;
        Ok(Self { window, model, objective, policy, settings, optimizer, accepted })
    }
    pub fn optimizer(&self) -> &LbfgsState { &self.optimizer }
    pub fn accepted(&self) -> &WindowEvaluation { &self.accepted }

    /// Additional accepted iterations; both optimizer and external window work
    /// allowances are cumulative. An interrupted line search restarts on retry,
    /// with all failed window evaluations/intervals still charged.
    pub fn run(&mut self, additional_iterations: usize, control: &mut WindowControl,
        cancelled: &mut impl FnMut()->bool) -> Result<LbfgsReport, WindowStudyError>
    {
        poll(cancelled).map_err(LbfgsError::Evaluation)?;
        let mut report = self.advance(0, control, cancelled)?;
        for _ in 0..additional_iterations {
            if report.reason != StopReason::IterationCap { break; }
            poll(cancelled).map_err(LbfgsError::Evaluation)?;
            report = self.advance(1, control, cancelled)?;
        }
        Ok(report)
    }
    fn advance(&mut self, iterations: usize, control: &mut WindowControl,
        cancelled: &mut impl FnMut()->bool) -> Result<LbfgsReport, WindowStudyError>
    {
        let (window, model, objective, policy) = (self.window, self.model, self.objective, &self.policy);
        let mut candidate = None;
        let result = self.optimizer.try_run(&mut |z| {
            let evaluation = window.evaluate_using(model, objective, z, policy, control, cancelled)?;
            let response = (evaluation.value, evaluation.gradient.clone()); candidate = Some(evaluation);
            Ok::<_, WindowError>(response)
        }, &StopRule::GradNorm(self.settings.gradient_tolerance), iterations, self.settings.max_evaluations);
        // At most one accepted step in this call. A complete but rejected trial
        // cannot overwrite the accepted field, even if a later callback failed.
        if let Some(candidate) = candidate {
            if candidate.controls == self.optimizer.x { self.accepted = candidate; }
        }
        result
    }
}

#[cfg(test)]
mod tests;

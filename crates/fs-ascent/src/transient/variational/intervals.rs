//! A complete discrete forecast and its matching adjoint for one window interval.
//!
//! This boundary lets physical integrators participate without impersonating an
//! ODE RHS or duplicating the variational objective/optimizer. Implementations
//! own time stepping, solver limits, replay checks and tape memory admission.
//! A tape must stay bound to the same immutable model and numerical map. The
//! outer window checks endpoint times, dimensions and finite accumulations;
//! those checks alone cannot establish that a supplied adjoint is correct.

use super::{IntervalPolicy, WindowError, poll, trajectory_error};
use fs_time::AdaptiveState;
use fs_time::adaptive::adjoint::OdeVjp;
use fs_time::adaptive::adjoint::trajectory::{
    RecordedRk45, RecordingStatus, TrajectoryGradient,
};

/// Checkpointed fine solver steps inside unchanged reconstruction knots.
pub mod substeps;

/// A completed interval map. No partial forward state may be exposed as a tape.
/// The transpose action differentiates this exact map, not the solver stopping
/// decisions. Derivatives of start/end times and numerical policy are excluded.
pub trait IntervalTape {
    fn endpoint(&self) -> &[f64];
    fn end_time(&self) -> f64;
    fn accepted_steps(&self) -> usize;
    /// Pull back the endpoint seed, adding explicit fixed-state parameter terms.
    /// Return n initial-state and p parameter partials, never partial gradients.
    fn pullback(&self, seed: &[f64], direct_parameters: &[f64],
        cancelled: &mut dyn FnMut() -> bool) -> Result<TrajectoryGradient, WindowError>;
}

/// Numerical policy for a fixed model type. GAT tapes borrow their model/policy,
/// preventing the window from replacing either during its forward/reverse pair.
/// A successful record owns or borrows everything needed to reverse that map.
/// Model dimensions and all scientific inputs must remain unchanged on retry.
pub trait IntervalScheme<M> {
    type Tape<'a>: IntervalTape where Self: 'a, M: 'a;
    fn dimension(&self, model: &M) -> usize;
    fn parameter_count(&self, model: &M) -> usize;
    /// Reject unsupported time grids or policies before trial admission/work.
    fn validate(&self, times: &[f64]) -> Result<(), WindowError>;
    fn record<'a>(&'a self, model: &'a M, interval: usize, start: f64, end: f64,
        initial: &[f64], cancelled: &mut dyn FnMut() -> bool)
        -> Result<Self::Tape<'a>, WindowError>;
}

/// Adapter for the existing RK45 recording. No numerical stages or reductions
/// are copied: the production recording and pullback perform all the work.
pub struct Rk45Interval<'a, M> {
    recording: RecordedRk45<'a, M>,
    policy: &'a IntervalPolicy,
    interval: usize,
}
impl<M: OdeVjp> IntervalScheme<M> for IntervalPolicy {
    type Tape<'a> = Rk45Interval<'a, M> where Self: 'a, M: 'a;
    fn dimension(&self, model: &M) -> usize { model.dimension() }
    fn parameter_count(&self, model: &M) -> usize { model.parameter_count() }
    fn validate(&self, _times: &[f64]) -> Result<(), WindowError> {
        if !self.initial_step.is_finite() || self.initial_step <= 0.0 {
            return Err(WindowError::Invalid("positive finite initial RK45 step required"));
        }
        Ok(())
    }
    fn record<'a>(&'a self, model: &'a M, interval: usize, start: f64, end: f64,
        initial: &[f64], cancelled: &mut dyn FnMut() -> bool)
        -> Result<Self::Tape<'a>, WindowError>
    {
        poll(cancelled)?;
        let mut config = self.recording.clone(); config.end = end;
        let initial = AdaptiveState::new(start, initial, self.initial_step);
        let mut recording = RecordedRk45::new(model, initial, config)
            .map_err(|e| trajectory_error(interval, e))?;
        let report = recording.advance(self.max_attempts, self.max_records, &mut || cancelled())
            .map_err(|e| trajectory_error(interval, e))?;
        if report.status == RecordingStatus::Cancelled { return Err(WindowError::Cancelled); }
        if report.status != RecordingStatus::ReachedEnd {
            return Err(WindowError::ForwardStopped { interval, status: report.status });
        }
        Ok(Rk45Interval { recording, policy: self, interval })
    }
}
impl<M: OdeVjp> IntervalTape for Rk45Interval<'_, M> {
    fn endpoint(&self) -> &[f64] { &self.recording.state().u }
    fn end_time(&self) -> f64 { self.recording.state().t }
    fn accepted_steps(&self) -> usize { self.recording.accepted_steps() }
    fn pullback(&self, seed: &[f64], direct_parameters: &[f64],
        cancelled: &mut dyn FnMut() -> bool) -> Result<TrajectoryGradient, WindowError>
    {
        self.recording.pullback(seed, direct_parameters, self.policy.replay, &mut || cancelled())
            .map_err(|e| trajectory_error(self.interval, e))
    }
}

#[cfg(test)]
mod tests;

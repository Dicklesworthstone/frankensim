//! Time-series objectives on one observation-aligned adaptive trajectory.
//!
//! Requested times become exact accepted endpoints of the production RK45
//! stepper. No linear/dense interpolation, nearest-step substitution, separate
//! solve per observation, or state-sized forward tangent per parameter is used.
//! The resulting discrete gradients hold this accepted mesh fixed. Adding an
//! observation can change the forward mesh; it is not a passive dense-output
//! query and is not differentiation of the controller or sampling times.

use super::*;

/// An additive scalar objective at one fixed observation time. `sample` indexes
/// the original immutable time array, including repeated times. Write every
/// state and direct-parameter partial, including zeros, holding time fixed.
/// Callbacks must be pure and deterministic: a failed reverse sweep can retry
/// them. They are called once per sample, in reverse declaration order, only
/// after the corresponding forward state has passed its replay check.
pub trait SampleObjective {
    fn evaluate(
        &self, sample: usize, time: f64, state: &[f64],
        state_bar: &mut [f64], parameter_bar: &mut [f64],
    ) -> Result<f64, String>;
}

#[derive(Debug, Clone, PartialEq)]
pub struct SampledGradient {
    /// Sum of the caller's observation objectives, not a confidence bound.
    pub value: f64,
    pub gradient: TrajectoryGradient,
    pub observations: usize,
}

impl<'a, M: OdeVjp> RecordedRk45<'a, M> {
    /// Retain a nondecreasing, finite, explicit observation timetable inside
    /// the integration interval. The cap is checked before copying/scanning.
    /// Repeated times are allowed for distinct sensors; their statistical
    /// dependence and weighting remain the objective's responsibility.
    /// An empty timetable gives exactly the ordinary recording path.
    pub fn new_sampled(
        model: &'a M, state: AdaptiveState, config: RecordingConfig,
        times: &[f64], max_samples: usize,
    ) -> Result<Self, TrajectoryError> {
        if times.len() > max_samples {
            return Err(AdjointError::InvalidInput("observation timetable exceeds sample cap").into());
        }
        if times.iter().any(|t| !t.is_finite() || *t < state.t || *t > config.end)
            || times.windows(2).any(|pair| pair[0] > pair[1])
        {
            return Err(AdjointError::InvalidInput("sample times must be ordered and inside the interval").into());
        }
        let mut recording = Self::new(model, state, config)?;
        recording.sample_times.try_reserve_exact(times.len()).map_err(|_| TrajectoryError::Allocation)?;
        recording.sample_times.extend_from_slice(times);
        Ok(recording)
    }

    pub fn sample_times(&self) -> &[f64] { &self.sample_times }

    /// Accumulate all observation partials in ONE checkpointed reverse sweep.
    /// No observed full states or Jacobian rows are retained. Per-boundary
    /// callback scratch is n+p and is dropped before the step VJP, within the
    /// existing 17n+2p workspace cap (the live sweep cotangent is separate).
    ///
    /// Initial-time terms act directly on the initial cotangent. Endpoint terms
    /// act before reversing their step. Multiple readings at one time are each
    /// evaluated once. Initial-condition parameter dependence must be contracted
    /// with the returned initial cotangent by the caller, just as for `pullback`.
    /// Errors/cancellation return neither partial objective nor partial gradient.
    pub fn pullback_samples<O: SampleObjective, Cancel: FnMut() -> bool>(
        &self, objective: &O, budget: ReplayBudget, cancelled: &mut Cancel,
    ) -> Result<SampledGradient, TrajectoryError> {
        self.model_shape()?;
        if self.state.t != self.config.end { return Err(TrajectoryError::Incomplete); }
        if self.sample_times.is_empty() {
            return Err(AdjointError::InvalidInput("a sampled objective requires observations").into());
        }
        workspace_size(self.initial.len(), self.parameters, self.config.max_workspace_components)?;
        let required = fs_ad::revolve::min_budget(self.records.len());
        if budget.checkpoints < required {
            return Err(TrajectoryError::CheckpointLimit { required, limit: budget.checkpoints });
        }
        poll(cancelled)?;
        let mut progress = Progress { replays: 0, peak: 0, budget };
        let mut value = 0.0;
        let mut observations = 0;
        let bar = Cotangent { initial: vec![0.0; self.initial.len()], parameters: vec![0.0; self.parameters] };
        let mut observe = |step: usize, state: &[f64], bar: &mut Cotangent, cancelled: &mut Cancel| {
            let time = self.records[step].end;
            let first = self.sample_times.partition_point(|t| *t < time);
            let last = self.sample_times.partition_point(|t| *t <= time);
            self.observe_range(objective, first..last, state, bar, &mut value, &mut observations, cancelled)
        };
        let mut bar = if self.records.is_empty() { bar } else {
            self.segment(&self.initial, 0, self.records.len(), bar, 1, &mut progress, &mut observe, cancelled)?
        };
        let start = self.records.first().map_or(self.state.t, |r| r.start);
        let initial_samples = self.sample_times.partition_point(|t| *t <= start);
        self.observe_range(objective, 0..initial_samples, &self.initial, &mut bar,
            &mut value, &mut observations, cancelled)?;
        if observations != self.sample_times.len() {
            return Err(AdjointError::InvalidInput("an observation has no recorded endpoint").into());
        }
        poll(cancelled)?;
        Ok(SampledGradient { value, observations, gradient: TrajectoryGradient {
            initial: bar.initial, parameters: bar.parameters,
            replayed_steps: progress.replays, peak_checkpoints: progress.peak,
        } })
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_range<O: SampleObjective, Cancel: FnMut() -> bool>(
        &self, objective: &O, range: std::ops::Range<usize>, state: &[f64],
        bar: &mut Cotangent, value: &mut f64, observations: &mut usize, cancelled: &mut Cancel,
    ) -> Result<(), TrajectoryError> {
        if range.is_empty() { return Ok(()); }
        let mut state_bar = vec![f64::NAN; self.initial.len()];
        let mut parameter_bar = vec![f64::NAN; self.parameters];
        for sample in range.rev() {
            poll(cancelled)?;
            state_bar.fill(f64::NAN);
            parameter_bar.fill(f64::NAN);
            let term = objective.evaluate(sample, self.sample_times[sample], state,
                &mut state_bar, &mut parameter_bar).map_err(TrajectoryError::Observation)?;
            poll(cancelled)?;
            if !term.is_finite() || !finite(&state_bar) || !finite(&parameter_bar) {
                return Err(TrajectoryError::Observation("non-finite or unwritten observation value/partials".into()));
            }
            *value += term;
            for (total, partial) in bar.initial.iter_mut().zip(&state_bar) { *total += partial; }
            for (total, partial) in bar.parameters.iter_mut().zip(&parameter_bar) { *total += partial; }
            if !value.is_finite() || !finite(&bar.initial) || !finite(&bar.parameters) {
                return Err(AdjointError::NonFiniteAccumulation.into());
            }
            *observations += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;

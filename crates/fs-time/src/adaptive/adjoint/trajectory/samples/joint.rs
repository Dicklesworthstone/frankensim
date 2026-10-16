//! Nonseparable objectives across fixed-time observations of one trajectory.
//!
//! Retain scalar predictions, not state histories. One checked forward replay
//! gathers them; a joint loss produces their cotangents; the existing sampled
//! reverse sweep propagates those cotangents. Both passes share ONE replay cap.
//! This covers cross-time correlations without pretending the loss is additive.

use super::*;

/// A scalar observation map and a joint loss over all its predictions. Every
/// callback must be deterministic and unchanged throughout this operation.
/// Observation partials hold time fixed. `loss` must overwrite all prediction
/// and direct-parameter partials, including zeros, and poll long-running work.
/// Statistical meaning, covariance validity and derivative correctness belong
/// to the caller; this interface does not confer those guarantees.
pub trait JointSampleObjective {
    fn observe(&self, sample: usize, time: f64, state: &[f64]) -> Result<f64, String>;
    fn loss<Cancel: FnMut() -> bool>(
        &self, predictions: &[f64], prediction_bar: &mut [f64],
        direct_parameter_bar: &mut [f64], cancelled: &mut Cancel,
    ) -> Result<f64, TrajectoryError>;
    #[allow(clippy::too_many_arguments)]
    fn observation_vjp(
        &self, sample: usize, time: f64, state: &[f64], seed: f64,
        state_bar: &mut [f64], parameter_bar: &mut [f64],
    ) -> Result<(), String>;
}

struct Seeded<'a, O> { objective: &'a O, predictions: &'a [f64], seeds: &'a [f64] }
impl<O: JointSampleObjective> SampleObjective for Seeded<'_, O> {
    fn evaluate(&self, sample: usize, time: f64, state: &[f64], x: &mut [f64], p: &mut [f64])
        -> Result<f64, String>
    {
        // The ODE endpoint hash alone cannot detect a changed sensor map.
        let replayed = self.objective.observe(sample, time, state)?;
        if replayed.to_bits() != self.predictions[sample].to_bits() {
            return Err(format!("joint observation replay mismatch at sample {sample}"));
        }
        self.objective.observation_vjp(sample, time, state, self.seeds[sample], x, p)?;
        // The nonseparable objective value is computed once, not once per row.
        Ok(0.0)
    }
}

fn allocated(count: usize) -> Result<Vec<f64>, TrajectoryError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| TrajectoryError::Allocation)?;
    values.resize(count, f64::NAN);
    Ok(values)
}

impl<M: OdeVjp> RecordedRk45<'_, M> {
    /// Differentiate a joint objective on a COMPLETE observation-aligned
    /// recording. The accepted mesh and observation times are held fixed.
    ///
    /// In addition to the existing step workspace and reverse checkpoints,
    /// this owns 2*m+p scalars (predictions, their seeds, direct partials), plus
    /// one live replay state while collecting predictions. `max_observations`
    /// caps m before allocation/callbacks; model/loss-owned memory is separate.
    /// No observed state histories or Jacobian matrices are retained.
    ///
    /// `gradient.replayed_steps` includes BOTH the collection pass and reverse
    /// recomputation. Cancellation, incomplete callbacks, changed predictions
    /// or any budget failure return no objective/gradient. Retry redoes the
    /// whole sweep without altering the retained forward recording.
    pub fn pullback_joint<O: JointSampleObjective, Cancel: FnMut() -> bool>(
        &self, objective: &O, budget: ReplayBudget, max_observations: usize,
        cancelled: &mut Cancel,
    ) -> Result<SampledGradient, TrajectoryError> {
        self.model_shape()?;
        if self.state.t != self.config.end { return Err(TrajectoryError::Incomplete); }
        let m = self.sample_times.len();
        if m == 0 || m > max_observations
            || m.checked_mul(2).and_then(|n| n.checked_add(self.parameters)).is_none()
        {
            return Err(AdjointError::InvalidInput("joint observations exceed the nonempty sample cap").into());
        }
        workspace_size(self.initial.len(), self.parameters, self.config.max_workspace_components)?;
        let required = fs_ad::revolve::min_budget(self.records.len());
        if budget.checkpoints < required {
            return Err(TrajectoryError::CheckpointLimit { required, limit: budget.checkpoints });
        }
        // Every leaf needs a second replay even before midpoint recomputation.
        if self.records.len() > budget.replayed_steps / 2 { return Err(TrajectoryError::ReplayLimit); }
        poll(cancelled)?;
        let mut predictions = allocated(m)?;
        let mut progress = Progress { replays: 0, peak: 0, budget };
        {
            let start = self.records.first().map_or(self.state.t, |r| r.start);
            let mut next_sample = 0;
            self.collect_joint(objective, start, &self.initial, &mut predictions, &mut next_sample, cancelled)?;
            let mut current = self.initial.clone();
            for step in 0..self.records.len() {
                current = self.checked_replay(step, &current, &mut progress, cancelled)?.next;
                self.collect_joint(objective, self.records[step].end, &current,
                    &mut predictions, &mut next_sample, cancelled)?;
            }
            if next_sample != m {
                return Err(AdjointError::InvalidInput("joint observation has no recorded endpoint").into());
            }
        }
        let mut seeds = allocated(m)?;
        let mut direct = allocated(self.parameters)?;
        poll(cancelled)?;
        let value = objective.loss(&predictions, &mut seeds, &mut direct, cancelled)?;
        poll(cancelled)?;
        if !value.is_finite() || !finite(&seeds) || !finite(&direct) {
            return Err(TrajectoryError::Observation("non-finite or unwritten joint loss/partials".into()));
        }
        let remaining = ReplayBudget { replayed_steps: budget.replayed_steps - progress.replays, ..budget };
        let seeded = Seeded { objective, predictions: &predictions, seeds: &seeds };
        let mut result = self.pullback_samples(&seeded, remaining, cancelled)?;
        for (total, partial) in result.gradient.parameters.iter_mut().zip(direct) { *total += partial; }
        if !finite(&result.gradient.parameters) { return Err(AdjointError::NonFiniteAccumulation.into()); }
        result.gradient.replayed_steps += progress.replays; // bounded by the shared cap
        result.value = value;
        poll(cancelled)?;
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_joint<O: JointSampleObjective, Cancel: FnMut() -> bool>(
        &self, objective: &O, time: f64, state: &[f64], predictions: &mut [f64],
        next: &mut usize, cancelled: &mut Cancel,
    ) -> Result<(), TrajectoryError> {
        while *next < predictions.len() && self.sample_times[*next] == time {
            poll(cancelled)?;
            let value = objective.observe(*next, time, state).map_err(TrajectoryError::Observation)?;
            poll(cancelled)?;
            if !value.is_finite() { return Err(TrajectoryError::Observation("non-finite joint observation".into())); }
            predictions[*next] = value;
            *next += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;

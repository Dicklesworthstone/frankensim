//! Additive observation objectives on one recorded IMEX trajectory.
//!
//! Observations name accepted endpoints by step number: zero is the initial
//! state, and N is the endpoint after N steps. Repeated indices permit several
//! sensors at one endpoint. This fixed-step interface does not interpolate
//! observations, round arbitrary times to steps, or change the forward mesh.

use super::{
    Cotangent, ImexReplayBudget, ImexTrajectoryError, ImexTrajectoryGradient, ImexVjp, Progress,
    RecordedImex2, accumulate, finite, imex_poll,
};
use fs_solver::FlexiblePreconditioner;

/// The same objective interface used by recorded RK45. Here `sample` indexes
/// the caller's endpoint-index array; `time` is that endpoint's recorded clock.
pub use crate::adaptive::adjoint::trajectory::samples::SampleObjective;

#[derive(Debug, Clone, PartialEq)]
pub struct ImexSampledGradient {
    /// Sum of the supplied scalar observation objectives.
    pub value: f64,
    pub gradient: ImexTrajectoryGradient,
    pub observations: usize,
}

impl<M: ImexVjp, P: FlexiblePreconditioner> RecordedImex2<'_, M, P> {
    /// Accumulate time-series objectives in one bounded checkpoint sweep.
    ///
    /// `indices` must be nonempty, nondecreasing, inside `0..=accepted_steps`,
    /// and fit `max_samples`. The cap is checked before scanning the indices.
    /// Indices, the step size and all model/solver policies are fixed during
    /// differentiation. Objective callbacks must be pure and overwrite all
    /// state and direct-parameter partials, including zeros. Noise weights,
    /// sensor calibration and dependence between readings belong to the caller.
    ///
    /// Each sample is evaluated once in reverse declaration order, after its
    /// endpoint has passed replay verification. Initial samples act directly
    /// on the initial cotangent. Observed states are not retained, and no
    /// per-parameter trajectory is computed. The callback uses n+p scalar
    /// scratch, dropped before the step pullback, within the existing per-step
    /// workspace ceiling. Checkpoints, records and the live cotangent remain
    /// separately accounted for as in `pullback`.
    ///
    /// The caller chains the returned initial cotangent through any parameter
    /// dependence of the initial state. Direct objective partials supplied by
    /// each callback are included once. Failure/cancellation returns neither
    /// a partial objective nor a partial gradient and leaves the recording
    /// retryable. Long objective callbacks must bound their own work.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn pullback_samples<O, Q, Cancel>(
        &self,
        indices: &[usize],
        max_samples: usize,
        objective: &O,
        adjoint_preconditioner: &Q,
        budget: ImexReplayBudget,
        cancelled: &mut Cancel,
    ) -> Result<ImexSampledGradient, ImexTrajectoryError>
    where
        O: SampleObjective + ?Sized,
        Q: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        self.model_shape()?;
        if self.records.len() != self.config.steps {
            return Err(ImexTrajectoryError::Incomplete);
        }
        if indices.is_empty() || indices.len() > max_samples {
            return Err(ImexTrajectoryError::InvalidInput(
                "nonempty observation indices must fit the sample cap",
            ));
        }
        imex_poll(cancelled)?;
        let mut previous = 0;
        for chunk in indices.chunks(256) {
            imex_poll(cancelled)?;
            for &index in chunk {
                if index < previous || index > self.records.len() {
                    return Err(ImexTrajectoryError::InvalidInput(
                        "observation indices must be ordered accepted endpoints",
                    ));
                }
                previous = index;
            }
        }
        let required = self.required_checkpoints();
        if budget.checkpoints < required {
            return Err(ImexTrajectoryError::CheckpointLimit {
                required,
                limit: budget.checkpoints,
            });
        }
        let mut progress = Progress {
            replays: 0,
            peak: 0,
            budget,
        };
        let mut value = 0.0;
        let mut observations = 0;
        let bar = Cotangent {
            initial: vec![0.0; self.initial.len()],
            parameters: vec![0.0; self.parameters],
        };
        let mut observe =
            |endpoint: usize, state: &[f64], bar: &mut Cotangent, cancelled: &mut Cancel| {
                let first = indices.partition_point(|index| *index < endpoint);
                let last = indices.partition_point(|index| *index <= endpoint);
                if first == last {
                    return Ok(());
                }
                let time = if endpoint == 0 {
                    self.initial_time
                } else {
                    self.records[endpoint - 1].time
                };
                let mut state_bar = vec![f64::NAN; self.initial.len()];
                let mut parameter_bar = vec![f64::NAN; self.parameters];
                for sample in (first..last).rev() {
                    imex_poll(cancelled)?;
                    state_bar.fill(f64::NAN);
                    parameter_bar.fill(f64::NAN);
                    let term =
                        objective.evaluate(sample, time, state, &mut state_bar, &mut parameter_bar);
                    imex_poll(cancelled)?;
                    let term = term.map_err(ImexTrajectoryError::Observation)?;
                    if !term.is_finite() || !finite(&state_bar) || !finite(&parameter_bar) {
                        return Err(ImexTrajectoryError::Observation(
                            "nonfinite or unwritten observation value/partials".into(),
                        ));
                    }
                    value += term;
                    if !value.is_finite() {
                        return Err(super::ImexAdjointError::NonFiniteAccumulation.into());
                    }
                    accumulate(&mut bar.initial, &state_bar)?;
                    accumulate(&mut bar.parameters, &parameter_bar)?;
                    observations += 1;
                }
                Ok(())
            };
        let mut bar = if self.records.is_empty() {
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
                &mut observe,
                cancelled,
            )?
        };
        observe(0, &self.initial, &mut bar, cancelled)?;
        debug_assert_eq!(observations, indices.len());
        imex_poll(cancelled)?;
        Ok(ImexSampledGradient {
            value,
            observations,
            gradient: ImexTrajectoryGradient {
                initial: bar.initial,
                parameters: bar.parameters,
                replayed_steps: progress.replays,
                peak_checkpoints: progress.peak,
            },
        })
    }
}

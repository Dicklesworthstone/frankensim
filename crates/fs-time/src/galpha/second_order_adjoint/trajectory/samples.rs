//! Additive displacement, velocity and acceleration observations on a recorded
//! structural trajectory. Endpoint zero is the supplied initial state; repeated
//! indices allow several sensors at one accepted endpoint.

use super::{
    Cotangent, Progress, RecordedStructural, SecondOrderAdjointError, SecondOrderState,
    StructuralReplayBudget, StructuralTrajectoryError, StructuralTrajectoryGradient,
    StructuralTrajectoryModel, accumulate, finite, structural_poll,
};
use fs_solver::FlexiblePreconditioner;

/// A scalar observation and its partial derivatives at one accepted endpoint.
/// The callback must be pure for a fixed model/parameter point and overwrite
/// every output component, including zeros. It must bound its own work.
pub trait StructuralSampleObjective {
    /// `sample` indexes the caller's observation-index array. `state` contains
    /// the actual accepted clock, absolute step counter and q/v/a, with empty
    /// solver history. Write the partials in displacement, velocity and
    /// acceleration order, holding the other fields fixed; `parameter_bar`
    /// contains only direct parameter partials at this fixed endpoint state.
    fn evaluate(
        &self,
        sample: usize,
        state: &SecondOrderState,
        state_bar: (&mut [f64], &mut [f64], &mut [f64]),
        parameter_bar: &mut [f64],
    ) -> Result<f64, String>;
}

/// A complete sampled objective and its discrete trajectory gradient.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralSampledGradient {
    /// Sum of all scalar observation terms, in reverse declaration order.
    pub value: f64,
    /// Initial q/v/a cotangents and total model/load/direct parameter partials.
    pub gradient: StructuralTrajectoryGradient,
    /// Number of observation callbacks successfully accumulated.
    pub observations: usize,
}

impl<M: StructuralTrajectoryModel + ?Sized> RecordedStructural<'_, M> {
    /// Accumulate structural measurements in one bounded checkpoint sweep.
    ///
    /// `indices` must be nonempty, nondecreasing, within `0..=accepted_steps`,
    /// and fit `max_samples`. Indices count steps relative to this recording,
    /// independently of the supplied state's absolute step counter. Repeats
    /// are allowed. The cap is checked before scanning the indices. Sampling
    /// does not interpolate states, round observation times or change steps.
    ///
    /// Each sample is evaluated once in reverse declaration order, after its
    /// endpoint passes replay verification. Initial observations contribute
    /// directly to the initial cotangent. Observed states are not retained;
    /// all parameters share one reverse sweep. Scratch is `3*n+p` scalars,
    /// dropped before the step pullback and within the admitted per-step
    /// workspace. Checkpoints, records and the live cotangent are separately
    /// accounted for as in `pullback`.
    ///
    /// Time, step size, observation indices and solver policies are fixed for
    /// differentiation. The caller chains initial q/v/a cotangents through any
    /// parameter-dependent initial state, including consistent acceleration.
    /// Each direct objective partial is included once. A refusal or cancellation
    /// returns no partial result and leaves the recording unchanged for retry.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn pullback_samples<O, P, Cancel>(
        &self,
        indices: &[usize],
        max_samples: usize,
        objective: &O,
        adjoint_preconditioner: &P,
        budget: StructuralReplayBudget,
        cancelled: &mut Cancel,
    ) -> Result<StructuralSampledGradient, StructuralTrajectoryError>
    where
        O: StructuralSampleObjective + ?Sized,
        P: FlexiblePreconditioner,
        Cancel: FnMut() -> bool,
    {
        self.model_shape()?;
        if self.records.len() != self.config.steps {
            return Err(StructuralTrajectoryError::Incomplete);
        }
        if indices.is_empty() || indices.len() > max_samples {
            return Err(StructuralTrajectoryError::InvalidInput(
                "nonempty observation indices must fit the sample cap",
            ));
        }
        structural_poll(cancelled)?;
        let mut previous = 0;
        for chunk in indices.chunks(256) {
            structural_poll(cancelled)?;
            for &index in chunk {
                if index < previous || index > self.records.len() {
                    return Err(StructuralTrajectoryError::InvalidInput(
                        "observation indices must be ordered accepted endpoints",
                    ));
                }
                previous = index;
            }
        }
        let required = self.required_checkpoints();
        if budget.checkpoints < required {
            return Err(StructuralTrajectoryError::CheckpointLimit {
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
        let n = self.initial.q.len();
        let bar = Cotangent {
            q: vec![0.0; n],
            v: vec![0.0; n],
            a: vec![0.0; n],
            parameters: vec![0.0; self.parameters],
        };
        let mut observe = |endpoint: usize,
                           state: &SecondOrderState,
                           bar: &mut Cotangent,
                           cancelled: &mut Cancel| {
            let first = indices.partition_point(|index| *index < endpoint);
            let last = indices.partition_point(|index| *index <= endpoint);
            if first == last {
                return Ok(());
            }
            let mut q_bar = vec![f64::NAN; n];
            let mut v_bar = vec![f64::NAN; n];
            let mut a_bar = vec![f64::NAN; n];
            let mut parameter_bar = vec![f64::NAN; self.parameters];
            for sample in (first..last).rev() {
                structural_poll(cancelled)?;
                q_bar.fill(f64::NAN);
                v_bar.fill(f64::NAN);
                a_bar.fill(f64::NAN);
                parameter_bar.fill(f64::NAN);
                let term = objective.evaluate(
                    sample,
                    state,
                    (&mut q_bar, &mut v_bar, &mut a_bar),
                    &mut parameter_bar,
                );
                structural_poll(cancelled)?;
                let term = term.map_err(StructuralTrajectoryError::Observation)?;
                if !term.is_finite()
                    || !finite(&q_bar)
                    || !finite(&v_bar)
                    || !finite(&a_bar)
                    || !finite(&parameter_bar)
                {
                    return Err(StructuralTrajectoryError::Observation(
                        "nonfinite or unwritten observation value/partials".into(),
                    ));
                }
                value += term;
                if !value.is_finite() {
                    return Err(SecondOrderAdjointError::NonFiniteAccumulation.into());
                }
                accumulate(&mut bar.q, &q_bar)?;
                accumulate(&mut bar.v, &v_bar)?;
                accumulate(&mut bar.a, &a_bar)?;
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
        structural_poll(cancelled)?;
        Ok(StructuralSampledGradient {
            value,
            observations,
            gradient: StructuralTrajectoryGradient {
                initial_q: bar.q,
                initial_v: bar.v,
                initial_a: bar.a,
                parameters: bar.parameters,
                replayed_steps: progress.replays,
                peak_checkpoints: progress.peak,
            },
        })
    }
}

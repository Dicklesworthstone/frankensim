//! Cycle-mean observations and their ordinary real Fourier-coordinate covectors.

use super::{
    ForcedSensitivityError, PreparedForcedSensitivity, Result, invalid, poll, positive,
    project_direction, scaled,
};
use crate::{AftSum, AftWorkspace, Packing, normalized_difference};
use fs_math::det;

/// Independent numerical acceptance for a differentiable cycle-mean objective.
#[derive(Clone, Copy, Debug)]
pub struct ForcedMeanObjectiveControls {
    /// Positive physical objective unit. It conditions only defect acceptance;
    /// the returned value and gradient retain original objective/state units.
    pub objective_scale: f64,
    /// Maximum absolute coarse/fine mean-value difference divided by objective_scale.
    pub maximum_value_sampling_defect: f64,
    /// Maximum L2 coarse/fine covector difference, with each original HB column
    /// scale divided by objective_scale. All real coordinates count.
    pub maximum_gradient_sampling_defect: f64,
}

/// Independently sampled cycle mean and its complete coefficient gradient.
#[derive(Clone, Debug, PartialEq)]
pub struct ForcedMeanObjective {
    /// Fine-grid mean of the original scalar observation, in physical units.
    pub value: f64,
    /// Ordinary real-coordinate gradient: all DC states, then all real states
    /// and all imaginary states for each positive harmonic. Pass this directly
    /// to [`PreparedForcedSensitivity::prepare_adjoint`]; no extra harmonic
    /// multiplicity, sign or conjugation applies.
    pub gradient: Vec<f64>,
    /// Absolute coarse/fine value difference divided by objective_scale.
    pub value_sampling_defect: f64,
    /// L2 difference of the full state/objective-scaled coefficient covectors.
    pub gradient_sampling_defect: f64,
}

impl PreparedForcedSensitivity {
    /// Evaluate `q = mean(f(t, x(t)))` and its coefficient gradient on the two
    /// original cached AFT grids, at fixed forcing frequency and phase origin.
    ///
    /// `observation(t, x, value, state_gradient)` must write the scalar value AND
    /// every entry of its derivative with respect to the full physical state,
    /// including explicit zeros. Both outputs are reset to NaN at every sample;
    /// false, incomplete or nonfinite output refuses the complete query.
    /// The callback is invoked once per sample, with the cached original state.
    /// No original source/tangent calls, orbit solves or factorizations occur.
    ///
    /// The shared centered AFT projection supplies the covector. For a real
    /// mean-square state, `q = X_0^2 + 2*sum(|X_n|^2)`, the gradient is `2*X_0`,
    /// `4*Re X_n`, `4*Im X_n`. Agreement of objective values alone is insufficient:
    /// the complete scaled coefficient gradients must independently agree.
    ///
    /// Explicit objective parameter, frequency and time-chart derivatives are
    /// separate partials for the adjoint's `explicit_derivative`; this operation
    /// does not infer them. Autonomous phase/amplitude closures are refused.
    /// The caller supplies a differentiable observation and must independently
    /// refine harmonic order; two-grid agreement is not an infinite-tail proof.
    ///
    /// Two O(N) vectors and the existing AFT scratch fit the original sensitivity
    /// payload/work allowance for one query. Callback-internal work is excluded.
    /// Cancellation and numerical refusal publish no partial objective.
    pub fn mean_objective(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        mut observation: impl FnMut(f64, &[f64], &mut f64, &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective> {
        poll(&mut cancelled)?;
        if self.pack.omega_slot.is_some()
            || self.phase.is_some()
            || self.amplitude.is_some()
            || !positive(controls.objective_scale)
            || !positive(controls.maximum_value_sampling_defect)
            || !positive(controls.maximum_gradient_sampling_defect)
        {
            return Err(invalid(
                "forced mean objective requires positive finite objective units and sampling tolerances",
            ));
        }
        let (coarse_value, coarse_gradient) = sample_objective(
            &self.pack,
            self.omega,
            &mut self.coarse,
            &mut observation,
            &mut cancelled,
        )?;
        let (value, gradient) = sample_objective(
            &self.pack,
            self.omega,
            &mut self.fine,
            &mut observation,
            &mut cancelled,
        )?;
        let value_sampling_defect =
            scaled(coarse_value - value, 1.0, controls.objective_scale)?.abs();
        let mut gradient_sampling_defect = 0.0;
        for ((coarse, fine), unit) in coarse_gradient
            .iter()
            .zip(&gradient)
            .zip(&self.scales.columns)
        {
            poll(&mut cancelled)?;
            gradient_sampling_defect = det::hypot(
                gradient_sampling_defect,
                scaled(coarse - fine, *unit, controls.objective_scale)?,
            );
        }
        if !gradient_sampling_defect.is_finite()
            || value_sampling_defect > controls.maximum_value_sampling_defect
            || gradient_sampling_defect > controls.maximum_gradient_sampling_defect
        {
            return Err(ForcedSensitivityError::ObjectiveSamplingUnresolved {
                value_defect: value_sampling_defect,
                gradient_defect: gradient_sampling_defect,
                value_limit: controls.maximum_value_sampling_defect,
                gradient_limit: controls.maximum_gradient_sampling_defect,
            });
        }
        poll(&mut cancelled)?;
        Ok(ForcedMeanObjective {
            value,
            gradient,
            value_sampling_defect,
            gradient_sampling_defect,
        })
    }
}

fn sample_objective(
    pack: &Packing,
    omega: f64,
    aft: &mut AftWorkspace,
    observation: &mut impl FnMut(f64, &[f64], &mut f64, &mut [f64]) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Vec<f64>)> {
    let samples = aft.samples;
    let mut reference = None;
    let mut mean = AftSum::new(0.0);
    // project_direction owns per-sample NaN initialization, state/phase access,
    // cancellation, centered projection and transform scratch invalidation.
    let mut gradient = project_direction(
        pack,
        omega,
        aft,
        &mut |time, state, derivative| {
            let mut value = f64::NAN;
            if !observation(time, state, &mut value, derivative) || !value.is_finite() {
                return false;
            }
            if let Some(first) = reference {
                mean.add(normalized_difference(value, first, samples));
            } else {
                reference = Some(value);
                mean = AftSum::new(value);
            }
            true
        },
        cancelled,
    )
    .map_err(|error| match error {
        ForcedSensitivityError::IslandDirectionRefused { sample, samples } => {
            ForcedSensitivityError::ObjectiveRefused { sample, samples }
        }
        other => other,
    })?;
    for (index, value) in gradient.iter_mut().enumerate() {
        // The residual projector returns -Fourier(f_x). Ordinary derivatives
        // of x = DC + 2*Re(X_n exp(i*n*theta)) undo that sign and double AC.
        *value *= if index < pack.d { -1.0 } else { -2.0 };
    }
    let value = mean.total();
    if !value.is_finite() || gradient.iter().any(|entry| !entry.is_finite()) {
        return Err(invalid(
            "finite mean objective and coefficient gradient required",
        ));
    }
    poll(cancelled)?;
    Ok((value, gradient))
}

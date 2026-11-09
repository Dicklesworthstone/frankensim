//! Cycle means with the autonomous period differentiated at fixed phase.

use super::super::objective::sample_objective;
use super::super::{
    ForcedMeanObjectiveControls, ForcedSensitivityError, Result, invalid, poll, positive, scaled,
};
use super::PreparedAutonomousSensitivity;
use crate::{AftSum, AftWorkspace, Packing, normalized_difference};
use fs_math::det;

/// A cycle mean and its complete coefficient/frequency covector.
#[derive(Clone, Debug, PartialEq)]
pub struct AutonomousMeanObjective {
    /// Fine-grid physical mean of the original observation.
    pub value: f64,
    /// Ordinary real coefficient gradient, excluding the frequency coordinate.
    /// Pass it unchanged to `PreparedAutonomousSensitivity::prepare_adjoint`.
    pub gradient: Vec<f64>,
    /// Partial with respect to the cycle's angular frequency, holding all
    /// Fourier coefficients and their phase origin fixed. Pass as the separate
    /// frequency argument of `prepare_adjoint`.
    pub angular_frequency_partial: f64,
    /// Absolute coarse/fine mean difference divided by objective_scale.
    pub value_sampling_defect: f64,
    /// L2 sampling change of the complete normalized coefficient AND frequency
    /// covector, using the original autonomous coordinate scales.
    pub gradient_sampling_defect: f64,
}

impl PreparedAutonomousSensitivity {
    /// Evaluate a differentiable cycle mean and its full autonomous covector.
    ///
    /// `observation(t, x, value, state_gradient, frequency_partial)` writes the
    /// scalar value, every physical state partial, and its COMPLETE frequency
    /// partial at fixed phase and state. Explicit zeros are required. All
    /// outputs are reset to NaN at each sample; incomplete/nonfinite declarations
    /// refuse. The callback is evaluated once on each cached trajectory sample.
    ///
    /// The mean is an integral over phase: `q = mean_theta f(theta/omega,x,omega)`.
    /// Thus the last output includes `f_omega - (t/omega)*f_t` when the original
    /// observation uses physical time. A state-only mean writes zero. For energy
    /// per cycle from power P(x), write `value = (TAU/omega)*P` and
    /// `frequency_partial = -value/omega`; its state gradient also includes the
    /// same period factor. No change of integration period is inferred silently.
    ///
    /// The shared centered AFT owner supplies the ordinary coefficient gradient.
    /// Independent value and full covector agreement includes the frequency
    /// coordinate; matching values or state partials alone are insufficient.
    /// The phase reference remains fixed, as in the autonomous adjoint. Direct
    /// physical-parameter dependence remains the adjoint's explicit_derivative.
    ///
    /// This reuses both original trajectories, scratch and physical scales,
    /// adding two coefficient vectors and bounded scalar accumulation within the
    /// existing query allowance. It calls no baseline source, tangent, frequency
    /// column, orbit solve or factorization. Callback work is excluded. Two-grid
    /// agreement does not bound omitted harmonics or establish stability.
    pub fn mean_objective(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        mut observation: impl FnMut(f64, &[f64], &mut f64, &mut [f64], &mut f64) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        poll(&mut cancelled)?;
        if !positive(controls.objective_scale)
            || !positive(controls.maximum_value_sampling_defect)
            || !positive(controls.maximum_gradient_sampling_defect)
        {
            return Err(invalid(
                "autonomous mean objective requires positive finite objective units and sampling tolerances",
            ));
        }
        let prepared = &mut self.prepared;
        let (coarse_value, coarse_gradient, coarse_frequency) = sample_with_clock(
            &prepared.pack,
            prepared.omega,
            &mut prepared.coarse,
            &mut observation,
            &mut cancelled,
        )?;
        let (value, gradient, angular_frequency_partial) = sample_with_clock(
            &prepared.pack,
            prepared.omega,
            &mut prepared.fine,
            &mut observation,
            &mut cancelled,
        )?;
        let value_sampling_defect =
            scaled(coarse_value - value, 1.0, controls.objective_scale)?.abs();
        let mut gradient_sampling_defect = scaled(
            coarse_frequency - angular_frequency_partial,
            prepared.scales.columns[prepared.pack.coeff_len()],
            controls.objective_scale,
        )?
        .abs();
        for ((coarse, fine), unit) in coarse_gradient
            .iter()
            .zip(&gradient)
            .zip(&prepared.scales.columns)
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
        Ok(AutonomousMeanObjective {
            value,
            gradient,
            angular_frequency_partial,
            value_sampling_defect,
            gradient_sampling_defect,
        })
    }
}

// The original scalar-objective owner retains value averaging, complete state
// derivative admission, centered Fourier projection and failure invalidation.
// Only the independently authored frequency partial needs another scalar mean.
fn sample_with_clock(
    pack: &Packing,
    omega: f64,
    aft: &mut AftWorkspace,
    observation: &mut impl FnMut(f64, &[f64], &mut f64, &mut [f64], &mut f64) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Vec<f64>, f64)> {
    let samples = aft.samples;
    let mut reference = None;
    let mut frequency_mean = AftSum::new(0.0);
    let (value, gradient) = sample_objective(
        pack,
        omega,
        aft,
        &mut |time, state, value, derivative| {
            let mut frequency_partial = f64::NAN;
            if !observation(time, state, value, derivative, &mut frequency_partial)
                || !frequency_partial.is_finite()
            {
                return false;
            }
            if let Some(first) = reference {
                frequency_mean.add(normalized_difference(frequency_partial, first, samples));
            } else {
                reference = Some(frequency_partial);
                frequency_mean = AftSum::new(frequency_partial);
            }
            true
        },
        cancelled,
    )?;
    let frequency = frequency_mean.total();
    if !frequency.is_finite() {
        return Err(invalid(
            "autonomous objective frequency partial must be finite",
        ));
    }
    Ok((value, gradient, frequency))
}

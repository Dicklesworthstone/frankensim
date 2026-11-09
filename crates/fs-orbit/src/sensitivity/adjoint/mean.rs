//! Direct objective partials sampled on the original cached physical cycles.

use super::{
    Dd, ForcedObjectiveGradient, ForcedParameterPoint, ForcedSensitivityDirection,
    ForcedSensitivityError, PreparedForcedAdjoint, Result, admit_objective_direction, invalid,
    poll, scaled,
};
use crate::{AftSum, AftWorkspace, TAU, normalized_difference};

impl PreparedForcedAdjoint<'_> {
    /// Differentiate a cycle mean with a sampled, explicit objective partial.
    ///
    /// `objective_partial(t, x, out)` writes the COMPLETE direct parameter
    /// partial of the scalar integrand at fixed Fourier coefficients and phase,
    /// including any authored clock/time-chart dependence. It must not include
    /// the state chain; the original objective adjoint supplies that once.
    /// Each output is reset to NaN, and false, incomplete or nonfinite output
    /// refuses. The island and port callbacks have the original split-direction
    /// contract, including all direct material, source and frequency changes.
    ///
    /// Direct partials use both unchanged cached trajectories and compensated
    /// scalar means. Their independent normalized sampling defect AND the
    /// complete `mean(q_p) - lambda^T R_p` defect must satisfy the prepared
    /// derivative allowance. Each contraction uses its own grid's direct mean;
    /// agreement of the residual column alone cannot hide objective aliasing.
    ///
    /// This adds two scalar sample passes, constant scratch and no baseline
    /// source/tangent call, state solve or factorization. Existing query bounds
    /// cover the numerical work; callback internals retain their exclusions.
    /// Two-grid agreement is not a harmonic-tail or stability certificate.
    pub fn gradient_mean_split_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        mut objective_partial: impl FnMut(f64, &[f64], &mut f64) -> bool,
        mut island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        mut port_direction: impl FnMut(ForcedParameterPoint<'_>, &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        admit_objective_direction(direction, 0.0, &mut cancelled)?;
        let explicit = self.sample_explicit(direction, &mut objective_partial, &mut cancelled)?;
        let [coarse, fine] = self.prepared.project_split_direction(
            direction,
            &mut island_direction,
            &mut port_direction,
            &mut cancelled,
        )?;
        self.contract_sampled(direction, explicit, &coarse, &fine, &mut cancelled)
    }

    /// Sample a complete explicit cycle-mean objective partial with a fixed
    /// linear port and prescribed clock. The state chain is supplied by the
    /// original adjoint. For a changing forced clock use the split entry above.
    /// Autonomous callers retain the solved clock through their bordered dual.
    /// Refusal and cancellation leave the original trajectories reusable.
    pub fn gradient_mean_island_direction(
        &mut self,
        direction: ForcedSensitivityDirection,
        objective_partial: impl FnMut(f64, &[f64], &mut f64) -> bool,
        island_direction: impl FnMut(f64, &[f64], &mut [f64]) -> bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        admit_objective_direction(direction, 0.0, &mut cancelled)?;
        if direction.angular_frequency_derivative != 0.0 {
            return Err(invalid(
                "island-only mean objective direction requires a fixed physical clock",
            ));
        }
        self.gradient_mean_split_direction(
            direction,
            objective_partial,
            island_direction,
            |_, out| {
                out.fill(0.0);
                true
            },
            cancelled,
        )
    }

    fn sample_explicit(
        &self,
        direction: ForcedSensitivityDirection,
        objective_partial: &mut impl FnMut(f64, &[f64], &mut f64) -> bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2]> {
        let p = &self.prepared;
        let coarse = sample_partial(p.pack.d, p.omega, &p.coarse, objective_partial, cancelled)?;
        let fine = sample_partial(p.pack.d, p.omega, &p.fine, objective_partial, cancelled)?;
        let defect = scaled(
            (Dd::from_f64(coarse) - Dd::from_f64(fine)).to_f64(),
            direction.parameter_scale,
            self.objective_scale,
        )?
        .abs();
        if defect > p.controls.maximum_derivative_residual {
            return Err(ForcedSensitivityError::DerivativeUnresolved {
                relative_residual: self.relative_residual,
                derivative_residual: defect,
            });
        }
        Ok([coarse, fine])
    }
}

fn sample_partial(
    dimension: usize,
    omega: f64,
    aft: &AftWorkspace,
    callback: &mut impl FnMut(f64, &[f64], &mut f64) -> bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut reference = None;
    let mut mean = AftSum::new(0.0);
    for sample in 0..aft.samples {
        poll(cancelled)?;
        let time = (TAU / omega) * sample as f64 / aft.samples as f64;
        let mut value = f64::NAN;
        if !callback(
            time,
            &aft.xs[sample * dimension..(sample + 1) * dimension],
            &mut value,
        ) || !value.is_finite()
        {
            return Err(ForcedSensitivityError::ObjectiveRefused {
                sample,
                samples: aft.samples,
            });
        }
        if let Some(first) = reference {
            mean.add(normalized_difference(value, first, aft.samples));
        } else {
            reference = Some(value);
            mean = AftSum::new(value);
        }
    }
    let value = mean.total();
    if !value.is_finite() {
        return Err(invalid("finite explicit cycle-mean derivative required"));
    }
    poll(cancelled)?;
    Ok(value)
}

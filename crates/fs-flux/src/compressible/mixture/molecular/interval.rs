//! Shared requested-duration control for the original molecular field owners.
//!
//! Geometry, source evaluations, SSP stages, conservative divergence and private
//! field/load buffers remain with the caller. This controller owns only time,
//! attempt bounds and accepted-stage extrema. A refused call never publishes an
//! owner; successful trials may replace only the caller's private scratch.

use fs_gas::multicomponent::dissipation::{
    DissipationControl, DissipationError, TransportFaceReport,
};

#[derive(Clone, Copy, Default)]
pub(in crate::compressible) struct Stage {
    pub rate: f64,
    pub gradient_residual: f64,
    pub mass_frame_residual: f64,
}
impl Stage {
    pub fn add(&mut self, report: TransportFaceReport) {
        self.rate = self.rate.max(report.diffusivity_m2_s);
        self.gradient_residual = self
            .gradient_residual
            .max(report.maximum_gradient_residual_per_m);
        self.mass_frame_residual = self
            .mass_frame_residual
            .max(report.mass_frame_residual_kg_m2_s.abs());
    }
    pub fn merge(self, other: Self) -> Self {
        Self {
            rate: self.rate.max(other.rate),
            gradient_residual: self.gradient_residual.max(other.gradient_residual),
            mass_frame_residual: self.mass_frame_residual.max(other.mass_frame_residual),
        }
    }
}

#[derive(Clone, Copy)]
pub(in crate::compressible) enum Action {
    /// Compute current-state faces once for both substep selection and stage one.
    First,
    /// Finish both SSP stages. `cached` means First's faces remain available;
    /// later retries must reconstruct them after stage two reused their buffers.
    /// Return combined first/second-stage extrema only after the private field
    /// and interval loads have been admitted. Never publish accepted owners here.
    Trial { dt: f64, first: Stage, cached: bool },
}

#[derive(Clone, Copy, Default)]
pub(in crate::compressible) struct Interval {
    pub elapsed: f64,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub maximum_parabolic_number: f64,
    pub maximum_gradient_residual_per_m: f64,
    pub maximum_mass_frame_residual_kg_m2_s: f64,
}

pub(in crate::compressible) fn advance<C: FnMut() -> bool>(
    duration: f64,
    controls: DissipationControl,
    cancelled: &mut C,
    mut evaluate: impl FnMut(Action, &mut C) -> Result<Stage, DissipationError>,
) -> Result<Interval, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.is_finite() || duration < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite molecular duration required",
        ));
    }
    let mut result = Interval::default();
    while result.elapsed < duration {
        poll(cancelled)?;
        work_limit(result, controls.maximum_steps)?;
        let first = evaluate(Action::First, cancelled)?;
        if !first.rate.is_finite() || first.rate < 0.0 {
            return Err(DissipationError::Invalid(
                "molecular parabolic scale overflow",
            ));
        }
        let remaining = duration - result.elapsed;
        let mut dt = remaining.min(controls.maximum_step_s);
        if first.rate > 0.0 {
            dt = dt.min(controls.parabolic_cfl / first.rate);
        }
        for halving in 0..=controls.maximum_halvings {
            poll(cancelled)?;
            work_limit(result, controls.maximum_steps)?;
            if !dt.is_finite()
                || dt <= 0.0
                || result.elapsed + dt <= result.elapsed
                || (dt < controls.minimum_step_s && dt != remaining)
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable molecular substep",
                ));
            }
            match evaluate(
                Action::Trial {
                    dt,
                    first,
                    cached: halving == 0,
                },
                cancelled,
            ) {
                Ok(stage) => {
                    result.maximum_parabolic_number =
                        result.maximum_parabolic_number.max(dt * stage.rate);
                    result.maximum_gradient_residual_per_m = result
                        .maximum_gradient_residual_per_m
                        .max(stage.gradient_residual);
                    result.maximum_mass_frame_residual_kg_m2_s = result
                        .maximum_mass_frame_residual_kg_m2_s
                        .max(stage.mass_frame_residual);
                    result.elapsed = if dt == remaining {
                        duration
                    } else {
                        result.elapsed + dt
                    };
                    result.accepted_steps += 1;
                    break;
                }
                Err(DissipationError::Cancelled) => return Err(DissipationError::Cancelled),
                Err(error) => {
                    result.rejected_steps += 1;
                    if halving == controls.maximum_halvings || dt * 0.5 < controls.minimum_step_s {
                        return Err(DissipationError::StepFailure {
                            elapsed_time_s: result.elapsed,
                            step_s: dt,
                            last: Box::new(error),
                        });
                    }
                    dt *= 0.5;
                }
            }
        }
    }
    poll(cancelled)?;
    Ok(result)
}

fn work_limit(result: Interval, maximum: usize) -> Result<(), DissipationError> {
    if result.accepted_steps + result.rejected_steps >= maximum {
        Err(DissipationError::WorkLimit {
            elapsed_time_s: result.elapsed,
            accepted_steps: result.accepted_steps,
            rejected_steps: result.rejected_steps,
        })
    } else {
        Ok(())
    }
}

pub(in crate::compressible) fn parabolic(
    observed: f64,
    maximum: f64,
) -> Result<(), DissipationError> {
    if !observed.is_finite() || observed > maximum * (1.0 + 32.0 * f64::EPSILON) {
        Err(DissipationError::ParabolicNumber { observed, maximum })
    } else {
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    if cancelled() {
        Err(DissipationError::Cancelled)
    } else {
        Ok(())
    }
}

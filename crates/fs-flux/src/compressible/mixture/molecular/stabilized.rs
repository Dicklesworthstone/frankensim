//! Shared bounded RKL2 interval control for original multidimensional fields.
//!
//! The original gas owner supplies the Legendre coefficients. Field owners
//! retain source evaluations, geometry, conservative divergence, caloric gates,
//! signed boundary quadrature and private trial buffers. This module selects
//! macrosteps and accounts for all work; it never evaluates substitute physics.

use super::interval::{Interval, Stage};
use super::{DissipationError, poll};
pub use fs_gas::multicomponent::dissipation::{Rkl2Scheme, StabilizedDissipationControl};

mod planar;
pub use planar::StabilizedMolecularReport;

#[derive(Clone, Copy)]
pub(in crate::compressible) enum Action {
    /// The controller has reserved one actual RHS evaluation. Prepare the
    /// current state and retain its original face fluxes as the trial anchor.
    First,
    /// Finish every original RKL2 stage and its endpoint admission. A retry has
    /// overwritten face scratch, so `cached=false` must reserve and reconstruct
    /// the anchor evaluation. All other stage/RHS calls are reserved by the
    /// field owner through `reserve_evaluation`, including failed prefixes.
    Trial {
        dt: f64,
        scheme: Rkl2Scheme,
        first: Stage,
        cached: bool,
    },
}

#[derive(Clone, Copy, Default)]
pub(in crate::compressible) struct StabilizedInterval {
    pub interval: Interval,
    pub operator_evaluations: usize,
    pub maximum_stages_used: usize,
    pub maximum_scaled_parabolic_number: f64,
}

/// Reserve before touching any source/constitutive workspace. A failed call
/// consumes its actual attempted evaluation; composing owners retain the same
/// counter across both molecular halves and discarded complete-interval probes.
pub(in crate::compressible) fn reserve_evaluation(
    evaluations: &mut usize,
    maximum: usize,
    elapsed: f64,
) -> Result<(), DissipationError> {
    if *evaluations >= maximum {
        return Err(DissipationError::OperatorWorkLimit {
            elapsed_time_s: elapsed,
            operator_evaluations: *evaluations,
        });
    }
    *evaluations += 1;
    Ok(())
}

pub(in crate::compressible) fn advance<C: FnMut() -> bool>(
    duration: f64,
    controls: StabilizedDissipationControl,
    evaluations: &mut usize,
    cancelled: &mut C,
    mut evaluate: impl FnMut(Action, f64, &mut usize, &mut C) -> Result<Stage, DissipationError>,
) -> Result<StabilizedInterval, DissipationError> {
    poll(cancelled)?;
    controls.validate()?;
    if !duration.is_finite() || duration < 0.0 {
        return Err(DissipationError::Invalid(
            "nonnegative finite stabilized molecular duration required",
        ));
    }
    let initial_evaluations = *evaluations;
    let base = controls.transport;
    let maximum_scheme = Rkl2Scheme::new(controls.maximum_stages)?;
    let mut result = StabilizedInterval::default();
    while result.interval.elapsed < duration {
        poll(cancelled)?;
        work_limit(result.interval, base.maximum_steps)?;
        reserve_evaluation(
            evaluations,
            controls.maximum_operator_evaluations,
            result.interval.elapsed,
        )?;
        let first = evaluate(
            Action::First,
            result.interval.elapsed,
            evaluations,
            cancelled,
        )?;
        if !first.rate.is_finite() || first.rate < 0.0 {
            return Err(DissipationError::Invalid(
                "stabilized molecular parabolic scale overflow",
            ));
        }
        let remaining = duration - result.interval.elapsed;
        let mut dt = remaining.min(base.maximum_step_s);
        if first.rate > 0.0 {
            dt = dt.min(base.parabolic_cfl * maximum_scheme.stability_extent() / first.rate);
        }
        for halving in 0..=base.maximum_halvings {
            poll(cancelled)?;
            work_limit(result.interval, base.maximum_steps)?;
            if !dt.is_finite()
                || dt <= 0.0
                || result.interval.elapsed + dt <= result.interval.elapsed
                || (dt < base.minimum_step_s && dt != remaining)
            {
                return Err(DissipationError::Invalid(
                    "unrepresentable stabilized molecular substep",
                ));
            }
            let mut scheme = maximum_scheme;
            for stages in 2..=controls.maximum_stages {
                let candidate = Rkl2Scheme::new(stages)?;
                if dt * first.rate <= base.parabolic_cfl * candidate.stability_extent() {
                    scheme = candidate;
                    break;
                }
            }
            match evaluate(
                Action::Trial {
                    dt,
                    scheme,
                    first,
                    cached: halving == 0,
                },
                result.interval.elapsed,
                evaluations,
                cancelled,
            ) {
                Ok(stage) => {
                    let number = dt * stage.rate;
                    super::interval::parabolic(
                        number,
                        base.parabolic_cfl * scheme.stability_extent(),
                    )?;
                    result.interval.maximum_parabolic_number =
                        result.interval.maximum_parabolic_number.max(number);
                    result.maximum_scaled_parabolic_number = result
                        .maximum_scaled_parabolic_number
                        .max(number / scheme.stability_extent());
                    result.maximum_stages_used =
                        result.maximum_stages_used.max(scheme.stage_count());
                    result.interval.maximum_gradient_residual_per_m = result
                        .interval
                        .maximum_gradient_residual_per_m
                        .max(stage.gradient_residual);
                    result.interval.maximum_mass_frame_residual_kg_m2_s = result
                        .interval
                        .maximum_mass_frame_residual_kg_m2_s
                        .max(stage.mass_frame_residual);
                    result.interval.elapsed = if dt == remaining {
                        duration
                    } else {
                        result.interval.elapsed + dt
                    };
                    result.interval.accepted_steps += 1;
                    break;
                }
                Err(
                    error @ (DissipationError::Cancelled
                    | DissipationError::OperatorWorkLimit { .. }
                    | DissipationError::WorkLimit { .. }),
                ) => return Err(error),
                Err(error) => {
                    result.interval.rejected_steps += 1;
                    if halving == base.maximum_halvings || dt * 0.5 < base.minimum_step_s {
                        return Err(DissipationError::StepFailure {
                            elapsed_time_s: result.interval.elapsed,
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
    result.operator_evaluations = *evaluations - initial_evaluations;
    Ok(result)
}

fn work_limit(interval: Interval, maximum: usize) -> Result<(), DissipationError> {
    if interval.accepted_steps + interval.rejected_steps >= maximum {
        Err(DissipationError::WorkLimit {
            elapsed_time_s: interval.elapsed,
            accepted_steps: interval.accepted_steps,
            rejected_steps: interval.rejected_steps,
        })
    } else {
        Ok(())
    }
}

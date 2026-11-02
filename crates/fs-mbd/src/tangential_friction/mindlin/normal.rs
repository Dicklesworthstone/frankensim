//! Monotone normal recoil around the original source-owned Hertz path average.

use super::*;

pub(super) fn solve_normal(
    source: &SourcePointContact,
    start: f64,
    free: f64,
    compliance: f64,
    dt: f64,
    controls: RigidMindlinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, SourceElasticPathAverage, usize), RigidMindlinError> {
    let envelope = source.envelope();
    let maximum_change =
        fs_math::next_down(finite(dt * envelope.maximum_absolute_rate_m_per_s)?).max(0.0);
    let mut lo = finite(start - maximum_change)?.max(0.0);
    let mut hi = finite(start + maximum_change)?.min(envelope.maximum_indentation_m);
    // The source checks the actual quotient, not a rounded product. Move a
    // rounded product and bracket endpoints inward if that quotient exceeds
    // its original rate bound. This changes only the numerical search interval.
    if (start - lo) / dt > envelope.maximum_absolute_rate_m_per_s {
        lo = fs_math::next_up(lo);
    }
    if (hi - start) / dt > envelope.maximum_absolute_rate_m_per_s {
        hi = fs_math::next_down(hi);
    }
    let evaluate = |indentation| -> Result<(SourceElasticPathAverage, f64), RigidMindlinError> {
        let path = source.elastic_path_average(start, indentation, dt)?;
        let residual =
            finite(indentation - start - free + finite(compliance * path.mean_force_n)?)?;
        Ok((path, residual))
    };
    for (indentation, lower) in [(lo, true), (hi, false)] {
        poll(cancelled)?;
        let (path, residual) = evaluate(indentation)?;
        if lower && residual >= 0.0 && indentation == 0.0 {
            return Err(RigidMindlinError::Separation);
        }
        if residual.abs() <= controls.displacement_tolerance_m && indentation > 0.0 {
            return Ok((indentation, path, 0));
        }
        if (lower && residual > 0.0) || (!lower && residual < 0.0) {
            return Err(RigidMindlinError::Unresolved(
                "normal recoil leaves the original rate/indentation envelope",
            ));
        }
    }
    let mut indentation = start.clamp(lo, hi);
    let mut residual = 0.0;
    for iteration in 1..=controls.maximum_normal_iterations {
        poll(cancelled)?;
        let (path, value) = evaluate(indentation)?;
        residual = value;
        if residual.abs() <= controls.displacement_tolerance_m && indentation > 0.0 {
            return Ok((indentation, path, iteration));
        }
        if residual < 0.0 {
            lo = indentation;
        } else {
            hi = indentation;
        }
        let derivative = finite(1.0 + finite(compliance * path.endpoint_derivative_n_per_m)?)?;
        let proposed = indentation - residual / derivative;
        indentation =
            if proposed.is_finite() && proposed > lo && proposed < hi && proposed != indentation {
                proposed
            } else {
                f64::midpoint(lo, hi)
            };
    }
    Err(RigidMindlinError::NormalSolve {
        residual_m: residual,
        iterations: controls.maximum_normal_iterations,
    })
}

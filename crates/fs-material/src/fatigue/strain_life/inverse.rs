//! Private bounded scalar life inversion shared by axial and shear calibrations.
//! This helper contains no material meaning or coefficient substitutions.

use crate::fatigue::FatigueError;
use fs_math::det;

pub(in crate::fatigue) fn invert_life(
    amplitude: f64,
    amplitude_domain: (f64, f64),
    cycle_domain: (f64, f64),
    log_reversal_domain: (f64, f64),
    evaluate: impl Fn(f64) -> (f64, f64),
) -> Result<Option<f64>, FatigueError> {
    let (minimum_amplitude, maximum_amplitude) = amplitude_domain;
    let (minimum_cycles, maximum_cycles) = cycle_domain;
    let (log_minimum_reversals, log_maximum_reversals) = log_reversal_domain;
    if !amplitude.is_finite() || amplitude < 0.0 {
        return Err(FatigueError::Arithmetic);
    }
    if amplitude == 0.0 {
        return Ok(None);
    }
    if !minimum_amplitude.is_finite()
        || minimum_amplitude <= 0.0
        || !maximum_amplitude.is_finite()
        || maximum_amplitude <= minimum_amplitude
    {
        return Err(FatigueError::Arithmetic);
    }
    if amplitude < minimum_amplitude || amplitude > maximum_amplitude {
        return Err(FatigueError::StrainAmplitudeOutsideCurve {
            amplitude,
            minimum: minimum_amplitude,
            maximum: maximum_amplitude,
        });
    }
    if amplitude == maximum_amplitude {
        return Ok(Some(minimum_cycles));
    }
    if amplitude == minimum_amplitude {
        return Ok(Some(maximum_cycles));
    }
    let mut lo = log_minimum_reversals;
    let mut hi = log_maximum_reversals;
    let mut x = lo + 0.5 * (hi - lo);
    for _ in 0..64 {
        let (value, derivative) = evaluate(x);
        let residual = value - amplitude;
        if !value.is_finite() {
            return Err(FatigueError::Arithmetic);
        }
        let next_newton = if derivative.is_finite() && derivative < 0.0 {
            x - residual / derivative
        } else {
            f64::NAN
        };
        // The local log-life step measures relative LIFE error; a small
        // strain residual alone would misjudge very shallow fatigue fits.
        if residual == 0.0
            || (next_newton.is_finite()
                && (next_newton - x).abs() <= 8.0 * f64::EPSILON * x.abs().max(1.0))
        {
            return life_from_log_reversals(x, cycle_domain);
        }
        if residual > 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let midpoint = lo + 0.5 * (hi - lo);
        if midpoint == lo || midpoint == hi {
            return life_from_log_reversals(midpoint, cycle_domain);
        }
        // Force geometric bracket shrinkage when Newton lands too near an
        // endpoint; this retains a bounded convergence guarantee.
        let margin = 0.05 * (hi - lo);
        x = if next_newton.is_finite() && next_newton > lo + margin && next_newton < hi - margin {
            next_newton
        } else {
            midpoint
        };
    }
    Err(FatigueError::Arithmetic)
}

fn life_from_log_reversals(x: f64, cycle_domain: (f64, f64)) -> Result<Option<f64>, FatigueError> {
    let cycles = det::exp(x - core::f64::consts::LN_2);
    if !cycles.is_finite() || cycles < cycle_domain.0 || cycles > cycle_domain.1 {
        return Err(FatigueError::Arithmetic);
    }
    Ok(Some(cycles))
}

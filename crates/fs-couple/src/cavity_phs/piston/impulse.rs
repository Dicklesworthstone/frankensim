//! The existing piston impulse solve, shared by single- and double-acting walls.
//!
//! This is a private coupling implementation detail, not a second public
//! nonlinear-solver API. Material proposals remain owned by the caller.

use super::{PistonError, PistonSolve, invalid, positive};

pub(super) struct Trial<T> {
    pub value: T,
    pub residual: f64,
    pub tangent: f64,
}

pub(super) struct Solution<T> {
    pub value: T,
    pub iterations: usize,
    pub proposals: usize,
}

/// Preserve the original evaluation order and bounded safeguarded Newton step.
/// The callback checks cancellation and never publishes a material proposal.
pub(super) fn solve<T>(
    initial_volume: f64,
    controls: PistonSolve,
    mut evaluate: impl FnMut(f64) -> Result<Trial<T>, PistonError>,
) -> Result<Solution<T>, PistonError> {
    let [mut lo, mut hi] = controls.volume_bracket_m3;
    if !positive(lo) || !positive(hi) || lo >= hi || !initial_volume.is_finite()
        || initial_volume < lo || initial_volume > hi
        || !positive(controls.momentum_tolerance_n_s)
        || !(1..=128).contains(&controls.maximum_iterations)
    {
        return Err(invalid("invalid piston impulse bracket, tolerance or work bound"));
    }
    let first = evaluate(initial_volume)?;
    if first.residual.abs() <= controls.momentum_tolerance_n_s {
        return Ok(Solution { value: first.value, iterations: 0, proposals: 1 });
    }
    let lower = evaluate(lo)?;
    let upper = evaluate(hi)?;
    if lower.residual.abs() <= controls.momentum_tolerance_n_s {
        return Ok(Solution { value: lower.value, iterations: 0, proposals: 3 });
    }
    if upper.residual.abs() <= controls.momentum_tolerance_n_s {
        return Ok(Solution { value: upper.value, iterations: 0, proposals: 3 });
    }
    if lower.residual > 0.0 || upper.residual < 0.0 || lower.residual > upper.residual {
        return Err(PistonError::OutsideBracket {
            lower_impulse_n_s: lower.residual, upper_impulse_n_s: upper.residual,
        });
    }
    // Release rejected source-bearing proposals before the interior loop.
    drop(lower);
    drop(upper);
    let proposed = initial_volume - first.residual / first.tangent;
    let mut volume = if proposed > lo && proposed < hi { proposed }
        else { lo + 0.5 * (hi - lo) };
    let mut residual = first.residual;
    drop(first);
    for iteration in 1..=controls.maximum_iterations {
        let trial = evaluate(volume)?;
        residual = trial.residual;
        if residual.abs() <= controls.momentum_tolerance_n_s {
            return Ok(Solution { value: trial.value, iterations: iteration, proposals: iteration + 3 });
        }
        if residual < 0.0 { lo = volume; } else { hi = volume; }
        let proposed = volume - residual / trial.tangent;
        let margin = 0.1 * (hi - lo);
        let next = if proposed > lo + margin && proposed < hi - margin { proposed }
            else { lo + 0.5 * (hi - lo) };
        if next <= lo || next >= hi || next == volume {
            return Err(PistonError::NotConverged {
                iterations: iteration, impulse_residual_n_s: residual,
            });
        }
        volume = next;
    }
    Err(PistonError::NotConverged {
        iterations: controls.maximum_iterations, impulse_residual_n_s: residual,
    })
}

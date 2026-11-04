//! One bounded whole-step energy fixed point for interface exchange and sources.
//! The physical owners return the actual residual and retain their own state.

use super::{HeatExchangeSolve, ThermalPistonError, positive};

/// Internal numerical failure; physical wrappers attach their own units.
pub(crate) enum ScalarError<E> {
    Evaluation(E),
    Invalid(&'static str),
    Cancelled,
    OutsideBracket { lower: f64, upper: f64 },
    NotConverged { iterations: usize, residual: f64 },
}

pub(crate) fn admit<E>(controls: HeatExchangeSolve) -> Result<(), ThermalPistonError<E>> {
    let [lo, hi] = controls.heat_bracket_j;
    if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo > 0.0 || hi < 0.0
        || !positive(controls.source_tolerance_j) || !positive(controls.energy_tolerance_j)
        || !(1..=128).contains(&controls.maximum_iterations)
    { return Err(ThermalPistonError::Invalid("invalid source/exchange energy bracket or work limits")); }
    Ok(())
}

pub(crate) fn solve<T, E, C: FnMut() -> bool>(
    controls: HeatExchangeSolve,
    evaluate: impl FnMut(f64, usize, &mut C) -> Result<(T, f64), ThermalPistonError<E>>,
    cancelled: &mut C,
) -> Result<T, ThermalPistonError<E>> {
    admit(controls)?;
    solve_scalar(controls.heat_bracket_j, controls.source_tolerance_j,
        controls.maximum_iterations, evaluate, cancelled).map_err(|error| match error {
        ScalarError::Evaluation(error) => error,
        ScalarError::Invalid(_) => ThermalPistonError::Invalid("nonfinite physical energy feedback residual"),
        ScalarError::Cancelled => ThermalPistonError::Cancelled,
        ScalarError::OutsideBracket { lower, upper } => ThermalPistonError::OutsideBracket {
            lower_residual_j: lower, upper_residual_j: upper,
        },
        ScalarError::NotConverged { iterations, residual } => ThermalPistonError::NotConverged {
            iterations, residual_j: residual,
        },
    })
}

/// One safeguarded secant/bisection owner shared by heat and reaction extent.
/// The caller supplies consistent units and admits the physical candidates.
/// A degenerate zero bracket can admit an idle state but cannot make progress.
pub(crate) fn solve_scalar<T, E, C: FnMut() -> bool>(
    bracket: [f64; 2], tolerance: f64, maximum_iterations: usize,
    mut evaluate: impl FnMut(f64, usize, &mut C) -> Result<(T, f64), E>,
    cancelled: &mut C,
) -> Result<T, ScalarError<E>> {
    let [mut lo, mut hi] = bracket;
    if !lo.is_finite() || !hi.is_finite() || lo > hi || lo > 0.0 || hi < 0.0
        || !positive(tolerance) || !(1..=128).contains(&maximum_iterations)
    { return Err(ScalarError::Invalid("invalid scalar bracket or work limits")); }
    let mut sample = |coordinate, trials| {
        if cancelled() { return Err(ScalarError::Cancelled); }
        let (candidate, residual) = evaluate(coordinate, trials, cancelled)
            .map_err(ScalarError::Evaluation)?;
        if !residual.is_finite() {
            return Err(ScalarError::Invalid("nonfinite physical feedback residual"));
        }
        if cancelled() { return Err(ScalarError::Cancelled); }
        Ok((candidate, residual))
    };
    let (first, f0) = sample(0.0, 1)?;
    if f0.abs() <= tolerance { return Ok(first); }
    // Select the side indicated by the actual zero-source residual. This is
    // not a uniqueness/monotonicity claim for an arbitrary constitutive law.
    let edge = if f0 < 0.0 { lo = 0.0; hi } else { hi = 0.0; lo };
    let (end, fe) = sample(edge, 2)?;
    if fe.abs() <= tolerance { return Ok(end); }
    let (mut fl, mut fu) = if f0 < 0.0 { (f0, fe) } else { (fe, f0) };
    if fl > 0.0 || fu < 0.0 || lo == hi {
        return Err(ScalarError::OutsideBracket { lower: fl, upper: fu });
    }
    let mut residual = f0;
    for iteration in 1..=maximum_iterations {
        let fraction = -fl / (fu - fl);
        let q = if fraction.is_finite() && fraction > 0.1 && fraction < 0.9 {
            (1.0 - fraction) * lo + fraction * hi
        } else { 0.5 * lo + 0.5 * hi };
        if q <= lo || q >= hi {
            return Err(ScalarError::NotConverged { iterations: iteration - 1, residual });
        }
        let (candidate, value) = sample(q, iteration + 2)?;
        residual = value;
        if residual.abs() <= tolerance { return Ok(candidate); }
        if residual < 0.0 { lo = q; fl = residual; } else { hi = q; fu = residual; }
    }
    Err(ScalarError::NotConverged { iterations: maximum_iterations, residual })
}


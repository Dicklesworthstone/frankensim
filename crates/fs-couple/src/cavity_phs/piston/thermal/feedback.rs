//! One bounded whole-step energy fixed point for interface exchange and sources.
//! The physical owners return the actual residual and retain their own state.

use super::{HeatExchangeSolve, ThermalPistonError, poll, positive};

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
    mut evaluate: impl FnMut(f64, usize, &mut C) -> Result<(T, f64), ThermalPistonError<E>>,
    cancelled: &mut C,
) -> Result<T, ThermalPistonError<E>> {
    admit(controls)?;
    let mut sample = |energy, trials| {
        poll(cancelled)?;
        let (candidate, residual) = evaluate(energy, trials, cancelled)?;
        if !residual.is_finite() {
            return Err(ThermalPistonError::Invalid("nonfinite physical energy feedback residual"));
        }
        poll(cancelled)?;
        Ok((candidate, residual))
    };
    let [mut lo, mut hi] = controls.heat_bracket_j;
    let (first, f0) = sample(0.0, 1)?;
    if f0.abs() <= controls.source_tolerance_j { return Ok(first); }
    // Select the side indicated by the actual zero-source residual. This is
    // not a uniqueness/monotonicity claim for an arbitrary constitutive law.
    let edge = if f0 < 0.0 { lo = 0.0; hi } else { hi = 0.0; lo };
    let (end, fe) = sample(edge, 2)?;
    if fe.abs() <= controls.source_tolerance_j { return Ok(end); }
    let (mut fl, mut fu) = if f0 < 0.0 { (f0, fe) } else { (fe, f0) };
    if fl > 0.0 || fu < 0.0 || lo == hi {
        return Err(ThermalPistonError::OutsideBracket { lower_residual_j: fl, upper_residual_j: fu });
    }
    let mut residual = f0;
    for iteration in 1..=controls.maximum_iterations {
        let fraction = -fl / (fu - fl);
        let q = if fraction.is_finite() && fraction > 0.1 && fraction < 0.9 {
            (1.0 - fraction) * lo + fraction * hi
        } else { 0.5 * lo + 0.5 * hi };
        if q <= lo || q >= hi {
            return Err(ThermalPistonError::NotConverged { iterations: iteration - 1, residual_j: residual });
        }
        let (candidate, value) = sample(q, iteration + 2)?;
        residual = value;
        if residual.abs() <= controls.source_tolerance_j { return Ok(candidate); }
        if residual < 0.0 { lo = q; fl = residual; } else { hi = q; fu = residual; }
    }
    Err(ThermalPistonError::NotConverged { iterations: controls.maximum_iterations, residual_j: residual })
}

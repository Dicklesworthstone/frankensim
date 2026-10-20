//! Exact schedule integration for declared domain-absorbed heating energy.
//!
//! Constant, ramp and linear/held table schedules are integrated piecewise,
//! including clamped tails. The returned value weights are seconds: E = w.P.
//! This is the integral of the authored input, not the endpoint quadrature a
//! thermal run actually applied. It is independent of source position/width
//! because the existing projection normalizes power absorbed IN the domain.

use super::{MovingHeatSource, invalid, poll};
use crate::ConductionError;
use fs_exec::Cx;
use fs_scenario::signal::{Interp, TimeSignal};

/// Joules in a fixed time window and their derivatives per authored watt.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorbedHeat {
    pub energy_j: f64,
    /// Constant: one entry; ramp: [from,to]; table: original value order.
    /// The derivative with respect to a knot time is not included.
    pub power_value_weights_s: Vec<f64>,
}

impl MovingHeatSource {
    /// Integrate the original power schedule over [start_s,end_s]. Breakpoints
    /// inside the interval are included without a caller-selected time grid.
    /// The existing clamped tails are retained; isolated jump values carry no
    /// duration. A zero-width interval is valid and has a zero gradient.
    ///
    /// # Errors
    /// Nonfinite/reversed/unrepresentable intervals, negative authored powers,
    /// insufficient value allowance, unsupported Chebyshev power, cancellation
    /// or nonfinite arithmetic. No source or physical trajectory is changed.
    pub fn integrate_absorbed_power(&self, cx: &Cx<'_>, start_s: f64,
        end_s: f64, max_values: usize) -> Result<AbsorbedHeat, ConductionError>
    {
        poll(cx, 0)?;
        let duration = end_s - start_s;
        if !start_s.is_finite() || !end_s.is_finite() || !duration.is_finite() || duration < 0.0 {
            return Err(invalid("absorbed-energy window needs finite ordered endpoints and representable duration"));
        }
        // This owner already enforces the representation and control cap before
        // copying the values. Chebyshev positivity needs a separate policy.
        let values = self.absorbed_power_values(max_values)?;
        for (i, &power) in values.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if !power.is_finite() || power < 0.0 {
                return Err(invalid("absorbed-energy integration requires nonnegative finite authored watts"));
            }
        }
        let mut weights = super::zeros(values.len())?;
        match &self.absorbed_power_w {
            TimeSignal::Constant(_) => weights[0] = duration,
            TimeSignal::Ramp { t_start, t_end, .. } => {
                tails(&mut weights, start_s, end_s, *t_start, *t_end)?;
                segment(&mut weights, 0, start_s, end_s, *t_start, *t_end, true)?;
            }
            TimeSignal::Table { times, interp, .. } => {
                tails(&mut weights, start_s, end_s, times[0], times[times.len()-1])?;
                for i in 0..times.len()-1 {
                    if i % 256 == 0 { poll(cx, i)?; }
                    segment(&mut weights, i, start_s, end_s, times[i], times[i+1], *interp == Interp::Linear)?;
                }
            }
            TimeSignal::Chebfun(_) => unreachable!("power-value owner refuses polynomial controls"),
        }
        let mut energy = 0.0;
        for (i, (&power, &weight)) in values.iter().zip(&weights).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            energy = checked(energy + checked(power * weight)?)?;
        }
        poll(cx, values.len())?;
        Ok(AbsorbedHeat { energy_j: energy, power_value_weights_s: weights })
    }
}

fn tails(weights: &mut [f64], start: f64, end: f64, first: f64, last: f64)
    -> Result<(), ConductionError>
{
    if start < first { weights[0] = checked(end.min(first) - start)?; }
    if end > last {
        let n = weights.len()-1;
        weights[n] = checked(weights[n] + (end - start.max(last)))?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn segment(weights: &mut [f64], index: usize, start: f64, end: f64,
    a: f64, b: f64, linear: bool) -> Result<(), ConductionError>
{
    let (lo, hi) = (start.max(a), end.min(b));
    if hi <= lo { return Ok(()); }
    let duration = hi-lo;
    if linear {
        // The mean of an affine function over an interval is its endpoint mean.
        // Average dimensionless fractions rather than large absolute times.
        let alpha = f64::midpoint(fraction(lo, a, b), fraction(hi, a, b));
        weights[index] = checked(weights[index] + duration * (1.0-alpha))?;
        weights[index+1] = checked(weights[index+1] + duration * alpha)?;
    } else {
        weights[index] = checked(weights[index] + duration)?;
    }
    Ok(())
}

fn fraction(t: f64, a: f64, b: f64) -> f64 {
    let width = b-a;
    if width.is_finite() { (t-a)/width }
    else {
        let scale = a.abs().max(b.abs()).max(t.abs());
        (t/scale-a/scale)/(b/scale-a/scale)
    }
}

fn checked(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(invalid("absorbed-energy integral or derivative is not representable")) }
}

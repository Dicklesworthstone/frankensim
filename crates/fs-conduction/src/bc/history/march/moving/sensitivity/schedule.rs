//! Fixed-time chain rules from sampled scan controls to authored history values.
//!
//! Vector entries follow the original signal representation: one constant,
//! [from, to] for a ramp, values in table order, or stored Chebyshev coefficients
//! (including the original half-weighted coefficient zero). Knot times, domains,
//! interpolation choices and the selected side of an event remain fixed.
//! This does not differentiate motion limits, source-support admission or optics.

use super::{MovingHeatGradient, MovingHeatSource, finite, poll, zeros};
use crate::ConductionError;
use fs_exec::Cx;
use fs_scenario::signal::{Interp, TimeSignal};

/// Objective derivatives in the SAME representation as the original source.
/// Entries differentiate values in SI units, not normalized optimizer variables.
#[derive(Debug, Clone, PartialEq)]
pub struct MovingHeatScheduleGradient {
    /// Per-metre derivatives of each center signal's values or coefficients.
    pub center_values_m: [Vec<f64>; 3],
    /// Per-metre derivatives of the three shared physical footprint half widths.
    pub half_width_m: [f64; 3],
    /// Per-watt derivatives of the absorbed-power signal's values or coefficients.
    pub absorbed_power_values_w: Vec<f64>,
}

impl MovingHeatSource {
    /// Total scalar design coordinates: four signals' values plus three widths.
    /// Time knots and Chebyshev domain endpoints are intentionally not controls.
    pub fn schedule_control_count(&self) -> Result<usize, ConductionError> {
        self.signals().try_fold(3usize, |total, signal| total.checked_add(count(signal))
            .ok_or_else(|| super::invalid("moving-source schedule control count overflows")))
    }

    /// Chain a SAME-source, SAME-time sampled control gradient to schedule values.
    /// The source's original right-continuous/clamped evaluation is validated.
    /// `max_controls` bounds all returned scalar entries before allocation.
    /// At zero power only feasible nonnegative-power directions are physical.
    pub fn pullback_schedule_at(&self, cx: &Cx<'_>, time_s: f64,
        gradient: MovingHeatGradient, max_controls: usize)
        -> Result<MovingHeatScheduleGradient, ConductionError>
    { self.schedule_pullback(cx, time_s, false, gradient, max_controls) }

    /// Same chain rule on the exact interval-ending side of a held-table knot.
    /// A heater's turn-off value never receives the preceding interval's gradient.
    pub fn pullback_schedule_before(&self, cx: &Cx<'_>, time_s: f64,
        gradient: MovingHeatGradient, max_controls: usize)
        -> Result<MovingHeatScheduleGradient, ConductionError>
    { self.schedule_pullback(cx, time_s, true, gradient, max_controls) }

    fn schedule_pullback(&self, cx: &Cx<'_>, time: f64, left: bool,
        gradient: MovingHeatGradient, max_controls: usize)
        -> Result<MovingHeatScheduleGradient, ConductionError>
    {
        poll(cx, 0)?;
        if self.schedule_control_count()? > max_controls {
            return Err(super::invalid("moving-source schedule exceeds its scalar-control allowance"));
        }
        // The existing owner checks every quantity, signal and evaluated value.
        self.values(time, left)?;
        for value in gradient.center_m.into_iter().chain(gradient.half_width_m)
            .chain(std::iter::once(gradient.absorbed_power_w)) { finite(value)?; }
        let center_values_m = [
            pull_signal(cx, &self.center_m[0], time, left, gradient.center_m[0])?,
            pull_signal(cx, &self.center_m[1], time, left, gradient.center_m[1])?,
            pull_signal(cx, &self.center_m[2], time, left, gradient.center_m[2])?,
        ];
        let absorbed_power_values_w = pull_signal(cx, &self.absorbed_power_w,
            time, left, gradient.absorbed_power_w)?;
        poll(cx, 0)?;
        Ok(MovingHeatScheduleGradient { center_values_m,
            half_width_m: gradient.half_width_m, absorbed_power_values_w })
    }
}

fn count(signal: &TimeSignal) -> usize {
    match signal {
        TimeSignal::Constant(_) => 1,
        TimeSignal::Ramp { .. } => 2,
        TimeSignal::Table { values, .. } => values.len(),
        TimeSignal::Chebfun(profile) => profile.cheb.coeffs().len(),
    }
}

// Derivative basis only: the primal is always evaluated by the original owner.
fn pull_signal(cx: &Cx<'_>, signal: &TimeSignal, time: f64, left: bool, seed: f64)
    -> Result<Vec<f64>, ConductionError>
{
    poll(cx, 0)?;
    let mut out = zeros(count(signal))?;
    match signal {
        TimeSignal::Constant(_) => out[0] = seed,
        TimeSignal::Ramp { t_start, t_end, .. } => {
            if time <= *t_start { out[0] = seed; }
            else if time >= *t_end { out[1] = seed; }
            else {
                let alpha = fraction(time, *t_start, *t_end);
                out[0] = finite(seed * (1.0 - alpha))?;
                out[1] = finite(seed * alpha)?;
            }
        }
        TimeSignal::Table { times, interp, .. } => {
            let before = times.partition_point(|&t| t < time);
            if left && *interp == Interp::Hold && before > 0
                && times.get(before).is_some_and(|&t| t == time)
            { out[before - 1] = seed; }
            else {
                match times.binary_search_by(|t| t.total_cmp(&time)) {
                    Ok(i) => out[i] = seed,
                    Err(0) => out[0] = seed,
                    Err(i) if i >= times.len() => out[times.len() - 1] = seed,
                    Err(i) => match interp {
                        Interp::Hold => out[i - 1] = seed,
                        Interp::Linear => {
                            let alpha = fraction(time, times[i - 1], times[i]);
                            out[i - 1] = finite(seed * (1.0 - alpha))?;
                            out[i] = finite(seed * alpha)?;
                        }
                    }
                }
            }
        }
        TimeSignal::Chebfun(profile) => {
            let (a, b) = profile.cheb.domain();
            let t = reference(time.clamp(a, b), a, b);
            for (index, (value, basis)) in out.iter_mut().zip(chebyshev_basis(t)).enumerate() {
                if index % 256 == 0 { poll(cx, index)?; }
                *value = finite(seed * basis)?;
            }
        }
    }
    poll(cx, out.len())?;
    Ok(out)
}

fn fraction(t: f64, a: f64, b: f64) -> f64 {
    let width = b - a;
    if width.is_finite() { (t - a) / width }
    else {
        let scale = a.abs().max(b.abs()).max(t.abs());
        ((t / scale) - (a / scale)) / ((b / scale) - (a / scale))
    }
}

// Match the original Chebyshev coordinate map, including its wide-domain path.
fn reference(t: f64, a: f64, b: f64) -> f64 {
    if t == a { return -1.0; }
    if t == b { return 1.0; }
    let width = b - a;
    let doubled = 2.0 * (t - a);
    if width.is_finite() && doubled.is_finite() { doubled / width - 1.0 }
    else { (t - f64::midpoint(a, b)) / (0.5 * b - 0.5 * a) }
}

// Derivatives with respect to fs-cheb's stored coefficients, not DCT samples.
fn chebyshev_basis(t: f64) -> impl Iterator<Item = f64> {
    let (mut previous, mut current) = (1.0, t);
    std::iter::once(0.5).chain(std::iter::from_fn(move || {
        let value = current;
        (previous, current) = (current, 2.0 * t * current - previous);
        Some(value)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_chebyshev_controls_use_half_weighted_constant_and_original_basis() {
        for t in [-1.0_f64, -0.75, 0.0, 0.3, 1.0] {
            let expected = [0.5, t, 2.0*t*t-1.0, 4.0*t*t*t-3.0*t,
                8.0*t.powi(4)-8.0*t*t+1.0];
            for (actual, expected) in chebyshev_basis(t).zip(expected) {
                assert!((actual - expected).abs() < 2e-15);
            }
        }
        assert_eq!(reference(0.0, -f64::MAX, f64::MAX), 0.0);
        assert_eq!(reference(-f64::MAX, -f64::MAX, f64::MAX), -1.0);
        assert_eq!(reference(f64::MAX, -f64::MAX, f64::MAX), 1.0);
        assert_eq!(fraction(0.0, -f64::MAX, f64::MAX), 0.5);
    }
}

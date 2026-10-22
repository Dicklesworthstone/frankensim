//! Resolve source travel independently of endpoint thermal agreement.
//!
//! A returning path can have equal coarse and half-step source samples while
//! heating different locations in between. Bound the L1 speed in ellipsoidal
//! coordinates instead: sum_i |dc_i/dt| / half_width_i. The allowed travel is
//! a fraction of those half widths, not a temperature or energy tolerance.
//! For a Chebyshev path use |T_n'| <= n^2 on its declared domain, including all
//! supplied coefficients. This is conservative in exact arithmetic but is not
//! an interval-arithmetic certificate or a spatial/temporal solution-error bound.

/// Physical path velocities and authored-coordinate derivatives, not timestep caps.
pub mod velocity;

use super::{MovingHeatSource, invalid, poll};
use crate::ConductionError;
use fs_exec::Cx;
use fs_scenario::signal::{Interp, TimeSignal};

/// Explicit geometric scan resolution, independent of the thermal error test.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanMotionLimit {
    /// Maximum sum of per-axis travel divided by the source's corresponding
    /// ellipsoid half width, in (0, 1]. Applied even during zero-power travel.
    pub max_normalized_travel: f64,
    /// Total coefficients admitted across this source's three Chebyshev axes.
    /// Zero permits only constant, ramp and tabulated paths. Lookup in a table
    /// is logarithmic; coefficient work is linear and checkpointed.
    pub max_chebyshev_coefficients: usize,
}

impl ScanMotionLimit {
    /// Admit the geometric allowance; the source checks its coefficient count.
    pub fn validate(self) -> Result<(), ConductionError> {
        if !self.max_normalized_travel.is_finite()
            || self.max_normalized_travel <= 0.0 || self.max_normalized_travel > 1.0
        { return Err(invalid("scan motion allowance must be a finite half-width fraction in (0,1]")); }
        Ok(())
    }
}

impl MovingHeatSource {
    /// Limit an interval by path breakpoints and footprint-relative speed.
    ///
    /// Position jumps belong to the interval after their knot; an interval
    /// ending there describes the pre-jump field. Linear pieces use their
    /// actual slope. Smooth paths use every coefficient, not endpoint chords.
    /// Outside a clamped path domain the speed is zero. The input end may
    /// already be limited by source power or another participant's events.
    ///
    /// # Errors
    /// Invalid times/policy, coefficient-budget exhaustion, cancellation, or
    /// an unrepresentable speed or positive interval. No source state changes.
    pub fn motion_limited_endpoint(&self, cx: &Cx<'_>, start: f64, end: f64,
        policy: ScanMotionLimit) -> Result<f64, ConductionError>
    {
        self.admit_motion_limit(cx, policy)?;
        if !start.is_finite() || !end.is_finite() || end <= start || !(end-start).is_finite() {
            return Err(invalid("scan motion requires a finite positive represented interval"));
        }
        let mut limited = end;
        let mut speed = 0.0;
        for (axis, signal) in self.center_m.iter().enumerate() {
            poll(cx, axis)?;
            let bound = axis_speed(cx, signal, self.half_width_m[axis], start, &mut limited)?;
            speed = upper(speed + bound)?;
        }
        if speed == 0.0 || (limited-start)*speed <= policy.max_normalized_travel {
            poll(cx, 3)?;
            return Ok(limited);
        }
        let dt = policy.max_normalized_travel / speed;
        limited = limited.min(start + dt);
        if (limited-start)*speed > policy.max_normalized_travel { limited = limited.next_down(); }
        if limited <= start || !limited.is_finite()
            || (limited-start)*speed > policy.max_normalized_travel
        { return Err(invalid("scan motion limit cannot advance the represented physical time")); }
        poll(cx, 3)?;
        Ok(limited)
    }

    pub(in crate::bc::history::march) fn admit_motion_limit(&self, cx: &Cx<'_>,
        policy: ScanMotionLimit) -> Result<(), ConductionError>
    {
        poll(cx, 0)?;
        policy.validate()?;
        let mut count = 0usize;
        for signal in &self.center_m {
            if let TimeSignal::Chebfun(profile) = signal {
                count = count.checked_add(profile.cheb.coeffs().len())
                    .filter(|&n| n <= policy.max_chebyshev_coefficients)
                    .ok_or_else(|| invalid("scan path exceeds the declared Chebyshev coefficient budget"))?;
            }
        }
        Ok(())
    }
}

fn axis_speed(cx: &Cx<'_>, signal: &TimeSignal, radius: f64, start: f64,
    end: &mut f64) -> Result<f64, ConductionError>
{
    let mut stop = |knot: f64| { if knot > start { *end = (*end).min(knot); } };
    match signal {
        TimeSignal::Constant(_) => Ok(0.0),
        TimeSignal::Ramp { t_start, t_end, from, to } => {
            if start < *t_start { stop(*t_start); return Ok(0.0); }
            if start >= *t_end { return Ok(0.0); }
            stop(*t_end);
            scaled_rate((to.value-from.value).abs(), *t_end-*t_start, radius)
        }
        TimeSignal::Table { times, values, interp, .. } => {
            // Construction validated the immutable, nonempty table.
            let next = times.partition_point(|&t| t <= start);
            if next == times.len() { return Ok(0.0); }
            stop(times[next]);
            if next == 0 || *interp == Interp::Hold { return Ok(0.0); }
            scaled_rate((values[next]-values[next-1]).abs(), times[next]-times[next-1], radius)
        }
        TimeSignal::Chebfun(profile) => {
            let (a, b) = profile.cheb.domain();
            if start < a { stop(a); return Ok(0.0); }
            if start >= b { return Ok(0.0); }
            stop(b);
            chebyshev_rate(profile.cheb.coeffs(), b-a, radius, |at| poll(cx, at))
        }
    }
}

fn chebyshev_rate(coefficients: &[f64], duration: f64, radius: f64,
    mut checkpoint: impl FnMut(usize) -> Result<(), ConductionError>) -> Result<f64, ConductionError>
{
    let mut variation = 0.0;
    for (n, &coefficient) in coefficients.iter().enumerate().skip(1) {
        if n % 256 == 0 { checkpoint(n)?; }
        let degree = n as f64;
        variation = upper(variation + upper(coefficient.abs()*degree*degree)?)?;
    }
    // t -> x in [-1,1] contributes 2/duration, not 1/duration.
    upper(scaled_rate(variation, duration, radius)? * 2.0)
}

fn scaled_rate(variation: f64, duration: f64, radius: f64) -> Result<f64, ConductionError> {
    if variation == 0.0 { return Ok(0.0); }
    if !duration.is_finite() || duration <= 0.0 {
        return Err(invalid("scan path duration is not representable"));
    }
    let rate = upper(upper(variation)? / duration)? / radius;
    if rate == 0.0 { return Err(invalid("nonzero scan motion rate underflowed")); }
    upper(rate)
}

fn upper(value: f64) -> Result<f64, ConductionError> {
    if value == 0.0 { return Ok(0.0); }
    let expanded = value.next_up();
    if !expanded.is_finite() || expanded <= 0.0 {
        return Err(invalid("scan path speed bound is not representable"));
    }
    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g1_returning_chebyshev_path_cannot_alias_into_zero_motion() {
        // x^2(1-x^2) = (1-T4(x))/8: zero at both ends and the midpoint.
        let coefficients = [0.125, 0.0, 0.0, 0.0, -0.125];
        let bound = chebyshev_rate(&coefficients, 2.0, 0.25, |_| Ok(())).unwrap();
        assert!(bound >= 8.0 && bound < 8.0*(1.0+1e-13));
        for i in 0..=200 {
            let x = -1.0 + 0.01*f64::from(i);
            assert!((2.0*x-4.0*x*x*x).abs()/0.25 <= bound);
        }
        assert_eq!(chebyshev_rate(&[17.0], 2.0, 0.25, |_| Ok(())).unwrap(), 0.0);
    }
    #[test]
    fn g4_spectral_motion_scans_poll_and_refuse_arithmetic_overflow() {
        let error = ConductionError::Cancelled { stage: "scan-motion-test", at: 256 };
        assert_eq!(chebyshev_rate(&vec![0.01; 300], 1.0, 1.0, |_| Err(error.clone())), Err(error));
        assert!(chebyshev_rate(&[0.0, f64::MAX], 1.0, 1.0, |_| Ok(())).is_err());
        assert!(scaled_rate(1.0, f64::INFINITY, 1.0).is_err());
    }
}

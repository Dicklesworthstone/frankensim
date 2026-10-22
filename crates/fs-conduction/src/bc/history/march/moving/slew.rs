//! Rates of authored continuous power schedules, independent of thermal sampling.
//!
//! A linear segment has dP/dt = (P1-P0)/(t1-t0) throughout its open interval.
//! Clamped tails are constant. Knots can have different one-sided slopes, but
//! cannot hide a jump. Held tables are deliberately unsupported: their finite
//! difference between timestamps is NOT the rate of the actual held signal.

use super::{MovingHeatSource, invalid, poll};
use crate::ConductionError;
use fs_exec::Cx;
use fs_scenario::signal::{Interp, TimeSignal};

/// One continuous ramp in its original physical units and authored ordering.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingHeatPowerSlew {
    pub start_time_s: f64,
    pub end_time_s: f64,
    /// Indices into the source's authored absorbed-power values.
    pub power_value_indices: [usize; 2],
    pub endpoint_powers_w: [f64; 2],
    /// Signed W/s; a negative value is a power decrease, not active cooling.
    pub rate_w_s: f64,
    /// Derivatives of the rate with respect to the two endpoint powers [1/s].
    pub power_value_weights_s_inv: [f64; 2],
}

impl MovingHeatSource {
    /// Enumerate EVERY continuous power segment, including an authored ramp
    /// to an unselected zero-power shutdown. No thermal solve or time sampling
    /// is performed. Constants return an empty vector; clamped tails have zero
    /// rate. The returned slopes also bound the one-sided rates at the knots.
    ///
    /// Supports constant, ramp and linearly interpolated table powers. Held
    /// tables and smooth function objects refuse instead of silently estimating
    /// a continuous rate from samples. `max_segments` bounds retained rows and
    /// traversal. Times are fixed: derivatives concern authored powers only.
    pub fn absorbed_power_slew_segments(&self, cx: &Cx<'_>, max_segments: usize)
        -> Result<Vec<MovingHeatPowerSlew>, ConductionError>
    {
        poll(cx, 0)?;
        let count = match &self.absorbed_power_w {
            TimeSignal::Constant(_) => 0,
            TimeSignal::Ramp { .. } => 1,
            TimeSignal::Table { values, interp: Interp::Linear, .. } => values.len().saturating_sub(1),
            _ => return Err(invalid("finite power slew requires a constant, ramp or linear table; held jumps and smooth profiles are not supported")),
        };
        if count > max_segments {
            return Err(invalid("absorbed-power slew segments exceed their explicit allowance"));
        }
        let mut out = Vec::new();
        out.try_reserve_exact(count).map_err(|_| invalid("power slew allocation refused"))?;
        let mut push = |index: usize, times: [f64; 2], powers: [f64; 2]| -> Result<(), ConductionError> {
            poll(cx, index)?;
            let dt = times[1] - times[0];
            if times.iter().any(|v| !v.is_finite()) || !dt.is_finite() || dt <= 0.0
                || powers.iter().any(|v| !v.is_finite() || *v < 0.0)
            { return Err(invalid("power slew needs ordered finite seconds and finite nonnegative watts")); }
            let inverse = 1.0 / dt;
            let rate = (powers[1] - powers[0]) / dt;
            if !inverse.is_finite() || !rate.is_finite() {
                return Err(invalid("power slew or its endpoint derivatives are not representable"));
            }
            out.push(MovingHeatPowerSlew { start_time_s: times[0], end_time_s: times[1],
                power_value_indices: [index, index + 1], endpoint_powers_w: powers,
                rate_w_s: rate, power_value_weights_s_inv: [-inverse, inverse] });
            Ok(())
        };
        match &self.absorbed_power_w {
            TimeSignal::Constant(value) => {
                if !value.value.is_finite() || value.value < 0.0 {
                    return Err(invalid("power slew requires finite nonnegative watts"));
                }
            }
            TimeSignal::Ramp { t_start, t_end, from, to } =>
                push(0, [*t_start, *t_end], [from.value, to.value])?,
            TimeSignal::Table { times, values, interp: Interp::Linear, .. } => {
                if times.len() != values.len() || values.is_empty() {
                    return Err(invalid("power slew table has inconsistent times and values"));
                }
                if values.len() == 1 && (!values[0].is_finite() || values[0] < 0.0) {
                    return Err(invalid("power slew requires finite nonnegative watts"));
                }
                for i in 0..count { push(i, [times[i], times[i+1]], [values[i], values[i+1]])?; }
            }
            _ => unreachable!("unsupported schedules refused before allocation"),
        }
        poll(cx, count)?;
        Ok(out)
    }
}

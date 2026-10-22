//! Exact segment kinematics of continuous, piecewise-linear source paths.
//!
//! The three position schedules may use different clocks. Their union partitions
//! time into intervals with one constant velocity vector, not three unrelated
//! axis maxima. No thermal sampling, power gating or footprint scaling is used.

use super::{MovingHeatSource, invalid, poll};
use crate::ConductionError;
use fs_exec::Cx;
use fs_scenario::signal::{Interp, TimeSignal};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MovingHeatPathBudget {
    /// Sum of authored position timestamps before deduplication across axes.
    pub max_knots: usize,
    /// Maximum intervals in the merged clock, including stationary intervals.
    pub max_segments: usize,
}

/// Velocity on one active authored axis segment. Constant/clamped axes have
/// no entry: translating their constant value cannot change their velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingHeatAxisRate {
    pub value_indices: [usize; 2],
    pub endpoint_positions_m: [f64; 2],
    /// Duration of this AXIS segment, not the shorter merged-clock interval.
    pub duration_s: f64,
    pub value_weights_s_inv: [f64; 2],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingHeatPathSegment {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub velocity_m_s: [f64; 3],
    pub axes: [Option<MovingHeatAxisRate>; 3],
}

impl MovingHeatSource {
    /// Enumerate the entire continuous path, even when the heater is off.
    /// Supports constants, ramps and linear tables. Held paths and Chebyshev
    /// profiles refuse; a timestamp chord is not their actual velocity.
    ///
    /// Tails outside all authored domains have zero velocity and are omitted.
    /// Velocity may jump at a knot: these are speed data, not acceleration or
    /// jerk bounds. Both one-sided velocities occur in the adjacent segments.
    /// Derivatives hold times fixed and concern the original authored values.
    pub fn path_velocity_segments(&self, cx: &Cx<'_>, budget: MovingHeatPathBudget)
        -> Result<Vec<MovingHeatPathSegment>, ConductionError>
    {
        poll(cx, 0)?;
        let mut count = 0usize;
        for signal in &self.center_m {
            let n = match signal {
                TimeSignal::Constant(_) => 0,
                TimeSignal::Ramp { .. } => 2,
                TimeSignal::Table { times, interp: Interp::Linear, .. } => times.len(),
                _ => return Err(invalid("path speed requires constant, ramp or linear-table positions; jumps and smooth profiles are unsupported")),
            };
            count = count.checked_add(n).filter(|&n| n <= budget.max_knots)
                .ok_or_else(|| invalid("path timestamps exceed their explicit allowance"))?;
        }
        let mut knots = Vec::new();
        knots.try_reserve_exact(count).map_err(|_| invalid("path clock allocation refused"))?;
        for signal in &self.center_m {
            match signal {
                TimeSignal::Constant(value) => {
                    if !value.value.is_finite() { return Err(invalid("path positions must be finite metres")); }
                }
                TimeSignal::Ramp { t_start, t_end, from, to } => {
                    rate(0, [*t_start, *t_end], [from.value, to.value])?;
                    knots.extend_from_slice(&[*t_start, *t_end]);
                }
                TimeSignal::Table { times, values, .. } => {
                    if times.is_empty() || times.len() != values.len() {
                        return Err(invalid("path table needs matching nonempty timestamps and values"));
                    }
                    for (i, (&time, &value)) in times.iter().zip(values).enumerate() {
                        if i % 256 == 0 { poll(cx, i)?; }
                        if !time.is_finite() || !value.is_finite() || (i > 0 && times[i-1] >= time) {
                            return Err(invalid("path table needs finite metres and strictly increasing finite seconds"));
                        }
                        knots.push(time);
                    }
                }
                _ => unreachable!("unsupported paths refused before allocation"),
            }
        }
        // Sorting is bounded by max_knots. Selection itself never samples a
        // midpoint, which could round onto the wrong side of a narrow interval.
        knots.sort_unstable_by(f64::total_cmp);
        knots.dedup_by(|a, b| *a == *b);
        poll(cx, count)?;
        let segments = knots.len().saturating_sub(1);
        if segments > budget.max_segments { return Err(invalid("merged path segments exceed their explicit allowance")); }
        let mut out = Vec::new();
        out.try_reserve_exact(segments).map_err(|_| invalid("path segment allocation refused"))?;
        for (i, pair) in knots.windows(2).enumerate() {
            poll(cx, i)?;
            if !(pair[1] - pair[0]).is_finite() {
                return Err(invalid("merged path duration is not representable"));
            }
            let mut axes = [None; 3];
            let mut velocity = [0.0; 3];
            for axis in 0..3 {
                axes[axis] = active_rate(&self.center_m[axis], pair[0])?;
                if let Some(segment) = axes[axis] {
                    let delta = segment.endpoint_positions_m[1] - segment.endpoint_positions_m[0];
                    velocity[axis] = delta / segment.duration_s;
                    if !velocity[axis].is_finite() || (delta != 0.0 && velocity[axis] == 0.0) {
                        return Err(invalid("path velocity is not representable"));
                    }
                }
            }
            out.push(MovingHeatPathSegment { start_time_s: pair[0], end_time_s: pair[1],
                velocity_m_s: velocity, axes });
        }
        poll(cx, segments)?;
        Ok(out)
    }
}

fn rate(index: usize, times: [f64; 2], values: [f64; 2])
    -> Result<MovingHeatAxisRate, ConductionError>
{
    let duration = times[1] - times[0];
    let inverse = 1.0 / duration;
    if times.iter().chain(&values).any(|x| !x.is_finite())
        || !duration.is_finite() || duration <= 0.0 || !inverse.is_finite() || inverse <= 0.0
    { return Err(invalid("path segment needs representable positive seconds and finite metre coordinates")); }
    Ok(MovingHeatAxisRate { value_indices: [index, index+1], endpoint_positions_m: values,
        duration_s: duration, value_weights_s_inv: [-inverse, inverse] })
}

fn active_rate(signal: &TimeSignal, start: f64) -> Result<Option<MovingHeatAxisRate>, ConductionError> {
    match signal {
        TimeSignal::Constant(_) => Ok(None),
        TimeSignal::Ramp { t_start, t_end, from, to } => {
            if start < *t_start || start >= *t_end { Ok(None) }
            else { rate(0, [*t_start, *t_end], [from.value, to.value]).map(Some) }
        }
        TimeSignal::Table { times, values, .. } => {
            let next = times.partition_point(|&time| time <= start);
            if next == 0 || next == times.len() { Ok(None) }
            else { rate(next-1, [times[next-1], times[next]], [values[next-1], values[next]]).map(Some) }
        }
        _ => unreachable!("unsupported paths refused before clock assembly"),
    }
}

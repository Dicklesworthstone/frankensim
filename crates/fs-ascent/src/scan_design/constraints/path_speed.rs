//! Three-dimensional actuator speed on the original authored source paths.
//!
//! One inequality per merged-clock interval uses (|v|^2 - limit^2)/scale^2.
//! This couples the axes, remains differentiable at rest, and includes fixed
//! return positions and heater-off travel. It is not the timestep motion limiter.

use super::{ConstrainedScanStudy, Limits, LimitAdmission, ScanControlBound,
    ScanEnergyCap, ScanStateBound, ScanPeakBound, ScanPowerSlewBound, ScanResponseBudget,
    ScanConstraintEvaluation, ScanGeometryControl, ScanPowerStudy, ScanObservationSpec,
    ScanSensitivityError, ScanningEnthalpyRun, MovingHeatSource, SqpState, SqpError,
    finite, invalid, poll};
use fs_conduction::bc::history::march::moving::motion::velocity::MovingHeatPathBudget;
use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanPathSpeedBound {
    pub source: usize,
    /// Nonnegative maximum Euclidean speed of the center, metres per second.
    pub maximum_m_s: f64,
    /// Positive speed scale. A scaled tolerance tau permits
    /// |v|^2 <= maximum_m_s^2 + tau * scale_m_s^2, not tau*scale speed error.
    pub scale_m_s: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScanPathSpeedBudget {
    pub max_knots_per_source: usize,
    /// Total retained intervals across all named sources, not per source.
    pub max_segments: usize,
}

/// Independent continuous power and path policies; either group can be empty.
/// Old constructors install no path-speed limits. Existing row/KKT budgets
/// still include every additional interval. Clock times are never variables.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ScanActuatorBounds {
    pub power_slew: Vec<ScanPowerSlewBound>,
    pub max_power_segments: usize,
    pub path_speed: Vec<ScanPathSpeedBound>,
    pub path_budget: ScanPathSpeedBudget,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanPathSpeedValue {
    pub source: usize,
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub velocity_m_s: [f64; 3],
    pub speed_m_s: f64,
}

#[derive(Clone, Copy)]
struct AxisRow {
    fixed_m: [f64; 2],
    duration_s: f64,
    /// Geometry-coordinate index and its metre-per-coordinate derivative.
    selected: [Option<(usize, f64)>; 2],
}
pub(super) struct PathSpeedRow {
    bound: ScanPathSpeedBound,
    times: [f64; 2],
    axes: [Option<AxisRow>; 3],
}

impl<F, O: ScanObservationSpec> ConstrainedScanStudy<F, O> {
    /// Couple vector path-speed and power-rate bounds to the same SQP problem
    /// as thermal responses. Requires continuous piecewise-linear positions on
    /// named sources. Acceleration, jerk, collision and hardware dynamics are
    /// not modeled. A zero speed cap is degenerate for a first-order KKT test;
    /// leave stationary path coordinates unselected when possible.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_actuator_bounds<'h, 'm, 'c>(cx: &Cx<'_>, mut study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>,
        peaks: Vec<ScanPeakBound>, response_budget: ScanResponseBudget, max_peak_samples: usize,
        actuators: ScanActuatorBounds, max_kkt_dimension: usize)
        -> Result<Self, SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        let limits = Limits::compile_with_actuators(cx, &study, bounds, energy, states, peaks,
            response_budget, max_peak_samples, actuators,
            LimitAdmission::Differentiated(max_kkt_dimension)).map_err(SqpError::Evaluation)?;
        let initial = study.initial_coordinates().to_vec();
        let state = SqpState::try_new(&initial, max_kkt_dimension,
            &mut |point| limits.sample(cx, &mut study, point), Some(cx))?;
        Ok(Self { study, limits, state })
    }
}

pub(super) fn compile<F, O: ScanObservationSpec>(cx: &Cx<'_>, study: &ScanPowerStudy<F, O>,
    bounds: &[ScanPathSpeedBound], budget: ScanPathSpeedBudget)
    -> Result<Vec<PathSpeedRow>, ScanSensitivityError>
{
    poll(cx)?;
    if bounds.len() > study.sources.len() { return Err(invalid("path speed needs at most one bound per source")); }
    let mut selected = BTreeMap::new();
    for (i, &control) in study.geometry_controls.iter().enumerate() {
        poll(cx)?;
        if let ScanGeometryControl::CenterValue { source, axis, value, scale_m, .. } = control {
            selected.insert((source, axis, value), (i, scale_m));
        }
    }
    let mut seen = BTreeSet::new();
    let mut rows = Vec::new();
    for &bound in bounds {
        poll(cx)?;
        if !seen.insert(bound.source) || !bound.maximum_m_s.is_finite() || bound.maximum_m_s < 0.0
            || !bound.scale_m_s.is_finite() || bound.scale_m_s <= 0.0
        { return Err(invalid("path speed needs distinct sources, nonnegative finite maxima and positive finite metre-per-second scales")); }
        normalized_square(bound.maximum_m_s, bound.scale_m_s)?;
        let source = study.sources.get(bound.source).ok_or_else(|| invalid("path speed names an absent source"))?;
        let segments = source.path_velocity_segments(cx, MovingHeatPathBudget {
            max_knots: budget.max_knots_per_source, max_segments: budget.max_segments - rows.len() })?;
        rows.try_reserve_exact(segments.len()).map_err(|_| invalid("path speed row allocation refused"))?;
        for segment in segments {
            poll(cx)?;
            rows.push(PathSpeedRow { bound, times: [segment.start_time_s, segment.end_time_s],
                axes: std::array::from_fn(|axis| segment.axes[axis].map(|rate| AxisRow {
                    fixed_m: rate.endpoint_positions_m, duration_s: rate.duration_s,
                    selected: rate.value_indices.map(|i| selected.get(&(bound.source, axis, i)).copied()),
                })) });
        }
    }
    poll(cx)?;
    Ok(rows)
}

fn normalized_square(value: f64, scale: f64) -> Result<f64, ScanSensitivityError> {
    let normalized = finite(value / scale)?;
    let square = finite(normalized * normalized)?;
    if value != 0.0 && square == 0.0 { return Err(invalid("scaled path speed squared underflowed")); }
    Ok(square)
}

pub(super) fn append<F, O: ScanObservationSpec>(cx: &Cx<'_>, study: &ScanPowerStudy<F, O>,
    rows: &[PathSpeedRow], geometry: &[f64], derivatives: bool, out: &mut ScanConstraintEvaluation)
    -> Result<(), ScanSensitivityError>
{
    let n = study.initial_coordinates().len();
    for row in rows {
        poll(cx)?;
        let mut velocity = [0.0; 3];
        let mut square = 0.0;
        for axis in 0..3 {
            if let Some(data) = row.axes[axis] {
                let positions: [f64; 2] = std::array::from_fn(|side|
                    data.selected[side].map_or(data.fixed_m[side], |(i, _)| geometry[i]));
                let delta = finite(positions[1] - positions[0])?;
                velocity[axis] = finite(delta / data.duration_s)?;
                if delta != 0.0 && velocity[axis] == 0.0 { return Err(invalid("path velocity underflowed")); }
            }
            square = finite(square + normalized_square(velocity[axis], row.bound.scale_m_s)?)?;
        }
        let maximum = velocity.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        let speed = if maximum == 0.0 { 0.0 } else {
            finite(maximum * velocity.iter().map(|v| (v / maximum).powi(2)).sum::<f64>().sqrt())?
        };
        if derivatives {
            let offset = out.residuals.len() * n;
            for (axis, data) in row.axes.iter().enumerate() {
                if let Some(data) = data {
                    for (side, sign) in [(0, -1.0), (1, 1.0)] {
                        if let Some((index, scale_m)) = data.selected[side] {
                            let slope = finite(finite(scale_m / data.duration_s)? / row.bound.scale_m_s)?;
                            let gradient = finite(2.0 * finite(velocity[axis] / row.bound.scale_m_s)? * sign * slope)?;
                            let slot = offset + study.controls.len() + index;
                            out.jacobian[slot] = finite(out.jacobian[slot] + gradient)?;
                        }
                    }
                }
            }
        }
        out.residuals.push(finite(square - normalized_square(row.bound.maximum_m_s, row.bound.scale_m_s)?)?);
        out.path_speed_values.push(ScanPathSpeedValue { source: row.bound.source,
            start_time_s: row.times[0], end_time_s: row.times[1], velocity_m_s: velocity, speed_m_s: speed });
    }
    poll(cx)?;
    Ok(())
}

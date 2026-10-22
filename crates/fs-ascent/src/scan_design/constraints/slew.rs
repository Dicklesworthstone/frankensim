//! Continuous absorbed-power rate inequalities on the authored schedules.
//!
//! Each original segment contributes (rate-rise)/scale <= 0 and
//! (-rate-fall)/scale <= 0. The rates are not estimated from thermal endpoints.
//! Frozen values, including zero shutdown values, stay in the physical rate;
//! only selected powers contribute derivatives through their actual coordinates.

use super::{ConstrainedScanStudy, Limits, LimitAdmission, ScanControlBound,
    ScanEnergyCap, ScanStateBound, ScanPeakBound, ScanResponseBudget,
    ScanConstraintEvaluation, ScanPowerStudy, ScanObservationSpec, ScanSensitivityError,
    ScanningEnthalpyRun, MovingHeatSource, SqpState, SqpError, finite, invalid, poll};
use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

/// Per-source continuous rate limits over its entire authored power schedule.
/// Constant/ramp/linear-table powers only. This limits absorbed power, not the
/// electrical input, optical dynamics or the mechanical speed of the source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanPowerSlewBound {
    pub source: usize,
    /// Maximum positive slope, W/s. Zero prohibits increases.
    pub maximum_rise_w_s: f64,
    /// Maximum magnitude of a negative slope, W/s. Zero prohibits decreases.
    pub maximum_fall_w_s: f64,
    /// Positive W/s scale; tolerance tau permits tau*scale_w_s violation.
    pub scale_w_s: f64,
}

pub(super) struct PowerSlewRow {
    bound: ScanPowerSlewBound,
    duration_s: f64,
    fixed_powers_w: [f64; 2],
    selected: [Option<usize>; 2],
}

impl<F, O: ScanObservationSpec> ConstrainedScanStudy<F, O> {
    /// Add actuator-rate inequalities to the original SQP model. All rows and
    /// their coordinate derivatives enter the same KKT problem as the thermal
    /// requirements; no after-the-fact clipping or surrogate is introduced.
    /// Original constructors retain their behavior when no rates are declared.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_power_slew_bounds<'h, 'm, 'c>(cx: &Cx<'_>, mut study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>,
        peaks: Vec<ScanPeakBound>, response_budget: ScanResponseBudget, max_peak_samples: usize,
        slew_bounds: Vec<ScanPowerSlewBound>, max_slew_segments: usize, max_kkt_dimension: usize)
        -> Result<Self, SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        let limits = Limits::compile_with_slew(cx, &study, bounds, energy, states, peaks,
            response_budget, max_peak_samples, slew_bounds, max_slew_segments,
            LimitAdmission::Differentiated(max_kkt_dimension)).map_err(SqpError::Evaluation)?;
        let initial = study.initial_coordinates().to_vec();
        let state = SqpState::try_new(&initial, max_kkt_dimension,
            &mut |point| limits.sample(cx, &mut study, point), Some(cx))?;
        Ok(Self { study, limits, state })
    }
}

pub(super) fn compile<F, O: ScanObservationSpec>(cx: &Cx<'_>, study: &ScanPowerStudy<F, O>,
    bounds: &[ScanPowerSlewBound], max_segments: usize)
    -> Result<Vec<PowerSlewRow>, ScanSensitivityError>
{
    poll(cx)?;
    if bounds.len() > study.sources.len() {
        return Err(invalid("power slew needs at most one declaration per existing source"));
    }
    let selected: BTreeMap<_, _> = study.controls.iter().enumerate()
        .map(|(i, c)| ((c.source, c.value), i)).collect();
    let mut seen = BTreeSet::new();
    let mut rows = Vec::new();
    for &bound in bounds {
        poll(cx)?;
        if !seen.insert(bound.source) || !bound.scale_w_s.is_finite() || bound.scale_w_s <= 0.0
            || [bound.maximum_rise_w_s, bound.maximum_fall_w_s].iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        { return Err(invalid("power slew needs distinct sources, nonnegative finite rates and positive finite W/s scales")); }
        let source = study.sources.get(bound.source)
            .ok_or_else(|| invalid("power slew names an absent source"))?;
        let segments = source.absorbed_power_slew_segments(cx, max_segments - rows.len())?;
        rows.try_reserve_exact(segments.len()).map_err(|_| invalid("power slew row allocation refused"))?;
        for segment in segments {
            poll(cx)?;
            rows.push(PowerSlewRow { bound,
                duration_s: segment.end_time_s - segment.start_time_s,
                fixed_powers_w: segment.endpoint_powers_w,
                selected: segment.power_value_indices.map(|i| selected.get(&(bound.source, i)).copied()) });
        }
    }
    poll(cx)?;
    Ok(rows)
}

pub(super) fn append<F, O: ScanObservationSpec>(cx: &Cx<'_>, study: &ScanPowerStudy<F, O>,
    rows: &[PowerSlewRow], powers: &[f64], derivatives: bool, out: &mut ScanConstraintEvaluation)
    -> Result<(), ScanSensitivityError>
{
    let n = study.initial_coordinates().len();
    for row in rows {
        poll(cx)?;
        let endpoints: [f64; 2] = std::array::from_fn(|side|
            row.selected[side].map_or(row.fixed_powers_w[side], |i| powers[i]));
        let rate = finite(finite(endpoints[1] - endpoints[0])? / row.duration_s)?;
        out.power_slew_values_w_s.push(rate);
        for (sign, limit) in [(1.0, row.bound.maximum_rise_w_s), (-1.0, row.bound.maximum_fall_w_s)] {
            if derivatives {
                let offset = out.residuals.len() * n;
                for (side, endpoint_sign) in [(0, -1.0), (1, 1.0)] {
                    if let Some(index) = row.selected[side] {
                        let derivative = finite(study.power_derivative(index, powers[index]) / row.duration_s)?;
                        out.jacobian[offset + index] = finite(sign * endpoint_sign * (derivative / row.bound.scale_w_s))?;
                    }
                }
            }
            out.residuals.push(finite(finite(sign * rate - limit)? / row.bound.scale_w_s)?);
        }
    }
    poll(cx)?;
    Ok(())
}

//! Opt-in event-aligned backward Euler, one accepted interval per call.
//!
//! A right-continuous held load belongs to the interval AFTER its knot. Merely
//! landing on a heater's turn-off knot and evaluating the new zero value would
//! erase the preceding pulse. This path samples the exact left limit at the
//! interval end, without subtracting an epsilon from physical time. Continuous
//! histories keep their original endpoint interpolation. Every supplied table
//! knot, ramp endpoint and clamped Chebyshev-domain endpoint is visited.
//!
//! This integrates held source work exactly up to arithmetic/solver error, but
//! does NOT integrate temperature-dependent boundary exchange exactly. Smooth
//! forcing remains backward-Euler quadrature. A prescribed-temperature jump
//! takes effect in the following finite step, with the original boundary-node
//! capacity reaction; no instantaneous impulse is invented. At a jump, accepted
//! fields and reports describe the interval-ending, pre-jump state.

use super::{HeatSourceHistory, ThermalHistoryMethod, ThermalHistoryRun,
    ThermalHistoryStep, invalid, poll};
use super::super::ThermalBoundaryHistory;
use crate::{ConductionError, ConductionProblem, ScalarField, ThermalBc,
    ThermalBoundary, ThermalBoundaryBuilder, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS};
use crate::transient::backward_euler::StepConfig;
use fs_exec::Cx;
use fs_qty::{Dims, QtyAny};
use fs_scenario::bc::{BcValue, BoundaryCondition};
use fs_scenario::signal::{Interp, TimeSignal};

/// One completed interval, which may stop before the requested target.
#[derive(Debug, Clone)]
pub struct HistoryEventStep<T> {
    /// Previously accepted physical time, seconds.
    pub start_time_s: f64,
    /// Next schedule knot, step-size limit or requested target, seconds.
    pub end_time_s: f64,
    /// The unchanged physical solver's complete response for this interval.
    pub report: T,
}

impl ThermalHistoryRun<'_, '_> {
    /// Advance ONE step, stopping at the next schedule knot or `max_step_s`.
    ///
    /// Repeat until `time_s() == target_s` to traverse a window. Each successful
    /// call consumes one existing lifetime-budget slot. Source and boundary
    /// histories use the exact interval-end left limit; [`Self::advance_to`]
    /// retains its original right-continuous endpoint behavior. No table is
    /// rescanned to build or sort a global event list. Cancellation is checked
    /// between histories; signal evaluation retains its existing scan behavior.
    ///
    /// # Errors
    /// Invalid/unrepresentable times, exhausted budget, or the original history,
    /// material, cancellation, residual or energy refusal. No failed call changes
    /// accepted state. Previously completed calls remain available for resume.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: StepConfig)
        -> Result<HistoryEventStep<ThermalHistoryStep>, ConductionError>
    {
        poll(cx, self.accepted_steps)?;
        if self.accepted_steps >= self.max_steps {
            return Err(invalid("thermal history step budget exhausted"));
        }
        let start_time_s = self.time_s;
        let end_time_s = next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start_time_s, target_s, max_step_s)?;
        let boundary = boundary_before(cx, self.problem.boundary, end_time_s)?;
        let source = source_before(&self.problem.source, end_time_s, self.temperature.len())?;
        let problem = ConductionProblem { mesh: self.problem.boundary.mesh(),
            boundary: &boundary, material: self.problem.material, source: &source,
            element_materials: self.problem.element_materials };
        let dt = end_time_s - start_time_s;
        let report = match self.method {
            ThermalHistoryMethod::Linear => ThermalHistoryStep::Linear(
                self.stepper.advance_prescribed(cx, problem, self.problem.interfaces,
                    &self.temperature, dt, config)?),
            ThermalHistoryMethod::Nonlinear(policy) => ThermalHistoryStep::Nonlinear(
                self.stepper.advance_nonlinear_prescribed(cx, problem, self.problem.interfaces,
                    &self.temperature, dt, config, policy)?),
        };
        self.publish(cx, end_time_s, report.solution())?;
        Ok(HistoryEventStep { start_time_s, end_time_s, report })
    }
}

// All histories are immutable after their original constructor admission.
// Table lookup is logarithmic; no event-count-sized allocation or epsilon merge.
pub(super) fn next_endpoint<'s>(cx: &Cx<'_>, history: &ThermalBoundaryHistory<'_>,
    source: &HeatSourceHistory, extra: impl IntoIterator<Item = &'s TimeSignal>,
    start: f64, target: f64, max_step: f64) -> Result<f64, ConductionError>
{
    poll(cx, 0)?;
    let remaining = target - start;
    if !start.is_finite() || !target.is_finite() || !remaining.is_finite()
        || remaining <= 0.0 || !max_step.is_finite() || max_step <= 0.0
    { return Err(invalid("event stepping requires a later finite target and a positive finite maximum step")); }
    let mut end = if remaining <= max_step { target } else { start + max_step };
    // Addition can round upward: the actual solver interval must respect the cap.
    if end - start > max_step { end = end.next_down(); }
    if !end.is_finite() || end <= start {
        return Err(invalid("event step is too small to advance the represented physical time"));
    }
    for (at, (_, row)) in history.rows.iter().enumerate() {
        poll(cx, at)?;
        if let Some(BcValue::Signal(signal)) = &row.row.value { limit(signal, start, &mut end); }
        if let Some(signal) = &row.reference_temperature { limit(signal, start, &mut end); }
    }
    match source {
        HeatSourceHistory::Fixed(_) => {}
        HeatSourceHistory::Uniform(signal) | HeatSourceHistory::Scaled { scale: signal, .. } => {
            limit(signal, start, &mut end);
        }
    }
    for (at, signal) in extra.into_iter().enumerate() {
        poll(cx, at)?;
        limit(signal, start, &mut end);
    }
    poll(cx, 0)?;
    Ok(end)
}

fn limit(signal: &TimeSignal, start: f64, end: &mut f64) {
    let mut consider = |time: f64| { if time > start && time < *end { *end = time; } };
    match signal {
        TimeSignal::Constant(_) => {}
        TimeSignal::Ramp { t_start, t_end, .. } => { consider(*t_start); consider(*t_end); }
        TimeSignal::Table { times, .. } => {
            if let Some(&time) = times.get(times.partition_point(|&t| t <= start)) { consider(time); }
        }
        TimeSignal::Chebfun(profile) => {
            let (a, b) = profile.cheb.domain();
            consider(a); consider(b);
        }
    }
}

pub(super) fn value_before(signal: &TimeSignal, time_s: f64) -> Result<QtyAny, ConductionError> {
    // The original owner still validates every supplied value and dimension.
    let value = signal.eval(time_s).map_err(|e| invalid(&format!("event history: {e}")))?;
    if let TimeSignal::Table { times, values, dims, interp: Interp::Hold } = signal {
        let index = times.partition_point(|&t| t < time_s);
        if index > 0 && times.get(index).is_some_and(|&t| t == time_s) {
            return Ok(QtyAny::new(values[index - 1], *dims));
        }
    }
    Ok(value)
}

fn quantity_before(signal: &TimeSignal, time: f64, dims: Dims) -> Result<f64, ConductionError> {
    let value = value_before(signal, time)?;
    if value.dims != dims {
        return Err(ConductionError::Dimensions { context: "event history".into(),
            expected: dims.0, found: value.dims.0 });
    }
    Ok(value.value)
}

pub(super) fn source_before(source: &HeatSourceHistory, time: f64, vertices: usize)
    -> Result<ScalarField, ConductionError>
{
    let sampled = match source {
        HeatSourceHistory::Fixed(field) => field.clone(),
        HeatSourceHistory::Uniform(signal) => {
            ScalarField::Uniform(quantity_before(signal, time, VOLUMETRIC_SOURCE_DIMS)?)
        }
        HeatSourceHistory::Scaled { field, scale } => {
            field.validate("heat source footprint", vertices)?;
            let factor = quantity_before(scale, time, Dims([0; 6]))?;
            match field {
                ScalarField::Uniform(v) => ScalarField::Uniform(v * factor),
                ScalarField::Nodal(v) => ScalarField::Nodal(v.iter().map(|v| v * factor).collect()),
            }
        }
    };
    sampled.validate("event heat source", vertices)?;
    Ok(sampled)
}

pub(super) fn boundary_before(cx: &Cx<'_>, history: &ThermalBoundaryHistory<'_>, time: f64)
    -> Result<ThermalBoundary, ConductionError>
{
    poll(cx, 0)?;
    let mut conditions = history.partition.conditions.clone();
    for (at, (index, scheduled)) in history.rows.iter().enumerate() {
        poll(cx, at)?;
        let reference = scheduled.reference_temperature.as_ref()
            .map(|s| quantity_before(s, time, TEMPERATURE_DIMS)).transpose()?;
        let row = &scheduled.row;
        let value = match &row.value {
            Some(BcValue::Signal(signal)) => Some(BcValue::Uniform(value_before(signal, time)?)),
            other => other.clone(),
        };
        let sampled = BoundaryCondition { region: row.region.clone(), physics: row.physics,
            kind: row.kind, value, compatibility: row.compatibility, frame: row.frame };
        conditions[*index] = ThermalBc::from_scenario_row_at(&sampled, &history.environment, reference, time)?;
    }
    let boundary = ThermalBoundaryBuilder { mesh: history.mesh(),
        names: history.partition.names.clone(), conditions,
        face_region: history.partition.face_region.clone(),
        adiabatic_remainder: history.partition.adiabatic_remainder != 0 }.finish()?;
    poll(cx, history.rows.len())?;
    Ok(boundary)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_adjacent_knots_are_distinct_and_left_limits_do_not_shift_continuous_signals() {
        let a = 1.0_f64;
        let b = a.next_up();
        let signal = TimeSignal::Table { times: vec![0.0, a, b], values: vec![0.0, 3.0, 7.0],
            dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold };
        let mut end = 2.0;
        limit(&signal, a, &mut end);
        assert_eq!(end, b);
        assert_eq!(value_before(&signal, a).unwrap().value, 0.0);
        assert_eq!(value_before(&signal, b).unwrap().value, 3.0);
        assert_eq!(value_before(&signal, b.next_up()).unwrap().value, 7.0);
        let ramp = TimeSignal::Ramp { t_start: 0.0, t_end: 2.0,
            from: QtyAny::new(0.0, TEMPERATURE_DIMS), to: QtyAny::new(20.0, TEMPERATURE_DIMS) };
        assert_eq!(value_before(&ramp, 1.0).unwrap(), ramp.eval(1.0).unwrap());
        assert_eq!(value_before(&ramp, 3.0).unwrap().value, 20.0);
        let malformed = TimeSignal::Table { times: vec![0.0, 1.0], values: vec![3.0],
            dims: TEMPERATURE_DIMS, interp: Interp::Hold };
        assert!(value_before(&malformed, 1.0).is_err());
    }
}

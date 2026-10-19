//! Least-squares inverse objectives over a complete fixed-grid heating history.
//!
//! Each endpoint runs the original nonlinear scan producer once. The reverse
//! sweep propagates the downstream enthalpy seed BEFORE adding the current
//! observation's seed, so early heating and intermediate measurements are not
//! lost. All sampled source gradients are contracted into the original authored
//! constant/ramp/table/Chebyshev values, on the same side of every event.
//!
//! This is a bounded full-storage evaluation, not checkpoint recomputation or
//! differentiation of adaptive grids, event times, geometry or material charts.

use super::{ScanSensitivityBudget, ScanSensitivityError, ScanStepLinearization};
use super::super::{ScanStep, ScanningEnthalpyRun};
use crate::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient;
use crate::{ConductionError, LinearConfig};
use crate::transient::enthalpy::EnthalpyStepConfig;
use fs_exec::Cx;

/// Physical quantity measured at a retained mesh vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanObservable {
    /// Absolute temperature [K], pulled back through the original chart dT/dh.
    Temperature,
    /// Specific enthalpy [J/kg], retaining sensitivity on latent plateaus.
    SpecificEnthalpy,
}

/// One scaled observation; duplicates add independent residual terms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanObservation {
    /// Zero-based index in the supplied endpoint array (not an initial sample).
    pub endpoint: usize,
    /// Vertex of the original, unchanged conduction mesh.
    pub vertex: usize,
    /// Selects the units of both target and scale.
    pub observable: ScanObservable,
    /// Target in kelvin or J/kg as selected above.
    pub target: f64,
    /// Positive residual scale in the same units. It need not be an uncertainty.
    pub scale: f64,
}

/// A fixed quadrature convention, held unchanged by every derivative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanObjectiveSampling {
    /// Same right-continuous endpoint rule as advance_linearized_to.
    RightEndpoint,
    /// Pre-switch values, requiring the supplied grid to include every event
    /// and selected motion stop. The evaluator refuses to alter the grid.
    EventLeft,
}

/// Total full-storage admission, separate from per-solve residual/work budgets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanObjectiveBudget {
    pub max_steps: usize,
    pub max_observations: usize,
    /// Endpoints times mesh vertices. Reports/tangents retain several arrays;
    /// this is a shape count, not a peak allocation or byte bound.
    pub max_retained_vertices: usize,
    /// Endpoints times sources times mesh vertices.
    pub max_retained_source_vertices: usize,
    /// Sum over sources of authored signal values plus their three half widths.
    pub max_schedule_controls: usize,
    /// Total radiation-feedback allowance divided among all retained endpoints.
    pub max_retained_radiation_entries: usize,
}

/// Complete objective and gradient, returned only after the reverse sweep passes.
#[derive(Debug)]
pub struct ScanObjectiveEvaluation {
    /// Dimensionless 0.5 * sum((prediction-target)/scale)^2.
    pub value: f64,
    /// Original observation order, in the selected physical units.
    pub predictions: Vec<f64>,
    pub standardized_residuals: Vec<f64>,
    /// Gradients in each source's original schedule representation and SI units.
    pub sources: Vec<MovingHeatScheduleGradient>,
    /// Sensitivity to the input run's initial specific enthalpy, in vertex order.
    pub initial_specific_enthalpy: Vec<f64>,
    /// Original last endpoint report, not a replacement or smoothed trajectory.
    pub final_step: ScanStep,
    /// Actual primal endpoints retained; solver iterations remain in their reports.
    pub endpoint_count: usize,
    /// Sum of reported iterations of the original physical transpose solves.
    pub transpose_iterations: usize,
    pub max_transpose_relative_residual: f64,
}

/// Evaluate a complete scan objective and its authored-control gradients.
///
/// The trial run is consumed: failed candidates cannot mutate any caller-owned
/// accepted run. Build a fresh run for each optimizer point. Source declarations
/// never change during evaluation, so sensitivities cannot be rebound to a
/// different path, power schedule, footprint, event side or material history.
/// Both nonlinear conductivity and attached radiation use the existing exact
/// endpoint tangents. Background/boundary data remain fixed design inputs.
///
/// Retains at most the explicitly admitted endpoint count and source-node shape
/// count; each endpoint's sparse tangent follows the original mesh and interface
/// topology. No global memory/continuum-error certificate is claimed. Classical
/// derivatives at phase/material corners are refused by the original owner.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_objective(
    cx: &Cx<'_>, mut run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
) -> Result<ScanObjectiveEvaluation, ScanSensitivityError> {
    poll(cx, 0)?;
    let count = endpoints.len();
    let vertices = run.conduction().temperature().len();
    let remaining = run.conduction().max_steps().checked_sub(run.conduction().accepted_steps());
    if count == 0 || count > budget.max_steps || remaining.is_none_or(|n| count > n)
        || observations.is_empty() || observations.len() > budget.max_observations
    { return Err(invalid("nonempty endpoints/observations must fit the objective and remaining physical-step budgets")); }
    let retained = count.checked_mul(vertices)
        .filter(|&n| n <= budget.max_retained_vertices)
        .ok_or_else(|| invalid("scan objective retained vertex allowance exceeded"))?;
    if retained.checked_mul(run.sources().len()).is_none_or(|n| n > budget.max_retained_source_vertices) {
        return Err(invalid("scan objective retained source-vertex allowance exceeded"));
    }
    let mut controls = 0usize;
    for (i, source) in run.sources().iter().enumerate() {
        poll(cx, i)?;
        controls = controls.checked_add(source.schedule_control_count()?)
            .filter(|&n| n <= budget.max_schedule_controls)
            .ok_or_else(|| invalid("scan objective authored-control allowance exceeded"))?;
    }
    let mut previous = run.conduction().time_s();
    for (i, &time) in endpoints.iter().enumerate() {
        poll(cx, i)?;
        if !time.is_finite() || !(time-previous).is_finite() || time <= previous {
            return Err(invalid("scan objective endpoints must increase from the input state with finite positive intervals"));
        }
        previous = time;
    }
    let mut at_endpoint = Vec::new();
    at_endpoint.try_reserve_exact(count).map_err(|_| invalid("scan observation index allocation refused"))?;
    at_endpoint.resize_with(count, Vec::new);
    for (i, observation) in observations.iter().enumerate() {
        poll(cx, i)?;
        if observation.endpoint >= count || observation.vertex >= vertices
            || !observation.target.is_finite() || !observation.scale.is_finite() || observation.scale <= 0.0
        { return Err(invalid("scan observation needs an existing endpoint/vertex, finite target and positive finite physical scale")); }
        at_endpoint[observation.endpoint].push(i);
    }
    let per_step = ScanSensitivityBudget {
        max_source_vertices: budget.max_retained_source_vertices / count,
        max_radiation_feedback_entries: budget.max_retained_radiation_entries / count,
    };
    let mut tape: Vec<ScanStepLinearization<'_>> = Vec::new();
    tape.try_reserve_exact(count).map_err(|_| invalid("scan objective endpoint allocation refused"))?;
    let mut predictions = zeros(observations.len())?;
    let mut residuals = zeros(observations.len())?;
    let mut value = 0.0;
    for (index, &time) in endpoints.iter().enumerate() {
        poll(cx, index)?;
        let binding = match sampling {
            ScanObjectiveSampling::RightEndpoint => run.advance_linearized_to(cx, time, step, per_step)?,
            ScanObjectiveSampling::EventLeft => run.advance_linearized_event_step(
                cx, time, time-run.conduction().time_s(), step, per_step)?,
        };
        if binding.end_time_s() != time {
            return Err(invalid("fixed scan objective grid omits a source/boundary event or selected motion stop"));
        }
        for &i in &at_endpoint[index] {
            poll(cx, i)?;
            let observation = observations[i];
            let solution = binding.report().solution();
            let prediction = match observation.observable {
                ScanObservable::Temperature => solution.temperature[observation.vertex],
                ScanObservable::SpecificEnthalpy => solution.specific_enthalpy_j_kg[observation.vertex],
            };
            let residual = finite((prediction-observation.target)/observation.scale)?;
            value = finite(value + 0.5 * residual * residual)?;
            predictions[i] = prediction;
            residuals[i] = residual;
        }
        tape.push(binding);
    }
    let mut seed = zeros(vertices)?;
    let mut accumulated: Option<Vec<MovingHeatScheduleGradient>> = None;
    let mut iterations = 0usize;
    let mut residual = 0.0_f64;
    for index in (0..count).rev() {
        poll(cx, index)?;
        let binding = &tape[index];
        let mut temperature_seed = zeros(vertices)?;
        for &i in &at_endpoint[index] {
            poll(cx, i)?;
            let observation = observations[i];
            let local = finite(residuals[i]/observation.scale)?;
            let target = match observation.observable {
                ScanObservable::Temperature => &mut temperature_seed[observation.vertex],
                ScanObservable::SpecificEnthalpy => &mut seed[observation.vertex],
            };
            *target = finite(*target + local)?;
        }
        let local_h = binding.transport().transport().temperature_pullback(cx, &temperature_seed)?;
        add(cx, &mut seed, &local_h)?;
        let gradient = binding.pullback(cx, &seed, transpose)?;
        iterations = iterations.checked_add(gradient.storage.transport.iterations)
            .ok_or_else(|| invalid("scan objective transpose work count overflow"))?;
        residual = residual.max(gradient.storage.transport.relative_residual);
        let mut schedules = Vec::with_capacity(gradient.sources.len());
        for (i, (source, gradient)) in run.sources().iter().zip(gradient.sources).enumerate() {
            poll(cx, i)?;
            schedules.push(match sampling {
                ScanObjectiveSampling::RightEndpoint => source.pullback_schedule_at(
                    cx, binding.end_time_s(), gradient, budget.max_schedule_controls)?,
                ScanObjectiveSampling::EventLeft => source.pullback_schedule_before(
                    cx, binding.end_time_s(), gradient, budget.max_schedule_controls)?,
            });
        }
        match &mut accumulated {
            None => accumulated = Some(schedules),
            Some(total) => {
                if total.len() != schedules.len() { return Err(invalid("scan source identity changed during reverse sweep")); }
                for (a, b) in total.iter_mut().zip(&schedules) {
                    for axis in 0..3 { add(cx, &mut a.center_values_m[axis], &b.center_values_m[axis])?; }
                    add(cx, &mut a.half_width_m, &b.half_width_m)?;
                    add(cx, &mut a.absorbed_power_values_w, &b.absorbed_power_values_w)?;
                }
            }
        }
        seed = gradient.storage.transport.previous_specific_enthalpy;
    }
    let final_step = tape.last().ok_or_else(|| invalid("missing scan endpoint"))?.report().clone();
    poll(cx, count)?;
    Ok(ScanObjectiveEvaluation { value, predictions, standardized_residuals: residuals,
        sources: accumulated.ok_or_else(|| invalid("missing scan reverse sweep"))?,
        initial_specific_enthalpy: seed, final_step, endpoint_count: count,
        transpose_iterations: iterations, max_transpose_relative_residual: residual })
}

fn add(cx: &Cx<'_>, target: &mut [f64], values: &[f64]) -> Result<(), ScanSensitivityError> {
    if target.len() != values.len() { return Err(invalid("scan objective gradient shape changed")); }
    for (i, (a, b)) in target.iter_mut().zip(values).enumerate() {
        if i % 256 == 0 { poll(cx, i)?; }
        *a = finite(*a + *b)?;
    }
    Ok(())
}
fn zeros(n: usize) -> Result<Vec<f64>, ScanSensitivityError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| invalid("scan objective vector allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}
fn finite(value: f64) -> Result<f64, ScanSensitivityError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("scan objective arithmetic is not representable")) }
}
fn invalid(what: &str) -> ScanSensitivityError {
    ConductionError::Config { parameter: "scan history objective", what: what.into() }.into()
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ScanSensitivityError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "scan-history-objective", at }.into())
}

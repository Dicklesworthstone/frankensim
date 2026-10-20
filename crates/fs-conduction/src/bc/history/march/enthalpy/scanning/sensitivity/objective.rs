//! Least-squares inverse objectives over a complete fixed-grid heating history.
//!
//! Each endpoint runs the original nonlinear scan producer once. The reverse
//! sweep propagates the downstream enthalpy seed BEFORE adding the current
//! observation's seed, so early heating and intermediate measurements are not
//! lost. All sampled source gradients are contracted into the original authored
//! constant/ramp/table/Chebyshev values, on the same side of every event.
//!
//! The default evaluator retains its full bounded tape; `checkpoint` adds an
//! explicit segmented replay policy. Neither differentiates adaptive grids,
//! event times, geometry or material charts.

/// Independent physical response rows sharing the original forward trajectory.
pub mod responses;
/// Short-tape reverse sweeps with bounded caloric-state checkpoints.
pub mod checkpoint;
/// Mesh-bound interior probes and finite spatial readouts.
pub mod spatial;

use super::{ScanSensitivityBudget, ScanSensitivityError, ScanStepLinearization};
use super::super::{ScanStep, ScanningEnthalpyRun};
use crate::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient;
use crate::{ConductionError, LinearConfig};
use crate::transient::enthalpy::EnthalpyStepConfig;
use fs_exec::Cx;
use spatial::ScanObservationSpec;

/// Physical quantity measured at a retained mesh vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanObservable {
    /// Absolute temperature [K], pulled back through the original chart dT/dh.
    Temperature,
    /// Specific enthalpy [J/kg], retaining sensitivity on latent plateaus.
    SpecificEnthalpy,
    /// Dimensionless liquid MASS fraction in [0,1], from the original per-node
    /// chart. Not percent, volume fraction or an inferred temperature threshold.
    LiquidMassFraction,
}

impl ScanObservable {
    // Only newly injected observation seeds are converted. Derivatives carried
    // from later timesteps already use enthalpy coordinates, even on a plateau.
    fn seed_in_enthalpy(self, cx: &Cx<'_>, binding: &ScanStepLinearization<'_>, seed: Vec<f64>)
        -> Result<Vec<f64>, ScanSensitivityError>
    {
        match self {
            Self::Temperature => Ok(binding.transport().transport().temperature_pullback(cx, &seed)?),
            Self::SpecificEnthalpy => Ok(seed),
            Self::LiquidMassFraction => binding.liquid_mass_fraction_pullback(cx, &seed),
        }
    }
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
    /// Target in kelvin, J/kg or dimensionless mass fraction as selected above.
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
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
) -> Result<ScanObjectiveEvaluation, ScanSensitivityError> {
    evaluate_with_tape(cx, run, endpoints, observations, sampling, step, transpose,
        budget, |_, _| Ok(())).map(|(objective, ())| objective)
}

// Keep one primal/least-squares owner. Additional response pullbacks borrow its
// checked endpoint tape before it is dropped; they cannot mutate the trajectory.
#[allow(clippy::too_many_arguments)]
fn evaluate_with_tape<R, O: ScanObservationSpec>(
    cx: &Cx<'_>, mut run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    responses: impl FnOnce(&ScanningEnthalpyRun<'_, '_, '_>, &[ScanStepLinearization<'_>])
        -> Result<R, ScanSensitivityError>,
) -> Result<(ScanObjectiveEvaluation, R), ScanSensitivityError> {
    let (at_endpoint, per_step) = admit_layout(cx, &run, endpoints, observations,
        budget, endpoints.len())?;
    let count = endpoints.len();
    let vertices = run.conduction().temperature().len();
    let mut tape: Vec<ScanStepLinearization<'_>> = Vec::new();
    tape.try_reserve_exact(count).map_err(|_| invalid("scan objective endpoint allocation refused"))?;
    let mut predictions = zeros(observations.len())?;
    let mut residuals = zeros(observations.len())?;
    let mut value = 0.0;
    for (index, &time) in endpoints.iter().enumerate() {
        poll(cx, index)?;
        let binding = advance_binding(cx, &mut run, time, sampling, step, per_step)?;
        for &i in &at_endpoint[index] {
            poll(cx, i)?;
            let observation = &observations[i];
            let prediction = observation.prediction(cx, binding.report())?;
            let residual = finite((prediction-observation.target())/observation.scale())?;
            value = finite(value + 0.5 * residual * residual)?;
            predictions[i] = prediction;
            residuals[i] = residual;
        }
        tape.push(binding);
    }
    let mut reverse = ReverseState::new(vertices)?;
    for index in (0..count).rev() {
        reverse.inject_observations(cx, &tape[index], observations, &at_endpoint[index], &residuals)?;
        reverse.pull(cx, &run, &tape[index], sampling, transpose, budget.max_schedule_controls)?;
    }
    let final_step = tape.last().ok_or_else(|| invalid("missing scan endpoint"))?.report().clone();
    poll(cx, count)?;
    let extra = responses(&run, &tape)?;
    poll(cx, count)?;
    Ok((ScanObjectiveEvaluation { value, predictions, standardized_residuals: residuals,
        sources: reverse.sources.ok_or_else(|| invalid("missing scan reverse sweep"))?,
        initial_specific_enthalpy: reverse.seed, final_step, endpoint_count: count,
        transpose_iterations: reverse.iterations, max_transpose_relative_residual: reverse.residual }, extra))
}

// Admission is shared by full-storage and checkpointed objectives. Only the
// number of simultaneously retained tangents changes; the physical grid does not.
fn admit_layout<O: ScanObservationSpec>(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    endpoints: &[f64], observations: &[O], budget: ScanObjectiveBudget,
    retained_steps: usize,
) -> Result<(Vec<Vec<usize>>, ScanSensitivityBudget), ScanSensitivityError> {
    poll(cx, 0)?;
    let count = endpoints.len();
    let vertices = run.conduction().temperature().len();
    let remaining = run.conduction().max_steps().checked_sub(run.conduction().accepted_steps());
    if count == 0 || count > budget.max_steps || remaining.is_none_or(|n| count > n)
        || observations.is_empty() || observations.len() > budget.max_observations
    { return Err(invalid("nonempty endpoints/observations must fit the objective and remaining physical-step budgets")); }
    if retained_steps == 0 || retained_steps > count {
        return Err(invalid("scan objective needs a nonempty retained tangent segment"));
    }
    let retained = retained_steps.checked_mul(vertices)
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
        if observation.endpoint() >= count
            || !observation.target().is_finite() || !observation.scale().is_finite() || observation.scale() <= 0.0
        { return Err(invalid("scan observation needs an existing endpoint/vertex, finite target and positive finite physical scale")); }
        observation.validate_readout(run.run.problem.boundary.mesh())?;
        if observation.observable() == ScanObservable::LiquidMassFraction
            && !(0.0..=1.0).contains(&observation.target())
        { return Err(invalid("liquid mass fraction targets must be dimensionless values in [0,1], not percentages")); }
        at_endpoint[observation.endpoint()].push(i);
    }
    let per_step = ScanSensitivityBudget {
        max_source_vertices: budget.max_retained_source_vertices / retained_steps,
        max_radiation_feedback_entries: budget.max_retained_radiation_entries / retained_steps,
    };
    Ok((at_endpoint, per_step))
}

fn advance_binding<'m>(cx: &Cx<'_>, run: &mut ScanningEnthalpyRun<'_, 'm, '_>,
    time: f64, sampling: ScanObjectiveSampling, step: EnthalpyStepConfig,
    per_step: ScanSensitivityBudget,
) -> Result<ScanStepLinearization<'m>, ScanSensitivityError> {
    let binding = match sampling {
        ScanObjectiveSampling::RightEndpoint => run.advance_linearized_to(cx, time, step, per_step)?,
        ScanObjectiveSampling::EventLeft => run.advance_linearized_event_step(
            cx, time, time-run.conduction().time_s(), step, per_step)?,
    };
    if binding.end_time_s() != time {
        return Err(invalid("fixed scan objective grid omits a source/boundary event or selected motion stop"));
    }
    Ok(binding)
}

fn observed(report: &ScanStep, observable: ScanObservable, vertex: usize) -> f64 {
    match observable {
        ScanObservable::Temperature => report.solution().temperature[vertex],
        ScanObservable::SpecificEnthalpy => report.solution().specific_enthalpy_j_kg[vertex],
        ScanObservable::LiquidMassFraction => report.solution().liquid_mass_fraction[vertex],
    }
}

// Only the original endpoint pullback does numerical differentiation. Reusing
// this accumulator keeps checkpoint boundaries invisible to summation order.
struct ReverseState {
    seed: Vec<f64>,
    sources: Option<Vec<MovingHeatScheduleGradient>>,
    iterations: usize,
    residual: f64,
}
impl ReverseState {
    fn new(vertices: usize) -> Result<Self, ScanSensitivityError> {
        Ok(Self { seed: zeros(vertices)?, sources: None, iterations: 0, residual: 0.0 })
    }

    fn inject_observations<O: ScanObservationSpec>(&mut self, cx: &Cx<'_>, binding: &ScanStepLinearization<'_>,
        observations: &[O], indices: &[usize], residuals: &[f64],
    ) -> Result<(), ScanSensitivityError> {
        let mut temperature_seed = zeros(self.seed.len())?;
        let mut fraction_seed = None;
        for &i in indices {
            poll(cx, i)?;
            let observation = &observations[i];
            let local = finite(residuals[i]/observation.scale())?;
            let target = match observation.observable() {
                ScanObservable::Temperature => &mut temperature_seed,
                ScanObservable::SpecificEnthalpy => &mut self.seed,
                ScanObservable::LiquidMassFraction => {
                    if fraction_seed.is_none() { fraction_seed = Some(zeros(self.seed.len())?); }
                    fraction_seed.as_mut().expect("fraction seed allocated")
                },
            };
            observation.add_seed(cx, target, local)?;
        }
        let local_h = binding.transport().transport().temperature_pullback(cx, &temperature_seed)?;
        add(cx, &mut self.seed, &local_h)?;
        if let Some(seed) = fraction_seed {
            let local_h = binding.liquid_mass_fraction_pullback(cx, &seed)?;
            add(cx, &mut self.seed, &local_h)?;
        }
        Ok(())
    }

    fn pull(&mut self, cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
        binding: &ScanStepLinearization<'_>, sampling: ScanObjectiveSampling,
        transpose: LinearConfig, max_controls: usize,
    ) -> Result<(), ScanSensitivityError> {
        poll(cx, 0)?;
        let gradient = binding.pullback(cx, &self.seed, transpose)?;
        self.iterations = self.iterations.checked_add(gradient.storage.transport.iterations)
            .ok_or_else(|| invalid("scan objective transpose work count overflow"))?;
        self.residual = self.residual.max(gradient.storage.transport.relative_residual);
        if gradient.sources.len() != run.sources().len() {
            return Err(invalid("scan source-gradient layout changed"));
        }
        let mut schedules = Vec::with_capacity(gradient.sources.len());
        for (i, (source, gradient)) in run.sources().iter().zip(gradient.sources).enumerate() {
            poll(cx, i)?;
            schedules.push(match sampling {
                ScanObjectiveSampling::RightEndpoint => source.pullback_schedule_at(
                    cx, binding.end_time_s(), gradient, max_controls)?,
                ScanObjectiveSampling::EventLeft => source.pullback_schedule_before(
                    cx, binding.end_time_s(), gradient, max_controls)?,
            });
        }
        match &mut self.sources {
            None => self.sources = Some(schedules),
            Some(total) => {
                if total.len() != schedules.len() { return Err(invalid("scan source identity changed during reverse sweep")); }
                for (a, b) in total.iter_mut().zip(&schedules) {
                    for axis in 0..3 { add(cx, &mut a.center_values_m[axis], &b.center_values_m[axis])?; }
                    add(cx, &mut a.half_width_m, &b.half_width_m)?;
                    add(cx, &mut a.absorbed_power_values_w, &b.absorbed_power_values_w)?;
                }
            }
        }
        self.seed = gradient.storage.transport.previous_specific_enthalpy;
        Ok(())
    }
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

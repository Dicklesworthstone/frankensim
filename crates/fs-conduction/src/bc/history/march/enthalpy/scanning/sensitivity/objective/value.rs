//! Value-only scan objectives for physical screening and derivative-free search.
//!
//! Every endpoint uses the original primal producer, including nonlinear
//! conductivity, latent heat, radiation, source normalization and energy gates.
//! No tangent is built, so a valid primal at a constitutive corner need not be
//! rejected merely because a classical derivative does not exist there.
//! No gradient, stationary-point, constraint-feasibility or optimality claim is
//! attached to the returned value. Only the fitting objective is evaluated.

/// Sampled state and hot-spot responses on the same value-only trajectory.
pub mod responses;

use super::{ScanObservable, ScanObjectiveSampling, ScanSensitivityError,
    ScanStep, ScanningEnthalpyRun, EnthalpyStepConfig, finite, invalid, poll, zeros};
use super::spatial::ScanObservationSpec;
use fs_exec::Cx;

/// Whole-trajectory admission without irrelevant retained-tangent allowances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanValueBudget {
    /// Maximum original physical endpoints, also limited by the run's allowance.
    pub max_steps: usize,
    /// Maximum physical readouts and retained scalar predictions/residuals.
    pub max_observations: usize,
    /// Combined authored source coordinates and half widths.
    pub max_schedule_controls: usize,
}

/// Original physical prediction and fitting loss, deliberately without gradients.
#[derive(Debug)]
pub struct ScanValueEvaluation {
    /// Sum of the declared residual losses, in original observation order per time.
    pub value: f64,
    /// Original input observation order, in each observation's physical units.
    pub predictions: Vec<f64>,
    /// Raw standardized residuals, including those downweighted by robust loss.
    pub standardized_residuals: Vec<f64>,
    /// Unmodified final physical report, including phase fields and energy checks.
    pub final_step: ScanStep,
    /// Original primal solves. There is no replay or transpose solve on this path.
    pub endpoint_count: usize,
}

/// Evaluate vertex or spatial measurements without requiring differentiability.
///
/// The trial run is consumed. Only its current state and one endpoint report
/// are retained; observations and the endpoint index are explicitly bounded.
/// Sampling uses the same fixed right-endpoint or event-left convention as the
/// adjoint evaluator, and omitted source/boundary/motion stops remain errors.
/// Loss is applied after spatial readout; no separate targets are invented for
/// the sensor's supporting vertices. Physical failures and cancellation propagate.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_values<O: ScanObservationSpec>(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling, step: EnthalpyStepConfig,
    budget: ScanValueBudget,
) -> Result<ScanValueEvaluation, ScanSensitivityError> {
    evaluate_with_observer(cx, run, endpoints, observations, sampling, step, budget,
        |_, _| Ok(()))
}

// Share the original primal and scoring loop. Response observers may inspect
// accepted reports, never alter the trajectory or replace the physical solve.
#[allow(clippy::too_many_arguments)]
fn evaluate_with_observer<O: ScanObservationSpec>(
    cx: &Cx<'_>, mut run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling, step: EnthalpyStepConfig,
    budget: ScanValueBudget,
    mut observe: impl FnMut(usize, &ScanStep) -> Result<(), ScanSensitivityError>,
) -> Result<ScanValueEvaluation, ScanSensitivityError> {
    let at_endpoint = admit_values(cx, &run, endpoints, observations, budget)?;
    let mut predictions = zeros(observations.len())?;
    let mut residuals = zeros(observations.len())?;
    let mut value = 0.0;
    let mut final_step = None;
    for (index, &time) in endpoints.iter().enumerate() {
        poll(cx, index)?;
        let report = match sampling {
            ScanObjectiveSampling::RightEndpoint => run.advance_to(cx, time, step)?,
            ScanObjectiveSampling::EventLeft => {
                let sampled = run.advance_event_step(cx, time, time - run.conduction().time_s(), step)?;
                if sampled.end_time_s != time {
                    return Err(invalid("fixed scan objective grid omits a source/boundary event or selected motion stop"));
                }
                sampled.report
            }
        };
        for &i in &at_endpoint[index] {
            poll(cx, i)?;
            let observation = &observations[i];
            let prediction = observation.prediction(cx, &report)?;
            let residual = finite((prediction - observation.target()) / observation.scale())?;
            value = finite(value + observation.loss().value_derivative(residual)?.0)?;
            predictions[i] = prediction;
            residuals[i] = residual;
        }
        observe(index, &report)?;
        if index + 1 == endpoints.len() { final_step = Some(report); }
    }
    poll(cx, endpoints.len())?;
    Ok(ScanValueEvaluation { value, predictions, standardized_residuals: residuals,
        final_step: final_step.ok_or_else(|| invalid("missing scan value endpoint"))?,
        endpoint_count: endpoints.len() })
}

// Shared with full/replayed adjoints. Only tangent storage admission differs.
pub(super) fn admit_values<O: ScanObservationSpec>(cx: &Cx<'_>,
    run: &ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64], observations: &[O],
    budget: ScanValueBudget,
) -> Result<Vec<Vec<usize>>, ScanSensitivityError> {
    poll(cx, 0)?;
    let count = endpoints.len();
    let remaining = run.conduction().max_steps().checked_sub(run.conduction().accepted_steps());
    if count == 0 || count > budget.max_steps || remaining.is_none_or(|n| count > n)
        || observations.is_empty() || observations.len() > budget.max_observations
    { return Err(invalid("nonempty endpoints/observations must fit the objective and remaining physical-step budgets")); }
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
        if !time.is_finite() || !(time - previous).is_finite() || time <= previous {
            return Err(invalid("scan objective endpoints must increase from the input state with finite positive intervals"));
        }
        previous = time;
    }
    let mut at_endpoint = Vec::new();
    at_endpoint.try_reserve_exact(count).map_err(|_| invalid("scan observation index allocation refused"))?;
    at_endpoint.resize_with(count, Vec::new);
    for (i, observation) in observations.iter().enumerate() {
        poll(cx, i)?;
        if observation.endpoint() >= count || !observation.target().is_finite()
            || !observation.scale().is_finite() || observation.scale() <= 0.0
        { return Err(invalid("scan observation needs an existing endpoint/vertex, finite target and positive finite physical scale")); }
        observation.validate_readout(run.run.problem.boundary.mesh())?;
        observation.loss().validate()?;
        if observation.observable() == ScanObservable::LiquidMassFraction
            && !(0.0..=1.0).contains(&observation.target())
        { return Err(invalid("liquid mass fraction targets must be dimensionless values in [0,1], not percentages")); }
        at_endpoint[observation.endpoint()].push(i);
    }
    poll(cx, count)?;
    Ok(at_endpoint)
}

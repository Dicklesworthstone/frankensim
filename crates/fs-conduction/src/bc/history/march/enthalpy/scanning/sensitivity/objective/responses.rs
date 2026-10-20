//! Independent state responses for constrained inverse design.
//!
//! A temperature or enthalpy limit needs the derivative of that actual response,
//! not the derivative of a squared fitting residual (which vanishes at its target).
//! Reuse the objective's forward tape, then perform one reverse sweep per response.
//! Earlier heat and the selected side of each schedule event remain represented.
//! Only supplied vertices and endpoints are observed: no continuous-time or
//! between-node maximum, extrema smoothing, or global safety bound is implied.

use super::{ScanObservable, ScanObjectiveBudget, ScanObjectiveEvaluation, ScanObservation,
    ScanObjectiveSampling, ScanSensitivityError, ScanStepLinearization, ScanningEnthalpyRun,
    MovingHeatScheduleGradient, EnthalpyStepConfig, LinearConfig, add, evaluate_with_tape,
    finite, invalid, poll, zeros};
use fs_exec::Cx;

/// A response in its original physical units at one supplied endpoint and vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanResponse {
    pub endpoint: usize,
    pub vertex: usize,
    pub observable: ScanObservable,
}

/// Separate allowance for extra reverse sweeps and returned response gradients.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanResponseBudget {
    pub max_responses: usize,
    /// Responses times (all authored source controls plus mesh vertices).
    /// A scalar shape cap, not a global allocator or byte bound.
    pub max_gradient_entries: usize,
    /// Objective endpoints plus sum(response.endpoint + 1), admitted before
    /// running physics. This bounds transpose calls, not Krylov iterations.
    pub max_transpose_solves: usize,
}

/// Physical response and derivatives, not a penalty or a squared residual.
#[derive(Debug)]
pub struct ScanResponseEvaluation {
    pub response: ScanResponse,
    /// Kelvin or J/kg, as selected by response.observable.
    pub value: f64,
    pub sources: Vec<MovingHeatScheduleGradient>,
    pub initial_specific_enthalpy: Vec<f64>,
    pub transpose_iterations: usize,
    pub max_transpose_relative_residual: f64,
}

/// Evaluate the original objective and independently seeded state responses.
/// There is ONE original forward simulation per candidate, not one per limit.
/// Every reverse solve retains the existing material/phase smoothness, residual,
/// energy, radiation and cancellation rules. Responses may precede the final
/// observation; subsequent source values then receive exactly zero sensitivity.
///
/// A failed response or objective returns neither partial result. The input run
/// is consumed, just as in evaluate_scan_objective. Times, event sides, mesh and
/// constitutive charts are fixed design inputs, not differentiated decisions.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_objective_with_responses(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    responses: &[ScanResponse], response_budget: ScanResponseBudget,
) -> Result<(ScanObjectiveEvaluation, Vec<ScanResponseEvaluation>), ScanSensitivityError> {
    poll(cx, 0)?;
    if responses.len() > response_budget.max_responses {
        return Err(invalid("scan response count exceeds its explicit allowance"));
    }
    let vertices = run.conduction().temperature().len();
    let mut entries_per_response = vertices;
    for (i, source) in run.sources().iter().enumerate() {
        poll(cx, i)?;
        entries_per_response = entries_per_response.checked_add(source.schedule_control_count()?)
            .ok_or_else(|| invalid("scan response gradient size overflow"))?;
    }
    if entries_per_response.checked_mul(responses.len())
        .is_none_or(|n| n > response_budget.max_gradient_entries)
    { return Err(invalid("scan response gradients exceed their retained scalar allowance")); }
    let mut solves = endpoints.len();
    for (i, response) in responses.iter().enumerate() {
        poll(cx, i)?;
        if response.endpoint >= endpoints.len() || response.vertex >= vertices {
            return Err(invalid("scan response must name an existing endpoint and mesh vertex"));
        }
        solves = solves.checked_add(response.endpoint + 1)
            .ok_or_else(|| invalid("scan response transpose count overflow"))?;
    }
    if solves > response_budget.max_transpose_solves {
        return Err(invalid("scan objective and response sweeps exceed their transpose allowance"));
    }
    evaluate_with_tape(cx, run, endpoints, observations, sampling, step, transpose, budget,
        |run, tape| {
            let mut output = Vec::new();
            output.try_reserve_exact(responses.len())
                .map_err(|_| invalid("scan response output allocation refused"))?;
            for &response in responses {
                poll(cx, output.len())?;
                output.push(pull_response(cx, run, tape, response, sampling, transpose,
                    budget.max_schedule_controls)?);
            }
            Ok(output)
        })
}

fn pull_response(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    tape: &[ScanStepLinearization<'_>], response: ScanResponse,
    sampling: ScanObjectiveSampling, transpose: LinearConfig, max_controls: usize,
) -> Result<ScanResponseEvaluation, ScanSensitivityError> {
    let solution = tape[response.endpoint].report().solution();
    let value = finite(match response.observable {
        ScanObservable::Temperature => solution.temperature[response.vertex],
        ScanObservable::SpecificEnthalpy => solution.specific_enthalpy_j_kg[response.vertex],
    })?;
    let mut seed = zeros(solution.temperature.len())?;
    seed[response.vertex] = 1.0;
    let mut accumulated: Option<Vec<MovingHeatScheduleGradient>> = None;
    let mut iterations = 0usize;
    let mut residual = 0.0_f64;
    for index in (0..=response.endpoint).rev() {
        poll(cx, index)?;
        let binding = &tape[index];
        // Only the local temperature seed is transformed by dT/dh. All earlier
        // seeds already represent derivatives with respect to stored enthalpy.
        let gradient = if index == response.endpoint && response.observable == ScanObservable::Temperature {
            binding.pullback_temperature(cx, &seed, transpose)?
        } else { binding.pullback(cx, &seed, transpose)? };
        iterations = iterations.checked_add(gradient.storage.transport.iterations)
            .ok_or_else(|| invalid("scan response transpose iteration count overflow"))?;
        residual = residual.max(gradient.storage.transport.relative_residual);
        if gradient.sources.len() != run.sources().len() {
            return Err(invalid("scan response source-gradient layout changed"));
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
        match &mut accumulated {
            None => accumulated = Some(schedules),
            Some(total) => {
                if total.len() != schedules.len() { return Err(invalid("scan response source identity changed")); }
                for (a, b) in total.iter_mut().zip(&schedules) {
                    for axis in 0..3 { add(cx, &mut a.center_values_m[axis], &b.center_values_m[axis])?; }
                    add(cx, &mut a.half_width_m, &b.half_width_m)?;
                    add(cx, &mut a.absorbed_power_values_w, &b.absorbed_power_values_w)?;
                }
            }
        }
        seed = gradient.storage.transport.previous_specific_enthalpy;
    }
    poll(cx, response.endpoint)?;
    Ok(ScanResponseEvaluation { response, value,
        sources: accumulated.ok_or_else(|| invalid("missing scan response sweep"))?,
        initial_specific_enthalpy: seed, transpose_iterations: iterations,
        max_transpose_relative_residual: residual })
}

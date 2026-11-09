//! Independent state responses for constrained inverse design.
//!
//! A temperature or enthalpy limit needs the derivative of that actual response,
//! not the derivative of a squared fitting residual (which vanishes at its target).
//! Reuse the objective's forward tape, then perform one reverse sweep per response.
//! Earlier heat and the selected side of each schedule event remain represented.
//! Only supplied vertices and endpoints are observed: no continuous-time or
//! between-node maximum, extrema smoothing, or global safety bound is implied.

/// One-sweep smooth hot-spot envelopes over explicit vertex/time selections.
pub mod peaks;

use super::{ScanObservable, ScanObjectiveBudget, ScanObjectiveEvaluation, ScanObservation,
    ScanObjectiveSampling, ScanSensitivityError, ScanStepLinearization, ScanningEnthalpyRun,
    MovingHeatScheduleGradient, ScanRadiationScheduleGradient, EnthalpyStepConfig,
    LinearConfig, ReverseState, evaluate_with_tape, finite, invalid, poll};
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
    /// Responses times (all authored source/radiation controls plus mesh vertices).
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
    /// Kelvin, J/kg or dimensionless mass fraction, as selected by the observable.
    pub value: f64,
    pub sources: Vec<MovingHeatScheduleGradient>,
    /// Original patch order; empty when no radiation is attached.
    pub radiation: Vec<ScanRadiationScheduleGradient>,
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
    let entries_per_response = vertices.checked_add(super::radiation::control_count(cx, &run)?)
        .ok_or_else(|| invalid("scan response gradient size overflow"))?;
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
        ScanObservable::LiquidMassFraction => solution.liquid_mass_fraction[response.vertex],
    })?;
    let mut reverse = ReverseState::new(solution.temperature.len())?;
    reverse.seed[response.vertex] = 1.0;
    // Convert only the local response seed, never the carried h derivative.
    reverse.seed = response.observable.seed_in_enthalpy(cx, &tape[response.endpoint], reverse.seed)?;
    for index in (0..=response.endpoint).rev() {
        poll(cx, index)?;
        reverse.pull(cx, run, &tape[index], sampling, transpose, max_controls)?;
    }
    poll(cx, response.endpoint)?;
    Ok(ScanResponseEvaluation { response, value,
        sources: reverse.sources.ok_or_else(|| invalid("missing scan response sweep"))?,
        radiation: reverse.radiation, initial_specific_enthalpy: reverse.seed,
        transpose_iterations: reverse.iterations, max_transpose_relative_residual: reverse.residual })
}

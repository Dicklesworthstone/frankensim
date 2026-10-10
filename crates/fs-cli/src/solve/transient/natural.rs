//! Natural-convection endpoints with immutable transient solid history.
//!
//! The existing card owner supplies h at the candidate area-mean wall
//! temperature. The shared backward-Euler evaluator independently reassembles
//! the complete physical endpoint, including k(T), contact and optional
//! radiation, before an h fixed point can publish a new temperature.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use fs_conduction::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, StepConfig, StepSolution,
};
use fs_conduction::{ConductionProblem, RobinFlux, ThermalBoundary, ThermalInterfaces};
use fs_exec::Cx;

use super::super::{conduction_error, natural as boundary_law, radiation, SolveRefusal};
use super::{bad, lower, nonlinear_row, poll, radiation_row, RadiativeEndpoint};

pub(super) struct NaturalStep {
    pub(super) endpoint: StepSolution,
    /// Final counted inner response; cumulative work is separate below.
    pub(super) nonlinear_evidence: String,
    pub(super) radiation_evidence: String,
    pub(super) radiation: Option<RadiativeEndpoint>,
    /// Natural coefficients actually used by that inner response.
    pub(super) applied_boundary: ThermalBoundary,
    pub(super) applied_coefficients: BTreeMap<String, f64>,
    /// Coefficients and heat evaluated at the accepted physical endpoint.
    pub(super) physical_coefficients: BTreeMap<String, f64>,
    pub(super) physical_convective_robin_fluxes: Vec<RobinFlux>,
    pub(super) receipt: String,
    pub(super) iterations: usize,
    pub(super) solid_solves: usize,
    pub(super) krylov_iterations: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
    pub(super) radiation_trials: usize,
    /// Direct physical and frozen-predictor audits by this outer producer.
    pub(super) endpoint_evaluations: usize,
    /// Inner radiative physical gates only, distinct from the natural gate.
    pub(super) maximum_physical_energy_residual_j: f64,
    pub(super) physical_residual_norm_j: f64,
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy_residual_j: f64,
    pub(super) physical_dirichlet_in_w: f64,
    pub(super) physical_convective_out_w: f64,
    pub(super) physical_radiative_out_w: f64,
}

fn finite(value: f64, what: &str) -> Result<f64, SolveRefusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!("transient natural-convection {what} is not finite")))
    }
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), SolveRefusal> {
    *total = total.checked_add(added)
        .filter(|value| *value <= cap)
        .ok_or_else(|| bad("transient natural convection exceeded its numerical work allowance"))?;
    Ok(())
}

fn tighter_relative(requested: f64, allowance_j: f64, initial_norm_j: f64)
    -> Result<f64, SolveRefusal>
{
    let tolerance = if initial_norm_j == 0.0 {
        requested
    } else {
        requested.min(allowance_j / initial_norm_j)
    };
    if tolerance.is_finite() && tolerance > 0.0 {
        Ok(tolerance)
    } else {
        Err(bad("transient natural convection cannot represent its inner residual tolerance"))
    }
}

/// Solve one bounded natural fixed point. The initial h guess comes from the
/// existing card owner, so a powered ambient start need not query an undefined
/// equilibrium law. Every candidate is checked against the ACTUAL card; exact
/// equilibrium and out-of-domain Rayleigh values keep their original refusals.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance(
    cx: &Cx<'_>,
    engine: &BackwardEuler<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    old: &[f64],
    dt_s: f64,
    step_config: StepConfig,
    nonlinear: Option<NonlinearStepConfig>,
    radiation: Option<&radiation::LoweredRadiation>,
    laws: &[boundary_law::NaturalLaw],
    pressure_pa: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<NaturalStep, SolveRefusal> {
    poll(cx, deadline)?;
    if laws.is_empty()
        || old.len() != problem.mesh.vertex_count()
        || !(dt_s.is_finite() && dt_s > 0.0
            && step_config.energy_tolerance_j.is_finite() && step_config.energy_tolerance_j > 0.0)
    {
        return Err(bad("transient natural convection requires laws, matching history and positive time/energy controls"));
    }
    let mut names = BTreeSet::new();
    let mut regions = Vec::with_capacity(laws.len());
    for law in laws {
        if !names.insert(law.target.as_str()) {
            return Err(bad("transient natural convection has duplicate targets"));
        }
        let region = problem.boundary.region_names().iter()
            .position(|name| name == &law.target)
            .ok_or_else(|| bad(format!(
                "natural-convection target '{}' is absent from the solid boundary", law.target,
            )))?;
        regions.push(region);
    }
    let mut coefficients = boundary_law::initial_coefficients(laws, pressure_pa)?;
    let mut predictor = old.to_vec();
    for &(vertex, prescribed) in problem.boundary.dirichlet() {
        *predictor.get_mut(vertex)
            .ok_or_else(|| bad("natural-convection prescribed vertex lies outside the history"))? = prescribed;
    }
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    if let Some(lowering) = radiation {
        lowering.config.validate().map_err(lower)?;
    }
    if let Some(policy) = nonlinear {
        policy.validate().map_err(lower)?;
    }
    let max_solid_solves = boundary_law::MAX_ITERATIONS.checked_mul(radiation_factor)
        .ok_or_else(|| bad("natural-convection/radiation work allowance overflows"))?;
    let max_krylov = max_solid_solves.checked_mul(step_config.linear.max_iterations)
        .ok_or_else(|| bad("natural-convection Krylov allowance overflows"))?;
    let max_updates = max_solid_solves.checked_mul(nonlinear.map_or(0, |p| p.max_iterations))
        .ok_or_else(|| bad("natural-convection Newton allowance overflows"))?;
    let max_backtracks = max_updates.checked_mul(nonlinear.map_or(0, |p| p.line_search.max_backtracks))
        .ok_or_else(|| bad("natural-convection backtrack allowance overflows"))?;
    let max_endpoint_evaluations = boundary_law::MAX_ITERATIONS.checked_mul(2)
        .ok_or_else(|| bad("natural-convection endpoint audit allowance overflows"))?;
    if let Some(policy) = nonlinear {
        policy.line_search.max_backtracks.checked_add(1)
            .and_then(|trials| trials.checked_mul(max_updates))
            .ok_or_else(|| bad("natural-convection Newton trial allowance overflows"))?;
    }
    let patches = radiation.map(|lowering| lowering.patches.as_slice());
    let mut solid_solves = 0;
    let mut krylov_iterations = 0;
    let mut nonlinear_updates = 0;
    let mut nonlinear_backtracks = 0;
    let mut radiation_trials = 0;
    let mut endpoint_evaluations = 0;
    let mut maximum_physical_energy_residual_j = 0.0_f64;
    let mut physical_target: Option<f64> = None;
    for iteration in 1..=boundary_law::MAX_ITERATIONS {
        poll(cx, deadline)?;
        let replacements = laws.iter().zip(&regions).map(|(law, &region)| {
            let h = *coefficients.get(&law.target)
                .ok_or_else(|| bad("natural-convection trial lost an applied coefficient"))?;
            Ok((region, h, law.ambient_k))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let boundary = problem.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
        let trial_problem = ConductionProblem { boundary: &boundary, ..problem };
        let mut inner_step = step_config;
        let mut inner_nonlinear = nonlinear;
        if let Some(target) = physical_target {
            // This is a residual of the supplied FROZEN h, not a query of the
            // natural law at old=ambient. The physical target never changes.
            count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
            let initial = engine.evaluate_prescribed_endpoint(
                cx, trial_problem, interfaces, old, &predictor, dt_s, step_config, patches,
            ).map_err(lower)?;
            if let Some(policy) = inner_nonlinear.as_mut() {
                policy.residual_atol_j = policy.residual_atol_j.min(target * 0.125);
                policy.residual_rtol = tighter_relative(
                    policy.residual_rtol, target * 0.125, initial.residual_norm_j,
                )?;
            } else {
                inner_step.linear.tolerance = tighter_relative(
                    step_config.linear.tolerance, target * 0.125, initial.residual_norm_j,
                )?;
                if radiation.is_some() {
                    // The radiation owner's linear physical gate includes a
                    // joule-derived rounding floor. Reserve space under the
                    // natural target without widening the final energy gate.
                    inner_step.energy_tolerance_j = step_config.energy_tolerance_j * 0.125;
                    if !(inner_step.energy_tolerance_j.is_finite() && inner_step.energy_tolerance_j > 0.0) {
                        return Err(bad("natural-convection inner energy allowance is not representable"));
                    }
                }
            }
        }
        let (endpoint, nonlinear_evidence, radiation_evidence, radiative) =
            if let Some(lowering) = radiation {
                let solved = engine.advance_prescribed_with_ambient_radiation(
                    cx, trial_problem, interfaces, old, dt_s, inner_step, inner_nonlinear,
                    &lowering.patches, lowering.config,
                ).map_err(lower)?;
                count(&mut solid_solves, solved.radiation.iterations, max_solid_solves)?;
                count(&mut radiation_trials, solved.radiation.iterations, max_solid_solves)?;
                count(&mut krylov_iterations, solved.radiation.krylov_iterations, max_krylov)?;
                count(&mut nonlinear_updates, solved.radiation.solid_iterations, max_updates)?;
                count(&mut nonlinear_backtracks, solved.nonlinear_backtracks, max_backtracks)?;
                maximum_physical_energy_residual_j = maximum_physical_energy_residual_j
                    .max(solved.physical_energy_residual_j.abs());
                let nonlinear_evidence = solved.nonlinear.as_ref().map(nonlinear_row)
                    .transpose()?.unwrap_or_else(|| "null".to_string());
                let radiation_evidence = radiation_row(&solved)?;
                let radiative = RadiativeEndpoint {
                    combined_boundary: solved.combined_boundary,
                    convective_robin_fluxes: solved.convective_robin_fluxes,
                    convective_out_w: solved.convective_out_w,
                    report: solved.radiation,
                };
                (solved.conduction, nonlinear_evidence, radiation_evidence, Some(radiative))
            } else if let Some(policy) = inner_nonlinear {
                let solved = engine.advance_nonlinear_prescribed(
                    cx, trial_problem, interfaces, old, dt_s, inner_step, policy,
                ).map_err(lower)?;
                count(&mut solid_solves, 1, max_solid_solves)?;
                count(&mut krylov_iterations, solved.step.krylov_iterations, max_krylov)?;
                count(&mut nonlinear_updates, solved.nonlinear_iterations, max_updates)?;
                count(&mut nonlinear_backtracks, solved.backtracks, max_backtracks)?;
                let row = nonlinear_row(&solved)?;
                (solved.step, row, "null".to_string(), None)
            } else {
                let endpoint = engine.advance_prescribed(
                    cx, trial_problem, interfaces, old, dt_s, inner_step,
                ).map_err(lower)?;
                count(&mut solid_solves, 1, max_solid_solves)?;
                count(&mut krylov_iterations, endpoint.krylov_iterations, max_krylov)?;
                (endpoint, "null".to_string(), "null".to_string(), None)
            };
        poll(cx, deadline)?;
        let fluxes = radiative.as_ref().map_or(
            endpoint.robin_fluxes.as_slice(), |result| result.convective_robin_fluxes.as_slice(),
        );
        let mut next = BTreeMap::new();
        let mut points = Vec::with_capacity(laws.len());
        let mut worst_change = 0.0_f64;
        let mut physical_rows = Vec::with_capacity(laws.len());
        for (law, &region) in laws.iter().zip(&regions) {
            let flux = fluxes.iter().find(|flux| flux.region == law.target)
                .ok_or_else(|| bad("natural-convection target has no solved convective trace"))?;
            let point = boundary_law::coefficient(
                law, finite(flux.mean_wall_temperature_k - law.ambient_k, "wall departure")?, pressure_pa,
            )?;
            let applied = *coefficients.get(&law.target)
                .ok_or_else(|| bad("natural-convection applied coefficient is missing"))?;
            worst_change = worst_change.max(finite(
                (point.htc_w_m2_k - applied).abs() / applied, "relative coefficient change",
            )?);
            next.insert(law.target.clone(), point.htc_w_m2_k);
            physical_rows.push((region, point.htc_w_m2_k, law.ambient_k));
            points.push(point);
        }
        let physical_boundary = problem.boundary.with_uniform_robin_replacements(&physical_rows)
            .map_err(lower)?;
        count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
        let physical = engine.evaluate_prescribed_endpoint(
            cx, ConductionProblem { boundary: &physical_boundary, ..problem },
            interfaces, old, &endpoint.temperature, dt_s, step_config, patches,
        ).map_err(lower)?;
        let target = match physical_target {
            Some(target) => target,
            None => {
                // A physical initial equilibrium coefficient is undefined by
                // this card. Freeze the scale at the first fully evaluated
                // candidate instead; neither later iteration nor time can
                // inflate it. The absolute floor is a part of the J budget.
                let absolute = 0.01 * step_config.energy_tolerance_j
                    / fs_math::det::sqrt(physical.free_dofs as f64);
                let relative = nonlinear.map_or(step_config.linear.tolerance, |p| p.residual_rtol);
                let target = finite(absolute + relative * physical.residual_norm_j, "physical residual target")?;
                if target <= 0.0 {
                    return Err(bad("natural-convection physical residual target is not representable"));
                }
                physical_target = Some(target);
                target
            }
        };
        poll(cx, deadline)?;
        if worst_change <= boundary_law::TOLERANCE_REL
            && physical.residual_norm_j <= target
            && physical.energy_residual_j.abs() <= step_config.energy_tolerance_j
        {
            let converged = laws.iter().zip(points).map(|(law, coefficient)| {
                let flux = physical.convective_robin_fluxes.iter()
                    .find(|flux| flux.region == law.target)
                    .ok_or_else(|| bad("physical natural-convection trace is missing"))?;
                Ok(boundary_law::Converged {
                    law: law.clone(), coefficient,
                    mean_wall_k: flux.mean_wall_temperature_k,
                    heat_rate_w: flux.heat_rate_w,
                })
            }).collect::<Result<Vec<_>, SolveRefusal>>()?;
            let receipt = boundary_law::receipt_fragment(&converged, iteration)?;
            return Ok(NaturalStep {
                endpoint, nonlinear_evidence, radiation_evidence, radiation: radiative,
                applied_boundary: boundary, applied_coefficients: coefficients,
                physical_coefficients: next,
                physical_convective_robin_fluxes: physical.convective_robin_fluxes,
                receipt, iterations: iteration,
                solid_solves, krylov_iterations, nonlinear_updates, nonlinear_backtracks,
                radiation_trials, endpoint_evaluations, maximum_physical_energy_residual_j,
                physical_residual_norm_j: physical.residual_norm_j,
                physical_residual_tolerance_j: target,
                physical_energy_residual_j: physical.energy_residual_j,
                physical_dirichlet_in_w: physical.dirichlet_in_w,
                physical_convective_out_w: physical.convective_out_w,
                physical_radiative_out_w: physical.radiation_out_w,
            });
        }
        if iteration == boundary_law::MAX_ITERATIONS {
            return Err(conduction_error(
                "cli-solve-conduction-transient-natural-unconverged",
                format!(
                    "natural convection exhausted {iteration} iterations: relative h change {worst_change:.4e}, physical residual {:.4e} J against {target:.4e} J, energy residual {:.4e} J against {:.4e} J",
                    physical.residual_norm_j, physical.energy_residual_j, step_config.energy_tolerance_j,
                ),
                "inspect the card regime and numerical tolerances; no partial thermal endpoint is published",
            ));
        }
        coefficients = next;
    }
    Err(bad("natural convection has no admitted iteration budget"))
}

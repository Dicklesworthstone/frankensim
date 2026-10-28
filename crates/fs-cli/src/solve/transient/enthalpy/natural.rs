//! Natural-convection fixed points over one immutable enthalpy time step.
//!
//! The existing card owns h(T_wall). The shared enthalpy endpoint evaluator
//! independently checks the returned candidate under that physical boundary,
//! including optional radiation, before any old enthalpy can be replaced.

use std::collections::{BTreeMap, BTreeSet};

use super::super::super::natural as boundary_law;
use super::*;

pub(super) struct NaturalStep {
    pub(super) endpoint: EnthalpyStepSolution,
    pub(super) radiation: Option<RadiativeEndpoint>,
    pub(super) enclosure: Option<radiation::EnclosureEndpoint>,
    pub(super) radiation_evidence: String,
    pub(super) applied_boundary: fs_conduction::ThermalBoundary,
    pub(super) applied_coefficients: BTreeMap<String, f64>,
    pub(super) physical_coefficients: BTreeMap<String, f64>,
    pub(super) physical_convective_robin_fluxes: Vec<fs_conduction::RobinFlux>,
    pub(super) receipt: String,
    pub(super) iterations: usize,
    pub(super) solid_solves: usize,
    pub(super) krylov_iterations: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) radiation_trials: usize,
    pub(super) endpoint_evaluations: usize,
    pub(super) maximum_radiation_energy_residual_j: f64,
    pub(super) physical_residual_norm_j: f64,
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy_residual_j: f64,
    pub(super) physical_convective_out_w: f64,
    pub(super) physical_radiative_out_w: f64,
    pub(super) physical_neumann_out_w: f64,
    pub(super) physical_dirichlet_in_w: f64,
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), SolveRefusal> {
    *total = total
        .checked_add(added)
        .filter(|value| *value <= cap)
        .ok_or_else(|| bad("enthalpy natural convection exceeded its numerical work allowance"))?;
    Ok(())
}

fn tighter_relative(
    requested: f64,
    allowance_j: f64,
    initial_norm_j: f64,
) -> Result<f64, SolveRefusal> {
    let tolerance = if initial_norm_j == 0.0 {
        requested
    } else {
        requested.min(allowance_j / initial_norm_j)
    };
    if tolerance.is_finite() && tolerance > 0.0 {
        Ok(tolerance)
    } else {
        Err(bad(
            "enthalpy natural convection cannot represent its inner residual tolerance",
        ))
    }
}

/// All counted inner solves use `old` unchanged. There is no apparent heat
/// capacity or latent account outside the shared enthalpy owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    old: &[f64],
    dt_s: f64,
    config: EnthalpyStepConfig,
    radiation: Option<&radiation::LoweredRadiation>,
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    laws: &[boundary_law::NaturalLaw],
    pressure_pa: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<NaturalStep, SolveRefusal> {
    poll(cx, deadline)?;
    if laws.is_empty()
        || old.len() != problem.mesh.vertex_count()
        || !(dt_s.is_finite() && dt_s > 0.0)
    {
        return Err(bad(
            "enthalpy natural convection requires laws, matching history and positive time",
        ));
    }
    let mut names = BTreeSet::new();
    let mut regions = Vec::with_capacity(laws.len());
    for law in laws {
        if !names.insert(law.target.as_str()) {
            return Err(bad("enthalpy natural convection has duplicate targets"));
        }
        regions.push(
            problem
                .boundary
                .region_names()
                .iter()
                .position(|name| name == &law.target)
                .ok_or_else(|| {
                    bad(format!(
                        "natural-convection target '{}' is absent from the enthalpy boundary",
                        law.target,
                    ))
                })?,
        );
    }
    let mut coefficients = boundary_law::initial_coefficients(laws, pressure_pa)?;
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    let max_solid_solves = boundary_law::MAX_ITERATIONS
        .checked_mul(radiation_factor)
        .ok_or_else(|| bad("natural-convection/radiation work allowance overflows"))?;
    let max_updates = max_solid_solves
        .checked_mul(config.max_newton_iterations)
        .ok_or_else(|| bad("natural-convection enthalpy Newton allowance overflows"))?;
    let max_krylov = max_updates
        .checked_mul(config.newton.linear_restart)
        .and_then(|value| value.checked_mul(config.newton.max_linear_cycles))
        .ok_or_else(|| bad("natural-convection enthalpy Krylov allowance overflows"))?;
    let max_endpoint_evaluations = boundary_law::MAX_ITERATIONS
        .checked_mul(2)
        .ok_or_else(|| bad("natural-convection endpoint allowance overflows"))?;
    let mut solid_solves = 0;
    let mut krylov_iterations = 0;
    let mut nonlinear_updates = 0;
    let mut radiation_trials = 0;
    let mut endpoint_evaluations = 0;
    let mut maximum_radiation_energy_residual_j = 0.0_f64;
    let mut physical_target = None;
    let initial_candidate = engine.prescribed_endpoint_initial_guess(cx, problem.boundary, old)
        .map_err(|error| match error {
            HeterogeneousEnthalpyError::Enthalpy(inner) => failed(inner),
            other => bad(format!("enthalpy initial bath state refused: {other}")),
        })?;
    for iteration in 1..=boundary_law::MAX_ITERATIONS {
        poll(cx, deadline)?;
        let replacements = laws
            .iter()
            .zip(&regions)
            .map(|(law, &region)| {
                let h = *coefficients
                    .get(&law.target)
                    .ok_or_else(|| bad("natural-convection trial lost an applied coefficient"))?;
                Ok((region, h, law.ambient_k))
            })
            .collect::<Result<Vec<_>, SolveRefusal>>()?;
        let boundary = problem
            .boundary
            .with_uniform_robin_replacements(&replacements)
            .map_err(lower)?;
        let trial_problem = ConductionProblem {
            boundary: &boundary,
            ..problem
        };
        let mut inner = config;
        if let Some(target) = physical_target {
            // This evaluates the current FROZEN natural coefficient, avoiding
            // an undefined natural-card query at an ambient initial state.
            count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
            let initial = radiative::evaluate(
                cx, engine, trial_problem, interfaces, old, &initial_candidate,
                dt_s, config, radiation, enclosure,
            )?;
            inner.newton.absolute_tolerance = config.newton.absolute_tolerance.min(target * 0.125);
            inner.newton.relative_tolerance = tighter_relative(
                config.newton.relative_tolerance,
                target * 0.125,
                initial.residual_norm_j,
            )?;
        }
        let (endpoint, radiation_evidence, ambient, enclosed) = if let Some(lowering) = radiation {
            let solved = radiative::advance(
                cx, engine, trial_problem, interfaces, old, dt_s, inner,
                lowering, enclosure, deadline,
            )?;
            count(&mut solid_solves, solved.solid_solves, max_solid_solves)?;
            count(&mut radiation_trials, solved.solid_solves, max_solid_solves)?;
            count(&mut krylov_iterations, solved.krylov_iterations, max_krylov)?;
            count(&mut nonlinear_updates, solved.updates, max_updates)?;
            maximum_radiation_energy_residual_j =
                maximum_radiation_energy_residual_j.max(solved.energy_residual_j.abs());
            (solved.endpoint, solved.evidence, solved.ambient, solved.enclosure)
        } else {
            let endpoint = engine
                .advance(cx, trial_problem, interfaces, old, dt_s, inner)
                .map_err(failed)?;
            count(&mut solid_solves, 1, max_solid_solves)?;
            count(
                &mut nonlinear_updates,
                endpoint.newton.iterations,
                max_updates,
            )?;
            for row in &endpoint.newton.history {
                count(&mut krylov_iterations, row.linear_iterations, max_krylov)?;
            }
            (endpoint, "null".to_string(), None, None)
        };
        poll(cx, deadline)?;
        let fluxes = radiative::fluxes(&endpoint, ambient.as_ref(), enclosed.as_ref());
        let mut next = BTreeMap::new();
        let mut points = Vec::with_capacity(laws.len());
        let mut worst_change = 0.0_f64;
        let mut physical_rows = Vec::with_capacity(laws.len());
        for (law, &region) in laws.iter().zip(&regions) {
            let flux = fluxes
                .iter()
                .find(|flux| flux.region == law.target)
                .ok_or_else(|| bad("natural-convection target has no solved convective trace"))?;
            let point = boundary_law::coefficient(
                law,
                finite(flux.mean_wall_temperature_k - law.ambient_k)?,
                pressure_pa,
            )?;
            let applied = *coefficients
                .get(&law.target)
                .ok_or_else(|| bad("natural-convection applied coefficient is missing"))?;
            worst_change = worst_change.max(finite((point.htc_w_m2_k - applied).abs() / applied)?);
            next.insert(law.target.clone(), point.htc_w_m2_k);
            physical_rows.push((region, point.htc_w_m2_k, law.ambient_k));
            points.push(point);
        }
        let physical_boundary = problem
            .boundary
            .with_uniform_robin_replacements(&physical_rows)
            .map_err(lower)?;
        count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
        let physical = radiative::evaluate(
            cx, engine, ConductionProblem { boundary: &physical_boundary, ..problem },
            interfaces, old, &endpoint.specific_enthalpy_j_kg, dt_s, config,
            radiation, enclosure,
        )?;
        let target = match physical_target {
            Some(target) => target,
            None => {
                // An actual initial equilibrium coefficient is undefined by
                // the card. Freeze the residual scale at this first physical
                // candidate; no later iterate may inflate the target.
                let target = finite(
                    config
                        .newton
                        .absolute_tolerance
                        .max(config.newton.relative_tolerance * physical.residual_norm_j),
                )?;
                if target <= 0.0 {
                    return Err(bad(
                        "natural-convection enthalpy residual target is not representable",
                    ));
                }
                physical_target = Some(target);
                target
            }
        };
        poll(cx, deadline)?;
        if worst_change <= boundary_law::TOLERANCE_REL
            && physical.residual_norm_j <= target
            && physical.energy_residual_j.abs() <= config.energy_tolerance_j
        {
            let converged = laws
                .iter()
                .zip(points)
                .map(|(law, coefficient)| {
                    let flux = physical
                        .convective_robin_fluxes
                        .iter()
                        .find(|flux| flux.region == law.target)
                        .ok_or_else(|| bad("physical natural-convection trace is missing"))?;
                    Ok(boundary_law::Converged {
                        law: law.clone(),
                        coefficient,
                        mean_wall_k: flux.mean_wall_temperature_k,
                        heat_rate_w: flux.heat_rate_w,
                    })
                })
                .collect::<Result<Vec<_>, SolveRefusal>>()?;
            let receipt = boundary_law::receipt_fragment(&converged, iteration)?;
            return Ok(NaturalStep {
                endpoint,
                radiation: ambient,
                enclosure: enclosed,
                radiation_evidence,
                applied_boundary: boundary,
                applied_coefficients: coefficients,
                physical_coefficients: next,
                physical_convective_robin_fluxes: physical.convective_robin_fluxes,
                receipt,
                iterations: iteration,
                solid_solves,
                krylov_iterations,
                nonlinear_updates,
                radiation_trials,
                endpoint_evaluations,
                maximum_radiation_energy_residual_j,
                physical_residual_norm_j: physical.residual_norm_j,
                physical_residual_tolerance_j: target,
                physical_energy_residual_j: physical.energy_residual_j,
                physical_convective_out_w: physical.convective_out_w,
                physical_radiative_out_w: physical.radiation_out_w,
                physical_neumann_out_w: physical.neumann_out_w,
                physical_dirichlet_in_w: physical.dirichlet_in_w,
            });
        }
        if iteration == boundary_law::MAX_ITERATIONS {
            return Err(conduction_error(
                "cli-solve-conduction-transient-natural-unconverged",
                format!(
                    "enthalpy natural convection exhausted {iteration} iterations: relative h change {worst_change:.4e}, physical residual {:.4e} J against {target:.4e} J, energy residual {:.4e} J against {:.4e} J",
                    physical.residual_norm_j, physical.energy_residual_j, config.energy_tolerance_j,
                ),
                "inspect the card regime and numerical tolerances; no partial thermal endpoint is published",
            ));
        }
        coefficients = next;
    }
    Err(bad(
        "enthalpy natural convection has no admitted iteration budget",
    ))
}


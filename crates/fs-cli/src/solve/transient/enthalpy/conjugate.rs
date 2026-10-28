//! Quasi-steady fan-network exchange over one immutable enthalpy time step.
//!
//! The shared branch-major IQN-ILS owner supplies every air update and its
//! independent branch heat gates. Applied Robin references stay attached to
//! the last counted solid response. A separate full endpoint evaluation uses
//! the air march's updated references, including nonlinear ambient radiation.
//! Each air callback can also close the existing natural-reservoir law; its
//! actual coefficients enter that same final audit, while its applied rows
//! stay attached to the counted solid response.

use std::collections::{BTreeMap, BTreeSet};

use fs_airflow::conjugate::SolidRegionState;
use fs_conduction::ThermalBoundary;

use super::super::super::conjugate as exchange;
use super::*;

pub(super) struct AirflowStep {
    pub(super) endpoint: EnthalpyStepSolution,
    pub(super) radiation: Option<RadiativeEndpoint>,
    pub(super) radiation_evidence: String,
    pub(super) natural: Option<NaturalEndpoint>,
    pub(super) natural_evidence: String,
    pub(super) natural_iterations: usize,
    pub(super) natural_solid_solves: usize,
    pub(super) natural_updates: usize,
    pub(super) natural_krylov_iterations: usize,
    pub(super) natural_endpoint_evaluations: usize,
    pub(super) outcome: exchange::ConjugateOutcome,
    pub(super) applied_boundary: ThermalBoundary,
    pub(super) applied_references: BTreeMap<String, f64>,
    pub(super) physical_references: BTreeMap<String, f64>,
    pub(super) physical_convective_robin_fluxes: Vec<fs_conduction::RobinFlux>,
    /// Every inner solid response; final-inner histories remain on endpoint.
    pub(super) solid_solves: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) krylov_iterations: usize,
    pub(super) radiation_trials: usize,
    /// Seed, old-state and final audits, plus every nested natural audit.
    pub(super) endpoint_evaluations: usize,
    pub(super) reference_tolerance_k: f64,
    pub(super) residual_reference: &'static str,
    pub(super) initial_seed_residual_norm_j: f64,
    pub(super) maximum_radiation_energy_residual_j: f64,
    pub(super) physical_residual_norm_j: f64,
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy_residual_j: f64,
    pub(super) physical_convective_out_w: f64,
    pub(super) physical_radiative_out_w: f64,
    pub(super) physical_dirichlet_in_w: f64,
    pub(super) air_heat_gain_w: f64,
    pub(super) off_path_convective_out_w: f64,
    pub(super) coupled_energy_residual_j: f64,
}

struct SolidTrial {
    endpoint: EnthalpyStepSolution,
    radiation: Option<RadiativeEndpoint>,
    radiation_evidence: String,
    natural: Option<NaturalEndpoint>,
    natural_evidence: String,
    physical_natural_boundary: Option<ThermalBoundary>,
    applied_boundary: ThermalBoundary,
    applied_references: BTreeMap<String, f64>,
    convective_out_w: f64,
    off_path_convective_out_w: f64,
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), SolveRefusal> {
    *total = total
        .checked_add(added)
        .filter(|value| *value <= cap)
        .ok_or_else(|| bad("enthalpy airflow exceeded its numerical work allowance"))?;
    Ok(())
}

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, SolveRefusal> {
    values
        .into_iter()
        .try_fold(0.0, |total, value| finite(total + value))
}

fn boundary_at(
    base: &ThermalBoundary,
    regions: &[(&str, usize, f64)],
    references: &BTreeMap<String, f64>,
) -> Result<ThermalBoundary, SolveRefusal> {
    if references.len() != regions.len() {
        return Err(bad("enthalpy airflow needs one reference per solid target"));
    }
    let replacements = regions
        .iter()
        .map(|&(target, region, coefficient)| {
            let reference = *references
                .get(target)
                .ok_or_else(|| bad("enthalpy airflow omitted a Robin reference"))?;
            Ok((region, coefficient, reference))
        })
        .collect::<Result<Vec<_>, SolveRefusal>>()?;
    base.with_uniform_robin_replacements(&replacements)
        .map_err(lower)
}

fn inner_config(
    mut config: EnthalpyStepConfig,
    allowance_j: f64,
    initial_norm_j: f64,
) -> Result<EnthalpyStepConfig, SolveRefusal> {
    config.newton.absolute_tolerance = config.newton.absolute_tolerance.min(allowance_j);
    if initial_norm_j > 0.0 {
        config.newton.relative_tolerance = config
            .newton
            .relative_tolerance
            .min(allowance_j / initial_norm_j);
    }
    // Both tolerances must remain strictly positive under the shared Newton
    // validator. Underflow is a refusal, never a hidden tolerance widening.
    if config.newton.absolute_tolerance <= 0.0 || config.newton.relative_tolerance <= 0.0 {
        return Err(bad("enthalpy airflow inner tolerance is not representable"));
    }
    config
        .newton
        .validate()
        .map_err(|error| failed(EnthalpyError::Newton(error)))?;
    Ok(config)
}

/// Every radiation, Newton and air trial receives exactly the same old h.
/// Air flow, inlet and transport properties remain the retained steady values.
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
    laws: &[buoyancy::NaturalLaw],
    pressure_pa: f64,
    path: &exchange::ConjugatePath,
    deadline: Option<(Instant, f64)>,
) -> Result<AirflowStep, SolveRefusal> {
    poll(cx, deadline)?;
    if path.segments.is_empty()
        || old.is_empty()
        || old.len() != problem.mesh.vertex_count()
        || config.max_newton_iterations == 0
        || !(dt_s.is_finite() && dt_s > 0.0)
        || !(config.energy_tolerance_j.is_finite() && config.energy_tolerance_j > 0.0)
    {
        return Err(bad(
            "enthalpy airflow requires segments, matching history and positive step controls",
        ));
    }
    config
        .newton
        .validate()
        .map_err(|error| failed(EnthalpyError::Newton(error)))?;
    if let Some(lowering) = radiation {
        lowering.config.validate().map_err(lower)?;
    }
    let mut targets = BTreeSet::new();
    let mut regions = Vec::with_capacity(path.segments.len());
    for (index, segment) in path.segments.iter().enumerate() {
        if index % 256 == 0 {
            poll(cx, deadline)?;
        }
        if !targets.insert(segment.target.as_str())
            || laws.iter().any(|law| law.target == segment.target)
            || !(segment.htc_w_m2_k.is_finite() && segment.htc_w_m2_k > 0.0)
            || !(segment.wetted_area_m2.is_finite() && segment.wetted_area_m2 > 0.0)
        {
            return Err(bad(
                "enthalpy airflow needs distinct positive-area, positive-coefficient targets, separate from natural reservoirs",
            ));
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == &segment.target)
            .ok_or_else(|| bad("enthalpy airflow target is absent from the solid boundary"))?;
        regions.push((segment.target.as_str(), region, segment.htc_w_m2_k));
    }
    let scaled_conductance_j_k = finite(
        dt_s
            * sum(
                path.segments
                    .iter()
                    .map(|segment| segment.htc_w_m2_k * segment.wetted_area_m2),
            )?,
    )?;
    if scaled_conductance_j_k <= 0.0 {
        return Err(bad(
            "enthalpy airflow step conductance must be representable and positive",
        ));
    }

    let mut exchange_config = exchange::exchange_config(false);
    let natural_factor = if laws.is_empty() { 1 } else { buoyancy::MAX_ITERATIONS };
    let max_natural_iterations = exchange_config.max_iterations.checked_mul(natural_factor)
        .ok_or_else(|| bad("enthalpy airflow/natural work allowance overflows"))?;
    let max_natural_endpoint_evaluations = max_natural_iterations.checked_mul(2)
        .ok_or_else(|| bad("enthalpy airflow/natural audit allowance overflows"))?;
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    let max_solid_solves = max_natural_iterations
        .checked_mul(radiation_factor)
        .ok_or_else(|| bad("enthalpy airflow/radiation work allowance overflows"))?;
    let max_updates = max_solid_solves
        .checked_mul(config.max_newton_iterations)
        .ok_or_else(|| bad("enthalpy airflow Newton allowance overflows"))?;
    let max_krylov = max_updates
        .checked_mul(config.newton.linear_restart)
        .and_then(|value| value.checked_mul(config.newton.max_linear_cycles))
        .ok_or_else(|| bad("enthalpy airflow Krylov allowance overflows"))?;
    let max_endpoint_evaluations = exchange_config.max_iterations.checked_add(2)
        .and_then(|air| if laws.is_empty() { Some(air) }
            else { air.checked_add(max_natural_endpoint_evaluations) })
        .ok_or_else(|| bad("enthalpy airflow endpoint allowance overflows"))?;

    // This is the same branch-inlet seed used by the shared exchange owner.
    // It is an explicit boundary guess, not an air march at the initial wall.
    // Natural rows here retain their initial applied coefficients; a powered
    // ambient start need not evaluate an undefined equilibrium correlation.
    let mut seed_references = BTreeMap::new();
    let mut seed_segments = 0;
    for branch in path.air_paths() {
        poll(cx, deadline)?;
        for segment in branch.segments() {
            let derivation = path
                .segments
                .get(seed_segments)
                .ok_or_else(|| bad("enthalpy airflow branch ordering has extra segments"))?;
            if segment.region() != derivation.target.as_str()
                || segment.area_m2() != derivation.wetted_area_m2
                || segment.htc_w_per_m2_k() != derivation.htc_w_m2_k
                || seed_references
                    .insert(segment.region().to_string(), branch.inlet_temperature_k())
                    .is_some()
            {
                return Err(bad(
                    "enthalpy airflow branch and solid segment declarations disagree",
                ));
            }
            count(&mut seed_segments, 1, path.segments.len())?;
        }
    }
    if seed_segments != path.segments.len() {
        return Err(bad("enthalpy airflow branch ordering omits solid segments"));
    }
    let patches = radiation.map(|lowering| lowering.patches.as_slice());
    let mut endpoint_evaluations = 0;
    // Air-reference updates replace Robin rows only; the prescribed trace is
    // unchanged across the seed and every inner coupling trial.
    let initial_candidate = engine.prescribed_endpoint_initial_guess(cx, problem.boundary, old)
        .map_err(|error| match error {
            HeterogeneousEnthalpyError::Enthalpy(inner) => failed(inner),
            other => bad(format!("enthalpy initial bath state refused: {other}")),
        })?;
    let initial_seed_residual_norm_j = {
        let seed_boundary = boundary_at(problem.boundary, &regions, &seed_references)?;
        count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
        engine
            .evaluate_endpoint(
                cx,
                ConductionProblem { boundary: &seed_boundary, ..problem },
                interfaces,
                old,
                &initial_candidate,
                dt_s,
                config,
                patches,
            )
            .map_err(failed)?
            .residual_norm_j
    };
    let physical_residual_tolerance_j = finite(
        config
            .newton
            .absolute_tolerance
            .max(config.newton.relative_tolerance * initial_seed_residual_norm_j),
    )?;
    let inner_allowance_j = finite(
        0.125
            * physical_residual_tolerance_j
                .min(config.energy_tolerance_j / fs_math::det::sqrt(old.len() as f64)),
    )?;
    // Uniform reference changes add dt*h*integral(N_i) times that change to
    // the nodal balance. Its Euclidean norm is bounded by dt*sum(hA)*max(dT).
    // Reserve room independently in both the full residual and joule budget.
    exchange_config.temperature_tolerance_k = exchange_config
        .temperature_tolerance_k
        .min(0.01 * config.energy_tolerance_j / scaled_conductance_j_k)
        .min(0.125 * physical_residual_tolerance_j / scaled_conductance_j_k);
    if !(inner_allowance_j.is_finite() && inner_allowance_j > 0.0)
        || !(exchange_config.temperature_tolerance_k.is_finite()
            && exchange_config.temperature_tolerance_k > 0.0)
    {
        return Err(bad(
            "enthalpy airflow cannot represent positive residual and reference tolerances",
        ));
    }

    let mut air_iterations = 0;
    let mut solid_solves = 0;
    let mut nonlinear_updates = 0;
    let mut krylov_iterations = 0;
    let mut radiation_trials = 0;
    let mut natural_iterations = 0;
    let mut natural_solid_solves = 0;
    let mut natural_updates = 0;
    let mut natural_krylov_iterations = 0;
    let mut natural_endpoint_evaluations = 0;
    let mut maximum_radiation_energy_residual_j = 0.0_f64;
    let mut last = None;
    let outcome = exchange::run_exchange_with_config(cx, path, exchange_config, |cx, references| {
        poll(cx, deadline)?;
        count(&mut air_iterations, 1, exchange_config.max_iterations)?;
        let mut boundary = boundary_at(problem.boundary, &regions, references)?;
        let trial_problem = ConductionProblem { boundary: &boundary, ..problem };
        count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
        let initial_norm_j = engine
            .evaluate_endpoint(
                cx, trial_problem, interfaces, old, &initial_candidate, dt_s, config, patches,
            )
            .map_err(failed)?
            .residual_norm_j;
        let inner = inner_config(config, inner_allowance_j, initial_norm_j)?;
        let mut natural_endpoint = None;
        let mut natural_evidence = "null".to_string();
        let mut physical_natural_boundary = None;
        let (endpoint, radiation_evidence, radiative) = if !laws.is_empty() {
            let mut natural_config = inner;
            // Reserve room below the original endpoint energy and residual
            // allowances for the remaining air-reference discrepancy.
            natural_config.energy_tolerance_j = finite(config.energy_tolerance_j * 0.125)?;
            let solved = natural::advance(
                cx, engine, trial_problem, interfaces, old, dt_s, natural_config,
                radiation, laws, pressure_pa, deadline,
            )?;
            natural_evidence = natural_step_row(&solved)?;
            count(&mut natural_iterations, solved.iterations, max_natural_iterations)?;
            count(&mut natural_solid_solves, solved.solid_solves, max_solid_solves)?;
            count(&mut natural_updates, solved.nonlinear_updates, max_updates)?;
            count(&mut natural_krylov_iterations, solved.krylov_iterations, max_krylov)?;
            count(&mut natural_endpoint_evaluations, solved.endpoint_evaluations,
                max_natural_endpoint_evaluations)?;
            count(&mut endpoint_evaluations, solved.endpoint_evaluations, max_endpoint_evaluations)?;
            count(&mut solid_solves, solved.solid_solves, max_solid_solves)?;
            count(&mut nonlinear_updates, solved.nonlinear_updates, max_updates)?;
            count(&mut krylov_iterations, solved.krylov_iterations, max_krylov)?;
            count(&mut radiation_trials, solved.radiation_trials, max_solid_solves)?;
            maximum_radiation_energy_residual_j = maximum_radiation_energy_residual_j
                .max(solved.maximum_radiation_energy_residual_j);
            let rows = laws.iter().map(|law| {
                let region = boundary.region_names().iter()
                    .position(|name| name == &law.target)
                    .ok_or_else(|| bad("enthalpy natural target disappeared"))?;
                let h = *solved.physical_coefficients.get(&law.target)
                    .ok_or_else(|| bad("enthalpy natural coefficient disappeared"))?;
                Ok((region, h, law.ambient_k))
            }).collect::<Result<Vec<_>, SolveRefusal>>()?;
            physical_natural_boundary = Some(
                boundary.with_uniform_robin_replacements(&rows).map_err(lower)?,
            );
            boundary = solved.applied_boundary;
            natural_endpoint = Some(NaturalEndpoint {
                receipt: solved.receipt,
                applied_coefficients: solved.applied_coefficients,
                boundary: boundary.clone(),
            });
            (solved.endpoint, solved.radiation_evidence, solved.radiation)
        } else if let Some(lowering) = radiation {
            let solved = engine
                .advance_with_ambient_radiation(
                    cx, trial_problem, interfaces, old, dt_s, inner,
                    &lowering.patches, lowering.config,
                )
                .map_err(failed)?;
            count(&mut solid_solves, solved.radiation.iterations, max_solid_solves)?;
            count(&mut radiation_trials, solved.radiation.iterations, max_solid_solves)?;
            count(&mut nonlinear_updates, solved.radiation.solid_iterations, max_updates)?;
            count(&mut krylov_iterations, solved.radiation.krylov_iterations, max_krylov)?;
            maximum_radiation_energy_residual_j =
                maximum_radiation_energy_residual_j.max(solved.physical_energy_residual_j.abs());
            let evidence = radiation_step_row(&solved)?;
            let radiative = RadiativeEndpoint {
                combined_boundary: solved.combined_boundary,
                convective_robin_fluxes: solved.convective_robin_fluxes,
                convective_out_w: solved.convective_out_w,
                report: solved.radiation,
            };
            (solved.conduction, evidence, Some(radiative))
        } else {
            let endpoint = engine
                .advance(cx, trial_problem, interfaces, old, dt_s, inner)
                .map_err(failed)?;
            count(&mut solid_solves, 1, max_solid_solves)?;
            count(&mut nonlinear_updates, endpoint.newton.iterations, max_updates)?;
            for row in &endpoint.newton.history {
                count(&mut krylov_iterations, row.linear_iterations, max_krylov)?;
            }
            (endpoint, "null".to_string(), None)
        };
        poll(cx, deadline)?;
        let (fluxes, convective_out_w) = radiative.as_ref().map_or(
            (endpoint.robin_fluxes.as_slice(), endpoint.robin_out_w),
            |result| (result.convective_robin_fluxes.as_slice(), result.convective_out_w),
        );
        // Only original convection reaches the air owner. Radiation on the
        // same trace exchanges with its own reservoir and is excluded here.
        let states = path
            .segments
            .iter()
            .map(|segment| {
                fluxes
                    .iter()
                    .find(|flux| flux.region == segment.target)
                    .map(SolidRegionState::from_robin_flux)
                    .ok_or_else(|| bad("enthalpy airflow target has no solved convective trace"))
            })
            .collect::<Result<Vec<_>, SolveRefusal>>()?;
        let off_path_convective_out_w = sum(
            fluxes
                .iter()
                .filter(|flux| !targets.contains(flux.region.as_str()))
                .map(|flux| flux.heat_rate_w),
        )?;
        last = Some(SolidTrial {
            endpoint,
            radiation: radiative,
            radiation_evidence,
            natural: natural_endpoint,
            natural_evidence,
            physical_natural_boundary,
            applied_boundary: boundary,
            applied_references: references.clone(),
            convective_out_w,
            off_path_convective_out_w,
        });
        Ok(states)
    })?;
    poll(cx, deadline)?;
    if outcome.solution.iterations != air_iterations {
        return Err(bad(
            "enthalpy airflow work does not match the shared exchange history",
        ));
    }
    let last = last.ok_or_else(|| bad("enthalpy airflow accepted no solid response"))?;
    // Standard branch receipts describe exactly the applied solid response.
    let outcome = exchange::cross_check_decomposition(
        outcome, last.convective_out_w, last.off_path_convective_out_w,
    )?;
    let mut physical_references = BTreeMap::new();
    for branch in &outcome.solution.branches {
        for segment in &branch.march.segments {
            if !targets.contains(segment.region.as_str())
                || physical_references
                    .insert(segment.region.clone(), segment.reference_temperature_k)
                    .is_some()
            {
                return Err(bad(
                    "enthalpy airflow physical references do not partition the targets",
                ));
            }
        }
    }
    let physical_boundary = boundary_at(
        last.physical_natural_boundary.as_ref().unwrap_or(problem.boundary),
        &regions, &physical_references,
    )?;
    count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
    let physical = engine
        .evaluate_endpoint(
            cx,
            ConductionProblem { boundary: &physical_boundary, ..problem },
            interfaces,
            old,
            &last.endpoint.specific_enthalpy_j_kg,
            dt_s,
            config,
            patches,
        )
        .map_err(failed)?;
    poll(cx, deadline)?;
    let off_path_convective_out_w = sum(
        physical
            .convective_robin_fluxes
            .iter()
            .filter(|flux| !targets.contains(flux.region.as_str()))
            .map(|flux| flux.heat_rate_w),
    )?;
    let air_heat_gain_w = sum(
        outcome.solution.branches.iter().map(|branch| branch.march.total_heat_rate_w),
    )?;
    let net_input_w = finite(
        physical.source_w + physical.dirichlet_in_w - physical.neumann_out_w - off_path_convective_out_w
            - air_heat_gain_w - physical.radiation_out_w,
    )?;
    let coupled_energy_residual_j =
        finite(physical.stored_energy_change_j - dt_s * net_input_w)?;
    if physical.residual_norm_j > physical_residual_tolerance_j
        || physical.energy_residual_j.abs() > config.energy_tolerance_j
        || coupled_energy_residual_j.abs() > config.energy_tolerance_j
    {
        return Err(conduction_error(
            "cli-solve-conduction-transient-airflow-energy",
            format!(
                "enthalpy airflow endpoint has physical residual {:.4e} J against {:.4e} J, physical energy residual {:.4e} J and actual-air energy residual {coupled_energy_residual_j:.4e} J against {:.4e} J",
                physical.residual_norm_j, physical_residual_tolerance_j,
                physical.energy_residual_j, config.energy_tolerance_j,
            ),
            "inspect the thermal tolerances, air paths and reservoir exchanges; no endpoint is published",
        ));
    }
    Ok(AirflowStep {
        endpoint: last.endpoint,
        radiation: last.radiation,
        radiation_evidence: last.radiation_evidence,
        natural: last.natural,
        natural_evidence: last.natural_evidence,
        natural_iterations, natural_solid_solves, natural_updates,
        natural_krylov_iterations, natural_endpoint_evaluations,
        outcome,
        applied_boundary: last.applied_boundary,
        applied_references: last.applied_references,
        physical_references,
        physical_convective_robin_fluxes: physical.convective_robin_fluxes,
        solid_solves,
        nonlinear_updates,
        krylov_iterations,
        radiation_trials,
        endpoint_evaluations,
        reference_tolerance_k: exchange_config.temperature_tolerance_k,
        residual_reference: "initial-enthalpy-with-branch-inlet-seed-boundary",
        initial_seed_residual_norm_j,
        maximum_radiation_energy_residual_j,
        physical_residual_norm_j: physical.residual_norm_j,
        physical_residual_tolerance_j,
        physical_energy_residual_j: physical.energy_residual_j,
        physical_convective_out_w: physical.convective_out_w,
        physical_radiative_out_w: physical.radiation_out_w,
        physical_dirichlet_in_w: physical.dirichlet_in_w,
        air_heat_gain_w,
        off_path_convective_out_w,
        coupled_energy_residual_j,
    })
}

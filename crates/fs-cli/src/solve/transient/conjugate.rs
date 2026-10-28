//! One transient solid endpoint coupled to frozen, quasi-steady air paths.
//!
//! Every air and radiation trial uses the same immutable physical history.
//! Only the last shared-solid response is retained; updated air references do
//! not trigger an extra, uncounted solve or replace the applied boundary.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use fs_airflow::conjugate::SolidRegionState;
use fs_conduction::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, StepConfig, StepSolution,
};
use fs_conduction::{ConductionProblem, ThermalBoundary, ThermalInterfaces};
use fs_exec::Cx;

use super::super::{SolveRefusal, conduction_error, conjugate as exchange, radiation};
use super::{PhysicalEnergy, RadiativeEndpoint, bad, lower, nonlinear_row, poll, radiative};

pub(super) struct AirflowStep {
    pub(super) endpoint: StepSolution,
    /// Final inner solve evidence; cumulative work is retained separately.
    pub(super) nonlinear_evidence: String,
    pub(super) radiation_evidence: String,
    pub(super) radiation: Option<RadiativeEndpoint>,
    pub(super) outcome: exchange::ConjugateOutcome,
    /// References the last solid actually used, before the air-law update.
    pub(super) applied_references: BTreeMap<String, f64>,
    pub(super) applied_boundary: ThermalBoundary,
    /// All BE solves, including every radiative trial of every air trial.
    pub(super) solid_solves: usize,
    pub(super) krylov_iterations: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
    pub(super) radiation_trials: usize,
    pub(super) reference_tolerance_k: f64,
    pub(super) maximum_physical_energy_residual_j: f64,
    /// Direct initial and accepted physical endpoint evaluations.
    pub(super) endpoint_evaluations: usize,
    pub(super) initial_physical_residual_norm_j: f64,
    pub(super) physical_residual_norm_j: f64,
    /// Fixed from the actual air law at the old free/new prescribed predictor.
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy_residual_j: f64,
    pub(super) physical_energy: PhysicalEnergy,
    pub(super) air_heat_gain_w: f64,
    pub(super) off_path_convective_out_w: f64,
    pub(super) radiation_out_w: f64,
    pub(super) physical_dirichlet_in_w: f64,
    pub(super) coupled_energy_residual_j: f64,
}

struct SolidTrial {
    endpoint: StepSolution,
    nonlinear_evidence: String,
    radiation_evidence: String,
    radiation: Option<RadiativeEndpoint>,
    applied_references: BTreeMap<String, f64>,
    applied_boundary: ThermalBoundary,
    convective_out_w: f64,
    off_path_convective_out_w: f64,
}

fn finite(value: f64, what: &str) -> Result<f64, SolveRefusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!("transient airflow {what} is not finite")))
    }
}

fn total(values: impl IntoIterator<Item = f64>, what: &str) -> Result<f64, SolveRefusal> {
    values
        .into_iter()
        .try_fold(0.0, |sum, value| finite(sum + value, what))
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), SolveRefusal> {
    *total = total
        .checked_add(added)
        .filter(|value| *value <= cap)
        .ok_or_else(|| bad("transient airflow exceeded its admitted numerical work allowance"))?;
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
            "transient airflow cannot represent its inner residual tolerance",
        ))
    }
}

/// Evaluate air at a supplied wall field without a solid solve or derivative.
/// Natural coefficients remain the outer owner's supplied frozen rows.
#[allow(clippy::too_many_arguments)]
fn predictor_air_boundary(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    temperature: &[f64],
    path: &exchange::ConjugatePath,
    regions: &[(&str, usize, f64)],
    deadline: Option<(Instant, f64)>,
) -> Result<(ThermalBoundary, Vec<f64>), SolveRefusal> {
    let mut areas = vec![0.0; regions.len()];
    let mut means = vec![0.0; regions.len()];
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 {
            poll(cx, deadline)?;
        }
        let Some(index) = problem
            .boundary
            .region_for(slot)
            .and_then(|region| regions.iter().position(|&(_, found, _)| found == region))
        else {
            continue;
        };
        // Reproduce the forward Robin report, including its fused area sum.
        let mut wall_sum = 0.0;
        for &vertex in &face.vertices {
            wall_sum = finite(
                wall_sum + temperature[vertex as usize],
                "predictor wall sum",
            )?;
        }
        areas[index] = finite(areas[index] + face.area, "predictor trace area")?;
        means[index] = finite(
            face.area.mul_add(wall_sum / 3.0, means[index]),
            "predictor wall integral",
        )?;
    }
    for (area, mean) in areas.iter().zip(&mut means) {
        if *area <= 0.0 {
            return Err(bad("transient airflow predictor has an empty wall trace"));
        }
        *mean = finite(*mean * (1.0 / *area), "predictor wall mean")?;
    }
    let mut replacements = Vec::with_capacity(regions.len());
    let mut start = 0usize;
    for branch in path.air_paths() {
        poll(cx, deadline)?;
        let end = start
            .checked_add(branch.segments().len())
            .filter(|&end| end <= regions.len())
            .ok_or_else(|| bad("transient airflow predictor has extra branch segments"))?;
        let march = branch.march(&means[start..end]).map_err(|error| {
            bad(format!(
                "transient airflow predictor march refused: {error}"
            ))
        })?;
        if march.segments.len() != branch.segments().len() {
            return Err(bad(
                "transient airflow predictor lost its marched segment layout",
            ));
        }
        for (local, (segment, state)) in branch.segments().iter().zip(&march.segments).enumerate() {
            let index = start + local;
            let (target, region, h) = regions[index];
            let area = segment.area_m2();
            if segment.region() != target
                || state.region != target
                || segment.htc_w_per_m2_k().to_bits() != h.to_bits()
                || (areas[index] - area).abs() > 128.0 * f64::EPSILON * areas[index].max(area)
            {
                return Err(bad(
                    "transient airflow predictor differs from the physical branch/solid partition",
                ));
            }
            replacements.push((region, h, state.reference_temperature_k));
        }
        start = end;
    }
    if start != regions.len() {
        return Err(bad("transient airflow predictor omits branch segments"));
    }
    let references = replacements
        .iter()
        .map(|&(_, _, reference)| reference)
        .collect();
    let boundary = problem
        .boundary
        .with_uniform_robin_replacements(&replacements)
        .map_err(lower)?;
    Ok((boundary, references))
}

/// Advance one common solid step and independently march every supplied air
/// branch. Air has no stored state here: its flow and inlet remain frozen.
/// Raw reference and branch-local watt gates come from the existing owner.
/// The unchanged endpoint must also pass its complete updated-reference nodal
/// residual and energy gates, with actual air enthalpy gain kept separate.
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
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    path: &exchange::ConjugatePath,
    deadline: Option<(Instant, f64)>,
) -> Result<AirflowStep, SolveRefusal> {
    poll(cx, deadline)?;
    if old.is_empty()
        || old.len() != problem.mesh.vertex_count()
        || !(dt_s.is_finite()
            && dt_s > 0.0
            && step_config.energy_tolerance_j.is_finite()
            && step_config.energy_tolerance_j > 0.0)
    {
        return Err(bad(
            "transient airflow needs matching history and positive finite step/energy controls",
        ));
    }
    if path.segments.is_empty() {
        return Err(bad("transient airflow has no lowered air segments"));
    }
    let mut targets = BTreeSet::new();
    let mut regions = Vec::with_capacity(path.segments.len());
    let conductance_w_k = total(
        path.segments
            .iter()
            .map(|segment| segment.htc_w_m2_k * segment.wetted_area_m2),
        "total interface conductance",
    )?;
    for segment in &path.segments {
        if !targets.insert(segment.target.as_str())
            || !(segment.htc_w_m2_k.is_finite()
                && segment.htc_w_m2_k > 0.0
                && segment.wetted_area_m2.is_finite()
                && segment.wetted_area_m2 > 0.0)
        {
            return Err(bad(
                "transient airflow needs distinct positive-area, positive-coefficient segments",
            ));
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == &segment.target)
            .ok_or_else(|| {
                bad(format!(
                    "transient airflow target '{}' is absent from the solid boundary",
                    segment.target,
                ))
            })?;
        regions.push((segment.target.as_str(), region, segment.htc_w_m2_k));
    }
    let scaled_conductance_j_k =
        finite(dt_s * conductance_w_k, "step-scaled interface conductance")?;
    let mut config = exchange::exchange_config(false);
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    if let Some(lowering) = radiation {
        lowering.config.validate().map_err(lower)?;
    }
    if let Some(policy) = nonlinear {
        policy.validate().map_err(lower)?;
    }
    let max_endpoint_evaluations = 2usize;
    let mut endpoint_evaluations = 0usize;
    let max_solid_solves = config
        .max_iterations
        .checked_mul(radiation_factor)
        .ok_or_else(|| bad("transient airflow/radiation trial allowance overflows"))?;
    let max_krylov = max_solid_solves
        .checked_mul(step_config.linear.max_iterations)
        .ok_or_else(|| bad("transient airflow Krylov allowance overflows"))?;
    let max_updates = max_solid_solves
        .checked_mul(nonlinear.map_or(0, |p| p.max_iterations))
        .ok_or_else(|| bad("transient airflow Newton allowance overflows"))?;
    let max_backtracks = max_updates
        .checked_mul(nonlinear.map_or(0, |p| p.line_search.max_backtracks))
        .ok_or_else(|| bad("transient airflow backtrack allowance overflows"))?;
    if let Some(policy) = nonlinear {
        policy
            .line_search
            .max_backtracks
            .checked_add(1)
            .and_then(|trials| trials.checked_mul(max_updates))
            .ok_or_else(|| bad("transient airflow Newton trial allowance overflows"))?;
    }
    let mut predictor = old.to_vec();
    for &(vertex, prescribed) in problem.boundary.dirichlet() {
        *predictor
            .get_mut(vertex)
            .ok_or_else(|| bad("transient airflow prescribed vertex lies outside its history"))? =
            prescribed;
    }
    // Unlike natural convection, the air law is defined at the old predictor,
    // including equilibrium. Freeze the complete physical scale here, before
    // any outer reference iterate can inflate it.
    let (initial_physical_residual_norm_j, free_dofs, predictor_references_k) = {
        let (boundary, references) =
            predictor_air_boundary(cx, problem, &predictor, path, &regions, deadline)?;
        count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
        let initial = radiative::evaluate(
            cx, engine, ConductionProblem { boundary: &boundary, ..problem },
            interfaces, old, &predictor, dt_s, step_config, radiation, enclosure, deadline,
        )?;
        (initial.residual_norm_j, initial.free_dofs, references)
    };
    let absolute = 0.01 * step_config.energy_tolerance_j / fs_math::det::sqrt(free_dofs as f64);
    let relative = nonlinear.map_or(step_config.linear.tolerance, |policy| policy.residual_rtol);
    let physical_residual_tolerance_j = finite(
        absolute + relative * initial_physical_residual_norm_j,
        "physical residual target",
    )?;
    if physical_residual_tolerance_j <= 0.0 {
        return Err(bad(
            "transient airflow physical residual target is not representable",
        ));
    }
    // The Euclidean nodal change from reference mismatch is bounded by
    // dt*sum(hA)*max(dT_ref). Keep one eighth of the fixed physical target and
    // at most one percent of the original joule allowance for this mismatch.
    config.temperature_tolerance_k = config
        .temperature_tolerance_k
        .min((0.01 * step_config.energy_tolerance_j) / scaled_conductance_j_k)
        .min((0.125 * physical_residual_tolerance_j) / scaled_conductance_j_k);
    if !(scaled_conductance_j_k > 0.0
        && config.temperature_tolerance_k.is_finite()
        && config.temperature_tolerance_k > 0.0)
    {
        return Err(bad(
            "transient airflow cannot represent a positive physical reference tolerance",
        ));
    }
    let mut solid_solves = 0;
    let mut krylov_iterations = 0;
    let mut nonlinear_updates = 0;
    let mut nonlinear_backtracks = 0;
    let mut radiation_trials = 0;
    let mut maximum_physical_energy_residual_j = 0.0_f64;
    let mut last = None;
    let outcome = exchange::run_exchange_with_config(cx, path, config, |cx, references| {
        poll(cx, deadline)?;
        let replacements = regions
            .iter()
            .map(|&(target, region, coefficient)| {
                let reference = *references
                    .get(target)
                    .ok_or_else(|| bad("transient airflow omitted an applied Robin reference"))?;
                Ok((region, coefficient, reference))
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
        // All material/Robin/radiation coefficients at this one predictor
        // remain fixed across air trials. Only original convection reference
        // loads change. Their Euclidean residual change is bounded by
        // dt*sum(hA)*max(dT_ref), including on radiating or prescribed traces.
        // Use the conservative bound to tighten relative solves; acceptance
        // still reassembles the actual final physical residual independently.
        let mut reference_change_k = 0.0_f64;
        for (&(target, _, _), &initial_reference) in regions.iter().zip(&predictor_references_k) {
            let reference = *references
                .get(target)
                .ok_or_else(|| bad("transient airflow omitted an inner predictor reference"))?;
            reference_change_k = reference_change_k.max(finite(
                (reference - initial_reference).abs(),
                "predictor reference change",
            )?);
        }
        let inner_initial_norm_j = finite(
            initial_physical_residual_norm_j + scaled_conductance_j_k * reference_change_k,
            "inner predictor residual bound",
        )?;
        let mut inner_step = step_config;
        let mut inner_nonlinear = nonlinear;
        if let Some(policy) = inner_nonlinear.as_mut() {
            policy.residual_atol_j = policy
                .residual_atol_j
                .min(0.125 * physical_residual_tolerance_j);
            policy.residual_rtol = tighter_relative(
                policy.residual_rtol,
                0.125 * physical_residual_tolerance_j,
                inner_initial_norm_j,
            )?;
        } else {
            inner_step.linear.tolerance = tighter_relative(
                step_config.linear.tolerance,
                0.125 * physical_residual_tolerance_j,
                inner_initial_norm_j,
            )?;
        }
        if radiation.is_some() {
            // The radiation owner's physical gate includes its own joule
            // rounding floor. Reserve space for it beneath this outer target.
            inner_step.energy_tolerance_j = 0.125 * step_config.energy_tolerance_j;
            if !(inner_step.energy_tolerance_j.is_finite() && inner_step.energy_tolerance_j > 0.0) {
                return Err(bad(
                    "transient airflow inner radiation allowance is not representable",
                ));
            }
        }
        let (
            endpoint,
            nonlinear_evidence,
            radiation_evidence,
            radiative,
            physical_energy_residual_j,
        ) = if let Some(lowering) = radiation {
            let solved = radiative::advance(
                cx, engine, trial_problem, interfaces, old, dt_s, inner_step,
                inner_nonlinear, lowering, enclosure, deadline,
            )?;
            count(&mut solid_solves, solved.solid_solves, max_solid_solves)?;
            count(&mut radiation_trials, solved.solid_solves, max_solid_solves)?;
            count(&mut krylov_iterations, solved.krylov_iterations, max_krylov)?;
            count(&mut nonlinear_updates, solved.nonlinear_updates, max_updates)?;
            count(&mut nonlinear_backtracks, solved.nonlinear_backtracks, max_backtracks)?;
            (
                solved.conduction,
                solved.nonlinear_evidence,
                solved.radiation_evidence,
                Some(solved.radiation),
                solved.physical_energy.energy_residual_j,
            )
        } else if let Some(policy) = inner_nonlinear {
            let solved = engine
                .advance_nonlinear_prescribed(
                    cx,
                    trial_problem,
                    interfaces,
                    old,
                    dt_s,
                    inner_step,
                    policy,
                )
                .map_err(lower)?;
            count(&mut solid_solves, 1, max_solid_solves)?;
            count(
                &mut krylov_iterations,
                solved.step.krylov_iterations,
                max_krylov,
            )?;
            count(
                &mut nonlinear_updates,
                solved.nonlinear_iterations,
                max_updates,
            )?;
            count(&mut nonlinear_backtracks, solved.backtracks, max_backtracks)?;
            let row = nonlinear_row(&solved)?;
            let energy = solved.step.energy_residual_j;
            (solved.step, row, "null".to_string(), None, energy)
        } else {
            let endpoint = engine
                .advance_prescribed(cx, trial_problem, interfaces, old, dt_s, inner_step)
                .map_err(lower)?;
            count(&mut solid_solves, 1, max_solid_solves)?;
            count(
                &mut krylov_iterations,
                endpoint.krylov_iterations,
                max_krylov,
            )?;
            let energy = endpoint.energy_residual_j;
            (
                endpoint,
                "null".to_string(),
                "null".to_string(),
                None,
                energy,
            )
        };
        poll(cx, deadline)?;
        maximum_physical_energy_residual_j =
            maximum_physical_energy_residual_j.max(physical_energy_residual_j.abs());
        let (fluxes, convective_out_w) = radiative.as_ref().map_or(
            (endpoint.robin_fluxes.as_slice(), endpoint.robin_out_w),
            |result| {
                (
                    result.convective_robin_fluxes(),
                    result.convective_out_w(),
                )
            },
        );
        // Only original convection reaches the air marcher. Ambient radiation
        // exchanges with its own reservoir even on a shared exterior trace.
        let states = path
            .segments
            .iter()
            .map(|segment| {
                fluxes
                    .iter()
                    .find(|flux| flux.region == segment.target)
                    .map(SolidRegionState::from_robin_flux)
                    .ok_or_else(|| {
                        bad(format!(
                            "transient airflow target '{}' has no solved convective flux",
                            segment.target,
                        ))
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let off_path_convective_out_w = total(
            fluxes
                .iter()
                .filter(|flux| !targets.contains(flux.region.as_str()))
                .map(|flux| flux.heat_rate_w),
            "off-path convective heat",
        )?;
        last = Some(SolidTrial {
            endpoint,
            nonlinear_evidence,
            radiation_evidence,
            radiation: radiative,
            applied_references: references.clone(),
            applied_boundary: boundary,
            convective_out_w,
            off_path_convective_out_w,
        });
        Ok(states)
    })?;
    poll(cx, deadline)?;
    let last = last.ok_or_else(|| bad("transient airflow accepted no solid response"))?;
    let outcome = exchange::cross_check_decomposition(
        outcome,
        last.convective_out_w,
        last.off_path_convective_out_w,
    )?;
    let references = &outcome.solution.reference_temperatures_k;
    if references.len() != regions.len() {
        return Err(bad(
            "transient airflow accepted references lost their physical segment layout",
        ));
    }
    let physical_rows = regions
        .iter()
        .zip(references)
        .map(|(&(_, region, h), &reference)| (region, h, reference))
        .collect::<Vec<_>>();
    let physical_boundary = problem
        .boundary
        .with_uniform_robin_replacements(&physical_rows)
        .map_err(lower)?;
    count(&mut endpoint_evaluations, 1, max_endpoint_evaluations)?;
    let physical = radiative::evaluate(
        cx, engine, ConductionProblem { boundary: &physical_boundary, ..problem },
        interfaces, old, &last.endpoint.temperature, dt_s, step_config,
        radiation, enclosure, deadline,
    )?;
    let off_path_convective_out_w = total(
        physical
            .convective_robin_fluxes
            .iter()
            .filter(|flux| !targets.contains(flux.region.as_str()))
            .map(|flux| flux.heat_rate_w),
        "physical off-path convection",
    )?;
    let air_heat_gain_w = total(
        outcome
            .solution
            .branches
            .iter()
            .map(|branch| branch.march.total_heat_rate_w),
        "air enthalpy gain",
    )?;
    let net_input_w = finite(
        physical.source_w + physical.dirichlet_in_w
            - physical.neumann_out_w
            - off_path_convective_out_w
            - air_heat_gain_w
            - physical.radiation_out_w,
        "coupled endpoint net input",
    )?;
    let coupled_energy_residual_j = finite(
        physical.stored_energy_change_j - dt_s * net_input_w,
        "coupled endpoint energy residual",
    )?;
    if physical.residual_norm_j > physical_residual_tolerance_j
        || physical.energy_residual_j.abs() > step_config.energy_tolerance_j
        || coupled_energy_residual_j.abs() > step_config.energy_tolerance_j
    {
        return Err(conduction_error(
            "cli-solve-conduction-transient-airflow-energy",
            format!(
                "transient airflow has physical residual {:.4e} J against {:.4e} J, physical energy residual {:.4e} J and actual-air energy residual {coupled_energy_residual_j:.4e} J against {:.4e} J",
                physical.residual_norm_j,
                physical_residual_tolerance_j,
                physical.energy_residual_j,
                step_config.energy_tolerance_j,
            ),
            "tighten the thermal solve and inspect the interface heat balance; no endpoint is published",
        ));
    }
    maximum_physical_energy_residual_j =
        maximum_physical_energy_residual_j.max(physical.energy_residual_j.abs());
    Ok(AirflowStep {
        endpoint: last.endpoint,
        nonlinear_evidence: last.nonlinear_evidence,
        radiation_evidence: last.radiation_evidence,
        radiation: last.radiation,
        outcome,
        applied_references: last.applied_references,
        applied_boundary: last.applied_boundary,
        solid_solves,
        krylov_iterations,
        nonlinear_updates,
        nonlinear_backtracks,
        radiation_trials,
        reference_tolerance_k: config.temperature_tolerance_k,
        maximum_physical_energy_residual_j,
        endpoint_evaluations,
        initial_physical_residual_norm_j,
        physical_residual_norm_j: physical.residual_norm_j,
        physical_residual_tolerance_j,
        physical_energy_residual_j: physical.energy_residual_j,
        physical_energy: PhysicalEnergy::physical(&physical),
        air_heat_gain_w,
        off_path_convective_out_w,
        radiation_out_w: physical.radiation_out_w,
        physical_dirichlet_in_w: physical.dirichlet_in_w,
        coupled_energy_residual_j,
    })
}


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

use super::super::{conduction_error, conjugate as exchange, radiation, SolveRefusal};
use super::{bad, lower, nonlinear_row, poll, radiation_row, RadiativeEndpoint};

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
    radiation_out_w: f64,
    physical_dirichlet_in_w: f64,
}

fn finite(value: f64, what: &str) -> Result<f64, SolveRefusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad(format!("transient airflow {what} is not finite")))
    }
}

fn total(values: impl IntoIterator<Item = f64>, what: &str) -> Result<f64, SolveRefusal> {
    values.into_iter().try_fold(0.0, |sum, value| finite(sum + value, what))
}

fn count(total: &mut usize, added: usize, cap: usize) -> Result<(), SolveRefusal> {
    *total = total.checked_add(added)
        .filter(|value| *value <= cap)
        .ok_or_else(|| bad("transient airflow exceeded its admitted numerical work allowance"))?;
    Ok(())
}

/// Advance one common solid step and independently march every supplied air
/// branch. Air has no stored state here: its flow and inlet remain frozen.
/// Raw reference and branch-local watt gates come from the existing owner;
/// the final storage gate additionally uses actual air enthalpy gain.
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
    path: &exchange::ConjugatePath,
    deadline: Option<(Instant, f64)>,
) -> Result<AirflowStep, SolveRefusal> {
    poll(cx, deadline)?;
    if !(dt_s.is_finite() && dt_s > 0.0
        && step_config.energy_tolerance_j.is_finite() && step_config.energy_tolerance_j > 0.0)
    {
        return Err(bad("transient airflow needs positive finite step time and energy tolerance"));
    }
    if path.segments.is_empty() {
        return Err(bad("transient airflow has no lowered air segments"));
    }
    let mut targets = BTreeSet::new();
    let mut regions = Vec::with_capacity(path.segments.len());
    let conductance_w_k = total(path.segments.iter().map(|segment| {
        segment.htc_w_m2_k * segment.wetted_area_m2
    }), "total interface conductance")?;
    for segment in &path.segments {
        if !targets.insert(segment.target.as_str())
            || !(segment.htc_w_m2_k.is_finite() && segment.htc_w_m2_k > 0.0
                && segment.wetted_area_m2.is_finite() && segment.wetted_area_m2 > 0.0)
        {
            return Err(bad("transient airflow needs distinct positive-area, positive-coefficient segments"));
        }
        let region = problem.boundary.region_names().iter()
            .position(|name| name == &segment.target)
            .ok_or_else(|| bad(format!(
                "transient airflow target '{}' is absent from the solid boundary", segment.target,
            )))?;
        regions.push((segment.target.as_str(), region, segment.htc_w_m2_k));
    }
    let scaled_conductance_j_k = finite(dt_s * conductance_w_k, "step-scaled interface conductance")?;
    let mut config = exchange::exchange_config(false);
    // Q_solid - Q_air = hA * (updated_reference - applied_reference).
    // Allocate at most one percent of the original joule allowance to this
    // mismatch, retaining the independent default watt gates without widening.
    let energy_reference_tolerance_k =
        (0.01 * step_config.energy_tolerance_j) / scaled_conductance_j_k;
    config.temperature_tolerance_k = config.temperature_tolerance_k.min(energy_reference_tolerance_k);
    if !(scaled_conductance_j_k > 0.0
        && config.temperature_tolerance_k.is_finite() && config.temperature_tolerance_k > 0.0)
    {
        return Err(bad("transient airflow cannot represent a positive energy-bounded reference tolerance"));
    }
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    if let Some(lowering) = radiation {
        lowering.config.validate().map_err(lower)?;
    }
    if let Some(policy) = nonlinear {
        policy.validate().map_err(lower)?;
    }
    let max_solid_solves = config.max_iterations.checked_mul(radiation_factor)
        .ok_or_else(|| bad("transient airflow/radiation trial allowance overflows"))?;
    let max_krylov = max_solid_solves.checked_mul(step_config.linear.max_iterations)
        .ok_or_else(|| bad("transient airflow Krylov allowance overflows"))?;
    let max_updates = max_solid_solves.checked_mul(nonlinear.map_or(0, |p| p.max_iterations))
        .ok_or_else(|| bad("transient airflow Newton allowance overflows"))?;
    let max_backtracks = max_updates.checked_mul(nonlinear.map_or(0, |p| p.line_search.max_backtracks))
        .ok_or_else(|| bad("transient airflow backtrack allowance overflows"))?;
    if let Some(policy) = nonlinear {
        policy.line_search.max_backtracks.checked_add(1)
            .and_then(|trials| trials.checked_mul(max_updates))
            .ok_or_else(|| bad("transient airflow Newton trial allowance overflows"))?;
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
        let replacements = regions.iter().map(|&(target, region, coefficient)| {
            let reference = *references.get(target)
                .ok_or_else(|| bad("transient airflow omitted an applied Robin reference"))?;
            Ok((region, coefficient, reference))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let boundary = problem.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
        let trial_problem = ConductionProblem { boundary: &boundary, ..problem };
        let (endpoint, nonlinear_evidence, radiation_evidence, radiative,
            physical_energy_residual_j, physical_dirichlet_in_w, radiation_out_w) =
            if let Some(lowering) = radiation {
                let solved = engine.advance_prescribed_with_ambient_radiation(
                    cx, trial_problem, interfaces, old, dt_s, step_config, nonlinear,
                    &lowering.patches, lowering.config,
                ).map_err(lower)?;
                count(&mut solid_solves, solved.radiation.iterations, max_solid_solves)?;
                count(&mut radiation_trials, solved.radiation.iterations, max_solid_solves)?;
                count(&mut krylov_iterations, solved.radiation.krylov_iterations, max_krylov)?;
                count(&mut nonlinear_updates, solved.radiation.solid_iterations, max_updates)?;
                count(&mut nonlinear_backtracks, solved.nonlinear_backtracks, max_backtracks)?;
                let nonlinear_evidence = solved.nonlinear.as_ref().map(nonlinear_row)
                    .transpose()?.unwrap_or_else(|| "null".to_string());
                let radiation_evidence = radiation_row(&solved)?;
                let radiation_out_w = solved.radiation.nonlinear_radiation_out_w;
                let radiative = RadiativeEndpoint {
                    combined_boundary: solved.combined_boundary,
                    convective_robin_fluxes: solved.convective_robin_fluxes,
                    convective_out_w: solved.convective_out_w,
                    report: solved.radiation,
                };
                (solved.conduction, nonlinear_evidence, radiation_evidence, Some(radiative),
                    solved.physical_energy_residual_j, solved.physical_dirichlet_in_w, radiation_out_w)
            } else if let Some(policy) = nonlinear {
                let solved = engine.advance_nonlinear_prescribed(
                    cx, trial_problem, interfaces, old, dt_s, step_config, policy,
                ).map_err(lower)?;
                count(&mut solid_solves, 1, max_solid_solves)?;
                count(&mut krylov_iterations, solved.step.krylov_iterations, max_krylov)?;
                count(&mut nonlinear_updates, solved.nonlinear_iterations, max_updates)?;
                count(&mut nonlinear_backtracks, solved.backtracks, max_backtracks)?;
                let row = nonlinear_row(&solved)?;
                let energy = solved.step.energy_residual_j;
                let reaction = solved.step.dirichlet_in_w;
                (solved.step, row, "null".to_string(), None, energy, reaction, 0.0)
            } else {
                let endpoint = engine.advance_prescribed(
                    cx, trial_problem, interfaces, old, dt_s, step_config,
                ).map_err(lower)?;
                count(&mut solid_solves, 1, max_solid_solves)?;
                count(&mut krylov_iterations, endpoint.krylov_iterations, max_krylov)?;
                let energy = endpoint.energy_residual_j;
                let reaction = endpoint.dirichlet_in_w;
                (endpoint, "null".to_string(), "null".to_string(), None, energy, reaction, 0.0)
            };
        poll(cx, deadline)?;
        maximum_physical_energy_residual_j =
            maximum_physical_energy_residual_j.max(physical_energy_residual_j.abs());
        let (fluxes, convective_out_w) = radiative.as_ref().map_or(
            (endpoint.robin_fluxes.as_slice(), endpoint.robin_out_w),
            |result| (result.convective_robin_fluxes.as_slice(), result.convective_out_w),
        );
        // Only original convection reaches the air marcher. Ambient radiation
        // exchanges with its own reservoir even on a shared exterior trace.
        let states = path.segments.iter().map(|segment| {
            fluxes.iter().find(|flux| flux.region == segment.target)
                .map(SolidRegionState::from_robin_flux)
                .ok_or_else(|| bad(format!(
                    "transient airflow target '{}' has no solved convective flux", segment.target,
                )))
        }).collect::<Result<Vec<_>, _>>()?;
        let off_path_convective_out_w = total(fluxes.iter()
            .filter(|flux| !targets.contains(flux.region.as_str()))
            .map(|flux| flux.heat_rate_w), "off-path convective heat")?;
        last = Some(SolidTrial {
            endpoint, nonlinear_evidence, radiation_evidence, radiation: radiative,
            applied_references: references.clone(), applied_boundary: boundary,
            convective_out_w, off_path_convective_out_w, radiation_out_w, physical_dirichlet_in_w,
        });
        Ok(states)
    })?;
    poll(cx, deadline)?;
    let last = last.ok_or_else(|| bad("transient airflow accepted no solid response"))?;
    let outcome = exchange::cross_check_decomposition(
        outcome, last.convective_out_w, last.off_path_convective_out_w,
    )?;
    let air_heat_gain_w = total(outcome.solution.branches.iter()
        .map(|branch| branch.march.total_heat_rate_w), "air enthalpy gain")?;
    let net_input_w = finite(
        last.endpoint.source_w + last.physical_dirichlet_in_w
            - last.endpoint.neumann_out_w - last.off_path_convective_out_w
            - air_heat_gain_w - last.radiation_out_w,
        "coupled endpoint net input",
    )?;
    let coupled_energy_residual_j = finite(
        last.endpoint.stored_energy_change_j - dt_s * net_input_w,
        "coupled endpoint energy residual",
    )?;
    if coupled_energy_residual_j.abs() > step_config.energy_tolerance_j {
        return Err(conduction_error(
            "cli-solve-conduction-transient-airflow-energy",
            format!(
                "transient storage minus actual air, radiation and boundary exchange is {coupled_energy_residual_j:.4e} J against the declared {:.4e} J gate",
                step_config.energy_tolerance_j,
            ),
            "tighten the thermal solve and inspect the interface heat balance; no endpoint is published",
        ));
    }
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
        air_heat_gain_w,
        off_path_convective_out_w: last.off_path_convective_out_w,
        radiation_out_w: last.radiation_out_w,
        physical_dirichlet_in_w: last.physical_dirichlet_in_w,
        coupled_energy_residual_j,
    })
}

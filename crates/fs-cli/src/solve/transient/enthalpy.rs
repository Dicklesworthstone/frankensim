//! Native fixed-mesh total-enthalpy trajectories through the shared P1 owner.
//!
//! Reference mass remains fixed. The accepted enthalpy, including latent heat,
//! is the only storage history; temperature is a constitutive observation.

use fs_conduction::transient::enthalpy::{
    EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError, ReferenceEnthalpyMaterial,
    },
    radiation::EnthalpyRadiationStepSolution,
};
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};

use super::*;

mod conductivity;
mod conjugate;
mod enclosure;
mod natural;
mod pointwise;
pub(in crate::solve::transient) use pointwise::evaluate as evaluate_pointwise_endpoint;
mod radiative;

const NO_CLAIM: &str = "Estimated stationary P1 total-enthalpy transport with frozen reference mass and caller-declared equilibrium phase charts. Latent heat is stored only in specific enthalpy; no apparent heat capacity or additional latent account is used. Existing sourced conductivity, finite-resistance contact and prescribed temperature, Neumann and Robin boundaries remain active. A prescribed temperature requires a unique enthalpy on every incident phase chart; a latent plateau cannot silently select a phase fraction. The initial physical enthalpy remains unchanged, and prescribed-temperature reactions include its complete storage jump and all physical boundary heat exactly once. Declared regional conductivity multipliers use each element's mean nodal liquid mass fraction and the shared exact enthalpy tangent, including phase feedback on latent plateaus. The base tensor retains its temperature dependence and source support; multipliers are declared phenomenological data, not a volume-fraction mixture or validated material claim. Declared ambient radiation defaults to sourced area-mean gray patches; every trial retains the same old enthalpy and the accepted endpoint passes independent physical residual and nonlinear-heat energy checks. Declared natural convection evaluates the existing card at each candidate wall temperature and requires the complete physical residual and energy balance to close; exact equilibrium and out-of-domain Rayleigh states retain their refusals. Declared quasi-steady airflow uses the existing branch heat-exchange owner and checks the complete endpoint at updated air references, with actual air enthalpy gain, off-path convection and nonlinear radiation counted separately. Natural reservoirs on distinct boundary targets compose with forced-air transport: every air trial encloses the existing natural/radiation endpoint solve while retaining the same accepted old enthalpy. The final physical audit evaluates updated natural coefficients and air references at the counted state; natural-reservoir and radiation heat stay separate from the forced-convective flux sent to air. Every workload switch is represented on both independent time grids. The final field comes only from the completed fine grid. Its forward sampled peak includes the original regional initial enthalpy temperature and every accepted fine endpoint in the same requirement region, retaining the earliest exact tie. This observation needs no adjoint and does not bound temperatures between samples. The temperature difference is an assumed-first-order temporal estimate, not a spatial or continuous-time maximum bound; agreement on a latent plateau does not prove convergence of enthalpy or phase fraction. Supplied material-card identities and source text are attribution, not verification of the declared curves. Requested final or sampled-peak history derivatives revalidate the accepted fine endpoints and carry the shared exact specific-enthalpy transpose through prescribed convection, contact, complete declared ambient-radiation feedback, and the full natural-convection and streamwise-air boundary derivatives. Reported controls are initial enthalpy, delivered regional power, fixed reference density, independently declared prescribed-bath temperatures, per-surface absolute emissivity and radiation reservoir temperature, natural-reservoir and air-inlet temperatures, and independently controlled absolute fan-bank speed ratios. Fan controls include the retained hydraulic, convection-card and downstream-air parameter chain; unsupported coordinates keep explicit reasons. Independently varied prescribed boundaries sharing vertices are unsupported; a peak selected at the initial state has zero future bath sensitivity. Constitutive calibration and time-grid derivatives are fixed. Radiative controls hold query temperature and selected claim support fixed, retain source identities and restrict emissivity increments to feasible directions. Geometry, reference density and energetic internal variables are fixed during each forward trajectory. Declared closed gray-diffuse enclosures compose with total-enthalpy storage and prescribed convection through the shared reflected-radiation owner; every spatial rung binds actual mesh areas and both time grids retain immutable old enthalpy across all trials. Fresh patch exchange is internal, and the complete physical residual and storage/heat balance must pass before history advances. Natural/airflow endpoint loops and requested history derivatives use the same fresh closed-exchange model and the existing complete boundary feedback; original convection alone heats transported air. Explicit pointwise radiation integrates local P1 fourth-power emission and its exact tangent inside the same enthalpy Newton solve. Closed enclosures retain patch-uniform reflected irradiation; original convection alone heats air. Natural and airflow endpoint audits use that same spatial radiation owner, including original prescribed flux and reservoir reactions. The pointwise path preserves source support and fixed reference mass. Its requested trajectory derivatives reuse the complete local-emission and reflected-irradiation transpose, including spatially resolved absolute emissivity and ambient-reservoir controls. No material motion, density-driven expansion, phase kinetics, airflow storage, time-varying fan dynamics is inferred.";

/// Storage-specific scope carried by the outer conduction report.
pub(super) fn no_claim() -> &'static str {
    NO_CLAIM
}

/// Check product combinations before the staged driver starts expensive work.
pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    let Some(policy) = policy(spec).filter(|time| time.enthalpy.is_some()) else {
        return Ok(());
    };
    if !policy.capacities.is_empty() {
        return Err(bad(
            "enthalpy and volumetric heat-capacity storage cannot both be active",
        ));
    }
    adjoint::requested(spec)?;
    Ok(())
}

fn failed(error: EnthalpyError) -> SolveRefusal {
    match error {
        EnthalpyError::Conduction(error @ ConductionError::Cancelled { .. }) => lower(error),
        other => bad(format!("total-enthalpy transport refused: {other}")),
    }
}

fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("enthalpy trajectory arithmetic is not finite"))
    }
}

fn numbers(
    cx: &Cx<'_>,
    values: &[f64],
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    let mut encoded = String::new();
    for (index, &value) in values.iter().enumerate() {
        if index % 256 == 0 {
            poll(cx, deadline)?;
        }
        if index != 0 {
            encoded.push(',');
        }
        encoded.push_str(&number(value)?);
    }
    Ok(encoded)
}

/// The legacy mean-patch ambient carrier is local to this enthalpy lane.
/// Capacity transients have their own enum and enclosure endpoint semantics.
struct RadiativeEndpoint {
    combined_boundary: fs_conduction::ThermalBoundary,
    convective_robin_fluxes: Vec<fs_conduction::RobinFlux>,
    convective_out_w: f64,
    report: fs_conduction::AmbientRadiationReport,
}

struct GridResult {
    endpoint: EnthalpyStepSolution,
    endpoint_dt_s: f64,
    threshold_j: f64,
    endpoint_residual_j: f64,
    endpoint_energy_residual_j: f64,
    endpoint_energy: PhysicalEnergy,
    rows: String,
    sampled_peak: SampledPeak,
    linear: Vec<LinearSolveEvidence>,
    residuals_w: Vec<f64>,
    updates: usize,
    stored_j: f64,
    net_input_j: f64,
    maximum_energy_residual_j: f64,
    radiation: Option<RadiativeEndpoint>,
    enclosure: Option<radiation::EnclosureEndpoint>,
    radiation_trials: usize,
    radiation_updates: usize,
    radiation_krylov_iterations: usize,
    maximum_radiation_energy_residual_j: f64,
    natural: Option<NaturalEndpoint>,
    natural_iterations: usize,
    natural_solid_solves: usize,
    natural_updates: usize,
    natural_krylov_iterations: usize,
    natural_endpoint_evaluations: usize,
    airflow: Option<AirflowEndpoint>,
    airflow_iterations: usize,
    airflow_solid_solves: usize,
    airflow_updates: usize,
    airflow_krylov_iterations: usize,
    airflow_endpoint_evaluations: usize,
    maximum_airflow_physical_energy_residual_j: f64,
}

fn airflow_step_row(
    path: &airflow::ConjugatePath,
    solved: &conjugate::AirflowStep,
) -> Result<String, SolveRefusal> {
    let references = |rows: &BTreeMap<String, f64>| -> Result<String, SolveRefusal> {
        rows.iter()
            .map(|(name, &value)| Ok(format!("{}:{}", json_string(name), number(value)?)))
            .collect::<Result<Vec<_>, SolveRefusal>>()
            .map(|rows| rows.join(","))
    };
    let fluxes = solved.physical_convective_robin_fluxes.iter().map(|row| {
        Ok(format!(
            "{{\"region\":{},\"area_m2\":{},\"coefficient_w_m2_k\":{},\"mean_wall_temperature_k\":{},\"reference_temperature_k\":{},\"heat_rate_w\":{}}}",
            json_string(&row.region), number(row.area_m2)?, number(row.mean_htc_w_per_m2_k)?,
            number(row.mean_wall_temperature_k)?, number(row.mean_reference_temperature_k)?,
            number(row.heat_rate_w)?,
        ))
    }).collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        "{{\"iterations\":{},\"solid_solves\":{},\"krylov_iterations\":{},\"nonlinear_iterations\":{},\"radiation_trials\":{},\"endpoint_evaluations\":{},\"reference_tolerance_k\":{},\"applied_references_k\":{{{}}},\"physical_references_k\":{{{}}},\"residual_reference\":{},\"initial_seed_residual_norm_j\":{},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"physical_convective_out_w\":{},\"radiative_out_w\":{},\"physical_dirichlet_in_w\":{},\"air_heat_gain_w\":{},\"off_path_convective_out_w\":{},\"coupled_energy_residual_j\":{},\"physical_convective_robin_fluxes\":[{}],\"exchange\":{}}}",
        solved.outcome.solution.iterations,
        solved.solid_solves,
        solved.krylov_iterations,
        solved.nonlinear_updates,
        solved.radiation_trials,
        solved.endpoint_evaluations,
        number(solved.reference_tolerance_k)?,
        references(&solved.applied_references)?,
        references(&solved.physical_references)?,
        json_string(solved.residual_reference),
        number(solved.initial_seed_residual_norm_j)?,
        number(solved.physical_residual_norm_j)?,
        number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?,
        number(solved.physical_convective_out_w)?,
        number(solved.physical_radiative_out_w)?,
        number(solved.physical_dirichlet_in_w)?,
        number(solved.air_heat_gain_w)?,
        number(solved.off_path_convective_out_w)?,
        number(solved.coupled_energy_residual_j)?,
        fluxes.join(","),
        airflow::receipt_fragment(path, &solved.outcome)?,
    ))
}

fn natural_step_row(solved: &natural::NaturalStep) -> Result<String, SolveRefusal> {
    let coefficients = |rows: &BTreeMap<String, f64>| -> Result<String, SolveRefusal> {
        rows.iter()
            .map(|(name, &value)| Ok(format!("{}:{}", json_string(name), number(value)?)))
            .collect::<Result<Vec<_>, SolveRefusal>>()
            .map(|rows| rows.join(","))
    };
    let fluxes = solved.physical_convective_robin_fluxes.iter().map(|row| {
        Ok(format!(
            "{{\"region\":{},\"area_m2\":{},\"coefficient_w_m2_k\":{},\"mean_wall_temperature_k\":{},\"reference_temperature_k\":{},\"heat_rate_w\":{}}}",
            json_string(&row.region), number(row.area_m2)?, number(row.mean_htc_w_per_m2_k)?,
            number(row.mean_wall_temperature_k)?, number(row.mean_reference_temperature_k)?,
            number(row.heat_rate_w)?,
        ))
    }).collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        "{{\"iterations\":{},\"solid_solves\":{},\"krylov_iterations\":{},\"nonlinear_iterations\":{},\"radiation_trials\":{},\"endpoint_evaluations\":{},\"applied_coefficients_w_m2_k\":{{{}}},\"physical_coefficients_w_m2_k\":{{{}}},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"physical_convective_out_w\":{},\"physical_radiative_out_w\":{},\"physical_dirichlet_in_w\":{},\"physical_convective_robin_fluxes\":[{}],\"exchange\":{}}}",
        solved.iterations,
        solved.solid_solves,
        solved.krylov_iterations,
        solved.nonlinear_updates,
        solved.radiation_trials,
        solved.endpoint_evaluations,
        coefficients(&solved.applied_coefficients)?,
        coefficients(&solved.physical_coefficients)?,
        number(solved.physical_residual_norm_j)?,
        number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?,
        number(solved.physical_convective_out_w)?,
        number(solved.physical_radiative_out_w)?,
        number(solved.physical_dirichlet_in_w)?,
        fluxes.join(","),
        solved.receipt,
    ))
}

fn radiation_step_row(solved: &EnthalpyRadiationStepSolution) -> Result<String, SolveRefusal> {
    let report = &solved.radiation;
    let patches = report.patches.iter().map(|row| {
        Ok(format!(
            "{{\"region\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"heat_mismatch_w\":{},\"heat_tolerance_w\":{}}}",
            json_string(row.patch.region()), number(row.area_m2)?,
            number(row.mean_surface_temperature_k)?, number(row.applied_heat_w)?,
            number(row.nonlinear_heat_w)?, number(row.heat_mismatch_w)?,
            number(row.heat_tolerance_w)?,
        ))
    }).collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        "{{\"iterations\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"max_temperature_change_k\":{},\"max_heat_mismatch_w\":{},\"decomposition_residual_w\":{},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"physical_dirichlet_in_w\":{},\"final_inner_residual_j\":{},\"patches\":[{}]}}",
        report.iterations,
        report.solid_iterations,
        report.krylov_iterations,
        number(solved.convective_out_w)?,
        number(report.applied_radiation_out_w)?,
        number(report.nonlinear_radiation_out_w)?,
        number(report.max_temperature_change_k)?,
        number(report.max_heat_mismatch_w)?,
        number(report.decomposition_residual_w)?,
        number(solved.physical_residual_norm_j)?,
        number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?,
        number(solved.physical_dirichlet_in_w)?,
        number(solved.conduction.newton.residual_norm)?,
        patches.join(","),
    ))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn march(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    initial: &[f64],
    initial_goal_k: f64,
    ends_s: &[f64],
    workload: Option<&Workload<'_>>,
    config: EnthalpyStepConfig,
    radiation: Option<&radiation::LoweredRadiation>,
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    airflow: Option<&airflow::ConjugatePath>,
    natural: &[buoyancy::NaturalLaw],
    pressure_pa: f64,
    labels: &[u32],
    goal_region: Option<u32>,
    mut tape: Option<&mut adjoint::EnthalpyTape>,
    deadline: Option<(Instant, f64)>,
) -> Result<GridResult, SolveRefusal> {
    let mut old = initial.to_vec();
    let mut sampled_peak = SampledPeak::new(initial_goal_k)?;
    let mut rows = String::new();
    let mut linear = Vec::new();
    let mut residuals_w = Vec::with_capacity(ends_s.len());
    let mut last = None;
    let mut time = 0.0;
    let mut endpoint_dt_s = 0.0;
    let mut threshold_j = 0.0;
    let mut endpoint_residual_j = 0.0;
    let mut endpoint_energy_residual_j = 0.0;
    let mut endpoint_energy = None;
    let mut updates = 0;
    let mut stored_j = 0.0;
    let mut net_input_j = 0.0;
    let mut maximum_energy_residual_j = 0.0_f64;
    let mut last_radiation = None;
    let mut last_enclosure = None;
    let mut radiation_trials = 0;
    let mut radiation_updates = 0;
    let mut radiation_krylov_iterations = 0;
    let mut maximum_radiation_energy_residual_j = 0.0_f64;
    let mut last_natural = None;
    let mut natural_iterations = 0;
    let mut natural_solid_solves = 0;
    let mut natural_updates = 0;
    let mut natural_krylov_iterations = 0;
    let mut natural_endpoint_evaluations = 0;
    let mut last_airflow = None;
    let mut airflow_iterations = 0;
    let mut airflow_solid_solves = 0;
    let mut airflow_updates = 0;
    let mut airflow_krylov_iterations = 0;
    let mut airflow_endpoint_evaluations = 0;
    let mut maximum_airflow_physical_energy_residual_j = 0.0_f64;
    for (index, &end) in ends_s.iter().enumerate() {
        poll(cx, deadline)?;
        let dt = end - time;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(bad(
                "every enthalpy time interval must have positive finite duration",
            ));
        }
        let source = workload
            .map(|history| history.source(cx, problem.mesh, labels, time, deadline))
            .transpose()?;
        let step_problem = ConductionProblem {
            source: source.as_ref().unwrap_or(problem.source),
            ..problem
        };
        let mut physical = None;
        let mut physical_neumann_out_w = None;
        let mut physical_dirichlet_in_w = None;
        let mut radiation_fragment = String::new();
        let mut natural_fragment = String::new();
        let mut airflow_fragment = String::new();
        let endpoint = if let Some(path) = airflow {
            let solved = conjugate::advance(
                cx, engine, step_problem, interfaces, &old, dt, config,
                radiation, enclosure, natural, pressure_pa, path, deadline,
            )?;
            physical = Some((
                finite(
                    solved.endpoint.source_w
                        + solved.physical_dirichlet_in_w
                        - solved.physical_neumann_out_w
                        - solved.off_path_convective_out_w
                        - solved.air_heat_gain_w
                        - solved.physical_radiative_out_w,
                )?,
                solved.physical_residual_norm_j,
                solved.physical_residual_tolerance_j,
                solved.coupled_energy_residual_j,
            ));
            physical_neumann_out_w = Some(solved.physical_neumann_out_w);
            physical_dirichlet_in_w = Some(solved.physical_dirichlet_in_w);
            airflow_fragment = format!(",\"conjugate\":{}", airflow_step_row(path, &solved)?);
            airflow_iterations += solved.outcome.solution.iterations;
            airflow_solid_solves += solved.solid_solves;
            airflow_updates += solved.nonlinear_updates;
            airflow_krylov_iterations += solved.krylov_iterations;
            airflow_endpoint_evaluations += solved.endpoint_evaluations;
            maximum_airflow_physical_energy_residual_j =
                maximum_airflow_physical_energy_residual_j
                    .max(solved.physical_energy_residual_j.abs());
            natural_iterations += solved.natural_iterations;
            natural_solid_solves += solved.natural_solid_solves;
            natural_updates += solved.natural_updates;
            natural_krylov_iterations += solved.natural_krylov_iterations;
            natural_endpoint_evaluations += solved.natural_endpoint_evaluations;
            if !natural.is_empty() {
                natural_fragment = format!(",\"natural\":{}", solved.natural_evidence);
            }
            last_natural = solved.natural;
            if radiation.is_some() {
                radiation_fragment = format!(",\"radiation\":{}", solved.radiation_evidence);
                radiation_trials += solved.radiation_trials;
                radiation_updates += solved.nonlinear_updates;
                radiation_krylov_iterations += solved.krylov_iterations;
                maximum_radiation_energy_residual_j = maximum_radiation_energy_residual_j
                    .max(solved.maximum_radiation_energy_residual_j);
            }
            last_radiation = solved.radiation;
            last_enclosure = solved.enclosure;
            last_airflow = Some(AirflowEndpoint {
                outcome: solved.outcome,
                applied_references: solved.applied_references,
                boundary: solved.applied_boundary,
            });
            solved.endpoint
        } else if !natural.is_empty() {
            let solved = natural::advance(
                cx,
                engine,
                step_problem,
                interfaces,
                &old,
                dt,
                config,
                radiation,
                enclosure,
                natural,
                pressure_pa,
                deadline,
            )?;
            physical = Some((
                finite(
                    solved.endpoint.source_w
                        + solved.physical_dirichlet_in_w
                        - solved.physical_neumann_out_w
                        - solved.physical_convective_out_w
                        - solved.physical_radiative_out_w,
                )?,
                solved.physical_residual_norm_j,
                solved.physical_residual_tolerance_j,
                solved.physical_energy_residual_j,
            ));
            physical_neumann_out_w = Some(solved.physical_neumann_out_w);
            physical_dirichlet_in_w = Some(solved.physical_dirichlet_in_w);
            natural_fragment = format!(",\"natural\":{}", natural_step_row(&solved)?);
            natural_iterations += solved.iterations;
            natural_solid_solves += solved.solid_solves;
            natural_updates += solved.nonlinear_updates;
            natural_krylov_iterations += solved.krylov_iterations;
            natural_endpoint_evaluations += solved.endpoint_evaluations;
            if radiation.is_some() {
                radiation_fragment = format!(",\"radiation\":{}", solved.radiation_evidence);
                radiation_trials += solved.radiation_trials;
                radiation_updates += solved.nonlinear_updates;
                radiation_krylov_iterations += solved.krylov_iterations;
                maximum_radiation_energy_residual_j = maximum_radiation_energy_residual_j
                    .max(solved.maximum_radiation_energy_residual_j);
            }
            last_radiation = solved.radiation;
            last_enclosure = solved.enclosure;
            last_natural = Some(NaturalEndpoint {
                receipt: solved.receipt,
                applied_coefficients: solved.applied_coefficients,
                boundary: solved.applied_boundary,
            });
            solved.endpoint
        } else if let Some(lowering) = radiation {
            let solved = radiative::advance(
                cx, engine, step_problem, interfaces, &old, dt, config,
                lowering, enclosure, deadline,
            )?;
            physical = Some((
                solved.net_input_w, solved.residual_norm_j,
                solved.residual_tolerance_j, solved.energy_residual_j,
            ));
            physical_neumann_out_w = Some(solved.neumann_out_w);
            physical_dirichlet_in_w = Some(solved.dirichlet_in_w);
            radiation_fragment = format!(",\"radiation\":{}", solved.evidence);
            radiation_trials += solved.solid_solves;
            radiation_updates += solved.updates;
            radiation_krylov_iterations += solved.krylov_iterations;
            maximum_radiation_energy_residual_j =
                maximum_radiation_energy_residual_j.max(solved.energy_residual_j.abs());
            last_radiation = solved.ambient;
            last_enclosure = solved.enclosure;
            solved.endpoint
        } else {
            engine
                .advance(cx, step_problem, interfaces, &old, dt, config)
                .map_err(failed)?
        };
        poll(cx, deadline)?;
        // The spatial assembler samples material at element means. Also check
        // every published nodal temperature against its incident source span.
        for (element, vertices) in problem.mesh.complex().tets.iter().enumerate() {
            if element % 256 == 0 {
                poll(cx, deadline)?;
            }
            let material = match problem.element_materials {
                Some(materials) => materials.model_for(element).map_err(lower)?,
                None => problem.material,
            };
            for &vertex in vertices {
                material
                    .temperature_span()
                    .check(endpoint.temperature[vertex as usize])
                    .map_err(lower)?;
            }
        }
        let initial_residual = endpoint
            .newton
            .history
            .first()
            .map_or(endpoint.newton.residual_norm, |row| row.residual_before);
        let inner_threshold_j = config
            .newton
            .absolute_tolerance
            .max(config.newton.relative_tolerance * initial_residual.max(f64::MIN_POSITIVE));
        let (net_w, residual_j, accepted_threshold_j, energy_residual_j) = match physical {
            Some(values) => values,
            None => (
                finite(endpoint.source_w + endpoint.dirichlet_in_w - endpoint.neumann_out_w - endpoint.robin_out_w)?,
                endpoint.newton.residual_norm,
                inner_threshold_j,
                endpoint.energy_residual_j,
            ),
        };
        threshold_j = accepted_threshold_j;
        endpoint_residual_j = residual_j;
        endpoint_energy_residual_j = energy_residual_j;
        let neumann_out_w = physical_neumann_out_w.unwrap_or(endpoint.neumann_out_w);
        let bath_w = physical_dirichlet_in_w.unwrap_or(endpoint.dirichlet_in_w);
        // Keep the fresh accepted heat law in the common report. The retained
        // endpoint still carries its last applied boundary for inner evidence.
        let robin_out_w = if physical.is_some() {
            finite(endpoint.source_w + bath_w - neumann_out_w - net_w)?
        } else {
            endpoint.robin_out_w
        };
        endpoint_energy = Some(PhysicalEnergy {
            neumann_out_w,
            robin_out_w,
            dirichlet_in_w: bath_w,
            net_input_w: net_w,
            energy_residual_j,
        });
        // Only final inner histories are retained by the boundary owners. The
        // conduction report counts precisely these retained Newton attempts;
        // all-trial work is reported separately in the boundary counters.
        let mut krylov_iterations = 0;
        for row in &endpoint.newton.history {
            residuals_w.push(finite(row.residual_before / dt)?);
            linear.push(LinearSolveEvidence {
                nonlinear_iteration: updates + row.iteration,
                method: if airflow.is_some() && !natural.is_empty() {
                    "fgmres-enthalpy-newton-final-mixed-convection-trial"
                } else if airflow.is_some() {
                    "fgmres-enthalpy-newton-final-airflow-trial"
                } else if !natural.is_empty() {
                    "fgmres-enthalpy-newton-final-natural-trial"
                } else if radiation.is_some() {
                    "fgmres-enthalpy-newton-final-radiation-trial"
                } else {
                    "fgmres-enthalpy-newton"
                },
                iterations: row.linear_iterations,
                reported: ResidualClaim::TrueEuclidean(row.linear_relative_residual),
                true_relative_residual: finite(row.linear_relative_residual)?,
                converged_true: row.linear_relative_residual <= row.forcing,
                stall: None,
            });
            krylov_iterations += row.linear_iterations;
        }
        updates += endpoint.newton.iterations;
        stored_j = finite(stored_j + endpoint.stored_energy_change_j)?;
        net_input_j = finite(net_input_j + dt * net_w)?;
        maximum_energy_residual_j = maximum_energy_residual_j.max(energy_residual_j.abs());
        let mut liquid_mass = 0.0;
        let mut total_mass = 0.0;
        for (vertex, (&mass, &fraction)) in engine
            .reference_nodal_masses_kg()
            .iter()
            .zip(&endpoint.liquid_mass_fraction)
            .enumerate()
        {
            if vertex % 256 == 0 {
                poll(cx, deadline)?;
            }
            liquid_mass = finite(liquid_mass + mass * fraction)?;
            total_mass = finite(total_mass + mass)?;
        }
        let qoi = ladder_functional(
            labels,
            &problem.mesh.complex().tets,
            &endpoint.temperature,
            goal_region,
        );
        sampled_peak.record(end, qoi)?;
        if index != 0 {
            rows.push(',');
        }
        let bath_fragment = if problem.boundary.dirichlet().is_empty() { String::new() } else {
            format!(",\"dirichlet_in_w\":{},\"physical_dirichlet_in_w\":{}",
                number(endpoint.dirichlet_in_w)?,
                number(physical_dirichlet_in_w.unwrap_or(endpoint.dirichlet_in_w))?)
        };
        rows.push_str(&format!(
            "{{\"step\":{},\"time_s\":{},\"dt_s\":{},\"final_region_max_k\":{},\"stored_energy_change_j\":{},\"net_input_w\":{},\"source_w\":{},\"energy_residual_j\":{},\"residual_j\":{},\"threshold_j\":{},\"newton_iterations\":{},\"krylov_iterations\":{krylov_iterations},\"liquid_mass_kg\":{},\"mean_liquid_mass_fraction\":{}{bath_fragment}{radiation_fragment}{natural_fragment}{airflow_fragment}}}",
            index + 1, number(end)?, number(dt)?, number(qoi)?,
            number(endpoint.stored_energy_change_j)?, number(net_w)?, number(endpoint.source_w)?,
            number(energy_residual_j)?, number(residual_j)?,
            number(threshold_j)?, endpoint.newton.iterations,
            number(liquid_mass)?, number(liquid_mass / total_mass)?,
        ));
        if let Some(history) = tape.as_deref_mut() {
            history.record(cx, &endpoint, threshold_j, deadline)?;
        }
        old.clone_from(&endpoint.specific_enthalpy_j_kg);
        last = Some(endpoint);
        endpoint_dt_s = dt;
        time = end;
    }
    Ok(GridResult {
        endpoint: last.ok_or_else(|| bad("enthalpy trajectory has no accepted time interval"))?,
        endpoint_dt_s,
        threshold_j,
        endpoint_residual_j,
        endpoint_energy_residual_j,
        endpoint_energy: endpoint_energy.ok_or_else(|| bad("enthalpy grid lost its final physical energy audit"))?,
        rows,
        sampled_peak,
        linear,
        residuals_w,
        updates,
        stored_j,
        net_input_j,
        maximum_energy_residual_j,
        radiation: last_radiation,
        enclosure: last_enclosure,
        radiation_trials,
        radiation_updates,
        radiation_krylov_iterations,
        maximum_radiation_energy_residual_j,
        natural: last_natural,
        natural_iterations,
        natural_solid_solves,
        natural_updates,
        natural_krylov_iterations,
        natural_endpoint_evaluations,
        airflow: last_airflow,
        airflow_iterations,
        airflow_solid_solves,
        airflow_updates,
        airflow_krylov_iterations,
        airflow_endpoint_evaluations,
        maximum_airflow_physical_energy_residual_j,
    })
}

/// Execute both time grids with the original physical problem and source data.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn solve(
    cx: &Cx<'_>,
    spec: &ProjectSpec,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    labels: &[u32],
    region_ids: &BTreeMap<String, u32>,
    radiation: Option<&radiation::LoweredRadiation>,
    airflow: Option<&airflow::ConjugatePath>,
    natural: &[buoyancy::NaturalLaw],
    pressure_pa: f64,
    linear: LinearConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<NativeTransientResult, SolveRefusal> {
    admit(spec)?;
    poll(cx, deadline)?;
    let policy = policy(spec).ok_or_else(|| bad("missing transient enthalpy declaration"))?;
    let storage = policy
        .enthalpy
        .as_ref()
        .ok_or_else(|| bad("missing enthalpy storage"))?;
    let objective = adjoint::requested(spec)?;
    let grids = TimeGrids::new(policy)?;
    let total_steps = grids
        .coarse
        .len()
        .checked_add(grids.fine.len())
        .ok_or_else(|| bad("enthalpy time work overflows"))?;
    let n = problem.mesh.vertex_count();
    let elements = problem.mesh.element_count();
    if n == 0 || elements == 0 || labels.len() != elements || linear.max_iterations == 0 {
        return Err(bad(
            "enthalpy needs a nonempty labeled mesh and a positive linear work cap",
        ));
    }
    let restart = n.min(24).min(linear.max_iterations);
    let mut minimum_step = 1.0;
    for _ in 0..storage.max_backtracks {
        minimum_step *= 0.5;
    }
    let newton = NewtonKrylovConfig {
        absolute_tolerance: 0.01 * policy.energy_tolerance.value / fs_math::det::sqrt(n as f64),
        relative_tolerance: linear.tolerance,
        linear_restart: restart,
        // Round down, so one Newton attempt never spends more than the
        // existing admitted linear-iteration allowance.
        max_linear_cycles: linear.max_iterations / restart,
        globalization: Globalization::LineSearch(LineSearchConfig {
            armijo: 1e-4,
            contraction: 0.5,
            minimum_step,
        }),
        ..NewtonKrylovConfig::default()
    };
    let config = EnthalpyStepConfig {
        newton,
        max_newton_iterations: storage.max_newton_iterations as usize,
        energy_tolerance_j: policy.energy_tolerance.value,
    };
    let natural_factor = if natural.is_empty() {
        1
    } else {
        buoyancy::MAX_ITERATIONS
    };
    let airflow_factor = airflow.map_or(1, |_| airflow::exchange_config(false).max_iterations);
    let radiation_factor = radiation.map_or(1, radiation::LoweredRadiation::max_solid_solves);
    let max_natural_iterations_per_step = natural_factor.checked_mul(airflow_factor)
        .ok_or_else(|| bad("combined enthalpy interface allowance overflows"))?;
    let max_solid_solves_per_step = max_natural_iterations_per_step
        .checked_mul(radiation_factor)
        .ok_or_else(|| bad("combined enthalpy boundary iteration allowance overflows"))?;
    let total_trial_limit = total_steps
        .checked_mul(max_solid_solves_per_step)
        .ok_or_else(|| bad("combined enthalpy radiation allowance overflows"))?;
    let max_newton_updates_per_step = max_solid_solves_per_step
        .checked_mul(config.max_newton_iterations)
        .ok_or_else(|| bad("per-step radiative enthalpy Newton allowance overflows"))?;
    let total_newton_limit = total_steps
        .checked_mul(max_newton_updates_per_step)
        .ok_or_else(|| bad("combined enthalpy Newton allowance overflows"))?;
    let max_krylov_iterations_per_step = max_newton_updates_per_step
        .checked_mul(restart)
        .and_then(|v| v.checked_mul(newton.max_linear_cycles))
        .ok_or_else(|| bad("per-step radiative enthalpy Krylov allowance overflows"))?;
    let total_krylov_limit = total_newton_limit
        .checked_mul(restart)
        .and_then(|v| v.checked_mul(newton.max_linear_cycles))
        .ok_or_else(|| bad("combined enthalpy Krylov allowance overflows"))?;
    let adjoint_krylov_limit = if objective.is_some() {
        grids.fine.len().checked_mul(linear.max_iterations)
            .and_then(|value| value.checked_mul(
                if radiation.is_some_and(|lowering| lowering.is_enclosure() && !lowering.is_pointwise()) { 2 } else { 1 }
            ))
            .ok_or_else(|| bad("enthalpy adjoint Krylov allowance overflows"))?
    } else {
        0
    };
    let total_forward_and_adjoint_krylov_limit = total_krylov_limit
        .checked_add(adjoint_krylov_limit)
        .ok_or_else(|| bad("combined enthalpy forward/adjoint work allowance overflows"))?;
    let total_backtrack_evaluation_limit = total_newton_limit
        .checked_mul(storage.max_backtracks as usize + 1)
        .ok_or_else(|| bad("combined enthalpy line-search allowance overflows"))?;
    let max_natural_endpoint_evaluations_per_step = if natural.is_empty() {
        0
    } else {
        max_natural_iterations_per_step.checked_mul(2)
            .ok_or_else(|| bad("enthalpy natural endpoint allowance overflows"))?
    };
    let total_natural_endpoint_evaluation_limit = total_steps
        .checked_mul(max_natural_endpoint_evaluations_per_step)
        .ok_or_else(|| bad("combined enthalpy natural endpoint allowance overflows"))?;
    let max_endpoint_evaluations_per_step = if airflow.is_some() {
        // Natural audits occur inside each air callback, alongside the seed,
        // per-callback old-state audit and final physical air audit.
        max_natural_endpoint_evaluations_per_step
            .checked_add(airflow_factor)
            .and_then(|value| value.checked_add(2))
            .ok_or_else(|| bad("enthalpy airflow endpoint allowance overflows"))?
    } else {
        max_natural_endpoint_evaluations_per_step
    };
    let total_endpoint_evaluation_limit = total_steps
        .checked_mul(max_endpoint_evaluations_per_step)
        .ok_or_else(|| bad("combined enthalpy endpoint allowance overflows"))?;
    let retained_newton_history_limit = total_steps.checked_mul(config.max_newton_iterations)
        .ok_or_else(|| bad("retained enthalpy Newton history allowance overflows"))?;
    let knot_count = storage.materials.iter().try_fold(0usize, |count, row| {
        count
            .checked_add(row.knots.len())
            .ok_or_else(|| bad("enthalpy knot count overflows"))
    })?;
    let radiation_bytes = if let Some(lowering) = radiation {
        n.checked_mul(256)
            .and_then(|v| elements.checked_mul(2048).and_then(|e| v.checked_add(e)))
            .and_then(|v| {
                lowering
                    .surface_count()
                    .checked_mul(total_steps)
                    .and_then(|p| p.checked_mul(1024))
                    .and_then(|p| v.checked_add(p))
            })
            .and_then(|v| {
                lowering.surface_count().checked_mul(lowering.surface_count())
                    .and_then(|entries| entries.checked_mul(128))
                    .and_then(|matrix| v.checked_add(matrix))
            })
            .ok_or_else(|| bad("enthalpy radiation boundary/receipt allowance overflows"))?
    } else {
        0
    };
    let natural_bytes = if natural.is_empty() {
        0
    } else {
        n.checked_mul(192)
            .and_then(|v| elements.checked_mul(1024).and_then(|e| v.checked_add(e)))
            .and_then(|v| {
                natural
                    .len()
                    .checked_mul(total_steps)
                    .and_then(|rows| rows.checked_mul(4096))
                    .and_then(|rows| v.checked_add(rows))
            })
            .ok_or_else(|| bad("enthalpy natural boundary/receipt allowance overflows"))?
    };
    let airflow_bytes = if let Some(path) = airflow {
        n.checked_mul(192)
            .and_then(|v| elements.checked_mul(1024).and_then(|e| v.checked_add(e)))
            .and_then(|v| {
                path.segments.len().checked_add(1)
                    .and_then(|rows| rows.checked_mul(total_steps))
                    .and_then(|rows| rows.checked_mul(4096))
                    .and_then(|rows| v.checked_add(rows))
            })
            .and_then(|v| {
                path.segments.len().checked_add(1)
                    .and_then(|rows| rows.checked_mul(airflow_factor))
                    .and_then(|rows| rows.checked_mul(512))
                    .and_then(|rows| v.checked_add(rows))
            })
            .ok_or_else(|| bad("enthalpy airflow boundary/receipt allowance overflows"))?
    } else {
        0
    };
    let adjoint_bytes = if objective.is_some() {
        adjoint::enthalpy_storage_bytes(
            n, elements, grids.fine.len(), storage.materials.len(),
            problem.boundary.conditions().len(), config.max_newton_iterations,
            linear.max_iterations,
        )?.checked_add(radiation.map(|lowering| {
            if lowering.is_pointwise() {
                adjoint::enthalpy_pointwise_storage_bytes(
                    n, problem.mesh.boundary().len(), lowering.surface_count(), grids.fine.len(),
                )
            } else {
                adjoint::enthalpy_radiation_storage_bytes(n, elements, lowering.surface_count())
            }
        }).transpose()?.unwrap_or(0))
            .and_then(|value| adjoint::enthalpy_convection_storage_bytes(
                n, elements, natural.len().checked_add(airflow.map_or(0, |path| path.segments.len()))?,
            ).ok().and_then(|extra| value.checked_add(extra)))
            .and_then(|value| adjoint::enthalpy_control_storage_bytes(
                spec, natural.len().checked_add(airflow.map_or(0, |path| path.segments.len()))?,
            ).ok().and_then(|extra| value.checked_add(extra)))
            .and_then(|value| adjoint::enthalpy_material_control_storage_bytes(
                n, elements, spec.materials.as_deref().unwrap_or(&[]).len(),
                spec.interface_cards.as_deref().unwrap_or(&[]).len(),
            ).ok().and_then(|extra| value.checked_add(extra)))
            .and_then(|value| adjoint::heat_flux_storage_bytes(spec, n)
                .ok().and_then(|extra| value.checked_add(extra)))
            .ok_or_else(|| bad("combined enthalpy forward/adjoint storage allowance overflows"))?
    } else {
        0
    };
    // Only the final inner Newton history of each accepted step is retained;
    // discarded air, natural and radiation trials reuse bounded workspace.
    // This is declared-memory admission, not allocator overhead or peak RSS.
    let bytes = n
        .checked_mul(8 * (2 * restart + 48) + 160)
        .and_then(|v| elements.checked_mul(4096).and_then(|e| v.checked_add(e)))
        .and_then(|v| knot_count.checked_mul(256).and_then(|k| v.checked_add(k)))
        .and_then(|v| total_steps.checked_mul(2048).and_then(|r| v.checked_add(r)))
        .and_then(|v| {
            retained_newton_history_limit
                .checked_mul(160)
                .and_then(|r| v.checked_add(r))
        })
        .and_then(|v| v.checked_add(radiation_bytes))
        .and_then(|v| v.checked_add(natural_bytes))
        .and_then(|v| v.checked_add(airflow_bytes))
        .and_then(|v| v.checked_add(adjoint_bytes))
        .ok_or_else(|| bad("enthalpy memory allowance overflows"))?;
    if u64::try_from(bytes).unwrap_or(u64::MAX)
        > spec
            .budgets
            .as_ref()
            .map_or(0, |budget| budget.memory_bytes)
    {
        return Err(bad(format!(
            "enthalpy numerical buffers and retained fields need an estimated {bytes} bytes, above the declared memory budget"
        )));
    }
    let mut curves = Vec::with_capacity(storage.materials.len());
    let mut by_region = BTreeMap::new();
    let mut material_rows = Vec::with_capacity(storage.materials.len());
    let mut initial_temperatures = Vec::with_capacity(storage.materials.len());
    for (index, row) in storage.materials.iter().enumerate() {
        poll(cx, deadline)?;
        let region = region_ids
            .get(&row.region)
            .ok_or_else(|| bad("enthalpy names an unknown mesh region"))?;
        if by_region.insert(*region, index).is_some() {
            return Err(bad("a mesh region has multiple enthalpy assignments"));
        }
        let curve = row
            .phase_curve()
            .map_err(|v| conduction_error(v.code, v.what, v.fix))?;
        let state = curve.state_at_specific_enthalpy(row.initial_specific_enthalpy.value)
            .map_err(|error| bad(format!("initial enthalpy for '{}' refuses: {error}", row.region)))?;
        let initial_temperature = state.temperature_k();
        if !(initial_temperature.is_finite() && initial_temperature > 0.0) {
            return Err(bad("regional initial enthalpy must resolve to positive absolute temperature"));
        }
        initial_temperatures.push(initial_temperature);
        let phase_conductivity = conductivity::receipt(row)?;
        material_rows.push(format!(
            "{{\"region\":{},\"material_index\":{index},\"material_card_identity\":{},\"phase_curve_identity\":{},\"source\":{},\"reference_density_kg_m3\":{},\"initial_specific_enthalpy_j_kg\":{},\"initial_temperature_k\":{}{phase_conductivity}}}",
            json_string(&row.region), json_string(&row.material_card_identity.to_hex()),
            json_string(&curve.identity().to_hex()), json_string(&row.source),
            number(row.reference_density.value)?, number(row.initial_specific_enthalpy.value)?,
            number(initial_temperature)?,
        ));
        curves.push(curve);
    }
    let materials = storage
        .materials
        .iter()
        .zip(&curves)
        .map(|(row, curve)| ReferenceEnthalpyMaterial {
            curve,
            reference_density_kg_m3: row.reference_density.value,
        })
        .collect::<Vec<_>>();
    let mut assignment = Vec::with_capacity(elements);
    for (element, label) in labels.iter().enumerate() {
        if element % 256 == 0 {
            poll(cx, deadline)?;
        }
        assignment.push(
            *by_region
                .get(label)
                .ok_or_else(|| bad("a retained element has no declared enthalpy chart"))?,
        );
    }
    // Validate the actual regional initial state against the sourced base
    // conductivity before any time integration or boundary iteration.
    for (element, &material_id) in assignment.iter().enumerate() {
        if element % 256 == 0 {
            poll(cx, deadline)?;
        }
        let material = match problem.element_materials {
            Some(materials) => materials.model_for(element).map_err(lower)?,
            None => problem.material,
        };
        material.temperature_span()
            .check(initial_temperatures[material_id])
            .map_err(lower)?;
    }
    let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, problem.mesh, &materials, &assignment,
        EnthalpyBudget { max_vertices: n, max_elements: elements })
        .map_err(|error| match error {
            HeterogeneousEnthalpyError::Enthalpy(inner) => failed(inner),
            other => bad(format!("enthalpy material assignment refused: {other}; distinct regional charts require distinct interface vertices and an explicit contact law")),
        })?;
    let engine = conductivity::bind(cx, engine, storage, &assignment, deadline)?;
    let mut initial = Vec::with_capacity(n);
    for (vertex, &id) in engine.vertex_material_ids().iter().enumerate() {
        if vertex % 256 == 0 {
            poll(cx, deadline)?;
        }
        initial.push(storage.materials[id].initial_specific_enthalpy.value);
    }
    let enclosure = radiation
        .map(|lowering| lowering.bind_enclosure(cx, problem))
        .transpose()
        .map_err(lower)?
        .flatten();
    let goal_region = ladder_target_region(spec, region_ids);
    if goal_region.is_none() {
        return Err(bad(
            "enthalpy temporal comparison needs a temperature-maximum requirement region",
        ));
    }
    let goal_region_name = spec.requirements.iter().flatten()
        .find(|row| region_ids.get(&row.region).copied() == goal_region)
        .map(|row| row.region.as_str())
        .ok_or_else(|| bad("transient trajectory lost its declared requirement region"))?;
    // A region owns one admitted phase chart and one original initial h.
    // Use that chart-resolved temperature, never a global uniform temperature
    // or the later prescribed endpoint, for sample zero.
    let initial_goal_k = goal_region.and_then(|id| by_region.get(&id))
        .and_then(|&index| initial_temperatures.get(index)).copied()
        .ok_or_else(|| bad("initial enthalpy peak has no selected regional chart"))?;
    let workload = Workload::bind(cx, spec, policy, problem.mesh, labels, region_ids, deadline)?;
    let coarse = march(
        cx,
        &engine,
        problem,
        interfaces,
        &initial,
        initial_goal_k,
        &grids.coarse,
        workload.as_ref(),
        config,
        radiation,
        enclosure.as_ref(),
        airflow,
        natural,
        pressure_pa,
        labels,
        goal_region,
        None,
        deadline,
    )?;
    let mut tape = if objective.is_some() {
        let mut temperatures = Vec::with_capacity(n);
        for (vertex, &id) in engine.vertex_material_ids().iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cx, deadline)?;
            }
            temperatures.push(initial_temperatures[id]);
        }
        Some(adjoint::EnthalpyTape::new(
            cx, &initial, &temperatures, grids.fine.len(), deadline,
        )?)
    } else {
        None
    };
    let fine = march(
        cx,
        &engine,
        problem,
        interfaces,
        &initial,
        initial_goal_k,
        &grids.fine,
        workload.as_ref(),
        config,
        radiation,
        enclosure.as_ref(),
        airflow,
        natural,
        pressure_pa,
        labels,
        goal_region,
        tape.as_mut(),
        deadline,
    )?;
    let adjoint_fragment = match (objective, tape) {
        (Some(objective), Some(tape)) => {
            let derivative = adjoint::finish_enthalpy(
                cx, spec, &engine, problem, interfaces, &grids.fine, workload.as_ref(),
                config, linear, labels, region_ids,
                goal_region.ok_or_else(|| bad("enthalpy adjoint lost its goal region"))?,
                storage, radiation, enclosure.as_ref(), airflow, natural, pressure_pa, tape, objective, deadline,
            )?;
            format!(",\"adjoint\":{derivative}")
        }
        (None, None) => String::new(),
        _ => return Err(bad("enthalpy adjoint request and retained history disagree")),
    };
    let coarse_qoi = ladder_functional(
        labels,
        &problem.mesh.complex().tets,
        &coarse.endpoint.temperature,
        goal_region,
    );
    let fine_qoi = ladder_functional(
        labels,
        &problem.mesh.complex().tets,
        &fine.endpoint.temperature,
        goal_region,
    );
    let difference = finite((fine_qoi - coarse_qoi).abs())?;
    let temporal_half_width_k = if difference > 0.0 {
        Some(finite(1.25 * difference)?)
    } else {
        None
    };
    let estimate = temporal_half_width_k
        .map(number)
        .transpose()?
        .unwrap_or_else(|| "null".into());
    let mut maximum_enthalpy_difference = 0.0_f64;
    let mut maximum_fraction_difference = 0.0_f64;
    for vertex in 0..n {
        if vertex % 256 == 0 {
            poll(cx, deadline)?;
        }
        maximum_enthalpy_difference = maximum_enthalpy_difference.max(finite(
            (fine.endpoint.specific_enthalpy_j_kg[vertex]
                - coarse.endpoint.specific_enthalpy_j_kg[vertex])
                .abs(),
        )?);
        maximum_fraction_difference = maximum_fraction_difference.max(
            (fine.endpoint.liquid_mass_fraction[vertex]
                - coarse.endpoint.liquid_mass_fraction[vertex])
                .abs(),
        );
    }
    let status = if temporal_half_width_k.is_some() {
        "estimated-order-one"
    } else {
        "no-resolved-time-difference"
    };
    let workload_receipt = workload
        .as_ref()
        .map_or("null", |history| history.receipt.as_str());
    let nonlinear_history_fragment = if airflow.is_some() {
        ",\"history_scope\":\"final-inner-newton-solve-of-each-accepted-airflow-endpoint\""
    } else if !natural.is_empty() {
        ",\"history_scope\":\"final-inner-newton-solve-of-each-accepted-natural-endpoint\""
    } else if radiation.is_some() {
        ",\"history_scope\":\"final-inner-newton-solve-of-each-accepted-time-step\""
    } else {
        ""
    };
    let radiation_receipt = if let Some(lowering) = radiation {
        let method = match (lowering.is_pointwise(), lowering.is_enclosure()) {
            (true, true) => "implicit-pointwise-p1-closed-gray-diffuse-enclosure",
            (true, false) => "implicit-pointwise-p1-gray-radiation",
            (false, true) => "implicit-closed-gray-diffuse-enclosure",
            (false, false) => "implicit-area-mean-gray-radiation",
        };
        format!(
            ",\"radiation\":{{\"method\":\"{method}\",\"max_trials_per_step\":{max_solid_solves_per_step},\"max_trials_per_natural_iteration\":{},\"max_trials_per_boundary_iteration\":{},\"max_newton_updates_per_trial\":{},\"max_newton_updates_per_step\":{max_newton_updates_per_step},\"max_krylov_iterations_per_newton_update\":{},\"max_krylov_iterations_per_step\":{max_krylov_iterations_per_step},\"total_trial_limit\":{total_trial_limit},\"total_backtrack_evaluation_limit\":{total_backtrack_evaluation_limit},\"coarse_trials\":{},\"fine_trials\":{},\"total_trials\":{},\"coarse_newton_updates\":{},\"fine_newton_updates\":{},\"total_newton_updates\":{},\"coarse_krylov_iterations\":{},\"fine_krylov_iterations\":{},\"total_krylov_iterations\":{},\"maximum_physical_energy_residual_j\":{},\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-trial\",\"nonlinear_detail\":\"step iteration fields, linear evidence and residual histories retain only the final inner solve; radiation counters include every trial across airflow and natural-convection iterations; step and final residuals and energy use the complete physical endpoint\"}}",
            lowering.max_solid_solves(),
            lowering.max_solid_solves(),
            config.max_newton_iterations,
            restart * newton.max_linear_cycles,
            coarse.radiation_trials,
            fine.radiation_trials,
            coarse.radiation_trials + fine.radiation_trials,
            coarse.radiation_updates,
            fine.radiation_updates,
            coarse.radiation_updates + fine.radiation_updates,
            coarse.radiation_krylov_iterations,
            fine.radiation_krylov_iterations,
            coarse.radiation_krylov_iterations + fine.radiation_krylov_iterations,
            number(
                coarse
                    .maximum_radiation_energy_residual_j
                    .max(fine.maximum_radiation_energy_residual_j)
            )?,
        )
    } else {
        String::new()
    };
    let natural_receipt_fragment = if natural.is_empty() {
        String::new()
    } else {
        format!(
            ",\"natural\":{{\"method\":\"implicit-area-mean-natural-convection\",\"max_iterations_per_step\":{max_natural_iterations_per_step},\"max_iterations_per_air_trial\":{natural_factor},\"max_solid_solves_per_step\":{max_solid_solves_per_step},\"max_newton_updates_per_step\":{max_newton_updates_per_step},\"max_krylov_iterations_per_step\":{max_krylov_iterations_per_step},\"max_endpoint_evaluations_per_step\":{max_natural_endpoint_evaluations_per_step},\"total_endpoint_evaluation_limit\":{total_natural_endpoint_evaluation_limit},\"coarse_iterations\":{},\"fine_iterations\":{},\"total_iterations\":{},\"coarse_solid_solves\":{},\"fine_solid_solves\":{},\"total_solid_solves\":{},\"coarse_newton_updates\":{},\"fine_newton_updates\":{},\"total_newton_updates\":{},\"coarse_krylov_iterations\":{},\"fine_krylov_iterations\":{},\"total_krylov_iterations\":{},\"total_endpoint_evaluations\":{},\"maximum_physical_energy_residual_j\":{},\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-air-natural-and-radiation-trial\",\"endpoint\":\"last-counted-solid-response-with-independent-full-physical-residual-and-energy-gates\",\"residual_reference\":\"first-physical-candidate-frozen-per-natural-solve\",\"equilibrium\":\"existing-card-domain-no-zero-delta-limiting-law\"}}",
            coarse.natural_iterations,
            fine.natural_iterations,
            coarse.natural_iterations + fine.natural_iterations,
            coarse.natural_solid_solves,
            fine.natural_solid_solves,
            coarse.natural_solid_solves + fine.natural_solid_solves,
            coarse.natural_updates,
            fine.natural_updates,
            coarse.natural_updates + fine.natural_updates,
            coarse.natural_krylov_iterations,
            fine.natural_krylov_iterations,
            coarse.natural_krylov_iterations + fine.natural_krylov_iterations,
            coarse.natural_endpoint_evaluations + fine.natural_endpoint_evaluations,
            number(if airflow.is_some() {
                coarse.maximum_airflow_physical_energy_residual_j
                    .max(fine.maximum_airflow_physical_energy_residual_j)
            } else {
                coarse.maximum_energy_residual_j.max(fine.maximum_energy_residual_j)
            })?,
        )
    };
    let airflow_receipt_fragment = if airflow.is_some() {
        format!(
            ",\"conjugate\":{{\"method\":\"implicit-total-enthalpy-quasi-steady-airflow\",\"max_iterations_per_step\":{airflow_factor},\"max_solid_solves_per_step\":{max_solid_solves_per_step},\"max_newton_updates_per_step\":{max_newton_updates_per_step},\"max_krylov_iterations_per_step\":{max_krylov_iterations_per_step},\"max_endpoint_evaluations_per_step\":{max_endpoint_evaluations_per_step},\"total_endpoint_evaluation_limit\":{total_endpoint_evaluation_limit},\"coarse_iterations\":{},\"fine_iterations\":{},\"total_iterations\":{},\"coarse_solid_solves\":{},\"fine_solid_solves\":{},\"total_solid_solves\":{},\"coarse_newton_updates\":{},\"fine_newton_updates\":{},\"total_newton_updates\":{},\"coarse_krylov_iterations\":{},\"fine_krylov_iterations\":{},\"total_krylov_iterations\":{},\"total_endpoint_evaluations\":{},\"maximum_physical_energy_residual_j\":{},\"maximum_coupled_energy_residual_j\":{},\"coupled_window_residual_j\":{},\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-airflow-natural-and-radiation-trial\",\"endpoint\":\"last-counted-solid-response-with-full-updated-air-reference-and-natural-coefficient-physical-residual-and-energy-gates\",\"residual_reference\":\"initial-enthalpy-with-branch-inlet-seed-boundary\",\"energy\":\"actual-air-enthalpy-gain-plus-off-path-convection-plus-nonlinear-radiation-counted-once\",\"nested_evidence\":\"step natural and radiation rows retain the final air trial; trajectory counters include every air trial\"}}",
            coarse.airflow_iterations,
            fine.airflow_iterations,
            coarse.airflow_iterations + fine.airflow_iterations,
            coarse.airflow_solid_solves,
            fine.airflow_solid_solves,
            coarse.airflow_solid_solves + fine.airflow_solid_solves,
            coarse.airflow_updates,
            fine.airflow_updates,
            coarse.airflow_updates + fine.airflow_updates,
            coarse.airflow_krylov_iterations,
            fine.airflow_krylov_iterations,
            coarse.airflow_krylov_iterations + fine.airflow_krylov_iterations,
            coarse.airflow_endpoint_evaluations + fine.airflow_endpoint_evaluations,
            number(coarse.maximum_airflow_physical_energy_residual_j
                .max(fine.maximum_airflow_physical_energy_residual_j))?,
            number(coarse.maximum_energy_residual_j.max(fine.maximum_energy_residual_j))?,
            number(fine.stored_j - fine.net_input_j)?,
        )
    } else {
        String::new()
    };
    let uniform_initial = policy.initial_temperature.uniform_temperature();
    let initial_temperature_receipt = uniform_initial
        .map(|quantity| number(quantity.value)).transpose()?
        .unwrap_or_else(|| "null".to_string());
    let initial_state = json_string(if uniform_initial.is_some() {
        "uniform-temperature"
    } else {
        "regional-enthalpy"
    });
    let sampled_peak_receipt = fine.sampled_peak.receipt(
        goal_region_name, grids.fine.len(), policy.horizon.value,
    )?;
    let receipt = format!(
        "{{\"schema\":\"fs-cli-transient-enthalpy-v1\",\"method\":\"backward-euler-total-enthalpy\",\"status\":\"completed\",\"qoi_time\":\"final\",\"initial_state\":{initial_state},\"initial_temperature_k\":{initial_temperature_receipt},\"final_time_s\":{},\"coarse_steps\":{},\"fine_steps\":{},\"total_steps\":{total_steps},\"max_steps\":{},\"coarse_max_step_s\":{},\"fine_max_step_s\":{},\"workload\":{workload_receipt},\"energy_tolerance_j\":{},\"materials\":[{}],\"nonlinear\":{{\"method\":\"newton-fgmres\",\"residual_units\":\"J\",\"report_residual_units\":\"W\",\"max_updates_per_step\":{},\"max_backtracks_per_update\":{},\"linear_restart\":{restart},\"max_linear_cycles\":{},\"total_newton_limit\":{total_newton_limit},\"retained_history_update_limit\":{retained_newton_history_limit},\"total_krylov_limit\":{total_krylov_limit},\"total_forward_and_adjoint_krylov_limit\":{total_forward_and_adjoint_krylov_limit},\"total_updates\":{},\"residual_atol_j\":{},\"residual_rtol\":{},\"armijo_c\":0.0001,\"line_search_shrink\":0.5{nonlinear_history_fragment}}}{radiation_receipt}{natural_receipt_fragment}{airflow_receipt_fragment}{adjoint_fragment},\"coarse\":[{}],\"fine\":[{}],\"sampled_peak\":{sampled_peak_receipt},\"temporal_error\":{{\"status\":{},\"qoi\":\"temperature-max\",\"coarse_final_k\":{},\"fine_final_k\":{},\"absolute_difference_k\":{},\"assumed_order\":1,\"safety_factor\":1.25,\"estimated_half_width_k\":{estimate},\"max_specific_enthalpy_difference_j_kg\":{},\"max_liquid_mass_fraction_difference\":{},\"spatial_error_measured\":false}},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{},\"maximum_step_residual_j\":{}}},\"endpoint\":{{\"vertex_order\":\"retained-conduction-mesh\",\"specific_enthalpy_j_kg\":[{}],\"liquid_mass_fraction\":[{}],\"reference_nodal_mass_kg\":[{}]}},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        number(policy.horizon.value)?,
        grids.coarse.len(),
        grids.fine.len(),
        policy.max_steps,
        number(TimeGrids::largest_step(&grids.coarse))?,
        number(TimeGrids::largest_step(&grids.fine))?,
        number(policy.energy_tolerance.value)?,
        material_rows.join(","),
        max_newton_updates_per_step,
        storage.max_backtracks,
        newton.max_linear_cycles,
        coarse.updates + fine.updates,
        number(newton.absolute_tolerance)?,
        number(newton.relative_tolerance)?,
        coarse.rows,
        fine.rows,
        json_string(status),
        number(coarse_qoi)?,
        number(fine_qoi)?,
        number(difference)?,
        number(maximum_enthalpy_difference)?,
        number(maximum_fraction_difference)?,
        number(fine.stored_j)?,
        number(fine.net_input_j)?,
        number(fine.stored_j - fine.net_input_j)?,
        number(fine.maximum_energy_residual_j)?,
        numbers(cx, &fine.endpoint.specific_enthalpy_j_kg, deadline)?,
        numbers(cx, &fine.endpoint.liquid_mass_fraction, deadline)?,
        numbers(cx, engine.reference_nodal_masses_kg(), deadline)?,
        json_string(NO_CLAIM),
    );
    let (conjugate_receipt, mut derived_boundary) = match (airflow, fine.airflow.as_ref()) {
        (Some(path), Some(accepted)) => {
            let derived = path.segments.iter().map(|segment| {
                let reference = accepted.applied_references.get(&segment.target)
                    .copied().ok_or_else(|| bad(
                        "the final enthalpy airflow endpoint lost an applied Robin reference",
                    ))?;
                Ok((segment.target.clone(), (segment.htc_w_m2_k, reference)))
            }).collect::<Result<BTreeMap<_, _>, SolveRefusal>>()?;
            (Some(airflow::receipt_fragment(path, &accepted.outcome)?), derived)
        }
        (None, None) => (None, BTreeMap::new()),
        _ => return Err(bad("the final enthalpy endpoint has inconsistent airflow evidence")),
    };
    let natural_receipt =
        match (natural.is_empty(), fine.natural.as_ref()) {
            (true, None) => None,
            (false, Some(accepted)) => {
                for law in natural {
                    let coefficient = *accepted.applied_coefficients.get(&law.target).ok_or_else(
                        || bad("the final natural enthalpy endpoint lost an applied coefficient"),
                    )?;
                    derived_boundary.insert(law.target.clone(), (coefficient, law.ambient_k));
                }
                Some(accepted.receipt.clone())
            }
            _ => {
                return Err(bad(
                    "the final enthalpy endpoint has inconsistent natural evidence",
                ));
            }
        };
    let endpoint = fine.endpoint;
    let energy = fine.endpoint_energy;
    let net_w = energy.net_input_w;
    let endpoint_storage_w = finite(endpoint.stored_energy_change_j / fine.endpoint_dt_s)?;
    let endpoint_energy_residual_w = finite(fine.endpoint_energy_residual_j / fine.endpoint_dt_s)?;
    let (provenance, receipts, assignment) = match problem.element_materials {
        Some(materials) => (
            materials.provenance(),
            materials.receipts().len(),
            Some(materials.identity()),
        ),
        None => (
            problem.material.provenance(),
            problem.material.receipts().len(),
            None,
        ),
    };
    let mut evidence = coarse.linear;
    for mut row in fine.linear {
        row.nonlinear_iteration += coarse.updates;
        evidence.push(row);
    }
    let mut residual_history = coarse.residuals_w;
    residual_history.extend(fine.residuals_w);
    let report = ConductionReport {
        iterations: coarse.updates + fine.updates,
        residual_history,
        final_residual: finite(fine.endpoint_residual_j / fine.endpoint_dt_s)?,
        residual_threshold: finite(fine.threshold_j / fine.endpoint_dt_s)?,
        stop_reason: StopReason::ResidualTolerance,
        linear: evidence,
        energy: EnergyBalance {
            source_w: endpoint.source_w,
            neumann_out_w: energy.neumann_out_w,
            robin_out_w: energy.robin_out_w,
            dirichlet_in_w: energy.dirichlet_in_w,
            closure_w: net_w,
            scale_w: endpoint
                .source_w
                .abs()
                .max(energy.dirichlet_in_w.abs())
                .max(energy.neumann_out_w.abs())
                .max(energy.robin_out_w.abs())
                .max(f64::MIN_POSITIVE),
        },
        material_provenance: provenance,
        material_receipts: receipts,
        element_material_identity: assignment,
        interface_fluxes: endpoint.contact_fluxes,
        robin_fluxes: endpoint.robin_fluxes,
        free_dofs: n - problem.boundary.dirichlet().len(),
        elements,
    };
    poll(cx, deadline)?;
    let conduction = ConductionSolution {
        temperature: endpoint.temperature,
        report,
    };
    let solution = if radiation.is_some_and(radiation::LoweredRadiation::is_pointwise) {
        if fine.enclosure.is_some() || fine.radiation.is_some() {
            return Err(bad("pointwise enthalpy endpoint retained a mean-patch radiation carrier"));
        }
        radiation::SolidSolution::from_pointwise_endpoint(
            conduction, endpoint.pointwise_radiation, endpoint.pointwise_enclosure_radiation,
        )?
    } else { match (enclosure, fine.enclosure, fine.radiation) {
        (Some(model), Some(endpoint), None) => {
            radiation::SolidSolution::from_enclosure_endpoint(conduction, model, endpoint)
        }
        (None, None, Some(endpoint)) => radiation::SolidSolution::from_radiative_endpoint(
            conduction,
            endpoint.combined_boundary,
            endpoint.convective_robin_fluxes,
            endpoint.convective_out_w,
            endpoint.report,
        ),
        (None, None, None) => radiation::SolidSolution::from_conduction(conduction),
        _ => return Err(bad("the final enthalpy endpoint has inconsistent radiation evidence")),
    }};
    Ok(NativeTransientResult {
        solution,
        receipt,
        conjugate_receipt,
        natural_receipt,
        derived_boundary,
        temporal_half_width_k,
        endpoint_storage_w,
        endpoint_energy_residual_w,
    })
}


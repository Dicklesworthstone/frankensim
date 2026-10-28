//! Native `.fsim` finite-time cooling through the shared backward-Euler kernel.
//!
//! Each grid advances immutable accepted history. The nested half grid starts
//! from the SAME declared initial field, and only its final field is published.
//! Coarse/fine work shares one explicit cap. Static regional capacities,
//! endpoint-evaluated conductivity, prescribed exterior laws and matching
//! contact are supported. Temperature-dependent conductivity uses the shared
//! Newton tangent, including k'(T), with immutable accepted physical history.
//! The retained final-QoI difference is Estimated temporal error; it supplies
//! neither a spatial bound nor a continuous-time peak claim.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_conduction::assemble::DofMap;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, NonlinearStepSolution, RadiationStepSolution, StepConfig,
    StepSolution,
};
use fs_conduction::{
    ConductionError, ConductionProblem, ConductionReport, ConductionSolution, EnergyBalance,
    LinearConfig, LinearSolveEvidence, ResidualClaim, StopReason, ThermalInterfaces,
};
use fs_exec::Cx;
use fs_project::{ConductionTransient, ProjectSpec, ThermalBoundaryCondition};

use super::{
    SolveRefusal, canonical_f64, conduction_error, json_string, ladder_functional,
    ladder_target_region, conjugate as airflow, natural as buoyancy, radiation,
};

mod adjoint;
mod conjugate;
mod enthalpy;
mod enclosure;
mod natural;
mod radiative;
use radiative::PhysicalEnergy;
mod sampled_peak;
use sampled_peak::SampledPeak;
pub(super) mod spatial;
mod workload;
use workload::{TimeGrids, Workload};

pub(super) const NO_CLAIM: &str = "Estimated finite-mesh final-time temperature from backward Euler with declared temperature-independent regional heat capacity and the bound conductivity evaluated at each endpoint. Temperature-dependent conductivity uses residual-gated Newton/FGMRES with the full k'(T) tangent and immutable physical history. Explicit regional power histories replace delivered watts on their named volume regions; both grids land on every workload switch. A forward sampled peak covers the original initial field and every accepted fine endpoint in the requirement region, retaining the earliest exact tie. It is a sampled-time observation and does not bound unsampled temperatures. The nested-grid final-maximum difference estimates temporal error at assumed order one; it is not a spatial error bound, observed-order proof, continuum or continuous-time maximum bound, validated heat capacity, or compliance certificate. Prescribed temperatures are endpoint data, including their discrete boundary-node storage reaction. Declared ambient radiation uses the shared area-mean gray-patch model. Declared closed gray-diffuse enclosures use the shared reflected-radiosity owner with genuine declared volumetric capacity, original convection and prescribed-flux heat, and an independent fresh endpoint residual and storage-energy gate. Every coupling trial retains the same physical history. Declared airflow uses quasi-steady branch transport at the retained flow-network operating point, with solid storage, a full nodal residual gate at actual marched references and a separate actual-air energy gate. Declared natural convection evaluates its admitted card at each endpoint and checks the full physical free residual and energy balance; exact equilibrium and out-of-domain Rayleigh states remain refused. Natural reservoirs and forced-air paths on distinct patches compose through the same stored solid: every natural update encloses a complete air solve, and the outer endpoint is independently checked under actual natural coefficients and actual air-law references with a separate air-enthalpy energy gate. The final air row describes its last natural iteration; natural and trajectory counters include all nested responses without advancing physical history. No airflow storage, time-varying fan dynamics, latent heat or undeclared transient adjoint is inferred. Endpoint heat storage is reported separately from numerical energy closure; transient energy_residual_j is the checked storage-minus-net-input balance.";

pub(super) struct NativeTransientResult {
    pub(super) solution: radiation::SolidSolution,
    pub(super) receipt: String,
    pub(super) conjugate_receipt: Option<String>,
    pub(super) natural_receipt: Option<String>,
    pub(super) derived_boundary: BTreeMap<String, (f64, f64)>,
    /// An assumed-first-order temporal estimate for the published fine QoI.
    /// Exact floating-point agreement does not establish zero temporal error.
    pub(super) temporal_half_width_k: Option<f64>,
    /// Actual final-step heat entering storage, ΔU/dt, W.
    pub(super) endpoint_storage_w: f64,
    /// Actual final-step (storage - net input) residual divided by dt, W.
    pub(super) endpoint_energy_residual_w: f64,
}

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-transient",
        message,
        "declare either regional capacities or enthalpy charts with admitted material and boundary laws; inspect the time/work/energy declaration",
    )
}

fn lower(error: ConductionError) -> SolveRefusal {
    match error {
        ConductionError::Cancelled { .. } => conduction_error(
            "cli-solve-cancelled",
            "transient conduction was cancelled before final publication",
            "resume the retained stage prefix; this conduction stage restarts from its declared initial state",
        ),
        other => bad(format!("backward-Euler thermal solve refused: {other}")),
    }
}

fn uniform_initial_temperature(policy: &ConductionTransient) -> Result<f64, SolveRefusal> {
    policy.initial_temperature.uniform_temperature()
        .map(|quantity| quantity.value)
        .ok_or_else(|| bad("regional initial temperatures require explicit total-enthalpy storage"))
}

fn number(value: f64) -> Result<String, SolveRefusal> {
    canonical_f64(value).ok_or_else(|| bad("transient evidence contains nonfinite arithmetic"))
}

fn poll(cx: &Cx<'_>, deadline: Option<(Instant, f64)>) -> Result<(), SolveRefusal> {
    cx.checkpoint().map_err(|_| {
        conduction_error(
            "cli-solve-cancelled",
            "transient conduction was cancelled",
            "resume the retained pipeline prefix",
        )
    })?;
    check_deadline(deadline)
}

pub(super) fn check_deadline(deadline: Option<(Instant, f64)>) -> Result<(), SolveRefusal> {
    if deadline.is_some_and(|(started, seconds)| started.elapsed().as_secs_f64() >= seconds) {
        return Err(conduction_error(
            "cli-solve-conduction-transient-wall-budget",
            "transient conduction exhausted the remaining wall budget between numerical operations",
            "increase the solve-time budget; preceding stages remain retained and this stage restarts",
        ));
    }
    Ok(())
}

fn policy(spec: &ProjectSpec) -> Option<&ConductionTransient> {
    spec.cooling
        .as_ref()?
        .conduction
        .as_ref()?
        .transient
        .as_ref()
}

pub(super) fn requested(spec: &ProjectSpec) -> bool {
    policy(spec).is_some()
}

pub(super) fn no_claim(spec: &ProjectSpec) -> &'static str {
    if policy(spec).is_some_and(|time| time.enthalpy.is_some()) {
        enthalpy::no_claim()
    } else {
        NO_CLAIM
    }
}

/// Refuse incompatible model combinations before meshing or solving.
pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    // Also reject a transient-adjoint output on a steady project; an ignored
    // report request must never look like a completed derivative calculation.
    adjoint::requested(spec)?;
    if radiation::pointwise_requested(spec) {
        if !policy(spec).is_some_and(|time| time.enthalpy.is_some()) {
            return Err(conduction_error(
                "cli-solve-conduction-pointwise-radiation",
                "pointwise radiation requires the native total-enthalpy transient owner",
                "declare total-enthalpy storage and finite-time integration, or select patch-mean radiation",
            ));
        }
        if spec.outputs.iter().flatten().any(|output| {
            output.name.contains("adjoint") || output.name == "temperature-volume-mean-bound"
        }) {
            return Err(conduction_error(
                "cli-solve-conduction-pointwise-radiation",
                "pointwise radiation does not yet expose the native trajectory adjoint or continuum mean bound",
                "request the forward temperature field or sampled peak; no mean-patch derivative is substituted",
            ));
        }
    }
    let Some(transient) = policy(spec) else {
        return Ok(());
    };
    if let Some(violation) = spec
        .validate()
        .into_iter()
        .find(|row| row.code.starts_with("project-conduction-transient"))
    {
        return Err(conduction_error(
            violation.code,
            violation.what,
            violation.fix,
        ));
    }
    if spec
        .solver
        .as_ref()
        .is_some_and(|solver| solver.fidelity == "adaptive")
    {
        return Err(bad(
            "transient adaptive refinement requires a trajectory goal estimator; use the uniform spatial ladder with its shared time-step allowance",
        ));
    }
    if spec
        .outputs
        .iter()
        .flatten()
        .any(|output| output.name.contains("adjoint") && !adjoint::is_output(&output.name))
    {
        return Err(bad(
            "steady nominal adjoints cannot differentiate a transient trajectory",
        ));
    }
    enthalpy::admit(spec)?;
    TimeGrids::new(transient)?;
    Ok(())
}

struct GridResult {
    endpoint: StepSolution,
    endpoint_dt_s: f64,
    endpoint_energy: PhysicalEnergy,
    rows: String,
    sampled_peak: SampledPeak,
    /// Present only on an explicitly requested fine-grid derivative path.
    tape: Option<adjoint::Tape>,
    linear: Vec<LinearSolveEvidence>,
    stored_j: f64,
    net_input_j: f64,
    maximum_energy_residual_j: f64,
    nonlinear_updates: usize,
    nonlinear_backtracks: usize,
    radiation: Option<RadiativeEndpoint>,
    radiation_trials: usize,
    maximum_physical_energy_residual_j: f64,
    airflow: Option<AirflowEndpoint>,
    airflow_iterations: usize,
    airflow_solid_solves: usize,
    airflow_endpoint_evaluations: usize,
    natural: Option<NaturalEndpoint>,
    natural_iterations: usize,
    natural_solid_solves: usize,
    natural_endpoint_evaluations: usize,
    coupled_net_input_j: f64,
    maximum_coupled_energy_residual_j: f64,
}

/// The accepted solid field belongs to these applied references, not the
/// updated air references returned by the final fixed-point evaluation.
struct AirflowEndpoint {
    outcome: airflow::ConjugateOutcome,
    applied_references: BTreeMap<String, f64>,
    boundary: fs_conduction::ThermalBoundary,
}

/// The final physical natural law is retained beside the coefficients the
/// counted solid response actually used.
struct NaturalEndpoint {
    receipt: String,
    applied_coefficients: BTreeMap<String, f64>,
    boundary: fs_conduction::ThermalBoundary,
}

/// Retained final physical boundary. Inner secant and nonlinear physical heat
/// remain distinct, so neither a frozen operator nor a prescribed reaction is
/// silently relabeled as the full radiative model.
enum RadiativeEndpoint {
    Ambient {
        combined_boundary: fs_conduction::ThermalBoundary,
        convective_robin_fluxes: Vec<fs_conduction::RobinFlux>,
        convective_out_w: f64,
        report: fs_conduction::AmbientRadiationReport,
    },
    Enclosure(radiation::EnclosureEndpoint),
}

impl RadiativeEndpoint {
    fn convective_robin_fluxes(&self) -> &[fs_conduction::RobinFlux] {
        match self {
            Self::Ambient { convective_robin_fluxes, .. } => convective_robin_fluxes,
            Self::Enclosure(endpoint) => &endpoint.convective_robin_fluxes,
        }
    }

    fn convective_out_w(&self) -> f64 {
        match self {
            Self::Ambient { convective_out_w, .. } => *convective_out_w,
            Self::Enclosure(endpoint) => endpoint.convective_out_w,
        }
    }
}

fn nonlinear_row(solved: &NonlinearStepSolution) -> Result<String, SolveRefusal> {
    Ok(format!(
        "{{\"updates\":{},\"backtracks\":{},\"initial_residual_j\":{},\"residual_j\":{},\"threshold_j\":{}}}",
        solved.nonlinear_iterations,
        solved.backtracks,
        number(solved.initial_residual_j)?,
        number(solved.residual_j)?,
        number(solved.threshold_j)?,
    ))
}

fn radiation_row(solved: &RadiationStepSolution) -> Result<String, SolveRefusal> {
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
        "{{\"iterations\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"max_temperature_change_k\":{},\"max_heat_mismatch_w\":{},\"decomposition_residual_w\":{},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"physical_dirichlet_in_w\":{},\"patches\":[{}]}}",
        report.iterations, report.solid_iterations, report.krylov_iterations,
        number(solved.convective_out_w)?, number(report.applied_radiation_out_w)?,
        number(report.nonlinear_radiation_out_w)?, number(report.max_temperature_change_k)?,
        number(report.max_heat_mismatch_w)?, number(report.decomposition_residual_w)?,
        number(solved.physical_residual_norm_j)?, number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?, number(solved.physical_dirichlet_in_w)?,
        patches.join(","),
    ))
}

fn airflow_row(
    path: &airflow::ConjugatePath,
    solved: &conjugate::AirflowStep,
) -> Result<String, SolveRefusal> {
    let references = solved.applied_references.iter()
        .map(|(target, temperature)| {
            Ok(format!("{}:{}", json_string(target), number(*temperature)?))
        })
        .collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        "{{\"air_heat_gain_w\":{},\"off_path_convective_out_w\":{},\"radiative_out_w\":{},\"physical_dirichlet_in_w\":{},\"coupled_energy_residual_j\":{},\"solid_solves\":{},\"krylov_iterations\":{},\"nonlinear_iterations\":{},\"nonlinear_backtracks\":{},\"radiation_trials\":{},\"reference_tolerance_k\":{},\"endpoint_evaluations\":{},\"initial_physical_residual_norm_j\":{},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"applied_references_k\":{{{}}},\"exchange\":{}}}",
        number(solved.air_heat_gain_w)?,
        number(solved.off_path_convective_out_w)?,
        number(solved.radiation_out_w)?,
        number(solved.physical_dirichlet_in_w)?,
        number(solved.coupled_energy_residual_j)?,
        solved.solid_solves, solved.krylov_iterations, solved.nonlinear_updates,
        solved.nonlinear_backtracks, solved.radiation_trials,
        number(solved.reference_tolerance_k)?, solved.endpoint_evaluations,
        number(solved.initial_physical_residual_norm_j)?,
        number(solved.physical_residual_norm_j)?, number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?, references.join(","),
        airflow::receipt_fragment(path, &solved.outcome)?,
    ))
}

fn natural_row(solved: &natural::NaturalStep) -> Result<String, SolveRefusal> {
    let coefficients = |values: &BTreeMap<String, f64>| {
        values.iter().map(|(target, value)| {
            Ok(format!("{}:{}", json_string(target), number(*value)?))
        }).collect::<Result<Vec<_>, SolveRefusal>>().map(|rows| rows.join(","))
    };
    let fluxes = solved.physical_convective_robin_fluxes.iter().map(|flux| {
        Ok(format!(
            "{{\"region\":{},\"area_m2\":{},\"heat_rate_w\":{}}}",
            json_string(&flux.region), number(flux.area_m2)?, number(flux.heat_rate_w)?,
        ))
    }).collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        "{{\"iterations\":{},\"solid_solves\":{},\"krylov_iterations\":{},\"nonlinear_iterations\":{},\"nonlinear_backtracks\":{},\"radiation_trials\":{},\"endpoint_evaluations\":{},\"airflow_iterations\":{},\"coupled_net_input_w\":{},\"coupled_energy_residual_j\":{},\"applied_coefficients_w_m2_k\":{{{}}},\"physical_coefficients_w_m2_k\":{{{}}},\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},\"physical_energy_residual_j\":{},\"physical_dirichlet_in_w\":{},\"physical_convective_out_w\":{},\"physical_radiative_out_w\":{},\"physical_convective_robin_fluxes\":[{}],\"exchange\":{}}}",
        solved.iterations, solved.solid_solves, solved.krylov_iterations,
        solved.nonlinear_updates, solved.nonlinear_backtracks, solved.radiation_trials,
        solved.endpoint_evaluations, solved.airflow_iterations,
        number(solved.coupled_net_input_w)?, number(solved.coupled_energy_residual_j)?,
        coefficients(&solved.applied_coefficients)?,
        coefficients(&solved.physical_coefficients)?,
        number(solved.physical_residual_norm_j)?, number(solved.physical_residual_tolerance_j)?,
        number(solved.physical_energy_residual_j)?, number(solved.physical_dirichlet_in_w)?,
        number(solved.physical_convective_out_w)?, number(solved.physical_radiative_out_w)?,
        fluxes.join(","), solved.receipt,
    ))
}

#[allow(clippy::too_many_arguments)]
fn march(
    cx: &Cx<'_>,
    engine: &BackwardEuler<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    policy: &ConductionTransient,
    ends_s: &[f64],
    workload: Option<&Workload<'_>>,
    linear: LinearConfig,
    nonlinear: Option<NonlinearStepConfig>,
    radiation: Option<&radiation::LoweredRadiation>,
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    airflow: Option<&airflow::ConjugatePath>,
    natural: &[buoyancy::NaturalLaw],
    pressure_pa: f64,
    labels: &[u32],
    goal_region: Option<u32>,
    retain_adjoint: bool,
    deadline: Option<(Instant, f64)>,
) -> Result<GridResult, SolveRefusal> {
    let steps = ends_s.len();
    let mut old = vec![uniform_initial_temperature(policy)?; problem.mesh.vertex_count()];
    let mut sampled_peak = SampledPeak::new(ladder_functional(
        labels, &problem.mesh.complex().tets, &old, goal_region,
    ))?;
    let mut tape = if retain_adjoint {
        Some(adjoint::Tape::new(cx, &old, steps, deadline)?)
    } else { None };
    let mut rows = Vec::with_capacity(steps);
    let mut evidence = Vec::with_capacity(steps);
    let mut last = None;
    let mut time = 0.0;
    let mut stored_j = 0.0;
    let mut net_input_j = 0.0;
    let mut maximum_energy_residual_j = 0.0_f64;
    let mut nonlinear_updates = 0;
    let mut nonlinear_backtracks = 0;
    let mut endpoint_dt_s = 0.0;
    let mut endpoint_energy = None;
    let mut last_radiation = None;
    let mut radiation_trials = 0;
    let mut maximum_physical_energy_residual_j = 0.0_f64;
    let mut last_airflow = None;
    let mut airflow_iterations = 0;
    let mut airflow_solid_solves = 0;
    let mut airflow_endpoint_evaluations = 0;
    let mut last_natural = None;
    let mut natural_iterations = 0;
    let mut natural_solid_solves = 0;
    let mut natural_endpoint_evaluations = 0;
    let mut coupled_net_input_j = 0.0;
    let mut maximum_coupled_energy_residual_j = 0.0_f64;
    for (index, &end) in ends_s.iter().enumerate() {
        let ordinal = index + 1;
        poll(cx, deadline)?;
        let dt = end - time;
        if !(end.is_finite() && dt.is_finite() && dt > 0.0) {
            return Err(bad(
                "the requested time grid cannot represent every positive step",
            ));
        }
        let config = StepConfig {
            linear,
            energy_tolerance_j: policy.energy_tolerance.value,
        };
        let source = workload
            .map(|history| history.source(cx, problem.mesh, labels, time, deadline))
            .transpose()?;
        let step_problem = ConductionProblem {
            source: source.as_ref().unwrap_or(problem.source),
            ..problem
        };
        let mut radiative_row = "null".to_string();
        let mut conjugate_row = "null".to_string();
        let mut natural_evidence = "null".to_string();
        let mut physical_target_j = None;
        let mut physical_energy = None;
        let (endpoint, nonlinear_evidence, krylov_iterations) = if !natural.is_empty() {
            let solved = natural::advance(
                cx, engine, step_problem, interfaces, &old, dt, config,
                nonlinear, radiation, enclosure, airflow, natural, pressure_pa, deadline,
            )?;
            natural_evidence = natural_row(&solved)?;
            physical_energy = Some(solved.physical_energy);
            // The outer natural owner accepted the complete physical model.
            // An inner radiation target belongs to its frozen natural trial.
            physical_target_j = Some(solved.physical_residual_tolerance_j);
            natural_iterations += solved.iterations;
            natural_solid_solves += solved.solid_solves;
            natural_endpoint_evaluations += solved.endpoint_evaluations;
            nonlinear_updates += solved.nonlinear_updates;
            nonlinear_backtracks += solved.nonlinear_backtracks;
            radiation_trials += solved.radiation_trials;
            maximum_physical_energy_residual_j = maximum_physical_energy_residual_j
                .max(solved.maximum_physical_energy_residual_j);
            maximum_coupled_energy_residual_j = maximum_coupled_energy_residual_j
                .max(solved.physical_energy_residual_j.abs())
                .max(solved.coupled_energy_residual_j.abs());
            coupled_net_input_j += dt * solved.coupled_net_input_w;
            conjugate_row = solved.airflow_evidence;
            if solved.airflow.is_some() {
                airflow_iterations += solved.airflow_iterations;
                // These are the same solves counted by natural_solid_solves,
                // presented under both owning interfaces, never summed twice.
                airflow_solid_solves += solved.solid_solves;
                airflow_endpoint_evaluations += solved.airflow_endpoint_evaluations;
            }
            last_airflow = solved.airflow;
            radiative_row = solved.radiation_evidence;
            last_radiation = solved.radiation;
            last_natural = Some(NaturalEndpoint {
                receipt: solved.receipt,
                applied_coefficients: solved.applied_coefficients,
                boundary: solved.applied_boundary,
            });
            (solved.endpoint, solved.nonlinear_evidence, solved.krylov_iterations)
        } else if let Some(path) = airflow {
            let solved = conjugate::advance(
                cx, engine, step_problem, interfaces, &old, dt, config,
                nonlinear, radiation, enclosure, path, deadline,
            )?;
            conjugate_row = airflow_row(path, &solved)?;
            physical_energy = Some(solved.physical_energy);
            physical_target_j = Some(solved.physical_residual_tolerance_j);
            airflow_iterations += solved.outcome.solution.iterations;
            airflow_solid_solves += solved.solid_solves;
            airflow_endpoint_evaluations += solved.endpoint_evaluations;
            nonlinear_updates += solved.nonlinear_updates;
            nonlinear_backtracks += solved.nonlinear_backtracks;
            radiation_trials += solved.radiation_trials;
            maximum_physical_energy_residual_j = maximum_physical_energy_residual_j
                .max(solved.maximum_physical_energy_residual_j);
            maximum_coupled_energy_residual_j = maximum_coupled_energy_residual_j
                .max(solved.coupled_energy_residual_j.abs());
            coupled_net_input_j += dt * (solved.endpoint.source_w
                + solved.physical_dirichlet_in_w - solved.physical_energy.neumann_out_w
                - solved.off_path_convective_out_w - solved.air_heat_gain_w
                - solved.radiation_out_w);
            radiative_row = solved.radiation_evidence;
            last_radiation = solved.radiation;
            last_airflow = Some(AirflowEndpoint {
                outcome: solved.outcome,
                applied_references: solved.applied_references,
                boundary: solved.applied_boundary,
            });
            (solved.endpoint, solved.nonlinear_evidence, solved.krylov_iterations)
        } else if let Some(radiation) = radiation {
            let solved = radiative::advance(
                cx, engine, step_problem, interfaces, &old, dt, config,
                nonlinear, radiation, enclosure, deadline,
            )?;
            radiative_row = solved.radiation_evidence;
            physical_target_j = Some(solved.physical_residual_tolerance_j);
            physical_energy = Some(solved.physical_energy);
            radiation_trials += solved.solid_solves;
            nonlinear_updates += solved.nonlinear_updates;
            nonlinear_backtracks += solved.nonlinear_backtracks;
            maximum_physical_energy_residual_j = maximum_physical_energy_residual_j
                .max(solved.physical_energy.energy_residual_j.abs());
            last_radiation = Some(solved.radiation);
            (solved.conduction, solved.nonlinear_evidence, solved.krylov_iterations)
        } else if let Some(nonlinear) = nonlinear {
            let solved = engine.advance_nonlinear_prescribed(
                cx, step_problem, interfaces, &old, dt, config, nonlinear,
            ).map_err(lower)?;
            nonlinear_updates += solved.nonlinear_iterations;
            nonlinear_backtracks += solved.backtracks;
            let row = nonlinear_row(&solved)?;
            let krylov_iterations = solved.step.krylov_iterations;
            (solved.step, row, krylov_iterations)
        } else {
            let endpoint = engine
                .advance_prescribed(cx, step_problem, interfaces, &old, dt, config)
                .map_err(lower)?;
            let krylov_iterations = endpoint.krylov_iterations;
            (endpoint, "null".to_string(), krylov_iterations)
        };
        // Conductivity was assembled at element means; keep the whole
        // published field inside every consuming material's retained span.
        for (element, tet) in problem.mesh.complex().tets.iter().enumerate() {
            if element % 256 == 0 {
                poll(cx, deadline)?;
            }
            let material = match problem.element_materials {
                Some(materials) => materials.model_for(element).map_err(lower)?,
                None => problem.material,
            };
            for &vertex in tet {
                let temperature = endpoint.temperature[vertex as usize];
                if !temperature.is_finite() || temperature <= 0.0 {
                    return Err(bad(
                        "a transient endpoint is not a positive finite absolute temperature",
                    ));
                }
                material
                    .temperature_span()
                    .check(temperature)
                    .map_err(lower)?;
            }
        }
        let qoi = ladder_functional(
            labels,
            &problem.mesh.complex().tets,
            &endpoint.temperature,
            goal_region,
        );
        sampled_peak.record(end, qoi)?;
        // Enclosure heat is internal exchange. Common trajectory energy must
        // use the fresh physical bath, original Neumann/convection and closed
        // reflected heat, while the inner StepSolution retains its applied law.
        let energy = if enclosure.is_some() {
            physical_energy.ok_or_else(|| bad("enclosure endpoint lost its physical energy audit"))?
        } else {
            PhysicalEnergy::inner(&endpoint)
        };
        let net_w = energy.net_input_w;
        stored_j += endpoint.stored_energy_change_j;
        net_input_j += dt * net_w;
        maximum_energy_residual_j = maximum_energy_residual_j.max(energy.energy_residual_j.abs());
        rows.push(format!("{{\"step\":{ordinal},\"time_s\":{},\"dt_s\":{},\"final_region_max_k\":{},\"stored_energy_change_j\":{},\"net_input_w\":{},\"source_w\":{},\"energy_residual_j\":{},\"relative_residual\":{},\"krylov_iterations\":{},\"nonlinear\":{nonlinear_evidence},\"radiation\":{radiative_row},\"conjugate\":{conjugate_row},\"natural\":{natural_evidence}}}",
            number(end)?, number(dt)?, number(qoi)?, number(endpoint.stored_energy_change_j)?,
            number(net_w)?, number(endpoint.source_w)?, number(energy.energy_residual_j)?, number(endpoint.relative_residual)?,
            krylov_iterations));
        evidence.push(LinearSolveEvidence {
            nonlinear_iteration: ordinal - 1,
            method: if !natural.is_empty() && airflow.is_some() {
                match (nonlinear.is_some(), radiation.is_some()) {
                    (true, true) => "fgmres-backward-euler-newton-radiation-mixed-convection",
                    (true, false) => "fgmres-backward-euler-newton-mixed-convection",
                    (false, true) => "pcg-backward-euler-radiation-mixed-convection",
                    (false, false) => "pcg-backward-euler-mixed-convection",
                }
            } else if !natural.is_empty() {
                match (nonlinear.is_some(), radiation.is_some()) {
                    (true, true) => "fgmres-backward-euler-newton-radiation-natural",
                    (true, false) => "fgmres-backward-euler-newton-natural",
                    (false, true) => "pcg-backward-euler-radiation-natural",
                    (false, false) => "pcg-backward-euler-natural",
                }
            } else { match (nonlinear.is_some(), radiation.is_some(), airflow.is_some()) {
                (true, true, true) => "fgmres-backward-euler-newton-radiation-airflow",
                (true, false, true) => "fgmres-backward-euler-newton-airflow",
                (false, true, true) => "pcg-backward-euler-radiation-airflow",
                (false, false, true) => "pcg-backward-euler-airflow",
                (true, true, false) => "fgmres-backward-euler-newton-radiation",
                (true, false, false) => "fgmres-backward-euler-newton",
                (false, true, false) => "pcg-backward-euler-radiation",
                (false, false, false) => "pcg-backward-euler",
            } },
            iterations: krylov_iterations,
            reported: ResidualClaim::TrueEuclidean(endpoint.relative_residual),
            true_relative_residual: endpoint.relative_residual,
            converged_true: endpoint.relative_residual < linear.tolerance,
            stall: None,
        });
        if let Some(tape) = &mut tape {
            match physical_target_j {
                Some(target) => tape.record_physical(cx, &endpoint.temperature, target, deadline)?,
                None => tape.record(cx, &endpoint.temperature, deadline)?,
            }
        }
        old.clone_from(&endpoint.temperature);
        endpoint_dt_s = dt;
        endpoint_energy = Some(energy);
        last = Some(endpoint);
        time = end;
    }
    poll(cx, deadline)?;
    Ok(GridResult {
        endpoint: last.ok_or_else(|| bad("the transient grid has no steps"))?,
        endpoint_dt_s,
        endpoint_energy: endpoint_energy.ok_or_else(|| bad("transient grid lost its final energy audit"))?,
        rows: rows.join(","),
        sampled_peak,
        tape,
        linear: evidence,
        stored_j,
        net_input_j,
        maximum_energy_residual_j,
        nonlinear_updates,
        nonlinear_backtracks,
        radiation: last_radiation,
        radiation_trials,
        maximum_physical_energy_residual_j,
        airflow: last_airflow,
        airflow_iterations,
        airflow_solid_solves,
        airflow_endpoint_evaluations,
        natural: last_natural,
        natural_iterations,
        natural_solid_solves,
        natural_endpoint_evaluations,
        coupled_net_input_j,
        maximum_coupled_energy_residual_j,
    })
}

/// Lower only declared capacities and execute the shared spatial problem.
#[allow(clippy::too_many_arguments)]
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
    let objective = adjoint::requested(spec)?;
    let policy = policy(spec).ok_or_else(|| bad("no transient declaration was supplied"))?;
    if labels.len() != problem.mesh.element_count() {
        return Err(bad(
            "transient regional capacities do not match the retained mesh labels",
        ));
    }
    if radiation.is_some()
        != spec
            .cooling
            .as_ref()
            .and_then(|cooling| cooling.conduction.as_ref())
            .is_some_and(|setup| setup.radiation.is_some() || setup.enclosure_radiation.is_some())
    {
        return Err(bad(
            "the transient radiation declaration has no matching material-card lowering",
        ));
    }
    if airflow.is_some()
        != spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
            .is_some_and(|setup| setup.boundaries.iter().any(|row| {
                matches!(row.condition, ThermalBoundaryCondition::AirflowConvection { .. })
            }))
    {
        return Err(bad(
            "the transient airflow declaration has no matching flow-network and wetted-area lowering",
        ));
    }
    let natural_requested = !natural.is_empty();
    if natural_requested
        != spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
            .is_some_and(|setup| setup.boundaries.iter().any(|row| {
                matches!(row.condition, ThermalBoundaryCondition::NaturalConvection { .. })
            }))
    {
        return Err(bad(
            "the transient natural-convection declaration has no matching card lowering",
        ));
    }
    if policy.enthalpy.is_some() {
        return enthalpy::solve(
            cx, spec, problem, interfaces, labels, region_ids, radiation,
            airflow, natural, pressure_pa, linear, deadline,
        );
    }
    let initial_temperature = uniform_initial_temperature(policy)?;
    let grids = TimeGrids::new(policy)?;
    let coarse_steps = grids.coarse.len();
    let fine_steps = grids.fine.len();
    let total_steps = coarse_steps
        .checked_add(fine_steps)
        .ok_or_else(|| bad("transient step count overflow"))?;
    if total_steps > policy.max_steps as usize || total_steps > 10_000 {
        return Err(bad(
            "coarse and fine grids exceed the admitted combined time-work cap",
        ));
    }
    // Includes endpoint/history/scratch vectors, retained time rows and an
    // element-scaled sparse assembly allowance. This is checked admission,
    // not a claim to meter allocator overhead or peak process RSS.
    let row_bytes = radiation
        .map_or(Some(2048_usize), |lowering| {
            lowering
                .surface_count()
                .checked_mul(1024)
                .and_then(|bytes| bytes.checked_add(4096))
        })
        .and_then(|row_bytes| match airflow {
            Some(path) => path.segments.len().checked_mul(4096)
                .and_then(|bytes| bytes.checked_add(4096))
                .and_then(|bytes| bytes.checked_add(row_bytes)),
            None => Some(row_bytes),
        })
        .and_then(|row_bytes| if natural_requested {
            natural.len().checked_mul(4096)
                .and_then(|bytes| bytes.checked_add(4096))
                .and_then(|bytes| bytes.checked_add(row_bytes))
        } else { Some(row_bytes) })
        .ok_or_else(|| bad("transient coupling row estimate overflows"))?;
    let bytes = problem
        .mesh
        .vertex_count()
        .checked_mul((if natural_requested { 192 } else { 128 })
            + if airflow.is_some() { 64 } else { 0 })
        .and_then(|n| {
            problem
                .mesh
                .element_count()
                .checked_mul(4096)
                .and_then(|e| n.checked_add(e))
        })
        .and_then(|n| total_steps.checked_mul(row_bytes).and_then(|s| n.checked_add(s)))
        .ok_or_else(|| bad("transient memory estimate overflows"))?;
    // The enclosure owner retains dense patch-space radiosity work, bounded
    // independently from the sparse solid. Include it in forward admission.
    let bytes = if let Some(lowering) = radiation {
        lowering.surface_count().checked_mul(lowering.surface_count())
            .and_then(|entries| entries.checked_mul(128))
            .and_then(|additional| bytes.checked_add(additional))
            .ok_or_else(|| bad("transient radiation workspace estimate overflows"))?
    } else { bytes };
    let bytes = if objective.is_some() {
        let mut additional = adjoint::storage_bytes(problem.mesh.vertex_count(), fine_steps,
            policy.capacities.len())?;
        if objective.is_some_and(adjoint::Objective::material_controls) {
            additional = additional.checked_add(adjoint::material_control_storage_bytes(
                spec, problem.mesh.vertex_count(), problem.mesh.element_count(),
            )?).ok_or_else(|| bad("combined material-control adjoint memory estimate overflows"))?;
        }
        let ports = problem.boundary.conditions().iter()
            .filter(|condition| matches!(condition, fs_conduction::ThermalBc::Robin { .. }))
            .count();
        if let Some(lowering) = radiation {
            additional = additional.checked_add(adjoint::radiation_storage_bytes(
                problem.mesh.vertex_count(), ports, lowering.surface_count(), lowering.is_enclosure(),
            )?).ok_or_else(|| bad("combined radiative adjoint memory estimate overflows"))?;
        }
        if natural_requested {
            additional = additional.checked_add(adjoint::natural_storage_bytes(
                problem.mesh.vertex_count(), ports, natural.len(),
            )?).ok_or_else(|| bad("combined natural-convection adjoint memory estimate overflows"))?;
        }
        if let Some(path) = airflow {
            additional = additional.checked_add(adjoint::airflow_storage_bytes(
                spec, problem.mesh.vertex_count(), ports, path,
            )?).ok_or_else(|| bad("combined airflow adjoint memory estimate overflows"))?;
        }
        bytes.checked_add(additional)
            .ok_or_else(|| bad("combined forward/adjoint memory estimate overflows"))?
    } else { bytes };
    if objective.is_some() {
        let reverse_solves = radiation.map_or(Some(1), |lowering| {
            lowering.surface_count().checked_add(2)
        }).ok_or_else(|| bad("reverse radiation work allowance overflows"))?;
        if fine_steps.checked_mul(reverse_solves)
            .and_then(|solves| solves.checked_mul(linear.max_iterations)).is_none()
        {
            return Err(bad("reverse reconstruction and transpose work allowances overflow"));
        }
    }
    if u64::try_from(bytes).unwrap_or(u64::MAX)
        > spec.budgets.as_ref().map_or(0, |b| b.memory_bytes)
    {
        return Err(bad(format!(
            "transient storage and retained step rows need an estimated {bytes} bytes, above the declared memory budget"
        )));
    }
    let mut by_region = BTreeMap::new();
    let mut capacity_rows = Vec::new();
    for row in &policy.capacities {
        let id = *region_ids
            .get(&row.region)
            .ok_or_else(|| bad("capacity names an unknown spatial region"))?;
        if by_region
            .insert(
                id,
                VolumetricHeatCapacity::declared(row.volumetric_heat_capacity.value)
                    .map_err(lower)?,
            )
            .is_some()
        {
            return Err(bad("a transient region has more than one capacity"));
        }
        capacity_rows.push(format!(
            "{{\"region\":{},\"volumetric_heat_capacity_j_m3_k\":{},\"source\":{}}}",
            json_string(&row.region),
            number(row.volumetric_heat_capacity.value)?,
            json_string(&row.source)
        ));
    }
    let capacities = labels
        .iter()
        .map(|id| {
            by_region
                .get(id)
                .copied()
                .ok_or_else(|| bad("a retained mesh region has no declared heat capacity"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut temperature_dependent = false;
    for element in 0..problem.mesh.element_count() {
        if element % 256 == 0 {
            poll(cx, deadline)?;
        }
        let material = match problem.element_materials {
            Some(materials) => materials.model_for(element).map_err(lower)?,
            None => problem.material,
        };
        temperature_dependent |= material.is_temperature_dependent();
        material
            .temperature_span()
            .check(initial_temperature)
            .map_err(lower)?;
    }
    // The source-bound material selects the numerical path. Newton tolerances
    // measure the actual transient residual in joules; absolute temperature
    // never supplies the residual scale. Reserve 99% of the energy allowance
    // for the independently checked relative-residual and rounding effects.
    // Both grids reuse these same bounded controls, retained below.
    let nonlinear = if temperature_dependent {
        let mut config = NonlinearStepConfig::default();
        if let Some(solver) = &spec.solver {
            config.residual_rtol = solver.tolerance_rel;
        }
        config.residual_atol_j = 0.01 * policy.energy_tolerance.value
            / fs_math::det::sqrt(problem.mesh.vertex_count().max(1) as f64);
        config.validate().map_err(lower)?;
        Some(config)
    } else {
        None
    };
    let goal_region = ladder_target_region(spec, region_ids);
    if goal_region.is_none() {
        return Err(bad(
            "transient final-time error comparison needs the declared temperature-maximum requirement region",
        ));
    }
    let goal_region_name = spec.requirements.iter().flatten()
        .find(|row| region_ids.get(&row.region).copied() == goal_region)
        .map(|row| row.region.as_str())
        .ok_or_else(|| bad("transient trajectory lost its declared requirement region"))?;
    // Either convection owner bounds interface evaluations, radiation bounds
    // repeated solid responses, and the linear allowance includes all
    // corrections in one response. Check complete both-grid products first.
    let airflow_factor = airflow.map_or(1, |_| airflow::exchange_config(false).max_iterations);
    let natural_factor = if natural_requested { buoyancy::MAX_ITERATIONS } else { 1 };
    let radiation_factor = radiation.map_or(1, |lowering| lowering.config.max_iterations);
    let max_endpoint_evaluations_per_step = if natural_requested {
        natural_factor.checked_mul(2)
            .ok_or_else(|| bad("natural endpoint evaluation allowance overflows"))?
    } else { 0 };
    let max_airflow_endpoint_evaluations_per_step = if airflow.is_some() {
        natural_factor.checked_mul(2)
            .ok_or_else(|| bad("airflow endpoint evaluation allowance overflows"))?
    } else { 0 };
    let max_airflow_iterations_per_step = airflow_factor.checked_mul(natural_factor)
        .ok_or_else(|| bad("coupled interface iteration allowance overflows"))?;
    let max_solid_solves_per_step = max_airflow_iterations_per_step.checked_mul(radiation_factor)
        .ok_or_else(|| bad("coupled solid response allowance overflows"))?;
    let max_krylov_per_step = linear
        .max_iterations
        .checked_mul(max_solid_solves_per_step)
        .ok_or_else(|| bad("coupled Krylov allowance overflows"))?;
    let max_updates_per_step = nonlinear.map_or(Ok(0), |config| {
        config
            .max_iterations
            .checked_mul(max_solid_solves_per_step)
            .ok_or_else(|| bad("coupled Newton allowance overflows"))
    })?;
    let max_backtracks_per_update = nonlinear.map_or(0, |config| config.line_search.max_backtracks);
    if max_solid_solves_per_step.checked_mul(total_steps).is_none()
        || max_endpoint_evaluations_per_step.checked_mul(total_steps).is_none()
        || max_airflow_endpoint_evaluations_per_step.checked_mul(total_steps).is_none()
        || max_krylov_per_step.checked_mul(total_steps).is_none()
        || max_updates_per_step
            .checked_mul(total_steps)
            .and_then(|updates| updates.checked_mul(max_backtracks_per_update))
            .is_none()
    {
        return Err(bad("combined transient numerical work allowance overflows"));
    }
    let engine = BackwardEuler::per_element(cx, problem.mesh, &capacities).map_err(lower)?;
    // One exact geometric binding per spatial rung is shared by both time
    // grids and every outer natural, air, radiation and reverse evaluation.
    let enclosure = radiation.map(|lowering| lowering.bind_enclosure(cx, problem))
        .transpose().map_err(lower)?.flatten();
    let workload = Workload::bind(cx, spec, policy, problem.mesh, labels, region_ids, deadline)?;
    let coarse = march(
        cx,
        &engine,
        problem,
        interfaces,
        policy,
        &grids.coarse,
        workload.as_ref(),
        linear,
        nonlinear,
        radiation,
        enclosure.as_ref(),
        airflow,
        natural,
        pressure_pa,
        labels,
        goal_region,
        false,
        deadline,
    )?;
    let fine = march(
        cx,
        &engine,
        problem,
        interfaces,
        policy,
        &grids.fine,
        workload.as_ref(),
        linear,
        nonlinear,
        radiation,
        enclosure.as_ref(),
        airflow,
        natural,
        pressure_pa,
        labels,
        goal_region,
        objective.is_some(),
        deadline,
    )?;
    let adjoint_fragment = match objective {
        Some(objective) => {
            let tape = fine.tape.as_ref().ok_or_else(|| bad("requested transient adjoint lost its fine history"))?;
            let natural_coupling = natural_requested.then_some(adjoint::NaturalCoupling {
                laws: natural, pressure_pa,
            });
            let derivative = adjoint::finish(cx, spec, &engine, problem, interfaces, &grids.fine,
                workload.as_ref(), StepConfig { linear, energy_tolerance_j: policy.energy_tolerance.value },
                nonlinear, labels, region_ids, goal_region.expect("admitted goal region"),
                radiation, enclosure.as_ref(), natural_coupling.as_ref(), airflow, tape, objective, deadline)?;
            format!(",\"adjoint\":{derivative}")
        }
        None => String::new(),
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
    let difference = (fine_qoi - coarse_qoi).abs();
    let temporal_half_width_k = (difference > 0.0).then_some(1.25 * difference);
    let estimate = temporal_half_width_k
        .map(number)
        .transpose()?
        .unwrap_or_else(|| "null".to_string());
    let temporal_status = if temporal_half_width_k.is_some() {
        "estimated-order-one"
    } else {
        "no-resolved-time-difference"
    };
    let nonlinear_updates = coarse.nonlinear_updates + fine.nonlinear_updates;
    let nonlinear_backtracks = coarse.nonlinear_backtracks + fine.nonlinear_backtracks;
    let nonlinear_receipt = if let Some(config) = nonlinear {
        format!(
            "{{\"method\":\"newton-fgmres\",\"conductivity\":\"endpoint-k(T)-with-k-prime-tangent\",\"max_updates_per_step\":{},\"max_krylov_iterations_per_step\":{},\"residual_rtol\":{},\"residual_atol_j\":{},\"armijo_c\":{},\"line_search_shrink\":{},\"max_backtracks_per_update\":{},\"coarse_updates\":{},\"fine_updates\":{},\"total_updates\":{nonlinear_updates},\"total_backtracks\":{nonlinear_backtracks}}}",
            max_updates_per_step,
            max_krylov_per_step,
            number(config.residual_rtol)?,
            number(config.residual_atol_j)?,
            number(config.line_search.armijo_c)?,
            number(config.line_search.shrink)?,
            config.line_search.max_backtracks,
            coarse.nonlinear_updates,
            fine.nonlinear_updates,
        )
    } else {
        "null".to_string()
    };
    let radiation_trials = coarse.radiation_trials + fine.radiation_trials;
    let radiation_receipt = if let Some(lowering) = radiation {
        let method = if enclosure.is_some() {
            "implicit-closed-gray-diffuse-enclosure"
        } else {
            "implicit-area-mean-gray-radiation"
        };
        format!(
            "{{\"method\":\"{method}\",\"max_trials_per_step\":{},\"max_trials_per_coupling_iteration\":{},\"max_krylov_iterations_per_solid_solve\":{},\"max_krylov_iterations_per_step\":{max_krylov_per_step},\"coarse_trials\":{},\"fine_trials\":{},\"total_trials\":{radiation_trials},\"maximum_physical_energy_residual_j\":{},\"history\":\"unchanged-physical-old-temperature-in-every-trial\",\"nonlinear_detail\":\"per-step nonlinear rows describe the final inner solve; radiation rows include that evaluation's trials; trajectory counters include every outer convection evaluation\"}}",
            max_solid_solves_per_step, lowering.config.max_iterations, linear.max_iterations,
            coarse.radiation_trials, fine.radiation_trials,
            number(coarse.maximum_physical_energy_residual_j
                .max(fine.maximum_physical_energy_residual_j))?,
        )
    } else {
        "null".to_string()
    };
    let airflow_receipt = if airflow.is_some() {
        let total_iterations = coarse.airflow_iterations + fine.airflow_iterations;
        let total_solid_solves = coarse.airflow_solid_solves + fine.airflow_solid_solves;
        let total_endpoint_evaluations =
            coarse.airflow_endpoint_evaluations + fine.airflow_endpoint_evaluations;
        format!(
            "{{\"method\":\"quasi-steady-branch-iqn-ils\",\"storage\":\"solid-only\",\"max_iterations_per_step\":{max_airflow_iterations_per_step},\"max_iterations_per_natural_iteration\":{airflow_factor},\"max_solid_solves_per_step\":{max_solid_solves_per_step},\"max_krylov_iterations_per_step\":{max_krylov_per_step},\"coarse_iterations\":{},\"fine_iterations\":{},\"total_iterations\":{total_iterations},\"coarse_solid_solves\":{},\"fine_solid_solves\":{},\"total_solid_solves\":{total_solid_solves},\"max_endpoint_evaluations_per_step\":{max_airflow_endpoint_evaluations_per_step},\"total_endpoint_evaluations\":{total_endpoint_evaluations},\"maximum_coupled_energy_residual_j\":{},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{}}},\"history\":\"unchanged-physical-old-temperature-in-every-natural-air-and-radiation-trial\",\"endpoint\":\"last-counted-solid-response-with-independent-actual-air-reference-physical-residual-and-energy-gates\",\"flow\":\"retained-steady-operating-point\"}}",
            coarse.airflow_iterations, fine.airflow_iterations,
            coarse.airflow_solid_solves, fine.airflow_solid_solves,
            number(coarse.maximum_coupled_energy_residual_j
                .max(fine.maximum_coupled_energy_residual_j))?,
            number(fine.stored_j)?, number(fine.coupled_net_input_j)?,
            number(fine.stored_j - fine.coupled_net_input_j)?,
        )
    } else {
        "null".to_string()
    };
    let natural_receipt = if natural_requested {
        let total_iterations = coarse.natural_iterations + fine.natural_iterations;
        let total_solid_solves = coarse.natural_solid_solves + fine.natural_solid_solves;
        let total_endpoint_evaluations =
            coarse.natural_endpoint_evaluations + fine.natural_endpoint_evaluations;
        format!(
            "{{\"method\":\"implicit-area-mean-natural-convection\",\"max_iterations_per_step\":{natural_factor},\"max_solid_solves_per_step\":{max_solid_solves_per_step},\"max_krylov_iterations_per_step\":{max_krylov_per_step},\"max_endpoint_evaluations_per_step\":{max_endpoint_evaluations_per_step},\"coarse_iterations\":{},\"fine_iterations\":{},\"total_iterations\":{total_iterations},\"coarse_solid_solves\":{},\"fine_solid_solves\":{},\"total_solid_solves\":{total_solid_solves},\"total_endpoint_evaluations\":{total_endpoint_evaluations},\"maximum_physical_energy_residual_j\":{},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{}}},\"history\":\"unchanged-physical-old-temperature-in-every-natural-air-and-radiation-trial\",\"endpoint\":\"last-counted-solid-response-with-independent-full-physical-residual-and-energy-gates\",\"equilibrium\":\"existing-card-domain-no-zero-delta-limiting-law\"}}",
            coarse.natural_iterations, fine.natural_iterations,
            coarse.natural_solid_solves, fine.natural_solid_solves,
            number(coarse.maximum_coupled_energy_residual_j
                .max(fine.maximum_coupled_energy_residual_j))?,
            number(fine.stored_j)?, number(fine.coupled_net_input_j)?,
            number(fine.stored_j - fine.coupled_net_input_j)?,
        )
    } else {
        "null".to_string()
    };
    let workload_receipt = workload
        .as_ref()
        .map_or("null", |history| history.receipt.as_str());
    let coarse_max_step = number(TimeGrids::largest_step(&grids.coarse))?;
    let fine_max_step = number(TimeGrids::largest_step(&grids.fine))?;
    let uniform_step = |steps: usize| {
        if workload.is_none() {
            number(policy.horizon.value / steps as f64)
        } else {
            Ok("null".to_string())
        }
    };
    let sampled_peak_receipt = fine.sampled_peak.receipt(
        goal_region_name, fine_steps, policy.horizon.value,
    )?;
    let receipt = format!(
        "{{\"schema\":\"fs-cli-transient-conduction-v1\",\"method\":\"backward-euler\",\"status\":\"completed\",\"qoi_time\":\"final\",\"initial_temperature_k\":{},\"final_time_s\":{},\"coarse_steps\":{coarse_steps},\"fine_steps\":{fine_steps},\"total_steps\":{total_steps},\"max_steps\":{},\"coarse_step_s\":{},\"fine_step_s\":{},\"coarse_max_step_s\":{coarse_max_step},\"fine_max_step_s\":{fine_max_step},\"workload\":{workload_receipt},\"energy_tolerance_j\":{},\"capacities\":[{}],\"nonlinear\":{nonlinear_receipt},\"radiation\":{radiation_receipt},\"conjugate\":{airflow_receipt},\"natural\":{natural_receipt},\"coarse\":[{}],\"fine\":[{}],\"sampled_peak\":{sampled_peak_receipt},\"temporal_error\":{{\"status\":{},\"qoi\":\"temperature-max\",\"region\":{},\"coarse_final_k\":{},\"fine_final_k\":{},\"absolute_difference_k\":{},\"assumed_order\":1,\"safety_factor\":1.25,\"estimated_half_width_k\":{estimate},\"spatial_error_measured\":false}},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{},\"maximum_step_residual_j\":{}}},\"authority\":\"Estimated\",\"no_claim\":{}{adjoint_fragment}}}",
        number(initial_temperature)?,
        number(policy.horizon.value)?,
        policy.max_steps,
        uniform_step(coarse_steps)?,
        uniform_step(fine_steps)?,
        number(policy.energy_tolerance.value)?,
        capacity_rows.join(","),
        coarse.rows,
        fine.rows,
        json_string(temporal_status),
        json_string(goal_region_name),
        number(coarse_qoi)?,
        number(fine_qoi)?,
        number(difference)?,
        number(fine.stored_j)?,
        number(fine.net_input_j)?,
        number(fine.stored_j - fine.net_input_j)?,
        number(fine.maximum_energy_residual_j)?,
        json_string(NO_CLAIM)
    );
    let (conjugate_receipt, mut derived_boundary) = match (airflow, fine.airflow.as_ref()) {
        (Some(path), Some(accepted)) => {
            let derived = path.segments.iter().map(|segment| {
                let reference = accepted.applied_references.get(&segment.target)
                    .copied().ok_or_else(|| bad(
                        "the final airflow endpoint lost an applied Robin reference",
                    ))?;
                Ok((segment.target.clone(), (segment.htc_w_m2_k, reference)))
            }).collect::<Result<BTreeMap<_, _>, SolveRefusal>>()?;
            (Some(airflow::receipt_fragment(path, &accepted.outcome)?), derived)
        }
        (None, None) => (None, BTreeMap::new()),
        _ => return Err(bad("the final transient endpoint has inconsistent airflow evidence")),
    };
    let natural_receipt = match (natural_requested, fine.natural.as_ref()) {
        (true, Some(accepted)) => {
            for law in natural {
                let coefficient = accepted.applied_coefficients.get(&law.target)
                    .copied().ok_or_else(|| bad(
                        "the final natural endpoint lost an applied Robin coefficient",
                    ))?;
                derived_boundary.insert(law.target.clone(), (coefficient, law.ambient_k));
            }
            Some(accepted.receipt.clone())
        }
        (false, None) => None,
        _ => return Err(bad("the final transient endpoint has inconsistent natural evidence")),
    };
    let endpoint = fine.endpoint;
    let endpoint_dt = fine.endpoint_dt_s;
    let endpoint_storage_w = endpoint.stored_energy_change_j / endpoint_dt;
    let energy = fine.endpoint_energy;
    let endpoint_energy_residual_w = energy.energy_residual_j / endpoint_dt;
    let net_w = energy.net_input_w;
    let scale_w = endpoint.source_w.abs()
        .max(energy.dirichlet_in_w.abs())
        .max(energy.neumann_out_w.abs())
        .max(energy.robin_out_w.abs())
        .max(f64::MIN_POSITIVE);
    let mut linear_evidence = coarse.linear;
    for mut row in fine.linear {
        row.nonlinear_iteration += coarse_steps;
        linear_evidence.push(row);
    }
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
    let report = ConductionReport {
        iterations: if nonlinear.is_some() {
            nonlinear_updates
        } else if natural_requested {
            coarse.natural_solid_solves + fine.natural_solid_solves
        } else if airflow.is_some() {
            coarse.airflow_solid_solves + fine.airflow_solid_solves
        } else if radiation.is_some() {
            radiation_trials
        } else {
            total_steps
        },
        residual_history: linear_evidence
            .iter()
            .map(|row| row.true_relative_residual)
            .collect(),
        final_residual: endpoint.relative_residual,
        residual_threshold: linear.tolerance,
        stop_reason: StopReason::ResidualTolerance,
        linear: linear_evidence,
        energy: EnergyBalance {
            source_w: endpoint.source_w,
            neumann_out_w: energy.neumann_out_w,
            robin_out_w: energy.robin_out_w,
            dirichlet_in_w: energy.dirichlet_in_w,
            closure_w: net_w,
            scale_w,
        },
        material_provenance: provenance,
        material_receipts: receipts,
        element_material_identity: assignment,
        interface_fluxes: endpoint.contact_fluxes,
        robin_fluxes: endpoint.robin_fluxes,
        free_dofs: DofMap::new(
            fine.airflow.as_ref().map(|accepted| &accepted.boundary)
                .or_else(|| fine.natural.as_ref().map(|accepted| &accepted.boundary))
                .unwrap_or(problem.boundary),
            problem.mesh.vertex_count(),
        ).map_err(lower)?.n(),
        elements: problem.mesh.element_count(),
    };
    poll(cx, deadline)?;
    let conduction = ConductionSolution {
        temperature: endpoint.temperature,
        report,
    };
    let solution = match fine.radiation {
        Some(RadiativeEndpoint::Ambient {
            combined_boundary, convective_robin_fluxes, convective_out_w, report,
        }) => radiation::SolidSolution::from_radiative_endpoint(
            conduction, combined_boundary, convective_robin_fluxes, convective_out_w, report,
        ),
        Some(RadiativeEndpoint::Enclosure(endpoint)) => radiation::SolidSolution::from_enclosure_endpoint(
            conduction, enclosure.ok_or_else(|| bad("final enclosure endpoint lost its bound model"))?,
            endpoint,
        ),
        None => radiation::SolidSolution::from_conduction(conduction),
    };
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



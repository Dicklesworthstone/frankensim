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
    BackwardEuler, NonlinearStepConfig, StepConfig, StepSolution,
};
use fs_conduction::{
    ConductionError, ConductionProblem, ConductionReport, ConductionSolution, EnergyBalance,
    LinearConfig, LinearSolveEvidence, ResidualClaim, StopReason, ThermalInterfaces,
};
use fs_exec::Cx;
use fs_project::{ConductionTransient, ProjectSpec, ThermalBoundaryCondition};

use super::{
    SolveRefusal, canonical_f64, conduction_error, json_string, ladder_functional,
    ladder_target_region,
};

mod workload;
use workload::{TimeGrids, Workload};

pub(super) const NO_CLAIM: &str = "Estimated finite-mesh final-time temperature from backward Euler with declared temperature-independent regional heat capacity and the bound conductivity evaluated at each endpoint. Temperature-dependent conductivity uses residual-gated Newton/FGMRES with the full k'(T) tangent and immutable physical history. Explicit regional power histories replace delivered watts on their named volume regions; both grids land on every workload switch. The nested-grid final-maximum difference estimates temporal error at assumed order one; it is not a spatial error bound, observed-order proof, continuum or continuous-time maximum bound, validated heat capacity, or compliance certificate. Prescribed temperatures are endpoint data, including their discrete boundary-node storage reaction. No airflow storage, natural convection, radiation, latent heat or transient adjoint is inferred. Endpoint heat storage is reported separately from numerical energy closure; transient energy_residual_j is the checked storage-minus-net-input balance.";

pub(super) struct NativeTransientResult {
    pub(super) solution: ConductionSolution,
    pub(super) receipt: String,
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
        "use explicit regional capacities and admitted conductivity curves with prescribed thermal laws; inspect the time/work/energy declaration",
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

/// Refuse incompatible model combinations before meshing or solving.
pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
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
    let setup = spec
        .cooling
        .as_ref()
        .and_then(|cooling| cooling.conduction.as_ref())
        .ok_or_else(|| bad("transient conduction has no spatial declaration"))?;
    if setup.radiation.is_some()
        || setup.boundaries.iter().any(|row| {
            matches!(
                row.condition,
                ThermalBoundaryCondition::AirflowConvection { .. }
                    | ThermalBoundaryCondition::NaturalConvection { .. }
            )
        })
    {
        return Err(bad(
            "native transient cooling currently admits static Dirichlet/Neumann/Robin boundaries and matching contact; coupled airflow, natural convection and radiation need their time-dependent producer",
        ));
    }
    if spec
        .solver
        .as_ref()
        .is_some_and(|solver| matches!(solver.fidelity.as_str(), "adaptive" | "ladder"))
    {
        return Err(bad(
            "the declared transient time-work cap applies to one fixed mesh; spatial ladder/adaptive runs are not yet admitted together with it",
        ));
    }
    if spec
        .outputs
        .iter()
        .flatten()
        .any(|output| output.name.contains("adjoint"))
    {
        return Err(bad(
            "steady nominal adjoints cannot differentiate a transient trajectory",
        ));
    }
    TimeGrids::new(transient)?;
    Ok(())
}

struct GridResult {
    endpoint: StepSolution,
    endpoint_dt_s: f64,
    rows: String,
    linear: Vec<LinearSolveEvidence>,
    stored_j: f64,
    net_input_j: f64,
    maximum_energy_residual_j: f64,
    nonlinear_updates: usize,
    nonlinear_backtracks: usize,
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
    labels: &[u32],
    goal_region: Option<u32>,
    deadline: Option<(Instant, f64)>,
) -> Result<GridResult, SolveRefusal> {
    let steps = ends_s.len();
    let mut old = vec![policy.initial_temperature.value; problem.mesh.vertex_count()];
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
        let (endpoint, nonlinear_row) = if let Some(nonlinear) = nonlinear {
            let solved = engine
                .advance_nonlinear_prescribed(
                    cx,
                    step_problem,
                    interfaces,
                    &old,
                    dt,
                    config,
                    nonlinear,
                )
                .map_err(lower)?;
            nonlinear_updates += solved.nonlinear_iterations;
            nonlinear_backtracks += solved.backtracks;
            let row = format!(
                "{{\"updates\":{},\"backtracks\":{},\"initial_residual_j\":{},\"residual_j\":{},\"threshold_j\":{}}}",
                solved.nonlinear_iterations,
                solved.backtracks,
                number(solved.initial_residual_j)?,
                number(solved.residual_j)?,
                number(solved.threshold_j)?,
            );
            (solved.step, row)
        } else {
            let endpoint = engine
                .advance_prescribed(cx, step_problem, interfaces, &old, dt, config)
                .map_err(lower)?;
            (endpoint, "null".to_string())
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
        let net_w = endpoint.source_w + endpoint.dirichlet_in_w
            - endpoint.neumann_out_w
            - endpoint.robin_out_w;
        stored_j += endpoint.stored_energy_change_j;
        net_input_j += dt * net_w;
        maximum_energy_residual_j = maximum_energy_residual_j.max(endpoint.energy_residual_j.abs());
        rows.push(format!("{{\"step\":{ordinal},\"time_s\":{},\"dt_s\":{},\"final_region_max_k\":{},\"stored_energy_change_j\":{},\"net_input_w\":{},\"source_w\":{},\"energy_residual_j\":{},\"relative_residual\":{},\"krylov_iterations\":{},\"nonlinear\":{nonlinear_row}}}",
            number(end)?, number(dt)?, number(qoi)?, number(endpoint.stored_energy_change_j)?,
            number(net_w)?, number(endpoint.source_w)?, number(endpoint.energy_residual_j)?, number(endpoint.relative_residual)?,
            endpoint.krylov_iterations));
        evidence.push(LinearSolveEvidence {
            nonlinear_iteration: ordinal - 1,
            method: if nonlinear.is_some() {
                "fgmres-backward-euler-newton"
            } else {
                "pcg-backward-euler"
            },
            iterations: endpoint.krylov_iterations,
            reported: ResidualClaim::TrueEuclidean(endpoint.relative_residual),
            true_relative_residual: endpoint.relative_residual,
            converged_true: endpoint.relative_residual < linear.tolerance,
            stall: None,
        });
        old.clone_from(&endpoint.temperature);
        endpoint_dt_s = dt;
        last = Some(endpoint);
        time = end;
    }
    poll(cx, deadline)?;
    Ok(GridResult {
        endpoint: last.ok_or_else(|| bad("the transient grid has no steps"))?,
        endpoint_dt_s,
        rows: rows.join(","),
        linear: evidence,
        stored_j,
        net_input_j,
        maximum_energy_residual_j,
        nonlinear_updates,
        nonlinear_backtracks,
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
    linear: LinearConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<NativeTransientResult, SolveRefusal> {
    admit(spec)?;
    poll(cx, deadline)?;
    let policy = policy(spec).ok_or_else(|| bad("no transient declaration was supplied"))?;
    if labels.len() != problem.mesh.element_count() {
        return Err(bad(
            "transient regional capacities do not match the retained mesh labels",
        ));
    }
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
    let bytes = problem
        .mesh
        .vertex_count()
        .checked_mul(128)
        .and_then(|n| {
            problem
                .mesh
                .element_count()
                .checked_mul(4096)
                .and_then(|e| n.checked_add(e))
        })
        .and_then(|n| total_steps.checked_mul(2048).and_then(|s| n.checked_add(s)))
        .ok_or_else(|| bad("transient memory estimate overflows"))?;
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
            .check(policy.initial_temperature.value)
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
    let engine = BackwardEuler::per_element(cx, problem.mesh, &capacities).map_err(lower)?;
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
        labels,
        goal_region,
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
        labels,
        goal_region,
        deadline,
    )?;
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
            config.max_iterations,
            linear.max_iterations,
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
    let receipt = format!(
        "{{\"schema\":\"fs-cli-transient-conduction-v1\",\"method\":\"backward-euler\",\"status\":\"completed\",\"qoi_time\":\"final\",\"initial_temperature_k\":{},\"final_time_s\":{},\"coarse_steps\":{coarse_steps},\"fine_steps\":{fine_steps},\"total_steps\":{total_steps},\"max_steps\":{},\"coarse_step_s\":{},\"fine_step_s\":{},\"coarse_max_step_s\":{coarse_max_step},\"fine_max_step_s\":{fine_max_step},\"workload\":{workload_receipt},\"energy_tolerance_j\":{},\"capacities\":[{}],\"nonlinear\":{nonlinear_receipt},\"coarse\":[{}],\"fine\":[{}],\"temporal_error\":{{\"status\":{},\"qoi\":\"temperature-max\",\"region\":{},\"coarse_final_k\":{},\"fine_final_k\":{},\"absolute_difference_k\":{},\"assumed_order\":1,\"safety_factor\":1.25,\"estimated_half_width_k\":{estimate},\"spatial_error_measured\":false}},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{},\"maximum_step_residual_j\":{}}},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        number(policy.initial_temperature.value)?,
        number(policy.horizon.value)?,
        policy.max_steps,
        uniform_step(coarse_steps)?,
        uniform_step(fine_steps)?,
        number(policy.energy_tolerance.value)?,
        capacity_rows.join(","),
        coarse.rows,
        fine.rows,
        json_string(temporal_status),
        json_string(
            spec.requirements
                .iter()
                .flatten()
                .find(|row| region_ids.get(&row.region).copied() == goal_region)
                .map_or("", |row| row.region.as_str())
        ),
        number(coarse_qoi)?,
        number(fine_qoi)?,
        number(difference)?,
        number(fine.stored_j)?,
        number(fine.net_input_j)?,
        number(fine.stored_j - fine.net_input_j)?,
        number(fine.maximum_energy_residual_j)?,
        json_string(NO_CLAIM)
    );
    let endpoint = fine.endpoint;
    let endpoint_dt = fine.endpoint_dt_s;
    let endpoint_storage_w = endpoint.stored_energy_change_j / endpoint_dt;
    let endpoint_energy_residual_w = endpoint.energy_residual_j / endpoint_dt;
    let net_w =
        endpoint.source_w + endpoint.dirichlet_in_w - endpoint.neumann_out_w - endpoint.robin_out_w;
    let scale_w = endpoint
        .source_w
        .abs()
        .max(endpoint.dirichlet_in_w.abs())
        .max(endpoint.neumann_out_w.abs())
        .max(endpoint.robin_out_w.abs())
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
            neumann_out_w: endpoint.neumann_out_w,
            robin_out_w: endpoint.robin_out_w,
            dirichlet_in_w: endpoint.dirichlet_in_w,
            closure_w: net_w,
            scale_w,
        },
        material_provenance: provenance,
        material_receipts: receipts,
        element_material_identity: assignment,
        interface_fluxes: endpoint.contact_fluxes,
        robin_fluxes: endpoint.robin_fluxes,
        free_dofs: DofMap::new(problem.boundary, problem.mesh.vertex_count())
            .map_err(lower)?
            .n(),
        elements: problem.mesh.element_count(),
    };
    poll(cx, deadline)?;
    Ok(NativeTransientResult {
        solution: ConductionSolution {
            temperature: endpoint.temperature,
            report,
        },
        receipt,
        temporal_half_width_k,
        endpoint_storage_w,
        endpoint_energy_residual_w,
    })
}

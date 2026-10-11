//! Native fixed-mesh total-enthalpy trajectories through the shared P1 owner.
//!
//! Reference mass remains fixed. The accepted enthalpy, including latent heat,
//! is the only storage history; temperature is a constitutive observation.

use fs_conduction::transient::enthalpy::{
    EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError, ReferenceEnthalpyMaterial,
    },
};
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};

use super::*;

const NO_CLAIM: &str = "Estimated stationary P1 total-enthalpy transport with frozen reference mass and caller-declared equilibrium phase charts. Latent heat is stored only in specific enthalpy; no apparent heat capacity or additional latent account is used. Existing sourced conductivity, finite-resistance contact and prescribed Neumann/Robin boundaries remain active. Every workload switch is represented on both independent time grids. The final field comes only from the completed fine grid. The temperature difference is an assumed-first-order temporal estimate, not a spatial or continuous-time maximum bound; agreement on a latent plateau does not prove convergence of enthalpy or phase fraction. Supplied material-card identities and source text are attribution, not verification of the declared curves. Geometry, reference density and energetic internal variables are fixed. No material motion, density-driven expansion, phase kinetics, phase-dependent conductivity, radiation, dynamic convection or trajectory adjoint is inferred.";

/// Storage-specific scope carried by the outer conduction report.
pub(super) fn no_claim() -> &'static str {
    NO_CLAIM
}

/// Check product combinations before the staged driver starts expensive work.
pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    let Some(policy) = policy(spec).filter(|time| time.enthalpy.is_some()) else {
        return Ok(());
    };
    let setup = spec
        .cooling
        .as_ref()
        .and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("enthalpy storage has no spatial conduction declaration"))?;
    if !policy.capacities.is_empty() {
        return Err(bad(
            "enthalpy and volumetric heat-capacity storage cannot both be active",
        ));
    }
    if setup.radiation.is_some() {
        return Err(bad(
            "native enthalpy currently requires non-radiating exterior laws",
        ));
    }
    if setup.boundaries.iter().any(|row| {
        matches!(
            row.condition,
            ThermalBoundaryCondition::FixedTemperature { .. }
                | ThermalBoundaryCondition::NaturalConvection { .. }
                | ThermalBoundaryCondition::AirflowConvection { .. }
        )
    }) {
        return Err(bad(
            "native enthalpy requires prescribed Neumann or Robin boundaries; temperature alone does not select latent enthalpy, and dynamic convection is not lowered to a fixed coefficient",
        ));
    }
    if adjoint::requested(spec)?.is_some() {
        return Err(bad(
            "native temperature-history adjoints do not differentiate enthalpy storage",
        ));
    }
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

struct GridResult {
    endpoint: EnthalpyStepSolution,
    endpoint_dt_s: f64,
    threshold_j: f64,
    rows: String,
    linear: Vec<LinearSolveEvidence>,
    residuals_w: Vec<f64>,
    updates: usize,
    stored_j: f64,
    net_input_j: f64,
    maximum_energy_residual_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn march(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    initial: &[f64],
    ends_s: &[f64],
    workload: Option<&Workload<'_>>,
    config: EnthalpyStepConfig,
    labels: &[u32],
    goal_region: Option<u32>,
    deadline: Option<(Instant, f64)>,
) -> Result<GridResult, SolveRefusal> {
    let mut old = initial.to_vec();
    let mut rows = String::new();
    let mut linear = Vec::new();
    let mut residuals_w = Vec::with_capacity(ends_s.len());
    let mut last = None;
    let mut time = 0.0;
    let mut endpoint_dt_s = 0.0;
    let mut threshold_j = 0.0;
    let mut updates = 0;
    let mut stored_j = 0.0;
    let mut net_input_j = 0.0;
    let mut maximum_energy_residual_j = 0.0_f64;
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
        let endpoint = engine
            .advance(cx, step_problem, interfaces, &old, dt, config)
            .map_err(failed)?;
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
        threshold_j = config
            .newton
            .absolute_tolerance
            .max(config.newton.relative_tolerance * initial_residual.max(f64::MIN_POSITIVE));
        let mut krylov_iterations = 0;
        for row in &endpoint.newton.history {
            residuals_w.push(finite(row.residual_before / dt)?);
            linear.push(LinearSolveEvidence {
                nonlinear_iteration: updates + row.iteration,
                method: "fgmres-enthalpy-newton",
                iterations: row.linear_iterations,
                reported: ResidualClaim::TrueEuclidean(row.linear_relative_residual),
                true_relative_residual: finite(row.linear_relative_residual)?,
                converged_true: row.linear_relative_residual <= row.forcing,
                stall: None,
            });
            krylov_iterations += row.linear_iterations;
        }
        updates += endpoint.newton.iterations;
        let net_w = finite(endpoint.source_w - endpoint.neumann_out_w - endpoint.robin_out_w)?;
        stored_j = finite(stored_j + endpoint.stored_energy_change_j)?;
        net_input_j = finite(net_input_j + dt * net_w)?;
        maximum_energy_residual_j = maximum_energy_residual_j.max(endpoint.energy_residual_j.abs());
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
        if index != 0 {
            rows.push(',');
        }
        rows.push_str(&format!(
            "{{\"step\":{},\"time_s\":{},\"dt_s\":{},\"final_region_max_k\":{},\"stored_energy_change_j\":{},\"net_input_w\":{},\"source_w\":{},\"energy_residual_j\":{},\"residual_j\":{},\"threshold_j\":{},\"newton_iterations\":{},\"krylov_iterations\":{krylov_iterations},\"liquid_mass_kg\":{},\"mean_liquid_mass_fraction\":{}}}",
            index + 1, number(end)?, number(dt)?, number(qoi)?,
            number(endpoint.stored_energy_change_j)?, number(net_w)?, number(endpoint.source_w)?,
            number(endpoint.energy_residual_j)?, number(endpoint.newton.residual_norm)?,
            number(threshold_j)?, endpoint.newton.iterations,
            number(liquid_mass)?, number(liquid_mass / total_mass)?,
        ));
        old.clone_from(&endpoint.specific_enthalpy_j_kg);
        last = Some(endpoint);
        endpoint_dt_s = dt;
        time = end;
    }
    Ok(GridResult {
        endpoint: last.ok_or_else(|| bad("enthalpy trajectory has no accepted time interval"))?,
        endpoint_dt_s,
        threshold_j,
        rows,
        linear,
        residuals_w,
        updates,
        stored_j,
        net_input_j,
        maximum_energy_residual_j,
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
    let total_newton_limit = total_steps
        .checked_mul(config.max_newton_iterations)
        .ok_or_else(|| bad("combined enthalpy Newton allowance overflows"))?;
    let total_krylov_limit = total_newton_limit
        .checked_mul(restart)
        .and_then(|v| v.checked_mul(newton.max_linear_cycles))
        .ok_or_else(|| bad("combined enthalpy Krylov allowance overflows"))?;
    total_newton_limit
        .checked_mul(storage.max_backtracks as usize + 1)
        .ok_or_else(|| bad("combined enthalpy line-search allowance overflows"))?;
    let knot_count = storage.materials.iter().try_fold(0usize, |count, row| {
        count
            .checked_add(row.knots.len())
            .ok_or_else(|| bad("enthalpy knot count overflows"))
    })?;
    // Same declared-memory admission convention as the existing native time
    // path: numerical buffers plus sparse staging and retained text/telemetry.
    // It is not a measurement of allocator overhead or process peak RSS.
    let bytes = n
        .checked_mul(8 * (2 * restart + 48) + 160)
        .and_then(|v| elements.checked_mul(4096).and_then(|e| v.checked_add(e)))
        .and_then(|v| knot_count.checked_mul(256).and_then(|k| v.checked_add(k)))
        .and_then(|v| total_steps.checked_mul(2048).and_then(|r| v.checked_add(r)))
        .and_then(|v| {
            total_newton_limit
                .checked_mul(160)
                .and_then(|r| v.checked_add(r))
        })
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
        material_rows.push(format!(
            "{{\"region\":{},\"material_index\":{index},\"material_card_identity\":{},\"phase_curve_identity\":{},\"source\":{},\"reference_density_kg_m3\":{},\"initial_specific_enthalpy_j_kg\":{}}}",
            json_string(&row.region), json_string(&row.material_card_identity.to_hex()),
            json_string(&curve.identity().to_hex()), json_string(&row.source),
            number(row.reference_density.value)?, number(row.initial_specific_enthalpy.value)?,
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
    let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, problem.mesh, &materials, &assignment,
        EnthalpyBudget { max_vertices: n, max_elements: elements })
        .map_err(|error| match error {
            HeterogeneousEnthalpyError::Enthalpy(inner) => failed(inner),
            other => bad(format!("enthalpy material assignment refused: {other}; distinct regional charts require distinct interface vertices and an explicit contact law")),
        })?;
    let mut initial = Vec::with_capacity(n);
    for (vertex, &id) in engine.vertex_material_ids().iter().enumerate() {
        if vertex % 256 == 0 {
            poll(cx, deadline)?;
        }
        initial.push(storage.materials[id].initial_specific_enthalpy.value);
    }
    let goal_region = ladder_target_region(spec, region_ids);
    if goal_region.is_none() {
        return Err(bad(
            "enthalpy temporal comparison needs a temperature-maximum requirement region",
        ));
    }
    let workload = Workload::bind(cx, spec, policy, problem.mesh, labels, region_ids, deadline)?;
    let coarse = march(
        cx,
        &engine,
        problem,
        interfaces,
        &initial,
        &grids.coarse,
        workload.as_ref(),
        config,
        labels,
        goal_region,
        deadline,
    )?;
    let fine = march(
        cx,
        &engine,
        problem,
        interfaces,
        &initial,
        &grids.fine,
        workload.as_ref(),
        config,
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
    let receipt = format!(
        "{{\"schema\":\"fs-cli-transient-enthalpy-v1\",\"method\":\"backward-euler-total-enthalpy\",\"status\":\"completed\",\"qoi_time\":\"final\",\"initial_temperature_k\":{},\"final_time_s\":{},\"coarse_steps\":{},\"fine_steps\":{},\"total_steps\":{total_steps},\"max_steps\":{},\"coarse_max_step_s\":{},\"fine_max_step_s\":{},\"workload\":{workload_receipt},\"energy_tolerance_j\":{},\"materials\":[{}],\"nonlinear\":{{\"method\":\"newton-fgmres\",\"residual_units\":\"J\",\"report_residual_units\":\"W\",\"max_updates_per_step\":{},\"max_backtracks_per_update\":{},\"linear_restart\":{restart},\"max_linear_cycles\":{},\"total_newton_limit\":{total_newton_limit},\"total_krylov_limit\":{total_krylov_limit},\"total_updates\":{},\"residual_atol_j\":{},\"residual_rtol\":{},\"armijo_c\":0.0001,\"line_search_shrink\":0.5}},\"coarse\":[{}],\"fine\":[{}],\"temporal_error\":{{\"status\":{},\"qoi\":\"temperature-max\",\"coarse_final_k\":{},\"fine_final_k\":{},\"absolute_difference_k\":{},\"assumed_order\":1,\"safety_factor\":1.25,\"estimated_half_width_k\":{estimate},\"max_specific_enthalpy_difference_j_kg\":{},\"max_liquid_mass_fraction_difference\":{},\"spatial_error_measured\":false}},\"energy\":{{\"stored_change_j\":{},\"integrated_net_input_j\":{},\"window_residual_j\":{},\"maximum_step_residual_j\":{}}},\"endpoint\":{{\"vertex_order\":\"retained-conduction-mesh\",\"specific_enthalpy_j_kg\":[{}],\"liquid_mass_fraction\":[{}],\"reference_nodal_mass_kg\":[{}]}},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        number(policy.initial_temperature.value)?,
        number(policy.horizon.value)?,
        grids.coarse.len(),
        grids.fine.len(),
        policy.max_steps,
        number(TimeGrids::largest_step(&grids.coarse))?,
        number(TimeGrids::largest_step(&grids.fine))?,
        number(policy.energy_tolerance.value)?,
        material_rows.join(","),
        config.max_newton_iterations,
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
    let endpoint = fine.endpoint;
    let net_w = finite(endpoint.source_w - endpoint.neumann_out_w - endpoint.robin_out_w)?;
    let endpoint_storage_w = finite(endpoint.stored_energy_change_j / fine.endpoint_dt_s)?;
    let endpoint_energy_residual_w = finite(endpoint.energy_residual_j / fine.endpoint_dt_s)?;
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
        final_residual: finite(endpoint.newton.residual_norm / fine.endpoint_dt_s)?,
        residual_threshold: finite(fine.threshold_j / fine.endpoint_dt_s)?,
        stop_reason: StopReason::ResidualTolerance,
        linear: evidence,
        energy: EnergyBalance {
            source_w: endpoint.source_w,
            neumann_out_w: endpoint.neumann_out_w,
            robin_out_w: endpoint.robin_out_w,
            dirichlet_in_w: 0.0,
            closure_w: net_w,
            scale_w: endpoint
                .source_w
                .abs()
                .max(endpoint.neumann_out_w.abs())
                .max(endpoint.robin_out_w.abs())
                .max(f64::MIN_POSITIVE),
        },
        material_provenance: provenance,
        material_receipts: receipts,
        element_material_identity: assignment,
        interface_fluxes: endpoint.contact_fluxes,
        robin_fluxes: endpoint.robin_fluxes,
        free_dofs: n,
        elements,
    };
    poll(cx, deadline)?;
    Ok(NativeTransientResult {
        solution: radiation::SolidSolution::from_conduction(ConductionSolution {
            temperature: endpoint.temperature,
            report,
        }),
        receipt,
        conjugate_receipt: None,
        natural_receipt: None,
        derived_boundary: BTreeMap::new(),
        temporal_half_width_k,
        endpoint_storage_w,
        endpoint_energy_residual_w,
    })
}

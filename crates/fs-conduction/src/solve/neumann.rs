//! Pure-Neumann conduction with a declared volume-mean temperature.
//!
//! The heat equation fixes temperature differences; a closed body's absolute
//! temperature requires an additional datum. This lane supplies the P1 volume
//! mean without introducing a thermostat, modifying the loads, or discarding a
//! physical equilibrium equation. It is a steady boundary-value solve, not a
//! prediction of the final temperature from an initial caloric state.
//!
//! K * 1 = 0 requires sum(load) = 0. The volume weights m_i = sum(V_e/4)
//! impose sum(m_i*T_i) = volume*mean without changing that original load.
//!
//! Constant anisotropic/heterogeneous bulk and sheet conductivity are admitted.
//! Temperature-dependent laws require a nonlinear constrained solve and refuse
//! here. Disconnected bodies require separate mean data and also refuse.

use fs_sparse::Coo;

use super::{
    CgState, ConductionError, ConductionProblem, ConductionReport, ConductionSolution,
    CsrOp, DofMap, ElementMaterials, LinearConfig, LinearSolveEvidence, StopReason,
    assemble_operator_scaled_with_interfaces, full_residual, norm2, physical_report,
    residual_roundoff_floor, spd_preconditioner, true_relative_residual,
};
use fs_exec::Cx;

/// Numerical policy for a balanced pure-Neumann equilibrium.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanTemperatureConfig {
    /// Relative full-field heat residual against the original load norm.
    pub residual_rtol: f64,
    /// Absolute full-field heat residual floor, W.
    pub residual_atol_w: f64,
    /// Maximum net supplied heat, W. Incompatible loads are never projected.
    pub net_heat_tolerance_w: f64,
    /// Recomputed relative residual required of the sparse difference solve.
    pub linear_tolerance: f64,
    /// Total Krylov iterations; cancellation is checked every 32 iterations.
    pub max_iterations: usize,
}

impl Default for MeanTemperatureConfig {
    fn default() -> Self {
        let linear = LinearConfig::default();
        Self {
            residual_rtol: 1e-10,
            residual_atol_w: 1e-30,
            net_heat_tolerance_w: 1e-10,
            linear_tolerance: linear.tolerance,
            max_iterations: linear.max_iterations,
        }
    }
}

/// Solve the original spatial heat equations with a specified P1 volume mean.
///
/// The original boundary must contain only prescribed fluxes and/or adiabatic
/// faces. Both the integrated source/boundary balance and the assembled load
/// sum must meet the explicit watt tolerance. One algebraic coordinate is
/// eliminated to solve temperature DIFFERENCES, then the volume mean is set.
/// Every original heat equation, including that eliminated row, is rechecked.
/// The returned report consequently contains all original free rows and zero
/// Dirichlet heat; there is no fictitious reservoir reaction.
///
/// The supplied mean is an independent physical datum. It equals a conserved
/// energy datum only under a separately justified caloric model (for example,
/// uniform constant volumetric heat capacity). No energy-to-temperature or
/// phase-equilibrium inference is made. The mean is imposed up to floating-point
/// rounding; the residual threshold reports the existing rowwise roundoff floor.
///
/// # Errors
/// Nonfinite/nonpositive temperatures, unsupported constitutive states,
/// temperature-dependent bulk/sheet laws, disconnected/orphan vertices,
/// anchored boundaries, incompatible heat, cancellation, and unresolved
/// linear/full-field residuals refuse without changing caller-owned state.
#[allow(clippy::too_many_lines)]
pub fn solve_mean_temperature(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    mean_temperature_k: f64,
    config: MeanTemperatureConfig,
) -> Result<ConductionSolution, ConductionError> {
    poll(cx, 0)?;
    if !mean_temperature_k.is_finite() || mean_temperature_k <= 0.0 {
        return Err(invalid("mean temperature must be finite and positive"));
    }
    if !config.residual_rtol.is_finite() || config.residual_rtol <= 0.0
        || !config.residual_atol_w.is_finite() || config.residual_atol_w < 0.0
        || !config.net_heat_tolerance_w.is_finite() || config.net_heat_tolerance_w < 0.0
        || !config.linear_tolerance.is_finite() || config.linear_tolerance <= 0.0
        || config.linear_tolerance >= 1.0 || config.max_iterations == 0
    {
        return Err(invalid("finite admissible heat tolerances and a positive Krylov budget required"));
    }
    if !problem.boundary.dirichlet().is_empty() || problem.boundary.has_robin()
        || problem.boundary.face_count() != problem.mesh.boundary().len()
    {
        return Err(invalid("mean constraint requires the original unanchored Neumann boundary"));
    }
    if problem.boundary.has_temperature_dependent_surface_conduction() {
        return Err(invalid("temperature-dependent sheets need a nonlinear mean-constrained solve"));
    }
    if let Some(assigned) = problem.element_materials {
        assigned.validate_for(problem.mesh)?;
    }
    let n = problem.mesh.vertex_count();
    let mut parents: Vec<_> = (0..n).collect();
    let mut weights = vec![0.0; n];
    for (e, tet) in problem.mesh.complex().tets.iter().enumerate() {
        if e % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx, e)?; }
        let model = match problem.element_materials {
            Some(assigned) => assigned.model_for(e)?,
            None => problem.material,
        };
        if model.is_temperature_dependent() {
            return Err(invalid("temperature-dependent bulk needs a nonlinear mean-constrained solve"));
        }
        let first = root(&mut parents, tet[0] as usize);
        for &vertex in tet {
            let vertex = vertex as usize;
            let other = root(&mut parents, vertex);
            parents[other] = first;
            weights[vertex] += problem.mesh.element_volume(e) / 4.0;
        }
    }
    let connected = root(&mut parents, 0);
    for (i, &weight) in weights.iter().enumerate() {
        if i % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx, i)?; }
        if !weight.is_finite() || weight <= 0.0 || root(&mut parents, i) != connected {
            return Err(invalid("one connected body with no orphan vertices requires one mean datum"));
        }
    }
    let volume = crate::require_finite("Neumann volume", compensated(weights.iter().copied()))?;
    if volume <= 0.0 { return Err(invalid("positive finite volume required")); }
    for weight in &mut weights { *weight /= volume; }

    let dofs = DofMap::new(problem.boundary, n)?;
    // Sample each constant coefficient inside its original local support;
    // different regions need not admit a common uniform query temperature.
    let initial = super::initial_free(&problem, &dofs, &super::InitialGuess::DirichletMean)?;
    let system = assemble_operator_scaled_with_interfaces(
        cx, problem.mesh, problem.boundary, problem.material, problem.source,
        &initial, None, None, problem.element_materials,
    )?;
    let (energy, _) = super::energy_balance(
        problem.mesh, problem.boundary, problem.source, &system, &dofs, &initial,
    );
    let net_heat = crate::require_finite("Neumann assembled heat", compensated(system.load.iter().copied()))?;
    crate::require_finite("Neumann integrated heat", energy.closure_w)?;
    if net_heat.abs() > config.net_heat_tolerance_w
        || energy.closure_w.abs() > config.net_heat_tolerance_w
    {
        return Err(invalid(&format!(
            "no steady heat balance: assembled {net_heat:e} W, integrated {:e} W; tolerance {:e} W",
            energy.closure_w, config.net_heat_tolerance_w,
        )));
    }
    let load_norm = crate::require_finite("Neumann heat norm", norm2(&system.load))?;
    let initial_residual = crate::require_finite("Neumann initial residual", norm2(&full_residual(&system, &initial)))?;

    // Zero is a temperature DIFFERENCE at one vertex, never an absolute wall
    // temperature. No load is changed, including the temporarily omitted row.
    let mut coo = Coo::new(n - 1, n - 1);
    for row in 1..n {
        if row % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx, row)?; }
        let (columns, values) = system.operator.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            crate::require_finite("Neumann operator", value)?;
            if column != 0 { coo.push(row - 1, column - 1, value); }
        }
    }
    let reduced = coo.assemble();
    let rhs = &system.load[1..];
    let preconditioner = spd_preconditioner(&reduced);
    let operator = CsrOp::symmetric(reduced.clone());
    let mut cg = CgState::new(&operator, &preconditioner, rhs);
    let linear_report = loop {
        poll(cx, cg.iters)?;
        let remaining = config.max_iterations - cg.iters;
        let report = cg.run(&operator, &preconditioner, config.linear_tolerance, remaining.min(32));
        if report.converged || !report.rel_residual.is_finite()
            || cg.iters >= config.max_iterations
        {
            break report;
        }
    };
    let truth = true_relative_residual(&reduced, rhs, &cg.x);
    if !truth.is_finite() || truth >= config.linear_tolerance {
        return Err(ConductionError::LinearSolveFailed {
            iteration: 0, krylov_iterations: linear_report.iters,
            true_relative_residual: truth, tolerance: config.linear_tolerance,
        });
    }
    let mean_difference = crate::require_finite(
        "Neumann mean difference",
        compensated(weights[1..].iter().zip(&cg.x).map(|(&w, &t)| w * t)),
    )?;
    let mut temperature = Vec::with_capacity(n);
    for (i, difference) in std::iter::once(0.0).chain(cg.x.iter().copied()).enumerate() {
        if i % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx, i)?; }
        let value = mean_temperature_k + (difference - mean_difference);
        if !value.is_finite() || value <= 0.0 {
            return Err(invalid("the requested mean and heat loads yield a nonpositive/nonfinite temperature"));
        }
        temperature.push(value);
    }
    // Requery every original material/table at the actual final field. A
    // constant table still has a validity span; setting a mean cannot waive it.
    let final_system = assemble_operator_scaled_with_interfaces(
        cx, problem.mesh, problem.boundary, problem.material, problem.source,
        &temperature, None, None, problem.element_materials,
    )?;
    let physical = physical_report(problem, None, &final_system, &dofs, &temperature)?;
    let threshold = crate::require_finite("Neumann residual threshold",
        config.residual_rtol.mul_add(load_norm, config.residual_atol_w)
            .max(residual_roundoff_floor(&final_system, &dofs, &temperature)))?;
    crate::require_finite("Neumann final residual", physical.final_residual)?;
    if physical.final_residual > threshold {
        return Err(ConductionError::NotConverged {
            iterations: 1, residual: physical.final_residual, threshold,
        });
    }
    poll(cx, n)?;
    Ok(ConductionSolution {
        temperature,
        report: ConductionReport {
            iterations: 1, residual_history: vec![initial_residual],
            final_residual: physical.final_residual, residual_threshold: threshold,
            stop_reason: StopReason::ResidualTolerance,
            linear: vec![LinearSolveEvidence {
                nonlinear_iteration: 0, method: "pcg", iterations: linear_report.iters,
                reported: linear_report.residual_claim(), true_relative_residual: truth,
                converged_true: true, stall: linear_report.diagnosis,
            }],
            energy: physical.energy,
            material_provenance: problem.element_materials.map_or_else(
                || problem.material.provenance(), ElementMaterials::provenance),
            material_receipts: problem.element_materials.map_or_else(
                || problem.material.receipts().len(), |assigned| assigned.receipts().len()),
            element_material_identity: problem.element_materials.map(ElementMaterials::identity),
            interface_fluxes: physical.interface_fluxes, robin_fluxes: physical.robin_fluxes,
            free_dofs: n, elements: problem.mesh.element_count(),
        },
    })
}

fn root(parents: &mut [usize], mut index: usize) -> usize {
    while parents[index] != index {
        parents[index] = parents[parents[index]];
        index = parents[index];
    }
    index
}

fn compensated(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "Neumann equilibrium",
        what: what.into(),
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "Neumann equilibrium",
        at,
    })
}

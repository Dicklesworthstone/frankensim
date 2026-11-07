//! Nonlinear pure-Neumann equilibrium with a prescribed P1 volume mean.
//!
//! The Jacobian's LEFT null vector is constant, but its right null vector
//! generally is not: shifting a nonlinear solution changes conductivity.
//! Two solves with one reduced Jacobian recover its reference-temperature
//! response and enforce the mean without adding a source or a reservoir.

use fs_exec::Cx;
use fs_solver::nonlinear::{FgmresState, SolverRunProgress};
use fs_solver::{CsrOp, norm2};
use fs_sparse::Coo;

use super::super::{
    AssembledSystem, ConductionError, ConductionProblem, ConductionReport, ConductionSolution,
    DofMap, ElementMaterials, InitialGuess, LineSearch, LinearConfig, LinearSolveEvidence,
    Nonlinearity, SolveConfig, StopReason, StopRule, ThermalInterfaces,
    assemble_jacobian_with_optional_interfaces, assemble_operator_scaled_with_interfaces,
    full_residual, physical_report, residual_roundoff_floor, spd_preconditioner,
    true_relative_residual,
};
use super::{compatible_heat, compensated, invalid, poll, require_connected, volume_weights};

/// Numerical policy for mean-constrained Newton conduction.
#[derive(Debug, Clone, PartialEq)]
pub struct MeanTemperatureNewtonConfig {
    /// Full original heat residual stop rule, scaled by the original load norm.
    pub stop: StopRule,
    /// Total inner Krylov iterations for EACH of the two solves per Newton
    /// iteration; restart length does not enlarge this budget.
    pub linear: LinearConfig,
    /// Armijo decrease and admissibility backtracking policy.
    pub line_search: LineSearch,
    /// Maximum net supplied heat, W; incompatible loads are never projected.
    pub net_heat_tolerance_w: f64,
    /// Optional nodal initial field. Its volume mean is shifted to the supplied
    /// datum BEFORE the first material query; the resulting field must be
    /// positive and in every original material's support. None starts uniformly.
    pub initial_temperature: Option<Vec<f64>>,
}

impl Default for MeanTemperatureNewtonConfig {
    fn default() -> Self {
        Self {
            stop: StopRule::default(),
            linear: LinearConfig::default(),
            line_search: LineSearch::default(),
            net_heat_tolerance_w: 1e-10,
            initial_temperature: None,
        }
    }
}

/// Solve temperature-dependent conduction without an artificial thermostat.
///
/// All original Neumann/adiabatic heat equations, including the algebraically
/// eliminated reference row, must pass the residual test. The volume mean is
/// maintained to floating-point accuracy at every accepted iterate. The native
/// bulk and sheet tangents retain piecewise material laws and their validity
/// spans, including heterogeneous and anisotropic assignments.
///
/// For a reference vertex 0, solve J_ff*u = -R_f and J_ff*z = -J_f0. Then
/// delta = (0,u) + a*(1,z), where a enforces the mean. Both right-hand sides
/// reuse the same sparse preconditioner. No dense mean row or Jacobian inverse
/// is formed. A singular reduced tangent or mean response refuses; no
/// regularization silently changes the problem.
///
/// The mean is an independent datum, not an inference from conserved energy.
/// This finds a discrete equilibrium within declared tolerances; it does not
/// establish uniqueness or a continuum error bound.
///
/// # Errors
/// Invalid policy/geometry, incompatible heat, material-domain violations,
/// cancellation, singular tangents, exhausted budgets and unresolved original
/// heat residuals return an error without modifying caller-owned state.
#[allow(clippy::too_many_lines)]
pub fn solve_mean_temperature_newton(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    mean_temperature_k: f64,
    config: MeanTemperatureNewtonConfig,
) -> Result<ConductionSolution, ConductionError> {
    solve_mean_temperature_newton_impl(cx, problem, None, mean_temperature_k, config)
}

/// Mean-constrained nonlinear conduction through a finite-resistance assembly.
///
/// The existing matching/nonmatching contact operator enters the residual,
/// exact tangent, both constrained Newton solves, and final heat report once.
/// Contact resistance remains fixed; bulk and tangential-sheet conductivity
/// may depend on the actual local temperature through their original laws.
/// Connectivity is checked on the complete admitted Picard operator, so one
/// mean is sufficient only when the contacts join every body. No artificial
/// thermostat, load projection or welded temperature continuity is introduced.
/// The numerical policy and refusals are those of [`solve_mean_temperature_newton`].
pub fn solve_mean_temperature_newton_with_interfaces(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: &ThermalInterfaces,
    mean_temperature_k: f64,
    config: MeanTemperatureNewtonConfig,
) -> Result<ConductionSolution, ConductionError> {
    solve_mean_temperature_newton_impl(cx, problem, Some(interfaces), mean_temperature_k, config)
}

#[allow(clippy::too_many_lines)]
fn solve_mean_temperature_newton_impl(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    mean_temperature_k: f64,
    mut config: MeanTemperatureNewtonConfig,
) -> Result<ConductionSolution, ConductionError> {
    let weights = volume_weights(cx, problem, mean_temperature_k)?;
    SolveConfig {
        nonlinearity: Nonlinearity::Newton {
            line_search: config.line_search,
        },
        stop: config.stop,
        linear: config.linear,
        initial: InitialGuess::Uniform(mean_temperature_k),
    }
    .validate()?;
    if !config.stop.residual_atol.is_finite()
        || config.stop.residual_atol < 0.0
        || !config.stop.step_atol.is_finite()
        || config.stop.step_atol < 0.0
        || !config.net_heat_tolerance_w.is_finite()
        || config.net_heat_tolerance_w < 0.0
        || config.linear.tolerance >= 1.0
        || config.linear.max_iterations == 0
    {
        return Err(invalid(
            "finite nonnegative heat/step tolerances and a positive Krylov budget required",
        ));
    }
    let n = weights.len();
    let dofs = DofMap::new(problem.boundary, n)?;
    let mut temperature = config
        .initial_temperature
        .take()
        .unwrap_or_else(|| vec![mean_temperature_k; n]);
    if temperature.len() != n {
        return Err(ConductionError::FieldLength {
            field: "Neumann initial temperature",
            expected: n,
            found: temperature.len(),
        });
    }
    for &value in &temperature {
        crate::require_finite("Neumann initial temperature", value)?;
    }
    let offset = mean_error(&weights, &temperature, mean_temperature_k);
    for value in &mut temperature {
        *value -= offset;
    }
    if temperature.iter().any(|t| !t.is_finite() || *t <= 0.0) {
        return Err(invalid(
            "mean-adjusted initial temperatures must be finite and positive",
        ));
    }
    let mut system = assemble(cx, problem, interfaces, &temperature)?;
    require_connected(cx, &system.operator)?;
    compatible_heat(
        problem,
        &system,
        &dofs,
        &temperature,
        config.net_heat_tolerance_w,
    )?;
    let heat_scale = crate::require_finite("Neumann heat norm", norm2(&system.load))?;
    let mut history = Vec::new();
    let mut linear = Vec::new();
    let mut last_step = f64::INFINITY;
    for iteration in 0..=config.stop.max_iterations {
        poll(cx, iteration)?;
        let residual = full_residual(&system, &temperature);
        let rnorm = crate::require_finite("Neumann residual", norm2(&residual))?;
        let threshold = crate::require_finite(
            "Neumann residual threshold",
            config
                .stop
                .residual_rtol
                .mul_add(heat_scale, config.stop.residual_atol)
                .max(residual_roundoff_floor(&system, &dofs, &temperature)),
        )?;
        let mean_residual = mean_error(&weights, &temperature, mean_temperature_k);
        let mean_floor = 32.0
            * f64::EPSILON
            * temperature
                .iter()
                .copied()
                .fold(mean_temperature_k, f64::max);
        if rnorm <= threshold && mean_residual.abs() <= mean_floor {
            compatible_heat(
                problem,
                &system,
                &dofs,
                &temperature,
                config.net_heat_tolerance_w,
            )?;
            let physical = physical_report(problem, interfaces, &system, &dofs, &temperature)?;
            poll(cx, iteration)?;
            return Ok(ConductionSolution {
                temperature,
                report: ConductionReport {
                    iterations: iteration,
                    residual_history: history,
                    final_residual: rnorm,
                    residual_threshold: threshold,
                    stop_reason: StopReason::ResidualTolerance,
                    linear,
                    energy: physical.energy,
                    material_provenance: problem.element_materials.map_or_else(
                        || problem.material.provenance(),
                        ElementMaterials::provenance,
                    ),
                    material_receipts: problem.element_materials.map_or_else(
                        || problem.material.receipts().len(),
                        |assigned| assigned.receipts().len(),
                    ),
                    element_material_identity: problem
                        .element_materials
                        .map(ElementMaterials::identity),
                    interface_fluxes: physical.interface_fluxes,
                    robin_fluxes: physical.robin_fluxes,
                    free_dofs: n,
                    elements: problem.mesh.element_count(),
                },
            });
        }
        if iteration == config.stop.max_iterations || last_step <= config.stop.step_atol {
            return Err(ConductionError::NotConverged {
                iterations: iteration,
                residual: rnorm,
                threshold,
            });
        }
        let direction = direction(
            cx,
            problem,
            interfaces,
            &temperature,
            &residual,
            &weights,
            mean_residual,
            &config.linear,
            iteration,
            &mut linear,
        )?;
        let mut alpha = 1.0_f64;
        let mut full_step_span = None;
        let mut accepted = None;
        for backtracks in 0..=config.line_search.max_backtracks {
            poll(cx, backtracks)?;
            let candidate: Vec<_> = temperature
                .iter()
                .zip(&direction)
                .map(|(&t, &d)| alpha.mul_add(d, t))
                .collect();
            if candidate.iter().all(|t| t.is_finite() && *t > 0.0) {
                match assemble(cx, problem, interfaces, &candidate) {
                    Ok(trial_system) => {
                        let trial = norm2(&full_residual(&trial_system, &candidate));
                        if trial.is_finite()
                            && trial <= config.line_search.armijo_c.mul_add(-alpha, 1.0) * rnorm
                        {
                            accepted = Some((candidate, trial_system));
                            break;
                        }
                    }
                    Err(span @ ConductionError::OutsideTemperatureSpan { .. }) => {
                        if backtracks == 0 {
                            full_step_span = Some(span);
                        }
                    }
                    Err(other) => return Err(other),
                }
            }
            if backtracks < config.line_search.max_backtracks {
                alpha *= config.line_search.shrink;
            }
        }
        let (next, next_system) = accepted.ok_or_else(|| {
            full_step_span.unwrap_or(ConductionError::LineSearchFailed {
                iteration,
                backtracks: config.line_search.max_backtracks,
                smallest_step: alpha,
            })
        })?;
        last_step = temperature
            .iter()
            .zip(&next)
            .map(|(&a, &b)| (b - a).abs())
            .fold(0.0_f64, f64::max);
        history.push(rnorm);
        temperature = next;
        system = next_system;
    }
    unreachable!("the declared final iteration returns convergence or refusal")
}

fn assemble(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature: &[f64],
) -> Result<AssembledSystem, ConductionError> {
    assemble_operator_scaled_with_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        problem.source,
        temperature,
        None,
        interfaces,
        problem.element_materials,
    )
}

fn mean_error(weights: &[f64], temperature: &[f64], mean: f64) -> f64 {
    // Subtract before summation so the absolute Kelvin baseline does not hide
    // the small correction needed to maintain the physical mean.
    compensated(
        weights
            .iter()
            .zip(temperature)
            .map(|(&w, &t)| w * (t - mean)),
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn direction(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature: &[f64],
    residual: &[f64],
    weights: &[f64],
    mean_residual: f64,
    config: &LinearConfig,
    iteration: usize,
    evidence: &mut Vec<LinearSolveEvidence>,
) -> Result<Vec<f64>, ConductionError> {
    let jacobian = assemble_jacobian_with_optional_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        temperature,
        interfaces,
        problem.element_materials,
    )?;
    let n = temperature.len();
    let mut coo = Coo::new(n - 1, n - 1);
    let mut reference_rhs = vec![0.0; n - 1];
    for row in 1..n {
        if row % crate::assemble::ASSEMBLY_TILE == 0 {
            poll(cx, row)?;
        }
        let (columns, values) = jacobian.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            crate::require_finite("Neumann Jacobian", value)?;
            if column == 0 {
                reference_rhs[row - 1] = -value;
            } else {
                coo.push(row - 1, column - 1, value);
            }
        }
    }
    let reduced = coo.assemble();
    let preconditioner = spd_preconditioner(&reduced);
    let operator = CsrOp::general(reduced.clone());
    let mut solve = |rhs: &[f64]| -> Result<Vec<f64>, ConductionError> {
        let mut gm = FgmresState::new(rhs, config.restart.min(config.max_iterations));
        let report = loop {
            poll(cx, gm.iters)?;
            // FGMRES budgets cycles; shorten the last cycle to honor the
            // caller's total inner-iteration cap exactly.
            gm.restart = config.restart.min(config.max_iterations - gm.iters);
            let before = gm.iters;
            let run = gm.run_cancellable(&operator, &preconditioner, rhs, config.tolerance, 1, cx);
            match run.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => {
                    return Err(ConductionError::Cancelled {
                        stage: "Neumann Newton Krylov",
                        at: gm.iters,
                    });
                }
                other => {
                    return Err(invalid(&format!(
                        "Neumann Krylov callback failed: {other:?}"
                    )));
                }
            }
            if run.report.converged
                || !run.report.rel_residual.is_finite()
                || gm.iters == before
                || gm.iters >= config.max_iterations
                || run.report.diagnosis == Some(super::super::StallDiagnosis::Breakdown)
            {
                break run.report;
            }
        };
        let truth = true_relative_residual(&reduced, rhs, &gm.x);
        if !truth.is_finite() || truth >= config.tolerance {
            return Err(ConductionError::LinearSolveFailed {
                iteration,
                krylov_iterations: report.iters,
                true_relative_residual: truth,
                tolerance: config.tolerance,
            });
        }
        evidence.push(LinearSolveEvidence {
            nonlinear_iteration: iteration,
            method: "fgmres",
            iterations: report.iters,
            reported: report.residual_claim(),
            true_relative_residual: truth,
            converged_true: true,
            stall: report.diagnosis,
        });
        Ok(gm.x)
    };
    let rhs: Vec<_> = residual[1..].iter().map(|r| -r).collect();
    let particular = solve(&rhs)?;
    let response = solve(&reference_rhs)?;
    let response_mean =
        weights[0] + compensated(weights[1..].iter().zip(&response).map(|(&w, &z)| w * z));
    let response_scale = weights[0]
        + compensated(
            weights[1..]
                .iter()
                .zip(&response)
                .map(|(&w, &z)| (w * z).abs()),
        );
    if !response_mean.is_finite() || response_mean.abs() <= 64.0 * f64::EPSILON * response_scale {
        return Err(invalid(
            "the nonlinear tangent cannot resolve the prescribed mean",
        ));
    }
    let correction_mean = compensated(weights[1..].iter().zip(&particular).map(|(&w, &u)| w * u));
    let reference_step = crate::require_finite(
        "Neumann mean correction",
        -(mean_residual + correction_mean) / response_mean,
    )?;
    let mut step = Vec::with_capacity(n);
    step.push(reference_step);
    step.extend(
        particular
            .iter()
            .zip(&response)
            .map(|(&u, &z)| reference_step.mul_add(z, u)),
    );
    Ok(step)
}

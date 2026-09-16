//! Nonlinear backward Euler with immutable physical history.
//!
//! Solve F(T) = C(T - T_old) + dt (A(T)T - b) = 0. The tangent is
//! C + dt J_steady(T), including the existing element K'(T) contribution.
//! Newton trials never update T_old, and every merit evaluation reassembles
//! the actual constitutive law. A small update is NOT a convergence test.

use super::*;
use crate::assemble::{AssembledSystem, assemble_jacobian_with_optional_interfaces};
use crate::solve::LineSearch;
use fs_solver::FgmresState;

/// Explicit nonlinear work, stopping and globalization policy for one endpoint.
/// Constant heat capacity and constant contact resistance remain prerequisites.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonlinearStepConfig {
    /// Maximum accepted Newton updates, not a convergence condition.
    pub max_iterations: usize,
    /// Relative tolerance against the INITIAL free transient residual norm.
    /// The scale is fixed during this solve; it is not the absolute-temperature
    /// load C*T_old, which could conceal a poor temperature correction.
    pub residual_rtol: f64,
    /// Absolute Euclidean free-residual floor in joules.
    pub residual_atol_j: f64,
    /// Armijo policy on the actual transient residual. Invalid material-range
    /// trials consume this same budget instead of extrapolating the curve.
    pub line_search: LineSearch,
}

impl Default for NonlinearStepConfig {
    fn default() -> Self {
        Self {
            max_iterations: 32,
            residual_rtol: 1.0e-10,
            residual_atol_j: 1.0e-10,
            line_search: LineSearch::default(),
        }
    }
}

impl NonlinearStepConfig {
    /// Admit the policy before allocation or assembly.
    ///
    /// # Errors
    /// Refuses zero update budget, non-finite/negative tolerances, relative
    /// tolerance outside (0,1), or invalid Armijo/shrink factors.
    pub fn validate(&self) -> Result<(), ConductionError> {
        if self.max_iterations == 0
            || !(self.residual_rtol.is_finite() && self.residual_rtol > 0.0 && self.residual_rtol < 1.0)
            || !(self.residual_atol_j.is_finite() && self.residual_atol_j >= 0.0)
            || !(self.line_search.armijo_c.is_finite()
                && self.line_search.armijo_c > 0.0 && self.line_search.armijo_c < 1.0)
            || !(self.line_search.shrink.is_finite()
                && self.line_search.shrink > 0.0 && self.line_search.shrink < 1.0)
        {
            return Err(invalid("nonlinear endpoint requires a positive Newton budget, rtol in (0,1), nonnegative joule tolerance and Armijo/shrink in (0,1)"));
        }
        Ok(())
    }
}

/// A residual- AND energy-admitted endpoint with its nonlinear solve evidence.
#[derive(Debug, Clone)]
pub struct NonlinearStepSolution {
    /// Same endpoint, flux and discrete-energy vocabulary as the linear step.
    /// Its Krylov count is the sum over all Newton corrections; its relative
    /// residual is the worst recomputed INNER linear-solve relative residual.
    pub step: StepSolution,
    /// Accepted Newton updates (zero if the initial state already passed).
    pub nonlinear_iterations: usize,
    /// Rejected Newton trials, including material-range refusals.
    pub backtracks: usize,
    /// Norm of F(T_old) on the free degrees of freedom, joules.
    pub initial_residual_j: f64,
    /// Recomputed norm of F(T_new) on the published field, joules.
    pub residual_j: f64,
    /// Fixed atol + rtol * initial_residual_j gate, joules.
    pub threshold_j: f64,
}

struct Evaluation {
    system: AssembledSystem,
    residual: Vec<f64>,
    norm: f64,
}

impl BackwardEuler<'_> {
    /// Advance one nonlinear endpoint from immutable old temperatures.
    ///
    /// Uses the same mesh, consistent capacity, heterogeneous material laws,
    /// boundaries and contact operators as `advance`. The NEW endpoint's k(T)
    /// is evaluated on every trial; conductivity is never frozen at T_old.
    /// Prescribed temperatures must already equal history. Constant models are
    /// also admitted, making the linear path an independent comparison seam.
    ///
    /// `config.linear.max_iterations` is a TOTAL inner-iteration cap across
    /// all Newton corrections of this endpoint. Restarts are capped at 64 and
    /// shortened for the remaining budget. Cancellation is checked at Newton,
    /// line-search and Krylov-cycle boundaries, as well as in assembly.
    ///
    /// # Errors
    /// Invalid data, sampled-material extrapolation, exhausted Newton/Krylov or
    /// line-search budgets, non-finite arithmetic, failed joule balance, or
    /// cancellation return no new physical state. The caller's history is never
    /// modified, even after several accepted internal Newton trials.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_nonlinear(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: NonlinearStepConfig,
    ) -> Result<NonlinearStepSolution, ConductionError> {
        poll(cx, 0)?;
        nonlinear.validate()?;
        if config.linear.restart == 0 {
            return Err(invalid("nonlinear endpoint requires a positive FGMRES restart"));
        }
        let dofs = self.admit_step(cx, problem, old, dt_s, config)?;
        let mut temperature = old.to_vec();
        let mut ev = self.evaluate_endpoint(cx, problem, interfaces, old, &temperature, dt_s, &dofs)?;
        let initial_residual_j = ev.norm;
        let threshold_j = finite(nonlinear.residual_atol_j + nonlinear.residual_rtol * initial_residual_j)?;
        let mut updates = 0;
        let mut work = 0;
        let mut backtracks = 0_usize;
        let mut worst_linear_residual = 0.0_f64;
        loop {
            poll(cx, updates)?;
            if ev.norm <= threshold_j {
                let residual_j = ev.norm;
                let step = self.finish_step(cx, problem, interfaces, old, dt_s, config,
                    &dofs, &ev.system, temperature, worst_linear_residual, work)?;
                return Ok(NonlinearStepSolution { step, nonlinear_iterations: updates,
                    backtracks, initial_residual_j, residual_j, threshold_j });
            }
            if updates == nonlinear.max_iterations {
                return Err(ConductionError::NotConverged {
                    iterations: updates, residual: ev.norm, threshold: threshold_j,
                });
            }
            if work == config.linear.max_iterations {
                return Err(ConductionError::LinearSolveFailed {
                    iteration: updates, krylov_iterations: work,
                    true_relative_residual: 1.0, tolerance: config.linear.tolerance,
                });
            }
            let jacobian = assemble_jacobian_with_optional_interfaces(cx, self.mesh,
                problem.boundary, problem.material, &temperature, interfaces, problem.element_materials)?;
            let full = axpy_csr(&self.capacity, 1.0, &jacobian, dt_s);
            // Newton increments on fixed nodes are ZERO: no absolute lift.
            let (matrix, _) = reduce_matrix_and_lift(&full, &dofs);
            // Precondition with the SPD Picard block, not with a nonsymmetric
            // K'(T) Jacobian masquerading as an SPD operator.
            let picard = axpy_csr(&self.capacity, 1.0, &ev.system.operator, dt_s);
            let (picard, _) = reduce_matrix_and_lift(&picard, &dofs);
            let rhs: Vec<_> = ev.residual.iter().map(|&r| -r).collect();
            let remaining = LinearConfig {
                max_iterations: config.linear.max_iterations - work, ..config.linear
            };
            let (direction, relative, iters) = solve_direction(cx, &matrix, &picard, &rhs, remaining, updates)?;
            work += iters;
            worst_linear_residual = worst_linear_residual.max(relative);
            let mut alpha = 1.0;
            let mut rejected = 0;
            loop {
                poll(cx, updates)?;
                let mut trial = temperature.clone();
                for (index, &vertex) in dofs.free().iter().enumerate() {
                    if index % 512 == 0 { poll(cx, index)?; }
                    trial[vertex] = temperature[vertex] + alpha * direction[index];
                }
                let candidate = if trial.iter().all(|t| t.is_finite()) {
                    match self.evaluate_endpoint(cx, problem, interfaces, old, &trial, dt_s, &dofs) {
                        Ok(candidate) => Some(candidate),
                        Err(ConductionError::OutsideTemperatureSpan { .. }) => None,
                        Err(error) => return Err(error),
                    }
                } else { None };
                if let Some(candidate) = candidate {
                    let sufficient = (1.0 - nonlinear.line_search.armijo_c * alpha) * ev.norm;
                    if candidate.norm <= threshold_j || candidate.norm <= sufficient {
                        temperature = trial;
                        ev = candidate;
                        break;
                    }
                }
                if rejected == nonlinear.line_search.max_backtracks {
                    return Err(ConductionError::LineSearchFailed {
                        iteration: updates, backtracks: rejected, smallest_step: alpha,
                    });
                }
                rejected += 1;
                backtracks = backtracks.checked_add(1)
                    .ok_or_else(|| invalid("nonlinear backtrack count overflow"))?;
                alpha *= nonlinear.line_search.shrink;
            }
            updates += 1;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], temperature: &[f64],
        dt: f64, dofs: &DofMap,
    ) -> Result<Evaluation, ConductionError> {
        poll(cx, 0)?;
        let system = assemble_operator_scaled_with_interfaces(cx, self.mesh, problem.boundary,
            problem.material, problem.source, temperature, None, interfaces, problem.element_materials)?;
        let delta = temperature.iter().zip(old).map(|(t, old)| finite(t - old))
            .collect::<Result<Vec<_>, _>>()?;
        let mut storage = vec![0.0; temperature.len()];
        self.capacity.spmv(&delta, &mut storage);
        poll(cx, 0)?;
        let mut flux = vec![0.0; temperature.len()];
        system.operator.spmv(temperature, &mut flux);
        poll(cx, 0)?;
        let residual = dofs.free().iter().map(|&v|
            finite(storage[v] + dt * finite(flux[v] - system.load[v])?))
            .collect::<Result<Vec<_>, _>>()?;
        let norm = finite(norm2(&residual))?;
        Ok(Evaluation { system, residual, norm })
    }
}

fn solve_direction(
    cx: &Cx<'_>, matrix: &Csr, picard: &Csr, rhs: &[f64],
    config: LinearConfig, iteration: usize,
) -> Result<(Vec<f64>, f64, usize), ConductionError> {
    poll(cx, iteration)?;
    let scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 { return Ok((vec![0.0; rhs.len()], 0.0, 0)); }
    let normalized: Vec<_> = rhs.iter().map(|v| v / scale).collect();
    let op = CsrOp::general(matrix.clone());
    let pre = spd_preconditioner(picard);
    let restart = config.restart.min(rhs.len()).min(64);
    let mut state = FgmresState::new(&normalized, restart);
    while state.rel_residual() >= config.tolerance && state.iters < config.max_iterations {
        poll(cx, state.iters)?;
        let before = state.iters;
        state.restart = restart.min(config.max_iterations - before);
        state.run(&op, &pre, &normalized, config.tolerance, 1);
        if state.iters == before { break; }
    }
    poll(cx, state.iters)?;
    let mut applied = vec![0.0; rhs.len()];
    matrix.spmv(&state.x, &mut applied);
    let residual = normalized.iter().zip(applied).map(|(b, a)| finite(b - a))
        .collect::<Result<Vec<_>, _>>()?;
    let relative = finite(norm2(&residual) / norm2(&normalized))?;
    if relative >= config.tolerance {
        return Err(ConductionError::LinearSolveFailed { iteration,
            krylov_iterations: state.iters, true_relative_residual: relative,
            tolerance: config.tolerance });
    }
    let direction = state.x.iter().map(|v| finite(v * scale)).collect::<Result<_, _>>()?;
    Ok((direction, relative, state.iters))
}

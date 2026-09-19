//! Checked scalar goals of the existing tetrahedral density-Poisson operator.
//! Both solves use the existing restarted GMRES implementation. Derivatives
//! use the owner's exact element pullback, never differentiated Krylov steps.

use super::DensityPoisson;
use fs_exec::Cx;
use fs_solver::{GmresState, LinearOp};
use std::cell::Cell;

/// Explicit per-sample dimension, solve-work and residual limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DensityGoalLimits {
    /// Maximum interior state coordinates (geometry is already constructed).
    pub max_dofs: usize,
    /// Maximum tetrahedral coefficient coordinates.
    pub max_cells: usize,
    /// Shared primal + adjoint operator-call cap, including fresh residuals.
    pub max_operator_applications: usize,
    /// Maximum Arnoldi length, in 1..=64. Shortened to fit remaining work.
    pub restart: usize,
    /// Maximum restart cycles per solve; cancellation is polled between cycles.
    pub max_cycles: usize,
    /// Strict relative Euclidean residual gate on the returned primal state.
    pub primal_tolerance: f64,
    /// Strict relative Euclidean residual gate on the returned adjoint state.
    pub adjoint_tolerance: f64,
}

impl DensityGoalLimits {
    /// Validate the solve envelope before allocating state or evaluating the operator.
    pub fn validate(self, dofs: usize, cells: usize) -> Result<(), DensityGoalError> {
        if dofs == 0 || dofs > self.max_dofs || cells > self.max_cells
            || !(1..=64).contains(&self.restart) || self.max_cycles == 0
            || dofs.checked_mul(self.restart + 8).is_none()
            || !self.primal_tolerance.is_finite() || !(0.0..1.0).contains(&self.primal_tolerance)
            || self.primal_tolerance == 0.0
            || !self.adjoint_tolerance.is_finite() || !(0.0..1.0).contains(&self.adjoint_tolerance)
            || self.adjoint_tolerance == 0.0 {
            return Err(DensityGoalError::Invalid("dimensions, solve limits or tolerances"));
        }
        Ok(())
    }
}

/// No result is returned when either solve or its derivative is unavailable.
#[derive(Debug, Clone, PartialEq)]
pub enum DensityGoalError {
    /// Invalid caller configuration, shape, or unrepresentable normalization.
    Invalid(&'static str),
    /// Another complete solve cycle could not be funded.
    Budget { phase: &'static str, applications: usize, limit: usize },
    /// The actual returned-state residual, not a Krylov recurrence estimate.
    NotConverged {
        phase: &'static str,
        iterations: usize,
        relative_residual: f64,
        tolerance: f64,
        applications: usize,
    },
    /// Finite inputs led to unrepresentable solver or derivative arithmetic.
    NonFinite { quantity: &'static str, bits: u64, applications: usize },
    /// No partial goal or gradient is published. Spent calls remain attributed.
    Cancelled { applications: usize },
}
impl core::fmt::Display for DensityGoalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "density-Poisson goal refused: {self:?}")
    }
}
impl std::error::Error for DensityGoalError {}

/// One coherent primal, adjoint, scalar residual and coefficient derivative.
/// Residual tolerances are numerical solve evidence, not a continuum error bound.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityGoalEvaluation {
    /// Interior finite-element state in the existing operator's ordering.
    pub state: Vec<f64>,
    /// Interior adjoint solving K(rho)^T lambda = goal_weights.
    pub adjoint: Vec<f64>,
    /// goal_weights^T state - target, before outer objective composition.
    pub residual: f64,
    /// Total derivative of the scalar residual, one value per tetrahedron.
    pub gradient: Vec<f64>,
    /// Freshly recomputed relative Euclidean primal residual.
    pub primal_relative_residual: f64,
    /// Freshly recomputed relative Euclidean adjoint residual.
    pub adjoint_relative_residual: f64,
    /// Arnoldi iterations used for the primal.
    pub primal_iterations: usize,
    /// Arnoldi iterations used for the adjoint.
    pub adjoint_iterations: usize,
    /// Actual physical operator calls across both solves and residual checks.
    pub operator_applications: usize,
}

fn poll(cx: Option<&Cx<'_>>, applications: &Cell<usize>) -> Result<(), DensityGoalError> {
    if let Some(cx) = cx {
        cx.checkpoint().map_err(|_| DensityGoalError::Cancelled { applications: applications.get() })?;
    }
    Ok(())
}
fn finite(value: f64, quantity: &'static str, applications: &Cell<usize>) -> Result<f64, DensityGoalError> {
    if value.is_finite() { Ok(value) } else {
        Err(DensityGoalError::NonFinite { quantity, bits: value.to_bits(), applications: applications.get() })
    }
}
fn norm(values: &[f64]) -> f64 {
    let scale = values.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
    if scale == 0.0 { return 0.0; }
    let sum: f64 = values.iter().map(|x| (x / scale) * (x / scale)).sum();
    scale * fs_math::det::sqrt(sum)
}

// An immutable coefficient view. Element assembly stays in DensityPoisson;
// its public signed-density action allocates one result vector per application.
struct GoalOp<'a, 'c> {
    model: &'a DensityPoisson<'c>,
    density: &'a [f64],
    applications: &'a Cell<usize>,
}
impl LinearOp for GoalOp<'_, '_> {
    fn n(&self) -> usize { self.model.n() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        // The driver funds the complete cycle (m+2 applies) plus a fresh check
        // before entering GMRES; this counter records actual, not reserved work.
        self.applications.set(self.applications.get() + 1);
        y.copy_from_slice(&self.model.apply_density(self.density, x));
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y); // this particular P1 Poisson operator IS symmetric
    }
}

fn solve(
    op: &GoalOp<'_, '_>, rhs: &[f64], tolerance: f64, phase: &'static str,
    limits: DensityGoalLimits, transposed: bool, cx: Option<&Cx<'_>>,
) -> Result<(Vec<f64>, usize, f64), DensityGoalError> {
    poll(cx, op.applications)?;
    let scale = rhs.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
    if scale == 0.0 { return Ok((vec![0.0; rhs.len()], 0, 0.0)); }
    if limits.max_operator_applications.saturating_sub(op.applications.get()) < 4 {
        return Err(DensityGoalError::Budget { phase, applications: op.applications.get(), limit: limits.max_operator_applications });
    }
    // Normalize the RHS before the existing Euclidean Krylov reductions. The
    // accepted-state check below uses the original RHS and rescaled solution.
    let scaled: Vec<f64> = rhs.iter().map(|x| x / scale).collect();
    if rhs.iter().zip(&scaled).any(|(a, b)| *a != 0.0 && *b == 0.0) {
        return Err(DensityGoalError::Invalid("RHS normalization erased a nonzero component"));
    }
    let denominator = norm(&scaled);
    let mut state = GmresState::new(&scaled, limits.restart);
    let mut relative = 1.0;
    for _ in 0..limits.max_cycles {
        poll(cx, op.applications)?;
        let remaining = limits.max_operator_applications.saturating_sub(op.applications.get());
        if remaining < 4 {
            return Err(DensityGoalError::Budget { phase, applications: op.applications.get(), limit: limits.max_operator_applications });
        }
        state.restart = limits.restart.min(remaining - 3);
        let _report = state.run(op, &scaled, tolerance, 1, transposed);
        poll(cx, op.applications)?;
        let mut candidate = Vec::with_capacity(rhs.len());
        for (i, value) in state.x.iter().enumerate() {
            if i % 256 == 0 { poll(cx, op.applications)?; }
            candidate.push(finite(value * scale, "rescaled solved state", op.applications)?);
        }
        let mut applied = vec![0.0; rhs.len()];
        if transposed { op.apply_transpose(&candidate, &mut applied); }
        else { op.apply(&candidate, &mut applied); }
        let mut residual = Vec::with_capacity(rhs.len());
        for (i, (b, ax)) in rhs.iter().zip(applied).enumerate() {
            if i % 256 == 0 { poll(cx, op.applications)?; }
            residual.push(finite((b - ax) / scale, "scaled true residual", op.applications)?);
        }
        relative = finite(norm(&residual) / denominator, "true relative residual", op.applications)?;
        poll(cx, op.applications)?;
        if relative < tolerance { return Ok((candidate, state.iters, relative)); }
    }
    Err(DensityGoalError::NotConverged { phase, iterations: state.iters,
        relative_residual: relative, tolerance, applications: op.applications.get() })
}

impl DensityPoisson<'_> {
    /// Number of coefficient coordinates, taken from the actual tetrahedra.
    #[must_use]
    pub fn density_count(&self) -> usize { self.complex.tets.len() }

    /// Solve K(rho)u=b and K(rho)^T lambda=w, then differentiate w^T u-target.
    /// Geometry and homogeneous unit-cube boundary conditions are exactly those
    /// of this existing operator. Density, load, goal weights and target are
    /// explicit; load/weights/target are independent of density in this API.
    ///
    /// The shared operator budget includes every GMRES call and every additional
    /// accepted-state residual check. Partial samples never escape. Cx polls at
    /// restart boundaries and around the whole element pullback; individual
    /// matrix-free applications and Arnoldi cycles remain indivisible phases.
    /// This is a fixed-mesh linear goal, not shape differentiation, nonlinear
    /// physics, a certified gradient-error bound or a multigrid performance claim.
    pub fn solve_goal(
        &self, density: &[f64], load: &[f64], weights: &[f64], target: f64,
        limits: DensityGoalLimits, cx: Option<&Cx<'_>>,
    ) -> Result<DensityGoalEvaluation, DensityGoalError> {
        let applications = Cell::new(0);
        poll(cx, &applications)?;
        limits.validate(self.n(), self.density_count())?;
        if density.len() != self.density_count() || load.len() != self.n() || weights.len() != self.n()
            || !target.is_finite() {
            return Err(DensityGoalError::Invalid("density, load, goal or target shape/value"));
        }
        for (i, rho) in density.iter().enumerate() {
            if i % 256 == 0 { poll(cx, &applications)?; }
            if !rho.is_finite() || *rho <= 0.0 { return Err(DensityGoalError::Invalid("density must be finite and positive")); }
        }
        if load.iter().chain(weights).any(|v| !v.is_finite()) {
            return Err(DensityGoalError::Invalid("load and goal weights must be finite"));
        }
        let op = GoalOp { model: self, density, applications: &applications };
        let (state, primal_iterations, primal_relative_residual) = solve(
            &op, load, limits.primal_tolerance, "primal", limits, false, cx)?;
        let (adjoint, adjoint_iterations, adjoint_relative_residual) = solve(
            &op, weights, limits.adjoint_tolerance, "adjoint", limits, true, cx)?;
        let residual = finite(fs_solver::dot(weights, &state) - target, "goal residual", &applications)?;
        poll(cx, &applications)?;
        let mut gradient = self.density_pullback(&adjoint, &state);
        for (i, derivative) in gradient.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cx, &applications)?; }
            *derivative = finite(-*derivative, "implicit coefficient derivative", &applications)?;
        }
        poll(cx, &applications)?;
        Ok(DensityGoalEvaluation { state, adjoint, residual, gradient,
            primal_relative_residual, adjoint_relative_residual, primal_iterations, adjoint_iterations,
            operator_applications: applications.get() })
    }
}

//! Retained, matrix-free compliance derivatives in element-modulus coordinates.
//!
//! For each independently applied dead load, `K(m) u = f`, with homogeneous
//! supports and `K(m) = sum m_e K_e`. Thus `g_e = -u^T K_e u`. A direction
//! requires ONE tangent solve `K du = -sum v_e K_e u`, then
//! `(H v)_e = -2 du^T K_e u`. Weighted load cases are solved independently;
//! summing their forces would introduce incorrect cross-load terms.
//!
//! The element blocks, supports and CG recurrence are the existing production
//! ones. No global stiffness, Jacobian or Hessian is assembled. Derivatives
//! are with respect to positive stiffness multipliers, NOT raw/filtered SIMP
//! densities. This module does not differentiate loads, shape or supports.

use super::DensityElasticity;
use crate::control::{EvaluationStop, SolveControl};
use crate::pipeline::LoadCase;
use fs_solver::LinearOp;
use fs_sparse::precond::Precond;

/// Relative true-residual gate and maximum Krylov iterations per solve.
#[derive(Debug, Clone, Copy)]
pub struct ComplianceSolveOptions {
    pub tolerance: f64,
    pub max_iterations: usize,
}
impl Default for ComplianceSolveOptions {
    fn default() -> Self { Self { tolerance: 1.0e-10, max_iterations: 50_000 } }
}

/// No partial derivative is returned on admission, numerical or work failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ComplianceError {
    InvalidInput(&'static str),
    NonFinite(&'static str),
    Stopped(EvaluationStop),
    Residual { stage: &'static str, relative: f64, tolerance: f64 },
}
impl std::fmt::Display for ComplianceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "elastic compliance evaluation refused: {self:?}")
    }
}
impl std::error::Error for ComplianceError {}
impl From<EvaluationStop> for ComplianceError {
    fn from(value: EvaluationStop) -> Self { Self::Stopped(value) }
}

/// An immutable physical linearization. Its borrow prevents mutation of the
/// element blocks/supports; a private modulus snapshot makes later actions
/// independent of the public operator's current `moduli` field.
///
/// `SolveControl` owns lifetime work, including rejected/cancelled calls.
/// Polling occurs between element tiles, linear-solver batches, and before
/// publication. As in the existing solver, one matrix application itself is
/// not preemptible. No continuum certificate or global byte cap is implied.
pub struct ComplianceGoal<'a> {
    elasticity: &'a DensityElasticity,
    moduli: Vec<f64>,
    diagonal: Vec<f64>,
    weights: Vec<f64>,
    displacements: Vec<Vec<f64>>,
    case_compliances: Vec<f64>,
    primal_residuals: Vec<f64>,
    gradient: Vec<f64>,
    value: f64,
    options: ComplianceSolveOptions,
}
impl std::fmt::Debug for ComplianceGoal<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ComplianceGoal").field("parameters", &self.moduli.len())
            .field("loads", &self.weights.len()).field("value", &self.value).finish_non_exhaustive()
    }
}

// Derivative operators have ZERO constrained rows. The solve operator has
// identity there, exactly like DensityElasticity::apply.
struct WeightedOp<'a> { elasticity: &'a DensityElasticity, weights: &'a [f64], identity: bool }
impl LinearOp for WeightedOp<'_> {
    fn n(&self) -> usize { self.elasticity.n }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let el = self.elasticity;
        y.fill(0.0);
        for ((k, tet), &weight) in el.ke.iter().zip(&el.tets).zip(self.weights) {
            let local = gather(el, tet, x);
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if el.free[dof] {
                        let row = 3 * a + component;
                        let mut value = 0.0;
                        for (j, &input) in local.iter().enumerate() {
                            value = k[row * 12 + j].mul_add(input, value);
                        }
                        y[dof] = (weight * value).mul_add(1.0, y[dof]);
                    }
                }
            }
        }
        if self.identity {
            for (dof, value) in y.iter_mut().enumerate() {
                if !el.free[dof] { *value = x[dof]; }
            }
        }
    }
}
struct Diagonal<'a>(&'a [f64]);
impl Precond for Diagonal<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((&r, &d), z) in r.iter().zip(self.0).zip(z) { *z = r / d; }
    }
}
fn gather(el: &DensityElasticity, tet: &[u32; 4], x: &[f64]) -> [f64; 12] {
    let mut local = [0.0; 12];
    for (a, &vertex) in tet.iter().enumerate() {
        for component in 0..3 {
            let dof = 3 * vertex as usize + component;
            if el.free[dof] { local[3 * a + component] = x[dof]; }
        }
    }
    local
}
fn finite(value: f64, what: &'static str) -> Result<f64, ComplianceError> {
    if value.is_finite() { Ok(value) } else { Err(ComplianceError::NonFinite(what)) }
}
fn norm(x: &[f64]) -> f64 {
    // Scaled sum of squares avoids overflow in the squares. The strict sqrt
    // keeps this acceptance gate independent of platform libm.
    if x.iter().any(|v| !v.is_finite()) { return f64::NAN; }
    let scale = x.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 || !scale.is_finite() { return scale; }
    let sum: f64 = x.iter().map(|v| { let q = v / scale; q * q }).sum();
    scale * fs_math::det::sqrt(sum)
}

impl<'a> ComplianceGoal<'a> {
    /// Prepare all original load cases at the exact supplied modulus point.
    /// Loads on fixed dofs are reactions and do no displacement work.
    pub fn new(
        elasticity: &'a DensityElasticity, moduli: &[f64], loads: &[LoadCase<'_>],
        options: ComplianceSolveOptions, control: &mut SolveControl<'_>,
    ) -> Result<Self, ComplianceError> {
        control.checkpoint("compliance-admission")?;
        if moduli.len() != elasticity.cells() || moduli.is_empty() {
            return Err(ComplianceError::InvalidInput("one modulus per nonempty element set is required"));
        }
        if !(options.tolerance.is_finite() && options.tolerance > 0.0 && options.tolerance < 1.0) {
            return Err(ComplianceError::InvalidInput("relative tolerance must lie in (0, 1)"));
        }
        if loads.is_empty() || !loads.iter().any(|case| case.weight > 0.0) {
            return Err(ComplianceError::InvalidInput("at least one positive-weight load is required"));
        }
        for (e, &m) in moduli.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("compliance-admission")?; }
            if !(m.is_finite() && m > 0.0) {
                return Err(ComplianceError::InvalidInput("moduli must be finite and positive"));
            }
        }
        for case in loads {
            if case.force.len() != elasticity.n || !(case.weight.is_finite() && case.weight >= 0.0) {
                return Err(ComplianceError::InvalidInput("load shape or weight is invalid"));
            }
            for (dof, &force) in case.force.iter().enumerate() {
                if dof % 256 == 0 { control.checkpoint("compliance-admission")?; }
                finite(force, "load")?;
            }
        }
        let mut diagonal = vec![0.0; elasticity.n];
        for (e, ((k, tet), &m)) in elasticity.ke.iter().zip(&elasticity.tets).zip(moduli).enumerate() {
            if e % 256 == 0 { control.checkpoint("compliance-diagonal")?; }
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if elasticity.free[dof] {
                        diagonal[dof] = m.mul_add(k[13 * (3 * a + component)], diagonal[dof]);
                    }
                }
            }
        }
        for (dof, d) in diagonal.iter_mut().enumerate() {
            if !elasticity.free[dof] { *d = 1.0; }
            if !(d.is_finite() && *d > 0.0) {
                return Err(ComplianceError::InvalidInput("stiffness diagonal must be finite and positive"));
            }
        }
        let mut goal = Self { elasticity, moduli: moduli.to_vec(), diagonal,
            weights: loads.iter().map(|case| case.weight).collect(), displacements: Vec::new(),
            case_compliances: Vec::new(), primal_residuals: Vec::new(),
            gradient: vec![0.0; moduli.len()], value: 0.0, options };
        for case in loads {
            let rhs: Vec<f64> = case.force.iter().zip(&elasticity.free)
                .map(|(&f, &free)| if free { f } else { 0.0 }).collect();
            let (u, residual) = goal.solve(&rhs, control, "compliance-primal")?;
            let value = finite(rhs.iter().zip(&u).map(|(f, u)| f * u).sum(), "compliance")?;
            goal.value = finite(goal.value + case.weight * value, "weighted compliance")?;
            let energy = goal.bilinear(&u, &u, control)?;
            for (g, e) in goal.gradient.iter_mut().zip(energy) {
                *g = finite(*g - case.weight * e, "compliance gradient")?;
            }
            goal.case_compliances.push(value);
            goal.primal_residuals.push(residual);
            goal.displacements.push(u);
        }
        control.checkpoint("compliance-publish")?;
        Ok(goal)
    }

    pub fn value(&self) -> f64 { self.value }
    pub fn gradient(&self) -> &[f64] { &self.gradient }
    pub fn displacements(&self) -> &[Vec<f64>] { &self.displacements }
    pub fn case_compliances(&self) -> &[f64] { &self.case_compliances }
    pub fn primal_relative_residuals(&self) -> &[f64] { &self.primal_residuals }

    fn bilinear(&self, a: &[f64], b: &[f64], control: &mut SolveControl<'_>)
        -> Result<Vec<f64>, ComplianceError>
    {
        let el = self.elasticity;
        let mut values = Vec::with_capacity(el.cells());
        for (e, (k, tet)) in el.ke.iter().zip(&el.tets).enumerate() {
            if e % 256 == 0 { control.checkpoint("compliance-contraction")?; }
            let (a, b) = (gather(el, tet, a), gather(el, tet, b));
            let mut value = 0.0;
            for i in 0..12 { for j in 0..12 { value += a[i] * k[12 * i + j] * b[j]; } }
            values.push(finite(value, "element bilinear form")?);
        }
        Ok(values)
    }

    fn solve(&self, rhs: &[f64], control: &mut SolveControl<'_>, stage: &'static str)
        -> Result<(Vec<f64>, f64), ComplianceError>
    {
        let op = WeightedOp { elasticity: self.elasticity, weights: &self.moduli, identity: true };
        let u = control.solve_preconditioned(&op, &Diagonal(&self.diagonal), rhs,
            self.options.tolerance, self.options.max_iterations, stage)?;
        control.checkpoint("compliance-true-residual")?;
        let mut residual = vec![0.0; rhs.len()];
        op.apply(&u, &mut residual);
        for (r, &f) in residual.iter_mut().zip(rhs) { *r = f - *r; }
        let relative = finite(norm(&residual) / norm(rhs).max(f64::MIN_POSITIVE), "true residual")?;
        if relative >= self.options.tolerance {
            return Err(ComplianceError::Residual { stage, relative, tolerance: self.options.tolerance });
        }
        control.checkpoint("compliance-solve-publish")?;
        Ok((u, relative))
    }

    /// One tangent solve per positive-weight load; zero directions cost none.
    /// The same original operator, preconditioner and true-residual gate apply.
    pub fn hessian_vector_product(&self, direction: &[f64], control: &mut SolveControl<'_>)
        -> Result<Vec<f64>, ComplianceError>
    {
        control.checkpoint("compliance-hessian")?;
        if direction.len() != self.moduli.len() {
            return Err(ComplianceError::InvalidInput("Hessian direction shape mismatch"));
        }
        for (e, &v) in direction.iter().enumerate() {
            if e % 256 == 0 { control.checkpoint("compliance-hessian")?; }
            finite(v, "Hessian direction")?;
        }
        let mut product = vec![0.0; direction.len()];
        if direction.iter().all(|&v| v == 0.0) { return Ok(product); }
        let derivative = WeightedOp { elasticity: self.elasticity, weights: direction, identity: false };
        for (&weight, u) in self.weights.iter().zip(&self.displacements) {
            if weight == 0.0 { continue; }
            control.checkpoint("compliance-tangent-rhs")?;
            let mut rhs = vec![0.0; self.elasticity.n];
            derivative.apply(u, &mut rhs);
            for value in &mut rhs { *value = -*value; }
            let (du, _) = self.solve(&rhs, control, "compliance-tangent")?;
            for (value, cross) in product.iter_mut().zip(self.bilinear(&du, u, control)?) {
                *value = finite(*value - 2.0 * weight * cross, "compliance Hessian action")?;
            }
        }
        control.checkpoint("compliance-hessian-publish")?;
        Ok(product)
    }
}

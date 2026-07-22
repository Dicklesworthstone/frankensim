//! The adjoint hook: `dJ/dρ` for a LINEAR conduction solve whose design
//! parameters are per-element conductivity multipliers.
//!
//! # What this establishes
//!
//! The design map is `K_e ← ρ_e · K` on the conduction block only.
//! With `R(T; ρ) = A(ρ)T − b(ρ) = 0` on the free dofs (`b` carries the
//! Dirichlet lift, which is itself `ρ`-dependent — hence the pullback
//! uses the FULL nodal temperature, not the free part), the implicit
//! function theorem gives
//!
//! ```text
//!   dJ/dρ_e = ∂J/∂ρ_e − λᵀ (∂R/∂ρ_e),   (∂R/∂T)ᵀ λ = ∂J/∂T
//!   ∂R/∂ρ_e = K_e · T_full |_free       (exact: A is LINEAR in ρ)
//! ```
//!
//! The adjoint solve runs through [`fs_adjoint::ift_gradient_matfree`],
//! i.e. ONE transposed solve on `fs-solver`'s stack — never a
//! differentiated Krylov iteration.
//!
//! # What this does NOT establish
//!
//! - Nothing about the NONLINEAR `k(T)` case. The construction refuses
//!   a temperature-dependent model rather than silently linearizing:
//!   the correct Jacobian there is [`crate::assemble_jacobian`]'s
//!   nonsymmetric operator, and wiring the nonlinear IFT path is a
//!   separate piece of work with its own evidence.
//! - Nothing about SHAPE derivatives. `ρ` is a coefficient, not a
//!   geometry; mesh-motion sensitivity is `fs-adjoint`'s Hadamard path
//!   and is not wired here.
//! - Nothing about goal-oriented ERROR. A gradient is not a DWR
//!   estimate; `fs-adjoint`'s `dwr-accept` feature is the place that
//!   claim lives, and this crate does not make it.

use fs_adjoint::{AdjointReport, ift_gradient_matfree};
use fs_exec::Cx;
use fs_solver::{CsrOp, norm2};

use crate::ConductionError;
use crate::assemble::{DofMap, assemble_operator_scaled, element_stiffness, reduce};
use crate::solve::{ConductionProblem, LinearConfig};

/// A conduction problem parameterized by per-element conductivity
/// multipliers, with an adjoint gradient.
pub struct ConductivityDesign<'m> {
    problem: ConductionProblem<'m>,
    dofs: DofMap,
    linear: LinearConfig,
    base_tensor: [[f64; 3]; 3],
}

/// A primal solve at one design point.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignSolution {
    /// Free-dof temperature, K.
    pub free_temperature: Vec<f64>,
    /// Full nodal temperature, K.
    pub temperature: Vec<f64>,
    /// Recomputed `‖b − Ax‖₂/‖b‖₂` of the primal solve.
    pub primal_relative_residual: f64,
    /// Krylov iterations spent on the primal solve.
    pub primal_iterations: usize,
}

impl<'m> ConductivityDesign<'m> {
    /// Bind a design problem.
    ///
    /// # Errors
    /// [`ConductionError::Conductivity`] when the material is
    /// temperature dependent (this hook is the LINEAR case only);
    /// [`ConductionError::NoFreeDofs`] / [`ConductionError::SingularPureNeumann`]
    /// from the degree-of-freedom map.
    pub fn new(
        problem: ConductionProblem<'m>,
        linear: LinearConfig,
    ) -> Result<ConductivityDesign<'m>, ConductionError> {
        if problem.material.is_temperature_dependent() {
            return Err(ConductionError::Conductivity {
                what: "the IFT hook in this crate covers the LINEAR case only; a k(T) model \
                       needs the nonsymmetric Newton Jacobian and its own verification"
                    .to_string(),
            });
        }
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        let base_tensor = problem.material.tensor_at(0.0)?;
        Ok(ConductivityDesign {
            problem,
            dofs,
            linear,
            base_tensor,
        })
    }

    /// One design parameter per element.
    #[must_use]
    pub fn parameter_count(&self) -> usize {
        self.problem.mesh.element_count()
    }

    /// The degree-of-freedom map.
    #[must_use]
    pub const fn dofs(&self) -> &DofMap {
        &self.dofs
    }

    fn check_rho(&self, rho: &[f64]) -> Result<(), ConductionError> {
        if rho.len() != self.parameter_count() {
            return Err(ConductionError::FieldLength {
                field: "conductivity design vector",
                expected: self.parameter_count(),
                found: rho.len(),
            });
        }
        for &r in rho {
            if !(r.is_finite() && r > 0.0) {
                return Err(ConductionError::Conductivity {
                    what: format!("design multiplier {r} must be finite and positive"),
                });
            }
        }
        Ok(())
    }

    /// Solve the primal problem at a design point.
    ///
    /// # Errors
    /// Assembly and linear-solve refusals.
    pub fn solve(&self, cx: &Cx<'_>, rho: &[f64]) -> Result<DesignSolution, ConductionError> {
        self.check_rho(rho)?;
        let (a_ff, b_f) = self.system(cx, rho)?;
        let op = CsrOp::symmetric(a_ff.clone());
        let precond = crate::solve::spd_preconditioner(&a_ff);
        let mut cg = fs_solver::CgState::new(&op, &precond, &b_f);
        let report = cg.run(
            &op,
            &precond,
            self.linear.tolerance,
            self.linear.max_iterations,
        );
        let (x, iters) = (cg.x, report.iters);
        let mut ax = vec![0.0f64; b_f.len()];
        a_ff.spmv(&x, &mut ax);
        let r: Vec<f64> = b_f.iter().zip(&ax).map(|(b, a)| b - a).collect();
        let denom = norm2(&b_f).max(f64::MIN_POSITIVE);
        let primal_relative_residual = norm2(&r) / denom;
        if primal_relative_residual >= self.linear.tolerance {
            return Err(ConductionError::LinearSolveFailed {
                iteration: 0,
                krylov_iterations: iters,
                true_relative_residual: primal_relative_residual,
                tolerance: self.linear.tolerance,
            });
        }
        Ok(DesignSolution {
            temperature: self.dofs.scatter(&x),
            free_temperature: x,
            primal_relative_residual,
            primal_iterations: iters,
        })
    }

    fn system(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
    ) -> Result<(fs_sparse::Csr, Vec<f64>), ConductionError> {
        let zero = vec![0.0f64; self.problem.mesh.vertex_count()];
        let system = assemble_operator_scaled(
            cx,
            self.problem.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            &zero,
            Some(rho),
        )?;
        Ok(reduce(&system, &self.dofs))
    }

    /// The linear quantity of interest `J(ρ) = Σ_i w_i T_i` over FREE
    /// dofs.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] for a mis-sized weight vector;
    /// every refusal [`ConductivityDesign::solve`] can produce.
    pub fn objective(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
        weights: &[f64],
    ) -> Result<f64, ConductionError> {
        if weights.len() != self.dofs.n() {
            return Err(ConductionError::FieldLength {
                field: "objective weights",
                expected: self.dofs.n(),
                found: weights.len(),
            });
        }
        let solution = self.solve(cx, rho)?;
        Ok(solution
            .free_temperature
            .iter()
            .zip(weights)
            .map(|(t, w)| t * w)
            .sum())
    }

    /// `dJ/dρ` by the IFT adjoint, plus the adjoint-solve report.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`] for a mis-sized weight vector;
    /// every refusal [`ConductivityDesign::solve`] can produce.
    pub fn gradient(
        &self,
        cx: &Cx<'_>,
        rho: &[f64],
        weights: &[f64],
    ) -> Result<(Vec<f64>, AdjointReport), ConductionError> {
        if weights.len() != self.dofs.n() {
            return Err(ConductionError::FieldLength {
                field: "objective weights",
                expected: self.dofs.n(),
                found: weights.len(),
            });
        }
        let solution = self.solve(cx, rho)?;
        let (a_ff, _) = self.system(cx, rho)?;
        let op = CsrOp::symmetric(a_ff);
        let temperature = solution.temperature.clone();
        let pullback =
            |lambda: &[f64]| -> Vec<f64> { self.parameter_pullback(lambda, &temperature) };
        let (gradient, report) = ift_gradient_matfree(
            &op,
            weights,
            &[],
            &pullback,
            self.linear.tolerance,
            self.linear
                .max_iterations
                .div_ceil(self.linear.restart.max(1))
                .max(1),
        );
        Ok((gradient, report))
    }

    /// `(∂R/∂ρ)ᵀ λ`, element by element: `λᵀ (K_e · T_full)|_free`.
    /// Exact, because the assembled operator is LINEAR in `ρ`.
    #[must_use]
    pub fn parameter_pullback(&self, lambda: &[f64], temperature: &[f64]) -> Vec<f64> {
        let mesh = self.problem.mesh;
        let mut out = vec![0.0f64; mesh.element_count()];
        for (e, slot) in out.iter_mut().enumerate() {
            let tet = mesh.complex().tets[e];
            let ke = element_stiffness(mesh, e, &self.base_tensor);
            let mut acc = 0.0f64;
            for a in 0..4 {
                let ia = self.dofs.slot_of(tet[a] as usize);
                let Some(ia) = ia else { continue };
                let mut row = 0.0f64;
                for b in 0..4 {
                    row = ke[a][b].mul_add(temperature[tet[b] as usize], row);
                }
                acc = lambda[ia].mul_add(row, acc);
            }
            *slot = acc;
        }
        out
    }
}

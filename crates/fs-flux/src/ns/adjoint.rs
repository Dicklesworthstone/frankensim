//! Checked discrete adjoints of the assembled velocity/pressure equations.
//!
//! Solve A^T lambda = dJ/dx using the same dense/sparse solver infrastructure
//! as the primal. The actual transpose includes strong boundary rows, the
//! pressure datum and nonsymmetric convection; it is never replaced by A.

use super::{FluxParams, FluxSystem};
use fs_sparse::{Csr, ops::transpose};

/// Converged multipliers for an assembled discrete objective gradient.
#[derive(Clone, Debug)]
pub struct FluxAdjointSolution {
    /// Multipliers in the original equation order: velocity, then pressure.
    pub multipliers: Vec<f64>,
    /// Krylov iterations; zero for dense factorization or an exact zero load.
    pub iters: usize,
    /// Independently recomputed Euclidean residual of A^T lambda = dJ/dx.
    pub rel_residual: f64,
    /// Independently recomputed residual with each transposed equation divided
    /// by its largest absolute matrix entry.
    pub scaled_rel_residual: f64,
}

/// Invalid input, a detected singular factorization, or an exhausted solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FluxAdjointError {
    /// Shape, coefficient, objective gradient or tolerance is invalid.
    Invalid(&'static str),
    /// The small direct solve detected a singular operator.
    SingularSystem,
    /// The requested tolerance was not attained; no multipliers are published.
    NonConvergence {
        /// Actual Krylov iterations consumed.
        iterations: usize,
        /// Original relative residual of the transposed equations.
        relative_residual: f64,
        /// Row-equilibrated relative residual of the transposed equations.
        scaled_relative_residual: f64,
    },
}

impl core::fmt::Display for FluxAdjointError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "flow adjoint: {reason}"),
            Self::SingularSystem => write!(f, "flow adjoint: singular operator"),
            Self::NonConvergence {
                iterations,
                relative_residual,
                scaled_relative_residual,
            } => write!(
                f,
                "flow adjoint: unconverged after {iterations} iterations (original residual {relative_residual}, scaled residual {scaled_relative_residual})"
            ),
        }
    }
}

impl std::error::Error for FluxAdjointError {}

impl FluxSystem<'_> {
    /// Solve the actual transposed assembled system with checked convergence.
    ///
    /// Sparse counting-sort transposition costs O(nnz + n); above 1500 unknowns
    /// the existing row-equilibrated block FGMRES builds its factors from this
    /// transpose. No dense matrix is created on that path. Small systems reuse
    /// dense LU with residual refinement. Both transposed residual criteria
    /// must meet `params.tol`, and `params.max_iters` bounds sparse iterations
    /// exactly, including zero and partial-restart budgets. Other fields of
    /// `params` do not modify the already assembled operator.
    ///
    /// For A(theta) x(theta) = b(theta), the returned multipliers give
    /// dJ/dtheta = partial_theta J + lambda^T (b_theta - A_theta x).
    /// Differentiate every parameter-dependent load and boundary term as well
    /// as the viscosity operator. Mesh and equation ordering must agree.
    ///
    /// This solves the supplied linearization. A nonlinear rheology or
    /// Navier--Stokes sensitivity requires its FULL residual Jacobian, including
    /// constitutive and convecting-velocity derivatives; a frozen Picard matrix
    /// only yields the sensitivity of that frozen linear problem. Residual
    /// convergence alone is not a condition-number or gradient-error bound.
    ///
    /// # Errors
    /// Refuses invalid inputs, detected singularity and failed convergence.
    pub fn solve_adjoint_checked(
        &self,
        matrix: &Csr,
        objective_gradient: &[f64],
        params: FluxParams,
    ) -> Result<FluxAdjointSolution, FluxAdjointError> {
        if matrix.nrows() != self.n
            || matrix.ncols() != self.n
            || objective_gradient.len() != self.n
        {
            return Err(FluxAdjointError::Invalid("dimension mismatch"));
        }
        if !params.tol.is_finite() || params.tol <= 0.0 {
            return Err(FluxAdjointError::Invalid(
                "tolerance must be finite and positive",
            ));
        }
        if objective_gradient.iter().any(|value| !value.is_finite())
            || (0..self.n).any(|row| matrix.row(row).1.iter().any(|value| !value.is_finite()))
        {
            return Err(FluxAdjointError::Invalid(
                "operator and objective gradient must be finite",
            ));
        }
        let transposed = transpose(matrix);
        let solution = self
            .solve_linear_result(&transposed, objective_gradient, params)
            .map_err(|_| FluxAdjointError::SingularSystem)?;
        if !solution.converged(params.tol) || solution.x.iter().any(|value| !value.is_finite()) {
            return Err(FluxAdjointError::NonConvergence {
                iterations: solution.iters,
                relative_residual: solution.rel_residual,
                scaled_relative_residual: solution.scaled_rel_residual,
            });
        }
        Ok(FluxAdjointSolution {
            multipliers: solution.x,
            iters: solution.iters,
            rel_residual: solution.rel_residual,
            scaled_rel_residual: solution.scaled_rel_residual,
        })
    }
}

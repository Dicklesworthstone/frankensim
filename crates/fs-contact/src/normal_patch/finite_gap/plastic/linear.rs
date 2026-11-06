//! Bounded equality solves for the pressure-saturation active set.
//!
//! The iterative lane applies the same symmetric free principal matrix as the
//! dense factorization, with scalar diagonal normalization. It stores the
//! translation-invariant kernel and linear-sized work vectors, never the dense
//! free matrix. A true displacement residual gates every returned free solve.

use core::cell::RefCell;

use fs_exec::Cx;
use fs_la::factor::cholesky_cancellable;
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, checked_cg, norm2};
use fs_sparse::precond::IdentityPrecond;

use super::{
    ComplianceKernel, FiniteGapPlasticError, FiniteGapPlasticModel,
    MAX_FINITE_GAP_PLASTIC_FREE_CELLS, checkpoint, factor_error, invalid,
};

const TILE: usize = 64;

/// Equality solver for the free pressures of the shared box active set.
///
/// These controls concern numerical work and acceptance, not a different
/// contact or hardness law. Iterative residuals are explicitly recomputed
/// floating-point residuals; they are not forward-error certificates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FiniteGapPlasticLinearSolver {
    /// Reference Cholesky solve, bounded by [`MAX_FINITE_GAP_PLASTIC_FREE_CELLS`].
    #[default]
    Dense,
    /// Matrix-free CG over the rectangular half-space compliance.
    ///
    /// This removes the dense free-cell storage limit. Each operator application
    /// still performs quadratic pair interactions; no runtime speedup is claimed.
    MatrixFree {
        /// Positive total CG iteration budget for each free equality solve.
        /// All defect corrections consume this same budget.
        maximum_iterations: usize,
        /// Maximum additional defect solves after the initial CG attempt.
        /// Zero disables correction but preserves recomputed-residual checks.
        maximum_corrections: usize,
    },
}

impl FiniteGapPlasticLinearSolver {
    pub(super) fn validate(self) -> Result<(), FiniteGapPlasticError> {
        if matches!(
            self,
            Self::MatrixFree {
                maximum_iterations: 0,
                ..
            }
        ) {
            return Err(invalid("positive matrix-free iteration budget"));
        }
        Ok(())
    }
}

impl FiniteGapPlasticModel {
    pub(super) fn solve_free(
        &self,
        free: &[usize],
        upper: &[usize],
        pressure: &[f64],
        target: &[f64],
        rank_one: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteGapPlasticError> {
        if self.linear_solver == FiniteGapPlasticLinearSolver::Dense
            && free.len() > MAX_FINITE_GAP_PLASTIC_FREE_CELLS
        {
            return Err(FiniteGapPlasticError::FreeSetTooLarge {
                free_cells: free.len(),
                maximum: MAX_FINITE_GAP_PLASTIC_FREE_CELLS,
            });
        }
        let upper_spring_opening = if rank_one == 0.0 {
            0.0
        } else {
            rank_one * upper.iter().map(|&i| pressure[i]).sum::<f64>()
        };
        if !upper_spring_opening.is_finite() {
            return Err(invalid("finite constrained spring opening"));
        }
        let mut rhs = Vec::with_capacity(free.len());
        for &i in free {
            checkpoint(cx)?;
            let mut value = target[i];
            for chunk in upper.chunks(TILE) {
                checkpoint(cx)?;
                for &j in chunk {
                    value = (-self.compliance.compliance(i, j)).mul_add(pressure[j], value);
                }
            }
            if rank_one != 0.0 {
                value -= upper_spring_opening;
            }
            if !value.is_finite() {
                return Err(invalid("finite constrained pressure right-hand side"));
            }
            rhs.push(value);
        }
        match self.linear_solver {
            FiniteGapPlasticLinearSolver::Dense => {
                let n = free.len();
                let mut matrix = vec![0.0; n * n];
                for (row, &i) in free.iter().enumerate() {
                    checkpoint(cx)?;
                    for (col, &j) in free.iter().enumerate().take(row + 1) {
                        let base = self.compliance.compliance(i, j);
                        let value = if rank_one == 0.0 {
                            base
                        } else {
                            base + rank_one
                        };
                        if !value.is_finite() {
                            return Err(invalid("finite free compliance"));
                        }
                        matrix[row * n + col] = value;
                        matrix[col * n + row] = value;
                    }
                }
                let factor = cholesky_cancellable(&matrix, n, || cx.checkpoint().is_err())
                    .map_err(factor_error)?;
                factor
                    .solve_cancellable(&mut rhs, || cx.checkpoint().is_err())
                    .map_err(factor_error)?;
                if rhs.iter().any(|p| !p.is_finite()) {
                    return Err(invalid("finite pressure solution"));
                }
                Ok(rhs)
            }
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations,
                maximum_corrections,
            } => self.solve_free_matrix_free(
                free,
                &rhs,
                rank_one,
                maximum_iterations,
                maximum_corrections,
                cx,
            ),
        }
    }

    fn solve_free_matrix_free(
        &self,
        free: &[usize],
        rhs: &[f64],
        rank_one: f64,
        maximum_iterations: usize,
        maximum_corrections: usize,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteGapPlasticError> {
        // Reserve displacement tolerance for rescaling, full-gradient work,
        // and the final physical/history check. Do not substitute a generic
        // relative Krylov tolerance for the caller's absolute contact budget.
        let allowance = 0.25 * self.config.complementarity_tolerance_m;
        if allowance <= 0.0 {
            return Err(invalid("representable matrix-free displacement allowance"));
        }
        let rhs_scale = rhs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        if rhs_scale <= allowance {
            return Ok(vec![0.0; free.len()]);
        }
        let rhs_shape: Vec<f64> = rhs.iter().map(|value| value / rhs_scale).collect();
        let relative_tolerance = ((allowance / rhs_scale) / norm2(&rhs_shape)).min(0.5);
        if !relative_tolerance.is_finite() || relative_tolerance <= 0.0 {
            return Err(invalid("representable matrix-free relative tolerance"));
        }
        let diagonal = self.compliance.compliance(0, 0) + rank_one;
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(invalid("finite positive free compliance diagonal"));
        }
        // Solve (C / d) q = rhs, q = d p. Normalizing the kernel once avoids
        // division inside every pair interaction, without storing a free matrix.
        let mut normalized = self.compliance.clone();
        for chunk in normalized.values.chunks_mut(TILE) {
            checkpoint(cx)?;
            for value in chunk {
                *value /= diagonal;
                if !value.is_finite() {
                    return Err(invalid("finite normalized free compliance"));
                }
            }
        }
        let operator = PrincipalCompliance {
            kernel: normalized,
            free,
            rank_one: rank_one / diagonal,
            cx,
            fault: RefCell::new(None),
        };
        let solved = checked_cg(
            &operator,
            &IdentityPrecond,
            rhs,
            CheckedCgConfig {
                tolerance: relative_tolerance,
                max_iterations: maximum_iterations,
                max_corrections: maximum_corrections,
            },
            |_| {
                operator.check_fault()?;
                checkpoint(cx)
            },
        )
        .map_err(|error| match error {
            CheckedCgError::Interrupted(error) => error,
            CheckedCgError::InvalidInput(_) => invalid("checked-CG free solve inputs"),
        })?;
        operator.check_fault()?;
        let pressure: Vec<f64> = solved.x.iter().map(|value| value / diagonal).collect();
        if pressure.iter().any(|value| !value.is_finite()) {
            return Err(invalid("finite rescaled free pressure"));
        }
        let residual_m = self.free_residual(free, &pressure, rhs, rank_one, cx)?;
        if !solved.report.converged_euclidean() || residual_m > allowance {
            return Err(FiniteGapPlasticError::LinearSolveDidNotConverge {
                free_cells: free.len(),
                iterations: solved.report.iters,
                residual_m,
            });
        }
        Ok(pressure)
    }

    fn free_residual(
        &self,
        free: &[usize],
        pressure: &[f64],
        rhs: &[f64],
        rank_one: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, FiniteGapPlasticError> {
        let common = if rank_one == 0.0 {
            0.0
        } else {
            rank_one * pressure.iter().sum::<f64>()
        };
        if !common.is_finite() {
            return Err(invalid("finite free spring displacement"));
        }
        let mut residual = 0.0_f64;
        for (row, &i) in free.iter().enumerate() {
            checkpoint(cx)?;
            let mut opening = 0.0;
            for (chunk_index, chunk) in free.chunks(TILE).enumerate() {
                checkpoint(cx)?;
                for (within_chunk, &j) in chunk.iter().enumerate() {
                    let col = chunk_index * TILE + within_chunk;
                    let value = self.compliance.compliance(i.max(j), i.min(j));
                    opening = value.mul_add(pressure[col], opening);
                }
            }
            opening += common;
            let defect = opening - rhs[row];
            if !defect.is_finite() {
                return Err(invalid("finite physical free residual"));
            }
            residual = residual.max(defect.abs());
        }
        Ok(residual)
    }
}

struct PrincipalCompliance<'a, 'scope> {
    kernel: ComplianceKernel,
    free: &'a [usize],
    rank_one: f64,
    cx: &'a Cx<'scope>,
    // LinearOp cannot return errors. Preserve the first original refusal for
    // the checked-CG callback, which runs before any iterate is accepted.
    fault: RefCell<Option<FiniteGapPlasticError>>,
}

impl PrincipalCompliance<'_, '_> {
    fn check_fault(&self) -> Result<(), FiniteGapPlasticError> {
        match self.fault.borrow().as_ref() {
            Some(error) => Err(error.clone()),
            None => Ok(()),
        }
    }

    fn try_apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), FiniteGapPlasticError> {
        checkpoint(self.cx)?;
        if x.len() != self.free.len() || y.len() != self.free.len() {
            return Err(invalid("free operator dimensions"));
        }
        let mut common = 0.0;
        if self.rank_one != 0.0 {
            for chunk in x.chunks(TILE) {
                checkpoint(self.cx)?;
                for &value in chunk {
                    common = self.rank_one.mul_add(value, common);
                }
            }
            if !common.is_finite() {
                return Err(invalid("finite normalized spring displacement"));
            }
        }
        y.fill(0.0);
        for (row_tile, rows) in self.free.chunks(TILE).enumerate() {
            for (col_tile, columns) in self.free.chunks(TILE).enumerate() {
                checkpoint(self.cx)?;
                for (row_offset, &i) in rows.iter().enumerate() {
                    let row = row_tile * TILE + row_offset;
                    for (col_offset, &j) in columns.iter().enumerate() {
                        let col = col_tile * TILE + col_offset;
                        // Exactly the lower-triangle mirroring used by dense
                        // assembly, despite signed-offset rounding differences.
                        let value = self.kernel.compliance(i.max(j), i.min(j));
                        y[row] = value.mul_add(x[col], y[row]);
                    }
                }
            }
            for row in &mut y[row_tile * TILE..row_tile * TILE + rows.len()] {
                *row += common;
                if !row.is_finite() {
                    return Err(invalid("finite normalized free opening"));
                }
            }
        }
        Ok(())
    }
}

impl LinearOp for PrincipalCompliance<'_, '_> {
    fn n(&self) -> usize {
        self.free.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if self.fault.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        if let Err(error) = self.try_apply(x, y) {
            *self.fault.borrow_mut() = Some(error);
            y.fill(f64::NAN);
        }
    }
}

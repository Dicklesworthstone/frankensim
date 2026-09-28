//! Small-system fallback when contact/Robin coefficients defeat comparison
//! dominance. Krylov columns propose R; only an outward I-A R defect proves it.

use super::{
    ConductionError, Cx, LinearConfig, LinearGoalAnalyzer, bounded_solve, invalid, map_enclosure,
    poll,
};

impl LinearGoalAnalyzer<'_> {
    /// Prepare interval-LDL inverse evidence for the stored solid matrix.
    ///
    /// The declared stability budget is the one authority for extra inverse
    /// work: a zero budget performs none, and interval updates are capped at
    /// the remaining iterations times `nnz(A)` (one iteration is one matrix
    /// pass), then charged back to the shared counter. Stored entries use the
    /// same memory-derived class as the residual limits. Only an established
    /// bound is retained, and nothing is solved.
    pub(super) fn prepare_sparse_inverse(&mut self, cx: &Cx<'_>) -> Result<(), ConductionError> {
        use fs_solver::goal::inverse::sparse::{SparseInverseLimits, VerifiedSparseInverse};
        poll(cx, self.stability_iterations)?;
        let remaining = self.config.max_stability_iterations
            .saturating_sub(self.stability_iterations);
        let pass = self.response.matrix.nnz().max(1);
        if remaining == 0 {
            return Ok(());
        }
        let limits = self.config.residual_limits;
        let sparse = SparseInverseLimits {
            max_rows: limits.max_rows,
            max_input_nonzeros: limits.max_nonzeros,
            max_entries: limits.max_nonzeros,
            max_updates: remaining.saturating_mul(pass),
        };
        let evidence = VerifiedSparseInverse::prepare(&self.response.matrix, sparse,
            || cx.checkpoint().is_ok()).map_err(map_enclosure)?;
        let paid = evidence.work().updates.div_ceil(pass).max(1);
        self.stability_iterations = self.stability_iterations.saturating_add(paid)
            .min(self.config.max_stability_iterations);
        if evidence.inverse_upper().is_some() {
            self.sparse_inverse = Some(evidence);
        }
        poll(cx, self.stability_iterations)?;
        Ok(())
    }

    pub(super) fn prepare_inverse_columns(
        &mut self,
        cx: &Cx<'_>,
        temperature: &[f64],
    ) -> Result<(), ConductionError> {
        poll(cx, self.stability_iterations)?;
        let matrix = &self.response.matrix;
        let n = matrix.nrows();
        let remaining = self
            .config
            .max_stability_iterations
            .saturating_sub(self.stability_iterations);
        // Explicit small dimension, dense storage and repeated sparse traversal
        // caps are checked before allocation. The existing residual entry cap
        // also bounds n*(n+nnz) verification visits, not merely the sparse input.
        let admitted = n > 0
            && n <= 256
            && n <= self.config.residual_limits.max_rows
            && remaining >= n
            && n.checked_mul(n)
                .is_some_and(|entries| entries <= self.config.residual_limits.max_nonzeros)
            && n.checked_add(matrix.nnz())
                .and_then(|visits| n.checked_mul(visits))
                .is_some_and(|visits| visits <= self.config.residual_limits.max_nonzeros);
        if !admitted {
            return Ok(());
        }
        let mut columns = Vec::new();
        columns
            .try_reserve_exact(n)
            .map_err(|_| invalid("inverse proposal allocation refused"))?;
        for column in 0..n {
            poll(cx, self.stability_iterations)?;
            let remaining = self
                .config
                .max_stability_iterations
                .saturating_sub(self.stability_iterations);
            if remaining == 0 {
                return Ok(());
            }
            let mut rhs = vec![0.0; n];
            rhs[column] = 1.0;
            let (values, _, iterations) = bounded_solve(
                cx,
                matrix,
                &rhs,
                LinearConfig {
                    tolerance: self.response.linear.tolerance.min(1e-10),
                    max_iterations: remaining,
                    restart: self.response.linear.restart,
                },
                false,
            )?;
            self.stability_iterations += iterations;
            if !values.iter().all(|value| value.is_finite()) {
                return Ok(());
            }
            columns.push(values);
        }
        let free_temperature = self.response.dofs.gather(temperature);
        let checked = fs_solver::goal::inverse::enclose_goal_error_with_inverse(
            matrix,
            &self.rhs,
            &free_temperature,
            &self.weights,
            &self.free_dual,
            self.stability_scaling.as_deref(),
            &columns,
            self.config.residual_limits,
            || cx.checkpoint().is_ok(),
        )
        .map_err(map_enclosure)?;
        // Retain only independently verified columns; analyze() repeats the
        // defect check against the exact immutable matrix on every field.
        if checked.inverse_infinity_upper().is_some() {
            self.inverse_columns = Some(columns);
        }
        poll(cx, self.stability_iterations)?;
        Ok(())
    }
}

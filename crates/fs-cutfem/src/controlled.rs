//! Cooperatively interruptible solves of the canonical elasticity operator.
//!
//! Assembly, boundary conditions, Jacobi scaling and the true-Euclidean
//! correction gate are the existing CutFEM mathematics. Only scheduling changes:
//! the caller can stop between bounded batches of CG iterations. An interruption
//! never returns a field masquerading as a converged solution. Assembly and one
//! matrix/vector operation are not preemptible; this is a work bound, not a
//! wall-clock latency guarantee.

use std::collections::BTreeMap;
use std::ops::ControlFlow;

use crate::elastic::CutElasticityOperator;
use crate::fem::{JacobiPrecond, recomputed_euclidean_residual_claim};
use crate::{CutFemError, NodeKey};
use fs_solver::krylov::{CgState, ResidualClaim};
use fs_solver::op::LinearOp;

/// Converged coefficients and nodal field from a controlled elasticity solve.
///
/// Construction is private: only an explicitly recomputed Euclidean residual
/// below the requested tolerance admits this value. Interrupted or exhausted
/// solves cannot construct one.
#[derive(Debug, Clone)]
pub struct ControlledElasticitySolution {
    coefficients: Vec<f64>,
    nodal: BTreeMap<NodeKey, [f64; 2]>,
    compliance: f64,
    iters: usize,
    residual_claim: ResidualClaim,
}

impl ControlledElasticitySolution {
    /// Displacement coefficients in the canonical operator's terminal order.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Reconstructed nodal displacements, including hanging-node constraints.
    #[must_use]
    pub fn nodal(&self) -> &BTreeMap<NodeKey, [f64; 2]> {
        &self.nodal
    }

    /// The admitted operator's assembled-load functional, `b^T x`.
    #[must_use]
    pub const fn compliance(&self) -> f64 {
        self.compliance
    }

    /// Aggregate CG iterations across all residual-correction passes.
    #[must_use]
    pub const fn iters(&self) -> usize {
        self.iters
    }

    /// Explicitly recomputed Euclidean residual provenance, never a recurrence.
    #[must_use]
    pub const fn residual_claim(&self) -> ResidualClaim {
        self.residual_claim
    }
}

impl CutElasticityOperator {
    /// Solve the admitted operator with caller-controlled interruption.
    ///
    /// `poll_iters` is a strictly positive bound on additional CG iterations
    /// between calls to `control`. The callback receives the cumulative count
    /// across all correction passes. `Break(reason)` returns that exact reason
    /// without publishing a solution; `Continue(())` permits more work. Checks
    /// also bracket setup, explicit residual verification, and publication.
    /// A callback must not mutate the operator or problem through side channels.
    ///
    /// An always-continuing callback follows the same arithmetic and correction
    /// order as the ordinary elasticity solve, independent of batch size. The
    /// global `max_iters` budget never resets at a batch or correction boundary.
    /// CG currently copies its diagnostic history at each batch boundary, so
    /// larger batches reduce that overhead at the cost of coarser interruption.
    ///
    /// # Errors
    /// Refuses invalid controls, numerical failure, or exhaustion before the
    /// *true* Euclidean residual meets `tol`. Cancellation is `Ok(Break(_))`,
    /// distinct from a failed PDE solve. This method does not assemble a mesh or
    /// promise interruption inside assembly, a sparse apply, or a reduction.
    pub fn solve_controlled<B>(
        &self,
        tol: f64,
        max_iters: usize,
        poll_iters: usize,
        mut control: impl FnMut(usize) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, ControlledElasticitySolution>, CutFemError> {
        if !(tol.is_finite() && tol > 0.0) || max_iters == 0 || poll_iters == 0 {
            return Err(CutFemError::InvalidElasticityInput {
                what: "controlled solve requires a finite positive tolerance and positive total/poll iteration budgets".into(),
            });
        }
        if let ControlFlow::Break(reason) = control(0) {
            return Ok(ControlFlow::Break(reason));
        }
        let preconditioner = JacobiPrecond::new(self.matrix());
        let rhs = self.rhs();
        let mut x = vec![0.0; rhs.len()];
        let mut total_iters = 0usize;

        loop {
            if let ControlFlow::Break(reason) = control(total_iters) {
                return Ok(ControlFlow::Break(reason));
            }
            let claim = recomputed_euclidean_residual_claim(self, &x, rhs);
            let residual = claim.euclidean().expect("explicit Euclidean residual");
            if let ControlFlow::Break(reason) = control(total_iters) {
                return Ok(ControlFlow::Break(reason));
            }
            if !residual.is_finite() {
                return Err(CutFemError::SolveNotConverged {
                    iters: total_iters,
                    rel_residual: residual,
                });
            }
            if residual < tol {
                let compliance = self.algebraic_compliance(&x)?;
                let nodal = self.nodal_values(&x);
                if let ControlFlow::Break(reason) = control(total_iters) {
                    return Ok(ControlFlow::Break(reason));
                }
                return Ok(ControlFlow::Continue(ControlledElasticitySolution {
                    coefficients: x,
                    nodal,
                    compliance,
                    iters: total_iters,
                    residual_claim: claim,
                }));
            }
            if total_iters >= max_iters {
                return Err(CutFemError::SolveNotConverged {
                    iters: total_iters,
                    rel_residual: residual,
                });
            }

            let mut applied = vec![0.0; rhs.len()];
            self.apply(&x, &mut applied);
            let correction_rhs: Vec<_> = rhs.iter().zip(applied).map(|(b, ax)| b - ax).collect();
            let mut correction = CgState::new(self, &preconditioner, &correction_rhs);
            let remaining = max_iters - total_iters;
            while correction.iters < remaining && !(correction.rel_residual() < tol) {
                if let ControlFlow::Break(reason) = control(total_iters + correction.iters) {
                    return Ok(ControlFlow::Break(reason));
                }
                let count = poll_iters.min(remaining - correction.iters);
                let before = correction.iters;
                let _ = correction.run(self, &preconditioner, tol, count);
                if let ControlFlow::Break(reason) = control(total_iters + correction.iters) {
                    return Ok(ControlFlow::Break(reason));
                }
                if correction.iters == before || !correction.rel_residual().is_finite() {
                    break;
                }
            }
            let completed = correction.iters;
            total_iters += completed;
            for (value, delta) in x.iter_mut().zip(correction.x) {
                *value += delta;
            }
            if completed == 0 {
                // A correction that performs no work cannot improve the
                // previously unmet explicit residual gate.
                let claim = recomputed_euclidean_residual_claim(self, &x, rhs);
                if let ControlFlow::Break(reason) = control(total_iters) {
                    return Ok(ControlFlow::Break(reason));
                }
                return Err(CutFemError::SolveNotConverged {
                    iters: total_iters,
                    rel_residual: claim.euclidean().expect("explicit Euclidean residual"),
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CutElasticity, CutStabilizationScaling, HalfPlane, Quadtree};
    use fs_material::IsotropicElastic;

    fn fixture(run: impl FnOnce(&CutElasticity<'_>)) {
        let grid = Quadtree::uniform(3);
        let sdf = HalfPlane { normal: [1.0, 0.0], offset: 2.0 };
        let material = IsotropicElastic::new(1.0, 0.3, 1.0).expect("material");
        let clamp = |x: f64, _: f64| x < 1e-9;
        run(&CutElasticity {
            grid: &grid,
            sdf: &sdf,
            material: &material,
            nitsche_beta: 20.0,
            ghost_gamma: 0.5,
            stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
            quad_depth: 2,
            clamp: Some(&clamp),
            boundary_traction: None,
            traction_free_interface: true,
            solver_tol: 1e-10,
            solver_max_iters: 2_000,
        });
    }

    #[test]
    fn batch_sizes_reproduce_the_canonical_solve_bits() {
        fixture(|solver| {
            let body = |_: f64, _: f64| [0.0, -1.0];
            let zero = |_: f64, _: f64| [0.0, 0.0];
            let reference = solver.solve(&body, &zero).expect("canonical solve");
            let operator = solver.assemble(&body, &zero).expect("canonical assembly");
            for batch in [1, 7, 64, usize::MAX] {
                let ControlFlow::Continue(solution) = operator.solve_controlled(
                    solver.solver_tol, solver.solver_max_iters, batch,
                    |_| ControlFlow::<()>::Continue(()),
                ).expect("controlled solve") else { panic!("unexpected interruption") };
                assert_eq!(
                    solution.coefficients().iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    reference.coefficients().iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                assert_eq!(solution.nodal(), reference.nodal());
                assert_eq!(solution.compliance().to_bits(), reference.compliance().to_bits());
                assert_eq!(solution.iters(), reference.iters);
                assert_eq!(solution.residual_claim(), reference.residual_claim());
            }
        });
    }

    #[test]
    fn interruption_inside_cg_returns_the_reason_not_a_solution() {
        fixture(|solver| {
            let operator = solver.assemble(&|_, _| [0.0, -1.0], &|_, _| [0.0, 0.0])
                .expect("assembly");
            let outcome = operator.solve_controlled(1e-10, 2_000, 2, |iters| {
                if iters >= 4 { ControlFlow::Break(("wall-budget", iters)) }
                else { ControlFlow::Continue(()) }
            }).expect("interrupted, not a PDE failure");
            assert!(matches!(outcome, ControlFlow::Break(("wall-budget", 4))));
            assert!(matches!(
                operator.solve_controlled(1e-10, 2_000, 2, |_| ControlFlow::Break("cancelled")),
                Ok(ControlFlow::Break("cancelled"))
            ));
        });
    }

    #[test]
    fn total_iteration_budget_does_not_reset_at_poll_boundaries() {
        fixture(|solver| {
            let operator = solver.assemble(&|_, _| [0.0, -1.0], &|_, _| [0.0, 0.0])
                .expect("assembly");
            let result = operator.solve_controlled(1e-14, 3, 1, |_| ControlFlow::<()>::Continue(()));
            assert!(matches!(result, Err(CutFemError::SolveNotConverged { iters: 3, .. })));
        });
    }

    #[test]
    fn controls_are_admitted_before_polling_and_zero_load_needs_no_cg() {
        fixture(|solver| {
            let operator = solver.assemble(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
                .expect("assembly");
            for (tol, total, batch) in [(f64::NAN, 2, 1), (0.0, 2, 1), (1e-10, 0, 1), (1e-10, 2, 0)] {
                let result = operator.solve_controlled(tol, total, batch, |_| -> ControlFlow<()> {
                    panic!("invalid controls must not start work")
                });
                assert!(matches!(result, Err(CutFemError::InvalidElasticityInput { .. })));
            }
            let ControlFlow::Continue(solution) = operator.solve_controlled(
                1e-10, 2, 1, |_| ControlFlow::<()>::Continue(()),
            ).expect("zero-load solve") else { panic!("unexpected interruption") };
            assert_eq!(solution.iters(), 0);
            assert_eq!(solution.compliance(), 0.0);
            assert_eq!(solution.residual_claim().euclidean(), Some(0.0));
        });
    }
}

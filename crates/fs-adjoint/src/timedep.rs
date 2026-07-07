//! Time-dependent adjoints: reverse sweeps under fs-ad's binomial
//! REVOLVE checkpointing (O(log N) memory instead of storing every
//! state — long horizons degrade to recomputation, not OOM; the
//! Ledger-backed SPILL tier is recorded follow-up scope). The
//! reference problem is backward-Euler heat conduction
//! (M + h·K)·u_{n+1} = M·u_n + h·b: each REVERSE step is a TRANSPOSED
//! solve through the same factorization the forward sweep used — the
//! fs-solver adjoint contract exercised through time.

use fs_ad::revolve::{checkpointed_adjoint, min_budget};
use fs_solver::{CgState, CsrOp};
use fs_sparse::Csr;
use fs_sparse::precond::IdentityPrecond;

/// Backward-Euler heat problem on assembled (M, K): the forward
/// stepper, the terminal-misfit objective, and its adjoint gradient
/// with respect to the INITIAL CONDITION.
pub struct HeatAdjoint {
    sys: CsrOp,
    mass: Csr,
    /// Time step.
    pub h: f64,
    /// Step count.
    pub steps: usize,
}

impl HeatAdjoint {
    /// Build from mass and stiffness (interior-reduced, SPD).
    #[must_use]
    pub fn new(mass: Csr, stiffness: &Csr, h: f64, steps: usize) -> HeatAdjoint {
        let n = mass.nrows();
        // A = M + h·K.
        let mut coo = fs_sparse::Coo::new(n, n);
        for r in 0..n {
            let (cols, vals) = mass.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, v);
            }
            let (cols, vals) = stiffness.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, h * v);
            }
        }
        HeatAdjoint {
            sys: CsrOp::symmetric(coo.assemble()),
            mass,
            h,
            steps,
        }
    }

    /// One forward step: solve (M + hK)·u⁺ = M·u.
    fn step_forward(&self, u: &[f64]) -> Vec<f64> {
        let n = u.len();
        let mut rhs = vec![0.0f64; n];
        self.mass.spmv(u, &mut rhs);
        let mut st = CgState::new(&self.sys, &IdentityPrecond, &rhs);
        let rep = st.run(&self.sys, &IdentityPrecond, 1e-13, 10_000);
        assert!(rep.converged, "forward heat step failed: {rep:?}");
        st.x
    }

    /// One reverse (adjoint) step: μ = M·(M + hK)⁻ᵀ·λ — the transpose
    /// of the forward map u ↦ (M + hK)⁻¹·M·u (symmetric factors, so
    /// the transposed solve reuses the same operator).
    fn step_reverse(&self, lambda: &[f64]) -> Vec<f64> {
        let n = lambda.len();
        let mut st = CgState::new(&self.sys, &IdentityPrecond, lambda);
        let rep = st.run(&self.sys, &IdentityPrecond, 1e-13, 10_000);
        assert!(rep.converged, "reverse heat step failed: {rep:?}");
        let mut out = vec![0.0f64; n];
        self.mass.spmv(&st.x, &mut out);
        out
    }

    /// Run the forward sweep to the terminal state.
    #[must_use]
    pub fn forward(&self, u0: &[f64]) -> Vec<f64> {
        let mut u = u0.to_vec();
        for _ in 0..self.steps {
            u = self.step_forward(&u);
        }
        u
    }
}

/// Gradient of the terminal misfit J = ½‖u_N − target‖² with respect
/// to u₀, by the revolve-checkpointed reverse sweep. Returns
/// (gradient, total forward evaluations) — the recompute count is the
/// price of O(log N) memory and is reported, not hidden.
#[must_use]
pub fn heat_initial_gradient(problem: &HeatAdjoint, u0: &[f64], target: &[f64]) -> (Vec<f64>, u64) {
    let steps = problem.steps;
    let budget = min_budget(steps);
    let forward = |_i: usize, u: &Vec<f64>| -> Vec<f64> { problem.step_forward(u) };
    let reverse =
        |_i: usize, _u: &Vec<f64>, bar: Vec<f64>| -> Vec<f64> { problem.step_reverse(&bar) };
    // Terminal cotangent: ∂J/∂u_N = u_N − target.
    let u_n = problem.forward(u0);
    let seed: Vec<f64> = u_n.iter().zip(target).map(|(a, b)| a - b).collect();
    let (bar, stats) = checkpointed_adjoint(&u0.to_vec(), steps, budget, &forward, &reverse, seed);
    (bar, stats.forward_steps)
}

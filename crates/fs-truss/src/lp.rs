//! The plastic-design layout LP with a first-order primal-dual solver.
//!
//! Variables: per member, split tension/compression forces
//! `q⁺, q⁻ ≥ 0`. Objective: material volume `Σ lᵢ(qᵢ⁺ + qᵢ⁻)/σ_y`.
//! Constraints: nodal equilibrium `B·(q⁺ − q⁻) = f` on FREE degrees of
//! freedom (supports drop their rows). Standard form:
//! `min cᵀx  s.t.  A x = b, x ≥ 0`.
//!
//! PDHG (Chambolle–Pock): `x ← Π₊(x − τ(c + Aᵀy))`,
//! `y ← y + σ(A(2x − x_prev) − b)`, with `τσ‖A‖² < 1` from a
//! power-iteration norm estimate. Sparse-matvec dominated (fs-sparse
//! CSR), bitwise deterministic, warm-startable across load cases. The
//! DUALITY GAP and KKT residuals are tracked every check interval —
//! the LP's own certificate of near-optimality (under this saddle the
//! dual objective is `−bᵀy`, feasible where `c + Aᵀy ≥ 0`; a uniform
//! shrink of y restores feasibility with the bound still certified).

use crate::ground::GroundStructure;
use fs_sparse::{Coo, Csr};
use std::fmt::Write as _;

/// PDHG controls.
#[derive(Debug, Clone, Copy)]
pub struct PdhgSettings {
    /// Iteration cap.
    pub max_iters: usize,
    /// Relative duality-gap target.
    pub gap_tol: f64,
    /// Check/ledger interval.
    pub check_every: usize,
}

impl Default for PdhgSettings {
    fn default() -> Self {
        PdhgSettings {
            max_iters: 200_000,
            gap_tol: 1e-6,
            check_every: 500,
        }
    }
}

/// Solve evidence.
#[derive(Debug, Clone, Default)]
pub struct PdhgReport {
    /// Iterations run.
    pub iters: usize,
    /// Final primal objective (volume).
    pub volume: f64,
    /// Final certified relative duality gap.
    pub gap: f64,
    /// Final equilibrium residual ‖Ax − b‖/‖b‖.
    pub eq_residual: f64,
    /// Gap trace (iteration, gap) at check intervals.
    pub trace: Vec<(usize, f64)>,
}

impl PdhgReport {
    /// Ledger row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"iters\":{},\"volume\":{:.8e},\"gap\":{:.3e},\"eq_residual\":{:.3e}}}",
            self.iters, self.volume, self.gap, self.eq_residual
        );
        s
    }
}

/// The assembled layout LP for one ground structure.
pub struct LayoutLp {
    /// Equilibrium matrix on free DOFs over split variables (n_free ×
    /// 2·members): columns `[q⁺ | q⁻]` with `B` and `−B` blocks.
    pub a: Csr,
    /// Aᵀ (materialized once; PDHG applies both directions).
    pub at: Csr,
    /// Cost per split variable (length/σ_y).
    pub c: Vec<f64>,
    /// Free-DOF load vector.
    pub b: Vec<f64>,
    /// Free-DOF index per (node, component); None = supported.
    pub dof_map: Vec<Option<usize>>,
    /// Estimated operator norm ‖A‖.
    pub norm_est: f64,
}

impl LayoutLp {
    /// Assemble from a ground structure, support predicate, nodal
    /// loads, and yield stress.
    #[must_use]
    pub fn assemble(
        gs: &GroundStructure,
        supported: &dyn Fn(usize, usize) -> bool,
        loads: &dyn Fn(usize) -> [f64; 2],
        sigma_y: f64,
    ) -> LayoutLp {
        let n = gs.nodes.len();
        let mut dof_map: Vec<Option<usize>> = Vec::with_capacity(2 * n);
        let mut nf = 0usize;
        for node in 0..n {
            for comp in 0..2 {
                if supported(node, comp) {
                    dof_map.push(None);
                } else {
                    dof_map.push(Some(nf));
                    nf += 1;
                }
            }
        }
        let m = gs.members.len();
        let mut coo = Coo::new(nf, 2 * m);
        for (k, &(a, b)) in gs.members.iter().enumerate() {
            let dx = (gs.nodes[b][0] - gs.nodes[a][0]) / gs.lengths[k];
            let dy = (gs.nodes[b][1] - gs.nodes[a][1]) / gs.lengths[k];
            // Unit tension in member k pulls node a toward b and b
            // toward a.
            let entries = [(2 * a, dx), (2 * a + 1, dy), (2 * b, -dx), (2 * b + 1, -dy)];
            for (dof, v) in entries {
                if let Some(row) = dof_map[dof] {
                    coo.push(row, k, v); // q⁺ column
                    coo.push(row, m + k, -v); // q⁻ column
                }
            }
        }
        let a_mat = coo.assemble();
        let at = fs_sparse::ops::transpose(&a_mat);
        let mut b_vec = vec![0.0f64; nf];
        for node in 0..n {
            let f = loads(node);
            for comp in 0..2 {
                if let Some(row) = dof_map[2 * node + comp] {
                    b_vec[row] = f[comp];
                }
            }
        }
        let mut c = Vec::with_capacity(2 * m);
        for &l in &gs.lengths {
            c.push(l / sigma_y);
        }
        for &l in &gs.lengths {
            c.push(l / sigma_y);
        }
        // Power iteration for ‖A‖ (deterministic start).
        let mut v: Vec<f64> = (0..2 * m).map(|i| 1.0 + ((i % 7) as f64) * 0.1).collect();
        let mut norm_est = 1.0;
        let mut av = vec![0.0f64; nf];
        for _ in 0..30 {
            a_mat.spmv(&v, &mut av);
            let mut atv = vec![0.0f64; 2 * m];
            at.spmv(&av, &mut atv);
            let nrm = atv.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
            norm_est = nrm.sqrt();
            for (vi, ai) in v.iter_mut().zip(&atv) {
                *vi = ai / nrm;
            }
        }
        LayoutLp {
            a: a_mat,
            at,
            c,
            b: b_vec,
            dof_map,
            norm_est,
        }
    }

    /// Run PDHG from a warm start (zeros for cold); returns the
    /// primal solution (split forces) and the report.
    #[must_use]
    #[allow(clippy::too_many_lines)] // one iteration loop with certificates
    pub fn solve(
        &self,
        warm_x: Option<Vec<f64>>,
        warm_y: Option<Vec<f64>>,
        settings: PdhgSettings,
    ) -> (Vec<f64>, Vec<f64>, PdhgReport) {
        let nvar = self.c.len();
        let nrow = self.b.len();
        let mut x = warm_x.unwrap_or_else(|| vec![0.0; nvar]);
        let mut y = warm_y.unwrap_or_else(|| vec![0.0; nrow]);
        let step = 0.95 / self.norm_est.max(1e-30);
        let (tau, sigma) = (step, step);
        let bnorm = self.b.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
        let mut report = PdhgReport::default();
        let mut aty = vec![0.0f64; nvar];
        let mut ax = vec![0.0f64; nrow];
        let mut x_prev = x.clone();
        for it in 0..settings.max_iters {
            // x ← Π₊(x − τ(c + Aᵀy))
            self.at.spmv(&y, &mut aty);
            x_prev.copy_from_slice(&x);
            for i in 0..nvar {
                x[i] = (x[i] - tau * (self.c[i] + aty[i])).max(0.0);
            }
            // y ← y + σ(A(2x − x_prev) − b)
            let xbar: Vec<f64> = x
                .iter()
                .zip(&x_prev)
                .map(|(xi, xp)| 2.0 * xi - xp)
                .collect();
            self.a.spmv(&xbar, &mut ax);
            for r in 0..nrow {
                y[r] += sigma * (ax[r] - self.b[r]);
            }
            if (it + 1) % settings.check_every == 0 || it + 1 == settings.max_iters {
                let (gap, eq_res, primal) = self.certificate(&x, &y, bnorm);
                report.trace.push((it + 1, gap));
                report.iters = it + 1;
                report.volume = primal;
                report.gap = gap;
                report.eq_residual = eq_res;
                if gap < settings.gap_tol && eq_res < settings.gap_tol {
                    break;
                }
            }
        }
        (x, y, report)
    }

    /// The certificate: (relative duality gap, equilibrium residual,
    /// primal objective). With the saddle `cᵀx + yᵀ(Ax − b)` this
    /// solver ascends, the dual objective is `−bᵀy` under feasibility
    /// `c + Aᵀy ≥ 0`; a uniform shrink λ·y restores feasibility with a
    /// still-certified bound `−λ·bᵀy` (c > 0 for every member).
    #[must_use]
    pub fn certificate(&self, x: &[f64], y: &[f64], bnorm: f64) -> (f64, f64, f64) {
        let primal: f64 = self.c.iter().zip(x).map(|(c, x)| c * x).sum();
        let mut aty = vec![0.0f64; self.c.len()];
        self.at.spmv(y, &mut aty);
        let mut scale = 1.0f64;
        for (a, c) in aty.iter().zip(&self.c) {
            // Violation where c + Aᵀy < 0, i.e. aty < −c.
            if *a < -c && *a < 0.0 {
                scale = scale.min(-c / a);
            }
        }
        let dual: f64 = -(y.iter().zip(&self.b).map(|(y, b)| y * b).sum::<f64>()) * scale.max(0.0);
        let mut ax = vec![0.0f64; self.b.len()];
        self.a.spmv(x, &mut ax);
        let eq_res = ax
            .iter()
            .zip(&self.b)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
            / bnorm;
        let gap = (primal - dual).abs() / primal.abs().max(1e-30);
        (gap, eq_res, primal)
    }
}

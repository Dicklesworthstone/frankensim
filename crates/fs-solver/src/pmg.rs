//! p-multigrid over the high-order tensor hierarchy (tfz.10 slice 2).
//!
//! The payoff of fs-feec's hierarchical Lobatto basis: the order-r′
//! space on the same mesh is a SUBSET of the order-r space (bubble k
//! is the same polynomial at every order), so prolongation is EXACT
//! INJECTION into matching lattice slots, restriction is its
//! transpose (entry picking), and the Galerkin coarse operator
//! PᵀAP IS the coarse-order Galerkin operator — every level applies
//! MATRIX-FREE through sum factorization (P6: only the r = 1 coarse
//! problem is ever assembled, and SA-AMG preconditions its CG solve).
//!
//! Smoothing: matrix-free Chebyshev on the Jacobi-scaled operator
//! (fixed-iteration power method for λ_max — deterministic), the
//! many-core-honest choice (no sequential dependencies).

use crate::op::LinearOp;
use fs_feec::TensorSpace;
use fs_sparse::precond::{Precond, SaAmg};
use fs_sparse::{Coo, Csr};

/// One p-MG level: the tensor space, its interior mask, Jacobi
/// diagonal, and Chebyshev bounds.
struct Level {
    space: TensorSpace,
    mask: Vec<bool>,
    diag: Vec<f64>,
    lambda_max: f64,
    /// 1D lattice injection map INTO the next-finer level (absent on
    /// the finest level).
    inject1: Option<Vec<usize>>,
}

/// Matrix-free p-multigrid V-cycle preconditioner for the Poisson
/// operator on the unit-cube tensor spaces (homogeneous Dirichlet).
pub struct PMultigrid {
    levels: Vec<Level>,
    coarse: Csr,
    coarse_amg: SaAmg,
    coarse_interior: Vec<usize>,
    coarse_slot: Vec<usize>,
    /// Chebyshev smoothing degree per pre/post sweep.
    pub smooth_degree: usize,
}

/// The 1D lattice injection from order rc into order rf (same mesh):
/// vertex j → j·rf/rc·rc… vertices map to vertices, bubble k keeps
/// its cell and k (bubbles are order-independent polynomials).
fn inject_1d(m: usize, rc: usize, rf: usize) -> Vec<usize> {
    let nc = m * rc + 1;
    let mut map = vec![0usize; nc];
    for j in 0..=m {
        map[j * rc] = j * rf;
    }
    for c in 0..m {
        for k in 2..=rc {
            map[c * rc + (k - 1)] = c * rf + (k - 1);
        }
    }
    map
}

/// Fixed-iteration power method for λ_max of the Jacobi-scaled masked
/// operator (deterministic start, 20 iterations, 10% safety margin).
fn lambda_max(space: &TensorSpace, mask: &[bool], diag: &[f64]) -> f64 {
    let n = space.ndof();
    let mut v: Vec<f64> = (0..n)
        .map(|i| if mask[i] { 1.0 + (i % 7) as f64 / 7.0 } else { 0.0 })
        .collect();
    let mut lam = 1.0f64;
    for _ in 0..20 {
        let mut av = space.apply_stiffness(&v);
        for i in 0..n {
            av[i] = if mask[i] { av[i] / diag[i] } else { 0.0 };
        }
        let norm = fs_math::det::sqrt(av.iter().map(|x| x * x).sum::<f64>());
        lam = norm;
        for (vi, ai) in v.iter_mut().zip(&av) {
            *vi = ai / norm;
        }
    }
    lam * 1.1
}

impl PMultigrid {
    /// Build the hierarchy for order `r` on an m³ grid: orders halve
    /// down to 1 (e.g. 4 → 2 → 1); the r = 1 level is assembled
    /// (interior Kronecker CSR) and preconditioned with SA-AMG.
    ///
    /// # Panics
    /// If `r < 2` (no hierarchy to build) or `m == 0`.
    #[must_use]
    pub fn new(m: usize, r: usize, smooth_degree: usize) -> PMultigrid {
        assert!(r >= 2, "p-MG needs r >= 2");
        let mut orders = vec![r];
        let mut cur = r;
        while cur > 1 {
            cur = (cur / 2).max(1);
            orders.push(cur);
        }
        let mut levels = Vec::new();
        for (li, &ord) in orders.iter().enumerate() {
            let space = TensorSpace::new(m, ord);
            let mask = space.interior_mask();
            let diag = space.stiffness_diagonal();
            let lmax = lambda_max(&space, &mask, &diag);
            let inject1 = if li == 0 {
                None
            } else {
                Some(inject_1d(m, ord, orders[li - 1]))
            };
            levels.push(Level { space, mask, diag, lambda_max: lmax, inject1 });
        }
        // Assemble the coarse (r = 1) interior operator: Kronecker of
        // assembled 1D matrices restricted to interior dofs.
        let coarse_space = &levels.last().expect("nonempty").space;
        let (m1, k1) = coarse_space.assembled_1d();
        let n1 = coarse_space.n1;
        let mask = &levels.last().expect("nonempty").mask;
        let n_full = coarse_space.ndof();
        let mut slot = vec![usize::MAX; n_full];
        let mut interior = Vec::new();
        for d in 0..n_full {
            if mask[d] {
                slot[d] = interior.len();
                interior.push(d);
            }
        }
        let mut coo = Coo::new(interior.len(), interior.len());
        let gid = |i: usize, j: usize, k: usize| (i * n1 + j) * n1 + k;
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    if !mask[gid(i, j, k)] {
                        continue;
                    }
                    for a in 0..n1 {
                        for b in 0..n1 {
                            for c in 0..n1 {
                                if !mask[gid(a, b, c)] {
                                    continue;
                                }
                                let v = k1[i * n1 + a] * m1[j * n1 + b] * m1[k * n1 + c]
                                    + m1[i * n1 + a] * k1[j * n1 + b] * m1[k * n1 + c]
                                    + m1[i * n1 + a] * m1[j * n1 + b] * k1[k * n1 + c];
                                if v != 0.0 {
                                    coo.push(slot[gid(i, j, k)], slot[gid(a, b, c)], v);
                                }
                            }
                        }
                    }
                }
            }
        }
        let coarse = coo.assemble();
        let coarse_amg = SaAmg::new(&coarse, 0.08, 2);
        PMultigrid {
            levels,
            coarse,
            coarse_amg,
            coarse_interior: interior,
            coarse_slot: slot,
            smooth_degree,
        }
    }

    /// Matrix-free Chebyshev smoothing on level `li`: x ← x + p(D⁻¹A)
    /// applied to the residual.
    fn smooth(&self, li: usize, x: &mut [f64], b: &[f64]) {
        let lv = &self.levels[li];
        let n = lv.space.ndof();
        // Smoothing band: with the hierarchical-injection coarse space
        // the complement modes spread well below λ_max/4, so the band
        // reaches down to λ_max/16 and the degree grows with r to keep
        // per-sweep damping strong (tuned on the ladder battery).
        let (lmax, lmin) = (lv.lambda_max, lv.lambda_max / 16.0);
        let theta = f64::midpoint(lmax, lmin);
        let delta = f64::midpoint(lmax, -lmin);
        let masked_scaled_residual = |x: &[f64]| -> Vec<f64> {
            let mut ax = lv.space.apply_stiffness(x);
            for i in 0..n {
                ax[i] = if lv.mask[i] { (b[i] - ax[i]) / lv.diag[i] } else { 0.0 };
            }
            ax
        };
        // Standard Chebyshev iteration on the scaled residual equation.
        let mut d = masked_scaled_residual(x);
        let mut alpha = 1.0 / theta;
        for di in &mut d {
            *di *= alpha;
        }
        for (xi, di) in x.iter_mut().zip(&d) {
            *xi += di;
        }
        let mut rho_prev = delta / theta;
        for _ in 1..self.smooth_degree {
            let rho = 1.0 / (2.0 * theta / delta - rho_prev);
            let r = masked_scaled_residual(x);
            alpha = 2.0 * rho / delta;
            let beta = rho * rho_prev;
            for i in 0..n {
                d[i] = alpha.mul_add(r[i], beta * d[i]);
            }
            for (xi, di) in x.iter_mut().zip(&d) {
                *xi += di;
            }
            rho_prev = rho;
        }
    }

    /// Restrict a residual from level `li` to `li + 1` (injection
    /// transpose: entry picking through the tensor of 1D maps).
    fn restrict(&self, li: usize, r: &[f64]) -> Vec<f64> {
        let fine = &self.levels[li];
        let coarse = &self.levels[li + 1];
        let map = coarse.inject1.as_ref().expect("coarse levels carry maps");
        let ncc = coarse.space.n1;
        let mut out = vec![0.0f64; coarse.space.ndof()];
        for i in 0..ncc {
            for j in 0..ncc {
                for k in 0..ncc {
                    out[coarse.space.gid(i, j, k)] =
                        r[fine.space.gid(map[i], map[j], map[k])];
                }
            }
        }
        out
    }

    /// Prolong-and-correct from level `li + 1` into `li` (injection).
    fn prolong_add(&self, li: usize, e_coarse: &[f64], x: &mut [f64]) {
        let fine = &self.levels[li];
        let coarse = &self.levels[li + 1];
        let map = coarse.inject1.as_ref().expect("coarse levels carry maps");
        let ncc = coarse.space.n1;
        for i in 0..ncc {
            for j in 0..ncc {
                for k in 0..ncc {
                    x[fine.space.gid(map[i], map[j], map[k])] +=
                        e_coarse[coarse.space.gid(i, j, k)];
                }
            }
        }
    }

    /// One V-cycle on level `li` for A·x = b (x in/out, full-space
    /// vectors with boundary dofs pinned to zero).
    fn vcycle(&self, li: usize, x: &mut [f64], b: &[f64]) {
        if li == self.levels.len() - 1 {
            // Coarse: assembled CSR + SA-AMG-preconditioned CG.
            let nred = self.coarse_interior.len();
            let rhs: Vec<f64> = self.coarse_interior.iter().map(|&d| b[d]).collect();
            let mut xi = vec![0.0f64; nred];
            let _ = fs_sparse::precond::pcg(&self.coarse, &rhs, &mut xi, &self.coarse_amg, 1e-10, 500);
            for (s, &d) in xi.iter().zip(&self.coarse_interior) {
                x[d] = *s;
            }
            let _ = &self.coarse_slot;
            return;
        }
        let lv = &self.levels[li];
        let n = lv.space.ndof();
        self.smooth(li, x, b);
        // Residual, restricted.
        let mut ax = lv.space.apply_stiffness(x);
        for i in 0..n {
            ax[i] = if lv.mask[i] { b[i] - ax[i] } else { 0.0 };
        }
        let r_coarse = self.restrict(li, &ax);
        let mut e_coarse = vec![0.0f64; self.levels[li + 1].space.ndof()];
        self.vcycle(li + 1, &mut e_coarse, &r_coarse);
        self.prolong_add(li, &e_coarse, x);
        self.smooth(li, x, b);
    }
}

impl Precond for PMultigrid {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        z.fill(0.0);
        // Mask the input (boundary dofs carry no residual).
        let masked: Vec<f64> = r
            .iter()
            .zip(&self.levels[0].mask)
            .map(|(v, &m)| if m { *v } else { 0.0 })
            .collect();
        self.vcycle(0, z, &masked);
        for (zi, &m) in z.iter_mut().zip(&self.levels[0].mask) {
            if !m {
                *zi = 0.0;
            }
        }
    }
}

/// The finest-level operator as a [`LinearOp`] (masked full-space
/// Poisson apply — what CG iterates on).
pub struct MaskedTensorOp {
    space: TensorSpace,
    mask: Vec<bool>,
}

impl MaskedTensorOp {
    /// Build for order r on an m³ grid.
    #[must_use]
    pub fn new(m: usize, r: usize) -> MaskedTensorOp {
        let space = TensorSpace::new(m, r);
        let mask = space.interior_mask();
        MaskedTensorOp { space, mask }
    }

    /// The interior mask.
    #[must_use]
    pub fn mask(&self) -> &[bool] {
        &self.mask
    }

    /// The underlying space.
    #[must_use]
    pub fn space(&self) -> &TensorSpace {
        &self.space
    }
}

impl LinearOp for MaskedTensorOp {
    fn n(&self) -> usize {
        self.space.ndof()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let ax = self.space.apply_stiffness(x);
        for i in 0..y.len() {
            y[i] = if self.mask[i] { ax[i] } else { x[i] };
        }
    }
}

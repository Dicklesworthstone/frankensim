//! Saddle-point block preconditioners (bead avuw): Stokes on the unit
//! cube over the FEEC tensor pair Q_r³ / P_{r−1}^disc (velocity
//! components in the H1 Lobatto tensor space; pressure in per-cell
//! TOTAL-DEGREE Legendre modes, whose L2 mass is DIAGONAL, so the
//! pressure-mass Schur approximation is trivially invertible — the
//! payoff the plan promised). MEASURED rejection: the full-tensor
//! Q_{r−1}^disc pressure (the bead's first reading) contains the
//! Q1/Q0 checkerboard pattern and iteration counts GREW with the mesh
//! (44 → 101 → 137 across m = 2..4 at r = 2); the total-degree
//! P_{r−1}^disc subset is the classic uniformly inf-sup-stable
//! discontinuous pair (Q2/P1disc and family) and the counts flatten.
//!
//! The system [[A, Bᵀ], [B, 0]] is driven by PRECONDITIONED MINRES
//! ([`crate::PminresState`]) with the block-diagonal SPD preconditioner
//! blockdiag(p-MG on each velocity component, pressure-mass inverse on
//! the Schur block) — the textbook-optimal Stokes combination
//! (Silvester–Wathen), h-robust because both blocks are.
//!
//! Enclosed flow fixes pressure only up to a constant: the operator
//! and rhs are PROJECTED against the constant-pressure null vector
//! (which the diagonal pressure mass preserves — all cells share one
//! h — so the projection commutes with the preconditioner and MINRES
//! stays in the orthogonal complement).

use crate::op::LinearOp;
use fs_feec::highorder::hex::TensorSpace;
use fs_feec::highorder::quad1d::{gauss_legendre, legendre, lobatto_shapes};
use fs_sparse::precond::Precond;
use fs_sparse::{Coo, Csr};

/// The assembled Stokes saddle system on an m³ mesh at order r
/// (velocity Q_r³ with homogeneous Dirichlet everywhere, pressure
/// P_{r−1}^disc), plus everything the solver and preconditioner need.
pub struct StokesSystem {
    /// Scalar velocity-component space (shared by all 3 components).
    pub space: TensorSpace,
    vmask: Vec<bool>,
    /// Scalar component dof count (velocity total = 3·nv).
    pub nv: usize,
    /// Pressure dof count: m³ cells times total-degree P_{r−1} modes.
    pub np: usize,
    /// Weak divergence, np × 3nv (boundary-velocity columns zeroed).
    b: Csr,
    /// Bᵀ, 3nv × np.
    bt: Csr,
    /// Diagonal pressure mass (Legendre orthogonality).
    pmass: Vec<f64>,
    /// Normalized constant-pressure null vector.
    pnull: Vec<f64>,
}

impl StokesSystem {
    /// Assemble the fixture.
    ///
    /// # Panics
    /// If `m < 2` (the p-MG preconditioner's smoother needs interior
    /// vertices) or `r < 2` (no p-hierarchy).
    #[must_use]
    pub fn new(m: usize, r: usize) -> StokesSystem {
        assert!(m >= 2 && r >= 2, "Stokes fixture needs m >= 2, r >= 2");
        let space = TensorSpace::new(m, r);
        let vmask = space.interior_mask();
        let nv = space.ndof();
        let h = 1.0 / m as f64;

        let modes = pressure_modes(r);
        let nmodes = modes.len();
        let ncells = m * m * m;
        let np = ncells * nmodes;
        let cell_id = |cx: usize, cy: usize, cz: usize| (cx * m + cy) * m + cz;

        // 1D couplings on one cell: V_e[k, l] = ∫ P_k N_l dx (physical),
        // W_e[k, l] = ∫ P_k N_l' dx (h factors cancel).
        let (qx, qw) = gauss_legendre(r + 3);
        let p1 = r + 1;
        let mut ve = vec![0.0f64; r * p1];
        let mut we = vec![0.0f64; r * p1];
        for (&x, &w) in qx.iter().zip(&qw) {
            let (nvals, nders) = lobatto_shapes(r, x);
            for k in 0..r {
                let (pk, _) = legendre(k, x);
                for l in 0..p1 {
                    ve[k * p1 + l] = (w * pk * h / 2.0).mul_add(nvals[l], ve[k * p1 + l]);
                    we[k * p1 + l] = (w * pk).mul_add(nders[l], we[k * p1 + l]);
                }
            }
        }
        // B: per-cell assembly (pressure modes never cross cells).
        let mut coo = Coo::new(np, 3 * nv);
        for cx in 0..m {
            for cy in 0..m {
                for cz in 0..m {
                    let cbase = cell_id(cx, cy, cz) * nmodes;
                    for (mi, &(ka, kb, kc)) in modes.iter().enumerate() {
                        let prow = cbase + mi;
                        for lx in 0..p1 {
                            let i = space.lat1(cx, lx);
                            for ly in 0..p1 {
                                let j = space.lat1(cy, ly);
                                for lz in 0..p1 {
                                    let k = space.lat1(cz, lz);
                                    let vg = space.gid(i, j, k);
                                    if !vmask[vg] {
                                        continue;
                                    }
                                    let ex = we[ka * p1 + lx] * ve[kb * p1 + ly] * ve[kc * p1 + lz];
                                    let ey = ve[ka * p1 + lx] * we[kb * p1 + ly] * ve[kc * p1 + lz];
                                    let ez = ve[ka * p1 + lx] * ve[kb * p1 + ly] * we[kc * p1 + lz];
                                    if ex != 0.0 {
                                        coo.push(prow, vg, ex);
                                    }
                                    if ey != 0.0 {
                                        coo.push(prow, nv + vg, ey);
                                    }
                                    if ez != 0.0 {
                                        coo.push(prow, 2 * nv + vg, ez);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let b_csr = coo.assemble();
        // Bᵀ by re-pushing transposed (deterministic order).
        let mut coot = Coo::new(3 * nv, np);
        for row in 0..np {
            let (cols, vals) = b_csr.row(row);
            for (&col, &val) in cols.iter().zip(vals) {
                coot.push(col, row, val);
            }
        }
        let bt = coot.assemble();
        // Diagonal pressure mass: ∏ h/(2k+1) per axis mode.
        let mut pmass = vec![0.0f64; np];
        for cell in 0..ncells {
            for (mi, &(ka, kb, kc)) in modes.iter().enumerate() {
                pmass[cell * nmodes + mi] =
                    h / (2 * ka + 1) as f64 * (h / (2 * kb + 1) as f64) * (h / (2 * kc + 1) as f64);
            }
        }
        // Constant pressure = mode 0 (P0) in every cell, equal weight.
        let mut pnull = vec![0.0f64; np];
        for cell in 0..ncells {
            pnull[cell * nmodes] = 1.0;
        }
        let nn = fs_math::det::sqrt(pnull.iter().map(|x| x * x).sum::<f64>());
        for v in &mut pnull {
            *v /= nn;
        }
        StokesSystem {
            space,
            vmask,
            nv,
            np,
            b: b_csr,
            bt,
            pmass,
            pnull,
        }
    }

    /// Total system size (3·nv + np).
    #[must_use]
    pub fn n(&self) -> usize {
        3 * self.nv + self.np
    }

    /// The interior mask of one velocity component.
    #[must_use]
    pub fn vmask(&self) -> &[bool] {
        &self.vmask
    }

    /// Project the constant-pressure mode out of the p-block (in place).
    pub fn project_pressure(&self, p: &mut [f64]) {
        let dot: f64 = p.iter().zip(&self.pnull).map(|(a, b)| a * b).sum();
        for (pi, ni) in p.iter_mut().zip(&self.pnull) {
            *pi = dot.mul_add(-ni, *pi);
        }
    }

    /// ‖B·u‖∞ of the velocity part of a solution vector — the
    /// divergence-freeness certificate.
    #[must_use]
    pub fn divergence_inf(&self, xu: &[f64]) -> f64 {
        let mut bu = vec![0.0f64; self.np];
        self.b.spmv(&xu[..3 * self.nv], &mut bu);
        self.project_pressure(&mut bu);
        bu.iter().fold(0.0f64, |m, v| m.max(v.abs()))
    }
}

fn pressure_modes(r: usize) -> Vec<(usize, usize, usize)> {
    let mut modes = Vec::new();
    for ka in 0..r {
        for kb in 0..r {
            for kc in 0..r {
                if ka + kb + kc < r {
                    modes.push((ka, kb, kc));
                }
            }
        }
    }
    modes
}

/// The saddle operator [[A, Bᵀ], [B, 0]] (velocity Laplacian per
/// component with identity on Dirichlet dofs; pressure block projected
/// against the constant mode so the operator is symmetric on the
/// working subspace).
pub struct StokesOp<'a> {
    sys: &'a StokesSystem,
}

impl<'a> StokesOp<'a> {
    /// Wrap a system.
    #[must_use]
    pub fn new(sys: &'a StokesSystem) -> StokesOp<'a> {
        StokesOp { sys }
    }
}

impl LinearOp for StokesOp<'_> {
    fn n(&self) -> usize {
        self.sys.n()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let sys = self.sys;
        let nv = sys.nv;
        // Velocity block: masked vector Laplacian.
        for comp in 0..3 {
            let u = &x[comp * nv..(comp + 1) * nv];
            let au = sys.space.apply_stiffness(u);
            for i in 0..nv {
                y[comp * nv + i] = if sys.vmask[i] { au[i] } else { u[i] };
            }
        }
        // + Bᵀ p (projected input; Bᵀ rows for Dirichlet dofs are empty
        // by construction, so identity rows stay clean).
        let mut p_proj = x[3 * nv..].to_vec();
        sys.project_pressure(&mut p_proj);
        let mut btp = vec![0.0f64; 3 * nv];
        sys.bt.spmv(&p_proj, &mut btp);
        for i in 0..3 * nv {
            y[i] += btp[i];
        }
        // Pressure block: B u, projected.
        let mut bu = vec![0.0f64; sys.np];
        sys.b.spmv(&x[..3 * nv], &mut bu);
        sys.project_pressure(&mut bu);
        y[3 * nv..].copy_from_slice(&bu);
    }
}

/// blockdiag(p-MG per velocity component, diagonal pressure-mass
/// inverse) — the Silvester–Wathen SPD Stokes preconditioner.
pub struct StokesBlockDiag<'a> {
    sys: &'a StokesSystem,
    pmg: crate::PMultigrid,
}

impl<'a> StokesBlockDiag<'a> {
    /// Build (one p-MG shared across the identical components).
    #[must_use]
    pub fn new(sys: &'a StokesSystem, smooth_degree: usize) -> StokesBlockDiag<'a> {
        StokesBlockDiag {
            sys,
            pmg: crate::PMultigrid::new(sys.space.m(), sys.space.r(), smooth_degree),
        }
    }
}

impl Precond for StokesBlockDiag<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let nv = self.sys.nv;
        for comp in 0..3 {
            let mut zc = vec![0.0f64; nv];
            self.pmg.apply(&r[comp * nv..(comp + 1) * nv], &mut zc);
            z[comp * nv..(comp + 1) * nv].copy_from_slice(&zc);
        }
        for i in 0..self.sys.np {
            z[3 * nv + i] = r[3 * nv + i] / self.sys.pmass[i];
        }
    }
}

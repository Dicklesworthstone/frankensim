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
//! Smoothing (bead x08j): Chebyshev accelerating VERTEX-CENTERED
//! overlapping-patch additive Schwarz. Each patch is the tensor window
//! of all dofs on the (up to) 8 elements around an interior vertex;
//! because the global operator is the Kronecker sum of assembled 1D
//! matrices, the patch operator is EXACTLY the Kronecker sum of 1D
//! windows, and local solves use FAST DIAGONALIZATION (per-axis
//! generalized eigenproblems K·v = λ·M·v, sizes ≤ 2r+1, shared across
//! the ≤ 3 distinct 1D trim signatures on the uniform grid) — exact
//! patch inverses, no (2r+1)³ dense factorization. This replaces the
//! Jacobi scaling inside the Chebyshev smoother and kills the
//! CBS-angle p-growth the pointwise smoother had (MEASURED history:
//! element-matrix EBE Schwarz 50 iters and Dirichlet-block ω-Schwarz
//! 41 iters where Jacobi-Chebyshev needed 8 — both in the tfz.10 bead
//! trail; Jacobi-Chebyshev itself grew mildly with r — the x08j
//! ladder gates the flat counts).

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};

use crate::op::LinearOp;
use fs_feec::TensorSpace;
use fs_la::eigen::jacobi_eigh;
use fs_sparse::precond::{Precond, SaAmg};
use fs_sparse::{Coo, Csr};

/// One p-MG level: the tensor space, its interior mask, the vertex
/// Schwarz smoother data, and the Chebyshev spectral bound (of the
/// Schwarz-preconditioned operator).
struct Level {
    space: TensorSpace,
    mask: Vec<bool>,
    schwarz: VertexSchwarz,
    lambda_max: f64,
    /// 1D lattice injection map INTO the next-finer level (absent on
    /// the finest level).
    inject1: Option<Vec<usize>>,
}

// ---------------------------------------------------------------------------
// Vertex-patch Schwarz via fast diagonalization.
// ---------------------------------------------------------------------------

/// One distinct 1D patch factor: the window indices, the generalized
/// eigenpairs (K·v = λ·M·v on the window) with V M-orthonormal
/// (VᵀMV = I, VᵀKV = Λ) — the axis ingredient of fast diagonalization.
struct AxisPatch {
    /// Window LENGTH (indices are per-vertex offsets; the eigendata is
    /// translation-invariant on the uniform grid, the window is not —
    /// sharing the representative's indices misdirected every second
    /// interior vertex at m ≥ 5 and the preconditioner went singular
    /// on the uncovered cells: MEASURED as 29 iters at m=5 and outright
    /// failure at m=6 before the fix).
    nloc: usize,
    /// Row-major nloc×nloc; COLUMNS are the eigenvectors.
    v: Vec<f64>,
    lam: Vec<f64>,
}

/// In-place lower-Cholesky of a small SPD row-major matrix.
fn cholesky(a: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..=i {
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                assert!(s > 0.0, "patch mass matrix must be SPD");
                a[i * n + i] = fs_math::det::sqrt(s);
            } else {
                a[i * n + j] = s / a[j * n + j];
            }
        }
    }
    // Zero the strict upper triangle (storage hygiene).
    for i in 0..n {
        for j in i + 1..n {
            a[i * n + j] = 0.0;
        }
    }
}

fn build_axis_patch(mass: &[f64], stiff: &[f64], n1: usize, idx: &[usize]) -> AxisPatch {
    let nl = idx.len();
    let mut ml = vec![0.0f64; nl * nl];
    let mut kl = vec![0.0f64; nl * nl];
    for (a, &i) in idx.iter().enumerate() {
        for (b, &j) in idx.iter().enumerate() {
            ml[a * nl + b] = mass[i * n1 + j];
            kl[a * nl + b] = stiff[i * n1 + j];
        }
    }
    // Generalized K v = λ M v via M = LLᵀ: C = L⁻¹ K L⁻ᵀ, eig(C) = (Λ, W),
    // V = L⁻ᵀ W (then VᵀMV = I, VᵀKV = Λ).
    let mut l = ml;
    cholesky(&mut l, nl);
    // C = L⁻¹ K L⁻ᵀ: forward-solve L·X = K (column-wise), then
    // C = X L⁻ᵀ i.e. solve L·Cᵀ = Xᵀ (C symmetric).
    let mut x = kl;
    // Forward substitution on each column of X (L X = K).
    for c in 0..nl {
        for i in 0..nl {
            let mut s = x[i * nl + c];
            for k in 0..i {
                s -= l[i * nl + k] * x[k * nl + c];
            }
            x[i * nl + c] = s / l[i * nl + i];
        }
    }
    // Now rows: C Lᵀ = X  ⇒  for each row of C: L (C row)ᵀ = (X row)ᵀ.
    let mut cmat = vec![0.0f64; nl * nl];
    for rrow in 0..nl {
        for i in 0..nl {
            let mut s = x[rrow * nl + i];
            for k in 0..i {
                s -= l[i * nl + k] * cmat[rrow * nl + k];
            }
            cmat[rrow * nl + i] = s / l[i * nl + i];
        }
    }
    let (lam, w) = jacobi_eigh(&cmat, nl);
    // V = L⁻ᵀ W: back-substitute Lᵀ V = W per column.
    let mut v = w;
    for c in 0..nl {
        for i in (0..nl).rev() {
            let mut s = v[i * nl + c];
            for k in i + 1..nl {
                s -= l[k * nl + i] * v[k * nl + c];
            }
            v[i * nl + c] = s / l[i * nl + i];
        }
    }
    AxisPatch { nloc: nl, v, lam }
}

/// Additive vertex-patch Schwarz over the (m−1)³ interior vertices,
/// exact per-patch inverses by fast diagonalization. SPD by
/// construction (Σ RᵀA_patch⁻¹R); used as the inner preconditioner of
/// the Chebyshev smoother, which normalizes its spectrum.
struct VertexSchwarz {
    /// Distinct 1D windows (interior / left-trimmed / right-trimmed).
    types: Vec<AxisPatch>,
    /// Per interior vertex (indexed v − 1, v ∈ 1..m): window start
    /// and eigendata type.
    windows: Vec<(usize, usize)>,
    /// Symmetrized partition-of-unity weights 1/√(patch multiplicity)
    /// per dof (0 on the boundary): S = D^½ Σ RᵀAᵢ⁻¹R D^½ stays SPD
    /// while normalizing the overlap counting — without it λ_max(S·A)
    /// JUMPS when the per-axis window multiplicity first reaches 3
    /// (m ≥ 4) and the m-ladder counts grew 8 → 13 (measured).
    dw: Vec<f64>,
}

impl VertexSchwarz {
    fn new(space: &TensorSpace) -> VertexSchwarz {
        let (m, r, n1) = (space.m(), space.r(), space.n1());
        assert!(m >= 2, "vertex-patch Schwarz needs m >= 2");
        let (mass, stiff) = space.assembled_1d();
        let mut types = Vec::new();
        let mut sig_of = std::collections::HashMap::new();
        let mut windows = Vec::with_capacity(m - 1);
        for v in 1..m {
            let lo = ((v - 1) * r).max(1);
            let hi = ((v + 1) * r).min(n1 - 2);
            let sig = (v == 1, v == m - 1);
            let t = *sig_of.entry(sig).or_insert_with(|| {
                let idx: Vec<usize> = (lo..=hi).collect();
                types.push(build_axis_patch(&mass, &stiff, n1, &idx));
                types.len() - 1
            });
            windows.push((lo, t));
        }
        // Per-dof multiplicity = product of per-axis 1D window counts.
        let mut c1 = vec![0u32; n1];
        for v in 1..m {
            let lo = ((v - 1) * r).max(1);
            let hi = ((v + 1) * r).min(n1 - 2);
            for c in c1.iter_mut().take(hi + 1).skip(lo) {
                *c += 1;
            }
        }
        let mut dw = vec![0.0f64; n1 * n1 * n1];
        for i in 0..n1 {
            for j in 0..n1 {
                for k in 0..n1 {
                    let mult = c1[i] * c1[j] * c1[k];
                    if mult > 0 {
                        dw[(i * n1 + j) * n1 + k] = 1.0 / fs_math::det::sqrt(f64::from(mult));
                    }
                }
            }
        }
        VertexSchwarz { types, windows, dw }
    }

    /// out = D^½ Σ_patches Rᵀ A_patch⁻¹ R D^½ · rhs (additive,
    /// PU-symmetrized; out overwritten).
    fn apply(&self, space: &TensorSpace, rhs: &[f64], out: &mut [f64]) {
        out.fill(0.0);
        let rhs: Vec<f64> = rhs.iter().zip(&self.dw).map(|(v, w)| v * w).collect();
        let rhs = rhs.as_slice();
        let m = space.m();
        for vx in 1..m {
            let (lox, tx) = self.windows[vx - 1];
            let ax = &self.types[tx];
            for vy in 1..m {
                let (loy, ty) = self.windows[vy - 1];
                let ay = &self.types[ty];
                for vz in 1..m {
                    let (loz, tz) = self.windows[vz - 1];
                    let az = &self.types[tz];
                    let (nx, ny, nz) = (ax.nloc, ay.nloc, az.nloc);
                    // Gather the local cube.
                    let mut loc = vec![0.0f64; nx * ny * nz];
                    for a in 0..nx {
                        for b in 0..ny {
                            for c in 0..nz {
                                loc[(a * ny + b) * nz + c] =
                                    rhs[space.gid(lox + a, loy + b, loz + c)];
                            }
                        }
                    }
                    // Fast diagonalization: U = (Vᵀ⊗Vᵀ⊗Vᵀ)·loc,
                    // U /= λ-sum, E = (V⊗V⊗V)·U.
                    let u = tensor3(&ax.v, nx, &loc, ny, nz, true);
                    let u = tensor3_mid(&ay.v, ny, &u, nx, nz, true);
                    let mut u = tensor3_last(&az.v, nz, &u, nx, ny, true);
                    for a in 0..nx {
                        for b in 0..ny {
                            for c in 0..nz {
                                u[(a * ny + b) * nz + c] /= ax.lam[a] + ay.lam[b] + az.lam[c];
                            }
                        }
                    }
                    let e = tensor3(&ax.v, nx, &u, ny, nz, false);
                    let e = tensor3_mid(&ay.v, ny, &e, nx, nz, false);
                    let e = tensor3_last(&az.v, nz, &e, nx, ny, false);
                    // Scatter-add.
                    for a in 0..nx {
                        for b in 0..ny {
                            for c in 0..nz {
                                out[space.gid(lox + a, loy + b, loz + c)] +=
                                    e[(a * ny + b) * nz + c];
                            }
                        }
                    }
                }
            }
        }
        for (o, w) in out.iter_mut().zip(&self.dw) {
            *o *= w;
        }
    }
}

/// Contract the square matrix `v` (row-major, columns = eigenvectors)
/// along the FIRST axis of an nx×ny×nz cube; `transpose` applies Vᵀ.
fn tensor3(v: &[f64], nx: usize, src: &[f64], ny: usize, nz: usize, transpose: bool) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny * nz];
    for a in 0..nx {
        for i in 0..nx {
            let w = if transpose {
                v[i * nx + a]
            } else {
                v[a * nx + i]
            };
            if w == 0.0 {
                continue;
            }
            let (dst, srcp) = (a * ny * nz, i * ny * nz);
            for t in 0..ny * nz {
                out[dst + t] = w.mul_add(src[srcp + t], out[dst + t]);
            }
        }
    }
    out
}

/// Contract along the MIDDLE axis of an nx×ny×nz cube.
fn tensor3_mid(
    v: &[f64],
    ny: usize,
    src: &[f64],
    nx: usize,
    nz: usize,
    transpose: bool,
) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny * nz];
    for x in 0..nx {
        for b in 0..ny {
            for j in 0..ny {
                let w = if transpose {
                    v[j * ny + b]
                } else {
                    v[b * ny + j]
                };
                if w == 0.0 {
                    continue;
                }
                let (dst, srcp) = ((x * ny + b) * nz, (x * ny + j) * nz);
                for t in 0..nz {
                    out[dst + t] = w.mul_add(src[srcp + t], out[dst + t]);
                }
            }
        }
    }
    out
}

/// Contract along the LAST axis of an nx×ny×nz cube.
fn tensor3_last(
    v: &[f64],
    nz: usize,
    src: &[f64],
    nx: usize,
    ny: usize,
    transpose: bool,
) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny * nz];
    for xy in 0..nx * ny {
        let base = xy * nz;
        for c in 0..nz {
            let mut s = 0.0f64;
            for k in 0..nz {
                let w = if transpose {
                    v[k * nz + c]
                } else {
                    v[c * nz + k]
                };
                s = w.mul_add(src[base + k], s);
            }
            out[base + c] = s;
        }
    }
    out
}

/// Matrix-free p-multigrid V-cycle preconditioner for the Poisson
/// operator on the unit-cube tensor spaces (homogeneous Dirichlet).
///
/// V-CYCLE FIXITY IS AN ASSUMPTION, NOT AN ESTABLISHED INVARIANT. Plain
/// CG requires a FIXED symmetric preconditioner, and this one is fixed
/// only while the r = 1 coarse level is solved near-exactly. That solve
/// is REQUESTED at 1e-13 within a 2000-iteration cap; if it ever exits
/// on the cap the V-cycle becomes an inexact, application-dependent
/// operator. The receipts used to be discarded (`let _ = pcg(..)`), so
/// nothing could say whether the assumption held. They are now retained:
/// see [`PMultigrid::coarse_solves_converged`].
pub struct PMultigrid {
    levels: Vec<Level>,
    coarse: Csr,
    coarse_amg: SaAmg,
    coarse_interior: Vec<usize>,
    coarse_slot: Vec<usize>,
    /// Chebyshev smoothing degree per pre/post sweep.
    pub smooth_degree: usize,
    /// Coarse-solve evidence. `Precond::apply` takes `&self`, so the
    /// receipts are accumulated through atomics rather than dropped.
    coarse_solves: AtomicUsize,
    coarse_all_converged: AtomicBool,
    /// Worst coarse relative residual seen, as `f64::to_bits`.
    worst_coarse_rel_residual: AtomicU64,
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

impl PMultigrid {
    /// Fold one coarse-solve receipt into the retained evidence.
    fn record_coarse_solve(&self, report: &fs_sparse::precond::PcgReport) {
        self.coarse_solves.fetch_add(1, Ordering::Relaxed);
        if !report.converged {
            self.coarse_all_converged.store(false, Ordering::Relaxed);
        }
        let candidate = if report.rel_residual.is_finite() {
            report.rel_residual
        } else {
            // A non-finite coarse residual is the worst possible outcome
            // and must not be silently dropped by a `<` comparison.
            f64::INFINITY
        };
        let mut current = self.worst_coarse_rel_residual.load(Ordering::Relaxed);
        loop {
            if f64::from_bits(current) >= candidate {
                break;
            }
            match self.worst_coarse_rel_residual.compare_exchange_weak(
                current,
                candidate.to_bits(),
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
    }

    /// Whether EVERY r = 1 coarse solve performed so far met its 1e-13
    /// request within the 2000-iteration cap.
    ///
    /// This is the evidence for the V-cycle's fixity assumption. `false`
    /// means at least one application solved the coarse level only to
    /// the cap, so the preconditioner varied between applications and a
    /// plain-CG convergence claim built on it must be demoted. A driver
    /// that never asks is asserting the invariant without evidence.
    #[must_use]
    pub fn coarse_solves_converged(&self) -> bool {
        self.coarse_all_converged.load(Ordering::Relaxed)
    }

    /// Worst relative residual any coarse solve finished at (`0.0`
    /// before the first application).
    #[must_use]
    pub fn worst_coarse_rel_residual(&self) -> f64 {
        f64::from_bits(self.worst_coarse_rel_residual.load(Ordering::Relaxed))
    }

    /// How many coarse solves have run (V-cycles plus the smoother's
    /// Pavarino coarse term).
    #[must_use]
    pub fn coarse_solves(&self) -> usize {
        self.coarse_solves.load(Ordering::Relaxed)
    }

    /// Drop the retained coarse-solve evidence so a driver can scope it
    /// to one outer solve.
    pub fn reset_coarse_evidence(&self) {
        self.coarse_solves.store(0, Ordering::Relaxed);
        self.coarse_all_converged.store(true, Ordering::Relaxed);
        self.worst_coarse_rel_residual
            .store(0.0f64.to_bits(), Ordering::Relaxed);
    }

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
        assert!(m >= 2, "vertex-patch Schwarz smoothing needs m >= 2");
        let mut levels = Vec::new();
        for (li, &ord) in orders.iter().enumerate() {
            let space = TensorSpace::new(m, ord);
            let mask = space.interior_mask();
            let schwarz = VertexSchwarz::new(&space);
            let inject1 = if li == 0 {
                None
            } else {
                Some(inject_1d(m, ord, orders[li - 1]))
            };
            levels.push(Level {
                space,
                mask,
                schwarz,
                lambda_max: 0.0, // measured below, once the coarse term exists
                inject1,
            });
        }
        // Assemble the coarse (r = 1) interior operator: Kronecker of
        // assembled 1D matrices restricted to interior dofs.
        let coarse_space = &levels.last().expect("nonempty").space;
        let (m1, k1) = coarse_space.assembled_1d();
        let n1 = coarse_space.n1();
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
        let mut pm = PMultigrid {
            levels,
            coarse,
            coarse_amg,
            coarse_interior: interior,
            coarse_slot: slot,
            smooth_degree,
            coarse_solves: AtomicUsize::new(0),
            coarse_all_converged: AtomicBool::new(true),
            worst_coarse_rel_residual: AtomicU64::new(0.0f64.to_bits()),
        };
        // Measure the Chebyshev bound of the FULL smoother preconditioner
        // (Schwarz + coarse term) per smoothing level.
        for li in 0..pm.levels.len() - 1 {
            pm.levels[li].lambda_max = pm.measure_lambda(li);
        }
        pm
    }

    /// The smoother's preconditioner: additive vertex-patch Schwarz PLUS
    /// the exact r = 1 coarse correction (Pavarino's combination — the
    /// coarse term is what makes the two-level bound independent of the
    /// MESH as well as the order; without it, counts were MEASURED to
    /// double from m = 3 to m = 4 at fixed r).
    fn precond_apply(&self, li: usize, r: &[f64], s: &mut [f64]) {
        let lv = &self.levels[li];
        lv.schwarz.apply(&lv.space, r, s);
        let ord = lv.space.r();
        let n1c = self.levels.last().expect("nonempty").space.n1();
        let nred = self.coarse_interior.len();
        let mut rc = vec![0.0f64; nred];
        for (sl, &dof) in self.coarse_interior.iter().enumerate() {
            let k = dof % n1c;
            let j = (dof / n1c) % n1c;
            let i = dof / (n1c * n1c);
            rc[sl] = r[lv.space.gid(i * ord, j * ord, k * ord)];
        }
        let mut xi = vec![0.0f64; nred];
        let coarse_report =
            fs_sparse::precond::pcg(&self.coarse, &rc, &mut xi, &self.coarse_amg, 1e-13, 2000);
        self.record_coarse_solve(&coarse_report);
        for (sl, &dof) in self.coarse_interior.iter().enumerate() {
            let k = dof % n1c;
            let j = (dof / n1c) % n1c;
            let i = dof / (n1c * n1c);
            s[lv.space.gid(i * ord, j * ord, k * ord)] += xi[sl];
        }
        for (si, &mk) in s.iter_mut().zip(&lv.mask) {
            if !mk {
                *si = 0.0;
            }
        }
    }

    /// Power method for λ_max of the smoother-preconditioned operator
    /// on level `li` (deterministic start, 40 iterations, 20% margin).
    fn measure_lambda(&self, li: usize) -> f64 {
        let lv = &self.levels[li];
        let n = lv.space.ndof();
        let mut v: Vec<f64> = (0..n)
            .map(|i| {
                if lv.mask[i] {
                    (if i % 2 == 0 { 1.0 } else { -1.0 }) * (1.0 + (i % 7) as f64 / 7.0)
                } else {
                    0.0
                }
            })
            .collect();
        let mut sav = vec![0.0f64; n];
        let mut lam = 1.0f64;
        for _ in 0..40 {
            let mut av = lv.space.apply_stiffness(&v);
            for (ai, &mk) in av.iter_mut().zip(&lv.mask) {
                if !mk {
                    *ai = 0.0;
                }
            }
            self.precond_apply(li, &av, &mut sav);
            let norm = fs_math::det::sqrt(sav.iter().map(|x| x * x).sum::<f64>());
            lam = norm;
            for (vi, ai) in v.iter_mut().zip(&sav) {
                *vi = ai / norm;
            }
        }
        lam * 1.2
    }

    /// Chebyshev smoothing on level `li` over the Schwarz-preconditioned
    /// operator S·A: band [λ_max/4, λ_max] — the exact patch solves
    /// keep the S·A spectrum TIGHT (vs /16 for Jacobi scaling), which
    /// is exactly the p-independence mechanism.
    fn smooth(&self, li: usize, x: &mut [f64], b: &[f64]) {
        let lv = &self.levels[li];
        let n = lv.space.ndof();
        let (lmax, lmin) = (lv.lambda_max, lv.lambda_max / 16.0);
        let theta = f64::midpoint(lmax, lmin);
        let delta = f64::midpoint(lmax, -lmin);
        let masked_scaled_residual = |x: &[f64]| -> Vec<f64> {
            let mut ax = lv.space.apply_stiffness(x);
            for i in 0..n {
                ax[i] = if lv.mask[i] { b[i] - ax[i] } else { 0.0 };
            }
            let mut s = vec![0.0f64; n];
            self.precond_apply(li, &ax, &mut s);
            s
        };
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
        let ncc = coarse.space.n1();
        let mut out = vec![0.0f64; coarse.space.ndof()];
        for i in 0..ncc {
            for j in 0..ncc {
                for k in 0..ncc {
                    out[coarse.space.gid(i, j, k)] = r[fine.space.gid(map[i], map[j], map[k])];
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
        let ncc = coarse.space.n1();
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
            // Near-exact coarse solve: a loosely-solved coarse level
            // makes the V-cycle a VARYING preconditioner and breaks
            // plain CG (observed as erratic residual histories). The
            // receipt is RETAINED, not discarded — it is the only
            // evidence that the fixity assumption held on this run.
            let coarse_report =
                fs_sparse::precond::pcg(&self.coarse, &rhs, &mut xi, &self.coarse_amg, 1e-13, 2000);
            self.record_coarse_solve(&coarse_report);
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

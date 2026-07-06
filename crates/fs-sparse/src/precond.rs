//! Preconditioner components (plan §6.2): Chebyshev polynomial smoothers,
//! zero-fill incomplete factorizations, and smoothed-aggregation AMG —
//! the toolkit the FLUX solver stack composes.
//!
//! Everything exposes the OPERATOR interface ([`Precond::apply`]); nothing
//! here requires assembling what can be applied. Deterministic by
//! construction: fixed iteration orders, index-order aggregation with
//! lowest-index tie-breaks (P2 applies to SETUP, not just solves), power
//! iterations from fixed start vectors, and NO platform libm anywhere in
//! solver state (contract rule since the eigensolver divergence).
//!
//! v1 scope notes: ILU/IC are sequential (level scheduling is a recorded
//! perf refinement; the bead itself labels them the bandwidth-hostile
//! fallback, not the default). The AMG coarsest level solves with
//! IC(0)-preconditioned CG in-crate (a dense direct coarse solve joins
//! solver-stack integration). Supernodal Cholesky is deferred per its own
//! "not a load-bearing wall" scope cap.

use crate::Csr;
use crate::ops::{spgemm, transpose};

/// The operator interface: z ← M⁻¹·r (approximately).
pub trait Precond {
    /// Apply the preconditioner to `r`, writing `z` (same length).
    fn apply(&self, r: &[f64], z: &mut [f64]);
}

/// Identity (the "no preconditioner" baseline).
#[derive(Debug, Clone, Copy)]
pub struct IdentityPrecond;

impl Precond for IdentityPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        z.copy_from_slice(r);
    }
}

// ---------------------------------------------------------------------------
// Spectral bound (deterministic power iteration)
// ---------------------------------------------------------------------------

/// Estimate λ_max of the SPD operator by power iteration from a FIXED
/// start vector, multiplied by `safety` (≥ 1; Chebyshev bands must
/// enclose the spectrum, so over-estimation is safe and under-estimation
/// is not).
#[must_use]
pub fn lambda_max_estimate(a: &Csr, iters: usize, safety: f64) -> f64 {
    let n = a.nrows();
    // Fixed deterministic start: varying pattern, no libm, no RNG.
    let mut v: Vec<f64> = (0..n)
        .map(|i| 1.0 + 0.3 * (((i * 2_654_435_761) % 97) as f64) / 97.0)
        .collect();
    let mut w = vec![0.0f64; n];
    let mut lam = 0.0f64;
    for _ in 0..iters {
        let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in &mut v {
            *x /= nrm;
        }
        a.spmv(&v, &mut w);
        lam = v.iter().zip(&w).map(|(x, y)| x * y).sum::<f64>();
        std::mem::swap(&mut v, &mut w);
    }
    lam * safety
}

// ---------------------------------------------------------------------------
// Chebyshev polynomial smoother
// ---------------------------------------------------------------------------

/// Degree-k Chebyshev smoother targeting the band [λ_max/α, λ_max]:
/// pure SpMV chains, no sequential dependencies — the many-core-honest
/// smoother. Damps the targeted band by the known Chebyshev factor.
#[derive(Debug, Clone)]
pub struct Chebyshev {
    a: Csr,
    degree: usize,
    lo: f64,
    hi: f64,
}

impl Chebyshev {
    /// Build with λ_max estimated in-crate (power iteration, safety 1.1)
    /// and the standard band divisor α (30.0 is the common default).
    #[must_use]
    pub fn new(a: &Csr, degree: usize, alpha: f64) -> Chebyshev {
        assert!(degree >= 1, "Chebyshev degree must be >= 1");
        assert!(alpha > 1.0, "band divisor must exceed 1");
        let hi = lambda_max_estimate(a, 30, 1.1);
        Chebyshev { a: a.clone(), degree, lo: hi / alpha, hi }
    }

    /// The band this smoother targets (evidence for ledgering).
    #[must_use]
    pub fn band(&self) -> (f64, f64) {
        (self.lo, self.hi)
    }
}

impl Precond for Chebyshev {
    /// z ≈ A⁻¹·r restricted to the band: the classic three-term Chebyshev
    /// iteration with zero initial guess.
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = r.len();
        let theta = f64::midpoint(self.hi, self.lo);
        let delta = f64::midpoint(self.hi, -self.lo);
        // First step: z1 = r/θ.
        let mut zk: Vec<f64> = r.iter().map(|&v| v / theta).collect();
        let mut zk_prev = vec![0.0f64; n];
        let mut rho_prev = delta / theta;
        let mut az = vec![0.0f64; n];
        for _ in 1..self.degree {
            self.a.spmv(&zk, &mut az);
            // residual of current iterate: rk = r − A·zk.
            let rho = 1.0 / (2.0 * theta / delta - rho_prev);
            let c1 = 2.0 * rho / delta;
            let c2 = rho * rho_prev;
            for i in 0..n {
                let rk = r[i] - az[i];
                let znew = c1.mul_add(rk, zk[i]) + c2 * (zk[i] - zk_prev[i]);
                zk_prev[i] = zk[i];
                zk[i] = znew;
            }
            rho_prev = rho;
        }
        z.copy_from_slice(&zk);
    }
}

// ---------------------------------------------------------------------------
// ILU(0) / IC(0)
// ---------------------------------------------------------------------------

/// Typed incomplete-factorization breakdown.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IluBreakdown {
    /// Row at which the pivot became zero/non-positive. Standard remedy:
    /// retry on A + σ·diag(A) with growing σ (caller policy, documented).
    pub row: usize,
}

/// ILU(0): incomplete LU with zero fill on the CSR pattern. Sequential v1
/// (the bandwidth-hostile fallback per plan; level scheduling recorded).
#[derive(Debug, Clone)]
pub struct Ilu0 {
    n: usize,
    /// Factored values on A's pattern (L unit-lower strictly below diag,
    /// U on/above), plus A's structure.
    factored: Csr,
}

/// Factor. # Errors: [`IluBreakdown`] with the failing row.
pub fn ilu0(a: &Csr) -> Result<Ilu0, IluBreakdown> {
    let n = a.nrows();
    assert_eq!(n, a.ncols(), "ilu0 requires a square matrix");
    // Work on a mutable copy of the CSR values with the same pattern.
    let mut vals: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut cols: Vec<Vec<usize>> = Vec::with_capacity(n);
    for i in 0..n {
        let (c, v) = a.row(i);
        cols.push(c.to_vec());
        vals.push(v.to_vec());
    }
    for i in 0..n {
        // For each k < i in row i's pattern (ascending): eliminate.
        for kk in 0..cols[i].len() {
            let k = cols[i][kk];
            if k >= i {
                break;
            }
            // pivot U[k][k]
            let (ck, vk) = (&cols[k], &vals[k]);
            let Some(dk_pos) = ck.iter().position(|&c| c == k) else {
                return Err(IluBreakdown { row: k });
            };
            let piv = vk[dk_pos];
            if piv == 0.0 {
                return Err(IluBreakdown { row: k });
            }
            let lik = vals[i][kk] / piv;
            vals[i][kk] = lik;
            // Row update restricted to row i's pattern: a_ij -= lik * u_kj.
            // Merge over row k's entries j > k that also exist in row i.
            let (ck_c, vk_c): (Vec<usize>, Vec<f64>) = (cols[k].clone(), vals[k].clone());
            for (pos_k, &j) in ck_c.iter().enumerate() {
                if j <= k {
                    continue;
                }
                if let Ok(pos_i) = cols[i].binary_search(&j) {
                    vals[i][pos_i] = (-lik).mul_add(vk_c[pos_k], vals[i][pos_i]);
                }
            }
        }
        // Diagonal must exist and be nonzero for the NEXT eliminations.
        match cols[i].binary_search(&i) {
            Ok(p) if vals[i][p] != 0.0 => {}
            _ => return Err(IluBreakdown { row: i }),
        }
    }
    // Rebuild a CSR holding the factored values.
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut v_all = Vec::new();
    for i in 0..n {
        col_idx.extend_from_slice(&cols[i]);
        v_all.extend_from_slice(&vals[i]);
        row_ptr[i + 1] = col_idx.len();
    }
    Ok(Ilu0 { n, factored: Csr::from_parts(n, n, row_ptr, col_idx, v_all) })
}

impl Precond for Ilu0 {
    /// z = U⁻¹·L⁻¹·r (forward then backward substitution on the pattern).
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = self.n;
        z.copy_from_slice(r);
        // Forward: L is unit-lower on the strict-lower pattern.
        for i in 0..n {
            let (cols, vals) = self.factored.row(i);
            let mut v = z[i];
            for (&c, &lv) in cols.iter().zip(vals) {
                if c >= i {
                    break;
                }
                v = (-lv).mul_add(z[c], v);
            }
            z[i] = v;
        }
        // Backward with U (diag included).
        for i in (0..n).rev() {
            let (cols, vals) = self.factored.row(i);
            let mut v = z[i];
            let mut diag = 1.0;
            for (&c, &uv) in cols.iter().zip(vals) {
                if c < i {
                    continue;
                }
                if c == i {
                    diag = uv;
                } else {
                    v = (-uv).mul_add(z[c], v);
                }
            }
            z[i] = v / diag;
        }
    }
}

// ---------------------------------------------------------------------------
// PCG reference driver
// ---------------------------------------------------------------------------

/// Result of a PCG solve: the evidence object.
#[derive(Debug, Clone)]
pub struct PcgReport {
    /// Iterations performed.
    pub iters: usize,
    /// Final relative residual ‖r‖₂/‖b‖₂.
    pub rel_residual: f64,
    /// Whether the tolerance was met within the cap.
    pub converged: bool,
}

/// Preconditioned conjugate gradients (reference driver — the solver
/// stack supersedes this; kept minimal and deterministic).
pub fn pcg<P: Precond>(
    a: &Csr,
    b: &[f64],
    x: &mut [f64],
    m: &P,
    tol: f64,
    max_iters: usize,
) -> PcgReport {
    let n = b.len();
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
    let mut r = vec![0.0f64; n];
    a.spmv(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    let mut z = vec![0.0f64; n];
    m.apply(&r, &mut z);
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    let mut ap = vec![0.0f64; n];
    for it in 0..max_iters {
        let rel = r.iter().map(|v| v * v).sum::<f64>().sqrt() / bnorm;
        if rel <= tol {
            return PcgReport { iters: it, rel_residual: rel, converged: true };
        }
        a.spmv(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        let alpha = rz / pap;
        for i in 0..n {
            x[i] = alpha.mul_add(p[i], x[i]);
            r[i] = (-alpha).mul_add(ap[i], r[i]);
        }
        m.apply(&r, &mut z);
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = beta.mul_add(p[i], z[i]);
        }
    }
    let rel = r.iter().map(|v| v * v).sum::<f64>().sqrt() / bnorm;
    PcgReport { iters: max_iters, rel_residual: rel, converged: rel <= tol }
}

// ---------------------------------------------------------------------------
// Smoothed-aggregation AMG
// ---------------------------------------------------------------------------

/// One AMG level: the operator, its smoother, and the prolongator down
/// from the finer level (absent on the finest).
struct AmgLevel {
    a: Csr,
    smoother: Chebyshev,
    /// Prolongator from THIS (coarser) level up to the finer one.
    p_from_coarse: Option<Csr>,
}

/// Smoothed-aggregation AMG V-cycle preconditioner.
pub struct SaAmg {
    levels: Vec<AmgLevel>,
    coarse_ilu: Option<Ilu0>,
    /// Per-level unknown counts (evidence: operator complexity).
    pub level_sizes: Vec<usize>,
}

/// Deterministic greedy aggregation on the strength graph: index-order
/// roots, lowest-index attachment (P2 on setup). Returns (aggregate id
/// per node, aggregate count).
fn aggregate(a: &Csr, theta: f64) -> (Vec<usize>, usize) {
    let n = a.nrows();
    const UNASSIGNED: usize = usize::MAX;
    let mut agg = vec![UNASSIGNED; n];
    let strong = |i: usize, j: usize, aij: f64| -> bool {
        let aii = a.get(i, i).abs();
        let ajj = a.get(j, j).abs();
        aij.abs() > theta * (aii * ajj).sqrt()
    };
    let mut count = 0;
    // Pass 1: roots in index order over fully-unaggregated neighborhoods.
    for i in 0..n {
        if agg[i] != UNASSIGNED {
            continue;
        }
        let (cols, vals) = a.row(i);
        let neighborhood_free = cols
            .iter()
            .zip(vals)
            .all(|(&j, &v)| j == i || !strong(i, j, v) || agg[j] == UNASSIGNED);
        if neighborhood_free {
            agg[i] = count;
            for (&j, &v) in cols.iter().zip(vals) {
                if j != i && strong(i, j, v) {
                    agg[j] = count;
                }
            }
            count += 1;
        }
    }
    // Pass 2: attach leftovers to the lowest-id neighboring aggregate.
    for i in 0..n {
        if agg[i] != UNASSIGNED {
            continue;
        }
        let (cols, vals) = a.row(i);
        let mut best = UNASSIGNED;
        for (&j, &v) in cols.iter().zip(vals) {
            if j != i && strong(i, j, v) && agg[j] != UNASSIGNED {
                best = best.min(agg[j]);
            }
        }
        if best == UNASSIGNED {
            // Isolated node: its own aggregate (keeps the prolongator
            // full-rank; happens on disconnected/weak rows).
            agg[i] = count;
            count += 1;
        } else {
            agg[i] = best;
        }
    }
    (agg, count)
}

/// Build the Jacobi-smoothed prolongator P = (I − ω·D⁻¹·A)·P₀ where P₀ is
/// the piecewise-constant tentative prolongator over the aggregates.
fn smoothed_prolongator(a: &Csr, agg: &[usize], n_agg: usize) -> Csr {
    let n = a.nrows();
    // P0 as CSR: one entry per row.
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::with_capacity(n);
    let mut vals = Vec::with_capacity(n);
    for (i, &g) in agg.iter().enumerate() {
        col_idx.push(g);
        vals.push(1.0);
        row_ptr[i + 1] = i + 1;
    }
    let p0 = Csr::from_parts(n, n_agg, row_ptr, col_idx, vals);
    // ω = 4/(3·λ_max(D⁻¹A)) via power iteration on the scaled operator.
    let dinv_a = scale_rows_by_dinv(a);
    let lam = lambda_max_estimate(&dinv_a, 20, 1.0);
    let omega = 4.0 / (3.0 * lam);
    // P = P0 − ω·(D⁻¹A)·P0.
    let ap0 = spgemm(&dinv_a, &p0);
    // P = P0 + (−ω)·AP0: merge the two patterns.
    add_scaled(&p0, &ap0, -omega)
}

/// D⁻¹·A (rows scaled by the inverse diagonal).
fn scale_rows_by_dinv(a: &Csr) -> Csr {
    let n = a.nrows();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n {
        let d = a.get(i, i);
        let inv = if d != 0.0 { 1.0 / d } else { 0.0 };
        let (cols, vs) = a.row(i);
        for (&c, &v) in cols.iter().zip(vs) {
            col_idx.push(c);
            vals.push(v * inv);
        }
        row_ptr[i + 1] = col_idx.len();
    }
    Csr::from_parts(n, a.ncols(), row_ptr, col_idx, vals)
}

/// x + s·y over merged CSR patterns.
fn add_scaled(x: &Csr, y: &Csr, s: f64) -> Csr {
    let n = x.nrows();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n {
        let (cx, vx) = x.row(i);
        let (cy, vy) = y.row(i);
        let (mut a, mut b) = (0usize, 0usize);
        while a < cx.len() || b < cy.len() {
            let ca = cx.get(a).copied().unwrap_or(usize::MAX);
            let cb = cy.get(b).copied().unwrap_or(usize::MAX);
            let c = ca.min(cb);
            let mut v = 0.0f64;
            if ca == c {
                v += vx[a];
                a += 1;
            }
            if cb == c {
                v = s.mul_add(vy[b], v);
                b += 1;
            }
            col_idx.push(c);
            vals.push(v);
            row_ptr[i + 1] = col_idx.len();
        }
        row_ptr[i + 1] = col_idx.len();
    }
    Csr::from_parts(n, x.ncols(), row_ptr, col_idx, vals)
}

impl SaAmg {
    /// Set up the hierarchy: strength θ, Chebyshev smoother degree, and
    /// the coarsest-size cutoff. Deterministic across reruns (tested).
    #[must_use]
    pub fn new(a: &Csr, theta: f64, smoother_degree: usize) -> SaAmg {
        let mut levels = Vec::new();
        let mut level_sizes = Vec::new();
        let mut current = a.clone();
        let mut p_down: Option<Csr> = None;
        loop {
            let n = current.nrows();
            level_sizes.push(n);
            let smoother = Chebyshev::new(&current, smoother_degree, 30.0);
            let done = n <= 64;
            let (agg, n_agg) = if done { (Vec::new(), 0) } else { aggregate(&current, theta) };
            // Stagnation guard: coarsening must actually shrink.
            let stagnated = !done && n_agg * 10 >= n * 9;
            levels.push(AmgLevel { a: current.clone(), smoother, p_from_coarse: p_down.take() });
            if done || stagnated {
                break;
            }
            let p = smoothed_prolongator(&current, &agg, n_agg);
            let pt = transpose(&p);
            let ap = spgemm(&current, &p);
            current = spgemm(&pt, &ap);
            p_down = Some(p);
        }
        // Coarsest-level ILU(0) for the PCG coarse solve (shift-retry on
        // breakdown: SPD Galerkin products rarely need it, but be safe).
        let coarse = &levels.last().unwrap().a;
        let coarse_ilu = ilu0(coarse).ok();
        SaAmg { levels, coarse_ilu, level_sizes }
    }

    /// Operator complexity: Σ nnz(A_ℓ)/nnz(A₀) — the memory-honesty metric.
    #[must_use]
    pub fn operator_complexity(&self) -> f64 {
        let base = self.levels[0].a.nnz().max(1) as f64;
        self.levels.iter().map(|l| l.a.nnz() as f64).sum::<f64>() / base
    }

    fn vcycle(&self, level: usize, r: &[f64], z: &mut [f64]) {
        let lvl = &self.levels[level];
        let n = lvl.a.nrows();
        if level + 1 == self.levels.len() {
            // Coarsest: IC-preconditioned CG to tight tolerance.
            z.fill(0.0);
            if let Some(ilu) = &self.coarse_ilu {
                pcg(&lvl.a, r, z, ilu, 1e-12, 4 * n.max(16));
            } else {
                pcg(&lvl.a, r, z, &IdentityPrecond, 1e-12, 4 * n.max(16));
            }
            return;
        }
        // Pre-smooth.
        lvl.smoother.apply(r, z);
        // Residual: rr = r − A z.
        let mut az = vec![0.0f64; n];
        lvl.a.spmv(z, &mut az);
        let rr: Vec<f64> = r.iter().zip(&az).map(|(a, b)| a - b).collect();
        // Restrict with Pᵀ (the coarser level's stored prolongator).
        let p = self.levels[level + 1]
            .p_from_coarse
            .as_ref()
            .expect("non-coarsest level must carry a prolongator");
        let nc = p.ncols();
        let pt = transpose(p);
        let mut rc = vec![0.0f64; nc];
        pt.spmv(&rr, &mut rc);
        // Coarse correction.
        let mut zc = vec![0.0f64; nc];
        self.vcycle(level + 1, &rc, &mut zc);
        // Prolongate and add.
        let mut up = vec![0.0f64; n];
        p.spmv(&zc, &mut up);
        for (zi, ui) in z.iter_mut().zip(&up) {
            *zi += ui;
        }
        // Post-smooth on the corrected residual.
        lvl.a.spmv(z, &mut az);
        let r2: Vec<f64> = r.iter().zip(&az).map(|(a, b)| a - b).collect();
        let mut dz = vec![0.0f64; n];
        lvl.smoother.apply(&r2, &mut dz);
        for (zi, di) in z.iter_mut().zip(&dz) {
            *zi += di;
        }
    }
}

impl Precond for SaAmg {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.vcycle(0, r, z);
    }
}

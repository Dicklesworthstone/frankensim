//! An in-house primal–dual interior-point solver for block semidefinite
//! programs with free variables (plan §9.8: the SDP engine under the
//! Lasserre/SOS layer).
//!
//! ## Problem form
//!
//! ```text
//! primal:  minimize   c_fᵀ x_f + Σ_b ⟨C_b, X_b⟩
//!          subject to A_f x_f + Σ_b 𝒜_b(X_b) = b,   X_b ⪰ 0,  x_f free
//! dual:    maximize   bᵀ y
//!          subject to A_fᵀ y = c_f,   Z_b = C_b − 𝒜_b*(y) ⪰ 0
//! ```
//!
//! Every constraint matrix is symmetric and given by its upper-triangle
//! entries `(block, i, j, v)` with `i ≤ j`; an off-diagonal entry stands for
//! `A_ij = A_ji = v`, so it contributes `2·v·X_ij` to `⟨A, X⟩`.
//!
//! ## Algorithm
//!
//! Infeasible-start primal–dual path following with the HKM
//! (Helmberg–Kojima–Monteiro) search direction and Mehrotra's
//! predictor–corrector, the scheme of CSDP/SDPT3. Free variables are handled
//! exactly through the reduced saddle system
//! `(A_fᵀ M⁻¹ A_f) Δx_f = A_fᵀ M⁻¹ h − r_f`, never by splitting. The Schur
//! complement `M_ij = tr(A_i X A_j Z⁻¹)` is formed directly from the sparse
//! entries (cost ∝ (Σ nnz)², i.e. O(N⁴) for an `N × N` Gram block).
//!
//! ## Honesty
//!
//! This solver produces FLOATING-POINT approximate solutions. It never issues
//! a certificate: `crate::verify` re-proves every claimed identity with
//! interval arithmetic and an interval Cholesky positivity test, so an
//! inaccurate or wrong SDP answer can only cost a refusal, never a false
//! claim. Infeasibility statuses are divergence heuristics (no Farkas
//! certificate is verified) and are reported as such. The plan names
//! Burer–Monteiro low-rank first-order methods for scale; this second-order
//! method is the small-problem engine and BM remains a recorded successor.
//!
//! Determinism: fixed iteration order, no threads, no randomness — results
//! are bit-identical for the same input on the same build/ISA.

#![allow(clippy::needless_range_loop)]

use crate::dense;

/// One equality constraint row: `Σ entries ⟨A, X⟩ + Σ free a_k x_k = rhs`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SdpRow {
    /// Upper-triangle matrix entries `(block, i, j, value)`, `i ≤ j`.
    pub entries: Vec<(usize, usize, usize, f64)>,
    /// Free-variable coefficients `(index, value)`.
    pub free: Vec<(usize, f64)>,
    /// Right-hand side.
    pub rhs: f64,
}

/// A block SDP with free variables (see the module docs for the form).
#[derive(Debug, Clone, PartialEq)]
pub struct SdpProblem {
    block_sizes: Vec<usize>,
    n_free: usize,
    rows: Vec<SdpRow>,
    c_entries: Vec<(usize, usize, usize, f64)>,
    c_free: Vec<f64>,
}

/// Structured refusal for malformed problems (P10: name the repair).
#[derive(Debug, Clone, PartialEq)]
pub enum SdpError {
    /// A block index, matrix index, or free index is out of range.
    IndexOutOfRange {
        /// What was out of range.
        what: String,
    },
    /// A coefficient or right-hand side is NaN or infinite.
    NonFinite {
        /// Which datum.
        what: String,
    },
    /// A PSD block was declared with size zero.
    EmptyBlock {
        /// The block.
        block: usize,
    },
    /// The problem has no constraints.
    NoConstraints,
}

impl core::fmt::Display for SdpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SdpError::IndexOutOfRange { what } => write!(f, "SDP index out of range: {what}"),
            SdpError::NonFinite { what } => write!(f, "SDP datum is not finite: {what}"),
            SdpError::EmptyBlock { block } => {
                write!(f, "PSD block {block} has size 0; drop it from the problem")
            }
            SdpError::NoConstraints => write!(f, "SDP has no equality constraints"),
        }
    }
}

impl std::error::Error for SdpError {}

impl SdpProblem {
    /// A problem with the given PSD block sizes and number of free variables,
    /// zero objective, and no rows yet.
    #[must_use]
    pub fn new(block_sizes: Vec<usize>, n_free: usize) -> SdpProblem {
        SdpProblem {
            block_sizes,
            n_free,
            rows: Vec::new(),
            c_entries: Vec::new(),
            c_free: vec![0.0; n_free],
        }
    }

    /// Append a constraint row; returns its index.
    pub fn add_row(&mut self, row: SdpRow) -> usize {
        self.rows.push(row);
        self.rows.len() - 1
    }

    /// Set the objective coefficient of free variable `k` (minimized).
    pub fn set_free_cost(&mut self, k: usize, c: f64) {
        if k < self.c_free.len() {
            self.c_free[k] = c;
        }
    }

    /// Add an objective entry `C_b[i][j] (= C_b[j][i]) += v`.
    pub fn add_cost_entry(&mut self, block: usize, i: usize, j: usize, v: f64) {
        self.c_entries.push((block, i.min(j), i.max(j), v));
    }

    /// PSD block sizes.
    #[must_use]
    pub fn block_sizes(&self) -> &[usize] {
        &self.block_sizes
    }

    /// Number of free variables.
    #[must_use]
    pub fn n_free(&self) -> usize {
        self.n_free
    }

    /// Number of constraint rows.
    #[must_use]
    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    /// The constraint rows.
    #[must_use]
    pub fn rows(&self) -> &[SdpRow] {
        &self.rows
    }

    fn validate(&self) -> Result<(), SdpError> {
        if self.rows.is_empty() {
            return Err(SdpError::NoConstraints);
        }
        for (b, &n) in self.block_sizes.iter().enumerate() {
            if n == 0 {
                return Err(SdpError::EmptyBlock { block: b });
            }
        }
        let check = |b: usize, i: usize, j: usize, v: f64, ctx: &str| -> Result<(), SdpError> {
            let n = *self
                .block_sizes
                .get(b)
                .ok_or_else(|| SdpError::IndexOutOfRange {
                    what: format!("{ctx}: block {b}"),
                })?;
            if i >= n || j >= n {
                return Err(SdpError::IndexOutOfRange {
                    what: format!("{ctx}: entry ({i},{j}) in block {b} of size {n}"),
                });
            }
            if !v.is_finite() {
                return Err(SdpError::NonFinite {
                    what: format!("{ctx}: entry ({b},{i},{j})"),
                });
            }
            Ok(())
        };
        for (r, row) in self.rows.iter().enumerate() {
            for &(b, i, j, v) in &row.entries {
                check(b, i, j, v, &format!("row {r}"))?;
            }
            for &(k, v) in &row.free {
                if k >= self.n_free {
                    return Err(SdpError::IndexOutOfRange {
                        what: format!("row {r}: free variable {k} of {}", self.n_free),
                    });
                }
                if !v.is_finite() {
                    return Err(SdpError::NonFinite {
                        what: format!("row {r}: free coefficient {k}"),
                    });
                }
            }
            if !row.rhs.is_finite() {
                return Err(SdpError::NonFinite {
                    what: format!("row {r}: rhs"),
                });
            }
        }
        for &(b, i, j, v) in &self.c_entries {
            check(b, i, j, v, "objective")?;
        }
        for (k, c) in self.c_free.iter().enumerate() {
            if !c.is_finite() {
                return Err(SdpError::NonFinite {
                    what: format!("objective free coefficient {k}"),
                });
            }
        }
        Ok(())
    }
}

/// Solver controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SdpSettings {
    /// Iteration cap.
    pub max_iter: usize,
    /// Relative primal/dual infeasibility tolerance.
    pub tol_feas: f64,
    /// Relative duality-gap tolerance.
    pub tol_gap: f64,
    /// Looser tolerance under which a run that cannot make further progress
    /// (stall, factorization breakdown, iteration cap) is reported
    /// [`SdpStatus::NearOptimal`] instead of failing.
    pub tol_inaccurate: f64,
    /// Fraction-to-the-boundary step factor in `(0, 1)`.
    pub step_fraction: f64,
}

impl Default for SdpSettings {
    fn default() -> Self {
        SdpSettings {
            max_iter: 120,
            tol_feas: 1e-8,
            tol_gap: 1e-8,
            tol_inaccurate: 1e-6,
            step_fraction: 0.95,
        }
    }
}

/// How the solve ended. Only `Optimal` means the tolerances were met.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdpStatus {
    /// Relative infeasibilities and gap below tolerance.
    Optimal,
    /// Progress stopped (stall, breakdown, or iteration cap) with every
    /// measure below `tol_inaccurate` — usable, but less accurate.
    NearOptimal,
    /// The dual objective diverged with a nearly feasible ray: the primal is
    /// (heuristically) infeasible.
    PrimalInfeasible,
    /// The primal objective diverged: the dual is (heuristically) infeasible
    /// / the primal is unbounded.
    DualInfeasible,
    /// The iteration cap was reached.
    MaxIterations,
    /// Step lengths collapsed before the tolerances were met.
    Stalled,
    /// A factorization failed beyond regularization.
    NumericalFailure,
}

/// The solver's final iterate and diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct SdpSolution {
    /// Termination status.
    pub status: SdpStatus,
    /// Primal PSD blocks, dense row-major.
    pub x: Vec<Vec<f64>>,
    /// Primal free variables.
    pub x_free: Vec<f64>,
    /// Dual multipliers, one per row.
    pub y: Vec<f64>,
    /// Dual slack blocks, dense row-major.
    pub z: Vec<Vec<f64>>,
    /// Primal objective value.
    pub primal_objective: f64,
    /// Dual objective value.
    pub dual_objective: f64,
    /// Relative primal infeasibility.
    pub primal_infeasibility: f64,
    /// Relative dual infeasibility.
    pub dual_infeasibility: f64,
    /// Relative duality gap.
    pub relative_gap: f64,
    /// Iterations performed.
    pub iterations: usize,
}

/// Full (both-triangle) positions of a row, grouped by block.
struct RowPositions {
    /// `(block, a, b, v)` with both `(i,j)` and `(j,i)` for off-diagonals.
    pos: Vec<(usize, usize, usize, f64)>,
}

fn expand(entries: &[(usize, usize, usize, f64)]) -> Vec<(usize, usize, usize, f64)> {
    let mut out = Vec::with_capacity(entries.len() * 2);
    for &(b, i, j, v) in entries {
        out.push((b, i, j, v));
        if i != j {
            out.push((b, j, i, v));
        }
    }
    out.sort_by_key(|x| (x.0, x.1, x.2));
    out
}

struct Work<'a> {
    p: &'a SdpProblem,
    rows: Vec<RowPositions>,
    c_pos: Vec<(usize, usize, usize, f64)>,
    af: Vec<f64>, // m × n_free dense
    m: usize,
    nf: usize,
}

impl Work<'_> {
    /// 𝒜(G) for arbitrary (possibly nonsymmetric) blocks.
    fn apply_a(&self, g: &[Vec<f64>]) -> Vec<f64> {
        self.rows
            .iter()
            .map(|r| {
                r.pos
                    .iter()
                    .map(|&(b, i, j, v)| v * g[b][i * self.p.block_sizes[b] + j])
                    .sum()
            })
            .collect()
    }

    /// 𝒜*(y) as dense symmetric blocks.
    fn apply_at(&self, y: &[f64]) -> Vec<Vec<f64>> {
        let mut out: Vec<Vec<f64>> = self
            .p
            .block_sizes
            .iter()
            .map(|&n| vec![0.0; n * n])
            .collect();
        for (r, row) in self.rows.iter().enumerate() {
            if y[r] == 0.0 {
                continue;
            }
            for &(b, i, j, v) in &row.pos {
                out[b][i * self.p.block_sizes[b] + j] += y[r] * v;
            }
        }
        out
    }

    fn c_dense(&self) -> Vec<Vec<f64>> {
        let mut out: Vec<Vec<f64>> = self
            .p
            .block_sizes
            .iter()
            .map(|&n| vec![0.0; n * n])
            .collect();
        for &(b, i, j, v) in &self.c_pos {
            out[b][i * self.p.block_sizes[b] + j] += v;
        }
        out
    }

    fn af_times(&self, x: &[f64]) -> Vec<f64> {
        (0..self.m)
            .map(|r| (0..self.nf).map(|k| self.af[r * self.nf + k] * x[k]).sum())
            .collect()
    }

    fn aft_times(&self, y: &[f64]) -> Vec<f64> {
        (0..self.nf)
            .map(|k| (0..self.m).map(|r| self.af[r * self.nf + k] * y[r]).sum())
            .collect()
    }

    /// Schur complement `M_ij = Σ_b tr(A_i X A_j Z⁻¹)` (symmetric PD).
    fn schur(&self, x: &[Vec<f64>], zinv: &[Vec<f64>]) -> Vec<f64> {
        let m = self.m;
        let mut s = vec![0.0; m * m];
        for i in 0..m {
            let pi = &self.rows[i].pos;
            for j in i..m {
                let pj = &self.rows[j].pos;
                let mut acc = 0.0;
                // Both lists are sorted by block; walk matching blocks.
                let (mut ia, mut ja) = (0, 0);
                while ia < pi.len() && ja < pj.len() {
                    let bi = pi[ia].0;
                    let bj = pj[ja].0;
                    if bi < bj {
                        ia += 1;
                        continue;
                    }
                    if bj < bi {
                        ja += 1;
                        continue;
                    }
                    let b = bi;
                    let n = self.p.block_sizes[b];
                    let ie = pi[ia..]
                        .iter()
                        .position(|e| e.0 != b)
                        .map_or(pi.len(), |o| ia + o);
                    let je = pj[ja..]
                        .iter()
                        .position(|e| e.0 != b)
                        .map_or(pj.len(), |o| ja + o);
                    let xb = &x[b];
                    let zb = &zinv[b];
                    for &(_, a, bb, v) in &pi[ia..ie] {
                        for &(_, c, d, w) in &pj[ja..je] {
                            acc += v * w * xb[bb * n + c] * zb[d * n + a];
                        }
                    }
                    ia = ie;
                    ja = je;
                }
                s[i * m + j] = acc;
                s[j * m + i] = acc;
            }
        }
        s
    }
}

/// Dense row-major PSD blocks.
type Blocks = Vec<Vec<f64>>;
/// A search direction `(ΔX, Δy, ΔZ, Δx_f)`.
type Direction = (Blocks, Vec<f64>, Blocks, Vec<f64>);

fn frob(a: &[f64]) -> f64 {
    a.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn inner(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Cholesky with escalating diagonal regularization.
fn chol_reg(a: &[f64], n: usize) -> Option<Vec<f64>> {
    if let Some(l) = dense::cholesky(a, n) {
        return Some(l);
    }
    let dmax = (0..n)
        .fold(0.0f64, |s, i| s.max(a[i * n + i].abs()))
        .max(1e-300);
    let mut delta = 1e-14;
    while delta <= 1e-4 {
        let mut b = a.to_vec();
        for i in 0..n {
            b[i * n + i] += delta * dmax;
        }
        if let Some(l) = dense::cholesky(&b, n) {
            return Some(l);
        }
        delta *= 100.0;
    }
    None
}

/// Largest `α` (capped at `cap`) keeping `X + α ΔX ⪰ 0`, given `chol(X)`.
fn max_step(lx: &[f64], dx: &[f64], n: usize, cap: f64) -> f64 {
    let t = dense::congruence_inv(lx, dx, n);
    let lmin = dense::sym_eigenvalues(&t, n)[0];
    if lmin >= 0.0 {
        cap
    } else {
        (-1.0 / lmin).min(cap)
    }
}

/// Solve the block SDP.
///
/// # Errors
/// [`SdpError`] when the problem data are malformed. Numerical outcomes
/// (including non-convergence) are reported in [`SdpSolution::status`].
#[allow(clippy::too_many_lines)] // one predictor–corrector loop, kept in one place
pub fn solve(problem: &SdpProblem, settings: &SdpSettings) -> Result<SdpSolution, SdpError> {
    problem.validate()?;
    let m = problem.rows.len();
    let nf = problem.n_free;
    let nb = problem.block_sizes.len();
    let mut af = vec![0.0; m * nf];
    for (r, row) in problem.rows.iter().enumerate() {
        for &(k, v) in &row.free {
            af[r * nf + k] += v;
        }
    }
    let w = Work {
        p: problem,
        rows: problem
            .rows
            .iter()
            .map(|r| RowPositions {
                pos: expand(&r.entries),
            })
            .collect(),
        c_pos: expand(&problem.c_entries),
        af,
        m,
        nf,
    };
    let b: Vec<f64> = problem.rows.iter().map(|r| r.rhs).collect();
    let c = w.c_dense();
    let cf = &problem.c_free;
    let n_tot: usize = problem.block_sizes.iter().sum();
    let norm_b = frob(&b);
    let norm_c = c.iter().map(|cb| frob(cb).powi(2)).sum::<f64>().sqrt() + frob(cf);

    // Initial point (SDPT3-style scaling).
    let row_norms: Vec<f64> = w
        .rows
        .iter()
        .map(|r| r.pos.iter().map(|e| e.3 * e.3).sum::<f64>().sqrt())
        .collect();
    let mut x: Vec<Vec<f64>> = Vec::with_capacity(nb);
    let mut z: Vec<Vec<f64>> = Vec::with_capacity(nb);
    for (bi, &n) in problem.block_sizes.iter().enumerate() {
        let sq = (n as f64).sqrt();
        let mut xi = 10.0f64.max(sq);
        for r in 0..m {
            xi = xi.max(sq * (1.0 + b[r].abs()) / (1.0 + row_norms[r]));
        }
        let cb = frob(&c[bi]);
        let mut eta = 10.0f64.max(sq).max(cb);
        for &rn in &row_norms {
            eta = eta.max(rn);
        }
        let mut xb = vec![0.0; n * n];
        let mut zb = vec![0.0; n * n];
        for i in 0..n {
            xb[i * n + i] = xi;
            zb[i * n + i] = eta;
        }
        x.push(xb);
        z.push(zb);
    }
    let mut y = vec![0.0; m];
    let mut xf = vec![0.0; nf];

    let status;
    let mut iterations = 0;
    let mut stalls = 0;
    let (mut pobj, mut dobj, mut pinf, mut dinf, mut gap);
    loop {
        // Residuals and measures at the current iterate.
        let ax = w.apply_a(&x);
        let afx = w.af_times(&xf);
        let rp: Vec<f64> = (0..m).map(|r| b[r] - ax[r] - afx[r]).collect();
        let aty = w.apply_at(&y);
        let rd: Vec<Vec<f64>> = (0..nb)
            .map(|bi| {
                c[bi]
                    .iter()
                    .zip(&aty[bi])
                    .zip(&z[bi])
                    .map(|((cv, av), zv)| cv - av - zv)
                    .collect()
            })
            .collect();
        let afty = w.aft_times(&y);
        let rf: Vec<f64> = (0..nf).map(|k| cf[k] - afty[k]).collect();
        pobj = inner(cf, &xf) + (0..nb).map(|bi| inner(&c[bi], &x[bi])).sum::<f64>();
        dobj = inner(&b, &y);
        let xz: f64 = (0..nb).map(|bi| inner(&x[bi], &z[bi])).sum();
        let mu = xz / n_tot as f64;
        pinf = frob(&rp) / (1.0 + norm_b);
        dinf = (rd.iter().map(|r| frob(r).powi(2)).sum::<f64>() + frob(&rf).powi(2)).sqrt()
            / (1.0 + norm_c);
        gap = (pobj - dobj).abs() / (1.0 + pobj.abs() + dobj.abs());
        if pinf < settings.tol_feas && dinf < settings.tol_feas && gap < settings.tol_gap {
            status = SdpStatus::Optimal;
            break;
        }
        if dobj > 1e10 * (1.0 + norm_c) && pinf > settings.tol_feas {
            status = SdpStatus::PrimalInfeasible;
            break;
        }
        if -pobj > 1e10 * (1.0 + norm_b) && dinf > settings.tol_feas {
            status = SdpStatus::DualInfeasible;
            break;
        }
        if iterations >= settings.max_iter {
            status = SdpStatus::MaxIterations;
            break;
        }
        iterations += 1;

        // Factorizations.
        let mut lx = Vec::with_capacity(nb);
        let mut zinv = Vec::with_capacity(nb);
        let mut lz = Vec::with_capacity(nb);
        let mut ok = true;
        for bi in 0..nb {
            let n = problem.block_sizes[bi];
            if let (Some(a), Some(bz)) = (dense::cholesky(&x[bi], n), dense::cholesky(&z[bi], n)) {
                zinv.push(dense::chol_inverse(&bz, n));
                lx.push(a);
                lz.push(bz);
            } else {
                ok = false;
                break;
            }
        }
        if !ok {
            status = SdpStatus::NumericalFailure;
            break;
        }
        let schur = w.schur(&x, &zinv);
        let Some(lm) = chol_reg(&schur, m) else {
            status = SdpStatus::NumericalFailure;
            break;
        };
        // Free-variable reduced system pieces.
        let mut minv_af = vec![0.0; m * nf];
        let mut ls = Vec::new();
        if nf > 0 {
            let mut col = vec![0.0; m];
            for k in 0..nf {
                for r in 0..m {
                    col[r] = w.af[r * nf + k];
                }
                dense::chol_solve(&lm, m, &mut col);
                for r in 0..m {
                    minv_af[r * nf + k] = col[r];
                }
            }
            let mut s = vec![0.0; nf * nf];
            for k in 0..nf {
                for l in 0..nf {
                    s[k * nf + l] = (0..m).map(|r| w.af[r * nf + k] * minv_af[r * nf + l]).sum();
                }
            }
            dense::symmetrize(&mut s, nf);
            if let Some(l) = chol_reg(&s, nf) {
                ls = l;
            } else {
                status = SdpStatus::NumericalFailure;
                break;
            }
        }
        let solve_dir = |h: &[f64]| -> (Vec<f64>, Vec<f64>) {
            let mut u = h.to_vec();
            dense::chol_solve(&lm, m, &mut u);
            if nf == 0 {
                return (u, Vec::new());
            }
            let atu = w.aft_times(&u);
            let mut dxf: Vec<f64> = (0..nf).map(|k| atu[k] - rf[k]).collect();
            dense::chol_solve(&ls, nf, &mut dxf);
            let dy: Vec<f64> = (0..m)
                .map(|r| u[r] - (0..nf).map(|k| minv_af[r * nf + k] * dxf[k]).sum::<f64>())
                .collect();
            (dy, dxf)
        };
        // Direction for a given centering target and second-order term.
        let direction = |sigma_mu: f64, corr: Option<(&Blocks, &Blocks)>| -> Direction {
            // T = σμ Z⁻¹ − X − X R_d Z⁻¹ − ΔXa ΔZa Z⁻¹
            let mut t: Vec<Vec<f64>> = Vec::with_capacity(nb);
            for bi in 0..nb {
                let n = problem.block_sizes[bi];
                let xr = dense::matmul(&x[bi], &rd[bi], n);
                let mut tb = dense::matmul(&xr, &zinv[bi], n);
                if let Some((dxa, dza)) = corr {
                    let p = dense::matmul(&dxa[bi], &dza[bi], n);
                    let pz = dense::matmul(&p, &zinv[bi], n);
                    for k in 0..n * n {
                        tb[k] += pz[k];
                    }
                }
                for k in 0..n * n {
                    tb[k] = sigma_mu * zinv[bi][k] - x[bi][k] - tb[k];
                }
                t.push(tb);
            }
            let at = w.apply_a(&t);
            let h: Vec<f64> = (0..m).map(|r| rp[r] - at[r]).collect();
            let (dy, dxf) = solve_dir(&h);
            let atdy = w.apply_at(&dy);
            let mut dz: Vec<Vec<f64>> = Vec::with_capacity(nb);
            let mut dx: Vec<Vec<f64>> = Vec::with_capacity(nb);
            for bi in 0..nb {
                let n = problem.block_sizes[bi];
                let dzb: Vec<f64> = rd[bi].iter().zip(&atdy[bi]).map(|(r, a)| r - a).collect();
                // ΔX = σμZ⁻¹ − X − X ΔZ Z⁻¹ [− ΔXa ΔZa Z⁻¹]
                let xdz = dense::matmul(&x[bi], &dzb, n);
                let xdzz = dense::matmul(&xdz, &zinv[bi], n);
                let mut dxb: Vec<f64> = vec![0.0; n * n];
                for k in 0..n * n {
                    dxb[k] = sigma_mu * zinv[bi][k] - x[bi][k] - xdzz[k];
                }
                if let Some((dxa, dza)) = corr {
                    let p = dense::matmul(&dxa[bi], &dza[bi], n);
                    let pz = dense::matmul(&p, &zinv[bi], n);
                    for k in 0..n * n {
                        dxb[k] -= pz[k];
                    }
                }
                dense::symmetrize(&mut dxb, n);
                dz.push(dzb);
                dx.push(dxb);
            }
            (dx, dy, dz, dxf)
        };
        // Predictor.
        let (dxa, _dya, dza, _dxfa) = direction(0.0, None);
        let mut ap = 1.0f64;
        let mut ad = 1.0f64;
        for bi in 0..nb {
            let n = problem.block_sizes[bi];
            ap = ap.min(max_step(&lx[bi], &dxa[bi], n, 1.0));
            ad = ad.min(max_step(&lz[bi], &dza[bi], n, 1.0));
        }
        let mut xz_aff = 0.0;
        for bi in 0..nb {
            for k in 0..x[bi].len() {
                xz_aff += (x[bi][k] + ap * dxa[bi][k]) * (z[bi][k] + ad * dza[bi][k]);
            }
        }
        let mu_aff = xz_aff / n_tot as f64;
        let sigma = if mu > 0.0 {
            (mu_aff / mu).clamp(0.0, 1.0).powi(3)
        } else {
            0.0
        };
        // Corrector.
        let (dx, dy, dz, dxf) = direction(sigma * mu, Some((&dxa, &dza)));
        let mut ap = f64::INFINITY;
        let mut ad = f64::INFINITY;
        for bi in 0..nb {
            let n = problem.block_sizes[bi];
            ap = ap.min(max_step(&lx[bi], &dx[bi], n, f64::INFINITY));
            ad = ad.min(max_step(&lz[bi], &dz[bi], n, f64::INFINITY));
        }
        let ap = (settings.step_fraction * ap).min(1.0);
        let ad = (settings.step_fraction * ad).min(1.0);
        if ap < 1e-10 && ad < 1e-10 {
            stalls += 1;
            if stalls >= 3 {
                status = SdpStatus::Stalled;
                break;
            }
        } else {
            stalls = 0;
        }
        for bi in 0..nb {
            for k in 0..x[bi].len() {
                x[bi][k] += ap * dx[bi][k];
                z[bi][k] += ad * dz[bi][k];
            }
            let n = problem.block_sizes[bi];
            dense::symmetrize(&mut x[bi], n);
            dense::symmetrize(&mut z[bi], n);
        }
        for k in 0..nf {
            xf[k] += ap * dxf[k];
        }
        for r in 0..m {
            y[r] += ad * dy[r];
        }
    }
    let status = match status {
        SdpStatus::Stalled | SdpStatus::NumericalFailure | SdpStatus::MaxIterations
            if pinf < settings.tol_inaccurate
                && dinf < settings.tol_inaccurate
                && gap < settings.tol_inaccurate =>
        {
            SdpStatus::NearOptimal
        }
        other => other,
    };
    Ok(SdpSolution {
        status,
        x,
        x_free: xf,
        y,
        z,
        primal_objective: pobj,
        dual_objective: dobj,
        primal_infeasibility: pinf,
        dual_infeasibility: dinf,
        relative_gap: gap,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// min ⟨C, X⟩ s.t. tr X = 1, X ⪰ 0 has value λ_min(C).
    #[test]
    fn min_eigenvalue_sdp() {
        let mut p = SdpProblem::new(vec![2], 0);
        p.add_cost_entry(0, 0, 0, 2.0);
        p.add_cost_entry(0, 0, 1, 1.0);
        p.add_cost_entry(0, 1, 1, 2.0);
        p.add_row(SdpRow {
            entries: vec![(0, 0, 0, 1.0), (0, 1, 1, 1.0)],
            free: vec![],
            rhs: 1.0,
        });
        let s = solve(&p, &SdpSettings::default()).unwrap();
        assert_eq!(s.status, SdpStatus::Optimal, "{s:?}");
        assert!(
            (s.primal_objective - 1.0).abs() < 1e-7,
            "{}",
            s.primal_objective
        );
        assert!((s.dual_objective - 1.0).abs() < 1e-7);
    }

    /// A free variable: min −t s.t. X = [[1, 0],[0, 1]] − t I ... expressed as
    /// X_00 + t = 1, X_11 + t = 3, X_01 = 0 → t* = 1.
    #[test]
    fn free_variable_is_exact() {
        let mut p = SdpProblem::new(vec![2], 1);
        p.set_free_cost(0, -1.0);
        p.add_row(SdpRow {
            entries: vec![(0, 0, 0, 1.0)],
            free: vec![(0, 1.0)],
            rhs: 1.0,
        });
        p.add_row(SdpRow {
            entries: vec![(0, 1, 1, 1.0)],
            free: vec![(0, 1.0)],
            rhs: 3.0,
        });
        p.add_row(SdpRow {
            entries: vec![(0, 0, 1, 1.0)],
            free: vec![],
            rhs: 0.0,
        });
        let s = solve(&p, &SdpSettings::default()).unwrap();
        assert_eq!(s.status, SdpStatus::Optimal, "{s:?}");
        assert!((s.x_free[0] - 1.0).abs() < 1e-7, "{:?}", s.x_free);
    }

    #[test]
    fn malformed_problems_are_refused() {
        let mut p = SdpProblem::new(vec![2], 0);
        assert_eq!(
            solve(&p, &SdpSettings::default()).unwrap_err(),
            SdpError::NoConstraints
        );
        p.add_row(SdpRow {
            entries: vec![(0, 0, 5, 1.0)],
            free: vec![],
            rhs: 1.0,
        });
        assert!(matches!(
            solve(&p, &SdpSettings::default()),
            Err(SdpError::IndexOutOfRange { .. })
        ));
    }

    #[test]
    fn infeasible_primal_is_flagged() {
        // X_00 = −1 with X ⪰ 0 is infeasible.
        let mut p = SdpProblem::new(vec![1], 0);
        p.add_row(SdpRow {
            entries: vec![(0, 0, 0, 1.0)],
            free: vec![],
            rhs: -1.0,
        });
        let s = solve(&p, &SdpSettings::default()).unwrap();
        assert_ne!(s.status, SdpStatus::Optimal);
    }
}

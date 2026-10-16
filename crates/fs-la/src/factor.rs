//! Dense factorizations (plan §6.1): blocked right-looking Cholesky,
//! partial-pivot LU, blocked Householder QR (compact-WY), deterministic
//! TSQR, one-sided Jacobi SVD, and 1-norm condition estimation.
//!
//! Determinism: fixed loop orders, fused `mul_add` accumulation, pivot
//! tie-breaks by LOWEST index (P2), TSQR tree shape a pure function of the
//! row-block count — bit-deterministic cross-ISA by construction,
//! golden-hashed. Blocked trailing updates go through [`crate::gemm`], so
//! its KC bit-contract is inherited here.
//!
//! Failure is DATA, not panic: singular pivots and non-SPD diagonals
//! surface as [`FactorError`] with the offending index — the structured
//! diagnostic the plan requires (rank decisions belong to callers).
//!
//! v1 is the single-threaded correctness core; fs-exec tile-parallel
//! panel/update drivers, arena packing, and the FrankenScipy cross-check
//! battery are recorded follow-up scope (bead comments).

/// FACTORIZATION BIT-SEMANTICS VERSION (bead y4pt): bump on ANY change
/// to qr/cholesky/lu/svd operation order or rounding that can move
/// downstream bits. Pinned by golden-couplings.json;
/// `cargo run -p xtask -- check-goldens` fails on drift until every
/// dependent golden is deliberately re-frozen.
pub const FACTOR_BIT_SEMANTICS_VERSION: u32 = 1;
use crate::gemm::gemm_f64;

/// Blocking width for panel factorizations (pre-autotuner default; bit-
/// relevant only through GEMM's KC contract, which NB does not reach).
const NB: usize = 32;

/// Compute a row-major matrix extent without allowing release-mode wrapping to
/// turn an impossible shape into an apparently valid slice length.
#[track_caller]
fn checked_matrix_entries(rows: usize, cols: usize, operation: &str) -> usize {
    rows.checked_mul(cols).unwrap_or_else(|| {
        panic!("{operation} matrix entry count overflows usize: {rows} * {cols}")
    })
}

/// Admit the complete storage proposition shared by the dense factorization
/// entry points before they read the matrix or allocate factor state.
#[track_caller]
fn admit_matrix_storage(a: &[f64], rows: usize, cols: usize, operation: &str) -> usize {
    let entries = checked_matrix_entries(rows, cols, operation);
    assert_eq!(
        a.len(),
        entries,
        "{operation} storage length mismatch: expected {rows} * {cols} = {entries} entries"
    );
    entries
}

/// Typed factorization failure: which structural assumption broke, where.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FactorError {
    /// A Cholesky pivot was ≤ 0 (matrix not positive definite) at `index`.
    NotSpd {
        /// Zero-based diagonal index of the failing pivot.
        index: usize,
    },
    /// An LU pivot column was exactly zero at elimination step `index`.
    Singular {
        /// Zero-based elimination step with no usable pivot.
        index: usize,
    },
}

impl core::fmt::Display for FactorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FactorError::NotSpd { index } => {
                write!(f, "matrix is not positive definite (pivot {index} <= 0)")
            }
            FactorError::Singular { index } => {
                write!(f, "matrix is singular (no nonzero pivot at step {index})")
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Cholesky
// ---------------------------------------------------------------------------

/// Lower Cholesky factor of an SPD matrix (row-major n×n): A = L·Lᵀ.
#[derive(Debug, Clone)]
pub struct Cholesky {
    n: usize,
    /// Row-major; lower triangle holds L, strict upper is untouched input.
    l: Vec<f64>,
}

/// Blocked right-looking Cholesky. Reads only the lower triangle of `a`.
///
/// # Errors
/// [`FactorError::NotSpd`] with the failing diagonal index.
pub fn cholesky(a: &[f64], n: usize) -> Result<Cholesky, FactorError> {
    admit_matrix_storage(a, n, n, "cholesky");
    let mut l = a.to_vec();
    let mut j = 0;
    while j < n {
        let nb = NB.min(n - j);
        // Unblocked panel factor on the diagonal block + rows below.
        for jj in j..j + nb {
            // d = a_jj − Σ_{k<jj, k≥j} l_jk² (columns < j already folded in
            // by earlier trailing updates).
            let mut d = l[jj * n + jj];
            for k in j..jj {
                d = (-l[jj * n + k]).mul_add(l[jj * n + k], d);
            }
            if d <= 0.0 || !d.is_finite() {
                return Err(FactorError::NotSpd { index: jj });
            }
            let piv = d.sqrt();
            l[jj * n + jj] = piv;
            for i in jj + 1..n {
                let mut v = l[i * n + jj];
                for k in j..jj {
                    v = (-l[i * n + k]).mul_add(l[jj * n + k], v);
                }
                l[i * n + jj] = v / piv;
            }
        }
        // Trailing update: A22 −= L21·L21ᵀ via GEMM (β = 1). Only the lower
        // triangle is meaningful; the full-block update also writes the
        // upper copy, which subsequent panels never read below their own
        // column range... they DO read their diagonal-block upper? No: the
        // panel loop reads l[i*n+k] with k < column — lower only. Safe.
        let rest = n - (j + nb);
        if rest > 0 {
            // L21 is rest×nb at rows j+nb.., cols j..j+nb — extract packed.
            let mut l21 = vec![0.0f64; rest * nb];
            for r in 0..rest {
                for cc in 0..nb {
                    l21[r * nb + cc] = l[(j + nb + r) * n + j + cc];
                }
            }
            let mut l21t = vec![0.0f64; nb * rest];
            for r in 0..rest {
                for cc in 0..nb {
                    l21t[cc * rest + r] = l21[r * nb + cc];
                }
            }
            let mut upd = vec![0.0f64; rest * rest];
            gemm_f64(rest, rest, nb, 1.0, &l21, &l21t, 0.0, &mut upd);
            for r in 0..rest {
                for cc in 0..=r {
                    l[(j + nb + r) * n + j + nb + cc] -= upd[r * rest + cc];
                }
            }
        }
        j += nb;
    }
    Ok(Cholesky { n, l })
}

impl Cholesky {
    /// L entry (lower triangle; upper queries return 0).
    #[must_use]
    pub fn l(&self, i: usize, j: usize) -> f64 {
        if j <= i { self.l[i * self.n + j] } else { 0.0 }
    }

    /// Solve A·x = b in place (forward + back substitution) — the
    /// refinement hook fs-la-mixed-precision consumes.
    pub fn solve(&self, b: &mut [f64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "b must have length {n}");
        for i in 0..n {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(i) {
                v = (-self.l[i * n + k]).mul_add(bk, v);
            }
            b[i] = v / self.l[i * n + i];
        }
        for i in (0..n).rev() {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(n).skip(i + 1) {
                v = (-self.l[k * n + i]).mul_add(bk, v);
            }
            b[i] = v / self.l[i * n + i];
        }
    }
}

// ---------------------------------------------------------------------------
// LU with partial pivoting
// ---------------------------------------------------------------------------

/// P·A = L·U, packed: unit-lower L below the diagonal, U on/above.
#[derive(Debug, Clone)]
pub struct Lu {
    n: usize,
    data: Vec<f64>,
    /// Row permutation: factored row i came from original row `perm[i]`.
    perm: Vec<usize>,
    /// Growth statistic max|U| / max|A| (pivot-growth ledger entry).
    pub growth: f64,
}

/// Partial-pivot LU, deterministic tie-break (lowest row index wins ties
/// in |pivot| — P2). Blocked with GEMM trailing updates.
///
/// # Errors
/// [`FactorError::Singular`] if an entire pivot column is zero.
pub fn lu(a: &[f64], n: usize) -> Result<Lu, FactorError> {
    admit_matrix_storage(a, n, n, "lu");
    let mut m = a.to_vec();
    let mut perm: Vec<usize> = (0..n).collect();
    let max_a = m.iter().fold(0.0f64, |acc, &v| acc.max(v.abs()));
    let mut k0 = 0;
    while k0 < n {
        let nb = NB.min(n - k0);
        // Unblocked panel factor over columns k0..k0+nb, pivoting whole rows.
        for k in k0..k0 + nb {
            // Pivot search: strict > keeps the LOWEST index on ties.
            let mut piv_row = k;
            let mut piv_val = m[k * n + k].abs();
            for r in k + 1..n {
                let v = m[r * n + k].abs();
                if v > piv_val {
                    piv_val = v;
                    piv_row = r;
                }
            }
            if piv_val == 0.0 {
                return Err(FactorError::Singular { index: k });
            }
            if piv_row != k {
                for c in 0..n {
                    m.swap(k * n + c, piv_row * n + c);
                }
                perm.swap(k, piv_row);
            }
            let piv = m[k * n + k];
            for r in k + 1..n {
                let mult = m[r * n + k] / piv;
                m[r * n + k] = mult;
                // Update only within the panel here; the trailing block is
                // GEMM-updated below.
                for c in k + 1..k0 + nb {
                    m[r * n + c] = (-mult).mul_add(m[k * n + c], m[r * n + c]);
                }
            }
        }
        // Block triangular solve: U12 = L11⁻¹·A12 (unit-lower forward).
        let rest = n - (k0 + nb);
        if rest > 0 {
            for c in k0 + nb..n {
                for i in k0..k0 + nb {
                    let mut v = m[i * n + c];
                    for k in k0..i {
                        v = (-m[i * n + k]).mul_add(m[k * n + c], v);
                    }
                    m[i * n + c] = v;
                }
            }
            // Trailing update A22 −= L21·U12 via GEMM.
            let mut l21 = vec![0.0f64; rest * nb];
            for r in 0..rest {
                for cc in 0..nb {
                    l21[r * nb + cc] = m[(k0 + nb + r) * n + k0 + cc];
                }
            }
            let mut u12 = vec![0.0f64; nb * rest];
            for rr in 0..nb {
                for cc in 0..rest {
                    u12[rr * rest + cc] = m[(k0 + rr) * n + k0 + nb + cc];
                }
            }
            let mut upd = vec![0.0f64; rest * rest];
            gemm_f64(rest, rest, nb, 1.0, &l21, &u12, 0.0, &mut upd);
            for r in 0..rest {
                for cc in 0..rest {
                    m[(k0 + nb + r) * n + k0 + nb + cc] -= upd[r * rest + cc];
                }
            }
        }
        k0 += nb;
    }
    let max_u = (0..n)
        .flat_map(|i| (i..n).map(move |j| (i, j)))
        .fold(0.0f64, |acc, (i, j)| acc.max(m[i * n + j].abs()));
    let growth = if max_a > 0.0 { max_u / max_a } else { 1.0 };
    Ok(Lu {
        n,
        data: m,
        perm,
        growth,
    })
}

impl Lu {
    /// The row permutation (factored row i ← original row `perm()[i]`).
    #[must_use]
    pub fn perm(&self) -> &[usize] {
        &self.perm
    }

    /// Solve A·x = b (applies P, forward, back).
    pub fn solve(&self, b: &mut [f64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "b must have length {n}");
        let pb: Vec<f64> = self.perm.iter().map(|&p| b[p]).collect();
        b.copy_from_slice(&pb);
        for i in 0..n {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(i) {
                v = (-self.data[i * n + k]).mul_add(bk, v);
            }
            b[i] = v;
        }
        for i in (0..n).rev() {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(n).skip(i + 1) {
                v = (-self.data[i * n + k]).mul_add(bk, v);
            }
            b[i] = v / self.data[i * n + i];
        }
    }

    /// Solve Aᵀ·x = b (needed by the 1-norm condition estimator).
    pub fn solve_transpose(&self, b: &mut [f64]) {
        let n = self.n;
        assert_eq!(b.len(), n, "b must have length {n}");
        // Aᵀ = Uᵀ·Lᵀ·P: forward with Uᵀ, back with Lᵀ, then P⁻¹.
        for i in 0..n {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(i) {
                v = (-self.data[k * n + i]).mul_add(bk, v);
            }
            b[i] = v / self.data[i * n + i];
        }
        for i in (0..n).rev() {
            let mut v = b[i];
            for (k, &bk) in b.iter().enumerate().take(n).skip(i + 1) {
                v = (-self.data[k * n + i]).mul_add(bk, v);
            }
            b[i] = v;
        }
        let mut out = vec![0.0; n];
        for (&p, &bv) in self.perm.iter().zip(b.iter()) {
            out[p] = bv;
        }
        b.copy_from_slice(&out);
    }

    /// 1-norm condition estimate κ₁ ≈ ‖A‖₁·‖A⁻¹‖₁ via the Hager/Higham
    /// power method on A⁻¹ (deterministic start; a lower bound in theory,
    /// usually within a small factor).
    #[must_use]
    pub fn condition_1(&self, a: &[f64]) -> f64 {
        let n = self.n;
        let norm_a = (0..n)
            .map(|j| (0..n).map(|i| a[i * n + j].abs()).sum::<f64>())
            .fold(0.0f64, f64::max);
        let mut x = vec![1.0 / n as f64; n];
        let mut est = 0.0f64;
        for _ in 0..5 {
            self.solve(&mut x);
            let new_est: f64 = x.iter().map(|v| v.abs()).sum();
            let z: Vec<f64> = x
                .iter()
                .map(|&v| if v >= 0.0 { 1.0 } else { -1.0 })
                .collect();
            let mut zt = z;
            self.solve_transpose(&mut zt);
            // Pick the max-magnitude unit vector (lowest index on ties).
            let (mut best, mut best_i) = (0.0f64, 0usize);
            for (i, &v) in zt.iter().enumerate() {
                if v.abs() > best {
                    best = v.abs();
                    best_i = i;
                }
            }
            if new_est <= est {
                break;
            }
            est = new_est;
            x = vec![0.0; n];
            x[best_i] = 1.0;
        }
        est * norm_a
    }
}

// ---------------------------------------------------------------------------
// Householder QR (blocked, compact-WY)
// ---------------------------------------------------------------------------

/// A = Q·R, m×n with m ≥ n: R in the upper triangle, Householder vectors
/// (unit leading 1 implicit) below the diagonal, scaling factors in `tau`.
#[derive(Debug, Clone)]
pub struct Qr {
    m: usize,
    n: usize,
    data: Vec<f64>,
    tau: Vec<f64>,
}

/// Blocked Householder QR with compact-WY trailing application.
#[must_use]
pub fn qr(a: &[f64], m: usize, n: usize) -> Qr {
    assert!(m >= n, "qr requires m >= n (got {m} x {n})");
    admit_matrix_storage(a, m, n, "qr");
    let mut q = a.to_vec();
    let mut tau = vec![0.0f64; n];
    let mut j0 = 0;
    while j0 < n {
        let nb = NB.min(n - j0);
        // Unblocked panel: columns j0..j0+nb.
        for j in j0..j0 + nb {
            // Householder for column j at rows j..m.
            let mut norm_sq = 0.0f64;
            for i in j..m {
                norm_sq = q[i * n + j].mul_add(q[i * n + j], norm_sq);
            }
            let alpha = q[j * n + j];
            let norm = norm_sq.sqrt();
            if norm == 0.0 {
                tau[j] = 0.0;
                continue;
            }
            let beta = if alpha >= 0.0 { -norm } else { norm };
            let v0 = alpha - beta;
            tau[j] = -v0 / beta; // = (beta - alpha)/beta, in (0, 2]
            let inv_v0 = 1.0 / v0;
            for i in j + 1..m {
                q[i * n + j] *= inv_v0;
            }
            q[j * n + j] = beta;
            // Apply (I − τ v vᵀ) to the remaining panel columns.
            for c in j + 1..j0 + nb {
                let mut dot = q[j * n + c];
                for i in j + 1..m {
                    dot = q[i * n + j].mul_add(q[i * n + c], dot);
                }
                let t = tau[j] * dot;
                q[j * n + c] -= t;
                for i in j + 1..m {
                    q[i * n + c] = (-t).mul_add(q[i * n + j], q[i * n + c]);
                }
            }
        }
        // Trailing block: apply the panel reflectors column-sequentially
        // (v1: sequential application IS the compact-WY product evaluated
        // factor by factor — bit-simpler, same math; the fused WY GEMM form
        // joins the perf lane).
        for j in j0..j0 + nb {
            if tau[j] == 0.0 {
                continue;
            }
            for c in j0 + nb..n {
                let mut dot = q[j * n + c];
                for i in j + 1..m {
                    dot = q[i * n + j].mul_add(q[i * n + c], dot);
                }
                let t = tau[j] * dot;
                q[j * n + c] -= t;
                for i in j + 1..m {
                    q[i * n + c] = (-t).mul_add(q[i * n + j], q[i * n + c]);
                }
            }
        }
        j0 += nb;
    }
    Qr { m, n, data: q, tau }
}

impl Qr {
    /// R entry (upper triangle).
    #[must_use]
    pub fn r(&self, i: usize, j: usize) -> f64 {
        if i <= j && i < self.n {
            self.data[i * self.n + j]
        } else {
            0.0
        }
    }

    /// Apply Qᵀ to a length-m vector in place (reflectors in order).
    pub fn apply_qt(&self, y: &mut [f64]) {
        assert_eq!(y.len(), self.m, "y must have length {}", self.m);
        for j in 0..self.n {
            if self.tau[j] == 0.0 {
                continue;
            }
            let mut dot = y[j];
            for (i, &yi) in y.iter().enumerate().take(self.m).skip(j + 1) {
                dot = self.data[i * self.n + j].mul_add(yi, dot);
            }
            let t = self.tau[j] * dot;
            y[j] -= t;
            for (i, ys) in y.iter_mut().enumerate().take(self.m).skip(j + 1) {
                *ys = (-t).mul_add(self.data[i * self.n + j], *ys);
            }
        }
    }

    /// Apply Q to a length-m vector in place (reflectors reversed).
    pub fn apply_q(&self, y: &mut [f64]) {
        assert_eq!(y.len(), self.m, "y must have length {}", self.m);
        for j in (0..self.n).rev() {
            if self.tau[j] == 0.0 {
                continue;
            }
            let mut dot = y[j];
            for (i, &yi) in y.iter().enumerate().take(self.m).skip(j + 1) {
                dot = self.data[i * self.n + j].mul_add(yi, dot);
            }
            let t = self.tau[j] * dot;
            y[j] -= t;
            for (i, ys) in y.iter_mut().enumerate().take(self.m).skip(j + 1) {
                *ys = (-t).mul_add(self.data[i * self.n + j], *ys);
            }
        }
    }

    /// Least-squares solve min‖A·x − b‖₂ (x returned; b length m).
    #[must_use]
    pub fn solve_ls(&self, b: &[f64]) -> Vec<f64> {
        let mut y = b.to_vec();
        self.apply_qt(&mut y);
        let n = self.n;
        let mut x = y[..n].to_vec();
        for i in (0..n).rev() {
            let mut v = x[i];
            for (k, &xk) in x.iter().enumerate().take(n).skip(i + 1) {
                v = (-self.data[i * n + k]).mul_add(xk, v);
            }
            x[i] = v / self.data[i * n + i];
        }
        x
    }
}

// ---------------------------------------------------------------------------
// TSQR — deterministic tree QR for tall-skinny blocks
// ---------------------------------------------------------------------------

/// TSQR: the R factor of a tall-skinny A (m×n, m ≫ n), computed by leaf
/// QRs over fixed row blocks and pairwise R-combines up a binary tree
/// whose SHAPE is a pure function of the block count — the reduction-tree
/// doctrine applied to orthogonalization. Strictly negative diagonal signs are
/// normalized by row flips. For a finite, full-column-rank matrix, the
/// exact-arithmetic QR theorem makes the resulting positive-diagonal R unique;
/// the finite full-rank regression battery checks agreement across the covered
/// tree shapes within tolerance. Rank-deficient inputs retain deterministic
/// fixed-tree behavior, but this implementation makes no uniqueness or
/// tree-independence claim for their R factors.
#[must_use]
pub fn tsqr_r(a: &[f64], m: usize, n: usize, row_block: usize) -> Vec<f64> {
    assert!(m >= n, "tsqr requires m >= n");
    assert!(row_block >= n, "row_block must be >= n");
    admit_matrix_storage(a, m, n, "tsqr");
    if n == 0 {
        return Vec::new();
    }
    let r_entries = checked_matrix_entries(n, n, "tsqr R");
    // Leaf QRs. Block boundaries are a pure function of (m, row_block, n):
    // a final fragment shorter than n rows is absorbed into the previous
    // block (a leaf QR needs rows >= n).
    let mut bounds: Vec<usize> = (0..m).step_by(row_block).collect();
    bounds.push(m);
    if bounds.len() > 2 && m - bounds[bounds.len() - 2] < n {
        bounds.remove(bounds.len() - 2);
    }
    let mut rs: Vec<Vec<f64>> = Vec::new();
    for w in bounds.windows(2) {
        let (r0, r1) = (w[0], w[1]);
        let start = checked_matrix_entries(r0, n, "tsqr leaf prefix");
        let end = checked_matrix_entries(r1, n, "tsqr leaf prefix");
        rs.push(r_of(&a[start..end], r1 - r0, n));
    }
    // Fixed binary combine tree over block indices: (0,1), (2,3), …
    while rs.len() > 1 {
        let mut next = Vec::with_capacity(rs.len().div_ceil(2));
        for pair in rs.chunks(2) {
            if pair.len() == 2 {
                let stacked_rows = n
                    .checked_mul(2)
                    .expect("tsqr stacked-row count overflows usize");
                let stacked_entries = r_entries
                    .checked_mul(2)
                    .expect("tsqr stacked-R entry count overflows usize");
                let mut stacked = vec![0.0f64; stacked_entries];
                stacked[..r_entries].copy_from_slice(&pair[0]);
                stacked[r_entries..].copy_from_slice(&pair[1]);
                next.push(r_of(&stacked, stacked_rows, n));
            } else {
                next.push(pair[0].clone());
            }
        }
        rs = next;
    }
    rs.pop().unwrap()
}

/// The strictly-negative-sign-normalized R of a QR factorization, returned
/// dense n×n. Zero diagonal sign bits are retained as produced by `qr`.
fn r_of(block: &[f64], rows: usize, n: usize) -> Vec<f64> {
    let f = qr(block, rows, n);
    let mut r = vec![0.0f64; checked_matrix_entries(n, n, "tsqr R")];
    for i in 0..n {
        let flip = if f.r(i, i) < 0.0 { -1.0 } else { 1.0 };
        for j in i..n {
            r[i * n + j] = flip * f.r(i, j);
        }
    }
    r
}

// ---------------------------------------------------------------------------
// One-sided Jacobi SVD
// ---------------------------------------------------------------------------

/// Thin SVD A = U·diag(σ)·Vᵀ (m×n, m ≥ n): superior RELATIVE accuracy for
/// small singular values (the certificate rank-decision property).
#[derive(Debug, Clone)]
pub struct Svd {
    /// Left vectors, m×n row-major.
    pub u: Vec<f64>,
    /// Singular values, descending.
    pub sigma: Vec<f64>,
    /// Right vectors, n×n row-major (columns are the right singular vecs).
    pub v: Vec<f64>,
}

/// One-sided Jacobi with cyclic-by-rows sweeps and a deterministic
/// rotation order. Zero columns yield zero singular values (sorted last).
#[must_use]
pub fn svd_jacobi(a: &[f64], m: usize, n: usize) -> Svd {
    assert!(m >= n, "svd requires m >= n");
    let matrix_entries = admit_matrix_storage(a, m, n, "svd");
    let right_entries = checked_matrix_entries(n, n, "svd right-vector matrix");
    let mut u = a.to_vec();
    let mut v = vec![0.0f64; right_entries];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    let tol = 1e-15;
    for _sweep in 0..60 {
        let mut rotated = false;
        for p in 0..n {
            for q in p + 1..n {
                let (mut app, mut aqq, mut apq) = (0.0f64, 0.0f64, 0.0f64);
                for i in 0..m {
                    let (x, y) = (u[i * n + p], u[i * n + q]);
                    app = x.mul_add(x, app);
                    aqq = y.mul_add(y, aqq);
                    apq = x.mul_add(y, apq);
                }
                if app == 0.0 || aqq == 0.0 {
                    continue;
                }
                if apq.abs() <= tol * (app * aqq).sqrt() {
                    continue;
                }
                rotated = true;
                // Rutishauser rotation.
                let zeta = (aqq - app) / (2.0 * apq);
                let t = if zeta >= 0.0 {
                    1.0 / (zeta + (1.0 + zeta * zeta).sqrt())
                } else {
                    -1.0 / (-zeta + (1.0 + zeta * zeta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let (x, y) = (u[i * n + p], u[i * n + q]);
                    u[i * n + p] = c.mul_add(x, -(s * y));
                    u[i * n + q] = s.mul_add(x, c * y);
                }
                for i in 0..n {
                    let (x, y) = (v[i * n + p], v[i * n + q]);
                    v[i * n + p] = c.mul_add(x, -(s * y));
                    v[i * n + q] = s.mul_add(x, c * y);
                }
            }
        }
        if !rotated {
            break;
        }
    }
    // Column norms → σ; normalize U columns; sort descending (stable:
    // lowest original index first on ties — deterministic).
    let mut order: Vec<usize> = (0..n).collect();
    let mut sig = vec![0.0f64; n];
    for (j, s) in sig.iter_mut().enumerate() {
        let mut nrm = 0.0f64;
        for i in 0..m {
            nrm = u[i * n + j].mul_add(u[i * n + j], nrm);
        }
        *s = nrm.sqrt();
    }
    order.sort_by(|&x, &y| sig[y].partial_cmp(&sig[x]).unwrap().then(x.cmp(&y)));
    let mut u_out = vec![0.0f64; matrix_entries];
    let mut v_out = vec![0.0f64; right_entries];
    let mut s_out = vec![0.0f64; n];
    for (slot, &src) in order.iter().enumerate() {
        s_out[slot] = sig[src];
        let inv = if sig[src] > 0.0 { 1.0 / sig[src] } else { 0.0 };
        for i in 0..m {
            u_out[i * n + slot] = u[i * n + src] * inv;
        }
        for i in 0..n {
            v_out[i * n + slot] = v[i * n + src];
        }
    }
    Svd {
        u: u_out,
        sigma: s_out,
        v: v_out,
    }
}

//! Complex nonsymmetric eigenvalues (bead urvw item 1): Givens-based
//! Hessenberg reduction + explicitly shifted QR iteration with Wilkinson
//! shifts and standard deflation — the general-spectrum solver the
//! Orr–Sommerfeld stack, colleague-matrix rootfinding, and stability
//! analysis all build on.
//!
//! v1 design decisions (recorded): EXPLICIT shifted QR via complex Givens
//! sweeps (numerically sound and simple to verify; the implicit-Q bulge
//! chase is a perf refinement), no balancing pass (moderate-condition
//! inputs; recorded no-claim), eigenvalues only (eigenvectors via inverse
//! iteration are follow-up). Deterministic: fixed sweep orders, exact
//! deflation thresholds, deterministic shift tie-breaks, canonical
//! (re, im)-lexicographic output ordering.

use fs_math::c64::C64;

mod cancellable;
pub use cancellable::{
    LuComplexError, lu_complex_blocked_with_cancel, lu_complex_blocked_workspace_bytes,
    lu_complex_with_cancel,
};

mod svd;
pub use svd::{ComplexSvd, ComplexSvdError, svd_complex_with_cancel};

/// Refusal of the explicitly bounded, cancellable general QR owner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigError {
    Invalid { what: &'static str },
    NotConverged(EigFailure),
    Cancelled,
}

impl core::fmt::Display for EigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complex eigensolve refused: {self:?}")
    }
}
impl std::error::Error for EigError {}

fn poll_eig(cancelled: &mut impl FnMut() -> bool) -> Result<(), EigError> {
    if cancelled() {
        Err(EigError::Cancelled)
    } else {
        Ok(())
    }
}

/// Typed convergence failure: which eigenvalue window stalled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EigFailure {
    /// Active window upper index at exhaustion.
    pub window_hi: usize,
}

impl core::fmt::Display for EigFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "QR iteration failed to converge (window ending {})",
            self.window_hi
        )
    }
}

/// Complex Givens pair (c real, s complex) zeroing `b` against `a`:
/// G·[a; b] = [r; 0] with c² + |s|² = 1.
fn givens(a: C64, b: C64) -> (f64, C64) {
    let bm = b.abs();
    if bm == 0.0 {
        return (1.0, C64::ZERO);
    }
    let am = a.abs();
    if am == 0.0 {
        // Pure swap-with-phase: r takes b's magnitude.
        return (0.0, (b.conj()).scale(1.0 / bm));
    }
    let t = fs_math::det::sqrt(am.mul_add(am, bm * bm));
    let c = am / t;
    // s = (a/|a|)·conj(b)/t — the phase choice that makes r = (a/|a|)·t.
    let phase = a.scale(1.0 / am);
    let s = phase * b.conj().scale(1.0 / t);
    (c, s)
}

/// Apply G to rows (i, j) of `m` over columns [col_lo, n):
/// row_i ← c·row_i + s·row_j; row_j ← −conj(s)·row_i_old + c·row_j.
fn rot_rows(m: &mut [C64], n: usize, i: usize, j: usize, c: f64, s: C64, col_lo: usize) {
    for k in col_lo..n {
        let a = m[i * n + k];
        let b = m[j * n + k];
        m[i * n + k] = a.scale(c) + s * b;
        m[j * n + k] = b.scale(c) - s.conj() * a;
    }
}

/// Apply Gᴴ to columns (i, j) over rows [0, row_hi).
fn rot_cols(m: &mut [C64], n: usize, i: usize, j: usize, c: f64, s: C64, row_hi: usize) {
    for k in 0..row_hi {
        let a = m[k * n + i];
        let b = m[k * n + j];
        m[k * n + i] = a.scale(c) + b * s.conj();
        m[k * n + j] = b.scale(c) - a * s;
    }
}

/// Reduce to upper Hessenberg by Givens similarity transforms
/// (deterministic column-by-column, bottom-up zeroing).
fn hessenberg(
    m: &mut [C64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), EigError> {
    for k in 0..n.saturating_sub(2) {
        for i in ((k + 2)..n).rev() {
            poll_eig(cancelled)?;
            let a = m[(i - 1) * n + k];
            let b = m[i * n + k];
            if b.abs() == 0.0 {
                continue;
            }
            let (c, s) = givens(a, b);
            rot_rows(m, n, i - 1, i, c, s, k);
            rot_cols(m, n, i - 1, i, c, s, n);
        }
    }
    Ok(())
}

/// Eigenvalues of the trailing 2×2 [[a,b],[c,d]] via the complex
/// quadratic formula; deterministic branch selection.
fn eig2(a: C64, b: C64, c: C64, d: C64) -> (C64, C64) {
    let t = (a + d).scale(0.5);
    let disc = ((a - d).scale(0.5) * (a - d).scale(0.5) + b * c).sqrt();
    (t + disc, t - disc)
}

/// All eigenvalues of a dense complex matrix (row-major n×n), sorted
/// canonically by (re, im). Deterministic across reruns and ISAs.
///
/// # Errors
/// [`EigFailure`] if a QR window exhausts its iteration cap (rare;
/// exceptional shifts are applied before giving up).
pub fn eig(a: &[C64], n: usize) -> Result<Vec<C64>, EigFailure> {
    assert_eq!(a.len(), n * n, "a must be n*n = {}", n * n);
    eig_with_cancel(a, n, 120, || false).map_err(|error| match error {
        EigError::NotConverged(failure) => failure,
        EigError::Invalid { .. } => EigFailure {
            window_hi: n.saturating_sub(1),
        },
        EigError::Cancelled => unreachable!("cancellation is disabled"),
    })
}

/// Same Hessenberg/QR owner as [`eig`], with an explicit per-window sweep cap
/// and cancellation at elimination/rotation boundaries. No partial spectrum
/// is returned. Callers admitting large matrices must also budget dense work.
pub fn eig_with_cancel(
    a: &[C64],
    n: usize,
    max_sweeps_per_window: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<C64>, EigError> {
    poll_eig(&mut cancelled)?;
    if n.checked_mul(n) != Some(a.len())
        || max_sweeps_per_window == 0
        || a.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(EigError::Invalid {
            what: "finite square matrix and positive sweep budget required",
        });
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut h = a.to_vec();
    hessenberg(&mut h, n, &mut cancelled)?;
    let mut eigs: Vec<C64> = Vec::with_capacity(n);
    let mut hi = n - 1; // active window is rows/cols [0, hi]
    let mut iters_here = 0usize;
    loop {
        poll_eig(&mut cancelled)?;
        // Deflate any negligible subdiagonals inside [0, hi].
        let mut k = hi;
        while k > 0 {
            let sub = h[k * n + (k - 1)].abs();
            let scale = h[k * n + k].abs() + h[(k - 1) * n + (k - 1)].abs();
            if sub <= 1e-15 * scale.max(f64::MIN_POSITIVE) {
                h[k * n + (k - 1)] = C64::ZERO;
            }
            k -= 1;
        }
        // Peel converged eigenvalues from the bottom of the window.
        loop {
            if hi == 0 {
                eigs.push(h[0]);
                eigs.sort_by(|x, y| x.re.total_cmp(&y.re).then_with(|| x.im.total_cmp(&y.im)));
                return Ok(eigs);
            }
            if h[hi * n + (hi - 1)].abs() == 0.0 {
                eigs.push(h[hi * n + hi]);
                hi -= 1;
                iters_here = 0;
                continue;
            }
            if hi >= 1 && (hi == 1 || h[(hi - 1) * n + (hi - 2)].abs() == 0.0) {
                // Isolated 2×2 block: closed form, deflate both.
                let (l1, l2) = eig2(
                    h[(hi - 1) * n + (hi - 1)],
                    h[(hi - 1) * n + hi],
                    h[hi * n + (hi - 1)],
                    h[hi * n + hi],
                );
                eigs.push(l1);
                eigs.push(l2);
                if hi == 1 {
                    eigs.sort_by(|x, y| x.re.total_cmp(&y.re).then_with(|| x.im.total_cmp(&y.im)));
                    return Ok(eigs);
                }
                hi -= 2;
                iters_here = 0;
                continue;
            }
            break;
        }
        // Shift: Wilkinson (trailing 2×2 eigenvalue nearest h[hi][hi]);
        // every 12th iteration an exceptional shift breaks cycles.
        iters_here += 1;
        if iters_here > max_sweeps_per_window {
            return Err(EigError::NotConverged(EigFailure { window_hi: hi }));
        }
        let hnn = h[hi * n + hi];
        let mu = if iters_here.is_multiple_of(12) {
            // Exceptional: perturb by the subdiagonal magnitude.
            hnn + C64::from_re(h[hi * n + (hi - 1)].abs())
        } else {
            let (l1, l2) = eig2(
                h[(hi - 1) * n + (hi - 1)],
                h[(hi - 1) * n + hi],
                h[hi * n + (hi - 1)],
                hnn,
            );
            let (d1, d2) = ((l1 - hnn).abs(), (l2 - hnn).abs());
            // Deterministic tie-break via total_cmp chains: distance to
            // h_nn, then |λ|, then (re, im) lexicographic.
            let pick_l1 = d1
                .total_cmp(&d2)
                .then_with(|| l1.abs().total_cmp(&l2.abs()))
                .then_with(|| l1.re.total_cmp(&l2.re))
                .then_with(|| l1.im.total_cmp(&l2.im))
                .is_le();
            if pick_l1 { l1 } else { l2 }
        };
        // Explicit shifted QR sweep on the window: H − μI = QR (Givens
        // zeroing the subdiagonal), then H ← RQ + μI (apply the same
        // rotations on the right).
        for d in 0..=hi {
            h[d * n + d] = h[d * n + d] - mu;
        }
        let mut rots: Vec<(usize, f64, C64)> = Vec::with_capacity(hi);
        for k in 0..hi {
            poll_eig(&mut cancelled)?;
            let aa = h[k * n + k];
            let bb = h[(k + 1) * n + k];
            let (c, s) = givens(aa, bb);
            rot_rows(&mut h, n, k, k + 1, c, s, k);
            rots.push((k, c, s));
        }
        for &(k, c, s) in &rots {
            poll_eig(&mut cancelled)?;
            rot_cols(&mut h, n, k, k + 1, c, s, (k + 2).min(hi + 1));
        }
        for d in 0..=hi {
            h[d * n + d] = h[d * n + d] + mu;
        }
        if h.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(EigError::Invalid {
                what: "QR arithmetic is not representable",
            });
        }
    }
}

/// Complex LU with partial pivoting (magnitude pivots, lowest index on
/// ties) — the solve kernel generalized-eigenproblem reductions need
/// (B⁻¹·A for the Orr–Sommerfeld stack).
#[derive(Debug, Clone)]
pub struct LuComplex {
    n: usize,
    data: Vec<C64>,
    perm: Vec<usize>,
}

std::thread_local! {
    static LU_WORKERS: core::cell::Cell<usize> = const { core::cell::Cell::new(1) };
}

/// Run `f` with `workers` (at least one) threads declared for every
/// [`lu_complex`] factorization it makes on this thread, then restore the
/// previous declaration, also on unwind. This is a scheduling budget only:
/// factors are bit-identical for every worker count.
pub fn with_lu_workers<R>(workers: usize, f: impl FnOnce() -> R) -> R {
    struct Restore(usize);
    impl Drop for Restore {
        fn drop(&mut self) {
            LU_WORKERS.with(|w| w.set(self.0));
        }
    }
    let _restore = Restore(LU_WORKERS.with(|w| w.replace(workers.max(1))));
    f()
}

/// Factor with the blocked elimination ([`lu_complex_blocked_with_cancel`])
/// on the workers declared by [`with_lu_workers`], one by default. Factors
/// are bit-identical to [`lu_complex_with_cancel`].
/// # Errors: [`EigFailure`] reused as a typed singularity signal
/// (`window_hi` = failing elimination step).
pub fn lu_complex(a: &[C64], n: usize) -> Result<LuComplex, EigFailure> {
    let workers = LU_WORKERS.with(core::cell::Cell::get);
    lu_complex_blocked_with_cancel(a, n, workers, || false).map_err(|error| match error {
        LuComplexError::Singular { index } => EigFailure { window_hi: index },
        LuComplexError::Cancelled => unreachable!("cancellation is disabled"),
    })
}

impl LuComplex {
    /// Solve A·x = b in place.
    pub fn solve(&self, b: &mut [C64]) {
        self.solve_with_cancel(b, || false)
            .expect("cancellation is disabled");
    }
}

/// Determinant via complex Gaussian elimination with partial pivoting
/// (test oracle for eigenvalue products; exposed because a complex det
/// has independent uses).
#[must_use]
pub fn det_complex(a: &[C64], n: usize) -> C64 {
    assert_eq!(a.len(), n * n);
    let mut m = a.to_vec();
    let mut det = C64::ONE;
    for k in 0..n {
        // Partial pivot by magnitude, lowest index on ties.
        let (mut best, mut best_i) = (m[k * n + k].abs(), k);
        for i in k + 1..n {
            let v = m[i * n + k].abs();
            if v > best {
                best = v;
                best_i = i;
            }
        }
        if best == 0.0 {
            return C64::ZERO;
        }
        if best_i != k {
            for c in 0..n {
                m.swap(k * n + c, best_i * n + c);
            }
            det = -det;
        }
        let piv = m[k * n + k];
        det = det * piv;
        for i in k + 1..n {
            let f = m[i * n + k] / piv;
            for c in k..n {
                let sub = f * m[k * n + c];
                m[i * n + c] = m[i * n + c] - sub;
            }
        }
    }
    det
}

//! Randomized NLA (plan §6.1 [F], bead 6ys.7): range finders, randomized
//! SVD, Nyström, sketch-and-precondition least squares, and
//! Hutchinson/Hutch++ trace estimation.
//!
//! REPLAYABILITY IS THE DESIGN CENTER: every random object draws from a
//! KEYED Philox stream (fs-rand's logical-identity doctrine), so a
//! "stochastic" estimate is a pure function of (seed, kernel, tile) —
//! bitwise reproducible, cross-ISA (the distributions themselves are
//! bit-deterministic, proven in fs-rand). Every result carries evidence
//! metadata (probe counts, error estimates, variance estimates) — the
//! Evidence<T> integration point.

use crate::factor::qr;
use fs_rand::{Stream, StreamKey};

/// Kernel ids for stream keying (stable registry — changing one is a
/// golden-evidence event).
const K_RANGE: u32 = 0x5A01;
const K_SKETCH: u32 = 0x5A02;
const K_TRACE: u32 = 0x5A03;

/// Evidence attached to a range approximation.
#[derive(Debug, Clone)]
pub struct RangeReport {
    /// Columns in the returned basis.
    pub rank: usize,
    /// Posterior estimate of ‖A − QQᵀA‖₂ (probabilistic upper indicator;
    /// the G0 battery validates its coverage empirically).
    pub est_error: f64,
    /// Random probes consumed by the estimator.
    pub probes: usize,
}

/// Orthonormal basis (m×k, row-major columns) approximating range(A),
/// with fixed rank + oversampling and `q` power iterations.
/// Deterministic given `seed`.
#[must_use]
pub fn range_finder(
    a: &[f64],
    m: usize,
    n: usize,
    rank: usize,
    oversample: usize,
    q_power: usize,
    seed: u64,
) -> (Vec<f64>, RangeReport) {
    assert_eq!(a.len(), m * n, "a must be m*n");
    let k = (rank + oversample).min(n).min(m);
    // Gaussian test matrix from the keyed stream.
    let mut s = StreamKey { seed, kernel: K_RANGE, tile: 0 }.stream();
    let mut y = vec![0.0f64; m * k]; // Y = A·Ω
    let mut omega_col = vec![0.0f64; n];
    let mut ycol = vec![0.0f64; m];
    for j in 0..k {
        for w in omega_col.iter_mut() {
            *w = s.next_normal();
        }
        matvec(a, m, n, &omega_col, &mut ycol);
        for i in 0..m {
            y[i * k + j] = ycol[i];
        }
    }
    // Power iterations: Y ← A·(Aᵀ·Y), re-orthonormalizing between steps
    // (numerical hygiene for slow spectral decay).
    let mut basis = orthonormal_columns(&y, m, k);
    for _ in 0..q_power {
        let mut z = vec![0.0f64; n * k];
        let mut zc = vec![0.0f64; n];
        let mut bc = vec![0.0f64; m];
        for j in 0..k {
            for i in 0..m {
                bc[i] = basis[i * k + j];
            }
            matvec_t(a, m, n, &bc, &mut zc);
            for i in 0..n {
                z[i * k + j] = zc[i];
            }
        }
        let zq = orthonormal_columns(&z, n, k);
        let mut y2 = vec![0.0f64; m * k];
        for j in 0..k {
            for i in 0..n {
                zc[i] = zq[i * k + j];
            }
            matvec(a, m, n, &zc, &mut bc);
            for i in 0..m {
                y2[i * k + j] = bc[i];
            }
        }
        basis = orthonormal_columns(&y2, m, k);
    }
    // Posterior error estimate: s probe vectors through the residual.
    const PROBES: usize = 8;
    let mut est = 0.0f64;
    let mut probe = vec![0.0f64; n];
    let mut aw = vec![0.0f64; m];
    for _ in 0..PROBES {
        for w in probe.iter_mut() {
            *w = s.next_normal();
        }
        matvec(a, m, n, &probe, &mut aw);
        project_out(&basis, m, k, &mut aw);
        let nrm = fs_math::det::sqrt(aw.iter().map(|t| t * t).sum::<f64>())
            / fs_math::det::sqrt(probe.iter().map(|t| t * t).sum::<f64>().max(f64::MIN_POSITIVE));
        est = est.max(nrm);
    }
    // Scale by the standard Gaussian max factor (10·sqrt(2/π) covers the
    // declared 1e-3-ish failure probability band for 8 probes).
    let report = RangeReport { rank: k, est_error: est * 10.0, probes: PROBES };
    (basis, report)
}

/// Randomized SVD: A ≈ U·diag(σ)·Vᵀ with rank-k factors.
/// Returns (u m×k, sigma k, v n×k, report).
#[must_use]
pub fn rsvd(
    a: &[f64],
    m: usize,
    n: usize,
    rank: usize,
    oversample: usize,
    q_power: usize,
    seed: u64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, RangeReport) {
    let (q, report) = range_finder(a, m, n, rank, oversample, q_power, seed);
    let k = report.rank;
    // B = Qᵀ·A is k×n; SVD via the landed one-sided Jacobi on Bᵀ (n×k,
    // n ≥ k): Bᵀ = W·Σ·Zᵀ  ⇒  B = Z·Σ·Wᵀ  ⇒  U = Q·Z, V = W.
    let mut bt = vec![0.0f64; n * k]; // Bᵀ, row-major n×k
    for j in 0..n {
        for i in 0..k {
            let mut acc = 0.0f64;
            for r in 0..m {
                acc = q[r * k + i].mul_add(a[r * n + j], acc);
            }
            bt[j * k + i] = acc;
        }
    }
    let svd = crate::factor::svd_jacobi(&bt, n, k);
    // U = Q·Z (Z = svd.v, k×k), keep the leading `rank` columns.
    let keep = rank.min(k);
    let mut u = vec![0.0f64; m * keep];
    for i in 0..m {
        for j in 0..keep {
            let mut acc = 0.0f64;
            for l in 0..k {
                acc = q[i * k + l].mul_add(svd.v[l * k + j], acc);
            }
            u[i * keep + j] = acc;
        }
    }
    let mut v = vec![0.0f64; n * keep];
    for i in 0..n {
        for j in 0..keep {
            v[i * keep + j] = svd.u[i * k + j];
        }
    }
    (u, svd.sigma[..keep].to_vec(), v, report)
}

/// Nyström approximation of a PSD matrix: Â = (AΩ)·(ΩᵀAΩ)⁺·(AΩ)ᵀ,
/// returned as a factor F (n×k) with Â = F·Fᵀ. Deterministic per seed.
#[must_use]
pub fn nystrom_psd(a: &[f64], n: usize, rank: usize, oversample: usize, seed: u64) -> Vec<f64> {
    assert_eq!(a.len(), n * n, "a must be n*n");
    let k = (rank + oversample).min(n);
    let mut s = StreamKey { seed, kernel: K_RANGE, tile: 1 }.stream();
    // Y = A·Ω (n×k), C = Ωᵀ·Y (k×k, PSD).
    let mut omega = vec![0.0f64; n * k];
    for w in omega.iter_mut() {
        *w = s.next_normal();
    }
    let mut y = vec![0.0f64; n * k];
    for j in 0..k {
        for i in 0..n {
            let mut acc = 0.0f64;
            for l in 0..n {
                acc = a[i * n + l].mul_add(omega[l * k + j], acc);
            }
            y[i * k + j] = acc;
        }
    }
    let mut c = vec![0.0f64; k * k];
    for i in 0..k {
        for j in 0..k {
            let mut acc = 0.0f64;
            for l in 0..n {
                acc = omega[l * k + i].mul_add(y[l * k + j], acc);
            }
            c[i * k + j] = acc;
        }
    }
    // Symmetrize + eigen; pseudo-inverse square root via jacobi_eigh.
    for i in 0..k {
        for j in i + 1..k {
            let avg = f64::midpoint(c[i * k + j], c[j * k + i]);
            c[i * k + j] = avg;
            c[j * k + i] = avg;
        }
    }
    let (vals, vecs) = crate::eigen::jacobi_eigh(&c, k);
    let vmax = vals.last().copied().unwrap_or(0.0).max(f64::MIN_POSITIVE);
    // F = Y·V·diag(1/sqrt(λ)) over the numerically nonzero spectrum.
    let mut f = vec![0.0f64; n * k];
    for j in 0..k {
        let lam = vals[j];
        if lam <= 1e-12 * vmax {
            continue; // pseudo-inverse: drop the null directions
        }
        let inv_sqrt = 1.0 / fs_math::det::sqrt(lam);
        for i in 0..n {
            let mut acc = 0.0f64;
            for l in 0..k {
                acc = y[i * k + l].mul_add(vecs[l * k + j], acc);
            }
            f[i * k + j] = acc * inv_sqrt;
        }
    }
    f
}

/// Sketch-and-precondition least squares: min ‖Ax − b‖₂ for tall dense A
/// (m ≫ n) via a sparse-sign sketch, QR of the sketch, and preconditioned
/// normal-equation CG. Returns (x, iterations).
#[must_use]
pub fn sketch_ls(a: &[f64], m: usize, n: usize, b: &[f64], seed: u64) -> (Vec<f64>, usize) {
    assert_eq!(a.len(), m * n, "a must be m*n");
    assert_eq!(b.len(), m, "b must be m");
    assert!(m >= 4 * n, "sketch LS wants m >= 4n (tall)");
    let ms = 4 * n; // sketch rows
    // Sparse sign sketch: 8 nonzeros per COLUMN of S (ms×m), ±1/√8,
    // positions/signs from the keyed stream (deterministic).
    let mut s = StreamKey { seed, kernel: K_SKETCH, tile: 0 }.stream();
    let mut sa = vec![0.0f64; ms * n];
    let mut sb = vec![0.0f64; ms];
    let scale = 1.0 / fs_math::det::sqrt(8.0);
    for col in 0..m {
        for _ in 0..8 {
            let row = s.next_below(ms as u64) as usize;
            let sign = if s.next_below(2) == 0 { scale } else { -scale };
            for j in 0..n {
                sa[row * n + j] = sign.mul_add(a[col * n + j], sa[row * n + j]);
            }
            sb[row] = sign.mul_add(b[col], sb[row]);
        }
    }
    // R from QR of the sketch is the preconditioner.
    let f = qr(&sa, ms, n);
    // Preconditioned normal equations: solve (R⁻ᵀAᵀAR⁻¹)·y = R⁻ᵀAᵀb by
    // CG, then x = R⁻¹y. Apply operators matrix-free.
    let rsolve = |v: &mut [f64]| {
        // R⁻¹·v (back substitution on the stored R).
        for i in (0..n).rev() {
            let mut acc = v[i];
            for k in i + 1..n {
                acc = (-f.r(i, k)).mul_add(v[k], acc);
            }
            v[i] = acc / f.r(i, i);
        }
    };
    let rtsolve = |v: &mut [f64]| {
        // R⁻ᵀ·v (forward substitution).
        for i in 0..n {
            let mut acc = v[i];
            for k in 0..i {
                acc = (-f.r(k, i)).mul_add(v[k], acc);
            }
            v[i] = acc / f.r(i, i);
        }
    };
    let op = |v: &[f64], out: &mut [f64]| {
        // out = R⁻ᵀ Aᵀ A R⁻¹ v.
        let mut t = v.to_vec();
        rsolve(&mut t);
        let mut at = vec![0.0f64; m];
        matvec(a, m, n, &t, &mut at);
        let mut ata = vec![0.0f64; n];
        matvec_t(a, m, n, &at, &mut ata);
        rtsolve(&mut ata);
        out.copy_from_slice(&ata);
    };
    let mut rhs = vec![0.0f64; n];
    matvec_t(a, m, n, b, &mut rhs);
    rtsolve(&mut rhs);
    // CG on the preconditioned system (condition ~O(1) by construction).
    let mut x = vec![0.0f64; n];
    let mut r = rhs.clone();
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|t| t * t).sum();
    let rhs_norm = fs_math::det::sqrt(rr).max(f64::MIN_POSITIVE);
    let mut ap = vec![0.0f64; n];
    let mut iters = 0usize;
    for it in 0..200 {
        if fs_math::det::sqrt(rr) <= 1e-13 * rhs_norm {
            iters = it;
            break;
        }
        op(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(u, w)| u * w).sum();
        let alpha = rr / pap;
        for i in 0..n {
            x[i] = alpha.mul_add(p[i], x[i]);
            r[i] = (-alpha).mul_add(ap[i], r[i]);
        }
        let rr_new: f64 = r.iter().map(|t| t * t).sum();
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..n {
            p[i] = beta.mul_add(p[i], r[i]);
        }
        iters = it + 1;
    }
    rsolve(&mut x);
    (x, iters)
}

/// Trace-estimate evidence.
#[derive(Debug, Clone)]
pub struct TraceReport {
    /// The estimate.
    pub estimate: f64,
    /// Probes consumed.
    pub probes: usize,
    /// Sample variance of the per-probe estimates (evidence metadata).
    pub variance_est: f64,
}

/// Hutchinson trace estimator with Rademacher probes from a keyed
/// stream: unbiased, replayable, evidence-carrying.
#[must_use]
pub fn hutchinson(a: &[f64], n: usize, probes: usize, seed: u64) -> TraceReport {
    assert_eq!(a.len(), n * n, "a must be n*n");
    let mut s = StreamKey { seed, kernel: K_TRACE, tile: 0 }.stream();
    let mut z = vec![0.0f64; n];
    let mut az = vec![0.0f64; n];
    let (mut m1, mut m2) = (0.0f64, 0.0f64);
    for _ in 0..probes {
        for zi in z.iter_mut() {
            *zi = if s.next_below(2) == 0 { 1.0 } else { -1.0 };
        }
        matvec(a, n, n, &z, &mut az);
        let est: f64 = z.iter().zip(&az).map(|(u, w)| u * w).sum();
        m1 += est;
        m2 += est * est;
    }
    let p = probes as f64;
    let mean = m1 / p;
    let var = (m2 / p - mean * mean).max(0.0) / p;
    TraceReport { estimate: mean, probes, variance_est: var }
}

/// Hutch++: exact trace on a rangefinder subspace plus Hutchinson on the
/// deflated remainder — provably lower variance at matched probe budgets
/// for decaying spectra (measured in the battery, not just cited).
#[must_use]
pub fn hutch_pp(a: &[f64], n: usize, probes: usize, seed: u64) -> TraceReport {
    assert_eq!(a.len(), n * n, "a must be n*n");
    let k = (probes / 3).max(1).min(n);
    // Subspace capture: Q from A·Ω with k columns.
    let (q, _) = range_finder(a, n, n, k, 0, 0, seed ^ 0x9E37_79B9);
    let kq = q.len() / n;
    // Exact part: tr(Qᵀ A Q).
    let mut exact = 0.0f64;
    let mut aq = vec![0.0f64; n];
    let mut qc = vec![0.0f64; n];
    for j in 0..kq {
        for i in 0..n {
            qc[i] = q[i * kq + j];
        }
        matvec(a, n, n, &qc, &mut aq);
        exact += qc.iter().zip(&aq).map(|(u, w)| u * w).sum::<f64>();
    }
    // Deflated Hutchinson on (I−QQᵀ)A(I−QQᵀ) with the remaining budget.
    let rem = probes.saturating_sub(kq).max(1);
    let mut s = StreamKey { seed, kernel: K_TRACE, tile: 1 }.stream();
    let mut z = vec![0.0f64; n];
    let mut az = vec![0.0f64; n];
    let (mut m1, mut m2) = (0.0f64, 0.0f64);
    for _ in 0..rem {
        for zi in z.iter_mut() {
            *zi = if s.next_below(2) == 0 { 1.0 } else { -1.0 };
        }
        project_out(&q, n, kq, &mut z);
        matvec(a, n, n, &z, &mut az);
        project_out(&q, n, kq, &mut az);
        let est: f64 = z.iter().zip(&az).map(|(u, w)| u * w).sum();
        m1 += est;
        m2 += est * est;
    }
    let p = rem as f64;
    let mean = m1 / p;
    let var = (m2 / p - mean * mean).max(0.0) / p;
    TraceReport { estimate: exact + mean, probes: kq + rem, variance_est: var }
}

// --- helpers ---

fn matvec(a: &[f64], m: usize, n: usize, x: &[f64], out: &mut [f64]) {
    for i in 0..m {
        let mut acc = 0.0f64;
        for j in 0..n {
            acc = a[i * n + j].mul_add(x[j], acc);
        }
        out[i] = acc;
    }
}

fn matvec_t(a: &[f64], m: usize, n: usize, x: &[f64], out: &mut [f64]) {
    for j in 0..n {
        out[j] = 0.0;
    }
    for i in 0..m {
        let xi = x[i];
        if xi == 0.0 {
            continue;
        }
        for j in 0..n {
            out[j] = a[i * n + j].mul_add(xi, out[j]);
        }
    }
}

/// Orthonormalize the k columns of a row-major m×k block via QR.
fn orthonormal_columns(block: &[f64], m: usize, k: usize) -> Vec<f64> {
    let f = qr(block, m, k);
    let mut q = vec![0.0f64; m * k];
    for j in 0..k {
        let mut e = vec![0.0f64; m];
        e[j] = 1.0;
        f.apply_q(&mut e);
        for i in 0..m {
            q[i * k + j] = e[i];
        }
    }
    q
}

/// v ← (I − QQᵀ)·v for a row-major m×k orthonormal Q.
fn project_out(q: &[f64], m: usize, k: usize, v: &mut [f64]) {
    for j in 0..k {
        let mut dot = 0.0f64;
        for i in 0..m {
            dot = q[i * k + j].mul_add(v[i], dot);
        }
        for i in 0..m {
            v[i] = (-dot).mul_add(q[i * k + j], v[i]);
        }
    }
}

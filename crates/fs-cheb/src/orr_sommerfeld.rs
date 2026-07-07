//! Orr–Sommerfeld spectral stability (bead urvw item 2): temporal
//! eigenvalues of plane Poiseuille flow U(y) = 1 − y² via Chebyshev
//! collocation with CLAMPED boundary conditions (φ(±1) = φ′(±1) = 0),
//! following the classical Trefethen construction: the clamped fourth-
//! derivative matrix comes from the substitution φ = (1 − y²)·u, i.e.
//! D4c = (diag(1−y²)·D⁴ − 8·diag(y)·D³ − 12·D²)·diag(1/(1−y²)) on the
//! interior nodes, with plain interior D² carrying the Dirichlet part.
//!
//! Operator (temporal form; eigenvalues λ with Re(λ) > 0 unstable):
//!   A = (D4c − 2α²·D2 + α⁴·I)/Re − 2iα·I − iα·diag(U)·(D2 − α²·I)
//!   B = D2 − α²·I,     A·φ = λ·B·φ.
//!
//! THE ACCEPTANCE TEST of the whole spectral stack: the neutral-stability
//! crossing at α = 1.02056 must reproduce the published critical Reynolds
//! number Re_c ≈ 5772.22.
//!
//! "Modal growth rates along a path" (the vessel flagship's first-class
//! query): [`growth_rates`] returns the k rightmost eigenvalues, sorted
//! by descending real part with deterministic tie-breaks.

use crate::{diff_matrix, lobatto_points};
use fs_la::eigen_complex::{EigFailure, eig, lu_complex};
use fs_math::c64::C64;

/// Real n×n matrix product (row-major helpers for the D powers).
fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                out[i * n + j] = aik.mul_add(b[k * n + j], out[i * n + j]);
            }
        }
    }
    out
}

/// The interior-node Orr–Sommerfeld pair (A, B) at collocation order `n`
/// (n+1 Lobatto points, n−1 interior nodes). Returns (a, b, dim).
fn os_matrices(re: f64, alpha: f64, n: usize) -> (Vec<C64>, Vec<C64>, usize) {
    assert!(n >= 8, "collocation order too small for OS");
    let m = n + 1;
    let x = lobatto_points(n);
    let d1 = diff_matrix(n);
    let d2 = matmul(&d1, &d1, m);
    let d3 = matmul(&d2, &d1, m);
    let d4 = matmul(&d3, &d1, m);
    // Clamped D4 on the full grid: (diag(1−x²)·D4 − 8·diag(x)·D3 − 12·D2)
    // · diag(1/(1−x²)) — the boundary columns are never used (interior
    // restriction below), so the 1/(1−x²) singularity at ±1 is moot.
    let ni = n - 1;
    let mut d4c = vec![0.0f64; ni * ni];
    let mut d2i = vec![0.0f64; ni * ni];
    for i in 0..ni {
        let gi = i + 1; // grid row
        let wi = 1.0 - x[gi] * x[gi];
        for j in 0..ni {
            let gj = j + 1; // grid col
            let raw = wi.mul_add(
                d4[gi * m + gj],
                (-8.0 * x[gi]).mul_add(d3[gi * m + gj], -12.0 * d2[gi * m + gj]),
            );
            d4c[i * ni + j] = raw / (1.0 - x[gj] * x[gj]);
            d2i[i * ni + j] = d2[gi * m + gj];
        }
    }
    let a2 = alpha * alpha;
    let mut a = vec![C64::ZERO; ni * ni];
    let mut b = vec![C64::ZERO; ni * ni];
    for i in 0..ni {
        let u = 1.0 - x[i + 1] * x[i + 1]; // U(y) = 1 − y²
        for j in 0..ni {
            let idelta = if i == j { 1.0 } else { 0.0 };
            // (D4c − 2α²·D2 + α⁴)/Re — the viscous block (real).
            let visc =
                (d4c[i * ni + j] - 2.0 * a2 * d2i[i * ni + j] + a2 * a2 * idelta) / re;
            // −2iα − iα·U·(D2 − α²): the inertial block (imaginary).
            let inert = -alpha * (2.0 * idelta + u * (d2i[i * ni + j] - a2 * idelta));
            a[i * ni + j] = C64::new(visc, inert);
            b[i * ni + j] = C64::from_re(d2i[i * ni + j] - a2 * idelta);
        }
    }
    (a, b, ni)
}

/// The k rightmost temporal eigenvalues (descending real part,
/// deterministic tie-breaks) — "modal growth rates σ₁..σ_k" as a
/// first-class query. Re(λ) > 0 means exponential growth.
///
/// # Errors
/// [`EigFailure`] if B is singular or QR exhausts (neither occurs on the
/// tested parameter ranges).
pub fn growth_rates(re: f64, alpha: f64, n: usize, k: usize) -> Result<Vec<C64>, EigFailure> {
    let (a, b, ni) = os_matrices(re, alpha, n);
    // Reduce the generalized problem to standard: M = B⁻¹·A.
    let fact = lu_complex(&b, ni)?;
    let mut m_mat = vec![C64::ZERO; ni * ni];
    let mut col = vec![C64::ZERO; ni];
    for j in 0..ni {
        for i in 0..ni {
            col[i] = a[i * ni + j];
        }
        fact.solve(&mut col);
        for i in 0..ni {
            m_mat[i * ni + j] = col[i];
        }
    }
    let mut eigs = eig(&m_mat, ni)?;
    eigs.sort_by(|p, q| {
        q.re.total_cmp(&p.re).then_with(|| q.im.total_cmp(&p.im))
    });
    eigs.truncate(k);
    Ok(eigs)
}

/// The rightmost growth rate max Re(λ) at (Re, α).
///
/// # Errors
/// Propagates [`EigFailure`] from the eigen solve.
pub fn max_growth(re: f64, alpha: f64, n: usize) -> Result<f64, EigFailure> {
    Ok(growth_rates(re, alpha, n, 1)?[0].re)
}

/// Critical Reynolds number at fixed α by deterministic bisection on the
/// sign of the rightmost growth rate. `lo`/`hi` must bracket the
/// crossing (checked, structured panic otherwise).
///
/// # Errors
/// Propagates [`EigFailure`] from the eigen solves.
pub fn critical_reynolds(
    alpha: f64,
    n: usize,
    mut lo: f64,
    mut hi: f64,
) -> Result<f64, EigFailure> {
    let glo = max_growth(lo, alpha, n)?;
    let ghi = max_growth(hi, alpha, n)?;
    assert!(
        glo < 0.0 && ghi > 0.0,
        "bracket must straddle neutral stability: growth({lo}) = {glo}, growth({hi}) = {ghi}"
    );
    for _ in 0..40 {
        let mid = f64::midpoint(lo, hi);
        if max_growth(mid, alpha, n)? > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < 0.005 {
            break;
        }
    }
    Ok(f64::midpoint(lo, hi))
}

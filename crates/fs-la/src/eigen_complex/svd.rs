//! Bounded one-sided complex Jacobi SVD.
//!
//! The working columns are orthogonalized directly, without forming A^H A.
//! Each rotation is J = [[c, s], [-conj(s), c]], as in the complex
//! one-sided Jacobi construction (LAPACK ZGESVJ/ZGSVJ0). This is a small
//! deterministic cyclic implementation, not the preconditioned LAPACK driver.
//! A global component scale and separately normalized column inner products
//! avoid squaring the input's scale or the ratio of its column norms.

use fs_math::{c64::C64, det};

/// A thin decomposition A = U diag(singular_values) V^H.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexSvd {
    /// Nonnegative singular values in descending order; no rank truncation.
    pub singular_values: Vec<f64>,
    /// Row-major m-by-n matrix; its columns are the left singular vectors.
    /// Exact zero singular values receive an orthonormal completion.
    pub u: Vec<C64>,
    /// Row-major n-by-n matrix; its columns are the right singular vectors.
    pub v: Vec<C64>,
    /// Completed cyclic sweeps (zero for initially orthogonal columns).
    pub sweeps: usize,
    /// Largest normalized Hermitian inner product between nonzero working
    /// columns at convergence, before their final normalization.
    pub max_correlation: f64,
}

/// Invalid/unsupported input, exhausted convergence, or cooperative cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum ComplexSvdError {
    /// Invalid shape/options, nonfinite data, or an unrepresentable numerical
    /// range. These cases never return a silently truncated decomposition.
    BadInput { what: &'static str },
    /// The requested relative column orthogonality was not reached.
    NotConverged { sweeps: usize, max_correlation: f64 },
    /// Cancellation was observed at a bounded row, pair, or column boundary.
    Cancelled,
}

impl core::fmt::Display for ComplexSvdError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::BadInput { what } => write!(f, "complex SVD: {what}"),
            Self::NotConverged {
                sweeps,
                max_correlation,
            } => write!(
                f,
                "complex SVD did not converge in {sweeps} sweeps (column correlation {max_correlation:e})"
            ),
            Self::Cancelled => write!(f, "complex SVD cancelled"),
        }
    }
}

impl std::error::Error for ComplexSvdError {}

fn bad(what: &'static str) -> ComplexSvdError {
    ComplexSvdError::BadInput { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ComplexSvdError> {
    if cancelled() {
        Err(ComplexSvdError::Cancelled)
    } else {
        Ok(())
    }
}

fn divide(z: C64, x: f64) -> C64 {
    // A reciprocal would overflow for representable subnormal x.
    C64::new(z.re / x, z.im / x)
}

/// Scaled norm of one strided column; no squared unscaled magnitudes.
fn column_norm(a: &[C64], n: usize, column: usize) -> f64 {
    let mut scale = 0.0_f64;
    for row in a.chunks_exact(n) {
        scale = scale.max(row[column].re.abs()).max(row[column].im.abs());
    }
    if scale == 0.0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for row in a.chunks_exact(n) {
        let z = divide(row[column], scale);
        sum += z.re.mul_add(z.re, z.im * z.im);
    }
    scale * det::sqrt(sum)
}

fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

/// Correlation is computed after normalizing EACH column, so weak independent
/// directions do not disappear when their squared norms would underflow.
fn correlation(a: &[C64], n: usize, p: usize, q: usize, ap: f64, aq: f64) -> C64 {
    if ap == 0.0 || aq == 0.0 {
        return C64::ZERO;
    }
    let (mut re, mut im, mut er, mut ei) = (0.0, 0.0, 0.0, 0.0);
    for row in a.chunks_exact(n) {
        let z = divide(row[p], ap).conj() * divide(row[q], aq);
        compensated_add(&mut re, &mut er, z.re);
        compensated_add(&mut im, &mut ei, z.im);
    }
    C64::new(re, im)
}

fn max_correlation(
    a: &[C64],
    norms: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ComplexSvdError> {
    let n = norms.len();
    let mut maximum = 0.0_f64;
    for p in 0..n {
        for q in (p + 1)..n {
            poll(cancelled)?;
            maximum = maximum.max(correlation(a, n, p, q, norms[p], norms[q]).abs());
        }
    }
    Ok(maximum)
}

fn rotate(a: &mut [C64], n: usize, p: usize, q: usize, c: f64, s: C64) {
    for row in a.chunks_exact_mut(n) {
        let (x, y) = (row[p], row[q]);
        row[p] = x.scale(c) - s.conj() * y;
        row[q] = s * x + y.scale(c);
    }
}

fn swap_columns(a: &mut [C64], n: usize, p: usize, q: usize) {
    for row in a.chunks_exact_mut(n) {
        row.swap(p, q);
    }
}

/// Fill only EXACT zero singular directions, without discarding weak modes.
/// Twice-projected coordinate vectors give a deterministic orthonormal
/// completion. The coordinate cursor only advances: at most m candidates.
fn complete_zero_columns(
    u: &mut [C64],
    m: usize,
    norms: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ComplexSvdError> {
    let n = norms.len();
    let mut coordinate = 0;
    for (j, &norm) in norms.iter().enumerate() {
        poll(cancelled)?;
        if norm != 0.0 {
            continue;
        }
        loop {
            poll(cancelled)?;
            if coordinate == m {
                return Err(bad("orthonormal completion lost numerical rank"));
            }
            for (i, row) in u.chunks_exact_mut(n).enumerate() {
                row[j] = if i == coordinate { C64::ONE } else { C64::ZERO };
            }
            coordinate += 1;
            for _ in 0..2 {
                for k in 0..j {
                    poll(cancelled)?;
                    let mut dot = C64::ZERO;
                    for row in u.chunks_exact(n) {
                        dot = dot + row[k].conj() * row[j];
                    }
                    for row in u.chunks_exact_mut(n) {
                        row[j] = row[j] - row[k] * dot;
                    }
                }
            }
            let length = column_norm(u, n, j);
            if length > 32.0 * m as f64 * f64::EPSILON {
                for row in u.chunks_exact_mut(n) {
                    row[j] = divide(row[j], length);
                }
                break;
            }
        }
    }
    Ok(())
}

/// Compute a thin complex SVD of row-major a, with 0 <= n <= m <= 1024.
///
/// The tolerance bounds |u_p^H u_q| for nonzero working columns and must be
/// finite, greater than machine epsilon, and less than one. The sweep cap must
/// be positive. A full convergence check follows each cyclic sweep; exhaustion
/// is an error, never a partial decomposition. The caller controls numerical
/// rank separately: this routine does not zero small singular values.
///
/// Input components are divided by a single maximum component magnitude.
/// Representability failures (including a nonzero component lost by scaling,
/// an underflowed required Jacobi angle, or an unrepresentable singular value)
/// return BadInput. This deliberately makes no full-dynamic-range claim for
/// matrices spanning more than the representable scaled range.
///
/// Arithmetic is deterministic: fixed pair order, descending singular-value
/// order with original-column tie breaks, and a fixed column phase. Storage
/// beyond borrowed a is two dense arrays of m*n + n*n complex entries plus
/// O(n) scalars. Work is O(max_sweeps * m * n^2), plus O(m^2*n) only when exact
/// zero columns require completion. Cancellation is checked at every pair and
/// throughout input/output preparation, with at most O(m+n) work between polls.
///
/// # Errors
/// Returns ComplexSvdError for invalid or unsupported inputs, failure to
/// converge within the explicit sweep cap, or cancellation.
pub fn svd_complex_with_cancel(
    a: &[C64],
    m: usize,
    n: usize,
    tolerance: f64,
    max_sweeps: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ComplexSvd, ComplexSvdError> {
    poll(&mut cancelled)?;
    if n > m || m > 1024 || n > 1024 || m.checked_mul(n) != Some(a.len()) {
        return Err(bad("expected row-major m*n input with 0 <= n <= m <= 1024"));
    }
    if !tolerance.is_finite() || tolerance <= f64::EPSILON || tolerance >= 1.0 {
        return Err(bad(
            "tolerance must be finite and between machine epsilon and one",
        ));
    }
    if max_sweeps == 0 {
        return Err(bad("max_sweeps must be positive"));
    }
    if n == 0 {
        return Ok(ComplexSvd {
            singular_values: Vec::new(),
            u: Vec::new(),
            v: Vec::new(),
            sweeps: 0,
            max_correlation: 0.0,
        });
    }
    let mut scale = 0.0_f64;
    for row in a.chunks_exact(n) {
        poll(&mut cancelled)?;
        for z in row {
            if !z.re.is_finite() || !z.im.is_finite() {
                return Err(bad("input contains a nonfinite component"));
            }
            scale = scale.max(z.re.abs()).max(z.im.abs());
        }
    }
    let mut b = Vec::with_capacity(a.len());
    for row in a.chunks_exact(n) {
        poll(&mut cancelled)?;
        for &z in row {
            let w = if scale == 0.0 { z } else { divide(z, scale) };
            if (z.re != 0.0 && w.re == 0.0) || (z.im != 0.0 && w.im == 0.0) {
                return Err(bad("global scaling would erase a nonzero input component"));
            }
            b.push(w);
        }
    }
    let mut v = vec![C64::ZERO; n * n];
    for j in 0..n {
        poll(&mut cancelled)?;
        v[j * n + j] = C64::ONE;
    }
    if scale == 0.0 {
        // The all-zero case has a canonical answer; avoid cubic completion.
        for j in 0..n {
            poll(&mut cancelled)?;
            b[j * n + j] = C64::ONE;
        }
        return Ok(ComplexSvd {
            singular_values: vec![0.0; n],
            u: b,
            v,
            sweeps: 0,
            max_correlation: 0.0,
        });
    }
    let mut norms = Vec::with_capacity(n);
    for j in 0..n {
        poll(&mut cancelled)?;
        norms.push(column_norm(&b, n, j));
    }
    let mut maximum = max_correlation(&b, &norms, &mut cancelled)?;
    let mut sweeps = 0;
    while maximum > tolerance && sweeps < max_sweeps {
        poll(&mut cancelled)?;
        for p in 0..n {
            for q in (p + 1)..n {
                poll(&mut cancelled)?;
                let (ap, aq) = (norms[p], norms[q]);
                let g = correlation(&b, n, p, q, ap, aq);
                let gamma = g.abs();
                if gamma <= tolerance {
                    continue;
                }
                // Stable version of tan(theta) for the 2x2 Hermitian Gram
                // block: delta = beta_q - beta_p, off = |p^H q|. Only
                // norm ratios are squared, never the input's absolute scale.
                let h = ap.max(aq);
                let (x, y) = (ap / h, aq / h);
                let delta = (y - x) * (y + x);
                let off = x * y * gamma;
                let denominator = delta.abs() + det::hypot(delta, 2.0 * off);
                if off == 0.0 || denominator == 0.0 {
                    return Err(bad(
                        "required Jacobi angle is outside the representable range",
                    ));
                }
                let sign = if delta < 0.0 { -1.0 } else { 1.0 };
                let t = sign * (2.0 * off / denominator);
                let c = 1.0 / det::sqrt(1.0 + t * t);
                let s = divide(g, gamma).scale(t * c);
                rotate(&mut b, n, p, q, c, s);
                rotate(&mut v, n, p, q, c, s);
                norms[p] = column_norm(&b, n, p);
                norms[q] = column_norm(&b, n, q);
            }
        }
        sweeps += 1;
        maximum = max_correlation(&b, &norms, &mut cancelled)?;
    }
    if maximum > tolerance {
        return Err(ComplexSvdError::NotConverged {
            sweeps,
            max_correlation: maximum,
        });
    }

    // Deterministic selection ordering keeps both dense arrays in place.
    let mut labels: Vec<usize> = (0..n).collect();
    for j in 0..n {
        poll(&mut cancelled)?;
        let mut best = j;
        for k in (j + 1)..n {
            if norms[k] > norms[best] || (norms[k] == norms[best] && labels[k] < labels[best]) {
                best = k;
            }
        }
        norms.swap(j, best);
        labels.swap(j, best);
        swap_columns(&mut b, n, j, best);
        swap_columns(&mut v, n, j, best);
    }
    let mut singular_values = Vec::with_capacity(n);
    for (j, &length) in norms.iter().enumerate() {
        poll(&mut cancelled)?;
        let sigma = length * scale;
        if !sigma.is_finite() || (length != 0.0 && sigma == 0.0) {
            return Err(bad("singular value is outside the representable range"));
        }
        singular_values.push(sigma);
        if length != 0.0 {
            for row in b.chunks_exact_mut(n) {
                row[j] = divide(row[j], length);
            }
        }
    }
    complete_zero_columns(&mut b, m, &norms, &mut cancelled)?;

    // Fix the phase using the largest left-vector coordinate; ties use the
    // first row. Apply the same phase to U and V to preserve U*S*V^H.
    for j in 0..n {
        poll(&mut cancelled)?;
        let mut pivot = b[j];
        for row in b.chunks_exact(n).skip(1) {
            if row[j].abs() > pivot.abs() {
                pivot = row[j];
            }
        }
        let phase = divide(pivot.conj(), pivot.abs());
        for row in b.chunks_exact_mut(n) {
            row[j] = row[j] * phase;
        }
        for row in v.chunks_exact_mut(n) {
            row[j] = row[j] * phase;
        }
    }
    poll(&mut cancelled)?;
    Ok(ComplexSvd {
        singular_values,
        u: b,
        v,
        sweeps,
        max_correlation: maximum,
    })
}

#[cfg(test)]
mod tests;

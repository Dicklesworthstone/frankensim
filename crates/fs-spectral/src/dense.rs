//! Complete nonsymmetric dense numerical spectra through existing L1 kernels.
//!
//! fs-la owns Hessenberg/QR and complex Jacobi SVD. This service budgets their
//! composition, retains EVERY QR root, obtains right vectors without squaring
//! a nearly singular matrix, and checks original-row residuals and the complete
//! eigenvector basis. A numerically unresolved basis refuses as a whole.
//! There is no interval enclosure, spectral separation proof or forward-error
//! guarantee for a nonnormal matrix. Callers supply dimensionless coordinates.

use fs_la::eigen_complex::{ComplexSvdError, EigError, eig_with_cancel, svd_complex_with_cancel};
use fs_math::c64::C64;

/// Explicit native iteration caps, numerical acceptance and aggregate resources.
#[derive(Clone, Copy, Debug)]
pub struct DenseSpectrumOptions {
    pub maximum_dimension: usize,
    /// Includes all returned eigenvectors and native scratch, excluding the
    /// borrowed matrix and allocator metadata. No partial spectrum is returned.
    pub maximum_dense_bytes: usize,
    /// Conservative native QR + at-most-n nullspace SVDs + one basis SVD work.
    /// O(qr_sweeps*n^4 + svd_sweeps*n^4); this is a cold dense analysis.
    pub maximum_work: u64,
    pub qr_max_sweeps_per_window: usize,
    pub svd_max_sweeps: usize,
    pub svd_tolerance: f64,
    /// Maximum original-row backward residual of A*v=lambda*v, in (0,1).
    pub maximum_relative_residual: f64,
    /// Minimum sigma_min(V)/sigma_max(V) of the COMPLETE unit-vector basis,
    /// in (0,1). Coordinate dependent; a numerical resolution gate, not a proof.
    pub minimum_eigenbasis_reciprocal_condition: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum DenseSpectrumError {
    Invalid { what: &'static str },
    Eigen(EigError),
    SingularVectors(ComplexSvdError),
    EigenpairResidual { index: usize, relative: f64 },
    UnresolvedEigenbasis { reciprocal_condition: f64 },
    Cancelled,
}
impl core::fmt::Display for DenseSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dense spectrum refused: {self:?}")
    }
}
impl std::error::Error for DenseSpectrumError {}

/// Complete dimension-many numerical eigenpairs in native canonical (re,im)
/// order. Repeated values retain algebraic multiplicity; no pole is merged.
#[derive(Debug, Clone)]
pub struct DenseSpectrum {
    pub values: Vec<C64>,
    /// One unit Euclidean right vector per value, in original supplied order.
    pub right_vectors: Vec<Vec<C64>>,
    pub relative_residuals: Vec<f64>,
    pub eigenbasis_reciprocal_condition: f64,
}

fn invalid(what: &'static str) -> DenseSpectrumError {
    DenseSpectrumError::Invalid { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DenseSpectrumError> {
    if cancelled() {
        Err(DenseSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Native QR discovery, native SVD right vectors and complete-basis admission.
/// No new eigenvalue, inverse-iteration or pole-refinement algorithm is used.
/// Exactly coincident values use distinct vectors from the SAME nullspace SVD;
/// unresolved near coincidences or defective roots are never silently pruned.
pub fn nonsymmetric_spectrum(
    matrix: &[C64],
    n: usize,
    options: DenseSpectrumOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DenseSpectrum, DenseSpectrumError> {
    poll(&mut cancelled)?;
    if n == 0
        || n > 1024
        || n > options.maximum_dimension
        || n.checked_mul(n) != Some(matrix.len())
        || matrix.iter().any(|&z| !finite(z))
        || options.qr_max_sweeps_per_window == 0
        || options.svd_max_sweeps == 0
        || !options.svd_tolerance.is_finite()
        || options.svd_tolerance <= f64::EPSILON
        || options.svd_tolerance >= 1.
        || !options.maximum_relative_residual.is_finite()
        || options.maximum_relative_residual <= 0.
        || options.maximum_relative_residual >= 1.
        || !options.minimum_eigenbasis_reciprocal_condition.is_finite()
        || options.minimum_eigenbasis_reciprocal_condition <= 0.
        || options.minimum_eigenbasis_reciprocal_condition >= 1.
    {
        return Err(invalid(
            "complete finite square matrix and valid numerical controls required",
        ));
    }
    let size = n as u128;
    let bytes = 16 * (18 * size * size + 24 * size);
    let work = 32 * options.qr_max_sweeps_per_window as u128 * size.pow(4)
        + 64 * options.svd_max_sweeps as u128 * (size + 1) * size.pow(3)
        + 256 * size.pow(3);
    if bytes > options.maximum_dense_bytes as u128 || work > options.maximum_work as u128 {
        return Err(invalid(
            "complete spectrum exceeds its declared dense payload or work budget",
        ));
    }
    let values = eig_with_cancel(matrix, n, options.qr_max_sweeps_per_window, &mut cancelled)
        .map_err(|e| {
            if e == EigError::Cancelled {
                DenseSpectrumError::Cancelled
            } else {
                DenseSpectrumError::Eigen(e)
            }
        })?;
    if values.len() != n || values.iter().any(|&z| !finite(z)) {
        return Err(invalid(
            "native QR did not return the complete finite spectrum",
        ));
    }
    let mut vectors = Vec::with_capacity(n);
    let mut residuals = Vec::with_capacity(n);
    let mut at = 0;
    while at < n {
        poll(&mut cancelled)?;
        let value = values[at];
        let mut end = at + 1;
        while end < n && values[end] == value {
            end += 1;
        }
        let mut shifted = matrix.to_vec();
        for i in 0..n {
            shifted[i * n + i] = shifted[i * n + i] - value;
        }
        let svd = svd_complex_with_cancel(
            &shifted,
            n,
            n,
            options.svd_tolerance,
            options.svd_max_sweeps,
            &mut cancelled,
        )
        .map_err(map_svd)?;
        for index in at..end {
            let column = n - 1 - (index - at);
            let vector: Vec<_> = (0..n).map(|i| svd.v[i * n + column]).collect();
            let rows = eigenpair_row_residuals(matrix, n, value, &vector, &mut cancelled)?;
            let relative = rows.into_iter().fold(0., f64::max);
            if relative > options.maximum_relative_residual {
                return Err(DenseSpectrumError::EigenpairResidual { index, relative });
            }
            vectors.push(vector);
            residuals.push(relative);
        }
        at = end;
    }
    // Individually tiny residuals do not establish distinct independent modes.
    // A near-defective or ambiguously recovered basis must not masquerade as
    // a complete diagonalization merely because every returned vector is close.
    let mut basis = vec![C64::ZERO; n * n];
    for i in 0..n {
        poll(&mut cancelled)?;
        for j in 0..n {
            basis[i * n + j] = vectors[j][i];
        }
    }
    let svd = svd_complex_with_cancel(
        &basis,
        n,
        n,
        options.svd_tolerance,
        options.svd_max_sweeps,
        &mut cancelled,
    )
    .map_err(map_svd)?;
    let reciprocal_condition = svd.singular_values[n - 1] / svd.singular_values[0];
    if !reciprocal_condition.is_finite()
        || reciprocal_condition < options.minimum_eigenbasis_reciprocal_condition
    {
        return Err(DenseSpectrumError::UnresolvedEigenbasis {
            reciprocal_condition,
        });
    }
    poll(&mut cancelled)?;
    Ok(DenseSpectrum {
        values,
        right_vectors: vectors,
        relative_residuals: residuals,
        eigenbasis_reciprocal_condition: reciprocal_condition,
    })
}

fn map_svd(error: ComplexSvdError) -> DenseSpectrumError {
    if error == ComplexSvdError::Cancelled {
        DenseSpectrumError::Cancelled
    } else {
        DenseSpectrumError::SingularVectors(error)
    }
}

/// Actual original-row backward residuals, also reusable after a domain's
/// coordinate/time back-transform. Row scaling avoids overflowing coefficients;
/// it does not replace the supplied equations or use a near-zero T(lambda) norm.
/// Each denominator is (sum_j |A_ij| + |lambda|)*||v||_infinity.
/// This remains meaningful for zero eigenvector components and dark coordinates;
/// it is a row-norm estimate, not a componentwise relative perturbation bound.
pub fn eigenpair_row_residuals(
    matrix: &[C64],
    n: usize,
    value: C64,
    vector: &[C64],
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<f64>, DenseSpectrumError> {
    poll(&mut cancelled)?;
    if n == 0
        || n.checked_mul(n) != Some(matrix.len())
        || vector.len() != n
        || !finite(value)
        || matrix.iter().chain(vector).any(|&z| !finite(z))
        || vector.iter().all(|&z| z == C64::ZERO)
    {
        return Err(invalid(
            "finite complete nonzero eigenvector and original matrix required",
        ));
    }
    let vector_norm = vector.iter().map(|z| z.abs()).fold(0., f64::max);
    if !vector_norm.is_finite() {
        return Err(invalid("eigenvector norm is unrepresentable"));
    }
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        poll(&mut cancelled)?;
        let row = &matrix[i * n..(i + 1) * n];
        let scale = row.iter().map(|z| z.abs()).fold(value.abs(), f64::max);
        if scale == 0. {
            result.push(0.);
            continue;
        }
        if !scale.is_finite() {
            return Err(invalid("original residual row scale is unrepresentable"));
        }
        let lambda = C64::new(value.re / scale, value.im / scale);
        let mut sum = -lambda * vector[i];
        let mut correction = C64::ZERO;
        let mut denominator = lambda.abs();
        for (&a, &x) in row.iter().zip(vector) {
            let a = C64::new(a.re / scale, a.im / scale);
            let term = a * x - correction;
            let next = sum + term;
            correction = (next - sum) - term;
            sum = next;
            denominator += a.abs();
        }
        denominator *= vector_norm;
        let relative = if denominator == 0. {
            if sum == C64::ZERO { 0. } else { f64::INFINITY }
        } else {
            sum.abs() / denominator
        };
        if !finite(sum) || !denominator.is_finite() || !relative.is_finite() {
            return Err(invalid("original eigenpair residual is unrepresentable"));
        }
        result.push(relative);
    }
    poll(&mut cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;

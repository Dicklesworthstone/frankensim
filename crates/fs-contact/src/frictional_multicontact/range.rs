//! A sufficient velocity-determinacy condition for redundant Coulomb contacts.
//!
//! For two law-satisfying impulses, normal complementarity gives
//! `delta p_n dot delta v_n <= 0`. Monotonicity of the Euclidean norm's
//! subgradient gives, contact by contact,
//! `delta p_t dot delta v_t <= mu * |delta p_n| * |delta v_t|`.
//! Consequently E = delta p^T W delta p satisfies
//! `E <= |D N delta p| * |T W delta p|`, where D = diag(mu), and N/T select
//! normal/tangential rows. If D N annihilates ker(W), the two norms are bounded
//! by `sqrt(||D N W^+ N^T D||) * sqrt(E)` and
//! `sqrt(||T W T^T||) * sqrt(E)`. A product below one implies E = 0, hence
//! unique W*p, but does NOT imply unique p. Symmetric absolute row sums bound
//! the two spectral norms. This is sufficient, not an existence classifier.
//!
//! Generic eigendecomposition remains in fs-la. Its range is used only for
//! this numerical witness: Newton, contact, and energy checks always use the
//! supplied original W. No eigenvalue shift, compliance, or truncated matrix
//! enters the physical solve. The small matrix cap bounds each unpolled call.

use super::{FrictionalContactError, PsdVelocityWitness, poll};
use fs_solver::orthant::active_set::SpdOrthantConfig;

pub(super) fn velocity_witness(
    original: &[f64],
    mu: &[f64],
    config: SpdOrthantConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<PsdVelocityWitness>, FrictionalContactError> {
    let n = mu.len();
    let d = 3 * n;
    poll(cancelled)?;
    if d == 0 || d > config.max_variables || d > 63 {
        return Ok(None);
    }
    // These checks are repeated because a failed Cholesky can precede the
    // original owner's full symmetry validation (e.g. a nonpositive diagonal).
    let scale = original
        .chunks_exact(d)
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return Ok(None);
    }
    let mut w = vec![0.0; d * d];
    for i in 0..d {
        poll(cancelled)?;
        if original[i * d + i] <= 0.0 {
            return Ok(None);
        }
        for j in 0..=i {
            let a = original[i * d + j];
            let b = original[j * d + i];
            let symmetry = (a / original[i * d + i].sqrt() / original[j * d + j].sqrt()
                - b / original[i * d + i].sqrt() / original[j * d + j].sqrt())
            .abs();
            if !symmetry.is_finite() || symmetry > config.symmetry_tolerance {
                return Ok(None);
            }
            let value = 0.5 * (a / scale) + 0.5 * (b / scale);
            w[i * d + j] = value;
            w[j * d + i] = value;
        }
    }
    poll(cancelled)?;
    let (eigenvalues, q) = fs_la::eigen::jacobi_eigh(&w, d);
    poll(cancelled)?;
    if eigenvalues.iter().chain(&q).any(|v| !v.is_finite()) {
        return Ok(None);
    }
    let zero_band = 64.0 * f64::EPSILON * d as f64;
    // An eigenvalue between numerical zero and a separated positive range is
    // rank-ambiguous. Refuse it rather than erasing a weak physical direction.
    let positive_band = (64.0 * zero_band).max(config.relative_pivot_tolerance);
    let mut spectral_residual = 0.0_f64;
    let mut orthogonality_residual = 0.0_f64;
    for i in 0..d {
        poll(cancelled)?;
        for j in 0..d {
            let mut image = 0.0;
            let mut gram = 0.0;
            for k in 0..d {
                // Check the original, unaveraged matrix, too. This does not
                // accept a spectral witness for a different physical matrix.
                image += (original[i * d + k] / scale) * q[k * d + j];
                gram += q[k * d + i] * q[k * d + j];
            }
            spectral_residual =
                spectral_residual.max((image - q[i * d + j] * eigenvalues[j]).abs());
            orthogonality_residual =
                orthogonality_residual.max((gram - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    if spectral_residual > zero_band || orthogonality_residual > zero_band {
        return Ok(None);
    }
    let mut retained = Vec::new();
    let mut null = Vec::new();
    for (i, &value) in eigenvalues.iter().enumerate() {
        if value.abs() <= zero_band {
            null.push(i);
        } else if value >= positive_band {
            retained.push(i);
        } else {
            // Includes a negative eigenvalue outside the roundoff band.
            return Ok(None);
        }
    }
    if retained.is_empty() || null.is_empty() {
        return Ok(None);
    }
    let mut null_residual = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        if mu[i] == 0.0 {
            continue;
        }
        let row_sum = null.iter().map(|&k| q[3 * i * d + k].abs()).sum::<f64>();
        // Any nonzero friction coefficient requires the corresponding normal
        // impulse to be range-observable. A tiny mu is not permission to treat
        // an unbounded normal null direction as a finite perturbation.
        if row_sum > zero_band {
            return Ok(None);
        }
        null_residual = null_residual.max(mu[i] * row_sum);
    }
    let mut normal_inverse_bound = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        let mut row_sum = 0.0;
        for j in 0..n {
            let mut value = 0.0;
            for &k in &retained {
                value += q[3 * i * d + k] * q[3 * j * d + k] / eigenvalues[k];
            }
            row_sum += (mu[i] * value * mu[j]).abs();
        }
        normal_inverse_bound = normal_inverse_bound.max(row_sum);
    }
    let mut tangent_bound = 0.0_f64;
    for i in (0..d).filter(|i| i % 3 != 0) {
        poll(cancelled)?;
        let row_sum = (0..d)
            .filter(|j| j % 3 != 0)
            .map(|j| w[i * d + j].abs())
            .sum::<f64>();
        tangent_bound = tangent_bound.max(row_sum);
    }
    // The common normalization of W cancels between the inverse and mobility.
    let bound =
        (normal_inverse_bound * tangent_bound).sqrt() * (1.0 + 512.0 * f64::EPSILON * d as f64);
    if !bound.is_finite() || bound >= 1.0 - 1e-8 || !null_residual.is_finite() {
        return Ok(None);
    }
    poll(cancelled)?;
    Ok(Some(PsdVelocityWitness {
        dimension: d,
        rank: retained.len(),
        friction_feedback_bound: bound,
        weighted_normal_nullspace_residual: null_residual,
        spectral_residual,
        orthogonality_residual,
        minimum_retained_eigenvalue: eigenvalues[retained[0]],
        maximum_null_eigenvalue_magnitude: null
            .iter()
            .map(|&k| eigenvalues[k].abs())
            .fold(0.0, f64::max),
    }))
}

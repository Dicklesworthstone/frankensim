//! Allocation-free, scaled Jacobi eigenvalues for one symmetric 3x3 matrix.
//!
//! Uses the same cyclic plane rotations as the general Jacobi owner, with
//! fixed storage and no eigenvectors. Scaling and a bounded rotation formula
//! avoid squaring dimensional entries or forming an unbounded Jacobi theta.
//! Unlike cubic roots, the rotations do not lose sqrt(epsilon) accuracy near
//! a repeated eigenvalue. These are floating-point estimates, not intervals.

/// Refusal from the bounded symmetric 3x3 eigenvalue kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymmetricEigenvalues3Error {
    /// Every matrix entry must be finite.
    NonFinite,
    /// Both triangles must describe the same symmetric matrix.
    NotSymmetric,
    /// A resulting eigenvalue is outside the finite f64 range.
    Unrepresentable,
    /// The fixed sweep budget did not reduce the off-diagonal residual.
    NoConvergence,
}

/// Eigenvalues in ascending order, without allocation or a dynamic work size.
///
/// At most 16 cyclic Jacobi sweeps are performed. Diagonal matrices retain
/// their entries exactly. Otherwise entries are scaled by their maximum
/// magnitude, and convergence requires off-diagonals <= 4*f64::EPSILON in
/// that scale. Deterministic pair order and strict fs-math square roots apply.
/// Relative accuracy of an eigenvalue close to zero is not guaranteed.
pub fn symmetric_eigenvalues3(
    matrix: &[[f64; 3]; 3],
) -> Result<[f64; 3], SymmetricEigenvalues3Error> {
    if matrix.iter().flatten().any(|x| !x.is_finite()) {
        return Err(SymmetricEigenvalues3Error::NonFinite);
    }
    if matrix[0][1] != matrix[1][0]
        || matrix[0][2] != matrix[2][0]
        || matrix[1][2] != matrix[2][1]
    {
        return Err(SymmetricEigenvalues3Error::NotSymmetric);
    }
    if matrix[0][1] == 0.0 && matrix[0][2] == 0.0 && matrix[1][2] == 0.0 {
        let mut values = [matrix[0][0], matrix[1][1], matrix[2][2]];
        values.sort_by(f64::total_cmp);
        return Ok(values);
    }
    let scale = matrix.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
    let mut a = matrix.map(|row| row.map(|x| x / scale));
    let tolerance = 4.0 * f64::EPSILON;
    for _ in 0..16 {
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            let b = a[p][q];
            if b.abs() <= tolerance {
                continue;
            }
            let delta = 0.5 * (a[q][q] - a[p][p]);
            let h = delta.abs().max(b.abs());
            let radius = h * fs_math::det::sqrt((delta / h).powi(2) + (b / h).powi(2));
            let t = b / (delta + radius.copysign(delta));
            let c = 1.0 / fs_math::det::sqrt(1.0 + t * t);
            let s = t * c;
            a[p][p] -= t * b;
            a[q][q] += t * b;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            // The one unaffected coordinate completes the symmetric rotation.
            let k = 3 - p - q;
            let (kp, kq) = (a[k][p], a[k][q]);
            a[k][p] = c * kp - s * kq;
            a[p][k] = a[k][p];
            a[k][q] = s * kp + c * kq;
            a[q][k] = a[k][q];
        }
        if a[0][1].abs().max(a[0][2].abs()).max(a[1][2].abs()) <= tolerance {
            let mut values = [a[0][0] * scale, a[1][1] * scale, a[2][2] * scale];
            if values.iter().any(|x| !x.is_finite()) {
                return Err(SymmetricEigenvalues3Error::Unrepresentable);
            }
            values.sort_by(f64::total_cmp);
            return Ok(values);
        }
    }
    Err(SymmetricEigenvalues3Error::NoConvergence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_retains_extreme_entries_exactly() {
        assert_eq!(
            symmetric_eigenvalues3(&[
                [f64::MAX, 0.0, 0.0],
                [0.0, f64::MIN_POSITIVE, 0.0],
                [0.0, 0.0, -f64::MAX],
            ]).unwrap(),
            [-f64::MAX, f64::MIN_POSITIVE, f64::MAX]
        );
    }

    #[test]
    fn shear_and_repeated_roots_remain_accurate_across_scales() {
        for scale in [1.0e-300, 1.0, 1.0e300] {
            let values = symmetric_eigenvalues3(&[
                [0.0, scale, 0.0],
                [scale, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]).unwrap();
            assert_eq!(values, [-scale, 0.0, scale]);
            // I - 3 n n^T, n=(1,1,1)/sqrt(3): eigenvalues -2, 1, 1.
            let values = symmetric_eigenvalues3(&[
                [0.0, -scale, -scale],
                [-scale, 0.0, -scale],
                [-scale, -scale, 0.0],
            ]).unwrap();
            for (actual, expected) in values.into_iter().zip([-2.0, 1.0, 1.0]) {
                assert!((actual / scale - expected).abs() < 3.0e-15);
            }
        }
    }

    #[test]
    fn finite_and_symmetric_admission_is_explicit() {
        assert_eq!(
            symmetric_eigenvalues3(&[[f64::NAN; 3]; 3]),
            Err(SymmetricEigenvalues3Error::NonFinite)
        );
        assert_eq!(
            symmetric_eigenvalues3(&[
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]),
            Err(SymmetricEigenvalues3Error::NotSymmetric)
        );
        assert_eq!(
            symmetric_eigenvalues3(&[[f64::MAX; 3]; 3]),
            Err(SymmetricEigenvalues3Error::Unrepresentable)
        );
    }
}

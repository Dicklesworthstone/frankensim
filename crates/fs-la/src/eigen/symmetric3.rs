//! Allocation-free Jacobi eigenvalues for one symmetric 3x3 matrix.
//!
//! Plane rotations use local scaling and fixed storage. Dimensional entries
//! remain intact, so a large compression cannot erase a small tensile block.
//! No dimensional squares or unbounded Jacobi theta are formed. Unlike cubic
//! roots, rotations avoid sqrt(epsilon) error near repeated eigenvalues.
//! These are floating-point estimates, not interval enclosures.

/// Refusal from the bounded symmetric 3x3 eigenvalue kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymmetricEigenvalues3Error {
    /// Every matrix entry must be finite.
    NonFinite,
    /// Both triangles must describe the same symmetric matrix.
    NotSymmetric,
    /// A required rotation or resulting eigenvalue is not representable.
    Unrepresentable,
    /// The fixed sweep budget did not reduce the off-diagonal residual.
    NoConvergence,
}

/// Eigenvalues in ascending order, without allocation or a dynamic work size.
///
/// At most 16 cyclic Jacobi sweeps are performed. Diagonal matrices retain
/// their entries exactly. Each rotation is scaled locally; convergence uses
/// |a_pq| <= 4*EPSILON*sqrt(|a_pp|)*sqrt(|a_qq|), never a global matrix norm.
/// A zero diagonal therefore cannot hide a nonzero shear. Deterministic pair
/// order and strict fs-math square roots apply. Relative accuracy of a nearly
/// zero eigenvalue produced by cancellation is not guaranteed.
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
    let mut a = *matrix;
    for _ in 0..16 {
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            let b = a[p][q];
            if negligible_pair(a[p][p], a[q][q], b) {
                continue;
            }
            let scale = a[p][p].abs().max(a[q][q].abs()).max(b.abs());
            let delta = 0.5 * (a[q][q] / scale - a[p][p] / scale);
            let scaled_b = b / scale;
            if scaled_b == 0.0 {
                return Err(SymmetricEigenvalues3Error::Unrepresentable);
            }
            let h = delta.abs().max(scaled_b.abs());
            let radius = h
                * fs_math::det::sqrt((delta / h).powi(2) + (scaled_b / h).powi(2));
            let t = scaled_b / (delta + radius.copysign(delta));
            // Apply the correction in physical units: normalizing the entire
            // matrix would underflow b^2/|a_pp| even when it is representable.
            let shift = t * b;
            if !t.is_finite()
                || t == 0.0
                || !shift.is_finite()
                || (shift == 0.0 && (a[p][p] == 0.0 || a[q][q] == 0.0))
            {
                // Do not silently turn a nascent positive eigenvalue into zero.
                // A subnormal correction to two nonzero diagonals can round to
                // zero without changing either sign; the rotation still proceeds.
                return Err(SymmetricEigenvalues3Error::Unrepresentable);
            }
            let c = 1.0 / fs_math::det::sqrt(1.0 + t * t);
            let s = t * c;
            a[p][p] -= shift;
            a[q][q] += shift;
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            // The one unaffected coordinate completes the symmetric rotation.
            let k = 3 - p - q;
            let (kp, kq) = (a[k][p], a[k][q]);
            a[k][p] = c * kp - s * kq;
            a[p][k] = a[k][p];
            a[k][q] = s * kp + c * kq;
            a[q][k] = a[k][q];
            if [a[p][p], a[q][q], a[k][p], a[k][q]]
                .iter()
                .any(|x| !x.is_finite())
            {
                return Err(SymmetricEigenvalues3Error::Unrepresentable);
            }
        }
        if negligible_pair(a[0][0], a[1][1], a[0][1])
            && negligible_pair(a[0][0], a[2][2], a[0][2])
            && negligible_pair(a[1][1], a[2][2], a[1][2])
        {
            let mut values = [a[0][0], a[1][1], a[2][2]];
            values.sort_by(f64::total_cmp);
            return Ok(values);
        }
    }
    Err(SymmetricEigenvalues3Error::NoConvergence)
}

fn negligible_pair(pp: f64, qq: f64, pq: f64) -> bool {
    pq == 0.0
        || pq.abs()
            <= (4.0 * f64::EPSILON * fs_math::det::sqrt(pp.abs()))
                * fs_math::det::sqrt(qq.abs())
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
    #[test]
    fn small_tensile_values_survive_large_compression() {
        for (compression, shear) in [(1.0e16, 1.0), (1.0e300, 1.0e-300)] {
            for [isolated, p, q] in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
                let mut matrix = [[0.0; 3]; 3];
                matrix[isolated][isolated] = -compression;
                matrix[p][q] = shear;
                matrix[q][p] = shear;
                assert_eq!(
                    symmetric_eigenvalues3(&matrix).unwrap(),
                    [-compression, -shear, shear]
                );
            }
        }
        // The positive root is b^2/|a| to far beyond f64 precision here.
        // Its normalized value underflows, but its dimensional value does not.
        let values = symmetric_eigenvalues3(&[
            [-1.0e300, 1.0e50, 0.0],
            [1.0e50, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]).unwrap();
        assert_eq!(values[0], -1.0e300);
        assert_eq!(values[1], 0.0);
        assert!((values[2] / 1.0e-200 - 1.0).abs() < 2.0e-15);
        assert_eq!(
            symmetric_eigenvalues3(&[
                [-f64::MAX, f64::MIN_POSITIVE, 0.0],
                [f64::MIN_POSITIVE, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]),
            Err(SymmetricEigenvalues3Error::Unrepresentable)
        );
    }

}

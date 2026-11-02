//! Minimum-change angular correction in the inertia's represented range.

use fs_la::factor::svd_jacobi;

/// The mass-weighted least-change velocity correction is `omega × r` about
/// the receiving centroid. Collinear cells have rank-two inertia: they carry
/// transverse angular momentum but cannot carry spin about their common line.
/// Scale before the native SVD so rank admission is independent of units.
pub(super) fn angular_correction(
    inertia: &[f64; 9],
    angular: [f64; 3],
    tolerance: f64,
) -> Option<[f64; 3]> {
    if !tolerance.is_finite()
        || tolerance < 0.0
        || inertia.iter().chain(&angular).any(|x| !x.is_finite())
    {
        return None;
    }
    let norm = |v: [f64; 3]| v[0].hypot(v[1]).hypot(v[2]);
    if norm(angular) <= tolerance {
        return Some([0.0; 3]);
    }
    let scale = inertia.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if scale == 0.0 {
        return None;
    }
    let a = inertia.map(|x| x / scale);
    let b = angular.map(|x| x / scale);
    if b.iter().any(|x| !x.is_finite()) {
        return None;
    }
    let svd = svd_jacobi(&a, 3, 3);
    let rank_floor = 64.0 * f64::EPSILON * svd.sigma[0];
    let mut omega = [0.0; 3];
    for mode in 0..3 {
        if svd.sigma[mode] > rank_floor {
            let projection = (0..3)
                .fold(0.0, |sum, row| svd.u[3 * row + mode].mul_add(b[row], sum))
                / svd.sigma[mode];
            for (row, value) in omega.iter_mut().enumerate() {
                *value = svd.v[3 * row + mode].mul_add(projection, *value);
            }
        }
    }
    if omega.iter().any(|x| !x.is_finite()) {
        return None;
    }
    // Reapply the original matrix, independently of the SVD decomposition.
    // A discarded direction is admissible only within the caller's existing
    // physical angular-momentum allowance, never a rank-dependent tolerance.
    let residual = core::array::from_fn(|row| {
        (0..3).fold(-b[row], |sum, col| {
            a[3 * row + col].mul_add(omega[col], sum)
        })
    });
    (norm(residual) <= tolerance / scale).then_some(omega)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
        for axis in 0..3 {
            assert!(
                (actual[axis] - expected[axis]).abs() <= tolerance,
                "actual={actual:?}, expected={expected:?}"
            );
        }
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    #[test]
    fn a_cell_row_carries_both_transverse_angular_components() {
        // Two unit masses at r=(±1,0,0). Their independently prescribed
        // velocities (0,±3,∓2) have P=0 and L=(0,4,6).
        let inertia = [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
        let omega = angular_correction(&inertia, [0.0, 4.0, 6.0], 1e-13).unwrap();
        close(cross(omega, [1.0, 0.0, 0.0]), [0.0, 3.0, -2.0], 1e-14);
        close(cross(omega, [-1.0, 0.0, 0.0]), [0.0, -3.0, 2.0], 1e-14);
        assert!(angular_correction(&inertia, [1.0, 4.0, 6.0], 1e-13).is_none());
    }

    #[test]
    fn an_oblique_row_preserves_real_velocities_without_null_spin() {
        // Two unit masses at ±(1,2,3). The target velocities ±(-3,-6,5)
        // independently give L=(56,-28,0), perpendicular to the row.
        let inertia = [26.0, -4.0, -6.0, -4.0, 20.0, -12.0, -6.0, -12.0, 10.0];
        let omega = angular_correction(&inertia, [56.0, -28.0, 0.0], 1e-12).unwrap();
        close(omega, [2.0, -1.0, 0.0], 1e-14);
        close(cross(omega, [1.0, 2.0, 3.0]), [-3.0, -6.0, 5.0], 3e-14);
        assert!(angular_correction(&inertia, [57.0, -26.0, 3.0], 1e-12).is_none());
    }

    #[test]
    fn density_and_length_units_do_not_change_the_resolved_rotation() {
        for mass in [1e-120, 1.0, 1e120] {
            for length in [1e-60, 1.0, 1e60] {
                let scale = mass * length * length;
                let inertia = [0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0].map(|x| x * scale);
                let omega =
                    angular_correction(&inertia, [0.0, 4.0 * scale, 6.0 * scale], 1e-13 * scale)
                        .unwrap();
                close(omega, [0.0, 2.0, 3.0], 2e-14);
                assert!(
                    angular_correction(&inertia, [scale, 4.0 * scale, 6.0 * scale], 1e-13 * scale,)
                        .is_none()
                );
            }
        }
    }

    #[test]
    fn full_inertia_and_zero_spin_keep_their_existing_solutions() {
        let inertia = [5.0, 1.0, 0.0, 1.0, 4.0, 0.5, 0.0, 0.5, 3.0];
        let omega = angular_correction(&inertia, [0.75, -1.25, 2.75], 1e-13).unwrap();
        close(omega, [0.25, -0.5, 1.0], 1e-14);
        assert_eq!(angular_correction(&[0.0; 9], [0.0; 3], 0.0), Some([0.0; 3]));
        assert!(angular_correction(&[0.0; 9], [0.0, 1.0, 0.0], 1e-13).is_none());
    }
}

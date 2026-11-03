//! Geometry-independent least-change rotation in the represented inertia range.

use fs_la::factor::svd_jacobi;

/// Solve `I * omega = angular` for the mass-weighted correction `omega × r`.
/// Collinear cells carry transverse angular momentum, but not axial spin.
/// Scaling before the native SVD makes rank admission independent of units;
/// the original physical residual, not rank truncation, decides acceptance.
pub(crate) fn angular_correction(
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
    let residual = core::array::from_fn(|row| {
        (0..3).fold(-b[row], |sum, col| {
            a[3 * row + col].mul_add(omega[col], sum)
        })
    });
    (norm(residual) <= tolerance / scale).then_some(omega)
}

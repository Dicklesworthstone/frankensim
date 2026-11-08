//! Scalar energetic Voce terms shared by the isotropic and Hill return maps.
//!
//! Callers admit Q >= 0 and b > 0, and retain their own tensor flow/history.
//! The series are needed when b*alpha or b*delta is too small for subtraction
//! of the hardening work and stored energy to retain their convex remainder.

use fs_math::det;

pub(crate) fn strength(saturation_pa: f64, rate_per_strain: f64, alpha: f64) -> f64 {
    saturation_pa * -det::expm1(-rate_per_strain * alpha)
}

pub(crate) fn slope(saturation_pa: f64, rate_per_strain: f64, alpha: f64) -> f64 {
    saturation_pa * rate_per_strain * det::exp(-rate_per_strain * alpha)
}

pub(crate) fn secant(
    saturation_pa: f64,
    rate_per_strain: f64,
    alpha: f64,
    delta: f64,
) -> f64 {
    let x = rate_per_strain * delta;
    let ratio = if x == 0.0 { 1.0 } else { -det::expm1(-x) / x };
    saturation_pa * rate_per_strain * det::exp(-rate_per_strain * alpha) * ratio
}

/// Integral_0^alpha Q (1-exp(-b p)) dp.
pub(crate) fn storage(saturation_pa: f64, rate_per_strain: f64, alpha: f64) -> f64 {
    let x = rate_per_strain * alpha;
    if x < 0.1 {
        let ratio = -1.0 / 362880.0 + x / 3628800.0;
        let ratio = 1.0 / 40320.0 + x * ratio;
        let ratio = -1.0 / 5040.0 + x * ratio;
        let ratio = 1.0 / 720.0 + x * ratio;
        let ratio = -1.0 / 120.0 + x * ratio;
        let ratio = 1.0 / 24.0 + x * ratio;
        let ratio = -1.0 / 6.0 + x * ratio;
        let ratio = 0.5 + x * ratio;
        (saturation_pa * rate_per_strain) * alpha * (alpha * ratio)
    } else {
        saturation_pa * (alpha * (1.0 + det::expm1(-x) / x))
    }
}

// [1-exp(-x)-x exp(-x)]/x, with its continuous zero limit.
fn loss_ratio(x: f64) -> f64 {
    if x < 0.1 {
        let ratio = -1.0 / 45360.0 + x / 403200.0;
        let ratio = 1.0 / 5760.0 + x * ratio;
        let ratio = -1.0 / 840.0 + x * ratio;
        let ratio = 1.0 / 144.0 + x * ratio;
        let ratio = -1.0 / 30.0 + x * ratio;
        let ratio = 1.0 / 8.0 + x * ratio;
        let ratio = -1.0 / 3.0 + x * ratio;
        x * (0.5 + x * ratio)
    } else {
        -det::expm1(-x) / x - det::exp(-x)
    }
}

/// R(alpha+delta)*delta - Integral_alpha^(alpha+delta) R(p) dp >= 0.
pub(crate) fn numerical_loss(
    saturation_pa: f64,
    rate_per_strain: f64,
    alpha: f64,
    delta: f64,
) -> f64 {
    if delta == 0.0 {
        return 0.0;
    }
    let x = rate_per_strain * delta;
    (saturation_pa * det::exp(-rate_per_strain * alpha)) * delta * loss_ratio(x)
}

/// R(alpha+delta)-R(alpha)-delta*R'(alpha+delta), in stress units.
/// `remaining_pa` is Q*exp(-b*alpha), already prepared by the return map.
pub(crate) fn slope_remainder(remaining_pa: f64, rate_per_strain: f64, delta: f64) -> f64 {
    let x = rate_per_strain * delta;
    if x < 0.1 {
        remaining_pa * x * loss_ratio(x)
    } else if x.is_infinite() {
        remaining_pa
    } else {
        remaining_pa * (-det::expm1(-x) - x * det::exp(-x))
    }
}

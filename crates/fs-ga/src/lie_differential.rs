//! Directional differentials of the canonical SO(3) exponential Jacobians.
//!
//! These are derivatives with respect to rotation-vector coordinates, NOT a
//! new pose or rotation authority. Inverse derivatives reuse the admitted
//! canonical inverse via d(J^-1) = -J^-1 (dJ) J^-1. Analytic coefficient series
//! avoid cancellation at zero. The underlying logarithm still has a branch
//! cut at pi; differentiating a selected Jacobian does not remove that cut.

use crate::{GaError, Mat3, So3Tangent, Vec3, so3_right_jacobian_inverse};
use fs_math::det;

/// d J_r(phi)[direction], in the existing right/body convention.
///
/// # Errors
/// Refuses nonfinite coordinates or unrepresentable intermediate/output data.
pub fn so3_right_jacobian_directional(
    phi: So3Tangent, direction: So3Tangent,
) -> Result<Mat3, GaError> {
    let p = phi.angular;
    let v = direction.angular;
    finite(&[p.x, p.y, p.z, v.x, v.y, v.z])?;
    let z = p.dot(p);
    let dz = 2.0 * p.dot(v);
    finite(&[z, dz])?;
    // J_r = I - a [phi] + b [phi]^2. Differentiate the complete
    // polynomial, not just its scalar coefficients. Series use z=theta^2.
    let (a, b, az, bz) = if z < 0.01 {
        (
            0.5 + z * (-1.0/24.0 + z * (1.0/720.0 + z * (-1.0/40320.0 + z/3628800.0))),
            1.0/6.0 + z * (-1.0/120.0 + z * (1.0/5040.0 + z * (-1.0/362880.0 + z/39916800.0))),
            -1.0/24.0 + z * (1.0/360.0 + z * (-1.0/13440.0 + z/907200.0)),
            -1.0/120.0 + z * (1.0/2520.0 + z * (-1.0/120960.0 + z/9979200.0)),
        )
    } else {
        let theta = det::sqrt(z);
        let a = (1.0 - det::cos(theta)) / z;
        let b = (theta - det::sin(theta)) / (z * theta);
        (a, b, (det::sin(theta)/(2.0*theta) - a)/z, (a - 3.0*b)/(2.0*z))
    };
    let h = hat(p);
    let dh = hat(v);
    let h2 = h.compose(h);
    let left = dh.compose(h);
    let right = h.compose(dh);
    let result = Mat3 { m: std::array::from_fn(|i|
        -az * dz * h.m[i] - a * dh.m[i] + bz * dz * h2.m[i]
            + b * (left.m[i] + right.m[i])) };
    finite(&result.m)?;
    Ok(result)
}

/// d J_l(phi)[direction] = d J_r(-phi)[-direction].
///
/// # Errors
/// Same admission as [`so3_right_jacobian_directional`].
pub fn so3_left_jacobian_directional(
    phi: So3Tangent, direction: So3Tangent,
) -> Result<Mat3, GaError> {
    so3_right_jacobian_directional(phi.scale(-1.0), direction.scale(-1.0))
}

/// d J_r(phi)^-1[direction], retaining the canonical singularity checks.
///
/// # Errors
/// Refuses nonfinite data, the original inverse-Jacobian singularities, or overflow.
pub fn so3_right_jacobian_inverse_directional(
    phi: So3Tangent, direction: So3Tangent,
) -> Result<Mat3, GaError> {
    let inverse = so3_right_jacobian_inverse(phi)?;
    let derivative = so3_right_jacobian_directional(phi, direction)?;
    let product = inverse.compose(derivative).compose(inverse);
    let result = Mat3 { m: product.m.map(|x| -x) };
    finite(&result.m)?;
    Ok(result)
}

/// d J_l(phi)^-1[direction], with no finite-difference probe.
///
/// # Errors
/// Same admission as [`so3_right_jacobian_inverse_directional`].
pub fn so3_left_jacobian_inverse_directional(
    phi: So3Tangent, direction: So3Tangent,
) -> Result<Mat3, GaError> {
    so3_right_jacobian_inverse_directional(phi.scale(-1.0), direction.scale(-1.0))
}

fn hat(v: Vec3) -> Mat3 {
    Mat3 { m: [0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0] }
}
fn finite(values: &[f64]) -> Result<(), GaError> {
    for (index, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(GaError::NonFinite { context: "SO(3) Jacobian differential", index });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{so3_right_jacobian, so3_left_jacobian_inverse};
    fn tangent(v: [f64; 3]) -> So3Tangent { So3Tangent::new(Vec3::new(v[0], v[1], v[2])) }

    #[test]
    fn g0_zero_angle_differentials_have_exact_analytic_limits() {
        let zero = tangent([0.0; 3]);
        let v = tangent([0.3, -0.4, 0.8]);
        let skew = hat(v.angular);
        let right = so3_right_jacobian_directional(zero, v).unwrap();
        let inverse = so3_right_jacobian_inverse_directional(zero, v).unwrap();
        for i in 0..9 {
            assert_eq!(right.m[i], -0.5 * skew.m[i]);
            assert_eq!(inverse.m[i], 0.5 * skew.m[i]);
        }
    }

    #[test]
    fn g3_jacobian_differentials_match_independent_centered_probes() {
        let direction = tangent([0.23, -0.37, 0.51]);
        for scale in [0.0, 1e-8, 1e-4, 0.02, 0.3, 1.0, 3.13] {
            let p = tangent([0.36*scale, -0.48*scale, 0.8*scale]);
            let h = 2e-5;
            let plus = So3Tangent::new(Vec3::new(p.angular.x + h*direction.angular.x,
                p.angular.y + h*direction.angular.y, p.angular.z + h*direction.angular.z));
            let minus = So3Tangent::new(Vec3::new(p.angular.x - h*direction.angular.x,
                p.angular.y - h*direction.angular.y, p.angular.z - h*direction.angular.z));
            for (exact, upper, lower) in [
                (so3_right_jacobian_directional(p, direction).unwrap(),
                    so3_right_jacobian(plus).unwrap(), so3_right_jacobian(minus).unwrap()),
                (so3_left_jacobian_inverse_directional(p, direction).unwrap(),
                    so3_left_jacobian_inverse(plus).unwrap(), so3_left_jacobian_inverse(minus).unwrap()),
            ] {
                for i in 0..9 {
                    let probe = (upper.m[i] - lower.m[i])/(2.0*h);
                    assert!((exact.m[i] - probe).abs() < 3e-6, "angle {scale}, entry {i}: {} != {probe}", exact.m[i]);
                }
            }
        }
    }

    #[test]
    fn g0_inverse_identity_and_refusals() {
        let p = tangent([0.2, -0.3, 0.4]);
        let v = tangent([0.4, 0.5, -0.1]);
        let j = so3_right_jacobian(p).unwrap();
        let ji = so3_right_jacobian_inverse(p).unwrap();
        let a = so3_right_jacobian_directional(p, v).unwrap().compose(ji);
        let b = j.compose(so3_right_jacobian_inverse_directional(p, v).unwrap());
        for i in 0..9 { assert!((a.m[i]+b.m[i]).abs() < 1e-14); }
        assert!(so3_right_jacobian_directional(p, tangent([f64::NAN,0.0,0.0])).is_err());
        assert!(so3_right_jacobian_directional(tangent([f64::MAX,0.0,0.0]), v).is_err());
        assert!(so3_right_jacobian_inverse_directional(tangent([core::f64::consts::TAU,0.0,0.0]), v).is_err());
    }
}

//! Façade doctrine (plan §7.7): GA is internal excellence — the public
//! boundary speaks Vec3 / quaternion / matrix so no caller pays a
//! formalism tax. Conversions are component relabelings where possible
//! (bitwise-exact) and tight-ULP products elsewhere; both directions are
//! conformance-tested.

use crate::GaError;
use crate::mv::Pga;
use crate::pga::{Motor, Point};
use fs_math::det;

const E12: usize = 0b0110;
const E13: usize = 0b1010;
const E23: usize = 0b1100;
const E01: usize = 0b0011;
const E02: usize = 0b0101;
const E03: usize = 0b1001;

/// A plain 3-vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl core::ops::Add for Vec3 {
    type Output = Vec3;
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl core::ops::Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Vec3 {
    /// Construct.
    #[must_use]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Uniform scale.
    #[must_use]
    pub fn scale(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product with both product roundoffs retained in each component.
    ///
    /// Near-parallel force and arm vectors can have a small, physically real
    /// moment even when their rounded products agree. Symmetric compensation
    /// retains that moment and antisymmetry, including `v.cross(v) == 0` for
    /// finite products. A one-sided FMA alone would create a spurious moment
    /// in that last case. No geometric tolerance or torque clamp is applied.
    /// Overflow and non-finite products retain ordinary IEEE propagation.
    #[must_use]
    pub fn cross(self, o: Vec3) -> Vec3 {
        #[inline]
        fn difference_of_products(a: f64, b: f64, c: f64, d: f64) -> f64 {
            let ab = a * b;
            let cd = c * d;
            let difference = ab - cd;
            if !difference.is_finite() {
                return difference;
            }
            // TwoProd tails, symmetrically: swapping the operands negates
            // both terms. When the products nearly cancel their subtraction
            // is exact (Sterbenz), so the tails recover the resolved moment.
            difference + (a.mul_add(b, -ab) - c.mul_add(d, -cd))
        }
        Vec3::new(
            difference_of_products(self.y, o.z, self.z, o.y),
            difference_of_products(self.z, o.x, self.x, o.z),
            difference_of_products(self.x, o.y, self.y, o.x),
        )
    }

    /// Euclidean norm.
    #[must_use]
    pub fn norm(self) -> f64 {
        det::sqrt(self.dot(self))
    }
}

/// A quaternion `w + x i + y j + z k` (the conventional rotation façade).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    /// Scalar part.
    pub w: f64,
    /// i coefficient.
    pub x: f64,
    /// j coefficient.
    pub y: f64,
    /// k coefficient.
    pub z: f64,
}

impl Quat {
    /// The identity rotation.
    #[must_use]
    pub fn identity() -> Self {
        Quat {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Rotation by `angle` about the unit `axis` (deterministic trig).
    #[must_use]
    pub fn from_axis_angle(axis: Vec3, angle: f64) -> Self {
        let (s, c) = (det::sin(angle / 2.0), det::cos(angle / 2.0));
        Quat {
            w: c,
            x: s * axis.x,
            y: s * axis.y,
            z: s * axis.z,
        }
    }

    /// Conjugate (inverse for unit quaternions).
    #[must_use]
    pub fn conjugate(self) -> Quat {
        Quat {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl core::ops::Mul for Quat {
    type Output = Quat;
    /// Hamilton product.
    fn mul(self, o: Quat) -> Quat {
        Quat {
            w: self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
            x: self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
            y: self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
            z: self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
        }
    }
}

impl Quat {
    /// Rotate a vector (the hand-written baseline path the GA motor is
    /// benchmarked against).
    ///
    /// Evaluate the projective action `q v q⁻¹`, including the norm-squared
    /// divisor. A finite, nonzero quaternion is required. This also prevents
    /// an admitted unit-norm roundoff defect from scaling physical vectors;
    /// the stored quaternion is never rescaled.
    #[must_use]
    pub fn rotate(self, v: Vec3) -> Vec3 {
        let q = Vec3::new(self.x, self.y, self.z);
        let norm_sq = self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z;
        let t = q.cross(v).scale(2.0 / norm_sq);
        v + t.scale(self.w) + q.cross(t)
    }

    /// Exact relabeling into the PGA rotor `w − x e23 − y e31 − z e12`.
    #[must_use]
    pub fn to_rotor(self) -> Motor {
        let mut v = Pga::zero();
        v.0[0] = self.w;
        v.0[E23] = -self.x;
        v.0[E13] = self.y; // −y e31 = +y e13
        v.0[E12] = -self.z;
        Motor(v)
    }

    /// Exact relabeling back from a rotor's scalar + Euclidean bivector
    /// components (inverse of [`Quat::to_rotor`], bitwise).
    #[must_use]
    pub fn from_rotor(m: &Motor) -> Quat {
        Quat {
            w: m.0.0[0],
            x: -m.0.0[E23],
            y: m.0.0[E13],
            z: -m.0.0[E12],
        }
    }
}

impl Motor {
    /// Build the motor for "rotate by `q`, then translate by `t`".
    #[must_use]
    pub fn from_parts(q: Quat, t: Vec3) -> Motor {
        Motor::translator(t.x, t.y, t.z).compose(&q.to_rotor())
    }

    /// Recover `(q, t)`. The quaternion is a bitwise-exact relabeling of
    /// the motor's scalar + Euclidean-bivector components; the
    /// translation is recovered by one motor product divided by the rotor's
    /// norm squared (tight ULP). This is the same projective scale cancellation
    /// as the point sandwich, not normalization of the stored motor.
    #[must_use]
    pub fn to_parts(&self) -> (Quat, Vec3) {
        let q = Quat::from_rotor(self);
        let translator = self.compose(&q.to_rotor().reverse());
        let scale = -2.0 / (q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
        let t = Vec3::new(
            scale * translator.0.0[E01],
            scale * translator.0.0[E02],
            scale * translator.0.0[E03],
        );
        (q, t)
    }
}

/// A 3×4 rigid-motion matrix (rotation columns + translation), the bulk
/// transform façade: compose in motor land, apply in matrix land.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat34 {
    /// Row-major 3×4 entries.
    pub m: [f64; 12],
}

impl Mat34 {
    /// Lower a motor to its rigid-motion matrix.
    ///
    /// # Errors
    /// Propagates [`GaError::IdealPoint`] for non-motor input.
    pub fn from_motor(motor: &Motor) -> Result<Mat34, GaError> {
        let o = motor.transform_point(Point::new(0.0, 0.0, 0.0))?;
        let cols = [
            motor.transform_point(Point::new(1.0, 0.0, 0.0))?,
            motor.transform_point(Point::new(0.0, 1.0, 0.0))?,
            motor.transform_point(Point::new(0.0, 0.0, 1.0))?,
        ];
        let mut m = [0.0f64; 12];
        for (c, p) in cols.iter().enumerate() {
            m[c] = p.x - o.x;
            m[4 + c] = p.y - o.y;
            m[8 + c] = p.z - o.z;
        }
        m[3] = o.x;
        m[7] = o.y;
        m[11] = o.z;
        Ok(Mat34 { m })
    }

    /// Apply to a point.
    #[must_use]
    pub fn apply(&self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.m[0] * v.x + self.m[1] * v.y + self.m[2] * v.z + self.m[3],
            self.m[4] * v.x + self.m[5] * v.y + self.m[6] * v.z + self.m[7],
            self.m[8] * v.x + self.m[9] * v.y + self.m[10] * v.z + self.m[11],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::Vec3;

    #[test]
    fn g0_force_moments_retain_product_tails_and_antisymmetry() {
        // Exact integer oracle, independent of floating-point products.
        // The central case is (2^27 + 1)(2^27 - 1) - 2^54 = -1 Nm.
        for exponent in [27, 32, 40] {
            let n = 1_i128 << exponent;
            for offset in [1_i128, 3, 7, 11] {
                let arm = [n + offset, n, n - offset];
                let force = [n, n - offset, n + 2 * offset];
                let expected = core::array::from_fn::<_, 3, _>(|i| {
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    (arm[j] * force[k] - arm[k] * force[j]) as f64
                });
                let r = Vec3::new(arm[0] as f64, arm[1] as f64, arm[2] as f64);
                let f = Vec3::new(force[0] as f64, force[1] as f64, force[2] as f64);
                let torque = r.cross(f);
                assert_eq!(torque, Vec3::new(expected[0], expected[1], expected[2]));
                assert_eq!(f.cross(r), torque.scale(-1.0));
                assert_eq!(r.cross(r), Vec3::default());
            }
        }
    }

    #[test]
    fn g3_torque_preserves_binary_unit_rescaling_and_force_reversal() {
        let n = 134_217_728.0;
        let arm = Vec3::new(n + 1.0, n, 0.0);
        let force = Vec3::new(n, n - 1.0, 0.0);
        for (length_scale, force_scale) in [(1.0, 1.0), (0.125, 32.0), (1024.0, 0.25)] {
            let torque = arm.scale(length_scale).cross(force.scale(force_scale));
            assert_eq!(torque, Vec3::new(0.0, 0.0, -length_scale * force_scale));
            assert_eq!(
                arm.scale(length_scale).cross(force.scale(-force_scale)),
                torque.scale(-1.0)
            );
        }
    }

    #[test]
    fn cross_preserves_ieee_nonfinite_propagation() {
        let x = Vec3::new(f64::MAX, 0.0, 0.0);
        let y = Vec3::new(0.0, 2.0, 0.0);
        assert_eq!(x.cross(y).z, f64::INFINITY);
        assert_eq!(y.cross(x).z, f64::NEG_INFINITY);
        assert!(Vec3::new(f64::NAN, 0.0, 0.0).cross(y).z.is_nan());
    }
}

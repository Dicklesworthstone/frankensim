//! Lie-group façades for the existing PGA kinematics owner.
//!
//! This module does not introduce another pose representation. [`So3`] is a
//! validated, canonical façade over [`Quat`], and [`Se3`] is a validated,
//! canonical façade over [`Motor`]. The fixed-size matrices below are boundary
//! views for differentials and frame changes only.
//!
//! Conventions are explicit throughout:
//! - twists use `[angular, linear]` ordering;
//! - wrenches use the dual `[torque, force]` ordering;
//! - `space_*` means a left perturbation, `Exp(delta) * group`;
//! - `body_*` means a right perturbation, `group * Exp(delta)`;
//! - a pose maps coordinates from its local/body frame into its parent/space
//!   frame, so its adjoint maps body-coordinate twists to space coordinates.

use crate::GaError;
use crate::facade::{Quat, Vec3};
use crate::mv::Pga;
use crate::pga::{
    EVEN_BLADES, Motor, Point, axis_bivector, exp_bivector, ideal_bivector, motor_log,
};
use fs_math::det;

/// Accepted defect for a quaternion or PGA motor supplied at the validated
/// boundary. Accepted values are never rescaled; this is a check, not hidden
/// normalization.
pub const UNIT_TOLERANCE: f64 = 4.0e-12;

const FORBIDDEN_COMPONENT_TOLERANCE: f64 = 0.0;
const PURE_ROTATION_TRANSLATION_TOLERANCE: f64 = 0.0;
const DEGENERATE_NORM_SQUARED: f64 = f64::MIN_POSITIVE;
const SMALL_ANGLE_SQUARED: f64 = 1.0e-12;
const JACOBIAN_SINGULAR_SINE: f64 = 1.0e-12;
const SE3_SERIES_MAX_TERMS: usize = 96;
const SE3_SERIES_TOLERANCE: f64 = 64.0 * f64::EPSILON;
const SE3_SERIES_MAX_AD_NORM: f64 = 32.0;

/// Row-major 3×3 boundary matrix.
///
/// This is intentionally not a rotation/pose authority. Use [`So3`] for group
/// operations and this type only for linear maps and interoperability.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    /// Row-major entries.
    pub m: [f64; 9],
}

impl Mat3 {
    /// Identity matrix.
    #[must_use]
    pub const fn identity() -> Self {
        Self {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        }
    }

    /// Zero matrix.
    #[must_use]
    pub const fn zero() -> Self {
        Self { m: [0.0; 9] }
    }

    /// Apply the linear map to a vector.
    #[must_use]
    pub fn apply(self, value: Vec3) -> Vec3 {
        Vec3::new(
            self.m[0] * value.x + self.m[1] * value.y + self.m[2] * value.z,
            self.m[3] * value.x + self.m[4] * value.y + self.m[5] * value.z,
            self.m[6] * value.x + self.m[7] * value.y + self.m[8] * value.z,
        )
    }

    /// Matrix product, with `rhs` applied first.
    #[must_use]
    pub fn compose(self, rhs: Self) -> Self {
        let mut out = [0.0; 9];
        let mut row = 0;
        while row < 3 {
            let mut col = 0;
            while col < 3 {
                let mut inner = 0;
                while inner < 3 {
                    out[row * 3 + col] += self.m[row * 3 + inner] * rhs.m[inner * 3 + col];
                    inner += 1;
                }
                col += 1;
            }
            row += 1;
        }
        Self { m: out }
    }

    /// Matrix transpose.
    #[must_use]
    pub fn transpose(self) -> Self {
        Self {
            m: [
                self.m[0], self.m[3], self.m[6], self.m[1], self.m[4], self.m[7], self.m[2],
                self.m[5], self.m[8],
            ],
        }
    }

    /// Largest absolute entry.
    #[must_use]
    pub fn max_abs(self) -> f64 {
        max_abs_slice(&self.m)
    }

    fn add_scaled(self, rhs: Self, scale: f64) -> Self {
        let mut out = self.m;
        let mut i = 0;
        while i < 9 {
            out[i] += scale * rhs.m[i];
            i += 1;
        }
        Self { m: out }
    }
}

/// Row-major 6×6 boundary matrix for twists and wrenches.
///
/// Coordinates are ordered `[angular, linear]` (or dually
/// `[torque, force]`). It is not a pose authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat6 {
    /// Row-major entries.
    pub m: [f64; 36],
}

impl Mat6 {
    /// Identity matrix.
    #[must_use]
    pub const fn identity() -> Self {
        let mut m = [0.0; 36];
        m[0] = 1.0;
        m[7] = 1.0;
        m[14] = 1.0;
        m[21] = 1.0;
        m[28] = 1.0;
        m[35] = 1.0;
        Self { m }
    }

    /// Zero matrix.
    #[must_use]
    pub const fn zero() -> Self {
        Self { m: [0.0; 36] }
    }

    /// Matrix product, with `rhs` applied first.
    #[must_use]
    pub fn compose(self, rhs: Self) -> Self {
        let mut out = [0.0; 36];
        let mut row = 0;
        while row < 6 {
            let mut col = 0;
            while col < 6 {
                let mut inner = 0;
                while inner < 6 {
                    out[row * 6 + col] += self.m[row * 6 + inner] * rhs.m[inner * 6 + col];
                    inner += 1;
                }
                col += 1;
            }
            row += 1;
        }
        Self { m: out }
    }

    /// Matrix transpose.
    #[must_use]
    pub fn transpose(self) -> Self {
        let mut out = [0.0; 36];
        let mut row = 0;
        while row < 6 {
            let mut col = 0;
            while col < 6 {
                out[row * 6 + col] = self.m[col * 6 + row];
                col += 1;
            }
            row += 1;
        }
        Self { m: out }
    }

    /// Apply to a twist.
    #[must_use]
    pub fn apply_twist(self, value: Twist) -> Twist {
        let out = self.apply_array(value.to_array());
        Twist::new(
            Vec3::new(out[0], out[1], out[2]),
            Vec3::new(out[3], out[4], out[5]),
        )
    }

    /// Apply to a wrench in `[torque, force]` ordering.
    #[must_use]
    pub fn apply_wrench(self, value: Wrench) -> Wrench {
        let out = self.apply_array(value.to_array());
        Wrench::new(
            Vec3::new(out[0], out[1], out[2]),
            Vec3::new(out[3], out[4], out[5]),
        )
    }

    /// Largest absolute entry.
    #[must_use]
    pub fn max_abs(self) -> f64 {
        max_abs_slice(&self.m)
    }

    fn apply_array(self, value: [f64; 6]) -> [f64; 6] {
        let mut out = [0.0; 6];
        let mut row = 0;
        while row < 6 {
            let mut col = 0;
            while col < 6 {
                out[row] += self.m[row * 6 + col] * value[col];
                col += 1;
            }
            row += 1;
        }
        out
    }

    fn add_scaled(self, rhs: Self, scale: f64) -> Self {
        let mut out = self.m;
        let mut i = 0;
        while i < 36 {
            out[i] += scale * rhs.m[i];
            i += 1;
        }
        Self { m: out }
    }

    fn norm_infinity(self) -> f64 {
        let mut maximum = 0.0_f64;
        let mut row = 0;
        while row < 6 {
            let mut sum = 0.0;
            let mut col = 0;
            while col < 6 {
                sum += self.m[row * 6 + col].abs();
                col += 1;
            }
            maximum = maximum.max(sum);
            row += 1;
        }
        maximum
    }
}

/// Tangent vector in `so(3)`, represented by an angular rotation vector in
/// radians.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct So3Tangent {
    /// Angular coordinates in radians.
    pub angular: Vec3,
}

impl So3Tangent {
    /// Construct an `so(3)` tangent.
    #[must_use]
    pub const fn new(angular: Vec3) -> Self {
        Self { angular }
    }

    /// Uniform scaling in the Lie algebra.
    #[must_use]
    pub fn scale(self, scale: f64) -> Self {
        Self::new(self.angular.scale(scale))
    }
}

/// A spatial/body twist in `[angular, linear]` ordering.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Twist {
    /// Angular velocity/rotation-vector part.
    pub angular: Vec3,
    /// Linear velocity/translation-generator part.
    pub linear: Vec3,
}

impl Twist {
    /// Construct a twist using the crate-wide `[angular, linear]` convention.
    #[must_use]
    pub const fn new(angular: Vec3, linear: Vec3) -> Self {
        Self { angular, linear }
    }

    /// Zero twist.
    #[must_use]
    pub const fn zero() -> Self {
        Self::new(
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        )
    }

    /// Uniform scaling in `se(3)`.
    #[must_use]
    pub fn scale(self, scale: f64) -> Self {
        Self::new(self.angular.scale(scale), self.linear.scale(scale))
    }

    /// Vector-space sum.
    #[must_use]
    pub fn plus(self, rhs: Self) -> Self {
        Self::new(self.angular + rhs.angular, self.linear + rhs.linear)
    }

    /// Vector-space difference.
    #[must_use]
    pub fn minus(self, rhs: Self) -> Self {
        Self::new(self.angular - rhs.angular, self.linear - rhs.linear)
    }

    /// Lie bracket `[self, rhs]` in `[angular, linear]` coordinates.
    #[must_use]
    pub fn bracket(self, rhs: Self) -> Self {
        Self::new(
            self.angular.cross(rhs.angular),
            self.angular.cross(rhs.linear) + self.linear.cross(rhs.angular),
        )
    }

    /// The adjoint-representation matrix `ad_self`, satisfying
    /// `ad_self * rhs = [self, rhs]`.
    #[must_use]
    pub fn ad(self) -> Mat6 {
        let angular_hat = hat(self.angular);
        let linear_hat = hat(self.linear);
        let mut out = Mat6::zero();
        set_mat3_block(&mut out, 0, 0, angular_hat);
        set_mat3_block(&mut out, 1, 0, linear_hat);
        set_mat3_block(&mut out, 1, 1, angular_hat);
        out
    }

    /// Transform body coordinates into space coordinates with a pose adjoint.
    #[must_use]
    pub fn transform_by(self, pose: &Se3) -> Self {
        pose.adjoint().apply_twist(self)
    }

    /// Coordinates as `[wx, wy, wz, vx, vy, vz]`.
    #[must_use]
    pub const fn to_array(self) -> [f64; 6] {
        [
            self.angular.x,
            self.angular.y,
            self.angular.z,
            self.linear.x,
            self.linear.y,
            self.linear.z,
        ]
    }
}

/// A wrench dual to [`Twist`], in `[torque, force]` ordering.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Wrench {
    /// Torque/moment part.
    pub torque: Vec3,
    /// Force part.
    pub force: Vec3,
}

impl Wrench {
    /// Construct a wrench using `[torque, force]` ordering.
    #[must_use]
    pub const fn new(torque: Vec3, force: Vec3) -> Self {
        Self { torque, force }
    }

    /// Dual pairing `torque·angular + force·linear` (instantaneous power
    /// when the coordinates carry physical rates and forces).
    #[must_use]
    pub fn pairing(self, twist: Twist) -> f64 {
        self.torque.dot(twist.angular) + self.force.dot(twist.linear)
    }

    /// Transform body coordinates into space coordinates with the coadjoint
    /// `Ad_pose^{-T}`, preserving [`Wrench::pairing`].
    #[must_use]
    pub fn transform_by(self, pose: &Se3) -> Self {
        pose.coadjoint().apply_wrench(self)
    }

    /// Coordinates as `[tx, ty, tz, fx, fy, fz]`.
    #[must_use]
    pub const fn to_array(self) -> [f64; 6] {
        [
            self.torque.x,
            self.torque.y,
            self.torque.z,
            self.force.x,
            self.force.y,
            self.force.z,
        ]
    }
}

/// Validated, canonical `SO(3)` façade backed by the existing [`Quat`].
///
/// Construction checks unit norm but never silently normalizes. The quaternion
/// double cover is canonicalized to the nonnegative-scalar hemisphere, with a
/// deterministic lexicographic tie-break at exactly π.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct So3 {
    quat: Quat,
}

impl So3 {
    /// Identity rotation.
    #[must_use]
    pub fn identity() -> Self {
        Self {
            quat: Quat::identity(),
        }
    }

    /// Validate and canonicalize an existing quaternion without rescaling it.
    ///
    /// # Errors
    /// Refuses non-finite, degenerate, or non-unit input.
    pub fn try_from_quat(quat: Quat) -> Result<Self, GaError> {
        validate_finite(&[quat.w, quat.x, quat.y, quat.z], "SO(3) quaternion")?;
        let norm_squared = quat.w * quat.w + quat.x * quat.x + quat.y * quat.y + quat.z * quat.z;
        if norm_squared <= DEGENERATE_NORM_SQUARED {
            return Err(GaError::DegenerateNorm {
                context: "SO(3) quaternion",
                norm_squared,
            });
        }
        let defect = (norm_squared - 1.0).abs();
        if defect > UNIT_TOLERANCE {
            return Err(GaError::NotUnit {
                context: "SO(3) quaternion",
                defect,
                tolerance: UNIT_TOLERANCE,
            });
        }
        Ok(Self {
            quat: canonical_quat(quat),
        })
    }

    /// Validate a PGA motor as a pure rotation.
    ///
    /// # Errors
    /// Refuses invalid motors and motors with non-negligible translation.
    pub fn try_from_motor(motor: Motor) -> Result<Self, GaError> {
        let pose = Se3::try_from_motor(motor)?;
        let translation = pose.translation();
        let defect = translation.norm();
        if defect > PURE_ROTATION_TRANSLATION_TOLERANCE {
            return Err(GaError::InvalidRepresentation {
                context: "SO(3) motor translation",
                defect,
                tolerance: PURE_ROTATION_TRANSLATION_TOLERANCE,
            });
        }
        Ok(pose.rotation())
    }

    /// Borrow the authoritative quaternion façade.
    #[must_use]
    pub const fn as_quat(&self) -> &Quat {
        &self.quat
    }

    /// Convert to the existing PGA rotor by bitwise component relabeling.
    #[must_use]
    pub fn to_motor(self) -> Motor {
        self.quat.to_rotor()
    }

    /// Lie exponential from an angular rotation vector.
    ///
    /// # Errors
    /// Refuses non-finite coordinates.
    pub fn exp(tangent: So3Tangent) -> Result<Self, GaError> {
        validate_vec3(tangent.angular, "SO(3) exponential")?;
        let theta_squared = tangent.angular.dot(tangent.angular);
        let (scalar, vector_scale) = if theta_squared < SMALL_ANGLE_SQUARED {
            let theta_fourth = theta_squared * theta_squared;
            (
                1.0 - theta_squared / 8.0 + theta_fourth / 384.0,
                0.5 - theta_squared / 48.0 + theta_fourth / 3840.0,
            )
        } else {
            let theta = det::sqrt(theta_squared);
            let half = 0.5 * theta;
            (det::cos(half), det::sin(half) / theta)
        };
        Self::try_from_quat(Quat {
            w: scalar,
            x: vector_scale * tangent.angular.x,
            y: vector_scale * tangent.angular.y,
            z: vector_scale * tangent.angular.z,
        })
    }

    /// Principal logarithm, with rotation-vector magnitude in `[0, π]`.
    #[must_use]
    pub fn log(self) -> So3Tangent {
        let vector = Vec3::new(self.quat.x, self.quat.y, self.quat.z);
        let vector_norm = vector.norm();
        let scale = if vector_norm < f64::EPSILON {
            // atan2(n, w) / n = 1/w - n²/(3w³) + O(n⁴).
            let inverse_w = 1.0 / self.quat.w;
            2.0 * (inverse_w - vector_norm * vector_norm * inverse_w.powi(3) / 3.0)
        } else {
            2.0 * det::atan2(vector_norm, self.quat.w) / vector_norm
        };
        So3Tangent::new(vector.scale(scale))
    }

    /// Compose rotations, applying `rhs` first.
    ///
    /// # Errors
    /// Refuses a result whose accumulated roundoff crosses the unit tolerance.
    pub fn compose(self, rhs: Self) -> Result<Self, GaError> {
        Self::try_from_quat(self.quat * rhs.quat)
    }

    /// Group inverse.
    #[must_use]
    pub fn inverse(self) -> Self {
        Self {
            quat: canonical_quat(self.quat.conjugate()),
        }
    }

    /// Rotate a finite vector.
    ///
    /// # Errors
    /// Refuses non-finite vector coordinates.
    pub fn rotate(self, value: Vec3) -> Result<Vec3, GaError> {
        validate_vec3(value, "SO(3) vector transform")?;
        Ok(self.quat.rotate(value))
    }

    /// Lower to a 3×3 matrix boundary view.
    #[must_use]
    pub fn matrix(self) -> Mat3 {
        rotation_matrix(self.quat)
    }

    /// Right/body perturbation: `self * Exp(delta_body)`.
    ///
    /// # Errors
    /// Propagates exponential or composition validation failures.
    pub fn body_plus(self, delta_body: So3Tangent) -> Result<Self, GaError> {
        self.compose(Self::exp(delta_body)?)
    }

    /// Left/space perturbation: `Exp(delta_space) * self`.
    ///
    /// # Errors
    /// Propagates exponential or composition validation failures.
    pub fn space_plus(self, delta_space: So3Tangent) -> Result<Self, GaError> {
        Self::exp(delta_space)?.compose(self)
    }

    /// Right/body difference `Log(reference^{-1} * self)`.
    ///
    /// # Errors
    /// Refuses accumulated group-invariant drift.
    pub fn body_minus(self, reference: Self) -> Result<So3Tangent, GaError> {
        Ok(reference.inverse().compose(self)?.log())
    }

    /// Left/space difference `Log(self * reference^{-1})`.
    ///
    /// # Errors
    /// Refuses accumulated group-invariant drift.
    pub fn space_minus(self, reference: Self) -> Result<So3Tangent, GaError> {
        Ok(self.compose(reference.inverse())?.log())
    }
}

impl TryFrom<Quat> for So3 {
    type Error = GaError;

    fn try_from(value: Quat) -> Result<Self, Self::Error> {
        Self::try_from_quat(value)
    }
}

/// Validated, canonical `SE(3)` façade backed by the existing PGA [`Motor`].
///
/// Construction checks even-grade and unit-motor invariants without silently
/// renormalizing. The motor double cover uses the same deterministic branch as
/// [`So3`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Se3 {
    motor: Motor,
}

impl Se3 {
    /// Identity rigid motion.
    #[must_use]
    pub fn identity() -> Self {
        Self {
            motor: Motor::identity(),
        }
    }

    /// Validate and canonicalize an existing PGA motor without rescaling it.
    ///
    /// # Errors
    /// Refuses non-finite, non-even, degenerate, or non-unit input.
    pub fn try_from_motor(motor: Motor) -> Result<Self, GaError> {
        validate_finite(&motor.0.0, "SE(3) motor")?;
        let mut forbidden_defect = 0.0_f64;
        let mut index = 0;
        while index < Pga::BLADES {
            if !EVEN_BLADES.contains(&index) {
                forbidden_defect = forbidden_defect.max(motor.0.0[index].abs());
            }
            index += 1;
        }
        if forbidden_defect > FORBIDDEN_COMPONENT_TOLERANCE {
            return Err(GaError::InvalidRepresentation {
                context: "SE(3) motor grade",
                defect: forbidden_defect,
                tolerance: FORBIDDEN_COMPONENT_TOLERANCE,
            });
        }
        let norm_squared = motor.0.gp(&motor.0.reverse()).scalar_part();
        if norm_squared <= DEGENERATE_NORM_SQUARED {
            return Err(GaError::DegenerateNorm {
                context: "SE(3) motor",
                norm_squared,
            });
        }
        let defect = motor.unit_defect();
        if !defect.is_finite() {
            return Err(GaError::NonFinite {
                context: "SE(3) motor unit defect",
                index: 0,
            });
        }
        if defect > UNIT_TOLERANCE {
            return Err(GaError::NotUnit {
                context: "SE(3) motor",
                defect,
                tolerance: UNIT_TOLERANCE,
            });
        }
        let quat = Quat::from_rotor(&motor);
        let mut motor = if canonical_quat_needs_flip(quat) {
            Motor(motor.0.scale(-1.0))
        } else {
            motor
        };
        for coordinate in &mut motor.0.0 {
            *coordinate = positive_zero(*coordinate);
        }
        Ok(Self { motor })
    }

    /// Construct from validated rotation and finite translation.
    ///
    /// # Errors
    /// Refuses non-finite translation or an invalid generated motor.
    pub fn from_parts(rotation: So3, translation: Vec3) -> Result<Self, GaError> {
        validate_vec3(translation, "SE(3) translation")?;
        Self::try_from_motor(Motor::from_parts(*rotation.as_quat(), translation))
    }

    /// Borrow the authoritative PGA motor.
    #[must_use]
    pub const fn as_motor(&self) -> &Motor {
        &self.motor
    }

    /// Validated rotation component.
    #[must_use]
    pub fn rotation(self) -> So3 {
        So3 {
            quat: canonical_quat(Quat::from_rotor(&self.motor)),
        }
    }

    /// Translation component in the parent/space frame.
    #[must_use]
    pub fn translation(self) -> Vec3 {
        self.motor.to_parts().1
    }

    /// Lie exponential from a twist in `[angular, linear]` coordinates.
    ///
    /// This is the existing PGA bivector exponential under the explicit
    /// identification `B = -½(axis(angular) + ideal(linear))`.
    ///
    /// # Errors
    /// Refuses non-finite coordinates or a generated motor outside the
    /// validated unit tolerance.
    pub fn exp(twist: Twist) -> Result<Self, GaError> {
        validate_twist(twist, "SE(3) exponential")?;
        let bivector = axis_bivector(twist.angular.x, twist.angular.y, twist.angular.z)
            .add(&ideal_bivector(
                twist.linear.x,
                twist.linear.y,
                twist.linear.z,
            ))
            .scale(-0.5);
        Self::try_from_motor(exp_bivector(&bivector))
    }

    /// Principal logarithm in `[angular, linear]` coordinates.
    ///
    /// The existing PGA `motor_log` chooses the canonical motor branch, so the
    /// angular magnitude is in `[0, π]`.
    #[must_use]
    pub fn log(self) -> Twist {
        let bivector = motor_log(&self.motor);
        Twist::new(
            Vec3::new(
                -2.0 * bivector.0[0b1100],
                2.0 * bivector.0[0b1010],
                -2.0 * bivector.0[0b0110],
            ),
            Vec3::new(
                -2.0 * bivector.0[0b0011],
                -2.0 * bivector.0[0b0101],
                -2.0 * bivector.0[0b1001],
            ),
        )
    }

    /// Compose poses, applying `rhs` first.
    ///
    /// # Errors
    /// Refuses a result whose accumulated roundoff crosses the unit tolerance.
    pub fn compose(self, rhs: Self) -> Result<Self, GaError> {
        Self::try_from_motor(self.motor.compose(&rhs.motor))
    }

    /// Group inverse.
    ///
    /// # Errors
    /// Refuses a result outside the validated group invariants.
    pub fn inverse(self) -> Result<Self, GaError> {
        Self::try_from_motor(self.motor.reverse())
    }

    /// Transform a point from body/local coordinates to parent/space
    /// coordinates.
    ///
    /// # Errors
    /// Propagates an impossible ideal-point result.
    pub fn transform_point(self, point: Vec3) -> Result<Vec3, GaError> {
        validate_vec3(point, "SE(3) point transform")?;
        let out = self
            .motor
            .transform_point(Point::new(point.x, point.y, point.z))?;
        Ok(Vec3::new(out.x, out.y, out.z))
    }

    /// Transform a free finite vector (rotation only).
    ///
    /// # Errors
    /// Refuses non-finite vector coordinates.
    pub fn transform_vector(self, vector: Vec3) -> Result<Vec3, GaError> {
        self.rotation().rotate(vector)
    }

    /// Group adjoint mapping body-coordinate twists to space coordinates.
    #[must_use]
    pub fn adjoint(self) -> Mat6 {
        let rotation = self.rotation().matrix();
        let translation_cross_rotation = hat(self.translation()).compose(rotation);
        let mut out = Mat6::zero();
        set_mat3_block(&mut out, 0, 0, rotation);
        set_mat3_block(&mut out, 1, 0, translation_cross_rotation);
        set_mat3_block(&mut out, 1, 1, rotation);
        out
    }

    /// Coadjoint `Ad_self^{-T}` mapping body-coordinate wrenches to space
    /// coordinates while preserving the twist/wrench dual pairing.
    #[must_use]
    pub fn coadjoint(self) -> Mat6 {
        let rotation = self.rotation().matrix();
        let translation_cross_rotation = hat(self.translation()).compose(rotation);
        let mut out = Mat6::zero();
        set_mat3_block(&mut out, 0, 0, rotation);
        set_mat3_block(&mut out, 0, 1, translation_cross_rotation);
        set_mat3_block(&mut out, 1, 1, rotation);
        out
    }

    /// Right/body perturbation: `self * Exp(delta_body)`.
    ///
    /// # Errors
    /// Propagates exponential or composition validation failures.
    pub fn body_plus(self, delta_body: Twist) -> Result<Self, GaError> {
        self.compose(Self::exp(delta_body)?)
    }

    /// Left/space perturbation: `Exp(delta_space) * self`.
    ///
    /// # Errors
    /// Propagates exponential or composition validation failures.
    pub fn space_plus(self, delta_space: Twist) -> Result<Self, GaError> {
        Self::exp(delta_space)?.compose(self)
    }

    /// Right/body difference `Log(reference^{-1} * self)`.
    ///
    /// # Errors
    /// Refuses accumulated group-invariant drift.
    pub fn body_minus(self, reference: Self) -> Result<Twist, GaError> {
        Ok(reference.inverse()?.compose(self)?.log())
    }

    /// Left/space difference `Log(self * reference^{-1})`.
    ///
    /// # Errors
    /// Refuses accumulated group-invariant drift.
    pub fn space_minus(self, reference: Self) -> Result<Twist, GaError> {
        Ok(self.compose(reference.inverse()?)?.log())
    }

    /// Left-trivialized differential of `Exp` at `twist`, returned with an
    /// analytic tail certificate for the deterministic `ad` series.
    ///
    /// # Errors
    /// Refuses non-finite or excessively ill-conditioned coordinates, or a
    /// series that cannot certify its tail within the bounded term budget.
    pub fn left_jacobian(twist: Twist) -> Result<Se3Jacobian, GaError> {
        se3_jacobian_series(twist, 1.0, "SE(3) left Jacobian")
    }

    /// Right-trivialized differential of `Exp` at `twist`.
    ///
    /// # Errors
    /// Has the same structured refusal policy as [`Se3::left_jacobian`].
    pub fn right_jacobian(twist: Twist) -> Result<Se3Jacobian, GaError> {
        se3_jacobian_series(twist, -1.0, "SE(3) right Jacobian")
    }
}

impl TryFrom<Motor> for Se3 {
    type Error = GaError;

    fn try_from(value: Motor) -> Result<Self, Self::Error> {
        Self::try_from_motor(value)
    }
}

/// Deterministic `SE(3)` exponential differential plus its convergence
/// receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Se3Jacobian {
    /// The 6×6 differential in `[angular, linear]` ordering.
    pub matrix: Mat6,
    /// Number of accumulated series terms, including the identity term.
    pub terms_used: usize,
    /// Analytic infinity-norm upper bound on the unaccumulated tail.
    pub tail_bound: f64,
    /// Infinity norm of the `ad` matrix used for conditioning/refusal.
    pub ad_norm: f64,
}

impl Se3Jacobian {
    /// Apply the differential to a twist.
    #[must_use]
    pub fn apply(self, twist: Twist) -> Twist {
        self.matrix.apply_twist(twist)
    }
}

/// `SO(3)` left Jacobian `J_l(phi)` with analytic small-angle limits.
///
/// # Errors
/// Refuses non-finite angular coordinates.
pub fn so3_left_jacobian(phi: So3Tangent) -> Result<Mat3, GaError> {
    so3_jacobian(phi.angular, 1.0)
}

/// `SO(3)` right Jacobian `J_r(phi) = J_l(-phi)`.
///
/// # Errors
/// Refuses non-finite angular coordinates.
pub fn so3_right_jacobian(phi: So3Tangent) -> Result<Mat3, GaError> {
    so3_jacobian(phi.angular, -1.0)
}

/// Inverse `SO(3)` left Jacobian with analytic zero-angle limit and an
/// explicit refusal at nonzero `2π` singularities.
///
/// # Errors
/// Refuses non-finite coordinates or a singular/ill-conditioned branch.
pub fn so3_left_jacobian_inverse(phi: So3Tangent) -> Result<Mat3, GaError> {
    so3_jacobian_inverse(phi.angular, 1.0)
}

/// Inverse `SO(3)` right Jacobian.
///
/// # Errors
/// Refuses non-finite coordinates or a singular/ill-conditioned branch.
pub fn so3_right_jacobian_inverse(phi: So3Tangent) -> Result<Mat3, GaError> {
    so3_jacobian_inverse(phi.angular, -1.0)
}

fn so3_jacobian(phi: Vec3, handedness: f64) -> Result<Mat3, GaError> {
    validate_vec3(phi, "SO(3) Jacobian")?;
    let theta_squared = phi.dot(phi);
    let (a, b) = if theta_squared < SMALL_ANGLE_SQUARED {
        let theta_fourth = theta_squared * theta_squared;
        (
            0.5 - theta_squared / 24.0 + theta_fourth / 720.0,
            1.0 / 6.0 - theta_squared / 120.0 + theta_fourth / 5040.0,
        )
    } else {
        let theta = det::sqrt(theta_squared);
        (
            (1.0 - det::cos(theta)) / theta_squared,
            (theta - det::sin(theta)) / (theta_squared * theta),
        )
    };
    let cross = hat(phi);
    Ok(Mat3::identity()
        .add_scaled(cross, handedness * a)
        .add_scaled(cross.compose(cross), b))
}

fn so3_jacobian_inverse(phi: Vec3, handedness: f64) -> Result<Mat3, GaError> {
    validate_vec3(phi, "SO(3) inverse Jacobian")?;
    let theta_squared = phi.dot(phi);
    let coefficient = if theta_squared < SMALL_ANGLE_SQUARED {
        let theta_fourth = theta_squared * theta_squared;
        1.0 / 12.0 + theta_squared / 720.0 + theta_fourth / 30_240.0
    } else {
        let theta = det::sqrt(theta_squared);
        let half = 0.5 * theta;
        let sine_half = det::sin(half);
        if sine_half.abs() < JACOBIAN_SINGULAR_SINE {
            return Err(GaError::IllConditioned {
                context: "SO(3) inverse Jacobian",
                measure: sine_half.abs(),
                limit: JACOBIAN_SINGULAR_SINE,
            });
        }
        (1.0 - half * det::cos(half) / sine_half) / theta_squared
    };
    let cross = hat(phi);
    Ok(Mat3::identity()
        .add_scaled(cross, -0.5 * handedness)
        .add_scaled(cross.compose(cross), coefficient))
}

fn se3_jacobian_series(
    twist: Twist,
    handedness: f64,
    context: &'static str,
) -> Result<Se3Jacobian, GaError> {
    validate_twist(twist, context)?;
    let ad = twist.ad();
    let mut signed_ad = Mat6::zero();
    let mut i = 0;
    while i < 36 {
        signed_ad.m[i] = handedness * ad.m[i];
        i += 1;
    }
    let ad_norm = signed_ad.norm_infinity();
    if ad_norm > SE3_SERIES_MAX_AD_NORM {
        return Err(GaError::IllConditioned {
            context,
            measure: ad_norm,
            limit: SE3_SERIES_MAX_AD_NORM,
        });
    }
    if ad_norm == 0.0 {
        return Ok(Se3Jacobian {
            matrix: Mat6::identity(),
            terms_used: 1,
            tail_bound: 0.0,
            ad_norm,
        });
    }

    // dexp = Σ ad^n/(n+1)!. The scalar majorant uses the submultiplicative
    // infinity norm. Once r = ||ad||/(n+2) < 1, the remaining tail is bounded
    // by next_term/(1-r), giving a deterministic convergence receipt.
    let mut sum = Mat6::identity();
    let mut power = Mat6::identity();
    let mut coefficient = 1.0;
    let mut term_bound = 1.0;
    let mut last_tail_bound = f64::INFINITY;
    let mut n = 1;
    while n < SE3_SERIES_MAX_TERMS {
        power = power.compose(signed_ad);
        coefficient /= (n + 1) as f64;
        term_bound *= ad_norm / (n + 1) as f64;
        sum = sum.add_scaled(power, coefficient);

        let ratio = ad_norm / (n + 2) as f64;
        if ratio < 1.0 {
            let next_term_bound = term_bound * ratio;
            last_tail_bound = next_term_bound / (1.0 - ratio);
            let target = SE3_SERIES_TOLERANCE * sum.norm_infinity().max(1.0);
            if last_tail_bound <= target {
                return Ok(Se3Jacobian {
                    matrix: sum,
                    terms_used: n + 1,
                    tail_bound: last_tail_bound,
                    ad_norm,
                });
            }
        }
        n += 1;
    }
    Err(GaError::SeriesDidNotConverge {
        context,
        terms: SE3_SERIES_MAX_TERMS,
        tail_bound: last_tail_bound,
    })
}

fn canonical_quat(quat: Quat) -> Quat {
    let canonical = if canonical_quat_needs_flip(quat) {
        Quat {
            w: -quat.w,
            x: -quat.x,
            y: -quat.y,
            z: -quat.z,
        }
    } else {
        quat
    };
    Quat {
        w: positive_zero(canonical.w),
        x: positive_zero(canonical.x),
        y: positive_zero(canonical.y),
        z: positive_zero(canonical.z),
    }
}

fn positive_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn canonical_quat_needs_flip(quat: Quat) -> bool {
    quat.w < 0.0
        || (quat.w == 0.0
            && (quat.x < 0.0
                || (quat.x == 0.0 && (quat.y < 0.0 || (quat.y == 0.0 && quat.z < 0.0)))))
}

fn rotation_matrix(quat: Quat) -> Mat3 {
    let xx = quat.x * quat.x;
    let yy = quat.y * quat.y;
    let zz = quat.z * quat.z;
    let xy = quat.x * quat.y;
    let xz = quat.x * quat.z;
    let yz = quat.y * quat.z;
    let wx = quat.w * quat.x;
    let wy = quat.w * quat.y;
    let wz = quat.w * quat.z;
    Mat3 {
        m: [
            1.0 - 2.0 * (yy + zz),
            2.0 * (xy - wz),
            2.0 * (xz + wy),
            2.0 * (xy + wz),
            1.0 - 2.0 * (xx + zz),
            2.0 * (yz - wx),
            2.0 * (xz - wy),
            2.0 * (yz + wx),
            1.0 - 2.0 * (xx + yy),
        ],
    }
}

fn hat(vector: Vec3) -> Mat3 {
    Mat3 {
        m: [
            0.0, -vector.z, vector.y, vector.z, 0.0, -vector.x, -vector.y, vector.x, 0.0,
        ],
    }
}

fn set_mat3_block(matrix: &mut Mat6, block_row: usize, block_col: usize, block: Mat3) {
    let mut row = 0;
    while row < 3 {
        let mut col = 0;
        while col < 3 {
            matrix.m[(block_row * 3 + row) * 6 + block_col * 3 + col] = block.m[row * 3 + col];
            col += 1;
        }
        row += 1;
    }
}

fn validate_vec3(value: Vec3, context: &'static str) -> Result<(), GaError> {
    validate_finite(&[value.x, value.y, value.z], context)
}

fn validate_twist(value: Twist, context: &'static str) -> Result<(), GaError> {
    validate_finite(&value.to_array(), context)
}

fn validate_finite(values: &[f64], context: &'static str) -> Result<(), GaError> {
    let mut index = 0;
    while index < values.len() {
        if !values[index].is_finite() {
            return Err(GaError::NonFinite { context, index });
        }
        index += 1;
    }
    Ok(())
}

fn max_abs_slice(values: &[f64]) -> f64 {
    let mut maximum = 0.0_f64;
    let mut index = 0;
    while index < values.len() {
        maximum = maximum.max(values[index].abs());
        index += 1;
    }
    maximum
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn vec_error(lhs: Vec3, rhs: Vec3) -> f64 {
        (lhs - rhs).norm()
    }

    fn twist_error(lhs: Twist, rhs: Twist) -> f64 {
        vec_error(lhs.angular, rhs.angular).max(vec_error(lhs.linear, rhs.linear))
    }

    fn matrix3_identity_error(matrix: Mat3) -> f64 {
        let mut residual = matrix.m;
        residual[0] -= 1.0;
        residual[4] -= 1.0;
        residual[8] -= 1.0;
        max_abs_slice(&residual)
    }

    fn matrix6_error(lhs: &Mat6, rhs: &Mat6) -> f64 {
        let mut residual = [0.0; 36];
        let mut i = 0;
        while i < 36 {
            residual[i] = lhs.m[i] - rhs.m[i];
            i += 1;
        }
        max_abs_slice(&residual)
    }

    #[test]
    fn g0_so3_exp_log_zero_tiny_general_and_near_pi() {
        let cases = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0e-13, -2.0e-13, 3.0e-13),
            Vec3::new(0.3, -0.7, 1.1),
            Vec3::new(PI - 1.0e-10, 0.0, 0.0),
        ];
        for angular in cases {
            let tangent = So3Tangent::new(angular);
            let recovered = So3::exp(tangent).unwrap().log();
            assert!(
                vec_error(recovered.angular, angular) < 2.0e-10,
                "SO(3) round trip failed for {angular:?}: {recovered:?}"
            );
        }
    }

    #[test]
    fn g0_se3_exp_log_zero_tiny_general_and_near_pi() {
        let cases = [
            Twist::zero(),
            Twist::new(
                Vec3::new(1.0e-12, -2.0e-12, 3.0e-12),
                Vec3::new(0.2, -0.1, 0.4),
            ),
            Twist::new(Vec3::new(0.3, -0.4, 0.7), Vec3::new(1.2, -0.8, 0.5)),
            Twist::new(Vec3::new(0.0, PI - 1.0e-9, 0.0), Vec3::new(-0.2, 0.4, 0.7)),
        ];
        for tangent in cases {
            let recovered = Se3::exp(tangent).unwrap().log();
            assert!(
                twist_error(recovered, tangent) < 3.0e-9,
                "SE(3) round trip failed for {tangent:?}: {recovered:?}"
            );
        }
    }

    #[test]
    fn g0_so3_jacobians_have_analytic_inverse() {
        for angular in [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0e-9, -2.0e-9, 3.0e-9),
            Vec3::new(0.4, -0.7, 1.2),
            Vec3::new(PI - 1.0e-8, 0.0, 0.0),
        ] {
            let tangent = So3Tangent::new(angular);
            let left = so3_left_jacobian(tangent).unwrap();
            let left_inverse = so3_left_jacobian_inverse(tangent).unwrap();
            let right = so3_right_jacobian(tangent).unwrap();
            let right_inverse = so3_right_jacobian_inverse(tangent).unwrap();
            assert!(matrix3_identity_error(left.compose(left_inverse)) < 2.0e-11);
            assert!(matrix3_identity_error(right.compose(right_inverse)) < 2.0e-11);
        }
    }

    #[test]
    fn g3_so3_and_se3_differentials_match_finite_difference() {
        let epsilon = 2.0e-7;
        let phi = So3Tangent::new(Vec3::new(0.4, -0.2, 0.7));
        let direction = So3Tangent::new(Vec3::new(-0.3, 0.5, 0.2));
        let base = So3::exp(phi).unwrap();
        let displaced = So3::exp(So3Tangent::new(
            phi.angular + direction.angular.scale(epsilon),
        ))
        .unwrap();
        let measured = displaced
            .space_minus(base)
            .unwrap()
            .angular
            .scale(1.0 / epsilon);
        let predicted = so3_left_jacobian(phi).unwrap().apply(direction.angular);
        assert!(vec_error(measured, predicted) < 2.0e-7);
        let measured_body = displaced
            .body_minus(base)
            .unwrap()
            .angular
            .scale(1.0 / epsilon);
        let predicted_body = so3_right_jacobian(phi).unwrap().apply(direction.angular);
        assert!(vec_error(measured_body, predicted_body) < 2.0e-7);

        let xi = Twist::new(Vec3::new(0.3, -0.4, 0.2), Vec3::new(0.8, -0.1, 0.5));
        let eta = Twist::new(Vec3::new(-0.2, 0.1, 0.4), Vec3::new(0.3, 0.7, -0.6));
        let base = Se3::exp(xi).unwrap();
        let displaced = Se3::exp(xi.plus(eta.scale(epsilon))).unwrap();
        let measured = displaced.space_minus(base).unwrap().scale(1.0 / epsilon);
        let differential = Se3::left_jacobian(xi).unwrap();
        assert!(differential.tail_bound < 1.0e-12);
        assert!(twist_error(measured, differential.apply(eta)) < 5.0e-7);
        let measured_body = displaced.body_minus(base).unwrap().scale(1.0 / epsilon);
        let differential_body = Se3::right_jacobian(xi).unwrap();
        assert!(differential_body.tail_bound < 1.0e-12);
        assert!(twist_error(measured_body, differential_body.apply(eta)) < 5.0e-7);
    }

    #[test]
    fn g0_adjoint_homomorphism_and_bracket_covariance() {
        let first = Se3::exp(Twist::new(
            Vec3::new(0.2, -0.5, 0.1),
            Vec3::new(0.7, 0.2, -0.4),
        ))
        .unwrap();
        let second = Se3::exp(Twist::new(
            Vec3::new(-0.3, 0.1, 0.6),
            Vec3::new(-0.2, 0.8, 0.3),
        ))
        .unwrap();
        let composed = first.compose(second).unwrap();
        assert!(
            matrix6_error(
                &composed.adjoint(),
                &first.adjoint().compose(second.adjoint())
            ) < 2.0e-12
        );

        let x = Twist::new(Vec3::new(0.4, 0.2, -0.1), Vec3::new(-0.7, 0.3, 0.5));
        let y = Twist::new(Vec3::new(-0.2, 0.6, 0.3), Vec3::new(0.1, -0.8, 0.4));
        let lhs = x.bracket(y).transform_by(&first);
        let rhs = x.transform_by(&first).bracket(y.transform_by(&first));
        assert!(twist_error(lhs, rhs) < 2.0e-12);
    }

    #[test]
    fn g0_coadjoint_preserves_dual_pairing() {
        let pose = Se3::exp(Twist::new(
            Vec3::new(0.3, -0.4, 0.2),
            Vec3::new(1.1, -0.7, 0.6),
        ))
        .unwrap();
        let twist = Twist::new(Vec3::new(-0.2, 0.5, 0.7), Vec3::new(0.4, -0.1, 0.8));
        let wrench = Wrench::new(Vec3::new(0.6, -0.3, 0.2), Vec3::new(-0.4, 0.9, 0.1));
        let before = wrench.pairing(twist);
        let after = wrench
            .transform_by(&pose)
            .pairing(twist.transform_by(&pose));
        assert!(
            (before - after).abs() < 2.0e-12,
            "pairing {before} != {after}"
        );
    }

    #[test]
    fn g0_left_and_right_perturbations_obey_adjoint_relation() {
        let pose = Se3::exp(Twist::new(
            Vec3::new(0.2, -0.1, 0.4),
            Vec3::new(0.6, 0.3, -0.5),
        ))
        .unwrap();
        let body_delta = Twist::new(Vec3::new(0.03, -0.02, 0.01), Vec3::new(-0.04, 0.05, 0.02));
        let space_delta = body_delta.transform_by(&pose);
        let from_body = pose.body_plus(body_delta).unwrap();
        let from_space = pose.space_plus(space_delta).unwrap();
        assert!(twist_error(from_body.body_minus(from_space).unwrap(), Twist::zero()) < 2.0e-12);
    }

    #[test]
    fn g5_replay_is_bit_identical_and_invalid_inputs_refuse() {
        let tangent = Twist::new(Vec3::new(0.31, -0.27, 0.83), Vec3::new(1.2, -0.4, 0.7));
        let first = Se3::exp(tangent).unwrap();
        let second = Se3::exp(tangent).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.log(), second.log());
        assert_eq!(
            first,
            Se3::try_from_motor(Motor(first.as_motor().0.scale(-1.0))).unwrap()
        );
        assert_eq!(
            Se3::left_jacobian(tangent).unwrap(),
            Se3::left_jacobian(tangent).unwrap()
        );

        let positive_pi = So3::try_from_quat(Quat {
            w: 0.0,
            x: 1.0,
            y: 0.0,
            z: 0.0,
        })
        .unwrap();
        let negative_pi = So3::try_from_quat(Quat {
            w: 0.0,
            x: -1.0,
            y: 0.0,
            z: 0.0,
        })
        .unwrap();
        assert_eq!(positive_pi, negative_pi);

        let signed_zero = So3::try_from_quat(Quat {
            w: 1.0,
            x: -0.0,
            y: 0.0,
            z: -0.0,
        })
        .unwrap();
        for coordinate in [
            signed_zero.as_quat().w,
            signed_zero.as_quat().x,
            signed_zero.as_quat().y,
            signed_zero.as_quat().z,
        ] {
            assert_ne!(coordinate.to_bits(), (-0.0f64).to_bits());
        }

        let mut signed_zero_motor = Motor::identity();
        for coordinate in &mut signed_zero_motor.0.0 {
            if *coordinate == 0.0 {
                *coordinate = -0.0;
            }
        }
        let canonical_motor = Se3::try_from_motor(signed_zero_motor).unwrap();
        assert_eq!(canonical_motor, Se3::identity());
        assert!(
            canonical_motor
                .as_motor()
                .0
                .0
                .iter()
                .all(|coordinate| coordinate.to_bits() != (-0.0f64).to_bits())
        );

        assert!(matches!(
            So3::try_from_quat(Quat {
                w: 2.0,
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }),
            Err(GaError::NotUnit { .. })
        ));
        assert!(matches!(
            So3::try_from_quat(Quat {
                w: 0.0,
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }),
            Err(GaError::DegenerateNorm { .. })
        ));
        assert!(matches!(
            Se3::exp(Twist::new(
                Vec3::new(f64::NAN, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
            )),
            Err(GaError::NonFinite { .. })
        ));
        assert!(matches!(
            Se3::left_jacobian(Twist::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(100.0, 0.0, 0.0),
            )),
            Err(GaError::IllConditioned { .. })
        ));
        let mut wrong_grade = Motor::identity();
        wrong_grade.0.0[1] = f64::EPSILON;
        assert!(matches!(
            Se3::try_from_motor(wrong_grade),
            Err(GaError::InvalidRepresentation { .. })
        ));
        assert!(matches!(
            So3::try_from_motor(Motor::translator(f64::EPSILON, 0.0, 0.0)),
            Err(GaError::InvalidRepresentation { .. })
        ));
        assert!(matches!(
            So3::identity().rotate(Vec3::new(0.0, f64::INFINITY, 0.0)),
            Err(GaError::NonFinite { .. })
        ));
    }
}

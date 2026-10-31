//! One potential, acceleration and spatial derivative for moving consumers.
//!
//! A deforming body samples this field at its actual quadrature positions.
//! Its gravity force and Newton tangent then differentiate the same potential;
//! density and reference mass remain with that body's original material owner.
//! No finite differences, eigensolve or allocation occur in these evaluations.

use super::{CentralGravity, GravityError, UniformGravity, finite, same_frame};
use crate::Vec3;
use fs_ga::Mat3;

/// A conservative point sample in one declared inertial coordinate frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointGravitySample {
    /// Acceleration in world coordinates, m/s².
    pub acceleration_world_m_s2: Vec3,
    /// Specific potential, J/kg; acceleration is its negative spatial gradient.
    pub specific_potential_j_kg: f64,
    /// Spatial acceleration derivative, s⁻², in row-major world coordinates.
    /// Entry `(i,j)` is `d acceleration_i / d position_j`, the negative
    /// Hessian of the specific potential. This is not a body inertia tensor.
    pub acceleration_gradient_s2: Mat3,
}

impl PointGravitySample {
    /// Check finite arithmetic and symmetry of a conservative-field derivative.
    ///
    /// Symmetry is checked relative to the largest entry, with a 64-epsilon
    /// roundoff allowance. This cannot establish that a user-supplied force,
    /// potential and derivative are derivatives of one another; that remains
    /// the field implementation's semantic obligation and verification task.
    pub fn validate(self) -> Result<Self, GravityError> {
        if !self.acceleration_world_m_s2.is_finite()
            || !self.specific_potential_j_kg.is_finite()
            || !self
                .acceleration_gradient_s2
                .m
                .iter()
                .all(|x| x.is_finite())
        {
            return Err(GravityError::Invalid("nonfinite point gravity sample"));
        }
        let matrix = self.acceleration_gradient_s2.m;
        let scale = matrix.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        if scale > 0.0 {
            for (a, b) in [(1, 3), (2, 6), (5, 7)] {
                if (matrix[a] / scale - matrix[b] / scale).abs() > 64.0 * f64::EPSILON {
                    return Err(GravityError::Invalid(
                        "conservative gravity acceleration gradient must be symmetric",
                    ));
                }
            }
        }
        Ok(self)
    }
}

/// A stationary conservative point field with its analytic spatial derivative.
///
/// Implementations return `g = -grad(phi)` and `Dg = -Hessian(phi)` in the same
/// inertial axes. The derivative enables consistent gravity tangents for
/// deformable bodies and variational/sensitivity equations. It does not include
/// centrifugal/Coriolis forces, a moving source's power, or finite-body moments.
/// A field may admit acceleration while refusing an unrepresentable derivative.
pub trait PointGravityField {
    /// Evaluate all three quantities at the actual requested world position.
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError>;
}

impl PointGravityField for UniformGravity {
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError> {
        same_frame(frame, self.frame)?;
        let displacement = position_world_m.sub(self.zero_potential_at_world_m);
        let acceleration = self.gravity.acceleration_world();
        PointGravitySample {
            acceleration_world_m_s2: acceleration,
            specific_potential_j_kg: -acceleration.dot(displacement),
            acceleration_gradient_s2: Mat3 { m: [0.0; 9] },
        }
        .validate()
    }
}

impl PointGravityField for CentralGravity {
    fn sample_point(
        &self,
        frame: u64,
        position_world_m: Vec3,
    ) -> Result<PointGravitySample, GravityError> {
        // Reuse the original scaled radius and source/frame admission. Repeated
        // division avoids forming r³, whose intermediate can overflow even
        // when every requested physical value is representable.
        let (direction, radius, potential) = self.radial(frame, position_world_m)?;
        let acceleration = direction.scale(finite(potential / radius)?);
        let mu_over_r3 = finite((-potential / radius) / radius)?;
        let n = [direction.x, direction.y, direction.z];
        let mut gradient = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                let diagonal = if i == j { 1.0 } else { 0.0 };
                gradient[3 * i + j] = mu_over_r3 * ((3.0 * n[i]) * n[j] - diagonal);
            }
        }
        PointGravitySample {
            acceleration_world_m_s2: acceleration,
            specific_potential_j_kg: potential,
            acceleration_gradient_s2: Mat3 { m: gradient },
        }
        .validate()
    }
}

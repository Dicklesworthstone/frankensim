//! Fixed material support and frame preparation for the common native solve.

use super::{PorousMomentumError, invalid, norm};
use fs_mbd::{UnitQuaternion, Vec3};

/// A centred rigid region carrying homogeneous porous resistance.
///
/// Dimensions are metres in the original body's material coordinates, not
/// grid coordinates. Fluid remains represented inside this region; the shape
/// neither excludes volume nor creates displaced mass or additional inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PorousBodySupport {
    /// Rotation-invariant centred sphere.
    Sphere { radius_m: f64 },
    /// Centred ellipsoid, with semiaxes along the body's x/y/z axes.
    Ellipsoid { semiaxes_m: [f64; 3] },
}

impl PorousBodySupport {
    pub(super) fn admit(self) -> Result<f64, PorousMomentumError> {
        let radii = match self {
            Self::Sphere { radius_m } => [radius_m; 3],
            Self::Ellipsoid { semiaxes_m } => semiaxes_m,
        };
        if radii.iter().any(|r| !r.is_finite() || *r <= 0.0) {
            return Err(invalid(
                "porous support requires finite positive dimensions",
            ));
        }
        Ok(radii.into_iter().fold(0.0_f64, f64::max))
    }
}

/// Homogeneous resistance in the body's fixed material axes.
///
/// Lattice `dx = dt = 1` units match the native Darcy--Forchheimer owner:
/// force density is `-R diag(darcy_principal) R^T (u-v_s)
/// - forchheimer |u-v_s| (u-v_s)`, with the actual rigid material frame R.
/// These coefficients are not accelerations or intrinsic permeabilities.
/// The linear axes may be unrelated to the body's full inertia eigenvectors.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PrincipalPorousResistance {
    pub darcy_principal: [f64; 3],
    pub forchheimer: f64,
}

impl PrincipalPorousResistance {
    /// Isotropic specialization, retaining the original native scalar law.
    pub const fn isotropic(darcy: f64, forchheimer: f64) -> Self {
        Self {
            darcy_principal: [darcy; 3],
            forchheimer,
        }
    }

    pub(super) fn admit(self) -> Result<(), PorousMomentumError> {
        if self
            .darcy_principal
            .iter()
            .chain(core::iter::once(&self.forchheimer))
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "principal porous resistance requires finite nonnegative coefficients",
            ));
        }
        Ok(())
    }
}

/// Prepared once per native trial, not per population or grid cell.
pub(super) struct PreparedSupport {
    shape: PorousBodySupport,
    /// Columns of the original proper material-to-world rotation.
    axes: [[f64; 3]; 3],
}

impl PreparedSupport {
    pub(super) fn new(shape: PorousBodySupport, orientation: UnitQuaternion) -> Self {
        Self {
            shape,
            axes: [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ]
            .map(|axis| {
                let v = orientation.rotate_body_to_world(axis);
                [v.x, v.y, v.z]
            }),
        }
    }

    /// Actual cell-centre support in this trial's midpoint material frame.
    pub(super) fn contains(&self, arm_world_m: Vec3) -> bool {
        match self.shape {
            PorousBodySupport::Sphere { radius_m } => norm(arm_world_m) <= radius_m,
            PorousBodySupport::Ellipsoid { semiaxes_m } => {
                let mut distance = 0.0_f64;
                for (axis, radius) in self.axes.iter().zip(semiaxes_m) {
                    let coordinate = axis[0].mul_add(
                        arm_world_m.x,
                        axis[1].mul_add(arm_world_m.y, axis[2] * arm_world_m.z),
                    );
                    // A box rejection avoids overflowing ratios far outside a
                    // thin region. Squared normalized coordinates are <= one.
                    if !coordinate.is_finite() || coordinate.abs() > radius {
                        return false;
                    }
                    let normalized = coordinate / radius;
                    distance = normalized.mul_add(normalized, distance);
                    if distance > 1.0 {
                        return false;
                    }
                }
                distance <= 1.0
            }
        }
    }

    pub(super) fn resistance_axes(&self, resistance: PrincipalPorousResistance) -> [[f64; 3]; 3] {
        // The original scalar lane keeps its identical native representation;
        // rotating an isotropic tensor would add arithmetic without physics.
        if resistance
            .darcy_principal
            .iter()
            .all(|&v| v == resistance.darcy_principal[0])
        {
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        } else {
            self.axes
        }
    }
}

/// Geodesic midpoint of the actual native attitudes, with quaternion sign
/// continuity across the canonical double-cover cut. A Cayley rotation is
/// shorter than pi; its normalized half-sum is the unique shortest midpoint.
pub(super) fn midpoint_orientation(
    before: UnitQuaternion,
    after: UnitQuaternion,
) -> Result<UnitQuaternion, PorousMomentumError> {
    if before == after {
        return Ok(before);
    }
    let a = before.components();
    let b = after.components();
    let product = a[0].mul_add(b[0], a[1].mul_add(b[1], a[2].mul_add(b[2], a[3] * b[3])));
    let sign = if product < 0.0 { -1.0 } else { 1.0 };
    Ok(UnitQuaternion::new(
        a[0] + sign * b[0],
        a[1] + sign * b[1],
        a[2] + sign * b[2],
        a[3] + sign * b[3],
    )?)
}

#[cfg(test)]
#[path = "support/tests.rs"]
mod tests;

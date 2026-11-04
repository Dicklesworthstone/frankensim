//! Authored anisotropic GGX slopes in a physical tangent frame.
//!
//! Reflection and transmission retain the original Walter half-vector,
//! Fresnel, masking-product and radiance Jacobians. Only D, directional G1
//! and the two-axis stretch in the original visible-normal warp differ.
//! See Walter et al., EGSR 2007, <https://doi.org/10.2312/EGWR/EGSR07/195-206>
//! and Heitz, JCGT 2018, <https://jcgt.org/published/0007/04/01/>.
//!
//! This is a single-scattering slope model. It does not infer grooves from
//! height roughness, add inter-facet multiple scattering, or deform geometry.

use super::{DielectricError, PI, scale_safe_norm, validate_roughness};
use fs_geom::Vec3;

/// Two admitted GGX slope scales and an explicit axis in the query frame.
/// The scalar scene tracer stores this axis in geometry-local coordinates.
/// At a contact its projection onto the tangent plane defines alpha-x;
/// alpha-y lies along `normal cross tangent`. A parallel axis refuses rather
/// than silently choosing a different finish at that contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnisotropicDielectricRoughness {
    alpha_x: f64,
    alpha_y: f64,
    tangent_axis: Vec3,
}

impl AnisotropicDielectricRoughness {
    pub fn try_new(
        alpha_x: f64,
        alpha_y: f64,
        tangent_axis: Vec3,
    ) -> Result<Self, DielectricError> {
        validate_roughness(alpha_x)?;
        validate_roughness(alpha_y)?;
        Ok(Self {
            alpha_x,
            alpha_y,
            tangent_axis: unit(tangent_axis)?,
        })
    }

    /// GGX slope scales `[alpha_x, alpha_y]`, with no perceptual remapping.
    pub const fn alphas(self) -> [f64; 2] {
        [self.alpha_x, self.alpha_y]
    }

    /// Original normalized authored axis, before contact projection.
    pub const fn tangent_axis(self) -> Vec3 {
        self.tangent_axis
    }

    /// Project an authored axis without the subtractive cancellation in
    /// `axis - normal*dot(axis, normal)`. Reversing the normal retains the same
    /// tangent and reverses only the bitangent, preserving the physical axes.
    pub(crate) fn tangent_frame(self, normal: Vec3) -> Result<(Vec3, Vec3), DielectricError> {
        let normal = unit(normal)?;
        let bitangent = unit(cross(normal, self.tangent_axis))?;
        let tangent = unit(cross(bitangent, normal))?;
        Ok((tangent, bitangent))
    }

    pub(super) fn prepare(self, normal: Vec3, alpha_x: f64) -> Result<Frame, DielectricError> {
        // The legacy interface accepts alpha-x separately. It cannot silently
        // override the complete law retained by an anisotropic surface.
        if alpha_x != self.alpha_x {
            return Err(DielectricError::InvalidRoughness);
        }
        let (tangent, bitangent) = self.tangent_frame(normal)?;
        Ok(Frame {
            normal,
            tangent,
            bitangent,
            law: self,
        })
    }
}

#[derive(Clone, Copy)]
pub(super) struct Frame {
    normal: Vec3,
    tangent: Vec3,
    bitangent: Vec3,
    law: AnisotropicDielectricRoughness,
}

impl Frame {
    pub(super) fn distribution(self, micro_normal: Vec3) -> f64 {
        let z = self.normal.dot(micro_normal);
        if z <= 0.0 {
            return 0.0;
        }
        let x = self.tangent.dot(micro_normal) / self.law.alpha_x;
        let y = self.bitangent.dot(micro_normal) / self.law.alpha_y;
        let denominator = x * x + y * y + z * z;
        1.0 / (PI * self.law.alpha_x * self.law.alpha_y * denominator * denominator)
    }

    pub(super) fn masking(self, direction: Vec3) -> f64 {
        let z = self.normal.dot(direction).abs().min(1.0);
        if z == 0.0 {
            return 0.0;
        }
        let x = self.law.alpha_x * self.tangent.dot(direction);
        let y = self.law.alpha_y * self.bitangent.dot(direction);
        2.0 * z / (z + (x * x + y * y + z * z).sqrt())
    }

    pub(super) fn sample(self, wo: Vec3, u: f64, v: f64) -> Result<Vec3, DielectricError> {
        super::sample_ggx_visible_normal_in_frame(
            self.normal,
            self.tangent,
            self.bitangent,
            wo,
            self.law.alpha_x,
            self.law.alpha_y,
            u,
            v,
        )
    }
}

fn unit(v: Vec3) -> Result<Vec3, DielectricError> {
    if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) {
        return Err(DielectricError::InvalidDirection);
    }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 {
        return Err(DielectricError::InvalidDirection);
    }
    let scaled = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
    let length = scale_safe_norm(scaled);
    if !length.is_finite() || length == 0.0 {
        return Err(DielectricError::InvalidDirection);
    }
    let normalized = scaled.scale(1.0 / length);
    Ok(Vec3::new(
        if normalized.x == 0.0 {
            0.0
        } else {
            normalized.x
        },
        if normalized.y == 0.0 {
            0.0
        } else {
            normalized.y
        },
        if normalized.z == 0.0 {
            0.0
        } else {
            normalized.z
        },
    ))
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

#[cfg(test)]
mod tests;

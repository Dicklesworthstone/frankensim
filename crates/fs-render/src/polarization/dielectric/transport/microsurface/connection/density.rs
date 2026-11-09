//! Deterministic density of the single-reflection subset of a Smith walk.
//!
//! This is the analytic height-integrated, unit-reflectance BRDF times the
//! outgoing cosine. Its integral is the probability of escaping after ONE
//! reflection; higher orders are missing. It is not a normalized full-walk
//! directional PDF. No stochastic response, optical evaluation, history or
//! guessed multiple-scattering albedo enters these constant-work queries.

use super::{SmithBeamError, local};
use crate::conductor::ConductorSurface;
use crate::conductor::microsurface::{SmithGgx, unit};
use crate::instances::RigidTransform;
use crate::polarization::dielectric::transport::connection::checked_product;
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, PolarizedBeam};
use crate::surface_detail::SurfaceDetail;
use fs_geom::{Point3, Vec3};

impl PolarizedBeam {
    /// One-reflection escape density in solid angle, with unit facet albedo.
    /// The outgoing direction points away from the surface. Macro-horizon
    /// clipping and the original normal-facing convention are retained.
    pub fn smith_conductor_first_reflection_density(
        self,
        normal: Vec3,
        surface: ConductorSurface,
        outgoing: Vec3,
    ) -> Result<f64, SmithBeamError> {
        let outgoing = unit(outgoing)?;
        let mut normal = unit(normal)?;
        if normal.dot(self.frame().direction()) > 0.0 {
            normal = normal.scale(-1.0);
        }
        if normal.dot(self.frame().direction()) >= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        if surface.is_delta() || normal.dot(outgoing) <= 0.0 {
            return Ok(0.0);
        }
        let reference = if normal.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let frame = PolarizationFrame::try_new(normal, reference)?;
        let alpha = surface.roughness_alpha();
        first_density(
            SmithGgx::try_new(alpha, alpha)?,
            [frame.horizontal(), frame.vertical(), normal],
            self.frame().direction(),
            outgoing,
        )
    }

    /// The same subdensity in the native body-locked brushed frame, including
    /// its anisotropy, authored mark and exact undefined-axis fallback.
    /// No radiance/cosine adjoint multiplier belongs in a sampling density:
    /// changing the local frame is a rotation of solid angle.
    pub fn instanced_smith_brushed_first_reflection_density(
        self,
        detail: &SurfaceDetail,
        transform: RigidTransform,
        local_point: Point3,
        world_normal: Vec3,
        outgoing: Vec3,
    ) -> Result<f64, SmithBeamError> {
        let outgoing = unit(outgoing)?;
        let Some(frame) =
            self.smith_instanced_brushed_frame(detail, transform, local_point, world_normal)?
        else {
            let surface = ConductorSurface::try_rough(detail.alpha_pair().0)
                .map_err(super::super::ConductorBeamError::from)?;
            return self.smith_conductor_first_reflection_density(world_normal, surface, outgoing);
        };
        if self
            .smith_brushed_endpoint_weight(&frame, outgoing, true)?
            .is_none()
        {
            return Ok(0.0);
        }
        let (ax, ay) = detail.alpha_pair_for_mark(frame.mark_covered);
        first_density(
            SmithGgx::try_new(ax, ay)?,
            [frame.tangent, frame.bitangent, frame.shading_normal],
            self.frame().direction(),
            outgoing,
        )
    }
}

fn first_density(
    model: SmithGgx,
    basis: [Vec3; 3],
    incoming: Vec3,
    outgoing: Vec3,
) -> Result<f64, SmithBeamError> {
    let outgoing = local(basis, outgoing);
    let Some(connection) = model.single_scattering_connection(local(basis, incoming), outgoing)?
    else {
        return Ok(0.0);
    };
    Ok(checked_product(connection.factor_per_sr, outgoing.z)?)
}

#[cfg(test)]
mod tests;

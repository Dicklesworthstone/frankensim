//! Evaluated rough metal connections using the original scalar GGX kernels.

use super::polarized_brushed_conductor::{admit_frame, instanced_frame, unit};
use crate::conductor::{ConductorOptics, ConductorSurface};
use crate::dielectric::DielectricEvent;
use crate::instances::RigidTransform;
use crate::polarization::Stokes;
use crate::polarization::dielectric::mirrors::ConductorBeamError;
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, PolarizedBeam};
use crate::polarization::dielectric::transport::connection::{
    SurfaceConnectionResponse, checked_product, normalized_camera, reciprocal_contraction,
};
use crate::surface_detail::SurfaceDetail;
use crate::tracer::{
    brushed_local, ggx_vndf_reflection_pdf, ggx_vndf_reflection_pdf_aniso,
    smith_g1, smith_g1_aniso,
};
use fs_geom::{Point3, Vec3};

impl PolarizedBeam {
    /// Evaluate continuous rough metal reflection and its unchanged VNDF PDF.
    ///
    /// Retains the original spectral complex IOR and authored coherent optical
    /// stack, including circular/elliptical polarization and the actual local
    /// phase index. `source_frame` travels into the vertex and its Stokes
    /// intensity is one. No sampled facet is reused or resampled. A direction
    /// below the geometric horizon and a smooth delta surface both have zero
    /// continuous response. See [`SurfaceConnectionResponse`] for radiometry.
    pub fn evaluate_rough_conductor_connection(
        self,
        normal: Vec3,
        optics: &ConductorOptics,
        surface: ConductorSurface,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, ConductorBeamError> {
        let camera = normalized_camera(self, normalized_source)?;
        let wo = camera.frame().direction().scale(-1.0);
        let wi = source_frame.direction().scale(-1.0);
        let mut normal = unit(normal)?;
        if normal.dot(wo) < 0.0 { normal = normal.scale(-1.0); }
        if normal.dot(wo) <= 0.0 { return Err(BeamError::InvalidFrame.into()); }
        if surface.is_delta() || normal.dot(wi) <= 0.0 {
            return Ok(SurfaceConnectionResponse::zero(DielectricEvent::Reflection));
        }
        let facet = unit(Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z))?;
        let alpha = surface.roughness_alpha();
        let pdf = ggx_vndf_reflection_pdf(alpha, normal, wo, facet);
        let mask = smith_g1(alpha, normal.dot(wi).clamp(0.0, 1.0));
        let reflected = camera.reflect_conductor(facet, optics)?.beam;
        let response = reciprocal_contraction(reflected, source_frame, normalized_source)?;
        // Native f = pdf * G1(wi) * F / geometric_cos_i. Retain this
        // factorization, avoiding division by a possibly tiny reflectance.
        let value = checked_product(pdf, mask)? / normal.dot(wi);
        let factor = checked_product(value, response)?;
        Ok(SurfaceConnectionResponse::checked(factor, pdf, DielectricEvent::Reflection)?)
    }

    /// Evaluate the original body-locked anisotropic metal for camera radiance.
    ///
    /// The physical receiving cosine belongs to the GEOMETRIC normal, while
    /// native GGX masking and proposal use the authored shading frame. The
    /// radiance law is D*G*F/(4 ns.wo ng.wi); the forward-power sampler's
    /// adjoint shading-normal correction must not be applied to this response.
    /// Zero detail and the undefined radial pole reuse isotropic evaluation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_instanced_brushed_conductor_connection(
        self,
        optics: &ConductorOptics,
        detail: &SurfaceDetail,
        transform: RigidTransform,
        local_point: Point3,
        world_normal: Vec3,
        source_frame: PolarizationFrame,
        normalized_source: Stokes,
    ) -> Result<SurfaceConnectionResponse, ConductorBeamError> {
        let camera = normalized_camera(self, normalized_source)?;
        let Some(frame) = instanced_frame(camera, detail, transform, local_point, world_normal)? else {
            return camera.evaluate_rough_conductor_connection(
                world_normal, optics, ConductorSurface::try_rough(detail.alpha_pair().0)?,
                source_frame, normalized_source,
            );
        };
        admit_frame(&frame)?;
        let wo = camera.frame().direction().scale(-1.0);
        let wi = source_frame.direction().scale(-1.0);
        let geometric_out = frame.geometric_normal.dot(wo);
        if geometric_out <= 0.0 { return Err(BeamError::InvalidFrame.into()); }
        let geometric_in = frame.geometric_normal.dot(wi);
        let wo_local = brushed_local(&frame, wo);
        let wi_local = brushed_local(&frame, wi);
        if geometric_in <= 0.0 || wo_local[2] <= 0.0 || wi_local[2] <= 0.0 {
            return Ok(SurfaceConnectionResponse::zero(DielectricEvent::Reflection));
        }
        let facet = unit(Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z))?;
        let (alpha_x, alpha_y) = detail.alpha_pair_for_mark(frame.mark_covered);
        let pdf = ggx_vndf_reflection_pdf_aniso(
            alpha_x, alpha_y, wo_local, brushed_local(&frame, facet),
        );
        let mask = smith_g1_aniso(alpha_x, alpha_y, wi_local);
        let reflected = camera.reflect_conductor(facet, optics)?.beam;
        let response = reciprocal_contraction(reflected, source_frame, normalized_source)?;
        let value = checked_product(pdf, mask)? / geometric_in;
        let factor = checked_product(value, response)?;
        Ok(SurfaceConnectionResponse::checked(factor, pdf, DielectricEvent::Reflection)?)
    }
}

#[cfg(test)]
mod tests;

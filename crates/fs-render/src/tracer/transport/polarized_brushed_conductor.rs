//! Polarized importance transport over the native body-locked brushed finish.
//!
//! The sampler, anisotropic distribution and procedural surface frame retain
//! their scalar owners. A tilted shading normal needs the *adjoint* of the
//! effective camera BSDF, not its radiance weight. See PBRT 3e, section 16.1.3:
//! https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation
//! This is the same shading approximation, not resolved groove geometry or an
//! energy-conservation guarantee for an arbitrarily normal-mapped surface.

use fs_geom::{Point3, Vec3};
use fs_math::det;

use crate::conductor::{ConductorError, ConductorOptics, ConductorSurface};
use crate::dielectric::DielectricEvent;
use crate::instances::RigidTransform;
use crate::polarization::Stokes;
use crate::polarization::dielectric::mirrors::ConductorBeamError;
use crate::polarization::dielectric::transport::{
    BeamError, BeamSample, PolarizationFrame, PolarizedBeam,
};
use crate::surface_detail::{ShadingFrame, SurfaceDetail, SurfaceDetailError, shading_frame};
use crate::tracer::{
    brushed_local, brushed_world, sample_ggx_vndf_reflection_aniso, smith_g1_aniso,
};

impl PolarizedBeam {
    /// Scatter in an explicitly supplied world-space brushed frame.
    ///
    /// The complete orthonormal frame must face the incoming ray according to
    /// its GEOMETRIC normal. A ray below the shading horizon, or a native
    /// reflected sample below either horizon, contributes zero (`Ok(None)`).
    /// Invalid frames, alpha values and uniforms are errors, not zero samples.
    /// No rejection is resampled or removed from the proposal normalization.
    ///
    /// Forward power uses `F_mueller * G1(wi) * C`, where
    /// `C = (ns.wo / ng.wo) * (ng.wi / ns.wi)`. This is the adjoint of the
    /// existing effective radiance BSDF. The small `ns.wi` cancels analytically
    /// against masking before division. Never clamp C to one: a valid weighted
    /// sample can exceed input power. This empirical shading-normal model is
    /// not a physical heat-deposition law or an inter-facet multiple-scattering
    /// model. With untilted normals C is one and the response is passive.
    ///
    /// Bare and coherent-coated optics, phase indices, polarization frames and
    /// reflection handedness use the existing conductor/beam owners. The PDF
    /// remains the native solid-angle density, independent of polarization.
    pub fn sample_brushed_conductor(
        self,
        optics: &ConductorOptics,
        detail: &SurfaceDetail,
        frame: &ShadingFrame,
        samples: [f64; 2],
    ) -> Result<Option<BeamSample>, ConductorBeamError> {
        admit_samples(samples)?;
        admit_detail(detail)?;
        admit_frame(frame)?;
        let wo = self.frame().direction().scale(-1.0);
        let geometric_out = frame.geometric_normal.dot(wo);
        if geometric_out <= 0.0 {
            return Err(BeamError::InvalidFrame.into());
        }
        let wo_local = brushed_local(frame, wo);
        if wo_local[2] <= 0.0 {
            return Ok(None);
        }
        let (alpha_x, alpha_y) = detail.alpha_pair_for_mark(frame.mark_covered);
        let Some((wi_local, pdf)) = sample_ggx_vndf_reflection_aniso(
            wo_local, alpha_x, alpha_y, samples[0], samples[1],
        ) else {
            return Ok(None);
        };
        let direction = brushed_world(frame, wi_local);
        let geometric_in = frame.geometric_normal.dot(direction);
        if geometric_in <= 0.0 {
            return Ok(None);
        }
        if !pdf.is_finite() || pdf <= 0.0 {
            return Err(BeamError::NumericalRange.into());
        }
        let facet = unit(Vec3::new(
            wo.x + direction.x, wo.y + direction.y, wo.z + direction.z,
        ))?;
        let reflected = self.reflect_conductor(facet, optics)?.beam;
        let weight = if frame.geometric_normal == frame.shading_normal {
            smith_g1_aniso(alpha_x, alpha_y, wi_local)
        } else {
            // G1(wi)/ns.wi = 2/(ns.wi + |(ax*wi.x,ay*wi.y,wi.z)|).
            // This avoids the removable shading-horizon singularity and the
            // huge D/pdf intermediates of a narrow lobe. The ng.wo horizon is
            // a real singularity of the authored shading approximation.
            let length = det::hypot(
                det::hypot(alpha_x * wi_local[0], alpha_y * wi_local[1]), wi_local[2],
            );
            (2.0 * geometric_in / (wi_local[2] + length)) * (wo_local[2] / geometric_out)
        };
        if !weight.is_finite() || weight < 0.0 {
            return Err(BeamError::NumericalRange.into());
        }
        let s = reflected.power();
        let power = Stokes { i: s.i * weight, q: s.q * weight, u: s.u * weight, v: s.v * weight };
        if ![power.i, power.q, power.u, power.v].iter().all(|v| v.is_finite()) {
            return Err(BeamError::NumericalRange.into());
        }
        let output_frame = PolarizationFrame::try_new(direction, reflected.frame().horizontal())?;
        let beam = Self::try_new(
            self.ray().origin, output_frame, self.wavelength_nm(), self.medium(), power,
        )?;
        Ok(Some(BeamSample { beam, event: DielectricEvent::Reflection, pdf, delta: false }))
    }

    /// Derive native body-locked brushing from an instance's local surface point.
    ///
    /// `transform` is the original proper-rigid body-to-world placement. The
    /// world geometric normal is normalized and face-forwarded BEFORE deriving
    /// the local shading frame, matching native scalar material evaluation.
    /// Local radial/circular grooves and the roughness mark therefore move with
    /// the object, not with the ray or polarization axis. This does not move the
    /// physical beam origin or change its surrounding medium.
    ///
    /// Exact zero detail and the chart's undefined radial pole reuse the
    /// original isotropic sampler (including its samples and PDF), not an
    /// invented tangent. A curved contact may supply its local geometric
    /// representative without snapping the integrated beam to that point.
    pub fn sample_instanced_brushed_conductor(
        self,
        optics: &ConductorOptics,
        detail: &SurfaceDetail,
        transform: RigidTransform,
        local_point: Point3,
        world_normal: Vec3,
        samples: [f64; 2],
    ) -> Result<Option<BeamSample>, ConductorBeamError> {
        admit_samples(samples)?;
        admit_detail(detail)?;
        if ![local_point.x, local_point.y, local_point.z].iter().all(|v| v.is_finite()) {
            return Err(BeamError::InvalidPosition.into());
        }
        let mut normal = unit(world_normal)?;
        if normal.dot(self.ray().dir) > 0.0 { normal = normal.scale(-1.0); }
        let fallback = || self.sample_rough_conductor(
            world_normal, optics, ConductorSurface::try_rough(detail.alpha_pair().0)?, samples,
        );
        if detail.is_zero_detail() { return fallback(); }
        let inverse = transform.inverse().map_err(|_| BeamError::InvalidFrame)?;
        let local_frame = match shading_frame(detail, local_point, inverse.transform_vector(normal)) {
            Ok(frame) => frame,
            Err(SurfaceDetailError::PoleTangentUndefined) => return fallback(),
            Err(_) => return Err(BeamError::InvalidFrame.into()),
        };
        let frame = ShadingFrame {
            geometric_normal: normal,
            shading_normal: transform.transform_vector(local_frame.shading_normal),
            tangent: transform.transform_vector(local_frame.tangent),
            bitangent: transform.transform_vector(local_frame.bitangent),
            mark_covered: local_frame.mark_covered,
        };
        self.sample_brushed_conductor(optics, detail, &frame, samples)
    }
}

fn admit_samples(samples: [f64; 2]) -> Result<(), ConductorBeamError> {
    if samples.iter().any(|v| !v.is_finite() || !(0.0..1.0).contains(v)) {
        return Err(ConductorError::InvalidSample.into());
    }
    Ok(())
}

fn admit_detail(detail: &SurfaceDetail) -> Result<(), ConductorBeamError> {
    // SurfaceDetail::none is const and intentionally does not validate alpha.
    // Validate the original pair before mark clamping could hide invalid data.
    let (x, y) = detail.alpha_pair();
    ConductorSurface::try_rough(x)?;
    ConductorSurface::try_rough(y)?;
    Ok(())
}

fn admit_frame(frame: &ShadingFrame) -> Result<(), BeamError> {
    let n = frame.geometric_normal;
    let s = frame.shading_normal;
    let t = frame.tangent;
    let b = frame.bitangent;
    const TOL: f64 = 2e-12;
    for v in [n, s, t, b] {
        if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) || (v.dot(v) - 1.0).abs() > TOL {
            return Err(BeamError::InvalidFrame);
        }
    }
    let cross = Vec3::new(t.y*b.z-t.z*b.y, t.z*b.x-t.x*b.z, t.x*b.y-t.y*b.x);
    if t.dot(b).abs() > TOL || s.dot(t).abs() > TOL || s.dot(b).abs() > TOL
        || n.dot(s) <= 0.0 || (cross.dot(s) - 1.0).abs() > TOL {
        return Err(BeamError::InvalidFrame);
    }
    Ok(())
}

fn unit(v: Vec3) -> Result<Vec3, BeamError> {
    if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) { return Err(BeamError::InvalidFrame); }
    let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
    if scale == 0.0 { return Err(BeamError::InvalidFrame); }
    let v = Vec3::new(v.x/scale, v.y/scale, v.z/scale);
    Ok(v.scale(1.0/det::sqrt(v.dot(v))))
}

#[cfg(test)]
mod tests;

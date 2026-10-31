//! Source-depth-conditioned incoherent imaging over admitted separated planes.
//!
//! Sample source direction independently of the eventual hit depth, then use
//! that source's PSF. Selecting a PSF from the unshifted destination hit would
//! apply the wrong optical kernel across a depth boundary. The renderer also
//! establishes full-pupil visibility; a first-hit depth check cannot do that.

use super::{
    AnnularPupil, ApertureKind, CameraError, CircularDefocus, LensSample, PhysicalCamera,
    PolygonPupil, SpectralCameraRay, annular, direct_unit,
};
use crate::charts::Ray;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

#[cfg(any(feature = "tracer", test))]
mod visibility;
#[cfg(feature = "tracer")]
pub(crate) use visibility::ProjectedPupilFootprint;

const PI: f64 = core::f64::consts::PI;
const MAX_LAYERED_PLANES: usize = 16;

/// Scalar incoherent imaging of several explicit, pupil-visible source planes.
///
/// Each plane is frontoparallel at its positive axial distance from the pupil.
/// The renderer admits separated static opaque diffuse meshes with black misses,
/// proving that no declared plane clips a ray from another plane to any pupil
/// point throughout the owning camera shot. Camera-only callers must establish
/// equivalent visibility and pupil-independent source radiance. This does not
/// model arbitrary occluding layers, coherent scene transport or inter-object
/// edge-wave scattering.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredDiffraction {
    depths_m: [f64; MAX_LAYERED_PLANES],
    count: usize,
    maximum_depth_phase_error_rad: f64,
}

impl LayeredDiffraction {
    /// Fixed bound on source planes and the pairwise pupil-clearance admission.
    pub const MAX_PLANES: usize = MAX_LAYERED_PLANES;

    /// Declare one to sixteen strictly increasing finite positive axial depths
    /// and a residual pupil-rim phase budget in `[0,0.1]` radians. Geometry must
    /// separately be coplanar to coordinate roundoff; this phase budget cannot
    /// admit an occluder at a nearby, undeclared depth.
    pub fn try_new(
        object_distances_m: &[f64],
        maximum_depth_phase_error_rad: f64,
    ) -> Result<Self, CameraError> {
        if object_distances_m.is_empty()
            || object_distances_m.len() > Self::MAX_PLANES
            || object_distances_m
                .iter()
                .any(|depth| !depth.is_finite() || *depth <= 0.0 || !(1.0 / depth).is_finite())
            || object_distances_m.windows(2).any(|pair| pair[0] >= pair[1])
            || !maximum_depth_phase_error_rad.is_finite()
            || !(0.0..=0.1).contains(&maximum_depth_phase_error_rad)
        {
            return Err(CameraError::InvalidDiffraction);
        }
        let mut depths_m = [0.0; Self::MAX_PLANES];
        depths_m[..object_distances_m.len()].copy_from_slice(object_distances_m);
        Ok(Self {
            depths_m,
            count: object_distances_m.len(),
            maximum_depth_phase_error_rad,
        })
    }

    /// Declared axial depths [m], ordered from nearest to farthest.
    #[must_use]
    pub fn object_distances_m(&self) -> &[f64] {
        &self.depths_m[..self.count]
    }

    /// Maximum phase discrepancy between an actual source and its declared plane.
    #[must_use]
    pub const fn maximum_depth_phase_error_rad(self) -> f64 {
        self.maximum_depth_phase_error_rad
    }
}

/// A sampled source ray whose spectral weights require its actual first hit.
///
/// This token borrows the exact camera that generated it, retaining the shared
/// proposal and wavelength packet until [`Self::resolve_source`] consumes it.
/// No provisional weights are exposed. Query the scene with [`Self::ray`], then
/// resolve that first hit, or `None` for the admitted black background.
#[must_use = "resolve the ray at its actual first visible source before accumulating radiance"]
pub struct LayeredCameraRay<'a, const N: usize> {
    camera: &'a PhysicalCamera,
    ray: Ray,
    q: [f64; 2],
    proposal: f64,
    reference_nm: f64,
    wavelengths_nm: [f64; N],
}

impl<const N: usize> LayeredCameraRay<'_, N> {
    /// Optical-centre ray to intersect against the admitted source geometry.
    #[must_use]
    pub const fn ray(&self) -> Ray {
        self.ray
    }

    /// Resolve the actual first source of this sampled ray. A black miss has
    /// zero weights. Finite hits retain their own depth and each wavelength's
    /// signed phase, rather than the depth at the original raster location.
    pub fn resolve_source(
        self,
        cx: &Cx<'_>,
        point: Option<Point3>,
    ) -> Result<SpectralCameraRay<N>, CameraError> {
        cx.checkpoint()?;
        let Some(point) = point else {
            return Ok(SpectralCameraRay {
                ray: self.ray,
                weights: [0.0; N],
            });
        };
        self.camera
            .admit_diffraction_depth(Some(point), &self.wavelengths_nm)?;
        let depth = point.delta_from(self.camera.eye).dot(self.camera.forward);
        let mut weights = [0.0; N];
        for (weight, &wavelength) in weights.iter_mut().zip(&self.wavelengths_nm) {
            cx.checkpoint()?;
            let ratio = self.reference_nm / wavelength;
            let phase = self.camera.pupil_defocus_phase(1.0 / depth, wavelength)?;
            *weight = ratio
                * ratio
                * self
                    .camera
                    .layered_pupil_density(self.q.map(|value| ratio * value), phase)?
                / self.proposal;
            if !weight.is_finite() || *weight < 0.0 {
                return Err(CameraError::InvalidDiffraction);
            }
        }
        Ok(SpectralCameraRay {
            ray: self.ray,
            weights,
        })
    }
}

impl PhysicalCamera {
    /// Select actual-source-depth diffraction over explicit separated planes.
    /// This replaces the focused or single-plane policy. Use the two-stage
    /// [`Self::sample_layered_diffraction_ray`] API; the older one-stage spectral
    /// method refuses because the first source depth is not yet known there.
    pub fn with_layered_diffraction(
        mut self,
        policy: LayeredDiffraction,
    ) -> Result<Self, CameraError> {
        if self.aperture.is_pinhole() {
            return Err(CameraError::InvalidDiffraction);
        }
        self.admit_diffraction_projection()?;
        self.layered_diffraction = Some(policy);
        self.diffraction = None;
        self.defocused_diffraction = None;
        Ok(self)
    }

    /// Selected source-depth-resolved imaging policy, if any.
    #[must_use]
    pub const fn layered_diffraction(&self) -> Option<LayeredDiffraction> {
        self.layered_diffraction
    }

    /// Sample a full-support direction independently of its unknown source depth.
    /// The proposal covers the widest declared defocus plus a narrow component;
    /// all admitted source depths share the same complete mixture density.
    /// Pixel, wavelength and shutter draws remain independent of this sample.
    pub fn sample_layered_diffraction_ray<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        diffraction_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<LayeredCameraRay<'_, N>, CameraError> {
        cx.checkpoint()?;
        if self.layered_diffraction.is_none() || !x_tan.is_finite() || !y_tan.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        let (minimum, maximum) = self.diffraction_wavelength_range(wavelengths_nm)?;
        let reference_nm = 0.5 * minimum + 0.5 * maximum;
        let phase = self.layered_maximum_phase(reference_nm)?;
        let (q, proposal) = match &self.aperture.kind {
            ApertureKind::RegularPolygon { vertices, .. } => {
                PolygonPupil::new(self.aperture.radius_m(), vertices)
                    .sample_phase_range([diffraction_sample.u(), diffraction_sample.v()], phase)
            }
            _ => {
                let (radius, proposal) = annular::sample_radius(
                    self.aperture.obscuration_ratio(),
                    phase,
                    diffraction_sample.u(),
                );
                let (sin, cos) = det::sin_cos(2.0 * PI * diffraction_sample.v());
                ([radius * cos, radius * sin], proposal)
            }
        };
        let scale = (reference_nm * 1.0e-9 / self.aperture.radius_m()) / (2.0 * PI);
        let x = x_tan + scale * q[0];
        let y = y_tan + scale * q[1];
        if !scale.is_finite()
            || scale <= 0.0
            || !proposal.is_finite()
            || proposal <= 0.0
            || !x.is_finite()
            || !y.is_finite()
        {
            return Err(CameraError::InvalidDiffraction);
        }
        let ray = Ray {
            origin: self.eye,
            dir: direct_unit(Vec3::new(
                self.forward.x + x * self.right.x + y * self.up.x,
                self.forward.y + x * self.right.y + y * self.up.y,
                self.forward.z + x * self.right.z + y * self.up.z,
            ))?,
        };
        Ok(LayeredCameraRay {
            camera: self,
            ray,
            q,
            proposal,
            reference_nm,
            wavelengths_nm: *wavelengths_nm,
        })
    }

    /// Angular-raster PSF density evaluated at an admitted actual source depth.
    /// Its integral is one for each depth and wavelength. Summing source
    /// contributions therefore preserves their image-integrated radiance.
    pub fn layered_diffraction_density(
        &self,
        tangent_offset: [f64; 2],
        source_depth_m: f64,
        wavelength_nm: f64,
    ) -> Result<f64, CameraError> {
        self.diffraction_wavelength_range(&[wavelength_nm])?;
        self.admit_layered_depth(source_depth_m, wavelength_nm)?;
        let scale = (wavelength_nm * 1.0e-9 / self.aperture.radius_m()) / (2.0 * PI);
        let q = tangent_offset.map(|value| value / scale);
        let phase = self.pupil_defocus_phase(1.0 / source_depth_m, wavelength_nm)?;
        if !scale.is_finite() || scale <= 0.0 || q.iter().any(|value| !value.is_finite()) {
            return Err(CameraError::InvalidDiffraction);
        }
        let density = self.layered_pupil_density(q, phase)? / scale / scale;
        if !density.is_finite() || density < 0.0 {
            return Err(CameraError::InvalidDiffraction);
        }
        Ok(density)
    }

    pub(super) fn layered_maximum_phase(&self, wavelength_nm: f64) -> Result<f64, CameraError> {
        let policy = self
            .layered_diffraction
            .ok_or(CameraError::InvalidDiffraction)?;
        // Phase is affine in inverse source depth. Its maximum magnitude on
        // this ordered set occurs at an endpoint, so packet admission and
        // proposal selection require two phase evaluations, not one per plane.
        let depths = policy.object_distances_m();
        let nearest = self.pupil_defocus_phase(1.0 / depths[0], wavelength_nm)?;
        let farthest = self.pupil_defocus_phase(1.0 / depths[depths.len() - 1], wavelength_nm)?;
        Ok(nearest.abs().max(farthest.abs()))
    }

    pub(super) fn admit_layered_depth(
        &self,
        depth: f64,
        wavelength_nm: f64,
    ) -> Result<(), CameraError> {
        let policy = self
            .layered_diffraction
            .ok_or(CameraError::InvalidDiffraction)?;
        if !depth.is_finite() || depth <= 0.0 {
            return Err(CameraError::DefocusOutsideDiffractionModel);
        }
        self.pupil_defocus_phase(1.0 / depth, wavelength_nm)?;
        let inverse_error = policy
            .object_distances_m()
            .iter()
            .map(|declared| (1.0 / depth - 1.0 / declared).abs())
            .fold(f64::INFINITY, f64::min);
        let phase_error = PI
            * (self.aperture.radius_m() / (wavelength_nm * 1.0e-9))
            * self.aperture.radius_m()
            * inverse_error;
        if !phase_error.is_finite() || phase_error > policy.maximum_depth_phase_error_rad {
            return Err(CameraError::DefocusOutsideDiffractionModel);
        }
        Ok(())
    }

    fn layered_pupil_density(&self, q: [f64; 2], phase: f64) -> Result<f64, CameraError> {
        let radius = det::hypot(q[0], q[1]);
        if !radius.is_finite() {
            return Err(CameraError::InvalidDiffraction);
        }
        match &self.aperture.kind {
            ApertureKind::Circular { .. } => Ok(CircularDefocus::new(phase).density(radius)),
            ApertureKind::Annular { .. } => {
                Ok(AnnularPupil::new(self.aperture.obscuration_ratio(), phase).density(radius))
            }
            ApertureKind::RegularPolygon { vertices, .. } => {
                PolygonPupil::new(self.aperture.radius_m(), vertices)
                    .defocused_density(q, phase)
                    .ok_or(CameraError::InvalidDiffraction)
            }
            ApertureKind::Pinhole => Err(CameraError::InvalidDiffraction),
        }
    }
}

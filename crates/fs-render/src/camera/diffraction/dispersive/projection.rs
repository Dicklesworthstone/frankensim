//! Exact local inverse of the original material-driven geometric thin lens.
//!
//! Fixed-size DD basis contractions retain the physical camera's authored
//! world frame. The final original ray is always checked against the target;
//! the supplied metre tolerance is a geometric residual, not a wave-optics
//! certificate, visibility test, or pupil-averaged image.

use super::{CameraError, PhysicalCamera};
use crate::camera::LensSample;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::{dd::Dd, det};

/// One wavelength- and pupil-specific geometric thin-lens projection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeometricSpectralProjection {
    tangent_xy: [f64; 2],
    ndc_xy: [f64; 2],
    depth_m: f64,
    residual_m: f64,
}

impl GeometricSpectralProjection {
    pub const fn tangent_xy(self) -> [f64; 2] {
        self.tangent_xy
    }
    /// Original raster convention: +x right, +y up; off-screen values remain.
    pub const fn ndc_xy(self) -> [f64; 2] {
        self.ndc_xy
    }
    pub const fn depth_m(self) -> f64 {
        self.depth_m
    }
    /// Closest forward distance between the original ray and target [m].
    pub const fn residual_m(self) -> f64 {
        self.residual_m
    }
}

impl PhysicalCamera {
    /// Invert one original geometric dispersive pupil realization analytically.
    ///
    /// The same exact phase-index laws, exterior medium, fixed sensor conjugate
    /// and aperture sampler determine both directions. This method performs no
    /// iterative optical solve and allocates no query storage. It verifies the
    /// result by generating one original spectral ray; source bands, optical
    /// opening and full-pupil paraxial limits remain with that native owner.
    ///
    /// The reference/achromatic limit retains the original world-space focus
    /// construction, including the actual rounded pupil origin. Positive target
    /// depth and a finite positive metre tolerance are required. Diffraction,
    /// sequential lenses, unsupported spectra and unresolved residuals refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn project_geometric_spectral_point(
        &self,
        cx: &Cx<'_>,
        target: Point3,
        aspect_ratio: f64,
        wavelength_nm: f64,
        lens_sample: LensSample,
        position_tolerance_m: f64,
    ) -> Result<GeometricSpectralProjection, CameraError> {
        cx.checkpoint()?;
        if !aspect_ratio.is_finite()
            || aspect_ratio <= 0.0
            || !position_tolerance_m.is_finite()
            || position_tolerance_m <= 0.0
        {
            return Err(CameraError::InvalidProjection);
        }
        if ![target.x, target.y, target.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(CameraError::NonFiniteInput);
        }
        if self.has_diffraction() || self.sequential_lens().is_some() {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let lens = self
            .dispersive_lens()
            .ok_or(CameraError::InvalidDispersiveLens)?;
        self.diffraction_wavelength_range(&[wavelength_nm])?;
        let difference = lens.power_difference_per_m(wavelength_nm)?;
        let [px, py] = self.aperture.sample(lens_sample);
        // Identical operation order to the original geometric ray origin.
        let origin = self
            .eye
            .offset(self.right.scale(px))
            .offset(self.up.scale(py));
        let depth = dot(difference3(target, self.eye), vector(self.forward)).to_f64();
        if !depth.is_finite() || depth <= 0.0 {
            return Err(CameraError::DispersiveLensOutsideModel);
        }

        // Reciprocal basis, not an assumed exactly orthogonal rounded frame.
        let [r, u, f] = [self.right, self.up, self.forward].map(vector);
        let dual = [cross(u, f), cross(f, r), cross(r, u)];
        let determinant = dot(r, dual[0]);
        if determinant == Dd::ZERO || !finite(determinant) {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let delta = difference3(target, origin);
        let coordinates = dual.map(|axis| dot(delta, axis) / determinant);
        if !coordinates.into_iter().all(finite) || !Dd::ZERO.lt(coordinates[2]) {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let tangent_xy = if difference == 0.0 {
            // The native achromatic limit targets a finite world focus point.
            // Intersect the target line with that same authored basis plane.
            let pupil = dual.map(|axis| dot(difference3(origin, self.eye), axis) / determinant);
            let focus = Dd::from_f64(self.focus_distance_m);
            let along = (focus - pupil[2]) / coordinates[2];
            [
                ((pupil[0] + along * coordinates[0]) / focus).to_f64(),
                ((pupil[1] + along * coordinates[1]) / focus).to_f64(),
            ]
        } else {
            let inverse_conjugate = Dd::from_f64(difference - 1.0 / self.focus_distance_m);
            [
                (coordinates[0] / coordinates[2] - Dd::from_f64(px) * inverse_conjugate).to_f64(),
                (coordinates[1] / coordinates[2] - Dd::from_f64(py) * inverse_conjugate).to_f64(),
            ]
        };
        if !tangent_xy.into_iter().all(f64::is_finite) {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let [ray] = self.generate_geometric_spectral_rays_from_tangent_offsets(
            cx,
            tangent_xy[0],
            tangent_xy[1],
            lens_sample,
            &[wavelength_nm],
        )?;
        let direction = vector(ray.dir);
        let delta = difference3(target, ray.origin);
        let along = dot(delta, direction) / dot(direction, direction);
        if !finite(along) || along.lt(Dd::ZERO) {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let miss: [f64; 3] = core::array::from_fn(|i| (delta[i] - along * direction[i]).to_f64());
        let residual_m = det::hypot(det::hypot(miss[0], miss[1]), miss[2]);
        if !residual_m.is_finite() || residual_m > position_tolerance_m {
            return Err(CameraError::DispersiveLensOutsideModel);
        }
        let half = self.raster_vertical_half_tan();
        let ndc_xy = [tangent_xy[0] / aspect_ratio / half, tangent_xy[1] / half];
        if !ndc_xy.into_iter().all(f64::is_finite) {
            return Err(CameraError::InvalidProjection);
        }
        cx.checkpoint()?;
        Ok(GeometricSpectralProjection {
            tangent_xy,
            ndc_xy,
            depth_m: depth,
            residual_m,
        })
    }
}

fn vector(v: Vec3) -> [Dd; 3] {
    [v.x, v.y, v.z].map(Dd::from_f64)
}
fn difference3(a: Point3, b: Point3) -> [Dd; 3] {
    [
        Dd::from_f64(a.x) - Dd::from_f64(b.x),
        Dd::from_f64(a.y) - Dd::from_f64(b.y),
        Dd::from_f64(a.z) - Dd::from_f64(b.z),
    ]
}
fn dot(a: [Dd; 3], b: [Dd; 3]) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

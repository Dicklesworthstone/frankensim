//! Finite rectangular pixels mapped from actual native detector contacts.
//!
//! The raster is centred on its local XY plane, with +X right and +Y up.
//! Row zero is the top row. Each point belongs to a half-open cell in image
//! coordinates: left/top edges are included and right/bottom edges excluded.
//! This is a numerical classification of the original surface witness, not a
//! virtual detector intersection or relocation of the propagated ray.
//!
//! Binding checks the detector's material and original geometry/frame identity.
//! Each observed contact must have native planar-patch support, an aligned
//! normal, and a plane residual within the declared registration tolerance,
//! numerical contact collar and transformation roundoff. This local admission
//! is not a global planarity certificate or a bound on pixel-assignment error.
//! Detector hits outside the raster remain explicit cropped contributions.

use super::{DetectorBinding, MeterError, validate_spectral_sample};
use crate::instances::RigidTransform;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::{CurvedPathResult, CurvedPolarizedScene};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

/// Maximum admitted number of pixels in one native detector raster.
pub const MAX_DETECTOR_PIXELS: usize = 1_048_576;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RasterError {
    InvalidLayout,
    PixelBudgetExceeded,
    Meter(MeterError),
    /// A contact has no native open planar-patch witness.
    NonPlanarContact,
    /// The native patch normal or position does not match the sensor plane.
    PlaneMismatch,
    /// A surviving limited path cannot be classified as a detector miss.
    UnresolvedPath(PathTermination),
    NumericalRange,
}
impl From<MeterError> for RasterError {
    fn from(error: MeterError) -> Self {
        Self::Meter(error)
    }
}
impl core::fmt::Display for RasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native detector raster: {self:?}")
    }
}
impl core::error::Error for RasterError {}

/// Physical rectangular sensor layout; pixels are area cells, not point probes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DetectorRaster {
    transform: RigidTransform,
    inverse_rotation: RigidTransform,
    width_m: f64,
    height_m: f64,
    cols: u32,
    rows: u32,
    len: usize,
    pixel_area_m2: f64,
    plane_tolerance_m: f64,
}
impl DetectorRaster {
    /// Admit a centred local XY rectangle with finite, representable pixel area.
    /// Registration tolerance starts at zero; native contact and floating-point
    /// evaluation collars still apply. Use `with_plane_tolerance_m` to declare
    /// an additional sensor-to-surface registration tolerance explicitly.
    pub fn try_new(
        transform: RigidTransform,
        width_m: f64,
        height_m: f64,
        width_pixels: u32,
        height_pixels: u32,
    ) -> Result<Self, RasterError> {
        if !width_m.is_finite()
            || width_m <= 0.0
            || !height_m.is_finite()
            || height_m <= 0.0
            || width_pixels == 0
            || height_pixels == 0
        {
            return Err(RasterError::InvalidLayout);
        }
        let count = u64::from(width_pixels) * u64::from(height_pixels);
        if count > MAX_DETECTOR_PIXELS as u64 {
            return Err(RasterError::PixelBudgetExceeded);
        }
        let pixel_area_m2 =
            (width_m / f64::from(width_pixels)) * (height_m / f64::from(height_pixels));
        if !pixel_area_m2.is_finite() || pixel_area_m2 <= 0.0 {
            return Err(RasterError::NumericalRange);
        }
        if 0.5 * width_m == 0.0 || 0.5 * height_m == 0.0 {
            return Err(RasterError::NumericalRange);
        }
        let [x, y, z, w] = transform.rotation_xyzw();
        let inverse_rotation = RigidTransform::try_new([-x, -y, -z, w], [0.0; 3])
            .map_err(|_| RasterError::NumericalRange)?;
        for x in [-0.5 * width_m, 0.5 * width_m] {
            for y in [-0.5 * height_m, 0.5 * height_m] {
                let p = transform.transform_point(Point3::new(x, y, 0.0));
                if !finite_point(p) {
                    return Err(RasterError::NumericalRange);
                }
            }
        }
        Ok(Self {
            transform,
            inverse_rotation,
            width_m,
            height_m,
            cols: width_pixels,
            rows: height_pixels,
            len: count as usize,
            pixel_area_m2,
            plane_tolerance_m: 0.0,
        })
    }

    pub fn with_plane_tolerance_m(mut self, tolerance_m: f64) -> Result<Self, RasterError> {
        if !tolerance_m.is_finite() || tolerance_m < 0.0 {
            return Err(RasterError::InvalidLayout);
        }
        self.plane_tolerance_m = tolerance_m;
        Ok(self)
    }

    pub const fn transform(self) -> RigidTransform {
        self.transform
    }
    pub const fn width_m(self) -> f64 {
        self.width_m
    }
    pub const fn height_m(self) -> f64 {
        self.height_m
    }
    pub const fn cols(self) -> u32 {
        self.cols
    }
    pub const fn rows(self) -> u32 {
        self.rows
    }
    pub const fn len(self) -> usize {
        self.len
    }
    pub const fn is_empty(self) -> bool {
        false
    }
    pub const fn pixel_area_m2(self) -> f64 {
        self.pixel_area_m2
    }
    pub const fn plane_tolerance_m(self) -> f64 {
        self.plane_tolerance_m
    }

    fn local_point(self, point: Point3) -> Result<Vec3, RasterError> {
        if !finite_point(point) {
            return Err(RasterError::NumericalRange);
        }
        let [x, y, z] = self.transform.translation_m();
        // Subtract before rotating, avoiding cancellation between two large
        // separately rotated world positions for a translated sensor.
        let local = self.inverse_rotation.transform_vector(Vec3::new(
            point.x - x,
            point.y - y,
            point.z - z,
        ));
        if !finite_vector(local) {
            return Err(RasterError::NumericalRange);
        }
        Ok(local)
    }

    fn classify(
        self,
        point: Point3,
        normal: Vec3,
        resolution_m: f64,
    ) -> Result<Option<usize>, RasterError> {
        let local = self.local_point(point)?;
        let normal = self.inverse_rotation.transform_vector(normal);
        if !finite_vector(normal) || !resolution_m.is_finite() || resolution_m < 0.0 {
            return Err(RasterError::NumericalRange);
        }
        if normal.x.hypot(normal.y) > 1e-8 || (normal.z.abs() - 1.0).abs() > 1e-8 {
            return Err(RasterError::PlaneMismatch);
        }
        let scale = local
            .x
            .abs()
            .max(local.y.abs())
            .max(local.z.abs())
            .max(self.width_m)
            .max(self.height_m);
        let tolerance = self.plane_tolerance_m + resolution_m + 64.0 * f64::EPSILON * scale;
        if !tolerance.is_finite() {
            return Err(RasterError::NumericalRange);
        }
        if local.z.abs() > tolerance {
            return Err(RasterError::PlaneMismatch);
        }
        // Test physical edges first: a rounded division must not pull an
        // outside point into the raster. Right/bottom belong to the crop.
        if local.x < -0.5 * self.width_m
            || local.x >= 0.5 * self.width_m
            || local.y > 0.5 * self.height_m
            || local.y <= -0.5 * self.height_m
        {
            return Ok(None);
        }
        let column = ((local.x / self.width_m + 0.5) * f64::from(self.cols)).floor();
        let row = ((0.5 - local.y / self.height_m) * f64::from(self.rows)).floor();
        if !column.is_finite() || !row.is_finite() || column < 0.0 || row < 0.0 {
            return Err(RasterError::NumericalRange);
        }
        // An interior coordinate can round to the upper normalized edge.
        // Its physical edge test above already establishes the final cell.
        let column = (column as usize).min(self.cols as usize - 1);
        let row = (row as usize).min(self.rows as usize - 1);
        Ok(Some(row * self.cols as usize + column))
    }
}

/// Immutable physical detector binding with a separate, explicit raster crop.
/// No per-bin accumulation storage or geometry copy is retained here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundDetectorRaster {
    detector: DetectorBinding,
    layout: DetectorRaster,
}

/// One original hit contribution, already expressed in the native detector's
/// reference frame. `pixel=None` is a hit outside the raster, not a missed ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RasterContribution {
    pub pixel: Option<usize>,
    pub stokes_w_nm: Stokes,
    pub contact_resolution_m: f64,
}

impl CurvedPolarizedScene<'_> {
    pub fn detector_raster(
        &self,
        object_id: u64,
        layout: DetectorRaster,
    ) -> Result<BoundDetectorRaster, RasterError> {
        Ok(BoundDetectorRaster {
            detector: DetectorBinding::for_detector(self.geometry, &self.bindings, object_id)?,
            layout,
        })
    }
}

impl BoundDetectorRaster {
    pub const fn layout(self) -> DetectorRaster {
        self.layout
    }
    pub const fn object_id(self) -> u64 {
        self.detector.object_id
    }
    pub const fn frame_identity(self) -> ContentHash {
        self.detector.frame_identity
    }

    /// Extract one receiving contribution without changing a counter or meter.
    /// All launches must still be included in the consumer's normalization.
    pub fn contribution(
        &self,
        bin: usize,
        path: &CurvedPathResult,
    ) -> Result<Option<RasterContribution>, RasterError> {
        validate_spectral_sample(bin, path.beam(), path.max_contact_resolution_m())?;
        if matches!(
            path.termination(),
            PathTermination::DistanceLimit
                | PathTermination::SurfaceLimit
                | PathTermination::VolumeLimit
        ) && path.beam().power().i > 0.0
        {
            return Err(RasterError::UnresolvedPath(path.termination()));
        }
        if path.termination() != PathTermination::Detected {
            return Ok(None);
        }
        let stokes = self.detector.detected_power(
            path.beam(),
            path.last_contact()
                .map(|contact| (contact.object_id, contact.frame_identity)),
        )?;
        let Some(stokes_w_nm) = stokes else {
            return Ok(None);
        };
        let contact = path
            .last_contact()
            .ok_or(MeterError::DetectorGeometryMismatch)?;
        if contact.geometry_identity != self.detector.geometry_identity {
            return Err(MeterError::DetectorGeometryMismatch.into());
        }
        if contact.witness.departure_radius_m.is_none() {
            return Err(RasterError::NonPlanarContact);
        }
        let pixel = self.layout.classify(
            contact.witness.surface_point,
            contact.witness.normal,
            contact.resolution_m,
        )?;
        Ok(Some(RasterContribution {
            pixel,
            stokes_w_nm,
            contact_resolution_m: path.max_contact_resolution_m(),
        }))
    }
}

fn finite_point(p: Point3) -> bool {
    [p.x, p.y, p.z].iter().all(|v| v.is_finite())
}
fn finite_vector(v: Vec3) -> bool {
    [v.x, v.y, v.z].iter().all(|v| v.is_finite())
}

#[cfg(test)]
mod tests;

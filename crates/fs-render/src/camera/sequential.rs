//! Exact sequential refraction through an authored, rotationally symmetric lens.
//!
//! The sensor is at z=0 and light is traced backwards toward positive z. Each
//! surface is a spherical cap (signed curvature) or a plane, with an explicit
//! clear aperture and the material on its scene side. Sensor and scene media
//! are index one. Only the direct transmitted image is modeled: reflected
//! ghosts, diffraction, aspheres and automatic focusing are separate models.
//!
//! We sample the actual rear cap, not a guessed exit-pupil bounding box. The
//! projected-area to sensor-solid-angle Jacobian retains obliquity and clipping
//! without dropping skew rays. See PBRT, *Realistic Cameras*, sections 6.4.2
//! and 6.4.7, https://pbr-book.org/3ed-2018/Camera_Models/Realistic_Cameras.

use std::sync::Arc;

use super::{Aperture, CameraError, CameraProjection, LensSample, PhysicalCamera};
use crate::charts::Ray;
use crate::dielectric::{DielectricGlass, refract};
use crate::polarization::Stokes;
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// A sphere cap or plane, ordered from sensor to scene in a lens prescription.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LensSurface {
    vertex_z_m: f64,
    curvature_per_m: f64,
    clear_radius_m: f64,
    medium_after: Option<DielectricGlass>,
}

impl LensSurface {
    /// Positive curvature places the sphere centre on the scene side of the
    /// vertex. Zero is a plane; a plane with unchanged material is an aperture
    /// stop. `None` explicitly denotes index-one vacuum, not a fitted air law.
    pub fn try_new(
        vertex_z_m: f64,
        curvature_per_m: f64,
        clear_radius_m: f64,
        medium_after: Option<DielectricGlass>,
    ) -> Result<Self, CameraError> {
        if ![vertex_z_m, curvature_per_m, clear_radius_m]
            .iter()
            .all(|v| v.is_finite())
            || vertex_z_m <= 0.0
            || clear_radius_m <= 0.0
            || !(clear_radius_m * clear_radius_m > 0.0)
            || !(clear_radius_m * clear_radius_m).is_finite()
            || (curvature_per_m * clear_radius_m).abs() >= 1.0
        {
            return Err(CameraError::InvalidSequentialLens);
        }
        let surface = Self {
            vertex_z_m,
            curvature_per_m,
            clear_radius_m,
            medium_after,
        };
        let edge = vertex_z_m + surface.sag(clear_radius_m * clear_radius_m);
        if !edge.is_finite() || edge <= 0.0 {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(surface)
    }

    pub const fn vertex_z_m(self) -> f64 {
        self.vertex_z_m
    }
    pub const fn curvature_per_m(self) -> f64 {
        self.curvature_per_m
    }
    pub const fn clear_radius_m(self) -> f64 {
        self.clear_radius_m
    }
    pub const fn medium_after(&self) -> Option<&DielectricGlass> {
        self.medium_after.as_ref()
    }

    fn sag(self, radius_squared: f64) -> f64 {
        let c = self.curvature_per_m;
        // Rationalized sphere sag retains the plane limit for weak curvature.
        let numerator = c * radius_squared;
        numerator / (1.0 + det::sqrt((1.0 - c * numerator).max(0.0)))
    }

    fn bounds(self) -> [f64; 2] {
        let edge = self.vertex_z_m + self.sag(self.clear_radius_m * self.clear_radius_m);
        [self.vertex_z_m.min(edge), self.vertex_z_m.max(edge)]
    }

    fn normal(self, p: Point3) -> Result<Vec3, CameraError> {
        super::scaled_unit(Vec3::new(
            self.curvature_per_m * p.x,
            self.curvature_per_m * p.y,
            self.curvature_per_m * (p.z - self.vertex_z_m) - 1.0,
        ))
    }

    fn intersect(self, ray: Ray) -> Result<Option<(Point3, f64)>, CameraError> {
        let o = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z - self.vertex_z_m);
        let c = self.curvature_per_m;
        let mut roots = [f64::INFINITY; 2];
        if c == 0.0 {
            if ray.dir.z == 0.0 {
                return Ok(None);
            }
            roots[0] = -o.z / ray.dir.z;
        } else {
            // c*(x²+y²+dz²)-2*dz=0 avoids a huge sphere centre/radius
            // subtraction when curvature tends to zero. Stable quadratic roots
            // retain the near cap; the opposite hemisphere is never a lens face.
            let a = c * ray.dir.dot(ray.dir);
            let b = c * o.dot(ray.dir) - ray.dir.z;
            let constant = c * o.dot(o) - 2.0 * o.z;
            let discriminant = b.mul_add(b, -a * constant);
            if ![a, b, constant, discriminant].iter().all(|v| v.is_finite()) {
                return Err(CameraError::SequentialLensOutsideModel);
            }
            if discriminant < 0.0 {
                return Ok(None);
            }
            let q = -b - det::sqrt(discriminant).copysign(b);
            if q == 0.0 {
                roots[0] = -b / a;
            } else {
                roots = [q / a, constant / q];
            }
        }
        let mut hit = None;
        for distance in roots {
            if !distance.is_finite() || distance <= 0.0 {
                continue;
            }
            let p = ray.origin.offset(ray.dir.scale(distance));
            if ![p.x, p.y, p.z].iter().all(|v| v.is_finite()) {
                return Err(CameraError::SequentialLensOutsideModel);
            }
            if c * (p.z - self.vertex_z_m) > 1.0 {
                continue;
            }
            if det::hypot(p.x, p.y) > self.clear_radius_m {
                continue;
            }
            if hit.is_none_or(|(_, previous)| distance < previous) {
                hit = Some((p, distance));
            }
        }
        Ok(hit)
    }
}

/// The direct transmitted ray, with power losses retained along its real path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialLensRay {
    pub ray: Ray,
    /// Unpolarized direct-image transmission. Polarization generated at one
    /// interface is retained through the others; Fresnel averages are not
    /// independently multiplied. Reciprocity gives the same unpolarized
    /// transmission from scene to sensor. Both exterior indices are one.
    pub transmittance: f64,
    /// Sensor-to-front-surface geometric path in metres.
    pub distance_m: f64,
    /// Integral of the phase index over that path, in metres.
    pub optical_path_m: f64,
}

/// One wavelength's outgoing camera ray and normalized radiance response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialCameraRay {
    pub ray: Ray,
    /// Direct sensor irradiance divided by `pi*rear_radius²/rear_vertex_z²`.
    /// Includes exact sensor/cap obliquity, Fresnel loss and bulk absorption.
    /// This explicit aperture-normalized convention matches the existing film's
    /// radiance scale; exposure duration, ISO and sensor quantum efficiency are
    /// not silently applied. Multiply by that solid-angle normalization to
    /// recover irradiance from a radiometric scene.
    pub weight: f64,
}

/// Immutable, bounded optical prescription; tracing performs no allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct SequentialLens {
    surfaces: Arc<[LensSurface]>,
    reference_wavelength_nm: f64,
    wavelength_bounds_nm: (f64, f64),
}

impl SequentialLens {
    pub const MAX_SURFACES: usize = 32;

    /// Admit ordered, disjoint axial cap bounds. This conservative condition
    /// makes surface order unambiguous without a general lens-solid Boolean
    /// model. Cemented elements share one interface with the next glass as its
    /// scene-side medium. The final surface must return to index one.
    pub fn try_new(
        reference_wavelength_nm: f64,
        surfaces: &[LensSurface],
    ) -> Result<Self, CameraError> {
        if surfaces.is_empty()
            || surfaces.len() > Self::MAX_SURFACES
            || !reference_wavelength_nm.is_finite()
            || reference_wavelength_nm <= 0.0
            || surfaces.last().is_some_and(|s| s.medium_after.is_some())
            || surfaces
                .windows(2)
                .any(|pair| pair[0].bounds()[1] >= pair[1].bounds()[0])
        {
            return Err(CameraError::InvalidSequentialLens);
        }
        let mut band = (0.0_f64, f64::INFINITY);
        for surface in surfaces {
            if let Some(glass) = surface.medium_after {
                let (lo, hi) = glass
                    .wavelength_bounds_nm()
                    .map_err(|_| CameraError::InvalidSequentialLens)?;
                band = (band.0.max(lo), band.1.min(hi));
                glass
                    .ior()
                    .eval(reference_wavelength_nm)
                    .map_err(|_| CameraError::InvalidSequentialLens)?;
                glass
                    .absorption_ref()
                    .extinction_per_m(reference_wavelength_nm)
                    .map_err(|_| CameraError::InvalidSequentialLens)?;
            }
        }
        if reference_wavelength_nm < band.0 || reference_wavelength_nm > band.1 {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(Self {
            surfaces: Arc::from(surfaces),
            reference_wavelength_nm,
            wavelength_bounds_nm: band,
        })
    }

    pub fn surfaces(&self) -> &[LensSurface] {
        &self.surfaces
    }
    pub const fn reference_wavelength_nm(&self) -> f64 {
        self.reference_wavelength_nm
    }
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }

    /// Paraxial reference focal length and axial object focus measured from
    /// the sensor. Exact rays do not use this approximation. The ABCD matrix
    /// propagates `(height, n*angle)` through real spacings and surface powers.
    /// Afocal/diverging systems or a virtual focus cannot attach to the
    /// current positive-focus camera, but remain usable through `trace_ray`.
    pub fn reference_focus(&self) -> Result<(f64, f64), CameraError> {
        let (mut a, mut b, mut c, mut d) = (1.0, 0.0, 0.0, 1.0);
        let (mut z, mut n) = (0.0, 1.0);
        for surface in self.surfaces() {
            let travel = (surface.vertex_z_m - z) / n;
            a += travel * c;
            b += travel * d;
            let next = surface
                .medium_after
                .map_or(Ok(1.0), |glass| {
                    glass.ior().eval(self.reference_wavelength_nm)
                })
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            let power = (n - next) * surface.curvature_per_m;
            c += power * a;
            d += power * b;
            z = surface.vertex_z_m;
            n = next;
        }
        let focal = -1.0 / c;
        let beyond_front = -b / d;
        let focus = z + beyond_front;
        if ![focal, beyond_front, focus]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        Ok((focal, focus))
    }

    /// Trace a ray from the sensor plane in local SI coordinates. Aperture
    /// misses, backward departures and TIR are `Ok(None)`. Invalid wavelengths,
    /// numerical range and cancellation remain errors, including for clipped
    /// rays. Every original material band is checked before traversal.
    pub fn trace_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
    ) -> Result<Option<SequentialLensRay>, CameraError> {
        cx.checkpoint()?;
        let (lo, hi) = self.wavelength_bounds_nm;
        if !wavelength_nm.is_finite()
            || wavelength_nm <= 0.0
            || wavelength_nm < lo
            || wavelength_nm > hi
            || ![ray.origin.x, ray.origin.y, ray.origin.z]
                .iter()
                .all(|v| v.is_finite())
            || ray.origin.z != 0.0
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let mut indices = [1.0; Self::MAX_SURFACES + 1];
        let mut extinction = [0.0; Self::MAX_SURFACES + 1];
        for (i, surface) in self.surfaces().iter().enumerate() {
            if let Some(glass) = surface.medium_after {
                indices[i + 1] = glass
                    .ior()
                    .eval(wavelength_nm)
                    .map_err(|_| CameraError::SequentialLensOutsideModel)?;
                extinction[i + 1] = glass
                    .absorption_ref()
                    .extinction_per_m(wavelength_nm)
                    .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            }
        }
        let mut ray = Ray {
            origin: ray.origin,
            dir: super::scaled_unit(ray.dir)?,
        };
        let axis = if ray.dir.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let mut frame = PolarizationFrame::try_new(ray.dir, axis)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        let mut power = Stokes::unpolarized(1.0);
        let (mut distance_m, mut optical_path_m) = (0.0, 0.0);
        for (i, surface) in self.surfaces().iter().enumerate() {
            cx.checkpoint()?;
            if ray.dir.z <= 0.0 {
                return Ok(None);
            }
            let Some((point, distance)) = surface.intersect(ray)? else {
                return Ok(None);
            };
            let normal = surface.normal(point)?;
            let cosine = (-normal.dot(ray.dir)).clamp(0.0, 1.0);
            if cosine == 0.0 {
                return Ok(None);
            }
            let response = fresnel_dielectric_polarized(cosine, indices[i], indices[i + 1])
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            if response.scalar.total_internal_reflection {
                return Ok(None);
            }
            let cross = super::cross(ray.dir, normal);
            let s = if cross.x == 0.0 && cross.y == 0.0 && cross.z == 0.0 {
                frame.horizontal()
            } else {
                super::scaled_unit(cross)?
            };
            let incident_frame = PolarizationFrame::try_new(ray.dir, s)
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            power = response
                .transmission
                .apply_diagonal_jones(frame.express_in(power, incident_frame));
            let attenuation = det::exp(-extinction[i] * distance);
            power = Stokes {
                i: power.i * attenuation,
                q: power.q * attenuation,
                u: power.u * attenuation,
                v: power.v * attenuation,
            };
            let direction = refract(
                ray.dir.scale(-1.0),
                normal,
                indices[i],
                indices[i + 1],
                response.scalar.transmitted_cosine,
            )
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            frame = PolarizationFrame::try_new(direction, s)
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            ray = Ray {
                origin: point,
                dir: direction,
            };
            distance_m += distance;
            optical_path_m += indices[i] * distance;
            if ![power.i, distance_m, optical_path_m]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
            {
                return Err(CameraError::SequentialLensOutsideModel);
            }
        }
        cx.checkpoint()?;
        Ok(Some(SequentialLensRay {
            ray,
            transmittance: power.i,
            distance_m,
            optical_path_m,
        }))
    }
}

impl PhysicalCamera {
    /// Replace the ideal lens with an explicitly positioned sequential train.
    /// The camera eye becomes the sensor centre. Sensor height and pose are
    /// retained; focal/focus diagnostics are derived from the real reference
    /// prescription. Changing focus requires a newly positioned prescription,
    /// rather than silently moving authored surfaces during an animation.
    pub fn with_sequential_lens(mut self, lens: SequentialLens) -> Result<Self, CameraError> {
        if self.has_diffraction() || self.dispersive_lens.is_some() {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let height = self
            .projection
            .sensor_height_m()
            .ok_or(CameraError::InvalidProjection)?;
        let (focal, focus) = lens.reference_focus()?;
        self.projection = CameraProjection::try_focal_sensor(focal, height)?;
        self.focus_distance_m = focus;
        self.aperture = Aperture::try_circular(lens.surfaces()[0].clear_radius_m)?;
        self.sequential_lens = Some(lens);
        Ok(self)
    }

    pub const fn sequential_lens(&self) -> Option<&SequentialLens> {
        self.sequential_lens.as_ref()
    }

    /// Generate independent wavelength paths using one uniform projected-disk
    /// sample of the actual rear cap. A rejected wavelength retains zero
    /// sample mass; callers must not resample it or copy another lane's ray.
    pub fn generate_sequential_spectral_rays_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<[Option<SequentialCameraRay>; N], CameraError> {
        cx.checkpoint()?;
        if N == 0 || !x_tan.is_finite() || !y_tan.is_finite() {
            return Err(CameraError::InvalidProjection);
        }
        let lens = self
            .sequential_lens
            .as_ref()
            .ok_or(CameraError::InvalidSequentialLens)?;
        let rear = lens.surfaces()[0];
        let sensor = Point3::new(-x_tan * rear.vertex_z_m, -y_tan * rear.vertex_z_m, 0.0);
        let [x, y] = self.aperture.sample(lens_sample);
        let target = Point3::new(x, y, rear.vertex_z_m + rear.sag(x * x + y * y));
        let displacement = target.delta_from(sensor);
        let distance = super::scaled_norm(displacement)?;
        let direction = super::scaled_unit(displacement)?;
        let normal = rear.normal(target)?;
        let projected_cosine = -normal.z;
        let incident_cosine = (-normal.dot(direction)).max(0.0);
        let geometry = direction.z * incident_cosine / projected_cosine
            * (rear.vertex_z_m / distance)
            * (rear.vertex_z_m / distance);
        if !geometry.is_finite() || geometry < 0.0 {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let mut result = [None; N];
        for (entry, &wavelength) in result.iter_mut().zip(wavelengths_nm) {
            if let Some(traced) = lens.trace_ray(
                cx,
                Ray {
                    origin: sensor,
                    dir: direction,
                },
                wavelength,
            )? {
                if geometry == 0.0 {
                    continue;
                }
                let p = traced.ray.origin;
                let d = traced.ray.dir;
                let origin = self
                    .eye
                    .offset(self.right.scale(p.x))
                    .offset(self.up.scale(p.y))
                    .offset(self.forward.scale(p.z));
                let direction = Vec3::new(
                    self.right.x * d.x + self.up.x * d.y + self.forward.x * d.z,
                    self.right.y * d.x + self.up.y * d.y + self.forward.y * d.z,
                    self.right.z * d.x + self.up.z * d.y + self.forward.z * d.z,
                );
                *entry = Some(SequentialCameraRay {
                    ray: Ray {
                        origin,
                        dir: super::scaled_unit(direction)?,
                    },
                    weight: geometry * traced.transmittance,
                });
            }
        }
        Ok(result)
    }
}

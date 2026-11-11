//! Exact sequential refraction through an authored, rotationally symmetric lens.
//!
//! The sensor is at z=0 and light is traced backwards toward positive z. Each
//! surface is a spherical cap, conic/even asphere or plane, with an explicit
//! clear aperture and the material on its scene side. Exterior media may be
//! authored explicitly; the original constructor retains index-one exteriors.
//! The direct transmitted image is the default. An explicit proposal also
//! admits incoherent two-reflection ghost paths through these same surfaces.
//! Higher-order ghosts, coherent bulk etalons, diffraction and scene-based
//! autofocus remain separate models.
//!
//! We sample the actual rear cap, not a guessed exit-pupil bounding box. The
//! projected-area to sensor-solid-angle Jacobian retains obliquity and clipping
//! without dropping skew rays. See PBRT, *Realistic Cameras*, sections 6.4.2
//! and 6.4.7, https://pbr-book.org/3ed-2018/Camera_Models/Realistic_Cameras.

use std::sync::Arc;

use super::{Aperture, CameraError, CameraProjection, LensSample, PhysicalCamera};
use crate::charts::Ray;
use crate::dielectric::{DielectricCoating, DielectricGlass, DielectricSurface, reflect, refract};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{PolarizationFrame, physical_roundoff};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::{dd::Dd, det};

const LIGHT_M_S: f64 = 299_792_458.0;

mod asphere;
mod coating;
pub use coating::{LensCoating, LensCoatingLayer};
mod ghost;
pub use ghost::LensGhostPath;
use ghost::{GhostSampling, PathSteps};
mod focus;
pub use focus::SequentialLensFocus;
mod paraxial;
mod projection;
pub use projection::{SequentialProjection, SequentialProjectionError, SequentialProjectionOptions};
pub use paraxial::ParaxialFocusSensitivity;

/// A rotationally symmetric surface, ordered from sensor to scene.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LensSurface {
    vertex_z_m: f64,
    curvature_per_m: f64,
    clear_radius_m: f64,
    medium_after: Option<DielectricGlass>,
    asphere: Option<asphere::AsphericSurface>,
    optical_surface: DielectricSurface,
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
            asphere: None,
            optical_surface: DielectricSurface::SMOOTH,
        };
        let edge = vertex_z_m + surface.sag(clear_radius_m * clear_radius_m);
        if !edge.is_finite() || edge <= 0.0 {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(surface)
    }

    /// An authored conic plus even-polynomial sag in SI metres:
    /// `z(r)-vertex = c*r²/(1+sqrt(1-(1+k)*c²*r²)) + A4*r⁴ + ... + A16*r¹⁶`.
    /// `c` is the vertex curvature [1/m], `k` the dimensionless conic constant,
    /// and the seven coefficients have units m^-3, m^-5, ..., m^-15.
    /// There is no quadratic departure term: the original paraxial reference
    /// continues to use the actual vertex curvature. A zero departure and
    /// zero conic constant retain the original sphere/plane arithmetic.
    ///
    /// Admission requires a finite single-valued conic over the clear disk
    /// and conservative axial bounds wholly in front of the sensor. Pure
    /// conics use analytic intersections. Polynomial rays use bounded,
    /// safeguarded iteration after a global slope bound proves a unique
    /// forward crossing; unqualified grazing/multiple-crossing cases refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_aspheric(
        vertex_z_m: f64,
        curvature_per_m: f64,
        conic_constant: f64,
        even_coefficients: [f64; 7],
        clear_radius_m: f64,
        medium_after: Option<DielectricGlass>,
    ) -> Result<Self, CameraError> {
        if conic_constant == 0.0 && even_coefficients == [0.0; 7] {
            return Self::try_new(vertex_z_m, curvature_per_m, clear_radius_m, medium_after);
        }
        if !vertex_z_m.is_finite() || vertex_z_m <= 0.0 {
            return Err(CameraError::InvalidSequentialLens);
        }
        let asphere = asphere::AsphericSurface::try_new(
            curvature_per_m,
            conic_constant,
            even_coefficients,
            clear_radius_m,
        )
        .ok_or(CameraError::InvalidSequentialLens)?;
        let surface = Self {
            vertex_z_m,
            curvature_per_m,
            clear_radius_m,
            medium_after,
            asphere: Some(asphere),
            optical_surface: DielectricSurface::SMOOTH,
        };
        let [low, high] = surface.bounds();
        if !low.is_finite() || !high.is_finite() || low <= 0.0 {
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
    /// Zero for the original spherical or planar surface.
    pub const fn conic_constant(&self) -> f64 {
        match &self.asphere {
            Some(profile) => profile.conic_constant(),
            None => 0.0,
        }
    }
    /// Original SI A4 through A16 coefficients; all zero for a sphere/plane.
    pub const fn even_aspheric_coefficients(&self) -> &[f64; 7] {
        match &self.asphere {
            Some(profile) => profile.even_coefficients(),
            None => &[0.0; 7],
        }
    }
    pub const fn medium_after(&self) -> Option<&DielectricGlass> {
        self.medium_after.as_ref()
    }

    /// Attach a passive coherent coating in **sensor-to-scene layer order**.
    /// Every lens surface uses this same convention, including the front face
    /// of an element. The original films' spectral interference, absorption,
    /// and s/p phase act on the complete reciprocal detector state. Snell
    /// geometry still uses the adjacent bulk indices; the locally parallel
    /// films are collapsed onto this face and do not displace the ray.
    ///
    /// Compact Cauchy coatings restrict the prescription to 380--780 nm.
    /// [`SequentialLens::with_spectral_coating`] can replace this coating with
    /// original Sellmeier, tabulated or absorbing films at the lens level.
    /// Nonzero film thickness requires steady generation: timed queries refuse
    /// because coherent coating group delay is not part of this lens model.
    #[must_use]
    pub const fn with_coating(mut self, coating: DielectricCoating) -> Self {
        self.optical_surface = self.optical_surface.with_coating(coating);
        self
    }

    /// Borrow this surface's compact coating, including its source identity.
    /// A lens-level spectral coating takes precedence during tracing.
    pub const fn coating(&self) -> Option<&DielectricCoating> {
        self.optical_surface.coating()
    }

    fn sag(&self, radius_squared: f64) -> f64 {
        if let Some(profile) = &self.asphere {
            return profile.sag(radius_squared);
        }
        let c = self.curvature_per_m;
        // Rationalized sphere sag retains the plane limit for weak curvature.
        let numerator = c * radius_squared;
        numerator / (1.0 + det::sqrt((1.0 - c * numerator).max(0.0)))
    }

    fn bounds(&self) -> [f64; 2] {
        if let Some(profile) = &self.asphere {
            let [low, high] = profile.sag_bounds();
            return [
                (self.vertex_z_m + low).next_down(),
                (self.vertex_z_m + high).next_up(),
            ];
        }
        let edge = self.vertex_z_m + self.sag(self.clear_radius_m * self.clear_radius_m);
        [self.vertex_z_m.min(edge), self.vertex_z_m.max(edge)]
    }

    fn normal(&self, p: Point3) -> Result<Vec3, CameraError> {
        if let Some(profile) = &self.asphere {
            let [x, y, z] = profile
                .normal([p.x, p.y])
                .ok_or(CameraError::SequentialLensOutsideModel)?;
            return super::scaled_unit(Vec3::new(x, y, z));
        }
        super::scaled_unit(Vec3::new(
            self.curvature_per_m * p.x,
            self.curvature_per_m * p.y,
            self.curvature_per_m * (p.z - self.vertex_z_m) - 1.0,
        ))
    }

    fn intersect(&self, ray: Ray) -> Result<Option<(Point3, f64)>, CameraError> {
        if let Some(profile) = &self.asphere {
            return profile
                .intersect(
                    self.vertex_z_m,
                    [ray.origin.x, ray.origin.y, ray.origin.z],
                    [ray.dir.x, ray.dir.y, ray.dir.z],
                )
                .map(|hit| hit.map(|(p, t)| (Point3::new(p[0], p[1], p[2]), t)))
                .map_err(|()| CameraError::SequentialLensOutsideModel);
        }
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

/// One transmitted image path, with power losses retained along its real path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialLensRay {
    pub ray: Ray,
    /// Unpolarized image-path transmission, including any prescribed ghost
    /// reflections. Polarization generated at one interface is retained through
    /// the others; Fresnel averages are not independently multiplied.
    /// Reciprocity gives the same unpolarized power transmission from scene to
    /// sensor. This excludes the endpoint
    /// radiance factor applied by the physical-camera generators.
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
    /// Authored local phase index required at the native scene launch.
    pub scene_phase_index: f64,
    /// Sampled sensor irradiance divided by `pi*rear_radius²/rear_vertex_z²`.
    /// Includes exact sensor/cap obliquity, Fresnel loss, bulk absorption and
    /// the sensor-to-scene phase-index-squared radiance factor, and the inverse
    /// selected path probability when two-reflection ghosts are enabled.
    /// This explicit aperture-normalized convention matches the existing film's
    /// radiance scale; exposure duration, ISO and sensor quantum efficiency are
    /// not silently applied. Multiply by that solid-angle normalization to
    /// recover irradiance from a radiometric scene.
    pub weight: f64,
}

/// A reciprocal detector response after transmission through the actual lens.
///
/// The frame points from sensor to scene. `power()` is the complete transported
/// Stokes detector weight, including Fresnel polarization and bulk absorption;
/// it is not an emitted physical source. Continue it through the reciprocal
/// scene transport, then multiply the measurement by `geometry_weight()` once.
/// The sensor-to-scene phase-index-squared radiance factor is already included.
/// Native scene continuation starts in [`Self::scene_phase_index`], so its
/// original endpoint conversion composes with this factor exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SequentialCameraBeam {
    origin: Point3,
    frame: PolarizationFrame,
    power: Stokes,
    wavelength_nm: f64,
    geometry_weight: f64,
    scene_phase_index: f64,
    group_delay_s: Option<Dd>,
}

impl SequentialCameraBeam {
    pub const fn ray(&self) -> Ray {
        Ray {
            origin: self.origin,
            dir: self.frame.direction(),
        }
    }
    pub const fn frame(&self) -> PolarizationFrame {
        self.frame
    }
    pub const fn power(&self) -> Stokes {
        self.power
    }
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    /// Authored scene-side phase index at this wavelength. A native scene
    /// continuation must admit this same local index at the actual launch.
    pub const fn scene_phase_index(&self) -> f64 {
        self.scene_phase_index
    }
    /// Sensor-to-front-surface group transit along the actual lens path,
    /// including repeated bulk segments of an admitted ghost.
    /// Only explicit timed generation evaluates group indices and returns
    /// `Some`; ordinary steady generation never requires their derivatives.
    pub const fn group_delay_s(&self) -> Option<Dd> {
        self.group_delay_s
    }
    /// Sensor/cap obliquity, projected-area sampling, and the inverse selected
    /// path probability when first-order ghost sampling is enabled. The same
    /// aperture normalization as [`SequentialCameraRay::weight`] is used;
    /// optical reflection/transmission and analyzer response are in [`Self::power`].
    pub const fn geometry_weight(&self) -> f64 {
        self.geometry_weight
    }
}

struct TracedLensBeam {
    ray: Ray,
    frame: PolarizationFrame,
    power: Stokes,
    distance_m: f64,
    optical_path_m: f64,
    sensor_phase_index: f64,
    scene_phase_index: f64,
    group_delay_s: Option<Dd>,
}

/// Immutable, bounded optical prescription; tracing performs no allocation.
#[derive(Clone, Debug, PartialEq)]
pub struct SequentialLens {
    surfaces: Arc<[LensSurface]>,
    spectral_coatings: Option<Arc<[Option<LensCoating>]>>,
    ghost_sampling: Option<GhostSampling>,
    sensor_medium: Option<DielectricGlass>,
    reference_wavelength_nm: f64,
    // Immutable bulk-only support allows replacing a visible compact coating
    // without accidentally retaining its old restriction in the new lens.
    bulk_wavelength_bounds_nm: (f64, f64),
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
        if surfaces.last().is_some_and(|s| s.medium_after.is_some()) {
            return Err(CameraError::InvalidSequentialLens);
        }
        Self::try_new_with_sensor_medium(reference_wavelength_nm, None, surfaces)
    }

    /// Admit an explicit medium from the sensor to the first surface. Each
    /// surface still declares its scene-side material; the final one is the
    /// scene exterior. `None` denotes vacuum independently on either side.
    /// All source bands intersect before admission. Steady optics requires
    /// phase index and absorption only, never a group-index derivative.
    pub fn try_new_with_sensor_medium(
        reference_wavelength_nm: f64,
        sensor_medium: Option<DielectricGlass>,
        surfaces: &[LensSurface],
    ) -> Result<Self, CameraError> {
        if surfaces.is_empty()
            || surfaces.len() > Self::MAX_SURFACES
            || !reference_wavelength_nm.is_finite()
            || reference_wavelength_nm <= 0.0
            || surfaces
                .windows(2)
                .any(|pair| pair[0].bounds()[1] >= pair[1].bounds()[0])
        {
            return Err(CameraError::InvalidSequentialLens);
        }
        let mut band = (0.0_f64, f64::INFINITY);
        for glass in sensor_medium
            .iter()
            .chain(surfaces.iter().filter_map(LensSurface::medium_after))
        {
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
        let bulk_wavelength_bounds_nm = band;
        if surfaces.iter().any(|surface| surface.coating().is_some()) {
            band = (
                band.0.max(crate::spectral::LAMBDA_MIN),
                band.1.min(crate::spectral::LAMBDA_MAX),
            );
        }
        if reference_wavelength_nm < band.0 || reference_wavelength_nm > band.1 {
            return Err(CameraError::InvalidSequentialLens);
        }
        Ok(Self {
            surfaces: Arc::from(surfaces),
            spectral_coatings: None,
            ghost_sampling: None,
            sensor_medium,
            reference_wavelength_nm,
            bulk_wavelength_bounds_nm,
            wavelength_bounds_nm: band,
        })
    }

    /// Attach original spectral films in sensor-to-scene order to one face.
    /// The immutable coating data are shared when a lens is cloned. This
    /// coating replaces any compact coating authored on the selected surface;
    /// it is not an additional interface or a second attenuation multiplier.
    /// The complete effective film/bulk support is recomputed, so replacing a
    /// compact coating can restore nonvisible wavelengths of the bulk laws.
    ///
    /// The face index and the reference wavelength must be admitted. Nonzero
    /// films require steady tracing; timed rays refuse before aperture clipping
    /// because coherent group residence delay is not part of this model.
    pub fn with_spectral_coating(
        mut self,
        surface_index: usize,
        coating: LensCoating,
    ) -> Result<Self, CameraError> {
        if surface_index >= self.surfaces.len() {
            return Err(CameraError::InvalidSequentialLens);
        }
        let mut band = self.bulk_wavelength_bounds_nm;
        for (index, surface) in self.surfaces.iter().enumerate() {
            let spectral = if index == surface_index {
                Some(&coating)
            } else {
                self.spectral_coating(index)
            };
            let support = match spectral {
                Some(coating) => coating.wavelength_bounds_nm(),
                None if surface.coating().is_some() => {
                    (crate::spectral::LAMBDA_MIN, crate::spectral::LAMBDA_MAX)
                }
                None => continue,
            };
            band = (band.0.max(support.0), band.1.min(support.1));
        }
        if self.reference_wavelength_nm < band.0 || self.reference_wavelength_nm > band.1 {
            return Err(CameraError::InvalidSequentialLens);
        }
        let surface_count = self.surfaces.len();
        let coatings = self
            .spectral_coatings
            .get_or_insert_with(|| vec![None; surface_count].into());
        Arc::make_mut(coatings)[surface_index] = Some(coating);
        self.wavelength_bounds_nm = band;
        Ok(self)
    }

    /// Borrow the full spectral coating that takes precedence over the
    /// surface's compact coating, if one was attached to this face.
    #[must_use]
    pub fn spectral_coating(&self, surface_index: usize) -> Option<&LensCoating> {
        self.spectral_coatings
            .as_ref()?
            .get(surface_index)?
            .as_ref()
    }

    /// Enable the incoherent sum of the direct image and every two-interface
    /// reflection ghost pair in this prescription. The direct probability is
    /// finite and strictly between zero and one; the remaining mass is uniform
    /// over all N*(N-1)/2 pairs. It controls sampling variance. Native Fresnel,
    /// coating and Beer laws determine each path's physical optical strength.
    ///
    /// Physical-camera generators partition and remap the first lens-sample
    /// coordinate into path choice plus a conditionally uniform pupil sample.
    /// Their measure weight includes the inverse path probability exactly once.
    /// Misses keep zero sample mass. A path visits at most 3N-2 surfaces, with
    /// no per-ray allocation or intensity cutoff. Higher-reflection paths and
    /// coherent interference between bulk paths are outside this model.
    pub fn with_two_reflection_ghosts(
        mut self,
        direct_probability: f64,
    ) -> Result<Self, CameraError> {
        self.ghost_sampling = Some(GhostSampling::try_new(
            self.surfaces.len(),
            direct_probability,
        )?);
        Ok(self)
    }

    /// Restore the original direct-only camera proposal.
    #[must_use]
    pub fn without_ghosts(mut self) -> Self {
        self.ghost_sampling = None;
        self
    }

    /// Direct-path probability of the explicit two-reflection proposal.
    /// `None` retains the original direct-only model and sampling coordinates.
    #[must_use]
    pub fn ghost_direct_probability(&self) -> Option<f64> {
        self.ghost_sampling.map(GhostSampling::direct_probability)
    }

    pub fn surfaces(&self) -> &[LensSurface] {
        &self.surfaces
    }
    pub const fn sensor_medium(&self) -> Option<&DielectricGlass> {
        self.sensor_medium.as_ref()
    }
    pub fn scene_medium(&self) -> Option<&DielectricGlass> {
        self.surfaces.last().and_then(LensSurface::medium_after)
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
    /// The effective focal length is on the sensor side; the final object
    /// distance uses the authored scene-side index.
    /// Afocal/diverging systems or a virtual focus cannot attach to the
    /// current positive-focus camera, but remain usable through `trace_ray`.
    pub fn reference_focus(&self) -> Result<(f64, f64), CameraError> {
        self.paraxial_focus_at(self.reference_wavelength_nm)
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
        Ok(self
            .trace_seeded_ray(
                cx,
                ray,
                wavelength_nm,
                None,
                Stokes::unpolarized(1.0),
                false,
            )?
            .map(|beam| SequentialLensRay {
                ray: beam.ray,
                transmittance: beam.power.i,
                distance_m: beam.distance_m,
                optical_path_m: beam.optical_path_m,
            }))
    }

    /// Trace one deterministic two-reflection ghost without a sampling PDF.
    /// Native spectral interface power and every actual bulk segment are
    /// retained. This diagnostic does not depend on the camera proposal;
    /// [`Self::trace_ray`] continues to return the direct transmitted path.
    pub fn trace_ghost_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        path: LensGhostPath,
    ) -> Result<Option<SequentialLensRay>, CameraError> {
        Ok(self
            .trace_seeded_path(
                cx,
                ray,
                wavelength_nm,
                None,
                Stokes::unpolarized(1.0),
                false,
                Some(path),
            )?
            .map(|beam| SequentialLensRay {
                ray: beam.ray,
                transmittance: beam.power.i,
                distance_m: beam.distance_m,
                optical_path_m: beam.optical_path_m,
            }))
    }

    fn trace_seeded_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        reference_axis: Option<Vec3>,
        power: Stokes,
        timed: bool,
    ) -> Result<Option<TracedLensBeam>, CameraError> {
        self.trace_seeded_path(cx, ray, wavelength_nm, reference_axis, power, timed, None)
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_seeded_path(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        reference_axis: Option<Vec3>,
        mut power: Stokes,
        timed: bool,
        ghost: Option<LensGhostPath>,
    ) -> Result<Option<TracedLensBeam>, CameraError> {
        cx.checkpoint()?;
        if ghost.is_some_and(|path| path.first_surface_index() >= self.surfaces.len()) {
            return Err(CameraError::InvalidSequentialLens);
        }
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
        // Coherent films can have dispersive residence times even when the
        // geometric bulk ray is unchanged. Do not report bulk-only transit as
        // the delay of a coated detector, including for a later clipped ray.
        if timed
            && self.surfaces().iter().enumerate().any(|(index, surface)| {
                self.spectral_coating(index).map_or_else(
                    || {
                        surface.coating().is_some_and(|coating| {
                            coating
                                .periodic()
                                .unit_cell()
                                .iter()
                                .any(|layer| layer.thickness_m() != 0.0)
                        })
                    },
                    LensCoating::is_active,
                )
            })
        {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let mut indices = [1.0; Self::MAX_SURFACES + 1];
        let mut extinction = [0.0; Self::MAX_SURFACES + 1];
        let mut group_indices = [1.0; Self::MAX_SURFACES];
        for (i, medium) in std::iter::once(self.sensor_medium())
            .chain(self.surfaces().iter().map(LensSurface::medium_after))
            .enumerate()
        {
            if let Some(glass) = medium {
                indices[i] = glass
                    .ior()
                    .eval(wavelength_nm)
                    .map_err(|_| CameraError::SequentialLensOutsideModel)?;
                extinction[i] = glass
                    .absorption_ref()
                    .extinction_per_m(wavelength_nm)
                    .map_err(|_| CameraError::SequentialLensOutsideModel)?;
                // Only finite sensor-to-surface segments belong to lens
                // transit. The exterior scene's group law is owned by its
                // subsequent native propagation, not this final interface.
                if timed && i < self.surfaces.len() {
                    let group = glass
                        .ior()
                        .group_index(wavelength_nm)
                        .map_err(|_| CameraError::SequentialLensOutsideModel)?;
                    if !group.is_finite() || group <= 0.0 {
                        return Err(CameraError::SequentialLensOutsideModel);
                    }
                    group_indices[i] = group;
                }
            }
        }
        let mut group_delay_s = timed.then_some(Dd::ZERO);
        let mut ray = Ray {
            origin: ray.origin,
            dir: super::scaled_unit(ray.dir)?,
        };
        let axis = reference_axis.unwrap_or_else(|| {
            if ray.dir.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            }
        });
        let mut frame = PolarizationFrame::try_new(ray.dir, axis)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        let (mut distance_m, mut optical_path_m) = (0.0, 0.0);
        for step in PathSteps::new(self.surfaces.len(), ghost) {
            cx.checkpoint()?;
            if (!step.reverse && ray.dir.z <= 0.0) || (step.reverse && ray.dir.z >= 0.0) {
                return Ok(None);
            }
            let i = step.surface;
            let surface = &self.surfaces[i];
            let Some((point, distance)) = surface.intersect(ray)? else {
                return Ok(None);
            };
            let mut normal = surface.normal(point)?;
            if step.reverse {
                normal = normal.scale(-1.0);
            }
            let (incident, transmitted) = if step.reverse { (i + 1, i) } else { (i, i + 1) };
            let cosine = (-normal.dot(ray.dir)).clamp(0.0, 1.0);
            if cosine == 0.0 {
                return Ok(None);
            }
            let response = match self.spectral_coating(i) {
                Some(coating) => coating.fresnel_polarized(
                    wavelength_nm,
                    step.reverse,
                    cosine,
                    indices[incident],
                    indices[transmitted],
                ),
                None => surface
                    .optical_surface
                    .interface(wavelength_nm, step.reverse)
                    .and_then(|interface| {
                        interface.fresnel_polarized(cosine, indices[incident], indices[transmitted])
                    })
                    .map_err(|_| CameraError::SequentialLensOutsideModel),
            }?;
            if response.scalar.total_internal_reflection {
                if !step.reflect {
                    return Ok(None);
                }
                // Steady reflection retains the native evanescent s/p phase.
                // Bulk transit alone cannot certify its dispersive delay.
                if timed {
                    return Err(CameraError::SequentialLensOutsideModel);
                }
            }
            let cross = super::cross(ray.dir, normal);
            let s = if cross.x == 0.0 && cross.y == 0.0 && cross.z == 0.0 {
                frame.horizontal()
            } else {
                super::scaled_unit(cross)?
            };
            let incident_frame = PolarizationFrame::try_new(ray.dir, s)
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            let input = frame.express_in(power, incident_frame);
            let matrix = if step.reflect {
                response.reflection
            } else {
                response.transmission
            };
            power = matrix
                .apply_diagonal_jones(input)
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            if step.reflect
                && input.i > 0.0
                && power.i == 0.0
                && ghost::reflection_has_zero_support(matrix, input)
            {
                // A matched interface or extinguished polarization is a real
                // zero-contribution branch, not Beer/throughput underflow.
                return Ok(None);
            }
            let attenuation = det::exp(-extinction[incident] * distance);
            power = Stokes {
                i: power.i * attenuation,
                q: power.q * attenuation,
                u: power.u * attenuation,
                v: power.v * attenuation,
            };
            let direction = if step.reflect {
                super::scaled_unit(reflect(ray.dir.scale(-1.0), normal))?
            } else {
                refract(
                    ray.dir.scale(-1.0),
                    normal,
                    indices[incident],
                    indices[transmitted],
                    response.scalar.transmitted_cosine,
                )
                .map_err(|_| CameraError::SequentialLensOutsideModel)?
            };
            frame = PolarizationFrame::try_new(direction, s)
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
            ray = Ray {
                origin: point,
                dir: direction,
            };
            distance_m += distance;
            optical_path_m += indices[incident] * distance;
            if let Some(delay) = &mut group_delay_s {
                *delay = *delay + segment_group_delay(group_indices[incident], distance)?;
                if !delay.hi.is_finite() || !delay.lo.is_finite() || delay.hi <= 0.0 {
                    return Err(CameraError::SequentialLensOutsideModel);
                }
            }
            if ![power.i, distance_m, optical_path_m]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0)
            {
                return Err(CameraError::SequentialLensOutsideModel);
            }
        }
        cx.checkpoint()?;
        // Bound only arithmetic drift in the physical cone, using the same
        // primitive as native scene transport. The scalar intensity is retained.
        let power = physical_roundoff(power, power.i)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        Ok(Some(TracedLensBeam {
            ray,
            frame,
            power,
            distance_m,
            optical_path_m,
            sensor_phase_index: indices[0],
            scene_phase_index: indices[self.surfaces.len()],
            group_delay_s,
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
        // Keep the scalar API's original zero-throughput behavior. Reciprocal
        // scene continuation instead refuses an underflowed non-clipped beam.
        Ok(self
            .generate_sequential_beams(cx, x_tan, y_tan, lens_sample, wavelengths_nm, None, false)?
            .map(|beam| {
                beam.map(|beam| SequentialCameraRay {
                    ray: beam.ray(),
                    scene_phase_index: beam.scene_phase_index,
                    weight: beam.geometry_weight * beam.power.i,
                })
            }))
    }

    /// Generate wavelength-specific reciprocal beams for an intensity detector.
    /// Retains lens-induced Q/U/V for subsequent polarizing scene interfaces.
    /// Aperture/TIR misses are `None`, with their original zero sample mass;
    /// finite optical attenuation that underflows to zero is an explicit error.
    pub fn generate_sequential_spectral_beams_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        admit_reciprocal_beams(self.generate_sequential_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            None,
            false,
        )?)
    }

    /// The same reciprocal rays with explicit sensor-to-front group transit.
    /// Each authored finite segment must supply a positive, defined original
    /// group index. Like source-band admission, this is checked before aperture
    /// rejection; the final scene exterior has no lens segment and is excluded.
    /// DD seconds retain each actual segment;
    /// phase optical path is never substituted for missing group information.
    /// No interface residence time or coherent pulse reshaping is modeled;
    /// active films and encountered TIR ghost reflections refuse timing.
    pub fn generate_sequential_spectral_timed_beams_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        admit_reciprocal_beams(self.generate_sequential_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            None,
            true,
        )?)
    }

    /// Place a passive detector row at the physical sensor, before tracing the
    /// lens. The row acts on arriving light along the negative sensor ray; its
    /// Q axis is the supplied world axis projected perpendicular to that ray.
    /// Require `0 < a0 <= 1` and `|(aQ,aU,aV)| <= min(a0,1-a0)`.
    ///
    /// Reciprocity launches `(a0,aQ,-aU,aV)` in the opposite-travel frame.
    /// Lens basis rotations, s/p Fresnel response and Beer loss act on that
    /// complete seed, preserving analyzer/lens/scene polarization correlation.
    /// An isotropic row has no distinguished axis; a nonzero finite reference
    /// is still required, but may be parallel to the ray in that case.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_sequential_spectral_beams_analyzed_from_tangent_offsets<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_row: Stokes,
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        self.generate_sequential_analyzed_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            reference_axis,
            detector_row,
            false,
        )
    }

    /// A sensor analyzer and explicit group transit on the same wavelength-
    /// specific ray. Arguments and reciprocal detector convention match the
    /// steady analyzed generator; only this query requires group derivatives.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets<
        const N: usize,
    >(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_row: Stokes,
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        self.generate_sequential_analyzed_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            reference_axis,
            detector_row,
            true,
        )
    }

    /// Install one absolute passive detector row for each original wavelength
    /// before the real lens. Rows share the supplied world reference axis and
    /// retain their throughput; the original ray, Fresnel, aperture and Beer
    /// owners are unchanged. Callers resolve their source-owned spectral law
    /// before this query; this method never interpolates or borrows support.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets<
        const N: usize,
    >(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_rows: &[Stokes; N],
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        self.generate_sequential_detector_rows(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            reference_axis,
            detector_rows,
            false,
        )
    }

    /// Per-wavelength detector rows with the same qualified sensor-to-front
    /// group delays as the original timed generator. The analyzer contributes
    /// no invented residence time, spectral phase or pulse reshaping.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_sequential_spectral_timed_beams_with_detector_rows_from_tangent_offsets<
        const N: usize,
    >(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_rows: &[Stokes; N],
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        self.generate_sequential_detector_rows(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            reference_axis,
            detector_rows,
            true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn generate_sequential_detector_rows<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_rows: &[Stokes; N],
        timed: bool,
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        cx.checkpoint()?;
        for &row in detector_rows {
            cx.checkpoint()?;
            admit_detector_row(row)?;
        }
        let axis = super::scaled_unit(reference_axis)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        admit_reciprocal_beams(self.generate_sequential_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            Some((axis, *detector_rows)),
            timed,
        )?)
    }

    #[allow(clippy::too_many_arguments)]
    fn generate_sequential_analyzed_beams<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        reference_axis: Vec3,
        detector_row: Stokes,
        timed: bool,
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        cx.checkpoint()?;
        admit_detector_row(detector_row)?;
        let axis = super::scaled_unit(reference_axis)
            .map_err(|_| CameraError::SequentialLensOutsideModel)?;
        admit_reciprocal_beams(self.generate_sequential_beams(
            cx,
            x_tan,
            y_tan,
            lens_sample,
            wavelengths_nm,
            Some((axis, [detector_row; N])),
            timed,
        )?)
    }

    #[allow(clippy::too_many_arguments)]
    fn generate_sequential_beams<const N: usize>(
        &self,
        cx: &Cx<'_>,
        x_tan: f64,
        y_tan: f64,
        lens_sample: LensSample,
        wavelengths_nm: &[f64; N],
        detector: Option<(Vec3, [Stokes; N])>,
        timed: bool,
    ) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
        cx.checkpoint()?;
        if N == 0 || !x_tan.is_finite() || !y_tan.is_finite() {
            return Err(CameraError::InvalidProjection);
        }
        let lens = self
            .sequential_lens
            .as_ref()
            .ok_or(CameraError::InvalidSequentialLens)?;
        let (ghost, lens_sample, path_probability) = match lens.ghost_sampling {
            Some(proposal) => proposal.sample(lens_sample)?,
            None => (None, lens_sample, 1.0),
        };
        let rear = &lens.surfaces()[0];
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
        let geometry = geometry / path_probability;
        if !geometry.is_finite() {
            return Err(CameraError::SequentialLensOutsideModel);
        }
        let detector = detector.map(|(axis, rows)| {
            let local_axis = Vec3::new(
                axis.dot(self.right),
                axis.dot(self.up),
                axis.dot(self.forward),
            );
            rows.map(|row| {
                let axis = if row.q == 0.0 && row.u == 0.0 && row.v == 0.0 {
                    None
                } else {
                    Some(local_axis)
                };
                // (right, up, forward) has determinant -1: world and local
                // right-handed frames share h but reverse v. H*A*row retains
                // both the geometric and reciprocal reversal conventions.
                (axis, Stokes { v: -row.v, ..row })
            })
        });
        let to_world = |v: Vec3| {
            Vec3::new(
                self.right.x * v.x + self.up.x * v.y + self.forward.x * v.z,
                self.right.y * v.x + self.up.y * v.y + self.forward.y * v.z,
                self.right.z * v.x + self.up.z * v.y + self.forward.z * v.z,
            )
        };
        let mut result = [None; N];
        for (lane, (entry, &wavelength)) in result.iter_mut().zip(wavelengths_nm).enumerate() {
            let (axis, power) = detector
                .as_ref()
                .map_or((None, Stokes::unpolarized(1.0)), |seeds| seeds[lane]);
            if let Some(traced) = lens.trace_seeded_path(
                cx,
                Ray {
                    origin: sensor,
                    dir: direction,
                },
                wavelength,
                axis,
                power,
                timed,
                ghost,
            )? {
                if geometry == 0.0 {
                    continue;
                }
                let p = traced.ray.origin;
                let origin = self
                    .eye
                    .offset(self.right.scale(p.x))
                    .offset(self.up.scale(p.y))
                    .offset(self.forward.scale(p.z));
                let frame = PolarizationFrame::try_new(
                    to_world(traced.ray.dir),
                    to_world(traced.frame.horizontal()),
                )
                .map_err(|_| CameraError::SequentialLensOutsideModel)?;
                let power = radiance_power(
                    traced.power,
                    traced.sensor_phase_index,
                    traced.scene_phase_index,
                )?;
                *entry = Some(SequentialCameraBeam {
                    origin,
                    frame,
                    power: Stokes {
                        u: -power.u,
                        v: -power.v,
                        ..power
                    },
                    wavelength_nm: wavelength,
                    geometry_weight: geometry,
                    scene_phase_index: traced.scene_phase_index,
                    group_delay_s: traced.group_delay_s,
                });
            }
        }
        Ok(result)
    }
}

fn admit_detector_row(row: Stokes) -> Result<(), CameraError> {
    if ![row.i, row.q, row.u, row.v].into_iter().all(f64::is_finite)
        || row.i <= 0.0
        || row.i > 1.0
        || det::hypot(det::hypot(row.q, row.u), row.v) > row.i.min(1.0 - row.i)
    {
        return Err(CameraError::SequentialLensOutsideModel);
    }
    Ok(())
}

fn segment_group_delay(group_index: f64, distance_m: f64) -> Result<Dd, CameraError> {
    let group_path = Dd::from_f64(group_index) * Dd::from_f64(distance_m);
    let delay = if group_path.hi.is_finite() && group_path.lo.is_finite() {
        group_path / Dd::from_f64(LIGHT_M_S)
    } else {
        // Retain a finite transit when the product only overflows before /c.
        let (large, small) = if group_index >= distance_m {
            (group_index, distance_m)
        } else {
            (distance_m, group_index)
        };
        (Dd::from_f64(large) / Dd::from_f64(LIGHT_M_S)) * Dd::from_f64(small)
    };
    if !delay.hi.is_finite() || !delay.lo.is_finite() || delay.hi <= 0.0 {
        return Err(CameraError::SequentialLensOutsideModel);
    }
    Ok(delay)
}

fn radiance_power(
    power: Stokes,
    sensor_index: f64,
    scene_index: f64,
) -> Result<Stokes, CameraError> {
    if sensor_index == scene_index {
        return Ok(power);
    }
    // L/n² is the lossless radiance invariant. The native scene subsequently
    // supplies n_scene²/n_source², so the lens factor must occur only here.
    let ratio = sensor_index / scene_index;
    let factor = ratio * ratio;
    if !factor.is_finite() || factor <= 0.0 {
        return Err(CameraError::SequentialLensOutsideModel);
    }
    let scaled = Stokes {
        i: power.i * factor,
        q: power.q * factor,
        u: power.u * factor,
        v: power.v * factor,
    };
    if power.i > 0.0 && scaled.i == 0.0 {
        return Err(CameraError::SequentialLensOutsideModel);
    }
    physical_roundoff(scaled, scaled.i).map_err(|_| CameraError::SequentialLensOutsideModel)
}

fn admit_reciprocal_beams<const N: usize>(
    beams: [Option<SequentialCameraBeam>; N],
) -> Result<[Option<SequentialCameraBeam>; N], CameraError> {
    if beams.iter().flatten().any(|beam| beam.power.i == 0.0) {
        return Err(CameraError::SequentialLensOutsideModel);
    }
    Ok(beams)
}

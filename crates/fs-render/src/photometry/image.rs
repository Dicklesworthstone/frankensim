//! Photopic surface and haze images through the existing physical camera.
//!
//! Pixels contain cd/m², before exposure or tone mapping. The implemented
//! transport consists of directly lit opaque surfaces plus fixture light
//! scattered once on the primary camera ray. It is not a global-illumination
//! solution: interreflection, scattering onto surfaces, and multiple scattering
//! remain unbounded omitted terms. Point emitters have no finite visible area.
//! Pixel centres are sampled once; no antialiasing or pixel-area error is claimed.

use crate::camera::{LensSample, PhysicalCamera};
use crate::charts::Ray;
use crate::tracer::{photometric_brdf, photometric_intersect};
use fs_exec::Cx;

use super::medium::{
    HomogeneousMedium, ScatterError, ScatterEstimate, ScatterSettings, single_scattered_luminance,
};
use super::source::{LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS, photopic_weight_nm, wavelength_nm};
use super::transport::gray_transmittance;
use super::{CompensatedSum, PhotometricScene, PhotometryError, finite_point, scaled_norm, unit};

/// Untonemapped luminance of one ray, with separate transport contributions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LuminanceEstimate {
    /// Surface plus single-scattered luminance, cd/m².
    pub luminance_cd_m2: f64,
    /// Directly lit surface radiance after extinction to the camera, cd/m².
    pub surface_cd_m2: f64,
    /// Haze contribution, including the declared quadrature estimate.
    pub single_scatter: ScatterEstimate,
    /// Unknown total error, including calibration, spectral discretization,
    /// geometry, pixel sampling, interreflection and multiple scattering.
    /// The haze refinement estimate must never be substituted for this bound.
    pub total_error_bound_cd_m2: Option<f64>,
}

/// The direct fixture-lit surface component after extinction to the camera.
/// This quantity excludes light scattered into the camera ray, medium thermal
/// emission and surface interreflection, even when those processes exist in
/// the borrowed physical scene. It is not a total image radiance estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceLuminanceEstimate {
    /// Unscattered directly reflected surface luminance, cd/m².
    pub surface_cd_m2: f64,
    /// Unknown error from source/observer calibration, geometry and spectra.
    pub total_error_bound_cd_m2: Option<f64>,
}

impl PhotometricScene<'_> {
    /// Trace one unit-direction camera ray over a finite declared range.
    /// Primary geometry terminates both visible surface and haze transport.
    /// Existing Lambertian, GGX and spectral conductor BRDFs are evaluated at
    /// the same wavelengths as the source's lumen normalization.
    pub fn luminance_along_ray(
        &self,
        cx: &Cx<'_>,
        ray: &Ray,
        distance_m: f64,
        scatter: ScatterSettings,
    ) -> Result<LuminanceEstimate, PhotometryError> {
        cx.checkpoint()?;
        if let Some(medium) = self.medium.homogeneous() {
            return self.homogeneous_luminance_along_ray(cx, ray, distance_m, scatter, medium);
        }
        if let Some(reason) = self.medium.image_refusal() {
            return Err(PhotometryError::UnsupportedMedium { reason });
        }
        // Passive absorbing grids have analytically zero haze and thermal
        // source terms. Their coloured attenuation still stays inside the
        // spectral BRDF integral on both physical connection legs.
        let single_scatter = single_scattered_luminance(
            &HomogeneousMedium::vacuum(),
            ray,
            distance_m,
            scatter,
            cx,
            |_, _| Ok(0.0),
        )
        .map_err(scatter_error)?;
        let surface = self.direct_surface_luminance_along_ray(cx, ray, distance_m)?;
        Ok(LuminanceEstimate {
            luminance_cd_m2: surface.surface_cd_m2,
            surface_cd_m2: surface.surface_cd_m2,
            single_scatter,
            total_error_bound_cd_m2: None,
        })
    }

    fn homogeneous_luminance_along_ray(
        &self,
        cx: &Cx<'_>,
        ray: &Ray,
        distance_m: f64,
        scatter: ScatterSettings,
        medium: HomogeneousMedium,
    ) -> Result<LuminanceEstimate, PhotometryError> {
        cx.checkpoint()?;
        validate_camera_ray(self, ray, distance_m)?;
        let primary = photometric_intersect(self.primitives, cx, ray, distance_m)?;
        let limit = primary.as_ref().map_or(distance_m, |(_, hit)| hit.t);
        if !limit.is_finite() || limit < 0.0 || limit > distance_m {
            return Err(PhotometryError::NonFiniteResult);
        }
        // Point-source irradiance is singular at the source. Midpoint samples
        // can miss that singularity entirely, so admit the whole visible
        // camera segment before estimating its in-scattered radiance.
        if medium.scattering_per_m() > 0.0 && limit > 0.0 {
            for (index, source) in self.sources.iter().enumerate() {
                cx.checkpoint()?;
                if source.flux_lm() == 0.0 {
                    continue;
                }
                let offset = source.position_m().delta_from(ray.origin);
                let projection = offset.dot(ray.dir) / ray.dir.dot(ray.dir);
                if !projection.is_finite() {
                    return Err(PhotometryError::NonFiniteResult);
                }
                let closest = ray.at(projection.clamp(0.0, limit));
                let separation = source.position_m().delta_from(closest);
                let distance = if separation.x == 0.0 && separation.y == 0.0 && separation.z == 0.0
                {
                    0.0
                } else {
                    scaled_norm(separation).ok_or(PhotometryError::NonFiniteResult)?
                };
                if distance <= self.options.ray_epsilon_m {
                    return Err(PhotometryError::SourceTooClose { source: index });
                }
            }
        }
        let mut surface_sum = CompensatedSum::default();
        if let Some((primitive_index, hit)) = primary {
            let outgoing = ray.dir.scale(-1.0);
            let mut normal =
                hit.normal
                    .and_then(unit)
                    .ok_or(PhotometryError::InvalidPrimitive {
                        primitive: primitive_index,
                    })?;
            // Same two-sided opaque convention as the existing spectral tracer.
            if normal.dot(outgoing) < 0.0 {
                normal = normal.scale(-1.0);
            }
            let material = &self.primitives[primitive_index].material;
            for (index, source) in self.sources.iter().enumerate() {
                cx.checkpoint()?;
                if source.flux_lm() == 0.0 {
                    continue;
                }
                let connection = self.connection(cx, source, index, hit.point)?;
                let cosine = normal.dot(connection.towards_source).max(0.0);
                if cosine == 0.0
                    || connection.illuminance_normal_lux == 0.0
                    || !self.visible(cx, source.position_m(), &connection)?
                {
                    continue;
                }
                let bins = source.spectrum().bins_w_per_lumen_nm();
                let mut response = CompensatedSum::default();
                for (bin, power_per_lumen_nm) in bins.iter().copied().enumerate() {
                    if bin % 16 == 0 {
                        cx.checkpoint()?;
                    }
                    let brdf = photometric_brdf(
                        material,
                        normal,
                        outgoing,
                        connection.towards_source,
                        wavelength_nm(bin),
                    )?;
                    if !brdf.is_finite() || brdf < 0.0 {
                        return Err(PhotometryError::NonFiniteResult);
                    }
                    response.add(power_per_lumen_nm * photopic_weight_nm(bin) * brdf);
                }
                debug_assert_eq!(bins.len(), SPECTRAL_BINS);
                let reflected = connection.illuminance_normal_lux
                    * cosine
                    * LUMINOUS_EFFICACY_LM_W
                    * response.finish()?;
                surface_sum.add(reflected);
            }
        }
        let surface_cd_m2 = surface_sum.finish()? * medium.transmittance(limit)?;
        let single_scatter = single_scattered_luminance(
            &medium,
            ray,
            limit,
            scatter,
            cx,
            |point, towards_camera| {
                let mut incident = CompensatedSum::default();
                for (index, source) in self.sources.iter().enumerate() {
                    cx.checkpoint()?;
                    if source.flux_lm() == 0.0 {
                        continue;
                    }
                    let connection = self.connection(cx, source, index, point)?;
                    if connection.illuminance_normal_lux == 0.0
                        || !self.visible(cx, source.position_m(), &connection)?
                    {
                        continue;
                    }
                    // Both vectors are photon propagation directions. The
                    // source connection stores the opposite (toward-light) one.
                    let cosine = -connection.towards_source.dot(towards_camera);
                    incident.add(connection.illuminance_normal_lux * medium.phase_per_sr(cosine)?);
                }
                incident.finish()
            },
        )
        .map_err(scatter_error)?;
        let luminance_cd_m2 = surface_cd_m2 + single_scatter.luminance_cd_m2;
        if !luminance_cd_m2.is_finite() || luminance_cd_m2 < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        cx.checkpoint()?;
        Ok(LuminanceEstimate {
            luminance_cd_m2,
            surface_cd_m2,
            single_scatter,
            total_error_bound_cd_m2: None,
        })
    }

    /// Evaluate only unscattered fixture light reflected once by the first
    /// opaque surface. Physical grid absorption and scattering both attenuate
    /// each wavelength on both legs; grid emission and in-scattered light are
    /// intentionally outside this explicitly named component.
    pub fn direct_surface_luminance_along_ray(
        &self,
        cx: &Cx<'_>,
        ray: &Ray,
        distance_m: f64,
    ) -> Result<SurfaceLuminanceEstimate, PhotometryError> {
        cx.checkpoint()?;
        validate_camera_ray(self, ray, distance_m)?;
        let primary = photometric_intersect(self.primitives, cx, ray, distance_m)?;
        let Some((primitive_index, hit)) = primary else {
            return Ok(SurfaceLuminanceEstimate {
                surface_cd_m2: 0.0,
                total_error_bound_cd_m2: None,
            });
        };
        if !hit.t.is_finite() || hit.t < 0.0 || hit.t > distance_m {
            return Err(PhotometryError::NonFiniteResult);
        }
        let outgoing = ray.dir.scale(-1.0);
        let mut normal = hit
            .normal
            .and_then(unit)
            .ok_or(PhotometryError::InvalidPrimitive {
                primitive: primitive_index,
            })?;
        if normal.dot(outgoing) < 0.0 {
            normal = normal.scale(-1.0);
        }
        let camera_transmission = self.medium.transmittance(cx, ray.origin, ray.dir, hit.t)?;
        let camera_gray = gray_transmittance(&camera_transmission);
        let material = &self.primitives[primitive_index].material;
        let mut total = CompensatedSum::default();
        for (index, source) in self.sources.iter().enumerate() {
            cx.checkpoint()?;
            if source.flux_lm() == 0.0 {
                continue;
            }
            let connection = self.unattenuated_connection(source, index, hit.point)?;
            let cosine = normal.dot(connection.towards_source).max(0.0);
            if cosine == 0.0
                || connection.illuminance_normal_lux == 0.0
                || !self.visible(cx, source.position_m(), &connection)?
            {
                continue;
            }
            let light_transmission = self.medium.transmittance(
                cx,
                source.position_m(),
                connection.towards_source.scale(-1.0),
                connection.distance_m,
            )?;
            let light_gray = gray_transmittance(&light_transmission);
            let mut response = CompensatedSum::default();
            for (bin, power) in source
                .spectrum()
                .bins_w_per_lumen_nm()
                .iter()
                .copied()
                .enumerate()
            {
                if bin % 16 == 0 {
                    cx.checkpoint()?;
                }
                let brdf = photometric_brdf(
                    material,
                    normal,
                    outgoing,
                    connection.towards_source,
                    wavelength_nm(bin),
                )?;
                if !brdf.is_finite() || brdf < 0.0 {
                    return Err(PhotometryError::NonFiniteResult);
                }
                let mut contribution = power * photopic_weight_nm(bin) * brdf;
                if light_gray.is_none() {
                    contribution *= light_transmission[bin];
                }
                if camera_gray.is_none() {
                    contribution *= camera_transmission[bin];
                }
                response.add(contribution);
            }
            // Gray legs retain the original multiplication order exactly.
            let incident = connection.illuminance_normal_lux * light_gray.unwrap_or(1.0);
            total.add(incident * cosine * LUMINOUS_EFFICACY_LM_W * response.finish()?);
        }
        let surface_cd_m2 = total.finish()? * camera_gray.unwrap_or(1.0);
        if !surface_cd_m2.is_finite() || surface_cd_m2 < 0.0 {
            return Err(PhotometryError::NonFiniteResult);
        }
        cx.checkpoint()?;
        Ok(SurfaceLuminanceEstimate {
            surface_cd_m2,
            total_error_bound_cd_m2: None,
        })
    }
}

fn scatter_error(error: ScatterError<PhotometryError>) -> PhotometryError {
    match error {
        ScatterError::Cancelled => PhotometryError::Cancelled,
        ScatterError::Incident(error) => error,
        ScatterError::Medium(error) => PhotometryError::Medium(error),
        ScatterError::Settings => PhotometryError::InvalidScatterSettings,
        ScatterError::Ray => PhotometryError::InvalidCameraRay,
        ScatterError::IncidentValue | ScatterError::Arithmetic => PhotometryError::NonFiniteResult,
    }
}

fn validate_distance(scene: &PhotometricScene<'_>, distance_m: f64) -> Result<(), PhotometryError> {
    if !distance_m.is_finite() || distance_m <= scene.options.ray_epsilon_m {
        return Err(PhotometryError::InvalidCameraRay);
    }
    if distance_m > scene.options.max_distance_m {
        return Err(PhotometryError::DistanceLimitExceeded);
    }
    Ok(())
}

fn validate_camera_ray(
    scene: &PhotometricScene<'_>,
    ray: &Ray,
    distance_m: f64,
) -> Result<(), PhotometryError> {
    validate_distance(scene, distance_m)?;
    if !finite_point(ray.origin)
        || !ray.dir.x.is_finite()
        || !ray.dir.y.is_finite()
        || !ray.dir.z.is_finite()
        || (ray.dir.norm() - 1.0).abs() > 1.0e-12
        || !finite_point(ray.at(distance_m))
    {
        return Err(PhotometryError::InvalidCameraRay);
    }
    Ok(())
}

/// Progressive row-major luminance image, bound to one immutable scene and
/// physical camera. It returns ordinary floating-point data with unit-bearing
/// field names. Exposure metadata is not applied to a luminance measurement.
pub struct LuminanceJob<'a, 'scene> {
    scene: &'a PhotometricScene<'scene>,
    camera: &'a PhysicalCamera,
    width: u32,
    height: u32,
    count: usize,
    distance_m: f64,
    scatter: ScatterSettings,
    values: Vec<LuminanceEstimate>,
}

impl<'a, 'scene> LuminanceJob<'a, 'scene> {
    /// Begin a centre-sampled pinhole image. Finite apertures explicitly refuse
    /// until a lens-area estimator and its uncertainty are supplied.
    pub fn new(
        scene: &'a PhotometricScene<'scene>,
        camera: &'a PhysicalCamera,
        width: u32,
        height: u32,
        distance_m: f64,
        scatter: ScatterSettings,
    ) -> Result<Self, PhotometryError> {
        validate_distance(scene, distance_m)?;
        if let Some(reason) = scene.medium.image_refusal() {
            return Err(PhotometryError::UnsupportedMedium { reason });
        }
        if width == 0 || height == 0 || !camera.aperture().is_pinhole() {
            return Err(PhotometryError::InvalidCameraRay);
        }
        let count = usize::try_from(u64::from(width) * u64::from(height))
            .map_err(|_| PhotometryError::OutputBudgetExceeded)?;
        if count > scene.options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(count)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        Ok(Self {
            scene,
            camera,
            width,
            height,
            count,
            distance_m,
            scatter,
            values,
        })
    }

    /// Append at most `max_pixels` whole pixels. Cancellation or a failed pixel
    /// preserves all previously accepted pixels and does not append that pixel.
    pub fn advance(&mut self, cx: &Cx<'_>, max_pixels: usize) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.values.len();
        let end = start.saturating_add(max_pixels).min(self.count);
        let half_tan = self.camera.projection().vertical_half_tan();
        let aspect = f64::from(self.width) / f64::from(self.height);
        for index in start..end {
            let x = (index % self.width as usize) as f64 + 0.5;
            let y = (index / self.width as usize) as f64 + 0.5;
            let x_tan = (2.0 * x / f64::from(self.width) - 1.0) * half_tan * aspect;
            let y_tan = (1.0 - 2.0 * y / f64::from(self.height)) * half_tan;
            let ray = self
                .camera
                .generate_ray_from_tangent_offsets(cx, x_tan, y_tan, LensSample::CENTER)
                .map_err(PhotometryError::from)?;
            let value = self
                .scene
                .luminance_along_ray(cx, &ray, self.distance_m, self.scatter)?;
            cx.checkpoint()?;
            self.values.push(value);
        }
        Ok(self.values.len() - start)
    }

    /// Successfully completed row-major prefix, with cd/m² and model metadata.
    #[must_use]
    pub fn results(&self) -> &[LuminanceEstimate] {
        &self.values
    }

    /// Requested image dimensions.
    #[must_use]
    pub const fn dimensions(&self) -> (u32, u32) {
        (self.width, self.height)
    }

    /// Whether the complete image is available.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.values.len() == self.count
    }
}

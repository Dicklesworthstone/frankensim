//! Native polarized environment next-event estimation and complementary MIS.
//!
//! One original map-direction proposal competes with the original continuous
//! receiver proposal. Its direct segment must stay in a declared homogeneous
//! nonabsorbing medium and encounter no other surface. Delta chains and GRIN
//! angular maps are not represented by this proposal; their camera paths keep
//! weight one. The native event owner proves visibility without a ray offset.
//! See PBRT 4e, section 13.4:
//! <https://pbr-book.org/4ed/Light_Transport_I_Surface_Reflection/A_Better_Path_Tracer>.

use super::super::area_lights::camera::{AreaLightCameraError, SurfaceReceiver};
use super::super::emission::SurfaceEmission;
use super::super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedCameraScene,
    EnvironmentRadiance, RadianceEndpoint, finite_nonnegative, phase_index,
};
use crate::charts::Ray;
use crate::gradient_index::event::GradientEventOptions;
use crate::polarization::dielectric::scene::curved::{
    CurvedPathResult, CurvedSceneError, CurvedTraceEvent, CurvedTraceSettings,
};
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

mod microsurface;
mod source;
mod visibility;
mod weights;

/// Environment connections on native molecular and retained particle walks.
#[cfg(feature = "volumes")]
pub mod participating;

use visibility::{clear_to_environment, supports_segment};
use weights::{positive_product, power_weight};

pub use microsurface::{SmithDielectricEnvironmentCameraScene, SmithEnvironmentCameraScene};

/// A fixed-direction physical response and deterministic comparison density.
/// For ordinary surfaces the latter is their exact proposal. The Smith view
/// supplies only the analytic first-reflection escape subdensity instead.
struct EnvironmentConnection {
    normal: fs_geom::Vec3,
    factor_per_sr: f64,
    medium: usize,
    comparison_density: f64,
}

/// Original camera limits and a separate bounded first-boundary shadow walk.
/// A distance limit is successful only when the original scene bounds also
/// prove that the remaining infinite halfline has no contact.
#[derive(Clone, Copy, Debug)]
pub struct EnvironmentCameraSettings {
    pub surfaces: CurvedTraceSettings,
    pub visibility: GradientEventOptions,
    pub visibility_distance_m: f64,
    pub max_environment_samples: u32,
    /// Cumulative geometry work across all environment shadow queries.
    pub max_query_work: u64,
}

impl Default for EnvironmentCameraSettings {
    fn default() -> Self {
        Self {
            surfaces: CurvedTraceSettings::default(),
            visibility: GradientEventOptions::default(),
            visibility_distance_m: 100.0,
            max_environment_samples: 64,
            max_query_work: 1_000_000,
        }
    }
}

impl EnvironmentCameraSettings {
    fn validate(self) -> Result<(), EnvironmentCameraError> {
        if !self.visibility_distance_m.is_finite()
            || self.visibility_distance_m <= 0.0
            || self.max_environment_samples == 0
            || self.max_query_work == 0
        {
            return Err(EnvironmentCameraError::InvalidSettings);
        }
        self.visibility
            .validate()
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum EnvironmentCameraError {
    Camera(CameraRadianceError),
    Receiver(AreaLightCameraError),
    MissingEnvironment,
    InvalidSettings,
    InvalidSample,
    SampleLimit,
    GeometryLimit,
    UnresolvedVisibility,
    NumericalRange,
}

impl From<CameraRadianceError> for EnvironmentCameraError {
    fn from(value: CameraRadianceError) -> Self {
        Self::Camera(value)
    }
}
impl From<AreaLightCameraError> for EnvironmentCameraError {
    fn from(value: AreaLightCameraError) -> Self {
        Self::Receiver(value)
    }
}
impl From<fs_exec::Cancelled> for EnvironmentCameraError {
    fn from(value: fs_exec::Cancelled) -> Self {
        Self::Camera(value.into())
    }
}
impl core::fmt::Display for EnvironmentCameraError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "native environment camera: {self:?}")
    }
}
impl core::error::Error for EnvironmentCameraError {}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Work {
    samples: u32,
    shadows: u32,
    fields: u64,
    geometry: u64,
}

/// The original native camera path plus the added direct environment term.
/// Transport weight, contacts and continuation draws retain their original
/// meaning. The endpoint weight partitions only escaped environment radiance;
/// already encountered surface and thermal emission are never scaled by it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnvironmentCameraSample {
    radiance: CameraRadianceSample,
    direct_environment: f64,
    endpoint_weight: f64,
    work: Work,
}

impl EnvironmentCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        &self.radiance
    }
    pub const fn environment_radiance_w_m2_sr_nm(&self) -> f64 {
        self.direct_environment
    }
    pub const fn environment_endpoint_weight(&self) -> f64 {
        self.endpoint_weight
    }
    pub const fn environment_samples(&self) -> u32 {
        self.work.samples
    }
    pub const fn visibility_queries(&self) -> u32 {
        self.work.shadows
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.work.fields
    }
    pub const fn query_work(&self) -> u64 {
        self.work.geometry
    }
}

/// An allocation-free borrowed view over the exact native camera and its
/// calibrated environment. Smooth delta receivers retain camera transport;
/// diffuse, coated, rough and brushed receivers share their original physical
/// polarized BSDF evaluator and solid-angle proposal.
pub struct CurvedEnvironmentCameraScene<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    environment: EnvironmentRadiance<'s>,
}

impl<'a, 's, 'm> CurvedEnvironmentCameraScene<'a, 's, 'm> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
    ) -> Result<Self, EnvironmentCameraError> {
        cx.checkpoint()?;
        let environment = camera
            .environment
            .ok_or(EnvironmentCameraError::MissingEnvironment)?;
        Ok(Self {
            camera,
            environment,
        })
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }

    /// Preserve the original camera Philox domain. Environment directions use
    /// their own domain keyed by the same native surface-event ordinal.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_ray_with_samples(
            cx,
            ray,
            wavelength_nm,
            initial_medium,
            settings,
            |event| counter_samples(seed, path_id, event, 0x706f_6c32),
            |event| {
                let samples = counter_samples(seed, path_id, event, 0x706f_6c39);
                [samples[0], samples[1]]
            },
        )
    }

    /// One environment sample at each nonblack continuous receiver. The two
    /// sample streams remain independent and each retains its own ordinal;
    /// no environment query resamples the camera's continuation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_with_launch(
            cx,
            || self.camera.launch(ray, wavelength_nm, initial_medium),
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
        )
    }

    /// Use the original passive analyzer launch before any camera or source
    /// interaction. Circular handedness and reciprocal frames stay with the
    /// native polarized transport owner, including rough receiver response.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_with_launch(
            cx,
            || analyzer.launch(self.camera, ray, wavelength_nm, initial_medium),
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_with_launch(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        settings.validate()?;
        let beam = launch()?;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut work = Work::default();
        let mut direct = 0.0;
        let mut continuous = false;
        let mut continuation_pdf = None;
        let mut source_error = None;
        let active = self.environment.scale_w_m2_sr_nm > 0.0 && !self.environment.map.is_black();
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings.surfaces,
            surface_samples,
            true,
            &self.camera.spatial_volume_emitters,
            |event| {
                let observed = match event {
                    CurvedTraceEvent::Arrival(path) => {
                        emission.observe(path)?;
                        (|| {
                            let receiver = SurfaceReceiver::at_arrival(self.camera, path)?;
                            continuous = active && receiver.is_some();
                            if let Some(receiver) = receiver.filter(|_| active) {
                                let value = self.add_surface(
                                    cx,
                                    receiver,
                                    path,
                                    initial_index,
                                    settings,
                                    &mut work,
                                    &mut environment_samples,
                                )?;
                                direct = finite_nonnegative(direct + value)?;
                            }
                            Ok(())
                        })()
                    }
                    CurvedTraceEvent::Departure {
                        path, pdf, delta, ..
                    } => {
                        emission.depart(path);
                        (|| {
                            continuation_pdf = None;
                            if continuous
                                && !delta
                                && supports_segment(
                                    self.camera.scene,
                                    path.medium_index(),
                                    path.beam().wavelength_nm(),
                                )?
                            {
                                if !pdf.is_finite() || pdf <= 0.0 {
                                    return Err(EnvironmentCameraError::NumericalRange);
                                }
                                continuation_pdf = Some(pdf);
                            }
                            Ok(())
                        })()
                    }
                };
                if let Err(error) = observed {
                    source_error = Some(error);
                    Err(CurvedSceneError::NumericalRange)
                } else {
                    Ok(())
                }
            },
        );
        if let Some(error) = source_error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        // Ask the original owner for only its endpoint term, preserving its
        // infinite-absorption rule, native analyzer and source contraction.
        let mut radiance = self.camera.finish(cx, path, initial_index, 0.0)?;
        let endpoint_weight = if radiance.endpoint() == RadianceEndpoint::Environment {
            if let Some(pdf) = continuation_pdf {
                let source = self
                    .environment
                    .map
                    .evaluate(path.beam().ray().dir)
                    .ok_or(EnvironmentCameraError::NumericalRange)?;
                power_weight(pdf, source.pdf_solid_angle)?
            } else {
                1.0
            }
        } else {
            1.0
        };
        let endpoint = positive_product(&[radiance.radiance_w_m2_sr_nm(), endpoint_weight], 1.0)?;
        radiance.radiance_w_m2_sr_nm = finite_nonnegative(emitted + direct + endpoint)?;
        cx.checkpoint()?;
        Ok(EnvironmentCameraSample {
            radiance,
            direct_environment: direct,
            endpoint_weight,
            work,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn add_surface(
        &self,
        cx: &Cx<'_>,
        receiver: SurfaceReceiver,
        path: &CurvedPathResult,
        initial_index: f64,
        settings: EnvironmentCameraSettings,
        work: &mut Work,
        samples: &mut impl FnMut(u32) -> [f64; 2],
    ) -> Result<f64, EnvironmentCameraError> {
        self.add_surface_with_response(cx, path, initial_index, settings, work, samples, |source| {
            let (response, medium) =
                receiver.directional_response(self.camera, source.frame, source.stokes)?;
            Ok(EnvironmentConnection {
                normal: response.normal,
                factor_per_sr: if receiver.is_rough() {
                    response.factor
                } else {
                    response.factor / core::f64::consts::PI
                },
                medium,
                comparison_density: response.pdf,
            })
        })
    }

    /// One shared calibrated proposal, native visibility query and physical
    /// radiometric accumulation. The response owner supplies a deterministic
    /// angular partition guide separately from its optical response.
    #[allow(clippy::too_many_arguments)]
    fn add_surface_with_response(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        initial_index: f64,
        settings: EnvironmentCameraSettings,
        work: &mut Work,
        samples: &mut impl FnMut(u32) -> [f64; 2],
        response: impl FnOnce(
            source::EnvironmentSourceSample,
        ) -> Result<EnvironmentConnection, EnvironmentCameraError>,
    ) -> Result<f64, EnvironmentCameraError> {
        if work.samples >= settings.max_environment_samples {
            return Err(EnvironmentCameraError::SampleLimit);
        }
        work.samples += 1;
        let ordinal = path
            .surface_events()
            .checked_sub(1)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        let Some(source) = self
            .environment
            .sample_direction(path.beam().wavelength_nm(), samples(ordinal))?
        else {
            return Ok(0.0);
        };
        let response = response(source)?;
        if response.factor_per_sr == 0.0
            || !supports_segment(
                self.camera.scene,
                response.medium,
                path.beam().wavelength_nm(),
            )?
        {
            return Ok(0.0);
        }
        let cosine = response.normal.dot(source.direction).abs().min(1.0);
        if cosine == 0.0 {
            return Ok(0.0);
        }
        let contact = *path
            .last_contact()
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        work.shadows = work
            .shadows
            .checked_add(1)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        let remaining = settings
            .max_query_work
            .checked_sub(work.geometry)
            .ok_or(EnvironmentCameraError::GeometryLimit)?;
        let visibility = clear_to_environment(
            cx,
            self.camera.scene,
            Some(contact),
            Ray {
                origin: path.beam().ray().origin,
                dir: source.direction,
            },
            source.frame.horizontal(),
            path.beam().wavelength_nm(),
            response.medium,
            settings,
            remaining,
        )?;
        work.fields = work
            .fields
            .checked_add(visibility.fields)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        work.geometry = work
            .geometry
            .checked_add(visibility.geometry)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        if !visibility.visible {
            return Ok(0.0);
        }
        let weight = power_weight(source.pdf, response.comparison_density)?;
        let local_index = phase_index(path.beam().medium(), path.beam().wavelength_nm())?;
        let ratio = initial_index / local_index;
        positive_product(
            &[
                path.beam().power().i,
                ratio,
                ratio,
                response.factor_per_sr,
                source.radiance,
                cosine,
                weight,
            ],
            source.pdf,
        )
    }
}

fn counter_samples(seed: u64, path_id: u64, event: u32, domain: u32) -> [f64; 3] {
    let block = philox4x32_10(
        [path_id as u32, (path_id >> 32) as u32, event, domain],
        [seed as u32, (seed >> 32) as u32],
    );
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}

#[cfg(test)]
mod tests;

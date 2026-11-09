//! Full Smith environment connections with a deterministic angular partition.
//!
//! The complete polarized BRDF is stochastic and its outgoing marginal PDF
//! is unavailable. Instead use the analytic ONE-reflection, unit-albedo escape
//! subdensity q1 and the environment proposal p only as deterministic guides:
//! w_light = p²/(p²+q1²), w_camera = q1²/(p²+q1²).
//! These are complementary functions of the macrostate and exit direction;
//! neither contains a sampled BRDF, hidden collision density or random facet.
//! Linearity gives integral(f*L*w_light) + integral(f*L*w_camera) = integral(f*L).
//! q1 is NOT substituted for the continuation PDF in its throughput.
//! This retains useful narrow-lobe continuation without a fitted MS albedo.
//!
//! Only direct segments in declared homogeneous, nonabsorbing media compete.
//! Intervening surfaces, delta chains and unsupported angular maps reset the
//! endpoint partition. Surface/LTE emission is never scaled by this weight.
//! See Heitz et al. (2016), and PBRT 4e §2.2.3 for partition conditions.

use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::{CurvedTransportEvent, SmithSurfaceSampler};
use crate::polarization::dielectric::scene::smith_draw;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;

mod dielectric;
pub use dielectric::SmithDielectricEnvironmentCameraScene;

/// Borrowed opt-in view: opaque rough and brushed receivers use complete
/// Smith reflection for both continuation and source connections. Other
/// materials keep the original native law and its original exact-PDF MIS.
/// Each light estimate has its own bounded, independent microfacet stream.
pub struct SmithEnvironmentCameraScene<'v, 'a, 's, 'm> {
    source: &'v CurvedEnvironmentCameraScene<'a, 's, 'm>,
    microsurface: SmithWalkSettings,
}

#[derive(Clone, Copy)]
enum EndpointPartition {
    Original,
    Ordinary(f64),
    Smith(SurfaceReceiver),
    Dielectric(SurfaceReceiver),
}

impl EndpointPartition {
    fn weight(
        self,
        source: &CurvedEnvironmentCameraScene<'_, '_, '_>,
        path: &CurvedPathResult,
    ) -> Result<f64, EnvironmentCameraError> {
        let density = match self {
            Self::Original => return Ok(1.0),
            Self::Ordinary(pdf) => pdf,
            Self::Smith(receiver) => {
                receiver.smith_first_reflection_density(source.camera, path.beam().ray().dir)?
            }
            Self::Dielectric(receiver) => {
                dielectric::angular_guide(receiver, source.camera, path.beam().ray().dir)?
            }
        };
        let environment = source
            .environment
            .map
            .evaluate(path.beam().ray().dir)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        // If the analytic subset has no support, the light estimator owns
        // that direction. An unsupported light direction stays camera-owned.
        if density == 0.0 {
            return Ok(if environment.pdf_solid_angle == 0.0 {
                1.0
            } else {
                0.0
            });
        }
        power_weight(density, environment.pdf_solid_angle)
    }
}

impl<'v, 'a, 's, 'm> SmithEnvironmentCameraScene<'v, 'a, 's, 'm> {
    pub const fn new(
        source: &'v CurvedEnvironmentCameraScene<'a, 's, 'm>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            source,
            microsurface,
        }
    }

    pub const fn source_environment_camera(&self) -> &'v CurvedEnvironmentCameraScene<'a, 's, 'm> {
        self.source
    }

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
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.source.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Preserve the resolved lens's reciprocal response; its caller still
    /// owns the original pupil/geometric raster weight.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| counter_samples(seed, path_id, event, 0x706f_6c32),
            |event| {
                let sample = counter_samples(seed, path_id, event, 0x706f_6c39);
                [sample[0], sample[1]]
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [(seed as u32) ^ 0x736d_656e, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
        )
    }

    /// Independent ordinary-surface, environment-direction, continuing-facet
    /// and fixed-direction-evaluation streams. All surface and facet ordinals
    /// use the original zero-based convention. Each facet stream remains
    /// independently indexed by (surface ordinal, facet ordinal).
    /// A refusal discards the whole estimate; averaging only successful
    /// histories after a facet/work limit would bias the image.
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
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
        )
    }

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
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.source.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
            None,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_microsurfaces(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut environment_samples: impl FnMut(u32) -> [f64; 2],
        mut facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut connection_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
        mut dielectric_connection_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        settings.validate()?;
        let camera = self.source.camera;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, initial_index, beam);
        let mut work = Work::default();
        let mut direct = 0.0;
        let mut receiver = None;
        let mut partition = EndpointPartition::Original;
        let mut source_error = None;
        let include_dielectric = dielectric_samples.is_some();
        let active = self.source.environment.scale_w_m2_sr_nm > 0.0
            && !self.source.environment.map.is_black();
        let result = camera
            .scene
            .trace_with_flights_events_and_smith(
                cx,
                beam,
                initial_medium,
                settings.surfaces,
                surface_samples,
                true,
                &camera.spatial_volume_emitters,
                |event| {
                    let observed = match event {
                        CurvedTransportEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            emission.observe(path)?;
                            (|| {
                                // A later surface ends the earlier direct segment.
                                // Only this newly arriving vertex may compete next.
                                partition = EndpointPartition::Original;
                                receiver =
                                    SurfaceReceiver::at_arrival(camera, path)?.filter(|_| active);
                                if let Some(receiving) = receiver {
                                    let value = if include_dielectric
                                        && receiving.uses_smith_dielectric(camera)?
                                    {
                                        // A rough equal-index interface can mix a
                                        // straight-through atom and a continuous
                                        // reflected law. Keep its complete native
                                        // walk and endpoint, without a mismatched
                                        // legacy single-scatter source estimate.
                                        if !receiving
                                            .supports_smith_dielectric_connections(camera)?
                                            || !dielectric::any_supported_medium(
                                                receiving,
                                                camera,
                                                path.beam().wavelength_nm(),
                                            )?
                                        {
                                            return Ok(());
                                        }
                                        self.source.add_surface_with_response(
                                            cx,
                                            path,
                                            initial_index,
                                            settings,
                                            &mut work,
                                            &mut environment_samples,
                                            |source| {
                                                let medium = receiving.directional_source_medium(
                                                    camera,
                                                    source.frame.direction(),
                                                )?;
                                                if !supports_segment(
                                                    camera.scene,
                                                    medium,
                                                    path.beam().wavelength_nm(),
                                                )? {
                                                    return Ok(EnvironmentConnection {
                                                        normal: source.direction,
                                                        factor_per_sr: 0.0,
                                                        medium,
                                                        comparison_density: 0.0,
                                                    });
                                                }
                                                let guide = dielectric::angular_guide(
                                                    receiving,
                                                    camera,
                                                    source.direction,
                                                )?;
                                                let ordinal =
                                                    path.surface_events().checked_sub(1).ok_or(
                                                        EnvironmentCameraError::NumericalRange,
                                                    )?;
                                                let draws =
                                                    dielectric_connection_samples.as_mut().ok_or(
                                                        EnvironmentCameraError::NumericalRange,
                                                    )?;
                                                let (normal, factor_per_sr, actual_medium) =
                                                    receiving
                                                        .directional_smith_dielectric_response(
                                                            cx,
                                                            camera,
                                                            source.frame,
                                                            source.stokes,
                                                            self.microsurface,
                                                            |facet| draws(ordinal, facet),
                                                        )?;
                                                if actual_medium != medium {
                                                    return Err(
                                                        EnvironmentCameraError::NumericalRange,
                                                    );
                                                }
                                                Ok(EnvironmentConnection {
                                                    normal,
                                                    factor_per_sr,
                                                    medium,
                                                    comparison_density: guide,
                                                })
                                            },
                                        )?
                                    } else if receiving.uses_smith(camera)? {
                                        // An opaque receiver's source leg has
                                        // this exact adjacent medium. Unsupported
                                        // angular maps stay continuation-owned
                                        // without consuming light/facet samples.
                                        if !supports_segment(
                                            camera.scene,
                                            path.medium_index(),
                                            path.beam().wavelength_nm(),
                                        )? {
                                            return Ok(());
                                        }
                                        self.source.add_surface_with_response(
                                            cx,
                                            path,
                                            initial_index,
                                            settings,
                                            &mut work,
                                            &mut environment_samples,
                                            |source| {
                                                let ordinal =
                                                    path.surface_events().checked_sub(1).ok_or(
                                                        EnvironmentCameraError::NumericalRange,
                                                    )?;
                                                let (normal, response, medium) = receiving
                                                    .directional_smith_response(
                                                        cx,
                                                        camera,
                                                        source.frame,
                                                        source.stokes,
                                                        self.microsurface,
                                                        |facet| connection_samples(ordinal, facet),
                                                    )?;
                                                Ok(EnvironmentConnection {
                                                    normal,
                                                    factor_per_sr: response.factor_per_sr(),
                                                    medium,
                                                    comparison_density: response
                                                        .first_reflection_density(),
                                                })
                                            },
                                        )?
                                    } else {
                                        self.source.add_surface(
                                            cx,
                                            receiving,
                                            path,
                                            initial_index,
                                            settings,
                                            &mut work,
                                            &mut environment_samples,
                                        )?
                                    };
                                    direct = finite_nonnegative(direct + value)?;
                                }
                                Ok(())
                            })()
                        }
                        CurvedTransportEvent::Surface(CurvedTraceEvent::Departure {
                            path,
                            pdf,
                            delta,
                            ..
                        }) => {
                            emission.depart(path);
                            (|| {
                                partition = EndpointPartition::Original;
                                if receiver.is_some()
                                    && !delta
                                    && supports_segment(
                                        camera.scene,
                                        path.medium_index(),
                                        path.beam().wavelength_nm(),
                                    )?
                                {
                                    if !pdf.is_finite() || pdf <= 0.0 {
                                        return Err(EnvironmentCameraError::NumericalRange);
                                    }
                                    partition = EndpointPartition::Ordinary(pdf);
                                }
                                Ok(())
                            })()
                        }
                        CurvedTransportEvent::MicrosurfaceDeparture(path) => {
                            emission.depart(path);
                            (|| {
                                partition = EndpointPartition::Original;
                                if let Some(receiving) = receiver {
                                    if supports_segment(
                                        camera.scene,
                                        path.medium_index(),
                                        path.beam().wavelength_nm(),
                                    )? {
                                        if include_dielectric
                                            && receiving.uses_smith_dielectric(camera)?
                                        {
                                            if receiving
                                                .supports_smith_dielectric_connections(camera)?
                                            {
                                                partition =
                                                    EndpointPartition::Dielectric(receiving);
                                            }
                                        } else if receiving.uses_smith(camera)? {
                                            partition = EndpointPartition::Smith(receiving);
                                        }
                                    }
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
                0,
                |_, _, _| Ok(None),
                |_, _, _, _| Err(CurvedSceneError::NumericalRange),
                |_, _, _| Ok(()),
                Some(SmithSurfaceSampler {
                    settings: self.microsurface,
                    samples: &mut facet_samples,
                    dielectric_samples: dielectric_samples
                        .as_mut()
                        .map(|samples| &mut **samples as &mut dyn FnMut(u32, u32) -> [f64; 4]),
                }),
            )
            .map(|(path, _)| path);
        if let Some(error) = source_error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        let mut radiance = camera.finish(cx, path, initial_index, 0.0)?;
        let endpoint_weight = if radiance.endpoint() == RadianceEndpoint::Environment {
            partition.weight(self.source, &path)?
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
}

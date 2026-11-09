//! Original molecular real-lens CIE/QE point histories and finite shutters.

use super::*;

impl CurvedCameraRaster<'_> {
    /// Observe each exact source at its own native group delay.
    /// The detector row precedes the physical lens; clipped lanes keep their mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_point_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighPointLightPixelSample<CurvedSequentialPixelSample>,
        TransientPointLightRasterError,
    > {
        let (observation_time_s, point) = self
            .sample_timed_sequential_rayleigh_point_lights_transport(
                cx,
                scene,
                initial_medium,
                Observation::Instant(observation_time_s),
                settings,
                seed,
                pixel,
                sample,
                |trace| {
                    self.sample_sequential_pixel_with_timing(
                        cx,
                        seed,
                        pixel,
                        sample,
                        scene.source_camera().analyzer_response(),
                        true,
                        |cx, beam, lane| trace(cx, *beam, lane),
                    )
                },
            )?;
        Ok(TransientRayleighPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// The detector row precedes the physical lens; clipped lanes keep their mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_point_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<ExposurePixelSample<CurvedSequentialPixelSample>>,
        TransientPointLightRasterError,
    > {
        let (_, point) = self.sample_timed_sequential_rayleigh_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_pixel_with_timing(
                    cx,
                    seed,
                    pixel,
                    sample,
                    scene.source_camera().analyzer_response(),
                    true,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )?;
        Ok(
            point
                .try_map_pixel(|mean| Ok::<_, TransientRasterError>(finish_pixel(window, mean)))?,
        )
    }

    /// Observe each exact source at its own native group delay.
    /// The detector row precedes the physical lens; clipped lanes keep their mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_point_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighPointLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        TransientPointLightRasterError,
    > {
        let (observation_time_s, point) = self
            .sample_timed_sequential_rayleigh_point_lights_transport(
                cx,
                scene,
                initial_medium,
                Observation::Instant(observation_time_s),
                settings,
                seed,
                pixel,
                sample,
                |trace| {
                    self.sample_sequential_sensor_with_timing(
                        cx,
                        sensor,
                        seed,
                        pixel,
                        sample,
                        scene.source_camera().analyzer_response(),
                        true,
                        |cx, beam, lane| trace(cx, *beam, lane),
                    )
                },
            )?;
        Ok(TransientRayleighPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// The detector row precedes the physical lens; clipped lanes keep their mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_point_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<
            ExposureSensorPixelSample<C, CurvedSequentialSensorPixelSample<C>>,
        >,
        TransientPointLightRasterError,
    > {
        let (_, point) = self.sample_timed_sequential_rayleigh_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_sensor_with_timing(
                    cx,
                    sensor,
                    seed,
                    pixel,
                    sample,
                    scene.source_camera().analyzer_response(),
                    true,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )?;
        Ok(point.try_map_pixel(|mean| {
            let energy = *mean.qe_weighted_radiance_w_m2_sr();
            let electrons = *mean.photoelectron_rate_per_m2_sr();
            finish_sensor_rates(window, mean, energy, electrons)
        })?)
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_timed_sequential_rayleigh_point_lights_transport<P>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                SequentialCameraBeam,
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<(Dd, RayleighPointLightPixelSample<P>), TransientPointLightRasterError> {
        // Admission precedes optional lens lanes, including completely clipped packets.
        let (clock, observation) = admit_observation(cx, observation)?;
        let point = self.sample_rayleigh_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, draws| match observation {
                Observation::Instant(time) => scene
                    .sample_prepared_sequential_ray_with_tracking_samples_at_time(
                        cx,
                        &beam,
                        initial_medium,
                        time,
                        settings,
                        |event| draws.surface(event),
                        |flight| draws.flight(flight),
                        |connection, flight| draws.source(connection, flight),
                    )
                    .map(|result| *result.point_light_sample()),
                Observation::Exposure(window) => scene
                    .sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
                        cx,
                        &beam,
                        initial_medium,
                        window,
                        settings,
                        |event| draws.surface(event),
                        |flight| draws.flight(flight),
                        |connection, flight| draws.source(connection, flight),
                    )
                    .map(|result| *result.mean_point_light_sample()),
            },
            integrate,
        )?;
        Ok((clock, point))
    }
}

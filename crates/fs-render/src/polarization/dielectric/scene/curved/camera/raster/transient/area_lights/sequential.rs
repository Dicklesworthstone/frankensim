//! Resolved-lens visible/QE pixels with retarded area-source NEE/MIS.

use super::super::super::CurvedSequentialPixelSample;
use super::super::super::sensor::CurvedSequentialSensorPixelSample;
use super::super::exposure::finish_sensor_rates;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::transient::admit_time;

impl CurvedCameraRaster<'_> {
    /// Trace the configured analyzer at the physical sensor, then the real
    /// lens, then the original finite-area camera estimator. Every wavelength
    /// retains its own full Stokes response and lens-plus-scene source clock.
    /// Clipped lanes retain zero contribution and their original proposal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_area_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientAreaLightPixelSample<CurvedSequentialPixelSample>,
        TransientAreaLightRasterError,
    > {
        let (observation_time_s, area) = self.sample_timed_sequential_area_transport(
            cx,
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            workspace,
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
        Ok(TransientAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// All authored QE channels share the same four optional real-lens paths
    /// and source images, including supported infrared bands. Lens geometry and
    /// QE/photon weights are applied once by the original sensor integrator.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_area_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        TransientAreaLightRasterError,
    > {
        let (observation_time_s, area) = self.sample_timed_sequential_area_transport(
            cx,
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            workspace,
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
        Ok(TransientAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// Frozen real-lens/source transport with each original history integrated
    /// over its own retarded sensor shutter. No time samples or retraces are
    /// introduced; per-lane NEE/MIS work remains in the returned area sample.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_area_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        AreaLightPixelSample<ExposurePixelSample<CurvedSequentialPixelSample>>,
        TransientAreaLightRasterError,
    > {
        let (_, area) = self.sample_timed_sequential_area_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            workspace,
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
        Ok(area.try_map_pixel(|mean| Ok::<_, TransientRasterError>(finish_pixel(window, mean)))?)
    }

    /// Original QE mean rates and sensor-shutter energy/electron integrals,
    /// before pixel etendue, detector gain, saturation or noise. Channels share
    /// the same transport and retain physically clipped wavelength samples.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_area_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        AreaLightPixelSample<ExposureSensorPixelSample<C, CurvedSequentialSensorPixelSample<C>>>,
        TransientAreaLightRasterError,
    > {
        let (_, area) = self.sample_timed_sequential_area_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            workspace,
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
        Ok(area.try_map_pixel(|mean| {
            let qe_rate = *mean.qe_weighted_radiance_w_m2_sr();
            let electron_rate = *mean.photoelectron_rate_per_m2_sr();
            finish_sensor_rates(window, mean, qe_rate, electron_rate)
        })?)
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_timed_sequential_area_transport<P>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
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
    ) -> Result<(Dd, AreaLightPixelSample<P>), TransientAreaLightRasterError> {
        // Admission is outside the surviving-lane callback: even a completely
        // clipped packet must honor cancellation and a valid sensor clock.
        cx.checkpoint()
            .map_err(|error| TransientAreaLightRasterError::Source(error.into()))?;
        let (clock, observation) = match observation {
            Observation::Instant(time) => {
                let time = admit_time(time).map_err(TransientAreaLightRasterError::Source)?;
                (time, Observation::Instant(time))
            }
            Observation::Exposure(window) => (window.open_s(), Observation::Exposure(window)),
        };
        let area = self.sample_area_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, samples| match observation {
                Observation::Instant(time) => scene
                    .sample_prepared_sequential_ray_with_samples_at_time(
                        cx,
                        &beam,
                        initial_medium,
                        time,
                        settings,
                        workspace,
                        |event| samples.surface(event),
                        |event, binding| samples.area(event, binding),
                    )
                    .map(|result| *result.area_sample()),
                Observation::Exposure(window) => scene
                    .sample_prepared_sequential_ray_with_samples_during_exposure(
                        cx,
                        &beam,
                        initial_medium,
                        window,
                        settings,
                        workspace,
                        |event| samples.surface(event),
                        |event, binding| samples.area(event, binding),
                    )
                    .map(|result| *result.mean_area_sample()),
            },
            integrate,
        )?;
        Ok((clock, area))
    }
}

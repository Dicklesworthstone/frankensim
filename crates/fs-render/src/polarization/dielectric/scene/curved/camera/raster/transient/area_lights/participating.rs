//! CIE/QE observations of timed real-lens molecular area-source transport.
//!
//! Native camera collisions, source survival and source-image work are shared
//! with steady molecular rendering. Only the physical sensor launch and each
//! finite surface source's own retarded history enter this detector adapter.

use super::super::super::CurvedSequentialPixelSample;
use super::super::super::area_lights::participating::RayleighAreaLightPixelSample;
use super::super::super::sensor::CurvedSequentialSensorPixelSample;
use super::super::exposure::finish_sensor_rates;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
    RayleighAreaLightCameraSettings, RayleighAreaLightCameraWorkspace,
};
use crate::polarization::dielectric::scene::curved::camera::transient::admit_time;
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::participating::TransientCurvedRayleighAreaLightCameraScene;

/// A complete sensor-clock molecular observation. All original camera and
/// independent source-work counters remain attached to their wavelength lane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighAreaLightPixelSample<P = CurvedSequentialPixelSample> {
    observation_time_s: Dd,
    area: RayleighAreaLightPixelSample<P>,
}

impl<P> TransientRayleighAreaLightPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        self.area.pixel_sample()
    }
    pub const fn area_sample(&self) -> &RayleighAreaLightPixelSample<P> {
        &self.area
    }
}

impl CurvedCameraRaster<'_> {
    /// Trace the configured analyzer at the physical sensor, then the real
    /// lens, then native molecular scattering and finite-area NEE/MIS. Every wavelength
    /// retains its full Stokes response and lens-plus-scene-plus-source clock.
    /// Clipped lanes retain zero contribution and their original proposal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_area_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighAreaLightPixelSample<CurvedSequentialPixelSample>,
        TransientAreaLightRasterError,
    > {
        let (observation_time_s, area) = self.sample_timed_sequential_rayleigh_area_transport(
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
        Ok(TransientRayleighAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// All authored QE channels share the same four optional real-lens paths
    /// and molecular/source histories, including supported infrared bands. Lens geometry and
    /// QE/photon weights are applied once by the original sensor integrator.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_area_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        TransientAreaLightRasterError,
    > {
        let (observation_time_s, area) = self.sample_timed_sequential_rayleigh_area_transport(
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
        Ok(TransientRayleighAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// Frozen real-lens/source transport with each original history integrated
    /// over its own retarded sensor shutter. No time samples or retraces are
    /// introduced; native molecular and source-survival work remains attached.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_area_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<ExposurePixelSample<CurvedSequentialPixelSample>>,
        TransientAreaLightRasterError,
    > {
        let (_, area) = self.sample_timed_sequential_rayleigh_area_transport(
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
    pub fn sample_pixel_sensor_sequential_rayleigh_area_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<
            ExposureSensorPixelSample<C, CurvedSequentialSensorPixelSample<C>>,
        >,
        TransientAreaLightRasterError,
    > {
        let (_, area) = self.sample_timed_sequential_rayleigh_area_transport(
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
    fn sample_timed_sequential_rayleigh_area_transport<P>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
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
    ) -> Result<(Dd, RayleighAreaLightPixelSample<P>), TransientAreaLightRasterError> {
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
        let area = self.sample_rayleigh_area_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, samples| match observation {
                Observation::Instant(time) => scene
                    .sample_prepared_sequential_ray_with_tracking_samples_at_time(
                        cx,
                        &beam,
                        initial_medium,
                        time,
                        settings,
                        workspace,
                        |event| samples.surface(event),
                        |flight| samples.flight(flight),
                        |event, binding| samples.area(event, binding),
                        |collision, binding| samples.volume_area(collision, binding),
                        |area, image, flight| samples.source(area, image, flight),
                    )
                    .map(|result| *result.area_sample()),
                Observation::Exposure(window) => scene
                    .sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
                        cx,
                        &beam,
                        initial_medium,
                        window,
                        settings,
                        workspace,
                        |event| samples.surface(event),
                        |flight| samples.flight(flight),
                        |event, binding| samples.area(event, binding),
                        |collision, binding| samples.volume_area(collision, binding),
                        |area, image, flight| samples.source(area, image, flight),
                    )
                    .map(|result| *result.mean_area_sample()),
            },
            integrate,
        )?;
        Ok((clock, area))
    }
}

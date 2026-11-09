//! Source histories on complete polarized Smith finite-area camera paths.
//!
//! The original estimator owns all facet walks, source proposals, endpoint
//! partitions, Stokes/eta transport and work. This observer changes only the
//! prescribed source radiance at that path's original retarded group clock.
//! A finite shutter uses each history's interval mean on one frozen transport.
//! The Smith boundary response is local: its dimensionless statistical height
//! is not a physical flight length or an inferred interface residence time.

use super::super::super::area_lights::camera::SmithDielectricAreaLightCameraScene;
use super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;

/// Allocation-free view of original histories and the complete opaque/glass
/// Smith source estimator. Preparation retains the native receiver/source
/// admission; no transport, wavelength or history is approximated anew.
pub struct TransientSmithDielectricAreaLightCameraScene<'v, 'a, 'h, 's, 'm, 'l> {
    source: &'v TransientCurvedAreaLightCameraScene<'a, 'h, 's, 'm, 'l>,
    lights: SmithDielectricAreaLightCameraScene<'a, 'h, 's, 'm, 'l>,
}

impl<'a, 'h, 's, 'm, 'l> TransientCurvedAreaLightCameraScene<'a, 'h, 's, 'm, 'l> {
    /// Prepare full two-sided Smith continuation and finite-area connection
    /// responses with these exact source histories. The original direct-source
    /// partition remains authoritative, including equal-index mixed measures.
    pub fn with_smith_dielectric<'v>(
        &'v self,
        cx: &Cx<'_>,
        microsurface: SmithWalkSettings,
    ) -> Result<
        TransientSmithDielectricAreaLightCameraScene<'v, 'a, 'h, 's, 'm, 'l>,
        TransientAreaLightCameraError,
    > {
        Ok(TransientSmithDielectricAreaLightCameraScene {
            source: self,
            lights: SmithDielectricAreaLightCameraScene::try_new(cx, self.lights, microsurface)?,
        })
    }
}

impl<'v, 'a, 'h, 's, 'm, 'l> TransientSmithDielectricAreaLightCameraScene<'v, 'a, 'h, 's, 'm, 'l> {
    pub const fn source_camera(&self) -> &'a TransientCurvedCameraScene<'h, 's, 'm> {
        self.source.camera
    }

    pub const fn source_area_camera(
        &self,
    ) -> &SmithDielectricAreaLightCameraScene<'a, 'h, 's, 'm, 'l> {
        &self.lights
    }

    /// Observe the original finite-area and continuing-emitter estimators at
    /// each source's own retarded time. All six native RNG domains are reused.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)
            .map_err(TransientCameraError::from)?;
        let area = self.observe(
            cx,
            beam,
            initial_medium,
            Observation::Instant(observation_time_s),
            Dd::ZERO,
            settings,
            workspace,
            seed,
            path_id,
        )?;
        Ok(TransientAreaLightCameraSample {
            observation_time_s,
            area,
        })
    }

    /// Return mean radiance and its physical shutter integral on one frozen
    /// Smith walk. Source histories supply their interval mean; no midpoint
    /// substitution, extra paths or temporal random draws are introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)
            .map_err(TransientCameraError::from)?;
        let area = self.observe(
            cx,
            beam,
            initial_medium,
            Observation::Exposure(window),
            Dd::ZERO,
            settings,
            workspace,
            seed,
            path_id,
        )?;
        exposure_sample(window, area)
    }

    /// Continue a qualified timed lens beam. Its analyzer must already precede
    /// the lens, and the original lens DD group delay is included exactly once
    /// in every source clock. The caller retains the separate pupil weight.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source.camera.analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_at_time(
            cx,
            optical,
            initial_medium,
            observation_time_s,
            settings,
            workspace,
            seed,
            path_id,
        )
    }

    /// Integrate each source over its own lens-plus-scene retarded interval.
    /// Geometry and optical coefficients remain frozen throughout the shutter.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source.camera.analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_during_exposure(
            cx,
            optical,
            initial_medium,
            window,
            settings,
            workspace,
            seed,
            path_id,
        )
    }

    // Only the native raster uses these seams after placing the configured
    // analyzer at the sensor, before generating the qualified timed lens beam.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let beam = self
            .source
            .camera
            .source_camera()
            .launch_sequential(optical, initial_medium)
            .map_err(TransientCameraError::from)?;
        let area = self.observe(
            cx,
            beam,
            initial_medium,
            Observation::Instant(observation_time_s),
            delay,
            settings,
            workspace,
            seed,
            path_id,
        )?;
        Ok(TransientAreaLightCameraSample {
            observation_time_s,
            area,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let delay = lens_delay(optical)?;
        let beam = self
            .source
            .camera
            .source_camera()
            .launch_sequential(optical, initial_medium)
            .map_err(TransientCameraError::from)?;
        let area = self.observe(
            cx,
            beam,
            initial_medium,
            Observation::Exposure(window),
            delay,
            settings,
            workspace,
            seed,
            path_id,
        )?;
        exposure_sample(window, area)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        observation: Observation,
        launch_delay_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, TransientAreaLightCameraError> {
        let mut source_error = None;
        let result = self.lights.sample_launched_with_modulation(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            seed,
            path_id,
            |path| {
                self.source
                    .scale(cx, observation, launch_delay_s, path)
                    .map_err(|error| {
                        source_error = Some(error);
                        AreaLightCameraError::NumericalRange
                    })
            },
        );
        if let Some(error) = source_error {
            return Err(error.into());
        }
        Ok(result?)
    }
}

fn exposure_sample(
    window: ExposureWindow,
    area: AreaLightCameraSample,
) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
    let integrated_radiance_j_m2_sr_nm = integrate_rate(
        area.radiance_sample().radiance_w_m2_sr_nm(),
        window.duration_s(),
    )
    .map_err(TransientCameraError::from)?;
    Ok(ExposureAreaLightCameraSample {
        window,
        area,
        integrated_radiance_j_m2_sr_nm,
    })
}

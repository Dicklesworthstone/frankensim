//! Finite shutter integration on frozen original refractive paths.
//!
//! Transport is linear in source radiance while geometry, media and surface
//! laws are fixed. One original path therefore suffices: each emitting hit
//! requests its source's mean over the shutter interval shifted by that hit's
//! own group delay. No temporal quadrature or extra transport variates occur.
//! The interval is rectangular at the original observer clock. Explicit timed
//! sequential beams include their lens transit before every source interval;
//! moving geometry and nonlinear sensor electronics remain separate models.

use super::*;

pub mod piecewise_linear;

/// Explicit finite rectangular shutter on the original physical clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureWindow {
    open_s: Dd,
    close_s: Dd,
    duration_s: Dd,
}

impl ExposureWindow {
    pub fn try_new(open_s: Dd, close_s: Dd) -> Result<Self, TransientCameraError> {
        let open_s = admit_time(open_s)?;
        let close_s = admit_time(close_s)?;
        let duration_s = close_s - open_s;
        if !open_s.lt(close_s) || !finite_dd(duration_s) || duration_s.hi <= 0.0 {
            return Err(TransientCameraError::InvalidExposureWindow);
        }
        Ok(Self {
            open_s,
            close_s,
            duration_s,
        })
    }

    pub const fn open_s(self) -> Dd {
        self.open_s
    }
    pub const fn close_s(self) -> Dd {
        self.close_s
    }
    pub const fn duration_s(self) -> Dd {
        self.duration_s
    }
}

/// The original emitting hit and both ends of its retarded source interval.
/// A source returns the mean multiplier, so a constant unit history retains
/// the original steady radiance exactly. Endpoints and duration remain DD.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceExposureQuery {
    pub opening: SurfaceEmissionQuery,
    pub observation_close_s: Dd,
    pub emission_close_s: Dd,
    pub duration_s: Dd,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureCameraSample {
    window: ExposureWindow,
    mean_radiance: CameraRadianceSample,
    integrated_radiance_j_m2_sr_nm: f64,
}

impl ExposureCameraSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    /// Mean spectral radiance; the retained path is the one original trace.
    pub const fn mean_radiance_sample(&self) -> &CameraRadianceSample {
        &self.mean_radiance
    }
    /// Time integral of spectral radiance, before any detector calibration.
    pub const fn integrated_radiance_j_m2_sr_nm(&self) -> f64 {
        self.integrated_radiance_j_m2_sr_nm
    }
}

impl TransientCurvedCameraScene<'_, '_, '_> {
    /// Integrate source histories over one shutter using original seeded work.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            Dd::ZERO,
            settings,
            seed,
            path_id,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        window: ExposureWindow,
        launch_delay_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut error = None;
        let result = self.camera.scene.trace_with_seed_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            true,
            &self.camera.spatial_volume_emitters,
            |event| {
                self.observe_exposure_with_delay(
                    cx,
                    window,
                    launch_delay_s,
                    event,
                    &mut emission,
                    &mut error,
                )
            },
        );
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        finish(
            window,
            self.camera.finish(cx, path, initial_index, emitted)?,
        )
    }

    /// One caller-sampled path; histories never consume event random variates.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            Dd::ZERO,
            settings,
            samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        window: ExposureWindow,
        launch_delay_s: Dd,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut error = None;
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            true,
            &self.camera.spatial_volume_emitters,
            |event| {
                self.observe_exposure_with_delay(
                    cx,
                    window,
                    launch_delay_s,
                    event,
                    &mut emission,
                    &mut error,
                )
            },
        );
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        finish(
            window,
            self.camera.finish(cx, path, initial_index, emitted)?,
        )
    }

    pub(super) fn observe_exposure(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        event: CurvedTraceEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        self.observe_exposure_with_delay(cx, window, Dd::ZERO, event, emission, error)
    }

    pub(super) fn observe_exposure_with_delay(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        launch_delay_s: Dd,
        event: CurvedTraceEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedTraceEvent::Arrival(path) => {
                let scale = self
                    .exposure_scale_with_delay(cx, window, launch_delay_s, path)
                    .map_err(|failure| {
                        *error = Some(failure);
                        CurvedSceneError::NumericalRange
                    })?;
                emission.observe_modulated(path, scale)
            }
            CurvedTraceEvent::Departure { path, .. } => {
                emission.depart(path);
                Ok(())
            }
        }
    }

    pub(super) fn exposure_scale(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        path: &CurvedPathResult,
    ) -> Result<f64, TransientCameraError> {
        self.exposure_scale_with_delay(cx, window, Dd::ZERO, path)
    }

    pub(super) fn exposure_scale_with_delay(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        launch_delay_s: Dd,
        path: &CurvedPathResult,
    ) -> Result<f64, TransientCameraError> {
        let Some((binding, opening)) =
            self.source_query_with_delay(cx, window.open_s, launch_delay_s, path)?
        else {
            return Ok(1.0);
        };
        average_history(cx, binding, opening, window)
    }
}

// Shared source-clock interval check; callers preserve their own source identity.
pub(super) fn checked_exposure_close(
    emission_open_s: Dd,
    group_delay_s: Dd,
    window: ExposureWindow,
) -> Option<Dd> {
    let emission_close_s = checked_retarded_time(window.close_s, group_delay_s)?;
    let duration_error = ((emission_close_s - emission_open_s) - window.duration_s).abs();
    if !emission_open_s.lt(emission_close_s)
        || (window.duration_s * Dd::from_f64(64.0 * f64::EPSILON)).lt(duration_error)
    {
        return None;
    }
    Some(emission_close_s)
}

pub(super) fn average_history(
    cx: &Cx<'_>,
    binding: SurfaceEmissionHistoryBinding<'_>,
    opening: SurfaceEmissionQuery,
    window: ExposureWindow,
) -> Result<f64, TransientCameraError> {
    let object_id = binding.object_id;
    let emission_close_s = checked_exposure_close(opening.emission_time_s, opening.group_delay_s, window)
        .ok_or(TransientCameraError::UnresolvedClock { object_id })?;
    let scale = binding
        .history
        .average(
            cx,
            SurfaceExposureQuery {
                opening,
                observation_close_s: window.close_s,
                emission_close_s,
                duration_s: window.duration_s,
            },
        )
        .map_err(|error| TransientCameraError::History { object_id, error })?;
    cx.checkpoint()?;
    if !scale.is_finite() || scale < 0.0 {
        return Err(TransientCameraError::InvalidHistoryScale { object_id });
    }
    Ok(scale)
}

pub(super) fn finish(
    window: ExposureWindow,
    mean_radiance: CameraRadianceSample,
) -> Result<ExposureCameraSample, TransientCameraError> {
    let integrated_radiance_j_m2_sr_nm =
        integrate_rate(mean_radiance.radiance_w_m2_sr_nm(), window.duration_s)?;
    Ok(ExposureCameraSample {
        window,
        mean_radiance,
        integrated_radiance_j_m2_sr_nm,
    })
}

/// Shared final unit conversion for radiance and already QE-weighted rates.
pub(crate) fn integrate_rate(rate: f64, duration_s: Dd) -> Result<f64, CameraRadianceError> {
    let integral = (Dd::from_f64(rate) * duration_s).to_f64();
    if !integral.is_finite() || integral < 0.0 || (rate > 0.0 && integral == 0.0) {
        Err(CameraRadianceError::NumericalRange)
    } else {
        Ok(integral)
    }
}

//! Timed real-lens seeds on the original refractive molecular history walk.

use super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;

mod exposure;

impl TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_> {
    /// Continue the complete lens detector state with the original independent
    /// surface and flight counters. Every source retains the supplied sensor
    /// clock and receives its own lens-plus-scene group delay exactly once.
    ///
    /// The beam must carry qualified lens timing. It also owns any analyzer
    /// installed at the sensor; an independently configured transient-scene
    /// analyzer refuses here, because applying it after the lens is incorrect.
    /// The sequential transient raster installs that analyzer before the lens.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_at_time(
            cx,
            beam,
            initial_medium,
            time,
            delay,
            settings,
            seed,
            path_id,
        )
    }

    /// Uniform-density flights keep their original [CDF,cosine,azimuth]
    /// variates. Timing changes neither collision proposals nor Mueller power.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    /// Uniform or heterogeneous flights retain the independent fourth
    /// acceptance variate. Null collisions preserve the original polarization,
    /// accumulated scene delay and work counts, and query no source history.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    // The raster installed this scene's analyzer at the sensor while generating
    // the timed lens beam. This seam admits the complete original seed without
    // applying that analyzer again at the lens exit.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let beam = self
            .camera
            .source_camera()
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    fn launch_timed_sequential(
        &self,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
    ) -> Result<(PolarizedBeam, Dd), TransientCameraError> {
        if self.camera.analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens);
        }
        let delay = lens_delay(optical)?;
        let beam = self
            .camera
            .source_camera()
            .launch_sequential(optical, initial_medium)?;
        Ok((beam, delay))
    }
}

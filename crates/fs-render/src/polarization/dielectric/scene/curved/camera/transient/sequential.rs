//! Sensor-clock histories continued from an opt-in timed real-lens beam.

use super::*;
use crate::camera::SequentialCameraBeam;

mod exposure;

impl TransientCurvedCameraScene<'_, '_, '_> {
    /// Continue a timed sequential beam through the original native history
    /// walk. Every source sees the supplied sensor observation clock and its
    /// own lens-plus-scene group delay. Lens Fresnel/Beer and endpoint radiance
    /// normalization are already in the reciprocal beam; geometry is separate.
    ///
    /// The supplied beam owns its detector response, including any analyzer
    /// applied at the physical sensor. An independently configured
    /// [`Self::with_analyzer`] is refused here: applying it at the lens exit
    /// would change the physics. The high-level sequential transient raster
    /// applies that configured analyzer before tracing the lens instead.
    /// Ordinary steady beams have no qualified group delay and are refused.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        if self.analyzer.is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens);
        }
        let delay = lens_delay(optical)?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
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

    /// Caller-supplied native facet/branch variates retain their original
    /// addresses. Source histories add no random draws or extra path traces.
    /// Analyzer and group-delay admission match [`Self::sample_sequential_ray_at_time`].
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        if self.analyzer.is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens);
        }
        self.sample_prepared_sequential_ray_with_samples_at_time(
            cx,
            optical,
            initial_medium,
            time,
            settings,
            samples,
        )
    }

    // The raster is responsible for installing this scene's analyzer at the
    // sensor during timed beam generation. This restricted seam never applies
    // it a second time after the lens, and is not a public raw-beam bypass.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn sample_prepared_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            delay,
            settings,
            samples,
        )
    }
}

pub(super) fn lens_delay(optical: &SequentialCameraBeam) -> Result<Dd, TransientCameraError> {
    optical
        .group_delay_s()
        .filter(|delay| finite_dd(*delay) && Dd::ZERO.lt(*delay))
        .ok_or(TransientCameraError::MissingLensGroupDelay)
}

#[cfg(test)]
mod tests;

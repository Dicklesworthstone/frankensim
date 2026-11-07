//! Exact finite-source shutter means at the physical sequential-lens sensor.

use super::super::exposure::{ExposureCameraSample, ExposureWindow};
use super::*;

impl TransientCurvedCameraScene<'_, '_, '_> {
    /// Trace one timed beam and integrate each source history over its own
    /// sensor shutter minus lens-plus-scene group delay. Transport, absorption
    /// and polarization stay frozen during this rectangular exposure.
    ///
    /// The beam already owns any analyzer applied at the physical sensor;
    /// an independently configured scene analyzer is refused. Histories must
    /// implement their original interval average: no point-time quadrature or
    /// extra random draws are substituted. Rear-cap geometry remains separate.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        if self.analyzer.is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens);
        }
        let delay = lens_delay(optical)?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            seed,
            path_id,
        )
    }

    /// The same finite shutter using caller-supplied native surface variates.
    /// Each emitting hit receives one source interval, without a retrace.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        if self.analyzer.is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens);
        }
        self.sample_prepared_sequential_ray_with_samples_during_exposure(
            cx,
            optical,
            initial_medium,
            window,
            settings,
            samples,
        )
    }

    // Raster generation installs the scene's analyzer before the lens. This
    // restricted entry point must not apply that analyzer a second time.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let delay = lens_delay(optical)?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            samples,
        )
    }
}

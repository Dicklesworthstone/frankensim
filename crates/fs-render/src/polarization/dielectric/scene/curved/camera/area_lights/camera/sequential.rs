//! Full reciprocal real-lens beams continued through native finite-area lighting.

use super::*;
use crate::camera::SequentialCameraBeam;

impl CurvedAreaLightCameraScene<'_, '_, '_, '_> {
    /// Continue the complete steady lens response through original finite-area
    /// next-event estimation and complementary emitter-hit MIS. Polarization,
    /// Fresnel/coating transmission and Beer attenuation remain in the beam.
    ///
    /// Group timing is not required: phase-qualified coated lenses and steady
    /// material tables retain their original domain. The beam already owns any
    /// sensor analyzer. Apply its separate geometric sampling weight once when
    /// integrating the resulting radiance over the physical pupil.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.sample_with_modulation(
            cx,
            || self.camera.launch_sequential(optical, initial_medium),
            initial_medium,
            settings,
            workspace,
            surface_samples,
            area_samples,
            |_| Ok(1.0),
        )
    }
}

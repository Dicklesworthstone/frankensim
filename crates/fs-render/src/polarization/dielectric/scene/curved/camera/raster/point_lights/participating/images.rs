//! Native spectral pixels from distinct refractive images of point fixtures.

use super::{CurvedCameraRaster, PointLightRasterError, RayleighPointLightPixelSample, unit};
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::RayleighPointLightCameraSettings;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::CurvedRayleighImageLightCameraScene;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

impl CurvedCameraRaster<'_> {
    /// Sum discovered regular optical images at each original camera vertex.
    /// The caller prepares one workspace per worker and reuses it across
    /// wavelength lanes, vertices and pixels. No scene cache or allocation is
    /// introduced by this raster operation. Any unresolved seed or lane refuses
    /// the whole pixel; finite authored seeds do not establish image completeness.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_image_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample, PointLightRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_pixel_rayleigh_point_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                scene.sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    settings,
                    workspace,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                    |connection, flight| {
                        let source_key = [key[0], key[1] ^ (0x6372_7360 + lane as u32)];
                        unit(philox4x32_10([pixel, sample, connection, flight], source_key)[0])
                    },
                )
            },
        )
    }
}

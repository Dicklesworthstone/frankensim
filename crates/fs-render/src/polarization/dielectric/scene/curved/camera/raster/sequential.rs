//! Resolved lens apertures and reciprocal polarization on native curved paths.

use super::*;

#[cfg(feature = "volumes")]
mod participating;
#[cfg(feature = "volumes")]
pub use participating::CurvedSequentialRayleighPixelSample;

/// One visible sample through a resolved sequential lens. A rejected optical
/// lane has no outgoing ray or scene path and remains a zero-contribution draw.
/// Average complete samples normally; never renormalize surviving wavelengths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedSequentialPixelSample {
    xyz: [f64; 3],
    packet: OptionalCameraPacket,
}

impl CurvedSequentialPixelSample {
    /// Equal-energy XYZ with the original four-wavelength integration measure.
    /// The camera's aperture-normalized radiance convention is retained.
    pub const fn xyz(&self) -> [f64; 3] {
        self.xyz
    }
    pub const fn camera_ray(&self) -> Option<Ray> {
        self.packet.rays[0]
    }
    pub const fn camera_rays(&self) -> &[Option<Ray>; CAMERA_WAVELENGTHS] {
        &self.packet.rays
    }
    pub const fn wavelengths_nm(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.packet.wavelengths_nm
    }
    /// Native paths already include the lens's reciprocal Fresnel/Beer seed.
    /// Clipped wavelengths remain `None`, not a fabricated black transport.
    pub const fn spectral_samples(&self) -> &[Option<CameraRadianceSample>; CAMERA_WAVELENGTHS] {
        &self.packet.spectral
    }
    /// Projected-cap geometry factors applied once to [`Self::xyz`]. Fresnel
    /// and lens absorption are already in the retained spectral radiances.
    pub const fn geometry_weights(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.packet.optical_weights
    }
}

impl CurvedCameraRaster<'_> {
    /// Trace a physical sequential lens into the original polarized curved
    /// scene. Each surviving wavelength starts at its actual front-surface
    /// point with the full reciprocal lens Stokes state. The scene exterior
    /// must have local phase index one at every launch; the original scene
    /// owner admits the selected medium, material spectra and path budgets.
    ///
    /// Aperture/TIR rejection keeps its original zero-contribution sample.
    /// Numerical, source-support and cancellation errors still abort the whole
    /// pixel. Use [`Self::sample_pixel`] for the existing ideal/wave result.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        self.sample_sequential_pixel(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            None,
        )
    }

    /// Place a passive analyzer at the physical sensor, before reciprocal
    /// lens propagation. Its wavelength-specific Stokes seed retains the
    /// correlation of the analyzer with both lens and scene polarization.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        self.sample_sequential_pixel(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(analyzer),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_pixel(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<CameraAnalyzer>,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sequential_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            analyzer,
            |cx, beam, lane| {
                scene.sample_sequential_ray_with_samples(
                    cx,
                    beam,
                    initial_medium,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_pixel_with_transport(
        &self,
        cx: &Cx<'_>,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<CameraAnalyzer>,
        mut trace: impl FnMut(
            &Cx<'_>,
            &SequentialCameraBeam,
            usize,
        ) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        if self.camera.sequential_lens().is_none() {
            return Err(CameraError::InvalidSequentialLens.into());
        }
        let packet = self.sample_optional_packet_with_transport(
            cx,
            seed,
            pixel,
            sample,
            analyzer,
            |u| {
                Ok(hero_wavelength_packet::<CAMERA_WAVELENGTHS>(
                    LAMBDA_MIN + u * (LAMBDA_MAX - LAMBDA_MIN),
                    LAMBDA_MIN,
                    LAMBDA_MAX,
                ))
            },
            |cx, launch, _, lane| {
                let CameraLaunch::Sequential(beam) = launch else {
                    return Err(CameraError::InvalidSequentialLens.into());
                };
                Ok(trace(cx, &beam, lane)?)
            },
        )?;
        let mut xyz = [0.0; 3];
        for ((lane, &wavelength), &geometry) in packet
            .spectral
            .iter()
            .zip(&packet.wavelengths_nm)
            .zip(&packet.optical_weights)
        {
            let Some(lane) = lane else {
                continue;
            };
            let weight = weighted_radiance(lane.radiance_w_m2_sr_nm(), geometry)? * self.xyz_scale;
            for (sum, observer) in
                xyz.iter_mut()
                    .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
            {
                *sum += weight * observer;
            }
        }
        if xyz.iter().any(|value| !value.is_finite()) {
            return Err(CurvedRasterError::NumericalRange);
        }
        cx.checkpoint().map_err(CameraError::from)?;
        Ok(CurvedSequentialPixelSample { xyz, packet })
    }
}

#[cfg(test)]
mod tests;

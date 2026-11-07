//! Source-bounded QE observations through the original resolved camera lens.

use super::*;

/// Detector integrals with actual wavelength-dependent optical rejection.
/// Optional paths carry the full lens/scene polarization and Beer response;
/// a clipped wavelength is retained in the proposal with zero contribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedSequentialSensorPixelSample<const CHANNELS: usize> {
    qe_weighted_radiance_w_m2_sr: [f64; CHANNELS],
    photoelectron_rate_per_m2_sr: [f64; CHANNELS],
    wavelengths: [DetectorWavelength<CHANNELS>; CAMERA_WAVELENGTHS],
    packet: OptionalCameraPacket,
}

impl<const C: usize> CurvedSequentialSensorPixelSample<C> {
    pub const fn qe_weighted_radiance_w_m2_sr(&self) -> &[f64; C] {
        &self.qe_weighted_radiance_w_m2_sr
    }
    /// QE*lambda_vac/(h*c) integral, before etendue/exposure/electronics.
    pub const fn photoelectron_rate_per_m2_sr(&self) -> &[f64; C] {
        &self.photoelectron_rate_per_m2_sr
    }
    pub const fn camera_ray(&self) -> Option<Ray> {
        self.packet.rays[0]
    }
    pub const fn camera_rays(&self) -> &[Option<Ray>; CAMERA_WAVELENGTHS] {
        &self.packet.rays
    }
    pub const fn detector_wavelengths(&self) -> &[DetectorWavelength<C>; CAMERA_WAVELENGTHS] {
        &self.wavelengths
    }
    /// Radiances already contain lens Fresnel/Beer and original scene transport.
    pub const fn spectral_samples(&self) -> &[Option<CameraRadianceSample>; CAMERA_WAVELENGTHS] {
        &self.packet.spectral
    }
    /// Cap-area geometry applied once to both channel integrals. This excludes
    /// the lens transmission already carried by each reciprocal native path.
    pub const fn geometry_weights(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.packet.optical_weights
    }
}

impl CurvedCameraRaster<'_> {
    /// Observe the resolved lens through the same sum-of-QE proposal as
    /// [`Self::sample_pixel_sensor`]. All channels reuse four independent
    /// native paths. Optical rejection keeps the original one-quarter packet
    /// measure; no retries, wavelength substitution or surviving-lane average.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialSensorPixelSample<C>, CurvedRasterError> {
        self.sample_sequential_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            None,
        )
    }

    /// The same QE observation with a passive analyzer at the original
    /// sensor. Lens and scene Stokes transport are composed before integrating
    /// any channel; no scalar analyzer factor is applied afterwards.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialSensorPixelSample<C>, CurvedRasterError> {
        self.sample_sequential_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(analyzer),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<CameraAnalyzer>,
    ) -> Result<CurvedSequentialSensorPixelSample<C>, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sequential_sensor_with_transport(
            cx,
            sensor,
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
    pub(in super::super) fn sample_sequential_sensor_with_transport<const C: usize>(
        &self,
        cx: &Cx<'_>,
        sensor: &QuantumEfficiencySensor<C>,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<CameraAnalyzer>,
        mut trace: impl FnMut(
            &Cx<'_>,
            &SequentialCameraBeam,
            usize,
        ) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<CurvedSequentialSensorPixelSample<C>, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        if self.camera.sequential_lens().is_none() {
            return Err(CameraError::InvalidSequentialLens.into());
        }
        let mut draws = None;
        let packet = self.sample_optional_packet_with_transport(
            cx,
            seed,
            pixel,
            sample,
            analyzer,
            |u| {
                let weights = detector_packet(sensor, u)?;
                draws = Some(weights);
                Ok(weights.map(|draw| draw.wavelength_nm()))
            },
            |cx, launch, _, lane| {
                let CameraLaunch::Sequential(beam) = launch else {
                    return Err(CameraError::InvalidSequentialLens.into());
                };
                Ok(trace(cx, &beam, lane)?)
            },
        )?;
        let wavelengths = draws.ok_or(CurvedRasterError::NumericalRange)?;
        let mut power = [Dd::ZERO; C];
        let mut electrons = [Dd::ZERO; C];
        for ((lane, draw), &geometry) in packet
            .spectral
            .iter()
            .zip(&wavelengths)
            .zip(&packet.optical_weights)
        {
            let Some(lane) = lane else {
                continue;
            };
            let radiance = weighted_radiance(lane.radiance_w_m2_sr_nm(), geometry)?;
            for channel in 0..C {
                power[channel] =
                    power[channel] + weighted_lane(radiance, draw.radiance_weights_nm()[channel]);
                electrons[channel] = electrons[channel]
                    + weighted_lane(radiance, draw.electron_weights_nm_per_joule()[channel]);
            }
        }
        let qe_weighted_radiance_w_m2_sr = power.map(Dd::to_f64);
        let photoelectron_rate_per_m2_sr = electrons.map(Dd::to_f64);
        if qe_weighted_radiance_w_m2_sr
            .iter()
            .chain(&photoelectron_rate_per_m2_sr)
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CurvedRasterError::NumericalRange);
        }
        cx.checkpoint().map_err(CameraError::from)?;
        Ok(CurvedSequentialSensorPixelSample {
            qe_weighted_radiance_w_m2_sr,
            photoelectron_rate_per_m2_sr,
            wavelengths,
            packet,
        })
    }
}

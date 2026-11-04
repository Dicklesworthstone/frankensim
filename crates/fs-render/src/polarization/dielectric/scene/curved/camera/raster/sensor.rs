//! Multichannel physical detector observations of original spectral ray paths.
//!
//! Four stratified quantiles of the sum-of-QE proposal independently trace the
//! native scene. Every channel reuses those same four paths with its own exact
//! importance weight. Preparation is shared, per-pixel storage is fixed, and
//! increasing the channel count does not increase geometry/transport work.
//! Absolute `(seed, pixel, sample)` identities retain scheduling/replay behavior.
//!
//! Outputs are integrals of incident spectral radiance against authored QE and
//! QE*lambda_vac/(h*c). They precede pixel etendue, exposure and electronics.
//! The ideal lens retains the original aperture-average convention; no sensor
//! irradiance, calibrated count, shot noise, gain or display color is inferred.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;
use crate::spectral::detector::{DetectorWavelength, QuantumEfficiencySensor};
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedSensorPixelSample<const CHANNELS: usize> {
    qe_weighted_radiance_w_m2_sr: [f64; CHANNELS],
    photoelectron_rate_per_m2_sr: [f64; CHANNELS],
    ray: Ray,
    wavelengths: [DetectorWavelength<CHANNELS>; CAMERA_WAVELENGTHS],
    spectral: [CameraRadianceSample; CAMERA_WAVELENGTHS],
}

impl<const C: usize> CurvedSensorPixelSample<C> {
    /// Monte Carlo contribution to each integral of L_lambda*QE(lambda) dλ.
    /// Average completed samples without a CIE/display or white normalization.
    pub const fn qe_weighted_radiance_w_m2_sr(&self) -> &[f64; C] {
        &self.qe_weighted_radiance_w_m2_sr
    }
    /// Contribution to each integral of L_lambda*QE(lambda)*lambda/(h*c) dλ.
    /// Collected-electron rate per m² sr, before etendue and exposure calibration.
    pub const fn photoelectron_rate_per_m2_sr(&self) -> &[f64; C] {
        &self.photoelectron_rate_per_m2_sr
    }
    pub const fn camera_ray(&self) -> Ray {
        self.ray
    }
    pub const fn detector_wavelengths(&self) -> &[DetectorWavelength<C>; CAMERA_WAVELENGTHS] {
        &self.wavelengths
    }
    pub const fn spectral_samples(&self) -> &[CameraRadianceSample; CAMERA_WAVELENGTHS] {
        &self.spectral
    }
}

/// The same detector observation through native molecular/particle scattering.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedRayleighSensorPixelSample<const CHANNELS: usize> {
    pixel: CurvedSensorPixelSample<CHANNELS>,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<const C: usize> CurvedRayleighSensorPixelSample<C> {
    pub const fn pixel_sample(&self) -> &CurvedSensorPixelSample<C> {
        &self.pixel
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    /// Includes rejected tracking candidates and surface-interrupted flights.
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

impl CurvedCameraRaster<'_> {
    /// Trace four original optical paths across an explicit detector band.
    /// Materials and emission retain their own source admission at every
    /// queried wavelength; unsupported spectra refuse, never become black.
    /// Each lane retains the native path budgets. Failure publishes no partial
    /// channel result and never retries or renormalizes surviving lanes.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                scene.sample_ray_with_samples(
                    cx,
                    ray,
                    wavelength,
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

    /// Observe the same source-bounded detector through the original scattering
    /// transport. Each wavelength retains its own source spectrum, phase law,
    /// geometry and free-flight counter. Zero scattering preserves the entire
    /// surface-only detector result; unresolved volume paths abort all channels.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedRayleighSensorPixelSample<C>, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let pixel = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let traced = scene.sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                )?;
                volume_events[lane] = traced.volume_events();
                free_flight_draws[lane] = traced.free_flight_draws();
                Ok(*traced.radiance_sample())
            },
        )?;
        Ok(CurvedRayleighSensorPixelSample {
            pixel,
            volume_events,
            free_flight_draws,
        })
    }

    fn sample_sensor_with_transport<const C: usize>(
        &self,
        cx: &Cx<'_>,
        sensor: &QuantumEfficiencySensor<C>,
        seed: u64,
        pixel: u32,
        sample: u32,
        trace: impl FnMut(&Cx<'_>, Ray, f64, usize) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        let mut draws = None;
        let (ray, _, spectral) = self.sample_packet_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |u| {
                let weights = [
                    sensor.sample_wavelength(u)?,
                    sensor.sample_wavelength((u + 0.25).fract())?,
                    sensor.sample_wavelength((u + 0.5).fract())?,
                    sensor.sample_wavelength((u + 0.75).fract())?,
                ];
                let wavelengths = weights.map(|draw| draw.wavelength_nm());
                draws = Some(weights);
                Ok(wavelengths)
            },
            trace,
        )?;
        let wavelengths = draws.ok_or(CurvedRasterError::NumericalRange)?;
        let mut power = [Dd::ZERO; C];
        let mut electrons = [Dd::ZERO; C];
        for (lane, draw) in spectral.iter().zip(&wavelengths) {
            let radiance = lane.radiance_w_m2_sr_nm();
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
        Ok(CurvedSensorPixelSample {
            qe_weighted_radiance_w_m2_sr,
            photoelectron_rate_per_m2_sr,
            ray,
            wavelengths,
            spectral,
        })
    }
}

// Divide the larger factor first, preserving tiny radiances multiplied by
// large photon weights and avoiding a needless large-product intermediate.
fn weighted_lane(radiance: f64, weight: f64) -> Dd {
    if radiance >= weight {
        Dd::from_f64(radiance / CAMERA_WAVELENGTHS as f64) * Dd::from_f64(weight)
    } else {
        Dd::from_f64(radiance) * Dd::from_f64(weight / CAMERA_WAVELENGTHS as f64)
    }
}

#[cfg(test)]
mod tests;

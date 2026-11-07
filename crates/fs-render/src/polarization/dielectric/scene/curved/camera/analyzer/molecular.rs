//! Passive reception through the original reciprocal molecular Mueller walk.
//!
//! Rayleigh and admitted depolarized-Rayleigh scattering obey the same
//! A=diag(1,1,-1,1) reciprocity as the surface owner. Seed the existing analyzer
//! once at the camera and retain its complete Stokes response through every
//! frame rotation, interface and volume event. Scattering survival belongs to
//! the original free-flight proposal; no extra scalar extinction or analyzer
//! throughput is applied after transport.

use super::super::CameraRadianceError;
use super::super::participating::{CurvedRayleighCameraScene, RayleighCameraSample};
use super::CameraAnalyzer;
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;
use fs_exec::Cx;

impl CurvedRayleighCameraScene<'_, '_, '_, '_> {
    /// Retain the original surface and molecular counter streams, thermal and
    /// surface sources, field identities and local radiance index conversion.
    /// The analyzer is the same passive launch-plane response as for clear
    /// camera rays. Authored surface polarization uses the same source-frame
    /// contraction at native hits; unresolved tails still refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.surface_camera(), ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Uniform-density flights use the original [CDF, cosine, azimuth]
    /// variates. The scalar angular proposal divides the full Stokes response
    /// once; it is not conditioned or resampled for the chosen analyzer.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.surface_camera(), ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    /// Uniform or heterogeneous density with original [majorant CDF, cosine,
    /// azimuth, acceptance] variates. Null candidates preserve the full Stokes
    /// state; all original tracking, volume, distance and work limits apply.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.surface_camera(), ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
        )
    }
}

#[cfg(test)]
mod tests;

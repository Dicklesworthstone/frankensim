//! Reciprocal camera-lens state continued through the original native scene.
//!
//! A sequential lens transmits a complete detector Stokes seed. Resetting it
//! to an unpolarized beam and multiplying two scalar transmissions would lose
//! polarization correlations between lens and scene interfaces. These entry
//! points preserve the outgoing frame, Q/U/V and absolute lens throughput.
//! Only the separate rear-cap sampling Jacobian remains for raster integration.

use super::{CameraRadianceError, CameraRadianceSample, CurvedCameraScene, phase_index};
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_exec::Cx;

impl CurvedCameraScene<'_, '_> {
    /// Continue a real-lens detector seed with the original surface counter
    /// stream, emission, local n² conversion and curved-path work budgets.
    /// The sample includes lens Fresnel/Beer losses, but not
    /// [`SequentialCameraBeam::geometry_weight`]. Its retained scene path starts
    /// at the lens's front surface, not at the sensor.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// The same optical continuation using native [facet_u,facet_v,branch]
    /// variates. Proposals remain independent of the detector state, retaining
    /// every native Mueller/proposal division exactly once.
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples(cx, beam, initial_medium, settings, samples)
    }

    fn launch_sequential(
        &self,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
    ) -> Result<PolarizedBeam, CameraRadianceError> {
        let beam = self.scene.launch_beam(
            optical.ray().origin,
            optical.frame(),
            optical.wavelength_nm(),
            initial_medium,
            optical.power(),
        )?;
        // Query the original local material, including its spectral support.
        // A coincident geometric boundary cannot stand in for this admission.
        if phase_index(beam.medium(), beam.wavelength_nm())? != 1.0 {
            return Err(CameraRadianceError::IncompatibleLensExterior);
        }
        Ok(beam)
    }
}

#[cfg(feature = "volumes")]
mod molecular {
    use super::*;
    use crate::polarization::dielectric::scene::curved::camera::participating::{
        CurvedRayleighCameraScene, RayleighCameraSample,
    };
    use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;

    impl CurvedRayleighCameraScene<'_, '_, '_, '_> {
        /// Retain lens polarization through the original molecular Mueller
        /// walk and its independent surface/free-flight counter streams.
        /// Lens transmission, scattering survival and analyzer response are
        /// each applied by their original owner, never as post-hoc scalars.
        #[allow(clippy::too_many_arguments)]
        pub fn sample_sequential_ray(
            &self,
            cx: &Cx<'_>,
            optical: &SequentialCameraBeam,
            initial_medium: usize,
            settings: CurvedRayleighTraceSettings,
            seed: u64,
            path_id: u64,
        ) -> Result<RayleighCameraSample, CameraRadianceError> {
            cx.checkpoint()?;
            let beam = self
                .surface_camera()
                .launch_sequential(optical, initial_medium)?;
            self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
        }

        /// Uniform-density flights retain [CDF,cosine,azimuth] variates and
        /// the original path/source refusal boundaries.
        #[allow(clippy::too_many_arguments)]
        pub fn sample_sequential_ray_with_samples(
            &self,
            cx: &Cx<'_>,
            optical: &SequentialCameraBeam,
            initial_medium: usize,
            settings: CurvedRayleighTraceSettings,
            surface_samples: impl FnMut(u32) -> [f64; 3],
            flight_samples: impl FnMut(u32) -> [f64; 3],
        ) -> Result<RayleighCameraSample, CameraRadianceError> {
            cx.checkpoint()?;
            let beam = self
                .surface_camera()
                .launch_sequential(optical, initial_medium)?;
            self.sample_launched_ray_with_samples(
                cx,
                beam,
                initial_medium,
                settings,
                surface_samples,
                flight_samples,
            )
        }

        /// Uniform or heterogeneous density with the original independent
        /// [majorant CDF,cosine,azimuth,acceptance] variates. Rejected tracking
        /// candidates preserve the entire lens-seeded Stokes state.
        #[allow(clippy::too_many_arguments)]
        pub fn sample_sequential_ray_with_tracking_samples(
            &self,
            cx: &Cx<'_>,
            optical: &SequentialCameraBeam,
            initial_medium: usize,
            settings: CurvedRayleighTraceSettings,
            surface_samples: impl FnMut(u32) -> [f64; 3],
            flight_samples: impl FnMut(u32) -> [f64; 4],
        ) -> Result<RayleighCameraSample, CameraRadianceError> {
            cx.checkpoint()?;
            let beam = self
                .surface_camera()
                .launch_sequential(optical, initial_medium)?;
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
}

#[cfg(test)]
mod tests;

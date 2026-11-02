//! Unpolarized camera radiance with molecular scattering on curved paths.
//!
//! Reduced radiance L/n² has the same local absorption and scattering rates as
//! the native reciprocal power walk. Its scattering-only free-flight density
//! already includes both scattering survival and the collision rate; native
//! propagation supplies absorption, and the full Rayleigh Mueller matrix is
//! already divided by its angular proposal. No extra extinction, albedo, phase
//! density, geometry Jacobian or index factor is inserted at a collision.
//!
//! The reciprocal Rayleigh/depolarized-Rayleigh matrix obeys reverse-transpose
//! symmetry with the native Stokes sign conventions. Keeping every intervening
//! Q/U/V component therefore retains reciprocal scalar M00 for unpolarized
//! emission and intensity-only reception, including multiple scattering. The
//! actual endpoint (n_camera/n_source)² factor is applied once by the existing
//! camera consumer. This is not a general polarized camera adjoint.
//!
//! References and endpoint radiometry are those of [`super::CurvedCameraScene`]; the molecular
//! matrix is owned by [`crate::volumes::transport::scattering::polarized`].
//! Only its admitted molecular phase laws are supported. No thermal volume
//! emission, scalar particle-to-Mueller inference, or moving field is added.
//!
//! Original volume, distance, integration and geometry budgets remain binding.
//! A positive global majorant does not prove vacuum escape when one density
//! sample is zero, including outside a finite grid. Only the native `Escaped`
//! termination admits environment evaluation; unresolved tails are errors.

use super::{CameraRadianceError, CameraRadianceSample, CurvedCameraScene, phase_index};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::SceneError;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighPath, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use fs_exec::Cx;

/// One completed camera sample and the original molecular traversal counters.
/// Rejected delta-tracking candidates consume free-flight draws but do not
/// increment real volume events. No partially resolved sample is returned.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighCameraSample {
    radiance: CameraRadianceSample,
    volume_events: u32,
    free_flight_draws: u32,
}

impl RayleighCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        &self.radiance
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
}

/// Borrowed composition of previously admitted emitters and molecular media.
/// Both must reference the exact same original surface scene; matching object
/// IDs or coincident geometry alone do not establish optical/source identity.
/// Preparation and sampling introduce no additional heap workspace.
pub struct CurvedRayleighCameraScene<'a, 's, 'm, 'o> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    scattering: &'a CurvedRayleighScene<'s, 'm, 'o>,
}

impl<'a, 's, 'm, 'o> CurvedRayleighCameraScene<'a, 's, 'm, 'o> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
        scattering: &'a CurvedRayleighScene<'s, 'm, 'o>,
    ) -> Result<Self, CameraRadianceError> {
        cx.checkpoint()?;
        if !core::ptr::eq(camera.source_scene(), scattering.surface_scene()) {
            return Err(CameraRadianceError::Scene(
                SceneError::InvalidBindings.into(),
            ));
        }
        Ok(Self { camera, scattering })
    }

    pub const fn surface_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }

    pub const fn participating_scene(&self) -> &'a CurvedRayleighScene<'s, 'm, 'o> {
        self.scattering
    }

    /// Retain the native independent surface/flight counter-based streams.
    /// Zero scattering consumes no flight draw and preserves surface-only
    /// results for the same seed/path identity, including every retained bit.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let path = self
            .scattering
            .for_camera()
            .trace(cx, beam, initial_medium, settings, seed, path_id)?;
        self.finish(cx, path, initial_index)
    }

    /// Uniform-density custom samples retain the native three-variate flight
    /// convention [distance CDF, phase cosine, phase azimuth]. Heterogeneous
    /// media require the tracking method's independent acceptance variate.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let path = self.scattering.for_camera().trace_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
        )?;
        self.finish(cx, path, initial_index)
    }

    /// Uniform or heterogeneous density with native independent variates:
    /// [majorant CDF, phase cosine, phase azimuth, delta acceptance]. Rejected
    /// candidates retain the same curved position, frame, power and medium.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let path = self.scattering.for_camera().trace_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
        )?;
        self.finish(cx, path, initial_index)
    }

    fn finish(
        &self,
        cx: &Cx<'_>,
        path: CurvedRayleighPath,
        initial_index: f64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        let radiance = self.camera.finish(cx, *path.path(), initial_index)?;
        Ok(RayleighCameraSample {
            radiance,
            volume_events: path.volume_events(),
            free_flight_draws: path.free_flight_draws(),
        })
    }
}

#[cfg(test)]
mod tests;

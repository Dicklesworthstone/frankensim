//! Reciprocal camera paths with all Smith/GGX conductor reflections.
//!
//! This is a path-only estimator: the original camera follows the complete
//! sampled microsurface path and observes actual surface/environment sources.
//! Reciprocity of the Smith model permits the same geometric walk from the
//! detector; each facet applies the original complex Fresnel/Mueller response.
//! The endpoint contracts the original source polarization, and applies the
//! local radiance n² conversion exactly once. No marginal directional PDF is
//! needed and the joint hidden-path density is NEVER used in MIS.
//!
//! `RoughConductor` and `BrushedConductor` select this law. Brushing retains
//! the original body-locked anisotropic frame, mark and geometric-horizon clip;
//! authored normal tilt remains an empirical shading approximation. Smooth
//! mirrors, dielectrics and diffuse surfaces keep their original laws.
//! The borrowed camera retains its source admission, physical geometry,
//! GRIN/Beer propagation, optical lengths, thermal emission and cancellation.
//! Microscopic statistical heights add no invented metric path or delay.
//!
//! Each geometric encounter has a finite facet budget. Exhaustion refuses the
//! complete sample, rather than returning black or accepting a truncated path.
//! This does not supply area-light/environment NEE, BDPT, molecular scattering,
//! or an unbiased early-stopping rule.
//! Small-source convergence remains that of a path-only estimator.
//!
//! Smith model and reciprocity: Heitz et al. (2016),
//! <https://doi.org/10.1145/2897824.2925943>.

use super::super::{CurvedSceneError, CurvedTransportEvent, SmithSurfaceSampler};
use super::emission::SurfaceEmission;
use super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedCameraScene,
    CurvedTraceSettings, phase_index,
};
use crate::camera::SequentialCameraBeam;
use crate::charts::Ray;
use crate::polarization::dielectric::scene::smith_draw;
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Explicit full rough-glass continuation, separate from opaque-light estimators.
pub mod dielectric;
pub use dielectric::SmithDielectricCameraScene;

/// Allocation-free opt-in view of an admitted physical source camera. Keeping
/// a distinct owner prevents a single-scattering NEE strategy from silently
/// receiving a different BSDF or an invented outgoing proposal probability.
pub struct SmithCameraScene<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    microsurface: SmithWalkSettings,
}

impl<'a, 's, 'm> SmithCameraScene<'a, 's, 'm> {
    pub const fn new(
        camera: &'a CurvedCameraScene<'s, 'm>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            camera,
            microsurface,
        }
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }

    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.microsurface
    }

    /// Trace one spectral radiance sample with independent counter streams for
    /// geometric surface decisions and microscopic free flights/facets.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Install the complete passive detector row before every microscopic
    /// Mueller event; there is no post-hoc scalar analyzer approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Continue the original resolved lens's reciprocal Stokes seed and local
    /// exterior index. The original raster applies its pupil geometry weight.
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
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray(cx, beam, initial_medium, settings, seed, path_id)
    }

    fn sample_launched_ray(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        self.sample_launched_ray_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                // Retain the original curved surface domain exactly.
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                    [seed as u32, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
        )
    }

    /// Ordinary surface variates remain [facet_u, facet_v, branch].
    /// Microscopic variates are independently indexed by (surface, free-flight)
    /// and follow [flight in (0,1), facet_u in [0,1), facet_v in [0,1)].
    /// A variable number of facets cannot advance a later surface's stream.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
        )
    }

    fn sample_launched_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        self.sample_launched_ray_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_with_microsurfaces(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        dielectric_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        self.sample_launched_with_emission(
            cx, beam, initial_medium, settings, surface_samples, facet_samples,
            dielectric_samples, |_, _| Ok(0.0),
        )
    }

    /// Explicit source adapters can add independently sampled physical radiance
    /// at original arrivals without duplicating the path or changing its law.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_with_emission(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
        mut extra_emission: impl FnMut(
            &super::super::CurvedPathResult,
            f64,
        ) -> Result<f64, CameraRadianceError>,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let result = self
            .camera
            .scene
            .trace_with_flights_events_and_smith(
                cx,
                beam,
                initial_medium,
                settings,
                surface_samples,
                true,
                &self.camera.spatial_volume_emitters,
                |event| match event {
                    CurvedTransportEvent::Surface(super::super::CurvedTraceEvent::Arrival(path)) => {
                        emission.observe(path)?;
                        emission.add_radiance(extra_emission(path, initial_index))
                    }
                    CurvedTransportEvent::Surface(event) => emission.observe_event(event),
                    CurvedTransportEvent::MicrosurfaceDeparture(path) => {
                        // LTE source accumulation needs the new Mueller prefix even
                        // though this departure has no marginal directional PDF.
                        emission.depart(path);
                        Ok(())
                    }
                },
                0,
                |_, _, _| Ok(None),
                |_, _, _, _| Err(CurvedSceneError::NumericalRange),
                |_, _, _| Ok(()),
                Some(SmithSurfaceSampler {
                    settings: self.microsurface,
                    samples: &mut facet_samples,
                    dielectric_samples: dielectric_samples.as_mut().map(|samples| {
                        &mut **samples as &mut dyn FnMut(u32, u32) -> [f64; 4]
                    }),
                }),
            )
            .map(|(path, _)| path);
        let (path, emitted) = emission.finish(result)?;
        self.camera.finish(cx, path, initial_index, emitted)
    }
}

#[cfg(test)]
mod tests;

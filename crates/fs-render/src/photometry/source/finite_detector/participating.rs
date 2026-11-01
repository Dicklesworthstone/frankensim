//! Calibrated finite receivers in refracting, Rayleigh-scattering media.
//!
//! The original source sampler, spectral packet calibration, native detector
//! contacts and statistical accumulators are shared with surface-only optical
//! measurements. The curved participating scene adds collisions along physical
//! arc length, including inside non-unit-index glass and gradient-index media.
//! Its phase Mueller response and native absorption already weight each packet;
//! this consumer adds no extra extinction, cosine or refractive-index factor.
//!
//! One source sample shares source direction, surface draws and volume-flight
//! draws across all 80 wavelengths. Each wavelength uses its own index and
//! scattering coefficient. Broadband moments retain the resulting correlations.
//! A surviving volume/surface/distance-limited path refuses the whole result.
//! Finite sample errors remain empirical estimates, not completeness bounds.

use super::{
    FiniteDetectorError, FiniteDetectorRequest, FiniteDetectorSpectrum, PhotometricSource,
    PreparedSource, SpectrumAccumulator,
    image::{DetectorImage, render_with_transport},
};
use crate::polarization::dielectric::scene::MAX_POLARIZED_PATH_EVENTS;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::scene::measurement::raster::DetectorRaster;
use crate::polarization::{Stokes, photometry::PolarizedSource};
use fs_exec::Cx;

/// Source/detector sampling controls and an independent per-path collision cap.
/// The scene's medium bindings determine where scattering is physically active.
#[derive(Clone, Copy, Debug)]
pub struct RayleighDetectorRequest {
    pub detector: FiniteDetectorRequest,
    pub max_volume_events: u32,
}

impl RayleighDetectorRequest {
    fn trace_settings(self) -> Result<CurvedRayleighTraceSettings, FiniteDetectorError> {
        if self.max_volume_events == 0 || self.max_volume_events > MAX_POLARIZED_PATH_EVENTS {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        Ok(CurvedRayleighTraceSettings {
            surfaces: self.detector.trace,
            max_volume_events: self.max_volume_events,
        })
    }
}

impl PhotometricSource {
    /// Measure spectral Stokes watts and lumens after native curved scattering
    /// and surface interactions. Escapes and blocked launches remain in the
    /// sample denominator; unresolved surviving paths are never counted as dark.
    pub fn measure_rayleigh_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighScene<'_, '_, '_>,
        request: RayleighDetectorRequest,
    ) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
        measure(self, cx, scene, request, |_| Stokes::unpolarized(1.0))
    }

    /// Rasterize actual detector contacts after refractive and volume transport.
    /// The result retains physical Stokes/XYZ units, crop conservation, sampling
    /// errors and the existing raw OpenEXR export.
    pub fn render_rayleigh_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighScene<'_, '_, '_>,
        request: RayleighDetectorRequest,
        raster: DetectorRaster,
    ) -> Result<DetectorImage, FiniteDetectorError> {
        render(self, cx, scene, request, raster, |_| {
            Stokes::unpolarized(1.0)
        })
    }
}

impl PolarizedSource {
    /// Retain the source's wavelength-dependent surviving polarization and
    /// upstream filtering through multiple native molecular scattering events.
    pub fn measure_rayleigh_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighScene<'_, '_, '_>,
        request: RayleighDetectorRequest,
    ) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
        measure(self.intensity_source(), cx, scene, request, |bin| {
            source_state(self, bin)
        })
    }

    /// Form a calibrated physical sensor image from the same surviving source
    /// spectrum, with correlated spectral and pixel uncertainty treatment.
    pub fn render_rayleigh_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighScene<'_, '_, '_>,
        request: RayleighDetectorRequest,
        raster: DetectorRaster,
    ) -> Result<DetectorImage, FiniteDetectorError> {
        render(self.intensity_source(), cx, scene, request, raster, |bin| {
            source_state(self, bin)
        })
    }
}

fn source_state(source: &PolarizedSource, bin: usize) -> Stokes {
    let state = source.relative_stokes()[bin];
    if state.i == 0.0 {
        Stokes::unpolarized(1.0)
    } else {
        Stokes {
            i: 1.0,
            q: state.q / state.i,
            u: state.u / state.i,
            v: state.v / state.i,
        }
    }
}

fn measure(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedRayleighScene<'_, '_, '_>,
    request: RayleighDetectorRequest,
    polarization: impl Fn(usize) -> Stokes,
) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
    let settings = request.trace_settings()?;
    let detector = request.detector;
    let surfaces = scene.surface_scene();
    let prepared = PreparedSource::new(source, cx, surfaces, detector, polarization)?;
    let mut accumulator = SpectrumAccumulator::new();
    for sample in 0..detector.samples {
        let meter = prepared.trace_sample_with(
            cx,
            surfaces,
            sample,
            |beam, path_id| {
                scene
                    .trace(
                        cx,
                        beam,
                        detector.initial_medium,
                        settings,
                        detector.seed,
                        path_id,
                    )
                    .map(|result| *result.path())
            },
            |_, _| Ok(()),
        )?;
        accumulator.add(&meter, sample)?;
    }
    cx.checkpoint()?;
    accumulator.finish(detector.samples)
}

fn render(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedRayleighScene<'_, '_, '_>,
    request: RayleighDetectorRequest,
    raster: DetectorRaster,
    polarization: impl Fn(usize) -> Stokes,
) -> Result<DetectorImage, FiniteDetectorError> {
    let settings = request.trace_settings()?;
    let detector = request.detector;
    render_with_transport(
        source,
        cx,
        scene.surface_scene(),
        detector,
        raster,
        polarization,
        |beam, path_id| {
            scene
                .trace(
                    cx,
                    beam,
                    detector.initial_medium,
                    settings,
                    detector.seed,
                    path_id,
                )
                .map(|result| *result.path())
        },
    )
}

#[cfg(test)]
mod tests;

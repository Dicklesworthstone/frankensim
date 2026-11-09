//! Steady molecular finite-area lighting through original polarized real lenses.
//!
//! All detector channels share each wavelength's native camera collisions,
//! source survival and NEE/MIS. Steady lens transport does not request group
//! timing from active coatings or nondifferentiable dispersion models.

use super::super::super::super::analyzer::{
    AnalyzerResponse, CameraAnalyzer, SpectralCameraAnalyzer,
};
use super::super::super::CurvedSequentialPixelSample;
use super::super::super::sensor::CurvedSequentialSensorPixelSample;
use super::*;
use crate::camera::SequentialCameraBeam;

impl CurvedCameraRaster<'_> {
    /// Integrate full-Stokes lens, molecular scattering and finite-area NEE/MIS paths
    /// into CIE XYZ. Active steady coatings need no group-delay model. Every
    /// clipped wavelength keeps zero contribution and its full proposal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_area_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighAreaLightPixelSample<CurvedSequentialPixelSample>, AreaLightRasterError>
    {
        self.sample_sequential_rayleigh_area_pixel(
            cx,
            scene,
            initial_medium,
            None,
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    /// Install the passive detector row at the physical sensor before the
    /// original lens, retaining its correlation with lens and scene optics.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_area_lights_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighAreaLightPixelSample<CurvedSequentialPixelSample>, AreaLightRasterError>
    {
        self.sample_sequential_rayleigh_area_pixel(
            cx,
            scene,
            initial_medium,
            Some(AnalyzerResponse::Constant(analyzer)),
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    /// Resolve the original spectral detector row over the complete visible
    /// band before tracing. Missing support refuses; it is not optical clipping.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_area_lights_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighAreaLightPixelSample<CurvedSequentialPixelSample>, AreaLightRasterError>
    {
        self.sample_sequential_rayleigh_area_pixel(
            cx,
            scene,
            initial_medium,
            Some(AnalyzerResponse::Spectral(analyzer)),
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_rayleigh_area_pixel(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighAreaLightPixelSample<CurvedSequentialPixelSample>, AreaLightRasterError>
    {
        self.sample_sequential_rayleigh_area_transport(
            scene,
            initial_medium,
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_pixel_with_timing(
                    cx,
                    seed,
                    pixel,
                    sample,
                    analyzer,
                    false,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )
    }

    /// Observe the same four optional lens/molecular/source paths through all authored
    /// QE channels. Original detector support and photon weights apply once,
    /// including nonvisible bands supported by the lens and scene sources.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_area_lights<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        AreaLightRasterError,
    > {
        self.sample_sequential_rayleigh_area_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            None,
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    /// The original QE observation with its passive analyzer before the lens,
    /// rather than an intensity factor after polarized source transport.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_area_lights_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        AreaLightRasterError,
    > {
        self.sample_sequential_rayleigh_area_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            Some(AnalyzerResponse::Constant(analyzer)),
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    /// Admit the spectral analyzer over the detector's original wavelength
    /// support, then retain each row through the complete real-lens path.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_area_lights_spectrally_analyzed<
        const C: usize,
    >(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        AreaLightRasterError,
    > {
        self.sample_sequential_rayleigh_area_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            Some(AnalyzerResponse::Spectral(analyzer)),
            settings,
            workspace,
            seed,
            pixel,
            sample,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_rayleigh_area_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighAreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        AreaLightRasterError,
    > {
        self.sample_sequential_rayleigh_area_transport(
            scene,
            initial_medium,
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_sensor_with_timing(
                    cx,
                    sensor,
                    seed,
                    pixel,
                    sample,
                    analyzer,
                    false,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_sequential_rayleigh_area_transport<P>(
        &self,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                SequentialCameraBeam,
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<RayleighAreaLightPixelSample<P>, AreaLightRasterError> {
        self.sample_rayleigh_area_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, samples| {
                scene.sample_sequential_ray_with_tracking_samples(
                    cx,
                    &beam,
                    initial_medium,
                    settings,
                    workspace,
                    |event| samples.surface(event),
                    |flight| samples.flight(flight),
                    |event, binding| samples.area(event, binding),
                    |collision, binding| samples.volume_area(collision, binding),
                    |area, image, flight| samples.source(area, image, flight),
                )
            },
            integrate,
        )
    }
}

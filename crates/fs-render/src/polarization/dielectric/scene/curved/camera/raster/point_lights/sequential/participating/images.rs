//! Steady real-lens observations of every native discovered fixture image.
use super::*;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::CurvedRayleighImageLightCameraScene;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;

impl CurvedCameraRaster<'_> {
    /// Integrate the original steady full-Stokes lens and calibrated fixture paths.
    /// Active coatings need no group-delay model; clipped lanes keep their proposal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_image_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError>
    {
        self.sample_sequential_rayleigh_image_pixel(
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

    /// Install the passive detector row before the native lens and source transport.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_image_lights_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError>
    {
        self.sample_sequential_rayleigh_image_pixel(
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

    /// Admit the spectral detector row over the original observation band, then
    /// retain its wavelength-specific polarization through the complete lens.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_image_lights_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError>
    {
        self.sample_sequential_rayleigh_image_pixel(
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
    fn sample_sequential_rayleigh_image_pixel(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError>
    {
        self.sample_sequential_rayleigh_image_transport(
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

    /// Integrate the original steady full-Stokes lens and calibrated fixture paths.
    /// Active coatings need no group-delay model; clipped lanes keep their proposal mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_image_lights<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        PointLightRasterError,
    > {
        self.sample_sequential_rayleigh_image_sensor(
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

    /// Install the passive detector row before the native lens and source transport.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_image_lights_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        PointLightRasterError,
    > {
        self.sample_sequential_rayleigh_image_sensor(
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

    /// Admit the spectral detector row over the original observation band, then
    /// retain its wavelength-specific polarization through the complete lens.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_image_lights_spectrally_analyzed<
        const C: usize,
    >(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        PointLightRasterError,
    > {
        self.sample_sequential_rayleigh_image_sensor(
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
    fn sample_sequential_rayleigh_image_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<CurvedSequentialSensorPixelSample<C>>,
        PointLightRasterError,
    > {
        self.sample_sequential_rayleigh_image_transport(
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
    fn sample_sequential_rayleigh_image_transport<P>(
        &self,
        scene: &CurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
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
    ) -> Result<RayleighPointLightPixelSample<P>, PointLightRasterError> {
        self.sample_rayleigh_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, draws| {
                scene.sample_sequential_ray_with_tracking_samples(
                    cx,
                    &beam,
                    initial_medium,
                    settings,
                    workspace,
                    |event| draws.surface(event),
                    |flight| draws.flight(flight),
                    |connection, flight| draws.source(connection, flight),
                )
            },
            integrate,
        )
    }
}

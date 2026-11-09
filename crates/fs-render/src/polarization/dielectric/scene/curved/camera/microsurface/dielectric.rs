//! Explicit path-only camera transport with two-sided rough-glass returns.
//!
//! This owner adds dielectric Smith walks to the original opaque Smith camera.
//! Ordinary point/environment/area NEE owners keep their existing dielectric
//! law: their single-scattering connection cannot be mixed with these paths.
//! The shared camera emission/endpoint code retains the reciprocal Stokes
//! launch and the single endpoint n² radiance conversion. Internal Fresnel
//! events carry power and the current half-space; no extra eta multiplier or
//! hidden-path density enters MIS. Brushed metals retain their original law.
//!
//! Sources, local GRIN indices, Beer absorption, geometry and optical distance
//! remain with the original camera. Microscopic heights have no metric length.

use super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedCameraScene,
    CurvedTraceSettings, Cx, PolarizedBeam, Ray, SequentialCameraBeam, SmithCameraScene,
    SmithWalkSettings, philox4x32_10, smith_draw,
};

/// Borrowed opt-in owner for full opaque AND dielectric Smith microsurfaces.
pub struct SmithDielectricCameraScene<'a, 's, 'm> {
    opaque: SmithCameraScene<'a, 's, 'm>,
}

impl<'a, 's, 'm> SmithDielectricCameraScene<'a, 's, 'm> {
    pub const fn new(camera: &'a CurvedCameraScene<'s, 'm>, settings: SmithWalkSettings) -> Self {
        Self {
            opaque: SmithCameraScene::new(camera, settings),
        }
    }
    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.opaque.source_camera()
    }
    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.opaque.microsurface_settings()
    }

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
        let beam = self
            .opaque
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

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
        let beam = analyzer.launch(self.opaque.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

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
        let beam = self
            .opaque
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                    [seed as u32, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [seed as u32 ^ 0x6469_656c, (seed >> 32) as u32 ^ 0x6563_7472],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
        )
    }

    /// Separate domains for ordinary surfaces, opaque microscopic facets and
    /// dielectric microscopic facets. The latter returns [flight,u,v,branch].
    /// Event/facet counters are logical identities, not a shared mutable stream.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
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
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.opaque.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<CameraRadianceSample, CameraRadianceError> {
        self.opaque.sample_launched_ray_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            Some(&mut dielectric),
        )
    }
}

//! Calibrated point fixtures on complete two-sided Smith glass paths.
//!
//! The same multiple-scattering law owns the receiver and continuation. Source
//! bundles retain their physical irradiance, native polarization and declared
//! optical branch. The local connection owns its net radiance eta factor; the
//! original camera retains its prefix conversion. No marginal PDF or MIS is
//! invented for a point source. Equal-index coated rough receivers explicitly
//! refuse the unsupported mixed-measure connection; they never use an older
//! single-scattering approximation or silently return black.

use super::*;

/// Borrowed full opaque/dielectric Smith view of an admitted point-light camera.
pub struct SmithDielectricPointLightCameraScene<'v, 'a, 's, 'm, 'l> {
    opaque: SmithPointLightCameraScene<'v, 'a, 's, 'm, 'l>,
}

impl<'v, 'a, 's, 'm, 'l> SmithDielectricPointLightCameraScene<'v, 'a, 's, 'm, 'l> {
    pub const fn new(
        source: &'v CurvedPointLightCameraScene<'a, 's, 'm, 'l>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            opaque: SmithPointLightCameraScene::new(source, microsurface),
        }
    }

    pub const fn source_point_camera(&self) -> &'v CurvedPointLightCameraScene<'a, 's, 'm, 'l> {
        self.opaque.source_point_camera()
    }

    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.opaque.microsurface
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .source
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
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(
            self.opaque.source.camera,
            ray,
            wavelength_nm,
            initial_medium,
        )?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// The original resolved lens owns its pupil weight and reciprocal seed.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .source
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
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
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
            |event, facet| draw(seed, path_id, event, facet, 0x7067_6c63, 0),
            |event, binding, facet| {
                let block = draw(seed, path_id, event, facet, 0x706d_6e65, binding as u32);
                [block[0], block[1], block[2]]
            },
            |event, binding, facet| draw(seed, path_id, event, facet, 0x7067_6e65, binding as u32),
        )
    }

    /// Independent streams for ordinary surfaces, opaque and glass continuation,
    /// then opaque and glass source evaluations. All surface/facet ordinals in
    /// this API are zero-based. Glass uses [flight,u,v,Fresnel branch].
    /// Budget or numerical refusal invalidates the entire camera sample.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        opaque_connections: impl FnMut(u32, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, u32) -> [f64; 4],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .source
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
            opaque_connections,
            dielectric_connections,
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
        settings: PointLightCameraSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        opaque_connections: impl FnMut(u32, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, u32) -> [f64; 4],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(
            self.opaque.source.camera,
            ray,
            wavelength_nm,
            initial_medium,
        )?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
            opaque_connections,
            dielectric_connections,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        opaque_connections: impl FnMut(u32, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, u32) -> [f64; 4],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .source
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
            opaque_connections,
            dielectric_connections,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric: impl FnMut(u32, u32) -> [f64; 4],
        mut opaque_connections: impl FnMut(u32, usize, u32) -> [f64; 3],
        mut dielectric_connections: impl FnMut(u32, usize, u32) -> [f64; 4],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        self.opaque.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            // The older opaque wrapper's source callback uses the post-arrival
            // count. This new API exposes the same zero-based identity as its
            // ordinary, opaque-continuation and all glass callbacks.
            |event, binding, facet| opaque_connections(event - 1, binding, facet),
            Some(&mut dielectric),
            Some(&mut dielectric_connections),
        )
    }
}

fn draw(seed: u64, path: u64, event: u32, facet: u32, domain: u32, binding: u32) -> [f64; 4] {
    let block = philox4x32_10(
        [path as u32, (path >> 32) as u32, event, facet],
        [(seed as u32) ^ domain, ((seed >> 32) as u32) ^ binding],
    );
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}

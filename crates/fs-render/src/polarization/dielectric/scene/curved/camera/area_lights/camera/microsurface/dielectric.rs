//! Complete polarized Smith glass under native finite-area illumination.
//!
//! A source may meet either adjacent half-space of the rough receiver. The
//! original source point, side/face selection, polarization channels, Beer
//! attenuation and native direct visibility retain their original owners.
//! A successful source connection uses the full two-sided dielectric response;
//! only a directly represented emitter hit after the same receiver is removed
//! from continuation. The final, not arriving, medium chooses that source set.
//!
//! Equal-index rough glass has a straight transmission atom. Its entire mixed
//! law remains on full path continuation with endpoint weight one and no source
//! draws; this owner never treats that atom as a solid-angle density. Smooth
//! interfaces and intervening scattering also keep their original endpoint term.
//! No stochastic BSDF value or hidden-path density is used in nonlinear MIS.

use super::*;

/// Borrow the original finite-source camera with full opaque and dielectric
/// Smith continuation. Each bound receiver is rough glass or an opaque Smith
/// conductor, and both receiver-adjacent media are homogeneous. Source chains
/// through other surfaces and GRIN source legs remain outside this direct
/// proposal. Its finite facet and native source/query budgets refuse unfinished
/// samples; never average only the histories that complete.
pub struct SmithDielectricAreaLightCameraScene<'v, 'a, 's, 'm, 'l> {
    shared: SmithAreaLightCameraScene<'v, 'a, 's, 'm, 'l>,
}

impl<'v, 'a, 's, 'm, 'l> SmithDielectricAreaLightCameraScene<'v, 'a, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        source: &'v CurvedAreaLightCameraScene<'a, 's, 'm, 'l>,
        microsurface: SmithWalkSettings,
    ) -> Result<Self, AreaLightCameraError> {
        cx.checkpoint()?;
        for (index, binding) in source.bindings.iter().enumerate() {
            cx.checkpoint()?;
            let invalid = || AreaLightCameraError::InvalidBinding { index };
            let scene = source.camera.scene;
            let slot = scene
                .bindings
                .binary_search_by_key(&binding.receiver_object_id, |native| native.object_id)
                .map_err(|_| invalid())?;
            let adjacent = match scene.bindings[slot].interaction {
                Interaction::RoughConductor {
                    medium, surface, ..
                } if !surface.is_delta() => [medium, medium],
                Interaction::BrushedConductor { medium, .. } => [medium, medium],
                Interaction::Dielectric {
                    outside,
                    inside,
                    surface,
                } if !surface.is_delta() => [outside, inside],
                _ => return Err(invalid()),
            };
            if !adjacent.contains(&binding.initial_medium)
                || adjacent.into_iter().any(|medium| {
                    !matches!(scene.media.get(medium), Some(SceneMedium::Homogeneous(_)))
                })
            {
                return Err(invalid());
            }
        }
        Ok(Self {
            shared: SmithAreaLightCameraScene {
                source,
                microsurface,
            },
        })
    }

    pub const fn source_area_camera(&self) -> &'v CurvedAreaLightCameraScene<'a, 's, 'm, 'l> {
        self.shared.source
    }
    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.shared.microsurface
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .shared
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(
            self.shared.source.camera,
            ray,
            wavelength_nm,
            initial_medium,
        )?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    /// Preserve the lens's original reciprocal Stokes seed. Its raster owner
    /// retains the separate pupil/geometric weight.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .shared
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
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
            |event, binding| {
                let block = philox4x32_10(
                    [
                        path_id as u32,
                        (path_id >> 32) as u32,
                        event,
                        binding as u32,
                    ],
                    [
                        (seed as u32) ^ 0x6172_6561,
                        ((seed >> 32) as u32) ^ 0x736d_6974,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, binding, channel, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [
                        (seed as u32) ^ 0x7361_7265,
                        ((seed >> 32) as u32) ^ (2 * binding + channel) as u32,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, binding, channel, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [
                        (seed as u32) ^ 0x6761_7265,
                        ((seed >> 32) as u32) ^ 0x6167_6c73 ^ (2 * binding + channel) as u32,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
        )
    }

    /// Independent streams for ordinary surfaces, opaque continuation,
    /// dielectric continuation [flight,u,v,branch], area points, opaque
    /// connections, and dielectric connections [flight,u,v,branch].
    /// Surface and facet ordinals start at zero. Each connection stream also
    /// carries its source binding and thermal polarization channel (0 or 1).
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        area: impl FnMut(u32, usize) -> [f64; 2],
        opaque_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 4],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .shared
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surfaces,
            opaque,
            dielectric,
            area,
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
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        area: impl FnMut(u32, usize) -> [f64; 2],
        opaque_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 4],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(
            self.shared.source.camera,
            ray,
            wavelength_nm,
            initial_medium,
        )?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surfaces,
            opaque,
            dielectric,
            area,
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
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        area: impl FnMut(u32, usize) -> [f64; 2],
        opaque_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        dielectric_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 4],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .shared
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surfaces,
            opaque,
            dielectric,
            area,
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
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric: impl FnMut(u32, u32) -> [f64; 4],
        area: impl FnMut(u32, usize) -> [f64; 2],
        opaque_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        mut dielectric_connections: impl FnMut(u32, usize, usize, u32) -> [f64; 4],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.shared.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surfaces,
            opaque,
            area,
            opaque_connections,
            Some(&mut dielectric),
            Some(&mut dielectric_connections),
        )
    }
}

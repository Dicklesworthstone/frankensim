//! Two-sided full Smith glass with environment next-event estimation.
//!
//! The original single-facet camera proposal is used only as a deterministic
//! angular guide q. Light and escaped-camera terms use p²/(p²+q²) and
//! q²/(p²+q²), respectively. The full-walk marginal PDF is neither known nor
//! substituted by q. The stochastic evaluated response appears only linearly;
//! no sampled facet or response changes these complementary weights.
//!
//! Rough equal-index interfaces remain complete path-only transport. A coating
//! can give them a straight-through atom alongside continuous reflection, so
//! treating their whole walk as an ordinary solid-angle proposal would lose
//! energy. Smooth surfaces, intervening contacts and unsupported infinity media
//! likewise retain the original endpoint rule. Each selected source direction
//! checks its own adjacent medium before consuming connection-facet samples.

use super::*;
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_geom::Vec3;

/// Borrowed opt-in view with complete opaque and two-sided dielectric Smith
/// transport. Native environment calibration, source polarization, visibility,
/// medium topology, surface/LTE emission and camera index conversion are shared
/// with the original camera owner. No synthetic source or receiving surface is
/// introduced.
pub struct SmithDielectricEnvironmentCameraScene<'v, 'a, 's, 'm> {
    opaque: SmithEnvironmentCameraScene<'v, 'a, 's, 'm>,
}

impl<'v, 'a, 's, 'm> SmithDielectricEnvironmentCameraScene<'v, 'a, 's, 'm> {
    pub const fn new(
        source: &'v CurvedEnvironmentCameraScene<'a, 's, 'm>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            opaque: SmithEnvironmentCameraScene::new(source, microsurface),
        }
    }

    pub const fn source_environment_camera(&self) -> &'v CurvedEnvironmentCameraScene<'a, 's, 'm> {
        self.opaque.source
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
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
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
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(
            self.opaque.source.camera,
            ray,
            wavelength_nm,
            initial_medium,
        )?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        let beam = self
            .opaque
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| counter_samples(seed, path_id, event, 0x706f_6c32),
            |event| {
                let sample = counter_samples(seed, path_id, event, 0x706f_6c39);
                [sample[0], sample[1]]
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [(seed as u32) ^ 0x736d_656e, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, facet| glass_draw(seed, path_id, event, facet, 0x6765_636f),
            |event, facet| glass_draw(seed, path_id, event, facet, 0x6765_6e65),
        )
    }

    /// Six independent streams: ordinary surface, environment direction,
    /// opaque continuation, opaque connection, dielectric continuation and
    /// dielectric connection. Surface and facet ordinals are zero based.
    /// Every glass quartet is [flight, facet u, facet v, Fresnel branch].
    /// A bounded-walk failure refuses the entire estimate; successful-history
    /// averaging after rejecting a failure is not a valid estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
        dielectric_samples: impl FnMut(u32, u32) -> [f64; 4],
        dielectric_connection_samples: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
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
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
            dielectric_samples,
            dielectric_connection_samples,
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
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
        dielectric_samples: impl FnMut(u32, u32) -> [f64; 4],
        dielectric_connection_samples: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
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
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
            dielectric_samples,
            dielectric_connection_samples,
        )
    }

    /// The resolved real-lens beam owns its reciprocal optical response. Its
    /// original pupil/geometric raster weight remains with the raster caller.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
        dielectric_samples: impl FnMut(u32, u32) -> [f64; 4],
        dielectric_connection_samples: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
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
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
            dielectric_samples,
            dielectric_connection_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        environment_samples: impl FnMut(u32) -> [f64; 2],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric_samples: impl FnMut(u32, u32) -> [f64; 4],
        mut dielectric_connection_samples: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<EnvironmentCameraSample, EnvironmentCameraError> {
        self.opaque.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            environment_samples,
            facet_samples,
            connection_samples,
            Some(&mut dielectric_samples),
            Some(&mut dielectric_connection_samples),
        )
    }
}

fn glass_draw(seed: u64, path_id: u64, event: u32, facet: u32, domain: u32) -> [f64; 4] {
    let block = philox4x32_10(
        [path_id as u32, (path_id >> 32) as u32, event, facet],
        [(seed as u32) ^ domain, ((seed >> 32) as u32) ^ 0x6469_656c],
    );
    block.map(|value| (f64::from(value) + 0.5) / 4_294_967_296.0)
}

// This guide uses only the authored surface, original arrival state and exact
// requested direction. Unit unpolarized source input makes the two calls
// independent of map polarization; the native proposal is polarization
// independent. Its legacy masking law is not the full Smith physical response.
pub(super) fn angular_guide(
    receiver: SurfaceReceiver,
    camera: &CurvedCameraScene<'_, '_>,
    outgoing: Vec3,
) -> Result<f64, EnvironmentCameraError> {
    let incoming = outgoing.scale(-1.0);
    let [x, y, z] = [incoming.x.abs(), incoming.y.abs(), incoming.z.abs()];
    let axis = if x <= y && x <= z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if y <= z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let frame = PolarizationFrame::try_new(incoming, axis)
        .map_err(CurvedSceneError::from)
        .map_err(CameraRadianceError::from)?;
    let (response, _) = receiver.directional_response(camera, frame, Stokes::unpolarized(1.0))?;
    if response.pdf.is_finite() && response.pdf >= 0.0 {
        Ok(response.pdf)
    } else {
        Err(EnvironmentCameraError::NumericalRange)
    }
}

pub(super) fn any_supported_medium(
    receiver: SurfaceReceiver,
    camera: &CurvedCameraScene<'_, '_>,
    wavelength_nm: f64,
) -> Result<bool, EnvironmentCameraError> {
    let media = receiver
        .smith_dielectric_media(camera)?
        .ok_or(EnvironmentCameraError::NumericalRange)?;
    Ok(supports_segment(camera.scene, media[0], wavelength_nm)?
        || supports_segment(camera.scene, media[1], wavelength_nm)?)
}

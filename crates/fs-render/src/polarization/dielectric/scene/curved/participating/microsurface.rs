//! Complete Smith rough surfaces on the original participating-medium walk.
//!
//! Physical volume free flights, null-collision acceptance, local polarized
//! phase laws, refractive fields, finite absorption cells and geometric event
//! limits are owned by the original scene. The sole opt-in is the surface law:
//! rough/brushed metal and rough glass use the existing bounded Smith walk.
//! Statistical facet height adds no physical distance, attenuation or delay.
//! Microscopic escape selects the actual final medium before the next volume
//! flight. Surface, physical-flight and facet streams have separate addresses.
//!
//! A complete walk has no supplied marginal directional density. Its hidden
//! path density is never passed off as a phase/BSDF PDF. Consumers that need NEE
//! must evaluate the matching complete surface response independently.

use super::{
    CurvedRayleighPath, CurvedRayleighScene, CurvedRayleighTraceSettings,
    CurvedRayleighTransportEvent, CurvedSceneError, PolarizedBeam, SmithSurfaceSampler,
};
use crate::polarization::dielectric::scene::smith_draw;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Borrow an admitted participating scene and select complete opaque and
/// dielectric microsurfaces. No geometry, optical law or density is copied.
pub struct SmithDielectricRayleighScene<'a, 's, 'm, 'o> {
    scene: &'a CurvedRayleighScene<'s, 'm, 'o>,
    microsurface: SmithWalkSettings,
}

impl<'a, 's, 'm, 'o> SmithDielectricRayleighScene<'a, 's, 'm, 'o> {
    pub const fn new(
        scene: &'a CurvedRayleighScene<'s, 'm, 'o>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            scene,
            microsurface,
        }
    }

    pub const fn source_scene(&self) -> &'a CurvedRayleighScene<'s, 'm, 'o> {
        self.scene
    }

    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.microsurface
    }

    /// Native source-power transport; no radiance index conversion is applied.
    /// The original forward brushing and volume Mueller conventions remain in
    /// force. A real volume cap remains an unresolved `VolumeLimit` path.
    #[allow(clippy::too_many_arguments)]
    pub fn trace(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            &[],
            |_| Ok(()),
        )
    }

    /// Uniform and heterogeneous media share native flight draws
    /// `[majorant CDF, phase cosine, phase azimuth, acceptance]`. Facet draws
    /// are `[flight,u,v]` for opaque and `[flight,u,v,branch]` for dielectric
    /// surfaces. Every ordinal starts at zero; rejected volume candidates use
    /// a flight draw but consume neither a surface nor a microscopic draw.
    #[allow(clippy::too_many_arguments)]
    pub fn trace_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.trace_with_tracking_samples_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            opaque,
            dielectric,
            &[],
            |_| Ok(()),
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn trace_events_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
        thermal_sources: &super::super::SpatialSources<'_>,
        observe: impl FnMut(CurvedRayleighTransportEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        let draws = |ordinal, domain| {
            philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, ordinal, domain],
                [seed as u32, (seed >> 32) as u32],
            )
            .map(|lane| (f64::from(lane) + 0.5) / 4_294_967_296.0)
        };
        self.trace_with_tracking_samples_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                let draw = draws(event, 0x706f_6c32);
                [draw[0], draw[1], draw[2]]
            },
            |flight| draws(flight, 0x706f_6c34),
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, facet| {
                philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [seed as u32 ^ 0x6469_656c, (seed >> 32) as u32 ^ 0x6563_7472],
                )
                .map(|lane| (f64::from(lane) + 0.5) / 4_294_967_296.0)
            },
            thermal_sources,
            observe,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn trace_with_tracking_samples_events_observed(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        mut opaque: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric: impl FnMut(u32, u32) -> [f64; 4],
        thermal_sources: &super::super::SpatialSources<'_>,
        observe: impl FnMut(CurvedRayleighTransportEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedRayleighPath, CurvedSceneError> {
        self.scene.trace_with_tracking_samples_events_and_smith(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            true,
            thermal_sources,
            observe,
            Some(SmithSurfaceSampler {
                settings: self.microsurface,
                samples: &mut opaque,
                dielectric_samples: Some(&mut dielectric),
            }),
        )
    }
}

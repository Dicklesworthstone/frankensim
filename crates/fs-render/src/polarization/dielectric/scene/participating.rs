//! Molecular scattering in the SAME native surface-to-detector traversal.
//!
//! One optional original `RayleighGrid` belongs to each explicit scene medium.
//! Its density/cell optics apply only while that medium is current. Thus air
//! inside a grid does not also scatter inside an overlapping glass object.
//! A bound carrier must be homogeneous, clear and index-one at the transported
//! wavelength: the grid owns extinction. Unbound media retain native absorption.
//! This straight-segment, source-power model does not add LTE emission, bent
//! rays within a scattering medium, particle polarization or optical heat.

use super::{Cx, PathResult, PolarizedBeam, PolarizedScene, SceneError, SegmentProgress,
    TraceSettings, MAX_POLARIZED_PATH_EVENTS, philox4x32_10};
use crate::volumes::transport::scattering::polarized::{
    PolarizedScatteringError, RayleighGrid, validate_carrier,
};
use crate::volumes::transport::VolumeTransportError;

impl From<PolarizedScatteringError> for SceneError {
    fn from(error: PolarizedScatteringError) -> Self {
        match error {
            PolarizedScatteringError::Volume(VolumeTransportError::Cancelled) => Self::Cancelled,
            other => Self::Scattering(other),
        }
    }
}

/// Independent physical-distance, surface-event and volume-collision limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighTraceSettings {
    pub surfaces: TraceSettings,
    pub max_volume_events: u32,
}

/// Original surface result, accepted by its spectral meter without an alternate
/// estimator or extra attenuation. A volume cap appears as `VolumeLimit` in path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighPath {
    path: PathResult,
    volume_events: u32,
    free_flight_draws: u32,
}
impl RayleighPath {
    pub const fn path(&self) -> &PathResult { &self.path }
    pub const fn volume_events(&self) -> u32 { self.volume_events }
    pub const fn free_flight_draws(&self) -> u32 { self.free_flight_draws }
}

/// Immutable borrowed binding; no geometry, material tables or density arrays
/// are copied. `volumes[i]` corresponds to the original scene's medium index i.
/// Preparation and tracing allocate no heap workspace.
pub struct RayleighScene<'a, 'scene, 'grid, 'data> {
    scene: &'a PolarizedScene<'scene>,
    volumes: &'a [Option<RayleighGrid<'grid, 'data>>],
}
impl<'a, 'scene, 'grid, 'data> RayleighScene<'a, 'scene, 'grid, 'data> {
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'a PolarizedScene<'scene>,
        volumes: &'a [Option<RayleighGrid<'grid, 'data>>],
    ) -> Result<Self, SceneError> {
        cx.checkpoint()?;
        if volumes.len() != scene.media.len() { return Err(SceneError::InvalidBindings); }
        Ok(Self { scene, volumes })
    }

    pub const fn surface_scene(&self) -> &'a PolarizedScene<'scene> { self.scene }

    /// Surface draws retain the original stream. Free flights use a separate
    /// domain and ordinal, including escaped segments. Replay is worker-order
    /// independent; repeated path IDs intentionally correlate wavelengths.
    #[allow(clippy::too_many_arguments)]
    pub fn trace(
        &self, cx: &Cx<'_>, beam: PolarizedBeam, initial_medium: usize,
        settings: RayleighTraceSettings, seed: u64, path_id: u64,
    ) -> Result<RayleighPath, SceneError> {
        let draws = |ordinal, domain| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, ordinal, domain],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        };
        self.trace_with_samples(cx, beam, initial_medium, settings,
            |event| draws(event, 0x706f_6c31), |flight| draws(flight, 0x706f_6c33))
    }

    /// Explicit independent sample streams. Surface draws have the original
    /// [facet_u, facet_v, event] semantics. Flight draws are [survival, cosine,
    /// azimuth], all in [0,1). Absorption-only grids consume no random draw.
    /// No phase, absorption, or geometry failure becomes a dark sample.
    #[allow(clippy::too_many_arguments)]
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: RayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<RayleighPath, SceneError> {
        cx.checkpoint()?;
        if settings.max_volume_events == 0 || settings.max_volume_events > MAX_POLARIZED_PATH_EVENTS {
            return Err(SceneError::InvalidSettings);
        }
        let mut volume_events = 0;
        let mut free_flight_draws = 0_u32;
        let path = self.scene.trace_segments(
            cx, beam, initial_medium, settings.surfaces, surface_samples,
            |beam, medium, distance| {
                let Some(grid) = self.volumes[medium] else {
                    return SegmentProgress::homogeneous(beam, distance);
                };
                // Refuse duplicate absorption/refraction before consuming the
                // caller's stream. These data were not inferred from the grid.
                validate_carrier(beam)?;
                let draw = if grid.grid().has_scattering() {
                    let draw = flight_samples(free_flight_draws);
                    free_flight_draws += 1;
                    draw
                } else { [0.5; 3] };
                let flight = grid.sample_beam_segment(cx, beam, distance, draw)?;
                let collision = flight.event.is_some();
                if collision { volume_events += 1; }
                Ok(SegmentProgress {
                    beam: flight.beam,
                    distance_m: flight.distance_m,
                    collision,
                    volume_limit: collision && volume_events >= settings.max_volume_events,
                })
            },
        )?;
        cx.checkpoint()?;
        Ok(RayleighPath { path, volume_events, free_flight_draws })
    }
}

#[cfg(test)]
mod tests;

//! Material exchange at actual smooth-mirror hits in the native scene path.

use super::{Cx, InstanceHit, PathResult, PolarizedBeam, PolarizedScene, SceneError, SegmentProgress, TraceSettings};
use super::super::mirrors::ConductorBeamBalance;

/// One actual native smooth-conductor hit and its optical/material exchange.
/// The hit retains the original object, geometry, frame and local coordinates
/// for assigning heat and force to the corresponding physical surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MirrorInteraction {
    pub hit: InstanceHit,
    pub balance: ConductorBeamBalance,
}

/// The original path plus its ordered smooth-mirror interactions.
///
/// This is not a whole-scene deposition report: bulk attenuation, dielectric,
/// diffuse, rough-metal, detector and opaque-termination exchanges are not
/// counted here. A capped path retains its surviving beam in `path`.
#[derive(Clone, Debug, PartialEq)]
pub struct MirrorPathResult {
    pub path: PathResult,
    pub mirrors: Vec<MirrorInteraction>,
}

impl PolarizedScene<'_> {
    /// Trace the original path and retain heat/force at every smooth metal hit.
    ///
    /// Uses the same `(seed, path_id, event)` draws, geometry, media, query
    /// offsets and reflected beam as `trace`. Exchange uses the actual beam
    /// after propagation to each hit, including absorption before the mirror
    /// and all earlier reflections. No extra cosine or roulette is applied.
    /// Smooth-mirror heat comes from native positive electromagnetic flux;
    /// force is available only in index-one hosts. The caller owns thermal and
    /// mechanical integration and wavelength/path ensemble normalization.
    ///
    /// Storage grows only with encountered smooth-mirror hits, bounded by the
    /// original event budget. A failure returns no partial report. Ordinary
    /// `trace` calls allocate no exchange storage and solve no absorption flux.
    pub fn trace_with_mirror_balances(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<MirrorPathResult, SceneError> {
        let mut mirrors = Vec::new();
        let path = self.trace_segments(
            cx,
            beam,
            initial_medium,
            settings,
            |event| super::path_draw(seed, path_id, event),
            |beam, _, distance| SegmentProgress::homogeneous(beam, distance),
            Some(&mut mirrors),
        )?;
        Ok(MirrorPathResult { path, mirrors })
    }
}

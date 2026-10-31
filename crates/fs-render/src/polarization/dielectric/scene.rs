//! Bounded forward polarized transport over the native immutable instance scene.
//!
//! Geometry, intersection authority, GGX sampling and material attenuation keep
//! their existing owners. Explicit media on both sides of each oriented surface
//! handle nested boundaries without guessing containment or cloning a medium
//! stack per ray. The encountered side must agree with the current medium ID.
//!
//! This is an opt-in light-path API, not a camera integrator or a topology proof.
//! A positive query offset suppresses self-hits: features/gaps below that offset
//! are unresolved. The offset is NOT subtracted from physical path length or
//! absorption. Detectors measure incident ray-power weights (no extra cosine or
//! inverse-square factor). Monte Carlo averages and spectral integration belong
//! to the caller. Limits return surviving beams, never invented absorption.

use fs_exec::Cx;
use fs_geom::Vec3;
use fs_rand::philox::philox4x32_10;

use super::transport::{BeamError, BeamMedium, PolarizedBeam};
use crate::dielectric::{DielectricEvent, DielectricGlass, DielectricSurface};
use crate::instances::{InstanceError, InstanceHit, InstanceScene};

/// Admission bound on scene bindings and media; not a performance claim.
pub const MAX_POLARIZED_SCENE_ITEMS: usize = 65_536;
/// Maximum admitted surface-event budget for one synchronous path.
pub const MAX_POLARIZED_PATH_EVENTS: u32 = 65_536;

/// Material or measurement attached to an existing stable geometry ID.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Interaction {
    /// The geometric normal points from `inside` to `outside`. Both indices
    /// refer to the scene's original media slice; equal-index glass may still
    /// have different absorption. Misoriented/overlapping media fail closed.
    Dielectric { outside: usize, inside: usize, surface: DielectricSurface },
    /// A two-sided ideal receiver. Positive Q uses the world reference axis
    /// projected transverse to each arriving ray. A parallel axis is refused.
    /// No lens/pupil transformation or grazing-angle polarizer model is inferred.
    Detector { reference_axis: Vec3 },
    /// Explicitly terminate at an opaque absorber, without invented reflection.
    Opaque,
}

/// Exactly one interaction is required for every native geometry instance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceBinding {
    pub object_id: u64,
    pub interaction: Interaction,
}

/// Explicit finite path and geometry-resolution limits, in SI metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceSettings {
    pub max_surface_events: u32,
    pub max_distance_m: f64,
    pub surface_epsilon_m: f64,
    /// Must be at least four times the intersection epsilon. The caller must
    /// choose it below the smallest feature/gap that needs to be resolved.
    pub query_offset_m: f64,
}

impl TraceSettings {
    fn validate(self) -> Result<(), SceneError> {
        if self.max_surface_events == 0 || self.max_surface_events > MAX_POLARIZED_PATH_EVENTS
            || !self.max_distance_m.is_finite() || self.max_distance_m <= 0.0
            || !self.surface_epsilon_m.is_finite() || self.surface_epsilon_m <= 0.0
            || !self.query_offset_m.is_finite()
            || self.query_offset_m / self.surface_epsilon_m < 4.0
            || self.query_offset_m >= self.max_distance_m {
            return Err(SceneError::InvalidSettings);
        }
        Ok(())
    }
}

/// Failures are not dark samples. Original geometry, media and input beams
/// remain unchanged. User sampling callbacks are not transactionally rolled back.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SceneError {
    InvalidBindings,
    InvalidSettings,
    InvalidInitialMedium,
    MediumMismatch { object_id: u64, expected: usize, actual: usize },
    InvalidHit,
    UnrepresentableProgress,
    Geometry(InstanceError),
    Beam(BeamError),
    Cancelled,
}
impl From<InstanceError> for SceneError {
    fn from(error: InstanceError) -> Self {
        match error { InstanceError::Cancelled => Self::Cancelled, other => Self::Geometry(other) }
    }
}
impl From<BeamError> for SceneError {
    fn from(error: BeamError) -> Self { Self::Beam(error) }
}
impl From<fs_exec::Cancelled> for SceneError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Cancelled }
}
impl core::fmt::Display for SceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized scene: {self:?}")
    }
}
impl core::error::Error for SceneError {}

/// Why a path stopped. A limited path is not an escaped or absorbed path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathTermination {
    Detected,
    Opaque,
    /// Native single-scattering GGX rejection, not a material heat measurement.
    RejectedMicrofacet,
    Extinguished,
    DistanceLimit,
    SurfaceLimit,
}

/// Immutable path result. The final native hit retains geometry authority and
/// local surface coordinates for downstream detector/pixel mapping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathResult {
    beam: PolarizedBeam,
    medium: usize,
    traveled_m: f64,
    surface_events: u32,
    termination: PathTermination,
    last_hit: Option<InstanceHit>,
}
impl PathResult {
    pub const fn beam(&self) -> &PolarizedBeam { &self.beam }
    pub const fn medium_index(&self) -> usize { self.medium }
    pub const fn traveled_m(&self) -> f64 { self.traveled_m }
    pub const fn surface_events(&self) -> u32 { self.surface_events }
    pub const fn termination(&self) -> PathTermination { self.termination }
    pub const fn last_hit(&self) -> Option<&InstanceHit> { self.last_hit.as_ref() }
}

/// Prepared bindings borrow the actual native scene and spectral materials.
/// Admission sorts once; paths have no owned heap workspace and do not clone
/// geometry. Traversal retains InstanceScene's per-object BLAS/chart behavior.
#[derive(Debug)]
pub struct PolarizedScene<'a> {
    geometry: &'a InstanceScene,
    media: &'a [DielectricGlass],
    bindings: Vec<SurfaceBinding>,
}
impl<'a> PolarizedScene<'a> {
    pub fn try_new(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        mut bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        if media.is_empty() || media.len() > MAX_POLARIZED_SCENE_ITEMS
            || bindings.len() != geometry.instances().len()
            || bindings.len() > MAX_POLARIZED_SCENE_ITEMS {
            return Err(SceneError::InvalidBindings);
        }
        bindings.sort_unstable_by_key(|b| b.object_id);
        for (binding, instance) in bindings.iter().zip(geometry.instances()) {
            if binding.object_id != instance.object_id() { return Err(SceneError::InvalidBindings); }
            match binding.interaction {
                Interaction::Dielectric { outside, inside, .. } => {
                    if outside >= media.len() || inside >= media.len() {
                        return Err(SceneError::InvalidBindings);
                    }
                }
                Interaction::Detector { reference_axis: a } => {
                    if ![a.x, a.y, a.z].iter().all(|x| x.is_finite())
                        || a.x.abs().max(a.y.abs()).max(a.z.abs()) == 0.0 {
                        return Err(SceneError::InvalidBindings);
                    }
                }
                Interaction::Opaque => {}
            }
        }
        Ok(Self { geometry, media, bindings })
    }

    /// Deterministic counter-based sampling keyed by `(seed, path_id, event)`.
    /// This stream is separate from scalar/camera renderers and worker order.
    /// Reusing a path ID intentionally correlates wavelengths or repeated runs.
    pub fn trace(
        &self, cx: &Cx<'_>, beam: PolarizedBeam, initial_medium: usize,
        settings: TraceSettings, seed: u64, path_id: u64,
    ) -> Result<PathResult, SceneError> {
        self.trace_with_samples(cx, beam, initial_medium, settings, |event| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c31],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        })
    }

    /// Trace with externally supplied `[microfacet_u, microfacet_v, event]`
    /// samples. The callback receives the zero-based surface-event index and
    /// is called only for dielectric hits, after geometry and medium admission.
    /// This scene uses piecewise-homogeneous straight segments. A sampled
    /// gradient-index medium refuses at launch, before any traversal or sample;
    /// its curved boundary intersections are not implemented by this owner.
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<PathResult, SceneError> {
        cx.checkpoint()?;
        settings.validate()?;
        let BeamMedium::Homogeneous(material) = beam.medium() else {
            return Err(SceneError::InvalidInitialMedium);
        };
        if self.media.get(initial_medium) != Some(&material) {
            return Err(SceneError::InvalidInitialMedium);
        }
        let mut result = PathResult {
            beam, medium: initial_medium, traveled_m: 0.0, surface_events: 0,
            termination: PathTermination::SurfaceLimit, last_hit: None,
        };
        let mut offset = 0.0;
        for event in 0..settings.max_surface_events {
            cx.checkpoint()?;
            if result.beam.power().i == 0.0 {
                result.termination = PathTermination::Extinguished;
                return finish(cx, result);
            }
            let remaining = settings.max_distance_m - result.traveled_m;
            if remaining <= offset {
                return distance_limit(cx, result, remaining, settings.max_distance_m);
            }
            let mut query = result.beam.ray();
            query.origin = query.at(offset);
            if ![query.origin.x, query.origin.y, query.origin.z].iter().all(|x| x.is_finite())
                || (offset > 0.0 && query.origin == result.beam.ray().origin) {
                return Err(SceneError::UnrepresentableProgress);
            }
            let Some(hit) = self.geometry.intersect(
                cx, &query, remaining - offset, settings.surface_epsilon_m,
            )? else {
                return distance_limit(cx, result, remaining, settings.max_distance_m);
            };
            if !hit.hit.t.is_finite() || hit.hit.t <= 0.0 || hit.hit.t > remaining - offset {
                return Err(SceneError::InvalidHit);
            }
            // Native traversal already admitted the query-space endpoint.
            // Adding the prefix back can round one ulp past the physical cap.
            let distance = (offset + hit.hit.t).min(remaining);
            let traveled = result.traveled_m + distance;
            if traveled <= result.traveled_m { return Err(SceneError::UnrepresentableProgress); }
            // Advance from the PHYSICAL origin, not the displaced query origin.
            result.beam = result.beam.advance(distance)?;
            result.traveled_m = traveled.min(settings.max_distance_m);
            result.surface_events += 1;
            result.last_hit = Some(hit);
            let index = self.bindings.binary_search_by_key(&hit.object_id, |b| b.object_id)
                .map_err(|_| SceneError::InvalidBindings)?;
            match self.bindings[index].interaction {
                Interaction::Detector { reference_axis } => {
                    result.beam = result.beam.with_reference_axis(reference_axis)?;
                    result.termination = PathTermination::Detected;
                    return finish(cx, result);
                }
                Interaction::Opaque => {
                    result.termination = PathTermination::Opaque;
                    return finish(cx, result);
                }
                Interaction::Dielectric { outside, inside, surface } => {
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    let cosine = normal.dot(result.beam.ray().dir);
                    if !cosine.is_finite() || cosine == 0.0 { return Err(SceneError::InvalidHit); }
                    let (expected, target) = if cosine < 0.0 { (outside, inside) } else { (inside, outside) };
                    if result.medium != expected {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected, actual: result.medium,
                        });
                    }
                    let Some(next) = result.beam.sample_interface(
                        normal, self.media[target], surface, samples(event),
                    )? else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        return finish(cx, result);
                    };
                    if next.event == DielectricEvent::Transmission { result.medium = target; }
                    result.beam = next.beam;
                    offset = settings.query_offset_m;
                }
            }
        }
        finish(cx, result)
    }
}

fn finish(cx: &Cx<'_>, result: PathResult) -> Result<PathResult, SceneError> {
    cx.checkpoint()?;
    Ok(result)
}
fn distance_limit(
    cx: &Cx<'_>, mut result: PathResult, remaining: f64, maximum: f64,
) -> Result<PathResult, SceneError> {
    result.beam = result.beam.advance(remaining.max(0.0))?;
    result.traveled_m = maximum;
    result.termination = PathTermination::DistanceLimit;
    finish(cx, result)
}

/// Spectral receiving power and ideal analyzer readouts on native CIE bins.
#[cfg(feature = "photometry")]
pub mod measurement;

#[cfg(test)]
mod tests;

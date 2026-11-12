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

use super::mirrors::ConductorBeamError;
use super::transport::microsurface::{SmithBeamError, SmithWalkSettings};
use super::transport::{BeamError, BeamMedium, BeamSample, PolarizedBeam};
use crate::conductor::{CoherentDielectricCoating, ConductorOptics, ConductorSurface};
use crate::dielectric::{DielectricCoatingStack, DielectricEvent, DielectricGlass, DielectricSurface};
use crate::instances::{InstanceError, InstanceHit, InstanceScene, RigidTransform};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, SpectralReflectance};

/// Admission bound on scene bindings and media; not a performance claim.
pub const MAX_POLARIZED_SCENE_ITEMS: usize = 65_536;
/// Maximum admitted surface-event budget for one synchronous path.
pub const MAX_POLARIZED_PATH_EVENTS: u32 = 65_536;

fn path_draw(seed: u64, path_id: u64, event: u32) -> [f64; 3] {
    let block = philox4x32_10(
        [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c31],
        [seed as u32, (seed >> 32) as u32],
    );
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}

// The extra microscopic coordinates have their own Philox key domain. Surface
// event and facet event are separate counters; no variable-length walk shifts
// another geometric hit's original random stream.
fn smith_draw(seed: u64, path_id: u64, surface_event: u32, facet_event: u32) -> [f64; 3] {
    let block = philox4x32_10(
        [path_id as u32, (path_id >> 32) as u32, surface_event, facet_event],
        [seed as u32 ^ 0x736d_6974, (seed >> 32) as u32 ^ 0x6867_6778],
    );
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}

/// Curved paths through immutable gradient-index fields and native geometry.
pub mod curved;

mod mirror_exchange;
pub use mirror_exchange::{MirrorInteraction, MirrorPathResult};

mod coating_stacks;
pub use coating_stacks::DielectricStackBinding;
use coating_stacks::admit_coating_stacks;

/// Molecular multiple scattering between the original native surface events.
#[cfg(feature = "volumes")]
pub mod participating;

/// Material or measurement attached to an existing stable geometry ID.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Interaction {
    /// The geometric normal points from `inside` to `outside`. Both indices
    /// refer to the scene's original media slice; equal-index glass may still
    /// have different absorption. Misoriented/overlapping media fail closed.
    /// An anisotropic surface axis is authored in geometry-local coordinates
    /// and follows the instance's original proper-rigid placement.
    Dielectric { outside: usize, inside: usize, surface: DielectricSurface },
    /// Smooth dielectric boundary carrying one original transparent film.
    /// `coating` indexes the borrowed coating table. The same physical film
    /// applies in both directions, with the actual incident/exit media.
    CoatedDielectric { outside: usize, inside: usize, coating: usize },
    /// Two-sided smooth opaque reflector in the declared surrounding medium.
    /// `optics` indexes the separately borrowed conductor table. The same law
    /// (including an explicitly supplied film) applies on both sides; distinct
    /// front/back coatings require distinct geometry bindings. No roughness or
    /// transmission into the opaque substrate is inferred.
    SmoothConductor { optics: usize, medium: usize },
    /// Two-sided isotropic GGX metal, single scattering by default. The explicit
    /// Smith multiple-scattering trace follows internal facet returns instead.
    /// Optical constants and
    /// coherent coatings remain in the shared table, independent of roughness.
    /// Native rejected microfacet samples terminate without resampling. This
    /// does not model anisotropic brushing or a tilted shading normal.
    RoughConductor { optics: usize, medium: usize, surface: ConductorSurface },
    /// Native body-locked radial/circular brushing, including its shading
    /// normal and radial roughness mark. Forward power uses the adjoint
    /// shading-normal law; no resolved groove geometry or heat law is inferred.
    BrushedConductor { optics: usize, medium: usize, detail: crate::surface_detail::SurfaceDetail },
    /// Two-sided ideal depolarizing Lambertian reflection in the declared
    /// surrounding medium. `reflectance` indexes the original borrowed
    /// spectral reflectance table; no emission or subsurface law is inferred.
    DepolarizingDiffuse { medium: usize, reflectance: usize },
    /// Reciprocal EON rough diffuse reflection. The original spectral table
    /// supplies diffuse-facet albedo, with wavelength-dependent interreflection.
    RoughDiffuse { medium: usize, reflectance: usize, surface: crate::rough_diffuse::RoughDiffuseSurface },
    /// The original tracer's lifted Lambertian spectrum, evaluated directly
    /// over its native visible band. No knot fitting or RGB conversion occurs.
    /// Reflection uses the same depolarizing cosine law as tabulated diffuse.
    DepolarizingLiftedDiffuse { medium: usize, reflectance: LiftedSpectrum },
    /// Smooth lossless dielectric over the original spectral diffuse base.
    /// Its constant exterior index must match the actual local medium. Mirror
    /// and diffuse events share the same surface; internal returns are summed.
    CoatedDiffuse { medium: usize, reflectance: usize, coating: crate::coated_diffuse::SmoothDiffuseCoating },
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
    AllocationRefused,
    InvalidSettings,
    InvalidInitialMedium,
    MediumMismatch { object_id: u64, expected: usize, actual: usize },
    InvalidHit,
    UnrepresentableProgress,
    Geometry(InstanceError),
    Beam(BeamError),
    Conductor(ConductorBeamError),
    /// A microscopic path failed or exceeded its explicit facet budget.
    Smith(SmithBeamError),
    /// Explicit two-sided dielectric microsurface path failure.
    SmithDielectric(super::transport::microsurface::dielectric::SmithDielectricError),
    #[cfg(feature = "volumes")]
    Scattering(crate::volumes::transport::scattering::polarized::PolarizedScatteringError),
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
impl From<ConductorBeamError> for SceneError {
    fn from(error: ConductorBeamError) -> Self { Self::Conductor(error) }
}
impl From<SmithBeamError> for SceneError {
    fn from(error: SmithBeamError) -> Self {
        match error { SmithBeamError::Cancelled => Self::Cancelled, other => Self::Smith(other) }
    }
}
impl From<super::transport::microsurface::dielectric::SmithDielectricError> for SceneError {
    fn from(error: super::transport::microsurface::dielectric::SmithDielectricError) -> Self {
        match error {
            super::transport::microsurface::dielectric::SmithDielectricError::Cancelled => {
                Self::Cancelled
            }
            other => Self::SmithDielectric(other),
        }
    }
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
    /// A straight continuation in a homogeneous medium provably misses every
    /// native scene surface. Surviving power is measured at the retained
    /// endpoint, not after hypothetical propagation to infinity.
    Escaped,
    DistanceLimit,
    SurfaceLimit,
    /// A surviving molecular path exhausted its explicit collision budget.
    VolumeLimit,
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

// Internal segment seam: native surface traversal/dispatch remains one owner.
// A segment either reaches its endpoint or supplies one volume continuation.
struct SegmentProgress {
    beam: PolarizedBeam,
    distance_m: f64,
    collision: bool,
    volume_limit: bool,
}
impl SegmentProgress {
    fn homogeneous(beam: PolarizedBeam, distance_m: f64) -> Result<Self, SceneError> {
        Ok(Self { beam: beam.advance(distance_m)?, distance_m,
            collision: false, volume_limit: false })
    }
}

/// Prepared bindings borrow the actual native scene and spectral materials.
/// Admission sorts once; paths have no owned heap workspace and do not clone
/// geometry. Traversal retains InstanceScene's per-object BLAS/chart behavior.
#[derive(Debug)]
pub struct PolarizedScene<'a> {
    geometry: &'a InstanceScene,
    media: &'a [DielectricGlass],
    conductors: &'a [ConductorOptics],
    coatings: &'a [CoherentDielectricCoating],
    reflectances: &'a [SpectralReflectance],
    coating_stacks: Vec<Option<DielectricCoatingStack<'a>>>,
    bindings: Vec<SurfaceBinding>,
}
impl<'a> PolarizedScene<'a> {
    pub fn try_new(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        Self::try_new_with_conductors(geometry, media, &[], bindings)
    }

    /// Bind native conductor optics once, by index, without copying large
    /// tables/coatings into every surface binding or altering dielectric paths.
    pub fn try_new_with_conductors(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        conductors: &'a [ConductorOptics],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        Self::try_new_with_optics(geometry, media, conductors, &[], bindings)
    }

    /// Borrow source-bound smooth dielectric films and metal optics alongside
    /// the original media/geometry. Coating data are not copied into every hit.
    pub fn try_new_with_optics(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        conductors: &'a [ConductorOptics],
        coatings: &'a [CoherentDielectricCoating],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        Self::try_new_with_materials(geometry, media, conductors, coatings, &[], bindings)
    }

    /// Borrow original optical and diffuse tables once. Surface records retain
    /// indices, keeping spectral tables out of every copied interaction.
    pub fn try_new_with_materials(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        conductors: &'a [ConductorOptics],
        coatings: &'a [CoherentDielectricCoating],
        reflectances: &'a [SpectralReflectance],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        let bindings = admit_bindings(
            geometry, media.len(), conductors.len(), coatings.len(), reflectances.len(), bindings,
        )?;
        Ok(Self { geometry, media, conductors, coatings, reflectances,
            coating_stacks: Vec::new(), bindings })
    }

    /// Deterministic counter-based sampling keyed by `(seed, path_id, event)`.
    /// This stream is separate from scalar/camera renderers and worker order.
    /// Reusing a path ID intentionally correlates wavelengths or repeated runs.
    pub fn trace(
        &self, cx: &Cx<'_>, beam: PolarizedBeam, initial_medium: usize,
        settings: TraceSettings, seed: u64, path_id: u64,
    ) -> Result<PathResult, SceneError> {
        self.trace_with_samples(cx, beam, initial_medium, settings,
            |event| path_draw(seed, path_id, event))
    }

    /// Forward polarized transport with full Smith/GGX internal reflections at
    /// every `RoughConductor` and `BrushedConductor` hit. Each geometric hit has
    /// an independent facet budget and random stream. The native body-locked
    /// anisotropic frame, mark and endpoint shading-normal adjoint are retained.
    /// This entry point uses homogeneous straight segments.
    ///
    /// Microscopic reflections retain local complex Fresnel phase and coatings.
    /// They do not add metric distance or claim resolved groove geometry. A
    /// facet limit returns an error, never fabricated absorption. This forward
    /// power sampler supplies no outgoing marginal PDF for camera/BDPT MIS.
    #[allow(clippy::too_many_arguments)]
    pub fn trace_with_smith_multiple_scattering(
        &self, cx: &Cx<'_>, beam: PolarizedBeam, initial_medium: usize,
        settings: TraceSettings, microsurface: SmithWalkSettings, seed: u64, path_id: u64,
    ) -> Result<PathResult, SceneError> {
        self.trace_segments_with_smith(cx, beam, initial_medium, settings,
            |event| path_draw(seed, path_id, event),
            |beam, _, distance| SegmentProgress::homogeneous(beam, distance),
            None, Some((microsurface, seed, path_id)))
    }

    /// Trace with externally supplied `[microfacet_u, microfacet_v, event]`
    /// samples. The callback receives the zero-based surface-event index and
    /// is called only for dielectric, rough or brushed conductor hits, after
    /// geometry and medium admission. Rough/brushed metals use two samples;
    /// smooth metals consume none.
    /// This scene uses piecewise-homogeneous straight segments. A sampled
    /// gradient-index medium refuses at launch, before any traversal or sample;
    /// its curved boundary intersections are not implemented by this owner.
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<PathResult, SceneError> {
        self.trace_segments(cx, beam, initial_medium, settings, samples,
            |beam, _, distance| SegmentProgress::homogeneous(beam, distance), None)
    }

    // A monomorphized segment continuation shares ALL geometry, surface laws,
    // media transitions and receiver semantics with the surface-only API.
    #[allow(clippy::too_many_arguments)]
    fn trace_segments(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
        propagate: impl FnMut(PolarizedBeam, usize, f64) -> Result<SegmentProgress, SceneError>,
        mirror_interactions: Option<&mut Vec<MirrorInteraction>>,
    ) -> Result<PathResult, SceneError> {
        self.trace_segments_with_smith(cx, beam, initial_medium, settings, samples,
            propagate, mirror_interactions, None)
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_segments_with_smith(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
        mut propagate: impl FnMut(PolarizedBeam, usize, f64) -> Result<SegmentProgress, SceneError>,
        mut mirror_interactions: Option<&mut Vec<MirrorInteraction>>,
        microsurface: Option<(SmithWalkSettings, u64, u64)>,
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
        while result.surface_events < settings.max_surface_events {
            cx.checkpoint()?;
            if result.beam.power().i == 0.0 {
                result.termination = PathTermination::Extinguished;
                return finish(cx, result);
            }
            let remaining = settings.max_distance_m - result.traveled_m;
            // Find a prospective surface ONCE, but let a physical volume event
            // preempt it. Extinction starts at the physical origin, including
            // any prefix suppressed only in the geometry query.
            let hit = if remaining <= offset {
                None
            } else {
                let mut query = result.beam.ray();
                query.origin = query.at(offset);
                if ![query.origin.x, query.origin.y, query.origin.z].iter().all(|x| x.is_finite())
                    || (offset > 0.0 && query.origin == result.beam.ray().origin) {
                    return Err(SceneError::UnrepresentableProgress);
                }
                self.geometry.intersect(
                    cx, &query, remaining - offset, settings.surface_epsilon_m,
                )?
            };
            let distance = if let Some(hit) = hit {
                if !hit.hit.t.is_finite() || hit.hit.t <= 0.0 || hit.hit.t > remaining - offset {
                    return Err(SceneError::InvalidHit);
                }
                (offset + hit.hit.t).min(remaining)
            } else { remaining.max(0.0) };
            let progress = propagate(result.beam, result.medium, distance)?;
            if !progress.distance_m.is_finite() || progress.distance_m < 0.0
                || progress.distance_m > distance {
                return Err(SceneError::UnrepresentableProgress);
            }
            let traveled = result.traveled_m + progress.distance_m;
            if progress.distance_m > 0.0 && traveled <= result.traveled_m {
                return Err(SceneError::UnrepresentableProgress);
            }
            result.beam = progress.beam;
            result.traveled_m = traveled.min(settings.max_distance_m);
            if progress.collision {
                if result.beam.power().i == 0.0 || progress.volume_limit {
                    result.termination = if result.beam.power().i == 0.0 {
                        PathTermination::Extinguished
                    } else { PathTermination::VolumeLimit };
                    return finish(cx, result);
                }
                // A volume collision is not a surface hit. Re-query the changed
                // direction from its true origin without a new surface offset.
                offset = 0.0;
                continue;
            }
            let Some(hit) = hit else {
                result.traveled_m = settings.max_distance_m;
                result.termination = PathTermination::DistanceLimit;
                return finish(cx, result);
            };
            let event = result.surface_events;
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
                interaction @ (Interaction::DepolarizingDiffuse { medium, .. }
                | Interaction::RoughDiffuse { medium, .. }
                | Interaction::DepolarizingLiftedDiffuse { medium, .. }) => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let rho = depolarizing_reflectance(
                        interaction, self.reflectances, result.beam.wavelength_nm(),
                    )?;
                    if rho == 0.0 {
                        result.termination = PathTermination::Opaque;
                        return finish(cx, result);
                    }
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    let draw = samples(event);
                    result.beam = match interaction {
                        Interaction::RoughDiffuse { surface, .. } => result.beam
                            .sample_rough_diffuse(normal, surface, rho, [draw[0], draw[1]])?,
                        _ => result.beam.sample_depolarizing_diffuse(normal, rho, [draw[0], draw[1]])?,
                    }.beam;
                    offset = settings.query_offset_m;
                }
                Interaction::CoatedDiffuse { medium, reflectance, coating } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let rho = self.reflectances[reflectance].at(result.beam.wavelength_nm())
                        .map_err(|_| BeamError::Dielectric(crate::dielectric::DielectricError::InvalidWavelength))?;
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    result.beam = result.beam.sample_coated_diffuse(normal, coating, rho, samples(event))?.beam;
                    offset = settings.query_offset_m;
                }
                Interaction::SmoothConductor { optics, medium } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    result.beam = if let Some(interactions) = mirror_interactions.as_mut() {
                        let balance = result.beam.reflect_conductor_with_balance(
                            normal, &self.conductors[optics],
                        )?;
                        let reflected = balance.reflection().beam;
                        interactions.push(MirrorInteraction {
                            hit, balance, propagation_distance_m: distance,
                            surface_event_index: event,
                        });
                        reflected
                    } else {
                        result.beam.reflect_conductor(normal, &self.conductors[optics])?.beam
                    };
                    // The beam remains in the same dielectric. Absorbed metal
                    // power is removed by the response, not a second roulette.
                    offset = settings.query_offset_m;
                }
                Interaction::RoughConductor { optics, medium, surface } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    result.beam = if let Some((limits, seed, path_id)) = microsurface {
                        result.beam.sample_smith_conductor(
                            cx, normal, &self.conductors[optics], surface, limits,
                            |facet| smith_draw(seed, path_id, event, facet),
                        )?.beam()
                    } else {
                        let draw = samples(event);
                        let Some(next) = result.beam.sample_rough_conductor(
                            normal, &self.conductors[optics], surface, [draw[0], draw[1]],
                        )? else {
                            result.termination = PathTermination::RejectedMicrofacet;
                            return finish(cx, result);
                        };
                        next.beam
                    };
                    offset = settings.query_offset_m;
                }
                Interaction::BrushedConductor { optics, medium, detail } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    let instance = &self.geometry.instances()[index];
                    result.beam = if let Some((limits, seed, path_id)) = microsurface {
                        let Some(next) = result.beam.sample_instanced_smith_brushed_conductor(
                            cx, &self.conductors[optics], &detail, instance.transform(),
                            hit.local_hit.point, normal, limits,
                            |facet| smith_draw(seed, path_id, event, facet),
                        )? else {
                            result.termination = PathTermination::RejectedMicrofacet;
                            return finish(cx, result);
                        };
                        next.beam()
                    } else {
                        let draw = samples(event);
                        let Some(next) = result.beam.sample_instanced_brushed_conductor(
                            &self.conductors[optics], &detail, instance.transform(),
                            hit.local_hit.point, normal, [draw[0], draw[1]],
                        )? else {
                            result.termination = PathTermination::RejectedMicrofacet;
                            return finish(cx, result);
                        };
                        next.beam
                    };
                    offset = settings.query_offset_m;
                }
                interaction @ (Interaction::Dielectric { outside, inside, .. }
                | Interaction::CoatedDielectric { outside, inside, .. }) => {
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    let cosine = normal.dot(result.beam.ray().dir);
                    if !cosine.is_finite() || cosine == 0.0 { return Err(SceneError::InvalidHit); }
                    let (expected, target) = if cosine < 0.0 { (outside, inside) } else { (inside, outside) };
                    if result.medium != expected {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected, actual: result.medium,
                        });
                    }
                    let Some(next) = sample_boundary(
                        result.beam, normal, self.media[target].into(), interaction,
                        self.geometry.instances()[index].transform(),
                        self.coatings, self.coating_stacks.get(index).copied().flatten(), samples(event),
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

fn admit_bindings(
    geometry: &InstanceScene,
    media_count: usize,
    conductor_count: usize,
    coating_count: usize,
    reflectance_count: usize,
    mut bindings: Vec<SurfaceBinding>,
) -> Result<Vec<SurfaceBinding>, SceneError> {
    if media_count == 0
        || media_count > MAX_POLARIZED_SCENE_ITEMS
        || conductor_count > MAX_POLARIZED_SCENE_ITEMS
        || coating_count > MAX_POLARIZED_SCENE_ITEMS
        || reflectance_count > MAX_POLARIZED_SCENE_ITEMS
        || bindings.len() != geometry.instances().len()
        || bindings.len() > MAX_POLARIZED_SCENE_ITEMS
    {
        return Err(SceneError::InvalidBindings);
    }
    bindings.sort_unstable_by_key(|b| b.object_id);
    for (binding, instance) in bindings.iter().zip(geometry.instances()) {
        if binding.object_id != instance.object_id() {
            return Err(SceneError::InvalidBindings);
        }
        match binding.interaction {
            Interaction::Dielectric { outside, inside, .. }
            | Interaction::CoatedDielectric { outside, inside, .. } => {
                if outside >= media_count || inside >= media_count {
                    return Err(SceneError::InvalidBindings);
                }
                if let Interaction::CoatedDielectric { coating, .. } = binding.interaction {
                    if coating >= coating_count { return Err(SceneError::InvalidBindings); }
                }
            }
            Interaction::Detector { reference_axis: a } => {
                if ![a.x, a.y, a.z].iter().all(|x| x.is_finite())
                    || a.x.abs().max(a.y.abs()).max(a.z.abs()) == 0.0
                {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::DepolarizingDiffuse { medium, reflectance }
            | Interaction::RoughDiffuse { medium, reflectance, .. }
            | Interaction::CoatedDiffuse { medium, reflectance, .. } => {
                if medium >= media_count || reflectance >= reflectance_count {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::DepolarizingLiftedDiffuse { medium, reflectance } => {
                if medium >= media_count || reflectance.c.iter().any(|c| !c.is_finite()) {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::BrushedConductor { optics, medium, detail } => {
                let (alpha_x, alpha_y) = detail.alpha_pair();
                if optics >= conductor_count || medium >= media_count
                    || ConductorSurface::try_rough(alpha_x).is_err()
                    || ConductorSurface::try_rough(alpha_y).is_err() {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::SmoothConductor { optics, medium }
            | Interaction::RoughConductor { optics, medium, .. } => {
                if optics >= conductor_count || medium >= media_count {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::Opaque => {}
        }
    }
    Ok(bindings)
}

// One exact reflectance evaluation for forward walks and reciprocal camera
// receivers. Tabulated laws retain their original support/error semantics;
// the tracer's lifted polynomial is never extrapolated outside its native band.
pub(crate) fn depolarizing_reflectance(
    interaction: Interaction, reflectances: &[SpectralReflectance], wavelength_nm: f64,
) -> Result<f64, BeamError> {
    let invalid_wavelength = || BeamError::Dielectric(crate::dielectric::DielectricError::InvalidWavelength);
    let rho = match interaction {
        Interaction::DepolarizingDiffuse { reflectance, .. }
        | Interaction::RoughDiffuse { reflectance, .. } => reflectances
            .get(reflectance).ok_or(BeamError::InvalidDiffuseSample)?
            .at(wavelength_nm).map_err(|_| invalid_wavelength())?,
        Interaction::DepolarizingLiftedDiffuse { reflectance, .. } => {
            if !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) { return Err(invalid_wavelength()); }
            reflectance.eval(wavelength_nm)
        }
        _ => return Err(BeamError::InvalidDiffuseSample),
    };
    if !(0.0..=1.0).contains(&rho) { return Err(BeamError::NumericalRange); }
    Ok(rho)
}

// Both straight and curved scenes use exactly the same boundary sampler.
#[allow(clippy::too_many_arguments)]
fn sample_boundary(
    beam: PolarizedBeam,
    normal: Vec3,
    target: BeamMedium,
    interaction: Interaction,
    transform: RigidTransform,
    coatings: &[CoherentDielectricCoating],
    stack: Option<DielectricCoatingStack<'_>>,
    samples: [f64; 3],
) -> Result<Option<BeamSample>, BeamError> {
    match interaction {
        Interaction::Dielectric { surface, .. } => {
            let surface = bind_dielectric_surface(surface, transform)?;
            match stack {
                Some(stack) => beam.sample_stack_interface(normal, target, surface, stack, samples),
                None => beam.sample_interface(normal, target, surface, samples),
            }
        }
        Interaction::CoatedDielectric { coating, .. } => {
            if samples.iter().any(|u| !u.is_finite() || !(0.0..1.0).contains(u)) {
                return Err(crate::dielectric::DielectricError::InvalidInterface.into());
            }
            let coating = coatings.get(coating)
                .ok_or(crate::dielectric::DielectricError::InvalidInterface)?;
            beam.sample_coated_interface(normal, target, coating, samples[2]).map(Some)
        }
        _ => Err(crate::dielectric::DielectricError::InvalidInterface.into()),
    }
}

// Sampling and reciprocal camera connections share this body-to-world binding.
// Isotropic paths return the unchanged owner, without rebuilding any law.
fn bind_dielectric_surface(
    surface: DielectricSurface,
    transform: RigidTransform,
) -> Result<DielectricSurface, BeamError> {
    match surface.anisotropy() {
        Some(law) => Ok(surface.with_tangent_axis(transform.transform_vector(law.tangent_axis()))?),
        None => Ok(surface),
    }
}

fn finish(cx: &Cx<'_>, result: PathResult) -> Result<PathResult, SceneError> {
    cx.checkpoint()?;
    Ok(result)
}
/// Spectral receiving power and ideal analyzer readouts on native CIE bins.
#[cfg(feature = "photometry")]
pub mod measurement;

#[cfg(test)]
mod tests;

#[cfg(test)]
mod anisotropic_tests;

#[cfg(test)]
mod brushed_tests;

#[cfg(test)]
mod microsurface_tests;

#[cfg(test)]
mod conductor_tests {
    use super::*;
    use super::super::transport::PolarizationFrame;
    use crate::charts::TriMesh;
    use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource};
    use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
    use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
    use crate::polarization::Stokes;
    use fs_geom::Point3;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| f(&Cx::new(&gate, arena,
            StreamKey { seed: 31, kernel_id: 7, tile: 0, iteration: 0 },
            asupersync::types::Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn glass(extinction: f64) -> DielectricGlass {
        DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(extinction).unwrap(), GlassProvenance::Custom)
    }
    fn metal() -> ConductorOptics {
        let source = ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.mirror-scene", b"synthetic"),
            ConductorDataStatus::Representative,
        ).unwrap();
        ConductorOptics::try_new(core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }), source).unwrap()
    }
    fn plane(id: u64, z: f64) -> GeometryInstance {
        GeometryInstance::try_new(id,
            fs_blake3::hash_domain("org.frankensim.test.mirror-plane", &id.to_le_bytes()),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-10.0,-10.0,z], [10.0,-10.0,z], [0.0,10.0,z]], vec![[0,1,2]],
            )), RigidTransform::identity()).unwrap()
    }
    fn geometry() -> InstanceScene { InstanceScene::try_new(vec![plane(1,1.0),plane(2,-1.0)]).unwrap() }
    fn bindings(medium: usize) -> Vec<SurfaceBinding> {
        vec![SurfaceBinding { object_id:1, interaction: Interaction::SmoothConductor { optics:0, medium } },
            SurfaceBinding { object_id:2, interaction: Interaction::Detector { reference_axis:Vec3::new(1.0,0.0,0.0) } }]
    }
    fn settings() -> TraceSettings {
        TraceSettings { max_surface_events:4, max_distance_m:10.0, surface_epsilon_m:1e-9, query_offset_m:1e-8 }
    }
    fn beam(medium: DielectricGlass, wavelength: f64) -> PolarizedBeam {
        PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
            PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0)).unwrap(),
            wavelength, medium, Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap()
    }
    fn close(a: f64, b: f64) { assert!((a-b).abs() <= 2e-11*a.abs().max(b.abs()).max(1.0), "{a} != {b}"); }

    #[test]
    fn g1_mirror_exchange_uses_actual_attenuated_hit_without_changing_scene_path() {
        with_cx(|cx| {
            let geometry = geometry();
            let media = [glass(0.1)];
            let metals = [metal()];
            let scene = PolarizedScene::try_new_with_conductors(
                &geometry, &media, &metals, bindings(0),
            ).unwrap();
            let source = beam(media[0], 550.0);
            let ordinary = scene.trace(cx, source, 0, settings(), 31, 7).unwrap();
            let exchanged = scene.trace_with_mirror_balances(cx, source, 0, settings(), 31, 7).unwrap();
            assert_eq!(exchanged.path, ordinary);
            assert_eq!(exchanged.mirrors.len(), 1);
            let mirror = &exchanged.mirrors[0];
            assert_eq!(mirror.hit.object_id, 1);
            let incident = fs_math::det::exp(-0.1);
            let r = (0.8 * 0.8 + 9.0) / (1.2 * 1.2 + 9.0);
            close(mirror.balance.absorbed_power_w(), incident * (1.0 - r));
            close(mirror.balance.reflection().beam.power().i, incident * r);
            close(mirror.balance.vacuum_force_n().unwrap().z * 299_792_458.0,
                incident * (1.0 + r));
            close(exchanged.path.beam().power().i, r * fs_math::det::exp(-0.3));
        });
    }

    #[test]
    fn g1_mirror_cavity_closes_heat_and_momentum_with_surviving_limited_power() {
        with_cx(|cx| {
            let geometry = geometry();
            let media = [glass(0.0)];
            let metals = [metal()];
            let bindings = [1, 2].into_iter().map(|object_id| SurfaceBinding {
                object_id,
                interaction: Interaction::SmoothConductor { optics: 0, medium: 0 },
            }).collect();
            let scene = PolarizedScene::try_new_with_conductors(
                &geometry, &media, &metals, bindings,
            ).unwrap();
            let mut limits = settings();
            limits.max_surface_events = 3;
            let source = beam(media[0], 550.0);
            let report = scene.trace_with_mirror_balances(cx, source, 0, limits, 31, 7).unwrap();
            assert_eq!(report.path, scene.trace(cx, source, 0, limits, 31, 7).unwrap());
            assert_eq!(report.path.termination(), PathTermination::SurfaceLimit);
            assert_eq!(report.mirrors.len(), 3);
            let r = (0.8 * 0.8 + 9.0) / (1.2 * 1.2 + 9.0);
            let mut incident = 1.0;
            let mut heat = 0.0;
            let mut force_z = 0.0;
            for (i, mirror) in report.mirrors.iter().enumerate() {
                assert_eq!(mirror.hit.object_id, if i % 2 == 0 { 1 } else { 2 });
                close(mirror.balance.absorbed_power_w(), incident * (1.0 - r));
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                let force = mirror.balance.vacuum_force_n().unwrap();
                close(force.z * 299_792_458.0, sign * incident * (1.0 + r));
                heat += mirror.balance.absorbed_power_w();
                force_z += force.z;
                incident *= r;
            }
            close(report.path.beam().power().i, incident);
            close(heat + report.path.beam().power().i, 1.0);
            // The unresolved surviving beam travels opposite to the source.
            close(force_z * 299_792_458.0, 1.0 + incident);
        });
    }

    #[test]
    fn g1_native_scene_mirror_detector_retains_phase_material_and_absorption() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.1)]; let metals = [metal()];
            let scene = PolarizedScene::try_new_with_conductors(&geometry, &media, &metals, bindings(0)).unwrap();
            let source = beam(media[0], 550.0);
            let out = scene.trace_with_samples(cx, source, 0, settings(), |_| panic!("smooth metal needs no random draw")).unwrap();
            assert_eq!(out.termination(), PathTermination::Detected);
            assert_eq!(out.medium_index(),0); assert_eq!(out.beam().medium(),BeamMedium::Homogeneous(media[0]));
            assert_eq!(out.last_hit().unwrap().object_id,2);
            assert_eq!(out.surface_events(),2); close(out.traveled_m(),3.0);
            let expected = metals[0].fresnel(550.0,1.0,1.0).unwrap()*fs_math::det::exp(-0.3);
            close(out.beam().power().i,expected); close(out.beam().power().v,-expected);
            assert_eq!(scene.trace(cx, source,0,settings(),1,2).unwrap(),
                scene.trace(cx,source,0,settings(),99,999).unwrap());
        });
    }

    #[test]
    fn g1_native_scene_uses_coating_and_original_spectral_index_at_the_hit() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.0)];
            let bare = metal();
            let film = CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(),125e-9,bare.source()).unwrap();
            let metals = [bare.with_coating(film)];
            let scene = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut values = Vec::new();
            for wavelength in [400.0,550.0,750.0] {
                let out = scene.trace(cx,beam(media[0],wavelength),0,settings(),1,2).unwrap();
                let expected = metals[0].fresnel(wavelength,1.0,1.0).unwrap();
                close(out.beam().power().i,expected); close(out.beam().power().v,-expected);
                values.push(out.beam().power().i);
            }
            assert!((values[0]-values[1]).abs()>1e-3);
        });
    }

    #[test]
    fn g0_scene_metal_bindings_medium_mismatch_and_limits_fail_explicitly() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.0),glass(0.5)]; let metals = [metal()];
            assert!(matches!(PolarizedScene::try_new(&geometry,&media,bindings(0)), Err(SceneError::InvalidBindings)));
            assert!(matches!(PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(2)), Err(SceneError::InvalidBindings)));
            let wrong = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(1)).unwrap();
            let source = beam(media[0],550.0); let original = source;
            assert_eq!(wrong.trace(cx,source,0,settings(),1,2),Err(SceneError::MediumMismatch {object_id:1,expected:1,actual:0}));
            assert_eq!(source,original);
            let scene = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut limits = settings(); limits.max_surface_events=1;
            let out = scene.trace(cx,source,0,limits,1,2).unwrap();
            assert_eq!(out.termination(),PathTermination::SurfaceLimit);
            assert_eq!(out.medium_index(),0);
            close(out.beam().power().i,metals[0].fresnel(550.0,1.0,1.0).unwrap());
        });
    }

    #[cfg(feature = "photometry")]
    #[test]
    fn g1_calibrated_source_mirror_polarimeter_keeps_circular_sample_power() {
        use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum, SPECTRAL_BINS};
        use crate::polarization::photometry::SpectralOptics;
        with_cx(|cx| {
            let geometry=geometry(); let media=[glass(0.0)]; let metals=[metal()];
            let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut meter=scene.power_meter(2).unwrap();
            let source=PhotometricSource::try_new(RigidTransform::identity(),1000.0,1.0,
                AngularDistribution::isotropic(),SourceSpectrum::try_from_bins([1.0;SPECTRAL_BINS]).unwrap()).unwrap()
                .with_unpolarized_optics()
                .through(&SpectralOptics::linear_polarizer(core::f64::consts::FRAC_PI_4).unwrap()).unwrap()
                .through(&SpectralOptics::linear_retarder(0.0,[core::f64::consts::FRAC_PI_2;SPECTRAL_BINS]).unwrap()).unwrap();
            for bin in 0..SPECTRAL_BINS {
                // One directional estimator sample per bin, not a claim that one
                // direction converges a full-sphere source's detector integral.
                let beam=source.emit_beam(Vec3::new(0.0,0.0,1.0),bin,4.0*core::f64::consts::PI,media[0]).unwrap();
                let path=scene.trace(cx,beam,0,settings(),31,bin as u64).unwrap();
                meter.record(bin,&path).unwrap();
            }
            let reading=meter.reading().unwrap().require_complete().unwrap();
            let expected=500.0*metals[0].fresnel(550.0,1.0,1.0).unwrap();
            close(reading.photopic_stokes_lm().i,expected);
            close(reading.photopic_stokes_lm().v,expected);
            close(reading.circular_analyzer(true).unwrap().luminous_lm,expected);
            close(reading.circular_analyzer(false).unwrap().luminous_lm,0.0);
        });
    }

    #[test]
    fn g1_native_rough_metal_scene_matches_explicit_beam_and_receiver_path() {
        with_cx(|cx| {
            let geometry=geometry(); let media=[glass(0.1)]; let metal=metal();
            let film=CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(),125e-9,metal.source()).unwrap();
            let metals=[metal.with_coating(film)];
            let surface=ConductorSurface::try_rough(0.4).unwrap();
            let mut bindings=bindings(0);
            bindings[0].interaction=Interaction::RoughConductor { optics:0,medium:0,surface };
            let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings).unwrap();
            let input=beam(media[0],550.0); let draw=[0.17,0.63];
            let reflected=input.advance(1.0).unwrap().sample_rough_conductor(
                Vec3::new(0.0,0.0,1.0),&metals[0],surface,draw,
            ).unwrap().unwrap().beam;
            let remaining=-2.0/reflected.ray().dir.z;
            let expected=reflected.advance(remaining).unwrap().with_reference_axis(Vec3::new(1.0,0.0,0.0)).unwrap();
            let mut calls=0;
            let out=scene.trace_with_samples(cx,input,0,settings(),|event| {
                assert_eq!(event,0); calls+=1; [draw[0],draw[1],0.5]
            }).unwrap();
            assert_eq!(calls,1); assert_eq!(out.termination(),PathTermination::Detected);
            assert_eq!(out.last_hit().unwrap().object_id,2);
            assert_eq!(out.medium_index(),0); close(out.traveled_m(),1.0+remaining);
            for (a,b) in [out.beam().power().i,out.beam().power().q,out.beam().power().u,out.beam().power().v]
                .into_iter().zip([expected.power().i,expected.power().q,expected.power().u,expected.power().v]) { close(a,b); }
            assert_eq!(scene.trace(cx,input,0,settings(),31,42).unwrap(),
                scene.trace(cx,input,0,settings(),31,42).unwrap());
        });
    }

    #[test]
    fn g0_scene_rough_metal_rejection_stops_without_retry_or_invented_absorption() {
        with_cx(|cx| {
            let geometry=geometry(); let media=[glass(0.1)]; let metals=[metal()];
            let mut bindings=bindings(0);
            // The bounded GGX reflection proposal (6a13950eae) rejects no facet at alpha 1 (every
            // cap draw reflects above the horizon) nor at any alpha at the normal view (the cap is
            // exact there). At alpha 0.5 and a 45-degree view a residual rejected set remains. The
            // beam runs along (1, 0, 1)/sqrt 2 onto plane 1, whose normal facing wo is (0, 0, -1);
            // basis_all_sphere gives the local wo (-1, 0, 1)/sqrt 2, and the draw [0.98, 0] reflects
            // to a local cosine of -0.301 (the warp of tracer.rs sample_ggx_bounded_normal, replayed).
            bindings[0].interaction=Interaction::RoughConductor {
                optics:0,medium:0,surface:ConductorSurface::try_rough(0.5).unwrap(),
            };
            assert!(matches!(PolarizedScene::try_new(&geometry,&media,bindings.clone()),Err(SceneError::InvalidBindings)));
            let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings).unwrap();
            let c=core::f64::consts::FRAC_1_SQRT_2;
            let tilted=PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
                PolarizationFrame::try_new(Vec3::new(c,0.0,c),Vec3::new(c,0.0,-c)).unwrap(),
                550.0,media[0],Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap();
            let mut calls=0;
            let out=scene.trace_with_samples(cx,tilted,0,settings(),|_| {
                calls+=1; [0.98,0.0,0.5]
            }).unwrap();
            assert_eq!(calls,1); assert_eq!(out.termination(),PathTermination::RejectedMicrofacet);
            assert_eq!(out.surface_events(),1); assert_eq!(out.medium_index(),0);
            close(out.traveled_m(),1.0/c);
            close(out.beam().power().i,fs_math::det::exp(-0.1/c));
        });
    }
}

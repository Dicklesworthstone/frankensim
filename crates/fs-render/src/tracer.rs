//! SPECTRAL PATH TRACER v1 (bead 872c, WS3; [F] — behind the `tracer`
//! feature): hero-wavelength Monte-Carlo transport with next-event
//! estimation + BSDF sampling combined by the crate's MIS balance
//! heuristic, over the certified chart/BVH backends, producing CIE XYZ
//! film and byte-exact EXR through fs-img.
//!
//! DETERMINISM: every random draw comes from a counter-based stream
//! keyed by `(pixel, sample, bounce-dimension)` — Philox 4×32-10 for
//! path decisions, optionally Owen-scrambled Sobol' for the pixel/
//! wavelength dimensions ([`Sampler::OwenSobol`], decorrelated across
//! pixels by a Philox-derived scramble seed). No draw depends on
//! scheduling, so images are bitwise invariant to tile traversal and
//! worker count, and a render RESUMED from an `spp` checkpoint equals
//! the straight-through render bitwise (the pause–serialize–resume
//! doctrine applied to images). All transcendentals in the radiance
//! path go through `fs_math::det` (goldens hash these bits — no
//! platform libm), and Fresnel/roughness powers are explicit
//! multiplications, never `powi` (the a55x/4xnt hazard class).
//!
//! The legacy one-rectangle/no-environment path retains its v1 random stream
//! and non-grazing arithmetic; lighting bit-semantics v2 makes its rectangle
//! forward/reverse support boundary identical. The opt-in lighting-v1 extension admits multiple static
//! rectangular area lights and one canonical lat-long environment, ordered and
//! importance-sampled independently of caller construction order. Rectangular
//! lights are also scene geometry so BSDF paths find them (MIS-weighted both
//! ways); materials are Lambertian, legacy Schlick-GGX, exact-Fresnel rough
//! conductors, and smooth/rough spectral dielectrics; no volumetric media (the
//! `volumes` module is separate); no Russian roulette (fixed depth keeps work
//! deterministic).

use crate::animated_instances::{AnimatedGeometryInstance, AnimatedInstanceError};
use crate::aov::{
    AdaptiveAovAccumulator, AdaptiveCinematicAovFilm, AlignedAovPrimary, AlignedAovSample,
    CinematicAovConfig, CinematicAovError, CinematicAovFilm, CinematicAovPalette,
    CinematicAovTileAccumulator, MAX_EXACT_F32_INTEGER, adaptive_render_binding, render_binding,
    validate_binding, validate_reference_times,
};
use crate::camera::{
    AnimatedCamera, CameraError, CameraExposure, CutSide, KeyframeFocus, LensSample,
    OpticalCenterProjection, PhysicalCamera,
};
use crate::charts::{Hit, Ray, TraceTermination, TriMesh, sphere_trace};
use crate::conductor::{
    CONDUCTOR_BSDF_SEMANTICS_VERSION, ConductorError, ConductorOptics, ConductorSurface,
};
use crate::dielectric::{
    BeerLambertParameters, DielectricError, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance, evaluate_rough_dielectric, fresnel_dielectric, sample_rough_dielectric,
    sample_smooth_dielectric,
};
use crate::instances::{
    GeometryInstance, InstanceError, InstanceHit, InstanceSurfaceFeature, SharedGeometry,
};
use crate::lighting::{AdmittedLighting, EnvironmentMap, LightSample, LightingError};
use crate::motion::{NormalizedShutterTime, ShutterInterval, TimedRay};
use crate::motion_vectors::{
    MotionEndpoint, MotionFrame, MotionVectorComputation, MotionVectorError, PrimarySurfaceSample,
    RasterSize, compute_motion_vectors,
};
use crate::spectral::{
    LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, cie_x, cie_y, cie_z, xyz_e_to_d65, xyz_to_linear_srgb,
    y_integral,
};
use crate::{balance_heuristic, hero_wavelengths};
use core::mem::size_of;
use fs_alloc::{AllocError, LeaseCharge, LeaseReceipt, LeaseRefusal, OperationMemoryLease};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{
    Budget, Cancelled, Cx, ExecMode, LocalTaskCaps, ParkedTilePool, PoolConfig, RunError, RunId,
    RunReport, TileKernel, TilePlan, TilePool,
};
use fs_geom::{Chart, Point3, Vec3};
use fs_math::det;
use fs_rand::philox::philox4x32_10;
use fs_rand::qmc::Sobol;
use std::collections::BTreeSet;
use std::num::NonZeroU32;
use std::ops::ControlFlow;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

mod checkpoint;
pub use checkpoint::{
    RENDER_CHECKPOINT_CONTENT_DOMAIN, RENDER_CHECKPOINT_EXECUTION_ENVIRONMENT_DOMAIN,
    RENDER_CHECKPOINT_JOB_DOMAIN, RENDER_CHECKPOINT_SCHEMA_VERSION, RenderCheckpointBinding,
    RenderCheckpointError, RenderCheckpointKind, RenderCheckpointReceipt,
    RenderCheckpointWriteError, adaptive_checkpoint_job_identity, uniform_checkpoint_job_identity,
};
mod sharding;
pub use sharding::{
    RENDER_SHARD_ARTIFACT_DOMAIN, RENDER_SHARD_IDENTITY_DOMAIN, RENDER_SHARD_SCHEMA_VERSION,
    RENDER_SHARD_SEMANTICS_VERSION, RenderShardError, RenderShardLimits, RenderShardMergeLimits,
    UniformRenderShardResult, UniformRenderShardSpec, merge_uniform_shards, render_cinematic_shard,
    render_motion_shard, render_static_shard,
};

/// Bit-affecting semantic surface version of the tracer (see
/// golden-couplings.json): the path-integrator estimator shape, the
/// Philox/Sobol stream keying, the BSDF forms, the CMF/adaptation
/// constants it inherits from `spectral`, and the EXR channel layout.
/// Bump ONLY with a semantic justification per docs/GOLDEN_POLICY.md.
pub const TRACER_BIT_SEMANTICS_VERSION: u32 = 2;

/// Bit-affecting semantics of the optional motion-time path. This is versioned
/// separately because the legacy static entry points do not draw a shutter
/// dimension and retain their existing image bits.
pub const MOTION_TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

/// Bit-affecting semantics of the opt-in cinematic-camera path. This is
/// versioned independently so adding lens and camera-trajectory dimensions
/// cannot silently change the legacy pinhole tracer's frozen stream.
pub const CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

/// Bit-affecting semantics of the opt-in spectral dielectric path. Existing
/// opaque materials retain tracer-v1 stream order and image bits.
pub const DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION: u32 = 5;

/// Domain for deterministic identities of complete tracer material values.
pub const MATERIAL_CONTENT_IDENTITY_DOMAIN: &str = "org.frankensim.render.material.v1";

/// Bit-affecting semantics of construction-order-independent multi-light and
/// environment sampling. The legacy one-rectangle/no-environment path remains
/// under [`TRACER_BIT_SEMANTICS_VERSION`] and keeps its frozen stream.
pub const LIGHTING_TRACER_BIT_SEMANTICS_VERSION: u32 = 2;

/// Dedicated Philox counter domain for the two lens coordinates. Lens draws
/// never advance [`PathRng`] and therefore cannot perturb light or BSDF draws.
const CAMERA_LENS_SAMPLE_DOMAIN_V1: u32 = 0x6c65_6e73;

const PI: f64 = core::f64::consts::PI;
/// Hero-wavelength packet width (the bead's 4-wavelength packets).
pub const PACKET: usize = 4;
#[cfg(test)]
const FORCED_TRANSMISSION_EVENT_SAMPLE: f64 = 1.0 - f64::EPSILON;
/// Self-intersection offset along the normal when spawning rays.
const RAY_EPS: f64 = 1e-6;
/// Sphere-trace surface tolerance.
const TRACE_EPS: f64 = 1e-7;
const MAX_MEDIUM_STACK_DEPTH: usize = 64;
const RECT_LIGHT_GEOMETRY_REL_TOLERANCE: f64 = 1.0e-10;
const SLAB_PARALLEL_COSINE_TOLERANCE: f64 = 2.0e-10;
const SLAB_CONNECTION_REL_TOLERANCE: f64 = 2.0e-8;
const SLAB_CONNECTION_BISECTION_STEPS: usize = 80;
const PATH_ANIMATED_INSTANCE_CACHE_CAPACITY: usize = 4;

struct PathTime {
    interval: ShutterInterval,
    normalized: NormalizedShutterTime,
    cached_animated: [Option<CachedAnimatedInstance>; PATH_ANIMATED_INSTANCE_CACHE_CAPACITY],
}

struct CachedAnimatedInstance {
    primitive_index: usize,
    instance: GeometryInstance,
}

#[derive(Clone, Copy)]
struct SceneIntersection {
    primitive_index: usize,
    hit: Hit,
    instance_hit: Option<InstanceHit>,
}

/// Accepted frontmost surface for one beauty-path sample.
///
/// `primitive_index` is scene-local diagnostic metadata, not a stable object
/// identity. Instanced geometry additionally carries [`PrimarySurfaceSample`],
/// whose object/geometry/material/feature tuple is stable across rigid poses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrimaryTraceHit {
    /// Index of the accepted primitive in [`Scene::primitives`].
    pub primitive_index: usize,
    /// Exact world-space hit shaded by the beauty path.
    pub hit: Hit,
    /// Face-forwarded world-space surface frame actually consumed by the
    /// current beauty integrator. This remains separate from the backend's
    /// pre-face-forward `Hit::shading_normal` diagnostic.
    pub beauty_shading_normal_world: Vec3,
    /// Deterministic identity of the complete material value on that primitive.
    pub material_identity: ContentHash,
    /// Stable local-space identity and correspondence for instanced geometry.
    /// Direct legacy mesh/chart primitives have no object identity and return
    /// `None` rather than inventing one.
    pub surface: Option<PrimarySurfaceSample>,
}

struct PathTraceSample {
    xyz: [f64; 3],
    primary: Option<PrimaryTraceHit>,
    absolute_time_s: Option<f64>,
    pixel_jitter: [f64; 2],
    contribution_split: Option<PathContributionSplit>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct PathContributionSplit {
    direct_xyz: [f64; 3],
    indirect_xyz: [f64; 3],
    emission_xyz: [f64; 3],
}

impl PathContributionSplit {
    const ZERO: Self = Self {
        direct_xyz: [0.0; 3],
        indirect_xyz: [0.0; 3],
        emission_xyz: [0.0; 3],
    };

    fn add_assign(&mut self, other: Self) {
        add_xyz(&mut self.direct_xyz, other.direct_xyz);
        add_xyz(&mut self.indirect_xyz, other.indirect_xyz);
        add_xyz(&mut self.emission_xyz, other.emission_xyz);
    }
}

fn add_xyz(target: &mut [f64; 3], value: [f64; 3]) {
    for channel in 0..3 {
        target[channel] += value[channel];
    }
}

#[derive(Clone, Copy)]
enum PathContributionClass {
    Direct,
    Indirect,
    Emission,
}

#[derive(Clone, Copy)]
struct PathContributionRadiance {
    direct: [f64; PACKET],
    indirect: [f64; PACKET],
    emission: [f64; PACKET],
}

impl PathContributionRadiance {
    const ZERO: Self = Self {
        direct: [0.0; PACKET],
        indirect: [0.0; PACKET],
        emission: [0.0; PACKET],
    };

    fn record(&mut self, class: PathContributionClass, lane: usize, value: f64) {
        let target = match class {
            PathContributionClass::Direct => &mut self.direct,
            PathContributionClass::Indirect => &mut self.indirect,
            PathContributionClass::Emission => &mut self.emission,
        };
        target[lane] += value;
    }
}

/// One deterministic cinematic beauty sample and its exactly aligned primary
/// surface record.
///
/// This is the lossless producer seam for depth, normal, identity, and motion
/// AOVs. It reports the accepted primary hit from the same path traversal that
/// produced `xyz`; it does not retrace the pixel. Visibility at a different
/// reference time still requires [`crate::motion_vectors::validate_reprojection`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CinematicPixelSample {
    /// Unaccumulated CIE XYZ beauty contribution for this logical sample.
    pub xyz: [f64; 3],
    /// Absolute shutter time in seconds used by camera and animated geometry.
    pub absolute_time_s: f64,
    /// Linear-light direct-illumination contribution in CIE XYZ.
    pub direct_xyz: [f64; 3],
    /// Linear-light multi-bounce contribution in CIE XYZ.
    pub indirect_xyz: [f64; 3],
    /// Camera-visible emitter/environment contribution in CIE XYZ.
    pub emission_xyz: [f64; 3],
    /// Frontmost accepted surface, or `None` for a background ray.
    pub primary: Option<PrimaryTraceHit>,
}

#[derive(Clone, Copy)]
struct SurfaceFrame {
    oriented: Vec3,
    geometric: Vec3,
    entering: bool,
}

#[derive(Clone, Copy)]
struct MediumEntry {
    boundary_primitive: usize,
    glass: DielectricGlass,
}

#[derive(Clone, Copy)]
struct MediumStack {
    entries: [Option<MediumEntry>; MAX_MEDIUM_STACK_DEPTH],
    len: usize,
}

impl MediumStack {
    const fn new() -> Self {
        Self {
            entries: [None; MAX_MEDIUM_STACK_DEPTH],
            len: 0,
        }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn last(&self) -> Option<&MediumEntry> {
        self.len
            .checked_sub(1)
            .and_then(|index| self.entries[index].as_ref())
    }

    fn get(&self, index: usize) -> Option<&MediumEntry> {
        (index < self.len)
            .then(|| self.entries[index].as_ref())
            .flatten()
    }

    fn iter(&self) -> impl Iterator<Item = &MediumEntry> {
        self.entries[..self.len]
            .iter()
            .map(|entry| entry.as_ref().expect("occupied medium-stack prefix"))
    }

    fn push(&mut self, entry: MediumEntry) -> Result<(), TracerError> {
        if self.len == MAX_MEDIUM_STACK_DEPTH {
            return Err(TracerError::MediumStackOverflow);
        }
        self.entries[self.len] = Some(entry);
        self.len += 1;
        Ok(())
    }

    fn pop(&mut self) {
        if let Some(index) = self.len.checked_sub(1) {
            self.entries[index] = None;
            self.len = index;
        }
    }
}

#[derive(Clone, Copy)]
enum MediumTransition {
    Enter(MediumEntry),
    Exit { boundary_primitive: usize },
}

#[derive(Clone, Copy)]
struct BoundaryMedia {
    incident: Option<DielectricGlass>,
    transmitted: Option<DielectricGlass>,
    transition: MediumTransition,
}

#[derive(Clone, Copy)]
struct PreviousBsdf {
    pdf: f64,
    delta: bool,
    /// Only opaque source vertices are covered by the parallel-slab NEE
    /// technique. Retaining the oriented geometric normal also lets reverse
    /// MIS replay the exact straight shadow-ray spawn used by forward NEE.
    /// `None` prevents a rough-dielectric vertex from acquiring a competitor
    /// that cannot evaluate its source BSDF.
    opaque_source_geometric_normal: Option<Vec3>,
    smooth_slab: Option<SmoothSlabPath>,
}

#[derive(Clone, Copy)]
enum SmoothSlabPath {
    Entered(SmoothSlabEntry),
    Exited(SmoothSlabExit),
}

#[derive(Clone, Copy)]
struct SmoothSlabEntry {
    source_origin: Point3,
    source_geometric_normal: Vec3,
    source_direction: Vec3,
    source_pdf_solid_angle: f64,
    slab: ParallelSlab,
    entry_transmission_probability: f64,
}

#[derive(Clone, Copy)]
struct SmoothSlabExit {
    source_origin: Point3,
    source_geometric_normal: Vec3,
    source_direction: Vec3,
    source_pdf_solid_angle: f64,
    slab: ParallelSlab,
    transmission_probability: f64,
}

#[derive(Clone, Copy)]
struct DielectricPathSample {
    direction: Vec3,
    event: DielectricEvent,
    pdf: f64,
    delta: bool,
    weights: [f64; PACKET],
}

#[derive(Clone, Copy)]
struct DielectricLaneContinuation {
    direction: Vec3,
    event: DielectricEvent,
    weight: f64,
    pdf: f64,
    delta: bool,
}

#[derive(Clone, Copy)]
struct SpectralPathState {
    ray: Ray,
    throughput: [f64; PACKET],
    previous_bsdf: Option<PreviousBsdf>,
    prev_origin: Point3,
    segment_origin: Point3,
    medium_stack: MediumStack,
    rng: PathRng,
    next_depth: u32,
    /// `None` is the unsplit four-wavelength packet. At the first wavelength-
    /// dependent dielectric boundary, each continuation owns exactly one
    /// active lane and can never split again.
    active_lane: Option<usize>,
}

#[derive(Clone, Copy)]
enum DirectLightTarget {
    Rectangle {
        primitive_index: usize,
        distance_m: f64,
        point: Point3,
        normal: Vec3,
    },
    Environment,
}

#[derive(Clone, Copy)]
struct PreparedDirectLight {
    direction: Vec3,
    emission: (LiftedSpectrum, f64),
    pdf_solid_angle: f64,
    target: DirectLightTarget,
}

#[derive(Clone, Copy)]
struct ParallelSlab {
    boundary_primitive: usize,
    glass: DielectricGlass,
    entry_reference: Point3,
    exit_reference: Point3,
    /// Unit normal from the incident ambient half-space through the slab.
    axis: Vec3,
    thickness_m: f64,
}

#[derive(Clone, Copy)]
struct SlabConnectionGeometry {
    incident_direction: Vec3,
    internal_direction: Vec3,
    /// `dA_light / dOmega_source` for a finite rectangle.  Infinite
    /// environment directions use the identity angular map instead.
    light_area_jacobian: Option<f64>,
}

#[derive(Clone, Copy)]
struct SlabDirectLane {
    incident_direction: Vec3,
    nee_pdf_solid_angle: f64,
    transmission_probability: f64,
    radiance_transport: f64,
    visible: bool,
}

#[derive(Clone, Copy)]
enum CameraPath<'a> {
    Legacy,
    Cinematic {
        camera: &'a AnimatedCamera,
        exposure: CameraExposure<'a>,
    },
}

/// The per-draw uniform stream: Philox keyed by (pixel, sample,
/// dimension). Counter-based — random access, no state shared between
/// pixels/samples/workers.
#[derive(Clone, Copy)]
struct PathRng {
    pixel: u32,
    sample: u32,
    dim: u32,
    key: [u32; 2],
}

impl PathRng {
    fn next2(&mut self) -> (f64, f64) {
        let out = philox4x32_10([self.pixel, self.sample, self.dim, 0x7261_7972], self.key);
        self.dim += 1;
        (u32_unit(out[0]), u32_unit(out[1]))
    }
}

fn u32_unit(x: u32) -> f64 {
    f64::from(x) / 4_294_967_296.0
}

fn prepare_direct_light(origin: Point3, sample: LightSample) -> Option<PreparedDirectLight> {
    let prepared = match sample {
        LightSample::Rectangle(sample) => {
            let displacement = sample.point.delta_from(origin);
            let distance_squared = displacement.dot(displacement);
            if !(distance_squared > 0.0 && distance_squared.is_finite()) {
                return None;
            }
            let distance_m = distance_squared.sqrt();
            PreparedDirectLight {
                direction: displacement.scale(1.0 / distance_m),
                emission: sample.emission,
                pdf_solid_angle: sample.pdf_solid_angle,
                target: DirectLightTarget::Rectangle {
                    primitive_index: sample.primitive_index,
                    distance_m,
                    point: sample.point,
                    normal: sample.normal,
                },
            }
        }
        LightSample::Environment(sample) => PreparedDirectLight {
            direction: sample.direction,
            emission: sample.emission,
            pdf_solid_angle: sample.pdf_solid_angle,
            target: DirectLightTarget::Environment,
        },
    };
    (prepared.pdf_solid_angle > 0.0
        && prepared.pdf_solid_angle.is_finite()
        && prepared.direction.x.is_finite()
        && prepared.direction.y.is_finite()
        && prepared.direction.z.is_finite())
    .then_some(prepared)
}

/// Pixel-space sampler for the (jitter-x, jitter-y, hero-λ) dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sampler {
    /// Independent Philox draws for every dimension.
    Iid,
    /// Owen-scrambled Sobol' over the three pixel dimensions,
    /// decorrelated across pixels by a Philox-derived scramble seed
    /// (the ambition-round upgrade; its equal-spp variance claim is
    /// measured in the battery, not assumed).
    OwenSobol,
}

/// How direct lighting is estimated — [`DirectStrategy::Mis`] is the
/// product setting; the single-technique modes exist so the battery
/// can MEASURE that MIS beats either alone (the bead's acceptance).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectStrategy {
    /// Next-event estimation only.
    NeeOnly,
    /// BSDF sampling only.
    BsdfOnly,
    /// Both, combined with the balance heuristic.
    Mis,
}

/// A surface material.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Material {
    /// Ideal diffuse with a spectral reflectance.
    Lambertian {
        /// Reflectance spectrum (bounded (0,1) by construction).
        reflectance: LiftedSpectrum,
    },
    /// GGX microfacet (Smith separable shadowing, Schlick Fresnel with
    /// the spectral reflectance as F0).
    Ggx {
        /// F0 reflectance spectrum.
        reflectance: LiftedSpectrum,
        /// Roughness α (GGX convention, > 0).
        alpha: f64,
    },
    /// Opaque spectral conductor using exact complex-IOR Fresnel and an
    /// isotropic single-scattering GGX microfacet surface.
    Conductor {
        /// Validated absolute complex-index table and honest provenance.
        optics: ConductorOptics,
        /// Validated isotropic GGX roughness.
        surface: ConductorSurface,
    },
    /// Homogeneous spectral dielectric boundary. Geometry must be a closed,
    /// consistently outward-oriented solid; encountered violations refuse
    /// through the path-local medium stack.
    Dielectric {
        /// Validated interior glass definition.
        glass: DielectricGlass,
        /// Smooth-delta or isotropic-GGX boundary treatment.
        surface: DielectricSurface,
    },
}

impl Material {
    /// Deterministic content identity of this complete tracer material value.
    ///
    /// This binds the numerical visual model, including declared glass
    /// provenance, but is not evidence that those parameters match a physical
    /// specimen. Emission remains separate primitive metadata.
    #[must_use]
    pub fn content_identity(self) -> ContentHash {
        let mut hasher = DomainHasher::new(MATERIAL_CONTENT_IDENTITY_DOMAIN);
        match self {
            Self::Lambertian { reflectance } => {
                hasher.update(&[0]);
                update_material_scalars(&mut hasher, reflectance.c);
            }
            Self::Ggx { reflectance, alpha } => {
                hasher.update(&[1]);
                update_material_scalars(&mut hasher, reflectance.c);
                update_material_scalar(&mut hasher, alpha);
            }
            Self::Dielectric { glass, surface } => {
                hasher.update(&[2]);
                update_material_scalars(&mut hasher, glass.ior().coefficients());
                match glass.absorption().parameters() {
                    BeerLambertParameters::Clear => hasher.update(&[0]),
                    BeerLambertParameters::Constant { extinction_per_m } => {
                        hasher.update(&[1]);
                        update_material_scalar(&mut hasher, extinction_per_m);
                    }
                    BeerLambertParameters::ReferenceRgb {
                        linear_rgb,
                        distance_m,
                    } => {
                        hasher.update(&[2]);
                        update_material_scalars(&mut hasher, linear_rgb);
                        update_material_scalar(&mut hasher, distance_m);
                    }
                }
                hasher.update(&[match glass.provenance() {
                    GlassProvenance::Custom => 0,
                    GlassProvenance::RepresentativeBorosilicateV1 => 1,
                    GlassProvenance::RepresentativeCrownV1 => 2,
                }]);
                match surface.roughness_alpha() {
                    None => hasher.update(&[0]),
                    Some(alpha) => {
                        hasher.update(&[1]);
                        update_material_scalar(&mut hasher, alpha);
                    }
                }
            }
            Self::Conductor { optics, surface } => {
                // Tags 0--2 are frozen legacy identities. Conductor is an
                // additive material family with independently versioned bits.
                hasher.update(&[3]);
                hasher.update(&CONDUCTOR_BSDF_SEMANTICS_VERSION.to_le_bytes());
                hasher.update(optics.content_identity().as_bytes());
                update_material_scalar(&mut hasher, surface.roughness_alpha());
            }
        }
        hasher.finalize()
    }
}

fn update_material_scalars(hasher: &mut DomainHasher, values: [f64; 3]) {
    for value in values {
        update_material_scalar(hasher, value);
    }
}

fn update_material_scalar(hasher: &mut DomainHasher, value: f64) {
    hasher.update(&value.to_bits().to_le_bytes());
}

/// Scene geometry: a triangle mesh (BVH) or any certified chart
/// (sphere-traced SDF/F-rep through the default [S] backend surface hardened by
/// bead 8ll9).
pub enum Shape {
    /// Triangle mesh over the deterministic median-split BVH.
    Mesh(TriMesh),
    /// A certified-Lipschitz chart, sphere-traced.
    Chart(Box<dyn Chart>),
    /// Shared immutable chart/mesh placed by a validated proper-rigid transform.
    Instance(GeometryInstance),
    /// Shared immutable chart/mesh evaluated from a rigid trajectory at the
    /// current path's absolute shutter time.
    AnimatedInstance(AnimatedGeometryInstance),
}

/// One scene object.
pub struct Primitive {
    /// Geometry.
    pub shape: Shape,
    /// Material (ignored for pure emitters in v1: lights do not
    /// reflect).
    pub material: Material,
    /// Emitted radiance: spectrum × scale (None = non-emissive).
    pub emission: Option<(LiftedSpectrum, f64)>,
}

/// Rectangular area-light metadata remains re-exported here for the existing
/// tracer API; validation and sampling live in [`crate::lighting`].
pub use crate::lighting::RectLight;

/// Pinhole camera. `half_tan` is tan(fov/2) supplied directly — the
/// library takes no trig on its API surface.
pub struct Camera {
    /// Eye point.
    pub eye: Point3,
    /// Unit view direction.
    pub forward: Vec3,
    /// Unit up (orthogonal to forward).
    pub up: Vec3,
    /// tan(vertical fov / 2).
    pub half_tan: f64,
}

/// A renderable scene.
pub struct Scene {
    /// Objects (lights included as emissive primitives).
    pub primitives: Vec<Primitive>,
    /// Static rectangular emitters eligible for next-event estimation.
    pub lights: Vec<RectLight>,
    /// Optional canonical lat-long environment emitter.
    pub environment: Option<EnvironmentMap>,
    /// Camera.
    pub camera: Camera,
}

/// Render settings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Settings {
    /// Image width (pixels).
    pub width: u32,
    /// Image height (pixels).
    pub height: u32,
    /// Samples per pixel.
    pub spp: u32,
    /// Maximum path depth (bounces).
    pub max_depth: u32,
    /// Pixel-dimension sampler.
    pub sampler: Sampler,
    /// Direct-lighting strategy.
    pub strategy: DirectStrategy,
    /// Stream seed (the replay identity).
    pub seed: u64,
}

/// Bit-affecting semantics of the deterministic adaptive stopping estimator.
///
/// Version 1 uses per-channel Welford means and second moments, checks only the
/// declared fixed sample checkpoints, and compares an IID standard-error
/// estimate against `absolute + relative * max(abs(mean), dark_floor)`. For
/// Owen-Sobol this same quantity is only a within-stream dispersion proxy; a
/// confidence interval would require independent scrambles.
pub const ADAPTIVE_SAMPLING_SEMANTICS_VERSION: u32 = 1;

/// Invalid adaptive policy or non-finite raw estimator state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveSamplingError {
    /// At least two samples are required before a sample variance exists.
    InvalidMinimumSamples,
    /// A decision batch must contain at least one sample.
    InvalidBatchSamples,
    /// [`Settings::spp`] was below the declared adaptive minimum.
    MaximumBelowMinimum,
    /// A tolerance was negative or non-finite.
    InvalidThreshold {
        /// Rejected policy field.
        field: &'static str,
    },
    /// A traced XYZ sample was non-finite.
    NonFiniteSample,
    /// A running sum, mean, or second moment became non-finite or invalid.
    InvalidMoment,
    /// The per-pixel sample counter overflowed.
    SampleCountOverflow,
}

impl core::fmt::Display for AdaptiveSamplingError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidMinimumSamples => {
                formatter.write_str("adaptive minimum samples must be at least two")
            }
            Self::InvalidBatchSamples => {
                formatter.write_str("adaptive decision batch must be nonzero")
            }
            Self::MaximumBelowMinimum => {
                formatter.write_str("render sample ceiling is below adaptive minimum")
            }
            Self::InvalidThreshold { field } => {
                write!(
                    formatter,
                    "adaptive threshold {field} must be finite and nonnegative"
                )
            }
            Self::NonFiniteSample => formatter.write_str("adaptive sample was non-finite"),
            Self::InvalidMoment => {
                formatter.write_str("adaptive running moment became non-finite or negative")
            }
            Self::SampleCountOverflow => {
                formatter.write_str("adaptive per-pixel sample count overflowed")
            }
        }
    }
}

impl core::error::Error for AdaptiveSamplingError {}

/// Deterministic raw-estimator stopping policy.
///
/// Decisions occur first at `minimum_samples`, then every `batch_samples`, and
/// finally at [`Settings::spp`] even when the ceiling is not batch-aligned.
/// Denoised pixels are deliberately absent from this API: only raw path
/// samples can satisfy the estimator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveSamplingConfig {
    minimum_samples: u32,
    batch_samples: u32,
    absolute_error: f64,
    relative_error: f64,
    dark_floor: f64,
}

impl AdaptiveSamplingConfig {
    /// Validate a deterministic adaptive policy.
    pub fn try_new(
        minimum_samples: u32,
        batch_samples: u32,
        absolute_error: f64,
        relative_error: f64,
        dark_floor: f64,
    ) -> Result<Self, AdaptiveSamplingError> {
        if minimum_samples < 2 {
            return Err(AdaptiveSamplingError::InvalidMinimumSamples);
        }
        if batch_samples == 0 {
            return Err(AdaptiveSamplingError::InvalidBatchSamples);
        }
        for (field, value) in [
            ("absolute_error", absolute_error),
            ("relative_error", relative_error),
            ("dark_floor", dark_floor),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(AdaptiveSamplingError::InvalidThreshold { field });
            }
        }
        let canonical_zero = |value: f64| if value == 0.0 { 0.0 } else { value };
        Ok(Self {
            minimum_samples,
            batch_samples,
            absolute_error: canonical_zero(absolute_error),
            relative_error: canonical_zero(relative_error),
            dark_floor: canonical_zero(dark_floor),
        })
    }

    /// First sample count at which convergence may be declared.
    #[must_use]
    pub const fn minimum_samples(self) -> u32 {
        self.minimum_samples
    }

    /// Spacing between deterministic decision checkpoints.
    #[must_use]
    pub const fn batch_samples(self) -> u32 {
        self.batch_samples
    }

    /// Per-channel absolute dispersion allowance in raw XYZ units.
    #[must_use]
    pub const fn absolute_error(self) -> f64 {
        self.absolute_error
    }

    /// Per-channel relative dispersion allowance.
    #[must_use]
    pub const fn relative_error(self) -> f64 {
        self.relative_error
    }

    /// Lower scale used by the relative term for dark channels.
    #[must_use]
    pub const fn dark_floor(self) -> f64 {
        self.dark_floor
    }

    fn validate_maximum(self, maximum_samples: u32) -> Result<(), AdaptiveSamplingError> {
        if maximum_samples < self.minimum_samples {
            Err(AdaptiveSamplingError::MaximumBelowMinimum)
        } else {
            Ok(())
        }
    }

    fn is_checkpoint(self, samples: u32, maximum_samples: u32) -> bool {
        samples == maximum_samples
            || (samples >= self.minimum_samples
                && (samples - self.minimum_samples).is_multiple_of(self.batch_samples))
    }

    fn accepts(self, dispersion: f64, mean: f64) -> bool {
        if !dispersion.is_finite() {
            return false;
        }
        if dispersion <= self.absolute_error {
            return true;
        }
        let scale = mean.abs().max(self.dark_floor);
        scale > 0.0 && (dispersion - self.absolute_error) / scale <= self.relative_error
    }
}

/// Why one adaptive pixel stopped consuming paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum AdaptiveDecision {
    /// All three raw XYZ standard-error proxies met the declared threshold at
    /// a deterministic decision checkpoint.
    ErrorThreshold,
    /// The hard [`Settings::spp`] ceiling was reached without satisfying the
    /// threshold. When both happen at the final checkpoint,
    /// [`Self::ErrorThreshold`] wins so the decision records that the target
    /// was achieved.
    MaximumSamples,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct AdaptivePixelAccumulator {
    sum_xyz: [f64; 3],
    mean_xyz: [f64; 3],
    m2_xyz: [f64; 3],
    samples: u32,
    decision: Option<AdaptiveDecision>,
}

impl AdaptivePixelAccumulator {
    const EMPTY: Self = Self {
        sum_xyz: [0.0; 3],
        mean_xyz: [0.0; 3],
        m2_xyz: [0.0; 3],
        samples: 0,
        decision: None,
    };

    fn push(&mut self, xyz: [f64; 3]) -> Result<(), AdaptiveSamplingError> {
        if xyz.iter().any(|value| !value.is_finite()) {
            return Err(AdaptiveSamplingError::NonFiniteSample);
        }
        let next_samples = self
            .samples
            .checked_add(1)
            .ok_or(AdaptiveSamplingError::SampleCountOverflow)?;
        let mut next_sum_xyz = self.sum_xyz;
        let mut next_mean_xyz = self.mean_xyz;
        let mut next_m2_xyz = self.m2_xyz;
        for channel in 0..3 {
            let next_sum = self.sum_xyz[channel] + xyz[channel];
            let delta = xyz[channel] - self.mean_xyz[channel];
            let next_mean = self.mean_xyz[channel] + delta / f64::from(next_samples);
            let next_m2 = delta.mul_add(xyz[channel] - next_mean, self.m2_xyz[channel]);
            if !next_sum.is_finite()
                || !next_mean.is_finite()
                || !next_m2.is_finite()
                || next_m2 < 0.0
            {
                return Err(AdaptiveSamplingError::InvalidMoment);
            }
            next_sum_xyz[channel] = next_sum;
            next_mean_xyz[channel] = next_mean;
            next_m2_xyz[channel] = if next_m2 == 0.0 { 0.0 } else { next_m2 };
        }
        self.sum_xyz = next_sum_xyz;
        self.mean_xyz = next_mean_xyz;
        self.m2_xyz = next_m2_xyz;
        self.samples = next_samples;
        Ok(())
    }

    fn mean_xyz(self) -> [f64; 3] {
        self.mean_xyz
    }

    fn sample_variance_xyz(self) -> [f64; 3] {
        if self.samples < 2 {
            return [f64::INFINITY; 3];
        }
        let inverse = 1.0 / f64::from(self.samples - 1);
        self.m2_xyz.map(|m2| m2 * inverse)
    }

    fn dispersion_proxy_xyz(self) -> [f64; 3] {
        let inverse_samples = 1.0 / f64::from(self.samples.max(1));
        self.sample_variance_xyz()
            .map(|variance| det::sqrt(variance * inverse_samples))
    }

    fn meets(self, policy: AdaptiveSamplingConfig) -> bool {
        let mean = self.mean_xyz();
        self.dispersion_proxy_xyz()
            .into_iter()
            .zip(mean)
            .all(|(dispersion, mean)| policy.accepts(dispersion, mean))
    }

    fn decision(
        self,
        policy: AdaptiveSamplingConfig,
        maximum_samples: u32,
    ) -> Option<AdaptiveDecision> {
        if !policy.is_checkpoint(self.samples, maximum_samples) {
            return None;
        }
        if self.meets(policy) {
            Some(AdaptiveDecision::ErrorThreshold)
        } else if self.samples == maximum_samples {
            Some(AdaptiveDecision::MaximumSamples)
        } else {
            None
        }
    }
}

/// Largest admitted edge of one logical image tile. This is an allocation
/// guard, not a performance recommendation; cinematic quality profiles use
/// 32 x 32 tiles today.
pub const MAX_RENDER_TILE_EDGE: u32 = 4_096;
/// Defensive ceiling on renderer-created worker threads.
pub const MAX_RENDER_WORKERS: usize = 256;
const RENDER_TILE_KERNEL: &str = "fs-render/spectral-film-tile-v1";
const CINEMATIC_AOV_TILE_KERNEL: &str = "fs-render/cinematic-aov-tile-v1";
/// Independent tile program identity for the adaptive cinematic AOV path.
///
/// The adaptive stop decision remains local to each pixel and uses only the
/// raw XYZ Welford state.  The co-staged AOV tile observes that same accepted
/// prefix, but is deliberately not part of this kernel's decision identity.
const ADAPTIVE_CINEMATIC_AOV_TILE_KERNEL: &str = "fs-render/adaptive-cinematic-aov-tile-v1";
const PENDING_ADAPTIVE_RENDER_TILE_KERNEL: &str =
    "fs-render/pending-adaptive-spectral-film-tile-v1";

/// Invalid explicit tile-render execution policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderExecutionConfigError {
    /// A tile edge was zero or exceeded [`MAX_RENDER_TILE_EDGE`].
    InvalidTileShape,
    /// Worker count was zero or exceeded [`MAX_RENDER_WORKERS`].
    InvalidWorkerCount,
    /// The operation memory ceiling was zero.
    InvalidMemoryLimit,
    /// Explicit worker weights did not provide one positive value per worker.
    InvalidWorkerWeights,
    /// A job requested a worker count, scheduling weights, or execution mode
    /// different from the already-parked render crew.
    ParkedCrewMismatch,
    /// A resumable job was retried under a different execution mode than the
    /// one bound when its private state was created.
    ResumeModeMismatch {
        /// Mode bound into the pending job.
        expected: ExecMode,
        /// Mode supplied by the retry context.
        actual: ExecMode,
    },
    /// A resumable job was retried under a different admitted compute budget.
    ResumeBudgetMismatch,
    /// Image dimensions were zero or exceeded the frozen `u32` pixel-identity
    /// domain used by tracer-v1 random streams.
    InvalidImageDimensions,
}

impl core::fmt::Display for RenderExecutionConfigError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidTileShape => formatter.write_str("invalid render tile shape"),
            Self::InvalidWorkerCount => formatter.write_str("invalid render worker count"),
            Self::InvalidMemoryLimit => {
                formatter.write_str("invalid render operation memory limit")
            }
            Self::InvalidWorkerWeights => formatter.write_str("invalid render worker weights"),
            Self::ParkedCrewMismatch => {
                formatter.write_str("render job does not match the parked worker crew")
            }
            Self::ResumeModeMismatch { expected, actual } => write!(
                formatter,
                "pending render mode mismatch: expected {}, got {}",
                expected.name(),
                actual.name()
            ),
            Self::ResumeBudgetMismatch => {
                formatter.write_str("pending render budget differs from its admitted budget")
            }
            Self::InvalidImageDimensions => {
                formatter.write_str("image exceeds the tracer-v1 pixel identity domain")
            }
        }
    }
}

impl core::error::Error for RenderExecutionConfigError {}

/// Explicit, replayable execution policy for the tile-parallel tracer.
///
/// Tile geometry is independent of worker count. Worker weights affect only
/// scheduling, never the pixel/sample/dimension random stream or arithmetic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RenderExecutionConfig {
    tile_width: u32,
    tile_height: u32,
    workers: usize,
    memory_limit_bytes: u64,
    run_id: RunId,
    quantum_weights: Vec<u32>,
}

impl RenderExecutionConfig {
    /// Validate a bounded deterministic execution policy. An empty weight
    /// vector selects equal weights; use [`Self::with_quantum_weights`] for an
    /// explicit heterogeneous schedule.
    pub fn try_new(
        tile_width: u32,
        tile_height: u32,
        workers: usize,
        memory_limit_bytes: u64,
        run_id: RunId,
    ) -> Result<Self, RenderExecutionConfigError> {
        if tile_width == 0
            || tile_height == 0
            || tile_width > MAX_RENDER_TILE_EDGE
            || tile_height > MAX_RENDER_TILE_EDGE
        {
            return Err(RenderExecutionConfigError::InvalidTileShape);
        }
        if workers == 0 || workers > MAX_RENDER_WORKERS {
            return Err(RenderExecutionConfigError::InvalidWorkerCount);
        }
        if memory_limit_bytes == 0 {
            return Err(RenderExecutionConfigError::InvalidMemoryLimit);
        }
        Ok(Self {
            tile_width,
            tile_height,
            workers,
            memory_limit_bytes,
            run_id,
            quantum_weights: Vec::new(),
        })
    }

    /// Select one positive initial-share weight per worker. This is useful for
    /// heterogeneous cores and for proving schedule independence.
    pub fn with_quantum_weights(
        mut self,
        quantum_weights: Vec<u32>,
    ) -> Result<Self, RenderExecutionConfigError> {
        if quantum_weights.len() != self.workers || quantum_weights.contains(&0) {
            return Err(RenderExecutionConfigError::InvalidWorkerWeights);
        }
        self.quantum_weights = quantum_weights;
        Ok(self)
    }

    /// Logical tile width in pixels.
    #[must_use]
    pub const fn tile_width(&self) -> u32 {
        self.tile_width
    }

    /// Logical tile height in pixels.
    #[must_use]
    pub const fn tile_height(&self) -> u32 {
        self.tile_height
    }

    /// Worker count used to construct the throughput lane.
    #[must_use]
    pub const fn workers(&self) -> usize {
        self.workers
    }

    /// Hard operation-memory ceiling, including retained input film (for a
    /// progressive append), staging film, tile scratch, and executor metadata.
    #[must_use]
    pub const fn memory_limit_bytes(&self) -> u64 {
        self.memory_limit_bytes
    }

    /// Caller-ledgered logical run identity.
    #[must_use]
    pub const fn run_id(&self) -> RunId {
        self.run_id
    }

    /// Explicit scheduling weights, or an empty equal-weight request.
    #[must_use]
    pub fn quantum_weights(&self) -> &[u32] {
        &self.quantum_weights
    }
}

fn equivalent_quantum_weights(left: &[u32], right: &[u32], workers: usize) -> bool {
    left == right
        || (left.is_empty() && right.len() == workers && right.iter().all(|weight| *weight == 1))
        || (right.is_empty() && left.len() == workers && left.iter().all(|weight| *weight == 1))
}

/// Reusable worker topology for a sequence of render jobs.
///
/// Call [`Self::with_parked_crew_local`] once around an animation or batch.
/// The callback receives a [`ParkedRenderScope`] whose jobs wake the same
/// worker crew instead of spawning and joining threads for every frame.
pub struct RenderWorkerPool {
    pool: TilePool,
    mode: ExecMode,
    workers: usize,
    quantum_weights: Vec<u32>,
}

impl RenderWorkerPool {
    /// Probe host placement once and prepare a reusable pool. The scheduler
    /// seed affects placement only; pixel/sample random streams remain keyed
    /// exclusively by [`Settings::seed`].
    #[must_use]
    pub fn new(execution: &RenderExecutionConfig, mode: ExecMode, scheduler_seed: u64) -> Self {
        let pool = build_render_pool(execution, mode, scheduler_seed);
        Self {
            pool,
            mode,
            workers: execution.workers,
            quantum_weights: execution.quantum_weights.clone(),
        }
    }

    /// Park the configured workers for the complete callback. All render jobs
    /// issued through the supplied scope must use the same worker count,
    /// scheduling weights, and execution mode, but may choose their own tile
    /// shape, memory ceiling, and logical run identity.
    pub fn with_parked_crew_local<R>(
        &self,
        operation: impl FnOnce(&ParkedRenderScope<'_>) -> R,
    ) -> R {
        self.pool.with_parked_crew_local(|parked| {
            let scope = ParkedRenderScope {
                pool: parked,
                mode: self.mode,
                workers: self.workers,
                quantum_weights: &self.quantum_weights,
            };
            operation(&scope)
        })
    }
}

/// Callback-scoped renderer backed by an already-parked worker crew.
///
/// This value cannot escape [`RenderWorkerPool::with_parked_crew_local`], so
/// worker lifetime remains structurally joined even on unwind.
pub struct ParkedRenderScope<'a> {
    pool: &'a ParkedTilePool<'a, LocalTaskCaps>,
    mode: ExecMode,
    workers: usize,
    quantum_weights: &'a [u32],
}

impl ParkedRenderScope<'_> {
    fn validate_job(
        &self,
        cx: &Cx<'_>,
        execution: &RenderExecutionConfig,
    ) -> Result<(), RenderExecutionError> {
        if cx.mode() != self.mode
            || execution.workers != self.workers
            || !equivalent_quantum_weights(
                &execution.quantum_weights,
                self.quantum_weights,
                self.workers,
            )
        {
            return Err(RenderExecutionError::Config(
                RenderExecutionConfigError::ParkedCrewMismatch,
            ));
        }
        Ok(())
    }

    /// Render a progressive static range on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_range(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        film: &mut Film,
        from: u32,
        to: u32,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionReport, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_range_with_execution_impl(
            scene,
            cx,
            settings,
            film,
            from,
            to,
            None,
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a progressive motion-blurred range on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_motion_range(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        film: &mut Film,
        from: u32,
        to: u32,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionReport, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_range_with_execution_impl(
            scene,
            cx,
            settings,
            film,
            from,
            to,
            Some(shutter),
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a progressive cinematic-camera range on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_cinematic_range(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: &Settings,
        film: &mut Film,
        from: u32,
        to: u32,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionReport, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        let exposure = camera
            .admit_shutter(cx, shutter, cut_side)
            .map_err(TracerError::from)?;
        render_range_with_execution_impl(
            scene,
            cx,
            settings,
            film,
            from,
            to,
            Some(shutter),
            CameraPath::Cinematic { camera, exposure },
            execution,
            self.pool,
        )
    }

    /// Render a fresh static film on the parked crew.
    pub fn render(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_fresh_with_execution_impl(
            scene,
            cx,
            settings,
            None,
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a fresh motion-blurred film on the parked crew.
    pub fn render_motion(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_fresh_with_execution_impl(
            scene,
            cx,
            settings,
            Some(shutter),
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a fresh cinematic-camera film on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_cinematic(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: &Settings,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<RenderExecutionOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        let exposure = camera
            .admit_shutter(cx, shutter, cut_side)
            .map_err(TracerError::from)?;
        render_fresh_with_execution_impl(
            scene,
            cx,
            settings,
            Some(shutter),
            CameraPath::Cinematic { camera, exposure },
            execution,
            self.pool,
        )
    }

    /// Render a fresh cinematic-camera film with exactly aligned AOVs on this
    /// already parked worker crew. The film remains private until all tiles,
    /// the executor drain, and the final cancellation checkpoint succeed.
    #[allow(clippy::too_many_arguments)]
    pub fn render_cinematic_with_aovs(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: &Settings,
        shutter: ShutterInterval,
        config: CinematicAovConfig,
        execution: &RenderExecutionConfig,
    ) -> Result<CinematicAovExecutionOutput, CinematicAovExecutionError> {
        self.validate_job(cx, execution)
            .map_err(CinematicAovExecutionError::Execution)?;
        render_cinematic_with_aovs_execution_impl(
            scene, camera, cut_side, cx, settings, shutter, config, execution, self.pool,
        )
    }

    /// Render a fresh adaptive cinematic-camera film with AOVs aligned to each
    /// pixel's exact accepted sample prefix on this already parked crew.
    ///
    /// As with the serial oracle, denoising and diagnostic AOVs observe the
    /// Welford-selected prefix but cannot affect its stopping decision.  No
    /// state becomes public unless every tile drains successfully.
    #[allow(clippy::too_many_arguments)]
    pub fn render_cinematic_adaptive_with_aovs(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: &Settings,
        policy: AdaptiveSamplingConfig,
        shutter: ShutterInterval,
        config: CinematicAovConfig,
        execution: &RenderExecutionConfig,
    ) -> Result<AdaptiveCinematicAovExecutionOutput, CinematicAovExecutionError> {
        self.validate_job(cx, execution)
            .map_err(CinematicAovExecutionError::Execution)?;
        render_cinematic_adaptive_with_aovs_execution_impl(
            scene, camera, cut_side, cx, settings, policy, shutter, config, execution, self.pool,
        )
    }

    /// Render a fresh static adaptive film on the parked crew.
    pub fn render_adaptive(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        policy: AdaptiveSamplingConfig,
        execution: &RenderExecutionConfig,
    ) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_adaptive_parallel_impl(
            scene,
            cx,
            settings,
            policy,
            None,
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a fresh legacy-camera motion adaptive film on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_motion_adaptive(
        &self,
        scene: &Scene,
        cx: &Cx<'_>,
        settings: &Settings,
        policy: AdaptiveSamplingConfig,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        render_adaptive_parallel_impl(
            scene,
            cx,
            settings,
            policy,
            Some(shutter),
            CameraPath::Legacy,
            execution,
            self.pool,
        )
    }

    /// Render a fresh cinematic-camera adaptive film on the parked crew.
    #[allow(clippy::too_many_arguments)]
    pub fn render_cinematic_adaptive(
        &self,
        scene: &Scene,
        camera: &AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: &Settings,
        policy: AdaptiveSamplingConfig,
        shutter: ShutterInterval,
        execution: &RenderExecutionConfig,
    ) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
        self.validate_job(cx, execution)?;
        let exposure = camera
            .admit_shutter(cx, shutter, cut_side)
            .map_err(TracerError::from)?;
        render_adaptive_parallel_impl(
            scene,
            cx,
            settings,
            policy,
            Some(shutter),
            CameraPath::Cinematic { camera, exposure },
            execution,
            self.pool,
        )
    }
}

/// Exact bounds of one row-major logical image tile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RenderTileBounds {
    /// Inclusive x origin.
    pub x: u32,
    /// Inclusive y origin.
    pub y: u32,
    /// Width, clipped at the image edge.
    pub width: u32,
    /// Height, clipped at the image edge.
    pub height: u32,
}

/// Validated fixed logical tile layout, independent of worker count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RenderTileLayout {
    image_width: u32,
    image_height: u32,
    tile_width: u32,
    tile_height: u32,
    tiles_x: u32,
    tiles_y: u32,
    tile_count: u64,
}

impl RenderTileLayout {
    /// Plan row-major tiles and reject images outside tracer-v1's frozen
    /// `u32` pixel identity domain.
    pub fn try_new(
        image_width: u32,
        image_height: u32,
        tile_width: u32,
        tile_height: u32,
    ) -> Result<Self, RenderExecutionConfigError> {
        if image_width == 0 || image_height == 0 || image_width.checked_mul(image_height).is_none()
        {
            return Err(RenderExecutionConfigError::InvalidImageDimensions);
        }
        if tile_width == 0
            || tile_height == 0
            || tile_width > MAX_RENDER_TILE_EDGE
            || tile_height > MAX_RENDER_TILE_EDGE
        {
            return Err(RenderExecutionConfigError::InvalidTileShape);
        }
        let tiles_x = 1 + (image_width - 1) / tile_width;
        let tiles_y = 1 + (image_height - 1) / tile_height;
        let tile_count = u64::from(tiles_x) * u64::from(tiles_y);
        Ok(Self {
            image_width,
            image_height,
            tile_width,
            tile_height,
            tiles_x,
            tiles_y,
            tile_count,
        })
    }

    /// Number of logical tiles.
    #[must_use]
    pub const fn tile_count(self) -> u64 {
        self.tile_count
    }

    /// Image width bound into this layout.
    #[must_use]
    pub const fn image_width(self) -> u32 {
        self.image_width
    }

    /// Image height bound into this layout.
    #[must_use]
    pub const fn image_height(self) -> u32 {
        self.image_height
    }

    /// Unclipped logical tile width.
    #[must_use]
    pub const fn tile_width(self) -> u32 {
        self.tile_width
    }

    /// Unclipped logical tile height.
    #[must_use]
    pub const fn tile_height(self) -> u32 {
        self.tile_height
    }

    /// Number of tile columns.
    #[must_use]
    pub const fn tiles_x(self) -> u32 {
        self.tiles_x
    }

    /// Number of tile rows.
    #[must_use]
    pub const fn tiles_y(self) -> u32 {
        self.tiles_y
    }

    /// Bounds for a logical row-major tile ID.
    #[must_use]
    pub fn bounds(self, tile: u64) -> Option<RenderTileBounds> {
        if tile >= self.tile_count {
            return None;
        }
        let tile = u32::try_from(tile).ok()?;
        let tile_x = tile % self.tiles_x;
        let tile_y = tile / self.tiles_x;
        let x = tile_x * self.tile_width;
        let y = tile_y * self.tile_height;
        Some(RenderTileBounds {
            x,
            y,
            width: self.tile_width.min(self.image_width - x),
            height: self.tile_height.min(self.image_height - y),
        })
    }

    /// Exact pixel count in a half-open contiguous row-major tile range.
    ///
    /// This is constant-time even for a layout containing billions of tiny
    /// tiles, so resource admission never has to enumerate untrusted work.
    #[must_use]
    pub fn pixel_count_in_range(self, start: u64, end: u64) -> Option<u64> {
        if start > end || end > self.tile_count {
            return None;
        }
        let start_pixels = self.pixel_prefix(start)?;
        let end_pixels = self.pixel_prefix(end)?;
        end_pixels.checked_sub(start_pixels)
    }

    fn pixel_prefix(self, tile: u64) -> Option<u64> {
        if tile > self.tile_count {
            return None;
        }
        let tiles_x = u64::from(self.tiles_x);
        let row = tile / tiles_x;
        let column = tile % tiles_x;
        let row_origin = row.checked_mul(u64::from(self.tile_height))?;
        let completed_height = row_origin.min(u64::from(self.image_height));
        let completed_rows = completed_height.checked_mul(u64::from(self.image_width))?;
        let current_height = u64::from(self.image_height)
            .saturating_sub(row_origin)
            .min(u64::from(self.tile_height));
        let current_width = column
            .checked_mul(u64::from(self.tile_width))?
            .min(u64::from(self.image_width));
        completed_rows.checked_add(current_height.checked_mul(current_width)?)
    }
}

/// Structured failure of an explicit tile-render operation. Tracer/domain
/// failures remain distinct from memory admission and executor containment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RenderExecutionError {
    /// Invalid explicit policy or image layout.
    Config(RenderExecutionConfigError),
    /// Geometry, lighting, camera, or cancellation failure from the tracer.
    Tracer(TracerError),
    /// The hard operation-memory lease refused retained or staging storage.
    Memory(LeaseRefusal),
    /// A tile-local admitted allocation failed.
    Allocation(AllocError),
    /// Adaptive-policy validation or raw moment accumulation refused.
    Adaptive(AdaptiveSamplingError),
    /// The throughput lane contained a panic, spawn failure, or internal
    /// execution refusal and drained all children before returning.
    Executor(RunError),
    /// A defensive renderer invariant failed without publishing a film.
    Internal(&'static str),
}

impl core::fmt::Display for RenderExecutionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Config(error) => write!(formatter, "tile render configuration refused: {error}"),
            Self::Tracer(error) => write!(formatter, "tile render refused: {error}"),
            Self::Memory(error) => write!(formatter, "tile render memory refused: {error}"),
            Self::Allocation(error) => write!(formatter, "tile render allocation refused: {error}"),
            Self::Adaptive(error) => write!(formatter, "adaptive render refused: {error}"),
            Self::Executor(error) => write!(formatter, "tile render execution failed: {error}"),
            Self::Internal(detail) => write!(formatter, "tile render invariant failed: {detail}"),
        }
    }
}

impl core::error::Error for RenderExecutionError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Config(error) => Some(error),
            Self::Tracer(error) => Some(error),
            Self::Allocation(error) => Some(error),
            Self::Adaptive(error) => Some(error),
            Self::Executor(error) => Some(error),
            Self::Memory(_) | Self::Internal(_) => None,
        }
    }
}

impl From<TracerError> for RenderExecutionError {
    fn from(error: TracerError) -> Self {
        Self::Tracer(error)
    }
}

impl From<AdaptiveSamplingError> for RenderExecutionError {
    fn from(error: AdaptiveSamplingError) -> Self {
        Self::Adaptive(error)
    }
}

impl From<Cancelled> for RenderExecutionError {
    fn from(_: Cancelled) -> Self {
        Self::Tracer(TracerError::Cancelled)
    }
}

/// Structured failure of deterministic tile-parallel cinematic AOV rendering.
/// AOV alignment and profile refusals remain distinct from executor, memory,
/// and tracer execution failures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicAovExecutionError {
    /// AOV admission, alignment, accumulation, or artifact-state refusal.
    Aov(CinematicAovError),
    /// Tile policy, operation-memory, tracer, or executor refusal.
    Execution(RenderExecutionError),
}

impl core::fmt::Display for CinematicAovExecutionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Aov(error) => write!(formatter, "parallel cinematic AOV render refused: {error}"),
            Self::Execution(error) => error.fmt(formatter),
        }
    }
}

impl core::error::Error for CinematicAovExecutionError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Aov(error) => Some(error),
            Self::Execution(error) => Some(error),
        }
    }
}

impl From<CinematicAovError> for CinematicAovExecutionError {
    fn from(error: CinematicAovError) -> Self {
        Self::Aov(error)
    }
}

impl From<RenderExecutionError> for CinematicAovExecutionError {
    fn from(error: RenderExecutionError) -> Self {
        Self::Execution(error)
    }
}

/// Actionable measured report for one tile-render run. Timing fields are
/// diagnostics only and never enter image semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RenderExecutionReport {
    /// Fixed logical tile layout.
    pub layout: RenderTileLayout,
    /// Caller-requested worker ceiling.
    pub requested_workers: usize,
    /// Workers actually admitted for this tile count. The throughput lane
    /// never launches more workers than logical tiles.
    pub workers: usize,
    /// One-based attempt number for resumable jobs. Compatibility one-shot
    /// runs use one; no-work reports use zero.
    pub attempt_index: u64,
    /// Bytes charged for an already-retained progressive input film. Zero for
    /// fresh rendering and empty progressive ranges.
    pub retained_film_bytes: u64,
    /// Bytes charged for private all-or-nothing output state. A fresh uniform
    /// render owns one XYZ film payload; an adaptive render additionally owns
    /// its Welford, count, and decision AOVs.
    pub staging_film_bytes: u64,
    /// Worst-case concurrent tile/row pixel-accumulator scratch reserved before
    /// dispatch.
    pub tile_scratch_envelope_bytes: u64,
    /// Shared Sobol direction-state payload charged for this run; zero for IID.
    pub sampler_state_bytes: u64,
    /// Retained per-tile row-prefix checkpoint payload. Zero for the
    /// compatibility all-or-nothing APIs.
    pub progress_state_bytes: u64,
    /// Setup and admission wall time. Resumable jobs retain their original
    /// admission time in every later attempt report.
    pub setup_ns: u64,
    /// Throughput-lane wall time, including drain on failure.
    pub traversal_ns: u64,
    /// Sum of per-tile compute time across workers.
    pub tile_compute_ns: u64,
    /// Sum of per-tile staging-copy time across workers.
    pub tile_merge_ns: u64,
    /// Final failure-free film publication time.
    pub publication_ns: u64,
    /// Conservative worker-capacity time not accounted for by tile compute or
    /// staging copies. Scheduling and measurement overhead are included.
    pub idle_worker_ns: u64,
    /// Executor scheduling, completion, and cancellation diagnostics.
    pub executor: RunReport,
    /// Exact operation-memory admission trace after transient charges release.
    /// For a resumable job, counters and peak are cumulative from admission
    /// through this attempt. Executor and traversal/compute/merge/publication
    /// timings describe this attempt; `setup_ns` remains the initial resumable
    /// job admission time.
    pub memory: LeaseReceipt,
}

/// Fresh film plus the execution evidence that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct RenderExecutionOutput {
    /// Fully published film.
    pub film: Film,
    /// Tile, scheduling, timing, and memory evidence.
    pub report: RenderExecutionReport,
}

/// Fully published aligned cinematic AOV film plus deterministic tile
/// execution evidence.
#[derive(Debug, PartialEq)]
pub struct CinematicAovExecutionOutput {
    /// Complete private-until-success beauty and aligned AOV state.
    pub film: CinematicAovFilm,
    /// Tile, scheduling, timing, and memory evidence.
    pub report: RenderExecutionReport,
}

/// Fully published adaptive cinematic AOV film plus the deterministic tile
/// execution evidence.  Each AOV sample is aligned to that pixel's terminal
/// Welford-selected prefix; `report` is observational and does not enter the
/// raw estimator or adaptive decision semantics.
#[derive(Debug, PartialEq)]
pub struct AdaptiveCinematicAovExecutionOutput {
    /// Complete private-until-success adaptive beauty and aligned AOV state.
    pub film: AdaptiveCinematicAovFilm,
    /// Tile, scheduling, timing, and memory evidence.
    pub report: RenderExecutionReport,
}

/// Adaptive film plus execution evidence and exact raw path-count summary.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveRenderOutput {
    /// Fully published adaptive film and statistical AOVs.
    pub film: AdaptiveFilm,
    /// Tile, scheduling, timing, and memory evidence.
    pub report: RenderExecutionReport,
}

impl AdaptiveRenderOutput {
    /// Deterministic aggregate of the authoritative film sample-count and
    /// decision maps.
    #[must_use]
    pub fn summary(&self) -> AdaptiveRenderSummary {
        self.film.summary()
    }
}

/// Public progress of one opaque in-memory render job. Counts never expose a
/// partially updated film.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RenderProgress {
    /// Tile rows committed atomically into the private eventual film.
    pub committed_tile_rows: u64,
    /// Total tile rows in the fixed logical layout.
    pub total_tile_rows: u64,
    /// Tiles whose every row is committed.
    pub completed_tiles: u64,
    /// Total logical tiles.
    pub total_tiles: u64,
    /// Execution API attempts already made, including failed, cancelled, and
    /// completed-job calls that require no worker dispatch.
    pub attempts: u64,
}

struct PendingRenderState {
    xyz: Vec<[f64; 3]>,
    next_row: Vec<u32>,
}

struct PendingAdaptiveRenderState {
    film: AdaptiveRenderState,
    next_row: Vec<u32>,
}

/// Opaque, single-buffer fresh render that can retain completed row prefixes
/// across cancellation or a contained worker failure.
///
/// The job borrows its exact scene/camera assets and owns its settings,
/// execution policy, mode, and budget, eventual film buffer, and row-prefix
/// checkpoints. Resume therefore cannot substitute a different scene, sample
/// range, layout, mode, budget, or `RunId`. The only public image is returned
/// after complete success.
#[must_use = "resume or retain the pending render; dropping it discards private progress"]
pub struct PendingRender<'assets> {
    scene: &'assets Scene,
    lighting: AdmittedLighting<'assets>,
    settings: Settings,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: RenderExecutionConfig,
    layout: RenderTileLayout,
    state: Mutex<PendingRenderState>,
    sobol: Option<Sobol>,
    lease: OperationMemoryLease,
    film_charge: Option<LeaseCharge>,
    progress_charge: Option<LeaseCharge>,
    sampler_charge: Option<LeaseCharge>,
    film_bytes: u64,
    progress_state_bytes: u64,
    sampler_state_bytes: u64,
    setup_ns: u64,
    attempts: u64,
}

/// Failed/cancelled attempt that retains the opaque job and its committed row
/// prefixes for an exact retry under a fresh `Cx` cancellation authority.
#[must_use = "inspect or resume the suspension; dropping it discards private progress"]
pub struct RenderSuspension<'assets> {
    work: PendingRender<'assets>,
    cause: RenderExecutionError,
    attempt: RenderExecutionReport,
}

/// Successful bounded render attempt that intentionally yielded at a
/// row-atomic safe point without publishing a partial film.
#[must_use = "checkpoint or resume the yielded render; dropping it discards private progress"]
pub struct RenderCheckpointYield<'assets> {
    work: PendingRender<'assets>,
    attempt: RenderExecutionReport,
}

/// Opaque adaptive render that retains only complete tile-row prefixes across
/// cancellation or contained worker failure.
///
/// The exact scene/camera assets, settings, adaptive policy, estimator
/// version, sampler stream, execution mode, budget, layout, and `RunId` are
/// bound at construction. No partial row or public film is observable.
#[must_use = "resume or retain the pending adaptive render; dropping it discards private progress"]
pub struct PendingAdaptiveRender<'assets> {
    scene: &'assets Scene,
    lighting: AdmittedLighting<'assets>,
    settings: Settings,
    policy: AdaptiveSamplingConfig,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    requested_mode: FilmTimeMode,
    execution_mode: ExecMode,
    execution_budget: Budget,
    execution: RenderExecutionConfig,
    layout: RenderTileLayout,
    state: Mutex<PendingAdaptiveRenderState>,
    sobol: Option<Sobol>,
    lease: OperationMemoryLease,
    state_charge: Option<LeaseCharge>,
    progress_charge: Option<LeaseCharge>,
    sampler_charge: Option<LeaseCharge>,
    state_bytes: u64,
    progress_state_bytes: u64,
    sampler_state_bytes: u64,
    setup_ns: u64,
    attempts: u64,
}

/// Failed/cancelled adaptive attempt retaining the exact private job for a
/// later retry.
#[must_use = "inspect or resume the adaptive suspension; dropping it discards private progress"]
pub struct AdaptiveRenderSuspension<'assets> {
    work: PendingAdaptiveRender<'assets>,
    cause: RenderExecutionError,
    attempt: RenderExecutionReport,
}

/// Successful bounded adaptive attempt that intentionally yielded at a
/// row-atomic safe point without publishing partial estimator AOVs.
#[must_use = "checkpoint or resume the yielded adaptive render; dropping it discards private progress"]
pub struct AdaptiveRenderCheckpointYield<'assets> {
    work: PendingAdaptiveRender<'assets>,
    attempt: RenderExecutionReport,
}

impl core::fmt::Debug for PendingRender<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("PendingRender")
            .field("settings", &self.settings)
            .field("execution", &self.execution)
            .field("layout", &self.layout)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for RenderSuspension<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("RenderSuspension")
            .field("cause", &self.cause)
            .field("attempt", &self.attempt)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for RenderCheckpointYield<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("RenderCheckpointYield")
            .field("attempt", &self.attempt)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for PendingAdaptiveRender<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("PendingAdaptiveRender")
            .field("settings", &self.settings)
            .field("policy", &self.policy)
            .field("execution", &self.execution)
            .field("layout", &self.layout)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for AdaptiveRenderSuspension<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("AdaptiveRenderSuspension")
            .field("cause", &self.cause)
            .field("attempt", &self.attempt)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl core::fmt::Debug for AdaptiveRenderCheckpointYield<'_> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_struct("AdaptiveRenderCheckpointYield")
            .field("attempt", &self.attempt)
            .field("progress", &self.progress())
            .finish_non_exhaustive()
    }
}

impl<'assets> RenderSuspension<'assets> {
    /// Structured cause of the most recent attempt.
    #[must_use]
    pub const fn cause(&self) -> &RenderExecutionError {
        &self.cause
    }

    /// Drain, timing, and progress evidence for the most recent attempt.
    /// Memory counters are job-cumulative through that attempt, and a
    /// suspended job legitimately retains nonzero lease bytes.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        &self.attempt
    }

    /// Current opaque progress without exposing film pixels.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.work.progress()
    }

    /// Recover the pending job for a later retry.
    #[must_use]
    pub fn into_pending(self) -> PendingRender<'assets> {
        self.work
    }
}

impl<'assets> RenderCheckpointYield<'assets> {
    /// Drain/timing evidence for the successful bounded attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        &self.attempt
    }

    /// Row-atomic progress retained by the opaque pending job.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.work.progress()
    }

    /// Recover the pending job for durable checkpointing or another attempt.
    #[must_use]
    pub fn into_pending(self) -> PendingRender<'assets> {
        self.work
    }
}

impl<'assets> AdaptiveRenderSuspension<'assets> {
    /// Structured cause of the most recent adaptive attempt.
    #[must_use]
    pub const fn cause(&self) -> &RenderExecutionError {
        &self.cause
    }

    /// Drain, timing, and cumulative memory evidence for the last attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        &self.attempt
    }

    /// Current opaque row/tile progress.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.work.progress()
    }

    /// Recover the exact pending adaptive job for retry.
    #[must_use]
    pub fn into_pending(self) -> PendingAdaptiveRender<'assets> {
        self.work
    }
}

impl<'assets> AdaptiveRenderCheckpointYield<'assets> {
    /// Drain/timing evidence for the successful bounded adaptive attempt.
    #[must_use]
    pub const fn attempt_report(&self) -> &RenderExecutionReport {
        &self.attempt
    }

    /// Row-atomic progress retained by the opaque adaptive job.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        self.work.progress()
    }

    /// Recover the pending job for durable checkpointing or another attempt.
    #[must_use]
    pub fn into_pending(self) -> PendingAdaptiveRender<'assets> {
        self.work
    }
}

/// Fail-closed spectral-tracer diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TracerError {
    /// The supplied execution context requested cancellation.
    Cancelled,
    /// A chart backend stopped in a state other than a clean miss or certified
    /// residual hit.
    BackendFailure(TraceTermination),
    /// A chart returned a terminal result without retaining its typed
    /// no-tunneling claim. Uncertified misses are not geometry absence.
    UncertifiedTrace,
    /// A progressive sample range had its exclusive end before its start.
    InvalidRange {
        /// Inclusive start of the rejected sample range.
        from: u32,
        /// Exclusive end of the rejected sample range.
        to: u32,
    },
    /// Render dimensions or a public film buffer were structurally invalid or
    /// could not be allocated without exceeding address-space bounds.
    InvalidInput,
    /// Shading requires a finite surface normal; no arbitrary fallback normal
    /// may be minted.
    MissingNormal,
    /// Instance construction, placement, hit data, or object IDs were invalid.
    InvalidInstance,
    /// The accepted instance hit could not form a valid stable surface record.
    MotionVector(MotionVectorError),
    /// An animated instance was supplied to a static render entry point.
    MissingRayTime,
    /// The resolved shutter extends beyond an animated instance trajectory.
    MotionOutsideTrajectory,
    /// A progressive append attempted to mix static, legacy-motion, or
    /// cinematic samples, or changed the admitted shutter, stream, or shot.
    ProgressiveTimeModeMismatch,
    /// The v1 NEE light is static metadata and cannot name animated geometry.
    AnimatedLightUnsupported,
    /// Rectangular light metadata did not name a matching emissive primitive.
    LightPrimitiveMismatch {
        /// Primitive index named by the rejected light.
        light_primitive: usize,
    },
    /// A cinematic camera was malformed, evaluated outside its admitted shot,
    /// or cancelled. The nested error retains ranked admission fixes.
    Camera(CameraError),
    /// Validated dielectric evaluation unexpectedly refused.
    Dielectric(DielectricError),
    /// Validated conductor evaluation unexpectedly refused.
    Conductor(ConductorError),
    /// Lighting artifact, selection, or probability admission refused.
    Lighting(LightingError),
    /// An encountered dielectric boundary violated strict LIFO nesting or was
    /// oriented as an exit while the path remained in ambient air.
    MediumStackMismatch {
        /// Boundary primitive being processed.
        boundary_primitive: usize,
        /// Active top boundary, if any.
        active_boundary: Option<usize>,
    },
    /// The defensive nested-medium ceiling was exceeded.
    MediumStackOverflow,
    /// A direct-light sample first encountered a dielectric, but the geometry
    /// or material was outside the exact smooth homogeneous parallel-slab
    /// estimator.  Refusing is required: treating an unsupported transparent
    /// blocker as either opaque or undeviating would bias the raw estimator.
    UnsupportedSlabNee {
        /// First dielectric primitive on the attempted connection.
        boundary_primitive: usize,
        /// Stable, actionable refusal class.
        reason: &'static str,
    },
    /// A ray missed all geometry while still inside a declared closed medium.
    UnclosedMedium {
        /// Active top boundary at the miss.
        boundary_primitive: usize,
        /// Integrator operation whose ray escaped the declared closed medium.
        context: &'static str,
        /// Film pixel's deterministic linear index.
        pixel: u32,
        /// Deterministic sample index within the pixel.
        sample: u32,
        /// Path depth at which the invalid miss occurred.
        depth: u32,
        /// Active spectral lane after a dispersive split, if any.
        active_lane: Option<usize>,
        /// Exact IEEE-754 bits of the escaped ray origin.
        origin_bits: [u64; 3],
        /// Exact IEEE-754 bits of the escaped ray direction.
        direction_bits: [u64; 3],
        /// Geometric surface normal at a next-event shadow-ray spawn, if any.
        geometric_normal_bits: Option<[u64; 3]>,
        /// Incident path direction at a next-event shadow-ray spawn, if any.
        incident_direction_bits: Option<[u64; 3]>,
        /// Oriented-normal dot shadow-direction bits, if applicable.
        shadow_side_dot_bits: Option<u64>,
    },
}

impl core::fmt::Display for TracerError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("spectral render cancelled"),
            Self::BackendFailure(termination) => {
                write!(formatter, "chart backend stopped with {termination:?}")
            }
            Self::UncertifiedTrace => {
                formatter.write_str("chart backend produced an uncertified trace result")
            }
            Self::InvalidRange { from, to } => {
                write!(formatter, "invalid progressive sample range {from}..{to}")
            }
            Self::InvalidInput => formatter.write_str("invalid spectral render input"),
            Self::MissingNormal => formatter.write_str("surface hit has no finite normal"),
            Self::InvalidInstance => formatter.write_str("invalid rigid render instance"),
            Self::MotionVector(error) => write!(
                formatter,
                "primary motion/AOV surface record refused: {error}"
            ),
            Self::MissingRayTime => {
                formatter.write_str("animated render instance requires explicit shutter time")
            }
            Self::MotionOutsideTrajectory => {
                formatter.write_str("render shutter lies outside an instance trajectory")
            }
            Self::ProgressiveTimeModeMismatch => {
                formatter.write_str("progressive film camera/time checkpoint mismatch")
            }
            Self::AnimatedLightUnsupported => {
                formatter.write_str("animated area-light geometry is unsupported")
            }
            Self::LightPrimitiveMismatch { light_primitive } => write!(
                formatter,
                "area light does not match emissive primitive {light_primitive}"
            ),
            Self::Camera(error) => write!(formatter, "cinematic camera refused: {error}"),
            Self::Dielectric(error) => write!(formatter, "dielectric transport refused: {error}"),
            Self::Conductor(error) => write!(formatter, "conductor transport refused: {error}"),
            Self::Lighting(error) => write!(formatter, "scene lighting refused: {error}"),
            Self::MediumStackMismatch {
                boundary_primitive,
                active_boundary,
            } => write!(
                formatter,
                "dielectric boundary {boundary_primitive} violated LIFO nesting; active boundary {active_boundary:?}"
            ),
            Self::MediumStackOverflow => formatter.write_str("dielectric medium stack overflow"),
            Self::UnsupportedSlabNee {
                boundary_primitive,
                reason,
            } => write!(
                formatter,
                "smooth parallel-slab next-event connection through dielectric boundary {boundary_primitive} refused: {reason}"
            ),
            Self::UnclosedMedium {
                boundary_primitive,
                context,
                pixel,
                sample,
                depth,
                active_lane,
                origin_bits,
                direction_bits,
                geometric_normal_bits,
                incident_direction_bits,
                shadow_side_dot_bits,
            } => {
                let origin = origin_bits.map(f64::from_bits);
                let direction = direction_bits.map(f64::from_bits);
                let geometric_normal = geometric_normal_bits.map(|bits| bits.map(f64::from_bits));
                let incident_direction =
                    incident_direction_bits.map(|bits| bits.map(f64::from_bits));
                let shadow_side_dot = shadow_side_dot_bits.map(f64::from_bits);
                write!(
                    formatter,
                    "ray missed while still inside dielectric boundary {boundary_primitive}; context={context} pixel={pixel} sample={sample} depth={depth} active_lane={active_lane:?} origin={origin:?} direction={direction:?} geometric_normal={geometric_normal:?} incident_direction={incident_direction:?} shadow_side_dot={shadow_side_dot:?}"
                )
            }
        }
    }
}

impl core::error::Error for TracerError {}

impl From<Cancelled> for TracerError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<InstanceError> for TracerError {
    fn from(error: InstanceError) -> Self {
        match error {
            InstanceError::Cancelled => Self::Cancelled,
            InstanceError::BackendFailure(termination) => Self::BackendFailure(termination),
            InstanceError::UncertifiedTrace => Self::UncertifiedTrace,
            InstanceError::MissingNormal => Self::MissingNormal,
            InstanceError::InvalidTransform
            | InstanceError::InvalidObjectId
            | InstanceError::InvalidGeometryIdentity
            | InstanceError::DuplicateObjectId
            | InstanceError::TooManyInstances
            | InstanceError::InvalidIntersectionInput
            | InstanceError::InvalidHit => Self::InvalidInstance,
        }
    }
}

impl From<AnimatedInstanceError> for TracerError {
    fn from(error: AnimatedInstanceError) -> Self {
        match error {
            AnimatedInstanceError::Cancelled => Self::Cancelled,
            AnimatedInstanceError::ShutterOutsideTrajectory
            | AnimatedInstanceError::Extrapolation => Self::MotionOutsideTrajectory,
            AnimatedInstanceError::Instance(error) => error.into(),
            AnimatedInstanceError::EmptyTrajectory
            | AnimatedInstanceError::InvalidKeyframeTime
            | AnimatedInstanceError::InvalidKeyframeVelocity
            | AnimatedInstanceError::NonIncreasingKeyframeTime
            | AnimatedInstanceError::InvalidEvaluationTime
            | AnimatedInstanceError::InvalidInterpolation => Self::InvalidInstance,
        }
    }
}

impl From<MotionVectorError> for TracerError {
    fn from(error: MotionVectorError) -> Self {
        match error {
            MotionVectorError::Cancelled => Self::Cancelled,
            other => Self::MotionVector(other),
        }
    }
}

impl From<CameraError> for TracerError {
    fn from(error: CameraError) -> Self {
        match error {
            CameraError::Cancelled => Self::Cancelled,
            other => Self::Camera(other),
        }
    }
}

impl From<DielectricError> for TracerError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}

impl From<ConductorError> for TracerError {
    fn from(error: ConductorError) -> Self {
        Self::Conductor(error)
    }
}

impl From<LightingError> for TracerError {
    fn from(error: LightingError) -> Self {
        match error {
            LightingError::Cancelled => Self::Cancelled,
            other => Self::Lighting(other),
        }
    }
}

/// Time semantics already committed to a progressive film.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilmTimeMode {
    /// No non-empty sample range has committed yet.
    Uninitialized,
    /// Samples came from the legacy static render path.
    Static,
    /// Samples came from one exact admitted shutter and render stream.
    Motion {
        /// Complete shutter definition used by every committed path.
        shutter: ShutterInterval,
        /// Render seed that domain-separated the shutter-time stream.
        stream_identity: u64,
    },
    /// Samples came from one admitted cinematic-camera shot. The stable shot
    /// ID prevents progressive appends from crossing a hard-cut side or
    /// silently switching between legacy and cinematic ray generation.
    Cinematic {
        /// Complete shutter definition used by every committed path.
        shutter: ShutterInterval,
        /// Render seed that domain-separated shutter and lens streams.
        stream_identity: u64,
        /// Stable nonzero identity of the shot that owns the exposure.
        shot_id: u64,
    },
}

/// Accumulated CIE XYZ film: `spp` samples summed per pixel (divide on
/// output). Checkpointable: rendering samples `[a, b)` then `[b, c)`
/// into the same film equals rendering `[a, c)` bitwise.
#[derive(Debug, Clone, PartialEq)]
pub struct Film {
    /// Width.
    pub width: u32,
    /// Height.
    pub height: u32,
    /// Row-major XYZ sums.
    pub xyz: Vec<[f64; 3]>,
    /// Samples accumulated so far.
    pub spp_done: u32,
    /// Camera path and exact-shutter provenance required for the next append.
    pub time_mode: FilmTimeMode,
}

impl Film {
    /// An empty film.
    #[must_use]
    pub fn new(width: u32, height: u32) -> Film {
        Self::try_new(width, height).expect("film dimensions must fit the address space")
    }

    /// Allocate an empty film without panicking on malformed or excessive
    /// dimensions.
    pub fn try_new(width: u32, height: u32) -> Result<Film, TracerError> {
        let len = checked_pixel_len(width, height)?;
        let mut xyz = Vec::new();
        xyz.try_reserve_exact(len)
            .map_err(|_| TracerError::InvalidInput)?;
        xyz.resize(len, [0.0; 3]);
        Ok(Film {
            width,
            height,
            xyz,
            spp_done: 0,
            time_mode: FilmTimeMode::Uninitialized,
        })
    }

    /// Linear-sRGB planes (R, G, B row-major), Bradford-adapted like
    /// the rest of the spectral pipeline; sums divided by `spp_done`.
    #[must_use]
    pub fn to_linear_srgb(&self) -> [Vec<f32>; 3] {
        let n = self.xyz.len();
        let mut planes = [vec![0.0f32; n], vec![0.0f32; n], vec![0.0f32; n]];
        let inv = if self.spp_done == 0 {
            0.0
        } else {
            1.0 / f64::from(self.spp_done)
        };
        for (i, xyz) in self.xyz.iter().enumerate() {
            let rgb = xyz_to_linear_srgb(xyz_e_to_d65([xyz[0] * inv, xyz[1] * inv, xyz[2] * inv]));
            for (p, v) in planes.iter_mut().zip(rgb) {
                #[allow(clippy::cast_possible_truncation)]
                {
                    p[i] = v as f32;
                }
            }
        }
        planes
    }
}

/// Raw adaptive film and its per-pixel statistical AOVs.
///
/// `xyz` remains an unnormalised path sum, exactly like [`Film::xyz`], but
/// every pixel has its own divisor in `sample_counts`. The second central
/// moment is retained so checkpointing and later adaptive batches never need
/// to infer variance from quantised output pixels. The retained sampler,
/// policy, seed, and time mode are diagnostic estimator provenance, not a
/// complete replay identity: this film does not bind the full [`Settings`] or
/// scene content.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveFilm {
    width: u32,
    height: u32,
    xyz: Vec<[f64; 3]>,
    mean_xyz_aov: Vec<[f64; 3]>,
    m2_xyz: Vec<[f64; 3]>,
    sample_counts: Vec<u32>,
    decisions: Vec<AdaptiveDecision>,
    maximum_samples: u32,
    policy: AdaptiveSamplingConfig,
    sampler: Sampler,
    stream_seed: u64,
    semantics_version: u32,
    time_mode: FilmTimeMode,
}

impl AdaptiveFilm {
    /// Width in pixels.
    #[must_use]
    pub const fn width(&self) -> u32 {
        self.width
    }

    /// Height in pixels.
    #[must_use]
    pub const fn height(&self) -> u32 {
        self.height
    }

    /// Row-major raw XYZ sums. These preserve the existing uniform tracer's
    /// per-pixel sequential addition order.
    #[must_use]
    pub fn xyz_sums(&self) -> &[[f64; 3]] {
        &self.xyz
    }

    /// Row-major Welford means retained for stable exact resume.
    #[must_use]
    pub fn running_means_xyz(&self) -> &[[f64; 3]] {
        &self.mean_xyz_aov
    }

    /// Row-major Welford second central moments.
    #[must_use]
    pub fn m2_xyz(&self) -> &[[f64; 3]] {
        &self.m2_xyz
    }

    /// Row-major number of raw path samples actually consumed.
    #[must_use]
    pub fn sample_counts(&self) -> &[u32] {
        &self.sample_counts
    }

    /// Row-major deterministic stopping decisions.
    #[must_use]
    pub fn decisions(&self) -> &[AdaptiveDecision] {
        &self.decisions
    }

    /// Hard sample ceiling supplied through [`Settings::spp`].
    #[must_use]
    pub const fn maximum_samples(&self) -> u32 {
        self.maximum_samples
    }

    /// Exact adaptive policy bound into this film's decisions.
    #[must_use]
    pub const fn policy(&self) -> AdaptiveSamplingConfig {
        self.policy
    }

    /// Raw sample-sequence family used by the estimator.
    #[must_use]
    pub const fn sampler(&self) -> Sampler {
        self.sampler
    }

    /// Replay seed for absolute `(pixel, sample, dimension)` streams.
    #[must_use]
    pub const fn stream_seed(&self) -> u64 {
        self.stream_seed
    }

    /// Bit-affecting adaptive estimator/stopping version.
    #[must_use]
    pub const fn semantics_version(&self) -> u32 {
        self.semantics_version
    }

    /// Camera path and exact-shutter provenance shared by every sample.
    #[must_use]
    pub const fn time_mode(&self) -> FilmTimeMode {
        self.time_mode
    }

    /// Per-pixel Welford running mean used by the stopping estimator.
    ///
    /// This value is retained independently from the raw beauty sum. Call
    /// [`Self::beauty_mean_xyz`] when converting or comparing rendered
    /// radiance so the existing uniform tracer's sequential summation order
    /// remains the oracle.
    #[must_use]
    pub fn estimator_mean_xyz(&self, pixel: usize) -> Option<[f64; 3]> {
        self.mean_xyz_aov.get(pixel).copied()
    }

    /// Per-pixel beauty mean computed from the raw sequential sum and the
    /// exact number of traced paths.
    #[must_use]
    pub fn beauty_mean_xyz(&self, pixel: usize) -> Option<[f64; 3]> {
        let sum = *self.xyz.get(pixel)?;
        let samples = *self.sample_counts.get(pixel)?;
        if samples == 0 {
            Some([0.0; 3])
        } else {
            let inverse = 1.0 / f64::from(samples);
            Some(sum.map(|value| value * inverse))
        }
    }

    /// Bessel-corrected per-channel sample variance for one pixel. At a fixed
    /// IID sample count this is unbiased; adaptive stopping removes that
    /// general guarantee. For Owen-Sobol it is only within-stream dispersion.
    /// Before two samples, the estimator is undefined and this returns
    /// positive infinity.
    #[must_use]
    pub fn sample_variance_xyz(&self, pixel: usize) -> Option<[f64; 3]> {
        let m2 = *self.m2_xyz.get(pixel)?;
        let samples = *self.sample_counts.get(pixel)?;
        if m2.iter().any(|value| !value.is_finite() || *value < 0.0) {
            return None;
        }
        if samples < 2 {
            Some([f64::INFINITY; 3])
        } else {
            let inverse = 1.0 / f64::from(samples - 1);
            Some(m2.map(|value| value * inverse))
        }
    }

    /// IID standard-error estimate, or a within-stream dispersion proxy for
    /// Owen-Sobol. This is not a confidence interval or image-error
    /// certificate.
    #[must_use]
    pub fn dispersion_proxy_xyz(&self, pixel: usize) -> Option<[f64; 3]> {
        let samples = *self.sample_counts.get(pixel)?;
        let inverse_samples = 1.0 / f64::from(samples.max(1));
        Some(
            self.sample_variance_xyz(pixel)?
                .map(|variance| det::sqrt(variance * inverse_samples)),
        )
    }

    /// Linear-sRGB planes normalised by each pixel's actual raw sample count.
    #[must_use]
    pub fn to_linear_srgb(&self) -> [Vec<f32>; 3] {
        let n = self.xyz.len();
        let mut planes = [vec![0.0f32; n], vec![0.0f32; n], vec![0.0f32; n]];
        for pixel in 0..n {
            let xyz = self
                .beauty_mean_xyz(pixel)
                .expect("adaptive film owns shape-matched private buffers");
            let rgb = xyz_to_linear_srgb(xyz_e_to_d65(xyz));
            for (plane, value) in planes.iter_mut().zip(rgb) {
                #[allow(clippy::cast_possible_truncation)]
                {
                    plane[pixel] = value as f32;
                }
            }
        }
        planes
    }

    /// Deterministic aggregate counters for progress, cost, and checkpoint
    /// metadata. This is a raw path-count summary, not a quality certificate.
    #[must_use]
    pub fn summary(&self) -> AdaptiveRenderSummary {
        let mut minimum_samples = u32::MAX;
        let mut maximum_samples = 0_u32;
        let mut total_samples = 0_u64;
        let mut converged_pixels = 0_u64;
        for pixel in 0..self.sample_counts.len() {
            let samples = self.sample_counts[pixel];
            let decision = self.decisions[pixel];
            minimum_samples = minimum_samples.min(samples);
            maximum_samples = maximum_samples.max(samples);
            total_samples += u64::from(samples);
            converged_pixels += u64::from(decision == AdaptiveDecision::ErrorThreshold);
        }
        if self.sample_counts.is_empty() {
            minimum_samples = 0;
        }
        AdaptiveRenderSummary {
            pixels: self.sample_counts.len() as u64,
            minimum_samples,
            maximum_samples,
            total_samples,
            converged_pixels,
            maximum_sample_pixels: self.sample_counts.len() as u64 - converged_pixels,
        }
    }

    /// Deterministic row-major aggregate for one logical tile. This reports a
    /// maximum dispersion proxy, never an average mislabeled as tile
    /// variance.
    #[must_use]
    pub fn tile_summary(&self, layout: RenderTileLayout, tile: u64) -> Option<AdaptiveTileSummary> {
        if (layout.image_width, layout.image_height) != (self.width, self.height) {
            return None;
        }
        let bounds = layout.bounds(tile)?;
        let mut minimum_samples = u32::MAX;
        let mut maximum_samples = 0_u32;
        let mut total_samples = 0_u64;
        let mut converged_pixels = 0_u64;
        let mut maximum_dispersion_xyz = [0.0_f64; 3];
        for y in bounds.y..bounds.y + bounds.height {
            for x in bounds.x..bounds.x + bounds.width {
                let pixel = y as usize * self.width as usize + x as usize;
                let samples = self.sample_counts[pixel];
                minimum_samples = minimum_samples.min(samples);
                maximum_samples = maximum_samples.max(samples);
                total_samples += u64::from(samples);
                converged_pixels +=
                    u64::from(self.decisions[pixel] == AdaptiveDecision::ErrorThreshold);
                let dispersion = self.dispersion_proxy_xyz(pixel)?;
                for channel in 0..3 {
                    maximum_dispersion_xyz[channel] =
                        maximum_dispersion_xyz[channel].max(dispersion[channel]);
                }
            }
        }
        let pixels = u64::from(bounds.width) * u64::from(bounds.height);
        Some(AdaptiveTileSummary {
            bounds,
            pixels,
            minimum_samples,
            maximum_samples,
            total_samples,
            converged_pixels,
            maximum_sample_pixels: pixels - converged_pixels,
            maximum_dispersion_xyz,
        })
    }
}

/// Exact aggregate of one adaptive film's decision map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdaptiveRenderSummary {
    /// Number of image pixels.
    pub pixels: u64,
    /// Smallest per-pixel sample count.
    pub minimum_samples: u32,
    /// Largest per-pixel sample count.
    pub maximum_samples: u32,
    /// Total raw paths consumed across all pixels.
    pub total_samples: u64,
    /// Pixels stopped by the declared error proxy.
    pub converged_pixels: u64,
    /// Pixels stopped only by the hard sample ceiling.
    pub maximum_sample_pixels: u64,
}

/// Deterministic raw-statistics aggregate for one logical image tile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveTileSummary {
    /// Exact tile bounds.
    pub bounds: RenderTileBounds,
    /// Number of pixels in the clipped tile.
    pub pixels: u64,
    /// Smallest per-pixel sample count in the tile.
    pub minimum_samples: u32,
    /// Largest per-pixel sample count in the tile.
    pub maximum_samples: u32,
    /// Exact raw paths consumed by the tile.
    pub total_samples: u64,
    /// Pixels stopped by the declared error proxy.
    pub converged_pixels: u64,
    /// Pixels stopped only by the hard ceiling.
    pub maximum_sample_pixels: u64,
    /// Per-channel maximum IID standard-error estimate or Owen-Sobol
    /// within-stream dispersion proxy.
    pub maximum_dispersion_xyz: [f64; 3],
}

/// Render samples `[from, to)` for every pixel into `film` (progressive
/// accumulation; `film.spp_done` must equal `from`).
///
/// This compatibility path is intentionally serial and serves as the bitwise
/// oracle for [`render_range_with_execution`]. Shape, checkpoint, and range
/// mismatches are returned as structured errors without changing `film`.
pub fn render_range(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    render_range_impl(scene, cx, s, film, from, to, None, CameraPath::Legacy)
}

/// Render a progressive sample range with one deterministic shutter-time draw
/// per camera path. Every secondary and shadow ray in that path retains the
/// same physical time. The operation is transactional just like
/// [`render_range`].
pub fn render_motion_range(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: ShutterInterval,
) -> Result<(), TracerError> {
    render_range_impl(
        scene,
        cx,
        s,
        film,
        from,
        to,
        Some(shutter),
        CameraPath::Legacy,
    )
}

/// Render a progressive range with a validated thin-lens/keyframed camera.
/// Camera and geometry share the one absolute shutter time drawn for each
/// path. The complete exposure must belong to one camera shot; a positive
/// shutter crossing a hard cut refuses before film publication.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_range(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    s: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: ShutterInterval,
) -> Result<(), TracerError> {
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    render_range_impl(
        scene,
        cx,
        s,
        film,
        from,
        to,
        Some(shutter),
        CameraPath::Cinematic { camera, exposure },
    )
}

#[allow(clippy::too_many_arguments)] // private seam keeps legacy public APIs unchanged
fn render_range_impl(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
) -> Result<(), TracerError> {
    let (lighting, requested_mode) =
        preflight_render(scene, cx, s, Some(film), from, to, shutter, camera_path)?;
    if to == from {
        return Ok(());
    }
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let sobol = (s.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, s.seed));
    let kn = 1.0 / y_integral();
    // Cancellation and backend refusals are transactional: a failed range
    // leaves both the accumulated sums and checkpoint unchanged, so retrying
    // cannot double-count a partially completed range.
    let mut staged_xyz = Vec::new();
    staged_xyz
        .try_reserve_exact(film.xyz.len())
        .map_err(|_| TracerError::InvalidInput)?;
    for chunk in film.xyz.chunks(4096) {
        cx.checkpoint()?;
        staged_xyz.extend_from_slice(chunk);
    }
    for py in 0..s.height {
        cx.checkpoint()?;
        for px in 0..s.width {
            let pixel = py * s.width + px;
            let slot = &mut staged_xyz[pixel as usize];
            for sample in from..to {
                cx.checkpoint()?;
                let xyz = trace_pixel_sample(
                    scene,
                    &lighting,
                    cx,
                    s,
                    kn,
                    sobol.as_ref(),
                    key,
                    pixel,
                    sample,
                    shutter,
                    camera_path,
                )?;
                slot[0] += xyz[0];
                slot[1] += xyz[1];
                slot[2] += xyz[2];
            }
        }
    }
    cx.checkpoint()?;
    film.xyz = staged_xyz;
    film.spp_done = to;
    if film.time_mode == FilmTimeMode::Uninitialized {
        film.time_mode = requested_mode;
    }
    Ok(())
}

#[derive(Debug, Clone)]
enum RenderTileFailure {
    Tracer(TracerError),
    Allocation(AllocError),
    Adaptive(AdaptiveSamplingError),
    Internal(&'static str),
}

struct ParallelRenderKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    settings: &'run Settings,
    base_xyz: Option<&'run [[f64; 3]]>,
    staging: &'run Mutex<Vec<[f64; 3]>>,
    failures: &'run Mutex<Option<(u64, RenderTileFailure)>>,
    layout: RenderTileLayout,
    from: u32,
    to: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl ParallelRenderKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: RenderTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("tile outside planned layout"),
            );
        };
        let Some(pixel_count) = bounds
            .width
            .checked_mul(bounds.height)
            .and_then(|count| usize::try_from(count).ok())
        else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("tile pixel count overflow"),
            );
        };
        let mut pixels = Vec::new();
        if pixels.try_reserve_exact(pixel_count).is_err() {
            return self.fail(
                tile,
                RenderTileFailure::Allocation(AllocError::OutOfMemory {
                    site: "render-tile-pixels",
                    requested_bytes: pixel_count.saturating_mul(size_of::<[f64; 3]>()),
                }),
            );
        }
        let compute_started = Instant::now();
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                return ControlFlow::Break(Cancelled);
            }
            for px in bounds.x..bounds.x + bounds.width {
                let Some(pixel) = py
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(px))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pixel identity overflow after preflight"),
                    );
                };
                let mut xyz = self.base_xyz.map_or([0.0; 3], |base| base[pixel as usize]);
                for sample in self.from..self.to {
                    if cx.checkpoint().is_err() {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    let sample_xyz = match trace_pixel_sample(
                        self.scene,
                        self.lighting,
                        cx,
                        self.settings,
                        kn,
                        self.sobol,
                        key,
                        pixel,
                        sample,
                        self.shutter,
                        self.camera_path,
                    ) {
                        Ok(sample_xyz) => sample_xyz,
                        Err(TracerError::Cancelled) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return ControlFlow::Break(Cancelled);
                        }
                        Err(error) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.fail(tile, RenderTileFailure::Tracer(error));
                        }
                    };
                    xyz[0] += sample_xyz[0];
                    xyz[1] += sample_xyz[1];
                    xyz[2] += sample_xyz[2];
                }
                pixels.push(xyz);
            }
        }
        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }

        let mut staging = self
            .staging
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let merge_started = Instant::now();
        let mut source_offset = 0usize;
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                return ControlFlow::Break(Cancelled);
            }
            let destination = (py as usize) * (self.settings.width as usize) + bounds.x as usize;
            let source_end = source_offset + bounds.width as usize;
            staging[destination..destination + bounds.width as usize]
                .copy_from_slice(&pixels[source_offset..source_end]);
            source_offset = source_end;
        }
        atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        ControlFlow::Continue(())
    }
}

impl TileKernel for ParallelRenderKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(RENDER_TILE_KERNEL, self.layout.tile_count())
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

#[derive(Debug, Clone)]
enum CinematicAovTileFailure {
    Aov(CinematicAovError),
    Internal(&'static str),
}

struct ParallelCinematicAovKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    camera: &'assets AnimatedCamera,
    exposure: CameraExposure<'assets>,
    cut_side: CutSide,
    settings: &'run Settings,
    config: CinematicAovConfig,
    palette: &'run CinematicAovPalette,
    albedo_cache: &'run AovAlbedoCache,
    staging: &'run Mutex<CinematicAovFilm>,
    failures: &'run Mutex<Option<(u64, CinematicAovTileFailure)>>,
    layout: RenderTileLayout,
    shutter: ShutterInterval,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl ParallelCinematicAovKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: CinematicAovTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn aov_error(&self, tile: u64, error: CinematicAovError) -> ControlFlow<Cancelled, ()> {
        if matches!(error, CinematicAovError::Tracer(TracerError::Cancelled)) {
            ControlFlow::Break(Cancelled)
        } else {
            self.fail(tile, CinematicAovTileFailure::Aov(error))
        }
    }

    #[allow(clippy::too_many_lines)]
    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("AOV tile outside planned layout"),
            );
        };
        let Some(pixel_count) = bounds
            .width
            .checked_mul(bounds.height)
            .and_then(|count| usize::try_from(count).ok())
        else {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("AOV tile pixel count overflow"),
            );
        };
        let mut pixels =
            match CinematicAovTileAccumulator::try_new(pixel_count, self.config.profile()) {
                Ok(pixels) => pixels,
                Err(error) => return self.aov_error(tile, error),
            };

        let compute_started = Instant::now();
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();
        let capture_primary = self.config.captures_primary();
        let capture_ids = self.config.captures_ids();
        let mut local_pixel = 0usize;
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                return ControlFlow::Break(Cancelled);
            }
            for px in bounds.x..bounds.x + bounds.width {
                let Some(pixel) = py
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(px))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        CinematicAovTileFailure::Internal(
                            "AOV pixel identity overflow after preflight",
                        ),
                    );
                };
                for sample in 0..self.settings.spp {
                    if cx.checkpoint().is_err() {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    let traced = match trace_pixel_sample_with_primary(
                        self.scene,
                        self.lighting,
                        cx,
                        self.settings,
                        kn,
                        self.sobol,
                        key,
                        pixel,
                        sample,
                        Some(self.shutter),
                        self.camera_path,
                        capture_primary,
                    ) {
                        Ok(traced) => traced,
                        Err(TracerError::Cancelled) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return ControlFlow::Break(Cancelled);
                        }
                        Err(error) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.aov_error(tile, CinematicAovError::Tracer(error));
                        }
                    };
                    let aligned = if capture_primary {
                        let Some(split) = traced.contribution_split else {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self
                                .aov_error(tile, CinematicAovError::SampleAlignmentMismatch);
                        };
                        let Some(absolute_time_s) = traced.absolute_time_s else {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.aov_error(
                                tile,
                                CinematicAovError::Tracer(TracerError::MissingRayTime),
                            );
                        };
                        let primary = match prepare_aligned_aov_primary(
                            self.scene,
                            self.camera,
                            self.exposure,
                            self.cut_side,
                            cx,
                            self.settings,
                            self.config.provenance(),
                            self.palette,
                            self.albedo_cache,
                            capture_ids,
                            absolute_time_s,
                            traced.primary.as_ref(),
                        ) {
                            Ok(primary) => primary,
                            Err(error) => {
                                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                                return self.aov_error(tile, error);
                            }
                        };
                        Some(AlignedAovSample {
                            beauty_xyz: traced.xyz,
                            direct_xyz: split.direct_xyz,
                            indirect_xyz: split.indirect_xyz,
                            emission_xyz: split.emission_xyz,
                            pixel_jitter: traced.pixel_jitter,
                            absolute_sample: sample,
                            primary,
                        })
                    } else {
                        None
                    };
                    if let Err(error) = pixels.push(local_pixel, traced.xyz, aligned) {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return self.aov_error(tile, error);
                    }
                }
                local_pixel += 1;
            }
        }
        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
        if local_pixel != pixel_count {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("AOV tile accumulator length mismatch"),
            );
        }
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }

        let mut staging = self
            .staging
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let merge_started = Instant::now();
        let copied = staging.copy_fresh_tile(
            bounds.x,
            bounds.y,
            bounds.width,
            bounds.height,
            &pixels,
            || cx.checkpoint().is_ok(),
        );
        atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        match copied {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => self.aov_error(tile, error),
        }
    }
}

impl TileKernel for ParallelCinematicAovKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(CINEMATIC_AOV_TILE_KERNEL, self.layout.tile_count())
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

struct AdaptiveRenderState {
    xyz: Vec<[f64; 3]>,
    mean_xyz: Vec<[f64; 3]>,
    m2_xyz: Vec<[f64; 3]>,
    sample_counts: Vec<u32>,
    decisions: Vec<AdaptiveDecision>,
}

impl AdaptiveRenderState {
    fn try_new(pixel_count: usize, state_bytes: u64) -> Result<Self, RenderExecutionError> {
        let allocation_error = |site| {
            RenderExecutionError::Allocation(AllocError::OutOfMemory {
                site,
                requested_bytes: usize::try_from(state_bytes).unwrap_or(usize::MAX),
            })
        };
        let mut xyz = Vec::new();
        xyz.try_reserve_exact(pixel_count)
            .map_err(|_| allocation_error("render-adaptive-xyz"))?;
        xyz.resize(pixel_count, [0.0; 3]);
        let mut mean_xyz = Vec::new();
        mean_xyz
            .try_reserve_exact(pixel_count)
            .map_err(|_| allocation_error("render-adaptive-mean"))?;
        mean_xyz.resize(pixel_count, [0.0; 3]);
        let mut m2_xyz = Vec::new();
        m2_xyz
            .try_reserve_exact(pixel_count)
            .map_err(|_| allocation_error("render-adaptive-m2"))?;
        m2_xyz.resize(pixel_count, [0.0; 3]);
        let mut sample_counts = Vec::new();
        sample_counts
            .try_reserve_exact(pixel_count)
            .map_err(|_| allocation_error("render-adaptive-sample-counts"))?;
        sample_counts.resize(pixel_count, 0);
        let mut decisions = Vec::new();
        decisions
            .try_reserve_exact(pixel_count)
            .map_err(|_| allocation_error("render-adaptive-decisions"))?;
        decisions.resize(pixel_count, AdaptiveDecision::MaximumSamples);
        Ok(Self {
            xyz,
            mean_xyz,
            m2_xyz,
            sample_counts,
            decisions,
        })
    }

    fn into_film(
        self,
        settings: &Settings,
        policy: AdaptiveSamplingConfig,
        time_mode: FilmTimeMode,
    ) -> AdaptiveFilm {
        AdaptiveFilm {
            width: settings.width,
            height: settings.height,
            xyz: self.xyz,
            mean_xyz_aov: self.mean_xyz,
            m2_xyz: self.m2_xyz,
            sample_counts: self.sample_counts,
            decisions: self.decisions,
            maximum_samples: settings.spp,
            policy,
            sampler: settings.sampler,
            stream_seed: settings.seed,
            semantics_version: ADAPTIVE_SAMPLING_SEMANTICS_VERSION,
            time_mode,
        }
    }
}

/// Private all-or-nothing state for one adaptive cinematic AOV render.
///
/// A single mutex protects the beauty Welford state and its corresponding AOV
/// planes at tile publication.  Tiles own disjoint pixels, so lock ordering
/// cannot alter arithmetic; the mutex only prevents a partially copied tile
/// from ever becoming a public result after cancellation or failure.
struct AdaptiveCinematicAovStaging {
    beauty: AdaptiveRenderState,
    aov: AdaptiveAovAccumulator,
}

struct ParallelAdaptiveCinematicAovKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    camera: &'assets AnimatedCamera,
    exposure: CameraExposure<'assets>,
    cut_side: CutSide,
    settings: &'run Settings,
    policy: AdaptiveSamplingConfig,
    config: CinematicAovConfig,
    palette: &'run CinematicAovPalette,
    albedo_cache: &'run AovAlbedoCache,
    staging: &'run Mutex<AdaptiveCinematicAovStaging>,
    failures: &'run Mutex<Option<(u64, CinematicAovTileFailure)>>,
    layout: RenderTileLayout,
    shutter: ShutterInterval,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl ParallelAdaptiveCinematicAovKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: CinematicAovTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn aov_error(&self, tile: u64, error: CinematicAovError) -> ControlFlow<Cancelled, ()> {
        if matches!(error, CinematicAovError::Tracer(TracerError::Cancelled)) {
            ControlFlow::Break(Cancelled)
        } else {
            self.fail(tile, CinematicAovTileFailure::Aov(error))
        }
    }

    #[allow(clippy::too_many_lines)]
    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("adaptive AOV tile outside planned layout"),
            );
        };
        let Some(pixel_count) = bounds
            .width
            .checked_mul(bounds.height)
            .and_then(|count| usize::try_from(count).ok())
        else {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("adaptive AOV tile pixel count overflow"),
            );
        };
        let mut aov_pixels =
            match CinematicAovTileAccumulator::try_new(pixel_count, self.config.profile()) {
                Ok(pixels) => pixels,
                Err(error) => return self.aov_error(tile, error),
            };
        let mut adaptive_pixels = Vec::new();
        if adaptive_pixels.try_reserve_exact(pixel_count).is_err() {
            return self.fail(
                tile,
                CinematicAovTileFailure::Aov(CinematicAovError::AllocationRefused),
            );
        }

        let compute_started = Instant::now();
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();
        let capture_primary = self.config.captures_primary();
        let capture_ids = self.config.captures_ids();
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                return ControlFlow::Break(Cancelled);
            }
            for px in bounds.x..bounds.x + bounds.width {
                let Some(pixel) = py
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(px))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        CinematicAovTileFailure::Internal(
                            "adaptive AOV pixel identity overflow after preflight",
                        ),
                    );
                };
                let local_pixel = adaptive_pixels.len();
                let mut accumulator = AdaptivePixelAccumulator::EMPTY;
                for sample in 0..self.settings.spp {
                    if cx.checkpoint().is_err() {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    let traced = match trace_pixel_sample_with_primary(
                        self.scene,
                        self.lighting,
                        cx,
                        self.settings,
                        kn,
                        self.sobol,
                        key,
                        pixel,
                        sample,
                        Some(self.shutter),
                        self.camera_path,
                        capture_primary,
                    ) {
                        Ok(traced) => traced,
                        Err(TracerError::Cancelled) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return ControlFlow::Break(Cancelled);
                        }
                        Err(error) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.aov_error(tile, CinematicAovError::Tracer(error));
                        }
                    };
                    if let Err(error) = accumulator.push(traced.xyz) {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return self.aov_error(tile, CinematicAovError::Adaptive(error));
                    }
                    let aligned = if capture_primary {
                        let Some(split) = traced.contribution_split else {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self
                                .aov_error(tile, CinematicAovError::SampleAlignmentMismatch);
                        };
                        let Some(absolute_time_s) = traced.absolute_time_s else {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.aov_error(
                                tile,
                                CinematicAovError::Tracer(TracerError::MissingRayTime),
                            );
                        };
                        let primary = match prepare_aligned_aov_primary(
                            self.scene,
                            self.camera,
                            self.exposure,
                            self.cut_side,
                            cx,
                            self.settings,
                            self.config.provenance(),
                            self.palette,
                            self.albedo_cache,
                            capture_ids,
                            absolute_time_s,
                            traced.primary.as_ref(),
                        ) {
                            Ok(primary) => primary,
                            Err(error) => {
                                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                                return self.aov_error(tile, error);
                            }
                        };
                        Some(AlignedAovSample {
                            beauty_xyz: traced.xyz,
                            direct_xyz: split.direct_xyz,
                            indirect_xyz: split.indirect_xyz,
                            emission_xyz: split.emission_xyz,
                            pixel_jitter: traced.pixel_jitter,
                            absolute_sample: sample,
                            primary,
                        })
                    } else {
                        None
                    };
                    if let Err(error) = aov_pixels.push(local_pixel, traced.xyz, aligned) {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return self.aov_error(tile, error);
                    }
                    if let Some(decision) = accumulator.decision(self.policy, self.settings.spp) {
                        accumulator.decision = Some(decision);
                        break;
                    }
                }
                if accumulator.decision.is_none() {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        CinematicAovTileFailure::Internal(
                            "adaptive AOV pixel had no terminal decision",
                        ),
                    );
                }
                adaptive_pixels.push(accumulator);
            }
        }
        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
        if adaptive_pixels.len() != pixel_count {
            return self.fail(
                tile,
                CinematicAovTileFailure::Internal("adaptive AOV tile accumulator length mismatch"),
            );
        }
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }

        let mut staging = self
            .staging
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let merge_started = Instant::now();
        let mut source_offset = 0usize;
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                return ControlFlow::Break(Cancelled);
            }
            let destination = py as usize * self.settings.width as usize + bounds.x as usize;
            for column in 0..bounds.width as usize {
                let pixel = adaptive_pixels[source_offset + column];
                let index = destination + column;
                staging.beauty.xyz[index] = pixel.sum_xyz;
                staging.beauty.mean_xyz[index] = pixel.mean_xyz;
                staging.beauty.m2_xyz[index] = pixel.m2_xyz;
                staging.beauty.sample_counts[index] = pixel.samples;
                let Some(decision) = pixel.decision else {
                    atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                    return self.fail(
                        tile,
                        CinematicAovTileFailure::Internal(
                            "adaptive AOV pixel lost terminal decision before merge",
                        ),
                    );
                };
                staging.beauty.decisions[index] = decision;
            }
            source_offset += bounds.width as usize;
        }
        let copied = staging.aov.copy_fresh_tile(
            bounds.x,
            bounds.y,
            bounds.width,
            bounds.height,
            &aov_pixels,
            || cx.checkpoint().is_ok(),
        );
        atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        match copied {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => self.aov_error(tile, error),
        }
    }
}

impl TileKernel for ParallelAdaptiveCinematicAovKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(ADAPTIVE_CINEMATIC_AOV_TILE_KERNEL, self.layout.tile_count())
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

fn adaptive_state_bytes(pixel_count: usize) -> Result<u64, RenderExecutionError> {
    let bytes_per_pixel = size_of::<[f64; 3]>()
        .checked_mul(3)
        .and_then(|bytes| bytes.checked_add(size_of::<u32>()))
        .and_then(|bytes| bytes.checked_add(size_of::<AdaptiveDecision>()))
        .ok_or(RenderExecutionError::Internal(
            "adaptive bytes-per-pixel overflow",
        ))?;
    u64::try_from(pixel_count)
        .ok()
        .and_then(|count| count.checked_mul(bytes_per_pixel as u64))
        .ok_or(RenderExecutionError::Config(
            RenderExecutionConfigError::InvalidImageDimensions,
        ))
}

#[allow(clippy::too_many_arguments)]
fn trace_adaptive_pixel(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    kn: f64,
    sobol: Option<&Sobol>,
    key: [u32; 2],
    pixel: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
) -> Result<AdaptivePixelAccumulator, RenderTileFailure> {
    let mut accumulator = AdaptivePixelAccumulator::EMPTY;
    for sample in 0..settings.spp {
        cx.checkpoint()
            .map_err(|_| RenderTileFailure::Tracer(TracerError::Cancelled))?;
        let xyz = trace_pixel_sample(
            scene,
            lighting,
            cx,
            settings,
            kn,
            sobol,
            key,
            pixel,
            sample,
            shutter,
            camera_path,
        )
        .map_err(RenderTileFailure::Tracer)?;
        accumulator.push(xyz).map_err(RenderTileFailure::Adaptive)?;
        if let Some(decision) = accumulator.decision(policy, settings.spp) {
            accumulator.decision = Some(decision);
            return Ok(accumulator);
        }
    }
    Err(RenderTileFailure::Internal(
        "adaptive pixel reached no terminal decision",
    ))
}

struct AdaptiveRenderKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    settings: &'run Settings,
    policy: AdaptiveSamplingConfig,
    state: &'run Mutex<AdaptiveRenderState>,
    failures: &'run Mutex<Option<(u64, RenderTileFailure)>>,
    layout: RenderTileLayout,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl AdaptiveRenderKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: RenderTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("adaptive tile outside planned layout"),
            );
        };
        let Some(pixel_count) = bounds
            .width
            .checked_mul(bounds.height)
            .and_then(|count| usize::try_from(count).ok())
        else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("adaptive tile pixel count overflow"),
            );
        };
        let mut pixels = Vec::new();
        if pixels.try_reserve_exact(pixel_count).is_err() {
            return self.fail(
                tile,
                RenderTileFailure::Allocation(AllocError::OutOfMemory {
                    site: "render-adaptive-tile-pixels",
                    requested_bytes: pixel_count
                        .saturating_mul(size_of::<AdaptivePixelAccumulator>()),
                }),
            );
        }
        let compute_started = Instant::now();
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                return ControlFlow::Break(Cancelled);
            }
            for px in bounds.x..bounds.x + bounds.width {
                let Some(pixel) = py
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(px))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("adaptive pixel identity overflow"),
                    );
                };
                match trace_adaptive_pixel(
                    self.scene,
                    self.lighting,
                    cx,
                    self.settings,
                    self.policy,
                    kn,
                    self.sobol,
                    key,
                    pixel,
                    self.shutter,
                    self.camera_path,
                ) {
                    Ok(pixel) => pixels.push(pixel),
                    Err(RenderTileFailure::Tracer(TracerError::Cancelled)) => {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    Err(error) => {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return self.fail(tile, error);
                    }
                }
            }
        }
        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }

        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let merge_started = Instant::now();
        let mut source_offset = 0usize;
        for py in bounds.y..bounds.y + bounds.height {
            if cx.checkpoint().is_err() {
                atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                return ControlFlow::Break(Cancelled);
            }
            let destination = py as usize * self.settings.width as usize + bounds.x as usize;
            for column in 0..bounds.width as usize {
                let pixel = pixels[source_offset + column];
                let index = destination + column;
                state.xyz[index] = pixel.sum_xyz;
                state.mean_xyz[index] = pixel.mean_xyz;
                state.m2_xyz[index] = pixel.m2_xyz;
                state.sample_counts[index] = pixel.samples;
                let Some(decision) = pixel.decision else {
                    atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("adaptive pixel had no terminal decision"),
                    );
                };
                state.decisions[index] = decision;
            }
            source_offset += bounds.width as usize;
        }
        atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        ControlFlow::Continue(())
    }
}

impl TileKernel for AdaptiveRenderKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(
            "fs-render/adaptive-spectral-film-tile-v1",
            self.layout.tile_count(),
        )
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

struct PendingAdaptiveRenderKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    settings: &'run Settings,
    policy: AdaptiveSamplingConfig,
    state: &'run Mutex<PendingAdaptiveRenderState>,
    failures: &'run Mutex<Option<(u64, RenderTileFailure)>>,
    layout: RenderTileLayout,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    row_quota: Option<NonZeroU32>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl PendingAdaptiveRenderKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: RenderTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("pending adaptive tile outside planned layout"),
            );
        };
        let initial_row = {
            let state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            state.next_row.get(tile as usize).copied()
        };
        let Some(initial_row) = initial_row else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("pending adaptive row state outside layout"),
            );
        };
        if initial_row >= bounds.height {
            return ControlFlow::Continue(());
        }
        let terminal_row = self.row_quota.map_or(bounds.height, |quota| {
            initial_row.saturating_add(quota.get()).min(bounds.height)
        });
        let row_pixels = bounds.width as usize;
        let mut pixels = Vec::new();
        if pixels.try_reserve_exact(row_pixels).is_err() {
            return self.fail(
                tile,
                RenderTileFailure::Allocation(AllocError::OutOfMemory {
                    site: "render-pending-adaptive-row-pixels",
                    requested_bytes: row_pixels
                        .saturating_mul(size_of::<AdaptivePixelAccumulator>()),
                }),
            );
        }
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();

        loop {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let row_offset = {
                let state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                let Some(row) = state.next_row.get(tile as usize) else {
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pending adaptive row state outside layout"),
                    );
                };
                *row
            };
            if row_offset >= terminal_row {
                return ControlFlow::Continue(());
            }

            let compute_started = Instant::now();
            pixels.clear();
            let y = bounds.y + row_offset;
            for column in 0..row_pixels {
                let x = bounds.x + column as u32;
                let Some(pixel) = y
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(x))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pending adaptive pixel identity overflow"),
                    );
                };
                match trace_adaptive_pixel(
                    self.scene,
                    self.lighting,
                    cx,
                    self.settings,
                    self.policy,
                    kn,
                    self.sobol,
                    key,
                    pixel,
                    self.shutter,
                    self.camera_path,
                ) {
                    Ok(pixel) => pixels.push(pixel),
                    Err(RenderTileFailure::Tracer(TracerError::Cancelled)) => {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    Err(error) => {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return self.fail(tile, error);
                    }
                }
            }
            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }

            let mut state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            let merge_started = Instant::now();
            if state.next_row[tile as usize] != row_offset {
                atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                return self.fail(
                    tile,
                    RenderTileFailure::Internal(
                        "pending adaptive row commit lost exclusive ownership",
                    ),
                );
            }
            let start = y as usize * self.settings.width as usize + bounds.x as usize;
            for (column, pixel) in pixels.iter().copied().enumerate() {
                let Some(decision) = pixel.decision else {
                    atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal(
                            "pending adaptive pixel had no terminal decision",
                        ),
                    );
                };
                let index = start + column;
                state.film.xyz[index] = pixel.sum_xyz;
                state.film.mean_xyz[index] = pixel.mean_xyz;
                state.film.m2_xyz[index] = pixel.m2_xyz;
                state.film.sample_counts[index] = pixel.samples;
                state.film.decisions[index] = decision;
            }
            state.next_row[tile as usize] = row_offset + 1;
            atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        }
    }
}

impl TileKernel for PendingAdaptiveRenderKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(
            PENDING_ADAPTIVE_RENDER_TILE_KERNEL,
            self.layout.tile_count(),
        )
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

struct PendingRenderKernel<'run, 'assets> {
    scene: &'assets Scene,
    lighting: &'run AdmittedLighting<'assets>,
    settings: &'run Settings,
    state: &'run Mutex<PendingRenderState>,
    failures: &'run Mutex<Option<(u64, RenderTileFailure)>>,
    layout: RenderTileLayout,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'assets>,
    sobol: Option<&'run Sobol>,
    row_quota: Option<NonZeroU32>,
    compute_ns: &'run AtomicU64,
    merge_ns: &'run AtomicU64,
}

impl PendingRenderKernel<'_, '_> {
    fn fail(&self, tile: u64, failure: RenderTileFailure) -> ControlFlow<Cancelled, ()> {
        let mut recorded = self
            .failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if recorded
            .as_ref()
            .is_none_or(|(recorded_tile, _)| tile < *recorded_tile)
        {
            *recorded = Some((tile, failure));
        }
        ControlFlow::Break(Cancelled)
    }

    fn run_tile(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let Some(bounds) = self.layout.bounds(tile) else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("pending tile outside planned layout"),
            );
        };
        let initial_row = {
            let state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            state.next_row.get(tile as usize).copied()
        };
        let Some(initial_row) = initial_row else {
            return self.fail(
                tile,
                RenderTileFailure::Internal("pending row state outside planned layout"),
            );
        };
        if initial_row >= bounds.height {
            return ControlFlow::Continue(());
        }
        let terminal_row = self.row_quota.map_or(bounds.height, |quota| {
            initial_row.saturating_add(quota.get()).min(bounds.height)
        });
        let row_pixels = bounds.width as usize;
        let mut pixels = Vec::new();
        if pixels.try_reserve_exact(row_pixels).is_err() {
            return self.fail(
                tile,
                RenderTileFailure::Allocation(AllocError::OutOfMemory {
                    site: "render-pending-row-pixels",
                    requested_bytes: row_pixels.saturating_mul(size_of::<[f64; 3]>()),
                }),
            );
        }
        let key = [
            (self.settings.seed & 0xffff_ffff) as u32,
            (self.settings.seed >> 32) as u32,
        ];
        let kn = 1.0 / y_integral();

        loop {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let row_offset = {
                let state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                let Some(row) = state.next_row.get(tile as usize) else {
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pending row state outside planned layout"),
                    );
                };
                *row
            };
            if row_offset >= terminal_row {
                return ControlFlow::Continue(());
            }

            let compute_started = Instant::now();
            pixels.clear();
            {
                let state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                if state.next_row[tile as usize] != row_offset {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pending tile row advanced concurrently"),
                    );
                }
                let y = bounds.y + row_offset;
                let start = y as usize * self.settings.width as usize + bounds.x as usize;
                pixels.extend_from_slice(&state.xyz[start..start + row_pixels]);
            }
            let y = bounds.y + row_offset;
            for (column, xyz) in pixels.iter_mut().enumerate() {
                let x = bounds.x + column as u32;
                let Some(pixel) = y
                    .checked_mul(self.settings.width)
                    .and_then(|row| row.checked_add(x))
                else {
                    atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                    return self.fail(
                        tile,
                        RenderTileFailure::Internal("pending pixel identity overflow"),
                    );
                };
                for sample in 0..self.settings.spp {
                    if cx.checkpoint().is_err() {
                        atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                        return ControlFlow::Break(Cancelled);
                    }
                    let sample_xyz = match trace_pixel_sample(
                        self.scene,
                        self.lighting,
                        cx,
                        self.settings,
                        kn,
                        self.sobol,
                        key,
                        pixel,
                        sample,
                        self.shutter,
                        self.camera_path,
                    ) {
                        Ok(sample_xyz) => sample_xyz,
                        Err(TracerError::Cancelled) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return ControlFlow::Break(Cancelled);
                        }
                        Err(error) => {
                            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
                            return self.fail(tile, RenderTileFailure::Tracer(error));
                        }
                    };
                    xyz[0] += sample_xyz[0];
                    xyz[1] += sample_xyz[1];
                    xyz[2] += sample_xyz[2];
                }
            }
            atomic_saturating_add(self.compute_ns, elapsed_ns(compute_started));
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }

            let mut state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            let merge_started = Instant::now();
            if state.next_row[tile as usize] != row_offset {
                atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
                return self.fail(
                    tile,
                    RenderTileFailure::Internal("pending row commit lost exclusive ownership"),
                );
            }
            let start = y as usize * self.settings.width as usize + bounds.x as usize;
            state.xyz[start..start + row_pixels].copy_from_slice(&pixels);
            state.next_row[tile as usize] = row_offset + 1;
            atomic_saturating_add(self.merge_ns, elapsed_ns(merge_started));
        }
    }
}

impl TileKernel for PendingRenderKernel<'_, '_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new(RENDER_TILE_KERNEL, self.layout.tile_count())
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        self.run_tile(tile, cx)
    }
}

trait RenderPoolRunner {
    fn run_render<K: TileKernel<Out = ()>>(
        &self,
        cx: &Cx<'_>,
        kernel: &K,
        run: RunId,
        lease: &OperationMemoryLease,
    ) -> (Result<(), RunError>, RunReport);
}

impl RenderPoolRunner for TilePool {
    fn run_render<K: TileKernel<Out = ()>>(
        &self,
        cx: &Cx<'_>,
        kernel: &K,
        run: RunId,
        lease: &OperationMemoryLease,
    ) -> (Result<(), RunError>, RunReport) {
        self.run_declared_leased_with_cx(cx, kernel, run, lease)
    }
}

impl RenderPoolRunner for ParkedTilePool<'_, LocalTaskCaps> {
    fn run_render<K: TileKernel<Out = ()>>(
        &self,
        cx: &Cx<'_>,
        kernel: &K,
        run: RunId,
        lease: &OperationMemoryLease,
    ) -> (Result<(), RunError>, RunReport) {
        self.run_declared_leased_with_cx(cx, kernel, run, lease)
    }
}

fn build_render_pool(
    execution: &RenderExecutionConfig,
    mode: ExecMode,
    scheduler_seed: u64,
) -> TilePool {
    let mut config = PoolConfig::for_host(execution.workers, scheduler_seed);
    config.mode = mode;
    config
        .quantum_weights
        .clone_from(&execution.quantum_weights);
    TilePool::new(config)
}

struct ParallelRenderResult {
    xyz: Option<Vec<[f64; 3]>>,
    requested_mode: FilmTimeMode,
    report: RenderExecutionReport,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn render_range_parallel_impl<R: RenderPoolRunner>(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    base_film: Option<&Film>,
    from: u32,
    to: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<ParallelRenderResult, RenderExecutionError> {
    let setup_started = Instant::now();
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        base_film,
        from,
        to,
        shutter,
        camera_path,
    )?;
    let layout = RenderTileLayout::try_new(
        settings.width,
        settings.height,
        execution.tile_width,
        execution.tile_height,
    )
    .map_err(RenderExecutionError::Config)?;
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);

    if to == from && base_film.is_some() {
        return Ok(ParallelRenderResult {
            xyz: None,
            requested_mode,
            report: empty_parallel_report(
                cx,
                layout,
                execution,
                RENDER_TILE_KERNEL,
                elapsed_ns(setup_started),
                0,
                0,
                lease.receipt(),
            ),
        });
    }

    let pixel_count = checked_pixel_len(settings.width, settings.height)?;
    let film_bytes = u64::try_from(pixel_count)
        .ok()
        .and_then(|count| count.checked_mul(size_of::<[f64; 3]>() as u64))
        .ok_or(RenderExecutionError::Config(
            RenderExecutionConfigError::InvalidImageDimensions,
        ))?;
    let retained_charge = base_film
        .map(|_| lease.reserve("render-retained-film", film_bytes))
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let staging_charge = lease
        .reserve("render-staging-film", film_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let mut staged_xyz = Vec::new();
    staged_xyz.try_reserve_exact(pixel_count).map_err(|_| {
        RenderExecutionError::Allocation(AllocError::OutOfMemory {
            site: "render-staging-film",
            requested_bytes: usize::try_from(film_bytes).unwrap_or(usize::MAX),
        })
    })?;
    staged_xyz.resize(pixel_count, [0.0; 3]);

    if to == from {
        drop(staging_charge);
        drop(retained_charge);
        let memory = lease.receipt();
        return Ok(ParallelRenderResult {
            xyz: Some(staged_xyz),
            requested_mode,
            report: empty_parallel_report(
                cx,
                layout,
                execution,
                RENDER_TILE_KERNEL,
                elapsed_ns(setup_started),
                0,
                film_bytes,
                memory,
            ),
        });
    }

    let max_tile_pixels = u64::from(execution.tile_width.min(settings.width))
        .checked_mul(u64::from(execution.tile_height.min(settings.height)))
        .ok_or(RenderExecutionError::Internal(
            "tile scratch pixel envelope overflow",
        ))?;
    let active_worker_ceiling = u64::try_from(execution.workers)
        .unwrap_or(u64::MAX)
        .min(layout.tile_count());
    let tile_scratch_envelope_bytes = max_tile_pixels
        .checked_mul(size_of::<[f64; 3]>() as u64)
        .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
        .ok_or(RenderExecutionError::Internal(
            "tile scratch byte envelope overflow",
        ))?;
    // Hold one worst-case charge across the run. Tile bodies allocate raw
    // fallible Vec storage inside this already-admitted aggregate envelope,
    // so scheduling overlap cannot turn memory admission into a race.
    let tile_scratch_charge = lease
        .reserve("render-tile-scratch-envelope", tile_scratch_envelope_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let sobol_bytes = if settings.sampler == Sampler::OwenSobol {
        3_u64 * size_of::<[u32; 32]>() as u64
    } else {
        0
    };
    let sobol_charge = (sobol_bytes != 0)
        .then(|| lease.reserve("render-sobol-directions", sobol_bytes))
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let staging = Mutex::new(staged_xyz);
    let failures = Mutex::new(None);
    let compute_ns = AtomicU64::new(0);
    let merge_ns = AtomicU64::new(0);
    let kernel = ParallelRenderKernel {
        scene,
        lighting: &lighting,
        settings,
        base_xyz: base_film.map(|film| film.xyz.as_slice()),
        staging: &staging,
        failures: &failures,
        layout,
        from,
        to,
        shutter,
        camera_path,
        sobol: sobol.as_ref(),
        compute_ns: &compute_ns,
        merge_ns: &merge_ns,
    };
    let setup_ns = elapsed_ns(setup_started);
    let traversal_started = Instant::now();
    let (outcome, executor) = runner.run_render(cx, &kernel, execution.run_id, &lease);
    let traversal_ns = elapsed_ns(traversal_started);
    let failure = failures
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone();
    match outcome {
        Err(RunError::Cancelled { .. }) => {
            if let Some((_tile, failure)) = failure {
                return Err(match failure {
                    RenderTileFailure::Tracer(error) => RenderExecutionError::Tracer(error),
                    RenderTileFailure::Allocation(error) => RenderExecutionError::Allocation(error),
                    RenderTileFailure::Adaptive(error) => RenderExecutionError::Adaptive(error),
                    RenderTileFailure::Internal(detail) => RenderExecutionError::Internal(detail),
                });
            }
            return Err(RenderExecutionError::Tracer(TracerError::Cancelled));
        }
        Err(error) => return Err(RenderExecutionError::Executor(error)),
        Ok(()) => {
            if failure.is_some() {
                return Err(RenderExecutionError::Internal(
                    "tile failure side channel disagreed with executor outcome",
                ));
            }
        }
    }
    cx.checkpoint()?;
    drop(kernel);
    let staged_xyz = staging
        .into_inner()
        .map_err(|_| RenderExecutionError::Internal("successful staging mutex was poisoned"))?;
    let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
    let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
    let active_workers = executor.tiles_by_worker.len();
    let worker_capacity_ns = traversal_ns.saturating_mul(active_workers as u64);
    let idle_worker_ns =
        worker_capacity_ns.saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
    drop(sobol);
    drop(sobol_charge);
    drop(tile_scratch_charge);
    drop(staging_charge);
    drop(retained_charge);
    let memory = lease.receipt();
    Ok(ParallelRenderResult {
        xyz: Some(staged_xyz),
        requested_mode,
        report: RenderExecutionReport {
            layout,
            requested_workers: execution.workers,
            workers: active_workers,
            attempt_index: 1,
            retained_film_bytes: if base_film.is_some() { film_bytes } else { 0 },
            staging_film_bytes: film_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: sobol_bytes,
            progress_state_bytes: 0,
            setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns: 0,
            idle_worker_ns,
            executor,
            memory,
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn render_cinematic_with_aovs_execution_impl<R: RenderPoolRunner>(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    config: CinematicAovConfig,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<CinematicAovExecutionOutput, CinematicAovExecutionError> {
    let setup_started = Instant::now();
    cx.checkpoint().map_err(RenderExecutionError::from)?;
    validate_reference_times(config, shutter)?;
    if config.captures_ids() && settings.spp > MAX_EXACT_F32_INTEGER {
        return Err(CinematicAovError::InexactSampleCount {
            samples: settings.spp,
        }
        .into());
    }
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)
        .map_err(RenderExecutionError::from)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        Some(shutter),
        camera_path,
    )
    .map_err(RenderExecutionError::from)?;
    let capture_ids = config.captures_ids();
    let palette = CinematicAovPalette::try_from_scene(scene, config.limits(), capture_ids, cx)?;
    let continuity_fingerprint = cinematic_input_continuity_fingerprint(scene, camera, cx)?;
    let layout = RenderTileLayout::try_new(
        settings.width,
        settings.height,
        execution.tile_width,
        execution.tile_height,
    )
    .map_err(RenderExecutionError::Config)?;
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
    let (pixel_count, staging_film_bytes) =
        CinematicAovFilm::admitted_retained_bytes(settings.width, settings.height, config)?;
    let staging_charge = lease
        .reserve("render-cinematic-aov-staging-film", staging_film_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let staged = CinematicAovFilm::try_new(settings.width, settings.height, config)?;
    validate_binding(
        &staged,
        *settings,
        shutter,
        exposure.shot_id(),
        cut_side,
        &palette,
        continuity_fingerprint,
    )?;

    if settings.spp == 0 {
        drop(staging_charge);
        let memory = lease.receipt();
        return Ok(CinematicAovExecutionOutput {
            film: staged,
            report: empty_parallel_report(
                cx,
                layout,
                execution,
                CINEMATIC_AOV_TILE_KERNEL,
                elapsed_ns(setup_started),
                0,
                staging_film_bytes,
                memory,
            ),
        });
    }

    let albedo_cache = AovAlbedoCache::try_new(scene.primitives.len(), config.captures_primary())?;
    let max_tile_pixels_u64 = u64::from(execution.tile_width.min(settings.width))
        .checked_mul(u64::from(execution.tile_height.min(settings.height)))
        .ok_or(RenderExecutionError::Internal(
            "cinematic AOV tile scratch pixel envelope overflow",
        ))?;
    let max_tile_pixels = usize::try_from(max_tile_pixels_u64).map_err(|_| {
        RenderExecutionError::Internal("cinematic AOV tile scratch length overflow")
    })?;
    let active_worker_ceiling = u64::try_from(execution.workers)
        .unwrap_or(u64::MAX)
        .min(layout.tile_count());
    let tile_scratch_envelope_bytes =
        CinematicAovTileAccumulator::retained_bytes(max_tile_pixels, config.profile())?
            .checked_mul(active_worker_ceiling)
            .ok_or(RenderExecutionError::Internal(
                "cinematic AOV tile scratch byte envelope overflow",
            ))?;
    let tile_scratch_charge = lease
        .reserve(
            "render-cinematic-aov-tile-scratch-envelope",
            tile_scratch_envelope_bytes,
        )
        .map_err(RenderExecutionError::Memory)?;
    let sobol_bytes = if settings.sampler == Sampler::OwenSobol {
        3_u64 * size_of::<[u32; 32]>() as u64
    } else {
        0
    };
    let sobol_charge = (sobol_bytes != 0)
        .then(|| lease.reserve("render-cinematic-aov-sobol-directions", sobol_bytes))
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let staging = Mutex::new(staged);
    let failures = Mutex::new(None);
    let compute_ns = AtomicU64::new(0);
    let merge_ns = AtomicU64::new(0);
    let kernel = ParallelCinematicAovKernel {
        scene,
        lighting: &lighting,
        camera,
        exposure,
        cut_side,
        settings,
        config,
        palette: &palette,
        albedo_cache: &albedo_cache,
        staging: &staging,
        failures: &failures,
        layout,
        shutter,
        camera_path,
        sobol: sobol.as_ref(),
        compute_ns: &compute_ns,
        merge_ns: &merge_ns,
    };
    let setup_ns = elapsed_ns(setup_started);
    let traversal_started = Instant::now();
    let (outcome, executor) = runner.run_render(cx, &kernel, execution.run_id, &lease);
    let traversal_ns = elapsed_ns(traversal_started);
    let failure = failures
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone();
    match outcome {
        Err(RunError::Cancelled { .. }) => {
            if let Some((_tile, failure)) = failure {
                return Err(match failure {
                    CinematicAovTileFailure::Aov(error) => CinematicAovExecutionError::Aov(error),
                    CinematicAovTileFailure::Internal(detail) => {
                        CinematicAovExecutionError::Execution(RenderExecutionError::Internal(
                            detail,
                        ))
                    }
                });
            }
            return Err(CinematicAovExecutionError::Execution(
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        Err(error) => {
            return Err(CinematicAovExecutionError::Execution(
                RenderExecutionError::Executor(error),
            ));
        }
        Ok(()) => {
            if failure.is_some() {
                return Err(CinematicAovExecutionError::Execution(
                    RenderExecutionError::Internal(
                        "cinematic AOV tile failure disagreed with executor outcome",
                    ),
                ));
            }
        }
    }
    cx.checkpoint().map_err(RenderExecutionError::from)?;
    drop(kernel);
    let mut film = staging.into_inner().map_err(|_| {
        CinematicAovExecutionError::Execution(RenderExecutionError::Internal(
            "successful cinematic AOV staging mutex was poisoned",
        ))
    })?;
    let publication_started = Instant::now();
    film.beauty_mut().spp_done = settings.spp;
    film.beauty_mut().time_mode = requested_mode;
    film.bind(render_binding(
        *settings,
        shutter,
        exposure.shot_id(),
        cut_side,
        palette,
        continuity_fingerprint,
    ));
    let publication_ns = elapsed_ns(publication_started);
    let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
    let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
    let active_workers = executor.tiles_by_worker.len();
    let worker_capacity_ns = traversal_ns.saturating_mul(active_workers as u64);
    let idle_worker_ns =
        worker_capacity_ns.saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
    drop(sobol);
    drop(sobol_charge);
    drop(tile_scratch_charge);
    drop(staging_charge);
    let memory = lease.receipt();
    debug_assert_eq!(film.beauty().xyz.len(), pixel_count);
    Ok(CinematicAovExecutionOutput {
        film,
        report: RenderExecutionReport {
            layout,
            requested_workers: execution.workers,
            workers: active_workers,
            attempt_index: 1,
            retained_film_bytes: 0,
            staging_film_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: sobol_bytes,
            progress_state_bytes: 0,
            setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns,
            idle_worker_ns,
            executor,
            memory,
        },
    })
}

/// Tile-parallel adaptive cinematic render with AOVs co-staged against the
/// exact per-pixel accepted path prefix.  The serial adaptive-AOV API remains
/// the arithmetic oracle; this function changes only execution topology.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn render_cinematic_adaptive_with_aovs_execution_impl<R: RenderPoolRunner>(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: ShutterInterval,
    config: CinematicAovConfig,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<AdaptiveCinematicAovExecutionOutput, CinematicAovExecutionError> {
    let setup_started = Instant::now();
    cx.checkpoint().map_err(RenderExecutionError::from)?;
    policy
        .validate_maximum(settings.spp)
        .map_err(CinematicAovError::Adaptive)?;
    validate_reference_times(config, shutter)?;
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)
        .map_err(RenderExecutionError::from)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        Some(shutter),
        camera_path,
    )
    .map_err(RenderExecutionError::from)?;
    let capture_ids = config.captures_ids();
    let palette = CinematicAovPalette::try_from_scene(scene, config.limits(), capture_ids, cx)?;
    let albedo_cache = AovAlbedoCache::try_new(scene.primitives.len(), config.captures_primary())?;
    let continuity_fingerprint = cinematic_input_continuity_fingerprint(scene, camera, cx)?;
    let layout = RenderTileLayout::try_new(
        settings.width,
        settings.height,
        execution.tile_width,
        execution.tile_height,
    )
    .map_err(RenderExecutionError::Config)?;
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
    let (pixel_count, staging_film_bytes) = AdaptiveAovAccumulator::admitted_retained_bytes(
        settings.width,
        settings.height,
        settings.spp,
        config,
    )?;
    let beauty_state_bytes = adaptive_state_bytes(pixel_count)?;
    let aov_only_bytes = staging_film_bytes.checked_sub(beauty_state_bytes).ok_or(
        RenderExecutionError::Internal("adaptive cinematic AOV staging byte accounting underflow"),
    )?;
    let staging_charge = lease
        .reserve(
            "render-adaptive-cinematic-aov-staging-film",
            staging_film_bytes,
        )
        .map_err(RenderExecutionError::Memory)?;
    let staging = AdaptiveCinematicAovStaging {
        beauty: AdaptiveRenderState::try_new(pixel_count, beauty_state_bytes)?,
        aov: AdaptiveAovAccumulator::try_new(
            settings.width,
            settings.height,
            settings.spp,
            config,
        )?,
    };

    let max_tile_pixels_u64 = u64::from(execution.tile_width.min(settings.width))
        .checked_mul(u64::from(execution.tile_height.min(settings.height)))
        .ok_or(RenderExecutionError::Internal(
            "adaptive cinematic AOV tile scratch pixel envelope overflow",
        ))?;
    let max_tile_pixels = usize::try_from(max_tile_pixels_u64).map_err(|_| {
        RenderExecutionError::Internal("adaptive cinematic AOV tile scratch length overflow")
    })?;
    let active_worker_ceiling = u64::try_from(execution.workers)
        .unwrap_or(u64::MAX)
        .min(layout.tile_count());
    let per_tile_scratch_bytes =
        CinematicAovTileAccumulator::retained_bytes(max_tile_pixels, config.profile())?
            .checked_add(
                u64::try_from(max_tile_pixels)
                    .ok()
                    .and_then(|pixels| {
                        pixels.checked_mul(size_of::<AdaptivePixelAccumulator>() as u64)
                    })
                    .ok_or(RenderExecutionError::Internal(
                        "adaptive cinematic AOV adaptive tile scratch byte overflow",
                    ))?,
            )
            .ok_or(RenderExecutionError::Internal(
                "adaptive cinematic AOV tile scratch byte envelope overflow",
            ))?;
    let tile_scratch_envelope_bytes = per_tile_scratch_bytes
        .checked_mul(active_worker_ceiling)
        .ok_or(RenderExecutionError::Internal(
            "adaptive cinematic AOV tile scratch worker envelope overflow",
        ))?;
    let tile_scratch_charge = lease
        .reserve(
            "render-adaptive-cinematic-aov-tile-scratch-envelope",
            tile_scratch_envelope_bytes,
        )
        .map_err(RenderExecutionError::Memory)?;
    let sobol_bytes = if settings.sampler == Sampler::OwenSobol {
        3_u64 * size_of::<[u32; 32]>() as u64
    } else {
        0
    };
    let sobol_charge = (sobol_bytes != 0)
        .then(|| {
            lease.reserve(
                "render-adaptive-cinematic-aov-sobol-directions",
                sobol_bytes,
            )
        })
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let staging = Mutex::new(staging);
    let failures = Mutex::new(None);
    let compute_ns = AtomicU64::new(0);
    let merge_ns = AtomicU64::new(0);
    let kernel = ParallelAdaptiveCinematicAovKernel {
        scene,
        lighting: &lighting,
        camera,
        exposure,
        cut_side,
        settings,
        policy,
        config,
        palette: &palette,
        albedo_cache: &albedo_cache,
        staging: &staging,
        failures: &failures,
        layout,
        shutter,
        camera_path,
        sobol: sobol.as_ref(),
        compute_ns: &compute_ns,
        merge_ns: &merge_ns,
    };
    let setup_ns = elapsed_ns(setup_started);
    let traversal_started = Instant::now();
    let (outcome, executor) = runner.run_render(cx, &kernel, execution.run_id, &lease);
    let traversal_ns = elapsed_ns(traversal_started);
    let failure = failures
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone();
    match outcome {
        Err(RunError::Cancelled { .. }) => {
            if let Some((_tile, failure)) = failure {
                return Err(match failure {
                    CinematicAovTileFailure::Aov(error) => CinematicAovExecutionError::Aov(error),
                    CinematicAovTileFailure::Internal(detail) => {
                        CinematicAovExecutionError::Execution(RenderExecutionError::Internal(
                            detail,
                        ))
                    }
                });
            }
            return Err(CinematicAovExecutionError::Execution(
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        Err(error) => {
            return Err(CinematicAovExecutionError::Execution(
                RenderExecutionError::Executor(error),
            ));
        }
        Ok(()) => {
            if failure.is_some() {
                return Err(CinematicAovExecutionError::Execution(
                    RenderExecutionError::Internal(
                        "adaptive cinematic AOV tile failure disagreed with executor outcome",
                    ),
                ));
            }
        }
    }
    cx.checkpoint().map_err(RenderExecutionError::from)?;
    drop(kernel);
    let staging = staging.into_inner().map_err(|_| {
        CinematicAovExecutionError::Execution(RenderExecutionError::Internal(
            "successful adaptive cinematic AOV staging mutex was poisoned",
        ))
    })?;
    let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
    let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
    let active_workers = executor.tiles_by_worker.len();
    let idle_worker_ns = traversal_ns
        .saturating_mul(active_workers as u64)
        .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
    drop(sobol);
    drop(sobol_charge);
    drop(tile_scratch_charge);
    // Publication transfers the adaptive beauty state and aligned AOV planes
    // to the returned film, so this operation lease no longer owns them.
    drop(staging_charge);
    let memory = lease.receipt();
    let publication_started = Instant::now();
    let beauty = staging.beauty.into_film(settings, policy, requested_mode);
    let film = staging.aov.publish(
        beauty,
        adaptive_render_binding(
            *settings,
            shutter,
            exposure.shot_id(),
            cut_side,
            palette,
            continuity_fingerprint,
            policy,
        ),
    )?;
    let publication_ns = elapsed_ns(publication_started);
    debug_assert_eq!(film.beauty().xyz_sums().len(), pixel_count);
    debug_assert_eq!(
        film.retained_bytes(),
        aov_only_bytes.saturating_add(beauty_state_bytes),
        "AOV retained-byte contract must continue to include adaptive beauty state"
    );
    Ok(AdaptiveCinematicAovExecutionOutput {
        film,
        report: RenderExecutionReport {
            layout,
            requested_workers: execution.workers,
            workers: active_workers,
            attempt_index: 1,
            retained_film_bytes: 0,
            staging_film_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: sobol_bytes,
            progress_state_bytes: 0,
            setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns,
            idle_worker_ns,
            executor,
            memory,
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn render_adaptive_parallel_impl<R: RenderPoolRunner>(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
    let setup_started = Instant::now();
    policy.validate_maximum(settings.spp)?;
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        shutter,
        camera_path,
    )?;
    let layout = RenderTileLayout::try_new(
        settings.width,
        settings.height,
        execution.tile_width,
        execution.tile_height,
    )
    .map_err(RenderExecutionError::Config)?;
    let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
    let pixel_count = checked_pixel_len(settings.width, settings.height)?;
    let state_bytes = adaptive_state_bytes(pixel_count)?;
    let state_charge = lease
        .reserve("render-adaptive-film", state_bytes)
        .map_err(RenderExecutionError::Memory)?;
    let state = AdaptiveRenderState::try_new(pixel_count, state_bytes)?;

    let max_tile_pixels = u64::from(execution.tile_width.min(settings.width))
        .checked_mul(u64::from(execution.tile_height.min(settings.height)))
        .ok_or(RenderExecutionError::Internal(
            "adaptive tile scratch pixel envelope overflow",
        ))?;
    let active_worker_ceiling = u64::try_from(execution.workers)
        .unwrap_or(u64::MAX)
        .min(layout.tile_count());
    let tile_scratch_envelope_bytes = max_tile_pixels
        .checked_mul(size_of::<AdaptivePixelAccumulator>() as u64)
        .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
        .ok_or(RenderExecutionError::Internal(
            "adaptive tile scratch byte envelope overflow",
        ))?;
    let tile_scratch_charge = lease
        .reserve(
            "render-adaptive-tile-scratch-envelope",
            tile_scratch_envelope_bytes,
        )
        .map_err(RenderExecutionError::Memory)?;
    let sampler_state_bytes = if settings.sampler == Sampler::OwenSobol {
        3_u64 * size_of::<[u32; 32]>() as u64
    } else {
        0
    };
    let sampler_charge = (sampler_state_bytes != 0)
        .then(|| lease.reserve("render-adaptive-sobol-directions", sampler_state_bytes))
        .transpose()
        .map_err(RenderExecutionError::Memory)?;
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let state = Mutex::new(state);
    let failures = Mutex::new(None);
    let compute_ns = AtomicU64::new(0);
    let merge_ns = AtomicU64::new(0);
    let kernel = AdaptiveRenderKernel {
        scene,
        lighting: &lighting,
        settings,
        policy,
        state: &state,
        failures: &failures,
        layout,
        shutter,
        camera_path,
        sobol: sobol.as_ref(),
        compute_ns: &compute_ns,
        merge_ns: &merge_ns,
    };
    let setup_ns = elapsed_ns(setup_started);
    let traversal_started = Instant::now();
    let (outcome, executor) = runner.run_render(cx, &kernel, execution.run_id, &lease);
    let traversal_ns = elapsed_ns(traversal_started);
    let failure = failures
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone();
    if let Some(error) = render_outcome_error(outcome, failure) {
        return Err(error);
    }
    cx.checkpoint()?;
    drop(kernel);
    let state = state
        .into_inner()
        .map_err(|_| RenderExecutionError::Internal("adaptive state mutex was poisoned"))?;
    let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
    let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
    let active_workers = executor.tiles_by_worker.len();
    let idle_worker_ns = traversal_ns
        .saturating_mul(active_workers as u64)
        .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
    drop(sobol);
    drop(sampler_charge);
    drop(tile_scratch_charge);
    // Publication transfers every adaptive AOV allocation to the returned
    // film, so the operation lease no longer owns those bytes.
    drop(state_charge);
    let memory = lease.receipt();
    let publication_started = Instant::now();
    let film = state.into_film(settings, policy, requested_mode);
    let publication_ns = elapsed_ns(publication_started);
    Ok(AdaptiveRenderOutput {
        film,
        report: RenderExecutionReport {
            layout,
            requested_workers: execution.workers,
            workers: active_workers,
            attempt_index: 1,
            retained_film_bytes: 0,
            staging_film_bytes: state_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes,
            progress_state_bytes: 0,
            setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns,
            idle_worker_ns,
            executor,
            memory,
        },
    })
}

fn empty_parallel_report(
    cx: &Cx<'_>,
    layout: RenderTileLayout,
    execution: &RenderExecutionConfig,
    kernel: &'static str,
    setup_ns: u64,
    retained_film_bytes: u64,
    staging_film_bytes: u64,
    memory: LeaseReceipt,
) -> RenderExecutionReport {
    RenderExecutionReport {
        layout,
        requested_workers: execution.workers,
        workers: 0,
        attempt_index: 0,
        retained_film_bytes,
        staging_film_bytes,
        tile_scratch_envelope_bytes: 0,
        sampler_state_bytes: 0,
        progress_state_bytes: 0,
        setup_ns,
        traversal_ns: 0,
        tile_compute_ns: 0,
        tile_merge_ns: 0,
        publication_ns: 0,
        idle_worker_ns: 0,
        executor: RunReport {
            kernel,
            mode: cx.mode().name(),
            declared_run: execution.run_id,
            completed: 0,
            total: 0,
            steals: 0,
            cross_ccd_steals: 0,
            cancel_latencies_ns: Vec::new(),
            tiles_by_worker: Vec::new(),
        },
        memory,
    }
}

fn elapsed_ns(started: Instant) -> u64 {
    u64::try_from(started.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

fn atomic_saturating_add(target: &AtomicU64, value: u64) {
    let _ = target.try_update(Ordering::Relaxed, Ordering::Relaxed, |current| {
        Some(current.saturating_add(value))
    });
}

impl<'assets> PendingRender<'assets> {
    /// Begin a fresh static render whose partial pixels remain opaque until a
    /// complete attempt publishes the final film.
    pub fn begin_static(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        Self::begin_impl(scene, cx, settings, None, CameraPath::Legacy, execution)
    }

    /// Begin a fresh legacy-camera motion render.
    pub fn begin_motion(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        shutter: ShutterInterval,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        Self::begin_impl(
            scene,
            cx,
            settings,
            Some(shutter),
            CameraPath::Legacy,
            execution,
        )
    }

    /// Begin a fresh cinematic-camera render with one admitted exposure.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_cinematic(
        scene: &'assets Scene,
        camera: &'assets AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: Settings,
        shutter: ShutterInterval,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        let exposure = camera
            .admit_shutter(cx, shutter, cut_side)
            .map_err(TracerError::from)?;
        Self::begin_impl(
            scene,
            cx,
            settings,
            Some(shutter),
            CameraPath::Cinematic { camera, exposure },
            execution,
        )
    }

    fn begin_impl(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        shutter: Option<ShutterInterval>,
        camera_path: CameraPath<'assets>,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        let setup_started = Instant::now();
        let (lighting, requested_mode) = preflight_render(
            scene,
            cx,
            &settings,
            None,
            0,
            settings.spp,
            shutter,
            camera_path,
        )?;
        let layout = RenderTileLayout::try_new(
            settings.width,
            settings.height,
            execution.tile_width,
            execution.tile_height,
        )
        .map_err(RenderExecutionError::Config)?;
        let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
        let pixel_count = checked_pixel_len(settings.width, settings.height)?;
        let film_bytes = u64::try_from(pixel_count)
            .ok()
            .and_then(|count| count.checked_mul(size_of::<[f64; 3]>() as u64))
            .ok_or(RenderExecutionError::Config(
                RenderExecutionConfigError::InvalidImageDimensions,
            ))?;
        let film_charge = lease
            .reserve("render-pending-film", film_bytes)
            .map_err(RenderExecutionError::Memory)?;
        let mut xyz = Vec::new();
        xyz.try_reserve_exact(pixel_count).map_err(|_| {
            RenderExecutionError::Allocation(AllocError::OutOfMemory {
                site: "render-pending-film",
                requested_bytes: usize::try_from(film_bytes).unwrap_or(usize::MAX),
            })
        })?;
        xyz.resize(pixel_count, [0.0; 3]);

        let tile_count = usize::try_from(layout.tile_count()).map_err(|_| {
            RenderExecutionError::Internal("pending tile count exceeds address space")
        })?;
        let progress_state_bytes = u64::try_from(tile_count)
            .ok()
            .and_then(|count| count.checked_mul(size_of::<u32>() as u64))
            .ok_or(RenderExecutionError::Internal(
                "pending row-prefix byte count overflow",
            ))?;
        let progress_charge = lease
            .reserve("render-pending-row-prefixes", progress_state_bytes)
            .map_err(RenderExecutionError::Memory)?;
        let mut next_row = Vec::new();
        next_row.try_reserve_exact(tile_count).map_err(|_| {
            RenderExecutionError::Allocation(AllocError::OutOfMemory {
                site: "render-pending-row-prefixes",
                requested_bytes: usize::try_from(progress_state_bytes).unwrap_or(usize::MAX),
            })
        })?;
        next_row.resize(tile_count, 0);
        if settings.spp == 0 {
            for (tile, row) in next_row.iter_mut().enumerate() {
                *row = layout
                    .bounds(tile as u64)
                    .ok_or(RenderExecutionError::Internal(
                        "pending zero-sample tile outside layout",
                    ))?
                    .height;
            }
        }

        let needs_sobol = settings.spp != 0 && settings.sampler == Sampler::OwenSobol;
        let sampler_state_bytes = if needs_sobol {
            3_u64 * size_of::<[u32; 32]>() as u64
        } else {
            0
        };
        let sampler_charge = (sampler_state_bytes != 0)
            .then(|| lease.reserve("render-pending-sobol-directions", sampler_state_bytes))
            .transpose()
            .map_err(RenderExecutionError::Memory)?;
        let sobol = needs_sobol.then(|| Sobol::scrambled(3, settings.seed));
        Ok(Self {
            scene,
            lighting,
            settings,
            shutter,
            camera_path,
            requested_mode,
            execution_mode: cx.mode(),
            execution_budget: cx.budget(),
            execution,
            layout,
            state: Mutex::new(PendingRenderState { xyz, next_row }),
            sobol,
            lease,
            film_charge: Some(film_charge),
            progress_charge: Some(progress_charge),
            sampler_charge,
            film_bytes,
            progress_state_bytes,
            sampler_state_bytes,
            setup_ns: elapsed_ns(setup_started),
            attempts: 0,
        })
    }

    /// Current row/tile progress without exposing partially accumulated pixels.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        let state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        pending_progress(self.layout, &state.next_row, self.attempts)
    }

    /// Resume using a one-shot scoped worker lane. On refusal, the returned
    /// suspension owns this exact job for retry under a fresh `Cx`.
    pub fn resume(self, cx: &Cx<'_>) -> Result<RenderExecutionOutput, RenderSuspension<'assets>> {
        let mut work = self;
        work.start_attempt();
        if cx.mode() != work.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: work.execution_mode,
                    actual: cx.mode(),
                });
            return Err(work.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != work.execution_budget {
            return Err(work.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        let pool = build_render_pool(&work.execution, cx.mode(), work.settings.seed);
        work.resume_with_runner(cx, &pool)
    }

    /// Resume on an already-parked animation/batch crew.
    pub fn resume_on_parked(
        mut self,
        renderer: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, RenderSuspension<'assets>> {
        self.start_attempt();
        if cx.mode() != self.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: self.execution_mode,
                    actual: cx.mode(),
                });
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != self.execution_budget {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        if let Err(cause) = renderer.validate_job(cx, &self.execution) {
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        self.resume_with_runner(cx, renderer.pool)
    }

    /// Advance every tile that was incomplete at attempt start by at most
    /// `rows_per_incomplete_tile` complete rows, then intentionally yield the
    /// still-opaque job at a durable-checkpoint safe point.
    ///
    /// This is a successful bounded attempt, not cancellation. Even when the
    /// quota completes the final rows, the film remains private until a later
    /// ordinary [`Self::resume`] publishes it.
    pub fn advance_to_safe_point(
        self,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<RenderCheckpointYield<'assets>, RenderSuspension<'assets>> {
        let mut work = self;
        work.start_attempt();
        if cx.mode() != work.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: work.execution_mode,
                    actual: cx.mode(),
                });
            return Err(work.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != work.execution_budget {
            return Err(work.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        let pool = build_render_pool(&work.execution, cx.mode(), work.settings.seed);
        work.advance_to_safe_point_with_runner(cx, &pool, rows_per_incomplete_tile)
    }

    /// Bounded safe-point attempt on an already parked animation/batch crew.
    pub fn advance_to_safe_point_on_parked(
        mut self,
        renderer: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<RenderCheckpointYield<'assets>, RenderSuspension<'assets>> {
        self.start_attempt();
        if cx.mode() != self.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: self.execution_mode,
                    actual: cx.mode(),
                });
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != self.execution_budget {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        if let Err(cause) = renderer.validate_job(cx, &self.execution) {
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        self.advance_to_safe_point_with_runner(cx, renderer.pool, rows_per_incomplete_tile)
    }

    fn advance_to_safe_point_with_runner<R: RenderPoolRunner>(
        self,
        cx: &Cx<'_>,
        runner: &R,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<RenderCheckpointYield<'assets>, RenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        if self.progress().completed_tiles == self.layout.tile_count() {
            return self.yield_without_dispatch(cx);
        }

        let row_pixels = u64::from(self.execution.tile_width.min(self.settings.width));
        let active_worker_ceiling = u64::try_from(self.execution.workers)
            .unwrap_or(u64::MAX)
            .min(self.layout.tile_count());
        let tile_scratch_envelope_bytes = row_pixels
            .checked_mul(size_of::<[f64; 3]>() as u64)
            .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
            .ok_or(RenderExecutionError::Internal(
                "pending row-scratch envelope overflow",
            ));
        let tile_scratch_envelope_bytes = match tile_scratch_envelope_bytes {
            Ok(bytes) => bytes,
            Err(cause) => return Err(self.suspend_without_dispatch(cx, cause)),
        };
        let tile_scratch_charge = match self.lease.reserve(
            "render-pending-row-scratch-envelope",
            tile_scratch_envelope_bytes,
        ) {
            Ok(charge) => charge,
            Err(error) => {
                return Err(self.suspend_without_dispatch(cx, RenderExecutionError::Memory(error)));
            }
        };

        let failures = Mutex::new(None);
        let compute_ns = AtomicU64::new(0);
        let merge_ns = AtomicU64::new(0);
        let kernel = PendingRenderKernel {
            scene: self.scene,
            lighting: &self.lighting,
            settings: &self.settings,
            state: &self.state,
            failures: &failures,
            layout: self.layout,
            shutter: self.shutter,
            camera_path: self.camera_path,
            sobol: self.sobol.as_ref(),
            row_quota: Some(rows_per_incomplete_tile),
            compute_ns: &compute_ns,
            merge_ns: &merge_ns,
        };
        let traversal_started = Instant::now();
        let (outcome, executor) =
            runner.run_render(cx, &kernel, self.execution.run_id, &self.lease);
        let traversal_ns = elapsed_ns(traversal_started);
        drop(tile_scratch_charge);
        let failure = failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone();
        let cause = render_outcome_error(outcome, failure);
        let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
        let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
        let active_workers = executor.tiles_by_worker.len();
        let idle_worker_ns = traversal_ns
            .saturating_mul(active_workers as u64)
            .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
        let mut report = RenderExecutionReport {
            layout: self.layout,
            requested_workers: self.execution.workers,
            workers: active_workers,
            attempt_index: self.attempts,
            retained_film_bytes: 0,
            staging_film_bytes: self.film_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: self.sampler_state_bytes,
            progress_state_bytes: self.progress_state_bytes,
            setup_ns: self.setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns: 0,
            idle_worker_ns,
            executor,
            memory: self.lease.receipt(),
        };
        if let Some(cause) = cause {
            report.memory = self.lease.receipt();
            return Err(RenderSuspension {
                work: self,
                cause,
                attempt: report,
            });
        }
        if cx.checkpoint().is_err() {
            report.memory = self.lease.receipt();
            return Err(RenderSuspension {
                work: self,
                cause: RenderExecutionError::Tracer(TracerError::Cancelled),
                attempt: report,
            });
        }
        Ok(RenderCheckpointYield {
            work: self,
            attempt: report,
        })
    }

    fn resume_with_runner<R: RenderPoolRunner>(
        self,
        cx: &Cx<'_>,
        runner: &R,
    ) -> Result<RenderExecutionOutput, RenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        if self.progress().completed_tiles == self.layout.tile_count() {
            return self.publish_without_dispatch(cx);
        }

        let row_pixels = u64::from(self.execution.tile_width.min(self.settings.width));
        let active_worker_ceiling = u64::try_from(self.execution.workers)
            .unwrap_or(u64::MAX)
            .min(self.layout.tile_count());
        let tile_scratch_envelope_bytes = row_pixels
            .checked_mul(size_of::<[f64; 3]>() as u64)
            .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
            .ok_or(RenderExecutionError::Internal(
                "pending row-scratch envelope overflow",
            ));
        let tile_scratch_envelope_bytes = match tile_scratch_envelope_bytes {
            Ok(bytes) => bytes,
            Err(cause) => return Err(self.suspend_without_dispatch(cx, cause)),
        };
        let tile_scratch_charge = match self.lease.reserve(
            "render-pending-row-scratch-envelope",
            tile_scratch_envelope_bytes,
        ) {
            Ok(charge) => charge,
            Err(error) => {
                return Err(self.suspend_without_dispatch(cx, RenderExecutionError::Memory(error)));
            }
        };

        let failures = Mutex::new(None);
        let compute_ns = AtomicU64::new(0);
        let merge_ns = AtomicU64::new(0);
        let kernel = PendingRenderKernel {
            scene: self.scene,
            lighting: &self.lighting,
            settings: &self.settings,
            state: &self.state,
            failures: &failures,
            layout: self.layout,
            shutter: self.shutter,
            camera_path: self.camera_path,
            sobol: self.sobol.as_ref(),
            row_quota: None,
            compute_ns: &compute_ns,
            merge_ns: &merge_ns,
        };
        let traversal_started = Instant::now();
        let (outcome, executor) =
            runner.run_render(cx, &kernel, self.execution.run_id, &self.lease);
        let traversal_ns = elapsed_ns(traversal_started);
        drop(tile_scratch_charge);
        let failure = failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone();
        let cause = render_outcome_error(outcome, failure);
        let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
        let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
        let active_workers = executor.tiles_by_worker.len();
        let idle_worker_ns = traversal_ns
            .saturating_mul(active_workers as u64)
            .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
        let mut report = RenderExecutionReport {
            layout: self.layout,
            requested_workers: self.execution.workers,
            workers: active_workers,
            attempt_index: self.attempts,
            retained_film_bytes: 0,
            staging_film_bytes: self.film_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: self.sampler_state_bytes,
            progress_state_bytes: self.progress_state_bytes,
            setup_ns: self.setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns: 0,
            idle_worker_ns,
            executor,
            memory: self.lease.receipt(),
        };
        if let Some(cause) = cause {
            report.memory = self.lease.receipt();
            return Err(RenderSuspension {
                work: self,
                cause,
                attempt: report,
            });
        }
        if self.progress().completed_tiles != self.layout.tile_count() {
            return Err(RenderSuspension {
                work: self,
                cause: RenderExecutionError::Internal(
                    "executor succeeded before every pending row committed",
                ),
                attempt: report,
            });
        }
        if cx.checkpoint().is_err() {
            report.memory = self.lease.receipt();
            return Err(RenderSuspension {
                work: self,
                cause: RenderExecutionError::Tracer(TracerError::Cancelled),
                attempt: report,
            });
        }
        Ok(self.publish(report))
    }

    fn start_attempt(&mut self) {
        self.attempts = self.attempts.saturating_add(1);
    }

    fn suspend_without_dispatch(
        self,
        cx: &Cx<'_>,
        cause: RenderExecutionError,
    ) -> RenderSuspension<'assets> {
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.film_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        RenderSuspension {
            work: self,
            cause,
            attempt: report,
        }
    }

    fn yield_without_dispatch(
        self,
        cx: &Cx<'_>,
    ) -> Result<RenderCheckpointYield<'assets>, RenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.film_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        Ok(RenderCheckpointYield {
            work: self,
            attempt: report,
        })
    }

    fn publish_without_dispatch(
        self,
        cx: &Cx<'_>,
    ) -> Result<RenderExecutionOutput, RenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.film_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        Ok(self.publish(report))
    }

    fn publish(mut self, mut report: RenderExecutionReport) -> RenderExecutionOutput {
        let publication_started = Instant::now();
        let state = self
            .state
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let PendingRenderState { xyz, next_row } = state;
        drop(next_row);
        drop(self.progress_charge.take());
        drop(self.sobol);
        drop(self.sampler_charge.take());
        // Publication transfers the film allocation to the returned value, so
        // it leaves the operation lease even though `xyz` remains live.
        drop(self.film_charge.take());
        report.memory = self.lease.receipt();
        let film = Film {
            width: self.settings.width,
            height: self.settings.height,
            xyz,
            spp_done: self.settings.spp,
            time_mode: if self.settings.spp == 0 {
                FilmTimeMode::Uninitialized
            } else {
                self.requested_mode
            },
        };
        report.publication_ns = elapsed_ns(publication_started);
        RenderExecutionOutput { film, report }
    }
}

impl<'assets> PendingAdaptiveRender<'assets> {
    /// Begin a fresh static adaptive render.
    pub fn begin_static(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        policy: AdaptiveSamplingConfig,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        Self::begin_impl(
            scene,
            cx,
            settings,
            policy,
            None,
            CameraPath::Legacy,
            execution,
        )
    }

    /// Begin a fresh legacy-camera motion adaptive render.
    pub fn begin_motion(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        policy: AdaptiveSamplingConfig,
        shutter: ShutterInterval,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        Self::begin_impl(
            scene,
            cx,
            settings,
            policy,
            Some(shutter),
            CameraPath::Legacy,
            execution,
        )
    }

    /// Begin a fresh cinematic-camera adaptive render.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_cinematic(
        scene: &'assets Scene,
        camera: &'assets AnimatedCamera,
        cut_side: CutSide,
        cx: &Cx<'_>,
        settings: Settings,
        policy: AdaptiveSamplingConfig,
        shutter: ShutterInterval,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        let exposure = camera
            .admit_shutter(cx, shutter, cut_side)
            .map_err(TracerError::from)?;
        Self::begin_impl(
            scene,
            cx,
            settings,
            policy,
            Some(shutter),
            CameraPath::Cinematic { camera, exposure },
            execution,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn begin_impl(
        scene: &'assets Scene,
        cx: &Cx<'_>,
        settings: Settings,
        policy: AdaptiveSamplingConfig,
        shutter: Option<ShutterInterval>,
        camera_path: CameraPath<'assets>,
        execution: RenderExecutionConfig,
    ) -> Result<Self, RenderExecutionError> {
        let setup_started = Instant::now();
        policy.validate_maximum(settings.spp)?;
        let (lighting, requested_mode) = preflight_render(
            scene,
            cx,
            &settings,
            None,
            0,
            settings.spp,
            shutter,
            camera_path,
        )?;
        let layout = RenderTileLayout::try_new(
            settings.width,
            settings.height,
            execution.tile_width,
            execution.tile_height,
        )
        .map_err(RenderExecutionError::Config)?;
        let lease = OperationMemoryLease::bounded(execution.memory_limit_bytes);
        let pixel_count = checked_pixel_len(settings.width, settings.height)?;
        let state_bytes = adaptive_state_bytes(pixel_count)?;
        let state_charge = lease
            .reserve("render-pending-adaptive-film", state_bytes)
            .map_err(RenderExecutionError::Memory)?;
        let film = AdaptiveRenderState::try_new(pixel_count, state_bytes)?;

        let tile_count = usize::try_from(layout.tile_count()).map_err(|_| {
            RenderExecutionError::Internal("pending adaptive tile count exceeds address space")
        })?;
        let progress_state_bytes = u64::try_from(tile_count)
            .ok()
            .and_then(|count| count.checked_mul(size_of::<u32>() as u64))
            .ok_or(RenderExecutionError::Internal(
                "pending adaptive row-prefix byte count overflow",
            ))?;
        let progress_charge = lease
            .reserve("render-pending-adaptive-row-prefixes", progress_state_bytes)
            .map_err(RenderExecutionError::Memory)?;
        let mut next_row = Vec::new();
        next_row.try_reserve_exact(tile_count).map_err(|_| {
            RenderExecutionError::Allocation(AllocError::OutOfMemory {
                site: "render-pending-adaptive-row-prefixes",
                requested_bytes: usize::try_from(progress_state_bytes).unwrap_or(usize::MAX),
            })
        })?;
        next_row.resize(tile_count, 0);

        let sampler_state_bytes = if settings.sampler == Sampler::OwenSobol {
            3_u64 * size_of::<[u32; 32]>() as u64
        } else {
            0
        };
        let sampler_charge = (sampler_state_bytes != 0)
            .then(|| {
                lease.reserve(
                    "render-pending-adaptive-sobol-directions",
                    sampler_state_bytes,
                )
            })
            .transpose()
            .map_err(RenderExecutionError::Memory)?;
        let sobol =
            (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
        Ok(Self {
            scene,
            lighting,
            settings,
            policy,
            shutter,
            camera_path,
            requested_mode,
            execution_mode: cx.mode(),
            execution_budget: cx.budget(),
            execution,
            layout,
            state: Mutex::new(PendingAdaptiveRenderState { film, next_row }),
            sobol,
            lease,
            state_charge: Some(state_charge),
            progress_charge: Some(progress_charge),
            sampler_charge,
            state_bytes,
            progress_state_bytes,
            sampler_state_bytes,
            setup_ns: elapsed_ns(setup_started),
            attempts: 0,
        })
    }

    /// Current row/tile progress without exposing partial film state.
    #[must_use]
    pub fn progress(&self) -> RenderProgress {
        let state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        pending_progress(self.layout, &state.next_row, self.attempts)
    }

    /// Resume on a one-shot scoped worker lane.
    pub fn resume(
        self,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, AdaptiveRenderSuspension<'assets>> {
        let mut work = self;
        work.start_attempt();
        if cx.mode() != work.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: work.execution_mode,
                    actual: cx.mode(),
                });
            return Err(work.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != work.execution_budget {
            return Err(work.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        let pool = build_render_pool(&work.execution, cx.mode(), work.settings.seed);
        work.resume_with_runner(cx, &pool)
    }

    /// Resume on an already-parked animation/batch crew.
    pub fn resume_on_parked(
        mut self,
        renderer: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, AdaptiveRenderSuspension<'assets>> {
        self.start_attempt();
        if cx.mode() != self.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: self.execution_mode,
                    actual: cx.mode(),
                });
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != self.execution_budget {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        if let Err(cause) = renderer.validate_job(cx, &self.execution) {
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        self.resume_with_runner(cx, renderer.pool)
    }

    /// Advance every adaptive tile that was incomplete at attempt start by at
    /// most `rows_per_incomplete_tile` complete rows, then intentionally yield
    /// the opaque sums, moments, counts, and decisions at a safe point.
    pub fn advance_to_safe_point(
        self,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<AdaptiveRenderCheckpointYield<'assets>, AdaptiveRenderSuspension<'assets>> {
        let mut work = self;
        work.start_attempt();
        if cx.mode() != work.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: work.execution_mode,
                    actual: cx.mode(),
                });
            return Err(work.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != work.execution_budget {
            return Err(work.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        let pool = build_render_pool(&work.execution, cx.mode(), work.settings.seed);
        work.advance_to_safe_point_with_runner(cx, &pool, rows_per_incomplete_tile)
    }

    /// Bounded adaptive safe-point attempt on a parked animation/batch crew.
    pub fn advance_to_safe_point_on_parked(
        mut self,
        renderer: &ParkedRenderScope<'_>,
        cx: &Cx<'_>,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<AdaptiveRenderCheckpointYield<'assets>, AdaptiveRenderSuspension<'assets>> {
        self.start_attempt();
        if cx.mode() != self.execution_mode {
            let cause =
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeModeMismatch {
                    expected: self.execution_mode,
                    actual: cx.mode(),
                });
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        if cx.budget() != self.execution_budget {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Config(RenderExecutionConfigError::ResumeBudgetMismatch),
            ));
        }
        if let Err(cause) = renderer.validate_job(cx, &self.execution) {
            return Err(self.suspend_without_dispatch(cx, cause));
        }
        self.advance_to_safe_point_with_runner(cx, renderer.pool, rows_per_incomplete_tile)
    }

    fn advance_to_safe_point_with_runner<R: RenderPoolRunner>(
        self,
        cx: &Cx<'_>,
        runner: &R,
        rows_per_incomplete_tile: NonZeroU32,
    ) -> Result<AdaptiveRenderCheckpointYield<'assets>, AdaptiveRenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        if self.progress().completed_tiles == self.layout.tile_count() {
            return self.yield_without_dispatch(cx);
        }

        let row_pixels = u64::from(self.execution.tile_width.min(self.settings.width));
        let active_worker_ceiling = u64::try_from(self.execution.workers)
            .unwrap_or(u64::MAX)
            .min(self.layout.tile_count());
        let tile_scratch_envelope_bytes = row_pixels
            .checked_mul(size_of::<AdaptivePixelAccumulator>() as u64)
            .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
            .ok_or(RenderExecutionError::Internal(
                "pending adaptive row-scratch envelope overflow",
            ));
        let tile_scratch_envelope_bytes = match tile_scratch_envelope_bytes {
            Ok(bytes) => bytes,
            Err(cause) => return Err(self.suspend_without_dispatch(cx, cause)),
        };
        let tile_scratch_charge = match self.lease.reserve(
            "render-pending-adaptive-row-scratch-envelope",
            tile_scratch_envelope_bytes,
        ) {
            Ok(charge) => charge,
            Err(error) => {
                return Err(self.suspend_without_dispatch(cx, RenderExecutionError::Memory(error)));
            }
        };

        let failures = Mutex::new(None);
        let compute_ns = AtomicU64::new(0);
        let merge_ns = AtomicU64::new(0);
        let kernel = PendingAdaptiveRenderKernel {
            scene: self.scene,
            lighting: &self.lighting,
            settings: &self.settings,
            policy: self.policy,
            state: &self.state,
            failures: &failures,
            layout: self.layout,
            shutter: self.shutter,
            camera_path: self.camera_path,
            sobol: self.sobol.as_ref(),
            row_quota: Some(rows_per_incomplete_tile),
            compute_ns: &compute_ns,
            merge_ns: &merge_ns,
        };
        let traversal_started = Instant::now();
        let (outcome, executor) =
            runner.run_render(cx, &kernel, self.execution.run_id, &self.lease);
        let traversal_ns = elapsed_ns(traversal_started);
        drop(tile_scratch_charge);
        let failure = failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone();
        let cause = render_outcome_error(outcome, failure);
        let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
        let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
        let active_workers = executor.tiles_by_worker.len();
        let idle_worker_ns = traversal_ns
            .saturating_mul(active_workers as u64)
            .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
        let mut report = RenderExecutionReport {
            layout: self.layout,
            requested_workers: self.execution.workers,
            workers: active_workers,
            attempt_index: self.attempts,
            retained_film_bytes: 0,
            staging_film_bytes: self.state_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: self.sampler_state_bytes,
            progress_state_bytes: self.progress_state_bytes,
            setup_ns: self.setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns: 0,
            idle_worker_ns,
            executor,
            memory: self.lease.receipt(),
        };
        if let Some(cause) = cause {
            report.memory = self.lease.receipt();
            return Err(AdaptiveRenderSuspension {
                work: self,
                cause,
                attempt: report,
            });
        }
        if cx.checkpoint().is_err() {
            report.memory = self.lease.receipt();
            return Err(AdaptiveRenderSuspension {
                work: self,
                cause: RenderExecutionError::Tracer(TracerError::Cancelled),
                attempt: report,
            });
        }
        Ok(AdaptiveRenderCheckpointYield {
            work: self,
            attempt: report,
        })
    }

    fn resume_with_runner<R: RenderPoolRunner>(
        self,
        cx: &Cx<'_>,
        runner: &R,
    ) -> Result<AdaptiveRenderOutput, AdaptiveRenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        if self.progress().completed_tiles == self.layout.tile_count() {
            return self.publish_without_dispatch(cx);
        }

        let row_pixels = u64::from(self.execution.tile_width.min(self.settings.width));
        let active_worker_ceiling = u64::try_from(self.execution.workers)
            .unwrap_or(u64::MAX)
            .min(self.layout.tile_count());
        let tile_scratch_envelope_bytes = row_pixels
            .checked_mul(size_of::<AdaptivePixelAccumulator>() as u64)
            .and_then(|bytes| bytes.checked_mul(active_worker_ceiling))
            .ok_or(RenderExecutionError::Internal(
                "pending adaptive row-scratch envelope overflow",
            ));
        let tile_scratch_envelope_bytes = match tile_scratch_envelope_bytes {
            Ok(bytes) => bytes,
            Err(cause) => return Err(self.suspend_without_dispatch(cx, cause)),
        };
        let tile_scratch_charge = match self.lease.reserve(
            "render-pending-adaptive-row-scratch-envelope",
            tile_scratch_envelope_bytes,
        ) {
            Ok(charge) => charge,
            Err(error) => {
                return Err(self.suspend_without_dispatch(cx, RenderExecutionError::Memory(error)));
            }
        };

        let failures = Mutex::new(None);
        let compute_ns = AtomicU64::new(0);
        let merge_ns = AtomicU64::new(0);
        let kernel = PendingAdaptiveRenderKernel {
            scene: self.scene,
            lighting: &self.lighting,
            settings: &self.settings,
            policy: self.policy,
            state: &self.state,
            failures: &failures,
            layout: self.layout,
            shutter: self.shutter,
            camera_path: self.camera_path,
            sobol: self.sobol.as_ref(),
            row_quota: None,
            compute_ns: &compute_ns,
            merge_ns: &merge_ns,
        };
        let traversal_started = Instant::now();
        let (outcome, executor) =
            runner.run_render(cx, &kernel, self.execution.run_id, &self.lease);
        let traversal_ns = elapsed_ns(traversal_started);
        drop(tile_scratch_charge);
        let failure = failures
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone();
        let cause = render_outcome_error(outcome, failure);
        let tile_compute_ns = compute_ns.load(Ordering::Relaxed);
        let tile_merge_ns = merge_ns.load(Ordering::Relaxed);
        let active_workers = executor.tiles_by_worker.len();
        let idle_worker_ns = traversal_ns
            .saturating_mul(active_workers as u64)
            .saturating_sub(tile_compute_ns.saturating_add(tile_merge_ns));
        let mut report = RenderExecutionReport {
            layout: self.layout,
            requested_workers: self.execution.workers,
            workers: active_workers,
            attempt_index: self.attempts,
            retained_film_bytes: 0,
            staging_film_bytes: self.state_bytes,
            tile_scratch_envelope_bytes,
            sampler_state_bytes: self.sampler_state_bytes,
            progress_state_bytes: self.progress_state_bytes,
            setup_ns: self.setup_ns,
            traversal_ns,
            tile_compute_ns,
            tile_merge_ns,
            publication_ns: 0,
            idle_worker_ns,
            executor,
            memory: self.lease.receipt(),
        };
        if let Some(cause) = cause {
            report.memory = self.lease.receipt();
            return Err(AdaptiveRenderSuspension {
                work: self,
                cause,
                attempt: report,
            });
        }
        if self.progress().completed_tiles != self.layout.tile_count() {
            return Err(AdaptiveRenderSuspension {
                work: self,
                cause: RenderExecutionError::Internal(
                    "executor succeeded before every pending adaptive row committed",
                ),
                attempt: report,
            });
        }
        if cx.checkpoint().is_err() {
            report.memory = self.lease.receipt();
            return Err(AdaptiveRenderSuspension {
                work: self,
                cause: RenderExecutionError::Tracer(TracerError::Cancelled),
                attempt: report,
            });
        }
        Ok(self.publish(report))
    }

    fn start_attempt(&mut self) {
        self.attempts = self.attempts.saturating_add(1);
    }

    fn suspend_without_dispatch(
        self,
        cx: &Cx<'_>,
        cause: RenderExecutionError,
    ) -> AdaptiveRenderSuspension<'assets> {
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            PENDING_ADAPTIVE_RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.state_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        AdaptiveRenderSuspension {
            work: self,
            cause,
            attempt: report,
        }
    }

    fn yield_without_dispatch(
        self,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderCheckpointYield<'assets>, AdaptiveRenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            PENDING_ADAPTIVE_RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.state_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        Ok(AdaptiveRenderCheckpointYield {
            work: self,
            attempt: report,
        })
    }

    fn publish_without_dispatch(
        self,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveRenderOutput, AdaptiveRenderSuspension<'assets>> {
        if cx.checkpoint().is_err() {
            return Err(self.suspend_without_dispatch(
                cx,
                RenderExecutionError::Tracer(TracerError::Cancelled),
            ));
        }
        let mut report = empty_parallel_report(
            cx,
            self.layout,
            &self.execution,
            PENDING_ADAPTIVE_RENDER_TILE_KERNEL,
            self.setup_ns,
            0,
            self.state_bytes,
            self.lease.receipt(),
        );
        report.attempt_index = self.attempts;
        report.sampler_state_bytes = self.sampler_state_bytes;
        report.progress_state_bytes = self.progress_state_bytes;
        Ok(self.publish(report))
    }

    fn publish(mut self, mut report: RenderExecutionReport) -> AdaptiveRenderOutput {
        let publication_started = Instant::now();
        let state = self
            .state
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let PendingAdaptiveRenderState { film, next_row } = state;
        drop(next_row);
        drop(self.progress_charge.take());
        drop(self.sobol);
        drop(self.sampler_charge.take());
        // Publication transfers all adaptive AOV allocations to the returned
        // film, so those bytes leave the operation lease while remaining live.
        drop(self.state_charge.take());
        report.memory = self.lease.receipt();
        let film = film.into_film(&self.settings, self.policy, self.requested_mode);
        report.publication_ns = elapsed_ns(publication_started);
        AdaptiveRenderOutput { film, report }
    }
}

fn pending_progress(layout: RenderTileLayout, next_row: &[u32], attempts: u64) -> RenderProgress {
    let mut committed_tile_rows = 0_u64;
    let mut total_tile_rows = 0_u64;
    let mut completed_tiles = 0_u64;
    for tile in 0..layout.tile_count() {
        let bounds = layout
            .bounds(tile)
            .expect("validated pending tile remains inside its layout");
        let committed = next_row[tile as usize].min(bounds.height);
        committed_tile_rows += u64::from(committed);
        total_tile_rows += u64::from(bounds.height);
        completed_tiles += u64::from(committed == bounds.height);
    }
    RenderProgress {
        committed_tile_rows,
        total_tile_rows,
        completed_tiles,
        total_tiles: layout.tile_count(),
        attempts,
    }
}

fn render_outcome_error(
    outcome: Result<(), RunError>,
    failure: Option<(u64, RenderTileFailure)>,
) -> Option<RenderExecutionError> {
    match outcome {
        Err(RunError::Cancelled { .. }) => Some(match failure {
            Some((_tile, RenderTileFailure::Tracer(error))) => RenderExecutionError::Tracer(error),
            Some((_tile, RenderTileFailure::Allocation(error))) => {
                RenderExecutionError::Allocation(error)
            }
            Some((_tile, RenderTileFailure::Adaptive(error))) => {
                RenderExecutionError::Adaptive(error)
            }
            Some((_tile, RenderTileFailure::Internal(detail))) => {
                RenderExecutionError::Internal(detail)
            }
            None => RenderExecutionError::Tracer(TracerError::Cancelled),
        }),
        Err(error) => Some(RenderExecutionError::Executor(error)),
        Ok(()) if failure.is_some() => Some(RenderExecutionError::Internal(
            "tile failure side channel disagreed with executor outcome",
        )),
        Ok(()) => None,
    }
}

/// Tile-parallel progressive static render under an explicit worker, tile,
/// memory, and run-identity policy. Publication is all-or-nothing.
pub fn render_range_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_range_with_execution_impl(
        scene,
        cx,
        settings,
        film,
        from,
        to,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Tile-parallel progressive motion render under an explicit execution policy.
#[allow(clippy::too_many_arguments)]
pub fn render_motion_range_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_range_with_execution_impl(
        scene,
        cx,
        settings,
        film,
        from,
        to,
        Some(shutter),
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Tile-parallel progressive cinematic-camera render under an explicit
/// execution policy.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_range_with_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    render_range_with_execution_impl(
        scene,
        cx,
        settings,
        film,
        from,
        to,
        Some(shutter),
        CameraPath::Cinematic { camera, exposure },
        execution,
        &pool,
    )
}

#[allow(clippy::too_many_arguments)]
fn render_range_with_execution_impl<R: RenderPoolRunner>(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    from: u32,
    to: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    let result = render_range_parallel_impl(
        scene,
        cx,
        settings,
        Some(film),
        from,
        to,
        shutter,
        camera_path,
        execution,
        runner,
    )?;
    let mut report = result.report;
    if let Some(xyz) = result.xyz {
        let publication_started = Instant::now();
        film.xyz = xyz;
        film.spp_done = to;
        if film.time_mode == FilmTimeMode::Uninitialized {
            film.time_mode = result.requested_mode;
        }
        report.publication_ns = elapsed_ns(publication_started);
    }
    Ok(report)
}

/// Render a fresh static film with one image-sized allocation plus bounded
/// per-tile scratch. Unlike progressive append, no second retained film is
/// required for transactional publication.
pub fn render_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_fresh_with_execution_impl(
        scene,
        cx,
        settings,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Render a fresh motion-blurred film with explicit tile execution.
pub fn render_motion_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_fresh_with_execution_impl(
        scene,
        cx,
        settings,
        Some(shutter),
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Render a fresh cinematic-camera film with explicit tile execution.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_with_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    render_fresh_with_execution_impl(
        scene,
        cx,
        settings,
        Some(shutter),
        CameraPath::Cinematic { camera, exposure },
        execution,
        &pool,
    )
}

/// Render a fresh cinematic-camera film with exactly aligned AOVs under an
/// explicit deterministic tile policy. The existing serial AOV API remains the
/// progressive oracle; this entry point owns one private staging film and
/// publishes it only after complete executor drain and final cancellation
/// admission.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_with_aovs_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    config: CinematicAovConfig,
    execution: &RenderExecutionConfig,
) -> Result<CinematicAovExecutionOutput, CinematicAovExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_cinematic_with_aovs_execution_impl(
        scene, camera, cut_side, cx, settings, shutter, config, execution, &pool,
    )
}

/// Render a fresh adaptive cinematic-camera film with exactly aligned AOVs
/// under an explicit deterministic tile policy.
///
/// This is the parallel counterpart of
/// [`render_cinematic_adaptive_with_aovs`].  It retains every pixel's exact
/// accepted `0..terminal_count` prefix, then publishes raw beauty and AOVs
/// only after all tiles have drained and the final cancellation checkpoint has
/// succeeded.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_adaptive_with_aovs_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: ShutterInterval,
    config: CinematicAovConfig,
    execution: &RenderExecutionConfig,
) -> Result<AdaptiveCinematicAovExecutionOutput, CinematicAovExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_cinematic_adaptive_with_aovs_execution_impl(
        scene, camera, cut_side, cx, settings, policy, shutter, config, execution, &pool,
    )
}

/// Render a fresh static adaptive film with explicit deterministic tile
/// execution. The existing uniform APIs remain the final-quality fallback and
/// bitwise oracle.
pub fn render_adaptive_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    execution: &RenderExecutionConfig,
) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_adaptive_parallel_impl(
        scene,
        cx,
        settings,
        policy,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Render a fresh motion-blurred adaptive film with explicit tile execution.
#[allow(clippy::too_many_arguments)]
pub fn render_motion_adaptive_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    render_adaptive_parallel_impl(
        scene,
        cx,
        settings,
        policy,
        Some(shutter),
        CameraPath::Legacy,
        execution,
        &pool,
    )
}

/// Render a fresh cinematic-camera adaptive film with explicit tile
/// execution.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_adaptive_with_execution(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: ShutterInterval,
    execution: &RenderExecutionConfig,
) -> Result<AdaptiveRenderOutput, RenderExecutionError> {
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    render_adaptive_parallel_impl(
        scene,
        cx,
        settings,
        policy,
        Some(shutter),
        CameraPath::Cinematic { camera, exposure },
        execution,
        &pool,
    )
}

fn render_fresh_with_execution_impl<R: RenderPoolRunner>(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
    execution: &RenderExecutionConfig,
    runner: &R,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    let result = render_range_parallel_impl(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        shutter,
        camera_path,
        execution,
        runner,
    )?;
    let publication_started = Instant::now();
    let film = Film {
        width: settings.width,
        height: settings.height,
        xyz: result.xyz.ok_or(RenderExecutionError::Internal(
            "fresh render omitted staging film",
        ))?,
        spp_done: settings.spp,
        time_mode: if settings.spp == 0 {
            FilmTimeMode::Uninitialized
        } else {
            result.requested_mode
        },
    };
    let mut report = result.report;
    report.publication_ns = elapsed_ns(publication_started);
    Ok(RenderExecutionOutput { film, report })
}

/// Render the full image (fresh film, samples `[0, spp)`).
pub fn render(scene: &Scene, cx: &Cx<'_>, s: &Settings) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = Film::try_new(s.width, s.height)?;
    render_range(scene, cx, s, &mut film, 0, s.spp)?;
    Ok(film)
}

/// Render a fresh film with deterministic motion blur over `shutter`.
pub fn render_motion(
    scene: &Scene,
    cx: &Cx<'_>,
    s: &Settings,
    shutter: ShutterInterval,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = Film::try_new(s.width, s.height)?;
    render_motion_range(scene, cx, s, &mut film, 0, s.spp, shutter)?;
    Ok(film)
}

/// Render a fresh film with a validated physical/keyframed camera.
pub fn render_cinematic(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    s: &Settings,
    shutter: ShutterInterval,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = Film::try_new(s.width, s.height)?;
    render_cinematic_range(scene, camera, cut_side, cx, s, &mut film, 0, s.spp, shutter)?;
    Ok(film)
}

/// Render a fresh cinematic beauty film with an exactly aligned opt-in AOV
/// profile. Legacy [`render_cinematic`] and [`film_to_exr`] behavior is not
/// changed by this entry point.
#[allow(clippy::too_many_arguments)]
pub fn render_cinematic_with_aovs(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    config: crate::aov::CinematicAovConfig,
) -> Result<CinematicAovFilm, CinematicAovError> {
    cx.checkpoint().map_err(TracerError::from)?;
    let mut film = CinematicAovFilm::try_new(settings.width, settings.height, config)?;
    render_cinematic_range_with_aovs(
        scene,
        camera,
        cut_side,
        cx,
        settings,
        &mut film,
        0,
        settings.spp,
        shutter,
    )?;
    Ok(film)
}

/// Render a fresh deterministic adaptive cinematic film with exactly aligned
/// denoising and diagnostic AOVs.
///
/// Adaptive decisions use only the existing raw XYZ Welford estimator. AOVs
/// observe, but never alter, the accepted `0..terminal_count` sample prefix for
/// each pixel. The complete result remains private until every pixel and a
/// final cancellation checkpoint succeed.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn render_cinematic_adaptive_with_aovs(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    policy: AdaptiveSamplingConfig,
    shutter: ShutterInterval,
    config: crate::aov::CinematicAovConfig,
) -> Result<AdaptiveCinematicAovFilm, CinematicAovError> {
    cx.checkpoint().map_err(TracerError::from)?;
    policy.validate_maximum(settings.spp)?;
    validate_reference_times(config, shutter)?;
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        None,
        0,
        settings.spp,
        Some(shutter),
        camera_path,
    )?;
    let capture_ids = config.captures_ids();
    let palette = CinematicAovPalette::try_from_scene(scene, config.limits(), capture_ids, cx)?;
    let albedo_cache = AovAlbedoCache::try_new(scene.primitives.len(), config.captures_primary())?;
    let continuity_fingerprint = cinematic_input_continuity_fingerprint(scene, camera, cx)?;
    let pixel_count = checked_pixel_len(settings.width, settings.height)?;
    // This constructor performs the complete AOV pixel and retained-memory
    // admission before the larger adaptive beauty allocation below. A request
    // outside its declared resource envelope therefore fails without first
    // consuming the very memory it refused.
    let mut aov =
        AdaptiveAovAccumulator::try_new(settings.width, settings.height, settings.spp, config)?;
    let state_bytes =
        adaptive_state_bytes(pixel_count).map_err(|_| CinematicAovError::AllocationRefused)?;
    let mut beauty_state = AdaptiveRenderState::try_new(pixel_count, state_bytes)
        .map_err(|_| CinematicAovError::AllocationRefused)?;
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let kn = 1.0 / y_integral();
    let capture_primary = config.captures_primary();

    for py in 0..settings.height {
        cx.checkpoint().map_err(TracerError::from)?;
        for px in 0..settings.width {
            let pixel = py
                .checked_mul(settings.width)
                .and_then(|row| row.checked_add(px))
                .ok_or(TracerError::InvalidInput)?;
            let mut accumulator = AdaptivePixelAccumulator::EMPTY;
            for sample in 0..settings.spp {
                cx.checkpoint().map_err(TracerError::from)?;
                let traced = trace_pixel_sample_with_primary(
                    scene,
                    &lighting,
                    cx,
                    settings,
                    kn,
                    sobol.as_ref(),
                    key,
                    pixel,
                    sample,
                    Some(shutter),
                    camera_path,
                    capture_primary,
                )?;
                accumulator.push(traced.xyz)?;
                if capture_primary {
                    let split = traced
                        .contribution_split
                        .ok_or(CinematicAovError::SampleAlignmentMismatch)?;
                    let absolute_time_s =
                        traced.absolute_time_s.ok_or(TracerError::MissingRayTime)?;
                    let primary = prepare_aligned_aov_primary(
                        scene,
                        camera,
                        exposure,
                        cut_side,
                        cx,
                        settings,
                        config.provenance(),
                        &palette,
                        &albedo_cache,
                        capture_ids,
                        absolute_time_s,
                        traced.primary.as_ref(),
                    )?;
                    aov.push(
                        pixel as usize,
                        AlignedAovSample {
                            beauty_xyz: traced.xyz,
                            direct_xyz: split.direct_xyz,
                            indirect_xyz: split.indirect_xyz,
                            emission_xyz: split.emission_xyz,
                            pixel_jitter: traced.pixel_jitter,
                            absolute_sample: sample,
                            primary,
                        },
                    )?;
                }
                if let Some(decision) = accumulator.decision(policy, settings.spp) {
                    accumulator.decision = Some(decision);
                    break;
                }
            }
            let Some(decision) = accumulator.decision else {
                return Err(CinematicAovError::SampleAlignmentMismatch);
            };
            let index = pixel as usize;
            beauty_state.xyz[index] = accumulator.sum_xyz;
            beauty_state.mean_xyz[index] = accumulator.mean_xyz;
            beauty_state.m2_xyz[index] = accumulator.m2_xyz;
            beauty_state.sample_counts[index] = accumulator.samples;
            beauty_state.decisions[index] = decision;
        }
    }
    cx.checkpoint().map_err(TracerError::from)?;
    let beauty = beauty_state.into_film(settings, policy, requested_mode);
    aov.publish(
        beauty,
        adaptive_render_binding(
            *settings,
            shutter,
            exposure.shot_id(),
            cut_side,
            palette,
            continuity_fingerprint,
            policy,
        ),
    )
}

/// Transactionally append one absolute sample range to aligned cinematic
/// beauty and AOV state.
///
/// The accepted primary, material, depth, normals, contribution split, and
/// shutter time all come from the same path traversal. Any cancellation,
/// geometry refusal, AOV refusal, or allocation failure leaves `film`
/// unchanged.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn render_cinematic_range_with_aovs(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut CinematicAovFilm,
    from: u32,
    to: u32,
    shutter: ShutterInterval,
) -> Result<(), CinematicAovError> {
    cx.checkpoint().map_err(TracerError::from)?;
    // A committed AOV prefix is resumable only inside the sample ceiling
    // carried by `Settings`. Refuse before scene/camera admission so a
    // successful public render cannot later become uncheckpointable.
    if to > settings.spp {
        return Err(TracerError::InvalidInput.into());
    }
    if film.config().captures_ids() && to > MAX_EXACT_F32_INTEGER {
        return Err(CinematicAovError::InexactSampleCount { samples: to });
    }
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, requested_mode) = preflight_render(
        scene,
        cx,
        settings,
        Some(film.beauty()),
        from,
        to,
        Some(shutter),
        camera_path,
    )?;
    let capture_ids = film.config().captures_ids();
    let palette =
        CinematicAovPalette::try_from_scene(scene, film.config().limits(), capture_ids, cx)?;
    let continuity_fingerprint = cinematic_input_continuity_fingerprint(scene, camera, cx)?;
    validate_binding(
        film,
        *settings,
        shutter,
        exposure.shot_id(),
        cut_side,
        &palette,
        continuity_fingerprint,
    )?;
    if to == from {
        return Ok(());
    }

    let albedo_cache =
        AovAlbedoCache::try_new(scene.primitives.len(), film.config().captures_primary())?;
    let mut staged = film.try_clone_for_stage(cx)?;
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let kn = 1.0 / y_integral();
    let capture_primary = film.config().captures_primary();
    for py in 0..settings.height {
        cx.checkpoint().map_err(TracerError::from)?;
        for px in 0..settings.width {
            let pixel = py * settings.width + px;
            for sample in from..to {
                cx.checkpoint().map_err(TracerError::from)?;
                let traced = trace_pixel_sample_with_primary(
                    scene,
                    &lighting,
                    cx,
                    settings,
                    kn,
                    sobol.as_ref(),
                    key,
                    pixel,
                    sample,
                    Some(shutter),
                    camera_path,
                    capture_primary,
                )?;
                let beauty = &mut staged.beauty_mut().xyz[pixel as usize];
                beauty[0] += traced.xyz[0];
                beauty[1] += traced.xyz[1];
                beauty[2] += traced.xyz[2];

                if capture_primary {
                    let split = traced
                        .contribution_split
                        .ok_or(CinematicAovError::SampleAlignmentMismatch)?;
                    let absolute_time_s =
                        traced.absolute_time_s.ok_or(TracerError::MissingRayTime)?;
                    let primary = prepare_aligned_aov_primary(
                        scene,
                        camera,
                        exposure,
                        cut_side,
                        cx,
                        settings,
                        film.config().provenance(),
                        &palette,
                        &albedo_cache,
                        capture_ids,
                        absolute_time_s,
                        traced.primary.as_ref(),
                    )?;
                    staged.push(
                        pixel as usize,
                        AlignedAovSample {
                            beauty_xyz: traced.xyz,
                            direct_xyz: split.direct_xyz,
                            indirect_xyz: split.indirect_xyz,
                            emission_xyz: split.emission_xyz,
                            pixel_jitter: traced.pixel_jitter,
                            absolute_sample: sample,
                            primary,
                        },
                    )?;
                }
            }
        }
    }
    cx.checkpoint().map_err(TracerError::from)?;
    staged.beauty_mut().spp_done = to;
    if staged.beauty().time_mode == FilmTimeMode::Uninitialized {
        staged.beauty_mut().time_mode = requested_mode;
    }
    if staged.binding().is_none() {
        staged.bind(render_binding(
            *settings,
            shutter,
            exposure.shot_id(),
            cut_side,
            palette,
            continuity_fingerprint,
        ));
    }
    *film = staged;
    Ok(())
}

type CachedAovAlbedo = Result<Option<[f64; 3]>, ()>;

/// Per-render lazy cache for the scene-linear albedo guide.
///
/// `LiftedSpectrum::rgb` performs the deterministic 80-bin spectral round
/// trip.  The material is immutable for a render, so repeating that work for
/// every primary sample only burns cycles; caching its exact `f64` result does
/// not change estimator or guide bits.  Lazy initialization preserves the
/// previous fail-closed behavior for invalid materials that are never hit.
struct AovAlbedoCache {
    by_primitive: Vec<OnceLock<CachedAovAlbedo>>,
}

impl AovAlbedoCache {
    fn try_new(primitive_count: usize, enabled: bool) -> Result<Self, CinematicAovError> {
        let count = if enabled { primitive_count } else { 0 };
        let mut by_primitive = Vec::new();
        by_primitive
            .try_reserve_exact(count)
            .map_err(|_| CinematicAovError::AllocationRefused)?;
        by_primitive.resize_with(count, OnceLock::new);
        Ok(Self { by_primitive })
    }

    fn get(
        &self,
        primitive_index: usize,
        material: Material,
    ) -> Result<Option<[f64; 3]>, CinematicAovError> {
        let slot = self
            .by_primitive
            .get(primitive_index)
            .ok_or(CinematicAovError::InvalidPrimary)?;
        match *slot.get_or_init(|| match material {
            Material::Lambertian { reflectance } | Material::Ggx { reflectance, .. } => {
                let rgb = reflectance.rgb();
                if rgb.iter().any(|value| !value.is_finite()) {
                    Err(())
                } else {
                    Ok(Some([
                        rgb[0].clamp(0.0, 1.0),
                        rgb[1].clamp(0.0, 1.0),
                        rgb[2].clamp(0.0, 1.0),
                    ]))
                }
            }
            Material::Conductor { .. } | Material::Dielectric { .. } => Ok(None),
        }) {
            Ok(albedo) => Ok(albedo),
            Err(()) => Err(CinematicAovError::InvalidPrimary),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_aligned_aov_primary(
    scene: &Scene,
    camera: &AnimatedCamera,
    exposure: CameraExposure<'_>,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    provenance: crate::aov::CinematicAovProvenance,
    palette: &CinematicAovPalette,
    albedo_cache: &AovAlbedoCache,
    capture_ids: bool,
    absolute_time_s: f64,
    primary: Option<&PrimaryTraceHit>,
) -> Result<Option<AlignedAovPrimary>, CinematicAovError> {
    let Some(primary) = primary else {
        return Ok(None);
    };
    let primitive = scene
        .primitives
        .get(primary.primitive_index)
        .ok_or(CinematicAovError::InvalidPrimary)?;
    let physical = camera
        .evaluate_exposure(cx, exposure, absolute_time_s)
        .map_err(TracerError::from)?;
    let depth_m = match physical
        .project_from_optical_center(primary.hit.point, settings.camera_aspect())
        .map_err(TracerError::from)?
    {
        OpticalCenterProjection::InFront { depth_m, .. } => depth_m,
        OpticalCenterProjection::BehindCamera { .. } => {
            return Err(CinematicAovError::InvalidPrimary);
        }
    };
    let geometric_normal_world = aov_unit_normal(primary.hit.normal)?;
    let shading_normal_world = aov_unit_normal(Some(primary.beauty_shading_normal_world))?;
    let albedo_linear_rgb = albedo_cache.get(primary.primitive_index, primitive.material)?;
    let material_palette_index = if capture_ids {
        palette.material_index(primary.material_identity)?
    } else {
        0
    };
    let object_palette_index = if capture_ids {
        primary
            .surface
            .map(|surface| palette.object_index(surface.identity().object_id()))
            .transpose()?
            .unwrap_or(0)
    } else {
        0
    };
    let previous_motion_pixels = primary
        .surface
        .map(|surface| {
            aligned_previous_motion(
                cx,
                camera,
                exposure,
                cut_side,
                &primitive.shape,
                surface,
                provenance,
                settings,
                absolute_time_s,
                physical.clone(),
            )
        })
        .transpose()?
        .flatten();
    Ok(Some(AlignedAovPrimary {
        primitive_index: primary.primitive_index,
        object_palette_index,
        material_palette_index,
        albedo_linear_rgb,
        geometric_normal_world,
        shading_normal_world,
        // The v1 beauty integrator shades with its face-forwarded geometric
        // frame, not the backend-authored shading normal. Keep the authored
        // validity bit clear until a transport-correct shading-normal design
        // is implemented.
        has_authored_shading_normal: false,
        depth_m,
        previous_motion_pixels,
    }))
}

#[allow(clippy::too_many_arguments)]
fn aligned_previous_motion(
    cx: &Cx<'_>,
    camera: &AnimatedCamera,
    exposure: CameraExposure<'_>,
    cut_side: CutSide,
    shape: &Shape,
    surface: PrimarySurfaceSample,
    provenance: crate::aov::CinematicAovProvenance,
    settings: &Settings,
    absolute_time_s: f64,
    current_camera: PhysicalCamera,
) -> Result<Option<[f64; 2]>, CinematicAovError> {
    let previous_sample_time_s =
        shutter_aligned_previous_reference_time(provenance, absolute_time_s);
    let Some(previous) =
        motion_frame_for_shape(cx, camera, cut_side, shape, previous_sample_time_s, None)?
    else {
        return Err(CinematicAovError::InvalidPrimary);
    };
    let Some(current) = motion_frame_for_shape(
        cx,
        camera,
        cut_side,
        shape,
        absolute_time_s,
        Some((exposure.shot_id(), current_camera)),
    )?
    else {
        return Err(CinematicAovError::InvalidPrimary);
    };
    // The aligned AOV exports only previous motion. Supplying the current
    // frame as the unused next endpoint avoids evaluating beyond a terminal
    // trajectory cut while preserving the low-level three-frame API.
    let next = current.clone();
    let raster =
        RasterSize::try_new(settings.width, settings.height).map_err(TracerError::MotionVector)?;
    match compute_motion_vectors(surface, &previous, &current, &next, raster)
        .map_err(TracerError::MotionVector)?
    {
        MotionVectorComputation::Available(sample) => match sample.previous {
            MotionEndpoint::Projected {
                displacement_pixels,
                ..
            } => Ok(Some(displacement_pixels)),
            MotionEndpoint::CameraCut { .. } | MotionEndpoint::BehindCamera { .. } => Ok(None),
        },
        MotionVectorComputation::Unavailable { .. } => Ok(None),
    }
}

/// Map one accepted shutter sample to the same shutter phase in the preceding
/// presentation interval. The provenance times are presentation timestamps,
/// not necessarily the mechanics times sampled by a front- or back-loaded
/// exposure. Using the timestamps themselves would therefore over- or
/// under-shoot the immediately preceding shutter-integrated image.
fn shutter_aligned_previous_reference_time(
    provenance: crate::aov::CinematicAovProvenance,
    absolute_time_s: f64,
) -> f64 {
    let previous_cadence_s = provenance.frame_time_s() - provenance.previous_frame_time_s();
    absolute_time_s - previous_cadence_s
}

fn motion_frame_for_shape(
    cx: &Cx<'_>,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    shape: &Shape,
    absolute_time_s: f64,
    current_camera: Option<(u64, PhysicalCamera)>,
) -> Result<Option<MotionFrame>, CinematicAovError> {
    let (shot_id, physical) = match current_camera {
        Some(current) => current,
        None => {
            let evaluated = camera
                .evaluate_with_shot(cx, absolute_time_s, cut_side)
                .map_err(TracerError::from)?;
            (evaluated.shot_id(), evaluated.into_camera())
        }
    };
    let frame = match shape {
        Shape::Instance(instance) => {
            MotionFrame::from_instance(absolute_time_s, shot_id, physical, instance)
                .map_err(TracerError::MotionVector)?
        }
        Shape::AnimatedInstance(instance) => {
            let evaluated = instance
                .instance_at(cx, absolute_time_s)
                .map_err(TracerError::from)?;
            MotionFrame::from_instance(absolute_time_s, shot_id, physical, &evaluated)
                .map_err(TracerError::MotionVector)?
        }
        Shape::Mesh(_) | Shape::Chart(_) => return Ok(None),
    };
    Ok(Some(frame))
}

fn aov_unit_normal(normal: Option<Vec3>) -> Result<[f64; 3], CinematicAovError> {
    let normal = normal.ok_or(CinematicAovError::InvalidPrimary)?;
    let scale = normal.x.abs().max(normal.y.abs()).max(normal.z.abs());
    if scale == 0.0 || !scale.is_finite() {
        return Err(CinematicAovError::InvalidPrimary);
    }
    let scaled = [normal.x / scale, normal.y / scale, normal.z / scale];
    let norm =
        scale * (scaled[0] * scaled[0] + scaled[1] * scaled[1] + scaled[2] * scaled[2]).sqrt();
    if norm == 0.0 || !norm.is_finite() {
        return Err(CinematicAovError::InvalidPrimary);
    }
    Ok([normal.x / norm, normal.y / norm, normal.z / norm])
}

const CINEMATIC_AOV_CONTINUITY_DOMAIN: &str =
    "org.frankensim.render.cinematic-aov-continuity-guard.v1";

/// Build a process-local progressive-continuity guard from every borrowed
/// scene and camera value that can affect cinematic path samples. The result is
/// deliberately not an authority-bearing content identity: direct chart
/// backends have no general content-addressing contract, so their allocation
/// address is included only to detect replacement during this process.
fn cinematic_input_continuity_fingerprint(
    scene: &Scene,
    camera: &AnimatedCamera,
    cx: &Cx<'_>,
) -> Result<ContentHash, CinematicAovError> {
    let mut hasher = DomainHasher::new(CINEMATIC_AOV_CONTINUITY_DOMAIN);
    hasher.update(&(scene.primitives.len() as u64).to_le_bytes());
    for (primitive_index, primitive) in scene.primitives.iter().enumerate() {
        if primitive_index.is_multiple_of(1_024) {
            cx.checkpoint().map_err(TracerError::from)?;
        }
        update_shape_continuity(&mut hasher, &primitive.shape, cx)?;
        hasher.update(primitive.material.content_identity().as_bytes());
        match primitive.emission {
            None => hasher.update(&[0]),
            Some((spectrum, scale)) => {
                hasher.update(&[1]);
                update_continuity_f64s(&mut hasher, spectrum.c);
                update_continuity_f64(&mut hasher, scale);
            }
        }
    }
    hasher.update(&(scene.lights.len() as u64).to_le_bytes());
    for light in &scene.lights {
        hasher.update(light.identity().as_bytes());
        hasher.update(&(light.prim as u64).to_le_bytes());
    }
    match &scene.environment {
        None => hasher.update(&[0]),
        Some(environment) => {
            hasher.update(&[1]);
            hasher.update(environment.semantic_hash().as_bytes());
        }
    }
    update_continuity_point(&mut hasher, scene.camera.eye);
    update_continuity_vec(&mut hasher, scene.camera.forward);
    update_continuity_vec(&mut hasher, scene.camera.up);
    update_continuity_f64(&mut hasher, scene.camera.half_tan);

    hasher.update(&(camera.shots().len() as u64).to_le_bytes());
    for shot in camera.shots() {
        hasher.update(&shot.shot_id().to_le_bytes());
        update_continuity_f64(&mut hasher, shot.start_s());
        update_continuity_f64(&mut hasher, shot.end_s());
        hasher.update(&(shot.keyframes().len() as u64).to_le_bytes());
        for keyframe in shot.keyframes() {
            update_continuity_f64(&mut hasher, keyframe.absolute_time_s());
            update_physical_camera_continuity(&mut hasher, keyframe.camera());
            match keyframe.focus() {
                KeyframeFocus::AxialDistance => hasher.update(&[0]),
                KeyframeFocus::WorldPoint(point) => {
                    hasher.update(&[1]);
                    update_continuity_point(&mut hasher, point);
                }
            }
        }
    }
    cx.checkpoint().map_err(TracerError::from)?;
    Ok(hasher.finalize())
}

fn update_shape_continuity(
    hasher: &mut DomainHasher,
    shape: &Shape,
    cx: &Cx<'_>,
) -> Result<(), CinematicAovError> {
    match shape {
        Shape::Mesh(mesh) => {
            hasher.update(&[0]);
            update_mesh_continuity(hasher, mesh, cx)?;
        }
        Shape::Chart(chart) => {
            hasher.update(&[1]);
            let address = std::ptr::from_ref::<dyn Chart>(&**chart).cast::<()>() as usize as u64;
            hasher.update(&address.to_le_bytes());
        }
        Shape::Instance(instance) => {
            hasher.update(&[2]);
            hasher.update(&instance.object_id().to_le_bytes());
            hasher.update(instance.geometry_identity().as_bytes());
            hasher.update(instance.transform().content_identity().as_bytes());
            update_shared_geometry_continuity(hasher, instance.geometry(), cx)?;
        }
        Shape::AnimatedInstance(instance) => {
            hasher.update(&[3]);
            hasher.update(&instance.object_id().to_le_bytes());
            hasher.update(instance.geometry_identity().as_bytes());
            update_shared_geometry_continuity(hasher, instance.geometry(), cx)?;
            hasher.update(&(instance.trajectory().keyframes().len() as u64).to_le_bytes());
            for keyframe in instance.trajectory().keyframes() {
                update_continuity_f64(&mut *hasher, keyframe.absolute_time_s());
                hasher.update(keyframe.transform().content_identity().as_bytes());
                update_continuity_f64s(hasher, keyframe.linear_velocity_m_per_s());
            }
        }
    }
    Ok(())
}

fn update_shared_geometry_continuity(
    hasher: &mut DomainHasher,
    geometry: &SharedGeometry,
    cx: &Cx<'_>,
) -> Result<(), CinematicAovError> {
    match geometry {
        SharedGeometry::Mesh(mesh) => {
            hasher.update(&[0]);
            update_mesh_continuity(hasher, mesh, cx)
        }
        SharedGeometry::Chart(chart) => {
            hasher.update(&[1]);
            let address = std::sync::Arc::as_ptr(chart).cast::<()>() as usize as u64;
            hasher.update(&address.to_le_bytes());
            Ok(())
        }
    }
}

fn update_mesh_continuity(
    hasher: &mut DomainHasher,
    mesh: &TriMesh,
    cx: &Cx<'_>,
) -> Result<(), CinematicAovError> {
    hasher.update(&(mesh.vertices.len() as u64).to_le_bytes());
    for (index, vertex) in mesh.vertices.iter().enumerate() {
        if index.is_multiple_of(1_024) {
            cx.checkpoint().map_err(TracerError::from)?;
        }
        update_continuity_f64s(hasher, *vertex);
    }
    hasher.update(&(mesh.triangles.len() as u64).to_le_bytes());
    for triangle in &mesh.triangles {
        for index in triangle {
            hasher.update(&index.to_le_bytes());
        }
    }
    Ok(())
}

fn update_physical_camera_continuity(hasher: &mut DomainHasher, camera: &PhysicalCamera) {
    update_continuity_point(hasher, camera.eye());
    update_continuity_vec(hasher, camera.forward());
    update_continuity_vec(hasher, camera.up());
    update_continuity_vec(hasher, camera.right());
    let projection = camera.projection();
    update_continuity_f64(hasher, projection.vertical_half_tan());
    match (
        projection.focal_length_m(),
        projection.sensor_height_m(),
        projection.vertical_fov_rad(),
    ) {
        (Some(focal), Some(sensor), None) => {
            hasher.update(&[0]);
            update_continuity_f64(hasher, focal);
            update_continuity_f64(hasher, sensor);
        }
        (None, None, Some(fov)) => {
            hasher.update(&[1]);
            update_continuity_f64(hasher, fov);
        }
        _ => hasher.update(&[2]),
    }
    update_continuity_f64(hasher, camera.focus_distance_m());
    let aperture = camera.aperture();
    let aperture_tag = if aperture.is_pinhole() {
        0
    } else if aperture.blades().is_some() {
        2
    } else {
        1
    };
    hasher.update(&[aperture_tag]);
    update_continuity_f64(hasher, aperture.radius_m());
    hasher.update(&[aperture.blades().unwrap_or(0)]);
    update_continuity_f64(hasher, aperture.rotation_rad().unwrap_or(0.0));
    let exposure = camera.exposure_metadata();
    update_continuity_f64(hasher, exposure.sensitivity_iso());
    update_continuity_f64(hasher, exposure.compensation_ev());
}

fn update_continuity_point(hasher: &mut DomainHasher, point: Point3) {
    update_continuity_f64s(hasher, [point.x, point.y, point.z]);
}

fn update_continuity_vec(hasher: &mut DomainHasher, vector: Vec3) {
    update_continuity_f64s(hasher, [vector.x, vector.y, vector.z]);
}

fn update_continuity_f64s<const N: usize>(hasher: &mut DomainHasher, values: [f64; N]) {
    for value in values {
        update_continuity_f64(hasher, value);
    }
}

fn update_continuity_f64(hasher: &mut DomainHasher, value: f64) {
    let value = if value == 0.0 { 0.0 } else { value };
    hasher.update(&value.to_bits().to_le_bytes());
}

/// Trace one absolute logical cinematic sample and retain the exact primary
/// hit accepted by the beauty path.
///
/// `pixel` is a row-major index into `settings.width × settings.height`, and
/// `sample` is the same absolute sample identity used by progressive renders.
/// The call performs the same scene/camera/shutter admission as
/// [`render_cinematic`] but does not allocate or mutate a [`Film`].
#[allow(clippy::too_many_arguments)]
pub fn trace_cinematic_pixel_sample(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    settings: &Settings,
    shutter: ShutterInterval,
    pixel: u32,
    sample: u32,
) -> Result<CinematicPixelSample, TracerError> {
    cx.checkpoint()?;
    let pixel_count = checked_pixel_len(settings.width, settings.height)?;
    if usize::try_from(pixel).map_or(true, |index| index >= pixel_count) {
        return Err(TracerError::InvalidInput);
    }
    let exposure = camera.admit_shutter(cx, shutter, cut_side)?;
    let camera_path = CameraPath::Cinematic { camera, exposure };
    let (lighting, _) =
        preflight_render(scene, cx, settings, None, 0, 1, Some(shutter), camera_path)?;
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let sobol =
        (settings.sampler == Sampler::OwenSobol).then(|| Sobol::scrambled(3, settings.seed));
    let traced = trace_pixel_sample_with_primary(
        scene,
        &lighting,
        cx,
        settings,
        1.0 / y_integral(),
        sobol.as_ref(),
        key,
        pixel,
        sample,
        Some(shutter),
        camera_path,
        true,
    )?;
    let split = traced.contribution_split.ok_or(TracerError::InvalidInput)?;
    Ok(CinematicPixelSample {
        xyz: traced.xyz,
        absolute_time_s: traced.absolute_time_s.ok_or(TracerError::MissingRayTime)?,
        direct_xyz: split.direct_xyz,
        indirect_xyz: split.indirect_xyz,
        emission_xyz: split.emission_xyz,
        primary: traced.primary,
    })
}

/// The (jitter-x, jitter-y, hero-λ) dimensions for one (pixel, sample).
#[allow(clippy::too_many_arguments)]
fn trace_pixel_sample(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    kn: f64,
    sobol: Option<&Sobol>,
    key: [u32; 2],
    pixel: u32,
    sample: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
) -> Result<[f64; 3], TracerError> {
    Ok(trace_pixel_sample_with_primary(
        scene,
        lighting,
        cx,
        settings,
        kn,
        sobol,
        key,
        pixel,
        sample,
        shutter,
        camera_path,
        false,
    )?
    .xyz)
}

#[allow(clippy::too_many_arguments)]
fn trace_pixel_sample_with_primary(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    kn: f64,
    sobol: Option<&Sobol>,
    key: [u32; 2],
    pixel: u32,
    sample: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
    capture_primary: bool,
) -> Result<PathTraceSample, TracerError> {
    // Pool seeds control scheduling/placement only. Bind downstream chart,
    // camera, and intersection work to the public render seed while retaining
    // the executor's logical tile/iteration and refusal routing.
    let render_cx = cx.with_stream_seed(settings.seed);
    let (jx, jy, ul) = pixel_dims(settings, sobol, key, pixel, sample)?;
    let ray_time = match shutter {
        Some(interval) => {
            let normalized =
                interval.sample_for_stream(settings.seed, u64::from(pixel), u64::from(sample));
            let absolute_time_s = interval.time_at(normalized);
            let mut cached_animated = std::array::from_fn(|_| None);
            let mut cached_count = 0;
            for (primitive_index, primitive) in scene.primitives.iter().enumerate() {
                if let Shape::AnimatedInstance(instance) = &primitive.shape {
                    let Some(slot) = cached_animated.get_mut(cached_count) else {
                        break;
                    };
                    *slot = Some(CachedAnimatedInstance {
                        primitive_index,
                        instance: instance.instance_at(&render_cx, absolute_time_s)?,
                    });
                    cached_count += 1;
                }
            }
            Some(PathTime {
                interval,
                normalized,
                cached_animated,
            })
        }
        None => None,
    };
    trace_path(
        scene,
        lighting,
        &render_cx,
        settings,
        kn,
        pixel,
        sample,
        jx,
        jy,
        ul,
        ray_time.as_ref(),
        camera_path,
        capture_primary,
        None,
    )
}

/// The (jitter-x, jitter-y, hero-λ) dimensions for one (pixel, sample).
fn pixel_dims(
    s: &Settings,
    sobol: Option<&Sobol>,
    key: [u32; 2],
    pixel: u32,
    sample: u32,
) -> Result<(f64, f64, f64), TracerError> {
    match s.sampler {
        Sampler::Iid => {
            let a = philox4x32_10([pixel, sample, 0xdead_0001, 0], key);
            Ok((u32_unit(a[0]), u32_unit(a[1]), u32_unit(a[2])))
        }
        Sampler::OwenSobol => {
            let sobol = sobol.ok_or(TracerError::InvalidInput)?;
            // One Sobol' point per sample index; Cranley–Patterson-free
            // decorrelation across pixels via a per-pixel Philox shift
            // of the SAMPLE INDEX ordering is not net-preserving, so
            // instead the scramble seed is shared and the pixel enters
            // through a Philox-derived toroidal shift of the point —
            // net-preserving per pixel, decorrelated across pixels.
            let mut pt = [0.0f64; 3];
            sobol.point(sample, &mut pt);
            let shift = philox4x32_10([pixel, 0x50b0_1000, 0, 0], key);
            let wrap = |x: f64, u: u32| {
                let v = x + u32_unit(u);
                if v >= 1.0 { v - 1.0 } else { v }
            };
            Ok((
                wrap(pt[0], shift[0]),
                wrap(pt[1], shift[1]),
                wrap(pt[2], shift[2]),
            ))
        }
    }
}

fn unclosed_medium_error(
    boundary_primitive: usize,
    context: &'static str,
    pixel: u32,
    sample: u32,
    depth: u32,
    active_lane: Option<usize>,
    ray: &Ray,
    surface_context: Option<(Vec3, Vec3, f64)>,
) -> TracerError {
    let (geometric_normal_bits, incident_direction_bits, shadow_side_dot_bits) = surface_context
        .map_or((None, None, None), |(normal, incident, side_dot)| {
            (
                Some([normal.x.to_bits(), normal.y.to_bits(), normal.z.to_bits()]),
                Some([
                    incident.x.to_bits(),
                    incident.y.to_bits(),
                    incident.z.to_bits(),
                ]),
                Some(side_dot.to_bits()),
            )
        });
    TracerError::UnclosedMedium {
        boundary_primitive,
        context,
        pixel,
        sample,
        depth,
        active_lane,
        origin_bits: [
            ray.origin.x.to_bits(),
            ray.origin.y.to_bits(),
            ray.origin.z.to_bits(),
        ],
        direction_bits: [
            ray.dir.x.to_bits(),
            ray.dir.y.to_bits(),
            ray.dir.z.to_bits(),
        ],
        geometric_normal_bits,
        incident_direction_bits,
        shadow_side_dot_bits,
    }
}

fn slab_nee_refusal(boundary_primitive: usize, reason: &'static str) -> TracerError {
    TracerError::UnsupportedSlabNee {
        boundary_primitive,
        reason,
    }
}

fn finite_unit_direction(direction: Vec3) -> bool {
    direction.x.is_finite()
        && direction.y.is_finite()
        && direction.z.is_finite()
        && (direction.norm() - 1.0).abs() <= 8.0e-12
}

fn vec_add(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(left.x + right.x, left.y + right.y, left.z + right.z)
}

fn vec_sub(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(left.x - right.x, left.y - right.y, left.z - right.z)
}

fn parallel_directions(left: Vec3, right: Vec3) -> bool {
    finite_unit_direction(left)
        && finite_unit_direction(right)
        && left.dot(right) >= 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE
}

fn mesh_vertex(vertex: [f64; 3]) -> Vec3 {
    Vec3::new(vertex[0], vertex[1], vertex[2])
}

fn normalized_cross(left: Vec3, right: Vec3) -> Option<Vec3> {
    let normal = cross(left, right);
    let norm = normal.norm();
    (norm.is_finite() && norm > 0.0).then(|| normal.scale(1.0 / norm))
}

fn mesh_projection_bounds(mesh: &TriMesh, axis: Vec3) -> Option<(f64, f64, f64)> {
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    let mut coordinate_scale = 1.0_f64;
    for &vertex in &mesh.vertices {
        let point = mesh_vertex(vertex);
        let projection = point.dot(axis);
        if !point.x.is_finite()
            || !point.y.is_finite()
            || !point.z.is_finite()
            || !projection.is_finite()
        {
            return None;
        }
        minimum = minimum.min(projection);
        maximum = maximum.max(projection);
        coordinate_scale = coordinate_scale
            .max(point.x.abs())
            .max(point.y.abs())
            .max(point.z.abs());
    }
    (minimum.is_finite() && maximum.is_finite() && maximum > minimum).then_some((
        minimum,
        maximum,
        coordinate_scale,
    ))
}

/// Projection bounds for two axes in one vertex-memory pass.
///
/// Each axis keeps the same row order and `min`/`max` arithmetic as two calls
/// to `mesh_projection_bounds`; only the redundant vertex load and unused
/// coordinate-scale reductions are removed.
fn mesh_projection_bounds_pair(
    mesh: &TriMesh,
    first_axis: Vec3,
    second_axis: Vec3,
) -> Option<((f64, f64), (f64, f64))> {
    let mut first_minimum = f64::INFINITY;
    let mut first_maximum = f64::NEG_INFINITY;
    let mut second_minimum = f64::INFINITY;
    let mut second_maximum = f64::NEG_INFINITY;
    for &vertex in &mesh.vertices {
        let point = mesh_vertex(vertex);
        let first_projection = point.dot(first_axis);
        let second_projection = point.dot(second_axis);
        if !point.x.is_finite()
            || !point.y.is_finite()
            || !point.z.is_finite()
            || !first_projection.is_finite()
            || !second_projection.is_finite()
        {
            return None;
        }
        first_minimum = first_minimum.min(first_projection);
        first_maximum = first_maximum.max(first_projection);
        second_minimum = second_minimum.min(second_projection);
        second_maximum = second_maximum.max(second_projection);
    }
    (first_minimum.is_finite()
        && first_maximum.is_finite()
        && first_maximum > first_minimum
        && second_minimum.is_finite()
        && second_maximum.is_finite()
        && second_maximum > second_minimum)
        .then_some((
            (first_minimum, first_maximum),
            (second_minimum, second_maximum),
        ))
}

fn triangle_unit_normal(mesh: &TriMesh, triangle_index: u32) -> Option<Vec3> {
    let triangle = *mesh.triangles.get(usize::try_from(triangle_index).ok()?)?;
    let a = mesh_vertex(*mesh.vertices.get(usize::try_from(triangle[0]).ok()?)?);
    let b = mesh_vertex(*mesh.vertices.get(usize::try_from(triangle[1]).ok()?)?);
    let c = mesh_vertex(*mesh.vertices.get(usize::try_from(triangle[2]).ok()?)?);
    normalized_cross(vec_sub(b, a), vec_sub(c, a))
}

fn triangle_is_on_mesh_support_plane(
    mesh: &TriMesh,
    triangle_index: u32,
    axis: Vec3,
    minimum: f64,
    maximum: f64,
    tolerance: f64,
) -> bool {
    let Some(triangle) = mesh
        .triangles
        .get(usize::try_from(triangle_index).unwrap_or(usize::MAX))
    else {
        return false;
    };
    let mut on_minimum = true;
    let mut on_maximum = true;
    for &vertex_index in triangle {
        let Some(&vertex) = mesh
            .vertices
            .get(usize::try_from(vertex_index).unwrap_or(usize::MAX))
        else {
            return false;
        };
        let projection = mesh_vertex(vertex).dot(axis);
        on_minimum &= (projection - minimum).abs() <= tolerance;
        on_maximum &= (projection - maximum).abs() <= tolerance;
    }
    on_minimum || on_maximum
}

fn mesh_face_has_admitted_thin_axis(
    mesh: &TriMesh,
    triangle_index: u32,
    local_geometric_normal: Vec3,
) -> bool {
    // Slab NEE deliberately admits only a support face whose through-mesh
    // extent is strictly smaller than both deterministic orthogonal extents.
    // This cheap O(vertices) test accepts the Euler plate's broad central
    // faces while refusing its vertical and chamfer faces; it is an admission
    // boundary, not a closed-solid or global thickness certificate.
    let Some(face_axis) = triangle_unit_normal(mesh, triangle_index) else {
        return false;
    };
    let local_normal_norm = local_geometric_normal.norm();
    if !local_normal_norm.is_finite() || local_normal_norm <= 0.0 {
        return false;
    }
    let local_geometric_normal = local_geometric_normal.scale(1.0 / local_normal_norm);
    if face_axis.dot(local_geometric_normal).abs() < 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE {
        return false;
    }
    let Some((minimum, maximum, coordinate_scale)) = mesh_projection_bounds(mesh, face_axis) else {
        return false;
    };
    let support_tolerance = 128.0 * f64::EPSILON * coordinate_scale;
    if !triangle_is_on_mesh_support_plane(
        mesh,
        triangle_index,
        face_axis,
        minimum,
        maximum,
        support_tolerance,
    ) {
        return false;
    }
    let candidate_width = maximum - minimum;
    let (tangent, bitangent) = basis_all_sphere(face_axis);
    let Some(((tangent_minimum, tangent_maximum), (bitangent_minimum, bitangent_maximum))) =
        mesh_projection_bounds_pair(mesh, tangent, bitangent)
    else {
        return false;
    };
    let tangent_width = tangent_maximum - tangent_minimum;
    let bitangent_width = bitangent_maximum - bitangent_minimum;
    let comparison_scale = candidate_width
        .max(tangent_width)
        .max(bitangent_width)
        .max(f64::MIN_POSITIVE);
    let comparison_tolerance = SLAB_CONNECTION_REL_TOLERANCE * comparison_scale;
    candidate_width + comparison_tolerance < tangent_width
        && candidate_width + comparison_tolerance < bitangent_width
}

fn require_instanced_mesh_face_witness(
    scene: &Scene,
    intersection: SceneIntersection,
    missing_reason: &'static str,
    non_slab_reason: &'static str,
) -> Result<(), TracerError> {
    let Some(primitive) = scene.primitives.get(intersection.primitive_index) else {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            missing_reason,
        ));
    };
    let geometry = match &primitive.shape {
        Shape::Instance(instance) => instance.geometry(),
        Shape::AnimatedInstance(instance) => instance.geometry(),
        Shape::Mesh(_) | Shape::Chart(_) => {
            return Err(slab_nee_refusal(
                intersection.primitive_index,
                missing_reason,
            ));
        }
    };
    let SharedGeometry::Mesh(mesh) = geometry else {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            missing_reason,
        ));
    };
    let Some(instance_hit) = intersection.instance_hit else {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            missing_reason,
        ));
    };
    let InstanceSurfaceFeature::MeshTriangle { triangle_index, .. } = instance_hit.surface_feature
    else {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            missing_reason,
        ));
    };
    let Some(local_geometric_normal) = instance_hit.local_hit.normal else {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            missing_reason,
        ));
    };
    if !mesh_face_has_admitted_thin_axis(mesh, triangle_index, local_geometric_normal) {
        return Err(slab_nee_refusal(
            intersection.primitive_index,
            non_slab_reason,
        ));
    }
    Ok(())
}

fn discover_parallel_slab(
    scene: &Scene,
    cx: &Cx<'_>,
    first_ray: &Ray,
    first_hit: SceneIntersection,
    ray_time: Option<&PathTime>,
) -> Result<ParallelSlab, TracerError> {
    let boundary_primitive = first_hit.primitive_index;
    require_instanced_mesh_face_witness(
        scene,
        first_hit,
        "first slab interface lacks an instanced triangle-mesh face witness",
        "first slab interface lacks admitted thin-axis support",
    )?;
    let (glass, surface) = match scene.primitives[boundary_primitive].material {
        Material::Dielectric { glass, surface } => (glass, surface),
        Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
            return Err(slab_nee_refusal(
                boundary_primitive,
                "first connection blocker is not dielectric",
            ));
        }
    };
    if !surface.is_delta() {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "first dielectric blocker is rough",
        ));
    }
    let entry_frame = surface_frame(&first_hit.hit, first_ray)?;
    if !entry_frame.entering {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "first interface is not an ambient-to-slab entry",
        ));
    }
    let axis = entry_frame.geometric.scale(-1.0);
    let probe = Ray {
        origin: dielectric_spawn_origin(first_hit.hit.point, entry_frame.geometric, axis),
        dir: axis,
    };
    let Some(opposite) = intersect(scene, cx, &probe, ray_time)? else {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "normal probe escaped before finding the opposite interface",
        ));
    };
    if opposite.primitive_index != boundary_primitive {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "normal probe encountered overlapping or nested geometry",
        ));
    }
    require_instanced_mesh_face_witness(
        scene,
        opposite,
        "opposite slab interface lacks an instanced triangle-mesh face witness",
        "opposite slab interface lacks admitted thin-axis support",
    )?;
    let exit_frame = surface_frame(&opposite.hit, &probe)?;
    if exit_frame.entering
        || entry_frame.geometric.dot(exit_frame.geometric) > -1.0 + SLAB_PARALLEL_COSINE_TOLERANCE
    {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "opposite interface is not parallel with reversed orientation",
        ));
    }
    let separation = opposite.hit.point.delta_from(first_hit.hit.point);
    let thickness_m = separation.dot(axis);
    let tangential = vec_sub(separation, axis.scale(thickness_m));
    let scale = thickness_m.abs().max(
        first_hit
            .hit
            .point
            .delta_from(Point3::new(0.0, 0.0, 0.0))
            .norm()
            * f64::EPSILON,
    );
    if !thickness_m.is_finite()
        || thickness_m <= 2.0 * RAY_EPS
        || tangential.norm() > SLAB_CONNECTION_REL_TOLERANCE * scale.max(RAY_EPS)
    {
        return Err(slab_nee_refusal(
            boundary_primitive,
            "normal probe did not establish a positive planar thickness",
        ));
    }
    Ok(ParallelSlab {
        boundary_primitive,
        glass,
        entry_reference: first_hit.hit.point,
        exit_reference: opposite.hit.point,
        axis,
        thickness_m,
    })
}

/// Solve the unique axisymmetric two-interface Snell connection to a finite
/// point.  A parallel slab returns the outgoing ray to the incident external
/// direction, but shifts its line laterally by a direction-dependent amount;
/// merely aiming the shadow ray at the light point is therefore wrong.
fn slab_connection_to_point(
    source: Point3,
    target: Point3,
    light_normal: Vec3,
    slab_axis: Vec3,
    thickness_m: f64,
    eta_ambient_over_glass: f64,
) -> Result<SlabConnectionGeometry, &'static str> {
    if !finite_unit_direction(slab_axis)
        || !finite_unit_direction(light_normal)
        || !thickness_m.is_finite()
        || thickness_m <= 0.0
        || !eta_ambient_over_glass.is_finite()
        || !(0.0..=1.0).contains(&eta_ambient_over_glass)
        || eta_ambient_over_glass == 0.0
    {
        return Err("invalid slab connection parameters");
    }
    let displacement = target.delta_from(source);
    let normal_distance = displacement.dot(slab_axis);
    let external_normal_distance = normal_distance - thickness_m;
    if !normal_distance.is_finite() || external_normal_distance <= 0.0 {
        return Err("source and finite light do not bracket the slab in ambient");
    }
    let tangent = vec_sub(displacement, slab_axis.scale(normal_distance));
    let tangent_distance = tangent.norm();
    if !tangent_distance.is_finite() {
        return Err("finite-light tangential displacement is invalid");
    }

    let ratio = eta_ambient_over_glass;
    let radial_residual = |sine: f64| {
        let cosine = (1.0 - sine * sine).max(0.0).sqrt();
        let internal_cosine = (1.0 - ratio * ratio * sine * sine).max(0.0).sqrt();
        external_normal_distance * sine / cosine + thickness_m * ratio * sine / internal_cosine
            - tangent_distance
    };
    let sine = if tangent_distance == 0.0 {
        0.0
    } else {
        let mut low = 0.0;
        let mut high = 1.0_f64.next_down();
        let high_residual = radial_residual(high);
        if !high_residual.is_finite() || high_residual <= 0.0 {
            return Err("finite-light Snell connection has no propagating root");
        }
        for _ in 0..SLAB_CONNECTION_BISECTION_STEPS {
            let middle = 0.5 * (low + high);
            if radial_residual(middle) < 0.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    };
    let cosine = (1.0 - sine * sine).max(0.0).sqrt();
    let internal_sine = ratio * sine;
    let internal_cosine = (1.0 - internal_sine * internal_sine).max(0.0).sqrt();
    if cosine <= 0.0 || internal_cosine <= 0.0 {
        return Err("finite-light Snell connection is grazing or evanescent");
    }
    let tangent_direction = if tangent_distance == 0.0 {
        basis_all_sphere(slab_axis).0
    } else {
        tangent.scale(1.0 / tangent_distance)
    };
    let incident_direction = vec_add(tangent_direction.scale(sine), slab_axis.scale(cosine));
    let internal_direction = vec_add(
        tangent_direction.scale(internal_sine),
        slab_axis.scale(internal_cosine),
    );

    // The outgoing external ray is parallel to `incident_direction` and can
    // be represented by a direction-dependent virtual-origin shift
    //   delta(w) = h [r/c_g - 1/c] w_t.
    // Differentiate that shift in an orthonormal tangent basis of the source
    // sphere.  Projection of the two ray differentials onto the emitter plane
    // gives dA/dOmega = |w.(a1 x a2)| / |n_L.w|.
    let tangent_component = vec_sub(
        incident_direction,
        slab_axis.scale(incident_cosine(incident_direction, slab_axis)),
    );
    let shift_coefficient = ratio / internal_cosine - 1.0 / cosine;
    let shift_derivative = 1.0 / (cosine * cosine)
        - ratio * ratio * ratio * cosine / (internal_cosine * internal_cosine * internal_cosine);
    let ray_parameter = normal_distance / cosine;
    let shift_differential = |sphere_tangent: Vec3| {
        let dc = slab_axis.dot(sphere_tangent);
        let dt = vec_sub(sphere_tangent, slab_axis.scale(dc));
        vec_add(
            dt.scale(thickness_m * shift_coefficient),
            tangent_component.scale(thickness_m * shift_derivative * dc),
        )
    };
    let (sphere_tangent, sphere_bitangent) = basis_all_sphere(incident_direction);
    let differential_u = vec_add(
        sphere_tangent.scale(ray_parameter),
        shift_differential(sphere_tangent),
    );
    let differential_v = vec_add(
        sphere_bitangent.scale(ray_parameter),
        shift_differential(sphere_bitangent),
    );
    let light_cosine = light_normal.dot(incident_direction).abs();
    if light_cosine <= 1.0e-12 {
        return Err("connected finite light is grazing");
    }
    let light_area_jacobian = incident_direction
        .dot(cross(differential_u, differential_v))
        .abs()
        / light_cosine;
    if !light_area_jacobian.is_finite() || light_area_jacobian <= 0.0 {
        return Err("finite-light slab Jacobian is singular");
    }

    let lateral_shift = tangent_component.scale(thickness_m * shift_coefficient);
    let reconstructed = source
        .offset(lateral_shift)
        .offset(incident_direction.scale(ray_parameter));
    let residual = reconstructed.delta_from(target).norm();
    let scale = normal_distance
        .abs()
        .max(tangent_distance)
        .max(thickness_m)
        .max(RAY_EPS);
    if !residual.is_finite() || residual > SLAB_CONNECTION_REL_TOLERANCE * scale {
        return Err("finite-light Snell solve failed its endpoint residual");
    }
    Ok(SlabConnectionGeometry {
        incident_direction,
        internal_direction,
        light_area_jacobian: Some(light_area_jacobian),
    })
}

fn incident_cosine(direction: Vec3, slab_axis: Vec3) -> f64 {
    direction.dot(slab_axis).clamp(0.0, 1.0)
}

fn slab_connection_to_environment(
    direction: Vec3,
    slab_axis: Vec3,
    eta_ambient_over_glass: f64,
) -> Result<SlabConnectionGeometry, &'static str> {
    if !finite_unit_direction(direction)
        || !finite_unit_direction(slab_axis)
        || !eta_ambient_over_glass.is_finite()
        || !(0.0..=1.0).contains(&eta_ambient_over_glass)
        || eta_ambient_over_glass == 0.0
    {
        return Err("invalid environment slab connection parameters");
    }
    let cosine = direction.dot(slab_axis);
    if cosine <= 0.0 {
        return Err("sampled environment direction does not cross the slab");
    }
    let tangent = vec_sub(direction, slab_axis.scale(cosine));
    let internal_tangent = tangent.scale(eta_ambient_over_glass);
    let internal_cosine = (1.0 - internal_tangent.dot(internal_tangent))
        .max(0.0)
        .sqrt();
    if internal_cosine <= 0.0 {
        return Err("sampled environment direction has no propagating slab path");
    }
    Ok(SlabConnectionGeometry {
        incident_direction: direction,
        internal_direction: vec_add(internal_tangent, slab_axis.scale(internal_cosine)),
        // A parallel slab embedded in the same ambient medium returns every
        // external direction unchanged, so its environment angular Jacobian
        // is exactly one.
        light_area_jacobian: None,
    })
}

fn deterministic_snell_transmission(
    normal_toward_incident: Vec3,
    travel_direction: Vec3,
    eta_i: f64,
    eta_t: f64,
) -> Result<(Vec3, f64), TracerError> {
    let incident_cosine = normal_toward_incident
        .dot(travel_direction.scale(-1.0))
        .clamp(0.0, 1.0);
    let fresnel = fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
    if fresnel.total_internal_reflection {
        return Err(TracerError::Dielectric(DielectricError::InvalidInterface));
    }
    let eta = eta_i / eta_t;
    let direction = vec_add(
        normal_toward_incident.scale(eta * incident_cosine - fresnel.transmitted_cosine),
        travel_direction.scale(eta),
    );
    let norm = direction.norm();
    if !norm.is_finite() || norm <= 0.0 {
        return Err(TracerError::Dielectric(DielectricError::InvalidDirection));
    }
    let direction = direction.scale(1.0 / norm);
    if normal_toward_incident.dot(direction) >= 0.0 {
        return Err(TracerError::Dielectric(DielectricError::InvalidDirection));
    }
    Ok((direction, 1.0 - fresnel.reflectance))
}

fn primitive_is_dielectric(scene: &Scene, primitive_index: usize) -> bool {
    scene
        .primitives
        .get(primitive_index)
        .is_some_and(|primitive| matches!(primitive.material, Material::Dielectric { .. }))
}

#[allow(clippy::too_many_arguments)]
fn trace_parallel_slab_direct_lane(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    source_point: Point3,
    source_geometric_normal: Vec3,
    slab: ParallelSlab,
    direct: PreparedDirectLight,
    wavelength_nm: f64,
) -> Result<SlabDirectLane, TracerError> {
    let eta_glass = slab.glass.ior().eval(wavelength_nm)?;
    let ratio = 1.0 / eta_glass;
    let geometry = match direct.target {
        DirectLightTarget::Rectangle { point, normal, .. } => slab_connection_to_point(
            source_point,
            point,
            normal,
            slab.axis,
            slab.thickness_m,
            ratio,
        )
        .map_err(|reason| slab_nee_refusal(slab.boundary_primitive, reason))?,
        DirectLightTarget::Environment => {
            slab_connection_to_environment(direct.direction, slab.axis, ratio)
                .map_err(|reason| slab_nee_refusal(slab.boundary_primitive, reason))?
        }
    };
    let nee_pdf_solid_angle = match direct.target {
        DirectLightTarget::Rectangle {
            distance_m, normal, ..
        } => {
            let direct_light_cosine = normal.dot(direct.direction).abs();
            let area_pdf = direct.pdf_solid_angle * direct_light_cosine / (distance_m * distance_m);
            area_pdf
                * geometry
                    .light_area_jacobian
                    .ok_or(TracerError::InvalidInput)?
        }
        DirectLightTarget::Environment => direct.pdf_solid_angle,
    };
    if !nee_pdf_solid_angle.is_finite() || nee_pdf_solid_angle <= 0.0 {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected light density is non-positive or non-finite",
        ));
    }

    let entry_ray = Ray {
        origin: dielectric_spawn_origin(
            source_point,
            source_geometric_normal,
            geometry.incident_direction,
        ),
        dir: geometry.incident_direction,
    };
    let Some(entry) = intersect(scene, cx, &entry_ray, ray_time)? else {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected source ray missed the discovered slab",
        ));
    };
    if entry.primitive_index != slab.boundary_primitive {
        if primitive_is_dielectric(scene, entry.primitive_index) {
            return Err(slab_nee_refusal(
                slab.boundary_primitive,
                "connected source ray encountered a nested dielectric",
            ));
        }
        return Ok(SlabDirectLane {
            incident_direction: geometry.incident_direction,
            nee_pdf_solid_angle,
            transmission_probability: 0.0,
            radiance_transport: 0.0,
            visible: false,
        });
    }
    require_instanced_mesh_face_witness(
        scene,
        entry,
        "connected slab entry lacks an instanced triangle-mesh face witness",
        "connected slab entry lacks admitted thin-axis support",
    )?;
    let entry_frame = surface_frame(&entry.hit, &entry_ray)?;
    if !entry_frame.entering
        || entry_frame.geometric.dot(slab.axis.scale(-1.0)) < 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE
    {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected entry lies on a nonparallel slab face",
        ));
    }
    let source_clearance = entry.hit.point.delta_from(source_point).dot(slab.axis);
    let target_clearance = match direct.target {
        DirectLightTarget::Rectangle { point, .. } => {
            point.delta_from(slab.exit_reference).dot(slab.axis)
        }
        DirectLightTarget::Environment => f64::INFINITY,
    };
    if source_clearance <= 0.0 || target_clearance <= 0.0 {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected source and target do not lie in opposite ambient half-spaces",
        ));
    }

    let (internal_direction, entry_transmission_probability) = deterministic_snell_transmission(
        entry_frame.geometric,
        geometry.incident_direction,
        1.0,
        eta_glass,
    )?;
    if !parallel_directions(internal_direction, geometry.internal_direction) {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected entry refraction disagrees with the slab solve",
        ));
    }
    let internal_ray = Ray {
        origin: dielectric_spawn_origin(entry.hit.point, entry_frame.geometric, internal_direction),
        dir: internal_direction,
    };
    let Some(exit) = intersect(scene, cx, &internal_ray, ray_time)? else {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected internal ray escaped the declared closed slab",
        ));
    };
    if exit.primitive_index != slab.boundary_primitive {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected internal ray encountered overlapping or nested geometry",
        ));
    }
    require_instanced_mesh_face_witness(
        scene,
        exit,
        "connected slab exit lacks an instanced triangle-mesh face witness",
        "connected slab exit lacks admitted thin-axis support",
    )?;
    let exit_frame = surface_frame(&exit.hit, &internal_ray)?;
    if exit_frame.entering
        || entry_frame.geometric.dot(exit_frame.geometric) > -1.0 + SLAB_PARALLEL_COSINE_TOLERANCE
    {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected exit is not the parallel opposite interface",
        ));
    }
    let internal_segment = exit.hit.point.delta_from(entry.hit.point);
    let actual_thickness = internal_segment.dot(slab.axis);
    let thickness_scale = slab.thickness_m.max(RAY_EPS);
    if (actual_thickness - slab.thickness_m).abs() > SLAB_CONNECTION_REL_TOLERANCE * thickness_scale
    {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "connected path observes a nonuniform slab thickness",
        ));
    }
    let (outgoing_direction, exit_transmission_probability) = deterministic_snell_transmission(
        exit_frame.geometric.scale(-1.0),
        internal_direction,
        eta_glass,
        1.0,
    )?;
    if !parallel_directions(outgoing_direction, geometry.incident_direction) {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "two-interface refraction did not restore the external direction",
        ));
    }
    let exit_ray = Ray {
        origin: dielectric_spawn_origin(exit.hit.point, exit_frame.geometric, outgoing_direction),
        dir: outgoing_direction,
    };
    let final_hit = intersect(scene, cx, &exit_ray, ray_time)?;
    let visible = match (direct.target, final_hit) {
        (
            DirectLightTarget::Rectangle {
                primitive_index,
                point,
                ..
            },
            Some(hit),
        ) if hit.primitive_index == primitive_index => {
            let point_error = hit.hit.point.delta_from(point).norm();
            let distance_scale = point.delta_from(exit.hit.point).norm().max(RAY_EPS);
            if point_error > 4.0 * RAY_EPS + SLAB_CONNECTION_REL_TOLERANCE * distance_scale {
                return Err(slab_nee_refusal(
                    slab.boundary_primitive,
                    "connected ray reached the sampled emitter at the wrong point",
                ));
            }
            true
        }
        (DirectLightTarget::Rectangle { .. }, Some(hit))
        | (DirectLightTarget::Environment, Some(hit)) => {
            if primitive_is_dielectric(scene, hit.primitive_index) {
                return Err(slab_nee_refusal(
                    slab.boundary_primitive,
                    "post-slab connection encountered another dielectric",
                ));
            }
            false
        }
        (DirectLightTarget::Rectangle { .. }, None) => false,
        (DirectLightTarget::Environment, None) => true,
    };
    let internal_distance_m = internal_segment.norm();
    let beer = medium_transmittance(Some(slab.glass), wavelength_nm, internal_distance_m)?;
    let transmission_probability = entry_transmission_probability * exit_transmission_probability;
    // Each delta BTDF carries its radiance-mode eta factor.  For the same
    // ambient on both sides the reciprocal factors cancel analytically, but
    // retaining both terms pins the transport convention and catches a future
    // one-interface omission.
    let entry_btdf = entry_transmission_probability * ratio * ratio;
    let exit_eta_ratio = eta_glass;
    let exit_btdf = exit_transmission_probability * exit_eta_ratio * exit_eta_ratio;
    let radiance_transport = entry_btdf * exit_btdf * beer;
    if !transmission_probability.is_finite()
        || transmission_probability < 0.0
        || !radiance_transport.is_finite()
        || radiance_transport < 0.0
    {
        return Err(slab_nee_refusal(
            slab.boundary_primitive,
            "slab Fresnel or Beer-Lambert transport became invalid",
        ));
    }
    Ok(SlabDirectLane {
        incident_direction: geometry.incident_direction,
        nee_pdf_solid_angle,
        transmission_probability,
        radiance_transport,
        visible,
    })
}

#[allow(clippy::too_many_arguments)]
fn previous_after_dielectric_sample(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    previous: Option<PreviousBsdf>,
    previous_origin: Point3,
    medium_stack: &MediumStack,
    boundary_primitive: usize,
    glass: DielectricGlass,
    surface: DielectricSurface,
    frame: SurfaceFrame,
    hit: Hit,
    instance_hit: Option<InstanceHit>,
    incident_ray: Ray,
    sampled_direction: Vec3,
    sampled_event: DielectricEvent,
    sampled_pdf: f64,
    sampled_delta: bool,
    wavelength_nm: f64,
) -> Result<PreviousBsdf, TracerError> {
    let mut next = PreviousBsdf {
        pdf: sampled_pdf,
        delta: sampled_delta,
        opaque_source_geometric_normal: None,
        smooth_slab: None,
    };
    if !surface.is_delta() || sampled_event != DielectricEvent::Transmission {
        return Ok(next);
    }
    let boundary = boundary_media(boundary_primitive, glass, frame.entering, medium_stack)?;
    let eta_i = medium_ior(boundary.incident, wavelength_nm)?;
    let eta_t = medium_ior(boundary.transmitted, wavelength_nm)?;
    let (_, transmission_probability) =
        deterministic_snell_transmission(frame.oriented, incident_ray.dir, eta_i, eta_t)?;

    if frame.entering && medium_stack.len() == 0 {
        let Some(source) = previous.filter(|source| {
            source.opaque_source_geometric_normal.is_some()
                && !source.delta
                && source.smooth_slab.is_none()
        }) else {
            return Ok(next);
        };
        let first_hit = SceneIntersection {
            primitive_index: boundary_primitive,
            hit,
            instance_hit,
        };
        let slab = match discover_parallel_slab(scene, cx, &incident_ray, first_hit, ray_time) {
            Ok(slab) => slab,
            // Slab NEE is optional proposal support, not a restriction on
            // ordinary BSDF transport. Retain the delta event without a slab
            // competitor so an eventual emissive hit receives unit MIS weight.
            Err(TracerError::UnsupportedSlabNee { .. }) => return Ok(next),
            Err(error) => return Err(error),
        };
        next.smooth_slab = Some(SmoothSlabPath::Entered(SmoothSlabEntry {
            source_origin: previous_origin,
            source_geometric_normal: source
                .opaque_source_geometric_normal
                .ok_or(TracerError::InvalidInput)?,
            source_direction: incident_ray.dir,
            source_pdf_solid_angle: source.pdf,
            slab,
            entry_transmission_probability: transmission_probability,
        }));
        return Ok(next);
    }

    if !frame.entering
        && medium_stack.len() == 1
        && let Some(SmoothSlabPath::Entered(entry)) = previous.and_then(|value| value.smooth_slab)
        && entry.slab.boundary_primitive == boundary_primitive
        && entry.slab.glass == glass
        && frame.geometric.dot(entry.slab.axis) >= 1.0 - SLAB_PARALLEL_COSINE_TOLERANCE
        && parallel_directions(sampled_direction, entry.source_direction)
    {
        require_instanced_mesh_face_witness(
            scene,
            SceneIntersection {
                primitive_index: boundary_primitive,
                hit,
                instance_hit,
            },
            "sampled slab exit lacks an instanced triangle-mesh face witness",
            "sampled slab exit lacks admitted thin-axis support",
        )?;
        let observed_thickness = hit
            .point
            .delta_from(entry.slab.entry_reference)
            .dot(entry.slab.axis);
        if (observed_thickness - entry.slab.thickness_m).abs()
            <= SLAB_CONNECTION_REL_TOLERANCE * entry.slab.thickness_m.max(RAY_EPS)
        {
            next.smooth_slab = Some(SmoothSlabPath::Exited(SmoothSlabExit {
                source_origin: entry.source_origin,
                source_geometric_normal: entry.source_geometric_normal,
                source_direction: entry.source_direction,
                source_pdf_solid_angle: entry.source_pdf_solid_angle,
                slab: entry.slab,
                transmission_probability: entry.entry_transmission_probability
                    * transmission_probability,
            }));
        }
    }
    Ok(next)
}

fn same_local_parallel_face_pair(left: ParallelSlab, right: ParallelSlab) -> bool {
    if left.boundary_primitive != right.boundary_primitive
        || left.glass != right.glass
        || !parallel_directions(left.axis, right.axis)
    {
        return false;
    }
    let scale = left
        .thickness_m
        .abs()
        .max(right.thickness_m.abs())
        .max(RAY_EPS);
    let tolerance = 4.0 * RAY_EPS + SLAB_CONNECTION_REL_TOLERANCE * scale;
    (left.thickness_m - right.thickness_m).abs() <= tolerance
        && right
            .entry_reference
            .delta_from(left.entry_reference)
            .dot(left.axis)
            .abs()
            <= tolerance
        && right
            .exit_reference
            .delta_from(left.exit_reference)
            .dot(left.axis)
            .abs()
            <= tolerance
}

/// Replay the exact straight shadow-ray eligibility test performed before the
/// forward slab-NEE branch. A completed BSDF slab path has no competing NEE
/// density when that ray misses the slab or first meets any other primitive.
/// If it first meets the same dielectric primitive, its local planar face pair
/// must still admit before the forward technique has a nonzero density. A
/// refusal therefore removes the competing NEE proposal instead of invalidating
/// the completed BSDF path.
fn replay_forward_slab_eligibility(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    slab_path: SmoothSlabExit,
    straight_direction: Vec3,
) -> Result<Option<ParallelSlab>, TracerError> {
    let source_normal = slab_path.source_geometric_normal;
    let source_normal_norm = source_normal.norm();
    if !source_normal.x.is_finite()
        || !source_normal.y.is_finite()
        || !source_normal.z.is_finite()
        || !source_normal_norm.is_finite()
        || source_normal_norm <= 0.0
        || !finite_unit_direction(straight_direction)
    {
        return Err(slab_nee_refusal(
            slab_path.slab.boundary_primitive,
            "reverse slab MIS retained an invalid source frame or direction",
        ));
    }
    if slab_path.source_geometric_normal.dot(straight_direction) <= 0.0 {
        return Ok(None);
    }
    let straight_shadow = Ray {
        origin: slab_path
            .source_origin
            .offset(slab_path.source_geometric_normal.scale(RAY_EPS)),
        dir: straight_direction,
    };
    let Some(first_blocker) = intersect(scene, cx, &straight_shadow, ray_time)? else {
        return Ok(None);
    };
    if first_blocker.primitive_index != slab_path.slab.boundary_primitive {
        return Ok(None);
    }
    let replayed =
        match discover_parallel_slab(scene, cx, &straight_shadow, first_blocker, ray_time) {
            Ok(replayed) => replayed,
            Err(TracerError::UnsupportedSlabNee { .. }) => return Ok(None),
            Err(error) => return Err(error),
        };
    if !same_local_parallel_face_pair(slab_path.slab, replayed) {
        return Err(slab_nee_refusal(
            slab_path.slab.boundary_primitive,
            "forward straight shadow and BSDF path do not cross the same local planar face pair",
        ));
    }
    Ok(Some(replayed))
}

fn completed_slab_rectangle_nee_pdf(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    slab_path: SmoothSlabExit,
    light_index: usize,
    light_point: Point3,
    wavelength_nm: f64,
) -> Result<f64, TracerError> {
    let light = scene
        .lights
        .get(light_index)
        .ok_or(TracerError::InvalidInput)?;
    let direct_pdf = lighting.rect_mixture_pdf(light_index, slab_path.source_origin, light_point);
    if direct_pdf == 0.0 {
        return Ok(0.0);
    }
    let displacement = light_point.delta_from(slab_path.source_origin);
    let distance_squared = displacement.dot(displacement);
    if !(distance_squared > 0.0 && distance_squared.is_finite()) {
        return Err(TracerError::InvalidInput);
    }
    let direct_direction = displacement.scale(1.0 / distance_squared.sqrt());
    let Some(replayed_slab) =
        replay_forward_slab_eligibility(scene, cx, ray_time, slab_path, direct_direction)?
    else {
        return Ok(0.0);
    };
    let distance_m = distance_squared.sqrt();
    let direct = PreparedDirectLight {
        direction: direct_direction,
        emission: light.emission,
        pdf_solid_angle: direct_pdf,
        target: DirectLightTarget::Rectangle {
            primitive_index: light.prim,
            distance_m,
            point: light_point,
            normal: light.normal(),
        },
    };
    let connection = match trace_parallel_slab_direct_lane(
        scene,
        cx,
        ray_time,
        slab_path.source_origin,
        slab_path.source_geometric_normal,
        replayed_slab,
        direct,
        wavelength_nm,
    ) {
        Ok(connection) => connection,
        Err(TracerError::UnsupportedSlabNee { .. }) => return Ok(0.0),
        Err(error) => return Err(error),
    };
    if !connection.visible {
        return Ok(0.0);
    }
    if !parallel_directions(connection.incident_direction, slab_path.source_direction) {
        return Err(slab_nee_refusal(
            slab_path.slab.boundary_primitive,
            "BSDF slab path disagrees with the corresponding NEE connection",
        ));
    }
    Ok(connection.nee_pdf_solid_angle)
}

#[allow(clippy::too_many_arguments)]
fn completed_slab_environment_nee_pdf(
    scene: &Scene,
    cx: &Cx<'_>,
    ray_time: Option<&PathTime>,
    slab_path: SmoothSlabExit,
    direction: Vec3,
    direct_pdf: f64,
    wavelength_nm: f64,
) -> Result<f64, TracerError> {
    if !direct_pdf.is_finite() || direct_pdf <= 0.0 {
        return Ok(0.0);
    }
    let Some(replayed_slab) =
        replay_forward_slab_eligibility(scene, cx, ray_time, slab_path, direction)?
    else {
        return Ok(0.0);
    };
    let direct = PreparedDirectLight {
        direction,
        // The transport validator does not consume radiance. Retain a finite
        // placeholder so reverse MIS can replay exactly the same geometric
        // support test as the forward environment proposal.
        emission: (LiftedSpectrum { c: [0.0; 3] }, 1.0),
        pdf_solid_angle: direct_pdf,
        target: DirectLightTarget::Environment,
    };
    let connection = match trace_parallel_slab_direct_lane(
        scene,
        cx,
        ray_time,
        slab_path.source_origin,
        slab_path.source_geometric_normal,
        replayed_slab,
        direct,
        wavelength_nm,
    ) {
        Ok(connection) => connection,
        Err(TracerError::UnsupportedSlabNee { .. }) => return Ok(0.0),
        Err(error) => return Err(error),
    };
    if !connection.visible {
        return Ok(0.0);
    }
    if !parallel_directions(connection.incident_direction, slab_path.source_direction) {
        return Err(slab_nee_refusal(
            replayed_slab.boundary_primitive,
            "BSDF environment slab path disagrees with the corresponding NEE connection",
        ));
    }
    Ok(connection.nee_pdf_solid_angle)
}

fn completed_slab_mis_weight(
    strategy: DirectStrategy,
    slab_path: SmoothSlabExit,
    nee_pdf_solid_angle: f64,
) -> f64 {
    match strategy {
        DirectStrategy::BsdfOnly => 1.0,
        DirectStrategy::NeeOnly => 0.0,
        DirectStrategy::Mis => balance_heuristic(
            1,
            slab_path.source_pdf_solid_angle * slab_path.transmission_probability,
            1,
            nee_pdf_solid_angle,
        ),
    }
}

fn slab_nee_mis_weight(
    strategy: DirectStrategy,
    nee_pdf_solid_angle: f64,
    source_bsdf_pdf_solid_angle: f64,
    transmission_probability: f64,
    competing_bsdf_path_is_evaluated: bool,
) -> Result<f64, TracerError> {
    match strategy {
        DirectStrategy::NeeOnly => Ok(1.0),
        DirectStrategy::BsdfOnly => Err(TracerError::InvalidInput),
        DirectStrategy::Mis if !competing_bsdf_path_is_evaluated => Ok(1.0),
        DirectStrategy::Mis => Ok(balance_heuristic(
            1,
            nee_pdf_solid_angle,
            1,
            source_bsdf_pdf_solid_angle * transmission_probability,
        )),
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one integrator, one story
fn trace_path(
    scene: &Scene,
    lighting: &AdmittedLighting<'_>,
    cx: &Cx<'_>,
    s: &Settings,
    kn: f64,
    pixel: u32,
    sample: u32,
    jx: f64,
    jy: f64,
    ul: f64,
    ray_time: Option<&PathTime>,
    camera_path: CameraPath<'_>,
    capture_primary: bool,
    resume: Option<SpectralPathState>,
) -> Result<PathTraceSample, TracerError> {
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let rng = PathRng {
        pixel,
        sample,
        dim: 1,
        key,
    };
    // Hero wavelengths: one stratified draw covers the packet.
    let hero = LAMBDA_MIN + ul * (LAMBDA_MAX - LAMBDA_MIN);
    let lambdas = hero_wavelengths(hero, PACKET, LAMBDA_MIN, LAMBDA_MAX);
    let state = if let Some(state) = resume {
        state
    } else {
        // Camera ray. Keep the legacy branch expression-for-expression compatible;
        // the opt-in cinematic branch owns separate lens dimensions and evaluates
        // at the same absolute time already carried by animated geometry.
        let px = pixel % s.width;
        let py = pixel / s.width;
        let (w, h) = (f64::from(s.width), f64::from(s.height));
        let ray = match camera_path {
            CameraPath::Legacy => {
                let ndc_x = (2.0 * (f64::from(px) + jx) / w - 1.0)
                    * s.camera_aspect()
                    * scene.camera.half_tan;
                let ndc_y = (1.0 - 2.0 * (f64::from(py) + jy) / h) * scene.camera.half_tan;
                legacy_camera_ray(&scene.camera, ndc_x, ndc_y)
            }
            CameraPath::Cinematic { camera, exposure } => {
                let time = ray_time.ok_or(TracerError::MissingRayTime)?;
                let physical = camera.evaluate_exposure(
                    cx,
                    exposure,
                    time.interval.time_at(time.normalized),
                )?;
                let half_tan = physical.projection().vertical_half_tan();
                let x_tan = (2.0 * (f64::from(px) + jx) / w - 1.0) * s.camera_aspect() * half_tan;
                let y_tan = (1.0 - 2.0 * (f64::from(py) + jy) / h) * half_tan;
                physical.generate_ray_from_tangent_offsets(
                    cx,
                    x_tan,
                    y_tan,
                    camera_lens_sample(key, pixel, sample)?,
                )?
            }
        };
        SpectralPathState {
            ray,
            throughput: [1.0; PACKET],
            previous_bsdf: None,
            prev_origin: ray.origin,
            segment_origin: ray.origin,
            medium_stack: MediumStack::new(),
            rng,
            next_depth: 0,
            active_lane: None,
        }
    };
    let mut radiance = [0.0f64; PACKET];
    let mut contribution_radiance = capture_primary.then_some(PathContributionRadiance::ZERO);
    let mut resumed_xyz = [0.0; 3];
    let mut resumed_contribution_split = capture_primary.then_some(PathContributionSplit::ZERO);
    let mut primary = None;
    let mut ray = state.ray;
    let mut throughput = state.throughput;
    let mut previous_bsdf = state.previous_bsdf;
    let mut prev_origin = state.prev_origin;
    let mut segment_origin = state.segment_origin;
    let mut medium_stack = state.medium_stack;
    let mut rng = state.rng;
    let active_lane = state.active_lane;
    for depth in state.next_depth..s.max_depth {
        cx.checkpoint()?;
        let Some(intersection) = intersect(scene, cx, &ray, ray_time)? else {
            if let Some(active) = medium_stack.last() {
                return Err(unclosed_medium_error(
                    active.boundary_primitive,
                    "path_miss",
                    pixel,
                    sample,
                    depth,
                    active_lane,
                    &ray,
                    None,
                ));
            }
            let completed_slab = previous_bsdf.and_then(|previous| match previous.smooth_slab {
                Some(SmoothSlabPath::Exited(slab)) => Some(slab),
                Some(SmoothSlabPath::Entered(_)) | None => None,
            });
            let environment_origin = completed_slab.map_or(prev_origin, |slab| slab.source_origin);
            if let Some(environment) = lighting.environment_evaluation(environment_origin, ray.dir)
            {
                let ordinary_weight = completed_slab.is_none().then(|| {
                    emissive_hit_weight(
                        s.strategy,
                        previous_bsdf,
                        previous_bsdf.map(|_| environment.pdf_solid_angle),
                    )
                });
                let (spectrum, scale) = environment.emission;
                for (lane, &lambda) in lambdas.iter().enumerate() {
                    let weight = if let Some(slab) = completed_slab {
                        let nee_pdf = completed_slab_environment_nee_pdf(
                            scene,
                            cx,
                            ray_time,
                            slab,
                            ray.dir,
                            environment.pdf_solid_angle,
                            lambda,
                        )?;
                        completed_slab_mis_weight(s.strategy, slab, nee_pdf)
                    } else {
                        ordinary_weight.expect("ordinary environment MIS weight")
                    };
                    let before = radiance[lane];
                    radiance[lane] += throughput[lane] * spectrum.eval(lambda) * scale * weight;
                    if let Some(contributions) = &mut contribution_radiance {
                        contributions.record(
                            emissive_contribution_class(depth),
                            lane,
                            radiance[lane] - before,
                        );
                    }
                }
            }
            break;
        };
        let SceneIntersection {
            primitive_index: prim_idx,
            hit,
            instance_hit,
        } = intersection;
        attenuate_segment(
            &mut throughput,
            &lambdas,
            &medium_stack,
            segment_origin,
            &hit,
        )?;
        let prim = &scene.primitives[prim_idx];
        let frame = surface_frame(&hit, &ray)?;
        if capture_primary && depth == 0 {
            let material_identity = prim.material.content_identity();
            let surface = instance_hit
                .as_ref()
                .map(|instance_hit| {
                    PrimarySurfaceSample::try_from_instance_hit(instance_hit, material_identity)
                })
                .transpose()?;
            primary = Some(PrimaryTraceHit {
                primitive_index: prim_idx,
                hit,
                beauty_shading_normal_world: frame.oriented,
                material_identity,
                surface,
            });
        }
        let n = frame.oriented;
        if let Some((spec, scale)) = &prim.emission {
            // MIS weight against NEE for this light, seen from the
            // previous vertex.
            let light_index = lighting.rect_index_for_primitive(prim_idx);
            let completed_slab = previous_bsdf.and_then(|previous| match previous.smooth_slab {
                Some(SmoothSlabPath::Exited(slab)) => Some(slab),
                Some(SmoothSlabPath::Entered(_)) | None => None,
            });
            let ordinary_nee_pdf = if completed_slab.is_none()
                && s.strategy == DirectStrategy::Mis
                && previous_bsdf.is_some()
            {
                light_index.map(|light_index| {
                    lighting.rect_mixture_pdf(light_index, prev_origin, hit.point)
                })
            } else {
                None
            };
            for (k, &l) in lambdas.iter().enumerate() {
                if completed_slab.is_some() && throughput[k].to_bits() == 0.0_f64.to_bits() {
                    continue;
                }
                let weight = if let (Some(slab), Some(light_index)) = (completed_slab, light_index)
                {
                    let nee_pdf = completed_slab_rectangle_nee_pdf(
                        scene,
                        lighting,
                        cx,
                        ray_time,
                        slab,
                        light_index,
                        hit.point,
                        l,
                    )?;
                    completed_slab_mis_weight(s.strategy, slab, nee_pdf)
                } else {
                    emissive_hit_weight(s.strategy, previous_bsdf, ordinary_nee_pdf)
                };
                let before = radiance[k];
                radiance[k] += throughput[k] * spec.eval(l) * scale * weight;
                if let Some(contributions) = &mut contribution_radiance {
                    contributions.record(
                        emissive_contribution_class(depth),
                        k,
                        radiance[k] - before,
                    );
                }
            }
            break; // v1: emitters do not reflect
        }

        let dielectric_boundary = match prim.material {
            Material::Dielectric { glass, .. } => Some(boundary_media(
                prim_idx,
                glass,
                frame.entering,
                &medium_stack,
            )?),
            Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => None,
        };

        // Next-event estimation.
        if s.strategy != DirectStrategy::BsdfOnly {
            match prim.material {
                Material::Dielectric { surface, .. } => {
                    if let Some(alpha) = surface.roughness_alpha() {
                        let boundary = dielectric_boundary.ok_or(TracerError::InvalidInput)?;
                        let (u1, u2) = rng.next2();
                        if let Some(direct) = lighting
                            .sample(hit.point, u1, u2)
                            .and_then(|sample| prepare_direct_light(hit.point, sample))
                        {
                            let wi = direct.direction;
                            let cos_s = n.dot(wi).abs();
                            let wo = ray.dir.scale(-1.0);
                            if cos_s > 0.0 {
                                let shadow = Ray {
                                    origin: dielectric_spawn_origin(hit.point, frame.geometric, wi),
                                    dir: wi,
                                };
                                let shadow_hit = intersect(scene, cx, &shadow, ray_time)?;
                                let visible = match (direct.target, shadow_hit) {
                                    (
                                        DirectLightTarget::Rectangle {
                                            primitive_index,
                                            distance_m,
                                            ..
                                        },
                                        Some(shadow_hit),
                                    ) => {
                                        shadow_hit.primitive_index == primitive_index
                                            && shadow_hit.hit.t > distance_m - 1.0e-4
                                    }
                                    (DirectLightTarget::Environment, None) => true,
                                    (DirectLightTarget::Rectangle { .. }, None)
                                    | (DirectLightTarget::Environment, Some(_)) => false,
                                };
                                if visible {
                                    let pdf_nee = direct.pdf_solid_angle;
                                    let shadow_medium = if n.dot(wi) > 0.0 {
                                        boundary.incident
                                    } else {
                                        boundary.transmitted
                                    };
                                    if matches!(direct.target, DirectLightTarget::Environment)
                                        && shadow_medium.is_some()
                                    {
                                        return Err(unclosed_medium_error(
                                            medium_stack.last().map_or(prim_idx, |active| {
                                                active.boundary_primitive
                                            }),
                                            "dielectric_nee_environment",
                                            pixel,
                                            sample,
                                            depth,
                                            active_lane,
                                            &shadow,
                                            Some((frame.geometric, ray.dir, n.dot(wi))),
                                        ));
                                    }
                                    let (emission, emission_scale) = direct.emission;
                                    for (lane, &lambda) in lambdas.iter().enumerate() {
                                        if throughput[lane].to_bits() == 0.0_f64.to_bits() {
                                            continue;
                                        }
                                        let eta_i = medium_ior(boundary.incident, lambda)?;
                                        let eta_t = medium_ior(boundary.transmitted, lambda)?;
                                        let evaluation = evaluate_rough_dielectric(
                                            n, wo, wi, eta_i, eta_t, alpha,
                                        )?;
                                        if evaluation.value <= 0.0 {
                                            continue;
                                        }
                                        // The first dispersive boundary fans out before
                                        // event selection, so every active wavelength is
                                        // sampled from its own complete BSDF. The competing
                                        // MIS density is therefore this lane's native PDF;
                                        // there is no hero-lane proposal ratio here.
                                        let competing_pdf = evaluation.pdf;
                                        let weight = match s.strategy {
                                            DirectStrategy::Mis
                                                if depth + 1 == s.max_depth
                                                    && !lighting.is_legacy_compatibility_path() =>
                                            {
                                                1.0
                                            }
                                            DirectStrategy::Mis => {
                                                balance_heuristic(1, pdf_nee, 1, competing_pdf)
                                            }
                                            DirectStrategy::NeeOnly => 1.0,
                                            DirectStrategy::BsdfOnly => {
                                                return Err(TracerError::InvalidInput);
                                            }
                                        };
                                        let attenuation = match direct.target {
                                            DirectLightTarget::Rectangle { distance_m, .. } => {
                                                medium_transmittance(
                                                    shadow_medium,
                                                    lambda,
                                                    distance_m,
                                                )?
                                            }
                                            DirectLightTarget::Environment => 1.0,
                                        };
                                        let before = radiance[lane];
                                        radiance[lane] += throughput[lane]
                                            * evaluation.value
                                            * cos_s
                                            * attenuation
                                            * emission.eval(lambda)
                                            * emission_scale
                                            / pdf_nee
                                            * weight;
                                        if let Some(contributions) = &mut contribution_radiance {
                                            contributions.record(
                                                direct_contribution_class(depth),
                                                lane,
                                                radiance[lane] - before,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
                    // Preserve the opaque tracer-v1 arithmetic and draw order
                    // expression-for-expression so the frozen Cornell path is
                    // unaffected by enabling dielectric support.
                    let (u1, u2) = rng.next2();
                    if let Some(direct) = lighting
                        .sample(hit.point, u1, u2)
                        .and_then(|sample| prepare_direct_light(hit.point, sample))
                    {
                        let wi = direct.direction;
                        let cos_s = n.dot(wi);
                        if cos_s > 0.0 {
                            let shadow = Ray {
                                origin: hit.point.offset(n.scale(RAY_EPS)),
                                dir: wi,
                            };
                            let shadow_hit = intersect(scene, cx, &shadow, ray_time)?;
                            let visible = match (direct.target, shadow_hit) {
                                (
                                    DirectLightTarget::Rectangle {
                                        primitive_index,
                                        distance_m,
                                        ..
                                    },
                                    Some(shadow_hit),
                                ) => {
                                    shadow_hit.primitive_index == primitive_index
                                        && shadow_hit.hit.t > distance_m - 1.0e-4
                                }
                                (DirectLightTarget::Environment, None) => true,
                                (DirectLightTarget::Rectangle { .. }, None)
                                | (DirectLightTarget::Environment, Some(_)) => false,
                            };
                            if visible {
                                if matches!(direct.target, DirectLightTarget::Environment)
                                    && let Some(active) = medium_stack.last()
                                {
                                    return Err(unclosed_medium_error(
                                        active.boundary_primitive,
                                        "opaque_nee_environment",
                                        pixel,
                                        sample,
                                        depth,
                                        active_lane,
                                        &shadow,
                                        Some((frame.geometric, ray.dir, n.dot(wi))),
                                    ));
                                }
                                let pdf_nee = direct.pdf_solid_angle;
                                let wo = ray.dir.scale(-1.0);
                                let bsdf_pdf = bsdf_pdf(&prim.material, n, wo, wi);
                                let weight = match s.strategy {
                                    DirectStrategy::Mis
                                        if depth + 1 == s.max_depth
                                            && !lighting.is_legacy_compatibility_path() =>
                                    {
                                        1.0
                                    }
                                    DirectStrategy::Mis => {
                                        balance_heuristic(1, pdf_nee, 1, bsdf_pdf)
                                    }
                                    _ => 1.0,
                                };
                                let (espec, escale) = direct.emission;
                                if let Some(active) = medium_stack.last() {
                                    if let DirectLightTarget::Rectangle { distance_m, .. } =
                                        direct.target
                                    {
                                        for (k, &l) in lambdas.iter().enumerate() {
                                            let f = opaque_bsdf_eval(
                                                &prim.material,
                                                n,
                                                wo,
                                                wi,
                                                l,
                                                Some(active.glass),
                                            )?;
                                            let attenuation = medium_transmittance(
                                                Some(active.glass),
                                                l,
                                                distance_m,
                                            )?;
                                            let before = radiance[k];
                                            radiance[k] += throughput[k]
                                                * f
                                                * cos_s
                                                * attenuation
                                                * espec.eval(l)
                                                * escale
                                                / pdf_nee
                                                * weight;
                                            if let Some(contributions) = &mut contribution_radiance
                                            {
                                                contributions.record(
                                                    direct_contribution_class(depth),
                                                    k,
                                                    radiance[k] - before,
                                                );
                                            }
                                        }
                                    }
                                } else {
                                    // Keep the ambient opaque tracer-v1 arithmetic
                                    // bit-for-bit identical to its frozen path.
                                    for (k, &l) in lambdas.iter().enumerate() {
                                        let f =
                                            opaque_bsdf_eval(&prim.material, n, wo, wi, l, None)?;
                                        let before = radiance[k];
                                        radiance[k] +=
                                            throughput[k] * f * cos_s * espec.eval(l) * escale
                                                / pdf_nee
                                                * weight;
                                        if let Some(contributions) = &mut contribution_radiance {
                                            contributions.record(
                                                direct_contribution_class(depth),
                                                k,
                                                radiance[k] - before,
                                            );
                                        }
                                    }
                                }
                            } else if let Some(blocker) = shadow_hit
                                && primitive_is_dielectric(scene, blocker.primitive_index)
                            {
                                'optional_slab_nee: {
                                    // The current slab connector models one isolated
                                    // dielectric layer in ambient. If the source is
                                    // already inside a medium, this would instead be a
                                    // nested-media connection, so this optional NEE
                                    // technique has no proposal. Ordinary BSDF transport
                                    // remains support-complete for the path.
                                    if medium_stack.len() != 0 {
                                        break 'optional_slab_nee;
                                    }
                                    let slab = match discover_parallel_slab(
                                        scene, cx, &shadow, blocker, ray_time,
                                    ) {
                                        Ok(slab) => slab,
                                        // A bevel, sidewall, rough boundary, or other
                                        // unsupported dielectric face simply has no
                                        // slab-NEE proposal. Ordinary BSDF transport
                                        // remains responsible for paths through it.
                                        Err(TracerError::UnsupportedSlabNee { .. }) => {
                                            break 'optional_slab_nee;
                                        }
                                        Err(error) => return Err(error),
                                    };
                                    let wo = ray.dir.scale(-1.0);
                                    let (espec, escale) = direct.emission;
                                    let competing_bsdf_path_is_evaluated = depth
                                        .checked_add(3)
                                        .is_some_and(|emitter_depth| emitter_depth < s.max_depth);
                                    for (lane, &lambda) in lambdas.iter().enumerate() {
                                        if throughput[lane].to_bits() == 0.0_f64.to_bits() {
                                            continue;
                                        }
                                        let connection = match trace_parallel_slab_direct_lane(
                                            scene,
                                            cx,
                                            ray_time,
                                            hit.point,
                                            frame.geometric,
                                            slab,
                                            direct,
                                            lambda,
                                        ) {
                                            Ok(connection) => connection,
                                            // The exact refracted connection can leave
                                            // the admitted broad face even when the
                                            // straight visibility probe entered it.
                                            // That wavelength has no slab-NEE proposal;
                                            // the BSDF random walk remains available.
                                            Err(TracerError::UnsupportedSlabNee { .. }) => continue,
                                            Err(error) => return Err(error),
                                        };
                                        let cos_connected = n.dot(connection.incident_direction);
                                        if !connection.visible
                                            || connection.radiance_transport == 0.0
                                            || cos_connected <= 0.0
                                        {
                                            continue;
                                        }
                                        let source_bsdf_pdf = bsdf_pdf(
                                            &prim.material,
                                            n,
                                            wo,
                                            connection.incident_direction,
                                        );
                                        let weight = slab_nee_mis_weight(
                                            s.strategy,
                                            connection.nee_pdf_solid_angle,
                                            source_bsdf_pdf,
                                            connection.transmission_probability,
                                            competing_bsdf_path_is_evaluated,
                                        )?;
                                        let f = opaque_bsdf_eval(
                                            &prim.material,
                                            n,
                                            wo,
                                            connection.incident_direction,
                                            lambda,
                                            None,
                                        )?;
                                        let before = radiance[lane];
                                        radiance[lane] += throughput[lane]
                                            * f
                                            * cos_connected
                                            * connection.radiance_transport
                                            * espec.eval(lambda)
                                            * escale
                                            / connection.nee_pdf_solid_angle
                                            * weight;
                                        if let Some(contributions) = &mut contribution_radiance {
                                            contributions.record(
                                                direct_contribution_class(depth),
                                                lane,
                                                radiance[lane] - before,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // BSDF sampling for the next bounce.
        let (u1, u2) = rng.next2();
        let wo = ray.dir.scale(-1.0);
        match prim.material {
            Material::Dielectric { glass, surface } => {
                let boundary = dielectric_boundary.ok_or(TracerError::InvalidInput)?;
                let proposal_lane = active_lane.unwrap_or(0);
                let event_sample = if surface.is_delta() {
                    u1
                } else {
                    rng.next2().0
                };
                if should_split_dispersive_boundary(active_lane, &boundary, &lambdas)? {
                    let continuations = sample_dispersive_dielectric_lanes(
                        n,
                        wo,
                        surface,
                        &boundary,
                        &lambdas,
                        u1,
                        u2,
                        event_sample,
                    )?;
                    // Split before conditioning on any one wavelength's event.
                    // This is support-complete at exit interfaces where a blue
                    // lane can undergo TIR while a red lane still transmits.
                    // Each suffix runs to completion before the next fixed lane,
                    // so at most two states are live. `active_lane` prevents a
                    // child from splitting again.
                    for (lane, continuation) in continuations.into_iter().enumerate() {
                        let Some(continuation) = continuation else {
                            continue;
                        };
                        if throughput[lane].to_bits() == 0.0_f64.to_bits()
                            || continuation.weight.to_bits() == 0.0_f64.to_bits()
                        {
                            continue;
                        }
                        let child_previous_bsdf = previous_after_dielectric_sample(
                            scene,
                            cx,
                            ray_time,
                            previous_bsdf,
                            prev_origin,
                            &medium_stack,
                            prim_idx,
                            glass,
                            surface,
                            frame,
                            hit,
                            instance_hit,
                            ray,
                            continuation.direction,
                            continuation.event,
                            continuation.pdf,
                            continuation.delta,
                            lambdas[lane],
                        )?;
                        let mut child_medium_stack = medium_stack;
                        if continuation.event == DielectricEvent::Transmission {
                            apply_medium_transition(&mut child_medium_stack, boundary.transition)?;
                        }
                        let mut child_throughput = [0.0; PACKET];
                        child_throughput[lane] = throughput[lane] * continuation.weight;
                        let child = trace_path(
                            scene,
                            lighting,
                            cx,
                            s,
                            kn,
                            pixel,
                            sample,
                            jx,
                            jy,
                            ul,
                            ray_time,
                            camera_path,
                            capture_primary,
                            Some(SpectralPathState {
                                ray: Ray {
                                    origin: dielectric_spawn_origin(
                                        hit.point,
                                        frame.geometric,
                                        continuation.direction,
                                    ),
                                    dir: continuation.direction,
                                },
                                throughput: child_throughput,
                                previous_bsdf: Some(child_previous_bsdf),
                                prev_origin: hit.point,
                                segment_origin: hit.point,
                                medium_stack: child_medium_stack,
                                rng,
                                next_depth: depth + 1,
                                active_lane: Some(lane),
                            }),
                        )?;
                        add_xyz(&mut resumed_xyz, child.xyz);
                        if let Some(total) = &mut resumed_contribution_split {
                            total.add_assign(
                                child.contribution_split.ok_or(TracerError::InvalidInput)?,
                            );
                        }
                    }
                    break;
                }
                let Some(sampled) = sample_dielectric_path(
                    n,
                    wo,
                    surface,
                    &boundary,
                    &lambdas,
                    proposal_lane,
                    u1,
                    u2,
                    event_sample,
                )?
                else {
                    break;
                };
                for (lane, weight) in sampled.weights.into_iter().enumerate() {
                    throughput[lane] *= weight;
                }
                let next_previous_bsdf = previous_after_dielectric_sample(
                    scene,
                    cx,
                    ray_time,
                    previous_bsdf,
                    prev_origin,
                    &medium_stack,
                    prim_idx,
                    glass,
                    surface,
                    frame,
                    hit,
                    instance_hit,
                    ray,
                    sampled.direction,
                    sampled.event,
                    sampled.pdf,
                    sampled.delta,
                    lambdas[proposal_lane],
                )?;
                if sampled.event == DielectricEvent::Transmission {
                    apply_medium_transition(&mut medium_stack, boundary.transition)?;
                }
                previous_bsdf = Some(next_previous_bsdf);
                prev_origin = hit.point;
                segment_origin = hit.point;
                ray = Ray {
                    origin: dielectric_spawn_origin(hit.point, frame.geometric, sampled.direction),
                    dir: sampled.direction,
                };
            }
            Material::Lambertian { .. } | Material::Ggx { .. } | Material::Conductor { .. } => {
                let Some((wi, pdf)) = bsdf_sample(&prim.material, n, wo, u1, u2) else {
                    break;
                };
                let cos_s = n.dot(wi).max(0.0);
                if pdf <= 0.0 || cos_s <= 0.0 {
                    break;
                }
                for (k, &l) in lambdas.iter().enumerate() {
                    throughput[k] *= opaque_bsdf_eval(
                        &prim.material,
                        n,
                        wo,
                        wi,
                        l,
                        medium_stack.last().map(|active| active.glass),
                    )? * cos_s
                        / pdf;
                }
                previous_bsdf = Some(PreviousBsdf {
                    pdf,
                    delta: false,
                    opaque_source_geometric_normal: Some(n),
                    smooth_slab: None,
                });
                prev_origin = hit.point;
                segment_origin = hit.point;
                ray = Ray {
                    origin: hit.point.offset(n.scale(RAY_EPS)),
                    dir: wi,
                };
            }
        }
    }
    // Hero-wavelength estimator → XYZ (same normalization convention
    // as `spectral::xyz_of_spectrum`: Y of unit radiance is 1).
    let range = LAMBDA_MAX - LAMBDA_MIN;
    let mut xyz = [0.0f64; 3];
    for (k, &l) in lambdas.iter().enumerate() {
        let w = radiance[k] * range / PACKET as f64 * kn;
        xyz[0] += w * cie_x(l);
        xyz[1] += w * cie_y(l);
        xyz[2] += w * cie_z(l);
    }
    add_xyz(&mut xyz, resumed_xyz);
    let contribution_split = contribution_radiance.map(|contributions| {
        let mut split = PathContributionSplit {
            direct_xyz: packet_radiance_to_xyz(contributions.direct, &lambdas, kn),
            indirect_xyz: packet_radiance_to_xyz(contributions.indirect, &lambdas, kn),
            emission_xyz: packet_radiance_to_xyz(contributions.emission, &lambdas, kn),
        };
        split
            .add_assign(resumed_contribution_split.expect("captured recursive contribution split"));
        split
    });
    Ok(PathTraceSample {
        xyz,
        primary,
        absolute_time_s: ray_time.map(|time| time.interval.time_at(time.normalized)),
        pixel_jitter: [jx, jy],
        contribution_split,
    })
}

const fn emissive_contribution_class(depth: u32) -> PathContributionClass {
    match depth {
        0 => PathContributionClass::Emission,
        1 => PathContributionClass::Direct,
        _ => PathContributionClass::Indirect,
    }
}

const fn direct_contribution_class(depth: u32) -> PathContributionClass {
    if depth == 0 {
        PathContributionClass::Direct
    } else {
        PathContributionClass::Indirect
    }
}

fn packet_radiance_to_xyz(radiance: [f64; PACKET], lambdas: &[f64], kn: f64) -> [f64; 3] {
    let range = LAMBDA_MAX - LAMBDA_MIN;
    let mut xyz = [0.0; 3];
    for (k, &lambda) in lambdas.iter().enumerate() {
        let weight = radiance[k] * range / PACKET as f64 * kn;
        xyz[0] += weight * cie_x(lambda);
        xyz[1] += weight * cie_y(lambda);
        xyz[2] += weight * cie_z(lambda);
    }
    xyz
}

/// Preserve the exact v1 camera arithmetic as one explicit compatibility
/// branch. In particular, a zero-aperture cinematic camera must call this
/// helper rather than construct and subtract a focus point: those expressions
/// are mathematically equivalent but are not generally bit-equivalent.
fn legacy_camera_ray(camera: &Camera, ndc_x: f64, ndc_y: f64) -> Ray {
    let right = cross(camera.forward, camera.up);
    let dir = unit(Vec3::new(
        camera.forward.x + ndc_x * right.x + ndc_y * camera.up.x,
        camera.forward.y + ndc_x * right.y + ndc_y * camera.up.y,
        camera.forward.z + ndc_x * right.z + ndc_y * camera.up.z,
    ));
    Ray {
        origin: camera.eye,
        dir,
    }
}

fn camera_lens_sample(key: [u32; 2], pixel: u32, sample: u32) -> Result<LensSample, TracerError> {
    let u = philox4x32_10([pixel, sample, CAMERA_LENS_SAMPLE_DOMAIN_V1, 1], key);
    let v = philox4x32_10([pixel, sample, CAMERA_LENS_SAMPLE_DOMAIN_V1, 2], key);
    Ok(LensSample::try_new(u32_unit(u[0]), u32_unit(v[0]))?)
}

fn surface_frame(hit: &Hit, ray: &Ray) -> Result<SurfaceFrame, TracerError> {
    let geometric = hit.normal.ok_or(TracerError::MissingNormal)?;
    if !geometric.x.is_finite()
        || !geometric.y.is_finite()
        || !geometric.z.is_finite()
        || geometric.norm() <= 0.0
    {
        return Err(TracerError::MissingNormal);
    }
    let norm = geometric.norm();
    let geometric_unit = geometric.scale(1.0 / norm);
    let entering = geometric.dot(ray.dir) <= 0.0;
    let oriented = if geometric.dot(ray.dir) > 0.0 {
        geometric.scale(-1.0)
    } else {
        geometric
    };
    Ok(SurfaceFrame {
        oriented,
        geometric: geometric_unit,
        entering,
    })
}

fn boundary_media(
    boundary_primitive: usize,
    glass: DielectricGlass,
    entering: bool,
    stack: &MediumStack,
) -> Result<BoundaryMedia, TracerError> {
    if entering {
        if stack.len() >= MAX_MEDIUM_STACK_DEPTH {
            return Err(TracerError::MediumStackOverflow);
        }
        if stack
            .iter()
            .any(|entry| entry.boundary_primitive == boundary_primitive)
        {
            return Err(TracerError::MediumStackMismatch {
                boundary_primitive,
                active_boundary: stack.last().map(|entry| entry.boundary_primitive),
            });
        }
        Ok(BoundaryMedia {
            incident: stack.last().map(|entry| entry.glass),
            transmitted: Some(glass),
            transition: MediumTransition::Enter(MediumEntry {
                boundary_primitive,
                glass,
            }),
        })
    } else {
        let Some(active) = stack.last() else {
            return Err(TracerError::MediumStackMismatch {
                boundary_primitive,
                active_boundary: None,
            });
        };
        if active.boundary_primitive != boundary_primitive || active.glass != glass {
            return Err(TracerError::MediumStackMismatch {
                boundary_primitive,
                active_boundary: Some(active.boundary_primitive),
            });
        }
        Ok(BoundaryMedia {
            incident: Some(active.glass),
            transmitted: stack
                .len()
                .checked_sub(2)
                .and_then(|index| stack.get(index))
                .map(|entry| entry.glass),
            transition: MediumTransition::Exit { boundary_primitive },
        })
    }
}

fn apply_medium_transition(
    stack: &mut MediumStack,
    transition: MediumTransition,
) -> Result<(), TracerError> {
    match transition {
        MediumTransition::Enter(entry) => {
            if stack.len() >= MAX_MEDIUM_STACK_DEPTH
                || stack
                    .iter()
                    .any(|active| active.boundary_primitive == entry.boundary_primitive)
            {
                return Err(if stack.len() >= MAX_MEDIUM_STACK_DEPTH {
                    TracerError::MediumStackOverflow
                } else {
                    TracerError::MediumStackMismatch {
                        boundary_primitive: entry.boundary_primitive,
                        active_boundary: stack.last().map(|active| active.boundary_primitive),
                    }
                });
            }
            stack.push(entry)?;
        }
        MediumTransition::Exit { boundary_primitive } => {
            if stack.last().map(|entry| entry.boundary_primitive) != Some(boundary_primitive) {
                return Err(TracerError::MediumStackMismatch {
                    boundary_primitive,
                    active_boundary: stack.last().map(|entry| entry.boundary_primitive),
                });
            }
            stack.pop();
        }
    }
    Ok(())
}

fn medium_ior(medium: Option<DielectricGlass>, wavelength_nm: f64) -> Result<f64, TracerError> {
    medium.map_or(Ok(1.0), |glass| {
        glass.ior().eval(wavelength_nm).map_err(Into::into)
    })
}

fn medium_transmittance(
    medium: Option<DielectricGlass>,
    wavelength_nm: f64,
    distance_m: f64,
) -> Result<f64, TracerError> {
    medium.map_or(Ok(1.0), |glass| {
        glass
            .absorption()
            .transmittance(wavelength_nm, distance_m)
            .map_err(Into::into)
    })
}

fn attenuate_segment(
    throughput: &mut [f64; PACKET],
    lambdas: &[f64],
    medium_stack: &MediumStack,
    physical_origin: Point3,
    hit: &Hit,
) -> Result<(), TracerError> {
    let Some(active) = medium_stack.last() else {
        return Ok(());
    };
    // Intersection rays are numerically offset from the preceding boundary;
    // absorption is measured from the unshifted physical vertex so the offset
    // cannot silently shorten every path through a thin medium.
    let distance_m = hit.point.delta_from(physical_origin).norm();
    if !distance_m.is_finite() || distance_m < 0.0 {
        return Err(TracerError::Dielectric(DielectricError::InvalidAbsorption));
    }
    for (lane, &lambda) in lambdas.iter().enumerate() {
        if throughput[lane].to_bits() != 0.0_f64.to_bits() {
            throughput[lane] *= active
                .glass
                .absorption()
                .transmittance(lambda, distance_m)?;
        }
    }
    Ok(())
}

fn boundary_direction_is_wavelength_dependent(
    boundary: &BoundaryMedia,
    lambdas: &[f64],
) -> Result<bool, TracerError> {
    let Some((&first, remaining)) = lambdas.split_first() else {
        return Err(TracerError::InvalidInput);
    };
    let first_ratio =
        medium_ior(boundary.incident, first)? / medium_ior(boundary.transmitted, first)?;
    for &lambda in remaining {
        let ratio =
            medium_ior(boundary.incident, lambda)? / medium_ior(boundary.transmitted, lambda)?;
        if ratio.to_bits() != first_ratio.to_bits() {
            return Ok(true);
        }
    }
    Ok(false)
}

fn should_split_dispersive_boundary(
    active_lane: Option<usize>,
    boundary: &BoundaryMedia,
    lambdas: &[f64],
) -> Result<bool, TracerError> {
    if active_lane.is_some() {
        return Ok(false);
    }
    boundary_direction_is_wavelength_dependent(boundary, lambdas)
}

#[cfg(test)]
const fn maximum_spectral_traversals(max_depth: u32) -> u64 {
    max_depth as u64 * PACKET as u64
}

/// Sample one complete dielectric continuation for every wavelength at the
/// first dispersive boundary. All lanes share the admitted microfacet and
/// event uniforms as a variance-reducing correlation, but each marginal uses
/// its own Fresnel event probability, Snell direction, throughput, and PDF.
/// Splitting before event conditioning is essential: at an exit boundary one
/// wavelength may undergo TIR while another still has nonzero BTDF support.
#[allow(clippy::too_many_arguments)]
fn sample_dispersive_dielectric_lanes(
    normal: Vec3,
    wo: Vec3,
    surface: DielectricSurface,
    boundary: &BoundaryMedia,
    lambdas: &[f64],
    microfacet_u: f64,
    microfacet_v: f64,
    event_sample: f64,
) -> Result<[Option<DielectricLaneContinuation>; PACKET], TracerError> {
    if lambdas.len() != PACKET || !event_sample.is_finite() || !(0.0..1.0).contains(&event_sample) {
        return Err(TracerError::InvalidInput);
    }
    let mut continuations = [None; PACKET];
    for (lane, &lambda) in lambdas.iter().enumerate() {
        let eta_i = medium_ior(boundary.incident, lambda)?;
        let eta_t = medium_ior(boundary.transmitted, lambda)?;
        let continuation = if let Some(alpha) = surface.roughness_alpha() {
            sample_rough_dielectric(
                normal,
                wo,
                eta_i,
                eta_t,
                alpha,
                microfacet_u,
                microfacet_v,
                event_sample,
            )?
            .map(|sample| DielectricLaneContinuation {
                direction: sample.direction,
                event: sample.event,
                weight: sample.radiance_weight,
                pdf: if sample.delta { 0.0 } else { sample.pdf },
                delta: sample.delta,
            })
        } else {
            let sample = sample_smooth_dielectric(normal, wo, eta_i, eta_t, event_sample)?;
            Some(DielectricLaneContinuation {
                direction: sample.direction,
                event: sample.event,
                weight: sample.radiance_weight,
                pdf: 0.0,
                delta: true,
            })
        };
        if let Some(continuation) = continuation {
            if !continuation.weight.is_finite()
                || continuation.weight < 0.0
                || !continuation.pdf.is_finite()
                || continuation.pdf < 0.0
            {
                return Err(TracerError::InvalidInput);
            }
        }
        continuations[lane] = continuation;
    }
    Ok(continuations)
}

#[allow(clippy::too_many_arguments)]
fn sample_dielectric_path(
    normal: Vec3,
    wo: Vec3,
    surface: DielectricSurface,
    boundary: &BoundaryMedia,
    lambdas: &[f64],
    proposal_lane: usize,
    microfacet_u: f64,
    microfacet_v: f64,
    event_sample: f64,
) -> Result<Option<DielectricPathSample>, TracerError> {
    let proposal_lambda = *lambdas
        .get(proposal_lane)
        .ok_or(TracerError::InvalidInput)?;
    let proposal_eta_i = medium_ior(boundary.incident, proposal_lambda)?;
    let proposal_eta_t = medium_ior(boundary.transmitted, proposal_lambda)?;
    if let Some(alpha) = surface.roughness_alpha() {
        let Some(sample) = sample_rough_dielectric(
            normal,
            wo,
            proposal_eta_i,
            proposal_eta_t,
            alpha,
            microfacet_u,
            microfacet_v,
            event_sample,
        )?
        else {
            return Ok(None);
        };
        let mut weights = [0.0; PACKET];
        if sample.delta {
            weights.fill(sample.radiance_weight);
        } else {
            let sampled_cosine = admitted_unit_dot(normal, sample.direction)
                .abs()
                .clamp(0.0, 1.0);
            for (lane, &lambda) in lambdas.iter().enumerate() {
                let eta_i = medium_ior(boundary.incident, lambda)?;
                let eta_t = medium_ior(boundary.transmitted, lambda)?;
                let evaluation =
                    evaluate_rough_dielectric(normal, wo, sample.direction, eta_i, eta_t, alpha)?;
                weights[lane] = evaluation.value * sampled_cosine / sample.pdf;
            }
        }
        return Ok(Some(DielectricPathSample {
            direction: sample.direction,
            event: sample.event,
            pdf: if sample.delta { 0.0 } else { sample.pdf },
            delta: sample.delta,
            weights,
        }));
    }

    let sample =
        sample_smooth_dielectric(normal, wo, proposal_eta_i, proposal_eta_t, event_sample)?;
    let mut weights = [0.0; PACKET];
    // The dielectric sampler has already admitted and normalized this frame.
    // Clamp the independently recomputed cosine so a final binary64 ulp above
    // one cannot make companion-wavelength Fresnel evaluation refuse.
    let incident_cosine = admitted_unit_dot(normal, wo).clamp(0.0, 1.0);
    for (lane, &lambda) in lambdas.iter().enumerate() {
        let eta_i = medium_ior(boundary.incident, lambda)?;
        let eta_t = medium_ior(boundary.transmitted, lambda)?;
        let fresnel = fresnel_dielectric(incident_cosine, eta_i, eta_t)?;
        weights[lane] = match sample.event {
            DielectricEvent::Reflection => fresnel.reflectance / sample.probability,
            DielectricEvent::Transmission => {
                let eta_ratio = eta_i / eta_t;
                (1.0 - fresnel.reflectance) * eta_ratio * eta_ratio / sample.probability
            }
        };
    }
    Ok(Some(DielectricPathSample {
        direction: sample.direction,
        event: sample.event,
        pdf: 0.0,
        delta: true,
        weights,
    }))
}

/// Reproduce the dielectric module's normalization after its sampler has
/// admitted the frame; this keeps packet-lane weights on the sampled BSDF's
/// cosine convention without adding a second public validation surface.
fn admitted_unit_dot(left: Vec3, right: Vec3) -> f64 {
    let left = left.scale(1.0 / left.norm());
    let right = right.scale(1.0 / right.norm());
    left.dot(right)
}

fn dielectric_spawn_origin(point: Point3, geometric_normal: Vec3, direction: Vec3) -> Point3 {
    let side = if geometric_normal.dot(direction) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    // Keep the metric displacement independent of world translation, then
    // step each affected coordinate outward so rounding cannot erase it.
    let offset = geometric_normal.scale(side * RAY_EPS);
    Point3::new(
        offset_component(point.x, offset.x),
        offset_component(point.y, offset.y),
        offset_component(point.z, offset.z),
    )
}

fn offset_component(coordinate: f64, offset: f64) -> f64 {
    let shifted = coordinate + offset;
    if offset > 0.0 {
        shifted.next_up()
    } else if offset < 0.0 {
        shifted.next_down()
    } else {
        coordinate
    }
}

fn emissive_hit_weight(
    strategy: DirectStrategy,
    previous_bsdf: Option<PreviousBsdf>,
    designated_light_nee_pdf: Option<f64>,
) -> f64 {
    let Some(previous) = previous_bsdf else {
        return 1.0;
    };
    if previous.delta {
        return 1.0;
    }
    match strategy {
        DirectStrategy::BsdfOnly => 1.0,
        DirectStrategy::NeeOnly => 0.0,
        DirectStrategy::Mis => designated_light_nee_pdf.map_or(1.0, |nee_pdf| {
            balance_heuristic(1, previous.pdf, 1, nee_pdf)
        }),
    }
}

impl Settings {
    fn camera_aspect(&self) -> f64 {
        f64::from(self.width) / f64::from(self.height)
    }
}

fn intersect(
    scene: &Scene,
    cx: &Cx<'_>,
    ray: &Ray,
    ray_time: Option<&PathTime>,
) -> Result<Option<SceneIntersection>, TracerError> {
    let mut best: Option<SceneIntersection> = None;
    for (i, prim) in scene.primitives.iter().enumerate() {
        cx.checkpoint()?;
        let candidate = match &prim.shape {
            Shape::Mesh(mesh) => mesh
                .intersect_with_cx(cx, ray)?
                .map(|hit| SceneIntersection {
                    primitive_index: i,
                    hit,
                    instance_hit: None,
                }),
            Shape::Chart(chart) => {
                let (hit, audit) = sphere_trace(chart.as_ref(), cx, ray, 1e4, TRACE_EPS, 1.0);
                if matches!(
                    audit.termination,
                    TraceTermination::Hit
                        | TraceTermination::ResidualLimit
                        | TraceTermination::Miss
                ) && !audit.certified
                {
                    return Err(TracerError::UncertifiedTrace);
                }
                match audit.termination {
                    TraceTermination::Cancelled => return Err(TracerError::Cancelled),
                    TraceTermination::Miss => None,
                    TraceTermination::Hit => Some(SceneIntersection {
                        primitive_index: i,
                        hit: hit.ok_or(TracerError::BackendFailure(TraceTermination::Hit))?,
                        instance_hit: None,
                    }),
                    termination => return Err(TracerError::BackendFailure(termination)),
                }
            }
            Shape::Instance(instance) => {
                instance
                    .intersect(cx, ray, 1e4, TRACE_EPS)?
                    .map(|instance_hit| SceneIntersection {
                        primitive_index: i,
                        hit: instance_hit.hit,
                        instance_hit: Some(instance_hit),
                    })
            }
            Shape::AnimatedInstance(instance) => {
                let time = ray_time.ok_or(TracerError::MissingRayTime)?;
                let instance_hit = if let Some(cached) = time
                    .cached_animated
                    .iter()
                    .flatten()
                    .find(|cached| cached.primitive_index == i)
                {
                    cached.instance.intersect(cx, ray, 1e4, TRACE_EPS)?
                } else {
                    let timed_ray = TimedRay::at_normalized(*ray, time.interval, time.normalized);
                    instance.intersect(cx, &timed_ray, 1e4, TRACE_EPS)?
                };
                instance_hit.map(|instance_hit| SceneIntersection {
                    primitive_index: i,
                    hit: instance_hit.hit,
                    instance_hit: Some(instance_hit),
                })
            }
        };
        if let Some(candidate) = candidate {
            let replace = match best.as_ref() {
                None => true,
                Some(best_hit) if candidate.hit.t < best_hit.hit.t => true,
                Some(best_hit)
                    if candidate.hit.t.total_cmp(&best_hit.hit.t) == core::cmp::Ordering::Equal =>
                {
                    instance_object_id(&prim.shape)
                        .zip(instance_object_id(
                            &scene.primitives[best_hit.primitive_index].shape,
                        ))
                        .is_some_and(|(candidate, current)| candidate < current)
                }
                Some(_) => false,
            };
            if replace {
                best = Some(candidate);
            }
        }
    }
    Ok(best)
}

fn validate_scene<'scene>(
    scene: &'scene Scene,
    cx: &Cx<'_>,
    shutter: Option<ShutterInterval>,
) -> Result<AdmittedLighting<'scene>, TracerError> {
    let lighting =
        AdmittedLighting::try_new_cancellable(cx, &scene.lights, scene.environment.as_ref())?;
    for light in &scene.lights {
        cx.checkpoint()?;
        let Some(light_primitive) = scene.primitives.get(light.prim) else {
            return Err(TracerError::LightPrimitiveMismatch {
                light_primitive: light.prim,
            });
        };
        if matches!(&light_primitive.shape, Shape::AnimatedInstance(_)) {
            return Err(TracerError::AnimatedLightUnsupported);
        }
        if !emissions_match(light_primitive.emission, Some(light.emission))
            || !light_geometry_matches(light, &light_primitive.shape)
        {
            return Err(TracerError::LightPrimitiveMismatch {
                light_primitive: light.prim,
            });
        }
    }
    cx.checkpoint()?;
    let mut object_ids = BTreeSet::new();
    for primitive in &scene.primitives {
        cx.checkpoint()?;
        let Some(object_id) = instance_object_id(&primitive.shape) else {
            continue;
        };
        if !object_ids.insert(object_id) {
            return Err(TracerError::InvalidInstance);
        }
        if let Shape::AnimatedInstance(instance) = &primitive.shape {
            let shutter = shutter.ok_or(TracerError::MissingRayTime)?;
            instance.trajectory().admit_shutter(shutter)?;
        }
    }
    cx.checkpoint()?;
    Ok(lighting)
}

#[allow(clippy::too_many_arguments)]
fn preflight_render<'a>(
    scene: &'a Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: Option<&Film>,
    from: u32,
    to: u32,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
) -> Result<(AdmittedLighting<'a>, FilmTimeMode), TracerError> {
    cx.checkpoint()?;
    let expected_len = checked_pixel_len(settings.width, settings.height)?;
    if to < from {
        return Err(TracerError::InvalidRange { from, to });
    }
    let requested_mode = requested_time_mode(shutter, settings.seed, camera_path)?;
    if let Some(film) = film {
        if (film.width, film.height) != (settings.width, settings.height)
            || film.xyz.len() != expected_len
            || film.spp_done != from
        {
            return Err(TracerError::InvalidInput);
        }
        validate_film_time_mode(film, requested_mode)?;
    } else if from != 0 {
        return Err(TracerError::InvalidInput);
    }
    let lighting = validate_scene(scene, cx, shutter)?;
    cx.checkpoint()?;
    Ok((lighting, requested_mode))
}

fn checked_pixel_len(width: u32, height: u32) -> Result<usize, TracerError> {
    let pixel_count = width.checked_mul(height).ok_or(TracerError::InvalidInput)?;
    if pixel_count == 0 {
        return Err(TracerError::InvalidInput);
    }
    usize::try_from(pixel_count).map_err(|_| TracerError::InvalidInput)
}

fn light_geometry_matches(light: &RectLight, shape: &Shape) -> bool {
    match shape {
        Shape::Mesh(mesh) => rectangle_mesh_matches(light, mesh, |point| point),
        Shape::Instance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => {
                let transform = instance.transform();
                rectangle_mesh_matches(light, mesh, |point| transform.transform_point(point))
            }
            SharedGeometry::Chart(_) => false,
        },
        Shape::Chart(_) | Shape::AnimatedInstance(_) => false,
    }
}

fn rectangle_mesh_matches(
    light: &RectLight,
    mesh: &TriMesh,
    to_world: impl Fn(Point3) -> Point3,
) -> bool {
    if mesh.vertices.len() != 4 || mesh.triangles.len() != 2 {
        return false;
    }
    let expected = [
        light.corner,
        light.corner.offset(light.edge_u),
        light.corner.offset(light.edge_u).offset(light.edge_v),
        light.corner.offset(light.edge_v),
    ];
    let edge_scale = light.edge_u.norm().max(light.edge_v.norm());
    let coordinate_scale = expected
        .iter()
        .flat_map(|point| [point.x.abs(), point.y.abs(), point.z.abs()])
        .fold(0.0_f64, f64::max);
    let tolerance =
        RECT_LIGHT_GEOMETRY_REL_TOLERANCE * edge_scale + 8.0 * f64::EPSILON * coordinate_scale;

    let mut actual_to_expected = [usize::MAX; 4];
    let mut expected_used = [false; 4];
    for (actual_index, vertex) in mesh.vertices.iter().enumerate() {
        let actual = to_world(Point3::new(vertex[0], vertex[1], vertex[2]));
        let mut matched = None;
        for (expected_index, expected_point) in expected.iter().enumerate() {
            if !expected_used[expected_index] && points_match(actual, *expected_point, tolerance) {
                if matched.is_some() {
                    return false;
                }
                matched = Some(expected_index);
            }
        }
        let Some(expected_index) = matched else {
            return false;
        };
        expected_used[expected_index] = true;
        actual_to_expected[actual_index] = expected_index;
    }

    let mut topology = [[usize::MAX; 3]; 2];
    for (triangle_index, triangle) in mesh.triangles.iter().enumerate() {
        if triangle
            .iter()
            .any(|index| *index as usize >= actual_to_expected.len())
        {
            return false;
        }
        topology[triangle_index] = triangle.map(|index| actual_to_expected[index as usize]);
        topology[triangle_index].sort_unstable();
        if topology[triangle_index][0] == topology[triangle_index][1]
            || topology[triangle_index][1] == topology[triangle_index][2]
        {
            return false;
        }
    }
    topology.sort_unstable();
    topology == [[0, 1, 2], [0, 2, 3]] || topology == [[0, 1, 3], [1, 2, 3]]
}

fn points_match(left: Point3, right: Point3, tolerance: f64) -> bool {
    [
        (left.x - right.x).abs(),
        (left.y - right.y).abs(),
        (left.z - right.z).abs(),
    ]
    .into_iter()
    .all(|difference| difference.is_finite() && difference <= tolerance)
}

fn emissions_match(
    left: Option<(LiftedSpectrum, f64)>,
    right: Option<(LiftedSpectrum, f64)>,
) -> bool {
    let canonical_bits = |value: f64| {
        if value == 0.0 {
            0.0_f64.to_bits()
        } else {
            value.to_bits()
        }
    };
    match (left, right) {
        (Some((left_spectrum, left_scale)), Some((right_spectrum, right_scale))) => {
            left_spectrum
                .c
                .into_iter()
                .zip(right_spectrum.c)
                .all(|(left, right)| canonical_bits(left) == canonical_bits(right))
                && canonical_bits(left_scale) == canonical_bits(right_scale)
        }
        (None, None) => true,
        (Some(_), None) | (None, Some(_)) => false,
    }
}

fn requested_time_mode(
    shutter: Option<ShutterInterval>,
    stream_identity: u64,
    camera_path: CameraPath<'_>,
) -> Result<FilmTimeMode, TracerError> {
    match (shutter, camera_path) {
        (None, CameraPath::Legacy) => Ok(FilmTimeMode::Static),
        (Some(shutter), CameraPath::Legacy) => Ok(FilmTimeMode::Motion {
            shutter,
            stream_identity,
        }),
        (Some(shutter), CameraPath::Cinematic { exposure, .. }) => Ok(FilmTimeMode::Cinematic {
            shutter,
            stream_identity,
            shot_id: exposure.shot_id(),
        }),
        (None, CameraPath::Cinematic { .. }) => Err(TracerError::InvalidInput),
    }
}

fn validate_film_time_mode(film: &Film, requested: FilmTimeMode) -> Result<(), TracerError> {
    match (film.spp_done, film.time_mode) {
        (0, FilmTimeMode::Uninitialized) => Ok(()),
        (
            0,
            FilmTimeMode::Static | FilmTimeMode::Motion { .. } | FilmTimeMode::Cinematic { .. },
        )
        | (1.., FilmTimeMode::Uninitialized) => Err(TracerError::InvalidInput),
        (_, accepted) if accepted == requested => Ok(()),
        _ => Err(TracerError::ProgressiveTimeModeMismatch),
    }
}

fn instance_object_id(shape: &Shape) -> Option<u64> {
    match shape {
        Shape::Instance(instance) => Some(instance.object_id()),
        Shape::AnimatedInstance(instance) => Some(instance.object_id()),
        Shape::Mesh(_) | Shape::Chart(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dielectric::{BeerLambertAbsorption, CauchyIor};
    use crate::instances::RigidTransform;
    use crate::spectral::lift_rgb;
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, StreamKey};
    use fs_geom::fixtures::SphereChart;

    fn with_test_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x534c_4142_4e45_4531,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn assert_near(actual: f64, expected: f64, tolerance: f64, context: &str) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{context}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g0_shutter_sample_previous_motion_preserves_phase_across_frame_cadence() {
        let identity = hash_domain("shutter-aligned-motion-reference-test", b"identity");
        let provenance = crate::aov::CinematicAovProvenance::try_new(
            1,
            1.0 / 24.0,
            0.0,
            2.0 / 24.0,
            identity,
            identity,
            identity,
        )
        .unwrap();
        // At 24 fps and a 15-degree back-loaded shutter, frame one covers
        // [47/576, 48/576]. Its midpoint must map to the same phase in frame
        // zero [23/576, 24/576].
        let current_midpoint_s = 47.5 / 576.0;
        let previous_s = shutter_aligned_previous_reference_time(provenance, current_midpoint_s);
        assert_near(
            previous_s,
            23.5 / 576.0,
            f64::EPSILON,
            "previous matched shutter phase",
        );

        let first_frame = crate::aov::CinematicAovProvenance::try_new(
            0,
            0.0,
            0.0,
            1.0 / 24.0,
            identity,
            identity,
            identity,
        )
        .unwrap();
        let first_midpoint_s = 23.5 / 576.0;
        let cut_previous_s = shutter_aligned_previous_reference_time(first_frame, first_midpoint_s);
        assert_eq!(cut_previous_s.to_bits(), first_midpoint_s.to_bits());
    }

    #[test]
    fn equal_ior_slab_connection_is_the_unoccluded_area_estimator_exactly() {
        let source = Point3::new(0.0, 0.0, 0.0);
        let target = Point3::new(0.4, -0.2, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let light_normal = axis;
        let connection =
            slab_connection_to_point(source, target, light_normal, axis, 0.5, 1.0).unwrap();
        let displacement = target.delta_from(source);
        let distance = displacement.norm();
        let direct = displacement.scale(1.0 / distance);
        assert!(parallel_directions(connection.incident_direction, direct));
        assert!(parallel_directions(connection.internal_direction, direct));
        assert_near(
            connection.light_area_jacobian.unwrap(),
            distance * distance / light_normal.dot(direct).abs(),
            2.0e-12,
            "equal-IOR dA/dOmega",
        );
    }

    #[test]
    fn normal_glass_slab_has_the_closed_form_paraxial_jacobian() {
        let source = Point3::new(0.0, 0.0, 0.0);
        let target = Point3::new(0.0, 0.0, 3.0);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let thickness = 0.6;
        let ratio = 2.0 / 3.0;
        let connection =
            slab_connection_to_point(source, target, axis, axis, thickness, ratio).unwrap();
        assert_eq!(connection.incident_direction, axis);
        assert_eq!(connection.internal_direction, axis);
        let effective_paraxial_distance = 3.0 - thickness + thickness * ratio;
        assert_near(
            connection.light_area_jacobian.unwrap(),
            effective_paraxial_distance * effective_paraxial_distance,
            3.0e-14,
            "normal-incidence slab Jacobian",
        );
    }

    fn slab_ray_plane_endpoint(
        source: Point3,
        direction: Vec3,
        axis: Vec3,
        thickness: f64,
        ratio: f64,
        plane_point: Point3,
        plane_normal: Vec3,
    ) -> Point3 {
        let cosine = axis.dot(direction);
        let tangent = vec_sub(direction, axis.scale(cosine));
        let internal_cosine = (1.0 - ratio * ratio * tangent.dot(tangent)).sqrt();
        let shift = tangent.scale(thickness * (ratio / internal_cosine - 1.0 / cosine));
        let shifted_origin = source.offset(shift);
        let parameter =
            plane_point.delta_from(shifted_origin).dot(plane_normal) / plane_normal.dot(direction);
        shifted_origin.offset(direction.scale(parameter))
    }

    #[test]
    fn finite_light_slab_jacobian_matches_an_independent_ray_differential() {
        let source = Point3::new(-0.2, 0.1, -0.7);
        let target = Point3::new(0.8, -0.35, 2.4);
        let axis = Vec3::new(0.0, 0.0, 1.0);
        let light_normal = unit(Vec3::new(0.15, -0.25, 1.0));
        let thickness = 0.45;
        let ratio = 1.0 / 1.52;
        let connection =
            slab_connection_to_point(source, target, light_normal, axis, thickness, ratio).unwrap();
        let (u, v) = basis_all_sphere(connection.incident_direction);
        let epsilon = 2.0e-6;
        let endpoint = |tangent: Vec3, sign: f64| {
            let perturbed = unit(vec_add(
                connection.incident_direction,
                tangent.scale(sign * epsilon),
            ));
            slab_ray_plane_endpoint(
                source,
                perturbed,
                axis,
                thickness,
                ratio,
                target,
                light_normal,
            )
        };
        let du = endpoint(u, 1.0)
            .delta_from(endpoint(u, -1.0))
            .scale(0.5 / epsilon);
        let dv = endpoint(v, 1.0)
            .delta_from(endpoint(v, -1.0))
            .scale(0.5 / epsilon);
        let finite_difference = cross(du, dv).norm();
        assert_near(
            connection.light_area_jacobian.unwrap(),
            finite_difference,
            2.0e-8 * finite_difference,
            "finite-difference dA/dOmega",
        );
    }

    #[test]
    fn slab_mis_uses_both_fresnel_event_probabilities_without_bias() {
        let nee_pdf = 0.73;
        let source_bsdf_pdf = 0.29;
        let two_interface_transmission_probability = 0.81 * 0.76;
        let nee_weight = slab_nee_mis_weight(
            DirectStrategy::Mis,
            nee_pdf,
            source_bsdf_pdf,
            two_interface_transmission_probability,
            true,
        )
        .unwrap();
        let bsdf_weight = completed_slab_mis_weight(
            DirectStrategy::Mis,
            SmoothSlabExit {
                source_origin: Point3::new(0.0, 0.0, 0.0),
                source_geometric_normal: Vec3::new(0.0, 0.0, 1.0),
                source_direction: Vec3::new(0.0, 0.0, 1.0),
                source_pdf_solid_angle: source_bsdf_pdf,
                slab: ParallelSlab {
                    boundary_primitive: 0,
                    glass: DielectricGlass::representative_crown(),
                    entry_reference: Point3::new(0.0, 0.0, 0.0),
                    exit_reference: Point3::new(0.0, 0.0, 0.1),
                    axis: Vec3::new(0.0, 0.0, 1.0),
                    thickness_m: 0.1,
                },
                transmission_probability: two_interface_transmission_probability,
            },
            nee_pdf,
        );
        assert_near(
            nee_weight + bsdf_weight,
            1.0,
            2.0 * f64::EPSILON,
            "paired slab MIS weights",
        );
    }

    fn closed_test_slab(z_bottom: f64, z_top: f64) -> TriMesh {
        closed_test_slab_with_extent(z_bottom, z_top, 4.0)
    }

    fn closed_test_slab_with_extent(z_bottom: f64, z_top: f64, half_extent: f64) -> TriMesh {
        let lo = -half_extent;
        let hi = half_extent;
        TriMesh::new(
            vec![
                [lo, lo, z_bottom],
                [hi, lo, z_bottom],
                [hi, hi, z_bottom],
                [lo, hi, z_bottom],
                [lo, lo, z_top],
                [hi, lo, z_top],
                [hi, hi, z_top],
                [lo, hi, z_top],
            ],
            vec![
                [0, 2, 1],
                [0, 3, 2],
                [4, 5, 6],
                [4, 6, 7],
                [0, 1, 5],
                [0, 5, 4],
                [1, 2, 6],
                [1, 6, 5],
                [2, 3, 7],
                [2, 7, 6],
                [3, 0, 4],
                [3, 4, 7],
            ],
        )
    }

    fn instanced_test_mesh(mesh: TriMesh, object_id: u64) -> Shape {
        let identity = hash_domain(
            "org.frankensim.fs-render.test.slab-mesh",
            &object_id.to_le_bytes(),
        );
        Shape::Instance(
            GeometryInstance::try_new(
                object_id,
                identity,
                SharedGeometry::mesh(mesh),
                RigidTransform::identity(),
            )
            .unwrap(),
        )
    }

    fn closed_test_wedge() -> TriMesh {
        let lo = -4.0;
        let hi = 4.0;
        TriMesh::new(
            vec![
                [lo, lo, 0.0],
                [hi, lo, 0.0],
                [hi, hi, 0.0],
                [lo, hi, 0.0],
                [lo, lo, 0.2],
                [hi, lo, 0.6],
                [hi, hi, 0.6],
                [lo, hi, 0.2],
            ],
            vec![
                [0, 2, 1],
                [0, 3, 2],
                [4, 5, 6],
                [4, 6, 7],
                [0, 1, 5],
                [0, 5, 4],
                [1, 2, 6],
                [1, 6, 5],
                [2, 3, 7],
                [2, 7, 6],
                [3, 0, 4],
                [3, 4, 7],
            ],
        )
    }

    fn chamfered_plate_face_mesh() -> TriMesh {
        // An open extrusion of the same eight-edge x/z profile used by the
        // Euler studio's beveled plate. The test intentionally claims only
        // local face geometry, not watertightness of this compact fixture.
        let profile = [
            [-0.95, 0.0],
            [0.95, 0.0],
            [1.0, 0.05],
            [1.0, 0.15],
            [0.95, 0.2],
            [-0.95, 0.2],
            [-1.0, 0.15],
            [-1.0, 0.05],
        ];
        let mut vertices = Vec::with_capacity(2 * profile.len());
        for [x, z] in profile {
            vertices.push([x, -1.0, z]);
            vertices.push([x, 1.0, z]);
        }
        let mut triangles = Vec::with_capacity(2 * profile.len());
        for index in 0..profile.len() {
            let next = (index + 1) % profile.len();
            let a = u32::try_from(2 * index).unwrap();
            let b = a + 1;
            let d = u32::try_from(2 * next).unwrap();
            let c = d + 1;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
        }
        TriMesh::new(vertices, triangles)
    }

    fn horizontal_quad(z: f64, half_extent: f64) -> TriMesh {
        horizontal_quad_at(z, 0.0, 0.0, half_extent)
    }

    fn horizontal_quad_at(z: f64, center_x: f64, center_y: f64, half_extent: f64) -> TriMesh {
        TriMesh::new(
            vec![
                [center_x - half_extent, center_y - half_extent, z],
                [center_x + half_extent, center_y - half_extent, z],
                [center_x + half_extent, center_y + half_extent, z],
                [center_x - half_extent, center_y + half_extent, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }

    fn reverse_finite_parallel_face_pair() -> TriMesh {
        // The refracted connection enters the narrow lower patch near x=.492,
        // exits the upper patch near x=.429, and the normal discovery probe
        // from the entry still reaches the upper patch. The corresponding
        // straight source-to-emitter ray crosses x=.505 at z=0 and x=.406 at
        // z=.2, so it genuinely misses both finite interfaces.
        TriMesh::new(
            vec![
                [-1.000, -0.2, 0.0],
                [0.500, -0.2, 0.0],
                [0.500, 0.2, 0.0],
                [-1.000, 0.2, 0.0],
                [0.420, -0.2, 0.2],
                [2.000, -0.2, 0.2],
                [2.000, 0.2, 0.2],
                [0.420, 0.2, 0.2],
            ],
            vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7]],
        )
    }

    fn slab_visibility_scene(with_blocker: bool) -> (Scene, PreparedDirectLight) {
        let white = lift_rgb([1.0, 1.0, 1.0]);
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.5).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let mut primitives = vec![Primitive {
            shape: instanced_test_mesh(closed_test_slab(0.0, 0.2), 1),
            material: Material::Dielectric {
                glass,
                surface: DielectricSurface::SMOOTH,
            },
            emission: None,
        }];
        if with_blocker {
            primitives.push(Primitive {
                shape: Shape::Mesh(horizontal_quad(1.0, 2.0)),
                material: Material::Lambertian { reflectance: white },
                emission: None,
            });
        }
        let light_primitive = primitives.len();
        let emission = (white, 2.0);
        primitives.push(Primitive {
            shape: Shape::Mesh(horizontal_quad(2.0, 2.0)),
            material: Material::Lambertian { reflectance: white },
            emission: Some(emission),
        });
        let source = Point3::new(0.0, 0.0, -1.0);
        let target = Point3::new(0.35, -0.15, 2.0);
        let displacement = target.delta_from(source);
        let distance = displacement.norm();
        let direct = PreparedDirectLight {
            direction: displacement.scale(1.0 / distance),
            emission,
            pdf_solid_angle: 0.8,
            target: DirectLightTarget::Rectangle {
                primitive_index: light_primitive,
                distance_m: distance,
                point: target,
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
        };
        (
            Scene {
                primitives,
                lights: Vec::new(),
                environment: None,
                camera: Camera {
                    eye: source,
                    forward: Vec3::new(0.0, 0.0, 1.0),
                    up: Vec3::new(0.0, 1.0, 0.0),
                    half_tan: 0.0,
                },
            },
            direct,
        )
    }

    #[test]
    fn slab_connection_traces_both_interfaces_and_respects_post_slab_visibility() {
        with_test_cx(|cx| {
            for (with_blocker, expected_visible) in [(false, true), (true, false)] {
                let (scene, direct) = slab_visibility_scene(with_blocker);
                let source = scene.camera.eye;
                let straight_ray = Ray {
                    origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                    dir: direct.direction,
                };
                let first = intersect(&scene, cx, &straight_ray, None)
                    .unwrap()
                    .expect("straight sample first reaches the slab");
                assert_eq!(first.primitive_index, 0);
                let slab = discover_parallel_slab(&scene, cx, &straight_ray, first, None).unwrap();
                let lane = trace_parallel_slab_direct_lane(
                    &scene,
                    cx,
                    None,
                    source,
                    Vec3::new(0.0, 0.0, 1.0),
                    slab,
                    direct,
                    550.0,
                )
                .unwrap();
                assert_eq!(lane.visible, expected_visible);
                assert!(lane.nee_pdf_solid_angle.is_finite());
                assert!(lane.nee_pdf_solid_angle > 0.0);
                assert_near(
                    lane.radiance_transport,
                    lane.transmission_probability,
                    2.0e-15,
                    "lossless reciprocal eta factors",
                );
            }
        });
    }

    #[test]
    fn environment_slab_connection_preserves_the_external_direction_and_pdf() {
        with_test_cx(|cx| {
            let (mut scene, _) = slab_visibility_scene(false);
            scene.primitives.truncate(1);
            let source = scene.camera.eye;
            let tangent = 0.3_f64;
            let direction = Vec3::new(tangent, 0.0, (1.0 - tangent * tangent).sqrt());
            let direct = PreparedDirectLight {
                direction,
                emission: (lift_rgb([1.0; 3]), 1.0),
                pdf_solid_angle: 0.37,
                target: DirectLightTarget::Environment,
            };
            let straight_ray = Ray {
                origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                dir: direction,
            };
            let first = intersect(&scene, cx, &straight_ray, None)
                .unwrap()
                .expect("environment ray first reaches the slab");
            let slab = discover_parallel_slab(&scene, cx, &straight_ray, first, None).unwrap();
            let lane = trace_parallel_slab_direct_lane(
                &scene,
                cx,
                None,
                source,
                Vec3::new(0.0, 0.0, 1.0),
                slab,
                direct,
                550.0,
            )
            .unwrap();
            assert!(lane.visible);
            assert!(parallel_directions(lane.incident_direction, direction));
            assert_eq!(
                lane.nee_pdf_solid_angle.to_bits(),
                direct.pdf_solid_angle.to_bits(),
                "parallel slab has an identity external angular map"
            );
        });
    }

    #[test]
    fn reverse_environment_slab_pdf_replays_the_forward_straight_shadow_gate() {
        with_test_cx(|cx| {
            let (mut scene, _) = slab_visibility_scene(false);
            scene.primitives.truncate(1);
            let source = scene.camera.eye;
            let tangent = 0.24_f64;
            let direction = Vec3::new(tangent, 0.0, (1.0 - tangent * tangent).sqrt());
            let straight_ray = Ray {
                origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                dir: direction,
            };
            let first = intersect(&scene, cx, &straight_ray, None)
                .unwrap()
                .expect("forward environment shadow must first meet the slab");
            let slab = discover_parallel_slab(&scene, cx, &straight_ray, first, None).unwrap();
            let retained_path = SmoothSlabExit {
                source_origin: source,
                source_geometric_normal: Vec3::new(0.0, 0.0, 1.0),
                source_direction: direction,
                source_pdf_solid_angle: 0.29,
                slab,
                transmission_probability: 0.74,
            };
            let direct_pdf = 0.37;
            assert_eq!(
                completed_slab_environment_nee_pdf(
                    &scene,
                    cx,
                    None,
                    retained_path,
                    direction,
                    direct_pdf,
                    550.0,
                )
                .unwrap()
                .to_bits(),
                direct_pdf.to_bits(),
                "reverse environment MIS must retain the forward-eligible NEE density"
            );
        });
    }

    #[test]
    fn slab_nee_refuses_rough_nonparallel_and_nested_dielectrics() {
        with_test_cx(|cx| {
            let (mut rough_scene, direct) = slab_visibility_scene(false);
            let Material::Dielectric { glass, .. } = rough_scene.primitives[0].material else {
                unreachable!();
            };
            rough_scene.primitives[0].material = Material::Dielectric {
                glass,
                surface: DielectricSurface::try_rough(0.15).unwrap(),
            };
            let source = rough_scene.camera.eye;
            let straight_ray = Ray {
                origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                dir: direct.direction,
            };
            let first = intersect(&rough_scene, cx, &straight_ray, None)
                .unwrap()
                .unwrap();
            assert!(matches!(
                discover_parallel_slab(&rough_scene, cx, &straight_ray, first, None),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "first dielectric blocker is rough",
                    ..
                })
            ));

            let (mut wedge_scene, direct) = slab_visibility_scene(false);
            wedge_scene.primitives[0].shape = instanced_test_mesh(closed_test_wedge(), 1);
            let straight_ray = Ray {
                origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                dir: direct.direction,
            };
            let first = intersect(&wedge_scene, cx, &straight_ray, None)
                .unwrap()
                .unwrap();
            assert!(matches!(
                discover_parallel_slab(&wedge_scene, cx, &straight_ray, first, None),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "opposite interface is not parallel with reversed orientation",
                    ..
                })
            ));

            let (mut nested_scene, direct) = slab_visibility_scene(true);
            nested_scene.primitives[1].material = Material::Dielectric {
                glass,
                surface: DielectricSurface::SMOOTH,
            };
            let straight_ray = Ray {
                origin: source.offset(Vec3::new(0.0, 0.0, RAY_EPS)),
                dir: direct.direction,
            };
            let first = intersect(&nested_scene, cx, &straight_ray, None)
                .unwrap()
                .unwrap();
            let slab = discover_parallel_slab(&nested_scene, cx, &straight_ray, first, None)
                .expect("first plate remains a supported slab");
            assert!(matches!(
                trace_parallel_slab_direct_lane(
                    &nested_scene,
                    cx,
                    None,
                    source,
                    Vec3::new(0.0, 0.0, 1.0),
                    slab,
                    direct,
                    550.0,
                ),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "post-slab connection encountered another dielectric",
                    ..
                })
            ));
        });
    }

    #[test]
    fn nested_medium_opaque_source_omits_optional_slab_nee_without_refusing_path() {
        with_test_cx(|cx| {
            let white = lift_rgb([1.0, 1.0, 1.0]);
            let glass = DielectricGlass::new(
                CauchyIor::try_constant(1.5).unwrap(),
                BeerLambertAbsorption::try_constant(0.0).unwrap(),
                GlassProvenance::Custom,
            );
            let emission = (white, 2.0);
            let light_primitive = 2;
            let scene = Scene {
                primitives: vec![
                    Primitive {
                        // Enclosing medium: the resumed path starts inside it.
                        shape: instanced_test_mesh(
                            closed_test_slab_with_extent(-1.0, 3.0, 8.0),
                            91,
                        ),
                        material: Material::Dielectric {
                            glass,
                            surface: DielectricSurface::SMOOTH,
                        },
                        emission: None,
                    },
                    Primitive {
                        // A second dielectric blocks the straight light sample.
                        shape: instanced_test_mesh(closed_test_slab(0.8, 1.0), 92),
                        material: Material::Dielectric {
                            glass,
                            surface: DielectricSurface::SMOOTH,
                        },
                        emission: None,
                    },
                    Primitive {
                        shape: Shape::Mesh(horizontal_quad(2.0, 2.0)),
                        material: Material::Lambertian { reflectance: white },
                        emission: Some(emission),
                    },
                    Primitive {
                        // The resumed camera ray reaches this opaque source surface
                        // without crossing either dielectric boundary.
                        shape: Shape::Mesh(horizontal_quad(0.0, 4.0)),
                        material: Material::Lambertian { reflectance: white },
                        emission: None,
                    },
                ],
                lights: vec![RectLight {
                    corner: Point3::new(-2.0, -2.0, 2.0),
                    edge_u: Vec3::new(4.0, 0.0, 0.0),
                    edge_v: Vec3::new(0.0, 4.0, 0.0),
                    prim: light_primitive,
                    emission,
                }],
                environment: None,
                camera: Camera {
                    eye: Point3::new(0.0, 0.0, 0.25),
                    forward: Vec3::new(0.0, 0.0, -1.0),
                    up: Vec3::new(0.0, 1.0, 0.0),
                    half_tan: 0.0,
                },
            };
            let lighting = AdmittedLighting::try_new(&scene.lights, None).unwrap();
            let settings = Settings {
                width: 1,
                height: 1,
                spp: 1,
                max_depth: 1,
                sampler: Sampler::Iid,
                strategy: DirectStrategy::Mis,
                seed: 0x4e45_5354_4544,
            };
            let mut medium_stack = MediumStack::new();
            medium_stack
                .push(MediumEntry {
                    boundary_primitive: 0,
                    glass,
                })
                .unwrap();
            let origin = Point3::new(0.0, 0.0, 0.25);
            let result = trace_path(
                &scene,
                &lighting,
                cx,
                &settings,
                1.0 / y_integral(),
                0,
                0,
                0.5,
                0.5,
                0.5,
                None,
                CameraPath::Legacy,
                false,
                Some(SpectralPathState {
                    ray: Ray {
                        origin,
                        dir: Vec3::new(0.0, 0.0, -1.0),
                    },
                    throughput: [1.0; PACKET],
                    previous_bsdf: None,
                    prev_origin: origin,
                    segment_origin: origin,
                    medium_stack,
                    rng: PathRng {
                        pixel: 0,
                        sample: 0,
                        dim: 1,
                        key: [0, 0],
                    },
                    next_depth: 0,
                    active_lane: None,
                }),
            );
            let sample = result.expect(
                "unsupported nested-media slab NEE must be omitted, not refuse the BSDF path",
            );
            assert!(sample.xyz.into_iter().all(f64::is_finite));
        });
    }

    #[test]
    fn slab_nee_refuses_a_smooth_chart_sphere_even_on_a_centered_chord() {
        with_test_cx(|cx| {
            let glass = DielectricGlass::representative_crown();
            let geometry_identity = hash_domain(
                "org.frankensim.fs-render.test.slab-chart",
                b"centered-sphere",
            );
            let sphere = GeometryInstance::try_new(
                41,
                geometry_identity,
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 0.5,
                }),
                RigidTransform::identity(),
            )
            .unwrap();
            let scene = Scene {
                primitives: vec![Primitive {
                    shape: Shape::Instance(sphere),
                    material: Material::Dielectric {
                        glass,
                        surface: DielectricSurface::SMOOTH,
                    },
                    emission: None,
                }],
                lights: Vec::new(),
                environment: None,
                camera: Camera {
                    eye: Point3::new(0.0, 0.0, -2.0),
                    forward: Vec3::new(0.0, 0.0, 1.0),
                    up: Vec3::new(0.0, 1.0, 0.0),
                    half_tan: 0.0,
                },
            };
            let ray = Ray {
                origin: scene.camera.eye,
                dir: scene.camera.forward,
            };
            let first = intersect(&scene, cx, &ray, None)
                .unwrap()
                .expect("centered ray must hit the certified sphere chart");
            assert!(matches!(
                discover_parallel_slab(&scene, cx, &ray, first, None),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "first slab interface lacks an instanced triangle-mesh face witness",
                    ..
                })
            ));
            let retained_path = SmoothSlabExit {
                source_origin: scene.camera.eye,
                source_geometric_normal: scene.camera.forward,
                source_direction: scene.camera.forward,
                source_pdf_solid_angle: 0.25,
                slab: ParallelSlab {
                    boundary_primitive: 0,
                    glass,
                    entry_reference: Point3::new(0.0, 0.0, -0.5),
                    exit_reference: Point3::new(0.0, 0.0, 0.5),
                    axis: Vec3::new(0.0, 0.0, 1.0),
                    thickness_m: 1.0,
                },
                transmission_probability: 0.8,
            };
            assert!(
                replay_forward_slab_eligibility(
                    &scene,
                    cx,
                    None,
                    retained_path,
                    scene.camera.forward,
                )
                .expect("an unsupported forward proposal is not a BSDF-path error")
                .is_none(),
                "the unsupported chart sphere must contribute zero competing NEE density"
            );
        });
    }

    #[test]
    fn paired_mesh_projection_scan_preserves_sequential_bound_bits() {
        let mesh = chamfered_plate_face_mesh();
        let first_axis = unit(Vec3::new(0.3, -0.4, 0.5));
        let second_axis = unit(Vec3::new(-0.2, 0.7, 0.1));
        let first = mesh_projection_bounds(&mesh, first_axis).unwrap();
        let second = mesh_projection_bounds(&mesh, second_axis).unwrap();
        let paired = mesh_projection_bounds_pair(&mesh, first_axis, second_axis).unwrap();

        assert_eq!(
            [paired.0.0.to_bits(), paired.0.1.to_bits()],
            [first.0.to_bits(), first.1.to_bits()]
        );
        assert_eq!(
            [paired.1.0.to_bits(), paired.1.1.to_bits()],
            [second.0.to_bits(), second.1.to_bits()]
        );
    }

    #[test]
    fn slab_nee_admits_chamfered_plate_faces_only_on_the_thin_axis() {
        with_test_cx(|cx| {
            let scene = Scene {
                primitives: vec![Primitive {
                    shape: instanced_test_mesh(chamfered_plate_face_mesh(), 57),
                    material: Material::Dielectric {
                        glass: DielectricGlass::representative_crown(),
                        surface: DielectricSurface::SMOOTH,
                    },
                    emission: None,
                }],
                lights: Vec::new(),
                environment: None,
                camera: Camera {
                    eye: Point3::new(0.0, 0.0, -1.0),
                    forward: Vec3::new(0.0, 0.0, 1.0),
                    up: Vec3::new(0.0, 1.0, 0.0),
                    half_tan: 0.0,
                },
            };
            let central_ray = Ray {
                origin: scene.camera.eye,
                dir: scene.camera.forward,
            };
            let central_hit = intersect(&scene, cx, &central_ray, None)
                .unwrap()
                .expect("central ray must hit the retained lower principal face");
            let slab = discover_parallel_slab(&scene, cx, &central_ray, central_hit, None)
                .expect("central principal faces are the plate's thin-axis slab pair");
            assert_near(slab.thickness_m, 0.2, 4.0 * f64::EPSILON, "plate thickness");

            let bevel_outward = unit(Vec3::new(1.0, 0.0, -1.0));
            let bevel_center = Point3::new(0.975, 0.0, 0.025);
            let bevel_ray = Ray {
                origin: bevel_center.offset(bevel_outward.scale(0.5)),
                dir: bevel_outward.scale(-1.0),
            };
            let bevel_hit = intersect(&scene, cx, &bevel_ray, None)
                .unwrap()
                .expect("normal ray must hit the lower-right chamfer face");
            assert!(matches!(
                discover_parallel_slab(&scene, cx, &bevel_ray, bevel_hit, None),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "first slab interface lacks admitted thin-axis support",
                    ..
                })
            ));

            let side_ray = Ray {
                origin: Point3::new(2.0, 0.0, 0.1),
                dir: Vec3::new(-1.0, 0.0, 0.0),
            };
            let side_hit = intersect(&scene, cx, &side_ray, None)
                .unwrap()
                .expect("normal ray must hit the retained vertical side face");
            assert!(matches!(
                discover_parallel_slab(&scene, cx, &side_ray, side_hit, None),
                Err(TracerError::UnsupportedSlabNee {
                    reason: "first slab interface lacks admitted thin-axis support",
                    ..
                })
            ));
        });
    }

    fn reverse_edge_fixture(
        cx: &Cx<'_>,
        with_opaque_straight_blocker: bool,
    ) -> (Scene, SmoothSlabExit, Point3) {
        let white = lift_rgb([1.0; 3]);
        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.5).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let source = Point3::new(1.0, 0.0, -1.0);
        let target = Point3::new(-0.485, 0.0, 2.0);
        let source_normal = Vec3::new(0.0, 0.0, 1.0);
        let slab_axis = source_normal;
        let connection =
            slab_connection_to_point(source, target, slab_axis, slab_axis, 0.2, 1.0 / 1.5).unwrap();
        let mut primitives = vec![Primitive {
            shape: instanced_test_mesh(reverse_finite_parallel_face_pair(), 71),
            material: Material::Dielectric {
                glass,
                surface: DielectricSurface::SMOOTH,
            },
            emission: None,
        }];
        if with_opaque_straight_blocker {
            let straight_halfway_x = source.x + (target.x - source.x) / 6.0;
            primitives.push(Primitive {
                shape: Shape::Mesh(horizontal_quad_at(-0.5, straight_halfway_x, 0.0, 1.0e-3)),
                material: Material::Lambertian { reflectance: white },
                emission: None,
            });
        }
        let light_primitive = primitives.len();
        let emission = (white, 1.0);
        primitives.push(Primitive {
            shape: Shape::Mesh(horizontal_quad(2.0, 2.0)),
            material: Material::Lambertian { reflectance: white },
            emission: Some(emission),
        });
        let light = RectLight {
            corner: Point3::new(-2.0, -2.0, 2.0),
            edge_u: Vec3::new(4.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 4.0, 0.0),
            prim: light_primitive,
            emission,
        };
        let scene = Scene {
            primitives,
            lights: vec![light],
            environment: None,
            camera: Camera {
                eye: source,
                forward: source_normal,
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.0,
            },
        };
        let connected_ray = Ray {
            origin: source.offset(source_normal.scale(RAY_EPS)),
            dir: connection.incident_direction,
        };
        let first = intersect(&scene, cx, &connected_ray, None)
            .unwrap()
            .expect("the physical refracted connection must enter the finite face pair");
        assert_eq!(first.primitive_index, 0);
        let slab = discover_parallel_slab(&scene, cx, &connected_ray, first, None).unwrap();
        (
            scene,
            SmoothSlabExit {
                source_origin: source,
                source_geometric_normal: source_normal,
                source_direction: connection.incident_direction,
                source_pdf_solid_angle: 0.31,
                slab,
                transmission_probability: 0.72,
            },
            target,
        )
    }

    #[test]
    fn reverse_slab_mis_has_no_competitor_when_the_straight_ray_misses_a_finite_plate() {
        with_test_cx(|cx| {
            let (scene, slab_path, target) = reverse_edge_fixture(cx, false);
            let lighting = AdmittedLighting::try_new(&scene.lights, None).unwrap();
            let straight = unit(target.delta_from(slab_path.source_origin));
            let first = intersect(
                &scene,
                cx,
                &Ray {
                    origin: slab_path
                        .source_origin
                        .offset(slab_path.source_geometric_normal.scale(RAY_EPS)),
                    dir: straight,
                },
                None,
            )
            .unwrap()
            .expect("the straight ray must reach the emitter after missing the finite plate");
            assert_eq!(first.primitive_index, scene.lights[0].prim);
            let nee_pdf = completed_slab_rectangle_nee_pdf(
                &scene, &lighting, cx, None, slab_path, 0, target, 550.0,
            )
            .unwrap();
            assert_eq!(nee_pdf.to_bits(), 0.0_f64.to_bits());
            assert_eq!(
                completed_slab_mis_weight(DirectStrategy::Mis, slab_path, nee_pdf).to_bits(),
                1.0_f64.to_bits(),
                "a forward-ineligible slab NEE technique must not downweight the BSDF path"
            );
        });
    }

    #[test]
    fn reverse_slab_mis_has_no_competitor_when_an_opaque_primitive_blocks_the_straight_ray() {
        with_test_cx(|cx| {
            let (scene, slab_path, target) = reverse_edge_fixture(cx, true);
            let lighting = AdmittedLighting::try_new(&scene.lights, None).unwrap();
            let straight = unit(target.delta_from(slab_path.source_origin));
            let first = intersect(
                &scene,
                cx,
                &Ray {
                    origin: slab_path
                        .source_origin
                        .offset(slab_path.source_geometric_normal.scale(RAY_EPS)),
                    dir: straight,
                },
                None,
            )
            .unwrap()
            .expect("the straight ray must meet the purpose-built blocker");
            assert_eq!(first.primitive_index, 1);
            assert!(!primitive_is_dielectric(&scene, first.primitive_index));
            let nee_pdf = completed_slab_rectangle_nee_pdf(
                &scene, &lighting, cx, None, slab_path, 0, target, 550.0,
            )
            .unwrap();
            assert_eq!(nee_pdf.to_bits(), 0.0_f64.to_bits());
            assert_eq!(
                completed_slab_mis_weight(DirectStrategy::Mis, slab_path, nee_pdf).to_bits(),
                1.0_f64.to_bits(),
                "an opaque-first straight shadow must remove the slab NEE competitor"
            );
        });
    }

    #[test]
    fn aov_albedo_cache_preserves_direct_spectral_round_trip_bits() {
        let reflectance = lift_rgb([0.2, 0.4, 0.8]);
        let material = Material::Lambertian { reflectance };
        let direct = reflectance.rgb().map(|value| value.clamp(0.0, 1.0));
        let cache = AovAlbedoCache::try_new(2, true).unwrap();
        let first = cache.get(0, material).unwrap().unwrap();
        let second = cache.get(0, material).unwrap().unwrap();

        assert_eq!(first.map(f64::to_bits), direct.map(f64::to_bits));
        assert_eq!(second.map(f64::to_bits), direct.map(f64::to_bits));
        assert_eq!(
            cache
                .get(
                    1,
                    Material::Dielectric {
                        glass: DielectricGlass::representative_crown(),
                        surface: DielectricSurface::SMOOTH,
                    },
                )
                .unwrap(),
            None
        );
    }

    #[test]
    fn material_content_identity_is_deterministic_and_value_complete() {
        let reflectance = lift_rgb([0.2, 0.4, 0.8]);
        let diffuse = Material::Lambertian { reflectance };
        let smoother = Material::Ggx {
            reflectance,
            alpha: 0.1,
        };
        let rougher = Material::Ggx {
            reflectance,
            alpha: 0.2,
        };

        assert_eq!(diffuse.content_identity(), diffuse.content_identity());
        assert_ne!(diffuse.content_identity(), smoother.content_identity());
        assert_ne!(smoother.content_identity(), rougher.content_identity());
    }

    #[test]
    fn conductor_material_identity_binds_optics_and_surface_independently() {
        let smooth = ConductorSurface::try_rough(0.08).unwrap();
        let rough = ConductorSurface::try_rough(0.24).unwrap();
        let tungsten = ConductorOptics::representative_tungsten();
        let stainless = ConductorOptics::representative_stainless_steel();
        let tungsten_smooth = Material::Conductor {
            optics: tungsten,
            surface: smooth,
        };
        let tungsten_rough = Material::Conductor {
            optics: tungsten,
            surface: rough,
        };
        let stainless_smooth = Material::Conductor {
            optics: stainless,
            surface: smooth,
        };

        assert_eq!(
            tungsten_smooth.content_identity(),
            tungsten_smooth.content_identity()
        );
        assert_ne!(
            tungsten_smooth.content_identity(),
            tungsten_rough.content_identity()
        );
        assert_ne!(
            tungsten_smooth.content_identity(),
            stainless_smooth.content_identity()
        );
    }

    #[test]
    fn ggx_visible_normal_sampling_matches_pdf_at_grazing_and_both_poles() {
        let materials = [
            Material::Ggx {
                reflectance: lift_rgb([0.72; 3]),
                alpha: 0.2,
            },
            Material::Conductor {
                optics: ConductorOptics::representative_tungsten(),
                surface: ConductorSurface::try_rough(0.2).unwrap(),
            },
        ];
        let near_south = Vec3::new(2.0e-5, -3.0e-5, -1.0);
        let near_south = near_south.scale(1.0 / near_south.norm());
        for material in materials {
            for normal in [
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -1.0),
                near_south,
            ] {
                let (tangent, _) = basis_all_sphere(normal);
                for cos_o in [1.0, 0.08] {
                    let sin_o = (1.0_f64 - cos_o * cos_o).sqrt();
                    let wo = Vec3::new(
                        tangent.x * sin_o + normal.x * cos_o,
                        tangent.y * sin_o + normal.y * cos_o,
                        tangent.z * sin_o + normal.z * cos_o,
                    );
                    let mut accepted = 0_usize;
                    for (u1, u2) in [
                        (0.0, 0.0),
                        (0.01, 0.13),
                        (0.08, 0.17),
                        (0.21, 0.43),
                        (0.55, 0.79),
                        (0.93, 0.97),
                    ] {
                        let Some((wi, sampled_pdf)) = bsdf_sample(&material, normal, wo, u1, u2)
                        else {
                            continue;
                        };
                        accepted += 1;
                        assert!(wi.x.is_finite() && wi.y.is_finite() && wi.z.is_finite());
                        assert!((wi.norm() - 1.0).abs() <= 3.0e-12);
                        assert!(normal.dot(wi) > 0.0);
                        let evaluated_pdf = bsdf_pdf(&material, normal, wo, wi);
                        let tolerance = 4.0e-12 * sampled_pdf.abs().max(1.0);
                        assert!(
                            (evaluated_pdf - sampled_pdf).abs() <= tolerance,
                            "GGX VNDF sample/PDF mismatch at normal={normal:?}, cos_o={cos_o}: sampled={sampled_pdf:.17e}, evaluated={evaluated_pdf:.17e}"
                        );
                    }
                    assert!(
                        accepted >= 4,
                        "too few admitted VNDF samples at normal={normal:?}, cos_o={cos_o}: {accepted}"
                    );
                }
            }
        }
    }

    #[test]
    fn ggx_visible_normal_directional_pdf_integrates_to_its_acceptance_mass() {
        const SIDE: u32 = 512;
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let cos_o = 0.08_f64;
        let wo = Vec3::new((1.0 - cos_o * cos_o).sqrt(), 0.0, cos_o);
        let material = Material::Ggx {
            reflectance: lift_rgb([0.72; 3]),
            alpha: 0.12,
        };

        let mut pdf_sum = 0.0;
        let mut accepted = 0_u64;
        for y in 0..SIDE {
            let z = (f64::from(y) + 0.5) / f64::from(SIDE);
            let radial = (1.0 - z * z).sqrt();
            for x in 0..SIDE {
                let azimuth = 2.0 * PI * (f64::from(x) + 0.5) / f64::from(SIDE);
                let wi = Vec3::new(radial * det::cos(azimuth), radial * det::sin(azimuth), z);
                pdf_sum += bsdf_pdf(&material, normal, wo, wi);

                let u1 = (f64::from(x) + 0.5) / f64::from(SIDE);
                let u2 = (f64::from(y) + 0.5) / f64::from(SIDE);
                if bsdf_sample(&material, normal, wo, u1, u2).is_some() {
                    accepted += 1;
                }
            }
        }
        let samples = f64::from(SIDE * SIDE);
        let integrated_mass = pdf_sum * (2.0 * PI) / samples;
        let sampled_mass = accepted as f64 / samples;
        assert!((0.0..=1.0).contains(&integrated_mass));
        assert!((0.0..=1.0).contains(&sampled_mass));
        assert!(
            (integrated_mass - sampled_mass).abs() <= 4.0e-3,
            "VNDF reflection density mass disagrees with sampler acceptance: integral={integrated_mass:.9}, sampled={sampled_mass:.9}"
        );
    }

    fn sample_ggx_ndf_reflection_for_comparison(
        n: Vec3,
        wo: Vec3,
        alpha: f64,
        u1: f64,
        u2: f64,
    ) -> Option<(Vec3, f64)> {
        let alpha_squared = alpha * alpha;
        let cos_m_squared = ((1.0 - u1) / (u1 * (alpha_squared - 1.0) + 1.0)).clamp(0.0, 1.0);
        let cos_m = cos_m_squared.sqrt();
        let sin_m = (1.0 - cos_m_squared).sqrt();
        let phi = 2.0 * PI * u2;
        let m = to_world_all_sphere(n, [sin_m * det::cos(phi), sin_m * det::sin(phi), cos_m]);
        let wo_dot_m = wo.dot(m);
        if wo_dot_m <= 0.0 {
            return None;
        }
        let wi = Vec3::new(
            2.0 * wo_dot_m * m.x - wo.x,
            2.0 * wo_dot_m * m.y - wo.y,
            2.0 * wo_dot_m * m.z - wo.z,
        );
        if n.dot(wi) <= 0.0 {
            return None;
        }
        let pdf = ggx_d(alpha, n.dot(m)) * n.dot(m).max(0.0) / (4.0 * wo_dot_m);
        (pdf > 0.0).then_some((wi, pdf))
    }

    fn conductor_furnace_sample(
        material: &Material,
        normal: Vec3,
        wo: Vec3,
        sample: Option<(Vec3, f64)>,
    ) -> f64 {
        let Some((wi, pdf)) = sample else {
            return 0.0;
        };
        opaque_bsdf_eval(material, normal, wo, wi, 550.0, None).unwrap() * normal.dot(wi) / pdf
    }

    #[test]
    fn ggx_vndf_reduces_equal_cost_grazing_conductor_variance_without_bias() {
        const SAMPLES: u32 = 1 << 16;
        const REFERENCE_SIDE: u32 = 512;
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let cos_o = 0.08_f64;
        let wo = Vec3::new((1.0 - cos_o * cos_o).sqrt(), 0.0, cos_o);
        let alpha = 0.12;
        let material = Material::Conductor {
            optics: ConductorOptics::representative_tungsten(),
            surface: ConductorSurface::try_rough(alpha).unwrap(),
        };
        let mut ndf = AdaptivePixelAccumulator::EMPTY;
        let mut vndf = AdaptivePixelAccumulator::EMPTY;
        for sample in 0..SAMPLES {
            let draws = philox4x32_10(
                [sample, 0x766e_6466, 0x6767_7802, 0],
                [0x4752_415a, 0x494e_475f],
            );
            let u1 = u32_unit(draws[0]);
            let u2 = u32_unit(draws[1]);
            let ndf_sample = conductor_furnace_sample(
                &material,
                normal,
                wo,
                sample_ggx_ndf_reflection_for_comparison(normal, wo, alpha, u1, u2),
            );
            let vndf_sample = conductor_furnace_sample(
                &material,
                normal,
                wo,
                bsdf_sample(&material, normal, wo, u1, u2),
            );
            ndf.push([ndf_sample; 3]).unwrap();
            vndf.push([vndf_sample; 3]).unwrap();
        }

        let mut reference = 0.0;
        for y in 0..REFERENCE_SIDE {
            for x in 0..REFERENCE_SIDE {
                let u1 = (f64::from(x) + 0.5) / f64::from(REFERENCE_SIDE);
                let u2 = (f64::from(y) + 0.5) / f64::from(REFERENCE_SIDE);
                let (wi, _) = cosine_sample(normal, u1, u2);
                reference += opaque_bsdf_eval(&material, normal, wo, wi, 550.0, None).unwrap() * PI;
            }
        }
        reference /= f64::from(REFERENCE_SIDE * REFERENCE_SIDE);

        let ndf_mean = ndf.mean_xyz()[0];
        let vndf_mean = vndf.mean_xyz()[0];
        let ndf_error = (ndf_mean - reference).abs();
        let vndf_error = (vndf_mean - reference).abs();
        let ndf_variance = ndf.sample_variance_xyz()[0];
        let vndf_variance = vndf.sample_variance_xyz()[0];
        assert!(
            ndf_error <= 1.5e-2 && vndf_error <= 1.5e-2,
            "GGX estimators disagree with the independent furnace integral: reference={reference:.9}, NDF={:.9}, VNDF={:.9}",
            ndf_mean,
            vndf_mean,
        );
        assert!(
            vndf_variance < 0.35 * ndf_variance,
            "VNDF did not materially reduce equal-cost grazing variance: NDF={ndf_variance:.9e}, VNDF={vndf_variance:.9e}"
        );
        eprintln!(
            "GGX grazing furnace: reference={reference:.9}, NDF_mean={:.9}, VNDF_mean={:.9}, NDF_variance={ndf_variance:.9e}, VNDF_variance={vndf_variance:.9e}",
            ndf_mean, vndf_mean,
        );
    }

    #[test]
    fn conductor_bsdf_is_reciprocal_medium_aware_and_energy_bounded() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo_raw = Vec3::new(0.23, -0.11, 0.97);
        let wi_raw = Vec3::new(-0.31, 0.19, 0.93);
        let wo = wo_raw.scale(1.0 / wo_raw.norm());
        let wi = wi_raw.scale(1.0 / wi_raw.norm());
        for optics in [
            ConductorOptics::representative_tungsten(),
            ConductorOptics::representative_stainless_steel(),
        ] {
            let material = Material::Conductor {
                optics,
                surface: ConductorSurface::try_rough(0.18).unwrap(),
            };
            let forward = opaque_bsdf_eval(&material, normal, wo, wi, 550.0, None).unwrap();
            let reverse = opaque_bsdf_eval(&material, normal, wi, wo, 550.0, None).unwrap();
            assert!(forward.is_finite() && forward >= 0.0);
            assert!((forward - reverse).abs() <= 2.0e-14 * forward.abs().max(1.0));

            let crown = opaque_bsdf_eval(
                &material,
                normal,
                wo,
                wi,
                550.0,
                Some(DielectricGlass::representative_crown()),
            )
            .unwrap();
            assert_ne!(
                forward.to_bits(),
                crown.to_bits(),
                "conductor Fresnel ignored the active incident medium"
            );

            let mut reflected_energy = 0.0;
            const SIDE: u32 = 48;
            for y in 0..SIDE {
                for x in 0..SIDE {
                    let u1 = (f64::from(x) + 0.5) / f64::from(SIDE);
                    let u2 = (f64::from(y) + 0.5) / f64::from(SIDE);
                    let (sampled_wi, _) = cosine_sample(normal, u1, u2);
                    reflected_energy +=
                        opaque_bsdf_eval(&material, normal, normal, sampled_wi, 550.0, None)
                            .unwrap()
                            * PI;
                }
            }
            reflected_energy /= f64::from(SIDE * SIDE);
            assert!(
                reflected_energy.is_finite() && (0.0..=1.0 + 1.0e-10).contains(&reflected_energy),
                "single-scattering conductor furnace escaped its energy bound: {reflected_energy:.17e}"
            );
        }
    }

    #[test]
    fn lighting_admission_cancellation_preserves_tracer_cancellation_authority() {
        assert_eq!(
            TracerError::from(LightingError::Cancelled),
            TracerError::Cancelled
        );
    }

    #[test]
    fn mis_weights_only_emitters_reachable_by_the_nee_technique() {
        let bsdf_pdf = 0.25;
        let nee_pdf = 0.75;
        let previous = Some(PreviousBsdf {
            pdf: bsdf_pdf,
            delta: false,
            opaque_source_geometric_normal: Some(Vec3::new(0.0, 0.0, 1.0)),
            smooth_slab: None,
        });
        assert_eq!(
            emissive_hit_weight(DirectStrategy::Mis, previous, None).to_bits(),
            1.0_f64.to_bits(),
            "an emitter absent from NEE has no competing technique"
        );
        assert_eq!(
            emissive_hit_weight(DirectStrategy::Mis, previous, Some(nee_pdf)).to_bits(),
            balance_heuristic(1, bsdf_pdf, 1, nee_pdf).to_bits()
        );
        assert_eq!(
            emissive_hit_weight(DirectStrategy::NeeOnly, previous, None).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            emissive_hit_weight(
                DirectStrategy::NeeOnly,
                Some(PreviousBsdf {
                    pdf: 0.0,
                    delta: true,
                    opaque_source_geometric_normal: None,
                    smooth_slab: None,
                }),
                Some(nee_pdf),
            )
            .to_bits(),
            1.0_f64.to_bits(),
            "a delta path has no competing solid-angle NEE technique"
        );
    }

    #[test]
    fn cinematic_lens_stream_v1_has_frozen_replay_bits_and_keying() {
        let sample = camera_lens_sample([0x89ab_cdef, 0x0123_4567], 17, 29).unwrap();
        assert_eq!(sample.u().to_bits(), 0x3f9f_babf_b800_0000);
        assert_eq!(sample.v().to_bits(), 0x3fe1_6175_5c60_0000);

        let reseeded = camera_lens_sample([0x89ab_cdee, 0x0123_4567], 17, 29).unwrap();
        assert_ne!(
            [sample.u().to_bits(), sample.v().to_bits()],
            [reseeded.u().to_bits(), reseeded.v().to_bits()],
            "lens stream ignored the explicit render key"
        );
    }

    #[test]
    fn dielectric_offset_is_independent_of_unrelated_world_translation() {
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let positive = Vec3::new(0.0, 0.0, 1.0);
        let negative = Vec3::new(0.0, 0.0, -1.0);
        let local = Point3::new(0.0, 0.0, 0.0);
        let translated = Point3::new(1.0e9, 0.0, 0.0);

        let local_positive = dielectric_spawn_origin(local, normal, positive);
        let translated_positive = dielectric_spawn_origin(translated, normal, positive);
        let local_negative = dielectric_spawn_origin(local, normal, negative);
        let translated_negative = dielectric_spawn_origin(translated, normal, negative);

        assert_eq!(translated_positive.x.to_bits(), translated.x.to_bits());
        assert_eq!(translated_negative.x.to_bits(), translated.x.to_bits());
        assert_eq!(translated_positive.z.to_bits(), local_positive.z.to_bits());
        assert_eq!(translated_negative.z.to_bits(), local_negative.z.to_bits());
        assert!(local_positive.z > RAY_EPS);
        assert!(local_negative.z < -RAY_EPS);
    }

    #[test]
    fn exiting_medium_selects_ambient_or_the_actual_enclosing_layer() {
        let glass = DielectricGlass::representative_crown();
        let inner = MediumEntry {
            boundary_primitive: 7,
            glass,
        };

        let mut single_layer = MediumStack::new();
        single_layer.push(inner).unwrap();
        let exit_to_ambient = boundary_media(7, glass, false, &single_layer).unwrap();
        assert_eq!(exit_to_ambient.incident, Some(glass));
        assert_eq!(exit_to_ambient.transmitted, None);

        let mut nested = MediumStack::new();
        nested
            .push(MediumEntry {
                boundary_primitive: 3,
                glass,
            })
            .unwrap();
        nested.push(inner).unwrap();
        let exit_to_enclosing_layer = boundary_media(7, glass, false, &nested).unwrap();
        assert_eq!(exit_to_enclosing_layer.incident, Some(glass));
        assert_eq!(exit_to_enclosing_layer.transmitted, Some(glass));
    }

    fn test_entry_boundary(glass: DielectricGlass) -> BoundaryMedia {
        BoundaryMedia {
            incident: None,
            transmitted: Some(glass),
            transition: MediumTransition::Enter(MediumEntry {
                boundary_primitive: 7,
                glass,
            }),
        }
    }

    fn test_exit_boundary(glass: DielectricGlass) -> BoundaryMedia {
        BoundaryMedia {
            incident: Some(glass),
            transmitted: None,
            transition: MediumTransition::Exit {
                boundary_primitive: 7,
            },
        }
    }

    fn test_oblique_wo() -> Vec3 {
        let raw = Vec3::new(0.6, 0.0, 0.8);
        raw.scale(1.0 / raw.norm())
    }

    fn assert_relative_close(observed: f64, expected: f64, tolerance: f64, context: &str) {
        assert!(
            (observed - expected).abs() <= tolerance * observed.abs().max(expected.abs()).max(1.0),
            "{context}: observed={observed:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn dispersive_split_is_unbiased_deterministic_and_uses_lane_snell_directions() {
        let boundary = test_entry_boundary(DielectricGlass::representative_crown());
        let lambdas = hero_wavelengths(411.25, PACKET, LAMBDA_MIN, LAMBDA_MAX);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = test_oblique_wo();
        let reflected = sample_dispersive_dielectric_lanes(
            normal,
            wo,
            DielectricSurface::SMOOTH,
            &boundary,
            &lambdas,
            0.2,
            0.7,
            0.0,
        )
        .unwrap();
        let transmitted = sample_dispersive_dielectric_lanes(
            normal,
            wo,
            DielectricSurface::SMOOTH,
            &boundary,
            &lambdas,
            0.2,
            0.7,
            FORCED_TRANSMISSION_EVENT_SAMPLE,
        )
        .unwrap();
        let replay = sample_dispersive_dielectric_lanes(
            normal,
            wo,
            DielectricSurface::SMOOTH,
            &boundary,
            &lambdas,
            0.2,
            0.7,
            FORCED_TRANSMISSION_EVENT_SAMPLE,
        )
        .unwrap();
        for lane in 0..PACKET {
            let reflection = reflected[lane].unwrap();
            let continuation = transmitted[lane].unwrap();
            let replayed = replay[lane].unwrap();
            assert_eq!(reflection.event, DielectricEvent::Reflection);
            assert_eq!(continuation.event, DielectricEvent::Transmission);
            assert_eq!(continuation.direction, replayed.direction);
            assert_eq!(continuation.weight.to_bits(), replayed.weight.to_bits());
            let eta_i = medium_ior(boundary.incident, lambdas[lane]).unwrap();
            let eta_t = medium_ior(boundary.transmitted, lambdas[lane]).unwrap();
            let fresnel = fresnel_dielectric(normal.dot(wo), eta_i, eta_t).unwrap();
            let expected = fresnel.reflectance
                + (1.0 - fresnel.reflectance) * (eta_i / eta_t) * (eta_i / eta_t);
            let observed = fresnel.reflectance * reflection.weight
                + (1.0 - fresnel.reflectance) * continuation.weight;
            assert_relative_close(observed, expected, 4.0e-14, "smooth split expectation");
        }
        assert_ne!(
            transmitted[0].unwrap().direction.x.to_bits(),
            transmitted[3].unwrap().direction.x.to_bits(),
            "dispersive companions followed one refracted direction"
        );
    }

    #[test]
    fn rough_split_proposal_replays_bsdf_and_children_cannot_resplit() {
        let boundary = test_entry_boundary(DielectricGlass::representative_crown());
        let lambdas = hero_wavelengths(437.0, PACKET, LAMBDA_MIN, LAMBDA_MAX);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = test_oblique_wo();
        let split = sample_dispersive_dielectric_lanes(
            normal,
            wo,
            DielectricSurface::POLISHED,
            &boundary,
            &lambdas,
            0.31,
            0.67,
            FORCED_TRANSMISSION_EVENT_SAMPLE,
        )
        .unwrap();
        assert!(should_split_dispersive_boundary(None, &boundary, &lambdas).unwrap());
        for lane in 0..PACKET {
            assert!(
                !should_split_dispersive_boundary(Some(lane), &boundary, &lambdas).unwrap(),
                "rough monochromatic child lane {lane} was allowed to split again"
            );
            let continuation = split[lane].unwrap();
            assert_eq!(continuation.event, DielectricEvent::Transmission);
            let eta_i = medium_ior(boundary.incident, lambdas[lane]).unwrap();
            let eta_t = medium_ior(boundary.transmitted, lambdas[lane]).unwrap();
            let evaluation = evaluate_rough_dielectric(
                normal,
                wo,
                continuation.direction,
                eta_i,
                eta_t,
                DielectricSurface::POLISHED.roughness_alpha().unwrap(),
            )
            .unwrap();
            assert_relative_close(
                evaluation.pdf,
                continuation.pdf,
                2.0e-12,
                "lane-native proposal PDF",
            );
            assert_relative_close(
                continuation.weight,
                evaluation.value * normal.dot(continuation.direction).abs() / continuation.pdf,
                2.0e-12,
                "rough split sample/PDF replay",
            );
        }
    }

    #[test]
    fn dispersive_exit_split_preserves_companion_transmission_across_tir_threshold() {
        let glass = DielectricGlass::representative_crown();
        let boundary = test_exit_boundary(glass);
        let lambdas = hero_wavelengths(411.25, PACKET, LAMBDA_MIN, LAMBDA_MAX);
        let mut high_ior_lane = 0;
        let mut low_ior_lane = 0;
        let mut high_ior = medium_ior(boundary.incident, lambdas[0]).unwrap();
        let mut low_ior = high_ior;
        for (lane, &lambda) in lambdas.iter().enumerate().skip(1) {
            let ior = medium_ior(boundary.incident, lambda).unwrap();
            if ior > high_ior {
                high_ior = ior;
                high_ior_lane = lane;
            }
            if ior < low_ior {
                low_ior = ior;
                low_ior_lane = lane;
            }
        }
        assert!(high_ior > low_ior);
        let incident_sine = 0.5 * (high_ior.recip() + low_ior.recip());
        let incident_cosine = (1.0 - incident_sine * incident_sine).sqrt();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = Vec3::new(incident_sine, 0.0, incident_cosine);
        assert!(should_split_dispersive_boundary(None, &boundary, &lambdas).unwrap());

        for surface in [DielectricSurface::SMOOTH, DielectricSurface::POLISHED] {
            // A nearly normal sampled microfacet keeps the rough case inside
            // the same interval between the two wavelength critical angles.
            let split = sample_dispersive_dielectric_lanes(
                normal,
                wo,
                surface,
                &boundary,
                &lambdas,
                1.0e-12,
                0.37,
                FORCED_TRANSMISSION_EVENT_SAMPLE,
            )
            .unwrap();
            let tir = split[high_ior_lane].unwrap();
            let transmitting = split[low_ior_lane].unwrap();
            assert_eq!(tir.event, DielectricEvent::Reflection);
            assert_eq!(transmitting.event, DielectricEvent::Transmission);
            assert!(tir.direction.z > 0.0);
            assert!(transmitting.direction.z < 0.0);
            assert!(transmitting.weight > 0.0);
        }
    }

    #[test]
    fn split_work_is_bounded_and_nondispersive_boundaries_retain_packet_bits() {
        for max_depth in 0..64_u32 {
            let bound = maximum_spectral_traversals(max_depth);
            assert!(u64::from(max_depth) <= bound);
            for split_depth in 0..max_depth {
                let prefix = u64::from(split_depth + 1);
                let suffix = u64::from(max_depth - split_depth - 1);
                for children in 0..=PACKET as u64 {
                    assert!(prefix + children * suffix <= bound);
                }
            }
        }

        let glass = DielectricGlass::new(
            CauchyIor::try_constant(1.52).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        );
        let boundary = test_entry_boundary(glass);
        let lambdas = hero_wavelengths(463.0, PACKET, LAMBDA_MIN, LAMBDA_MAX);
        assert!(!boundary_direction_is_wavelength_dependent(&boundary, &lambdas).unwrap());
        assert!(!should_split_dispersive_boundary(None, &boundary, &lambdas).unwrap());
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = test_oblique_wo();
        let first = sample_dielectric_path(
            normal,
            wo,
            DielectricSurface::POLISHED,
            &boundary,
            &lambdas,
            0,
            0.23,
            0.71,
            FORCED_TRANSMISSION_EVENT_SAMPLE,
        )
        .unwrap()
        .unwrap();
        let last = sample_dielectric_path(
            normal,
            wo,
            DielectricSurface::POLISHED,
            &boundary,
            &lambdas,
            PACKET - 1,
            0.23,
            0.71,
            FORCED_TRANSMISSION_EVENT_SAMPLE,
        )
        .unwrap()
        .unwrap();
        assert_eq!(first.direction, last.direction);
        assert_eq!(first.pdf.to_bits(), last.pdf.to_bits());
        assert_eq!(
            first.weights.map(f64::to_bits),
            last.weights.map(f64::to_bits)
        );
    }

    #[test]
    fn four_lane_split_reduces_equal_sample_chroma_variance_without_mean_shift() {
        const SAMPLES: u32 = 4_096;
        let boundary = test_entry_boundary(DielectricGlass::representative_crown());
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let wo = test_oblique_wo();
        let mut collapsed = AdaptivePixelAccumulator::EMPTY;
        let mut split_accumulator = AdaptivePixelAccumulator::EMPTY;
        for sample_index in 0..SAMPLES {
            let hero = LAMBDA_MIN
                + (f64::from(sample_index) + 0.5) / f64::from(SAMPLES) * (LAMBDA_MAX - LAMBDA_MIN);
            let lambdas = hero_wavelengths(hero, PACKET, LAMBDA_MIN, LAMBDA_MAX);
            let draw = philox4x32_10(
                [sample_index, 0x7370_6c69, 0x745f_7633, 0],
                [0x4c41_4e45, 0x5350_4c49],
            );
            let sampled = sample_dielectric_path(
                normal,
                wo,
                DielectricSurface::SMOOTH,
                &boundary,
                &lambdas,
                0,
                0.2,
                0.7,
                u32_unit(draw[0]),
            )
            .unwrap()
            .unwrap();
            let mut old_xyz = [0.0; 3];
            let mut split_xyz = [0.0; 3];
            if sampled.event == DielectricEvent::Reflection {
                for lane in 0..PACKET {
                    let basis = [
                        cie_x(lambdas[lane]),
                        cie_y(lambdas[lane]),
                        cie_z(lambdas[lane]),
                    ];
                    for channel in 0..3 {
                        old_xyz[channel] += sampled.weights[lane] * basis[channel];
                    }
                }
            } else {
                old_xyz = [
                    PACKET as f64 * sampled.weights[0] * cie_x(lambdas[0]),
                    PACKET as f64 * sampled.weights[0] * cie_y(lambdas[0]),
                    PACKET as f64 * sampled.weights[0] * cie_z(lambdas[0]),
                ];
            }
            let continuations = sample_dispersive_dielectric_lanes(
                normal,
                wo,
                DielectricSurface::SMOOTH,
                &boundary,
                &lambdas,
                0.2,
                0.7,
                u32_unit(draw[0]),
            )
            .unwrap();
            for lane in 0..PACKET {
                let weight = continuations[lane].unwrap().weight;
                split_xyz[0] += weight * cie_x(lambdas[lane]);
                split_xyz[1] += weight * cie_y(lambdas[lane]);
                split_xyz[2] += weight * cie_z(lambdas[lane]);
            }
            collapsed.push(old_xyz).unwrap();
            split_accumulator.push(split_xyz).unwrap();
        }
        for channel in 0..3 {
            let mean_gap =
                (collapsed.mean_xyz()[channel] - split_accumulator.mean_xyz()[channel]).abs();
            let conservative_standard_error = ((collapsed.sample_variance_xyz()[channel]
                + split_accumulator.sample_variance_xyz()[channel])
                / f64::from(SAMPLES))
            .sqrt();
            assert!(mean_gap <= 5.0 * conservative_standard_error);
        }
        let old_chroma_variance =
            collapsed.sample_variance_xyz()[0] + collapsed.sample_variance_xyz()[2];
        let split_chroma_variance =
            split_accumulator.sample_variance_xyz()[0] + split_accumulator.sample_variance_xyz()[2];
        assert!(
            split_chroma_variance < 0.4 * old_chroma_variance,
            "split chroma variance was not materially lower: old={old_chroma_variance:.9e}, split={split_chroma_variance:.9e}"
        );
    }

    #[test]
    fn adaptive_policy_admission_canonicalizes_zero_and_fixes_checkpoints() {
        assert_eq!(
            AdaptiveSamplingConfig::try_new(1, 1, 0.0, 0.0, 0.0),
            Err(AdaptiveSamplingError::InvalidMinimumSamples)
        );
        assert_eq!(
            AdaptiveSamplingConfig::try_new(2, 0, 0.0, 0.0, 0.0),
            Err(AdaptiveSamplingError::InvalidBatchSamples)
        );
        for (field, value) in [
            ("absolute_error", -1.0),
            ("relative_error", f64::NAN),
            ("dark_floor", f64::INFINITY),
        ] {
            let values = match field {
                "absolute_error" => (value, 0.0, 0.0),
                "relative_error" => (0.0, value, 0.0),
                _ => (0.0, 0.0, value),
            };
            assert_eq!(
                AdaptiveSamplingConfig::try_new(2, 1, values.0, values.1, values.2),
                Err(AdaptiveSamplingError::InvalidThreshold { field })
            );
        }

        let policy = AdaptiveSamplingConfig::try_new(3, 4, -0.0, -0.0, -0.0)
            .expect("signed zero is a valid nonnegative threshold");
        assert_eq!(policy.absolute_error().to_bits(), 0.0_f64.to_bits());
        assert_eq!(policy.relative_error().to_bits(), 0.0_f64.to_bits());
        assert_eq!(policy.dark_floor().to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            (0..=12)
                .filter(|samples| policy.is_checkpoint(*samples, 12))
                .collect::<Vec<_>>(),
            vec![3, 7, 11, 12],
            "the hard ceiling is a final truncated checkpoint"
        );
        assert_eq!(
            policy.validate_maximum(2),
            Err(AdaptiveSamplingError::MaximumBelowMinimum)
        );
    }

    #[test]
    fn adaptive_welford_boundary_and_final_checkpoint_precedence_are_exact() {
        let mut pixel = AdaptivePixelAccumulator::EMPTY;
        pixel.push([1.0; 3]).unwrap();
        pixel.push([3.0; 3]).unwrap();
        assert_eq!(pixel.sum_xyz.map(f64::to_bits), [4.0; 3].map(f64::to_bits));
        assert_eq!(
            pixel.mean_xyz().map(f64::to_bits),
            [2.0; 3].map(f64::to_bits)
        );
        assert_eq!(
            pixel.sample_variance_xyz().map(f64::to_bits),
            [2.0; 3].map(f64::to_bits)
        );
        assert_eq!(
            pixel.dispersion_proxy_xyz().map(f64::to_bits),
            [1.0; 3].map(f64::to_bits)
        );

        let exact = AdaptiveSamplingConfig::try_new(2, 1, 1.0, 0.0, 0.0).unwrap();
        assert_eq!(
            pixel.decision(exact, 2),
            Some(AdaptiveDecision::ErrorThreshold),
            "threshold equality wins, including at the hard ceiling"
        );
        let next_down = f64::from_bits(1.0_f64.to_bits() - 1);
        let strict = AdaptiveSamplingConfig::try_new(2, 1, next_down, 0.0, 0.0).unwrap();
        assert_eq!(
            pixel.decision(strict, 2),
            Some(AdaptiveDecision::MaximumSamples)
        );
    }

    #[test]
    fn adaptive_accumulator_failures_are_transactional() {
        let mut pixel = AdaptivePixelAccumulator::EMPTY;
        pixel.push([1.0, 2.0, 3.0]).unwrap();
        let before_nan = pixel;
        assert_eq!(
            pixel.push([4.0, f64::NAN, 6.0]),
            Err(AdaptiveSamplingError::NonFiniteSample)
        );
        assert_eq!(pixel, before_nan);

        let mut later_channel_overflow = AdaptivePixelAccumulator {
            sum_xyz: [1.0, f64::MAX, 1.0],
            mean_xyz: [1.0, f64::MAX, 1.0],
            m2_xyz: [0.0; 3],
            samples: 1,
            decision: None,
        };
        let before_overflow = later_channel_overflow;
        assert_eq!(
            later_channel_overflow.push([2.0, f64::MAX, 2.0]),
            Err(AdaptiveSamplingError::InvalidMoment)
        );
        assert_eq!(later_channel_overflow, before_overflow);

        let mut count_overflow = AdaptivePixelAccumulator {
            samples: u32::MAX,
            ..AdaptivePixelAccumulator::EMPTY
        };
        let before_count = count_overflow;
        assert_eq!(
            count_overflow.push([0.0; 3]),
            Err(AdaptiveSamplingError::SampleCountOverflow)
        );
        assert_eq!(count_overflow, before_count);
    }

    #[test]
    fn adaptive_hdr_moments_and_power_of_two_scaling_remain_stable() {
        let mut hdr = AdaptivePixelAccumulator::EMPTY;
        for value in [1.0e12 + 1.0, 1.0e12 + 2.0, 1.0e12 + 3.0, 1.0e12 + 4.0] {
            hdr.push([value; 3]).unwrap();
        }
        let expected_variance = 5.0 / 3.0;
        for variance in hdr.sample_variance_xyz() {
            assert!(
                (variance - expected_variance).abs() <= 2.0 * f64::EPSILON,
                "centered HDR variance drifted: actual={variance:?} expected={expected_variance:?}"
            );
        }

        let mut base = AdaptivePixelAccumulator::EMPTY;
        let mut scaled = AdaptivePixelAccumulator::EMPTY;
        for sample in [[1.0, 2.0, 4.0], [3.0, 6.0, 8.0], [5.0, 10.0, 16.0]] {
            base.push(sample).unwrap();
            scaled.push(sample.map(|value| value * 8.0)).unwrap();
        }
        assert_eq!(
            scaled.sum_xyz.map(f64::to_bits),
            base.sum_xyz.map(|value| value * 8.0).map(f64::to_bits)
        );
        assert_eq!(
            scaled.mean_xyz.map(f64::to_bits),
            base.mean_xyz.map(|value| value * 8.0).map(f64::to_bits)
        );
        assert_eq!(
            scaled.m2_xyz.map(f64::to_bits),
            base.m2_xyz.map(|value| value * 64.0).map(f64::to_bits)
        );
        let base_policy = AdaptiveSamplingConfig::try_new(2, 1, 0.25, 0.1, 0.5).unwrap();
        let scaled_policy = AdaptiveSamplingConfig::try_new(2, 1, 2.0, 0.1, 4.0).unwrap();
        assert_eq!(base.meets(base_policy), scaled.meets(scaled_policy));
    }

    #[test]
    fn adaptive_dark_floor_is_per_channel_and_does_not_hide_a_noisy_channel() {
        let mut pixel = AdaptivePixelAccumulator::EMPTY;
        pixel.push([0.0, 0.0, 0.0]).unwrap();
        pixel.push([0.02, 0.0, 0.0]).unwrap();
        assert_eq!(
            pixel.dispersion_proxy_xyz().map(f64::to_bits),
            [0.01, 0.0, 0.0].map(f64::to_bits)
        );

        let without_floor = AdaptiveSamplingConfig::try_new(2, 1, 0.0, 0.01, 0.0).unwrap();
        assert!(
            !pixel.meets(without_floor),
            "one unresolved channel must keep the whole pixel active"
        );
        let with_floor = AdaptiveSamplingConfig::try_new(2, 1, 0.0, 0.01, 1.0).unwrap();
        assert!(
            pixel.meets(with_floor),
            "dark-floor equality should admit the declared channel threshold"
        );
    }
}

// ---- BSDF machinery --------------------------------------------------

/// Deterministic orthonormal basis from a unit normal (Frisvad's
/// branch on the pole, fixed arithmetic).
fn basis(n: Vec3) -> (Vec3, Vec3) {
    if n.z < -0.999_999_9 {
        return (Vec3::new(0.0, -1.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
    }
    let a = 1.0 / (1.0 + n.z);
    let b = -n.x * n.y * a;
    (
        Vec3::new(1.0 - n.x * n.x * a, b, -n.x),
        Vec3::new(b, 1.0 - n.y * n.y * a, -n.y),
    )
}

fn to_world(n: Vec3, local: [f64; 3]) -> Vec3 {
    let (t, b) = basis(n);
    Vec3::new(
        t.x * local[0] + b.x * local[1] + n.x * local[2],
        t.y * local[0] + b.y * local[1] + n.y * local[2],
        t.z * local[0] + b.z * local[1] + n.z * local[2],
    )
}

/// Cosine-weighted hemisphere sample around `n` using `det` trig (this
/// path feeds the frozen goldens; the crate-root helper uses platform
/// trig and stays for the un-hashed v0 batteries).
fn cosine_sample(n: Vec3, u1: f64, u2: f64) -> (Vec3, f64) {
    let r = u1.sqrt();
    let phi = 2.0 * PI * u2;
    let (sp, cp) = (det::sin(phi), det::cos(phi));
    let z = (1.0 - u1).max(0.0).sqrt();
    (to_world(n, [r * cp, r * sp, z]), z / PI)
}

fn ggx_d(alpha: f64, cos_m: f64) -> f64 {
    if cos_m <= 0.0 {
        return 0.0;
    }
    let a2 = alpha * alpha;
    let c2 = cos_m * cos_m;
    let t = c2 * (a2 - 1.0) + 1.0;
    a2 / (PI * t * t)
}

fn smith_g1(alpha: f64, cos_v: f64) -> f64 {
    let a2 = alpha * alpha;
    2.0 * cos_v / (cos_v + (a2 + (1.0 - a2) * cos_v * cos_v).sqrt())
}

fn schlick(f0: f64, cos_i: f64) -> f64 {
    let m = (1.0 - cos_i).clamp(0.0, 1.0);
    let m2 = m * m;
    let m5 = m2 * m2 * m; // explicit powers — never powi (hazard class)
    f0 + (1.0 - f0) * m5
}

fn bsdf_eval(mat: &Material, n: Vec3, wo: Vec3, wi: Vec3, lambda: f64) -> f64 {
    let (cos_o, cos_i) = (n.dot(wo), n.dot(wi));
    if cos_o <= 0.0 || cos_i <= 0.0 {
        return 0.0;
    }
    match mat {
        Material::Lambertian { reflectance } => reflectance.eval(lambda) / PI,
        Material::Ggx { reflectance, alpha } => {
            let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
            let hn = hsum.norm();
            if hn < 1e-12 {
                return 0.0;
            }
            let m = hsum.scale(1.0 / hn);
            let d = ggx_d(*alpha, n.dot(m));
            let g = smith_g1(*alpha, cos_o) * smith_g1(*alpha, cos_i);
            let f = schlick(reflectance.eval(lambda), wo.dot(m).max(0.0));
            d * g * f / (4.0 * cos_o * cos_i)
        }
        Material::Conductor { .. } => 0.0,
        Material::Dielectric { .. } => 0.0,
    }
}

fn opaque_bsdf_eval(
    mat: &Material,
    n: Vec3,
    wo: Vec3,
    wi: Vec3,
    lambda: f64,
    incident_medium: Option<DielectricGlass>,
) -> Result<f64, TracerError> {
    let Material::Conductor { optics, surface } = mat else {
        // Keep the frozen Lambertian/GGX arithmetic in `bsdf_eval`; merely
        // adding the conductor family must not perturb legacy image bits.
        return Ok(bsdf_eval(mat, n, wo, wi, lambda));
    };
    let (cos_o, cos_i) = (n.dot(wo), n.dot(wi));
    if cos_o <= 0.0 || cos_i <= 0.0 {
        return Ok(0.0);
    }
    let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
    let hn = hsum.norm();
    if hn < 1.0e-12 {
        return Ok(0.0);
    }
    let microfacet_normal = hsum.scale(1.0 / hn);
    let alpha = surface.roughness_alpha();
    let d = ggx_d(alpha, n.dot(microfacet_normal));
    let g = smith_g1(alpha, cos_o) * smith_g1(alpha, cos_i);
    let incident_eta = medium_ior(incident_medium, lambda)?;
    let fresnel = optics.fresnel(
        lambda,
        incident_eta,
        wo.dot(microfacet_normal).clamp(0.0, 1.0),
    )?;
    Ok(d * g * fresnel / (4.0 * cos_o * cos_i))
}

fn bsdf_pdf(mat: &Material, n: Vec3, wo: Vec3, wi: Vec3) -> f64 {
    let cos_i = n.dot(wi);
    if cos_i <= 0.0 || n.dot(wo) <= 0.0 {
        return 0.0;
    }
    match mat {
        Material::Lambertian { .. } => cos_i / PI,
        Material::Ggx { alpha, .. } => {
            let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
            let hn = hsum.norm();
            if hn < 1e-12 {
                return 0.0;
            }
            let m = hsum.scale(1.0 / hn);
            ggx_vndf_reflection_pdf(*alpha, n, wo, m)
        }
        Material::Conductor { surface, .. } => {
            let hsum = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
            let hn = hsum.norm();
            if hn < 1e-12 {
                return 0.0;
            }
            let m = hsum.scale(1.0 / hn);
            ggx_vndf_reflection_pdf(surface.roughness_alpha(), n, wo, m)
        }
        Material::Dielectric { .. } => 0.0,
    }
}

fn bsdf_sample(mat: &Material, n: Vec3, wo: Vec3, u1: f64, u2: f64) -> Option<(Vec3, f64)> {
    match mat {
        Material::Lambertian { .. } => {
            let (wi, pdf) = cosine_sample(n, u1, u2);
            (pdf > 0.0).then_some((wi, pdf))
        }
        Material::Ggx { alpha, .. } => sample_ggx_vndf_reflection(n, wo, *alpha, u1, u2),
        Material::Conductor { surface, .. } => {
            sample_ggx_vndf_reflection(n, wo, surface.roughness_alpha(), u1, u2)
        }
        Material::Dielectric { .. } => None,
    }
}

/// Duff et al.'s cancellation-free all-sphere basis. Unlike the legacy
/// Lambertian basis, this remains well conditioned at both poles and is used
/// by the view-dependent GGX sampler for every reflective microfacet material.
fn basis_all_sphere(normal: Vec3) -> (Vec3, Vec3) {
    let sign = if normal.z < 0.0 { -1.0 } else { 1.0 };
    let a = -1.0 / (sign + normal.z);
    let b = normal.x * normal.y * a;
    (
        Vec3::new(
            1.0 + sign * normal.x * normal.x * a,
            sign * b,
            -sign * normal.x,
        ),
        Vec3::new(b, sign + normal.y * normal.y * a, -normal.y),
    )
}

fn to_world_all_sphere(normal: Vec3, local: [f64; 3]) -> Vec3 {
    let (tangent, bitangent) = basis_all_sphere(normal);
    Vec3::new(
        tangent.x * local[0] + bitangent.x * local[1] + normal.x * local[2],
        tangent.y * local[0] + bitangent.y * local[1] + normal.y * local[2],
        tangent.z * local[0] + bitangent.z * local[1] + normal.z * local[2],
    )
}

/// Directional density induced by isotropic GGX visible-normal sampling.
///
/// The sampled microfacet density is
/// `D(m) G1(wo) (wo·m) / (n·wo)`. Reflection contributes the Jacobian
/// `1 / (4 wo·m)`, so the two dot products cancel. Samples whose reflected
/// direction falls below the macrosurface are rejected by the caller; the
/// remaining directional density is therefore intentionally sub-normalized by
/// exactly that rejection probability.
fn ggx_vndf_reflection_pdf(alpha: f64, n: Vec3, wo: Vec3, m: Vec3) -> f64 {
    let cos_o = n.dot(wo);
    if cos_o <= 0.0 || wo.dot(m) <= 0.0 {
        return 0.0;
    }
    ggx_d(alpha, n.dot(m)) * smith_g1(alpha, cos_o) / (4.0 * cos_o)
}

/// Sample Heitz's isotropic GGX distribution of visible normals and reflect
/// `wo` about the selected microfacet ("Sampling the GGX Distribution of
/// Visible Normals", JCGT 2018). This uses exactly the existing two BSDF
/// uniform dimensions, so determinism and cancellation checkpoints are
/// unchanged while the estimator's sample/PDF bits deliberately change.
fn sample_ggx_vndf_reflection(
    n: Vec3,
    wo: Vec3,
    alpha: f64,
    u1: f64,
    u2: f64,
) -> Option<(Vec3, f64)> {
    let cos_o = n.dot(wo);
    if !alpha.is_finite()
        || alpha <= 0.0
        || cos_o <= 0.0
        || !u1.is_finite()
        || !u2.is_finite()
        || !(0.0..=1.0).contains(&u1)
        || !(0.0..=1.0).contains(&u2)
    {
        return None;
    }

    let (tangent, bitangent) = basis_all_sphere(n);
    let local_wo = [tangent.dot(wo), bitangent.dot(wo), cos_o];

    // Stretch the view so isotropic GGX becomes a unit-roughness problem.
    let stretched = [alpha * local_wo[0], alpha * local_wo[1], local_wo[2]];
    let inverse_stretched_norm = 1.0
        / (stretched[0] * stretched[0] + stretched[1] * stretched[1] + stretched[2] * stretched[2])
            .sqrt();
    let view = stretched.map(|component| component * inverse_stretched_norm);

    // Orthonormal frame around the stretched view direction.
    let lensq = view[0] * view[0] + view[1] * view[1];
    let frame_1 = if lensq > 0.0 {
        let inverse_len = 1.0 / lensq.sqrt();
        [-view[1] * inverse_len, view[0] * inverse_len, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let frame_2 = [
        view[1] * frame_1[2] - view[2] * frame_1[1],
        view[2] * frame_1[0] - view[0] * frame_1[2],
        view[0] * frame_1[1] - view[1] * frame_1[0],
    ];

    // Uniform disk sample warped onto the visible projected GGX hemisphere.
    let radius = u1.sqrt();
    let phi = 2.0 * PI * u2;
    let disk_1 = radius * det::cos(phi);
    let mut disk_2 = radius * det::sin(phi);
    let blend = 0.5 * (1.0 + view[2]);
    disk_2 = (1.0 - blend) * (1.0 - disk_1 * disk_1).max(0.0).sqrt() + blend * disk_2;
    let projected = (1.0 - disk_1 * disk_1 - disk_2 * disk_2).max(0.0).sqrt();
    let stretched_normal = [
        disk_1 * frame_1[0] + disk_2 * frame_2[0] + projected * view[0],
        disk_1 * frame_1[1] + disk_2 * frame_2[1] + projected * view[1],
        disk_1 * frame_1[2] + disk_2 * frame_2[2] + projected * view[2],
    ];

    // Undo the stretch and normalize back into the macrosurface frame.
    let unstretched = [
        alpha * stretched_normal[0],
        alpha * stretched_normal[1],
        stretched_normal[2].max(0.0),
    ];
    let inverse_normal_len = 1.0
        / (unstretched[0] * unstretched[0]
            + unstretched[1] * unstretched[1]
            + unstretched[2] * unstretched[2])
            .sqrt();
    let local_m = unstretched.map(|component| component * inverse_normal_len);
    let m = to_world_all_sphere(n, local_m);
    let wo_dot_m = wo.dot(m);
    if wo_dot_m <= 0.0 {
        return None;
    }
    let wi = Vec3::new(
        2.0 * wo_dot_m * m.x - wo.x,
        2.0 * wo_dot_m * m.y - wo.y,
        2.0 * wo_dot_m * m.z - wo.z,
    );
    if n.dot(wi) <= 0.0 {
        return None;
    }
    let pdf = ggx_vndf_reflection_pdf(alpha, n, wo, m);
    (pdf > 0.0 && pdf.is_finite()).then_some((wi, pdf))
}

// ---- vector helpers (fs-geom's Vec3 has no cross) ---------------------

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}

/// Encode a film as a linear-sRGB float EXR (channels R, G, B) —
/// byte-exact through fs-img's writer.
///
/// # Errors
/// Propagates [`fs_img::ImgError`] on shape defects.
pub fn film_to_exr(film: &Film) -> Result<Vec<u8>, fs_img::ImgError> {
    let [r, g, b] = film.to_linear_srgb();
    let ch = |name: &str, data: Vec<f32>| fs_img::Channel {
        name: name.to_string(),
        ty: fs_img::PixelType::Float,
        data,
    };
    fs_img::write_exr(
        film.width,
        film.height,
        &[ch("R", r), ch("G", g), ch("B", b)],
    )
}

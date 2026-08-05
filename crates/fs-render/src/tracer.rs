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
//! and arithmetic. The opt-in lighting-v1 extension admits multiple static
//! rectangular area lights and one canonical lat-long environment, ordered and
//! importance-sampled independently of caller construction order. Rectangular
//! lights are also scene geometry so BSDF paths find them (MIS-weighted both
//! ways); materials are Lambertian and
//! GGX (Smith separable G, Schlick Fresnel with the spectral
//! reflectance as F0); no volumetric media (the `volumes` module is
//! separate); no Russian roulette (fixed depth keeps work deterministic).

use crate::animated_instances::{AnimatedGeometryInstance, AnimatedInstanceError};
use crate::camera::{AnimatedCamera, CameraError, CameraExposure, CutSide, LensSample};
use crate::charts::{Hit, Ray, TraceTermination, TriMesh, sphere_trace};
use crate::dielectric::{
    DielectricError, DielectricEvent, DielectricGlass, DielectricSurface,
    evaluate_rough_dielectric, fresnel_dielectric, sample_rough_dielectric,
    sample_smooth_dielectric,
};
use crate::instances::{GeometryInstance, InstanceError, SharedGeometry};
use crate::lighting::{AdmittedLighting, EnvironmentMap, LightSample, LightingError};
use crate::motion::{NormalizedShutterTime, ShutterInterval, TimedRay};
use crate::spectral::{
    LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, cie_x, cie_y, cie_z, xyz_e_to_d65, xyz_to_linear_srgb,
    y_integral,
};
use crate::{balance_heuristic, hero_wavelengths};
use core::mem::size_of;
use fs_alloc::{AllocError, LeaseCharge, LeaseReceipt, LeaseRefusal, OperationMemoryLease};
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
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
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
pub const TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

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
pub const DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

/// Bit-affecting semantics of construction-order-independent multi-light and
/// environment sampling. The legacy one-rectangle/no-environment path remains
/// under [`TRACER_BIT_SEMANTICS_VERSION`] and keeps its frozen stream.
pub const LIGHTING_TRACER_BIT_SEMANTICS_VERSION: u32 = 1;

/// Dedicated Philox counter domain for the two lens coordinates. Lens draws
/// never advance [`PathRng`] and therefore cannot perturb light or BSDF draws.
const CAMERA_LENS_SAMPLE_DOMAIN_V1: u32 = 0x6c65_6e73;

const PI: f64 = core::f64::consts::PI;
/// Hero-wavelength packet width (the bead's 4-wavelength packets).
pub const PACKET: usize = 4;
/// Self-intersection offset along the normal when spawning rays.
const RAY_EPS: f64 = 1e-6;
/// Sphere-trace surface tolerance.
const TRACE_EPS: f64 = 1e-7;
const MAX_MEDIUM_STACK_DEPTH: usize = 64;
const RECT_LIGHT_GEOMETRY_REL_TOLERANCE: f64 = 1.0e-10;

#[derive(Clone, Copy)]
struct PathTime {
    interval: ShutterInterval,
    normalized: NormalizedShutterTime,
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
enum DirectLightTarget {
    Rectangle {
        primitive_index: usize,
        distance_m: f64,
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
    /// A ray missed all geometry while still inside a declared closed medium.
    UnclosedMedium {
        /// Active top boundary at the miss.
        boundary_primitive: usize,
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
            Self::Lighting(error) => write!(formatter, "scene lighting refused: {error}"),
            Self::MediumStackMismatch {
                boundary_primitive,
                active_boundary,
            } => write!(
                formatter,
                "dielectric boundary {boundary_primitive} violated LIFO nesting; active boundary {active_boundary:?}"
            ),
            Self::MediumStackOverflow => formatter.write_str("dielectric medium stack overflow"),
            Self::UnclosedMedium { boundary_primitive } => write!(
                formatter,
                "ray missed while still inside dielectric boundary {boundary_primitive}"
            ),
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
    let exposure = camera.admit_shutter(cx, shutter, cut_side)?;
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
    // Pool seeds control scheduling/placement only. Bind downstream chart,
    // camera, and intersection work to the public render seed while retaining
    // the executor's logical tile/iteration and refusal routing.
    let render_cx = cx.with_stream_seed(settings.seed);
    let (jx, jy, ul) = pixel_dims(settings, sobol, key, pixel, sample)?;
    let ray_time = shutter.map(|interval| PathTime {
        interval,
        normalized: interval.sample_for_stream(settings.seed, u64::from(pixel), u64::from(sample)),
    });
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
        ray_time,
        camera_path,
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
    ray_time: Option<PathTime>,
    camera_path: CameraPath<'_>,
) -> Result<[f64; 3], TracerError> {
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let mut rng = PathRng {
        pixel,
        sample,
        dim: 1,
        key,
    };
    // Hero wavelengths: one stratified draw covers the packet.
    let hero = LAMBDA_MIN + ul * (LAMBDA_MAX - LAMBDA_MIN);
    let lambdas = hero_wavelengths(hero, PACKET, LAMBDA_MIN, LAMBDA_MAX);
    // Camera ray. Keep the legacy branch expression-for-expression compatible;
    // the opt-in cinematic branch owns separate lens dimensions and evaluates
    // at the same absolute time already carried by animated geometry.
    let px = pixel % s.width;
    let py = pixel / s.width;
    let (w, h) = (f64::from(s.width), f64::from(s.height));
    let mut ray = match camera_path {
        CameraPath::Legacy => {
            let ndc_x =
                (2.0 * (f64::from(px) + jx) / w - 1.0) * s.camera_aspect() * scene.camera.half_tan;
            let ndc_y = (1.0 - 2.0 * (f64::from(py) + jy) / h) * scene.camera.half_tan;
            legacy_camera_ray(&scene.camera, ndc_x, ndc_y)
        }
        CameraPath::Cinematic { camera, exposure } => {
            let time = ray_time.ok_or(TracerError::MissingRayTime)?;
            let physical =
                camera.evaluate_exposure(cx, exposure, time.interval.time_at(time.normalized))?;
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
    let mut throughput = [1.0f64; PACKET];
    let mut radiance = [0.0f64; PACKET];
    let mut previous_bsdf: Option<PreviousBsdf> = None;
    let mut prev_origin = ray.origin;
    let mut segment_origin = ray.origin;
    let mut medium_stack = MediumStack::new();
    let mut packet_collapsed = false;
    for depth in 0..s.max_depth {
        cx.checkpoint()?;
        let Some((prim_idx, hit)) = intersect(scene, cx, &ray, ray_time)? else {
            if let Some(active) = medium_stack.last() {
                return Err(TracerError::UnclosedMedium {
                    boundary_primitive: active.boundary_primitive,
                });
            }
            if let Some(environment) = lighting.environment_evaluation(prev_origin, ray.dir) {
                let competing_pdf = previous_bsdf.map(|_| environment.pdf_solid_angle);
                let weight = emissive_hit_weight(s.strategy, previous_bsdf, competing_pdf);
                let (spectrum, scale) = environment.emission;
                for (lane, &lambda) in lambdas.iter().enumerate() {
                    radiance[lane] += throughput[lane] * spectrum.eval(lambda) * scale * weight;
                }
            }
            break;
        };
        attenuate_segment(
            &mut throughput,
            &lambdas,
            &medium_stack,
            segment_origin,
            &hit,
        )?;
        let prim = &scene.primitives[prim_idx];
        let frame = surface_frame(&hit, &ray)?;
        let n = frame.oriented;
        if let Some((spec, scale)) = &prim.emission {
            // MIS weight against NEE for this light, seen from the
            // previous vertex.
            let nee_pdf = if s.strategy == DirectStrategy::Mis && previous_bsdf.is_some() {
                lighting
                    .rect_index_for_primitive(prim_idx)
                    .map(|light_index| {
                        lighting.rect_mixture_pdf(light_index, prev_origin, hit.point)
                    })
            } else {
                None
            };
            let weight = emissive_hit_weight(s.strategy, previous_bsdf, nee_pdf);
            for (k, &l) in lambdas.iter().enumerate() {
                radiance[k] += throughput[k] * spec.eval(l) * scale * weight;
            }
            break; // v1: emitters do not reflect
        }

        let dielectric_boundary = match prim.material {
            Material::Dielectric { glass, .. } => {
                collapse_dispersive_packet(&mut throughput, &mut packet_collapsed, glass);
                Some(boundary_media(
                    prim_idx,
                    glass,
                    frame.entering,
                    &medium_stack,
                )?)
            }
            Material::Lambertian { .. } | Material::Ggx { .. } => None,
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
                            let eta_i = medium_ior(boundary.incident, lambdas[0])?;
                            let eta_t = medium_ior(boundary.transmitted, lambdas[0])?;
                            let evaluation =
                                evaluate_rough_dielectric(n, wo, wi, eta_i, eta_t, alpha)?;
                            if evaluation.value > 0.0 && evaluation.pdf > 0.0 && cos_s > 0.0 {
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
                                        },
                                        Some((index, shadow_hit)),
                                    ) => {
                                        index == primitive_index
                                            && shadow_hit.t > distance_m - 1.0e-4
                                    }
                                    (DirectLightTarget::Environment, None) => true,
                                    (DirectLightTarget::Rectangle { .. }, None)
                                    | (DirectLightTarget::Environment, Some(_)) => false,
                                };
                                if visible {
                                    let pdf_nee = direct.pdf_solid_angle;
                                    let weight = match s.strategy {
                                        DirectStrategy::Mis
                                            if depth + 1 == s.max_depth
                                                && !lighting.is_legacy_compatibility_path() =>
                                        {
                                            1.0
                                        }
                                        DirectStrategy::Mis => {
                                            balance_heuristic(1, pdf_nee, 1, evaluation.pdf)
                                        }
                                        DirectStrategy::NeeOnly => 1.0,
                                        DirectStrategy::BsdfOnly => {
                                            return Err(TracerError::InvalidInput);
                                        }
                                    };
                                    let shadow_medium = match evaluation.event {
                                        DielectricEvent::Reflection => boundary.incident,
                                        DielectricEvent::Transmission => boundary.transmitted,
                                    };
                                    if matches!(direct.target, DirectLightTarget::Environment)
                                        && shadow_medium.is_some()
                                    {
                                        return Err(TracerError::UnclosedMedium {
                                            boundary_primitive: medium_stack
                                                .last()
                                                .map_or(prim_idx, |active| {
                                                    active.boundary_primitive
                                                }),
                                        });
                                    }
                                    let (emission, emission_scale) = direct.emission;
                                    for (lane, &lambda) in lambdas.iter().enumerate() {
                                        if throughput[lane].to_bits() == 0.0_f64.to_bits() {
                                            continue;
                                        }
                                        let eta_i = medium_ior(boundary.incident, lambda)?;
                                        let eta_t = medium_ior(boundary.transmitted, lambda)?;
                                        let value = evaluate_rough_dielectric(
                                            n, wo, wi, eta_i, eta_t, alpha,
                                        )?
                                        .value;
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
                                        radiance[lane] += throughput[lane]
                                            * value
                                            * cos_s
                                            * attenuation
                                            * emission.eval(lambda)
                                            * emission_scale
                                            / pdf_nee
                                            * weight;
                                    }
                                }
                            }
                        }
                    }
                }
                Material::Lambertian { .. } | Material::Ggx { .. } => {
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
                                    },
                                    Some((index, shadow_hit)),
                                ) => index == primitive_index && shadow_hit.t > distance_m - 1.0e-4,
                                (DirectLightTarget::Environment, None) => true,
                                (DirectLightTarget::Rectangle { .. }, None)
                                | (DirectLightTarget::Environment, Some(_)) => false,
                            };
                            if visible {
                                if matches!(direct.target, DirectLightTarget::Environment)
                                    && let Some(active) = medium_stack.last()
                                {
                                    return Err(TracerError::UnclosedMedium {
                                        boundary_primitive: active.boundary_primitive,
                                    });
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
                                            let f = bsdf_eval(&prim.material, n, wo, wi, l);
                                            let attenuation = medium_transmittance(
                                                Some(active.glass),
                                                l,
                                                distance_m,
                                            )?;
                                            radiance[k] += throughput[k]
                                                * f
                                                * cos_s
                                                * attenuation
                                                * espec.eval(l)
                                                * escale
                                                / pdf_nee
                                                * weight;
                                        }
                                    }
                                } else {
                                    // Keep the ambient opaque tracer-v1 arithmetic
                                    // bit-for-bit identical to its frozen path.
                                    for (k, &l) in lambdas.iter().enumerate() {
                                        let f = bsdf_eval(&prim.material, n, wo, wi, l);
                                        radiance[k] +=
                                            throughput[k] * f * cos_s * espec.eval(l) * escale
                                                / pdf_nee
                                                * weight;
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
            Material::Dielectric { surface, .. } => {
                let boundary = dielectric_boundary.ok_or(TracerError::InvalidInput)?;
                let event_sample = if surface.is_delta() {
                    u1
                } else {
                    rng.next2().0
                };
                let Some(sampled) = sample_dielectric_path(
                    n,
                    wo,
                    surface,
                    &boundary,
                    &lambdas,
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
                if sampled.event == DielectricEvent::Transmission {
                    apply_medium_transition(&mut medium_stack, boundary.transition)?;
                }
                previous_bsdf = Some(PreviousBsdf {
                    pdf: sampled.pdf,
                    delta: sampled.delta,
                });
                prev_origin = hit.point;
                segment_origin = hit.point;
                ray = Ray {
                    origin: dielectric_spawn_origin(hit.point, frame.geometric, sampled.direction),
                    dir: sampled.direction,
                };
            }
            Material::Lambertian { .. } | Material::Ggx { .. } => {
                let Some((wi, pdf)) = bsdf_sample(&prim.material, n, wo, u1, u2) else {
                    break;
                };
                let cos_s = n.dot(wi).max(0.0);
                if pdf <= 0.0 || cos_s <= 0.0 {
                    break;
                }
                for (k, &l) in lambdas.iter().enumerate() {
                    throughput[k] *= bsdf_eval(&prim.material, n, wo, wi, l) * cos_s / pdf;
                }
                previous_bsdf = Some(PreviousBsdf { pdf, delta: false });
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
    Ok(xyz)
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
                .get(stack.len().saturating_sub(2))
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

fn collapse_dispersive_packet(
    throughput: &mut [f64; PACKET],
    packet_collapsed: &mut bool,
    glass: DielectricGlass,
) {
    if glass.ior().is_dispersive() && !*packet_collapsed {
        throughput[0] *= PACKET as f64;
        for lane in &mut throughput[1..] {
            *lane = 0.0;
        }
        *packet_collapsed = true;
    }
}

#[allow(clippy::too_many_arguments)]
fn sample_dielectric_path(
    normal: Vec3,
    wo: Vec3,
    surface: DielectricSurface,
    boundary: &BoundaryMedia,
    lambdas: &[f64],
    microfacet_u: f64,
    microfacet_v: f64,
    event_sample: f64,
) -> Result<Option<DielectricPathSample>, TracerError> {
    let hero_eta_i = medium_ior(boundary.incident, lambdas[0])?;
    let hero_eta_t = medium_ior(boundary.transmitted, lambdas[0])?;
    if let Some(alpha) = surface.roughness_alpha() {
        let Some(sample) = sample_rough_dielectric(
            normal,
            wo,
            hero_eta_i,
            hero_eta_t,
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

    let sample = sample_smooth_dielectric(normal, wo, hero_eta_i, hero_eta_t, event_sample)?;
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
    ray_time: Option<PathTime>,
) -> Result<Option<(usize, Hit)>, TracerError> {
    let mut best: Option<(usize, Hit)> = None;
    for (i, prim) in scene.primitives.iter().enumerate() {
        cx.checkpoint()?;
        let hit = match &prim.shape {
            Shape::Mesh(mesh) => mesh.intersect_with_cx(cx, ray)?,
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
                    TraceTermination::Hit => {
                        Some(hit.ok_or(TracerError::BackendFailure(TraceTermination::Hit))?)
                    }
                    termination => return Err(TracerError::BackendFailure(termination)),
                }
            }
            Shape::Instance(instance) => instance
                .intersect(cx, ray, 1e4, TRACE_EPS)?
                .map(|instance_hit| instance_hit.hit),
            Shape::AnimatedInstance(instance) => {
                let time = ray_time.ok_or(TracerError::MissingRayTime)?;
                let timed_ray = TimedRay::at_normalized(*ray, time.interval, time.normalized);
                instance
                    .intersect(cx, &timed_ray, 1e4, TRACE_EPS)?
                    .map(|instance_hit| instance_hit.hit)
            }
        };
        if let Some(h) = hit {
            let replace = match best.as_ref() {
                None => true,
                Some((_, best_hit)) if h.t < best_hit.t => true,
                Some((best_index, best_hit))
                    if h.t.total_cmp(&best_hit.t) == core::cmp::Ordering::Equal =>
                {
                    instance_object_id(&prim.shape)
                        .zip(instance_object_id(&scene.primitives[*best_index].shape))
                        .is_some_and(|(candidate, current)| candidate < current)
                }
                Some(_) => false,
            };
            if replace {
                best = Some((i, h));
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
    fn dispersive_packet_collapses_to_the_weighted_hero_lane_exactly_once() {
        let glass = DielectricGlass::representative_crown();
        assert!(glass.ior().is_dispersive());
        let mut throughput = [1.0, 2.0, 3.0, 4.0];
        let mut collapsed = false;

        collapse_dispersive_packet(&mut throughput, &mut collapsed, glass);
        assert!(collapsed);
        assert_eq!(
            throughput.map(f64::to_bits),
            [4.0, 0.0, 0.0, 0.0].map(f64::to_bits)
        );

        collapse_dispersive_packet(&mut throughput, &mut collapsed, glass);
        assert_eq!(
            throughput.map(f64::to_bits),
            [4.0, 0.0, 0.0, 0.0].map(f64::to_bits),
            "a later dispersive boundary must not apply the estimator weight twice"
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
        Material::Dielectric { .. } => 0.0,
    }
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
            let wom = wo.dot(m);
            if wom <= 0.0 {
                return 0.0;
            }
            ggx_d(*alpha, n.dot(m)) * n.dot(m).max(0.0) / (4.0 * wom)
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
        Material::Ggx { alpha, .. } => {
            // Sample the half-vector from D(m)·cos m (standard GGX NDF
            // sampling; VNDF is a recorded follow-up).
            let a2 = alpha * alpha;
            let cos_m2 = ((1.0 - u1) / (u1 * (a2 - 1.0) + 1.0)).clamp(0.0, 1.0);
            let cos_m = cos_m2.sqrt();
            let sin_m = (1.0 - cos_m2).max(0.0).sqrt();
            let phi = 2.0 * PI * u2;
            let m = to_world(n, [sin_m * det::cos(phi), sin_m * det::sin(phi), cos_m]);
            let wom = wo.dot(m);
            if wom <= 0.0 {
                return None;
            }
            let wi = Vec3::new(
                2.0 * wom * m.x - wo.x,
                2.0 * wom * m.y - wo.y,
                2.0 * wom * m.z - wo.z,
            );
            if n.dot(wi) <= 0.0 {
                return None;
            }
            let pdf = ggx_d(*alpha, n.dot(m)) * n.dot(m).max(0.0) / (4.0 * wom);
            (pdf > 0.0).then_some((wi, pdf))
        }
        Material::Dielectric { .. } => None,
    }
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

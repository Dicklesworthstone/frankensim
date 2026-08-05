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
//! v1 scope (documented, falsifiable): single rectangular area light
//! per scene for NEE; lights are also scene geometry so BSDF paths
//! find them (MIS-weighted both ways); materials are Lambertian and
//! GGX (Smith separable G, Schlick Fresnel with the spectral
//! reflectance as F0); no volumetric media (the `volumes` module is
//! separate); no environment light; no Russian roulette (fixed depth
//! keeps work deterministic).

use crate::animated_instances::{AnimatedGeometryInstance, AnimatedInstanceError};
use crate::camera::{AnimatedCamera, CameraError, CameraExposure, CutSide, LensSample};
use crate::charts::{Hit, Ray, TraceTermination, TriMesh, sphere_trace};
use crate::dielectric::{
    DielectricError, DielectricEvent, DielectricGlass, DielectricSurface,
    evaluate_rough_dielectric, fresnel_dielectric, sample_rough_dielectric,
    sample_smooth_dielectric,
};
use crate::instances::{GeometryInstance, InstanceError};
use crate::motion::{NormalizedShutterTime, ShutterInterval, TimedRay};
use crate::spectral::{
    LAMBDA_MAX, LAMBDA_MIN, LiftedSpectrum, cie_x, cie_y, cie_z, xyz_e_to_d65, xyz_to_linear_srgb,
    y_integral,
};
use crate::{balance_heuristic, hero_wavelengths};
use fs_exec::{Cancelled, Cx};
use fs_geom::{Chart, Point3, Vec3};
use fs_math::det;
use fs_rand::philox::philox4x32_10;
use fs_rand::qmc::Sobol;
use std::collections::BTreeSet;

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

/// The single rectangular area light (v1) used by next-event
/// estimation. The SAME rectangle must also be present as an emissive
/// mesh primitive (index `prim`) so BSDF-sampled paths hit it.
pub struct RectLight {
    /// One corner.
    pub corner: Point3,
    /// First edge.
    pub edge_u: Vec3,
    /// Second edge.
    pub edge_v: Vec3,
    /// Index of the emissive primitive this light corresponds to.
    pub prim: usize,
    /// Emitted radiance spectrum × scale (must match the primitive's).
    pub emission: (LiftedSpectrum, f64),
}

impl RectLight {
    fn area(&self) -> f64 {
        cross(self.edge_u, self.edge_v).norm()
    }

    fn normal(&self) -> Vec3 {
        let n = cross(self.edge_u, self.edge_v);
        n.scale(1.0 / n.norm())
    }
}

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
    /// The NEE light (v1: exactly one).
    pub light: RectLight,
    /// Camera.
    pub camera: Camera,
}

/// Render settings.
#[derive(Debug, Clone, Copy)]
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
    /// A cinematic camera was malformed, evaluated outside its admitted shot,
    /// or cancelled. The nested error retains ranked admission fixes.
    Camera(CameraError),
    /// Validated dielectric evaluation unexpectedly refused.
    Dielectric(DielectricError),
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
            Self::Camera(error) => write!(formatter, "cinematic camera refused: {error}"),
            Self::Dielectric(error) => write!(formatter, "dielectric transport refused: {error}"),
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
        if width == 0 || height == 0 {
            return Err(TracerError::InvalidInput);
        }
        let len = (width as usize)
            .checked_mul(height as usize)
            .ok_or(TracerError::InvalidInput)?;
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

/// Render samples `[from, to)` for every pixel into `film` (progressive
/// accumulation; `film.spp_done` must equal `from`).
///
/// # Panics
/// If the film dimensions or checkpoint do not match. A malformed public XYZ
/// buffer or invalid range is returned as a structured error.
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
    assert_eq!((film.width, film.height), (s.width, s.height), "film shape");
    if film.width == 0 || film.height == 0 || s.width == 0 || s.height == 0 {
        return Err(TracerError::InvalidInput);
    }
    let expected_len = (film.width as usize)
        .checked_mul(film.height as usize)
        .ok_or(TracerError::InvalidInput)?;
    if film.xyz.len() != expected_len {
        return Err(TracerError::InvalidInput);
    }
    if to < from {
        return Err(TracerError::InvalidRange { from, to });
    }
    validate_film_time_mode(film, shutter, s.seed, camera_path)?;
    validate_instance_ids(scene, shutter)?;
    assert_eq!(film.spp_done, from, "progressive checkpoint mismatch");
    cx.checkpoint()?;
    if to == from {
        return Ok(());
    }
    let key = [(s.seed & 0xffff_ffff) as u32, (s.seed >> 32) as u32];
    let sobol = Sobol::scrambled(3, s.seed);
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
                let (jx, jy, ul) = pixel_dims(s, &sobol, key, pixel, sample);
                let ray_time = shutter.map(|interval| PathTime {
                    interval,
                    normalized: interval.sample_for_stream(
                        s.seed,
                        u64::from(pixel),
                        u64::from(sample),
                    ),
                });
                let xyz = trace_path(
                    scene,
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
        film.time_mode = requested_time_mode(shutter, s.seed, camera_path)?;
    }
    Ok(())
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
fn pixel_dims(
    s: &Settings,
    sobol: &Sobol,
    key: [u32; 2],
    pixel: u32,
    sample: u32,
) -> (f64, f64, f64) {
    match s.sampler {
        Sampler::Iid => {
            let a = philox4x32_10([pixel, sample, 0xdead_0001, 0], key);
            (u32_unit(a[0]), u32_unit(a[1]), u32_unit(a[2]))
        }
        Sampler::OwenSobol => {
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
            (
                wrap(pt[0], shift[0]),
                wrap(pt[1], shift[1]),
                wrap(pt[2], shift[2]),
            )
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one integrator, one story
fn trace_path(
    scene: &Scene,
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
    let mut medium_stack = Vec::<MediumEntry>::new();
    let mut packet_collapsed = false;
    for _depth in 0..s.max_depth {
        cx.checkpoint()?;
        let Some((prim_idx, hit)) = intersect(scene, cx, &ray, ray_time)? else {
            if let Some(active) = medium_stack.last() {
                return Err(TracerError::UnclosedMedium {
                    boundary_primitive: active.boundary_primitive,
                });
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
            let nee_pdf = if s.strategy == DirectStrategy::Mis
                && previous_bsdf.is_some()
                && prim_idx == scene.light.prim
            {
                Some({
                    let d = hit.point.delta_from(prev_origin);
                    let d2 = d.dot(d);
                    let cos_l = scene.light.normal().dot(unit(d)).abs();
                    if cos_l > 1e-12 {
                        d2 / (cos_l * scene.light.area())
                    } else {
                        0.0
                    }
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
                        let q = scene
                            .light
                            .corner
                            .offset(scene.light.edge_u.scale(u1))
                            .offset(scene.light.edge_v.scale(u2));
                        let to_light = q.delta_from(hit.point);
                        let d2 = to_light.dot(to_light);
                        if d2 > 0.0 && d2.is_finite() {
                            let dist = d2.sqrt();
                            let wi = to_light.scale(1.0 / dist);
                            let cos_s = n.dot(wi).abs();
                            let cos_l = scene.light.normal().dot(wi).abs();
                            let wo = ray.dir.scale(-1.0);
                            let eta_i = medium_ior(boundary.incident, lambdas[0])?;
                            let eta_t = medium_ior(boundary.transmitted, lambdas[0])?;
                            let evaluation =
                                evaluate_rough_dielectric(n, wo, wi, eta_i, eta_t, alpha)?;
                            if evaluation.value > 0.0
                                && evaluation.pdf > 0.0
                                && cos_s > 0.0
                                && cos_l > 1.0e-9
                            {
                                let shadow = Ray {
                                    origin: dielectric_spawn_origin(hit.point, frame.geometric, wi),
                                    dir: wi,
                                };
                                let visible = match intersect(scene, cx, &shadow, ray_time)? {
                                    Some((index, shadow_hit)) => {
                                        index == scene.light.prim && shadow_hit.t > dist - 1.0e-4
                                    }
                                    None => false,
                                };
                                if visible {
                                    let pdf_nee = d2 / (cos_l * scene.light.area());
                                    let weight = match s.strategy {
                                        DirectStrategy::Mis => {
                                            balance_heuristic(1, pdf_nee, 1, evaluation.pdf)
                                        }
                                        DirectStrategy::NeeOnly => 1.0,
                                        DirectStrategy::BsdfOnly => unreachable!(),
                                    };
                                    let shadow_medium = match evaluation.event {
                                        DielectricEvent::Reflection => boundary.incident,
                                        DielectricEvent::Transmission => boundary.transmitted,
                                    };
                                    let (emission, emission_scale) = &scene.light.emission;
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
                                        let attenuation =
                                            medium_transmittance(shadow_medium, lambda, dist)?;
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
                    let q = scene
                        .light
                        .corner
                        .offset(scene.light.edge_u.scale(u1))
                        .offset(scene.light.edge_v.scale(u2));
                    let to_light = q.delta_from(hit.point);
                    let d2 = to_light.dot(to_light);
                    let dist = d2.sqrt();
                    let wi = to_light.scale(1.0 / dist);
                    let cos_s = n.dot(wi);
                    let cos_l = scene.light.normal().dot(wi).abs();
                    if cos_s > 0.0 && cos_l > 1e-9 {
                        let shadow = Ray {
                            origin: hit.point.offset(n.scale(RAY_EPS)),
                            dir: wi,
                        };
                        let vis = match intersect(scene, cx, &shadow, ray_time)? {
                            Some((i, h)) => i == scene.light.prim && h.t > dist - 1e-4,
                            None => false,
                        };
                        if vis {
                            let pdf_nee = d2 / (cos_l * scene.light.area());
                            let wo = ray.dir.scale(-1.0);
                            let bsdf_pdf = bsdf_pdf(&prim.material, n, wo, wi);
                            let weight = match s.strategy {
                                DirectStrategy::Mis => balance_heuristic(1, pdf_nee, 1, bsdf_pdf),
                                _ => 1.0,
                            };
                            let (espec, escale) = &scene.light.emission;
                            if let Some(active) = medium_stack.last() {
                                for (k, &l) in lambdas.iter().enumerate() {
                                    let f = bsdf_eval(&prim.material, n, wo, wi, l);
                                    let attenuation =
                                        medium_transmittance(Some(active.glass), l, dist)?;
                                    radiance[k] += throughput[k]
                                        * f
                                        * cos_s
                                        * attenuation
                                        * espec.eval(l)
                                        * escale
                                        / pdf_nee
                                        * weight;
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
                    boundary,
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
    stack: &[MediumEntry],
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
    stack: &mut Vec<MediumEntry>,
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
            stack.push(entry);
        }
        MediumTransition::Exit { boundary_primitive } => {
            if stack.last().map(|entry| entry.boundary_primitive) != Some(boundary_primitive) {
                return Err(TracerError::MediumStackMismatch {
                    boundary_primitive,
                    active_boundary: stack.last().map(|entry| entry.boundary_primitive),
                });
            }
            let _ = stack.pop();
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
    medium_stack: &[MediumEntry],
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
    boundary: BoundaryMedia,
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
            for (lane, &lambda) in lambdas.iter().enumerate() {
                let eta_i = medium_ior(boundary.incident, lambda)?;
                let eta_t = medium_ior(boundary.transmitted, lambda)?;
                let evaluation =
                    evaluate_rough_dielectric(normal, wo, sample.direction, eta_i, eta_t, alpha)?;
                weights[lane] = evaluation.value * normal.dot(sample.direction).abs() / sample.pdf;
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
    for (lane, &lambda) in lambdas.iter().enumerate() {
        let eta_i = medium_ior(boundary.incident, lambda)?;
        let eta_t = medium_ior(boundary.transmitted, lambda)?;
        let fresnel = fresnel_dielectric(normal.dot(wo), eta_i, eta_t)?;
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

fn dielectric_spawn_origin(point: Point3, geometric_normal: Vec3, direction: Vec3) -> Point3 {
    let scale = point.x.abs().max(point.y.abs()).max(point.z.abs()).max(1.0);
    let side = if geometric_normal.dot(direction) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    point.offset(geometric_normal.scale(side * RAY_EPS * scale))
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

fn validate_instance_ids(
    scene: &Scene,
    shutter: Option<ShutterInterval>,
) -> Result<(), TracerError> {
    let light_primitive = scene
        .primitives
        .get(scene.light.prim)
        .ok_or(TracerError::InvalidInput)?;
    if matches!(&light_primitive.shape, Shape::AnimatedInstance(_)) {
        return Err(TracerError::AnimatedLightUnsupported);
    }
    let mut object_ids = BTreeSet::new();
    for primitive in &scene.primitives {
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
    Ok(())
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

fn validate_film_time_mode(
    film: &Film,
    shutter: Option<ShutterInterval>,
    stream_identity: u64,
    camera_path: CameraPath<'_>,
) -> Result<(), TracerError> {
    let requested = requested_time_mode(shutter, stream_identity, camera_path)?;
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

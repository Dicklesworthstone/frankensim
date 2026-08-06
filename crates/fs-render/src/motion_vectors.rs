//! Stable primary-hit correspondence and exact geometric motion vectors.
//!
//! Motion here means the screen displacement of one admitted local material
//! point under explicit rigid object and camera transforms. It is a rendering
//! derivative, not a measured mechanical velocity. Projection uses the camera
//! optical centre so aperture samples cannot masquerade as object motion.

use core::fmt;

use fs_blake3::ContentHash;
use fs_exec::{Cancelled, Cx};
use fs_geom::{Point3, Vec3};

use crate::animated_instances::{AnimatedGeometryInstance, AnimatedInstanceError};
use crate::camera::{
    AnimatedCamera, CameraError, CutSide, OpticalCenterProjection, PhysicalCamera,
};
use crate::instances::{GeometryInstance, InstanceHit, InstanceSurfaceFeature, RigidTransform};

/// Bit-affecting semantics for feature selection, projection, and motion-vector
/// direction. Increment when any result bit or categorical outcome can change.
pub const MOTION_VECTOR_SEMANTICS_VERSION: u32 = 1;

const BARYCENTRIC_SUM_TOLERANCE: f64 = 16.0 * f64::EPSILON;

/// Validated raster extent shared by beauty and motion/AOV consumers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RasterSize {
    width: u32,
    height: u32,
}

impl RasterSize {
    /// Admit nonzero dimensions.
    pub fn try_new(width: u32, height: u32) -> Result<Self, MotionVectorError> {
        if width == 0 || height == 0 {
            return Err(MotionVectorError::InvalidRasterSize);
        }
        Ok(Self { width, height })
    }

    /// Width in pixels.
    #[must_use]
    pub const fn width(self) -> u32 {
        self.width
    }

    /// Height in pixels.
    #[must_use]
    pub const fn height(self) -> u32 {
        self.height
    }

    /// Exact render aspect ratio `width / height`.
    #[must_use]
    pub fn aspect_ratio(self) -> f64 {
        f64::from(self.width) / f64::from(self.height)
    }

    /// Row-major index when `(x, y)` is inside this raster.
    #[must_use]
    pub fn linear_index(self, x: u32, y: u32) -> Option<u64> {
        (x < self.width && y < self.height)
            .then(|| u64::from(y) * u64::from(self.width) + u64::from(x))
    }

    /// Pixel count without allocating an AOV.
    #[must_use]
    pub const fn pixel_count(self) -> u64 {
        self.width as u64 * self.height as u64
    }
}

/// Backend feature component of a stable primary-hit identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StableFeatureIdentity {
    /// Original triangle index in an immutable ordered mesh artifact.
    MeshTriangle(u32),
    /// The chart supplies no admitted stable feature parameter.
    ChartUnavailable,
}

/// Stable categorical identity retained with a primary surface hit.
///
/// The frame/pose identity is intentionally absent: it changes while the same
/// object, geometry, material, and feature move between frames.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StableHitIdentity {
    object_id: u64,
    geometry_identity: ContentHash,
    material_identity: ContentHash,
    feature: StableFeatureIdentity,
}

impl StableHitIdentity {
    /// Stable object identity.
    #[must_use]
    pub const fn object_id(self) -> u64 {
        self.object_id
    }

    /// Caller-supplied immutable geometry identity.
    #[must_use]
    pub const fn geometry_identity(self) -> ContentHash {
        self.geometry_identity
    }

    /// Caller-supplied material content identity.
    #[must_use]
    pub const fn material_identity(self) -> ContentHash {
        self.material_identity
    }

    /// Stable backend feature, or an explicit chart refusal.
    #[must_use]
    pub const fn feature(self) -> StableFeatureIdentity {
        self.feature
    }
}

/// Reconstructible point on one immutable mesh triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeshSurfacePoint {
    triangle_index: u32,
    barycentric: [f64; 3],
    local_point: Point3,
}

impl MeshSurfacePoint {
    /// Original triangle index.
    #[must_use]
    pub const fn triangle_index(self) -> u32 {
        self.triangle_index
    }

    /// Barycentric coordinates ordered like the triangle's vertex indices.
    #[must_use]
    pub const fn barycentric(self) -> [f64; 3] {
        self.barycentric
    }

    /// Accepted local-space hit point. The barycentrics are retained as the
    /// stronger feature-local witness for target-frame validation.
    #[must_use]
    pub const fn local_point(self) -> Point3 {
        self.local_point
    }
}

/// Why a visible primary hit cannot be mapped between frames.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CorrespondenceUnavailable {
    /// Generic charts currently expose no stable surface parameter. A local
    /// point is not silently promoted to a material coordinate.
    ChartHasNoStableParameter,
}

/// Backend-specific local correspondence retained at the accepted primary hit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SurfaceCorrespondence {
    /// Exact rigid-mesh correspondence.
    Mesh(MeshSurfacePoint),
    /// Fail-closed no-correspondence outcome.
    Unavailable(CorrespondenceUnavailable),
}

/// Primary surface data shared by motion, depth, and normal AOV extraction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrimarySurfaceSample {
    identity: StableHitIdentity,
    source_frame_identity: ContentHash,
    correspondence: SurfaceCorrespondence,
    geometric_normal_local: Vec3,
    shading_normal_local: Option<Vec3>,
}

impl PrimarySurfaceSample {
    /// Validate and retain correspondence directly from the accepted instance
    /// hit. `material_identity` is supplied by scene composition because a
    /// material value alone is not an authoritative asset identity.
    pub fn try_from_instance_hit(
        hit: &InstanceHit,
        material_identity: ContentHash,
    ) -> Result<Self, MotionVectorError> {
        if is_zero_hash(material_identity) {
            return Err(MotionVectorError::InvalidMaterialIdentity);
        }
        let local_point = hit.local_hit.point;
        ensure_point_finite(local_point)?;
        let geometric_normal_local = hit
            .local_hit
            .normal
            .ok_or(MotionVectorError::InvalidSurfaceWitness)?;
        ensure_direction_finite(geometric_normal_local)?;
        if let Some(normal) = hit.local_hit.shading_normal {
            ensure_direction_finite(normal)?;
        }

        let (feature, correspondence) = match hit.surface_feature {
            InstanceSurfaceFeature::ChartUnavailable => (
                StableFeatureIdentity::ChartUnavailable,
                SurfaceCorrespondence::Unavailable(
                    CorrespondenceUnavailable::ChartHasNoStableParameter,
                ),
            ),
            InstanceSurfaceFeature::MeshTriangle {
                triangle_index,
                barycentric,
            } => {
                validate_barycentric(barycentric)?;
                (
                    StableFeatureIdentity::MeshTriangle(triangle_index),
                    SurfaceCorrespondence::Mesh(MeshSurfacePoint {
                        triangle_index,
                        barycentric,
                        local_point,
                    }),
                )
            }
        };

        Ok(Self {
            identity: StableHitIdentity {
                object_id: hit.object_id,
                geometry_identity: hit.geometry_identity,
                material_identity,
                feature,
            },
            source_frame_identity: hit.frame_identity,
            correspondence,
            geometric_normal_local,
            shading_normal_local: hit.local_hit.shading_normal,
        })
    }

    /// Stable categorical identity.
    #[must_use]
    pub const fn identity(self) -> StableHitIdentity {
        self.identity
    }

    /// Object/geometry/pose identity at the accepted primary hit.
    #[must_use]
    pub const fn source_frame_identity(self) -> ContentHash {
        self.source_frame_identity
    }

    /// Reconstructible correspondence or explicit refusal.
    #[must_use]
    pub const fn correspondence(self) -> SurfaceCorrespondence {
        self.correspondence
    }

    /// Local geometric normal.
    #[must_use]
    pub const fn geometric_normal_local(self) -> Vec3 {
        self.geometric_normal_local
    }

    /// Local shading normal, kept distinct from geometry.
    #[must_use]
    pub const fn shading_normal_local(self) -> Option<Vec3> {
        self.shading_normal_local
    }
}

/// Exact camera and object state used at one reference frame.
#[derive(Clone, Debug, PartialEq)]
pub struct MotionFrame {
    absolute_time_s: f64,
    shot_id: u64,
    object_id: u64,
    geometry_identity: ContentHash,
    frame_identity: ContentHash,
    object_to_world: RigidTransform,
    camera: PhysicalCamera,
}

impl MotionFrame {
    /// Bind a static placed instance to an already evaluated camera.
    pub fn from_instance(
        absolute_time_s: f64,
        shot_id: u64,
        camera: PhysicalCamera,
        instance: &GeometryInstance,
    ) -> Result<Self, MotionVectorError> {
        if !absolute_time_s.is_finite() || shot_id == 0 {
            return Err(MotionVectorError::InvalidFrame);
        }
        Ok(Self {
            absolute_time_s: canonical_zero(absolute_time_s),
            shot_id,
            object_id: instance.object_id(),
            geometry_identity: instance.geometry_identity(),
            frame_identity: instance.frame_identity(),
            object_to_world: instance.transform(),
            camera,
        })
    }

    /// Evaluate camera and object at exactly the same absolute time. Camera
    /// shot identity is retained so hard cuts become categorical refusals.
    pub fn from_animated(
        cx: &Cx<'_>,
        absolute_time_s: f64,
        cut_side: CutSide,
        camera: &AnimatedCamera,
        instance: &AnimatedGeometryInstance,
    ) -> Result<Self, MotionVectorError> {
        cx.checkpoint()?;
        let evaluated_camera = camera.evaluate_with_shot(cx, absolute_time_s, cut_side)?;
        let evaluated_instance = instance.instance_at(cx, absolute_time_s)?;
        let frame = Self::from_instance(
            absolute_time_s,
            evaluated_camera.shot_id(),
            evaluated_camera.into_camera(),
            &evaluated_instance,
        )?;
        cx.checkpoint()?;
        Ok(frame)
    }

    /// Absolute reference time in seconds.
    #[must_use]
    pub const fn absolute_time_s(&self) -> f64 {
        self.absolute_time_s
    }

    /// Owning continuous-shot identity.
    #[must_use]
    pub const fn shot_id(&self) -> u64 {
        self.shot_id
    }

    /// Stable object identity.
    #[must_use]
    pub const fn object_id(&self) -> u64 {
        self.object_id
    }

    /// Immutable geometry identity.
    #[must_use]
    pub const fn geometry_identity(&self) -> ContentHash {
        self.geometry_identity
    }

    /// Object/geometry/pose identity at this reference time.
    #[must_use]
    pub const fn frame_identity(&self) -> ContentHash {
        self.frame_identity
    }

    /// Body-to-world transform.
    #[must_use]
    pub const fn object_to_world(&self) -> RigidTransform {
        self.object_to_world
    }

    /// Exact physical camera.
    #[must_use]
    pub const fn camera(&self) -> &PhysicalCamera {
        &self.camera
    }
}

/// Finite perspective projection in both NDC and row-major pixel coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RasterProjection {
    ndc_xy: [f64; 2],
    pixel_xy: [f64; 2],
    depth_m: f64,
    in_frame: bool,
}

impl RasterProjection {
    /// NDC coordinate, `+x` right and `+y` up.
    #[must_use]
    pub const fn ndc_xy(self) -> [f64; 2] {
        self.ndc_xy
    }

    /// Continuous pixel-edge coordinate, `+x` right and `+y` down. Pixel
    /// centres occupy `(x + 0.5, y + 0.5)` under the tracer convention.
    #[must_use]
    pub const fn pixel_xy(self) -> [f64; 2] {
        self.pixel_xy
    }

    /// Positive axial camera depth in metres.
    #[must_use]
    pub const fn depth_m(self) -> f64 {
        self.depth_m
    }

    /// Whether the continuous coordinate lies in the half-open raster extent.
    #[must_use]
    pub const fn in_frame(self) -> bool {
        self.in_frame
    }
}

/// Previous/next projection relative to the accepted current primary hit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MotionEndpoint {
    /// Finite projection. Off-screen endpoints retain their exact displacement
    /// but are rejected by reprojection validation.
    Projected {
        /// Target-frame projection.
        target: RasterProjection,
        /// `target - current` in NDC, with `+y` up.
        displacement_ndc: [f64; 2],
        /// `target - current` in raster pixels, with `+y` down.
        displacement_pixels: [f64; 2],
    },
    /// Reference frames belong to different continuous shots.
    CameraCut {
        /// Shot containing the current primary hit.
        current_shot_id: u64,
        /// Shot containing the target reference frame.
        target_shot_id: u64,
    },
    /// The corresponding point lies on or behind the target lens plane.
    BehindCamera {
        /// Nonpositive signed axial depth in metres.
        signed_depth_m: f64,
    },
}

/// Available three-frame motion and aligned normal data for one primary hit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionVectorSample {
    /// Stable categorical identity for AOVs and target validation.
    pub identity: StableHitIdentity,
    /// Current center-projection in the exact beauty raster convention.
    pub current: RasterProjection,
    /// Current-to-previous projection.
    pub previous: MotionEndpoint,
    /// Current-to-next projection.
    pub next: MotionEndpoint,
    /// Current world-space geometric normal.
    pub geometric_normal_world: Vec3,
    /// Current world-space shading normal, kept separate.
    pub shading_normal_world: Option<Vec3>,
}

/// Result of attempting motion correspondence for one accepted primary hit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MotionVectorComputation {
    /// Exact rigid local correspondence was available.
    Available(MotionVectorSample),
    /// Object/material AOV identity is still usable, but motion is not.
    Unavailable {
        /// Stable categorical identity retained from the hit.
        identity: StableHitIdentity,
        /// Explicit refusal reason.
        reason: CorrespondenceUnavailable,
    },
}

/// Compute previous/current/next geometric motion from one accepted primary
/// surface sample.
///
/// Reference times must be nondecreasing. The current frame must exactly match
/// the source pose identity of the primary hit; all three frames must name the
/// same object and immutable geometry. A hard camera cut is reported per
/// endpoint instead of producing a numeric vector.
pub fn compute_motion_vectors(
    sample: PrimarySurfaceSample,
    previous: &MotionFrame,
    current: &MotionFrame,
    next: &MotionFrame,
    raster: RasterSize,
) -> Result<MotionVectorComputation, MotionVectorError> {
    if previous.absolute_time_s > current.absolute_time_s
        || current.absolute_time_s > next.absolute_time_s
    {
        return Err(MotionVectorError::NonMonotonicReferenceTimes);
    }
    for frame in [previous, current, next] {
        if frame.object_id != sample.identity.object_id
            || frame.geometry_identity != sample.identity.geometry_identity
        {
            return Err(MotionVectorError::ForeignObjectFrame);
        }
    }
    if current.frame_identity != sample.source_frame_identity {
        return Err(MotionVectorError::ForeignCurrentPose);
    }

    let surface = match sample.correspondence {
        SurfaceCorrespondence::Mesh(surface) => surface,
        SurfaceCorrespondence::Unavailable(reason) => {
            return Ok(MotionVectorComputation::Unavailable {
                identity: sample.identity,
                reason,
            });
        }
    };

    let current_world = current.object_to_world.transform_point(surface.local_point);
    let current_projection = project_world(current.camera(), current_world, raster)?;
    let ProjectedWorldPoint::InFront(current_projection) = current_projection else {
        return Err(MotionVectorError::CurrentPointBehindCamera);
    };

    let previous_endpoint =
        project_endpoint(surface, previous, current, current_projection, raster)?;
    let next_endpoint = project_endpoint(surface, next, current, current_projection, raster)?;
    let geometric_normal_world = current
        .object_to_world
        .transform_vector(sample.geometric_normal_local);
    ensure_direction_finite(geometric_normal_world)?;
    let shading_normal_world = sample
        .shading_normal_local
        .map(|normal| current.object_to_world.transform_vector(normal));
    if let Some(normal) = shading_normal_world {
        ensure_direction_finite(normal)?;
    }

    Ok(MotionVectorComputation::Available(MotionVectorSample {
        identity: sample.identity,
        current: current_projection,
        previous: previous_endpoint,
        next: next_endpoint,
        geometric_normal_world,
        shading_normal_world,
    }))
}

fn project_endpoint(
    surface: MeshSurfacePoint,
    target_frame: &MotionFrame,
    current_frame: &MotionFrame,
    current_projection: RasterProjection,
    raster: RasterSize,
) -> Result<MotionEndpoint, MotionVectorError> {
    if target_frame.shot_id != current_frame.shot_id {
        return Ok(MotionEndpoint::CameraCut {
            current_shot_id: current_frame.shot_id,
            target_shot_id: target_frame.shot_id,
        });
    }
    let target_world = target_frame
        .object_to_world
        .transform_point(surface.local_point);
    match project_world(target_frame.camera(), target_world, raster)? {
        ProjectedWorldPoint::BehindCamera { signed_depth_m } => {
            Ok(MotionEndpoint::BehindCamera { signed_depth_m })
        }
        ProjectedWorldPoint::InFront(target) => Ok(MotionEndpoint::Projected {
            target,
            displacement_ndc: subtract2(target.ndc_xy, current_projection.ndc_xy),
            displacement_pixels: subtract2(target.pixel_xy, current_projection.pixel_xy),
        }),
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ProjectedWorldPoint {
    InFront(RasterProjection),
    BehindCamera { signed_depth_m: f64 },
}

fn project_world(
    camera: &PhysicalCamera,
    world_point: Point3,
    raster: RasterSize,
) -> Result<ProjectedWorldPoint, MotionVectorError> {
    match camera.project_from_optical_center(world_point, raster.aspect_ratio())? {
        OpticalCenterProjection::BehindCamera { signed_depth_m } => {
            Ok(ProjectedWorldPoint::BehindCamera { signed_depth_m })
        }
        OpticalCenterProjection::InFront { ndc_xy, depth_m } => {
            let width = f64::from(raster.width);
            let height = f64::from(raster.height);
            let pixel_xy = [
                (ndc_xy[0] + 1.0) * 0.5 * width,
                (1.0 - ndc_xy[1]) * 0.5 * height,
            ];
            if pixel_xy.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(MotionVectorError::InvalidProjection);
            }
            let in_frame = pixel_xy[0] >= 0.0
                && pixel_xy[0] < width
                && pixel_xy[1] >= 0.0
                && pixel_xy[1] < height;
            Ok(ProjectedWorldPoint::InFront(RasterProjection {
                ndc_xy,
                pixel_xy,
                depth_m,
                in_frame,
            }))
        }
    }
}

/// Target-frame frontmost observation at the reprojected coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReprojectionObservation {
    /// No target-frame primary surface was present.
    Background,
    /// One frontmost target-frame surface sample.
    Surface {
        /// Stable object/geometry/material/feature identity.
        identity: StableHitIdentity,
        /// Positive target-camera axial depth in metres.
        depth_m: f64,
        /// Target mesh barycentrics when the feature supports them.
        barycentric: Option<[f64; 3]>,
    },
}

impl ReprojectionObservation {
    /// Construct a finite positive-depth surface observation.
    pub fn try_surface(
        identity: StableHitIdentity,
        depth_m: f64,
        barycentric: Option<[f64; 3]>,
    ) -> Result<Self, MotionVectorError> {
        if !depth_m.is_finite() || depth_m <= 0.0 {
            return Err(MotionVectorError::InvalidObservation);
        }
        if let Some(barycentric) = barycentric {
            validate_barycentric(barycentric)?;
        }
        Ok(Self::Surface {
            identity,
            depth_m,
            barycentric,
        })
    }
}

/// Explicit tolerances for target-frame depth and local-coordinate agreement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReprojectionTolerance {
    depth_absolute_m: f64,
    depth_relative: f64,
    barycentric_absolute: f64,
}

impl ReprojectionTolerance {
    /// Admit finite nonnegative tolerances.
    pub fn try_new(
        depth_absolute_m: f64,
        depth_relative: f64,
        barycentric_absolute: f64,
    ) -> Result<Self, MotionVectorError> {
        if [depth_absolute_m, depth_relative, barycentric_absolute]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(MotionVectorError::InvalidReprojectionTolerance);
        }
        Ok(Self {
            depth_absolute_m,
            depth_relative,
            barycentric_absolute,
        })
    }

    fn depth_band(self, expected: f64, observed: f64) -> f64 {
        self.depth_absolute_m + self.depth_relative * expected.abs().max(observed.abs())
    }
}

/// Visibility/correspondence verdict for one target endpoint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReprojectionStatus {
    /// Identity, depth, and feature-local coordinates agree.
    VisibleAtTarget,
    /// A different frontmost surface is measurably nearer.
    OccludedAtTarget,
    /// The expected surface is absent or a different surface is farther away.
    DisoccludedAtTarget,
    /// The finite target projection lies outside the raster.
    TargetOffScreen,
    /// The target point is on or behind the lens plane.
    BehindTargetCamera,
    /// The endpoint crosses a hard camera cut.
    CameraCut,
    /// Different categorical identities occurred at indistinguishable depth.
    IdentityMismatch,
    /// The same identity appeared at an incompatible depth.
    DepthMismatch,
    /// The same feature/depth appeared at incompatible local coordinates.
    TopologyAmbiguous,
    /// The target AOV omitted the mesh-local witness needed for validation.
    MissingSurfaceWitness,
}

/// Validate one projected endpoint against the frontmost target-frame AOV
/// observation. Projection alone never claims visibility.
#[must_use]
pub fn validate_reprojection(
    sample: PrimarySurfaceSample,
    endpoint: MotionEndpoint,
    observation: ReprojectionObservation,
    tolerance: ReprojectionTolerance,
) -> ReprojectionStatus {
    let MotionEndpoint::Projected { target, .. } = endpoint else {
        return match endpoint {
            MotionEndpoint::CameraCut { .. } => ReprojectionStatus::CameraCut,
            MotionEndpoint::BehindCamera { .. } => ReprojectionStatus::BehindTargetCamera,
            MotionEndpoint::Projected { .. } => unreachable!(),
        };
    };
    if !target.in_frame {
        return ReprojectionStatus::TargetOffScreen;
    }
    let ReprojectionObservation::Surface {
        identity,
        depth_m,
        barycentric,
    } = observation
    else {
        return ReprojectionStatus::DisoccludedAtTarget;
    };

    let depth_band = tolerance.depth_band(target.depth_m, depth_m);
    if identity != sample.identity {
        if depth_m + depth_band < target.depth_m {
            return ReprojectionStatus::OccludedAtTarget;
        }
        if target.depth_m + depth_band < depth_m {
            return ReprojectionStatus::DisoccludedAtTarget;
        }
        return ReprojectionStatus::IdentityMismatch;
    }
    if (depth_m - target.depth_m).abs() > depth_band {
        return ReprojectionStatus::DepthMismatch;
    }
    let SurfaceCorrespondence::Mesh(expected) = sample.correspondence else {
        return ReprojectionStatus::MissingSurfaceWitness;
    };
    let Some(observed) = barycentric else {
        return ReprojectionStatus::MissingSurfaceWitness;
    };
    if expected
        .barycentric
        .into_iter()
        .zip(observed)
        .any(|(left, right)| (left - right).abs() > tolerance.barycentric_absolute)
    {
        return ReprojectionStatus::TopologyAmbiguous;
    }
    ReprojectionStatus::VisibleAtTarget
}

/// Fail-closed construction or projection error.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MotionVectorError {
    /// Width and height must both be nonzero.
    InvalidRasterSize,
    /// Frame time or shot identity is invalid.
    InvalidFrame,
    /// Material identity is the reserved all-zero digest.
    InvalidMaterialIdentity,
    /// Local point, normal, or barycentric witness is malformed.
    InvalidSurfaceWitness,
    /// Reference times are not ordered previous <= current <= next.
    NonMonotonicReferenceTimes,
    /// A reference frame names another object or immutable geometry.
    ForeignObjectFrame,
    /// The current frame pose does not match the accepted primary hit.
    ForeignCurrentPose,
    /// The accepted current primary point does not project in front of camera.
    CurrentPointBehindCamera,
    /// Projection arithmetic did not remain finite.
    InvalidProjection,
    /// Reprojection tolerances are negative or non-finite.
    InvalidReprojectionTolerance,
    /// Target observation depth or local witness is malformed.
    InvalidObservation,
    /// Camera admission/evaluation/projection refused.
    Camera(CameraError),
    /// Animated instance evaluation refused.
    AnimatedInstance(AnimatedInstanceError),
    /// Execution scope requested cancellation.
    Cancelled,
}

impl From<CameraError> for MotionVectorError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}

impl From<AnimatedInstanceError> for MotionVectorError {
    fn from(error: AnimatedInstanceError) -> Self {
        Self::AnimatedInstance(error)
    }
}

impl From<Cancelled> for MotionVectorError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl fmt::Display for MotionVectorError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "motion-vector evaluation refused: {self:?}")
    }
}

impl core::error::Error for MotionVectorError {}

fn validate_barycentric(barycentric: [f64; 3]) -> Result<(), MotionVectorError> {
    if barycentric
        .iter()
        .any(|weight| !weight.is_finite() || *weight < 0.0 || *weight > 1.0)
        || (barycentric.into_iter().sum::<f64>() - 1.0).abs() > BARYCENTRIC_SUM_TOLERANCE
    {
        return Err(MotionVectorError::InvalidSurfaceWitness);
    }
    Ok(())
}

fn ensure_point_finite(point: Point3) -> Result<(), MotionVectorError> {
    if [point.x, point.y, point.z]
        .iter()
        .all(|value| value.is_finite())
    {
        Ok(())
    } else {
        Err(MotionVectorError::InvalidSurfaceWitness)
    }
}

fn ensure_direction_finite(vector: Vec3) -> Result<(), MotionVectorError> {
    let components = [vector.x, vector.y, vector.z];
    let scale = components.into_iter().map(f64::abs).fold(0.0, f64::max);
    if components.iter().all(|value| value.is_finite()) && scale > 0.0 {
        Ok(())
    } else {
        Err(MotionVectorError::InvalidSurfaceWitness)
    }
}

fn is_zero_hash(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn subtract2(left: [f64; 2], right: [f64; 2]) -> [f64; 2] {
    [
        canonical_zero(left[0] - right[0]),
        canonical_zero(left[1] - right[1]),
    ]
}

#[cfg(test)]
mod tests {
    use asupersync::types::Budget;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::fixtures::SphereChart;

    use super::*;
    use crate::animated_instances::{RigidTransformTrajectory, TransformKeyframe};
    use crate::camera::{Aperture, CameraKeyframe, CameraProjection, CameraShot};
    use crate::charts::{Ray, TriMesh};
    use crate::instances::SharedGeometry;

    const TOL: f64 = 2.0e-12;

    fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x6d6f_7469_6f6e,
                    kernel_id: 19,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn identity(label: &str) -> ContentHash {
        hash_domain("org.frankensim.test.motion-vector", label.as_bytes())
    }

    fn triangle() -> TriMesh {
        TriMesh::new(
            vec![[0.0, -1.0, 0.0], [2.0, -1.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
    }

    fn transform_z(angle: f64, translation: [f64; 3]) -> RigidTransform {
        let half = 0.5 * angle;
        RigidTransform::try_new([0.0, 0.0, half.sin(), half.cos()], translation).unwrap()
    }

    fn camera_at(x: f64, aperture_radius: f64) -> PhysicalCamera {
        PhysicalCamera::try_look_at(
            Point3::new(x, 0.0, 0.0),
            Point3::new(x, 0.0, -1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.5).unwrap(),
            5.0,
            Aperture::try_circular(aperture_radius).unwrap(),
        )
        .unwrap()
    }

    fn instance(transform: RigidTransform) -> GeometryInstance {
        GeometryInstance::try_new(
            41,
            identity("mesh"),
            SharedGeometry::mesh(triangle()),
            transform,
        )
        .unwrap()
    }

    fn hit_sample(cx: &Cx<'_>, instance: &GeometryInstance) -> PrimarySurfaceSample {
        let target = instance
            .transform()
            .transform_point(Point3::new(1.0, 0.0, 0.0));
        let origin = camera_at(0.0, 0.0).eye();
        let ray = Ray {
            origin,
            dir: target.delta_from(origin),
        };
        let hit = instance.intersect(cx, &ray, 2.0, 1.0e-9).unwrap().unwrap();
        PrimarySurfaceSample::try_from_instance_hit(&hit, identity("steel")).unwrap()
    }

    fn frame(
        time: f64,
        shot_id: u64,
        camera: PhysicalCamera,
        transform: RigidTransform,
    ) -> MotionFrame {
        MotionFrame::from_instance(time, shot_id, camera, &instance(transform)).unwrap()
    }

    fn available(result: MotionVectorComputation) -> MotionVectorSample {
        let MotionVectorComputation::Available(sample) = result else {
            panic!("expected available motion");
        };
        sample
    }

    fn projected(endpoint: MotionEndpoint) -> (RasterProjection, [f64; 2], [f64; 2]) {
        let MotionEndpoint::Projected {
            target,
            displacement_ndc,
            displacement_pixels,
        } = endpoint
        else {
            panic!("expected projected endpoint");
        };
        (target, displacement_ndc, displacement_pixels)
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= TOL * actual.abs().max(expected.abs()).max(1.0),
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn mesh_hit_retains_barycentrics_and_lowest_triangle_tie() {
        let mesh = TriMesh::new(
            vec![[0.0, -1.0, 0.0], [2.0, -1.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[2, 1, 0], [0, 1, 2]],
        );
        let ray = Ray {
            origin: Point3::new(1.0, 0.0, 2.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let bvh = mesh.intersect_surface(&ray).unwrap();
        let brute = mesh.intersect_surface_bruteforce(&ray).unwrap();
        assert_eq!(bvh.triangle_index, 0);
        assert_eq!(brute.triangle_index, 0);
        assert_eq!(bvh.barycentric, brute.barycentric);
        assert_close(bvh.barycentric.into_iter().sum(), 1.0);
        let indices = mesh.triangles[bvh.triangle_index as usize];
        let mut reconstructed = [0.0; 3];
        for (weight, vertex_index) in bvh.barycentric.into_iter().zip(indices) {
            for (coordinate, vertex) in reconstructed
                .iter_mut()
                .zip(mesh.vertices[vertex_index as usize])
            {
                *coordinate += weight * vertex;
            }
        }
        assert_close(reconstructed[0], bvh.hit.point.x);
        assert_close(reconstructed[1], bvh.hit.point.y);
        assert_close(reconstructed[2], bvh.hit.point.z);
    }

    #[test]
    fn translating_and_rotating_object_match_analytic_projection() {
        with_cx(|cx| {
            let current_transform = transform_z(0.0, [0.0, 0.0, -5.0]);
            let sample = hit_sample(cx, &instance(current_transform));
            let raster = RasterSize::try_new(200, 100).unwrap();
            let previous = frame(
                -1.0,
                7,
                camera_at(0.0, 0.0),
                transform_z(0.0, [-1.0, 0.0, -5.0]),
            );
            let current = frame(0.0, 7, camera_at(0.0, 0.0), current_transform);
            let next = frame(
                1.0,
                7,
                camera_at(0.0, 0.0),
                transform_z(core::f64::consts::FRAC_PI_2, [0.0, 0.0, -5.0]),
            );
            let motion = available(
                compute_motion_vectors(sample, &previous, &current, &next, raster).unwrap(),
            );
            let (_, previous_ndc, previous_pixels) = projected(motion.previous);
            let (_, next_ndc, next_pixels) = projected(motion.next);
            // Current local point is (1,0,0). Previous translation places it
            // at x=0; next quarter-turn places it at y=1.
            assert_close(previous_ndc[0], -0.2);
            assert_close(previous_ndc[1], 0.0);
            assert_close(previous_pixels[0], -20.0);
            assert_close(previous_pixels[1], 0.0);
            assert_close(next_ndc[0], -0.2);
            assert_close(next_ndc[1], 0.4);
            assert_close(next_pixels[0], -20.0);
            assert_close(next_pixels[1], -20.0);
        });
    }

    #[test]
    fn camera_pan_and_common_camera_object_translation_have_expected_motion() {
        with_cx(|cx| {
            let current_transform = transform_z(0.0, [0.0, 0.0, -5.0]);
            let sample = hit_sample(cx, &instance(current_transform));
            let raster = RasterSize::try_new(200, 100).unwrap();
            let previous = frame(-1.0, 9, camera_at(0.0, 0.0), current_transform);
            let current = frame(0.0, 9, camera_at(0.0, 0.0), current_transform);
            let camera_pan = frame(1.0, 9, camera_at(1.0, 0.0), current_transform);
            let panned = available(
                compute_motion_vectors(sample, &previous, &current, &camera_pan, raster).unwrap(),
            );
            let (_, pan_ndc, _) = projected(panned.next);
            assert_close(pan_ndc[0], -0.2);

            let common = frame(
                1.0,
                9,
                camera_at(1.0, 0.0),
                transform_z(0.0, [1.0, 0.0, -5.0]),
            );
            let common_motion = available(
                compute_motion_vectors(sample, &previous, &current, &common, raster).unwrap(),
            );
            let (_, common_ndc, common_pixels) = projected(common_motion.next);
            assert_eq!(common_ndc.map(f64::to_bits), [0.0f64.to_bits(); 2]);
            assert_eq!(common_pixels.map(f64::to_bits), [0.0f64.to_bits(); 2]);
        });
    }

    #[test]
    fn aperture_is_irrelevant_and_quaternion_double_cover_is_identical() {
        with_cx(|cx| {
            let positive = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -5.0]).unwrap();
            let negative =
                RigidTransform::try_new([0.0, 0.0, 0.0, -1.0], [0.0, 0.0, -5.0]).unwrap();
            let sample = hit_sample(cx, &instance(positive));
            let raster = RasterSize::try_new(3840, 2160).unwrap();
            let previous = frame(-1.0, 11, camera_at(0.0, 0.025), negative);
            let current = frame(0.0, 11, camera_at(0.0, 0.0), positive);
            let next = frame(1.0, 11, camera_at(0.0, 0.050), negative);
            assert_eq!(previous.frame_identity(), current.frame_identity());
            let motion = available(
                compute_motion_vectors(sample, &previous, &current, &next, raster).unwrap(),
            );
            for endpoint in [motion.previous, motion.next] {
                let (_, ndc, pixels) = projected(endpoint);
                assert_eq!(ndc.map(f64::to_bits), [0.0f64.to_bits(); 2]);
                assert_eq!(pixels.map(f64::to_bits), [0.0f64.to_bits(); 2]);
            }
        });
    }

    #[test]
    fn cuts_offscreen_and_behind_camera_are_explicit() {
        with_cx(|cx| {
            let current_transform = transform_z(0.0, [0.0, 0.0, -5.0]);
            let sample = hit_sample(cx, &instance(current_transform));
            let raster = RasterSize::try_new(200, 100).unwrap();
            let previous = frame(-1.0, 3, camera_at(0.0, 0.0), current_transform);
            let current = frame(0.0, 4, camera_at(0.0, 0.0), current_transform);
            let offscreen = frame(
                1.0,
                4,
                camera_at(0.0, 0.0),
                transform_z(0.0, [100.0, 0.0, -5.0]),
            );
            let motion = available(
                compute_motion_vectors(sample, &previous, &current, &offscreen, raster).unwrap(),
            );
            assert!(matches!(motion.previous, MotionEndpoint::CameraCut { .. }));
            let (target, _, _) = projected(motion.next);
            assert!(!target.in_frame());
            let tolerance = ReprojectionTolerance::try_new(1.0e-6, 1.0e-6, 1.0e-6).unwrap();
            assert_eq!(
                validate_reprojection(
                    sample,
                    motion.next,
                    ReprojectionObservation::Background,
                    tolerance,
                ),
                ReprojectionStatus::TargetOffScreen
            );

            let behind = frame(
                1.0,
                4,
                camera_at(0.0, 0.0),
                transform_z(0.0, [0.0, 0.0, 1.0]),
            );
            let behind_motion = available(
                compute_motion_vectors(sample, &current, &current, &behind, raster).unwrap(),
            );
            assert!(matches!(
                behind_motion.next,
                MotionEndpoint::BehindCamera { .. }
            ));
        });
    }

    #[test]
    fn chart_correspondence_refuses_without_fabricating_motion() {
        with_cx(|cx| {
            let instance = GeometryInstance::try_new(
                91,
                identity("sphere"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                }),
                RigidTransform::identity(),
            )
            .unwrap();
            let hit = instance
                .intersect(
                    cx,
                    &Ray {
                        origin: Point3::new(0.0, 0.0, 3.0),
                        dir: Vec3::new(0.0, 0.0, -1.0),
                    },
                    8.0,
                    1.0e-8,
                )
                .unwrap()
                .unwrap();
            let sample =
                PrimarySurfaceSample::try_from_instance_hit(&hit, identity("mat")).unwrap();
            assert_eq!(
                sample.identity().feature(),
                StableFeatureIdentity::ChartUnavailable
            );
            let frame = MotionFrame::from_instance(0.0, 1, camera_at(0.0, 0.0), &instance).unwrap();
            assert_eq!(
                compute_motion_vectors(
                    sample,
                    &frame,
                    &frame,
                    &frame,
                    RasterSize::try_new(64, 64).unwrap(),
                )
                .unwrap(),
                MotionVectorComputation::Unavailable {
                    identity: sample.identity(),
                    reason: CorrespondenceUnavailable::ChartHasNoStableParameter,
                }
            );
        });
    }

    #[test]
    fn reprojection_distinguishes_visibility_occlusion_and_topology() {
        with_cx(|cx| {
            let transform = transform_z(0.0, [0.0, 0.0, -5.0]);
            let sample = hit_sample(cx, &instance(transform));
            let frame = frame(0.0, 5, camera_at(0.0, 0.0), transform);
            let motion = available(
                compute_motion_vectors(
                    sample,
                    &frame,
                    &frame,
                    &frame,
                    RasterSize::try_new(200, 100).unwrap(),
                )
                .unwrap(),
            );
            let SurfaceCorrespondence::Mesh(surface) = sample.correspondence() else {
                unreachable!();
            };
            let tolerance = ReprojectionTolerance::try_new(1.0e-8, 1.0e-8, 1.0e-8).unwrap();
            let visible = ReprojectionObservation::try_surface(
                sample.identity(),
                motion.current.depth_m(),
                Some(surface.barycentric()),
            )
            .unwrap();
            assert_eq!(
                validate_reprojection(sample, motion.next, visible, tolerance),
                ReprojectionStatus::VisibleAtTarget
            );
            assert_eq!(
                validate_reprojection(
                    sample,
                    motion.next,
                    ReprojectionObservation::Background,
                    tolerance,
                ),
                ReprojectionStatus::DisoccludedAtTarget
            );
            let foreign_identity = StableHitIdentity {
                object_id: 999,
                ..sample.identity()
            };
            let occluder = ReprojectionObservation::try_surface(
                foreign_identity,
                motion.current.depth_m() - 1.0,
                Some(surface.barycentric()),
            )
            .unwrap();
            assert_eq!(
                validate_reprojection(sample, motion.next, occluder, tolerance),
                ReprojectionStatus::OccludedAtTarget
            );
            let mut wrong_barycentric = surface.barycentric();
            wrong_barycentric.swap(0, 1);
            let ambiguous = ReprojectionObservation::try_surface(
                sample.identity(),
                motion.current.depth_m(),
                Some(wrong_barycentric),
            )
            .unwrap();
            assert_eq!(
                validate_reprojection(sample, motion.next, ambiguous, tolerance),
                ReprojectionStatus::TopologyAmbiguous
            );
        });
    }

    #[test]
    fn animated_frames_bind_exact_shot_and_pose_and_replay_deterministically() {
        with_cx(|cx| {
            let start = transform_z(0.0, [0.0, 0.0, -5.0]);
            let end = transform_z(core::f64::consts::FRAC_PI_2, [1.0, 0.0, -5.0]);
            let trajectory = RigidTransformTrajectory::try_new(vec![
                TransformKeyframe::try_new(0.0, start, [1.0, 0.0, 0.0]).unwrap(),
                TransformKeyframe::try_new(1.0, end, [1.0, 0.0, 0.0]).unwrap(),
            ])
            .unwrap();
            let animated = AnimatedGeometryInstance::try_new(
                41,
                identity("mesh"),
                SharedGeometry::mesh(triangle()),
                trajectory,
            )
            .unwrap();
            let physical = camera_at(0.0, 0.0);
            let camera = AnimatedCamera::try_new(vec![
                CameraShot::try_new(
                    23,
                    0.0,
                    1.0,
                    vec![CameraKeyframe::try_new(0.0, physical.clone()).unwrap()],
                )
                .unwrap(),
            ])
            .unwrap();
            let current_instance = animated.instance_at(cx, 0.0).unwrap();
            let sample = hit_sample(cx, &current_instance);
            let previous =
                MotionFrame::from_animated(cx, 0.0, CutSide::After, &camera, &animated).unwrap();
            let current = previous.clone();
            let next =
                MotionFrame::from_animated(cx, 1.0, CutSide::Before, &camera, &animated).unwrap();
            assert_eq!(current.shot_id(), 23);
            assert_ne!(current.frame_identity(), next.frame_identity());
            let raster = RasterSize::try_new(3840, 2160).unwrap();
            let first = compute_motion_vectors(sample, &previous, &current, &next, raster).unwrap();
            let second =
                compute_motion_vectors(sample, &previous, &current, &next, raster).unwrap();
            assert_eq!(first, second);
            assert_eq!(raster.pixel_count(), 8_294_400);
            assert_eq!(raster.linear_index(3839, 2159), Some(8_294_399));
            assert_eq!(raster.linear_index(3840, 0), None);
        });
    }

    #[test]
    fn invalid_inputs_and_foreign_current_pose_fail_closed() {
        with_cx(|cx| {
            assert_eq!(
                RasterSize::try_new(0, 10),
                Err(MotionVectorError::InvalidRasterSize)
            );
            assert!(ReprojectionTolerance::try_new(-1.0, 0.0, 0.0).is_err());
            let transform = transform_z(0.0, [0.0, 0.0, -5.0]);
            let source = instance(transform);
            let hit = source
                .intersect(
                    cx,
                    &Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(1.0, 0.0, -5.0),
                    },
                    2.0,
                    1.0e-9,
                )
                .unwrap()
                .unwrap();
            assert_eq!(
                PrimarySurfaceSample::try_from_instance_hit(&hit, ContentHash([0; 32])),
                Err(MotionVectorError::InvalidMaterialIdentity)
            );
            let sample =
                PrimarySurfaceSample::try_from_instance_hit(&hit, identity("mat")).unwrap();
            let current = frame(
                0.0,
                2,
                camera_at(0.0, 0.0),
                transform_z(0.0, [0.5, 0.0, -5.0]),
            );
            assert_eq!(
                compute_motion_vectors(
                    sample,
                    &current,
                    &current,
                    &current,
                    RasterSize::try_new(100, 100).unwrap(),
                ),
                Err(MotionVectorError::ForeignCurrentPose)
            );
        });
    }
}

//! Time-varying proper-rigid instances evaluated at absolute shutter time.
//!
//! Translation uses cubic Hermite interpolation with producer-supplied endpoint
//! velocities. Rotation uses shortest-arc quaternion interpolation. This module
//! reconstructs only the admitted trajectory; it neither extrapolates nor adds
//! mechanical bandwidth to its keyframes.

use core::fmt;

use fs_blake3::ContentHash;
use fs_exec::{Cancelled, Cx};
use fs_math::det;

use crate::charts::Ray;
use crate::instances::{
    GeometryInstance, InstanceError, InstanceHit, RigidTransform, SharedGeometry,
};
use crate::motion::{ShutterInterval, TimedRay};

/// One admitted pose sample and its world-space translational velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransformKeyframe {
    absolute_time_s: f64,
    transform: RigidTransform,
    linear_velocity_m_per_s: [f64; 3],
}

impl TransformKeyframe {
    /// Admit a finite absolute time and finite world-space velocity.
    pub fn try_new(
        absolute_time_s: f64,
        transform: RigidTransform,
        linear_velocity_m_per_s: [f64; 3],
    ) -> Result<Self, AnimatedInstanceError> {
        if !absolute_time_s.is_finite() {
            return Err(AnimatedInstanceError::InvalidKeyframeTime);
        }
        if linear_velocity_m_per_s
            .iter()
            .any(|component| !component.is_finite())
        {
            return Err(AnimatedInstanceError::InvalidKeyframeVelocity);
        }
        Ok(Self {
            // Numeric ordering treats both signs of zero as one instant. Store
            // that instant canonically so the total-order binary search below
            // cannot distinguish a query solely by its zero sign.
            absolute_time_s: canonical_time(absolute_time_s),
            transform,
            linear_velocity_m_per_s,
        })
    }

    /// Absolute sample time [s].
    #[must_use]
    pub const fn absolute_time_s(self) -> f64 {
        self.absolute_time_s
    }

    /// Admitted body-to-world transform.
    #[must_use]
    pub const fn transform(self) -> RigidTransform {
        self.transform
    }

    /// World-space translational velocity [m/s].
    #[must_use]
    pub const fn linear_velocity_m_per_s(self) -> [f64; 3] {
        self.linear_velocity_m_per_s
    }
}

/// A reconstructed transform and translation derivative at one exact time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaluatedTransform {
    /// Requested absolute time [s].
    pub absolute_time_s: f64,
    /// Interpolated proper-rigid body-to-world transform.
    pub transform: RigidTransform,
    /// Derivative of the Hermite translation [m/s].
    pub linear_velocity_m_per_s: [f64; 3],
}

/// Strictly time-ordered proper-rigid transform samples.
#[derive(Clone, Debug)]
pub struct RigidTransformTrajectory {
    keyframes: Vec<TransformKeyframe>,
}

impl RigidTransformTrajectory {
    /// Admit one or more keyframes in strictly increasing absolute-time order.
    ///
    /// A single keyframe is a static trajectory valid at that exact time.
    pub fn try_new(keyframes: Vec<TransformKeyframe>) -> Result<Self, AnimatedInstanceError> {
        if keyframes.is_empty() {
            return Err(AnimatedInstanceError::EmptyTrajectory);
        }
        for pair in keyframes.windows(2) {
            let duration_s = pair[1].absolute_time_s - pair[0].absolute_time_s;
            if !duration_s.is_finite() || duration_s <= 0.0 {
                return Err(AnimatedInstanceError::NonIncreasingKeyframeTime);
            }
            for velocity in [
                pair[0].linear_velocity_m_per_s,
                pair[1].linear_velocity_m_per_s,
            ] {
                if velocity
                    .iter()
                    .any(|component| !(duration_s * component).is_finite())
                {
                    return Err(AnimatedInstanceError::InvalidInterpolation);
                }
            }
        }
        Ok(Self { keyframes })
    }

    /// Ordered trajectory keyframes.
    #[must_use]
    pub fn keyframes(&self) -> &[TransformKeyframe] {
        &self.keyframes
    }

    /// Earliest admitted absolute time [s].
    #[must_use]
    pub fn start_time_s(&self) -> f64 {
        self.keyframes[0].absolute_time_s
    }

    /// Latest admitted absolute time [s].
    #[must_use]
    pub fn end_time_s(&self) -> f64 {
        self.keyframes[self.keyframes.len() - 1].absolute_time_s
    }

    /// Require the complete resolved exposure to lie inside this trajectory.
    pub fn admit_shutter(&self, shutter: ShutterInterval) -> Result<(), AnimatedInstanceError> {
        if shutter.open_s() < self.start_time_s() || shutter.close_s() > self.end_time_s() {
            return Err(AnimatedInstanceError::ShutterOutsideTrajectory);
        }
        Ok(())
    }

    /// Evaluate without extrapolation at one finite absolute time.
    pub fn evaluate(
        &self,
        absolute_time_s: f64,
    ) -> Result<EvaluatedTransform, AnimatedInstanceError> {
        if !absolute_time_s.is_finite() {
            return Err(AnimatedInstanceError::InvalidEvaluationTime);
        }
        let absolute_time_s = canonical_time(absolute_time_s);
        if absolute_time_s < self.start_time_s() || absolute_time_s > self.end_time_s() {
            return Err(AnimatedInstanceError::Extrapolation);
        }
        match self
            .keyframes
            .binary_search_by(|keyframe| keyframe.absolute_time_s.total_cmp(&absolute_time_s))
        {
            Ok(index) => {
                let keyframe = self.keyframes[index];
                Ok(EvaluatedTransform {
                    absolute_time_s,
                    transform: keyframe.transform,
                    linear_velocity_m_per_s: keyframe.linear_velocity_m_per_s,
                })
            }
            Err(right_index) => {
                let left = self.keyframes[right_index - 1];
                let right = self.keyframes[right_index];
                interpolate_keyframes(left, right, absolute_time_s)
            }
        }
    }
}

fn canonical_time(time_s: f64) -> f64 {
    if time_s == 0.0 { 0.0 } else { time_s }
}

/// Immutable local geometry paired with a time-varying rigid placement.
#[derive(Clone, Debug)]
pub struct AnimatedGeometryInstance {
    prototype: GeometryInstance,
    trajectory: RigidTransformTrajectory,
}

impl AnimatedGeometryInstance {
    /// Bind stable object/local-geometry identities to one admitted trajectory.
    pub fn try_new(
        object_id: u64,
        geometry_identity: ContentHash,
        geometry: SharedGeometry,
        trajectory: RigidTransformTrajectory,
    ) -> Result<Self, AnimatedInstanceError> {
        let prototype = GeometryInstance::try_new(
            object_id,
            geometry_identity,
            geometry,
            trajectory.keyframes[0].transform,
        )?;
        Ok(Self {
            prototype,
            trajectory,
        })
    }

    /// Stable nonzero object ID retained at every ray time.
    #[must_use]
    pub const fn object_id(&self) -> u64 {
        self.prototype.object_id()
    }

    /// Immutable local-geometry identity retained at every ray time.
    #[must_use]
    pub const fn geometry_identity(&self) -> ContentHash {
        self.prototype.geometry_identity()
    }

    /// Shared immutable local geometry.
    #[must_use]
    pub const fn geometry(&self) -> &SharedGeometry {
        self.prototype.geometry()
    }

    /// Admitted transform trajectory.
    #[must_use]
    pub const fn trajectory(&self) -> &RigidTransformTrajectory {
        &self.trajectory
    }

    /// Materialize the static instance used for one absolute-time query.
    pub fn instance_at(
        &self,
        cx: &Cx<'_>,
        absolute_time_s: f64,
    ) -> Result<GeometryInstance, AnimatedInstanceError> {
        cx.checkpoint()?;
        let evaluated = self.trajectory.evaluate(absolute_time_s)?;
        let instance = GeometryInstance::try_new(
            self.object_id(),
            self.geometry_identity(),
            self.geometry().clone(),
            evaluated.transform,
        )?;
        cx.checkpoint()?;
        Ok(instance)
    }

    /// Intersect at `timed_ray.absolute_time_s()` with no pose extrapolation.
    ///
    /// The evaluated static instance retains the same object ID and immutable
    /// local-geometry identity. Its frame identity binds the interpolated pose.
    pub fn intersect(
        &self,
        cx: &Cx<'_>,
        timed_ray: &TimedRay<Ray>,
        t_max: f64,
        eps: f64,
    ) -> Result<Option<InstanceHit>, AnimatedInstanceError> {
        cx.checkpoint()?;
        let instance = self.instance_at(cx, timed_ray.absolute_time_s())?;
        let hit = instance.intersect(cx, timed_ray.spatial(), t_max, eps)?;
        cx.checkpoint()?;
        Ok(hit)
    }
}

/// Fail-closed trajectory admission, evaluation, and intersection errors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnimatedInstanceError {
    /// No pose samples were supplied.
    EmptyTrajectory,
    /// A keyframe time was NaN or infinite.
    InvalidKeyframeTime,
    /// A keyframe translational velocity was NaN or infinite.
    InvalidKeyframeVelocity,
    /// Keyframe times were duplicated, decreasing, or had no finite interval.
    NonIncreasingKeyframeTime,
    /// The resolved shutter is not wholly covered by trajectory samples.
    ShutterOutsideTrajectory,
    /// The requested evaluation time was NaN or infinite.
    InvalidEvaluationTime,
    /// The requested time lies outside the admitted keyframe interval.
    Extrapolation,
    /// Interpolation produced an unrepresentable transform or derivative.
    InvalidInterpolation,
    /// Static instance admission or intersection refused.
    Instance(InstanceError),
    /// Execution was cancelled at a bounded evaluation/intersection boundary.
    Cancelled,
}

impl From<Cancelled> for AnimatedInstanceError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<InstanceError> for AnimatedInstanceError {
    fn from(error: InstanceError) -> Self {
        match error {
            InstanceError::Cancelled => Self::Cancelled,
            other => Self::Instance(other),
        }
    }
}

impl fmt::Display for AnimatedInstanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "invalid animated render instance: {self:?}")
    }
}

impl std::error::Error for AnimatedInstanceError {}

fn interpolate_keyframes(
    left: TransformKeyframe,
    right: TransformKeyframe,
    absolute_time_s: f64,
) -> Result<EvaluatedTransform, AnimatedInstanceError> {
    let duration_s = right.absolute_time_s - left.absolute_time_s;
    let alpha = (absolute_time_s - left.absolute_time_s) / duration_s;
    let left_translation = left.transform.translation_m();
    let right_translation = right.transform.translation_m();
    let translation_m = core::array::from_fn(|axis| {
        hermite_scalar(
            left_translation[axis],
            left.linear_velocity_m_per_s[axis],
            right_translation[axis],
            right.linear_velocity_m_per_s[axis],
            duration_s,
            alpha,
        )
    });
    let linear_velocity_m_per_s = core::array::from_fn(|axis| {
        hermite_scalar_derivative(
            left_translation[axis],
            left.linear_velocity_m_per_s[axis],
            right_translation[axis],
            right.linear_velocity_m_per_s[axis],
            duration_s,
            alpha,
        )
    });
    if translation_m
        .iter()
        .chain(linear_velocity_m_per_s.iter())
        .any(|component| !component.is_finite())
    {
        return Err(AnimatedInstanceError::InvalidInterpolation);
    }
    let rotation_xyzw = slerp_shortest(
        left.transform.rotation_xyzw(),
        right.transform.rotation_xyzw(),
        alpha,
    )?;
    let transform = RigidTransform::try_new(rotation_xyzw, translation_m)
        .map_err(|_| AnimatedInstanceError::InvalidInterpolation)?;
    Ok(EvaluatedTransform {
        absolute_time_s,
        transform,
        linear_velocity_m_per_s,
    })
}

fn slerp_shortest(
    left: [f64; 4],
    mut right: [f64; 4],
    alpha: f64,
) -> Result<[f64; 4], AnimatedInstanceError> {
    let mut dot = left
        .iter()
        .zip(right)
        .map(|(first, second)| first * second)
        .sum::<f64>();
    if dot < 0.0 {
        for component in &mut right {
            *component = -*component;
        }
        dot = -dot;
    }
    dot = dot.clamp(-1.0, 1.0);
    let mut interpolated: [f64; 4] = if dot > 1.0 - 1.0e-12 {
        core::array::from_fn(|index| left[index] + alpha * (right[index] - left[index]))
    } else {
        let angle = det::acos(dot);
        let denominator = det::sin(angle);
        let left_weight = det::sin((1.0 - alpha) * angle) / denominator;
        let right_weight = det::sin(alpha * angle) / denominator;
        core::array::from_fn(|index| left_weight.mul_add(left[index], right_weight * right[index]))
    };
    let norm = interpolated
        .iter()
        .map(|component| component * component)
        .sum::<f64>()
        .sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return Err(AnimatedInstanceError::InvalidInterpolation);
    }
    for component in &mut interpolated {
        *component /= norm;
    }
    Ok(interpolated)
}

fn hermite_scalar(
    left_position: f64,
    left_velocity: f64,
    right_position: f64,
    right_velocity: f64,
    duration_s: f64,
    alpha: f64,
) -> f64 {
    let alpha_squared = alpha * alpha;
    let alpha_cubed = alpha_squared * alpha;
    let left_position_weight = 2.0 * alpha_cubed - 3.0 * alpha_squared + 1.0;
    let left_velocity_weight = alpha_cubed - 2.0 * alpha_squared + alpha;
    let right_position_weight = -2.0 * alpha_cubed + 3.0 * alpha_squared;
    let right_velocity_weight = alpha_cubed - alpha_squared;
    left_position_weight * left_position
        + left_velocity_weight * duration_s * left_velocity
        + right_position_weight * right_position
        + right_velocity_weight * duration_s * right_velocity
}

fn hermite_scalar_derivative(
    left_position: f64,
    left_velocity: f64,
    right_position: f64,
    right_velocity: f64,
    duration_s: f64,
    alpha: f64,
) -> f64 {
    let alpha_squared = alpha * alpha;
    ((6.0 * alpha_squared - 6.0 * alpha) * left_position
        + (-6.0 * alpha_squared + 6.0 * alpha) * right_position)
        / duration_s
        + (3.0 * alpha_squared - 4.0 * alpha + 1.0) * left_velocity
        + (3.0 * alpha_squared - 2.0 * alpha) * right_velocity
}

#[cfg(test)]
mod tests {
    use asupersync::types::Budget;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};

    use super::*;
    use crate::charts::TriMesh;
    use crate::instances::InstanceBackendAudit;
    use crate::motion::{
        NormalizedShutterTime, ShotTimeBounds, ShutterConvention, ShutterDistribution,
    };

    fn keyframe(
        time_s: f64,
        rotation_xyzw: [f64; 4],
        translation_m: [f64; 3],
        velocity_m_per_s: [f64; 3],
    ) -> TransformKeyframe {
        TransformKeyframe::try_new(
            time_s,
            RigidTransform::try_new(rotation_xyzw, translation_m).unwrap(),
            velocity_m_per_s,
        )
        .unwrap()
    }

    fn trajectory() -> RigidTransformTrajectory {
        RigidTransformTrajectory::try_new(vec![
            keyframe(2.0, [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            keyframe(4.0, [0.0, 0.0, 0.0, 1.0], [4.0, 0.0, 0.0], [3.0, 0.0, 0.0]),
        ])
        .unwrap()
    }

    fn shutter(open_s: f64, duration_s: f64) -> ShutterInterval {
        ShutterInterval::resolve(
            open_s,
            duration_s,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(-10.0, 10.0).unwrap(),
        )
        .unwrap()
    }

    fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 53,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&gate, &cx)
        })
    }

    #[test]
    fn admission_and_evaluation_refuse_malformed_or_uncovered_times() {
        assert_eq!(
            TransformKeyframe::try_new(f64::NAN, RigidTransform::identity(), [0.0; 3]),
            Err(AnimatedInstanceError::InvalidKeyframeTime)
        );
        assert_eq!(
            TransformKeyframe::try_new(0.0, RigidTransform::identity(), [f64::INFINITY, 0.0, 0.0]),
            Err(AnimatedInstanceError::InvalidKeyframeVelocity)
        );
        assert_eq!(
            RigidTransformTrajectory::try_new(Vec::new()).unwrap_err(),
            AnimatedInstanceError::EmptyTrajectory
        );
        let duplicate = vec![
            keyframe(1.0, [0.0, 0.0, 0.0, 1.0], [0.0; 3], [0.0; 3]),
            keyframe(1.0, [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0; 3]),
        ];
        assert_eq!(
            RigidTransformTrajectory::try_new(duplicate).unwrap_err(),
            AnimatedInstanceError::NonIncreasingKeyframeTime
        );

        let trajectory = trajectory();
        assert_eq!(
            trajectory.evaluate(f64::NAN),
            Err(AnimatedInstanceError::InvalidEvaluationTime)
        );
        assert_eq!(
            trajectory.evaluate(1.999),
            Err(AnimatedInstanceError::Extrapolation)
        );
        assert_eq!(
            trajectory.evaluate(4.001),
            Err(AnimatedInstanceError::Extrapolation)
        );
        assert_eq!(trajectory.admit_shutter(shutter(2.0, 2.0)), Ok(()));
        assert_eq!(
            trajectory.admit_shutter(shutter(1.5, 1.0)),
            Err(AnimatedInstanceError::ShutterOutsideTrajectory)
        );
    }

    #[test]
    fn signed_zero_times_are_one_canonical_instant_without_index_underflow() {
        let singleton = RigidTransformTrajectory::try_new(vec![keyframe(
            0.0,
            [0.0, 0.0, 0.0, 1.0],
            [3.0, 0.0, 0.0],
            [0.0; 3],
        )])
        .unwrap();
        let singleton_at_negative_zero = singleton.evaluate(-0.0).unwrap();
        assert_eq!(
            singleton_at_negative_zero.absolute_time_s.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            singleton_at_negative_zero.transform.translation_m()[0].to_bits(),
            3.0_f64.to_bits()
        );

        let trajectory = RigidTransformTrajectory::try_new(vec![
            keyframe(-0.0, [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0; 3]),
            keyframe(1.0, [0.0, 0.0, 0.0, 1.0], [2.0, 0.0, 0.0], [0.0; 3]),
        ])
        .unwrap();
        assert_eq!(trajectory.start_time_s().to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            trajectory.evaluate(0.0).unwrap().transform.translation_m()[0].to_bits(),
            1.0_f64.to_bits()
        );
    }

    #[test]
    fn hermite_translation_retains_endpoint_velocities_and_smooth_midpoint() {
        let trajectory = trajectory();
        let left = trajectory.evaluate(2.0).unwrap();
        let middle = trajectory.evaluate(3.0).unwrap();
        let right = trajectory.evaluate(4.0).unwrap();
        assert_eq!(
            left.transform.translation_m().map(f64::to_bits),
            [0.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            left.linear_velocity_m_per_s.map(f64::to_bits),
            [1.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            middle.transform.translation_m().map(f64::to_bits),
            [1.5, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            middle.linear_velocity_m_per_s.map(f64::to_bits),
            [2.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            right.transform.translation_m().map(f64::to_bits),
            [4.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            right.linear_velocity_m_per_s.map(f64::to_bits),
            [3.0, 0.0, 0.0].map(f64::to_bits)
        );
    }

    #[test]
    fn quaternion_interpolation_uses_the_shortest_rotation_arc() {
        let half_270 = 3.0 * core::f64::consts::FRAC_PI_4;
        let trajectory = RigidTransformTrajectory::try_new(vec![
            keyframe(0.0, [0.0, 0.0, 0.0, 1.0], [0.0; 3], [0.0; 3]),
            keyframe(
                2.0,
                [0.0, 0.0, half_270.sin(), half_270.cos()],
                [0.0; 3],
                [0.0; 3],
            ),
        ])
        .unwrap();
        let rotated = trajectory
            .evaluate(1.0)
            .unwrap()
            .transform
            .transform_vector(Vec3::new(1.0, 0.0, 0.0));
        let expected = core::f64::consts::FRAC_1_SQRT_2;
        assert!((rotated.x - expected).abs() < 2.0e-12);
        assert!((rotated.y + expected).abs() < 2.0e-12);
        assert!(rotated.z.abs() < 2.0e-12);
    }

    #[test]
    fn timed_intersection_uses_absolute_time_and_preserves_local_identity() {
        with_cx(|gate, cx| {
            let mesh = TriMesh::new(
                vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[0, 1, 2]],
            );
            let bvh_fingerprint = mesh.bvh_fingerprint();
            let geometry = SharedGeometry::mesh(mesh);
            let identity = hash_domain("org.frankensim.test.animated-instance", b"triangle");
            let trajectory = RigidTransformTrajectory::try_new(vec![
                keyframe(0.0, [0.0, 0.0, 0.0, 1.0], [0.0; 3], [1.0, 0.0, 0.0]),
                keyframe(2.0, [0.0, 0.0, 0.0, 1.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            ])
            .unwrap();
            let animated =
                AnimatedGeometryInstance::try_new(17, identity, geometry, trajectory).unwrap();
            let exposure = shutter(0.0, 2.0);
            animated.trajectory().admit_shutter(exposure).unwrap();
            let ray = TimedRay::at_normalized(
                Ray {
                    origin: Point3::new(1.0, 0.0, 2.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                },
                exposure,
                NormalizedShutterTime::try_new(0.5).unwrap(),
            );
            let hit = animated.intersect(cx, &ray, 4.0, 1.0e-9).unwrap().unwrap();
            assert_eq!(hit.object_id, 17);
            assert_eq!(hit.geometry_identity, identity);
            assert_eq!(
                hit.backend_audit,
                InstanceBackendAudit::Mesh { bvh_fingerprint }
            );
            assert!((hit.hit.t - 2.0).abs() < 2.0e-12);
            assert!((hit.hit.point.x - 1.0).abs() < 2.0e-12);
            let timed_instance = animated.instance_at(cx, 1.0).unwrap();
            assert!(timed_instance.geometry().ptr_eq(animated.geometry()));
            assert_eq!(hit.frame_identity, timed_instance.frame_identity());

            let outside = TimedRay::at_normalized(
                *ray.spatial(),
                shutter(3.0, 1.0),
                NormalizedShutterTime::try_new(0.0).unwrap(),
            );
            assert_eq!(
                animated.intersect(cx, &outside, 4.0, 1.0e-9),
                Err(AnimatedInstanceError::Extrapolation)
            );

            gate.request();
            assert_eq!(
                animated.intersect(cx, &ray, 4.0, 1.0e-9),
                Err(AnimatedInstanceError::Cancelled)
            );
        });
    }

    #[test]
    fn equal_camera_and_object_translation_preserves_relative_intersection() {
        with_cx(|_, cx| {
            let mesh = TriMesh::new(
                vec![[-0.25, -0.25, 0.0], [0.25, -0.25, 0.0], [0.0, 0.25, 0.0]],
                vec![[0, 1, 2]],
            );
            let trajectory = RigidTransformTrajectory::try_new(vec![
                keyframe(0.0, [0.0, 0.0, 0.0, 1.0], [0.0; 3], [1.0, 0.0, 0.0]),
                keyframe(1.0, [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            ])
            .unwrap();
            let animated = AnimatedGeometryInstance::try_new(
                29,
                hash_domain(
                    "org.frankensim.test.animated-relative-motion",
                    b"small-triangle",
                ),
                SharedGeometry::mesh(mesh),
                trajectory,
            )
            .unwrap();
            let exposure = shutter(0.0, 1.0);
            let start_ray = TimedRay::at_normalized(
                Ray {
                    origin: Point3::new(0.0, 0.0, 2.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                },
                exposure,
                NormalizedShutterTime::try_new(0.0).unwrap(),
            );
            let co_moving_end_ray = TimedRay::at_normalized(
                Ray {
                    origin: Point3::new(1.0, 0.0, 2.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                },
                exposure,
                NormalizedShutterTime::try_new(1.0).unwrap(),
            );
            let fixed_camera_end_ray = TimedRay::at_normalized(
                *start_ray.spatial(),
                exposure,
                NormalizedShutterTime::try_new(1.0).unwrap(),
            );

            let start_hit = animated
                .intersect(cx, &start_ray, 4.0, 1.0e-9)
                .unwrap()
                .unwrap();
            let co_moving_hit = animated
                .intersect(cx, &co_moving_end_ray, 4.0, 1.0e-9)
                .unwrap()
                .unwrap();
            assert_eq!(start_hit.hit.t.to_bits(), co_moving_hit.hit.t.to_bits());
            assert_eq!(start_hit.hit.point.x.to_bits(), 0.0_f64.to_bits());
            assert_eq!(co_moving_hit.hit.point.x.to_bits(), 1.0_f64.to_bits());
            assert_eq!(
                animated
                    .intersect(cx, &fixed_camera_end_ray, 4.0, 1.0e-9)
                    .unwrap(),
                None,
                "object motion must still matter when the camera stays fixed",
            );
        });
    }
}

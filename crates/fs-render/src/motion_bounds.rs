//! Conservative world-space bounds for one rigid-motion segment.
//!
//! This leaf deliberately does not sample quaternion interpolation.  For a
//! segment whose orientation may vary, every local point is enclosed by an
//! outward-rounded sphere about the body origin. Linear motion may use a tight
//! endpoint translation envelope; runtime cubic-Hermite trajectories use the
//! convex hull of their equivalent Bezier translation controls. The result
//! therefore contains translational overshoot and arbitrary proper rotation,
//! including additional full spins that endpoint quaternions cannot encode.
//! A separately declared constant identity rotation uses the exact translated
//! local-box envelope.

use core::fmt;

use fs_geom::{Aabb, Point3};

use crate::animated_instances::RigidTransformTrajectory;
use crate::instances::RigidTransform;
use crate::motion::ShutterInterval;

/// Typed refusal from conservative rigid-motion bounding.
#[derive(Debug, Clone, PartialEq)]
pub enum MotionBoundsError {
    /// A local bound coordinate is not finite.
    NonFiniteLocalAabb {
        /// Coordinate that failed validation.
        field: &'static str,
    },
    /// A local bound axis has `min > max`.
    InvertedLocalAabb {
        /// Axis that failed validation.
        axis: &'static str,
    },
    /// Segment times are non-finite or do not define a positive interval.
    InvalidTimeInterval,
    /// Constant rotation was declared, but the canonical endpoint rotations differ.
    ConstantRotationMismatch,
    /// The requested shutter is not wholly covered by the trajectory.
    ShutterOutsideTrajectory,
    /// Finite inputs require a non-finite derived bound in binary64.
    Unrepresentable {
        /// Derived quantity that could not be represented finitely.
        field: &'static str,
    },
}

impl fmt::Display for MotionBoundsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for MotionBoundsError {}

/// A finite, ordered local-space box with a retained conservative radius about
/// the body origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteLocalAabb {
    bounds: Aabb,
    origin_radius_upper_m: f64,
}

impl FiniteLocalAabb {
    /// Validate a local box without normalizing malformed public fields.
    ///
    /// Degenerate axes and point boxes are admitted.  Unbounded boxes and
    /// boxes whose body-origin radius overflows are refused because this leaf
    /// promises a finite world AABB.
    pub fn try_new(bounds: Aabb) -> Result<Self, MotionBoundsError> {
        for (field, value) in [
            ("min.x", bounds.min.x),
            ("min.y", bounds.min.y),
            ("min.z", bounds.min.z),
            ("max.x", bounds.max.x),
            ("max.y", bounds.max.y),
            ("max.z", bounds.max.z),
        ] {
            if !value.is_finite() {
                return Err(MotionBoundsError::NonFiniteLocalAabb { field });
            }
        }
        for (axis, minimum, maximum) in [
            ("x", bounds.min.x, bounds.max.x),
            ("y", bounds.min.y, bounds.max.y),
            ("z", bounds.min.z, bounds.max.z),
        ] {
            if minimum > maximum {
                return Err(MotionBoundsError::InvertedLocalAabb { axis });
            }
        }
        let axis_radius = [
            bounds.min.x.abs().max(bounds.max.x.abs()),
            bounds.min.y.abs().max(bounds.max.y.abs()),
            bounds.min.z.abs().max(bounds.max.z.abs()),
        ];
        let origin_radius_upper_m =
            upper_euclidean_norm(axis_radius).ok_or(MotionBoundsError::Unrepresentable {
                field: "local_origin_radius_m",
            })?;
        Ok(Self {
            bounds,
            origin_radius_upper_m,
        })
    }

    /// Retained finite local-space box.
    #[must_use]
    pub const fn bounds(self) -> Aabb {
        self.bounds
    }

    /// Outward-rounded radius of a body-origin sphere containing the box [m].
    #[must_use]
    pub const fn origin_radius_upper_m(self) -> f64 {
        self.origin_radius_upper_m
    }
}

/// Rotation knowledge supplied by the motion owner.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotationSweep {
    /// Endpoint rotations must be identical and orientation remains constant.
    Constant,
    /// Orientation may follow any proper-rotation path, including extra spins.
    Arbitrary,
}

/// Two validated rigid poses over a positive time interval.
///
/// Translation is defined to interpolate linearly between endpoint body-origin
/// translations.  [`RotationSweep::Arbitrary`] makes no interpolation or
/// shortest-arc assumption about orientation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidMotionSegment {
    start: RigidTransform,
    end: RigidTransform,
    start_time_s: f64,
    end_time_s: f64,
    rotation_sweep: RotationSweep,
}

impl RigidMotionSegment {
    /// Admit a finite positive interval and a rotation declaration consistent
    /// with its canonical endpoint transforms.
    pub fn try_new(
        start: RigidTransform,
        end: RigidTransform,
        start_time_s: f64,
        end_time_s: f64,
        rotation_sweep: RotationSweep,
    ) -> Result<Self, MotionBoundsError> {
        if !start_time_s.is_finite() || !end_time_s.is_finite() || start_time_s >= end_time_s {
            return Err(MotionBoundsError::InvalidTimeInterval);
        }
        if rotation_sweep == RotationSweep::Constant
            && !same_rotation_bits(start.rotation_xyzw(), end.rotation_xyzw())
        {
            return Err(MotionBoundsError::ConstantRotationMismatch);
        }
        Ok(Self {
            start,
            end,
            start_time_s,
            end_time_s,
            rotation_sweep,
        })
    }

    /// Initial body-to-world transform.
    #[must_use]
    pub const fn start(self) -> RigidTransform {
        self.start
    }

    /// Final body-to-world transform.
    #[must_use]
    pub const fn end(self) -> RigidTransform {
        self.end
    }

    /// Closed motion time interval `[start, end]` [s].
    #[must_use]
    pub const fn time_interval_s(self) -> [f64; 2] {
        [self.start_time_s, self.end_time_s]
    }

    /// Retained orientation-path declaration.
    #[must_use]
    pub const fn rotation_sweep(self) -> RotationSweep {
        self.rotation_sweep
    }
}

/// Compute a finite conservative world AABB for the complete motion segment.
///
/// Arbitrary rotation is bounded without sampling: the linearly translated
/// body origin is inflated by a sphere containing every local-box point.
/// Constant identity rotation uses the tighter analytic translation envelope.
/// The rotation sphere includes an operation-count margin for the renderer's
/// finite-precision quaternion normalization and vector rotation; final world
/// sums are additionally rounded outward by one binary64 step.
pub fn conservative_world_swept_aabb(
    local: FiniteLocalAabb,
    motion: RigidMotionSegment,
) -> Result<Aabb, MotionBoundsError> {
    let start_translation = motion.start.translation_m();
    let end_translation = motion.end.translation_m();
    let translation_min = [
        start_translation[0].min(end_translation[0]),
        start_translation[1].min(end_translation[1]),
        start_translation[2].min(end_translation[2]),
    ];
    let translation_max = [
        start_translation[0].max(end_translation[0]),
        start_translation[1].max(end_translation[1]),
        start_translation[2].max(end_translation[2]),
    ];

    let (lower_terms, upper_terms) = if motion.rotation_sweep == RotationSweep::Constant
        && same_rotation_bits(
            motion.start.rotation_xyzw(),
            RigidTransform::identity().rotation_xyzw(),
        ) {
        (
            [local.bounds.min.x, local.bounds.min.y, local.bounds.min.z],
            [local.bounds.max.x, local.bounds.max.y, local.bounds.max.z],
        )
    } else {
        let radius = rotation_evaluation_radius(local.origin_radius_upper_m)?;
        ([-radius; 3], [radius; 3])
    };

    world_aabb_from_translation_envelope(translation_min, translation_max, lower_terms, upper_terms)
}

/// Bound the runtime cubic-Hermite translation and arbitrary proper rotation
/// over one admitted shutter.
///
/// Every overlapping keyframe segment is converted to its equivalent cubic
/// Bezier translation controls. A Bezier curve lies in the convex hull of its
/// controls, so this envelope includes velocity-driven interior overshoot that
/// endpoint-only motion bounds miss. Boundary segments are deliberately kept
/// whole when the shutter clips them: the result may be loose, but never relies
/// on a numerically fragile root solve or extrapolation. Rotation is enclosed
/// by the local body-origin sphere.
pub fn conservative_trajectory_swept_aabb(
    local: FiniteLocalAabb,
    trajectory: &RigidTransformTrajectory,
    shutter: ShutterInterval,
) -> Result<Aabb, MotionBoundsError> {
    if shutter.open_s() < trajectory.start_time_s() || shutter.close_s() > trajectory.end_time_s() {
        return Err(MotionBoundsError::ShutterOutsideTrajectory);
    }

    let keyframes = trajectory.keyframes();
    let mut translation_min = [f64::INFINITY; 3];
    let mut translation_max = [f64::NEG_INFINITY; 3];
    let mut scale = [0.0_f64; 3];
    if keyframes.len() == 1 {
        include_translation_control(
            keyframes[0].transform().translation_m(),
            &mut translation_min,
            &mut translation_max,
            &mut scale,
        )?;
    } else {
        for pair in keyframes.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if right.absolute_time_s() < shutter.open_s()
                || left.absolute_time_s() > shutter.close_s()
            {
                continue;
            }
            let duration_s = right.absolute_time_s() - left.absolute_time_s();
            let left_translation = left.transform().translation_m();
            let right_translation = right.transform().translation_m();
            let left_velocity = left.linear_velocity_m_per_s();
            let right_velocity = right.linear_velocity_m_per_s();
            include_translation_control(
                left_translation,
                &mut translation_min,
                &mut translation_max,
                &mut scale,
            )?;
            include_translation_control(
                right_translation,
                &mut translation_min,
                &mut translation_max,
                &mut scale,
            )?;
            // Multiply duration before dividing by three. Dividing a subnormal
            // velocity first can underflow to zero even though the admitted
            // duration-times-velocity product, and therefore the runtime
            // Hermite displacement, is representable.
            let left_tangent = core::array::from_fn(|axis| {
                left_translation[axis] + (duration_s * left_velocity[axis]) / 3.0
            });
            let right_tangent = core::array::from_fn(|axis| {
                right_translation[axis] - (duration_s * right_velocity[axis]) / 3.0
            });
            include_translation_control(
                left_tangent,
                &mut translation_min,
                &mut translation_max,
                &mut scale,
            )?;
            include_translation_control(
                right_tangent,
                &mut translation_min,
                &mut translation_max,
                &mut scale,
            )?;
        }
    }
    if translation_min[0].is_infinite() && translation_min[0].is_sign_positive() {
        return Err(MotionBoundsError::ShutterOutsideTrajectory);
    }

    // The runtime Hermite basis and the equivalent Bezier controls use
    // different floating-point evaluation orders. Expand the exact-real convex
    // hull by a conservative operation-count margin before adding geometry.
    for axis in 0..3 {
        let margin = (64.0 * f64::EPSILON * scale[axis]).max(64.0 * f64::from_bits(1));
        let lower = translation_min[axis] - margin;
        let upper = translation_max[axis] + margin;
        if !lower.is_finite() || !upper.is_finite() {
            return Err(MotionBoundsError::Unrepresentable {
                field: "hermite_translation_envelope",
            });
        }
        translation_min[axis] = next_down(lower);
        translation_max[axis] = next_up(upper);
    }

    let radius = rotation_evaluation_radius(local.origin_radius_upper_m)?;
    world_aabb_from_translation_envelope(
        translation_min,
        translation_max,
        [-radius; 3],
        [radius; 3],
    )
}

fn rotation_evaluation_radius(origin_radius_upper_m: f64) -> Result<f64, MotionBoundsError> {
    if origin_radius_upper_m.to_bits() == 0 {
        return Ok(0.0);
    }
    // RigidTransform normalizes an admitted quaternion and evaluates q*v*q^-1
    // with ordinary binary64 arithmetic. The real rotation preserves norm, but
    // that evaluation may enlarge one coordinate by several ulps. This
    // operation-count allowance is deliberately much larger than the concrete
    // multiply/add path and retains a subnormal floor.
    let margin = (64.0 * f64::EPSILON * origin_radius_upper_m).max(64.0 * f64::from_bits(1));
    let inflated = origin_radius_upper_m + margin;
    if !inflated.is_finite() {
        return Err(MotionBoundsError::Unrepresentable {
            field: "rotation_evaluation_radius_m",
        });
    }
    Ok(next_up(inflated))
}

fn include_translation_control(
    control: [f64; 3],
    minimum: &mut [f64; 3],
    maximum: &mut [f64; 3],
    scale: &mut [f64; 3],
) -> Result<(), MotionBoundsError> {
    for axis in 0..3 {
        if !control[axis].is_finite() {
            return Err(MotionBoundsError::Unrepresentable {
                field: "hermite_bezier_control",
            });
        }
        minimum[axis] = minimum[axis].min(control[axis]);
        maximum[axis] = maximum[axis].max(control[axis]);
        scale[axis] = scale[axis].max(control[axis].abs());
    }
    Ok(())
}

fn world_aabb_from_translation_envelope(
    translation_min: [f64; 3],
    translation_max: [f64; 3],
    lower_terms: [f64; 3],
    upper_terms: [f64; 3],
) -> Result<Aabb, MotionBoundsError> {
    let mut lower = [0.0; 3];
    let mut upper = [0.0; 3];
    for axis in 0..3 {
        lower[axis] = outward_sum_lower(translation_min[axis], lower_terms[axis]).ok_or(
            MotionBoundsError::Unrepresentable {
                field: "world_swept_aabb.min",
            },
        )?;
        upper[axis] = outward_sum_upper(translation_max[axis], upper_terms[axis]).ok_or(
            MotionBoundsError::Unrepresentable {
                field: "world_swept_aabb.max",
            },
        )?;
    }
    Ok(Aabb {
        min: Point3::new(lower[0], lower[1], lower[2]),
        max: Point3::new(upper[0], upper[1], upper[2]),
    })
}

fn upper_euclidean_norm(values: [f64; 3]) -> Option<f64> {
    let scale = values[0].max(values[1]).max(values[2]);
    if scale.to_bits() == 0 {
        return Some(0.0);
    }
    let mut squared_sum = 0.0;
    for value in values {
        let ratio = next_up(value / scale);
        let squared = next_up(ratio * ratio);
        squared_sum = next_up(squared_sum + squared);
    }
    let norm = next_up(scale * next_up(squared_sum.sqrt()));
    norm.is_finite().then_some(norm)
}

fn outward_sum_lower(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    if !sum.is_finite() {
        return None;
    }
    Some(if right.to_bits() << 1 == 0 {
        sum
    } else {
        next_down(sum)
    })
}

fn outward_sum_upper(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    if !sum.is_finite() {
        return None;
    }
    Some(if right.to_bits() << 1 == 0 {
        sum
    } else {
        next_up(sum)
    })
}

fn next_up(value: f64) -> f64 {
    if value.is_infinite() && value.is_sign_positive() {
        return value;
    }
    if value.to_bits() << 1 == 0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(value: f64) -> f64 {
    if value.is_infinite() && value.is_sign_negative() {
        return value;
    }
    if value.to_bits() << 1 == 0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    f64::from_bits(if value > 0.0 { bits - 1 } else { bits + 1 })
}

fn same_rotation_bits(left: [f64; 4], right: [f64; 4]) -> bool {
    left.into_iter()
        .zip(right)
        .all(|(left, right)| left.to_bits() == right.to_bits())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_geom::Vec3;

    use crate::animated_instances::TransformKeyframe;
    use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};

    fn transform(rotation_xyzw: [f64; 4], translation_m: [f64; 3]) -> RigidTransform {
        RigidTransform::try_new(rotation_xyzw, translation_m).expect("valid test transform")
    }

    fn rotation_z(angle: f64, translation_m: [f64; 3]) -> RigidTransform {
        let half = 0.5 * angle;
        transform([0.0, 0.0, half.sin(), half.cos()], translation_m)
    }

    fn corners(bounds: Aabb) -> [Point3; 8] {
        [
            Point3::new(bounds.min.x, bounds.min.y, bounds.min.z),
            Point3::new(bounds.min.x, bounds.min.y, bounds.max.z),
            Point3::new(bounds.min.x, bounds.max.y, bounds.min.z),
            Point3::new(bounds.min.x, bounds.max.y, bounds.max.z),
            Point3::new(bounds.max.x, bounds.min.y, bounds.min.z),
            Point3::new(bounds.max.x, bounds.min.y, bounds.max.z),
            Point3::new(bounds.max.x, bounds.max.y, bounds.min.z),
            Point3::new(bounds.max.x, bounds.max.y, bounds.max.z),
        ]
    }

    fn assert_contains_transform(bounds: Aabb, local: Aabb, transform: RigidTransform) {
        for corner in corners(local) {
            let world = transform.transform_point(corner);
            assert!(
                bounds.contains(world),
                "world corner {world:?} escaped {bounds:?}"
            );
        }
    }

    fn shutter(open_s: f64, duration_s: f64, shot_end_s: f64) -> ShutterInterval {
        ShutterInterval::resolve(
            open_s,
            duration_s,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(0.0, shot_end_s).expect("shot"),
        )
        .expect("shutter")
    }

    #[test]
    fn g0_endpoints_are_contained_for_arbitrary_rotation() {
        let local_bounds = Aabb::new(Point3::new(-2.0, -0.5, -1.0), Point3::new(1.0, 3.0, 0.25));
        let local = FiniteLocalAabb::try_new(local_bounds).expect("finite local bounds");
        let start = RigidTransform::identity();
        let end = rotation_z(core::f64::consts::PI - 1.0e-12, [4.0, -2.0, 0.5]);
        let segment = RigidMotionSegment::try_new(start, end, 2.0, 3.0, RotationSweep::Arbitrary)
            .expect("valid arbitrary-rotation segment");
        let swept = conservative_world_swept_aabb(local, segment).expect("finite swept bounds");
        assert_contains_transform(swept, local_bounds, start);
        assert_contains_transform(swept, local_bounds, end);
    }

    #[test]
    fn g0_interior_high_spin_and_near_pi_samples_are_contained() {
        let local_bounds = Aabb::new(Point3::new(0.25, -1.5, -0.2), Point3::new(2.0, 0.75, 0.8));
        let local = FiniteLocalAabb::try_new(local_bounds).expect("finite local bounds");
        let start = RigidTransform::identity();
        let end = rotation_z(core::f64::consts::PI - 1.0e-13, [3.0, -1.0, 2.0]);
        let segment = RigidMotionSegment::try_new(start, end, 0.0, 1.0, RotationSweep::Arbitrary)
            .expect("valid arbitrary-rotation segment");
        let swept = conservative_world_swept_aabb(local, segment).expect("finite swept bounds");
        for (fraction, angle) in [
            (0.125, 4.0 * core::f64::consts::PI),
            (0.5, 9.0 * core::f64::consts::PI),
            (0.875, core::f64::consts::PI - 1.0e-14),
        ] {
            let translation = [3.0 * fraction, -fraction, 2.0 * fraction];
            assert_contains_transform(swept, local_bounds, rotation_z(angle, translation));
        }
    }

    #[test]
    fn g0_quaternion_roundoff_cannot_escape_arbitrary_rotation_sphere() {
        let local_bounds = Aabb::new(
            Point3::new(
                -30.293_556_256_793_288,
                29.793_706_838_168_61,
                -17.989_728_812_426_2,
            ),
            Point3::new(
                -30.293_556_256_793_288,
                29.793_706_838_168_61,
                -17.989_728_812_426_2,
            ),
        );
        let local = FiniteLocalAabb::try_new(local_bounds).expect("finite point bound");
        let rotation = RigidTransform::try_new(
            [
                -0.0,
                0.470_420_147_449_030_55,
                0.779_086_783_908_792_3,
                -0.414_401_578_197_630_7,
            ],
            [0.0; 3],
        )
        .expect("admitted near-unit quaternion");
        let motion =
            RigidMotionSegment::try_new(rotation, rotation, 0.0, 1.0, RotationSweep::Arbitrary)
                .expect("arbitrary rotation path");
        let swept = conservative_world_swept_aabb(local, motion).expect("swept bound");

        assert_contains_transform(swept, local_bounds, rotation);
    }

    #[test]
    fn g0_pure_translation_matches_the_analytic_envelope() {
        let local_bounds = Aabb::new(Point3::new(-1.0, 2.0, -3.0), Point3::new(4.0, 5.0, 6.0));
        let local = FiniteLocalAabb::try_new(local_bounds).expect("finite local bounds");
        let start = transform([0.0, 0.0, 0.0, 1.0], [2.0, -4.0, 1.0]);
        let end = transform([0.0, 0.0, 0.0, 1.0], [7.0, 3.0, -2.0]);
        let segment = RigidMotionSegment::try_new(start, end, 5.0, 9.0, RotationSweep::Constant)
            .expect("pure translation segment");
        let swept = conservative_world_swept_aabb(local, segment).expect("finite swept bounds");
        let analytic_min = Point3::new(1.0, -2.0, -5.0);
        let analytic_max = Point3::new(11.0, 8.0, 7.0);
        assert!(swept.min.x <= analytic_min.x && swept.min.y <= analytic_min.y);
        assert!(swept.min.z <= analytic_min.z);
        assert!(swept.max.x >= analytic_max.x && swept.max.y >= analytic_max.y);
        assert!(swept.max.z >= analytic_max.z);
        assert!(analytic_min.x - swept.min.x <= f64::EPSILON);
        assert!(swept.max.x - analytic_max.x <= 2.0 * f64::EPSILON * analytic_max.x);
    }

    #[test]
    fn g0_runtime_hermite_overshoot_is_inside_trajectory_bound() {
        let local_bounds = Aabb::new(Point3::new(-0.25, -0.5, -0.1), Point3::new(0.75, 0.25, 0.2));
        let local = FiniteLocalAabb::try_new(local_bounds).expect("finite local bounds");
        let start = TransformKeyframe::try_new(0.0, RigidTransform::identity(), [10.0, -2.0, 0.0])
            .expect("start");
        let end = TransformKeyframe::try_new(1.0, RigidTransform::identity(), [-10.0, 2.0, 0.0])
            .expect("end");
        let trajectory = RigidTransformTrajectory::try_new(vec![start, end]).expect("trajectory");
        let swept = conservative_trajectory_swept_aabb(local, &trajectory, shutter(0.0, 1.0, 1.0))
            .expect("Hermite swept bound");

        let midpoint = trajectory.evaluate(0.5).expect("midpoint").transform;
        assert_eq!(midpoint.translation_m()[0].to_bits(), 2.5_f64.to_bits());
        assert!(
            swept.max.x > 2.5,
            "endpoint-only translation envelope was retained: {swept:?}"
        );
        for sample in 0..=256 {
            let time_s = f64::from(sample) / 256.0;
            assert_contains_transform(
                swept,
                local_bounds,
                trajectory.evaluate(time_s).expect("sample").transform,
            );
        }
    }

    #[test]
    fn g0_subnormal_velocity_with_huge_duration_remains_inside_trajectory_bound() {
        let local_bounds = Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0));
        let local = FiniteLocalAabb::try_new(local_bounds).expect("point bound");
        let duration_s = 1.0e300;
        let start = TransformKeyframe::try_new(
            0.0,
            RigidTransform::identity(),
            [f64::from_bits(1), 0.0, 0.0],
        )
        .expect("start");
        let end = TransformKeyframe::try_new(duration_s, RigidTransform::identity(), [0.0; 3])
            .expect("end");
        let trajectory = RigidTransformTrajectory::try_new(vec![start, end]).expect("trajectory");
        let swept = conservative_trajectory_swept_aabb(
            local,
            &trajectory,
            shutter(0.0, duration_s, duration_s),
        )
        .expect("Hermite swept bound");
        let midpoint = trajectory.evaluate(duration_s / 2.0).expect("midpoint");
        let displacement = midpoint.transform.translation_m()[0];

        assert!(
            displacement > 0.0,
            "regression setup must survive runtime arithmetic"
        );
        assert!(
            swept.contains(
                midpoint
                    .transform
                    .transform_point(Point3::new(0.0, 0.0, 0.0))
            ),
            "subnormal-velocity Hermite sample {displacement:e} escaped {swept:?}"
        );
    }

    #[test]
    fn g0_trajectory_bound_refuses_a_shutter_outside_keyframes() {
        let local = FiniteLocalAabb::try_new(Aabb::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ))
        .expect("point bound");
        let keyframes = vec![
            TransformKeyframe::try_new(0.0, RigidTransform::identity(), [0.0; 3]).expect("start"),
            TransformKeyframe::try_new(1.0, RigidTransform::identity(), [0.0; 3]).expect("end"),
        ];
        let trajectory = RigidTransformTrajectory::try_new(keyframes).expect("trajectory");
        assert_eq!(
            conservative_trajectory_swept_aabb(local, &trajectory, shutter(1.0, 1.0, 2.0)),
            Err(MotionBoundsError::ShutterOutsideTrajectory)
        );
    }

    #[test]
    fn g0_degenerate_boxes_remain_finite_and_conservative() {
        let origin = FiniteLocalAabb::try_new(Aabb::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ))
        .expect("origin point box");
        assert_eq!(origin.origin_radius_upper_m().to_bits(), 0.0_f64.to_bits());
        let start = RigidTransform::identity();
        let end = transform([0.0, 0.0, 0.0, 1.0], [1.0, 2.0, 3.0]);
        let segment = RigidMotionSegment::try_new(start, end, 0.0, 1.0, RotationSweep::Arbitrary)
            .expect("point motion");
        let swept = conservative_world_swept_aabb(origin, segment).expect("finite point sweep");
        assert_eq!(
            [swept.min.x, swept.min.y, swept.min.z].map(f64::to_bits),
            [0.0, 0.0, 0.0].map(f64::to_bits)
        );
        assert_eq!(
            [swept.max.x, swept.max.y, swept.max.z].map(f64::to_bits),
            [1.0, 2.0, 3.0].map(f64::to_bits)
        );

        let offset_point = Aabb::new(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0));
        let local = FiniteLocalAabb::try_new(offset_point).expect("offset point box");
        let swept = conservative_world_swept_aabb(local, segment).expect("offset point sweep");
        assert_contains_transform(
            swept,
            offset_point,
            rotation_z(core::f64::consts::FRAC_PI_2, [0.5, 1.0, 1.5]),
        );
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_refuse_explicitly() {
        let inverted = Aabb {
            min: Point3::new(1.0, 0.0, 0.0),
            max: Point3::new(-1.0, 0.0, 0.0),
        };
        assert_eq!(
            FiniteLocalAabb::try_new(inverted),
            Err(MotionBoundsError::InvertedLocalAabb { axis: "x" })
        );
        let non_finite = Aabb {
            min: Point3::new(f64::NAN, 0.0, 0.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        assert_eq!(
            FiniteLocalAabb::try_new(non_finite),
            Err(MotionBoundsError::NonFiniteLocalAabb { field: "min.x" })
        );
        let huge = Aabb::new(
            Point3::new(-f64::MAX, -f64::MAX, -f64::MAX),
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
        );
        assert_eq!(
            FiniteLocalAabb::try_new(huge),
            Err(MotionBoundsError::Unrepresentable {
                field: "local_origin_radius_m"
            })
        );

        let identity = RigidTransform::identity();
        assert_eq!(
            RigidMotionSegment::try_new(identity, identity, 1.0, 1.0, RotationSweep::Constant),
            Err(MotionBoundsError::InvalidTimeInterval)
        );
        let rotated = rotation_z(core::f64::consts::FRAC_PI_2, [0.0; 3]);
        assert_eq!(
            RigidMotionSegment::try_new(identity, rotated, 0.0, 1.0, RotationSweep::Constant),
            Err(MotionBoundsError::ConstantRotationMismatch)
        );
        assert!(RigidTransform::try_new([0.0; 4], [0.0; 3]).is_err());
    }

    #[test]
    fn g0_arbitrary_rotation_bound_is_deterministic() {
        let local = FiniteLocalAabb::try_new(Aabb::new(
            Point3::new(-1.0, -2.0, -3.0),
            Point3::new(4.0, 5.0, 6.0),
        ))
        .expect("finite local bounds");
        let motion = RigidMotionSegment::try_new(
            RigidTransform::identity(),
            rotation_z(1.75, [1.0, 2.0, 3.0]),
            0.0,
            0.25,
            RotationSweep::Arbitrary,
        )
        .expect("valid motion");
        let first = conservative_world_swept_aabb(local, motion).expect("first result");
        let second = conservative_world_swept_aabb(local, motion).expect("second result");
        assert_eq!(first, second);
    }

    #[test]
    fn helper_rotation_is_proper() {
        let rotated = rotation_z(core::f64::consts::FRAC_PI_2, [0.0; 3])
            .transform_vector(Vec3::new(1.0, 0.0, 0.0));
        assert!((rotated.x.abs()) <= 4.0 * f64::EPSILON);
        assert!((rotated.y - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert_eq!(rotated.z.to_bits(), 0.0_f64.to_bits());
    }
}

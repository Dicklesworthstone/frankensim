//! Analytic rigid-frame velocities, distinct from rotating vector components.

use super::{
    DISPLACEMENT_DIMS, FrameId, FrameMotion, FrameTree, RATE_DIMS, WORLD, axis_is_unit,
    motor_is_finite, quat_is_unit, rotation_about, vec3_is_finite,
};
use crate::ScenarioError;
use crate::signal::{Interp, TimeSignal};
use fs_cheb::{Cheb1, ChebError};
use fs_ga::{Motor, Quat, Vec3};
use fs_qty::Dims;
use std::ops::ControlFlow;

/// Explicit storage/work admission for a reusable frame path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrameKinematicsBudget {
    /// Maximum rows in the supplied frame tree, bounding path lookup work.
    pub max_frames: usize,
    /// Maximum dynamic signal scalars on the selected path, including retained
    /// Chebyshev derivative coefficients. Borrowed tables are counted once.
    pub max_history_scalars: usize,
}

impl Default for FrameKinematicsBudget {
    fn default() -> Self {
        Self {
            max_frames: 256,
            max_history_scalars: 1 << 20,
        }
    }
}

/// Refusal before publishing a prepared path or physical velocity sample.
#[derive(Debug, Clone, PartialEq)]
pub enum FrameKinematicsError {
    /// Invalid frame topology, motion data, signal structure or dimensions.
    Scenario(ScenarioError),
    /// The original Chebyshev derivative owner refused a physical rate.
    Derivative(ChebError),
    /// An explicit admission limit would be exceeded.
    Budget {
        /// Resource being admitted.
        resource: &'static str,
        /// Required count.
        needed: usize,
        /// Caller-supplied maximum.
        maximum: usize,
    },
    /// An allocation or count cannot be represented/admitted.
    Storage(&'static str),
    /// The clamped schedule has a corner or jump at the requested time.
    Nondifferentiable {
        /// Frame carrying the schedule.
        frame: FrameId,
        /// Physical query time in seconds.
        time_s: f64,
    },
    /// An input or computed physical quantity is nonfinite.
    NonFinite(&'static str),
    /// The caller cancelled; no partial path/sample is returned.
    Cancelled,
}

impl From<ScenarioError> for FrameKinematicsError {
    fn from(error: ScenarioError) -> Self {
        Self::Scenario(error)
    }
}

impl core::fmt::Display for FrameKinematicsError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "frame kinematics refused: {self:?}")
    }
}

impl std::error::Error for FrameKinematicsError {}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), FrameKinematicsError> {
    if checkpoint().is_break() {
        Err(FrameKinematicsError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite_vector(vector: Vec3, what: &'static str) -> Result<Vec3, FrameKinematicsError> {
    if vec3_is_finite(vector) {
        Ok(vector)
    } else {
        Err(FrameKinematicsError::NonFinite(what))
    }
}

fn slope(a: f64, b: f64, start: f64, end: f64) -> Result<f64, FrameKinematicsError> {
    let difference = b - a;
    let width = end - start;
    let rate = if difference.is_finite() && width.is_finite() {
        difference / width
    } else if width.is_finite() {
        let scale = a.abs().max(b.abs());
        (b / scale - a / scale) / (width / scale)
    } else {
        let scale = start.abs().max(end.abs());
        let numerator = if difference.is_finite() {
            difference / scale
        } else {
            b / scale - a / scale
        };
        numerator / (end / scale - start / scale)
    };
    if rate.is_finite() {
        Ok(rate)
    } else {
        Err(FrameKinematicsError::NonFinite("physical schedule rate"))
    }
}

#[derive(Debug)]
struct PreparedSignal<'a> {
    source: &'a TimeSignal,
    derivative: Option<Cheb1>,
    second_derivative: Option<Cheb1>,
}

impl<'a> PreparedSignal<'a> {
    fn new(
        source: &'a TimeSignal,
        expected: Dims,
        budget: FrameKinematicsBudget,
        used_scalars: &mut usize,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, FrameKinematicsError> {
        if source.dims() != expected {
            return Err(ScenarioError::Dimensions {
                context: "frame motion history".into(),
                expected: expected.0,
                got: source.dims().0,
            }
            .into());
        }
        let dynamic = source
            .validation_dynamic_scalars()
            .ok_or(FrameKinematicsError::Storage(
                "signal scalar count overflow",
            ))?;
        let (derivative_count, second_derivative_count) = match source {
            TimeSignal::Chebfun(profile) => {
                let n = profile.cheb.coeffs().len();
                (n.saturating_sub(1).max(1), n.saturating_sub(2).max(1))
            }
            _ => (0, 0),
        };
        let needed = used_scalars
            .checked_add(dynamic)
            .and_then(|count| count.checked_add(derivative_count))
            .and_then(|count| count.checked_add(second_derivative_count))
            .ok_or(FrameKinematicsError::Storage(
                "signal scalar count overflow",
            ))?;
        if needed > budget.max_history_scalars {
            return Err(FrameKinematicsError::Budget {
                resource: "history scalars",
                needed,
                maximum: budget.max_history_scalars,
            });
        }
        let mut violations = Vec::new();
        source.check_with_checkpoint("frame motion history", &mut violations, &mut |_| {
            poll(checkpoint)
        })?;
        if let Some(violation) = violations.into_iter().next() {
            return Err(ScenarioError::Evaluate {
                what: violation.what,
            }
            .into());
        }
        let derivative = match source {
            TimeSignal::Chebfun(profile) => Some(
                profile
                    .cheb
                    .try_differentiate_with_checkpoint(&mut || {
                        if checkpoint().is_break() {
                            Err(ChebError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })
                    .map_err(|error| match error {
                        ChebError::Cancelled => FrameKinematicsError::Cancelled,
                        other => FrameKinematicsError::Derivative(other),
                    })?,
            ),
            _ => None,
        };
        let second_derivative = match derivative.as_ref() {
            Some(first) => Some(
                first
                    .try_differentiate_with_checkpoint(&mut || {
                        if checkpoint().is_break() {
                            Err(ChebError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })
                    .map_err(|error| match error {
                        ChebError::Cancelled => FrameKinematicsError::Cancelled,
                        other => FrameKinematicsError::Derivative(other),
                    })?,
            ),
            None => None,
        };
        *used_scalars = needed;
        Ok(Self {
            source,
            derivative,
            second_derivative,
        })
    }

    fn sample(
        &self,
        frame: FrameId,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(f64, f64, f64), FrameKinematicsError> {
        let value = self
            .source
            .eval_prevalidated_with_checkpoint(time_s, &mut || poll(checkpoint))??
            .value;
        let corner = || FrameKinematicsError::Nondifferentiable { frame, time_s };
        let (rate, acceleration) = match self.source {
            TimeSignal::Constant(_) => (0.0, 0.0),
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => {
                if from.value == to.value || time_s < *t_start || time_s > *t_end {
                    (0.0, 0.0)
                } else if time_s == *t_start || time_s == *t_end {
                    return Err(corner());
                } else {
                    (slope(from.value, to.value, *t_start, *t_end)?, 0.0)
                }
            }
            TimeSignal::Table {
                times,
                values,
                interp,
                ..
            } => {
                let found = times.binary_search_by(|probe| {
                    if *probe == time_s {
                        std::cmp::Ordering::Equal
                    } else {
                        probe.total_cmp(&time_s)
                    }
                });
                match (*interp, found) {
                    (Interp::Hold, Ok(index))
                        if index > 0 && values[index] != values[index - 1] =>
                    {
                        return Err(corner());
                    }
                    (Interp::Hold, _) => (0.0, 0.0),
                    (Interp::Linear, Ok(index)) => {
                        let before = if index == 0 {
                            0.0
                        } else {
                            slope(
                                values[index - 1],
                                values[index],
                                times[index - 1],
                                times[index],
                            )?
                        };
                        let after = if index + 1 == times.len() {
                            0.0
                        } else {
                            slope(
                                values[index],
                                values[index + 1],
                                times[index],
                                times[index + 1],
                            )?
                        };
                        if before != after {
                            return Err(corner());
                        }
                        (before, 0.0)
                    }
                    (Interp::Linear, Err(index)) if index == 0 || index == times.len() => {
                        (0.0, 0.0)
                    }
                    (Interp::Linear, Err(index)) => (
                        slope(
                            values[index - 1],
                            values[index],
                            times[index - 1],
                            times[index],
                        )?,
                        0.0,
                    ),
                }
            }
            TimeSignal::Chebfun(profile) => {
                let (start, end) = profile.cheb.domain();
                if time_s < start || time_s > end {
                    (0.0, 0.0)
                } else {
                    let derivative = self.derivative.as_ref().expect("prepared Chebyshev rate");
                    let second = self
                        .second_derivative
                        .as_ref()
                        .expect("prepared Chebyshev acceleration");
                    let rate = derivative.eval_with_checkpoint(time_s, &mut || poll(checkpoint))?;
                    let acceleration =
                        second.eval_with_checkpoint(time_s, &mut || poll(checkpoint))?;
                    if (time_s == start || time_s == end)
                        && (rate != 0.0 || acceleration != 0.0)
                    {
                        return Err(corner());
                    }
                    (rate, acceleration)
                }
            }
        };
        if !rate.is_finite() {
            return Err(FrameKinematicsError::NonFinite("physical schedule rate"));
        }
        if !acceleration.is_finite() {
            return Err(FrameKinematicsError::NonFinite(
                "physical schedule acceleration",
            ));
        }
        Ok((value, rate, acceleration))
    }
}

#[derive(Debug)]
enum PreparedMotion<'a> {
    Fixed {
        pose: Motor,
        translation: Vec3,
    },
    Rotating {
        axis: [f64; 3],
        center: Vec3,
        rate: f64,
    },
    Tilt {
        axis: [f64; 3],
        center: Vec3,
        angle: PreparedSignal<'a>,
    },
    Translating {
        direction: Vec3,
        displacement: PreparedSignal<'a>,
    },
}

/// An immutable admitted path with borrowed histories and prepared analytic rates.
///
/// Tables are not cloned. Scalar history validation and Chebyshev
/// differentiation run once; samples do no allocation or frame-tree lookup.
/// The source borrow prevents mutation of a schedule while this path is used.
#[derive(Debug)]
pub struct PreparedFrameKinematics<'a> {
    frame: FrameId,
    path: Vec<(FrameId, PreparedMotion<'a>)>,
}

impl FrameTree {
    /// Prepare analytic rigid motion through arbitrary admitted ancestor chains.
    ///
    /// Supports fixed, rotating, translating and scheduled-tilt declarations.
    /// No finite differencing or constant-screw approximation is used. Budgets
    /// apply before history work/allocation; cancellation publishes no path.
    pub fn prepare_kinematics(
        &self,
        frame: FrameId,
        budget: FrameKinematicsBudget,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedFrameKinematics<'_>, FrameKinematicsError> {
        poll(checkpoint)?;
        if self.frames.len() > budget.max_frames {
            return Err(FrameKinematicsError::Budget {
                resource: "frame rows",
                needed: self.frames.len(),
                maximum: budget.max_frames,
            });
        }
        let mut path = Vec::new();
        path.try_reserve_exact(self.frames.len())
            .map_err(|_| FrameKinematicsError::Storage("frame path allocation failed"))?;
        let mut current = frame;
        let mut used_scalars = 0;
        while current != WORLD {
            poll(checkpoint)?;
            if path.len() == self.frames.len() {
                return Err(ScenarioError::Frame {
                    what: format!("cyclic frame chain reached from id {}", frame.0),
                }
                .into());
            }
            let node = self.find_unique(current)?;
            let invalid = || ScenarioError::Frame {
                what: format!("frame {} has invalid finite motion coordinates", current.0),
            };
            let prepared = match &node.motion {
                FrameMotion::Fixed {
                    orientation,
                    translation,
                } => {
                    if !quat_is_unit(*orientation) || !vec3_is_finite(*translation) {
                        return Err(invalid().into());
                    }
                    let pose = Motor::from_parts(*orientation, *translation);
                    if !motor_is_finite(&pose) {
                        return Err(invalid().into());
                    }
                    PreparedMotion::Fixed {
                        pose,
                        translation: *translation,
                    }
                }
                FrameMotion::Rotating { axis, center, rate } => {
                    if !axis_is_unit(axis) || !vec3_is_finite(*center) || !rate.value.is_finite() {
                        return Err(invalid().into());
                    }
                    if rate.dims != RATE_DIMS {
                        return Err(ScenarioError::Dimensions {
                            context: "frame angular rate".into(),
                            expected: RATE_DIMS.0,
                            got: rate.dims.0,
                        }
                        .into());
                    }
                    PreparedMotion::Rotating {
                        axis: *axis,
                        center: *center,
                        rate: rate.value,
                    }
                }
                FrameMotion::Tilt {
                    axis,
                    center,
                    angle,
                } => {
                    if !axis_is_unit(axis) || !vec3_is_finite(*center) {
                        return Err(invalid().into());
                    }
                    PreparedMotion::Tilt {
                        axis: *axis,
                        center: *center,
                        angle: PreparedSignal::new(
                            angle,
                            Dims::NONE,
                            budget,
                            &mut used_scalars,
                            checkpoint,
                        )?,
                    }
                }
                FrameMotion::Translating {
                    direction,
                    displacement,
                } => {
                    if !axis_is_unit(direction) {
                        return Err(invalid().into());
                    }
                    PreparedMotion::Translating {
                        direction: Vec3::new(direction[0], direction[1], direction[2]),
                        displacement: PreparedSignal::new(
                            displacement,
                            DISPLACEMENT_DIMS,
                            budget,
                            &mut used_scalars,
                            checkpoint,
                        )?,
                    }
                }
            };
            path.push((current, prepared));
            current = node.parent;
        }
        poll(checkpoint)?;
        Ok(PreparedFrameKinematics { frame, path })
    }
}

/// One rigid velocity field at a physical time, reusable for many points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameKinematics {
    frame: FrameId,
    time_s: f64,
    pose: Motor,
    origin: Vec3,
    linear: Vec3,
    angular: Vec3,
    linear_acceleration: Vec3,
    angular_acceleration: Vec3,
}

/// Absolute WORLD position and velocity of an explicitly specified local point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FramePointKinematics {
    /// Position in meters.
    pub position_world: Vec3,
    /// Inertial velocity in meters/second, including frame transport.
    pub velocity_world: Vec3,
}

/// Position, velocity and acceleration of a point moving relative to the frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FramePointMotion {
    pub position_world: Vec3,
    pub velocity_world: Vec3,
    pub acceleration_world: Vec3,
}

/// Point coordinates and their first two derivatives in the sampled source frame.
///
/// Relative velocity and acceleration are derivatives of moving-frame
/// coordinates, not merely inertial WORLD vectors expressed in different axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameRelativePointMotion {
    /// Position relative to the frame origin, in frame axes [m].
    pub position_frame: Vec3,
    /// Time derivative of `position_frame`, in frame axes [m/s].
    pub velocity_relative: Vec3,
    /// Time derivative of `velocity_relative`, in frame axes [m/s²].
    pub acceleration_relative: Vec3,
}

impl PreparedFrameKinematics<'_> {
    /// Source frame whose admitted path this object borrows.
    #[must_use]
    pub const fn source_frame(&self) -> FrameId {
        self.frame
    }

    /// Sample pose, frame-origin velocity/acceleration and angular
    /// velocity/acceleration once at `time_s`.
    ///
    /// A clamped schedule's corners/jumps refuse rather than selecting an
    /// undeclared temporal side. At a linear table knot the represented
    /// one-sided rates must agree; constant exterior intervals have zero rate.
    /// Queries allocate nothing, poll per frame and Chebyshev coefficient,
    /// and never clone/revalidate a history or advance any physical solver.
    pub fn sample(
        &self,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<FrameKinematics, FrameKinematicsError> {
        poll(checkpoint)?;
        if !time_s.is_finite() {
            return Err(FrameKinematicsError::NonFinite("physical time"));
        }
        let mut pose = Motor::identity();
        let mut origin = Vec3::default();
        let mut linear = Vec3::default();
        let mut angular = Vec3::default();
        let mut linear_acceleration = Vec3::default();
        let mut angular_acceleration = Vec3::default();
        for (frame, motion) in &self.path {
            poll(checkpoint)?;
            match motion {
                PreparedMotion::Fixed {
                    pose: local,
                    translation,
                } => {
                    let rotation = Quat::from_rotor(local);
                    origin = rotation.rotate(origin) + *translation;
                    linear = rotation.rotate(linear);
                    angular = rotation.rotate(angular);
                    linear_acceleration = rotation.rotate(linear_acceleration);
                    angular_acceleration = rotation.rotate(angular_acceleration);
                    pose = local.compose(&pose);
                }
                PreparedMotion::Translating {
                    direction,
                    displacement,
                } => {
                    let (value, rate, acceleration) =
                        displacement.sample(*frame, time_s, checkpoint)?;
                    let offset = direction.scale(value);
                    origin = origin + offset;
                    linear = linear + direction.scale(rate);
                    linear_acceleration =
                        linear_acceleration + direction.scale(acceleration);
                    pose = Motor::translator(offset.x, offset.y, offset.z).compose(&pose);
                }
                _ => {
                    let (axis, center, angle, rate, rate_acceleration) = match motion {
                        PreparedMotion::Rotating { axis, center, rate } => {
                            (*axis, *center, rate * time_s, *rate, 0.0)
                        }
                        PreparedMotion::Tilt {
                            axis,
                            center,
                            angle,
                        } => {
                            let (value, rate, acceleration) =
                                angle.sample(*frame, time_s, checkpoint)?;
                            (*axis, *center, value, rate, acceleration)
                        }
                        _ => unreachable!("fixed and translating motions handled above"),
                    };
                    if !angle.is_finite() {
                        return Err(FrameKinematicsError::NonFinite("rotation angle"));
                    }
                    let local = rotation_about(axis, center, angle);
                    let rotation = Quat::from_rotor(&local);
                    let offset = rotation.rotate(origin - center);
                    let spin_axis = Vec3::new(axis[0], axis[1], axis[2]);
                    let spin = spin_axis.scale(rate);
                    let spin_acceleration = spin_axis.scale(rate_acceleration);
                    let rotated_linear = rotation.rotate(linear);
                    let rotated_angular = rotation.rotate(angular);
                    linear_acceleration = rotation.rotate(linear_acceleration)
                        + spin.scale(2.0).cross(rotated_linear)
                        + spin_acceleration.cross(offset)
                        + spin.cross(spin.cross(offset));
                    angular_acceleration = rotation.rotate(angular_acceleration)
                        + spin.cross(rotated_angular)
                        + spin_acceleration;
                    origin = center + offset;
                    linear = rotated_linear + spin.cross(offset);
                    angular = rotated_angular + spin;
                    pose = local.compose(&pose);
                }
            }
            finite_vector(origin, "frame origin")?;
            finite_vector(linear, "frame-origin velocity")?;
            finite_vector(angular, "frame angular velocity")?;
            finite_vector(linear_acceleration, "frame-origin acceleration")?;
            finite_vector(angular_acceleration, "frame angular acceleration")?;
            if !motor_is_finite(&pose) {
                return Err(FrameKinematicsError::NonFinite("frame pose"));
            }
        }
        poll(checkpoint)?;
        Ok(FrameKinematics {
            frame: self.frame,
            time_s,
            pose,
            origin,
            linear,
            angular,
            linear_acceleration,
            angular_acceleration,
        })
    }
}

impl FrameKinematics {
    /// Source frame of this sample.
    #[must_use]
    pub const fn source_frame(&self) -> FrameId {
        self.frame
    }

    /// Physical sample time in seconds.
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Original native motor mapping this frame into WORLD.
    #[must_use]
    pub const fn world_from_frame(&self) -> Motor {
        self.pose
    }

    /// Frame-origin position in WORLD meters.
    #[must_use]
    pub const fn origin_world(&self) -> Vec3 {
        self.origin
    }

    /// Inertial velocity of the frame origin in WORLD meters/second.
    #[must_use]
    pub const fn origin_velocity_world(&self) -> Vec3 {
        self.linear
    }

    /// Angular velocity in WORLD radians/second.
    #[must_use]
    pub const fn angular_velocity_world(&self) -> Vec3 {
        self.angular
    }

    /// Inertial acceleration of the frame origin in WORLD meters/second².
    #[must_use]
    pub const fn origin_acceleration_world(&self) -> Vec3 {
        self.linear_acceleration
    }

    /// Angular acceleration in WORLD radians/second².
    #[must_use]
    pub const fn angular_acceleration_world(&self) -> Vec3 {
        self.angular_acceleration
    }

    /// Instantaneous acceleration of the rigid carrier material point that is
    /// currently at `point` in WORLD coordinates.
    pub fn material_acceleration_at_world_point(
        &self,
        point: Vec3,
    ) -> Result<Vec3, FrameKinematicsError> {
        finite_vector(point, "WORLD query point")?;
        let offset = finite_vector(point - self.origin, "WORLD point offset")?;
        finite_vector(
            self.linear_acceleration
                + self.angular_acceleration.cross(offset)
                + self.angular.cross(self.angular.cross(offset)),
            "WORLD point acceleration",
        )
    }

    /// Evaluate the rigid carrier velocity at a fixed WORLD query position.
    ///
    /// This is the Eulerian velocity field `v_origin + omega × (x - origin)`.
    /// The query position is not transported from an earlier material point.
    /// Geometry ownership/contact and representability of the supplied point
    /// remain with the caller; this does not move a wall or a mesh.
    pub fn velocity_at_world_point(&self, point: Vec3) -> Result<Vec3, FrameKinematicsError> {
        finite_vector(point, "WORLD query point")?;
        let offset = finite_vector(point - self.origin, "WORLD point offset")?;
        finite_vector(
            self.linear + self.angular.cross(offset),
            "WORLD point velocity",
        )
    }

    /// Transform a local point with explicitly supplied relative velocity and
    /// acceleration, retaining rigid-frame Euler, centripetal and Coriolis terms.
    pub fn point_motion_in_world(
        &self,
        point: Vec3,
        relative_velocity: Vec3,
        relative_acceleration: Vec3,
    ) -> Result<FramePointMotion, FrameKinematicsError> {
        finite_vector(point, "local point")?;
        finite_vector(relative_velocity, "local relative velocity")?;
        finite_vector(relative_acceleration, "local relative acceleration")?;
        let rotation = Quat::from_rotor(&self.pose);
        let offset = finite_vector(rotation.rotate(point), "rotated local point")?;
        let relative_velocity_world =
            finite_vector(rotation.rotate(relative_velocity), "rotated relative velocity")?;
        let position_world = finite_vector(self.origin + offset, "WORLD point position")?;
        let velocity_world = finite_vector(
            self.linear + self.angular.cross(offset) + relative_velocity_world,
            "WORLD point velocity",
        )?;
        let acceleration_world = finite_vector(
            self.linear_acceleration
                + self.angular_acceleration.cross(offset)
                + self.angular.cross(self.angular.cross(offset))
                + self.angular.scale(2.0).cross(relative_velocity_world)
                + rotation.rotate(relative_acceleration),
            "WORLD point acceleration",
        )?;
        Ok(FramePointMotion {
            position_world,
            velocity_world,
            acceleration_world,
        })
    }

    /// Express inertial WORLD point motion relative to this sampled frame.
    ///
    /// This is the inverse of [`Self::point_motion_in_world`] at the same time.
    /// With WORLD offset `r = point_world - origin` and relative velocity
    /// `u = velocity_world - V - omega cross r` expressed in WORLD axes,
    /// the returned acceleration is
    /// `R^T (acceleration_world - A - alpha cross r
    ///       - omega cross (omega cross r) - 2 omega cross u)`.
    ///
    /// Translation, Euler, centripetal and Coriolis transport are removed
    /// exactly once. No force, trajectory or time step is inferred. This
    /// constant-work query allocates nothing; nonfinite inputs or results
    /// refuse. Coordinate information already rounded out of an absolute
    /// WORLD position cannot be recovered.
    pub fn point_motion_in_frame(
        &self,
        point_world: Vec3,
        velocity_world: Vec3,
        acceleration_world: Vec3,
    ) -> Result<FrameRelativePointMotion, FrameKinematicsError> {
        finite_vector(point_world, "WORLD point position")?;
        finite_vector(velocity_world, "WORLD point velocity")?;
        finite_vector(acceleration_world, "WORLD point acceleration")?;
        let offset = finite_vector(point_world - self.origin, "WORLD point offset")?;
        let relative_velocity_world = finite_vector(
            velocity_world - self.linear - self.angular.cross(offset),
            "WORLD relative velocity",
        )?;
        let relative_acceleration_world = finite_vector(
            acceleration_world
                - self.linear_acceleration
                - self.angular_acceleration.cross(offset)
                - self.angular.cross(self.angular.cross(offset))
                - self.angular.scale(2.0).cross(relative_velocity_world),
            "WORLD relative acceleration",
        )?;
        let inverse_rotation = Quat::from_rotor(&self.pose).conjugate();
        Ok(FrameRelativePointMotion {
            position_frame: finite_vector(inverse_rotation.rotate(offset), "local point")?,
            velocity_relative: finite_vector(
                inverse_rotation.rotate(relative_velocity_world),
                "local relative velocity",
            )?,
            acceleration_relative: finite_vector(
                inverse_rotation.rotate(relative_acceleration_world),
                "local relative acceleration",
            )?,
        })
    }

    /// Transform local point coordinates and their explicitly relative velocity.
    ///
    /// Adds origin transport, angular transport and the rotated relative rate:
    /// `v_world = v_origin + omega_world × R point + R v_relative`.
    /// Existing absolute-vector evaluation remains a rotation-only operation.
    pub fn point_in_world(
        &self,
        point: Vec3,
        relative_velocity: Vec3,
    ) -> Result<FramePointKinematics, FrameKinematicsError> {
        let motion =
            self.point_motion_in_world(point, relative_velocity, Vec3::default())?;
        Ok(FramePointKinematics {
            position_world: motion.position_world,
            velocity_world: motion.velocity_world,
        })
    }
}

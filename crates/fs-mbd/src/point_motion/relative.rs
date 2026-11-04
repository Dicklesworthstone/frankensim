//! A moving point's local coordinate derivatives, transported by native motion.
//!
//! These are prescribed kinematic queries, not extra material mass, a new joint,
//! or an integration/constraint law. The carrier's original dynamics still owns
//! its pose, angular velocity and physical acceleration.

use super::{DynamicsError, PreparedRigidPointMotion, RigidPointMotion, Vec3};
use crate::articulated::{ArticulatedError, ArticulatedPointMotion, PreparedLinkPointMotion};
use fs_ga::Vec3 as SpatialVec3;

impl PreparedRigidPointMotion {
    /// Observe a point whose body coordinates change at this same instant.
    ///
    /// `arm_body` is COM-relative. The other inputs are its first and second
    /// derivatives in the rotating BODY coordinates, not world derivatives
    /// expressed in body axes. With `u = R * arm_dot` and `b = R * arm_ddot`,
    /// add `u` to the original rigid-point velocity and `2 omega cross u + b`
    /// to its physical acceleration. The native tangential and centripetal
    /// terms are retained exactly once.
    ///
    /// This can describe a sliding probe or a migrating source location. It
    /// does not supply the relative-motion law, reaction loads, or acceleration
    /// through an impact. Carrier angular fields still describe the BODY.
    /// Zero relative derivatives return the original fixed-point result exactly.
    pub fn relative_point_motion(
        &self,
        arm_body: Vec3,
        velocity_body: Vec3,
        acceleration_body: Vec3,
    ) -> Result<RigidPointMotion, DynamicsError> {
        if !velocity_body.is_finite() || !acceleration_body.is_finite() {
            return Err(DynamicsError::NonFinite("relative_point_coordinate_derivatives"));
        }
        let mut point = self.point_motion(arm_body)?;
        if velocity_body == Vec3::ZERO && acceleration_body == Vec3::ZERO {
            return Ok(point);
        }
        let rotation = self.pose().orientation();
        let slip = rotation.rotate_body_to_world(velocity_body);
        let relative_acceleration = rotation.rotate_body_to_world(acceleration_body);
        let coriolis = point.kinematics.angular_velocity_world.cross(slip).scale(2.0);
        point.kinematics.point_velocity_world = point.kinematics.point_velocity_world.add(slip);
        point.point_acceleration_world = point
            .point_acceleration_world
            .add(coriolis)
            .add(relative_acceleration);
        if !point.kinematics.point_velocity_world.is_finite()
            || !point.point_acceleration_world.is_finite()
        {
            return Err(DynamicsError::NonFinite("relative_point_world_motion"));
        }
        Ok(point)
    }
}

impl PreparedLinkPointMotion {
    /// Observe a point moving relative to the original articulated LINK origin.
    ///
    /// The inputs are the point's coordinates and their first/second derivatives
    /// in that rotating link frame. The native receipt already includes base,
    /// joint and gravity effects. Only `R * point_dot` and
    /// `2 omega cross (R * point_dot) + R * point_ddot` are added here; no spatial
    /// acceleration is reinterpreted and no gravity is added a second time.
    ///
    /// Query the actual child link instead when it already owns this motion:
    /// supplying that child's joint motion again as relative motion double-counts
    /// it. A migrating contact location is not the velocity of either material
    /// surface at contact and must not replace those velocities in friction laws.
    /// Carrier angular fields retain the link's rotation. This read-only query
    /// allocates nothing. Zero derivatives preserve the fixed-point result.
    pub fn relative_point_motion(
        &self,
        point_body: SpatialVec3,
        velocity_body: SpatialVec3,
        acceleration_body: SpatialVec3,
    ) -> Result<ArticulatedPointMotion, ArticulatedError> {
        let rotation = self.world_from_link().rotation();
        // The original checked rotation validates even relative inputs that
        // would otherwise be inactive on a stationary carrier.
        let slip = rotation.rotate(velocity_body)?;
        let relative_acceleration = rotation.rotate(acceleration_body)?;
        let mut point = self.point_motion(point_body)?;
        if velocity_body == SpatialVec3::default()
            && acceleration_body == SpatialVec3::default()
        {
            return Ok(point);
        }
        point.point_velocity_world = point.point_velocity_world + slip;
        point.point_acceleration_world = point.point_acceleration_world
            + point.angular_velocity_world.cross(slip).scale(2.0)
            + relative_acceleration;
        for (field, value) in [
            ("relative_point.velocity_world", point.point_velocity_world),
            ("relative_point.acceleration_world", point.point_acceleration_world),
        ] {
            for (index, component) in [value.x, value.y, value.z].iter().enumerate() {
                if !component.is_finite() {
                    return Err(ArticulatedError::NonFinite { field, index });
                }
            }
        }
        Ok(point)
    }
}

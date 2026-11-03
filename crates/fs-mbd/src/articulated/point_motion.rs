//! Physical material-point acceleration from the original articulated receipts.
//!
//! Prescribed-base RNEA/ABA stores gravity-shifted spatial acceleration, while
//! floating ABA stores physical spatial acceleration. Typed entry methods keep
//! those conventions separate before one shared Cartesian point-motion map.

use super::{
    ArticulatedError, ForwardDynamics, FreeFloatingForwardDynamics, InverseDynamics, Kinematics,
    origin_linear_acceleration_body, validate_vec3,
};
use fs_ga::{Se3, Twist, Vec3};

/// Instantaneous physical motion of a material point fixed in an articulated link.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedPointMotion {
    /// Point position in world coordinates, in metres.
    pub point_world: Vec3,
    /// Cartesian point velocity in world coordinates, in m/s.
    pub point_velocity_world: Vec3,
    /// Cartesian point acceleration in world coordinates, in m/s².
    /// This includes gravity and is not an ideal accelerometer's specific force.
    pub point_acceleration_world: Vec3,
    /// Link angular velocity in world coordinates, in rad/s.
    pub angular_velocity_world: Vec3,
    /// Link angular acceleration in world coordinates, in rad/s².
    pub angular_acceleration_world: Vec3,
}

impl ArticulatedPointMotion {
    /// Ideal accelerometer specific force `a_point - gravity`, in world m/s².
    ///
    /// Supply the same physical gravity used by the dynamics. A freely falling
    /// nonrotating point reads zero; a stationary supported point reads minus
    /// gravity. Rotating into the sensor's own axes remains the sensor owner's
    /// task. This adds no sensor bias, bandwidth, contact or measurement model.
    pub fn specific_acceleration_world(
        &self,
        gravity_world: Vec3,
    ) -> Result<Vec3, ArticulatedError> {
        validate_vec3(gravity_world, "point_motion.gravity_world")?;
        let result = self.point_acceleration_world - gravity_world;
        validate_vec3(result, "point_motion.specific_acceleration_world")?;
        Ok(result)
    }
}

/// Shared instantaneous link motion for any number of fixed material points.
///
/// Preparation and each query use bounded, allocation-free vector operations.
/// The captured poses, twists and accelerations must describe the same original
/// model, state and instant. Preparation checks lengths, the selected index and
/// selected vectors; it cannot establish that separately supplied receipts have
/// matching state provenance. Reprepare after any state or load change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedLinkPointMotion {
    world_from_link: Se3,
    origin_velocity_world: Vec3,
    origin_acceleration_world: Vec3,
    angular_velocity_world: Vec3,
    angular_acceleration_world: Vec3,
}

impl InverseDynamics {
    /// Prepare physical link motion from prescribed-base inverse dynamics.
    ///
    /// `gravity_world` must be the original RNEA gravity. The receipt stores
    /// apparent spatial acceleration with that gravity subtracted; this method
    /// restores it exactly once after the moving-frame Cartesian conversion.
    /// `kinematics` must use the same base, joint positions and joint velocities.
    pub fn prepare_point_motion(
        &self,
        kinematics: &Kinematics,
        link: usize,
        gravity_world: Vec3,
    ) -> Result<PreparedLinkPointMotion, ArticulatedError> {
        prepare(kinematics, &self.body_acceleration, link, gravity_world)
    }
}

impl ForwardDynamics {
    /// Prepare physical link motion from prescribed-base forward dynamics.
    ///
    /// Supply the original ABA gravity and kinematics at the same state. As in
    /// RNEA, `body_acceleration` has gravity subtracted; this method restores it
    /// once. It neither re-solves joint acceleration nor infers support reactions.
    pub fn prepare_point_motion(
        &self,
        kinematics: &Kinematics,
        link: usize,
        gravity_world: Vec3,
    ) -> Result<PreparedLinkPointMotion, ArticulatedError> {
        prepare(kinematics, &self.body_acceleration, link, gravity_world)
    }
}

impl FreeFloatingForwardDynamics {
    /// Prepare physical link motion from free-floating forward dynamics.
    ///
    /// This receipt already includes physical gravity, so no extra gravity
    /// input or addition is appropriate. Supply kinematics from the same base
    /// pose/twist and original joint state. The spatial linear component still
    /// needs `omega cross v` before it is a Cartesian origin acceleration.
    pub fn prepare_point_motion(
        &self,
        kinematics: &Kinematics,
        link: usize,
    ) -> Result<PreparedLinkPointMotion, ArticulatedError> {
        prepare(
            kinematics,
            &self.body_spatial_acceleration,
            link,
            Vec3::default(),
        )
    }
}

fn prepare(
    kinematics: &Kinematics,
    acceleration: &[Twist],
    link: usize,
    gravity_correction_world: Vec3,
) -> Result<PreparedLinkPointMotion, ArticulatedError> {
    let link_count = kinematics.world_from_link.len();
    for (field, got) in [
        ("point_motion.body_twist", kinematics.body_twist.len()),
        ("point_motion.body_acceleration", acceleration.len()),
    ] {
        if got != link_count {
            return Err(ArticulatedError::LengthMismatch {
                field,
                expected: link_count,
                got,
            });
        }
    }
    if link >= link_count {
        return Err(ArticulatedError::InvalidLinkIndex { link, link_count });
    }
    validate_vec3(gravity_correction_world, "point_motion.gravity_world")?;
    let world_from_link = kinematics.world_from_link[link];
    let rotation = world_from_link.rotation();
    let velocity = kinematics.body_twist[link];
    let spatial = acceleration[link];
    let origin_acceleration_body = origin_linear_acceleration_body(spatial, velocity)?;
    let origin_acceleration_world =
        rotation.rotate(origin_acceleration_body)? + gravity_correction_world;
    validate_vec3(
        origin_acceleration_world,
        "point_motion.origin_acceleration_world",
    )?;
    Ok(PreparedLinkPointMotion {
        world_from_link,
        origin_velocity_world: rotation.rotate(velocity.linear)?,
        origin_acceleration_world,
        angular_velocity_world: rotation.rotate(velocity.angular)?,
        angular_acceleration_world: rotation.rotate(spatial.angular)?,
    })
}

impl PreparedLinkPointMotion {
    /// Original world-from-link pose captured during preparation.
    #[must_use]
    pub const fn world_from_link(&self) -> Se3 {
        self.world_from_link
    }

    /// Evaluate a body-fixed point relative to the original link origin.
    ///
    /// For world arm `r`, `v = v_origin + omega cross r` and
    /// `a = a_origin + alpha cross r + omega cross (omega cross r)`.
    /// Joint Coriolis and centrifugal terms are already in the native RNEA/ABA
    /// acceleration and twist. They are retained, not separately reconstructed.
    /// A point sliding within the queried link needs its own relative motion;
    /// instead query the actual moving child link when it is in the tree.
    /// This is a smooth instantaneous query, not an integration or impact law.
    pub fn point_motion(
        &self,
        point_body: Vec3,
    ) -> Result<ArticulatedPointMotion, ArticulatedError> {
        validate_vec3(point_body, "point_motion.point_body")?;
        let arm_world = self.world_from_link.rotation().rotate(point_body)?;
        let rotational_velocity = self.angular_velocity_world.cross(arm_world);
        let point_velocity_world = self.origin_velocity_world + rotational_velocity;
        let point_acceleration_world = self.origin_acceleration_world
            + self.angular_acceleration_world.cross(arm_world)
            + self.angular_velocity_world.cross(rotational_velocity);
        validate_vec3(point_velocity_world, "point_motion.point_velocity_world")?;
        validate_vec3(
            point_acceleration_world,
            "point_motion.point_acceleration_world",
        )?;
        Ok(ArticulatedPointMotion {
            point_world: self.world_from_link.transform_point(point_body)?,
            point_velocity_world,
            point_acceleration_world,
            angular_velocity_world: self.angular_velocity_world,
            angular_acceleration_world: self.angular_acceleration_world,
        })
    }
}

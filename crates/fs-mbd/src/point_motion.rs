//! Instantaneous smooth motion of body-fixed material points.
//!
//! The original native pose, momentum, full inertia and configured uniform
//! gravity determine these queries. Preparation shares angular and centre-of-mass
//! quantities across any number of attachments without advancing the body.

use crate::{
    DynamicsError, MassProperties, PointKinematics, Pose, RigidBodyIntegrator, RigidBodyState,
    Vec3, Wrench,
};
use fs_math::dd::Dd;

mod relative;

/// Immutable instantaneous motion shared by points on one native rigid body.
///
/// Preparation evaluates the original full inertia once for angular velocity
/// and once for angular acceleration. Each subsequent point query uses a fixed
/// number of vector operations and allocates no storage. The captured state,
/// mass and inertia remain fixed; a new accepted state or wrench needs a new
/// preparation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedRigidPointMotion {
    pose: Pose,
    center_of_mass_velocity_world: Vec3,
    angular_velocity_body: Vec3,
    angular_velocity_world: Vec3,
    center_of_mass_acceleration_world: Vec3,
    angular_acceleration_body: Vec3,
    angular_acceleration_world: Vec3,
}

/// Kinematics and instantaneous acceleration of a body-fixed material point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidPointMotion {
    /// Native point position and velocity, including both angular-velocity frames.
    pub kinematics: PointKinematics,
    /// Centre-of-mass acceleration in m/s², including configured gravity.
    pub center_of_mass_acceleration_world: Vec3,
    /// Angular acceleration in rad/s², expressed in the actual body frame.
    pub angular_acceleration_body: Vec3,
    /// The same angular acceleration in the world frame.
    pub angular_acceleration_world: Vec3,
    /// Material-point acceleration in m/s² in the world frame.
    pub point_acceleration_world: Vec3,
}

impl RigidPointMotion {
    /// Ideal accelerometer specific force `a_point - gravity`, in world m/s².
    ///
    /// Supply the physical gravity used by the original dynamics, even when
    /// its weight was already included in an external wrench. A nonrotating
    /// freely falling point reads zero; a stationary supported point reads
    /// minus gravity. Sensor-axis rotation, bias and bandwidth belong to the
    /// sensor consumer. This is the same convention as articulated-point motion.
    pub fn specific_acceleration_world(self, gravity_world: Vec3) -> Result<Vec3, DynamicsError> {
        if !gravity_world.is_finite() || !self.point_acceleration_world.is_finite() {
            return Err(DynamicsError::NonFinite("specific_acceleration_world"));
        }
        let result = self.point_acceleration_world.sub(gravity_world);
        if !result.is_finite() {
            return Err(DynamicsError::NonFinite("specific_acceleration_world"));
        }
        Ok(result)
    }
}

impl RigidBodyIntegrator {
    /// Captures smooth instantaneous motion under the supplied native wrench.
    ///
    /// `wrench` has exactly the convention of [`Self::step`]: it excludes this
    /// integrator's configured gravity. It must contain all other actual forces
    /// and torques at this state and instant, including any reaction forces
    /// needed to describe the accepted motion. State-dependent loads must be
    /// evaluated at this pose before preparation.
    ///
    /// For constant mass and body-frame inertia, Euler's equation gives
    /// `alpha_body = I_body^-1 (torque_body - omega_body cross L_body)`.
    /// The original inertia owner retains nonprincipal body axes. This is a
    /// read-only derivative query, not a time step, constraint solve or impulse
    /// derivative; it does not assign a finite acceleration to an impact.
    /// Nonfinite inputs and unrepresentable derived vectors return an error.
    pub fn prepare_point_motion(
        self,
        state: RigidBodyState,
        properties: MassProperties,
        wrench: Wrench,
    ) -> Result<PreparedRigidPointMotion, DynamicsError> {
        self.gravity().validate()?;
        wrench.validate()?;
        let center = state.point_kinematics(properties, Vec3::ZERO)?;
        let angular_momentum_rate_body = wrench.torque_body.sub(
            center
                .angular_velocity_body
                .cross(state.angular_momentum_body()),
        );
        if !angular_momentum_rate_body.is_finite() {
            return Err(DynamicsError::NonFinite("angular_momentum_rate_body"));
        }
        // The checked native momentum-to-velocity map is the same full inverse
        // inertia action required by Euler's equation; no new tensor inverse.
        let angular_acceleration_body =
            properties.angular_velocity_body_checked(angular_momentum_rate_body)?;
        let angular_acceleration_world = state
            .pose()
            .orientation()
            .rotate_body_to_world(angular_acceleration_body);
        if !angular_acceleration_world.is_finite() {
            return Err(DynamicsError::NonFinite("angular_acceleration_world"));
        }
        let force = [
            wrench.force_world.x,
            wrench.force_world.y,
            wrench.force_world.z,
        ];
        let gravity = crate::verlet::lifted(self.gravity().acceleration_world());
        // Keep F/m's division remainder until gravity has been added. This is
        // the derivative of the steppers' compensated F + m*g resultant, while
        // avoiding an overflowing m*g when the acceleration itself is finite.
        let center_of_mass_acceleration_world = crate::verlet::rounded(
            core::array::from_fn(|i| {
                let quotient = force[i] / properties.mass();
                let remainder = (-quotient).mul_add(properties.mass(), force[i]);
                Dd::from_pair(quotient, remainder / properties.mass()) + gravity[i]
            }),
            "center_of_mass_acceleration_world",
        )?;
        if !center_of_mass_acceleration_world.is_finite() {
            return Err(DynamicsError::NonFinite(
                "center_of_mass_acceleration_world",
            ));
        }
        Ok(PreparedRigidPointMotion {
            pose: state.pose(),
            center_of_mass_velocity_world: center.center_of_mass_velocity_world,
            angular_velocity_body: center.angular_velocity_body,
            angular_velocity_world: center.angular_velocity_world,
            center_of_mass_acceleration_world,
            angular_acceleration_body,
            angular_acceleration_world,
        })
    }
}

impl PreparedRigidPointMotion {
    /// Returns the captured native centre-of-mass pose.
    #[must_use]
    pub const fn pose(&self) -> Pose {
        self.pose
    }

    /// Returns the captured angular velocity in the actual body frame.
    #[must_use]
    pub const fn angular_velocity_body(&self) -> Vec3 {
        self.angular_velocity_body
    }

    /// Evaluates a centre-of-mass-relative, body-fixed material point.
    ///
    /// With world arm `r`, velocity is `v_com + omega cross r`; acceleration is
    /// `a_com + alpha cross r + omega cross (omega cross r)`. The last term is
    /// present even during steady rotation. A sliding or deforming attachment
    /// needs its own relative motion and is outside this rigid-point query.
    pub fn point_motion(&self, arm_body: Vec3) -> Result<RigidPointMotion, DynamicsError> {
        if !arm_body.is_finite() {
            return Err(DynamicsError::NonFinite("arm_body"));
        }
        let arm_world = self.pose.orientation().rotate_body_to_world(arm_body);
        if !arm_world.is_finite() {
            return Err(DynamicsError::NonFinite("arm_world"));
        }
        let point_world = self.pose.point_world_from_body(arm_body)?;
        let rotational_velocity_world = self.angular_velocity_world.cross(arm_world);
        let point_velocity_world = self
            .center_of_mass_velocity_world
            .add(rotational_velocity_world);
        if !point_velocity_world.is_finite() {
            return Err(DynamicsError::NonFinite("point_velocity_world"));
        }
        let point_acceleration_world = self
            .center_of_mass_acceleration_world
            .add(self.angular_acceleration_world.cross(arm_world))
            .add(self.angular_velocity_world.cross(rotational_velocity_world));
        if !point_acceleration_world.is_finite() {
            return Err(DynamicsError::NonFinite("point_acceleration_world"));
        }
        Ok(RigidPointMotion {
            kinematics: PointKinematics {
                arm_body,
                arm_world,
                point_world,
                center_of_mass_velocity_world: self.center_of_mass_velocity_world,
                angular_velocity_body: self.angular_velocity_body,
                angular_velocity_world: self.angular_velocity_world,
                relative_point_velocity_world: Vec3::ZERO,
                point_velocity_world,
            },
            center_of_mass_acceleration_world: self.center_of_mass_acceleration_world,
            angular_acceleration_body: self.angular_acceleration_body,
            angular_acceleration_world: self.angular_acceleration_world,
            point_acceleration_world,
        })
    }
}

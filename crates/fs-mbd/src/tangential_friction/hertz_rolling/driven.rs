//! Actual frozen-body drive response and work for the joint Hertz/rolling phase.

use super::*;

#[derive(Clone, Debug)]
pub(super) struct HeldPairDrive {
    external: [Wrench; 2],
    pub acceleration: [f64; 5],
}

impl HeldPairDrive {
    pub fn new(
        core: &PreparedRigidHertzPair,
        radius: f64,
        external: [Wrench; 2],
    ) -> Result<Self, RigidHertzRollingError> {
        let mut point = [Vec3::ZERO; 2];
        let mut angular = [Vec3::ZERO; 2];
        for i in 0..2 {
            external[i].validate()?;
            let property = core.pair_properties[i];
            angular[i] = core.before[i].pose().orientation().rotate_body_to_world(
                property.angular_velocity_body_checked(external[i].torque_body)?,
            );
            point[i] = external[i]
                .force_world
                .scale(property.mass().recip())
                .add(angular[i].cross(core.points_before[i].arm_world));
            if !point[i].is_finite() || !angular[i].is_finite() {
                return Err(invalid("unrepresentable external drive response"));
            }
        }
        let point = point[0].sub(point[1]);
        let angular = angular[0].sub(angular[1]);
        let axes = core.axes;
        let acceleration = [
            point.dot(axes[0]),
            point.dot(axes[1]),
            core.contact_radius_m() * angular.dot(axes[2]),
            radius * angular.dot(axes[0]),
            radius * angular.dot(axes[1]),
        ];
        if acceleration.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable scaled external drive"));
        }
        Ok(Self {
            external,
            acceleration,
        })
    }

    /// One summed resultant update per actual body. The zero-drive path retains
    /// the original native publication arithmetic exactly.
    pub fn apply(
        &self,
        core: &PreparedRigidHertzPair,
        impulse: Vec3,
        couple: Vec3,
        gate: &CancelGate,
    ) -> Result<[RigidBodyState; 2], RigidHertzRollingError> {
        if self.external == [Wrench::ZERO; 2] {
            return Ok(apply_pair_wrench(
                core.before,
                core.pair_properties,
                core.arms,
                impulse,
                couple,
                gate,
            )?);
        }
        let mut after = core.before;
        let dt = core.duration_s();
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            poll(gate)?;
            let orientation = core.before[i].pose().orientation();
            let point_impulse = impulse.scale(sign);
            let linear = point_impulse.add(self.external[i].force_world.scale(dt));
            let angular = core.arms[i]
                .cross(orientation.rotate_world_to_body(point_impulse))
                .add(orientation.rotate_world_to_body(couple.scale(sign)))
                .add(self.external[i].torque_body.scale(dt));
            after[i] = core.before[i]
                .apply_momentum_impulse(core.pair_properties[i], linear, angular)?
                .state_after;
        }
        Ok(after)
    }

    /// Exact work of the declared constant external wrenches along the affine
    /// momentum path of this held-effort phase, including all common motion.
    pub fn work(
        &self,
        core: &PreparedRigidHertzPair,
        after: [PointKinematics; 2],
    ) -> Result<f64, RigidHertzRollingError> {
        let dt = core.duration_s();
        let mut work = 0.0;
        for i in 0..2 {
            let mean_v = core.points_before[i]
                .center_of_mass_velocity_world
                .scale(0.5)
                .add(after[i].center_of_mass_velocity_world.scale(0.5));
            let mean_w = core.points_before[i]
                .angular_velocity_body
                .scale(0.5)
                .add(after[i].angular_velocity_body.scale(0.5));
            work = finite(
                work + self.external[i].force_world.scale(dt).dot(mean_v)
                    + self.external[i].torque_body.scale(dt).dot(mean_w),
            )?;
        }
        Ok(work)
    }

    pub fn impulses(
        &self,
        core: &PreparedRigidHertzPair,
    ) -> Result<(Vec3, Vec3), RigidHertzRollingError> {
        let mut linear = Vec3::ZERO;
        let mut angular = Vec3::ZERO;
        let dt = core.duration_s();
        for i in 0..2 {
            let force_impulse = self.external[i].force_world.scale(dt);
            linear = linear.add(force_impulse);
            angular = angular
                .add(core.before[i].pose().position_world().cross(force_impulse))
                .add(
                    core.before[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(self.external[i].torque_body.scale(dt)),
                );
        }
        if !linear.is_finite() || !angular.is_finite() {
            return Err(invalid("unrepresentable external momentum transfer"));
        }
        Ok((linear, angular))
    }
}

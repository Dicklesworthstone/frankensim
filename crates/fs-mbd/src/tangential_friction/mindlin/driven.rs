//! Held native wrenches, simultaneous resultant impulses and independent work.

use super::*;

pub(super) struct HeldDrive {
    external: [Wrench; 2],
    pub acceleration: [f64; 2],
}

pub(super) struct NativeResponse {
    pub bodies: [RigidBodyState; 2],
    pub points: [PointKinematics; 2],
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub kinetic_energy_change_j: f64,
    pub contact_work_j: f64,
    pub external_work_j: f64,
    pub work_energy_residual_j: f64,
    pub linear_change_ns: Vec3,
    pub external_linear_ns: Vec3,
    pub external_angular_nms: Vec3,
}

impl HeldDrive {
    pub fn new(
        bodies: [RigidBodyState; 2],
        properties: [MassProperties; 2],
        points: [PointKinematics; 2],
        axes: [Vec3; 2],
        external: [Wrench; 2],
    ) -> Result<Self, RigidMindlinError> {
        let mut acceleration = [Vec3::ZERO; 2];
        for i in 0..2 {
            external[i].validate()?;
            let angular = bodies[i].pose().orientation().rotate_body_to_world(
                properties[i].angular_velocity_body_checked(external[i].torque_body)?,
            );
            acceleration[i] = external[i]
                .force_world
                .scale(properties[i].mass().recip())
                .add(angular.cross(points[i].arm_world));
            norm(acceleration[i])?;
        }
        let relative = acceleration[0].sub(acceleration[1]);
        Ok(Self {
            external,
            acceleration: [
                finite(relative.dot(axes[0]))?,
                finite(relative.dot(axes[1]))?,
            ],
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn apply(
        &self,
        before: [RigidBodyState; 2],
        properties: [MassProperties; 2],
        arms: [Vec3; 2],
        points: [PointKinematics; 2],
        impulse_on_b: Vec3,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NativeResponse, RigidMindlinError> {
        if self.external == [Wrench::ZERO; 2] {
            // Retain the established unforced arithmetic, including its
            // independent per-body work receipts and summed energy change.
            let native = apply_equal_and_opposite_impulse_at_body_points(
                before[0],
                properties[0],
                arms[0],
                before[1],
                properties[1],
                arms[1],
                impulse_on_b,
            )?;
            return Ok(NativeResponse {
                bodies: [native.first.state_after, native.second.state_after],
                points: [
                    native.first.point_kinematics_after,
                    native.second.point_kinematics_after,
                ],
                kinetic_before_j: finite(
                    native.first.work.kinetic_energy_before_j
                        + native.second.work.kinetic_energy_before_j,
                )?,
                kinetic_after_j: finite(
                    native.first.work.kinetic_energy_after_j
                        + native.second.work.kinetic_energy_after_j,
                )?,
                kinetic_energy_change_j: native.kinetic_energy_change_j,
                contact_work_j: native.impulse_work_j,
                external_work_j: 0.0,
                work_energy_residual_j: native.work_energy_residual_j,
                linear_change_ns: native.linear_momentum_change_world,
                external_linear_ns: Vec3::ZERO,
                external_angular_nms: Vec3::ZERO,
            });
        }
        let mut after = before;
        let mut after_points = points;
        let (mut kinetic_before, mut kinetic_after, mut kinetic_change, mut total_work) =
            (0.0, 0.0, 0.0, 0.0);
        let (mut contact_work, mut external_work) = (0.0, 0.0);
        let (mut external_linear, mut external_angular) = (Vec3::ZERO, Vec3::ZERO);
        for (i, sign) in [(0, -1.0), (1, 1.0)] {
            poll(cancelled)?;
            let orientation = before[i].pose().orientation();
            let contact_impulse = impulse_on_b.scale(sign);
            let force_impulse = self.external[i].force_world.scale(dt);
            let torque_impulse = self.external[i].torque_body.scale(dt);
            let angular = arms[i]
                .cross(orientation.rotate_world_to_body(contact_impulse))
                .add(torque_impulse);
            let native = before[i].apply_momentum_impulse(
                properties[i],
                contact_impulse.add(force_impulse),
                angular,
            )?;
            after[i] = native.state_after;
            after_points[i] = after[i].point_kinematics(properties[i], arms[i])?;
            let mean_point = points[i]
                .point_velocity_world
                .scale(0.5)
                .add(after_points[i].point_velocity_world.scale(0.5));
            let mean_com = points[i]
                .center_of_mass_velocity_world
                .scale(0.5)
                .add(after_points[i].center_of_mass_velocity_world.scale(0.5));
            let mean_omega = points[i]
                .angular_velocity_body
                .scale(0.5)
                .add(after_points[i].angular_velocity_body.scale(0.5));
            contact_work = finite(contact_work + contact_impulse.dot(mean_point))?;
            external_work = finite(
                external_work + force_impulse.dot(mean_com) + torque_impulse.dot(mean_omega),
            )?;
            kinetic_before = finite(kinetic_before + native.work.kinetic_energy_before_j)?;
            kinetic_after = finite(kinetic_after + native.work.kinetic_energy_after_j)?;
            kinetic_change = finite(kinetic_change + native.work.kinetic_energy_change_j)?;
            total_work = finite(total_work + native.work.impulse_work_j)?;
            external_linear = external_linear.add(force_impulse);
            external_angular = external_angular
                .add(before[i].pose().position_world().cross(force_impulse))
                .add(orientation.rotate_body_to_world(torque_impulse));
        }
        norm(external_linear)?;
        norm(external_angular)?;
        let linear_change = after[0]
            .linear_momentum_world()
            .add(after[1].linear_momentum_world())
            .sub(
                before[0]
                    .linear_momentum_world()
                    .add(before[1].linear_momentum_world()),
            );
        norm(linear_change)?;
        Ok(NativeResponse {
            bodies: after,
            points: after_points,
            kinetic_before_j: kinetic_before,
            kinetic_after_j: kinetic_after,
            kinetic_energy_change_j: kinetic_change,
            contact_work_j: contact_work,
            external_work_j: external_work,
            work_energy_residual_j: finite(kinetic_change - total_work)?,
            linear_change_ns: linear_change,
            external_linear_ns: external_linear,
            external_angular_nms: external_angular,
        })
    }
}

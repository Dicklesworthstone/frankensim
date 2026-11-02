//! One held wrench per native body, shared by all original contact rows.

use super::*;

#[derive(Clone, Debug)]
pub(super) struct HeldIslandDrive {
    external: Vec<Wrench>,
    /// COM and angular acceleration in world coordinates, in captured order.
    acceleration: Vec<[Vec3; 2]>,
    active: bool,
}

impl HeldIslandDrive {
    pub fn new(
        whole_body_count: usize,
        indices: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        external: Option<&[Wrench]>,
        gate: &CancelGate,
    ) -> Result<Self, RigidHertzRollingIslandError> {
        if let Some(wrenches) = external {
            if wrenches.len() != whole_body_count {
                return Err(invalid("one held wrench per original global body required"));
            }
            for (i, wrench) in wrenches.iter().enumerate() {
                poll(gate)?;
                wrench.validate()?;
                if *wrench != Wrench::ZERO && indices.binary_search(&i).is_err() {
                    return Err(invalid(
                        "nonparticipating body cannot receive drive in this contact assembly",
                    ));
                }
            }
        }
        let external: Vec<_> = indices
            .iter()
            .map(|&i| external.map_or(Wrench::ZERO, |w| w[i]))
            .collect();
        let active = external.iter().any(|w| *w != Wrench::ZERO);
        let mut acceleration = vec![[Vec3::ZERO; 2]; bodies.len()];
        if active {
            for i in 0..bodies.len() {
                poll(gate)?;
                acceleration[i] = [
                    external[i].force_world.scale(properties[i].mass().recip()),
                    bodies[i].pose().orientation().rotate_body_to_world(
                        properties[i].angular_velocity_body_checked(external[i].torque_body)?,
                    ),
                ];
                if acceleration[i].iter().any(|v| !v.is_finite()) {
                    return Err(invalid("unrepresentable native external acceleration"));
                }
            }
        }
        Ok(Self {
            external,
            acceleration,
            active,
        })
    }

    pub fn is_active(&self) -> bool {
        self.active
    }

    pub fn patch_acceleration(
        &self,
        bodies: [usize; 2],
        points: [PointKinematics; 2],
        axes: [Vec3; 3],
        contact_radius: f64,
        sphere_radius: f64,
    ) -> Result<[f64; 5], RigidHertzRollingIslandError> {
        if !self.active {
            return Ok([0.0; 5]);
        }
        let point: [Vec3; 2] = std::array::from_fn(|i| {
            let [linear, angular] = self.acceleration[bodies[i]];
            linear.add(angular.cross(points[i].arm_world))
        });
        let relative = point[0].sub(point[1]);
        let angular = self.acceleration[bodies[0]][1].sub(self.acceleration[bodies[1]][1]);
        let result = [
            relative.dot(axes[0]),
            relative.dot(axes[1]),
            contact_radius * angular.dot(axes[2]),
            sphere_radius * angular.dot(axes[0]),
            sphere_radius * angular.dot(axes[1]),
        ];
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "unrepresentable five-coordinate external acceleration",
            ));
        }
        Ok(result)
    }

    pub fn add_impulses(
        &self,
        linear: &mut [Vec3],
        angular: &mut [Vec3],
        dt: f64,
        gate: &CancelGate,
    ) -> Result<(), RigidHertzRollingIslandError> {
        if !self.active {
            return Ok(());
        }
        for i in 0..self.external.len() {
            poll(gate)?;
            linear[i] = linear[i].add(self.external[i].force_world.scale(dt));
            angular[i] = angular[i].add(self.external[i].torque_body.scale(dt));
            if !linear[i].is_finite() || !angular[i].is_finite() {
                return Err(invalid("unrepresentable summed external/contact resultant"));
            }
        }
        Ok(())
    }

    /// Contract each original body's wrench with its actual common midpoint
    /// velocity once, including motion that cancels from every contact row.
    pub fn work_and_impulses(
        &self,
        before: &[RigidBodyState],
        after: &[RigidBodyState],
        properties: &[MassProperties],
        dt: f64,
        gate: &CancelGate,
    ) -> Result<(f64, Vec3, Vec3), RigidHertzRollingIslandError> {
        if !self.active {
            return Ok((0.0, Vec3::ZERO, Vec3::ZERO));
        }
        let mut work = 0.0;
        let mut linear = Vec3::ZERO;
        let mut angular = Vec3::ZERO;
        for i in 0..before.len() {
            poll(gate)?;
            let mean_v = before[i]
                .center_of_mass_velocity_world(properties[i])?
                .scale(0.5)
                .add(
                    after[i]
                        .center_of_mass_velocity_world(properties[i])?
                        .scale(0.5),
                );
            let mean_omega = properties[i]
                .angular_velocity_body_checked(before[i].angular_momentum_body())?
                .scale(0.5)
                .add(
                    properties[i]
                        .angular_velocity_body_checked(after[i].angular_momentum_body())?
                        .scale(0.5),
                );
            let force_impulse = self.external[i].force_world.scale(dt);
            let torque_impulse = self.external[i].torque_body.scale(dt);
            work = finite(work + force_impulse.dot(mean_v) + torque_impulse.dot(mean_omega))?;
            linear = linear.add(force_impulse);
            angular = angular
                .add(before[i].pose().position_world().cross(force_impulse))
                .add(
                    before[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(torque_impulse),
                );
        }
        if !linear.is_finite() || !angular.is_finite() {
            return Err(invalid("unrepresentable external momentum transfer"));
        }
        Ok((work, linear, angular))
    }
}

//! Held-wrench mechanics and scalar events for the original collinear law.

use super::*;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Departure {
    external: [Wrench; 2],
    capacity: f64,
    gamma: f64,
    pub sign: f64,
}

#[derive(Clone, Debug)]
pub(super) struct HeldDrive {
    pub external: [Wrench; 2],
    pub acceleration: [f64; 2],
    equivalent_force: [f64; 2],
}
impl HeldDrive {
    pub fn new(
        pair: &PreparedRigidHertzRollingPair<'_>,
        response: [[f64; 2]; 2],
        external: [Wrench; 2],
    ) -> Result<Self, RigidHertzRollingCollinearError> {
        if external == [Wrench::ZERO; 2] {
            return Ok(Self {
                external,
                acceleration: [0.0; 2],
                equivalent_force: [0.0; 2],
            });
        }
        let mut point_acceleration = [Vec3::ZERO; 2];
        let mut angular_acceleration = [Vec3::ZERO; 2];
        for i in 0..2 {
            external[i].validate()?;
            let property = pair.core.pair_properties[i];
            angular_acceleration[i] = pair.core.before[i]
                .pose()
                .orientation()
                .rotate_body_to_world(
                    property.angular_velocity_body_checked(external[i].torque_body)?,
                );
            point_acceleration[i] = external[i]
                .force_world
                .scale(property.mass().recip())
                .add(angular_acceleration[i].cross(pair.core.points_before[i].arm_world));
            if !point_acceleration[i].is_finite() || !angular_acceleration[i].is_finite() {
                return Err(invalid("unrepresentable held-wrench acceleration"));
            }
        }
        let point = point_acceleration[0].sub(point_acceleration[1]);
        let angular = angular_acceleration[0].sub(angular_acceleration[1]);
        let axes = pair.core.axes;
        let point_roundoff = 128.0
            * f64::EPSILON
            * (magnitude(point_acceleration[0]) + magnitude(point_acceleration[1]));
        let angular_roundoff = 128.0
            * f64::EPSILON
            * (magnitude(angular_acceleration[0]) + magnitude(angular_acceleration[1]));
        if point.dot(axes[1]).abs() > point_roundoff
            || point.dot(axes[2]).abs() > point_roundoff
            || angular.dot(axes[0]).abs() > angular_roundoff
            || angular.dot(axes[2]).abs() > angular_roundoff
        {
            return Err(invalid(
                "held drive leaves the admitted shear/rolling subspace",
            ));
        }
        let acceleration = [
            finite(point.dot(axes[0]))?,
            finite(pair.radius_m * angular.dot(axes[1]))?,
        ];
        let [[a, b], [_b, c]] = response;
        let q1 = finite((acceleration[1] - (b / a) * acceleration[0]) / (c - b * (b / a)))?;
        let q0 = finite((acceleration[0] - b * q1) / a)?;
        Ok(Self {
            external,
            acceleration,
            equivalent_force: [q0, q1],
        })
    }

    pub fn apply(
        &self,
        before: [RigidBodyState; 2],
        duration: f64,
    ) -> Result<[RigidBodyState; 2], RigidHertzRollingCollinearError> {
        if self.external == [Wrench::ZERO; 2] {
            return Ok(before);
        }
        let mut after = before;
        for i in 0..2 {
            after[i] = RigidBodyState::new(
                before[i].pose(),
                before[i]
                    .linear_momentum_world()
                    .add(self.external[i].force_world.scale(duration)),
                before[i]
                    .angular_momentum_body()
                    .add(self.external[i].torque_body.scale(duration)),
            )?;
        }
        Ok(after)
    }

    /// Orthogonal mass-metric decomposition: free full-body work minus free
    /// relative work plus drive along the actual resolved contact trajectory.
    pub fn work(
        &self,
        before: [RigidBodyState; 2],
        properties: [MassProperties; 2],
        duration: f64,
        velocity: [f64; 2],
        displacement: [f64; 2],
    ) -> Result<f64, RigidHertzRollingCollinearError> {
        if self.external == [Wrench::ZERO; 2] {
            return Ok(0.0);
        }
        let mut free_work = 0.0;
        for i in 0..2 {
            let dp = self.external[i].force_world.scale(duration);
            let dl = self.external[i].torque_body.scale(duration);
            let mean_v = before[i]
                .linear_momentum_world()
                .add(dp.scale(0.5))
                .scale(properties[i].mass().recip());
            let mean_w = properties[i].angular_velocity_body_checked(
                before[i].angular_momentum_body().add(dl.scale(0.5)),
            )?;
            free_work = finite(free_work + dp.dot(mean_v) + dl.dot(mean_w))?;
        }
        let relative_correction = (0..2)
            .map(|i| {
                self.equivalent_force[i]
                    * (displacement[i]
                        - duration * (velocity[i] + 0.5 * duration * self.acceleration[i]))
            })
            .sum::<f64>();
        finite(free_work + relative_correction).map_err(Into::into)
    }

    pub fn impulses(&self, before: [RigidBodyState; 2], duration: f64) -> (Vec3, Vec3) {
        let mut linear = Vec3::ZERO;
        let mut angular = Vec3::ZERO;
        for (body, wrench) in before.into_iter().zip(self.external) {
            let dp = wrench.force_world.scale(duration);
            linear = linear.add(dp);
            angular = angular.add(body.pose().position_world().cross(dp)).add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(wrench.torque_body.scale(duration)),
            );
        }
        (linear, angular)
    }
}

pub(super) struct StaticExit {
    pub root: StopRoot,
    pub sign: f64,
}

impl PreparedRigidHertzRollingCollinear<'_> {
    pub(super) fn departure(&self, capacity: f64, sign: f64) -> Departure {
        Departure {
            external: self.drive.external,
            capacity,
            gamma: self.gamma,
            sign,
        }
    }

    pub(super) fn retained_departure(&self, capacity: f64) -> Option<f64> {
        let previous = self.history.stopped.as_ref()?.departure.as_ref()?;
        (previous.external == self.drive.external
            && previous.capacity == capacity
            && previous.gamma == self.gamma)
            .then_some(previous.sign)
    }

    pub(super) fn first_static_exit(
        &self,
        r: f64,
        capacity: f64,
        horizon: f64,
        gate: &CancelGate,
    ) -> Result<Option<StaticExit>, RigidHertzRollingCollinearError> {
        let [[a, b], [_b, c]] = self.response;
        let [du, dr] = self.drive.acceleration;
        let force_at = |t| -> Result<f64, RigidHertzRollingCollinearError> {
            let path = ExponentialRatePath::new(
                r,
                dr - (b / a) * du,
                self.gamma * (c - b * (b / a)),
                t,
                || gate.is_requested(),
            )?;
            Ok(finite(((b * self.gamma) * path.endpoint_rate() - du) / a)?)
        };
        let end = force_at(horizon)?;
        if end.abs() <= capacity {
            return Ok(None);
        }
        let direction = end.signum();
        let initial = finite((b * self.gamma * r - du) / a)?;
        if direction * initial >= capacity {
            return Ok(Some(StaticExit {
                root: StopRoot {
                    time: 0.0,
                    iterations: 0,
                    bracket: 0.0,
                },
                sign: -direction,
            }));
        }
        let tolerance = self
            .controls
            .mechanical
            .mechanical
            .contact
            .force_tolerance_n
            .min(128.0 * f64::EPSILON * (initial.abs() + end.abs() + capacity));
        let root = self.scalar_crossing(
            0.0,
            horizon,
            self.controls.maximum_root_iterations,
            |t| Ok((direction * force_at(t)? - capacity, tolerance)),
            gate,
        )?;
        Ok(Some(StaticExit {
            root,
            sign: -direction,
        }))
    }

    /// The sliding velocity derivative is monotone. Splitting at its one
    /// possible zero exposes a stop even when the interval endpoints share a sign.
    pub(super) fn first_driven_stop(
        &self,
        u: f64,
        r: f64,
        force: f64,
        sign: f64,
        horizon: f64,
        gate: &CancelGate,
    ) -> Result<Option<StopRoot>, RigidHertzRollingCollinearError> {
        let initial_derivative = self.slip_at(u, r, force, 0.0, gate)?.1;
        // The zero-viscosity law is linear, including arbitrarily small real
        // incoming slip. Its exact stop needs no time-bracket iterations.
        if self.gamma == 0.0 {
            let stop = -u / initial_derivative;
            return Ok(
                (u != 0.0 && stop > 0.0 && stop <= horizon).then_some(StopRoot {
                    time: stop,
                    iterations: 0,
                    bracket: 0.0,
                }),
            );
        }
        let end_derivative = self.slip_at(u, r, force, horizon, gate)?.1;
        let mut cuts = vec![0.0];
        let mut iterations = 0;
        if initial_derivative != 0.0
            && end_derivative != 0.0
            && initial_derivative.signum() != end_derivative.signum()
        {
            let derivative_sign = end_derivative.signum();
            let tolerance =
                128.0 * f64::EPSILON * (initial_derivative.abs() + end_derivative.abs());
            let extremum = self.scalar_crossing(
                0.0,
                horizon,
                self.controls.maximum_root_iterations,
                |t| {
                    Ok((
                        derivative_sign * self.slip_at(u, r, force, t, gate)?.1,
                        tolerance,
                    ))
                },
                gate,
            )?;
            iterations = extremum.iterations;
            cuts.push(extremum.time);
        }
        cuts.push(horizon);
        for pair in cuts.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            let (end, _, _) = self.slip_at(u, r, force, hi, gate)?;
            if sign * end > 0.0 {
                continue;
            }
            if lo == 0.0 && u == 0.0 {
                if sign * initial_derivative < 0.0 || sign * end < 0.0 {
                    return Err(invalid("unresolved outgoing driven slip direction"));
                }
                continue;
            }
            if end == 0.0 {
                return Ok(Some(StopRoot {
                    time: hi,
                    iterations,
                    bracket: 0.0,
                }));
            }
            let mut root = self.scalar_crossing(
                lo,
                hi,
                self.controls.maximum_root_iterations - iterations,
                |t| {
                    let (velocity, _, tolerance) = self.slip_at(u, r, force, t, gate)?;
                    Ok((-sign * velocity, tolerance))
                },
                gate,
            )?;
            root.iterations += iterations;
            return Ok(Some(root));
        }
        Ok(None)
    }

    /// A monotone first crossing is bracketed with f(lo)<=0<=f(hi).
    fn scalar_crossing(
        &self,
        mut lo: f64,
        mut hi: f64,
        maximum_iterations: usize,
        mut value: impl FnMut(f64) -> Result<(f64, f64), RigidHertzRollingCollinearError>,
        gate: &CancelGate,
    ) -> Result<StopRoot, RigidHertzRollingCollinearError> {
        for iteration in 1..=maximum_iterations {
            poll(gate)?;
            let t = lo + 0.5 * (hi - lo);
            let (f, tolerance) = value(t)?;
            if f == 0.0 || (hi - lo <= self.controls.time_tolerance_s && f.abs() <= tolerance) {
                return Ok(StopRoot {
                    time: t,
                    iterations: iteration,
                    bracket: hi - lo,
                });
            }
            if t == lo || t == hi {
                return Err(RigidHertzRollingCollinearError::NoProgress);
            }
            if f < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
        }
        Err(RigidHertzRollingCollinearError::RootBudget)
    }
}

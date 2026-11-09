//! Rigid tangential motion through the original halfway reflection.

use super::super::pressure::{GasPressure3, GasPressureField3};
use super::super::transport::MassTransport3;
use super::*;
use fs_math::dd::Dd;

mod pair;
pub use pair::PlanarRigidWallPairReport3;

/// Prescribed rigid material motion of a fixed, completely flooded domain plane.
///
/// Coordinates use cell centers at `index + 0.5`; domain planes are at zero or
/// the corresponding cell count. Angular speed is right handed about the
/// POSITIVE coordinate axis normal to the plane, on either domain side.
/// The origin is a fixed reference point, not an inferred body center.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarRigidWallMotion3 {
    /// Domain plane; no normal motion, tilting or cover/vacate occurs.
    pub wall: PlanarWall3,
    /// Reference point for velocity and every reported moment [lattice length].
    pub origin: [f64; 3],
    /// Tangential velocity at that reference point [lattice length/tick].
    pub velocity_at_origin: [f64; 3],
    /// Angular velocity about the positive plane-normal axis [rad/tick].
    pub angular_speed: f64,
}

impl PlanarRigidWallMotion3 {
    fn validate(self) -> Result<(), FreeSurfaceError3> {
        if !self.origin.into_iter().all(f64::is_finite)
            || !self.velocity_at_origin.into_iter().all(f64::is_finite)
            || !self.angular_speed.is_finite()
            || self.velocity_at_origin[self.wall.axis()] != 0.0
        {
            return Err(refused(
                None,
                "rigid planar motion requires finite tangential translation and axial spin",
            ));
        }
        Ok(())
    }

    fn at_arm(self, arm: [f64; 3]) -> [f64; 3] {
        // Retain exact old uniform-velocity population arithmetic at zero spin.
        if self.angular_speed == 0.0 {
            return self.velocity_at_origin;
        }
        let a = (self.wall.axis() + 1) % 3;
        let b = (self.wall.axis() + 2) % 3;
        let mut velocity = self.velocity_at_origin;
        velocity[a] -= self.angular_speed * arm[b];
        velocity[b] += self.angular_speed * arm[a];
        velocity
    }

    /// Evaluate the original rigid velocity at a declared lattice point.
    /// Refuses invalid motion, nonfinite positions or unrepresentable velocity.
    pub fn velocity_at(self, point: [f64; 3]) -> Result<[f64; 3], FreeSurfaceError3> {
        self.validate()?;
        let arm = core::array::from_fn(|d| point[d] - self.origin[d]);
        let velocity = self.at_arm(arm);
        if !arm.into_iter().chain(velocity).all(f64::is_finite) {
            return Err(refused(
                None,
                "rigid wall position or local velocity is unrepresentable",
            ));
        }
        Ok(velocity)
    }

    /// Maximum local speed over the finite rectangular domain plane.
    /// An affine velocity's norm attains its maximum at a corner. This is a
    /// geometric bound, not a claim about speeds on an infinite spinning plate.
    pub fn maximum_speed(self, cells: [usize; 3]) -> Result<f64, FreeSurfaceError3> {
        self.validate()?;
        if cells.contains(&0) {
            return Err(refused(None, "rigid wall requires a nonempty domain"));
        }
        let axis = self.wall.axis();
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let mut maximum = 0.0_f64;
        for corner in 0..4 {
            let mut point = [0.0; 3];
            point[axis] = if self.wall.upper() {
                cells[axis] as f64
            } else {
                0.0
            };
            point[a] = if corner & 1 != 0 {
                cells[a] as f64
            } else {
                0.0
            };
            point[b] = if corner & 2 != 0 {
                cells[b] as f64
            } else {
                0.0
            };
            maximum = maximum.max(speed(self.velocity_at(point)?));
        }
        if !maximum.is_finite() {
            return Err(refused(None, "rigid wall speed bound is unrepresentable"));
        }
        Ok(maximum)
    }
}

/// Actual link wrench, work and affine response of one original VOF timestep.
///
/// Raw forces and moments include numerical density-one pressure. Subtract the
/// separately reported reference wrench before interpreting physical loads.
/// Work is boundary mechanical exchange, not isothermal liquid energy or heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarRigidWallReport3 {
    /// Prescribed motion and common moment origin.
    pub motion: PlanarRigidWallMotion3,
    /// Actual raw population impulse delivered to liquid.
    pub fluid_impulse: [f64; 3],
    /// Exactly reciprocal raw wall impulse.
    pub wall_impulse: [f64; 3],
    /// Sum of actual halfway-link arms crossed with liquid impulses.
    pub fluid_angular_impulse: [f64; 3],
    /// Exactly reciprocal angular impulse about `motion.origin`.
    pub wall_angular_impulse: [f64; 3],
    /// Density-one pressure impulse into the liquid.
    pub reference_pressure_fluid_impulse: [f64; 3],
    /// Moment of that same reference pressure about `motion.origin`.
    pub reference_pressure_fluid_angular_impulse: [f64; 3],
    /// Derivative of `(J_a, J_b, L_axis)` into liquid with respect to
    /// `(V_a, V_b, angular_speed)`, at the SAME fixed origin. Here
    /// `a=(axis+1)%3`, `b=(axis+2)%3`. Accumulated from each original reflected
    /// link's `2 w rho / cs²`; symmetric positive semidefinite up to roundoff.
    /// It predicts the raw wrench for another admitted motion from THIS input
    /// fluid only, not later steps or different origins. No fluid trial runs
    /// during response evaluation, and native replay still owns acceptance.
    pub fluid_response: [[f64; 3]; 3],
    /// Sum of local rigid velocity dotted with actual liquid impulse.
    pub fluid_work: f64,
    /// Exactly reciprocal actuator work.
    pub wall_work: f64,
    /// Complete flooded plane area in lattice spacings squared.
    pub wetted_area: f64,
    /// Geometric centroid of that complete plane, in lattice coordinates.
    pub wetted_centroid: [f64; 3],
    /// Original liquid mass before and after the candidate.
    pub mass_before: f64,
    pub mass_after: f64,
    /// Absolute relative native liquid mass defect.
    pub relative_mass_defect: f64,
    /// Maximum physical wet-cell speed at either endpoint.
    pub maximum_speed: f64,
    /// Maximum prescribed local wall speed over the finite plane.
    pub maximum_wall_speed: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) struct RigidWallLink3 {
    pub(super) velocity: [f64; 3],
    arm: [f64; 3],
}

pub(super) struct RigidWallExchange3 {
    motion: PlanarRigidWallMotion3,
    angular: [Dd; 3],
    work: Dd,
    response: [[Dd; 3]; 3],
}

impl RigidWallExchange3 {
    fn new(motion: PlanarRigidWallMotion3) -> Self {
        Self {
            motion,
            angular: [Dd::ZERO; 3],
            work: Dd::ZERO,
            response: [[Dd::ZERO; 3]; 3],
        }
    }

    pub(super) fn link(
        &self,
        source: [i64; 3],
        direction: [f64; 3],
    ) -> Result<RigidWallLink3, FreeSurfaceError3> {
        let arm = core::array::from_fn(|d| {
            source[d] as f64 + 0.5 + 0.5 * direction[d] - self.motion.origin[d]
        });
        let velocity = self.motion.at_arm(arm);
        if !arm.into_iter().chain(velocity).all(f64::is_finite) {
            return Err(refused(
                None,
                "rigid wall link position or speed overflowed",
            ));
        }
        Ok(RigidWallLink3 { velocity, arm })
    }

    pub(super) fn record(
        &mut self,
        link: RigidWallLink3,
        direction: [f64; 3],
        strength: f64,
        derivative: f64,
    ) -> Result<(), FreeSurfaceError3> {
        let moment = cross(link.arm, direction);
        let axis = self.motion.wall.axis();
        let generalized = [
            direction[(axis + 1) % 3],
            direction[(axis + 2) % 3],
            moment[axis],
        ];
        for d in 0..3 {
            self.angular[d] = self.angular[d] + Dd::from_f64(moment[d]) * Dd::from_f64(strength);
            self.work = self.work
                + Dd::from_f64(link.velocity[d])
                    * Dd::from_f64(direction[d])
                    * Dd::from_f64(strength);
            for j in d..3 {
                self.response[d][j] = self.response[d][j]
                    + Dd::from_f64(derivative)
                        * Dd::from_f64(generalized[d])
                        * Dd::from_f64(generalized[j]);
            }
        }
        if !self
            .angular
            .into_iter()
            .chain([self.work])
            .chain(self.response.into_iter().flatten())
            .all(|value| value.to_f64().is_finite())
        {
            return Err(refused(
                None,
                "rigid wall angular impulse, work or response overflowed",
            ));
        }
        Ok(())
    }
}

struct PreparedRigidWall3 {
    exchange: PlanarWallStep3,
    wetted_area: f64,
    mass_before: f64,
    before_speed: f64,
    maximum_wall_speed: f64,
}

impl FreeSurface3 {
    /// Propose spatial rigid wall motion through one original liquid step.
    ///
    /// The full plane must remain flooded. Only tangential translation and
    /// normal-axis spin are supported; the fixed geometry cannot tilt, change
    /// coverage or move normally. Every actual reflection uses the local
    /// velocity, and the same link supplies reciprocal impulse, moment and work.
    /// All trials begin from this unchanged liquid. Cancellation or any native,
    /// local speed, population, mass or representability refusal publishes none.
    pub fn planar_rigid_wall_motion_candidate(
        &self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarRigidWallReport3), FreeSurfaceError3> {
        self.planar_rigid_wall_candidate_with_acceleration(
            motion,
            maximum_relative_mass_defect,
            None,
            gate,
        )
    }

    /// Propose rigid wall motion and a body-acceleration endpoint in ONE native
    /// collision, gas reconstruction, moving reflection and conversion step.
    ///
    /// First bind the initial uniform/affine acceleration explicitly. Collision
    /// uses its retained local force `F_n`; the accepted physical velocity and
    /// gas pressure use `F_next = rho_next * next(x)`, including newly wet cells.
    /// No force reframe discards the changing load's physical half impulse.
    /// The original reflection alone supplies the reported wall impulse,
    /// angular impulse and work; body loading is not charged to the wall.
    ///
    /// The motion is the prescribed interval wall velocity; its geometric plane
    /// remains fixed and completely flooded. Translation must be tangential
    /// and rotation axial. Supply the appropriate interval/midpoint velocity
    /// when the actuator accelerates. This does not evolve a finite-inertia
    /// body or supply an isothermal energy/heat balance. The actual final wet
    /// velocities, including the next half source, retain native Mach limits.
    /// Any source, boundary, conversion or cancellation refusal changes no
    /// accepted population, acceleration history or mass.
    pub fn planar_rigid_wall_motion_with_body_acceleration_candidate(
        &self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        next: BodyAccelerationField3,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarRigidWallReport3), FreeSurfaceError3> {
        self.planar_rigid_wall_candidate_with_acceleration(
            motion,
            maximum_relative_mass_defect,
            Some(next),
            gate,
        )
    }

    fn planar_rigid_wall_candidate_with_acceleration(
        &self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        next: Option<BodyAccelerationField3>,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarRigidWallReport3), FreeSurfaceError3> {
        let mut prepared =
            self.prepare_planar_rigid_wall_motion(motion, maximum_relative_mass_defect, gate)?;
        let mut candidate = match next {
            Some(next) => self.body_acceleration_endpoint_candidate(next, gate)?,
            None => self.clone(),
        };
        candidate.step_inner_with_wall(gate, None, None, None, Some(&mut prepared.exchange))?;
        let report = candidate.finish_planar_rigid_wall_motion(
            motion,
            maximum_relative_mass_defect,
            prepared,
            gate,
        )?;
        Ok((candidate, report))
    }

    /// The thermal owner supplies its one private candidate. The same native
    /// rigid reflection records the actual collision/conversion mass routes.
    /// Gas traction, when supplied, acts on its own missing populations in the
    /// same tick and does not become a solid-wall impulse or heat source.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_planar_rigid_wall_transport(
        &mut self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        maximum_routes: usize,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
        gate: &CancelGate,
    ) -> Result<(MassTransport3, PlanarRigidWallReport3), FreeSurfaceError3> {
        check_cancel(gate)?;
        let pressure = gas_pressure
            .map(|(bound, field)| GasPressure3::new(bound, field))
            .transpose()?;
        if let Some(forces) = forces {
            self.rebind_force_field(forces, gate)?;
        }
        let mut prepared =
            self.prepare_planar_rigid_wall_motion(motion, maximum_relative_mass_defect, gate)?;
        let mut trace = MassTransport3::new(self.mass_snapshot(), maximum_routes);
        self.step_inner_with_wall(
            gate,
            Some(&mut trace),
            Some(tau),
            pressure.as_ref(),
            Some(&mut prepared.exchange),
        )?;
        trace.after = self.mass_snapshot();
        let report = self.finish_planar_rigid_wall_motion(
            motion,
            maximum_relative_mass_defect,
            prepared,
            gate,
        )?;
        Ok((trace, report))
    }

    fn prepare_planar_rigid_wall_motion(
        &self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PreparedRigidWall3, FreeSurfaceError3> {
        check_cancel(gate)?;
        let cells = self.domain_cells();
        let maximum_wall_speed = motion.maximum_speed(cells)?;
        if maximum_wall_speed / CS2_3.sqrt() >= MACH_LIMIT {
            return Err(refused(
                None,
                "local rigid wall speed exceeds the lattice Mach domain",
            ));
        }
        if !maximum_relative_mass_defect.is_finite() || maximum_relative_mass_defect < 0.0 {
            return Err(refused(
                None,
                "rigid wall relative mass tolerance must be finite and nonnegative",
            ));
        }
        let (wetted_area, before_speed) = self.admit_planar_rigid_wall_state(motion, gate)?;
        let mass_before = self.ledger_mass();
        if !mass_before.is_finite() || mass_before <= 0.0 {
            return Err(refused(
                None,
                "rigid wall requires finite positive liquid mass",
            ));
        }
        Ok(PreparedRigidWall3 {
            exchange: PlanarWallStep3 {
                wall: motion.wall,
                extent: cells[motion.wall.axis()] as i64,
                velocity: motion.velocity_at_origin,
                impulse: [0.0; 3],
                rigid: Some(RigidWallExchange3::new(motion)),
                field: None,
            },
            wetted_area,
            mass_before,
            before_speed,
            maximum_wall_speed,
        })
    }

    fn admit_planar_rigid_wall_state(
        &self,
        motion: PlanarRigidWallMotion3,
        gate: &CancelGate,
    ) -> Result<(f64, f64), FreeSurfaceError3> {
        let admit = |p: [usize; 3], liquid: [f64; 3]| {
            for q in 1..Q3 {
                let direction = [E3[q].0 as f64, E3[q].1 as f64, E3[q].2 as f64];
                let normal = direction[motion.wall.axis()];
                if normal != if motion.wall.upper() { -1.0 } else { 1.0 } {
                    continue;
                }
                let point = core::array::from_fn(|d| p[d] as f64 + 0.5 - 0.5 * direction[d]);
                let velocity = motion.velocity_at(point)?;
                let relative = core::array::from_fn(|d| liquid[d] - velocity[d]);
                if !speed(relative).is_finite() || speed(relative) / CS2_3.sqrt() >= MACH_LIMIT {
                    return Err(refused(
                        Some(p),
                        "local liquid/rigid-wall speed exceeds the lattice Mach domain",
                    ));
                }
            }
            Ok(())
        };
        self.admit_planar_wall_state_with(motion.wall, gate, admit)
    }

    /// Coupled thermal force/wetting refresh must admit the actual final state.
    pub(crate) fn refresh_planar_rigid_wall_endpoint(
        &self,
        report: &mut PlanarRigidWallReport3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        let (_, after_speed) = self.admit_planar_rigid_wall_state(report.motion, gate)?;
        let mass_before = report.mass_before;
        let mass_after = self.ledger_mass();
        let relative_mass_defect = ((mass_after - mass_before) / mass_before).abs();
        if !mass_after.is_finite()
            || !relative_mass_defect.is_finite()
            || relative_mass_defect > maximum_relative_mass_defect
        {
            return Err(refused(
                None,
                "rigid wall exceeds its liquid mass-conservation tolerance",
            ));
        }
        report.mass_after = mass_after;
        report.relative_mass_defect = relative_mass_defect;
        report.maximum_speed = report.maximum_speed.max(after_speed);
        check_cancel(gate)?;
        Ok(())
    }

    fn finish_planar_rigid_wall_motion(
        &self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        prepared: PreparedRigidWall3,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallReport3, FreeSurfaceError3> {
        let PreparedRigidWall3 {
            exchange,
            wetted_area,
            mass_before,
            before_speed,
            maximum_wall_speed,
        } = prepared;
        let cells = self.domain_cells();
        let axis = motion.wall.axis();
        let mut reference_pressure_fluid_impulse = [0.0; 3];
        reference_pressure_fluid_impulse[axis] =
            (if motion.wall.upper() { -CS2_3 } else { CS2_3 }) * wetted_area;
        let mut wetted_centroid = cells.map(|n| n as f64 * 0.5);
        wetted_centroid[axis] = if motion.wall.upper() {
            cells[axis] as f64
        } else {
            0.0
        };
        let reference_pressure_fluid_angular_impulse = cross(
            core::array::from_fn(|d| wetted_centroid[d] - motion.origin[d]),
            reference_pressure_fluid_impulse,
        );
        let rigid = exchange
            .rigid
            .expect("rigid exchange retained by original streaming");
        let fluid_angular_impulse = rigid.angular.map(Dd::to_f64);
        let fluid_response = core::array::from_fn(|i| {
            core::array::from_fn(|j| rigid.response[i.min(j)][i.max(j)].to_f64())
        });
        let fluid_work = rigid.work.to_f64();
        if !reference_pressure_fluid_angular_impulse
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(refused(
                None,
                "rigid wall reference-pressure moment overflowed",
            ));
        }
        let mut report = PlanarRigidWallReport3 {
            motion,
            fluid_impulse: exchange.impulse,
            wall_impulse: exchange.impulse.map(|v| -v),
            fluid_angular_impulse,
            wall_angular_impulse: fluid_angular_impulse.map(|v| -v),
            reference_pressure_fluid_impulse,
            reference_pressure_fluid_angular_impulse,
            fluid_response,
            fluid_work,
            wall_work: -fluid_work,
            wetted_area,
            wetted_centroid,
            mass_before,
            mass_after: mass_before,
            relative_mass_defect: 0.0,
            maximum_speed: before_speed,
            maximum_wall_speed,
        };
        self.refresh_planar_rigid_wall_endpoint(&mut report, maximum_relative_mass_defect, gate)?;
        Ok(report)
    }

    /// Publish the same fully admitted rigid-wall candidate atomically.
    pub fn step_with_planar_rigid_wall_motion(
        &mut self,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallReport3, FreeSurfaceError3> {
        let (candidate, report) =
            self.planar_rigid_wall_motion_candidate(motion, maximum_relative_mass_defect, gate)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }

    /// Publish the fully admitted rigid-wall/body-acceleration candidate.
    /// Equal acceleration endpoints preserve the original held-field step.
    pub fn advance_body_acceleration_with_planar_rigid_wall_motion(
        &mut self,
        next: BodyAccelerationField3,
        motion: PlanarRigidWallMotion3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallReport3, FreeSurfaceError3> {
        let (candidate, report) = self.planar_rigid_wall_motion_with_body_acceleration_candidate(
            motion,
            maximum_relative_mass_defect,
            next,
            gate,
        )?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;

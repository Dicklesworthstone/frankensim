//! Joint original-body wrench equilibrium and independent physical publication gates.
use super::*;
use fs_mbd::StepReceipt;

#[derive(Clone, Copy)]
pub(super) struct Work {
    evaluations: usize,
    material_evaluations: usize,
}

#[derive(Clone, Copy)]
pub(super) struct Evaluation {
    pub(super) receipts: [StepReceipt; 2],
    residual: [f64; 12],
    expected: [f64; 12],
    contact: [RigidCohesiveLoad; 2],
    cayley: [Vec3; 2],
    storage: f64,
    fracture: f64,
    viscous: f64,
    spatial_work: f64,
}

impl RigidCohesiveSystem<'_> {
    pub(super) fn solve(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        gate: Option<&CancelGate>,
    ) -> Result<(Evaluation, [Wrench; 2], usize, Work), RigidCohesiveError> {
        let mut work = Work {
            evaluations: 0,
            material_evaluations: 0,
        };
        let mut x = pack(self.last_wrenches);
        let mut value = self.evaluate(x, loads, &mut work, gate)?;
        let mut matrix = [0.0; 144];
        let mut rhs = [0.0; 12];
        let mut correction = [0.0; 12];
        let scale: [f64; 12] = core::array::from_fn(|j| {
            let (absolute, maximum) = if j % 6 < 3 {
                (
                    self.config.force_absolute_tolerance_n,
                    self.config.maximum_total_force_n,
                )
            } else {
                (
                    self.config.torque_absolute_tolerance_n_m,
                    self.config.maximum_total_torque_n_m,
                )
            };
            (absolute / self.config.relative_tolerance)
                .min(maximum)
                .max(absolute)
        });
        for iteration in 0..=self.config.maximum_iterations {
            poll(gate)?;
            let error = self.scaled_residual(x, value)?;
            if error <= 1.0 {
                return Ok((value, unpack(x), iteration, work));
            }
            if iteration == self.config.maximum_iterations {
                return Err(RigidCohesiveError::NotConverged {
                    scaled_residual: error,
                    iterations: iteration,
                });
            }
            for j in 0..12 {
                // cbrt(binary64 epsilon): central differentiation balances
                // truncation and subtraction without a physical stiffness.
                let increment = 6.055454452393343e-6 * x[j].abs().max(scale[j]);
                let mut high_x = x;
                let mut low_x = x;
                high_x[j] += increment;
                low_x[j] -= increment;
                let high = self.evaluate(high_x, loads, &mut work, gate);
                if let Err(ref refusal) = high {
                    terminal(refusal)?;
                }
                let low = self.evaluate(low_x, loads, &mut work, gate);
                if let Err(ref refusal) = low {
                    terminal(refusal)?;
                }
                for i in 0..12 {
                    let derivative = match (&high, &low) {
                        (Ok(hi), Ok(lo)) => {
                            (hi.residual[i] - lo.residual[i]) / (high_x[j] - low_x[j])
                        }
                        (Ok(hi), Err(_)) => {
                            (hi.residual[i] - value.residual[i]) / (high_x[j] - x[j])
                        }
                        (Err(_), Ok(lo)) => {
                            (value.residual[i] - lo.residual[i]) / (x[j] - low_x[j])
                        }
                        (Err(error), Err(_)) => return Err(error.clone()),
                    };
                    matrix[12 * i + j] = finite(derivative * scale[j] / scale[i])?;
                }
            }
            for i in 0..12 {
                rhs[i] = -value.residual[i] / scale[i];
            }
            self.linear_solve
                .solve_into(&matrix, &rhs, &mut correction)?;
            let mut fraction = 1.0;
            let mut accepted = false;
            for _ in 0..self.config.maximum_backtracks {
                let trial_x = core::array::from_fn(|i| x[i] + fraction * scale[i] * correction[i]);
                match self.evaluate(trial_x, loads, &mut work, gate) {
                    Ok(trial) if self.scaled_residual(trial_x, trial)? < error => {
                        x = trial_x;
                        value = trial;
                        accepted = true;
                        break;
                    }
                    Err(refusal) => {
                        terminal(&refusal)?;
                    }
                    _ => {}
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(RigidCohesiveError::NotConverged {
                    scaled_residual: error,
                    iterations: iteration + 1,
                });
            }
        }
        unreachable!("iteration limit returns before the loop ends")
    }

    fn scaled_residual(&self, x: [f64; 12], value: Evaluation) -> Result<f64, RigidCohesiveError> {
        let mut worst: f64 = 0.0;
        for i in 0..12 {
            let absolute = if i % 6 < 3 {
                self.config.force_absolute_tolerance_n
            } else {
                self.config.torque_absolute_tolerance_n_m
            };
            let tolerance = finite(
                absolute + self.config.relative_tolerance * x[i].abs().max(value.expected[i].abs()),
            )?;
            worst = worst.max(finite(value.residual[i].abs() / tolerance)?);
        }
        Ok(worst)
    }

    fn evaluate(
        &mut self,
        x: [f64; 12],
        loads: [RigidCohesiveLoad; 2],
        work: &mut Work,
        gate: Option<&CancelGate>,
    ) -> Result<Evaluation, RigidCohesiveError> {
        poll(gate)?;
        if work.evaluations >= self.config.maximum_evaluations {
            return Err(RigidCohesiveError::Budget(
                "complete rigid patch evaluations",
            ));
        }
        work.evaluations += 1;
        for value in x {
            finite(value)?;
        }
        let wrenches = unpack(x);
        let receipts = [
            self.integrator.step_cayley(
                self.bodies[0].state,
                self.bodies[0].properties,
                wrenches[0],
                self.dt_s,
            )?,
            self.integrator.step_cayley(
                self.bodies[1].state,
                self.bodies[1].properties,
                wrenches[1],
                self.dt_s,
            )?,
        ];
        let before = self.bodies.each_ref().map(|body| body.state);
        let after = receipts.map(|receipt| receipt.state_after);
        let mut cayley = [Vec3::ZERO; 2];
        for body in 0..2 {
            cayley[body] = cayley_increment(
                before[body].pose().orientation(),
                after[body].pose().orientation(),
                self.config.objective.maximum_cayley_increment,
            )?;
        }
        let mut contact = [RigidCohesiveLoad::ZERO; 2];
        let mut storage = 0.0;
        let mut fracture = 0.0;
        let mut viscous = 0.0;
        let mut spatial_work = 0.0;
        for (index, site) in self.sites.iter().enumerate() {
            poll(gate)?;
            if work.material_evaluations >= self.config.maximum_material_evaluations {
                return Err(RigidCohesiveError::Budget(
                    "original cohesive material evaluations",
                ));
            }
            work.material_evaluations += 1;
            let frame = compose(after[0].pose().orientation(), site.frame_in_a)?;
            let trial = site.point.trial(
                frame,
                separation(&after, site.arms_body_m)?,
                self.dt_s,
                self.config.objective,
            )?;
            let numerical = trial.numerical();
            let force = numerical.discrete_wrench.force_on_a_n;
            let couple = numerical.discrete_wrench.frame_couple_n_m;
            for body in 0..2 {
                let arm = before[body]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(site.arms_body_m[body])
                    .scale(0.5)
                    .add(
                        after[body]
                            .pose()
                            .orientation()
                            .rotate_body_to_world(site.arms_body_m[body])
                            .scale(0.5),
                    );
                let action = if body == 0 { force } else { force.scale(-1.0) };
                contact[body].force_world_n = contact[body].force_world_n.add(action);
                contact[body].moment_world_n_m =
                    contact[body].moment_world_n_m.add(arm.cross(action));
            }
            contact[0].moment_world_n_m = contact[0].moment_world_n_m.add(couple);
            let area = numerical.point.area_m2();
            storage += area * numerical.point.material().stored_energy_j_m2();
            fracture += area * numerical.material.equilibrium.energy.dissipated;
            viscous += area * numerical.material.viscous_dissipated_j_m2;
            spatial_work += numerical.spatial_work_j;
            self.staged_sites[index] = Some(trial);
        }
        let mut expected = [Wrench::ZERO; 2];
        for body in 0..2 {
            expected[body].force_world = loads[body].force_world_n.add(contact[body].force_world_n);
            let moment = loads[body]
                .moment_world_n_m
                .add(contact[body].moment_world_n_m);
            // Rbar = R0 (I - [c_body]/2)^-1, so Rbar^-1 M equals
            // R0^T (M - c_world x M / 2). No normalized midpoint rotation.
            expected[body].torque_body = before[body]
                .pose()
                .orientation()
                .rotate_world_to_body(moment.sub(cayley[body].cross(moment).scale(0.5)));
        }
        let expected = pack(expected);
        let residual = core::array::from_fn(|i| x[i] - expected[i]);
        for value in
            expected
                .into_iter()
                .chain(residual)
                .chain([storage, fracture, viscous, spatial_work])
        {
            finite(value)?;
        }
        Ok(Evaluation {
            receipts,
            residual,
            expected,
            contact,
            cayley,
            storage,
            fracture,
            viscous,
            spatial_work,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn audit(
        &self,
        trial: Evaluation,
        wrenches: [Wrench; 2],
        loads: [RigidCohesiveLoad; 2],
        step: u64,
        time_s: f64,
        iterations: usize,
        work: Work,
        gate: Option<&CancelGate>,
    ) -> Result<RigidCohesiveFrame, RigidCohesiveError> {
        let mut mechanical = 0.0;
        let mut external_work = 0.0;
        let mut old_p = Vec3::ZERO;
        let mut new_p = Vec3::ZERO;
        let mut old_l = Vec3::ZERO;
        let mut new_l = Vec3::ZERO;
        let mut external_impulse = Vec3::ZERO;
        let mut external_angular_impulse = Vec3::ZERO;
        let mut force_error: f64 = 0.0;
        let mut torque_error: f64 = 0.0;
        for body in 0..2 {
            poll(gate)?;
            let receipt = trial.receipts[body];
            let old = receipt.state_before;
            let new = receipt.state_after;
            self.bodies[body].check(new)?;
            bound(
                "total body force",
                norm(wrenches[body].force_world),
                self.config.maximum_total_force_n,
            )?;
            bound(
                "total body torque",
                norm(wrenches[body].torque_body),
                self.config.maximum_total_torque_n_m,
            )?;
            let lmid = old
                .angular_momentum_body()
                .scale(0.5)
                .add(new.angular_momentum_body().scale(0.5));
            let native_cayley = old.pose().orientation().rotate_body_to_world(
                self.bodies[body]
                    .properties
                    .angular_velocity_body_checked(lmid)?
                    .scale(self.dt_s),
            );
            balance(
                "original Cayley kinematics",
                norm(native_cayley.sub(trial.cayley[body])),
                self.config.frame_tolerance,
            )?;
            let dx = new.pose().position_world().sub(old.pose().position_world());
            let xmid = old
                .pose()
                .position_world()
                .scale(0.5)
                .add(new.pose().position_world().scale(0.5));
            external_work += loads[body].force_world_n.dot(dx)
                + loads[body].moment_world_n_m.dot(trial.cayley[body]);
            mechanical += receipt.diagnostics_after.mechanical_energy;
            old_p = old_p.add(receipt.diagnostics_before.linear_momentum_world);
            new_p = new_p.add(receipt.diagnostics_after.linear_momentum_world);
            old_l = old_l.add(receipt.diagnostics_before.angular_momentum_world);
            new_l = new_l.add(receipt.diagnostics_after.angular_momentum_world);
            let total_external = loads[body].force_world_n.add(
                self.integrator
                    .gravity()
                    .acceleration_world()
                    .scale(self.bodies[body].properties.mass()),
            );
            external_impulse = external_impulse.add(total_external.scale(self.dt_s));
            external_angular_impulse = external_angular_impulse.add(
                xmid.cross(total_external)
                    .add(loads[body].moment_world_n_m)
                    .scale(self.dt_s),
            );
            for axis in 0..3 {
                force_error = force_error.max(trial.residual[6 * body + axis].abs());
                torque_error = torque_error.max(trial.residual[6 * body + 3 + axis].abs());
            }
        }
        for (site, staged) in self.sites.iter().zip(&self.staged_sites) {
            poll(gate)?;
            let value = staged
                .as_ref()
                .ok_or(RigidCohesiveError::Invalid("incomplete staged patch"))?
                .numerical();
            bound(
                "site mean force",
                norm(value.discrete_wrench.force_on_a_n),
                self.config.maximum_site_force_n,
            )?;
            let previous = site.point.numerical().material().jump_m();
            let next = value.point.material().jump_m();
            for axis in 0..3 {
                bound(
                    "site jump increment",
                    (next[axis] - previous[axis]).abs(),
                    self.config.maximum_jump_increment_m,
                )?;
            }
        }
        finite(mechanical)?;
        finite(external_work)?;
        let old_energy = self.total_energy_j()?;
        let energy_residual = finite(
            mechanical + trial.storage - old_energy + trial.fracture + trial.viscous
                - external_work,
        )?;
        let energy_tolerance = finite(
            self.config.energy_absolute_tolerance_j
                + self.config.energy_relative_tolerance
                    * (old_energy.abs()
                        + mechanical.abs()
                        + trial.storage.abs()
                        + external_work.abs()
                        + trial.fracture
                        + trial.viscous),
        )?;
        balance(
            "actual rigid/cohesive energy",
            energy_residual,
            energy_tolerance,
        )?;
        let p_residual = new_p.sub(old_p).sub(external_impulse);
        let l_residual = new_l.sub(old_l).sub(external_angular_impulse);
        let p_tolerance = finite(
            self.config.linear_momentum_absolute_tolerance_n_s
                + self.config.momentum_relative_tolerance
                    * (norm(old_p) + norm(new_p) + norm(external_impulse)),
        )?;
        let l_tolerance = finite(
            self.config.angular_momentum_absolute_tolerance_n_m_s
                + self.config.momentum_relative_tolerance
                    * (norm(old_l) + norm(new_l) + norm(external_angular_impulse)),
        )?;
        balance(
            "actual total linear momentum",
            norm(p_residual),
            p_tolerance,
        )?;
        balance(
            "actual total angular momentum",
            norm(l_residual),
            l_tolerance,
        )?;
        poll(gate)?;
        Ok(RigidCohesiveFrame {
            step,
            time_s,
            mechanical_energy_j: mechanical,
            cohesive_stored_energy_j: trial.storage,
            external_work_j: external_work,
            fracture_dissipation_j: trial.fracture,
            viscous_dissipation_j: trial.viscous,
            cohesive_spatial_work_j: trial.spatial_work,
            energy_residual_j: energy_residual,
            energy_tolerance_j: energy_tolerance,
            linear_momentum_residual_n_s: p_residual,
            linear_momentum_tolerance_n_s: p_tolerance,
            angular_momentum_residual_n_m_s: l_residual,
            angular_momentum_tolerance_n_m_s: l_tolerance,
            contact_loads: trial.contact,
            native_wrenches: wrenches,
            maximum_force_residual_n: force_error,
            maximum_torque_residual_n_m: torque_error,
            iterations,
            evaluations: work.evaluations,
            material_evaluations: work.material_evaluations,
        })
    }
}

fn pack(wrenches: [Wrench; 2]) -> [f64; 12] {
    let mut packed = [0.0; 12];
    for (body, wrench) in wrenches.into_iter().enumerate() {
        packed[6 * body..6 * body + 6].copy_from_slice(&[
            wrench.force_world.x,
            wrench.force_world.y,
            wrench.force_world.z,
            wrench.torque_body.x,
            wrench.torque_body.y,
            wrench.torque_body.z,
        ]);
    }
    packed
}
fn unpack(x: [f64; 12]) -> [Wrench; 2] {
    core::array::from_fn(|body| Wrench {
        force_world: Vec3::new(x[6 * body], x[6 * body + 1], x[6 * body + 2]),
        torque_body: Vec3::new(x[6 * body + 3], x[6 * body + 4], x[6 * body + 5]),
    })
}
fn terminal(error: &RigidCohesiveError) -> Result<(), RigidCohesiveError> {
    match error {
        RigidCohesiveError::Cancelled | RigidCohesiveError::Budget(_) => Err(error.clone()),
        _ => Ok(()),
    }
}

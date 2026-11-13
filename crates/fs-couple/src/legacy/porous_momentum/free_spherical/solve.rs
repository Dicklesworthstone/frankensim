//! Six-component residuals evaluated through the original native owners.

use super::*;
use fs_lbm::d3q19::LinearResistance3;
use fs_mbd::{Gravity, RigidBodyIntegrator, Wrench};

impl Step<'_> {
    pub(super) fn solve(self, gate: &CancelGate) -> Result<FreePorousStep, PorousMomentumError> {
        let c = self.controls.flow;
        let scratch = RefCell::new(self.motion.fluid.clone());
        let trials = Cell::new(0usize);
        let trial_cells = Cell::new(0usize);
        let evaluate = |x: &[f64]| -> Result<Trial, PorousMomentumError> {
            poll(gate)?;
            let count = trials.get() + 1;
            let work = trial_cells
                .get()
                .checked_add(self.motion.cells)
                .ok_or(invalid("free porous trial work overflow"))?;
            if count > c.maximum_native_trials || work > c.maximum_trial_cells {
                return Err(PorousMomentumError::Budget {
                    trials: count,
                    trial_cells: work,
                });
            }
            trials.set(count);
            trial_cells.set(work);
            self.evaluate(x, &mut scratch.borrow_mut(), gate)
        };
        let problem = FallibleFiniteDifference::new(6, c.difference_step, |x: &[f64]| {
            evaluate(x).map(|t| t.residual)
        })
        .map_err(difference_error)?;
        // A momentum-converged trial must also resolve its interface work.
        // The speed bound converts the two scaled impulse defects to joules;
        // keep most of the physical budget for independent native defects.
        let work_scale = finite(
            c.maximum_speed_m_s * self.impulse_scales[0]
                + (c.maximum_speed_m_s / self.motion.radius_m) * self.impulse_scales[1],
        )?;
        let tolerance = 0.25
            * (c.linear_tolerance_kg_m_s / self.impulse_scales[0])
                .min(c.angular_tolerance_kg_m2_s / self.impulse_scales[1])
                .min(c.work_tolerance_j / work_scale);
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(invalid("free porous solve tolerance is not representable"));
        }
        let config = NewtonKrylovConfig {
            absolute_tolerance: tolerance,
            relative_tolerance: 0.0,
            linear_restart: 6,
            max_linear_cycles: 4,
            forcing_minimum: 1e-3,
            forcing_maximum: 0.1,
            ..Default::default()
        };
        let construction = NewtonKrylovState::new(&problem, vec![0.0; 6], config);
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = construction.map_err(PorousMomentumError::Solver)?;
        let result = solver.run_controlled(&problem, c.maximum_iterations, &mut || {
            problem.has_failed() || gate.is_requested()
        });
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        poll(gate)?;
        if result.progress != SolverRunProgress::Complete {
            return Err(PorousMomentumError::SolverProgress(result.progress));
        }
        if !result.report.converged {
            return Err(PorousMomentumError::Unresolved(result.report));
        }
        // This replay is authoritative: its fluid remains in the same scratch
        // owner that is returned only after the independent physical checks.
        let trial = evaluate(&solver.x)?;
        let mut nonlinear = result.report;
        nonlinear.residual_norm = fs_solver::norm2(&trial.residual);
        if !nonlinear.residual_norm.is_finite() || nonlinear.residual_norm > tolerance {
            nonlinear.converged = false;
            return Err(PorousMomentumError::Unresolved(nonlinear));
        }
        self.finish(
            scratch.into_inner(),
            trial,
            nonlinear,
            trials.get(),
            trial_cells.get(),
            gate,
        )
    }

    fn evaluate(
        &self,
        x: &[f64],
        fluid: &mut BoundaryGrid3,
        gate: &CancelGate,
    ) -> Result<Trial, PorousMomentumError> {
        let m = self.motion;
        let c = self.controls.flow;
        let dt = m.units.timestep_s;
        let load = |offset: usize, scale: f64| {
            v(std::array::from_fn(|i| {
                (dd(x[offset + i]) * dd(scale) / dd(dt)).to_f64()
            }))
        };
        let mechanics = RigidBodyIntegrator::new(Gravity::ZERO).step_cayley(
            m.body,
            m.properties,
            Wrench {
                force_world: load(0, self.impulse_scales[0]),
                torque_body: load(3, self.impulse_scales[1]),
            },
            dt,
        )?;
        let body = mechanics.state_after;
        contain_sphere(
            body.pose().position_world(),
            m.radius_m,
            m.units,
            m.fluid.dimensions(),
        )?;
        admit_speed(
            body.center_of_mass_velocity_world(m.properties)?,
            angular_velocity(body, m.properties)?,
            m.radius_m,
            c.maximum_speed_m_s,
        )?;
        let centre = average(m.body.pose().position_world(), body.pose().position_world());
        let velocity = average(
            m.body.center_of_mass_velocity_world(m.properties)?,
            body.center_of_mass_velocity_world(m.properties)?,
        );
        let omega_body = m.properties.angular_velocity_body_checked(average(
            m.body.angular_momentum_body(),
            body.angular_momentum_body(),
        ))?;
        let omega_world = m.body.pose().orientation().rotate_body_to_world(omega_body);
        admit_speed(velocity, omega_world, m.radius_m, c.maximum_speed_m_s)?;
        if norm(omega_world) * dt > self.controls.maximum_rotation_rad {
            return Err(invalid(
                "free porous trial rotation exceeds the declared domain",
            ));
        }
        let linear_change =
            difference(body.linear_momentum_world(), m.body.linear_momentum_world());
        let spin_change = difference(world_spin(body), world_spin(m.body));
        let lattice_origin = v(std::array::from_fn(|i| {
            ((dd(a(centre)[i]) - dd(a(m.units.grid_origin_m)[i])) / dd(m.units.cell_size_m))
                .to_f64()
        }));
        if !lattice_origin.is_finite() {
            return Err(invalid("free porous midpoint origin is not representable"));
        }
        // One prepared material frame per native trial. The Cayley endpoint
        // rotation is shorter than pi; its geodesic half carries actual fixed
        // body geometry and anisotropy without changing the work velocity.
        let frame = midpoint_orientation(m.body.pose().orientation(), body.pose().orientation())?;
        let support = PreparedSupport::new(m.support, frame);
        let axes = support.resistance_axes(self.resistance);
        let [nx, ny, _] = m.fluid.dimensions();
        fluid.clone_from(&m.fluid);
        let exchange = fluid.step_with_darcy_forchheimer_and_wall_momentum(
            |cell @ [x, y, z]| {
                let arm = v(std::array::from_fn(|i| {
                    (dd(a(m.units.grid_origin_m)[i]) - dd(a(centre)[i])
                        + (dd(cell[i] as f64) + dd(0.5)) * dd(m.units.cell_size_m))
                    .to_f64()
                }));
                let inside = support.contains(arm);
                let law = if inside {
                    self.resistance
                } else {
                    PrincipalPorousResistance::default()
                };
                let skeleton_velocity = if inside {
                    a(velocity.add(omega_world.cross(arm))).map(|u| u / self.units.velocity)
                } else {
                    [0.0; 3]
                };
                DarcyForchheimer3 {
                    linear: LinearResistance3 {
                        principal: law.darcy_principal,
                        axes,
                        skeleton_velocity,
                        body_force: self.forcing[(z * ny + y) * nx + x],
                    },
                    quadratic: law.forchheimer,
                }
            },
            &m.links,
            a(lattice_origin),
            gate,
        )?;
        let impulse = v(exchange.skeleton_impulse).scale(self.units.impulse);
        let moment = v(exchange.skeleton_angular_impulse).scale(self.units.angular);
        let mut residual = Vec::with_capacity(6);
        for (change, drag, external, scale) in [
            (
                linear_change,
                impulse,
                self.external_linear,
                self.impulse_scales[0],
            ),
            (
                spin_change,
                moment,
                self.external_couple,
                self.impulse_scales[1],
            ),
        ] {
            for i in 0..3 {
                residual.push(finite(
                    ((dd(a(change)[i]) - dd(a(drag)[i]) - dd(a(external)[i])) / dd(scale)).to_f64(),
                )?);
            }
        }
        Ok(Trial {
            mechanics,
            centre,
            velocity,
            omega_body,
            omega_world,
            linear_change,
            spin_change,
            exchange,
            residual,
        })
    }
}

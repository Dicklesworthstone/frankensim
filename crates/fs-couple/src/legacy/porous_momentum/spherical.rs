//! Full rotational recoil of finite spherical porous support.
//!
//! Only the COM is constrained. All three components of the measured fluid
//! couple act on the ORIGINAL rigid inertia, including its gyroscopic motion.
//! Radial isotropic coefficients and spherical support stay invariant under
//! every attitude; this does not stand in for a rotating general solid mask.
//!
//! The original rigid Cayley step supplies the exact discrete spatial impulse
//! and work map. Its conjugate angular velocity is `R0 * omega_body_mid`, not
//! the average of endpoint world angular velocities or a normalized midpoint
//! quaternion. Actual nonlinear fluid trials share that same velocity field.

use super::rotor::{PorousRotorControls, PorousRotorResistance};
use super::*;
use fs_lbm::d3q19::{BoundaryLink3, LinearResistance3};
use fs_mbd::{Gravity, RigidBodyIntegrator, StepReceipt, Wrench};

/// Original fluid and finite dry inertia, free to rotate about a fixed COM.
/// Radius bounds the porous support; native mass/inertia remain authored data,
/// not inferred from lattice populations or replaced by apparent fluid inertia.
pub struct SphericalPorousBody {
    fluid: BoundaryGrid3,
    body: RigidBodyState,
    properties: MassProperties,
    units: PorousLatticeUnits3,
    radius_m: f64,
    links: Vec<BoundaryLink3>,
    cells: usize,
    time_s: f64,
    accepted_steps: u64,
}

/// Accepted original owners and actual finite-spin fluid/body exchange.
pub struct SphericalPorousStep {
    pub motion: SphericalPorousBody,
    /// Native lattice-unit exchange, with moments about `exchange_origin_m`.
    pub exchange: ResistanceMomentumExchange3,
    pub exchange_origin_m: Vec3,
    pub mechanics: StepReceipt,
    pub nonlinear: NewtonReport,
    pub native_trials: usize,
    pub trial_cells: usize,
    /// Native Cayley work-conjugate velocity `R0 * I^-1 (L0+L1)/2`.
    pub angular_velocity_world_rad_s: Vec3,
    /// Bearing force impulse ON the body at its COM. It supplies no couple
    /// and does no work; ALL of the distributed drag moment reaches the body.
    pub guide_impulse_kg_m_s: Vec3,
    pub balance_origin_m: Vec3,
    pub fluid_linear_change_kg_m_s: Vec3,
    pub fluid_angular_change_kg_m2_s: Vec3,
    pub body_angular_change_kg_m2_s: Vec3,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    /// Participant norms are summed before acceptance, preventing cancellation.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub fluid_mass_change_kg: f64,
    pub fluid_kinetic_change_j: f64,
    pub body_kinetic_change_j: f64,
    pub body_external_work_j: f64,
    pub fluid_external_work_j: f64,
    pub fluid_drag_work_j: f64,
    pub relative_dissipation_j: f64,
    /// Native fluid transport/pressure/collision remainder, not physical heat
    /// or a complete-fluid energy error certificate.
    pub other_fluid_kinetic_change_j: f64,
    pub mechanical_work_residual_j: f64,
}

struct Sample {
    law: PorousRotorResistance,
    force: [f64; 3],
    /// Physical lever arm divided by lattice velocity scale, zero off support.
    arm_per_velocity: Vec3,
}
struct Trial {
    mechanics: StepReceipt,
    omega_body: Vec3,
    omega_world: Vec3,
    spin_change: Vec3,
    exchange: ResistanceMomentumExchange3,
    residual: Vec<f64>,
}

impl SphericalPorousBody {
    /// Bind bounded original populations and dry native inertia. The whole
    /// sphere must fit inside six stationary no-slip walls, without voxel
    /// solids. Initial linear momentum is zero; initial spin and principal
    /// frame may be arbitrary. The bearing constrains translation only.
    pub fn new(
        fluid: BoundaryGrid3,
        body: RigidBodyState,
        properties: MassProperties,
        units: PorousLatticeUnits3,
        radius_m: f64,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PorousMomentumError> {
        poll(gate)?;
        units.admit()?;
        body.kinetic_energy(properties)?;
        if body.linear_momentum_world() != Vec3::ZERO || !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(invalid(
                "spherical porous support requires a fixed COM and positive finite radius",
            ));
        }
        let centre = body.pose().position_world();
        let dimensions = fluid.dimensions();
        for i in 0..3 {
            let lo = a(units.grid_origin_m)[i];
            let hi = finite(lo + dimensions[i] as f64 * units.cell_size_m)?;
            if a(centre)[i] - radius_m < lo || a(centre)[i] + radius_m > hi {
                return Err(invalid(
                    "complete spherical support must fit inside the box",
                ));
            }
        }
        let (cells, links) = closed_domain(&fluid, maximum_cells, gate)?;
        Ok(Self {
            fluid,
            body,
            properties,
            units,
            radius_m,
            links,
            cells,
            time_s: 0.0,
            accepted_steps: 0,
        })
    }
    pub const fn fluid(&self) -> &BoundaryGrid3 {
        &self.fluid
    }
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    pub const fn units(&self) -> PorousLatticeUnits3 {
        self.units
    }
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Actual endpoint angular velocity, resolved in world coordinates.
    pub fn angular_velocity_world_rad_s(&self) -> Result<Vec3, PorousMomentumError> {
        Ok(self.body.pose().orientation().rotate_body_to_world(
            self.properties
                .angular_velocity_body_checked(self.body.angular_momentum_body())?,
        ))
    }

    /// Advance one original lattice tick with reciprocal three-axis rotation.
    /// `resistance(radius_m)` supplies isotropic native lattice coefficients
    /// inside the sphere. Outside, drag and imposed skeleton velocity are zero.
    /// The independent `force(cell)` is a native lattice force density over the
    /// complete box. Sample each callback once per eligible cell and hold it at
    /// this tick's midpoint. `external_couple_nm` is a world-resolved couple.
    ///
    /// Each three-coordinate Newton trial uses the ORIGINAL Euler-top/Cayley
    /// step and original nonlinear fluid collision/stream/endpoint-force step.
    /// With `Rbar=(R0+R1)/2`, the native spatial impulse is `h Rbar tau_body`,
    /// paired exactly with `Omega=R0 omega_body_mid` in work. Solving actual
    /// endpoint spatial momentum avoids assuming a fixed axis or adding an
    /// unphysical transverse guide couple. Native final replay is authoritative.
    ///
    /// `controls.flow.velocity_scale_m_s` sets a rim-speed impulse scale through
    /// the largest original principal inertia. Its midpoint tolerance bounds
    /// the native Cayley velocity identity. `maximum_rotation_rad` bounds
    /// `h |Omega|` (actual Cayley angle is `2 atan(h |Omega|/2)`). Independent
    /// time/space refinement remains necessary; source-work closure is not a
    /// fluid-energy or attitude-accuracy theorem. No candidate publishes on error.
    pub fn step_candidate(
        &self,
        mut resistance: impl FnMut(f64) -> PorousRotorResistance,
        mut force: impl FnMut([usize; 3]) -> [f64; 3],
        external_couple_nm: Vec3,
        controls: PorousRotorControls,
        gate: &CancelGate,
    ) -> Result<SphericalPorousStep, PorousMomentumError> {
        poll(gate)?;
        let c = controls.flow;
        let units = self.units.admit()?;
        c.admit(units, self.cells)?;
        if !external_couple_nm.is_finite()
            || !controls.maximum_rotation_rad.is_finite()
            || controls.maximum_rotation_rad <= 0.0
            || norm(self.angular_velocity_world_rad_s()?) * self.radius_m > c.maximum_speed_m_s
        {
            return Err(invalid(
                "finite spatial couple, bounded native rotation and rim speed required",
            ));
        }
        let dt = self.units.timestep_s;
        let external = external_couple_nm.scale(dt);
        let next_time = finite(self.time_s + dt)?;
        let next_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or(invalid("step counter overflow"))?;
        if next_time <= self.time_s || !external.is_finite() {
            return Err(invalid(
                "physical clock or external angular impulse is not representable",
            ));
        }
        let largest_inertia = a(self.properties.principal_inertia_body())
            .into_iter()
            .fold(0.0_f64, f64::max);
        let impulse_scale = finite(2.0 * largest_inertia * (c.velocity_scale_m_s / self.radius_m))?;
        if impulse_scale <= 0.0 {
            return Err(invalid("rotational impulse scale is not representable"));
        }
        let origin = self.body.pose().position_world();
        let lattice_origin = origin
            .sub(self.units.grid_origin_m)
            .scale(self.units.cell_size_m.recip());
        if !lattice_origin.is_finite() {
            return Err(invalid("spherical moment origin is not representable"));
        }
        let _ = fluid_change(&self.fluid, &self.fluid, self.units, c, gate)?;
        let [nx, ny, nz] = self.fluid.dimensions();
        let mut field = Vec::with_capacity(self.cells);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x % 256 == 0 {
                        poll(gate)?;
                    }
                    let cell = [x, y, z];
                    let arm = v(std::array::from_fn(|i| {
                        (dd(a(self.units.grid_origin_m)[i]) - dd(a(origin)[i])
                            + (dd(cell[i] as f64) + dd(0.5)) * dd(self.units.cell_size_m))
                        .to_f64()
                    }));
                    let radial = finite(norm(arm))?;
                    let inside = radial <= self.radius_m;
                    let law = if inside {
                        resistance(radial)
                    } else {
                        PorousRotorResistance::default()
                    };
                    let force = force(cell);
                    if [law.darcy, law.forchheimer]
                        .iter()
                        .any(|x| !x.is_finite() || *x < 0.0)
                        || force.iter().any(|x| !x.is_finite())
                    {
                        return Err(invalid(
                            "radial native resistance must be finite nonnegative; force finite",
                        ));
                    }
                    let arm_per_velocity = if inside {
                        arm.scale(units.velocity.recip())
                    } else {
                        Vec3::ZERO
                    };
                    if !arm_per_velocity.is_finite() {
                        return Err(invalid("spherical lever arm is not representable"));
                    }
                    field.push(Sample {
                        law,
                        force,
                        arm_per_velocity,
                    });
                }
            }
        }
        let scratch = RefCell::new(self.fluid.clone());
        let trials = Cell::new(0usize);
        let trial_cells = Cell::new(0usize);
        let evaluate = |x: &[f64]| -> Result<Trial, PorousMomentumError> {
            poll(gate)?;
            let count = trials.get() + 1;
            let work = trial_cells
                .get()
                .checked_add(self.cells)
                .ok_or(invalid("trial work overflow"))?;
            if count > c.maximum_native_trials || work > c.maximum_trial_cells {
                return Err(PorousMomentumError::Budget {
                    trials: count,
                    trial_cells: work,
                });
            }
            trials.set(count);
            trial_cells.set(work);
            let torque_body = v(std::array::from_fn(|i| {
                (dd(x[i]) * dd(impulse_scale) / dd(dt)).to_f64()
            }));
            let mechanics = RigidBodyIntegrator::new(Gravity::ZERO).step_cayley(
                self.body,
                self.properties,
                Wrench {
                    force_world: Vec3::ZERO,
                    torque_body,
                },
                dt,
            )?;
            let body = mechanics.state_after;
            let midpoint_spin = self
                .body
                .angular_momentum_body()
                .scale(0.5)
                .add(body.angular_momentum_body().scale(0.5));
            let omega_body = self
                .properties
                .angular_velocity_body_checked(midpoint_spin)?;
            let omega_world = self
                .body
                .pose()
                .orientation()
                .rotate_body_to_world(omega_body);
            if norm(omega_world) * self.radius_m > c.maximum_speed_m_s
                || norm(omega_world) * dt > controls.maximum_rotation_rad
            {
                return Err(invalid(
                    "trial spherical rim speed or rotation exceeds the declared domain",
                ));
            }
            let spin_change = spin_change(self.body, body);
            let mut fluid = scratch.borrow_mut();
            fluid.clone_from(&self.fluid);
            let exchange = fluid.step_with_darcy_forchheimer_and_wall_momentum(
                |[x, y, z]| {
                    let sample = &field[(z * ny + y) * nx + x];
                    DarcyForchheimer3 {
                        linear: LinearResistance3 {
                            principal: [sample.law.darcy; 3],
                            skeleton_velocity: a(omega_world.cross(sample.arm_per_velocity)),
                            body_force: sample.force,
                            ..Default::default()
                        },
                        quadratic: sample.law.forchheimer,
                    }
                },
                &self.links,
                a(lattice_origin),
                gate,
            )?;
            let moment = v(exchange.skeleton_angular_impulse).scale(units.angular);
            let residual = (0..3)
                .map(|i| {
                    finite(
                        ((dd(a(spin_change)[i]) - dd(a(moment)[i]) - dd(a(external)[i]))
                            / dd(impulse_scale))
                        .to_f64(),
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(Trial {
                mechanics,
                omega_body,
                omega_world,
                spin_change,
                exchange,
                residual,
            })
        };
        let problem = FallibleFiniteDifference::new(3, c.difference_step, |x: &[f64]| {
            evaluate(x).map(|t| t.residual)
        })
        .map_err(difference_error)?;
        let tolerance = 0.25 * c.angular_tolerance_kg_m2_s / impulse_scale;
        let config = NewtonKrylovConfig {
            absolute_tolerance: tolerance,
            relative_tolerance: 0.0,
            linear_restart: 3,
            max_linear_cycles: 4,
            forcing_minimum: 1e-3,
            forcing_maximum: 0.1,
            ..Default::default()
        };
        let construction = NewtonKrylovState::new(&problem, vec![0.0; 3], config);
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
            external,
            controls,
            units,
            nonlinear,
            trials.get(),
            trial_cells.get(),
            next_time,
            next_steps,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        fluid: BoundaryGrid3,
        trial: Trial,
        external: Vec3,
        controls: PorousRotorControls,
        units: Units,
        nonlinear: NewtonReport,
        native_trials: usize,
        trial_cells: usize,
        time_s: f64,
        accepted_steps: u64,
        gate: &CancelGate,
    ) -> Result<SphericalPorousStep, PorousMomentumError> {
        let c = controls.flow;
        let body = trial.mechanics.state_after;
        let origin = self.body.pose().position_world();
        if body.pose().position_world() != origin || body.linear_momentum_world() != Vec3::ZERO {
            return Err(invalid("native spherical motion changed its fixed COM"));
        }
        let omega_after = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        if norm(omega_after) * self.radius_m > c.maximum_speed_m_s {
            return Err(invalid(
                "accepted spherical rim speed exceeds the declared domain",
            ));
        }
        check(
            "native Cayley work-velocity identity",
            norm(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(trial.omega_body)
                    .sub(trial.omega_world),
            ) * self.radius_m,
            c.midpoint_tolerance_m_s,
        )?;
        let exchange = trial.exchange;
        let impulse = v(exchange.skeleton_impulse).scale(units.impulse);
        let moment = v(exchange.skeleton_angular_impulse).scale(units.angular);
        let guide = impulse.scale(-1.0);
        let arm = origin.sub(c.balance_origin_m);
        let shift = |moment: Vec3, impulse: Vec3| {
            moment.add(add_cross(arm, impulse, Vec3::ZERO, Vec3::ZERO))
        };
        let drag_angular = shift(moment, impulse);
        let guide_angular = shift(Vec3::ZERO, guide);
        let body_angular = v(std::array::from_fn(|i| {
            (dd(a(trial.spin_change)[i])
                - dd(a(drag_angular)[i])
                - dd(a(guide_angular)[i])
                - dd(a(external)[i]))
            .to_f64()
        }));
        let body_linear = impulse.add(guide).scale(-1.0);
        let change = fluid_change(&self.fluid, &fluid, self.units, c, gate)?;
        let fluid_external = v(exchange.external_impulse).scale(units.impulse);
        let fluid_external_angular = shift(
            v(exchange.external_angular_impulse).scale(units.angular),
            fluid_external,
        );
        let wall = v(exchange.walls.fluid_population_impulse).scale(units.impulse);
        let wall_angular = shift(
            v(exchange.walls.fluid_population_angular_impulse).scale(units.angular),
            wall,
        );
        let fluid_linear = v(std::array::from_fn(|i| {
            (dd(a(change.linear)[i]) - dd(a(fluid_external)[i]) + dd(a(impulse)[i])
                - dd(a(wall)[i]))
            .to_f64()
        }));
        let fluid_angular = v(std::array::from_fn(|i| {
            (dd(a(change.angular)[i]) - dd(a(fluid_external_angular)[i]) + dd(a(drag_angular)[i])
                - dd(a(wall_angular)[i]))
            .to_f64()
        }));
        let absolute_linear = norm(body_linear) + norm(fluid_linear);
        let absolute_angular = norm(body_angular) + norm(fluid_angular);
        check(
            "spherical fluid/body linear momentum including bearing and walls",
            absolute_linear,
            c.linear_tolerance_kg_m_s,
        )?;
        check(
            "spherical fluid/body angular momentum without a bearing couple",
            absolute_angular,
            c.angular_tolerance_kg_m2_s,
        )?;
        let kinetic = (0..3)
            .fold(Dd::ZERO, |sum, i| {
                sum + dd(a(trial.omega_body)[i])
                    * (dd(a(body.angular_momentum_body())[i])
                        - dd(a(self.body.angular_momentum_body())[i]))
            })
            .to_f64();
        let fluid_work = finite(exchange.external_work * units.energy)?;
        let drag_work = finite(exchange.fluid_drag_work * units.energy)?;
        let skeleton_work = finite(exchange.skeleton_work * units.energy)?;
        let loss = finite(exchange.relative_dissipation * units.energy)?;
        if loss < 0.0 {
            return Err(invalid("native porous dissipation is negative"));
        }
        let external_work = dot(trial.omega_world, external);
        let body_drag_work = dot(trial.omega_world, moment);
        let constitutive = (dd(drag_work) + dd(skeleton_work) + dd(loss)).to_f64();
        let interface = (dd(body_drag_work) - dd(skeleton_work)).to_f64();
        let rigid_work = (dd(kinetic) - dd(body_drag_work) - dd(external_work)).to_f64();
        check(
            "noncancelling spherical constitutive/interface/native work",
            constitutive.abs() + interface.abs() + rigid_work.abs(),
            c.work_tolerance_j,
        )?;
        let mechanical = (dd(kinetic) + dd(drag_work) + dd(loss) - dd(external_work)).to_f64();
        let remainder = (dd(change.kinetic) - dd(fluid_work) - dd(drag_work)).to_f64();
        let linear_residual = body_linear.add(fluid_linear);
        let angular_residual = body_angular.add(fluid_angular);
        if [
            impulse,
            moment,
            guide,
            trial.spin_change,
            linear_residual,
            angular_residual,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || [kinetic, external_work, mechanical, remainder]
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "combined spherical momentum/work account is not representable",
            ));
        }
        poll(gate)?;
        Ok(SphericalPorousStep {
            motion: Self {
                fluid,
                body,
                properties: self.properties,
                units: self.units,
                radius_m: self.radius_m,
                links: self.links.clone(),
                cells: self.cells,
                time_s,
                accepted_steps,
            },
            exchange,
            exchange_origin_m: origin,
            mechanics: trial.mechanics,
            nonlinear,
            native_trials,
            trial_cells,
            angular_velocity_world_rad_s: trial.omega_world,
            guide_impulse_kg_m_s: guide,
            balance_origin_m: c.balance_origin_m,
            fluid_linear_change_kg_m_s: change.linear,
            fluid_angular_change_kg_m2_s: change.angular,
            body_angular_change_kg_m2_s: trial.spin_change,
            linear_residual_kg_m_s: linear_residual,
            angular_residual_kg_m2_s: angular_residual,
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            fluid_mass_change_kg: change.mass,
            fluid_kinetic_change_j: change.kinetic,
            body_kinetic_change_j: kinetic,
            body_external_work_j: external_work,
            fluid_external_work_j: fluid_work,
            fluid_drag_work_j: drag_work,
            relative_dissipation_j: loss,
            other_fluid_kinetic_change_j: remainder,
            mechanical_work_residual_j: mechanical,
        })
    }
}

fn spin_change(before: RigidBodyState, after: RigidBodyState) -> Vec3 {
    let spin = |body: RigidBodyState| {
        body.pose()
            .orientation()
            .rotate_body_to_world(body.angular_momentum_body())
    };
    let p0 = spin(before);
    let p1 = spin(after);
    v(std::array::from_fn(|i| {
        (dd(a(p1)[i]) - dd(a(p0)[i])).to_f64()
    }))
}

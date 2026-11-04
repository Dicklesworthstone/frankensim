//! Finite axial inertia and actual attitude for an axisymmetric porous rotor.
//!
//! A finite cylinder inside six stationary walls has fixed COM and a fixed
//! native principal axis. Isotropic Darcy/Forchheimer coefficients depend only
//! on cylindrical radius and axial position; their support is invariant under
//! the actual rotation. This does not approximate a rotating nonaxisymmetric
//! mask, a translated body, or a rotating velocity field across periodic cuts.
//! The original lattice law and rigid midpoint integrator own every candidate.

use super::*;
use fs_lbm::d3q19::{BoundaryLink3, FaceBoundary3, LinearResistance3};
use fs_mbd::{Gravity, RigidBodyIntegrator, StepReceipt, Wrench};

/// Physical, rotation-invariant support, centred at the native body COM.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousRotorCylinder {
    /// Signed unit world axis, also a principal axis of the native dry inertia.
    pub axis_world: Vec3,
    pub radius_m: f64,
    pub half_height_m: f64,
}

/// Isotropic coefficients in the ORIGINAL lattice force-density units.
/// `F_drag = -darcy * slip - forchheimer * |slip| * slip`.
/// These are not intrinsic permeabilities or a porosity/added-inertia model.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PorousRotorResistance {
    pub darcy: f64,
    pub forchheimer: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct PorousRotorControls {
    /// Original flow/work/budget controls. Velocity scale and midpoint
    /// tolerance apply to rim speed; maximum speed bounds fluid and rotor rim
    /// at both endpoints and every trial. No new angular-speed unit is implied.
    pub flow: PorousMomentumControls,
    /// Absolute native midpoint angle during this tick; refine phase in time.
    pub maximum_rotation_rad: f64,
}

/// Native fluid and finite dry rotor. The bearing admits axial rotation only.
pub struct AxisymmetricPorousRotor {
    fluid: BoundaryGrid3,
    body: RigidBodyState,
    properties: MassProperties,
    units: PorousLatticeUnits3,
    cylinder: PorousRotorCylinder,
    axis_body: Vec3,
    inertia_kg_m2: f64,
    links: Vec<BoundaryLink3>,
    cells: usize,
    time_s: f64,
    accepted_steps: u64,
}

pub struct PorousRotorStep {
    pub motion: AxisymmetricPorousRotor,
    /// Original lattice receipt, all moments about the fixed body COM.
    pub exchange: ResistanceMomentumExchange3,
    pub exchange_origin_m: Vec3,
    /// Actual native attitude/momentum step; no apparent inertia replaces it.
    pub mechanics: StepReceipt,
    pub nonlinear: NewtonReport,
    pub native_trials: usize,
    pub trial_cells: usize,
    pub midpoint_angular_speed_rad_s: f64,
    pub midpoint_residual_rad_s: f64,
    /// Bearing impulses ON the rotor. The couple is about its fixed COM.
    pub guide_impulse_kg_m_s: Vec3,
    pub guide_couple_kg_m2_s: Vec3,
    /// Actual measured guide work; zero within the native principal-axis gate.
    pub guide_work_j: f64,
    pub balance_origin_m: Vec3,
    pub fluid_linear_change_kg_m_s: Vec3,
    pub fluid_angular_change_kg_m2_s: Vec3,
    pub rotor_angular_change_kg_m2_s: Vec3,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    /// Sums of participant norms: opposing fluid/body defects cannot cancel.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub fluid_mass_change_kg: f64,
    pub fluid_kinetic_change_j: f64,
    pub rotor_kinetic_change_j: f64,
    pub rotor_external_work_j: f64,
    pub fluid_external_work_j: f64,
    pub fluid_drag_work_j: f64,
    pub relative_dissipation_j: f64,
    /// Transport/pressure/collision remainder, not asserted to be heat/error.
    pub other_fluid_kinetic_change_j: f64,
    pub mechanical_work_residual_j: f64,
}

struct Sample {
    law: PorousRotorResistance,
    force: [f64; 3],
    /// Lattice skeleton velocity per physical rad/s; zero outside support.
    tangent: [f64; 3],
}
struct RotorTrial {
    exchange: ResistanceMomentumExchange3,
    omega: f64,
    residual: Vec<f64>,
}

impl AxisymmetricPorousRotor {
    /// Bind original populations and dry inertia. The complete cylinder must
    /// fit inside the closed box, including for an oblique axis. Voxel solids,
    /// open/periodic/moving walls and initial transverse motion refuse.
    /// Cell work is bounded before wall-link or trial storage is allocated.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fluid: BoundaryGrid3,
        body: RigidBodyState,
        properties: MassProperties,
        units: PorousLatticeUnits3,
        cylinder: PorousRotorCylinder,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PorousMomentumError> {
        poll(gate)?;
        units.admit()?;
        body.kinetic_energy(properties)?;
        let cells = fluid
            .dimensions()
            .into_iter()
            .try_fold(1usize, |n, d| n.checked_mul(d))
            .ok_or(invalid("rotor domain cell count overflow"))?;
        if cells == 0
            || cells > maximum_cells
            || fluid.boundaries() != BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6])
            || body.linear_momentum_world() != Vec3::ZERO
            || !cylinder.axis_world.is_finite()
            || (norm(cylinder.axis_world) - 1.0).abs() > 64.0 * f64::EPSILON
            || [cylinder.radius_m, cylinder.half_height_m]
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "bounded stationary closed box, fixed COM and finite unit-axis cylinder required",
            ));
        }
        let axis_body = body
            .pose()
            .orientation()
            .rotate_world_to_body(cylinder.axis_world);
        let unit_spin = properties.angular_momentum_body_checked(axis_body)?;
        let inertia = finite(dot(axis_body, unit_spin))?;
        let spin = body.angular_momentum_body();
        if inertia <= 0.0
            || norm(unit_spin.sub(axis_body.scale(inertia))) > 128.0 * f64::EPSILON * inertia
            || norm(spin.sub(axis_body.scale(dot(axis_body, spin))))
                > 128.0 * f64::EPSILON * norm(spin)
        {
            return Err(invalid(
                "rotor requires one native principal axis and axial initial spin",
            ));
        }
        let centre = body.pose().position_world();
        let dimensions = fluid.dimensions();
        for i in 0..3 {
            let ai = a(cylinder.axis_world)[i];
            let extent = finite(
                ai.abs() * cylinder.half_height_m
                    + cylinder.radius_m * (1.0 - ai * ai).max(0.0).sqrt(),
            )?;
            let lo = a(units.grid_origin_m)[i];
            let hi = finite(lo + dimensions[i] as f64 * units.cell_size_m)?;
            if a(centre)[i] - extent < lo || a(centre)[i] + extent > hi {
                return Err(invalid(
                    "entire rotation-invariant cylinder must fit inside the box",
                ));
            }
        }
        let [nx, ny, nz] = dimensions;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x % 256 == 0 {
                        poll(gate)?;
                    }
                    if fluid.is_solid(x, y, z) {
                        return Err(invalid("porous rotor excludes voxel solids"));
                    }
                }
            }
        }
        let links = fluid.boundary_links();
        poll(gate)?;
        Ok(Self {
            fluid,
            body,
            properties,
            units,
            cylinder,
            axis_body,
            inertia_kg_m2: inertia,
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
    pub const fn cylinder(&self) -> PorousRotorCylinder {
        self.cylinder
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }
    pub fn angular_speed_rad_s(&self) -> Result<f64, PorousMomentumError> {
        self.body_speed(self.body)
    }
    fn body_speed(&self, body: RigidBodyState) -> Result<f64, PorousMomentumError> {
        finite(dot(
            self.properties
                .angular_velocity_body_checked(body.angular_momentum_body())?,
            self.axis_body,
        ))
    }

    /// Advance one lattice tick and the actual axial rotor attitude together.
    /// `resistance(r_m, z_m)` authors isotropic coefficients on the declared
    /// cylinder; outside it both coefficients and skeleton velocity are zero.
    /// `force(cell)` is a separate native lattice body-force density everywhere.
    /// Each callback is sampled once per eligible cell, before any trial, and
    /// held at this tick's midpoint. Author time-varying data at that midpoint.
    /// `external_torque_nm` is a signed couple along the fixed world axis.
    ///
    /// Each scalar nonlinear residual executes the complete original fluid
    /// operator with `v_s = omega_mid axis x (x-COM)`. There is no affine pilot
    /// approximation. The solved midpoint is shared with the native rigid
    /// integrator. Source-work closure does not certify fluid kinetic-energy
    /// accuracy or rotor phase; refine space/time and the native relaxation.
    /// Failure or cancellation publishes neither candidate state.
    pub fn step_candidate(
        &self,
        mut resistance: impl FnMut(f64, f64) -> PorousRotorResistance,
        mut force: impl FnMut([usize; 3]) -> [f64; 3],
        external_torque_nm: f64,
        controls: PorousRotorControls,
        gate: &CancelGate,
    ) -> Result<PorousRotorStep, PorousMomentumError> {
        poll(gate)?;
        let c = controls.flow;
        let units = self.units.admit()?;
        let radius = self.cylinder.radius_m;
        let omega0 = self.angular_speed_rad_s()?;
        if !external_torque_nm.is_finite()
            || !c.balance_origin_m.is_finite()
            || !(1..=64).contains(&c.maximum_iterations)
            || !(2..=16384).contains(&c.maximum_native_trials)
            || c.maximum_trial_cells < self.cells
            || [
                c.velocity_scale_m_s,
                c.maximum_speed_m_s,
                c.difference_step,
                c.linear_tolerance_kg_m_s,
                c.angular_tolerance_kg_m2_s,
                c.midpoint_tolerance_m_s,
                c.position_tolerance_m,
                c.work_tolerance_j,
                controls.maximum_rotation_rad,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || !c.relative_mass_tolerance.is_finite()
            || !(0.0..1.0).contains(&c.relative_mass_tolerance)
            || c.maximum_speed_m_s / units.velocity > fs_lbm::MACH_LIMIT * fs_lbm::CS2.sqrt()
            || omega0.abs() * radius > c.maximum_speed_m_s
        {
            return Err(invalid(
                "finite rotor drive, bounded rim speed, solve and physical tolerances required",
            ));
        }
        let dt = self.units.timestep_s;
        let external = finite(external_torque_nm * dt)?;
        let next_time = finite(self.time_s + dt)?;
        let next_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or(invalid("step counter overflow"))?;
        if next_time <= self.time_s {
            return Err(invalid("physical clock increment is not representable"));
        }
        let omega_scale = finite(c.velocity_scale_m_s / radius)?;
        let impulse_scale = finite(2.0 * self.inertia_kg_m2 * omega_scale)?;
        if omega_scale <= 0.0 || impulse_scale <= 0.0 {
            return Err(invalid("rotor solve scale is not representable"));
        }
        let origin = self.body.pose().position_world();
        let lattice_origin = origin
            .sub(self.units.grid_origin_m)
            .scale(self.units.cell_size_m.recip());
        if !lattice_origin.is_finite() {
            return Err(invalid("rotor moment origin is not representable"));
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
                    let axial = finite(dot(self.cylinder.axis_world, arm))?;
                    let tangent = self.cylinder.axis_world.cross(arm);
                    let radial = finite(norm(tangent))?;
                    let inside = radial <= radius && axial.abs() <= self.cylinder.half_height_m;
                    let law = if inside {
                        resistance(radial, axial)
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
                            "native rotor coefficients must be finite nonnegative; force finite",
                        ));
                    }
                    let tangent = if inside {
                        a(tangent.scale(units.velocity.recip()))
                    } else {
                        [0.0; 3]
                    };
                    if tangent.iter().any(|x| !x.is_finite()) {
                        return Err(invalid("rotor lever arm is not representable"));
                    }
                    field.push(Sample {
                        law,
                        force,
                        tangent,
                    });
                }
            }
        }
        let scratch = RefCell::new(self.fluid.clone());
        let trials = Cell::new(0usize);
        let trial_cells = Cell::new(0usize);
        let evaluate = |x: &[f64]| -> Result<RotorTrial, PorousMomentumError> {
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
            let omega = finite((dd(omega0) + dd(x[0]) * dd(omega_scale)).to_f64())?;
            if omega.abs() * radius > c.maximum_speed_m_s
                || omega.abs() * dt > controls.maximum_rotation_rad
            {
                return Err(invalid(
                    "trial rotor speed or rotation exceeds the declared domain",
                ));
            }
            let mut fluid = scratch.borrow_mut();
            fluid.clone_from(&self.fluid);
            trials.set(count);
            trial_cells.set(work);
            let exchange = fluid.step_with_darcy_forchheimer_and_wall_momentum(
                |[x, y, z]| {
                    let s = &field[(z * ny + y) * nx + x];
                    DarcyForchheimer3 {
                        linear: LinearResistance3 {
                            principal: [s.law.darcy; 3],
                            skeleton_velocity: s.tangent.map(|t| omega * t),
                            body_force: s.force,
                            ..Default::default()
                        },
                        quadratic: s.law.forchheimer,
                    }
                },
                &self.links,
                a(lattice_origin),
                gate,
            )?;
            let axial = dot(
                self.cylinder.axis_world,
                v(exchange.skeleton_angular_impulse).scale(units.angular),
            );
            let residual = ((dd(2.0) * dd(self.inertia_kg_m2) * (dd(omega) - dd(omega0))
                - dd(axial)
                - dd(external))
                / dd(impulse_scale))
            .to_f64();
            Ok(RotorTrial {
                exchange,
                omega,
                residual: vec![residual],
            })
        };
        let problem = FallibleFiniteDifference::new(1, c.difference_step, |x: &[f64]| {
            evaluate(x).map(|t| t.residual)
        })
        .map_err(difference_error)?;
        let tolerance = 0.25
            * (c.angular_tolerance_kg_m2_s / impulse_scale)
                .min(c.midpoint_tolerance_m_s / c.velocity_scale_m_s);
        let config = NewtonKrylovConfig {
            absolute_tolerance: tolerance,
            relative_tolerance: 0.0,
            linear_restart: 1,
            max_linear_cycles: 4,
            forcing_minimum: 1e-3,
            forcing_maximum: 0.1,
            ..Default::default()
        };
        let construction = NewtonKrylovState::new(&problem, vec![0.0], config);
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
        // Charge the final replay, and keep its exact native populations.
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
        trial: RotorTrial,
        external: f64,
        controls: PorousRotorControls,
        units: Units,
        nonlinear: NewtonReport,
        native_trials: usize,
        trial_cells: usize,
        time_s: f64,
        accepted_steps: u64,
        gate: &CancelGate,
    ) -> Result<PorousRotorStep, PorousMomentumError> {
        let c = controls.flow;
        let exchange = trial.exchange;
        let axis = self.cylinder.axis_world;
        let impulse = v(exchange.skeleton_impulse).scale(units.impulse);
        let moment = v(exchange.skeleton_angular_impulse).scale(units.angular);
        let axial = finite(dot(axis, moment))?;
        let total = dd(axial) + dd(external);
        let torque = finite((total / dd(self.units.timestep_s)).to_f64())?;
        let mechanics = RigidBodyIntegrator::new(Gravity::ZERO).step(
            self.body,
            self.properties,
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: self.axis_body.scale(torque),
            },
            self.units.timestep_s,
        )?;
        let body = mechanics.state_after;
        let origin = self.body.pose().position_world();
        if body.pose().position_world() != origin
            || body.linear_momentum_world() != Vec3::ZERO
            || norm(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(self.axis_body)
                    .sub(axis),
            ) > 256.0 * f64::EPSILON
        {
            return Err(invalid(
                "native rotor departed from its fixed COM/principal axis",
            ));
        }
        let omega0 = self
            .properties
            .angular_velocity_body_checked(self.body.angular_momentum_body())?;
        let omega1 = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        let before_world = self.body.pose().orientation().rotate_body_to_world(omega0);
        let after_world = body.pose().orientation().rotate_body_to_world(omega1);
        let midpoint = v(std::array::from_fn(|i| {
            ((dd(a(before_world)[i]) + dd(a(after_world)[i])) * dd(0.5)).to_f64()
        }));
        let midpoint_error = midpoint.sub(axis.scale(trial.omega));
        check(
            "rotor shared midpoint rim speed",
            norm(midpoint_error) * self.cylinder.radius_m,
            c.midpoint_tolerance_m_s,
        )?;
        if norm(after_world) * self.cylinder.radius_m > c.maximum_speed_m_s
            || norm(midpoint) * self.units.timestep_s > controls.maximum_rotation_rad
        {
            return Err(invalid(
                "accepted rotor speed or rotation exceeds the declared domain",
            ));
        }
        let before_spin = self
            .body
            .pose()
            .orientation()
            .rotate_body_to_world(self.body.angular_momentum_body());
        let after_spin = body
            .pose()
            .orientation()
            .rotate_body_to_world(body.angular_momentum_body());
        let rotor_change = v(std::array::from_fn(|i| {
            (dd(a(after_spin)[i]) - dd(a(before_spin)[i])).to_f64()
        }));
        let kinetic = (0..3)
            .fold(Dd::ZERO, |sum, i| {
                sum + dd(0.5)
                    * (dd(a(omega0)[i]) + dd(a(omega1)[i]))
                    * (dd(a(body.angular_momentum_body())[i])
                        - dd(a(self.body.angular_momentum_body())[i]))
            })
            .to_f64();
        // Bearing force and transverse couple follow the measured drag wrench,
        // independently of the measured rotor increment or balance residual.
        let guide = impulse.scale(-1.0);
        let guide_couple = axis.scale(axial).sub(moment);
        let arm = origin.sub(c.balance_origin_m);
        let shift = |moment: Vec3, impulse: Vec3| {
            moment.add(add_cross(arm, impulse, Vec3::ZERO, Vec3::ZERO))
        };
        let drag_angular = shift(moment, impulse);
        let guide_angular = shift(guide_couple, guide);
        let external_angular = axis.scale(external);
        let body_residual = v(std::array::from_fn(|i| {
            (dd(a(rotor_change)[i])
                - dd(a(drag_angular)[i])
                - dd(a(guide_angular)[i])
                - dd(a(external_angular)[i]))
            .to_f64()
        }));
        let change = fluid_change(&self.fluid, &fluid, self.units, c, gate)?;
        let fluid_external = v(exchange.external_impulse).scale(units.impulse);
        let fluid_external_angular = shift(
            v(exchange.external_angular_impulse).scale(units.angular),
            fluid_external,
        );
        let wall_fluid = v(exchange.walls.fluid_population_impulse).scale(units.impulse);
        let wall_angular = shift(
            v(exchange.walls.fluid_population_angular_impulse).scale(units.angular),
            wall_fluid,
        );
        let fluid_linear = v(std::array::from_fn(|i| {
            (dd(a(change.linear)[i]) - dd(a(fluid_external)[i]) + dd(a(impulse)[i])
                - dd(a(wall_fluid)[i]))
            .to_f64()
        }));
        let fluid_angular = v(std::array::from_fn(|i| {
            (dd(a(change.angular)[i]) - dd(a(fluid_external_angular)[i]) + dd(a(drag_angular)[i])
                - dd(a(wall_angular)[i]))
            .to_f64()
        }));
        // The native rotor's unchanged linear momentum equals drag + bearing.
        let body_linear = impulse.add(guide).scale(-1.0);
        let absolute_linear = norm(body_linear) + norm(fluid_linear);
        let absolute_angular = norm(body_residual) + norm(fluid_angular);
        check(
            "rotor/fluid linear momentum including bearing and walls",
            absolute_linear,
            c.linear_tolerance_kg_m_s,
        )?;
        check(
            "rotor/fluid angular momentum including bearing and walls",
            absolute_angular,
            c.angular_tolerance_kg_m2_s,
        )?;
        let fluid_work = finite(exchange.external_work * units.energy)?;
        let drag_work = finite(exchange.fluid_drag_work * units.energy)?;
        let skeleton_work = finite(exchange.skeleton_work * units.energy)?;
        let loss = finite(exchange.relative_dissipation * units.energy)?;
        if loss < 0.0 {
            return Err(invalid("native porous dissipation is negative"));
        }
        let external_work = dot(midpoint, external_angular);
        let guide_work = dot(midpoint, guide_couple);
        let rotor_drag_work = dot(midpoint, moment);
        let constitutive = (dd(drag_work) + dd(skeleton_work) + dd(loss)).to_f64();
        let interface = (dd(rotor_drag_work) - dd(skeleton_work)).to_f64();
        let rigid_work =
            (dd(kinetic) - dd(rotor_drag_work) - dd(guide_work) - dd(external_work)).to_f64();
        check(
            "noncancelling rotor constitutive/interface/native work",
            constitutive.abs() + interface.abs() + rigid_work.abs() + guide_work.abs(),
            c.work_tolerance_j,
        )?;
        let mechanical =
            (dd(kinetic) + dd(drag_work) + dd(loss) - dd(external_work) - dd(guide_work)).to_f64();
        let remainder = (dd(change.kinetic) - dd(fluid_work) - dd(drag_work)).to_f64();
        let linear_residual = body_linear.add(fluid_linear);
        let angular_residual = body_residual.add(fluid_angular);
        if [
            impulse,
            moment,
            guide,
            guide_couple,
            rotor_change,
            linear_residual,
            angular_residual,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || [kinetic, external_work, guide_work, mechanical, remainder]
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(invalid("combined rotor account is not representable"));
        }
        poll(gate)?;
        Ok(PorousRotorStep {
            motion: Self {
                fluid,
                body,
                properties: self.properties,
                units: self.units,
                cylinder: self.cylinder,
                axis_body: self.axis_body,
                inertia_kg_m2: self.inertia_kg_m2,
                links: self.links.clone(),
                cells: self.cells,
                time_s,
                accepted_steps,
            },
            exchange,
            exchange_origin_m: origin,
            mechanics,
            nonlinear,
            native_trials,
            trial_cells,
            midpoint_angular_speed_rad_s: trial.omega,
            midpoint_residual_rad_s: dot(axis, midpoint_error),
            guide_impulse_kg_m_s: guide,
            guide_couple_kg_m2_s: guide_couple,
            guide_work_j: guide_work,
            balance_origin_m: c.balance_origin_m,
            fluid_linear_change_kg_m_s: change.linear,
            fluid_angular_change_kg_m2_s: change.angular,
            rotor_angular_change_kg_m2_s: rotor_change,
            linear_residual_kg_m_s: linear_residual,
            angular_residual_kg_m2_s: angular_residual,
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            fluid_mass_change_kg: change.mass,
            fluid_kinetic_change_j: change.kinetic,
            rotor_kinetic_change_j: kinetic,
            rotor_external_work_j: external_work,
            fluid_external_work_j: fluid_work,
            fluid_drag_work_j: drag_work,
            relative_dissipation_j: loss,
            other_fluid_kinetic_change_j: remainder,
            mechanical_work_residual_j: mechanical,
        })
    }
}

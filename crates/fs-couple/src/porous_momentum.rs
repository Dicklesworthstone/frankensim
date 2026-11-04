//! Finite dry-skeleton recoil under the original periodic Darcy--Forchheimer flow.
//!
//! One homogeneous porous continuum translates through a representative periodic
//! domain. Its dry mass excludes the fluid already represented by populations.
//! Uniform coefficients and fixed material axes make support invariant under
//! translation. An ideal rotational guide holds orientation and receives the
//! actual distributed drag couple. This is not a moving bounded porous object,
//! a porosity/added-mass model, or a substitute for native immersed-body motion.
//!
//! Three midpoint velocities are solved with the original Newton--Krylov owner.
//! Every residual executes the original nonlinear fluid step from the accepted
//! populations; transport-dependent density is never replaced by an affine
//! response. Native drag quadrature is retained, including its stiff-step slip
//! reversal. Constitutive work closes separately from transported fluid kinetic
//! energy. No thermal store, heat partition, global energy theorem or trajectory
//! error certificate is inferred. Resolve spatial and relaxation timescales.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, BoundaryStepError3, DarcyForchheimer3,
    ResistanceMomentumExchange3,
};
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::nonlinear::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Unit correspondence of the EXISTING lattice and the native SI rigid state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousLatticeUnits3 {
    pub cell_size_m: f64,
    pub timestep_s: f64,
    /// Fluid mass in one cell at lattice density one, kg. The caller retains
    /// the original density/volume convention; no porosity factor is inferred.
    pub mass_per_density_cell_kg: f64,
    /// World position of grid corner [0,0,0]; grid and world axes coincide.
    pub grid_origin_m: Vec3,
}

#[derive(Clone, Copy)]
struct Units {
    velocity: f64,
    impulse: f64,
    angular: f64,
    energy: f64,
}
impl PorousLatticeUnits3 {
    fn admit(self) -> Result<Units, PorousMomentumError> {
        if !self.grid_origin_m.is_finite()
            || [
                self.cell_size_m,
                self.timestep_s,
                self.mass_per_density_cell_kg,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid("positive finite original lattice scales required"));
        }
        let velocity = self.cell_size_m / self.timestep_s;
        let impulse = self.mass_per_density_cell_kg * velocity;
        let angular = impulse * self.cell_size_m;
        let energy = impulse * velocity;
        if [velocity, impulse, angular, energy]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid("lattice correspondence is not representable"));
        }
        Ok(Units {
            velocity,
            impulse,
            angular,
            energy,
        })
    }
}

/// Bounded native work and independent physical acceptance tolerances.
#[derive(Clone, Copy, Debug)]
pub struct PorousMomentumControls {
    /// Scale of the three Newton velocity-increment coordinates, m/s.
    pub velocity_scale_m_s: f64,
    /// Applies to initial/final fluid and body velocities and every trial body
    /// velocity. Must also fit the original lattice Mach limit.
    pub maximum_speed_m_s: f64,
    pub difference_step: f64,
    /// Newton attempts in 1..=64; every stencil/trial is additionally charged.
    pub maximum_iterations: usize,
    pub maximum_native_trials: usize,
    pub maximum_trial_cells: usize,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub midpoint_tolerance_m_s: f64,
    pub position_tolerance_m: f64,
    pub work_tolerance_j: f64,
    pub relative_mass_tolerance: f64,
    pub balance_origin_m: Vec3,
}

#[derive(Debug)]
pub enum PorousMomentumError {
    Invalid(&'static str),
    Fluid(BoundaryStepError3),
    Mechanics(DynamicsError),
    Solver(NewtonError),
    SolverProgress(SolverRunProgress),
    Unresolved(NewtonReport),
    Budget {
        trials: usize,
        trial_cells: usize,
    },
    Balance {
        account: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for PorousMomentumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "periodic porous recoil: {self:?}")
    }
}
impl std::error::Error for PorousMomentumError {}
impl From<BoundaryStepError3> for PorousMomentumError {
    fn from(e: BoundaryStepError3) -> Self {
        Self::Fluid(e)
    }
}
impl From<DynamicsError> for PorousMomentumError {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanics(e)
    }
}

/// Original periodic fluid and one nonrotating finite dry-skeleton momentum.
pub struct PeriodicPorousSkeleton {
    fluid: BoundaryGrid3,
    body: RigidBodyState,
    properties: MassProperties,
    units: PorousLatticeUnits3,
    cells: usize,
    time_s: f64,
    accepted_steps: u64,
}

/// Actual accepted native candidates, exchange and independent physical checks.
pub struct PorousMomentumStep {
    pub motion: PeriodicPorousSkeleton,
    /// Original lattice-unit receipt about `exchange_origin_m`.
    pub exchange: ResistanceMomentumExchange3,
    pub exchange_origin_m: Vec3,
    pub momentum: MomentumImpulseReceipt,
    pub nonlinear: NewtonReport,
    pub native_trials: usize,
    pub trial_cells: usize,
    pub midpoint_velocity_m_s: Vec3,
    pub midpoint_residual_m_s: Vec3,
    /// Ideal guide couple ON the skeleton; zero work because spin stays zero.
    pub guide_angular_impulse_kg_m2_s: Vec3,
    pub periodic_angular_transport_kg_m2_s: Vec3,
    pub balance_origin_m: Vec3,
    pub fluid_linear_change_kg_m_s: Vec3,
    pub fluid_angular_change_kg_m2_s: Vec3,
    pub skeleton_angular_change_kg_m2_s: Vec3,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    /// Sum of participant residual norms, so cancellation cannot pass the gate.
    pub absolute_linear_residual_kg_m_s: f64,
    pub absolute_angular_residual_kg_m2_s: f64,
    pub fluid_mass_change_kg: f64,
    pub fluid_kinetic_change_j: f64,
    pub skeleton_kinetic_change_j: f64,
    pub skeleton_external_work_j: f64,
    pub fluid_external_work_j: f64,
    pub fluid_drag_work_j: f64,
    pub relative_dissipation_j: f64,
    /// Fluid kinetic change minus native independent/drag source work. Includes
    /// pressure, transport and collision; not reclassified as heat or error.
    pub other_fluid_kinetic_change_j: f64,
    /// Skeleton kinetic change + fluid drag work + dissipation - skeleton drive.
    pub mechanical_work_residual_j: f64,
}

struct Trial {
    exchange: ResistanceMomentumExchange3,
    velocity: Vec3,
    residual: Vec<f64>,
}

impl PeriodicPorousSkeleton {
    /// Bind one homogeneous periodic representative domain. The source's actual
    /// dry mass and inertia remain native; density-one cell mass belongs solely
    /// to the original fluid. No body occupancy or displaced fluid is invented.
    /// Domain size is bounded and solids/open/wall faces refuse before cloning.
    pub fn new(
        fluid: BoundaryGrid3,
        body: RigidBodyState,
        properties: MassProperties,
        units: PorousLatticeUnits3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PorousMomentumError> {
        poll(gate)?;
        units.admit()?;
        body.center_of_mass_velocity_world(properties)?;
        let cells = fluid
            .dimensions()
            .into_iter()
            .try_fold(1usize, |a, b| a.checked_mul(b))
            .ok_or(invalid("domain cell count overflow"))?;
        if cells == 0
            || cells > maximum_cells
            || fluid.boundaries() != BoundarySpec3::periodic()
            || body.angular_momentum_body() != Vec3::ZERO
        {
            return Err(invalid(
                "bounded fully periodic fluid and zero skeleton spin required",
            ));
        }
        let [nx, ny, nz] = fluid.dimensions();
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if x % 256 == 0 {
                        poll(gate)?;
                    }
                    if fluid.is_solid(x, y, z) {
                        return Err(invalid("homogeneous porous support excludes voxel solids"));
                    }
                }
            }
        }
        Ok(Self {
            fluid,
            body,
            properties,
            units,
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
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Execute one original lattice tick with simultaneous finite-body recoil.
    /// `law` is homogeneous and uses original lattice units; its prescribed
    /// skeleton velocity must be zero because this owner solves that velocity.
    /// `external_force_n` is a held force on the dry skeleton at its COM.
    ///
    /// Coefficients/axes and fluid body force are held at the step midpoint;
    /// support is translationally invariant. The original endpoint-averaged
    /// fluid law uses the solved constant skeleton midpoint velocity. This is
    /// second-order coupling in smooth resolved regimes, not an exact nonlinear
    /// relaxation or a joint endpoint evaluation of changing skeleton speeds.
    /// Failed native trials, solves, balance checks and cancellation publish no
    /// state. One scratch grid reuses the original native clone_from buffers.
    pub fn step_candidate(
        &self,
        law: DarcyForchheimer3,
        external_force_n: Vec3,
        c: PorousMomentumControls,
        gate: &CancelGate,
    ) -> Result<PorousMomentumStep, PorousMomentumError> {
        poll(gate)?;
        let units = self.units.admit()?;
        let before_v = self.body.center_of_mass_velocity_world(self.properties)?;
        if law.linear.skeleton_velocity != [0.0; 3]
            || !external_force_n.is_finite()
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
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || !c.relative_mass_tolerance.is_finite()
            || !(0.0..1.0).contains(&c.relative_mass_tolerance)
            || c.maximum_speed_m_s / units.velocity > fs_lbm::MACH_LIMIT * fs_lbm::CS2.sqrt()
            || norm(before_v) > c.maximum_speed_m_s
        {
            return Err(invalid(
                "finite source, bounded solve, explicit units and physical tolerances required",
            ));
        }
        let dt = self.units.timestep_s;
        let next_time = finite(self.time_s + dt)?;
        let next_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or(invalid("step counter overflow"))?;
        if next_time <= self.time_s {
            return Err(invalid("physical clock increment is not representable"));
        }
        let external = external_force_n.scale(dt);
        let origin = self
            .body
            .pose()
            .position_world()
            .add(before_v.scale(0.5 * dt));
        let lattice_origin = origin
            .sub(self.units.grid_origin_m)
            .scale(self.units.cell_size_m.recip());
        let momentum_scale = finite(2.0 * self.properties.mass() * c.velocity_scale_m_s)?;
        if !external.is_finite()
            || !origin.is_finite()
            || !lattice_origin.is_finite()
            || momentum_scale <= 0.0
        {
            return Err(invalid(
                "native impulse correspondence is not representable",
            ));
        }
        // Also validate the accepted physical fluid before creating trial storage.
        let _ = fluid_change(&self.fluid, &self.fluid, self.units, c, gate)?;
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
            let velocity = before_v.add(v([x[0], x[1], x[2]]).scale(c.velocity_scale_m_s));
            if !velocity.is_finite() || norm(velocity) > c.maximum_speed_m_s {
                return Err(invalid(
                    "trial skeleton velocity exceeds the declared domain",
                ));
            }
            let mut selected = law;
            selected.linear.skeleton_velocity = a(velocity.scale(units.velocity.recip()));
            let mut fluid = scratch.borrow_mut();
            fluid.clone_from(&self.fluid);
            trials.set(count);
            trial_cells.set(work);
            let exchange = fluid.step_with_darcy_forchheimer_and_wall_momentum(
                |_| selected,
                &[],
                a(lattice_origin),
                gate,
            )?;
            let impulse = v(exchange.skeleton_impulse).scale(units.impulse);
            let residual = (0..3)
                .map(|i| {
                    ((dd(2.0)
                        * dd(self.properties.mass())
                        * (dd(a(velocity)[i]) - dd(a(before_v)[i]))
                        - dd(a(impulse)[i])
                        - dd(a(external)[i]))
                        / dd(momentum_scale))
                    .to_f64()
                })
                .collect();
            Ok(Trial {
                exchange,
                velocity,
                residual,
            })
        };
        let problem = FallibleFiniteDifference::new(3, c.difference_step, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual)
        })
        .map_err(difference_error)?;
        let tolerance = 0.25
            * (c.linear_tolerance_kg_m_s / momentum_scale)
                .min(c.midpoint_tolerance_m_s / c.velocity_scale_m_s);
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
        // This charged final native replay retains the exact accepted populations.
        let trial = evaluate(&solver.x)?;
        let mut nonlinear = result.report;
        nonlinear.residual_norm = fs_solver::norm2(&trial.residual);
        if !nonlinear.residual_norm.is_finite() || nonlinear.residual_norm > tolerance {
            nonlinear.converged = false;
            return Err(PorousMomentumError::Unresolved(nonlinear));
        }
        let fluid = scratch.into_inner();
        self.finish(
            fluid,
            trial,
            external,
            origin,
            units,
            c,
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
        origin: Vec3,
        units: Units,
        c: PorousMomentumControls,
        nonlinear: NewtonReport,
        native_trials: usize,
        trial_cells: usize,
        time_s: f64,
        accepted_steps: u64,
        gate: &CancelGate,
    ) -> Result<PorousMomentumStep, PorousMomentumError> {
        let exchange = trial.exchange;
        let impulse = v(exchange.skeleton_impulse).scale(units.impulse);
        let momentum =
            self.body
                .apply_momentum_impulse(self.properties, impulse.add(external), Vec3::ZERO)?;
        let before_v = self.body.center_of_mass_velocity_world(self.properties)?;
        let after_v = momentum
            .state_after
            .center_of_mass_velocity_world(self.properties)?;
        let midpoint = before_v.scale(0.5).add(after_v.scale(0.5));
        let midpoint_residual = midpoint.sub(trial.velocity);
        check(
            "shared midpoint",
            norm(midpoint_residual),
            c.midpoint_tolerance_m_s,
        )?;
        if norm(after_v) > c.maximum_speed_m_s {
            return Err(invalid(
                "accepted skeleton speed exceeds the declared domain",
            ));
        }
        let displacement = midpoint.scale(self.units.timestep_s);
        let old_position = self.body.pose().position_world();
        let pose = Pose::new(
            old_position.add(displacement),
            self.body.pose().orientation(),
        )?;
        let represented = pose.position_world().sub(old_position);
        check(
            "skeleton displacement",
            norm(represented.sub(displacement)),
            c.position_tolerance_m,
        )?;
        let body = RigidBodyState::new(
            pose,
            momentum.state_after.linear_momentum_world(),
            Vec3::ZERO,
        )?;
        let midpoint_position = old_position.add(represented.scale(0.5));
        let source_moment = v(exchange.skeleton_angular_impulse).scale(units.angular);
        // The native receipt may use any fixed moment origin. The physical
        // rotational guide cancels drag about the ACTUAL midpoint COM; the
        // independent drive acts there too, even when it is not parallel to drag.
        let guide = source_moment
            .add(add_cross(
                origin.sub(midpoint_position),
                impulse,
                Vec3::ZERO,
                Vec3::ZERO,
            ))
            .scale(-1.0);
        let periodic = v(exchange.periodic.angular_momentum_change).scale(units.angular);
        let body_angular = add_cross(
            old_position.sub(c.balance_origin_m),
            momentum.linear_increment_world,
            represented,
            body.linear_momentum_world(),
        );
        let arm = origin.sub(c.balance_origin_m);
        let drag_angular = source_moment.add(add_cross(arm, impulse, Vec3::ZERO, Vec3::ZERO));
        let applied_angular = drag_angular.add(guide).add(add_cross(
            midpoint_position.sub(c.balance_origin_m),
            external,
            Vec3::ZERO,
            Vec3::ZERO,
        ));
        let body_linear_residual = v(std::array::from_fn(|i| {
            (dd(a(momentum.linear_increment_world)[i]) - dd(a(impulse)[i]) - dd(a(external)[i]))
                .to_f64()
        }));
        let body_angular_residual = body_angular.sub(applied_angular);
        let change = fluid_change(&self.fluid, &fluid, self.units, c, gate)?;
        let fluid_external = v(exchange.external_impulse).scale(units.impulse);
        let fluid_external_angular = v(exchange.external_angular_impulse)
            .scale(units.angular)
            .add(arm.cross(fluid_external));
        let fluid_linear_residual = v(std::array::from_fn(|i| {
            (dd(a(change.linear)[i]) - dd(a(fluid_external)[i]) + dd(a(impulse)[i])).to_f64()
        }));
        let fluid_angular_residual = change
            .angular
            .sub(fluid_external_angular)
            .add(drag_angular)
            .sub(periodic);
        let absolute_linear = norm(body_linear_residual) + norm(fluid_linear_residual);
        let absolute_angular = norm(body_angular_residual) + norm(fluid_angular_residual);
        check(
            "pair linear momentum",
            absolute_linear,
            c.linear_tolerance_kg_m_s,
        )?;
        check(
            "pair angular momentum including guide and periodic cuts",
            absolute_angular,
            c.angular_tolerance_kg_m2_s,
        )?;
        let external_work = dot(midpoint, external);
        let fluid_work = finite(exchange.external_work * units.energy)?;
        let drag_work = finite(exchange.fluid_drag_work * units.energy)?;
        let skeleton_work = finite(exchange.skeleton_work * units.energy)?;
        let loss = finite(exchange.relative_dissipation * units.energy)?;
        if loss < 0.0 {
            return Err(invalid("native porous dissipation is negative"));
        }
        let constitutive = (dd(drag_work) + dd(skeleton_work) + dd(loss)).to_f64();
        let interface = (dd(dot(midpoint, impulse)) - dd(skeleton_work)).to_f64();
        let native_work = momentum.work.work_energy_residual_j;
        let summed_impulse_work = (0..3)
            .fold(Dd::ZERO, |sum, i| {
                sum + dd(a(midpoint)[i])
                    * (dd(a(momentum.linear_impulse_world)[i])
                        - dd(a(impulse)[i])
                        - dd(a(external)[i]))
            })
            .to_f64();
        check(
            "noncancelling constitutive/rigid interface work",
            constitutive.abs() + interface.abs() + native_work.abs() + summed_impulse_work.abs(),
            c.work_tolerance_j,
        )?;
        let mechanical =
            (dd(constitutive) + dd(interface) + dd(native_work) + dd(summed_impulse_work)).to_f64();
        let report = PorousMomentumStep {
            motion: Self {
                fluid,
                body,
                properties: self.properties,
                units: self.units,
                cells: self.cells,
                time_s,
                accepted_steps,
            },
            exchange,
            exchange_origin_m: origin,
            momentum,
            nonlinear,
            native_trials,
            trial_cells,
            midpoint_velocity_m_s: trial.velocity,
            midpoint_residual_m_s: midpoint_residual,
            guide_angular_impulse_kg_m2_s: guide,
            periodic_angular_transport_kg_m2_s: periodic,
            balance_origin_m: c.balance_origin_m,
            fluid_linear_change_kg_m_s: change.linear,
            fluid_angular_change_kg_m2_s: change.angular,
            skeleton_angular_change_kg_m2_s: body_angular,
            linear_residual_kg_m_s: body_linear_residual.add(fluid_linear_residual),
            angular_residual_kg_m2_s: body_angular_residual.add(fluid_angular_residual),
            absolute_linear_residual_kg_m_s: absolute_linear,
            absolute_angular_residual_kg_m2_s: absolute_angular,
            fluid_mass_change_kg: change.mass,
            fluid_kinetic_change_j: change.kinetic,
            skeleton_kinetic_change_j: momentum.work.kinetic_energy_change_j,
            skeleton_external_work_j: external_work,
            fluid_external_work_j: fluid_work,
            fluid_drag_work_j: drag_work,
            relative_dissipation_j: loss,
            other_fluid_kinetic_change_j: (dd(change.kinetic) - dd(fluid_work) - dd(drag_work))
                .to_f64(),
            mechanical_work_residual_j: mechanical,
        };
        if !guide.is_finite()
            || !periodic.is_finite()
            || !report.other_fluid_kinetic_change_j.is_finite()
            || !external_work.is_finite()
            || !mechanical.is_finite()
        {
            return Err(invalid("combined porous balance is not representable"));
        }
        poll(gate)?;
        Ok(report)
    }
}

struct FluidChange {
    mass: f64,
    linear: Vec3,
    angular: Vec3,
    kinetic: f64,
}
fn fluid_change(
    before: &BoundaryGrid3,
    after: &BoundaryGrid3,
    units: PorousLatticeUnits3,
    c: PorousMomentumControls,
    gate: &CancelGate,
) -> Result<FluidChange, PorousMomentumError> {
    let scale = units.admit()?;
    let mut mass0 = Dd::ZERO;
    let mut mass = Dd::ZERO;
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    let mut kinetic = Dd::ZERO;
    let [nx, ny, nz] = before.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if x % 256 == 0 {
                    poll(gate)?;
                }
                let r0 = before.density(x, y, z);
                let r1 = after.density(x, y, z);
                let u0 = before.velocity(x, y, z);
                let u1 = after.velocity(x, y, z);
                if ![r0, r1].iter().all(|r| r.is_finite() && *r > 0.0)
                    || u0.iter().chain(&u1).any(|v| !v.is_finite())
                    || norm(v(u0)) * scale.velocity > c.maximum_speed_m_s
                    || norm(v(u1)) * scale.velocity > c.maximum_speed_m_s
                {
                    return Err(invalid(
                        "accepted fluid density/velocity is outside the declared domain",
                    ));
                }
                mass0 = mass0 + dd(r0);
                mass = mass + dd(r1) - dd(r0);
                let dp: [Dd; 3] = std::array::from_fn(|i| {
                    (dd(r1) * dd(u1[i]) - dd(r0) * dd(u0[i])) * dd(scale.impulse)
                });
                let cell = [x, y, z];
                let arm: [Dd; 3] = std::array::from_fn(|i| {
                    dd(a(units.grid_origin_m)[i]) - dd(a(c.balance_origin_m)[i])
                        + (dd(cell[i] as f64) + dd(0.5)) * dd(units.cell_size_m)
                });
                let moment = cross(arm, dp);
                for i in 0..3 {
                    linear[i] = linear[i] + dp[i];
                    angular[i] = angular[i] + moment[i];
                    kinetic = kinetic
                        + dd(0.5)
                            * dd(scale.energy)
                            * (dd(r1) * dd(u1[i]) * dd(u1[i]) - dd(r0) * dd(u0[i]) * dd(u0[i]));
                }
            }
        }
    }
    let mass_change = (mass * dd(units.mass_per_density_cell_kg)).to_f64();
    check(
        "periodic fluid mass",
        mass_change.abs(),
        (mass0 * dd(units.mass_per_density_cell_kg) * dd(c.relative_mass_tolerance)).to_f64(),
    )?;
    let result = FluidChange {
        mass: mass_change,
        linear: v(linear.map(Dd::to_f64)),
        angular: v(angular.map(Dd::to_f64)),
        kinetic: kinetic.to_f64(),
    };
    if !result.linear.is_finite() || !result.angular.is_finite() || !result.kinetic.is_finite() {
        return Err(invalid(
            "actual fluid momentum/kinetic change is not representable",
        ));
    }
    Ok(result)
}

fn add_cross(a0: Vec3, b0: Vec3, a1: Vec3, b1: Vec3) -> Vec3 {
    let c0 = cross(a(a0).map(dd), a(b0).map(dd));
    let c1 = cross(a(a1).map(dd), a(b1).map(dd));
    v(std::array::from_fn(|i| (c0[i] + c1[i]).to_f64()))
}
fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(x: Vec3, y: Vec3) -> f64 {
    (dd(x.x) * dd(y.x) + dd(x.y) * dd(y.y) + dd(x.z) * dd(y.z)).to_f64()
}
fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn a(x: Vec3) -> [f64; 3] {
    [x.x, x.y, x.z]
}
fn v(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn norm(x: Vec3) -> f64 {
    x.x.hypot(x.y).hypot(x.z)
}
fn invalid(s: &'static str) -> PorousMomentumError {
    PorousMomentumError::Invalid(s)
}
fn finite(x: f64) -> Result<f64, PorousMomentumError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("nonfinite coupled porous arithmetic"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), PorousMomentumError> {
    if gate.is_requested() {
        Err(PorousMomentumError::Cancelled)
    } else {
        Ok(())
    }
}
fn check(account: &'static str, residual: f64, allowance: f64) -> Result<(), PorousMomentumError> {
    if !residual.is_finite() || !allowance.is_finite() || residual > allowance {
        Err(PorousMomentumError::Balance {
            account,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}
fn difference_error(e: FiniteDifferenceError<PorousMomentumError>) -> PorousMomentumError {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => invalid(s),
    }
}

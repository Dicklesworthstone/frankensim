//! Finite cohesive patches evolving two original rigid bodies and their bonds.
//!
//! Each source-bound material point keeps its original area, damage, viscosity,
//! body attachment arms and a frame fixed in body A. All sites participate in
//! one twelve-variable wrench solve through the original `fs-mbd` inertia and
//! Cayley midpoint integrator. No independently moving quadrature masses or
//! substitute contact springs are introduced. Work per nonlinear evaluation is
//! linear in the number of sites; the prepared dense solve remains 12 by 12.
//!
//! The original objective constitutive owner supplies the actual local-path
//! mean force F and frame couple M. With midpoint spatial arms r_a, r_b, the
//! bodies receive `(F, r_a x F + M)` and `(-F, -r_b x F)`. The frame couple is
//! essential: these internal forces have zero total spatial moment even for a
//! finite separation. A spatial moment is pulled back by the inverse ARITHMETIC
//! rotation average, not a midpoint quaternion. This makes both native angular
//! impulse and torque work agree with the original objective work ports.
//!
//! Accepted poses, momenta, every material point, observations and the clock
//! publish together after independent original energy/momentum checks. Failure
//! modifies only reusable scratch. Source conditions and attachment topology
//! stay fixed. Complete failure preserves the fragments and irreversible
//! history; the original normal compression penalty allows frictionless
//! reclosure. No bond creation, healing, collision search, damaged-area friction,
//! global softening-branch uniqueness or exact intrastep trajectory is claimed.

use super::objective::{
    ObjectiveCohesiveConfig, ObjectiveCohesiveError, cayley_increment,
    source::{
        SourceBoundObjectiveCohesiveError, SourceBoundObjectiveCohesivePoint,
        SourceBoundObjectiveCohesiveStep,
    },
};
use fs_exec::CancelGate;
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_matdb::MaterialStateId;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, UnitQuaternion,
    Vec3, Wrench,
};

mod solve;

/// The same physical state bounds apply at admission and after each real step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveBodyLimits {
    pub maximum_linear_speed_m_s: f64,
    pub maximum_angular_speed_rad_s: f64,
    pub maximum_kinetic_energy_j: f64,
    pub maximum_step_displacement_m: f64,
}

/// Original native body state and full inertia, with an explicit material owner.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidCohesiveBody {
    material: MaterialStateId,
    properties: MassProperties,
    state: RigidBodyState,
    limits: RigidCohesiveBodyLimits,
}
impl RigidCohesiveBody {
    pub fn new(
        material: MaterialStateId,
        properties: MassProperties,
        state: RigidBodyState,
        limits: RigidCohesiveBodyLimits,
    ) -> Result<Self, RigidCohesiveError> {
        for value in [
            limits.maximum_linear_speed_m_s,
            limits.maximum_angular_speed_rad_s,
            limits.maximum_kinetic_energy_j,
            limits.maximum_step_displacement_m,
        ] {
            positive(value)?;
        }
        let body = Self {
            material,
            properties,
            state,
            limits,
        };
        body.check(state)?;
        Ok(body)
    }
    pub const fn material(&self) -> &MaterialStateId {
        &self.material
    }
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    pub const fn state(&self) -> RigidBodyState {
        self.state
    }
    pub const fn limits(&self) -> RigidCohesiveBodyLimits {
        self.limits
    }

    fn check(&self, state: RigidBodyState) -> Result<(), RigidCohesiveError> {
        bound(
            "body linear speed",
            norm(state.center_of_mass_velocity_world(self.properties)?),
            self.limits.maximum_linear_speed_m_s,
        )?;
        bound(
            "body angular speed",
            norm(
                self.properties
                    .angular_velocity_body_checked(state.angular_momentum_body())?,
            ),
            self.limits.maximum_angular_speed_rad_s,
        )?;
        bound(
            "body kinetic energy",
            state.kinetic_energy(self.properties)?,
            self.limits.maximum_kinetic_energy_j,
        )?;
        bound(
            "body step displacement",
            norm(
                state
                    .pose()
                    .position_world()
                    .sub(self.state.pose().position_world()),
            ),
            self.limits.maximum_step_displacement_m,
        )
    }
}

/// One existing original bond at two material attachment points.
///
/// `frame_in_a` maps local `[normal, shear_1, shear_2]` into body A. Construction
/// of the system verifies that this frame and both actual body attachments
/// represent the supplied point's accepted geometry; no history is reset.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveSite<'a> {
    point: SourceBoundObjectiveCohesivePoint<'a>,
    arms_body_m: [Vec3; 2],
    frame_in_a: UnitQuaternion,
}
impl<'a> RigidCohesiveSite<'a> {
    pub fn new(
        point: SourceBoundObjectiveCohesivePoint<'a>,
        arms_body_m: [Vec3; 2],
        frame_in_a: UnitQuaternion,
    ) -> Result<Self, RigidCohesiveError> {
        for arm in arms_body_m {
            finite_vector(arm)?;
        }
        Ok(Self {
            point,
            arms_body_m,
            frame_in_a,
        })
    }
    pub const fn point(&self) -> SourceBoundObjectiveCohesivePoint<'a> {
        self.point
    }
    pub const fn arms_body_m(&self) -> [Vec3; 2] {
        self.arms_body_m
    }
    pub const fn frame_in_a(&self) -> UnitQuaternion {
        self.frame_in_a
    }
}

/// Explicit physical domain, solve work and independent balance gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveConfig {
    /// At most 64 original quadrature points; not 64 independent mechanics.
    pub maximum_sites: usize,
    pub maximum_iterations: usize,
    /// Complete two-body/patch evaluations, including all Jacobian and line-search calls.
    pub maximum_evaluations: usize,
    /// Each original material trial consumes one, including rejected candidates.
    pub maximum_material_evaluations: usize,
    pub maximum_backtracks: usize,
    pub maximum_step_s: f64,
    pub maximum_site_force_n: f64,
    pub maximum_total_force_n: f64,
    pub maximum_total_torque_n_m: f64,
    pub maximum_jump_increment_m: f64,
    pub force_absolute_tolerance_n: f64,
    pub torque_absolute_tolerance_n_m: f64,
    pub relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
    pub linear_momentum_absolute_tolerance_n_s: f64,
    pub angular_momentum_absolute_tolerance_n_m_s: f64,
    pub momentum_relative_tolerance: f64,
    /// Initial frame agreement and native Cayley-kinematics gate [dimensionless].
    pub frame_tolerance: f64,
    /// Original constitutive geometry, work and finite-rotation chart bounds.
    pub objective: ObjectiveCohesiveConfig,
}

/// World force and a spatial moment conjugate to the finite Cayley increment.
///
/// They are held over one step. The force is applied at the centre of mass;
/// the moment is an independent couple. Actual angular impulse is dt * moment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RigidCohesiveLoad {
    pub force_world_n: Vec3,
    pub moment_world_n_m: Vec3,
}
impl RigidCohesiveLoad {
    pub const ZERO: Self = Self {
        force_world_n: Vec3::ZERO,
        moment_world_n_m: Vec3::ZERO,
    };
}

/// Accepted per-site physics; subsequent failed trials cannot replace it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveSiteFrame {
    pub jump_m: [f64; 3],
    pub damage: f64,
    pub discrete_force_on_a_n: Vec3,
    pub frame_couple_on_a_n_m: Vec3,
    pub stored_energy_j: f64,
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub spatial_work_j: f64,
}

/// Independent work, storage and momentum checks on actual native successors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidCohesiveFrame {
    pub step: u64,
    pub time_s: f64,
    /// Kinetic plus the original uniform-gravity potential of both bodies.
    pub mechanical_energy_j: f64,
    pub cohesive_stored_energy_j: f64,
    pub external_work_j: f64,
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub cohesive_spatial_work_j: f64,
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
    pub linear_momentum_residual_n_s: Vec3,
    pub linear_momentum_tolerance_n_s: f64,
    pub angular_momentum_residual_n_m_s: Vec3,
    pub angular_momentum_tolerance_n_m_s: f64,
    pub contact_loads: [RigidCohesiveLoad; 2],
    pub native_wrenches: [Wrench; 2],
    pub maximum_force_residual_n: f64,
    pub maximum_torque_residual_n_m: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub material_evaluations: usize,
}

/// Refusals leave every accepted body, bond, observation and clock unchanged.
#[derive(Debug, Clone, PartialEq)]
pub enum RigidCohesiveError {
    Invalid(&'static str),
    Dynamics(DynamicsError),
    Source(SourceBoundObjectiveCohesiveError),
    Objective(ObjectiveCohesiveError),
    Factor(LuWorkspaceError),
    Cancelled,
    Budget(&'static str),
    Range {
        what: &'static str,
        value: f64,
        maximum: f64,
    },
    Balance {
        what: &'static str,
        residual: f64,
        tolerance: f64,
    },
    NotConverged {
        scaled_residual: f64,
        iterations: usize,
    },
}
impl core::fmt::Display for RigidCohesiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rigid cohesive trial refused: {self:?}")
    }
}
impl std::error::Error for RigidCohesiveError {}
impl From<DynamicsError> for RigidCohesiveError {
    fn from(value: DynamicsError) -> Self {
        Self::Dynamics(value)
    }
}
impl From<SourceBoundObjectiveCohesiveError> for RigidCohesiveError {
    fn from(value: SourceBoundObjectiveCohesiveError) -> Self {
        Self::Source(value)
    }
}
impl From<ObjectiveCohesiveError> for RigidCohesiveError {
    fn from(value: ObjectiveCohesiveError) -> Self {
        Self::Objective(value)
    }
}
impl From<LuWorkspaceError> for RigidCohesiveError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::Factor(value)
    }
}

/// One accepted physical clock and reusable scratch for the complete pair.
pub struct RigidCohesiveSystem<'a> {
    bodies: [RigidCohesiveBody; 2],
    sites: Vec<RigidCohesiveSite<'a>>,
    integrator: RigidBodyIntegrator,
    dt_s: f64,
    config: RigidCohesiveConfig,
    steps: u64,
    time_s: f64,
    linear_solve: LuWorkspace,
    staged_sites: Vec<Option<SourceBoundObjectiveCohesiveStep<'a>>>,
    site_frames: Vec<Option<RigidCohesiveSiteFrame>>,
    last_wrenches: [Wrench; 2],
    last_frame: Option<RigidCohesiveFrame>,
}
impl<'a> RigidCohesiveSystem<'a> {
    pub fn new(
        bodies: [RigidCohesiveBody; 2],
        sites: Vec<RigidCohesiveSite<'a>>,
        gravity: Gravity,
        dt_s: f64,
        config: RigidCohesiveConfig,
        gate: &CancelGate,
    ) -> Result<Self, RigidCohesiveError> {
        poll(Some(gate))?;
        validate_config(config)?;
        positive(dt_s)?;
        bound("step duration", dt_s, config.maximum_step_s)?;
        if sites.is_empty() || sites.len() > config.maximum_sites {
            return Err(RigidCohesiveError::Budget("finite cohesive sites"));
        }
        let integrator = RigidBodyIntegrator::new(gravity);
        for body in &bodies {
            body.check(body.state)?;
            integrator.diagnostics(body.state, body.properties)?;
        }
        for site in &sites {
            poll(Some(gate))?;
            site.point
                .source()
                .require_surface_materials(&bodies[0].material, &bodies[1].material)
                .map_err(SourceBoundObjectiveCohesiveError::MaterialAssignment)?;
            let frame = compose(bodies[0].state.pose().orientation(), site.frame_in_a)?;
            let frame_error = cayley_increment(
                frame,
                site.point.numerical().frame(),
                config.objective.maximum_cayley_increment,
            )?;
            balance(
                "initial material frame",
                norm(frame_error),
                config.frame_tolerance,
            )?;
            let separation = separation(&bodies.each_ref().map(|b| b.state), site.arms_body_m)?;
            balance(
                "initial attachment separation",
                norm(separation.sub(site.point.numerical().separation_world_m())),
                config.objective.kinematic_tolerance_m,
            )?;
        }
        let mut staged_sites = Vec::new();
        staged_sites
            .try_reserve_exact(sites.len())
            .map_err(|_| RigidCohesiveError::Budget("site workspace allocation"))?;
        staged_sites.resize_with(sites.len(), || None);
        let mut site_frames = Vec::new();
        site_frames
            .try_reserve_exact(sites.len())
            .map_err(|_| RigidCohesiveError::Budget("site observation allocation"))?;
        site_frames.resize(sites.len(), None);
        let system = Self {
            bodies,
            sites,
            integrator,
            dt_s,
            config,
            steps: 0,
            time_s: 0.0,
            linear_solve: LuWorkspace::new(12)?,
            staged_sites,
            site_frames,
            last_wrenches: [Wrench::ZERO; 2],
            last_frame: None,
        };
        finite(system.total_energy_j()?)?;
        poll(Some(gate))?;
        Ok(system)
    }

    pub const fn bodies(&self) -> &[RigidCohesiveBody; 2] {
        &self.bodies
    }
    pub fn sites(&self) -> &[RigidCohesiveSite<'a>] {
        &self.sites
    }
    pub fn site_frames(&self) -> &[Option<RigidCohesiveSiteFrame>] {
        &self.site_frames
    }
    pub const fn config(&self) -> RigidCohesiveConfig {
        self.config
    }
    pub const fn step_duration_s(&self) -> f64 {
        self.dt_s
    }
    pub const fn steps(&self) -> u64 {
        self.steps
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn last_frame(&self) -> Option<&RigidCohesiveFrame> {
        self.last_frame.as_ref()
    }
    pub const fn gravity(&self) -> Gravity {
        self.integrator.gravity()
    }

    pub fn total_energy_j(&self) -> Result<f64, RigidCohesiveError> {
        let mut total = 0.0;
        for body in &self.bodies {
            total += self
                .integrator
                .diagnostics(body.state, body.properties)?
                .mechanical_energy;
        }
        for site in &self.sites {
            total += site.point.numerical().area_m2()
                * site.point.numerical().material().stored_energy_j_m2();
        }
        finite(total)
    }

    pub fn step(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_inner(loads, None)
    }
    pub fn step_under_gate(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        gate: &CancelGate,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        self.step_inner(loads, Some(gate))
    }

    fn step_inner(
        &mut self,
        loads: [RigidCohesiveLoad; 2],
        gate: Option<&CancelGate>,
    ) -> Result<&RigidCohesiveFrame, RigidCohesiveError> {
        poll(gate)?;
        for load in loads {
            finite_vector(load.force_world_n)?;
            finite_vector(load.moment_world_n_m)?;
        }
        let next_step = self
            .steps
            .checked_add(1)
            .ok_or(RigidCohesiveError::Budget("clock"))?;
        let next_time = finite(next_step as f64 * self.dt_s)?;
        if next_time <= self.time_s {
            return Err(RigidCohesiveError::Budget("clock resolution"));
        }
        let (candidate, wrenches, iterations, work) = self.solve(loads, gate)?;
        let frame = self.audit(
            candidate, wrenches, loads, next_step, next_time, iterations, work, gate,
        )?;
        if self.staged_sites.iter().any(Option::is_none) {
            return Err(RigidCohesiveError::Invalid("incomplete staged patch"));
        }
        poll(gate)?;
        // All fallible native, source, geometry and work checks precede publication.
        for body in 0..2 {
            self.bodies[body].state = candidate.receipts[body].state_after;
        }
        for (index, site) in self.sites.iter_mut().enumerate() {
            let step = self.staged_sites[index]
                .as_ref()
                .expect("complete staged patch checked");
            let numerical = step.numerical();
            site.point = step.point();
            self.site_frames[index] = Some(RigidCohesiveSiteFrame {
                jump_m: numerical.point.material().jump_m(),
                damage: numerical.point.material().damage(),
                discrete_force_on_a_n: numerical.discrete_wrench.force_on_a_n,
                frame_couple_on_a_n_m: numerical.discrete_wrench.frame_couple_n_m,
                stored_energy_j: numerical.point.area_m2()
                    * numerical.point.material().stored_energy_j_m2(),
                fracture_dissipation_j: numerical.point.area_m2()
                    * numerical.material.equilibrium.energy.dissipated,
                viscous_dissipation_j: numerical.point.area_m2()
                    * numerical.material.viscous_dissipated_j_m2,
                spatial_work_j: numerical.spatial_work_j,
            });
        }
        self.steps = next_step;
        self.time_s = next_time;
        self.last_wrenches = wrenches;
        self.last_frame = Some(frame);
        Ok(self.last_frame.as_ref().expect("frame just published"))
    }
}

fn validate_config(config: RigidCohesiveConfig) -> Result<(), RigidCohesiveError> {
    if !(1..=64).contains(&config.maximum_sites)
        || !(1..=128).contains(&config.maximum_iterations)
        || config.maximum_evaluations == 0
        || config.maximum_material_evaluations == 0
        || !(1..=32).contains(&config.maximum_backtracks)
    {
        return Err(RigidCohesiveError::Invalid(
            "bounded site and nonlinear work budgets required",
        ));
    }
    for value in [
        config.maximum_step_s,
        config.maximum_site_force_n,
        config.maximum_total_force_n,
        config.maximum_total_torque_n_m,
        config.maximum_jump_increment_m,
        config.force_absolute_tolerance_n,
        config.torque_absolute_tolerance_n_m,
        config.energy_absolute_tolerance_j,
        config.linear_momentum_absolute_tolerance_n_s,
        config.angular_momentum_absolute_tolerance_n_m_s,
        config.frame_tolerance,
        config.objective.maximum_cayley_increment,
        config.objective.kinematic_tolerance_m,
        config.objective.material_energy_tolerance_j_m2,
        config.objective.spatial_work_tolerance_j,
    ] {
        positive(value)?;
    }
    for value in [
        config.relative_tolerance,
        config.energy_relative_tolerance,
        config.momentum_relative_tolerance,
        config.frame_tolerance,
    ] {
        if !value.is_finite() || value <= 0.0 || value >= 1.0 {
            return Err(RigidCohesiveError::Invalid(
                "relative and frame tolerances must be in (0,1)",
            ));
        }
    }
    Ok(())
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> Result<UnitQuaternion, DynamicsError> {
    let [aw, ax, ay, az] = a.components();
    let [bw, bx, by, bz] = b.components();
    UnitQuaternion::new(
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    )
}
fn separation(states: &[RigidBodyState; 2], arms: [Vec3; 2]) -> Result<Vec3, DynamicsError> {
    Ok(states[1]
        .pose()
        .point_world_from_body(arms[1])?
        .sub(states[0].pose().point_world_from_body(arms[0])?))
}
fn finite(x: f64) -> Result<f64, RigidCohesiveError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RigidCohesiveError::Invalid("nonfinite derived value"))
    }
}
fn finite_vector(x: Vec3) -> Result<(), RigidCohesiveError> {
    if x.is_finite() {
        Ok(())
    } else {
        Err(RigidCohesiveError::Invalid("nonfinite vector"))
    }
}
fn positive(x: f64) -> Result<f64, RigidCohesiveError> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(RigidCohesiveError::Invalid(
            "positive finite physical bounds required",
        ))
    }
}
fn norm(x: Vec3) -> f64 {
    x.x.hypot(x.y).hypot(x.z)
}
fn bound(what: &'static str, value: f64, maximum: f64) -> Result<(), RigidCohesiveError> {
    finite(value)?;
    if value > maximum {
        Err(RigidCohesiveError::Range {
            what,
            value,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn balance(what: &'static str, residual: f64, tolerance: f64) -> Result<(), RigidCohesiveError> {
    finite(residual)?;
    positive(tolerance)?;
    if residual.abs() > tolerance {
        Err(RigidCohesiveError::Balance {
            what,
            residual,
            tolerance,
        })
    } else {
        Ok(())
    }
}
fn poll(gate: Option<&CancelGate>) -> Result<(), RigidCohesiveError> {
    if gate.is_some_and(CancelGate::is_requested) {
        Err(RigidCohesiveError::Cancelled)
    } else {
        Ok(())
    }
}

//! Implicit pressure/viscous loading of a guided or planar rigid body.
//!
//! [`planar`] extends the same native coupling to both translations and rotation.
//! The guided owner below retains its scalar, three-solve path.
//!
//! The original material BDF1 Stokes solver supplies the affine load response
//! `F(v) = F(0) + Z v`. Its scalar Schur complement closes body inertia, a linear
//! spring and a dashpot at the SAME endpoint velocity. A final native flow solve
//! checks that endpoint before publishing fluid history, rigid momentum, pose
//! and time together. No empirical added mass or drag coefficient is introduced.
//!
//! This is a fixed-reference SMALL-DISPLACEMENT model. Boundary velocities act
//! on the original mesh; its geometry does not follow the reported body pose.
//! The caller declares a displacement bound appropriate to that approximation,
//! which every accepted endpoint must satisfy. There is no ALE, transported
//! fluid volume, convection, thermal pressure work or finite-motion claim.
//! Incompressible pressure gauges remain explicit. The source material and its
//! frozen Newtonian transport state are retained by the original binding.
//! Explicit [`feedback`] steps instead sample source viscosity from the original
//! thermal cells for each interval while retaining that fixed density and mesh.
//!
//! Three native solves per step (zero motion, homogeneous unit motion, endpoint)
//! avoid a partitioned added-mass fixed-point iteration. Source and exterior
//! boundary callbacks must be pure and constant during those solves. Cancellation
//! is polled between solves and before publication; the synchronous flow solver
//! itself has no new cancellation guarantee.

use fs_flux::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxError, MaterialFluxNumerics, MaterialFluxSystem,
    MaterialStokesBoundary, MaterialStokesSolution,
};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};

pub mod planar;

/// Physical bulk viscous heat on the original finite-cell thermal timeline.
pub mod heating;

/// Native heat powers drive source-bound instantaneous liquid expansion.
pub mod thermal_dilation;

/// Source viscosity sampled from the original finite thermal cells each interval.
pub mod feedback;

/// Finite guided translation on the original nonlinear moving-domain flow.
pub mod translating;

/// Finite planar translation and rotation through the original rigid ALE flow.
pub mod rigid;

mod search;

/// Fixed physical and approximation data; all mechanical values cover the
/// explicit out-of-plane depth, rather than being per-unit-depth quantities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GuidedStokesConfig {
    /// Unit direction in the mesh's physical x/y frame.
    pub direction: [f64; 2],
    /// Extrusion depth [m] converting original 2-D stresses into forces [N].
    pub depth_m: f64,
    /// Maximum absolute translation from the construction pose [m]. This must
    /// be small compared with the caller's relevant gaps/geometry lengths.
    pub max_displacement_m: f64,
    /// Spring stiffness [N/m], with rest position at the construction pose.
    pub stiffness_n_m: f64,
    /// Mechanical dashpot coefficient [N s/m]; excludes resolved fluid stress.
    pub damping_n_s_m: f64,
    /// Allowed mismatch of strong normal velocity moments [m/s].
    pub velocity_tolerance_m_s: f64,
}

/// Explicit acceptance tolerances for the coupled endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GuidedStokesControls {
    /// Absolute residual in the guided body momentum equation [N].
    pub force_absolute_n: f64,
    /// Relative residual, scaled by the largest term in that equation.
    pub force_relative: f64,
    /// Absolute residual of the body's discrete mechanical-work identity [J].
    pub work_absolute_j: f64,
    /// Relative residual, scaled by the largest work/storage term.
    pub work_relative: f64,
}

impl Default for GuidedStokesControls {
    fn default() -> Self {
        Self {
            force_absolute_n: 1e-9,
            force_relative: 1e-8,
            work_absolute_j: 1e-10,
            work_relative: 1e-8,
        }
    }
}

/// All mutable physical history, committed only after a complete accepted step.
#[derive(Clone, Debug, PartialEq)]
pub struct GuidedStokesState {
    body: RigidBodyState,
    velocity: Vec<f64>,
    displacement_m: f64,
    time_s: f64,
}

impl GuidedStokesState {
    /// Original rigid-body momentum and accepted small-displacement pose.
    #[must_use]
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Accepted BDM1 velocity moments, in m/s.
    #[must_use]
    pub fn fluid_velocity(&self) -> &[f64] {
        &self.velocity
    }
    /// Guided translation relative to the reference pose [m].
    #[must_use]
    pub const fn displacement_m(&self) -> f64 {
        self.displacement_m
    }
    /// Physical accepted time [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
}

/// Rejection leaves the accepted physical state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum GuidedStokesError {
    /// Invalid dimensions, units, guide, boundary ownership or arithmetic.
    Invalid(&'static str),
    /// Original material/flow admission or solver failure.
    Flow(MaterialFluxError),
    /// Original rigid-body validation or impulse failure.
    Mechanics(DynamicsError),
    /// A strong normal interface velocity does not match the guided body.
    InterfaceVelocity { edge: usize, residual_m_s: f64 },
    /// The proposed endpoint is outside the declared linearized geometry domain.
    DisplacementLimit { proposed_m: f64, limit_m: f64 },
    /// The scalar coupled momentum operator is nonpositive or unrepresentable.
    SingularCoupling { effective_coefficient_n_s_m: f64 },
    /// Actual final fluid/body coupling fails its force or mechanical-work gate.
    Balance {
        force_residual_n: f64,
        mechanical_work_residual_j: f64,
    },
    /// Cancellation at a composition boundary, before any physical publication.
    Cancelled,
}

impl core::fmt::Display for GuidedStokesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "guided Stokes motion refused: {self:?}")
    }
}
impl std::error::Error for GuidedStokesError {}
impl From<MaterialFluxError> for GuidedStokesError {
    fn from(value: MaterialFluxError) -> Self {
        Self::Flow(value)
    }
}
impl From<DynamicsError> for GuidedStokesError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}

/// Accepted original flow and explicit mechanical/interface accounts.
pub struct GuidedStokesStep {
    /// Native material flow at exactly the accepted body endpoint speed.
    pub flow: MaterialStokesSolution,
    /// Original pressure and viscous reaction on selected body edges, still per
    /// metre depth, about the fixed reference COM. Advective transport is absent.
    pub boundary_load: MaterialBoundaryLoad,
    /// Original fixed-pose momentum update. `body_after` includes the subsequent
    /// backward-Euler drift; no angular impulse or rotational drift is applied.
    pub momentum: MomentumImpulseReceipt,
    /// Accepted translated body state.
    pub body_after: RigidBodyState,
    /// New time [s].
    pub time_s: f64,
    /// Endpoint speed along the guide [m/s].
    pub speed_m_s: f64,
    /// Actual generalized pressure/viscous force on the body [N].
    pub fluid_force_n: f64,
    /// Homogeneous derivative dF_fluid/dv [N s/m], computed from the native flow.
    pub fluid_force_slope_n_s_m: f64,
    /// Force exerted by the ideal guide to cancel the transverse fluid load [N].
    pub guide_force_n: [f64; 2],
    /// Guide torque canceling the fluid moment about the reference COM [N m].
    pub guide_torque_n_m: f64,
    /// Fluid-on-body work `dt F_fluid v_new` [J]. The opposite sign is the
    /// conjugate interface transfer. This is not a whole-fluid energy audit.
    pub interface_work_j: f64,
    /// Work by the authored endpoint external force [J].
    pub external_work_j: f64,
    /// Change in original body kinetic energy [J].
    pub body_kinetic_change_j: f64,
    /// Change in spring storage [J].
    pub spring_energy_change_j: f64,
    /// Physical mechanical dashpot loss [J], excluding fluid viscous dissipation.
    pub damping_loss_j: f64,
    /// Backward-Euler body inertia and spring loss [J]; numerical, not heat.
    pub mechanical_time_discretization_loss_j: f64,
    /// Change in original fluid kinetic energy over the declared depth [J].
    /// This remains separate from physical and numerical dissipation.
    pub fluid_kinetic_change_j: f64,
    /// Bulk fluid viscous heat [J], using BDF1 endpoint quadrature
    /// `dt * depth * integral(2 mu D:D)`. This excludes wall slip and numerical
    /// losses, is not an exact trajectory integral, and updates no thermal state.
    pub bulk_viscous_dissipation_j: f64,
    /// Actual residual of the endpoint guided momentum equation [N].
    pub force_residual_n: f64,
    /// Body storage plus physical/numerical losses minus applied/interface work.
    pub mechanical_work_residual_j: f64,
}

/// The original fixed mesh and guided rigid body share one accepted history.
/// Ordinary methods retain the frozen construction material. Explicit thermal
/// feedback methods resolve a checked source for their requested interval,
/// preserving the original density, mesh and mechanical reference pose.
pub struct GuidedStokesMotion<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'m>,
    fluid: MaterialFluxSystem<'s, 'm, 'a, F>,
    properties: MassProperties,
    reference_pose: Pose,
    config: GuidedStokesConfig,
    edges: Vec<usize>,
    state: GuidedStokesState,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> GuidedStokesMotion<'s, 'm, 'a, F> {
    /// Bind existing native physical owners. Body edges are the complete fluid
    /// interface assigned to this one body; remaining edges are supplied at each
    /// step. Initial normal moments must agree with its initial guided velocity.
    /// The guide admits translation only; initial angular and transverse
    /// momentum must vanish (the transverse test permits arithmetic roundoff).
    ///
    /// # Errors
    /// Refuses invalid physical inputs, empty/duplicate/nonboundary selections,
    /// off-guide initial motion, invalid fluid history or unavailable transport.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_fluid_velocity: Vec<f64>,
        body_edges: &[usize],
        config: GuidedStokesConfig,
    ) -> Result<Self, GuidedStokesError> {
        let norm = config.direction[0].hypot(config.direction[1]);
        if !norm.is_finite()
            || (norm - 1.0).abs() > 64.0 * f64::EPSILON
            || !positive(config.depth_m)
            || !positive(config.max_displacement_m)
            || !nonnegative(config.stiffness_n_m)
            || !nonnegative(config.damping_n_s_m)
            || !positive(config.velocity_tolerance_m_s)
            || initial_fluid_velocity.len() != system.n_u
            || !initial_fluid_velocity.iter().all(|v| v.is_finite())
            || body_edges.is_empty()
        {
            return Err(GuidedStokesError::Invalid(
                "finite physical data, unit guide and nonempty matching history/interface required",
            ));
        }
        let mut edges = body_edges.to_vec();
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1])
            || edges
                .iter()
                .any(|&edge| edge >= system.mesh.edges.len() || !system.mesh.is_boundary(edge))
        {
            return Err(GuidedStokesError::Invalid(
                "body interface must contain unique boundary edges",
            ));
        }
        let axis = Vec3::new(config.direction[0], config.direction[1], 0.0);
        let momentum = body.linear_momentum_world();
        let transverse = momentum.sub(axis.scale(momentum.dot(axis)));
        let momentum_norm = momentum.x.hypot(momentum.y).hypot(momentum.z);
        let transverse_norm = transverse.x.hypot(transverse.y).hypot(transverse.z);
        if body.angular_momentum_body() != Vec3::ZERO
            || !transverse_norm.is_finite()
            || transverse_norm > 128.0 * f64::EPSILON * momentum_norm.max(f64::MIN_POSITIVE)
        {
            return Err(GuidedStokesError::Invalid(
                "initial body momentum must follow the translation guide without rotation",
            ));
        }
        let fluid = MaterialFluxSystem::new(system, material, numerics)?;
        let result = Self {
            system,
            fluid,
            properties,
            reference_pose: body.pose(),
            config,
            edges,
            state: GuidedStokesState {
                body,
                velocity: initial_fluid_velocity,
                displacement_m: 0.0,
                time_s: 0.0,
            },
        };
        result.check_interface(&result.state.velocity, result.speed(body)?)?;
        Ok(result)
    }

    /// Accepted physical state; failed steps never mutate it.
    #[must_use]
    pub const fn state(&self) -> &GuidedStokesState {
        &self.state
    }
    /// The original source-bearing material used by ordinary frozen-material steps.
    /// Temperature-feedback receipts retain the source resolved for their interval.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.fluid.material()
    }

    fn axis(&self) -> Vec3 {
        Vec3::new(self.config.direction[0], self.config.direction[1], 0.0)
    }
    fn speed(&self, body: RigidBodyState) -> Result<f64, GuidedStokesError> {
        Ok(body
            .center_of_mass_velocity_world(self.properties)?
            .dot(self.axis()))
    }
    fn check_interface(&self, velocity: &[f64], speed: f64) -> Result<(), GuidedStokesError> {
        for &edge in &self.edges {
            let normal = self.system.mesh.edges[edge].normal;
            let expected = speed
                * (normal[0] * self.config.direction[0] + normal[1] * self.config.direction[1]);
            let residual = (velocity[2 * edge] - expected)
                .abs()
                .max(velocity[2 * edge + 1].abs());
            if !residual.is_finite() || residual > self.config.velocity_tolerance_m_s {
                return Err(GuidedStokesError::InterfaceVelocity {
                    edge,
                    residual_m_s: residual,
                });
            }
        }
        Ok(())
    }
    fn load(
        &self,
        flow: &MaterialStokesSolution,
        offsets: &[Option<f64>],
    ) -> Result<(MaterialBoundaryLoad, f64), GuidedStokesError> {
        let reference = self.reference_pose.position_world();
        let load = flow.boundary_loads(&self.edges, [reference.x, reference.y], offsets)?;
        let force = self.config.depth_m
            * (load.stress.force_n_per_m[0] * self.config.direction[0]
                + load.stress.force_n_per_m[1] * self.config.direction[1]);
        if !force.is_finite() {
            return Err(GuidedStokesError::Invalid("physical fluid force overflow"));
        }
        Ok((load, force))
    }
    fn solve(
        &self,
        speed: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        previous: &[f64],
        dt: f64,
    ) -> Result<MaterialStokesSolution, GuidedStokesError> {
        let motion = |_: [f64; 2]| self.config.direction.map(|d| d * speed);
        Ok(solve_body_velocity(
            &self.fluid,
            &self.edges,
            &motion,
            force_density,
            exterior,
            previous,
            dt,
        )?)
    }

    /// Close one implicit fluid/body step with endpoint force and boundary data.
    /// `external_force_n` acts along the guide and must include any exterior
    /// pressure on the other side of the body; that pressure is never inferred.
    /// `exterior` owns every mesh boundary edge not assigned to this body.
    /// Mean-zero fluid components require explicit pressure offsets [Pa], while
    /// traction-datum components require None, exactly as in the native load API.
    /// Callbacks must be pure and fixed across all three linear solves.
    ///
    /// The mechanical work gate is a body/spring/dashpot BDF1 identity. Fluid
    /// kinetic change is separately reported; neither this identity nor opposite
    /// interface work signs prove global fluid/body passivity. Weak tangential
    /// SIP slip and the original flow's viscous/numerical terms remain its own.
    ///
    /// # Errors
    /// All original solver failures, cancellation, invalid timestep/controls,
    /// singular coupling, displacement-domain escape or failed endpoint balances
    /// preserve the entire accepted state.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesStep, GuidedStokesError> {
        let (state, step) = self.propose_step(
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            controls,
            cancelled,
        )?;
        self.state = state;
        Ok(step)
    }

    // Both ordinary stepping and thermal composition use this same native
    // proposal. Neither fluid nor body history changes before acceptance.
    #[allow(clippy::too_many_arguments)]
    fn propose_step(
        &self,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(GuidedStokesState, GuidedStokesStep), GuidedStokesError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        let time = self.state.time_s + dt;
        if !positive(dt)
            || !positive(dt.recip())
            || !time.is_finite()
            || time <= self.state.time_s
            || !external_force_n.is_finite()
            || ![
                controls.force_absolute_n,
                controls.force_relative,
                controls.work_absolute_j,
                controls.work_relative,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(GuidedStokesError::Invalid(
                "representable positive timestep and tolerances, finite applied force required",
            ));
        }
        let old_speed = self.speed(self.state.body)?;
        let q0 = self.state.displacement_m;
        let base = self.solve(0.0, force_density, exterior, &self.state.velocity, dt)?;
        let (_, base_force) = self.load(&base, pressure_offsets_pa)?;
        drop(base);
        poll(&mut cancelled)?;
        // Solve the homogeneous response directly. Subtracting two loads with
        // large absolute pressure would lose the small dynamic impedance.
        let zero = |_: [f64; 2]| [0.0; 2];
        let zero_traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
        let homogeneous: Vec<_> = exterior
            .iter()
            .map(|boundary| match *boundary {
                MaterialStokesBoundary::Velocity { edge, .. } => {
                    MaterialStokesBoundary::Velocity { edge, value: &zero }
                }
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    tangential_resistance_pa_s_per_m,
                    ..
                } => {
                    // Differentiate wall data while retaining the physical
                    // Robin operator in the native homogeneous response.
                    MaterialStokesBoundary::NavierSlip {
                        edge,
                        value: &zero,
                        tangential_resistance_pa_s_per_m,
                    }
                }
                MaterialStokesBoundary::TotalTraction { edge, .. } => {
                    MaterialStokesBoundary::TotalTraction {
                        edge,
                        value: &zero_traction,
                    }
                }
            })
            .collect();
        let zero_history = vec![0.0; self.system.n_u];
        let response = self.solve(1.0, &zero, &homogeneous, &zero_history, dt)?;
        let response_offsets: Vec<_> = pressure_offsets_pa.iter().map(|p| p.map(|_| 0.0)).collect();
        let (_, slope) = self.load(&response, &response_offsets)?;
        drop(response);
        poll(&mut cancelled)?;
        let mass = self.properties.mass();
        let stiffness = self.config.stiffness_n_m;
        let damping = self.config.damping_n_s_m;
        let coefficient = mass / dt + damping + stiffness * dt - slope;
        if !positive(coefficient) {
            return Err(GuidedStokesError::SingularCoupling {
                effective_coefficient_n_s_m: coefficient,
            });
        }
        let proposed_speed =
            (mass / dt * old_speed + external_force_n + base_force - stiffness * q0) / coefficient;
        if !proposed_speed.is_finite() {
            return Err(GuidedStokesError::Invalid("unrepresentable coupled speed"));
        }
        // Use the original rigid impulse owner, then solve the final fluid at
        // its actual representable speed, not at an uncommitted scalar guess.
        let momentum = self.state.body.apply_momentum_impulse(
            self.properties,
            self.axis().scale(mass * (proposed_speed - old_speed)),
            Vec3::ZERO,
        )?;
        let speed = self.speed(momentum.state_after)?;
        let displacement = dt * speed;
        let q1 = q0 + displacement;
        if !q1.is_finite() || q1.abs() > self.config.max_displacement_m {
            return Err(GuidedStokesError::DisplacementLimit {
                proposed_m: q1,
                limit_m: self.config.max_displacement_m,
            });
        }
        let reference = self.reference_pose;
        let body = RigidBodyState::new(
            Pose::new(
                reference.position_world().add(self.axis().scale(q1)),
                reference.orientation(),
            )?,
            momentum.state_after.linear_momentum_world(),
            Vec3::ZERO,
        )?;
        let actual_displacement = q1 - q0;
        let actual_drift = body
            .pose()
            .position_world()
            .sub(self.state.body.pose().position_world());
        let expected_drift = self.axis().scale(displacement);
        let drift_error = actual_drift.sub(expected_drift);
        let position_tolerance = self.config.velocity_tolerance_m_s * dt;
        if (actual_displacement - displacement).abs() > position_tolerance
            || drift_error.x.hypot(drift_error.y).hypot(drift_error.z) > position_tolerance
        {
            return Err(GuidedStokesError::Invalid(
                "body displacement is not representable at the declared velocity tolerance",
            ));
        }
        let flow = self.solve(speed, force_density, exterior, &self.state.velocity, dt)?;
        let (load, force) = self.load(&flow, pressure_offsets_pa)?;
        self.check_interface(flow.velocity_dofs(), speed)?;
        poll(&mut cancelled)?;
        let dv = speed - old_speed;
        // Velocity division can round independently of the stored momentum.
        // Audit the native increment, including any rounded or lost small kick.
        let inertia_force = momentum.linear_increment_world.dot(self.axis()) / dt;
        let spring_force = stiffness * q1;
        let damping_force = damping * speed;
        let force_residual =
            inertia_force - (force + external_force_n - spring_force - damping_force);
        let interface_work = force * displacement;
        let external_work = external_force_n * displacement;
        let kinetic_change = momentum.work.kinetic_energy_change_j;
        let spring_change = 0.5 * stiffness * (q1 + q0) * actual_displacement;
        let damping_loss = damping_force * displacement;
        let temporal_loss =
            0.5 * mass * dv * dv + 0.5 * stiffness * actual_displacement * actual_displacement;
        let work_residual = kinetic_change + spring_change + damping_loss + temporal_loss
            - interface_work
            - external_work;
        let force_scale = [
            inertia_force,
            force,
            external_force_n,
            spring_force,
            damping_force,
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
        let work_scale = [
            kinetic_change,
            spring_change,
            damping_loss,
            temporal_loss,
            interface_work,
            external_work,
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
        if !force_residual.is_finite()
            || !work_residual.is_finite()
            || force_residual.abs()
                > controls.force_absolute_n + controls.force_relative * force_scale
            || work_residual.abs() > controls.work_absolute_j + controls.work_relative * work_scale
        {
            return Err(GuidedStokesError::Balance {
                force_residual_n: force_residual,
                mechanical_work_residual_j: work_residual,
            });
        }
        let full_force = load.stress.force_n_per_m.map(|v| v * self.config.depth_m);
        let guide_force =
            core::array::from_fn(|axis| force * self.config.direction[axis] - full_force[axis]);
        let guide_torque = -load.stress.moment_n * self.config.depth_m;
        let fluid_energy_change = self.config.depth_m
            * (flow.kinetic_energy_j_per_m()
                - self.fluid.transport().density_kg_m3()
                    * self.system.kinetic_energy(&self.state.velocity));
        let viscous_heat =
            dt * self.config.depth_m * flow.viscous_dissipation().total_power_w_per_m();
        if !guide_force.iter().all(|v| v.is_finite())
            || !guide_torque.is_finite()
            || !fluid_energy_change.is_finite()
            || !viscous_heat.is_finite()
        {
            return Err(GuidedStokesError::Invalid(
                "physical effort or energy overflow",
            ));
        }
        let accepted_velocity = flow.velocity_dofs().to_vec();
        poll(&mut cancelled)?;
        let state = GuidedStokesState {
            body,
            velocity: accepted_velocity,
            displacement_m: q1,
            time_s: time,
        };
        Ok((
            state,
            GuidedStokesStep {
                flow,
                boundary_load: load,
                momentum,
                body_after: body,
                time_s: time,
                speed_m_s: speed,
                fluid_force_n: force,
                fluid_force_slope_n_s_m: slope,
                guide_force_n: guide_force,
                guide_torque_n_m: guide_torque,
                interface_work_j: interface_work,
                external_work_j: external_work,
                body_kinetic_change_j: kinetic_change,
                spring_energy_change_j: spring_change,
                damping_loss_j: damping_loss,
                mechanical_time_discretization_loss_j: temporal_loss,
                fluid_kinetic_change_j: fluid_energy_change,
                bulk_viscous_dissipation_j: viscous_heat,
                force_residual_n: force_residual,
                mechanical_work_residual_j: work_residual,
            },
        ))
    }
}

/// Keep both rigid lanes on the same original mixed-boundary flow entry.
fn solve_body_velocity<F: FluidTransport + ?Sized>(
    fluid: &MaterialFluxSystem<'_, '_, '_, F>,
    edges: &[usize],
    motion: &dyn Fn([f64; 2]) -> [f64; 2],
    force_density: &dyn Fn([f64; 2]) -> [f64; 2],
    exterior: &[MaterialStokesBoundary<'_>],
    previous: &[f64],
    dt: f64,
) -> Result<MaterialStokesSolution, MaterialFluxError> {
    let mut boundaries = Vec::with_capacity(exterior.len() + edges.len());
    boundaries.extend_from_slice(exterior);
    boundaries.extend(edges.iter().map(|&edge| MaterialStokesBoundary::Velocity {
        edge,
        value: motion,
    }));
    fluid.bdf1_stokes_mixed_step(force_density, &boundaries, previous, dt)
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GuidedStokesError> {
    if cancelled() {
        Err(GuidedStokesError::Cancelled)
    } else {
        Ok(())
    }
}

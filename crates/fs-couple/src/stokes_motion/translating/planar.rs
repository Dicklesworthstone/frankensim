//! Finite translation with TWO reciprocally coupled momentum components.
//!
//! The same nonlinear translating ALE owner supplies every trial stress load.
//! A bounded damped Newton search measures the full two-component response;
//! x/y equations are never solved as independent scalar guides. Orientation is
//! fixed, with an explicit zero-work restraint couple about the moving COM.
//! The entire domain translates: exterior callbacks live on its moving boundary.
//! This does not move exterior obstacles independently or rotate/deform the mesh.

use super::super::search::{EndpointSearch, Observation, SearchError};
use super::{impulse_time_loss, maximum, nonnegative, positive, sum};
use fs_flux::ale::translating::{
    TranslatingAleError, TranslatingAleFlow, TranslatingAleState, TranslatingAleStep,
};
use fs_flux::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxNumerics, MaterialNavierStokesBoundary,
    NavierStokesControls,
};
use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use std::cell::Cell;

/// Whole-body coefficients in fixed world x/y coordinates. Independent springs
/// rest at the construction pose; fluid cross-coordinate response is retained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarTranslationConfig {
    pub depth_m: f64,
    pub stiffness_n_m: [f64; 2],
    pub damping_n_s_m: [f64; 2],
    /// Body/mesh/drift and original unrelaxed normal-trace tolerance [m/s].
    pub velocity_tolerance_m_s: f64,
}

/// Authored finite search domain and physical acceptance, not a root-existence
/// guarantee. The initial trial is the old velocity clipped into this box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarTranslationControls {
    /// Ascending [minimum, maximum] endpoint velocities for each world axis.
    pub velocity_bounds_m_s: [[f64; 2]; 2],
    /// Positive one-sided response-probe sizes [m/s], shortened at a bound.
    /// Choose resolved physical scales relative to native flow accuracy.
    pub jacobian_step_m_s: [f64; 2],
    /// Includes base, derivative probes and rejected line-search flow solves.
    pub maximum_flow_solves: usize,
    pub force_absolute_n: f64,
    pub force_relative: f64,
    pub work_absolute_j: f64,
    pub work_relative: f64,
}

#[derive(Debug)]
pub enum PlanarTranslationError {
    Invalid(&'static str),
    Flow(TranslatingAleError),
    Mechanics(DynamicsError),
    LinearSolve(LuWorkspaceError),
    SolveLimit {
        flow_solves: usize,
        residual_n: [f64; 2],
    },
    /// Bounds, native velocity quantization or line search prevent progress.
    SearchStalled {
        residual_n: [f64; 2],
    },
    WorkBalance {
        residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for PlanarTranslationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar translating fluid/body coupling: {self:?}")
    }
}
impl std::error::Error for PlanarTranslationError {}
impl From<TranslatingAleError> for PlanarTranslationError {
    fn from(value: TranslatingAleError) -> Self {
        Self::Flow(value)
    }
}
impl From<DynamicsError> for PlanarTranslationError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<LuWorkspaceError> for PlanarTranslationError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::LinearSolve(value)
    }
}

impl From<SearchError<2>> for PlanarTranslationError {
    fn from(value: SearchError<2>) -> Self {
        match value {
            SearchError::Invalid(message) => Self::Invalid(message),
            SearchError::Stalled { residual } => Self::SearchStalled {
                residual_n: residual,
            },
            SearchError::LinearSolve(error) => Self::LinearSolve(error),
        }
    }
}

/// Accepted original flow, physical rigid endpoint and separate work accounts.
/// Body work alone is not a complete fluid/body energy or passivity theorem.
pub struct PlanarTranslationStep {
    pub flow: TranslatingAleStep,
    pub boundary_load: MaterialBoundaryLoad,
    pub momentum: MomentumImpulseReceipt,
    pub body_after: RigidBodyState,
    pub velocity_m_s: [f64; 2],
    pub translation_m: [f64; 2],
    pub flow_solves: usize,
    pub fluid_force_n: [f64; 2],
    /// Couple enforcing fixed orientation [N m]; its angular velocity is zero.
    pub restraint_torque_n_m: f64,
    pub interface_work_j: f64,
    /// Physical slip dissipation on the body [J], kept separately from body
    /// work, bulk viscosity and backward-Euler loss. No heat deposit is made.
    pub interface_slip_dissipation_j: f64,
    pub external_work_j: f64,
    pub body_kinetic_change_j: f64,
    pub spring_energy_change_j: f64,
    pub damping_loss_j: f64,
    /// Backward-Euler inertia/spring loss, never silently deposited as heat.
    pub mechanical_time_discretization_loss_j: f64,
    pub force_residual_n: [f64; 2],
    pub mechanical_work_residual_j: f64,
}

/// One original material/mesh owner, one translating native body, one clock.
pub struct PlanarTranslationMotion<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    flow: TranslatingAleFlow<'s, 'm, 'a, F>,
    properties: MassProperties,
    reference_pose: Pose,
    config: PlanarTranslationConfig,
    edges: Vec<usize>,
    tangential_resistance_pa_s_per_m: Option<f64>,
    body: RigidBodyState,
    schur: LuWorkspace,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> PlanarTranslationMotion<'s, 'm, 'a, F> {
    /// Body edges are impermeable no-slip walls. Initial BDM normal moments
    /// must agree with both native body velocity components. Rotation and
    /// out-of-plane momentum refuse rather than being projected away.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: PlanarTranslationConfig,
    ) -> Result<Self, PlanarTranslationError> {
        Self::with_wall_law(
            system,
            material,
            numerics,
            properties,
            body,
            initial_world_velocity,
            body_edges,
            config,
            None,
        )
    }

    /// The same freely translating body with impermeable Navier-slip edges.
    /// Normal momentum remains strongly coupled to the native fluid. Tangential
    /// fluid-on-body traction is beta*(u-body)_t [Pa], with the opposite fluid
    /// effort assembled by the original Robin owner. Finite beta >= 0 has units
    /// [Pa s/m]; zero gives perfect slip. Use new for exact no-slip.
    /// This coefficient is fixed on all selected body edges during this owner.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_slip(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: PlanarTranslationConfig,
        tangential_resistance_pa_s_per_m: f64,
    ) -> Result<Self, PlanarTranslationError> {
        Self::with_wall_law(
            system,
            material,
            numerics,
            properties,
            body,
            initial_world_velocity,
            body_edges,
            config,
            Some(tangential_resistance_pa_s_per_m),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn with_wall_law(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: PlanarTranslationConfig,
        tangential_resistance_pa_s_per_m: Option<f64>,
    ) -> Result<Self, PlanarTranslationError> {
        if tangential_resistance_pa_s_per_m.is_some_and(|beta| !nonnegative(beta)) {
            return Err(invalid(
                "finite nonnegative physical wall resistance required",
            ));
        }
        if !positive(config.depth_m)
            || !positive(config.velocity_tolerance_m_s)
            || !config
                .stiffness_n_m
                .into_iter()
                .chain(config.damping_n_s_m)
                .all(nonnegative)
            || initial_world_velocity.len() != system.n_u
            || !initial_world_velocity.iter().all(|v| v.is_finite())
            || body_edges.is_empty()
            || body.linear_momentum_world().z != 0.0
            || body.angular_momentum_body() != Vec3::ZERO
        {
            return Err(invalid(
                "finite planar coefficients, spin-free momentum and complete history required",
            ));
        }
        let mut edges = body_edges.to_vec();
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1])
            || edges
                .iter()
                .any(|&e| e >= system.mesh.edges.len() || !system.mesh.is_boundary(e))
        {
            return Err(invalid(
                "body interface requires unique original boundary edges",
            ));
        }
        let v = body.center_of_mass_velocity_world(properties)?;
        let velocity = [v.x, v.y];
        for &edge in &edges {
            let normal = system.mesh.edges[edge].normal;
            let expected = sum([velocity[0] * normal[0], velocity[1] * normal[1]]);
            if (initial_world_velocity[2 * edge] - expected).abs() > config.velocity_tolerance_m_s
                || initial_world_velocity[2 * edge + 1].abs() > config.velocity_tolerance_m_s
            {
                return Err(invalid("initial body and fluid normal velocities disagree"));
            }
        }
        Ok(Self {
            flow: TranslatingAleFlow::new(
                system,
                material,
                numerics,
                velocity,
                initial_world_velocity,
            )?,
            properties,
            reference_pose: body.pose(),
            config,
            edges,
            tangential_resistance_pa_s_per_m,
            body,
            schur: LuWorkspace::new(2)?,
        })
    }

    #[must_use]
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    #[must_use]
    pub const fn flow_state(&self) -> &TranslatingAleState {
        self.flow.state()
    }
    #[must_use]
    pub const fn material(&self) -> &F {
        self.flow.material()
    }

    /// Close both endpoint momentum equations from the SAME accepted histories.
    /// External force acts through the COM; pressure on unmeshed body surfaces
    /// must be included explicitly. All source/exterior callbacks must be pure
    /// across trials, and receive endpoint WORLD coordinates from native ALE.
    ///
    /// Two one-sided flow probes measure the full response. The original
    /// reusable LU owner solves a scaled 2-by-2 correction; backtracking admits
    /// only actual residual decrease within the velocity box. This is a local
    /// nonlinear method, not a uniqueness or global convergence claim. Every
    /// probe retains native nonlinear acceptance. Failure or final cancellation
    /// publishes neither body, fluid history, translation nor time.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialNavierStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        fluid_controls: NavierStokesControls,
        controls: PlanarTranslationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarTranslationStep, PlanarTranslationError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        if !positive(dt)
            || !positive(dt.recip())
            || !external_force_n.into_iter().all(f64::is_finite)
            || controls.maximum_flow_solves == 0
            || !controls.jacobian_step_m_s.into_iter().all(positive)
            || ![
                controls.force_absolute_n,
                controls.force_relative,
                controls.work_absolute_j,
                controls.work_relative,
            ]
            .into_iter()
            .all(positive)
            || controls
                .velocity_bounds_m_s
                .into_iter()
                .any(|[lo, hi]| !lo.is_finite() || !hi.is_finite() || !positive(hi - lo))
        {
            return Err(invalid(
                "finite timestep, velocity box, probe sizes, solve budget and positive tolerances required",
            ));
        }
        let cfg = self.config;
        let body0 = self.body;
        let original_velocity = body0.center_of_mass_velocity_world(self.properties)?;
        let v0 = [original_velocity.x, original_velocity.y];
        let q0 = self.flow.state().translation_m();
        let mass = self.properties.mass();
        let wall_velocity = Cell::new(v0);
        let wall = |_: [f64; 2]| wall_velocity.get();
        let mut boundaries = Vec::with_capacity(exterior.len() + self.edges.len());
        boundaries.extend_from_slice(exterior);
        boundaries.extend(self.edges.iter().map(
            |&edge| match self.tangential_resistance_pa_s_per_m {
                Some(beta) => MaterialNavierStokesBoundary::NavierSlip {
                    edge,
                    value: &wall,
                    tangential_resistance_pa_s_per_m: beta,
                },
                None => MaterialNavierStokesBoundary::Velocity { edge, value: &wall },
            },
        ));
        let mut search =
            EndpointSearch::new(v0, controls.velocity_bounds_m_s, controls.jacobian_step_m_s);
        let mut last_residual = [f64::NAN; 2];
        for index in 0..controls.maximum_flow_solves {
            poll(&mut cancelled)?;
            let requested = search.requested;
            let impulse = core::array::from_fn::<_, 2, _>(|i| mass * (requested[i] - v0[i]));
            let momentum = body0.apply_momentum_impulse(
                self.properties,
                Vec3::new(impulse[0], impulse[1], 0.0),
                Vec3::ZERO,
            )?;
            let v = momentum
                .state_after
                .center_of_mass_velocity_world(self.properties)?;
            let velocity = [v.x, v.y];
            if (0..2).any(|i| {
                velocity[i] < controls.velocity_bounds_m_s[i][0]
                    || velocity[i] > controls.velocity_bounds_m_s[i][1]
            }) {
                return Err(invalid(
                    "represented native velocity escapes the authored search box",
                ));
            }
            wall_velocity.set(velocity);
            let candidate = self.flow.candidate(
                dt,
                velocity,
                cfg.velocity_tolerance_m_s,
                force_world,
                &boundaries,
                fluid_controls,
                cfg.velocity_tolerance_m_s,
                &mut cancelled,
            )?;
            let flow = candidate.step();
            let q1 = flow.translation_m;
            let dx: [f64; 2] = core::array::from_fn(|i| q1[i] - q0[i]);
            let body = RigidBodyState::new(
                Pose::new(
                    self.reference_pose
                        .position_world()
                        .add(Vec3::new(q1[0], q1[1], 0.0)),
                    self.reference_pose.orientation(),
                )?,
                momentum.state_after.linear_momentum_world(),
                Vec3::ZERO,
            )?;
            let drift = body
                .pose()
                .position_world()
                .sub(body0.pose().position_world());
            if (0..2).any(|i| {
                let actual = [drift.x, drift.y][i] / dt;
                !actual.is_finite()
                    || (actual - velocity[i]).abs() > cfg.velocity_tolerance_m_s
                    || (dx[i] / dt - velocity[i]).abs() > cfg.velocity_tolerance_m_s
                    || (flow.mesh_velocity_m_s[i] - velocity[i]).abs() > cfg.velocity_tolerance_m_s
            }) {
                return Err(invalid("body and mesh drift are not jointly representable"));
            }
            let center = body.pose().position_world();
            let load =
                flow.boundary_loads_world(&self.edges, [center.x, center.y], pressure_offsets_pa)?;
            // Only original pressure/viscous stress acts on this impermeable
            // wall. Outward advective momentum is a separate native account.
            let force = load.stress.force_n_per_m.map(|f| cfg.depth_m * f);
            let increment = momentum.linear_increment_world;
            let inertia = [increment.x / dt, increment.y / dt];
            let spring: [f64; 2] = core::array::from_fn(|i| cfg.stiffness_n_m[i] * q1[i]);
            let damping: [f64; 2] = core::array::from_fn(|i| cfg.damping_n_s_m[i] * velocity[i]);
            let residual = core::array::from_fn(|i| {
                sum([
                    inertia[i],
                    -force[i],
                    -external_force_n[i],
                    spring[i],
                    damping[i],
                ])
            });
            let limits = core::array::from_fn(|i| {
                controls.force_absolute_n
                    + controls.force_relative
                        * maximum([
                            inertia[i],
                            force[i],
                            external_force_n[i],
                            spring[i],
                            damping[i],
                        ])
            });
            if !force
                .into_iter()
                .chain(inertia)
                .chain(spring)
                .chain(damping)
                .chain(residual)
                .chain(limits)
                .all(f64::is_finite)
            {
                return Err(invalid(
                    "unrepresentable physical momentum equation or allowance",
                ));
            }
            last_residual = residual;
            if (0..2).all(|i| residual[i].abs() <= limits[i]) {
                let interface_work = sum(core::array::from_fn::<_, 2, _>(|i| force[i] * dx[i]));
                let slip_loss = dt * cfg.depth_m * load.slip_dissipation_w_per_m;
                let external_work = sum(core::array::from_fn::<_, 2, _>(|i| {
                    external_force_n[i] * dx[i]
                }));
                let kinetic_change = momentum.work.kinetic_energy_change_j;
                let spring_change = sum(core::array::from_fn::<_, 2, _>(|i| {
                    0.5 * cfg.stiffness_n_m[i] * (q1[i] + q0[i]) * dx[i]
                }));
                let damping_loss = sum(core::array::from_fn::<_, 2, _>(|i| {
                    damping[i] * dt * velocity[i]
                }));
                let numerical_loss = sum([
                    impulse_time_loss(increment, mass),
                    0.5 * cfg.stiffness_n_m[0] * dx[0] * dx[0],
                    0.5 * cfg.stiffness_n_m[1] * dx[1] * dx[1],
                ]);
                let terms = [
                    kinetic_change,
                    spring_change,
                    damping_loss,
                    numerical_loss,
                    -interface_work,
                    -external_work,
                ];
                let work = sum(terms);
                let work_limit = controls.work_absolute_j + controls.work_relative * maximum(terms);
                let torque = -cfg.depth_m * load.stress.moment_n;
                if !terms
                    .into_iter()
                    .chain([work, work_limit, torque, slip_loss])
                    .all(f64::is_finite)
                {
                    return Err(invalid(
                        "unrepresentable mechanical work or restraint couple",
                    ));
                }
                if work.abs() > work_limit {
                    return Err(PlanarTranslationError::WorkBalance { residual_j: work });
                }
                poll(&mut cancelled)?;
                let flow = candidate.commit();
                self.body = body;
                return Ok(PlanarTranslationStep {
                    flow,
                    boundary_load: load,
                    momentum,
                    body_after: body,
                    velocity_m_s: velocity,
                    translation_m: q1,
                    flow_solves: index + 1,
                    fluid_force_n: force,
                    restraint_torque_n_m: torque,
                    interface_work_j: interface_work,
                    interface_slip_dissipation_j: slip_loss,
                    external_work_j: external_work,
                    body_kinetic_change_j: kinetic_change,
                    spring_energy_change_j: spring_change,
                    damping_loss_j: damping_loss,
                    mechanical_time_discretization_loss_j: numerical_loss,
                    force_residual_n: residual,
                    mechanical_work_residual_j: work,
                });
            }
            drop(candidate);
            if index + 1 < controls.maximum_flow_solves {
                search.observe(
                    Observation {
                        velocity,
                        residual,
                        limits,
                    },
                    &mut self.schur,
                )?;
            }
        }
        Err(PlanarTranslationError::SolveLimit {
            flow_solves: controls.maximum_flow_solves,
            residual_n: last_residual,
        })
    }
}

fn invalid(message: &'static str) -> PlanarTranslationError {
    PlanarTranslationError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlanarTranslationError> {
    if cancelled() {
        Err(PlanarTranslationError::Cancelled)
    } else {
        Ok(())
    }
}

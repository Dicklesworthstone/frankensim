//! Finite planar rigid momentum coupled to the original rotating ALE fluid.
//!
//! One bounded three-coordinate endpoint solve closes actual x/y momentum and
//! principal-z spin using the SAME native nonlinear flow history. Fluid force
//! and torque come from native pressure/viscous reactions about the moving COM.
//! The whole domain follows this body's translation and rotation. This is not
//! independent obstacle motion, deformation, or a complete coupled-energy theorem.

mod constraints;
mod rheology;

pub use self::constraints::{PlanarRigidConstraints, PlanarRigidSupport};

use self::rheology::{FlowTrial, Newtonian};
use super::search::{EndpointSearch, Observation, SearchError, maximum, sum};
use super::{nonnegative, positive};
use fs_flux::ale::rigid::{
    RigidAleError, RigidAleFlow, RigidAleMotion, RigidAleMotionTolerance, RigidAleState,
    RigidAleStep,
};
use fs_flux::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxNumerics, MaterialNavierStokesBoundary,
    NavierStokesControls,
};
use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use std::cell::Cell;

/// Whole-depth body coefficients. Springs rest at the construction pose; x/y
/// springs use fixed world axes and the torsional spring uses unwrapped angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarRigidConfig {
    pub depth_m: f64,
    pub translation_stiffness_n_m: [f64; 2],
    pub rotation_stiffness_n_m_rad: f64,
    pub translation_damping_n_s_m: [f64; 2],
    pub rotation_damping_n_m_s_rad: f64,
    /// Native body, mesh and strong normal-trace agreement [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Native spin, represented body/mesh angle and angular rate agreement [rad/s].
    pub angular_velocity_tolerance_rad_s: f64,
}

/// Explicit search region and separate dimensional acceptance for force/torque.
/// These are local nonlinear-search controls, not a root-existence guarantee.
/// A locked coordinate may use zero-width bounds containing zero and a zero
/// probe step. Free coordinates require a positive-width interval and probe.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarRigidControls {
    pub velocity_bounds_m_s: [[f64; 2]; 2],
    pub angular_velocity_bounds_rad_s: [f64; 2],
    pub jacobian_step_m_s: [f64; 2],
    pub jacobian_step_rad_s: f64,
    /// Includes all base, derivative and rejected backtracking flow solves.
    pub maximum_flow_solves: usize,
    pub force_absolute_n: f64,
    pub torque_absolute_n_m: f64,
    pub momentum_relative: f64,
    pub work_absolute_j: f64,
    pub work_relative: f64,
}

#[derive(Debug)]
pub enum PlanarRigidError {
    Invalid(&'static str),
    Flow(RigidAleError),
    Mechanics(DynamicsError),
    LinearSolve(LuWorkspaceError),
    InterfaceVelocity {
        edge: usize,
        residual_m_s: f64,
    },
    SolveLimit {
        flow_solves: usize,
        force_residual_n: [f64; 2],
        torque_residual_n_m: f64,
    },
    SearchStalled {
        force_residual_n: [f64; 2],
        torque_residual_n_m: f64,
    },
    WorkBalance {
        residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for PlanarRigidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar rigid ALE coupling: {self:?}")
    }
}
impl std::error::Error for PlanarRigidError {}
impl From<RigidAleError> for PlanarRigidError {
    fn from(value: RigidAleError) -> Self {
        Self::Flow(value)
    }
}
impl From<DynamicsError> for PlanarRigidError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<LuWorkspaceError> for PlanarRigidError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::LinearSolve(value)
    }
}
impl From<SearchError<3>> for PlanarRigidError {
    fn from(value: SearchError<3>) -> Self {
        match value {
            SearchError::Invalid(message) => Self::Invalid(message),
            SearchError::Stalled { residual } => Self::SearchStalled {
                force_residual_n: [residual[0], residual[1]],
                torque_residual_n_m: residual[2],
            },
            SearchError::LinearSolve(error) => Self::LinearSolve(error),
        }
    }
}

/// One jointly accepted body, original flow, moving geometry and physical clock.
pub struct PlanarRigidStep<G = ()> {
    /// Native flow and, for a rheological step, its original endpoint material owners.
    pub flow: RigidAleStep<G>,
    /// Original physical stress and transport, per unit depth, about the new COM.
    pub boundary_load: MaterialBoundaryLoad,
    /// Actual native fixed-pose momentum kick, before the subsequent rigid drift.
    pub momentum: MomentumImpulseReceipt,
    pub body_after: RigidBodyState,
    pub velocity_m_s: [f64; 2],
    pub angular_velocity_rad_s: f64,
    pub translation_m: [f64; 2],
    /// Signed unwrapped angle relative to the construction pose [rad].
    pub rotation_rad: f64,
    pub flow_solves: usize,
    pub fluid_force_n: [f64; 2],
    pub fluid_torque_n_m: f64,
    /// Physical reaction/impulse from ideal stationary supports, about the COM.
    pub support: PlanarRigidSupport,
    /// Endpoint stress paired with actual represented COM/angular drift [J].
    pub interface_work_j: f64,
    /// Physical endpoint wall-slip power times duration/depth, separate from work.
    pub interface_slip_dissipation_j: f64,
    /// Bulk fluid viscous heat [J], using BDF1 endpoint quadrature
    /// `dt * depth * integral(2 mu D:D)`. This uses the accepted rheology, excludes
    /// wall slip and numerical losses, and updates no thermal state. It is not
    /// an exact integral over the unobserved within-step trajectory.
    pub bulk_viscous_dissipation_j: f64,
    pub external_work_j: f64,
    pub body_kinetic_change_j: f64,
    pub spring_energy_change_j: f64,
    pub damping_loss_j: f64,
    /// Backward-Euler inertia/spring loss; no heat recipient is inferred.
    pub mechanical_time_discretization_loss_j: f64,
    pub force_residual_n: [f64; 2],
    pub torque_residual_n_m: f64,
    pub mechanical_work_residual_j: f64,
}

/// Finite x/y/angle motion, using a world-aligned principal body z axis.
/// Fluid geometry, boundary moments and torque arms follow the accepted pose.
pub struct PlanarRigidMotion<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    flow: RigidAleFlow<'s, 'm, 'a, F>,
    properties: MassProperties,
    masses: [f64; 3],
    reference_pose: Pose,
    config: PlanarRigidConfig,
    edges: Vec<usize>,
    tangential_resistance_pa_s_per_m: Option<f64>,
    constraints: PlanarRigidConstraints,
    body: RigidBodyState,
    rotation_rad: f64,
    schur: LuWorkspace,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> PlanarRigidMotion<'s, 'm, 'a, F> {
    /// Bind no-slip body edges on the initial world mesh. Initial strong normal
    /// means AND linear moments must match full rigid velocity, including spin.
    /// The body z axis must be principal and aligned with world z. No transverse
    /// angular momentum is projected away.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: PlanarRigidConfig,
    ) -> Result<Self, PlanarRigidError> {
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

    /// Bind the same impermeable body with finite tangential resistance [Pa s/m].
    /// Zero gives perfect slip; exact no-slip remains new(). The physical
    /// coefficient is fixed on the selected body edges, and loss stays explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_slip(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_world_velocity: &[f64],
        body_edges: &[usize],
        config: PlanarRigidConfig,
        tangential_resistance_pa_s_per_m: f64,
    ) -> Result<Self, PlanarRigidError> {
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
        config: PlanarRigidConfig,
        tangential_resistance_pa_s_per_m: Option<f64>,
    ) -> Result<Self, PlanarRigidError> {
        if !positive(config.depth_m)
            || !positive(config.velocity_tolerance_m_s)
            || !positive(config.angular_velocity_tolerance_rad_s)
            || !config
                .translation_stiffness_n_m
                .into_iter()
                .chain(config.translation_damping_n_s_m)
                .chain([
                    config.rotation_stiffness_n_m_rad,
                    config.rotation_damping_n_m_s_rad,
                ])
                .all(nonnegative)
            || tangential_resistance_pa_s_per_m.is_some_and(|beta| !nonnegative(beta))
            || initial_world_velocity.len() != system.n_u
            || !initial_world_velocity.iter().all(|v| v.is_finite())
            || body_edges.is_empty()
        {
            return Err(invalid(
                "finite planar coefficients, wall law and complete fluid history required",
            ));
        }
        let orientation = body.pose().orientation().components();
        let inertia_z = properties.angular_momentum_body_checked(Vec3::new(0.0, 0.0, 1.0))?;
        let angular = body.angular_momentum_body();
        if orientation[1] != 0.0
            || orientation[2] != 0.0
            || inertia_z.x != 0.0
            || inertia_z.y != 0.0
            || !positive(inertia_z.z)
            || body.linear_momentum_world().z != 0.0
            || angular.x != 0.0
            || angular.y != 0.0
        {
            return Err(invalid(
                "planar momenta and a world-aligned principal body z axis required",
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
        let rates = rates(body, properties)?;
        let center = body.pose().position_world();
        for &edge in &edges {
            let g = &system.mesh.edges[edge];
            let v = wall_velocity(rates, [center.x, center.y], g.mid);
            let a = system.mesh.verts[g.verts.0];
            let b = system.mesh.verts[g.verts.1];
            let linear =
                rates[2] * (-(b[1] - a[1]) * g.normal[0] + (b[0] - a[0]) * g.normal[1]) / 12.0;
            let errors = [
                initial_world_velocity[2 * edge] - sum([v[0] * g.normal[0], v[1] * g.normal[1]]),
                initial_world_velocity[2 * edge + 1] - linear,
            ];
            let residual = maximum(errors);
            if !errors.into_iter().all(f64::is_finite) || residual > config.velocity_tolerance_m_s {
                return Err(PlanarRigidError::InterfaceVelocity {
                    edge,
                    residual_m_s: residual,
                });
            }
        }
        Ok(Self {
            flow: RigidAleFlow::new(
                system,
                material,
                numerics,
                [center.x, center.y],
                motion(rates),
                initial_world_velocity,
            )?,
            properties,
            masses: [properties.mass(), properties.mass(), inertia_z.z],
            reference_pose: body.pose(),
            config,
            edges,
            tangential_resistance_pa_s_per_m,
            constraints: PlanarRigidConstraints::FREE,
            body,
            rotation_rad: 0.0,
            schur: LuWorkspace::new(3)?,
        })
    }

    #[must_use]
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    #[must_use]
    pub const fn flow_state(&self) -> &RigidAleState {
        self.flow.state()
    }
    #[must_use]
    pub const fn rotation_rad(&self) -> f64 {
        self.rotation_rad
    }
    #[must_use]
    pub const fn material(&self) -> &F {
        self.flow.material()
    }

    /// Solve actual body momentum from the same accepted histories.
    ///
    /// Stationary coordinates return their physical support reactions. The
    /// shared search probes only free coordinates (three for a free body)
    /// with the original reusable LU. Force and torque have separate
    /// dimensional allowances; descent uses fixed-base physical scaling.
    /// Every trial is an actual original nonlinear solve. Exterior callbacks
    /// receive endpoint world points/normals on the rigidly moving domain.
    /// External force acts through the COM; external torque acts about world z.
    ///
    /// Native kinematics preview supplies the exact admitted pivot/rates to
    /// wall callbacks, so an impermeable body does not acquire a fictitious
    /// reservoir flux. Body/mesh pose, actual momentum, work and cancellation
    /// must all pass before the unique flow candidate commits. Each rotation
    /// increment must stay below pi, as required by the native ALE owner.
    ///
    /// The body work gate is not a whole-fluid/body passivity or exact finite-step
    /// angular-conservation theorem. Original covariant BDF1 geometry/time defects,
    /// viscous dissipation, external streams and thermal ownership remain explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialNavierStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        fluid_controls: NavierStokesControls,
        controls: PlanarRigidControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarRigidStep, PlanarRigidError> {
        self.step_impl(
            timestep_s,
            force_world,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            external_torque_n_m,
            fluid_controls,
            controls,
            Newtonian,
            cancelled,
        )
    }

    // Constitutive dispatch changes only the original native fluid trial.
    // Geometry, actual body momentum, stress loads, work and atomic commit
    // retain one implementation for Newtonian and rheological materials.
    #[allow(clippy::too_many_arguments)]
    fn step_impl<G>(
        &mut self,
        timestep_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialNavierStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        fluid_controls: NavierStokesControls,
        controls: PlanarRigidControls,
        mut flow_trial: impl FlowTrial<F, G>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarRigidStep<G>, PlanarRigidError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        let free_axes = self.constraints.free_axes();
        let bounds = [
            controls.velocity_bounds_m_s[0],
            controls.velocity_bounds_m_s[1],
            controls.angular_velocity_bounds_rad_s,
        ];
        let steps = [
            controls.jacobian_step_m_s[0],
            controls.jacobian_step_m_s[1],
            controls.jacobian_step_rad_s,
        ];
        let external = [
            external_force_n[0],
            external_force_n[1],
            external_torque_n_m,
        ];
        if !positive(dt)
            || !positive(dt.recip())
            || !external.into_iter().all(f64::is_finite)
            || controls.maximum_flow_solves == 0
            || steps.into_iter().enumerate().any(|(i, h)| {
                if free_axes[i] {
                    !positive(h)
                } else {
                    !nonnegative(h)
                }
            })
            || ![
                controls.force_absolute_n,
                controls.torque_absolute_n_m,
                controls.work_absolute_j,
            ]
            .into_iter()
            .all(positive)
            || ![controls.momentum_relative, controls.work_relative]
                .into_iter()
                .all(nonnegative)
            || bounds.into_iter().enumerate().any(|(i, [lo, hi])| {
                !lo.is_finite()
                    || !hi.is_finite()
                    || if free_axes[i] {
                        !positive(hi - lo)
                    } else {
                        lo > 0.0 || hi < 0.0
                    }
            })
        {
            return Err(invalid(
                "finite time/tolerances and bounds/probes for free axes or stationary zero required",
            ));
        }
        let cfg = self.config;
        let body0 = self.body;
        let v0 = rates(body0, self.properties)?;
        let pivot0 = self.flow.state().pivot_world_m();
        let reference = self.reference_pose;
        let center0 = reference.position_world();
        let q0 = [
            pivot0[0] - center0.x,
            pivot0[1] - center0.y,
            self.rotation_rad,
        ];
        let stiffness = [
            cfg.translation_stiffness_n_m[0],
            cfg.translation_stiffness_n_m[1],
            cfg.rotation_stiffness_n_m_rad,
        ];
        let damping = [
            cfg.translation_damping_n_s_m[0],
            cfg.translation_damping_n_s_m[1],
            cfg.rotation_damping_n_m_s_rad,
        ];
        let motion_tolerance = RigidAleMotionTolerance {
            linear_velocity_m_s: cfg.velocity_tolerance_m_s,
            angular_velocity_rad_s: cfg.angular_velocity_tolerance_rad_s,
        };
        let wall_data = Cell::new((v0, pivot0));
        let wall = |point| {
            let (v, center) = wall_data.get();
            wall_velocity(v, center, point)
        };
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
        let mut search = EndpointSearch::with_free_axes(v0, bounds, steps, free_axes);
        let mut last_residual = [f64::NAN; 3];
        for index in 0..controls.maximum_flow_solves {
            poll(&mut cancelled)?;
            let impulse: [f64; 3] =
                core::array::from_fn(|i| self.masses[i] * (search.requested[i] - v0[i]));
            let momentum = body0.apply_momentum_impulse(
                self.properties,
                Vec3::new(impulse[0], impulse[1], 0.0),
                Vec3::new(0.0, 0.0, impulse[2]),
            )?;
            let velocity = rates(momentum.state_after, self.properties)?;
            if (0..3).any(|i| !free_axes[i] && velocity[i] != 0.0) {
                return Err(invalid(
                    "stationary support acquired nonzero native velocity",
                ));
            }
            if (0..3).any(|i| velocity[i] < bounds[i][0] || velocity[i] > bounds[i][1]) {
                return Err(invalid(
                    "represented native rates escape the authored search box",
                ));
            }
            let preview = self
                .flow
                .kinematics(dt, motion(velocity), motion_tolerance)?;
            let mesh = preview.mesh_velocity;
            wall_data.set((
                [
                    mesh.linear_velocity_world_m_s[0],
                    mesh.linear_velocity_world_m_s[1],
                    mesh.angular_velocity_rad_s,
                ],
                preview.pivot_world_m,
            ));
            let angle = preview.orientation_cos_sin[1].atan2(preview.orientation_cos_sin[0]);
            let pose = Pose::new(
                Vec3::new(
                    preview.pivot_world_m[0],
                    preview.pivot_world_m[1],
                    center0.z,
                ),
                reference
                    .orientation()
                    .right_exp(Vec3::new(0.0, 0.0, angle))?,
            )?;
            let body = RigidBodyState::new(
                pose,
                momentum.state_after.linear_momentum_world(),
                momentum.state_after.angular_momentum_body(),
            )?;
            let q1 = [
                preview.pivot_world_m[0] - center0.x,
                preview.pivot_world_m[1] - center0.y,
                self.rotation_rad + preview.rotation_increment_rad,
            ];
            let drift: [f64; 3] = core::array::from_fn(|i| q1[i] - q0[i]);
            let world_drift = pose.position_world().sub(body0.pose().position_world());
            let body_angle = relative_angle(body0.pose(), pose);
            let linear_errors = [
                drift[0] / dt - velocity[0],
                drift[1] / dt - velocity[1],
                world_drift.x / dt - velocity[0],
                world_drift.y / dt - velocity[1],
                mesh.linear_velocity_world_m_s[0] - velocity[0],
                mesh.linear_velocity_world_m_s[1] - velocity[1],
            ];
            let angular_errors = [
                drift[2] / dt - velocity[2],
                body_angle / dt - velocity[2],
                mesh.angular_velocity_rad_s - velocity[2],
            ];
            if !q1
                .into_iter()
                .chain(drift)
                .chain(linear_errors)
                .chain(angular_errors)
                .all(f64::is_finite)
                || maximum(linear_errors) > cfg.velocity_tolerance_m_s
                || maximum(angular_errors) > cfg.angular_velocity_tolerance_rad_s
            {
                return Err(invalid(
                    "body and mesh pose/rates are not jointly representable",
                ));
            }
            let candidate = flow_trial.candidate(
                &mut self.flow,
                dt,
                motion(velocity),
                motion_tolerance,
                force_world,
                &boundaries,
                fluid_controls,
                cfg.velocity_tolerance_m_s,
                &mut cancelled,
            )?;
            let flow = candidate.step();
            if flow.pivot_world_m != preview.pivot_world_m
                || flow.orientation_cos_sin != preview.orientation_cos_sin
                || flow.mesh_velocity != preview.mesh_velocity
            {
                return Err(invalid(
                    "native flow candidate does not retain its admitted kinematics",
                ));
            }
            let load =
                flow.boundary_loads_world(&self.edges, preview.pivot_world_m, pressure_offsets_pa)?;
            // Only native stress acts on the impermeable solid. Advected
            // momentum and slip dissipation retain separate original owners.
            let wrench = [
                load.stress.force_n_per_m[0] * cfg.depth_m,
                load.stress.force_n_per_m[1] * cfg.depth_m,
                load.stress.moment_n * cfg.depth_m,
            ];
            let dp = momentum.linear_increment_world;
            let dh = body.angular_momentum_body().z - body0.angular_momentum_body().z;
            let increments = [dp.x, dp.y, dh];
            let inertia = increments.map(|v| v / dt);
            let spring: [f64; 3] = core::array::from_fn(|i| stiffness[i] * q1[i]);
            let dashpot: [f64; 3] = core::array::from_fn(|i| damping[i] * velocity[i]);
            let physical_residual = core::array::from_fn(|i| {
                sum([inertia[i], -wrench[i], -external[i], spring[i], dashpot[i]])
            });
            let support = PlanarRigidSupport::resolve(free_axes, physical_residual, drift, dt)?;
            // A stationary support carries the entire missing load in its
            // locked coordinate. Free coordinates retain the original equation.
            let residual = core::array::from_fn(|i| {
                if free_axes[i] {
                    physical_residual[i]
                } else {
                    0.0
                }
            });
            let limits = core::array::from_fn(|i| {
                let absolute = if i == 2 {
                    controls.torque_absolute_n_m
                } else {
                    controls.force_absolute_n
                };
                absolute
                    + controls.momentum_relative
                        * maximum([inertia[i], wrench[i], external[i], spring[i], dashpot[i]])
            });
            if !wrench
                .into_iter()
                .chain(inertia)
                .chain(spring)
                .chain(dashpot)
                .chain(residual)
                .chain(limits)
                .all(f64::is_finite)
            {
                return Err(invalid(
                    "unrepresentable physical momentum equation or allowance",
                ));
            }
            last_residual = residual;
            if (0..3).all(|i| residual[i].abs() <= limits[i]) {
                let interface_work = sum(core::array::from_fn::<_, 3, _>(|i| wrench[i] * drift[i]));
                let external_work =
                    sum(core::array::from_fn::<_, 3, _>(|i| external[i] * drift[i]));
                let kinetic_change = momentum.work.kinetic_energy_change_j;
                let spring_change = sum(core::array::from_fn::<_, 3, _>(|i| {
                    0.5 * stiffness[i] * (q1[i] + q0[i]) * drift[i]
                }));
                let damping_loss = sum(core::array::from_fn::<_, 3, _>(|i| {
                    dashpot[i] * dt * velocity[i]
                }));
                let numerical_loss = sum(core::array::from_fn::<_, 3, _>(|i| {
                    let normalized = increments[i] / self.masses[i].sqrt();
                    0.5 * normalized * normalized + 0.5 * stiffness[i] * drift[i] * drift[i]
                }));
                let slip_loss = dt * cfg.depth_m * load.slip_dissipation_w_per_m;
                let viscous_heat = dt
                    * cfg.depth_m
                    * flow
                        .flow
                        .solution
                        .viscous_dissipation()
                        .total_power_w_per_m();
                let terms = [
                    kinetic_change,
                    spring_change,
                    damping_loss,
                    numerical_loss,
                    -interface_work,
                    -external_work,
                ];
                let work = sum(terms);
                let allowance = controls.work_absolute_j + controls.work_relative * maximum(terms);
                if !terms
                    .into_iter()
                    .chain([slip_loss, viscous_heat, work, allowance])
                    .all(f64::is_finite)
                {
                    return Err(invalid("unrepresentable physical or temporal work"));
                }
                if work.abs() > allowance {
                    return Err(PlanarRigidError::WorkBalance { residual_j: work });
                }
                poll(&mut cancelled)?;
                let flow = candidate.commit();
                self.body = body;
                self.rotation_rad = q1[2];
                return Ok(PlanarRigidStep {
                    flow,
                    boundary_load: load,
                    momentum,
                    body_after: body,
                    velocity_m_s: [velocity[0], velocity[1]],
                    angular_velocity_rad_s: velocity[2],
                    translation_m: [q1[0], q1[1]],
                    rotation_rad: q1[2],
                    flow_solves: index + 1,
                    fluid_force_n: [wrench[0], wrench[1]],
                    fluid_torque_n_m: wrench[2],
                    support,
                    interface_work_j: interface_work,
                    interface_slip_dissipation_j: slip_loss,
                    bulk_viscous_dissipation_j: viscous_heat,
                    external_work_j: external_work,
                    body_kinetic_change_j: kinetic_change,
                    spring_energy_change_j: spring_change,
                    damping_loss_j: damping_loss,
                    mechanical_time_discretization_loss_j: numerical_loss,
                    force_residual_n: [residual[0], residual[1]],
                    torque_residual_n_m: residual[2],
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
        Err(PlanarRigidError::SolveLimit {
            flow_solves: controls.maximum_flow_solves,
            force_residual_n: [last_residual[0], last_residual[1]],
            torque_residual_n_m: last_residual[2],
        })
    }
}

fn rates(body: RigidBodyState, properties: MassProperties) -> Result<[f64; 3], PlanarRigidError> {
    let v = body.center_of_mass_velocity_world(properties)?;
    let omega = properties.angular_velocity_body_checked(body.angular_momentum_body())?;
    if v.z != 0.0 || omega.x != 0.0 || omega.y != 0.0 {
        return Err(invalid("native body rates left the planar subspace"));
    }
    Ok([v.x, v.y, omega.z])
}
fn motion(rates: [f64; 3]) -> RigidAleMotion {
    RigidAleMotion {
        linear_velocity_world_m_s: [rates[0], rates[1]],
        angular_velocity_rad_s: rates[2],
    }
}
fn wall_velocity(rates: [f64; 3], center: [f64; 2], point: [f64; 2]) -> [f64; 2] {
    [
        (-rates[2]).mul_add(point[1] - center[1], rates[0]),
        rates[2].mul_add(point[0] - center[0], rates[1]),
    ]
}
fn relative_angle(reference: Pose, pose: Pose) -> f64 {
    let a = reference.orientation().components();
    let b = pose.orientation().components();
    let mut angle = 2.0 * (a[0] * b[3] - a[3] * b[0]).atan2(a[0] * b[0] + a[3] * b[3]);
    if angle > core::f64::consts::PI {
        angle -= core::f64::consts::TAU;
    }
    if angle < -core::f64::consts::PI {
        angle += core::f64::consts::TAU;
    }
    angle
}
fn invalid(message: &'static str) -> PlanarRigidError {
    PlanarRigidError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlanarRigidError> {
    if cancelled() {
        Err(PlanarRigidError::Cancelled)
    } else {
        Ok(())
    }
}

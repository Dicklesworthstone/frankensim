//! Implicit native Stokes force AND torque driving planar rigid momentum.
//!
//! The complete three-by-three fluid impedance couples x/y translation and
//! counterclockwise rotation about the reference COM. One loaded zero-motion
//! solve, three homogeneous unit-motion solves and a final endpoint solve use
//! the original material BDF1 flow. No prescribed body trajectory, fitted drag,
//! rotational added-inertia model or second fluid solver enters the update.
//!
//! Geometry and moment arms remain at the construction pose. Explicit bounds on
//! translation and angle limit this linearized, fixed-reference approximation;
//! the reported native body pose does not move the mesh. This is not ALE or a
//! finite-motion model. Body-frame z must be a principal axis aligned with world
//! z. Pressure datums, extrusion depth, source ownership and the native weak
//! tangential SIP treatment retain their original meanings.

mod feedback;

use super::heating::{
    PlanarStokesHeatError, PlanarStokesHeatStep, PlanarStokesObservedHeatError, StokesHeatBinding,
    StokesSlipHeatRoute, StokesThermalCells, unobserved_planar_error,
};
use super::{nonnegative, positive, solve_body_velocity};
use fs_flux::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxError, MaterialFluxNumerics, MaterialFluxSystem,
    MaterialStokesBoundary, MaterialStokesSolution,
};
use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, Pose, RigidBodyState, Vec3};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

/// Physical whole-body coefficients, including the declared extrusion depth.
/// Springs are independent in the fixed x/y/angle coordinates with rest at the
/// construction pose. Fluid forces retain all cross-coordinate coupling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarStokesConfig {
    /// Physical extrusion depth [m]; stress resultants otherwise remain per m.
    pub depth_m: f64,
    /// Maximum Euclidean COM translation from the construction pose [m].
    pub max_translation_m: f64,
    /// Maximum absolute angle [rad], positive and less than pi/2. The caller
    /// must choose a SMALL angle appropriate to the fixed reference geometry.
    pub max_rotation_rad: f64,
    /// Independent x/y spring stiffnesses [N/m].
    pub translation_stiffness_n_m: [f64; 2],
    /// Torsional stiffness [N m/rad].
    pub rotation_stiffness_n_m_rad: f64,
    /// Independent x/y mechanical dashpots [N s/m], excluding fluid viscosity.
    pub translation_damping_n_s_m: [f64; 2],
    /// Mechanical rotational dashpot [N m s/rad], excluding fluid viscosity.
    pub rotation_damping_n_m_s_rad: f64,
    /// Strong normal velocity and represented translation-rate tolerance [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Represented quaternion angular-drift tolerance [rad/s].
    pub angular_velocity_tolerance_rad_s: f64,
}

/// Separate force, torque and work tolerances; unlike units are never compared
/// against one absolute residual threshold.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarStokesControls {
    /// Absolute tolerance for each linear momentum equation [N].
    pub force_absolute_n: f64,
    /// Absolute tolerance for angular momentum [N m].
    pub torque_absolute_n_m: f64,
    /// Relative residual against each equation's largest physical term.
    pub momentum_relative: f64,
    /// Absolute tolerance for the body/spring/dashpot work identity [J].
    pub work_absolute_j: f64,
    /// Relative tolerance against the largest mechanical-work/storage term.
    pub work_relative: f64,
}

impl Default for PlanarStokesControls {
    fn default() -> Self {
        Self {
            force_absolute_n: 1e-9,
            torque_absolute_n_m: 1e-9,
            momentum_relative: 1e-8,
            work_absolute_j: 1e-10,
            work_relative: 1e-8,
        }
    }
}

/// Accepted physical history; unsuccessful steps never publish a partial state.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarStokesState {
    body: RigidBodyState,
    velocity: Vec<f64>,
    coordinates: [f64; 3],
    time_s: f64,
}

impl PlanarStokesState {
    /// Original momentum-form body with the accepted translated/rotated pose.
    #[must_use]
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Accepted native BDM1 velocity moments [m/s].
    #[must_use]
    pub fn fluid_velocity(&self) -> &[f64] {
        &self.velocity
    }
    /// Represented COM translation from the reference pose [m].
    #[must_use]
    pub const fn translation_m(&self) -> [f64; 2] {
        [self.coordinates[0], self.coordinates[1]]
    }
    /// Represented counterclockwise rotation from the reference pose [rad].
    #[must_use]
    pub const fn rotation_rad(&self) -> f64 {
        self.coordinates[2]
    }
    /// Accepted physical time [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
}

/// A refusal preserves body, fluid history, pose and clock together.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarStokesError {
    /// Invalid units, planar-axis admission, boundary ownership or arithmetic.
    Invalid(&'static str),
    /// Original material/flow refusal.
    Flow(MaterialFluxError),
    /// Original rigid momentum or pose refusal.
    Mechanics(DynamicsError),
    /// Original dense-solve refusal, including a singular coupled operator.
    LinearSolve(LuWorkspaceError),
    /// A strong normal mean or linear moment disagrees with rigid motion.
    InterfaceVelocity { edge: usize, residual_m_s: f64 },
    /// Translation escapes the declared fixed-reference approximation domain.
    TranslationLimit { proposed_m: f64, limit_m: f64 },
    /// Rotation escapes the declared fixed-reference approximation domain.
    RotationLimit { proposed_rad: f64, limit_rad: f64 },
    /// Actual endpoint force, torque or body-work identity fails its gate.
    Balance {
        force_residual_n: [f64; 2],
        torque_residual_n_m: f64,
        mechanical_work_residual_j: f64,
    },
    /// Cancellation between native solves or before physical publication.
    Cancelled,
}

impl core::fmt::Display for PlanarStokesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar Stokes motion refused: {self:?}")
    }
}
impl std::error::Error for PlanarStokesError {}
impl From<MaterialFluxError> for PlanarStokesError {
    fn from(value: MaterialFluxError) -> Self {
        Self::Flow(value)
    }
}
impl From<DynamicsError> for PlanarStokesError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}
impl From<LuWorkspaceError> for PlanarStokesError {
    fn from(value: LuWorkspaceError) -> Self {
        Self::LinearSolve(value)
    }
}

/// Accepted original flow, rigid momentum and explicit physical work accounts.
pub struct PlanarStokesStep {
    /// Native flow at the actual accepted rigid endpoint rates.
    pub flow: MaterialStokesSolution,
    /// Native pressure/viscous resultant per metre depth, about the reference COM.
    pub boundary_load: MaterialBoundaryLoad,
    /// Native fixed-pose impulse receipt, before the subsequent pose drift.
    pub momentum: MomentumImpulseReceipt,
    /// Accepted native body, including translation and quaternion drift.
    pub body_after: RigidBodyState,
    /// New time [s].
    pub time_s: f64,
    /// Endpoint COM velocity in fixed world x/y coordinates [m/s].
    pub velocity_m_s: [f64; 2],
    /// Endpoint counterclockwise angular velocity [rad/s].
    pub angular_velocity_rad_s: f64,
    /// Actual whole-depth fluid-on-body force [N].
    pub fluid_force_n: [f64; 2],
    /// Actual whole-depth fluid-on-body torque about the reference COM [N m].
    pub fluid_torque_n_m: f64,
    /// Full derivative of [Fx, Fy, Tz] with respect to [vx, vy, omega]. Rows
    /// have units [N, N, N m]; columns divide by [m/s, m/s, rad/s]. Computed
    /// homogeneously, without subtracting large absolute-pressure loads.
    pub fluid_impedance: [[f64; 3]; 3],
    /// `dt * (F dot v + T omega)` [J], fluid on body. Opposite interface signs
    /// do not establish a complete fluid/body energy or passivity audit.
    pub interface_work_j: f64,
    /// Endpoint authored external force/torque work [J].
    pub external_work_j: f64,
    /// Native body kinetic-energy change [J].
    pub body_kinetic_change_j: f64,
    /// Spring storage change using represented pose coordinates [J].
    pub spring_energy_change_j: f64,
    /// Physical mechanical dashpot loss [J], excluding resolved fluid viscosity.
    pub damping_loss_j: f64,
    /// Backward-Euler mechanical inertia/spring loss [J]; numerical, not heat.
    pub mechanical_time_discretization_loss_j: f64,
    /// Native fluid kinetic-energy change over the declared depth [J].
    pub fluid_kinetic_change_j: f64,
    /// Bulk fluid viscous heat [J], using BDF1 endpoint quadrature
    /// `dt * depth * integral(2 mu D:D)`. This excludes wall slip and numerical
    /// losses and is not an exact trajectory integral. Standalone `step`
    /// updates no thermal state; `step_with_heating` and `step_with_slip_heating`
    /// deposit this heat through the admitted finite-cell receiver map.
    pub bulk_viscous_dissipation_j: f64,
    /// Actual represented linear momentum residuals [N].
    pub force_residual_n: [f64; 2],
    /// Actual represented angular momentum residual [N m].
    pub torque_residual_n_m: f64,
    /// Body storage plus physical/numerical loss minus applied/interface work [J].
    pub mechanical_work_residual_j: f64,
}

/// One fixed native mesh/material, one planar rigid body and one accepted clock.
pub struct PlanarStokesMotion<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'m>,
    fluid: MaterialFluxSystem<'s, 'm, 'a, F>,
    properties: MassProperties,
    masses: [f64; 3],
    reference_pose: Pose,
    config: PlanarStokesConfig,
    edges: Vec<usize>,
    schur: LuWorkspace,
    state: PlanarStokesState,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> PlanarStokesMotion<'s, 'm, 'a, F> {
    /// Bind the original owners without reinterpreting native body coordinates.
    /// `MassProperties` admits only a COM reference; the pose position is that
    /// COM. The body-frame z axis must exactly align with world z and be a
    /// principal inertia axis. Initial out-of-plane momenta are refused.
    /// Initial strong normal moments, including the nonzero linear moment of
    /// rotation, must agree with the body. Tangential velocity remains weak SIP.
    ///
    /// # Errors
    /// Invalid planar data, history, selected edges or original material/rigid
    /// admission. A boundary selection must permit every rigid unit motion
    /// together with the exterior data; native incompressibility checks this
    /// during stepping (e.g. a sealed piston cannot translate normally).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        properties: MassProperties,
        body: RigidBodyState,
        initial_fluid_velocity: Vec<f64>,
        body_edges: &[usize],
        config: PlanarStokesConfig,
    ) -> Result<Self, PlanarStokesError> {
        if ![
            config.depth_m,
            config.max_translation_m,
            config.max_rotation_rad,
            config.velocity_tolerance_m_s,
            config.angular_velocity_tolerance_rad_s,
        ]
        .into_iter()
        .all(positive)
            || config.max_rotation_rad >= core::f64::consts::FRAC_PI_2
            || ![
                config.translation_stiffness_n_m[0],
                config.translation_stiffness_n_m[1],
                config.rotation_stiffness_n_m_rad,
                config.translation_damping_n_s_m[0],
                config.translation_damping_n_s_m[1],
                config.rotation_damping_n_m_s_rad,
            ]
            .into_iter()
            .all(nonnegative)
            || initial_fluid_velocity.len() != system.n_u
            || !initial_fluid_velocity.iter().all(|v| v.is_finite())
            || body_edges.is_empty()
        {
            return Err(PlanarStokesError::Invalid(
                "finite physical data, small-motion bounds and matching history required",
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
            return Err(PlanarStokesError::Invalid(
                "planar momenta and a world-aligned principal body z axis required",
            ));
        }
        let mut edges = body_edges.to_vec();
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1])
            || edges
                .iter()
                .any(|&edge| edge >= system.mesh.edges.len() || !system.mesh.is_boundary(edge))
        {
            return Err(PlanarStokesError::Invalid(
                "body interface must contain unique boundary edges",
            ));
        }
        let result = Self {
            system,
            fluid: MaterialFluxSystem::new(system, material, numerics)?,
            properties,
            masses: [properties.mass(), properties.mass(), inertia_z.z],
            reference_pose: body.pose(),
            config,
            edges,
            schur: LuWorkspace::new(3)?,
            state: PlanarStokesState {
                body,
                velocity: initial_fluid_velocity,
                coordinates: [0.0; 3],
                time_s: 0.0,
            },
        };
        result.check_interface(&result.state.velocity, result.rates(body)?)?;
        Ok(result)
    }

    /// Entire accepted physical history, unchanged by every failed step.
    #[must_use]
    pub const fn state(&self) -> &PlanarStokesState {
        &self.state
    }
    /// Original construction source. Ordinary methods retain its frozen transport;
    /// temperature-feedback steps retain their resolved source in the returned receipt.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.fluid.material()
    }

    fn rates(&self, body: RigidBodyState) -> Result<[f64; 3], PlanarStokesError> {
        let v = body.center_of_mass_velocity_world(self.properties)?;
        let omega = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        if v.z != 0.0 || omega.x != 0.0 || omega.y != 0.0 {
            return Err(PlanarStokesError::Invalid(
                "native body rates left the planar subspace",
            ));
        }
        Ok([v.x, v.y, omega.z])
    }

    fn boundary_velocity(&self, rates: [f64; 3], point: [f64; 2]) -> [f64; 2] {
        let origin = self.reference_pose.position_world();
        [
            rates[0] - rates[2] * (point[1] - origin.y),
            rates[1] + rates[2] * (point[0] - origin.x),
        ]
    }

    fn check_interface(&self, velocity: &[f64], rates: [f64; 3]) -> Result<(), PlanarStokesError> {
        for &edge in &self.edges {
            let geometry = &self.system.mesh.edges[edge];
            let n = geometry.normal;
            let u = self.boundary_velocity(rates, geometry.mid);
            let a = self.system.mesh.verts[geometry.verts.0];
            let b = self.system.mesh.verts[geometry.verts.1];
            // BDM1 stores means against 1 and normalized s in [-1/2,1/2].
            // Rigid rotation's linear normal moment does NOT vanish.
            let linear = rates[2] * (-(b[1] - a[1]) * n[0] + (b[0] - a[0]) * n[1]) / 12.0;
            let mean_error = (velocity[2 * edge] - (u[0] * n[0] + u[1] * n[1])).abs();
            let linear_error = (velocity[2 * edge + 1] - linear).abs();
            let residual = mean_error.max(linear_error);
            if !mean_error.is_finite()
                || !linear_error.is_finite()
                || residual > self.config.velocity_tolerance_m_s
            {
                return Err(PlanarStokesError::InterfaceVelocity {
                    edge,
                    residual_m_s: residual,
                });
            }
        }
        Ok(())
    }

    fn solve(
        &self,
        rates: [f64; 3],
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        previous: &[f64],
        dt: f64,
    ) -> Result<MaterialStokesSolution, PlanarStokesError> {
        let motion = |point| self.boundary_velocity(rates, point);
        Ok(solve_body_velocity(
            &self.fluid,
            &self.edges,
            &motion,
            force,
            exterior,
            previous,
            dt,
        )?)
    }

    fn load(
        &self,
        flow: &MaterialStokesSolution,
        offsets: &[Option<f64>],
    ) -> Result<(MaterialBoundaryLoad, [f64; 3]), PlanarStokesError> {
        let origin = self.reference_pose.position_world();
        let load = flow.boundary_loads(&self.edges, [origin.x, origin.y], offsets)?;
        let wrench = [
            load.stress.force_n_per_m[0],
            load.stress.force_n_per_m[1],
            load.stress.moment_n,
        ]
        .map(|v| self.config.depth_m * v);
        if !wrench.iter().all(|v| v.is_finite()) {
            return Err(PlanarStokesError::Invalid("physical fluid wrench overflow"));
        }
        Ok((load, wrench))
    }

    fn check_limits(&self, q: [f64; 3]) -> Result<(), PlanarStokesError> {
        let translation = q[0].hypot(q[1]);
        if !translation.is_finite() || translation > self.config.max_translation_m {
            return Err(PlanarStokesError::TranslationLimit {
                proposed_m: translation,
                limit_m: self.config.max_translation_m,
            });
        }
        if !q[2].is_finite() || q[2].abs() > self.config.max_rotation_rad {
            return Err(PlanarStokesError::RotationLimit {
                proposed_rad: q[2],
                limit_rad: self.config.max_rotation_rad,
            });
        }
        Ok(())
    }

    /// Close all three body momentum equations at the same native fluid endpoint.
    /// Authored external force [N] and torque [N m] are about the reference COM
    /// and must include pressure on any unmodeled side of the body. Remaining
    /// boundary edges belong to `exterior`. Pressure offsets [Pa] obey the native
    /// per-component mean-zero/traction-datum contract.
    ///
    /// Five native solves per step retain full translation/rotation impedance.
    /// Callbacks must be pure and fixed throughout. Cancellation is checked
    /// between solves and before publication, not inside the synchronous solver.
    /// The body work identity is independent of a whole-fluid energy audit;
    /// resolved viscous and weak SIP terms retain their native scope.
    ///
    /// # Errors
    /// Invalid time/controls, native solver errors, geometry-domain escape,
    /// unrepresentable pose/rates, failed balances or cancellation leave the
    /// entire accepted physical history unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesStep, PlanarStokesError> {
        let (state, step) = self.propose_step(
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            external_torque_n_m,
            controls,
            cancelled,
        )?;
        self.state = state;
        Ok(step)
    }

    // The native momentum solve is shared with thermal composition. Only the
    // reusable linear-solve workspace changes before physical acceptance.
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    fn propose_step(
        &mut self,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(PlanarStokesState, PlanarStokesStep), PlanarStokesError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        let time = self.state.time_s + dt;
        let external = [
            external_force_n[0],
            external_force_n[1],
            external_torque_n_m,
        ];
        if !positive(dt)
            || !positive(dt.recip())
            || !time.is_finite()
            || time <= self.state.time_s
            || !external.iter().all(|v| v.is_finite())
            || ![
                controls.force_absolute_n,
                controls.torque_absolute_n_m,
                controls.momentum_relative,
                controls.work_absolute_j,
                controls.work_relative,
            ]
            .into_iter()
            .all(positive)
        {
            return Err(PlanarStokesError::Invalid(
                "representable positive timestep/tolerances and finite external wrench required",
            ));
        }
        let v0 = self.rates(self.state.body)?;
        let q0 = self.state.coordinates;
        let base = self.solve([0.0; 3], force_density, exterior, &self.state.velocity, dt)?;
        let (_, base_wrench) = self.load(&base, pressure_offsets_pa)?;
        drop(base);
        poll(&mut cancelled)?;
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
        let response_offsets: Vec<_> = pressure_offsets_pa.iter().map(|p| p.map(|_| 0.0)).collect();
        let mut impedance = [[0.0; 3]; 3];
        for col in 0..3 {
            let mut unit = [0.0; 3];
            unit[col] = 1.0;
            let response = self.solve(unit, &zero, &homogeneous, &zero_history, dt)?;
            let (_, wrench) = self.load(&response, &response_offsets)?;
            for row in 0..3 {
                impedance[row][col] = wrench[row];
            }
            drop(response);
            poll(&mut cancelled)?;
        }
        let stiffness = [
            self.config.translation_stiffness_n_m[0],
            self.config.translation_stiffness_n_m[1],
            self.config.rotation_stiffness_n_m_rad,
        ];
        let damping = [
            self.config.translation_damping_n_s_m[0],
            self.config.translation_damping_n_s_m[1],
            self.config.rotation_damping_n_m_s_rad,
        ];
        // M/dt + C + dt*K - Z closes the complete fluid/body Schur complement.
        // Congruence with sqrt(M) gives all three unknowns common kinetic units
        // and avoids comparing unscaled metres with radians in pivot selection.
        let roots = self.masses.map(f64::sqrt);
        let mut matrix = [0.0; 9];
        let mut rhs = [0.0; 3];
        for row in 0..3 {
            for col in 0..3 {
                matrix[3 * row + col] = -impedance[row][col] / roots[row] / roots[col];
            }
            matrix[4 * row] += dt.recip()
                + damping[row] / self.masses[row]
                + dt * stiffness[row] / self.masses[row];
            rhs[row] = roots[row] * v0[row] / dt
                + (external[row] + base_wrench[row] - stiffness[row] * q0[row]) / roots[row];
        }
        let mut scaled_velocity = [0.0; 3];
        self.schur.solve_into(&matrix, &rhs, &mut scaled_velocity)?;
        let proposed: [f64; 3] = core::array::from_fn(|i| scaled_velocity[i] / roots[i]);
        let impulse: [f64; 3] = core::array::from_fn(|i| self.masses[i] * (proposed[i] - v0[i]));
        let momentum = self.state.body.apply_momentum_impulse(
            self.properties,
            Vec3::new(impulse[0], impulse[1], 0.0),
            Vec3::new(0.0, 0.0, impulse[2]),
        )?;
        let velocity = self.rates(momentum.state_after)?;
        let drift = velocity.map(|v| dt * v);
        let proposed_q = core::array::from_fn(|i| q0[i] + drift[i]);
        self.check_limits(proposed_q)?;
        let reference = self.reference_pose;
        let pose = Pose::new(
            reference
                .position_world()
                .add(Vec3::new(proposed_q[0], proposed_q[1], 0.0)),
            reference
                .orientation()
                .right_exp(Vec3::new(0.0, 0.0, proposed_q[2]))?,
        )?;
        let position = pose.position_world().sub(reference.position_world());
        let q1 = [position.x, position.y, relative_angle(reference, pose)];
        self.check_limits(q1)?;
        let actual_drift: [f64; 3] = core::array::from_fn(|i| q1[i] - q0[i]);
        let world_drift = pose
            .position_world()
            .sub(self.state.body.pose().position_world());
        let position_tolerance = self.config.velocity_tolerance_m_s * dt;
        if (actual_drift[0] - drift[0]).hypot(actual_drift[1] - drift[1]) > position_tolerance
            || (world_drift.x - drift[0]).hypot(world_drift.y - drift[1]) > position_tolerance
            || (actual_drift[2] - drift[2]).abs()
                > self.config.angular_velocity_tolerance_rad_s * dt
        {
            return Err(PlanarStokesError::Invalid(
                "body pose drift is not representable at the declared rate tolerances",
            ));
        }
        let body = RigidBodyState::new(
            pose,
            momentum.state_after.linear_momentum_world(),
            momentum.state_after.angular_momentum_body(),
        )?;
        let flow = self.solve(velocity, force_density, exterior, &self.state.velocity, dt)?;
        let (load, wrench) = self.load(&flow, pressure_offsets_pa)?;
        self.check_interface(flow.velocity_dofs(), velocity)?;
        poll(&mut cancelled)?;
        let dp = body
            .linear_momentum_world()
            .sub(self.state.body.linear_momentum_world());
        let dl = body.angular_momentum_body().z - self.state.body.angular_momentum_body().z;
        let inertia = [dp.x / dt, dp.y / dt, dl / dt];
        let mut residual = [0.0; 3];
        let mut momentum_ok = true;
        let mut spring_change = 0.0;
        let mut damping_loss = 0.0;
        let mut temporal_loss = 0.0;
        let mut interface_work = 0.0;
        let mut external_work = 0.0;
        for i in 0..3 {
            let spring = stiffness[i] * q1[i];
            let dashpot = damping[i] * velocity[i];
            residual[i] = inertia[i] - (wrench[i] + external[i] - spring - dashpot);
            let scale = [inertia[i], wrench[i], external[i], spring, dashpot]
                .into_iter()
                .map(f64::abs)
                .fold(0.0, f64::max);
            let absolute = if i == 2 {
                controls.torque_absolute_n_m
            } else {
                controls.force_absolute_n
            };
            momentum_ok &= residual[i].is_finite()
                && residual[i].abs() <= absolute + controls.momentum_relative * scale;
            let dv = velocity[i] - v0[i];
            spring_change += 0.5 * stiffness[i] * (q1[i] + q0[i]) * actual_drift[i];
            damping_loss += dashpot * drift[i];
            temporal_loss += 0.5 * self.masses[i] * dv * dv
                + 0.5 * stiffness[i] * actual_drift[i] * actual_drift[i];
            interface_work += wrench[i] * drift[i];
            external_work += external[i] * drift[i];
        }
        let kinetic_change = momentum.work.kinetic_energy_change_j;
        let work_residual = kinetic_change + spring_change + damping_loss + temporal_loss
            - interface_work
            - external_work;
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
        if !momentum_ok
            || !work_residual.is_finite()
            || work_residual.abs() > controls.work_absolute_j + controls.work_relative * work_scale
        {
            return Err(PlanarStokesError::Balance {
                force_residual_n: [residual[0], residual[1]],
                torque_residual_n_m: residual[2],
                mechanical_work_residual_j: work_residual,
            });
        }
        let fluid_energy_change = self.config.depth_m
            * (flow.kinetic_energy_j_per_m()
                - self.fluid.transport().density_kg_m3()
                    * self.system.kinetic_energy(&self.state.velocity));
        let viscous_heat =
            dt * self.config.depth_m * flow.viscous_dissipation().total_power_w_per_m();
        if !fluid_energy_change.is_finite() || !viscous_heat.is_finite() {
            return Err(PlanarStokesError::Invalid("physical fluid energy overflow"));
        }
        let accepted_velocity = flow.velocity_dofs().to_vec();
        poll(&mut cancelled)?;
        let state = PlanarStokesState {
            body,
            velocity: accepted_velocity,
            coordinates: q1,
            time_s: time,
        };
        let step = PlanarStokesStep {
            flow,
            boundary_load: load,
            momentum,
            body_after: body,
            time_s: time,
            velocity_m_s: [velocity[0], velocity[1]],
            angular_velocity_rad_s: velocity[2],
            fluid_force_n: [wrench[0], wrench[1]],
            fluid_torque_n_m: wrench[2],
            fluid_impedance: impedance,
            interface_work_j: interface_work,
            external_work_j: external_work,
            body_kinetic_change_j: kinetic_change,
            spring_energy_change_j: spring_change,
            damping_loss_j: damping_loss,
            mechanical_time_discretization_loss_j: temporal_loss,
            fluid_kinetic_change_j: fluid_energy_change,
            bulk_viscous_dissipation_j: viscous_heat,
            force_residual_n: [residual[0], residual[1]],
            torque_residual_n_m: residual[2],
            mechanical_work_residual_j: work_residual,
        };
        Ok((state, step))
    }

    /// Advance planar flow/body motion and deposit its bulk strain heat.
    ///
    /// The same fixed-reference triangular-prism map as the guided path binds
    /// the original finite thermal cells: geometry, density-derived mass and
    /// clocks must match. Each cell retains its authored caloric law, sources,
    /// reservoirs and contacts. The caller must bind the same physical fluid's
    /// caloric law; a frozen transport snapshot alone cannot establish that.
    ///
    /// Heat is `dt * depth * integral(2 mu D:D)` from the actual accepted
    /// endpoint, held during native thermal substeps. This is a first-order,
    /// one-way split without thermal advection or viscosity feedback. Resolve
    /// the approximation by refining both timesteps. Body dashpots, mechanical
    /// work and numerical losses are not thermal sources.
    ///
    /// # Errors
    /// Admission, motion, caloric-domain, budget and cancellation failures leave
    /// fluid velocity, body momentum/pose, thermal state and both clocks intact.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heating(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesHeatStep, PlanarStokesHeatError> {
        self.step_heating(
            thermal, &mut (), cells, None, timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, external_torque_n_m, controls,
            thermal_config, cancelled,
        )
        .map_err(unobserved_planar_error)
    }

    /// Also deposit native Navier-slip loss with an explicit fluid/wall split.
    ///
    /// Every exterior Navier-slip edge, including a perfect-slip edge, needs
    /// exactly one route to an active wall cell outside the fluid map. Moving
    /// body edges remain no-slip and cannot be routed. Bulk and slip heating use
    /// the same accepted flow; material and boundary callbacks are not queried
    /// again. All admission and approximation limits of
    /// [`Self::step_with_heating`] apply.
    ///
    /// # Errors
    /// Invalid/missing routes, original solver or thermal failures and
    /// cancellation preserve every accepted owner, including late refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_slip_heating(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesHeatStep, PlanarStokesHeatError> {
        self.step_heating(
            thermal, &mut (), cells, Some(routes), timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, external_torque_n_m, controls,
            thermal_config, cancelled,
        )
        .map_err(unobserved_planar_error)
    }

    /// Advance planar viscous heating with another owner on the common thermal clock.
    ///
    /// The same observer composition as the guided path forwards every source,
    /// transition, interval, conductance, exchange and timestep hook. Additional
    /// heat enters the native audit once; the separate Stokes receipts retain
    /// only the original bulk and interfacial dissipation.
    ///
    /// The observer's `Clone` must isolate mutable physical history according to
    /// `DepositionObserver`. Motion, thermal state and accepted observer history
    /// publish together after all acceptance checks and the final cancellation
    /// poll. Additional thermal physics does not change the fixed reference mesh
    /// or the frozen fluid transport used by the native planar motion solve.
    /// The time-refinement requirements of [`Self::step_with_heating`] apply.
    ///
    /// # Errors
    /// Native map, motion, thermal, observer or cancellation failures preserve
    /// all three inputs, including refusals after successful thermal substeps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heating_observed<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesHeatStep, PlanarStokesObservedHeatError<O::Error>> {
        self.step_heating(
            thermal, observer, cells, None, timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, external_torque_n_m, controls,
            thermal_config, cancelled,
        )
    }

    /// Also deposit routed Navier-slip heat on that same observed timeline.
    ///
    /// The route and partition requirements of [`Self::step_with_slip_heating`]
    /// and transaction requirements of [`Self::step_with_heating_observed`] apply.
    /// Bulk and slip heating share the same accepted native flow and heat owner.
    ///
    /// # Errors
    /// Route admission and native motion, thermal, observer or cancellation
    /// failures preserve every accepted owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_slip_heating_observed<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesHeatStep, PlanarStokesObservedHeatError<O::Error>> {
        self.step_heating(
            thermal, observer, cells, Some(routes), timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, external_torque_n_m, controls,
            thermal_config, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_heating<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        routes: Option<&[StokesSlipHeatRoute]>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarStokesHeatStep, PlanarStokesObservedHeatError<O::Error>> {
        poll(&mut cancelled)?;
        let binding = StokesHeatBinding {
            mesh: self.system.mesh,
            density_kg_m3: self.fluid.transport().density_kg_m3(),
            depth_m: self.config.depth_m,
            time_s: self.state.time_s,
        };
        binding.check_thermal_cells(thermal, cells, &mut cancelled)?;
        let routes = match routes {
            Some(routes) => binding.check_slip_heat_routes(
                thermal, cells, &self.edges, exterior, routes, &mut cancelled,
            )?,
            None => std::collections::BTreeMap::new(),
        };
        let (state, motion) = self.propose_step(
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            external_torque_n_m,
            controls,
            &mut cancelled,
        )?;
        let (candidate, accepted_observer, audit, heat_j) = binding.advance_heat_observed(
            thermal,
            observer,
            cells,
            &routes,
            &motion.flow,
            pressure_offsets_pa,
            timestep_s,
            state.time_s,
            motion.bulk_viscous_dissipation_j,
            thermal_config,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        self.state = state;
        *thermal = candidate;
        *observer = accepted_observer;
        Ok(PlanarStokesHeatStep {
            motion,
            thermal: audit,
            deposited_viscous_heat_j: heat_j[0],
            deposited_slip_heat_j: heat_j[1],
        })
    }
}

fn relative_angle(reference: Pose, pose: Pose) -> f64 {
    let a = reference.orientation().components();
    let b = pose.orientation().components();
    let mut angle = 2.0 * (a[0] * b[3] - a[3] * b[0]).atan2(a[0] * b[0] + a[3] * b[3]);
    // Canonical quaternion sign can flip when the reference lies near pi.
    if angle > core::f64::consts::PI {
        angle -= core::f64::consts::TAU;
    }
    if angle < -core::f64::consts::PI {
        angle += core::f64::consts::TAU;
    }
    angle
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlanarStokesError> {
    if cancelled() {
        Err(PlanarStokesError::Cancelled)
    } else {
        Ok(())
    }
}

//! Native incompressible flow on globally rigid translating/rotating domains.
//!
//! x = X + R(y-c), with fixed reference pivot c and proper planar rotation R.
//! The unknown a(y) = R^T u_world(x) is ABSOLUTE fluid velocity in mesh axes.
//! Relative advection is a - R^T V - omega J(y-c), while the carried velocity,
//! viscous stress, normal boundary data and physical source all use a.
//!
//! BDF1 integrates (a_new - R_new^T R_old a_old)/dt against the original
//! cell bases. The old field is rotated inside the mass quadrature: rotating
//! global H(div) edge coefficients would not preserve normal continuity.
//! Constant world velocity is therefore preserved without a Coriolis proxy.
//! The isometry retains the original metrics, area and zero-divergence mesh
//! velocity. No moved mesh, remap or second spatial operator is introduced.
//!
//! Geometry uses endpoint pivot velocity and a bounded exponential rotation.
//! This is a first-order covariant time discretization. Angular inventory and
//! boundary torque are physical, but finite-step moving-geometry angular or
//! energy conservation is not an exact theorem of this BDF1 map. Refine time.
//! Only prescribed global rigid mesh motion and a fixed Newtonian material
//! snapshot are owned here; mechanics, deformation, independent moving
//! obstacles, topology change and thermodynamic work/heat remain external.

use crate::{
    FluxSystem, MaterialBoundaryLoad, MaterialBoundaryTransport, MaterialBoundaryWrench,
    MaterialFluxError, MaterialFluxNumerics, MaterialFluxSystem, MaterialNavierStokesBoundary,
    MaterialNavierStokesStep, NavierStokesControls, OpenBoundaryPolicy,
};
use fs_material::fluid_transport::FluidTransport;

/// Physical motion of the reference pivot and the entire rigid mesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAleMotion {
    /// Velocity of the mesh pivot in world axes [m/s].
    pub linear_velocity_world_m_s: [f64; 2],
    /// Positive counterclockwise angular rate [rad/s].
    pub angular_velocity_rad_s: f64,
}
impl RigidAleMotion {
    /// Stationary physical mesh.
    pub const ZERO: Self = Self {
        linear_velocity_world_m_s: [0.0; 2],
        angular_velocity_rad_s: 0.0,
    };
    fn finite(self) -> bool {
        self.linear_velocity_world_m_s.iter().all(|v| v.is_finite())
            && self.angular_velocity_rad_s.is_finite()
    }
}

/// Absolute agreement required between requested and represented mesh rates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAleMotionTolerance {
    /// Maximum pivot-rate component error [m/s], positive and finite.
    pub linear_velocity_m_s: f64,
    /// Maximum angular-rate error [rad/s], positive and finite.
    pub angular_velocity_rad_s: f64,
}

/// Constant-work pose/rate preview, reused by the actual flow candidate.
///
/// Author moving-wall callbacks from this pivot and these admitted rates.
/// Full vertex, boundary and native-flow admission remains with candidate().
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAleKinematics {
    /// New physical reference-pivot position [m].
    pub pivot_world_m: [f64; 2],
    /// Proper reference-to-world planar rotation as [cosine, sine].
    pub orientation_cos_sin: [f64; 2],
    /// Represented rates used by the SAME native relative advector.
    pub mesh_velocity: RigidAleMotion,
    /// Signed represented rotation since the accepted state, in (-pi, pi).
    pub rotation_increment_rad: f64,
    /// New physical time [s].
    pub time_s: f64,
}

/// Accepted physical geometry, clock and original absolute-velocity history.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAleState {
    kinematics: RigidAleKinematics,
    absolute_velocity: Vec<f64>,
}
impl RigidAleState {
    /// Physical pivot position [m].
    #[must_use]
    pub const fn pivot_world_m(&self) -> [f64; 2] {
        self.kinematics.pivot_world_m
    }
    /// Current reference-to-world rotation [cosine, sine].
    #[must_use]
    pub const fn orientation_cos_sin(&self) -> [f64; 2] {
        self.kinematics.orientation_cos_sin
    }
    /// Accepted physical rigid mesh rates.
    #[must_use]
    pub const fn mesh_velocity(&self) -> RigidAleMotion {
        self.kinematics.mesh_velocity
    }
    /// Absolute fluid velocity expressed in the current mesh axes.
    ///
    /// These are the original reference BDM normal moments. They are also the
    /// physical normal moments on the rotated edges; they are NOT world-vector
    /// coefficients on the unrotated reference mesh.
    #[must_use]
    pub fn absolute_velocity_dofs(&self) -> &[f64] {
        &self.absolute_velocity
    }
    /// Accepted time [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.kinematics.time_s
    }
}

/// An invalid trial never changes accepted geometry, flow or time.
#[derive(Clone, Debug, PartialEq)]
pub enum RigidAleError {
    /// Motion, dimensions, unsupported law or arithmetic outside the domain.
    Invalid(&'static str),
    /// Original material, boundary, nonlinear or linear-solver refusal.
    Flow(MaterialFluxError),
    /// Cancellation before or after the synchronous native solve.
    Cancelled,
}
impl core::fmt::Display for RigidAleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid ALE flow refused: {self:?}")
    }
}
impl std::error::Error for RigidAleError {}
impl From<MaterialFluxError> for RigidAleError {
    fn from(value: MaterialFluxError) -> Self {
        Self::Flow(value)
    }
}

/// One solved rigid moving-domain endpoint, in coherent SI per unit depth.
pub struct RigidAleStep {
    /// Original nonlinear material solve. Its absolute velocity and vector
    /// reactions use current MESH axes; pressure and kinetic energy are physical.
    pub flow: MaterialNavierStokesStep,
    /// Fixed pivot in the unchanged reference mesh [m].
    pub reference_pivot_m: [f64; 2],
    /// Physical pivot at the endpoint [m].
    pub pivot_world_m: [f64; 2],
    /// Proper reference-to-world rotation [cosine, sine].
    pub orientation_cos_sin: [f64; 2],
    /// Actual rates used by relative convection and boundary transport.
    pub mesh_velocity: RigidAleMotion,
    /// Signed represented rotation over this interval [rad].
    pub rotation_increment_rad: f64,
    /// New physical time [s].
    pub time_s: f64,
    /// Constant source density times unchanged physical area [kg/m].
    pub mass_kg_per_m: f64,
    /// Integral rho*u_world over the physical domain [N s/m].
    pub world_momentum_n_s_per_m: [f64; 2],
    /// Integral rho*(x-X) cross u_world, per depth [(N m s)/m = N s].
    pub world_angular_momentum_about_pivot_n_s: f64,
    /// Original native kinetic quadrature, rotation invariant [J/m].
    pub world_kinetic_energy_j_per_m: f64,
}
impl RigidAleStep {
    /// Convert a reference-chart point to its physical endpoint position.
    pub fn point_world(&self, point_reference_m: [f64; 2]) -> Result<[f64; 2], RigidAleError> {
        let point = mapped_point(
            self.reference_pivot_m,
            self.pivot_world_m,
            self.orientation_cos_sin,
            point_reference_m,
        );
        if !point_reference_m
            .iter()
            .chain(&point)
            .all(|x| x.is_finite())
        {
            return Err(invalid("unrepresentable world point"));
        }
        Ok(point)
    }

    /// Actual angular inventory about an authored WORLD origin [N s per depth].
    /// This is a storage query, not a finite-step angular-conservation claim.
    pub fn world_angular_momentum(&self, origin_world_m: [f64; 2]) -> Result<f64, RigidAleError> {
        let arm = sub(self.pivot_world_m, origin_world_m);
        let angular =
            self.world_angular_momentum_about_pivot_n_s + cross(arm, self.world_momentum_n_s_per_m);
        if !origin_world_m.iter().chain(&arm).all(|v| v.is_finite()) || !angular.is_finite() {
            return Err(invalid("unrepresentable world angular momentum"));
        }
        Ok(angular)
    }

    /// Physical fluid-on-exterior stress and outward WORLD momentum.
    ///
    /// Native retained advection is relative to the MOVING control surface.
    /// Carried velocity is already absolute, so only vector axes and the moment
    /// origin transform; no mesh velocity or pressure force is added twice.
    /// Scalar planar torque and slip dissipation are rotation invariant.
    pub fn boundary_loads_world(
        &self,
        edges: &[usize],
        origin_world_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, RigidAleError> {
        let origin = self.reference_origin(origin_world_m)?;
        let mut load = self
            .flow
            .solution
            .boundary_loads(edges, origin, pressure_offsets_pa)?;
        for wrench in [
            &mut load.stress,
            &mut load.imposed,
            &mut load.backflow,
            &mut load.outward_momentum,
        ] {
            rotate_wrench(wrench, self.orientation_cos_sin)?;
        }
        load.origin_m = origin_world_m;
        Ok(load)
    }

    /// Physical mass, linear/angular momentum and carried kinetic transport
    /// through the rigid moving surface, using the SAME admitted upwind trace.
    /// Stress work, enthalpy and heat remain separately owned.
    pub fn boundary_transport_world(
        &self,
        edges: &[usize],
        origin_world_m: [f64; 2],
    ) -> Result<MaterialBoundaryTransport, RigidAleError> {
        let origin = self.reference_origin(origin_world_m)?;
        let mut transport = self
            .flow
            .solution
            .boundary_transport(edges, origin, [0.0; 2])?;
        rotate_wrench(&mut transport.outward_momentum, self.orientation_cos_sin)?;
        transport.origin_m = origin_world_m;
        Ok(transport)
    }

    fn reference_origin(&self, origin: [f64; 2]) -> Result<[f64; 2], RigidAleError> {
        let local = add(
            self.reference_pivot_m,
            inverse_rotate(self.orientation_cos_sin, sub(origin, self.pivot_world_m)),
        );
        if !origin.iter().chain(&local).all(|v| v.is_finite()) {
            return Err(invalid("unrepresentable world load origin"));
        }
        Ok(local)
    }
}

/// Original material/operator owner on one globally rigid moving mesh.
///
/// Reference bases, topology and metrics are prepared once. Each trial needs
/// only an affine mesh-advection vector and the native nonlinear workspace;
/// no deformed element preparation, global history projection or remap occurs.
pub struct RigidAleFlow<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'m>,
    fluid: MaterialFluxSystem<'s, 'm, 'a, F>,
    reference_pivot_m: [f64; 2],
    state: RigidAleState,
}

/// Solved trial borrowing its original owner until atomic acceptance.
pub struct RigidAleCandidate<'o, 's, 'm, 'a, F: FluidTransport + ?Sized> {
    owner: &'o mut RigidAleFlow<'s, 'm, 'a, F>,
    state: RigidAleState,
    step: RigidAleStep,
}
impl<F: FluidTransport + ?Sized> RigidAleCandidate<'_, '_, '_, '_, F> {
    /// Solved physical endpoint and native boundary queries.
    #[must_use]
    pub const fn step(&self) -> &RigidAleStep {
        &self.step
    }
    /// Proposed complete history; accepted owner is still unchanged.
    #[must_use]
    pub const fn state(&self) -> &RigidAleState {
        &self.state
    }
    /// Publish only after every mechanical/thermal consumer has accepted.
    /// A dropped candidate leaves the original owner unchanged.
    #[must_use]
    pub fn commit(self) -> RigidAleStep {
        self.owner.state = self.state;
        self.step
    }
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> RigidAleFlow<'s, 'm, 'a, F> {
    /// Bind a material snapshot and WORLD velocity on the initial mesh.
    /// Initially R=I, X=reference_pivot and time=0, so geometry is unchanged.
    /// Initial mesh rates are recorded, but absolute history needs no invented
    /// acceleration source when those rates change.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        reference_pivot_m: [f64; 2],
        initial_motion: RigidAleMotion,
        initial_world_velocity: &[f64],
    ) -> Result<Self, RigidAleError> {
        if initial_world_velocity.len() != system.n_u
            || !reference_pivot_m
                .iter()
                .chain(initial_world_velocity)
                .all(|v| v.is_finite())
            || !initial_motion.finite()
        {
            return Err(invalid(
                "finite pivot, motion and full physical velocity history required",
            ));
        }
        let fluid = MaterialFluxSystem::new(system, material, numerics)?;
        Ok(Self {
            system,
            fluid,
            reference_pivot_m,
            state: RigidAleState {
                kinematics: RigidAleKinematics {
                    pivot_world_m: reference_pivot_m,
                    orientation_cos_sin: [1.0, 0.0],
                    mesh_velocity: initial_motion,
                    rotation_increment_rad: 0.0,
                    time_s: 0.0,
                },
                absolute_velocity: initial_world_velocity.to_vec(),
            },
        })
    }
    /// Complete accepted state, preserved by every refusal and dropped trial.
    #[must_use]
    pub const fn state(&self) -> &RigidAleState {
        &self.state
    }
    /// The original retained source owner, without a repeated material query.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.fluid.material()
    }
    /// Fixed physical location of the pivot in the initial/reference mesh.
    #[must_use]
    pub const fn reference_pivot_m(&self) -> [f64; 2] {
        self.reference_pivot_m
    }

    /// Preview the SAME represented pose/rates used by candidate().
    ///
    /// X_new=X_old+dt*V_requested. Rotation is exp(dt*omega_requested);
    /// each requested angle must lie strictly in (-pi,pi), so the represented
    /// signed increment has an unambiguous rate. Larger turns need substeps.
    /// Normalization removes only floating multiplication drift in SO(2).
    /// Requested and represented rates must pass their authored absolute gates.
    /// This constant-work call does not invoke material or boundary callbacks.
    pub fn kinematics(
        &self,
        timestep_s: f64,
        requested: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
    ) -> Result<RigidAleKinematics, RigidAleError> {
        let dt = timestep_s;
        let time = self.state.time_s() + dt;
        let angle = dt * requested.angular_velocity_rad_s;
        if !dt.is_finite()
            || dt <= 0.0
            || !dt.recip().is_finite()
            || !time.is_finite()
            || time <= self.state.time_s()
            || !requested.finite()
            || !angle.is_finite()
            || angle.abs() >= core::f64::consts::PI
            || !positive(tolerance.linear_velocity_m_s)
            || !positive(tolerance.angular_velocity_rad_s)
        {
            return Err(invalid(
                "finite positive time/tolerances and rigid rotation increment below pi required",
            ));
        }
        let old = self.state.pivot_world_m();
        let pivot = core::array::from_fn(|i| old[i] + dt * requested.linear_velocity_world_m_s[i]);
        let linear = core::array::from_fn(|i| (pivot[i] - old[i]) / dt);
        let [c0, s0] = self.state.orientation_cos_sin();
        let (s, c) = angle.sin_cos();
        let unscaled = [c0.mul_add(c, -s0 * s), s0.mul_add(c, c0 * s)];
        let norm = unscaled[0].hypot(unscaled[1]);
        let rotation = if angle == 0.0 {
            [c0, s0]
        } else {
            unscaled.map(|v| v / norm)
        };
        let delta = relative_rotation(rotation, [c0, s0]);
        // relative_rotation returns old-to-new vector coordinates, so negate
        // its sine to recover the physical reference-to-world increment.
        let increment = (-delta[1]).atan2(delta[0]);
        let angular = increment / dt;
        if !pivot
            .iter()
            .chain(&linear)
            .chain(&rotation)
            .all(|v| v.is_finite())
            || !angular.is_finite()
            || (0..2).any(|i| {
                (linear[i] - requested.linear_velocity_world_m_s[i]).abs()
                    > tolerance.linear_velocity_m_s
            })
            || (angular - requested.angular_velocity_rad_s).abs() > tolerance.angular_velocity_rad_s
        {
            return Err(invalid(
                "rigid mesh pose and rates are not jointly representable",
            ));
        }
        Ok(RigidAleKinematics {
            pivot_world_m: pivot,
            orientation_cos_sin: rotation,
            mesh_velocity: RigidAleMotion {
                linear_velocity_world_m_s: linear,
                angular_velocity_rad_s: angular,
            },
            rotation_increment_rad: increment,
            time_s: time,
        })
    }

    /// Solve and commit one prescribed rigid-domain physical endpoint.
    /// Direct stepping retains the native relative-only correction criterion.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        requested: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RigidAleStep, RigidAleError> {
        Ok(self
            .candidate(
                timestep_s,
                requested,
                tolerance,
                force_world,
                boundaries_world,
                controls,
                0.0,
                cancelled,
            )?
            .commit())
    }

    /// Solve the original nonlinear equations without publishing any state.
    ///
    /// Callbacks receive endpoint WORLD positions/normals and provide physical
    /// force density [N/m³], velocity [m/s] and total traction [Pa].
    /// Existing velocity, Navier-slip and OpenReservoir laws remain intact.
    /// Fixed-frame DissipativeBackflow refuses because substituting relative
    /// normal flow would change that authored stress law.
    ///
    /// The absolute affine normal-velocity correction bound supplements only
    /// the native relative correction criterion; both nonlinear residual gates
    /// remain mandatory. Zero retains relative-only correction admission.
    /// Cancellation polls coordinate/storage work and surrounds the synchronous
    /// native solve; it does not interrupt its current Picard/linear internals.
    #[allow(clippy::too_many_arguments)]
    pub fn candidate(
        &mut self,
        timestep_s: f64,
        requested: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidAleCandidate<'_, 's, 'm, 'a, F>, RigidAleError> {
        poll(&mut cancelled)?;
        let motion = self.kinematics(timestep_s, requested, tolerance)?;
        if !absolute_velocity_tolerance_m_s.is_finite()
            || absolute_velocity_tolerance_m_s < 0.0
            || boundaries_world.iter().any(|b| {
                matches!(
                    b,
                    MaterialNavierStokesBoundary::Open {
                        backflow: OpenBoundaryPolicy::DissipativeBackflow,
                        ..
                    }
                )
            })
        {
            return Err(invalid(
                "finite nonnegative correction bound and literal physical boundary laws required",
            ));
        }
        let pivot_ref = self.reference_pivot_m;
        let rotation = motion.orientation_cos_sin;
        let pivot = motion.pivot_world_m;
        for vertices in self.system.mesh.verts.chunks(64) {
            poll(&mut cancelled)?;
            if vertices.iter().any(|&p| {
                !mapped_point(pivot_ref, pivot, rotation, p)
                    .iter()
                    .all(|v| v.is_finite())
            }) {
                return Err(invalid("rigid world geometry overflow"));
            }
        }
        let local_translation =
            inverse_rotate(rotation, motion.mesh_velocity.linear_velocity_world_m_s);
        let omega = motion.mesh_velocity.angular_velocity_rad_s;
        let mut mesh_velocity = vec![0.0; self.system.n_u];
        // Exact moments of the globally affine, divergence-free rigid rate.
        // No quadrature interpolation or per-cell discontinuous grid field.
        for (edge, geometry) in self.system.mesh.edges.iter().enumerate() {
            if edge % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let arm = sub(geometry.mid, pivot_ref);
            let delta = sub(
                self.system.mesh.verts[geometry.verts.1],
                self.system.mesh.verts[geometry.verts.0],
            );
            let mean = add(local_translation, [-omega * arm[1], omega * arm[0]]);
            mesh_velocity[2 * edge] = dot(mean, geometry.normal);
            mesh_velocity[2 * edge + 1] =
                dot([-omega * delta[1], omega * delta[0]], geometry.normal) / 12.0;
        }
        if mesh_velocity.iter().any(|v| !v.is_finite()) {
            return Err(invalid("rigid mesh velocity moments overflow"));
        }
        let history_rotation = relative_rotation(rotation, self.state.orientation_cos_sin());
        let world_point = |p| mapped_point(pivot_ref, pivot, rotation, p);
        let force = |p| inverse_rotate(rotation, force_world(world_point(p)));
        let velocities: Vec<_> = boundaries_world
            .iter()
            .map(|boundary| {
                move |point| match boundary {
                    MaterialNavierStokesBoundary::Velocity { value, .. }
                    | MaterialNavierStokesBoundary::NavierSlip { value, .. }
                    | MaterialNavierStokesBoundary::OpenReservoir {
                        velocity: value, ..
                    } => inverse_rotate(rotation, value(world_point(point))),
                    MaterialNavierStokesBoundary::Open { .. } => {
                        unreachable!("only velocity traces request velocity samples")
                    }
                }
            })
            .collect();
        let tractions: Vec<_> = boundaries_world
            .iter()
            .map(|boundary| {
                move |point, normal| match boundary {
                    MaterialNavierStokesBoundary::Open { traction, .. }
                    | MaterialNavierStokesBoundary::OpenReservoir { traction, .. } => {
                        inverse_rotate(
                            rotation,
                            traction(world_point(point), rotate(rotation, normal)),
                        )
                    }
                    MaterialNavierStokesBoundary::Velocity { .. }
                    | MaterialNavierStokesBoundary::NavierSlip { .. } => {
                        unreachable!("only open boundaries request traction samples")
                    }
                }
            })
            .collect();
        let boundaries: Vec<_> = boundaries_world
            .iter()
            .enumerate()
            .map(|(i, boundary)| match *boundary {
                MaterialNavierStokesBoundary::Velocity { edge, .. } => {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &velocities[i],
                    }
                }
                MaterialNavierStokesBoundary::NavierSlip {
                    edge,
                    tangential_resistance_pa_s_per_m,
                    ..
                } => MaterialNavierStokesBoundary::NavierSlip {
                    edge,
                    value: &velocities[i],
                    tangential_resistance_pa_s_per_m,
                },
                MaterialNavierStokesBoundary::Open { edge, backflow, .. } => {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &tractions[i],
                        backflow,
                    }
                }
                MaterialNavierStokesBoundary::OpenReservoir { edge, .. } => {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &tractions[i],
                        velocity: &velocities[i],
                    }
                }
            })
            .collect();
        poll(&mut cancelled)?;
        let flow = self.fluid.bdf1_navier_stokes_mixed_rigid_ale_step(
            &force,
            &boundaries,
            &self.state.absolute_velocity,
            history_rotation,
            &mesh_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
        )?;
        poll(&mut cancelled)?;
        let velocity = flow.solution.velocity_dofs();
        let rho = self.fluid.transport().density_kg_m3();
        let mut storage = Sums::<4>::new();
        for cell in 0..self.system.mesh.tris.len() {
            if cell % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let area = self.system.mesh.areas[cell];
            let vertices = self.system.mesh.tris[cell].map(|i| self.system.mesh.verts[i]);
            // All four integrands have degree at most two. These are original
            // physical field integrals, not a balance inferred from input loads.
            for (point, weight) in crate::bdm::tri_quad(vertices, area) {
                let mut value = [0.0; 2];
                for i in 0..6 {
                    let edge = self.system.mesh.tri_edges[cell][i / 2].0;
                    let coefficient = velocity[2 * edge + i % 2];
                    let basis = crate::bdm::eval_basis(&self.system.bases[cell], i, point);
                    value[0] += coefficient * basis[0];
                    value[1] += coefficient * basis[1];
                }
                let world_velocity = rotate(rotation, value);
                let arm = rotate(rotation, sub(point, pivot_ref));
                let mass = rho * weight;
                storage.add([
                    mass,
                    mass * world_velocity[0],
                    mass * world_velocity[1],
                    mass * cross(arm, world_velocity),
                ]);
            }
        }
        let [mass, px, py, angular] = storage.values();
        let kinetic = flow.solution.kinetic_energy_j_per_m();
        if !positive(mass) || ![px, py, angular, kinetic].iter().all(|v| v.is_finite()) {
            return Err(invalid("rigid world storage overflow"));
        }
        let absolute_velocity = velocity.to_vec();
        let state = RigidAleState {
            kinematics: motion,
            absolute_velocity,
        };
        let step = RigidAleStep {
            flow,
            reference_pivot_m: pivot_ref,
            pivot_world_m: pivot,
            orientation_cos_sin: rotation,
            mesh_velocity: motion.mesh_velocity,
            rotation_increment_rad: motion.rotation_increment_rad,
            time_s: motion.time_s,
            mass_kg_per_m: mass,
            world_momentum_n_s_per_m: [px, py],
            world_angular_momentum_about_pivot_n_s: angular,
            world_kinetic_energy_j_per_m: kinetic,
        };
        poll(&mut cancelled)?;
        Ok(RigidAleCandidate {
            owner: self,
            state,
            step,
        })
    }
}

fn rotate([c, s]: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    [c.mul_add(v[0], -s * v[1]), s.mul_add(v[0], c * v[1])]
}
fn inverse_rotate([c, s]: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    rotate([c, -s], v)
}
fn relative_rotation([c1, s1]: [f64; 2], [c0, s0]: [f64; 2]) -> [f64; 2] {
    if [c1, s1] == [c0, s0] {
        return [1.0, 0.0];
    }
    [c1.mul_add(c0, s1 * s0), c1.mul_add(s0, -s1 * c0)]
}
fn mapped_point(reference: [f64; 2], pivot: [f64; 2], rotation: [f64; 2], p: [f64; 2]) -> [f64; 2] {
    add(pivot, rotate(rotation, sub(p, reference)))
}
fn add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[0], a[1] * b[1])
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn invalid(message: &'static str) -> RigidAleError {
    RigidAleError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RigidAleError> {
    if cancelled() {
        Err(RigidAleError::Cancelled)
    } else {
        Ok(())
    }
}
fn rotate_wrench(
    value: &mut MaterialBoundaryWrench,
    rotation: [f64; 2],
) -> Result<(), RigidAleError> {
    value.force_n_per_m = rotate(rotation, value.force_n_per_m);
    if value.force_n_per_m.iter().any(|v| !v.is_finite()) {
        return Err(invalid("world boundary force overflow"));
    }
    Ok(())
}

struct Sums<const N: usize> {
    values: [f64; N],
    correction: [f64; N],
}
impl<const N: usize> Sums<N> {
    const fn new() -> Self {
        Self {
            values: [0.0; N],
            correction: [0.0; N],
        }
    }
    fn add(&mut self, values: [f64; N]) {
        for (i, value) in values.into_iter().enumerate() {
            let old = self.values[i];
            let next = old + value;
            self.correction[i] += if old.abs() >= value.abs() {
                (old - next) + value
            } else {
                (value - next) + old
            };
            self.values[i] = next;
        }
    }
    fn values(self) -> [f64; N] {
        core::array::from_fn(|i| self.values[i] + self.correction[i])
    }
}

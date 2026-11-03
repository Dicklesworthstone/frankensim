//! Incompressible nonlinear flow on an exactly translating physical domain.
//!
//! With x = X(t) + y and u = W(t) + v, incompressibility gives
//! dv/dt + div(v tensor v) - div(2 nu D(v)) + grad(p/rho)
//! = f_world(y + X)/rho - dW/dt. Translation has identity deformation
//! gradient and unit Jacobian: the original BDM moments, metrics and fully
//! nonlinear relative convection are used without a deformation/remap model.
//! BDF1 advances X by dt*W at the endpoint, and its same represented W enters
//! the acceleration source, boundary conversion and world momentum reports.
//!
//! All quantities are coherent SI, per metre depth where appropriate. This
//! owner prescribes mesh translation; it does not solve rigid mechanics,
//! rotation, deformation, changing topology, or a whole-domain energy audit.

use crate::{
    FluxSystem, MaterialBoundaryLoad, MaterialFluxError, MaterialFluxNumerics,
    MaterialFluxSystem, MaterialNavierStokesBoundary, MaterialNavierStokesStep,
    NavierStokesControls, OpenBoundaryPolicy,
};
use fs_material::fluid_transport::FluidTransport;

/// Accepted physical translation, clock and relative-frame velocity history.
#[derive(Clone, Debug, PartialEq)]
pub struct TranslatingAleState {
    translation_m: [f64; 2],
    mesh_velocity_m_s: [f64; 2],
    relative_velocity: Vec<f64>,
    time_s: f64,
}

impl TranslatingAleState {
    /// Physical domain is reference mesh plus this translation [m].
    #[must_use]
    pub const fn translation_m(&self) -> [f64; 2] {
        self.translation_m
    }

    /// Accepted represented endpoint mesh velocity [m/s].
    #[must_use]
    pub const fn mesh_velocity_m_s(&self) -> [f64; 2] {
        self.mesh_velocity_m_s
    }

    /// BDM1 moments of u_world - mesh_velocity, used by relative convection.
    #[must_use]
    pub fn relative_velocity_dofs(&self) -> &[f64] {
        &self.relative_velocity
    }

    /// Physical accepted time [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
}

/// A refusal leaves the complete physical state unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum TranslatingAleError {
    /// Invalid motion, dimensions, boundary law or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original nonlinear, boundary, material or linear-solver refusal.
    Flow(MaterialFluxError),
    /// Cancellation before or after the synchronous native flow solve.
    Cancelled,
}

impl core::fmt::Display for TranslatingAleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "translating ALE flow refused: {self:?}")
    }
}
impl std::error::Error for TranslatingAleError {}
impl From<MaterialFluxError> for TranslatingAleError {
    fn from(value: MaterialFluxError) -> Self {
        Self::Flow(value)
    }
}

/// One accepted moving-domain result, with explicit frame distinctions.
pub struct TranslatingAleStep {
    /// Original nonlinear solve in relative velocity coordinates. Its pressure
    /// and stress are physical; its kinetic energy and advected momentum use v.
    pub relative_flow: MaterialNavierStokesStep,
    /// BDM1 moments of physical world velocity u = v + W.
    pub world_velocity_dofs: Vec<f64>,
    /// Physical domain translation at this endpoint [m].
    pub translation_m: [f64; 2],
    /// Actual represented (X1-X0)/dt, admitted against the requested speed.
    pub mesh_velocity_m_s: [f64; 2],
    /// The acceleration used in the inertial source, (W1-W0)/dt [m/s²].
    pub mesh_acceleration_m_s2: [f64; 2],
    /// New physical time [s].
    pub time_s: f64,
    /// Constant reference density times unchanged physical area [kg/m].
    pub mass_kg_per_m: f64,
    /// Physical integral rho*u over the moving domain [N s/m].
    pub world_momentum_n_s_per_m: [f64; 2],
    /// Original kinetic-energy quadrature applied to u, rather than v [J/m].
    pub world_kinetic_energy_j_per_m: f64,
}

impl TranslatingAleStep {
    /// Physical stress loads and outward WORLD momentum across a MOVING
    /// control surface. Advective flux is rho*u_upwind*((w_world-W).n),
    /// using reservoir world velocity on incoming OpenReservoir traces.
    /// Origins and offsets are physical world metres and Pascals.
    ///
    /// The native relative-flux load is corrected by W times its SAME retained
    /// advector's mass flux and first spatial moment. Stress/pressure are invariant
    /// under translation; the moment origin is translated before evaluation.
    /// Native per-component pressure datum admission remains authoritative.
    /// No stress/advective double counting or power/passivity claim is implied.
    ///
    /// # Errors
    /// Native invalid edge/reference selections and nonfinite conversions refuse.
    pub fn boundary_loads_world(
        &self,
        edges: &[usize],
        origin_world_m: [f64; 2],
        pressure_offsets_pa: &[Option<f64>],
    ) -> Result<MaterialBoundaryLoad, TranslatingAleError> {
        let origin = [
            origin_world_m[0] - self.translation_m[0],
            origin_world_m[1] - self.translation_m[1],
        ];
        if !origin_world_m.iter().chain(&origin).all(|x| x.is_finite()) {
            return Err(TranslatingAleError::Invalid("nonfinite world load origin"));
        }
        let mut load = self.relative_flow.solution.boundary_loads(
            edges,
            origin,
            pressure_offsets_pa,
        )?;
        let (mass_flux, first_moment) = self.relative_flow.solution
            .outward_mass_flux_moments(edges, origin)?;
        for axis in 0..2 {
            load.outward_momentum.force_n_per_m[axis] +=
                self.mesh_velocity_m_s[axis] * mass_flux;
        }
        load.outward_momentum.moment_n +=
            first_moment[0] * self.mesh_velocity_m_s[1]
                - first_moment[1] * self.mesh_velocity_m_s[0];
        load.origin_m = origin_world_m;
        if !load.outward_momentum.force_n_per_m.iter().all(|v| v.is_finite())
            || !load.outward_momentum.moment_n.is_finite()
        {
            return Err(TranslatingAleError::Invalid("world momentum-flux conversion overflow"));
        }
        Ok(load)
    }
}

/// One material snapshot and unchanged reference operator on a translating
/// physical domain. Reference bases and topology are prepared once; no moved
/// mesh or new element basis is assembled during stepping.
pub struct TranslatingAleFlow<'s, 'm, 'a, F: FluidTransport + ?Sized> {
    system: &'s FluxSystem<'m>,
    fluid: MaterialFluxSystem<'s, 'm, 'a, F>,
    state: TranslatingAleState,
}

impl<'s, 'm, 'a, F: FluidTransport + ?Sized> TranslatingAleFlow<'s, 'm, 'a, F> {
    /// Bind the original material/mesh and a physical WORLD velocity history.
    /// Initial translation and time are zero; the reference mesh is the actual
    /// initial domain. A uniform frame velocity changes only the mean normal
    /// BDM moment, never the linear moment.
    ///
    /// # Errors
    /// Refuses invalid initial rates/history, source data or numerical controls.
    pub fn new(
        system: &'s FluxSystem<'m>,
        material: &'a F,
        numerics: MaterialFluxNumerics,
        initial_mesh_velocity_m_s: [f64; 2],
        initial_world_velocity: &[f64],
    ) -> Result<Self, TranslatingAleError> {
        let relative_velocity = shifted_velocity(
            system,
            initial_world_velocity,
            initial_mesh_velocity_m_s.map(|v| -v),
        )?;
        let fluid = MaterialFluxSystem::new(system, material, numerics)?;
        Ok(Self {
            system,
            fluid,
            state: TranslatingAleState {
                translation_m: [0.0; 2],
                mesh_velocity_m_s: initial_mesh_velocity_m_s,
                relative_velocity,
                time_s: 0.0,
            },
        })
    }

    /// Complete accepted physical state, unchanged by every refused trial.
    #[must_use]
    pub const fn state(&self) -> &TranslatingAleState {
        &self.state
    }

    /// Original immutable constitutive owner, retained rather than re-queried
    /// or replaced when the moving-frame pressure changes.
    #[must_use]
    pub const fn material(&self) -> &F {
        self.fluid.material()
    }

    /// Advance nonlinear NS on a translating domain, with physical endpoint
    /// forces [N/m³], world velocity [m/s] and total traction [Pa] callbacks.
    /// Every callback receives the ENDPOINT WORLD position.
    ///
    /// Backward Euler uses X1 = X0 + dt*requested_W1. The represented
    /// W1 = (X1-X0)/dt must agree within mesh_velocity_tolerance_m_s and is
    /// used everywhere thereafter, including -rho*(W1-W0)/dt. Relative
    /// convection is fully implicit through the original nonlinear solver.
    ///
    /// NaturalTraction and velocity boundaries retain their physical laws.
    /// The original DissipativeBackflow law is explicitly fixed-frame and
    /// refuses; silently applying it to relative velocity would change that
    /// boundary law. Pressure uses the original component-local gauges.
    ///
    /// Cancellation is polled around coordinate work and the synchronous
    /// native solve, not inside its current Picard/linear implementation.
    /// Geometry, frame rate, flow history and clock publish together only after
    /// all physical conversions and the final poll have succeeded.
    ///
    /// # Errors
    /// Invalid/unrepresentable motion, unsupported boundaries, original solve
    /// refusal, overflow or cancellation leaves accepted state unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        timestep_s: f64,
        requested_mesh_velocity_m_s: [f64; 2],
        mesh_velocity_tolerance_m_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TranslatingAleStep, TranslatingAleError> {
        poll(&mut cancelled)?;
        let dt = timestep_s;
        let time = self.state.time_s + dt;
        if !dt.is_finite()
            || dt <= 0.0
            || !dt.recip().is_finite()
            || !time.is_finite()
            || time <= self.state.time_s
            || !mesh_velocity_tolerance_m_s.is_finite()
            || mesh_velocity_tolerance_m_s <= 0.0
            || !requested_mesh_velocity_m_s.iter().all(|v| v.is_finite())
            || boundaries_world.iter().any(|b| matches!(
                b,
                MaterialNavierStokesBoundary::Open {
                    backflow: OpenBoundaryPolicy::DissipativeBackflow,
                    ..
                }
            ))
        {
            return Err(TranslatingAleError::Invalid(
                "finite representable timestep/motion and literal-traction or velocity boundaries required",
            ));
        }
        let translation: [f64; 2] = core::array::from_fn(|i| {
            self.state.translation_m[i] + dt * requested_mesh_velocity_m_s[i]
        });
        let velocity: [f64; 2] = core::array::from_fn(|i| {
            (translation[i] - self.state.translation_m[i]) / dt
        });
        let acceleration: [f64; 2] = core::array::from_fn(|i| {
            (velocity[i] - self.state.mesh_velocity_m_s[i]) / dt
        });
        if !translation.iter().chain(&velocity).chain(&acceleration).all(|v| v.is_finite())
            || (0..2).any(|i| {
                (velocity[i] - requested_mesh_velocity_m_s[i]).abs()
                    > mesh_velocity_tolerance_m_s
            })
        {
            return Err(TranslatingAleError::Invalid("mesh displacement and velocity are not jointly representable"));
        }
        // Convex cells and their quadrature points lie within these translated
        // vertex bounds. The metric is retained in the reference coordinate
        // chart; absolute-coordinate subtraction does not rebuild a new mesh.
        for vertices in self.system.mesh.verts.chunks(64) {
            poll(&mut cancelled)?;
            if vertices.iter().any(|point| {
                !(point[0] + translation[0]).is_finite()
                    || !(point[1] + translation[1]).is_finite()
            }) {
                return Err(TranslatingAleError::Invalid("translated world coordinates overflow"));
            }
        }
        let rho = self.fluid.transport().density_kg_m3();
        let force = |point: [f64; 2]| {
            let f = force_world([point[0] + translation[0], point[1] + translation[1]]);
            [f[0] - rho * acceleration[0], f[1] - rho * acceleration[1]]
        };
        // Homogeneous closure vectors avoid a box allocation per boundary.
        let velocities: Vec<_> = boundaries_world.iter().map(|boundary| {
            move |point: [f64; 2]| match boundary {
                MaterialNavierStokesBoundary::Velocity { value, .. }
                | MaterialNavierStokesBoundary::OpenReservoir { velocity: value, .. } => {
                    let u = value([point[0] + translation[0], point[1] + translation[1]]);
                    [u[0] - velocity[0], u[1] - velocity[1]]
                }
                MaterialNavierStokesBoundary::Open { .. } => {
                    unreachable!("only velocity edges request velocity samples")
                }
            }
        }).collect();
        let tractions: Vec<_> = boundaries_world.iter().map(|boundary| {
            move |point: [f64; 2], normal| match boundary {
                MaterialNavierStokesBoundary::Open { traction, .. }
                | MaterialNavierStokesBoundary::OpenReservoir { traction, .. } => {
                    traction([point[0] + translation[0], point[1] + translation[1]], normal)
                }
                MaterialNavierStokesBoundary::Velocity { .. } => {
                    unreachable!("only open edges request traction samples")
                }
            }
        }).collect();
        let boundaries: Vec<_> = boundaries_world.iter().enumerate().map(|(i, boundary)| {
            match *boundary {
                MaterialNavierStokesBoundary::Velocity { edge, .. } => {
                    MaterialNavierStokesBoundary::Velocity { edge, value: &velocities[i] }
                }
                MaterialNavierStokesBoundary::Open { edge, backflow, .. } => {
                    MaterialNavierStokesBoundary::Open { edge, traction: &tractions[i], backflow }
                }
                MaterialNavierStokesBoundary::OpenReservoir { edge, .. } => {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &tractions[i],
                        velocity: &velocities[i],
                    }
                }
            }
        }).collect();
        poll(&mut cancelled)?;
        let flow = self.fluid.bdf1_navier_stokes_mixed_step(
            &force,
            &boundaries,
            &self.state.relative_velocity,
            dt,
            controls,
        )?;
        poll(&mut cancelled)?;
        let world = shifted_velocity(self.system, flow.solution.velocity_dofs(), velocity)?;
        let kinetic = rho * self.system.kinetic_energy(&world);
        let mut mass = 0.0;
        let mut momentum = [0.0; 2];
        // BDM1 is affine in each cell: centroid quadrature integrates momentum
        // exactly. Energy retains the original higher-order native quadrature.
        for cell in 0..self.system.mesh.tris.len() {
            if cell % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let point = self.system.mesh.centroids[cell];
            let cell_mass = rho * self.system.mesh.areas[cell];
            mass += cell_mass;
            for i in 0..6 {
                let edge = self.system.mesh.tri_edges[cell][i / 2].0;
                let dof = 2 * edge + i % 2;
                let basis = crate::bdm::eval_basis(&self.system.bases[cell], i, point);
                for axis in 0..2 {
                    momentum[axis] += cell_mass * world[dof] * basis[axis];
                }
            }
        }
        if !mass.is_finite() || !kinetic.is_finite() || !momentum.iter().all(|v| v.is_finite()) {
            return Err(TranslatingAleError::Invalid("world storage conversion overflow"));
        }
        let relative_velocity = flow.solution.velocity_dofs().to_vec();
        poll(&mut cancelled)?;
        self.state = TranslatingAleState {
            translation_m: translation,
            mesh_velocity_m_s: velocity,
            relative_velocity,
            time_s: time,
        };
        Ok(TranslatingAleStep {
            relative_flow: flow,
            world_velocity_dofs: world,
            translation_m: translation,
            mesh_velocity_m_s: velocity,
            mesh_acceleration_m_s2: acceleration,
            time_s: time,
            mass_kg_per_m: mass,
            world_momentum_n_s_per_m: momentum,
            world_kinetic_energy_j_per_m: kinetic,
        })
    }
}

fn shifted_velocity(
    system: &FluxSystem<'_>,
    original: &[f64],
    shift: [f64; 2],
) -> Result<Vec<f64>, TranslatingAleError> {
    if original.len() != system.n_u
        || !original.iter().chain(&shift).all(|v| v.is_finite())
    {
        return Err(TranslatingAleError::Invalid("finite full velocity history and frame rate required"));
    }
    let mut shifted = original.to_vec();
    for (edge, geometry) in system.mesh.edges.iter().enumerate() {
        shifted[2 * edge] += shift[0] * geometry.normal[0] + shift[1] * geometry.normal[1];
    }
    if !shifted.iter().all(|v| v.is_finite()) {
        return Err(TranslatingAleError::Invalid("velocity-frame conversion overflow"));
    }
    Ok(shifted)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), TranslatingAleError> {
    if cancelled() {
        Err(TranslatingAleError::Cancelled)
    } else {
        Ok(())
    }
}

//! Original material inertia and units on the simultaneous physical-dry solver.

use super::super::settling::{SubmergedEllipsoid, SubmergedSphere};
use super::super::*;
use fs_lbm::d3q19::{DryBodiesExchange3, DryBodyCouplingConfig3};

mod planes;
pub use planes::{PlaneFilmSi, ResolvedImmersedPlanes};

struct StepLoads<const N: usize> {
    next: u64,
    impulses: [[f64; 6]; N],
    origin: [f64; 3],
}

/// One to eight ordinary material bodies sharing one resolved fluid state.
///
/// Light, neutral and heavy bodies retain their original positive dry mass and
/// full inertia. Reduced-pressure gravity changes only the applied load:
/// `(M_dry-rho_reference*V_original)*g`. The same reference fluid supplies bulk
/// and spherical-film viscosity. No rigid-interior approximation or excess
/// inertia is substituted for the native measured interior-momentum correction.
///
/// Construct with [`Self::spheres`] or [`Self::ellipsoids`]. Fluid starts at rest;
/// supplied body motion remains intact. All native solve, clearance, Mach,
/// quadrature and volume-change limits remain authoritative. Positive inertia
/// is not a guarantee of convergence for every density ratio or timestep.
pub struct ResolvedImmersedBodies<'a, F: FluidTransport + ?Sized, const N: usize> {
    native: RetainedFluidBodies<'a, F, N, &'a ResolvedBody>,
    volumes_m3: [f64; N],
    displaced_masses_kg: [f64; N],
    reduced_gravity_forces_n: [[f64; 3]; N],
    gravity_m_s2: [f64; 3],
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> ResolvedImmersedBodies<'a, F, N> {
    /// Bind original spherical material regions to the simultaneous dry-body
    /// and missing-film owner. Normal and optional tangent/spin film resistance
    /// participate inside every coupled candidate, never as a later impulse.
    /// Geometric film cutoffs are explicit; viscosity comes from `fluid` only.
    #[allow(clippy::too_many_arguments)]
    pub fn spheres(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [SubmergedSphere<'a>; N],
        initial: [BodyStateSi; N],
        gravity_m_s2: [f64; 3],
        film: FilmGeometry,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        Self::with_surfaces(
            fluid,
            domain,
            sources.map(SubmergedSphere::material),
            initial,
            gravity_m_s2,
            |scaling| {
                let dx = scaling.cell_size_m();
                let mut model = SphereLubrication3::new(
                    scaling.flow().viscosity,
                    film.cutoff_gap_m / dx,
                    film.roughness_gap_m / dx,
                    film.maximum_pairs,
                )?
                .with_resistance_mode_limit(film.maximum_resistance_modes)?;
                if let Some([cutoff, roughness]) = film.tangential_gaps_m {
                    model = model.with_tangential_pairs(cutoff / dx, roughness / dx)?;
                }
                Ok(CarrierSurfaces::Spheres {
                    radii: sources.map(|source| source.radius_m() / dx),
                    film: model,
                })
            },
            gate,
        )
    }

    /// Bind original ellipsoidal regions, including their complete rotated
    /// material-frame inertia and shape axes. This uses the resolved fluid law;
    /// enclosing-sphere clearance remains required and no ellipsoid film law is
    /// inferred from a representative radius.
    #[allow(clippy::too_many_arguments)]
    pub fn ellipsoids(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [SubmergedEllipsoid<'a>; N],
        initial: [BodyStateSi; N],
        gravity_m_s2: [f64; 3],
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        Self::with_surfaces(
            fluid,
            domain,
            sources.map(SubmergedEllipsoid::material),
            initial,
            gravity_m_s2,
            |scaling| {
                let mut shapes = Vec::with_capacity(N);
                for source in sources {
                    poll(gate)?;
                    shapes.push(
                        AdditionalMaterialEllipsoid::new(
                            source.material(),
                            source.semiaxes_m(),
                            source.body_from_shape(),
                        )?
                        .lattice_shape(scaling.cell_size_m())?,
                    );
                }
                Ok(CarrierSurfaces::Ellipsoids(
                    shapes
                        .try_into()
                        .map_err(|_| Error::Invalid("dry ellipsoid count"))?,
                ))
            },
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn with_surfaces(
        fluid: &'a F,
        domain: ImmersedDomain,
        materials: [&'a ResolvedBody; N],
        initial: [BodyStateSi; N],
        gravity_m_s2: [f64; 3],
        build: impl FnOnce(&MaterialScaling<'a, F>) -> Result<CarrierSurfaces<N>, Error>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        if !(1..=8).contains(&N) {
            return Err(Error::Invalid(
                "physical dry assembly requires one through eight bodies",
            ));
        }
        let gravity_m_s2 = convert(gravity_m_s2, 1.0)?;
        // The source descriptors already require exactly one matching uniform
        // region. The original source volume, never sampled occupancy, owns
        // the reduced-pressure load. Native construction retains dry inertia.
        let volumes_m3 = materials.map(|source| source.regions()[0].geometry().volume_m3());
        let native = RetainedFluidBodies::new_with_surfaces(
            fluid, domain, materials, materials, initial, build, gate,
        )?;
        let density = native.scaling.transport().density_kg_m3();
        let mut displaced_masses_kg = [0.0; N];
        let mut reduced_gravity_forces_n = [[0.0; 3]; N];
        for i in 0..N {
            poll(gate)?;
            displaced_masses_kg[i] = positive(density * positive(volumes_m3[i])?)?;
            let reduced_mass = finite(materials[i].inertia().mass() - displaced_masses_kg[i])?;
            reduced_gravity_forces_n[i] = convert(gravity_m_s2, reduced_mass)?;
        }
        Ok(Self {
            native,
            volumes_m3,
            displaced_masses_kg,
            reduced_gravity_forces_n,
            gravity_m_s2,
        })
    }

    /// Original immutable material bodies in the native stable body order.
    pub const fn materials(&self) -> &[&'a ResolvedBody; N] {
        self.native.sources()
    }
    /// One source-owned fluid snapshot and lattice/SI scaling.
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        self.native.scaling()
    }
    /// Original material volumes, not occupancy-quadrature estimates [m^3].
    pub const fn volumes_m3(&self) -> &[f64; N] {
        &self.volumes_m3
    }
    /// Reference fluid mass in each original material volume [kg].
    pub const fn displaced_masses_kg(&self) -> &[f64; N] {
        &self.displaced_masses_kg
    }
    /// Gravity held at construction; the reduced-pressure grid has no gravity source.
    pub const fn gravity_m_s2(&self) -> [f64; 3] {
        self.gravity_m_s2
    }
    /// Signed reduced weight, including exact neutral cancellation when source values agree.
    pub const fn reduced_gravity_forces_n(&self) -> &[[f64; 3]; N] {
        &self.reduced_gravity_forces_n
    }
    /// Accepted authoritative fluid; arbitrary independent population writes are excluded.
    pub const fn grid(&self) -> &BoundaryGrid3 {
        self.native.grid()
    }

    /// Prepare physical initial fluid flow [m/s] before the first accepted step.
    /// Original dry-body initial momenta and poses remain intact.
    pub fn initialize_uniform_flow_si(
        &mut self,
        velocity_m_s: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.native.initialize_uniform_flow_si(velocity_m_s, gate)
    }

    /// Hold an external fluid volume force [N/m^3] without an instantaneous
    /// momentum/stress jump. The simultaneous native dry correction owns each
    /// interior-fluid force subtraction; do not repeat it in authored solid loads.
    pub fn rebind_uniform_fluid_force_si(
        &mut self,
        force_density_n_m3: [f64; 3],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.native
            .rebind_uniform_fluid_force_si(force_density_n_m3, gate)
    }

    /// Hold a dense x-fastest spatial fluid volume force [N/m^3], added to the
    /// uniform force. All body states and accepted time remain intact; the
    /// simultaneous dry solve subtracts each body's interior source once.
    pub fn rebind_fluid_force_field_si(
        &mut self,
        force_density_n_m3: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<(), Error> {
        self.native
            .rebind_fluid_force_field_si(force_density_n_m3, gate)
    }

    /// Held external fluid volume force [N/m^3].
    pub fn uniform_fluid_force_density_n_m3(&self) -> Result<[f64; 3], Error> {
        self.native.uniform_fluid_force_density_n_m3()
    }
    /// Original dry inertia and accepted body kinematics in lattice units.
    pub const fn bodies_lattice(&self) -> &[MovingPartialSaturationBody3; N] {
        self.native.bodies_lattice()
    }
    /// Number of fully accepted physical intervals.
    pub const fn accepted_steps(&self) -> u64 {
        self.native.accepted_steps()
    }
    /// Accepted pose, velocity and original dry spin momentum in SI units.
    pub fn state_si(&self, index: usize) -> Result<BodyStateSi, Error> {
        self.native.state_si(index)
    }

    /// Refresh resolved-flow and any film viscosity from one new source while
    /// retaining all accepted mechanical/fluid history. The original owner
    /// refuses density changes, which would require conservative mass/momentum
    /// transfer and new units. No thermal storage is supplied by this update.
    pub fn rebind_fluid(&mut self, fluid: &'a F, gate: &CancelGate) -> Result<(), Error> {
        self.native.rebind_fluid(fluid, gate)
    }

    /// Advance all bodies and fluid together under authored WORLD center forces
    /// [N] and center couples [N m]. Reduced weight is added once; authored loads
    /// must not include weight or buoyancy again. Body torques are world vectors.
    ///
    /// This returns the complete original dry-assembly receipt IN LATTICE UNITS.
    /// Use explicit conversion methods for observations. Interior correction,
    /// resolved-fluid loads, films and external loads retain their native separate
    /// accounts. No nominal displaced-inertia acceleration replaces measured
    /// fluid feedback and no film loss is assigned a thermal destination.
    ///
    /// All conversions and the next step count are admitted before the native
    /// atomic publication. Any refusal preserves fluid, bodies, sources and
    /// accepted time; only reusable sampler/solver scratch may change.
    pub fn step_forces_si(
        &mut self,
        authored: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        controls: DryBodyCouplingConfig3,
        gate: &CancelGate,
    ) -> Result<DryBodiesExchange3<N>, Error> {
        let loads = self.prepare_step(authored, moment_origin_m, gate)?;
        let report = match &self.native.surfaces {
            CarrierSurfaces::Spheres { radii, film } => {
                self.native.grid.step_with_dry_lubricated_spheres(
                    &mut self.native.samplers,
                    *radii,
                    &mut self.native.bodies,
                    loads.impulses,
                    *film,
                    loads.origin,
                    controls,
                    gate,
                )?
            }
            CarrierSurfaces::Ellipsoids(shapes) => {
                self.native.grid.step_with_dry_partial_saturation_bodies(
                    &mut self.native.samplers,
                    core::array::from_fn(|i| &shapes[i] as &dyn BodyFrameShape3),
                    &mut self.native.bodies,
                    loads.impulses,
                    loads.origin,
                    controls,
                    gate,
                )?
            }
        };
        self.native.accepted_steps = loads.next;
        Ok(report)
    }

    // One conversion/load authority for sphere-only, ellipsoid and plane-film
    // consumers. No fallible physical-unit operation may follow publication.
    fn prepare_step(
        &self,
        authored: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        gate: &CancelGate,
    ) -> Result<StepLoads<N>, Error> {
        poll(gate)?;
        let next = self
            .native
            .accepted_steps
            .checked_add(1)
            .ok_or(Error::Invalid("step count overflowed"))?;
        let mut impulses = [[0.0; 6]; N];
        for i in 0..N {
            poll(gate)?;
            for axis in 0..6 {
                let load = if axis < 3 {
                    finite(authored[i][axis] + self.reduced_gravity_forces_n[i][axis])?
                } else {
                    finite(authored[i][axis])?
                };
                let unit = if axis < 3 {
                    self.native.impulse_unit
                } else {
                    self.native.angular_impulse_unit
                };
                impulses[i][axis] = finite(load * (self.native.scaling.timestep_s() / unit))?;
            }
        }
        let origin = convert(moment_origin_m, 1.0 / self.native.scaling.cell_size_m())?;
        Ok(StepLoads {
            next,
            impulses,
            origin,
        })
    }

    /// Convert native mass to kg without changing accepted state.
    pub fn mass_kg(&self, lattice: f64) -> Result<f64, Error> {
        self.native.mass_kg(lattice)
    }
    /// Convert native tensor components to kg m^2.
    pub fn inertia_kg_m2(&self, lattice: f64) -> Result<f64, Error> {
        self.native.inertia_kg_m2(lattice)
    }
    /// Convert signed native linear impulse/momentum to N s.
    pub fn impulse_ns(&self, lattice: f64) -> Result<f64, Error> {
        self.native.impulse_ns(lattice)
    }
    /// Convert signed native angular impulse/momentum to N m s.
    pub fn angular_impulse_nms(&self, lattice: f64) -> Result<f64, Error> {
        self.native.angular_impulse_nms(lattice)
    }
    /// Convert signed native energy/work to J without relabeling its physical role.
    pub fn energy_j(&self, lattice: f64) -> Result<f64, Error> {
        self.native.energy_j(lattice)
    }
}

//! One material source for resolved flow and subgrid films around moving spheres.
//!
//! The original D3Q19 owner advances populations, poses and momenta together.
//! This adapter only binds source material, physical units and geometry. Its
//! native partial-saturation model retains fluid INSIDE each sphere: the sourced
//! carrier inertia is explicitly ADDITIONAL to that fluid. It is not an ordinary
//! displaced-fluid rigid solid, and no density subtraction is inferred here.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, FaceBoundary3,
    MovingPartialSaturationBodiesExchange3, MovingPartialSaturationBody3, PartialSaturationError3,
    SphereLubrication3, TILE,
};
use fs_lbm::material::{MaterialLbmError, MaterialScaling, plan_scaling_from_material};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd_material::{MaterialInertiaError, ResolvedBody, VolumeMoments};
use fs_time::lie::gyrostat::InertiaTensor3;
use fs_time::{So3, Vec3};

/// Source, physical-unit, geometry or native step refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialImmersedError {
    /// Unsupported controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original fluid/scaling owner refused.
    Scaling(MaterialLbmError),
    /// Original material/geometry owner refused.
    Geometry(MaterialInertiaError),
    /// Original simultaneous flow/body owner refused, including cancellation.
    Native(PartialSaturationError3),
}
impl core::fmt::Display for MaterialImmersedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material immersed assembly: {self:?}")
    }
}
impl std::error::Error for MaterialImmersedError {}
impl From<MaterialLbmError> for MaterialImmersedError {
    fn from(e: MaterialLbmError) -> Self {
        Self::Scaling(e)
    }
}
impl From<MaterialInertiaError> for MaterialImmersedError {
    fn from(e: MaterialInertiaError) -> Self {
        Self::Geometry(e)
    }
}
impl From<PartialSaturationError3> for MaterialImmersedError {
    fn from(e: PartialSaturationError3) -> Self {
        Self::Native(e)
    }
}
type Error = MaterialImmersedError;

/// Fixed incompressible domain and numerical resolution. Walls are stationary;
/// paired periodic faces are selected per axis. Bodies cannot cross the halo.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedDomain {
    /// Cells on each axis; positive multiples of the native tile size.
    pub cells: [usize; 3],
    /// Each true entry pairs both faces of its axis periodically.
    pub periodic: [bool; 3],
    /// Uniform physical cell edge [m].
    pub cell_size_m: f64,
    /// One native tick [s].
    pub timestep_s: f64,
    /// Scenario speed bound for the original lattice-scaling admission. The
    /// native owner also checks actual point and fluid speeds against its Mach
    /// bound. Neither check establishes physical gas incompressibility.
    pub maximum_speed_m_s: f64,
    /// Maximum total cell count, before grid and per-body sampling allocation.
    pub maximum_cells: usize,
    /// Native midpoint volume quadrature, in 1..=8 samples along each cell axis.
    pub samples_per_axis: usize,
}

/// Explicit film geometry/calibration; viscosity is deliberately absent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilmGeometry {
    /// Gap above which the extra normal resistance is zero [m].
    pub cutoff_gap_m: f64,
    /// Positive singularity regularization, no larger than the cutoff [m].
    pub roughness_gap_m: f64,
    /// Optional reciprocal tangent/spin cutoff and regularization lengths [m].
    /// The native owner requires this cutoff <= 0.02 of each sphere radius.
    pub tangential_gaps_m: Option<[f64; 2]>,
    /// Maximum active sphere-pair films in the original simultaneous solve.
    pub maximum_pairs: usize,
    /// Maximum scalar resistance modes, including tangent/spin terms.
    pub maximum_resistance_modes: usize,
}

/// Source-backed additional inertia with the same declared spherical geometry.
/// This initial lane admits one uniform spherical region per carrier. Moment
/// matching is a consistency check, not proof of a source's geometric accuracy.
#[derive(Clone, Copy, Debug)]
pub struct AdditionalMaterialSphere<'a> {
    material: &'a ResolvedBody,
    radius_m: f64,
}
impl<'a> AdditionalMaterialSphere<'a> {
    /// Bind original material inertia to its matching spherical geometry.
    pub fn new(material: &'a ResolvedBody, radius_m: f64) -> Result<Self, Error> {
        let geometry = VolumeMoments::sphere(radius_m)?;
        if material.regions().len() != 1 || material.regions()[0].geometry() != geometry {
            return Err(Error::Invalid(
                "additional material must be one matching centered uniform spherical region",
            ));
        }
        Ok(Self { material, radius_m })
    }
    /// Original body with all density/query receipts retained.
    pub const fn material(self) -> &'a ResolvedBody {
        self.material
    }
    /// Radius of both the material region and native sampled surface [m].
    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }
}

/// Physical kinematics. Angular momentum belongs to the ADDITIONAL carrier;
/// resolved interior-fluid momentum remains in the native populations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereStateSi {
    /// World center of the spherical region [m].
    pub center_m: [f64; 3],
    /// Existing proper rotation from body coordinates to world coordinates.
    pub orientation: So3,
    /// World center velocity [m/s].
    pub velocity_m_s: [f64; 3],
    /// Additional carrier's spin momentum in body coordinates [kg m2/s].
    pub angular_momentum_body_kg_m2_s: [f64; 3],
}

/// Frozen Newtonian fluid and sourced carriers on ONE native motion solver.
///
/// The reference density maps to lattice density one. Thus the physical mass
/// unit is rho*dx^3 and mu_lattice = mu*dt/(rho*dx^2) = nu_lattice. This equality
/// is derived from the shared snapshot; dynamic and kinematic viscosity are
/// not interchangeable in general. No separate film-fluid scalar is accepted.
pub struct RetainedFluidSpheres<'a, F: FluidTransport + ?Sized, const N: usize> {
    scaling: MaterialScaling<'a, F>,
    sources: [AdditionalMaterialSphere<'a>; N],
    grid: BoundaryGrid3,
    samplers: [BodyFrameSampler3; N],
    bodies: [MovingPartialSaturationBody3; N],
    radii: [f64; N],
    film: SphereLubrication3,
    mass_unit: f64,
    inertia_unit: f64,
    impulse_unit: f64,
    angular_impulse_unit: f64,
    energy_unit: f64,
    accepted_steps: u64,
}

impl<'a, F: FluidTransport + ?Sized, const N: usize> RetainedFluidSpheres<'a, F, N> {
    /// Construct a quiescent fluid and its material-derived moving carriers.
    /// Source replacement builds a NEW model, never silently changes the mass
    /// or momentum of an existing trajectory. Density, temperature and transport
    /// stay frozen; this is not a thermal/phase or variable-density FSI closure.
    pub fn new(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [AdditionalMaterialSphere<'a>; N],
        initial: [SphereStateSi; N],
        film: FilmGeometry,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let count = domain
            .cells
            .into_iter()
            .try_fold(1usize, usize::checked_mul)
            .ok_or(Error::Invalid("grid cell count overflowed"))?;
        if !(1..=64).contains(&N)
            || !(1..=8).contains(&domain.samples_per_axis)
            || domain
                .cells
                .into_iter()
                .any(|n| n == 0 || !n.is_multiple_of(TILE))
            || count > domain.maximum_cells
            || count > isize::MAX as usize / (3 * 19 * core::mem::size_of::<f64>())
        {
            return Err(Error::Invalid(
                "invalid grid/body/quadrature allocation budget",
            ));
        }
        let scaling = plan_scaling_from_material(
            fluid,
            domain.cell_size_m,
            domain.timestep_s,
            domain.maximum_speed_m_s,
        )?;
        let dx = domain.cell_size_m;
        let dt = domain.timestep_s;
        let mass_unit = positive(scaling.transport().density_kg_m3() * dx * dx * dx)?;
        let inertia_unit = positive(mass_unit * dx * dx)?;
        let impulse_unit = positive(mass_unit * (dx / dt))?;
        let angular_impulse_unit = positive(impulse_unit * dx)?;
        let energy_unit = positive(impulse_unit * (dx / dt))?;
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
        let mut bodies = Vec::with_capacity(N);
        for (source, state) in sources.iter().zip(initial) {
            poll(gate)?;
            let mut tensor = source.material.inertia().inertia_com();
            for entry in &mut tensor.m {
                *entry = finite(*entry / inertia_unit)?;
            }
            let inertia = InertiaTensor3::new(tensor)
                .map_err(|source| Error::Native(PartialSaturationError3::Motion { source }))?;
            bodies.push(MovingPartialSaturationBody3 {
                center: convert(state.center_m, 1.0 / dx)?,
                orientation: state.orientation,
                mass: positive(source.material.inertia().mass() / mass_unit)?,
                inertia,
                velocity: convert(state.velocity_m_s, dt / dx)?,
                angular_momentum_body: vector(convert(
                    state.angular_momentum_body_kg_m2_s,
                    1.0 / angular_impulse_unit,
                )?),
            });
        }
        let boundaries = BoundarySpec3::new(core::array::from_fn(|i| {
            if domain.periodic[i / 2] {
                FaceBoundary3::Periodic
            } else {
                FaceBoundary3::stationary_wall()
            }
        }));
        let grid = BoundaryGrid3::new(
            domain.cells[0],
            domain.cells[1],
            domain.cells[2],
            scaling.flow().tau,
            [0.0; 3],
            boundaries,
        );
        let mut samplers = Vec::with_capacity(N);
        for _ in 0..N {
            poll(gate)?;
            samplers.push(BodyFrameSampler3::new(&grid, domain.samples_per_axis)?);
        }
        poll(gate)?;
        Ok(Self {
            scaling,
            sources,
            grid,
            samplers: samplers
                .try_into()
                .map_err(|_| Error::Invalid("sampler count"))?,
            bodies: bodies
                .try_into()
                .map_err(|_| Error::Invalid("body count"))?,
            radii: sources.map(|s| s.radius_m / dx),
            film: model,
            mass_unit,
            inertia_unit,
            impulse_unit,
            angular_impulse_unit,
            energy_unit,
            accepted_steps: 0,
        })
    }

    /// Original source, receipts and one frozen transport snapshot.
    pub const fn scaling(&self) -> &MaterialScaling<'a, F> {
        &self.scaling
    }
    /// Original source-backed additional carriers, in stable body order.
    pub const fn sources(&self) -> &[AdditionalMaterialSphere<'a>; N] {
        &self.sources
    }
    /// Read-only native fields prevent independent changes to viscosity,
    /// density normalization, carrier mass, or the authoritative populations.
    pub const fn grid(&self) -> &BoundaryGrid3 {
        &self.grid
    }
    /// Authoritative accepted native body states, in lattice units.
    pub const fn bodies_lattice(&self) -> &[MovingPartialSaturationBody3; N] {
        &self.bodies
    }
    /// Number of fully accepted physical intervals; refusals do not increment it.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Convert one accepted native carrier to physical units without changing it.
    pub fn state_si(&self, index: usize) -> Result<SphereStateSi, Error> {
        let body = self
            .bodies
            .get(index)
            .ok_or(Error::Invalid("body index out of range"))?;
        Ok(SphereStateSi {
            center_m: convert(body.center, self.scaling.cell_size_m())?,
            orientation: body.orientation,
            velocity_m_s: convert(
                body.velocity,
                self.scaling.cell_size_m() / self.scaling.timestep_s(),
            )?,
            angular_momentum_body_kg_m2_s: convert(
                array(body.angular_momentum_body),
                self.angular_impulse_unit,
            )?,
        })
    }

    /// Advance one physical dt with held WORLD forces [N] and center couples
    /// [N m]. Actual evolving midpoint gaps drive the original matched-film
    /// resistance: the cutoff contribution is subtracted, not a second full
    /// exterior drag law. All bodies and populations advance simultaneously.
    ///
    /// Returns the unchanged native receipt IN LATTICE UNITS; use the conversion
    /// methods below for momentum/work. Film loss has no finite heat destination
    /// in this isothermal model. Native scratch may change on refusal, but fields,
    /// body states and accepted_steps do not. No fallible operation follows the
    /// native atomic publication.
    pub fn step_forces_si(
        &mut self,
        forces_and_center_couples: [[f64; 6]; N],
        moment_origin_m: [f64; 3],
        gate: &CancelGate,
    ) -> Result<MovingPartialSaturationBodiesExchange3<N>, Error> {
        poll(gate)?;
        let next = self
            .accepted_steps
            .checked_add(1)
            .ok_or(Error::Invalid("step count overflowed"))?;
        let mut external = [[0.0; 6]; N];
        for a in 0..N {
            for axis in 0..6 {
                let unit = if axis < 3 {
                    self.impulse_unit
                } else {
                    self.angular_impulse_unit
                };
                external[a][axis] = finite(
                    forces_and_center_couples[a][axis] * (self.scaling.timestep_s() / unit),
                )?;
            }
        }
        let origin = convert(moment_origin_m, 1.0 / self.scaling.cell_size_m())?;
        let report = self.grid.step_with_lubricated_spheres(
            &mut self.samplers,
            self.radii,
            &mut self.bodies,
            external,
            self.film,
            origin,
            gate,
        )?;
        self.accepted_steps = next;
        Ok(report)
    }

    /// Convert signed native mass/change to kg, including density normalization.
    pub fn mass_kg(&self, lattice: f64) -> Result<f64, Error> {
        finite(lattice * self.mass_unit)
    }
    /// Convert a signed native inertia tensor component to kg m2.
    pub fn inertia_kg_m2(&self, lattice: f64) -> Result<f64, Error> {
        finite(lattice * self.inertia_unit)
    }
    /// Convert a signed native linear impulse/momentum component to N s.
    pub fn impulse_ns(&self, lattice: f64) -> Result<f64, Error> {
        finite(lattice * self.impulse_unit)
    }
    /// Convert a signed native angular impulse/momentum component to N m s.
    pub fn angular_impulse_nms(&self, lattice: f64) -> Result<f64, Error> {
        finite(lattice * self.angular_impulse_unit)
    }
    /// Convert a signed native work/storage/loss account to J without relabeling it.
    pub fn energy_j(&self, lattice: f64) -> Result<f64, Error> {
        finite(lattice * self.energy_unit)
    }
}

fn positive(x: f64) -> Result<f64, Error> {
    if x.is_finite() && x > 0.0 {
        Ok(x)
    } else {
        Err(Error::Invalid(
            "nonpositive or unrepresentable physical scale",
        ))
    }
}
fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Invalid("physical conversion is not finite"))
    }
}
fn convert<const N: usize>(values: [f64; N], scale: f64) -> Result<[f64; N], Error> {
    let converted = values.map(|x| x * scale);
    if converted.into_iter().all(f64::is_finite) {
        Ok(converted)
    } else {
        Err(Error::Invalid("physical vector conversion is not finite"))
    }
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(Error::Native(PartialSaturationError3::Cancelled))
    } else {
        Ok(())
    }
}

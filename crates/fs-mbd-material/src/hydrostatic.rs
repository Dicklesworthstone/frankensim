//! Material-driven hydrostatics on moving, possibly heterogeneous rigid bodies.
//!
//! Integrate displaced volume below a planar free surface, then apply
//! `F = -rho_fluid V_sub g` at its centroid. Body mass and inertia come from
//! the existing material compiler; fluid density comes from the same transport
//! snapshot used by heat/flow consumers. The two densities are never conflated.
//!
//! `HydrostaticMedium` has uniform gravity and homogeneous constant-density
//! fluid, negligible density above the surface, prescribed hydrostatic reservoir.
//! `LayeredHydrostaticMedium` also admits an overlying fluid and stable density
//! layers, integrating displaced fluid mass and its centroid across interfaces.
//! In the single-medium model gauge pressure is zero
//! on the plane; uniform exterior pressure cancels on a closed envelope. No
//! waves, drag, added mass, flooding, capillarity, compressible stratification
//! or finite-container reaction is modeled. A gas snapshot is admissible only
//! as a caller-chosen local constant-density approximation, not its full EOS.
//!
//! Reference: R. Fitzpatrick, Fluid Mechanics, Hydrostatics, sections
//! "Buoyancy" and "Equilibrium of Floating Bodies" (University of Texas).
//! <https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node27.html>

mod geometry;
mod layered;
mod step;
pub use step::{HydrostaticStep, HydrostaticStepControls, HydrostaticStepError};
pub use geometry::{DisplacedVolume, DisplacementMesh, MAX_DISPLACEMENT_CELLS};
pub use layered::{HydrostaticLayer, LayeredHydrostaticMedium, MAX_HYDROSTATIC_LAYERS};

use fs_ga::{GaError, Se3, Vec3};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_mbd::{DynamicsError, Gravity, RigidBodyState, Wrench};
use crate::MaterialInertiaError;
use crate::principal::{PrincipalBody, PrincipalFrame, PrincipalFrameError};

/// A failed evaluation returns no partial load, state or displaced volume.
#[derive(Clone, Debug, PartialEq)]
pub enum HydrostaticError {
    /// Invalid controls, geometry or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The original shared fluid owner refused its snapshot.
    Material(MaterialError),
    /// The original volume-moment owner refused a displacement cell.
    Geometry(MaterialInertiaError),
    /// The principal-frame owner refused the state/geometry transformation.
    Principal(PrincipalFrameError),
    /// The existing rigid-body dynamics owner refused an operation.
    Dynamics(DynamicsError),
    /// Cancellation observed between bounded cell operations or before return.
    Cancelled,
}
impl core::fmt::Display for HydrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material hydrostatics: {self:?}")
    }
}
impl std::error::Error for HydrostaticError {}
impl From<MaterialError> for HydrostaticError { fn from(e: MaterialError) -> Self { Self::Material(e) } }
impl From<DynamicsError> for HydrostaticError { fn from(e: DynamicsError) -> Self { Self::Dynamics(e) } }
impl From<PrincipalFrameError> for HydrostaticError { fn from(e: PrincipalFrameError) -> Self { Self::Principal(e) } }
impl From<GaError> for HydrostaticError { fn from(e: GaError) -> Self { Self::Principal(PrincipalFrameError::Frame(e)) } }

/// Conservative hydrostatic loading on an admitted material body frame.
///
/// Uniform and layered reservoirs share this boundary so thermal expansion
/// changes the same displacement geometry, fluid mass, force and pressure
/// potential consumed by ordinary rigid motion. A load excludes body weight;
/// its pressure-potential pose derivative must agree with its wrench.
/// Reservoir sources and interfaces stay frozen during an accounting window.
/// Cancellation remains generic, preserving the native per-cell callback path.
pub trait HydrostaticReservoir {
    /// Gravity supplied exactly once to the native rigid-body integrator.
    fn gravity(&self) -> Gravity;

    /// Evaluate on the original or thermally expanded principal frame.
    /// No failed source, geometry or cancellation may publish a partial load.
    fn load_frame(&self, mesh: &DisplacementMesh, body: &PrincipalFrame,
        state: RigidBodyState, cancelled: impl FnMut() -> bool)
        -> Result<HydrostaticLoad, HydrostaticError>;

    /// Reuse the native endpoint-load Verlet step and absolute energy-defect gate.
    /// Both endpoint loads use the supplied frame and actual displacement mesh.
    fn step_frame(&self, mesh: &DisplacementMesh, body: &PrincipalFrame,
        state: RigidBodyState, controls: HydrostaticStepControls,
        mut cancelled: impl FnMut() -> bool)
        -> Result<HydrostaticStep, HydrostaticStepError> {
        step::step_with_load(self.gravity(), body, state, controls, &mut cancelled,
            |state, cancel| self.load_frame(mesh, body, state, cancel))
    }
}

/// One frozen fluid evaluation, retaining its source owner. Temperature or
/// composition changes require evaluating that source and making a new medium;
/// no coefficient changes during a mechanics step or hidden water default.
#[derive(Debug)]
pub struct HydrostaticMedium<'a, F: FluidTransport + ?Sized> {
    source: &'a F,
    transport: FluidTransportState,
    gravity: Gravity,
    gravity_magnitude: f64,
    down_world: Vec3,
    surface_point_world: Vec3,
}
impl<'a, F: FluidTransport + ?Sized> HydrostaticMedium<'a, F> {
    /// Resolve the shared transport view exactly once. The horizontal free
    /// surface passes through `surface_point_world`, perpendicular to nonzero
    /// gravity. Both the gravity and the geometry use the same world frame.
    pub fn new(source: &'a F, gravity: Gravity, surface_point_world: Vec3) -> Result<Self, HydrostaticError> {
        finite_vec(surface_point_world)?;
        let (down_world, gravity_magnitude) = gravity_direction(gravity)?;
        let transport = source.transport_state()?;
        Ok(Self { source, transport, gravity, gravity_magnitude, down_world, surface_point_world })
    }
    /// Original source and its receipts, not just a copied scalar density.
    pub const fn source(&self) -> &'a F { self.source }
    /// The exact shared snapshot consumed by this frozen hydrostatic model.
    pub const fn transport(&self) -> FluidTransportState { self.transport }
    /// Gravity to pass ONCE to the existing dynamics integrator.
    pub const fn gravity(&self) -> Gravity { self.gravity }
    /// A point on the prescribed planar free surface, in world metres.
    pub const fn surface_point_world(&self) -> Vec3 { self.surface_point_world }

    /// Clip the displacement envelope at its current ORIGINAL body-frame pose.
    /// All cut cases (0..=4 wet tetrahedron vertices) integrate analytically;
    /// this is not a sampled wet-fraction estimate. Polls once per cell.
    pub fn displacement(&self, mesh: &DisplacementMesh, world_from_body: Se3,
        mut cancelled: impl FnMut() -> bool) -> Result<DisplacedVolume, HydrostaticError> {
        mesh.submerged(world_from_body, self.down_world, self.surface_point_world, &mut cancelled)
    }

    /// Buoyant force and torque in the existing free-body integrator's frames.
    /// Body mass/inertia come from `body`; displaced volume comes from `mesh`.
    /// They may differ for a sealed hollow body. The caller must bind the
    /// correct envelope to this specimen and keep the state on its principal
    /// frame; no moment/mass comparison can prove envelope coverage.
    ///
    /// The returned wrench EXCLUDES body weight: use `self.gravity()` on the
    /// integrator to avoid double counting. Pressure potential is
    /// `rho |g| integral max(depth,0) dV`; its negative pose derivative is this
    /// buoyancy wrench. It enables an independent mechanical energy audit.
    pub fn load(&self, mesh: &DisplacementMesh, body: &PrincipalBody<'_>, state: RigidBodyState,
        mut cancelled: impl FnMut() -> bool) -> Result<HydrostaticLoad, HydrostaticError> {
        self.load_frame(mesh, body.frame(), state, &mut cancelled)
    }

    /// Shared load assembly for the already admitted original or expanded frame.
    pub(crate) fn load_frame(&self, mesh: &DisplacementMesh, body: &PrincipalFrame,
        state: RigidBodyState, mut cancelled: impl FnMut() -> bool)
        -> Result<HydrostaticLoad, HydrostaticError> {
        poll(&mut cancelled)?;
        let displaced = self.displacement(mesh, body.world_from_body(state)?, &mut cancelled)?;
        let mass = finite(self.transport.density_kg_m3()*displaced.volume_m3())?;
        if displaced.volume_m3() > 0.0 && mass <= 0.0 {
            return Err(HydrostaticError::Invalid("displaced fluid mass is unrepresentable"));
        }
        let pressure_potential_j = finite(self.transport.density_kg_m3()
            * self.gravity_magnitude * displaced.depth_integral_m4)?;
        let load = HydrostaticLoad::from_mass_moments(displaced, mass,
            displaced.centroid_body_m(), pressure_potential_j, self.gravity, body, state)?;
        poll(&mut cancelled)?;
        Ok(load)
    }
}

impl<F: FluidTransport + ?Sized> HydrostaticReservoir for HydrostaticMedium<'_, F> {
    fn gravity(&self) -> Gravity { HydrostaticMedium::gravity(self) }

    fn load_frame(&self, mesh: &DisplacementMesh, body: &PrincipalFrame,
        state: RigidBodyState, cancelled: impl FnMut() -> bool)
        -> Result<HydrostaticLoad, HydrostaticError> {
        HydrostaticMedium::load_frame(self, mesh, body, state, cancelled)
    }
}

/// Immutable hydrostatic load. Its center of buoyancy is weighted by the
/// displaced fluid mass; this equals the volume centroid in a uniform medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrostaticLoad {
    displaced: DisplacedVolume,
    displaced_mass_kg: f64,
    center_world: Option<Vec3>,
    pressure_potential_j: f64,
    wrench: Wrench,
}
impl HydrostaticLoad {
    fn from_mass_moments(displaced: DisplacedVolume, mass: f64,
        center_body: Option<Vec3>, pressure_potential_j: f64, gravity: Gravity,
        body: &PrincipalFrame, state: RigidBodyState) -> Result<Self, HydrostaticError> {
        let force = gravity.acceleration_world().scale(-mass);
        finite_vec(ga(force))?;
        let (center_world, torque) = match center_body {
            None => (None, fs_mbd::Vec3::ZERO),
            Some(center) => {
                let arm = body.arm_for_body_point(center)?;
                let force_principal = state.pose().orientation().rotate_world_to_body(force);
                let torque = arm.cross(force_principal);
                finite_vec(ga(torque))?;
                (Some(ga(state.pose().point_world_from_body(arm)?)), torque)
            }
        };
        Ok(Self { displaced, displaced_mass_kg: mass, center_world, pressure_potential_j,
            wrench: Wrench { force_world: force, torque_body: torque } })
    }
    /// Integrated immersed geometry in the original body frame.
    pub const fn displaced(&self) -> DisplacedVolume { self.displaced }
    /// Actual displaced fluid mass [kg], including all declared density layers.
    pub const fn displaced_mass_kg(&self) -> f64 { self.displaced_mass_kg }
    /// World-space center of buoyancy [m], or None when completely dry.
    pub const fn center_of_buoyancy_world(&self) -> Option<Vec3> { self.center_world }
    /// Buoyancy alone, in the exact frames accepted by `RigidBodyIntegrator`.
    pub const fn wrench(&self) -> Wrench { self.wrench }
    /// Hydrostatic reservoir displacement potential [J]; zero when dry.
    /// With an overlying fluid, the reference plane fixes a signed datum.
    pub const fn pressure_potential_j(&self) -> f64 { self.pressure_potential_j }
}

fn gravity_direction(gravity: Gravity) -> Result<(Vec3, f64), HydrostaticError> {
    let g = ga(gravity.acceleration_world());
    let scale = g.x.abs().max(g.y.abs()).max(g.z.abs());
    if scale == 0.0 { return Err(HydrostaticError::Invalid("hydrostatic surface requires nonzero gravity")); }
    let scaled = g.scale(scale.recip());
    let norm = scaled.dot(scaled).sqrt();
    let magnitude = finite(scale * norm)?;
    if magnitude <= 0.0 { return Err(HydrostaticError::Invalid("unrepresentable gravity magnitude")); }
    Ok((finite_vec(scaled.scale(norm.recip()))?, magnitude))
}

fn finite(x: f64) -> Result<f64, HydrostaticError> {
    if x.is_finite() { Ok(x) } else { Err(HydrostaticError::Invalid("nonfinite hydrostatic arithmetic")) }
}
fn finite_vec(v: Vec3) -> Result<Vec3, HydrostaticError> {
    finite(v.x)?; finite(v.y)?; finite(v.z)?; Ok(v)
}
fn ga(v: fs_mbd::Vec3) -> Vec3 { Vec3::new(v.x,v.y,v.z) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HydrostaticError> {
    if cancelled() { Err(HydrostaticError::Cancelled) } else { Ok(()) }
}

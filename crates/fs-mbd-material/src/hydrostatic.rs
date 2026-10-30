//! Material-driven hydrostatics on moving, possibly heterogeneous rigid bodies.
//!
//! Integrate displaced volume below a planar free surface, then apply
//! `F = -rho_fluid V_sub g` at its centroid. Body mass and inertia come from
//! the existing material compiler; fluid density comes from the same transport
//! snapshot used by heat/flow consumers. The two densities are never conflated.
//!
//! Uniform gravity and homogeneous constant-density fluid, negligible density
//! above the surface, prescribed hydrostatic reservoir. Gauge pressure is zero
//! on the plane; uniform exterior pressure cancels on a closed envelope. No
//! waves, drag, added mass, flooding, capillarity, compressible stratification
//! or finite-container reaction is modeled. A gas snapshot is admissible only
//! as a caller-chosen local constant-density approximation, not its full EOS.
//!
//! Reference: R. Fitzpatrick, Fluid Mechanics, Hydrostatics, sections
//! "Buoyancy" and "Equilibrium of Floating Bodies" (University of Texas).
//! <https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node27.html>

mod geometry;
mod step;
pub use step::{HydrostaticStep, HydrostaticStepControls, HydrostaticStepError};
pub use geometry::{DisplacedVolume, DisplacementMesh, MAX_DISPLACEMENT_CELLS};

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
        let g = ga(gravity.acceleration_world());
        let scale = g.x.abs().max(g.y.abs()).max(g.z.abs());
        if scale == 0.0 { return Err(HydrostaticError::Invalid("hydrostatic surface requires nonzero gravity")); }
        let scaled = g.scale(scale.recip());
        let norm = scaled.dot(scaled).sqrt();
        let gravity_magnitude = finite(scale * norm)?;
        if gravity_magnitude <= 0.0 { return Err(HydrostaticError::Invalid("unrepresentable gravity magnitude")); }
        let down_world = finite_vec(scaled.scale(norm.recip()))?;
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
        let force = self.gravity.acceleration_world().scale(-mass);
        finite_vec(ga(force))?;
        let (center_world, torque) = match displaced.centroid_body_m() {
            None => (None, fs_mbd::Vec3::ZERO),
            Some(center) => {
                let arm = body.arm_for_body_point(center)?;
                let force_principal = state.pose().orientation().rotate_world_to_body(force);
                let torque = arm.cross(force_principal);
                finite_vec(ga(torque))?;
                (Some(ga(state.pose().point_world_from_body(arm)?)), torque)
            }
        };
        poll(&mut cancelled)?;
        Ok(HydrostaticLoad { displaced, center_world, pressure_potential_j,
            wrench: Wrench { force_world: force, torque_body: torque } })
    }
}

/// Immutable load: center of buoyancy is volume-weighted, not mass-weighted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrostaticLoad {
    displaced: DisplacedVolume,
    center_world: Option<Vec3>,
    pressure_potential_j: f64,
    wrench: Wrench,
}
impl HydrostaticLoad {
    /// Integrated immersed geometry in the original body frame.
    pub const fn displaced(&self) -> DisplacedVolume { self.displaced }
    /// World-space center of buoyancy [m], or None when completely dry.
    pub const fn center_of_buoyancy_world(&self) -> Option<Vec3> { self.center_world }
    /// Buoyancy alone, in the exact frames accepted by `RigidBodyIntegrator`.
    pub const fn wrench(&self) -> Wrench { self.wrench }
    /// Hydrostatic reservoir displacement potential [J]; zero when dry.
    pub const fn pressure_potential_j(&self) -> f64 { self.pressure_potential_j }
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

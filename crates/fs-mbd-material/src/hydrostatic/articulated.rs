//! Hydrostatic pressure on the original material-defined articulated mechanism.
//!
//! Every link keeps its own displaced-fluid volume and mass centroid. The
//! existing coadjoint moves each resultant from its principal COM frame to the
//! authored link origin before native articulated dynamics consumes it. Passing
//! a COM torque unchanged to that owner would lose the joint's buoyancy moment.

use super::{
    DensityProfileHydrostaticMedium, DisplacementMesh, HydrostaticError, HydrostaticLoad,
    HydrostaticMedium, HydrostaticReservoir, LayeredHydrostaticMedium,
    LinearDensityHydrostaticMedium, MAX_DISPLACEMENT_CELLS, finite, finite_vec, ga, poll,
};
use crate::articulated::ResolvedArticulated;
use crate::principal::PrincipalFrame;
use fs_ga::{Se3, Twist, Wrench};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::RigidBodyState;
use fs_mbd::articulated::{ArticulatedError, BaseState, forward_kinematics};

mod finite;
pub use finite::ArticulatedFiniteTankLoad;

mod sealed {
    pub trait Prescribed {}
}

/// A prescribed reservoir whose load and pressure potential are additive over
/// nonoverlapping displacement envelopes. This trait is sealed: the existing
/// uniform, layered and continuous-profile media satisfy this contract.
///
/// A finite tank is intentionally excluded. Its interfaces depend on ALL
/// displaced bodies and its liquid potential must be counted once; separate
/// per-link tank queries would change the inventory and duplicate its energy.
pub trait PrescribedHydrostaticReservoir: HydrostaticReservoir + sealed::Prescribed {}

impl<F: FluidTransport + ?Sized> sealed::Prescribed for HydrostaticMedium<'_, F> {}
impl<F: FluidTransport + ?Sized> PrescribedHydrostaticReservoir for HydrostaticMedium<'_, F> {}
impl<F: FluidTransport + ?Sized> sealed::Prescribed for LinearDensityHydrostaticMedium<'_, F> {}
impl<F: FluidTransport + ?Sized> PrescribedHydrostaticReservoir
    for LinearDensityHydrostaticMedium<'_, F>
{
}
impl sealed::Prescribed for LayeredHydrostaticMedium<'_> {}
impl PrescribedHydrostaticReservoir for LayeredHydrostaticMedium<'_> {}
impl sealed::Prescribed for DensityProfileHydrostaticMedium<'_> {}
impl PrescribedHydrostaticReservoir for DensityProfileHydrostaticMedium<'_> {}

/// A failed load query publishes neither partial link loads nor a changed
/// mechanism. Original source, geometry, frame and kinematic errors are retained.
#[derive(Clone, Debug, PartialEq)]
pub enum ArticulatedHydrostaticError {
    Physics(HydrostaticError),
    Kinematics(ArticulatedError),
}

impl core::fmt::Display for ArticulatedHydrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated hydrostatics: {self:?}")
    }
}
impl std::error::Error for ArticulatedHydrostaticError {}
impl From<HydrostaticError> for ArticulatedHydrostaticError {
    fn from(error: HydrostaticError) -> Self {
        Self::Physics(error)
    }
}
impl From<ArticulatedError> for ArticulatedHydrostaticError {
    fn from(error: ArticulatedError) -> Self {
        Self::Kinematics(error)
    }
}

/// Buoyancy-only loading in the exact link order and wrench convention of the
/// borrowed articulated model. No body weight is included: supply the SAME
/// reservoir gravity once to the original forward/inverse dynamics owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedHydrostaticLoad {
    /// Native `external_body_wrench`: force [N] and torque [N m], both in each
    /// authored link's axes, with torque about that link's origin (not its COM).
    pub external_body_wrenches: Vec<Wrench>,
    /// Original per-link displaced mass, geometry, buoyancy center and potential.
    /// Its individual native rigid-body wrenches retain their principal frames.
    pub link_loads: Vec<HydrostaticLoad>,
    /// Sum of the prescribed reservoir's displacement potentials [J]. Its
    /// negative generalized-coordinate gradient is the corresponding buoyant
    /// joint effort; body gravitational energy is a separate native account.
    pub pressure_potential_j: f64,
}

/// Prepared hydrostatics for a material-defined articulated mechanism.
///
/// The immutable borrowed source fixes topology, inertias and link frames.
/// Displacement envelopes are supplied separately in those same link frames:
/// a sealed hollow link excludes its outer envelope, including sealed voids,
/// while flooded voids must be omitted. Material volume is never substituted.
/// Envelopes must remain nonoverlapping at the evaluated configuration. This
/// owner does not perform contact detection or certify the mesh partition.
///
/// Principal frames are prepared once. A query performs one native kinematic
/// traversal and reuses the original clipped hydrostatic integration per link;
/// no material queries, new eigensolves, force law or dynamics solver are added.
/// The native kinematic traversal and returned vectors allocate O(links) storage.
/// Gravity, fluid density, partial immersion and stratification affect actual
/// native joint/base dynamics through the returned external wrenches.
#[derive(Debug)]
pub struct PreparedArticulatedHydrostatics<'a> {
    source: &'a ResolvedArticulated,
    displacement: &'a [DisplacementMesh],
    frames: Vec<PrincipalFrame>,
    cell_count: usize,
}

impl<'a> PreparedArticulatedHydrostatics<'a> {
    /// Bind exactly one envelope per native link and enforce a GLOBAL clipping
    /// cell budget, at most `MAX_DISPLACEMENT_CELLS`. Tensor reconstruction uses
    /// the original principal-frame tolerance contract, in (0, 1e-8].
    /// Cancellation is polled between links and by each bounded frame operation.
    pub fn new(
        source: &'a ResolvedArticulated,
        displacement: &'a [DisplacementMesh],
        max_cells: usize,
        relative_tensor_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HydrostaticError> {
        poll(&mut cancelled)?;
        if displacement.len() != source.model().link_count()
            || max_cells == 0
            || max_cells > MAX_DISPLACEMENT_CELLS
        {
            return Err(HydrostaticError::Invalid(
                "one displacement envelope per link and a bounded global cell budget required",
            ));
        }
        let mut cell_count = 0usize;
        for mesh in displacement {
            poll(&mut cancelled)?;
            cell_count = cell_count
                .checked_add(mesh.cells().len())
                .filter(|&count| count <= max_cells)
                .ok_or(HydrostaticError::Invalid(
                    "articulated displacement exceeds the global cell budget",
                ))?;
        }
        let mut frames = Vec::with_capacity(displacement.len());
        for body in source.bodies() {
            poll(&mut cancelled)?;
            frames.push(PrincipalFrame::from_inertia(
                body.inertia(),
                relative_tensor_tolerance,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            displacement,
            frames,
            cell_count,
        })
    }

    /// Original material-defined mechanism accepted by native dynamics.
    pub const fn source(&self) -> &'a ResolvedArticulated {
        self.source
    }

    /// Total displacement tetrahedra across every link, not a per-link budget.
    pub const fn cell_count(&self) -> usize {
        self.cell_count
    }

    /// Evaluate every link at its actual joint configuration and base pose.
    /// Loads depend on configuration only; velocities are neither sampled nor
    /// used to invent drag, added mass or wave radiation. Both fixed/prescribed
    /// and free-floating native solvers consume the same returned wrench array.
    ///
    /// Joint coordinates and limits are validated by the original kinematics
    /// owner. Cancellation surrounds that bounded traversal (at most the source
    /// model's 4096 links), then each link and native clipping cell. A late error
    /// or cancellation returns no load collection; the borrowed state is intact.
    pub fn load<R: PrescribedHydrostaticReservoir + ?Sized>(
        &self,
        reservoir: &R,
        world_from_base: Se3,
        joint_position: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedHydrostaticLoad, ArticulatedHydrostaticError> {
        let poses = self.link_poses(world_from_base, joint_position, &mut cancelled)?;
        let mut external_body_wrenches = Vec::with_capacity(self.frames.len());
        let mut link_loads = Vec::with_capacity(self.frames.len());
        let mut pressure_potential_j = 0.0;
        for ((frame, mesh), pose) in self.frames.iter().zip(self.displacement).zip(&poses) {
            poll(&mut cancelled)?;
            let state = frame
                .state_from_kinematics(*pose, Twist::zero())
                .map_err(HydrostaticError::from)?;
            let load = reservoir.load_frame(mesh, frame, state, &mut cancelled)?;
            let wrench = link_wrench(frame, state, load.wrench())?;
            pressure_potential_j = finite(pressure_potential_j + load.pressure_potential_j())?;
            external_body_wrenches.push(wrench);
            link_loads.push(load);
        }
        poll(&mut cancelled)?;
        Ok(ArticulatedHydrostaticLoad {
            external_body_wrenches,
            link_loads,
            pressure_potential_j,
        })
    }

    fn link_poses(
        &self,
        world_from_base: Se3,
        joint_position: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Se3>, ArticulatedHydrostaticError> {
        poll(cancelled)?;
        let kinematics = forward_kinematics(
            self.source.model(),
            BaseState::stationary(world_from_base),
            joint_position,
            &vec![0.0; self.source.model().dof_count()],
        )?;
        poll(cancelled)?;
        Ok(kinematics.world_from_link)
    }
}

/// Preserve the native authored-link origin as well as its axes.
fn link_wrench(
    frame: &PrincipalFrame,
    state: RigidBodyState,
    rigid: fs_mbd::Wrench,
) -> Result<Wrench, HydrostaticError> {
    let force_principal = state
        .pose()
        .orientation()
        .rotate_world_to_body(rigid.force_world);
    let wrench = Wrench::new(ga(rigid.torque_body), ga(force_principal))
        .transform_by(&frame.body_from_principal());
    finite_vec(wrench.force)?;
    finite_vec(wrench.torque)?;
    Ok(wrench)
}

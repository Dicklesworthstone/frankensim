//! Typed SI velocity histories on the shared D3Q19/D3Q27 reservoir owner.
//!
//! Scenario sampling, frame semantics, scaling, and complete batch preparation
//! are shared with the specialized D3Q19 adapter. The original generic owner
//! retains collision selection, Mach limits, force correction, and publication.

use super::{ScenarioVelocityError, VelocitySample3, prepare_velocity_faces3};
use fs_exec::CancelGate;
use fs_lbm::d3q19::Face3;
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{EquilibriumLattice3, GridError3};
use fs_lbm::material::MaterialScaling;
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::{BoundaryCondition, FrameTree};

/// Refused preparation or native publication, retaining the original error.
#[derive(Debug)]
pub enum Error<E> {
    /// Typed history, frame, scaling, speed bound, or cancellation refusal.
    Scenario(ScenarioVelocityError),
    /// Original grid admission error, including the collision's own Mach bound.
    Boundary(GridError3<E>),
}

impl<E: core::fmt::Debug> core::fmt::Display for Error<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario lattice velocity boundary refused: {self:?}")
    }
}

impl<E: core::fmt::Debug + 'static> std::error::Error for Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Scenario(error) => Some(error),
            Self::Boundary(error) => Some(error),
        }
    }
}

/// Apply one typed velocity history through the shared native reservoir owner.
///
/// The caller binds the scenario region to `face` and supplies the material
/// scaling used for the lattice. This uses [`rebind_velocity_faces3`] with one
/// selected slot; it performs no physics step or population initialization.
pub fn rebind_velocity_face3<L, F, const Q: usize>(
    grid: &mut OpenBoundaryGrid3<L, Q>,
    face: Face3,
    boundary: &BoundaryCondition,
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<VelocitySample3, Error<L::Error>>
where
    L: EquilibriumLattice3<Q>,
    F: FluidTransport + ?Sized,
{
    let mut boundaries = [None; 6];
    boundaries[face as usize] = Some(boundary);
    let samples = rebind_velocity_faces3(grid, boundaries, frames, time_s, scaling, gate)?;
    Ok(samples[face as usize].expect("selected velocity face has a sample"))
}

/// Atomically publish a batch of typed velocity histories on either quadrature.
///
/// Slots follow [`Face3::ALL`]; absent slots retain their current targets. All
/// selected histories use the same physical time and the exact preparation
/// path of [`super::rebind_velocity_faces3`]. Velocities remain absolute inertial
/// values expressed in their named axes; this does not add carrier motion.
///
/// One native rebind follows complete preparation. The native owner then checks
/// boundary kinds and both grid/collision speed envelopes before publication.
/// Selected velocity profiles refuse even if the dormant uniform target would
/// be unchanged; all unselected profiles retain their native ownership.
/// Invalid input or cancellation retains every
/// target, population, completed step and force.
///
/// No equilibrium, collision, geometry or scale is inferred or replaced. The
/// following native step consumes the held targets and reports its own measured
/// reservoir exchange. This update neither enforces mass-flow balance nor
/// advances physical time. Storage is fixed at six faces without field copies.
pub fn rebind_velocity_faces3<L, F, const Q: usize>(
    grid: &mut OpenBoundaryGrid3<L, Q>,
    boundaries: [Option<&BoundaryCondition>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<[Option<VelocitySample3>; 6], Error<L::Error>>
where
    L: EquilibriumLattice3<Q>,
    F: FluidTransport + ?Sized,
{
    let prepared = prepare_velocity_faces3(
        grid.boundaries(),
        boundaries,
        frames,
        time_s,
        scaling,
        gate,
        |face| grid.velocity_profile(face).is_some(),
    )
    .map_err(Error::Scenario)?;
    grid.rebind_open_boundaries(prepared.targets, gate)
        .map_err(Error::Boundary)?;
    Ok(prepared.samples)
}

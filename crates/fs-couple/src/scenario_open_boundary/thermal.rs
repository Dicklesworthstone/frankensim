//! Typed reservoir targets with original native momentum and D3Q7 heat.
//!
//! Sampling does not rebind the flow. Proposed targets enter the original
//! coupled transaction and publish with both physical fields and source epochs.

use super::{
    BoundarySample3, PreparedOpenFaces3, ScenarioBoundary3, ScenarioBoundaryError,
    admit_step_frequencies, prepare_open_faces3,
};
use fs_exec::Cx;
use fs_lbm::d3q19::{BoundarySpec3, Face3};
use fs_lbm::lattice3::EquilibriumLattice3;
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
pub use fs_lbm::lattice3::open_boundaries::thermal::OpenThermalDiffusion3;
use fs_lbm::lattice3::open_boundaries::thermal::{OpenThermalError3, OpenThermalStep3};
use fs_lbm::material::MaterialScaling;
use fs_lbm::thermal::open::ThermalOpenField3;
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::FrameTree;
use std::ops::ControlFlow;

/// Accepted scenario controls and independent native momentum/heat accounts.
#[derive(Debug, Clone)]
pub struct Step3 {
    pub samples: [Option<BoundarySample3>; 6],
    pub transport: OpenThermalStep3,
}

/// Refusal preserves targets, flow, force, heat, finite walls and time counters.
#[derive(Debug)]
pub enum Error<E> {
    Scenario(ScenarioBoundaryError),
    Native(OpenThermalError3<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for Error<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario flow/heat transaction refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for Error<E> {}

#[allow(clippy::too_many_arguments)]
fn prepare<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    cx: &Cx<'_>,
    has_velocity_profile: impl FnMut(Face3) -> bool,
    has_pressure_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedOpenFaces3, ScenarioBoundaryError> {
    let mut checkpoint = || {
        if cx.checkpoint().is_err() {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    admit_step_frequencies(boundaries, time_s, scaling.timestep_s(), &mut checkpoint)?;
    prepare_open_faces3(
        before,
        boundaries,
        frames,
        time_s,
        scaling,
        &mut checkpoint,
        has_velocity_profile,
        has_pressure_profile,
    )
}

/// Sample typed velocity/pressure controls and advance both native fields.
///
/// The original material scale converts momentum boundary data only. Heat
/// sources, reservoir temperatures, wall laws and scalar/tensor diffusion keep
/// the existing thermal carrier's explicit lattice/reference-capacity units.
/// Opening heat transport uses accepted physical velocities; the closing
/// source half-kick uses the newly admitted flow. Mechanical work is not heat.
/// Selected spatial-profile conflicts and unresolved sources refuse before
/// native work; a later thermal refusal also preserves the old control targets.
#[allow(clippy::too_many_arguments)]
pub fn step_serial<L, F, const Q: usize>(
    grid: &mut OpenBoundaryGrid3<L, Q>,
    thermal: &mut ThermalOpenField3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    cx: &Cx<'_>,
    collisions: Option<&[L::Collision]>,
    heat_per_step: &[f64],
    diffusion: OpenThermalDiffusion3<'_>,
    origin: [f64; 3],
) -> Result<Step3, Error<L::Error>>
where
    L: EquilibriumLattice3<Q>,
    F: FluidTransport + ?Sized,
{
    let prepared = prepare(
        grid.boundaries(),
        boundaries,
        frames,
        time_s,
        scaling,
        cx,
        |face| grid.velocity_profile(face).is_some(),
        |face| grid.pressure_profile(face).is_some(),
    )
    .map_err(Error::Scenario)?;
    let transport = grid
        .step_serial_with_thermal_and_open_boundaries(
            thermal,
            cx,
            prepared.targets,
            collisions,
            heat_per_step,
            diffusion,
            origin,
        )
        .map_err(Error::Native)?;
    Ok(Step3 {
        samples: prepared.samples,
        transport,
    })
}

/// The same transaction with retained spatial or authored endpoint forces.
///
/// None holds the current force without copying. Some supplies the next force
/// density in native x-fast lattice units and receives the old buffer only on
/// acceptance. Native half-force velocities and trapezoidal body impulse retain
/// their original owner. Refusal leaves the caller's endpoint buffer unchanged.
/// No force-dependent thermal source or buoyancy feedback is inferred.
#[allow(clippy::too_many_arguments)]
pub fn step_forced_serial<L, F, const Q: usize>(
    grid: &mut ForcedOpenBoundaryGrid3<L, Q>,
    thermal: &mut ThermalOpenField3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    cx: &Cx<'_>,
    collisions: Option<&[L::Collision]>,
    forces_at_end: Option<&mut Vec<[f64; 3]>>,
    heat_per_step: &[f64],
    diffusion: OpenThermalDiffusion3<'_>,
    origin: [f64; 3],
) -> Result<Step3, Error<L::Error>>
where
    L: EquilibriumLattice3<Q>,
    F: FluidTransport + ?Sized,
{
    let prepared = prepare(
        grid.boundaries(),
        boundaries,
        frames,
        time_s,
        scaling,
        cx,
        |face| grid.velocity_profile(face).is_some(),
        |face| grid.pressure_profile(face).is_some(),
    )
    .map_err(Error::Scenario)?;
    let transport = grid
        .step_serial_with_thermal_and_open_boundaries(
            thermal,
            cx,
            prepared.targets,
            collisions,
            forces_at_end,
            heat_per_step,
            diffusion,
            origin,
        )
        .map_err(Error::Native)?;
    Ok(Step3 {
        samples: prepared.samples,
        transport,
    })
}

//! One control-time transaction for typed velocity inlets and pressure outlets.
//!
//! This composes the existing scenario adapters before one native boundary
//! publication. It does not introduce another fluid solver, interpolation law,
//! mass-balance certificate, or time integrator.

use crate::scenario_pressure::{
    PressureBoundary3, PressureSample3, ScenarioPressureError, prepare_pressure_faces3,
};
use crate::scenario_velocity::{
    ScenarioVelocityError, VelocitySample3, prepare_velocity_faces3,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, Face3};
use fs_lbm::material::MaterialScaling;
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::{BoundaryCondition, FrameTree};

/// One explicit scenario-region binding on an already configured open face.
///
/// A slot has one owner: pressure and velocity cannot both prescribe that face.
/// The native owner retains its original boundary kind and geometry.
#[derive(Debug, Clone, Copy)]
pub enum ScenarioBoundary3<'a> {
    /// Absolute inertial velocity, expressed in the condition's named axes.
    Velocity(&'a BoundaryCondition),
    /// Scalar pressure with an explicitly declared gauge or absolute reference.
    Pressure(PressureBoundary3<'a>),
}

/// The SI/lattice sample actually installed for a selected face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundarySample3 {
    /// Original velocity adapter's sample, including its WORLD rotation.
    Velocity(VelocitySample3),
    /// Original pressure adapter's sample, including its pressure reference.
    Pressure(PressureSample3),
}

/// Refused preparation or native publication; no selected target is changed.
#[derive(Debug)]
pub enum ScenarioBoundaryError {
    /// Original velocity sampling, scaling or selected-profile refusal.
    Velocity(ScenarioVelocityError),
    /// Original pressure sampling, scaling or selected-profile refusal.
    Pressure(ScenarioPressureError),
    /// Original native boundary admission or final cancellation refusal.
    Boundary(BoundaryForcingError3),
}

impl core::fmt::Display for ScenarioBoundaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario open-boundary transaction refused: {self:?}")
    }
}

impl std::error::Error for ScenarioBoundaryError {}

/// Sample all selected inlet/outlet controls and publish them together.
///
/// Slots follow Face3::ALL. None preserves a face's complete existing control,
/// including any spatial profile. All selected histories use the same physical
/// time and supplied material scaling; explicit pressure references and moving
/// velocity frames retain their original adapter semantics.
///
/// Both preparations finish before the sole native rebind. A later invalid
/// pressure, velocity, incompatible native face, selected spatial profile or
/// cancellation cannot partially advance an earlier target. An empty selection
/// still requires finite time and no cancellation.
///
/// Populations, forces and completed steps are unchanged by this operation.
/// The next native physics step consumes these held targets. Callers still own
/// control-time scheduling, region-to-face geometry and the physical regime.
pub fn rebind_open_faces3<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<[Option<BoundarySample3>; 6], ScenarioBoundaryError> {
    let prepared = prepare_open_faces3(
        grid.boundaries(),
        boundaries,
        frames,
        time_s,
        scaling,
        gate,
        |face| grid.velocity_profile(face).is_some(),
        |face| grid.pressure_profile(face).is_some(),
    )?;
    grid.rebind_open_boundaries(prepared.targets, gate)
        .map_err(ScenarioBoundaryError::Boundary)?;
    Ok(prepared.samples)
}

struct PreparedOpenFaces3 {
    targets: BoundarySpec3,
    samples: [Option<BoundarySample3>; 6],
}

fn prepare_open_faces3<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
    has_velocity_profile: impl FnMut(Face3) -> bool,
    has_pressure_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedOpenFaces3, ScenarioBoundaryError> {
    let velocities = boundaries.map(|boundary| match boundary {
        Some(ScenarioBoundary3::Velocity(value)) => Some(value),
        _ => None,
    });
    let pressures = boundaries.map(|boundary| match boundary {
        Some(ScenarioBoundary3::Pressure(value)) => Some(value),
        _ => None,
    });
    let velocity = prepare_velocity_faces3(
        before,
        velocities,
        frames,
        time_s,
        scaling,
        gate,
        has_velocity_profile,
    )
    .map_err(ScenarioBoundaryError::Velocity)?;
    let pressure = prepare_pressure_faces3(
        velocity.targets,
        pressures,
        frames,
        time_s,
        scaling,
        gate,
        has_pressure_profile,
    )
    .map_err(ScenarioBoundaryError::Pressure)?;
    let samples = std::array::from_fn(|index| {
        velocity.samples[index]
            .map(BoundarySample3::Velocity)
            .or(pressure.samples[index].map(BoundarySample3::Pressure))
    });
    Ok(PreparedOpenFaces3 {
        targets: pressure.targets,
        samples,
    })
}

/// The same mixed transaction on the shared native D3Q19/D3Q27 owner.
pub mod lattice3 {
    use super::{
        BoundarySample3, ScenarioBoundary3, ScenarioBoundaryError, prepare_open_faces3,
    };
    use fs_exec::CancelGate;
    use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
    use fs_lbm::lattice3::{EquilibriumLattice3, GridError3};
    use fs_lbm::material::MaterialScaling;
    use fs_material::fluid_transport::FluidTransport;
    use fs_scenario::FrameTree;

    /// Preserve both adapter diagnostics and the selected collision's admission.
    #[derive(Debug)]
    pub enum Error<E> {
        /// Typed scenario, material-scaling or cancellation refusal.
        Scenario(ScenarioBoundaryError),
        /// Native boundary-kind, grid/collision envelope or cancellation refusal.
        Boundary(GridError3<E>),
    }

    impl<E: core::fmt::Debug> core::fmt::Display for Error<E> {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            write!(f, "scenario lattice open-boundary transaction refused: {self:?}")
        }
    }

    impl<E: core::fmt::Debug + 'static> std::error::Error for Error<E> {}

    /// Prepare both control families before one native admission/publication.
    ///
    /// The grid's own equilibrium, collision Mach bound, force correction and
    /// accepted state remain authoritative. Other semantics match the
    /// specialized owner's rebind_open_faces3.
    pub fn rebind_open_faces3<L, F, const Q: usize>(
        grid: &mut OpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        gate: &CancelGate,
    ) -> Result<[Option<BoundarySample3>; 6], Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let prepared = prepare_open_faces3(
            grid.boundaries(),
            boundaries,
            frames,
            time_s,
            scaling,
            gate,
            |_| false,
            |face| grid.pressure_profile(face).is_some(),
        )
        .map_err(Error::Scenario)?;
        grid.rebind_open_boundaries(prepared.targets, gate)
            .map_err(Error::Boundary)?;
        Ok(prepared.samples)
    }
}

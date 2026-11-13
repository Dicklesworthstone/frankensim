//! Typed velocity inlets and pressure outlets on the original fluid owners.
//!
//! This composes the existing scenario adapters before one native boundary
//! publication. The generic lattice owner also accepts these controls in the
//! same transaction as collision, streaming and measured reservoir transport.
//! Interpolation, physical admission and fluid stepping retain their owners.

/// Scenario reservoirs in the native joint flow/heat transaction.
pub mod thermal;

use crate::scenario_pressure::{
    PressureBoundary3, PressureSample3, ScenarioPressureError,
    prepare_pressure_faces3_with_checkpoint,
};
use crate::scenario_velocity::{
    ScenarioVelocityError, VelocitySample3, prepare_velocity_faces3_with_checkpoint,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, Face3};
use fs_lbm::material::MaterialScaling;
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::signal::DEFAULT_MAXIMUM_CYCLES_PER_STEP;
use fs_scenario::{BcValue, BoundaryCondition, FrameTree, TimeSignal};
use std::ops::ControlFlow;

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
    /// The fixed material timestep cannot resolve a selected pressure source.
    /// Rebuild the material scaling and grid with a smaller physical timestep.
    FrequencyResolution {
        /// Face carrying the unresolved pressure history.
        face: Face3,
        /// Caller-selected physical sample time [s].
        time_s: f64,
        /// Frozen physical duration of one native lattice step [s].
        timestep_s: f64,
    },
}

impl core::fmt::Display for ScenarioBoundaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if let Self::FrequencyResolution {
            face,
            time_s,
            timestep_s,
        } = self
        {
            return write!(
                f,
                "pressure source on {face:?} at {time_s} s is unresolved by the fixed lattice timestep {timestep_s} s; rebuild material scaling and the grid with a smaller timestep"
            );
        }
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
        &mut || gate_checkpoint(gate),
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

fn gate_checkpoint(gate: &CancelGate) -> ControlFlow<()> {
    if gate.is_requested() {
        ControlFlow::Break(())
    } else {
        ControlFlow::Continue(())
    }
}

fn admit_step_frequencies(
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    time_s: f64,
    timestep_s: f64,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ScenarioBoundaryError> {
    for (face, boundary) in Face3::ALL.into_iter().zip(boundaries) {
        let Some(ScenarioBoundary3::Pressure(input)) = boundary else {
            continue;
        };
        let Some(BcValue::Signal(TimeSignal::Chirp(chirp))) = &input.boundary.value else {
            continue;
        };
        if checkpoint().is_break() {
            return Err(ScenarioBoundaryError::Pressure(
                ScenarioPressureError::Cancelled,
            ));
        }
        let admitted = chirp
            .admits_sampled_step(time_s, timestep_s, DEFAULT_MAXIMUM_CYCLES_PER_STEP)
            .map_err(|error| {
                ScenarioBoundaryError::Pressure(ScenarioPressureError::Scenario(error))
            })?;
        if !admitted {
            return Err(ScenarioBoundaryError::FrequencyResolution {
                face,
                time_s,
                timestep_s,
            });
        }
    }
    Ok(())
}

fn prepare_open_faces3<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<ScenarioBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
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
    let velocity = prepare_velocity_faces3_with_checkpoint(
        before,
        velocities,
        frames,
        time_s,
        scaling,
        checkpoint,
        has_velocity_profile,
    )
    .map_err(ScenarioBoundaryError::Velocity)?;
    let pressure = prepare_pressure_faces3_with_checkpoint(
        velocity.targets,
        pressures,
        frames,
        time_s,
        scaling,
        checkpoint,
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
        BoundarySample3, ScenarioBoundary3, ScenarioBoundaryError, admit_step_frequencies,
        gate_checkpoint, prepare_open_faces3,
    };
    use fs_exec::{CancelGate, Cx, KernelRunner};
    use fs_lbm::d3q19::OpenBoundaryMomentumTransport3;
    use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
    use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
    use fs_lbm::lattice3::{EquilibriumLattice3, GridError3};
    use fs_lbm::material::MaterialScaling;
    use fs_material::fluid_transport::FluidTransport;
    use fs_scenario::FrameTree;
    use std::ops::ControlFlow;

    /// Accepted physical controls and their original native reservoir exchange.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BoundaryStep3 {
        /// Samples used by this step, in the original Face3::ALL ordering.
        pub samples: [Option<BoundarySample3>; 6],
        /// Measured population mass, momentum and angular momentum transport
        /// in lattice units; this is not a thermal or pressure-work balance.
        pub transport: OpenBoundaryMomentumTransport3,
    }

    /// Preserve both adapter diagnostics and the selected collision's admission.
    /// Refusal preserves all accepted targets, populations and completed steps.
    #[derive(Debug)]
    pub enum Error<E> {
        /// Typed scenario, material-scaling or cancellation refusal.
        Scenario(ScenarioBoundaryError),
        /// Native boundary-kind, grid/collision envelope or cancellation refusal.
        Boundary(GridError3<E>),
    }

    impl<E: core::fmt::Debug> core::fmt::Display for Error<E> {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            if let Self::Scenario(error) = self {
                return write!(
                    f,
                    "scenario lattice open-boundary transaction refused: {error}"
                );
            }
            write!(
                f,
                "scenario lattice open-boundary transaction refused: {self:?}"
            )
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
            &mut || gate_checkpoint(gate),
            |face| grid.velocity_profile(face).is_some(),
            |face| grid.pressure_profile(face).is_some(),
        )
        .map_err(Error::Scenario)?;
        grid.rebind_open_boundaries(prepared.targets, gate)
            .map_err(Error::Boundary)?;
        Ok(prepared.samples)
    }

    /// Sample the selected controls and advance one native serial fluid step.
    ///
    /// The explicit sample time may be an interval endpoint or midpoint, as
    /// chosen by the caller's time integrator. Both histories share that time,
    /// material scale and cancellation context. None preserves the complete
    /// held control, including an unselected spatial pressure profile; a
    /// selected profiled face refuses even if its uniform target is unchanged.
    ///
    /// Targets, populations and completed-step count publish only after the
    /// original collision, streaming and reservoir measurement all succeed.
    /// The returned samples therefore describe the accepted physical step.
    /// Native pressure chirps must fit the source cycle limit over any interval
    /// containing the selected sample. An unresolved source requires a smaller
    /// material-scaling timestep; this operation never changes lattice scaling.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial<L, F, const Q: usize>(
        grid: &mut OpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<BoundaryStep3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let mut checkpoint = || {
            if cx.checkpoint().is_err() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        admit_step_frequencies(boundaries, time_s, scaling.timestep_s(), &mut checkpoint)
            .map_err(Error::Scenario)?;
        let prepared = prepare_open_faces3(
            grid.boundaries(),
            boundaries,
            frames,
            time_s,
            scaling,
            &mut checkpoint,
            |face| grid.velocity_profile(face).is_some(),
            |face| grid.pressure_profile(face).is_some(),
        )
        .map_err(Error::Scenario)?;
        let transport = grid
            .step_serial_with_open_boundaries(cx, prepared.targets, collisions, origin)
            .map_err(Error::Boundary)?;
        Ok(BoundaryStep3 {
            samples: prepared.samples,
            transport,
        })
    }

    /// Executor counterpart of [`step_serial`], sharing one cancellation gate
    /// across scenario sampling and the original tiled fluid transaction.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled<L, F, P, const Q: usize>(
        grid: &mut OpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        pool: &P,
        gate: &CancelGate,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<BoundaryStep3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
        P: KernelRunner,
    {
        let mut checkpoint = || gate_checkpoint(gate);
        admit_step_frequencies(boundaries, time_s, scaling.timestep_s(), &mut checkpoint)
            .map_err(Error::Scenario)?;
        let prepared = prepare_open_faces3(
            grid.boundaries(),
            boundaries,
            frames,
            time_s,
            scaling,
            &mut checkpoint,
            |face| grid.velocity_profile(face).is_some(),
            |face| grid.pressure_profile(face).is_some(),
        )
        .map_err(Error::Scenario)?;
        let transport = grid
            .step_pooled_with_open_boundaries(pool, gate, prepared.targets, collisions, origin)
            .map_err(Error::Boundary)?;
        Ok(BoundaryStep3 {
            samples: prepared.samples,
            transport,
        })
    }

    fn prepare_forced_step<L, F, const Q: usize>(
        grid: &ForcedOpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        cx: &Cx<'_>,
    ) -> Result<super::PreparedOpenFaces3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let mut checkpoint = || {
            if cx.checkpoint().is_err() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        admit_step_frequencies(boundaries, time_s, scaling.timestep_s(), &mut checkpoint)
            .map_err(Error::Scenario)?;
        prepare_open_faces3(
            grid.boundaries(),
            boundaries,
            frames,
            time_s,
            scaling,
            &mut checkpoint,
            |face| grid.velocity_profile(face).is_some(),
            |face| grid.pressure_profile(face).is_some(),
        )
        .map_err(Error::Scenario)
    }

    /// Advance typed reservoirs with the native owner's retained spatial force.
    ///
    /// The complete force field stays attached to its populations. Scenario
    /// sampling, selected-profile refusal and pressure-source resolution use
    /// the same owners as the uniform-force path. New targets, populations,
    /// wall exchange and step count publish together after the original native
    /// force-aware reconstruction and reservoir measurement succeed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_forced_serial<L, F, const Q: usize>(
        grid: &mut ForcedOpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        origin: [f64; 3],
    ) -> Result<BoundaryStep3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let prepared = prepare_forced_step(grid, boundaries, frames, time_s, scaling, cx)?;
        let transport = grid
            .step_serial_with_open_boundaries(cx, prepared.targets, collisions, origin)
            .map_err(Error::Boundary)?;
        Ok(BoundaryStep3 {
            samples: prepared.samples,
            transport,
        })
    }

    /// Advance typed reservoirs and an explicit native endpoint force field.
    ///
    /// Force densities use the original grid's lattice units and x-fast order;
    /// conversion or material-dependent source evaluation remains with their
    /// existing owner. Collision uses the opening force, while reconstructed
    /// donors and targets use their own endpoint forces. The native trapezoidal
    /// body impulse is separate from the returned reservoir transport. A late
    /// refusal preserves every accepted source, control and fluid epoch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_forced_serial<L, F, const Q: usize>(
        grid: &mut ForcedOpenBoundaryGrid3<L, Q>,
        boundaries: [Option<ScenarioBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        forces_at_end_lu: Vec<[f64; 3]>,
        origin: [f64; 3],
    ) -> Result<BoundaryStep3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let prepared = prepare_forced_step(grid, boundaries, frames, time_s, scaling, cx)?;
        let transport = grid
            .advance_serial_with_open_boundaries(
                cx,
                prepared.targets,
                collisions,
                forces_at_end_lu,
                origin,
            )
            .map_err(Error::Boundary)?;
        Ok(BoundaryStep3 {
            samples: prepared.samples,
            transport,
        })
    }
}

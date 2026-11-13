//! Typed scenario pressure histories on the original open-boundary owners.
//!
//! Pressure is a scalar: named frames are checked but never rotate its value.
//! The caller explicitly declares gauge data or an absolute reference. Both
//! use the original material scaling's incompressible pressure perturbation;
//! no physical gas equation of state, finite reservoir or pressure-work model
//! is inferred. The caller binds each scenario region to a native lattice face.

/// Source-owned compressible reservoir exchange from pressure/temperature histories.
#[cfg(any(feature = "experimental-shock-loading", feature = "nonlinear-brass"))]
pub mod compressible;

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3};
use fs_lbm::material::{MaterialLbmError, MaterialScaling};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Pressure;
use fs_scenario::scalar::ScalarEvaluationError;
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameId, FrameTree, Physics, ScenarioError};
use std::ops::ControlFlow;

/// How the declared scalar pressure relates to the lattice's reference density.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureReference {
    /// The declared value is already a signed pressure perturbation [Pa].
    Gauge,
    /// The declared value is positive absolute pressure [Pa]. Subtract this
    /// explicitly supplied positive reference before converting to density.
    /// The caller owns its correspondence to the material reference state.
    Absolute {
        /// Absolute pressure [Pa] represented by native lattice density one.
        reference_pa: f64,
    },
}

/// One borrowed scenario boundary and its explicit pressure reference.
#[derive(Debug, Clone, Copy)]
pub struct PressureBoundary3<'a> {
    /// An IncompressibleFlow/PressureOutlet boundary on a named region.
    pub boundary: &'a BoundaryCondition,
    /// No gauge or absolute convention is guessed from the magnitude or units.
    pub reference: PressureReference,
}

/// The physical sample and the native target produced by the original scaling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PressureSample3 {
    /// Physical time at which the history and declared frame were evaluated [s].
    pub time_s: f64,
    /// Value in the caller's declared gauge or absolute convention [Pa].
    pub supplied_pa: f64,
    /// Signed perturbation about the reference-one lattice density [Pa].
    pub perturbation_pa: f64,
    /// Positive native density target, in the grid's existing lattice units.
    pub density_lattice: f64,
}

/// No native target is published after a sampling or admission refusal.
#[derive(Debug)]
pub enum ScenarioPressureError {
    /// Unsupported role, value carrier, reference, frame size or physical value.
    Invalid(&'static str),
    /// Original signal, dimensions or frame validation refused the sample.
    Scenario(ScenarioError),
    /// Original typed scalar source, physical semantics or frame refused.
    Sample(ScalarEvaluationError),
    /// The original SI pressure-to-density mapping refused the value.
    Scaling(MaterialLbmError),
    /// Original specialized-grid boundary admission refused publication.
    Boundary(BoundaryForcingError3),
    /// Caller cancellation observed before publication.
    Cancelled,
}

impl core::fmt::Display for ScenarioPressureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario pressure boundary refused: {self:?}")
    }
}
impl std::error::Error for ScenarioPressureError {}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ScenarioPressureError> {
    if checkpoint().is_break() {
        Err(ScenarioPressureError::Cancelled)
    } else {
        Ok(())
    }
}

/// Sample a legacy scalar/signal or canonical typed scalar in SI pressure units.
///
/// Signal interpolation, endpoint clamping and complete history validation
/// remain owned by TimeSignal. Validation and Chebyshev evaluation use that
/// owner's checkpoints. Frame evaluation retains the existing indivisible
/// world_pose operation and the velocity adapter's 256-frame-row limit.
/// A scalar does not acquire dynamic, hydrostatic or carrier-motion corrections.
///
/// Typed scalars use the original scalar evaluator, including explicit table
/// tails and static/instantaneous Pressure semantics. Stress and acoustic
/// amplitude contracts refuse despite sharing Pa dimensions. Spatial profiles
/// need a geometry/parameter binding and refuse here. No source history or
/// native population is copied.
pub fn sample_pressure3<F: FluidTransport + ?Sized>(
    input: PressureBoundary3<'_>,
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<PressureSample3, ScenarioPressureError> {
    use ScenarioPressureError as E;
    poll(checkpoint)?;
    if !time_s.is_finite() {
        return Err(E::Invalid("pressure sample time must be finite"));
    }
    let boundary = input.boundary;
    if boundary.physics != Physics::IncompressibleFlow || boundary.kind != BcKind::PressureOutlet {
        return Err(E::Invalid("pressure history requires IncompressibleFlow/PressureOutlet"));
    }
    if boundary.region.is_empty() || boundary.compatibility.is_some() {
        return Err(E::Invalid("pressure history requires a named patch and no mass-flow compatibility tag"));
    }
    if let PressureReference::Absolute { reference_pa } = input.reference
        && (!reference_pa.is_finite() || reference_pa <= 0.0)
    {
        return Err(E::Invalid("absolute pressure reference must be finite and positive"));
    }
    if frames.frames.len() > 256 {
        return Err(E::Invalid("pressure sampling admits at most 256 frame rows"));
    }
    frames
        .world_pose(FrameId(boundary.frame), time_s)
        .map_err(E::Scenario)?;
    poll(checkpoint)?;
    let value = match &boundary.value {
        Some(BcValue::Uniform(value)) => *value,
        Some(BcValue::Signal(signal)) => signal
            .eval_with_checkpoint(time_s, &mut |_| {
                if checkpoint().is_break() { Err(()) } else { Ok(()) }
            })
            .map_err(|()| E::Cancelled)?
            .map_err(E::Scenario)?,
        Some(BcValue::Profile(_)) => {
            return Err(E::Invalid("spatial pressure profiles require an explicit geometry binding"));
        }
        Some(BcValue::Typed(_)) => {
            let sample = boundary.pressure_at(frames, time_s, checkpoint)
                .map_err(|error| match error {
                    ScalarEvaluationError::Cancelled => E::Cancelled,
                    error => E::Sample(error),
                })?;
            fs_qty::QtyAny::new(sample.value, Pressure::DIMS)
        }
        None => return Err(E::Invalid("pressure outlet has no declared value")),
    };
    if value.dims != Pressure::DIMS {
        return Err(E::Scenario(ScenarioError::Dimensions {
            context: format!("pressure outlet on {:?}", boundary.region),
            expected: Pressure::DIMS.0,
            got: value.dims.0,
        }));
    }
    if !value.value.is_finite() {
        return Err(E::Invalid("supplied pressure must be finite"));
    }
    let perturbation_pa = match input.reference {
        PressureReference::Gauge => value.value,
        PressureReference::Absolute { reference_pa } => {
            if value.value <= 0.0 {
                return Err(E::Invalid("supplied absolute pressure must be positive"));
            }
            value.value - reference_pa
        }
    };
    let density_lattice = scaling
        .density_lattice_from_pressure_perturbation(perturbation_pa)
        .map_err(E::Scaling)?;
    poll(checkpoint)?;
    Ok(PressureSample3 {
        time_s,
        supplied_pa: value.value,
        perturbation_pa,
        density_lattice,
    })
}

pub(crate) struct PreparedPressureFaces3 {
    pub(crate) targets: BoundarySpec3,
    pub(crate) samples: [Option<PressureSample3>; 6],
}

/// Shared preparation never mutates a grid. Each selected face may use its
/// own explicit input convention; every result uses the same material scale.
pub(crate) fn prepare_pressure_faces3<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<PressureBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
    has_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedPressureFaces3, ScenarioPressureError> {
    prepare_pressure_faces3_with_checkpoint(
        before,
        boundaries,
        frames,
        time_s,
        scaling,
        &mut || {
            if gate.is_requested() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        has_profile,
    )
}

/// Share sampling with native steps using the caller's original context.
pub(crate) fn prepare_pressure_faces3_with_checkpoint<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<PressureBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    mut has_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedPressureFaces3, ScenarioPressureError> {
    poll(checkpoint)?;
    if !time_s.is_finite() {
        return Err(ScenarioPressureError::Invalid("pressure sample time must be finite"));
    }
    let mut targets = Face3::ALL.map(|face| before.face(face));
    let mut samples = [None; 6];
    for (index, input) in boundaries.into_iter().enumerate() {
        let Some(input) = input else { continue; };
        let sample = sample_pressure3(input, frames, time_s, scaling, checkpoint)?;
        // Even an unchanged dormant uniform target must not claim to have
        // applied a history while an installed spatial profile still wins.
        if has_profile(Face3::ALL[index]) {
            return Err(ScenarioPressureError::Boundary(BoundaryForcingError3::Invalid(
                "clear the selected pressure profile before applying a uniform scenario target",
            )));
        }
        targets[index] = FaceBoundary3::Pressure { density: sample.density_lattice };
        samples[index] = Some(sample);
    }
    poll(checkpoint)?;
    Ok(PreparedPressureFaces3 { targets: BoundarySpec3::new(targets), samples })
}

/// Atomically update selected existing pressure faces on the specialized D3Q19 owner.
///
/// Slots follow Face3::ALL. All selected histories sample at the same physical
/// time before one original native rebind. None retains the face's complete
/// control, including spatial profiles. A selected installed profile refuses,
/// even if its dormant uniform density equals the sampled density.
///
/// Invalid histories, references, scaling, native face kinds or cancellation
/// preserve all targets, populations, forcing and time. This updates held
/// controls only; the next native step owns actual reservoir transport.
pub fn rebind_pressure_faces3<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    boundaries: [Option<PressureBoundary3<'_>>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<[Option<PressureSample3>; 6], ScenarioPressureError> {
    let prepared = prepare_pressure_faces3(
        grid.boundaries(), boundaries, frames, time_s, scaling, gate,
        |face| grid.pressure_profile(face).is_some(),
    )?;
    grid.rebind_open_boundaries(prepared.targets, gate)
        .map_err(ScenarioPressureError::Boundary)?;
    Ok(prepared.samples)
}

/// Single-face form of rebind_pressure_faces3 with identical admission.
pub fn rebind_pressure_face3<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    face: Face3,
    input: PressureBoundary3<'_>,
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<PressureSample3, ScenarioPressureError> {
    let mut boundaries = [None; 6];
    boundaries[face as usize] = Some(input);
    let samples = rebind_pressure_faces3(grid, boundaries, frames, time_s, scaling, gate)?;
    Ok(samples[face as usize].expect("selected pressure face has a sample"))
}

/// The same pressure histories on the shared D3Q19/D3Q27 reservoir owner.
pub mod lattice3 {
    use super::{PressureBoundary3, PressureSample3, ScenarioPressureError, prepare_pressure_faces3};
    use fs_exec::CancelGate;
    use fs_lbm::d3q19::Face3;
    use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
    use fs_lbm::lattice3::{EquilibriumLattice3, GridError3};
    use fs_lbm::material::MaterialScaling;
    use fs_material::fluid_transport::FluidTransport;
    use fs_scenario::FrameTree;

    /// Refused preparation or original generic-grid publication.
    #[derive(Debug)]
    pub enum Error<E> {
        /// Original scenario/scaling refusal.
        Scenario(ScenarioPressureError),
        /// Original grid admission error.
        Boundary(GridError3<E>),
    }
    impl<E: core::fmt::Debug> core::fmt::Display for Error<E> {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            write!(f, "scenario lattice pressure boundary refused: {self:?}")
        }
    }
    impl<E: core::fmt::Debug + 'static> std::error::Error for Error<E> {}

    /// Prepare every selected pressure history, then publish once through the
    /// original generic owner. Its collision, force correction, face-kind
    /// restrictions, reservoir accounting and physical time stay authoritative.
    /// A selected spatial pressure profile refuses even when the sampled
    /// density equals its dormant uniform target. Unselected profiles remain
    /// authoritative and unchanged, as on the specialized D3Q19 owner.
    /// No equilibrium reconstruction or pressure mapping is duplicated here.
    pub fn rebind_pressure_faces3<L, F, const Q: usize>(
        grid: &mut OpenBoundaryGrid3<L, Q>,
        boundaries: [Option<PressureBoundary3<'_>>; 6],
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        gate: &CancelGate,
    ) -> Result<[Option<PressureSample3>; 6], Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let prepared = prepare_pressure_faces3(
            grid.boundaries(), boundaries, frames, time_s, scaling, gate,
            |face| grid.pressure_profile(face).is_some(),
        ).map_err(Error::Scenario)?;
        grid.rebind_open_boundaries(prepared.targets, gate).map_err(Error::Boundary)?;
        Ok(prepared.samples)
    }

    /// Single-face form of the same generic native pressure update.
    pub fn rebind_pressure_face3<L, F, const Q: usize>(
        grid: &mut OpenBoundaryGrid3<L, Q>,
        face: Face3,
        input: PressureBoundary3<'_>,
        frames: &FrameTree,
        time_s: f64,
        scaling: &MaterialScaling<'_, F>,
        gate: &CancelGate,
    ) -> Result<PressureSample3, Error<L::Error>>
    where
        L: EquilibriumLattice3<Q>,
        F: FluidTransport + ?Sized,
    {
        let mut boundaries = [None; 6];
        boundaries[face as usize] = Some(input);
        let samples = rebind_pressure_faces3(grid, boundaries, frames, time_s, scaling, gate)?;
        Ok(samples[face as usize].expect("selected pressure face has a sample"))
    }
}

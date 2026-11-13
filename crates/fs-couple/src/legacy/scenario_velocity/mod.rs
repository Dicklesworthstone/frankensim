//! Typed scenario velocities on the native three-dimensional open-boundary owners.
//!
//! The caller explicitly maps the scenario region to a lattice face. Lattice
//! axes are WORLD Cartesian axes, and velocities are absolute inertial values
//! expressed in the scenario's named axes. Moving coordinate frames rotate the
//! components; they do not supply an unrequested carrier velocity. This adapter
//! does not move geometry or infer a spatial profile from a uniform vector.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryForcingError3, BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3};
use fs_lbm::material::{MaterialLbmError, MaterialScaling};
use fs_material::fluid_transport::FluidTransport;
use fs_scenario::vector::VectorEvaluationError;
use fs_scenario::{BoundaryCondition, FrameTree};
use std::ops::ControlFlow;

/// Typed velocity controls for the shared native D3Q19/D3Q27 grid.
pub mod lattice3;

/// One checked time sample, retaining coherent SI and native lattice units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocitySample3 {
    /// Physical sample time [s], including moving-frame evaluation.
    pub time_s: f64,
    /// Signed inertial velocity in WORLD Cartesian components [m/s].
    pub world_m_s: [f64; 3],
    /// The same vector in the supplied material scaling's lattice units.
    pub lattice: [f64; 3],
}

/// A refused sample or target update leaves all native flow state unchanged.
#[derive(Debug)]
pub enum ScenarioVelocityError {
    /// Invalid boundary role, vector metadata, time, frame or cancellation.
    Scenario(VectorEvaluationError),
    /// Physical-to-lattice conversion failed.
    Scaling(MaterialLbmError),
    /// The full vector exceeds the scenario speed admitted by the scaling.
    SpeedBound {
        /// Measured Euclidean speed in lattice units.
        actual: f64,
        /// Previously admitted maximum in those same units.
        maximum: f64,
    },
    /// The original boundary owner refused the target or cancellation.
    Boundary(BoundaryForcingError3),
}

impl core::fmt::Display for ScenarioVelocityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scenario velocity boundary refused: {self:?}")
    }
}
impl std::error::Error for ScenarioVelocityError {}

/// Sample a typed velocity history, rotate it once, and convert all components.
///
/// The full Euclidean speed must fit the existing material scaling's scenario
/// bound; checking each component independently would admit a diagonal vector
/// up to sqrt(3) times too fast. Uses no allocation or material re-evaluation.
/// Legacy directionless scalar conditions cannot execute through this path.
pub fn sample_velocity3<F: FluidTransport + ?Sized>(
    boundary: &BoundaryCondition,
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<VelocitySample3, ScenarioVelocityError> {
    let world_m_s = boundary
        .velocity_world_at(frames, time_s, checkpoint)
        .map_err(ScenarioVelocityError::Scenario)?;
    let mut lattice = [0.0; 3];
    for (component, value) in lattice.iter_mut().zip(world_m_s) {
        *component = scaling
            .velocity_lattice(value)
            .map_err(ScenarioVelocityError::Scaling)?;
    }
    let actual = lattice[0].hypot(lattice[1]).hypot(lattice[2]);
    let maximum = scaling.flow().u_lattice;
    if actual > maximum {
        return Err(ScenarioVelocityError::SpeedBound { actual, maximum });
    }
    if checkpoint().is_break() {
        return Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Cancelled,
        ));
    }
    Ok(VelocitySample3 {
        time_s,
        world_m_s,
        lattice,
    })
}

/// Apply a sampled velocity to one existing open velocity face transactionally.
///
/// Sampling, unit conversion and all native admission complete before a target
/// is published. Other faces, populations, force frames and physical time are
/// retained. The next native step applies the held target using its original
/// force-corrected reconstruction. Call at each desired control time; this is
/// not a second integrator or an implicit continuous-time boundary scheme.
///
/// An installed velocity profile on the selected face refuses even when the
/// sample equals that face's stored uniform target. The original owner refuses
/// wall/pressure/periodic conversion. The explicit `face` is the caller's
/// geometry binding for `boundary.region`; no region-to-face match is inferred.
pub fn rebind_velocity_face3<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    face: Face3,
    boundary: &BoundaryCondition,
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<VelocitySample3, ScenarioVelocityError> {
    let mut boundaries = [None; 6];
    boundaries[face as usize] = Some(boundary);
    let samples = rebind_velocity_faces3(grid, boundaries, frames, time_s, scaling, gate)?;
    // Every selected slot is sampled before the batch can publish successfully.
    Ok(samples[face as usize].expect("selected velocity face has a sample"))
}

/// Sample and publish several typed velocity histories as one native update.
///
/// Input and output slots follow [`Face3::ALL`]: x-min, x-max, y-min, y-max,
/// z-min, z-max. A selected slot explicitly binds that boundary's region to
/// its face. `None` retains the face's complete existing control, including any
/// installed spatial profile. Each selected slot is sampled at the same
/// physical time, rotated in its declared frame, and admitted against the
/// original material scaling's full-vector speed bound.
///
/// Sampling and profile checks finish for every selected face before one call
/// to the original open-boundary owner publishes any target. A later invalid
/// history, bound violation, incompatible face, or cancellation therefore
/// cannot leave an earlier face at a new time. Selected spatial profiles refuse
/// even if the new uniform target equals their dormant uniform value.
///
/// This uses fixed six-slot storage, without copying histories or populations.
/// It does not advance time or enforce a mass-flow balance between faces. The
/// native restriction to open faces on one axis remains authoritative. Empty
/// selections preserve all controls but still require finite time and no
/// cancellation. Other semantics are those of [`rebind_velocity_face3`].
pub fn rebind_velocity_faces3<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    boundaries: [Option<&BoundaryCondition>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
) -> Result<[Option<VelocitySample3>; 6], ScenarioVelocityError> {
    let prepared = prepare_velocity_faces3(
        grid.boundaries(),
        boundaries,
        frames,
        time_s,
        scaling,
        gate,
        |face| grid.velocity_profile(face).is_some(),
    )?;
    grid.rebind_open_boundaries(prepared.targets, gate)
        .map_err(ScenarioVelocityError::Boundary)?;
    Ok(prepared.samples)
}

pub(crate) struct PreparedVelocityFaces3 {
    pub(crate) targets: BoundarySpec3,
    pub(crate) samples: [Option<VelocitySample3>; 6],
}

/// Both native owners share this preparation; only profile ownership and final
/// native admission differ. Nothing is published by this fixed-size operation.
pub(crate) fn prepare_velocity_faces3<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<&BoundaryCondition>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    gate: &CancelGate,
    has_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedVelocityFaces3, ScenarioVelocityError> {
    prepare_velocity_faces3_with_checkpoint(
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
pub(crate) fn prepare_velocity_faces3_with_checkpoint<F: FluidTransport + ?Sized>(
    before: BoundarySpec3,
    boundaries: [Option<&BoundaryCondition>; 6],
    frames: &FrameTree,
    time_s: f64,
    scaling: &MaterialScaling<'_, F>,
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    mut has_profile: impl FnMut(Face3) -> bool,
) -> Result<PreparedVelocityFaces3, ScenarioVelocityError> {
    if checkpoint().is_break() {
        return Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::Cancelled,
        ));
    }
    if !time_s.is_finite() {
        return Err(ScenarioVelocityError::Scenario(
            VectorEvaluationError::InvalidTime,
        ));
    }
    let mut targets = Face3::ALL.map(|face| before.face(face));
    let mut samples = [None; 6];
    for (index, boundary) in boundaries.into_iter().enumerate() {
        let Some(boundary) = boundary else {
            continue;
        };
        let sample = sample_velocity3(boundary, frames, time_s, scaling, checkpoint)?;
        // Unchanged native profiles remain authoritative. A selected scenario
        // target must never report success while an installed profile hides it.
        if has_profile(Face3::ALL[index]) {
            return Err(ScenarioVelocityError::Boundary(
                BoundaryForcingError3::Invalid(
                    "clear the selected velocity profile before applying a uniform scenario target",
                ),
            ));
        }
        targets[index] = FaceBoundary3::Velocity {
            velocity: sample.lattice,
        };
        samples[index] = Some(sample);
    }
    Ok(PreparedVelocityFaces3 {
        targets: BoundarySpec3::new(targets),
        samples,
    })
}

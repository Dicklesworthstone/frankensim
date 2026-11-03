//! Original spatial loads and finite attachments on a prescribed base.
//!
//! The force, attachment, storage and loss laws are shared with the free-flight
//! path. Only the mechanical boundary changes: native prescribed-base ABA and
//! RNEA retain the actual support reaction and its signed work.

/// Bounded implicit midpoint for stiff supported connections and smooth feedback.
pub mod implicit;

use super::super::prescribed::{
    EffortLaw, EvaluatedPrescribedStep, PrescribedBaseMotionError, PrescribedBaseMotionState,
    PrescribedBaseMotionStep, PrescribedBasePath, PrescribedStage, external_sources,
    step_evaluated,
};
use super::super::{FreeFloatingBaseState, MomentumMap, poll, validate_slice};
use super::bushing::{self, BushingLaw, SpatialBushing};
use super::{SpatialMotionLoads, separation, storage};
use crate::articulated::{ArticulatedModel, Kinematics};

/// Accepted supported motion and original finite-attachment energy accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedSpatialMotionStep {
    /// Native supported trajectory, including actuator, external, internal and
    /// support work. Connection forces enter both dynamics and support recovery.
    pub motion: PrescribedBaseMotionStep,
    /// Original axial and bushing spring storage at the supplied state, in joules.
    pub stored_energy_before: f64,
    /// The same storage at the accepted endpoint, in joules.
    pub stored_energy_after: f64,
    /// Nonnegative midpoint damping-loss quadrature, in joules.
    /// No thermal receiver or continuous-time passivity is inferred.
    pub dissipated_energy: f64,
    /// Internal connection work + actual storage change + damping loss.
    /// This discretization defect is not an additional source of heat.
    pub connection_work_residual: f64,
    /// Kinetic + connection storage change + loss, minus actuator, external
    /// and prescribed-support work. This is a numerical balance diagnostic.
    pub mechanical_work_residual: f64,
}

/// Advance a supported mechanism with original world loads and axial links.
///
/// Point forces retain their world direction at both actual stages. Finite
/// spring/damper connections use the existing material-point geometry and
/// relative rates, including prescribed base rotation. Their wrenches enter
/// the original ABA and RNEA together, so the support reacts to the complete
/// mechanism. Their storage and loss never become external forcing.
///
/// Uses the same explicit midpoint step as the held-load supported entry.
/// The duration must resolve both joint dynamics and connection modes; this
/// is not a stiff connector solver, contact law or constraint projection.
/// The caller supplies a physically consistent base path. All state and load
/// inputs are borrowed; refusal or cancellation returns no partial endpoint.
pub fn step_prescribed_base_spatial(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    step_prescribed_base_bushings(model, state, path, loads, &[], duration, cancelled)
}

/// Advance original spatial loads with joint efforts evaluated at both stages.
///
/// The effort slice is reset to the held generalized forces before each
/// callback. The callback receives elapsed time and the actual prescribed-base
/// stage, with its absolute clock, joint positions and rates. Its forces or
/// torques therefore act on the same motion as every spatial attachment.
///
/// See the supported held-load effort entry for smoothness, units and callback
/// cancellation conventions. Caller-owned callback side effects cannot be
/// rolled back, although no native state is published on any refusal.
pub fn step_prescribed_base_spatial_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    step_prescribed_base_bushings_with_effort(
        model,
        state,
        path,
        loads,
        &[],
        duration,
        effort,
        cancelled,
    )
}

/// Advance original spatial loads and finite six-axis bushings on a support.
///
/// The existing diagonal anisotropic law owns translational and principal-log
/// rotational strain, the exact SO(3) force pullback, stored energy and damping.
/// All internal wrenches enter the same native dynamics and support reaction;
/// their work, storage and loss remain separate from external and support work.
/// Empty bushings retain the ordinary supported-spatial trajectory exactly.
///
/// The explicit midpoint step must resolve the fastest attachment mode. Active
/// rotational chart admission is checked at start, midpoint and endpoint; these
/// sampled checks do not enclose the continuous path or permit crossing pi.
/// No rigid constraint, implicit solve or thermal destination is introduced.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    step_with_bushings_and_effort(
        model, state, path, loads, bushings, duration, None, cancelled,
    )
}

/// Advance the same supported bushings with stage-dependent joint effort.
///
/// The callback sees the actual complete stage and its clock, with its output
/// reset to held effort. Actuation and all attachment laws use the same stage;
/// each retains its own work/storage/loss account. The original coefficient, effort,
/// geometry, chart and cancellation refusals precede endpoint publication.
/// Callback-owned side effects cannot be rolled back.
#[allow(clippy::too_many_arguments)]
pub fn step_prescribed_base_bushings_with_effort(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[SpatialBushing],
    duration: f64,
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    step_with_bushings_and_effort(
        model,
        state,
        path,
        loads,
        bushings,
        duration,
        Some(&mut effort),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn step_with_bushings_and_effort<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    path: PrescribedBasePath,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    duration: f64,
    effort: Option<&mut EffortLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    poll(&mut cancelled)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        path,
        loads.held,
        duration,
        &mut cancelled,
        |stage, kinematics, cancelled| {
            evaluate_shared(model, stage, kinematics, loads, bushings, cancelled)
        },
        effort,
    )?;
    finish(model, evaluated, loads, bushings, &mut cancelled)
}

fn evaluate_shared<'a, B: BushingLaw>(
    model: &ArticulatedModel,
    stage: &PrescribedBaseMotionState,
    kinematics: &Kinematics,
    loads: SpatialMotionLoads<'a>,
    bushings: &[B],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedStage<'a>, PrescribedBaseMotionError> {
    // Only the original FK-based geometry map is needed. Do not reconstruct a
    // free base, solve a floating inertia, or require a fixed root joint.
    let map = MomentumMap::from_joint_state(
        model,
        &stage.joint_position,
        &stage.joint_velocity,
        cancelled,
    )?;
    let base = FreeFloatingBaseState::new(stage.base.world_from_base, stage.base.twist_body);
    let evaluated = bushing::evaluate(&map, base, loads, bushings, cancelled)?;
    // Preserve original source loads; internal connection rounding must never
    // become an external impulse in the supported motion receipt.
    let mut external = external_sources(
        model,
        stage.base.world_from_base,
        kinematics,
        &map.shape,
        loads.held,
        cancelled,
    )?;
    for (i, applied) in loads.world_points.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        external.add_world_wrench(
            applied.force_world,
            applied.couple_world,
            &[
                stage.base.world_from_base.translation(),
                stage
                    .base
                    .world_from_base
                    .rotation()
                    .rotate(map.shape.world_from_link[applied.point.link].translation())?,
                kinematics.world_from_link[applied.point.link]
                    .rotation()
                    .rotate(applied.point.point_body)?,
            ],
        );
    }
    Ok(PrescribedStage {
        body_wrench: evaluated.body_wrench,
        external,
        external_power: evaluated.load.body_wrench_power,
        internal_power: evaluated.load.internal_body_power,
        stored_energy: evaluated.load.stored_energy,
        dissipation_power: evaluated.load.dissipation_power,
    })
}

fn finish<B: BushingLaw>(
    model: &ArticulatedModel,
    evaluated: EvaluatedPrescribedStep,
    loads: SpatialMotionLoads<'_>,
    bushings: &[B],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    let mut stored_energy_after = 0.0;
    if !loads.connections.is_empty() || !bushings.is_empty() {
        let state = &evaluated.step.state_after;
        let map = MomentumMap::from_joint_state(
            model,
            &state.joint_position,
            &state.joint_velocity,
            cancelled,
        )?;
        for (index, connection) in loads.connections.iter().enumerate() {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            let (_, length) = separation(&map, *connection)?;
            stored_energy_after += storage(*connection, length);
        }
        stored_energy_after += bushing::stored_energy(&map, bushings, cancelled)?;
    }
    let connection_work_residual = evaluated.step.internal_body_work + stored_energy_after
        - evaluated.stored_energy_before
        + evaluated.dissipated_energy;
    let mechanical_work_residual = evaluated.step.kinetic_work_residual + connection_work_residual;
    validate_slice(
        &[
            stored_energy_after,
            connection_work_residual,
            mechanical_work_residual,
        ],
        "prescribed_spatial_motion.final_energy",
    )?;
    poll(cancelled)?;
    Ok(PrescribedSpatialMotionStep {
        motion: evaluated.step,
        stored_energy_before: evaluated.stored_energy_before,
        stored_energy_after,
        dissipated_energy: evaluated.dissipated_energy,
        connection_work_residual,
        mechanical_work_residual,
    })
}

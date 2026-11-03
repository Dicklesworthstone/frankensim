//! Original spatial loads and finite axial connections on a prescribed base.
//!
//! The force, attachment, storage and loss laws are shared with the free-flight
//! path. Only the mechanical boundary changes: native prescribed-base ABA and
//! RNEA retain the actual support reaction and its signed work.

use super::super::prescribed::{
    EvaluatedPrescribedStep, PrescribedBaseMotionError, PrescribedBaseMotionState,
    PrescribedBaseMotionStep, PrescribedBasePath, PrescribedStage, step_evaluated,
};
use super::super::{FreeFloatingBaseState, MomentumMap, poll, validate_slice};
use super::{SpatialMotionLoads, admit, evaluate, separation, storage};
use crate::articulated::ArticulatedModel;

/// Accepted supported motion and original finite-attachment energy accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedSpatialMotionStep {
    /// Native supported trajectory, including actuator, external, internal and
    /// support work. Connection forces enter both dynamics and support recovery.
    pub motion: PrescribedBaseMotionStep,
    /// Original axial spring storage at the supplied state, in joules.
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
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    poll(&mut cancelled)?;
    admit(model, loads, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        path,
        loads.held,
        duration,
        &mut cancelled,
        |stage, _, cancelled| evaluate_shared(model, stage, loads, cancelled),
        None,
    )?;
    finish(model, evaluated, loads, &mut cancelled)
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
    mut effort: impl FnMut(
        f64,
        &PrescribedBaseMotionState,
        &mut [f64],
    ) -> Result<(), PrescribedBaseMotionError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    poll(&mut cancelled)?;
    admit(model, loads, &mut cancelled)?;
    let evaluated = step_evaluated(
        model,
        state,
        path,
        loads.held,
        duration,
        &mut cancelled,
        |stage, _, cancelled| evaluate_shared(model, stage, loads, cancelled),
        Some(&mut effort),
    )?;
    finish(model, evaluated, loads, &mut cancelled)
}

fn evaluate_shared<'a>(
    model: &ArticulatedModel,
    stage: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'a>,
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
    let evaluated = evaluate(&map, base, loads, cancelled)?;
    Ok(PrescribedStage {
        body_wrench: evaluated.body_wrench,
        external_power: evaluated.load.body_wrench_power,
        internal_power: evaluated.load.internal_body_power,
        stored_energy: evaluated.load.stored_energy,
        dissipation_power: evaluated.load.dissipation_power,
    })
}

fn finish(
    model: &ArticulatedModel,
    evaluated: EvaluatedPrescribedStep,
    loads: SpatialMotionLoads<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedSpatialMotionStep, PrescribedBaseMotionError> {
    let mut stored_energy_after = 0.0;
    if !loads.connections.is_empty() {
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

//! Actual spatial force moments on the original finite mechanism graph.

use super::super::super::free_motion::SpatialFlow;
use super::{
    ArticulatedGraphStep, ArticulatedJoint, DistanceGraphControls, DistanceGraphError,
    DistanceStepError, DrivenJoint, DrivenJointGraphStep, DrivenJointStep, JointReactionImpulse,
    JointRows, MAX_DISTANCE_ROWS, MassProperties, RigidBodyIntegrator, RigidBodyState, Row, Rows,
    advance_between_with, driven_joints, limit, planned_bytes_for_motion, poll, reserved, revolute,
};
use crate::spatial_loads::{PreparedSpatialLoads, SpatialStepReceipt};

/// Advance a finite rod/ball/hinge mechanism under prepared spatial loads.
///
/// Every trial uses the original native spatial-load integrator. World-held
/// forces remain attached to their rotating material arms; body/world couples
/// retain their declared frames. The complete native source receipt is kept
/// for each free step, separate from the original joint reaction account.
///
/// Force/couple cross responses, nonlinear trials, endpoint impulses and rank
/// tests use the SAME graph owner as [`super::advance_articulated_graph`].
/// This adds one admitted torque-bound scalar per body, prepared once from the
/// constant-size load representation. Trial buffers are reused and no force
/// set is traversed during a trial. The original held-wrench entry keeps its
/// original free-motion map, receipt type, arithmetic and payload.
///
/// The loads are held for one interval. This is not a state-dependent force
/// callback, motor law, contact integrator or exact energy certificate. Native
/// source and reaction defects remain visible. Local travel bounds can require
/// timestep refinement, including near cancellation of large source moments.
/// Cancellation or any late refusal publishes no body or partial mechanism.
///
/// # Errors
/// Original graph/source shape, resource, rank, local-branch, numerical and
/// cancellation refusals; unrepresentable prepared torque bounds also refuse.
#[allow(clippy::too_many_arguments)]
pub fn advance_spatial_articulated_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    integrator: RigidBodyIntegrator,
    loads: &[PreparedSpatialLoads],
    duration_s: f64,
    controls: DistanceGraphControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedGraphStep<SpatialStepReceipt>, DistanceGraphError> {
    poll(&mut cancelled)?;
    if states.len() != properties.len() || states.len() != loads.len() {
        return Err(DistanceStepError::Invalid("spatial mechanism body/load shapes").into());
    }
    limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    let maximum = controls.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut rows = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut row_count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let width = joint.row_count();
        limit("constraints", (row_count + width) as u128, maximum as u128)?;
        joint.write_rows(&mut rows[row_count..row_count + width]);
        row_count += width;
    }
    let extra = joints.len() as u128 * core::mem::size_of::<JointRows>() as u128
        + loads.len() as u128 * core::mem::size_of::<f64>() as u128;
    let bytes = planned_bytes_for_motion::<JointReactionImpulse, SpatialStepReceipt>(
        states.len(),
        row_count,
    ) + extra;
    let bounds = prepare_bounds(loads, bytes, controls, &mut cancelled)?;
    let mut native_controls = controls;
    native_controls.maximum_live_bytes -= extra as usize;
    let (step, _) = advance_between_with(
        states,
        properties,
        Rows::Mixed(&rows[..row_count]),
        integrator,
        SpatialFlow {
            loads,
            bounds: &bounds,
        },
        duration_s,
        native_controls,
        None,
        &mut cancelled,
    )?;
    revolute::finish(step, joints, bytes, &mut cancelled)
}

/// Advance the same spatial loads through a finite actuated mechanism graph.
///
/// This composes [`super::advance_driven_joint_graph`]'s original quadratic
/// command and work owner with [`advance_spatial_articulated_graph`]'s native
/// spatial free motion. Source forces/couples, radial/angular actuator work and
/// passive constraint defects remain distinct in the returned receipts. Every row is
/// solved simultaneously, with finite recoil for every body.
///
/// The extra live payload is one prepared torque bound per body. The original
/// target bank and continuation vector are shared with the held-wrench driven
/// path, with no duplicate solve, command preparation or per-trial allocation.
/// Held loads, positive quadratic ranges, resolved local motion and the same
/// smooth-joint limitations of both composing owners remain required.
///
/// # Errors
/// Original graph, prescribed-command and spatial-load errors. Cancellation or
/// any late numerical/resource refusal returns no partial state or continuation.
#[allow(clippy::too_many_arguments)]
pub fn advance_spatial_driven_joint_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[DrivenJoint],
    integrator: RigidBodyIntegrator,
    loads: &[PreparedSpatialLoads],
    duration_s: f64,
    controls: DistanceGraphControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DrivenJointGraphStep<SpatialStepReceipt>, DistanceGraphError> {
    poll(&mut cancelled)?;
    if states.len() != properties.len() || states.len() != loads.len() {
        return Err(DistanceStepError::Invalid("spatial driven body/load shapes").into());
    }
    limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    let maximum = controls.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut row_count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        row_count += joint.row_count();
        limit("constraints", row_count as u128, maximum as u128)?;
    }
    let extra = loads.len() as u128 * core::mem::size_of::<f64>() as u128;
    let bytes = planned_bytes_for_motion::<JointReactionImpulse, SpatialStepReceipt>(
        states.len(),
        row_count,
    ) + joints.len() as u128 * core::mem::size_of::<DrivenJointStep>() as u128
        + extra;
    let bounds = prepare_bounds(loads, bytes, controls, &mut cancelled)?;
    driven_joints::advance_with(
        states,
        properties,
        joints,
        integrator,
        SpatialFlow {
            loads,
            bounds: &bounds,
        },
        duration_s,
        controls,
        extra,
        &mut cancelled,
    )
}

fn prepare_bounds(
    loads: &[PreparedSpatialLoads],
    bytes: u128,
    controls: DistanceGraphControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, DistanceGraphError> {
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut bounds = reserved(loads.len())?;
    for load in loads {
        poll(cancelled)?;
        bounds.push(load.torque_norm_bound()?);
    }
    Ok(bounds)
}

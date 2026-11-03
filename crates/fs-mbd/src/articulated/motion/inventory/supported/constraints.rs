//! Ideal affine joint locks and gears during supported rigid-payload capture.
//!
//! The loaded inertia and incoming momentum mismatch remain owned by the
//! original capture lifecycle. Native prescribed-base inverse-mass responses
//! form only the bounded row mobility. One simultaneous final impulse solve
//! applies both the mismatch and A^T lambda, retaining the original support
//! reaction, endpoint speed limits and representability checks.
//!
//! Capture's nonnegative relative-motion quadratic Q describes the complete
//! projection. With finite admitted row residuals its identity is
//! delta_K = W_support + (A^T lambda) dot qdot_final - Q.
//! The final-rate pairing remains visible in the original bounded energy
//! defect; midpoint reaction work is reported and capped separately. Neither
//! is subtracted from Q or assigned to heat.

use super::{
    SupportedPayloadCapture, SupportedPayloadControls, capture_evaluated,
};
use crate::articulated::motion::inventory::{PayloadError, PayloadState, finite, poll};
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
    validate_condition_floor,
};
use crate::articulated::{
    ArticulatedModel, PrescribedBaseImpulseResponse, apply_prescribed_base_impulses,
    prescribed_base_impulse_delta,
};
use fs_ga::Wrench;
use fs_math::dd::Dd;

/// Independent-row conditioning and actual ideal-reaction work acceptance.
#[derive(Clone, Copy, Debug)]
pub struct SupportedPayloadConstraintControls {
    /// Native normalized row-mobility reciprocal-condition floor in
    /// [4096*EPSILON, 1]. Dependent or unresolved rows are never regularized.
    pub minimum_reciprocal_condition: f64,
    /// Nonnegative cap on absolute actual midpoint joint-reaction work, J.
    /// This is separate from the capture projection-loss budget.
    pub maximum_reaction_work_j: f64,
}

/// Simultaneous ideal reactions and actual accepted affine/tangent residuals.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedPayloadConstraintReaction {
    /// Row impulses in input order, N m s for Angle and N s for Length.
    pub multipliers: Vec<f64>,
    /// Applied A^T lambda in native compact joint order, N m s or N s.
    /// Empty when there are no rows; no generalized joint impulse is applied.
    pub joint_impulse: Vec<f64>,
    /// Actual A q - b; configuration is unchanged during capture.
    pub position_residuals: Vec<f64>,
    /// Actual A qdot in row units/second before the event.
    pub velocity_residuals_before: Vec<f64>,
    /// Actual A qdot in row units/second at the returned endpoint.
    pub velocity_residuals_after: Vec<f64>,
    /// Applied joint impulse paired with the actual midpoint joint rates, J.
    /// Exact tangent endpoints give zero; finite admitted defects retain work.
    pub midpoint_work_j: f64,
    /// Applied joint impulse paired with the actual final joint rates, J.
    /// This remains in capture.balance.energy_residual_j, under its tolerance;
    /// it is not subtracted from the nonnegative capture-projection quadratic.
    pub final_rate_pairing_j: f64,
    /// Reciprocal condition of the normalized loaded-inertia row mobility.
    /// Equal to one for the empty row set.
    pub reciprocal_condition: f64,
}

impl SupportedPayloadConstraintReaction {
    fn empty() -> Self {
        Self {
            multipliers: Vec::new(),
            joint_impulse: Vec::new(),
            position_residuals: Vec::new(),
            velocity_residuals_before: Vec::new(),
            velocity_residuals_after: Vec::new(),
            midpoint_work_j: 0.0,
            final_rate_pairing_j: 0.0,
            reciprocal_condition: 1.0,
        }
    }
}

/// Original accepted capture plus the ideal joint impulse acting in that event.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstrainedSupportedPayloadCapture {
    /// Original inventory token and actual finite-payload/support accounts.
    /// Its joint momentum residual subtracts reaction.joint_impulse.
    /// With finite row defects, removed_relative_energy_j is the complete
    /// nonnegative projection quadratic, not an inferred thermal deposit.
    pub capture: SupportedPayloadCapture,
    /// Independent reaction work, row residuals and row-mobility diagnostic.
    pub reaction: SupportedPayloadConstraintReaction,
}

/// Capture a finite rigid payload while retaining existing ideal gear/lock rows.
///
/// Initial configuration and rates must already satisfy every constant
/// A q = b and A qdot = 0 within their declared tolerances. The incoming payload
/// changes the complete loaded inertia before the simultaneous row response is
/// solved. Mobile scalar roots and arbitrary native branches remain supported.
/// Numerical inverse-mass probes are not physical endpoint velocities; only the
/// completed native event is tested against speed limits. Ideal impulses do not
/// pass through actuator force/torque effort limits.
///
/// At most 32 independent rows and 4096 links are admitted by the original
/// owners. Only the row-space mobility is dense. Empty rows return exactly the
/// original capture and an empty reaction. No contact discovery, constraint
/// activation, varying gear ratio, backlash or nonlinear loop is inferred.
/// Cancellation or any refusal publishes no partial inventory/state/account.
#[allow(clippy::too_many_arguments)]
pub fn capture_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target_link: usize,
    payload: PayloadState,
    rows: &[JointConstraint],
    controls: SupportedPayloadControls<'_>,
    constraint_controls: SupportedPayloadConstraintControls,
    cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedSupportedPayloadCapture, PayloadError> {
    capture_evaluated(
        model,
        state,
        target_link,
        payload,
        controls,
        Some((rows, constraint_controls)),
        cancelled,
    )
}

pub(super) fn capture_response(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    zero: &[f64],
    body_impulses: &[Wrench],
    constraints: Option<(&[JointConstraint], SupportedPayloadConstraintControls)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(PrescribedBaseImpulseResponse, SupportedPayloadConstraintReaction), PayloadError> {
    let Some((rows, controls)) = constraints else {
        return Ok((
            apply_prescribed_base_impulses(
                model, state.base, &state.joint_position, &state.joint_velocity,
                zero, body_impulses,
            )?,
            SupportedPayloadConstraintReaction::empty(),
        ));
    };
    validate_condition_floor(controls.minimum_reciprocal_condition)?;
    if !controls.maximum_reaction_work_j.is_finite() || controls.maximum_reaction_work_j < 0.0 {
        return Err(ConstraintError::Invalid("supported capture reaction-work budget").into());
    }
    let initial = admit_constraints(
        rows, &state.joint_position, &state.joint_velocity, cancelled,
    )?;
    if rows.is_empty() {
        return Ok((
            apply_prescribed_base_impulses(
                model, state.base, &state.joint_position, &state.joint_velocity,
                zero, body_impulses,
            )?,
            SupportedPayloadConstraintReaction::empty(),
        ));
    }
    let mut prepared = PreparedJointConstraints::prescribed_trial(
        model, &state.joint_position, rows,
        controls.minimum_reciprocal_condition, cancelled,
    )?;
    poll(cancelled)?;
    let free = prescribed_base_impulse_delta(
        model, &state.joint_position, zero, body_impulses,
    )?;
    poll(cancelled)?;
    let mut candidate = Vec::with_capacity(model.dof_count());
    for (index, (&old, &delta)) in state.joint_velocity.iter().zip(&free.velocity.joint).enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        candidate.push(finite(old + delta, "supported capture free joint rate")?);
    }
    let rhs = apply_rows(rows, &candidate, cancelled)?
        .into_iter().map(|value| -value).collect::<Vec<_>>();
    let reaction = prepared.solve_response(&rhs, cancelled)?;
    poll(cancelled)?;
    // Reapply the simultaneous physical inputs through the original event
    // owner. Its native reaction and actual endpoint guards remain independent
    // of cached columns, rather than accepting a projected/clipped velocity.
    let response = apply_prescribed_base_impulses(
        model, state.base, &state.joint_position, &state.joint_velocity,
        &reaction.joint_effort, body_impulses,
    )?;
    poll(cancelled)?;
    let final_rows = admit_constraints(
        rows, &state.joint_position, &response.joint_velocity_after, cancelled,
    )?;
    let mut midpoint = Dd::ZERO;
    let mut final_pairing = Dd::ZERO;
    for (index, ((&applied, &old), &new)) in reaction.joint_effort.iter()
        .zip(&state.joint_velocity).zip(&response.joint_velocity_after).enumerate()
    {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let applied = Dd::from_f64(applied);
        midpoint = midpoint + applied
            * (Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new));
        final_pairing = final_pairing + applied * Dd::from_f64(new);
    }
    let midpoint_work_j = finite(midpoint.to_f64(), "supported capture reaction work")?;
    let final_rate_pairing_j = finite(
        final_pairing.to_f64(), "supported capture final reaction pairing",
    )?;
    if midpoint_work_j.abs() > controls.maximum_reaction_work_j {
        return Err(PayloadError::ConstraintWork {
            work_j: midpoint_work_j,
            maximum_j: controls.maximum_reaction_work_j,
        });
    }
    poll(cancelled)?;
    Ok((response, SupportedPayloadConstraintReaction {
        multipliers: reaction.multipliers,
        joint_impulse: reaction.joint_effort,
        position_residuals: initial.position,
        velocity_residuals_before: initial.velocity,
        velocity_residuals_after: final_rows.velocity,
        midpoint_work_j,
        final_rate_pairing_j,
        reciprocal_condition: prepared.reciprocal_condition,
    }))
}

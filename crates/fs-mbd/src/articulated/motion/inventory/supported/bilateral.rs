//! Supported payload inventory retaining material closures and ideal joint rows.
//!
//! The original capture/release lifecycle owns inventory and actual motion.
//! One native prescribed-base event transmits the payload mismatch or recoil
//! through distance, spherical and gear/lock rows simultaneously. Its support
//! reaction remains the only external momentum exchange. Full body and joint
//! reactions are projected and paired in the inventory's common energy frame.

use super::constraints::{SupportedPayloadConstraintControls, admit_controls};
use super::release::{
    SupportedPayloadRelease, SupportedPayloadReleaseControls, SupportedPayloadSeparationControls,
    release_evaluated,
};
use super::{SupportedPayloadCapture, SupportedPayloadControls, capture_evaluated};
use crate::articulated::motion::inventory::bilateral::{
    constraint_work, project_constraint_impulses,
};
use crate::articulated::motion::inventory::release::SeparationImpulse;
use crate::articulated::motion::inventory::{
    PayloadAttachment, PayloadError, PayloadState, finite, poll, twist_in_frame,
};
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::motion::reconfigure::constraints::{ConstraintError, JointConstraint};
use crate::articulated::motion::spatial_loads::implicit::distance::prescribed::{
    PrescribedBilateralImpulseResponse, apply_impulses,
};
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop,
};
use crate::articulated::motion::sum_products;
use crate::articulated::{ArticulatedModel, Kinematics};
use fs_ga::{Vec3, Wrench};
use fs_math::dd::Dd;

/// Actual mixed constraint reaction and its complete inventory-frame account.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedBilateralPayloadReaction {
    /// Native event on the loaded capture inertia or restored release inertia.
    /// Its source/reaction work is in the lab frame. For capture, its initial
    /// state carries the payload rigidly and is an auxiliary mismatch event,
    /// not the physical tree-plus-free-payload inventory before capture.
    pub impulse: PrescribedBilateralImpulseResponse,
    /// Zero-rate RNEA projection of actual body plus joint constraint impulses.
    /// Ordered as base angular/linear six-vector, then compact joint impulses.
    /// Base entries use the actual root joint configuration, including mobile
    /// roots. These entries are physical row covectors, not support impulses.
    pub generalized_constraint_impulse: Vec<f64>,
    /// Same constant translating energy frame as the inventory balance, m/s.
    pub accounting_boost_world: Vec3,
    /// Complete body-plus-joint impulse paired with actual mean rates, J.
    /// Release subtracts this work separately from support and source work.
    pub midpoint_work_j: f64,
    /// Complete constraint impulse paired with actual final rates, J.
    /// Capture retains this finite tangency defect in its bounded energy
    /// residual; it never subtracts midpoint work from capture loss.
    pub final_rate_pairing_j: f64,
}

/// Accepted supported capture with the material and joint rows it preserved.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedBilateralPayloadCapture {
    /// Original inventory identity, actual post-event state and physical audits.
    pub capture: SupportedPayloadCapture,
    /// Simultaneous mixed reaction, condition and finite constraint work.
    pub reaction: SupportedBilateralPayloadReaction,
}

/// Accepted supported release with its retained-tree mixed constraint event.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedBilateralPayloadRelease {
    /// Restored inertia, actual tree/payload recoil, support and source work.
    pub release: SupportedPayloadRelease,
    /// Zero for passive release; actual mixed reaction for driven separation.
    pub reaction: SupportedBilateralPayloadReaction,
}

/// Capture while preserving existing distance, spherical and affine joint rows.
///
/// All rows must already meet their position and velocity tolerances. The
/// original lifecycle combines the full rotated payload inertia, then the
/// native prescribed-base owner applies mismatch and all reactions together.
/// Fixed and driven supports and mobile scalar root joints are supported.
/// Only the final combined physical endpoint faces native speed limits.
///
/// The nonnegative capture quadratic Q obeys
/// delta_K - W_support + Q = J_constraint dot v_final in the common frame,
/// up to the independently bounded numerical defect. Actual midpoint reaction
/// work is capped separately. No constraint engagement or heat is inferred.
/// Empty rows still use the same native event and common-frame audit.
/// Cancellation or refusal publishes no partial inventory, state or account.
#[allow(clippy::too_many_arguments)]
pub fn capture_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target_link: usize,
    payload: PayloadState,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: SupportedPayloadControls<'_>,
    constraint_controls: SupportedPayloadConstraintControls,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<SupportedBilateralPayloadCapture, PayloadError> {
    let (result, reaction) = capture_evaluated(
        model,
        state,
        target_link,
        payload,
        controls,
        Some((rows, constraint_controls)),
        Some((loops, config)),
        cancelled,
    )?;
    Ok(SupportedBilateralPayloadCapture {
        capture: result.capture,
        reaction: reaction.ok_or(PayloadError::Invalid(
            "missing supported mixed capture reaction",
        ))?,
    })
}

/// Detach passively while admitting material and joint rows on retained inertia.
///
/// The original state is copied exactly. Geometry, tangency and full mixed-row
/// mobility are admitted, but a finite accepted rate defect is never repaired
/// by an invented impulse. Reaction and work remain zero.
#[allow(clippy::too_many_arguments)]
pub fn release_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: SupportedPayloadReleaseControls<'_>,
    constraint_controls: SupportedPayloadConstraintControls,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<SupportedBilateralPayloadRelease, PayloadError> {
    let (result, reaction) = release_evaluated(
        model,
        state,
        attachment,
        None,
        controls,
        0.0,
        Some((rows, constraint_controls)),
        Some((loops, config)),
        cancelled,
    )?;
    Ok(SupportedBilateralPayloadRelease {
        release: result.release,
        reaction: reaction.ok_or(PayloadError::Invalid(
            "missing supported mixed release reaction",
        ))?,
    })
}

/// Separate by one internal point impulse and free couple, retaining all rows.
///
/// The restored tree and full payload tensor respond to equal/opposite actions.
/// Material and gear reactions share one retained-inertia solve with recoil;
/// the prescribed support retains its independent six-component reaction and
/// signed work. Actual kinetic change subtracts support work, source work and
/// full common-frame midpoint constraint work separately. Numerical inventory
/// and work defects remain visible under the original acceptance tolerances.
///
/// The supported source policy is unchanged: negative separation work refuses,
/// and nonnegative supplied work must fit its explicit budget. Zero impulse
/// follows the exact passive state branch with retained-row admission.
#[allow(clippy::too_many_arguments)]
pub fn release_payload_with_impulse(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    impulse: SeparationImpulse,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: SupportedPayloadSeparationControls<'_>,
    constraint_controls: SupportedPayloadConstraintControls,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<SupportedBilateralPayloadRelease, PayloadError> {
    let (result, reaction) = release_evaluated(
        model,
        state,
        attachment,
        Some(impulse),
        controls.balance,
        controls.maximum_separation_work_j,
        Some((rows, constraint_controls)),
        Some((loops, config)),
        cancelled,
    )?;
    Ok(SupportedBilateralPayloadRelease {
        release: result.release,
        reaction: reaction.ok_or(PayloadError::Invalid(
            "missing supported mixed release reaction",
        ))?,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn native_response(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    body_impulses: &[Wrench],
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: SupportedPayloadConstraintControls,
    config: BilateralImpulseConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PrescribedBilateralImpulseResponse, PayloadError> {
    admit_controls(controls)?;
    let zero = vec![0.0; model.dof_count()];
    let response = apply_impulses(
        model,
        state,
        &zero,
        body_impulses,
        loops,
        rows,
        config,
        cancelled,
    )?;
    if !response.reciprocal_condition.is_finite() {
        return Err(PayloadError::Invalid(
            "supported mixed reciprocal condition",
        ));
    }
    if response.reciprocal_condition < controls.minimum_reciprocal_condition {
        return Err(ConstraintError::IllConditioned {
            reciprocal_condition: response.reciprocal_condition,
            minimum: controls.minimum_reciprocal_condition,
        }
        .into());
    }
    Ok(response)
}

/// Complete the receipt from the actual inventory endpoints, retaining Dd work
/// until release subtracts it from its support/source/kinetic balance.
#[allow(clippy::too_many_arguments)]
pub(super) fn finish_reaction(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    before: &Kinematics,
    after: &Kinematics,
    joint_velocity_after: &[f64],
    boost: Vec3,
    impulse: PrescribedBilateralImpulseResponse,
    controls: SupportedPayloadConstraintControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(SupportedBilateralPayloadReaction, Dd), PayloadError> {
    let generalized_constraint_impulse = project_constraint_impulses(
        model,
        state.base.world_from_base,
        &state.joint_position,
        &impulse.constraint_body_impulses,
        &impulse.constraint_joint_impulses,
        cancelled,
    )?;
    let work = constraint_work(
        &impulse.constraint_body_impulses,
        &impulse.constraint_joint_impulses,
        before,
        after,
        &state.joint_velocity,
        joint_velocity_after,
        boost,
        cancelled,
    )?;
    let midpoint_work_j = finite(work.to_f64(), "supported mixed constraint work")?;
    if midpoint_work_j.abs() > controls.maximum_reaction_work_j {
        return Err(PayloadError::ConstraintWork {
            work_j: midpoint_work_j,
            maximum_j: controls.maximum_reaction_work_j,
        });
    }
    let mut final_pairing = sum_products(
        impulse.constraint_joint_impulses.iter().copied(),
        joint_velocity_after.iter().copied(),
    );
    for (index, body) in impulse.constraint_body_impulses.iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        final_pairing = final_pairing
            + sum_products(
                body.to_array(),
                twist_in_frame(after.body_twist[index], after.world_from_link[index], boost)?
                    .to_array(),
            );
    }
    let final_rate_pairing_j = finite(final_pairing.to_f64(), "supported mixed final pairing")?;
    poll(cancelled)?;
    Ok((
        SupportedBilateralPayloadReaction {
            impulse,
            generalized_constraint_impulse,
            accounting_boost_world: boost,
            midpoint_work_j,
            final_rate_pairing_j,
        },
        work,
    ))
}

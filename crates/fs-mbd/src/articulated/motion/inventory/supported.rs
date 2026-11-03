//! Instantaneous rigid payload capture with a fixed or driven support.
//!
//! The original inventory owner combines the full rotated payload inertia.
//! The original prescribed-base impulse owner absorbs the incoming payload's
//! momentum mismatch through that loaded inertia, including mobile root joints.
//! Only the support exchanges external momentum; a driven support also supplies
//! signed work. Capture loss is the nonnegative relative-velocity quadratic,
//! never a raw decrease of lab-frame kinetic energy.
//!
//! This is a perfectly inelastic rigid attachment at the actual supplied pose.
//! It does not discover contact, impose extra gear/loop rows, resolve a finite
//! grasp duration, or assign the removed mechanical energy to a thermal body.

/// Passive or internally driven release of an exact captured inventory.
#[path = "supported/release.rs"]
pub mod release;

use super::super::prescribed::PrescribedBaseMotionState;
use super::super::reconfigure::{self, ReconfigurationControls, ReconfigurationError};
use super::super::{
    MAX_FREE_FLOATING_STEP_LINKS, SpatialSum, kinetic_link_increment, sum_products,
};
use super::{
    PayloadAttachment, PayloadError, PayloadState, account_boost, add_wrench, combine_inertias,
    energy_in_frame, finite, payload_at_link, poll, replace_link_inertia, twist_in_frame,
};
use crate::articulated::{
    ArticulatedModel, ArticulatedMomentum, BaseState, apply_prescribed_base_impulses,
    forward_kinematics, inverse_dynamics, momentum_from_kinematics, validate_twist,
    validate_wrench,
};
use fs_ga::{Vec3, Wrench};
use fs_math::dd::Dd;

/// Absolute acceptance controls for the actual returned inventory event.
#[derive(Clone, Copy, Debug)]
pub struct SupportedPayloadControls<'a> {
    /// Positive per-component world angular impulse residual allowance, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Positive per-component world linear impulse residual allowance, N s.
    pub linear_momentum_tolerance_ns: f64,
    /// Positive allowances for each native joint momentum residual, N m s or N s.
    pub joint_momentum_tolerances: &'a [f64],
    /// Nonnegative absolute common-frame kinetic/work residual allowance, J.
    pub energy_tolerance_j: f64,
    /// Nonnegative maximum perfectly inelastic relative-motion energy loss, J.
    pub maximum_removed_relative_energy_j: f64,
}

/// Complete-system momentum and work accounting across a supported inventory event.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedPayloadBalance {
    /// Complete system before the event: tree and any incoming free payload.
    pub momentum_before: ArticulatedMomentum,
    /// Complete system after the event: tree and any released free payload.
    pub momentum_after: ArticulatedMomentum,
    /// Actual linear momentum change, before rounding whole-system totals, N s.
    pub linear_momentum_change_world: Vec3,
    /// Actual spin plus orbital increment about the world origin, N m s.
    pub angular_momentum_change_world: Vec3,
    /// Independent native ABA support impulse ON the complete system.
    /// Passive release has no impulse; driven separation retains the reaction.
    /// Expressed about the base origin in base axes, in [N m s, N s].
    pub support_impulse_body: Wrench,
    /// The same reaction about the world origin, in world axes.
    pub support_impulse_world: Wrench,
    /// Actual linear increment minus the independent support impulse, N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual world-origin angular increment minus support impulse, N m s.
    pub angular_momentum_residual_world: Vec3,
    /// Native RNEA projection of actual body increments or release decomposition.
    /// Capture and equal/opposite separation have no external relative-joint impulse.
    pub joint_momentum_residual: Vec<f64>,
    /// Largest absolute residual divided by its declared allowance, at most one.
    pub maximum_normalized_momentum_residual: f64,
    /// Actual lab-frame kinetic change or passive-release decomposition defect, J.
    pub kinetic_energy_change_j: f64,
    /// Support impulse paired with the unchanged prescribed base twist, J.
    pub support_work_j: f64,
    /// Constant world velocity of the inertial energy-accounting frame, m/s.
    pub accounting_boost_world: Vec3,
    /// Complete-system kinetic energy in that frame before the event, J.
    pub kinetic_before_j: f64,
    /// Complete-system kinetic energy in the same frame after the event, J.
    pub kinetic_after_j: f64,
    /// Actual common-frame kinetic change or passive-release decomposition defect, J.
    pub accounting_kinetic_energy_change_j: f64,
    /// Support work paired with the base twist in that same frame, J.
    pub accounting_support_work_j: f64,
    /// Common-frame kinetic change minus support/separation work plus capture loss, J.
    /// This finite numerical defect is never assigned to removed energy.
    pub energy_residual_j: f64,
}

/// Accepted inventory and supported motion state, ready for continuation.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedPayloadCapture {
    /// Unchanged tree topology and joints with the target's loaded full inertia.
    pub model_after: ArticulatedModel,
    /// Unchanged time, support boundary and coordinates; simultaneous new rates.
    pub state_after: PrescribedBaseMotionState,
    /// Original inventory identity, including the retained target inertia.
    pub attachment: PayloadAttachment,
    /// Actual attached payload pose and velocity in the returned state.
    pub payload_after: PayloadState,
    /// Sum of nonnegative old-body and payload relative-motion quadratics, J.
    pub removed_relative_energy_j: f64,
    /// Independent impulse, generalized momentum and support-work audits.
    pub balance: SupportedPayloadBalance,
}

/// Capture a free rigid payload on any link of a prescribed-base tree.
///
/// At the unchanged configuration, the loaded tree initially assigns the
/// payload the target link's velocity. Its difference from the real incoming
/// velocity is an impulse on the loaded inertia. Native prescribed-base ABA
/// maps that mismatch to the joint-rate jump and the independent support
/// reaction. This is the same impulse map used for other supported events;
/// no large surrogate support mass, dense generalized matrix or rate clipping
/// is introduced. Fixed, revolute, prismatic and helical root joints retain
/// their original geometry and endpoint speed limits.
///
/// The support pose and twist cannot jump during this event. Its declared
/// smooth acceleration is retained but contributes no finite impulse. Actual
/// full-system momentum increments include the payload's before/after motion.
/// Native zero-bias RNEA independently checks relative-joint momentum. Energy
/// is audited in one translating inertial frame, retaining support work and
/// the separate nonnegative capture loss.
///
/// All inputs remain borrowed. Cancellation is checked between bounded native
/// calls and every 64 links in additional passes; at most 4096 links are
/// admitted. Any refusal publishes no model, state, attachment or spent budget.
pub fn capture_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target_link: usize,
    payload: PayloadState,
    controls: SupportedPayloadControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SupportedPayloadCapture, PayloadError> {
    poll(&mut cancelled)?;
    admit(model, state, target_link, controls)?;
    let before = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    validate_twist(payload.twist_body, "supported_payload.incoming_twist")?;
    let link_from_payload = before.world_from_link[target_link]
        .inverse()?
        .compose(payload.world_from_body)?;
    let loaded_inertia = combine_inertias(
        model.links()[target_link].inertia(),
        payload.inertia,
        link_from_payload,
    )?;
    let model_after = replace_link_inertia(model, target_link, loaded_inertia)?;
    let attached_before =
        payload_at_link(&before, target_link, link_from_payload, payload.inertia)?;
    let mismatch = payload
        .inertia
        .momentum(payload.twist_body.minus(attached_before.twist_body))?
        .transform_by(&link_from_payload);
    validate_wrench(mismatch, "supported_payload.momentum_mismatch")?;
    let mut body_impulses = vec![Wrench::default(); model.link_count()];
    body_impulses[target_link] = mismatch;
    let zero = vec![0.0; model.dof_count()];
    poll(&mut cancelled)?;
    let response = apply_prescribed_base_impulses(
        &model_after,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &zero,
        &body_impulses,
    )?;
    poll(&mut cancelled)?;
    let state_after = PrescribedBaseMotionState {
        time_s: state.time_s,
        base: state.base,
        joint_position: state.joint_position.clone(),
        joint_velocity: response.joint_velocity_after,
    };
    let after = forward_kinematics(
        &model_after,
        state_after.base,
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    let payload_after = payload_at_link(&after, target_link, link_from_payload, payload.inertia)?;
    let boost = account_boost(model, &before)?;
    let mut momentum_before = momentum_from_kinematics(model, &before)?;
    let incoming = payload
        .inertia
        .momentum(payload.twist_body)?
        .transform_by(&payload.world_from_body);
    momentum_before.linear_momentum_world = momentum_before.linear_momentum_world + incoming.force;
    momentum_before.angular_momentum_world =
        momentum_before.angular_momentum_world + incoming.torque;
    momentum_before.kinetic_energy = finite(
        momentum_before.kinetic_energy + payload.inertia.kinetic_energy(payload.twist_body)?,
        "supported_payload.kinetic_before",
    )?;
    validate_wrench(
        Wrench::new(
            momentum_before.angular_momentum_world,
            momentum_before.linear_momentum_world,
        ),
        "supported_payload.momentum_before",
    )?;
    poll(&mut cancelled)?;

    let mut change = SpatialSum::default();
    let mut kinetic_change = Dd::ZERO;
    let mut centered_change = Dd::ZERO;
    let mut removed = Dd::ZERO;
    // Negated physical increments are external inputs to zero-bias RNEA,
    // which subtracts them and therefore projects the positive increments.
    let mut projected = Vec::with_capacity(model.link_count());
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let pose = before.world_from_link[index];
        let old = before.body_twist[index];
        let new = after.body_twist[index];
        let delta = new.minus(old);
        let impulse = link.inertia().momentum(delta)?;
        change.add_world_wrench(
            pose.rotation().rotate(impulse.force)?,
            pose.rotation().rotate(impulse.torque)?,
            &[pose.translation()],
        );
        projected.push(negate(impulse));
        kinetic_change =
            kinetic_change + kinetic_link_increment(link.inertia(), pose, old, pose, new)?;
        centered_change = centered_change
            + kinetic_link_increment(
                link.inertia(),
                pose,
                twist_in_frame(old, pose, boost)?,
                pose,
                twist_in_frame(new, pose, boost)?,
            )?;
        removed = removed + Dd::from_f64(link.inertia().kinetic_energy(delta)?);
    }
    let payload_delta = payload_after.twist_body.minus(payload.twist_body);
    let payload_impulse = payload.inertia.momentum(payload_delta)?;
    change.add_world_wrench(
        payload
            .world_from_body
            .rotation()
            .rotate(payload_impulse.force)?,
        payload
            .world_from_body
            .rotation()
            .rotate(payload_impulse.torque)?,
        &[payload.world_from_body.translation()],
    );
    projected[target_link] = add_wrench(
        projected[target_link],
        negate(payload_impulse.transform_by(&link_from_payload)),
    );
    kinetic_change = kinetic_change
        + kinetic_link_increment(
            payload.inertia,
            payload.world_from_body,
            payload.twist_body,
            payload.world_from_body,
            payload_after.twist_body,
        )?;
    centered_change = centered_change
        + kinetic_link_increment(
            payload.inertia,
            payload.world_from_body,
            twist_in_frame(payload.twist_body, payload.world_from_body, boost)?,
            payload.world_from_body,
            twist_in_frame(payload_after.twist_body, payload.world_from_body, boost)?,
        )?;
    removed = removed + Dd::from_f64(payload.inertia.kinetic_energy(payload_delta)?);
    poll(&mut cancelled)?;
    let joint_momentum_residual = inverse_dynamics(
        model,
        BaseState::stationary(state.base.world_from_base),
        &state.joint_position,
        &zero,
        &zero,
        Vec3::default(),
        &projected,
    )?
    .generalized_force;
    poll(&mut cancelled)?;
    let mut support = SpatialSum::default();
    support.add_world_wrench(
        state
            .base
            .world_from_base
            .rotation()
            .rotate(response.support_impulse_body.force)?,
        state
            .base
            .world_from_base
            .rotation()
            .rotate(response.support_impulse_body.torque)?,
        &[state.base.world_from_base.translation()],
    );
    let residual = change
        .subtract(support)
        .checked_wrench("supported_payload.momentum_residual")?;
    let actual = change.checked_wrench("supported_payload.momentum_change")?;
    let maximum_normalized_momentum_residual =
        audit_momentum(residual, &joint_momentum_residual, controls)?;
    let support_work = sum_products(
        response.support_impulse_body.to_array(),
        twist_in_frame(state.base.twist_body, state.base.world_from_base, boost)?.to_array(),
    );
    let removed_relative_energy_j = finite(removed.to_f64(), "supported_payload.capture_loss")?;
    if removed_relative_energy_j < 0.0 {
        return Err(PayloadError::Invalid("negative supported capture loss"));
    }
    if removed_relative_energy_j > controls.maximum_removed_relative_energy_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: removed_relative_energy_j,
            maximum_j: controls.maximum_removed_relative_energy_j,
        });
    }
    let energy_residual_j = finite(
        (centered_change - support_work + removed).to_f64(),
        "supported_payload.energy_residual",
    )?;
    if energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(PayloadError::EnergyResidual {
            residual_j: energy_residual_j,
            maximum_j: controls.energy_tolerance_j,
        });
    }
    let kinetic_before_j = energy_in_frame(model, &before, Some(&payload), boost, &mut cancelled)?;
    let kinetic_after_j = energy_in_frame(&model_after, &after, None, boost, &mut cancelled)?;
    let balance = SupportedPayloadBalance {
        momentum_before,
        momentum_after: response.momentum_after,
        linear_momentum_change_world: actual.force,
        angular_momentum_change_world: actual.torque,
        support_impulse_body: response.support_impulse_body,
        support_impulse_world: response.support_impulse_world,
        linear_momentum_residual_world: residual.force,
        angular_momentum_residual_world: residual.torque,
        joint_momentum_residual,
        maximum_normalized_momentum_residual,
        kinetic_energy_change_j: finite(
            kinetic_change.to_f64(),
            "supported_payload.kinetic_change",
        )?,
        support_work_j: response.support_work,
        accounting_boost_world: boost,
        kinetic_before_j,
        kinetic_after_j,
        accounting_kinetic_energy_change_j: finite(
            centered_change.to_f64(),
            "supported_payload.centered_change",
        )?,
        accounting_support_work_j: finite(
            support_work.to_f64(),
            "supported_payload.centered_support_work",
        )?,
        energy_residual_j,
    };
    poll(&mut cancelled)?;
    Ok(SupportedPayloadCapture {
        model_after,
        state_after,
        attachment: PayloadAttachment {
            target_link,
            retained_link: model.links()[target_link].clone(),
            loaded_inertia,
            link_from_payload,
            payload_inertia: payload.inertia,
        },
        payload_after,
        removed_relative_energy_j,
        balance,
    })
}

fn admit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target: usize,
    controls: SupportedPayloadControls<'_>,
) -> Result<(), PayloadError> {
    if target >= model.link_count() || !state.time_s.is_finite() {
        return Err(PayloadError::Invalid("supported payload target or clock"));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(ReconfigurationError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    reconfigure::validate_controls(
        model,
        ReconfigurationControls {
            relative_momentum_tolerance: 0.0,
            angular_momentum_tolerance_nms: controls.angular_momentum_tolerance_nms,
            linear_momentum_tolerance_ns: controls.linear_momentum_tolerance_ns,
            joint_momentum_tolerances: controls.joint_momentum_tolerances,
            maximum_parameter_work_j: controls.maximum_removed_relative_energy_j,
        },
    )?;
    if !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j < 0.0 {
        return Err(PayloadError::Invalid("supported payload energy tolerance"));
    }
    Ok(())
}

fn negate(value: Wrench) -> Wrench {
    Wrench::new(value.torque.scale(-1.0), value.force.scale(-1.0))
}

fn audit_momentum(
    physical: Wrench,
    joint: &[f64],
    controls: SupportedPayloadControls<'_>,
) -> Result<f64, PayloadError> {
    if joint.len() != controls.joint_momentum_tolerances.len() {
        return Err(PayloadError::Invalid(
            "supported payload joint residual dimensions",
        ));
    }
    let mut maximum: f64 = 0.0;
    for (index, value) in physical.to_array().into_iter().enumerate() {
        let allowance = if index < 3 {
            controls.angular_momentum_tolerance_nms
        } else {
            controls.linear_momentum_tolerance_ns
        };
        maximum = maximum.max(finite(
            (value / allowance).abs(),
            "supported payload normalized residual",
        )?);
    }
    for (&value, &allowance) in joint.iter().zip(controls.joint_momentum_tolerances) {
        maximum = maximum.max(finite(
            (value / allowance).abs(),
            "supported payload joint normalized residual",
        )?);
    }
    if maximum > 1.0 {
        return Err(ReconfigurationError::MomentumResidual {
            maximum_normalized: maximum,
        }
        .into());
    }
    Ok(maximum)
}

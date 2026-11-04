//! Rigid payload release with fixed or driven articulated support.
//!
//! The original inventory restores the retained tensor and current rigid
//! payload velocity. Optional point impulses and free couples act equally and
//! oppositely on both surviving systems through their native inverse inertias.
//! Support impulse and signed work remain external. Release work comes from the
//! supplied internal actuator; numerical decomposition defects are not recoil,
//! mechanical loss, or an additional source of energy.

use super::super::release::{
    SeparationImpulse, admit_attachment, admit_separation_impulse, separation_impulses,
};
use super::super::{
    PayloadAttachment, PayloadError, PayloadState, account_boost, energy_in_frame, finite,
    payload_at_link, poll, replace_link_inertia, twist_in_frame,
};
use super::bilateral::{self, SupportedBilateralPayloadReaction};
use super::constraints::{
    ConstrainedSupportedPayloadRelease, SupportedPayloadConstraintControls,
    SupportedPayloadConstraintReaction, impulse_response, passive_reaction,
};
use super::{SupportedPayloadBalance, SupportedPayloadControls, admit, audit_momentum, negate};
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::motion::reconfigure::constraints::JointConstraint;
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop,
};
use crate::articulated::motion::{SpatialSum, kinetic_link_increment, sum_products};
use crate::articulated::{
    ArticulatedModel, BaseState, forward_kinematics, inverse_dynamics, momentum_from_kinematics,
    solve_floating_base_system, validate_twist, validate_wrench,
};
use fs_ga::{Vec3, Wrench};
use fs_math::dd::Dd;

/// Independent absolute allowances for the complete inventory event.
#[derive(Clone, Copy, Debug)]
pub struct SupportedPayloadReleaseControls<'a> {
    /// Positive world angular momentum residual allowance per axis, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Positive world linear momentum residual allowance per axis, N s.
    pub linear_momentum_tolerance_ns: f64,
    /// Positive allowances for native joint momentum residuals, N m s or N s.
    pub joint_momentum_tolerances: &'a [f64],
    /// Nonnegative common-frame kinetic/work residual allowance, J.
    pub energy_tolerance_j: f64,
}

/// Balance admission and explicit internal energy source for driven separation.
#[derive(Clone, Copy, Debug)]
pub struct SupportedPayloadSeparationControls<'a> {
    /// Independent world/joint momentum and common-frame work allowances.
    pub balance: SupportedPayloadReleaseControls<'a>,
    /// Nonnegative maximum mechanical work supplied by the impulse pair, J.
    /// The composing actuator must debit this actual work from its own source.
    pub maximum_separation_work_j: f64,
}

/// Restored supported inventory and the actual newly independent rigid body.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedPayloadRelease {
    /// Same topology and joints; exact retained target inertia from capture.
    pub model_after: ArticulatedModel,
    /// Original time, base boundary and coordinates, with native recoiled rates.
    /// Passive or zero-impulse release preserves every state bit.
    pub state_after: PrescribedBaseMotionState,
    /// Original payload tensor at the current pose with its actual outgoing twist.
    pub payload_after: PayloadState,
    /// Actual nonnegative work supplied by the simultaneous internal impulse pair.
    /// This excludes support work and is exactly zero for passive release.
    pub supplied_work_j: f64,
    /// Original impulse descriptor, including its shared application point.
    pub applied_impulse: Option<SeparationImpulse>,
    /// Applied impulse about the payload origin, expressed in payload axes.
    pub impulse_on_payload_body: Wrench,
    /// Equal/opposite applied impulse about the target link origin in link axes.
    pub impulse_on_retained_link: Wrench,
    /// Actual complete-system momentum and common-frame work diagnostics.
    /// Energy residual subtracts support work and supplied_work_j, plus any
    /// full ideal-reaction midpoint work from a constrained wrapper.
    pub balance: SupportedPayloadBalance,
}

/// Detach a captured payload without applying a separation impulse.
///
/// The original inventory validator requires the target's exact loaded inertia,
/// name, parent, joint and reference transform. The current coordinates and
/// rates may have evolved since capture; other links retain their current
/// inertias. Nested captures on one link therefore detach in reverse order,
/// and a repeated release against the unloaded model refuses.
///
/// The retained state is copied exactly. The released body inherits the native
/// current rigid velocity field, including angular velocity and offset-origin
/// motion. Its inertia is never reconstructed by subtracting loaded tensors.
/// Physical impulse, support work and energy loss are zero. Independent
/// decomposition defects remain visible under the declared tolerances.
///
/// This only removes the ideal rigid attachment. Collision discovery, ongoing
/// contact, fracture energy and finite-duration separation remain separate
/// models. Refusal or cancellation publishes no partial state.
pub fn release_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    controls: SupportedPayloadReleaseControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<SupportedPayloadRelease, PayloadError> {
    release_evaluated(
        model, state, attachment, None, controls, 0.0, None, None, cancelled,
    )
    .map(|(result, _)| result.release)
}

/// Detach and drive separation by one internal point impulse and free couple.
///
/// Uses the original free-release SeparationImpulse convention: all vectors
/// are in the current target link axes and both actions share one physical
/// point. The retained tree receives the opposite action through native
/// prescribed-base ABA; the full payload tensor uses the original six-component
/// inverse inertia. All retained joint rates respond simultaneously. No dense
/// generalized mass matrix, surrogate support mass or force-duration conversion
/// is introduced.
///
/// The support pose/twist cannot jump. Its independent native reaction may
/// exchange both momentum and signed work. The internal pair begins at one
/// rigid velocity field, so its initial work cancels; actual supplied work is
/// half each impulse paired with its resulting velocity jump. This positive
/// quadratic is evaluated without subtracting large lab-frame energies.
/// Whole-system and joint momentum, support work, and actual kinetic change
/// are audited independently, including original inventory decomposition error.
///
/// Original endpoint speed limits apply; impulses are not force/torque effort.
/// A completely erased outgoing payload jump refuses. Zero impulse preserves
/// the exact passive state and accounts. No contact law, extra gear/loop row,
/// actuator energy reservoir or thermal destination is inferred.
pub fn release_payload_with_impulse(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    impulse: SeparationImpulse,
    controls: SupportedPayloadSeparationControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<SupportedPayloadRelease, PayloadError> {
    release_evaluated(
        model,
        state,
        attachment,
        Some(impulse),
        controls.balance,
        controls.maximum_separation_work_j,
        None,
        None,
        cancelled,
    )
    .map(|(result, _)| result.release)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn release_evaluated(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    impulse: Option<SeparationImpulse>,
    controls: SupportedPayloadReleaseControls<'_>,
    maximum_work_j: f64,
    constraints: Option<(&[JointConstraint], SupportedPayloadConstraintControls)>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (
        ConstrainedSupportedPayloadRelease,
        Option<SupportedBilateralPayloadReaction>,
    ),
    PayloadError,
> {
    poll(&mut cancelled)?;
    admit_attachment(model, attachment)?;
    if !maximum_work_j.is_finite() || maximum_work_j < 0.0 {
        return Err(PayloadError::Invalid(
            "supported payload separation work budget",
        ));
    }
    let shared_controls = SupportedPayloadControls {
        angular_momentum_tolerance_nms: controls.angular_momentum_tolerance_nms,
        linear_momentum_tolerance_ns: controls.linear_momentum_tolerance_ns,
        joint_momentum_tolerances: controls.joint_momentum_tolerances,
        energy_tolerance_j: controls.energy_tolerance_j,
        maximum_removed_relative_energy_j: 0.0,
    };
    let target = attachment.target_link;
    admit(model, state, target, shared_controls)?;
    let prescribed = impulse.unwrap_or_default();
    admit_separation_impulse(prescribed)?;
    let (impulse_on_retained_link, impulse_on_payload_body) =
        separation_impulses(prescribed, attachment.link_from_payload)?;
    let no_impulse = prescribed.impulse_on_payload_link_ns == Vec3::default()
        && prescribed.couple_on_payload_link_nms == Vec3::default();
    let kin = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    let payload_before = payload_at_link(
        &kin,
        target,
        attachment.link_from_payload,
        attachment.payload_inertia,
    )?;
    let model_after = replace_link_inertia(model, target, attachment.retained_link.inertia())?;
    let momentum_before = momentum_from_kinematics(model, &kin)?;
    poll(&mut cancelled)?;
    let zero = vec![0.0; model.dof_count()];
    let mut state_after = state.clone();
    let mut payload_after = payload_before;
    let mut bilateral_event = None;
    let (after_owned, support_impulse_body, reaction) = if no_impulse && bilateral_rows.is_none() {
        // Admit retained-inertia rows without repairing their finite residuals.
        // Passive release retains the exact original state and its FK.
        let reaction = passive_reaction(&model_after, state, constraints, &mut cancelled)?;
        (None, Wrench::default(), reaction)
    } else {
        let mut body_impulses = vec![Wrench::default(); model.link_count()];
        body_impulses[target] = impulse_on_retained_link;
        let (response, reaction) = if let Some((loops, config)) = bilateral_rows {
            let (rows, constraint_controls) = constraints.ok_or(PayloadError::Invalid(
                "missing supported mixed constraint controls",
            ))?;
            let event = bilateral::native_response(
                &model_after,
                state,
                &body_impulses,
                loops,
                rows,
                constraint_controls,
                config,
                &mut cancelled,
            )?;
            let response = event.impulse.clone();
            bilateral_event = Some(event);
            (response, SupportedPayloadConstraintReaction::empty())
        } else {
            impulse_response(
                &model_after,
                state,
                &zero,
                &body_impulses,
                constraints,
                &mut cancelled,
            )?
        };
        poll(&mut cancelled)?;
        if no_impulse {
            // Full mixed mobility was admitted above, without projecting the
            // finite tangent residual. Preserve the passive state's bits.
            (None, Wrench::default(), reaction)
        } else {
            state_after.joint_velocity = response.joint_velocity_after;
            let delta = solve_floating_base_system(
                &payload_before.inertia.matrix()?,
                impulse_on_payload_body,
            )?;
            payload_after.twist_body = payload_before.twist_body.plus(delta);
            validate_twist(
                payload_after.twist_body,
                "supported_payload_release.payload_velocity",
            )?;
            if payload_after.twist_body == payload_before.twist_body {
                return Err(PayloadError::Invalid(
                    "supported payload velocity jump erased",
                ));
            }
            poll(&mut cancelled)?;
            let after = forward_kinematics(
                &model_after,
                state_after.base,
                &state_after.joint_position,
                &state_after.joint_velocity,
            )?;
            poll(&mut cancelled)?;
            (Some(after), response.support_impulse_body, reaction)
        }
    };
    let after = after_owned.as_ref().unwrap_or(&kin);
    let mut momentum_after = momentum_from_kinematics(&model_after, after)?;
    let payload_momentum = payload_after.inertia.momentum(payload_after.twist_body)?;
    let payload_world = payload_momentum.transform_by(&payload_after.world_from_body);
    momentum_after.linear_momentum_world =
        momentum_after.linear_momentum_world + payload_world.force;
    momentum_after.angular_momentum_world =
        momentum_after.angular_momentum_world + payload_world.torque;
    momentum_after.kinetic_energy = finite(
        momentum_after.kinetic_energy
            + payload_after
                .inertia
                .kinetic_energy(payload_after.twist_body)?,
        "supported_payload_release.kinetic_after",
    )?;
    validate_wrench(
        Wrench::new(
            momentum_after.angular_momentum_world,
            momentum_after.linear_momentum_world,
        ),
        "supported_payload_release.momentum_after",
    )?;
    poll(&mut cancelled)?;

    // Preserve the initial loaded-to-separated decomposition defect before any
    // motion response. It is neither an applied impulse nor repaired inventory.
    let target_pose = kin.world_from_link[target];
    let target_twist = kin.body_twist[target];
    let loaded = model.links()[target].inertia();
    let retained = model_after.links()[target].inertia();
    let loaded_momentum = loaded.momentum(target_twist)?;
    let retained_momentum = retained.momentum(target_twist)?;
    let payload_initial_momentum = payload_before.inertia.momentum(payload_before.twist_body)?;
    let mut change = SpatialSum::default();
    for (value, pose) in [
        (negate(loaded_momentum), target_pose),
        (retained_momentum, target_pose),
        (payload_initial_momentum, payload_before.world_from_body),
    ] {
        change.add_world_wrench(
            pose.rotation().rotate(value.force)?,
            pose.rotation().rotate(value.torque)?,
            &[pose.translation()],
        );
    }
    let incoming = payload_initial_momentum
        .transform_by(&attachment.link_from_payload)
        .to_array();
    let old = loaded_momentum.to_array();
    let new = retained_momentum.to_array();
    let mut local: [Dd; 6] = core::array::from_fn(|i| {
        Dd::from_f64(new[i]) + Dd::from_f64(incoming[i]) - Dd::from_f64(old[i])
    });
    let mut projected = vec![Wrench::default(); model.link_count()];
    let boost = account_boost(model, &kin)?;
    let (bilateral_reaction, bilateral_work) = if let Some(event) = bilateral_event {
        let (_, constraint_controls) = constraints.ok_or(PayloadError::Invalid(
            "missing supported mixed constraint controls",
        ))?;
        let (reaction, work) = bilateral::finish_reaction(
            &model_after,
            state,
            &kin,
            after,
            &state_after.joint_velocity,
            boost,
            event,
            constraint_controls,
            &mut cancelled,
        )?;
        (Some(reaction), Some(work))
    } else {
        (None, None)
    };
    let centered_target = twist_in_frame(target_twist, target_pose, boost)?;
    let centered_payload = twist_in_frame(
        payload_before.twist_body,
        payload_before.world_from_body,
        boost,
    )?;
    let decomposition = |target_velocity, payload_velocity| -> Result<Dd, PayloadError> {
        let change = Dd::from_f64(retained.kinetic_energy(target_velocity)?)
            + Dd::from_f64(payload_before.inertia.kinetic_energy(payload_velocity)?)
            - Dd::from_f64(loaded.kinetic_energy(target_velocity)?);
        finite(
            change.to_f64(),
            "supported_payload_release.energy_decomposition",
        )?;
        Ok(change)
    };
    let mut kinetic_change = decomposition(target_twist, payload_before.twist_body)?;
    let mut centered_change = decomposition(centered_target, centered_payload)?;
    let mut supplied_work = Dd::ZERO;
    if !no_impulse {
        for (index, link) in model_after.links().iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let pose = kin.world_from_link[index];
            let old = kin.body_twist[index];
            let new = after.body_twist[index];
            let actual = link.inertia().momentum(new.minus(old))?;
            change.add_world_wrench(
                pose.rotation().rotate(actual.force)?,
                pose.rotation().rotate(actual.torque)?,
                &[pose.translation()],
            );
            if index == target {
                for (sum, value) in local.iter_mut().zip(actual.to_array()) {
                    *sum = *sum + Dd::from_f64(value);
                }
            } else {
                projected[index] = negate(actual);
            }
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
        }
        let payload_delta = payload_after.twist_body.minus(payload_before.twist_body);
        let actual = payload_before.inertia.momentum(payload_delta)?;
        let pose = payload_before.world_from_body;
        change.add_world_wrench(
            pose.rotation().rotate(actual.force)?,
            pose.rotation().rotate(actual.torque)?,
            &[pose.translation()],
        );
        let target_actual = actual.transform_by(&attachment.link_from_payload);
        for (sum, value) in local.iter_mut().zip(target_actual.to_array()) {
            *sum = *sum + Dd::from_f64(value);
        }
        kinetic_change = kinetic_change
            + kinetic_link_increment(
                payload_before.inertia,
                pose,
                payload_before.twist_body,
                pose,
                payload_after.twist_body,
            )?;
        centered_change = centered_change
            + kinetic_link_increment(
                payload_before.inertia,
                pose,
                centered_payload,
                pose,
                twist_in_frame(payload_after.twist_body, pose, boost)?,
            )?;
        supplied_work =
            (sum_products(
                impulse_on_retained_link.to_array(),
                after.body_twist[target].minus(target_twist).to_array(),
            ) + sum_products(impulse_on_payload_body.to_array(), payload_delta.to_array()))
                * Dd::from_f64(0.5);
    }
    let local = local.map(Dd::to_f64);
    let local = Wrench::new(
        Vec3::new(local[0], local[1], local[2]),
        Vec3::new(local[3], local[4], local[5]),
    );
    validate_wrench(local, "supported_payload_release.local_decomposition")?;
    projected[target] = negate(local);
    poll(&mut cancelled)?;
    let mut joint_momentum_residual = inverse_dynamics(
        &model_after,
        BaseState::stationary(state.base.world_from_base),
        &state.joint_position,
        &zero,
        &zero,
        Vec3::default(),
        &projected,
    )?
    .generalized_force;
    let constraint_joint_impulse = bilateral_reaction
        .as_ref()
        .map_or(reaction.joint_impulse.as_slice(), |reaction| {
            &reaction.generalized_constraint_impulse[6..]
        });
    for (index, (residual, &applied)) in joint_momentum_residual
        .iter_mut()
        .zip(constraint_joint_impulse)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        *residual = finite(
            (Dd::from_f64(*residual) - Dd::from_f64(applied)).to_f64(),
            "supported_payload_release.constraint_momentum_residual",
        )?;
    }
    poll(&mut cancelled)?;
    let mut support = SpatialSum::default();
    support.add_world_wrench(
        state
            .base
            .world_from_base
            .rotation()
            .rotate(support_impulse_body.force)?,
        state
            .base
            .world_from_base
            .rotation()
            .rotate(support_impulse_body.torque)?,
        &[state.base.world_from_base.translation()],
    );
    let physical = change.checked_wrench("supported_payload_release.change")?;
    let residual = change
        .subtract(support)
        .checked_wrench("supported_payload_release.residual")?;
    let support_impulse_world = support.checked_wrench("supported_payload_release.support")?;
    let maximum_normalized_momentum_residual =
        audit_momentum(residual, &joint_momentum_residual, shared_controls)?;
    let supplied_work_j = finite(
        supplied_work.to_f64(),
        "supported_payload_release.supplied_work",
    )?;
    if supplied_work_j < 0.0 {
        return Err(PayloadError::Invalid("negative supported separation work"));
    }
    if supplied_work_j > maximum_work_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: supplied_work_j,
            maximum_j: maximum_work_j,
        });
    }
    let support_work = sum_products(
        support_impulse_body.to_array(),
        state.base.twist_body.to_array(),
    );
    let centered_support_work = sum_products(
        support_impulse_body.to_array(),
        twist_in_frame(state.base.twist_body, state.base.world_from_base, boost)?.to_array(),
    );
    // Retain the native reaction pairing before rounding its public receipt.
    // Small signed work must survive cancellation with support/source work.
    let constraint_work = bilateral_work.unwrap_or_else(|| {
        (sum_products(
            reaction.joint_impulse.iter().copied(),
            state.joint_velocity.iter().copied(),
        ) + sum_products(
            reaction.joint_impulse.iter().copied(),
            state_after.joint_velocity.iter().copied(),
        )) * Dd::from_f64(0.5)
    });
    let energy_residual_j = finite(
        (centered_change - centered_support_work - supplied_work - constraint_work).to_f64(),
        "supported_payload_release.energy_residual",
    )?;
    if energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(PayloadError::EnergyResidual {
            residual_j: energy_residual_j,
            maximum_j: controls.energy_tolerance_j,
        });
    }
    let kinetic_before_j = energy_in_frame(model, &kin, None, boost, &mut cancelled)?;
    let kinetic_after_j = energy_in_frame(
        &model_after,
        after,
        Some(&payload_after),
        boost,
        &mut cancelled,
    )?;
    let balance = SupportedPayloadBalance {
        momentum_before,
        momentum_after,
        linear_momentum_change_world: physical.force,
        angular_momentum_change_world: physical.torque,
        support_impulse_body,
        support_impulse_world,
        linear_momentum_residual_world: residual.force,
        angular_momentum_residual_world: residual.torque,
        joint_momentum_residual,
        maximum_normalized_momentum_residual,
        kinetic_energy_change_j: finite(
            kinetic_change.to_f64(),
            "supported_payload_release.kinetic_change",
        )?,
        support_work_j: finite(
            support_work.to_f64(),
            "supported_payload_release.support_work",
        )?,
        accounting_boost_world: boost,
        kinetic_before_j,
        kinetic_after_j,
        accounting_kinetic_energy_change_j: finite(
            centered_change.to_f64(),
            "supported_payload_release.centered_change",
        )?,
        accounting_support_work_j: finite(
            centered_support_work.to_f64(),
            "supported_payload_release.centered_support_work",
        )?,
        energy_residual_j,
    };
    poll(&mut cancelled)?;
    Ok((
        ConstrainedSupportedPayloadRelease {
            release: SupportedPayloadRelease {
                model_after,
                state_after,
                payload_after,
                supplied_work_j,
                applied_impulse: impulse,
                impulse_on_payload_body,
                impulse_on_retained_link,
                balance,
            },
            reaction,
        },
        bilateral_reaction,
    ))
}

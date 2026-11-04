//! Simultaneous inventory separation with one retained-mechanism response.
//!
//! Handle order validates inventory only. Every outgoing body inherits the
//! original loaded velocity field, and every retained-side impulse acts on the
//! same final retained inertia and velocity. Shared recoil work is not split
//! into a sequence of physically different releases.

use super::super::{
    PayloadAttachment, PayloadBalance, PayloadControls, PayloadError, PayloadState, account_boost,
    bilateral, energy_in_frame, energy_in_frame_many, finish_common_frame_balance, finite,
    payload_at_link, poll, system_totals, system_totals_many, validate_inputs,
};
use super::{
    SeparationImpulse, add_joint_delta, admit_attachment_link, admit_separation_impulse, centered,
    separation_impulses,
};
use crate::articulated::motion::reconfigure::constraints::{
    JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
};
use crate::articulated::motion::reconfigure::state_kinematics;
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralImpulseResponse, BilateralLoop, joint,
};
use crate::articulated::motion::{FreeFloatingMotionState, SpatialSum};
use crate::articulated::{
    ArticulatedModel, free_floating_impulse_delta, solve_floating_base_system, validate_twist,
};
use fs_ga::{Twist, Vec3, Wrench};
use fs_math::dd::Dd;

/// One captured inventory to detach in this simultaneous event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PayloadReleaseInput<'a> {
    /// Existing identity from capture. On a repeated target, consume handles
    /// in reverse capture order; independent targets may be interleaved.
    pub attachment: &'a PayloadAttachment,
    /// Optional point impulse and free couple delivered to this payload.
    pub impulse: Option<SeparationImpulse>,
}

/// One detached body and its contribution to the simultaneous source work.
#[derive(Debug, Clone, PartialEq)]
pub struct ReleasedPayload {
    /// Actual independent body, retaining its original frame and full inertia.
    pub payload_after: PayloadState,
    /// Original prescribed impulse, including its shared physical point.
    pub applied_impulse: Option<SeparationImpulse>,
    /// Applied impulse about the payload origin, expressed in payload axes.
    pub impulse_on_payload_body: Wrench,
    /// Equal/opposite impulse about its retained target link origin.
    pub impulse_on_retained_link: Wrench,
    /// Signed pair work against the ONE final retained state, in the common
    /// translating frame [J]. Shared recoil cross terms are included.
    pub supplied_work_j: f64,
}

/// All surviving bodies after one native simultaneous separation response.
#[derive(Debug)]
pub struct PayloadBatchRelease {
    /// Original topology and exact retained inertias after every detachment.
    pub model_after: ArticulatedModel,
    /// One retained-tree state after all sources and constraint reactions.
    pub state_after: FreeFloatingMotionState,
    /// Detached bodies and source accounts in input order.
    pub payloads_after: Vec<ReleasedPayload>,
    /// Compensated signed work of all separation pairs [J]. The existing
    /// absolute source-work cap applies to this total, not each rounded entry.
    pub supplied_work_j: f64,
    /// Independently capped actual ideal-constraint work [J], never heat.
    pub constraint_work_j: f64,
    /// Complete physical and projected momentum, common-frame energy and work.
    pub balance: PayloadBalance,
}

/// Detach a nonempty collection of captured payloads at one physical instant.
///
/// Every payload starts from the original loaded tree's actual velocity field.
/// All retained inertias are restored before the accumulated separation
/// impulses and existing gear/lock rows receive ONE native coupled response.
/// Each detached body then responds with its own original full spatial inertia.
/// Only the complete physical endpoint receives native speed admission.
///
/// On each repeated target, handles must appear in reverse capture order. This
/// validates exact inventory prefixes without solving intermediate motion;
/// independent targets may be interleaved or permuted. Unknown, repeated,
/// skipped or stale handles refuse without publishing a partial result.
///
/// Passive releases preserve the retained state bits after row/rank admission.
/// Exactly cancelling retained-side sources also leave that tree passive while
/// their individual payload impulses still act. Signed pair work uses the same
/// final retained velocity, preserving shared recoil cross terms. The existing
/// absolute total source-work and separate constraint-work limits still apply.
/// No capture loss, fracture law, contact discovery or thermal credit is inferred.
///
/// Assembly and accounting outside the original constraint solve use
/// O(links + payloads) work/storage, with cancellation polls every 64 entries
/// and between bounded native solves. All inputs remain unchanged on refusal.
///
/// # Errors
/// Empty or invalid inventory, supplied impulse, state, constraint regularity,
/// native endpoint, work/energy/momentum admission, or cancellation refusal.
pub fn release_payloads(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadReleaseInput<'_>],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PayloadBatchRelease, PayloadError> {
    release_evaluated(model, state, inputs, constraints, controls, None, cancelled)
        .map(|(release, _)| release)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(in crate::articulated::motion::inventory) fn release_evaluated(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadReleaseInput<'_>],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PayloadBatchRelease, Option<BilateralImpulseResponse>), PayloadError> {
    poll(&mut cancelled)?;
    let first = inputs
        .first()
        .ok_or(PayloadError::Invalid("empty payload release"))?;
    let before = validate_inputs(
        model,
        state,
        first.attachment.target_link,
        constraints,
        controls,
        &mut cancelled,
    )?;
    let mut links = model.links().to_vec();
    let mut payloads_before = Vec::with_capacity(inputs.len());
    let mut payloads_after = Vec::with_capacity(inputs.len());
    let mut sums = Vec::with_capacity(model.link_count());
    for index in 0..model.link_count() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        sums.push(SpatialSum::default());
    }
    for (index, input) in inputs.iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let attachment = input.attachment;
        let target = attachment.target_link;
        let link = links
            .get(target)
            .ok_or(PayloadError::Invalid("release target link"))?;
        admit_attachment_link(link, attachment)?;
        let prescribed = input.impulse.unwrap_or_default();
        admit_separation_impulse(prescribed)?;
        let payload = payload_at_link(
            &before,
            target,
            attachment.link_from_payload,
            attachment.payload_inertia,
        )?;
        let (retained, detached) = separation_impulses(prescribed, attachment.link_from_payload)?;
        sums[target].add_world_wrench(retained.force, retained.torque, &[]);
        payloads_before.push(payload);
        payloads_after.push(ReleasedPayload {
            payload_after: payload,
            applied_impulse: input.impulse,
            impulse_on_payload_body: detached,
            impulse_on_retained_link: retained,
            supplied_work_j: 0.0,
        });
        // Inventory-only fold. No intermediate model or motion is admitted.
        links[target] = attachment.retained_link.clone();
    }
    let model_after = ArticulatedModel::new(links)?;
    poll(&mut cancelled)?;
    let mut body_impulses = Vec::with_capacity(sums.len());
    let mut tree_has_impulse = false;
    for (index, sum) in sums.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let impulse = sum.checked_wrench("simultaneous release link impulse")?;
        tree_has_impulse |= impulse != Wrench::default();
        body_impulses.push(impulse);
    }
    let before_totals = system_totals(model, state, &before, None, &mut cancelled)?;
    let boost = account_boost(model, &before)?;
    let kinetic_before_j = energy_in_frame(model, &before, None, boost, &mut cancelled)?;
    let mut state_after = state.clone();
    let mut multipliers = vec![0.0; constraints.len()];
    let mut generalized_constraint_impulse = vec![0.0; 6 + model.dof_count()];
    let mut bilateral_reaction = None;
    if let Some((loops, config)) = bilateral_rows {
        let reaction = joint::apply_impulses(
            &model_after,
            state,
            &vec![0.0; model.dof_count()],
            &body_impulses,
            loops,
            constraints,
            config,
            &mut cancelled,
        )?;
        bilateral::admit_condition(reaction.reciprocal_condition, controls)?;
        if tree_has_impulse {
            state_after.base = reaction.impulse.base_after;
            state_after.joint_velocity = reaction.impulse.joint_velocity_after.clone();
        }
        generalized_constraint_impulse = bilateral::project_constraint_impulses(
            &model_after,
            state,
            &reaction.constraint_body_impulses,
            &reaction.constraint_joint_impulses,
            &mut cancelled,
        )?;
        multipliers = reaction.joint_multipliers.clone();
        bilateral_reaction = Some(reaction);
    } else {
        let mut prepared = if constraints.is_empty() {
            None
        } else {
            Some(PreparedJointConstraints::new(
                &model_after,
                &state.joint_position,
                constraints,
                controls.minimum_reciprocal_condition,
                &mut cancelled,
            )?)
        };
        if tree_has_impulse {
            let delta = free_floating_impulse_delta(
                &model_after,
                &state.joint_position,
                &vec![0.0; model.dof_count()],
                &body_impulses,
            )?;
            poll(&mut cancelled)?;
            state_after.base.twist_body = state_after.base.twist_body.plus(delta.base);
            add_joint_delta(
                &mut state_after.joint_velocity,
                &delta.joint,
                &mut cancelled,
            )?;
            if let Some(mobility) = &mut prepared {
                let rhs = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?
                    .into_iter()
                    .map(|rate| -rate)
                    .collect::<Vec<_>>();
                let reaction = mobility.solve_response(&rhs, &mut cancelled)?;
                state_after.base.twist_body = state_after.base.twist_body.plus(reaction.base_delta);
                add_joint_delta(
                    &mut state_after.joint_velocity,
                    &reaction.joint_delta,
                    &mut cancelled,
                )?;
                generalized_constraint_impulse[6..].copy_from_slice(&reaction.joint_effort);
                multipliers = reaction.multipliers;
            }
        }
    }
    for (index, (input, output)) in inputs.iter().zip(&mut payloads_after).enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let prescribed = input.impulse.unwrap_or_default();
        if prescribed.impulse_on_payload_link_ns != Vec3::default()
            || prescribed.couple_on_payload_link_nms != Vec3::default()
        {
            let delta = solve_floating_base_system(
                &output.payload_after.inertia.matrix()?,
                output.impulse_on_payload_body,
            )?;
            poll(&mut cancelled)?;
            output.payload_after.twist_body = output.payload_after.twist_body.plus(delta);
            validate_twist(
                output.payload_after.twist_body,
                "payload_release.payload_velocity",
            )?;
        }
    }
    admit_constraints(
        constraints,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut cancelled,
    )?;
    let after = state_kinematics(&model_after, &state_after)?;
    poll(&mut cancelled)?;
    let after_totals = system_totals_many(
        &model_after,
        &state_after,
        &after,
        inputs
            .iter()
            .zip(&payloads_after)
            .map(|(input, output)| (input.attachment.target_link, &output.payload_after)),
        &mut cancelled,
    )?;
    let kinetic_after_j = energy_in_frame_many(
        &model_after,
        &after,
        payloads_after.iter().map(|output| &output.payload_after),
        boost,
        &mut cancelled,
    )?;
    let mut source_work = Dd::ZERO;
    for (index, ((input, payload_before), output)) in inputs
        .iter()
        .zip(&payloads_before)
        .zip(&mut payloads_after)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let target = input.attachment.target_link;
        let pair = mean_work(
            output.impulse_on_retained_link,
            centered(
                before.body_twist[target],
                before.world_from_link[target],
                boost,
            )?,
            centered(
                after.body_twist[target],
                after.world_from_link[target],
                boost,
            )?,
        ) + mean_work(
            output.impulse_on_payload_body,
            centered(
                payload_before.twist_body,
                payload_before.world_from_body,
                boost,
            )?,
            centered(
                output.payload_after.twist_body,
                output.payload_after.world_from_body,
                boost,
            )?,
        );
        output.supplied_work_j = finite(pair.to_f64(), "release pair midpoint work")?;
        source_work = source_work + pair;
    }
    let supplied_work_j = finite(source_work.to_f64(), "simultaneous release midpoint work")?;
    if supplied_work_j.abs() > controls.maximum_energy_exchange_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: supplied_work_j,
            maximum_j: controls.maximum_energy_exchange_j,
        });
    }
    let constraint_work_j = if let Some(reaction) = &bilateral_reaction {
        bilateral::constraint_work(
            &reaction.constraint_body_impulses,
            &reaction.constraint_joint_impulses,
            &before,
            &after,
            &state.joint_velocity,
            &state_after.joint_velocity,
            boost,
            &mut cancelled,
        )?
    } else {
        let rates_before = apply_rows(constraints, &state.joint_velocity, &mut cancelled)?;
        let rates_after = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?;
        let mut work = 0.0;
        for ((&multiplier, &initial), &final_rate) in
            multipliers.iter().zip(&rates_before).zip(&rates_after)
        {
            work = finite(
                multiplier.mul_add(0.5 * initial + 0.5 * final_rate, work),
                "release constraint work",
            )?;
        }
        work
    };
    let energy_residual_j = finite(
        kinetic_after_j - kinetic_before_j - supplied_work_j - constraint_work_j,
        "release common-frame energy residual",
    )?;
    let balance = finish_common_frame_balance(
        before_totals,
        after_totals,
        generalized_constraint_impulse,
        multipliers,
        constraint_work_j,
        boost,
        kinetic_before_j,
        kinetic_after_j,
        energy_residual_j,
        controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok((
        PayloadBatchRelease {
            model_after,
            state_after,
            payloads_after,
            supplied_work_j,
            constraint_work_j,
            balance,
        },
        bilateral_reaction,
    ))
}

fn mean_work(impulse: Wrench, before: Twist, after: Twist) -> Dd {
    let mut work = Dd::ZERO;
    for ((j, old), new) in impulse
        .to_array()
        .into_iter()
        .zip(before.to_array())
        .zip(after.to_array())
    {
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        work = work + Dd::from_f64(j) * mean;
    }
    work
}

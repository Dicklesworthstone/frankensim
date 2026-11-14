//! Simultaneous perfectly inelastic capture on one original articulated tree.
//!
//! The inventory folds below construct release identities only. No intermediate
//! motion is solved or accepted: every arriving body contributes to the same
//! loaded-inertia response and the same final relative-energy quadratic.

use super::{
    PayloadAttachment, PayloadBalance, PayloadControls, PayloadError, PayloadState, account_boost,
    accumulate_body_change, attachment_decomposition, bilateral, combine_inertias, energy_in_frame,
    energy_in_frame_many, finish_common_frame_balance, finite, payload_at_link, poll,
    system_totals, system_totals_many, twist_in_frame, validate_inputs,
};
use crate::MomentumWork;
use crate::articulated::motion::reconfigure::constraints::{
    JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
};
use crate::articulated::motion::reconfigure::state_kinematics;
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralLoop, joint,
};
use crate::articulated::motion::{FreeFloatingMotionState, SpatialSum};
use crate::articulated::{
    ArticulatedModel, FreeFloatingBaseState, Link, free_floating_impulse_delta, validate_twist,
};
use fs_math::dd::Dd;

/// One independent arriving body and its link in the original target tree.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PayloadCaptureInput {
    /// Existing link index; repeated targets are admitted.
    pub target_link: usize,
    /// Actual independent body pose, full inertia and velocity at this event.
    pub payload: PayloadState,
}

/// One simultaneous capture, with every result in the authored input order.
#[derive(Debug, Clone, PartialEq)]
pub struct PayloadBatchCapture {
    /// Original topology and joints with all target inertias loaded together.
    pub model_after: ArticulatedModel,
    /// One coupled native recoil solution at the unchanged configuration.
    pub state_after: FreeFloatingMotionState,
    /// Existing release identities in input order. For repeated target links,
    /// release in REVERSE input order on that link; different links are independent.
    /// The prefix inertias bind inventory only, never intermediate physical states.
    pub attachments: Vec<PayloadAttachment>,
    /// Each captured body's actual rigid motion in the returned tree state.
    pub payloads_after: Vec<PayloadState>,
    /// Nonnegative total relative-motion loss, including every shared recoil [J].
    pub removed_relative_energy_j: f64,
    /// Whole-system physical and native projected-momentum/work accounting.
    pub balance: PayloadBalance,
}

/// Capture a nonempty collection of finite rigid bodies simultaneously.
///
/// All attachments use the original accepted tree pose. Their full centroidal
/// tensors, orientations, COM offsets and momenta enter ONE loaded model and
/// ONE native free-floating inverse-mass response. Existing constant affine
/// gear/lock rows may react on that same loaded inertia before speed admission.
/// No sequence of capture solves or dense generalized mass matrix is introduced.
///
/// Repeated targets build a stack of exact existing release handles in input
/// order; only their inventory identities are sequential. Captures are physical
/// simultaneous events, and their final rates/loss agree under permutations to
/// floating-point accuracy. Releases on one target must consume its handles in
/// reverse input order. Unrelated links may release independently.
///
/// Incoming and retained link momenta use the shared compensated spatial sum.
/// Solving in a common translating inertial frame retains small recoil under a
/// large common velocity. Removed energy is the sum of original-body velocity-
/// change quadratics evaluated against the final coupled state. The common-frame
/// energy audit remains separate; no residual is assigned to heat.
///
/// Work/storage outside original constraint handling is O(links + payloads).
/// One final native model construction replaces repeated full tree copies.
/// Extra loops poll every 64 entries; native passes retain the 4096-link cap.
/// Any geometry, inventory, momentum, work, speed or cancellation refusal leaves
/// all inputs unchanged and publishes no model, state or partial attachment.
///
/// Contact discovery, finite grasp duration, overlapping-volume validity,
/// independent moving articulations and thermal conversion are caller-owned.
pub fn capture_payloads(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadCaptureInput],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PayloadBatchCapture, PayloadError> {
    capture_evaluated(model, state, inputs, constraints, controls, None, cancelled)
        .map(|result| result.capture)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn capture_evaluated(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadCaptureInput],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<bilateral::BilateralPayloadBatchCapture, PayloadError> {
    poll(&mut cancelled)?;
    let first = inputs
        .first()
        .ok_or(PayloadError::Invalid("empty payload capture"))?;
    let before_kin = validate_inputs(
        model,
        state,
        first.target_link,
        constraints,
        controls,
        &mut cancelled,
    )?;
    let mut links = model.links().to_vec();
    let mut attachments = Vec::with_capacity(inputs.len());
    for (index, input) in inputs.iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let target = input.target_link;
        let retained = links
            .get(target)
            .ok_or(PayloadError::Invalid("payload target link index"))?;
        validate_twist(input.payload.twist_body, "payload.twist_body")?;
        let link_from_payload = before_kin.world_from_link[target]
            .inverse()?
            .compose(input.payload.world_from_body)?;
        let loaded_inertia =
            combine_inertias(retained.inertia(), input.payload.inertia, link_from_payload)?;
        attachments.push(PayloadAttachment {
            target_link: target,
            retained_link: retained.clone(),
            loaded_inertia,
            link_from_payload,
            payload_inertia: input.payload.inertia,
        });
        links[target] = Link::new(
            retained.name(),
            retained.parent(),
            retained.parent_from_child_zero(),
            retained.joint(),
            loaded_inertia,
        );
    }
    let model_after = ArticulatedModel::new(links)?;
    poll(&mut cancelled)?;
    let before = system_totals_many(
        model,
        state,
        &before_kin,
        inputs
            .iter()
            .map(|input| (input.target_link, &input.payload)),
        &mut cancelled,
    )?;
    let boost = account_boost(model, &before_kin)?;
    let mut sums: Vec<_> = (0..model.link_count())
        .map(|_| SpatialSum::default())
        .collect();
    for (index, (input, attachment)) in inputs.iter().zip(&attachments).enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let payload = input.payload;
        let carried = payload_at_link(
            &before_kin,
            input.target_link,
            attachment.link_from_payload,
            payload.inertia,
        )?;
        // Resolve each weak mismatch before applying inertia. The same
        // payload axes and boost subtraction preserve an exactly carried
        // input as zero even when finite tangent defects are admitted.
        let incoming = twist_in_frame(payload.twist_body, payload.world_from_body, boost)?;
        let carried = twist_in_frame(carried.twist_body, payload.world_from_body, boost)?;
        let mismatch = payload
            .inertia
            .momentum(incoming.minus(carried))?
            .transform_by(&attachment.link_from_payload);
        sums[input.target_link].add_world_wrench(mismatch.force, mismatch.torque, &[]);
    }
    let mut body_impulses = Vec::with_capacity(sums.len());
    for (index, sum) in sums.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        body_impulses.push(sum.checked_wrench("simultaneous capture momentum mismatch")?);
    }
    let (state_after, impulse, multipliers, reaction) = if let Some((loops, config)) =
        bilateral_rows
    {
        poll(&mut cancelled)?;
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
        poll(&mut cancelled)?;
        let state_after = FreeFloatingMotionState {
            base: reaction.impulse.base_after,
            joint_position: state.joint_position.clone(),
            joint_velocity: reaction.impulse.joint_velocity_after.clone(),
        };
        let impulse = bilateral::project_constraint_impulses(
            &model_after,
            state.base.world_from_base,
            &state.joint_position,
            &reaction.constraint_body_impulses,
            &reaction.constraint_joint_impulses,
            &mut cancelled,
        )?;
        let multipliers = reaction.joint_multipliers.clone();
        (state_after, impulse, multipliers, Some(reaction))
    } else {
        // The loaded tree takes only the mismatch impulse, from the state's
        // own velocity: an exactly carried payload changes no velocity, and the
        // finite assembly defect stays in the inventory rather than becoming a
        // rounded velocity jump that would be counted as exchanged energy.
        let solved = free_floating_impulse_delta(
            &model_after,
            &state.joint_position,
            &vec![0.0; model.dof_count()],
            &body_impulses,
        )?;
        poll(&mut cancelled)?;
        let mut joint_velocity = Vec::with_capacity(solved.joint.len());
        for (index, (rate, delta)) in state.joint_velocity.iter().zip(&solved.joint).enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            joint_velocity.push(finite(rate + delta, "payload joint rate")?);
        }
        let mut state_after = FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                state.base.world_from_base,
                state.base.twist_body.plus(solved.base),
            ),
            joint_position: state.joint_position.clone(),
            joint_velocity,
        };
        let mut impulse = vec![0.0; 6 + model.dof_count()];
        let mut multipliers = Vec::new();
        if !constraints.is_empty() {
            let mut prepared = PreparedJointConstraints::new(
                &model_after,
                &state.joint_position,
                constraints,
                controls.minimum_reciprocal_condition,
                &mut cancelled,
            )?;
            let rhs: Vec<_> = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?
                .into_iter()
                .map(|value| -value)
                .collect();
            let response = prepared.solve_response(&rhs, &mut cancelled)?;
            state_after.base.twist_body = state_after.base.twist_body.plus(response.base_delta);
            for (index, delta) in response.joint_delta.iter().enumerate() {
                if index % 64 == 0 {
                    poll(&mut cancelled)?;
                }
                state_after.joint_velocity[index] = finite(
                    state_after.joint_velocity[index] + delta,
                    "payload final joint rate",
                )?;
            }
            impulse[6..].copy_from_slice(&response.joint_effort);
            multipliers = response.multipliers;
        }
        (state_after, impulse, multipliers, None)
    };
    let final_rows = admit_constraints(
        constraints,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut cancelled,
    )?;
    let after_kin = state_kinematics(&model_after, &state_after)?;
    poll(&mut cancelled)?;
    let mut payloads_after = Vec::with_capacity(inputs.len());
    for (index, attachment) in attachments.iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        payloads_after.push(payload_at_link(
            &after_kin,
            attachment.target_link,
            attachment.link_from_payload,
            attachment.payload_inertia,
        )?);
    }
    let after = system_totals(&model_after, &state_after, &after_kin, None, &mut cancelled)?;
    let mut removed = Dd::ZERO;
    let mut kinetic_change = MomentumWork::default();
    let mut decomposition = MomentumWork::default();
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        removed =
            removed
                + Dd::from_f64(link.inertia().kinetic_energy(
                    before_kin.body_twist[index].minus(after_kin.body_twist[index]),
                )?);
        accumulate_body_change(
            &mut kinetic_change,
            link.inertia(),
            twist_in_frame(
                before_kin.body_twist[index],
                before_kin.world_from_link[index],
                boost,
            )?,
            twist_in_frame(
                after_kin.body_twist[index],
                after_kin.world_from_link[index],
                boost,
            )?,
        );
    }
    for (index, ((input, attached), attachment)) in inputs
        .iter()
        .zip(&payloads_after)
        .zip(&attachments)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        removed = removed
            + Dd::from_f64(
                input
                    .payload
                    .inertia
                    .kinetic_energy(input.payload.twist_body.minus(attached.twist_body))?,
            );
        let centered_after = twist_in_frame(attached.twist_body, attached.world_from_body, boost)?;
        accumulate_body_change(
            &mut kinetic_change,
            input.payload.inertia,
            twist_in_frame(
                input.payload.twist_body,
                input.payload.world_from_body,
                boost,
            )?,
            centered_after,
        );
        decomposition.merge(attachment_decomposition(
            attachment,
            twist_in_frame(
                after_kin.body_twist[input.target_link],
                after_kin.world_from_link[input.target_link],
                boost,
            )?,
            centered_after,
        )?);
    }
    // The inventory folds telescope at ONE accepted endpoint, including
    // repeated targets. Their finite assembly defect is separate from loss.
    kinetic_change.subtract(decomposition);
    let removed_relative_energy_j = finite(removed.to_f64(), "payload removed relative energy")?;
    if removed_relative_energy_j < 0.0 {
        return Err(PayloadError::Invalid("negative capture relative energy"));
    }
    if removed_relative_energy_j > controls.maximum_energy_exchange_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: removed_relative_energy_j,
            maximum_j: controls.maximum_energy_exchange_j,
        });
    }
    let reaction_work_j = if let Some(reaction) = &reaction {
        bilateral::constraint_work(
            &reaction.constraint_body_impulses,
            &reaction.constraint_joint_impulses,
            &before_kin,
            &after_kin,
            &state.joint_velocity,
            &state_after.joint_velocity,
            boost,
            &mut cancelled,
        )?
        .to_f64()
    } else {
        let initial_rows = apply_rows(constraints, &state.joint_velocity, &mut cancelled)?;
        let mut work = 0.0;
        for (index, multiplier) in multipliers.iter().enumerate() {
            poll(&mut cancelled)?;
            work = finite(
                work + multiplier * (0.5 * initial_rows[index] + 0.5 * final_rows.velocity[index]),
                "payload reaction work",
            )?;
        }
        work
    };
    let kinetic_before_j = energy_in_frame_many(
        model,
        &before_kin,
        inputs.iter().map(|input| &input.payload),
        boost,
        &mut cancelled,
    )?;
    let kinetic_after_j = energy_in_frame(&model_after, &after_kin, None, boost, &mut cancelled)?;
    // The loss is the complete passive projection onto the loaded tree and
    // admitted gear rows. Retain finite gear work separately, never as heat.
    let energy_residual_j = finite(
        (kinetic_change.total() + removed).to_f64(),
        "payload common-frame capture energy residual",
    )?;
    let balance = finish_common_frame_balance(
        before,
        after,
        impulse,
        multipliers,
        reaction_work_j,
        boost,
        kinetic_before_j,
        kinetic_after_j,
        kinetic_change.total().to_f64(),
        -decomposition.total().to_f64(),
        energy_residual_j,
        controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(bilateral::BilateralPayloadBatchCapture {
        capture: PayloadBatchCapture {
            model_after,
            state_after,
            attachments,
            payloads_after,
            removed_relative_energy_j,
            balance,
        },
        reaction,
    })
}

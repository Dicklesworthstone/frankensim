//! Release a retained payload with simultaneous articulated and gear recoil.
//!
//! Geometry is unchanged at the event. The detached body inherits the actual
//! link velocity field; an optional internal impulse pair then changes both
//! surviving systems through their native inverse inertias. No contact,
//! fracture, collision geometry, or destination for mechanical work is inferred.

use super::{
    PayloadAttachment, PayloadBalance, PayloadControls, PayloadError, PayloadState, account_boost,
    energy_in_frame, finish_common_frame_balance, finite, payload_at_link, poll,
    replace_link_inertia, system_totals, validate_inputs,
};
use crate::articulated::motion::FreeFloatingMotionState;
use crate::articulated::motion::reconfigure::constraints::{
    JointConstraint, PreparedJointConstraints, admit_constraints, apply_rows,
};
use crate::articulated::motion::reconfigure::state_kinematics;
use crate::articulated::{
    ArticulatedModel, free_floating_impulse_delta, solve_floating_base_system, validate_twist,
    validate_vec3, validate_wrench,
};
use fs_ga::{Se3, Twist, Vec3, Wrench};

/// One equal-and-opposite internal impulse at a shared physical point.
///
/// Every vector is expressed in the target link's axes at the release event.
/// The retained link receives the exact opposite force impulse and free couple.
/// The point-force moment is formed about the link origin and transformed to
/// the payload origin, retaining all orbital and full-tensor contributions.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SeparationImpulse {
    /// Application point measured from the target link origin [m].
    pub point_link_m: Vec3,
    /// Force impulse delivered to the payload [N s].
    pub impulse_on_payload_link_ns: Vec3,
    /// Additional free angular impulse on the payload [N m s].
    pub couple_on_payload_link_nms: Vec3,
}

/// Both actual surviving states and their common-frame physical accounting.
#[derive(Debug)]
pub struct PayloadRelease {
    /// Original tree topology with the exact retained link inertia restored.
    pub model_after: ArticulatedModel,
    /// Accepted retained-tree motion after both separation and gear reactions.
    pub state_after: FreeFloatingMotionState,
    /// Independent released payload, with its original body frame and inertia.
    pub payload_after: PayloadState,
    /// Actual midpoint work of the complete internal separation impulse pair [J].
    /// The pair is evaluated together, including its quadratic recoil work.
    pub supplied_work_j: f64,
    /// Actual finite midpoint work of the simultaneous ideal joint reactions [J].
    /// This is also retained in `balance.reaction_work_j`; no heat is inferred.
    pub constraint_work_j: f64,
    /// Original prescribed impulse, including its shared physical point.
    pub applied_impulse: Option<SeparationImpulse>,
    /// Applied impulse about the payload's origin, expressed in its own axes.
    pub impulse_on_payload_body: Wrench,
    /// Equal/opposite applied impulse about the retained target link origin.
    pub impulse_on_retained_link: Wrench,
    /// Actual world and projected-total momentum, common-frame energy and work.
    pub balance: PayloadBalance,
}

/// Detach an actual captured inventory, optionally supplying internal recoil work.
///
/// The attachment identifies the exact loaded target link. Its name, parent,
/// reference transform, joint and loaded inertia must still match; changed
/// inventory or a repeated release against the unloaded model is refused.
/// Motion and joint coordinates may have evolved since capture. Other links
/// retain their current inertias. Nested captures on the same link therefore
/// release in reverse inventory order; no tensor subtraction is guessed.
///
/// The payload inherits the current rigid velocity field at its actual pose.
/// A prescribed point impulse and free couple act equally and oppositely on
/// payload and retained link. The original articulated ABA supplies all coupled
/// base/joint recoil, while the native six-component solve supplies payload
/// recoil. Optional constant ideal joint constraints react under the RETAINED
/// inertia, before final native speed admission. An inadmissible free candidate
/// cannot reject a physically admissible constrained endpoint.
///
/// No impulse preserves the retained state bits, after row/rank admission. The
/// returned payload is always a separate native rigid state. Work and energy
/// are evaluated on actual outputs in the declared common translating inertial
/// frame; a large uniform velocity is not subtracted through lab energy totals.
/// Total projected momentum pulls the released payload covector back to the
/// target's virtual link at this fixed geometry. The paired separation impulse
/// cancels there; only the reported gear reaction changes that covector.
///
/// `maximum_energy_exchange_j` caps the absolute supplied separation work.
/// `maximum_constraint_work_j` separately caps actual ideal-reaction work.
/// The caller chooses the physical release and its work source; finite-duration
/// separation, fracture, ongoing contact, mass flux and heat partition are
/// outside this instantaneous event.
///
/// # Errors
/// Changed attachment inventory, invalid controls/state/impulse, unresolved
/// constraints, native inertia or endpoint refusal, momentum/energy/work
/// admission failure, or cancellation. No input or partial result is published.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn release_payload(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    attachment: &PayloadAttachment,
    impulse: Option<SeparationImpulse>,
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PayloadRelease, PayloadError> {
    poll(&mut cancelled)?;
    let target = attachment.target_link;
    admit_attachment(model, attachment)?;
    let retained_link = &attachment.retained_link;
    let before = validate_inputs(model, state, target, constraints, controls, &mut cancelled)?;
    let prescribed = impulse.unwrap_or_default();
    validate_vec3(prescribed.point_link_m, "payload_release.point")?;
    validate_vec3(
        prescribed.impulse_on_payload_link_ns,
        "payload_release.force_impulse",
    )?;
    validate_vec3(
        prescribed.couple_on_payload_link_nms,
        "payload_release.free_couple",
    )?;
    let no_impulse = prescribed.impulse_on_payload_link_ns == Vec3::default()
        && prescribed.couple_on_payload_link_nms == Vec3::default();
    let payload_before = payload_at_link(
        &before,
        target,
        attachment.link_from_payload,
        attachment.payload_inertia,
    )?;
    let model_after = replace_link_inertia(model, target, retained_link.inertia())?;
    let before_totals = system_totals(model, state, &before, None, &mut cancelled)?;
    let boost = account_boost(model, &before)?;
    let kinetic_before_j = energy_in_frame(model, &before, None, boost, &mut cancelled)?;
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
    let impulse_at_link = Wrench::new(
        prescribed
            .point_link_m
            .cross(prescribed.impulse_on_payload_link_ns)
            + prescribed.couple_on_payload_link_nms,
        prescribed.impulse_on_payload_link_ns,
    );
    validate_wrench(impulse_at_link, "payload_release.link_impulse")?;
    let impulse_on_retained_link = Wrench::new(
        impulse_at_link.torque.scale(-1.0),
        impulse_at_link.force.scale(-1.0),
    );
    let impulse_on_payload_body =
        impulse_at_link.transform_by(&attachment.link_from_payload.inverse()?);
    validate_wrench(impulse_on_payload_body, "payload_release.payload_impulse")?;

    let mut state_after = state.clone();
    let mut payload_after = payload_before;
    let mut multipliers = vec![0.0; constraints.len()];
    let mut generalized_constraint_impulse = vec![0.0; 6 + model.dof_count()];
    if !no_impulse {
        let mut body_impulses = vec![Wrench::default(); model.link_count()];
        body_impulses[target] = impulse_on_retained_link;
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
        let payload_delta =
            solve_floating_base_system(&payload_after.inertia.matrix()?, impulse_on_payload_body)?;
        poll(&mut cancelled)?;
        payload_after.twist_body = payload_after.twist_body.plus(payload_delta);
        validate_twist(payload_after.twist_body, "payload_release.payload_velocity")?;
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
    admit_constraints(
        constraints,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut cancelled,
    )?;
    let after = state_kinematics(&model_after, &state_after)?;
    poll(&mut cancelled)?;
    let after_totals = system_totals(
        &model_after,
        &state_after,
        &after,
        Some((target, &payload_after)),
        &mut cancelled,
    )?;
    let kinetic_after_j = energy_in_frame(
        &model_after,
        &after,
        Some(&payload_after),
        boost,
        &mut cancelled,
    )?;
    let supplied_work_j = if no_impulse {
        0.0
    } else {
        let midpoint_link = centered(
            before.body_twist[target],
            before.world_from_link[target],
            boost,
        )?
        .scale(0.5)
        .plus(
            centered(
                after.body_twist[target],
                after.world_from_link[target],
                boost,
            )?
            .scale(0.5),
        );
        let midpoint_payload = centered(
            payload_before.twist_body,
            payload_before.world_from_body,
            boost,
        )?
        .scale(0.5)
        .plus(
            centered(
                payload_after.twist_body,
                payload_after.world_from_body,
                boost,
            )?
            .scale(0.5),
        );
        finite(
            impulse_on_retained_link.pairing(midpoint_link)
                + impulse_on_payload_body.pairing(midpoint_payload),
            "release midpoint work",
        )?
    };
    if supplied_work_j.abs() > controls.maximum_energy_exchange_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: supplied_work_j,
            maximum_j: controls.maximum_energy_exchange_j,
        });
    }
    let rates_before = apply_rows(constraints, &state.joint_velocity, &mut cancelled)?;
    let rates_after = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?;
    let mut constraint_work_j = 0.0;
    for ((&multiplier, &initial), &final_rate) in
        multipliers.iter().zip(&rates_before).zip(&rates_after)
    {
        constraint_work_j = finite(
            multiplier.mul_add(0.5 * initial + 0.5 * final_rate, constraint_work_j),
            "release constraint work",
        )?;
    }
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
    Ok(PayloadRelease {
        model_after,
        state_after,
        payload_after,
        supplied_work_j,
        constraint_work_j,
        applied_impulse: impulse,
        impulse_on_payload_body,
        impulse_on_retained_link,
        balance,
    })
}

fn centered(velocity: Twist, pose: Se3, boost: Vec3) -> Result<Twist, PayloadError> {
    let local_boost = pose.rotation().inverse().rotate(boost)?;
    let centered = Twist::new(velocity.angular, velocity.linear - local_boost);
    validate_twist(centered, "payload_release.centered_velocity")?;
    Ok(centered)
}

fn add_joint_delta(
    rates: &mut [f64],
    delta: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PayloadError> {
    if rates.len() != delta.len() {
        return Err(PayloadError::Invalid("release joint response dimension"));
    }
    for (index, (rate, change)) in rates.iter_mut().zip(delta).enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        *rate = finite(*rate + change, "release joint velocity")?;
    }
    Ok(())
}

pub(super) fn admit_attachment(
    model: &ArticulatedModel,
    attachment: &PayloadAttachment,
) -> Result<(), PayloadError> {
    let loaded_link = model
        .links()
        .get(attachment.target_link)
        .ok_or(PayloadError::Invalid("release target link"))?;
    let retained_link = &attachment.retained_link;
    if loaded_link.name() != retained_link.name()
        || loaded_link.parent() != retained_link.parent()
        || loaded_link.parent_from_child_zero() != retained_link.parent_from_child_zero()
        || loaded_link.joint() != retained_link.joint()
        || loaded_link.inertia() != attachment.loaded_inertia
    {
        return Err(PayloadError::Invalid(
            "release attachment inventory changed",
        ));
    }
    Ok(())
}

//! Instantaneous response with the original prescribed base held fixed.
//!
//! The ABA backward pass condenses the joint impulse equations to
//! `I_base * delta_v_base + p_base = impulse_support`. With zero prescribed
//! base velocity jump, its existing bias is the support reaction. No large
//! surrogate base mass, force-duration conversion, or velocity projection is
//! needed. A moving prescribed support can exchange signed impulse work.

use super::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, Dd, GaError,
    ImpulseVelocityChange, PRESCRIBED_ABA_FIELDS, Se3, Twist, Vec3, Wrench, articulated_body_pass,
    ensure_finite_scalar, forward_kinematics, momentum_from_kinematics,
    momentum_increment_from_kinematics, propagate_aba_accelerations, require_len, validate_vec3,
    validate_wrench, wrench_from_array,
};

/// Fixed-configuration articulated response with no jump in the base twist.
///
/// The support impulse is exerted **on** the retained tree, including any
/// root-joint actuator reaction. Its opposite acts on the prescribed support.
/// A stationary base does no work; an externally driven base supplies the
/// explicitly reported signed work. This is an instantaneous inertia map,
/// not a contact law, a finite-duration integration, or a finite support model.
#[derive(Debug, Clone, PartialEq)]
pub struct PrescribedBaseImpulseResponse {
    /// Unchanged base pose, twist, and smooth acceleration declaration.
    pub base_after: BaseState,
    /// Actual updated generalized velocities, in compact joint order.
    pub joint_velocity_after: Vec<f64>,
    /// Physical momentum and kinetic energy before the event.
    pub momentum_before: ArticulatedMomentum,
    /// Physical momentum and kinetic energy at the actual returned endpoint.
    pub momentum_after: ArticulatedMomentum,
    /// Actual world linear increment, without subtracting large endpoint totals.
    pub linear_momentum_change_world: Vec3,
    /// Actual world-origin angular increment, including all retained links.
    pub angular_momentum_change_world: Vec3,
    /// Polarized kinetic-energy increment from actual link velocity changes, J.
    pub kinetic_energy_change: f64,
    /// Commanded external body impulse about the world origin, excluding support.
    /// Force entries are N s; torque entries are N m s.
    pub external_impulse_world: Wrench,
    /// Support reaction about the base origin, expressed in base coordinates.
    /// Force entries are N s; torque entries are N m s.
    pub support_impulse_body: Wrench,
    /// The same support reaction about the world origin, in world coordinates.
    pub support_impulse_world: Wrench,
    /// External body and relative joint impulses paired with midpoint rates, J.
    pub applied_impulse_work: f64,
    /// Support reaction paired with the unchanged prescribed base twist, J.
    pub support_work: f64,
    /// Total applied plus support impulse work, J.
    pub impulse_work: f64,
    /// Actual kinetic increment minus total impulse work, J.
    pub kinetic_work_residual: f64,
    /// Actual linear increment minus external and support impulses, N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual angular increment minus external and support impulses, N m s.
    pub angular_momentum_residual_world: Vec3,
}

/// Apply external spatial and joint impulses with a prescribed base boundary.
///
/// `body_impulses` are about each link origin in its own coordinates, in
/// `[N m s, N s]`; a point impulse J at local r is `Wrench::new(r.cross(J), J)`.
/// Joint impulses are N m s for revolute and N s for prismatic coordinates.
/// They are never compared with force/torque effort limits. Configuration and
/// actual endpoint velocity limits remain enforced, including at a scalar root
/// joint. The base pose and twist remain unchanged; its smooth acceleration does
/// not enter this instantaneous map. Gravity contributes no finite impulse.
///
/// The existing ABA recursions solve the constrained inverse mass in O(n)
/// work/storage without a dense generalized matrix or a floating-base solve.
/// The support reaction is independently obtained from the condensed dynamics;
/// actual endpoint momentum and kinetic work provide residuals, not repaired
/// state. A nonzero solved joint jump completely lost in endpoint rounding
/// refuses, including a changed joint rate erased by the unchanged parent
/// twist in native link kinematics. Constrained zero jumps retain their
/// legitimate support reactions.
/// All inputs are borrowed/copied; refusal publishes no partial update.
pub fn apply_prescribed_base_impulses(
    model: &ArticulatedModel,
    base: BaseState,
    joint_position: &[f64],
    joint_velocity: &[f64],
    joint_impulses: &[f64],
    body_impulses: &[Wrench],
) -> Result<PrescribedBaseImpulseResponse, ArticulatedError> {
    let before = forward_kinematics(model, base, joint_position, joint_velocity)?;
    let delta =
        prescribed_base_impulse_delta(model, joint_position, joint_impulses, body_impulses)?;
    let joint_velocity_after: Vec<_> = joint_velocity
        .iter()
        .zip(&delta.velocity.joint)
        .map(|(rate, jump)| {
            let endpoint = rate + jump;
            // Exact no-change detection, not an adjustable accuracy tolerance.
            #[allow(clippy::float_cmp)]
            let swallowed = *jump != 0.0 && endpoint == *rate;
            if swallowed {
                Err(ArticulatedError::Geometry(GaError::InvalidRepresentation {
                    context: "prescribed_impulse.joint_velocity_jump",
                    defect: jump.abs(),
                    tolerance: 0.0,
                }))
            } else {
                Ok(endpoint)
            }
        })
        .collect::<Result<_, _>>()?;
    let after = forward_kinematics(model, base, joint_position, &joint_velocity_after)?;
    for (link_index, link) in model.links.iter().enumerate() {
        if let Some(dof) = model.dof_index[link_index] {
            let unchanged_parent = link
                .parent
                .is_none_or(|parent| before.body_twist[parent] == after.body_twist[parent]);
            // With an unchanged parent and fixed configuration, a nonzero
            // scalar-joint jump MUST change this link twist. The joint rate
            // can be representable while its addition to a large base/parent
            // twist loses the whole physical response. Refuse that exact
            // no-change case without confusing deliberate cancellation by a
            // moving parent with an unrepresentable event.
            if delta.velocity.joint[dof] != 0.0
                && unchanged_parent
                && before.body_twist[link_index] == after.body_twist[link_index]
            {
                return Err(ArticulatedError::Geometry(GaError::InvalidRepresentation {
                    context: "prescribed_impulse.link_velocity_jump",
                    defect: delta.velocity.joint[dof].abs(),
                    tolerance: 0.0,
                }));
            }
        }
    }
    let momentum_before = momentum_from_kinematics(model, &before)?;
    let momentum_after = momentum_from_kinematics(model, &after)?;
    let change = momentum_increment_from_kinematics(model, &before, &after)?;

    let mut external = [Dd::ZERO; 6];
    let mut work = Dd::ZERO;
    for (index, impulse) in body_impulses.iter().enumerate() {
        add_world_impulse(&mut external, *impulse, before.world_from_link[index])?;
        for ((old, new), applied) in before.body_twist[index]
            .to_array()
            .into_iter()
            .zip(after.body_twist[index].to_array())
            .zip(impulse.to_array())
        {
            work =
                work + (Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new)) * Dd::from_f64(applied);
        }
    }
    for (index, impulse) in joint_impulses.iter().enumerate() {
        work = work
            + Dd::from_f64(*impulse)
                * (Dd::from_f64(0.5 * joint_velocity[index])
                    + Dd::from_f64(0.5 * joint_velocity_after[index]));
    }
    let external_impulse_world = wrench_from_array(external.map(Dd::to_f64));
    validate_wrench(external_impulse_world, "prescribed_impulse.external_world")?;
    let support_impulse_body = delta.support_body;
    let mut support = [Dd::ZERO; 6];
    add_world_impulse(&mut support, support_impulse_body, base.world_from_base)?;
    let support_impulse_world = wrench_from_array(support.map(Dd::to_f64));
    validate_wrench(support_impulse_world, "prescribed_impulse.support_world")?;
    let support_work = pairing(support_impulse_body, base.twist_body);
    let total_work = work + support_work;
    let residual = Dd::from_f64(change.kinetic_energy) - total_work;
    for value in [work, support_work, total_work, residual] {
        ensure_finite_scalar(value.to_f64(), "prescribed_impulse.work", 0)?;
    }
    let actual = Wrench::new(change.angular_momentum_world, change.linear_momentum_world);
    let actual = actual.to_array();
    let momentum_residual = wrench_from_array(std::array::from_fn(|index| {
        (Dd::from_f64(actual[index]) - external[index] - support[index]).to_f64()
    }));
    validate_vec3(
        momentum_residual.force,
        "prescribed_impulse.linear_residual",
    )?;
    validate_vec3(
        momentum_residual.torque,
        "prescribed_impulse.angular_residual",
    )?;
    Ok(PrescribedBaseImpulseResponse {
        base_after: base,
        joint_velocity_after,
        momentum_before,
        momentum_after,
        linear_momentum_change_world: change.linear_momentum_world,
        angular_momentum_change_world: change.angular_momentum_world,
        kinetic_energy_change: change.kinetic_energy,
        external_impulse_world,
        support_impulse_body,
        support_impulse_world,
        applied_impulse_work: work.to_f64(),
        support_work: support_work.to_f64(),
        impulse_work: total_work.to_f64(),
        kinetic_work_residual: residual.to_f64(),
        linear_momentum_residual_world: momentum_residual.force,
        angular_momentum_residual_world: momentum_residual.torque,
    })
}

pub(crate) struct PrescribedImpulseDelta {
    pub(crate) velocity: ImpulseVelocityChange,
    pub(crate) support_body: Wrench,
}

/// Inverse-mass probes do not apply an event or enforce endpoint speed limits.
pub(crate) fn prescribed_base_impulse_delta(
    model: &ArticulatedModel,
    joint_position: &[f64],
    joint_impulses: &[f64],
    body_impulses: &[Wrench],
) -> Result<PrescribedImpulseDelta, ArticulatedError> {
    model.validate_configuration(joint_position)?;
    model.validate_generalized("joint_impulses", joint_impulses)?;
    require_len("body_impulses", body_impulses.len(), model.link_count())?;
    for impulse in body_impulses {
        validate_wrench(*impulse, "body_impulses")?;
    }
    let zero = vec![0.0; model.dof_count()];
    let pass = articulated_body_pass(
        model,
        Twist::zero(),
        joint_position,
        &zero,
        joint_impulses,
        body_impulses,
        None,
        true,
        PRESCRIBED_ABA_FIELDS,
    )?;
    let (joint, body) = propagate_aba_accelerations(
        model,
        &pass,
        Twist::zero(),
        "prescribed_impulse.joint_velocity_change",
        "prescribed_impulse.body_velocity_change",
    )?;
    Ok(PrescribedImpulseDelta {
        velocity: ImpulseVelocityChange {
            base: Twist::zero(),
            joint,
            body,
        },
        support_body: pass.base_articulated_bias,
    })
}

fn pairing(impulse: Wrench, velocity: Twist) -> Dd {
    impulse
        .to_array()
        .into_iter()
        .zip(velocity.to_array())
        .fold(Dd::ZERO, |sum, (p, v)| {
            sum + Dd::from_f64(p) * Dd::from_f64(v)
        })
}

// Preserve local torque until the complete signed world-origin reduction.
// Combining it with a large orbital term inside each transformed wrench can
// erase a real couple when opposite forces cancel elsewhere in the tree.
fn add_world_impulse(
    sum: &mut [Dd; 6],
    impulse: Wrench,
    world_from_body: Se3,
) -> Result<(), ArticulatedError> {
    let force = world_from_body.rotation().rotate(impulse.force)?;
    let torque = world_from_body.rotation().rotate(impulse.torque)?;
    let orbital = world_from_body.translation().cross(force);
    for (total, value) in sum.iter_mut().zip(Wrench::new(torque, force).to_array()) {
        *total = *total + Dd::from_f64(value);
    }
    for (total, value) in sum[..3].iter_mut().zip([orbital.x, orbital.y, orbital.z]) {
        *total = *total + Dd::from_f64(value);
    }
    validate_wrench(
        wrench_from_array(sum.map(Dd::to_f64)),
        "prescribed_impulse.world_resultant",
    )
}

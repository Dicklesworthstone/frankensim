//! Continuous-velocity release of a supported tree's captured rigid inventory.
//!
//! The exact retained inertia comes from the original attachment token. Native
//! FK supplies the payload's current pose and rigid velocity field. Removing a
//! weld without an impulse changes no velocity, so the prescribed support has
//! no impulse or work. Independent decomposition residuals remain visible;
//! they are never repaired by a velocity jump, fictitious recoil or energy loss.

use super::super::release::admit_attachment;
use super::super::{
    PayloadAttachment, PayloadError, PayloadState, account_boost, energy_in_frame, finite,
    payload_at_link, poll, replace_link_inertia, twist_in_frame,
};
use super::{
    SupportedPayloadBalance, SupportedPayloadControls, admit, audit_momentum, negate,
};
use crate::articulated::motion::SpatialSum;
use crate::articulated::motion::prescribed::PrescribedBaseMotionState;
use crate::articulated::{
    ArticulatedModel, BaseState, forward_kinematics, inverse_dynamics, momentum_from_kinematics,
    validate_wrench,
};
use fs_ga::{Vec3, Wrench};
use fs_math::dd::Dd;

/// Independent absolute allowances for a zero-impulse inventory decomposition.
#[derive(Clone, Copy, Debug)]
pub struct SupportedPayloadReleaseControls<'a> {
    /// Positive world angular momentum decomposition allowance per axis, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Positive world linear momentum decomposition allowance per axis, N s.
    pub linear_momentum_tolerance_ns: f64,
    /// Positive allowances for native joint momentum decomposition, N m s or N s.
    pub joint_momentum_tolerances: &'a [f64],
    /// Nonnegative common-frame energy decomposition defect allowance, J.
    pub energy_tolerance_j: f64,
}

/// Restored supported inventory and the actual newly independent rigid body.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportedPayloadRelease {
    /// Same topology and joints; exact retained target inertia from capture.
    pub model_after: ArticulatedModel,
    /// Bit-preserved time, base pose/twist/acceleration, coordinates and rates.
    pub state_after: PrescribedBaseMotionState,
    /// Original payload inertia with its current pose and inherited body twist.
    pub payload_after: PayloadState,
    /// Actual loaded-to-separated momentum/energy decomposition diagnostics.
    /// Support impulses and work are zero; any nonzero residual is rounding.
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
/// Native FK/RNEA and independent body sums audit the restored tree plus free
/// payload against the loaded tree. Physical impulse, support work and energy
/// loss are all zero, with finite decomposition defects reported separately.
///
/// This only removes the ideal rigid attachment. Collision discovery, ongoing
/// contact, fracture energy, finite-duration separation and ejection actuators
/// remain separate models. Refusal or cancellation publishes no partial state.
pub fn release_payload(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    attachment: &PayloadAttachment,
    controls: SupportedPayloadReleaseControls<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SupportedPayloadRelease, PayloadError> {
    poll(&mut cancelled)?;
    admit_attachment(model, attachment)?;
    // Shared admission/audit fields; this release has no loss/exchange budget.
    let shared_controls = SupportedPayloadControls {
        angular_momentum_tolerance_nms: controls.angular_momentum_tolerance_nms,
        linear_momentum_tolerance_ns: controls.linear_momentum_tolerance_ns,
        joint_momentum_tolerances: controls.joint_momentum_tolerances,
        energy_tolerance_j: controls.energy_tolerance_j,
        maximum_removed_relative_energy_j: 0.0,
    };
    let target = attachment.target_link;
    admit(model, state, target, shared_controls)?;
    let kin = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    let payload_after = payload_at_link(
        &kin,
        target,
        attachment.link_from_payload,
        attachment.payload_inertia,
    )?;
    let model_after = replace_link_inertia(model, target, attachment.retained_link.inertia())?;
    // FK is inertia-independent, and neither geometry nor any state bit changes.
    let momentum_before = momentum_from_kinematics(model, &kin)?;
    poll(&mut cancelled)?;
    let mut momentum_after = momentum_from_kinematics(&model_after, &kin)?;
    let payload_momentum = payload_after.inertia.momentum(payload_after.twist_body)?;
    let payload_world = payload_momentum.transform_by(&payload_after.world_from_body);
    momentum_after.linear_momentum_world =
        momentum_after.linear_momentum_world + payload_world.force;
    momentum_after.angular_momentum_world =
        momentum_after.angular_momentum_world + payload_world.torque;
    momentum_after.kinetic_energy = finite(
        momentum_after.kinetic_energy
            + payload_after.inertia.kinetic_energy(payload_after.twist_body)?,
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

    let target_pose = kin.world_from_link[target];
    let target_twist = kin.body_twist[target];
    let loaded = model.links()[target].inertia();
    let retained = model_after.links()[target].inertia();
    let loaded_momentum = loaded.momentum(target_twist)?;
    let retained_momentum = retained.momentum(target_twist)?;
    let mut change = SpatialSum::default();
    for (value, pose) in [
        (negate(loaded_momentum), target_pose),
        (retained_momentum, target_pose),
        (payload_momentum, payload_after.world_from_body),
    ] {
        change.add_world_wrench(
            pose.rotation().rotate(value.force)?,
            pose.rotation().rotate(value.torque)?,
            &[pose.translation()],
        );
    }
    let physical = change.checked_wrench("supported_payload_release.decomposition")?;
    // Only the target inertia changed. Project its independently reconstructed
    // difference through the original joint screws; no generalized mass solve.
    let incoming = payload_momentum
        .transform_by(&attachment.link_from_payload)
        .to_array();
    let old = loaded_momentum.to_array();
    let new = retained_momentum.to_array();
    let local: [f64; 6] = core::array::from_fn(|i| {
        (Dd::from_f64(new[i]) + Dd::from_f64(incoming[i]) - Dd::from_f64(old[i]))
            .to_f64()
    });
    let local = Wrench::new(
        Vec3::new(local[0], local[1], local[2]),
        Vec3::new(local[3], local[4], local[5]),
    );
    validate_wrench(local, "supported_payload_release.local_decomposition")?;
    let mut projected = vec![Wrench::default(); model.link_count()];
    projected[target] = negate(local);
    let zero = vec![0.0; model.dof_count()];
    poll(&mut cancelled)?;
    let joint_momentum_residual = inverse_dynamics(
        &model_after,
        BaseState::stationary(state.base.world_from_base),
        &state.joint_position,
        &zero,
        &zero,
        Vec3::default(),
        &projected,
    )?
    .generalized_force;
    poll(&mut cancelled)?;
    let maximum_normalized_momentum_residual =
        audit_momentum(physical, &joint_momentum_residual, shared_controls)?;

    let boost = account_boost(model, &kin)?;
    let centered_target = twist_in_frame(target_twist, target_pose, boost)?;
    let centered_payload = twist_in_frame(
        payload_after.twist_body,
        payload_after.world_from_body,
        boost,
    )?;
    let decomposition = |target_velocity, payload_velocity| -> Result<f64, PayloadError> {
        finite(
            (Dd::from_f64(retained.kinetic_energy(target_velocity)?)
                + Dd::from_f64(payload_after.inertia.kinetic_energy(payload_velocity)?)
                - Dd::from_f64(loaded.kinetic_energy(target_velocity)?))
            .to_f64(),
            "supported_payload_release.energy_decomposition",
        )
    };
    let kinetic_energy_change_j = decomposition(target_twist, payload_after.twist_body)?;
    let energy_residual_j = decomposition(centered_target, centered_payload)?;
    if energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(PayloadError::EnergyResidual {
            residual_j: energy_residual_j,
            maximum_j: controls.energy_tolerance_j,
        });
    }
    let kinetic_before_j = energy_in_frame(model, &kin, None, boost, &mut cancelled)?;
    let kinetic_after_j = energy_in_frame(
        &model_after,
        &kin,
        Some(&payload_after),
        boost,
        &mut cancelled,
    )?;
    let balance = SupportedPayloadBalance {
        momentum_before,
        momentum_after,
        linear_momentum_change_world: physical.force,
        angular_momentum_change_world: physical.torque,
        support_impulse_body: Wrench::default(),
        support_impulse_world: Wrench::default(),
        linear_momentum_residual_world: physical.force,
        angular_momentum_residual_world: physical.torque,
        joint_momentum_residual,
        maximum_normalized_momentum_residual,
        kinetic_energy_change_j,
        support_work_j: 0.0,
        accounting_boost_world: boost,
        kinetic_before_j,
        kinetic_after_j,
        accounting_kinetic_energy_change_j: energy_residual_j,
        accounting_support_work_j: 0.0,
        energy_residual_j,
    };
    poll(&mut cancelled)?;
    Ok(SupportedPayloadRelease {
        model_after,
        state_after: state.clone(),
        payload_after,
        balance,
    })
}

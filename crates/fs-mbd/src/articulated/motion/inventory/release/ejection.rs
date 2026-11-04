//! Final-relative-twist ejection through the original articulated inventory.
//!
//! A commanded payload twist is relative to its carrier's FINAL velocity field
//! at the payload origin. Loaded-tree momentum determines the shared recoil;
//! the original retained-tree release owns all physical endpoints and work.

use super::super::{
    PayloadAttachment, PayloadControls, PayloadError, account_boost, bilateral, finite, poll,
    validate_inputs,
};
use super::{
    PayloadBatchRelease, PayloadReleaseInput, SeparationImpulse, admit_attachment_link,
    separation_impulses,
};
use crate::articulated::motion::FreeFloatingMotionState;
use crate::articulated::motion::SpatialSum;
use crate::articulated::motion::reconfigure::constraints::{
    JointConstraint, PreparedJointConstraints, apply_rows,
};
use crate::articulated::motion::reconfigure::state_kinematics;
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralImpulseResponse, BilateralLoop, joint,
};
use crate::articulated::{
    ArticulatedModel, ImpulseVelocityChange, free_floating_impulse_delta, validate_twist,
};
use fs_ga::{Twist, Vec3, Wrench};
use fs_math::dd::Dd;

/// Maximum commanded payloads, admitted before any per-event staging.
pub const MAX_PAYLOAD_EJECTION_TARGETS: usize = 4096;

/// One exact captured payload and its final motion relative to its carrier.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PayloadEjectionInput<'a> {
    /// Existing inventory identity; repeated targets consume handles in reverse
    /// capture order, as required by the original native batch release.
    pub attachment: &'a PayloadAttachment,
    /// Angular velocity difference [rad/s] and ORIGIN linear velocity difference
    /// [m/s], both in the original payload axes. The carrier velocity is evaluated
    /// at this SAME payload origin in its FINAL state; this is not a difference
    /// between velocities at two different frame origins or two different COMs.
    pub relative_twist_payload: Twist,
}

/// Original inventory/work limits and actual endpoint velocity tolerances.
#[derive(Debug, Clone, Copy)]
pub struct PayloadEjectionControls<'a> {
    /// All native momentum, common-frame energy and separate source/reaction work caps.
    pub payload: PayloadControls<'a>,
    /// Finite nonnegative maximum absolute component of linear target error [m/s].
    pub relative_velocity_tolerance_m_s: f64,
    /// Finite nonnegative maximum absolute component of spin target error [rad/s].
    pub relative_angular_velocity_tolerance_rad_s: f64,
}

/// One accepted native release realizing every commanded relative twist.
#[derive(Debug)]
pub struct PayloadEjection {
    /// Original authoritative retained/payload states, derived impulses and work.
    pub release: PayloadBatchRelease,
    /// Largest component of actual returned linear target error [m/s].
    pub maximum_relative_velocity_residual_m_s: f64,
    /// Largest component of actual returned spin target error [rad/s].
    pub maximum_relative_angular_velocity_residual_rad_s: f64,
}

/// Native refusal or failure of the actual returned velocity boundary condition.
#[derive(Debug, Clone, PartialEq)]
pub enum PayloadEjectionError {
    /// Original inventory, native physics, work, budget or cancellation refusal.
    Payload(PayloadError),
    /// An actual payload-origin velocity does not realize its command.
    RelativeVelocity {
        /// Index in the authored input order.
        payload_index: usize,
        /// Actual minus requested relative velocity, in payload axes [m/s].
        residual_payload_m_s: Vec3,
        /// Declared maximum absolute component [m/s].
        tolerance_m_s: f64,
    },
    /// An actual relative angular velocity does not realize its command.
    RelativeAngularVelocity {
        /// Index in the authored input order.
        payload_index: usize,
        /// Actual minus requested relative angular velocity, payload axes [rad/s].
        residual_payload_rad_s: Vec3,
        /// Declared maximum absolute component [rad/s].
        tolerance_rad_s: f64,
    },
}
impl From<PayloadError> for PayloadEjectionError {
    fn from(error: PayloadError) -> Self {
        Self::Payload(error)
    }
}
impl core::fmt::Display for PayloadEjectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "relative payload ejection: {self:?}")
    }
}
impl std::error::Error for PayloadEjectionError {}

/// Eject a nonempty collection at commanded final relative linear/spin rates.
///
/// Let S_i map generalized carrier velocity into payload-origin coordinates and
/// let I_i be the original full payload spatial inertia. Desired relative twist
/// d_i gives payload velocity S_i*q_after+d_i. Conservation therefore uses the
/// ORIGINAL LOADED inertia:
///
/// M_loaded*(q_after-q_before) = -sum(S_i^T*I_i*d_i) + J_constraint.
///
/// One loaded native response determines common recoil. Each required payload
/// impulse is I_i*(S_i*delta_q+d_i), evaluated from velocity CHANGES so a large
/// common translation never enters this weak impulse calculation. Force acts
/// at the original payload origin and its angular impulse is represented as a
/// free couple. This wrench decomposition does not identify real actuators.
///
/// Existing joint reactions use the loaded inertia and retain correction of
/// admitted nonzero initial tangent defects. The original batch release then
/// applies the derived impulses under the actual RETAINED inertia, reconstructs
/// every detached full-inertia state and owns all work and momentum gates.
/// Auxiliary loaded-response work is neither returned nor charged as physical
/// separation work. At an exact endpoint, each pair's actual work is
/// 0.5*d_i.pairing(payload_impulse_i), including shared recoil cross terms.
///
/// All-zero commands and exactly cancelling loaded-side sources keep the
/// original passive state semantics; there is no unsolicited tangent repair.
/// A zero command in an active batch may require an impulse to follow shared
/// carrier recoil. Input count is bounded before staging. Outside the original
/// row response, work/storage is O(links + payloads), with no dense 6N mobility.
///
/// Actual returned relative rates must satisfy their declared tolerances before
/// publication. This admits representable endpoints, not ideal rates lost under
/// a huge common velocity. Native source/constraint work limits remain separate.
/// Contact selection, actuator geometry, finite-time force laws, fracture,
/// supported-base ejection and thermal/source-energy evolution remain external.
///
/// # Errors
/// Empty/oversized or stale inventory, malformed target/tolerance, native
/// conditioning/constraint/endpoint refusal, physical work or balance refusal,
/// unrealized final relative rates, or cancellation. No input is mutated.
pub fn eject_payloads(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadEjectionInput<'_>],
    constraints: &[JointConstraint],
    controls: PayloadEjectionControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PayloadEjection, PayloadEjectionError> {
    eject_evaluated(model, state, inputs, constraints, controls, None, cancelled)
        .map(|(ejection, _)| ejection)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::articulated::motion::inventory) fn eject_evaluated(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadEjectionInput<'_>],
    constraints: &[JointConstraint],
    controls: PayloadEjectionControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PayloadEjection, Option<BilateralImpulseResponse>), PayloadEjectionError> {
    poll(&mut cancelled)?;
    if [
        controls.relative_velocity_tolerance_m_s,
        controls.relative_angular_velocity_tolerance_rad_s,
    ]
    .into_iter()
    .any(|x| !x.is_finite() || x < 0.0)
    {
        return Err(PayloadError::Invalid(
            "finite nonnegative ejection velocity tolerances required",
        )
        .into());
    }
    let prescribed = derive_impulses(
        model,
        state,
        inputs,
        constraints,
        controls.payload,
        bilateral_rows,
        &mut cancelled,
    )?;
    let (release, reaction) = super::release_payloads_evaluated(
        model,
        state,
        &prescribed,
        constraints,
        controls.payload,
        bilateral_rows,
        &mut cancelled,
    )?;
    let after =
        state_kinematics(&release.model_after, &release.state_after).map_err(PayloadError::from)?;
    let mut maximum_relative_velocity_residual_m_s: f64 = 0.0;
    let mut maximum_relative_angular_velocity_residual_rad_s: f64 = 0.0;
    for (index, (input, output)) in inputs.iter().zip(&release.payloads_after).enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        let payload_from_link = input
            .attachment
            .link_from_payload
            .inverse()
            .map_err(PayloadError::from)?;
        let carried =
            after.body_twist[input.attachment.target_link].transform_by(&payload_from_link);
        let residual = residual_twist(
            output.payload_after.twist_body,
            carried,
            input.relative_twist_payload,
        )?;
        let linear = maximum_component(residual.linear);
        let angular = maximum_component(residual.angular);
        if linear > controls.relative_velocity_tolerance_m_s {
            return Err(PayloadEjectionError::RelativeVelocity {
                payload_index: index,
                residual_payload_m_s: residual.linear,
                tolerance_m_s: controls.relative_velocity_tolerance_m_s,
            });
        }
        if angular > controls.relative_angular_velocity_tolerance_rad_s {
            return Err(PayloadEjectionError::RelativeAngularVelocity {
                payload_index: index,
                residual_payload_rad_s: residual.angular,
                tolerance_rad_s: controls.relative_angular_velocity_tolerance_rad_s,
            });
        }
        maximum_relative_velocity_residual_m_s = maximum_relative_velocity_residual_m_s.max(linear);
        maximum_relative_angular_velocity_residual_rad_s =
            maximum_relative_angular_velocity_residual_rad_s.max(angular);
    }
    poll(&mut cancelled)?;
    Ok((
        PayloadEjection {
            release,
            maximum_relative_velocity_residual_m_s,
            maximum_relative_angular_velocity_residual_rad_s,
        },
        reaction,
    ))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn derive_impulses<'a>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadEjectionInput<'a>],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<PayloadReleaseInput<'a>>, PayloadError> {
    poll(cancelled)?;
    let first = inputs
        .first()
        .ok_or(PayloadError::Invalid("empty payload ejection"))?;
    if inputs.len() > MAX_PAYLOAD_EJECTION_TARGETS {
        return Err(PayloadError::Invalid(
            "payload ejection exceeds simultaneous body limit",
        ));
    }
    let before = validate_inputs(
        model,
        state,
        first.attachment.target_link,
        constraints,
        controls,
        cancelled,
    )?;
    // References validate inventory prefixes without cloning names/inertias or
    // constructing any intermediate physical model.
    let mut prefixes: Vec<_> = model.links().iter().collect();
    let mut sums = Vec::with_capacity(model.link_count());
    for index in 0..model.link_count() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        sums.push(SpatialSum::default());
    }
    let mut passive = true;
    for (index, input) in inputs.iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let attachment = input.attachment;
        let target = attachment.target_link;
        let prefix = prefixes
            .get_mut(target)
            .ok_or(PayloadError::Invalid("ejection target link"))?;
        admit_attachment_link(*prefix, attachment)?;
        *prefix = &attachment.retained_link;
        validate_twist(
            input.relative_twist_payload,
            "payload_ejection.relative_twist",
        )?;
        passive &= input.relative_twist_payload == Twist::zero();
        if let Some(prescribed) = prescribe(attachment, input.relative_twist_payload)? {
            // Use the original separation decomposition for this X*I*d action.
            // If sources cancel exactly, this is the SAME rounded link sum
            // later seen by native release; no tiny source activates projection.
            let (retained, _) = separation_impulses(prescribed, attachment.link_from_payload)?;
            sums[target].add_world_wrench(retained.force, retained.torque, &[]);
        }
    }
    let mut source = Vec::with_capacity(sums.len());
    let mut nonzero_source = false;
    for (index, sum) in sums.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let impulse = sum.checked_wrench("relative payload loaded source")?;
        nonzero_source |= impulse != Wrench::default();
        source.push(impulse);
    }
    // A passive command must retain even finite admitted initial row defects.
    // Exactly cancelling nonzero targets also need no loaded-carrier response.
    let delta = if passive || !nonzero_source {
        None
    } else {
        Some(loaded_response(
            model,
            state,
            &before,
            &source,
            constraints,
            controls,
            bilateral_rows,
            cancelled,
        )?)
    };
    let mut prescribed = Vec::with_capacity(inputs.len());
    for (index, input) in inputs.iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let attachment = input.attachment;
        let change = if let Some(delta) = &delta {
            delta.body[attachment.target_link]
                .transform_by(&attachment.link_from_payload.inverse()?)
                .plus(input.relative_twist_payload)
        } else {
            input.relative_twist_payload
        };
        prescribed.push(PayloadReleaseInput {
            attachment,
            impulse: prescribe(attachment, change)?,
        });
    }
    poll(cancelled)?;
    Ok(prescribed)
}

fn prescribe(
    attachment: &PayloadAttachment,
    change: Twist,
) -> Result<Option<SeparationImpulse>, PayloadError> {
    validate_twist(change, "payload_ejection.payload_change")?;
    let impulse = attachment.payload_inertia.momentum(change)?;
    let rotation = attachment.link_from_payload.rotation();
    let separation = SeparationImpulse {
        point_link_m: attachment.link_from_payload.translation(),
        impulse_on_payload_link_ns: rotation.rotate(impulse.force)?,
        couple_on_payload_link_nms: rotation.rotate(impulse.torque)?,
    };
    Ok(
        if separation.impulse_on_payload_link_ns == Vec3::default()
            && separation.couple_on_payload_link_nms == Vec3::default()
        {
            None
        } else {
            Some(separation)
        },
    )
}

/// Auxiliary native loaded response; its work is not a physical energy debit.
#[allow(clippy::too_many_arguments)]
fn loaded_response(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    before: &crate::articulated::Kinematics,
    source: &[Wrench],
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ImpulseVelocityChange, PayloadError> {
    let zero = vec![0.0; model.dof_count()];
    if let Some((loops, config)) = bilateral_rows {
        // Remove only a uniform translation. Original angular/joint rates and
        // finite admitted tangent defects remain in the native mixed response.
        let boost = account_boost(model, before)?;
        let local_boost = state
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)?;
        let mut centered = state.clone();
        centered.base.twist_body.linear = centered.base.twist_body.linear - local_boost;
        let response = joint::apply_impulses(
            model,
            &centered,
            &zero,
            source,
            loops,
            constraints,
            config,
            &mut *cancelled,
        )?;
        bilateral::admit_condition(response.reciprocal_condition, controls)?;
        let mut combined = Vec::with_capacity(source.len());
        for (index, (source, reaction)) in source
            .iter()
            .zip(&response.constraint_body_impulses)
            .enumerate()
        {
            if index % 64 == 0 {
                poll(cancelled)?;
            }
            // Match the original native source/reaction addition exactly.
            combined.push(Wrench::new(
                source.torque + reaction.torque,
                source.force + reaction.force,
            ));
        }
        poll(cancelled)?;
        let delta = free_floating_impulse_delta(
            model,
            &state.joint_position,
            &response.constraint_joint_impulses,
            &combined,
        )?;
        poll(cancelled)?;
        return Ok(delta);
    }
    let source_delta = free_floating_impulse_delta(model, &state.joint_position, &zero, source)?;
    poll(cancelled)?;
    if constraints.is_empty() {
        return Ok(source_delta);
    }
    let mut prepared = PreparedJointConstraints::new(
        model,
        &state.joint_position,
        constraints,
        controls.minimum_reciprocal_condition,
        cancelled,
    )?;
    let mut candidate_rates = Vec::with_capacity(state.joint_velocity.len());
    for (index, (&old, &change)) in state
        .joint_velocity
        .iter()
        .zip(&source_delta.joint)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        candidate_rates.push(finite(old + change, "ejection auxiliary joint rate")?);
    }
    let rhs = apply_rows(constraints, &candidate_rates, cancelled)?
        .into_iter()
        .map(|rate| -rate)
        .collect::<Vec<_>>();
    let reaction = prepared.solve_response(&rhs, cancelled)?;
    // Reuse the native full-body propagation without subtracting endpoint
    // velocities or applying endpoint speed limits to an inverse-mass column.
    let delta =
        free_floating_impulse_delta(model, &state.joint_position, &reaction.joint_effort, source)?;
    poll(cancelled)?;
    Ok(delta)
}

fn residual_twist(actual: Twist, carried: Twist, target: Twist) -> Result<Twist, PayloadError> {
    let mut values = [0.0; 6];
    for (value, ((actual, carried), target)) in values.iter_mut().zip(
        actual
            .to_array()
            .into_iter()
            .zip(carried.to_array())
            .zip(target.to_array()),
    ) {
        *value = finite(
            (Dd::from_f64(actual) - Dd::from_f64(carried) - Dd::from_f64(target)).to_f64(),
            "ejection actual relative velocity residual",
        )?;
    }
    Ok(Twist::new(
        Vec3::new(values[0], values[1], values[2]),
        Vec3::new(values[3], values[4], values[5]),
    ))
}

fn maximum_component(vector: Vec3) -> f64 {
    vector.x.abs().max(vector.y.abs()).max(vector.z.abs())
}

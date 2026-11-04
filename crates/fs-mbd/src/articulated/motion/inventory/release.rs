//! Release a retained payload with simultaneous articulated constraint recoil.
//!
//! Geometry is unchanged at the event. The detached body inherits the actual
//! link velocity field; an optional internal impulse pair then changes both
//! surviving systems through their native inverse inertias. No contact,
//! fracture, collision geometry, or destination for mechanical work is inferred.

/// Simultaneous release with shared retained-mechanism recoil.
#[path = "release/batch.rs"]
mod batch;
pub(super) use batch::release_evaluated as release_payloads_evaluated;
pub use batch::{PayloadBatchRelease, PayloadReleaseInput, ReleasedPayload, release_payloads};

/// Prescribed final payload-origin velocity and spin relative to retained links.
#[path = "release/ejection.rs"]
pub mod ejection;

use super::{
    PayloadAttachment, PayloadBalance, PayloadControls, PayloadError, PayloadState, finite, poll,
};
use crate::articulated::motion::FreeFloatingMotionState;
use crate::articulated::motion::reconfigure::constraints::JointConstraint;
use crate::articulated::motion::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralImpulseResponse, BilateralLoop,
};
use crate::articulated::{ArticulatedModel, Link, validate_twist, validate_vec3, validate_wrench};
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
    /// Accepted retained-tree motion after separation and ideal-constraint reactions.
    pub state_after: FreeFloatingMotionState,
    /// Independent released payload, with its original body frame and inertia.
    pub payload_after: PayloadState,
    /// Actual midpoint work of the complete internal separation impulse pair [J].
    /// The pair is evaluated together, including its quadratic recoil work.
    pub supplied_work_j: f64,
    /// Actual finite midpoint work of the simultaneous ideal-constraint reactions [J].
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
    cancelled: impl FnMut() -> bool,
) -> Result<PayloadRelease, PayloadError> {
    release_evaluated(
        model,
        state,
        attachment,
        impulse,
        constraints,
        controls,
        None,
        cancelled,
    )
    .map(|(release, _)| release)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn release_evaluated(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    attachment: &PayloadAttachment,
    impulse: Option<SeparationImpulse>,
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    bilateral_rows: Option<(&[BilateralLoop], BilateralImpulseConfig)>,
    cancelled: impl FnMut() -> bool,
) -> Result<(PayloadRelease, Option<BilateralImpulseResponse>), PayloadError> {
    let (release, reaction) = release_payloads_evaluated(
        model,
        state,
        &[PayloadReleaseInput {
            attachment,
            impulse,
        }],
        constraints,
        controls,
        bilateral_rows,
        cancelled,
    )?;
    let PayloadBatchRelease {
        model_after,
        state_after,
        mut payloads_after,
        supplied_work_j,
        constraint_work_j,
        balance,
    } = release;
    // One shared inventory/physics owner; unwrap only the singleton containers.
    let payload = payloads_after.remove(0);
    Ok((
        PayloadRelease {
            model_after,
            state_after,
            payload_after: payload.payload_after,
            supplied_work_j,
            constraint_work_j,
            applied_impulse: payload.applied_impulse,
            impulse_on_payload_body: payload.impulse_on_payload_body,
            impulse_on_retained_link: payload.impulse_on_retained_link,
            balance,
        },
        reaction,
    ))
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
    admit_attachment_link(loaded_link, attachment)
}

/// Check one exact inventory prefix without constructing an intermediate model.
fn admit_attachment_link(
    loaded_link: &Link,
    attachment: &PayloadAttachment,
) -> Result<(), PayloadError> {
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

/// Original point/couple admission, shared by free and prescribed releases.
pub(super) fn admit_separation_impulse(impulse: SeparationImpulse) -> Result<(), PayloadError> {
    validate_vec3(impulse.point_link_m, "payload_release.point")?;
    validate_vec3(
        impulse.impulse_on_payload_link_ns,
        "payload_release.force_impulse",
    )?;
    validate_vec3(
        impulse.couple_on_payload_link_nms,
        "payload_release.free_couple",
    )?;
    Ok(())
}

/// Equal/opposite impulses about the original link and payload frame origins.
/// Both actions retain one shared physical application point.
pub(super) fn separation_impulses(
    impulse: SeparationImpulse,
    link_from_payload: Se3,
) -> Result<(Wrench, Wrench), PayloadError> {
    let at_link = Wrench::new(
        impulse
            .point_link_m
            .cross(impulse.impulse_on_payload_link_ns)
            + impulse.couple_on_payload_link_nms,
        impulse.impulse_on_payload_link_ns,
    );
    validate_wrench(at_link, "payload_release.link_impulse")?;
    let on_retained = Wrench::new(at_link.torque.scale(-1.0), at_link.force.scale(-1.0));
    let on_payload = at_link.transform_by(&link_from_payload.inverse()?);
    validate_wrench(on_payload, "payload_release.payload_impulse")?;
    Ok((on_retained, on_payload))
}

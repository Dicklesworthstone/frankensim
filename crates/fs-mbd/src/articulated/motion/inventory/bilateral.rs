//! Payload capture and release retaining the original mixed rigid-constraint response.
//!
//! Loaded inertia, payload mismatch and all existing material/joint rows react
//! together. The original inventory owner retains attachment identity, physical
//! momentum audits, common-frame energy and nonnegative capture quadratics.
//! No constraint engagement, pose projection or thermal destination is inferred.

/// Payload release and ejection retaining material closures and ideal joint rows.
#[path = "bilateral/release.rs"]
pub mod release;

use super::super::spatial_loads::implicit::distance::{
    BilateralImpulseConfig, BilateralImpulseResponse, BilateralLoop,
};
use super::*;
use fs_math::dd::Dd;

/// Original capture with the native event that transmitted its mismatch.
#[derive(Debug, Clone, PartialEq)]
pub struct BilateralPayloadCapture {
    /// Actual two-body-before/loaded-tree-after inventory and common-frame audit.
    pub capture: PayloadCapture,
    /// Auxiliary loaded-tree mismatch event, absent for exact legacy delegation
    /// when no material loops were supplied. Its initial state carries the
    /// payload rigidly; it is not the physical pre-capture inventory.
    /// Native source/reaction work here is in the lab frame. The complete
    /// inventory reaction work uses capture.balance.accounting_boost_world.
    pub reaction: Option<BilateralImpulseResponse>,
}

/// One simultaneous inventory capture with its shared mixed-constraint event.
#[derive(Debug, Clone, PartialEq)]
pub struct BilateralPayloadBatchCapture {
    /// Complete many-body-before/loaded-tree-after inventory and physical audit.
    /// Attachment handles and payload states retain the authored input order.
    pub capture: PayloadBatchCapture,
    /// Auxiliary loaded-tree mismatch event, absent for exact empty-loop
    /// delegation. Its initial state carries every payload rigidly, and its
    /// lab-frame work is distinct from the complete inventory's energy account.
    pub reaction: Option<BilateralImpulseResponse>,
}

/// Capture while retaining distance, spherical and constant gear/lock rows.
///
/// The original incoming payload pose fixes its attachment. One native loaded
/// impulse solve transmits its momentum mismatch through all simultaneous rows,
/// with full inertia and carrier recoil. Only the combined physical endpoint
/// receives native speed admission. Both initial and final constraint states
/// must meet their declared tolerances; a capture does not engage a new row.
///
/// The original capture quadratic remains nonnegative. Its common-frame
/// identity is delta_K + Q = J_constraint dot v_after; finite tangent defects
/// remain in that measured residual. Actual mean-velocity constraint work is
/// separately capped and never subtracted from Q or counted as damping heat.
///
/// The mixed event's scaled-pivot admission and the existing controls' actual
/// reciprocal-condition floor both apply. Empty material loops delegate to the
/// original capture exactly, with its original gear admission and controls.
///
/// # Errors
/// Original inventory, native geometry/state/limit, mixed-row regularity,
/// momentum/energy/work budget or cancellation refusal. Inputs remain unchanged.
#[allow(clippy::too_many_arguments)]
pub fn capture_payload(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target_link: usize,
    payload: PayloadState,
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: PayloadControls<'_>,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralPayloadCapture, PayloadError> {
    let result = capture_payloads(
        model,
        state,
        &[PayloadCaptureInput {
            target_link,
            payload,
        }],
        loops,
        rows,
        controls,
        config,
        cancelled,
    )?;
    Ok(BilateralPayloadCapture {
        capture: super::single_capture(result.capture),
        reaction: result.reaction,
    })
}

/// Capture several finite bodies in one closed-mechanism momentum event.
///
/// All arriving momentum mismatches are accumulated on the original accepted
/// geometry. One fully loaded inertia and one native mixed impulse solve then
/// determine the complete carrier/link recoil, retaining cross responses
/// between payloads, material closures and gears. Only that combined endpoint
/// receives physical speed admission. The original batch owner keeps the
/// inventory stack and the complete nonnegative capture-loss quadratic.
///
/// Repeated targets are allowed; their handles must later be released in reverse
/// input order on each link. Every actual carried mismatch is formed in the same
/// payload axes and translating frame. A zero total source preserves admitted
/// finite tangency defects after combined geometry/rank/condition admission.
/// Source-free capture may still remove incoming relative energy when opposing
/// payload mismatches cancel; that loss remains in the physical inventory.
///
/// Empty material loops reproduce the original batch capture exactly. The
/// existing capture and constraint-work budgets apply to the whole event.
///
/// # Errors
/// Invalid/empty inventory, native state or limits, unresolved mixed rows,
/// momentum/energy/work admission or cancellation leaves all inputs unchanged.
#[allow(clippy::too_many_arguments)]
pub fn capture_payloads(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadCaptureInput],
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
    controls: PayloadControls<'_>,
    config: BilateralImpulseConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralPayloadBatchCapture, PayloadError> {
    super::batch::capture_evaluated(
        model,
        state,
        inputs,
        rows,
        controls,
        (!loops.is_empty()).then_some((loops, config)),
        cancelled,
    )
}

/// The inventory's existing condition policy also bounds the complete mixed row
/// mobility. A Cholesky pivot floor alone is not this infinity-norm condition.
pub(super) fn admit_condition(
    reciprocal_condition: f64,
    controls: PayloadControls<'_>,
) -> Result<(), PayloadError> {
    if !reciprocal_condition.is_finite() {
        return Err(PayloadError::Invalid("payload mixed reciprocal condition"));
    }
    if reciprocal_condition < controls.minimum_reciprocal_condition {
        return Err(ConstraintError::IllConditioned {
            reciprocal_condition,
            minimum: controls.minimum_reciprocal_condition,
        }
        .into());
    }
    Ok(())
}

/// Project only physical constraint covectors through the original zero-rate
/// RNEA. No dense joint Jacobian or independent kinematics owner is introduced.
pub(super) fn project_constraint_impulses(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    body_impulses: &[Wrench],
    joint_impulses: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, PayloadError> {
    poll(cancelled)?;
    if body_impulses.len() != model.link_count() || joint_impulses.len() != model.dof_count() {
        return Err(PayloadError::Invalid("payload mixed impulse dimensions"));
    }
    let mut external = Vec::with_capacity(body_impulses.len());
    for (index, impulse) in body_impulses.iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        external.push(Wrench::new(
            impulse.torque.scale(-1.0),
            impulse.force.scale(-1.0),
        ));
    }
    let zero = vec![0.0; model.dof_count()];
    let evaluated = inverse_dynamics(
        model,
        BaseState::prescribed(state.base.world_from_base, Twist::zero(), Twist::zero()),
        &state.joint_position,
        &zero,
        &zero,
        Vec3::default(),
        &external,
    )?;
    poll(cancelled)?;
    let base = evaluated.body_wrench[0].transform_by(&model.links()[0].parent_from_child_zero());
    let mut generalized = base.to_array().to_vec();
    for (index, (&body, &joint)) in evaluated
        .generalized_force
        .iter()
        .zip(joint_impulses)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        generalized.push(finite(body + joint, "payload mixed generalized impulse")?);
    }
    validate_slice(&generalized, "payload mixed generalized impulse")?;
    poll(cancelled)?;
    Ok(generalized)
}

/// Constraint work in the same common translating frame as the inventory's
/// source and energy accounts. The native auxiliary event retains lab-frame
/// work separately, including any rounded resultant of internal body rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn constraint_work(
    body_impulses: &[Wrench],
    joint_impulses: &[f64],
    before: &Kinematics,
    after: &Kinematics,
    joint_velocity_before: &[f64],
    joint_velocity_after: &[f64],
    boost: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, PayloadError> {
    poll(cancelled)?;
    if body_impulses.len() != before.body_twist.len()
        || body_impulses.len() != after.body_twist.len()
        || body_impulses.len() != before.world_from_link.len()
        || body_impulses.len() != after.world_from_link.len()
        || joint_impulses.len() != joint_velocity_before.len()
        || joint_impulses.len() != joint_velocity_after.len()
    {
        return Err(PayloadError::Invalid("payload mixed work dimensions"));
    }
    let mut work = Dd::ZERO;
    for (index, impulse) in body_impulses.iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let old = twist_in_frame(
            before.body_twist[index],
            before.world_from_link[index],
            boost,
        )?;
        let new = twist_in_frame(after.body_twist[index], after.world_from_link[index], boost)?;
        for ((impulse, old), new) in impulse
            .to_array()
            .into_iter()
            .zip(old.to_array())
            .zip(new.to_array())
        {
            let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
            work = work + Dd::from_f64(impulse) * mean;
        }
    }
    for (index, ((&impulse, &old), &new)) in joint_impulses
        .iter()
        .zip(joint_velocity_before)
        .zip(joint_velocity_after)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        work = work + Dd::from_f64(impulse) * mean;
    }
    poll(cancelled)?;
    finite(work.to_f64(), "payload common-frame mixed constraint work")
}

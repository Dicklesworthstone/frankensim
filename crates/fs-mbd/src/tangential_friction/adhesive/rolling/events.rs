//! Attachment, loaded separation and recontact of a native rolling pair.
//!
//! One prescribed monotone indentation stage precedes the original held-axis
//! midpoint rolling phase. A snap releases the source-owned normal energy and
//! any spring energy remaining at the actual normal fold; it supplies no body
//! impulse. Held COM forces and body-frame torques can drive attached rolling
//! or detached acceleration through the original native momentum owners. Their
//! actual work is separate from prescribed normal work and contact loss. Free
//! angular travel cannot preload a subsequently formed junction. Poses and the frame
//! remain held. This split does not infer normal motion, finite rotation,
//! rolling-induced area relaxation or a destination for the released energy.

use super::{
    HeldGeometry, RigidAdhesiveRollingControls, RigidAdhesiveRollingError,
    RigidAdhesiveRollingReport, SourceJkrRollingContact, moment, poll, vector,
};
use crate::tangential_friction::hertz::{points, relative_spin};
use crate::tangential_friction::{
    RigidContactBinding, add_components, magnitude, rounded_components,
};
use crate::{MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::rolling::events::{
    JkrRollingEventReport, JkrRollingEventState, SourceJkrRollingEvents,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Original source and the physical pair; longitudinal is the rolling axis.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveRollingEventContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceJkrRollingEvents,
}

/// Source-owned branch, accumulated releases and spring history at held geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingEventState {
    source: JkrRollingEventState,
    geometry: HeldGeometry,
}
impl RigidAdhesiveRollingEventState {
    pub const fn source_state(&self) -> &JkrRollingEventState {
        &self.source
    }
}

/// Separate normal-stage releases and actual subsequent native recoil.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingEventReport {
    pub normal: JkrRollingEventReport,
    /// None means no junction exists during the mechanical phase.
    pub mechanical: Option<RigidAdhesiveRollingReport>,
    /// Original native resultant-impulse receipts when no junction exists.
    /// COM forces add no body-frame lever-arm torque.
    pub detached_impulses: Option<[MomentumImpulseReceipt; 2]>,
    /// Held loads in the original binding's A/B order.
    pub external_wrenches: [Wrench; 2],
    /// Actual native kinetic change, including detached acceleration [J].
    pub kinetic_energy_change_j: f64,
    /// Body force/torque work at actual midpoint velocities [J].
    pub external_work_j: f64,
    /// Actual pair increment minus the applied external impulse [N s].
    pub linear_momentum_residual_ns: Vec3,
    /// Orbital-plus-spin increment minus the external moment impulse [N m s].
    pub angular_momentum_residual_nms: Vec3,
    /// Actual midpoint rolling travel, or free travel while detached [rad].
    pub angle_increment_rad: f64,
    /// Native kinetic plus contact storage increment [J].
    pub represented_storage_change_j: f64,
    /// Prescribed normal, signed stiffness and applied body-load work [J].
    pub applied_work_j: f64,
    /// Original normal snap plus remaining rolling spring release [J].
    pub released_energy_j: f64,
    /// Rolling plastic loss only; not the outgoing snap release [J].
    pub irreversible_rolling_loss_j: f64,
    pub energy_residual_j: f64,
    /// Sum of independent stage defects, so opposite errors cannot cancel.
    pub absolute_stage_residual_j: f64,
    /// Allowance only for observing differences of rounded absolute stores [J].
    pub storage_observation_roundoff_j: f64,
}

/// Retryable proposal; no fallible operation remains after final publication.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingEventCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidAdhesiveRollingEventState,
    pub report: RigidAdhesiveRollingEventReport,
}

impl RigidAdhesiveRollingEventContact<'_> {
    /// Bind a newly formed stress-free junction to the actual body geometry.
    pub fn initial_attached_state(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        indentation_m: f64,
        controls: RigidAdhesiveRollingControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveRollingEventState, RigidAdhesiveRollingError> {
        let native = moment::capture::<SourceJkrRollingContact>(
            self.binding,
            self.frame,
            vector(self.frame.longitudinal()),
            bodies,
            properties,
            controls,
            gate,
        )?;
        let source = self.source.initial_attached_state(indentation_m)?;
        poll(gate)?;
        Ok(RigidAdhesiveRollingEventState {
            source,
            geometry: native.geometry,
        })
    }

    /// Bind a separated pair; later attachment starts a fresh rolling spring.
    pub fn initial_detached_state(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        indentation_m: f64,
        controls: RigidAdhesiveRollingControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveRollingEventState, RigidAdhesiveRollingError> {
        let native = moment::capture::<SourceJkrRollingContact>(
            self.binding,
            self.frame,
            vector(self.frame.longitudinal()),
            bodies,
            properties,
            controls,
            gate,
        )?;
        let source = self.source.initial_detached_state(indentation_m)?;
        poll(gate)?;
        Ok(RigidAdhesiveRollingEventState {
            source,
            geometry: native.geometry,
        })
    }
}

/// Prescribe normal loading, then advance unforced rolling for `dt_s`.
///
/// Normal-stage motion has its own work supply and is not reconstructed from
/// body positions. All source, event, geometry, inertia, momentum and work
/// checks finish before any accepted state can change.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_rolling_events(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingEventState,
    contact: RigidAdhesiveRollingEventContact<'_>,
    indentation_m: f64,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingEventCandidate, RigidAdhesiveRollingError> {
    prepare_rigid_adhesive_rolling_events_with_wrenches(
        bodies,
        properties,
        state,
        contact,
        indentation_m,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        gate,
    )
}

/// Prescribe normal loading, then apply held COM forces and body-frame torques
/// to the actual attached or detached body pair. All loads are in A/B order.
/// Attachment uses the original rolling solver; detached motion uses original
/// native resultant impulses and carries no spring preload into recontact.
/// Poses remain held and a driven attached reversal requires time refinement.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_rolling_events_with_wrenches(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingEventState,
    contact: RigidAdhesiveRollingEventContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingEventCandidate, RigidAdhesiveRollingError> {
    let native = moment::capture::<SourceJkrRollingContact>(
        contact.binding,
        contact.frame,
        vector(contact.frame.longitudinal()),
        bodies,
        properties,
        controls,
        gate,
    )?;
    if native.geometry != state.geometry {
        return Err(RigidAdhesiveRollingError::GeometryMismatch);
    }
    if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > controls.maximum_step_s {
        return Err(RigidAdhesiveRollingError::StepLimit);
    }
    let normal_step = contact
        .source
        .advance_normal(&state.source, indentation_m)?;
    let (normal_state, normal_report) = normal_step.into_parts();
    let observed_normal = Dd::from_f64(normal_state.normal_response().reversible_energy_j)
        - Dd::from_f64(state.source.normal_response().reversible_energy_j);
    let observed_rolling = Dd::from_f64(normal_state.rolling_storage_j())
        - Dd::from_f64(state.source.rolling_storage_j());
    let observation_defect = (observed_normal
        - (Dd::from_f64(normal_report.normal_contact_work_j)
            - Dd::from_f64(normal_report.normal_released_energy_j)))
    .abs()
        + (observed_rolling
            - (Dd::from_f64(normal_report.rolling_parameter_work_j)
                - Dd::from_f64(normal_report.rolling_released_energy_j)))
        .abs();
    // These absolute stores are independently rounded observations. Preserve
    // the exact incremental primitives below and allow only their observation
    // gate a cancellation-scale allowance, as in the attached normal owner.
    let observation_roundoff = [
        state.source.normal_response().reversible_energy_j,
        normal_state.normal_response().reversible_energy_j,
        state.source.rolling_storage_j(),
        normal_state.rolling_storage_j(),
    ]
    .into_iter()
    .fold(Dd::ZERO, |sum, energy| {
        sum + Dd::from_f64(64.0 * f64::EPSILON * energy.abs())
    })
    .to_f64();
    let mut after = [
        bodies[contact.binding.body_a],
        bodies[contact.binding.body_b],
    ];
    let mut detached_impulses = None;
    let (external_work, linear_residual, angular_residual);
    let (source_after, mechanical, angle, kinetic, rolling_storage, loss, mechanical_defect) =
        if let Some((source, memory)) = normal_state.attached() {
            let candidate = moment::prepare(
                native,
                state.geometry,
                source,
                &memory,
                wrenches,
                dt_s,
                controls,
                gate,
            )?;
            let report: RigidAdhesiveRollingReport = candidate.report.into();
            // Replay one accepted scalar angle through the original lifecycle
            // owner. This preserves its sealed cumulative ledgers and never
            // exposes a writable history-adoption interface.
            let replay = contact
                .source
                .advance_rolling(&normal_state, report.angle_increment_rad)?;
            if replay.report().rolling != Some(report.contact)
                || replay.state().attached().map(|(_, s)| s) != Some(candidate.source_after)
            {
                return Err(RigidAdhesiveRollingError::Invalid(
                    "original lifecycle and native rolling endpoints disagree",
                ));
            }
            after = candidate.bodies_after;
            let angle = report.angle_increment_rad;
            let kinetic = report.kinetic_energy_change_j;
            let storage = report.contact.rolling().storage_change_j;
            let loss = report.contact.rolling().irreversible_loss_j;
            let defect = report.total_energy_residual_j.abs();
            external_work = report.external_work_j;
            linear_residual = report.linear_momentum_residual_ns;
            angular_residual = report.angular_momentum_residual_nms;
            let (successor, _) = replay.into_parts();
            (
                successor,
                Some(report),
                angle,
                kinetic,
                storage,
                loss,
                defect,
            )
        } else {
            let detached = prepare_detached(
                after,
                [
                    properties[contact.binding.body_a],
                    properties[contact.binding.body_b],
                ],
                [contact.binding.arm_a_body, contact.binding.arm_b_body],
                vector(contact.frame.longitudinal()),
                wrenches,
                dt_s,
                controls,
                gate,
            )?;
            after = detached.impulses.map(|r| r.state_after);
            let angle = detached.angle;
            detached_impulses = Some(detached.impulses);
            external_work = detached.work;
            linear_residual = detached.linear_residual;
            angular_residual = detached.angular_residual;
            let (successor, _) = contact
                .source
                .advance_rolling(&normal_state, angle)?
                .into_parts();
            (
                successor,
                None,
                angle,
                detached.kinetic,
                0.0,
                0.0,
                detached.defect,
            )
        };
    // Keep original incremental channels separate until the final rounding.
    // In particular, weak spring release must survive a large normal snap.
    let prescribed = Dd::from_f64(normal_report.normal_contact_work_j)
        + Dd::from_f64(normal_report.rolling_parameter_work_j);
    let supplied = prescribed + Dd::from_f64(external_work);
    let released = Dd::from_f64(normal_report.normal_released_energy_j)
        + Dd::from_f64(normal_report.rolling_released_energy_j);
    let storage = prescribed - released + Dd::from_f64(kinetic) + Dd::from_f64(rolling_storage);
    let residual = (storage + released + Dd::from_f64(loss) - supplied).to_f64();
    let absolute_stages = (observation_defect
        + Dd::from_f64(normal_report.energy_residual_j.abs())
        + Dd::from_f64(mechanical_defect))
    .to_f64();
    if [
        supplied.to_f64(),
        released.to_f64(),
        storage.to_f64(),
        residual,
        absolute_stages,
        observation_roundoff,
    ]
    .iter()
    .any(|value| !value.is_finite())
        || residual.abs() > controls.energy_tolerance_j
        || absolute_stages > controls.energy_tolerance_j + observation_roundoff
    {
        return Err(RigidAdhesiveRollingError::Invalid(
            "normal event and native rolling balance exceeds tolerance",
        ));
    }
    poll(gate)?;
    Ok(RigidAdhesiveRollingEventCandidate {
        bodies_after: after,
        state_after: RigidAdhesiveRollingEventState {
            source: source_after,
            geometry: state.geometry,
        },
        report: RigidAdhesiveRollingEventReport {
            normal: normal_report,
            mechanical,
            detached_impulses,
            external_wrenches: wrenches,
            kinetic_energy_change_j: kinetic,
            external_work_j: external_work,
            linear_momentum_residual_ns: linear_residual,
            angular_momentum_residual_nms: angular_residual,
            angle_increment_rad: angle,
            represented_storage_change_j: storage.to_f64(),
            applied_work_j: supplied.to_f64(),
            released_energy_j: released.to_f64(),
            irreversible_rolling_loss_j: loss,
            energy_residual_j: residual,
            absolute_stage_residual_j: absolute_stages,
            storage_observation_roundoff_j: observation_roundoff,
        },
    })
}

struct DetachedPhase {
    impulses: [MomentumImpulseReceipt; 2],
    angle: f64,
    kinetic: f64,
    work: f64,
    defect: f64,
    linear_residual: Vec3,
    angular_residual: Vec3,
}

#[allow(clippy::too_many_arguments)]
fn prepare_detached(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axis: Vec3,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<DetachedPhase, RigidAdhesiveRollingError> {
    let apply = |i: usize| -> Result<_, RigidAdhesiveRollingError> {
        poll(gate)?;
        wrenches[i].validate()?;
        Ok(before[i].apply_momentum_impulse(
            properties[i],
            wrenches[i].force_world.scale(dt_s),
            wrenches[i].torque_body.scale(dt_s),
        )?)
    };
    // Each original native receipt independently retains actual momentum and
    // kinetic change, commanded midpoint work and their compensated defect.
    let impulses = [apply(0)?, apply(1)?];
    let mut kinetic = Dd::ZERO;
    let mut work = Dd::ZERO;
    let mut defect = Dd::ZERO;
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    for receipt in impulses {
        if magnitude(receipt.linear_residual_world) > controls.linear_momentum_tolerance_ns
            || magnitude(receipt.angular_residual_body) > controls.angular_momentum_tolerance_nms
        {
            return Err(RigidAdhesiveRollingError::Invalid(
                "detached applied impulse exceeds native momentum tolerance",
            ));
        }
        kinetic = kinetic + Dd::from_f64(receipt.work.kinetic_energy_change_j);
        work = work + Dd::from_f64(receipt.work.impulse_work_j);
        defect = defect + Dd::from_f64(receipt.work.work_energy_residual_j.abs());
        add_components(&mut linear, receipt.linear_residual_world);
        let pose = receipt.state_before.pose();
        add_components(
            &mut angular,
            pose.position_world().cross(receipt.linear_residual_world),
        );
        add_components(
            &mut angular,
            pose.orientation()
                .rotate_body_to_world(receipt.angular_residual_body),
        );
    }
    let linear_residual = rounded_components(linear);
    let angular_residual = rounded_components(angular);
    // Keep independently recovered per-body work defects as well as the summed
    // reported work balance; neither opposite body defects nor large totals can
    // hide a rejected impulse. The original owner computes all kinetic work.
    defect = defect + (kinetic - work).abs();
    if [kinetic.to_f64(), work.to_f64(), defect.to_f64()]
        .iter()
        .any(|v| !v.is_finite())
        || defect.to_f64() > controls.energy_tolerance_j
        || !linear_residual.is_finite()
        || !angular_residual.is_finite()
        || magnitude(linear_residual) > controls.linear_momentum_tolerance_ns
        || magnitude(angular_residual) > controls.angular_momentum_tolerance_nms
    {
        return Err(RigidAdhesiveRollingError::Invalid(
            "detached applied-load work or momentum balance exceeds tolerance",
        ));
    }
    let initial = points(before, properties, arms)?;
    let final_points = points(impulses.map(|r| r.state_after), properties, arms)?;
    let rate_before = relative_spin(initial).dot(axis);
    let rate_after = relative_spin(final_points).dot(axis);
    let angle = ((Dd::from_f64(rate_before) + Dd::from_f64(rate_after)) * Dd::from_f64(0.5 * dt_s))
        .to_f64();
    if !angle.is_finite() || angle.abs() > controls.maximum_angle_increment_rad {
        return Err(RigidAdhesiveRollingError::StepLimit);
    }
    if (rate_before > 0.0 && rate_after < 0.0) || (rate_before < 0.0 && rate_after > 0.0) {
        // Constant held acceleration can turn between equal-angle endpoints.
        // Admit the actual quadratic turning excursion, not just net travel.
        let scale = rate_before.abs().max(rate_after.abs());
        let first = rate_before.abs() / scale;
        let fraction = first / (first + rate_after.abs() / scale);
        let excursion =
            (Dd::from_f64(rate_before) * Dd::from_f64(0.5 * dt_s) * Dd::from_f64(fraction))
                .to_f64();
        if !excursion.is_finite() || excursion.abs() > controls.maximum_angle_increment_rad {
            return Err(RigidAdhesiveRollingError::StepLimit);
        }
    }
    poll(gate)?;
    Ok(DetachedPhase {
        impulses,
        angle,
        kinetic: kinetic.to_f64(),
        work: work.to_f64(),
        defect: defect.to_f64(),
        linear_residual,
        angular_residual,
    })
}

/// Publish both native momenta and the complete source-owned lifecycle together.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling_events(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingEventState,
    contact: RigidAdhesiveRollingEventContact<'_>,
    indentation_m: f64,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingEventReport, RigidAdhesiveRollingError> {
    apply_rigid_adhesive_rolling_events_with_wrenches(
        bodies,
        properties,
        state,
        contact,
        indentation_m,
        [Wrench::ZERO; 2],
        dt_s,
        controls,
        gate,
    )
}

/// Publish driven body momenta and the original lifecycle only after admission.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling_events_with_wrenches(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingEventState,
    contact: RigidAdhesiveRollingEventContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingEventReport, RigidAdhesiveRollingError> {
    let candidate = prepare_rigid_adhesive_rolling_events_with_wrenches(
        bodies,
        properties,
        state,
        contact,
        indentation_m,
        wrenches,
        dt_s,
        controls,
        gate,
    )?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *state = candidate.state_after;
    Ok(candidate.report)
}

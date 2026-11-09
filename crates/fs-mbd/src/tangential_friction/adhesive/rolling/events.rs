//! Attachment, loaded separation and recontact of an unforced native rolling pair.
//!
//! One prescribed monotone indentation stage precedes the original held-axis
//! midpoint rolling phase. A snap releases the source-owned normal energy and
//! any spring energy remaining at the actual normal fold; it supplies no body
//! impulse. Detached bodies retain their actual momenta, and their free angular
//! travel cannot preload a subsequently formed junction. Poses and the frame
//! remain held. This split does not infer normal motion, finite rotation,
//! rolling-induced area relaxation or a destination for the released energy.

use super::{
    HeldGeometry, RigidAdhesiveRollingControls, RigidAdhesiveRollingError,
    RigidAdhesiveRollingReport, SourceJkrRollingContact, moment, poll, vector,
};
use crate::tangential_friction::RigidContactBinding;
use crate::tangential_friction::hertz::{points, relative_spin};
use crate::{MassProperties, RigidBodyState, Wrench};
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
    /// Actual midpoint rolling travel, or free travel while detached [rad].
    pub angle_increment_rad: f64,
    /// Native kinetic plus contact storage increment [J].
    pub represented_storage_change_j: f64,
    /// Prescribed normal work plus signed rolling stiffness work [J].
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
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
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
    let (source_after, mechanical, angle, kinetic, rolling_storage, loss, mechanical_defect) =
        if let Some((source, memory)) = normal_state.attached() {
            let candidate = moment::prepare(
                native,
                state.geometry,
                source,
                &memory,
                [Wrench::ZERO; 2],
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
            // In the absence of an attached junction and external wrenches,
            // actual momenta are constant. Record free travel, but do not
            // carry it into the elastic memory of a future contact.
            let kinematics = points(
                after,
                [
                    properties[contact.binding.body_a],
                    properties[contact.binding.body_b],
                ],
                [contact.binding.arm_a_body, contact.binding.arm_b_body],
            )?;
            let rate = relative_spin(kinematics).dot(vector(contact.frame.longitudinal()));
            let angle = (Dd::from_f64(dt_s) * Dd::from_f64(rate)).to_f64();
            if !angle.is_finite() || angle.abs() > controls.maximum_angle_increment_rad {
                return Err(RigidAdhesiveRollingError::StepLimit);
            }
            let (successor, _) = contact
                .source
                .advance_rolling(&normal_state, angle)?
                .into_parts();
            (successor, None, angle, 0.0, 0.0, 0.0, 0.0)
        };
    // Keep original incremental channels separate until the final rounding.
    // In particular, weak spring release must survive a large normal snap.
    let supplied = Dd::from_f64(normal_report.normal_contact_work_j)
        + Dd::from_f64(normal_report.rolling_parameter_work_j);
    let released = Dd::from_f64(normal_report.normal_released_energy_j)
        + Dd::from_f64(normal_report.rolling_released_energy_j);
    let storage = supplied - released + Dd::from_f64(kinetic) + Dd::from_f64(rolling_storage);
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
    let candidate = prepare_rigid_adhesive_rolling_events(
        bodies,
        properties,
        state,
        contact,
        indentation_m,
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

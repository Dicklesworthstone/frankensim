//! Native bonded twist through normal-driven separation and fresh recontact.
//!
//! A prescribed monotone normal stage precedes the original held-axis torsion
//! solve. At the original mode-I fold, remaining twist storage and normal snap
//! release have distinct accounts and supply no invented body impulse. Actual
//! COM forces and body-frame torques act on either branch; detached angular
//! travel cannot preload a future junction. All owners publish together.
//!
//! This retains the source's homogeneous bonded-twist approximation and its
//! mode-III/adhesion ceiling. Poses/frame stay held. Normal body motion,
//! reciprocal mixed-mode area stability, annular recruitment and heat transfer
//! require their own physical owners.

use super::{
    HeldGeometry, RigidAdhesiveTorsionControls, RigidAdhesiveTorsionError,
    RigidAdhesiveTorsionReport, SourceJkrTorsionContact, moment, poll, vector,
};
use crate::tangential_friction::RigidContactBinding;
use crate::{MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3, Wrench};
use fs_contact::normal_patch::adhesive::torsion::events::{
    JkrTorsionEventReport, JkrTorsionEventState, SourceJkrTorsionEvents,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Original physical pair and source; the contact normal is the twist axis.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveTorsionEventContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceJkrTorsionEvents,
}

/// Sealed branch, twist history and cumulative losses at the original geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionEventState {
    source: JkrTorsionEventState,
    geometry: HeldGeometry,
}
impl RigidAdhesiveTorsionEventState {
    pub const fn source_state(&self) -> &JkrTorsionEventState {
        &self.source
    }
}

/// Original source releases, actual body work and subsequent calibrated twist.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionEventReport {
    pub normal: JkrTorsionEventReport,
    /// None when the normal stage leaves no attached junction.
    pub mechanical: Option<RigidAdhesiveTorsionReport>,
    /// Original resultant-impulse receipts for a separated body pair.
    pub detached_impulses: Option<[MomentumImpulseReceipt; 2]>,
    pub external_wrenches: [Wrench; 2],
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    /// Actual midpoint twist travel, or unloaded travel when detached [rad].
    pub angle_increment_rad: f64,
    /// Native kinetic plus normal and twist stored-energy increment [J].
    pub represented_storage_change_j: f64,
    /// Prescribed normal, signed twist-stiffness and applied body work [J].
    pub applied_work_j: f64,
    /// Original normal snap plus twist spring remaining at the actual fold [J].
    pub released_energy_j: f64,
    /// Plastic torsional loss only; releases above remain separate [J].
    pub irreversible_torsion_loss_j: f64,
    pub energy_residual_j: f64,
    /// Independent normal, parameter and mechanical defects cannot cancel [J].
    pub absolute_stage_residual_j: f64,
    /// Allowance only for differences of rounded absolute contact stores [J].
    pub storage_observation_roundoff_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionEventCandidate {
    /// Original binding A/B order.
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidAdhesiveTorsionEventState,
    pub report: RigidAdhesiveTorsionEventReport,
}

impl RigidAdhesiveTorsionEventContact<'_> {
    /// Bind a new, unstrained junction to the actual pair and held geometry.
    pub fn initial_attached_state(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        indentation_m: f64,
        controls: RigidAdhesiveTorsionControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTorsionEventState, RigidAdhesiveTorsionError> {
        self.initial_state(bodies, properties, indentation_m, true, controls, gate)
    }

    /// Bind a separated pair. Future attachment starts with zero stored twist.
    pub fn initial_detached_state(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        indentation_m: f64,
        controls: RigidAdhesiveTorsionControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTorsionEventState, RigidAdhesiveTorsionError> {
        self.initial_state(bodies, properties, indentation_m, false, controls, gate)
    }

    fn initial_state(
        self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        indentation_m: f64,
        attached: bool,
        controls: RigidAdhesiveTorsionControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveTorsionEventState, RigidAdhesiveTorsionError> {
        let native = moment::capture::<SourceJkrTorsionContact>(
            self.binding,
            self.frame,
            vector(self.frame.normal()),
            bodies,
            properties,
            controls,
            gate,
        )?;
        let source = if attached {
            self.source.initial_attached_state(indentation_m)
        } else {
            self.source.initial_detached_state(indentation_m)
        }
        .map_err(RigidAdhesiveTorsionError::Source)?;
        poll(gate)?;
        Ok(RigidAdhesiveTorsionEventState {
            source,
            geometry: native.geometry,
        })
    }
}

/// Prescribe a normal stage, then advance the original unforced native twist.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_torsion_events(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveTorsionEventState,
    contact: RigidAdhesiveTorsionEventContact<'_>,
    indentation_m: f64,
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionEventCandidate, RigidAdhesiveTorsionError> {
    prepare_rigid_adhesive_torsion_events_with_wrenches(
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

/// Stage normal-driven branch changes followed by held-force/torque motion.
/// Native attached and detached phases share their existing momentum/work
/// owners with rolling. The source alone supplies branch events, radius,
/// retained twist and loss; public edited observations are never adopted.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_torsion_events_with_wrenches(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveTorsionEventState,
    contact: RigidAdhesiveTorsionEventContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionEventCandidate, RigidAdhesiveTorsionError> {
    let native = moment::capture::<SourceJkrTorsionContact>(
        contact.binding,
        contact.frame,
        vector(contact.frame.normal()),
        bodies,
        properties,
        controls,
        gate,
    )?;
    if native.geometry != state.geometry {
        return Err(RigidAdhesiveTorsionError::GeometryMismatch);
    }
    if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > controls.maximum_step_s {
        return Err(RigidAdhesiveTorsionError::StepLimit);
    }
    let (normal_state, normal) = contact
        .source
        .advance_normal(&state.source, indentation_m)
        .map_err(RigidAdhesiveTorsionError::Source)?
        .into_parts();
    let observed_normal = Dd::from_f64(normal_state.normal_response().reversible_energy_j)
        - Dd::from_f64(state.source.normal_response().reversible_energy_j);
    let observed_twist = Dd::from_f64(normal_state.torsion_storage_j())
        - Dd::from_f64(state.source.torsion_storage_j());
    let observation_defect = (observed_normal
        - (Dd::from_f64(normal.normal_contact_work_j)
            - Dd::from_f64(normal.normal_released_energy_j)))
    .abs()
        + (observed_twist
            - (Dd::from_f64(normal.torsion_parameter_work_j)
                - Dd::from_f64(normal.torsion_released_energy_j)))
        .abs();
    let observation_roundoff = [
        state.source.normal_response().reversible_energy_j,
        normal_state.normal_response().reversible_energy_j,
        state.source.torsion_storage_j(),
        normal_state.torsion_storage_j(),
    ]
    .into_iter()
    .fold(Dd::ZERO, |sum, energy| {
        sum + Dd::from_f64(64.0 * f64::EPSILON * energy.abs())
    })
    .to_f64();
    let (
        after,
        successor,
        mechanical,
        detached_impulses,
        angle,
        kinetic,
        external_work,
        twist_storage,
        loss,
        mechanical_defect,
        linear_residual,
        angular_residual,
    ) = if let Some((source, memory)) = normal_state.attached() {
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
        let report: RigidAdhesiveTorsionReport = candidate.report.into();
        let replay = contact
            .source
            .advance_torsion(&normal_state, report.angle_increment_rad)
            .map_err(RigidAdhesiveTorsionError::Source)?;
        if replay.report().torsion != Some(report.contact)
            || replay.state().attached().map(|(_, s)| s) != Some(candidate.source_after)
        {
            return Err(RigidAdhesiveTorsionError::Invalid(
                "original lifecycle and native torsion endpoints disagree",
            ));
        }
        let angle = report.angle_increment_rad;
        let kinetic = report.kinetic_energy_change_j;
        let external_work = report.external_work_j;
        let storage = report.contact.torsion().storage_change_j;
        let loss = report.contact.torsion().irreversible_loss_j;
        let defect = report.total_energy_residual_j.abs();
        let linear = report.linear_momentum_residual_ns;
        let angular = report.angular_momentum_residual_nms;
        (
            candidate.bodies_after,
            replay.into_parts().0,
            Some(report),
            None,
            angle,
            kinetic,
            external_work,
            storage,
            loss,
            defect,
            linear,
            angular,
        )
    } else {
        let detached = moment::prepare_detached::<SourceJkrTorsionContact>(
            [
                bodies[contact.binding.body_a],
                bodies[contact.binding.body_b],
            ],
            [
                properties[contact.binding.body_a],
                properties[contact.binding.body_b],
            ],
            [contact.binding.arm_a_body, contact.binding.arm_b_body],
            vector(contact.frame.normal()),
            wrenches,
            dt_s,
            controls,
            gate,
        )?;
        let successor = contact
            .source
            .advance_torsion(&normal_state, detached.angle)
            .map_err(RigidAdhesiveTorsionError::Source)?
            .into_parts()
            .0;
        (
            detached.impulses.map(|r| r.state_after),
            successor,
            None,
            Some(detached.impulses),
            detached.angle,
            detached.kinetic,
            detached.work,
            0.0,
            0.0,
            detached.defect,
            detached.linear_residual,
            detached.angular_residual,
        )
    };
    // Keep original incremental normal and torsional channels separate until
    // the final observation. Weak twist release survives large normal work.
    let prescribed =
        Dd::from_f64(normal.normal_contact_work_j) + Dd::from_f64(normal.torsion_parameter_work_j);
    let supplied = prescribed + Dd::from_f64(external_work);
    let released = Dd::from_f64(normal.normal_released_energy_j)
        + Dd::from_f64(normal.torsion_released_energy_j);
    let storage = prescribed - released + Dd::from_f64(kinetic) + Dd::from_f64(twist_storage);
    let residual = (storage + released + Dd::from_f64(loss) - supplied).to_f64();
    let absolute_stages = (observation_defect
        + Dd::from_f64(normal.energy_residual_j.abs())
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
    .any(|v| !v.is_finite())
        || residual.abs() > controls.energy_tolerance_j
        || absolute_stages > controls.energy_tolerance_j + observation_roundoff
    {
        return Err(RigidAdhesiveTorsionError::Invalid(
            "normal event and native torsion balance exceeds tolerance",
        ));
    }
    poll(gate)?;
    Ok(RigidAdhesiveTorsionEventCandidate {
        bodies_after: after,
        state_after: RigidAdhesiveTorsionEventState {
            source: successor,
            geometry: state.geometry,
        },
        report: RigidAdhesiveTorsionEventReport {
            normal,
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
            irreversible_torsion_loss_j: loss,
            energy_residual_j: residual,
            absolute_stage_residual_j: absolute_stages,
            storage_observation_roundoff_j: observation_roundoff,
        },
    })
}

/// Publish both actual native momenta and the complete unforced lifecycle.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_torsion_events(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveTorsionEventState,
    contact: RigidAdhesiveTorsionEventContact<'_>,
    indentation_m: f64,
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionEventReport, RigidAdhesiveTorsionError> {
    apply_rigid_adhesive_torsion_events_with_wrenches(
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

/// Publish driven momenta and branch/twist history after all fallible work.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_torsion_events_with_wrenches(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveTorsionEventState,
    contact: RigidAdhesiveTorsionEventContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionEventReport, RigidAdhesiveTorsionError> {
    let candidate = prepare_rigid_adhesive_torsion_events_with_wrenches(
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

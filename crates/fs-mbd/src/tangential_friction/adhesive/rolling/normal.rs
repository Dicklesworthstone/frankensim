//! Prescribed JKR normal loading followed by the original native rolling phase.
//!
//! The first stage holds rolling angle while the original JKR branch changes
//! its attached radius. Normal contact work and the resulting change of rolling
//! spring energy are separate external supplies. The second stage advances the
//! same retained elastic/plastic memory with the new stiffness and actual body
//! inertias. Neither stage turns recoverable energy into frictional heat.
//!
//! This is an explicit operator split with prescribed local indentation and
//! held body poses/frame. It does not solve normal body motion, supply its
//! impulse, or include rolling-induced changes of the JKR equilibrium area.
//! Detachment, frame transport and changed material/calibration still refuse.

use super::{
    RigidAdhesiveRollingContact, RigidAdhesiveRollingControls, RigidAdhesiveRollingError,
    RigidAdhesiveRollingReport, RigidAdhesiveRollingState, SourceJkrRollingContact, moment, poll,
    vector,
};
use crate::tangential_friction::RigidContactBinding;
use crate::{MassProperties, RigidBodyState, Wrench};
use fs_contact::normal_patch::adhesive::JkrResponse;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Independent contact, rolling-memory and actual native-body work accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingNormalReport {
    pub normal_before: JkrResponse,
    pub normal_after: JkrResponse,
    /// Exact original attached-branch integral F d(indentation), J.
    pub normal_contact_work_j: f64,
    /// Signed supply changing rolling stiffness at fixed elastic shift, J.
    /// This is not normal JKR storage, irreversible rolling loss, or heat.
    pub rolling_parameter_work_j: f64,
    /// Original midpoint recoil, actual body work and irreversible rolling loss.
    pub mechanical: RigidAdhesiveRollingReport,
    /// Native kinetic plus normal and rolling stored-energy increment, J.
    pub represented_storage_change_j: f64,
    /// Normal work + rolling parameter work + applied body-wrench work, J.
    pub applied_work_j: f64,
    /// Storage change + rolling loss - all external work, J.
    pub energy_residual_j: f64,
    /// Sum of absolute normal, area-refresh and native mechanical defects.
    /// Opposite numerical errors cannot cancel to pass the acceptance gate.
    pub absolute_stage_residual_j: f64,
    /// Roundoff allowance for subtracting the rounded absolute contact stores.
    /// Incremental storage/work above retains the original factored primitives.
    pub storage_observation_roundoff_j: f64,
}

/// Retryable original-source, history and finite-body proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingNormalCandidate {
    /// Body endpoints in the original binding's A/B order.
    pub bodies_after: [RigidBodyState; 2],
    pub source_after: SourceJkrRollingContact,
    pub state_after: RigidAdhesiveRollingState,
    pub report: RigidAdhesiveRollingNormalReport,
}

/// Continue an attached normal load without resetting a loaded rolling spring.
///
/// The original source admits the complete monotone attached indentation path.
/// Its material, curvature, calibrated critical shift and small-shift limit stay
/// unchanged. Full native inertia, geometry and momentum/work checks then apply
/// to the ordinary rolling solve. A refusal publishes no source, history or body.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_rolling_with_normal_update(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingNormalCandidate, RigidAdhesiveRollingError> {
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
    let update = contact.source.advance_normal(&state.source, indentation_m)?;
    let normal_before = *contact.source.normal_response();
    let normal_after = *update.source().normal_response();
    let normal_work = update.normal_contact_work_j();
    let rolling_work = update.rolling_parameter_work_j();
    let observed_normal_storage = Dd::from_f64(normal_after.reversible_energy_j)
        - Dd::from_f64(normal_before.reversible_energy_j);
    let observed_rolling_storage = Dd::from_f64(update.rolling_storage_after_j())
        - Dd::from_f64(update.rolling_storage_before_j());
    // Absolute stores are observations rounded independently. They cannot
    // replace the original exact-path increments, particularly for weak work
    // beside large preload. Still cross-check them with a cancellation-scale
    // allowance; the native mechanical work keeps its original tighter gate.
    let observation_roundoff = [
        normal_before.reversible_energy_j,
        normal_after.reversible_energy_j,
        update.rolling_storage_before_j(),
        update.rolling_storage_after_j(),
    ]
    .into_iter()
    .fold(Dd::ZERO, |sum, value| {
        sum + Dd::from_f64(64.0 * f64::EPSILON * value.abs())
    })
    .to_f64();
    let (source_after, refreshed) = update.into_parts();
    let mechanical = moment::prepare(
        native,
        state.geometry,
        &source_after,
        &refreshed,
        wrenches,
        dt_s,
        controls,
        gate,
    )?;
    let report: RigidAdhesiveRollingReport = mechanical.report.into();
    let storage = Dd::from_f64(report.kinetic_energy_change_j)
        + Dd::from_f64(normal_work)
        + Dd::from_f64(rolling_work)
        + Dd::from_f64(report.contact.rolling().storage_change_j);
    let work = Dd::from_f64(normal_work)
        + Dd::from_f64(rolling_work)
        + Dd::from_f64(report.external_work_j);
    let residual =
        (storage + Dd::from_f64(report.contact.rolling().irreversible_loss_j) - work).to_f64();
    let absolute_stages = ((observed_normal_storage - Dd::from_f64(normal_work)).abs()
        + (observed_rolling_storage - Dd::from_f64(rolling_work)).abs()
        + Dd::from_f64(report.total_energy_residual_j.abs()))
    .to_f64();
    if [
        storage.to_f64(),
        work.to_f64(),
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
            "prescribed normal and rolling-memory energy balance exceeds tolerance",
        ));
    }
    poll(gate)?;
    Ok(RigidAdhesiveRollingNormalCandidate {
        bodies_after: mechanical.bodies_after,
        source_after,
        state_after: RigidAdhesiveRollingState {
            source: mechanical.source_after,
            geometry: state.geometry,
        },
        report: RigidAdhesiveRollingNormalReport {
            normal_before,
            normal_after,
            normal_contact_work_j: normal_work,
            rolling_parameter_work_j: rolling_work,
            mechanical: report,
            represented_storage_change_j: storage.to_f64(),
            applied_work_j: work.to_f64(),
            energy_residual_j: residual,
            absolute_stage_residual_j: absolute_stages,
            storage_observation_roundoff_j: observation_roundoff,
        },
    })
}

/// Publish the original source, rolling history and both body momenta together.
/// No normal displacement/impulse or caloric state is silently assigned.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling_with_normal_update(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingState,
    source: &mut SourceJkrRollingContact,
    binding: RigidContactBinding,
    frame: TangentFrame,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingNormalReport, RigidAdhesiveRollingError> {
    let candidate = prepare_rigid_adhesive_rolling_with_normal_update(
        bodies,
        properties,
        state,
        RigidAdhesiveRollingContact {
            binding,
            frame,
            source,
        },
        indentation_m,
        wrenches,
        dt_s,
        controls,
        gate,
    )?;
    poll(gate)?;
    bodies[binding.body_a] = candidate.bodies_after[0];
    bodies[binding.body_b] = candidate.bodies_after[1];
    *source = candidate.source_after;
    *state = candidate.state_after;
    Ok(candidate.report)
}

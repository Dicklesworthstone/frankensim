//! Prescribed normal loading of a retained, calibrated bonded-twist junction.
//!
//! The original JKR owner first changes area at fixed elastic twist. Its normal
//! work and the signed torsional-stiffness work are separate supplies. The
//! original moment solver then advances the loaded spring against actual body
//! inertias and applied wrenches. Every owner publishes together on acceptance.
//!
//! This held-pose split prescribes the normal motion and homogeneous bonded
//! twist; it does not infer normal impulse, stress-free annular recruitment,
//! detachment, frame transport, or reciprocal mixed-mode area relaxation.

use super::{
    RigidAdhesiveTorsionContact, RigidAdhesiveTorsionControls, RigidAdhesiveTorsionError,
    RigidAdhesiveTorsionReport, RigidAdhesiveTorsionState, SourceJkrTorsionContact, moment, poll,
    vector,
};
use crate::tangential_friction::RigidContactBinding;
use crate::{MassProperties, RigidBodyState, Wrench};
use fs_contact::normal_patch::adhesive::JkrResponse;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Original normal work, signed area work, and the actual native twist response.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionNormalReport {
    pub normal_before: JkrResponse,
    pub normal_after: JkrResponse,
    /// Original attached-branch integral of normal force against indentation [J].
    pub normal_contact_work_j: f64,
    /// Signed supply changing bonded-twist stiffness at fixed elastic twist [J].
    /// This is recoverable parameter work, not plastic loss or generated heat.
    pub torsion_parameter_work_j: f64,
    pub mechanical: RigidAdhesiveTorsionReport,
    /// Actual kinetic change plus original normal/torsional storage changes [J].
    pub represented_storage_change_j: f64,
    /// Normal work, torsional parameter work and applied body-wrench work [J].
    pub applied_work_j: f64,
    pub energy_residual_j: f64,
    /// Independent normal, parameter and native defects cannot cancel [J].
    pub absolute_stage_residual_j: f64,
    /// Only the independently rounded absolute-store observations use this [J].
    pub storage_observation_roundoff_j: f64,
}

/// Fully checked endpoints; bodies are in the binding's original A/B order.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveTorsionNormalCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub source_after: SourceJkrTorsionContact,
    pub state_after: RigidAdhesiveTorsionState,
    pub report: RigidAdhesiveTorsionNormalReport,
}

/// Change attached normal area, then advance retained twist under actual loads.
///
/// The original source checks the whole monotone normal path, material and
/// calibration identity, and its mode-III/adhesion approximation ceiling.
/// Changing radius never initializes a new spring or imports caller-authored
/// stiffness. The native solve retains its original finite-inertia, hidden-turn,
/// momentum, bounded-effort and work checks.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_torsion_with_normal_update(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveTorsionState,
    contact: RigidAdhesiveTorsionContact<'_>,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionNormalCandidate, RigidAdhesiveTorsionError> {
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
    let update = contact
        .source
        .advance_normal(&state.source, indentation_m)
        .map_err(RigidAdhesiveTorsionError::Source)?;
    let normal_before = *contact.source.normal_response();
    let normal_after = *update.source().normal_response();
    let normal_work = update.normal_contact_work_j();
    let torsion_work = update.torsion_parameter_work_j();
    let observed_normal_storage = Dd::from_f64(normal_after.reversible_energy_j)
        - Dd::from_f64(normal_before.reversible_energy_j);
    let observed_torsion_storage = Dd::from_f64(update.torsion_storage_after_j())
        - Dd::from_f64(update.torsion_storage_before_j());
    // Retain the original factored incremental work beside large preloads.
    // Independently rounded absolute endpoints remain a separate observation,
    // whose subtraction alone receives a scale-dependent roundoff allowance.
    let observation_roundoff = [
        normal_before.reversible_energy_j,
        normal_after.reversible_energy_j,
        update.torsion_storage_before_j(),
        update.torsion_storage_after_j(),
    ]
    .into_iter()
    .fold(Dd::ZERO, |sum, value| {
        sum + Dd::from_f64(64.0 * f64::EPSILON * value.abs())
    })
    .to_f64();
    let (source_after, refreshed) = update.into_parts();
    let native_after = moment::prepare(
        native,
        state.geometry,
        &source_after,
        &refreshed,
        wrenches,
        dt_s,
        controls,
        gate,
    )?;
    let mechanical: RigidAdhesiveTorsionReport = native_after.report.into();
    let storage = Dd::from_f64(mechanical.kinetic_energy_change_j)
        + Dd::from_f64(normal_work)
        + Dd::from_f64(torsion_work)
        + Dd::from_f64(mechanical.contact.torsion().storage_change_j);
    let work = Dd::from_f64(normal_work)
        + Dd::from_f64(torsion_work)
        + Dd::from_f64(mechanical.external_work_j);
    let residual =
        (storage + Dd::from_f64(mechanical.contact.torsion().irreversible_loss_j) - work).to_f64();
    let absolute_stages = ((observed_normal_storage - Dd::from_f64(normal_work)).abs()
        + (observed_torsion_storage - Dd::from_f64(torsion_work)).abs()
        + Dd::from_f64(mechanical.total_energy_residual_j.abs()))
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
        return Err(RigidAdhesiveTorsionError::Invalid(
            "normal-area update and native torsion failed independent work/storage gates",
        ));
    }
    poll(gate)?;
    Ok(RigidAdhesiveTorsionNormalCandidate {
        bodies_after: native_after.bodies_after,
        source_after,
        state_after: RigidAdhesiveTorsionState {
            source: native_after.source_after,
            geometry: state.geometry,
        },
        report: RigidAdhesiveTorsionNormalReport {
            normal_before,
            normal_after,
            normal_contact_work_j: normal_work,
            torsion_parameter_work_j: torsion_work,
            mechanical,
            represented_storage_change_j: storage.to_f64(),
            applied_work_j: work.to_f64(),
            energy_residual_j: residual,
            absolute_stage_residual_j: absolute_stages,
            storage_observation_roundoff_j: observation_roundoff,
        },
    })
}

/// Publish the original normal source, torsional memory and both body momenta.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_torsion_with_normal_update(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveTorsionState,
    source: &mut SourceJkrTorsionContact,
    binding: RigidContactBinding,
    frame: TangentFrame,
    indentation_m: f64,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveTorsionControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveTorsionNormalReport, RigidAdhesiveTorsionError> {
    let candidate = prepare_rigid_adhesive_torsion_with_normal_update(
        bodies,
        properties,
        state,
        RigidAdhesiveTorsionContact {
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

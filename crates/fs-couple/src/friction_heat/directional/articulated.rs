//! Original source friction in articulated impact and loaded contact motion.
//!
//! The original articulated owner retains joint constraints, floating-base
//! recoil, prescribed support reactions and every cross-contact influence.
//! These adapters sample the same ordered interface claims as the independent
//! rigid-body path and return a complete mechanical/source candidate without
//! mutating any accepted state.

use super::impact::{SourceRigidImpactBinding, SourceRigidImpactError};
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseSample;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_exec::CancelGate;
use fs_mbd::articulated::FreeFloatingImpulseResponse;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated_contacts::friction::ArticulatedImpactKinematics;
use fs_mbd::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticContactResult, ArticulatedEllipticSourceError, ArticulatedImpactGeometry,
    resolve_articulated_elliptic_contacts_with,
};
use fs_mbd::articulated_contacts::mixed::{
    MixedArticulatedContactSystem, MixedArticulatedImpulseResponse,
    resolve_mixed_articulated_elliptic_contacts_with,
};
use fs_mbd::articulated_contacts::step::mixed::MixedArticulatedContactPath;
use fs_mbd::articulated_contacts::step::mixed::elliptic::{
    MixedArticulatedEllipticContactStepResult, MixedArticulatedEllipticSourceStepError,
    step_mixed_articulated_elliptic_contacts_with,
};
use fs_mbd::articulated_contacts::step::{
    ArticulatedContactStepConfig, ArticulatedContactStepError,
};
use fs_mbd::articulated_contacts::trajectory::mixed::MixedArticulatedTrajectorySystem;
use fs_mbd::articulated_contacts::{ArticulatedContactError, ArticulatedContactSystem};

/// Actual incoming contact geometry and the original claims sampled there.
#[derive(Clone, Debug)]
pub struct SourceArticulatedImpactReceipt {
    pub kinematics: ArticulatedImpactKinematics,
    /// Includes the complete candidate history and selected source receipts.
    pub sample: SourceImpulseSample,
}

/// Complete, unpublished response of the original articulated systems.
///
/// The surrounding event owner must publish every returned native system
/// response and active source candidate together. Contact positions and joint
/// coordinates remain unchanged. Inactive histories are absent and stay as
/// accepted; each source receipt identifies its original input index.
#[derive(Clone, Debug)]
pub struct SourceArticulatedImpactCandidate<R = FreeFloatingImpulseResponse> {
    pub mechanical: ArticulatedEllipticContactResult<R>,
    /// Geometrically active contacts in ascending stable contact-ID order.
    pub sources: Vec<SourceArticulatedImpactReceipt>,
}

/// Resolve original material friction through one joint-constrained event.
///
/// Contacts, histories and bindings have the same input order. Sources are
/// sampled once per geometrically active contact, before normal impulses are
/// known, using the actual preimpact tangent speed of the original links.
/// The explicit coefficient role, claim identities, accepted path coordinates,
/// texture and full static/kinetic continuation history retain the existing
/// source sampler's semantics. A separating active row still needs its source:
/// another contact can close it through the original joint/body mobility.
///
/// The conditions callback receives the assembled normalized frame and both
/// physical witness positions. It owns local temperature, pressure and every
/// external query coordinate. No pressure is inferred from an impulse. Texture
/// history must already be transported to the supplied tangent frame using the
/// existing objective transport owner when that frame has changed.
///
/// Every source and mechanical result remains private through cancellation,
/// solver admission and the final native endpoint/momentum/work audit. The
/// callback must stage any external state itself. Returning a candidate neither
/// advances time nor allocates total impact loss to heat, wear or finite slip.
#[allow(clippy::too_many_arguments)]
pub fn resolve_source_articulated_elliptic_contacts<E>(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    config: FrictionalContactConfig,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &ArticulatedImpactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    gate: &CancelGate,
) -> Result<
    SourceArticulatedImpactCandidate,
    ArticulatedEllipticSourceError<SourceRigidImpactError<E>>,
> {
    if histories.len() != contacts.len() || bindings.len() != contacts.len() {
        return Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::Invalid(
                "one original source history and binding per articulated contact required",
            ),
        ));
    }
    // Unit markers let the native callback protocol stage its own transaction
    // without cloning every accepted material history or the returned sample.
    // The complete source candidate is owned once by its final receipt.
    let markers = vec![(); contacts.len()];
    let mut sources = Vec::new();
    let result = resolve_articulated_elliptic_contacts_with(
        systems,
        contacts,
        &markers,
        config,
        |_, kinematics| {
            sample_source(
                histories,
                bindings,
                &mut conditions,
                &mut sources,
                kinematics,
                gate,
            )
        },
        || gate.is_requested(),
    )?;
    // No fallible work follows the original owner's final physical acceptance.
    Ok(SourceArticulatedImpactCandidate {
        mechanical: result.mechanical,
        sources,
    })
}

/// Resolve original material friction across free and prescribed-base systems.
///
/// The same source sampler and history/condition semantics as
/// [`resolve_source_articulated_elliptic_contacts`] apply. Native constrained
/// mobility supplies prescribed support reactions and work without an effective
/// support mass. The returned mechanical `boundary_work` includes both plane
/// and base-support work; under nonunique impulses individual source-work
/// contributions remain representative, as admitted by the original contact law.
/// Source candidates and all native responses remain unpublished until the
/// complete event passes the original endpoint, momentum and work gates.
#[allow(clippy::too_many_arguments)]
pub fn resolve_source_mixed_articulated_elliptic_contacts<E>(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    config: FrictionalContactConfig,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &ArticulatedImpactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    gate: &CancelGate,
) -> Result<
    SourceArticulatedImpactCandidate<MixedArticulatedImpulseResponse>,
    ArticulatedEllipticSourceError<SourceRigidImpactError<E>>,
> {
    if histories.len() != contacts.len() || bindings.len() != contacts.len() {
        return Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::Invalid(
                "one original source history and binding per articulated contact required",
            ),
        ));
    }
    let markers = vec![(); contacts.len()];
    let mut sources = Vec::new();
    let result = resolve_mixed_articulated_elliptic_contacts_with(
        systems,
        contacts,
        &markers,
        config,
        |_, kinematics| {
            sample_source(
                histories,
                bindings,
                &mut conditions,
                &mut sources,
                kinematics,
                gate,
            )
        },
        || gate.is_requested(),
    )?;
    Ok(SourceArticulatedImpactCandidate {
        mechanical: result.mechanical,
        sources,
    })
}

/// Complete material and mechanical candidate for a loaded contact interval.
///
/// The native step retains every load, contact and prescribed-support reaction,
/// the actual endpoint states and clocks, and independent impulse/work accounts.
/// Source candidates stay private through native drift, full-path admission and
/// final cancellation. The caller must adopt all endpoints and source histories
/// together.
pub struct SourceMixedArticulatedContactStepCandidate<P> {
    /// Original native directional force/contact/drift result and path admission.
    pub mechanical: MixedArticulatedEllipticContactStepResult<P>,
    /// Original active material samples in native canonical contact-ID order.
    pub sources: Vec<SourceArticulatedImpactReceipt>,
}

/// Evolve loaded free and supported mechanisms with original material friction.
///
/// Each active interface uses the existing source sampler exactly once, at the
/// actual incoming tangent speed AFTER the held-load kick. The original native
/// simultaneous contact solve, driven support paths and full-path geometry owner
/// determine the entire candidate; no source history escapes on a later refusal.
/// The full native boundary work includes plane and prescribed-base work.
///
/// Conditions, support-path and geometry callbacks must stage any external state
/// until success. Geometry must admit the complete frozen-witness path, not just
/// its endpoints. The source coefficient role is explicit and frozen over this
/// first-order interval. Sampling retains existing travel/texture/history: this
/// adapter does not integrate finite slip, infer pressure from impulse, or
/// deposit contact/splitting loss into heat or wear.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn step_source_mixed_articulated_elliptic_contacts<P, E, SE>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    start_time_s: f64,
    config: ArticulatedContactStepConfig,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &ArticulatedImpactKinematics,
    ) -> Result<SourceFrictionConditions, SE>,
    support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    admit_path: impl FnOnce(&MixedArticulatedContactPath<'_, ()>) -> Result<P, E>,
    gate: &CancelGate,
) -> Result<
    SourceMixedArticulatedContactStepCandidate<P>,
    MixedArticulatedEllipticSourceStepError<E, SourceRigidImpactError<SE>>,
> {
    if histories.len() != contacts.len() || bindings.len() != contacts.len() {
        return Err(MixedArticulatedEllipticSourceStepError::Step(
            ArticulatedContactStepError::Invalid(
                "one original source history and binding per articulated contact required",
            ),
        ));
    }
    let markers = vec![(); contacts.len()];
    let mut sources = Vec::new();
    let result = step_mixed_articulated_elliptic_contacts_with(
        systems,
        contacts,
        &markers,
        start_time_s,
        config,
        |_, kinematics| {
            sample_source(
                histories,
                bindings,
                &mut conditions,
                &mut sources,
                kinematics,
                gate,
            )
        },
        support_path,
        admit_path,
        || gate.is_requested(),
    )?;
    Ok(SourceMixedArticulatedContactStepCandidate {
        mechanical: result.mechanical,
        sources,
    })
}

fn sample_source<E>(
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    conditions: &mut impl FnMut(
        &SourceThermalEllipticState,
        &ArticulatedImpactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    sources: &mut Vec<SourceArticulatedImpactReceipt>,
    kinematics: &ArticulatedImpactKinematics,
    gate: &CancelGate,
) -> Result<(EllipticFriction, ()), SourceRigidImpactError<E>> {
    let index = kinematics.original_index;
    let accepted = &histories[index];
    let binding = bindings[index];
    let boundary = conditions(accepted, kinematics).map_err(SourceRigidImpactError::Conditions)?;
    let sample = binding
        .law
        .sample_preimpact_coefficients(
            accepted,
            [
                kinematics.relative_velocity_m_s[1],
                kinematics.relative_velocity_m_s[2],
            ],
            &boundary,
            binding.coefficient_role,
            &mut || gate.is_requested(),
        )
        .map_err(SourceRigidImpactError::Coefficients)?;
    let friction = EllipticFriction {
        coefficients_mu: sample.coefficients_mu,
        principal_axis: sample.state.constitutive_state().principal_axis(),
    };
    sources.push(SourceArticulatedImpactReceipt {
        kinematics: *kinematics,
        sample,
    });
    Ok((friction, ()))
}

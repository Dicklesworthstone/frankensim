//! Original source friction through simultaneous articulated impact response.
//!
//! The original articulated owner retains joint constraints, floating-base
//! recoil and every cross-contact influence. This adapter samples the same
//! ordered interface claims as the independent rigid-body path and returns one
//! complete mechanical/source candidate, without mutating any accepted state.

use super::impact::{SourceRigidImpactBinding, SourceRigidImpactError};
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseSample;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_exec::CancelGate;
use fs_mbd::articulated_contacts::friction::ArticulatedImpactKinematics;
use fs_mbd::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticContactResult, ArticulatedEllipticSourceError, ArticulatedImpactGeometry,
    resolve_articulated_elliptic_contacts_with,
};
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
pub struct SourceArticulatedImpactCandidate {
    pub mechanical: ArticulatedEllipticContactResult,
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
            let index = kinematics.original_index;
            let accepted = &histories[index];
            let binding = bindings[index];
            let boundary =
                conditions(accepted, kinematics).map_err(SourceRigidImpactError::Conditions)?;
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
        },
        || gate.is_requested(),
    )?;
    // No fallible work follows the original owner's final physical acceptance.
    Ok(SourceArticulatedImpactCandidate {
        mechanical: result.mechanical,
        sources,
    })
}

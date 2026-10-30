//! Original material sources through simultaneous rigid contact impulses.
//!
//! Every active interface samples its incoming physical slip before the one
//! coupled mechanical solve. The original rigid owner assembles all body and
//! support influences, retains normal/tangential coupling and admits actual
//! momentum and work before publishing bodies and source histories together.
//! Collision loss remains unallocated; no finite-area heat route is inferred.

use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceImpulseCoefficientRole, SourceImpulseError, SourceImpulseSample,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::elliptic::{
    RigidEllipticImpactReport, RigidEllipticSourceError, apply_rigid_elliptic_impacts_with,
};
use fs_mbd::frictional_impact::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidImpactControls, RigidImpactKinematics,
};
use fs_mbd::{MassProperties, RigidBodyState};

/// One original ordered interface, in the corresponding input contact's frame.
#[derive(Clone, Copy, Debug)]
pub struct SourceRigidImpactBinding<'a> {
    pub law: &'a SourceThermalEllipticCoulomb,
    /// Explicit source role for the entire instantaneous contact event.
    pub coefficient_role: SourceImpulseCoefficientRole,
}

/// Incoming physical geometry and original claims used by one active contact.
#[derive(Clone, Debug)]
pub struct SourceRigidImpactReceipt {
    pub kinematics: RigidImpactKinematics,
    pub sample: SourceImpulseSample,
}

/// All body and material receipts belong to the same accepted coupled event.
#[derive(Clone, Debug)]
pub struct SourceRigidImpactReport {
    pub mechanical: RigidEllipticImpactReport,
    /// Geometrically active contacts only, in canonical contact-ID order.
    /// `kinematics.original_index` identifies the original binding/history.
    /// Separating active contacts remain present: another contact can close
    /// them through the shared body mobility, so their laws must be sampled.
    pub sources: Vec<SourceRigidImpactReceipt>,
}

#[derive(Debug)]
pub enum SourceRigidImpactError<E> {
    Conditions(E),
    Coefficients(SourceImpulseError),
}
impl<E: core::fmt::Debug> core::fmt::Display for SourceRigidImpactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source simultaneous rigid impact: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SourceRigidImpactError<E> {}

/// Resolve source-held directional friction through the full rigid contact law.
///
/// Contacts, bindings and histories use the same original input order. The
/// existing rigid owner prepares actual material points and normalized frames,
/// queries each geometrically active contact once in canonical ID order, and
/// applies all normal and tangent impulses in one admitted event. Source
/// coefficients are sampled before normal impulses are known. This retains
/// every normal/tangential and cross-contact mobility term; no successive
/// normal solve or independent tangent clip is introduced.
///
/// `conditions` supplies the complete local material query and explicit
/// temperature/pressure. Its point coordinates, incoming velocity and frame
/// come from the original mechanical assembly. No contact pressure is derived
/// from impulse, and no timestep, slip distance or heat destination is invented.
/// Geometrically inactive contacts neither query sources nor change histories.
///
/// Histories retain their texture in each supplied contact tangent frame.
/// A frame change between calls requires the existing objective history
/// transport. Source laws own material claims, while collision geometry and
/// normal restitution remain explicit inputs from their respective owners.
///
/// Query callbacks must keep external state private. Any source/native refusal,
/// work or momentum failure, or cancellation publishes neither body changes
/// nor source history. The resulting total impact loss is an unallocated
/// mechanical receipt; transformed numerical witnesses retain their explicitly
/// declared coordinate basis in the native contact report.
#[allow(clippy::too_many_arguments)]
pub fn apply_source_rigid_elliptic_impacts<E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidFrictionalContact<()>],
    histories: &mut [SourceThermalEllipticState],
    bindings: &[SourceRigidImpactBinding<'_>],
    controls: RigidImpactControls,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &RigidImpactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    gate: &CancelGate,
) -> Result<SourceRigidImpactReport, RigidEllipticSourceError<SourceRigidImpactError<E>>> {
    if bindings.len() != contacts.len() {
        return Err(RigidFrictionalImpactError::Invalid(
            "one original interface binding per contact required",
        ));
    }
    let mut sources = Vec::new();
    let mechanical = apply_rigid_elliptic_impacts_with(
        bodies,
        properties,
        contacts,
        histories,
        controls,
        |accepted, kinematics| {
            let binding = bindings[kinematics.original_index];
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
            let next = sample.state.clone();
            sources.push(SourceRigidImpactReceipt {
                kinematics: *kinematics,
                sample,
            });
            Ok((friction, next))
        },
        gate,
    )?;
    // The original owner has atomically admitted and published both arrays.
    // Constructing this returned receipt performs no further fallible work.
    Ok(SourceRigidImpactReport {
        mechanical,
        sources,
    })
}

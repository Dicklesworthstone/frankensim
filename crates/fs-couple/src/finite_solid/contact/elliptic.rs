//! Original directional interface claims through coupled rigid/solid impulses.
//!
//! Source coefficients are sampled at the actual incoming relative point speed,
//! in the exact contact frame used by the complete simultaneous solve. The
//! original fs-contact owner performs its work-conjugate ellipse transformation;
//! only its physical impulses enter the rigid and consistent solid mass owners.
//! Native endpoint, spin/orbital momentum, support and work gates are shared
//! with circular contact. No source state is published by a rejected event.

use super::{
    Cx, FiniteSolidBody, FiniteSolidRigidContact, FiniteSolidRigidContactControls,
    FiniteSolidRigidContactError, FiniteSolidRigidContactKinematics, FiniteSolidRigidContactReport,
    ImpactResponse, PhysicalImpulseView,
};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactError,
    EllipticFrictionalContactResult, solve_elliptic_frictional_contacts,
};
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceImpulseError, SourceImpulseSample,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_mbd::{MassProperties, RigidBodyState};

/// Reuse the original ordered-interface binding and explicit impact role.
pub use crate::friction_heat::directional::impact::SourceRigidImpactBinding as SourceFiniteSolidRigidContactBinding;

/// Geometry and normal law before the original source resolves the ellipse.
pub type FiniteSolidRigidContactGeometry = FiniteSolidRigidContact<()>;
pub type FiniteSolidRigidEllipticContact = FiniteSolidRigidContact<EllipticFriction>;
pub type FiniteSolidRigidEllipticContactReport =
    FiniteSolidRigidContactReport<EllipticFrictionalContactResult>;
pub type FiniteSolidRigidEllipticContactError =
    FiniteSolidRigidContactError<EllipticFrictionalContactResult, EllipticFrictionalContactError>;

/// Actual incoming geometry and complete original coefficient/state receipts.
#[derive(Clone, Debug)]
pub struct SourceFiniteSolidRigidContactReceipt<K = FiniteSolidRigidContactKinematics> {
    pub kinematics: K,
    /// `state` is a candidate in the contact frame, ready only with this event.
    pub sample: SourceImpulseSample,
}

/// One admitted mechanical event and the source candidates that produced it.
#[derive(Clone, Debug)]
pub struct SourceFiniteSolidRigidContactReport<K = FiniteSolidRigidContactKinematics> {
    pub mechanical: FiniteSolidRigidEllipticContactReport,
    /// Geometrically active contacts, in canonical ID order. Each original
    /// history is replaced by `sample.state` at `kinematics.original_index` when
    /// the caller adopts both body candidates. Inactive histories stay unchanged.
    pub sources: Vec<SourceFiniteSolidRigidContactReceipt<K>>,
}

/// Original source or condition refusal with the caller's contact identity.
#[derive(Debug)]
pub enum SourceFiniteSolidRigidContactFailure<E> {
    Conditions {
        contact_id: u64,
        original_index: usize,
        error: E,
    },
    Coefficients {
        contact_id: u64,
        original_index: usize,
        error: SourceImpulseError,
    },
    Contact(EllipticFrictionalContactError),
    Cancelled,
}
pub type SourceFiniteSolidRigidContactError<E> = FiniteSolidRigidContactError<
    EllipticFrictionalContactResult,
    SourceFiniteSolidRigidContactFailure<E>,
>;

impl ImpactResponse for EllipticFrictionalContactResult {
    fn physical(&self) -> PhysicalImpulseView<'_> {
        PhysicalImpulseView {
            impulses: &self.impulses,
            relative_velocities_after: &self.relative_velocities_after,
            relative_work: self.relative_work,
            dissipated_energy: self.dissipated_energy,
        }
    }
}

impl<'a> FiniteSolidBody<'a> {
    /// Resolve one directional impact with all original rigid/solid mobility.
    ///
    /// The ellipse is expressed in the normalized contact tangent frame.
    /// Both zero coefficients use the native frictionless lane; a single zero
    /// principal coefficient remains the contact owner's singular-law refusal.
    /// Its transformed determinacy witness stays labeled in the contact report.
    /// A unique velocity does not promote representative redundant impulses
    /// into unique per-contact work or heat partitions.
    #[allow(clippy::too_many_arguments)]
    pub fn rigid_elliptic_contact_candidate(
        &self,
        rigid_bodies: &[RigidBodyState],
        rigid_properties: &[MassProperties],
        contacts: &[FiniteSolidRigidEllipticContact],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            FiniteSolidRigidEllipticContactReport,
        ),
        FiniteSolidRigidEllipticContactError,
    > {
        self.rigid_contact_candidate_with(
            rigid_bodies,
            rigid_properties,
            contacts,
            origin_m,
            controls,
            cx,
            |kinematics, matrix| {
                let laws: Vec<_> = kinematics
                    .iter()
                    .map(|query| {
                        let c = contacts[query.original_index];
                        EllipticFrictionalContact {
                            id: c.id,
                            gap: query.gap_m,
                            activation_tolerance: c.activation_tolerance_m,
                            relative_velocity: query.relative_velocity_m_s,
                            friction: c.friction_coefficient,
                            normal_law: c.normal_law,
                        }
                    })
                    .collect();
                solve_elliptic_frictional_contacts(&laws, matrix, controls.contact, || {
                    cx.checkpoint().is_err()
                })
            },
        )
    }

    /// Resolve ordered material-interface sources at actual preimpact slip.
    ///
    /// Contacts, histories and bindings share the original input order. Each
    /// geometrically active source is queried once, in canonical contact-ID
    /// order, after native point/frame assembly. A separating active contact
    /// is still sampled because another contact may close it through the shared
    /// mobility. Inactive contacts query no source and retain their history.
    ///
    /// `conditions` supplies the complete admitted local query, temperature and
    /// optional pressure. The source owner replaces declared speed placeholders
    /// with the actual incoming slip norm and retains every accepted path value.
    /// Normal impulse defines neither pressure nor impact duration. This is a
    /// declared fixed preimpact-sampled ellipse, not an implicit speed-dependent
    /// impact law or automatic static/kinetic switching during the event.
    ///
    /// Each accepted history's texture must be expressed in the exact supplied
    /// tangent frame. Changing that frame requires the existing objective
    /// history transport; opaque source texture names do not define rotations.
    /// Source claims, ordered surface identity, authority, domains and context
    /// remain owned and checked by `SourceThermalEllipticCoulomb`.
    ///
    /// The returned source states and both mechanical candidates belong to one
    /// admitted event. No borrowed state changes on source failure, cancellation
    /// or late native audit rejection. Callback-owned external side effects are
    /// outside this transaction. Total impact loss stays unallocated; no heat,
    /// wear increment, elapsed time or trajectory accuracy is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn rigid_source_elliptic_contact_candidate<E>(
        &self,
        rigid_bodies: &[RigidBodyState],
        rigid_properties: &[MassProperties],
        contacts: &[FiniteSolidRigidContactGeometry],
        histories: &[SourceThermalEllipticState],
        bindings: &[SourceFiniteSolidRigidContactBinding<'_>],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        mut conditions: impl FnMut(
            &SourceThermalEllipticState,
            &FiniteSolidRigidContactKinematics,
        ) -> Result<SourceFrictionConditions, E>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            SourceFiniteSolidRigidContactReport,
        ),
        SourceFiniteSolidRigidContactError<E>,
    > {
        if contacts.len() > controls.contact.maximum_contacts
            || controls.contact.maximum_contacts > 21
        {
            return Err(FiniteSolidRigidContactError::Invalid(
                "bounded contact count required",
            ));
        }
        let mixed: Vec<_> = contacts
            .iter()
            .copied()
            .map(super::mixed::FiniteSolidMixedContact::from)
            .collect();
        let convert = |q: super::mixed::FiniteSolidMixedContactKinematics| {
            FiniteSolidRigidContactKinematics {
                original_index: q.original_index,
                contact_id: q.contact_id,
                axes_world: q.axes_world,
                rigid_point_world: q.first_point_world,
                solid_point_world: q.second_point_world,
                gap_m: q.gap_m,
                relative_velocity_m_s: q.relative_velocity_m_s,
            }
        };
        let (solid, rigid, report) = self.mixed_source_elliptic_contact_candidate(
            rigid_bodies,
            rigid_properties,
            &mixed,
            histories,
            bindings,
            origin_m,
            controls,
            |history, query| conditions(history, &convert(*query)),
            cx,
        )?;
        Ok((
            solid,
            rigid,
            SourceFiniteSolidRigidContactReport {
                mechanical: report.mechanical,
                sources: report
                    .sources
                    .into_iter()
                    .map(|receipt| SourceFiniteSolidRigidContactReceipt {
                        kinematics: convert(receipt.kinematics),
                        sample: receipt.sample,
                    })
                    .collect(),
            },
        ))
    }
}

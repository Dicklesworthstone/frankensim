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
pub struct SourceFiniteSolidRigidContactReceipt {
    pub kinematics: FiniteSolidRigidContactKinematics,
    /// `state` is a candidate in the contact frame, ready only with this event.
    pub sample: SourceImpulseSample,
}

/// One admitted mechanical event and the source candidates that produced it.
#[derive(Clone, Debug)]
pub struct SourceFiniteSolidRigidContactReport {
    pub mechanical: FiniteSolidRigidEllipticContactReport,
    /// Geometrically active contacts, in canonical ID order. Each original
    /// history is replaced by `sample.state` at `kinematics.original_index` when
    /// the caller adopts both body candidates. Inactive histories stay unchanged.
    pub sources: Vec<SourceFiniteSolidRigidContactReceipt>,
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
        use SourceFiniteSolidRigidContactFailure as Failure;
        if bindings.len() != contacts.len() || histories.len() != contacts.len() {
            return Err(FiniteSolidRigidContactError::Invalid(
                "one original interface binding and accepted history per contact required",
            ));
        }
        let mut sources = Vec::new();
        let (solid, rigid, mechanical) = self.rigid_contact_candidate_with(
            rigid_bodies,
            rigid_properties,
            contacts,
            origin_m,
            controls,
            cx,
            |kinematics, matrix| {
                let mut laws = Vec::with_capacity(kinematics.len());
                for query in kinematics {
                    if cx.checkpoint().is_err() {
                        return Err(Failure::Cancelled);
                    }
                    let c = contacts[query.original_index];
                    let friction = if query.gap_m > c.activation_tolerance_m {
                        EllipticFriction {
                            coefficients_mu: [0.0; 2],
                            principal_axis: [1.0, 0.0],
                        }
                    } else {
                        let accepted = &histories[query.original_index];
                        let binding = bindings[query.original_index];
                        let boundary =
                            conditions(accepted, query).map_err(|error| Failure::Conditions {
                                contact_id: query.contact_id,
                                original_index: query.original_index,
                                error,
                            })?;
                        let sample = binding
                            .law
                            .sample_preimpact_coefficients(
                                accepted,
                                [
                                    query.relative_velocity_m_s[1],
                                    query.relative_velocity_m_s[2],
                                ],
                                &boundary,
                                binding.coefficient_role,
                                &mut || cx.checkpoint().is_err(),
                            )
                            .map_err(|error| Failure::Coefficients {
                                contact_id: query.contact_id,
                                original_index: query.original_index,
                                error,
                            })?;
                        let friction = EllipticFriction {
                            coefficients_mu: sample.coefficients_mu,
                            principal_axis: sample.state.constitutive_state().principal_axis(),
                        };
                        sources.push(SourceFiniteSolidRigidContactReceipt {
                            kinematics: *query,
                            sample,
                        });
                        friction
                    };
                    laws.push(EllipticFrictionalContact {
                        id: c.id,
                        gap: query.gap_m,
                        activation_tolerance: c.activation_tolerance_m,
                        relative_velocity: query.relative_velocity_m_s,
                        friction,
                        normal_law: c.normal_law,
                    });
                }
                solve_elliptic_frictional_contacts(&laws, matrix, controls.contact, || {
                    cx.checkpoint().is_err()
                })
                .map_err(Failure::Contact)
            },
        )?;
        Ok((
            solid,
            rigid,
            SourceFiniteSolidRigidContactReport {
                mechanical,
                sources,
            },
        ))
    }
}

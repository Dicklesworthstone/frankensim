//! One contact system for rigid/rigid and rigid/deformable mechanical events.
//!
//! A row is the velocity of its first rigid material point minus the velocity
//! of its rigid or solid counterpart. Both rigid incidences and every consistent
//! solid cross-point block enter ONE matrix and ONE original contact solve.
//! Separate pairwise solves cannot reproduce that coupled response.

use super::{
    Cx, FiniteSolidBody, FiniteSolidRigidContact, FiniteSolidRigidContactControls,
    FiniteSolidRigidContactError, FiniteSolidRigidContactReport,
};
use crate::finite_solid::impulse::FiniteSolidPoint;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactError,
    EllipticFrictionalContactResult, solve_elliptic_frictional_contacts,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactError, FrictionalContactResult, NormalVelocityLaw,
    solve_frictional_contacts,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState,
};
use fs_mbd::{MassProperties, RigidBodyState, Vec3};

/// Original ordered source binding and source refusal identity are shared.
pub use super::elliptic::{
    SourceFiniteSolidRigidContactBinding as SourceFiniteSolidMixedContactBinding,
    SourceFiniteSolidRigidContactError as SourceFiniteSolidMixedContactError,
    SourceFiniteSolidRigidContactFailure as SourceFiniteSolidMixedContactFailure,
};
/// Same body-point representation used by native rigid contact geometry.
pub use fs_mbd::frictional_impact::RigidImpactPoint;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FiniteSolidContactCounterpart {
    Solid(FiniteSolidPoint),
    Rigid(RigidImpactPoint),
}

/// Physical paired witnesses; the first point receives +J, its counterpart -J.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidMixedContact<F = f64> {
    pub id: u64,
    pub point: RigidImpactPoint,
    pub counterpart: FiniteSolidContactCounterpart,
    /// Points from the counterpart toward the first point.
    pub normal_world: Vec3,
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: F,
    pub normal_law: NormalVelocityLaw,
}
impl<F> FiniteSolidMixedContact<F> {
    #[must_use]
    pub fn with_friction<T>(self, friction_coefficient: T) -> FiniteSolidMixedContact<T> {
        FiniteSolidMixedContact {
            id: self.id,
            point: self.point,
            counterpart: self.counterpart,
            normal_world: self.normal_world,
            tangent_world: self.tangent_world,
            activation_tolerance_m: self.activation_tolerance_m,
            friction_coefficient,
            normal_law: self.normal_law,
        }
    }
}
impl<F> From<FiniteSolidRigidContact<F>> for FiniteSolidMixedContact<F> {
    fn from(c: FiniteSolidRigidContact<F>) -> Self {
        Self {
            id: c.id,
            point: RigidImpactPoint {
                body: c.rigid_body,
                arm_body: c.arm_body_m,
            },
            counterpart: FiniteSolidContactCounterpart::Solid(c.solid_point),
            normal_world: c.normal_world,
            tangent_world: c.tangent_world,
            activation_tolerance_m: c.activation_tolerance_m,
            friction_coefficient: c.friction_coefficient,
            normal_law: c.normal_law,
        }
    }
}

/// Original native material-point kinematics, independent of counterpart kind.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidMixedContactKinematics {
    pub original_index: usize,
    pub contact_id: u64,
    /// Normal, tangent, bitangent; exactly the frame used by the matrix/solve.
    pub axes_world: [Vec3; 3],
    pub first_point_world: Vec3,
    pub second_point_world: Vec3,
    pub gap_m: f64,
    /// First minus second material-point velocity [m/s].
    pub relative_velocity_m_s: [f64; 3],
}

pub type FiniteSolidMixedContactReport<R = FrictionalContactResult> =
    FiniteSolidRigidContactReport<R>;
pub type FiniteSolidMixedContactError<R = FrictionalContactResult, E = FrictionalContactError> =
    FiniteSolidRigidContactError<R, E>;
pub type FiniteSolidMixedEllipticContact = FiniteSolidMixedContact<EllipticFriction>;
pub type FiniteSolidMixedEllipticContactReport =
    FiniteSolidRigidContactReport<EllipticFrictionalContactResult>;
pub type FiniteSolidMixedEllipticContactError =
    FiniteSolidRigidContactError<EllipticFrictionalContactResult, EllipticFrictionalContactError>;
pub type FiniteSolidMixedContactGeometry = FiniteSolidMixedContact<()>;
pub type SourceFiniteSolidMixedContactReceipt =
    super::elliptic::SourceFiniteSolidRigidContactReceipt<FiniteSolidMixedContactKinematics>;
pub type SourceFiniteSolidMixedContactReport =
    super::elliptic::SourceFiniteSolidRigidContactReport<FiniteSolidMixedContactKinematics>;

impl<'a> FiniteSolidBody<'a> {
    /// Simultaneous circular Coulomb response across both contact pair kinds.
    /// Counterpart kinds do not partition the solve. Original source/material
    /// points, full inertia, solid mass, support and work gates remain shared.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_contact_candidate(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[FiniteSolidMixedContact],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        (Self, Vec<RigidBodyState>, FiniteSolidMixedContactReport),
        FiniteSolidMixedContactError,
    > {
        self.mixed_contact_candidate_with(
            bodies,
            properties,
            contacts,
            origin_m,
            controls,
            cx,
            |queries, matrix| {
                let laws: Vec<_> = queries
                    .iter()
                    .map(|q| {
                        let c = contacts[q.original_index];
                        FrictionalContact {
                            id: c.id,
                            gap: q.gap_m,
                            activation_tolerance: c.activation_tolerance_m,
                            relative_velocity: q.relative_velocity_m_s,
                            friction_coefficient: c.friction_coefficient,
                            normal_law: c.normal_law,
                        }
                    })
                    .collect();
                solve_frictional_contacts(&laws, matrix, controls.contact, || {
                    cx.checkpoint().is_err()
                })
            },
        )
    }

    /// Simultaneous directional response with the original physical ellipse law.
    /// The original transformed determinacy receipt stays in the result; only
    /// physical impulses and velocities are applied to original native bodies.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_elliptic_contact_candidate(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[FiniteSolidMixedEllipticContact],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            FiniteSolidMixedEllipticContactReport,
        ),
        FiniteSolidMixedEllipticContactError,
    > {
        self.mixed_contact_candidate_with(
            bodies,
            properties,
            contacts,
            origin_m,
            controls,
            cx,
            |queries, matrix| {
                let laws: Vec<_> = queries
                    .iter()
                    .map(|q| {
                        let c = contacts[q.original_index];
                        EllipticFrictionalContact {
                            id: c.id,
                            gap: q.gap_m,
                            activation_tolerance: c.activation_tolerance_m,
                            relative_velocity: q.relative_velocity_m_s,
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

    /// Sample original ordered interface sources for one mixed mechanical event.
    ///
    /// Contacts, accepted histories and bindings share the caller's original
    /// order. Every geometrically active contact is sampled once in canonical
    /// ID order, including separating contacts that another row may close.
    /// The actual incoming tangential point speed is supplied to the original
    /// source owner. A normal impulse supplies neither pressure nor duration.
    /// Inactive histories are untouched; accepted texture frames must already
    /// match each query's exact tangent frame through objective transport.
    ///
    /// Returned source candidates belong to the same admitted mechanical event;
    /// source failure, cancellation or native audit rejection publishes none.
    /// Callback-owned side effects are outside this transaction. No automatic
    /// heat partition, slip-distance increment or implicit friction law is made.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_source_elliptic_contact_candidate<E>(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        contacts: &[FiniteSolidMixedContactGeometry],
        histories: &[SourceThermalEllipticState],
        bindings: &[SourceFiniteSolidMixedContactBinding<'_>],
        origin_m: [f64; 3],
        controls: FiniteSolidRigidContactControls,
        mut conditions: impl FnMut(
            &SourceThermalEllipticState,
            &FiniteSolidMixedContactKinematics,
        ) -> Result<SourceFrictionConditions, E>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            SourceFiniteSolidMixedContactReport,
        ),
        SourceFiniteSolidMixedContactError<E>,
    > {
        use SourceFiniteSolidMixedContactFailure as Failure;
        if bindings.len() != contacts.len() || histories.len() != contacts.len() {
            return Err(FiniteSolidRigidContactError::Invalid(
                "one original interface binding and accepted history per contact required",
            ));
        }
        let mut sources = Vec::new();
        let (solid, rigid, mechanical) = self.mixed_contact_candidate_with(
            bodies,
            properties,
            contacts,
            origin_m,
            controls,
            cx,
            |queries, matrix| {
                let mut laws = Vec::with_capacity(queries.len());
                for query in queries {
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
                        sources.push(SourceFiniteSolidMixedContactReceipt {
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
            SourceFiniteSolidMixedContactReport {
                mechanical,
                sources,
            },
        ))
    }
}

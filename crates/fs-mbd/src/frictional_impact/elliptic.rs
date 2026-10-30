//! Directional simultaneous impulses through the original rigid-body owner.
//!
//! An ellipse is specified in each original contact's normalized tangent frame.
//! Normal/tangent and intercontact mobility blocks remain fully coupled. The
//! contact core's congruence solve supplies physical impulses and velocities;
//! its numerical witness stays explicitly in transformed coordinates. Original
//! body inertia, point impulses, support work and all global audits are shared
//! with the circular-friction entry point.

use super::{
    ImpactResponse, PhysicalImpulseView, RigidFrictionalContact, RigidFrictionalImpactError,
    RigidFrictionalImpactReport, RigidImpactControls, RigidImpactKinematics, apply_impacts_with,
};
use crate::{MassProperties, RigidBodyState};
use fs_contact::frictional_multicontact::FrictionalContactDeterminacy;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactError,
    EllipticFrictionalContactResult, solve_elliptic_frictional_contacts,
};
use fs_exec::CancelGate;

/// Directional coefficients replace the scalar coefficient; every original
/// body/point/geometry/normal-law field retains its original meaning.
pub type RigidEllipticContact = RigidFrictionalContact<EllipticFriction>;
/// Geometry and normal law for a source-owned coefficient query.
pub type RigidImpactGeometry = RigidFrictionalContact<()>;
/// Native body audit with a physical ellipse receipt and its labeled transform.
pub type RigidEllipticImpactReport = RigidFrictionalImpactReport<EllipticFrictionalContactResult>;
pub type RigidEllipticImpactError =
    RigidFrictionalImpactError<EllipticFrictionalContactResult, EllipticFrictionalContactError>;

/// Source refusal preserves both contact identity and original binding index.
#[derive(Debug)]
pub enum RigidEllipticSourceFailure<E> {
    Source {
        contact_id: u64,
        original_index: usize,
        source: E,
    },
    Contact(EllipticFrictionalContactError),
    Cancelled,
}
pub type RigidEllipticSourceError<E> =
    RigidFrictionalImpactError<EllipticFrictionalContactResult, RigidEllipticSourceFailure<E>>;

impl ImpactResponse for EllipticFrictionalContactResult {
    fn physical(&self) -> PhysicalImpulseView<'_> {
        PhysicalImpulseView {
            impulses: &self.impulses,
            relative_velocities_after: &self.relative_velocities_after,
            relative_work: self.relative_work,
            dissipated_energy: self.dissipated_energy,
            // An invertible tangent congruence preserves the distinction
            // between unique impulses and a determined mobility image.
            velocity_only: matches!(
                self.transformed.determinacy,
                FrictionalContactDeterminacy::UniqueVelocity(_)
            ),
        }
    }
}

/// Apply one simultaneous anisotropic-friction event at the original poses.
///
/// Principal coefficients and axes are admitted by the original contact core.
/// Both zero coefficients use its frictionless lane; a single zero principal
/// coefficient is an unsupported singular ellipse and refuses. Dissipated
/// impact energy remains total mechanical loss, without an inferred thermal
/// partition, contact pressure or event duration.
pub fn apply_rigid_elliptic_impacts(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidEllipticContact],
    controls: RigidImpactControls,
    gate: &CancelGate,
) -> Result<RigidEllipticImpactReport, RigidEllipticImpactError> {
    apply_impacts_with(
        bodies,
        properties,
        contacts,
        controls,
        gate,
        |kinematics, delassus| {
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
            solve_elliptic_frictional_contacts(&laws, delassus, controls.contact, || {
                gate.is_requested()
            })
        },
    )
}

/// Resolve one material coefficient pair per active contact, then solve the complete
/// original body event and publish source histories with the accepted bodies.
///
/// Geometry and `histories` share the caller's original input order. Source
/// callbacks execute once per geometrically active contact in canonical ID order, after native
/// kinematic/geometry assembly and before the simultaneous solve. Each query
/// carries the actual relative point velocity and the exact normal-first world
/// frame used by that solve. It supplies no guessed normal load, pressure or
/// impact duration. Sources requiring these conditions must obtain independent
/// caller-owned evidence or refuse; preimpact sampling is not an inner
/// velocity-dependent impact law.
/// Geometrically inactive contacts retain their original histories without a
/// query; their actual endpoint velocities still include all other contacts.
///
/// A callback returns its candidate history by value. All candidates remain
/// private through the constitutive solve and the complete native-body audit.
/// Source failure, cancellation, unresolved support work or a late global
/// balance failure discards them without modifying any accepted history/body.
/// Callback-owned external side effects are outside this transaction.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_elliptic_impacts_with<S, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidImpactGeometry],
    histories: &mut [S],
    controls: RigidImpactControls,
    mut resolve: impl FnMut(&S, &RigidImpactKinematics) -> Result<(EllipticFriction, S), E>,
    gate: &CancelGate,
) -> Result<RigidEllipticImpactReport, RigidEllipticSourceError<E>> {
    if histories.len() != contacts.len() {
        return Err(RigidFrictionalImpactError::Invalid(
            "one accepted source history per original contact required",
        ));
    }
    // States are returned by the source, so no clone of potentially large
    // accepted material history is needed merely to stage the transaction.
    let mut candidates = Vec::new();
    let report =
        apply_impacts_with(
            bodies,
            properties,
            contacts,
            controls,
            gate,
            |kinematics, delassus| {
                let mut laws = Vec::with_capacity(kinematics.len());
                candidates.reserve(kinematics.len());
                for query in kinematics {
                    if gate.is_requested() {
                        return Err(RigidEllipticSourceFailure::Cancelled);
                    }
                    let c = contacts[query.original_index];
                    if query.gap_m > c.activation_tolerance_m {
                        laws.push(EllipticFrictionalContact {
                            id: c.id,
                            gap: query.gap_m,
                            activation_tolerance: c.activation_tolerance_m,
                            relative_velocity: query.relative_velocity_m_s,
                            friction: EllipticFriction {
                                coefficients_mu: [0.0; 2],
                                principal_axis: [1.0, 0.0],
                            },
                            normal_law: c.normal_law,
                        });
                        continue;
                    }
                    let (friction, next) = resolve(&histories[query.original_index], query)
                        .map_err(|source| RigidEllipticSourceFailure::Source {
                            contact_id: query.contact_id,
                            original_index: query.original_index,
                            source,
                        })?;
                    laws.push(EllipticFrictionalContact {
                        id: c.id,
                        gap: query.gap_m,
                        activation_tolerance: c.activation_tolerance_m,
                        relative_velocity: query.relative_velocity_m_s,
                        friction,
                        normal_law: c.normal_law,
                    });
                    candidates.push((query.original_index, next));
                }
                solve_elliptic_frictional_contacts(&laws, delassus, controls.contact, || {
                    gate.is_requested()
                })
                .map_err(RigidEllipticSourceFailure::Contact)
            },
        )?;
    // The original kernel's final cancellation poll precedes body publication.
    // From here to publication of every history there are no fallible calls.
    for (index, mut next) in candidates {
        core::mem::swap(&mut histories[index], &mut next);
    }
    Ok(report)
}

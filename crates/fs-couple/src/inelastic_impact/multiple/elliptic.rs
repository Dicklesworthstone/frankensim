//! Directional friction with original plastic, Maxwell and thermal histories.
//!
//! The original elliptic contact owner maps the COMPLETE coupled mobility into
//! work-conjugate coordinates, solves the same non-associated Coulomb law, then
//! returns physical impulses. The native mixed mechanical transaction applies
//! those impulses once per body and checks original endpoint momentum/work.
//! No principal direction is independently clipped and no material state is
//! rebuilt. Cross-point and rigid/rigid contact coupling remain in the matrix.
//!
//! The authored ellipse is expressed in each contact's normalized tangent frame.
//! Equal principal coefficients recover circular friction exactly; a zero
//! principal coefficient retains exact free slip without epsilon traction.
//! These are fixed-pose impact coefficients, not an inferred pressure, time,
//! source texture evolution, friction-heat partition or constitutive update.

use super::*;
use crate::finite_solid::contact::mixed::FiniteSolidMixedContactKinematics;
pub use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFrictionalContact, EllipticFrictionalContactError, EllipticFrictionalContactResult,
    solve_elliptic_frictional_contacts,
};

pub type InelasticMixedEllipticContact<P = TetImpulsePoint> =
    FiniteSolidMixedContact<EllipticFriction, P>;
pub type InelasticMixedEllipticContactReport =
    FiniteSolidRigidContactReport<EllipticFrictionalContactResult, SolidImpulseReport>;
pub type InelasticMixedEllipticContactError<E> = FiniteSolidRigidContactError<
    EllipticFrictionalContactResult,
    EllipticFrictionalContactError,
    E,
    SolidImpulseReport,
>;
pub type InelasticMixedEllipticCandidate<M> =
    InelasticMixedContactCandidate<M, EllipticFrictionalContactResult>;

pub(super) fn solve<P: Copy>(
    contacts: &[InelasticMixedEllipticContact<P>],
    queries: &[FiniteSolidMixedContactKinematics],
    matrix: &[f64],
    controls: FrictionalContactConfig,
    cx: &Cx<'_>,
) -> Result<EllipticFrictionalContactResult, EllipticFrictionalContactError> {
    let laws: Vec<_> = queries
        .iter()
        .map(|query| {
            let contact = contacts[query.original_index];
            EllipticFrictionalContact {
                id: contact.id,
                gap: query.gap_m,
                activation_tolerance: contact.activation_tolerance_m,
                relative_velocity: query.relative_velocity_m_s,
                normal_law: contact.normal_law,
                friction: contact.friction_coefficient,
            }
        })
        .collect();
    solve_elliptic_frictional_contacts(&laws, matrix, controls, || cx.checkpoint().is_err())
}

/// Directional friction with the original quadratic mixed Maxwell solid.
/// Continuous pressure and all material histories survive the single coupled
/// event; signed P2 cross-point response enters the complete physical matrix.
/// The finite-compliance formulation is required by this mass-only entry.
#[allow(clippy::too_many_arguments)]
pub fn mixed_visco_mixed_elliptic_contact_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact<Tet10ImpulsePoint>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedEllipticCandidate<MixedViscoTetMotion<'a>>,
    InelasticMixedEllipticContactError<MixedViscoTetError>,
> {
    native_candidate!(
        @law solid, held, bodies, properties, contacts, origin, controls, cx,
        Some(solid.velocity_mps()), motion,
        |queries, matrix| solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Coupled directional impacts retaining finite J2 deformation and return state.
#[allow(clippy::too_many_arguments)]
pub fn plastic_mixed_elliptic_contact_candidate<'a, L: FiniteJ2Material>(
    solid: &FinitePlasticTetMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedEllipticCandidate<FinitePlasticTetMotion<'a, L>>,
    InelasticMixedEllipticContactError<FinitePlasticTetError>,
> {
    native_candidate!(
        @law solid, held, bodies, properties, contacts, origin, controls, cx,
        Some(solid.velocity_mps()), motion,
        |queries, matrix| solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Same directional law with all accepted finite-strain Maxwell branches.
#[allow(clippy::too_many_arguments)]
pub fn visco_mixed_elliptic_contact_candidate<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedEllipticCandidate<FiniteViscoTetMotion<'a>>,
    InelasticMixedEllipticContactError<FiniteViscoTetError>,
> {
    native_candidate!(
        @law solid, held, bodies, properties, contacts, origin, controls, cx,
        Some(solid.velocity_mps()), motion,
        |queries, matrix| solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Thermal Tet4 event retaining temperature, free/internal energy and plasticity.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_tet_mixed_elliptic_contact_candidate<'a, L: FiniteThermalJ2Material>(
    solid: &FiniteThermoplasticTetSolid<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedEllipticCandidate<FiniteThermoplasticTetSolid<'a, L>>,
    InelasticMixedEllipticContactError<FiniteThermoplasticTetError>,
> {
    native_candidate!(
        @law solid, held, bodies, properties, contacts, origin, controls, cx,
        solid.velocity_mps(), solid,
        |queries, matrix| solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Thermal Hex8 event retaining full/mean-dilatation history and native inertia.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_hex_mixed_elliptic_contact_candidate<'a, L: FiniteThermalJ2Material>(
    solid: &FiniteThermoplasticHexMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact<HexImpulsePoint>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedEllipticCandidate<FiniteThermoplasticHexMotion<'a, L>>,
    InelasticMixedEllipticContactError<FiniteThermoplasticHexError>,
> {
    native_candidate!(
        @law solid, held, bodies, properties, contacts, origin, controls, cx,
        Some(solid.velocity_mps()), motion,
        |queries, matrix| solve(contacts, queries, matrix, controls.contact, cx)
    )
}

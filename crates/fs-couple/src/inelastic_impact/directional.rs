//! Original directional interface sources through retained-history impacts.
//!
//! Each geometrically active interface samples its complete original claims at
//! the actual incoming tangential point speed. One coupled elliptic contact
//! solve then uses the original consistent solid mass and full rigid inertia.
//! Body candidates and source candidates are returned together only after the
//! existing native momentum, support-work and determinacy checks succeed.
//!
//! The explicit coefficient role defines the entire instantaneous ellipse.
//! Source pressure and temperature come from caller-owned conditions; an
//! impulse supplies neither pressure nor duration. Native material memory,
//! accumulated slip/path and texture do not advance. Accepted textures must
//! already be expressed in the exact contact tangent frame. Contact loss stays
//! unallocated; this event neither deposits heat nor invents a slip increment.

use super::multiple::elliptic::InelasticMixedEllipticCandidate;
use super::multiple::{InelasticMixedContactControls, native_candidate};
use crate::finite_solid::contact::FiniteSolidRigidContactError;
use crate::finite_solid::contact::mixed::FiniteSolidMixedContact;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactResult,
    solve_elliptic_frictional_contacts,
};
use fs_exec::Cx;
use fs_material::plastic::finite::{FiniteJ2Material, thermal::FiniteThermalJ2Material};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::{
    plastic::{
        FinitePlasticTetError,
        dynamics::FinitePlasticTetMotion,
        hex::{
            FiniteThermoplasticHexError,
            dynamics::{FiniteThermoplasticHexMotion, HexImpulsePoint},
        },
        thermal::{FiniteThermoplasticTetError, FiniteThermoplasticTetSolid},
    },
    visco::{
        FiniteViscoTetError,
        dynamics::{FiniteViscoTetMotion, SolidImpulseReport, TetImpulsePoint},
    },
};

/// Reuse the original source bindings, identity-bearing receipts and refusals.
pub use crate::finite_solid::contact::mixed::{
    FiniteSolidMixedContactKinematics,
    SourceFiniteSolidMixedContactBinding as SourceInelasticMixedContactBinding,
    SourceFiniteSolidMixedContactFailure as SourceInelasticMixedContactFailure,
    SourceFiniteSolidMixedContactReceipt as SourceInelasticMixedContactReceipt,
};
pub use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticState, impulse::SourceImpulseCoefficientRole,
};

/// Caller-owned physical witnesses and normal law before resolving friction.
pub type SourceInelasticMixedContact<P = TetImpulsePoint> = FiniteSolidMixedContact<(), P>;
/// Both the original source failure and native material-owner failure survive.
pub type SourceInelasticMixedContactError<S, E> = FiniteSolidRigidContactError<
    EllipticFrictionalContactResult,
    SourceInelasticMixedContactFailure<E>,
    S,
    SolidImpulseReport,
>;

#[derive(Debug)]
pub struct SourceInelasticMixedContactCandidate<M> {
    pub mechanical: InelasticMixedEllipticCandidate<M>,
    /// Active contacts in canonical ID order. Adopt `sample.state` at the
    /// original input index only together with the mechanical candidates.
    /// Inactive histories are unchanged. Callback-owned external side effects
    /// remain the caller's responsibility; no accepted input is mutated here.
    pub sources: Vec<SourceInelasticMixedContactReceipt>,
}

#[allow(clippy::too_many_arguments)]
fn sample_and_solve<P: Copy, E>(
    contacts: &[SourceInelasticMixedContact<P>],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceInelasticMixedContactBinding<'_>],
    queries: &[FiniteSolidMixedContactKinematics],
    matrix: &[f64],
    controls: InelasticMixedContactControls,
    conditions: &mut impl FnMut(
        &SourceThermalEllipticState,
        &FiniteSolidMixedContactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    cx: &Cx<'_>,
) -> Result<
    (
        EllipticFrictionalContactResult,
        Vec<SourceInelasticMixedContactReceipt>,
    ),
    SourceInelasticMixedContactFailure<E>,
> {
    use SourceInelasticMixedContactFailure as Failure;
    let mut laws = Vec::with_capacity(queries.len());
    let mut sources = Vec::with_capacity(queries.len());
    for query in queries {
        cx.checkpoint().map_err(|_| Failure::Cancelled)?;
        let contact = contacts[query.original_index];
        let friction = if query.gap_m > contact.activation_tolerance_m {
            EllipticFriction {
                coefficients_mu: [0.0; 2],
                principal_axis: [1.0, 0.0],
            }
        } else {
            let accepted = &histories[query.original_index];
            let binding = bindings[query.original_index];
            let boundary = conditions(accepted, query).map_err(|error| Failure::Conditions {
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
            sources.push(SourceInelasticMixedContactReceipt {
                kinematics: *query,
                sample,
            });
            friction
        };
        laws.push(EllipticFrictionalContact {
            id: contact.id,
            gap: query.gap_m,
            activation_tolerance: contact.activation_tolerance_m,
            relative_velocity: query.relative_velocity_m_s,
            friction,
            normal_law: contact.normal_law,
        });
    }
    let result = solve_elliptic_frictional_contacts(&laws, matrix, controls.contact, || {
        cx.checkpoint().is_err()
    })
    .map_err(Failure::Contact)?;
    Ok((result, sources))
}

macro_rules! source_candidate {
    ($solid:expr, $held:expr, $bodies:expr, $properties:expr, $contacts:expr,
     $histories:expr, $bindings:expr, $origin:expr, $controls:expr, $conditions:expr,
     $cx:expr, $velocity:expr, $candidate_field:ident) => {{
        if $histories.len() != $contacts.len() || $bindings.len() != $contacts.len() {
            return Err(FiniteSolidRigidContactError::Invalid(
                "one original interface binding and accepted history per contact required",
            ));
        }
        let mut sources = Vec::new();
        let mechanical: Result<_, SourceInelasticMixedContactError<_, _>> = (|| {
            native_candidate!(
                @law $solid, $held, $bodies, $properties, $contacts,
                $origin, $controls, $cx, $velocity, $candidate_field,
                |queries, matrix| {
                    let (result, samples) = sample_and_solve(
                        $contacts, $histories, $bindings, queries, matrix,
                        $controls, &mut $conditions, $cx,
                    )?;
                    sources = samples;
                    Ok(result)
                }
            )
        })();
        Ok(SourceInelasticMixedContactCandidate { mechanical: mechanical?, sources })
    }};
}

/// Sample original interface sources and retain every finite J2 return state.
#[allow(clippy::too_many_arguments)]
pub fn plastic_source_elliptic_contact_candidate<'a, L: FiniteJ2Material, E>(
    solid: &FinitePlasticTetMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[SourceInelasticMixedContact],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceInelasticMixedContactBinding<'_>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &FiniteSolidMixedContactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    cx: &Cx<'_>,
) -> Result<
    SourceInelasticMixedContactCandidate<FinitePlasticTetMotion<'a, L>>,
    SourceInelasticMixedContactError<FinitePlasticTetError, E>,
> {
    source_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        histories,
        bindings,
        origin,
        controls,
        conditions,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

/// Source-derived simultaneous friction without advancing Maxwell memory/time.
#[allow(clippy::too_many_arguments)]
pub fn visco_source_elliptic_contact_candidate<'a, E>(
    solid: &FiniteViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[SourceInelasticMixedContact],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceInelasticMixedContactBinding<'_>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &FiniteSolidMixedContactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    cx: &Cx<'_>,
) -> Result<
    SourceInelasticMixedContactCandidate<FiniteViscoTetMotion<'a>>,
    SourceInelasticMixedContactError<FiniteViscoTetError, E>,
> {
    source_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        histories,
        bindings,
        origin,
        controls,
        conditions,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

/// Source-derived Tet4 contact retaining plastic/thermal states and storage.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_tet_source_elliptic_contact_candidate<'a, L: FiniteThermalJ2Material, E>(
    solid: &FiniteThermoplasticTetSolid<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[SourceInelasticMixedContact],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceInelasticMixedContactBinding<'_>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &FiniteSolidMixedContactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    cx: &Cx<'_>,
) -> Result<
    SourceInelasticMixedContactCandidate<FiniteThermoplasticTetSolid<'a, L>>,
    SourceInelasticMixedContactError<FiniteThermoplasticTetError, E>,
> {
    source_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        histories,
        bindings,
        origin,
        controls,
        conditions,
        cx,
        solid.velocity_mps(),
        solid
    )
}

/// Source-derived Hex8 contact through the original trilinear physical points.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_hex_source_elliptic_contact_candidate<'a, L: FiniteThermalJ2Material, E>(
    solid: &FiniteThermoplasticHexMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[SourceInelasticMixedContact<HexImpulsePoint>],
    histories: &[SourceThermalEllipticState],
    bindings: &[SourceInelasticMixedContactBinding<'_>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    mut conditions: impl FnMut(
        &SourceThermalEllipticState,
        &FiniteSolidMixedContactKinematics,
    ) -> Result<SourceFrictionConditions, E>,
    cx: &Cx<'_>,
) -> Result<
    SourceInelasticMixedContactCandidate<FiniteThermoplasticHexMotion<'a, L>>,
    SourceInelasticMixedContactError<FiniteThermoplasticHexError, E>,
> {
    source_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        histories,
        bindings,
        origin,
        controls,
        conditions,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

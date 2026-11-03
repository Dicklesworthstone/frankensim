//! Simultaneous Coulomb impacts with complete native material history.
//!
//! Rigid/rigid and rigid/solid witnesses share the existing coupled mechanical
//! owner: one full cross-point mobility, one original fs-contact solve, and one
//! summed native impulse update per body. The solid's accepted consistent mass,
//! graded density, deformation, internal variables and temperature survive the
//! event. Mobility borrows them and does not clone constitutive history.
//!
//! Held solid DOFs retain their current velocity. Their external work is
//! reported separately from contact dissipation; a moving support requires
//! unique impulses because velocity determinacy alone does not establish that
//! work. Stationary supports retain the contact owner's velocity-only lane.
//! This fixed-pose event supplies no detection, time stepping, penetration
//! repair, persistent constraint or inferred heat deposition.

use crate::finite_solid::contact::{
    FiniteSolidRigidContactControls, FiniteSolidRigidContactError, FiniteSolidRigidContactReport,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResult,
    solve_frictional_contacts,
};
use fs_exec::Cx;
use fs_material::plastic::finite::{FiniteJ2Material, thermal::FiniteThermalJ2Material};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::{
    mixed::visco::{
        MixedViscoTetError,
        dynamics::{MixedViscoTetMotion, Tet10ImpulsePoint},
    },
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
        dynamics::{
            FiniteViscoTetMotion, SolidImpulseControls, SolidImpulseReport, TetImpulsePoint,
        },
    },
};

/// Shared paired-witness representation, with the original native solid point.
pub use crate::finite_solid::contact::mixed::{
    FiniteSolidContactCounterpart, FiniteSolidMixedContact, RigidImpactPoint,
};
pub type InelasticMixedContact<P = TetImpulsePoint> = FiniteSolidMixedContact<f64, P>;
pub type InelasticMixedContactControls = FiniteSolidRigidContactControls<SolidImpulseControls>;
pub type InelasticMixedContactReport =
    FiniteSolidRigidContactReport<FrictionalContactResult, SolidImpulseReport>;
pub type InelasticMixedContactError<E> = FiniteSolidRigidContactError<
    FrictionalContactResult,
    FrictionalContactError,
    E,
    SolidImpulseReport,
>;

/// Native candidates are published together only after all contact/work gates.
#[derive(Debug)]
pub struct InelasticMixedContactCandidate<M, R = FrictionalContactResult, I = SolidImpulseReport> {
    pub solid: M,
    pub rigid_bodies: Vec<RigidBodyState>,
    pub report: FiniteSolidRigidContactReport<R, I>,
}

pub(super) fn solve<P: Copy>(
    contacts: &[InelasticMixedContact<P>],
    queries: &[crate::finite_solid::contact::mixed::FiniteSolidMixedContactKinematics],
    matrix: &[f64],
    controls: FrictionalContactConfig,
    cx: &Cx<'_>,
) -> Result<FrictionalContactResult, FrictionalContactError> {
    let laws: Vec<_> = queries
        .iter()
        .map(|query| {
            let contact = contacts[query.original_index];
            FrictionalContact {
                id: contact.id,
                gap: query.gap_m,
                activation_tolerance: contact.activation_tolerance_m,
                relative_velocity: query.relative_velocity_m_s,
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            }
        })
        .collect();
    solve_frictional_contacts(&laws, matrix, controls, || cx.checkpoint().is_err())
}

pub(super) fn moving_supports(velocity: Option<&[[f64; 3]]>, held: &[usize]) -> bool {
    // Invalid/unbound supports still reach their original native error below.
    velocity.is_some_and(|velocity| {
        held.iter().any(|&dof| {
            velocity
                .get(dof / 3)
                .is_some_and(|node| node[dof % 3] != 0.0)
        })
    })
}

// Only owner-specific point calls and candidate field names vary. All geometry,
// rigid incidence, full inertia, Coulomb solving and audits remain in the
// existing mechanical core; no pairwise collision sweep is introduced here.
macro_rules! native_candidate {
    ($solid:expr, $held:expr, $bodies:expr, $properties:expr, $contacts:expr,
     $origin:expr, $controls:expr, $cx:expr, $velocity:expr, $candidate_field:ident) => {{
        native_candidate!(
            @law $solid, $held, $bodies, $properties, $contacts,
            $origin, $controls, $cx, $velocity, $candidate_field,
            |queries, matrix| $crate::inelastic_impact::multiple::solve(
                $contacts, queries, matrix, $controls.contact, $cx
            )
        )
    }};
    (@law $solid:expr, $held:expr, $bodies:expr, $properties:expr, $contacts:expr,
     $origin:expr, $controls:expr, $cx:expr, $velocity:expr, $candidate_field:ident,
     $solve:expr) => {{
        let (solid, rigid_bodies, report) = $crate::finite_solid::contact::mixed_contact_candidate_with(
            $bodies,
            $properties,
            $contacts,
            $origin,
            $controls,
            $cx,
            $crate::inelastic_impact::multiple::moving_supports($velocity, $held),
            |point| $solid.point_kinematics(point, $cx),
            |rows, maximum_rows| {
                $solid.point_mobility_matrix_per_kg(
                    $held,
                    rows,
                    maximum_rows,
                    $controls.solid_impulse,
                    $cx,
                )
            },
            |impulses| {
                let candidate = $solid.point_impulses_candidate(
                    $held,
                    impulses,
                    $origin,
                    $controls.solid_impulse,
                    $cx,
                )?;
                Ok((candidate.$candidate_field, candidate.impulse))
            },
            |candidate, point| {
                candidate
                    .point_kinematics(point, $cx)
                    .map(|trace| trace.velocity_mps)
            },
            $solve,
        )?;
        Ok($crate::inelastic_impact::multiple::InelasticMixedContactCandidate {
            solid,
            rigid_bodies,
            report,
        })
    }};
}
pub(super) use native_candidate;

/// Directional friction ellipses through these same native material owners.
pub mod elliptic;

/// Pressure-aware impacts for exact Maxwell and conservative quadratic solids.
pub mod volume;

/// Simultaneous Coulomb impacts on the original quadratic mixed Maxwell solid.
///
/// All ten signed shape values, full cross-point response and moving-support
/// work participate in the same coupled solve as the rigid/rigid rows. Physical
/// pressure, deformation and every relaxation history remain unchanged. This
/// entry requires the original finite-compliance volumetric formulation.
#[allow(clippy::too_many_arguments)]
pub fn mixed_visco_mixed_contact_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact<Tet10ImpulsePoint>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedContactCandidate<MixedViscoTetMotion<'a>>,
    InelasticMixedContactError<MixedViscoTetError>,
> {
    native_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

/// Resolve a coupled impact while retaining every accepted plastic state.
#[allow(clippy::too_many_arguments)]
pub fn plastic_mixed_contact_candidate<'a, L: FiniteJ2Material>(
    solid: &FinitePlasticTetMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedContactCandidate<FinitePlasticTetMotion<'a, L>>,
    InelasticMixedContactError<FinitePlasticTetError>,
> {
    native_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

/// Resolve a coupled impact without advancing or resetting any Maxwell branch.
#[allow(clippy::too_many_arguments)]
pub fn visco_mixed_contact_candidate<'a>(
    solid: &FiniteViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedContactCandidate<FiniteViscoTetMotion<'a>>,
    InelasticMixedContactError<FiniteViscoTetError>,
> {
    native_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

/// Native thermal Tet4 event; temperature and free/internal energy stay fixed.
/// The input must already own its native consistent mass and velocity field.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_tet_mixed_contact_candidate<'a, L: FiniteThermalJ2Material>(
    solid: &FiniteThermoplasticTetSolid<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedContactCandidate<FiniteThermoplasticTetSolid<'a, L>>,
    InelasticMixedContactError<FiniteThermoplasticTetError>,
> {
    native_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        solid.velocity_mps(),
        solid
    )
}

/// Native thermal Hex8 event using its original trilinear point trace and mass.
/// Full/mean-dilatation quadrature histories and caloric storage remain intact.
#[allow(clippy::too_many_arguments)]
pub fn thermoplastic_hex_mixed_contact_candidate<'a, L: FiniteThermalJ2Material>(
    solid: &FiniteThermoplasticHexMotion<'a, L>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact<HexImpulsePoint>],
    origin: [f64; 3],
    controls: InelasticMixedContactControls,
    cx: &Cx<'_>,
) -> Result<
    InelasticMixedContactCandidate<FiniteThermoplasticHexMotion<'a, L>>,
    InelasticMixedContactError<FiniteThermoplasticHexError>,
> {
    native_candidate!(
        solid,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        Some(solid.velocity_mps()),
        motion
    )
}

//! Quadratic contact with native pressure constraints and conservative history.
//!
//! Each event prepares the original supported solid response once. Exact
//! incompressibility retains the native mass/pressure factor across every
//! contact direction and the final summed impulse. Conservative finite-volume
//! formulations instead use their original consistent mass. The same coupled
//! mechanical owner supplies rigid incidence, friction and endpoint audits.
//!
//! Pressure multipliers remain impulses [Pa s], separate from retained pressure
//! [Pa]. Their resultant, torque and work defects enter the absolute coupled
//! error account, so they cannot cancel another body's numerical error.

use super::elliptic::InelasticMixedEllipticContact;
use super::{
    FiniteSolidRigidContactControls, FiniteSolidRigidContactError, FiniteSolidRigidContactReport,
    FrictionalContactError, FrictionalContactResult, InelasticMixedContact,
    InelasticMixedContactCandidate, MassProperties, MixedViscoTetError, MixedViscoTetMotion,
    RigidBodyState, Tet10ImpulsePoint,
};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFrictionalContactError, EllipticFrictionalContactResult,
};
use fs_exec::Cx;
use fs_solid::hyper3::mixed::{
    conservative::{ConservativeTetError, ConservativeTetSolid},
    impulse::{MixedImpulseControls, MixedImpulseReport},
};

pub type VolumeContactControls = FiniteSolidRigidContactControls<MixedImpulseControls>;
pub type VolumeContactReport<R = FrictionalContactResult> =
    FiniteSolidRigidContactReport<R, MixedImpulseReport>;
pub type VolumeContactCandidate<M, R = FrictionalContactResult> =
    InelasticMixedContactCandidate<M, R, MixedImpulseReport>;
pub type VolumeContactError<S, R = FrictionalContactResult, E = FrictionalContactError> =
    FiniteSolidRigidContactError<R, E, S, MixedImpulseReport>;

// Only the original motion and candidate field vary. No new contact, pressure,
// support or work solver is introduced by these material-owner adapters.
macro_rules! prepared_candidate {
    ($solid:expr, $response:expr, $held:expr, $bodies:expr, $properties:expr,
     $contacts:expr, $origin:expr, $controls:expr, $cx:expr, $field:ident, $solve:expr) => {{
        let (solid, rigid_bodies, report) =
            crate::finite_solid::contact::mixed_contact_candidate_with(
                $bodies,
                $properties,
                $contacts,
                $origin,
                $controls,
                $cx,
                super::moving_supports(Some($solid.velocity_mps()), $held),
                |point| $solid.point_kinematics(point, $cx),
                |rows, maximum_rows| {
                    $response.point_mobility_matrix_per_kg(rows, maximum_rows, $cx)
                },
                |impulses| {
                    let candidate = $response.point_impulses_candidate(impulses, $origin, $cx)?;
                    Ok((candidate.$field, candidate.impulse))
                },
                |candidate, point| {
                    candidate
                        .point_kinematics(point, $cx)
                        .map(|trace| trace.velocity_mps)
                },
                $solve,
            )?;
        Ok(InelasticMixedContactCandidate {
            solid,
            rigid_bodies,
            report,
        })
    }};
}

/// Simultaneous Coulomb impacts on an exactly incompressible Maxwell solid.
///
/// The original accepted velocity must meet its weak-volume-rate allowance.
/// No projection is used to hide an incompatible pre-event state. Physical
/// pressure, relaxation histories, geometry and time remain unchanged.
#[allow(clippy::too_many_arguments)]
pub fn mixed_visco_contact_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact<Tet10ImpulsePoint>],
    pressure_gauge_node: Option<usize>,
    origin: [f64; 3],
    controls: VolumeContactControls,
    cx: &Cx<'_>,
) -> Result<VolumeContactCandidate<MixedViscoTetMotion<'a>>, VolumeContactError<MixedViscoTetError>>
{
    let response = solid
        .constrained_impulse_response(held, pressure_gauge_node, controls.solid_impulse, cx)
        .map_err(VolumeContactError::<MixedViscoTetError>::Solid)?;
    prepared_candidate!(
        solid,
        response,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        motion,
        |queries, matrix| super::solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Directional Coulomb response through the same exact pressure factor.
/// The original ellipse owner retains cross-axis coupling and zero-traction
/// principal directions; no pressure penalty or independent-axis solve is used.
#[allow(clippy::too_many_arguments)]
pub fn mixed_visco_elliptic_contact_candidate<'a>(
    solid: &MixedViscoTetMotion<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact<Tet10ImpulsePoint>],
    pressure_gauge_node: Option<usize>,
    origin: [f64; 3],
    controls: VolumeContactControls,
    cx: &Cx<'_>,
) -> Result<
    VolumeContactCandidate<MixedViscoTetMotion<'a>, EllipticFrictionalContactResult>,
    VolumeContactError<
        MixedViscoTetError,
        EllipticFrictionalContactResult,
        EllipticFrictionalContactError,
    >,
> {
    let response = solid
        .constrained_impulse_response(held, pressure_gauge_node, controls.solid_impulse, cx)
        .map_err(
            VolumeContactError::<
                MixedViscoTetError,
                EllipticFrictionalContactResult,
                EllipticFrictionalContactError,
            >::Solid,
        )?;
    prepared_candidate!(
        solid,
        response,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        motion,
        |queries, matrix| super::elliptic::solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Couple the original discrete-gradient solid to a simultaneous rigid manifold.
/// All three native formulations retain their original response: consistent
/// mass for finite compliance and the exact pressure KKT for incompressibility.
/// Event work enters the native cumulative account before any result is returned.
#[allow(clippy::too_many_arguments)]
pub fn conservative_contact_candidate<'a>(
    solid: &ConservativeTetSolid<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedContact<Tet10ImpulsePoint>],
    pressure_gauge_node: Option<usize>,
    origin: [f64; 3],
    controls: VolumeContactControls,
    cx: &Cx<'_>,
) -> Result<
    VolumeContactCandidate<ConservativeTetSolid<'a>>,
    VolumeContactError<ConservativeTetError>,
> {
    let response = solid
        .point_impulse_response(held, pressure_gauge_node, controls.solid_impulse, cx)
        .map_err(VolumeContactError::<ConservativeTetError>::Solid)?;
    prepared_candidate!(
        solid,
        response,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        solid,
        |queries, matrix| super::solve(contacts, queries, matrix, controls.contact, cx)
    )
}

/// Original directional friction with conservative pressure and work history.
#[allow(clippy::too_many_arguments)]
pub fn conservative_elliptic_contact_candidate<'a>(
    solid: &ConservativeTetSolid<'a>,
    held: &[usize],
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[InelasticMixedEllipticContact<Tet10ImpulsePoint>],
    pressure_gauge_node: Option<usize>,
    origin: [f64; 3],
    controls: VolumeContactControls,
    cx: &Cx<'_>,
) -> Result<
    VolumeContactCandidate<ConservativeTetSolid<'a>, EllipticFrictionalContactResult>,
    VolumeContactError<
        ConservativeTetError,
        EllipticFrictionalContactResult,
        EllipticFrictionalContactError,
    >,
> {
    let response = solid
        .point_impulse_response(held, pressure_gauge_node, controls.solid_impulse, cx)
        .map_err(
            VolumeContactError::<
                ConservativeTetError,
                EllipticFrictionalContactResult,
                EllipticFrictionalContactError,
            >::Solid,
        )?;
    prepared_candidate!(
        solid,
        response,
        held,
        bodies,
        properties,
        contacts,
        origin,
        controls,
        cx,
        solid,
        |queries, matrix| super::elliptic::solve(contacts, queries, matrix, controls.contact, cx)
    )
}

//! Directional friction through the original force/contact/drift composition.
//!
//! The original body owner supplies the complete coupled impulse response.
//! Coefficients are held for this first-order step; the geometry owner still
//! admits the complete native drift and every feature's applicability. Source
//! candidates cannot publish until that path and all physical audits succeed.

use super::{
    RigidContactPath, RigidContactStepControls, RigidContactStepError, RigidContactStepReport,
    step_with,
};
use crate::frictional_impact::RigidImpactKinematics;
use crate::frictional_impact::elliptic::{
    RigidEllipticContact, RigidEllipticSourceFailure, RigidImpactGeometry,
    apply_rigid_elliptic_impacts, apply_rigid_elliptic_impacts_with,
};
use crate::{Gravity, MassProperties, RigidBodyState, Wrench};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactError, EllipticFrictionalContactResult,
};
use fs_exec::CancelGate;

pub type RigidEllipticContactStepReport<P> =
    RigidContactStepReport<P, EllipticFrictionalContactResult>;
pub type RigidEllipticContactStepError<E> =
    RigidContactStepError<E, EllipticFrictionalContactResult, EllipticFrictionalContactError>;
pub type RigidEllipticSourceContactStepError<E, SE> =
    RigidContactStepError<E, EllipticFrictionalContactResult, RigidEllipticSourceFailure<SE>>;

/// Advance native poses under held loads and a coupled directional friction law.
///
/// The normal law, oriented principal ellipse and contact geometry are explicit
/// caller inputs. All body/normal/tangent couplings are retained in the original
/// simultaneous impulse solve. Equal principal coefficients recover the scalar
/// entry point. The same native force work, midpoint drift, moving-support work,
/// complete path admission and global balances determine step acceptance.
///
/// This retains the original first-order kick/contact/drift model. Its signed
/// splitting defect and total contact loss are not physical friction heat.
/// Geometry and timestep refinement remain caller-owned; no complete Moreau-Jean
/// method, collision discovery or automatic static/kinetic switch is inferred.
#[allow(clippy::too_many_arguments)]
pub fn step_rigid_elliptic_contacts<P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    gravity: Gravity,
    contacts: &[RigidEllipticContact],
    dt_s: f64,
    controls: RigidContactStepControls,
    admit_path: impl FnOnce(&RigidContactPath<'_, EllipticFriction>) -> Result<P, E>,
    gate: &CancelGate,
) -> Result<RigidEllipticContactStepReport<P>, RigidEllipticContactStepError<E>> {
    step_with(
        bodies,
        properties,
        wrenches,
        gravity,
        contacts,
        dt_s,
        controls,
        apply_rigid_elliptic_impacts,
        admit_path,
        gate,
    )
}

/// Resolve active source coefficients at the native post-kick contact velocity.
///
/// Histories and contact geometries retain caller input order. The existing
/// rigid contact owner calls `resolve` once per geometrically active contact in
/// stable contact-ID order, with actual point positions and the exact normalized
/// tangent frame. Coefficients are sampled before the coupled impulses are
/// known. This is a held-coefficient discretization, not an implicit constitutive
/// solve over the final slip or a resolved continuous friction trajectory.
///
/// Candidate histories stay private through impulse response, native drift,
/// endpoint gaps, complete geometric path admission and final cancellation.
/// A late refusal leaves every original body and history unchanged. No cloning
/// of accepted histories is required. Inactive histories remain untouched.
/// Callbacks must stage any state outside these two supplied arrays themselves.
///
/// The callback owns source conditions and any explicit history-update model.
/// This owner does not infer contact pressure, finite sliding distance or heat
/// partition from an impulse. A changed tangent frame requires objective history
/// transport before resolving its directional source law.
#[allow(clippy::too_many_arguments)]
pub fn step_rigid_elliptic_contacts_with<S, P, E, SE>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    gravity: Gravity,
    contacts: &[RigidImpactGeometry],
    histories: &mut [S],
    dt_s: f64,
    controls: RigidContactStepControls,
    mut resolve: impl FnMut(&S, &RigidImpactKinematics) -> Result<(EllipticFriction, S), SE>,
    admit_path: impl FnOnce(&RigidContactPath<'_, ()>) -> Result<P, E>,
    gate: &CancelGate,
) -> Result<RigidEllipticContactStepReport<P>, RigidEllipticSourceContactStepError<E, SE>> {
    super::poll(gate)?;
    if histories.len() != contacts.len() || contacts.len() > controls.impact.maximum_contacts {
        return Err(RigidContactStepError::Invalid(
            "one bounded accepted source history per original contact required",
        ));
    }
    // Unit markers reuse the original source-query/impulse owner without
    // publishing source state at its earlier, impulse-only boundary. Actual
    // source candidates belong to the complete-step transaction below.
    let mut markers = vec![(); contacts.len()];
    let mut candidates = Vec::new();
    let report = step_with(
        bodies,
        properties,
        wrenches,
        gravity,
        contacts,
        dt_s,
        controls,
        |staged, properties, contacts, controls, gate| {
            apply_rigid_elliptic_impacts_with(
                staged,
                properties,
                contacts,
                &mut markers,
                controls,
                |_, query| {
                    let (friction, next) = resolve(&histories[query.original_index], query)?;
                    candidates.push((query.original_index, next));
                    Ok((friction, ()))
                },
                gate,
            )
        },
        admit_path,
        gate,
    )?;
    // The shared step has completed every fallible operation and published
    // its native endpoint. No callback or cancellation poll splits this commit.
    // Keep replaced states alive until every slot is committed, so their
    // destructors cannot interrupt publication between two source histories.
    for (index, next) in &mut candidates {
        core::mem::swap(&mut histories[*index], next);
    }
    Ok(report)
}

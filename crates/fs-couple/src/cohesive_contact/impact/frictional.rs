//! Simultaneous Signorini-Coulomb impacts on retained cohesive bodies.
//!
//! The caller supplies sparse velocity traces in an orthonormal normal-first
//! frame at each actual contact. The original normal and both tangential laws
//! enter one native simultaneous solve with the full contact mobility. The
//! resulting summed body impulses use the same history-preserving owner as the
//! normal-only event. No finite-duration fracture or friction heat is invented.

use fs_exec::Cx;
use fs_solid::cohesive::dynamics::CohesiveMotion;
use fs_solid::cohesive::dynamics::impulse::{ImpulseControls, ImpulseReport};
use fs_sparse::precond::Precond;

use super::{ImpactError, admit_motions, poll, stage_impulses};

mod response;
pub use fs_contact::frictional_multicontact::NormalVelocityLaw;
pub use response::{FrictionalContact, FrictionalResponseControls, FrictionalResponseReport};

#[derive(Clone, Copy, Debug)]
pub struct FrictionalImpactControls {
    pub response: FrictionalResponseControls,
    /// Original native body acceptance, applied once to its entire impulse.
    pub impulse: ImpulseControls,
}

#[derive(Debug)]
pub struct FrictionalImpactReport {
    pub response: FrictionalResponseReport,
    /// Original body order; zero summed impulses preserve the accepted state.
    pub bodies: Vec<Option<ImpulseReport>>,
}

pub struct FrictionalImpactCandidate<'m, 'a> {
    pub motions: Vec<CohesiveMotion<'m, 'a>>,
    pub report: FrictionalImpactReport,
}

/// Resolve normal restitution and tangential maximum dissipation together.
///
/// Friction is the explicit coefficient of this instantaneous Coulomb law;
/// the adapter does not infer it from cohesive strength or damage. Prescribed
/// relative motion supplies its own signed work. All native determinacy,
/// passivity and unresolved-work refusals retain their original meaning.
pub fn impacts_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    preconditioners: &[&dyn Precond],
    controls: FrictionalImpactControls,
    cx: &Cx<'_>,
) -> Result<FrictionalImpactCandidate<'m, 'a>, ImpactError> {
    impacts_candidate_controlled(motions, contacts, preconditioners, controls, &mut || {
        cx.checkpoint().is_err()
    })
}

/// A failure on any body, contact or assembled audit publishes no successor.
/// Original physical clocks must agree; source/material data may be shared,
/// but each input denotes a distinct accepted mechanical body.
pub fn impacts_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[FrictionalContact<'_>],
    preconditioners: &[&dyn Precond],
    controls: FrictionalImpactControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FrictionalImpactCandidate<'m, 'a>, ImpactError> {
    admit_motions(
        motions,
        preconditioners,
        controls.response.maximum_bodies,
        controls.response.maximum_workspace_entries,
        cancelled,
    )?;
    let bodies: Vec<_> = motions
        .iter()
        .zip(preconditioners)
        .map(|(motion, preconditioner)| response::Body {
            mass: motion.mass_operator(),
            velocity: motion.velocity_m_s(),
            preconditioner: *preconditioner,
        })
        .collect();
    let plan = response::prepare(&bodies, contacts, controls.response, cancelled)?;
    drop(bodies);
    let staged = stage_impulses(
        motions,
        &plan.body_impulses,
        preconditioners,
        controls.impulse,
        cancelled,
    )?;
    let final_velocities: Vec<_> = staged
        .motions
        .iter()
        .map(CohesiveMotion::velocity_m_s)
        .collect();
    let report = response::finish(
        plan,
        contacts,
        &final_velocities,
        &staged.kinetic_changes,
        &staged.body_work_residuals,
        controls.response,
        cancelled,
    )?;
    poll(cancelled)?;
    Ok(FrictionalImpactCandidate {
        motions: staged.motions,
        report: FrictionalImpactReport {
            response: report,
            bodies: staged.receipts,
        },
    })
}

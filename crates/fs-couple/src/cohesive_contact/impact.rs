//! Simultaneous impacts between retained cohesive assemblies.
//!
//! Contact discovery supplies signed velocity traces in each original reduced
//! coordinate system. Their transpose scatters the impulses, preserving work
//! conjugacy. Original mass solves retain every cross-contact mobility entry;
//! the existing fs-contact owner resolves the complete event once. Each body
//! receives one summed impulse through its existing history-preserving owner.
//!
//! These are frozen-geometry event traces, not collision detection or trajectory
//! integration. The caller owns witness geometry, a common inertial frame and
//! the event time. Cartesian angular-momentum claims need geometric point maps;
//! arbitrary reduced traces supply only their stated generalized balance.

use fs_exec::Cx;
use fs_solid::cohesive::dynamics::CohesiveMotion;
use fs_solid::cohesive::dynamics::impulse::{ImpulseControls, ImpulseError, ImpulseReport};
use fs_sparse::precond::Precond;

/// Coupled normal and tangential Coulomb impact on the same retained bodies.
pub mod frictional;
mod mechanics;
mod response;
pub use response::{
    MassSolveReport, NormalContact, ResponseControls, ResponseError, ResponseReport, TraceTerm,
};

/// Native contact/mobility controls and each body's original impulse gates.
#[derive(Clone, Copy, Debug)]
pub struct ImpactControls {
    pub response: ResponseControls,
    /// Applied once to each body's summed event impulse, after contact response
    /// is resolved. A refusal rejects the entire multi-body candidate.
    pub impulse: ImpulseControls,
}

/// Original contact receipt and the actual native body momentum/work receipts.
#[derive(Debug)]
pub struct ImpactReport {
    pub response: ResponseReport,
    /// Original body order; `None` means the summed generalized impulse was
    /// exactly zero and the accepted motion was retained without a mass solve.
    pub bodies: Vec<Option<ImpulseReport>>,
}

/// Every retained body after one admitted event, in original body order.
/// Nothing is published until all native body and assembled work gates pass.
pub struct ImpactCandidate<'m, 'a> {
    pub motions: Vec<CohesiveMotion<'m, 'a>>,
    pub report: ImpactReport,
}

#[derive(Debug)]
pub enum ImpactError {
    Invalid(&'static str),
    Cancelled,
    Response(ResponseError),
    Body { body: usize, error: ImpulseError },
}
impl core::fmt::Display for ImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive simultaneous impact refused: {self:?}")
    }
}
impl std::error::Error for ImpactError {}
impl From<ResponseError> for ImpactError {
    fn from(error: ResponseError) -> Self {
        Self::Response(error)
    }
}

/// Resolve all declared normal contacts in one native simultaneous solve.
///
/// Each input denotes a distinct mechanical body with its own accepted state;
/// immutable mass/material data may be shared. Every body's physical clock must
/// agree. Gaps and signed contact traces come from the caller's actual event
/// geometry. An authored restitution coefficient is a system-level impact law,
/// not inferred from cohesive stiffness, strength or fracture energy.
pub fn normal_impacts_candidate<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[NormalContact<'_>],
    preconditioners: &[&dyn Precond],
    controls: ImpactControls,
    cx: &Cx<'_>,
) -> Result<ImpactCandidate<'m, 'a>, ImpactError> {
    normal_impacts_candidate_controlled(motions, contacts, preconditioners, controls, &mut || {
        cx.checkpoint().is_err()
    })
}

/// Predicate counterpart with the same native mechanics and whole-event
/// publication. A refusal on any body leaves every accepted input unchanged.
pub fn normal_impacts_candidate_controlled<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    contacts: &[NormalContact<'_>],
    preconditioners: &[&dyn Precond],
    controls: ImpactControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ImpactCandidate<'m, 'a>, ImpactError> {
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
    Ok(ImpactCandidate {
        motions: staged.motions,
        report: ImpactReport {
            response: report,
            bodies: staged.receipts,
        },
    })
}

pub(super) fn admit_motions(
    motions: &[&CohesiveMotion<'_, '_>],
    preconditioners: &[&dyn Precond],
    maximum_bodies: usize,
    maximum_workspace_entries: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ImpactError> {
    poll(cancelled)?;
    if motions.len() != preconditioners.len() {
        return Err(ImpactError::Invalid(
            "one native preconditioner per original body required",
        ));
    }
    admit_states(
        motions,
        maximum_bodies,
        maximum_workspace_entries,
        cancelled,
    )
}

// Shared accepted-state admission also serves held-force contact, whose native
// midpoint owner retains its own preconditioning and needs no impulse adapter.
pub(super) fn admit_states(
    motions: &[&CohesiveMotion<'_, '_>],
    maximum_bodies: usize,
    maximum_workspace_entries: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ImpactError> {
    poll(cancelled)?;
    if motions.is_empty() || motions.len() > maximum_bodies {
        return Err(ImpactError::Invalid(
            "bounded nonempty original bodies required",
        ));
    }
    if motions
        .len()
        .checked_mul(8)
        .is_none_or(|entries| entries > maximum_workspace_entries)
    {
        return Err(ImpactError::Invalid(
            "body table exceeds response workspace",
        ));
    }
    let time_s = motions[0].time_s();
    for (index, motion) in motions.iter().enumerate() {
        poll(cancelled)?;
        if motion.time_s() != time_s {
            return Err(ImpactError::Invalid(
                "all impacting bodies must share physical time",
            ));
        }
        if motions[..index]
            .iter()
            .any(|previous| core::ptr::eq(*previous, *motion))
        {
            return Err(ImpactError::Invalid(
                "the same accepted body was supplied more than once",
            ));
        }
    }
    Ok(())
}

struct StagedMotions<'m, 'a> {
    motions: Vec<CohesiveMotion<'m, 'a>>,
    receipts: Vec<Option<ImpulseReport>>,
    kinetic_changes: Vec<f64>,
    body_work_residuals: Vec<f64>,
}

// Both contact laws enter the original owner once per body, with one complete
// summed impulse. Material/contact history and whole-event refusal are shared.
fn stage_impulses<'m, 'a>(
    motions: &[&CohesiveMotion<'m, 'a>],
    impulses: &[Vec<f64>],
    preconditioners: &[&dyn Precond],
    controls: ImpulseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StagedMotions<'m, 'a>, ImpactError> {
    if motions.len() != impulses.len() {
        return Err(ImpactError::Invalid(
            "one summed impulse per original body required",
        ));
    }
    let mut candidates = Vec::with_capacity(motions.len());
    let mut receipts = Vec::with_capacity(motions.len());
    let mut kinetic_changes = Vec::with_capacity(motions.len());
    let mut body_work_residuals = Vec::with_capacity(motions.len());
    for (body, (motion, impulse)) in motions.iter().zip(impulses).enumerate() {
        poll(cancelled)?;
        if impulse.iter().all(|value| *value == 0.0) {
            candidates.push((**motion).clone());
            receipts.push(None);
            kinetic_changes.push(0.0);
            body_work_residuals.push(0.0);
            continue;
        }
        let accepted = motion
            .impulse_candidate_controlled(
                impulse,
                &BorrowedPrecond(preconditioners[body]),
                controls,
                cancelled,
            )
            .map_err(|error| ImpactError::Body { body, error })?;
        kinetic_changes.push(accepted.report.kinetic_energy_change_j);
        body_work_residuals.push(accepted.report.work_residual_j);
        receipts.push(Some(accepted.report));
        candidates.push(accepted.motion);
    }
    Ok(StagedMotions {
        motions: candidates,
        receipts,
        kinetic_changes,
        body_work_residuals,
    })
}

struct BorrowedPrecond<'p>(&'p dyn Precond);
impl Precond for BorrowedPrecond<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.0.apply(r, z);
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ImpactError> {
    if cancelled() {
        Err(ImpactError::Cancelled)
    } else {
        Ok(())
    }
}

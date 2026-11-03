//! Original reduced inertia and one simultaneous normal-contact law.
//!
//! Sparse mass response, conjugate distribution and actual-state work are
//! shared with the Coulomb adapter. The original normal law remains unchanged.

use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactResult, NormalContactImpact, solve_normal_impacts,
};
use fs_solver::checked_cg::CheckedCgConfig;

pub(super) use super::mechanics::Body;
use super::mechanics::{
    self, Mobility, MobilityControls, TraceRow, add, mul, nonnegative, poll, positive,
};
pub use super::mechanics::{MassSolveReport, ResponseError, TraceTerm};

/// Caller-owned frozen contact geometry. Positive velocity separates the sides.
/// The trace is strictly ordered by `(body, dof)` without duplicate or zero
/// entries. Arbitrary reduced traces do not certify Cartesian momentum/torque.
#[derive(Clone, Copy, Debug)]
pub struct NormalContact<'t> {
    pub id: u64,
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    pub restitution: f64,
    /// Prescribed contribution `c` in `relative_velocity = B v + c` [m/s].
    /// Its work on the retained bodies is `-c * normal_impulse`.
    pub prescribed_normal_velocity_m_s: f64,
    pub trace: &'t [TraceTerm],
}

#[derive(Clone, Copy, Debug)]
pub struct ResponseControls {
    pub mass_solve: CheckedCgConfig,
    pub contact_solve: MultiImpactConfig,
    pub maximum_bodies: usize,
    pub maximum_trace_terms: usize,
    /// Conservative scalar-equivalent numerical workspace ceiling, including
    /// body tables, native CG histories, contact-solver scratch and reports.
    /// Excludes borrowed operators/preconditioners and accepted-state clones.
    pub maximum_workspace_entries: usize,
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug)]
pub struct ResponseReport {
    pub contact_ids: Vec<u64>,
    /// Complete row-major contact mobility [1/kg], retained without averaging.
    pub contact_mobility_per_kg: Vec<f64>,
    pub mass_solves: Vec<MassSolveReport>,
    /// Native determinacy/passivity receipt; redundant impulses retain the
    /// original representative-only meaning, never individual heat authority.
    pub contact: MultiImpactResult,
    pub relative_velocities_before_m_s: Vec<f64>,
    /// Re-sampled from the actual returned body velocities, including `c`.
    pub relative_velocities_after_m_s: Vec<f64>,
    pub kinetic_change_j: f64,
    pub prescribed_work_j: f64,
    pub relative_impulse_work_j: f64,
    /// `sum Delta K - prescribed_work - actual relative impulse work`.
    pub work_energy_residual_j: f64,
    /// `sum Delta K - prescribed_work + native impact loss`. Any numerical
    /// positive work remains in this physical residual rather than subtracted.
    pub dissipation_residual_j: f64,
    pub maximum_velocity_residual_m_s: f64,
    /// Sum of ABSOLUTE native body event-work defects: opposing errors cannot
    /// cancel inside the whole-assembly energy account.
    pub absolute_body_work_residual_j: f64,
}

pub(super) struct Plan {
    pub body_impulses: Vec<Vec<f64>>,
    contact_ids: Vec<u64>,
    mobility: Mobility,
    contact: MultiImpactResult,
}

pub(super) fn prepare(
    bodies: &[Body<'_>],
    contacts: &[NormalContact<'_>],
    controls: ResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Plan, ResponseError> {
    let contact_workspace_entries = admit_contacts(contacts, controls, cancelled)?;
    let rows = rows(contacts);
    let mobility = mechanics::assemble(
        bodies,
        &rows,
        MobilityControls {
            mass_solve: controls.mass_solve,
            maximum_bodies: controls.maximum_bodies,
            maximum_trace_terms: controls.maximum_trace_terms,
            maximum_workspace_entries: controls.maximum_workspace_entries,
            contact_workspace_entries,
        },
        cancelled,
    )?;
    let native: Vec<_> = contacts
        .iter()
        .zip(&mobility.before)
        .map(|(row, &speed)| NormalContactImpact {
            gap: row.gap_m,
            activation_tolerance: row.activation_tolerance_m,
            relative_normal_speed: speed,
            restitution: row.restitution,
        })
        .collect();
    let contact = solve_normal_impacts(
        &native,
        &mobility.matrix,
        controls.contact_solve,
        &mut *cancelled,
    )
    .map_err(ResponseError::Contact)?;
    let body_impulses = mechanics::scatter(bodies, &rows, &contact.normal_impulses, cancelled)?;
    poll(cancelled)?;
    Ok(Plan {
        body_impulses,
        mobility,
        contact,
        contact_ids: contacts.iter().map(|row| row.id).collect(),
    })
}

/// Audit original physical successors. The enclosing owner publishes no
/// candidate until every body and the complete event account are admitted.
pub(super) fn finish(
    plan: Plan,
    contacts: &[NormalContact<'_>],
    final_velocities: &[&[f64]],
    kinetic_changes_j: &[f64],
    body_work_residuals_j: &[f64],
    controls: ResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResponseReport, ResponseError> {
    poll(cancelled)?;
    if contacts.len() != plan.mobility.before.len()
        || contacts
            .iter()
            .zip(&plan.contact_ids)
            .any(|(row, &id)| row.id != id)
    {
        return Err(ResponseError::Invalid(
            "original body/contact audit ordering required",
        ));
    }
    let audit = mechanics::audit(
        &rows(contacts),
        &plan.body_impulses,
        &plan.mobility.before,
        &plan.contact.relative_normal_speeds_after,
        &plan.contact.normal_impulses,
        plan.contact.dissipated_energy,
        final_velocities,
        kinetic_changes_j,
        body_work_residuals_j,
        controls.velocity_tolerance_m_s,
        controls.energy_tolerance_j,
        cancelled,
    )?;
    Ok(ResponseReport {
        contact_ids: plan.contact_ids,
        contact_mobility_per_kg: plan.mobility.matrix,
        mass_solves: plan.mobility.mass_solves,
        contact: plan.contact,
        relative_velocities_before_m_s: plan.mobility.before,
        relative_velocities_after_m_s: audit.after,
        kinetic_change_j: audit.kinetic_change_j,
        prescribed_work_j: audit.prescribed_work_j,
        relative_impulse_work_j: audit.relative_impulse_work_j,
        work_energy_residual_j: audit.work_energy_residual_j,
        dissipation_residual_j: audit.dissipation_residual_j,
        maximum_velocity_residual_m_s: audit.maximum_velocity_residual_m_s,
        absolute_body_work_residual_j: audit.absolute_body_work_residual_j,
    })
}

fn rows<'t>(contacts: &[NormalContact<'t>]) -> Vec<TraceRow<'t>> {
    contacts
        .iter()
        .enumerate()
        .map(|(contact, row)| TraceRow {
            contact,
            component: 0,
            terms: row.trace,
            prescribed_velocity: row.prescribed_normal_velocity_m_s,
        })
        .collect()
}

fn admit_contacts(
    contacts: &[NormalContact<'_>],
    controls: ResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, ResponseError> {
    poll(cancelled)?;
    let config = controls.contact_solve.solver;
    if contacts.len() > 64
        || contacts.len() > config.max_variables
        || !(1..=64).contains(&config.max_variables)
        || !(1..=4096).contains(&config.max_iterations)
        || !positive(config.gradient_tolerance)
        || !positive(config.symmetry_tolerance)
        || config.symmetry_tolerance > 1e-6
        || !(64.0 * f64::EPSILON..=1e-4).contains(&config.relative_pivot_tolerance)
        || !nonnegative(controls.contact_solve.energy_tolerance)
        || !nonnegative(controls.velocity_tolerance_m_s)
        || !nonnegative(controls.energy_tolerance_j)
    {
        return Err(ResponseError::Invalid(
            "bounded dimensions and finite numerical controls required",
        ));
    }
    for (i, row) in contacts.iter().enumerate() {
        poll(cancelled)?;
        if contacts[..i].iter().any(|previous| previous.id == row.id)
            || !nonnegative(row.gap_m)
            || !nonnegative(row.activation_tolerance_m)
            || !row.restitution.is_finite()
            || !(0.0..=1.0).contains(&row.restitution)
            || !row.prescribed_normal_velocity_m_s.is_finite()
        {
            return Err(ResponseError::Invalid(
                "unique IDs and finite nonpenetrating contact data required",
            ));
        }
    }
    let n = contacts.len();
    let mut workspace = add(mul(mul(n, n)?, 32)?, mul(n, 64)?)?;
    workspace = add(workspace, mechanics::packed_factor_workspace(n))?;
    if workspace > controls.maximum_workspace_entries {
        return Err(ResponseError::Budget("complete impact numerical workspace"));
    }
    Ok(workspace)
}

#[cfg(test)]
#[path = "response/tests.rs"]
mod tests;

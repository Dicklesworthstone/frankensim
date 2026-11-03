//! Complete coupled Coulomb response through the original reduced inertia.
//!
//! Caller-owned normal-first orthonormal frames define three conjugate sparse
//! rows per contact. One native solve retains every cross-axis and cross-contact
//! mobility entry; no successive normal solve or tangential clipping is used.

use fs_contact::frictional_multicontact::{
    FrictionalContact as NativeContact, FrictionalContactConfig, FrictionalContactResult,
    NormalVelocityLaw, solve_frictional_contacts,
};
use fs_solver::checked_cg::CheckedCgConfig;

pub(super) use super::super::mechanics::Body;
use super::super::mechanics::{
    self, Mobility, MobilityControls, TraceRow, add, mul, nonnegative, poll, positive,
};
pub use super::super::mechanics::{MassSolveReport, ResponseError, TraceTerm};

/// Frozen contact geometry in a caller-owned orthonormal normal-first frame.
/// Positive normal velocity separates the sides. Each sparse row is strictly
/// sorted by (body, dof), with no duplicate or zero coefficients. The trace
/// transpose distributes work-conjugate impulses to the original body DOFs.
/// Cartesian force/torque claims require the caller's actual geometric trace.
#[derive(Clone, Copy, Debug)]
pub struct FrictionalContact<'t> {
    pub id: u64,
    pub gap_m: f64,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
    /// Prescribed contribution c in relative_velocity = B v + c [m/s],
    /// in normal, tangent 1, tangent 2 order. Its body work is -c dot p.
    pub prescribed_velocity_m_s: [f64; 3],
    pub trace: [&'t [TraceTerm]; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct FrictionalResponseControls {
    pub mass_solve: CheckedCgConfig,
    pub contact_solve: FrictionalContactConfig,
    pub maximum_bodies: usize,
    pub maximum_trace_terms: usize,
    /// Conservative scalar-equivalent numerical workspace ceiling, including
    /// body tables, native CG/Newton/FGMRES histories, dense contact admission,
    /// and packed factorization scratch. Borrowed operators/preconditioners
    /// and accepted-state clones belong to their existing physical owner.
    pub maximum_workspace_entries: usize,
    pub velocity_tolerance_m_s: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug)]
pub struct FrictionalResponseReport {
    pub contact_ids: Vec<u64>,
    /// Full row-major 3m by 3m mobility [1/kg], normal-first within each contact.
    /// Off-diagonal normal/tangent and cross-contact responses remain intact.
    pub contact_mobility_per_kg: Vec<f64>,
    pub mass_solves: Vec<MassSolveReport>,
    /// Original law/determinacy/work receipt. A velocity-only result retains
    /// representative-only impulses; it does not assign individual heat/loads.
    pub contact: FrictionalContactResult,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    /// Re-sampled from actual returned body velocities, including prescribed c.
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    pub kinetic_change_j: f64,
    pub prescribed_work_j: f64,
    pub relative_impulse_work_j: f64,
    /// sum Delta K - prescribed_work - actual relative impulse work.
    pub work_energy_residual_j: f64,
    /// sum Delta K - prescribed_work + native impact loss. Native positive
    /// work allowances remain visible in this physical conservation residual.
    pub dissipation_residual_j: f64,
    pub maximum_velocity_residual_m_s: f64,
    /// Sum of absolute native body work defects, before any cancellation.
    pub absolute_body_work_residual_j: f64,
}

pub(super) struct Plan {
    pub body_impulses: Vec<Vec<f64>>,
    contact_ids: Vec<u64>,
    mobility: Mobility,
    contact: FrictionalContactResult,
}

pub(super) fn prepare(
    bodies: &[Body<'_>],
    contacts: &[FrictionalContact<'_>],
    controls: FrictionalResponseControls,
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
        .enumerate()
        .map(|(i, row)| NativeContact {
            id: row.id,
            gap: row.gap_m,
            activation_tolerance: row.activation_tolerance_m,
            relative_velocity: [
                mobility.before[3 * i],
                mobility.before[3 * i + 1],
                mobility.before[3 * i + 2],
            ],
            friction_coefficient: row.friction_coefficient,
            normal_law: row.normal_law,
        })
        .collect();
    let contact = solve_frictional_contacts(
        &native,
        &mobility.matrix,
        controls.contact_solve,
        &mut *cancelled,
    )
    .map_err(ResponseError::FrictionalContact)?;
    let impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let body_impulses = mechanics::scatter(bodies, &rows, &impulses, cancelled)?;
    poll(cancelled)?;
    Ok(Plan {
        body_impulses,
        mobility,
        contact,
        contact_ids: contacts.iter().map(|row| row.id).collect(),
    })
}

/// Admit only the actual physical successors produced by the original body
/// impulse owner. Failure leaves the enclosing owner with no candidate.
pub(super) fn finish(
    plan: Plan,
    contacts: &[FrictionalContact<'_>],
    final_velocities: &[&[f64]],
    kinetic_changes_j: &[f64],
    body_work_residuals_j: &[f64],
    controls: FrictionalResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FrictionalResponseReport, ResponseError> {
    poll(cancelled)?;
    if contacts.len() != plan.contact_ids.len()
        || contacts
            .iter()
            .zip(&plan.contact_ids)
            .any(|(row, &id)| row.id != id)
    {
        return Err(ResponseError::Invalid(
            "original body/contact audit ordering required",
        ));
    }
    let impulses: Vec<_> = plan.contact.impulses.iter().flatten().copied().collect();
    let predicted: Vec<_> = plan
        .contact
        .relative_velocities_after
        .iter()
        .flatten()
        .copied()
        .collect();
    let audit = mechanics::audit(
        &rows(contacts),
        &plan.body_impulses,
        &plan.mobility.before,
        &predicted,
        &impulses,
        plan.contact.dissipated_energy,
        final_velocities,
        kinetic_changes_j,
        body_work_residuals_j,
        controls.velocity_tolerance_m_s,
        controls.energy_tolerance_j,
        cancelled,
    )?;
    Ok(FrictionalResponseReport {
        contact_ids: plan.contact_ids,
        contact_mobility_per_kg: plan.mobility.matrix,
        mass_solves: plan.mobility.mass_solves,
        contact: plan.contact,
        relative_velocities_before_m_s: triples(&plan.mobility.before),
        relative_velocities_after_m_s: triples(&audit.after),
        kinetic_change_j: audit.kinetic_change_j,
        prescribed_work_j: audit.prescribed_work_j,
        relative_impulse_work_j: audit.relative_impulse_work_j,
        work_energy_residual_j: audit.work_energy_residual_j,
        dissipation_residual_j: audit.dissipation_residual_j,
        maximum_velocity_residual_m_s: audit.maximum_velocity_residual_m_s,
        absolute_body_work_residual_j: audit.absolute_body_work_residual_j,
    })
}

fn rows<'t>(contacts: &[FrictionalContact<'t>]) -> Vec<TraceRow<'t>> {
    contacts
        .iter()
        .enumerate()
        .flat_map(|(contact, row)| {
            (0..3).map(move |component| TraceRow {
                contact,
                component,
                terms: row.trace[component],
                prescribed_velocity: row.prescribed_velocity_m_s[component],
            })
        })
        .collect()
}

fn triples(values: &[f64]) -> Vec<[f64; 3]> {
    values.chunks_exact(3).map(|v| [v[0], v[1], v[2]]).collect()
}

fn admit_contacts(
    contacts: &[FrictionalContact<'_>],
    controls: FrictionalResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, ResponseError> {
    poll(cancelled)?;
    let config = controls.contact_solve;
    let matrix = config.matrix;
    if !(1..=21).contains(&config.maximum_contacts)
        || contacts.len() > config.maximum_contacts
        || !(1..=1_000_000).contains(&config.maximum_evaluations)
        || !(1..=256).contains(&config.maximum_newton_iterations)
        || !positive(config.velocity_tolerance)
        || !positive(config.impulse_tolerance)
        || !positive(config.complementarity_tolerance)
        || !nonnegative(config.energy_tolerance)
        || config.newton.linear_restart > 64
        || config.newton.max_linear_cycles > 128
        || !(1..=64).contains(&matrix.max_variables)
        || !(1..=4096).contains(&matrix.max_iterations)
        || 3 * contacts.len() > matrix.max_variables
        || !positive(matrix.gradient_tolerance)
        || !positive(matrix.symmetry_tolerance)
        || matrix.symmetry_tolerance > 1e-6
        || !(64.0 * f64::EPSILON..=1e-4).contains(&matrix.relative_pivot_tolerance)
        || !nonnegative(controls.velocity_tolerance_m_s)
        || !nonnegative(controls.energy_tolerance_j)
    {
        return Err(ResponseError::Invalid(
            "bounded dimensions and finite Coulomb controls required",
        ));
    }
    config.newton.validate().map_err(|error| {
        ResponseError::FrictionalContact(
            fs_contact::frictional_multicontact::FrictionalContactError::Newton(error),
        )
    })?;
    for (i, row) in contacts.iter().enumerate() {
        poll(cancelled)?;
        let law_valid = match row.normal_law {
            NormalVelocityLaw::NewtonRestitution { coefficient } => {
                coefficient.is_finite() && (0.0..=1.0).contains(&coefficient)
            }
            NormalVelocityLaw::Target { speed } => nonnegative(speed),
        };
        if contacts[..i].iter().any(|previous| previous.id == row.id)
            || !nonnegative(row.gap_m)
            || !nonnegative(row.activation_tolerance_m)
            || !nonnegative(row.friction_coefficient)
            || !law_valid
            || row.prescribed_velocity_m_s.iter().any(|v| !v.is_finite())
        {
            return Err(ResponseError::Invalid(
                "unique IDs and finite nonpenetrating Coulomb data required",
            ));
        }
    }
    let d = mul(contacts.len(), 3)?;
    // Account for retained/temporary dense range, SVD, energy and Newton
    // matrices; packed Cholesky scratch is separate. FGMRES restart is the
    // caller's actual restart, which its native owner does NOT clamp to d.
    let restart = config.newton.linear_restart;
    let mut workspace = add(mul(mul(d, d)?, 96)?, mul(d, 128)?)?;
    for entries in [
        mul(mul(add(restart, 1)?, d)?, 2)?,
        mul(add(restart, 1)?, restart)?,
        mul(restart, 4)?,
        mul(config.maximum_newton_iterations, 64)?,
        mul(config.newton.max_linear_cycles, 8)?,
        mechanics::packed_factor_workspace(d),
    ] {
        workspace = add(workspace, entries)?;
    }
    if workspace > controls.maximum_workspace_entries {
        return Err(ResponseError::Budget(
            "complete Coulomb numerical workspace",
        ));
    }
    Ok(workspace)
}

#[cfg(test)]
#[path = "response/tests.rs"]
mod tests;

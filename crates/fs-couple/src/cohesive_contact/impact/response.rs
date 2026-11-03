//! Original reduced inertia and one simultaneous normal-contact law.
//!
//! Signed traces pair velocity sampling with transpose impulse distribution.
//! Every cross-contact mobility entry is assembled from the retained mass; no
//! diagonal approximation, matrix averaging or pairwise impact sweep is used.
//! Sparse geometry and SPD mass/preconditioner admission belong to their owners.

use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::{LinearOp, SolveReport};
use fs_sparse::precond::Precond;

/// One signed dimensionless entry in a work-conjugate normal-velocity trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceTerm {
    pub body: usize,
    pub dof: usize,
    pub coefficient: f64,
}

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

/// Original checked-CG receipt for one nonzero body/contact mobility column.
#[derive(Clone, Debug)]
pub struct MassSolveReport {
    pub body: usize,
    /// Index into the caller's contact order, also used by the report arrays.
    pub contact: usize,
    pub solve: SolveReport,
    pub corrections: usize,
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

#[derive(Debug)]
pub enum ResponseError {
    Invalid(&'static str),
    Budget(&'static str),
    NonFinite(&'static str),
    Cancelled,
    MassSolve {
        body: usize,
        contact: usize,
        report: SolveReport,
    },
    Contact(MultiImpactError),
    Audit {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for ResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive normal impact refused: {self:?}")
    }
}
impl std::error::Error for ResponseError {}

pub(super) struct Body<'b> {
    pub mass: &'b dyn LinearOp,
    pub velocity: &'b [f64],
    pub preconditioner: &'b dyn Precond,
}

pub(super) struct Plan {
    pub body_impulses: Vec<Vec<f64>>,
    contact_ids: Vec<u64>,
    mobility: Vec<f64>,
    mass_solves: Vec<MassSolveReport>,
    before: Vec<f64>,
    contact: MultiImpactResult,
}

struct Mass<'a>(&'a dyn LinearOp);
impl LinearOp for Mass<'_> {
    fn n(&self) -> usize {
        self.0.n()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply_transpose(x, y);
    }
}
struct Preconditioner<'a>(&'a dyn Precond);
impl Precond for Preconditioner<'_> {
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.apply(x, y);
    }
}

pub(super) fn prepare(
    bodies: &[Body<'_>],
    contacts: &[NormalContact<'_>],
    controls: ResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Plan, ResponseError> {
    let solves = admit(bodies, contacts, controls, cancelled)?;
    let m = contacts.len();
    let mut mobility = vec![Sum::default(); m * m];
    let mut mass_solves = Vec::with_capacity(solves);
    for (body_index, body) in bodies.iter().enumerate() {
        poll(cancelled)?;
        if contacts
            .iter()
            .all(|row| body_terms(row, body_index).is_empty())
        {
            continue;
        }
        let mut rhs = vec![0.0; body.velocity.len()];
        for (j, contact) in contacts.iter().enumerate() {
            let terms = body_terms(contact, body_index);
            if terms.is_empty() {
                continue;
            }
            poll(cancelled)?;
            for chunk in rhs.chunks_mut(256) {
                poll(cancelled)?;
                chunk.fill(0.0);
            }
            for (k, term) in terms.iter().enumerate() {
                if k % 256 == 0 {
                    poll(cancelled)?;
                }
                rhs[term.dof] = term.coefficient;
            }
            let solution = checked_cg(
                &Mass(body.mass),
                &Preconditioner(body.preconditioner),
                &rhs,
                controls.mass_solve,
                |_| poll(cancelled),
            )
            .map_err(|error| match error {
                CheckedCgError::InvalidInput(reason) => ResponseError::Invalid(reason),
                CheckedCgError::Interrupted(error) => error,
            })?;
            if !solution.report.converged_euclidean() {
                return Err(ResponseError::MassSolve {
                    body: body_index,
                    contact: j,
                    report: solution.report,
                });
            }
            for (i, row) in contacts.iter().enumerate() {
                poll(cancelled)?;
                let mut entry = Sum::default();
                for (k, term) in body_terms(row, body_index).iter().enumerate() {
                    if k % 256 == 0 {
                        poll(cancelled)?;
                    }
                    entry.add(term.coefficient * solution.x[term.dof])?;
                }
                mobility[i * m + j].add(entry.total()?)?;
            }
            mass_solves.push(MassSolveReport {
                body: body_index,
                contact: j,
                solve: solution.report,
                corrections: solution.corrections,
            });
        }
    }
    let mobility = mobility
        .into_iter()
        .map(Sum::total)
        .collect::<Result<Vec<_>, _>>()?;
    let mut before = Vec::with_capacity(m);
    let mut native = Vec::with_capacity(m);
    for contact in contacts {
        let speed = sample(
            contact,
            |term| bodies[term.body].velocity[term.dof],
            cancelled,
        )?;
        before.push(speed);
        native.push(NormalContactImpact {
            gap: contact.gap_m,
            activation_tolerance: contact.activation_tolerance_m,
            relative_normal_speed: speed,
            restitution: contact.restitution,
        });
    }
    let contact = solve_normal_impacts(&native, &mobility, controls.contact_solve, &mut *cancelled)
        .map_err(ResponseError::Contact)?;
    let mut impulses: Vec<Vec<Sum>> = bodies
        .iter()
        .map(|body| vec![Sum::default(); body.velocity.len()])
        .collect();
    for (row, &impulse) in contacts.iter().zip(&contact.normal_impulses) {
        for (k, term) in row.trace.iter().enumerate() {
            if k % 256 == 0 {
                poll(cancelled)?;
            }
            impulses[term.body][term.dof].add(term.coefficient * impulse)?;
        }
    }
    let mut body_impulses = Vec::with_capacity(bodies.len());
    for body in impulses {
        poll(cancelled)?;
        let mut values = Vec::with_capacity(body.len());
        for (i, value) in body.into_iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            values.push(value.total()?);
        }
        body_impulses.push(values);
    }
    poll(cancelled)?;
    Ok(Plan {
        body_impulses,
        mobility,
        mass_solves,
        before,
        contact,
        contact_ids: contacts.iter().map(|row| row.id).collect(),
    })
}

/// Audit original physical successors. All arrays keep the original body/contact
/// order; the enclosing owner publishes no candidate until this succeeds.
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
    let n = plan.body_impulses.len();
    if final_velocities.len() != n
        || kinetic_changes_j.len() != n
        || body_work_residuals_j.len() != n
        || contacts.len() != plan.before.len()
        || contacts
            .iter()
            .zip(&plan.contact_ids)
            .any(|(row, &id)| row.id != id)
    {
        return Err(ResponseError::Invalid(
            "original body/contact audit ordering required",
        ));
    }
    let mut kinetic = Sum::default();
    let mut absolute_body_defect = Sum::default();
    let mut closure = Sum::default();
    let mut dissipation = Sum::default();
    for (b, &velocity) in final_velocities.iter().enumerate() {
        if velocity.len() != plan.body_impulses[b].len() {
            return Err(ResponseError::Invalid("returned body velocity dimension"));
        }
        for (i, &value) in velocity.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            finite(value)?;
        }
        kinetic.add(kinetic_changes_j[b])?;
        closure.add(kinetic_changes_j[b])?;
        dissipation.add(kinetic_changes_j[b])?;
        absolute_body_defect.add(body_work_residuals_j[b].abs())?;
    }
    let mut after = Vec::with_capacity(contacts.len());
    let mut relative_work = Sum::default();
    let mut prescribed_work = Sum::default();
    let mut velocity_defect = 0.0_f64;
    for (i, row) in contacts.iter().enumerate() {
        let speed = sample(row, |term| final_velocities[term.body][term.dof], cancelled)?;
        let residual = finite(speed - plan.contact.relative_normal_speeds_after[i])?;
        velocity_defect = velocity_defect.max(residual.abs());
        let impulse = plan.contact.normal_impulses[i];
        let work = finite(f64::midpoint(plan.before[i], speed) * impulse)?;
        let supplied = finite(-row.prescribed_normal_velocity_m_s * impulse)?;
        relative_work.add(work)?;
        prescribed_work.add(supplied)?;
        closure.add(-supplied)?;
        closure.add(-work)?;
        dissipation.add(-supplied)?;
        after.push(speed);
    }
    let relative_work = relative_work.total()?;
    dissipation.add(plan.contact.dissipated_energy)?;
    let report = ResponseReport {
        contact_ids: plan.contact_ids,
        contact_mobility_per_kg: plan.mobility,
        mass_solves: plan.mass_solves,
        contact: plan.contact,
        relative_velocities_before_m_s: plan.before,
        relative_velocities_after_m_s: after,
        kinetic_change_j: kinetic.total()?,
        prescribed_work_j: prescribed_work.total()?,
        relative_impulse_work_j: relative_work,
        work_energy_residual_j: closure.total()?,
        dissipation_residual_j: dissipation.total()?,
        maximum_velocity_residual_m_s: velocity_defect,
        absolute_body_work_residual_j: absolute_body_defect.total()?,
    };
    for (quantity, residual, tolerance) in [
        (
            "native endpoint velocity [m/s]",
            report.maximum_velocity_residual_m_s,
            controls.velocity_tolerance_m_s,
        ),
        (
            "absolute native body work [J]",
            report.absolute_body_work_residual_j,
            controls.energy_tolerance_j,
        ),
        (
            "assembly work closure [J]",
            report.work_energy_residual_j.abs(),
            controls.energy_tolerance_j,
        ),
        (
            "native contact dissipation [J]",
            report.dissipation_residual_j.abs(),
            controls.energy_tolerance_j,
        ),
    ] {
        if residual > tolerance {
            return Err(ResponseError::Audit {
                quantity,
                residual,
                tolerance,
            });
        }
    }
    poll(cancelled)?;
    Ok(report)
}

fn admit(
    bodies: &[Body<'_>],
    contacts: &[NormalContact<'_>],
    controls: ResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, ResponseError> {
    poll(cancelled)?;
    let config = controls.contact_solve.solver;
    if bodies.is_empty()
        || bodies.len() > controls.maximum_bodies
        || contacts.len() > 64
        || contacts.len() > config.max_variables
        || !(1..=64).contains(&config.max_variables)
        || !(1..=4096).contains(&config.max_iterations)
        || !positive(config.gradient_tolerance)
        || !positive(config.symmetry_tolerance)
        || config.symmetry_tolerance > 1e-6
        || !(64.0 * f64::EPSILON..=1e-4).contains(&config.relative_pivot_tolerance)
        || !controls.mass_solve.tolerance.is_finite()
        || controls.mass_solve.tolerance <= 0.0
        || controls.mass_solve.tolerance >= 1.0
        || !nonnegative(controls.contact_solve.energy_tolerance)
        || !nonnegative(controls.velocity_tolerance_m_s)
        || !nonnegative(controls.energy_tolerance_j)
    {
        return Err(ResponseError::Invalid(
            "bounded dimensions and finite numerical controls required",
        ));
    }
    let mut total_dofs = 0usize;
    let mut largest_body = 0usize;
    for body in bodies {
        poll(cancelled)?;
        if body.velocity.is_empty() || body.mass.n() != body.velocity.len() {
            return Err(ResponseError::Invalid(
                "original nonempty mass/velocity dimensions required",
            ));
        }
        total_dofs = add(total_dofs, body.velocity.len())?;
        largest_body = largest_body.max(body.velocity.len());
        for (i, &v) in body.velocity.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            finite(v)?;
        }
    }
    let mut trace_count = 0usize;
    let mut solves = 0usize;
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
        trace_count = add(trace_count, row.trace.len())?;
        if trace_count > controls.maximum_trace_terms {
            return Err(ResponseError::Budget("contact trace terms"));
        }
        let mut previous = None;
        for (k, term) in row.trace.iter().enumerate() {
            if k % 256 == 0 {
                poll(cancelled)?;
            }
            let key = (term.body, term.dof);
            if term.body >= bodies.len()
                || term.dof >= bodies[term.body].velocity.len()
                || !term.coefficient.is_finite()
                || term.coefficient == 0.0
                || previous.is_some_and(|old| old >= key)
            {
                return Err(ResponseError::Invalid(
                    "strictly sorted unique finite nonzero trace entries required",
                ));
            }
            if previous.is_none_or(|(body, _)| body != term.body) {
                solves = add(solves, 1)?;
            }
            previous = Some(key);
        }
    }
    // Sixteen largest-body vectors cover the original checked CG and its RHS.
    // Returned reports retain only outer-correction histories; native inner CG
    // also grows and clones its history transiently. The 32 contact matrices
    // reserve the bounded native SPD/PSD spectral/Gram solve, not dense masses.
    let iterations = controls.mass_solve.max_iterations;
    let corrections = controls.mass_solve.max_corrections.min(iterations);
    let retained = mul(solves, add(mul(corrections, 2)?, 24)?)?;
    let mut workspace = add(mul(largest_body, 16)?, mul(total_dofs, 3)?)?;
    for entries in [
        retained,
        mul(iterations, 4)?,
        mul(corrections, 2)?,
        mul(mul(contacts.len(), contacts.len())?, 32)?,
        mul(contacts.len(), 64)?,
        mul(bodies.len(), 24)?,
        64,
    ] {
        workspace = add(workspace, entries)?;
    }
    if contacts.len() > 32 {
        // The original blocked Cholesky uses the original serial GEMM's
        // MC*KC + KC*NC packed buffers, even for a small trailing block.
        workspace = add(workspace, 128 * 256 + 256 * 512)?;
    }
    if workspace > controls.maximum_workspace_entries {
        return Err(ResponseError::Budget("complete impact numerical workspace"));
    }
    Ok(solves)
}

fn body_terms<'r>(row: &'r NormalContact<'_>, body: usize) -> &'r [TraceTerm] {
    let start = row.trace.partition_point(|term| term.body < body);
    let end = row.trace.partition_point(|term| term.body <= body);
    &row.trace[start..end]
}

fn sample(
    row: &NormalContact<'_>,
    value: impl Fn(&TraceTerm) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ResponseError> {
    poll(cancelled)?;
    let mut sum = Sum::default();
    sum.add(row.prescribed_normal_velocity_m_s)?;
    for (k, term) in row.trace.iter().enumerate() {
        if k % 256 == 0 {
            poll(cancelled)?;
        }
        sum.add(term.coefficient * value(term))?;
    }
    sum.total()
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), ResponseError> {
        finite(value)?;
        let next = finite(self.value + value)?;
        let residual = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + residual)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, ResponseError> {
        finite(self.value + self.correction)
    }
}
fn finite(value: f64) -> Result<f64, ResponseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ResponseError::NonFinite("derived impact value"))
    }
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn add(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_add(b)
        .ok_or(ResponseError::Budget("workspace/count overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_mul(b)
        .ok_or(ResponseError::Budget("workspace/count overflow"))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResponseError> {
    if cancelled() {
        Err(ResponseError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "response/tests.rs"]
mod tests;

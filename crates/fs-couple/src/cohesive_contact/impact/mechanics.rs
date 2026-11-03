//! Shared mechanical response for normal and Coulomb impact laws.
//!
//! The physical law chooses contact impulses. This module supplies only the
//! original mass mobility, conjugate distribution, and actual-state work audit.

use fs_contact::frictional_multicontact::FrictionalContactError;
use fs_contact::multicontact::MultiImpactError;
use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::{LinearOp, SolveReport};
use fs_sparse::precond::Precond;

/// One signed dimensionless entry in a work-conjugate velocity trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceTerm {
    pub body: usize,
    pub dof: usize,
    pub coefficient: f64,
}

/// Original checked-CG receipt for one nonzero body/contact mobility column.
#[derive(Clone, Debug)]
pub struct MassSolveReport {
    pub body: usize,
    /// Index into the caller's original contact order.
    pub contact: usize,
    /// Zero for normal contact; normal/tangent 1/tangent 2 for Coulomb contact.
    pub component: usize,
    pub solve: SolveReport,
    pub corrections: usize,
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
    FrictionalContact(FrictionalContactError),
    Audit {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for ResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive impact response refused: {self:?}")
    }
}
impl std::error::Error for ResponseError {}

pub(super) struct Body<'b> {
    pub mass: &'b dyn LinearOp,
    pub velocity: &'b [f64],
    pub preconditioner: &'b dyn Precond,
}

/// A mechanical coordinate, with no assumed normal or friction law.
#[derive(Clone, Copy)]
pub(super) struct TraceRow<'t> {
    pub contact: usize,
    pub component: usize,
    pub terms: &'t [TraceTerm],
    pub prescribed_velocity: f64,
}

#[derive(Clone, Copy)]
pub(super) struct MobilityControls {
    pub mass_solve: CheckedCgConfig,
    pub maximum_bodies: usize,
    pub maximum_trace_terms: usize,
    pub maximum_workspace_entries: usize,
    /// Additional scratch retained by the actual contact-law numerical owner.
    pub contact_workspace_entries: usize,
}

pub(super) struct Mobility {
    pub matrix: Vec<f64>,
    pub mass_solves: Vec<MassSolveReport>,
    pub before: Vec<f64>,
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

/// Assemble all entries without diagonalizing or averaging the response.
/// There is one native mass solve per nonzero body/coordinate column.
pub(super) fn assemble(
    bodies: &[Body<'_>],
    rows: &[TraceRow<'_>],
    controls: MobilityControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Mobility, ResponseError> {
    let solves = admit(bodies, rows, controls, cancelled)?;
    let n = rows.len();
    let mut matrix = vec![Sum::default(); n * n];
    let mut mass_solves = Vec::with_capacity(solves);
    for (body_index, body) in bodies.iter().enumerate() {
        poll(cancelled)?;
        if rows
            .iter()
            .all(|row| body_terms(row, body_index).is_empty())
        {
            continue;
        }
        let mut rhs = vec![0.0; body.velocity.len()];
        for (column, row) in rows.iter().enumerate() {
            let terms = body_terms(row, body_index);
            if terms.is_empty() {
                continue;
            }
            for chunk in rhs.chunks_mut(256) {
                poll(cancelled)?;
                chunk.fill(0.0);
            }
            for (i, term) in terms.iter().enumerate() {
                if i % 256 == 0 {
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
                    contact: row.contact,
                    report: solution.report,
                });
            }
            for (i, sampled) in rows.iter().enumerate() {
                poll(cancelled)?;
                let mut entry = Sum::default();
                for (k, term) in body_terms(sampled, body_index).iter().enumerate() {
                    if k % 256 == 0 {
                        poll(cancelled)?;
                    }
                    entry.add(term.coefficient * solution.x[term.dof])?;
                }
                matrix[i * n + column].add(entry.total()?)?;
            }
            mass_solves.push(MassSolveReport {
                body: body_index,
                contact: row.contact,
                component: row.component,
                solve: solution.report,
                corrections: solution.corrections,
            });
        }
    }
    let matrix = matrix
        .into_iter()
        .map(Sum::total)
        .collect::<Result<Vec<_>, _>>()?;
    let mut before = Vec::with_capacity(rows.len());
    for row in rows {
        before.push(sample(
            row,
            |term| bodies[term.body].velocity[term.dof],
            cancelled,
        )?);
    }
    poll(cancelled)?;
    Ok(Mobility {
        matrix,
        mass_solves,
        before,
    })
}

pub(super) fn scatter(
    bodies: &[Body<'_>],
    rows: &[TraceRow<'_>],
    impulses: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Vec<f64>>, ResponseError> {
    if impulses.len() != rows.len() {
        return Err(ResponseError::Invalid(
            "complete contact impulse coordinates required",
        ));
    }
    let mut sums: Vec<Vec<Sum>> = bodies
        .iter()
        .map(|body| vec![Sum::default(); body.velocity.len()])
        .collect();
    for (row, &impulse) in rows.iter().zip(impulses) {
        finite(impulse)?;
        for (i, term) in row.terms.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            sums[term.body][term.dof].add(term.coefficient * impulse)?;
        }
    }
    let mut result = Vec::with_capacity(bodies.len());
    for body in sums {
        poll(cancelled)?;
        let mut values = Vec::with_capacity(body.len());
        for (i, value) in body.into_iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            values.push(value.total()?);
        }
        result.push(values);
    }
    poll(cancelled)?;
    Ok(result)
}

pub(super) struct WorkAudit {
    pub after: Vec<f64>,
    pub kinetic_change_j: f64,
    pub prescribed_work_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub dissipation_residual_j: f64,
    pub maximum_velocity_residual_m_s: f64,
    pub absolute_body_work_residual_j: f64,
}

/// Audit physical successors through the original velocity coordinates.
/// The total dissipation account retains any allowed native positive work.
#[allow(clippy::too_many_arguments)]
pub(super) fn audit(
    rows: &[TraceRow<'_>],
    body_impulses: &[Vec<f64>],
    before: &[f64],
    predicted: &[f64],
    impulses: &[f64],
    native_loss_j: f64,
    final_velocities: &[&[f64]],
    kinetic_changes_j: &[f64],
    body_work_residuals_j: &[f64],
    velocity_tolerance_m_s: f64,
    energy_tolerance_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WorkAudit, ResponseError> {
    poll(cancelled)?;
    let n = body_impulses.len();
    if final_velocities.len() != n
        || kinetic_changes_j.len() != n
        || body_work_residuals_j.len() != n
        || before.len() != rows.len()
        || predicted.len() != rows.len()
        || impulses.len() != rows.len()
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
        if velocity.len() != body_impulses[b].len() {
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
    let mut after = Vec::with_capacity(rows.len());
    let mut relative_work = Sum::default();
    let mut prescribed_work = Sum::default();
    let mut velocity_defect = 0.0_f64;
    for (i, row) in rows.iter().enumerate() {
        let speed = sample(row, |term| final_velocities[term.body][term.dof], cancelled)?;
        let residual = finite(speed - predicted[i])?;
        velocity_defect = velocity_defect.max(residual.abs());
        let work = finite(f64::midpoint(before[i], speed) * impulses[i])?;
        let supplied = finite(-row.prescribed_velocity * impulses[i])?;
        relative_work.add(work)?;
        prescribed_work.add(supplied)?;
        closure.add(-supplied)?;
        closure.add(-work)?;
        dissipation.add(-supplied)?;
        after.push(speed);
    }
    dissipation.add(native_loss_j)?;
    let report = WorkAudit {
        after,
        kinetic_change_j: kinetic.total()?,
        prescribed_work_j: prescribed_work.total()?,
        relative_impulse_work_j: relative_work.total()?,
        work_energy_residual_j: closure.total()?,
        dissipation_residual_j: dissipation.total()?,
        maximum_velocity_residual_m_s: velocity_defect,
        absolute_body_work_residual_j: absolute_body_defect.total()?,
    };
    for (quantity, residual, tolerance) in [
        (
            "native endpoint velocity [m/s]",
            report.maximum_velocity_residual_m_s,
            velocity_tolerance_m_s,
        ),
        (
            "absolute native body work [J]",
            report.absolute_body_work_residual_j,
            energy_tolerance_j,
        ),
        (
            "assembly work closure [J]",
            report.work_energy_residual_j.abs(),
            energy_tolerance_j,
        ),
        (
            "native contact dissipation [J]",
            report.dissipation_residual_j.abs(),
            energy_tolerance_j,
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
    rows: &[TraceRow<'_>],
    controls: MobilityControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize, ResponseError> {
    poll(cancelled)?;
    if bodies.is_empty()
        || bodies.len() > controls.maximum_bodies
        || rows.len() > 64
        || !controls.mass_solve.tolerance.is_finite()
        || controls.mass_solve.tolerance <= 0.0
        || controls.mass_solve.tolerance >= 1.0
    {
        return Err(ResponseError::Invalid(
            "bounded dimensions and finite mass controls required",
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
    for row in rows {
        poll(cancelled)?;
        finite(row.prescribed_velocity)?;
        trace_count = add(trace_count, row.terms.len())?;
        if trace_count > controls.maximum_trace_terms {
            return Err(ResponseError::Budget("contact trace terms"));
        }
        let mut previous = None;
        for (i, term) in row.terms.iter().enumerate() {
            if i % 256 == 0 {
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
    // Native mass vectors and both growing/copying inner and retained outer
    // histories. Law-specific factorizations/Newton scratch are added by owner.
    let iterations = controls.mass_solve.max_iterations;
    let corrections = controls.mass_solve.max_corrections.min(iterations);
    let retained = mul(solves, add(mul(corrections, 2)?, 24)?)?;
    let mut workspace = add(mul(largest_body, 16)?, mul(total_dofs, 3)?)?;
    for entries in [
        retained,
        mul(iterations, 4)?,
        mul(corrections, 2)?,
        mul(bodies.len(), 24)?,
        controls.contact_workspace_entries,
        64,
    ] {
        workspace = add(workspace, entries)?;
    }
    if workspace > controls.maximum_workspace_entries {
        return Err(ResponseError::Budget("complete impact numerical workspace"));
    }
    Ok(solves)
}

fn body_terms<'r>(row: &'r TraceRow<'_>, body: usize) -> &'r [TraceTerm] {
    let start = row.terms.partition_point(|term| term.body < body);
    let end = row.terms.partition_point(|term| term.body <= body);
    &row.terms[start..end]
}
fn sample(
    row: &TraceRow<'_>,
    value: impl Fn(&TraceTerm) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ResponseError> {
    poll(cancelled)?;
    let mut sum = Sum::default();
    sum.add(row.prescribed_velocity)?;
    for (i, term) in row.terms.iter().enumerate() {
        if i % 256 == 0 {
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
pub(super) fn finite(value: f64) -> Result<f64, ResponseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ResponseError::NonFinite("derived impact value"))
    }
}
pub(super) fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
pub(super) fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
pub(super) fn add(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_add(b)
        .ok_or(ResponseError::Budget("workspace/count overflow"))
}
pub(super) fn mul(a: usize, b: usize) -> Result<usize, ResponseError> {
    a.checked_mul(b)
        .ok_or(ResponseError::Budget("workspace/count overflow"))
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResponseError> {
    if cancelled() {
        Err(ResponseError::Cancelled)
    } else {
        Ok(())
    }
}

/// Original serial packed GEMM scratch used by blocked Cholesky above 32 rows.
pub(super) fn packed_factor_workspace(rows: usize) -> usize {
    if rows > 32 { 128 * 256 + 256 * 512 } else { 0 }
}

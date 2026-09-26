//! An outward error bound for the *coupled* stored system
//! `(A - B C) x = b + B d`. Neither the interface variables nor their feedback
//! are frozen. A small-gain check uses residual-verified columns of A^-1 B;
//! approximate response solves never become exact inverse authority.
//!
//! If q >= ||A^-1 B C||_inf and q < 1, the Neumann argument gives
//! ||(A-B C)^-1||_inf <= ||A^-1||_inf / (1-q). A failed sufficient check
//! means unknown, NOT that the coupled system is singular or unstable.

use fs_sparse::Csr;
use super::{
    GoalResidualError, GoalResidualLimits, ScalarEnclosure, Work, enclose_goal_error,
    scratch,
};
use super::arithmetic::{add_up, down, mul_up, up};

/// Structural limits for the coupled residual and optional response checks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FeedbackResidualLimits {
    /// Square solid operator and each of its residual passes.
    pub solid: GoalResidualLimits,
    /// Number of interface coordinates.
    pub max_ports: usize,
    /// Combined stored entries of B and C, including explicit zeros.
    pub max_transfer_nonzeros: usize,
    /// Total scalar entries of caller-supplied response columns.
    pub max_response_entries: usize,
    /// Solid row/entry visits across all response checks, including the
    /// initial inverse check. Checked BEFORE repeated sparse traversal.
    pub max_verification_entries: usize,
}

/// Why a whole coupled-system error bound is present or absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeedbackBoundStatus {
    /// The solid inverse and a strict small-gain inequality were checked.
    Enclosed,
    /// The solid operator has no admitted inverse-infinity bound.
    SolidInverseUnavailable,
    /// The computed gain bound is at least one. No ill-posedness is inferred.
    ContractionNotEstablished,
    /// A finite error/inverse/gain bound could not be represented.
    BoundNotRepresentable,
}

/// Immutable evidence for the exact stored A, B, C, b and d, not assembly,
/// continuum, nonlinear-physics or empirical-validation authority.
#[derive(Debug, Clone, PartialEq)]
pub struct FeedbackResidualReport {
    status: FeedbackBoundStatus,
    residual_infinity_upper: f64,
    solid_inverse_infinity_upper: Option<f64>,
    gain_infinity_upper: Option<f64>,
    coupled_inverse_infinity_upper: Option<f64>,
    state_error_infinity_upper: Option<f64>,
    response_residual_infinity_upper: Vec<f64>,
}

impl FeedbackResidualReport {
    /// Whole-system enclosure disposition.
    #[must_use]
    pub const fn status(&self) -> FeedbackBoundStatus { self.status }
    /// ||b+B(d+C x)-A x||_inf, with every evaluation rounded outwards.
    #[must_use]
    pub const fn residual_infinity_upper(&self) -> f64 { self.residual_infinity_upper }
    /// Checked bound for the original solid inverse.
    #[must_use]
    pub const fn solid_inverse_infinity_upper(&self) -> Option<f64> { self.solid_inverse_infinity_upper }
    /// Verified upper bound for ||A^-1 B C||_inf, when representable.
    #[must_use]
    pub const fn gain_infinity_upper(&self) -> Option<f64> { self.gain_infinity_upper }
    /// ||(A-B C)^-1||_inf, only after strict contraction was established.
    #[must_use]
    pub const fn coupled_inverse_infinity_upper(&self) -> Option<f64> { self.coupled_inverse_infinity_upper }
    /// Whole-field error. A regional maximum has this same error allowance
    /// even when its active vertex moves. None is never zero.
    #[must_use]
    pub const fn state_error_infinity_upper(&self) -> Option<f64> { self.state_error_infinity_upper }
    /// One outward ||B[:,j]-A response[j]||_inf per supplied column.
    /// Empty when the norm-only route was used or no solid inverse exists.
    #[must_use]
    pub fn response_residual_infinity_upper(&self) -> &[f64] { &self.response_residual_infinity_upper }
}

fn limit(field: &'static str, required: usize, allowed: usize) -> Result<(), GoalResidualError> {
    if required > allowed { Err(GoalResidualError::Limit { field, required, allowed }) } else { Ok(()) }
}
fn scalar_vector<F: FnMut() -> bool>(n: usize, work: &mut Work<F>) -> Result<Vec<f64>, GoalResidualError> {
    n.checked_mul(std::mem::size_of::<f64>()).ok_or(GoalResidualError::Allocation)?;
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| GoalResidualError::Allocation)?;
    for _ in 0..n { work.tick()?; out.push(0.0); }
    Ok(out)
}
fn vector<F: FnMut() -> bool>(
    field: &'static str, values: &[f64], n: usize, work: &mut Work<F>,
) -> Result<(), GoalResidualError> {
    if values.len() != n { return Err(GoalResidualError::Length { field, expected: n, found: values.len() }); }
    for (index, value) in values.iter().enumerate() {
        work.tick()?;
        if !value.is_finite() { return Err(GoalResidualError::NonFinite { field, index }); }
    }
    Ok(())
}
fn transfer<F: FnMut() -> bool>(
    matrix: &Csr, rows: usize, cols: usize, work: &mut Work<F>,
) -> Result<(), GoalResidualError> {
    if matrix.nrows() != rows || matrix.ncols() != cols {
        return Err(GoalResidualError::Shape { rows: matrix.nrows(), columns: matrix.ncols() });
    }
    for i in 0..rows {
        work.tick()?;
        for &a in matrix.row(i).1 {
            work.tick()?;
            if !a.is_finite() { return Err(GoalResidualError::NonFinite { field: "feedback transfer", index: i }); }
        }
    }
    Ok(())
}

/// Enclose a coupled residual and, when possible, its full solution error.
///
/// B has shape n x p, C has shape p x n, and d has length p. They may have
/// signed coefficients. Optional `responses[j]` approximate A^-1 B[:,j],
/// in *column* order; every supplied column is independently residual-checked.
/// Without responses, a cheaper norm-only perturbation bound is attempted.
/// The better of that bound and the response-based bound is used.
///
/// Storage is O(n+p), excluding borrowed matrices/responses. Response checks
/// cost O(p*(n+nnz(A))) and are explicitly capped before any allocation.
/// No dense n x n coupled matrix or transpose is formed. Checkpoints occur
/// at most 512 visits apart and before publication; cancellation is an error.
///
/// # Errors
/// Malformed/nonfinite data, exceeded structural/work limits, allocation,
/// unrepresentable residual evaluation and cancellation refuse. Missing
/// inverse or contraction evidence remains an explicit non-success report.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn enclose_affine_feedback_error(
    matrix: &Csr, rhs: &[f64], primal: &[f64],
    injection: &Csr, feedback: &Csr, offset: &[f64],
    responses: Option<&[Vec<f64>]>, scaling: Option<&[f64]>,
    limits: FeedbackResidualLimits, checkpoint: impl FnMut() -> bool,
) -> Result<FeedbackResidualReport, GoalResidualError> {
    let mut work = Work { checkpoint, left: 512 };
    work.poll()?;
    let n = matrix.nrows();
    let p = injection.ncols();
    if n == 0 || matrix.ncols() != n || p == 0 {
        return Err(GoalResidualError::Shape { rows: n, columns: matrix.ncols() });
    }
    limit("rows", n, limits.solid.max_rows)?;
    limit("nonzeros", matrix.nnz(), limits.solid.max_nonzeros)?;
    limit("ports", p, limits.max_ports)?;
    let entries = injection.nnz().checked_add(feedback.nnz()).ok_or(GoalResidualError::Allocation)?;
    limit("transfer nonzeros", entries, limits.max_transfer_nonzeros)?;
    let passes = if let Some(columns) = responses {
        if columns.len() != p { return Err(GoalResidualError::Length { field: "response columns", expected: p, found: columns.len() }); }
        let entries = n.checked_mul(p).ok_or(GoalResidualError::Allocation)?;
        limit("response entries", entries, limits.max_response_entries)?;
        p.checked_add(1).ok_or(GoalResidualError::Allocation)?
    } else { 1 };
    let visits = n.checked_add(matrix.nnz()).and_then(|v| v.checked_mul(passes))
        .ok_or(GoalResidualError::Allocation)?;
    limit("verification entries", visits, limits.max_verification_entries)?;
    vector("rhs", rhs, n, &mut work)?;
    vector("primal", primal, n, &mut work)?;
    vector("reference offset", offset, p, &mut work)?;
    transfer(injection, n, p, &mut work)?;
    transfer(feedback, p, n, &mut work)?;
    if let Some(columns) = responses {
        for column in columns { vector("response column", column, n, &mut work)?; }
    }
    let zeros = scalar_vector(n, &mut work)?;
    let solid = enclose_goal_error(matrix, rhs, primal, &zeros, &zeros, scaling,
        limits.solid, || (work.checkpoint)())?;
    let mut ports = scratch(p, &mut work)?;
    let mut c_norm = scalar_vector(p, &mut work)?;
    let mut gain_range = false;
    for j in 0..p {
        work.tick()?;
        ports[j] = ScalarEnclosure::point(offset[j]);
        let (indices, values) = feedback.row(j);
        for (&i, &value) in indices.iter().zip(values) {
            work.tick()?;
            ports[j] = ports[j].add_scaled(ScalarEnclosure::point(primal[i]), value)?;
            match add_up(c_norm[j], value.abs()) {
                Ok(sum) => c_norm[j] = sum,
                Err(_) => gain_range = true,
            }
        }
    }
    let mut residual_inf = 0.0_f64;
    let mut perturbation_inf = 0.0_f64;
    for i in 0..n {
        work.tick()?;
        let mut residual = ScalarEnclosure::point(rhs[i]);
        let (indices, values) = matrix.row(i);
        for (&j, &value) in indices.iter().zip(values) {
            work.tick()?;
            residual = residual.add_scaled(ScalarEnclosure::point(primal[j]), -value)?;
        }
        let (indices, values) = injection.row(i);
        let mut row_norm = 0.0;
        for (&j, &value) in indices.iter().zip(values) {
            work.tick()?;
            residual = residual.add_scaled(ports[j], value)?;
            match mul_up(value.abs(), c_norm[j]).and_then(|term| add_up(row_norm, term)) {
                Ok(sum) => row_norm = sum,
                Err(_) => gain_range = true,
            }
        }
        residual_inf = residual_inf.max(residual.magnitude_upper());
        perturbation_inf = perturbation_inf.max(row_norm);
    }
    let inverse = solid.inverse_infinity_upper();
    let mut report = FeedbackResidualReport {
        status: FeedbackBoundStatus::SolidInverseUnavailable,
        residual_infinity_upper: residual_inf, solid_inverse_infinity_upper: inverse,
        gain_infinity_upper: None, coupled_inverse_infinity_upper: None,
        state_error_infinity_upper: None, response_residual_infinity_upper: Vec::new(),
    };
    if let Some(inverse) = inverse {
        let mut gain = if gain_range { None } else { mul_up(inverse, perturbation_inf).ok() };
        if let Some(columns) = responses {
            let mut accumulated = scalar_vector(n, &mut work)?;
            let mut column_rhs = scalar_vector(n, &mut work)?;
            let mut response_range = gain_range;
            report.response_residual_infinity_upper.try_reserve_exact(p)
                .map_err(|_| GoalResidualError::Allocation)?;
            for (j, column) in columns.iter().enumerate() {
                for (i, out) in column_rhs.iter_mut().enumerate() {
                    work.tick()?;
                    // CSR lookup is logarithmic in this row's stored port count.
                    *out = injection.get(i, j);
                }
                let checked = enclose_goal_error(matrix, &column_rhs, column, &zeros, &zeros,
                    scaling, limits.solid, || (work.checkpoint)())?;
                let residual = checked.primal_residual_infinity_upper();
                report.response_residual_infinity_upper.push(residual);
                let error = mul_up(inverse, residual);
                match error {
                    Ok(error) => for (i, &value) in column.iter().enumerate() {
                        work.tick()?;
                        let term = add_up(value.abs(), error).and_then(|v| mul_up(v, c_norm[j]));
                        match term.and_then(|v| add_up(accumulated[i], v)) {
                            Ok(sum) => accumulated[i] = sum,
                            Err(_) => response_range = true,
                        }
                    },
                    Err(_) => response_range = true,
                }
            }
            if !response_range {
                let mut q = 0.0_f64;
                for value in accumulated { work.tick()?; q = q.max(value); }
                gain = Some(gain.map_or(q, |old| old.min(q)));
            }
        }
        report.gain_infinity_upper = gain;
        report.status = match gain {
            None => FeedbackBoundStatus::BoundNotRepresentable,
            Some(q) if q >= 1.0 => FeedbackBoundStatus::ContractionNotEstablished,
            Some(q) => {
                let bound = (|| {
                    let denominator = if q == 0.0 { 1.0 } else { down(1.0 - q)? };
                    if denominator <= 0.0 { return Err(GoalResidualError::ArithmeticRange); }
                    let coupled = if q == 0.0 { inverse } else { up(inverse / denominator)? };
                    Ok::<_, GoalResidualError>((coupled, mul_up(coupled, residual_inf)?))
                })();
                match bound {
                    Ok((coupled, error)) => {
                        report.coupled_inverse_infinity_upper = Some(coupled);
                        report.state_error_infinity_upper = Some(error);
                        FeedbackBoundStatus::Enclosed
                    }
                    Err(_) => FeedbackBoundStatus::BoundNotRepresentable,
                }
            }
        };
    }
    work.poll()?;
    Ok(report)
}

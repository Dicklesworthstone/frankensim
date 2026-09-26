//! Verify the small port system instead of requiring whole-state contraction.
//! H = A^-1 B, S = I-C H. For invertible S,
//! (A-B C)^-1 = (I+H S^-1 C) A^-1. Each interval entry of S includes the
//! residual-checked response-column error. Strict interval row dominance
//! then proves nonsingularity; no approximate Schur inverse is trusted.

use fs_sparse::Csr;

use super::{
    FeedbackBoundStatus, FeedbackInverseMethod, FeedbackResidualLimits,
    FeedbackResidualReport, GoalResidualError, ScalarEnclosure, Work,
    enclose_affine_feedback_error, scalar_vector,
};
use super::super::arithmetic::{add_up, down, mul_up, up};

/// Assess the same coupled system as [`enclose_affine_feedback_error`], with
/// an optional port-Schur fallback when its whole-state contraction check
/// cannot establish a bound. Already enclosed results are returned unchanged.
///
/// Supplied columns approximate H=A^-1 B. The original evaluator checks their
/// residuals, giving |H_ij-V_ij| <= alpha*||B_j-A V_j||_inf with a proved
/// alpha >= ||A^-1||_inf. These errors widen every entry of S=I-C H before
/// checking strict diagonal dominance. If beta bounds ||S^-1||_inf, then
/// alpha*(1+||H||_inf*beta*||C||_inf) bounds the coupled inverse. This proves
/// algebraic nonsingularity only, NOT stability of fixed-point iteration or
/// time evolution. In particular, a gain >= 1 may coexist with an enclosure.
///
/// The extra row/entry visits spend the REMAINING `max_verification_entries`
/// allowance. Insufficient remaining work or an inconclusive finite-range /
/// dominance check leaves the original no-bound report unchanged. No new
/// solve is performed, no dense state or port matrix is allocated, and scratch
/// storage is O(ports). All traversal, including long rows, polls cancellation.
///
/// # Errors
/// The original evaluator's input/work refusals, allocation failure, and
/// cancellation. No partially verified inverse is published.
#[allow(clippy::too_many_arguments)]
pub fn enclose_affine_feedback_error_with_schur(
    matrix: &Csr,
    rhs: &[f64],
    primal: &[f64],
    injection: &Csr,
    feedback: &Csr,
    offset: &[f64],
    responses: Option<&[Vec<f64>]>,
    scaling: Option<&[f64]>,
    limits: FeedbackResidualLimits,
    mut checkpoint: impl FnMut() -> bool,
) -> Result<FeedbackResidualReport, GoalResidualError> {
    let mut report = enclose_affine_feedback_error(
        matrix, rhs, primal, injection, feedback, offset, responses, scaling,
        limits, &mut checkpoint,
    )?;
    if report.status() == FeedbackBoundStatus::Enclosed {
        return Ok(report);
    }
    let mut work = Work { checkpoint, left: 512 };
    work.poll()?;
    let Some(columns) = responses else { return Ok(report) };
    let Some(alpha) = report.solid_inverse_infinity_upper() else { return Ok(report) };
    let n = matrix.nrows();
    let p = columns.len();
    // Counts are deterministic structural proxies, consistent with the base
    // evaluator. Include recomputed C norms, each interval Schur entry, the
    // response norm and per-port error/margin scratch work before allocation.
    let total = (|| {
        let base = n.checked_add(matrix.nnz())?.checked_mul(p.checked_add(1)?)?;
        let extra = p.checked_add(1)?.checked_mul(feedback.nnz())?
            .checked_add(n.checked_mul(p)?)?
            .checked_add(p.checked_mul(p)?)?
            .checked_add(p.checked_mul(3)?)?;
        base.checked_add(extra)
    })();
    if total.is_none_or(|needed| needed > limits.max_verification_entries) {
        return Ok(report);
    }
    let checked = schur_bound(alpha, feedback, columns,
        report.response_residual_infinity_upper(), n, &mut work);
    match checked {
        Ok(Some((schur, inverse))) => {
            if let Ok(error) = mul_up(inverse, report.residual_infinity_upper()) {
                report.status = FeedbackBoundStatus::Enclosed;
                report.coupled_inverse_infinity_upper = Some(inverse);
                report.state_error_infinity_upper = Some(error);
                report.inverse_method = Some(FeedbackInverseMethod::PortSchurDominance);
                report.schur_inverse_infinity_upper = Some(schur);
            }
        }
        Ok(None) | Err(GoalResidualError::ArithmeticRange) => {},
        Err(error) => return Err(error),
    }
    work.poll()?;
    Ok(report)
}

fn schur_bound<F: FnMut() -> bool>(
    alpha: f64,
    feedback: &Csr,
    columns: &[Vec<f64>],
    residuals: &[f64],
    n: usize,
    work: &mut Work<F>,
) -> Result<Option<(f64, f64)>, GoalResidualError> {
    let p = columns.len();
    // The base evaluator always checks all columns when the solid inverse is
    // present. Never treat an absent residual as an exact response solve.
    if residuals.len() != p { return Ok(None); }
    let mut errors = scalar_vector(p, work)?;
    let mut c_norms = scalar_vector(p, work)?;
    let mut c_norm = 0.0_f64;
    for (j, error) in errors.iter_mut().enumerate() {
        work.tick()?;
        *error = mul_up(alpha, residuals[j])?;
        for &value in feedback.row(j).1 {
            work.tick()?;
            c_norms[j] = add_up(c_norms[j], value.abs())?;
        }
        c_norm = c_norm.max(c_norms[j]);
    }
    let mut h_norm = 0.0_f64;
    for i in 0..n {
        let mut row = 0.0;
        for (column, &error) in columns.iter().zip(&errors) {
            work.tick()?;
            row = add_up(row, add_up(column[i].abs(), error)?)?;
        }
        h_norm = h_norm.max(row);
    }
    let mut minimum_margin = f64::INFINITY;
    for (i, &row_norm) in c_norms.iter().enumerate() {
        let mut diagonal_lower = 0.0;
        let mut off_diagonal_upper = 0.0;
        let (indices, values) = feedback.row(i);
        for (j, column) in columns.iter().enumerate() {
            work.tick()?;
            let mut entry = ScalarEnclosure::point(if i == j { 1.0 } else { 0.0 });
            for (&k, &value) in indices.iter().zip(values) {
                work.tick()?;
                entry = entry.add_scaled(ScalarEnclosure::point(column[k]), -value)?;
            }
            entry = entry.widen(mul_up(row_norm, errors[j])?)?;
            if i == j {
                diagonal_lower = if entry.lower() > 0.0 { entry.lower() }
                    else if entry.upper() < 0.0 { -entry.upper() }
                    else { 0.0 };
            } else {
                off_diagonal_upper = add_up(off_diagonal_upper, entry.magnitude_upper())?;
            }
        }
        let margin = down(diagonal_lower - off_diagonal_upper)?;
        if margin <= 0.0 { return Ok(None); }
        minimum_margin = minimum_margin.min(margin);
    }
    let beta = up(1.0 / minimum_margin)?;
    let amplification = add_up(1.0, mul_up(mul_up(h_norm, beta)?, c_norm)?)?;
    Ok(Some((beta, mul_up(alpha, amplification)?)))
}

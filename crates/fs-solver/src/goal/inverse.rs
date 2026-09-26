//! A bounded inverse check for matrices whose comparison matrix is not
//! strictly dominant. The supplied columns are numerical proposals, never
//! inverse authority: every entry of `I - A R` is enclosed independently.

use super::arithmetic::{add_up, down, mul_up, up};
use super::{
    GoalBoundStatus, GoalResidualError, GoalResidualLimits, GoalResidualReport, ScalarEnclosure,
    Work, enclose_goal_error,
};
use fs_sparse::Csr;

/// Enclose a discrete goal error, optionally proving an inverse bound from
/// column-major `columns[j][i] = R_ij`. The ordinary residual checker runs
/// first. Its established bound is preserved after validating the proposal.
/// Otherwise, outward evaluation proves `delta >= ||I - A R||_infinity`;
/// if `delta < 1`, then `||A^-1||_infinity <= ||R||_infinity / (1 - delta)`.
/// No symmetry, definiteness, exact inverse solve, or caller-supplied inverse
/// constant is assumed. Dual error is retained in the final goal enclosure.
///
/// Dimension is capped at 256. Both the dense entry count `n*n` and the
/// verification visits `n*(n + nnz)` must fit `limits.max_nonzeros`, in
/// addition to the ordinary sparse-system limits. No dense scratch is
/// allocated. Checkpoints cover input validation and every verification
/// entry, with at most 512 visits between polls and a final publication poll.
///
/// # Errors
/// Returns the ordinary checker refusals, plus typed proposal length,
/// finiteness, and structural-budget refusals. An inaccurate or singular
/// inverse proposal, or an unrepresentable inverse/remainder calculation,
/// leaves the full goal bound absent while retaining the residual enclosure.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn enclose_goal_error_with_inverse(
    matrix: &Csr,
    rhs: &[f64],
    primal: &[f64],
    goal: &[f64],
    dual: &[f64],
    scaling: Option<&[f64]>,
    columns: &[Vec<f64>],
    limits: GoalResidualLimits,
    mut checkpoint: impl FnMut() -> bool,
) -> Result<GoalResidualReport, GoalResidualError> {
    let mut report = enclose_goal_error(
        matrix,
        rhs,
        primal,
        goal,
        dual,
        scaling,
        limits,
        &mut checkpoint,
    )?;
    let mut work = Work {
        checkpoint,
        left: 512,
    };
    work.poll()?;
    let n = matrix.nrows();
    let visits = n
        .checked_add(matrix.nnz())
        .and_then(|row| n.checked_mul(row));
    for (field, required, allowed) in [
        ("inverse rows", n, 256),
        (
            "inverse entries",
            n.checked_mul(n).unwrap_or(usize::MAX),
            limits.max_nonzeros,
        ),
        (
            "inverse verification entries",
            visits.unwrap_or(usize::MAX),
            limits.max_nonzeros,
        ),
    ] {
        if required > allowed {
            return Err(GoalResidualError::Limit {
                field,
                required,
                allowed,
            });
        }
    }
    if columns.len() != n {
        return Err(GoalResidualError::Length {
            field: "inverse columns",
            expected: n,
            found: columns.len(),
        });
    }
    for column in columns {
        work.tick()?;
        if column.len() != n {
            return Err(GoalResidualError::Length {
                field: "inverse column",
                expected: n,
                found: column.len(),
            });
        }
        for (index, value) in column.iter().enumerate() {
            work.tick()?;
            if !value.is_finite() {
                return Err(GoalResidualError::NonFinite {
                    field: "inverse column",
                    index,
                });
            }
        }
    }
    if report.inverse_infinity_upper.is_some() {
        work.poll()?;
        return Ok(report);
    }
    let checked = (|| {
        let mut defect = 0.0_f64;
        let mut inverse_norm = 0.0_f64;
        for row in 0..n {
            let (indices, coefficients) = matrix.row(row);
            let mut defect_row = 0.0;
            let mut inverse_row = 0.0;
            for (j, column) in columns.iter().enumerate() {
                work.tick()?;
                inverse_row = add_up(inverse_row, column[row].abs())?;
                let mut entry = ScalarEnclosure::point(if row == j { 1.0 } else { 0.0 });
                for (&k, &a) in indices.iter().zip(coefficients) {
                    work.tick()?;
                    entry = entry.add_scaled(ScalarEnclosure::point(column[k]), -a)?;
                }
                defect_row = add_up(defect_row, entry.magnitude_upper())?;
            }
            defect = defect.max(defect_row);
            inverse_norm = inverse_norm.max(inverse_row);
        }
        if defect >= 1.0 {
            return Ok(None);
        }
        let norm = up(inverse_norm / down(1.0 - defect)?)?;
        let remainder = mul_up(
            mul_up(
                report.dual_residual_one_upper,
                report.primal_residual_infinity_upper,
            )?,
            norm,
        )?;
        Ok(Some((
            norm,
            remainder,
            report.weighted_residual.widen(remainder)?,
        )))
    })();
    match checked {
        Ok(Some((norm, remainder, error))) => {
            report.inverse_infinity_upper = Some(norm);
            report.dual_error_upper = Some(remainder);
            report.goal_error = Some(error);
            report.status = GoalBoundStatus::Enclosed;
        }
        Ok(None) => report.status = GoalBoundStatus::InverseBoundUnavailable,
        Err(GoalResidualError::ArithmeticRange) => {
            report.status = GoalBoundStatus::BoundNotRepresentable
        }
        Err(error) => return Err(error),
    }
    work.poll()?;
    Ok(report)
}

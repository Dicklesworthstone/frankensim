//! Storage increments through the original sparse and full element operators.

use super::*;
use fs_math::quadratic::{self, QuadraticError};

pub(super) struct StepEnergy {
    pub kinetic: f64,
    pub strain: f64,
    pub kinetic_change: Dd,
    pub strain_change: Dd,
    pub mechanical: Dd,
    pub support: Dd,
    pub thermal: Dd,
    pub damping: Dd,
}

pub(super) fn retained(value: Dd) -> Result<Dd, TransientError3> {
    number(value.hi)?;
    number(value.lo)?;
    number(value.to_f64())?;
    Ok(value)
}

fn failure(error: QuadraticError<TransientError3>) -> TransientError3 {
    match error {
        QuadraticError::Cancelled => TransientError3::Cancelled,
        QuadraticError::Operator(error) => error,
        _ => TransientError3::Invalid("invalid original storage increment"),
    }
}

pub(super) fn sparse_change(
    matrix: &Csr, before: &[f64], after: &[f64], cx: &Cx<'_>,
) -> Result<Dd, TransientError3> {
    quadratic::increment(before, after, |x, y, _| {
        // Exactly the original sparse row action; no reconstructed or lumped mass.
        Operator { m: matrix, k: matrix, cm: 1.0, ck: 0.0, cx }.apply(x, y);
        check(cx)
    }, &mut || cx.checkpoint().is_err())
        .map(|value| value.change).map_err(failure)
}

pub(super) fn full_change(
    problem: &TetLinearElasticProblem<'_>, mass: bool,
    before: &[f64], after: &[f64], nodal_reference_density: Option<&[f64]>, cx: &Cx<'_>,
) -> Result<Dd, TransientError3> {
    quadratic::increment(before, after, |x, y, _| {
        problem.apply_full_storage(x, y, mass, nodal_reference_density, cx)
            .map_err(TransientError3::from)
    }, &mut || cx.checkpoint().is_err())
        .map(|value| value.change).map_err(failure)
}

pub(super) fn source_work(
    before: &[f64], after: &[f64], force_before: &[f64], force_after: &[f64], cx: &Cx<'_>,
) -> Result<Dd, TransientError3> {
    quadratic::endpoint_work(before, after, force_before, force_after,
        &mut || cx.checkpoint().is_err()).map_err(|error| match error {
            QuadraticError::Cancelled => TransientError3::Cancelled,
            _ => TransientError3::Invalid("invalid endpoint force work"),
        })
}

/// Independent changes from the original symmetric metric coefficients:
/// [1/2 x1' delta(A) x1, 1/2 x0' delta(A) x1]. Mass uses velocities;
/// stiffness uses displacements. Differences precede the coordinate products.
/// The first changes endpoint storage; the second is parameter work. They are
/// different terms and neither is manufactured from an energy residual.
pub(super) fn coefficient_change(
    before: &Csr, after: &Csr, q0: &[f64], q1: &[f64], cx: &Cx<'_>,
) -> Result<[Dd; 2], TransientError3> {
    let mut result = [Dd::default(); 2];
    for row in 0..q0.len() {
        check(cx)?;
        let (old_cols, old_values) = before.row(row);
        let (new_cols, new_values) = after.row(row);
        let (mut old, mut new) = (0, 0);
        let mut action = Dd::default();
        while old < old_cols.len() || new < new_cols.len() {
            if (old + new) % 256 == 0 { check(cx)?; }
            let a = old_cols.get(old).copied().unwrap_or(usize::MAX);
            let b = new_cols.get(new).copied().unwrap_or(usize::MAX);
            let column = a.min(b);
            let mut coefficient = Dd::default();
            if b == column { coefficient = Dd::from_f64(new_values[new]); new += 1; }
            if a == column { coefficient = coefficient - Dd::from_f64(old_values[old]); old += 1; }
            action = action + coefficient * Dd::from_f64(q1[column]);
        }
        result[0] = result[0] + Dd::from_f64(0.5 * q1[row]) * action;
        result[1] = result[1] + Dd::from_f64(0.5 * q0[row]) * action;
    }
    Ok([retained(result[0])?, retained(result[1])?])
}

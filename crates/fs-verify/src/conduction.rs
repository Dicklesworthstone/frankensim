//! Solver-independent retention of verified whole-volume mean-temperature fields.
//!
//! The native FEM consumers are `fs_conduction::verification::solve_with_mean_bound`
//! and `bound_temperature_mean`, behind that crate's `thermal-verification`
//! feature. They admit the actual material and boundary model, solve the
//! same-operator unit-source homogeneous dual and optionally the primal, then
//! retain the original solution reports through this seam. The solver dependency
//! stays in that direction to avoid a cycle through fs-adjoint.
//!
//! The enclosure has exactly the scope of [`tet::tensor_mean_bound`]: the nominal
//! admitted linear tensor PDE on the declared conforming polyhedral domain.
//! It does not certify a point maximum, CAD fidelity, or parameter uncertainty.

use crate::tet::{self, FluxBudget, MeanBound, TensorTetProblem, TetError};

/// The caller's original primal and dual solutions, retained with their mean
/// enclosure. Solver reports and provenance are not replaced by this verifier.
#[derive(Debug)]
pub struct MeanTemperatureSolution<S> {
    pub primal: S,
    pub dual: S,
    pub bound: MeanBound,
}

/// A bound on a supplied primal field and the caller's actual dual solution.
/// No primal solver report is invented for a supplied field.
#[derive(Debug)]
pub struct MeanFieldBound<S> {
    pub dual: S,
    pub bound: MeanBound,
}

/// Enclose a supplied P1 temperature field using the actual dual field in the
/// retained solution. Fields must use the declared problem's vertex order;
/// `dual_temperature` extracts that field without copying or substituting a
/// solver convergence tolerance for the verifier's residual correction.
///
/// The dual has the same tensor operator, homogeneous boundary data and unit
/// source. Both discretization and algebraic error remain covered for inexact
/// fields. A refused model, exhausted flux budget or cancellation yields no
/// result. Exact tensor coefficients are passed unchanged to the verifier.
pub fn bound_mean_temperature<S>(
    problem: &TensorTetProblem<'_>,
    temperature: &[f64],
    dual: S,
    dual_temperature: impl for<'a> Fn(&'a S) -> &'a [f64],
    budget: FluxBudget,
    keep_going: impl FnMut() -> bool,
) -> Result<MeanFieldBound<S>, TetError> {
    let bound = tet::tensor_mean_bound(
        problem,
        temperature,
        dual_temperature(&dual),
        budget,
        keep_going,
    )?;
    Ok(MeanFieldBound { dual, bound })
}

#[cfg(test)]
#[path = "conduction/tensor_tests.rs"]
mod tensor_tests;
#[cfg(test)]
mod tests;

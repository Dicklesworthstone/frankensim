//! Continuous original-storage tangents for the shared harmonic-balance solver.

use super::{
    HbBudget, HbRefinementBudget, HbScales, PeriodicPortInput, PhsOrbitError, PhsOrbitProblem, bad,
    forced_problem,
};
use core::cell::RefCell;
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, ScaledHbOrbit, solve_hb_adaptive_scaled_seeded,
    solve_hb_scaled_seeded,
};
use fs_phs::PortHamiltonian;
use std::borrow::Cow;

/// Complete continuous Hessian action of the original Hamiltonian.
///
/// Given physical state `x` and physical coordinate direction `v`, write
/// `D(gradient H)(x) v` into every output component and return true. This is
/// the derivative of the storage's actual effort law, including every explicit
/// material-memory coordinate, rather than a discrete-step Jacobian or a
/// quadratic approximation to a nonlinear Hamiltonian. The callback must be
/// deterministic and read-only. False or unwritten/nonfinite outputs refuse
/// the solve; no finite-difference replacement is used for a declared action.
pub type PhsStorageHessianAction<'a> = dyn Fn(&[f64], &[f64], &mut [f64]) -> bool + 'a;

pub(super) enum StorageTangent<'a> {
    Unavailable,
    Quadratic {
        matrix: Cow<'a, [f64]>,
        prepared: RefCell<Option<Vec<f64>>>,
    },
    Hessian {
        action: &'a PhsStorageHessianAction<'a>,
        scratch: RefCell<(Vec<f64>, Vec<f64>)>,
    },
}

impl<'a> StorageTangent<'a> {
    pub(super) fn from_storage(system: &'a PortHamiltonian) -> Self {
        let matrix = if let Some(matrix) = system.storage().exact_quadratic_matrix() {
            Cow::Borrowed(matrix)
        } else if let Some(part) = system.storage().quadratic_part().filter(|part| part.exact) {
            // Composed and relaxation storage expose their complete original
            // quadratic form through owned metadata, including memory rows.
            // A merely quadratic PART of a nonlinear law is never sufficient.
            Cow::Owned(part.matrix)
        } else {
            return Self::Unavailable;
        };
        Self::Quadratic {
            matrix,
            prepared: RefCell::new(None),
        }
    }

    pub(super) fn hessian(n: usize, action: &'a PhsStorageHessianAction<'a>) -> Self {
        Self::Hessian {
            action,
            scratch: RefCell::new((vec![0.0; n], vec![0.0; n])),
        }
    }

    pub(super) fn available(&self) -> bool {
        !matches!(self, Self::Unavailable)
    }
}

impl PhsOrbitProblem<'_> {
    pub(super) fn state_jacobian(&self, _time: f64, state: &[f64], out: &mut [f64]) -> bool {
        let n = self.system.state_dim();
        let Some(square) = n.checked_mul(n) else {
            return false;
        };
        if state.len() != n || out.len() != square {
            return false;
        }
        let (j, r, _) = self.system.structure();
        match &self.tangent {
            StorageTangent::Unavailable => false,
            StorageTangent::Quadratic { matrix, prepared } => {
                if matrix.len() != square || matrix.iter().any(|v| !v.is_finite()) {
                    return false;
                }
                let mut prepared = prepared.borrow_mut();
                let tangent = prepared.get_or_insert_with(|| {
                    let mut result = vec![0.0; square];
                    // Constant original vector-field tangent: prepare only once,
                    // and skip exactly absent J-R entries without thresholding.
                    for row in 0..n {
                        for column in 0..n {
                            let mut value = 0.0;
                            let mut correction = 0.0;
                            for k in 0..n {
                                let flow = j[row * n + k] - r[row * n + k];
                                if flow == 0.0 {
                                    continue;
                                }
                                let term = flow * matrix[k * n + column];
                                let adjusted = term - correction;
                                let next = value + adjusted;
                                correction = (next - value) - adjusted;
                                value = next;
                            }
                            result[row * n + column] = value;
                        }
                    }
                    result
                });
                out.copy_from_slice(tangent);
                out.iter().all(|v| v.is_finite())
            }
            StorageTangent::Hessian { action, scratch } => {
                let mut scratch = scratch.borrow_mut();
                let (direction, effort) = &mut *scratch;
                for column in 0..n {
                    direction.fill(0.0);
                    direction[column] = 1.0;
                    effort.fill(f64::NAN);
                    if !action(state, direction, effort) || effort.iter().any(|v| !v.is_finite()) {
                        return false;
                    }
                    for row in 0..n {
                        let mut value = 0.0;
                        let mut correction = 0.0;
                        for k in 0..n {
                            let flow = j[row * n + k] - r[row * n + k];
                            if flow == 0.0 {
                                continue;
                            }
                            let adjusted = flow * effort[k] - correction;
                            let next = value + adjusted;
                            correction = (next - value) - adjusted;
                            value = next;
                        }
                        out[row * n + column] = value;
                    }
                }
                out.iter().all(|v| v.is_finite())
            }
        }
    }
}

/// Solve a forced periodic response using the original nonlinear storage's
/// explicitly supplied continuous Hessian action and physical Newton scales.
///
/// The original vector field `(J-R) gradient H + G u(t)` supplies every residual;
/// its tangent `(J-R) D(gradient H)` is assembled once per AFT sample by fs-orbit.
/// All Fourier coefficients, phases, memory coordinates, sampling checks and
/// source refusals keep the semantics of the ordinary scaled seeded wrapper.
/// Exact quadratic storage is already accelerated automatically by that wrapper;
/// this entry admits an explicit action for nonlinear storage.
///
/// # Errors
/// Invalid forcing/scales/seeds, source evaluation failure, declared tangent
/// refusal or fs-orbit's ordinary numerical refusals. A tangent is required only
/// when Newton needs a step; an already converged independently checked seed
/// need not evaluate it.
pub fn solve_forced_harmonic_balance_analytic_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: Vec<f64>,
    budget: &HbBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<ScaledHbOrbit, PhsOrbitError> {
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    let mut problem = forced_problem(system, input)?;
    problem.tangent = StorageTangent::hessian(system.state_dim(), hessian);
    let result = solve_hb_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        scales,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

/// The same original-storage analytic solve with adaptive harmonic refinement.
///
/// The supplied Hessian action is used at every admitted order. All independent
/// finite-band residual and waveform checks remain owned by fs-orbit; an exact
/// tangent does not certify an omitted infinite harmonic tail or orbit stability.
///
/// # Errors
/// The fixed-order analytic wrapper's refusals and unresolved adaptive budgets.
pub fn solve_forced_harmonic_balance_analytic_adaptive_scaled_seeded(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    seed: Vec<f64>,
    budget: &HbBudget,
    refinement: &HbRefinementBudget,
    scales: &HbScales,
    hessian: &PhsStorageHessianAction<'_>,
) -> Result<AdaptiveHbOrbit, PhsOrbitError> {
    if input.highest_nonzero_harmonic() > budget.harmonics {
        return Err(bad(
            "nonzero forcing harmonic exceeds the retained HB truncation",
        ));
    }
    let mut problem = forced_problem(system, input)?;
    problem.tangent = StorageTangent::hessian(system.state_dim(), hessian);
    let result = solve_hb_adaptive_scaled_seeded(
        &problem,
        HbAnchor::Forced {
            omega: input.omega_rad_s(),
        },
        seed,
        budget,
        refinement,
        scales,
    );
    if problem.source_refused() {
        return Err(PhsOrbitError::SourceRefused);
    }
    result.map_err(PhsOrbitError::Orbit)
}

#[cfg(test)]
mod tests;

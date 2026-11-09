//! Bounded Newton solve with the original aqueous tridiagonal operator.
//!
//! `exchange(c)` is the actual represented immobile amount change and its
//! nonnegative concentration tangent. Native adjacent-pivot LU is reused for
//! every Jacobian. Box-preserving backtracking decreases the original amount
//! residual; it supplies neither a global convergence nor a time-error claim.

use super::{DiluteSoluteError, Result, finite, nonnegative, poll, positive};
use fs_la::tridiagonal::TridiagonalWorkspace;

#[derive(Clone, Copy)]
pub(in super::super) struct Controls {
    pub tolerance_mol_m2: f64,
    pub maximum_iterations: usize,
    pub maximum_line_searches: usize,
    pub remaining_linear_solves: usize,
}

#[allow(clippy::too_many_arguments)]
fn evaluate(
    concentration: &[f64],
    lower: &[f64],
    diagonal: &[f64],
    upper: &[f64],
    rhs: &[f64],
    residual: &mut [f64],
    tangent: &mut [f64],
    response: &mut impl FnMut(usize, f64) -> Result<(f64, f64)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut maximum = 0.0_f64;
    for i in 0..concentration.len() {
        if i % 128 == 0 {
            poll(cancelled)?;
        }
        let (exchange, derivative) = response(i, concentration[i])?;
        if !exchange.is_finite() || !nonnegative(derivative) {
            return Err(DiluteSoluteError::Invalid(
                "invalid nonlinear storage response",
            ));
        }
        // Neumaier retains low terms when transport and inventory magnitudes
        // differ. The native owner independently reconstructs physical faces
        // and cell/whole-column changes after this solve.
        let terms = [
            finite(diagonal[i] * concentration[i])?,
            if i == 0 {
                0.0
            } else {
                finite(lower[i - 1] * concentration[i - 1])?
            },
            if i + 1 == concentration.len() {
                0.0
            } else {
                finite(upper[i] * concentration[i + 1])?
            },
            exchange,
            -rhs[i],
        ];
        let (mut total, mut correction) = (0.0_f64, 0.0_f64);
        for x in terms {
            let next = finite(total + x)?;
            correction += if total.abs() >= x.abs() {
                (total - next) + x
            } else {
                (x - next) + total
            };
            total = next;
        }
        residual[i] = finite(total + correction)?;
        tangent[i] = derivative;
        maximum = maximum.max(residual[i].abs());
    }
    Ok(maximum)
}

/// Eight cell arrays are allocated once per constituent; all trial evaluations
/// are charged by the source response closure, including rejected proposals.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(in super::super) fn solve(
    species: usize,
    lower: &[f64],
    diagonal: &[f64],
    upper: &[f64],
    rhs: &mut [f64],
    initial: &[f64],
    limits: impl Fn(usize) -> [f64; 2],
    solver: &mut TridiagonalWorkspace,
    controls: Controls,
    mut response: impl FnMut(usize, f64) -> Result<(f64, f64)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<usize> {
    let n = rhs.len();
    if initial.len() != n || !positive(controls.tolerance_mol_m2) {
        return Err(DiluteSoluteError::Invalid(
            "invalid nonlinear storage solve extent or tolerance",
        ));
    }
    let mut current = initial.to_vec();
    let mut residual = vec![0.0; n];
    let mut tangent = vec![0.0; n];
    let mut increment = vec![0.0; n];
    let mut jacobian = vec![0.0; n];
    let mut trial = vec![0.0; n];
    let mut trial_residual = vec![0.0; n];
    let mut trial_tangent = vec![0.0; n];
    for (i, &c) in current.iter().enumerate() {
        let [low, high] = limits(i);
        if !nonnegative(low)
            || !nonnegative(high)
            || low > high
            || !c.is_finite()
            || c < low
            || c > high
        {
            return Err(DiluteSoluteError::Invalid(
                "nonlinear initial concentration outside original source support",
            ));
        }
    }
    let mut maximum = evaluate(
        &current,
        lower,
        diagonal,
        upper,
        rhs,
        &mut residual,
        &mut tangent,
        &mut response,
        cancelled,
    )?;
    for solves in 0..=controls.maximum_iterations {
        poll(cancelled)?;
        if maximum <= controls.tolerance_mol_m2 {
            rhs.copy_from_slice(&current);
            return Ok(solves);
        }
        if solves == controls.maximum_iterations {
            break;
        }
        if solves >= controls.remaining_linear_solves {
            return Err(DiluteSoluteError::NonlinearWork {
                required_linear_solves: solves + 1,
                maximum_linear_solves: controls.remaining_linear_solves,
            });
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            jacobian[i] = finite(diagonal[i] + tangent[i])?;
            increment[i] = -residual[i];
        }
        solver.factor(lower, &jacobian, upper, &mut *cancelled)?;
        solver.solve_factored_in_place(&mut increment, &mut *cancelled)?;
        let mut fraction = 1.0_f64;
        for i in 0..n {
            let [low, high] = limits(i);
            if increment[i] > 0.0 {
                fraction = fraction.min((high - current[i]) / increment[i]);
            } else if increment[i] < 0.0 {
                fraction = fraction.min((low - current[i]) / increment[i]);
            }
        }
        if fraction < 1.0 {
            fraction *= 0.5;
        }
        let mut accepted = false;
        for _ in 0..controls.maximum_line_searches {
            poll(cancelled)?;
            if !positive(fraction) {
                break;
            }
            let mut changed = false;
            let mut inside = true;
            for i in 0..n {
                if i % 128 == 0 {
                    poll(cancelled)?;
                }
                trial[i] = finite(fraction.mul_add(increment[i], current[i]))?;
                let [low, high] = limits(i);
                inside &= trial[i] >= low && trial[i] <= high;
                changed |= trial[i] != current[i];
            }
            if !changed {
                break;
            }
            if inside {
                let next = evaluate(
                    &trial,
                    lower,
                    diagonal,
                    upper,
                    rhs,
                    &mut trial_residual,
                    &mut trial_tangent,
                    &mut response,
                    cancelled,
                )?;
                if next <= controls.tolerance_mol_m2 || next < maximum * (1.0 - 1e-4 * fraction) {
                    core::mem::swap(&mut current, &mut trial);
                    core::mem::swap(&mut residual, &mut trial_residual);
                    core::mem::swap(&mut tangent, &mut trial_tangent);
                    maximum = next;
                    accepted = true;
                    break;
                }
            }
            fraction *= 0.5;
        }
        if !accepted {
            return Err(DiluteSoluteError::Nonlinear {
                species,
                linear_solves: solves + 1,
                maximum_residual_mol_m2: maximum,
            });
        }
    }
    Err(DiluteSoluteError::Nonlinear {
        species,
        linear_solves: controls.maximum_iterations,
        maximum_residual_mol_m2: maximum,
    })
}

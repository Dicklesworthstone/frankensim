//! Retained pressure constraints for a rational terminal with exactly zero D.
//!
//! Original states obey x' = A x + F q + b u, H x = 0, where q is outward
//! terminal flow. Differentiating pressure continuity once gives
//! (-H F) q = H A x + H b u. Only these zero-storage flow multipliers are
//! eliminated; no DG or rational state is dropped or assigned invented storage.
//! A singular differentiated constraint needs a higher-index treatment and
//! refuses. This is the index-two Hessenberg case, not an epsilon resistance.

use super::super::terminal::TerminalSolve;
use super::{AcousticRealizeError, invalid};
use fs_exec::CancelGate;

pub(super) struct ConstrainedTerminal {
    pub generator: Vec<f64>,
    pub inlet_column: Vec<f64>,
    pub flow_columns: Vec<Vec<f64>>,
    pub inlet_flow_factors: Vec<f64>,
    /// Row-equilibrated original pressure constraints, used independently at
    /// initial and accepted candidate states. Scaling does not change H x = 0.
    pub constraints: Vec<Vec<f64>>,
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

/// Project the original vector field onto the original pressure manifold.
/// The admitted initial state must already satisfy H x = 0; it is never reset.
pub(super) fn eliminate(
    source: &[f64],
    inputs: &[Vec<f64>],
    constraints: &[Vec<f64>],
    inlet: &[f64],
    gate: &CancelGate,
) -> Result<ConstrainedTerminal, AcousticRealizeError> {
    poll(gate)?;
    let n = inlet.len();
    let p = constraints.len();
    if n == 0
        || n > 256
        || p == 0
        || p > n
        || source.len() != n * n
        || inputs.len() != p
        || inputs.iter().chain(constraints).any(|row| row.len() != n)
        || source
            .iter()
            .chain(inlet)
            .chain(inputs.iter().flatten())
            .chain(constraints.iter().flatten())
            .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "zero-direct terminal needs complete bounded finite original matrices",
        ));
    }
    let mut h = constraints.to_vec();
    for row in &mut h {
        let scale = row.iter().map(|v| v.abs()).fold(0.0, f64::max);
        if scale == 0.0 {
            return Err(invalid(
                "zero-direct terminal has an unresolved zero pressure constraint",
            ));
        }
        for v in row {
            *v /= scale;
        }
    }
    let mut constraint_mass = vec![0.0; p * p];
    for i in 0..p {
        poll(gate)?;
        for j in 0..p {
            constraint_mass[i * p + j] =
                -h[i].iter().zip(&inputs[j]).map(|(a, b)| a * b).sum::<f64>();
        }
    }
    let solve = TerminalSolve::new(&constraint_mass, p).map_err(|_| invalid(
        "zero-direct terminal differentiated pressure constraints are singular, nonpositive or ill-conditioned; higher-index elimination is unsupported",
    ))?;
    let mut flow_columns = vec![vec![0.0; n]; p];
    for col in 0..n {
        poll(gate)?;
        let rhs: Vec<_> = h
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(i, h)| h * source[i * n + col])
                    .sum()
            })
            .collect();
        for (row, value) in solve.apply(&rhs)?.into_iter().enumerate() {
            flow_columns[row][col] = value;
        }
    }
    let rhs: Vec<_> = h
        .iter()
        .map(|row| row.iter().zip(inlet).map(|(h, b)| h * b).sum())
        .collect();
    let inlet_flow_factors = solve.apply(&rhs)?;
    let mut generator = source.to_vec();
    let mut inlet_column = inlet.to_vec();
    for i in 0..n {
        poll(gate)?;
        for k in 0..p {
            inlet_column[i] += inputs[k][i] * inlet_flow_factors[k];
            for j in 0..n {
                generator[i * n + j] += inputs[k][i] * flow_columns[k][j];
            }
        }
    }
    if generator
        .iter()
        .chain(&inlet_column)
        .chain(flow_columns.iter().flatten())
        .chain(&inlet_flow_factors)
        .any(|v| !v.is_finite())
    {
        return Err(invalid(
            "zero-direct terminal constraint elimination is not representable",
        ));
    }
    // Independently check both differentiated identities, rather than infer
    // consistency from successful factorization of -H F.
    let tolerance = 256.0 * f64::EPSILON * n as f64;
    for row in &h {
        poll(gate)?;
        for j in 0..=n {
            let mut residual = 0.0;
            let mut scale = 0.0;
            for i in 0..n {
                residual += row[i]
                    * if j == n {
                        inlet_column[i]
                    } else {
                        generator[i * n + j]
                    };
                // Use the original equation's terms, before the cancellation
                // that enforces zero pressure. A scale formed only from the
                // projected row can collapse to its own roundoff residual.
                scale += (row[i] * if j == n { inlet[i] } else { source[i * n + j] }).abs();
                for k in 0..p {
                    scale += (row[i]
                        * inputs[k][i]
                        * if j == n {
                            inlet_flow_factors[k]
                        } else {
                            flow_columns[k][j]
                        })
                    .abs();
                }
            }
            if !residual.is_finite()
                || !scale.is_finite()
                || residual.abs() > tolerance * scale.max(f64::MIN_POSITIVE)
            {
                return Err(invalid(
                    "zero-direct terminal fails its differentiated pressure constraint",
                ));
            }
        }
    }
    Ok(ConstrainedTerminal {
        generator,
        inlet_column,
        flow_columns,
        inlet_flow_factors,
        constraints: h,
    })
}

#[cfg(test)]
mod tests;

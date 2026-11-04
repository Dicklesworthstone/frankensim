//! Finish only the part held back by cycles; the original DAG path stays exact.
//!
//! At cell i, `sum(inlet flux + q_ji Y_j) = q_out,i Y_i`. Row scaling by
//! positive incoming flow gives diagonal `q_out/q_in` and nonpositive
//! off-diagonals `-q_ji/q_in` (the balanced-flow `I - P` form).
//! Inlet reachability supplies the composition boundary condition. The same LU
//! solves every species and one positive comparison vector; no artificial leak,
//! concentration floor, least-squares regularization or new flow law is added.

use std::collections::{BTreeMap, BTreeSet};

use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};

use super::{
    CellId, DepositionError, MAX_RECIRCULATION_CELLS, RECIRCULATION_FRACTION_TOLERANCE,
    close_generated_simplex, finite, invalid,
};

fn numerical(error: LuWorkspaceError) -> DepositionError {
    match error {
        LuWorkspaceError::Capacity => DepositionError::Budget("stationary recirculation workspace"),
        LuWorkspaceError::Singular { .. } => {
            invalid("stationary recirculation has no unique composition")
        }
        _ => DepositionError::Arithmetic,
    }
}

/// A positive v with Av > 0 makes this Z-matrix a nonsingular M-matrix, so
/// `||A^-1||_inf <= max(v) / min(Av)`. This follows by multiplying Av >= m*1
/// by its nonnegative inverse. See Higham, "What Is an M-Matrix?" (2021),
/// https://nhigham.com/2021/03/16/what-is-an-m-matrix/ . The finite-f64 margins
/// below are conservative admission checks, not directed-rounding certificates.
fn inverse_norm_bound(
    lu: &mut LuWorkspace,
    matrix: &[f64],
    n: usize,
) -> Result<f64, DepositionError> {
    let mut positive = vec![0.0; n];
    lu.solve_factored_into(&vec![1.0; n], &mut positive)
        .map_err(numerical)?;
    if positive.iter().any(|v| *v <= 0.0) {
        return Err(invalid(
            "stationary recirculation is singular or ill-conditioned",
        ));
    }
    let mut lower = f64::INFINITY;
    for row in matrix.chunks_exact(n) {
        let mut product = 0.0;
        let mut magnitude = 0.0;
        for (&a, &v) in row.iter().zip(&positive) {
            product = finite(a.mul_add(v, product))?;
            magnitude = finite(magnitude + finite(a * v)?.abs())?;
        }
        let allowance = finite(8.0 * f64::EPSILON * (n + 2) as f64 * magnitude)?;
        lower = lower.min(product - allowance);
    }
    if lower <= 0.0 {
        return Err(invalid(
            "stationary recirculation is singular or ill-conditioned",
        ));
    }
    finite(positive.iter().copied().fold(0.0, f64::max) / lower)
}

pub(super) fn complete(
    rates: &BTreeMap<CellId, [f64; 2]>,
    edges: &BTreeMap<CellId, Vec<(CellId, f64)>>,
    incoming_species: &BTreeMap<CellId, Vec<f64>>,
    mut uniform: BTreeMap<CellId, Option<Vec<f64>>>,
    cells: &mut BTreeMap<CellId, Vec<f64>>,
    count: usize,
) -> Result<f64, DepositionError> {
    let ids: Vec<_> = rates
        .keys()
        .copied()
        .filter(|id| !cells.contains_key(id))
        .collect();
    let n = ids.len();
    if n > MAX_RECIRCULATION_CELLS {
        return Err(DepositionError::Budget("stationary recirculation cells"));
    }
    let indices: BTreeMap<_, _> = ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();

    // Propagate inlet identity separately from floating arithmetic. Each cell
    // changes at most twice: unreachable -> uniform -> mixed. Thus uniform
    // supplied fractions remain bit-identical even after recirculation.
    let mut ready: BTreeSet<_> = uniform.keys().copied().collect();
    while let Some(up) = ready.pop_first() {
        let source = uniform[&up].clone();
        for &(down, _) in edges.get(&up).map_or(&[][..], Vec::as_slice) {
            if !indices.contains_key(&down) {
                continue;
            }
            let changed = match uniform.entry(down) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(source.clone());
                    true
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    if entry.get().is_some() && entry.get() != &source {
                        entry.insert(None);
                        true
                    } else {
                        false
                    }
                }
            };
            if changed {
                ready.insert(down);
            }
        }
    }
    if ids.iter().any(|id| !uniform.contains_key(id)) {
        return Err(invalid("stationary species component has no inlet"));
    }

    let mut matrix = vec![0.0; n * n];
    let mut rhs = vec![vec![0.0; n]; count];
    for (i, id) in ids.iter().enumerate() {
        matrix[i * n + i] = finite(rates[id][1] / rates[id][0])?;
        if let Some(fluxes) = incoming_species.get(id) {
            for (s, &flux) in fluxes.iter().enumerate() {
                rhs[s][i] = finite(flux / rates[id][0])?;
            }
        }
    }
    for (&up, outgoing) in edges {
        let Some(&j) = indices.get(&up) else {
            continue;
        };
        for &(down, rate) in outgoing {
            let Some(&i) = indices.get(&down) else {
                continue;
            };
            let fraction = finite(rate / rates[&down][0])?;
            if fraction == 0.0 {
                return Err(invalid("stationary recirculation rate is unrepresentable"));
            }
            matrix[i * n + j] = finite(matrix[i * n + j] - fraction)?;
        }
    }
    let mut lu = LuWorkspace::new(n).map_err(numerical)?;
    lu.factor(&matrix).map_err(numerical)?;
    let inverse_bound = inverse_norm_bound(&mut lu, &matrix, n)?;
    let mut fractions = vec![vec![0.0; count]; n];
    let mut solution = vec![0.0; n];
    for (s, source) in rhs.iter().enumerate() {
        lu.solve_factored_into(source, &mut solution)
            .map_err(numerical)?;
        for (row, &value) in fractions.iter_mut().zip(&solution) {
            row[s] = value;
        }
    }
    for (i, id) in ids.iter().enumerate() {
        if let Some(exact) = &uniform[id] {
            fractions[i].clone_from(exact);
        } else {
            close_generated_simplex(&mut fractions[i])?;
        }
    }

    // Check the actual retained simplex, after roundoff-only closure, against
    // the original scaled equations. A tiny residual alone is insufficient
    // near a sealed loop: its amplification must also stay bounded.
    let mut residual_with_roundoff: f64 = 0.0;
    for i in 0..n {
        for (s, source) in rhs.iter().enumerate() {
            let mut residual = -source[i];
            let mut magnitude = source[i].abs();
            for j in 0..n {
                let term = finite(matrix[i * n + j] * fractions[j][s])?;
                residual = finite(residual + term)?;
                magnitude = finite(magnitude + term.abs())?;
            }
            let allowance = finite(16.0 * f64::EPSILON * (n + count + 2) as f64 * magnitude)?;
            residual_with_roundoff = residual_with_roundoff.max(residual.abs() + allowance);
        }
    }
    if finite(inverse_bound * residual_with_roundoff)? > RECIRCULATION_FRACTION_TOLERANCE {
        return Err(invalid(
            "stationary recirculation exceeds the mass-fraction error allowance",
        ));
    }

    // Independently reconstruct kg/s from the SAME unscaled connections used
    // by enthalpy. No algebraic solver residual substitutes for species balance.
    let mut fluxes: Vec<_> = ids
        .iter()
        .map(|id| {
            incoming_species
                .get(id)
                .cloned()
                .unwrap_or_else(|| vec![0.0; count])
        })
        .collect();
    for (&up, outgoing) in edges {
        let Some(&j) = indices.get(&up) else {
            continue;
        };
        for &(down, rate) in outgoing {
            let Some(&i) = indices.get(&down) else {
                continue;
            };
            for (flux, &fraction) in fluxes[i].iter_mut().zip(&fractions[j]) {
                *flux = finite(*flux + finite(rate * fraction)?)?;
            }
        }
    }
    let mut maximum_residual_kg_s: f64 = 0.0;
    for (i, id) in ids.iter().enumerate() {
        let [incoming, outgoing] = rates[id];
        for (&flux, &fraction) in fluxes[i].iter().zip(&fractions[i]) {
            let residual = finite(flux - finite(outgoing * fraction)?)?.abs();
            if residual > 128.0 * f64::EPSILON * count as f64 * incoming.max(outgoing) {
                return Err(invalid("stationary species flux does not close"));
            }
            maximum_residual_kg_s = maximum_residual_kg_s.max(residual);
        }
    }
    for (id, row) in ids.into_iter().zip(fractions) {
        cells.insert(id, row);
    }
    Ok(maximum_residual_kg_s)
}

//! Optional exact elimination of an acyclic thermal Jacobian's off-diagonal
//! support. Opposite entries need not be equal, nonzero or the same sign.
const NONE: usize = usize::MAX;

#[derive(Debug)]
pub(super) enum FactorError<E> {
    Capacity,
    Checkpoint(E),
}

struct Elimination {
    cell: usize,
    parent: usize,
    toward_parent: f64,
    multiplier: f64,
}

pub(super) struct Forest {
    pivots: Vec<f64>,
    order: Vec<Elimination>,
}

fn filled<T: Clone, E>(n: usize, value: T) -> Result<Vec<T>, FactorError<E>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| FactorError::Capacity)?;
    values.resize(n, value);
    Ok(values)
}

impl Forest {
    /// Build in O(cells + contributions), including parallel interfaces and
    /// directed-transfer entries. A cycle or unsafe pivot selects the original
    /// diagonal preconditioner; neither changes nor rejects the physical matrix.
    pub(super) fn factor<E>(
        diagonal: &[f64],
        entries: &[(usize, usize, f64)],
        checkpoint: &mut impl FnMut() -> Result<(), E>,
    ) -> Result<Option<Self>, FactorError<E>> {
        checkpoint().map_err(FactorError::Checkpoint)?;
        if entries.is_empty() {
            return Ok(None);
        }
        let n = diagonal.len();
        let mut head = filled(n, NONE)?;
        let mut next = filled(entries.len(), [NONE; 2])?;
        for (index, &(a, b, value)) in entries.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint().map_err(FactorError::Checkpoint)?;
            }
            if a >= n || b >= n || a == b || !value.is_finite() {
                return Ok(None);
            }
            next[index] = [head[a], head[b]];
            head[a] = index;
            head[b] = index;
        }

        // Both endpoint incidence lists encounter parallel contributions in
        // the same reverse insertion order. Their first contribution is thus
        // a shared edge identity. Marker arrays combine duplicates without a
        // map, sort, per-cell allocation or assumption about cell-ID ordering.
        let mut coefficients = filled(entries.len(), [0.0; 2])?;
        let mut seen = filled(n, NONE)?;
        let mut first = filled(n, NONE)?;
        let mut degree = filled(n, 0_usize)?;
        let mut incident = filled(n, 0_usize)?;
        let mut neighbours = Vec::new();
        neighbours
            .try_reserve_exact(n.min(entries.len()))
            .map_err(|_| FactorError::Capacity)?;
        let mut visited = 0_usize;
        for cell in 0..n {
            if cell % 256 == 0 {
                checkpoint().map_err(FactorError::Checkpoint)?;
            }
            neighbours.clear();
            let mut index = head[cell];
            while index != NONE {
                if visited % 256 == 0 {
                    checkpoint().map_err(FactorError::Checkpoint)?;
                }
                visited += 1;
                let (a, b, value) = entries[index];
                let side = usize::from(cell != a);
                let other = if side == 0 { b } else { a };
                if seen[other] != cell {
                    seen[other] = cell;
                    first[other] = index;
                    neighbours.push(index);
                }
                // The lower-index endpoint assembles both directional values
                // before the upper endpoint reads their shared representative.
                if cell < other {
                    let coefficient = &mut coefficients[first[other]][usize::from(a > b)];
                    *coefficient += value;
                    if !coefficient.is_finite() {
                        return Ok(None);
                    }
                }
                index = next[index][side];
            }
            for (i, &edge) in neighbours.iter().enumerate() {
                if i % 256 == 0 {
                    checkpoint().map_err(FactorError::Checkpoint)?;
                }
                if coefficients[edge] != [0.0; 2] {
                    degree[cell] += 1;
                    incident[cell] ^= edge;
                }
            }
        }

        let mut pivots = filled(n, 0.0)?;
        let mut leaves = Vec::new();
        let mut order = Vec::new();
        leaves
            .try_reserve_exact(n)
            .map_err(|_| FactorError::Capacity)?;
        order
            .try_reserve_exact(n)
            .map_err(|_| FactorError::Capacity)?;
        for (cell, &value) in diagonal.iter().enumerate() {
            if cell % 256 == 0 {
                checkpoint().map_err(FactorError::Checkpoint)?;
            }
            pivots[cell] = value;
            if degree[cell] <= 1 {
                leaves.push(cell);
            }
        }
        let mut cursor = 0;
        while cursor < leaves.len() {
            if cursor % 256 == 0 {
                checkpoint().map_err(FactorError::Checkpoint)?;
            }
            let cell = leaves[cursor];
            cursor += 1;
            let pivot = pivots[cell];
            if !pivot.is_finite() || pivot.abs() < f64::MIN_POSITIVE {
                return Ok(None);
            }
            if degree[cell] == 0 {
                order.push(Elimination {
                    cell,
                    parent: NONE,
                    toward_parent: 0.0,
                    multiplier: 0.0,
                });
                continue;
            }
            // XOR of the remaining edge identities identifies a leaf's sole
            // neighbour in constant time, even at a very high-degree hub.
            let edge = incident[cell];
            let (a, b, _) = entries[edge];
            let parent = if a == cell { b } else { a };
            let side = usize::from(cell > parent);
            let toward_parent = coefficients[edge][side];
            let multiplier = coefficients[edge][1 - side] / pivot;
            pivots[parent] = (-multiplier).mul_add(toward_parent, pivots[parent]);
            if !multiplier.is_finite() || !pivots[parent].is_finite() {
                return Ok(None);
            }
            order.push(Elimination {
                cell,
                parent,
                toward_parent,
                multiplier,
            });
            degree[cell] = 0;
            degree[parent] -= 1;
            incident[parent] ^= edge;
            if degree[parent] == 1 {
                leaves.push(parent);
            }
        }
        Ok((order.len() == n).then_some(Self { pivots, order }))
    }

    /// Apply without scratch allocation. False means cancellation or unsafe
    /// intermediate arithmetic; the caller can use diagonal scaling in the
    /// latter case. FGMRES still checks the original matrix's true residual.
    pub(super) fn apply(
        &self,
        residual: &[f64],
        output: &mut [f64],
        mut poll: impl FnMut() -> bool,
    ) -> bool {
        for (i, (out, &value)) in output.iter_mut().zip(residual).enumerate() {
            if i % 256 == 0 && !poll() {
                return false;
            }
            *out = value;
        }
        for (i, step) in self.order.iter().enumerate() {
            if i % 256 == 0 && !poll() {
                return false;
            }
            if step.parent != NONE {
                let value = (-step.multiplier).mul_add(output[step.cell], output[step.parent]);
                if !value.is_finite() {
                    return false;
                }
                output[step.parent] = value;
            }
        }
        for (i, step) in self.order.iter().rev().enumerate() {
            if i % 256 == 0 && !poll() {
                return false;
            }
            let value = if step.parent == NONE {
                output[step.cell]
            } else {
                (-step.toward_parent).mul_add(output[step.parent], output[step.cell])
            } / self.pivots[step.cell];
            if !value.is_finite() {
                return false;
            }
            output[step.cell] = value;
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stiff_unequal_capacity_tree_preserves_original_heat_balance() {
        let n = 49;
        let dt = 1000.0;
        let mass: Vec<_> = (0..n).map(|i| 0.5 + 0.25 * ((5 * i) % 11) as f64).collect();
        let target: Vec<_> = (0..n)
            .map(|i| 400.0 + (((7 * i) % 13) as f64 - 6.0) / dt)
            .collect();
        let mut links = Vec::new();
        for i in 1..n {
            links.push(((i - 1) / 2, i, 1.0 / (1.0 + 0.25 * (i % 5) as f64)));
            if i % 7 == 0 {
                links.push((i, (i - 1) / 2, 1.0 / 1.5));
            }
        }
        let mut initial = target.clone();
        for &(a, b, conductance) in &links {
            let heat = dt * conductance * (target[a] - target[b]);
            initial[a] += heat / mass[a];
            initial[b] -= heat / mass[b];
        }
        let mut diagonal = vec![1.0; n];
        let mut entries = Vec::new();
        let mut rhs = vec![0.0; n];
        for &(a, b, conductance) in &links {
            let dg = dt * conductance;
            diagonal[a] += dg / mass[a];
            diagonal[b] += dg / mass[b];
            entries.push((a, b, -dg / mass[b]));
            entries.push((b, a, -dg / mass[a]));
            let heat = dg * (initial[b] - initial[a]);
            rhs[a] += heat;
            rhs[b] -= heat;
        }
        let action = |x: &[f64]| {
            let mut result: Vec<_> = diagonal.iter().zip(x).map(|(d, x)| d * x).collect();
            for &(a, b, value) in &entries {
                result[a] += value * x[b];
            }
            result
        };
        let norm = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>().sqrt();
        // The best single Krylov direction with diagonal preconditioning is
        // alpha * D^-1 b. Its independently minimized true residual remains
        // large on this stiff tree; this is not a substituted solver backend.
        let diagonal_direction: Vec<_> = rhs.iter().zip(&diagonal).map(|(b, d)| b / d).collect();
        let applied = action(&diagonal_direction);
        let alpha = applied.iter().zip(&rhs).map(|(a, b)| a * b).sum::<f64>()
            / applied.iter().map(|a| a * a).sum::<f64>();
        let diagonal_residual: Vec<_> = applied
            .iter()
            .zip(&rhs)
            .map(|(a, b)| alpha * a - b)
            .collect();
        assert!(norm(&diagonal_residual) / norm(&rhs) > 1e-2);

        let forest = Forest::factor(&diagonal, &entries, &mut || Ok::<(), ()>(()))
            .unwrap()
            .unwrap();
        let mut energy = vec![0.0; n];
        assert!(forest.apply(&rhs, &mut energy, || true));
        for i in 0..n {
            assert!((initial[i] + energy[i] / mass[i] - target[i]).abs() < 2e-9);
        }
        let residual: Vec<_> = action(&energy)
            .iter()
            .zip(&rhs)
            .map(|(a, b)| a - b)
            .collect();
        assert!(norm(&residual) / norm(&rhs) < 1e-12);
        assert!(energy.iter().sum::<f64>().abs() < 1e-8);
    }

    #[test]
    fn nonsymmetric_parallel_and_one_sided_forest_solves_original_matrix() {
        let diagonal = [4.0, 3.0, 2.0, 5.0, 6.0, 7.0];
        let entries = [
            (0, 2, -1.5),
            (2, 0, -0.25),
            (0, 2, 0.5),
            (2, 0, -0.5),
            (2, 1, 0.0),
            (1, 2, -0.2),
            (4, 2, 2.0),
            (2, 4, -1.25),
            (4, 5, 0.75),
            (5, 4, -0.25),
            // A cancelled contribution is not an active edge of a cycle.
            (0, 1, 2.0),
            (0, 1, -2.0),
            (1, 0, 0.0),
        ];
        let expected = [2.0, -1.0, 3.0, -2.0, 0.5, 1.75];
        let mut rhs = [0.0; 6];
        for i in 0..6 {
            rhs[i] = diagonal[i] * expected[i];
        }
        for &(row, column, value) in &entries {
            rhs[row] += value * expected[column];
        }
        let forest = Forest::factor(&diagonal, &entries, &mut || Ok::<(), ()>(()))
            .unwrap()
            .unwrap();
        let mut solution = [0.0; 6];
        assert!(forest.apply(&rhs, &mut solution, || true));
        let mut residual = rhs;
        for i in 0..6 {
            assert!((solution[i] - expected[i]).abs() < 2e-14);
            residual[i] -= diagonal[i] * solution[i];
        }
        for &(row, column, value) in &entries {
            residual[row] -= value * solution[column];
        }
        assert!(residual.iter().all(|value| value.abs() < 2e-14));
    }

    #[test]
    fn cycles_unsafe_pivots_and_cancelled_factorization_do_not_claim_an_inverse() {
        assert!(
            Forest::factor(
                &[4.0; 3],
                &[(0, 1, -1.0), (1, 2, -1.0), (2, 0, -1.0)],
                &mut || Ok::<(), ()>(())
            )
            .unwrap()
            .is_none()
        );
        // This matrix is invertible but requires pivoting; forest elimination
        // must defer to FGMRES rather than label the physical system singular.
        assert!(
            Forest::factor(
                &[0.0, 2.0],
                &[(0, 1, 1.0), (1, 0, 1.0)],
                &mut || Ok::<(), ()>(())
            )
            .unwrap()
            .is_none()
        );
        let mut polls = 0;
        let result = Forest::factor(&[2.0; 3], &[(0, 1, -1.0), (1, 2, -1.0)], &mut || {
            polls += 1;
            if polls == 3 { Err(()) } else { Ok(()) }
        });
        assert!(matches!(result, Err(FactorError::Checkpoint(()))));
    }
}

//! The existing zero-fill ILU algorithm on flat canonical CSR storage.
//! Kept independent of matrix wrappers so the numerical kernel is testable
//! without the executor/runtime dependency graph.

/// Typed incomplete-factorization breakdown.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IluBreakdown {
    /// Row with a missing/zero diagonal or nonfinite factor arithmetic.
    /// A shifted retry is caller policy; this factor never perturbs the matrix.
    pub row: usize,
}

/// Incomplete-factorization failure with the caller's original checkpoint error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IluControlledError<E> {
    /// A missing/zero pivot or non-finite coefficient/intermediate.
    Breakdown(IluBreakdown),
    /// Sparse storage could not be reserved.
    Capacity,
    /// Cancellation or another caller-owned work limit.
    Checkpoint(E),
}

/// ILU(0): incomplete LU with zero fill on the original canonical CSR pattern.
/// Sequential elimination retains ascending row/column arithmetic order.
#[derive(Debug, Clone)]
pub struct Ilu0 {
    rows: Vec<usize>,
    columns: Vec<usize>,
    values: Vec<f64>,
    diagonal: Vec<usize>,
}

fn tick<E>(work: &mut usize, checkpoint: &mut impl FnMut() -> Result<(), E>) -> Result<(), E> {
    *work += 1;
    if *work == 256 {
        *work = 0;
        checkpoint()?;
    }
    Ok(())
}

fn reserved<T, E>(capacity: usize) -> Result<Vec<T>, IluControlledError<E>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| IluControlledError::Capacity)?;
    Ok(result)
}

fn copy<E, T: Copy>(
    source: &[T],
    work: &mut usize,
    checkpoint: &mut impl FnMut() -> Result<(), E>,
) -> Result<Vec<T>, IluControlledError<E>> {
    let mut result = reserved(source.len())?;
    for &value in source {
        tick(work, checkpoint).map_err(IluControlledError::Checkpoint)?;
        result.push(value);
    }
    Ok(result)
}

/// The public Csr adapter guarantees canonical square input.
pub(super) fn factor<E>(
    rows: &[usize],
    columns: &[usize],
    source_values: &[f64],
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<Ilu0, IluControlledError<E>> {
    checkpoint().map_err(IluControlledError::Checkpoint)?;
    let n = rows.len() - 1;
    let mut work = 0;
    let rows = copy(rows, &mut work, &mut checkpoint)?;
    let columns = copy(columns, &mut work, &mut checkpoint)?;
    let mut values = reserved(source_values.len())?;
    let mut diagonal = reserved(n)?;
    let mut position = reserved(n)?;
    // Validate every source coefficient before elimination, retaining the
    // controlled API's original offending-row detail and failure precedence.
    for i in 0..n {
        tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
        position.push(usize::MAX);
        let mut pivot = usize::MAX;
        for p in rows[i]..rows[i + 1] {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            let value = source_values[p];
            if !value.is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
            }
            if columns[p] == i {
                pivot = p;
            }
            values.push(value);
        }
        diagonal.push(pivot);
    }
    for i in 0..n {
        tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
        let start = rows[i];
        let end = rows[i + 1];
        for p in start..end {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            position[columns[p]] = p;
        }
        // Same ascending k and j order as the original row-vector ILU(0).
        // A position belongs to this row iff it lies in [start, end), so no
        // per-pivot row clone, binary search or marker-array reset is needed.
        for p in start..end {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            let k = columns[p];
            if k >= i {
                break;
            }
            let dk = diagonal[k];
            if dk == usize::MAX || values[dk] == 0.0 || !values[dk].is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row: k }));
            }
            let multiplier = values[p] / values[dk];
            if !multiplier.is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
            }
            values[p] = multiplier;
            for q in dk + 1..rows[k + 1] {
                tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
                let target = position[columns[q]];
                if (start..end).contains(&target) {
                    let updated = (-multiplier).mul_add(values[q], values[target]);
                    if !updated.is_finite() {
                        return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
                    }
                    values[target] = updated;
                }
            }
        }
        let di = diagonal[i];
        if di == usize::MAX || values[di] == 0.0 || !values[di].is_finite() {
            return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
        }
    }
    Ok(Ilu0 {
        rows,
        columns,
        values,
        diagonal,
    })
}

impl Ilu0 {
    /// Apply U^-1 L^-1 with checkpoints at bounded numerical-work boundaries.
    /// No scratch allocation or matrix symmetry assumption is needed.
    ///
    /// Numerical breakdown retains its row; checkpoint errors are unchanged.
    /// On either failure output is scratch and must be discarded or overwritten.
    /// The factor remains unchanged.
    /// Input and output lengths must equal the original matrix dimension.
    pub fn apply_controlled<E>(
        &self,
        residual: &[f64],
        output: &mut [f64],
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), IluControlledError<E>> {
        let n = self.diagonal.len();
        assert_eq!(residual.len(), n, "ilu0 residual length");
        assert_eq!(output.len(), n, "ilu0 output length");
        checkpoint().map_err(IluControlledError::Checkpoint)?;
        let mut work = 0;
        for (row, (out, &value)) in output.iter_mut().zip(residual).enumerate() {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            if !value.is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row }));
            }
            *out = value;
        }
        // Retain the established controlled traversal and checkpoint cadence,
        // including diagonal visits and the skipped lower entries in U.
        for i in 0..n {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            let mut value = output[i];
            for p in self.rows[i]..self.rows[i + 1] {
                tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
                let column = self.columns[p];
                if column >= i {
                    break;
                }
                value = (-self.values[p]).mul_add(output[column], value);
            }
            if !value.is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
            }
            output[i] = value;
        }
        for i in (0..n).rev() {
            tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
            let mut value = output[i];
            let mut diagonal = 1.0;
            for p in self.rows[i]..self.rows[i + 1] {
                tick(&mut work, &mut checkpoint).map_err(IluControlledError::Checkpoint)?;
                let column = self.columns[p];
                if column < i {
                    continue;
                }
                if column == i {
                    diagonal = self.values[p];
                } else {
                    value = (-self.values[p]).mul_add(output[column], value);
                }
            }
            value /= diagonal;
            output[i] = value;
            if !value.is_finite() {
                return Err(IluControlledError::Breakdown(IluBreakdown { row: i }));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonsymmetric_cycles_and_negative_pivots_keep_the_original_matrix() {
        let rows = [0, 3, 6, 9, 10];
        let columns = [0, 1, 2, 0, 1, 2, 0, 1, 2, 3];
        let values = [4.0, -1.0, 0.5, -2.0, 5.0, -1.0, 0.25, -1.0, 3.0, -2.0];
        let rhs = [10.5, -12.0, 10.5, -8.0];
        let expected = [2.0, -1.0, 3.0, 4.0];
        let ilu = factor(&rows, &columns, &values, || Ok::<(), ()>(())).unwrap();
        let mut solution = [0.0; 4];
        ilu.apply_controlled(&rhs, &mut solution, || Ok::<(), ()>(()))
            .unwrap();
        for i in 0..4 {
            assert!((solution[i] - expected[i]).abs() < 2e-14);
            let applied: f64 = (rows[i]..rows[i + 1])
                .map(|p| values[p] * solution[columns[p]])
                .sum();
            assert!((applied - rhs[i]).abs() < 2e-14);
        }
        assert_eq!(ilu.rows, rows);
        assert_eq!(ilu.columns, columns);
        assert_eq!(ilu.values.len(), values.len());
    }

    #[test]
    fn stiff_cyclic_heat_networks_improve_the_original_true_residual() {
        for grid in [false, true] {
            let n = 49;
            let dt = 1000.0;
            let mass: Vec<_> = (0..n).map(|i| 0.5 + 0.25 * ((5 * i) % 11) as f64).collect();
            let target: Vec<_> = (0..n)
                .map(|i| 400.0 + (((7 * i) % 13) as f64 - 6.0) / dt)
                .collect();
            let mut links = Vec::new();
            for i in 0..n {
                if grid {
                    if i % 7 < 6 {
                        links.push((i, i + 1, 1.0 / (1.0 + 0.25 * (i % 5) as f64)));
                    }
                    if i + 7 < n {
                        links.push((i, i + 7, 1.0 / (1.0 + 0.25 * ((i + 2) % 5) as f64)));
                    }
                } else {
                    links.push((i, (i + 1) % n, 1.0 / (1.0 + 0.25 * (i % 5) as f64)));
                }
                if i % 7 == 0 {
                    links.push(((i + 1) % n, i, 1.0 / 1.5));
                }
            }
            let mut initial = target.clone();
            for &(a, b, conductance) in &links {
                let heat = dt * conductance * (target[a] - target[b]);
                initial[a] += heat / mass[a];
                initial[b] -= heat / mass[b];
            }
            assert!(initial.iter().all(|h| (100.0..1500.0).contains(h)));
            // Independent dense test oracle for energy-coordinate heat
            // balances. Production storage and factors remain sparse.
            let mut matrix = vec![vec![0.0; n]; n];
            let mut rhs = vec![0.0; n];
            for (i, row) in matrix.iter_mut().enumerate() {
                row[i] = 1.0;
            }
            for &(a, b, conductance) in &links {
                let dg = dt * conductance;
                matrix[a][a] += dg / mass[a];
                matrix[b][b] += dg / mass[b];
                matrix[a][b] -= dg / mass[b];
                matrix[b][a] -= dg / mass[a];
                let heat = dg * (initial[b] - initial[a]);
                rhs[a] += heat;
                rhs[b] -= heat;
            }
            for column in 0..n {
                let sum: f64 = matrix.iter().map(|row| row[column]).sum();
                assert!((sum - 1.0).abs() < 2e-11);
            }
            let mut rows = vec![0];
            let mut columns = Vec::new();
            let mut values = Vec::new();
            for row in &matrix {
                for (column, &value) in row.iter().enumerate() {
                    if value != 0.0 {
                        columns.push(column);
                        values.push(value);
                    }
                }
                rows.push(values.len());
            }
            let ilu = factor(&rows, &columns, &values, || Ok::<(), ()>(())).unwrap();
            let mut correction = vec![0.0; n];
            ilu.apply_controlled(&rhs, &mut correction, || Ok::<(), ()>(()))
                .unwrap();
            let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
            // Independent minimization over the first preconditioned Krylov
            // direction. This measures the real A, not the incomplete factors.
            let relative_residual = |direction: &[f64]| {
                let applied: Vec<_> = matrix.iter().map(|row| dot(row, direction)).collect();
                let alpha = dot(&applied, &rhs) / dot(&applied, &applied);
                let residual: Vec<_> = applied
                    .iter()
                    .zip(&rhs)
                    .map(|(a, b)| alpha * a - b)
                    .collect();
                (dot(&residual, &residual) / dot(&rhs, &rhs)).sqrt()
            };
            let diagonal: Vec<_> = (0..n).map(|i| rhs[i] / matrix[i][i]).collect();
            let before = relative_residual(&diagonal);
            let after = relative_residual(&correction);
            assert!(before > 0.05);
            assert!(after * 2.0 < before, "grid={grid}: {before} -> {after}");
            assert_eq!(ilu.values.len(), values.len());
            eprintln!(
                "cyclic heat network grid={grid}: diagonal residual {before:.9e}, ILU(0) residual {after:.9e}"
            );
        }
    }

    #[test]
    fn unsafe_factor_and_substitution_arithmetic_refuse_without_shifts() {
        for value in [0.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                factor(&[0, 1], &[0], &[value], || Ok::<(), ()>(())),
                Err(IluControlledError::Breakdown(IluBreakdown { row: 0 }))
            ));
        }
        assert!(matches!(
            factor(&[0, 1, 2], &[0, 0], &[1.0, 1.0], || Ok::<(), ()>(())),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 1 }))
        ));
        assert!(matches!(
            factor(&[0, 2, 4], &[0, 1, 0, 1], &[1.0; 4], || Ok::<(), ()>(())),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 1 }))
        ));
        assert!(matches!(
            factor(
                &[0, 2, 4],
                &[0, 1, 0, 1],
                &[1e-300, 1e300, 1e300, 1.0],
                || Ok::<(), ()>(())
            ),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 1 }))
        ));
        let tiny = f64::from_bits(1);
        let ilu = factor(&[0, 1], &[0], &[tiny], || Ok::<(), ()>(())).unwrap();
        let mut output = [0.0];
        assert_eq!(
            ilu.apply_controlled(&[1.0], &mut output, || Ok::<(), ()>(())),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 0 }))
        );
        ilu.apply_controlled(&[tiny], &mut output, || Ok::<(), ()>(()))
            .unwrap();
        assert_eq!(output, [1.0]);
    }

    #[test]
    fn cancelled_setup_and_substitution_preserve_the_original_error_and_factor() {
        let n = 257;
        let mut rows = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            if i > 0 {
                columns.push(i - 1);
                values.push(-1.0);
            }
            columns.push(i);
            values.push(3.0);
            if i + 1 < n {
                columns.push(i + 1);
                values.push(-0.5);
            }
            rows.push(values.len());
        }
        let mut setup_polls = 0;
        let ilu = factor(&rows, &columns, &values, || {
            setup_polls += 1;
            Ok::<(), usize>(())
        })
        .unwrap();
        for stop in [1, setup_polls - 1] {
            let mut polls = 0;
            let result = factor(&rows, &columns, &values, || {
                polls += 1;
                if polls == stop { Err(stop) } else { Ok(()) }
            });
            assert!(matches!(result, Err(IluControlledError::Checkpoint(value)) if value == stop));
        }
        let rhs = vec![1.0; n];
        let mut expected = vec![0.0; n];
        let mut apply_polls = 0;
        ilu.apply_controlled(&rhs, &mut expected, || {
            apply_polls += 1;
            Ok::<(), usize>(())
        })
        .unwrap();
        let mut output = vec![0.0; n];
        for stop in [1, apply_polls - 1] {
            let mut polls = 0;
            assert_eq!(
                ilu.apply_controlled(&rhs, &mut output, || {
                    polls += 1;
                    if polls == stop { Err(stop) } else { Ok(()) }
                }),
                Err(IluControlledError::Checkpoint(stop))
            );
        }
        ilu.apply_controlled(&rhs, &mut output, || Ok::<(), ()>(()))
            .unwrap();
        assert_eq!(output, expected);
    }

    #[test]
    fn controlled_capacity_row_details_and_long_row_cancellation_are_preserved() {
        assert!(matches!(
            reserved::<f64, ()>(usize::MAX),
            Err(IluControlledError::Capacity)
        ));
        // All source values are validated before eliminating an earlier
        // zero pivot, preserving the established offending-source row.
        assert!(matches!(
            factor(&[0, 1, 2], &[0, 1], &[0.0, f64::NAN], || Ok::<(), ()>(())),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 1 }))
        ));
        let n = 700;
        let mut rows = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        for row in 0..n {
            if row == 0 {
                for column in 0..n {
                    columns.push(column);
                    values.push(if column == 0 { 2.0 } else { -0.001 });
                }
            } else {
                columns.push(row);
                values.push(2.0);
            }
            rows.push(values.len());
        }
        let ilu = factor(&rows, &columns, &values, || Ok::<(), ()>(())).unwrap();
        let rhs = vec![1.0; n];
        let mut output = vec![0.0; n];
        let mut polls = 0;
        // The current public API regression cancels at this checkpoint,
        // inside the final high-degree row of the backward substitution.
        assert_eq!(
            ilu.apply_controlled(&rhs, &mut output, || {
                polls += 1;
                if polls == 17 {
                    Err("cancel solve")
                } else {
                    Ok(())
                }
            }),
            Err(IluControlledError::Checkpoint("cancel solve"))
        );
        let mut invalid = rhs;
        invalid[13] = f64::NAN;
        assert_eq!(
            ilu.apply_controlled(&invalid, &mut output, || Ok::<(), ()>(())),
            Err(IluControlledError::Breakdown(IluBreakdown { row: 13 }))
        );
    }
}

//! Sparse identification of ideal equipotential coordinates.
//!
//! This is the original difference-form operator under u = P q, not a
//! terminal-resistor fit. Contributions internal to one group vanish before
//! summation. Reconstructing each diagonal from external couplings avoids
//! subtracting large internal conductances to recover a small physical link.

#[derive(Debug, PartialEq)]
pub(super) enum GroupedError {
    Invalid(&'static str),
    Cancelled,
}

pub(super) struct GroupedMatrix {
    pub row_ptr: Vec<usize>,
    pub columns: Vec<usize>,
    pub values: Vec<f64>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GroupedError> {
    if cancelled() {
        Err(GroupedError::Cancelled)
    } else {
        Ok(())
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), GroupedError> {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    let remainder = (next - *sum) - adjusted;
    if !value.is_finite() || !adjusted.is_finite() || !next.is_finite() || !remainder.is_finite() {
        return Err(GroupedError::Invalid(
            "grouped potential arithmetic is not representable",
        ));
    }
    *sum = next;
    *correction = remainder;
    Ok(())
}

/// Deterministic sparse P^T K P for an already admitted constant-nullspace
/// operator. Storage is O(vertices + original nonzeros); only touched columns
/// are sorted per grouped row. Original vertex and row order fix reductions.
pub(super) fn condense<'a>(
    groups: &[usize],
    count: usize,
    row: impl Fn(usize) -> (&'a [usize], &'a [f64]),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GroupedMatrix, GroupedError> {
    poll(cancelled)?;
    if count == 0 || count > groups.len() || groups.iter().any(|&g| g >= count) {
        return Err(GroupedError::Invalid(
            "complete bounded potential groups required",
        ));
    }
    let mut starts = vec![0usize; count + 1];
    for &group in groups {
        poll(cancelled)?;
        starts[group + 1] += 1;
    }
    if starts[1..].contains(&0) {
        return Err(GroupedError::Invalid(
            "potential groups must all be populated",
        ));
    }
    for group in 0..count {
        starts[group + 1] += starts[group];
    }
    let mut cursor = starts[..count].to_vec();
    let mut vertices = vec![0usize; groups.len()];
    for (vertex, &group) in groups.iter().enumerate() {
        poll(cancelled)?;
        vertices[cursor[group]] = vertex;
        cursor[group] += 1;
    }
    let mut sums = vec![0.0; count];
    let mut corrections = vec![0.0; count];
    let mut marks = vec![usize::MAX; count];
    let mut touched = Vec::new();
    let mut result = GroupedMatrix {
        row_ptr: vec![0],
        columns: Vec::new(),
        values: Vec::new(),
    };
    for group in 0..count {
        poll(cancelled)?;
        touched.clear();
        for &vertex in &vertices[starts[group]..starts[group + 1]] {
            let (columns, values) = row(vertex);
            if columns.len() != values.len() {
                return Err(GroupedError::Invalid(
                    "matching sparse columns and values required",
                ));
            }
            for (&column, &value) in columns.iter().zip(values) {
                poll(cancelled)?;
                if column >= groups.len() || !value.is_finite() {
                    return Err(GroupedError::Invalid(
                        "finite in-range sparse entries required",
                    ));
                }
                let other = groups[column];
                if other == group || value == 0.0 {
                    continue;
                }
                if marks[other] != group {
                    marks[other] = group;
                    touched.push(other);
                    sums[other] = 0.0;
                    corrections[other] = 0.0;
                }
                add(&mut sums[other], &mut corrections[other], value)?;
            }
        }
        // Include the diagonal in canonical CSR order, including an isolated
        // prescribed coordinate. Its zero diagonal is never inverted by CG.
        touched.push(group);
        touched.sort_unstable();
        let (mut diagonal, mut correction) = (0.0, 0.0);
        for &other in &touched {
            poll(cancelled)?;
            if other != group {
                add(&mut diagonal, &mut correction, -sums[other])?;
            }
        }
        for &other in &touched {
            poll(cancelled)?;
            let value = if other == group {
                diagonal
            } else {
                sums[other]
            };
            if other == group || value != 0.0 {
                result.columns.push(other);
                result.values.push(value);
            }
        }
        result.row_ptr.push(result.columns.len());
    }
    poll(cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense(result: &GroupedMatrix) -> Vec<Vec<f64>> {
        let n = result.row_ptr.len() - 1;
        let mut matrix = vec![vec![0.0; n]; n];
        for (i, row) in matrix.iter_mut().enumerate() {
            for k in result.row_ptr[i]..result.row_ptr[i + 1] {
                row[result.columns[k]] = result.values[k];
            }
        }
        matrix
    }

    // G0/G1: the shared ideal middle electrode joins two physical pieces.
    // Bulk links 1 S and 1/3 S at q=[0,2,8] carry [-2,0,2] A and 16 W.
    #[test]
    fn floating_join_preserves_reactions_and_energy() {
        let columns = [vec![0, 1], vec![0, 1], vec![2, 3], vec![2, 3]];
        let values = [
            vec![1.0, -1.0],
            vec![-1.0, 1.0],
            vec![1.0 / 3.0, -1.0 / 3.0],
            vec![-1.0 / 3.0, 1.0 / 3.0],
        ];
        let result = condense(&[0, 1, 1, 2], 3, |v| (&columns[v], &values[v]), &mut || {
            false
        })
        .unwrap();
        let matrix = dense(&result);
        let q = [0.0, 2.0, 8.0];
        let current: Vec<f64> = matrix
            .iter()
            .map(|r| r.iter().zip(q).map(|(a, b)| a * b).sum())
            .collect();
        for (actual, expected) in current.iter().zip([-2.0, 0.0, 2.0]) {
            assert!((actual - expected).abs() < 1e-14);
        }
        let power: f64 = current.iter().zip(q).map(|(a, b)| a * b).sum();
        assert!((power - 16.0).abs() < 1e-14);
        assert_eq!(matrix[0][1], matrix[1][0]);
        assert_eq!(matrix[1][2], matrix[2][1]);
    }

    // G1: identifying an ideal electrode removes its internal 1e30 S link
    // exactly; a remaining 3 S physical connection must not disappear.
    #[test]
    fn strong_internal_links_do_not_erase_small_external_conductance() {
        let columns = [vec![0, 1, 2], vec![0, 1, 2], vec![0, 1, 2]];
        let values = [
            vec![1e30, -1e30, -1.0],
            vec![-1e30, 1e30, -2.0],
            vec![-1.0, -2.0, 3.0],
        ];
        let result = condense(&[0, 0, 1], 2, |v| (&columns[v], &values[v]), &mut || false).unwrap();
        assert_eq!(dense(&result), vec![vec![3.0, -3.0], vec![-3.0, 3.0]]);
    }

    // G0/G3: signed off-diagonal entries occur in anisotropic/obtuse elements.
    // Grouping must preserve their signed contribution, not clamp to a graph.
    #[test]
    fn signed_couplings_and_nontrivial_group_order_are_preserved() {
        let columns = [vec![0, 1, 2], vec![0, 1, 2], vec![0, 1, 2]];
        let values = [
            vec![3.0, -4.0, 1.0],
            vec![-4.0, 6.0, -2.0],
            vec![1.0, -2.0, 1.0],
        ];
        let result = condense(&[1, 0, 1], 2, |v| (&columns[v], &values[v]), &mut || false).unwrap();
        assert_eq!(dense(&result), vec![vec![6.0, -6.0], vec![-6.0, 6.0]]);
        let all = condense(&[0, 0, 0], 1, |v| (&columns[v], &values[v]), &mut || false).unwrap();
        assert_eq!(dense(&all), vec![vec![0.0]]);
    }

    // G4: stop inside assembly, then repeat the original inputs bit-for-bit.
    #[test]
    fn cancellation_does_not_change_retry_or_sparse_inputs() {
        let columns = [vec![0, 1], vec![0, 1]];
        let values = [vec![2.0, -2.0], vec![-2.0, 2.0]];
        let run = || condense(&[0, 1], 2, |v| (&columns[v], &values[v]), &mut || false).unwrap();
        let expected = run();
        let mut calls = 0;
        let cancelled = condense(&[0, 1], 2, |v| (&columns[v], &values[v]), &mut || {
            calls += 1;
            calls >= 10
        });
        assert!(matches!(cancelled, Err(GroupedError::Cancelled)));
        let retry = run();
        assert_eq!(retry.row_ptr, expected.row_ptr);
        assert_eq!(retry.columns, expected.columns);
        assert_eq!(retry.values, expected.values);
        assert_eq!(values[0], vec![2.0, -2.0]);
    }

    #[test]
    fn malformed_and_nonrepresentable_operators_refuse() {
        let columns = [vec![0, 1], vec![0, 1]];
        let values = [vec![1.0, -1.0], vec![-1.0, 1.0]];
        assert!(condense(&[0, 2], 2, |v| (&columns[v], &values[v]), &mut || false).is_err());
        assert!(condense(&[0, 0], 2, |v| (&columns[v], &values[v]), &mut || false).is_err());
        let bad = [vec![f64::NAN, -1.0], vec![-1.0, 1.0]];
        assert!(condense(&[0, 1], 2, |v| (&columns[v], &bad[v]), &mut || false).is_err());
        let columns = [vec![2], vec![2], vec![0, 1]];
        let values = [vec![-f64::MAX], vec![-f64::MAX], vec![-f64::MAX, -f64::MAX]];
        assert!(condense(&[0, 0, 1], 2, |v| (&columns[v], &values[v]), &mut || false).is_err());
    }
}

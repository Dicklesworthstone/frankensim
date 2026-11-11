//! Dense linear solve for a block's interface Jacobian.
//!
//! Interface systems of lumped blocks are small (one unknown per scalar port
//! side in a loop), so the driver uses plain Gaussian elimination with partial
//! pivoting. Large or sparse blocks belong to the planner and to fs-solver.
//!
//! The factors are kept apart from the solve, so a block whose Jacobian does
//! not change from step to step is factored once.

/// Factor a row-major `n x n` matrix in place as `P A = L U` with partial
/// pivoting. `a` then holds `U` on and above the diagonal and the multipliers
/// of the unit lower triangle `L` below it; `pivots[k]` is the row exchanged
/// with row `k`.
///
/// Returns `false` if a pivot is zero or not finite; `a` is then partly
/// eliminated and must not be used.
pub(crate) fn factor_dense(a: &mut [f64], pivots: &mut [usize], n: usize) -> bool {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(pivots.len(), n);
    for column in 0..n {
        let mut pivot_row = column;
        let mut pivot = a[column * n + column].abs();
        for row in column + 1..n {
            let candidate = a[row * n + column].abs();
            if candidate > pivot {
                pivot = candidate;
                pivot_row = row;
            }
        }
        if !(pivot.is_finite() && pivot > 0.0) {
            return false;
        }
        pivots[column] = pivot_row;
        if pivot_row != column {
            // Whole rows, so the multipliers already stored move with them.
            for k in 0..n {
                a.swap(column * n + k, pivot_row * n + k);
            }
        }
        let diagonal = a[column * n + column];
        for row in column + 1..n {
            let factor = a[row * n + column] / diagonal;
            a[row * n + column] = factor;
            if factor != 0.0 {
                for k in column + 1..n {
                    a[row * n + k] -= factor * a[column * n + k];
                }
            }
        }
    }
    true
}

/// Solve `A x = b` with the factors of [`factor_dense`]; `b` becomes `x`.
///
/// Returns `false` if the solution is not finite.
pub(crate) fn solve_factored(a: &[f64], pivots: &[usize], b: &mut [f64], n: usize) -> bool {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(pivots.len(), n);
    debug_assert_eq!(b.len(), n);
    for (row, &exchanged) in pivots.iter().enumerate() {
        b.swap(row, exchanged);
    }
    for column in 0..n {
        for row in column + 1..n {
            let factor = a[row * n + column];
            if factor != 0.0 {
                b[row] -= factor * b[column];
            }
        }
    }
    for row in (0..n).rev() {
        let mut sum = b[row];
        for k in row + 1..n {
            sum -= a[row * n + k] * b[k];
        }
        b[row] = sum / a[row * n + row];
    }
    b.iter().all(|value| value.is_finite())
}

#[cfg(test)]
mod tests {
    use super::{factor_dense, solve_factored};

    /// Factor and solve in one call; `a` is replaced by its factors.
    fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
        let mut pivots = vec![0; n];
        factor_dense(a, &mut pivots, n) && solve_factored(a, &pivots, b, n)
    }

    /// Elimination applied to the right-hand side as it goes, with no kept
    /// factors: the solve this module had before the factors were kept.
    fn eliminate(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
        for column in 0..n {
            let mut pivot_row = column;
            let mut pivot = a[column * n + column].abs();
            for row in column + 1..n {
                let candidate = a[row * n + column].abs();
                if candidate > pivot {
                    pivot = candidate;
                    pivot_row = row;
                }
            }
            if !(pivot.is_finite() && pivot > 0.0) {
                return false;
            }
            if pivot_row != column {
                for k in 0..n {
                    a.swap(column * n + k, pivot_row * n + k);
                }
                b.swap(column, pivot_row);
            }
            let diagonal = a[column * n + column];
            for row in column + 1..n {
                let factor = a[row * n + column] / diagonal;
                if factor != 0.0 {
                    for k in column..n {
                        a[row * n + k] -= factor * a[column * n + k];
                    }
                    b[row] -= factor * b[column];
                }
            }
        }
        for row in (0..n).rev() {
            let mut sum = b[row];
            for k in row + 1..n {
                sum -= a[row * n + k] * b[k];
            }
            b[row] = sum / a[row * n + row];
        }
        b.iter().all(|value| value.is_finite())
    }

    #[test]
    fn solves_a_system_that_needs_pivoting() {
        // First pivot is zero without a row exchange.
        let mut a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0];
        // b = a * [1, 2, 1].
        let mut b = vec![5.0, 3.0, 4.0];
        assert!(solve_dense(&mut a, &mut b, 3));
        let expected = [1.0, 2.0, 1.0];
        for (got, want) in b.iter().zip(expected) {
            assert!((got - want).abs() <= 1.0e-14, "{got} vs {want}");
        }
    }

    #[test]
    fn refuses_a_singular_matrix() {
        let mut a = vec![1.0, 2.0, 2.0, 4.0];
        let mut b = vec![1.0, 2.0];
        assert!(!solve_dense(&mut a, &mut b, 2));
    }

    #[test]
    fn an_empty_system_is_solved() {
        assert!(solve_dense(&mut [], &mut [], 0));
    }

    #[test]
    fn kept_factors_solve_every_right_hand_side_to_the_bit_of_one_elimination() {
        // Matrices that exchange rows at several columns, with zeros that
        // make some multipliers vanish. Deterministic pseudo-random entries.
        let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        };
        for n in [1_usize, 2, 5, 12] {
            let matrix: Vec<f64> = (0..n * n)
                .map(|entry| {
                    let value = next();
                    // A seventh of the entries are zero (never a whole
                    // column at these sizes), the diagonal small.
                    if entry % 7 == 3 {
                        0.0
                    } else if entry % (n + 1) == 0 {
                        1.0e-3 * value
                    } else {
                        value
                    }
                })
                .collect();
            let mut factors = matrix.clone();
            let mut pivots = vec![0; n];
            assert!(factor_dense(&mut factors, &mut pivots, n));
            assert!(n < 5 || pivots.iter().enumerate().any(|(row, &to)| row != to));
            for _ in 0..3 {
                let right: Vec<f64> = (0..n).map(|_| next()).collect();
                let mut kept = right.clone();
                assert!(solve_factored(&factors, &pivots, &mut kept, n));
                let mut once = right.clone();
                assert!(eliminate(&mut matrix.clone(), &mut once, n));
                for (a, b) in kept.iter().zip(&once) {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
                // And it is a solution.
                for row in 0..n {
                    let product: f64 = (0..n).map(|k| matrix[row * n + k] * kept[k]).sum();
                    assert!((product - right[row]).abs() <= 1.0e-9, "{product} vs {}", right[row]);
                }
            }
        }
    }
}

//! Dense linear solve for a block's interface Jacobian.
//!
//! Interface systems of lumped blocks are small (one unknown per scalar port
//! side in a loop), so the driver uses plain Gaussian elimination with partial
//! pivoting. Large or sparse blocks belong to the planner and to fs-solver.

/// Solve `a x = b` in place for a row-major `n x n` matrix; `b` becomes `x`.
///
/// Returns `false` if a pivot is zero or not finite, or if the solution is
/// not finite. `a` is destroyed.
pub(crate) fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(b.len(), n);
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

#[cfg(test)]
mod tests {
    use super::solve_dense;

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
}

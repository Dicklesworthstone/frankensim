//! Complex tridiagonal solve with the ORIGINAL operation graph of the dense owner.
//!
//! The dense elimination ([`super::lu_complex_with_cancel`]) updates element (i, c) once per earlier pivot k, in
//! ascending k, as `difference(a, product(f_ik, u_kc))`, picks the pivot of column k by magnitude (lowest index on
//! ties), swaps WHOLE rows (earlier multipliers travel with the row) and divides by the pivot. On a tridiagonal
//! matrix almost every one of those updates has a structural zero as `f` or as `u`, and `x - (+-0)` returns `x`
//! exactly, so only these survive: at step k the candidates for the pivot are the rows at positions k and k + 1
//! (every lower row has a zero in column k), the pivot row has at most the three entries (k, k), (k, k + 1),
//! (k, k + 2) (a carried row has (k, k) and (k, k + 1), the untouched row k + 1 has all three), and exactly one
//! row, at position k + 1, is eliminated, getting the multiplier `f = b0 / a0`. This routine performs just those
//! operations, in the same order, with the same operators; the multiplier is applied to every right-hand side at
//! once, at the step it is formed, which is the dense forward substitution in the same ascending-k order for the
//! row that carries it (rows swapped at a later step carry their updated right-hand side along, as the dense
//! solve's permuted copy does). The back substitution is the dense one on the at most two nonzero entries to the
//! right of the diagonal. Every nonzero value is therefore bit-identical to `lu_complex(..).solve(..)`; only the
//! sign of a zero can differ, which no later nonzero value or magnitude comparison reads.
//!
//! Cost: O(n) per matrix and O(n) per right-hand side, against O(n^3) and O(n^2); no n x n storage.
use super::{C64, LuComplexError, difference, product};

/// Solve `T X = B` in place for the complex tridiagonal matrix `T` with sub-diagonal `lower` (`T[i + 1][i]`),
/// diagonal `diag` and super-diagonal `upper` (`T[i][i + 1]`), and `columns` right-hand sides stored one after the
/// other in `rhs` (column `j` is `rhs[j * n..(j + 1) * n]`, `n = diag.len()`).
///
/// Same magnitude pivots (lowest index on ties), whole-row swaps, singular index and scalar operation order as
/// [`super::lu_complex_with_cancel`] followed by its `solve`, so every nonzero result is bit-identical to a dense
/// factorization of the same matrix with its structural zeros stored; see the module header. `rhs` is unchanged
/// when the matrix is singular (the error names the elimination step, as the dense error does).
///
/// # Errors
/// [`LuComplexError::Singular`] when the magnitude pivot of a column is zero.
/// # Panics
/// If `lower` and `upper` do not have `n - 1` entries or `rhs` does not have `n * columns` entries.
pub fn solve_tridiagonal_complex(
    lower: &[C64], diag: &[C64], upper: &[C64], rhs: &mut [C64], columns: usize,
) -> Result<(), LuComplexError> {
    let n = diag.len();
    assert_eq!(lower.len(), n.saturating_sub(1), "lower must have n - 1 entries");
    assert_eq!(upper.len(), n.saturating_sub(1), "upper must have n - 1 entries");
    assert_eq!(rhs.len(), n * columns, "rhs must hold n * columns entries");
    if n == 0 {
        return Ok(());
    }
    // The pivot row of every step: entries (k, k), (k, k + 1), (k, k + 2).
    let mut pivot_rows = vec![[C64::ZERO; 3]; n];
    // Right-hand sides are updated in place only after the elimination cannot fail, so stage the swaps and
    // multipliers: (swap at step k, multiplier of step k). Singular matrices then leave `rhs` untouched.
    let mut steps: Vec<(bool, C64)> = Vec::with_capacity(n);
    // The row at position k: columns k, k + 1, k + 2.
    let mut a = [diag[0], if n > 1 { upper[0] } else { C64::ZERO }, C64::ZERO];
    for k in 0..n {
        if k + 1 == n {
            if a[0].abs() == 0.0 {
                return Err(LuComplexError::Singular { index: k });
            }
            pivot_rows[k] = a;
            break;
        }
        // The untouched row k + 1.
        let mut b = [lower[k], diag[k + 1], if k + 2 < n { upper[k + 1] } else { C64::ZERO }];
        let mut best = a[0].abs();
        let candidate = b[0].abs();
        let swap = candidate > best;
        if swap {
            best = candidate;
        }
        if best == 0.0 {
            return Err(LuComplexError::Singular { index: k });
        }
        if swap {
            core::mem::swap(&mut a, &mut b);
        }
        let f = b[0] / a[0];
        pivot_rows[k] = a;
        steps.push((swap, f));
        a = [difference(b[1], product(f, a[1])), difference(b[2], product(f, a[2])), C64::ZERO];
    }
    for column in rhs.chunks_exact_mut(n) {
        for (k, &(swap, f)) in steps.iter().enumerate() {
            if swap {
                column.swap(k, k + 1);
            }
            let pivot_value = column[k];
            column[k + 1] = difference(column[k + 1], product(f, pivot_value));
        }
        for i in (0..n).rev() {
            let row = &pivot_rows[i];
            let mut v = column[i];
            if i + 1 < n {
                v = difference(v, product(row[1], column[i + 1]));
            }
            if i + 2 < n {
                v = difference(v, product(row[2], column[i + 2]));
            }
            column[i] = v / row[0];
        }
    }
    Ok(())
}

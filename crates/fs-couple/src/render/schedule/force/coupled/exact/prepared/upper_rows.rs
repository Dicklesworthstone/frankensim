//! The quadratic form `z^T W z` of an exactly symmetric `W`, summed by rows over the upper triangle.
//!
//! Row `i` contributes `z_i (W_ii z_i + 2 s_i)`, where `s_i` is the sum of `W_ij z_j` over the columns `j > i`, taken in
//! column order. That is half the terms of the full square. Eight consecutive rows advance together, each with its own
//! sum, so one row's additions do not wait for another row's: a row first takes its columns inside its block of eight
//! and then the columns every row of the block shares. The rows are then added in row order. The magnitude is the same
//! sum over the absolute terms. The order is a function of the side alone and takes no thread count.
//!
//! The lower triangle is never read, so the form is `z^T W z` only for a `W` that equals its transpose:
//! [`bitwise_symmetric`] is the check a caller makes once before it relies on this order.
//!
//! This file names nothing of its crate. `tests/upper_row_forms.rs` compiles it a second time and holds it to a plain
//! loop of the same order, bit for bit.

/// Rows that advance together.
const BLOCK: usize = 8;

/// Whether the row-major square `matrix` of side `side` equals its transpose bit for bit.
pub(super) fn bitwise_symmetric(matrix: &[f64], side: usize) -> bool {
    matrix.len() == side * side
        && (0..side).all(|i| {
            (i + 1..side).all(|j| matrix[i * side + j].to_bits() == matrix[j * side + i].to_bits())
        })
}

/// The sum and the magnitude of `z^T W z` for the row-major symmetric `matrix` of side `z.len()`, in the order of the
/// module's description.
pub(super) fn upper_row_form(matrix: &[f64], z: &[f64]) -> (f64, f64) {
    let side = z.len();
    let (mut sum, mut magnitude) = (0.0, 0.0);
    let mut first = 0;
    while first < side {
        let last = (first + BLOCK).min(side);
        let (mut sums, mut sizes) = ([0.0_f64; BLOCK], [0.0_f64; BLOCK]);
        // The columns inside the block: row i takes i + 1 .. last.
        for (k, i) in (first..last).enumerate() {
            let row = &matrix[i * side + i + 1..i * side + last];
            for (value, right) in row.iter().zip(&z[i + 1..last]) {
                let term = value * right;
                sums[k] += term;
                sizes[k] += term.abs();
            }
        }
        // The columns every row of a whole block shares.
        if last - first == BLOCK {
            let rows: [&[f64]; BLOCK] = core::array::from_fn(|k| {
                &matrix[(first + k) * side + last..(first + k + 1) * side]
            });
            for (j, right) in z[last..].iter().enumerate() {
                for ((partial, size), row) in sums.iter_mut().zip(&mut sizes).zip(&rows) {
                    let term = row[j] * right;
                    *partial += term;
                    *size += term.abs();
                }
            }
        }
        for (k, i) in (first..last).enumerate() {
            let diagonal = matrix[i * side + i] * z[i];
            sum += z[i] * (diagonal + 2.0 * sums[k]);
            magnitude += z[i].abs() * (diagonal.abs() + 2.0 * sizes[k]);
        }
        first = last;
    }
    (sum, magnitude)
}

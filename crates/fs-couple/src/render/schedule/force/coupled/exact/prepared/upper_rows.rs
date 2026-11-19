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
//! [`row_storage`] is the other sum this module orders by rows: the storage of unit-mass oscillators under a stiffness
//! that need not equal its transpose, so it reads the whole square.
//!
//! [`identity_rows`] is a check of another kind: whether the last rows of a square are the identity's, so that a
//! product by that square need not compute them.
//!
//! [`product_from`] is a product of whole rows by `[state | force]` that keeps each row's sum over the state's columns
//! ([`StateSums`]): the next product on the same state adds only the force's columns, the same additions in the same
//! order as one left-to-right sum.
//!
//! [`upper_row_form_from`] is the form of `[state | force]` that keeps the same kind of sums ([`StateForms`]): a row
//! takes the state's columns first, so the next form on the same state continues every row over the force's columns.
//!
//! This file names nothing of its crate. `tests/upper_row_forms.rs` compiles it a second time and holds each sum to a
//! plain loop of the same order, bit for bit.

/// Rows that advance together.
const BLOCK: usize = 8;

/// Whether the row-major square `matrix` of side `side` equals its transpose bit for bit.
pub(super) fn bitwise_symmetric(matrix: &[f64], side: usize) -> bool {
    matrix.len() == side * side
        && (0..side).all(|i| {
            (i + 1..side).all(|j| matrix[i * side + j].to_bits() == matrix[j * side + i].to_bits())
        })
}

/// Whether the rows `from..side` of the row-major square `matrix` of side `side` are the identity's bit for bit: 1.0
/// on the row's own column and +0.0 on every other. No row before `from` is read.
pub(super) fn identity_rows(matrix: &[f64], side: usize, from: usize) -> bool {
    matrix.len() == side * side
        && (from..side).all(|i| {
            let row = &matrix[i * side..(i + 1) * side];
            (0..side).all(|j| row[j].to_bits() == if i == j { 1.0_f64 } else { 0.0 }.to_bits())
        })
}

/// The sum and the magnitude of `z^T W z` for the row-major symmetric `matrix` of side `z.len()`, in the order of the
/// module's description.
pub(super) fn upper_row_form(matrix: &[f64], z: &[f64]) -> (f64, f64) {
    upper_row_form_from(&mut None, matrix, z, 0)
}

/// Each row's sum (`.1`) and magnitude (`.2`) over a state's columns right of its diagonal, with the bits of the
/// state's entries (`.0`).
pub(super) struct StateForms(Vec<u64>, Vec<f64>, Vec<f64>);

/// [`upper_row_form`] of `matrix` and `z`, whose first `state` entries are a state and whose others are a held force.
/// A row's columns run from left to right, so its sum over the state's columns comes first: those sums are kept in
/// `slot`, a call whose state has the kept bits in every entry continues them over the force's columns, and any other
/// call computes and keeps its own. Either way a row adds its terms in the order of the module's description, and
/// the rows are added in row order.
pub(super) fn upper_row_form_from(
    slot: &mut Option<StateForms>,
    matrix: &[f64],
    z: &[f64],
    state: usize,
) -> (f64, f64) {
    let side = z.len();
    let held = |kept: &StateForms| {
        let same = |(bits, entry): (&u64, &f64)| *bits == entry.to_bits();
        kept.1.len() == side && kept.0.len() == state && kept.0.iter().zip(z).all(same)
    };
    let kept = slot.take().filter(held).unwrap_or_else(|| {
        let key = z[..state].iter().map(|entry| entry.to_bits()).collect();
        let (mut sums, mut sizes) = (vec![0.0; side], vec![0.0; side]);
        rows_extended(matrix, z, 0..state, &mut sums, &mut sizes);
        StateForms(key, sums, sizes)
    });
    let (mut sums, mut sizes) = (kept.1.clone(), kept.2.clone());
    rows_extended(matrix, z, state..side, &mut sums, &mut sizes);
    *slot = Some(kept);
    let (mut sum, mut magnitude) = (0.0, 0.0);
    for (i, (partial, size)) in sums.iter().zip(&sizes).enumerate() {
        let diagonal = matrix[i * side + i] * z[i];
        sum += z[i] * (diagonal + 2.0 * partial);
        magnitude += z[i].abs() * (diagonal.abs() + 2.0 * size);
    }
    (sum, magnitude)
}

/// Add to each row's sum and magnitude its terms `W_ij z_j` for the columns of `columns` right of its diagonal, in
/// the order of the module's description: the columns inside its block of eight, then the columns every row of a
/// whole block shares. Two calls on adjoining ranges add what one call on their union adds, in its order.
fn rows_extended(
    matrix: &[f64],
    z: &[f64],
    columns: core::ops::Range<usize>,
    sums: &mut [f64],
    sizes: &mut [f64],
) {
    let side = z.len();
    for first in (0..side).step_by(BLOCK) {
        let last = (first + BLOCK).min(side);
        let (mut partial, mut size) = ([0.0_f64; BLOCK], [0.0_f64; BLOCK]);
        partial[..last - first].copy_from_slice(&sums[first..last]);
        size[..last - first].copy_from_slice(&sizes[first..last]);
        for (k, i) in (first..last).enumerate() {
            let from = (i + 1).max(columns.start);
            let inside = from..last.min(columns.end).max(from);
            let row = &matrix[i * side..][inside.clone()];
            for (value, right) in row.iter().zip(&z[inside]) {
                let term = value * right;
                partial[k] += term;
                size[k] += term.abs();
            }
        }
        let from = last.max(columns.start);
        let shared = from..columns.end.max(from);
        if last - first == BLOCK {
            let row = |k: usize| &matrix[(first + k) * side..][shared.clone()];
            let rows: [&[f64]; BLOCK] = core::array::from_fn(row);
            for (j, right) in z[shared.clone()].iter().enumerate() {
                for ((sum, magnitude), row) in partial.iter_mut().zip(&mut size).zip(&rows) {
                    let term = row[j] * right;
                    *sum += term;
                    *magnitude += term.abs();
                }
            }
        }
        sums[first..last].copy_from_slice(&partial[..last - first]);
        sizes[first..last].copy_from_slice(&size[..last - first]);
    }
}

/// The storage `sum_i (v_i^2 + q_i s_i) / 2` of the state `modes = [q_0, v_0, q_1, v_1, ..]` under the row-major
/// `stiffness`, where `s_i` is the sum of `K_ij q_j` over every column in column order. Eight consecutive rows advance
/// together, each with its own sum, and the rows left after the last whole block run one at a time. The rows are then
/// added in row order, each as `0.5 v_i v_i + 0.5 q_i s_i`. Every entry of the square is read once.
pub(super) fn row_storage(stiffness: &[f64], modes: &[f64]) -> f64 {
    let positions: Vec<f64> = modes.iter().step_by(2).copied().collect();
    let side = positions.len().max(1);
    let mut sums = Vec::with_capacity(positions.len());
    let mut blocks = stiffness.chunks_exact(BLOCK * side);
    for block in blocks.by_ref() {
        let rows: [&[f64]; BLOCK] = core::array::from_fn(|k| &block[k * side..(k + 1) * side]);
        let mut partial = [0.0_f64; BLOCK];
        for (j, right) in positions.iter().enumerate() {
            for (sum, row) in partial.iter_mut().zip(&rows) {
                *sum += row[j] * right;
            }
        }
        sums.extend(partial);
    }
    for row in blocks.remainder().chunks_exact(side) {
        let terms = row.iter().zip(&positions);
        sums.push(terms.fold(0.0, |sum, (value, right)| sum + value * right));
    }
    let velocities = modes.iter().skip(1).step_by(2);
    let mut energy = 0.0;
    for ((position, sum), velocity) in positions.iter().zip(&sums).zip(velocities) {
        energy += 0.5 * velocity * velocity + 0.5 * position * sum;
    }
    energy
}

/// Each row's sum over the columns of a state, in column order from 0.0, with the bits of the state's entries.
pub(super) struct StateSums {
    key: Vec<u64>,
    sums: Vec<f64>,
}

/// The product of the whole rows `matrix` (row-major, each of length `z.len()`) by `z`, whose first `state` entries are a
/// state and whose others are a held force. Each entry is its row's sum of `row[j] * z[j]` in column order from 0.0.
/// The sums over the state's columns are kept in `slot`: a call whose state has the kept bits in every entry continues
/// from them, and any other call computes and keeps its own. Either way a row adds the terms of one left-to-right sum
/// in that sum's order, so each entry has its bits.
pub(super) fn product_from(
    slot: &mut Option<StateSums>,
    matrix: &[f64],
    z: &[f64],
    state: usize,
) -> Vec<f64> {
    let (head, tail) = z.split_at(state);
    let rows = matrix.len() / z.len().max(1);
    let held = |kept: &StateSums| {
        let same = |(bits, entry): (&u64, &f64)| *bits == entry.to_bits();
        kept.sums.len() == rows && kept.key.len() == state && kept.key.iter().zip(head).all(same)
    };
    let kept = match slot.take().filter(held) {
        Some(kept) => kept,
        None => StateSums {
            key: head.iter().map(|entry| entry.to_bits()).collect(),
            sums: row_sums(matrix, z.len(), 0, head, &vec![0.0; rows]),
        },
    };
    let product = row_sums(matrix, z.len(), state, tail, &kept.sums);
    *slot = Some(kept);
    product
}

/// For each whole row of `matrix` (rows of length `side`): its entry of `start` plus `row[first + j] * values[j]` in
/// column order. Eight rows advance together, each with its own sum; the rows after the last whole block run alone.
fn row_sums(matrix: &[f64], side: usize, first: usize, values: &[f64], start: &[f64]) -> Vec<f64> {
    let mut sums = start.to_vec();
    let columns = first..first + values.len();
    let mut blocks = matrix.chunks_exact(BLOCK * side);
    let mut done = 0;
    for block in blocks.by_ref() {
        let rows: [&[f64]; BLOCK] = core::array::from_fn(|k| &block[k * side..][columns.clone()]);
        let together = &mut sums[done..done + BLOCK];
        for (j, value) in values.iter().enumerate() {
            for (sum, row) in together.iter_mut().zip(&rows) {
                *sum += row[j] * value;
            }
        }
        done += BLOCK;
    }
    for (sum, row) in sums[done..]
        .iter_mut()
        .zip(blocks.remainder().chunks_exact(side))
    {
        *sum = row[columns.clone()]
            .iter()
            .zip(values)
            .fold(*sum, |s, (a, b)| s + a * b);
    }
    sums
}

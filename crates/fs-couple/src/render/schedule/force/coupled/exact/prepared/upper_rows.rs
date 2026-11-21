//! The quadratic form `z^T W z` of an exactly symmetric `W`, summed by rows over the upper triangle.
//!
//! Row `i` contributes `z_i (W_ii z_i + 2 s_i)`, where `s_i` is the sum of `W_ij z_j` over the columns `j > i`, taken in
//! column order. That is half the terms of the full square. Eight consecutive rows advance together, each with its own
//! sum, so one row's additions do not wait for another row's: a row first takes its columns inside its block of eight
//! and then the columns every row of the block shares. The rows are then added in row order. The magnitude is the same
//! sum over the absolute terms. The order is a function of the side alone and takes no thread count.
//!
//! The form reads `W` from [`UpperBlocks`], a copy made once that holds the upper triangle only, by blocks of eight
//! rows: the eight entries of a column that the block's rows share lie together, so a column is a few vector loads
//! where eight rows of the square would need eight scalar ones. The lower triangle is never copied, so the form is
//! `z^T W z` only for a `W` that equals its transpose: [`bitwise_symmetric`] is the check a caller makes once before it
//! relies on this order.
//!
//! [`row_storage`] is the other sum this module orders by rows: the storage of unit-mass oscillators under a stiffness
//! that need not equal its transpose, so it reads the whole square.
//!
//! [`identity_rows`] is a check of another kind: whether the last rows of a square are the identity's, so that a
//! product by that square need not compute them. [`held_force_defect`] is another: how far two forces differ past
//! their rounding bounds, with the bounds computed only when the forces differ past the tolerance without them.
//!
//! [`product_from`] is a product of whole rows by `[state | force]` that keeps each row's sum over the state's columns
//! ([`StateSums`]): the next product on the same state adds only the force's columns, the same additions in the same
//! order as one left-to-right sum. It reads the rows from [`RowBlocks`], a copy by blocks of eight rows in which a
//! column's eight entries that a block's rows share lie together, as in [`UpperBlocks`], and shares the blocks out to
//! threads. The sums are kept a row, so a product on another thread count reads the same ones. [`state_sums`] gives a
//! caller the kept sums of a state before its product, computing and keeping them first if need be, and [`row_from`]
//! continues one row from them: the same sum as that row's entry of the product.
//!
//! [`upper_row_form_from`] is the form of `[state | force]` that keeps the same kind of sums ([`StateForms`]): a row
//! takes the state's columns first, so the next form on the same state continues every row over the force's columns.
//! [`upper_row_forms_from`] computes several such forms of one `[state | force]` with their rows shared out to threads
//! by blocks of eight. Which thread computes a row changes no addition, so the bits take no thread count.
//!
//! This file names nothing of its crate. `tests/upper_row_forms.rs` compiles it a second time and holds each sum to a
//! plain loop of the same order over the row-major square, bit for bit.

use core::ops::Range;

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

/// How far the forces `a` and `b` differ past each entry's rounding bound: `max(0, max_i(|a_i - b_i| - r_i))`, a caller
/// refusing them when that exceeds `tolerance`. `bounds` gives the `r`, each at least zero, and runs only when the plain
/// defect `max(0, max_i |a_i - b_i|)` exceeds `tolerance`; otherwise the plain defect is returned. For `r >= 0` the
/// rounded `|a_i - b_i| - r_i` is at most `|a_i - b_i|`, so a plain defect within the tolerance leaves the bounded one
/// within it too: the caller refuses the same forces as if `bounds` always ran, except that an error `bounds` would
/// return is met only past the tolerance.
pub(super) fn held_force_defect<E>(
    a: &[f64],
    b: &[f64],
    tolerance: f64,
    bounds: impl FnOnce() -> Result<Vec<f64>, E>,
) -> Result<f64, E> {
    let plain = a
        .iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    if plain <= tolerance {
        return Ok(plain);
    }
    let bounds = bounds()?;
    Ok(a.iter()
        .zip(b)
        .zip(&bounds)
        .fold(0.0_f64, |m, ((x, y), r)| m.max((x - y).abs() - r)))
}

/// The upper triangle of a row-major square, by blocks of eight rows. Each block holds its own square of eight rows
/// and columns, row by row, with +0.0 below its diagonal and past the side; then, for a whole block, each column
/// right of the block with the block's eight entries together.
pub(super) struct UpperBlocks {
    starts: Vec<usize>,
    entries: Vec<f64>,
}

impl UpperBlocks {
    /// The copy of the row-major square `matrix` of side `side`. Its lower triangle is not read.
    pub(super) fn new(matrix: &[f64], side: usize) -> Self {
        let (mut starts, mut entries) = (Vec::new(), Vec::new());
        for first in (0..side).step_by(BLOCK) {
            let last = (first + BLOCK).min(side);
            starts.push(entries.len());
            for i in first..first + BLOCK {
                let upper = |j: usize| (i <= j && j < last).then(|| matrix[i * side + j]);
                entries.extend((first..first + BLOCK).map(|j| upper(j).unwrap_or(0.0)));
            }
            if last - first == BLOCK {
                for j in last..side {
                    entries.extend((first..last).map(|i| matrix[i * side + j]));
                }
            }
        }
        Self { starts, entries }
    }

    /// Row `i`'s entry on the diagonal.
    fn diagonal(&self, i: usize) -> f64 {
        self.entries[self.starts[i / BLOCK] + (i % BLOCK) * (BLOCK + 1)]
    }
}

/// The sum and the magnitude of `z^T W z`, `W` the symmetric square `upper` holds, in the module's order.
#[allow(dead_code, reason = "only tests/upper_row_forms.rs calls this one")]
pub(super) fn upper_row_form(upper: &UpperBlocks, z: &[f64]) -> (f64, f64) {
    upper_row_form_from(&mut None, upper, z, 0)
}

/// The bits of a state's entries (`key`) and four vectors of rows (`rows`): each row's sum and magnitude over the
/// state's columns right of its diagonal, then the same over all its columns right of its diagonal as the last form
/// left them, which are kept only so that a form allocates no vector.
#[derive(Default)]
pub(super) struct StateForms {
    key: Vec<u64>,
    rows: [Vec<f64>; 4],
}

/// [`upper_row_form`] of `upper` and `z`, whose first `state` entries are a state and whose others are a held force.
/// A row's columns run from left to right, so its sum over the state's columns comes first: those sums are kept in
/// `slot`, a call whose state has the kept bits in every entry continues them over the force's columns, and any other
/// call computes and keeps its own. Either way a row adds its terms in the module's order, and rows add in row order.
#[allow(dead_code, reason = "only tests/upper_row_forms.rs calls this one")]
pub(super) fn upper_row_form_from(
    slot: &mut Option<StateForms>,
    upper: &UpperBlocks,
    z: &[f64],
    state: usize,
) -> (f64, f64) {
    upper_row_forms_from(&mut [(upper, slot)], z, state, 1)[0]
}

/// [`upper_row_form_from`] of each `(upper, slot)` of `forms`, with the same `z` and `state`, on `threads` threads (the
/// calling thread among them). Each thread takes a contiguous range of whole blocks of eight rows, cut where the upper
/// triangle's entries before it reach an equal share, and runs every form over its rows. A row adds its terms in the
/// order of the module's description whichever thread runs it, and the calling thread adds the rows in row order, so
/// each form has the bits of [`upper_row_form_from`] for every thread count. One set of threads serves every form.
pub(super) fn upper_row_forms_from(
    forms: &mut [(&UpperBlocks, &mut Option<StateForms>)],
    z: &[f64],
    state: usize,
    threads: usize,
) -> Vec<(f64, f64)> {
    let side = z.len();
    let mut kept: Vec<_> = forms
        .iter_mut()
        .map(|(_, slot)| slot.take().unwrap_or_default())
        .collect();
    let held: Vec<_> = kept.iter_mut().map(|kept| kept.ready(z, state)).collect();
    let (starts, total) = forms.first().map_or((&[][..], 0), |(upper, _)| {
        (&upper.starts[..], upper.entries.len())
    });
    let threads = threads.clamp(1, starts.len().max(1));
    let cuts: Vec<_> = (0..=threads)
        .map(|t| starts.partition_point(|s| s * threads < total * t))
        .collect();
    let mut work: Vec<Vec<_>> = (0..threads).map(|_| Vec::new()).collect();
    for (((upper, _), kept), held) in forms.iter().zip(&mut kept).zip(&held) {
        let mut rest = kept.rows.each_mut().map(|row| &mut row[..]);
        for (part, cut) in work.iter_mut().zip(cuts.windows(2)) {
            let rows = (cut[1] * BLOCK).min(side) - (cut[0] * BLOCK).min(side);
            let mine = rest.each_mut().map(|row| {
                let (mine, others) = core::mem::take(row).split_at_mut(rows);
                *row = others;
                mine
            });
            part.push((*upper, *held, mine));
        }
    }
    let run = |part: Vec<(&UpperBlocks, bool, [&mut [f64]; 4])>, blocks: Range<usize>| {
        for (upper, held, [sums, sizes, rows, row_sizes]) in part {
            if !held {
                rows_extended(upper, z, 0..state, blocks.clone(), sums, sizes);
            }
            rows.copy_from_slice(sums);
            row_sizes.copy_from_slice(sizes);
            rows_extended(upper, z, state..side, blocks.clone(), rows, row_sizes);
        }
    };
    // The scope joins every thread it started before it returns, and panics if one of them did.
    let mut parts = work
        .into_iter()
        .zip(cuts.windows(2).map(|cut| cut[0]..cut[1]));
    let (mine, blocks) = parts.next().unwrap_or_default();
    if threads > 1 {
        std::thread::scope(|scope| {
            for (part, blocks) in parts {
                scope.spawn(move || run(part, blocks));
            }
            run(mine, blocks);
        });
    } else {
        drop(parts);
        run(mine, blocks);
    }
    let mut values = Vec::with_capacity(forms.len());
    for ((upper, slot), kept) in forms.iter_mut().zip(kept) {
        let (mut sum, mut magnitude) = (0.0, 0.0);
        for (i, (partial, size)) in kept.rows[2].iter().zip(&kept.rows[3]).enumerate() {
            let diagonal = upper.diagonal(i) * z[i];
            sum += z[i] * (diagonal + 2.0 * partial);
            magnitude += z[i].abs() * (diagonal.abs() + 2.0 * size);
        }
        **slot = Some(kept);
        values.push((sum, magnitude));
    }
    values
}

impl StateForms {
    /// Make the four vectors of rows `z.len()` long, and say whether the kept sums are over `state` entries with the
    /// bits of `z`'s; when they are not, they become `z`'s, set to +0.0 and still to be computed.
    fn ready(&mut self, z: &[f64], state: usize) -> bool {
        let same = |(bits, entry): (&u64, &f64)| *bits == entry.to_bits();
        let held = self.rows[0].len() == z.len()
            && self.key.len() == state
            && self.key.iter().zip(z).all(same);
        for row in &mut self.rows {
            row.resize(z.len(), 0.0);
        }
        if !held {
            self.key.clear();
            self.key.extend(z[..state].iter().map(|v| v.to_bits()));
            self.rows[0].fill(0.0);
            self.rows[1].fill(0.0);
        }
        held
    }
}

/// Add to the sum and magnitude of each row of the blocks `blocks` its terms `W_ij z_j` for the columns of `columns`
/// right of its diagonal, in the order of the module's description: the columns inside its block of eight, then the
/// columns every row of a whole block shares. `sums` and `sizes` hold those rows only, from the first block's first
/// row. Two calls on adjoining ranges of columns add what one call on their union adds, in its order.
fn rows_extended(
    upper: &UpperBlocks,
    z: &[f64],
    columns: Range<usize>,
    blocks: Range<usize>,
    sums: &mut [f64],
    sizes: &mut [f64],
) {
    let side = z.len();
    let offset = blocks.start * BLOCK;
    for (start, first) in upper.starts[blocks]
        .iter()
        .zip((offset..side).step_by(BLOCK))
    {
        let last = (first + BLOCK).min(side);
        let (square, shared) = upper.entries[*start..].split_at(BLOCK * BLOCK);
        let (mut partial, mut size) = ([0.0_f64; BLOCK], [0.0_f64; BLOCK]);
        partial[..last - first].copy_from_slice(&sums[first - offset..last - offset]);
        size[..last - first].copy_from_slice(&sizes[first - offset..last - offset]);
        for (k, i) in (first..last).enumerate() {
            let inside = (i + 1).max(columns.start)..last.min(columns.end);
            if inside.is_empty() {
                continue;
            }
            let row = &square[k * BLOCK..][inside.start - first..inside.end - first];
            for (value, right) in row.iter().zip(&z[inside]) {
                let term = value * right;
                partial[k] += term;
                size[k] += term.abs();
            }
        }
        let across = last.max(columns.start)..columns.end;
        if last - first == BLOCK && !across.is_empty() {
            let cells = &shared[(across.start - last) * BLOCK..(across.end - last) * BLOCK];
            for (cell, right) in cells.as_chunks::<BLOCK>().0.iter().zip(&z[across]) {
                for ((sum, magnitude), value) in partial.iter_mut().zip(&mut size).zip(cell) {
                    let term = value * right;
                    *sum += term;
                    *magnitude += term.abs();
                }
            }
        }
        sums[first - offset..last - offset].copy_from_slice(&partial[..last - first]);
        sizes[first - offset..last - offset].copy_from_slice(&size[..last - first]);
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

/// The whole rows of a row-major matrix by blocks of eight: each whole block's columns in order, with the block's
/// eight entries of a column together; then the rows after the last whole block as they are.
pub(super) struct RowBlocks {
    side: usize,
    rows: usize,
    entries: Vec<f64>,
}

impl RowBlocks {
    /// The copy of the whole rows `matrix`, each of `side` entries.
    pub(super) fn new(matrix: &[f64], side: usize) -> Self {
        let rows = matrix.len() / side.max(1);
        let whole = rows / BLOCK * BLOCK;
        let mut entries = Vec::with_capacity(rows * side);
        for first in (0..whole).step_by(BLOCK) {
            for j in 0..side {
                entries.extend((first..first + BLOCK).map(|i| matrix[i * side + j]));
            }
        }
        entries.extend_from_slice(&matrix[whole * side..rows * side]);
        Self {
            side,
            rows,
            entries,
        }
    }
}

/// Each row's sum over the columns of a state, in column order from 0.0, with the bits of the state's entries.
#[derive(Default)]
pub(super) struct StateSums {
    key: Vec<u64>,
    sums: Vec<f64>,
}

impl StateSums {
    /// Whether these are the sums of `rows` rows over a state with the bits of `head` in every entry. If not, they
    /// become that state's: every sum +0.0 and still to be computed.
    fn hold(&mut self, rows: usize, head: &[f64]) -> bool {
        let same = |(bits, entry): (&u64, &f64)| *bits == entry.to_bits();
        let held = self.sums.len() == rows
            && self.key.len() == head.len()
            && self.key.iter().zip(head).all(same);
        if !held {
            self.key.clear();
            self.key.extend(head.iter().map(|v| v.to_bits()));
            self.sums.clear();
            self.sums.resize(rows, 0.0);
        }
        held
    }
}

/// The product of the whole rows that `blocks` holds (each of length `z.len()`) by `z`, whose first `state` entries are
/// a state and whose others are a held force, on `threads` threads (the calling thread among them). Each entry is its
/// row's sum of `row[j] * z[j]` in column order from 0.0. The sums over the state's columns are kept in `slot`, one a
/// row: a call whose state has the kept bits in every entry continues from them, and any other call computes and keeps
/// every row's. Each thread takes a contiguous range of whole blocks of eight rows (the last range also the rows after
/// the last whole block) and reads and writes only its own rows' sums, so the split may change from one call to the
/// next. Either way a row adds the terms of one left-to-right sum in that sum's order, so each entry has its bits for
/// every thread count.
pub(super) fn product_from(
    slot: &mut Option<StateSums>,
    blocks: &RowBlocks,
    z: &[f64],
    state: usize,
    threads: usize,
) -> Vec<f64> {
    let (head, tail) = z.split_at(state);
    let mut kept = slot.take().unwrap_or_default();
    let held = kept.hold(blocks.rows, head);
    let mut product = vec![0.0; blocks.rows];
    let run = |start: usize, sums: &mut [f64], product: &mut [f64]| {
        if !held {
            row_sums(blocks, start, 0, head, sums);
        }
        product.copy_from_slice(sums);
        row_sums(blocks, start, state, tail, product);
    };
    let per = (blocks.rows.div_ceil(BLOCK).div_ceil(threads.max(1)) * BLOCK).max(1);
    let rows = kept.sums.chunks_mut(per).zip(product.chunks_mut(per));
    let mut parts = (0..).step_by(per).zip(rows);
    let first = parts.next();
    // The scope joins every thread it started before it returns, and panics if one of them did. One thread starts none.
    std::thread::scope(|scope| {
        for (start, (sums, product)) in parts {
            scope.spawn(move || run(start, sums, product));
        }
        if let Some((start, (sums, product))) = first {
            run(start, sums, product);
        }
    });
    *slot = Some(kept);
    product
}

/// Every row's sum over the columns of the state `head`, in column order from 0.0: the kept ones when `slot` holds
/// that state's, else computed on `threads` threads, which share the rows as [`product_from`] does, and kept, so a
/// product on the same state continues from them. Each sum has the same bits either way and for every thread count.
pub(super) fn state_sums<'a>(
    slot: &'a mut Option<StateSums>,
    blocks: &RowBlocks,
    head: &[f64],
    threads: usize,
) -> &'a [f64] {
    let mut kept = slot.take().unwrap_or_default();
    if !kept.hold(blocks.rows, head) {
        let per = (blocks.rows.div_ceil(BLOCK).div_ceil(threads.max(1)) * BLOCK).max(1);
        let mut parts = (0..).step_by(per).zip(kept.sums.chunks_mut(per));
        let first = parts.next();
        // The scope joins every thread it started before it returns, and panics if one of them did. One thread starts none.
        std::thread::scope(|scope| {
            for (start, sums) in parts {
                scope.spawn(move || row_sums(blocks, start, 0, head, sums));
            }
            if let Some((start, sums)) = first {
                row_sums(blocks, start, 0, head, sums);
            }
        });
    }
    &slot.insert(kept).sums
}

/// Row `index`'s sum of `row[j] * z[j]` over every column in order from 0.0, where `row` is that row of the matrix
/// whose sums over the first `state` columns `sums` holds ([`state_sums`]): its kept sum continued over the others.
pub(super) fn row_from(sums: &[f64], index: usize, row: &[f64], z: &[f64], state: usize) -> f64 {
    let terms = row[state..].iter().zip(&z[state..]);
    terms.fold(sums[index], |s, (a, b)| s + a * b)
}

/// Add to each entry of `sums`, the rows from `start` (a multiple of eight) on, its row's `row[first + j] * values[j]` in
/// column order. Eight rows advance together, each with its own sum; the rows after the last whole block run alone.
fn row_sums(blocks: &RowBlocks, start: usize, first: usize, values: &[f64], sums: &mut [f64]) {
    let side = blocks.side;
    let columns = first * BLOCK..(first + values.len()) * BLOCK;
    let (together, alone) = sums.as_chunks_mut::<BLOCK>();
    let whole = blocks.entries[start * side..].chunks_exact(BLOCK * side.max(1));
    for (block, sums) in whole.zip(together) {
        let mut partial = *sums;
        let cells = block[columns.clone()].as_chunks::<BLOCK>().0;
        for (cell, value) in cells.iter().zip(values) {
            for (sum, entry) in partial.iter_mut().zip(cell) {
                *sum += entry * value;
            }
        }
        *sums = partial;
    }
    let rest = &blocks.entries[blocks.rows / BLOCK * BLOCK * side..];
    for (sum, row) in alone.iter_mut().zip(rest.chunks_exact(side.max(1))) {
        *sum = row[first..first + values.len()]
            .iter()
            .zip(values)
            .fold(*sum, |s, (a, b)| s + a * b);
    }
}

//! Blocked, thread-parallel complex elimination with the ORIGINAL operation
//! graph. The unblocked owner updates element (i, c) once per earlier pivot k,
//! in ascending k, as `difference(a, product(f_ik, u_kc))`. This module only
//! reorders WHICH elements are updated when; every element still receives that
//! exact sequence. Factors, permutations and singular indices are therefore
//! bit-identical to [`super::lu_complex_with_cancel`] for every block size and
//! thread count.
//!
//! Per block of `block` pivot columns [k0, k1):
//! 1. Panel: the original magnitude pivots, lowest-index ties and full-row
//!    swaps, eliminating only inside the block's columns.
//! 2. Upper block U12 (rows k0..k1, columns k1..n): each row receives the
//!    block's earlier pivot rows in ascending k. Threads own column chunks.
//! 3. Trailing update (rows k1..n, columns k1..n): each element receives the
//!    block's pivots in ascending k. Threads own row chunks; columns are tiled
//!    so the pivot rows' tile stays cache-resident.
//!
//! A swap at pivot k exchanges rows whose columns past the block carry the
//! same earlier updates in both rows; the delayed updates use each row's own
//! (swapped) multipliers. Relabelling rows commutes with the delayed updates.
use super::{C64, LuComplexError, difference, fma, poll, product};

/// Pivot columns per block. The trailing update reads the block's pivot rows
/// once per row tile instead of once per pivot.
pub(super) const BLOCK: usize = 64;
/// Trailing-update column tile: BLOCK rows by TILE columns of pivot data is
/// 256 KiB, inside a per-core L2.
const TILE: usize = 256;

/// Panel elimination of pivots [k0, k1): `factor_body` with its column range
/// cut to the block. Polls at every pivot.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
pub(super) fn panel_body(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    k0: usize,
    k1: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    for k in k0..k1 {
        poll(cancelled)?;
        let (mut best, mut best_i) = (m[k * n + k].abs(), k);
        for i in k + 1..n {
            let v = m[i * n + k].abs();
            if v > best {
                best = v;
                best_i = i;
            }
        }
        if best == 0.0 {
            return Err(LuComplexError::Singular { index: k });
        }
        if best_i != k {
            for c in 0..n {
                m.swap(k * n + c, best_i * n + c);
            }
            perm.swap(k, best_i);
        }
        let piv = m[k * n + k];
        for i in k + 1..n {
            let f = m[i * n + k] / piv;
            m[i * n + k] = f;
            let (upper, lower) = m.split_at_mut(i * n);
            let pivot_row = &upper[k * n + k + 1..k * n + k1];
            for (target, &source) in lower[k + 1..k1].iter_mut().zip(pivot_row) {
                *target = difference(*target, product(f, source));
            }
        }
    }
    Ok(())
}

/// One column chunk of U12. `rows[r]` is (row k0+r's multipliers in columns
/// k0..k1, that row's chunk). Row r receives rows 0..r in ascending order.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
pub(super) fn upper_body(rows: &mut [(&[C64], &mut [C64])]) {
    for r in 1..rows.len() {
        let (done, rest) = rows.split_at_mut(r);
        let (multipliers, target) = &mut rest[0];
        for (j, (_, pivot)) in done.iter().enumerate() {
            let f = multipliers[j];
            for (t, &s) in target.iter_mut().zip(pivot.iter()) {
                *t = difference(*t, product(f, s));
            }
        }
    }
}

/// One row chunk of the trailing update. `rows` holds whole rows from k1 on;
/// `above` holds rows 0..k1, so pivot row k starts at `k * n`.
#[inline(always)]
#[allow(clippy::inline_always)] // must inline INTO the fma capsule's target_feature fns
pub(super) fn trailing_body(rows: &mut [C64], above: &[C64], n: usize, k0: usize, k1: usize) {
    let mut start = k1;
    while start < n {
        let end = (start + TILE).min(n);
        for row in rows.chunks_exact_mut(n) {
            let (left, right) = row.split_at_mut(k1);
            let target = &mut right[start - k1..end - k1];
            for k in k0..k1 {
                let f = left[k];
                let source = &above[k * n + start..k * n + end];
                for (t, &s) in target.iter_mut().zip(source) {
                    *t = difference(*t, product(f, s));
                }
            }
        }
        start = end;
    }
}

/// Contiguous `parts`-way split of `0..len`; part p is `bounds[p]..bounds[p+1]`.
fn bounds(len: usize, parts: usize) -> Vec<usize> {
    (0..=parts).map(|p| p * len / parts).collect()
}

/// Blocked elimination. `threads` workers (at least one) share each U12 and
/// trailing phase; the panel runs on the calling thread. Polls at every pivot
/// and before each U12 and trailing phase, so between polls at most one
/// block's U12 or trailing update runs.
pub(super) fn factor_blocked(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    block: usize,
    threads: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    let block = block.max(1);
    let threads = threads.max(1);
    let mut k0 = 0;
    while k0 < n {
        let k1 = (k0 + block).min(n);
        fma::panel(m, perm, n, k0, k1, cancelled)?;
        if k1 < n {
            poll(cancelled)?;
            let width = n - k1;
            let parts = threads.min(width);
            let edges = bounds(width, parts);
            let (top, bottom) = m.split_at_mut(k1 * n);
            if k1 - k0 > 1 {
                let mut chunks: Vec<Vec<(&[C64], &mut [C64])>> =
                    (0..parts).map(|_| Vec::with_capacity(k1 - k0)).collect();
                for row in top[k0 * n..].chunks_exact_mut(n) {
                    let (left, mut right) = row.split_at_mut(k1);
                    let left: &[C64] = left;
                    let multipliers = &left[k0..k1];
                    for (p, chunk) in chunks.iter_mut().enumerate() {
                        let (head, tail) =
                            core::mem::take(&mut right).split_at_mut(edges[p + 1] - edges[p]);
                        chunk.push((multipliers, head));
                        right = tail;
                    }
                }
                if parts == 1 {
                    fma::upper(&mut chunks[0]);
                } else {
                    std::thread::scope(|scope| {
                        for chunk in &mut chunks {
                            scope.spawn(move || fma::upper(chunk));
                        }
                    });
                }
            }
            poll(cancelled)?;
            let above: &[C64] = top;
            if parts == 1 {
                fma::trailing(bottom, above, n, k0, k1);
            } else {
                let mut rest = bottom;
                std::thread::scope(|scope| {
                    for p in 0..parts {
                        let (head, tail) =
                            core::mem::take(&mut rest).split_at_mut((edges[p + 1] - edges[p]) * n);
                        rest = tail;
                        scope.spawn(move || fma::trailing(head, above, n, k0, k1));
                    }
                });
            }
        }
        k0 = k1;
    }
    Ok(())
}

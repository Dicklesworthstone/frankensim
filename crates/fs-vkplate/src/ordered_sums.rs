//! Row sums that keep `Iterator::sum`'s order, eight rows at a time.
//!
//! `sums_of_products` gives each `out[r]` the bits of
//! `row(r).iter().zip(x).map(|(a, b)| a * b).sum::<f64>()`: the products are
//! added one at a time in ascending index, starting from the value an empty
//! `f64` sum returns (read from the toolchain, not written here). One such
//! chain waits on each addition before the next; eight chains run together so
//! that their additions overlap. No chain's order changes, so no bit does.
//! Rows of unequal length within a block of eight are summed one at a time.
//! `sums_of_products_from` is the same with each chain starting from a given value.
//! The file names nothing of its crate, so a test can include it by `#[path]`.

/// The value `Iterator::sum` returns for an empty `f64` iterator, which each of its sums starts from.
#[inline]
pub(crate) fn empty_sum() -> f64 {
    std::iter::empty::<f64>().sum()
}

/// One row's sum from `init`, its products added in ascending index.
#[inline]
fn one_row(init: f64, row: &[f64], x: &[f64]) -> f64 {
    let mut acc = init;
    for (a, b) in row.iter().zip(x) {
        acc += a * b;
    }
    acc
}

/// `out[r] = row(r) . x` for `r < out.len()`, each with the bits of
/// `row(r).iter().zip(x).map(|(a, b)| a * b).sum::<f64>()`.
pub(crate) fn sums_of_products<'a>(row: impl Fn(usize) -> &'a [f64], x: &[f64], out: &mut [f64]) {
    sums_of_products_from(empty_sum(), row, x, out);
}

/// `out[r]` for `r < out.len()` is `init` with each `row(r)[i] * x[i]` added in ascending `i`, the bits of
/// `row(r).iter().zip(x).fold(init, |acc, (a, b)| acc + a * b)`.
pub(crate) fn sums_of_products_from<'a>(init: f64, row: impl Fn(usize) -> &'a [f64], x: &[f64], out: &mut [f64]) {
    let blocks = out.len() / 8;
    for block in 0..blocks {
        let first = 8 * block;
        let rows: [&[f64]; 8] = std::array::from_fn(|j| row(first + j));
        let n = rows[0].len().min(x.len());
        if rows.iter().any(|r| r.len().min(x.len()) != n) {
            for (j, r) in rows.iter().enumerate() {
                out[first + j] = one_row(init, r, x);
            }
            continue;
        }
        let rows: [&[f64]; 8] = rows.map(|r| &r[..n]);
        let mut acc = [init; 8];
        for (i, xi) in x[..n].iter().enumerate() {
            for (a, r) in acc.iter_mut().zip(&rows) {
                *a += r[i] * xi;
            }
        }
        out[first..first + 8].copy_from_slice(&acc);
    }
    for (r, o) in out.iter_mut().enumerate().skip(8 * blocks) {
        *o = one_row(init, row(r), x);
    }
}

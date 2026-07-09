//! fs-sparse perf lane (bead wsbf, segment 1): compact-index CSR,
//! sharded parallel SpMV, and tiled parallel assembly — every kernel
//! BITWISE-equal to its serial/wide twin by construction, so the perf
//! program never trades away the determinism contract.
//!
//! - [`CsrCompact`]: u32 column indices halve the per-nnz index
//!   traffic that bounds CSR SpMV (the values stream is irreducible;
//!   the index stream is not). `row_ptr` stays `usize` — nnz may
//!   exceed u32 at production scale. Accumulation order is IDENTICAL
//!   to [`Csr::spmv`] (ascending-column `mul_add` from +0.0), so
//!   equality is bitwise, gated in conformance.
//! - [`CsrCompact::spmv_sharded`]: disjoint contiguous row ranges
//!   balanced by nnz prefix, scoped threads, each thread owning its
//!   `y` slice via `split_at_mut` (write-side first-touch for NUMA).
//!   Per-row accumulation is untouched, so the result is bitwise equal
//!   to the serial kernel at EVERY thread count — the bead's item-2
//!   determinism constraint holds by construction, no golden bump.
//! - [`Coo::assemble_parallel`]: row-range tiles bucketed in one
//!   serial pass (preserving GLOBAL insertion order within each tile),
//!   per-tile stable sort + accumulation on scoped threads, tiles
//!   concatenated in order. Rows never span tiles, so every (row, col)
//!   duplicate chain accumulates in exactly the serial order — bitwise
//!   equal to [`Coo::assemble`] for ANY thread count (the G5
//!   criterion, gated).

use crate::{Coo, Csr};

/// CSR with compact u32 column indices (the SpMV bandwidth diet).
#[derive(Debug, Clone, PartialEq)]
pub struct CsrCompact {
    nrows: usize,
    ncols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
    vals: Vec<f64>,
}

impl CsrCompact {
    /// Convert from canonical CSR. Panics (structured) when `ncols`
    /// exceeds the u32 index space — the caller keeps the wide format
    /// there; silently truncating an index would be corruption.
    #[must_use]
    pub fn from_csr(a: &Csr) -> CsrCompact {
        assert!(
            a.ncols() <= u32::MAX as usize,
            "compact CSR needs ncols <= u32::MAX, got {}",
            a.ncols()
        );
        let mut row_ptr = Vec::with_capacity(a.nrows() + 1);
        let mut col_idx = Vec::with_capacity(a.nnz());
        let mut vals = Vec::with_capacity(a.nnz());
        row_ptr.push(0usize);
        for r in 0..a.nrows() {
            let (cols, v) = a.row(r);
            for &c in cols {
                col_idx.push(u32::try_from(c).expect("checked ncols bound"));
            }
            vals.extend_from_slice(v);
            row_ptr.push(col_idx.len());
        }
        CsrCompact {
            nrows: a.nrows(),
            ncols: a.ncols(),
            row_ptr,
            col_idx,
            vals,
        }
    }

    /// Row count.
    #[must_use]
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Column count.
    #[must_use]
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Stored nonzeros.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.vals.len()
    }

    /// Serial SpMV — the same ascending-column `mul_add` accumulation
    /// from +0.0 as [`Csr::spmv`]; bitwise equality is gated.
    pub fn spmv(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        for (r, out) in y.iter_mut().enumerate() {
            let lo = self.row_ptr[r];
            let hi = self.row_ptr[r + 1];
            let mut acc = 0.0f64;
            // Slice windows: checked-free iteration (indexed loops were
            // MEASURED slower than the wide kernel from bounds checks).
            for (&c, &v) in self.col_idx[lo..hi].iter().zip(&self.vals[lo..hi]) {
                acc = v.mul_add(x[c as usize], acc);
            }
            *out = acc;
        }
    }

    /// nnz-balanced contiguous row shards for `threads` workers.
    fn shard_bounds(&self, threads: usize) -> Vec<usize> {
        let t = threads.max(1);
        let total = self.nnz();
        let mut bounds = Vec::with_capacity(t + 1);
        bounds.push(0usize);
        let mut next_target = 1usize;
        for r in 0..self.nrows {
            let filled = self.row_ptr[r + 1];
            while next_target < t && filled * t >= next_target * total {
                bounds.push(r + 1);
                next_target += 1;
            }
        }
        while bounds.len() < t {
            bounds.push(self.nrows);
        }
        bounds.push(self.nrows);
        bounds
    }

    /// Sharded parallel SpMV: bitwise equal to [`Self::spmv`] at every
    /// thread count (disjoint row ranges; per-row order untouched).
    pub fn spmv_sharded(&self, x: &[f64], y: &mut [f64], threads: usize) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        let bounds = self.shard_bounds(threads);
        std::thread::scope(|scope| {
            let mut rest = y;
            let mut offset = 0usize;
            for w in bounds.windows(2) {
                let (lo, hi) = (w[0], w[1]);
                let (mine, tail) = rest.split_at_mut(hi - offset);
                rest = tail;
                offset = hi;
                if lo == hi {
                    continue;
                }
                scope.spawn(move || {
                    for r in lo..hi {
                        let a = self.row_ptr[r];
                        let b = self.row_ptr[r + 1];
                        let mut acc = 0.0f64;
                        for (&c, &v) in self.col_idx[a..b].iter().zip(&self.vals[a..b]) {
                            acc = v.mul_add(x[c as usize], acc);
                        }
                        mine[r - lo] = acc;
                    }
                });
            }
        });
    }
}

impl Coo {
    /// Tiled PARALLEL assembly, bitwise equal to [`Coo::assemble`] for
    /// any `threads` (rows never span tiles; each duplicate chain
    /// accumulates in the global insertion order the serial path uses).
    #[must_use]
    pub fn assemble_parallel(&self, threads: usize) -> Csr {
        let t = threads.max(1);
        let nrows = self.nrows;
        // Row-range tiles balanced by STAGED triplet count.
        let mut per_row = vec![0usize; nrows + 1];
        for &r in &self.rows {
            per_row[r + 1] += 1;
        }
        for r in 0..nrows {
            per_row[r + 1] += per_row[r];
        }
        let total = self.vals.len();
        let mut tile_of_row = vec![0usize; nrows];
        {
            let mut tile = 0usize;
            for (r, slot) in tile_of_row.iter_mut().enumerate() {
                while tile + 1 < t && per_row[r + 1] * t > (tile + 1) * total {
                    tile += 1;
                }
                *slot = tile;
            }
        }
        // One serial bucketing pass: global order preserved per tile.
        let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); t];
        for (i, &r) in self.rows.iter().enumerate() {
            buckets[tile_of_row[r]].push(i);
        }
        // Per-tile canonical accumulation on scoped threads.
        struct TileOut {
            row_counts: Vec<usize>,
            col_idx: Vec<usize>,
            vals: Vec<f64>,
            row_lo: usize,
        }
        let mut outs: Vec<Option<TileOut>> = Vec::new();
        for _ in 0..t {
            outs.push(None);
        }
        // Per-tile contiguous row ranges (tile_of_row is monotone).
        let ranges: Vec<(usize, usize)> = (0..t)
            .map(|tile| {
                let lo = tile_of_row.iter().position(|&x| x == tile);
                match lo {
                    Some(lo) => {
                        let hi = tile_of_row
                            .iter()
                            .rposition(|&x| x == tile)
                            .map_or(lo, |p| p + 1);
                        (lo, hi)
                    }
                    None => (0, 0),
                }
            })
            .collect();
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (tile, idxs) in buckets.iter().enumerate() {
                let (row_lo, row_hi) = ranges[tile];
                handles.push((tile, scope.spawn(move || {
                    let mut order: Vec<usize> = idxs.clone();
                    order.sort_by_key(|&i| (self.rows[i], self.cols[i]));
                    let mut row_counts = vec![0usize; row_hi.saturating_sub(row_lo)];
                    let mut col_idx = Vec::new();
                    let mut vals = Vec::new();
                    let mut i = 0usize;
                    while i < order.len() {
                        let (r, c) = (self.rows[order[i]], self.cols[order[i]]);
                        let mut acc = self.vals[order[i]];
                        i += 1;
                        while i < order.len()
                            && self.rows[order[i]] == r
                            && self.cols[order[i]] == c
                        {
                            acc += self.vals[order[i]];
                            i += 1;
                        }
                        row_counts[r - row_lo] += 1;
                        col_idx.push(c);
                        vals.push(acc);
                    }
                    TileOut {
                        row_counts,
                        col_idx,
                        vals,
                        row_lo,
                    }
                })));
            }
            for (tile, h) in handles {
                outs[tile] = Some(h.join().expect("tile worker"));
            }
        });
        // Concatenate tiles in order (rows are contiguous and disjoint).
        let mut row_ptr = vec![0usize; nrows + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        for out in outs.into_iter().flatten() {
            for (dr, &cnt) in out.row_counts.iter().enumerate() {
                row_ptr[out.row_lo + dr + 1] = cnt;
            }
            col_idx.extend_from_slice(&out.col_idx);
            vals.extend_from_slice(&out.vals);
        }
        for r in 0..nrows {
            row_ptr[r + 1] += row_ptr[r];
        }
        Csr::from_parts(nrows, self.ncols, row_ptr, col_idx, vals)
    }
}

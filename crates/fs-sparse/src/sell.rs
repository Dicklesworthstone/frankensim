//! SELL-C-σ — sliced ELL with row sorting (plan §6.2): rows are grouped into
//! chunks of C, and within a sorting window of σ rows, sorted by descending
//! length so chunk widths hug the real work. The layout is column-of-chunk
//! major (lane-fastest), which is exactly what a C-wide SIMD kernel wants on
//! both NEON and AVX-512.
//!
//! v1 determinism decision: per-row lengths are STORED and the scalar kernel
//! iterates each lane to its true length — pad slots exist physically (for
//! the future SIMD lane) but are never read, so SpMV is bit-identical to CSR
//! by construction instead of "identical modulo signed-zero pad accidents".
//! The SIMD kernel (perf follow-up bead) must preserve this equality or
//! justify a golden-hash bump.

use crate::Csr;

/// Structured refusal from [`Sell::from_csr`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SellError {
    /// A SELL chunk must contain at least one lane.
    InvalidChunkHeight {
        /// The refused chunk height.
        chunk_height: usize,
    },
    /// The requested layout cannot be represented by `usize` storage offsets.
    GeometryOverflow {
        /// Layout operation whose result was not representable.
        operation: &'static str,
        /// Left operand of the checked operation.
        lhs: usize,
        /// Right operand of the checked operation.
        rhs: usize,
    },
}

impl std::fmt::Display for SellError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SellError::InvalidChunkHeight { chunk_height } => {
                write!(
                    f,
                    "SELL chunk height must be at least 1, got {chunk_height}"
                )
            }
            SellError::GeometryOverflow {
                operation,
                lhs,
                rhs,
            } => write!(
                f,
                "SELL geometry overflow while {operation}: operands {lhs} and {rhs}"
            ),
        }
    }
}

impl std::error::Error for SellError {}

fn checked_geometry_add(
    operation: &'static str,
    lhs: usize,
    rhs: usize,
) -> Result<usize, SellError> {
    lhs.checked_add(rhs).ok_or(SellError::GeometryOverflow {
        operation,
        lhs,
        rhs,
    })
}

fn checked_geometry_mul(
    operation: &'static str,
    lhs: usize,
    rhs: usize,
) -> Result<usize, SellError> {
    lhs.checked_mul(rhs).ok_or(SellError::GeometryOverflow {
        operation,
        lhs,
        rhs,
    })
}

fn bounded_worker_count(requested: usize, useful_chunks: usize) -> usize {
    let host_parallelism = std::thread::available_parallelism().map_or(1, std::num::NonZero::get);
    requested
        .max(1)
        .min(useful_chunks.max(1))
        .min(host_parallelism)
}

/// Sliced ELL with sorting window. Row `perm[r]`'s data lives in lane
/// `r % C` of chunk `r / C` (r is the SORTED position; `perm` maps sorted
/// position → original row).
#[derive(Debug, Clone, PartialEq)]
pub struct Sell {
    c: usize,
    sigma: usize,
    nrows: usize,
    ncols: usize,
    /// Sorted position → original row index.
    perm: Vec<usize>,
    /// Per sorted-position true row length.
    row_len: Vec<usize>,
    /// Chunk start offsets into `col_idx`/`vals` (len = nchunks + 1).
    chunk_ptr: Vec<usize>,
    /// Lane-fastest storage: entry k of lane l in chunk ch lives at
    /// `chunk_ptr[ch] + k * C + l`. Pad slots hold (0, 0.0) and are inert.
    col_idx: Vec<usize>,
    vals: Vec<f64>,
}

impl Sell {
    /// Convert from canonical CSR. `c` is the chunk height (≥ 1); `sigma`
    /// the sorting window in rows, rounded up to a multiple of `c`
    /// internally. Sorting is STABLE by descending length, so equal-length
    /// rows keep matrix order — the permutation is a pure function of the
    /// structure (deterministic).
    ///
    /// # Errors
    ///
    /// Returns [`SellError::InvalidChunkHeight`] when `c == 0`, or
    /// [`SellError::GeometryOverflow`] when the rounded sorting window or
    /// lane-fastest storage geometry cannot be represented by `usize`.
    pub fn from_csr(a: &Csr, c: usize, sigma: usize) -> Result<Sell, SellError> {
        if c == 0 {
            return Err(SellError::InvalidChunkHeight { chunk_height: c });
        }
        let sigma_chunks = sigma.max(c).div_ceil(c);
        let sigma = checked_geometry_mul("rounding sigma to chunk height", sigma_chunks, c)?;
        let nrows = a.nrows();
        let mut perm: Vec<usize> = (0..nrows).collect();
        for window in perm.chunks_mut(sigma) {
            window.sort_by_key(|&r| std::cmp::Reverse(a.row(r).0.len()));
        }
        let nchunks = nrows.div_ceil(c);
        let row_len: Vec<usize> = perm.iter().map(|&r| a.row(r).0.len()).collect();
        let chunk_ptr_len = checked_geometry_add("sizing the chunk pointer table", nchunks, 1)?;
        let mut chunk_ptr = vec![0usize; chunk_ptr_len];
        for ch in 0..nchunks {
            let row0 = checked_geometry_mul("locating a chunk row start", ch, c)?;
            let next_chunk = checked_geometry_add("locating the next chunk", ch, 1)?;
            let row1 = checked_geometry_mul("locating a chunk row end", next_chunk, c)?.min(nrows);
            let width = (row0..row1).map(|p| row_len[p]).max().unwrap_or(0);
            let chunk_slots = checked_geometry_mul("sizing a chunk", width, c)?;
            chunk_ptr[next_chunk] =
                checked_geometry_add("accumulating chunk storage", chunk_ptr[ch], chunk_slots)?;
        }
        let total = chunk_ptr[nchunks];
        let mut col_idx = vec![0usize; total];
        let mut vals = vec![0.0f64; total];
        for (pos, &orig) in perm.iter().enumerate() {
            let (ch, lane) = (pos / c, pos % c);
            let (cols, values) = a.row(orig);
            for (k, (&cc, &vv)) in cols.iter().zip(values).enumerate() {
                let stride = checked_geometry_mul("locating a SELL slot", k, c)?;
                let slot = checked_geometry_add("offsetting a SELL slot", chunk_ptr[ch], stride)?;
                let slot = checked_geometry_add("selecting a SELL lane", slot, lane)?;
                col_idx[slot] = cc;
                vals[slot] = vv;
            }
        }
        Ok(Sell {
            c,
            sigma,
            nrows,
            ncols: a.ncols(),
            perm,
            row_len,
            chunk_ptr,
            col_idx,
            vals,
        })
    }

    /// Exact (bitwise-lossless) expansion back to CSR: true row lengths are
    /// stored, so stored zeros and all values survive unchanged.
    #[must_use]
    pub fn to_csr(&self) -> Csr {
        let mut per_row: Vec<(Vec<usize>, Vec<f64>)> = vec![(Vec::new(), Vec::new()); self.nrows];
        for (pos, &orig) in self.perm.iter().enumerate() {
            let (ch, lane) = (pos / self.c, pos % self.c);
            let (cols, values) = &mut per_row[orig];
            for k in 0..self.row_len[pos] {
                cols.push(self.col_idx[self.chunk_ptr[ch] + k * self.c + lane]);
                values.push(self.vals[self.chunk_ptr[ch] + k * self.c + lane]);
            }
        }
        let mut row_ptr = vec![0usize; self.nrows + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        for (r, (cols, values)) in per_row.into_iter().enumerate() {
            col_idx.extend_from_slice(&cols);
            vals.extend_from_slice(&values);
            row_ptr[r + 1] = col_idx.len();
        }
        Csr::from_parts(self.nrows, self.ncols, row_ptr, col_idx, vals)
    }

    /// Rows.
    #[must_use]
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Columns.
    #[must_use]
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Physical slots including padding (the quantity the perf bead's
    /// padding-overhead metric is computed from).
    #[must_use]
    pub fn physical_slots(&self) -> usize {
        self.vals.len()
    }

    /// Chunk height C.
    #[must_use]
    pub fn chunk_height(&self) -> usize {
        self.c
    }

    /// Sorting window σ (as rounded internally).
    #[must_use]
    pub fn sorting_window(&self) -> usize {
        self.sigma
    }

    /// y = A·x, bit-identical to CSR SpMV: each lane accumulates its row's
    /// entries in ascending-column order with fused mul_add, iterating to the
    /// TRUE row length (pads never read).
    pub fn spmv(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(
            x.len(),
            self.ncols,
            "spmv: x length must equal ncols {}",
            self.ncols
        );
        assert_eq!(
            y.len(),
            self.nrows,
            "spmv: y length must equal nrows {}",
            self.nrows
        );
        crate::fma::sell_spmv_dispatch(self, x, y);
    }

    /// Y = A·X for row-major X (ncols × nrhs) and Y (nrows × nrhs).
    /// Reads the true row lengths directly from SELL storage, never padding;
    /// sorted rows are scattered back to their original positions. Each
    /// right-hand side follows CSR's ascending-column FMA accumulation.
    /// Uses eight stack accumulators and the existing runtime FMA dispatch.
    pub fn spmm_blocked(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        let x_len = self.ncols.checked_mul(nrhs).expect("spmm: X size overflow");
        let y_len = self.nrows.checked_mul(nrhs).expect("spmm: Y size overflow");
        assert_eq!(x.len(), x_len, "spmm: X is ncols x nrhs");
        assert_eq!(y.len(), y_len, "spmm: Y is nrows x nrhs");
        crate::fma::sell_spmm_blocked_dispatch(self, x, nrhs, y);
    }

    /// RHS-tiled SELL body, inlined into the existing FMA capsule.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_blocked_body(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        if nrhs == 0 || self.nrows == 0 {
            return;
        }
        const NB: usize = 8;
        let mut acc = [0.0f64; NB];
        for j0 in (0..nrhs).step_by(NB) {
            let width = NB.min(nrhs - j0);
            for (pos, &orig) in self.perm.iter().enumerate() {
                let (ch, lane) = (pos / self.c, pos % self.c);
                acc[..width].fill(0.0);
                for k in 0..self.row_len[pos] {
                    let slot = self.chunk_ptr[ch] + k * self.c + lane;
                    let x_offset = self.col_idx[slot] * nrhs + j0;
                    for (a, &v) in acc[..width]
                        .iter_mut()
                        .zip(&x[x_offset..x_offset + width])
                    {
                        *a = self.vals[slot].mul_add(v, *a);
                    }
                }
                let dst = orig * nrhs + j0;
                y[dst..dst + width].copy_from_slice(&acc[..width]);
            }
        }
    }

    /// Y = A * X for row-major X (ncols * nrhs) and Y (nrows * nrhs).
    ///
    /// Traverses SELL storage directly and writes rows in original matrix order.
    /// Each RHS follows [`Self::spmv`]'s ascending-column FMA sequence; neither
    /// ragged lanes nor row padding read X. No conversion or scratch allocation
    /// is performed. Shape/overflow errors panic before changing Y. Zero RHS
    /// count requires empty buffers; empty rows produce positive zero.
    pub fn spmm(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        let x_len = self
            .ncols
            .checked_mul(nrhs)
            .expect("spmm: RHS dimensions ncols*nrhs overflow usize");
        let y_len = self
            .nrows
            .checked_mul(nrhs)
            .expect("spmm: output dimensions nrows*nrhs overflow usize");
        assert_eq!(x.len(), x_len, "spmm: X must be ncols*nrhs");
        assert_eq!(y.len(), y_len, "spmm: Y must be nrows*nrhs");
        if nrhs == 0 {
            return;
        }
        crate::fma::sell_spmm_dispatch(self, x, nrhs, y);
    }

    /// Direct batched SELL body, inlined into the existing FMA capsule.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_body(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        for (pos, &original_row) in self.perm.iter().enumerate() {
            let (chunk, lane) = (pos / self.c, pos % self.c);
            let out = &mut y[original_row * nrhs..][..nrhs];
            out.fill(0.0);
            for entry in 0..self.row_len[pos] {
                let index = self.chunk_ptr[chunk] + entry * self.c + lane;
                let value = self.vals[index];
                let xrow = &x[self.col_idx[index] * nrhs..][..nrhs];
                for (output, &input) in out.iter_mut().zip(xrow) {
                    *output = value.mul_add(input, *output);
                }
            }
        }
    }

    /// The spmv loop body — recompiled by the x86 FMA-codegen capsule
    /// (bead nabk); MUST stay `inline(always)`.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmv_body(&self, x: &[f64], y: &mut [f64]) {
        for (pos, &orig) in self.perm.iter().enumerate() {
            let (ch, lane) = (pos / self.c, pos % self.c);
            let mut acc = 0.0f64;
            for k in 0..self.row_len[pos] {
                let idx = self.chunk_ptr[ch] + k * self.c + lane;
                acc = self.vals[idx].mul_add(x[self.col_idx[idx]], acc);
            }
            y[orig] = acc;
        }
    }
}

impl Sell {
    /// CHUNK-MAJOR SpMV (bead wsbf segment 2): processes each chunk's
    /// C lanes together — per k, the C (col, val) entries are
    /// CONTIGUOUS (lane-fastest storage), so the inner loop is the
    /// SIMD shape fs-tilelang's lane width names, and the C
    /// independent accumulator chains break the per-row FMA latency
    /// bound the roofline lane measured on the CSR kernels.
    ///
    /// DETERMINISM: per-lane accumulation is k-ascending `mul_add`
    /// from +0.0 — the same order as [`Sell::spmv`] — and each lane
    /// stops at its stored logical row length. Physical pad slots are
    /// never read, including after an accumulation has produced −0.0.
    /// Bitwise equality with the row-major kernel is GATED in
    /// conformance.
    pub fn spmv_chunked(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        crate::fma::sell_spmv_chunked_dispatch(self, x, y);
    }

    /// The chunk-major loop body — recompiled by the x86 FMA-codegen
    /// capsule (bead nabk); MUST stay `inline(always)`.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmv_chunked_body(&self, x: &[f64], y: &mut [f64]) {
        let c = self.c;
        let nchunks = self.chunk_ptr.len() - 1;
        if nchunks == 0 {
            return;
        }
        let mut acc = vec![0.0f64; c.min(self.nrows)];
        for ch in 0..nchunks {
            let base = self.chunk_ptr[ch];
            let kmax = (self.chunk_ptr[ch + 1] - base) / c;
            let row0 = ch.saturating_mul(c).min(self.nrows);
            let live = self.nrows.saturating_sub(row0).min(c);
            acc.fill(0.0);
            for k in 0..kmax {
                let off = base + k * c;
                let cols = &self.col_idx[off..off + c];
                let vals = &self.vals[off..off + c];
                for l in 0..live {
                    if k < self.row_len[row0 + l] {
                        acc[l] = vals[l].mul_add(x[cols[l]], acc[l]);
                    }
                }
            }
            for (l, &a) in acc.iter().enumerate().take(live) {
                y[self.perm[row0 + l]] = a;
            }
        }
    }

    /// The per-shard chunk-range body (staging pairs) — recompiled by
    /// the x86 FMA-codegen capsule (bead nabk); MUST stay
    /// `inline(always)`.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn shard_body(&self, x: &[f64], lo: usize, hi: usize) -> Vec<(usize, f64)> {
        if lo >= hi {
            return Vec::new();
        }
        let c = self.c;
        let row_lo = lo.saturating_mul(c).min(self.nrows);
        let row_hi = hi.saturating_mul(c).min(self.nrows);
        let mut out = Vec::with_capacity(row_hi - row_lo);
        let mut acc = vec![0.0f64; self.nrows.saturating_sub(row_lo).min(c)];
        for ch in lo..hi {
            let base = self.chunk_ptr[ch];
            let kmax = (self.chunk_ptr[ch + 1] - base) / c;
            let row0 = ch.saturating_mul(c).min(self.nrows);
            let live = self.nrows.saturating_sub(row0).min(c);
            acc.fill(0.0);
            for k in 0..kmax {
                let off = base + k * c;
                let cols = &self.col_idx[off..off + c];
                let vals = &self.vals[off..off + c];
                for l in 0..live {
                    if k < self.row_len[row0 + l] {
                        acc[l] = vals[l].mul_add(x[cols[l]], acc[l]);
                    }
                }
            }
            for (l, &a) in acc.iter().enumerate().take(live) {
                out.push((self.perm[row0 + l], a));
            }
        }
        out
    }

    /// Sharded chunk-major SpMV: disjoint chunk ranges on scoped
    /// threads; each thread writes only its own rows (perm within a
    /// chunk is thread-exclusive). Bitwise equal to [`Self::spmv_chunked`]
    /// at every thread count.
    pub fn spmv_chunked_sharded(&self, x: &[f64], y: &mut [f64], threads: usize) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        let nchunks = self.chunk_ptr.len() - 1;
        if nchunks == 0 {
            return;
        }
        let t = bounded_worker_count(threads, nchunks);
        let per = nchunks.div_ceil(t);
        // Each thread owns disjoint ROWS (chunks partition sorted rows;
        // perm is a bijection), so unsynchronized writes through a raw
        // pointer wrapper would be needed... instead: write into a
        // per-thread staging of (orig_row, value) pairs and apply
        // serially — the apply order is chunk-ascending, deterministic.
        let mut stages: Vec<Vec<(usize, f64)>> = Vec::with_capacity(t);
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(t);
            let mut lo = 0usize;
            while lo < nchunks {
                let hi = lo.saturating_add(per).min(nchunks);
                handles.push(scope.spawn(move || {
                    // Closures are separate codegen units: target_feature
                    // does NOT propagate in, so the shard body dispatches
                    // HERE (one detection per shard — bead nabk).
                    crate::fma::sell_shard_dispatch(self, x, lo, hi)
                }));
                lo = hi;
            }
            for h in handles {
                stages.push(h.join().expect("chunk worker"));
            }
        });
        for stage in stages {
            for (r, v) in stage {
                y[r] = v;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::{laplacian_2d, random_coo};

    #[test]
    fn round_trip_is_bitwise_lossless() {
        for (n, c, sigma) in [(50usize, 4usize, 16usize), (64, 8, 8), (7, 4, 32)] {
            let a = random_coo(n, n, 5, 0x5E11 + n as u64).assemble();
            let s = Sell::from_csr(&a, c, sigma).expect("test SELL geometry is representable");
            let back = s.to_csr();
            assert_eq!(a.nnz(), back.nnz());
            for r in 0..n {
                let (c1, v1) = a.row(r);
                let (c2, v2) = back.row(r);
                assert_eq!(c1, c2, "row {r} pattern changed (C={c}, sigma={sigma})");
                assert!(v1.iter().zip(v2).all(|(x, y)| x.to_bits() == y.to_bits()));
            }
        }
    }

    #[test]
    fn spmv_bitwise_equals_csr() {
        let mut seed = 31u64;
        let mut lcg = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        };
        // Mixed shapes including a ragged final chunk and heavy row-length
        // skew (dense row + empty rows) to stress the permutation.
        for (n, c, sigma) in [(81usize, 4usize, 16usize), (64, 8, 64), (33, 16, 16)] {
            let mut coo = random_coo(n, n, 4, 0xACE + n as u64);
            for col in 0..n {
                coo.push(n / 2, col, 0.25); // dense row
            }
            let a = coo.assemble();
            let s = Sell::from_csr(&a, c, sigma).expect("test SELL geometry is representable");
            let x: Vec<f64> = (0..n).map(|_| lcg()).collect();
            let (mut y1, mut y2) = (vec![0.0; n], vec![0.0; n]);
            a.spmv(&x, &mut y1);
            s.spmv(&x, &mut y2);
            for r in 0..n {
                assert_eq!(
                    y1[r].to_bits(),
                    y2[r].to_bits(),
                    "SELL(C={c},sigma={sigma}) diverged from CSR at row {r}"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"sell-bitwise\",\"verdict\":\"pass\",\"detail\":\"SELL spmv == CSR spmv bitwise, 3 configs incl ragged+skewed\"}}"
        );
    }

    #[test]
    fn spmm_matches_csr_bitwise_with_sorting_and_padding() {
        for (n, c, sigma, nrhs) in [
            (81usize, 4usize, 16usize, 1usize),
            (64, 8, 64, 7),
            (33, 16, 16, 13),
            (7, 4, 32, 0),
        ] {
            let mut coo = random_coo(n, n, 4, 0x5E11 + n as u64);
            for col in 0..n {
                coo.push(n / 2, col, 0.25);
            }
            let a = coo.assemble();
            let s = Sell::from_csr(&a, c, sigma).expect("valid SELL geometry");
            let x: Vec<f64> = (0..n * nrhs)
                .map(|k| (k as f64 * 0.013).cos())
                .collect();
            let (mut expected, mut actual) =
                (vec![0.0; n * nrhs], vec![0.0; n * nrhs]);
            a.spmm_blocked(&x, nrhs, &mut expected);
            s.spmm_blocked(&x, nrhs, &mut actual);
            assert!(
                actual.iter().zip(&expected).all(|(x, y)| x.to_bits() == y.to_bits()),
                "SELL C={c} sigma={sigma} nrhs={nrhs} must agree bitwise with CSR"
            );
        }
    }

    #[test]
    fn g5_spmm_matches_csr_and_portable_after_sorting() {
        let mut coo = crate::Coo::new(7, 9);
        coo.push(0, 7, 2.0);
        coo.push(2, 4, 3.0);
        coo.push(6, 0, -1.0);
        for col in 0..9 {
            coo.push(3, col, col as f64 * 0.125 - 0.5);
        }
        let a = coo.assemble();
        for (c, sigma) in [(1, 7), (4, 8), (16, 32)] {
            let s = Sell::from_csr(&a, c, sigma).unwrap();
            for nrhs in [0usize, 1, 2, 5, 17] {
                let x: Vec<f64> = (0..a.ncols() * nrhs)
                    .map(|i| (i % 19) as f64 * 0.125 - 1.0)
                    .collect();
                let mut expected = vec![0.0; a.nrows() * nrhs];
                let mut actual = vec![f64::NAN; expected.len()];
                let mut portable = actual.clone();
                let mut blocked = actual.clone();
                let mut blocked_portable = actual.clone();
                a.spmm(&x, nrhs, &mut expected);
                s.spmm(&x, nrhs, &mut actual);
                s.spmm_body(&x, nrhs, &mut portable);
                s.spmm_blocked(&x, nrhs, &mut blocked);
                s.spmm_blocked_body(&x, nrhs, &mut blocked_portable);
                for output in [&actual, &portable, &blocked, &blocked_portable] {
                    for (got, want) in output.iter().zip(&expected) {
                        assert_eq!(got.to_bits(), want.to_bits(), "SELL C={c}, RHS {nrhs}");
                    }
                }
            }
        }
    }

    #[test]
    fn g0_spmm_does_not_evaluate_padding() {
        let a = Csr::from_parts(
            3, 4, vec![0, 1, 1, 3], vec![0, 0, 3], vec![-f64::from_bits(1), 0.0, 2.0],
        );
        let s = Sell::from_csr(&a, 4, 8).unwrap();
        let x = [0.25, 1.0, f64::INFINITY, f64::NAN, f64::NAN, f64::NEG_INFINITY, 2.0, -3.0];
        let mut expected = [0.0; 6];
        let mut actual = [f64::NAN; 6];
        a.spmm(&x, 2, &mut expected);
        s.spmm(&x, 2, &mut actual);
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert_eq!(actual[0].to_bits(), (-0.0f64).to_bits());
        actual.fill(f64::NAN);
        s.spmm_blocked(&x, 2, &mut actual);
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }

    #[test]
    fn g0_spmm_handles_zero_axes() {
        let empty_columns = Csr::from_parts(3, 0, vec![0; 4], vec![], vec![]);
        let mut y = [f64::NAN; 9];
        Sell::from_csr(&empty_columns, 4, 8).unwrap().spmm(&[], 3, &mut y);
        assert!(y.iter().all(|v| v.to_bits() == 0));
        y.fill(f64::NAN);
        Sell::from_csr(&empty_columns, 4, 8).unwrap().spmm_blocked(&[], 3, &mut y);
        assert!(y.iter().all(|v| v.to_bits() == 0));
        let empty_rows = Csr::from_parts(0, 6, vec![0], vec![], vec![]);
        Sell::from_csr(&empty_rows, 4, 8).unwrap().spmm(&[f64::NAN; 18], 3, &mut []);
        Sell::from_csr(&empty_rows, 4, 8).unwrap().spmm_blocked(&[f64::NAN; 18], 3, &mut []);
        let empty = Sell::from_csr(&Csr::identity(0), 4, 8).unwrap();
        empty.spmm(&[], usize::MAX, &mut []);
        empty.spmm_blocked(&[], usize::MAX, &mut []);
        empty.spmm_blocked_body(&[], usize::MAX, &mut []);
    }

    #[test]
    fn g0_spmm_shape_refusals_leave_output_unchanged() {
        let s = Sell::from_csr(&Csr::identity(2), 4, 8).unwrap();
        for (x_len, nrhs, y_len) in [(3, 2, 4), (4, 2, 3), (0, usize::MAX, 4), (0, 0, 1)] {
            let x = vec![1.0; x_len];
            let mut y = vec![42.0; y_len];
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                s.spmm(&x, nrhs, &mut y);
            }));
            assert!(result.is_err());
            assert_eq!(y, vec![42.0; y_len]);
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                s.spmm_blocked(&x, nrhs, &mut y);
            }));
            assert!(result.is_err());
            assert_eq!(y, vec![42.0; y_len]);
        }
        let zero_columns = Csr::from_parts(2, 0, vec![0; 3], vec![], vec![]);
        let s = Sell::from_csr(&zero_columns, 4, 8).unwrap();
        let mut y = [42.0];
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            s.spmm(&[], usize::MAX, &mut y);
        }));
        assert!(result.is_err(), "output dimension overflow must be rejected");
        assert_eq!(y, [42.0]);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            s.spmm_blocked(&[], usize::MAX, &mut y);
        }));
        assert!(result.is_err(), "blocked output dimension overflow must be rejected");
        assert_eq!(y, [42.0]);
    }

    #[test]
    fn sorting_reduces_padding_on_skewed_rows() {
        // The reason sigma exists: with one dense row per window, sorting
        // within the window bounds the wide chunk count.
        let n = 64;
        let mut coo = random_coo(n, n, 2, 7);
        for col in 0..n {
            coo.push(0, col, 1.0);
        }
        let a = coo.assemble();
        let unsorted = Sell::from_csr(&a, 8, 8).expect("test SELL geometry is representable"); // sigma == C → no sorting effect
        let sorted = Sell::from_csr(&a, 8, 64).expect("test SELL geometry is representable");
        assert!(
            sorted.physical_slots() <= unsorted.physical_slots(),
            "sorting must not increase padding: {} vs {}",
            sorted.physical_slots(),
            unsorted.physical_slots()
        );
        // Laplacian is near-uniform: padding overhead should be tiny.
        let lap = laplacian_2d(16);
        let s = Sell::from_csr(&lap, 8, 64).expect("test SELL geometry is representable");
        let overhead = s.physical_slots() as f64 / lap.nnz() as f64;
        assert!(
            overhead < 1.35,
            "Laplacian SELL overhead {overhead} unexpectedly high"
        );
    }
}

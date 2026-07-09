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
    #[must_use]
    pub fn from_csr(a: &Csr, c: usize, sigma: usize) -> Sell {
        assert!(c >= 1, "chunk height must be at least 1");
        let sigma = sigma.max(c).div_ceil(c) * c;
        let nrows = a.nrows();
        let mut perm: Vec<usize> = (0..nrows).collect();
        for window in perm.chunks_mut(sigma) {
            window.sort_by_key(|&r| std::cmp::Reverse(a.row(r).0.len()));
        }
        let nchunks = nrows.div_ceil(c);
        let row_len: Vec<usize> = perm.iter().map(|&r| a.row(r).0.len()).collect();
        let mut chunk_ptr = vec![0usize; nchunks + 1];
        for ch in 0..nchunks {
            let width = (ch * c..((ch + 1) * c).min(nrows))
                .map(|p| row_len[p])
                .max()
                .unwrap_or(0);
            chunk_ptr[ch + 1] = chunk_ptr[ch] + width * c;
        }
        let total = chunk_ptr[nchunks];
        let mut col_idx = vec![0usize; total];
        let mut vals = vec![0.0f64; total];
        for (pos, &orig) in perm.iter().enumerate() {
            let (ch, lane) = (pos / c, pos % c);
            let (cols, values) = a.row(orig);
            for (k, (&cc, &vv)) in cols.iter().zip(values).enumerate() {
                col_idx[chunk_ptr[ch] + k * c + lane] = cc;
                vals[chunk_ptr[ch] + k * c + lane] = vv;
            }
        }
        Sell {
            c,
            sigma,
            nrows,
            ncols: a.ncols(),
            perm,
            row_len,
            chunk_ptr,
            col_idx,
            vals,
        }
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
    /// from +0.0 — the same order as [`Sell::spmv`] — and lanes
    /// shorter than the chunk read their PAD slots (col 0, val +0.0):
    /// `mul_add(+0.0, x[0], acc)` is exactly `acc` when `acc` is not
    /// −0.0, and acc starts +0.0 and can only stay +0.0 under
    /// (+0.0) + (±0.0) — the sell.rs signed-zero argument, inherited
    /// here because the kernel reads pads. Bitwise equality with the
    /// row-major kernel is GATED in conformance.
    pub fn spmv_chunked(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        let c = self.c;
        let nchunks = self.chunk_ptr.len() - 1;
        let mut acc = vec![0.0f64; c];
        for ch in 0..nchunks {
            let base = self.chunk_ptr[ch];
            let kmax = (self.chunk_ptr[ch + 1] - base) / c;
            acc.fill(0.0);
            for k in 0..kmax {
                let off = base + k * c;
                let cols = &self.col_idx[off..off + c];
                let vals = &self.vals[off..off + c];
                for l in 0..c {
                    acc[l] = vals[l].mul_add(x[cols[l]], acc[l]);
                }
            }
            let row0 = ch * c;
            for (l, &a) in acc
                .iter()
                .enumerate()
                .take(self.nrows.saturating_sub(row0).min(c))
            {
                y[self.perm[row0 + l]] = a;
            }
        }
    }

    /// Sharded chunk-major SpMV: disjoint chunk ranges on scoped
    /// threads; each thread writes only its own rows (perm within a
    /// chunk is thread-exclusive). Bitwise equal to [`Self::spmv_chunked`]
    /// at every thread count.
    pub fn spmv_chunked_sharded(&self, x: &[f64], y: &mut [f64], threads: usize) {
        assert_eq!(x.len(), self.ncols, "spmv: x length");
        assert_eq!(y.len(), self.nrows, "spmv: y length");
        let t = threads.max(1);
        let nchunks = self.chunk_ptr.len() - 1;
        let per = nchunks.div_ceil(t);
        // Each thread owns disjoint ROWS (chunks partition sorted rows;
        // perm is a bijection), so unsynchronized writes through a raw
        // pointer wrapper would be needed... instead: write into a
        // per-thread staging of (orig_row, value) pairs and apply
        // serially — the apply order is chunk-ascending, deterministic.
        let mut stages: Vec<Vec<(usize, f64)>> = Vec::new();
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for w in 0..t {
                let lo = (w * per).min(nchunks);
                let hi = ((w + 1) * per).min(nchunks);
                handles.push(scope.spawn(move || {
                    let c = self.c;
                    let mut out = Vec::with_capacity((hi - lo) * c);
                    let mut acc = vec![0.0f64; c];
                    for ch in lo..hi {
                        let base = self.chunk_ptr[ch];
                        let kmax = (self.chunk_ptr[ch + 1] - base) / c;
                        acc.fill(0.0);
                        for k in 0..kmax {
                            let off = base + k * c;
                            let cols = &self.col_idx[off..off + c];
                            let vals = &self.vals[off..off + c];
                            for l in 0..c {
                                acc[l] = vals[l].mul_add(x[cols[l]], acc[l]);
                            }
                        }
                        let row0 = ch * c;
                        let live = self.nrows.saturating_sub(row0).min(c);
                        for (l, &a) in acc.iter().enumerate().take(live) {
                            out.push((self.perm[row0 + l], a));
                        }
                    }
                    out
                }));
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
            let s = Sell::from_csr(&a, c, sigma);
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
            let s = Sell::from_csr(&a, c, sigma);
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
    fn sorting_reduces_padding_on_skewed_rows() {
        // The reason sigma exists: with one dense row per window, sorting
        // within the window bounds the wide chunk count.
        let n = 64;
        let mut coo = random_coo(n, n, 2, 7);
        for col in 0..n {
            coo.push(0, col, 1.0);
        }
        let a = coo.assemble();
        let unsorted = Sell::from_csr(&a, 8, 8); // sigma == C → no sorting effect
        let sorted = Sell::from_csr(&a, 8, 64);
        assert!(
            sorted.physical_slots() <= unsorted.physical_slots(),
            "sorting must not increase padding: {} vs {}",
            sorted.physical_slots(),
            unsorted.physical_slots()
        );
        // Laplacian is near-uniform: padding overhead should be tiny.
        let lap = laplacian_2d(16);
        let s = Sell::from_csr(&lap, 8, 64);
        let overhead = s.physical_slots() as f64 / lap.nnz() as f64;
        assert!(
            overhead < 1.35,
            "Laplacian SELL overhead {overhead} unexpectedly high"
        );
    }
}

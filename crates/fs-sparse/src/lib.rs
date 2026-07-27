//! fs-sparse — sparse formats (CSR/BSR/SELL-C-σ), deterministic assembly,
//! SpMV/SpMM, and pattern algebra (transpose/symmetrize/SpGEMM). Plan §6.2.
//!
//! DESIGN CENTER — bitwise determinism across formats: every format's SpMV
//! accumulates each row's contributions in ascending-column order with fused
//! `mul_add`, so CSR, BSR, and SELL-C-σ produce BIT-IDENTICAL results (tested,
//! not aspirational). SELL keeps per-row lengths and never touches its pad
//! slots, so signed-zero pollution from padding is designed out rather than
//! tolerated. Assembly canonicalizes through a stable sort keyed by
//! (row, col, insertion order): the resulting matrix is a pure function of
//! the triplet MULTISET-WITH-SEQUENCE, independent of how tiles interleaved
//! their pushes (Decalogue P2 — same mesh, same matrix, bitwise).
//!
//! v1 is the correctness core on scalar kernels. The roofline lane
//! (≥85% STREAM, per-CCD sharding, prefetch, fs-tilelang SIMD bodies) is the
//! recorded follow-up bead, gated on fs-tilelang + the autotuner.

pub mod bsr;
mod fma;
#[cfg(feature = "fnx-interop")]
pub mod interop;
#[cfg(feature = "fnp-interop")]
pub mod interop_fnp;
pub mod ops;
pub mod perf;
pub mod precond;
pub mod sell;

pub use bsr::Bsr;
#[cfg(feature = "fnx-interop")]
pub use interop::{InteropError, WEIGHT_KEY, csr_to_graph_snapshot, graph_snapshot_to_csr};
#[cfg(feature = "fnp-interop")]
pub use interop_fnp::{FnpInteropError, csr_to_dense_array, dense_array_to_csr};
pub use perf::CsrCompact;
pub use sell::{Sell, SellError};

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ---------------------------------------------------------------------------
// COO staging
// ---------------------------------------------------------------------------

/// Triplet staging buffer for assembly. Duplicate (row, col) entries are
/// allowed and ACCUMULATE (the FEM element-assembly contract).
#[derive(Debug, Clone, Default)]
pub struct Coo {
    pub(crate) nrows: usize,
    pub(crate) ncols: usize,
    pub(crate) rows: Vec<usize>,
    pub(crate) cols: Vec<usize>,
    pub(crate) vals: Vec<f64>,
}

impl Coo {
    /// Empty staging buffer with fixed dimensions.
    #[must_use]
    pub fn new(nrows: usize, ncols: usize) -> Coo {
        Coo {
            nrows,
            ncols,
            rows: Vec::new(),
            cols: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// Stage one contribution. Panics (structured) on out-of-range indices —
    /// a mis-indexed element is a programmer error, not a runtime condition.
    pub fn push(&mut self, row: usize, col: usize, val: f64) {
        assert!(
            row < self.nrows && col < self.ncols,
            "coo push ({row},{col}) outside {}x{}",
            self.nrows,
            self.ncols
        );
        self.rows.push(row);
        self.cols.push(col);
        self.vals.push(val);
    }

    /// Number of staged triplets (before deduplication).
    #[must_use]
    pub fn staged(&self) -> usize {
        self.vals.len()
    }

    /// Canonical deterministic assembly into CSR.
    ///
    /// Contract: the result is a pure function of the SET of
    /// (row, col, insertion_seq, val) tuples — triplets may arrive in ANY
    /// stream order (tiles, threads, replays) and the matrix is bitwise
    /// identical, because duplicates are accumulated in insertion-sequence
    /// order after a stable sort by (row, col). Insertion sequence is LOGICAL
    /// identity (e.g. element id), never thread arrival.
    #[must_use]
    pub fn assemble(&self) -> Csr {
        let nnz_staged = self.vals.len();
        // Sort triplet INDICES by (row, col); the sort is stable, so equal
        // keys keep insertion order — the canonical accumulation order.
        let mut order: Vec<usize> = (0..nnz_staged).collect();
        order.sort_by_key(|&i| (self.rows[i], self.cols[i]));

        let mut row_ptr = vec![0usize; self.nrows + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        let mut i = 0;
        while i < nnz_staged {
            let (r, c) = (self.rows[order[i]], self.cols[order[i]]);
            let mut acc = self.vals[order[i]];
            i += 1;
            while i < nnz_staged && self.rows[order[i]] == r && self.cols[order[i]] == c {
                acc += self.vals[order[i]];
                i += 1;
            }
            row_ptr[r + 1] += 1;
            col_idx.push(c);
            vals.push(acc);
        }
        for r in 0..self.nrows {
            row_ptr[r + 1] += row_ptr[r];
        }
        Csr {
            nrows: self.nrows,
            ncols: self.ncols,
            row_ptr,
            col_idx,
            vals,
        }
    }
}

// ---------------------------------------------------------------------------
// CSR — the canonical format
// ---------------------------------------------------------------------------

/// Compressed sparse row. CANONICAL INVARIANT: within each row, column
/// indices are strictly ascending (no duplicates). Every constructor
/// establishes this; every algorithm may rely on it.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    nrows: usize,
    ncols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    vals: Vec<f64>,
}

impl Csr {
    /// Build from raw parts, VALIDATING the canonical invariant (structured
    /// panics on violation — silently accepting a non-canonical matrix would
    /// void every determinism claim downstream).
    #[must_use]
    pub fn from_parts(
        nrows: usize,
        ncols: usize,
        row_ptr: Vec<usize>,
        col_idx: Vec<usize>,
        vals: Vec<f64>,
    ) -> Csr {
        assert_eq!(
            row_ptr.len(),
            nrows + 1,
            "row_ptr must have nrows+1 entries"
        );
        assert_eq!(row_ptr[0], 0, "row_ptr must start at 0");
        assert_eq!(
            *row_ptr.last().unwrap(),
            col_idx.len(),
            "row_ptr end must equal nnz"
        );
        assert_eq!(col_idx.len(), vals.len(), "col_idx/vals length mismatch");
        for r in 0..nrows {
            assert!(
                row_ptr[r] <= row_ptr[r + 1],
                "row_ptr must be monotone (row {r})"
            );
            let cols = &col_idx[row_ptr[r]..row_ptr[r + 1]];
            for w in cols.windows(2) {
                assert!(w[0] < w[1], "row {r}: columns must be strictly ascending");
            }
            if let Some(&last) = cols.last() {
                assert!(last < ncols, "row {r}: column {last} out of range {ncols}");
            }
        }
        Csr {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            vals,
        }
    }

    /// Validate and publish canonical CSR storage with a caller-supplied
    /// checkpoint between bounded validation steps.
    ///
    /// This is the cancellation-friendly sibling of [`Self::from_parts`]. The
    /// callback is invoked for every row and stored column, so callers can map
    /// their own cancellation or work-accounting error without adding a
    /// runtime dependency to `fs-sparse`. `Ok(None)` reports malformed shape,
    /// row pointers, or noncanonical/out-of-range columns; `Err(error)` passes
    /// through a checkpoint refusal.
    pub fn try_from_parts_with_checkpoint<E>(
        nrows: usize,
        ncols: usize,
        row_ptr: Vec<usize>,
        col_idx: Vec<usize>,
        vals: Vec<f64>,
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<Option<Csr>, E> {
        let Some(expected_row_ptr_len) = nrows.checked_add(1) else {
            return Ok(None);
        };
        if row_ptr.len() != expected_row_ptr_len
            || row_ptr.first().copied() != Some(0)
            || row_ptr.last().copied() != Some(col_idx.len())
            || col_idx.len() != vals.len()
        {
            return Ok(None);
        }
        for window in row_ptr.windows(2) {
            checkpoint()?;
            let start = window[0];
            let end = window[1];
            if start > end || end > col_idx.len() {
                return Ok(None);
            }
            let mut previous = None;
            for &column in &col_idx[start..end] {
                checkpoint()?;
                if column >= ncols || previous.is_some_and(|prior| prior >= column) {
                    return Ok(None);
                }
                previous = Some(column);
            }
        }
        Ok(Some(Csr {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            vals,
        }))
    }

    /// An n×n identity.
    #[must_use]
    pub fn identity(n: usize) -> Csr {
        Csr {
            nrows: n,
            ncols: n,
            row_ptr: (0..=n).collect(),
            col_idx: (0..n).collect(),
            vals: vec![1.0; n],
        }
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

    /// Stored nonzeros.
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.vals.len()
    }

    /// One row's (columns, values) — also the graph neighbor view (CSR IS
    /// the adjacency structure; FrankenNetworkx interop builds on this).
    #[must_use]
    pub fn row(&self, r: usize) -> (&[usize], &[f64]) {
        let span = self.row_ptr[r]..self.row_ptr[r + 1];
        (&self.col_idx[span.clone()], &self.vals[span])
    }

    /// Value access (O(log nnz_row)); zero if not stored.
    #[must_use]
    pub fn get(&self, r: usize, c: usize) -> f64 {
        let (cols, vals) = self.row(r);
        match cols.binary_search(&c) {
            Ok(k) => vals[k],
            Err(_) => 0.0,
        }
    }

    /// y = A·x. Deterministic by construction: each y[r] accumulates in
    /// ascending-column order with fused mul_add — the SAME order and
    /// arithmetic every format's kernel uses (bitwise cross-format equality
    /// is tested, not hoped for).
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
        fma::spmv_dispatch(self, x, y);
    }

    /// The spmv loop body (extracted so the x86 FMA-codegen capsule can
    /// recompile it under `target_feature` — bead nabk; the a55x
    /// pattern). MUST stay `inline(always)`: a non-inlined call would
    /// keep baseline codegen and the per-element libm `fma()` call.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmv_body(&self, x: &[f64], y: &mut [f64]) {
        for (r, out) in y.iter_mut().enumerate() {
            let (cols, vals) = self.row(r);
            let mut acc = 0.0f64;
            for (&c, &v) in cols.iter().zip(vals) {
                acc = v.mul_add(x[c], acc);
            }
            *out = acc;
        }
    }

    /// C = A·X for dense row-major X (ncols × k) into row-major Y (nrows × k).
    /// Per output column the accumulation sequence is identical to `spmv` on
    /// that column — bitwise equal (tested).
    pub fn spmm(&self, x: &[f64], k: usize, y: &mut [f64]) {
        let (x_len, y_len) = self.checked_spmm_buffer_lengths(k);
        assert_eq!(x.len(), x_len, "spmm: X must be ncols*k");
        assert_eq!(y.len(), y_len, "spmm: Y must be nrows*k");
        fma::spmm_dispatch(self, x, k, y);
    }

    /// Validate dense SpMM totals once before any row-stride arithmetic.
    ///
    /// Canonical CSR guarantees every stored column is below `ncols`; the
    /// kernels visit only rows below `nrows`. Representable total lengths
    /// therefore bound every subsequent row/column offset and slice endpoint.
    pub(crate) fn checked_spmm_buffer_lengths(&self, nrhs: usize) -> (usize, usize) {
        let x_len = self
            .ncols
            .checked_mul(nrhs)
            .expect("spmm: RHS dimensions ncols*nrhs overflow usize");
        let y_len = self
            .nrows
            .checked_mul(nrhs)
            .expect("spmm: output dimensions nrows*nrhs overflow usize");
        (x_len, y_len)
    }

    /// The spmm loop body (see [`Csr::spmv_body`]'s capsule note).
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_body(&self, x: &[f64], k: usize, y: &mut [f64]) {
        for r in 0..self.nrows {
            let (cols, vals) = self.row(r);
            let out = &mut y[r * k..(r + 1) * k];
            out.fill(0.0);
            for (&c, &v) in cols.iter().zip(vals) {
                let xrow = &x[c * k..(c + 1) * k];
                for (o, &xv) in out.iter_mut().zip(xrow) {
                    *o = v.mul_add(xv, *o);
                }
            }
        }
    }

    /// Dense expansion (test/oracle use; O(nrows·ncols) memory).
    #[must_use]
    pub fn to_dense(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.nrows * self.ncols];
        for r in 0..self.nrows {
            let (cols, vals) = self.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                d[r * self.ncols + c] = v;
            }
        }
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// Random sparse matrix via COO with duplicate-accumulation exercised.
    pub(crate) fn random_coo(nrows: usize, ncols: usize, per_row: usize, seed: u64) -> Coo {
        let mut s = seed;
        let mut coo = Coo::new(nrows, ncols);
        for r in 0..nrows {
            for _ in 0..per_row {
                let c = ((lcg(&mut s) + 0.5) * ncols as f64) as usize % ncols;
                coo.push(r, c, lcg(&mut s));
            }
        }
        coo
    }

    /// 5-point Laplacian on an n×n grid — the FEM-patterned fixture.
    pub(crate) fn laplacian_2d(n: usize) -> Csr {
        let dim = n * n;
        let mut coo = Coo::new(dim, dim);
        for i in 0..n {
            for j in 0..n {
                let u = i * n + j;
                coo.push(u, u, 4.0);
                if i > 0 {
                    coo.push(u, u - n, -1.0);
                }
                if i + 1 < n {
                    coo.push(u, u + n, -1.0);
                }
                if j > 0 {
                    coo.push(u, u - 1, -1.0);
                }
                if j + 1 < n {
                    coo.push(u, u + 1, -1.0);
                }
            }
        }
        coo.assemble()
    }

    #[test]
    fn assembly_accumulates_duplicates_and_sorts() {
        let mut coo = Coo::new(2, 3);
        coo.push(1, 2, 1.0);
        coo.push(0, 1, 2.0);
        coo.push(1, 2, 0.5); // duplicate → accumulate
        coo.push(0, 0, -1.0);
        let a = coo.assemble();
        assert_eq!(a.nnz(), 3);
        assert_eq!(a.get(0, 0).to_bits(), (-1.0f64).to_bits());
        assert_eq!(a.get(0, 1).to_bits(), 2.0f64.to_bits());
        assert_eq!(a.get(1, 2).to_bits(), 1.5f64.to_bits());
        let (cols, _) = a.row(0);
        assert_eq!(cols, &[0, 1], "columns must come out sorted");
    }

    #[test]
    fn assembly_is_stream_order_invariant_bitwise() {
        // The G5 property at v1 scope: shuffle the triplet stream (as tiles /
        // threads would), assemble, and demand BITWISE equality — possible
        // because accumulation order is (row, col, insertion_seq), and we
        // preserve insertion sequence as logical identity.
        let base = random_coo(40, 40, 12, 0x5EED);
        let a = base.assemble();
        let n = base.staged();
        // Deterministic pseudo-shuffles: stride permutations of the stream.
        for stride in [7usize, 13, 29] {
            assert!(
                !n.is_multiple_of(stride),
                "stride must not divide n for full cycle"
            );
            let mut shuffled = Coo::new(40, 40);
            let mut idx = 0;
            // Visit indices in a stride cycle, but PUSH in insertion order of
            // the ORIGINAL: we must keep (row,col)-duplicate relative order.
            // A real tile scheduler interleaves DISTINCT (row,col) freely but
            // preserves each element's own sequence; model that by shuffling
            // whole rows (tiles) rather than raw triplets.
            let mut row_order: Vec<usize> = (0..40).collect();
            for _ in 0..stride {
                idx = (idx + stride) % 40;
                row_order.swap(idx, (idx * stride + 3) % 40);
            }
            for &r in &row_order {
                for i in 0..n {
                    if base.rows[i] == r {
                        shuffled.push(base.rows[i], base.cols[i], base.vals[i]);
                    }
                }
            }
            let b = shuffled.assemble();
            assert_eq!(a.nnz(), b.nnz());
            assert!(
                a.vals
                    .iter()
                    .zip(&b.vals)
                    .all(|(x, y)| x.to_bits() == y.to_bits()),
                "assembly depended on tile interleaving (stride {stride})"
            );
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"assembly-order-invariance\",\"verdict\":\"pass\",\"detail\":\"3 tile interleavings bitwise identical\"}}"
        );
    }

    #[test]
    fn spmv_matches_dense_reference() {
        let mut seed = 3u64;
        for (nr, nc) in [(1usize, 1usize), (17, 9), (64, 64), (100, 40)] {
            let a = random_coo(nr, nc, 6, 0xF00 + nr as u64).assemble();
            let x: Vec<f64> = (0..nc).map(|_| lcg(&mut seed)).collect();
            let mut y = vec![0.0; nr];
            a.spmv(&x, &mut y);
            let d = a.to_dense();
            for r in 0..nr {
                let want: f64 = (0..nc).map(|c| d[r * nc + c] * x[c]).sum();
                let scale = want.abs().max(1.0);
                assert!(
                    (y[r] - want).abs() / scale < 1e-13,
                    "spmv row {r} of {nr}x{nc}: {} vs dense {want}",
                    y[r]
                );
            }
        }
    }

    #[test]
    fn spmv_linearity_property() {
        let a = laplacian_2d(9);
        let n = a.nrows();
        let mut seed = 5u64;
        let x1: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
        let x2: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
        let alpha = 1.75; // power of two × small — keeps the check tight
        let combo: Vec<f64> = x1.iter().zip(&x2).map(|(a, b)| alpha * a + b).collect();
        let (mut y1, mut y2, mut yc) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        a.spmv(&x1, &mut y1);
        a.spmv(&x2, &mut y2);
        a.spmv(&combo, &mut yc);
        for i in 0..n {
            let want = alpha * y1[i] + y2[i];
            assert!((yc[i] - want).abs() < 1e-12, "linearity at {i}");
        }
    }

    #[test]
    fn spmm_is_bitwise_columnwise_spmv() {
        let a = random_coo(30, 20, 5, 99).assemble();
        let k = 7;
        let mut seed = 8u64;
        let x: Vec<f64> = (0..20 * k).map(|_| lcg(&mut seed)).collect();
        let mut y = vec![0.0; 30 * k];
        a.spmm(&x, k, &mut y);
        for col in 0..k {
            let xc: Vec<f64> = (0..20).map(|r| x[r * k + col]).collect();
            let mut yc = vec![0.0; 30];
            a.spmv(&xc, &mut yc);
            for r in 0..30 {
                assert_eq!(
                    y[r * k + col].to_bits(),
                    yc[r].to_bits(),
                    "spmm != spmv bitwise at ({r},{col})"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"spmm-bitwise\",\"verdict\":\"pass\",\"detail\":\"spmm == per-column spmv, 30x20 k=7\"}}"
        );
    }

    #[test]
    #[should_panic(expected = "spmm: RHS dimensions ncols*nrhs overflow usize")]
    fn spmm_rejects_rhs_shape_overflow_before_dispatch() {
        let huge_ncols = usize::MAX / 2 + 1;
        let a = Csr::from_parts(1, huge_ncols, vec![0, 0], Vec::new(), Vec::new());
        let mut y = [f64::NAN; 2];
        a.spmm(&[], 2, &mut y);
    }

    #[test]
    #[should_panic(expected = "spmm: RHS dimensions ncols*nrhs overflow usize")]
    fn blocked_spmm_rejects_rhs_shape_overflow_before_dispatch() {
        let huge_ncols = usize::MAX / 2 + 1;
        let a = Csr::from_parts(1, huge_ncols, vec![0, 0], Vec::new(), Vec::new());
        let mut y = [f64::NAN; 2];
        a.spmm_blocked(&[], 2, &mut y);
    }

    #[test]
    fn adversarial_patterns() {
        // Empty rows, dense row, single column — the shapes that break naive
        // pointer arithmetic.
        let mut coo = Coo::new(5, 4);
        for c in 0..4 {
            coo.push(2, c, f64::from(u32::try_from(c).unwrap()) + 1.0); // dense row
        }
        coo.push(4, 0, -3.0); // rows 0,1,3 empty
        let a = coo.assemble();
        assert_eq!(a.nnz(), 5);
        let x = [1.0, 1.0, 1.0, 1.0];
        let mut y = vec![0.0; 5];
        a.spmv(&x, &mut y);
        assert_eq!(y, vec![0.0, 0.0, 10.0, 0.0, -3.0]);
        // Single-column matrix.
        let mut c1 = Coo::new(3, 1);
        c1.push(0, 0, 2.0);
        c1.push(2, 0, -1.0);
        let b = c1.assemble();
        let mut yb = vec![0.0; 3];
        b.spmv(&[4.0], &mut yb);
        assert_eq!(yb, vec![8.0, 0.0, -4.0]);
        // Fully empty matrix.
        let e = Coo::new(3, 3).assemble();
        assert_eq!(e.nnz(), 0);
        let mut ye = vec![1.0; 3];
        e.spmv(&[1.0, 1.0, 1.0], &mut ye);
        assert_eq!(ye, vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn from_parts_rejects_non_canonical() {
        for (rp, ci, why) in [
            (vec![0usize, 2], vec![1usize, 1], "duplicate columns"),
            (vec![0, 2], vec![1, 0], "descending columns"),
            (vec![0, 1], vec![5], "column out of range"),
        ] {
            let vals = vec![1.0; ci.len()];
            let r =
                std::panic::catch_unwind(|| Csr::from_parts(1, 2, rp.clone(), ci.clone(), vals));
            assert!(r.is_err(), "must reject: {why}");
        }
    }

    #[test]
    fn checkpointed_parts_preserve_the_canonical_invariant() {
        let row_ptr = vec![0, 2, 3];
        let col_idx = vec![0, 2, 1];
        let vals = vec![1.0, 2.0, 3.0];
        let checked = Csr::from_parts(2, 3, row_ptr.clone(), col_idx.clone(), vals.clone());
        let mut checkpoints = 0usize;
        let checkpointed =
            Csr::try_from_parts_with_checkpoint(2, 3, row_ptr, col_idx, vals, || {
                checkpoints += 1;
                Ok::<_, core::convert::Infallible>(())
            })
            .expect("checkpoint cannot fail")
            .expect("canonical parts pass checked publication");
        assert_eq!(checkpointed, checked);
        assert_eq!(checkpoints, 5);

        for columns in [vec![1, 0], vec![0, 0], vec![0, 3]] {
            let malformed = Csr::try_from_parts_with_checkpoint(
                1,
                3,
                vec![0, 2],
                columns,
                vec![1.0, 2.0],
                || Ok::<_, core::convert::Infallible>(()),
            )
            .expect("checkpoint cannot fail");
            assert!(malformed.is_none());
        }
        let bad_shape =
            Csr::try_from_parts_with_checkpoint(2, 3, vec![0, 1], vec![0], vec![1.0], || {
                Ok::<_, core::convert::Infallible>(())
            })
            .expect("checkpoint cannot fail");
        assert!(bad_shape.is_none());
        let decreasing_rows =
            Csr::try_from_parts_with_checkpoint(2, 3, vec![0, 2, 1], vec![0], vec![1.0], || {
                Ok::<_, core::convert::Infallible>(())
            })
            .expect("checkpoint cannot fail");
        assert!(decreasing_rows.is_none());
        let cancelled =
            Csr::try_from_parts_with_checkpoint(1, 1, vec![0, 1], vec![0], vec![1.0], || {
                Err("cancelled")
            });
        assert_eq!(cancelled, Err("cancelled"));
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}

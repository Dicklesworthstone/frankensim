//! Block CSR — FEM vector-valued unknowns make dense r×c micro-blocks;
//! BSR stores them contiguously so the future SIMD lane gets
//! batched-small-dense-adjacent work (plan §6.2).
//!
//! Bitwise note: BSR visits stored entries in ascending global column
//! order, using the same fused accumulation as CSR. An occupancy mask
//! distinguishes padding from explicitly stored zeros: evaluating padding
//! would read unrelated NaNs/infinities and can change a negative zero
//! produced by underflow. Padding never participates in arithmetic.

use crate::Csr;

/// Block compressed sparse row with fixed r×c blocks (row-major within a
/// block). Block columns are strictly ascending within each block row.
#[derive(Debug, Clone, PartialEq)]
pub struct Bsr {
    br: usize,
    bc: usize,
    nrows: usize,
    ncols: usize,
    brow_ptr: Vec<usize>,
    bcol_idx: Vec<usize>,
    blocks: Vec<f64>,
    /// True only for entries present in the source CSR, including +/-0.0.
    stored: Vec<bool>,
}

impl Bsr {
    /// Convert from canonical CSR. Panics (structured) unless dimensions are
    /// divisible by the block shape — padding a matrix is a MODELING decision
    /// the caller must make explicitly, not something a format conversion
    /// invents.
    #[must_use]
    pub fn from_csr(a: &Csr, br: usize, bc: usize) -> Bsr {
        assert!(br >= 1 && bc >= 1, "block shape must be at least 1x1");
        assert!(
            a.nrows().is_multiple_of(br) && a.ncols().is_multiple_of(bc),
            "matrix {}x{} not divisible by block {}x{}",
            a.nrows(),
            a.ncols(),
            br,
            bc
        );
        let block_area = br.checked_mul(bc).expect("bsr: block area overflows usize");
        let nbrows = a.nrows() / br;
        let mut brow_ptr = vec![0usize; nbrows + 1];
        let mut bcol_idx: Vec<usize> = Vec::new();
        let mut blocks: Vec<f64> = Vec::new();
        let mut stored: Vec<bool> = Vec::new();
        for rb in 0..nbrows {
            // Which block columns appear in this block row? Merge the rows'
            // column sets; ascending order falls out of a simple scan.
            let mut present: Vec<usize> = Vec::new();
            for i in 0..br {
                let (cols, _) = a.row(rb * br + i);
                for &c in cols {
                    let cb = c / bc;
                    if let Err(pos) = present.binary_search(&cb) {
                        present.insert(pos, cb);
                    }
                }
            }
            let base_block = blocks.len();
            let new_slots = present
                .len()
                .checked_mul(block_area)
                .and_then(|count| base_block.checked_add(count))
                .expect("bsr: block storage overflows usize");
            blocks.resize(new_slots, 0.0);
            stored.resize(new_slots, false);
            for i in 0..br {
                let (cols, vals) = a.row(rb * br + i);
                for (&c, &v) in cols.iter().zip(vals) {
                    let cb = c / bc;
                    let slot = present.binary_search(&cb).expect("built above");
                    let index = base_block + slot * block_area + i * bc + (c % bc);
                    blocks[index] = v;
                    stored[index] = true;
                }
            }
            brow_ptr[rb + 1] = brow_ptr[rb] + present.len();
            bcol_idx.extend_from_slice(&present);
        }
        Bsr {
            br,
            bc,
            nrows: a.nrows(),
            ncols: a.ncols(),
            brow_ptr,
            bcol_idx,
            blocks,
            stored,
        }
    }

    /// Expand back to CSR, preserving the original pattern and value bits,
    /// including explicit signed zeros and non-finite values. Only padding
    /// introduced by blocking is dropped; occupancy is not inferred from value.
    #[must_use]
    pub fn to_csr(&self) -> Csr {
        let mut row_ptr = vec![0usize; self.nrows + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        for rb in 0..self.brow_ptr.len() - 1 {
            for i in 0..self.br {
                let r = rb * self.br + i;
                for slot in self.brow_ptr[rb]..self.brow_ptr[rb + 1] {
                    let cb = self.bcol_idx[slot];
                    for j in 0..self.bc {
                        let index = slot * self.br * self.bc + i * self.bc + j;
                        if self.stored[index] {
                            col_idx.push(cb * self.bc + j);
                            vals.push(self.blocks[index]);
                        }
                    }
                }
                row_ptr[r + 1] = col_idx.len();
            }
        }
        Csr::from_parts(self.nrows, self.ncols, row_ptr, col_idx, vals)
    }

    /// Rows (scalar).
    #[must_use]
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Columns (scalar).
    #[must_use]
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Stored blocks.
    #[must_use]
    pub fn nblocks(&self) -> usize {
        self.bcol_idx.len()
    }

    /// y = A·x, bit-identical to CSR SpMV (ascending-column accumulation with
    /// fused mul_add; padding is skipped before reading x).
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
        crate::fma::bsr_spmv_dispatch(self, x, y);
    }

    /// Y = A·X for row-major X (ncols × nrhs) and Y (nrows × nrhs).
    /// Traverses stored blocks directly, retaining scalar ascending-column
    /// accumulation and bitwise agreement with CSR's blocked SpMM. Padding
    /// is skipped before reading X. Uses eight stack accumulators and the
    /// same runtime FMA dispatch as the unblocked kernel.
    pub fn spmm_blocked(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        let x_len = self.ncols.checked_mul(nrhs).expect("spmm: X size overflow");
        let y_len = self.nrows.checked_mul(nrhs).expect("spmm: Y size overflow");
        assert_eq!(x.len(), x_len, "spmm: X is ncols x nrhs");
        assert_eq!(y.len(), y_len, "spmm: Y is nrows x nrhs");
        crate::fma::bsr_spmm_blocked_dispatch(self, x, nrhs, y);
    }

    /// RHS-tiled block-row body, inlined into the existing FMA capsule.
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
            for rb in 0..self.brow_ptr.len() - 1 {
                for i in 0..self.br {
                    acc[..width].fill(0.0);
                    for slot in self.brow_ptr[rb]..self.brow_ptr[rb + 1] {
                        let cb = self.bcol_idx[slot];
                        let base = slot * self.br * self.bc + i * self.bc;
                        let block_row = &self.blocks[base..][..self.bc];
                        for (k, &value) in block_row.iter().enumerate() {
                            if !self.stored[base + k] {
                                continue;
                            }
                            let col_offset = (cb * self.bc + k) * nrhs + j0;
                            for (a, &v) in acc[..width]
                                .iter_mut()
                                .zip(&x[col_offset..col_offset + width])
                            {
                                *a = value.mul_add(v, *a);
                            }
                        }
                    }
                    let dst = (rb * self.br + i) * nrhs + j0;
                    y[dst..dst + width].copy_from_slice(&acc[..width]);
                }
            }
        }
    }

    /// Y = A * X for row-major X (ncols * nrhs) and Y (nrows * nrhs).
    ///
    /// Uses the block storage directly without conversion or scratch allocation.
    /// Each RHS follows the same stored-column FMA sequence as [`Self::spmv`].
    /// Explicit zeros participate; padding never reads X. Shape/overflow
    /// errors panic before changing Y. Zero RHS count is admitted with empty
    /// buffers, and empty rows overwrite their output with positive zero.
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
        crate::fma::bsr_spmm_dispatch(self, x, nrhs, y);
    }

    /// Direct batched block-row body, inlined into the existing FMA capsule.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_body(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
        for rb in 0..self.brow_ptr.len() - 1 {
            for i in 0..self.br {
                let row = rb * self.br + i;
                let out = &mut y[row * nrhs..][..nrhs];
                out.fill(0.0);
                for slot in self.brow_ptr[rb]..self.brow_ptr[rb + 1] {
                    let column0 = self.bcol_idx[slot] * self.bc;
                    let base = slot * self.br * self.bc + i * self.bc;
                    for j in 0..self.bc {
                        if self.stored[base + j] {
                            let value = self.blocks[base + j];
                            let xrow = &x[(column0 + j) * nrhs..][..nrhs];
                            for (output, &input) in out.iter_mut().zip(xrow) {
                                *output = value.mul_add(input, *output);
                            }
                        }
                    }
                }
            }
        }
    }

    /// The spmv loop body, extracted so the x86 FMA-codegen capsule can
    /// recompile it under `target_feature` (bead nabk). MUST stay
    /// `inline(always)`: a non-inlined call keeps baseline codegen and
    /// the per-element libm `fma()` call.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmv_body(&self, x: &[f64], y: &mut [f64]) {
        for rb in 0..self.brow_ptr.len() - 1 {
            for i in 0..self.br {
                let mut acc = 0.0f64;
                for slot in self.brow_ptr[rb]..self.brow_ptr[rb + 1] {
                    let cb = self.bcol_idx[slot];
                    let base = slot * self.br * self.bc + i * self.bc;
                    let block_row = &self.blocks[base..][..self.bc];
                    for (j, &v) in block_row.iter().enumerate() {
                        if self.stored[base + j] {
                            acc = v.mul_add(x[cb * self.bc + j], acc);
                        }
                    }
                }
                y[rb * self.br + i] = acc;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::{laplacian_2d, random_coo};

    #[test]
    fn round_trip_dense_bitwise_and_structural() {
        // Semantically lossless: dense expansion identical (bitwise).
        let a = laplacian_2d(8); // 64x64, 2x2-blockable
        let b = Bsr::from_csr(&a, 2, 2);
        let back = b.to_csr();
        let (da, db) = (a.to_dense(), back.to_dense());
        assert!(
            da.iter().zip(&db).all(|(x, y)| x.to_bits() == y.to_bits()),
            "dense expansion changed through BSR"
        );
        // Structurally lossless when the pattern is block-aligned: a random
        // block-diagonal matrix reproduces the exact CSR.
        let mut coo = crate::Coo::new(12, 12);
        for blk in 0..4 {
            for i in 0..3 {
                for j in 0..3 {
                    coo.push(
                        blk * 3 + i,
                        blk * 3 + j,
                        f64::from(u32::try_from(blk * 9 + i * 3 + j).unwrap()) + 1.0,
                    );
                }
            }
        }
        let c = coo.assemble();
        let rt = Bsr::from_csr(&c, 3, 3).to_csr();
        assert_eq!(c.nnz(), rt.nnz());
        for r in 0..12 {
            let (c1, v1) = c.row(r);
            let (c2, v2) = rt.row(r);
            assert_eq!(c1, c2);
            assert!(v1.iter().zip(v2).all(|(a, b)| a.to_bits() == b.to_bits()));
        }
    }

    #[test]
    fn spmv_bitwise_equals_csr() {
        let mut seed = 21u64;
        let mut lcg = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
        };
        for (n, br, bc) in [(64usize, 2usize, 2usize), (64, 4, 4), (60, 3, 5)] {
            let a = if br == bc && n == 64 {
                laplacian_2d(8)
            } else {
                random_coo(n, n, 7, 0xB5 + n as u64).assemble()
            };
            let b = Bsr::from_csr(&a, br, bc);
            let x: Vec<f64> = (0..a.ncols()).map(|_| lcg()).collect();
            let (mut y1, mut y2) = (vec![0.0; a.nrows()], vec![0.0; a.nrows()]);
            a.spmv(&x, &mut y1);
            b.spmv(&x, &mut y2);
            for r in 0..a.nrows() {
                assert_eq!(
                    y1[r].to_bits(),
                    y2[r].to_bits(),
                    "BSR({br}x{bc}) diverged from CSR at row {r}"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"bsr-bitwise\",\"verdict\":\"pass\",\"detail\":\"BSR spmv == CSR spmv bitwise, 3 block shapes\"}}"
        );
    }

    #[test]
    fn spmm_matches_csr_bitwise() {
        for (n, br, bc, nrhs) in [
            (60usize, 3usize, 5usize, 1usize),
            (60, 4, 3, 7),
            (64, 2, 2, 13),
            (64, 4, 4, 0),
        ] {
            let a = random_coo(n, n, 7, 0xB57 + n as u64).assemble();
            let b = Bsr::from_csr(&a, br, bc);
            let x: Vec<f64> = (0..n * nrhs)
                .map(|k| (k as f64 * 0.013).sin())
                .collect();
            let (mut expected, mut actual) =
                (vec![0.0; n * nrhs], vec![0.0; n * nrhs]);
            a.spmm_blocked(&x, nrhs, &mut expected);
            b.spmm_blocked(&x, nrhs, &mut actual);
            assert!(
                actual.iter().zip(&expected).all(|(x, y)| x.to_bits() == y.to_bits()),
                "BSR {br}x{bc} nrhs={nrhs} must agree bitwise with CSR"
            );
        }
    }

    #[test]
    fn g0_padding_does_not_read_unrelated_nonfinite_rhs() {
        let a = Csr::from_parts(
            4, 6, vec![0, 1, 1, 2, 3], vec![1, 4, 5], vec![2.0, -3.0, 4.0],
        );
        let b = Bsr::from_csr(&a, 2, 3);
        let x = [f64::NAN, 3.0, f64::INFINITY, f64::NEG_INFINITY, 2.0, 1.0];
        let mut expected = [0.0; 4];
        let mut actual = [f64::NAN; 4];
        a.spmv(&x, &mut expected);
        b.spmv(&x, &mut actual);
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert_eq!(actual, [6.0, 0.0, -6.0, 4.0]);
    }

    #[test]
    fn g0_padding_preserves_negative_zero_from_underflow() {
        let a = Csr::from_parts(1, 2, vec![0, 1], vec![0], vec![-f64::from_bits(1)]);
        let b = Bsr::from_csr(&a, 1, 2);
        let mut expected = [1.0];
        let mut actual = [1.0];
        a.spmv(&[0.25, 1.0], &mut expected);
        b.spmv(&[0.25, 1.0], &mut actual);
        assert_eq!(expected[0].to_bits(), (-0.0f64).to_bits());
        assert_eq!(actual[0].to_bits(), expected[0].to_bits());
    }

    #[test]
    fn g0_round_trip_preserves_explicit_zeros_and_nonfinite_payloads() {
        let a = Csr::from_parts(
            2, 4, vec![0, 2, 4], vec![0, 3, 1, 2],
            vec![-0.0, f64::from_bits(0x7ff8_0000_0000_0042), 0.0, f64::INFINITY],
        );
        let back = Bsr::from_csr(&a, 2, 2).to_csr();
        assert_eq!(back.nnz(), a.nnz());
        for row in 0..a.nrows() {
            let (expected_cols, expected_vals) = a.row(row);
            let (actual_cols, actual_vals) = back.row(row);
            assert_eq!(actual_cols, expected_cols);
            for (actual, expected) in actual_vals.iter().zip(expected_vals) {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
        }
    }

    #[test]
    fn g0_explicit_zero_still_participates_in_arithmetic() {
        let a = Csr::from_parts(2, 2, vec![0, 1, 1], vec![1], vec![0.0]);
        let b = Bsr::from_csr(&a, 2, 2);
        let mut y = [1.0; 2];
        b.spmv(&[1.0, f64::INFINITY], &mut y);
        assert!(y[0].is_nan(), "stored 0 * infinity must not be skipped");
        assert_eq!(y[1].to_bits(), 0.0f64.to_bits(), "empty row is +0.0");
    }

    #[test]
    fn g5_spmm_matches_csr_and_portable_body() {
        let a = random_coo(12, 18, 5, 0xB5_5A).assemble();
        for (br, bc) in [(1, 1), (2, 3), (3, 2), (4, 6)] {
            let b = Bsr::from_csr(&a, br, bc);
            for nrhs in [0usize, 1, 3, 9, 17] {
                let x: Vec<f64> = (0..a.ncols() * nrhs)
                    .map(|i| (i % 23) as f64 * 0.125 - 1.0)
                    .collect();
                let mut expected = vec![0.0; a.nrows() * nrhs];
                let mut actual = vec![f64::NAN; expected.len()];
                let mut portable = actual.clone();
                let mut blocked = actual.clone();
                let mut blocked_portable = actual.clone();
                a.spmm(&x, nrhs, &mut expected);
                b.spmm(&x, nrhs, &mut actual);
                b.spmm_body(&x, nrhs, &mut portable);
                b.spmm_blocked(&x, nrhs, &mut blocked);
                b.spmm_blocked_body(&x, nrhs, &mut blocked_portable);
                for output in [&actual, &portable, &blocked, &blocked_portable] {
                    for (got, want) in output.iter().zip(&expected) {
                        assert_eq!(got.to_bits(), want.to_bits(), "BSR {br}x{bc}, RHS {nrhs}");
                    }
                }
            }
        }
    }

    #[test]
    fn g0_spmm_skips_padding_but_not_stored_zeros() {
        let a = Csr::from_parts(
            2, 4, vec![0, 1, 2], vec![0, 3], vec![-f64::from_bits(1), 2.0],
        );
        let b = Bsr::from_csr(&a, 2, 2);
        let x = [0.25, 1.0, f64::INFINITY, f64::NAN, f64::NAN, f64::NEG_INFINITY, 2.0, -3.0];
        let mut expected = [0.0; 4];
        let mut actual = [f64::NAN; 4];
        a.spmm(&x, 2, &mut expected);
        b.spmm(&x, 2, &mut actual);
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        assert_eq!(actual[0].to_bits(), (-0.0f64).to_bits());
        actual.fill(f64::NAN);
        b.spmm_blocked(&x, 2, &mut actual);
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));

        let zeros = Csr::from_parts(2, 2, vec![0, 1, 1], vec![0], vec![0.0]);
        Bsr::from_csr(&zeros, 2, 2).spmm(&[f64::INFINITY, 2.0, f64::NAN, f64::NAN], 2, &mut actual);
        assert!(actual[0].is_nan());
        assert_eq!(actual[1..].iter().map(|v| v.to_bits()).collect::<Vec<_>>(), vec![0; 3]);
        Bsr::from_csr(&zeros, 2, 2).spmm_blocked(
            &[f64::INFINITY, 2.0, f64::NAN, f64::NAN], 2, &mut actual,
        );
        assert!(actual[0].is_nan());
        assert!(actual[1..].iter().all(|v| v.to_bits() == 0));
    }

    #[test]
    fn g0_spmm_handles_zero_axes() {
        let empty_columns = Csr::from_parts(3, 0, vec![0; 4], vec![], vec![]);
        let mut y = [f64::NAN; 9];
        Bsr::from_csr(&empty_columns, 3, 1).spmm(&[], 3, &mut y);
        assert!(y.iter().all(|v| v.to_bits() == 0));
        y.fill(f64::NAN);
        Bsr::from_csr(&empty_columns, 3, 1).spmm_blocked(&[], 3, &mut y);
        assert!(y.iter().all(|v| v.to_bits() == 0));
        let empty_rows = Csr::from_parts(0, 6, vec![0], vec![], vec![]);
        Bsr::from_csr(&empty_rows, 2, 3).spmm(&[f64::NAN; 18], 3, &mut []);
        Bsr::from_csr(&empty_rows, 2, 3).spmm_blocked(&[f64::NAN; 18], 3, &mut []);
        let empty = Bsr::from_csr(&Csr::identity(0), 1, 1);
        empty.spmm(&[], usize::MAX, &mut []);
        empty.spmm_blocked(&[], usize::MAX, &mut []);
        empty.spmm_blocked_body(&[], usize::MAX, &mut []);
    }

    #[test]
    fn g0_spmm_shape_refusals_leave_output_unchanged() {
        let b = Bsr::from_csr(&Csr::identity(2), 1, 1);
        for (x_len, nrhs, y_len) in [(3, 2, 4), (4, 2, 3), (0, usize::MAX, 4), (0, 0, 1)] {
            let x = vec![1.0; x_len];
            let mut y = vec![42.0; y_len];
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                b.spmm(&x, nrhs, &mut y);
            }));
            assert!(result.is_err());
            assert_eq!(y, vec![42.0; y_len]);
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                b.spmm_blocked(&x, nrhs, &mut y);
            }));
            assert!(result.is_err());
            assert_eq!(y, vec![42.0; y_len]);
        }
        let zero_columns = Csr::from_parts(2, 0, vec![0; 3], vec![], vec![]);
        let b = Bsr::from_csr(&zero_columns, 1, 1);
        let mut y = [42.0];
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            b.spmm(&[], usize::MAX, &mut y);
        }));
        assert!(result.is_err(), "output dimension overflow must be rejected");
        assert_eq!(y, [42.0]);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            b.spmm_blocked(&[], usize::MAX, &mut y);
        }));
        assert!(result.is_err(), "blocked output dimension overflow must be rejected");
        assert_eq!(y, [42.0]);
    }

    #[test]
    fn rejects_indivisible_dimensions() {
        let a = laplacian_2d(3); // 9x9
        let r = std::panic::catch_unwind(|| Bsr::from_csr(&a, 2, 2));
        assert!(r.is_err(), "9x9 with 2x2 blocks must be refused");
    }
}

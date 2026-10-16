//! Block CSR — FEM vector-valued unknowns make dense r×c micro-blocks;
//! BSR stores them contiguously so the future SIMD lane gets
//! batched-small-dense-adjacent work (plan §6.2).
//!
//! Bitwise note: BSR SpMV visits each scalar row's entries in ascending
//! global column order (block columns ascend, columns within a block
//! ascend), and fill zeros provably cannot perturb bits (a fused
//! `0·x + acc` is exactly `acc`, and `acc` can never become −0.0 from a
//! +0.0 start under round-to-nearest) — so BSR SpMV is bit-identical to
//! CSR SpMV. Tested, not just argued.

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
        let nbrows = a.nrows() / br;
        let mut brow_ptr = vec![0usize; nbrows + 1];
        let mut bcol_idx: Vec<usize> = Vec::new();
        let mut blocks: Vec<f64> = Vec::new();
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
            blocks.resize(base_block + present.len() * br * bc, 0.0);
            for i in 0..br {
                let (cols, vals) = a.row(rb * br + i);
                for (&c, &v) in cols.iter().zip(vals) {
                    let cb = c / bc;
                    let slot = present.binary_search(&cb).expect("built above");
                    blocks[base_block + slot * br * bc + i * bc + (c % bc)] = v;
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
        }
    }

    /// Expand back to CSR. Fill zeros introduced by blocking are DROPPED
    /// (exact-0.0 test), so `from_csr(a).to_csr()` reproduces `a` whenever
    /// `a` stores no explicit zeros; the semantically lossless claim
    /// (identical dense expansion) holds unconditionally and is tested.
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
                        let v = self.blocks[slot * self.br * self.bc + i * self.bc + j];
                        if v != 0.0 {
                            col_idx.push(cb * self.bc + j);
                            vals.push(v);
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
    /// fused mul_add; fill zeros provably inert — see module docs).
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
                    let block_row =
                        &self.blocks[slot * self.br * self.bc + i * self.bc..][..self.bc];
                    for (j, &v) in block_row.iter().enumerate() {
                        acc = v.mul_add(x[cb * self.bc + j], acc);
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
    fn rejects_indivisible_dimensions() {
        let a = laplacian_2d(3); // 9x9
        let r = std::panic::catch_unwind(|| Bsr::from_csr(&a, 2, 2));
        assert!(r.is_err(), "9x9 with 2x2 blocks must be refused");
    }
}

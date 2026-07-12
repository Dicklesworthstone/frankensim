//! Pattern algebra on canonical CSR: transpose, symmetrization, and
//! Gustavson SpGEMM (the kernel AMG's Galerkin triple product R·A·P is built
//! from later). All deterministic by construction: iteration order is a pure
//! function of the (canonical) structure.

use crate::Csr;

/// Aᵀ via counting sort — O(nnz + nrows + ncols), values preserved bitwise,
/// output canonical (within a column of A, rows are visited ascending, which
/// becomes ascending columns of Aᵀ).
#[must_use]
pub fn transpose(a: &Csr) -> Csr {
    let (nr, nc) = (a.nrows(), a.ncols());
    let mut counts = vec![0usize; nc + 1];
    for r in 0..nr {
        for &c in a.row(r).0 {
            counts[c + 1] += 1;
        }
    }
    for c in 0..nc {
        counts[c + 1] += counts[c];
    }
    let row_ptr = counts.clone();
    let nnz = a.nnz();
    let mut col_idx = vec![0usize; nnz];
    let mut vals = vec![0.0f64; nnz];
    let mut cursor = counts;
    for r in 0..nr {
        let (cols, values) = a.row(r);
        for (&c, &v) in cols.iter().zip(values) {
            col_idx[cursor[c]] = r;
            vals[cursor[c]] = v;
            cursor[c] += 1;
        }
    }
    Csr::from_parts(nc, nr, row_ptr, col_idx, vals)
}

/// S = (A + Aᵀ)/2 on the pattern UNION. Requires square A. The per-entry
/// value is `0.5 * (a + at)`: IEEE addition is commutative, so
/// S[i][j] and S[j][i] compute bit-identical values — symmetry holds
/// BITWISE (tested), which downstream symmetric-solver checks rely on.
#[must_use]
pub fn symmetrize(a: &Csr) -> Csr {
    assert_eq!(a.nrows(), a.ncols(), "symmetrize requires a square matrix");
    let at = transpose(a);
    let n = a.nrows();
    let mut row_ptr = vec![0usize; n + 1];
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for r in 0..n {
        let (ca, va) = a.row(r);
        let (cb, vb) = at.row(r);
        // Two-pointer merge over sorted column lists; a missing side
        // contributes exactly 0.0. `midpoint` is commutative in IEEE
        // arithmetic, so the bitwise-symmetry argument is preserved.
        let (mut i, mut j) = (0usize, 0usize);
        while i < ca.len() || j < cb.len() {
            let next_a = ca.get(i).copied().unwrap_or(usize::MAX);
            let next_b = cb.get(j).copied().unwrap_or(usize::MAX);
            let c = next_a.min(next_b);
            let mut x = 0.0f64;
            let mut y = 0.0f64;
            if next_a == c {
                x = va[i];
                i += 1;
            }
            if next_b == c {
                y = vb[j];
                j += 1;
            }
            col_idx.push(c);
            vals.push(f64::midpoint(x, y));
        }
        row_ptr[r + 1] = col_idx.len();
    }
    Csr::from_parts(n, n, row_ptr, col_idx, vals)
}

/// C = A·B, Gustavson's algorithm with a dense sparse-accumulator (SPA) per
/// row. Deterministic: contributions to C[i][j] accumulate in ascending-k
/// order (A's column order), and the output pattern is sorted ascending.
/// Complexity O(Σᵢ Σ_{k∈row i} nnz(B row k)) — the right shape for the AMG
/// triple product.
#[must_use]
pub fn spgemm(a: &Csr, b: &Csr) -> Csr {
    assert_eq!(
        a.ncols(),
        b.nrows(),
        "spgemm dimension mismatch: {}x{} times {}x{}",
        a.nrows(),
        a.ncols(),
        b.nrows(),
        b.ncols()
    );
    crate::fma::spgemm_dispatch(a, b)
}

/// C = A·B with a SPARSE accumulator (bead wsbf item 6): a per-row
/// BTreeMap replaces the dense O(ncols) SPA scratch — the win for VERY
/// WIDE B (a 10⁶-column product allocates nothing proportional to
/// ncols). DETERMINISM: per (row, col) the mul_add chain accumulates
/// in exactly the dense SPA's traversal order, and BTreeMap iteration
/// is ascending-column — bitwise equality with [`spgemm`] is GATED.
/// (BTree, not hash: deterministic order needs no post-sort and no
/// seeded hasher — the "hash-SPA" of the bead, honestly named.)
#[must_use]
pub fn spgemm_sparse_spa(a: &Csr, b: &Csr) -> Csr {
    assert_eq!(
        a.ncols(),
        b.nrows(),
        "spgemm dimension mismatch: {}x{} times {}x{}",
        a.nrows(),
        a.ncols(),
        b.nrows(),
        b.ncols()
    );
    crate::fma::spgemm_sparse_spa_dispatch(a, b)
}

impl crate::Csr {
    /// The [`spgemm`] loop body (dense SPA), extracted so the x86
    /// FMA-codegen capsule can recompile it under `target_feature`
    /// (bead nabk). MUST stay `inline(always)`: a non-inlined call
    /// keeps baseline codegen and the per-element libm `fma()` call.
    /// Accumulation order and output pattern untouched: pure codegen.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spgemm_body(&self, b: &Csr) -> Csr {
        let (m, n) = (self.nrows(), b.ncols());
        let mut spa = vec![0.0f64; n];
        let mut occupied = vec![false; n];
        let mut touched: Vec<usize> = Vec::new();
        let mut row_ptr = vec![0usize; m + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        for i in 0..m {
            let (acols, avals) = self.row(i);
            for (&k, &aik) in acols.iter().zip(avals) {
                let (bcols, bvals) = b.row(k);
                for (&j, &bkj) in bcols.iter().zip(bvals) {
                    spa[j] = aik.mul_add(bkj, spa[j]);
                    if !occupied[j] {
                        occupied[j] = true;
                        touched.push(j);
                    }
                }
            }
            touched.sort_unstable();
            for &j in &touched {
                col_idx.push(j);
                vals.push(spa[j]);
                spa[j] = 0.0;
                occupied[j] = false;
            }
            row_ptr[i + 1] = col_idx.len();
            touched.clear();
        }
        Csr::from_parts(m, n, row_ptr, col_idx, vals)
    }

    /// The [`spgemm_sparse_spa`] loop body (see `spgemm_body`'s
    /// extraction note).
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spgemm_sparse_spa_body(&self, b: &Csr) -> Csr {
        let (m, n) = (self.nrows(), b.ncols());
        let mut row_ptr = vec![0usize; m + 1];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        let mut spa: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
        for i in 0..m {
            let (acols, avals) = self.row(i);
            for (&k, &aik) in acols.iter().zip(avals) {
                let (bcols, bvals) = b.row(k);
                for (&j, &bkj) in bcols.iter().zip(bvals) {
                    let slot = spa.entry(j).or_insert(0.0);
                    *slot = aik.mul_add(bkj, *slot);
                }
            }
            for (&j, &v) in &spa {
                col_idx.push(j);
                vals.push(v);
            }
            row_ptr[i + 1] = col_idx.len();
            spa.clear();
        }
        Csr::from_parts(m, n, row_ptr, col_idx, vals)
    }

    /// Blocked SpMM (bead wsbf segment 2): Y = A · B for row-major
    /// dense B (ncols × nrhs) into row-major Y (nrows × nrhs). The
    /// rhs columns are processed in blocks so each A traversal feeds
    /// `NB` accumulators (block-Krylov/POD shapes reuse the matrix
    /// stream instead of re-reading it per vector). DETERMINISM:
    /// column j of Y accumulates in exactly [`crate::Csr::spmv`]'s
    /// k-ascending order — bitwise equality with per-column SpMV is
    /// GATED in conformance.
    pub fn spmm_blocked(&self, b: &[f64], nrhs: usize, y: &mut [f64]) {
        assert_eq!(b.len(), self.ncols() * nrhs, "spmm: B is ncols x nrhs");
        assert_eq!(y.len(), self.nrows() * nrhs, "spmm: Y is nrows x nrhs");
        crate::fma::spmm_blocked_dispatch(self, b, nrhs, y);
    }

    /// The [`Csr::spmm_blocked`] loop body (see `spgemm_body`'s
    /// extraction note).
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_blocked_body(&self, b: &[f64], nrhs: usize, y: &mut [f64]) {
        const NB: usize = 8;
        let mut acc = [0.0f64; NB];
        for j0 in (0..nrhs).step_by(NB) {
            let jw = NB.min(nrhs - j0);
            for r in 0..self.nrows() {
                let (cols, vals) = self.row(r);
                acc[..jw].fill(0.0);
                for (&c, &v) in cols.iter().zip(vals) {
                    let brow = &b[c * nrhs + j0..c * nrhs + j0 + jw];
                    for (a, &bv) in acc[..jw].iter_mut().zip(brow) {
                        *a = v.mul_add(bv, *a);
                    }
                }
                y[r * nrhs + j0..r * nrhs + j0 + jw].copy_from_slice(&acc[..jw]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Coo;
    use crate::tests::{laplacian_2d, random_coo};

    fn bitwise_eq(a: &Csr, b: &Csr) -> bool {
        if a.nrows() != b.nrows() || a.ncols() != b.ncols() || a.nnz() != b.nnz() {
            return false;
        }
        (0..a.nrows()).all(|r| {
            let (c1, v1) = a.row(r);
            let (c2, v2) = b.row(r);
            c1 == c2 && v1.iter().zip(v2).all(|(x, y)| x.to_bits() == y.to_bits())
        })
    }

    #[test]
    fn transpose_involution_is_bitwise() {
        for n in [1usize, 13, 50] {
            let a = random_coo(n, n + 3, 5, 0x7A + n as u64).assemble();
            let att = transpose(&transpose(&a));
            assert!(bitwise_eq(&a, &att), "(A^T)^T != A at n={n}");
        }
        // Rectangular shape checks too.
        let a = random_coo(20, 7, 3, 9).assemble();
        let at = transpose(&a);
        assert_eq!(at.nrows(), 7);
        assert_eq!(at.ncols(), 20);
        for r in 0..20 {
            let (cols, vals) = a.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                assert_eq!(at.get(c, r).to_bits(), v.to_bits());
            }
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"transpose\",\"verdict\":\"pass\",\"detail\":\"involution bitwise, rectangular entries verified\"}}"
        );
    }

    #[test]
    fn symmetrize_is_bitwise_symmetric() {
        let a = random_coo(40, 40, 6, 0x51).assemble();
        let s = symmetrize(&a);
        let st = transpose(&s);
        assert!(
            bitwise_eq(&s, &st),
            "symmetrize output not bitwise symmetric"
        );
        // Already-symmetric input is a fixed point (Laplacian).
        let lap = laplacian_2d(6);
        let sl = symmetrize(&lap);
        assert!(
            bitwise_eq(&lap, &sl),
            "symmetric input must be a fixed point"
        );
    }

    #[test]
    fn spgemm_identity_is_bitwise_and_matches_dense() {
        let a = random_coo(25, 30, 5, 0x6E).assemble();
        // A·I == A and I·A == A, bitwise (single exact contribution per slot).
        let right = spgemm(&a, &Csr::identity(30));
        let left = spgemm(&Csr::identity(25), &a);
        assert!(bitwise_eq(&a, &right), "A*I != A");
        assert!(bitwise_eq(&a, &left), "I*A != A");
        // General product vs dense oracle.
        let b = random_coo(30, 18, 4, 0xBEE).assemble();
        let c = spgemm(&a, &b);
        let (da, db, dc) = (a.to_dense(), b.to_dense(), c.to_dense());
        for i in 0..25 {
            for j in 0..18 {
                let want: f64 = (0..30).map(|k| da[i * 30 + k] * db[k * 18 + j]).sum();
                let got = dc[i * 18 + j];
                let scale = want.abs().max(1.0);
                assert!(
                    (got - want).abs() / scale < 1e-13,
                    "spgemm ({i},{j}): {got} vs dense {want}"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-sparse\",\"case\":\"spgemm\",\"verdict\":\"pass\",\"detail\":\"identity bitwise + 25x30x18 dense oracle\"}}"
        );
    }

    #[test]
    fn spgemm_laplacian_square_pattern_sanity() {
        // A² of the 5-point Laplacian has the 9-point-plus pattern: bandwidth
        // doubles; diagonal entry = 4²+ Σ neighbors² = 20 in the interior.
        let lap = laplacian_2d(8);
        let sq = spgemm(&lap, &lap);
        let interior = 3 * 8 + 3; // (3,3) — fully interior node
        assert_eq!(sq.get(interior, interior).to_bits(), 20.0f64.to_bits());
        assert!(sq.nnz() > lap.nnz(), "A^2 must widen the pattern");
    }

    #[test]
    fn spgemm_dimension_mismatch_is_refused() {
        let a = Coo::new(3, 4).assemble();
        let b = Coo::new(5, 2).assemble();
        let r = std::panic::catch_unwind(|| spgemm(&a, &b));
        assert!(r.is_err(), "3x4 times 5x2 must be refused");
    }
}

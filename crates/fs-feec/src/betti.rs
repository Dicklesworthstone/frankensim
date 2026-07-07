//! Exact integer rank of incidence operators (fraction-free Bareiss
//! elimination in i128 — no floating point, no tolerance knobs) and
//! the rank–nullity Betti bookkeeping b_k = dim C_k − rank d_k −
//! rank d_{k−1}. Fixture-scale only: entries grow during Bareiss, and
//! the point is CERTIFYING small complexes, not persistent homology
//! at scale (that is fs-topo's territory).

use fs_rep_mesh::TetComplex;

/// Exact rank over ℚ of an integer incidence operator given as rows
/// of (column, ±1) pairs.
///
/// # Panics
/// If an intermediate exceeds i128 (fixture far too large for this
/// certifier — use a bigger tool, don't loosen this one).
#[must_use]
pub fn integer_rank(rows: &[Vec<(usize, i8)>], cols: usize) -> usize {
    let nr = rows.len();
    if nr == 0 || cols == 0 {
        return 0;
    }
    let mut m = vec![0i128; nr * cols];
    for (r, row) in rows.iter().enumerate() {
        for &(c, s) in row {
            m[r * cols + c] = i128::from(s);
        }
    }
    // Fraction-free Gaussian elimination (Bareiss): division by the
    // previous pivot is exact.
    let mut rank = 0usize;
    let mut prev_pivot = 1i128;
    let mut row = 0usize;
    for col in 0..cols {
        // Deterministic pivot: lowest row index with a nonzero entry.
        let Some(p) = (row..nr).find(|&r| m[r * cols + col] != 0) else {
            continue;
        };
        if p != row {
            for c in 0..cols {
                m.swap(row * cols + c, p * cols + c);
            }
        }
        let pivot = m[row * cols + col];
        for r in (row + 1)..nr {
            let head = m[r * cols + col];
            for c in 0..cols {
                let a = m[r * cols + c]
                    .checked_mul(pivot)
                    .expect("Bareiss overflow: fixture too large for i128 certifier");
                let b = m[row * cols + c]
                    .checked_mul(head)
                    .expect("Bareiss overflow: fixture too large for i128 certifier");
                m[r * cols + c] = (a - b) / prev_pivot;
            }
        }
        prev_pivot = pivot;
        row += 1;
        rank += 1;
        if row == nr {
            break;
        }
    }
    rank
}

/// Betti numbers (b₀, b₁, b₂, b₃) of the complex by rank–nullity over
/// the exact integer incidence operators. For a solid ball fixture:
/// (1, 0, 0, 0).
#[must_use]
pub fn betti_numbers(complex: &TetComplex) -> [usize; 4] {
    let d0 = complex.d0();
    let d1 = complex.d1();
    let d2 = complex.d2();
    let (nv, ne, nf, nt) = (
        complex.vertex_count,
        complex.edges.len(),
        complex.faces.len(),
        complex.tets.len(),
    );
    let r0 = integer_rank(&d0.rows, d0.cols);
    let r1 = integer_rank(&d1.rows, d1.cols);
    let r2 = integer_rank(&d2.rows, d2.cols);
    [nv - r0, ne - r0 - r1, nf - r1 - r2, nt - r2]
}

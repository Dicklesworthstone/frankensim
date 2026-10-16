//! NEON complex-transpose capsule (bead 27d3): the 8×8-tiled
//! out-of-place transpose of an n₁×n₁ interleaved-complex matrix, the
//! pass that dominates fs-fft's six-step decomposition.
//!
//! One interleaved complex f64 IS one 128-bit q-register, so the tile
//! move needs no shuffles at all: `vld1q_f64` from the strided source
//! column, `vst1q_f64` to the sequential destination row. The win over
//! the scalar twin is eliminating the per-element bounds checks and
//! moving 16 bytes per instruction; the access PATTERN (sequential
//! writes, 8-line-bounded strided reads) is the twin's exactly.
//!
//! BITWISE contract: pure exact element moves — no floating-point
//! arithmetic anywhere, so tier equivalence is equality of moves, gated
//! bitwise against the scalar twin in the lib battery.

use core::arch::aarch64::{vld1q_f64, vst1q_f64};

/// Tile edge: 8 complex = 128 B per row segment (one Apple cache line).
const TILE: usize = 8;

/// `dst[i·n1 + j] = src[j·n1 + i]` over n₁×n₁ complex elements stored
/// interleaved (`[re, im]` per element; slice length `2·n1²`).
///
/// # Panics
/// Structured panics when either slice length is not `2·n1²`.
pub fn trn1c64(src: &[f64], dst: &mut [f64], n1: usize) {
    let need = crate::checked_trn1c64_len(n1)
        .unwrap_or_else(|| panic!("trn1c64 extent overflow (programmer error)"));
    assert_eq!(src.len(), need, "trn1c64 src length (programmer error)");
    assert_eq!(dst.len(), need, "trn1c64 dst length (programmer error)");
    let sp = src.as_ptr();
    let dp = dst.as_mut_ptr();
    let mut bi = 0;
    while bi < n1 {
        let i_end = (bi + TILE).min(n1);
        let mut bj = 0;
        while bj < n1 {
            let j_end = (bj + TILE).min(n1);
            for i in bi..i_end {
                // SAFETY: i < n1 and j < n1 throughout, so both element
                // indices j·n1+i and i·n1+j are < n1², and the byte
                // offsets 2·(idx)+1 are inside the asserted 2·n1² slice
                // lengths. `src` and `dst` are distinct slices (shared +
                // exclusive borrows), so reads and writes cannot alias.
                unsafe {
                    for j in bj..j_end {
                        let v = vld1q_f64(sp.add(2 * (j * n1 + i)));
                        vst1q_f64(dp.add(2 * (i * n1 + j)), v);
                    }
                }
            }
            bj += TILE;
        }
        bi += TILE;
    }
}

/// Gather columns `g..g+8` of an n₁×n₁ interleaved-complex matrix into
/// a dense 8×n₁ buffer: `bufs[c·n1 + i] = src[i·n1 + g + c]` (complex
/// indices; slices are f64 with lengths `2·n1²` and `2·8·n1`). One
/// 128-byte line of the source feeds all eight buffer streams, so the
/// strided side still touches every line exactly once. Pure exact
/// moves — BITWISE across tiers by construction (bead 27d3, the fused
/// six-step's inner loop).
///
/// # Panics
/// Structured panics on length/geometry mismatches.
pub fn gath8c64(src: &[f64], bufs: &mut [f64], n1: usize, g: usize) {
    let need = crate::checked_trn1c64_len(n1)
        .unwrap_or_else(|| panic!("gath8c64 extent overflow (programmer error)"));
    assert_eq!(src.len(), need, "gath8c64 src length (programmer error)");
    assert_eq!(
        bufs.len(),
        16 * n1,
        "gath8c64 bufs length (programmer error)"
    );
    let group_end = g
        .checked_add(TILE)
        .unwrap_or_else(|| panic!("gath8c64 column group out of range (programmer error)"));
    assert!(
        group_end <= n1,
        "gath8c64 column group out of range (programmer error)"
    );
    let sp = src.as_ptr();
    let bp = bufs.as_mut_ptr();
    for i in 0..n1 {
        // SAFETY: i < n1 and g + c < n1 (asserted), so the source
        // element index i·n1 + g + c is < n1² and its two f64s are
        // inside the asserted 2·n1² slice; the buffer index c·n1 + i is
        // < 8·n1 inside the asserted 16·n1 f64 slice. src and bufs are
        // distinct slices (shared + exclusive borrows) — no aliasing.
        unsafe {
            let line = sp.add(2 * (i * n1 + g));
            for c in 0..8 {
                let v = vld1q_f64(line.add(2 * c));
                vst1q_f64(bp.add(2 * (c * n1 + i)), v);
            }
        }
    }
}

/// Scatter a dense 8×n₁ buffer into columns `g..g+8` of an n₁×n₁
/// interleaved-complex matrix: `dst[k·n1 + g + c] = bufs[c·n1 + k]` —
/// the exact inverse of [`gath8c64`], and also the fused six-step's
/// row→column output move (eight contiguous rows ARE a dense buffer).
/// Pure exact moves — BITWISE across tiers by construction.
///
/// # Panics
/// Structured panics on length/geometry mismatches.
pub fn scat8c64(bufs: &[f64], dst: &mut [f64], n1: usize, g: usize) {
    let need = crate::checked_trn1c64_len(n1)
        .unwrap_or_else(|| panic!("scat8c64 extent overflow (programmer error)"));
    assert_eq!(dst.len(), need, "scat8c64 dst length (programmer error)");
    assert_eq!(
        bufs.len(),
        16 * n1,
        "scat8c64 bufs length (programmer error)"
    );
    let group_end = g
        .checked_add(TILE)
        .unwrap_or_else(|| panic!("scat8c64 column group out of range (programmer error)"));
    assert!(
        group_end <= n1,
        "scat8c64 column group out of range (programmer error)"
    );
    let bp = bufs.as_ptr();
    let dp = dst.as_mut_ptr();
    for k in 0..n1 {
        // SAFETY: mirror of gath8c64 — every index is bounded by the
        // asserts above; bufs and dst cannot alias (shared + exclusive
        // borrows).
        unsafe {
            let line = dp.add(2 * (k * n1 + g));
            for c in 0..8 {
                let v = vld1q_f64(bp.add(2 * (c * n1 + k)));
                vst1q_f64(line.add(2 * c), v);
            }
        }
    }
}

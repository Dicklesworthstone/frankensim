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
    let need = 2 * n1 * n1;
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

//! NEON f32 batched-GEMM capsule (bead 9ekv scope e): the PACKED
//! l-contiguous 4×4 tile kernel at f32 — FOUR lanes per vector
//! register instead of two, the single biggest per-register win
//! available on NEON-128 for the batched-lane layout (see the bead's
//! ceiling analysis: batched lanes cannot broadcast, so throughput
//! scales with lane width). Registered in unsafe-capsules.json;
//! SAFETY.md beside this file. Same layout contract as the f64
//! `btile4x4p_f64`: A packed i-major, B packed j-major, both walks
//! stride `mb` per l; per element identical to the scalar twin —
//! BITWISE per tier.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use core::arch::aarch64::{vdupq_n_f32, vfmaq_f32, vld1q_f32, vst1q_f32};

/// PACKED f32 batched-GEMM 4×4 tile (see module docs). Lane counts not
/// divisible by 4 take the scalar twin whole.
#[allow(clippy::too_many_arguments)] // packed-layout bundle (see fs-la::batched)
pub fn btile4x4pf32(
    a: &[f32],
    b: &[f32],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f32],
) {
    assert!(
        k >= 1
            && (i0 + 3) * k * mb + k * mb <= a.len()
            && (j0 + 3) * k * mb + k * mb <= b.len()
            && dst.len() >= 16 * mb,
        "btile4x4pf32 packed bounds (programmer error)"
    );
    if mb % 4 != 0 {
        crate::scalar::btile4x4pf32(a, b, i0, j0, k, mb, dst);
        return;
    }
    let quads = mb / 4;
    // SAFETY: every pointer is base(t) + l·mb + 4·q with base(t) =
    // ((i0+t)·k)·mb (a) or ((j0+t)·k)·mb (b); the maximal dereferenced
    // offset over t ≤ 3, l ≤ k−1, 4q ≤ mb−4 is inside the extents
    // asserted above. Every access is exactly 4 f32; f32 has no invalid
    // bit patterns; unaligned access is permitted. The per-quad rewind
    // (−k·mb + 4) never leaves the borrowed allocations.
    unsafe {
        let (mut a0p, mut a1p, mut a2p, mut a3p) = (
            a.as_ptr().add(i0 * k * mb),
            a.as_ptr().add((i0 + 1) * k * mb),
            a.as_ptr().add((i0 + 2) * k * mb),
            a.as_ptr().add((i0 + 3) * k * mb),
        );
        let (mut b0p, mut b1p, mut b2p, mut b3p) = (
            b.as_ptr().add(j0 * k * mb),
            b.as_ptr().add((j0 + 1) * k * mb),
            b.as_ptr().add((j0 + 2) * k * mb),
            b.as_ptr().add((j0 + 3) * k * mb),
        );
        let op = dst.as_mut_ptr();
        for q in 0..quads {
            let mut acc = [vdupq_n_f32(0.0); 16];
            for _l in 0..k {
                let a0 = vld1q_f32(a0p);
                let a1 = vld1q_f32(a1p);
                let a2 = vld1q_f32(a2p);
                let a3 = vld1q_f32(a3p);
                let b0 = vld1q_f32(b0p);
                let b1 = vld1q_f32(b1p);
                let b2 = vld1q_f32(b2p);
                let b3 = vld1q_f32(b3p);
                acc[0] = vfmaq_f32(acc[0], a0, b0);
                acc[1] = vfmaq_f32(acc[1], a0, b1);
                acc[2] = vfmaq_f32(acc[2], a0, b2);
                acc[3] = vfmaq_f32(acc[3], a0, b3);
                acc[4] = vfmaq_f32(acc[4], a1, b0);
                acc[5] = vfmaq_f32(acc[5], a1, b1);
                acc[6] = vfmaq_f32(acc[6], a1, b2);
                acc[7] = vfmaq_f32(acc[7], a1, b3);
                acc[8] = vfmaq_f32(acc[8], a2, b0);
                acc[9] = vfmaq_f32(acc[9], a2, b1);
                acc[10] = vfmaq_f32(acc[10], a2, b2);
                acc[11] = vfmaq_f32(acc[11], a2, b3);
                acc[12] = vfmaq_f32(acc[12], a3, b0);
                acc[13] = vfmaq_f32(acc[13], a3, b1);
                acc[14] = vfmaq_f32(acc[14], a3, b2);
                acc[15] = vfmaq_f32(acc[15], a3, b3);
                a0p = a0p.add(mb);
                a1p = a1p.add(mb);
                a2p = a2p.add(mb);
                a3p = a3p.add(mb);
                b0p = b0p.add(mb);
                b1p = b1p.add(mb);
                b2p = b2p.add(mb);
                b3p = b3p.add(mb);
            }
            let dp = op.add(4 * q);
            for (t, &v) in acc.iter().enumerate() {
                vst1q_f32(dp.add(t * mb), v);
            }
            let rewind = k * mb - 4;
            a0p = a0p.sub(rewind);
            a1p = a1p.sub(rewind);
            a2p = a2p.sub(rewind);
            a3p = a3p.sub(rewind);
            b0p = b0p.sub(rewind);
            b1p = b1p.sub(rewind);
            b2p = b2p.sub(rewind);
            b3p = b3p.sub(rewind);
        }
    }
}

//! NEON GEMM microkernel capsule (split from neon/mod.rs under the
//! 300-line capsule cap, bead 8nfp): the packed-panel 8x4 kernel
//! (bead xdgf) and the batched plane-SoA 4x4 entry-tile kernel (bead
//! 9ekv). Registered in unsafe-capsules.json; SAFETY.md beside this
//! file. Bitwise contracts and the tier battery are unchanged — this
//! is a pure file move.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use core::arch::aarch64::{
    float64x2_t, vdupq_n_f64, vfmaq_f64, vfmaq_laneq_f64, vld1q_f64, vst1q_f64,
};

/// The 8×4 f64 GEMM register microkernel: 16 `float64x2` accumulators
/// (8 rows × 2 column pairs) resident across the whole k loop, k
/// ascending, `vfmaq_laneq` broadcasting each packed A lane. Per
/// element this is exactly `acc[r][s] = fma(a[r], b[s], acc[r][s])` in
/// the scalar twin's order — BITWISE-identical, so fs-la's GEMM golden
/// is tier-invariant.
pub fn mk8x4_f64(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    assert!(
        a_panel.len() >= kc * 8 && b_panel.len() >= kc * 4,
        "mk8x4 panel length mismatch (programmer error)"
    );
    // SAFETY: every vld1q/vst1q reads or writes exactly 2 f64 at offsets
    // kept in bounds by the assert above (a: kk·8+6+2 ≤ kc·8; b: kk·4+2+2
    // ≤ kc·4) and by `acc`'s [[f64; 4]; 8] type. f64 has no invalid bit
    // patterns; vld1q/vst1q tolerate unaligned addresses.
    unsafe {
        let ap = a_panel.as_ptr();
        let bp = b_panel.as_ptr();
        let mut va: [[float64x2_t; 2]; 8] = [[vdupq_n_f64(0.0); 2]; 8];
        for (r, v) in va.iter_mut().enumerate() {
            v[0] = vld1q_f64(acc[r].as_ptr());
            v[1] = vld1q_f64(acc[r].as_ptr().add(2));
        }
        for kk in 0..kc {
            let b0 = vld1q_f64(bp.add(kk * 4));
            let b1 = vld1q_f64(bp.add(kk * 4 + 2));
            let a01 = vld1q_f64(ap.add(kk * 8));
            let a23 = vld1q_f64(ap.add(kk * 8 + 2));
            let a45 = vld1q_f64(ap.add(kk * 8 + 4));
            let a67 = vld1q_f64(ap.add(kk * 8 + 6));
            va[0][0] = vfmaq_laneq_f64::<0>(va[0][0], b0, a01);
            va[0][1] = vfmaq_laneq_f64::<0>(va[0][1], b1, a01);
            va[1][0] = vfmaq_laneq_f64::<1>(va[1][0], b0, a01);
            va[1][1] = vfmaq_laneq_f64::<1>(va[1][1], b1, a01);
            va[2][0] = vfmaq_laneq_f64::<0>(va[2][0], b0, a23);
            va[2][1] = vfmaq_laneq_f64::<0>(va[2][1], b1, a23);
            va[3][0] = vfmaq_laneq_f64::<1>(va[3][0], b0, a23);
            va[3][1] = vfmaq_laneq_f64::<1>(va[3][1], b1, a23);
            va[4][0] = vfmaq_laneq_f64::<0>(va[4][0], b0, a45);
            va[4][1] = vfmaq_laneq_f64::<0>(va[4][1], b1, a45);
            va[5][0] = vfmaq_laneq_f64::<1>(va[5][0], b0, a45);
            va[5][1] = vfmaq_laneq_f64::<1>(va[5][1], b1, a45);
            va[6][0] = vfmaq_laneq_f64::<0>(va[6][0], b0, a67);
            va[6][1] = vfmaq_laneq_f64::<0>(va[6][1], b1, a67);
            va[7][0] = vfmaq_laneq_f64::<1>(va[7][0], b0, a67);
            va[7][1] = vfmaq_laneq_f64::<1>(va[7][1], b1, a67);
        }
        for (r, v) in va.iter().enumerate() {
            vst1q_f64(acc[r].as_mut_ptr(), v[0]);
            vst1q_f64(acc[r].as_mut_ptr().add(2), v[1]);
        }
    }
}

/// Batched-GEMM 4×4 entry-tile microkernel (bead 9ekv): 16 resident
/// float64x2 accumulators (one per tile entry, 2 batch lanes), plane
/// pointers advancing by `stride` per l — 8 loads : 16 FMAs, so the
/// kernel is FMA-bound instead of drowning in accumulator round-trips.
/// Per element identical to the scalar twin (zero start, l-ascending
/// fused accumulate): BITWISE, and batch lanes are independent
/// matrices. Odd-tail lanes (mb % 2) go through the twin.
#[allow(clippy::too_many_arguments)] // plane-SoA layout bundle (see fs-la::batched)
pub fn btile4x4_f64(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    stride: usize,
    k: usize,
    m0: usize,
    mb: usize,
    dst: &mut [f64],
) {
    assert!(
        k >= 1
            && ((i0 + 3) * k + (k - 1)) * stride + m0 + mb <= a.len()
            && ((k - 1) * k + j0 + 3) * stride + m0 + mb <= b.len()
            && dst.len() >= 16 * mb,
        "btile4x4 plane bounds (programmer error)"
    );
    let pairs = mb / 2;
    // SAFETY: every pointer below is a_base(ti) + l·stride + 2·p (resp.
    // b), whose maximum over ti ≤ 3, l ≤ k−1, 2p ≤ mb−2 is inside the
    // extents asserted above; each vld1q/vst1q touches exactly 2 f64;
    // f64 has no invalid bit patterns; unaligned access is permitted.
    unsafe {
        // Tile bases hoisted out of the pair loop; per pair the eight
        // stream pointers are one add each and REWIND by k·stride (a)
        // / k²·stride (b) after the l walk.
        let mut ab = [core::ptr::null::<f64>(); 4];
        let mut bb = [core::ptr::null::<f64>(); 4];
        for t in 0..4 {
            ab[t] = a.as_ptr().add(((i0 + t) * k) * stride + m0);
            bb[t] = b.as_ptr().add((j0 + t) * stride + m0);
        }
        for p in 0..pairs {
            let mut acc = [vdupq_n_f64(0.0); 16];
            let mut l = 0;
            // l-unroll ×2: same per-element order (l ascending), half
            // the loop control.
            while l + 2 <= k {
                for step in 0..2 {
                    let _ = step;
                    let a0 = vld1q_f64(ab[0]);
                    let a1 = vld1q_f64(ab[1]);
                    let a2 = vld1q_f64(ab[2]);
                    let a3 = vld1q_f64(ab[3]);
                    let b0 = vld1q_f64(bb[0]);
                    let b1 = vld1q_f64(bb[1]);
                    let b2 = vld1q_f64(bb[2]);
                    let b3 = vld1q_f64(bb[3]);
                    acc[0] = vfmaq_f64(acc[0], a0, b0);
                    acc[1] = vfmaq_f64(acc[1], a0, b1);
                    acc[2] = vfmaq_f64(acc[2], a0, b2);
                    acc[3] = vfmaq_f64(acc[3], a0, b3);
                    acc[4] = vfmaq_f64(acc[4], a1, b0);
                    acc[5] = vfmaq_f64(acc[5], a1, b1);
                    acc[6] = vfmaq_f64(acc[6], a1, b2);
                    acc[7] = vfmaq_f64(acc[7], a1, b3);
                    acc[8] = vfmaq_f64(acc[8], a2, b0);
                    acc[9] = vfmaq_f64(acc[9], a2, b1);
                    acc[10] = vfmaq_f64(acc[10], a2, b2);
                    acc[11] = vfmaq_f64(acc[11], a2, b3);
                    acc[12] = vfmaq_f64(acc[12], a3, b0);
                    acc[13] = vfmaq_f64(acc[13], a3, b1);
                    acc[14] = vfmaq_f64(acc[14], a3, b2);
                    acc[15] = vfmaq_f64(acc[15], a3, b3);
                    for t in 0..4 {
                        ab[t] = ab[t].add(stride);
                        bb[t] = bb[t].add(k * stride);
                    }
                }
                l += 2;
            }
            if l < k {
                let a0 = vld1q_f64(ab[0]);
                let a1 = vld1q_f64(ab[1]);
                let a2 = vld1q_f64(ab[2]);
                let a3 = vld1q_f64(ab[3]);
                let b0 = vld1q_f64(bb[0]);
                let b1 = vld1q_f64(bb[1]);
                let b2 = vld1q_f64(bb[2]);
                let b3 = vld1q_f64(bb[3]);
                acc[0] = vfmaq_f64(acc[0], a0, b0);
                acc[1] = vfmaq_f64(acc[1], a0, b1);
                acc[2] = vfmaq_f64(acc[2], a0, b2);
                acc[3] = vfmaq_f64(acc[3], a0, b3);
                acc[4] = vfmaq_f64(acc[4], a1, b0);
                acc[5] = vfmaq_f64(acc[5], a1, b1);
                acc[6] = vfmaq_f64(acc[6], a1, b2);
                acc[7] = vfmaq_f64(acc[7], a1, b3);
                acc[8] = vfmaq_f64(acc[8], a2, b0);
                acc[9] = vfmaq_f64(acc[9], a2, b1);
                acc[10] = vfmaq_f64(acc[10], a2, b2);
                acc[11] = vfmaq_f64(acc[11], a2, b3);
                acc[12] = vfmaq_f64(acc[12], a3, b0);
                acc[13] = vfmaq_f64(acc[13], a3, b1);
                acc[14] = vfmaq_f64(acc[14], a3, b2);
                acc[15] = vfmaq_f64(acc[15], a3, b3);
                for t in 0..4 {
                    ab[t] = ab[t].add(stride);
                    bb[t] = bb[t].add(k * stride);
                }
            }
            let dp = dst.as_mut_ptr().add(2 * p);
            for (t, &v) in acc.iter().enumerate() {
                vst1q_f64(dp.add(t * mb), v);
            }
            // Rewind to l = 0 and advance two batch lanes.
            for t in 0..4 {
                ab[t] = ab[t].sub(k * stride).add(2);
                bb[t] = bb[t].sub(k * k * stride).add(2);
            }
        }
    }
    // Odd batch-lane tail: the scalar twin on the last lane.
    if mb % 2 == 1 {
        let mut tail = vec![0.0f64; 16];
        crate::scalar::btile4x4_f64(a, b, i0, j0, stride, k, m0 + mb - 1, 1, &mut tail);
        for t in 0..16 {
            dst[t * mb + mb - 1] = tail[t];
        }
    }
}

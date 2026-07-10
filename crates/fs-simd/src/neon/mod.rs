//! NEON capsule (aarch64): registered unsafe capsule per unsafe-capsules.json
//! and SAFETY.md in this directory. THE exemplar capsule (unsafe-safety-cases
//! bead): safe façade, <300 lines, scalar-twin equivalence property-tested.
//!
//! Every public function here is SAFE TO CALL: NEON is architecturally
//! guaranteed on aarch64 (no runtime-detection precondition), and all
//! pointer arithmetic derives from `as_chunks::<N>()` fixed-size arrays
//! whose bounds the type system already proved. Tails (the `as_chunks`
//! remainders) are handled by the scalar twin INSIDE each function, so
//! callers never see a partial contract.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use core::arch::aarch64::{
    float64x2_t, float64x2x2_t, vaddq_f64, vaddvq_f64, vdupq_n_f64, vfmaq_f64, vfmaq_laneq_f64,
    vld1q_f64, vld2q_f64, vmulq_f64, vnegq_f64, vst1q_f64, vst2q_f64, vsubq_f64,
};

const LANES: usize = 2; // float64x2_t

/// y[i] = a * x[i] + y[i] (fused, matching the scalar twin's mul_add).
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len(), "axpy length mismatch (programmer error)");
    let (xc, xt) = x.as_chunks::<LANES>();
    let (yc, yt) = y.as_chunks_mut::<LANES>();
    // SAFETY: vld1q/vst1q read/write exactly LANES f64 at addresses of
    // [f64; LANES] arrays produced by as_chunks — inside the allocation,
    // correctly typed. f64 has no invalid bit patterns; vld1q/vst1q do not
    // require alignment.
    unsafe {
        let va = vdupq_n_f64(a);
        for (xk, yk) in xc.iter().zip(yc) {
            let vx = vld1q_f64(xk.as_ptr());
            let vy = vld1q_f64(yk.as_ptr());
            vst1q_f64(yk.as_mut_ptr(), vfmaq_f64(vy, va, vx));
        }
    }
    crate::scalar::axpy(a, xt, yt);
}

/// x[i] *= a.
pub fn scale(a: f64, x: &mut [f64]) {
    let (xc, xt) = x.as_chunks_mut::<LANES>();
    // SAFETY: as in `axpy` — chunk-array pointers, exact LANES extents.
    unsafe {
        let va = vdupq_n_f64(a);
        for xk in xc {
            let vx = vld1q_f64(xk.as_ptr());
            vst1q_f64(xk.as_mut_ptr(), vmulq_f64(vx, va));
        }
    }
    crate::scalar::scale(a, xt);
}

/// out[i] = a[i] * b[i].
pub fn mul_elem(a: &[f64], b: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "mul_elem length mismatch");
    assert_eq!(a.len(), out.len(), "mul_elem length mismatch");
    let (ac, at) = a.as_chunks::<LANES>();
    let (bc, bt) = b.as_chunks::<LANES>();
    let (oc, ot) = out.as_chunks_mut::<LANES>();
    // SAFETY: as in `axpy`.
    unsafe {
        for ((ak, bk), ok) in ac.iter().zip(bc).zip(oc) {
            let va = vld1q_f64(ak.as_ptr());
            let vb = vld1q_f64(bk.as_ptr());
            vst1q_f64(ok.as_mut_ptr(), vmulq_f64(va, vb));
        }
    }
    crate::scalar::mul_elem(at, bt, ot);
}

/// out[i] = a[i] * b[i] + c[i] (fused).
pub fn fma3(a: &[f64], b: &[f64], c: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "fma3 length mismatch");
    assert_eq!(a.len(), c.len(), "fma3 length mismatch");
    assert_eq!(a.len(), out.len(), "fma3 length mismatch");
    let (ac, at) = a.as_chunks::<LANES>();
    let (bc, bt) = b.as_chunks::<LANES>();
    let (cc, ct) = c.as_chunks::<LANES>();
    let (oc, ot) = out.as_chunks_mut::<LANES>();
    // SAFETY: as in `axpy`.
    unsafe {
        for (((ak, bk), ck), ok) in ac.iter().zip(bc).zip(cc).zip(oc) {
            let va = vld1q_f64(ak.as_ptr());
            let vb = vld1q_f64(bk.as_ptr());
            let vc = vld1q_f64(ck.as_ptr());
            vst1q_f64(ok.as_mut_ptr(), vfmaq_f64(vc, va, vb));
        }
    }
    crate::scalar::fma3(at, bt, ct, ot);
}

/// Σ x[i]·y[i]. FIXED reduction shape for this tier: two 2-lane fused
/// accumulators filled in index order over 4-wide blocks (acc0 ← low half,
/// acc1 ← high half), combined as (acc0 + acc1) then lane-summed low-to-high,
/// then the remainder appended via the scalar twin. Same input → same bits.
#[must_use]
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "dot length mismatch");
    let (xc, xt) = x.as_chunks::<{ 2 * LANES }>();
    let (yc, yt) = y.as_chunks::<{ 2 * LANES }>();
    // SAFETY: pointers into [f64; 4] arrays; `.add(LANES)` stays inside the
    // same 4-element array. Exact LANES extents per load.
    let vec_part = unsafe {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        for (xk, yk) in xc.iter().zip(yc) {
            acc0 = vfmaq_f64(acc0, vld1q_f64(xk.as_ptr()), vld1q_f64(yk.as_ptr()));
            acc1 = vfmaq_f64(
                acc1,
                vld1q_f64(xk.as_ptr().add(LANES)),
                vld1q_f64(yk.as_ptr().add(LANES)),
            );
        }
        vaddvq_f64(vaddq_f64(acc0, acc1))
    };
    vec_part + crate::scalar::dot(xt, yt)
}

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

/// Σ x[i]; same fixed two-accumulator shape as [`dot`].
#[must_use]
pub fn sum(x: &[f64]) -> f64 {
    let (xc, xt) = x.as_chunks::<{ 2 * LANES }>();
    // SAFETY: as in `dot`.
    let vec_part = unsafe {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        for xk in xc {
            acc0 = vaddq_f64(acc0, vld1q_f64(xk.as_ptr()));
            acc1 = vaddq_f64(acc1, vld1q_f64(xk.as_ptr().add(LANES)));
        }
        vaddvq_f64(vaddq_f64(acc0, acc1))
    };
    vec_part + crate::scalar::sum(xt)
}

/// Radix-4 Stockham q-run butterfly (bead 27d3): `vld2q`/`vst2q`
/// deinterleave two complex elements per iteration into (re, im) SoA
/// vregs, so every add/sub/mul_add below is the scalar twin's exact
/// per-element operation on two lanes at once — BITWISE. Runs whose
/// length is not a multiple of 4 f64 delegate to the twin whole.
pub fn r4qrun_f64(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
    out: &mut [f64],
    w: &[f64; 6],
    inverse: bool,
) {
    let s2 = a.len();
    assert!(
        s2.is_multiple_of(2)
            && b.len() == s2
            && c.len() == s2
            && d.len() == s2
            && out.len() == 4 * s2,
        "r4qrun run-length mismatch (programmer error)"
    );
    if !s2.is_multiple_of(4) {
        crate::scalar::r4qrun_f64(a, b, c, d, out, w, inverse);
        return;
    }
    // SAFETY: every vld2q/vst2q touches exactly 4 f64 at offset 4·q2
    // with 4·q2 + 4 ≤ s2 (asserted above, s2 % 4 == 0); the four output
    // rows live at disjoint offsets j·s2 within `out` (len 4·s2). f64
    // has no invalid bit patterns; unaligned access is permitted.
    unsafe {
        let (v1re, v1im) = (vdupq_n_f64(w[0]), vdupq_n_f64(w[1]));
        let (v2re, v2im) = (vdupq_n_f64(w[2]), vdupq_n_f64(w[3]));
        let (v3re, v3im) = (vdupq_n_f64(w[4]), vdupq_n_f64(w[5]));
        let op = out.as_mut_ptr();
        for q2 in 0..s2 / 4 {
            let o = 4 * q2;
            let av = vld2q_f64(a.as_ptr().add(o));
            let bv = vld2q_f64(b.as_ptr().add(o));
            let cv = vld2q_f64(c.as_ptr().add(o));
            let dv = vld2q_f64(d.as_ptr().add(o));
            let t0re = vaddq_f64(av.0, cv.0);
            let t0im = vaddq_f64(av.1, cv.1);
            let t1re = vsubq_f64(av.0, cv.0);
            let t1im = vsubq_f64(av.1, cv.1);
            let t2re = vaddq_f64(bv.0, dv.0);
            let t2im = vaddq_f64(bv.1, dv.1);
            let t3re = vsubq_f64(bv.0, dv.0);
            let t3im = vsubq_f64(bv.1, dv.1);
            let (t3ire, t3iim) = if inverse {
                (vnegq_f64(t3im), t3re)
            } else {
                (t3im, vnegq_f64(t3re))
            };
            vst2q_f64(
                op.add(o),
                float64x2x2_t(vaddq_f64(t0re, t2re), vaddq_f64(t0im, t2im)),
            );
            let u1re = vaddq_f64(t1re, t3ire);
            let u1im = vaddq_f64(t1im, t3iim);
            vst2q_f64(
                op.add(s2 + o),
                float64x2x2_t(
                    vfmaq_f64(vnegq_f64(vmulq_f64(u1im, v1im)), u1re, v1re),
                    vfmaq_f64(vmulq_f64(u1im, v1re), u1re, v1im),
                ),
            );
            let u2re = vsubq_f64(t0re, t2re);
            let u2im = vsubq_f64(t0im, t2im);
            vst2q_f64(
                op.add(2 * s2 + o),
                float64x2x2_t(
                    vfmaq_f64(vnegq_f64(vmulq_f64(u2im, v2im)), u2re, v2re),
                    vfmaq_f64(vmulq_f64(u2im, v2re), u2re, v2im),
                ),
            );
            let u3re = vsubq_f64(t1re, t3ire);
            let u3im = vsubq_f64(t1im, t3iim);
            vst2q_f64(
                op.add(3 * s2 + o),
                float64x2x2_t(
                    vfmaq_f64(vnegq_f64(vmulq_f64(u3im, v3im)), u3re, v3re),
                    vfmaq_f64(vmulq_f64(u3im, v3re), u3re, v3im),
                ),
            );
        }
    }
}

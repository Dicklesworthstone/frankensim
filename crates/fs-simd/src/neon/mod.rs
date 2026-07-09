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
    float64x2_t, vaddq_f64, vaddvq_f64, vdupq_n_f64, vfmaq_f64, vfmaq_laneq_f64, vld1q_f64,
    vmulq_f64, vst1q_f64,
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

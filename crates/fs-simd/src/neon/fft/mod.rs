//! NEON FFT stage-kernel capsule (split from neon/mod.rs under the
//! 300-line capsule cap, bead 8nfp): the radix-4 Stockham q-run
//! butterfly (bead 27d3). Registered in unsafe-capsules.json;
//! SAFETY.md beside this file. Bitwise contract and the tier battery
//! are unchanged — this is a pure file move.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use core::arch::aarch64::{
    float64x2x2_t, vaddq_f64, vdupq_n_f64, vfmaq_f64, vld2q_f64, vmulq_f64, vnegq_f64, vst2q_f64,
    vsubq_f64,
};

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

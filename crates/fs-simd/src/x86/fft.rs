//! x86-64 FFT stage-kernel capsule (bead 27d3): the radix-4 Stockham
//! q-run butterfly, AVX2+FMA variant, twin of [`crate::scalar::r4qrun_f64`].
//! Registered in unsafe-capsules.json; SAFETY.md beside this file.
//!
//! Feature-gating contract (identical to the sibling `x86/mod.rs` capsule):
//! the `#[target_feature]` inner function is reached ONLY through the safe
//! façade below, which re-checks avx2+fma at runtime and falls back to the
//! scalar twin otherwise — so the façade is unconditionally safe to call.
//!
//! Bitwise contract: the AVX2 path deinterleaves each 4-complex (8 f64)
//! chunk into (re, im) SoA lanes with `unpack`/`permute4x64` (pure data
//! movement, lossless), performs EXACTLY the scalar twin's per-element
//! add/sub/mul_add composition on four lanes at once — fused re part via
//! `_mm256_fmadd_pd`, the separate `im·w` product via `_mm256_mul_pd`, the
//! `-(…)` negation via a sign-bit XOR (bit-exact, matches Rust unary `-`
//! and NEON `vneg` on ±0) — then re-interleaves and stores. Runs whose
//! length is not a multiple of 8 f64 (4 complex) delegate to the twin
//! whole, exactly as the NEON capsule delegates the non-multiple-of-4 tail.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::{
    __m256d, _mm256_add_pd, _mm256_fmadd_pd, _mm256_loadu_pd, _mm256_mul_pd, _mm256_permute4x64_pd,
    _mm256_set1_pd, _mm256_storeu_pd, _mm256_sub_pd, _mm256_unpackhi_pd, _mm256_unpacklo_pd,
    _mm256_xor_pd,
};

/// Safe façade: AVX2+FMA radix-4 q-run butterfly, else the scalar twin.
/// Unconditionally safe — the feature is re-checked here at runtime.
pub fn r4qrun_f64(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
    out: &mut [f64],
    w: &[f64; 6],
    inverse: bool,
) {
    assert_r4qrun_bounds(a, b, c, d, out);
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified on this CPU immediately above; the
            // inner body's loads/stores are bounds-argued in its own block.
            return unsafe { r4qrun_256(a, b, c, d, out, w, inverse) };
        }
    }
    crate::scalar::r4qrun_f64(a, b, c, d, out, w, inverse);
}

fn assert_r4qrun_bounds(a: &[f64], b: &[f64], c: &[f64], d: &[f64], out: &[f64]) {
    let s2 = a.len();
    let out_len = crate::checked_r4qrun_output_len(s2);
    assert!(
        s2.is_multiple_of(2)
            && b.len() == s2
            && c.len() == s2
            && d.len() == s2
            && out_len == Some(out.len()),
        "r4qrun run-length mismatch (programmer error)"
    );
}

/// AVX2+FMA body: four complex elements (8 f64) per iteration.
///
/// # Safety
/// Requires avx2+fma and slice geometry established by the safe façade. Every
/// `loadu`/`storeu` touches
/// exactly 4 f64 at an offset `o` (resp. `o + 4`) with `o + 8 ≤ s2`
/// (loop bound `q8 < s2/8`, `s2 % 8 == 0`); the four output rows live at
/// disjoint offsets `j·s2` within `out` (len `4·s2`). f64 has no invalid
/// bit patterns and unaligned access is permitted.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn r4qrun_256(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
    out: &mut [f64],
    w: &[f64; 6],
    inverse: bool,
) {
    let s2 = a.len();
    if !s2.is_multiple_of(8) {
        crate::scalar::r4qrun_f64(a, b, c, d, out, w, inverse);
        return;
    }
    // SAFETY: all intrinsics below run under the verified avx2+fma feature.
    unsafe {
        let sign = _mm256_set1_pd(-0.0);
        let neg = |x: __m256d| _mm256_xor_pd(x, sign);
        let (v1re, v1im) = (_mm256_set1_pd(w[0]), _mm256_set1_pd(w[1]));
        let (v2re, v2im) = (_mm256_set1_pd(w[2]), _mm256_set1_pd(w[3]));
        let (v3re, v3im) = (_mm256_set1_pd(w[4]), _mm256_set1_pd(w[5]));
        // Deinterleave [re0,im0,re1,im1,re2,im2,re3,im3] -> (re,im) SoA.
        let deint = |p: *const f64| -> (__m256d, __m256d) {
            let x0 = _mm256_loadu_pd(p);
            let x1 = _mm256_loadu_pd(p.add(4));
            let lo = _mm256_unpacklo_pd(x0, x1); // [re0,re2,re1,re3]
            let hi = _mm256_unpackhi_pd(x0, x1); // [im0,im2,im1,im3]
            // permute [0,2,1,3] (imm 0xD8) -> ascending element order.
            (
                _mm256_permute4x64_pd::<0xD8>(lo),
                _mm256_permute4x64_pd::<0xD8>(hi),
            )
        };
        // Re-interleave (re,im) SoA -> AoS and store the two 256-bit halves.
        let inter_store = |p: *mut f64, re: __m256d, im: __m256d| {
            let tre = _mm256_permute4x64_pd::<0xD8>(re); // [r0,r2,r1,r3]
            let tim = _mm256_permute4x64_pd::<0xD8>(im); // [i0,i2,i1,i3]
            _mm256_storeu_pd(p, _mm256_unpacklo_pd(tre, tim)); // [r0,i0,r1,i1]
            _mm256_storeu_pd(p.add(4), _mm256_unpackhi_pd(tre, tim)); // [r2,i2,r3,i3]
        };
        let op = out.as_mut_ptr();
        for q8 in 0..s2 / 8 {
            let o = 8 * q8;
            let (are, aim) = deint(a.as_ptr().add(o));
            let (bre, bim) = deint(b.as_ptr().add(o));
            let (cre, cim) = deint(c.as_ptr().add(o));
            let (dre, dim) = deint(d.as_ptr().add(o));
            let t0re = _mm256_add_pd(are, cre);
            let t0im = _mm256_add_pd(aim, cim);
            let t1re = _mm256_sub_pd(are, cre);
            let t1im = _mm256_sub_pd(aim, cim);
            let t2re = _mm256_add_pd(bre, dre);
            let t2im = _mm256_add_pd(bim, dim);
            let t3re = _mm256_sub_pd(bre, dre);
            let t3im = _mm256_sub_pd(bim, dim);
            let (t3ire, t3iim) = if inverse {
                (neg(t3im), t3re)
            } else {
                (t3im, neg(t3re))
            };
            // o0 = t0 + t2 (no twiddle)
            inter_store(
                op.add(o),
                _mm256_add_pd(t0re, t2re),
                _mm256_add_pd(t0im, t2im),
            );
            // o1 = (t1 + t3i) · w1
            let u1re = _mm256_add_pd(t1re, t3ire);
            let u1im = _mm256_add_pd(t1im, t3iim);
            inter_store(
                op.add(s2 + o),
                _mm256_fmadd_pd(u1re, v1re, neg(_mm256_mul_pd(u1im, v1im))),
                _mm256_fmadd_pd(u1re, v1im, _mm256_mul_pd(u1im, v1re)),
            );
            // o2 = (t0 - t2) · w2
            let u2re = _mm256_sub_pd(t0re, t2re);
            let u2im = _mm256_sub_pd(t0im, t2im);
            inter_store(
                op.add(2 * s2 + o),
                _mm256_fmadd_pd(u2re, v2re, neg(_mm256_mul_pd(u2im, v2im))),
                _mm256_fmadd_pd(u2re, v2im, _mm256_mul_pd(u2im, v2re)),
            );
            // o3 = (t1 - t3i) · w3
            let u3re = _mm256_sub_pd(t1re, t3ire);
            let u3im = _mm256_sub_pd(t1im, t3iim);
            inter_store(
                op.add(3 * s2 + o),
                _mm256_fmadd_pd(u3re, v3re, neg(_mm256_mul_pd(u3im, v3im))),
                _mm256_fmadd_pd(u3re, v3im, _mm256_mul_pd(u3im, v3re)),
            );
        }
    }
}

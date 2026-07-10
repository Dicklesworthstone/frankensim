//! x86-64 batched-GEMM tile capsule (bead 9ekv): the PACKED 4×4 batched-GEMM
//! tile microkernel, AVX2+FMA variant, twin of [`crate::scalar::btile4x4p_f64`].
//! Registered in unsafe-capsules.json; SAFETY.md beside this file.
//!
//! Feature-gating contract (identical to the sibling `x86/mod.rs` capsule):
//! the `#[target_feature]` inner function is reached ONLY through the safe
//! façade below, which re-checks avx2+fma at runtime and falls back to the
//! scalar twin otherwise — so the façade is unconditionally safe to call.
//!
//! Bitwise contract: over the lane (`mb`) dimension the AVX2 body accumulates
//! each of the 16 output tile elements in a 4-lane `__m256d` starting from
//! `_mm256_setzero_pd()` (+0.0) with `_mm256_fmadd_pd` in l-ascending order —
//! EXACTLY the scalar twin's per-lane `mul_add` from a +0.0 start — so every
//! lane is bit-identical. Lanes past the last full group of 4 (`mb % 4`) run
//! the scalar per-lane loop. The 4×4 tile keeps 16 live accumulators (> the 16
//! YMM registers), so LLVM spills; this is a correctness-first vectorization of
//! the lane dimension, not a register-optimal kernel.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::{
    _mm256_fmadd_pd, _mm256_loadu_pd, _mm256_setzero_pd, _mm256_storeu_pd,
};

/// Safe façade: AVX2+FMA packed 4×4 batched-GEMM tile, else the scalar twin.
/// Unconditionally safe — the feature is re-checked here at runtime.
#[allow(clippy::too_many_arguments)] // packed-layout bundle (matches the twin)
pub fn btile4x4p_f64(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f64],
) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified on this CPU immediately above; the
            // inner body's loads/stores are bounds-argued in its own block.
            return unsafe { btile4x4p_256(a, b, i0, j0, k, mb, dst) };
        }
    }
    crate::scalar::btile4x4p_f64(a, b, i0, j0, k, mb, dst);
}

/// AVX2+FMA body: 4 lanes (`__m256d`) of the packed batched tile per block.
///
/// # Safety
/// Requires avx2+fma (façade-verified). Every pointer is `base(t) + l·mb + m`
/// with `base(t)` = `((i0+t)·k)·mb` (a) or `((j0+t)·k)·mb` (b); the maximal
/// dereferenced offset over `t ≤ 3`, `l ≤ k−1`, `m ≤ mb−4` (vector) or
/// `≤ mb−1` (scalar tail) is inside the extents asserted by
/// `checked_btile4x4p_lengths`. Each vector access is 4 f64; f64 has no
/// invalid bit patterns and unaligned access is permitted; the 16 output rows
/// live at disjoint offsets `(ti·4+tj)·mb` within `dst`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
#[allow(clippy::too_many_arguments)]
unsafe fn btile4x4p_256(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f64],
) {
    let bounds = crate::checked_btile4x4p_lengths(i0, j0, k, mb);
    assert!(
        matches!(
            bounds,
            Some((a_len, b_len, dst_len))
                if a_len <= a.len() && b_len <= b.len() && dst_len <= dst.len()
        ),
        "btile4x4p packed bounds (programmer error)"
    );
    let full = (mb / 4) * 4;
    // SAFETY: all pointer arithmetic and vector ops below run under the
    // verified avx2+fma feature; offsets are bounded as argued on the fn.
    unsafe {
        let ap0 = a.as_ptr();
        let bp0 = b.as_ptr();
        let op = dst.as_mut_ptr();
        let a_base = [
            ap0.add(i0 * k * mb),
            ap0.add((i0 + 1) * k * mb),
            ap0.add((i0 + 2) * k * mb),
            ap0.add((i0 + 3) * k * mb),
        ];
        let b_base = [
            bp0.add(j0 * k * mb),
            bp0.add((j0 + 1) * k * mb),
            bp0.add((j0 + 2) * k * mb),
            bp0.add((j0 + 3) * k * mb),
        ];
        let mut m = 0;
        while m < full {
            let mut acc = [_mm256_setzero_pd(); 16];
            let mut ap = [
                a_base[0].add(m),
                a_base[1].add(m),
                a_base[2].add(m),
                a_base[3].add(m),
            ];
            let mut bp = [
                b_base[0].add(m),
                b_base[1].add(m),
                b_base[2].add(m),
                b_base[3].add(m),
            ];
            for _l in 0..k {
                let av = [
                    _mm256_loadu_pd(ap[0]),
                    _mm256_loadu_pd(ap[1]),
                    _mm256_loadu_pd(ap[2]),
                    _mm256_loadu_pd(ap[3]),
                ];
                let bv = [
                    _mm256_loadu_pd(bp[0]),
                    _mm256_loadu_pd(bp[1]),
                    _mm256_loadu_pd(bp[2]),
                    _mm256_loadu_pd(bp[3]),
                ];
                for ti in 0..4 {
                    for tj in 0..4 {
                        acc[ti * 4 + tj] = _mm256_fmadd_pd(av[ti], bv[tj], acc[ti * 4 + tj]);
                    }
                }
                for t in 0..4 {
                    ap[t] = ap[t].add(mb);
                    bp[t] = bp[t].add(mb);
                }
            }
            for ti in 0..4 {
                for tj in 0..4 {
                    _mm256_storeu_pd(op.add((ti * 4 + tj) * mb + m), acc[ti * 4 + tj]);
                }
            }
            m += 4;
        }
        // Scalar tail for the mb % 4 lanes past the last full group — the
        // scalar twin's exact per-lane l-ascending fused accumulation.
        for ti in 0..4 {
            for tj in 0..4 {
                let out_base = (ti * 4 + tj) * mb;
                for lane in full..mb {
                    let mut s = 0.0f64;
                    for l in 0..k {
                        let am = *a_base[ti].add(l * mb + lane);
                        let bm = *b_base[tj].add(l * mb + lane);
                        s = am.mul_add(bm, s);
                    }
                    *op.add(out_base + lane) = s;
                }
            }
        }
    }
}

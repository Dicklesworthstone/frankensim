//! x86-64 GEMM capsules: the dense 8×4 register microkernel and the packed
//! 4×4 batched-GEMM tile, both AVX2+FMA twins of their scalar definitions.
//! Registered in unsafe-capsules.json; SAFETY.md beside this file.
//!
//! Feature-gating contract: public safe façades re-check AVX2+FMA on every
//! call. The process-wide function table instead calls [`select_mk8x4_f64`]
//! once; only that selector can name and return the private unchecked thunk.
//!
//! Bitwise contract: over the lane (`mb`) dimension the AVX2 body accumulates
//! each of the 16 output tile elements in a 4-lane `__m256d` starting from
//! `_mm256_setzero_pd()` (+0.0) with `_mm256_fmadd_pd` in l-ascending order —
//! EXACTLY the scalar twin's per-lane `mul_add` from a +0.0 start — so every
//! lane is bit-identical. Lanes past the last full group of 4 (`mb % 4`) run
//! the scalar per-lane loop. The 4×4 tile runs as two ti-half passes of 8
//! accumulators so the live set (8 acc + 4 b + 2 a vectors) fits the 16-entry
//! YMM file; the second pass re-reads the L1-resident B stream, and the split
//! is bit-neutral because each output element's l-ascending fused chain is
//! untouched.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::{
    __m256d, _mm256_fmadd_pd, _mm256_loadu_pd, _mm256_set1_pd, _mm256_setzero_pd, _mm256_storeu_pd,
};
use fs_substrate::SimdTier;

/// Resolve the dense GEMM microkernel and its truthful operation tier once.
///
/// The returned function pointer is safe because the only unchecked thunk is
/// private to this module and is returned only after the required CPU features
/// are observed. The scalar pointer is returned for every other combination.
pub(crate) fn select_mk8x4_f64(global_tier: SimdTier) -> (crate::Mk8x4, SimdTier) {
    let avx2_available = std::arch::is_x86_feature_detected!("avx2");
    let fma_available = std::arch::is_x86_feature_detected!("fma");
    let tier = crate::mk8x4_f64_tier_for(global_tier, avx2_available, fma_available);
    if tier == SimdTier::Avx2 {
        (mk8x4_f64_selected, tier)
    } else {
        (crate::scalar::mk8x4_f64, SimdTier::Scalar)
    }
}

/// Safe façade for the 8×4 microkernel: AVX2+FMA body, scalar twin fallback.
pub fn mk8x4_f64(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    assert_mk8x4_bounds(a_panel, b_panel, kc);
    if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma") {
        return mk8x4_f64_selected(a_panel, b_panel, kc, acc);
    }
    crate::scalar::mk8x4_f64(a_panel, b_panel, kc, acc);
}

fn mk8x4_f64_selected(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    assert_mk8x4_bounds(a_panel, b_panel, kc);
    // SAFETY: only `select_mk8x4_f64`, after one-shot AVX2+FMA detection, can
    // export this private thunk; the public façade checks on every call.
    unsafe { mk8x4_f64_256(a_panel, b_panel, kc, acc) }
}

fn assert_mk8x4_bounds(a_panel: &[f64], b_panel: &[f64], kc: usize) {
    let lengths = crate::checked_mk8x4_lengths(kc);
    assert!(
        matches!(lengths, Some((a_len, b_len)) if a_len <= a_panel.len() && b_len <= b_panel.len()),
        "mk8x4 panel length mismatch (programmer error)"
    );
}

/// 8×4 f64 GEMM register microkernel, AVX2+FMA body with k-ascending fused
/// accumulation exactly matching the scalar twin.
///
/// # Safety
/// Requires AVX2+FMA and panel bounds established by `assert_mk8x4_bounds`.
#[target_feature(enable = "avx2,fma")]
unsafe fn mk8x4_f64_256(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    // SAFETY: load/store extents are discharged by the checked panel bounds.
    unsafe {
        let ap = a_panel.as_ptr();
        let bp = b_panel.as_ptr();
        let mut va: [__m256d; 8] = std::array::from_fn(|r| _mm256_loadu_pd(acc[r].as_ptr()));
        for kk in 0..kc {
            let b = _mm256_loadu_pd(bp.add(kk * 4));
            for (r, v) in va.iter_mut().enumerate() {
                let ar = _mm256_set1_pd(*ap.add(kk * 8 + r));
                *v = _mm256_fmadd_pd(ar, b, *v);
            }
        }
        for (r, v) in va.iter().enumerate() {
            _mm256_storeu_pd(acc[r].as_mut_ptr(), *v);
        }
    }
}

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
    assert_btile4x4p_bounds(a, b, i0, j0, k, mb, dst);
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

#[allow(clippy::too_many_arguments)]
fn assert_btile4x4p_bounds(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &[f64],
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
}

/// AVX2+FMA body: 4 lanes (`__m256d`) of the packed batched tile per block.
///
/// # Safety
/// Requires avx2+fma and slice geometry established by the safe façade. Every
/// pointer is `base(t) + l·mb + m`
/// with `base(t)` = `((i0+t)·k)·mb` (a) or `((j0+t)·k)·mb` (b); the maximal
/// dereferenced offset over `t ≤ 3`, `l ≤ k−1`, `m ≤ mb−4` (vector) or
/// `≤ mb−1` (scalar tail) is inside the extents asserted by
/// `checked_btile4x4p_lengths`. Each vector access is 4 f64; f64 has no
/// invalid bit patterns and unaligned access is permitted; the 16 output rows
/// live at disjoint offsets `(ti·4+tj)·mb` within `dst`. Each lane block
/// starts from validated bases and advances its cursors only when another
/// reduction load remains, so it never forms one-past plus `m`.
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
            // Two ti-half passes of 8 accumulators each: a single 16-vector
            // accumulator block plus 8 operand vectors exceeds the 16-entry
            // YMM file and spills every l iteration (measured ~0.27-0.34
            // attainment ceiling on Zen 3). 8 acc + 4 bv + 2 av = 14 live
            // vectors fits; the second pass re-reads the L1-resident B
            // stream. Per-element accumulation order is UNCHANGED (each
            // output's chain stays l-ascending fused), so the split is
            // bit-neutral by construction.
            for half in 0..2 {
                let ti0 = half * 2;
                let mut acc = [_mm256_setzero_pd(); 8];
                let mut ap = [a_base[ti0].add(m), a_base[ti0 + 1].add(m)];
                let mut bp = [
                    b_base[0].add(m),
                    b_base[1].add(m),
                    b_base[2].add(m),
                    b_base[3].add(m),
                ];
                for l in 0..k {
                    let av = [_mm256_loadu_pd(ap[0]), _mm256_loadu_pd(ap[1])];
                    let bv = [
                        _mm256_loadu_pd(bp[0]),
                        _mm256_loadu_pd(bp[1]),
                        _mm256_loadu_pd(bp[2]),
                        _mm256_loadu_pd(bp[3]),
                    ];
                    for ti in 0..2 {
                        for tj in 0..4 {
                            acc[ti * 4 + tj] = _mm256_fmadd_pd(av[ti], bv[tj], acc[ti * 4 + tj]);
                        }
                    }
                    if l + 1 < k {
                        for ptr in &mut ap {
                            *ptr = (*ptr).add(mb);
                        }
                        for ptr in &mut bp {
                            *ptr = (*ptr).add(mb);
                        }
                    }
                }
                for ti in 0..2 {
                    for tj in 0..4 {
                        _mm256_storeu_pd(op.add(((ti0 + ti) * 4 + tj) * mb + m), acc[ti * 4 + tj]);
                    }
                }
            }
            m += 4;
        }
        // Scalar tail for the mb % 4 lanes past the last full group — the
        // scalar twin's exact per-lane l-ascending fused accumulation.
        for (ti, &a_ptr) in a_base.iter().enumerate() {
            for (tj, &b_ptr) in b_base.iter().enumerate() {
                let out_base = (ti * 4 + tj) * mb;
                for lane in full..mb {
                    let mut s = 0.0f64;
                    for l in 0..k {
                        let am = *a_ptr.add(l * mb + lane);
                        let bm = *b_ptr.add(l * mb + lane);
                        s = am.mul_add(bm, s);
                    }
                    *op.add(out_base + lane) = s;
                }
            }
        }
    }
}

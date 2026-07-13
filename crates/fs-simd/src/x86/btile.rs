//! AVX2+FMA packed batched f32 tile (bead 9ekv, the x86 row's f32
//! sibling of `x86::gemm::btile4x4p_f64`): eight lanes per `__m256`,
//! the same two ti-half register-blocking (8 accumulators + 6 operand
//! vectors fits the 16-entry YMM file; one 16-accumulator block would
//! spill every reduction step), and the same scalar tail. Without this
//! capsule the f32 path dispatched to the scalar twin, whose
//! `f32::mul_add` lowers to a per-element libm call on baseline
//! x86-64 (the fz2.2/a55x hole).
//!
//! BITWISE contract: lanes are independent matrices; per lane every
//! path is the scalar twin's zero-start, l-ascending fused chain, so
//! vector width and the ti-half split cannot move bits. Registered in
//! unsafe-capsules.json; safety argument in the shared
//! `crates/fs-simd/src/x86/SAFETY.md`.
#![allow(unsafe_code)] // registered capsule — see x86/SAFETY.md

use core::arch::x86_64::{_mm256_fmadd_ps, _mm256_loadu_ps, _mm256_setzero_ps, _mm256_storeu_ps};

/// PACKED batched 4×4 f32 tile over l-contiguous operands (A i-major,
/// B j-major, both walks stride `mb`): safe façade — AVX2+FMA when the
/// host has it, the scalar twin otherwise. BITWISE across paths.
///
/// # Panics
/// Structured panics on packed-extent mismatches (the twin's contract).
pub fn btile4x4pf32(
    a: &[f32],
    b: &[f32],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f32],
) {
    let bounds = crate::checked_btile4x4p_lengths(i0, j0, k, mb);
    assert!(
        matches!(
            bounds,
            Some((a_len, b_len, dst_len))
                if a_len <= a.len() && b_len <= b.len() && dst_len <= dst.len()
        ),
        "btile4x4pf32 packed bounds (programmer error)"
    );
    if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma") {
        // SAFETY: feature availability re-verified immediately above.
        unsafe { btile4x4pf32_256(a, b, i0, j0, k, mb, dst) };
        return;
    }
    crate::scalar::btile4x4pf32(a, b, i0, j0, k, mb, dst);
}

/// AVX2+FMA body: 8 f32 lanes per `__m256` block, two ti-half passes.
///
/// # Safety
/// Requires avx2+fma and the slice geometry asserted by the façade.
/// Every pointer is `base(t) + l·mb + m` with `base(t)` =
/// `((i0+t)·k)·mb` (a) or `((j0+t)·k)·mb` (b); the maximal dereferenced
/// offset over `t ≤ 3`, `l ≤ k−1`, `m ≤ mb−8` (vector) or `≤ mb−1`
/// (scalar tail) is inside the asserted extents. Each vector access is
/// 8 f32, unaligned-tolerant (`loadu`/`storeu`); f32 has no invalid bit
/// patterns; the 16 output rows live at disjoint `(ti·4+tj)·mb`
/// offsets. Cursors advance only while another reduction load remains.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
#[allow(clippy::too_many_arguments)]
unsafe fn btile4x4pf32_256(
    a: &[f32],
    b: &[f32],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f32],
) {
    let full = (mb / 8) * 8;
    // SAFETY: all pointer arithmetic below is bounded as argued on the fn.
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
            // Two ti-half passes of 8 accumulators each (the f64
            // kernel's measured register-blocking, same YMM budget:
            // 8 acc + 4 bv + 2 av = 14 live vectors). Per-element
            // accumulation order is UNCHANGED — bit-neutral split.
            for half in 0..2 {
                let ti0 = half * 2;
                let mut acc = [_mm256_setzero_ps(); 8];
                let mut ap = [a_base[ti0].add(m), a_base[ti0 + 1].add(m)];
                let mut bp = [
                    b_base[0].add(m),
                    b_base[1].add(m),
                    b_base[2].add(m),
                    b_base[3].add(m),
                ];
                for l in 0..k {
                    let av = [_mm256_loadu_ps(ap[0]), _mm256_loadu_ps(ap[1])];
                    let bv = [
                        _mm256_loadu_ps(bp[0]),
                        _mm256_loadu_ps(bp[1]),
                        _mm256_loadu_ps(bp[2]),
                        _mm256_loadu_ps(bp[3]),
                    ];
                    for ti in 0..2 {
                        for tj in 0..4 {
                            acc[ti * 4 + tj] = _mm256_fmadd_ps(av[ti], bv[tj], acc[ti * 4 + tj]);
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
                        _mm256_storeu_ps(op.add(((ti0 + ti) * 4 + tj) * mb + m), acc[ti * 4 + tj]);
                    }
                }
            }
            m += 8;
        }
        // Scalar tail for the mb % 8 lanes past the last full group —
        // the scalar twin's exact per-lane l-ascending fused chain.
        for (ti, &a_ptr) in a_base.iter().enumerate() {
            for (tj, &b_ptr) in b_base.iter().enumerate() {
                let out_base = (ti * 4 + tj) * mb;
                for lane in full..mb {
                    let mut s = 0.0f32;
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

//! x86-64 elementwise capsule (bead fz2.2 tier audit): the fused
//! `fma3` vector path. Registered in unsafe-capsules.json; SAFETY.md
//! beside this file (the sibling capsules' contract applies verbatim).
//!
//! WHY THIS EXISTS: the tier audit measured `fma3` dispatched to the
//! scalar twin on x86 — and baseline x86-64 has no compile-time FMA,
//! so `f64::mul_add` lowers to a per-element libm CALL (the same
//! hazard the fs-roofline peak probe hit). `_mm256_fmadd_pd` is the
//! honest fused path.
//!
//! Bitwise contract: element-wise fused multiply-add — each lane is
//! EXACTLY `a[i].mul_add(b[i], c[i])` (both are a single IEEE fused
//! operation), so the vector path is bit-identical to the scalar twin
//! per element; the non-multiple-of-4 tail delegates to the twin.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::{_mm256_fmadd_pd, _mm256_loadu_pd, _mm256_storeu_pd};

/// Safe façade: AVX2+FMA fused `out[i] = a[i]·b[i] + c[i]`, else the
/// scalar twin. Unconditionally safe — features re-checked at runtime.
pub fn fma3(a: &[f64], b: &[f64], c: &[f64], out: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified on this CPU immediately above.
            return unsafe { fma3_256(a, b, c, out) };
        }
    }
    crate::scalar::fma3(a, b, c, out);
}

/// AVX2+FMA body: 4 lanes per iteration, tails to the scalar twin.
///
/// # Safety
/// Requires avx2+fma (verified by the façade). All loads/stores use
/// chunk-array pointers with exact 4-element extents.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn fma3_256(a: &[f64], b: &[f64], c: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "fma3 length mismatch");
    assert_eq!(a.len(), c.len(), "fma3 length mismatch");
    assert_eq!(a.len(), out.len(), "fma3 length mismatch");
    let (ac, at) = a.as_chunks::<4>();
    let (bc, bt) = b.as_chunks::<4>();
    let (cc, ct) = c.as_chunks::<4>();
    let (oc, ot) = out.as_chunks_mut::<4>();
    // SAFETY: chunk-array pointers, exact 4-lane extents per chunk.
    unsafe {
        for (((ak, bk), ck), ok) in ac.iter().zip(bc).zip(cc).zip(oc) {
            let va = _mm256_loadu_pd(ak.as_ptr());
            let vb = _mm256_loadu_pd(bk.as_ptr());
            let vc = _mm256_loadu_pd(ck.as_ptr());
            _mm256_storeu_pd(ok.as_mut_ptr(), _mm256_fmadd_pd(va, vb, vc));
        }
    }
    crate::scalar::fma3(at, bt, ct, ot);
}

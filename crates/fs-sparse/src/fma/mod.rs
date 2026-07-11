//! x86 FMA-codegen capsule for the CSR hot loops (bead nabk; the a55x
//! pattern): baseline x86-64 lowers every `f64::mul_add` in the spmv /
//! spmm reduction chains to a PER-ELEMENT libm `fma()` CALL (no
//! compile-time FMA) — the hazard class measured at 14–28x deficits on
//! the 5995WX across feec/batch-gemm/batch-lu. This capsule
//! re-compiles the SAME `#[inline(always)]` bodies under
//! `#[target_feature(enable = "avx2,fma")]`, where `mul_add` becomes a
//! native fused instruction. One correctly-rounded fused IEEE op per
//! element either way: BIT-IDENTICAL to the portable path (the crate's
//! cross-format bitwise suites gate it), and the reduction SHAPE —
//! ascending-column fused chain from +0.0 — is untouched: this is pure
//! codegen, never reordering. Registered capsule; SAFETY.md beside.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use crate::Csr;

/// Run spmv with FMA codegen when the CPU has it. Safe to call.
#[inline]
pub(crate) fn spmv_dispatch(csr: &Csr, x: &[f64], y: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { spmv_x86(csr, x, y) };
        }
    }
    csr.spmv_body(x, y);
}

/// Run spmm with FMA codegen when the CPU has it. Safe to call.
#[inline]
pub(crate) fn spmm_dispatch(csr: &Csr, x: &[f64], k: usize, y: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { spmm_x86(csr, x, k, y) };
        }
    }
    csr.spmm_body(x, k, y);
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before the
/// call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn spmv_x86(csr: &Csr, x: &[f64], y: &mut [f64]) {
    csr.spmv_body(x, y);
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before the
/// call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn spmm_x86(csr: &Csr, x: &[f64], k: usize, y: &mut [f64]) {
    csr.spmm_body(x, k, y);
}

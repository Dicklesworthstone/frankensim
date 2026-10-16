//! x86 FMA-codegen capsule for the fs-cheb hot loops (bead nabk; the
//! a55x pattern): baseline x86-64 lowers every `f64::mul_add` in the
//! Clenshaw / collocation reduction chains to a PER-ELEMENT libm
//! `fma()` CALL (no compile-time FMA) — the hazard class measured at
//! 14–28x deficits on the 5995WX across feec/batch-gemm/batch-lu.
//! This capsule re-compiles the SAME `#[inline(always)]` bodies under
//! `#[target_feature(enable = "avx2,fma")]`, where `mul_add` becomes a
//! native fused instruction. One correctly-rounded fused IEEE op per
//! element either way: BIT-IDENTICAL to the portable path (the crate's
//! batteries and the workspace goldens gate it), and every chain shape
//! — Clenshaw's backward recurrence, the j-inner dot chains, the
//! k-outer skip-zero saxpy — is untouched: pure codegen, never
//! reordering.
//!
//! Covered sites (the loop-hot subset of the nabk census): `Cheb1::
//! eval` (Clenshaw, called per evaluation point by rootfinding /
//! quadrature / cheb2), `dirichlet_laplace_eigs`' D·D product and
//! Rayleigh matvec, and `orr_sommerfeld`'s D-power matmul. The census'
//! remaining fs-cheb sites are deliberately NOT capsuled: `Fourier
//! Series::eval` is trig-call-bound (two `det::` calls per element
//! dominate its one `mul_add`), `Cheb1::differentiate` is a once-per-
//! call coefficient recurrence behind a Vec allocation, and the
//! `os_matrices` D4-clamp assembly is one O(n²) pass beside three
//! O(n³) matmuls.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use crate::Cheb1;

/// Clenshaw evaluation with FMA codegen when the CPU has it. Safe to
/// call; bit-identical either way.
#[inline]
pub(crate) fn cheb_eval_dispatch(c: &Cheb1, t: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { cheb_eval_x86(c, t) };
        }
    }
    c.eval_body(t)
}

/// D·D product with FMA codegen when the CPU has it. Safe to call;
/// bit-identical either way.
#[inline]
pub(crate) fn dsq_into_dispatch(d: &[f64], m: usize, d2: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { dsq_into_x86(d, m, d2) };
        }
    }
    crate::dsq_into_body(d, m, d2);
}

/// Dense matvec with FMA codegen when the CPU has it. Safe to call;
/// bit-identical either way.
#[inline]
pub(crate) fn matvec_into_dispatch(a: &[f64], v: &[f64], ni: usize, av: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { matvec_into_x86(a, v, ni, av) };
        }
    }
    crate::matvec_into_body(a, v, ni, av);
}

/// Dense n×n matmul with FMA codegen when the CPU has it. Safe to
/// call; bit-identical either way.
#[inline]
pub(crate) fn os_matmul_dispatch(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified immediately above; the inlined
            // body is pure safe slice arithmetic.
            return unsafe { os_matmul_x86(a, b, n) };
        }
    }
    crate::orr_sommerfeld::matmul_body(a, b, n)
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn cheb_eval_x86(c: &Cheb1, t: f64) -> f64 {
    c.eval_body(t)
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn dsq_into_x86(d: &[f64], m: usize, d2: &mut [f64]) {
    crate::dsq_into_body(d, m, d2);
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn matvec_into_x86(a: &[f64], v: &[f64], ni: usize, av: &mut [f64]) {
    crate::matvec_into_body(a, v, ni, av);
}

/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn os_matmul_x86(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    crate::orr_sommerfeld::matmul_body(a, b, n)
}

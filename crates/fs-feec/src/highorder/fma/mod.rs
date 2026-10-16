//! x86 FMA-codegen capsule for the sum-factorized apply (bead a55x):
//! baseline x86-64 lowers every `f64::mul_add` in the const-P element
//! loop to a PER-ELEMENT libm `fma()` CALL (no compile-time FMA),
//! measured at 0.026 attainment on the 5995WX vs 0.372 on aarch64 —
//! the same hazard class as the fs-roofline peak probe and the
//! fs-simd fma3 dispatch hole. This capsule re-compiles the SAME
//! `#[inline(always)]` kernel body under
//! `#[target_feature(enable = "avx2,fma")]`, where `mul_add` becomes a
//! native fused instruction. One fused IEEE op per element either way:
//! BIT-IDENTICAL to the portable path (gated by the crate's sf-kron
//! golden, which does not move). Registered capsule; SAFETY.md beside.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use super::hex::TensorSpace;

/// Run the const-P element loop with FMA codegen when the CPU has it,
/// portable codegen otherwise. Unconditionally safe to call.
#[inline]
pub(crate) fn apply_mono_dispatch<const P: usize>(space: &TensorSpace, u: &[f64], y: &mut [f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") && std::arch::is_x86_feature_detected!("fma")
        {
            // SAFETY: avx2+fma verified on this CPU immediately above;
            // the inlined body is pure safe slice arithmetic.
            return unsafe { apply_mono_x86::<P>(space, u, y) };
        }
    }
    space.apply_mono_body::<P>(u, y);
}

/// The identical body compiled with FMA enabled (the `inline(always)`
/// on `apply_mono_body` makes its code generate HERE, under these
/// target features — a non-inlined call would keep baseline codegen).
///
/// # Safety
/// Requires avx2+fma, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn apply_mono_x86<const P: usize>(space: &TensorSpace, u: &[f64], y: &mut [f64]) {
    space.apply_mono_body::<P>(u, y);
}

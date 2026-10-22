//! x86 AVX2-codegen capsule for the complex FMM's M2L dot product. Baseline
//! x86-64 runs `dot_body`'s four partial-sum lanes two at a time (SSE2,
//! `mulpd`/`addpd`). This capsule recompiles the SAME `#[inline(always)]`
//! body under `#[target_feature(enable = "avx2")]`, so the lanes run four wide.
//! The body's products are unfused (AVX2 does not enable FMA, and Rust never
//! contracts a separate multiply and add) and each lane keeps its operation
//! order, so the results are BIT-IDENTICAL to the portable path. The crate's
//! tests pin that against the portable body. Registered capsule; SAFETY.md
//! beside. Other CPUs run the portable body.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use super::{C64, dot_body};

/// Whether this CPU runs the AVX2-compiled body.
#[cfg(test)]
pub(super) fn accelerated() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// The M2L dot product with AVX2 codegen when the CPU has it, portable
/// codegen otherwise. Unconditionally safe to call.
#[inline]
pub(super) fn dot(a: &[C64], b: &[C64]) -> C64 {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: avx2 verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            return unsafe { dot_avx2(a, b) };
        }
    }
    dot_body(a, b)
}

/// The identical dot body compiled with AVX2 enabled (its `inline(always)`
/// makes the code generate HERE, under this feature).
///
/// # Safety
/// Requires the `avx2` feature, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dot_avx2(a: &[C64], b: &[C64]) -> C64 {
    dot_body(a, b)
}

//! x86 FMA-codegen capsule for the complex LU (WhiteAnchor mail 5404). Baseline
//! x86-64 lowers every `f64::mul_add` (two per complex product) to a call of the
//! libm/compiler-builtins `fma()`: a 2080-panel Burton-Miller sample spent about
//! 700 s in this factorization on trj. This capsule recompiles the SAME
//! `#[inline(always)]` factor and solve bodies under
//! `#[target_feature(enable = "fma")]`, where `mul_add` is one native fused
//! instruction. IEEE-754 fma is correctly rounded under both codegens, so the
//! results are BIT-IDENTICAL to the portable path, pinned by the crate's tests
//! against the frozen original operation graph. The blocked elimination's
//! panel, upper-block and trailing bodies are recompiled the same way.
//! Registered capsule; SAFETY.md beside. aarch64 and x86-64 CPUs without FMA
//! run the portable bodies.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use super::blocked::{panel_body, trailing_body, upper_body};
use super::{C64, LuComplex, LuComplexError, factor_body, solve_body};

/// Whether this CPU runs the FMA-compiled bodies.
#[cfg(test)]
pub(super) fn accelerated() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("fma")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Eliminate with FMA codegen when the CPU has it, portable codegen otherwise.
/// Unconditionally safe to call.
#[inline]
pub(super) fn factor(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("fma") {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            return unsafe { factor_fma(m, perm, n, cancelled) };
        }
    }
    factor_body(m, perm, n, cancelled)
}

/// Substitute with FMA codegen when the CPU has it, portable codegen otherwise.
/// Unconditionally safe to call.
#[inline]
pub(super) fn solve(
    lu: &LuComplex,
    b: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, LuComplexError> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("fma") {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            return unsafe { solve_fma(lu, b, cancelled) };
        }
    }
    solve_body(lu, b, cancelled)
}

/// Blocked panel with FMA codegen when the CPU has it. Unconditionally safe.
#[inline]
pub(super) fn panel(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    k0: usize,
    k1: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("fma") {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            return unsafe { panel_fma(m, perm, n, k0, k1, cancelled) };
        }
    }
    panel_body(m, perm, n, k0, k1, cancelled)
}

/// One U12 column chunk with FMA codegen when the CPU has it. Safe to call.
#[inline]
pub(super) fn upper(rows: &mut [(&[C64], &mut [C64])]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("fma") {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            unsafe { upper_fma(rows) };
            return;
        }
    }
    upper_body(rows);
}

/// One trailing row chunk with FMA codegen when the CPU has it. Safe to call.
#[inline]
pub(super) fn trailing(rows: &mut [C64], above: &[C64], n: usize, k0: usize, k1: usize) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("fma") {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe slice arithmetic.
            unsafe { trailing_fma(rows, above, n, k0, k1) };
            return;
        }
    }
    trailing_body(rows, above, n, k0, k1);
}

/// The identical elimination body compiled with FMA enabled (its
/// `inline(always)` makes the code generate HERE, under this feature).
///
/// # Safety
/// Requires the `fma` feature, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn factor_fma(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    factor_body(m, perm, n, cancelled)
}

/// The identical substitution body compiled with FMA enabled.
///
/// # Safety
/// Requires the `fma` feature, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn solve_fma(
    lu: &LuComplex,
    b: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, LuComplexError> {
    solve_body(lu, b, cancelled)
}

/// The blocked panel body compiled with FMA enabled.
///
/// # Safety
/// Requires the `fma` feature, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn panel_fma(
    m: &mut [C64],
    perm: &mut [usize],
    n: usize,
    k0: usize,
    k1: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), LuComplexError> {
    panel_body(m, perm, n, k0, k1, cancelled)
}

/// The U12 chunk body compiled with FMA enabled.
///
/// # Safety
/// As for [`panel_fma`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn upper_fma(rows: &mut [(&[C64], &mut [C64])]) {
    upper_body(rows);
}

/// The trailing chunk body compiled with FMA enabled.
///
/// # Safety
/// As for [`panel_fma`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn trailing_fma(rows: &mut [C64], above: &[C64], n: usize, k0: usize, k1: usize) {
    trailing_body(rows, above, n, k0, k1);
}

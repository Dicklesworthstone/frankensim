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

use crate::bsr::Bsr;
use crate::perf::CsrCompact;
use crate::sell::Sell;

macro_rules! fma_dispatch {
    ($name:ident, $ty:ty, $body:ident, ($($arg:ident: $t:ty),*) $(-> $ret:ty)?) => {
        /// FMA-codegen dispatch (bead nabk): native fused instruction
        /// when the CPU has avx2+fma, the identical portable body
        /// otherwise. Safe to call; bit-identical either way.
        #[inline]
        pub(crate) fn $name(this: &$ty $(, $arg: $t)*) $(-> $ret)? {
            #[cfg(target_arch = "x86_64")]
            {
                if std::arch::is_x86_feature_detected!("avx2")
                    && std::arch::is_x86_feature_detected!("fma")
                {
                    // SAFETY: avx2+fma verified immediately above; the
                    // inlined body is pure safe slice arithmetic.
                    return unsafe { paste_x86::$name(this $(, $arg)*) };
                }
            }
            this.$body($($arg),*)
        }
    };
}

fma_dispatch!(bsr_spmv_dispatch, Bsr, spmv_body, (x: &[f64], y: &mut [f64]));
fma_dispatch!(sell_spmv_dispatch, Sell, spmv_body, (x: &[f64], y: &mut [f64]));
fma_dispatch!(sell_spmv_chunked_dispatch, Sell, spmv_chunked_body, (x: &[f64], y: &mut [f64]));
fma_dispatch!(sell_shard_dispatch, Sell, shard_body, (x: &[f64], lo: usize, hi: usize) -> Vec<(usize, f64)>);
fma_dispatch!(compact_spmv_dispatch, CsrCompact, spmv_body, (x: &[f64], y: &mut [f64]));
fma_dispatch!(compact_shard_dispatch, CsrCompact, shard_body, (x: &[f64], lo: usize, hi: usize, mine: &mut [f64]));
fma_dispatch!(spgemm_dispatch, Csr, spgemm_body, (b: &Csr) -> Csr);
fma_dispatch!(spgemm_sparse_spa_dispatch, Csr, spgemm_sparse_spa_body, (b: &Csr) -> Csr);
fma_dispatch!(spmm_blocked_dispatch, Csr, spmm_blocked_body, (b: &[f64], nrhs: usize, y: &mut [f64]));

/// The `target_feature`-recompiled bodies (one per dispatcher above).
#[cfg(target_arch = "x86_64")]
mod paste_x86 {
    use super::{Bsr, Csr, CsrCompact, Sell};

    macro_rules! fma_body {
        ($name:ident, $ty:ty, $body:ident, ($($arg:ident: $t:ty),*) $(-> $ret:ty)?) => {
            /// # Safety
            /// Requires avx2+fma, verified by the dispatcher immediately
            /// before the call. The body itself is safe code.
            #[target_feature(enable = "avx2,fma")]
            pub(super) unsafe fn $name(this: &$ty $(, $arg: $t)*) $(-> $ret)? {
                this.$body($($arg),*)
            }
        };
    }

    fma_body!(bsr_spmv_dispatch, Bsr, spmv_body, (x: &[f64], y: &mut [f64]));
    fma_body!(sell_spmv_dispatch, Sell, spmv_body, (x: &[f64], y: &mut [f64]));
    fma_body!(sell_spmv_chunked_dispatch, Sell, spmv_chunked_body, (x: &[f64], y: &mut [f64]));
    fma_body!(sell_shard_dispatch, Sell, shard_body, (x: &[f64], lo: usize, hi: usize) -> Vec<(usize, f64)>);
    fma_body!(compact_spmv_dispatch, CsrCompact, spmv_body, (x: &[f64], y: &mut [f64]));
    fma_body!(compact_shard_dispatch, CsrCompact, shard_body, (x: &[f64], lo: usize, hi: usize, mine: &mut [f64]));
    fma_body!(spgemm_dispatch, Csr, spgemm_body, (b: &Csr) -> Csr);
    fma_body!(spgemm_sparse_spa_dispatch, Csr, spgemm_sparse_spa_body, (b: &Csr) -> Csr);
    fma_body!(spmm_blocked_dispatch, Csr, spmm_blocked_body, (b: &[f64], nrhs: usize, y: &mut [f64]));
}

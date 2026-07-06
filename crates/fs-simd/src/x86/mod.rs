//! x86-64 capsule: AVX2+FMA and AVX-512F variants of the primitive set.
//! Registered unsafe capsule — see SAFETY.md beside this file.
//!
//! Feature-gating contract: the `#[target_feature]` inner functions are
//! reached ONLY through the safe façades below, each of which re-checks the
//! CPU feature at runtime and falls back to the scalar twin otherwise —
//! the façades are therefore unconditionally safe to call, and the dispatch
//! table's tier choice is an optimization, not a soundness precondition.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

// Only the intrinsics the WIRED ops (axpy/dot/sum) use; mul intrinsics
// return here when scale/mul_elem get vector paths (caught by the CI
// both-ISA clippy gate — unused imports never compile on local aarch64).
#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::{
    __m256d, _mm256_add_pd, _mm256_fmadd_pd, _mm256_loadu_pd, _mm256_set1_pd, _mm256_storeu_pd,
    _mm512_add_pd, _mm512_fmadd_pd, _mm512_loadu_pd, _mm512_reduce_add_pd, _mm512_set1_pd,
    _mm512_storeu_pd,
};

/// Horizontal sum of a __m256d, fixed low-to-high lane order.
///
/// # Safety
/// Caller must ensure AVX is available (enforced by the façades' checks).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn hsum256(v: __m256d) -> f64 {
    let mut lanes = [0.0f64; 4];
    // SAFETY: storeu writes exactly 4 f64 into the local array.
    unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), v) };
    ((lanes[0] + lanes[1]) + lanes[2]) + lanes[3]
}

macro_rules! facade {
    ($name:ident, $avx512:ident, $avx2:ident, ($($arg:ident : $ty:ty),*) -> $ret:ty, $scalar:expr) => {
        /// Safe façade: AVX-512 → AVX2+FMA → scalar twin, checked at runtime.
        pub fn $name($($arg: $ty),*) -> $ret {
            #[cfg(target_arch = "x86_64")]
            {
                if std::arch::is_x86_feature_detected!("avx512f") {
                    // SAFETY: feature verified on this CPU immediately above.
                    return unsafe { $avx512($($arg),*) };
                }
                if std::arch::is_x86_feature_detected!("avx2")
                    && std::arch::is_x86_feature_detected!("fma")
                {
                    // SAFETY: features verified on this CPU immediately above.
                    return unsafe { $avx2($($arg),*) };
                }
            }
            #[allow(clippy::redundant_closure_call)]
            ($scalar)($($arg),*)
        }
    };
}

facade!(axpy, axpy_512, axpy_256, (a: f64, x: &[f64], y: &mut [f64]) -> (), crate::scalar::axpy);
facade!(dot, dot_512, dot_256, (x: &[f64], y: &[f64]) -> f64, crate::scalar::dot);
facade!(sum, sum_512, sum_256, (x: &[f64]) -> f64, crate::scalar::sum);

/// y[i] = a*x[i] + y[i], AVX2+FMA lane body, scalar tail.
///
/// # Safety
/// Requires avx2+fma (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn axpy_256(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len(), "axpy length mismatch (programmer error)");
    let (xc, xt) = x.as_chunks::<4>();
    let (yc, yt) = y.as_chunks_mut::<4>();
    // SAFETY: loadu/storeu access exactly 4 f64 at [f64; 4] chunk-array
    // addresses inside live slices; unaligned access is supported.
    unsafe {
        let va = _mm256_set1_pd(a);
        for (xk, yk) in xc.iter().zip(yc) {
            let vx = _mm256_loadu_pd(xk.as_ptr());
            let vy = _mm256_loadu_pd(yk.as_ptr());
            _mm256_storeu_pd(yk.as_mut_ptr(), _mm256_fmadd_pd(va, vx, vy));
        }
    }
    crate::scalar::axpy(a, xt, yt);
}

/// y[i] = a*x[i] + y[i], AVX-512 lane body, scalar tail.
///
/// # Safety
/// Requires avx512f (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn axpy_512(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len(), "axpy length mismatch (programmer error)");
    let (xc, xt) = x.as_chunks::<8>();
    let (yc, yt) = y.as_chunks_mut::<8>();
    // SAFETY: as in axpy_256, 8-lane extents.
    unsafe {
        let va = _mm512_set1_pd(a);
        for (xk, yk) in xc.iter().zip(yc) {
            let vx = _mm512_loadu_pd(xk.as_ptr());
            let vy = _mm512_loadu_pd(yk.as_ptr());
            _mm512_storeu_pd(yk.as_mut_ptr(), _mm512_fmadd_pd(va, vx, vy));
        }
    }
    crate::scalar::axpy(a, xt, yt);
}

/// dot, AVX2+FMA: two 4-lane fused accumulators (even/odd chunks), combined
/// then lane-summed low-to-high, scalar tail appended. Fixed shape per tier.
///
/// # Safety
/// Requires avx2+fma (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn dot_256(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "dot length mismatch");
    let (xc, xt) = x.as_chunks::<8>();
    let (yc, yt) = y.as_chunks::<8>();
    // SAFETY: pointers into [f64; 8] arrays; `.add(4)` stays inside the array.
    let vec_part = unsafe {
        let mut acc0 = _mm256_set1_pd(0.0);
        let mut acc1 = _mm256_set1_pd(0.0);
        for (xk, yk) in xc.iter().zip(yc) {
            acc0 = _mm256_fmadd_pd(
                _mm256_loadu_pd(xk.as_ptr()),
                _mm256_loadu_pd(yk.as_ptr()),
                acc0,
            );
            acc1 = _mm256_fmadd_pd(
                _mm256_loadu_pd(xk.as_ptr().add(4)),
                _mm256_loadu_pd(yk.as_ptr().add(4)),
                acc1,
            );
        }
        hsum256(_mm256_add_pd(acc0, acc1))
    };
    vec_part + crate::scalar::dot(xt, yt)
}

/// dot, AVX-512: two 8-lane fused accumulators, `_mm512_reduce_add_pd`
/// combine, scalar tail. Fixed shape per tier.
///
/// # Safety
/// Requires avx512f (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn dot_512(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "dot length mismatch");
    let (xc, xt) = x.as_chunks::<16>();
    let (yc, yt) = y.as_chunks::<16>();
    // SAFETY: pointers into [f64; 16] arrays; `.add(8)` stays inside.
    let vec_part = unsafe {
        let mut acc0 = _mm512_set1_pd(0.0);
        let mut acc1 = _mm512_set1_pd(0.0);
        for (xk, yk) in xc.iter().zip(yc) {
            acc0 = _mm512_fmadd_pd(
                _mm512_loadu_pd(xk.as_ptr()),
                _mm512_loadu_pd(yk.as_ptr()),
                acc0,
            );
            acc1 = _mm512_fmadd_pd(
                _mm512_loadu_pd(xk.as_ptr().add(8)),
                _mm512_loadu_pd(yk.as_ptr().add(8)),
                acc1,
            );
        }
        _mm512_reduce_add_pd(_mm512_add_pd(acc0, acc1))
    };
    vec_part + crate::scalar::dot(xt, yt)
}

/// sum, AVX2: two-accumulator shape as in dot.
///
/// # Safety
/// Requires avx2+fma (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn sum_256(x: &[f64]) -> f64 {
    let (xc, xt) = x.as_chunks::<8>();
    // SAFETY: pointers into [f64; 8] arrays; `.add(4)` stays inside.
    let vec_part = unsafe {
        let mut acc0 = _mm256_set1_pd(0.0);
        let mut acc1 = _mm256_set1_pd(0.0);
        for xk in xc {
            acc0 = _mm256_add_pd(acc0, _mm256_loadu_pd(xk.as_ptr()));
            acc1 = _mm256_add_pd(acc1, _mm256_loadu_pd(xk.as_ptr().add(4)));
        }
        hsum256(_mm256_add_pd(acc0, acc1))
    };
    vec_part + crate::scalar::sum(xt)
}

/// sum, AVX-512: two-accumulator shape as in dot.
///
/// # Safety
/// Requires avx512f (façade-verified).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn sum_512(x: &[f64]) -> f64 {
    let (xc, xt) = x.as_chunks::<16>();
    // SAFETY: pointers into [f64; 16] arrays; `.add(8)` stays inside.
    let vec_part = unsafe {
        let mut acc0 = _mm512_set1_pd(0.0);
        let mut acc1 = _mm512_set1_pd(0.0);
        for xk in xc {
            acc0 = _mm512_add_pd(acc0, _mm512_loadu_pd(xk.as_ptr()));
            acc1 = _mm512_add_pd(acc1, _mm512_loadu_pd(xk.as_ptr().add(8)));
        }
        _mm512_reduce_add_pd(_mm512_add_pd(acc0, acc1))
    };
    vec_part + crate::scalar::sum(xt)
}

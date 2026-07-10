//! f64 views of C64 slices — the bridge to fs-simd's interleaved-complex
//! stage kernels (bead 27d3). Registered unsafe capsule (SAFETY.md beside
//! this file): two cast helpers, nothing else.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use crate::C64;

/// Reinterpret a complex slice as its interleaved (re, im) f64 bytes.
pub(crate) fn as_f64(v: &[C64]) -> &[f64] {
    // SAFETY: C64 is #[repr(C)] { re: f64, im: f64 } — size 16, align 8,
    // no padding — so `v`'s allocation is exactly `2·len` valid f64s; the
    // shared borrow and lifetime carry over unchanged.
    unsafe { core::slice::from_raw_parts(v.as_ptr().cast::<f64>(), v.len() * 2) }
}

/// Mutable variant of [`as_f64`].
pub(crate) fn as_f64_mut(v: &mut [C64]) -> &mut [f64] {
    // SAFETY: as in `as_f64`; the exclusive borrow carries over, so no
    // aliasing is introduced.
    unsafe { core::slice::from_raw_parts_mut(v.as_mut_ptr().cast::<f64>(), v.len() * 2) }
}

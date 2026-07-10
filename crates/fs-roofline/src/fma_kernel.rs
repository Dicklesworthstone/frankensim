//! The FMA chain kernel behind the roofline compute axis — split from
//! axes.rs under the 300-line capsule cap. The ONLY unsafe here is the
//! x86 `target_feature` body (baseline x86-64 lowers `f64::mul_add` to
//! a libm CALL, which measured 1.0 GFLOP/s on a Threadripper and
//! inflated every attainment 40× — the xlvx finding). Registered in
//! unsafe-capsules.json; SAFETY.md beside this file.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use super::{LANES, STEPS};

fn fma_pass_portable(acc: &mut [f64; LANES], m: f64, a: f64) {
    for _ in 0..STEPS {
        for lane in acc.iter_mut() {
            *lane = lane.mul_add(m, a);
        }
    }
}

/// x86 body with REAL fused-multiply-add codegen (see module docs).
///
/// # Safety
/// Requires avx+fma, verified by the dispatcher immediately before the
/// call. The body is pure safe arithmetic on a caller-owned array.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,fma")]
unsafe fn fma_pass_x86(acc: &mut [f64; LANES], m: f64, a: f64) {
    for _ in 0..STEPS {
        for lane in acc.iter_mut() {
            *lane = lane.mul_add(m, a);
        }
    }
}

/// One chain pass, dispatched to the fastest sound body for this host.
pub(super) fn fma_pass(acc: &mut [f64; LANES], m: f64, a: f64) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx") && std::arch::is_x86_feature_detected!("fma") {
        // SAFETY: features verified on this CPU immediately above.
        return unsafe { fma_pass_x86(acc, m, a) };
    }
    fma_pass_portable(acc, m, a);
}

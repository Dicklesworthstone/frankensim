//! x86 FMA-codegen capsule for the D3Q19 cell collision (pexfp.138, ruling
//! 9653). Baseline x86-64 lowers every `f64::mul_add` to a call of the
//! compiler-builtins `fma()` trampoline. A moment-space cell collision makes
//! about 800 of them (two dense 19 x 19 central-moment transforms, the
//! relaxation and the binomial shift back). This capsule recompiles the SAME
//! `#[inline(always)]` collision body under `#[target_feature(enable =
//! "fma")]`, where `mul_add` is one native fused instruction and no call.
//! IEEE-754 fma is correctly rounded under both codegens and Rust never
//! contracts a separate multiply and add, so the results are BIT-IDENTICAL to
//! the portable path, pinned by `tests.rs` beside this file. Registered
//! capsule; SAFETY.md beside. aarch64 and x86-64 CPUs without FMA run the
//! portable body.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

use super::{CollisionError3, CollisionModel3, Q3, collide_cell3_with_projection};

#[cfg(test)]
mod tests;

/// Test-only switch that sends every collision down the portable path, so a
/// test can run the same grid both ways on one FMA host.
#[cfg(test)]
static PORTABLE_ONLY: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Whether this call runs the FMA-compiled body.
#[cfg(target_arch = "x86_64")]
#[inline]
fn fma_selected() -> bool {
    #[cfg(test)]
    {
        if PORTABLE_ONLY.load(std::sync::atomic::Ordering::SeqCst) {
            return false;
        }
    }
    std::arch::is_x86_feature_detected!("fma")
}

/// Collide one cell with FMA codegen when the CPU has it, portable codegen
/// otherwise. Unconditionally safe to call; the result is bitwise the same.
#[inline]
pub(super) fn collide_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    #[cfg(target_arch = "x86_64")]
    {
        if fma_selected() {
            // SAFETY: fma verified on this CPU immediately above; the inlined
            // body is ordinary safe array arithmetic.
            return unsafe { collide_cell3_fma(populations, model, force) };
        }
    }
    collide_cell3_with_projection(populations, model, force, false)
}

/// The identical collision body compiled with FMA enabled (its
/// `inline(always)` makes the code generate HERE, under this feature).
///
/// # Safety
/// Requires the `fma` feature, verified by the dispatcher immediately before
/// the call. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn collide_cell3_fma(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, force, false)
}

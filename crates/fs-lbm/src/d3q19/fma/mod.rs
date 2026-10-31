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

/// One cell's collision, held by a hot loop for the length of the loop.
pub(crate) type CollisionKernel3 =
    fn([f64; Q3], CollisionModel3, [f64; 3]) -> Result<[f64; Q3], CollisionError3>;

/// The collision kernel for this CPU: FMA codegen when the CPU has it,
/// portable codegen otherwise, bitwise the same. Feature detection runs once
/// per process. Call this before a hot loop, not inside it.
pub(crate) fn collision_kernel3() -> CollisionKernel3 {
    #[cfg(test)]
    {
        if PORTABLE_ONLY.load(std::sync::atomic::Ordering::SeqCst) {
            return portable;
        }
    }
    #[cfg(target_arch = "x86_64")]
    {
        static KERNEL: std::sync::OnceLock<CollisionKernel3> = std::sync::OnceLock::new();
        *KERNEL.get_or_init(|| {
            if std::arch::is_x86_feature_detected!("fma") {
                fma_checked
            } else {
                portable
            }
        })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        portable
    }
}

/// Collide one cell through [`collision_kernel3`], for single calls.
/// Unconditionally safe to call; the result is bitwise the same on any CPU.
#[inline]
pub(super) fn collide_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collision_kernel3()(populations, model, force)
}

fn portable(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, force, false)
}

/// Safe thunk over [`collide_cell3_fma`]. Private, and named only by
/// [`collision_kernel3`], after detection.
#[cfg(target_arch = "x86_64")]
fn fma_checked(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    // SAFETY: `collision_kernel3` is the only code that names this thunk, and
    // it does so only after detecting fma on this CPU. The inlined body is
    // ordinary safe array arithmetic.
    unsafe { collide_cell3_fma(populations, model, force) }
}

/// The identical collision body compiled with FMA enabled (its
/// `inline(always)` makes the code generate HERE, under this feature).
///
/// # Safety
/// Requires the `fma` feature, verified once by [`collision_kernel3`], the
/// only code that admits [`fma_checked`]. The body itself is safe code.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn collide_cell3_fma(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, force, false)
}

//! See SAFETY.md. wasm32 only; an empty crate everywhere else.
#![cfg_attr(not(target_arch = "wasm32"), allow(unused))]

/// Which fused multiply-add is active.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FmaPath {
    /// The CPU's FMA through `f64x2.relaxed_madd`, proven fused by the probe.
    Hardware,
    /// Correctly rounded software FMA (`libm::fma`): the host does not fuse.
    Software,
}

#[cfg(target_arch = "wasm32")]
mod wasm {
    use super::FmaPath;
    use core::arch::wasm32::{f64x2_extract_lane, f64x2_relaxed_madd, f64x2_splat};
    use core::sync::atomic::{AtomicU8, Ordering};

    const UNKNOWN: u8 = 0;
    const HARDWARE: u8 = 1;
    const SOFTWARE: u8 = 2;
    static PATH: AtomicU8 = AtomicU8::new(UNKNOWN);

    #[inline]
    #[target_feature(enable = "relaxed-simd")]
    fn relaxed(a: f64, b: f64, c: f64) -> f64 {
        f64x2_extract_lane::<0>(f64x2_relaxed_madd(f64x2_splat(a), f64x2_splat(b), f64x2_splat(c)))
    }

    #[inline]
    fn hardware(a: f64, b: f64, c: f64) -> f64 {
        // SAFETY: a module containing relaxed-simd opcodes only instantiates
        // on engines that validated them; see SAFETY.md.
        unsafe { relaxed(a, b, c) }
    }

    #[cold]
    fn probe() -> u8 {
        // (1+e)^2 = 1 + 2e + e^2 with e = 2^-30. Rounded first, the product is
        // 1 + 2e and the sum is exactly 0. Fused, the sum is exactly e^2 = 2^-60.
        let e = f64::from_bits(0x3E10_0000_0000_0000); // 2^-30
        let x = core::hint::black_box(1.0 + e);
        let c = core::hint::black_box(-(1.0 + 2.0 * e));
        let fused = hardware(x, x, c) == e * e;
        // A second, sign-and-magnitude-different witness: 3 * (1/3) - 1 is
        // exactly -2^-54 when fused and 0 when the product rounds to 1 first.
        let third = core::hint::black_box(1.0_f64 / 3.0);
        let fused = fused && hardware(core::hint::black_box(3.0), third, -1.0) == -f64::from_bits(0x3C90_0000_0000_0000);
        let path = if fused { HARDWARE } else { SOFTWARE };
        PATH.store(path, Ordering::Relaxed);
        path
    }

    /// The C `fma` the compiler calls for `f64::mul_add` on wasm32.
    #[unsafe(no_mangle)]
    pub extern "C" fn fma(a: f64, b: f64, c: f64) -> f64 {
        let mut path = PATH.load(Ordering::Relaxed);
        if path == UNKNOWN {
            path = probe();
        }
        if path == HARDWARE { hardware(a, b, c) } else { libm::fma(a, b, c) }
    }

    pub fn active() -> FmaPath {
        let mut path = PATH.load(Ordering::Relaxed);
        if path == UNKNOWN {
            path = probe();
        }
        if path == HARDWARE { FmaPath::Hardware } else { FmaPath::Software }
    }

    /// Both implementations on the same operands, for bit-comparison in tests.
    pub fn both(a: f64, b: f64, c: f64) -> (f64, f64) {
        (hardware(a, b, c), libm::fma(a, b, c))
    }
}

/// The active path. Referencing this from the engine also guarantees the
/// object that defines `fma` is linked.
#[cfg(target_arch = "wasm32")]
#[must_use]
pub fn active() -> FmaPath {
    wasm::active()
}

/// (relaxed madd, software fma) of the same operands.
#[cfg(target_arch = "wasm32")]
#[must_use]
pub fn both(a: f64, b: f64, c: f64) -> (f64, f64) {
    wasm::both(a, b, c)
}

//! Error-free transformations (EFT) — the exact building blocks of the
//! precision ladder (plan §6.1/§6.4). Relocated here (L0) from the fs-la
//! mixed-precision bead's original scope so BOTH consumers — fs-ivl's
//! double-double oracle and fs-la's iterative refinement — share one
//! implementation (recorded on beads 6ys.8 / 6ys.12).
//!
//! "Error-free" is literal: `two_sum(a, b)` returns `(s, e)` with
//! `s = fl(a + b)` and `s + e == a + b` EXACTLY (as real numbers), for any
//! finite inputs. These identities are bitwise-testable algebra, not
//! approximations — the G0 suite exploits that.
//!
//! Determinism: pure +,−,× and `mul_add`; cross-ISA bit-deterministic by
//! construction like everything else in this crate.

/// Knuth's TwoSum: `(s, e)` with `s = fl(a+b)` and `s + e = a + b` exactly.
/// No ordering precondition. 6 flops.
#[must_use]
pub fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let e = (a - (s - bb)) + (b - bb);
    (s, e)
}

/// Dekker's FastTwoSum: same contract as [`two_sum`] but REQUIRES
/// `|a| >= |b|` (or a = 0). 3 flops. Debug-asserted, not checked in
/// release — callers on hot paths own the precondition.
#[must_use]
pub fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    debug_assert!(
        a == 0.0 || b == 0.0 || a.abs() >= b.abs() || a.is_nan() || b.is_nan(),
        "quick_two_sum precondition |a| >= |b| violated: a={a}, b={b}"
    );
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

/// TwoProd via FMA: `(p, e)` with `p = fl(a·b)` and `p + e = a·b` exactly
/// (finite, non-overflowing inputs). 2 flops — the FMA IS the splitter
/// (no Dekker splitting needed on FrankenSim's baseline ISAs).
#[must_use]
pub fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let e = a.mul_add(b, -p);
    (p, e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// The EFT identity is EXACT — verify with dd-style reconstruction over
    /// scale-mixed random pairs, plus the classic adversarial cases.
    #[test]
    fn two_sum_identity_is_exact() {
        let mut seed = 0xEF7_u64;
        for _ in 0..100_000 {
            let scale_a = f64::exp2(f64::from((seed % 80) as u32) - 40.0);
            let a = lcg(&mut seed) * scale_a;
            let b = lcg(&mut seed);
            let (s, e) = two_sum(a, b);
            // Exactness check without higher precision: s+e must equal a+b
            // as an UNEVALUATED pair. Reconstruct by subtracting exactly:
            // two_sum(s, e) must return (s, e) itself (already normalized),
            // and s must equal fl(a+b).
            assert_eq!(s.to_bits(), (a + b).to_bits());
            let (s2, e2) = two_sum(s, e);
            assert_eq!(s2.to_bits(), s.to_bits(), "not normalized at a={a}, b={b}");
            assert_eq!(e2.to_bits(), e.to_bits());
            // And the error term is genuinely the rounding error: for cases
            // where the sum is exact, e must be exactly zero.
            let (s3, e3) = two_sum(a, a); // 2a is exact (no rounding)
            assert_eq!(s3.to_bits(), (2.0 * a).to_bits());
            assert_eq!(e3.to_bits(), 0.0f64.to_bits(), "2a is exact; e must be 0");
        }
        // Catastrophic-cancellation classic: (1e16 + 1) - 1e16.
        let (s, e) = two_sum(1e16, 1.0);
        assert_eq!(s.to_bits(), 1e16f64.to_bits(), "1.0 is absorbed");
        assert_eq!(
            e.to_bits(),
            1.0f64.to_bits(),
            "two_sum must recover the absorbed 1.0 exactly"
        );
    }

    #[test]
    fn quick_two_sum_matches_two_sum_when_ordered() {
        let mut seed = 0xFA57_u64;
        for _ in 0..100_000 {
            let x = lcg(&mut seed) * 1e8;
            let y = lcg(&mut seed);
            let (a, b) = if x.abs() >= y.abs() { (x, y) } else { (y, x) };
            let (s1, e1) = two_sum(a, b);
            let (s2, e2) = quick_two_sum(a, b);
            assert_eq!(s1.to_bits(), s2.to_bits());
            assert_eq!(e1.to_bits(), e2.to_bits(), "at a={a}, b={b}");
        }
    }

    #[test]
    fn two_prod_identity_is_exact() {
        let mut seed = 0x920D_u64;
        for _ in 0..100_000 {
            let a = lcg(&mut seed) * 1e5;
            let b = lcg(&mut seed) * 1e-3;
            let (p, e) = two_prod(a, b);
            assert_eq!(p.to_bits(), (a * b).to_bits());
            // fma(a, b, -p) is BY DEFINITION the single-rounding residual;
            // verify the reconstruction is normalized and that exact products
            // give e = 0.
            let (p2, e2) = two_sum(p, e);
            assert_eq!(p2.to_bits(), p.to_bits());
            assert_eq!(e2.to_bits(), e.to_bits());
        }
        // Exact product: powers of two.
        let (p, e) = two_prod(3.0, 0.5);
        assert_eq!(p.to_bits(), 1.5f64.to_bits());
        assert_eq!(e.to_bits(), 0.0f64.to_bits());
        // Inexact product with known residual: (1 + 2^-52)² = 1 + 2^-51 +
        // 2^-104; the tail 2^-104 cannot fit and is exactly the residual.
        let x = 1.0 + f64::EPSILON;
        let (_, e) = two_prod(x, x);
        assert_eq!(
            e.to_bits(),
            (f64::EPSILON * f64::EPSILON).to_bits(),
            "residual of (1+2^-52)^2 is 2^-104"
        );
    }
}

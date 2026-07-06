//! Double-double arithmetic (~106-bit significand) on the EFT primitives —
//! the precision ladder's second rung and fs-ivl's high-precision oracle
//! (plan §6.1/§6.4; relocated to L0 per beads 6ys.8/6ys.12 so fs-la's
//! iterative refinement and fs-ivl share ONE implementation).
//!
//! Representation: `Dd { hi, lo }` with `hi = fl(hi + lo)` and
//! `|lo| ≤ ½ ulp(hi)` (normalized). Documented error bounds (QD-library
//! classical results): add/sub/mul relative error ≤ 2⁻¹⁰⁴, div/sqrt
//! ≤ 2⁻¹⁰³, for finite, non-over/underflowing operands.
//!
//! Determinism: pure +,−,×,÷, sqrt, mul_add — cross-ISA bit-deterministic
//! by construction. Quad-double (the ladder's top rung) is deliberately
//! deferred: dd covers every current oracle need (recorded on 6ys.12).

use crate::eft::{quick_two_sum, two_prod, two_sum};

/// A double-double value: unevaluated sum `hi + lo`, normalized.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Dd {
    /// Leading component (the f64 closest to the represented value).
    pub hi: f64,
    /// Trailing component (the residual; `|lo| ≤ ½ ulp(hi)`).
    pub lo: f64,
}

impl Dd {
    /// Zero.
    pub const ZERO: Dd = Dd { hi: 0.0, lo: 0.0 };
    /// One.
    pub const ONE: Dd = Dd { hi: 1.0, lo: 0.0 };

    /// Lift an f64 (exact).
    #[must_use]
    pub const fn from_f64(x: f64) -> Dd {
        Dd { hi: x, lo: 0.0 }
    }

    /// Construct from an unevaluated pair (renormalizes).
    #[must_use]
    pub fn from_pair(hi: f64, lo: f64) -> Dd {
        let (s, e) = two_sum(hi, lo);
        Dd { hi: s, lo: e }
    }

    /// Round to nearest f64 (the leading component, by the normalization
    /// invariant).
    #[must_use]
    pub const fn to_f64(self) -> f64 {
        self.hi
    }

    /// Negation (exact).
    #[must_use]
    pub const fn neg(self) -> Dd {
        Dd {
            hi: -self.hi,
            lo: -self.lo,
        }
    }

    /// Absolute value (exact).
    #[must_use]
    pub fn abs(self) -> Dd {
        if self.hi < 0.0 || (self.hi == 0.0 && self.lo < 0.0) {
            self.neg()
        } else {
            self
        }
    }

    /// Addition (Knuth accurate variant; relative error ≤ 2⁻¹⁰⁴).
    #[must_use]
    pub fn add(self, o: Dd) -> Dd {
        let (s1, s2) = two_sum(self.hi, o.hi);
        let (t1, t2) = two_sum(self.lo, o.lo);
        let (s1, s2) = quick_two_sum(s1, s2 + t1);
        let (hi, lo) = quick_two_sum(s1, s2 + t2);
        Dd { hi, lo }
    }

    /// Subtraction.
    #[must_use]
    pub fn sub(self, o: Dd) -> Dd {
        self.add(o.neg())
    }

    /// Multiplication (relative error ≤ 2⁻¹⁰⁴).
    #[must_use]
    pub fn mul(self, o: Dd) -> Dd {
        let (p1, p2) = two_prod(self.hi, o.hi);
        let p2 = p2 + (self.hi * o.lo + self.lo * o.hi);
        let (hi, lo) = quick_two_sum(p1, p2);
        Dd { hi, lo }
    }

    /// Division (long-division refinement; relative error ≤ 2⁻¹⁰³).
    #[must_use]
    pub fn div(self, o: Dd) -> Dd {
        let q1 = self.hi / o.hi;
        let r = self.sub(o.mul(Dd::from_f64(q1)));
        let q2 = r.hi / o.hi;
        let r = r.sub(o.mul(Dd::from_f64(q2)));
        let q3 = r.hi / o.hi;
        let (s, e) = quick_two_sum(q1, q2);
        Dd::from_pair(s, e + q3)
    }

    /// Square root (Karp's method; relative error ≤ 2⁻¹⁰³). Negative input
    /// yields NaN components (IEEE convention), zero is exact.
    #[must_use]
    pub fn sqrt(self) -> Dd {
        if self.hi == 0.0 && self.lo == 0.0 {
            return Dd::ZERO;
        }
        let s = self.hi.sqrt(); // IEEE-correctly-rounded (0 ULP)
        let sd = Dd::from_f64(s);
        let e = self.sub(sd.mul(sd)).hi / (2.0 * s);
        Dd::from_pair(s, e)
    }

    /// Compare by value (total on non-NaN).
    #[must_use]
    pub fn lt(self, o: Dd) -> bool {
        self.hi < o.hi || (self.hi == o.hi && self.lo < o.lo)
    }
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

    /// |a − b| as dd, compared against a relative bound.
    fn rel_err_le(a: Dd, b: Dd, bound: f64) -> bool {
        let d = a.sub(b).abs();
        let scale = b.abs().hi.max(f64::MIN_POSITIVE);
        d.hi <= bound * scale
    }

    #[test]
    fn known_values() {
        // 1/3 to dd precision: hi = fl(1/3), lo = the residual; verify by
        // multiplying back: 3 · (1/3) must equal 1 to ≤ 2⁻¹⁰³ relative.
        let third = Dd::ONE.div(Dd::from_f64(3.0));
        let back = third.mul(Dd::from_f64(3.0));
        assert!(rel_err_le(back, Dd::ONE, 1e-30), "3·(1/3) = {back:?}");
        // sqrt(2)² = 2.
        let r2 = Dd::from_f64(2.0).sqrt();
        let two = r2.mul(r2);
        assert!(
            rel_err_le(two, Dd::from_f64(2.0), 1e-30),
            "sqrt(2)^2 = {two:?}"
        );
        // The dd value of sqrt(2) must match the known decimal expansion:
        // 1.41421356237309504880168872420969807856967...
        // hi + lo reconstructed against 1.4142135623730950488016887242097 via
        // a dd computed from two exact halves.
        let known = Dd::from_pair(1.4142135623730951, -9.667293313452913e-17);
        assert!(rel_err_le(r2, known, 1e-31), "sqrt(2) dd = {r2:?}");
    }

    #[test]
    fn addition_recovers_cancellation() {
        // (1e16 + 1) − 1e16 in f64 loses the 1; in dd it survives exactly.
        let a = Dd::from_f64(1e16).add(Dd::ONE);
        let b = a.sub(Dd::from_f64(1e16));
        assert_eq!(b.hi, 1.0);
        assert_eq!(b.lo, 0.0);
    }

    #[test]
    fn normalization_invariant_holds() {
        let mut seed = 0xDD_u64;
        for _ in 0..200_000 {
            let a = Dd::from_f64(lcg(&mut seed) * 1e10).add(Dd::from_f64(lcg(&mut seed)));
            let b = Dd::from_f64(lcg(&mut seed) * 1e-6).add(Dd::from_f64(lcg(&mut seed)));
            for v in [a.add(b), a.sub(b), a.mul(b), a.div(b)] {
                if v.hi.is_finite() && v.hi != 0.0 {
                    assert_eq!(
                        (v.hi + v.lo).to_bits(),
                        v.hi.to_bits(),
                        "not normalized: {v:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn field_properties_to_dd_precision() {
        let mut seed = 0xF1E1D_u64;
        for _ in 0..50_000 {
            let a = Dd::from_f64(lcg(&mut seed) * 100.0).add(Dd::from_f64(lcg(&mut seed) * 1e-14));
            let b = Dd::from_f64(lcg(&mut seed) * 100.0).add(Dd::from_f64(lcg(&mut seed) * 1e-14));
            let c = Dd::from_f64(lcg(&mut seed) * 100.0);
            // a + b − b ≈ a: the round trip loses ~2⁻¹⁰⁵ of the LARGEST
            // intermediate (|a+b| or |b|), not of |a| — near-cancelling
            // pairs (a ≈ −b, or tiny a against b ~ 50) otherwise blow up
            // the relative-to-a error by |b|/|a|.
            let scale = a.abs().hi.max(b.abs().hi).max(f64::MIN_POSITIVE);
            let round_trip = a.add(b).sub(b).sub(a).abs();
            assert!(round_trip.hi <= 1e-29 * scale, "add/sub round trip");
            // Distributivity residual, scaled by the products' magnitude
            // (b + c may itself nearly cancel).
            let lhs = a.mul(b.add(c));
            let rhs = a.mul(b).add(a.mul(c));
            let pscale = (a.abs().hi * (b.abs().hi + c.abs().hi)).max(f64::MIN_POSITIVE);
            assert!(
                lhs.sub(rhs).abs().hi <= 1e-28 * pscale,
                "distributivity residual too big"
            );
            // Division inverse.
            if b.hi.abs() > 1e-3 {
                assert!(rel_err_le(a.div(b).mul(b), a, 1e-29), "a/b*b round trip");
            }
        }
        println!(
            "{{\"suite\":\"fs-math\",\"case\":\"dd-properties\",\"verdict\":\"pass\",\"detail\":\"50k random field-property checks at dd precision\"}}"
        );
    }
}

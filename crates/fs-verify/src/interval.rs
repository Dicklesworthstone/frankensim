//! Rigorous OUTWARD-ROUNDED interval arithmetic (the verifier's
//! foundation): every operation widens its result by one ulp in each
//! direction, so enclosures survive floating-point rounding — no
//! "to-nearest plus slack" caveats. Small on purpose: the verifier
//! needs sums, products, squares, square roots, and division by
//! positive constants; unification with fs-ivl's forms is a CONTRACT
//! no-claim.

/// Nudge one ulp toward −∞.
#[must_use]
pub fn down(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return x;
    }
    let bits = x.to_bits();
    let next = if x > 0.0 {
        bits - 1
    } else if x < 0.0 {
        bits + 1
    } else {
        // ±0 → smallest negative subnormal.
        (1u64 << 63) | 1
    };
    f64::from_bits(next)
}

/// Nudge one ulp toward +∞.
#[must_use]
pub fn up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    let bits = x.to_bits();
    let next = if x > 0.0 {
        bits + 1
    } else if x < 0.0 {
        bits - 1
    } else {
        1u64 // smallest positive subnormal
    };
    f64::from_bits(next)
}

/// A closed interval with outward-rounded arithmetic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Iv {
    /// Lower end.
    pub lo: f64,
    /// Upper end.
    pub hi: f64,
}

impl Iv {
    fn entire() -> Iv {
        Iv {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        }
    }

    fn is_valid_finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }

    /// The point interval.
    #[must_use]
    pub fn point(v: f64) -> Iv {
        Iv { lo: v, hi: v }
    }

    /// Zero.
    #[must_use]
    pub fn zero() -> Iv {
        Iv { lo: 0.0, hi: 0.0 }
    }

    /// True when the interval is non-finite or reversed (the FAIL-CLOSED
    /// trigger).
    #[must_use]
    pub fn is_unbounded(&self) -> bool {
        !self.is_valid_finite()
    }

    /// Outward-rounded sum.
    #[must_use]
    #[allow(clippy::should_implement_trait)] // deliberate: no operator sugar for rigor-bearing ops
    pub fn add(self, o: Iv) -> Iv {
        if !self.is_valid_finite() || !o.is_valid_finite() {
            return Iv::entire();
        }
        Iv {
            lo: down(self.lo + o.lo),
            hi: up(self.hi + o.hi),
        }
    }

    /// Outward-rounded difference.
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, o: Iv) -> Iv {
        if !self.is_valid_finite() || !o.is_valid_finite() {
            return Iv::entire();
        }
        Iv {
            lo: down(self.lo - o.hi),
            hi: up(self.hi - o.lo),
        }
    }

    /// Outward-rounded product.
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn mul(self, o: Iv) -> Iv {
        if !self.is_valid_finite() || !o.is_valid_finite() {
            return Iv::entire();
        }
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Iv {
            lo: down(c.iter().copied().fold(f64::INFINITY, f64::min)),
            hi: up(c.iter().copied().fold(f64::NEG_INFINITY, f64::max)),
        }
    }

    /// Outward-rounded division by a strictly positive interval.
    ///
    /// An invalid, non-finite, or zero-containing divisor returns the entire
    /// real line. Callers must fail closed on that unbounded result.
    #[must_use]
    pub fn div_pos(self, divisor: Iv) -> Iv {
        if !self.lo.is_finite()
            || !self.hi.is_finite()
            || self.lo > self.hi
            || !divisor.lo.is_finite()
            || !divisor.hi.is_finite()
            || divisor.lo <= 0.0
            || divisor.lo > divisor.hi
        {
            return Iv::entire();
        }
        self.mul(Iv {
            lo: down(1.0 / divisor.hi),
            hi: up(1.0 / divisor.lo),
        })
    }

    /// Outward-rounded square (dependency-aware: never negative).
    #[must_use]
    pub fn sq(self) -> Iv {
        let m = self.mul(self);
        Iv {
            lo: m.lo.max(0.0),
            hi: m.hi,
        }
    }

    /// Outward-rounded square root. A wholly negative, reversed, or
    /// non-finite input fails closed; a valid interval crossing zero clamps
    /// its lower endpoint to zero.
    #[must_use]
    pub fn sqrt(self) -> Iv {
        if !self.is_valid_finite() || self.hi < 0.0 {
            return Iv::entire();
        }
        Iv {
            lo: down(self.lo.max(0.0).sqrt()).max(0.0),
            hi: up(self.hi.max(0.0).sqrt()),
        }
    }

    /// Outward-rounded scale by a positive point constant.
    #[must_use]
    pub fn scale_pos(self, s: f64) -> Iv {
        if !self.is_valid_finite() || !s.is_finite() || s <= 0.0 {
            return Iv::entire();
        }
        Iv {
            lo: down(self.lo * s),
            hi: up(self.hi * s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn outward_rounding_contains_reals() {
        // 0.1 + 0.2 is not representable; the enclosure must contain
        // the REAL 0.3 (which lies strictly between the neighbors of
        // the f64 nearest).
        let s = Iv::point(0.1).add(Iv::point(0.2));
        assert!(s.lo < 0.3 && 0.3 < s.hi || (s.lo <= 0.3 && 0.3 <= s.hi));
        // Repeated products stay enclosing.
        let mut p = Iv::point(1.0);
        for _ in 0..40 {
            p = p.mul(Iv::point(1.1));
        }
        // det::powi: const-base f64::powi(40) is exactly the release
        // const-fold divergence case (bead 4xnt).
        let truth = fs_math::det::powi(1.1f64, 40);
        assert!(p.lo <= truth && truth <= p.hi);
    }

    #[test]
    fn nudges_move_strictly() {
        assert!(down(1.0) < 1.0 && up(1.0) > 1.0);
        assert!(down(0.0) < 0.0 && up(0.0) > 0.0);
        assert!(down(-2.5) < -2.5 && up(-2.5) > -2.5);
    }

    /// Direction-aware containment of the EXACT real truth `T = hi + resid`:
    /// `resid > 0` ⇒ `T` exceeds `hi`, so a rigorous enclosure must reach
    /// `up(hi)`; `resid < 0` ⇒ it must reach `down(hi)`. This is what actually
    /// GATES the outward nudge — a plain `lo <= hi <= up_end` passes an
    /// enclosure one ulp too narrow, i.e. an op that DROPPED its widen, which
    /// is the verifier's entire rigor.
    fn contains_truth(iv: Iv, hi: f64, resid: f64) -> bool {
        let hi_ok = if resid > 0.0 {
            iv.hi >= up(hi)
        } else {
            iv.hi >= hi
        };
        let lo_ok = if resid < 0.0 {
            iv.lo <= down(hi)
        } else {
            iv.lo <= hi
        };
        lo_ok && hi_ok
    }

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    #[test]
    fn every_op_rounds_strictly_outward_past_the_true_result() {
        // The two tests above gate only `add`'s down-nudge and the `up`/`down`
        // primitives; dropping the widen from sub/mul/sq/sqrt/scale_pos would
        // keep them green while making the verifier's bounds UNSOUND by a ulp.
        // This battery gates EVERY op in BOTH directions against an exact
        // oracle: add/sub/mul/sq/scale_pos truths are two-sum / two-prod exact,
        // so the residual's sign pins which neighbour the enclosure must reach;
        // sqrt is checked by the exact sign of `root*root - p` (real sqrt(p) is
        // irrational in general and must be strictly straddled).
        use fs_math::dd::Dd;
        let mut seed = 0x0051_A9EF_u64;
        let mut checked = 0u64;
        for _ in 0..20_000 {
            let a = lcg(&mut seed) * 40.0;
            let b = lcg(&mut seed) * 40.0;
            let (da, db) = (Dd::from_f64(a), Dd::from_f64(b));
            for (iv, dd) in [
                (Iv::point(a).add(Iv::point(b)), da + db),
                (Iv::point(a).sub(Iv::point(b)), da - db),
                (Iv::point(a).mul(Iv::point(b)), da * db),
                (Iv::point(a).sq(), da * da),
            ] {
                assert!(
                    contains_truth(iv, dd.hi, dd.lo),
                    "outward rounding lost the truth: {iv:?} vs {dd:?}"
                );
                checked += 1;
            }
            let divisor = b.abs() + 1e-3;
            let quotient = a / divisor;
            // For moderate finite inputs, FMA gives the sign of q*d-a without
            // a second product rounding. That sign places the exact quotient
            // strictly below or above q; unlike Dd division, this is not an
            // approximate oracle being treated as exact.
            let residual = quotient.mul_add(divisor, -a);
            let quotient_iv = Iv::point(a).div_pos(Iv::point(divisor));
            let division_ok = if residual > 0.0 {
                quotient_iv.lo <= down(quotient)
            } else if residual < 0.0 {
                quotient_iv.hi >= up(quotient)
            } else {
                quotient_iv.lo <= quotient && quotient <= quotient_iv.hi
            };
            assert!(division_ok, "positive division lost the truth");
            checked += 1;
            let s = lcg(&mut seed).abs() + 1e-3;
            let sc = da * Dd::from_f64(s);
            assert!(
                contains_truth(Iv::point(a).scale_pos(s), sc.hi, sc.lo),
                "scale_pos lost the truth"
            );
            checked += 1;
            // sqrt: `root = p.sqrt()` is correctly rounded, so real sqrt(p) sits
            // within half a ulp of it and must be strictly straddled. `root² - p`
            // fused (single rounding ⇒ exact sign) is OPPOSITE in sign to
            // `real_sqrt(p) - root`.
            let p = a.abs();
            let root = p.sqrt();
            let d = root.mul_add(root, -p);
            let root_iv = Iv::point(p).sqrt();
            let sqrt_ok = if d > 0.0 {
                root_iv.lo <= down(root) // root overshoots ⇒ truth below root
            } else if d < 0.0 {
                root_iv.hi >= up(root) // root undershoots ⇒ truth above root
            } else {
                root_iv.lo <= root && root <= root_iv.hi // p is an exact square
            };
            assert!(
                sqrt_ok,
                "sqrt did not straddle sqrt({p}): {root_iv:?} root {root}"
            );
            checked += 1;
        }
        assert_eq!(checked, 140_000);
    }

    #[test]
    fn positive_division_refuses_invalid_divisors() {
        for divisor in [
            Iv::point(0.0),
            Iv { lo: -1.0, hi: 1.0 },
            Iv::point(f64::INFINITY),
            Iv {
                lo: f64::NAN,
                hi: 1.0,
            },
        ] {
            assert!(Iv::point(1.0).div_pos(divisor).is_unbounded());
        }
    }

    #[test]
    fn invalid_public_domains_fail_closed_without_panicking() {
        let reversed = Iv { lo: 2.0, hi: 1.0 };
        assert!(reversed.is_unbounded());
        assert!(reversed.add(Iv::point(1.0)).is_unbounded());
        assert!(Iv::point(1.0).sub(reversed).is_unbounded());
        assert!(reversed.mul(Iv::point(2.0)).is_unbounded());
        assert!(reversed.sq().is_unbounded());
        assert!(reversed.sqrt().is_unbounded());
        assert!(Iv { lo: -2.0, hi: -1.0 }.sqrt().is_unbounded());
        assert!(!Iv { lo: -1.0, hi: 4.0 }.sqrt().is_unbounded());

        for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Iv::point(1.0).scale_pos(scale).is_unbounded());
        }
        assert!(reversed.scale_pos(2.0).is_unbounded());
    }
}

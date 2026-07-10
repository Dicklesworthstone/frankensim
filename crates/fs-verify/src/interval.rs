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

    /// True when either end is non-finite (the FAIL-CLOSED trigger).
    #[must_use]
    pub fn is_unbounded(&self) -> bool {
        !self.lo.is_finite() || !self.hi.is_finite()
    }

    /// Outward-rounded sum.
    #[must_use]
    #[allow(clippy::should_implement_trait)] // deliberate: no operator sugar for rigor-bearing ops
    pub fn add(self, o: Iv) -> Iv {
        Iv {
            lo: down(self.lo + o.lo),
            hi: up(self.hi + o.hi),
        }
    }

    /// Outward-rounded difference.
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn sub(self, o: Iv) -> Iv {
        Iv {
            lo: down(self.lo - o.hi),
            hi: up(self.hi - o.lo),
        }
    }

    /// Outward-rounded product.
    #[must_use]
    #[allow(clippy::should_implement_trait)]
    pub fn mul(self, o: Iv) -> Iv {
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

    /// Outward-rounded square (dependency-aware: never negative).
    #[must_use]
    pub fn sq(self) -> Iv {
        let m = self.mul(self);
        Iv {
            lo: m.lo.max(0.0),
            hi: m.hi,
        }
    }

    /// Outward-rounded square root (requires `lo ≥ 0`; clamps tiny
    /// negative rounding residue at zero).
    #[must_use]
    pub fn sqrt(self) -> Iv {
        Iv {
            lo: down(self.lo.max(0.0).sqrt()).max(0.0),
            hi: up(self.hi.max(0.0).sqrt()),
        }
    }

    /// Outward-rounded scale by a positive point constant.
    #[must_use]
    pub fn scale_pos(self, s: f64) -> Iv {
        debug_assert!(s > 0.0);
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
}

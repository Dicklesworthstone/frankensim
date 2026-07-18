//! Outward-rounded interval arithmetic (plan §6.4). The postcondition of
//! EVERY operation is ENCLOSURE: the result interval contains the true
//! image of the inputs. That single law is what the word "certified" means
//! everywhere else in FrankenSim, and it is the G0 property the test
//! battery hammers.
//!
//! Rounding strategy: no global rounding-mode changes (Rust-safe,
//! thread-safe, SIMD-mixing-safe). Basic ops are IEEE correctly rounded
//! (error ≤ ½ ulp), so nudging the computed endpoint ONE step outward is
//! rigorous. Elementary functions come from fs-math's det module with
//! DECLARED ULP budgets (see fs-math CONTRACT.md); endpoints are nudged
//! outward by exactly that budget. If fs-math's budgets hold (empirically
//! enforced by its 200k-sample battery), these enclosures are rigorous.

use fs_math::det;

/// A closed interval [lo, hi] of f64 (±∞ allowed as endpoints; NaN never —
/// constructors reject it). INVARIANT: lo ≤ hi.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

/// Nudge a FINITE value k ULP steps down; ±∞ passes through (an infinite
/// bound is already maximally pessimistic).
fn down_k(mut x: f64, k: u64) -> f64 {
    if x.is_finite() {
        for _ in 0..k {
            x = fs_math::next_down(x);
        }
    }
    x
}

/// Nudge k steps up; ±∞ passes through.
fn up_k(mut x: f64, k: u64) -> f64 {
    if x.is_finite() {
        for _ in 0..k {
            x = fs_math::next_up(x);
        }
    }
    x
}

/// Enclose the exact real result represented by one rounded binary endpoint
/// operation.
///
/// A finite IEEE operation normally needs only the usual one-ULP outward
/// nudge. The exceptional case is overflow from two finite operands: a rounded
/// `+inf` stands for some finite real value above `f64::MAX`, and a rounded
/// `-inf` stands for one below `-f64::MAX`. Treating either as a singleton
/// infinity would exclude that finite truth. Actual infinity operands retain
/// their extended-real singleton semantics.
fn enclose_rounded_binary(value: f64, finite_operands: bool) -> Option<Interval> {
    if value.is_nan() {
        return None;
    }
    if finite_operands && value == f64::INFINITY {
        return Some(Interval {
            lo: f64::MAX,
            hi: f64::INFINITY,
        });
    }
    if finite_operands && value == f64::NEG_INFINITY {
        return Some(Interval {
            lo: f64::NEG_INFINITY,
            hi: -f64::MAX,
        });
    }
    if value.is_finite() {
        return Some(Interval {
            lo: down_k(value, 1),
            hi: up_k(value, 1),
        });
    }
    Some(Interval {
        lo: value,
        hi: value,
    })
}

impl Interval {
    /// The whole extended real line — the "no useful enclosure" answer.
    pub const WHOLE: Interval = Interval {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    /// Construct from endpoints. Panics (structured) on NaN or inverted
    /// endpoints — an invalid interval silently propagating would void
    /// every certificate downstream.
    #[must_use]
    pub fn new(lo: f64, hi: f64) -> Interval {
        assert!(
            !lo.is_nan() && !hi.is_nan(),
            "interval endpoints must not be NaN"
        );
        assert!(lo <= hi, "inverted interval [{lo}, {hi}]");
        Interval { lo, hi }
    }

    /// The degenerate interval [x, x] (exact — a point is its own enclosure).
    #[must_use]
    pub fn point(x: f64) -> Interval {
        Interval::new(x, x)
    }

    /// Lower endpoint.
    #[must_use]
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// Upper endpoint.
    #[must_use]
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// Width (hi − lo, rounded up — a certified width never understates).
    #[must_use]
    pub fn width(self) -> f64 {
        up_k(self.hi - self.lo, 1)
    }

    /// Midpoint (a representative point, NOT certified to be central).
    #[must_use]
    pub fn midpoint(self) -> f64 {
        let m = f64::midpoint(self.lo, self.hi);
        if m.is_finite() { m } else { 0.0 }
    }

    /// Does the interval contain `x`?
    #[must_use]
    pub fn contains(self, x: f64) -> bool {
        self.lo <= x && x <= self.hi
    }

    /// Does the interval contain zero (division hazard predicate)?
    #[must_use]
    pub fn contains_zero(self) -> bool {
        self.contains(0.0)
    }

    /// Is `self` a superset of `other`? (The containment-law comparator.)
    #[must_use]
    pub fn encloses(self, other: Interval) -> bool {
        self.lo <= other.lo && other.hi <= self.hi
    }

    /// Convex hull of two intervals.
    #[must_use]
    pub fn hull(self, o: Interval) -> Interval {
        Interval {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    /// Intersection; `None` when disjoint.
    #[must_use]
    pub fn intersect(self, o: Interval) -> Option<Interval> {
        let lo = self.lo.max(o.lo);
        let hi = self.hi.min(o.hi);
        (lo <= hi).then_some(Interval { lo, hi })
    }

    /// Absolute value (exact endpoint selection).
    #[must_use]
    pub fn abs(self) -> Interval {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            -self
        } else {
            Interval {
                lo: 0.0,
                hi: self.hi.max(-self.lo),
            }
        }
    }

    /// Square root. Requires hi ≥ 0; the sub-zero part of a straddling
    /// interval is clipped (sqrt's real domain), lo < 0 ≤ hi is legal.
    /// Panics (structured) if the entire interval is negative.
    #[must_use]
    pub fn sqrt(self) -> Interval {
        assert!(
            self.hi >= 0.0,
            "sqrt of entirely negative interval [{}, {}]",
            self.lo,
            self.hi
        );
        let lo = if self.lo <= 0.0 {
            0.0
        } else {
            down_k(self.lo.sqrt(), 1).max(0.0)
        };
        Interval {
            lo,
            hi: up_k(self.hi.sqrt(), 1),
        }
    }

    /// e^x (monotone; endpoints nudged by the declared budget). The exact
    /// lower bound 0 is preserved (exp > 0: never nudged below zero).
    #[must_use]
    pub fn exp(self) -> Interval {
        let elo = det::exp(self.lo);
        // Overflow: `exp` of a FINITE argument that rounds to +∞ is a finite
        // real above `f64::MAX`, so the largest finite value is a valid lower
        // bound — leaving `lo = +∞` would make [∞, ∞] EXCLUDE every true value
        // (bead wa8i V2). The upper bound may legitimately be +∞.
        let lo = if elo.is_infinite() && self.lo.is_finite() {
            f64::MAX
        } else {
            down_k(elo, det::EXP_ULP_BUDGET).max(0.0)
        };
        Interval {
            lo,
            hi: up_k(det::exp(self.hi), det::EXP_ULP_BUDGET),
        }
    }

    /// ln(x). Requires hi > 0; a zero/negative-straddling interval gets
    /// lower bound −∞ (enclosure of ln over the positive part). Panics if
    /// hi ≤ 0.
    #[must_use]
    pub fn ln(self) -> Interval {
        assert!(
            self.hi > 0.0,
            "ln of non-positive interval [{}, {}]",
            self.lo,
            self.hi
        );
        let lo = if self.lo <= 0.0 {
            f64::NEG_INFINITY
        } else {
            down_k(det::ln(self.lo), det::LN_ULP_BUDGET)
        };
        Interval {
            lo,
            hi: up_k(det::ln(self.hi), det::LN_ULP_BUDGET),
        }
    }

    /// tanh (monotone, odd; budget 5).
    #[must_use]
    pub fn tanh(self) -> Interval {
        Interval {
            lo: down_k(det::tanh(self.lo), det::TANH_ULP_BUDGET).max(-1.0),
            hi: up_k(det::tanh(self.hi), det::TANH_ULP_BUDGET).min(1.0),
        }
    }

    /// sin over the interval, rigorous under fs-math's declared trig budget
    /// (valid within |x| ≤ TRIG_DOMAIN; wider inputs → [-1, 1]).
    #[must_use]
    pub fn sin(self) -> Interval {
        trig_enclosure(self, det::sin, det::SIN_ULP_BUDGET, 0.5)
    }

    /// cos over the interval (critical points at kπ).
    #[must_use]
    pub fn cos(self) -> Interval {
        trig_enclosure(self, det::cos, det::COS_ULP_BUDGET, 0.0)
    }
}

impl core::ops::Neg for Interval {
    type Output = Interval;
    /// Negation (exact: IEEE negation is sign-bit flip).
    fn neg(self) -> Interval {
        Interval {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

impl core::ops::Add for Interval {
    type Output = Interval;
    /// Addition, outward-rounded.
    fn add(self, o: Interval) -> Interval {
        let lo = enclose_rounded_binary(self.lo + o.lo, self.lo.is_finite() && o.lo.is_finite());
        let hi = enclose_rounded_binary(self.hi + o.hi, self.hi.is_finite() && o.hi.is_finite());
        // ∞ + (−∞) = NaN; a NaN endpoint would bypass `new`'s guard and give a
        // non-enclosing interval — fall back to WHOLE (bead wa8i V3, mirroring
        // the Mul/Div guard).
        match (lo, hi) {
            (Some(lo), Some(hi)) => Interval {
                lo: lo.lo,
                hi: hi.hi,
            },
            _ => Interval::WHOLE,
        }
    }
}

impl core::ops::Sub for Interval {
    type Output = Interval;
    /// Subtraction, outward-rounded.
    fn sub(self, o: Interval) -> Interval {
        let lo = enclose_rounded_binary(self.lo - o.hi, self.lo.is_finite() && o.hi.is_finite());
        let hi = enclose_rounded_binary(self.hi - o.lo, self.hi.is_finite() && o.lo.is_finite());
        // ∞ − ∞ = NaN → WHOLE (bead wa8i V3), as in Add/Mul/Div.
        match (lo, hi) {
            (Some(lo), Some(hi)) => Interval {
                lo: lo.lo,
                hi: hi.hi,
            },
            _ => Interval::WHOLE,
        }
    }
}

impl core::ops::Mul for Interval {
    type Output = Interval;
    /// Multiplication, outward-rounded (min/max over endpoint products).
    /// 0·∞ ambiguities (NaN products) fall back to [`Interval::WHOLE`] —
    /// conservative, never wrong.
    fn mul(self, o: Interval) -> Interval {
        let pairs = [
            (self.lo, o.lo),
            (self.lo, o.hi),
            (self.hi, o.lo),
            (self.hi, o.hi),
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (left, right) in pairs {
            let Some(product) =
                enclose_rounded_binary(left * right, left.is_finite() && right.is_finite())
            else {
                return Interval::WHOLE;
            };
            lo = lo.min(product.lo);
            hi = hi.max(product.hi);
        }
        Interval { lo, hi }
    }
}

impl core::ops::Div for Interval {
    type Output = Interval;
    /// Division, outward-rounded. A zero-containing divisor yields
    /// [`Interval::WHOLE`] (documented: "no useful enclosure" is an answer,
    /// not a panic — certified callers branch on it).
    fn div(self, o: Interval) -> Interval {
        if o.contains_zero() {
            return Interval::WHOLE;
        }
        let pairs = [
            (self.lo, o.lo),
            (self.lo, o.hi),
            (self.hi, o.lo),
            (self.hi, o.hi),
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (numerator, denominator) in pairs {
            let Some(quotient) = enclose_rounded_binary(
                numerator / denominator,
                numerator.is_finite() && denominator.is_finite(),
            ) else {
                return Interval::WHOLE;
            };
            lo = lo.min(quotient.lo);
            hi = hi.max(quotient.hi);
        }
        Interval { lo, hi }
    }
}

/// Shared sin/cos enclosure. Critical points sit at (k + offset)·π where
/// extremum sign alternates with k parity: +1 at even k, −1 at odd k
/// (offset 0.5 → sin's π/2 + kπ; offset 0 → cos's kπ). Candidate critical
/// points are tested against a CONSERVATIVE interval enclosure of π
/// (machine π is correctly rounded, so real π ∈ [next_down(π), next_up(π)]);
/// possible containment counts as containment — over-approximation is
/// rigorous, under-approximation would be a lie.
fn trig_enclosure(x: Interval, f: fn(f64) -> f64, budget: u64, offset: f64) -> Interval {
    const UNIT: Interval = Interval { lo: -1.0, hi: 1.0 };
    let (a, b) = (x.lo(), x.hi());
    // Outside the declared trig domain (or unbounded): the honest answer.
    if !a.is_finite() || !b.is_finite() || a.abs() > det::TRIG_DOMAIN || b.abs() > det::TRIG_DOMAIN
    {
        return UNIT;
    }
    // Width ≥ 2π certainly covers a full period.
    if b - a >= 6.3 {
        return UNIT;
    }
    let pi_lo = fs_math::next_down(std::f64::consts::PI);
    let pi_hi = fs_math::next_up(std::f64::consts::PI);
    // Endpoint images, budget-nudged (clamped to the range of sin/cos —
    // nudging must never claim values outside [-1, 1]).
    let fa_lo = down_k(f(a), budget).clamp(-1.0, 1.0);
    let fa_hi = up_k(f(a), budget).clamp(-1.0, 1.0);
    let fb_lo = down_k(f(b), budget).clamp(-1.0, 1.0);
    let fb_hi = up_k(f(b), budget).clamp(-1.0, 1.0);
    let mut out = Interval::new(fa_lo.min(fb_lo), fa_hi.max(fb_hi));
    // Candidate k window: (a/π − offset) − 2 .. (b/π − offset) + 2 is
    // generously conservative for the ≤ 2π-wide interval.
    let k_min = (a / pi_hi - offset).floor() as i64 - 2;
    let k_max = (b / pi_lo - offset).ceil() as i64 + 2;
    for k in k_min..=k_max {
        let m = k as f64 + offset; // exact for |k| ≤ 2^51
        // Conservative enclosure of m·π.
        let (c_lo, c_hi) = if m >= 0.0 {
            (down_k(m * pi_lo, 1), up_k(m * pi_hi, 1))
        } else {
            (down_k(m * pi_hi, 1), up_k(m * pi_lo, 1))
        };
        // Might the critical point lie inside [a, b]?
        if c_hi >= a && c_lo <= b {
            let extremum = if k.rem_euclid(2) == 0 { 1.0 } else { -1.0 };
            out = out.hull(Interval::point(extremum));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::dd::Dd;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// RIGOROUS containment of the EXACT real truth `T = hi + resid` by `iv`,
    /// used to gate the ±1 ulp outward nudge that makes every enclosure sound.
    ///
    /// For `+`, `−`, `×` of f64-exact operands the dd oracle is exact
    /// (two-sum / two-prod), so `resid`'s SIGN pins which side of `hi` the truth
    /// lies on: if `resid > 0` the truth exceeds `hi` and the enclosure must
    /// reach `next_up(hi)` (symmetrically below). A plain `iv.contains(hi)` — or
    /// the old `|| iv.contains(hi + 2·resid)` grace — only demands the endpoint
    /// reach `hi`, so it PASSES an enclosure that is one ulp too narrow: exactly
    /// the enclosure a DROPPED nudge produces, which excludes the true value and
    /// is unsound. That blind spot fires on 18,959 of this battery's 300k checks
    /// (measured), so the grace form gated nothing. dd division carries only a
    /// ~2⁻¹⁰⁴ error, far below a ulp, so the same directed bound applies.
    fn contains_truth(iv: Interval, hi: f64, resid: f64) -> bool {
        let hi_ok = if resid > 0.0 {
            iv.hi() >= fs_math::next_up(hi)
        } else {
            iv.hi() >= hi
        };
        let lo_ok = if resid < 0.0 {
            iv.lo() <= fs_math::next_down(hi)
        } else {
            iv.lo() <= hi
        };
        lo_ok && hi_ok
    }

    #[test]
    fn arithmetic_containment_vs_dd_oracle() {
        // The G0 law, gated RIGOROUSLY (see `contains_truth`): for random
        // interval pairs and random interior points, the interval result must
        // contain the EXACT point truth — reaching past the nearest float when
        // the true value falls between machine numbers, which is what actually
        // exercises the outward nudge.
        let mut seed = 0x1712_u64;
        let mut cases = 0u64;
        for _ in 0..20_000 {
            let c1 = lcg(&mut seed) * 100.0;
            let w1 = (lcg(&mut seed) + 0.5).abs() * 3.0;
            let c2 = lcg(&mut seed) * 100.0;
            let w2 = (lcg(&mut seed) + 0.5).abs() * 3.0;
            let x = Interval::new(c1 - w1, c1 + w1);
            let y = Interval::new(c2 - w2, c2 + w2);
            for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                // Clamp: the affine sample can round OUTSIDE the box at the
                // ends (t = 1), and the containment law only speaks for
                // points inside the inputs.
                let px = (x.lo() + t * (x.hi() - x.lo())).clamp(x.lo(), x.hi());
                let py = (y.lo() + (1.0 - t) * (y.hi() - y.lo())).clamp(y.lo(), y.hi());
                let (dx, dy) = (Dd::from_f64(px), Dd::from_f64(py));
                for (iv, dd) in [(x + y, dx + dy), (x - y, dx - dy), (x * y, dx * dy)] {
                    assert!(
                        contains_truth(iv, dd.hi, dd.lo),
                        "containment violated: {iv:?} excludes truth {}+{} (oracle {dd:?})",
                        dd.hi,
                        dd.lo
                    );
                    cases += 1;
                }
                if !y.contains_zero() {
                    let (iv, dd) = (x / y, dx / dy);
                    assert!(
                        contains_truth(iv, dd.hi, dd.lo),
                        "div containment: {iv:?} excludes truth (oracle {dd:?})"
                    );
                    cases += 1;
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-ivl\",\"case\":\"containment-arith\",\"verdict\":\"pass\",\"detail\":\"{cases} point-oracle checks, 0 violations\"}}"
        );
    }

    #[test]
    fn rewrite_pairs_both_contain_truth() {
        // x·(y+z) and x·y + x·z are the same real function; BOTH interval
        // evaluations must contain the true value (they may differ in
        // tightness — the dependency problem — but never in truth).
        let mut seed = 0x2E11_u64;
        for _ in 0..20_000 {
            let x = Interval::new(-2.0 + lcg(&mut seed), -1.0 + lcg(&mut seed) + 3.0);
            let y = Interval::new(lcg(&mut seed) * 4.0 - 2.5, lcg(&mut seed).abs() * 3.0 + 2.0);
            let z = Interval::new(lcg(&mut seed) - 1.5, lcg(&mut seed) + 1.5);
            let form_a = x * (y + z);
            let form_b = x * y + x * z;
            let (px, py, pz) = (x.midpoint(), y.midpoint(), z.midpoint());
            let truth = Dd::from_f64(px) * (Dd::from_f64(py) + Dd::from_f64(pz));
            assert!(form_a.contains(truth.hi), "factored form lost truth");
            assert!(form_b.contains(truth.hi), "distributed form lost truth");
            // Sub-distributivity: the factored form is never wider than the
            // distributed one on these shapes... NOT asserted (not a law);
            // what IS a law: they intersect (both contain truth).
            assert!(
                form_a.intersect(form_b).is_some(),
                "rewrites disjoint — impossible"
            );
        }
    }

    #[test]
    fn elementary_functions_contain_platform_oracle() {
        // Platform libm (different implementation) as an independent point
        // oracle: its value is within a few ULP of truth, so the enclosure
        // must contain it after accounting a 1-ulp oracle slack — in
        // practice enclosures are ≥ 3 ulp wide, so direct containment.
        let mut seed = 0xE1E_u64;
        for _ in 0..50_000 {
            let c = lcg(&mut seed) * 20.0;
            let w = (lcg(&mut seed) + 0.5).abs() * 2.0;
            let x = Interval::new(c - w, c + w);
            let p = x.lo() + (lcg(&mut seed) + 0.5).clamp(0.0, 1.0) * (x.hi() - x.lo());
            assert!(x.exp().contains(p.exp()), "exp lost {p}");
            assert!(x.sin().contains(p.sin()), "sin lost {p} in {x:?}");
            assert!(x.cos().contains(p.cos()), "cos lost {p} in {x:?}");
            assert!(x.tanh().contains(p.tanh()), "tanh lost {p}");
            if x.lo() > 0.0 {
                assert!(x.ln().contains(p.ln()), "ln lost {p}");
                assert!(x.sqrt().contains(p.sqrt()), "sqrt lost {p}");
            }
        }
        println!(
            "{{\"suite\":\"fs-ivl\",\"case\":\"containment-elem\",\"verdict\":\"pass\",\"detail\":\"50k intervals x 6 functions vs libm point oracle\"}}"
        );
    }

    #[test]
    fn trig_extrema_are_captured() {
        // An interval straddling π/2 must have sin upper bound exactly 1.
        let x = Interval::new(1.4, 1.8);
        let s = x.sin();
        assert!(
            (s.hi() - 1.0).abs() < 1e-15 && s.hi() >= 1.0,
            "sin peak missed: {s:?}"
        );
        // Straddling π: cos lower bound −1.
        let y = Interval::new(3.0, 3.3);
        let c = y.cos();
        assert!(c.lo() <= -1.0 + 1e-15, "cos trough missed: {c:?}");
        // Width ≥ 2π: the unit interval.
        let wide = Interval::new(0.0, 7.0);
        assert_eq!(wide.sin(), Interval::new(-1.0, 1.0));
        // Narrow monotone stretch must NOT include ±1.
        let narrow = Interval::new(0.1, 0.2);
        let sn = narrow.sin();
        assert!(
            sn.hi() < 0.21 && sn.lo() > 0.09,
            "narrow sin too loose: {sn:?}"
        );
    }

    #[test]
    fn division_by_zero_straddling_interval_is_whole_line() {
        let x = Interval::new(1.0, 2.0);
        let y = Interval::new(-0.5, 0.5);
        assert_eq!(x / y, Interval::WHOLE);
        // But a signed non-zero divisor works.
        let z = Interval::new(0.25, 0.5);
        let q = x / z;
        assert!(q.contains(4.0) && q.contains(8.0) && q.lo() >= 1.9);
    }

    #[test]
    fn exp_overflow_lower_bound_stays_finite_and_encloses() {
        // exp of a huge-but-FINITE argument overflows f64 to +∞; the LOWER
        // bound must stay finite (f64::MAX ≤ every true value, which are all
        // above f64::MAX) so the interval still encloses them — bead wa8i V2.
        let r = Interval::new(710.0, 720.0).exp();
        assert!(
            r.lo().is_finite() && r.lo() > 1e308,
            "exp overflow lower bound must be finite (f64::MAX), got {}",
            r.lo()
        );
        assert!(r.hi().is_infinite() && r.hi() > 0.0);
        // A non-overflowing exp still tightly encloses the truth.
        let e = Interval::new(0.0, 1.0).exp();
        assert!(e.lo() <= 1.0 && e.hi() >= std::f64::consts::E && e.lo() >= 0.0);
    }

    #[test]
    fn add_sub_opposing_infinities_are_whole_not_nan() {
        // ∞ + (−∞) and ∞ − ∞ yield NaN endpoints; the ops must fall back to
        // WHOLE rather than emit a NaN interval that encloses nothing (and
        // would bypass `new`'s NaN guard) — bead wa8i V3.
        let pinf = Interval::new(f64::INFINITY, f64::INFINITY);
        let neg = Interval::new(f64::NEG_INFINITY, 0.0);
        assert_eq!(pinf + neg, Interval::WHOLE);
        let wide = Interval::new(0.0, f64::INFINITY);
        assert_eq!(pinf - wide, Interval::WHOLE);
        // Finite arithmetic is unaffected (outward-rounded enclosure).
        let s = Interval::new(1.0, 2.0) + Interval::new(3.0, 4.0);
        assert!(s.contains(4.0) && s.contains(6.0) && s.lo() <= 4.0 && s.hi() >= 6.0);
    }

    #[test]
    fn finite_basic_operation_overflow_is_one_sided_not_singleton_infinity() {
        let max = Interval::point(f64::MAX);
        let neg_max = Interval::point(-f64::MAX);
        let two = Interval::point(2.0);
        let min_positive = Interval::point(f64::MIN_POSITIVE);

        for positive in [max + max, max - neg_max, max * two, max / min_positive] {
            assert_eq!(
                positive.lo().to_bits(),
                f64::MAX.to_bits(),
                "finite positive overflow lost its finite lower bound: {positive:?}"
            );
            assert_eq!(positive.hi(), f64::INFINITY);
        }

        for negative in [
            neg_max + neg_max,
            neg_max - max,
            neg_max * two,
            neg_max / min_positive,
        ] {
            assert_eq!(negative.lo(), f64::NEG_INFINITY);
            assert_eq!(
                negative.hi().to_bits(),
                (-f64::MAX).to_bits(),
                "finite negative overflow lost its finite upper bound: {negative:?}"
            );
        }

        // If both signs can overflow, the honest hull is the whole line.
        assert_eq!(Interval::new(-f64::MAX, f64::MAX) * two, Interval::WHOLE);
    }

    #[test]
    fn actual_extended_infinities_remain_distinct_from_finite_overflow() {
        let pinf = Interval::point(f64::INFINITY);
        let ninf = Interval::point(f64::NEG_INFINITY);
        let two = Interval::point(2.0);

        assert_eq!(pinf + Interval::point(1.0), pinf);
        assert_eq!(pinf * two, pinf);
        assert_eq!(ninf - Interval::point(1.0), ninf);
        assert_eq!(ninf * two, ninf);
        assert_eq!(Interval::point(0.0) * pinf, Interval::WHOLE);
    }

    #[test]
    fn mixed_finite_and_infinite_endpoint_hulls_keep_every_limit_value() {
        let positive_tail = Interval::new(f64::MAX, f64::INFINITY);
        let negative_tail = Interval::new(f64::NEG_INFINITY, -f64::MAX);
        let one = Interval::point(1.0);

        let shifted_positive = positive_tail + one;
        assert_eq!(shifted_positive.lo().to_bits(), f64::MAX.to_bits());
        assert_eq!(shifted_positive.hi(), f64::INFINITY);

        let shifted_negative = negative_tail - one;
        assert_eq!(shifted_negative.lo(), f64::NEG_INFINITY);
        assert_eq!(shifted_negative.hi().to_bits(), (-f64::MAX).to_bits());

        // The quotient image approaches zero at the infinite endpoint and
        // reaches 2/MAX at the finite endpoint.  Outward rounding may include
        // one negative subnormal around the exact zero limit, but it must not
        // omit either limiting value.
        let reciprocal_tail = Interval::new(1.0, 2.0) / positive_tail;
        assert!(reciprocal_tail.contains(0.0));
        assert!(reciprocal_tail.lo() <= 0.0);
        assert!(reciprocal_tail.hi() >= 2.0 / f64::MAX);
    }

    #[test]
    fn underflow_and_extreme_cancellation_still_enclose_the_real_result() {
        let min_subnormal = f64::from_bits(1);
        let underflow = Interval::point(min_subnormal) * Interval::point(0.5);
        assert!(
            underflow.lo() < 0.0 && underflow.hi() > 0.0,
            "half a minimum subnormal must lie strictly inside {underflow:?}"
        );

        let cancellation = Interval::point(f64::MAX) + Interval::point(-f64::MAX);
        assert!(cancellation.contains(0.0));
        assert!(cancellation.lo() < 0.0 && cancellation.hi() > 0.0);
    }

    #[test]
    fn rumps_polynomial_contains_exact_value() {
        // Rump (1988): f(77617, 33096) with
        // f = 333.75 y⁶ + x²(11x²y² − y⁶ − 121y⁴ − 2) + 5.5 y⁸ + x/(2y).
        // Naive f64 evaluation confidently reports ~1.18e21; the true value
        // is −54767/66192 ≈ −0.8274. The certified evaluation may be WIDE —
        // that's the honest outcome under this cancellation — but it must
        // CONTAIN the truth, which is exactly the lesson: wide-and-right
        // beats narrow-and-wrong.
        let x = Interval::point(77617.0);
        let y = Interval::point(33096.0);
        let y2 = y * y;
        let y4 = y2 * y2;
        let y6 = y4 * y2;
        let y8 = y4 * y4;
        let x2 = x * x;
        let t1 = Interval::point(333.75) * y6;
        let inner = Interval::point(11.0) * x2 * y2
            - y6
            - Interval::point(121.0) * y4
            - Interval::point(2.0);
        let t2 = x2 * inner;
        let t3 = Interval::point(5.5) * y8;
        let t4 = x / (Interval::point(2.0) * y);
        let f = t1 + t2 + t3 + t4;
        let truth = -54767.0 / 66192.0;
        assert!(
            f.contains(truth),
            "Rump containment violated: {f:?} vs {truth}"
        );
        // And the naive f64 answer is NOT presented as certain: the interval
        // is honestly wide (width >> 1).
        assert!(
            f.width() > 1.0,
            "suspiciously narrow for this cancellation: {f:?}"
        );
        println!(
            "{{\"suite\":\"fs-ivl\",\"case\":\"rump\",\"verdict\":\"pass\",\"detail\":\"enclosure [{:.3e}, {:.3e}] contains -54767/66192\"}}",
            f.lo(),
            f.hi()
        );
    }

    #[test]
    fn set_operations_and_constructors() {
        let a = Interval::new(1.0, 3.0);
        let b = Interval::new(2.0, 5.0);
        assert_eq!(a.hull(b), Interval::new(1.0, 5.0));
        assert_eq!(a.intersect(b), Some(Interval::new(2.0, 3.0)));
        assert_eq!(a.intersect(Interval::new(4.0, 5.0)), None);
        assert!(Interval::new(0.0, 4.0).encloses(a));
        assert!(!a.encloses(b));
        assert_eq!(-a, Interval::new(-3.0, -1.0));
        assert_eq!(Interval::new(-2.0, 1.0).abs(), Interval::new(0.0, 2.0));
        for bad in [(f64::NAN, 1.0), (2.0, 1.0)] {
            let r = std::panic::catch_unwind(|| Interval::new(bad.0, bad.1));
            assert!(r.is_err(), "must reject {bad:?}");
        }
    }
}

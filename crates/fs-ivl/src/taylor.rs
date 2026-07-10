//! Taylor models (plan §6.4): a polynomial part plus a RIGOROUS interval
//! remainder — functional enclosures whose width shrinks like O(wⁿ⁺¹)
//! under subdivision where plain interval arithmetic manages O(w). The
//! containment law extends from values to FUNCTIONS: for every x in the
//! domain, f(x) ∈ P(x−c) + remainder.
//!
//! Rounding rigor follows the affine-module pattern: coefficient
//! arithmetic runs through [`Interval`], the midpoint is stored, and the
//! enclosure width is absorbed into the remainder — every absorption
//! outward-rounded. Elementary compositions carry LAGRANGE remainders
//! with derivative bounds from interval evaluation.

use crate::Interval;

/// A univariate Taylor model of fixed order: `f(x) ∈ Σ pₖ·(x−c)ᵏ + rem`
/// for all x in `domain`.
#[derive(Debug, Clone)]
pub struct TaylorModel1 {
    /// Expansion center.
    c: f64,
    /// Domain of validity.
    domain: Interval,
    /// Polynomial coefficients p₀..p_n (in powers of x − c).
    poly: Vec<f64>,
    /// Rigorous remainder: encloses truncation + all rounding.
    rem: Interval,
}

/// Absorb an interval's deviation from its midpoint into a remainder and
/// return the midpoint (the store-midpoint/absorb-width primitive).
fn split_mid(iv: Interval, rem: &mut Interval) -> f64 {
    let m = iv.midpoint();
    let dev = iv - Interval::point(m);
    *rem = *rem + dev;
    m
}

impl TaylorModel1 {
    /// The identity variable x as a Taylor model on `domain` centered at
    /// its midpoint: P(x) = c + (x−c), remainder 0.
    #[must_use]
    pub fn variable(domain: Interval, order: usize) -> TaylorModel1 {
        assert!(order >= 1, "order must be at least 1");
        let c = domain.midpoint();
        let mut poly = vec![0.0; order + 1];
        poly[0] = c;
        poly[1] = 1.0;
        TaylorModel1 {
            c,
            domain,
            poly,
            rem: Interval::point(0.0),
        }
    }

    /// A constant.
    #[must_use]
    pub fn constant(v: f64, domain: Interval, order: usize) -> TaylorModel1 {
        let mut poly = vec![0.0; order + 1];
        poly[0] = v;
        TaylorModel1 {
            c: domain.midpoint(),
            domain,
            poly,
            rem: Interval::point(0.0),
        }
    }

    /// Order (polynomial degree bound).
    #[must_use]
    pub fn order(&self) -> usize {
        self.poly.len() - 1
    }

    /// The domain.
    #[must_use]
    pub fn domain(&self) -> Interval {
        self.domain
    }

    /// The remainder interval (evidence of tightness).
    #[must_use]
    pub fn remainder(&self) -> Interval {
        self.rem
    }

    /// Interval enclosure of the model over a subdomain (must be ⊆ the
    /// model's domain): Horner in interval arithmetic + remainder.
    #[must_use]
    pub fn eval_interval(&self, x: Interval) -> Interval {
        assert!(
            self.domain.encloses(x),
            "eval domain {x:?} outside model domain {:?}",
            self.domain
        );
        let h = x - Interval::point(self.c);
        let mut acc = Interval::point(0.0);
        for &p in self.poly.iter().rev() {
            acc = acc * h + Interval::point(p);
        }
        acc + self.rem
    }

    /// Enclosure over the full domain.
    #[must_use]
    pub fn bound(&self) -> Interval {
        self.eval_interval(self.domain)
    }

    /// Scale by a constant.
    #[must_use]
    pub fn scale(&self, k: f64) -> TaylorModel1 {
        let mut rem = self.rem * Interval::point(k);
        let mut poly = vec![0.0; self.poly.len()];
        for (slot, &p) in poly.iter_mut().zip(&self.poly) {
            *slot = split_mid(Interval::point(p) * Interval::point(k), &mut rem);
        }
        TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        }
    }

    /// exp ∘ self with a Lagrange remainder: writes exp(m + g) =
    /// exp(m)·Σ gᵏ/k! + exp(sup)·|g|ⁿ⁺¹/(n+1)! where g = self − m is the
    /// centered part and the sup runs over the model's range.
    #[must_use]
    pub fn exp(&self) -> TaylorModel1 {
        let order = self.order();
        let range = self.bound();
        let m = range.midpoint();
        let em = Interval::point(fs_math::det::exp(m)) * Interval::new(1.0 - 1e-15, 1.0 + 1e-15); // covers the 3-ulp exp budget
        // g = self − m (a TM with small range).
        let g = self - &TaylorModel1::constant(m, self.domain, order);
        // Σ gᵏ/k! via Horner-free accumulation of powers.
        let mut sum = TaylorModel1::constant(1.0, self.domain, order);
        let mut gk = TaylorModel1::constant(1.0, self.domain, order);
        let mut fact = 1.0f64;
        for k in 1..=order {
            gk = &gk * &g;
            fact *= k as f64;
            sum = &sum + &gk.scale(1.0 / fact);
        }
        let mut out = sum;
        // Multiply by exp(m) rigorously (interval scalar).
        let mut rem = out.rem * em;
        let mut poly = vec![0.0; out.poly.len()];
        for (slot, &p) in poly.iter_mut().zip(&out.poly) {
            *slot = split_mid(Interval::point(p) * em, &mut rem);
        }
        // Lagrange remainder: exp(sup(range)) · |g|ⁿ⁺¹/(n+1)!
        let gmag = g.bound().abs();
        let sup = fs_math::det::exp(range.hi());
        let lag = Interval::new(-1.0, 1.0)
            * Interval::point(sup * 1.000_000_000_000_1)
            * gmag.powi(order + 1) // det-ok: Interval::powi, pinned sequential product (4xnt)
            * Interval::point(1.0 / factorial(order + 1));
        rem = rem + lag;
        out.poly = poly;
        out.rem = rem;
        out
    }

    /// sin ∘ self with the universal Lagrange bound |R| ≤ |g|ⁿ⁺¹/(n+1)!
    /// (all sine derivatives are bounded by 1).
    #[must_use]
    pub fn sin(&self) -> TaylorModel1 {
        let order = self.order();
        let range = self.bound();
        let m = range.midpoint();
        let (sm, cm) = (fs_math::det::sin(m), fs_math::det::cos(m));
        let g = self - &TaylorModel1::constant(m, self.domain, order);
        // sin(m+g) = Σ terms with derivatives cycling sin/cos at m.
        let mut sum = TaylorModel1::constant(0.0, self.domain, order);
        let mut gk = TaylorModel1::constant(1.0, self.domain, order);
        let mut fact = 1.0f64;
        for k in 0..=order {
            if k > 0 {
                gk = &gk * &g;
                fact *= k as f64;
            }
            // k-th derivative of sin at m: cycles sm, cm, −sm, −cm.
            let dk = match k % 4 {
                0 => sm,
                1 => cm,
                2 => -sm,
                _ => -cm,
            };
            // Budget slack on the strict sin/cos values (3 ulp declared).
            let dki = Interval::point(dk) + Interval::new(-2e-15, 2e-15);
            let term = {
                let mut rem = gk.rem * dki * Interval::point(1.0 / fact);
                let mut poly = vec![0.0; gk.poly.len()];
                for (slot, &p) in poly.iter_mut().zip(&gk.poly) {
                    *slot = split_mid(
                        Interval::point(p) * dki * Interval::point(1.0 / fact),
                        &mut rem,
                    );
                }
                TaylorModel1 {
                    c: self.c,
                    domain: self.domain,
                    poly,
                    rem,
                }
            };
            sum = &sum + &term;
        }
        let gmag = g.bound().abs();
        let lag = Interval::new(-1.0, 1.0)
            * gmag.powi(order + 1) // det-ok: Interval::powi, pinned sequential product (4xnt)
            * Interval::point(1.0 / factorial(order + 1));
        sum.rem = sum.rem + lag;
        sum
    }

    /// Enclosure of the polynomial part alone over the domain.
    fn poly_bound(&self) -> Interval {
        let h = self.domain - Interval::point(self.c);
        let mut acc = Interval::point(0.0);
        for &p in self.poly.iter().rev() {
            acc = acc * h + Interval::point(p);
        }
        acc
    }

    fn check_compatible(&self, o: &TaylorModel1) {
        assert!(
            self.c.to_bits() == o.c.to_bits()
                && self.domain.lo().to_bits() == o.domain.lo().to_bits()
                && self.domain.hi().to_bits() == o.domain.hi().to_bits(),
            "Taylor models must share center and domain"
        );
    }
}

impl core::ops::Add<&TaylorModel1> for &TaylorModel1 {
    type Output = TaylorModel1;
    /// Sum (same domain/center by construction discipline).
    fn add(self, o: &TaylorModel1) -> TaylorModel1 {
        self.check_compatible(o);
        let n = self.poly.len().max(o.poly.len());
        let mut rem = self.rem + o.rem;
        let mut poly = vec![0.0; n];
        for (i, slot) in poly.iter_mut().enumerate() {
            let a = self.poly.get(i).copied().unwrap_or(0.0);
            let b = o.poly.get(i).copied().unwrap_or(0.0);
            *slot = split_mid(Interval::point(a) + Interval::point(b), &mut rem);
        }
        TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        }
    }
}

impl core::ops::Sub<&TaylorModel1> for &TaylorModel1 {
    type Output = TaylorModel1;
    /// Difference.
    fn sub(self, o: &TaylorModel1) -> TaylorModel1 {
        self.check_compatible(o);
        let n = self.poly.len().max(o.poly.len());
        let mut rem = self.rem - o.rem;
        let mut poly = vec![0.0; n];
        for (i, slot) in poly.iter_mut().enumerate() {
            let a = self.poly.get(i).copied().unwrap_or(0.0);
            let b = o.poly.get(i).copied().unwrap_or(0.0);
            *slot = split_mid(Interval::point(a) - Interval::point(b), &mut rem);
        }
        TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        }
    }
}

impl core::ops::Mul<&TaylorModel1> for &TaylorModel1 {
    type Output = TaylorModel1;
    /// Product, truncated at the common order: kept cross terms absorb
    /// their rounding; DROPPED high-order terms are bounded over the
    /// centered domain and folded into the remainder; remainder cross
    /// terms use full polynomial enclosures (all outward-rounded).
    fn mul(self, o: &TaylorModel1) -> TaylorModel1 {
        self.check_compatible(o);
        let order = self.order().min(o.order());
        let h = self.domain - Interval::point(self.c);
        let mut rem = Interval::point(0.0);
        // Polynomial product with truncation.
        let mut acc: Vec<Interval> = vec![Interval::point(0.0); order + 1];
        for (i, &a) in self.poly.iter().enumerate() {
            for (j, &b) in o.poly.iter().enumerate() {
                let prod = Interval::point(a) * Interval::point(b);
                if i + j <= order {
                    acc[i + j] = acc[i + j] + prod;
                } else {
                    // Bound the dropped term a·b·h^(i+j) over the domain.
                    // det-ok: Interval::powi, pinned sequential product (4xnt)
                    rem = rem + prod * h.powi(i + j);
                }
            }
        }
        let mut poly = vec![0.0; order + 1];
        for (slot, iv) in poly.iter_mut().zip(&acc) {
            *slot = split_mid(*iv, &mut rem);
        }
        // Remainder cross terms: P₁·R₂ + P₂·R₁ + R₁·R₂ over the domain.
        let b1 = self.poly_bound();
        let b2 = o.poly_bound();
        rem = rem + b1 * o.rem + b2 * self.rem + self.rem * o.rem;
        TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        }
    }
}

impl Interval {
    /// Integer power with outward rounding (helper for remainder bounds;
    /// exact even-power tightening deliberately omitted — conservative).
    #[must_use]
    pub fn powi(self, k: usize) -> Interval {
        let mut acc = Interval::point(1.0);
        for _ in 0..k {
            acc = acc * self;
        }
        acc
    }

    /// The magnitude interval [0, max|self|].
    #[must_use]
    pub fn abs_bound(self) -> Interval {
        Interval::new(0.0, self.lo().abs().max(self.hi().abs()))
    }
}

fn factorial(k: usize) -> f64 {
    let mut f = 1.0f64;
    for i in 2..=k {
        f *= i as f64;
    }
    f
}

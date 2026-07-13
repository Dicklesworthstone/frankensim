//! Taylor models (plan §6.4): a polynomial part plus a RIGOROUS interval
//! remainder — functional enclosures whose truncation component shrinks like
//! O(wⁿ⁺¹) under subdivision until coefficient roundoff dominates, where plain
//! interval arithmetic manages O(w). The
//! containment law extends from values to FUNCTIONS: for every x in the
//! domain, f(x) ∈ P(x−c) + remainder.
//!
//! Rounding rigor follows the affine-module pattern: coefficient
//! arithmetic runs through [`Interval`], the midpoint is stored, and the
//! enclosure width is absorbed into the remainder — every absorption
//! outward-rounded. Elementary compositions carry LAGRANGE remainders
//! with derivative bounds from interval evaluation.

use crate::Interval;

/// Largest supported Taylor order.
///
/// Elementary-function remainders use `1 / (order + 1)!`. The order-169
/// remainder coefficient `1 / 170!` is normal in binary64, while `1 / 171!`
/// is subnormal. The cap also bounds each multiply to 28,900 coefficient pairs
/// and each elementary composition to about 4.9 million pairs.
pub const MAX_TAYLOR_ORDER: usize = 169;

/// Structured refusal from bounded Taylor-model construction or arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaylorModelError {
    /// The identity variable needs at least a linear coefficient.
    VariableOrderTooSmall {
        /// Order supplied by the caller.
        requested: usize,
        /// Smallest order that can represent the identity variable.
        minimum: usize,
    },
    /// The request exceeds the numerically and operationally admitted order.
    OrderTooLarge {
        /// Order supplied by the caller.
        requested: usize,
        /// Largest admitted order.
        maximum: usize,
    },
    /// Taylor domains must have finite endpoints.
    NonFiniteDomain,
    /// Constant values must be finite before allocation begins.
    NonFiniteConstant,
    /// Scale factors must be finite before arithmetic or allocation begins.
    NonFiniteScaleFactor,
    /// The requested bounded coefficient storage could not be reserved.
    AllocationFailed {
        /// Number of coefficient slots whose reservation failed.
        coefficients: usize,
    },
    /// Binary arithmetic requires exactly matching centers and domains.
    IncompatibleModels,
}

impl core::fmt::Display for TaylorModelError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::VariableOrderTooSmall { requested, minimum } => {
                write!(
                    f,
                    "Taylor variable order {requested} is below minimum {minimum}"
                )
            }
            Self::OrderTooLarge { requested, maximum } => {
                write!(f, "Taylor order {requested} exceeds maximum {maximum}")
            }
            Self::NonFiniteDomain => write!(f, "Taylor domain is not finite"),
            Self::NonFiniteConstant => write!(f, "Taylor constant is not finite"),
            Self::NonFiniteScaleFactor => write!(f, "Taylor scale factor is not finite"),
            Self::AllocationFailed { coefficients } => {
                write!(f, "could not reserve {coefficients} Taylor coefficients")
            }
            Self::IncompatibleModels => write!(f, "Taylor models have incompatible domains"),
        }
    }
}

impl std::error::Error for TaylorModelError {}

fn validate_order(order: usize, minimum: usize) -> Result<(), TaylorModelError> {
    if order < minimum {
        return Err(TaylorModelError::VariableOrderTooSmall {
            requested: order,
            minimum,
        });
    }
    if order > MAX_TAYLOR_ORDER {
        return Err(TaylorModelError::OrderTooLarge {
            requested: order,
            maximum: MAX_TAYLOR_ORDER,
        });
    }
    Ok(())
}

fn validate_domain(domain: Interval) -> Result<(), TaylorModelError> {
    if !domain.lo().is_finite() || !domain.hi().is_finite() {
        return Err(TaylorModelError::NonFiniteDomain);
    }
    Ok(())
}

fn zero_coefficients(order: usize) -> Result<Vec<f64>, TaylorModelError> {
    let coefficients = order + 1;
    let mut poly = Vec::new();
    poly.try_reserve_exact(coefficients)
        .map_err(|_| TaylorModelError::AllocationFailed { coefficients })?;
    poly.resize(coefficients, 0.0);
    Ok(poly)
}

fn zero_interval_coefficients(order: usize) -> Result<Vec<Interval>, TaylorModelError> {
    let coefficients = order + 1;
    let mut poly = Vec::new();
    poly.try_reserve_exact(coefficients)
        .map_err(|_| TaylorModelError::AllocationFailed { coefficients })?;
    poly.resize(coefficients, Interval::point(0.0));
    Ok(poly)
}

/// Preserve certificate soundness when finite inputs overflow interval
/// arithmetic. A one-sided infinite interval such as `[+inf, +inf]` does not
/// enclose the finite real value that overflowed binary64, so the only honest
/// representable result is the whole line.
fn finite_or_whole(iv: Interval) -> Interval {
    if iv.lo().is_finite() && iv.hi().is_finite() {
        iv
    } else {
        Interval::WHOLE
    }
}

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

/// Absorb an interval coefficient's deviation from its midpoint at the
/// coefficient's centered-domain power and return the midpoint.
fn split_mid(iv: Interval, h_power: Interval, rem: &mut Interval) -> f64 {
    if !iv.lo().is_finite()
        || !iv.hi().is_finite()
        || !h_power.lo().is_finite()
        || !h_power.hi().is_finite()
    {
        *rem = Interval::WHOLE;
        return 0.0;
    }
    let m = iv.midpoint();
    let dev = iv - Interval::point(m);
    *rem = finite_or_whole(*rem + dev * h_power);
    m
}

impl TaylorModel1 {
    /// The identity variable x as a Taylor model on `domain` centered at
    /// its midpoint: P(x) = c + (x−c), remainder 0.
    pub fn variable(domain: Interval, order: usize) -> Result<TaylorModel1, TaylorModelError> {
        validate_order(order, 1)?;
        validate_domain(domain)?;
        let c = domain.midpoint();
        let mut poly = zero_coefficients(order)?;
        poly[0] = c;
        poly[1] = 1.0;
        Ok(TaylorModel1 {
            c,
            domain,
            poly,
            rem: Interval::point(0.0),
        })
    }

    /// A constant.
    pub fn constant(
        v: f64,
        domain: Interval,
        order: usize,
    ) -> Result<TaylorModel1, TaylorModelError> {
        validate_order(order, 0)?;
        validate_domain(domain)?;
        if !v.is_finite() {
            return Err(TaylorModelError::NonFiniteConstant);
        }
        let mut poly = zero_coefficients(order)?;
        poly[0] = v;
        Ok(TaylorModel1 {
            c: domain.midpoint(),
            domain,
            poly,
            rem: Interval::point(0.0),
        })
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
        finite_or_whole(self.rem)
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
            acc = finite_or_whole(acc * h + Interval::point(p));
        }
        finite_or_whole(acc + self.rem)
    }

    /// Enclosure over the full domain.
    #[must_use]
    pub fn bound(&self) -> Interval {
        self.eval_interval(self.domain)
    }

    /// Scale by a constant.
    pub fn scale(&self, k: f64) -> Result<TaylorModel1, TaylorModelError> {
        if !k.is_finite() {
            return Err(TaylorModelError::NonFiniteScaleFactor);
        }
        self.scale_interval(Interval::point(k))
    }

    /// Scale by an interval known to enclose the mathematical scalar. This is
    /// private because public callers supply an exact finite binary64 scalar;
    /// elementary compositions use it for outward-rounded reciprocal
    /// factorials.
    fn scale_interval(&self, k: Interval) -> Result<TaylorModel1, TaylorModelError> {
        let mut rem = finite_or_whole(self.rem * k);
        let mut poly = zero_coefficients(self.order())?;
        let h = self.domain - Interval::point(self.c);
        let mut h_power = Interval::point(1.0);
        for (slot, &p) in poly.iter_mut().zip(&self.poly) {
            *slot = split_mid(Interval::point(p) * k, h_power, &mut rem);
            h_power = h_power * h;
        }
        Ok(TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        })
    }

    /// exp ∘ self with a Lagrange remainder: writes exp(m + g) =
    /// exp(m)·Σ gᵏ/k! + exp(sup)·|g|ⁿ⁺¹/(n+1)! where g = self − m is the
    /// centered part and the sup runs over the model's range.
    pub fn exp(&self) -> Result<TaylorModel1, TaylorModelError> {
        let order = self.order();
        let range = self.bound();
        let m = range.midpoint();
        let em = Interval::point(m).exp();
        // g = self − m (a TM with small range).
        let g = self.try_sub(&TaylorModel1::constant(m, self.domain, order)?)?;
        // Σ gᵏ/k! via Horner-free accumulation of powers.
        let mut sum = TaylorModel1::constant(1.0, self.domain, order)?;
        let mut gk = TaylorModel1::constant(1.0, self.domain, order)?;
        let mut inv_fact = Interval::point(1.0);
        for k in 1..=order {
            gk = gk.try_mul(&g)?;
            inv_fact = inv_fact / Interval::point(k as f64);
            sum = sum.try_add(&gk.scale_interval(inv_fact)?)?;
        }
        let mut out = sum;
        // Multiply by exp(m) rigorously (interval scalar).
        let mut rem = finite_or_whole(out.rem * em);
        let mut poly = zero_coefficients(order)?;
        let h = self.domain - Interval::point(self.c);
        let mut h_power = Interval::point(1.0);
        for (slot, &p) in poly.iter_mut().zip(&out.poly) {
            *slot = split_mid(Interval::point(p) * em, h_power, &mut rem);
            h_power = h_power * h;
        }
        // Lagrange remainder: exp(sup(range)) · |g|ⁿ⁺¹/(n+1)!
        let gmag = g.bound().abs();
        let exp_range = range.exp();
        let derivative_enclosure = Interval::new(-exp_range.hi(), exp_range.hi());
        let lag = derivative_enclosure
            * gmag.powi(order + 1) // det-ok: Interval::powi, pinned sequential product (4xnt)
            * reciprocal_factorial(order + 1);
        rem = finite_or_whole(rem + lag);
        out.poly = poly;
        out.rem = rem;
        Ok(out)
    }

    /// sin ∘ self with the universal Lagrange bound |R| ≤ |g|ⁿ⁺¹/(n+1)!
    /// (all sine derivatives are bounded by 1).
    pub fn sin(&self) -> Result<TaylorModel1, TaylorModelError> {
        let order = self.order();
        let range = self.bound();
        let m = range.midpoint();
        let (sm, cm) = (fs_math::det::sin(m), fs_math::det::cos(m));
        let g = self.try_sub(&TaylorModel1::constant(m, self.domain, order)?)?;
        // sin(m+g) = Σ terms with derivatives cycling sin/cos at m.
        let mut sum = TaylorModel1::constant(0.0, self.domain, order)?;
        let mut gk = TaylorModel1::constant(1.0, self.domain, order)?;
        let mut inv_fact = Interval::point(1.0);
        for k in 0..=order {
            if k > 0 {
                gk = gk.try_mul(&g)?;
                inv_fact = inv_fact / Interval::point(k as f64);
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
                let scalar = dki * inv_fact;
                let mut rem = finite_or_whole(gk.rem * scalar);
                let mut poly = zero_coefficients(order)?;
                let h = self.domain - Interval::point(self.c);
                let mut h_power = Interval::point(1.0);
                for (slot, &p) in poly.iter_mut().zip(&gk.poly) {
                    *slot = split_mid(Interval::point(p) * scalar, h_power, &mut rem);
                    h_power = h_power * h;
                }
                TaylorModel1 {
                    c: self.c,
                    domain: self.domain,
                    poly,
                    rem,
                }
            };
            sum = sum.try_add(&term)?;
        }
        let gmag = g.bound().abs();
        let lag = Interval::new(-1.0, 1.0)
            * gmag.powi(order + 1) // det-ok: Interval::powi, pinned sequential product (4xnt)
            * reciprocal_factorial(order + 1);
        sum.rem = finite_or_whole(sum.rem + lag);
        Ok(sum)
    }

    /// Enclosure of the polynomial part alone over the domain.
    fn poly_bound(&self) -> Interval {
        let h = self.domain - Interval::point(self.c);
        let mut acc = Interval::point(0.0);
        for &p in self.poly.iter().rev() {
            acc = finite_or_whole(acc * h + Interval::point(p));
        }
        finite_or_whole(acc)
    }

    /// Fallible sum. Both models must have bit-identical domains and centers.
    pub fn try_add(&self, o: &TaylorModel1) -> Result<TaylorModel1, TaylorModelError> {
        self.check_compatible(o)?;
        let order = self.order().max(o.order());
        let mut rem = finite_or_whole(self.rem + o.rem);
        let mut poly = zero_coefficients(order)?;
        let h = self.domain - Interval::point(self.c);
        let mut h_power = Interval::point(1.0);
        for (i, slot) in poly.iter_mut().enumerate() {
            let a = self.poly.get(i).copied().unwrap_or(0.0);
            let b = o.poly.get(i).copied().unwrap_or(0.0);
            *slot = split_mid(Interval::point(a) + Interval::point(b), h_power, &mut rem);
            h_power = h_power * h;
        }
        Ok(TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        })
    }

    /// Fallible difference. Both models must have bit-identical domains and centers.
    pub fn try_sub(&self, o: &TaylorModel1) -> Result<TaylorModel1, TaylorModelError> {
        self.check_compatible(o)?;
        let order = self.order().max(o.order());
        let mut rem = finite_or_whole(self.rem - o.rem);
        let mut poly = zero_coefficients(order)?;
        let h = self.domain - Interval::point(self.c);
        let mut h_power = Interval::point(1.0);
        for (i, slot) in poly.iter_mut().enumerate() {
            let a = self.poly.get(i).copied().unwrap_or(0.0);
            let b = o.poly.get(i).copied().unwrap_or(0.0);
            *slot = split_mid(Interval::point(a) - Interval::point(b), h_power, &mut rem);
            h_power = h_power * h;
        }
        Ok(TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        })
    }

    /// Fallible product, truncated at the common order.
    pub fn try_mul(&self, o: &TaylorModel1) -> Result<TaylorModel1, TaylorModelError> {
        self.check_compatible(o)?;
        let order = self.order().min(o.order());
        let h = self.domain - Interval::point(self.c);
        let mut rem = Interval::point(0.0);
        let mut acc = zero_interval_coefficients(order)?;
        let max_power = self.order() + o.order();
        let mut h_powers = zero_interval_coefficients(max_power)?;
        h_powers[0] = Interval::point(1.0);
        for power in 1..=max_power {
            h_powers[power] = h_powers[power - 1] * h;
        }
        for (i, &a) in self.poly.iter().enumerate() {
            for (j, &b) in o.poly.iter().enumerate() {
                let prod = Interval::point(a) * Interval::point(b);
                if i + j <= order {
                    acc[i + j] = acc[i + j] + prod;
                } else {
                    rem = finite_or_whole(rem + prod * h_powers[i + j]);
                }
            }
        }
        let mut poly = zero_coefficients(order)?;
        let mut h_power = Interval::point(1.0);
        for (slot, iv) in poly.iter_mut().zip(&acc) {
            *slot = split_mid(*iv, h_power, &mut rem);
            h_power = h_power * h;
        }
        let b1 = self.poly_bound();
        let b2 = o.poly_bound();
        rem = finite_or_whole(rem + b1 * o.rem + b2 * self.rem + self.rem * o.rem);
        Ok(TaylorModel1 {
            c: self.c,
            domain: self.domain,
            poly,
            rem,
        })
    }

    fn check_compatible(&self, o: &TaylorModel1) -> Result<(), TaylorModelError> {
        if self.c.to_bits() != o.c.to_bits()
            || self.domain.lo().to_bits() != o.domain.lo().to_bits()
            || self.domain.hi().to_bits() != o.domain.hi().to_bits()
        {
            return Err(TaylorModelError::IncompatibleModels);
        }
        Ok(())
    }
}

impl core::ops::Add<&TaylorModel1> for &TaylorModel1 {
    type Output = Result<TaylorModel1, TaylorModelError>;
    fn add(self, o: &TaylorModel1) -> Self::Output {
        self.try_add(o)
    }
}

impl core::ops::Sub<&TaylorModel1> for &TaylorModel1 {
    type Output = Result<TaylorModel1, TaylorModelError>;
    fn sub(self, o: &TaylorModel1) -> Self::Output {
        self.try_sub(o)
    }
}

impl core::ops::Mul<&TaylorModel1> for &TaylorModel1 {
    type Output = Result<TaylorModel1, TaylorModelError>;
    fn mul(self, o: &TaylorModel1) -> Self::Output {
        self.try_mul(o)
    }
}

impl Interval {
    /// Integer power with outward rounding (helper for remainder bounds;
    /// exact even-power tightening deliberately omitted — conservative).
    #[must_use]
    pub fn powi(self, k: usize) -> Interval {
        let mut acc = Interval::point(1.0);
        for _ in 0..k {
            acc = finite_or_whole(acc * self);
        }
        acc
    }

    /// The magnitude interval [0, max|self|].
    #[must_use]
    pub fn abs_bound(self) -> Interval {
        Interval::new(0.0, self.lo().abs().max(self.hi().abs()))
    }
}

fn reciprocal_factorial(k: usize) -> Interval {
    debug_assert!(k <= MAX_TAYLOR_ORDER + 1);
    let mut inv = Interval::point(1.0);
    for i in 2..=k {
        inv = inv / Interval::point(i as f64);
    }
    debug_assert!(inv.lo().is_finite() && inv.lo() > 0.0);
    debug_assert!(inv.hi().is_finite() && inv.hi() > 0.0);
    inv
}

#[cfg(test)]
mod admission_tests {
    use super::{MAX_TAYLOR_ORDER, reciprocal_factorial};

    #[test]
    fn maximum_remainder_reciprocal_is_normal_and_outward() {
        let admitted = reciprocal_factorial(MAX_TAYLOR_ORDER + 1);
        assert!(admitted.lo().is_normal());
        assert!(admitted.hi().is_normal());
        assert!(admitted.lo() < admitted.hi());

        let first_refused = admitted / crate::Interval::point((MAX_TAYLOR_ORDER + 2) as f64);
        assert!(first_refused.lo() > 0.0 && !first_refused.lo().is_normal());
        assert!(first_refused.hi() > 0.0 && !first_refused.hi().is_normal());

        // The old point reciprocal is known to miss the exact value for some
        // admitted orders. The interval implementation must retain width.
        let order_58 = reciprocal_factorial(58);
        assert!(order_58.lo() < order_58.hi());
    }
}

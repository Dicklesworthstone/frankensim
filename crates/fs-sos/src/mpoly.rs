//! Multivariate polynomials over `f64` ([`MPoly`]) and over outward-rounded
//! intervals ([`IPoly`]).
//!
//! `MPoly` is the modelling type: SOS programs are stated with it and the SDP
//! solver consumes its coefficients. `IPoly` is the certification type: every
//! quantity that a certificate depends on and that is *derived* by floating
//! arithmetic (a Lie derivative, a product of a multiplier with a decision
//! polynomial, a residual) is recomputed coefficient-by-coefficient with
//! `fs-ivl` intervals, so the identity a certificate proves is the identity of
//! the exact real polynomials, not of their rounded images.
//!
//! Storage is a `BTreeMap` keyed by exponent vectors, so iteration order (and
//! therefore every downstream floating reduction) is deterministic. Exact zero
//! coefficients are dropped; tiny nonzero coefficients are NEVER trimmed,
//! because a certificate must account for every term.

use std::collections::BTreeMap;

use fs_ivl::Interval;

/// The exponent vector of a monomial `x₀^e₀ · x₁^e₁ ⋯`. Its length is the
/// number of variables of the ambient ring.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Monomial(Vec<u32>);

impl Monomial {
    /// A monomial from its exponent vector.
    #[must_use]
    pub fn new(exponents: Vec<u32>) -> Monomial {
        Monomial(exponents)
    }

    /// The constant monomial `1` in `nvars` variables.
    #[must_use]
    pub fn one(nvars: usize) -> Monomial {
        Monomial(vec![0; nvars])
    }

    /// The degree-one monomial `x_k` in `nvars` variables.
    ///
    /// # Panics
    /// If `k >= nvars`.
    #[must_use]
    pub fn var(nvars: usize, k: usize) -> Monomial {
        assert!(
            k < nvars,
            "variable index {k} out of range for {nvars} variables"
        );
        let mut e = vec![0; nvars];
        e[k] = 1;
        Monomial(e)
    }

    /// The exponent vector.
    #[must_use]
    pub fn exponents(&self) -> &[u32] {
        &self.0
    }

    /// Number of variables of the ambient ring.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.0.len()
    }

    /// Total degree.
    #[must_use]
    pub fn degree(&self) -> u32 {
        self.0.iter().sum()
    }

    /// The product monomial (exponents add).
    ///
    /// # Panics
    /// If the variable counts differ.
    #[must_use]
    pub fn mul(&self, other: &Monomial) -> Monomial {
        assert_eq!(self.0.len(), other.0.len(), "monomial arity mismatch");
        Monomial(self.0.iter().zip(&other.0).map(|(a, b)| a + b).collect())
    }

    /// `self / other` when `other` divides `self`, else `None`.
    #[must_use]
    pub fn checked_div(&self, other: &Monomial) -> Option<Monomial> {
        if self.0.len() != other.0.len() {
            return None;
        }
        let mut out = Vec::with_capacity(self.0.len());
        for (a, b) in self.0.iter().zip(&other.0) {
            out.push(a.checked_sub(*b)?);
        }
        Some(Monomial(out))
    }

    /// Evaluate at a point (`x.len() == nvars`).
    #[must_use]
    pub fn eval(&self, x: &[f64]) -> f64 {
        self.0
            .iter()
            .zip(x)
            .map(|(&e, &xi)| xi.powi(e as i32))
            .product()
    }

    /// Enclose the value over an interval box: repeated multiplication so
    /// every factor is outward rounded. Even powers use the sign-aware square
    /// so `[−1, 1]²` encloses as `[0, 1]`, not `[−1, 1]`.
    #[must_use]
    pub fn eval_interval(&self, x: &[Interval]) -> Interval {
        let mut acc = Interval::point(1.0);
        for (&e, &xi) in self.0.iter().zip(x) {
            acc = acc * ipow(xi, e);
        }
        acc
    }

    /// The deterministic GRADED order used for bases: lower total degree
    /// first, then reverse-lexicographic on exponents, which lists
    /// `1, x₀, x₁, …, x₀², x₀x₁, …`.
    #[must_use]
    pub fn graded_cmp(&self, other: &Monomial) -> std::cmp::Ordering {
        self.degree()
            .cmp(&other.degree())
            .then_with(|| other.0.cmp(&self.0))
    }
}

/// Sign-aware interval square.
#[must_use]
pub(crate) fn isqr(x: Interval) -> Interval {
    let p = x * x;
    if x.contains_zero() {
        Interval::new(0.0, p.hi().max(0.0))
    } else {
        Interval::new(p.lo().max(0.0), p.hi())
    }
}

/// Interval integer power by repeated (sign-aware) squaring.
#[must_use]
pub(crate) fn ipow(x: Interval, e: u32) -> Interval {
    match e {
        0 => Interval::point(1.0),
        1 => x,
        _ => {
            let half = ipow(x, e / 2);
            let sq = isqr(half);
            if e.is_multiple_of(2) { sq } else { sq * x }
        }
    }
}

/// All monomials in `nvars` variables with total degree in
/// `min_degree..=max_degree`, in graded order.
#[must_use]
pub fn monomials_in_degree_range(nvars: usize, min_degree: u32, max_degree: u32) -> Vec<Monomial> {
    fn rec(k: usize, left: u32, cur: &mut Vec<u32>, out: &mut Vec<Monomial>, lo: u32, hi: u32) {
        if k == cur.len() {
            let d: u32 = cur.iter().sum();
            if d >= lo && d <= hi {
                out.push(Monomial(cur.clone()));
            }
            return;
        }
        for e in 0..=left {
            cur[k] = e;
            rec(k + 1, left - e, cur, out, lo, hi);
        }
        cur[k] = 0;
    }
    let mut out = Vec::new();
    let mut cur = vec![0u32; nvars];
    if nvars == 0 {
        if min_degree == 0 {
            out.push(Monomial(Vec::new()));
        }
        return out;
    }
    rec(0, max_degree, &mut cur, &mut out, min_degree, max_degree);
    out.sort_by(Monomial::graded_cmp);
    out
}

/// A multivariate polynomial with `f64` coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct MPoly {
    nvars: usize,
    terms: BTreeMap<Monomial, f64>,
}

impl MPoly {
    /// The zero polynomial in `nvars` variables.
    #[must_use]
    pub fn zero(nvars: usize) -> MPoly {
        MPoly {
            nvars,
            terms: BTreeMap::new(),
        }
    }

    /// A constant polynomial.
    #[must_use]
    pub fn constant(nvars: usize, c: f64) -> MPoly {
        MPoly::monomial(Monomial::one(nvars), c)
    }

    /// The coordinate polynomial `x_k`.
    #[must_use]
    pub fn var(nvars: usize, k: usize) -> MPoly {
        MPoly::monomial(Monomial::var(nvars, k), 1.0)
    }

    /// `c · m`.
    #[must_use]
    pub fn monomial(m: Monomial, c: f64) -> MPoly {
        let mut p = MPoly::zero(m.nvars());
        p.add_term(m, c);
        p
    }

    /// Build from `(exponents, coefficient)` pairs; repeated monomials add.
    ///
    /// # Panics
    /// If an exponent vector has the wrong length.
    #[must_use]
    pub fn from_terms<I>(nvars: usize, terms: I) -> MPoly
    where
        I: IntoIterator<Item = (Vec<u32>, f64)>,
    {
        let mut p = MPoly::zero(nvars);
        for (e, c) in terms {
            assert_eq!(e.len(), nvars, "exponent vector length != nvars");
            p.add_term(Monomial(e), c);
        }
        p
    }

    /// Add `c · m` in place (exact zero results are removed).
    ///
    /// # Panics
    /// If the monomial arity differs from the polynomial's.
    pub fn add_term(&mut self, m: Monomial, c: f64) {
        assert_eq!(m.nvars(), self.nvars, "monomial arity mismatch");
        if c == 0.0 {
            return;
        }
        match self.terms.entry(m) {
            std::collections::btree_map::Entry::Occupied(mut o) => {
                *o.get_mut() += c;
                // Remove exact cancellation so supports stay minimal.
                if *o.get() == 0.0 {
                    o.remove();
                }
            }
            std::collections::btree_map::Entry::Vacant(v) => {
                v.insert(c);
            }
        }
    }

    /// Number of variables.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// The nonzero terms in deterministic order.
    pub fn terms(&self) -> impl Iterator<Item = (&Monomial, f64)> {
        self.terms.iter().map(|(m, c)| (m, *c))
    }

    /// Number of nonzero terms.
    #[must_use]
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// Is this the zero polynomial?
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Coefficient of `m` (0 when absent).
    #[must_use]
    pub fn coeff(&self, m: &Monomial) -> f64 {
        self.terms.get(m).copied().unwrap_or(0.0)
    }

    /// Total degree (0 for the zero polynomial).
    #[must_use]
    pub fn degree(&self) -> u32 {
        self.terms.keys().map(Monomial::degree).max().unwrap_or(0)
    }

    /// Lowest total degree among nonzero terms (0 for the zero polynomial).
    #[must_use]
    pub fn min_degree(&self) -> u32 {
        self.terms.keys().map(Monomial::degree).min().unwrap_or(0)
    }

    /// Largest coefficient magnitude.
    #[must_use]
    pub fn max_abs_coeff(&self) -> f64 {
        self.terms.values().fold(0.0, |m, c| m.max(c.abs()))
    }

    /// Sum.
    #[must_use]
    pub fn add(&self, other: &MPoly) -> MPoly {
        let mut out = self.clone();
        for (m, c) in other.terms() {
            out.add_term(m.clone(), c);
        }
        out
    }

    /// Difference.
    #[must_use]
    pub fn sub(&self, other: &MPoly) -> MPoly {
        let mut out = self.clone();
        for (m, c) in other.terms() {
            out.add_term(m.clone(), -c);
        }
        out
    }

    /// Product.
    #[must_use]
    pub fn mul(&self, other: &MPoly) -> MPoly {
        assert_eq!(self.nvars, other.nvars, "polynomial arity mismatch");
        let mut out = MPoly::zero(self.nvars);
        for (ma, ca) in self.terms() {
            for (mb, cb) in other.terms() {
                out.add_term(ma.mul(mb), ca * cb);
            }
        }
        out
    }

    /// `s · self`.
    #[must_use]
    pub fn scale(&self, s: f64) -> MPoly {
        let mut out = MPoly::zero(self.nvars);
        for (m, c) in self.terms() {
            out.add_term(m.clone(), s * c);
        }
        out
    }

    /// Integer power.
    #[must_use]
    pub fn pow(&self, e: u32) -> MPoly {
        let mut out = MPoly::constant(self.nvars, 1.0);
        for _ in 0..e {
            out = out.mul(self);
        }
        out
    }

    /// Partial derivative `∂/∂x_k`.
    #[must_use]
    pub fn derivative(&self, k: usize) -> MPoly {
        let mut out = MPoly::zero(self.nvars);
        for (m, c) in self.terms() {
            let e = m.0[k];
            if e > 0 {
                let mut d = m.0.clone();
                d[k] -= 1;
                out.add_term(Monomial(d), c * f64::from(e));
            }
        }
        out
    }

    /// Evaluate at a point.
    #[must_use]
    pub fn eval(&self, x: &[f64]) -> f64 {
        assert_eq!(x.len(), self.nvars, "point dimension mismatch");
        self.terms().map(|(m, c)| c * m.eval(x)).sum()
    }

    /// Enclose the range over an interval box (or the exact value at a point
    /// box). Sound but not tight: dependency effects are not removed.
    #[must_use]
    pub fn eval_interval(&self, x: &[Interval]) -> Interval {
        assert_eq!(x.len(), self.nvars, "box dimension mismatch");
        self.terms().fold(Interval::point(0.0), |acc, (m, c)| {
            acc + Interval::point(c) * m.eval_interval(x)
        })
    }

    /// The same polynomial with point-interval coefficients.
    #[must_use]
    pub fn to_interval(&self) -> IPoly {
        IPoly {
            nvars: self.nvars,
            terms: self
                .terms
                .iter()
                .map(|(m, c)| (m.clone(), Interval::point(*c)))
                .collect(),
        }
    }
}

/// A multivariate polynomial with interval coefficients: each coefficient is
/// an enclosure of the exact real coefficient of the polynomial it models.
#[derive(Debug, Clone, PartialEq)]
pub struct IPoly {
    nvars: usize,
    terms: BTreeMap<Monomial, Interval>,
}

impl IPoly {
    /// The zero polynomial.
    #[must_use]
    pub fn zero(nvars: usize) -> IPoly {
        IPoly {
            nvars,
            terms: BTreeMap::new(),
        }
    }

    /// Number of variables.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// Add an enclosed term in place.
    ///
    /// # Panics
    /// If the monomial arity differs.
    pub fn add_term(&mut self, m: Monomial, c: Interval) {
        assert_eq!(m.nvars(), self.nvars, "monomial arity mismatch");
        if c.lo() == 0.0 && c.hi() == 0.0 {
            return;
        }
        let e = self.terms.entry(m).or_insert(Interval::point(0.0));
        *e = *e + c;
    }

    /// The terms in deterministic order (point-zero terms never stored).
    pub fn terms(&self) -> impl Iterator<Item = (&Monomial, Interval)> {
        self.terms.iter().map(|(m, c)| (m, *c))
    }

    /// Enclosure of the coefficient of `m`.
    #[must_use]
    pub fn coeff(&self, m: &Monomial) -> Interval {
        self.terms.get(m).copied().unwrap_or(Interval::point(0.0))
    }

    /// Sum.
    #[must_use]
    pub fn add(&self, other: &IPoly) -> IPoly {
        let mut out = self.clone();
        for (m, c) in other.terms() {
            out.add_term(m.clone(), c);
        }
        out
    }

    /// Difference.
    #[must_use]
    pub fn sub(&self, other: &IPoly) -> IPoly {
        let mut out = self.clone();
        for (m, c) in other.terms() {
            out.add_term(m.clone(), -c);
        }
        out
    }

    /// Product.
    #[must_use]
    pub fn mul(&self, other: &IPoly) -> IPoly {
        assert_eq!(self.nvars, other.nvars, "polynomial arity mismatch");
        let mut out = IPoly::zero(self.nvars);
        for (ma, ca) in self.terms() {
            for (mb, cb) in other.terms() {
                out.add_term(ma.mul(mb), ca * cb);
            }
        }
        out
    }

    /// Scale by an interval.
    #[must_use]
    pub fn scale(&self, s: Interval) -> IPoly {
        let mut out = IPoly::zero(self.nvars);
        for (m, c) in self.terms() {
            out.add_term(m.clone(), s * c);
        }
        out
    }

    /// Midpoint polynomial (a representative, NOT an enclosure).
    #[must_use]
    pub fn midpoint(&self) -> MPoly {
        let mut out = MPoly::zero(self.nvars);
        for (m, c) in self.terms() {
            out.add_term(m.clone(), c.midpoint());
        }
        out
    }

    /// Largest coefficient radius (0 for exact polynomials).
    #[must_use]
    pub fn max_radius(&self) -> f64 {
        self.terms.values().fold(0.0, |m, c| m.max(c.rad()))
    }

    /// Enclose the value over an interval box.
    #[must_use]
    pub fn eval_interval(&self, x: &[Interval]) -> Interval {
        assert_eq!(x.len(), self.nvars, "box dimension mismatch");
        self.terms().fold(Interval::point(0.0), |acc, (m, c)| {
            acc + c * m.eval_interval(x)
        })
    }
}

/// Enclosure of the Lie derivative `∇V · f` of `v` along the vector field
/// `f`, computed with interval arithmetic from the exact `f64` coefficients
/// of both.
///
/// # Panics
/// If `f.len() != v.nvars()` or arities disagree.
#[must_use]
pub fn lie_derivative(v: &MPoly, f: &[MPoly]) -> IPoly {
    assert_eq!(f.len(), v.nvars(), "vector field dimension != nvars");
    let mut out = IPoly::zero(v.nvars());
    for (k, fk) in f.iter().enumerate() {
        // ∂V/∂x_k: an integer multiple of a float coefficient — enclosed.
        let mut dv = IPoly::zero(v.nvars());
        for (m, c) in v.terms() {
            let e = m.exponents()[k];
            if e > 0 {
                let mut d = m.exponents().to_vec();
                d[k] -= 1;
                dv.add_term(
                    Monomial::new(d),
                    Interval::point(c) * Interval::point(f64::from(e)),
                );
            }
        }
        out = out.add(&dv.mul(&fk.to_interval()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn graded_basis_order_and_counts() {
        let b = monomials_in_degree_range(2, 0, 2);
        let e: Vec<Vec<u32>> = b.iter().map(|m| m.exponents().to_vec()).collect();
        assert_eq!(
            e,
            vec![
                vec![0, 0],
                vec![1, 0],
                vec![0, 1],
                vec![2, 0],
                vec![1, 1],
                vec![0, 2]
            ]
        );
        // C(3+3, 3) = 20 monomials of degree <= 3 in 3 variables.
        assert_eq!(monomials_in_degree_range(3, 0, 3).len(), 20);
        assert_eq!(monomials_in_degree_range(3, 1, 1).len(), 3);
    }

    #[test]
    fn arithmetic_and_derivative() {
        let x = MPoly::var(2, 0);
        let y = MPoly::var(2, 1);
        let p = x
            .mul(&x)
            .add(&x.mul(&y).scale(3.0))
            .sub(&MPoly::constant(2, 1.0));
        assert!((p.eval(&[2.0, -1.0]) - (4.0 - 6.0 - 1.0)).abs() < 1e-15);
        let dpx = p.derivative(0); // 2x + 3y
        assert!((dpx.eval(&[2.0, -1.0]) - 1.0).abs() < 1e-15);
        let q = x.add(&y).pow(2);
        assert_eq!(
            q.coeff(&Monomial::new(vec![1, 1])).to_bits(),
            2.0f64.to_bits()
        );
        assert!(x.sub(&x).is_zero());
    }

    #[test]
    fn interval_eval_encloses_and_squares_are_sign_aware() {
        let x = MPoly::var(1, 0);
        let p = x.mul(&x);
        let r = p.eval_interval(&[Interval::new(-1.0, 1.0)]);
        // Sign-aware: the enclosure of x² over [−1, 1] stays (essentially) nonnegative
        // — only the one-ulp outward nudge of the exact 0 can dip below.
        assert!(r.lo() > -1e-300 && r.hi() >= 1.0);
        let lie = lie_derivative(&p, &[x.scale(-1.0)]); // d/dt x² along ẋ=−x = −2x²
        assert!(lie.coeff(&Monomial::new(vec![2])).contains(-2.0));
    }
}

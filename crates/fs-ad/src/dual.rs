//! Forward-mode dual numbers, generic over their scalar: `Dual<T, N>` carries
//! a value and an N-wide derivative vector. Because `Dual<T, N>: Real`
//! whenever `T: Real`, NESTING gives higher orders from one implementation —
//! `Dual<Dual<f64, 1>, 1>` yields exact second directional derivatives (the
//! Hessian-vector sanity check the gradient gate uses, plan §6.6/§9.2).
//!
//! Elementary functions route through [`crate::Real`], whose `f64`
//! implementation calls fs-math's STRICT functions — so code made generic
//! over `Real` keeps cross-ISA bit-determinism on both the primal and the
//! dual path (dual arithmetic is +,−,×,÷,mul_add on top: deterministic by
//! construction).

use crate::Real;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// A dual number: value `re` plus N derivative components `eps`.
///
/// Comparison convention (standard forward AD): `PartialEq`/`PartialOrd`
/// compare the PRIMAL value only — generic code branches on values, and two
/// duals at the same point with different sensitivities must take the same
/// branch. Use field access when derivative equality matters (tests do).
#[derive(Debug, Clone, Copy)]
pub struct Dual<T, const N: usize> {
    /// The value (primal) component.
    pub re: T,
    /// Derivative components (one per seeded direction).
    pub eps: [T; N],
}

impl<T: Real, const N: usize> PartialEq for Dual<T, N> {
    fn eq(&self, other: &Self) -> bool {
        self.re == other.re
    }
}

impl<T: Real, const N: usize> PartialOrd for Dual<T, N> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.re.partial_cmp(&other.re)
    }
}

/// First-order duals over f64 — the workhorse (`Dual64<4>` etc.).
pub type Dual64<const N: usize> = Dual<f64, N>;

impl<T: Real, const N: usize> Dual<T, N> {
    /// A constant (zero derivative).
    #[must_use]
    pub fn constant(value: T) -> Self {
        Dual {
            re: value,
            eps: [T::zero(); N],
        }
    }

    /// A seeded variable: derivative 1 in direction `dir`, 0 elsewhere.
    #[must_use]
    pub fn variable(value: T, dir: usize) -> Self {
        let mut eps = [T::zero(); N];
        eps[dir] = T::one();
        Dual { re: value, eps }
    }

    /// Map over all components (value + derivative factor form): given the
    /// primal image `f` and the derivative factor `df`, produce
    /// `(f, df·eps)` — the chain rule as a helper.
    #[must_use]
    fn chain(self, f: T, df: T) -> Self {
        let mut eps = self.eps;
        for e in &mut eps {
            *e = *e * df;
        }
        Dual { re: f, eps }
    }
}

impl<T: Real, const N: usize> Add for Dual<T, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut eps = self.eps;
        for (a, b) in eps.iter_mut().zip(rhs.eps) {
            *a = *a + b;
        }
        Dual {
            re: self.re + rhs.re,
            eps,
        }
    }
}

impl<T: Real, const N: usize> Sub for Dual<T, N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut eps = self.eps;
        for (a, b) in eps.iter_mut().zip(rhs.eps) {
            *a = *a - b;
        }
        Dual {
            re: self.re - rhs.re,
            eps,
        }
    }
}

impl<T: Real, const N: usize> Mul for Dual<T, N> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // (a + ε u)(b + ε v) = ab + ε (av + bu); fused where T provides it.
        let mut eps = [T::zero(); N];
        for ((e, &u), &v) in eps.iter_mut().zip(self.eps.iter()).zip(rhs.eps.iter()) {
            *e = self.re.mul_add(v, rhs.re * u);
        }
        Dual {
            re: self.re * rhs.re,
            eps,
        }
    }
}

impl<T: Real, const N: usize> Div for Dual<T, N> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        // (u/v)' = (u' v − u v') / v² = (u' − (u/v) v') / v.
        let q = self.re / rhs.re;
        let inv = T::one() / rhs.re;
        let mut eps = [T::zero(); N];
        for ((e, &u), &v) in eps.iter_mut().zip(self.eps.iter()).zip(rhs.eps.iter()) {
            *e = (u - q * v) * inv;
        }
        Dual { re: q, eps }
    }
}

impl<T: Real, const N: usize> Neg for Dual<T, N> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut eps = self.eps;
        for e in &mut eps {
            *e = -*e;
        }
        Dual { re: -self.re, eps }
    }
}

impl<T: Real, const N: usize> Real for Dual<T, N> {
    fn zero() -> Self {
        Dual::constant(T::zero())
    }

    fn one() -> Self {
        Dual::constant(T::one())
    }

    fn from_f64(v: f64) -> Self {
        Dual::constant(T::from_f64(v))
    }

    fn value(self) -> f64 {
        self.re.value()
    }

    fn mul_add(self, a: Self, b: Self) -> Self {
        // The PRIMAL must stay fused (one rounding) so dual evaluation is
        // bit-identical to the scalar path — desugaring to mul+add here broke
        // the primal-bitwise test and would have corrupted every gradient
        // check that compares against fused kernels (fs-simd's FMA policy).
        let re = self.re.mul_add(a.re, b.re);
        let mut eps = [T::zero(); N];
        // d(x·a + b) = x'·a + x·a' + b', fused pairwise.
        for (i, e) in eps.iter_mut().enumerate() {
            *e = self
                .re
                .mul_add(a.eps[i], a.re.mul_add(self.eps[i], b.eps[i]));
        }
        Dual { re, eps }
    }

    fn recip(self) -> Self {
        let f = T::one() / self.re;
        self.chain(f, -(f * f))
    }

    fn sqrt(self) -> Self {
        let s = self.re.sqrt();
        // d√x = 1/(2√x); at x = 0 the derivative is unbounded — the +∞/NaN
        // that results is the HONEST answer (silently clamping would corrupt
        // gradient checks; documented convention).
        let df = T::one() / (T::from_f64(2.0) * s);
        self.chain(s, df)
    }

    fn abs(self) -> Self {
        // Subgradient convention at 0: derivative 0 (documented; matches the
        // symmetric choice and keeps optimizers stable at kinks).
        if self.re > T::zero() {
            self
        } else if self.re < T::zero() {
            -self
        } else {
            Dual::constant(self.re.abs())
        }
    }

    fn exp(self) -> Self {
        let f = self.re.exp();
        self.chain(f, f)
    }

    fn ln(self) -> Self {
        self.chain(self.re.ln(), T::one() / self.re)
    }

    fn sin(self) -> Self {
        self.chain(self.re.sin(), self.re.cos())
    }

    fn cos(self) -> Self {
        self.chain(self.re.cos(), -self.re.sin())
    }

    fn tanh(self) -> Self {
        let t = self.re.tanh();
        self.chain(t, T::one() - t * t)
    }

    fn powi(self, n: i32) -> Self {
        if n == 0 {
            return Self::one();
        }
        // d(xⁿ) = n·xⁿ⁻¹; computed via T's powi so the primal matches the
        // scalar path bitwise.
        let f = self.re.powi(n);
        let df = T::from_f64(f64::from(n)) * self.re.powi(n - 1);
        self.chain(f, df)
    }
}

// ---------------------------------------------------------------------------
// Seeding and extraction helpers.
// ---------------------------------------------------------------------------

/// Seed a full gradient computation: variables x[0..N] each seeded in their
/// own direction. `f` must be generic over `Real`; returns (value, gradient).
pub fn gradient<const N: usize>(
    x: [f64; N],
    f: impl Fn([Dual64<N>; N]) -> Dual64<N>,
) -> (f64, [f64; N]) {
    let mut vars = [Dual64::<N>::constant(0.0); N];
    for i in 0..N {
        vars[i] = Dual64::variable(x[i], i);
    }
    let out = f(vars);
    (out.re, out.eps)
}

/// Jacobian-vector product: seed all variables along ONE direction `v`
/// (a single dual lane regardless of dimension — the JVP shape).
pub fn jvp<const M: usize>(
    x: [f64; M],
    v: [f64; M],
    f: impl Fn([Dual64<1>; M]) -> Dual64<1>,
) -> (f64, f64) {
    let mut vars = [Dual64::<1>::constant(0.0); M];
    for i in 0..M {
        vars[i] = Dual {
            re: x[i],
            eps: [v[i]],
        };
    }
    let out = f(vars);
    (out.re, out.eps[0])
}

/// Exact second directional derivative d²/dt² f(x + t·v) via NESTED duals —
/// the Hessian-vector sanity primitive (vᵀ H v).
pub fn second_directional<const M: usize>(
    x: [f64; M],
    v: [f64; M],
    f: impl Fn([Dual<Dual64<1>, 1>; M]) -> Dual<Dual64<1>, 1>,
) -> (f64, f64, f64) {
    let mut vars = [Dual::<Dual64<1>, 1>::constant(Dual64::constant(0.0)); M];
    for i in 0..M {
        // Inner dual seeds direction v (first derivative); outer dual seeds
        // the SAME direction (differentiating the derivative).
        vars[i] = Dual {
            re: Dual {
                re: x[i],
                eps: [v[i]],
            },
            eps: [Dual {
                re: v[i],
                eps: [0.0],
            }],
        };
    }
    let out = f(vars);
    // out.re.re = f; out.re.eps[0] = ∇f·v; out.eps[0].eps[0] = vᵀHv.
    (out.re.re, out.re.eps[0], out.eps[0].eps[0])
}

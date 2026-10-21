//! Law-aware tensor polynomial chaos; no new regression or sampling engine.
//!
//! Hermite polynomials are orthonormal for N(0,1). Shifted Legendre
//! polynomials sqrt(2k+1) P_k(2u-1) are orthonormal for U(0,1).
//! Coefficient moments and support-based ANOVA therefore use the same algebra.
//! Bead: `frankensim-extreal-program-f85xj.6.7`.

use crate::pce::{PceModel, PceSensitivityError, SobolIndices, hermite_orthonormal};

/// Independent reference measure and its orthonormal polynomial family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PceBasis {
    /// Probabilists' Hermite polynomials divided by sqrt(k!), under N(0,1).
    StandardNormal,
    /// Shifted Legendre polynomials, under U(0,1); inputs must lie in [0,1].
    Uniform01,
}

impl PceBasis {
    pub(crate) fn value(self, degree: usize, x: f64) -> f64 {
        match self {
            Self::StandardNormal => hermite_orthonormal(degree, x),
            Self::Uniform01 => shifted_legendre(degree, x).0,
        }
    }

    fn derivative(self, degree: usize, x: f64) -> f64 {
        match self {
            Self::StandardNormal if degree == 0 => 0.0,
            Self::StandardNormal => {
                fs_math::det::sqrt(degree as f64) * hermite_orthonormal(degree - 1, x)
            }
            Self::Uniform01 => shifted_legendre(degree, x).1,
        }
    }
}

// Differentiate the recurrence itself, not P'_n=n(P_{n-1}-x P_n)/(1-x^2).
// This remains finite at both uniform endpoints and at every polynomial root.
fn shifted_legendre(degree: usize, u: f64) -> (f64, f64) {
    if degree == 0 {
        return (1.0, 0.0);
    }
    let x = 2.0 * u - 1.0;
    let (mut previous, mut value) = (1.0, x);
    let (mut d_previous, mut derivative) = (0.0, 2.0);
    for j in 1..degree {
        let a = (2 * j + 1) as f64;
        let b = j as f64;
        let divisor = (j + 1) as f64;
        let next = (a * x * value - b * previous) / divisor;
        let d_next = (a * (2.0 * value + x * derivative) - b * d_previous) / divisor;
        previous = value;
        value = next;
        d_previous = derivative;
        derivative = d_next;
    }
    let scale = fs_math::det::sqrt((2 * degree + 1) as f64);
    (scale * value, scale * derivative)
}

/// Immutable surrogate whose polynomial families cannot be detached from its
/// coefficients. Coordinates are independent germs, not necessarily physical
/// parameters (a correlated Gaussian transformation mixes those coordinates).
///
/// This is an algebraic surrogate, not a physical error bound or certificate.
/// The legacy Hermite model remains unchanged. Its coefficient storage and
/// family-independent ANOVA reduction are reused privately; all polynomial
/// evaluation here uses the retained families.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductPceModel {
    expansion: PceModel,
    basis: Vec<PceBasis>,
}

impl ProductPceModel {
    /// Admit 1..=16 germs, 1..=128 unique terms of total degree <=8, finite
    /// coefficients, and exactly one constant term at index zero. Sparse bases
    /// are supported. Bounds precede duplicate scans and polynomial work.
    pub fn new(
        basis: Vec<PceBasis>,
        indices: Vec<Vec<usize>>,
        coefficients: Vec<f64>,
    ) -> Result<Self, &'static str> {
        let dim = basis.len();
        if dim == 0 || dim > 16 || indices.is_empty() || indices.len() > 128 {
            return Err("product PCE needs 1..=16 germs and 1..=128 terms");
        }
        if coefficients.len() != indices.len() || coefficients.iter().any(|x| !x.is_finite()) {
            return Err("product PCE needs one finite coefficient per basis term");
        }
        for (i, alpha) in indices.iter().enumerate() {
            if alpha.len() != dim || alpha.iter().any(|&degree| degree > 8)
                || alpha.iter().sum::<usize>() > 8
            {
                return Err("product PCE basis dimension or total degree is invalid");
            }
            if indices[..i].contains(alpha) {
                return Err("product PCE basis terms must be unique");
            }
        }
        if indices[0].iter().any(|&degree| degree != 0) {
            return Err("product PCE needs its constant term first");
        }
        Ok(Self { expansion: PceModel { indices, coefficients, dim }, basis })
    }

    /// Reference measures in original germ order.
    #[must_use]
    pub fn basis(&self) -> &[PceBasis] { &self.basis }

    /// Read-only tensor multi-indices, aligned with the coefficients.
    #[must_use]
    pub fn indices(&self) -> &[Vec<usize>] { &self.expansion.indices }

    /// Read-only coefficients, beginning with the constant term.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] { &self.expansion.coefficients }

    /// Independent germ dimension.
    #[must_use]
    pub const fn dimension(&self) -> usize { self.expansion.dim }

    fn check_germ(&self, germ: &[f64]) -> Result<(), &'static str> {
        if germ.len() != self.dimension() || germ.iter().any(|x| !x.is_finite()) {
            return Err("product PCE needs finite germs of the declared dimension");
        }
        if self.basis.iter().zip(germ).any(|(&basis, &x)| {
            basis == PceBasis::Uniform01 && !(0.0..=1.0).contains(&x)
        }) {
            return Err("uniform PCE germs must lie in [0,1]");
        }
        Ok(())
    }

    /// Evaluate with typed refusal for malformed germs or unrepresentable sums.
    /// Gaussian-only arithmetic retains the legacy term/product/sum order.
    pub fn try_eval(&self, germ: &[f64]) -> Result<f64, &'static str> {
        self.check_germ(germ)?;
        let value = self.indices().iter().zip(self.coefficients()).map(|(alpha, c)| {
            let mut phi = 1.0;
            for (j, &x) in germ.iter().enumerate() {
                phi *= self.basis[j].value(alpha[j], x);
            }
            c * phi
        }).sum::<f64>();
        if value.is_finite() { Ok(value) } else { Err("product PCE evaluation overflowed") }
    }

    /// Convenience evaluation with the legacy model's scalar calling convention.
    /// # Panics
    /// Panics on the refusals documented by `try_eval`.
    #[must_use]
    pub fn eval(&self, germ: &[f64]) -> f64 {
        self.try_eval(germ).expect("product PCE evaluation refused")
    }

    /// Differentiate the surrogate in its retained germ coordinates. Prefix and
    /// suffix products avoid division by a zero polynomial value. Shifted
    /// Legendre derivatives include d(2u-1)/du=2, also at u=0 and u=1.
    pub fn try_eval_with_gradient(&self, germ: &[f64]) -> Result<(f64, Vec<f64>), &'static str> {
        self.check_germ(germ)?;
        let dim = self.dimension();
        let mut value = 0.0;
        let mut gradient = vec![0.0; dim];
        let mut factors = [1.0; 16];
        let mut prefix = [1.0; 17];
        for (alpha, &coefficient) in self.indices().iter().zip(self.coefficients()) {
            if coefficient == 0.0 { continue; }
            for j in 0..dim {
                factors[j] = self.basis[j].value(alpha[j], germ[j]);
                prefix[j + 1] = prefix[j] * factors[j];
            }
            value += coefficient * prefix[dim];
            let mut suffix = 1.0;
            for j in (0..dim).rev() {
                if alpha[j] != 0 {
                    let derivative = self.basis[j].derivative(alpha[j], germ[j]);
                    gradient[j] += coefficient * prefix[j] * derivative * suffix;
                }
                suffix *= factors[j];
            }
        }
        if value.is_finite() && gradient.iter().all(|g| g.is_finite()) {
            Ok((value, gradient))
        } else {
            Err("product PCE value or gradient overflowed")
        }
    }

    /// Convenience analytic germ gradient; does not differentiate surrogate error.
    /// # Panics
    /// Panics on the refusals documented by `try_eval_with_gradient`.
    #[must_use]
    pub fn eval_with_gradient(&self, germ: &[f64]) -> (f64, Vec<f64>) {
        self.try_eval_with_gradient(germ).expect("product PCE gradient refused")
    }

    /// Exact coefficient identity for the surrogate expectation (up to rounding).
    #[must_use]
    pub fn mean(&self) -> f64 { self.expansion.mean() }

    /// Sum of squared nonconstant coefficients; may overflow or underflow f64.
    /// The product execution report separately handles unrepresentable variance.
    #[must_use]
    pub fn variance(&self) -> f64 { self.expansion.variance() }

    /// Algebraic ANOVA under these independent germ laws, using the canonical
    /// scaled coefficient reduction. This is not physical-model validation.
    /// Under joint Gaussian inputs the effects belong to independent germs,
    /// not to the correlated physical coordinates. Constant models refuse.
    pub fn sobol_indices(&self) -> Result<SobolIndices, PceSensitivityError> {
        self.expansion.sobol_indices()
    }
}

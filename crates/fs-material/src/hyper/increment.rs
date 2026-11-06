//! Finite changes of the original elastic potentials, independent of their
//! absolute storage. This is two-point arithmetic, not infinitesimal AD.
//! Products propagate `d(ab)=mean(a) db+mean(b) da`; logarithms/exponentials
//! retain their small changes with log1p/expm1 identities. The shared DD
//! primitives retain invariant cancellations before the final f64 rounding.
//!
//! The existing small-strain expressions are reused only when BOTH endpoints
//! admit their series. Fibre activation is applied separately at each end.
//! Ogden's general spectral lift is only a second-order AD lift, so finite
//! changes instead use the exact two-eigensystem trace identity documented
//! below. Every branch still evaluates the original constitutive potential.

use super::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_ad::Real;
use fs_math::{dd::Dd, det};

#[derive(Clone, Copy, Debug)]
struct Pair {
    before: Dd,
    after: Dd,
    change: Dd,
}

impl PartialEq for Pair {
    fn eq(&self, other: &Self) -> bool {
        self.before.to_f64() == other.before.to_f64()
    }
}

impl Pair {
    fn coefficient(before: f64, after: f64) -> Self {
        Self::endpoints(Dd::from_f64(before), Dd::from_f64(after))
    }

    fn endpoints(before: Dd, after: Dd) -> Self {
        Self {
            before,
            after,
            change: after - before,
        }
    }

    fn mean(self) -> Dd {
        self.before * Dd::from_f64(0.5) + self.after * Dd::from_f64(0.5)
    }

    fn both_small(self, bound: f64) -> bool {
        self.before.abs().to_f64() <= bound && self.after.abs().to_f64() <= bound
    }

    fn positive_part(self) -> Self {
        if !self.before.lt(Dd::ZERO) && !self.after.lt(Dd::ZERO) {
            self
        } else {
            Self::endpoints(
                if self.before.lt(Dd::ZERO) {
                    Dd::ZERO
                } else {
                    self.before
                },
                if self.after.lt(Dd::ZERO) {
                    Dd::ZERO
                } else {
                    self.after
                },
            )
        }
    }

    // These elementary functions are required by Real but not used by the
    // elastic potentials. Their ordinary endpoint semantics are unchanged.
    fn endpoint_function(self, f: impl Fn(f64) -> f64) -> Self {
        Self::endpoints(
            Dd::from_f64(f(self.before.to_f64())),
            Dd::from_f64(f(self.after.to_f64())),
        )
    }
}

impl PartialOrd for Pair {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        self.before.to_f64().partial_cmp(&other.before.to_f64())
    }
}

impl core::ops::Add for Pair {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            before: self.before + other.before,
            after: self.after + other.after,
            change: self.change + other.change,
        }
    }
}

impl core::ops::Neg for Pair {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            before: -self.before,
            after: -self.after,
            change: -self.change,
        }
    }
}

impl core::ops::Sub for Pair {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + -other
    }
}

impl core::ops::Mul for Pair {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            before: self.before * other.before,
            after: self.after * other.after,
            change: self.mean() * other.change + other.mean() * self.change,
        }
    }
}

impl core::ops::Div for Pair {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let before = self.before / other.before;
        let after = self.after / other.after;
        let mean = before * Dd::from_f64(0.5) + after * Dd::from_f64(0.5);
        Self {
            before,
            after,
            change: (self.change - mean * other.change) / other.mean(),
        }
    }
}

// 2 atanh(z), for |z| <= 1/8. The 36th-order tail is below 2e-33;
// this is an absolute function bound, not an interval certificate.
fn twice_atanh(z: Dd) -> Dd {
    let square = z * z;
    let mut term = z;
    let mut sum = z;
    for n in 1..18 {
        term = term * square;
        sum = sum + term / Dd::from_f64(f64::from(2 * n + 1));
    }
    Dd::from_f64(2.0) * sum
}

fn log_dd(x: Dd) -> Dd {
    let remainder = Dd::from_f64(x.lo) / Dd::from_f64(x.hi);
    Dd::from_f64(det::ln(x.hi)) + twice_atanh(remainder / (Dd::from_f64(2.0) + remainder))
}

fn log_change(before: Dd, after: Dd, change: Dd) -> Dd {
    let z = change / (before + after);
    if z.abs().to_f64() <= 0.125 {
        twice_atanh(z)
    } else {
        log_dd(after) - log_dd(before)
    }
}

fn exprel(x: Dd) -> Dd {
    if x.abs().to_f64() > 0.125 {
        return expm1(x) / x;
    }
    let mut sum = Dd::ONE;
    let mut term = Dd::ONE;
    for n in 1..20 {
        term = term * x / Dd::from_f64(f64::from(n + 1));
        sum = sum + term;
    }
    sum
}

fn expm1(x: Dd) -> Dd {
    if x.abs().to_f64() <= 0.125 {
        x * exprel(x)
    } else {
        Dd::from_f64(det::expm1(x.hi)) + Dd::from_f64(det::exp(x.hi)) * Dd::from_f64(x.lo)
    }
}

fn exp_dd(x: Dd) -> Dd {
    Dd::from_f64(det::exp(x.hi)) * (Dd::ONE + expm1(Dd::from_f64(x.lo)))
}

impl Real for Pair {
    fn zero() -> Self {
        Self::from_f64(0.0)
    }
    fn one() -> Self {
        Self::from_f64(1.0)
    }
    fn from_f64(value: f64) -> Self {
        Self {
            before: Dd::from_f64(value),
            after: Dd::from_f64(value),
            change: Dd::ZERO,
        }
    }
    fn value(self) -> f64 {
        self.before.to_f64()
    }
    fn mul_add(self, a: Self, b: Self) -> Self {
        self * a + b
    }
    fn recip(self) -> Self {
        Self::one() / self
    }
    fn sqrt(self) -> Self {
        let before = self.before.sqrt();
        let after = self.after.sqrt();
        if before == Dd::ZERO && after == Dd::ZERO && self.change == Dd::ZERO {
            return Self::zero();
        }
        Self {
            before,
            after,
            change: self.change / (before + after),
        }
    }
    fn abs(self) -> Self {
        if self.before.lt(Dd::ZERO) && self.after.lt(Dd::ZERO) {
            -self
        } else if !self.before.lt(Dd::ZERO) && !self.after.lt(Dd::ZERO) {
            self
        } else {
            Self::endpoints(self.before.abs(), self.after.abs())
        }
    }
    fn exp(self) -> Self {
        let before = exp_dd(self.before);
        let after = exp_dd(self.after);
        let change = if self.change.lt(Dd::ZERO) {
            before * expm1(self.change)
        } else {
            -after * expm1(-self.change)
        };
        Self {
            before,
            after,
            change,
        }
    }
    fn ln(self) -> Self {
        Self {
            before: log_dd(self.before),
            after: log_dd(self.after),
            change: log_change(self.before, self.after, self.change),
        }
    }
    fn sin(self) -> Self {
        self.endpoint_function(det::sin)
    }
    fn cos(self) -> Self {
        self.endpoint_function(det::cos)
    }
    fn tanh(self) -> Self {
        self.endpoint_function(det::tanh)
    }
    fn asin(self) -> Self {
        self.endpoint_function(det::asin)
    }
    fn acos(self) -> Self {
        self.endpoint_function(det::acos)
    }
    fn atan(self) -> Self {
        self.endpoint_function(det::atan)
    }
    fn atan2(self, x: Self) -> Self {
        Self::endpoints(
            Dd::from_f64(det::atan2(self.before.to_f64(), x.before.to_f64())),
            Dd::from_f64(det::atan2(self.after.to_f64(), x.after.to_f64())),
        )
    }
    fn powi(self, n: i32) -> Self {
        let mut exponent = n.unsigned_abs();
        let mut base = self;
        let mut result = Self::one();
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = result * base;
            }
            exponent >>= 1;
            if exponent != 0 {
                base = base * base;
            }
        }
        if n < 0 { result.recip() } else { result }
    }
}

fn bulk(volumetric: Volumetric, j: Pair) -> Pair {
    if let Volumetric::SimoTaylor { kappa } = volumetric
        && !(j - Pair::one()).both_small(0.125)
    {
        let x = j - Pair::one();
        return Pair::from_f64(kappa)
            * (Pair::from_f64(0.25) * x * (j + Pair::one()) - Pair::from_f64(0.5) * j.ln());
    }
    volumetric.energy(j)
}

fn fibre(strain: Pair, k1: Pair, k2: f64) -> Pair {
    if k1.before == Dd::ZERO && k1.after == Dd::ZERO {
        return Pair::zero();
    }
    let mut active = strain.positive_part();
    if k1.before == Dd::ZERO || k1.after == Dd::ZERO {
        // An absent family stores zero regardless of its unused exponent.
        // Keep that endpoint inactive before exponentiation, as its source
        // does; zero times an overflowing unused exponential is not storage.
        active = Pair::endpoints(
            if k1.before == Dd::ZERO {
                Dd::ZERO
            } else {
                active.before
            },
            if k1.after == Dd::ZERO {
                Dd::ZERO
            } else {
                active.after
            },
        );
    }
    let exponent = Pair::from_f64(k2) * active * active;
    if exponent.both_small(0.125) {
        Pair::from_f64(0.5) * k1 * active * active * super::small_exp_relative(exponent)
    } else {
        (k1 / Pair::from_f64(2.0 * k2)) * (exponent.exp() - Pair::one())
    }
}

// Subtract the common linear part BEFORE contraction. Near Cbar=I this
// preserves quadratic Ogden work rather than cancelling O(strain) terms.
fn power_secant_minus_linear(a: f64, b: f64, r: f64) -> Dd {
    let (a, b) = (a.min(b), a.max(b));
    let a = Dd::from_f64(a);
    let b = Dd::from_f64(b);
    let logarithm = log_change(a, b, b - a);
    let power_gap = Dd::from_f64(r) * logarithm;
    if power_gap.abs().to_f64() > 0.125 {
        // Endpoint powers may both be finite although their RATIO overflows.
        // Anchor at the larger power and use only a nonpositive exponent.
        let change = if r > 0.0 {
            -exp_dd(Dd::from_f64(r) * log_dd(b)) * expm1(-power_gap)
        } else {
            exp_dd(Dd::from_f64(r) * log_dd(a)) * expm1(power_gap)
        };
        return change / (b - a) - Dd::from_f64(r);
    }
    let log_relative =
        Dd::from_f64(r - 1.0) * log_dd(a) + log_dd(exprel(power_gap)) - log_dd(exprel(logarithm));
    Dd::from_f64(r) * expm1(log_relative)
}

fn ogden(c: &[Pair; 9], j: Pair, terms: &[OgdenTerm], after_terms: &[OgdenTerm]) -> Pair {
    let strain: [Pair; 9] =
        core::array::from_fn(|i| c[i] - Pair::from_f64(if i % 4 == 0 { 1.0 } else { 0.0 }));
    // This existing series is bounded by the DEVIATORIC strain, even when
    // the common dilation is far from one. Its derivation remains C=m(I+B).
    if super::spectral::small_ogden_energy(&strain.map(|x| x.before.to_f64()), terms).is_some()
        && super::spectral::small_ogden_energy(&strain.map(|x| x.after.to_f64()), after_terms)
            .is_some()
    {
        let value = if terms == after_terms {
            super::spectral::small_ogden_energy(&strain, terms)
        } else {
            super::spectral::small_ogden_energy_with_amplitudes(&strain, terms, |i| {
                Pair::from_f64(0.25)
                    * Pair::coefficient(terms[i].mu_pa, after_terms[i].mu_pa)
                    * Pair::from_f64(terms[i].alpha)
            })
        };
        if let Some(value) = value {
            return value;
        }
    }
    let scale = (Pair::from_f64(-2.0 / 3.0) * j.ln()).exp();
    let c = c.map(|x| scale * x);
    let (d0, v0) = fs_la::eigen::jacobi_eigh(&c.map(|x| x.before.to_f64()), 3);
    let (d1, v1) = fs_la::eigen::jacobi_eigh(&c.map(|x| x.after.to_f64()), 3);
    let trace = c[0] + c[4] + c[8];
    let mut result = Pair::zero();
    for (term, after_term) in terms.iter().zip(after_terms) {
        let r = 0.5 * term.alpha;
        let mut power = trace;
        if r != 1.0 {
            power = Pair::endpoints(
                d0.iter().fold(Dd::ZERO, |sum, d| {
                    sum + exp_dd(Dd::from_f64(r) * log_dd(Dd::from_f64(*d)))
                }),
                d1.iter().fold(Dd::ZERO, |sum, d| {
                    sum + exp_dd(Dd::from_f64(r) * log_dd(Dd::from_f64(*d)))
                }),
            );
            power.change = Dd::from_f64(r) * trace.change;
            // Exact finite spectral identity for symmetric A,B:
            // tr g(B)-tr g(A) = sum_ij g[a_i,b_j]
            //   (u_i^T (B-A) v_j) (u_i^T v_j).
            // It follows from the two orthonormal resolutions of identity;
            // no eigenvector derivative or subtraction of close eigenvalues
            // is used for B-A. Coincident roots use the continuous secant.
            for i in 0..3 {
                for k in 0..3 {
                    let mut overlap = Dd::ZERO;
                    let mut change = Dd::ZERO;
                    for p in 0..3 {
                        overlap =
                            overlap + Dd::from_f64(v0[3 * p + i]) * Dd::from_f64(v1[3 * p + k]);
                        for q in 0..3 {
                            change = change
                                + Dd::from_f64(v0[3 * p + i])
                                    * c[3 * p + q].change
                                    * Dd::from_f64(v1[3 * q + k]);
                        }
                    }
                    power.change = power.change
                        + power_secant_minus_linear(d0[i], d1[k], r) * overlap * change;
                }
            }
        }
        let coefficient = if term.mu_pa == after_term.mu_pa {
            Pair::from_f64(term.mu_pa / term.alpha)
        } else {
            Pair::coefficient(term.mu_pa, after_term.mu_pa) / Pair::from_f64(term.alpha)
        };
        result = result + coefficient * (power - Pair::from_f64(3.0));
    }
    result
}

// The same native isochoric series also applies near a uniform dilation:
// its logarithm depends on B=C/mean(C)-I, not on the absolute dilation.
// Bound BOTH endpoint B tensors before using it. This preserves the absolute
// strain invariants multiplying a changed source amplitude, not just Delta C.
fn small_isochoric_input(c: &[Pair; 9]) -> Option<[Pair; 9]> {
    let a = core::array::from_fn(|i| c[i] - Pair::from_f64(if i % 4 == 0 { 1.0 } else { 0.0 }));
    let trace = c[0] + c[4] + c[8];
    let mean = trace / Pair::from_f64(3.0);
    let mut b = a.map(|x| x / mean);
    b[0] = ((c[0] - c[4]) + (c[0] - c[8])) / trace;
    b[4] = ((c[4] - c[0]) + (c[4] - c[8])) / trace;
    b[8] = -b[0] - b[4];
    let retained_mean = Pair::one() + (a[0] + a[4] + a[8]) / Pair::from_f64(3.0);
    (retained_mean.before.to_f64() > 0.0
        && retained_mean.after.to_f64() > 0.0
        && b.iter().all(|x| x.both_small(1.0 / 64.0)))
    .then_some(a)
}

pub(super) fn energy_change(
    law: &Hyperelastic,
    before: &[f64; 9],
    after: &[f64; 9],
    include_bulk: bool,
) -> f64 {
    energy_change_between(law, before, law, after, include_bulk)
}

pub(super) fn energy_change_between(
    law: &Hyperelastic,
    before: &[f64; 9],
    after_law: &Hyperelastic,
    after: &[f64; 9],
    include_bulk: bool,
) -> f64 {
    let energy = energy_pair(law, before, after_law, after, include_bulk);
    if !energy.before.to_f64().is_finite() || !energy.after.to_f64().is_finite() {
        return f64::NAN;
    }
    if energy.change.to_f64().is_finite() {
        return energy.change.to_f64();
    }
    // Widely opposed source amplitudes and strains can overflow symmetric
    // cross-products although each original energy remains representable.
    // Their independently retained endpoints provide the bounded fallback;
    // finite direct changes keep the cancellation-resistant operation above.
    (energy.after - energy.before).to_f64()
}

fn energy_pair(
    law: &Hyperelastic,
    before: &[f64; 9],
    after_law: &Hyperelastic,
    after: &[f64; 9],
    include_bulk: bool,
) -> Pair {
    let f =
        core::array::from_fn(|i| Pair::endpoints(Dd::from_f64(before[i]), Dd::from_f64(after[i])));
    let j = super::det3(&f);
    let small = if super::small_objective_strain(after, super::det3(after)).is_some() {
        super::small_objective_strain(&f, j)
    } else {
        None
    };
    let c: [Pair; 9] = core::array::from_fn(|ik| {
        (0..3).fold(Pair::zero(), |sum, r| {
            sum + f[3 * r + ik / 3] * f[3 * r + ik % 3]
        })
    });
    let trace = c[0] + c[4] + c[8];
    let scale = (Pair::from_f64(-2.0 / 3.0) * j.ln()).exp();
    let isochoric_input = small.or_else(|| match &law.model {
        // NH uses its own coupled strain expression; Ogden already owns a
        // bounded deviatoric logarithm. Neither consumes these invariants.
        HyperelasticModel::NeoHookean { .. } | HyperelasticModel::Ogden { .. } => None,
        _ => small_isochoric_input(&c),
    });
    let bar = isochoric_input.as_ref().map(super::small_isochoric_strain);
    let trace_shift = bar.as_ref().map_or_else(
        || scale * trace - Pair::from_f64(3.0),
        |strain| strain.i1_minus_three,
    );
    let remainder = match (&law.model, &after_law.model) {
        (
            HyperelasticModel::NeoHookean { mu, lambda },
            HyperelasticModel::NeoHookean {
                mu: after_mu,
                lambda: after_lambda,
            },
        ) => {
            let mu = Pair::coefficient(*mu, *after_mu);
            let lambda = if include_bulk {
                Pair::coefficient(*lambda, *after_lambda)
            } else {
                Pair::zero()
            };
            if let Some(a) = &small {
                // The original objective NH expression, with its two source
                // amplitudes carried in the same finite-pair scalar.
                let quadratic = a[0] * a[4] + a[0] * a[8] + a[4] * a[8]
                    - a[1] * a[3]
                    - a[2] * a[6]
                    - a[5] * a[7];
                let cubic = super::det3(a);
                let d = a[0] + a[4] + a[8] + quadratic + cubic;
                let (log_c, correction) = super::small_log1p_and_remainder(d);
                Pair::from_f64(0.5) * mu * (-quadratic - cubic + correction)
                    + Pair::from_f64(0.125) * lambda * log_c * log_c
            } else {
                let logarithm = j.ln();
                Pair::from_f64(0.5) * mu * (trace - Pair::from_f64(3.0)) - mu * logarithm
                    + Pair::from_f64(0.5) * lambda * logarithm * logarithm
            }
        }
        (
            HyperelasticModel::MooneyRivlin { c10, c01, .. },
            HyperelasticModel::MooneyRivlin {
                c10: after_c10,
                c01: after_c01,
                ..
            },
        ) => {
            let second_shift = if let Some(strain) = &bar {
                strain.i2_minus_three
            } else {
                // I2=|cof(F)|^2 avoids cancellation of trace(C)^2-trace(C^2).
                let mut second = Pair::zero();
                for (a, b) in [(0, 1), (0, 2), (1, 2)] {
                    for (p, q) in [(0, 1), (0, 2), (1, 2)] {
                        let minor = f[3 * a + p] * f[3 * b + q] - f[3 * a + q] * f[3 * b + p];
                        second = second + minor * minor;
                    }
                }
                scale * scale * second - Pair::from_f64(3.0)
            };
            Pair::coefficient(*c10, *after_c10) * trace_shift
                + Pair::coefficient(*c01, *after_c01) * second_shift
        }
        (
            HyperelasticModel::Gent {
                mu,
                locking_invariant,
                ..
            },
            HyperelasticModel::Gent {
                mu: after_mu,
                locking_invariant: after_locking,
                ..
            },
        ) if locking_invariant == after_locking => {
            let x = trace_shift / Pair::from_f64(*locking_invariant);
            let half_mu = Pair::from_f64(0.5) * Pair::coefficient(*mu, *after_mu);
            if x.both_small(0.125) {
                super::gent_energy_amplitude(trace_shift, half_mu, *locking_invariant)
            } else {
                -half_mu * (Pair::one() - x).ln() * Pair::from_f64(*locking_invariant)
            }
        }
        (
            HyperelasticModel::Ogden { terms, .. },
            HyperelasticModel::Ogden {
                terms: after_terms, ..
            },
        ) if terms.len() == after_terms.len()
            && terms
                .iter()
                .zip(after_terms)
                .all(|(a, b)| a.alpha == b.alpha) =>
        {
            ogden(&c, j, terms, after_terms)
        }
        (
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion,
                fibres,
                invariants,
                ..
            },
            HyperelasticModel::HolzapfelGasserOgden {
                c10: after_c10,
                k1: after_k1,
                k2: after_k2,
                dispersion: after_dispersion,
                fibres: after_fibres,
                invariants: after_invariants,
                ..
            },
        ) if k2 == after_k2
            && dispersion == after_dispersion
            && fibres == after_fibres
            && invariants == after_invariants =>
        {
            let mut energy = Pair::coefficient(*c10, *after_c10) * trace_shift;
            if *k1 != 0.0 || *after_k1 != 0.0 {
                let structural_trace = match invariants {
                    FibreInvariants::Isochoric => trace_shift,
                    FibreInvariants::Full => trace - Pair::from_f64(3.0),
                };
                let finite_bar = bar.as_ref().filter(|bar| {
                    [bar.full_quadratic, bar.full_cubic]
                        .iter()
                        .all(|x| x.before.to_f64().is_finite() && x.after.to_f64().is_finite())
                });
                for direction in fibres {
                    let shift = if let (Some(a), Some(bar)) = (&isochoric_input, finite_bar) {
                        let shift = super::structural_invariant_shift(a, direction);
                        match invariants {
                            FibreInvariants::Isochoric => {
                                super::isochoric_structural_shift(shift, bar)
                            }
                            FibreInvariants::Full => shift,
                        }
                    } else {
                        let mut invariant = Pair::zero();
                        for r in 0..3 {
                            let component = (0..3).fold(Pair::zero(), |sum, k| {
                                sum + f[3 * r + k] * Pair::from_f64(direction[k])
                            });
                            invariant = invariant + component * component;
                        }
                        (if *invariants == FibreInvariants::Isochoric {
                            scale * invariant
                        } else {
                            invariant
                        }) - Pair::one()
                    };
                    let strain = Pair::from_f64(*dispersion) * structural_trace
                        + Pair::from_f64(1.0 - 3.0 * dispersion) * shift;
                    energy = energy + fibre(strain, Pair::coefficient(*k1, *after_k1), *k2);
                }
            }
            energy
        }
        _ => {
            // Distinct nonlinear shapes have distinct endpoint potentials.
            // Evaluate each only under its own source and domain, and keep the
            // bulk separate so a common large preload cannot erase this change.
            let original = energy_pair(
                law,
                before,
                law,
                before,
                include_bulk && matches!(law.model, HyperelasticModel::NeoHookean { .. }),
            );
            let changed = energy_pair(
                after_law,
                after,
                after_law,
                after,
                include_bulk && matches!(after_law.model, HyperelasticModel::NeoHookean { .. }),
            );
            Pair::endpoints(original.before, changed.after)
        }
    };
    if include_bulk {
        remainder + bulk_sources(original_bulk(law), original_bulk(after_law), j)
    } else {
        remainder
    }
}

fn original_bulk(law: &Hyperelastic) -> Option<Volumetric> {
    match &law.model {
        HyperelasticModel::NeoHookean { .. } => None,
        HyperelasticModel::MooneyRivlin { kappa, .. } => {
            Some(Volumetric::Quadratic { kappa: *kappa })
        }
        HyperelasticModel::Gent { volumetric, .. }
        | HyperelasticModel::Ogden { volumetric, .. }
        | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => Some(*volumetric),
    }
}

fn bulk_sources(before: Option<Volumetric>, after: Option<Volumetric>, j: Pair) -> Pair {
    if before == after {
        return before.map_or_else(Pair::zero, |law| bulk(law, j));
    }
    if let (Some(before), Some(after)) = (before, after)
        && core::mem::discriminant(&before) == core::mem::discriminant(&after)
    {
        let normalized = match before {
            Volumetric::Quadratic { .. } => Volumetric::Quadratic { kappa: 1.0 },
            Volumetric::LogSquared { .. } => Volumetric::LogSquared { kappa: 1.0 },
            Volumetric::SimoTaylor { .. } => Volumetric::SimoTaylor { kappa: 1.0 },
        };
        return Pair::coefficient(before.kappa(), after.kappa()) * bulk(normalized, j);
    }
    Pair::endpoints(
        before.map_or(Dd::ZERO, |law| bulk(law, j).before),
        after.map_or(Dd::ZERO, |law| bulk(law, j).after),
    )
}

pub(super) fn determinant_change(before: &[f64; 9], after: &[f64; 9]) -> f64 {
    let f =
        core::array::from_fn(|i| Pair::endpoints(Dd::from_f64(before[i]), Dd::from_f64(after[i])));
    super::det3(&f).change.to_f64()
}

/// Original-volume bulk storage across two prescribed source snapshots.
/// The retained total volumes and reference volumes are independent inputs:
/// changing the natural volume must not reset the actual mixed volume.
pub(super) fn referenced_bulk_change(
    before: Volumetric,
    reference_before: f64,
    volume_before: f64,
    after: Volumetric,
    reference_after: f64,
    volume_after: f64,
) -> f64 {
    let pair = |a, b| Pair::endpoints(Dd::from_f64(a), Dd::from_f64(b));
    let reference = pair(reference_before, reference_after);
    let elastic = pair(volume_before, volume_after) / reference;
    let normalized = match (before, after) {
        (Volumetric::Quadratic { .. }, Volumetric::Quadratic { .. }) => {
            Some(Volumetric::Quadratic { kappa: 1.0 })
        }
        (Volumetric::LogSquared { .. }, Volumetric::LogSquared { .. }) => {
            Some(Volumetric::LogSquared { kappa: 1.0 })
        }
        (Volumetric::SimoTaylor { .. }, Volumetric::SimoTaylor { .. }) => {
            Some(Volumetric::SimoTaylor { kappa: 1.0 })
        }
        _ => None,
    };
    let endpoint_difference = || {
        // Each admitted endpoint has finite relaxed and original storage.
        // Evaluate those products separately when unlike source scales would
        // overflow the otherwise useful symmetric cross-products.
        (reference * bulk(after, elastic)).after - (reference * bulk(before, elastic)).before
    };
    let change = if let Some(normalized) = normalized {
        // All three laws are linear in kappa. Carry its finite change alongside
        // Delta(J/Jg), retaining cancellation between reference and modulus work.
        let change =
            (reference * pair(before.kappa(), after.kappa()) * bulk(normalized, elastic)).change;
        if change.to_f64().is_finite() {
            change
        } else {
            endpoint_difference()
        }
    } else {
        // A prescribed change of constitutive family has two distinct original
        // potentials. Evaluate both in DD; no interpolated surrogate is implied.
        endpoint_difference()
    };
    change.to_f64()
}

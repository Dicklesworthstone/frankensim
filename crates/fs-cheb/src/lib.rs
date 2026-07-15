//! fs-cheb — compute with FUNCTIONS as values (plan §6.5): smooth 1D
//! functions as adaptively truncated Chebyshev expansions with automatic
//! near-machine-precision degree selection, plus spectral collocation
//! differentiation matrices.
//!
//! Representation: coefficients over FIRST-KIND Chebyshev points (the
//! roots grid) — chosen deliberately so values ↔ coefficients is exactly
//! fs-fft's DCT-II/III pair (cross-ISA bit-deterministic by construction).
//! The 2D low-rank, Fourier-periodic, colleague-matrix root, and
//! Orr–Sommerfeld complex eigenproblem paths are implemented at v1
//! fixture scale. The Lobatto/DCT-I flavor and 3D low-rank functions are
//! recorded follow-up scope.
//!
//! Determinism: sampling grids, plateau detection, Clenshaw evaluation,
//! and rootfinding subdivision are all fixed-order arithmetic on strict
//! kernels — NO platform libm in any path that feeds function state
//! (workspace contract rule).

pub mod budget;
pub mod cheb2;
pub mod colleague;
pub(crate) mod fma;
pub mod fourier;
pub mod orr_sommerfeld;

pub use budget::{
    BuildRun, CHEB_BUDGET_SCHEMA_VERSION, ChebAdmission, ChebBudget, ChebError, EigsRun,
    WorkReceipt, admit_adaptive_build, admit_dirichlet_eigs, admit_root_scan,
    dirichlet_laplace_eigs_budgeted, try_build_budgeted,
};

use fs_fft::{dct2, dct3};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A smooth function on [a, b] as a truncated Chebyshev series
/// f(x) ≈ Σ' cₖ·Tₖ(t(x)) with t the affine pullback to [−1, 1]
/// (the k = 0 term is halved — the DCT-II convention).
/// Equality is exact bitwise structural equality (domain + coefficient
/// bits) — the replay/parity notion the budgeted paths certify.
#[derive(Debug, Clone, PartialEq)]
pub struct Cheb1 {
    a: f64,
    b: f64,
    /// Chebyshev coefficients, c[0] stored UN-halved (Clenshaw applies
    /// the ½ convention at evaluation).
    coeffs: Vec<f64>,
}

/// Relative plateau threshold for adaptive truncation. Sits ABOVE the
/// DCT rounding floor (~n·eps effects at large n): 10·2⁻⁵² — chasing the
/// floor itself inflates degrees ~20× for oscillatory functions (measured
/// during bring-up: sin(20x) resolved at 1090 instead of ~45).
const PLATEAU_REL: f64 = 2.2e-15;

/// Map a physical coordinate to the reference interval without forming an
/// overflowing `b - a`. Exact endpoints remain exact.
pub(crate) fn affine_to_reference(x: f64, a: f64, b: f64) -> f64 {
    if x == a {
        return -1.0;
    }
    if x == b {
        return 1.0;
    }
    let width = b - a;
    if width.is_finite() {
        // Preserve the established rounding path whenever its doubled
        // numerator is representable. A finite width alone is insufficient:
        // on [0, MAX], 2*(0.75*MAX) overflows although t=0.5 is representable.
        let doubled_delta = 2.0 * (x - a);
        if doubled_delta.is_finite() {
            return doubled_delta / width - 1.0;
        }
    }
    // The center/radius form handles either an overflowing width or an
    // overflowing doubled delta. It also preserves tiny offsets around the
    // center: x=1 on [-MAX, MAX] maps to 1/MAX instead of being absorbed while
    // forming a fractional coordinate.
    let center = f64::midpoint(a, b);
    (x - center) / half_width(a, b)
}

/// Map a reference coordinate back to `[a, b]` with a stable center/radius
/// form, avoiding the overflowing width used by the textbook affine formula.
pub(crate) fn affine_from_reference(t: f64, a: f64, b: f64) -> f64 {
    if t == -1.0 {
        return a;
    }
    if t == 1.0 {
        return b;
    }
    let center = f64::midpoint(a, b);
    let width = b - a;
    let radius = if width.is_finite() {
        width / 2.0
    } else {
        half_width(a, b)
    };
    // This is the pre-existing sampling expression whenever b-a is finite,
    // preserving ordinary-domain golden bits. Unlike convex weights, it does
    // not lose a tiny t in the rounded expressions 1±t on extreme domains.
    let mapped = center + t * radius;
    if (-1.0..=1.0).contains(&t) {
        mapped.clamp(a, b)
    } else {
        mapped
    }
}

fn half_width(a: f64, b: f64) -> f64 {
    let width = b - a;
    if width.is_finite() {
        width / 2.0
    } else {
        b / 2.0 - a / 2.0
    }
}

/// Largest power of two no greater than a positive finite value.
pub(crate) fn normalization_power_of_two(value: f64) -> f64 {
    debug_assert!(value.is_finite() && value > 0.0);
    let bits = value.to_bits();
    let exponent = (bits >> 52) & 0x7ff;
    if exponent != 0 {
        f64::from_bits(exponent << 52)
    } else {
        let mantissa = bits & ((1_u64 << 52) - 1);
        let highest_bit = 63_u32 - mantissa.leading_zeros();
        f64::from_bits(1_u64 << highest_bit)
    }
}

/// Normalize finite coefficients by a common exact power of two. Refuse an
/// exponent range for which the scaled f64 would lose information.
pub(crate) fn normalize_coefficients_exact(coeffs: &mut [f64], operation: &str) -> f64 {
    let cmax = coeffs
        .iter()
        .fold(0.0f64, |scale, coefficient| scale.max(coefficient.abs()));
    assert!(
        cmax.is_finite() && cmax > 0.0,
        "finite non-zero coefficients required for {operation}"
    );
    let scale = normalization_power_of_two(cmax);
    for coefficient in coeffs {
        let original = *coefficient;
        let normalized = original / scale;
        assert!(
            normalized.is_finite() && normalized * scale == original,
            "{operation} cannot normalize this coefficient exponent range without information loss"
        );
        *coefficient = normalized;
    }
    scale
}

/// Expansion summation of finite terms. The error-free `two_sum` residuals are
/// retained as non-overlapping partials, then collapsed with the final-rounding
/// step used by robust floating-point summation algorithms. `None` means an
/// intermediate overflow requires a common power-of-two rescale.
fn expansion_sum<F>(len: usize, term: &F, scale: f64, exact_scale: bool) -> Option<f64>
where
    F: Fn(usize) -> f64,
{
    let mut partials: Vec<f64> = Vec::with_capacity(len);
    for index in 0..len {
        let value = term(index);
        assert!(value.is_finite(), "summation terms must be finite");
        let mut x = value / scale;
        if exact_scale {
            assert!(
                x.is_finite() && x * scale == value,
                "summation cannot normalize its term exponent range without information loss"
            );
        }
        let mut write = 0usize;
        for read in 0..partials.len() {
            let mut y = partials[read];
            if x.abs() < y.abs() {
                core::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            if !hi.is_finite() {
                return None;
            }
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[write] = lo;
                write += 1;
            }
            x = hi;
        }
        partials.truncate(write);
        partials.push(x);
    }

    let mut count = partials.len();
    if count == 0 {
        return Some(0.0);
    }
    count -= 1;
    let mut hi = partials[count];
    let mut lo = 0.0;
    while count > 0 {
        let x = hi;
        count -= 1;
        let y = partials[count];
        hi = x + y;
        if !hi.is_finite() {
            return None;
        }
        let rounded_y = hi - x;
        lo = y - rounded_y;
        if lo != 0.0 {
            break;
        }
    }
    // Preserve the last half-even rounding correction when the remaining
    // partial has the same sign as the low residual.
    if count > 0
        && ((lo < 0.0 && partials[count - 1] < 0.0) || (lo > 0.0 && partials[count - 1] > 0.0))
    {
        let doubled = lo * 2.0;
        let corrected = hi + doubled;
        if corrected.is_finite() && corrected - hi == doubled {
            hi = corrected;
        }
    }
    Some(hi)
}

/// Accumulate finite terms without losing a representable cancellation
/// residual merely because every naive prefix happened to stay finite. A raw
/// expansion is tried first; only a prefix overflow triggers an exact common
/// power-of-two rescale.
pub(crate) fn stable_finite_sum<F>(len: usize, term: F, operation: &str) -> f64
where
    F: Fn(usize) -> f64,
{
    if let Some(result) = expansion_sum(len, &term, 1.0, false) {
        return result;
    }

    let mut largest = 0.0f64;
    for index in 0..len {
        let value = term(index);
        assert!(value.is_finite(), "{operation} terms must be finite");
        largest = largest.max(value.abs());
    }
    assert!(
        largest > 0.0,
        "{operation} overflowed with no non-zero term"
    );
    let scale = normalization_power_of_two(largest);
    let normalized = expansion_sum(len, &term, scale, true)
        .expect("power-of-two-normalized expansion must not overflow");
    let result = normalized * scale;
    assert!(result.is_finite(), "{operation} result must be finite");
    result
}

/// Multiply three finite factors while avoiding a lossy subnormal,
/// overflowing, or underflowing first pair whenever another pairing keeps the
/// intermediate normal. Ordinary callers retain their established first
/// pairing.
pub(crate) fn stable_finite_product3(a: f64, b: f64, c: f64, operation: &str) -> f64 {
    assert!(
        a.is_finite() && b.is_finite() && c.is_finite(),
        "{operation} factors must be finite"
    );
    let factors = [(a, b, c), (a, c, b), (b, c, a)];
    let exact_zero = a == 0.0 || b == 0.0 || c == 0.0;
    let mut finite_fallback = None;
    for (left, right, last) in factors {
        let intermediate = left * right;
        let candidate = intermediate * last;
        if !candidate.is_finite() {
            continue;
        }
        if candidate != 0.0 && finite_fallback.is_none() {
            finite_fallback = Some(candidate);
        }
        if (intermediate.is_normal() || (intermediate == 0.0 && exact_zero))
            && (candidate != 0.0 || exact_zero)
        {
            return candidate;
        }
    }
    if let Some(candidate) = finite_fallback {
        return candidate;
    }
    let direct = (a * b) * c;
    assert!(
        direct.is_finite(),
        "{operation} is not representable as finite f64"
    );
    direct
}

fn scaled_derivative_term(coefficient: f64, degree: usize, a: f64, b: f64) -> f64 {
    // The physical derivative recurrence contributes
    //     (2*k*c_k) * 2/(b-a) = 4*k*c_k/(b-a).
    // Combine the domain scale before either 2*k*c_k or b-a can overflow.
    let factor = 4.0 * degree as f64;
    let width = b - a;
    if width.is_finite() {
        let scaled_first = coefficient * factor;
        if scaled_first.is_finite() && (scaled_first != 0.0 || coefficient == 0.0) {
            scaled_first / width
        } else {
            (coefficient / width) * factor
        }
    } else {
        let scale = a.abs().max(b.abs());
        let normalized_width = (b / scale) - (a / scale);
        let scaled_first = coefficient * factor;
        if scaled_first.is_finite() && (scaled_first != 0.0 || coefficient == 0.0) {
            (scaled_first / scale) / normalized_width
        } else {
            ((coefficient / scale) * factor) / normalized_width
        }
    }
}

fn scaled_integral_term(coefficient: f64, reference_weight: f64, a: f64, b: f64) -> f64 {
    let width = b - a;
    let radius = half_width(a, b);
    if width.is_finite() && radius == 0.0 && width > 0.0 {
        // Preserve a subnormal physical width by applying the factor 1/2 to
        // the coefficient/weight side first.
        return (coefficient * (reference_weight / 2.0)) * width;
    }

    // Try the geometric scale first. If it overflows, apply the (at most
    // unit-magnitude) Chebyshev integration weight before the large radius.
    // If it underflows, the alternate ordering can still retain a subnormal.
    let geometric_first = coefficient * radius;
    let weighted = if geometric_first.is_finite()
        && (geometric_first != 0.0 || coefficient == 0.0 || radius == 0.0)
    {
        geometric_first * reference_weight
    } else {
        (coefficient * reference_weight) * radius
    };
    assert!(
        weighted.is_finite(),
        "physical Chebyshev integral term is not representable as finite f64"
    );
    weighted
}

impl Cheb1 {
    /// Build adaptively from a scalar function on [a, b]: sample at
    /// first-kind Chebyshev grids of doubling size until the trailing
    /// quarter of coefficients sits at the machine-precision plateau,
    /// then truncate. Panics (structured) if `max_degree` cannot resolve
    /// the function (non-smooth input is a modeling error here).
    #[must_use]
    pub fn build<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, max_degree: usize) -> Cheb1 {
        assert!(
            a.is_finite() && b.is_finite() && a < b,
            "domain must be finite and satisfy a < b (got [{a}, {b}])"
        );
        let degree_cap = max_degree.max(16);
        let mut n = 16usize;
        loop {
            let coeffs = Self::coeffs_at(f, a, b, n);
            let maxc = coeffs
                .iter()
                .fold(0.0f64, |m, &c| m.max(c.abs()))
                .max(f64::MIN_POSITIVE);
            let tail = &coeffs[3 * n / 4..];
            if tail.iter().all(|&c| c.abs() <= PLATEAU_REL * maxc) {
                // Truncate at the last coefficient above the plateau.
                let keep = coeffs
                    .iter()
                    .rposition(|&c| c.abs() > PLATEAU_REL * maxc)
                    .map_or(1, |p| p + 1);
                return Cheb1 {
                    a,
                    b,
                    coeffs: coeffs[..keep].to_vec(),
                };
            }
            n *= 2;
            assert!(
                n <= degree_cap,
                "function not resolved at degree {max_degree} on [{a}, {b}] \
                 (non-smooth or too oscillatory; raise max_degree or split the domain)"
            );
        }
    }

    /// Coefficients from n samples at first-kind points via DCT-II:
    /// cⱼ = (2/n)·Σₖ f(xₖ)·cos(πj(2k+1)/(2n)).
    fn coeffs_at<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> Vec<f64> {
        let vals = sample_first_kind(f, a, b, n);
        let mut c = dct2(&vals);
        let scale = 2.0 / n as f64;
        for v in &mut c {
            *v *= scale;
        }
        assert!(
            c.iter().all(|coefficient| coefficient.is_finite()),
            "Chebyshev transform coefficients must be representable as finite f64"
        );
        c
    }

    /// Construct directly from coefficients (c[0] un-halved convention).
    #[must_use]
    pub fn from_coeffs(a: f64, b: f64, coeffs: Vec<f64>) -> Cheb1 {
        assert!(
            a.is_finite() && b.is_finite() && a < b,
            "domain must be finite and satisfy a < b"
        );
        assert!(!coeffs.is_empty(), "need at least one coefficient");
        assert!(
            coeffs.iter().all(|c| c.is_finite()),
            "Cheb1 coefficients must be finite"
        );
        Cheb1 { a, b, coeffs }
    }

    /// Degree (number of retained coefficients − 1).
    #[must_use]
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// The domain.
    #[must_use]
    pub fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    /// Coefficient view (c[0] un-halved).
    #[must_use]
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// Evaluate by Clenshaw recurrence (fixed order, fused).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let t = affine_to_reference(x, self.a, self.b);
        fma::cheb_eval_dispatch(self, t)
    }

    /// The Clenshaw loop body, extracted so the x86 FMA-codegen capsule
    /// can recompile it under `target_feature` (bead nabk). MUST stay
    /// `inline(always)`: a non-inlined call keeps baseline codegen and
    /// the per-element libm `fma()` call.
    #[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
    #[inline(always)]
    pub(crate) fn eval_body(&self, t: f64) -> f64 {
        let (mut b1, mut b2) = (0.0f64, 0.0f64);
        for &c in self.coeffs.iter().skip(1).rev() {
            let b0 = (2.0 * t).mul_add(b1, c - b2);
            b2 = b1;
            b1 = b0;
        }
        // Σ' convention: half the k = 0 coefficient.
        t.mul_add(b1, 0.5f64.mul_add(self.coeffs[0], -b2))
    }

    /// Derivative as a new Chebyshev object (coefficient recurrence with
    /// the domain chain rule).
    #[must_use]
    pub fn differentiate(&self) -> Cheb1 {
        let n = self.coeffs.len();
        if n == 1 {
            return Cheb1 {
                a: self.a,
                b: self.b,
                coeffs: vec![0.0],
            };
        }
        // Preserve the established arithmetic/golden path whenever every
        // traditional intermediate is finite. Extreme domains and coefficient
        // scales fall through to the jointly-scaled recurrence below.
        let width = self.b - self.a;
        let traditional_scale = 2.0 / width;
        if width.is_finite() && traditional_scale.is_finite() && traditional_scale != 0.0 {
            let mut reference = vec![0.0f64; n];
            let mut finite = true;
            for k in (1..n).rev() {
                let above = if k + 2 < n { reference[k + 1] } else { 0.0 };
                reference[k - 1] = (2.0 * k as f64).mul_add(self.coeffs[k], above);
                finite &= reference[k - 1].is_finite();
            }
            if finite {
                let traditional: Vec<f64> = reference[..n - 1]
                    .iter()
                    .map(|value| value * traditional_scale)
                    .collect();
                if traditional.iter().all(|value| value.is_finite()) {
                    return Cheb1 {
                        a: self.a,
                        b: self.b,
                        coeffs: traditional,
                    };
                }
            }
        }
        // Run the standard recurrence after folding in the physical chain
        // rule. This avoids constructing an overflowing reference derivative
        // only to multiply it by a tiny domain scale afterward.
        let mut out = vec![0.0f64; n - 1];
        for k in (1..n).rev() {
            let above = if k + 1 < out.len() { out[k + 1] } else { 0.0 };
            out[k - 1] = scaled_derivative_term(self.coeffs[k], k, self.a, self.b) + above;
            assert!(
                out[k - 1].is_finite(),
                "physical Chebyshev derivative is not representable as finite f64 coefficients"
            );
        }
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs: out,
        }
    }

    /// Definite integral over the whole domain: only even coefficients
    /// contribute (∫₋₁¹ Tₖ = 2/(1−k²) for even k, else 0).
    #[must_use]
    pub fn integral(&self) -> f64 {
        let mut terms = Vec::with_capacity((self.coeffs.len() + 1) / 2);
        terms.push(self.coeffs[0]); // (½·c0)·2 = c0 with stored-un-halved
        let mut direct_terms_are_finite = true;
        for (k, &c) in self.coeffs.iter().enumerate().skip(2).step_by(2) {
            let term = 2.0 * c / (1.0 - (k as f64) * (k as f64));
            if !term.is_finite() {
                direct_terms_are_finite = false;
                break;
            }
            terms.push(term);
        }
        if direct_terms_are_finite {
            let reference_sum = stable_finite_sum(
                terms.len(),
                |index| terms[index],
                "Chebyshev integral accumulation",
            );
            return scaled_integral_term(reference_sum, 1.0, self.a, self.b);
        }

        // Scaling the reference-coordinate sum only after accumulation can
        // overflow even when the final physical integral is finite. Normalize
        // every coefficient by one exact power of two, sum at that bounded
        // scale, then combine the common coefficient and domain scales.
        let mut normalized = self.coeffs.clone();
        let coefficient_scale =
            normalize_coefficients_exact(&mut normalized, "Chebyshev integral fallback");
        let mut normalized_terms = Vec::with_capacity((normalized.len() + 1) / 2);
        normalized_terms.push(normalized[0]);
        for (k, &coefficient) in normalized.iter().enumerate().skip(2).step_by(2) {
            let degree = k as f64;
            let reference_weight = 2.0 / (1.0 - degree * degree);
            let term = coefficient * reference_weight;
            assert!(term.is_finite(), "normalized integral term must be finite");
            normalized_terms.push(term);
        }
        let normalized_sum = stable_finite_sum(
            normalized_terms.len(),
            |index| normalized_terms[index],
            "normalized Chebyshev integral accumulation",
        );
        let width = self.b - self.a;
        if width.is_finite() && width / 2.0 == 0.0 && width > 0.0 {
            stable_finite_product3(
                normalized_sum / 2.0,
                coefficient_scale,
                width,
                "physical Chebyshev integral",
            )
        } else {
            stable_finite_product3(
                normalized_sum,
                coefficient_scale,
                half_width(self.a, self.b),
                "physical Chebyshev integral",
            )
        }
    }

    /// Sum of two functions on the same domain.
    #[must_use]
    pub fn add(&self, o: &Cheb1) -> Cheb1 {
        assert!(self.a == o.a && self.b == o.b, "domain mismatch");
        let n = self.coeffs.len().max(o.coeffs.len());
        let mut coeffs = vec![0.0f64; n];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs.get(i).copied().unwrap_or(0.0)
                + o.coeffs.get(i).copied().unwrap_or(0.0);
            assert!(
                c.is_finite(),
                "Chebyshev sum coefficient is not representable as finite f64"
            );
        }
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs,
        }
    }

    /// Product via resampling at a grid resolving the sum of degrees.
    #[must_use]
    pub fn mul(&self, o: &Cheb1) -> Cheb1 {
        assert!(self.a == o.a && self.b == o.b, "domain mismatch");
        let n = (self.coeffs.len() + o.coeffs.len())
            .next_power_of_two()
            .max(16);
        let f = |x: f64| self.eval(x) * o.eval(x);
        let coeffs = Cheb1::coeffs_at(&f, self.a, self.b, n);
        // Truncate at plateau.
        let maxc = coeffs
            .iter()
            .fold(0.0f64, |m, &c| m.max(c.abs()))
            .max(f64::MIN_POSITIVE);
        let keep = coeffs
            .iter()
            .rposition(|&c| c.abs() > PLATEAU_REL * maxc)
            .map_or(1, |p| p + 1);
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs: coeffs[..keep].to_vec(),
        }
    }

    /// Sign-changing, numerically isolated roots on a fixed reference grid,
    /// refined by safeguarded bisection/Newton in reference coordinates.
    ///
    /// v1 limitations (documented): even-multiplicity roots and an even number
    /// of roots inside one scan cell are not found, and the returned vector is
    /// not a complete root-set certificate. A detected candidate with a small
    /// reference-space slope is refused by a conditioning heuristic; returned
    /// physical coordinates remain uncertified `f64` candidates. Use the
    /// colleague/certified APIs for candidate generation and per-box evidence.
    #[must_use]
    pub fn roots(&self) -> Vec<f64> {
        assert!(
            self.coeffs.iter().any(|&coefficient| coefficient != 0.0),
            "the identically zero polynomial has a continuum of roots"
        );
        // Common exact normalization prevents coefficient/domain scale from
        // contaminating sign and conditioning decisions. Root refinement then
        // stays in t-space, where an extreme physical domain cannot underflow
        // the derivative or magnify a premature rounded zero.
        let mut reference_coeffs = self.coeffs.clone();
        normalize_coefficients_exact(&mut reference_coeffs, "Chebyshev root scan");
        let reference = Cheb1 {
            a: -1.0,
            b: 1.0,
            coeffs: reference_coeffs,
        };
        let derivative = reference.differentiate();

        let mut roots_t = Vec::new();
        let samples = self.coeffs.len().saturating_mul(8).max(64);
        let mut prev_t = -1.0;
        let mut prev_v = reference.eval_reference(prev_t);
        assert!(prev_v.is_finite(), "root scan evaluation became non-finite");
        for k in 1..=samples {
            let t = 2.0 * (k as f64) / (samples as f64) - 1.0;
            let v = reference.eval_reference(t);
            assert!(v.is_finite(), "root scan evaluation became non-finite");
            if prev_v == 0.0 {
                reference.assert_resolvable_simple_root(&derivative, prev_t);
                roots_t.push(prev_t);
            } else if v != 0.0 && prev_v.is_sign_negative() != v.is_sign_negative() {
                roots_t.push(reference.bisect_newton_reference(&derivative, prev_t, t));
            }
            prev_t = t;
            prev_v = v;
        }
        if prev_v == 0.0 {
            reference.assert_resolvable_simple_root(&derivative, prev_t);
            roots_t.push(prev_t);
        }
        roots_t
            .into_iter()
            .map(|t| affine_from_reference(t, self.a, self.b))
            .collect()
    }

    fn eval_reference(&self, t: f64) -> f64 {
        fma::cheb_eval_dispatch(self, t)
    }

    fn assert_resolvable_simple_root(&self, derivative: &Cheb1, t: f64) {
        let slope = derivative.eval_reference(t);
        let degree_scale = self.degree().max(1) as f64;
        let slope_floor = 64.0 * 1.490_116_119_384_765_6e-8 * degree_scale;
        assert!(
            slope.is_finite() && slope.abs() > slope_floor,
            "fixed-grid root scan cannot resolve a multiple or ill-conditioned root; use colleague/certified root evidence"
        );
    }

    fn bisect_newton_reference(&self, derivative: &Cheb1, mut lo: f64, mut hi: f64) -> f64 {
        for _ in 0..40 {
            let mid = f64::midpoint(lo, hi);
            let v = self.eval_reference(mid);
            assert!(
                v.is_finite(),
                "root refinement evaluation became non-finite"
            );
            if v == 0.0 {
                self.assert_resolvable_simple_root(derivative, mid);
                return mid;
            }
            let lo_value = self.eval_reference(lo);
            assert!(
                lo_value.is_finite(),
                "root refinement evaluation became non-finite"
            );
            if lo_value != 0.0 && lo_value.is_sign_negative() != v.is_sign_negative() {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        // Newton polish from the bisection estimate.
        let mut x = f64::midpoint(lo, hi);
        for _ in 0..4 {
            let value = self.eval_reference(x);
            let dv = derivative.eval_reference(x);
            assert!(
                value.is_finite() && dv.is_finite(),
                "root refinement evaluation became non-finite"
            );
            if dv == 0.0 {
                break;
            }
            let step = value / dv;
            if !step.is_finite() {
                break;
            }
            let candidate = x - step;
            if !candidate.is_finite() || candidate < lo || candidate > hi {
                break;
            }
            x = candidate;
        }
        self.assert_resolvable_simple_root(derivative, x);
        x
    }
}

/// Sample f at the n first-kind Chebyshev points mapped to [a, b],
/// ordered k = 0..n (xₖ = cos(π(k+½)/n) descending in t).
fn sample_first_kind<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (k as f64 + 0.5) / (n as f64);
            let t = fs_math::det::cos(theta);
            let x = affine_from_reference(t, a, b);
            let y = f(x);
            assert!(y.is_finite(), "Cheb1 samples must be finite");
            y
        })
        .collect()
}

/// Synthesis: values at the n first-kind points from coefficients
/// (inverse of the analysis map; used by tests and resampling).
#[must_use]
pub fn values_from_coeffs(coeffs: &[f64], n: usize) -> Vec<f64> {
    // DCT-III with the k = 0 halving convention: dct3 already applies it.
    let mut padded = coeffs.to_vec();
    padded.resize(n, 0.0);
    dct3(&padded)
}

// ---------------------------------------------------------------------------
// Spectral collocation (Chebyshev–Lobatto differentiation matrices)
// ---------------------------------------------------------------------------

/// The n+1 Chebyshev–Lobatto points on [−1, 1], DESCENDING (x₀ = 1),
/// the classical collocation ordering.
#[must_use]
pub fn lobatto_points(n: usize) -> Vec<f64> {
    (0..=n)
        .map(|j| fs_math::det::cos(std::f64::consts::PI * (j as f64) / (n as f64)))
        .collect()
}

/// The (n+1)×(n+1) first-derivative collocation matrix on the Lobatto
/// grid (Trefethen's construction, with the NEGATIVE-SUM TRICK on the
/// diagonal — the classic accuracy fix: rows must sum to zero exactly
/// because differentiation annihilates constants).
#[must_use]
pub fn diff_matrix(n: usize) -> Vec<f64> {
    assert!(n >= 1, "need at least two points");
    let x = lobatto_points(n);
    let m = n + 1;
    let c = |i: usize| -> f64 {
        let ci = if i == 0 || i == n { 2.0 } else { 1.0 };
        if i.is_multiple_of(2) { ci } else { -ci }
    };
    let mut d = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            if i != j {
                d[i * m + j] = (c(i) / c(j)) / (x[i] - x[j]);
            }
        }
    }
    // Negative-sum trick: D[i][i] = −Σ_{j≠i} D[i][j].
    for i in 0..m {
        let mut s = 0.0f64;
        for j in 0..m {
            if j != i {
                s += d[i * m + j];
            }
        }
        d[i * m + i] = -s;
    }
    d
}

/// The D·D product body (j-inner fused dot chains) and the Rayleigh
/// matvec body, extracted so the x86 FMA-codegen capsule can recompile
/// them under `target_feature` (bead nabk). MUST stay `inline(always)`
/// — see `Cheb1::eval_body`. Chain shapes untouched: pure codegen.
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
pub(crate) fn dsq_into_body(d: &[f64], m: usize, d2: &mut [f64]) {
    for i in 0..m {
        for j in 0..m {
            let mut acc = 0.0f64;
            for l in 0..m {
                acc = d[i * m + l].mul_add(d[l * m + j], acc);
            }
            d2[i * m + j] = acc;
        }
    }
}

/// The Rayleigh matvec body (see `dsq_into_body`'s extraction note).
#[allow(clippy::inline_always)] // required to inherit the target-feature FMA capsule
#[inline(always)]
pub(crate) fn matvec_into_body(a: &[f64], v: &[f64], ni: usize, av: &mut [f64]) {
    for i in 0..ni {
        let mut acc = 0.0f64;
        for j in 0..ni {
            acc = a[i * ni + j].mul_add(v[j], acc);
        }
        av[i] = acc;
    }
}

/// Smallest `k` eigenvalues of the Dirichlet problem −u″ = λu on [−1, 1]
/// by collocation: interior block of −D², solved by SHIFT-INVERTED power
/// iteration. Shifts come from a coarse FINITE-DIFFERENCE surrogate (a
/// symmetric tridiagonal whose spectrum `fs_la::eigen::jacobi_eigh`
/// handles) — deterministic, independent of the analytic answer, and it
/// sidesteps the missing general nonsymmetric eigensolver (that solver
/// is the Orr–Sommerfeld follow-up's first deliverable).
#[must_use]
pub fn dirichlet_laplace_eigs(n: usize, k: usize) -> Vec<f64> {
    let m = n + 1;
    let d = diff_matrix(n);
    // D² then take the interior (n−1)×(n−1) block, negated.
    let mut d2 = vec![0.0f64; m * m];
    fma::dsq_into_dispatch(&d, m, &mut d2);
    let ni = n - 1;
    let mut a = vec![0.0f64; ni * ni];
    for i in 0..ni {
        for j in 0..ni {
            a[i * ni + j] = -d2[(i + 1) * m + (j + 1)];
        }
    }
    // FD surrogate on a uniform interior grid: (-1, 2, -1)/h² tridiag —
    // symmetric, so the landed dense Jacobi handles it. Its k smallest
    // eigenvalues approximate the true ones well enough to be shifts.
    let nf = 64usize;
    let h = 2.0 / (nf as f64 + 1.0);
    let mut fd = vec![0.0f64; nf * nf];
    for i in 0..nf {
        fd[i * nf + i] = 2.0 / (h * h);
        if i + 1 < nf {
            fd[i * nf + i + 1] = -1.0 / (h * h);
            fd[(i + 1) * nf + i] = -1.0 / (h * h);
        }
    }
    let (fd_eigs, _) = fs_la::eigen::jacobi_eigh(&fd, nf);
    let mut eigs = Vec::with_capacity(k);
    let mut shifted = vec![0.0f64; a.len()];
    for &fd_est in fd_eigs.iter().take(k) {
        // Shift slightly BELOW the surrogate estimate (FD underestimates
        // continuum eigenvalues; the offset keeps the shifted matrix
        // definite and the iteration locked to the intended eigenvalue).
        let mu = fd_est * 0.95;
        shifted.copy_from_slice(&a);
        for i in 0..ni {
            shifted[i * ni + i] -= mu;
        }
        let lu =
            fs_la::factor::lu(&shifted, ni).expect("shifted collocation operator is nonsingular");
        let mut v: Vec<f64> = (0..ni)
            .map(|i| 1.0 + 0.25 * (((i * 7 + 3) % 11) as f64))
            .collect();
        for _ in 0..100 {
            let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in &mut v {
                *x /= nrm;
            }
            lu.solve(&mut v);
        }
        // Rayleigh quotient λ = vᵀAv / vᵀv on the UNSHIFTED operator.
        let nrm2: f64 = v.iter().map(|x| x * x).sum();
        let mut av = vec![0.0f64; ni];
        fma::matvec_into_dispatch(&a, &v, ni, &mut av);
        eigs.push(v.iter().zip(&av).map(|(x, y)| x * y).sum::<f64>() / nrm2);
    }
    eigs
}

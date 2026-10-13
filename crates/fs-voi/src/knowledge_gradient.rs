//! One-measurement decision value for an arbitrary finite Gaussian menu.
//!
//! Conditional posterior means have the form `a[i] + b[i] Z`, `Z ~ N(0,1)`.
//! For MINIMIZATION, the knowledge gradient is
//! `min(a) - E[min_i(a[i] + b[i] Z)]`. All alternatives participate. This is
//! expected improvement of the *decision*, not variance reduction, EI against
//! a known best truth, or the full value of perfect information (EVPI).
//!
//! We integrate the lower envelope of lines. At each finite crossing `c`, its
//! slope decreases by `db > 0`; the contribution is
//! `db * (phi(|c|) - |c| Phi(-|c|))`. This sum of nonnegative hinge expectations
//! avoids subtracting two large objective levels. See Frazier, Powell & Dayanik,
//! "The Knowledge-Gradient Policy for Correlated Normal Beliefs" (2009),
//! doi:10.1287/ijoc.1080.0314, with signs reversed for minimization.
//!
//! Algebraic Gaussian integration does NOT certify a covariance model or a
//! floating-point error bound. The far-tail factor uses fixed positive
//! quadrature. Log values retain useful rankings when linear values underflow.

use fs_math::det;

/// Independent version: existing EVPI and action policies are unchanged.
pub const KNOWLEDGE_GRADIENT_VERSION: u32 = 1;
/// Bounded leaf work; callers with larger menus must explicitly decompose them.
pub const MAX_ALTERNATIVES: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnowledgeGradientError {
    InvalidShape,
    NonFiniteInput,
    NumericalRange,
    AllocationRefused,
    Interrupted,
}
impl core::fmt::Display for KnowledgeGradientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Gaussian knowledge gradient: {self:?}")
    }
}
impl std::error::Error for KnowledgeGradientError {}

/// Numerical one-step value, not a certificate or a global stopping criterion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KnowledgeGradient {
    expected_improvement: f64,
    log_expected_improvement: f64,
    envelope_segments: usize,
}
impl KnowledgeGradient {
    /// Objective-difference units. May round to zero for a positive far-tail value.
    #[must_use]
    pub const fn expected_improvement(&self) -> f64 { self.expected_improvement }
    /// Logarithm in the caller's declared units; negative infinity means exact
    /// zero for this computed envelope. A finite log with a zero linear value
    /// means underflow, NOT no information. Never compare logs of different units.
    #[must_use]
    pub const fn log_expected_improvement(&self) -> f64 { self.log_expected_improvement }
    #[must_use]
    pub const fn envelope_segments(&self) -> usize { self.envelope_segments }
    #[must_use]
    pub fn has_information(&self) -> bool { self.log_expected_improvement.is_finite() }
}

/// Bounded convenience leaf. Execution-aware callers use the controlled form.
pub fn expected_minimum_improvement(means: &[f64], slopes: &[f64]) -> Result<KnowledgeGradient, KnowledgeGradientError> {
    expected_minimum_improvement_controlled(means, slopes, || true)
}

/// Integrate the full menu's affine-normal lower envelope. Equal slopes retain
/// only their lowest intercept, not an arbitrary top-two-by-mean truncation.
/// O(n log n) envelope construction, O(n) bounded tail integrations and storage.
/// `keep_going` is called during admission, envelope construction and each tail
/// tile; returning false publishes no partial result. Sorting is bounded to 256
/// entries. Caller-provided mean/slopes are never changed.
pub fn expected_minimum_improvement_controlled(
    means: &[f64], slopes: &[f64], mut keep_going: impl FnMut() -> bool,
) -> Result<KnowledgeGradient, KnowledgeGradientError> {
    if means.is_empty() || means.len() != slopes.len() || means.len() > MAX_ALTERNATIVES {
        return Err(KnowledgeGradientError::InvalidShape);
    }
    check(&mut keep_going)?;
    let mut lines = Vec::new();
    lines.try_reserve_exact(means.len()).map_err(|_| KnowledgeGradientError::AllocationRefused)?;
    for (&a, &b) in means.iter().zip(slopes) {
        check(&mut keep_going)?;
        if !a.is_finite() || !b.is_finite() { return Err(KnowledgeGradientError::NonFiniteInput); }
        lines.push(Line { a: zero(a), b: zero(b), start: f64::NEG_INFINITY });
    }
    lines.sort_unstable_by(|x, y| y.b.total_cmp(&x.b).then_with(|| x.a.total_cmp(&y.a)));
    let mut hull: Vec<Line> = Vec::new();
    hull.try_reserve_exact(lines.len()).map_err(|_| KnowledgeGradientError::AllocationRefused)?;
    for mut line in lines {
        check(&mut keep_going)?;
        // Intercept-ascending tie order makes every later equal-slope line dominated.
        if hull.last().is_some_and(|previous| previous.b == line.b) { continue; }
        while let Some(previous) = hull.last() {
            check(&mut keep_going)?;
            let da = finite(line.a - previous.a)?;
            let db = finite(previous.b - line.b)?;
            if db <= 0.0 { return Err(KnowledgeGradientError::NumericalRange); }
            line.start = finite(da / db)?;
            if line.start > previous.start { break; }
            hull.pop();
        }
        if hull.is_empty() { line.start = f64::NEG_INFINITY; }
        hull.push(line);
    }
    let mut log_sum = f64::NEG_INFINITY;
    for pair in hull.windows(2) {
        check(&mut keep_going)?;
        let db = finite(pair[0].b - pair[1].b)?;
        let term = finite(det::ln(db) + log_normal_loss(pair[1].start.abs(), &mut keep_going)?)?;
        log_sum = if log_sum == f64::NEG_INFINITY { term } else {
            let hi = log_sum.max(term);
            hi + det::ln(1.0 + det::exp(log_sum.min(term) - hi))
        };
    }
    let expected_improvement = if log_sum == f64::NEG_INFINITY { 0.0 } else { finite(det::exp(log_sum))? };
    check(&mut keep_going)?;
    Ok(KnowledgeGradient { expected_improvement, log_expected_improvement: log_sum, envelope_segments: hull.len() })
}

#[derive(Clone, Copy)]
struct Line { a: f64, b: f64, start: f64 }
fn zero(value: f64) -> f64 { if value == 0.0 { 0.0 } else { value } }
fn finite(value: f64) -> Result<f64, KnowledgeGradientError> {
    if value.is_finite() { Ok(value) } else { Err(KnowledgeGradientError::NumericalRange) }
}
fn check(keep_going: &mut impl FnMut() -> bool) -> Result<(), KnowledgeGradientError> {
    if keep_going() { Ok(()) } else { Err(KnowledgeGradientError::Interrupted) }
}

fn log_normal_loss(t: f64, keep_going: &mut impl FnMut() -> bool) -> Result<f64, KnowledgeGradientError> {
    let log_phi = finite((-0.5 * t) * t - 0.5 * det::ln(std::f64::consts::TAU))?;
    if t <= 8.0 {
        let loss = det::exp(log_phi) - t * 0.5 * det::erfc(t / std::f64::consts::SQRT_2);
        if loss <= 0.0 { return Err(KnowledgeGradientError::NumericalRange); }
        return finite(det::ln(loss));
    }
    // phi(t)-t*Phi(-t) = phi(t)/(t+1)^2 * integral_0^inf
    // v*exp(-t/(t+1)*v - v^2/(2*(t+1)^2)) dv. The integrand is positive.
    // For t>8, truncating at 48 omits <2e-17 of this dimensionless integral.
    // GL32 on [0,48] avoids cancellation and underflow of phi(t) itself.
    let scale = t + 1.0;
    let alpha = t / scale;
    let mut integral = 0.0;
    for (node, weight) in GL32_POSITIVE {
        check(keep_going)?;
        for sign in [-1.0, 1.0] {
            let v = 24.0 * (1.0 + sign * node);
            let ratio = v / scale;
            integral += 24.0 * weight * v * det::exp(-alpha * v - 0.5 * ratio * ratio);
        }
    }
    if integral <= 0.0 { return Err(KnowledgeGradientError::NumericalRange); }
    finite(log_phi - 2.0 * det::ln(scale) + det::ln(integral))
}

// Positive half of the 32-point Gauss-Legendre rule, paired symmetrically.
#[allow(clippy::excessive_precision)]
const GL32_POSITIVE: [(f64, f64); 16] = [
    (4.83076656877383243e-02, 9.65400885147276594e-02),
    (1.44471961582796488e-01, 9.56387200792747083e-02),
    (2.39287362252137065e-01, 9.38443990808045109e-02),
    (3.31868602282127667e-01, 9.11738786957637798e-02),
    (4.21351276130635333e-01, 8.76520930044037833e-02),
    (5.06899908932229359e-01, 8.33119242269467070e-02),
    (5.87715757240762304e-01, 7.81938957870702278e-02),
    (6.63044266930215231e-01, 7.23457941088483381e-02),
    (7.32182118740289711e-01, 6.58222227763616829e-02),
    (7.94483795967942386e-01, 5.86840934785355650e-02),
    (8.49367613732569970e-01, 5.09980592623760914e-02),
    (8.96321155766052091e-01, 4.28358980222268357e-02),
    (9.34906075937739667e-01, 3.42738629130217645e-02),
    (9.64762255587506390e-01, 2.53920653092620241e-02),
    (9.85611511545268382e-01, 1.62743947309057432e-02),
    (9.97263861849481570e-01, 7.01861000947050576e-03),
];

#[cfg(test)]
mod tests;

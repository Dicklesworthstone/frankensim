//! Global positive-real certificates for stored scalar rational coefficients.
//!
//! For stable poles Re Z(i*omega)=P(x)/Q(x), x=(omega/scale)^2,
//! with Q(x)>0 for every x>=0. We build the compactified numerator directly in
//! Bernstein form on t=x/(1+x). Outward-rounded coefficients and de Casteljau
//! subdivision prove nonnegativity on the ENTIRE interval, including both tails.
//! When d=0 the common (1-t) factor is removed algebraically, never by a loss
//! floor. The separately admitted e*s is lossless and e>=0.
//!
//! The convex-hull certificate is sufficient, not complete: touching zeros,
//! ill-conditioned coefficients and exhausted budgets return Unresolved.
//! References: Bertot/Guilhot/Mahboubi, doi:10.1017/S0960129511000090;
//! Harris/Parrilo, https://arxiv.org/abs/2309.10675; improper decomposition:
//! Cherifi/Gernandt/Hinsen, doi:10.1007/s00498-023-00373-2, Lemma 19.
//! Arithmetic enclosure relies only on fs-ivl's basic IEEE operations, not on
//! transcendental-function error budgets or a floating eigensolve.

use super::{PassivityError, build_grid, qp_step, residue_coords};
use crate::model::{PoleTerm, RationalModel};
use fs_ivl::interval::Interval;

pub(crate) const MAX_GLOBAL_ORDER: usize = 64;

/// Cold-path degree and subdivision allowances. Effective caps are 64, 65536,
/// and 52 respectively. Zero subdivisions/depth still allows a direct proof.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GlobalPassivityBudget {
    /// Maximum real model order / polynomial degree.
    pub max_order: usize,
    /// Maximum number of interval bisections.
    pub max_subdivisions: usize,
    /// Maximum bisection depth.
    pub max_depth: u32,
}

impl Default for GlobalPassivityBudget {
    fn default() -> Self {
        Self {
            max_order: 64,
            max_subdivisions: 4096,
            max_depth: 48,
        }
    }
}

/// A complete enclosure proof, a proved violation, or an honest no-claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlobalPassivityVerdict {
    /// Nonnegative resistance is proved at all frequencies.
    Certified,
    /// A negative numerator value is proved somewhere in the domain.
    Violated,
    /// The finite calculation cannot establish either conclusion.
    Unresolved,
}

/// Why a global proof could not be completed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlobalPassivityUnresolved {
    /// The model/polynomial exceeds the admitted degree.
    OrderLimit,
    /// Finite scaling or coefficient enclosures cannot be represented.
    NumericRange,
    /// The subdivision allowance was exhausted.
    SubdivisionLimit,
    /// A cell reached the depth allowance, including unresolved touching zeros.
    DepthLimit,
}

/// Global result for the exact real values of the stored coefficients.
/// The producing checker determines whether the domain is the continuous axis
/// or the digital frequency circle. This does not certify rounded runtime
/// stepping, a coupled solver, or the physical accuracy of the represented law.
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalPassivityReport {
    /// Certified, proved active, or unresolved.
    pub verdict: GlobalPassivityVerdict,
    /// Bisections actually performed.
    pub subdivisions: usize,
    /// Finite interval-rechecked negative-resistance witness when representable.
    /// A negative asymptotic numerator can prove a violation without such a
    /// representable witness; in that case this is None.
    pub violation_frequency_rad_s: Option<f64>,
    /// Present exactly when the verdict is Unresolved.
    pub unresolved: Option<GlobalPassivityUnresolved>,
}

impl GlobalPassivityReport {
    fn unresolved(reason: GlobalPassivityUnresolved) -> Self {
        Self {
            verdict: GlobalPassivityVerdict::Unresolved,
            subdivisions: 0,
            violation_frequency_rad_s: None,
            unresolved: Some(reason),
        }
    }

    fn certified() -> Self {
        Self {
            verdict: GlobalPassivityVerdict::Certified,
            subdivisions: 0,
            violation_frequency_rad_s: None,
            unresolved: None,
        }
    }
}

/// Repair diagnostics; only a completed global proof permits success.
#[derive(Debug, Clone)]
pub struct GlobalRepairReport {
    /// Residue-QP steps used; zero leaves the original model unchanged.
    pub rounds: usize,
    /// L2 residue perturbation relative to the original residue norm.
    pub relative_perturbation: f64,
    /// Final QP stationarity residual (diagnostic, not the certificate).
    pub kkt_residual: f64,
    /// Global certificate of the returned model.
    pub certificate: GlobalPassivityReport,
}

/// Shared polynomial result: witness_parameter belongs to the unit interval.
#[derive(Debug, Clone)]
pub(crate) struct BernsteinReport {
    pub(crate) verdict: GlobalPassivityVerdict,
    pub(crate) subdivisions: usize,
    pub(crate) witness_parameter: Option<f64>,
    pub(crate) unresolved: Option<GlobalPassivityUnresolved>,
}

fn point(x: f64) -> Interval {
    Interval::point(x)
}
fn is_point(x: Interval, value: f64) -> bool {
    x.lo() == value && x.hi() == value
}
fn finite(coefficients: &[Interval]) -> bool {
    coefficients
        .iter()
        .all(|x| x.lo().is_finite() && x.hi().is_finite())
}

// Preserve only algebraically exact identities. An uncertain interval is never
// snapped to zero, including when two intervals happen to have equal endpoints.
pub(crate) fn interval_add(a: Interval, b: Interval) -> Interval {
    if is_point(a, 0.0) {
        return b;
    }
    if is_point(b, 0.0) {
        return a;
    }
    if a.lo() == a.hi() && b.lo() == b.hi() && a.lo() == -b.lo() {
        return point(0.0);
    }
    a + b
}

pub(crate) fn interval_mul(a: Interval, b: Interval) -> Interval {
    if is_point(a, 0.0) || is_point(b, 0.0) {
        return point(0.0);
    }
    if is_point(a, 1.0) {
        return b;
    }
    if is_point(b, 1.0) {
        return a;
    }
    if is_point(a, -1.0) {
        return -b;
    }
    if is_point(b, -1.0) {
        return -a;
    }
    a * b
}

pub(crate) fn interval_div(a: Interval, b: Interval) -> Interval {
    if !b.contains_zero() {
        if is_point(a, 0.0) {
            return point(0.0);
        }
        if is_point(b, 1.0) {
            return a;
        }
        if is_point(b, -1.0) {
            return -a;
        }
        if a.lo() == a.hi() && b.lo() == b.hi() {
            if a.lo() == b.lo() {
                return point(1.0);
            }
            if a.lo() == -b.lo() {
                return point(-1.0);
            }
        }
    }
    a / b
}

fn binomial_row(n: usize) -> Vec<u128> {
    let mut row = vec![1_u128; n + 1];
    for k in 1..=n {
        row[k] = row[k - 1] * (n + 1 - k) as u128 / k as u128;
    }
    row
}

fn integer(n: u128) -> Interval {
    let x = n as f64;
    if n <= (1_u128 << 53) {
        point(x)
    } else {
        Interval::new(fs_math::next_down(x), fs_math::next_up(x))
    }
}

fn rational(n: u128, d: u128) -> Interval {
    if n == d {
        point(1.0)
    } else {
        interval_div(integer(n), integer(d))
    }
}

/// Bernstein product; all binomial weights enclose exact integer ratios.
/// Above the hard degree cap, WHOLE forces an unresolved result without a
/// large allocation. The empty slice represents the zero polynomial.
pub(crate) fn bernstein_mul(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let (n, m) = (a.len() - 1, b.len() - 1);
    if n > MAX_GLOBAL_ORDER || m > MAX_GLOBAL_ORDER - n {
        return vec![Interval::WHOLE];
    }
    let (ca, cb, cc) = (binomial_row(n), binomial_row(m), binomial_row(n + m));
    let mut out = vec![point(0.0); n + m + 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            let weight = rational(ca[i] * cb[j], cc[i + j]);
            out[i + j] = interval_add(out[i + j], interval_mul(interval_mul(x, y), weight));
        }
    }
    out
}

/// Degree elevation preserves the same polynomial, with outward enclosures.
pub(crate) fn bernstein_elevate(a: &[Interval], degree: usize) -> Vec<Interval> {
    if degree > MAX_GLOBAL_ORDER || a.len() > degree + 1 {
        return vec![Interval::WHOLE];
    }
    if a.is_empty() {
        return vec![point(0.0); degree + 1];
    }
    let mut out = a.to_vec();
    while out.len() <= degree {
        let n = out.len();
        let mut next = vec![point(0.0); n + 1];
        next[0] = out[0];
        next[n] = out[n - 1];
        for i in 1..n {
            next[i] = interval_add(
                interval_mul(out[i - 1], rational(i as u128, n as u128)),
                interval_mul(out[i], rational((n - i) as u128, n as u128)),
            );
        }
        out = next;
    }
    out
}

/// Sum in a common Bernstein degree, preserving exact zero coefficients.
pub(crate) fn bernstein_add(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    if a.is_empty() {
        return b.to_vec();
    }
    if b.is_empty() {
        return a.to_vec();
    }
    let degree = a.len().max(b.len()) - 1;
    if degree > MAX_GLOBAL_ORDER {
        return vec![Interval::WHOLE];
    }
    let (a, b) = (bernstein_elevate(a, degree), bernstein_elevate(b, degree));
    a.into_iter()
        .zip(b)
        .map(|(x, y)| interval_add(x, y))
        .collect()
}

fn bisect(coefficients: &[Interval]) -> (Vec<Interval>, Vec<Interval>) {
    let n = coefficients.len() - 1;
    let mut work = coefficients.to_vec();
    let mut left = vec![point(0.0); n + 1];
    let mut right = left.clone();
    left[0] = work[0];
    right[n] = work[n];
    for level in 1..=n {
        for i in 0..=n - level {
            work[i] = interval_add(
                interval_mul(point(0.5), work[i]),
                interval_mul(point(0.5), work[i + 1]),
            );
        }
        left[level] = work[0];
        right[n - level] = work[n - level];
    }
    (left, right)
}

/// Sufficient nonnegativity proof by the Bernstein convex-hull property.
/// Only a proved negative endpoint reports Violated; negative INTERIOR control
/// points alone never do. Depth-first storage is O(degree * max_depth).
pub(crate) fn check_bernstein_nonnegative(
    coefficients: Vec<Interval>,
    budget: GlobalPassivityBudget,
) -> BernsteinReport {
    let report = |verdict, subdivisions, witness_parameter, unresolved| BernsteinReport {
        verdict,
        subdivisions,
        witness_parameter,
        unresolved,
    };
    let order = coefficients.len().saturating_sub(1);
    if order > budget.max_order.min(MAX_GLOBAL_ORDER) {
        return report(
            GlobalPassivityVerdict::Unresolved,
            0,
            None,
            Some(GlobalPassivityUnresolved::OrderLimit),
        );
    }
    if coefficients.is_empty() {
        return report(GlobalPassivityVerdict::Certified, 0, None, None);
    }
    let mut stack = vec![(coefficients, 0.0_f64, 1.0_f64, 0_u32)];
    let (mut subdivisions, mut unresolved) = (0, None);
    while let Some((coefficients, lo, hi, depth)) = stack.pop() {
        if !finite(&coefficients) {
            unresolved.get_or_insert(GlobalPassivityUnresolved::NumericRange);
            continue;
        }
        if coefficients[0].hi() < 0.0 {
            return report(
                GlobalPassivityVerdict::Violated,
                subdivisions,
                Some(lo),
                None,
            );
        }
        if coefficients[coefficients.len() - 1].hi() < 0.0 {
            return report(
                GlobalPassivityVerdict::Violated,
                subdivisions,
                Some(hi),
                None,
            );
        }
        if coefficients.iter().all(|x| x.lo() >= 0.0) {
            continue;
        }
        if coefficients.len() == 1 {
            unresolved.get_or_insert(GlobalPassivityUnresolved::NumericRange);
            continue;
        }
        if depth >= budget.max_depth.min(52) {
            unresolved.get_or_insert(GlobalPassivityUnresolved::DepthLimit);
            continue;
        }
        if subdivisions >= budget.max_subdivisions.min(65_536) {
            unresolved.get_or_insert(GlobalPassivityUnresolved::SubdivisionLimit);
            continue;
        }
        let mid = lo + 0.5 * (hi - lo);
        let (left, right) = bisect(&coefficients);
        subdivisions += 1;
        stack.push((right, mid, hi, depth + 1));
        stack.push((left, lo, mid, depth + 1));
    }
    report(
        if unresolved.is_some() {
            GlobalPassivityVerdict::Unresolved
        } else {
            GlobalPassivityVerdict::Certified
        },
        subdivisions,
        None,
        unresolved,
    )
}

pub(crate) fn value(coefficients: &[Interval], t: Interval) -> Interval {
    if coefficients.is_empty() {
        return point(0.0);
    }
    let mut work = coefficients.to_vec();
    let complement = interval_add(point(1.0), -t);
    for remaining in (1..coefficients.len()).rev() {
        for i in 0..remaining {
            work[i] = interval_add(
                interval_mul(complement, work[i]),
                interval_mul(t, work[i + 1]),
            );
        }
    }
    work[0]
}

fn finite_witness(coefficients: &[Interval], scale: f64, t: f64) -> Option<f64> {
    let candidate = |t: f64| {
        if !(0.0..1.0).contains(&t) {
            return None;
        }
        let omega = scale * fs_math::det::sqrt(t / (1.0 - t));
        if !omega.is_finite() {
            return None;
        }
        // Recheck at the ROUNDED returned frequency; the square root is only a
        // candidate generator and its accuracy is never part of the proof.
        let w = interval_div(point(omega), point(scale));
        let x = interval_mul(w, w);
        let u = interval_div(x, interval_add(point(1.0), x));
        let u = u.intersect(Interval::new(0.0, 1.0))?;
        (value(coefficients, u).hi() < 0.0).then_some(omega)
    };
    if let Some(w) = candidate(t) {
        return Some(w);
    }
    if t == 1.0 {
        // A negative leading coefficient proves a negative far tail. Seek a
        // finite repair constraint, without ever evaluating an infinite s.
        for depth in 1..=52 {
            if let Some(w) = candidate(1.0 - 2.0_f64.powi(-depth)) {
                return Some(w);
            }
        }
    }
    None
}

/// Prove scalar positive realness without a frequency grid or an eigensolve.
/// Stable denominators are positive by construction; all polynomial arithmetic
/// encloses the exact real values of the supplied f64 coefficients.
/// No damping, residues, poles, d or e are changed.
///
/// # Errors
/// Nonfinite/malformed models, unstable poles, or negative d/e refuse.
pub fn check_global_passivity(
    model: &RationalModel,
    budget: GlobalPassivityBudget,
) -> Result<GlobalPassivityReport, PassivityError> {
    if !model.d.is_finite() || !model.e.is_finite() {
        return Err(PassivityError::InvalidModel {
            what: "nonfinite d or e",
        });
    }
    if model.d < 0.0 || model.e < 0.0 {
        return Err(PassivityError::NegativeAsymptote);
    }
    if model.terms.len() > budget.max_order.min(MAX_GLOBAL_ORDER)
        || model.order() > budget.max_order.min(MAX_GLOBAL_ORDER)
    {
        return Ok(GlobalPassivityReport::unresolved(
            GlobalPassivityUnresolved::OrderLimit,
        ));
    }
    let mut scale = 0.0_f64;
    let mut nonzero_residue = false;
    for term in &model.terms {
        let (re, im, rr, ri) = match *term {
            PoleTerm::Real { pole, residue } => (pole, 0.0, residue, 0.0),
            PoleTerm::Pair { pole, residue } => {
                if !pole.im.is_finite() || pole.im <= 0.0 {
                    return Err(PassivityError::InvalidModel {
                        what: "pair pole needs finite positive imaginary part",
                    });
                }
                (pole.re, pole.im, residue.re, residue.im)
            }
        };
        if ![re, im, rr, ri].into_iter().all(f64::is_finite) {
            return Err(PassivityError::InvalidModel {
                what: "nonfinite pole or residue",
            });
        }
        if re >= 0.0 {
            return Err(PassivityError::Unstable);
        }
        scale = scale.max(re.abs()).max(im.abs());
        nonzero_residue |= rr != 0.0 || ri != 0.0;
    }
    if !nonzero_residue {
        return Ok(GlobalPassivityReport::certified());
    }
    // Normalize frequency and impedance by positive finite numbers. These are
    // numerical coordinate changes, not frequency/loss floors or model edits.
    let mut amplitude = model.d;
    for term in &model.terms {
        match *term {
            PoleTerm::Real { residue, .. } => amplitude = amplitude.max(residue.abs() / scale),
            PoleTerm::Pair { residue, .. } => {
                amplitude = amplitude
                    .max(residue.re.abs() / scale)
                    .max(residue.im.abs() / scale);
            }
        }
    }
    if !scale.is_finite() || scale <= 0.0 || !amplitude.is_finite() || amplitude <= 0.0 {
        return Ok(GlobalPassivityReport::unresolved(
            GlobalPassivityUnresolved::NumericRange,
        ));
    }
    let normalize_pole = |p: f64| interval_div(point(p), point(scale));
    let normalize_residue =
        |r: f64| interval_div(interval_div(point(r), point(scale)), point(amplitude));
    let mut factors = Vec::with_capacity(model.terms.len());
    let mut positive_terms = true;
    for term in &model.terms {
        let (denominator, numerator) = match *term {
            PoleTerm::Real { pole, residue } => {
                let (a, r) = (normalize_pole(pole), normalize_residue(residue));
                (
                    vec![interval_mul(a, a), point(1.0)],
                    vec![-interval_mul(a, r)],
                )
            }
            PoleTerm::Pair { pole, residue } => {
                let (a, b) = (normalize_pole(pole.re), normalize_pole(pole.im));
                let (r, i) = (normalize_residue(residue.re), normalize_residue(residue.im));
                let (a2, b2) = (interval_mul(a, a), interval_mul(b, b));
                let norm = interval_add(a2, b2);
                let (ar, bi) = (interval_mul(a, r), interval_mul(b, i));
                (
                    vec![interval_mul(norm, norm), interval_add(a2, -b2), point(1.0)],
                    vec![
                        -interval_mul(point(2.0), interval_mul(interval_add(ar, bi), norm)),
                        interval_mul(point(2.0), interval_add(bi, -ar)),
                    ],
                )
            }
        };
        if !finite(&denominator) || !finite(&numerator) {
            return Ok(GlobalPassivityReport::unresolved(
                GlobalPassivityUnresolved::NumericRange,
            ));
        }
        positive_terms &= numerator.iter().all(|x| x.lo() >= 0.0);
        factors.push((denominator, numerator));
    }
    // A sum of individually nonnegative rational terms also has a global
    // proof; this avoids unnecessary cancellation around narrow resonances.
    if positive_terms {
        return Ok(GlobalPassivityReport::certified());
    }
    let mut denominator = vec![point(1.0)];
    let mut numerator = if model.d == 0.0 {
        Vec::new()
    } else {
        vec![interval_div(point(model.d), point(amplitude))]
    };
    for (q, mut p) in factors {
        if model.d != 0.0 {
            p = bernstein_mul(&p, &[point(1.0), point(0.0)]);
        }
        numerator = bernstein_add(
            &bernstein_mul(&numerator, &q),
            &bernstein_mul(&denominator, &p),
        );
        denominator = bernstein_mul(&denominator, &q);
    }
    // With d=0, numerator has degree order-1: the exact common (1-t)
    // was removed before construction, so strictly proper tails remain visible.
    let proof = check_bernstein_nonnegative(numerator.clone(), budget);
    let witness = proof
        .witness_parameter
        .and_then(|t| finite_witness(&numerator, scale, t));
    Ok(GlobalPassivityReport {
        verdict: proof.verdict,
        subdivisions: proof.subdivisions,
        violation_frequency_rad_s: witness,
        unresolved: proof.unresolved,
    })
}

fn relative_perturbation(original: &RationalModel, repaired: &RationalModel) -> f64 {
    let (a, b) = (residue_coords(original), residue_coords(repaired));
    let scale = a.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 {
        return 0.0;
    } // A zero-residue original cannot be active.
    let base = fs_math::det::sqrt(a.iter().map(|x| (x / scale).powi(2)).sum());
    let delta: Vec<_> = a
        .iter()
        .zip(&b)
        .map(|(x, y)| y / scale - x / scale)
        .collect();
    let change_scale = delta.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if change_scale == 0.0 {
        return 0.0;
    }
    if !change_scale.is_finite() {
        return f64::INFINITY;
    }
    let change = fs_math::det::sqrt(delta.iter().map(|x| (x / change_scale).powi(2)).sum());
    change_scale * (change / base)
}

/// Residue-only repair with global re-certification after every QP step.
/// Poles, d and e remain fixed. A proved passive original is returned unchanged;
/// an unresolved proof refuses instead of adding an artificial loss floor.
/// The existing sampled/eigen diagnostic repair remains available separately.
///
/// # Errors
/// Model preconditions as in check_global_passivity, invalid fit band,
/// GlobalUnresolved for an incomplete proof, or RepairFailed after 12 steps.
pub fn repair_passivity_global(
    model: &RationalModel,
    band: (f64, f64),
    budget: GlobalPassivityBudget,
) -> Result<(RationalModel, GlobalRepairReport), PassivityError> {
    if !band.0.is_finite()
        || !band.1.is_finite()
        || band.0 < 0.0
        || band.1 <= band.0
        || !(band.1 * 10.0).is_finite()
    {
        return Err(PassivityError::InvalidModel {
            what: "global repair needs a finite nonnegative increasing band",
        });
    }
    let mut current = model.clone();
    let mut constraints = Vec::new();
    let mut last_kkt = 0.0;
    for round in 0..=12 {
        let certificate = check_global_passivity(&current, budget)?;
        match certificate.verdict {
            GlobalPassivityVerdict::Certified => {
                let relative_perturbation = relative_perturbation(model, &current);
                if !relative_perturbation.is_finite() {
                    return Err(PassivityError::GlobalUnresolved {
                        reason: GlobalPassivityUnresolved::NumericRange,
                    });
                }
                return Ok((
                    current,
                    GlobalRepairReport {
                        rounds: round,
                        relative_perturbation,
                        kkt_residual: last_kkt,
                        certificate,
                    },
                ));
            }
            GlobalPassivityVerdict::Unresolved => {
                return Err(PassivityError::GlobalUnresolved {
                    reason: certificate
                        .unresolved
                        .unwrap_or(GlobalPassivityUnresolved::NumericRange),
                });
            }
            GlobalPassivityVerdict::Violated => {}
        }
        if round == 12 {
            break;
        }
        let Some(witness) = certificate.violation_frequency_rad_s else {
            break;
        };
        // The global witness is indispensable: unlike the sampled band it
        // includes DC, narrow dips and both unbounded frequency tails.
        if !constraints.contains(&witness) {
            constraints.push(witness);
        }
        let mut worst = (0.0, 0.0);
        for w in build_grid(&current, band) {
            let re = current.eval_iw(w).re;
            if re.is_finite() && re < worst.0 {
                worst = (re, w);
            }
        }
        if worst.0 < 0.0 && !constraints.contains(&worst.1) {
            constraints.push(worst.1);
        }
        let dip = constraints
            .iter()
            .map(|&w| model.eval_iw(w).re)
            .fold(f64::INFINITY, f64::min);
        if !dip.is_finite() || dip >= 0.0 {
            break;
        }
        let margin = -dip * 0.05 * (round as f64 + 1.0);
        let (candidate, kkt) = qp_step(model, &constraints, margin);
        current = candidate;
        last_kkt = kkt;
    }
    Err(PassivityError::RepairFailed { rounds: 12 })
}

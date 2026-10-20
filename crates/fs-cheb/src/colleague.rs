//! Colleague-matrix rootfinding (bead kw89): the Chebyshev companion.
//! The v1 subdivision scanner (`Cheb1::roots`) has no generic
//! even-multiplicity guarantee (there may be no sign change to see, and the
//! retained fixture is missed). In exact algebra, roots of the represented polynomial are
//! eigenvalues of its colleague matrix; this finite-precision implementation
//! returns approximate candidates after a documented filtering and clustering
//! policy, so it does not claim complete enumeration or multiplicity recovery.
//! Eigenvalues come from the fs-la complex nonsymmetric stack. Independently,
//! `certified_roots_bounded` searches the full real domain using interval
//! Newton and retains exhaustion, unresolved regions, and polynomial class.
//! A local unique-root proof does not establish complete enumeration.

use crate::budget::{ChebBudget, ChebError, admit_root_certification};
#[cfg(test)]
use crate::normalize_coefficients_exact;
use crate::{Cheb1, affine_from_reference};
use fs_ivl::{
    DEFAULT_MAX_ROOT_BOXES, Interval, RootBox, RootSearchConfig, RootSearchTermination,
    krawczyk_step, newton_roots_with_checkpoint,
};
use fs_la::eigen_complex::eig;
use fs_math::c64::C64;

/// Root-filter policy (all tolerances RELATIVE to the domain size
/// where lengths are involved).
#[derive(Debug, Clone, Copy)]
pub struct ColleaguePolicy {
    /// Trailing coefficients below `trim_rel`·max|c| are trimmed
    /// before building the matrix (degree deflation).
    pub trim_rel: f64,
    /// Eigenvalues with |Im| above this are discarded as complex.
    pub im_tol: f64,
    /// Reference-domain slack outside [−1, 1] still accepted (roots
    /// on the boundary jitter by rounding), clamped back in.
    pub domain_slack: f64,
    /// Roots closer than this (reference domain) merge into one
    /// reported position (multiple eigenvalues of a multiple root).
    pub cluster_tol: f64,
}

impl Default for ColleaguePolicy {
    fn default() -> Self {
        ColleaguePolicy {
            trim_rel: 1e-13,
            im_tol: 1e-8,
            domain_slack: 1e-8,
            // A double root's eigenvalue pair splits at the √ε scale
            // (measured 5e-9 on the battery's fixture); the cluster
            // width must sit above it or the pair reports twice.
            cluster_tol: 1e-6,
        }
    }
}

fn assert_policy(policy: ColleaguePolicy) {
    assert!(
        policy.trim_rel.is_finite() && policy.trim_rel >= 0.0,
        "colleague trim tolerance must be finite and non-negative"
    );
    assert!(
        policy.im_tol.is_finite() && policy.im_tol >= 0.0,
        "colleague imaginary tolerance must be finite and non-negative"
    );
    assert!(
        policy.domain_slack.is_finite() && policy.domain_slack >= 0.0,
        "colleague domain slack must be finite and non-negative"
    );
    assert!(
        policy.cluster_tol.is_finite() && policy.cluster_tol > 0.0,
        "colleague cluster tolerance must be finite and positive"
    );
}

/// Evaluate `-numerator / (2*denominator)` without committing to an
/// overflowing/underflowing doubled denominator or to an overflowing full
/// ratio whose half is representable. The historical expression remains the
/// first path so ordinary inputs keep their established bits.
fn negative_half_ratio(numerator: f64, denominator: f64) -> f64 {
    assert!(
        numerator.is_finite() && denominator.is_finite() && denominator != 0.0,
        "finite numerator and finite non-zero denominator required"
    );
    let mut finite_zero_fallback = None;
    let doubled_denominator = 2.0 * denominator;
    if doubled_denominator.is_finite() && doubled_denominator != 0.0 {
        let ordinary = -numerator / doubled_denominator;
        if ordinary.is_finite() && (ordinary != 0.0 || numerator == 0.0) {
            return ordinary;
        }
        if ordinary.is_finite() {
            finite_zero_fallback = Some(ordinary);
        }
    }

    let full_ratio = numerator / denominator;
    if full_ratio.is_finite() {
        let ratio_first = -0.5 * full_ratio;
        if ratio_first.is_finite() && (ratio_first != 0.0 || numerator == 0.0) {
            return ratio_first;
        }
        if ratio_first.is_finite() {
            finite_zero_fallback = Some(ratio_first);
        }
    }

    let half_numerator = 0.5 * numerator;
    if half_numerator != 0.0 || numerator == 0.0 {
        let numerator_first = -half_numerator / denominator;
        if numerator_first.is_finite() && (numerator_first != 0.0 || numerator == 0.0) {
            return numerator_first;
        }
        if numerator_first.is_finite() {
            finite_zero_fallback = Some(numerator_first);
        }
    }
    if finite_zero_fallback.is_some() {
        panic!("colleague matrix half-ratio underflows a non-zero coefficient");
    }
    panic!("colleague matrix half-ratio is not representable as finite f64")
}

/// Greatest finite power of two no larger than `value`.
///
/// Dividing by this scale is exact whenever the result remains representable;
/// the caller verifies that property coefficient by coefficient.  An arbitrary
/// `cmax` divisor can round a non-zero coefficient to zero and thereby change
/// the represented polynomial before the eigensolve even begins.
fn power_of_two_scale(value: f64) -> f64 {
    assert!(value.is_finite() && value > 0.0);
    let bits = value.to_bits();
    let exponent = (bits >> 52) & 0x7ff;
    if exponent == 0 {
        let mantissa = bits & ((1_u64 << 52) - 1);
        debug_assert_ne!(mantissa, 0);
        let highest_bit = 63_u32 - mantissa.leading_zeros();
        f64::from_bits(1_u64 << highest_bit)
    } else {
        f64::from_bits(exponent << 52)
    }
}

/// Policy-filtered real-root candidates from the colleague matrix,
/// deduplicated per the policy. Even-multiplicity candidates can be found, but
/// the imaginary/domain filters and `cluster_tol` deliberately make this an
/// approximate candidate API rather than a completeness theorem.
///
/// # Panics
/// If the trimmed polynomial is constant (no roots to define), if an exact
/// power-of-two normalization or matrix half-ratio cannot preserve every
/// non-zero coefficient in finite `f64`, or if the eigensolver fails to
/// converge (typed failure surfaced as a panic at fixture scale).
#[must_use]
pub fn colleague_roots(p: &Cheb1, policy: ColleaguePolicy) -> Vec<f64> {
    assert_policy(policy);
    // Normalize the stored coefficients before halving c₀. Halving a
    // subnormal stored c₀ first can round it to zero and change the roots.
    // The normalization scale MUST be a power of two and every non-zero
    // coefficient MUST survive an exact round trip.  Refusing an exponent span
    // that cannot fit in one f64 matrix is preferable to silently solving a
    // different polynomial.
    // Cheb1 stores the Σ′ convention (f = c₀/2 + Σ c_k T_k), while
    // the colleague algebra wants the true constant term.
    let mut coeffs = p.coeffs().to_vec();
    // The policy is defined on the mathematical (unprimed) series, whose
    // constant coefficient is c0/2. Compute that scale before normalization;
    // using the stored c0 would move the trailing-trim boundary by a factor of
    // two whenever the constant term dominates.
    let cmax = coeffs.iter().enumerate().fold(0.0f64, |m, (index, &c)| {
        m.max(if index == 0 { 0.5 * c.abs() } else { c.abs() })
    });
    assert!(
        cmax.is_finite() && cmax > 0.0,
        "finite non-zero polynomial required for colleague matrix"
    );
    let scale = power_of_two_scale(cmax);
    for coefficient in &mut coeffs {
        let original = *coefficient;
        let scaled = original / scale;
        assert!(
            scaled.is_finite() && (original == 0.0 || scaled != 0.0) && scaled * scale == original,
            "colleague normalization cannot preserve the coefficient exponent span"
        );
        *coefficient = scaled;
    }
    let normalized_cmax = cmax / scale;
    let stored_c0 = coeffs[0];
    coeffs[0] *= 0.5;
    assert!(
        stored_c0 == 0.0 || coeffs[0] != 0.0,
        "colleague normalization cannot represent the unprimed constant coefficient"
    );
    let mut n = coeffs.len() - 1;
    while n > 0 && coeffs[n].abs() / normalized_cmax <= policy.trim_rel {
        n -= 1;
    }
    assert!(n >= 1, "constant polynomial has no roots to define");
    // Colleague matrix (n×n) for Σ_{k=0}^{n} a_k T_k:
    // row 0:      x T_0 = T_1                → [0, 1, 0, ...]
    // row k:      x T_k = (T_{k−1}+T_{k+1})/2 → [.., ½, 0, ½, ..]
    // row n−1:    coefficient-loaded: −a_j/(2 a_n) + ½·δ_{j,n−2}.
    let mut m = vec![C64::ZERO; n * n];
    let set = |m: &mut Vec<C64>, r: usize, c: usize, v: f64| {
        m[r * n + c] = C64::new(v, 0.0);
    };
    if n == 1 {
        // a_0 + a_1 T_1 = 0 → x = −a_0/a_1.
        let entry = -coeffs[0] / coeffs[1];
        assert!(
            entry.is_finite(),
            "linear colleague matrix entry is not representable after coefficient normalization"
        );
        set(&mut m, 0, 0, entry);
    } else {
        set(&mut m, 0, 1, 1.0);
        for r in 1..n - 1 {
            set(&mut m, r, r - 1, 0.5);
            set(&mut m, r, r + 1, 0.5);
        }
        let an = coeffs[n];
        for (j, &coeff) in coeffs.iter().enumerate().take(n) {
            let mut v = negative_half_ratio(coeff, an);
            if j == n - 2 {
                let combined = v + 0.5;
                assert!(
                    v == 0.0 || combined != 0.5,
                    "colleague recurrence addition absorbed a non-zero coefficient term"
                );
                assert!(
                    combined != v,
                    "colleague recurrence addition absorbed its non-zero half term"
                );
                v = combined;
            }
            assert!(
                v.is_finite(),
                "colleague matrix entry is not representable after coefficient normalization"
            );
            set(&mut m, n - 1, j, v);
        }
    }
    let eigs = eig(&m, n).expect("colleague eigensolve converges");
    let mut roots: Vec<f64> = eigs
        .into_iter()
        .filter(|l| l.im.abs() <= policy.im_tol)
        .map(|l| l.re)
        .filter(|&t| t >= -1.0 - policy.domain_slack && t <= 1.0 + policy.domain_slack)
        .map(|t| t.clamp(-1.0, 1.0))
        .collect();
    roots.sort_by(f64::total_cmp);
    // Cluster dedupe (multiple roots surface as eigenvalue clusters).
    let mut out: Vec<f64> = Vec::new();
    for t in roots {
        if out.last().is_none_or(|&prev| t - prev > policy.cluster_tol) {
            out.push(t);
        }
    }
    // Map reference → domain.
    let (a, b) = p.domain();
    out.iter()
        .map(|&t| affine_from_reference(t, a, b))
        .collect()
}

/// Clenshaw evaluation of a Chebyshev series in INTERVAL arithmetic —
/// a rigorous enclosure of p over the reference-domain interval `t`.
#[must_use]
pub fn clenshaw_interval(coeffs: &[f64], t: Interval) -> Interval {
    assert!(
        !coeffs.is_empty() && coeffs.iter().all(|c| c.is_finite()),
        "interval Clenshaw needs finite coefficients"
    );
    let two_t = Interval::point(2.0) * t;
    let mut b1 = Interval::point(0.0);
    let mut b2 = Interval::point(0.0);
    for &coefficient in coeffs.iter().rev() {
        let b0 = Interval::point(coefficient) + two_t * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    b1 - t * b2
}

/// Evaluate a plain-coefficient Chebyshev series and its exact-real
/// derivative together in interval arithmetic. Differentiating rounded f64
/// coefficients first would make interval Newton certify a nearby polynomial
/// instead of the polynomial supplied by the caller.
#[cfg(test)]
fn clenshaw_interval_with_derivative(coeffs: &[Interval], t: Interval) -> (Interval, Interval) {
    clenshaw_interval_checkpointed(coeffs, t, &mut || Ok::<(), core::convert::Infallible>(()))
        .unwrap()
}

/// Certified root enclosures over the interpolant's domain via
/// interval Newton (fs-ivl): every returned [`RootBox::Certified`]
/// PROVES a unique root inside; [`RootBox::Possible`] boxes are
/// honest ambiguity (multiple roots, endpoint roots, `min_width` termination,
/// or the interval search's finite box budget can all land here). `min_width`
/// is a dimensionless width in the reference coordinate, not a physical
/// length.
///
/// This list-only convenience API refuses box-budget exhaustion and the zero
/// polynomial. Use [`certified_roots_bounded`] or [`certified_roots_budgeted`]
/// to retain explicit classification, work, and partial-search results.
#[must_use]
pub fn certified_roots(p: &Cheb1, min_width: f64) -> Vec<RootBox> {
    let report = certified_roots_bounded(
        p,
        RootSearchConfig {
            min_width,
            max_boxes: DEFAULT_MAX_ROOT_BOXES,
        },
        &ChebBudget::default(),
    )
    .unwrap_or_else(|error| panic!("invalid certified root search: {error}"));
    assert!(
        report.polynomial_class != RootPolynomialClass::IdenticallyZero,
        "the identically zero polynomial has a continuum of roots; use certified_roots_bounded"
    );
    assert!(
        report.domain_exhausted,
        "certified root search exhausted its box budget; use certified_roots_bounded to retain unresolved regions"
    );
    report.roots
}

/// Exact classification of the stored finite binary64 coefficients, including
/// trailing-zero degree drops. This does not classify a sampled source function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootPolynomialClass {
    /// Cancellation occurred before coefficient classification finished.
    Unclassified,
    /// Every point of the closed physical domain is a root.
    IdenticallyZero,
    /// The represented polynomial has no roots.
    NonzeroConstant,
    /// A nonconstant real polynomial with the given effective degree.
    Nonconstant {
        /// Highest index with a nonzero stored coefficient.
        degree: usize,
    },
}

/// Full-domain real root-search result. Physical enclosures and reference
/// widths have deliberately separate fields. No complex-root or algebraic
/// multiplicity claim is made.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedRootReport {
    /// Closed physical domain, with both endpoints included.
    pub domain: Interval,
    /// Classification of the exact represented polynomial.
    pub polynomial_class: RootPolynomialClass,
    /// Disjoint physical root clusters in increasing order. Overlapping
    /// affine images are merged and downgraded to `Possible`.
    pub roots: Vec<RootBox>,
    /// Number of reference-domain subdivision boxes evaluated.
    pub boxes_examined: usize,
    /// Whether subdivision finished without budget exhaustion or cancellation.
    /// `true` can still accompany `Possible` boxes and does not establish counts.
    pub domain_exhausted: bool,
    /// Evaluated box widths in dimensionless reference coordinates, in visit order.
    pub reference_width_trajectory: Vec<f64>,
}

impl CertifiedRootReport {
    /// A checked finite enumeration, available only after a completed search
    /// whose disjoint boxes each prove existence and uniqueness. An empty slice
    /// then proves no real roots in the declared domain. The zero polynomial
    /// and any unresolved region cannot pass this projection.
    #[must_use]
    pub fn isolated_unique_roots(&self) -> Option<&[RootBox]> {
        (self.domain_exhausted
            && matches!(
                self.polynomial_class,
                RootPolynomialClass::NonzeroConstant | RootPolynomialClass::Nonconstant { .. }
            )
            && self.roots.iter().all(RootBox::is_certified))
        .then_some(self.roots.as_slice())
    }
}

/// Root enclosures and the reason the bounded computation stopped.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedRootOutcome<E> {
    /// Physical root coverage, including unresolved work.
    pub report: CertifiedRootReport,
    /// Exhaustion and cancellation are distinct from a completed subdivision.
    pub termination: RootSearchTermination<E>,
}

/// Search the full real domain with caller-supplied box, work, and memory caps.
/// `config.min_width` is in the dimensionless reference coordinate.
/// Unlike [`certified_roots`], this API represents zero polynomials and partial
/// searches without a panic or a list-only projection.
pub fn certified_roots_bounded(
    p: &Cheb1,
    config: RootSearchConfig,
    budget: &ChebBudget,
) -> Result<CertifiedRootReport, ChebError> {
    let outcome = certified_roots_with_checkpoint(p, config, budget, &mut || {
        Ok::<(), core::convert::Infallible>(())
    })?;
    Ok(outcome.report)
}

/// The explicit [`fs_exec::Cx`] entry point for bounded root certification.
pub fn certified_roots_budgeted(
    p: &Cheb1,
    config: RootSearchConfig,
    budget: &ChebBudget,
    cx: &fs_exec::Cx<'_>,
) -> Result<CertifiedRootOutcome<fs_exec::Cancelled>, ChebError> {
    certified_roots_with_checkpoint(p, config, budget, &mut || cx.checkpoint())
}

/// Cancellable interval-Newton search. Polls occur during coefficient admission
/// and every 64 Clenshaw coefficients, before Newton evaluations, and before
/// publication. Cancellation retains conservative physical coverage, never an
/// apparently complete prefix of roots. `Err` is an admission/numerical refusal;
/// checkpoint errors are retained in the returned `Cancelled` outcome.
pub fn certified_roots_with_checkpoint<E, C>(
    p: &Cheb1,
    config: RootSearchConfig,
    budget: &ChebBudget,
    checkpoint: &mut C,
) -> Result<CertifiedRootOutcome<E>, ChebError>
where
    C: FnMut() -> Result<(), E>,
{
    if !config.min_width.is_finite() || config.min_width <= 0.0 {
        return Err(ChebError::Shape {
            what: "certified root reference min_width must be finite and positive",
        });
    }
    admit_root_certification(p.coeffs().len(), config.max_boxes, budget)?;
    let (a, b) = p.domain();
    let domain = Interval::new(a, b);
    let mut report = CertifiedRootReport {
        domain,
        polynomial_class: RootPolynomialClass::Unclassified,
        roots: vec![RootBox::Possible(domain)],
        boxes_examined: 0,
        domain_exhausted: false,
        reference_width_trajectory: Vec::new(),
    };
    let mut last_nonzero = None;
    let mut cmax = 0.0_f64;
    for (index, &coefficient) in p.coeffs().iter().enumerate() {
        if index.is_multiple_of(64) {
            if let Err(error) = checkpoint() {
                return Ok(CertifiedRootOutcome {
                    report,
                    termination: RootSearchTermination::Cancelled(error),
                });
            }
        }
        if coefficient != 0.0 {
            last_nonzero = Some(index);
        }
        cmax = cmax.max(coefficient.abs());
    }
    report.polynomial_class = match last_nonzero {
        None => RootPolynomialClass::IdenticallyZero,
        Some(0) => RootPolynomialClass::NonzeroConstant,
        Some(degree) => RootPolynomialClass::Nonconstant { degree },
    };
    if !matches!(
        report.polynomial_class,
        RootPolynomialClass::Nonconstant { .. }
    ) {
        if let Err(error) = checkpoint() {
            return Ok(CertifiedRootOutcome {
                report,
                termination: RootSearchTermination::Cancelled(error),
            });
        }
        if report.polynomial_class == RootPolynomialClass::NonzeroConstant {
            report.roots.clear();
        }
        report.domain_exhausted = true;
        return Ok(CertifiedRootOutcome {
            report,
            termination: RootSearchTermination::Complete,
        });
    }
    let scale = crate::normalization_power_of_two(cmax);
    let mut coeffs = Vec::new();
    coeffs
        .try_reserve_exact(p.coeffs().len())
        .map_err(|_| ChebError::Numerical {
            what: "root interval coefficient allocation failed",
        })?;
    for (index, &original) in p.coeffs().iter().enumerate() {
        if index.is_multiple_of(64) {
            if let Err(error) = checkpoint() {
                return Ok(CertifiedRootOutcome {
                    report,
                    termination: RootSearchTermination::Cancelled(error),
                });
            }
        }
        let normalized = original / scale;
        if !normalized.is_finite() || normalized * scale != original {
            return Err(ChebError::Numerical {
                what: "root certification cannot normalize coefficients without information loss",
            });
        }
        let value = Interval::point(normalized);
        coeffs.push(if index == 0 {
            Interval::point(0.5) * value
        } else {
            value
        });
    }
    // The interval solver's callbacks return enclosures, so a checkpoint
    // failure inside Clenshaw supplies WHOLE temporarily. The next solver
    // checkpoint consumes the original error before accepting any contraction.
    let callback = core::cell::RefCell::new(checkpoint);
    let pending = core::cell::RefCell::new(None);
    let poll = || match pending.borrow_mut().take() {
        Some(error) => Err(error),
        None => (*callback.borrow_mut())(),
    };
    let evaluate = |t| match clenshaw_interval_checkpointed(&coeffs, t, &mut || poll()) {
        Ok(value) => value,
        Err(error) => {
            *pending.borrow_mut() = Some(error);
            (Interval::WHOLE, Interval::WHOLE)
        }
    };
    let f = |t| evaluate(t).0;
    let fp = |t| evaluate(t).1;
    let outcome =
        newton_roots_with_checkpoint(&f, &fp, Interval::new(-1.0, 1.0), config, &mut || poll())
            .map_err(|_| ChebError::Shape {
                what: "invalid interval root-search configuration",
            })?;
    report.boxes_examined = outcome.report.boxes_examined;
    report.reference_width_trajectory = outcome.report.width_trajectory;
    report.domain_exhausted = outcome.report.complete;
    let mut termination = outcome.termination;
    let mut roots = outcome.report.roots;
    for (index, root) in roots.iter_mut().enumerate() {
        if index.is_multiple_of(64)
            && !matches!(termination, RootSearchTermination::Cancelled(_))
            && let Err(error) = poll()
        {
            termination = RootSearchTermination::Cancelled(error);
            report.domain_exhausted = false;
        }
        let mapped = map_interval_to_domain(root.interval(), a, b);
        let mut certified = false;
        if root.is_certified() && !matches!(termination, RootSearchTermination::Cancelled(_)) {
            if let Err(error) = poll() {
                termination = RootSearchTermination::Cancelled(error);
                report.domain_exhausted = false;
            } else {
                let back = map_interval_to_reference(mapped, a, b);
                certified = back.lo().is_finite()
                    && back.hi().is_finite()
                    && krawczyk_step(&f, &fp, back).is_some_and(|(_, strict)| strict);
                if let Some(error) = pending.borrow_mut().take() {
                    termination = RootSearchTermination::Cancelled(error);
                    report.domain_exhausted = false;
                    certified = false;
                }
            }
        }
        *root = if certified {
            RootBox::Certified(mapped)
        } else {
            RootBox::Possible(mapped)
        };
    }
    // Outward affine rounding can make distinct reference boxes overlap,
    // especially on a domain containing only a few representable numbers.
    // Such a physical cluster has no unique-root count until independently resolved.
    roots.sort_by(|x, y| {
        x.interval()
            .lo()
            .total_cmp(&y.interval().lo())
            .then_with(|| x.interval().hi().total_cmp(&y.interval().hi()))
    });
    let mut kept = 0;
    for index in 0..roots.len() {
        let root = roots[index];
        if kept > 0 && roots[kept - 1].interval().hi() >= root.interval().lo() {
            roots[kept - 1] = RootBox::Possible(roots[kept - 1].interval().hull(root.interval()));
        } else {
            roots[kept] = root;
            kept += 1;
        }
    }
    roots.truncate(kept);
    report.roots = roots;
    if !matches!(termination, RootSearchTermination::Cancelled(_)) {
        if let Err(error) = poll() {
            termination = RootSearchTermination::Cancelled(error);
            report.domain_exhausted = false;
        }
    }
    Ok(CertifiedRootOutcome {
        report,
        termination,
    })
}

fn clenshaw_interval_checkpointed<E, C>(
    coeffs: &[Interval],
    t: Interval,
    checkpoint: &mut C,
) -> Result<(Interval, Interval), E>
where
    C: FnMut() -> Result<(), E>,
{
    let two = Interval::point(2.0);
    let two_t = two * t;
    let mut b1 = Interval::point(0.0);
    let mut b2 = Interval::point(0.0);
    let mut db1 = Interval::point(0.0);
    let mut db2 = Interval::point(0.0);
    for (index, &coefficient) in coeffs.iter().rev().enumerate() {
        if index.is_multiple_of(64) {
            checkpoint()?;
        }
        let b0 = coefficient + two_t * b1 - b2;
        let db0 = two * b1 + two_t * db1 - db2;
        b2 = b1;
        b1 = b0;
        db2 = db1;
        db1 = db0;
    }
    Ok((b1 - t * b2, db1 - b2 - t * db2))
}

fn map_interval_to_domain(iv: Interval, a: f64, b: f64) -> Interval {
    let lo_image = map_reference_point_to_domain(iv.lo(), a, b);
    let hi_image = map_reference_point_to_domain(iv.hi(), a, b);
    Interval::new(lo_image.lo().max(a), hi_image.hi().min(b))
}

fn affine_center_radius_intervals(a: f64, b: f64) -> (Interval, Interval) {
    let half = Interval::point(0.5);
    let a_interval = Interval::point(a);
    let b_interval = Interval::point(b);
    let raw_center = half * a_interval + half * b_interval;
    let raw_radius = half * b_interval - half * a_interval;
    // Mathematical center and radius of two finite ordered endpoints remain
    // inside these hard bounds. Intersecting with them prevents an outward
    // nudge at f64::MAX from turning an otherwise finite affine certificate
    // into an artificial infinity.
    let center = Interval::new(raw_center.lo().max(a), raw_center.hi().min(b));
    let radius = Interval::new(raw_radius.lo().max(0.0), raw_radius.hi().min(f64::MAX));
    (center, radius)
}

fn map_reference_point_to_domain(t: f64, a: f64, b: f64) -> Interval {
    if t == -1.0 {
        return Interval::point(a);
    }
    if t == 1.0 {
        return Interval::point(b);
    }
    let (center, radius) = affine_center_radius_intervals(a, b);
    center + Interval::point(t) * radius
}

fn map_domain_point_to_reference(x: f64, a: f64, b: f64) -> Interval {
    if x == a {
        return Interval::point(-1.0);
    }
    if x == b {
        return Interval::point(1.0);
    }
    let (center, radius) = affine_center_radius_intervals(a, b);
    (Interval::point(x) - center) / radius
}

fn map_interval_to_reference(iv: Interval, a: f64, b: f64) -> Interval {
    let lo_image = map_domain_point_to_reference(iv.lo(), a, b);
    let hi_image = map_domain_point_to_reference(iv.hi(), a, b);
    Interval::new(lo_image.lo().max(-1.0), hi_image.hi().min(1.0))
}

#[cfg(test)]
mod affine_certificate_tests {
    use super::*;

    #[test]
    fn interval_clenshaw_checks_inside_long_coefficient_arrays() {
        let coefficients = vec![Interval::point(1.0); 129];
        let mut polls = 0;
        let result =
            clenshaw_interval_checkpointed(&coefficients, Interval::point(0.5), &mut || {
                polls += 1;
                if polls == 2 {
                    Err("coefficient tile cancelled")
                } else {
                    Ok(())
                }
            });
        assert_eq!(result, Err("coefficient tile cancelled"));
        assert_eq!(polls, 2);
    }

    #[test]
    fn colleague_half_ratio_avoids_representable_intermediate_overflow() {
        let denominator = f64::MIN_POSITIVE / 7.0;
        assert!((1.0 / denominator).is_infinite());
        let ratio = negative_half_ratio(1.0, denominator);
        assert!(ratio.is_finite() && ratio < 0.0);
        assert!((ratio.abs() * denominator - 0.5).abs() <= f64::EPSILON);

        let min_subnormal = f64::from_bits(1);
        let tiny_numerator = (min_subnormal * f64::MAX) * 1.2;
        assert_eq!((tiny_numerator / f64::MAX).to_bits(), 1);
        assert_eq!((-0.5 * (tiny_numerator / f64::MAX)).to_bits(), 1_u64 << 63);
        assert_eq!(
            negative_half_ratio(tiny_numerator, f64::MAX).to_bits(),
            (1_u64 << 63) | 1,
            "numerator-first ordering must retain a representable negative min-subnormal"
        );

        let unrepresentable = std::panic::catch_unwind(|| {
            negative_half_ratio(f64::from_bits(1), 1.0);
        });
        assert!(
            unrepresentable.is_err(),
            "a non-zero matrix coefficient must never be silently rounded to zero"
        );
    }

    #[test]
    fn colleague_recurrence_refuses_absorbed_nonzero_term() {
        let polynomial = Cheb1::from_coeffs(
            -1.0,
            1.0,
            vec![0.0, fs_math::det::powi(2.0, -100), 0.0, 1.0],
        );
        let result = std::panic::catch_unwind(|| {
            let _ = colleague_roots(
                &polynomial,
                ColleaguePolicy {
                    trim_rel: 0.0,
                    ..ColleaguePolicy::default()
                },
            );
        });
        assert!(
            result.is_err(),
            "a recurrence row must not silently discard a non-zero coefficient"
        );
    }

    #[test]
    fn colleague_refuses_an_exponent_span_that_would_drop_a_term() {
        let min_subnormal = f64::from_bits(1);
        let amplitude = f64::MAX / 2.0;
        // f(t) = min_subnormal + amplitude*t^3.  Expressing t^3 in the
        // Chebyshev basis produces a huge middle/high pair while the constant
        // remains non-zero.  Scaling by the former arbitrary cmax rounded that
        // constant to zero and changed the polynomial.  No single f64
        // colleague matrix can preserve this exponent span, so the approximate
        // API must refuse rather than return roots of the truncated polynomial.
        let p = Cheb1::from_coeffs(
            -1.0,
            1.0,
            vec![2.0 * min_subnormal, 0.75 * amplitude, 0.0, 0.25 * amplitude],
        );
        let refused = std::panic::catch_unwind(|| colleague_roots(&p, ColleaguePolicy::default()));
        assert!(
            refused.is_err(),
            "unsupported exponent span must fail closed"
        );
    }

    #[test]
    fn directed_affine_image_contains_adversarial_exact_value() {
        let a = -1.278_585_201_636_465_8e17;
        let b = 2.708_877_815_465_865_6e17;
        let t = -0.327_713_252_609_748_8;
        let image = map_reference_point_to_domain(t, a, b);
        // Nearest f64 to the independently evaluated exact-real affine image
        // 6177406941685617.699570... . The former scalar + one-ULP scheme
        // ended at 6177406941685617 and excluded the real value.
        assert!(image.contains(6_177_406_941_685_618.0), "{image:?}");
        assert!(image.lo().is_finite() && image.hi().is_finite());
    }

    #[test]
    fn finite_domain_images_never_expand_to_infinity() {
        for t in [-1.0, -0.5, 0.0, 0.5, 1.0] {
            let image = map_interval_to_domain(Interval::point(t), -f64::MAX, f64::MAX);
            assert!(image.lo().is_finite() && image.hi().is_finite());
            if t > -1.0 && t < 1.0 {
                assert!(image.lo() > -f64::MAX && image.hi() < f64::MAX);
                assert!(image.contains(t * f64::MAX));
            }
        }

        let asymmetric = map_interval_to_domain(Interval::point(0.25), -f64::MAX, f64::MAX / 2.0);
        assert!(asymmetric.lo().is_finite() && asymmetric.hi().is_finite());
        assert!(asymmetric.contains(affine_from_reference(0.25, -f64::MAX, f64::MAX / 2.0,)));

        let adjacent_hi = f64::from_bits(1.0_f64.to_bits() + 1);
        let adjacent = map_interval_to_domain(Interval::point(0.0), 1.0, adjacent_hi);
        assert!(adjacent.lo().is_finite() && adjacent.hi().is_finite());
        assert!(adjacent.contains(f64::midpoint(1.0, adjacent_hi)));
    }

    #[test]
    fn interval_clenshaw_derivative_encloses_the_exact_series_derivative() {
        // p(t) = 3 + 2 T1(t) + 4 T2(t), so p'(t) = 2 + 16t.
        let coefficients = [3.0, 2.0, 4.0].map(Interval::point);
        let point = Interval::point(0.25);
        let (value, derivative) = clenshaw_interval_with_derivative(&coefficients, point);
        assert!(value.contains(0.0));
        assert!(derivative.contains(6.0));
    }

    #[test]
    fn exact_normalization_refuses_silent_subnormal_information_loss() {
        let mut coefficients = [f64::from_bits(3), 2.0];
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            normalize_coefficients_exact(&mut coefficients, "test normalization");
        }));
        assert!(result.is_err());
    }
}

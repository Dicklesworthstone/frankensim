//! The scalar abstraction and B-spline basis machinery, written ONCE and
//! instantiated at both `f64` (fast path) and [`crate::Rat`] (the exact
//! path the refinement-exactness claims are proved in).

use crate::NurbsError;
use crate::rat::Rat;

/// Defensive ceiling on Cox-de Boor triangular work in the legacy unbudgeted
/// basis API. Typed caller budgets and cancellation belong to its successor.
pub const BASIS_MAX_WORK_UNITS: u128 = 16_777_216;

// Conservative price for finite/order/run/multiplicity/clamping validation of
// one public knot entry. This intentionally overcounts the simple comparisons:
// admission must happen before any full scan, not after three nominally cheap
// passes over caller-mutable storage.
const KNOT_VALIDATION_WORK_PER_ENTRY: u128 = 16;

/// The field the spline algebra runs over.
pub trait Scalar:
    Copy
    + PartialEq
    + PartialOrd
    + core::fmt::Debug
    + core::ops::Add<Output = Self>
    + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self>
    + core::ops::Div<Output = Self>
    + core::ops::Neg<Output = Self>
{
    /// Additive identity.
    fn zero() -> Self;
    /// Multiplicative identity.
    fn one() -> Self;
    /// Lift a small integer.
    fn from_int(v: i64) -> Self;
    /// Whether this value belongs to the finite numeric domain admitted by
    /// spline structure. Exact scalar domains return `true`; floating and dual
    /// domains must reject NaN and infinities.
    fn is_finite(self) -> bool;
    /// Whether a positive rational weight is numerically representable without
    /// an immediate zero-denominator hazard. Exact domains may accept every
    /// positive value. Floating domains must reject subnormal weights because
    /// multiplying them by an ordinary basis value can underflow to zero even
    /// when every source value is finite.
    fn is_admissible_weight(self) -> bool {
        self.is_finite() && self > Self::zero()
    }
    /// Whether dividing a homogeneous numerator by an admitted weight stays in
    /// this scalar's finite Cartesian domain. Exact domains can answer without
    /// performing a potentially huge intermediate division.
    fn quotient_is_finite(self, denominator: Self) -> bool {
        (self / denominator).is_finite()
    }
}

impl Scalar for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn from_int(v: i64) -> Self {
        #[allow(clippy::cast_precision_loss)]
        {
            v as f64
        }
    }
    fn is_finite(self) -> bool {
        self.is_finite()
    }
    fn is_admissible_weight(self) -> bool {
        self.is_normal() && self > 0.0
    }
}

impl Scalar for Rat {
    fn zero() -> Self {
        Rat::int(0)
    }
    fn one() -> Self {
        Rat::int(1)
    }
    fn from_int(v: i64) -> Self {
        Rat::int(v)
    }
    fn is_finite(self) -> bool {
        true
    }
    fn quotient_is_finite(self, _denominator: Self) -> bool {
        true
    }
}

/// A clamped knot vector for degree-p splines.
#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector<S: Scalar> {
    /// Non-decreasing knots (first/last with multiplicity p+1).
    pub knots: Vec<S>,
    /// Polynomial degree.
    pub degree: usize,
}

impl<S: Scalar> KnotVector<S> {
    fn admitted_scan_work(&self, include_span_search: bool) -> Result<u128, NurbsError> {
        let knot_count = self.knots.len() as u128;
        let mut work = knot_count
            .checked_mul(KNOT_VALIDATION_WORK_PER_ENTRY)
            .and_then(|work| work.checked_add(self.degree as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "knot-scan work accounting overflows u128".to_string(),
            })?;
        if include_span_search {
            work = work
                .checked_add(self.control_count() as u128)
                .ok_or_else(|| NurbsError::Domain {
                    what: "span-search work accounting overflows u128".to_string(),
                })?;
        }
        Ok(work)
    }

    fn enforce_work(units: u128, operation: &str) -> Result<(), NurbsError> {
        if units > BASIS_MAX_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "{operation} requests {units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
                ),
            });
        }
        Ok(())
    }

    fn validated_domain(&self) -> (S, S) {
        (
            self.knots[self.degree],
            self.knots[self.knots.len() - 1 - self.degree],
        )
    }

    fn span_after_validation(&self, t: S) -> Result<usize, NurbsError> {
        let (lo, hi) = self.validated_domain();
        if !t.is_finite() || t < lo || t > hi {
            return Err(NurbsError::Domain {
                what: format!("parameter {t:?} outside {lo:?}..{hi:?}"),
            });
        }
        let n = self.control_count() - 1;
        if t == hi {
            // Validation guarantees at least one non-empty span, so this walk
            // cannot underflow.
            let mut s = n;
            while self.knots[s] == self.knots[s + 1] {
                s -= 1;
            }
            return Ok(s);
        }
        let mut span = self.degree;
        while span < n && self.knots[span + 1] <= t {
            span += 1;
        }
        Ok(span)
    }

    /// Revalidate the current public fields before any indexing algorithm uses
    /// them. This is intentionally allocation-free: callers can mutate the
    /// early-stage public representation after construction, and safe public
    /// input must become a structured refusal rather than an index panic.
    pub(crate) fn validate_live(&self) -> Result<(), NurbsError> {
        let endpoint_multiplicity =
            self.degree
                .checked_add(1)
                .ok_or_else(|| NurbsError::Structure {
                    what: format!("degree {} overflows knot-count arithmetic", self.degree),
                })?;
        let minimum_knots =
            endpoint_multiplicity
                .checked_mul(2)
                .ok_or_else(|| NurbsError::Structure {
                    what: format!("degree {} overflows knot-count arithmetic", self.degree),
                })?;
        if self.degree == 0 || self.knots.len() < minimum_knots {
            return Err(NurbsError::Structure {
                what: format!(
                    "degree {} needs at least {minimum_knots} knots, got {}",
                    self.degree,
                    self.knots.len()
                ),
            });
        }
        if self.knots.iter().copied().any(|knot| !knot.is_finite()) {
            return Err(NurbsError::Structure {
                what: "knots must be finite".to_string(),
            });
        }
        if self.knots.windows(2).any(|window| window[1] < window[0]) {
            return Err(NurbsError::Structure {
                what: "knots must be non-decreasing".to_string(),
            });
        }
        let mut run_start = 0usize;
        while run_start < self.knots.len() {
            let mut run_end = run_start + 1;
            while run_end < self.knots.len() && self.knots[run_end] == self.knots[run_start] {
                run_end += 1;
            }
            let multiplicity = run_end - run_start;
            let endpoint = run_start == 0 || run_end == self.knots.len();
            if (endpoint && multiplicity != endpoint_multiplicity)
                || (!endpoint && multiplicity > endpoint_multiplicity)
            {
                return Err(NurbsError::Structure {
                    what: format!(
                        "knot multiplicity {multiplicity} is invalid for degree {}",
                        self.degree
                    ),
                });
            }
            run_start = run_end;
        }
        for offset in 0..self.degree {
            if self.knots[offset + 1] != self.knots[0]
                || self.knots[self.knots.len() - 2 - offset] != self.knots[self.knots.len() - 1]
            {
                return Err(NurbsError::Structure {
                    what: "knot vector must be clamped (end multiplicity degree+1)".to_string(),
                });
            }
        }
        if self.knots[self.degree] == self.knots[self.knots.len() - 1 - self.degree] {
            return Err(NurbsError::Structure {
                what: "knot vector has an empty parametric domain (lo == hi)".to_string(),
            });
        }
        Ok(())
    }

    /// Validate and construct.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on ordering/clamping defects.
    pub fn new(knots: Vec<S>, degree: usize) -> Result<Self, NurbsError> {
        let endpoint_multiplicity = degree.checked_add(1).ok_or_else(|| NurbsError::Structure {
            what: format!("degree {degree} overflows knot-count arithmetic"),
        })?;
        let minimum_knots =
            endpoint_multiplicity
                .checked_mul(2)
                .ok_or_else(|| NurbsError::Structure {
                    what: format!("degree {degree} overflows knot-count arithmetic"),
                })?;
        if degree == 0 || knots.len() < minimum_knots {
            return Err(NurbsError::Structure {
                what: format!(
                    "degree {degree} needs at least {} knots, got {}",
                    minimum_knots,
                    knots.len()
                ),
            });
        }
        if knots.iter().copied().any(|knot| !knot.is_finite()) {
            return Err(NurbsError::Structure {
                what: "knots must be finite".to_string(),
            });
        }
        if knots.windows(2).any(|w| w[1] < w[0]) {
            return Err(NurbsError::Structure {
                what: "knots must be non-decreasing".to_string(),
            });
        }
        let mut run_start = 0usize;
        while run_start < knots.len() {
            let mut run_end = run_start + 1;
            while run_end < knots.len() && knots[run_end] == knots[run_start] {
                run_end += 1;
            }
            let multiplicity = run_end - run_start;
            let endpoint = run_start == 0 || run_end == knots.len();
            if (endpoint && multiplicity != endpoint_multiplicity)
                || (!endpoint && multiplicity > endpoint_multiplicity)
            {
                return Err(NurbsError::Structure {
                    what: format!(
                        "knot multiplicity {multiplicity} is invalid for degree {degree}"
                    ),
                });
            }
            run_start = run_end;
        }
        for k in 0..degree {
            if knots[k + 1] != knots[0] || knots[knots.len() - 2 - k] != knots[knots.len() - 1] {
                return Err(NurbsError::Structure {
                    what: "knot vector must be clamped (end multiplicity degree+1)".to_string(),
                });
            }
        }
        // The parametric domain [knots[degree], knots[len-1-degree]] must be
        // non-empty. An all-equal (zero-width) knot vector passes every check
        // above but has lo == hi, and `span(hi)`'s degenerate-span walk-back
        // would decrement past 0 (usize underflow → panic).
        if knots[degree] == knots[knots.len() - 1 - degree] {
            return Err(NurbsError::Structure {
                what: "knot vector has an empty parametric domain (lo == hi)".to_string(),
            });
        }
        Ok(KnotVector { knots, degree })
    }

    /// Number of basis functions / control points.
    #[must_use]
    pub fn control_count(&self) -> usize {
        self.knots
            .len()
            .checked_sub(self.degree)
            .and_then(|count| count.checked_sub(1))
            .unwrap_or(0)
    }

    /// The parametric domain `[u_min, u_max]`, after revalidating the mutable
    /// early-stage public representation.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the knot vector was mutated into an
    /// invalid shape; [`NurbsError::Domain`] when the defensive live-scan work
    /// ceiling is exceeded.
    pub fn domain(&self) -> Result<(S, S), NurbsError> {
        Self::enforce_work(self.admitted_scan_work(false)?, "knot-domain validation")?;
        self.validate_live()?;
        Ok(self.validated_domain())
    }

    /// The knot span index containing `t` (Piegl–Tiller A2.1 semantics;
    /// the end parameter maps into the last non-empty span).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or when defensive
    /// live-validation/span-search work admission refuses the request.
    pub fn span(&self, t: S) -> Result<usize, NurbsError> {
        Self::enforce_work(self.admitted_scan_work(true)?, "knot-span evaluation")?;
        self.validate_live()?;
        self.span_after_validation(t)
    }

    /// All nonzero basis-function values at `t` (Cox–de Boor triangle,
    /// Piegl–Tiller A2.2): `N_{span-p..=span, p}(t)`.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or when defensive
    /// validation, span-search, triangular-work, or allocation admission
    /// refuses the request.
    pub fn basis(&self, t: S) -> Result<(usize, Vec<S>), NurbsError> {
        let p = self.degree;
        let order = p.checked_add(1).ok_or_else(|| NurbsError::Domain {
            what: "basis order overflows usize".to_string(),
        })?;
        let triangular_work = (p as u128)
            .checked_mul(order as u128)
            .map(|product| product / 2)
            .and_then(|work| work.checked_add(order as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "basis-work accounting overflows u128".to_string(),
            })?;
        let total_work = self
            .admitted_scan_work(true)?
            .checked_add(triangular_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "basis total-work accounting overflows u128".to_string(),
            })?;
        Self::enforce_work(total_work, "basis evaluation")?;
        self.validate_live()?;
        let span = self.span_after_validation(t)?;
        let mut n = Vec::new();
        let mut left = Vec::new();
        let mut right = Vec::new();
        for (buffer, stage) in [
            (&mut n, "values"),
            (&mut left, "left workspace"),
            (&mut right, "right workspace"),
        ] {
            buffer
                .try_reserve_exact(order)
                .map_err(|_| NurbsError::Domain {
                    what: format!("basis {stage} allocation was refused"),
                })?;
            buffer.resize(order, S::zero());
        }
        n[0] = S::one();
        for j in 1..=p {
            left[j] = t - self.knots[span + 1 - j];
            right[j] = self.knots[span + j] - t;
            let mut saved = S::zero();
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = n[r] / denom;
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        if n.iter().copied().any(|value| !value.is_finite()) {
            return Err(NurbsError::Domain {
                what: format!("basis evaluation at {t:?} left the finite numeric domain"),
            });
        }
        Ok((span, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_domain_knot_vector_is_rejected_not_paniced() {
        // Regression: an all-equal knot vector passes the count / monotone /
        // clamped checks but has an empty domain (lo == hi). `span(hi)` then
        // underflowed its degenerate-span walk-back (usize `0 - 1`). Must refuse
        // at construction instead.
        assert!(KnotVector::new(vec![5.0f64; 6], 2).is_err());
        assert!(KnotVector::new(vec![0.0f64, 0.0, 0.0, 0.0], 1).is_err());
        // A proper clamped vector with a real domain builds and resolves the
        // upper-endpoint span without panicking.
        let kv = KnotVector::new(vec![0.0f64, 0.0, 0.0, 1.0, 1.0, 1.0], 2).expect("valid");
        assert_eq!(kv.span(1.0).expect("hi is in domain"), 2);
    }

    #[test]
    fn excessive_endpoint_and_interior_multiplicity_are_rejected() {
        assert!(KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0], 1).is_err());
        assert!(KnotVector::new(vec![0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0], 1).is_err());
    }

    #[test]
    fn non_finite_query_parameter_is_rejected() {
        let kv = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("valid line knots");
        assert!(kv.span(f64::NAN).is_err());
        assert!(kv.basis(f64::INFINITY).is_err());
    }

    #[test]
    fn domain_and_basis_fail_closed_on_mutation_and_quadratic_work() {
        let mut malformed = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("valid line knots");
        malformed.knots.clear();
        assert!(
            malformed.domain().is_err(),
            "safe public mutation must not turn domain access into an indexing panic"
        );

        let degree = 6_000usize;
        let mut knots = vec![0.0; degree + 1];
        knots.extend(vec![1.0; degree + 1]);
        let high_degree = KnotVector::new(knots, degree).expect("large but structurally valid");
        assert!(
            high_degree.basis(0.5).is_err(),
            "quadratic Cox-de Boor work must be refused before entering billions of iterations"
        );

        let interior_count = 1_000_000usize;
        let mut many_knots = Vec::with_capacity(interior_count + 4);
        many_knots.extend([0.0, 0.0]);
        for index in 1..=interior_count {
            #[allow(clippy::cast_precision_loss)]
            many_knots.push(index as f64 / (interior_count + 1) as f64);
        }
        many_knots.extend([1.0, 1.0]);
        let low_degree_many_spans = KnotVector {
            knots: many_knots,
            degree: 1,
        };
        assert!(
            low_degree_many_spans.basis(0.5).is_err(),
            "low polynomial degree must not bypass full knot-scan admission"
        );
        assert!(
            low_degree_many_spans.span(0.5).is_err(),
            "the public span search must share the defensive scan ceiling"
        );
    }
}

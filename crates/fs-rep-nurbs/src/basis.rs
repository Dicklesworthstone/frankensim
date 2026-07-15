//! The scalar abstraction and B-spline basis machinery, written ONCE and
//! instantiated at both `f64` (fast path) and [`crate::Rat`] (the exact
//! path the refinement-exactness claims are proved in).

use crate::NurbsError;
use crate::rat::Rat;
use fs_exec::Cx;

/// Defensive ceiling on Cox-de Boor triangular work in the measured basis
/// APIs. Typed caller budgets and cancellation-aware owning admission remain
/// successor work.
pub const BASIS_MAX_WORK_UNITS: u128 = 16_777_216;

// An admitted basis evaluation polls after at most this many logical span,
// initialization, triangle, or finite-check operations. The caller still owns
// request -> drain -> finalize; this primitive only observes the shared gate.
const BASIS_CANCELLATION_STRIDE: usize = 64;

// Conservative price for finite/order/run/multiplicity/clamping validation of
// one public knot entry. This intentionally overcounts the simple comparisons:
// admission must happen before any full scan, not after three nominally cheap
// passes over untrusted storage.
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
///
/// The representation is sealed after construction. Callers can inspect it
/// through [`Self::knots`] and [`Self::degree`], but cannot mutate around a
/// successful validation:
///
/// ```compile_fail
/// use fs_rep_nurbs::KnotVector;
/// let mut knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
/// knots.knots.clear();
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector<S: Scalar> {
    /// Non-decreasing knots (first/last with multiplicity p+1).
    pub(crate) knots: Vec<S>,
    /// Polynomial degree.
    pub(crate) degree: usize,
}

/// A validate-once borrow of one exact immutable knot-vector snapshot.
///
/// The borrow is the authority: safe Rust cannot mutate or replace the source
/// while this view is live, so no content hash or recomputed token is needed to
/// detect stale structure.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedKnotVector<'a, S: Scalar> {
    inner: &'a KnotVector<S>,
}

/// Transactional terminal state of one cancellation-aware basis evaluation.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum BasisRun<S: Scalar> {
    /// The complete nonzero basis row is safe to publish.
    Complete {
        /// Knot span containing the requested parameter.
        span: usize,
        /// `degree + 1` nonzero basis values in ascending control index.
        values: Vec<S>,
    },
    /// Cancellation was observed at a bounded poll; no partial row escapes.
    Cancelled,
}

fn basis_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < BASIS_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

impl<S: Scalar> KnotVector<S> {
    fn validation_work_for(knot_count: usize, degree: usize) -> Result<u128, NurbsError> {
        (knot_count as u128)
            .checked_mul(KNOT_VALIDATION_WORK_PER_ENTRY)
            .and_then(|work| work.checked_add(degree as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "knot-scan work accounting overflows u128".to_string(),
            })
    }

    pub(crate) fn validation_work(&self) -> Result<u128, NurbsError> {
        Self::validation_work_for(self.knots.len(), self.degree)
    }

    fn span_search_work(&self) -> u128 {
        self.control_count() as u128
    }

    fn basis_operation_work(&self) -> Result<u128, NurbsError> {
        let order = self
            .degree
            .checked_add(1)
            .ok_or_else(|| NurbsError::Domain {
                what: "basis order overflows usize".to_string(),
            })?;
        (self.degree as u128)
            .checked_mul(order as u128)
            .map(|product| product / 2)
            .and_then(|work| work.checked_add(order as u128))
            .and_then(|work| work.checked_add(self.span_search_work()))
            .ok_or_else(|| NurbsError::Domain {
                what: "basis-work accounting overflows u128".to_string(),
            })
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

    /// Validate the sealed fields before any indexing algorithm uses them.
    /// This remains allocation-free defense for crate-internal construction;
    /// public callers cannot mutate the representation after construction.
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
    /// [`NurbsError::Structure`] on ordering/clamping defects, or
    /// [`NurbsError::Domain`] when validation work exceeds the defensive cap.
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
        let validation_work = Self::validation_work_for(knots.len(), degree)?;
        Self::enforce_work(validation_work, "knot-vector construction")?;
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

    /// Borrow the immutable knot entries.
    #[must_use]
    pub fn knots(&self) -> &[S] {
        &self.knots
    }

    /// Polynomial degree.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.degree
    }

    /// Refuse an invalid parameter in constant time before any structural scan
    /// or allocation-bearing basis work. Public sealed owners establish these
    /// endpoint invariants at construction; the checked indexing keeps this
    /// defensive for crate-internal candidates as well.
    pub(crate) fn preflight_parameter(
        &self,
        parameter: S,
        operation: &str,
    ) -> Result<(), NurbsError> {
        if !parameter.is_finite() {
            return Err(NurbsError::Domain {
                what: format!("{operation} parameter must be finite"),
            });
        }
        let Some(&lo) = self.knots.get(self.degree) else {
            return Err(NurbsError::Structure {
                what: format!("{operation} knot vector has no lower domain endpoint"),
            });
        };
        let Some(hi_index) = self
            .knots
            .len()
            .checked_sub(1)
            .and_then(|last| last.checked_sub(self.degree))
        else {
            return Err(NurbsError::Structure {
                what: format!("{operation} knot vector has no upper domain endpoint"),
            });
        };
        let Some(&hi) = self.knots.get(hi_index) else {
            return Err(NurbsError::Structure {
                what: format!("{operation} knot vector has no upper domain endpoint"),
            });
        };
        if !lo.is_finite() || !hi.is_finite() || lo >= hi {
            return Err(NurbsError::Structure {
                what: format!("{operation} knot vector has an invalid parametric domain"),
            });
        }
        if parameter < lo || parameter > hi {
            return Err(NurbsError::Domain {
                what: format!("{operation} parameter {parameter:?} outside {lo:?}..{hi:?}"),
            });
        }
        Ok(())
    }

    /// Fallibly copy this sealed knot vector without revalidating unchanged
    /// entries.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when the destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let mut knots = Vec::new();
        knots
            .try_reserve_exact(self.knots.len())
            .map_err(|_| NurbsError::Domain {
                what: "knot-vector copy allocation was refused".to_string(),
            })?;
        knots.extend_from_slice(&self.knots);
        Ok(KnotVector {
            knots,
            degree: self.degree,
        })
    }

    /// Validate this exact immutable snapshot once and bind the proof to its
    /// borrow lifetime.
    ///
    /// # Errors
    /// Returns a structured refusal when validation work exceeds the defensive
    /// ceiling or the representation is malformed.
    pub fn admit(&self) -> Result<AdmittedKnotVector<'_, S>, NurbsError> {
        Self::enforce_work(self.validation_work()?, "knot-vector admission")?;
        self.validate_live()?;
        Ok(self.admitted_after_validation())
    }

    pub(crate) const fn admitted_after_validation(&self) -> AdmittedKnotVector<'_, S> {
        AdmittedKnotVector { inner: self }
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

    /// The parametric domain `[u_min, u_max]`, after structural admission.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the knot vector was mutated into an
    /// invalid shape; [`NurbsError::Domain`] when the defensive live-scan work
    /// ceiling is exceeded.
    pub fn domain(&self) -> Result<(S, S), NurbsError> {
        self.admit().map(|admitted| admitted.domain())
    }

    /// The knot span index containing `t` (Piegl–Tiller A2.1 semantics;
    /// the end parameter maps into the last non-empty span).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or when defensive
    /// live-validation/span-search work admission refuses the request.
    pub fn span(&self, t: S) -> Result<usize, NurbsError> {
        let total_work = self
            .validation_work()?
            .checked_add(self.span_search_work())
            .ok_or_else(|| NurbsError::Domain {
                what: "knot-span work accounting overflows u128".to_string(),
            })?;
        Self::enforce_work(total_work, "knot-span evaluation")?;
        self.validate_live()?;
        self.admitted_after_validation().span_after_preflight(t)
    }

    /// All nonzero basis-function values at `t` (Cox–de Boor triangle,
    /// Piegl–Tiller A2.2): `N_{span-p..=span, p}(t)`.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or when defensive
    /// validation, span-search, triangular-work, or allocation admission
    /// refuses the request.
    pub fn basis(&self, t: S) -> Result<(usize, Vec<S>), NurbsError> {
        let total_work = self
            .validation_work()?
            .checked_add(self.basis_operation_work()?)
            .ok_or_else(|| NurbsError::Domain {
                what: "basis total-work accounting overflows u128".to_string(),
            })?;
        Self::enforce_work(total_work, "basis evaluation")?;
        self.validate_live()?;
        self.admitted_after_validation().basis_after_preflight(t)
    }
}

impl<'a, S: Scalar> AdmittedKnotVector<'a, S> {
    /// The exact immutable source bound to this view.
    #[must_use]
    pub const fn source(&self) -> &'a KnotVector<S> {
        self.inner
    }

    /// Borrow the validated knot entries.
    #[must_use]
    pub fn knots(&self) -> &'a [S] {
        self.inner.knots()
    }

    /// Polynomial degree.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.inner.degree()
    }

    /// Number of basis functions / control points.
    #[must_use]
    pub fn control_count(&self) -> usize {
        self.inner.control_count()
    }

    /// The already-validated parametric domain.
    #[must_use]
    pub fn domain(&self) -> (S, S) {
        self.inner.validated_domain()
    }

    /// Resolve a knot span without rescanning structure.
    ///
    /// # Errors
    /// Returns a structured refusal for out-of-domain parameters or excessive
    /// span-search work.
    pub fn span(&self, t: S) -> Result<usize, NurbsError> {
        KnotVector::<S>::enforce_work(
            self.inner.span_search_work(),
            "admitted knot-span evaluation",
        )?;
        self.span_after_preflight(t)
    }

    fn span_after_preflight(&self, t: S) -> Result<usize, NurbsError> {
        self.inner.span_after_validation(t)
    }

    /// Evaluate all nonzero basis values without rescanning the sealed source.
    ///
    /// # Errors
    /// Returns a structured refusal for domain, work, allocation, or finite
    /// arithmetic failures.
    pub fn basis(&self, t: S) -> Result<(usize, Vec<S>), NurbsError> {
        KnotVector::<S>::enforce_work(
            self.inner.basis_operation_work()?,
            "admitted basis evaluation",
        )?;
        self.basis_after_preflight(t)
    }

    /// Evaluate all nonzero basis values with bounded cancellation polling.
    ///
    /// This method is deliberately available only on an admitted immutable
    /// snapshot: it does not imply that the owning [`KnotVector::admit`]
    /// structural scan is cancellation-aware. Cancellation is checked after
    /// cheap request/work admission, at fixed work strides, before each
    /// allocation, and immediately before publication. [`BasisRun::Cancelled`]
    /// carries no partial basis row. The caller remains responsible for
    /// draining and finalizing any surrounding executor scope; `Cx` budgets are
    /// not consumed by this measured primitive.
    ///
    /// # Errors
    /// Returns a structured refusal for domain, work, allocation, or finite
    /// arithmetic failures that precede successful publication.
    pub fn basis_with_cx(&self, t: S, cx: &Cx<'_>) -> Result<BasisRun<S>, NurbsError> {
        KnotVector::<S>::enforce_work(
            self.inner.basis_operation_work()?,
            "admitted basis evaluation",
        )?;
        self.basis_after_preflight_with_poll(t, || cx.checkpoint().is_err())
    }

    fn basis_after_preflight(&self, t: S) -> Result<(usize, Vec<S>), NurbsError> {
        match self.basis_after_preflight_with_poll(t, || false)? {
            BasisRun::Complete { span, values } => Ok((span, values)),
            BasisRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling basis evaluation observed cancellation".to_string(),
            }),
        }
    }

    fn basis_after_preflight_with_poll(
        &self,
        t: S,
        mut should_cancel: impl FnMut() -> bool,
    ) -> Result<BasisRun<S>, NurbsError> {
        let inner = self.inner;
        let p = inner.degree;
        let order = p.checked_add(1).ok_or_else(|| NurbsError::Domain {
            what: "basis order overflows usize".to_string(),
        })?;
        let (lo, hi) = inner.validated_domain();
        if !t.is_finite() || t < lo || t > hi {
            return Err(NurbsError::Domain {
                what: format!("parameter {t:?} outside {lo:?}..{hi:?}"),
            });
        }
        if should_cancel() {
            return Ok(BasisRun::Cancelled);
        }

        let mut operations_since_poll = 0usize;
        let n_last = inner.control_count() - 1;
        let span = if t == hi {
            // Admission guarantees a non-empty span, so this cannot underflow.
            let mut span = n_last;
            while inner.knots[span] == inner.knots[span + 1] {
                span -= 1;
                if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(BasisRun::Cancelled);
                }
            }
            span
        } else {
            let mut span = p;
            while span < n_last && inner.knots[span + 1] <= t {
                span += 1;
                if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(BasisRun::Cancelled);
                }
            }
            span
        };

        let mut n = Vec::new();
        let mut left = Vec::new();
        let mut right = Vec::new();
        for (buffer, stage) in [
            (&mut n, "values"),
            (&mut left, "left workspace"),
            (&mut right, "right workspace"),
        ] {
            if should_cancel() {
                return Ok(BasisRun::Cancelled);
            }
            operations_since_poll = 0;
            buffer
                .try_reserve_exact(order)
                .map_err(|_| NurbsError::Domain {
                    what: format!("basis {stage} allocation was refused"),
                })?;
            for _ in 0..order {
                buffer.push(S::zero());
                if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(BasisRun::Cancelled);
                }
            }
        }
        n[0] = S::one();
        for j in 1..=p {
            left[j] = t - inner.knots[span + 1 - j];
            right[j] = inner.knots[span + j] - t;
            if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(BasisRun::Cancelled);
            }
            let mut saved = S::zero();
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = n[r] / denom;
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
                if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(BasisRun::Cancelled);
                }
            }
            n[j] = saved;
        }
        for value in &n {
            if !value.is_finite() {
                return Err(NurbsError::Domain {
                    what: format!("basis evaluation at {t:?} left the finite numeric domain"),
                });
            }
            if basis_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(BasisRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(BasisRun::Cancelled);
        }
        Ok(BasisRun::Complete { span, values: n })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_basis_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0xB4515,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn cancellation_fixture() -> KnotVector<f64> {
        let degree = 16;
        let mut knots = vec![0.0; degree + 1];
        knots.extend(vec![1.0; degree + 1]);
        KnotVector::new(knots, degree).expect("cancellation fixture")
    }

    #[test]
    fn construction_admits_work_before_the_first_knot_scan() {
        let exact_cap_count = 1_048_575usize;
        assert_eq!(
            KnotVector::<f64>::validation_work_for(exact_cap_count, 16).expect("exact-cap work"),
            BASIS_MAX_WORK_UNITS
        );
        assert_eq!(
            KnotVector::<f64>::validation_work_for(exact_cap_count, 17).expect("cap-plus-one work"),
            BASIS_MAX_WORK_UNITS + 1
        );

        let over_cap = KnotVector::new(vec![f64::NAN; exact_cap_count], 17)
            .expect_err("cap-plus-one construction must be refused");
        assert!(
            matches!(over_cap, NurbsError::Domain { .. }),
            "work refusal must precede the non-finite scalar scan"
        );

        let exact_cap = KnotVector::new(vec![f64::NAN; exact_cap_count], 16)
            .expect_err("the exact-cap request reaches finite-value validation");
        assert!(
            matches!(exact_cap, NurbsError::Structure { .. }),
            "an exact-cap request must reach semantic validation"
        );
    }

    #[test]
    fn admitted_basis_cancellation_is_transactional_and_preserves_request_precedence() {
        let knots = cancellation_fixture();
        let admitted = knots.admit().expect("admitted fixture");
        with_basis_cx(true, |cx| {
            assert_eq!(
                admitted
                    .basis_with_cx(0.5, cx)
                    .expect("valid request reaches cancellation"),
                BasisRun::Cancelled
            );
            assert!(
                matches!(
                    admitted.basis_with_cx(f64::NAN, cx),
                    Err(NurbsError::Domain { .. })
                ),
                "constant-time request validation must precede cancellation"
            );
        });
    }

    #[test]
    fn admitted_basis_cancellation_replays_at_a_fixed_poll_ordinal() {
        let knots = cancellation_fixture();
        let admitted = knots.admit().expect("admitted fixture");
        let run = || {
            let mut polls = 0usize;
            let outcome = admitted
                .basis_after_preflight_with_poll(0.5, || {
                    polls += 1;
                    polls == 6
                })
                .expect("valid basis arithmetic");
            (outcome, polls)
        };
        let first = run();
        let second = run();
        assert_eq!(first, second);
        assert_eq!(first, (BasisRun::Cancelled, 6));
    }

    #[test]
    fn admitted_basis_with_cx_matches_legacy_exactly() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.25, 0.75, 1.0, 1.0, 1.0], 2)
            .expect("quadratic multispan knots");
        let admitted = knots.admit().expect("admitted fixture");
        let legacy = admitted.basis(0.6).expect("legacy basis");
        with_basis_cx(false, |cx| {
            let run = admitted.basis_with_cx(0.6, cx).expect("cancellable basis");
            assert_eq!(
                run,
                BasisRun::Complete {
                    span: legacy.0,
                    values: legacy.1,
                }
            );
        });
    }

    #[test]
    fn admitted_basis_final_checkpoint_gates_publication() {
        let knots = cancellation_fixture();
        let admitted = knots.admit().expect("admitted fixture");
        let mut total_polls = 0usize;
        let complete = admitted
            .basis_after_preflight_with_poll(0.5, || {
                total_polls += 1;
                false
            })
            .expect("healthy basis run");
        assert!(matches!(complete, BasisRun::Complete { .. }));

        let mut replay_polls = 0usize;
        let cancelled = admitted
            .basis_after_preflight_with_poll(0.5, || {
                replay_polls += 1;
                replay_polls == total_polls
            })
            .expect("final-checkpoint cancellation");
        assert_eq!(cancelled, BasisRun::Cancelled);
        assert_eq!(replay_polls, total_polls);
    }

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
    fn domain_and_basis_fail_closed_on_internal_corruption_and_quadratic_work() {
        let mut malformed = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("valid line knots");
        malformed.knots.clear();
        assert!(
            malformed.domain().is_err(),
            "crate-internal corruption must not turn domain access into an indexing panic"
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

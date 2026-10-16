//! Decision-contextual useful-bound outcomes (E09.2).
//!
//! A rigorous enclosure and a useful engineering bound are different claims.
//! This module preserves the enclosure while making failure of a
//! caller-declared usefulness criterion a typed, absorbing result.

use core::fmt;

use crate::ClaimClass;

/// Semantic version of the typed useful-bound outcome.
pub const USEFUL_BOUND_SCHEMA_VERSION: u16 = 1;
/// Maximum UTF-8 bytes in a decision-context or unit label.
pub const MAX_USEFUL_BOUND_FIELD_BYTES: usize = 512;

/// Closed reason that a valid enclosure cannot support the requested decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NoUsefulBoundCause {
    /// The declared prediction horizon exceeds the admitted useful horizon.
    HorizonTooLong,
    /// Propagated sensitivity or interval width blew up.
    LipschitzBlowup,
    /// The computation left the declared model or chart domain.
    DomainExit,
    /// The work budget ended before a useful enclosure was established.
    BudgetExhausted,
}

impl NoUsefulBoundCause {
    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::HorizonTooLong => "horizon-too-long",
            Self::LipschitzBlowup => "lipschitz-blowup",
            Self::DomainExit => "domain-exit",
            Self::BudgetExhausted => "budget-exhausted",
        }
    }
}

/// Ordered extended-real interval retained by the useful-bound protocol.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundInterval {
    lower: f64,
    upper: f64,
}

impl BoundInterval {
    /// Validate an enclosure. Ordered infinite endpoints are allowed, but a
    /// degenerate point at either infinity is not an enclosure of a real value.
    pub fn try_new(lower: f64, upper: f64) -> Result<Self, UsefulBoundError> {
        if lower.is_nan() || upper.is_nan() {
            return Err(UsefulBoundError::NaNInterval);
        }
        if lower > upper {
            return Err(UsefulBoundError::InvertedInterval { lower, upper });
        }
        if lower.is_infinite() && lower.to_bits() == upper.to_bits() {
            return Err(UsefulBoundError::PointAtInfinity);
        }
        Ok(Self { lower, upper })
    }

    /// Lower endpoint.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Upper endpoint.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    /// Conservative width. Any non-finite endpoint yields positive infinity.
    #[must_use]
    pub fn width(self) -> f64 {
        if self.lower.is_finite() && self.upper.is_finite() {
            let width = self.upper - self.lower;
            if width.is_finite() {
                next_up_nonnegative(width)
            } else {
                f64::INFINITY
            }
        } else {
            f64::INFINITY
        }
    }
}

/// Caller-declared decision context and maximum useful enclosure width.
#[derive(Debug, Clone, PartialEq)]
pub struct UsefulnessCriterion {
    decision_context: String,
    unit: String,
    max_width: f64,
}

impl UsefulnessCriterion {
    /// Construct one bounded, finite criterion.
    pub fn try_new(
        decision_context: impl Into<String>,
        unit: impl Into<String>,
        max_width: f64,
    ) -> Result<Self, UsefulBoundError> {
        let decision_context = validated_field("decision_context", decision_context.into())?;
        let unit = validated_field("unit", unit.into())?;
        if !(max_width.is_finite() && max_width > 0.0) {
            return Err(UsefulBoundError::InvalidMaximumWidth { max_width });
        }
        Ok(Self {
            decision_context,
            unit,
            max_width,
        })
    }

    /// Stable human or machine decision-context label.
    #[must_use]
    pub fn decision_context(&self) -> &str {
        &self.decision_context
    }

    /// Unit of the interval width and threshold.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// Largest width the caller declared useful.
    #[must_use]
    pub const fn max_width(&self) -> f64 {
        self.max_width
    }
}

/// Enclosure that satisfies its retained caller-declared criterion.
#[derive(Debug, Clone, PartialEq)]
pub struct Bound {
    interval: BoundInterval,
    criterion: UsefulnessCriterion,
}

impl Bound {
    /// Retained enclosure.
    #[must_use]
    pub const fn interval(&self) -> BoundInterval {
        self.interval
    }

    /// Exact criterion under which the enclosure was useful.
    #[must_use]
    pub const fn criterion(&self) -> &UsefulnessCriterion {
        &self.criterion
    }
}

/// Valid enclosure that cannot support the caller's requested decision.
#[derive(Debug, Clone, PartialEq)]
pub struct NoUsefulBound {
    interval: BoundInterval,
    width_achieved: f64,
    criterion: UsefulnessCriterion,
    cause: NoUsefulBoundCause,
    suggested_reformulation: ClaimClass,
}

impl NoUsefulBound {
    /// Enclosure obtained before the usefulness refusal.
    #[must_use]
    pub const fn interval(&self) -> BoundInterval {
        self.interval
    }

    /// Conservative width actually achieved.
    #[must_use]
    pub const fn width_achieved(&self) -> f64 {
        self.width_achieved
    }

    /// Caller-declared criterion that failed.
    #[must_use]
    pub const fn criterion(&self) -> &UsefulnessCriterion {
        &self.criterion
    }

    /// Closed failure cause.
    #[must_use]
    pub const fn cause(&self) -> NoUsefulBoundCause {
        self.cause
    }

    /// E09 claim class the caller should consider instead.
    #[must_use]
    pub const fn suggested_reformulation(&self) -> ClaimClass {
        self.suggested_reformulation
    }

    /// Deterministic reviewer-facing rendering. This never uses certificate,
    /// verified, compliant, or non-compliant language.
    #[must_use]
    pub fn render_report(&self) -> String {
        format!(
            "outcome=no-useful-bound cause={} interval=[{},{}] width-achieved={} usefulness-threshold={} unit={} decision-context={} suggested-reformulation={}\n",
            self.cause.code(),
            self.interval.lower(),
            self.interval.upper(),
            self.width_achieved,
            self.criterion.max_width(),
            self.criterion.unit(),
            self.criterion.decision_context(),
            self.suggested_reformulation.code(),
        )
    }
}

/// Typed useful enclosure or typed honest refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundOutcome {
    /// The interval satisfies its retained criterion.
    Bound(Bound),
    /// The interval is valid but cannot support the retained criterion.
    NoUsefulBound(NoUsefulBound),
}

impl BoundOutcome {
    /// Classify an enclosure against one caller-declared criterion.
    #[must_use]
    pub fn classify(
        interval: BoundInterval,
        criterion: UsefulnessCriterion,
        cause_if_too_wide: NoUsefulBoundCause,
        suggested_reformulation: ClaimClass,
    ) -> Self {
        if interval.width() <= criterion.max_width {
            Self::Bound(Bound {
                interval,
                criterion,
            })
        } else {
            Self::NoUsefulBound(NoUsefulBound {
                interval,
                width_achieved: interval.width(),
                criterion,
                cause: cause_if_too_wide,
                suggested_reformulation,
            })
        }
    }

    /// Produce an explicit refusal even when the retained interval is narrow.
    ///
    /// Narrowness cannot erase budget exhaustion, domain exit, or an
    /// already-known horizon violation.
    #[must_use]
    pub fn refuse(
        interval: BoundInterval,
        criterion: UsefulnessCriterion,
        cause: NoUsefulBoundCause,
        suggested_reformulation: ClaimClass,
    ) -> Self {
        Self::NoUsefulBound(NoUsefulBound {
            interval,
            width_achieved: interval.width(),
            criterion,
            cause,
            suggested_reformulation,
        })
    }

    /// The useful bound, if and only if this result is decision-usable.
    #[must_use]
    pub const fn bound(&self) -> Option<&Bound> {
        match self {
            Self::Bound(bound) => Some(bound),
            Self::NoUsefulBound(_) => None,
        }
    }

    /// The typed refusal, if present.
    #[must_use]
    pub const fn no_useful_bound(&self) -> Option<&NoUsefulBound> {
        match self {
            Self::Bound(_) => None,
            Self::NoUsefulBound(refusal) => Some(refusal),
        }
    }

    /// Retained interval regardless of usefulness.
    #[must_use]
    pub const fn interval(&self) -> BoundInterval {
        match self {
            Self::Bound(bound) => bound.interval,
            Self::NoUsefulBound(refusal) => refusal.interval,
        }
    }

    /// Retained decision criterion.
    #[must_use]
    pub const fn criterion(&self) -> &UsefulnessCriterion {
        match self {
            Self::Bound(bound) => &bound.criterion,
            Self::NoUsefulBound(refusal) => &refusal.criterion,
        }
    }

    /// Conservatively compose two outcomes. A refusal is absorbing and the
    /// combining closure is never called when either input is refused.
    pub fn compose_absorbing<F>(
        &self,
        other: &Self,
        cause_if_too_wide: NoUsefulBoundCause,
        suggested_reformulation: ClaimClass,
        combine: F,
    ) -> Result<Self, UsefulBoundError>
    where
        F: FnOnce(BoundInterval, BoundInterval) -> Result<BoundInterval, UsefulBoundError>,
    {
        if self.criterion() != other.criterion() {
            return Err(UsefulBoundError::IncompatibleCriteria);
        }
        if let Self::NoUsefulBound(refusal) = self {
            return Ok(Self::NoUsefulBound(refusal.clone()));
        }
        if let Self::NoUsefulBound(refusal) = other {
            return Ok(Self::NoUsefulBound(refusal.clone()));
        }
        let interval = combine(self.interval(), other.interval())?;
        Ok(Self::classify(
            interval,
            self.criterion().clone(),
            cause_if_too_wide,
            suggested_reformulation,
        ))
    }

    /// Deterministic reviewer-facing explanation.
    #[must_use]
    pub fn render_report(&self) -> String {
        match self {
            Self::Bound(bound) => format!(
                "outcome=bound interval=[{},{}] width={} usefulness-threshold={} unit={} decision-context={}\n",
                bound.interval.lower(),
                bound.interval.upper(),
                bound.interval.width(),
                bound.criterion.max_width(),
                bound.criterion.unit(),
                bound.criterion.decision_context(),
            ),
            Self::NoUsefulBound(refusal) => refusal.render_report(),
        }
    }
}

/// Refusal to construct or compose a useful-bound value.
#[derive(Debug, Clone, PartialEq)]
pub enum UsefulBoundError {
    /// An interval endpoint was NaN.
    NaNInterval,
    /// The interval endpoints were out of order.
    InvertedInterval {
        /// Presented lower endpoint.
        lower: f64,
        /// Presented upper endpoint.
        upper: f64,
    },
    /// A point interval at infinity encloses no real value.
    PointAtInfinity,
    /// A decision context or unit field was empty or too large.
    InvalidField {
        /// Field name.
        field: &'static str,
        /// Presented UTF-8 byte count.
        bytes: usize,
    },
    /// The usefulness threshold was not finite and strictly positive.
    InvalidMaximumWidth {
        /// Presented threshold.
        max_width: f64,
    },
    /// Composition attempted to mix different decision contexts or thresholds.
    IncompatibleCriteria,
}

impl fmt::Display for UsefulBoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NaNInterval => write!(f, "bound interval endpoints must not be NaN"),
            Self::InvertedInterval { lower, upper } => {
                write!(f, "bound interval is inverted: [{lower}, {upper}]")
            }
            Self::PointAtInfinity => {
                write!(f, "a point interval at infinity encloses no real value")
            }
            Self::InvalidField { field, bytes } => write!(
                f,
                "{field} must be nonempty and at most {MAX_USEFUL_BOUND_FIELD_BYTES} bytes; found {bytes}"
            ),
            Self::InvalidMaximumWidth { max_width } => write!(
                f,
                "maximum useful width must be finite and positive; found {max_width}"
            ),
            Self::IncompatibleCriteria => {
                write!(
                    f,
                    "cannot compose bounds from different usefulness criteria"
                )
            }
        }
    }
}

impl std::error::Error for UsefulBoundError {}

fn validated_field(field: &'static str, value: String) -> Result<String, UsefulBoundError> {
    let bytes = value.len();
    if value.trim().is_empty()
        || bytes > MAX_USEFUL_BOUND_FIELD_BYTES
        || value.chars().any(char::is_control)
    {
        return Err(UsefulBoundError::InvalidField { field, bytes });
    }
    Ok(value)
}

fn next_up_nonnegative(value: f64) -> f64 {
    if value == 0.0 {
        f64::from_bits(1)
    } else if value.is_finite() {
        f64::from_bits(value.to_bits() + 1)
    } else {
        value
    }
}

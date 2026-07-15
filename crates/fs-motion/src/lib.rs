//! fs-motion — certified rigid motion for the MORPH layer.
//!
//! `fs-ga` owns instantaneous SE(3) motor algebra; nothing previously
//! bound a motor PATH to a chart. This crate provides
//! [`CertifiedMotorTube`] (piecewise Taylor-model enclosures of the
//! motor components with rigorously measured versor-defect bounds),
//! [`MotorPath`] (the point-evaluation view), [`SpacetimeChart`]
//! (moving geometry with frozen-time snapshots and certified
//! time-span field enclosures), [`SweptChart`] (certified implicit
//! infimum bounds), fail-closed [`EnvelopeChart`] characteristic
//! classification, analytic screw and Wankel-pose constructors, and
//! the [`LowerToMotorTube`] builder contract that lets higher layers
//! lower their motions here without upward dependencies.
//!
//! See `CONTRACT.md` for invariants, determinism class, and no-claim
//! boundaries. Beads: `frankensim-ext-motion-motor-tube-c70j` and
//! `frankensim-ext-motion-swept-envelope-c58q`.

#![forbid(unsafe_code)]

pub mod algebra;
pub mod analytic;
pub mod spacetime;
pub mod swept;
pub mod tube;

pub use analytic::{ScrewParams, WankelParams, screw_tube, wankel_tube};
pub use spacetime::{FieldEnclosure, MotionSnapshot, SpacetimeChart};
pub use swept::{
    EnvelopeBranch, EnvelopeBranchClass, EnvelopeChart, EnvelopeConfig, EnvelopeDecision,
    EnvelopeEvidence, EnvelopeOracle, EnvelopeTraceReceipt, EnvelopeTraceStats, ProofState,
    SweepDecision, SweepReceipt, SweptChart, SweptConfig, WankelApexPoint, WankelSealCircle,
    classify_envelope_branch, envelope,
};
pub use tube::{
    BoxActionEnclosure, CertifiedMotorTube, EnclosureClass, LowerToMotorTube, MotorPath,
    MotorTubeSegment, PathSample, PointActionEnclosure,
};

use fs_ivl::TaylorModelError;

/// Typed refusals for motion construction and evaluation.
#[derive(Debug, Clone, PartialEq)]
pub enum MotionError {
    /// A parameter was NaN or infinite.
    NonFiniteInput {
        /// Which parameter family refused.
        what: &'static str,
    },
    /// The time domain is empty, inverted, or non-finite.
    EmptyTimeDomain,
    /// Zero segments requested.
    InvalidSegments,
    /// A component model does not share the multivector's domain and
    /// order.
    MixedModelShape {
        /// Blade index of the offending component.
        blade: usize,
        /// The multivector's order.
        expected_order: usize,
        /// The offered model's order.
        got_order: usize,
    },
    /// Propagated fs-ivl Taylor-model refusal.
    Taylor(TaylorModelError),
    /// The homogeneous weight enclosure contains zero.
    DegenerateWeight {
        /// Weight lower bound.
        lo: f64,
        /// Weight upper bound.
        hi: f64,
    },
    /// Every component midpoint at the sign anchor is below tolerance;
    /// the double-cover branch cannot be fixed deterministically.
    DoubleCoverAmbiguous {
        /// The anchor time.
        at: f64,
    },
    /// Adjacent segments fail the transition test (enclosure overlap
    /// plus positive representative dot product) at a boundary.
    ChartTransition {
        /// The boundary time.
        at: f64,
        /// The representative dot product (NaN when the domains do not
        /// abut or the enclosures do not overlap).
        dot: f64,
    },
    /// A query left the tube's time domain.
    OutOfDomain {
        /// Query lower bound.
        lo: f64,
        /// Query upper bound.
        hi: f64,
        /// Domain lower bound.
        domain_lo: f64,
        /// Domain upper bound.
        domain_hi: f64,
    },
    /// `eval_over` requires an `ExactDistance` base chart.
    UnsupportedBaseClaim,
    /// The base chart's sample certificate is not a rigorous
    /// enclosure.
    UncertifiedBaseSample,
    /// A finite support enclosure is required by the requested operation.
    UnboundedSupport,
    /// A caller-supplied accuracy or work configuration is invalid.
    InvalidConfiguration {
        /// The rejected condition.
        what: &'static str,
    },
    /// Independently certified lower/upper bounds contradicted one another.
    InconsistentEnclosure {
        /// Purported lower bound.
        lower: f64,
        /// Purported upper bound.
        upper: f64,
    },
    /// Declared machine geometry violates a construction precondition.
    InvalidGeometry {
        /// The rejected condition.
        what: &'static str,
    },
    /// A finite PGA point action unexpectedly produced an ideal point.
    PointActionFailed,
    /// Cooperative cancellation was observed.
    Cancelled,
}

impl std::fmt::Display for MotionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MotionError::NonFiniteInput { what } => {
                write!(f, "non-finite {what}")
            }
            MotionError::EmptyTimeDomain => {
                write!(f, "empty, inverted, or non-finite time domain")
            }
            MotionError::InvalidSegments => write!(f, "segment count must be positive"),
            MotionError::MixedModelShape {
                blade,
                expected_order,
                got_order,
            } => write!(
                f,
                "component model at blade {blade} has order {got_order}, expected \
                 {expected_order} on the shared domain"
            ),
            MotionError::Taylor(e) => write!(f, "taylor model refusal: {e}"),
            MotionError::DegenerateWeight { lo, hi } => {
                write!(f, "homogeneous weight enclosure [{lo}, {hi}] contains zero")
            }
            MotionError::DoubleCoverAmbiguous { at } => write!(
                f,
                "double-cover sign is ambiguous at anchor time {at}: every component \
                 midpoint is below tolerance"
            ),
            MotionError::ChartTransition { at, dot } => write!(
                f,
                "chart transition at t = {at} refused (representative dot product {dot}); \
                 adjacent segments must abut, overlap, and agree in double-cover sign"
            ),
            MotionError::OutOfDomain {
                lo,
                hi,
                domain_lo,
                domain_hi,
            } => write!(
                f,
                "query span [{lo}, {hi}] leaves the tube domain [{domain_lo}, {domain_hi}]"
            ),
            MotionError::UnsupportedBaseClaim => write!(
                f,
                "eval_over requires a base chart claiming ExactDistance; other claims \
                 refuse instead of guessing"
            ),
            MotionError::UncertifiedBaseSample => write!(
                f,
                "base chart sample certificate is not a rigorous enclosure"
            ),
            MotionError::UnboundedSupport => {
                write!(f, "operation requires a finite base support enclosure")
            }
            MotionError::InvalidConfiguration { what } => {
                write!(f, "invalid motion configuration: {what}")
            }
            MotionError::InconsistentEnclosure { lower, upper } => write!(
                f,
                "certified infimum bounds are inconsistent: lower {lower} exceeds upper {upper}"
            ),
            MotionError::InvalidGeometry { what } => {
                write!(f, "invalid machine geometry: {what}")
            }
            MotionError::PointActionFailed => {
                write!(f, "finite motor action produced no finite point")
            }
            MotionError::Cancelled => write!(f, "cancelled at a tile boundary"),
        }
    }
}

impl std::error::Error for MotionError {}

impl From<TaylorModelError> for MotionError {
    fn from(e: TaylorModelError) -> Self {
        MotionError::Taylor(e)
    }
}

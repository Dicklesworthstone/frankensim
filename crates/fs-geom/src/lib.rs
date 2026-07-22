//! fs-geom — the Region/Chart abstraction (plan §7.1): the geometry
//! kernel's founding move. Layer: L2.
//!
//! An abstract `Region` — semantically a measurable subset of ℝ³ with
//! piecewise-smooth boundary — is NEVER stored directly. It is PRESENTED
//! through [`Chart`]s: concrete representations answering signed-distance
//! queries with a value, a gradient where one exists, a certified local
//! Lipschitz bound, and a DECLARED error model relative to the abstract
//! region (an [`fs_evidence::NumericalCertificate`]).
//!
//! Three consequences no conventional kernel delivers (plan §7.1):
//! - "the same shape held three ways" is a normal, coherent state:
//!   [`Region`] holds multiple charts with provenance, and AGREEMENT
//!   BETWEEN CHARTS IS A CHECKABLE PROPOSITION
//!   ([`Region::check_agreement`]) with localized diagnostics, not an
//!   assumption;
//! - every conversion's error and authority are explicit and composable:
//!   [`Convert`] always returns evidence feeding the Error Ledger, while only
//!   rigorous global conversions may promote that evidence to
//!   [`fs_evidence::Certified`]
//!   (Decalogue P4);
//! - no chart type is privileged, ever — the Rep Router (a later bead)
//!   picks per OPERATION from declared capabilities.
//!
//! Cancellation: chart evaluation and sampling/build kernels take
//! `&fs_exec::Cx` and poll at bounded work units. Some finite algebraic
//! diagnostics remain synchronous legacy APIs without a `Cx`; their contracts
//! identify that no-claim boundary rather than implying universal P7 coverage.
//!
//! Object safety note: plan Appendix B sketches `Chart { type Param; ... }`.
//! `Region` must hold heterogeneous charts (`Arc<dyn Chart>`), so the
//! design-lever handle lives on the [`DesignChart`] subtrait instead —
//! same contract, object-safe core (fs-xform builds on `DesignChart`).

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;

mod convert;
#[cfg(feature = "derived-geometry")]
pub mod derived;
#[cfg(feature = "derived-geometry")]
pub mod derived_morphism;
#[cfg(feature = "semantic-diff")]
pub mod diff;
#[cfg(feature = "derived-geometry")]
pub mod exit_path;
pub mod fixtures;
pub mod ident;
mod region;
pub mod router;
pub mod sheaf;
#[cfg(feature = "sheaf-merge")]
pub mod sheaf_merge;
#[cfg(feature = "sheaf-repair")]
pub mod sheaf_repair;

pub use convert::{Convert, ConvertDiag, ErrBudget, SampledSdf};
pub use ident::{EntityId, IdTransform, IdentityMap};
pub use region::{
    AgreementConfig, AgreementReport, AgreementScope, AgreementStatus, AgreementUnknown,
    AgreementUnknownReason, Disagreement, Region, RegionChart,
};
pub use sheaf::{
    AdmittedSheafComplex, Interface, InterfaceBound, InterfaceSample, OUTSIDE_RAY_MAX_EVALUATIONS,
    OutsideRaySampleError, OutsideRaySampleReport, RayEndpoint, SHEAF_MAX_CHARTS,
    SHEAF_MAX_INTERFACE_EVALUATIONS, SHEAF_MAX_PAIR_CANDIDATES,
    SHEAF_MAX_RETAINED_INTERFACE_SAMPLES, SHEAF_MAX_TRIPLE_CANDIDATES, SheafAlgebraError,
    SheafBuildError, SheafBuildProgressUnit, SheafComplex, SheafVerdict, TripleCell,
    validate_outside_ray_samples,
};

pub use router::{
    Binding, ChainOutcome, ConverterSpec, CostOracle, CostOracleError, EdgeEvidenceClass,
    EdgeOutcome, EdgeOutcomeError, EdgeRunner, ErrorModel, ExecuteError, ExecuteErrorKind,
    MAX_MEMORY_ORACLE_EDGES, MAX_ROUTER_ID_BYTES, MemoryCostOracle, RouteCandidate,
    RouteExplanation, RoutePlan, RoutePlanError, RouteRefusal, RouteRequest, Router, RouterError,
    ValidatedEdgeObservation,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A point in ℝ³. Minimal geometry-local type (fs-la owns real linear
/// algebra; these exist so charts need no L1 dependency).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
    /// z coordinate.
    pub z: f64,
}

impl Point3 {
    /// Construct a point.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Point3 { x, y, z }
    }

    /// Difference vector `self - other`.
    #[must_use]
    pub fn delta_from(self, other: Point3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    /// Translate by a vector.
    #[must_use]
    pub fn offset(self, v: Vec3) -> Point3 {
        Point3 {
            x: self.x + v.x,
            y: self.y + v.y,
            z: self.z + v.z,
        }
    }
}

/// A vector in ℝ³ (see [`Point3`]'s scope note).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct a vector.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Euclidean norm.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Scale by a scalar.
    #[must_use]
    pub fn scale(self, s: f64) -> Vec3 {
        Vec3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Componentwise minimum corner.
    pub min: Point3,
    /// Componentwise maximum corner.
    pub max: Point3,
}

impl Aabb {
    /// The whole extended Euclidean space. Infinite endpoints are an honest
    /// support declaration for an unbounded region; samplers must resolve such
    /// a support through [`SamplingDomain`] before doing span arithmetic.
    pub const WHOLE_SPACE: Self = Self {
        min: Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        max: Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
    };

    /// Box from corners (normalized componentwise, total function).
    #[must_use]
    pub fn new(a: Point3, b: Point3) -> Self {
        let (min_x, max_x) = normalize_axis(a.x, b.x);
        let (min_y, max_y) = normalize_axis(a.y, b.y);
        let (min_z, max_z) = normalize_axis(a.z, b.z);
        Aabb {
            min: Point3::new(min_x, min_y, min_z),
            max: Point3::new(max_x, max_y, max_z),
        }
    }

    /// True when `p` lies inside or on the boundary.
    #[must_use]
    pub fn contains(&self, p: Point3) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// Smallest box containing both.
    #[must_use]
    pub fn union(&self, other: &Aabb) -> Aabb {
        // Set operations must not turn malformed public AABB fields into a
        // plausible finite support. Preserve the first invalid operand so the
        // shared sampling-domain gate can report its exact axis and bits.
        if !self.is_well_formed() {
            return *self;
        }
        if !other.is_well_formed() {
            return *other;
        }
        Aabb::new(
            Point3::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z),
            ),
            Point3::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z),
            ),
        )
    }

    /// Grow outward by `pad` on every side.
    #[must_use]
    pub fn inflate(&self, pad: f64) -> Aabb {
        Aabb::new(
            Point3::new(self.min.x - pad, self.min.y - pad, self.min.z - pad),
            Point3::new(self.max.x + pad, self.max.y + pad, self.max.z + pad),
        )
    }

    /// Whether this is a valid extended AABB. Correctly oriented infinite
    /// endpoints are allowed; NaNs, inverted axes, `+inf` minima, and `-inf`
    /// maxima are not.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        axis_bounds(*self).into_iter().all(|(_, lo, hi)| {
            !lo.is_nan()
                && !hi.is_nan()
                && lo <= hi
                && lo != f64::INFINITY
                && hi != f64::NEG_INFINITY
        })
    }

    /// Whether every endpoint is finite as well as well formed.
    ///
    /// This deliberately does not promise that `max - min` is finite; use
    /// [`SamplingDomain::resolve`] before midpoint, span, diagonal, or count
    /// arithmetic.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.is_well_formed()
            && axis_bounds(*self)
                .into_iter()
                .all(|(_, lo, hi)| lo.is_finite() && hi.is_finite())
    }

    /// Closed intersection of two well-formed extended boxes. Touching boxes
    /// produce a degenerate intersection; finite-volume samplers subsequently
    /// reject that through [`SamplingDomainError::DegenerateDomain`].
    #[must_use]
    pub fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        if !self.is_well_formed() || !other.is_well_formed() {
            return None;
        }
        let min = Point3::new(
            self.min.x.max(other.min.x),
            self.min.y.max(other.min.y),
            self.min.z.max(other.min.z),
        );
        let max = Point3::new(
            self.max.x.min(other.max.x),
            self.max.y.min(other.max.y),
            self.max.z.min(other.max.z),
        );
        (min.x <= max.x && min.y <= max.y && min.z <= max.z).then_some(Aabb { min, max })
    }
}

fn normalize_axis(a: f64, b: f64) -> (f64, f64) {
    // `f64::min`/`max` deliberately select the numeric operand when the other
    // is NaN. That behavior would launder malformed support before admission.
    if a.is_nan() || b.is_nan() || a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Coordinate axis used by structured sampling-domain refusals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// x axis.
    X,
    /// y axis.
    Y,
    /// z axis.
    Z,
}

impl Axis {
    fn name(self) -> &'static str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
        }
    }
}

fn axis_bounds(box_: Aabb) -> [(Axis, f64, f64); 3] {
    [
        (Axis::X, box_.min.x, box_.max.x),
        (Axis::Y, box_.min.y, box_.max.y),
        (Axis::Z, box_.min.z, box_.max.z),
    ]
}

/// Why an extended chart support could not be admitted as a finite sampling
/// domain. Every variant is detected before evaluation or allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SamplingDomainError {
    /// Chart support contained NaN, inverted bounds, or wrongly oriented
    /// infinity.
    InvalidSupport {
        /// First invalid axis in x/y/z order.
        axis: Axis,
        /// Exact lower-endpoint bits.
        min_bits: u64,
        /// Exact upper-endpoint bits.
        max_bits: u64,
    },
    /// The caller's explicit clip was not a finite, ordered AABB.
    InvalidClip {
        /// First invalid axis in x/y/z order.
        axis: Axis,
        /// Exact lower-endpoint bits.
        min_bits: u64,
        /// Exact upper-endpoint bits.
        max_bits: u64,
    },
    /// No clip resolved an infinite support axis.
    UnboundedSupport {
        /// First unresolved axis in x/y/z order.
        axis: Axis,
    },
    /// Support and explicit clip have no closed intersection.
    EmptyIntersection,
    /// The admitted box has zero width on an axis and cannot drive a 3-D
    /// midpoint/span/count sampler.
    DegenerateDomain {
        /// First zero-width axis in x/y/z order.
        axis: Axis,
    },
    /// Finite endpoints nevertheless overflowed their subtraction.
    NonFiniteSpan {
        /// First overflowing axis in x/y/z order.
        axis: Axis,
        /// Exact lower-endpoint bits.
        min_bits: u64,
        /// Exact upper-endpoint bits.
        max_bits: u64,
    },
    /// The three finite spans have no representable finite Euclidean diagonal.
    NonFiniteDiagonal,
}

impl core::fmt::Display for SamplingDomainError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::InvalidSupport { axis, .. } => {
                write!(
                    f,
                    "sampling refused: invalid chart support on {} axis",
                    axis.name()
                )
            }
            Self::InvalidClip { axis, .. } => write!(
                f,
                "sampling refused: explicit clip must have finite ordered bounds on {} axis",
                axis.name()
            ),
            Self::UnboundedSupport { axis } => write!(
                f,
                "sampling refused: unbounded support on {} axis; provide an explicit finite clip AABB",
                axis.name()
            ),
            Self::EmptyIntersection => {
                write!(
                    f,
                    "sampling refused: explicit clip does not intersect chart support"
                )
            }
            Self::DegenerateDomain { axis } => write!(
                f,
                "sampling refused: admitted domain has zero width on {} axis",
                axis.name()
            ),
            Self::NonFiniteSpan { axis, .. } => write!(
                f,
                "sampling refused: finite {} endpoints overflow their span",
                axis.name()
            ),
            Self::NonFiniteDiagonal => write!(
                f,
                "sampling refused: admitted finite spans have no finite representable diagonal"
            ),
        }
    }
}

impl core::error::Error for SamplingDomainError {}

/// A validated finite, strictly three-dimensional domain. Its private fields
/// ensure midpoint/span/count consumers cannot accidentally accept an extended
/// support without first resolving it against an explicit finite clip.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SamplingDomain {
    bounds: Aabb,
    spans: Vec3,
    diagonal: f64,
}

impl SamplingDomain {
    /// Validate an extended support without requiring it to be finite. This is
    /// useful before union/intersection/inflation, whose arithmetic must not
    /// launder a malformed raw chart support.
    pub fn validate_support(support: Aabb) -> Result<(), SamplingDomainError> {
        validate_support(support)
    }

    /// Resolve `support`, optionally through a finite clip, before any sampling
    /// arithmetic. The effective domain is `support intersection clip` when a
    /// clip is supplied.
    ///
    /// # Errors
    /// [`SamplingDomainError`] identifies the first deterministic admission
    /// failure. No chart evaluation or allocation is performed here.
    pub fn admit(support: Aabb, explicit_clip: Option<Aabb>) -> Result<Self, SamplingDomainError> {
        validate_support(support)?;
        let bounds = if let Some(clip) = explicit_clip {
            validate_clip(clip)?;
            support
                .intersection(&clip)
                .ok_or(SamplingDomainError::EmptyIntersection)?
        } else {
            support
        };

        let mut span_values = [0.0; 3];
        for (index, (axis, lo, hi)) in axis_bounds(bounds).into_iter().enumerate() {
            if !lo.is_finite() || !hi.is_finite() {
                return Err(SamplingDomainError::UnboundedSupport { axis });
            }
            let span = hi - lo;
            if !span.is_finite() {
                return Err(SamplingDomainError::NonFiniteSpan {
                    axis,
                    min_bits: lo.to_bits(),
                    max_bits: hi.to_bits(),
                });
            }
            if span <= 0.0 {
                return Err(SamplingDomainError::DegenerateDomain { axis });
            }
            span_values[index] = span;
        }
        let spans = Vec3::new(span_values[0], span_values[1], span_values[2]);
        let scale = spans.x.max(spans.y).max(spans.z);
        let diagonal = scale
            * ((spans.x / scale) * (spans.x / scale)
                + (spans.y / scale) * (spans.y / scale)
                + (spans.z / scale) * (spans.z / scale))
                .sqrt();
        if !diagonal.is_finite() {
            return Err(SamplingDomainError::NonFiniteDiagonal);
        }
        Ok(Self {
            bounds,
            spans,
            diagonal,
        })
    }

    /// Alias emphasizing that an extended support is being resolved into a
    /// finite domain.
    pub fn resolve(
        support: Aabb,
        explicit_clip: Option<Aabb>,
    ) -> Result<Self, SamplingDomainError> {
        Self::admit(support, explicit_clip)
    }

    /// The admitted finite bounds.
    #[must_use]
    pub const fn bounds(&self) -> Aabb {
        self.bounds
    }

    /// Finite positive per-axis spans.
    #[must_use]
    pub const fn spans(&self) -> Vec3 {
        self.spans
    }

    /// Finite Euclidean diagonal.
    #[must_use]
    pub const fn diagonal(&self) -> f64 {
        self.diagonal
    }

    /// Overflow-safe midpoint of the admitted bounds.
    #[must_use]
    pub fn midpoint(&self) -> Point3 {
        Point3::new(
            self.bounds.min.x + 0.5 * self.spans.x,
            self.bounds.min.y + 0.5 * self.spans.y,
            self.bounds.min.z + 0.5 * self.spans.z,
        )
    }

    /// Largest finite axis span.
    #[must_use]
    pub fn max_span(&self) -> f64 {
        self.spans.x.max(self.spans.y).max(self.spans.z)
    }
}

fn validate_support(support: Aabb) -> Result<(), SamplingDomainError> {
    for (axis, lo, hi) in axis_bounds(support) {
        if lo.is_nan() || hi.is_nan() || lo > hi || lo == f64::INFINITY || hi == f64::NEG_INFINITY {
            return Err(SamplingDomainError::InvalidSupport {
                axis,
                min_bits: lo.to_bits(),
                max_bits: hi.to_bits(),
            });
        }
    }
    Ok(())
}

fn validate_clip(clip: Aabb) -> Result<(), SamplingDomainError> {
    for (axis, lo, hi) in axis_bounds(clip) {
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(SamplingDomainError::InvalidClip {
                axis,
                min_bits: lo.to_bits(),
                max_bits: hi.to_bits(),
            });
        }
    }
    Ok(())
}

/// Betti-number bounds `(lower, upper)` per dimension — the topology hint
/// charts may advertise without proving (certificates are a later bead;
/// `unknown()` is the honest default).
///
/// The numbers describe the SOLID REGION the chart presents, i.e. the set
/// where the signed distance is negative — not its boundary surface. This
/// is the convention `fs_topo::cubical::verify_topology` measures against
/// (a ball is `(1, 0, 0)`, a solid torus `(1, 1, 0)`, a hollow ball
/// `(1, 0, 1)`), so a hint and a measurement of the same chart are
/// directly comparable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BettiBounds {
    /// Connected components of the solid region.
    pub b0: (u32, u32),
    /// Tunnels/handles of the solid region.
    pub b1: (u32, u32),
    /// Enclosed voids (cavities) INSIDE the solid region. A solid ball,
    /// box or torus encloses none; only a hollow shell does.
    pub b2: (u32, u32),
}

impl BettiBounds {
    /// No topology claim at all.
    #[must_use]
    pub const fn unknown() -> Self {
        BettiBounds {
            b0: (0, u32::MAX),
            b1: (0, u32::MAX),
            b2: (0, u32::MAX),
        }
    }

    /// Exact known Betti numbers.
    #[must_use]
    pub const fn exact(b0: u32, b1: u32, b2: u32) -> Self {
        BettiBounds {
            b0: (b0, b0),
            b1: (b1, b1),
            b2: (b2, b2),
        }
    }
}

/// Differentiability class a chart advertises for its signed distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Differentiability {
    /// No continuity claim. Sampling, finite-budget search, or discrete
    /// selection may change the returned value discontinuously.
    Unknown,
    /// Continuous only (gradients may be `None` anywhere).
    C0,
    /// Continuously differentiable away from the medial axis.
    C1,
    /// Smooth away from the medial axis.
    Smooth,
}

/// The theorem a chart exposes to ray steppers. A finite Lipschitz number by
/// itself is not enough: the field must also have a certified relationship to
/// the represented boundary. The default is deliberately no-claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceStepClaim {
    /// No generic no-tunneling theorem. Callers may offer an explicitly
    /// uncertified preview, but production render paths must fail closed.
    NoClaim,
    /// The represented real field is the exact signed distance. Each sample's
    /// `error` is either a genuinely exact singleton or a rigorous enclosure
    /// of its rounded evaluation; steppers use the enclosure endpoint closest
    /// to zero as the no-tunneling radius.
    ExactDistance,
    /// The field has the exact sign and zero set of the represented region;
    /// the represented real field is continuous on every finite line segment;
    /// each sample's positive finite Lipschitz bound is certified over the
    /// entire closed `|f| / L` step ball. [`Chart::trace_value_enclosure`]
    /// encloses the real implicit-field evaluation used for that step, making
    /// the radius safe even when the magnitude is not the exact distance. The
    /// continuity theorem also lets a consumer turn rigorously opposite signs
    /// at the ends of a short segment into existence of a zero inside it; the
    /// Lipschitz bound alone still supplies no upper proximity bound.
    LipschitzImplicit,
}

/// One signed-field query's answer (plan Appendix B: value + gradient +
/// certified Lipschitz data + the declared abstract-distance error model).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChartSample {
    /// Signed scalar representative of the region (negative inside). It is the
    /// Euclidean signed distance only under [`TraceStepClaim::ExactDistance`];
    /// a Lipschitz-implicit chart preserves sign and the zero set instead.
    pub signed_distance: f64,
    /// Gradient of the reported scalar field where it exists (`None` on medial
    /// axes/edges or for C0 charts).
    pub gradient: Option<Vec3>,
    /// Certified LOCAL Lipschitz bound for the reported scalar field near the
    /// query (sphere-tracing fuel; `None` = no claim). A chart opting into
    /// [`TraceStepClaim::LipschitzImplicit`] strengthens this to validity over
    /// the entire closed step ball specified by that claim.
    pub lipschitz: Option<f64>,
    /// Declared error of `signed_distance` relative to the ABSTRACT region
    /// (fs-evidence certificate: by-construction identities may say Exact,
    /// rounded analytic evaluations say Enclosure, and heuristics say Estimate).
    pub error: NumericalCertificate,
}

/// The chart contract: a concrete presentation of an abstract region.
/// Object-safe (see module docs for the `Param` note).
pub trait Chart: Send + Sync {
    /// Answer a signed-distance query. Implementations poll
    /// `cx.checkpoint()` at bounded strides inside expensive evaluations.
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample;

    /// A box guaranteed to contain the region (queries outside are
    /// positive-distance by definition).
    fn support(&self) -> Aabb;

    /// State the certified relationship that makes this chart safe for a ray
    /// stepper. Implementations must opt in explicitly; a `Some(lipschitz)`
    /// sample does not upgrade the default no-claim.
    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::NoClaim
    }

    /// Rigorous enclosure of the real scalar field used by the typed trace
    /// theorem at `x`. This is distinct from [`ChartSample::error`], which is
    /// relative to the abstract region's signed distance and therefore may be
    /// only an `Estimate` for a non-distance implicit field.
    ///
    /// Exact-distance charts inherit their sample's distance certificate.
    /// Lipschitz-implicit charts must override this method; the default refuses
    /// to promote a rounded field value into a certified trace step.
    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if self.trace_step_claim() == TraceStepClaim::ExactDistance {
            sample.error
        } else {
            NumericalCertificate::no_claim()
        }
    }

    /// Topology bounds this chart is willing to state.
    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    /// Stable chart-kind name (provenance, reports, router tables).
    fn name(&self) -> &'static str;

    /// Advertised differentiability class of the signed distance.
    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }

    /// Convenience: strict inside test via the SDF convention.
    fn inside(&self, x: Point3, cx: &Cx<'_>) -> bool {
        self.eval(x, cx).signed_distance < 0.0
    }
}

/// A finite geometric intersection between an arbitrary source chart and an
/// explicit clip AABB. Unlike merely replacing `support()`, this wrapper also
/// intersects the represented negative set: its field is
/// `max(source_field, finite_box_sdf)`, so it is negative exactly where both
/// the source and clip are inside.
///
/// The hard maximum and box edges make the wrapper C0. Its sign and support
/// are honest, but v1 deliberately makes no abstract-distance or certified ray
/// step claim for the composite magnitude.
pub struct ClippedChart<'a> {
    source: &'a dyn Chart,
    clip: Aabb,
    domain: SamplingDomain,
}

impl<'a> ClippedChart<'a> {
    /// Intersect `source` with a finite clip before any evaluation.
    ///
    /// # Errors
    /// [`SamplingDomainError`] when source support or clip cannot produce a
    /// finite, positive-volume domain.
    pub fn new(source: &'a dyn Chart, clip: Aabb) -> Result<Self, SamplingDomainError> {
        let domain = SamplingDomain::resolve(source.support(), Some(clip))?;
        Ok(Self {
            source,
            clip,
            domain,
        })
    }

    /// The admitted finite intersection domain.
    #[must_use]
    pub fn domain(&self) -> SamplingDomain {
        self.domain
    }

    /// The caller-supplied finite box whose exact SDF participates in the
    /// composite field.
    #[must_use]
    pub fn clip(&self) -> Aabb {
        self.clip
    }

    /// The wrapped source chart.
    #[must_use]
    pub fn source(&self) -> &'a dyn Chart {
        self.source
    }
}

impl Chart for ClippedChart<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let source = self.source.eval(x, cx);
        let clip_distance = signed_distance_to_box(x, self.clip);
        let (signed_distance, gradient) = if source.signed_distance > clip_distance {
            (source.signed_distance, source.gradient)
        } else if clip_distance > source.signed_distance {
            (clip_distance, None)
        } else {
            (source.signed_distance, None)
        };
        ChartSample {
            signed_distance,
            gradient,
            lipschitz: source
                .lipschitz
                .filter(|bound| bound.is_finite() && *bound >= 0.0)
                .map(|bound| bound.max(1.0)),
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        self.domain.bounds()
    }

    fn name(&self) -> &'static str {
        "geom/clipped"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

fn signed_distance_to_box(p: Point3, box_: Aabb) -> f64 {
    let outside = Vec3::new(
        if p.x < box_.min.x {
            box_.min.x - p.x
        } else if p.x > box_.max.x {
            p.x - box_.max.x
        } else {
            0.0
        },
        if p.y < box_.min.y {
            box_.min.y - p.y
        } else if p.y > box_.max.y {
            p.y - box_.max.y
        } else {
            0.0
        },
        if p.z < box_.min.z {
            box_.min.z - p.z
        } else if p.z > box_.max.z {
            p.z - box_.max.z
        } else {
            0.0
        },
    );
    if outside.x > 0.0 || outside.y > 0.0 || outside.z > 0.0 {
        let scale = outside.x.max(outside.y).max(outside.z);
        return scale
            * ((outside.x / scale) * (outside.x / scale)
                + (outside.y / scale) * (outside.y / scale)
                + (outside.z / scale) * (outside.z / scale))
                .sqrt();
    }
    let face_distance = (p.x - box_.min.x)
        .min(box_.max.x - p.x)
        .min(p.y - box_.min.y)
        .min(box_.max.y - p.y)
        .min(p.z - box_.min.z)
        .min(box_.max.z - p.z);
    -face_distance
}

/// A chart with design levers: the differentiable map θ → Region handle
/// (plan §7.6; fs-xform builds the parameterization zoo on this).
pub trait DesignChart: Chart {
    /// The design-lever handle.
    type Param;

    /// Report the current lever value.
    fn param(&self) -> &Self::Param;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn aabb_laws() {
        let a = Aabb::new(Point3::new(1.0, 2.0, 3.0), Point3::new(-1.0, 0.0, 5.0));
        assert_eq!(a.min, Point3::new(-1.0, 0.0, 3.0), "corners normalize");
        assert!(a.contains(Point3::new(0.0, 1.0, 4.0)));
        assert!(!a.contains(Point3::new(2.0, 1.0, 4.0)));
        let b = Aabb::new(Point3::new(5.0, 5.0, 5.0), Point3::new(6.0, 6.0, 6.0));
        let u = a.union(&b);
        assert!(u.contains(Point3::new(5.5, 5.5, 5.5)) && u.contains(Point3::new(0.0, 1.0, 4.0)));
        assert!(a.inflate(1.0).contains(Point3::new(1.5, 0.5, 4.0)));
    }

    #[test]
    fn sampling_domain_refuses_unresolved_and_malformed_support() {
        let malformed = Aabb::new(
            Point3::new(f64::NAN, -1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
        );
        assert!(malformed.min.x.is_nan(), "Aabb::new must not launder NaN");
        assert!(matches!(
            SamplingDomain::resolve(malformed, None),
            Err(SamplingDomainError::InvalidSupport { axis: Axis::X, .. })
        ));
        assert!(matches!(
            SamplingDomain::resolve(Aabb::WHOLE_SPACE, None),
            Err(SamplingDomainError::UnboundedSupport { axis: Axis::X })
        ));
        let huge = Aabb::new(
            Point3::new(-f64::MAX, -1.0, -1.0),
            Point3::new(f64::MAX, 1.0, 1.0),
        );
        assert!(matches!(
            SamplingDomain::resolve(huge, None),
            Err(SamplingDomainError::NonFiniteSpan { axis: Axis::X, .. })
        ));
    }

    #[test]
    fn finite_clip_resolves_whole_space_without_unstable_arithmetic() {
        let clip = Aabb::new(Point3::new(-4.0, -3.0, -2.0), Point3::new(6.0, 5.0, 4.0));
        let domain = SamplingDomain::resolve(Aabb::WHOLE_SPACE, Some(clip))
            .expect("finite clip admits whole-space support");
        assert_eq!(domain.bounds(), clip);
        assert_eq!(domain.spans(), Vec3::new(10.0, 8.0, 6.0));
        assert_eq!(domain.midpoint(), Point3::new(1.0, 1.0, 1.0));
        assert!(domain.diagonal().is_finite());
        let bounded = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
        assert_eq!(Aabb::WHOLE_SPACE.intersection(&bounded), Some(bounded));
    }

    #[test]
    fn aabb_union_preserves_malformed_support_for_admission() {
        let malformed = Aabb {
            min: Point3::new(f64::NAN, -1.0, -1.0),
            max: Point3::new(1.0, 1.0, 1.0),
        };
        let bounded = Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0));
        for union in [malformed.union(&bounded), bounded.union(&malformed)] {
            assert!(union.min.x.is_nan());
            assert!(matches!(
                SamplingDomain::resolve(union, None),
                Err(SamplingDomainError::InvalidSupport { axis: Axis::X, .. })
            ));
        }
    }

    #[test]
    fn vec_ops_are_the_usual_ones() {
        let v = Point3::new(1.0, 2.0, 2.0).delta_from(Point3::new(0.0, 0.0, 0.0));
        assert!((v.norm() - 3.0).abs() < 1e-12);
        assert!((v.dot(Vec3::new(1.0, 0.0, 0.0)) - 1.0).abs() < 1e-12);
        assert_eq!(v.scale(2.0), Vec3::new(2.0, 4.0, 4.0));
    }

    #[test]
    fn betti_bounds_default_to_no_claim() {
        let u = BettiBounds::unknown();
        assert_eq!(u.b0, (0, u32::MAX));
        assert_eq!(BettiBounds::exact(1, 0, 1).b2, (1, 1));
    }
}

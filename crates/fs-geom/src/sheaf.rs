//! SHEAF CERTIFICATES (plan §7.3, Bet 11): a finite, constant-scalar
//! graph-gauge substrate inspired by cellular sheaves, with sampled
//! interface-agreement evidence for multi-representation models. The current
//! implementation retains patch-adjacency incidence and independently sampled
//! interface rows; it does not yet implement general stalk spaces or admitted
//! trace/conversion restriction maps and therefore is not a full cellular-sheaf
//! model. The current positive certificate is an INTERVAL-VERIFIED bound
//! `‖δs‖∞ ≤ tol` on sampled interface mismatches. When a mismatch
//! enclosure lies entirely above tolerance, the offending interface cells are
//! reported as proven interface violations. This base certificate does not
//! establish between-sample coverage, continuum watertightness, cocycle
//! membership, or non-exactness and therefore makes no global or H¹ claim.
//!
//! The construction is finite linear algebra: the edge-level
//! [`SheafComplex::delta0_edges`] and [`SheafComplex::delta1`] maps validate and
//! cap requests fallibly before the current infallible `Coo` staging/assembly
//! produces sparse matrices with entries in {−1, 0, +1}. Thus their
//! `δ¹·δ⁰ = 0` identity holds BITWISE — small-integer f64 arithmetic is exact. The separate
//! sample-row restriction incidence is [`SheafComplex::delta0`]. The
//! least-squares section solve (per-patch gauge offsets over the adjacency
//! Laplacian) reports the fractional reduction in uncentered sample-level
//! midpoint-mismatch mean-square energy. That graph-gauge diagnostic is not a cohomology
//! certificate; the feature-gated repair classifier owns exact/coexact/harmonic
//! claims.

use crate::{Aabb, Chart, ChartSample, Point3, SamplingDomain, SamplingDomainError};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, NumericalKind, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate,
};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_sparse::{Coo, Csr};
use std::fmt::Write as _;

/// Samples drawn per pairwise interface.
pub const SAMPLES_PER_INTERFACE: usize = 32;

/// Zero-band half-width as a fraction of the overlap-box diagonal:
/// a point belongs to the shared surface region when BOTH charts place
/// it within this band of their zero set.
pub const BAND_FRACTION: f64 = 0.05;

/// Maximum number of chart evaluations admitted by one outside-ray sample run.
///
/// The sign-sequence replay diagnostic is deliberately bounded rather than
/// hiding an unbounded marching workload behind a sample-count argument.
pub const OUTSIDE_RAY_MAX_EVALUATIONS: usize = 1_048_576;

/// Maximum number of neighbor-membership probes admitted while discovering
/// fully connected triples during one sheaf build.
pub const SHEAF_MAX_TRIPLE_CANDIDATES: usize = 1_048_576;

/// Maximum charts whose supports may enter one sampled-interface build.
pub const SHEAF_MAX_CHARTS: usize = 4_096;

/// Maximum chart-pair support-overlap probes admitted before sampling.
pub const SHEAF_MAX_PAIR_CANDIDATES: usize = 1_048_576;

/// Maximum worst-case chart evaluations admitted across all overlapping pairs.
pub const SHEAF_MAX_INTERFACE_EVALUATIONS: usize = 16_777_216;

/// Maximum retained interface samples a build may allocate.
pub const SHEAF_MAX_RETAINED_INTERFACE_SAMPLES: usize = 131_072;

/// Allocation-free writer for the legacy FNV provenance stream.
///
/// `watertightness` historically materialized its complete canonical transcript
/// in one `String` before hashing it.  A public complex can contain many samples,
/// so that duplicated all evidence bytes without an admission budget.  Streaming
/// the identical bytes preserves the legacy fingerprint while keeping auxiliary
/// memory constant.  This remains a non-cryptographic fingerprint; strong
/// identity migration is tracked separately.
struct LegacyProvenanceWriter(u64);

impl LegacyProvenanceWriter {
    const fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    const fn finish(self) -> u64 {
        self.0
    }
}

impl std::fmt::Write for LegacyProvenanceWriter {
    fn write_str(&mut self, value: &str) -> std::fmt::Result {
        for byte in value.bytes() {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
        Ok(())
    }
}

/// Endpoint of a ray named by a structured sampling refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RayEndpoint {
    /// Segment start.
    Start,
    /// Segment end.
    End,
}

/// Why outside-to-outside ray sample validation could not return a report.
#[derive(Debug, Clone, PartialEq)]
pub enum OutsideRaySampleError {
    /// A union-of-charts model needs at least one presentation.
    EmptyCharts,
    /// A diagnostic run with no rays gathers no evidence.
    EmptyRays,
    /// At least one interval is required to define a segment march.
    InvalidSteps {
        /// Offending interval count.
        steps: usize,
    },
    /// The requested diagnostic exceeds its public deterministic work cap.
    WorkLimitExceeded {
        /// Requested chart evaluations.
        requested: u128,
        /// Public chart-evaluation cap.
        cap: usize,
    },
    /// A ray endpoint was NaN or infinite.
    NonFiniteEndpoint {
        /// Ray index.
        ray: usize,
        /// Which endpoint was invalid.
        endpoint: RayEndpoint,
        /// Offending point.
        point: Point3,
    },
    /// Finite endpoints did not yield a representable convex interpolation.
    NonRepresentableSamplePoint {
        /// Ray index.
        ray: usize,
        /// Sample index in `0..=steps`.
        step: usize,
        /// Offending interpolated point.
        point: Point3,
    },
    /// A chart returned a NaN or infinite nominal field value.
    NonFiniteSample {
        /// Ray index.
        ray: usize,
        /// Sample index in `0..=steps`.
        step: usize,
        /// Chart index.
        chart: usize,
        /// Offending value.
        value: f64,
    },
    /// The validator requires both endpoints to be strictly outside the union
    /// model so transition telemetry has an unambiguous replay contract.
    EndpointNotOutside {
        /// Ray index.
        ray: usize,
        /// Which endpoint violated the precondition.
        endpoint: RayEndpoint,
        /// Minimum signed field across the charts.
        min_signed_distance: f64,
    },
    /// An endpoint is nominally outside, but at least one chart neither excludes
    /// it through its support nor supplies a rigorous positive distance
    /// enclosure. Nominal sign alone cannot establish the endpoint precondition.
    EndpointOutsideUnproven {
        /// Ray index.
        ray: usize,
        /// Which endpoint lacked proof.
        endpoint: RayEndpoint,
        /// Chart that lacked outside authority.
        chart: usize,
        /// Nominal value retained for diagnosis.
        nominal: f64,
        /// Declared numerical certificate retained for diagnosis.
        certificate: NumericalCertificate,
    },
    /// Cancellation was observed before a report could be published.
    Cancelled {
        /// Rays fully classified before cancellation.
        completed_rays: usize,
        /// Ray points fully evaluated before cancellation.
        completed_points: usize,
        /// Individual chart evaluations completed before cancellation.
        completed_chart_evaluations: usize,
    },
}

impl core::fmt::Display for OutsideRaySampleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyCharts => write!(f, "outside-ray sampling refused: the chart set is empty"),
            Self::EmptyRays => write!(f, "outside-ray sampling refused: the ray set is empty"),
            Self::InvalidSteps { steps } => {
                write!(
                    f,
                    "outside-ray sampling refused: steps must be positive, got {steps}"
                )
            }
            Self::WorkLimitExceeded { requested, cap } => write!(
                f,
                "outside-ray sampling refused: {requested} chart evaluations exceed the public cap {cap}"
            ),
            Self::NonFiniteEndpoint {
                ray,
                endpoint,
                point,
            } => write!(
                f,
                "outside-ray sampling refused: ray {ray} {endpoint:?} endpoint is non-finite: {point:?}"
            ),
            Self::NonRepresentableSamplePoint { ray, step, point } => write!(
                f,
                "outside-ray sampling refused: ray {ray} sample {step} is not representable: {point:?}"
            ),
            Self::NonFiniteSample {
                ray,
                step,
                chart,
                value,
            } => write!(
                f,
                "outside-ray sampling refused: chart {chart} returned non-finite value {value} at ray {ray} sample {step}"
            ),
            Self::EndpointNotOutside {
                ray,
                endpoint,
                min_signed_distance,
            } => write!(
                f,
                "outside-ray sampling refused: ray {ray} {endpoint:?} endpoint is not nominally outside (minimum field {min_signed_distance})"
            ),
            Self::EndpointOutsideUnproven {
                ray,
                endpoint,
                chart,
                nominal,
                certificate,
            } => write!(
                f,
                "outside-ray sampling refused: ray {ray} {endpoint:?} endpoint is nominally outside at chart {chart} ({nominal}) but lacks a rigorous positive enclosure or excluding support ({certificate:?})"
            ),
            Self::Cancelled {
                completed_rays,
                completed_points,
                completed_chart_evaluations,
            } => write!(
                f,
                "outside-ray sampling cancelled after {completed_rays} rays, {completed_points} points, and {completed_chart_evaluations} chart evaluations"
            ),
        }
    }
}

impl core::error::Error for OutsideRaySampleError {}

/// Work and nominal transition counts retained by a successful proven-outside
/// to proven-outside ray sample validation. Only endpoint outside status is
/// certificate-checked; interior toggles come from finite nominal minimum-field
/// signs and may include samples whose certificates are `NoClaim`. `toggles`
/// is necessarily even on every admitted ray; it is replay telemetry, not
/// topology evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OutsideRaySampleReport {
    /// Rays fully validated.
    pub rays: usize,
    /// Ray sample points fully evaluated.
    pub sample_points: usize,
    /// Individual chart evaluations performed.
    pub chart_evaluations: usize,
    /// Total nominal minimum-field sign transitions across all rays.
    pub toggles: usize,
}

/// One shared interface sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceSample {
    /// The world-space point.
    pub point: Point3,
    /// Stored intervals for the two charts' signed distances. In an
    /// [`AdmittedSheafComplex`] these were derived from the charts' rigorous
    /// error certificates; values in a raw public [`SheafComplex`] are
    /// caller-supplied diagnostics and carry no such authority.
    pub values: [Interval; 2],
}

/// One interface (edge of the adjacency complex).
#[derive(Debug, Clone)]
pub struct Interface {
    /// Patch indices (u < v; the edge is oriented u → v).
    pub patches: (usize, usize),
    /// Shared samples.
    pub samples: Vec<InterfaceSample>,
}

/// An unverified pairwise-interface clique completion (candidate 2-cell).
/// Pairwise sampled overlaps do not by themselves prove a common triple
/// overlap or aligned restriction samples, so this carries no Čech/topology
/// authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TripleCell {
    /// Patch indices (sorted).
    pub patches: (usize, usize, usize),
    /// Minimum of the three independent pairwise sample counts. This is not a
    /// count of common aligned triple samples.
    pub samples: usize,
}

/// The patch-adjacency sheaf complex.
#[derive(Debug)]
pub struct SheafComplex {
    /// Explicit finite scope used to gather this evidence. `None` means the
    /// admitted global supports were sampled.
    pub sampling_clip: Option<Aabb>,
    /// Patch count.
    pub n_patches: usize,
    /// Pairwise interfaces (sorted by patch pair — deterministic).
    pub interfaces: Vec<Interface>,
    /// Unverified pairwise-interface clique completions.
    pub triples: Vec<TripleCell>,
}

/// Immutable complex produced by the chart-sampling admission path. Only this
/// wrapper may publish positive or negative sampled-interface evidence. It
/// dereferences immutably for incidence algebra and diagnostics, but exposes no
/// `DerefMut` or ownership escape that could mutate retained evidence after
/// admission.
///
/// ```compile_fail
/// fn cannot_mutate(
///     admitted: &mut fs_geom::AdmittedSheafComplex,
///     fabricated: fs_geom::Interface,
/// ) {
///     admitted.interfaces.push(fabricated);
/// }
/// ```
#[derive(Debug)]
pub struct AdmittedSheafComplex {
    inner: SheafComplex,
}

impl core::ops::Deref for AdmittedSheafComplex {
    type Target = SheafComplex;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl AdmittedSheafComplex {
    /// Per-interface sampled mismatch bounds from immutable builder-retained
    /// chart evidence. Raw public complexes intentionally have no originless
    /// method returning this positive/negative evidence type.
    #[must_use]
    pub fn mismatch_bounds(&self) -> Result<Vec<InterfaceBound>, SheafAlgebraError> {
        self.inner.mismatch_bounds()
    }

    /// Assess builder-retained sampled interfaces. This authority says only
    /// that the declared chart certificates were sampled through the admitted
    /// path and then kept immutable; it does not independently prove a chart
    /// implementation truthful.
    #[must_use]
    pub fn watertightness(&self, tol: f64) -> Evidence<SheafVerdict> {
        self.inner.assess_sampled_agreement(tol, true)
    }
}

/// One admitted interface's context-free numeric mismatch enclosure.
/// Construction is private so raw caller-fabricated complexes cannot present
/// originless chart samples as builder-retained diagnostics. Tolerance
/// predicates deliberately take the tolerance at use time: no detached boolean
/// can be reused across tolerances.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceBound {
    /// The patch pair.
    patches: (usize, usize),
    /// The interface has ordered in-range patch indices and at least one sample,
    /// and every sample point and mismatch-interval endpoint is finite.
    determinate: bool,
    /// At least one sample on an otherwise valid interface has a finite
    /// mismatch enclosure. This is enough for an existential localized leak,
    /// even when another sample is indeterminate; it is not enough for PASS.
    has_determinate_sample: bool,
    /// Reported worst lower bound.
    lo_report: f64,
    /// Reported worst upper bound.
    hi_report: f64,
}

/// Structured refusal for algebra over public raw sheaf-complex parts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafAlgebraError {
    /// Patch/interface/triple ordering, indices, or sample data are malformed.
    InvalidStructure,
    /// A raw diagnostic request exceeds a deterministic legacy-path ceiling.
    WorkLimit {
        /// Stable operation stage.
        stage: &'static str,
        /// Requested items.
        requested: u128,
        /// Public ceiling.
        cap: usize,
    },
    /// A bounded diagnostic allocation could not be reserved.
    ResourceExhausted {
        /// Stable allocation stage.
        stage: &'static str,
    },
    /// A section diagnostic requires finite interval midpoints, but a retained
    /// sample was unbounded or otherwise indeterminate.
    IndeterminateSampleValue {
        /// Interface index in canonical edge order.
        interface: usize,
        /// Sample index within that interface.
        sample: usize,
    },
    /// Finite inputs overflowed a section-diagnostic intermediate.
    NumericalOverflow {
        /// Stable arithmetic stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for SheafAlgebraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidStructure => write!(f, "raw sheaf algebra requires a valid complex"),
            Self::WorkLimit {
                stage,
                requested,
                cap,
            } => write!(
                f,
                "raw sheaf algebra stage {stage} requests {requested} items above cap {cap}"
            ),
            Self::ResourceExhausted { stage } => {
                write!(f, "raw sheaf algebra could not reserve storage for {stage}")
            }
            Self::IndeterminateSampleValue { interface, sample } => write!(
                f,
                "section diagnostic requires a finite midpoint at interface {interface}, sample {sample}"
            ),
            Self::NumericalOverflow { stage } => {
                write!(f, "section diagnostic arithmetic overflowed during {stage}")
            }
        }
    }
}

impl core::error::Error for SheafAlgebraError {}

impl InterfaceBound {
    /// Patch pair.
    #[must_use]
    pub const fn patches(self) -> (usize, usize) {
        self.patches
    }

    /// Whether every retained mismatch enclosure lies within `tolerance`.
    /// Invalid tolerances and indeterminate bounds fail closed.
    #[must_use]
    pub fn all_within(self, tolerance: f64) -> bool {
        self.determinate && tolerance.is_finite() && tolerance >= 0.0 && self.hi_report <= tolerance
    }

    /// Whether at least one retained mismatch enclosure lies wholly above
    /// `tolerance`. Invalid tolerances fail closed.
    #[must_use]
    pub fn proven_leak(self, tolerance: f64) -> bool {
        self.has_determinate_sample
            && tolerance.is_finite()
            && tolerance >= 0.0
            && self.lo_report > tolerance
    }

    /// Whether structure, tolerance, points, and interval endpoints were
    /// determinate for this interface.
    #[must_use]
    pub const fn determinate(self) -> bool {
        self.determinate
    }

    /// Reported worst lower bound.
    #[must_use]
    pub const fn lo_report(self) -> f64 {
        self.lo_report
    }

    /// Reported worst upper bound.
    #[must_use]
    pub const fn hi_report(self) -> f64 {
        self.hi_report
    }
}

/// The certificate verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum SheafVerdict {
    /// Every retained sample proves `‖δs‖∞ ≤ tol`, with per-interface
    /// margins attached. This is sampled agreement, not by itself a continuum
    /// covering or global watertightness theorem.
    Pass {
        /// Upper bound of the worst interface mismatch.
        worst_mismatch: f64,
        /// Per-interface (patch pair, mismatch upper bound).
        margins: Vec<((usize, usize), f64)>,
    },
    /// Some interface's mismatch enclosure lies entirely above tolerance —
    /// an interval-proven interface violation, localized without a topology
    /// claim.
    Fail {
        /// Offending interfaces: (patch pair, mismatch lower bound).
        interface_violations: Vec<((usize, usize), f64)>,
        /// Fractional reduction in uncentered midpoint-mismatch mean-square
        /// energy from per-patch graph gauge offsets, in `[0, 1]`. Near 1
        /// means a constant re-gauge fits the sampled edge means; it does not
        /// prove exactness or classify the residual topologically. `None`
        /// refuses the diagnostic when its least-squares arithmetic is not
        /// representable.
        gauge_fit_share: Option<f64>,
    },
    /// No sound aggregate claim: enclosures may straddle the tolerance, or the
    /// retained structure/scope/tolerance may be absent or malformed.
    Unknown {
        /// Reported non-authoritative interface bounds: (patch pair, lower,
        /// upper). For admitted evidence these are indeterminate/straddling
        /// bounds; for raw public parts every available diagnostic bound is
        /// returned because raw origin itself forces `Unknown`. This may be
        /// empty when the structure or interface set is absent.
        reported_bounds: Vec<((usize, usize), f64, f64)>,
    },
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// Geometry-derived deterministic seed: identical boxes give identical
/// samples regardless of patch indexing (re-index invariance is exact).
fn box_seed(b: &Aabb) -> u64 {
    let mut bytes = Vec::with_capacity(48);
    for v in [b.min.x, b.min.y, b.min.z, b.max.x, b.max.y, b.max.z] {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    fnv(&bytes)
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

/// Unit carried by cooperative build-progress diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafBuildProgressUnit {
    /// Admission checkpoints.
    AdmissionChecks,
    /// Chart supports snapshotted.
    Charts,
    /// Candidate chart pairs inspected.
    PairCandidates,
    /// Rejection-sampling draws completed for one chart pair.
    InterfaceDraws,
    /// Builder edges scanned while constructing triple adjacency.
    Edges,
    /// Neighbor-membership probes completed during triple discovery.
    NeighborProbes,
    /// Triple cells retained before final publication.
    TripleCells,
}

/// Why sheaf interface discovery could not safely sample a chart pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafBuildError {
    /// Cooperative cancellation was observed before publishing a complex.
    Cancelled {
        /// Stable build stage at which cancellation was observed.
        stage: &'static str,
        /// Pair being sampled, when cancellation occurred inside an interface.
        patches: Option<(usize, usize)>,
        /// Number of fully completed units.
        completed_work: usize,
        /// Exact unit measured by `completed_work`.
        unit: SheafBuildProgressUnit,
    },
    /// The requested finite sampling clip itself was not admissible. This is
    /// checked even when there are fewer than two charts or every chart pair
    /// is disjoint, so malformed caller input cannot silently succeed.
    SamplingClip {
        /// Exact clip admission failure.
        source: SamplingDomainError,
    },
    /// A chart's support was malformed before pair discovery began.
    ChartSupport {
        /// Chart index.
        chart: usize,
        /// Exact support validation failure.
        source: SamplingDomainError,
    },
    /// Pair-attributed finite-domain admission refusal.
    SamplingDomain {
        /// Chart indices in deterministic ascending order.
        patches: (usize, usize),
        /// Exact shared admission failure.
        source: SamplingDomainError,
    },
    /// A chart returned a non-finite signed-distance sample while an
    /// interface was being discovered. Such a producer cannot be treated as
    /// merely outside the sampled zero band.
    NonFiniteSample {
        /// Chart pair being sampled.
        patches: (usize, usize),
        /// Index of the chart that returned the malformed value.
        chart: usize,
        /// Exact sampled point coordinates, encoded as IEEE-754 bits.
        point_bits: [u64; 3],
        /// Exact malformed signed-distance bits.
        value_bits: u64,
        /// Pair draws fully evaluated before the refusal.
        completed_draws: usize,
    },
    /// A preflighted build stage exceeds its deterministic work/memory ceiling.
    BuildWorkLimit {
        /// Stable preflight stage.
        stage: &'static str,
        /// Requested work units or retained items.
        requested: u128,
        /// Public ceiling.
        cap: usize,
    },
    /// A bounded build allocation could not be reserved.
    ResourceExhausted {
        /// Stable allocation stage.
        stage: &'static str,
    },
    /// A builder-owned invariant was violated before publication.
    InternalInvariant {
        /// Stable invariant stage.
        stage: &'static str,
    },
    /// Triple discovery exceeded its deterministic membership-probe cap.
    TripleWorkLimit {
        /// Neighbor membership probes completed before refusal.
        completed_probes: usize,
        /// Public work cap.
        cap: usize,
    },
}

impl core::fmt::Display for SheafBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                stage,
                patches,
                completed_work,
                unit,
            } => write!(
                f,
                "sheaf build cancelled during {stage} for patches {patches:?} after \
                 {completed_work} completed {unit:?}; no complex was published"
            ),
            Self::SamplingClip { source } => {
                write!(f, "sheaf explicit sampling clip {source}")
            }
            Self::ChartSupport { chart, source } => {
                write!(f, "sheaf chart {chart} support {source}")
            }
            Self::SamplingDomain { patches, source } => {
                write!(f, "sheaf interface ({}, {}) {source}", patches.0, patches.1)
            }
            Self::NonFiniteSample {
                patches,
                chart,
                point_bits,
                value_bits,
                completed_draws,
            } => write!(
                f,
                "sheaf interface ({}, {}) chart {chart} returned non-finite signed-distance \
                 bits {value_bits:#018x} at point bits [{:#018x}, {:#018x}, {:#018x}] after \
                 {completed_draws} completed draws",
                patches.0, patches.1, point_bits[0], point_bits[1], point_bits[2]
            ),
            Self::BuildWorkLimit {
                stage,
                requested,
                cap,
            } => write!(
                f,
                "sheaf build stage {stage} requests {requested} work units/items; the deterministic cap is {cap}"
            ),
            Self::ResourceExhausted { stage } => {
                write!(
                    f,
                    "sheaf build could not reserve bounded storage for {stage}"
                )
            }
            Self::InternalInvariant { stage } => {
                write!(f, "sheaf builder invariant failed during {stage}")
            }
            Self::TripleWorkLimit {
                completed_probes,
                cap,
            } => write!(
                f,
                "sheaf triple discovery refused after {completed_probes} completed neighbor probes; the deterministic cap is {cap}"
            ),
        }
    }
}

impl core::error::Error for SheafBuildError {}

/// Outward enclosure of a chart sample's signed distance from its own error
/// certificate. Only well-formed rigorous claims are usable: estimates,
/// no-claims, and malformed rigorous certificates poison that sample into the
/// whole extended real line. It cannot contribute positive authority, though
/// it may coexist with an independently proven violation and aggregate `Fail`.
fn sample_interval(s: &ChartSample) -> Interval {
    match s.error.kind {
        NumericalKind::Exact
            if s.signed_distance.is_finite()
                && s.error.lo.is_finite()
                && s.error.hi.is_finite()
                && s.error.lo.to_bits() == s.signed_distance.to_bits()
                && s.error.hi.to_bits() == s.signed_distance.to_bits() =>
        {
            Interval::point(s.signed_distance)
        }
        NumericalKind::Enclosure
            if s.signed_distance.is_finite()
                && s.error.lo.is_finite()
                && s.error.hi.is_finite()
                && s.error.lo <= s.signed_distance
                && s.signed_distance <= s.error.hi =>
        {
            Interval::new(s.error.lo, s.error.hi)
        }
        NumericalKind::Exact
        | NumericalKind::Enclosure
        | NumericalKind::Estimate
        | NumericalKind::NoClaim => Interval::WHOLE,
    }
}

fn discover_triples(
    n_patches: usize,
    interfaces: &[Interface],
    cx: &Cx<'_>,
) -> Result<Vec<TripleCell>, SheafBuildError> {
    let mut degrees = Vec::new();
    degrees
        .try_reserve_exact(n_patches)
        .map_err(|_| SheafBuildError::ResourceExhausted {
            stage: "triple-degree-table",
        })?;
    degrees.resize(n_patches, 0usize);
    for (completed_edges, interface) in interfaces.iter().enumerate() {
        if completed_edges.is_multiple_of(256) {
            cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                stage: "triple-discovery",
                patches: None,
                completed_work: completed_edges,
                unit: SheafBuildProgressUnit::Edges,
            })?;
        }
        let (a, b) = interface.patches;
        if a >= n_patches || b >= n_patches || a >= b {
            return Err(SheafBuildError::InternalInvariant {
                stage: "triple-invalid-builder-edge",
            });
        }
        degrees[a] = degrees[a].saturating_add(1);
        degrees[b] = degrees[b].saturating_add(1);
    }
    let mut adjacency = Vec::new();
    adjacency
        .try_reserve_exact(n_patches)
        .map_err(|_| SheafBuildError::ResourceExhausted {
            stage: "triple-adjacency-table",
        })?;
    for degree in degrees {
        let mut neighbors = Vec::new();
        neighbors
            .try_reserve_exact(degree)
            .map_err(|_| SheafBuildError::ResourceExhausted {
                stage: "triple-adjacency-row",
            })?;
        adjacency.push(neighbors);
    }
    for interface in interfaces {
        let (a, b) = interface.patches;
        adjacency[a].push(b);
        adjacency[b].push(a);
    }

    let mut triples = Vec::new();
    let mut inspected = 0usize;
    for interface in interfaces {
        let (a, b) = interface.patches;
        let ab_samples = interface.samples.len();
        let (a_neighbors, b_neighbors) = (&adjacency[a], &adjacency[b]);
        // Iterate the smaller neighbor set explicitly so every membership
        // probe is counted and cancellation-checkable. A hidden set-
        // intersection iterator would conceal work spent skipping non-common
        // neighbors in a dense triangle-free graph.
        let (probe, lookup) = if a_neighbors.len() <= b_neighbors.len() {
            (a_neighbors, b_neighbors)
        } else {
            (b_neighbors, a_neighbors)
        };
        for &c in probe {
            if inspected >= SHEAF_MAX_TRIPLE_CANDIDATES {
                return Err(SheafBuildError::TripleWorkLimit {
                    completed_probes: inspected,
                    cap: SHEAF_MAX_TRIPLE_CANDIDATES,
                });
            }
            if inspected.is_multiple_of(256) {
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "triple-discovery",
                    patches: None,
                    completed_work: inspected,
                    unit: SheafBuildProgressUnit::NeighborProbes,
                })?;
            }
            inspected += 1;
            if c <= b || lookup.binary_search(&c).is_err() {
                continue;
            }
            let Ok(bc_index) = interfaces.binary_search_by_key(&(b, c), |edge| edge.patches) else {
                continue;
            };
            let Ok(ac_index) = interfaces.binary_search_by_key(&(a, c), |edge| edge.patches) else {
                continue;
            };
            let bc_samples = interfaces[bc_index].samples.len();
            let ac_samples = interfaces[ac_index].samples.len();
            triples
                .try_reserve(1)
                .map_err(|_| SheafBuildError::ResourceExhausted {
                    stage: "triple-cells",
                })?;
            triples.push(TripleCell {
                patches: (a, b, c),
                samples: ab_samples.min(bc_samples).min(ac_samples),
            });
        }
    }
    cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
        stage: "triple-discovery",
        patches: None,
        completed_work: inspected,
        unit: SheafBuildProgressUnit::NeighborProbes,
    })?;
    Ok(triples)
}

impl SheafComplex {
    /// Build the complex from charts: interface discovery via
    /// support-overlap, shared-surface samples via the zero band of BOTH
    /// charts, plus deterministic pairwise-interface clique completions. The
    /// latter are algebraic candidate cells, not verified common triple
    /// overlaps.
    pub fn from_charts(
        charts: &[&dyn Chart],
        cx: &Cx<'_>,
    ) -> Result<AdmittedSheafComplex, SheafBuildError> {
        Self::from_charts_with_clip(charts, None, cx)
    }

    /// Build interface evidence inside an explicit finite clip. Pairs outside
    /// the clip are skipped; invalid supports and unresolved unbounded shared
    /// domains are structured refusals.
    pub fn from_charts_clipped(
        charts: &[&dyn Chart],
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<AdmittedSheafComplex, SheafBuildError> {
        Self::from_charts_with_clip(charts, Some(clip), cx)
    }

    #[allow(clippy::too_many_lines)] // One ordered discovery, sampling, cancellation, and finalize transaction.
    fn from_charts_with_clip(
        charts: &[&dyn Chart],
        clip: Option<Aabb>,
        cx: &Cx<'_>,
    ) -> Result<AdmittedSheafComplex, SheafBuildError> {
        cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
            stage: "admission",
            patches: None,
            completed_work: 0,
            unit: SheafBuildProgressUnit::AdmissionChecks,
        })?;
        let n = charts.len();
        if n > SHEAF_MAX_CHARTS {
            return Err(SheafBuildError::BuildWorkLimit {
                stage: "chart-support-snapshot",
                requested: n as u128,
                cap: SHEAF_MAX_CHARTS,
            });
        }
        let pair_candidates = (n as u128).saturating_mul(n.saturating_sub(1) as u128) / 2;
        if pair_candidates > SHEAF_MAX_PAIR_CANDIDATES as u128 {
            return Err(SheafBuildError::BuildWorkLimit {
                stage: "pair-support-preflight",
                requested: pair_candidates,
                cap: SHEAF_MAX_PAIR_CANDIDATES,
            });
        }
        if let Some(explicit_clip) = clip {
            SamplingDomain::resolve(Aabb::WHOLE_SPACE, Some(explicit_clip))
                .map_err(|source| SheafBuildError::SamplingClip { source })?;
        }

        let mut supports = Vec::new();
        supports
            .try_reserve_exact(n)
            .map_err(|_| SheafBuildError::ResourceExhausted {
                stage: "chart-support-snapshot",
            })?;
        for (chart, source) in charts.iter().enumerate() {
            cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                stage: "chart-support-snapshot",
                patches: None,
                completed_work: chart,
                unit: SheafBuildProgressUnit::Charts,
            })?;
            let support = source.support();
            SamplingDomain::validate_support(support)
                .map_err(|source| SheafBuildError::ChartSupport { chart, source })?;
            supports.push(support);
        }

        let evaluations_per_pair = (SAMPLES_PER_INTERFACE as u128)
            .saturating_mul(64)
            .saturating_mul(2);
        let maximum_overlap_pairs = (SHEAF_MAX_INTERFACE_EVALUATIONS as u128
            / evaluations_per_pair)
            .min(SHEAF_MAX_RETAINED_INTERFACE_SAMPLES as u128 / SAMPLES_PER_INTERFACE as u128);
        let domain_capacity =
            usize::try_from(pair_candidates.min(maximum_overlap_pairs)).map_err(|_| {
                SheafBuildError::BuildWorkLimit {
                    stage: "pair-domain-preflight",
                    requested: pair_candidates.min(maximum_overlap_pairs),
                    cap: usize::MAX,
                }
            })?;
        let mut domains = Vec::new();
        domains.try_reserve_exact(domain_capacity).map_err(|_| {
            SheafBuildError::ResourceExhausted {
                stage: "pair-domain-preflight",
            }
        })?;
        let mut pair_probes = 0usize;
        for u in 0..n {
            for v in (u + 1)..n {
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "pair-domain-preflight",
                    patches: Some((u, v)),
                    completed_work: pair_probes,
                    unit: SheafBuildProgressUnit::PairCandidates,
                })?;
                pair_probes += 1;
                let Some(shared_support) = supports[u].intersection(&supports[v]) else {
                    continue;
                };
                let domain = match SamplingDomain::resolve(shared_support, clip) {
                    Ok(domain) => domain,
                    Err(
                        SamplingDomainError::EmptyIntersection
                        | SamplingDomainError::DegenerateDomain { .. },
                    ) => continue,
                    Err(source) => {
                        return Err(SheafBuildError::SamplingDomain {
                            patches: (u, v),
                            source,
                        });
                    }
                };
                let requested_pairs = domains.len() as u128 + 1;
                let requested_evaluations = requested_pairs.saturating_mul(evaluations_per_pair);
                if requested_evaluations > SHEAF_MAX_INTERFACE_EVALUATIONS as u128 {
                    return Err(SheafBuildError::BuildWorkLimit {
                        stage: "interface-sampling-evaluations",
                        requested: requested_evaluations,
                        cap: SHEAF_MAX_INTERFACE_EVALUATIONS,
                    });
                }
                let requested_samples =
                    requested_pairs.saturating_mul(SAMPLES_PER_INTERFACE as u128);
                if requested_samples > SHEAF_MAX_RETAINED_INTERFACE_SAMPLES as u128 {
                    return Err(SheafBuildError::BuildWorkLimit {
                        stage: "retained-interface-samples",
                        requested: requested_samples,
                        cap: SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
                    });
                }
                domains.push((u, v, domain));
            }
        }

        let mut interfaces = Vec::new();
        interfaces.try_reserve_exact(domains.len()).map_err(|_| {
            SheafBuildError::ResourceExhausted {
                stage: "retained-interfaces",
            }
        })?;
        for (u, v, domain) in domains {
            let shared = domain.bounds();
            let spans = domain.spans();
            let diag = domain.diagonal();
            let band = BAND_FRACTION * diag;
            let mut state = box_seed(&shared);
            let mut samples = Vec::new();
            samples
                .try_reserve_exact(SAMPLES_PER_INTERFACE)
                .map_err(|_| SheafBuildError::ResourceExhausted {
                    stage: "interface-samples",
                })?;
            // Rejection-sample the shared zero band (bounded draws).
            for draw_index in 0..SAMPLES_PER_INTERFACE * 64 {
                if samples.len() >= SAMPLES_PER_INTERFACE {
                    break;
                }
                let completed_draws = draw_index;
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "interface-sampling",
                    patches: Some((u, v)),
                    completed_work: completed_draws,
                    unit: SheafBuildProgressUnit::InterfaceDraws,
                })?;
                let p = Point3::new(
                    shared.min.x + lcg(&mut state) * spans.x,
                    shared.min.y + lcg(&mut state) * spans.y,
                    shared.min.z + lcg(&mut state) * spans.z,
                );
                let su = charts[u].eval(p, cx);
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "interface-sampling",
                    patches: Some((u, v)),
                    completed_work: completed_draws,
                    unit: SheafBuildProgressUnit::InterfaceDraws,
                })?;
                if !su.signed_distance.is_finite() {
                    return Err(SheafBuildError::NonFiniteSample {
                        patches: (u, v),
                        chart: u,
                        point_bits: [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()],
                        value_bits: su.signed_distance.to_bits(),
                        completed_draws,
                    });
                }
                let sv = charts[v].eval(p, cx);
                let completed_draws = draw_index + 1;
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "interface-sampling",
                    patches: Some((u, v)),
                    completed_work: completed_draws,
                    unit: SheafBuildProgressUnit::InterfaceDraws,
                })?;
                if !sv.signed_distance.is_finite() {
                    return Err(SheafBuildError::NonFiniteSample {
                        patches: (u, v),
                        chart: v,
                        point_bits: [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()],
                        value_bits: sv.signed_distance.to_bits(),
                        completed_draws,
                    });
                }
                if su.signed_distance.abs() <= band && sv.signed_distance.abs() <= band {
                    samples.push(InterfaceSample {
                        point: p,
                        values: [sample_interval(&su), sample_interval(&sv)],
                    });
                }
            }
            if !samples.is_empty() {
                interfaces.push(Interface {
                    patches: (u, v),
                    samples,
                });
            }
        }
        let triples = discover_triples(n, &interfaces, cx)?;
        cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
            stage: "finalize",
            patches: None,
            completed_work: triples.len(),
            unit: SheafBuildProgressUnit::TripleCells,
        })?;
        Ok(AdmittedSheafComplex {
            inner: SheafComplex {
                sampling_clip: clip,
                n_patches: n,
                interfaces,
                triples,
            },
        })
    }

    /// Assemble the sampled restriction incidence (samples × patches) with
    /// ±1 entries: one row per interface sample, `+1` on patch v's slot and
    /// `−1` on patch u's. These sample rows are not dimension-compatible with
    /// [`Self::delta1`]; use [`Self::delta0_edges`] for the edge-level cochain
    /// map in the bitwise `δ¹δ⁰ = 0` identity.
    pub fn delta0(&self) -> Result<Csr, SheafAlgebraError> {
        let rows = self.validate_algebra_shape()?;
        let mut coo = Coo::new(rows, self.n_patches);
        let mut r = 0usize;
        for iface in &self.interfaces {
            for _ in &iface.samples {
                coo.push(r, iface.patches.0, -1.0);
                coo.push(r, iface.patches.1, 1.0);
                r += 1;
            }
        }
        Ok(coo.assemble())
    }

    /// Assemble the edge-level δ¹ (triples × interfaces) with ±1 entries
    /// per the oriented
    /// triangle boundary: for triple (a,b,c) with edges e_ab, e_bc, e_ac:
    /// `+e_ab + e_bc − e_ac` (edge-level stalks: one column per edge).
    /// Malformed raw complexes return [`SheafAlgebraError`] before sparse
    /// assembly; admitted complexes have all three indexed edges by construction.
    pub fn delta1(&self) -> Result<Csr, SheafAlgebraError> {
        let _ = self.validate_algebra_shape()?;
        let mut coo = Coo::new(self.triples.len(), self.interfaces.len());
        for (t, triple) in self.triples.iter().enumerate() {
            let (a, b, c) = triple.patches;
            let indices = [
                self.interfaces
                    .binary_search_by_key(&(a.min(b), a.max(b)), |edge| edge.patches)
                    .ok(),
                self.interfaces
                    .binary_search_by_key(&(b.min(c), b.max(c)), |edge| edge.patches)
                    .ok(),
                self.interfaces
                    .binary_search_by_key(&(a.min(c), a.max(c)), |edge| edge.patches)
                    .ok(),
            ];
            if let [Some(ab), Some(bc), Some(ac)] = indices {
                coo.push(t, ab, 1.0);
                coo.push(t, bc, 1.0);
                coo.push(t, ac, -1.0);
            }
        }
        Ok(coo.assemble())
    }

    /// Edge-level δ⁰ (edges × patches, one row per INTERFACE): the
    /// companion of [`Self::delta1`] for the bitwise δδ = 0 identity.
    pub fn delta0_edges(&self) -> Result<Csr, SheafAlgebraError> {
        let _ = self.validate_algebra_shape()?;
        let mut coo = Coo::new(self.interfaces.len(), self.n_patches);
        for (r, iface) in self.interfaces.iter().enumerate() {
            coo.push(r, iface.patches.0, -1.0);
            coo.push(r, iface.patches.1, 1.0);
        }
        Ok(coo.assemble())
    }

    fn validate_algebra_shape(&self) -> Result<usize, SheafAlgebraError> {
        for (stage, requested, cap) in [
            ("patches", self.n_patches as u128, SHEAF_MAX_CHARTS),
            (
                "interfaces",
                self.interfaces.len() as u128,
                SHEAF_MAX_PAIR_CANDIDATES,
            ),
            (
                "triples",
                self.triples.len() as u128,
                SHEAF_MAX_TRIPLE_CANDIDATES,
            ),
        ] {
            if requested > cap as u128 {
                return Err(SheafAlgebraError::WorkLimit {
                    stage,
                    requested,
                    cap,
                });
            }
        }
        let mut samples = 0usize;
        for interface in &self.interfaces {
            samples = samples.checked_add(interface.samples.len()).ok_or(
                SheafAlgebraError::WorkLimit {
                    stage: "interface-samples",
                    requested: u128::MAX,
                    cap: SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
                },
            )?;
            if samples > SHEAF_MAX_RETAINED_INTERFACE_SAMPLES {
                return Err(SheafAlgebraError::WorkLimit {
                    stage: "interface-samples",
                    requested: samples as u128,
                    cap: SHEAF_MAX_RETAINED_INTERFACE_SAMPLES,
                });
            }
        }
        if !self.structure_is_valid() {
            return Err(SheafAlgebraError::InvalidStructure);
        }
        Ok(samples)
    }

    pub(crate) fn structure_is_valid(&self) -> bool {
        if self.n_patches == 0 {
            return false;
        }
        if self
            .sampling_clip
            .is_some_and(|clip| SamplingDomain::resolve(Aabb::WHOLE_SPACE, Some(clip)).is_err())
        {
            return false;
        }
        let mut previous_edge = None;
        for interface in &self.interfaces {
            let (u, v) = interface.patches;
            if u >= v
                || v >= self.n_patches
                || interface.samples.is_empty()
                || previous_edge.is_some_and(|previous| previous >= (u, v))
                || interface.samples.iter().any(|sample| {
                    !finite_point(sample.point)
                        || self
                            .sampling_clip
                            .is_some_and(|clip| !clip.contains(sample.point))
                })
            {
                return false;
            }
            previous_edge = Some((u, v));
        }

        let mut previous_triple = None;
        for triple in &self.triples {
            let (a, b, c) = triple.patches;
            if a >= b
                || b >= c
                || c >= self.n_patches
                || triple.samples == 0
                || previous_triple.is_some_and(|previous| previous >= (a, b, c))
            {
                return false;
            }
            // The interface loop above has already established strict edge
            // ordering, so each lookup is logarithmic rather than rescanning
            // every interface three times for every triple.
            let edge_samples = |u: usize, v: usize| {
                self.interfaces
                    .binary_search_by_key(&(u, v), |interface| interface.patches)
                    .ok()
                    .map(|index| self.interfaces[index].samples.len())
            };
            let expected = match (edge_samples(a, b), edge_samples(a, c), edge_samples(b, c)) {
                (Some(ab), Some(ac), Some(bc)) => ab.min(ac).min(bc),
                _ => return false,
            };
            if triple.samples != expected {
                return false;
            }
            previous_triple = Some((a, b, c));
        }
        true
    }

    /// Per-interface mismatch assessment. The VERDICT bits come from
    /// fs-ivl's sound predicates (`encloses`/`contains`). Reported magnitudes
    /// are aggregated directly from the intervals' outward endpoints; an
    /// indeterminate interval retains an infinite upper report and cannot
    /// authorize numerical evidence.
    #[must_use]
    fn mismatch_bounds(&self) -> Result<Vec<InterfaceBound>, SheafAlgebraError> {
        let mut bounds = Vec::new();
        bounds
            .try_reserve_exact(self.interfaces.len())
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "interface-mismatch-bounds",
            })?;
        for iface in &self.interfaces {
            let (u, v) = iface.patches;
            let valid_interface = u < v
                && v < self.n_patches
                && !iface.samples.is_empty()
                && iface
                    .samples
                    .iter()
                    .all(|sample| finite_point(sample.point));
            let mut determinate = valid_interface;
            let mut has_determinate_sample = false;
            let mut lo = 0.0f64;
            let mut hi = if valid_interface { 0.0 } else { f64::INFINITY };
            for s in &iface.samples {
                let d = (s.values[1] - s.values[0]).abs();
                if !(d.lo().is_finite() && d.hi().is_finite()) {
                    determinate = false;
                    hi = f64::INFINITY;
                    continue;
                }
                has_determinate_sample = true;
                lo = lo.max(d.lo());
                hi = hi.max(d.hi());
            }
            bounds.push(InterfaceBound {
                patches: iface.patches,
                determinate,
                has_determinate_sample: valid_interface && has_determinate_sample,
                lo_report: lo,
                hi_report: hi,
            });
        }
        Ok(bounds)
    }

    /// Least-squares section solve: per-patch gauge offsets minimizing the
    /// mean-square mismatch (graph-Laplacian normal equations, one smallest-
    /// index gauge root pinned in every connected component; deterministic
    /// Gauss–Seidel). Returns (offsets, raw ms mismatch, residual ms mismatch).
    pub fn section_solve(&self) -> Result<(Vec<f64>, f64, f64), SheafAlgebraError> {
        let _ = self.validate_algebra_shape()?;
        let n = self.n_patches;
        let mut offsets = Vec::new();
        offsets
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-offsets",
            })?;
        offsets.resize(n, 0.0f64);
        // Edge means of the midpoint mismatch.
        let mut degrees = Vec::new();
        degrees
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-degrees",
            })?;
        degrees.resize(n, 0usize);
        for iface in &self.interfaces {
            degrees[iface.patches.0] = degrees[iface.patches.0].saturating_add(1);
            degrees[iface.patches.1] = degrees[iface.patches.1].saturating_add(1);
        }
        let mut incident = Vec::new();
        incident
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-incidence",
            })?;
        for degree in degrees {
            let mut row = Vec::new();
            row.try_reserve_exact(degree)
                .map_err(|_| SheafAlgebraError::ResourceExhausted {
                    stage: "section-incidence-row",
                })?;
            incident.push(row);
        }
        for (interface_index, iface) in self.interfaces.iter().enumerate() {
            let (u, v) = iface.patches;
            if u >= n || v >= n || u == v || iface.samples.is_empty() {
                continue;
            }
            let mut sum = 0.0;
            for (sample_index, s) in iface.samples.iter().enumerate() {
                if s.values.iter().any(|value| {
                    !(value.lo().is_finite() && value.hi().is_finite() && value.lo() <= value.hi())
                }) {
                    return Err(SheafAlgebraError::IndeterminateSampleValue {
                        interface: interface_index,
                        sample: sample_index,
                    });
                }
                let left = s.values[0].midpoint();
                let right = s.values[1].midpoint();
                if !(left.is_finite() && right.is_finite()) {
                    return Err(SheafAlgebraError::IndeterminateSampleValue {
                        interface: interface_index,
                        sample: sample_index,
                    });
                }
                let mismatch = right - left;
                if !mismatch.is_finite() {
                    return Err(SheafAlgebraError::NumericalOverflow {
                        stage: "section-edge-mismatch",
                    });
                }
                sum += mismatch;
                if !sum.is_finite() {
                    return Err(SheafAlgebraError::NumericalOverflow {
                        stage: "section-edge-sum",
                    });
                }
            }
            let count = s_len(iface);
            incident[u].push((v, sum, count));
            incident[v].push((u, -sum, count));
        }
        // Fix one deterministic gauge root (the smallest patch index) in every
        // connected component, including isolated patches. This removes the
        // otherwise implicit iteration/order-dependent null mode on components
        // not containing patch 0.
        let mut pinned = Vec::new();
        pinned
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-component-roots",
            })?;
        pinned.resize(n, false);
        let mut visited = Vec::new();
        visited
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-component-visited",
            })?;
        visited.resize(n, false);
        let mut stack = Vec::new();
        stack
            .try_reserve_exact(n)
            .map_err(|_| SheafAlgebraError::ResourceExhausted {
                stage: "section-component-stack",
            })?;
        for root in 0..n {
            if visited[root] {
                continue;
            }
            pinned[root] = true;
            visited[root] = true;
            stack.push(root);
            while let Some(patch) = stack.pop() {
                for (neighbor, _, _) in &incident[patch] {
                    if !visited[*neighbor] {
                        visited[*neighbor] = true;
                        stack.push(*neighbor);
                    }
                }
            }
        }
        let raw_ms = sample_mean_square(&self.interfaces, &offsets)?;
        for _ in 0..200 {
            for p in 0..n {
                if pinned[p] {
                    continue;
                }
                // Optimal c_p given the rest: weighted average balance.
                let mut num = 0.0f64;
                let mut den = 0.0f64;
                for (neighbor, sum, count) in &incident[p] {
                    #[allow(clippy::cast_precision_loss)]
                    let w = *count as f64;
                    num += sum + w * offsets[*neighbor];
                    den += w;
                    if !(num.is_finite() && den.is_finite()) {
                        return Err(SheafAlgebraError::NumericalOverflow {
                            stage: "section-gauss-seidel",
                        });
                    }
                }
                if den > 0.0 {
                    offsets[p] = num / den;
                    if !offsets[p].is_finite() {
                        return Err(SheafAlgebraError::NumericalOverflow {
                            stage: "section-offset-update",
                        });
                    }
                }
            }
        }
        let residual_ms = sample_mean_square(&self.interfaces, &offsets)?;
        Ok((offsets, raw_ms, residual_ms))
    }

    /// Assess raw public parts as non-authoritative diagnostics. Raw complexes
    /// can exercise incidence and mismatch algebra, but callers can construct
    /// every field themselves, so neither `Pass` nor `Fail` would be evidence.
    /// Use [`Self::from_charts`] to obtain an immutable
    /// [`AdmittedSheafComplex`] when sampled-interface authority is required.
    #[must_use]
    pub fn watertightness(&self, tol: f64) -> Evidence<SheafVerdict> {
        self.assess_sampled_agreement(tol, false)
    }

    fn assess_sampled_agreement(
        &self,
        tol: f64,
        admitted_builder_origin: bool,
    ) -> Evidence<SheafVerdict> {
        let shape_validation = self.validate_algebra_shape();
        let bounded_shape = !matches!(shape_validation, Err(SheafAlgebraError::WorkLimit { .. }));
        let bounds_result = if bounded_shape {
            self.mismatch_bounds()
        } else {
            Ok(Vec::new())
        };
        let bounds_available = bounds_result.is_ok();
        let bounds = bounds_result.unwrap_or_default();
        let valid_tolerance = tol.is_finite() && tol >= 0.0;
        let structure_is_valid = shape_validation.is_ok() && bounds_available;
        let worst_hi = if !structure_is_valid || bounds.is_empty() {
            // No sampled interface means no finite mismatch quantity was
            // measured. Infinity is the fail-closed QoI sentinel; retaining
            // zero here makes an Unknown/NoClaim verdict look numerically clean
            // to consumers that inspect only the scalar. Malformed public
            // structure likewise cannot publish a finite clean-looking QoI
            // merely because some fabricated bounds happened to be parseable.
            f64::INFINITY
        } else {
            bounds.iter().map(|b| b.hi_report).fold(0.0f64, f64::max)
        };
        let worst_lo = bounds.iter().map(|b| b.lo_report).fold(0.0f64, f64::max);
        let all_determinate = valid_tolerance
            && structure_is_valid
            && !bounds.is_empty()
            && bounds.iter().all(|bound| bound.determinate);
        // PASS requires at least one DISCOVERED interface whose samples all lie
        // inside [0, tol]. An empty `bounds` means NO interface was found (charts
        // are disjoint/gapped/near-tangent, or the geometry is empty) — the
        // interface-agreement check gathered NO evidence, so the honest verdict
        // is `Unknown`, not a positive `worst_mismatch = 0` PASS. `all()` on an
        // empty set is vacuously true; guarding on non-emptiness closes the
        // vacuous-truth false certificate (bead obnw).
        let all_pass = all_determinate && bounds.iter().all(|b| b.all_within(tol));
        let interface_violations: Vec<((usize, usize), f64)> = bounds
            .iter()
            .filter(|b| b.proven_leak(tol))
            .map(|b| (b.patches, b.lo_report))
            .collect();
        let verdict = if !admitted_builder_origin {
            SheafVerdict::Unknown {
                reported_bounds: bounds
                    .iter()
                    .map(|b| (b.patches, b.lo_report, b.hi_report))
                    .collect(),
            }
        } else if all_pass {
            SheafVerdict::Pass {
                worst_mismatch: worst_hi,
                margins: bounds.iter().map(|b| (b.patches, b.hi_report)).collect(),
            }
        } else if structure_is_valid && !interface_violations.is_empty() {
            let share = if all_determinate {
                match self.section_solve() {
                    Ok((_, raw, residual))
                        if raw.is_finite()
                            && residual.is_finite()
                            && raw > 0.0
                            && residual >= 0.0 =>
                    {
                        Some((1.0 - residual / raw).clamp(0.0, 1.0))
                    }
                    Ok((_, raw, residual)) if raw == 0.0 && residual == 0.0 => Some(0.0),
                    Ok(_) | Err(_) => None,
                }
            } else {
                None
            };
            SheafVerdict::Fail {
                interface_violations,
                gauge_fit_share: share,
            }
        } else {
            SheafVerdict::Unknown {
                reported_bounds: bounds
                    .iter()
                    .filter(|b| !b.determinate || !b.all_within(tol))
                    .map(|b| (b.patches, b.lo_report, b.hi_report))
                    .collect(),
            }
        };
        let mut canon = LegacyProvenanceWriter::new();
        let _ = write!(
            canon,
            "sheaf-sampled-agreement;schema=4;origin={};patches={};interfaces={};triples={};tol={:016x};structure_valid={structure_is_valid}",
            if admitted_builder_origin {
                "chart-sampling-builder"
            } else {
                "raw-public-parts"
            },
            self.n_patches,
            self.interfaces.len(),
            self.triples.len(),
            tol.to_bits(),
        );
        match self.sampling_clip {
            None => {
                let _ = canon.write_str(";sampling_clip=none");
            }
            Some(clip) => {
                let _ = write!(
                    canon,
                    ";sampling_clip=some:{:016x},{:016x},{:016x},{:016x},{:016x},{:016x}",
                    clip.min.x.to_bits(),
                    clip.min.y.to_bits(),
                    clip.min.z.to_bits(),
                    clip.max.x.to_bits(),
                    clip.max.y.to_bits(),
                    clip.max.z.to_bits()
                );
            }
        }
        if bounded_shape {
            for interface in &self.interfaces {
                let _ = write!(
                    canon,
                    ";interface={}-{};samples={}",
                    interface.patches.0,
                    interface.patches.1,
                    interface.samples.len()
                );
                for sample in &interface.samples {
                    let _ = write!(
                        canon,
                        ";sample={:016x},{:016x},{:016x}:{:016x},{:016x}:{:016x},{:016x}",
                        sample.point.x.to_bits(),
                        sample.point.y.to_bits(),
                        sample.point.z.to_bits(),
                        sample.values[0].lo().to_bits(),
                        sample.values[0].hi().to_bits(),
                        sample.values[1].lo().to_bits(),
                        sample.values[1].hi().to_bits(),
                    );
                }
            }
            for triple in &self.triples {
                let _ = write!(
                    canon,
                    ";triple={}-{}-{}:{}",
                    triple.patches.0, triple.patches.1, triple.patches.2, triple.samples
                );
            }
        } else {
            let _ = canon.write_str(";raw-payload=omitted-over-work-limit");
        }
        for b in &bounds {
            let _ = write!(
                canon,
                ";bound={}-{}:{:016x}:{:016x}:within={}:leak={}:determinate={}",
                b.patches.0,
                b.patches.1,
                b.lo_report.to_bits(),
                b.hi_report.to_bits(),
                b.all_within(tol),
                b.proven_leak(tol),
                b.determinate,
            );
        }
        match &verdict {
            SheafVerdict::Pass {
                worst_mismatch,
                margins,
            } => {
                let _ = write!(canon, ";verdict=pass:{:016x}", worst_mismatch.to_bits());
                for (patches, margin) in margins {
                    let _ = write!(
                        canon,
                        ";margin={}-{}:{:016x}",
                        patches.0,
                        patches.1,
                        margin.to_bits()
                    );
                }
            }
            SheafVerdict::Fail {
                interface_violations,
                gauge_fit_share,
            } => {
                let _ = canon.write_str(";verdict=fail");
                for (patches, lower) in interface_violations {
                    let _ = write!(
                        canon,
                        ";violation={}-{}:{:016x}",
                        patches.0,
                        patches.1,
                        lower.to_bits()
                    );
                }
                match gauge_fit_share {
                    Some(share) => {
                        let _ = write!(canon, ";gauge_fit_share={:016x}", share.to_bits());
                    }
                    None => {
                        let _ = canon.write_str(";gauge_fit_share=none");
                    }
                }
            }
            SheafVerdict::Unknown { reported_bounds } => {
                let _ = canon.write_str(";verdict=unknown");
                for (patches, lower, upper) in reported_bounds {
                    let _ = write!(
                        canon,
                        ";reported-bound={}-{}:{:016x}:{:016x}",
                        patches.0,
                        patches.1,
                        lower.to_bits(),
                        upper.to_bits()
                    );
                }
            }
        }
        let numerical = if admitted_builder_origin
            && all_determinate
            && !matches!(&verdict, SheafVerdict::Unknown { .. })
            && worst_lo.is_finite()
            && worst_hi.is_finite()
        {
            NumericalCertificate::enclosure(worst_lo, worst_hi)
        } else {
            NumericalCertificate::no_claim()
        };
        Evidence {
            qoi: worst_hi,
            numerical,
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash(canon.finish()),
            adjoint_ref: None,
            value: verdict,
        }
    }
}

fn s_len(iface: &Interface) -> usize {
    iface.samples.len()
}

fn sample_mean_square(interfaces: &[Interface], offsets: &[f64]) -> Result<f64, SheafAlgebraError> {
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for (interface_index, interface) in interfaces.iter().enumerate() {
        let (u, v) = interface.patches;
        if u >= offsets.len() || v >= offsets.len() || u == v {
            continue;
        }
        for (sample_index, sample) in interface.samples.iter().enumerate() {
            if sample.values.iter().any(|value| {
                !(value.lo().is_finite() && value.hi().is_finite() && value.lo() <= value.hi())
            }) {
                return Err(SheafAlgebraError::IndeterminateSampleValue {
                    interface: interface_index,
                    sample: sample_index,
                });
            }
            let mismatch = sample.values[1].midpoint() - sample.values[0].midpoint();
            let gauged = mismatch + offsets[v] - offsets[u];
            let square = gauged * gauged;
            if !(mismatch.is_finite() && gauged.is_finite() && square.is_finite()) {
                return Err(SheafAlgebraError::NumericalOverflow {
                    stage: "section-mean-square",
                });
            }
            num += square;
            den += 1.0;
            if !(num.is_finite() && den.is_finite()) {
                return Err(SheafAlgebraError::NumericalOverflow {
                    stage: "section-mean-square-accumulation",
                });
            }
        }
    }
    Ok(if den > 0.0 { num / den } else { 0.0 })
}

/// Proven-outside to proven-outside sign-sequence diagnostic for input
/// validation and replay.
/// It is NOT an independent topology falsifier: because both endpoints are
/// required to be strictly outside, the sampled boolean sequence begins and
/// ends with the same sign and therefore has an even number of toggles by
/// construction. Authentic cross-examination needs certified oriented
/// intersections or winding/degree evidence. The sample count remains
/// work-capped, every endpoint and nominal field sample is validated, each
/// endpoint is excluded by every chart either through its declared support or a
/// rigorous positive signed-distance enclosure, convex
/// interpolation avoids an overflowing endpoint subtraction, and cancellation
/// is observed after every chart evaluation.
///
/// # Errors
/// [`OutsideRaySampleError`] when the inputs do not satisfy the finite outside-ray
/// preconditions, the public work cap would be exceeded, a chart produces a
/// non-finite value, or cancellation is observed.
pub fn validate_outside_ray_samples(
    charts: &[&dyn Chart],
    rays: &[(Point3, Point3)],
    steps: usize,
    cx: &Cx<'_>,
) -> Result<OutsideRaySampleReport, OutsideRaySampleError> {
    if charts.is_empty() {
        return Err(OutsideRaySampleError::EmptyCharts);
    }
    if rays.is_empty() {
        return Err(OutsideRaySampleError::EmptyRays);
    }
    if steps == 0 {
        return Err(OutsideRaySampleError::InvalidSteps { steps });
    }
    let requested = (rays.len() as u128)
        .saturating_mul((steps as u128).saturating_add(1))
        .saturating_mul(charts.len() as u128);
    if requested > OUTSIDE_RAY_MAX_EVALUATIONS as u128 {
        return Err(OutsideRaySampleError::WorkLimitExceeded {
            requested,
            cap: OUTSIDE_RAY_MAX_EVALUATIONS,
        });
    }

    for (ray, (start, end)) in rays.iter().copied().enumerate() {
        for (endpoint, point) in [(RayEndpoint::Start, start), (RayEndpoint::End, end)] {
            if !finite_point(point) {
                return Err(OutsideRaySampleError::NonFiniteEndpoint {
                    ray,
                    endpoint,
                    point,
                });
            }
        }
    }

    let mut completed_points = 0usize;
    let mut completed_chart_evaluations = 0usize;
    let mut total_toggles = 0usize;
    for (ri, (a, b)) in rays.iter().enumerate() {
        let mut crossings = 0usize;
        let mut prev_sign = None;
        for k in 0..=steps {
            checkpoint_outside_ray_samples(cx, ri, completed_points, completed_chart_evaluations)?;
            let p = ray_sample_point(*a, *b, k, steps);
            if !finite_point(p) {
                return Err(OutsideRaySampleError::NonRepresentableSamplePoint {
                    ray: ri,
                    step: k,
                    point: p,
                });
            }
            let mut d = f64::INFINITY;
            let mut unproven_endpoint = None;
            for (chart_index, chart) in charts.iter().enumerate() {
                let sample = chart.eval(p, cx);
                let value = sample.signed_distance;
                completed_chart_evaluations = completed_chart_evaluations.saturating_add(1);
                checkpoint_outside_ray_samples(
                    cx,
                    ri,
                    completed_points,
                    completed_chart_evaluations,
                )?;
                if !value.is_finite() {
                    return Err(OutsideRaySampleError::NonFiniteSample {
                        ray: ri,
                        step: k,
                        chart: chart_index,
                        value,
                    });
                }
                d = d.min(value);
                if k == 0 || k == steps {
                    let support = chart.support();
                    let excluded_by_support = support.is_well_formed() && !support.contains(p);
                    let interval = sample_interval(&sample);
                    if !excluded_by_support && !(interval.lo().is_finite() && interval.lo() > 0.0) {
                        unproven_endpoint.get_or_insert((chart_index, sample));
                    }
                }
            }
            if k == 0 && d <= 0.0 {
                return Err(OutsideRaySampleError::EndpointNotOutside {
                    ray: ri,
                    endpoint: RayEndpoint::Start,
                    min_signed_distance: d,
                });
            }
            if k == steps && d <= 0.0 {
                return Err(OutsideRaySampleError::EndpointNotOutside {
                    ray: ri,
                    endpoint: RayEndpoint::End,
                    min_signed_distance: d,
                });
            }
            if let Some((chart, sample)) = unproven_endpoint {
                return Err(OutsideRaySampleError::EndpointOutsideUnproven {
                    ray: ri,
                    endpoint: if k == 0 {
                        RayEndpoint::Start
                    } else {
                        RayEndpoint::End
                    },
                    chart,
                    nominal: sample.signed_distance,
                    certificate: sample.error,
                });
            }
            let sign = d < 0.0;
            if let Some(ps) = prev_sign
                && ps != sign
            {
                crossings += 1;
            }
            prev_sign = Some(sign);
            completed_points = completed_points.saturating_add(1);
        }
        debug_assert!(crossings.is_multiple_of(2));
        total_toggles = total_toggles.saturating_add(crossings);
    }
    checkpoint_outside_ray_samples(
        cx,
        rays.len(),
        completed_points,
        completed_chart_evaluations,
    )?;
    Ok(OutsideRaySampleReport {
        rays: rays.len(),
        sample_points: completed_points,
        chart_evaluations: completed_chart_evaluations,
        toggles: total_toggles,
    })
}

fn checkpoint_outside_ray_samples(
    cx: &Cx<'_>,
    completed_rays: usize,
    completed_points: usize,
    completed_chart_evaluations: usize,
) -> Result<(), OutsideRaySampleError> {
    cx.checkpoint()
        .map_err(|_| OutsideRaySampleError::Cancelled {
            completed_rays,
            completed_points,
            completed_chart_evaluations,
        })
}

fn finite_point(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn ray_sample_point(start: Point3, end: Point3, step: usize, steps: usize) -> Point3 {
    if step == 0 {
        return start;
    }
    if step == steps {
        return end;
    }
    #[allow(clippy::cast_precision_loss)]
    let t = step as f64 / steps as f64;
    let lerp = |a: f64, b: f64| a.mul_add(1.0 - t, b * t);
    Point3::new(
        lerp(start.x, end.x),
        lerp(start.y, end.y),
        lerp(start.z, end.z),
    )
}

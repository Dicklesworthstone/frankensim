//! SHEAF CERTIFICATES (plan §7.3, Bet 11): the cellular-sheaf
//! watertightness certificate for multi-representation models. A model
//! whose patches live in different charts is a cellular sheaf over the
//! patch-adjacency complex: stalks are per-patch sample spaces,
//! restriction maps are selection operators onto shared interface
//! samples, and GLOBAL CONSISTENCY is the existence of a section. The
//! watertightness certificate is an INTERVAL-VERIFIED bound
//! `‖δs‖∞ ≤ tol` on the interface mismatch cocycle; when the coboundary
//! cannot be driven below tolerance, the H¹ obstruction is reported WITH
//! THE OFFENDING INTERFACE CELLS ATTACHED — exactly the diagnostic an
//! agent needs to fix a leaky model.
//!
//! The construction is finite linear algebra: δ⁰ and δ¹ assemble as
//! sparse matrices with entries in {−1, 0, +1} (restrictions are point
//! samplers), so `δ¹·δ⁰ = 0` holds BITWISE — small-integer f64
//! arithmetic is exact. The least-squares section solve (per-patch gauge
//! offsets over the adjacency Laplacian) splits the mismatch into a
//! reconcilable coboundary part and the structural residual — the same
//! split Proposal 10's merge semantics reuses unmodified.

use crate::{Aabb, Chart, ChartSample, Point3, SamplingDomain, SamplingDomainError};
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, NumericalKind, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate,
};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_sparse::{Coo, Csr};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Samples drawn per pairwise interface.
pub const SAMPLES_PER_INTERFACE: usize = 32;

/// Zero-band half-width as a fraction of the overlap-box diagonal:
/// a point belongs to the shared surface region when BOTH charts place
/// it within this band of their zero set.
pub const BAND_FRACTION: f64 = 0.05;

/// Maximum number of chart evaluations admitted by one ray-parity probe.
///
/// The falsifier is deliberately a bounded diagnostic, not an unbounded
/// marching workload hidden behind a sample-count argument.
pub const RAY_PARITY_MAX_EVALUATIONS: usize = 1_048_576;

/// Maximum number of fully connected triple candidates admitted during one
/// sheaf build.
pub const SHEAF_MAX_TRIPLE_CANDIDATES: usize = 1_048_576;

/// Endpoint of a ray named by a structured parity refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RayEndpoint {
    /// Segment start.
    Start,
    /// Segment end.
    End,
}

/// Why the independent ray-parity falsifier could not return a verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum RayParityError {
    /// A union-of-charts model needs at least one presentation.
    EmptyCharts,
    /// A falsification run with no rays gathers no evidence.
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
    /// The parity theorem requires both endpoints to be strictly outside the
    /// union model.
    EndpointNotOutside {
        /// Ray index.
        ray: usize,
        /// Which endpoint violated the precondition.
        endpoint: RayEndpoint,
        /// Minimum signed field across the charts.
        min_signed_distance: f64,
    },
    /// Cancellation was observed before a verdict could be published.
    Cancelled {
        /// Rays fully classified before cancellation.
        completed_rays: usize,
        /// Ray points fully evaluated before cancellation.
        completed_points: usize,
        /// Individual chart evaluations completed before cancellation.
        completed_chart_evaluations: usize,
    },
}

impl core::fmt::Display for RayParityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyCharts => write!(f, "ray parity refused: the chart set is empty"),
            Self::EmptyRays => write!(f, "ray parity refused: the ray set is empty"),
            Self::InvalidSteps { steps } => {
                write!(f, "ray parity refused: steps must be positive, got {steps}")
            }
            Self::WorkLimitExceeded { requested, cap } => write!(
                f,
                "ray parity refused: {requested} chart evaluations exceed the public cap {cap}"
            ),
            Self::NonFiniteEndpoint {
                ray,
                endpoint,
                point,
            } => write!(
                f,
                "ray parity refused: ray {ray} {endpoint:?} endpoint is non-finite: {point:?}"
            ),
            Self::NonRepresentableSamplePoint { ray, step, point } => write!(
                f,
                "ray parity refused: ray {ray} sample {step} is not representable: {point:?}"
            ),
            Self::NonFiniteSample {
                ray,
                step,
                chart,
                value,
            } => write!(
                f,
                "ray parity refused: chart {chart} returned non-finite value {value} at ray {ray} sample {step}"
            ),
            Self::EndpointNotOutside {
                ray,
                endpoint,
                min_signed_distance,
            } => write!(
                f,
                "ray parity refused: ray {ray} {endpoint:?} endpoint is not strictly outside (minimum field {min_signed_distance})"
            ),
            Self::Cancelled {
                completed_rays,
                completed_points,
                completed_chart_evaluations,
            } => write!(
                f,
                "ray parity cancelled after {completed_rays} rays, {completed_points} points, and {completed_chart_evaluations} chart evaluations"
            ),
        }
    }
}

impl core::error::Error for RayParityError {}

/// One shared interface sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceSample {
    /// The world-space point.
    pub point: Point3,
    /// Enclosures of the two charts' signed distances at the point
    /// (outward bounds from each chart's own error certificate).
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

/// A triple junction (2-cell): three patches with a common overlap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TripleCell {
    /// Patch indices (sorted).
    pub patches: (usize, usize, usize),
    /// Sample count shared by all three pairwise interfaces.
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
    /// Triple junctions.
    pub triples: Vec<TripleCell>,
}

/// One interface's assessed mismatch. Verdict bits are predicate-sound and
/// reported magnitudes come directly from outward interval endpoints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceBound {
    /// The patch pair.
    pub patches: (usize, usize),
    /// Every sample's |mismatch| enclosure lies inside `[0, tol]`.
    pub all_within: bool,
    /// Some sample's enclosure lies ENTIRELY above tol (proven leak).
    pub proven_leak: bool,
    /// Every mismatch interval had finite outward endpoints and the supplied
    /// tolerance was a finite non-negative value.
    pub determinate: bool,
    /// Reported worst lower bound.
    pub lo_report: f64,
    /// Reported worst upper bound.
    pub hi_report: f64,
}

/// The certificate verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum SheafVerdict {
    /// `‖δs‖∞ ≤ tol` with the enclosure upper bound as the margin —
    /// per-interface margins attached.
    Pass {
        /// Upper bound of the worst interface mismatch.
        worst_mismatch: f64,
        /// Per-interface (patch pair, mismatch upper bound).
        margins: Vec<((usize, usize), f64)>,
    },
    /// The H¹ obstruction: some interface's mismatch enclosure lies
    /// ENTIRELY above tolerance — a proven leak, localized.
    Fail {
        /// Offending interfaces: (patch pair, mismatch lower bound).
        obstruction: Vec<((usize, usize), f64)>,
        /// The reconcilable (coboundary) share of the raw mismatch in
        /// `[0, 1]`: how much a re-gauge would fix (the merge-semantics
        /// split — near 1 means gauge drift, near 0 means structure). `None`
        /// refuses the diagnostic when its least-squares arithmetic is not
        /// representable.
        coboundary_share: Option<f64>,
    },
    /// Enclosures straddle the tolerance: no sound claim either way
    /// (tighten chart certificates or the band and re-run).
    Unknown {
        /// Straddling interfaces: (patch pair, lower, upper).
        straddling: Vec<((usize, usize), f64, f64)>,
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

/// Why sheaf interface discovery could not safely sample a chart pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafBuildError {
    /// Cooperative cancellation was observed before publishing a complex.
    Cancelled {
        /// Stable build stage at which cancellation was observed.
        stage: &'static str,
        /// Pair being sampled, when cancellation occurred inside an interface.
        patches: Option<(usize, usize)>,
        /// Candidate points fully evaluated for that pair, or connected
        /// triple candidates inspected when `patches` is `None`.
        completed_draws: usize,
    },
    /// The requested finite sampling clip itself was not admissible. This is
    /// checked even when there are fewer than two charts or every chart pair
    /// is disjoint, so malformed caller input cannot silently succeed.
    SamplingClip {
        /// Exact clip admission failure.
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
    /// Fully connected triple discovery exceeded its deterministic work cap.
    TripleWorkLimit {
        /// Connected triple candidates encountered.
        candidates: usize,
        /// Public candidate cap.
        cap: usize,
    },
}

impl core::fmt::Display for SheafBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                stage,
                patches,
                completed_draws,
            } => write!(
                f,
                "sheaf build cancelled during {stage} for patches {patches:?} after \
                 {completed_draws} completed interface draws; no complex was published"
            ),
            Self::SamplingClip { source } => {
                write!(f, "sheaf explicit sampling clip {source}")
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
            Self::TripleWorkLimit { candidates, cap } => write!(
                f,
                "sheaf triple discovery refused {candidates} connected candidates; the deterministic cap is {cap}"
            ),
        }
    }
}

impl core::error::Error for SheafBuildError {}

/// Outward enclosure of a chart sample's signed distance from its own error
/// certificate. Only well-formed rigorous claims are usable: estimates,
/// no-claims, and malformed rigorous certificates poison the interface into
/// the whole extended real line, which can only yield `Unknown`.
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
    interfaces: &[Interface],
    cx: &Cx<'_>,
) -> Result<Vec<TripleCell>, SheafBuildError> {
    let mut adjacency: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    let mut edge_samples = BTreeMap::new();
    for (completed_edges, interface) in interfaces.iter().enumerate() {
        if completed_edges.is_multiple_of(256) {
            cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                stage: "triple-discovery",
                patches: None,
                completed_draws: completed_edges,
            })?;
        }
        let (a, b) = interface.patches;
        adjacency.entry(a).or_default().insert(b);
        adjacency.entry(b).or_default().insert(a);
        edge_samples.insert((a, b), interface.samples.len());
    }

    let mut triples = Vec::new();
    let mut candidates = 0usize;
    for (&(a, b), &ab_samples) in &edge_samples {
        let (Some(a_neighbors), Some(b_neighbors)) = (adjacency.get(&a), adjacency.get(&b)) else {
            continue;
        };
        for &c in a_neighbors.intersection(b_neighbors) {
            if c <= b {
                continue;
            }
            if candidates >= SHEAF_MAX_TRIPLE_CANDIDATES {
                return Err(SheafBuildError::TripleWorkLimit {
                    candidates: candidates.saturating_add(1),
                    cap: SHEAF_MAX_TRIPLE_CANDIDATES,
                });
            }
            if candidates.is_multiple_of(256) {
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "triple-discovery",
                    patches: None,
                    completed_draws: candidates,
                })?;
            }
            let Some(&bc_samples) = edge_samples.get(&(b, c)) else {
                continue;
            };
            let Some(&ac_samples) = edge_samples.get(&(a, c)) else {
                continue;
            };
            triples.push(TripleCell {
                patches: (a, b, c),
                samples: ab_samples.min(bc_samples).min(ac_samples),
            });
            candidates += 1;
        }
    }
    cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
        stage: "triple-discovery",
        patches: None,
        completed_draws: candidates,
    })?;
    Ok(triples)
}

impl SheafComplex {
    /// Build the complex from charts: interface discovery via
    /// support-overlap, shared-surface samples via the zero band of BOTH
    /// charts, triple junctions via triple overlaps.
    pub fn from_charts(
        charts: &[&dyn Chart],
        cx: &Cx<'_>,
    ) -> Result<SheafComplex, SheafBuildError> {
        Self::from_charts_with_clip(charts, None, cx)
    }

    /// Build interface evidence inside an explicit finite clip. Pairs outside
    /// the clip are skipped; invalid supports and unresolved unbounded shared
    /// domains are structured refusals.
    pub fn from_charts_clipped(
        charts: &[&dyn Chart],
        clip: Aabb,
        cx: &Cx<'_>,
    ) -> Result<SheafComplex, SheafBuildError> {
        Self::from_charts_with_clip(charts, Some(clip), cx)
    }

    #[allow(clippy::too_many_lines)] // One ordered discovery, sampling, cancellation, and finalize transaction.
    fn from_charts_with_clip(
        charts: &[&dyn Chart],
        clip: Option<Aabb>,
        cx: &Cx<'_>,
    ) -> Result<SheafComplex, SheafBuildError> {
        cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
            stage: "admission",
            patches: None,
            completed_draws: 0,
        })?;
        let n = charts.len();
        if let Some(explicit_clip) = clip {
            SamplingDomain::resolve(Aabb::WHOLE_SPACE, Some(explicit_clip))
                .map_err(|source| SheafBuildError::SamplingClip { source })?;
        }
        let mut interfaces = Vec::new();
        for u in 0..n {
            for v in (u + 1)..n {
                cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                    stage: "interface-discovery",
                    patches: Some((u, v)),
                    completed_draws: 0,
                })?;
                let support_u = charts[u].support();
                let support_v = charts[v].support();
                SamplingDomain::validate_support(support_u).map_err(|source| {
                    SheafBuildError::SamplingDomain {
                        patches: (u, v),
                        source,
                    }
                })?;
                SamplingDomain::validate_support(support_v).map_err(|source| {
                    SheafBuildError::SamplingDomain {
                        patches: (u, v),
                        source,
                    }
                })?;
                let Some(shared_support) = support_u.intersection(&support_v) else {
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
                let shared = domain.bounds();
                let spans = domain.spans();
                let diag = domain.diagonal();
                let band = BAND_FRACTION * diag;
                let mut state = box_seed(&shared);
                let mut samples = Vec::new();
                // Rejection-sample the shared zero band (bounded draws).
                for draw_index in 0..SAMPLES_PER_INTERFACE * 64 {
                    if samples.len() >= SAMPLES_PER_INTERFACE {
                        break;
                    }
                    let completed_draws = draw_index;
                    cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
                        stage: "interface-sampling",
                        patches: Some((u, v)),
                        completed_draws,
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
                        completed_draws,
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
                        completed_draws,
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
        }
        let triples = discover_triples(&interfaces, cx)?;
        cx.checkpoint().map_err(|_| SheafBuildError::Cancelled {
            stage: "finalize",
            patches: None,
            completed_draws: triples.len(),
        })?;
        Ok(SheafComplex {
            sampling_clip: clip,
            n_patches: n,
            interfaces,
            triples,
        })
    }

    /// Assemble δ⁰ (edges × patches) with ±1 entries: one row per
    /// interface SAMPLE, `+1` on patch v's slot, `−1` on patch u's.
    /// (Per-sample rows: the stalk of an edge is its sample space.)
    #[must_use]
    pub fn delta0(&self) -> Csr {
        let rows: usize = self.interfaces.iter().map(|i| i.samples.len()).sum();
        let mut coo = Coo::new(rows, self.n_patches);
        let mut r = 0usize;
        for iface in &self.interfaces {
            for _ in &iface.samples {
                coo.push(r, iface.patches.0, -1.0);
                coo.push(r, iface.patches.1, 1.0);
                r += 1;
            }
        }
        coo.assemble()
    }

    /// Assemble δ¹ (triples × edges) with ±1 entries per the oriented
    /// triangle boundary: for triple (a,b,c) with edges e_ab, e_bc, e_ac:
    /// `+e_ab + e_bc − e_ac` (edge-level stalks: one column per edge).
    #[must_use]
    pub fn delta1(&self) -> Csr {
        let edge_index = |a: usize, b: usize| {
            self.interfaces
                .iter()
                .position(|i| i.patches == (a.min(b), a.max(b)))
                .expect("triple implies edges")
        };
        let mut coo = Coo::new(self.triples.len(), self.interfaces.len());
        for (t, triple) in self.triples.iter().enumerate() {
            let (a, b, c) = triple.patches;
            coo.push(t, edge_index(a, b), 1.0);
            coo.push(t, edge_index(b, c), 1.0);
            coo.push(t, edge_index(a, c), -1.0);
        }
        coo.assemble()
    }

    /// Edge-level δ⁰ (edges × patches, one row per INTERFACE): the
    /// companion of [`Self::delta1`] for the bitwise δδ = 0 identity.
    #[must_use]
    pub fn delta0_edges(&self) -> Csr {
        let mut coo = Coo::new(self.interfaces.len(), self.n_patches);
        for (r, iface) in self.interfaces.iter().enumerate() {
            coo.push(r, iface.patches.0, -1.0);
            coo.push(r, iface.patches.1, 1.0);
        }
        coo.assemble()
    }

    /// Per-interface mismatch assessment. The VERDICT bits come from
    /// fs-ivl's sound predicates (`encloses`/`contains`). Reported magnitudes
    /// are aggregated directly from the intervals' outward endpoints; an
    /// indeterminate interval retains an infinite upper report and cannot
    /// authorize numerical evidence.
    #[must_use]
    pub fn mismatch_bounds(&self, tol: f64) -> Vec<InterfaceBound> {
        let valid_tolerance = tol.is_finite() && tol >= 0.0;
        let ok_band = valid_tolerance.then(|| Interval::new(0.0, tol));
        self.interfaces
            .iter()
            .map(|iface| {
                let mut all_within = valid_tolerance;
                let mut proven_leak = false;
                let mut determinate = valid_tolerance;
                let mut lo = 0.0f64;
                let mut hi = 0.0f64;
                for s in &iface.samples {
                    let d = (s.values[1] - s.values[0]).abs();
                    if !(d.lo().is_finite() && d.hi().is_finite()) {
                        determinate = false;
                        all_within = false;
                        hi = f64::INFINITY;
                        continue;
                    }
                    let within = ok_band.is_some_and(|band| band.encloses(d));
                    all_within &= within;
                    // |mismatch| enclosure entirely above tol: a proven
                    // violation (sound: the true value is inside d).
                    if valid_tolerance && d.lo() > tol {
                        proven_leak = true;
                    }
                    lo = lo.max(d.lo());
                    hi = hi.max(d.hi());
                }
                InterfaceBound {
                    patches: iface.patches,
                    all_within,
                    proven_leak,
                    determinate,
                    lo_report: lo,
                    hi_report: hi,
                }
            })
            .collect()
    }

    /// Least-squares section solve: per-patch gauge offsets minimizing
    /// the mean-square mismatch (graph-Laplacian normal equations, gauge
    /// fixed by pinning patch 0; deterministic Gauss–Seidel). Returns
    /// (offsets, raw ms mismatch, residual ms mismatch).
    #[must_use]
    pub fn section_solve(&self) -> (Vec<f64>, f64, f64) {
        let n = self.n_patches;
        let mut offsets = vec![0.0f64; n];
        // Edge means of the midpoint mismatch.
        let edges: Vec<((usize, usize), f64, usize)> = self
            .interfaces
            .iter()
            .map(|iface| {
                let mut sum = 0.0;
                for s in &iface.samples {
                    sum += s.values[1].midpoint() - s.values[0].midpoint();
                }
                (iface.patches, sum, s_len(iface))
            })
            .collect();
        let raw_ms = mean_square(&edges, &offsets);
        for _ in 0..200 {
            for p in 1..n {
                // Optimal c_p given the rest: weighted average balance.
                let mut num = 0.0f64;
                let mut den = 0.0f64;
                for ((u, v), sum, count) in &edges {
                    #[allow(clippy::cast_precision_loss)]
                    let w = *count as f64;
                    if *u == p {
                        num += sum + w * offsets[*v];
                        den += w;
                    } else if *v == p {
                        num += -sum + w * offsets[*u];
                        den += w;
                    }
                }
                if den > 0.0 {
                    offsets[p] = num / den;
                }
            }
        }
        let residual_ms = mean_square(&edges, &offsets);
        (offsets, raw_ms, residual_ms)
    }

    /// The watertightness certificate: interval-verified verdict as
    /// Evidence (enclosure numerics; content-addressed provenance).
    /// PASS requires every sample's enclosure INSIDE `[0, tol]` (sound);
    /// FAIL requires a proven-above-tolerance interface; anything else
    /// is an honest Unknown.
    #[must_use]
    pub fn watertightness(&self, tol: f64) -> Evidence<SheafVerdict> {
        let bounds = self.mismatch_bounds(tol);
        let worst_hi = bounds.iter().map(|b| b.hi_report).fold(0.0f64, f64::max);
        let worst_lo = bounds.iter().map(|b| b.lo_report).fold(0.0f64, f64::max);
        let all_determinate = !bounds.is_empty() && bounds.iter().all(|bound| bound.determinate);
        // PASS requires at least one DISCOVERED interface whose samples all lie
        // inside [0, tol]. An empty `bounds` means NO interface was found (charts
        // are disjoint/gapped/near-tangent, or the geometry is empty) — the
        // interface-agreement check gathered NO evidence, so the honest verdict
        // is `Unknown`, not a positive `worst_mismatch = 0` PASS. `all()` on an
        // empty set is vacuously true; guarding on non-emptiness closes the
        // vacuous-truth false certificate (bead obnw).
        let all_pass = all_determinate && bounds.iter().all(|b| b.all_within);
        let obstruction: Vec<((usize, usize), f64)> = bounds
            .iter()
            .filter(|b| b.proven_leak)
            .map(|b| (b.patches, b.lo_report))
            .collect();
        let verdict = if all_pass {
            SheafVerdict::Pass {
                worst_mismatch: worst_hi,
                margins: bounds.iter().map(|b| (b.patches, b.hi_report)).collect(),
            }
        } else if all_determinate && !obstruction.is_empty() {
            let (_, raw, residual) = self.section_solve();
            let share = if raw.is_finite() && residual.is_finite() && raw > 0.0 && residual >= 0.0 {
                Some((1.0 - residual / raw).clamp(0.0, 1.0))
            } else if raw == 0.0 && residual == 0.0 {
                Some(0.0)
            } else {
                None
            };
            SheafVerdict::Fail {
                obstruction,
                coboundary_share: share,
            }
        } else {
            SheafVerdict::Unknown {
                straddling: bounds
                    .iter()
                    .filter(|b| !b.determinate || !b.all_within)
                    .map(|b| (b.patches, b.lo_report, b.hi_report))
                    .collect(),
            }
        };
        let mut canon = format!(
            "sheaf-watertightness;patches={};interfaces={};tol={tol}",
            self.n_patches,
            self.interfaces.len()
        );
        match self.sampling_clip {
            None => canon.push_str(";sampling_clip=none"),
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
        for b in &bounds {
            let _ = write!(
                canon,
                ";{}-{}:{}:{}",
                b.patches.0, b.patches.1, b.lo_report, b.hi_report
            );
        }
        let numerical = if all_determinate
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
            provenance: ProvenanceHash::of_bytes(canon.as_bytes()),
            adjoint_ref: None,
            value: verdict,
        }
    }
}

fn s_len(iface: &Interface) -> usize {
    iface.samples.len()
}

fn mean_square(edges: &[((usize, usize), f64, usize)], offsets: &[f64]) -> f64 {
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for ((u, v), sum, count) in edges {
        #[allow(clippy::cast_precision_loss)]
        let w = *count as f64;
        // Mean mismatch on the edge after gauging.
        let gauged = sum / w + offsets[*v] - offsets[*u];
        num += w * gauged * gauged;
        den += w;
    }
    if den > 0.0 { num / den } else { 0.0 }
}

/// THE INDEPENDENT FALSIFIER (registry pairing: watertightness →
/// ray-parity): a different algorithm on a different code path. March a
/// segment through the union-of-charts model counting sign changes of
/// the min-SDF; a closed (watertight) model yields an EVEN count on
/// segments with both endpoints strictly outside. Returns the violating ray
/// index, if any. The sample count is work-capped, every endpoint and field
/// sample is validated, convex interpolation avoids an overflowing endpoint
/// subtraction, and cancellation is observed after every chart evaluation.
///
/// # Errors
/// [`RayParityError`] when the inputs do not satisfy the finite outside-ray
/// preconditions, the public work cap would be exceeded, a chart produces a
/// non-finite value, or cancellation is observed.
pub fn ray_parity_falsifier(
    charts: &[&dyn Chart],
    rays: &[(Point3, Point3)],
    steps: usize,
    cx: &Cx<'_>,
) -> Result<Option<usize>, RayParityError> {
    if charts.is_empty() {
        return Err(RayParityError::EmptyCharts);
    }
    if rays.is_empty() {
        return Err(RayParityError::EmptyRays);
    }
    if steps == 0 {
        return Err(RayParityError::InvalidSteps { steps });
    }
    let requested = (rays.len() as u128)
        .saturating_mul((steps as u128).saturating_add(1))
        .saturating_mul(charts.len() as u128);
    if requested > RAY_PARITY_MAX_EVALUATIONS as u128 {
        return Err(RayParityError::WorkLimitExceeded {
            requested,
            cap: RAY_PARITY_MAX_EVALUATIONS,
        });
    }

    for (ray, (start, end)) in rays.iter().copied().enumerate() {
        for (endpoint, point) in [(RayEndpoint::Start, start), (RayEndpoint::End, end)] {
            if !finite_point(point) {
                return Err(RayParityError::NonFiniteEndpoint {
                    ray,
                    endpoint,
                    point,
                });
            }
        }
    }

    let mut completed_points = 0usize;
    let mut completed_chart_evaluations = 0usize;
    for (ri, (a, b)) in rays.iter().enumerate() {
        let mut crossings = 0usize;
        let mut prev_sign = None;
        for k in 0..=steps {
            checkpoint_ray_parity(cx, ri, completed_points, completed_chart_evaluations)?;
            let p = ray_sample_point(*a, *b, k, steps);
            if !finite_point(p) {
                return Err(RayParityError::NonRepresentableSamplePoint {
                    ray: ri,
                    step: k,
                    point: p,
                });
            }
            let mut d = f64::INFINITY;
            for (chart_index, chart) in charts.iter().enumerate() {
                let value = chart.eval(p, cx).signed_distance;
                completed_chart_evaluations = completed_chart_evaluations.saturating_add(1);
                checkpoint_ray_parity(cx, ri, completed_points, completed_chart_evaluations)?;
                if !value.is_finite() {
                    return Err(RayParityError::NonFiniteSample {
                        ray: ri,
                        step: k,
                        chart: chart_index,
                        value,
                    });
                }
                d = d.min(value);
            }
            if k == 0 && d <= 0.0 {
                return Err(RayParityError::EndpointNotOutside {
                    ray: ri,
                    endpoint: RayEndpoint::Start,
                    min_signed_distance: d,
                });
            }
            if k == steps && d <= 0.0 {
                return Err(RayParityError::EndpointNotOutside {
                    ray: ri,
                    endpoint: RayEndpoint::End,
                    min_signed_distance: d,
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
        if !crossings.is_multiple_of(2) {
            return Ok(Some(ri));
        }
    }
    checkpoint_ray_parity(
        cx,
        rays.len(),
        completed_points,
        completed_chart_evaluations,
    )?;
    Ok(None)
}

fn checkpoint_ray_parity(
    cx: &Cx<'_>,
    completed_rays: usize,
    completed_points: usize,
    completed_chart_evaluations: usize,
) -> Result<(), RayParityError> {
    cx.checkpoint().map_err(|_| RayParityError::Cancelled {
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

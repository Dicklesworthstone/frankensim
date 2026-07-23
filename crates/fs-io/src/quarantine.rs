//! THE QUARANTINE PRINCIPLE: imports land as [`Quarantined`] (raw value,
//! source receipt, detected defects) and promote to `Evidence` ONLY
//! after the repair suite runs and validity checks pass. Unrepaired
//! defects BLOCK promotion with actionable diagnostics. The receipt is
//! the ledger `imports` row payload (HELM writes it; L2 emits it).

use crate::IoError;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate,
};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_rep_mesh::{HalfEdgeMesh, Soup, repair};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Semantic identity of the tolerance-aware import census and promotion receipt.
pub const IMPORT_CENSUS_SEMANTICS_VERSION: &str = "fs-io-import-census-v1";

/// A defect found at import time (pre-repair census).
#[derive(Debug, Clone, PartialEq)]
pub struct ImportDefect {
    /// Defect class ("degenerate-face", "duplicate-face",
    /// "unreferenced-vertex", "non-manifold-or-open").
    pub class: &'static str,
    /// How many instances.
    pub count: usize,
}

/// How the quadratic shell-overlap/self-intersection check spends its pair budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectionInspection {
    /// Visit triangle pairs in lexicographic order until all pairs or the budget is exhausted.
    ExhaustiveF64 {
        /// Maximum raw triangle-pair positions to visit.
        max_pair_tests: usize,
    },
    /// Visit an evenly spaced deterministic sample of the raw triangle-pair index space.
    DeterministicSampleF64 {
        /// Maximum raw triangle-pair positions to sample.
        sample_count: usize,
    },
}

impl IntersectionInspection {
    const fn budget(self) -> usize {
        match self {
            Self::ExhaustiveF64 { max_pair_tests } => max_pair_tests,
            Self::DeterministicSampleF64 { sample_count } => sample_count,
        }
    }

    const fn requested_level(self) -> &'static str {
        match self {
            Self::ExhaustiveF64 { .. } => "exhaustive-f64-filter",
            Self::DeterministicSampleF64 { .. } => "deterministic-even-sample-f64-filter",
        }
    }
}

/// Explicit tolerance, pair budget, and cancellation envelope for an extended census.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImportCensusPolicy {
    model_tolerance: f64,
    intersection: IntersectionInspection,
    poll_stride: usize,
}

impl ImportCensusPolicy {
    /// Admit a census policy.
    ///
    /// # Errors
    /// Returns a structured refusal for a non-positive/non-finite tolerance, a zero
    /// intersection budget, or a zero cancellation-poll stride.
    pub fn try_new(
        model_tolerance: f64,
        intersection: IntersectionInspection,
        poll_stride: usize,
    ) -> Result<Self, CensusRefusal> {
        if !model_tolerance.is_finite() || model_tolerance <= 0.0 {
            return Err(CensusRefusal::invalid(
                "model_tolerance",
                "must be finite and strictly positive",
            ));
        }
        if intersection.budget() == 0 {
            return Err(CensusRefusal::invalid(
                "intersection budget",
                "must admit at least one pair test or sample",
            ));
        }
        if poll_stride == 0 {
            return Err(CensusRefusal::invalid(
                "poll_stride",
                "must be at least one",
            ));
        }
        Ok(Self {
            model_tolerance,
            intersection,
            poll_stride,
        })
    }

    /// Declared model tolerance in the caller's geometry length unit.
    #[must_use]
    pub const fn model_tolerance(self) -> f64 {
        self.model_tolerance
    }

    /// Requested intersection-inspection strategy.
    #[must_use]
    pub const fn intersection(self) -> IntersectionInspection {
        self.intersection
    }

    /// Maximum loop iterations between explicit `Cx` cancellation polls.
    #[must_use]
    pub const fn poll_stride(self) -> usize {
        self.poll_stride
    }
}

/// A structured census admission or cancellation refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CensusRefusal {
    /// Stable stage or invalid policy field.
    pub stage: &'static str,
    /// Deterministic position within the stage, when applicable.
    pub at: usize,
    /// Actionable explanation.
    pub reason: String,
    /// True only when `Cx` cancellation was observed.
    pub cancelled: bool,
}

impl CensusRefusal {
    fn invalid(stage: &'static str, reason: &'static str) -> Self {
        Self {
            stage,
            at: 0,
            reason: reason.to_string(),
            cancelled: false,
        }
    }

    fn cancelled(stage: &'static str, at: usize) -> Self {
        Self {
            stage,
            at,
            reason: "cancellation requested before census publication".to_string(),
            cancelled: true,
        }
    }
}

impl core::fmt::Display for CensusRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "import census refused at {} position {}: {}",
            self.stage, self.at, self.reason
        )
    }
}

impl core::error::Error for CensusRefusal {}

/// What the overlap/self-intersection phase actually inspected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntersectionCensusEvidence {
    /// Requested inspection level.
    pub requested_level: &'static str,
    /// Pair-test/sample budget from policy.
    pub pair_budget: usize,
    /// Raw triangle-pair positions visited or sampled.
    pub visited_raw_pairs: usize,
    /// Non-adjacent triangle pairs actually tested.
    pub inspected_pairs: usize,
    /// Raw sampled pairs skipped because they shared an indexed vertex.
    pub shared_vertex_pairs_skipped: usize,
    /// Intersecting non-adjacent pairs found.
    pub intersecting_pairs: usize,
    /// True only when every raw triangle pair was visited.
    pub complete: bool,
    /// Explicit numerical authority boundary.
    pub authority: &'static str,
}

/// Diagnostic input for E08 geometry budgeting.
///
/// This is deliberately not called an error bound: small features and a sampled
/// intersection census do not prove a Hausdorff enclosure.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportGeometryBudgetInput {
    /// Declared model tolerance in the caller's geometry length unit.
    pub model_tolerance: f64,
    /// Residual small-edge, sliver, and near-loop-gap finding count.
    pub tolerance_sensitive_residuals: usize,
    /// Largest detected near-loop gap, if any.
    pub largest_detected_boundary_gap: Option<f64>,
    /// Whether the intersection inspection covered every raw triangle pair.
    pub intersection_inspection_complete: bool,
    /// Explicit no-claim boundary for downstream budget composition.
    pub authority: &'static str,
}

/// Complete tolerance-aware diagnostic census.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportCensusReport {
    /// Stable semantic version.
    pub semantics_version: &'static str,
    /// Counted defect classes in deterministic class order.
    pub findings: Vec<ImportDefect>,
    /// Number of simple closed boundary-loop components.
    pub boundary_loops: usize,
    /// Smallest valid indexed edge length observed.
    pub smallest_edge: Option<f64>,
    /// Smallest non-degenerate triangle altitude observed.
    pub smallest_altitude: Option<f64>,
    /// Largest detected gap no greater than the declared tolerance.
    pub largest_detected_boundary_gap: Option<f64>,
    /// Pair-inspection evidence and authority.
    pub intersection: IntersectionCensusEvidence,
    /// Explicit downstream geometry-budget input.
    pub geometry_budget: ImportGeometryBudgetInput,
    /// Explicit cancellation poll stride exercised by the census.
    pub cancellation_poll_stride: usize,
}

impl ImportCensusReport {
    /// Count one class, returning zero when it is absent.
    #[must_use]
    pub fn count(&self, class: &str) -> usize {
        self.findings
            .iter()
            .find(|finding| finding.class == class)
            .map_or(0, |finding| finding.count)
    }

    /// Canonical JSON suitable for a supplier-corpus row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut output = format!(
            "{{\"semantics\":\"{}\",\"model_tolerance\":{},\
             \"model_tolerance_bits\":\"{:016x}\",\"findings\":[",
            self.semantics_version,
            canonical_f64(self.geometry_budget.model_tolerance),
            self.geometry_budget.model_tolerance.to_bits()
        );
        append_findings_json(&mut output, &self.findings);
        let _ = write!(
            output,
            "],\"boundary_loops\":{},\"smallest_edge\":{},\"smallest_altitude\":{},\
             \"largest_detected_boundary_gap\":{},\"intersection\":{{\
             \"requested_level\":\"{}\",\"pair_budget\":{},\"visited_raw_pairs\":{},\
             \"inspected_pairs\":{},\
             \"shared_vertex_pairs_skipped\":{},\"intersecting_pairs\":{},\"complete\":{},\
             \"authority\":\"{}\"}},\"cancellation\":{{\"kind\":\"cx-polled\",\
             \"poll_stride\":{}}},\"geometry_budget\":{{\"model_tolerance\":{},\
             \"tolerance_sensitive_residuals\":{},\"largest_detected_boundary_gap\":{},\
             \"intersection_inspection_complete\":{},\"authority\":\"{}\"}}}}",
            self.boundary_loops,
            optional_f64_json(self.smallest_edge),
            optional_f64_json(self.smallest_altitude),
            optional_f64_json(self.largest_detected_boundary_gap),
            self.intersection.requested_level,
            self.intersection.pair_budget,
            self.intersection.visited_raw_pairs,
            self.intersection.inspected_pairs,
            self.intersection.shared_vertex_pairs_skipped,
            self.intersection.intersecting_pairs,
            self.intersection.complete,
            self.intersection.authority,
            self.cancellation_poll_stride,
            canonical_f64(self.geometry_budget.model_tolerance),
            self.geometry_budget.tolerance_sensitive_residuals,
            optional_f64_json(self.geometry_budget.largest_detected_boundary_gap),
            self.geometry_budget.intersection_inspection_complete,
            self.geometry_budget.authority
        );
        output
    }
}

/// Maximum accepted residual count per defect class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImportRefusalThresholds {
    /// Degenerate faces accepted after repair.
    pub degenerate_faces: usize,
    /// Duplicate faces accepted after repair.
    pub duplicate_faces: usize,
    /// Unreferenced vertices accepted after repair.
    pub unreferenced_vertices: usize,
    /// Invalid face indices accepted after repair.
    pub invalid_face_indices: usize,
    /// Non-finite vertices accepted after repair.
    pub non_finite_vertices: usize,
    /// Non-manifold/open edges accepted after repair.
    pub non_manifold_or_open: usize,
    /// Tolerance-relative small edges accepted after repair.
    pub small_edges: usize,
    /// Tolerance-relative sliver faces accepted after repair.
    pub sliver_faces: usize,
    /// Near-adjacent boundary-loop gaps accepted after repair.
    pub boundary_loop_gaps: usize,
    /// Non-adjacent triangle intersections accepted after repair.
    pub shell_intersections: usize,
    /// Require every raw triangle pair to be visited.
    pub require_complete_intersection_census: bool,
}

impl ImportRefusalThresholds {
    /// Validation-grade defaults: no residual defect is accepted except
    /// unreferenced vertices, matching the legacy cosmetic exception.
    #[must_use]
    pub const fn validation_grade() -> Self {
        Self {
            degenerate_faces: 0,
            duplicate_faces: 0,
            unreferenced_vertices: usize::MAX,
            invalid_face_indices: 0,
            non_finite_vertices: 0,
            non_manifold_or_open: 0,
            small_edges: 0,
            sliver_faces: 0,
            boundary_loop_gaps: 0,
            shell_intersections: 0,
            require_complete_intersection_census: true,
        }
    }

    fn maximum_for(self, class: &str) -> usize {
        match class {
            "degenerate-face" => self.degenerate_faces,
            "duplicate-face" => self.duplicate_faces,
            "unreferenced-vertex" => self.unreferenced_vertices,
            "invalid-face-index" => self.invalid_face_indices,
            "non-finite-vertex" => self.non_finite_vertices,
            "non-manifold-or-open" => self.non_manifold_or_open,
            "small-edge" => self.small_edges,
            "sliver-face" => self.sliver_faces,
            "near-boundary-loop-gap" => self.boundary_loop_gaps,
            "shell-overlap-or-self-intersection" => self.shell_intersections,
            _ => 0,
        }
    }

    fn to_json(self) -> String {
        format!(
            "{{\"degenerate-face\":{},\"duplicate-face\":{},\"unreferenced-vertex\":{},\
             \"invalid-face-index\":{},\"non-finite-vertex\":{},\
             \"non-manifold-or-open\":{},\"small-edge\":{},\"sliver-face\":{},\
             \"near-boundary-loop-gap\":{},\"shell-overlap-or-self-intersection\":{},\
             \"require_complete_intersection_census\":{}}}",
            self.degenerate_faces,
            self.duplicate_faces,
            self.unreferenced_vertices,
            self.invalid_face_indices,
            self.non_finite_vertices,
            self.non_manifold_or_open,
            self.small_edges,
            self.sliver_faces,
            self.boundary_loop_gaps,
            self.shell_intersections,
            self.require_complete_intersection_census
        )
    }
}

/// Project-selected healing and refusal policy. Its profile and every threshold
/// are embedded in the promotion receipt.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImportPromotionPolicy {
    /// Stable project policy/profile identifier.
    pub profile: &'static str,
    /// Maximum boundary-loop size passed to the existing fan-fill repair.
    pub max_hole_edges: usize,
    /// Extended census envelope.
    pub census: ImportCensusPolicy,
    /// Residual thresholds.
    pub thresholds: ImportRefusalThresholds,
}

/// A typed promotion receipt retaining pre-repair, repair, and residual evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportPromotionReceipt {
    /// Source receipt.
    pub source: ImportReceipt,
    /// Project-selected policy.
    pub policy: ImportPromotionPolicy,
    /// Pre-repair census.
    pub before: ImportCensusReport,
    /// Structured repair operations from `fs-rep-mesh`.
    pub repairs: Vec<fs_rep_mesh::RepairReceipt>,
    /// Post-repair residual census.
    pub after: ImportCensusReport,
    /// `promoted` or `refused`.
    pub trust: &'static str,
}

impl ImportPromotionReceipt {
    /// Canonical ledger/corpus JSON, including class deltas and no-claim boundaries.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut output = format!(
            "{{\"kind\":\"import-promotion-receipt\",\"semantics\":\"{}\",\
             \"source\":{{\"format\":\"{}\",\"source_hash\":\"{:016x}\",\
             \"parser\":\"{}\",\"vertices\":{},\"triangles\":{}}},\
             \"policy\":{{\"profile\":\"{}\",\"max_hole_edges\":{},\
             \"model_tolerance\":{},\"model_tolerance_bits\":\"{:016x}\",\
             \"intersection_level\":\"{}\",\"intersection_budget\":{},\
             \"poll_stride\":{},\"thresholds\":{}}},\"before\":{},\
             \"repair_history\":{{\"operations\":[",
            IMPORT_CENSUS_SEMANTICS_VERSION,
            json_escape(self.source.format),
            self.source.source_hash,
            json_escape(self.source.parser_version),
            self.source.parsed.0,
            self.source.parsed.1,
            json_escape(self.policy.profile),
            self.policy.max_hole_edges,
            canonical_f64(self.policy.census.model_tolerance),
            self.policy.census.model_tolerance.to_bits(),
            self.policy.census.intersection.requested_level(),
            self.policy.census.intersection.budget(),
            self.policy.census.poll_stride,
            self.policy.thresholds.to_json(),
            self.before.to_json()
        );
        for (index, repair) in self.repairs.iter().enumerate() {
            if index > 0 {
                output.push(',');
            }
            output.push_str(&repair.to_json());
        }
        output.push_str("],\"class_deltas\":[");
        let classes = finding_classes(&self.before.findings, &self.after.findings);
        for (index, class) in classes.iter().enumerate() {
            if index > 0 {
                output.push(',');
            }
            let _ = write!(
                output,
                "{{\"class\":\"{}\",\"before\":{},\"after\":{}}}",
                class,
                self.before.count(class),
                self.after.count(class)
            );
        }
        let _ = write!(
            output,
            "],\"tolerance_consumed\":0.0,\"tolerance_based_repair_performed\":false}},\
             \"after\":{},\"geometry_budget\":{{\"model_tolerance\":{},\
             \"tolerance_sensitive_residuals\":{},\"largest_detected_boundary_gap\":{},\
             \"intersection_inspection_complete\":{},\"authority\":\"{}\"}},\
             \"trust\":\"{}\"}}",
            self.after.to_json(),
            canonical_f64(self.after.geometry_budget.model_tolerance),
            self.after.geometry_budget.tolerance_sensitive_residuals,
            optional_f64_json(self.after.geometry_budget.largest_detected_boundary_gap),
            self.after.geometry_budget.intersection_inspection_complete,
            self.after.geometry_budget.authority,
            self.trust
        );
        output
    }
}

/// Failure before or after repair under an explicit promotion policy.
#[derive(Debug)]
pub enum ImportPromotionError {
    /// Policy admission or `Cx` cancellation.
    Census(CensusRefusal),
    /// Completed census exceeded the receipted project thresholds.
    Refused(Box<PromotionRefusal>),
}

impl core::fmt::Display for ImportPromotionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Census(error) => error.fmt(f),
            Self::Refused(error) => {
                write!(f, "import promotion refused: {}", error.blocking.join(", "))
            }
        }
    }
}

impl core::error::Error for ImportPromotionError {}

/// The certification receipt: full provenance for "where did this
/// geometry come from and what did we fix."
#[derive(Debug, Clone, PartialEq)]
pub struct ImportReceipt {
    /// Declared source format ("stl", "obj", "ply", …).
    pub format: &'static str,
    /// Content hash of the raw source bytes (FNV-1a; BLAKE3-class hash
    /// upgrades HELM-side with the same field).
    pub source_hash: u64,
    /// Parser version (this crate's version).
    pub parser_version: &'static str,
    /// Element counts as parsed (vertices, triangles).
    pub parsed: (usize, usize),
}

impl ImportReceipt {
    /// Canonical JSON (the ledger `imports` row payload).
    #[must_use]
    pub fn to_json(&self, defects: &[ImportDefect], trust: &str) -> String {
        let mut s = format!(
            "{{\"kind\":\"import-receipt\",\"format\":\"{}\",\"source_hash\":\"{:016x}\",\
             \"parser\":\"{}\",\"vertices\":{},\"triangles\":{},\"defects\":[",
            self.format, self.source_hash, self.parser_version, self.parsed.0, self.parsed.1
        );
        for (i, d) in defects.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{{\"class\":\"{}\",\"count\":{}}}", d.class, d.count);
        }
        let _ = write!(s, "],\"trust\":\"{trust}\"}}");
        s
    }
}

/// An imported value in quarantine: NOT trusted, carries its census.
#[derive(Debug)]
pub struct Quarantined<T> {
    /// The raw parsed value.
    pub raw: T,
    /// Source provenance.
    pub source_receipt: ImportReceipt,
    /// Defects detected at import.
    pub defects: Vec<ImportDefect>,
}

/// A structured promotion refusal: what blocked, and what to do.
#[derive(Debug, Clone, PartialEq)]
pub struct PromotionRefusal {
    /// The blocking defect classes after repair.
    pub blocking: Vec<String>,
    /// Actionable guidance.
    pub fixes: Vec<String>,
    /// The post-repair receipt JSON (for the ledger row, trust=refused).
    pub receipt_json: String,
}

/// Census a soup for import defects (pre-repair, read-only).
#[must_use]
pub fn census(soup: &Soup) -> Vec<ImportDefect> {
    let mut defects = Vec::new();
    let mut degenerate = 0usize;
    let mut seen = std::collections::BTreeSet::new();
    let mut duplicate = 0usize;
    for t in &soup.triangles {
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            degenerate += 1;
        }
        let mut key = *t;
        key.sort_unstable();
        if !seen.insert(key) {
            duplicate += 1;
        }
    }
    if degenerate > 0 {
        defects.push(ImportDefect {
            class: "degenerate-face",
            count: degenerate,
        });
    }
    if duplicate > 0 {
        defects.push(ImportDefect {
            class: "duplicate-face",
            count: duplicate,
        });
    }
    let mut referenced = vec![false; soup.positions.len()];
    for t in &soup.triangles {
        for &i in t {
            if let Some(slot) = referenced.get_mut(i as usize) {
                *slot = true;
            }
        }
    }
    let unreferenced = referenced.iter().filter(|&&r| !r).count();
    if unreferenced > 0 {
        defects.push(ImportDefect {
            class: "unreferenced-vertex",
            count: unreferenced,
        });
    }
    // Closedness + edge-manifoldness by direct edge counting: every
    // undirected edge of a watertight 2-manifold appears exactly twice.
    // (The half-edge builder alone is insufficient: it legally accepts
    // open boundaries.)
    let mut edge_counts: std::collections::BTreeMap<(u32, u32), u32> =
        std::collections::BTreeMap::new();
    for t in &soup.triangles {
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            continue; // degenerates counted above
        }
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = (a.min(b), a.max(b));
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }
    let bad_edges = edge_counts.values().filter(|&&c| c != 2).count();
    let vertex_nonmanifold =
        HalfEdgeMesh::from_triangles(soup.positions.clone(), &soup.triangles).is_err();
    if bad_edges > 0 || vertex_nonmanifold {
        defects.push(ImportDefect {
            class: "non-manifold-or-open",
            count: bad_edges.max(1),
        });
    }
    defects
}

/// Run the tolerance-aware census under an explicit cost and cancellation policy.
///
/// Small edges are unique indexed edges whose length is no greater than the
/// declared model tolerance. A sliver is a non-degenerate triangle whose
/// altitude to its longest edge is no greater than that tolerance. Boundary
/// gaps are counted once per pair of simple boundary loops whose closest indexed
/// vertices are within tolerance.
///
/// The intersection pass tests only triangle pairs that do not share an indexed
/// vertex. It is an f64 filter, not an exact-predicate certificate; the receipt
/// states whether the raw pair space was exhaustive or sampled.
///
/// # Errors
/// [`CensusRefusal`] when cancellation is observed.
pub fn census_with_policy(
    soup: &Soup,
    policy: ImportCensusPolicy,
    cx: &Cx<'_>,
) -> Result<ImportCensusReport, CensusRefusal> {
    checkpoint(cx, "census-start", 0)?;

    let mut findings = census(soup);
    checkpoint(cx, "basic-topology-census", soup.triangles.len())?;

    let non_finite_vertices = soup
        .positions
        .iter()
        .filter(|point| !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite())
        .count();
    let invalid_face_indices = soup
        .triangles
        .iter()
        .filter(|triangle| {
            triangle
                .iter()
                .any(|&vertex| vertex as usize >= soup.positions.len())
        })
        .count();
    push_finding(&mut findings, "invalid-face-index", invalid_face_indices);
    push_finding(&mut findings, "non-finite-vertex", non_finite_vertices);

    let mut unique_edges = BTreeSet::new();
    let mut smallest_edge: Option<f64> = None;
    let mut smallest_altitude: Option<f64> = None;
    let mut small_edges = 0usize;
    let mut sliver_faces = 0usize;
    for (face_index, triangle) in soup.triangles.iter().enumerate() {
        poll(
            cx,
            policy.poll_stride,
            "tolerance-relative-feature-census",
            face_index,
        )?;
        if triangle
            .iter()
            .any(|&vertex| vertex as usize >= soup.positions.len())
        {
            continue;
        }
        let [a, b, c] = triangle.map(|vertex| soup.positions[vertex as usize]);
        if !points_are_finite([a, b, c]) {
            continue;
        }
        let lengths = [distance(a, b), distance(b, c), distance(c, a)];
        for ((u, v), length) in [
            ((triangle[0], triangle[1]), lengths[0]),
            ((triangle[1], triangle[2]), lengths[1]),
            ((triangle[2], triangle[0]), lengths[2]),
        ] {
            let key = (u.min(v), u.max(v));
            if unique_edges.insert(key) {
                smallest_edge = Some(smallest_edge.map_or(length, |old| old.min(length)));
                if length <= policy.model_tolerance {
                    small_edges += 1;
                }
            }
        }
        let longest = lengths.into_iter().fold(0.0_f64, f64::max);
        let twice_area = norm(cross(subtract(b, a), subtract(c, a)));
        if longest > 0.0 && twice_area > 0.0 {
            let altitude = twice_area / longest;
            smallest_altitude = Some(smallest_altitude.map_or(altitude, |old| old.min(altitude)));
            if altitude <= policy.model_tolerance {
                sliver_faces += 1;
            }
        }
    }
    push_finding(&mut findings, "small-edge", small_edges);
    push_finding(&mut findings, "sliver-face", sliver_faces);

    let boundary_loops = simple_boundary_loops(soup, policy.poll_stride, cx)?;
    let (boundary_loop_gaps, largest_detected_boundary_gap) =
        boundary_gap_census(soup, &boundary_loops, policy, cx)?;
    push_finding(&mut findings, "near-boundary-loop-gap", boundary_loop_gaps);

    let intersection = intersection_census(soup, policy, cx)?;
    push_finding(
        &mut findings,
        "shell-overlap-or-self-intersection",
        intersection.intersecting_pairs,
    );
    findings.sort_by_key(|finding| finding_order(finding.class));

    let tolerance_sensitive_residuals = small_edges + sliver_faces + boundary_loop_gaps;
    let geometry_budget = ImportGeometryBudgetInput {
        model_tolerance: policy.model_tolerance,
        tolerance_sensitive_residuals,
        largest_detected_boundary_gap,
        intersection_inspection_complete: intersection.complete,
        authority: "diagnostic-input-not-a-spatial-error-bound",
    };
    checkpoint(cx, "census-publication", soup.triangles.len())?;
    Ok(ImportCensusReport {
        semantics_version: IMPORT_CENSUS_SEMANTICS_VERSION,
        findings,
        boundary_loops: boundary_loops.len(),
        smallest_edge,
        smallest_altitude,
        largest_detected_boundary_gap,
        intersection,
        geometry_budget,
        cancellation_poll_stride: policy.poll_stride,
    })
}

fn push_finding(findings: &mut Vec<ImportDefect>, class: &'static str, count: usize) {
    if count == 0 {
        return;
    }
    if let Some(existing) = findings.iter_mut().find(|finding| finding.class == class) {
        existing.count = count;
    } else {
        findings.push(ImportDefect { class, count });
    }
}

fn finding_order(class: &str) -> (usize, &str) {
    let order = match class {
        "invalid-face-index" => 0,
        "non-finite-vertex" => 1,
        "degenerate-face" => 2,
        "duplicate-face" => 3,
        "unreferenced-vertex" => 4,
        "non-manifold-or-open" => 5,
        "small-edge" => 6,
        "sliver-face" => 7,
        "near-boundary-loop-gap" => 8,
        "shell-overlap-or-self-intersection" => 9,
        _ => usize::MAX,
    };
    (order, class)
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), CensusRefusal> {
    cx.checkpoint()
        .map_err(|_| CensusRefusal::cancelled(stage, at))
}

fn poll(
    cx: &Cx<'_>,
    poll_stride: usize,
    stage: &'static str,
    at: usize,
) -> Result<(), CensusRefusal> {
    if at.is_multiple_of(poll_stride) {
        checkpoint(cx, stage, at)?;
    }
    Ok(())
}

fn points_are_finite(points: [Point3; 3]) -> bool {
    points
        .iter()
        .all(|point| point.x.is_finite() && point.y.is_finite() && point.z.is_finite())
}

fn subtract(a: Point3, b: Point3) -> [f64; 3] {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -(a[2] * b[1])),
        a[2].mul_add(b[0], -(a[0] * b[2])),
        a[0].mul_add(b[1], -(a[1] * b[0])),
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn distance(a: Point3, b: Point3) -> f64 {
    norm(subtract(a, b))
}

fn simple_boundary_loops(
    soup: &Soup,
    poll_stride: usize,
    cx: &Cx<'_>,
) -> Result<Vec<Vec<u32>>, CensusRefusal> {
    let mut edge_counts: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for (face_index, triangle) in soup.triangles.iter().enumerate() {
        poll(cx, poll_stride, "boundary-loop-edge-census", face_index)?;
        if triangle
            .iter()
            .any(|&vertex| vertex as usize >= soup.positions.len())
        {
            continue;
        }
        for (a, b) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            *edge_counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut adjacency: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (edge_index, (&(a, b), &count)) in edge_counts.iter().enumerate() {
        poll(cx, poll_stride, "boundary-loop-adjacency", edge_index)?;
        if count == 1 && a != b {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
    }
    for (vertex_index, neighbors) in adjacency.values_mut().enumerate() {
        poll(
            cx,
            poll_stride,
            "boundary-loop-neighbor-order",
            vertex_index,
        )?;
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    let mut visited = BTreeSet::new();
    let mut loops = Vec::new();
    let mut traversal_visits = 0usize;
    for (component_index, &start) in adjacency.keys().enumerate() {
        poll(cx, poll_stride, "boundary-loop-component", component_index)?;
        if visited.contains(&start) {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        visited.insert(start);
        while let Some(vertex) = stack.pop() {
            poll(cx, poll_stride, "boundary-loop-traversal", traversal_visits)?;
            traversal_visits = traversal_visits.saturating_add(1);
            component.push(vertex);
            if let Some(neighbors) = adjacency.get(&vertex) {
                for &neighbor in neighbors {
                    if visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
        component.sort_unstable();
        if component.len() >= 3
            && component.iter().all(|vertex| {
                adjacency
                    .get(vertex)
                    .is_some_and(|neighbors| neighbors.len() == 2)
            })
        {
            loops.push(component);
        }
    }
    loops.sort();
    Ok(loops)
}

fn boundary_gap_census(
    soup: &Soup,
    loops: &[Vec<u32>],
    policy: ImportCensusPolicy,
    cx: &Cx<'_>,
) -> Result<(usize, Option<f64>), CensusRefusal> {
    let mut gaps = 0usize;
    let mut largest_gap: Option<f64> = None;
    let mut visits = 0usize;
    for left in 0..loops.len() {
        for right in (left + 1)..loops.len() {
            let mut closest = f64::INFINITY;
            for &a in &loops[left] {
                for &b in &loops[right] {
                    poll(cx, policy.poll_stride, "boundary-loop-gap-census", visits)?;
                    visits = visits.saturating_add(1);
                    closest = closest.min(distance(
                        soup.positions[a as usize],
                        soup.positions[b as usize],
                    ));
                }
            }
            if closest.is_finite() && closest <= policy.model_tolerance {
                gaps += 1;
                largest_gap = Some(largest_gap.map_or(closest, |old| old.max(closest)));
            }
        }
    }
    Ok((gaps, largest_gap))
}

fn intersection_census(
    soup: &Soup,
    policy: ImportCensusPolicy,
    cx: &Cx<'_>,
) -> Result<IntersectionCensusEvidence, CensusRefusal> {
    let budget = policy.intersection.budget();
    let mut visited_raw_pairs = 0usize;
    let mut inspected_pairs = 0usize;
    let mut shared_vertex_pairs_skipped = 0usize;
    let mut intersecting_pairs = 0usize;
    let mut complete = true;
    match policy.intersection {
        IntersectionInspection::ExhaustiveF64 { .. } => {
            'pairs: for left in 0..soup.triangles.len() {
                for right in (left + 1)..soup.triangles.len() {
                    let raw_pair =
                        pair_prefix(soup.triangles.len(), left).saturating_add(right - left - 1);
                    if visited_raw_pairs == budget {
                        complete = false;
                        break 'pairs;
                    }
                    poll(
                        cx,
                        policy.poll_stride,
                        "shell-intersection-census",
                        raw_pair,
                    )?;
                    visited_raw_pairs += 1;
                    if triangles_share_vertex(soup.triangles[left], soup.triangles[right]) {
                        shared_vertex_pairs_skipped += 1;
                        continue;
                    }
                    inspected_pairs += 1;
                    if indexed_triangles_intersect(soup, left, right) {
                        intersecting_pairs += 1;
                    }
                }
            }
        }
        IntersectionInspection::DeterministicSampleF64 { sample_count } => {
            let raw_pairs = soup
                .triangles
                .len()
                .saturating_mul(soup.triangles.len().saturating_sub(1))
                / 2;
            complete = sample_count >= raw_pairs;
            let samples = sample_count.min(raw_pairs);
            visited_raw_pairs = samples;
            for sample in 0..samples {
                poll(cx, policy.poll_stride, "shell-intersection-sample", sample)?;
                let target = if samples == raw_pairs {
                    sample
                } else {
                    sample.saturating_mul(raw_pairs) / samples
                };
                let (left, right) = unrank_pair(soup.triangles.len(), target);
                if triangles_share_vertex(soup.triangles[left], soup.triangles[right]) {
                    shared_vertex_pairs_skipped += 1;
                    continue;
                }
                inspected_pairs += 1;
                if indexed_triangles_intersect(soup, left, right) {
                    intersecting_pairs += 1;
                }
            }
        }
    }
    Ok(IntersectionCensusEvidence {
        requested_level: policy.intersection.requested_level(),
        pair_budget: budget,
        visited_raw_pairs,
        inspected_pairs,
        shared_vertex_pairs_skipped,
        intersecting_pairs,
        complete,
        authority: "f64-geometric-filter-no-exact-predicate-certificate",
    })
}

fn pair_prefix(n: usize, left: usize) -> usize {
    left.saturating_mul(n.saturating_mul(2).saturating_sub(left).saturating_sub(1)) / 2
}

fn unrank_pair(n: usize, rank: usize) -> (usize, usize) {
    let mut low = 0usize;
    let mut high = n.saturating_sub(1);
    while low < high {
        let middle = low + (high - low).div_ceil(2);
        if pair_prefix(n, middle) <= rank {
            low = middle;
        } else {
            high = middle - 1;
        }
    }
    let left = low.min(n.saturating_sub(2));
    (left, left + 1 + rank.saturating_sub(pair_prefix(n, left)))
}

fn triangles_share_vertex(a: [u32; 3], b: [u32; 3]) -> bool {
    a.iter().any(|vertex| b.contains(vertex))
}

fn indexed_triangles_intersect(soup: &Soup, left: usize, right: usize) -> bool {
    let left_indices = soup.triangles[left];
    let right_indices = soup.triangles[right];
    if left_indices
        .iter()
        .chain(right_indices.iter())
        .any(|&vertex| vertex as usize >= soup.positions.len())
    {
        return false;
    }
    let a = left_indices.map(|vertex| soup.positions[vertex as usize]);
    let b = right_indices.map(|vertex| soup.positions[vertex as usize]);
    points_are_finite(a) && points_are_finite(b) && triangles_intersect(a, b)
}

fn triangles_intersect(a: [Point3; 3], b: [Point3; 3]) -> bool {
    if !aabb_overlaps(a, b) {
        return false;
    }
    let normal_a = cross(subtract(a[1], a[0]), subtract(a[2], a[0]));
    let normal_b = cross(subtract(b[1], b[0]), subtract(b[2], b[0]));
    let norm_a = norm(normal_a);
    let norm_b = norm(normal_b);
    if norm_a == 0.0 || norm_b == 0.0 {
        return false;
    }
    let scale = a
        .iter()
        .chain(b.iter())
        .map(|point| point.x.abs().max(point.y.abs()).max(point.z.abs()))
        .fold(1.0_f64, f64::max);
    let epsilon = 128.0 * f64::EPSILON * scale * scale;
    let normal_cross = norm(cross(normal_a, normal_b));
    let coplanar = normal_cross <= epsilon * norm_a.max(norm_b)
        && a.iter()
            .all(|point| dot(normal_b, subtract(*point, b[0])).abs() <= epsilon)
        && b.iter()
            .all(|point| dot(normal_a, subtract(*point, a[0])).abs() <= epsilon);
    if coplanar {
        return coplanar_triangles_intersect(a, b, normal_a, epsilon);
    }
    triangle_edges(a)
        .into_iter()
        .any(|(start, end)| segment_intersects_triangle(start, end, b, epsilon))
        || triangle_edges(b)
            .into_iter()
            .any(|(start, end)| segment_intersects_triangle(start, end, a, epsilon))
}

fn aabb_overlaps(a: [Point3; 3], b: [Point3; 3]) -> bool {
    for axis in 0..3 {
        let coordinate = |point: Point3| match axis {
            0 => point.x,
            1 => point.y,
            _ => point.z,
        };
        let a_min = a.into_iter().map(coordinate).fold(f64::INFINITY, f64::min);
        let a_max = a
            .into_iter()
            .map(coordinate)
            .fold(f64::NEG_INFINITY, f64::max);
        let b_min = b.into_iter().map(coordinate).fold(f64::INFINITY, f64::min);
        let b_max = b
            .into_iter()
            .map(coordinate)
            .fold(f64::NEG_INFINITY, f64::max);
        if a_max < b_min || b_max < a_min {
            return false;
        }
    }
    true
}

fn triangle_edges(triangle: [Point3; 3]) -> [(Point3, Point3); 3] {
    [
        (triangle[0], triangle[1]),
        (triangle[1], triangle[2]),
        (triangle[2], triangle[0]),
    ]
}

fn segment_intersects_triangle(
    start: Point3,
    end: Point3,
    triangle: [Point3; 3],
    epsilon: f64,
) -> bool {
    let direction = subtract(end, start);
    let edge_1 = subtract(triangle[1], triangle[0]);
    let edge_2 = subtract(triangle[2], triangle[0]);
    let p = cross(direction, edge_2);
    let determinant = dot(edge_1, p);
    if determinant.abs() <= epsilon {
        return false;
    }
    let inverse = determinant.recip();
    let t_vector = subtract(start, triangle[0]);
    let u = dot(t_vector, p) * inverse;
    if u < -epsilon || u > 1.0 + epsilon {
        return false;
    }
    let q = cross(t_vector, edge_1);
    let v = dot(direction, q) * inverse;
    if v < -epsilon || u + v > 1.0 + epsilon {
        return false;
    }
    let segment_t = dot(edge_2, q) * inverse;
    segment_t >= -epsilon && segment_t <= 1.0 + epsilon
}

fn coplanar_triangles_intersect(
    a: [Point3; 3],
    b: [Point3; 3],
    normal: [f64; 3],
    epsilon: f64,
) -> bool {
    let drop_axis = if normal[0].abs() >= normal[1].abs() && normal[0].abs() >= normal[2].abs() {
        0
    } else if normal[1].abs() >= normal[2].abs() {
        1
    } else {
        2
    };
    let project = |point: Point3| match drop_axis {
        0 => [point.y, point.z],
        1 => [point.x, point.z],
        _ => [point.x, point.y],
    };
    let a2 = a.map(project);
    let b2 = b.map(project);
    for (a_start, a_end) in edges_2d(a2) {
        for (b_start, b_end) in edges_2d(b2) {
            if segments_intersect_2d(a_start, a_end, b_start, b_end, epsilon) {
                return true;
            }
        }
    }
    point_in_triangle_2d(a2[0], b2, epsilon) || point_in_triangle_2d(b2[0], a2, epsilon)
}

fn edges_2d(triangle: [[f64; 2]; 3]) -> [([f64; 2], [f64; 2]); 3] {
    [
        (triangle[0], triangle[1]),
        (triangle[1], triangle[2]),
        (triangle[2], triangle[0]),
    ]
}

fn orient_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]).mul_add(c[1] - a[1], -((b[1] - a[1]) * (c[0] - a[0])))
}

fn segments_intersect_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2], epsilon: f64) -> bool {
    let ab_c = orient_2d(a, b, c);
    let ab_d = orient_2d(a, b, d);
    let cd_a = orient_2d(c, d, a);
    let cd_b = orient_2d(c, d, b);
    (ab_c <= epsilon && ab_d >= -epsilon || ab_d <= epsilon && ab_c >= -epsilon)
        && (cd_a <= epsilon && cd_b >= -epsilon || cd_b <= epsilon && cd_a >= -epsilon)
        && ranges_overlap(a[0], b[0], c[0], d[0], epsilon)
        && ranges_overlap(a[1], b[1], c[1], d[1], epsilon)
}

fn ranges_overlap(a: f64, b: f64, c: f64, d: f64, epsilon: f64) -> bool {
    a.min(b) <= c.max(d) + epsilon && c.min(d) <= a.max(b) + epsilon
}

fn point_in_triangle_2d(point: [f64; 2], triangle: [[f64; 2]; 3], epsilon: f64) -> bool {
    let signs = [
        orient_2d(triangle[0], triangle[1], point),
        orient_2d(triangle[1], triangle[2], point),
        orient_2d(triangle[2], triangle[0], point),
    ];
    signs.iter().all(|&sign| sign >= -epsilon) || signs.iter().all(|&sign| sign <= epsilon)
}

fn append_findings_json(output: &mut String, findings: &[ImportDefect]) {
    for (index, finding) in findings.iter().enumerate() {
        if index > 0 {
            output.push(',');
        }
        let _ = write!(
            output,
            "{{\"class\":\"{}\",\"count\":{}}}",
            finding.class, finding.count
        );
    }
}

fn finding_classes<'a>(before: &'a [ImportDefect], after: &'a [ImportDefect]) -> Vec<&'a str> {
    let mut classes = BTreeSet::new();
    classes.extend(before.iter().map(|finding| finding.class));
    classes.extend(after.iter().map(|finding| finding.class));
    classes.into_iter().collect()
}

fn canonical_f64(value: f64) -> String {
    format!("{value:.17e}")
}

fn optional_f64_json(value: Option<f64>) -> String {
    value.map_or_else(|| "null".to_string(), canonical_f64)
}

fn json_escape(value: &str) -> String {
    let mut output = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(output, "\\u{:04x}", character as u32);
            }
            character => output.push(character),
        }
    }
    output
}

/// Wrap a parsed soup in quarantine with its census and receipt.
#[must_use]
pub fn quarantine(soup: Soup, format: &'static str, source: &[u8]) -> Quarantined<Soup> {
    let defects = census(&soup);
    Quarantined {
        source_receipt: ImportReceipt {
            format,
            source_hash: fs_obs::fnv1a64(source),
            parser_version: crate::VERSION,
            parsed: (soup.positions.len(), soup.triangles.len()),
        },
        raw: soup,
        defects,
    }
}

/// Promote a quarantined mesh: run the fs-rep-mesh repair suite, re-run
/// the validity census, and REFUSE if blocking defects remain. On
/// success the value becomes `Evidence<Soup>` with the full receipt in
/// its provenance chain.
///
/// # Errors
/// [`PromotionRefusal`] with blocking defects + actionable fixes.
pub fn promote(
    q: Quarantined<Soup>,
    max_hole_edges: usize,
) -> Result<(Evidence<Soup>, String), Box<PromotionRefusal>> {
    let outcome = repair(q.raw, max_hole_edges);
    let post = census(&outcome.soup);
    let blocking: Vec<String> = post
        .iter()
        .filter(|d| d.class != "unreferenced-vertex") // cosmetic, repair keeps welded verts
        .map(|d| format!("{} x{}", d.class, d.count))
        .collect();
    if !blocking.is_empty() {
        let receipt_json = q.source_receipt.to_json(&post, "refused");
        return Err(Box::new(PromotionRefusal {
            fixes: post
                .iter()
                .map(|d| match d.class {
                    "non-manifold-or-open" => format!(
                        "{} unrepaired: increase max_hole_edges (currently {max_hole_edges}) \
                         or route through the SDF re-mesh pipeline",
                        d.class
                    ),
                    other => format!("{other} survived repair: report a repair-suite gap"),
                })
                .collect(),
            blocking,
            receipt_json,
        }));
    }
    // Trusted: exact numerics (the mesh IS the value), receipt-chained
    // provenance, no model claims (geometry, not physics).
    let receipt_json = q.source_receipt.to_json(&q.defects, "promoted");
    let mut canon = receipt_json.clone();
    let _ = write!(canon, ";repairs={}", outcome.receipts_json());
    let provenance = ProvenanceHash::of_bytes(canon.as_bytes());
    let n_tris = outcome.soup.triangles.len();
    #[allow(clippy::cast_precision_loss)]
    let qoi = n_tris as f64;
    Ok((
        Evidence {
            value: outcome.soup,
            qoi,
            numerical: NumericalCertificate::exact(qoi),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        },
        receipt_json,
    ))
}

/// Promote under an explicit tolerance, intersection-cost, cancellation, and
/// residual-defect policy.
///
/// Unlike [`promote`], the success receipt is typed and includes pre/post
/// censuses, every repair-suite action, per-class before/after counts, the
/// project threshold profile, and the diagnostic geometry-budget input.
///
/// # Errors
/// [`ImportPromotionError::Census`] for policy/cancellation failure, or
/// [`ImportPromotionError::Refused`] when residuals exceed the receipted
/// thresholds.
pub fn promote_with_policy(
    q: Quarantined<Soup>,
    policy: ImportPromotionPolicy,
    cx: &Cx<'_>,
) -> Result<(Evidence<Soup>, ImportPromotionReceipt), ImportPromotionError> {
    let before =
        census_with_policy(&q.raw, policy.census, cx).map_err(ImportPromotionError::Census)?;

    // The current fs-rep-mesh repair code indexes every face. Refuse unsafe
    // source structure before entering it; a threshold cannot opt into a panic.
    if before.count("invalid-face-index") > 0 || before.count("non-finite-vertex") > 0 {
        let receipt = ImportPromotionReceipt {
            source: q.source_receipt,
            policy,
            before: before.clone(),
            repairs: Vec::new(),
            after: before,
            trust: "refused",
        };
        let receipt_json = receipt.to_json();
        return Err(ImportPromotionError::Refused(Box::new(PromotionRefusal {
            blocking: vec!["unsafe source structure cannot enter the repair suite".to_string()],
            fixes: vec![
                "remove invalid face indices and non-finite vertices at the source adapter"
                    .to_string(),
            ],
            receipt_json,
        })));
    }

    checkpoint(cx, "pre-repair", 0).map_err(ImportPromotionError::Census)?;
    let outcome = repair(q.raw, policy.max_hole_edges);
    checkpoint(cx, "post-repair", outcome.soup.triangles.len())
        .map_err(ImportPromotionError::Census)?;
    let after = census_with_policy(&outcome.soup, policy.census, cx)
        .map_err(ImportPromotionError::Census)?;

    let mut blocking = Vec::new();
    let mut fixes = Vec::new();
    for finding in &after.findings {
        let accepted = policy.thresholds.maximum_for(finding.class);
        if finding.count > accepted {
            blocking.push(format!(
                "{} x{} exceeds accepted maximum {}",
                finding.class, finding.count, accepted
            ));
            fixes.push(fix_for_finding(
                finding.class,
                policy.max_hole_edges,
                accepted,
            ));
        }
    }
    if policy.thresholds.require_complete_intersection_census && !after.intersection.complete {
        blocking.push(format!(
            "intersection-census-incomplete: visited {} raw pairs and inspected {} non-adjacent \
             pairs under budget {}",
            after.intersection.visited_raw_pairs,
            after.intersection.inspected_pairs,
            after.intersection.pair_budget
        ));
        fixes.push(
            "raise the exhaustive pair budget to cover the mesh or select a scoping profile that \
             explicitly accepts sampled intersection evidence"
                .to_string(),
        );
    }

    let mut receipt = ImportPromotionReceipt {
        source: q.source_receipt,
        policy,
        before,
        repairs: outcome.receipts.clone(),
        after,
        trust: if blocking.is_empty() {
            "promoted"
        } else {
            "refused"
        },
    };
    let receipt_json = receipt.to_json();
    if !blocking.is_empty() {
        return Err(ImportPromotionError::Refused(Box::new(PromotionRefusal {
            blocking,
            fixes,
            receipt_json,
        })));
    }

    // The complete promotion receipt, rather than merely the source census,
    // is the provenance preimage.
    let provenance = ProvenanceHash::of_bytes(receipt_json.as_bytes());
    let n_tris = outcome.soup.triangles.len();
    #[allow(clippy::cast_precision_loss)]
    let qoi = n_tris as f64;
    receipt.trust = "promoted";
    Ok((
        Evidence {
            value: outcome.soup,
            qoi,
            numerical: NumericalCertificate::exact(qoi),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        },
        receipt,
    ))
}

fn fix_for_finding(class: &str, max_hole_edges: usize, accepted: usize) -> String {
    match class {
        "non-manifold-or-open" | "near-boundary-loop-gap" => format!(
            "{class} exceeds the project maximum {accepted}: increase max_hole_edges from \
             {max_hole_edges} only for a justified local fill, or route through the SDF re-mesh \
             pipeline"
        ),
        "small-edge" | "sliver-face" => format!(
            "{class} exceeds the project maximum {accepted}: heal under an explicit \
             tolerance-consuming weld/collapse operation or select a receipted scoping threshold"
        ),
        "shell-overlap-or-self-intersection" => format!(
            "{class} exceeds the project maximum {accepted}: resolve the intersecting shells or \
             use the SDF re-mesh pipeline; sampled evidence is not a clean certificate"
        ),
        "invalid-face-index" | "non-finite-vertex" => {
            format!("{class} is structurally unsafe: repair the source adapter")
        }
        other => format!(
            "{other} exceeds the project maximum {accepted}: repair the source or report a \
             repair-suite gap"
        ),
    }
}

/// Convenience: parse-with-format + quarantine in one step.
///
/// # Errors
/// [`IoError`] from the parser (quarantine itself cannot fail).
pub fn import_mesh(bytes: &[u8], format: &'static str) -> Result<Quarantined<Soup>, IoError> {
    let soup = match format {
        "stl" => crate::stl::read_stl(bytes)?,
        "obj" => {
            let text = core::str::from_utf8(bytes).map_err(|e| IoError::Malformed {
                at: e.valid_up_to(),
                what: "OBJ must be UTF-8".to_string(),
            })?;
            crate::obj::read_obj(text)?
        }
        "ply" => crate::ply::read_ply(bytes)?,
        other => {
            return Err(IoError::Unsupported {
                what: format!("format {other:?} (stl/obj/ply)"),
            });
        }
    };
    Ok(quarantine(soup, format, bytes))
}

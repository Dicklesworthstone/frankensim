//! Retained real-supplier CAD corpus contract and deterministic import scorecard.
//!
//! The corpus is an external-reality boundary, not a fixture generator. Every
//! row binds exact source bytes to a pinned upstream revision and object
//! identity, a redistribution record, a quality tier, and an expected outcome. An
//! annotation affects lane authority only after an identified human locks it;
//! proposed annotations remain visible but make the standing lane fail closed.
//!
//! Mesh rows exercise the existing parser -> quarantine -> tolerance-aware
//! census -> repair -> promotion path. STEP rows first pass the bounded Part-21
//! parser and strict triangular `FACETED_BREP` decoder, then enter that same
//! quarantine path. A refusal before materialization is a legitimate observed
//! outcome with a stable stage code; it is never relabeled as a clean import.
//!
//! Determinism class: exact input bytes, manifest bytes, policy, and `Cx` mode
//! determine a canonical JSON artifact. Cancellation at case boundaries or in
//! the existing import kernels refuses the entire scorecard publication.
//!
//! No-claim boundary: a sampled intersection census is standing diagnostic
//! evidence, not a proof of self-intersection freedom. Corpus pass rates measure
//! only the retained population and must not be generalized to all supplier CAD.

use crate::quarantine::{
    CensusRefusal, ImportCensusPolicy, ImportPromotionError, ImportPromotionPolicy,
    ImportPromotionReceipt, ImportRefusalThresholds, IntersectionInspection, Quarantined,
    import_mesh, promote_with_policy, quarantine,
};
use crate::step_faceted::{StepFacetedLimits, StepFacetedRefusal};
use crate::{IoError, decode_faceted_brep_with_limits, parse_step};
use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Component, Path};

/// Stable semantic identity of the retained supplier-corpus contract.
pub const SUPPLIER_CORPUS_SEMANTICS_VERSION: &str = "fs-io-supplier-corpus-v1";
/// Stable semantic identity of the compact dashboard-source projection.
pub const SUPPLIER_CORPUS_SUMMARY_SEMANTICS_VERSION: &str = "fs-io-supplier-corpus-summary-v1";

/// Minimum real-file population required before the standing metric can pass.
pub const MINIMUM_RETAINED_SUPPLIER_CASES: usize = 20;

/// Exact tab-separated manifest header.
pub const SUPPLIER_CORPUS_MANIFEST_HEADER: &str = concat!(
    "case_id\trelative_path\tformat\tquality_tier\tsource_kind\tsource_origin\t",
    "source_revision\tsource_path\tsource_object_identity\tlicense_spdx\t",
    "license_url\tcontent_blake3\texpected_outcome\texpected_detail\treview_state\t",
    "reviewer\treviewed_at\tannotation_revision\tjustification"
);

const SCORECARD_IDENTITY_DOMAIN: &str = "fs-io supplier corpus scorecard v1";
/// Domain separator for the exact manifest-byte identity.
pub const SUPPLIER_CORPUS_MANIFEST_IDENTITY_DOMAIN: &str =
    "fs-io supplier corpus manifest bytes v1";
const MAX_SUPPLIER_CORPUS_MANIFEST_BYTES: usize = 8 * 1024 * 1024;
const MAX_SUPPLIER_CORPUS_CASES: usize = 10_000;

/// File format admitted by the retained corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SupplierCadFormat {
    /// Binary or ASCII stereolithography.
    Stl,
    /// Wavefront polygon mesh.
    Obj,
    /// ASCII polygon file format.
    Ply,
    /// ISO 10303-21 source, restricted to the native triangular faceted subset.
    Step,
}

impl SupplierCadFormat {
    /// Stable lowercase manifest spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Stl => "stl",
            Self::Obj => "obj",
            Self::Ply => "ply",
            Self::Step => "step",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "stl" => Some(Self::Stl),
            "obj" => Some(Self::Obj),
            "ply" => Some(Self::Ply),
            "step" => Some(Self::Step),
            _ => None,
        }
    }
}

/// Upstream revision system used to reconcile retained bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CorpusSourceKind {
    /// A file pinned by repository commit and Git blob SHA-1.
    Git,
    /// A versioned HTTPS resource pinned by provider revision and SHA-256.
    HttpSnapshot,
}

impl CorpusSourceKind {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Git => "git",
            Self::HttpSnapshot => "http-snapshot",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "git" => Some(Self::Git),
            "http-snapshot" => Some(Self::HttpSnapshot),
            _ => None,
        }
    }
}

/// Review-assigned source-quality stratum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CorpusQualityTier {
    /// Direct parametric-system export with no known defect.
    CleanParametricExport,
    /// High-detail tessellation exported by a real design workflow.
    TessellationHeavyExport,
    /// Scan-derived or externally repaired mesh.
    ScannedOrRepairedMesh,
    /// Source retained specifically because it is known to refuse or regress.
    KnownBroken,
}

impl CorpusQualityTier {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::CleanParametricExport => "clean-parametric-export",
            Self::TessellationHeavyExport => "tessellation-heavy-export",
            Self::ScannedOrRepairedMesh => "scanned-or-repaired-mesh",
            Self::KnownBroken => "known-broken",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "clean-parametric-export" => Some(Self::CleanParametricExport),
            "tessellation-heavy-export" => Some(Self::TessellationHeavyExport),
            "scanned-or-repaired-mesh" => Some(Self::ScannedOrRepairedMesh),
            "known-broken" => Some(Self::KnownBroken),
            _ => None,
        }
    }
}

/// Counted diagnostic expected from the pre-repair census or repair history.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CorpusExpectedFinding {
    /// Stable namespace-qualified class, such as `census:duplicate-face` or
    /// `repair:flipped-patch`.
    pub class: String,
    /// Exact expected count.
    pub count: usize,
}

/// Locked or proposed expected outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorpusExpectedOutcome {
    /// Parser, census, and promotion succeed without findings or repair actions.
    ImportClean,
    /// Promotion succeeds after the exact named diagnostic/repair multiset.
    Repair {
        /// Exact sorted finding multiset.
        findings: Vec<CorpusExpectedFinding>,
    },
    /// Import refuses at the named stable stage code.
    Refuse {
        /// Stable refusal code, deliberately excluding volatile diagnostics.
        code: String,
    },
}

impl CorpusExpectedOutcome {
    fn kind(&self) -> &'static str {
        match self {
            Self::ImportClean => "clean",
            Self::Repair { .. } => "repair",
            Self::Refuse { .. } => "refuse",
        }
    }

    fn detail(&self) -> String {
        match self {
            Self::ImportClean => String::new(),
            Self::Repair { findings } => findings
                .iter()
                .map(|finding| format!("{}={}", finding.class, finding.count))
                .collect::<Vec<_>>()
                .join(";"),
            Self::Refuse { code } => code.clone(),
        }
    }
}

/// Review authority attached to one source-tier and expected-outcome annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorpusAnnotationAuthority {
    /// Machine- or agent-proposed baseline. It may be inspected, but cannot
    /// make the standing lane green.
    Proposed,
    /// Human-reviewed baseline locked at an explicit revision.
    HumanLocked {
        /// Identified reviewer.
        reviewer: String,
        /// Review date in `YYYY-MM-DD` form.
        reviewed_at: String,
        /// Positive golden annotation revision.
        revision: u32,
    },
}

impl CorpusAnnotationAuthority {
    /// Whether this annotation may participate in pass/fail agreement.
    #[must_use]
    pub const fn is_locked(&self) -> bool {
        matches!(self, Self::HumanLocked { .. })
    }
}

/// One exact retained supplier input and its annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupplierCorpusCase {
    /// Stable corpus-local identifier.
    pub case_id: String,
    /// Normalized path relative to the corpus root.
    pub relative_path: String,
    /// Declared parser format.
    pub format: SupplierCadFormat,
    /// Proposed or human-reviewed source-quality stratum.
    pub quality_tier: CorpusQualityTier,
    /// Upstream revision system.
    pub source_kind: CorpusSourceKind,
    /// Canonical upstream repository or provider URL.
    pub source_origin: String,
    /// Exact repository commit or provider revision identifier.
    pub source_revision: String,
    /// Upstream path at the pinned revision.
    pub source_path: String,
    /// Exact upstream object identity: Git blob SHA-1 for `git`, SHA-256 for
    /// `http-snapshot`.
    pub source_object_identity: String,
    /// SPDX license identifier or explicit permission token.
    pub license_spdx: String,
    /// Immutable or repository-scoped license/permission record.
    pub license_url: String,
    /// Collision-resistant identity of the retained local bytes.
    pub content_blake3: ContentHash,
    /// Proposed or locked expected result.
    pub expected: CorpusExpectedOutcome,
    /// Review authority governing both the quality tier and expected outcome.
    pub annotation_authority: CorpusAnnotationAuthority,
    /// Human-readable baseline rationale. This is never used for comparison.
    pub justification: String,
}

/// Parsed exact manifest plus its source-byte identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusManifest {
    /// Collision-resistant identity of the exact manifest bytes.
    pub manifest_identity: ContentHash,
    /// Cases in required increasing `case_id` order.
    pub cases: Vec<SupplierCorpusCase>,
}

impl CorpusManifest {
    /// True only when every row has human-locked annotation authority.
    #[must_use]
    pub fn annotations_locked(&self) -> bool {
        self.cases
            .iter()
            .all(|case| case.annotation_authority.is_locked())
    }

    /// True only when the retained population meets the standing minimum.
    #[must_use]
    pub fn meets_minimum_population(&self) -> bool {
        self.cases.len() >= MINIMUM_RETAINED_SUPPLIER_CASES
    }

    /// True only when STEP and at least one supported mesh format are present.
    #[must_use]
    pub fn covers_required_formats(&self) -> bool {
        let has_step = self
            .cases
            .iter()
            .any(|case| case.format == SupplierCadFormat::Step);
        let has_mesh = self.cases.iter().any(|case| {
            matches!(
                case.format,
                SupplierCadFormat::Stl | SupplierCadFormat::Obj | SupplierCadFormat::Ply
            )
        });
        has_step && has_mesh
    }

    /// True only when every declared real-supplier quality stratum is present.
    #[must_use]
    pub fn covers_required_quality_tiers(&self) -> bool {
        [
            CorpusQualityTier::CleanParametricExport,
            CorpusQualityTier::TessellationHeavyExport,
            CorpusQualityTier::ScannedOrRepairedMesh,
            CorpusQualityTier::KnownBroken,
        ]
        .into_iter()
        .all(|tier| self.cases.iter().any(|case| case.quality_tier == tier))
    }
}

/// Structured strict-manifest refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusManifestError {
    /// One-based source line, or zero for whole-document admission.
    pub line: usize,
    /// Stable field or document stage.
    pub field: &'static str,
    /// Actionable diagnosis.
    pub reason: String,
}

impl std::fmt::Display for CorpusManifestError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.line == 0 {
            write!(
                formatter,
                "supplier corpus manifest refused at {}: {}",
                self.field, self.reason
            )
        } else {
            write!(
                formatter,
                "supplier corpus manifest line {} field {} refused: {}",
                self.line, self.field, self.reason
            )
        }
    }
}

impl std::error::Error for CorpusManifestError {}

/// Standing scorecard policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupplierCorpusPolicy {
    /// Mesh quarantine/promotion envelope.
    pub promotion: ImportPromotionPolicy,
    /// Native strict-faceted STEP decoder bounds.
    pub step_limits: StepFacetedLimits,
}

impl SupplierCorpusPolicy {
    /// Construct the standing diagnostic policy.
    ///
    /// The deterministic sample prevents quadratic work on real tessellations.
    /// Consequently this policy does not require a complete intersection census,
    /// and every receipt preserves that no-claim boundary.
    ///
    /// # Errors
    /// Returns a structured refusal only if the compile-time standing constants
    /// cease to satisfy `ImportCensusPolicy` admission.
    pub fn try_standing_lane() -> Result<Self, CensusRefusal> {
        let census = ImportCensusPolicy::try_new(
            1.0e-6,
            IntersectionInspection::DeterministicSampleF64 {
                sample_count: 4_096,
            },
            1_024,
        )?;
        let mut thresholds = ImportRefusalThresholds::validation_grade();
        thresholds.require_complete_intersection_census = false;
        Ok(Self {
            promotion: ImportPromotionPolicy {
                profile: "supplier-corpus-standing-v1",
                max_hole_edges: 32,
                census,
                thresholds,
            },
            step_limits: StepFacetedLimits::default(),
        })
    }
}

/// Coarse standing metric verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorpusVerdict {
    /// Promoted with no pre-repair findings and no repair actions.
    Clean,
    /// Promoted after at least one diagnostic finding or repair action.
    Repaired,
    /// Refused before trusted promotion.
    Refused,
}

impl CorpusVerdict {
    /// Stable lowercase scorecard spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Clean => "clean",
            Self::Repaired => "repaired",
            Self::Refused => "refused",
        }
    }
}

/// One case's actual import evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusObservation {
    /// Coarse verdict.
    pub verdict: CorpusVerdict,
    /// Exact sorted pre-census and repair-history counts.
    pub findings: Vec<CorpusExpectedFinding>,
    /// Stable refusal stage code, present only for refusals.
    pub refusal_code: Option<String>,
    /// Canonical underlying import/promotion receipt or refusal receipt.
    pub receipt_json: String,
    /// Additional deterministic diagnostic; comparison never depends on it.
    pub diagnostic: String,
}

/// Relationship between one observation and its proposed/locked annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnnotationAgreement {
    /// Human-locked annotation matched exactly.
    Match,
    /// Human-locked annotation disagreed with the observation.
    Mismatch {
        /// Deterministic explanation.
        reason: String,
    },
    /// Annotation is still proposed. `proposed_matches` is advisory only.
    Unreviewed {
        /// Whether the proposal currently agrees with the observation.
        proposed_matches: bool,
        /// Deterministic advisory explanation.
        reason: String,
    },
}

impl AnnotationAgreement {
    fn as_str(&self) -> &'static str {
        match self {
            Self::Match => "match",
            Self::Mismatch { .. } => "mismatch",
            Self::Unreviewed { .. } => "unreviewed",
        }
    }

    fn reason(&self) -> &str {
        match self {
            Self::Match => "",
            Self::Mismatch { reason } | Self::Unreviewed { reason, .. } => reason,
        }
    }
}

/// One deterministic scorecard row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusScorecardRow {
    /// Stable case identifier.
    pub case_id: String,
    /// Normalized retained path relative to the corpus root.
    pub relative_path: String,
    /// Exact retained-byte identity.
    pub content_blake3: ContentHash,
    /// Declared parser format.
    pub format: SupplierCadFormat,
    /// Upstream revision system.
    pub source_kind: CorpusSourceKind,
    /// Canonical upstream repository or provider URL.
    pub source_origin: String,
    /// Exact upstream revision identifier.
    pub source_revision: String,
    /// Upstream path at the pinned revision.
    pub source_path: String,
    /// Exact upstream object identity under the declared source kind.
    pub source_object_identity: String,
    /// SPDX license identifier or permission token.
    pub license_spdx: String,
    /// Immutable or repository-scoped license/permission record.
    pub license_url: String,
    /// Source quality stratum.
    pub quality_tier: CorpusQualityTier,
    /// Expected outcome.
    pub expected: CorpusExpectedOutcome,
    /// Review authority governing both the quality tier and expected outcome.
    pub annotation_authority: CorpusAnnotationAuthority,
    /// Actual observation.
    pub observed: CorpusObservation,
    /// Locked/proposed comparison.
    pub agreement: AnnotationAgreement,
}

/// Authority-aware counts exported to program-level metric consumers.
///
/// The outcome counts include only rows whose annotation is human-locked.
/// This makes `reviewed` the only valid denominator for clean, repaired, and
/// refused rates. Proposed annotations remain visible in the full
/// [`CorpusScorecard`], but cannot manufacture a program metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImportScorecardMetrics {
    total: usize,
    reviewed: usize,
    clean: usize,
    repaired: usize,
    refused: usize,
    annotation_mismatches: usize,
}

impl ImportScorecardMetrics {
    /// Total retained files, regardless of annotation authority.
    #[must_use]
    pub const fn total(self) -> usize {
        self.total
    }

    /// Files whose expected outcome has human-locked review authority.
    #[must_use]
    pub const fn reviewed(self) -> usize {
        self.reviewed
    }

    /// Human-reviewed files observed to import cleanly.
    #[must_use]
    pub const fn clean(self) -> usize {
        self.clean
    }

    /// Human-reviewed files observed to import after repair.
    #[must_use]
    pub const fn repaired(self) -> usize {
        self.repaired
    }

    /// Human-reviewed files observed to be refused.
    #[must_use]
    pub const fn refused(self) -> usize {
        self.refused
    }

    /// Observations disagreeing with a human-locked expected outcome.
    #[must_use]
    pub const fn annotation_mismatches(self) -> usize {
        self.annotation_mismatches
    }
}

/// Deterministic aggregate scorecard.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusScorecard {
    /// Exact manifest-byte identity.
    pub manifest_identity: ContentHash,
    /// Per-file rows in manifest order.
    pub rows: Vec<CorpusScorecardRow>,
    /// Clean import count.
    pub clean: usize,
    /// Repaired import count.
    pub repaired: usize,
    /// Refusal count.
    pub refused: usize,
    /// Human-locked annotation mismatch count.
    pub mismatches: usize,
    /// Advisory mismatch count among annotations that still lack human authority.
    pub proposed_mismatches: usize,
    /// Proposed/unreviewed annotation count.
    pub unreviewed: usize,
}

impl CorpusScorecard {
    /// Project authority-aware counts for program-level metrics.
    ///
    /// The full scorecard counts every observation for corpus diagnostics.
    /// This projection deliberately filters outcome counts to human-locked
    /// rows so an agent-proposed baseline can never become a rate denominator.
    #[must_use]
    pub fn import_metrics(&self) -> ImportScorecardMetrics {
        let mut clean = 0;
        let mut repaired = 0;
        let mut refused = 0;
        for row in self
            .rows
            .iter()
            .filter(|row| row.annotation_authority.is_locked())
        {
            match row.observed.verdict {
                CorpusVerdict::Clean => clean += 1,
                CorpusVerdict::Repaired => repaired += 1,
                CorpusVerdict::Refused => refused += 1,
            }
        }
        ImportScorecardMetrics {
            total: self.rows.len(),
            reviewed: clean + repaired + refused,
            clean,
            repaired,
            refused,
            annotation_mismatches: self.mismatches,
        }
    }

    /// Canonical compact projection for a tracked dashboard source.
    ///
    /// This binds the exact manifest and full per-file scorecard identities.
    /// Whole-population observations remain available for diagnostics, while
    /// the `reviewed` object is the only authority permitted as a rate source.
    #[must_use]
    pub fn summary_json(&self) -> String {
        let metrics = self.import_metrics();
        format!(
            "{{\"schema\":1,\"semantics\":\"{}\",\
             \"manifest_identity_domain\":\"{}\",\"manifest_identity\":\"{}\",\
             \"scorecard_identity\":\"{}\",\"population\":{{\"total\":{},\"clean\":{},\
             \"repaired\":{},\"refused\":{},\"unreviewed\":{},\
             \"annotation_mismatch\":{},\"proposal_mismatch\":{}}},\
             \"reviewed\":{{\"total\":{},\"clean\":{},\"repaired\":{},\"refused\":{},\
             \"annotation_mismatch\":{}}},\
             \"authority\":\"human-locked-only-dashboard-denominator\"}}",
            SUPPLIER_CORPUS_SUMMARY_SEMANTICS_VERSION,
            SUPPLIER_CORPUS_MANIFEST_IDENTITY_DOMAIN,
            self.manifest_identity,
            self.artifact_identity(),
            self.rows.len(),
            self.clean,
            self.repaired,
            self.refused,
            self.unreviewed,
            self.mismatches,
            self.proposed_mismatches,
            metrics.reviewed(),
            metrics.clean(),
            metrics.repaired(),
            metrics.refused(),
            metrics.annotation_mismatches(),
        )
    }

    /// True only when population, review authority, and exact agreement all pass.
    #[must_use]
    pub fn lane_passes(&self) -> bool {
        self.rows.len() >= MINIMUM_RETAINED_SUPPLIER_CASES
            && self.covers_required_formats()
            && self.covers_required_quality_tiers()
            && self.mismatches == 0
            && self.unreviewed == 0
    }

    /// True only when STEP and at least one supported mesh format are present.
    #[must_use]
    pub fn covers_required_formats(&self) -> bool {
        let has_step = self
            .rows
            .iter()
            .any(|row| row.format == SupplierCadFormat::Step);
        let has_mesh = self.rows.iter().any(|row| {
            matches!(
                row.format,
                SupplierCadFormat::Stl | SupplierCadFormat::Obj | SupplierCadFormat::Ply
            )
        });
        has_step && has_mesh
    }

    /// True only when every declared real-supplier quality stratum is present.
    #[must_use]
    pub fn covers_required_quality_tiers(&self) -> bool {
        [
            CorpusQualityTier::CleanParametricExport,
            CorpusQualityTier::TessellationHeavyExport,
            CorpusQualityTier::ScannedOrRepairedMesh,
            CorpusQualityTier::KnownBroken,
        ]
        .into_iter()
        .all(|tier| self.rows.iter().any(|row| row.quality_tier == tier))
    }

    /// Domain-separated identity of the canonical payload (excluding the
    /// self-referential artifact-identity field).
    #[must_use]
    pub fn artifact_identity(&self) -> ContentHash {
        hash_domain(SCORECARD_IDENTITY_DOMAIN, self.payload_json().as_bytes())
    }

    /// Canonical machine-readable scorecard.
    #[must_use]
    pub fn to_json(&self) -> String {
        let payload = self.payload_json();
        format!(
            "{{\"kind\":\"supplier-cad-import-scorecard\",\"semantics\":\"{}\",\
             \"artifact_identity\":\"{}\",\"payload\":{}}}",
            SUPPLIER_CORPUS_SEMANTICS_VERSION,
            self.artifact_identity(),
            payload
        )
    }

    fn payload_json(&self) -> String {
        let total = self.rows.len();
        let mut output = format!(
            "{{\"manifest_identity\":\"{}\",\"minimum_population\":{},\
             \"population\":{},\"lane_passes\":{},\
             \"coverage\":{{\"minimum_population\":{},\"required_formats\":{},\
             \"required_quality_tiers\":{}}},\"rates\":{{\
             \"clean\":{},\"repaired\":{},\"refused\":{},\"annotation_mismatch\":{},\
             \"proposal_mismatch\":{},\"unreviewed\":{}}},\"rows\":[",
            self.manifest_identity,
            MINIMUM_RETAINED_SUPPLIER_CASES,
            total,
            self.lane_passes(),
            total >= MINIMUM_RETAINED_SUPPLIER_CASES,
            self.covers_required_formats(),
            self.covers_required_quality_tiers(),
            rate_json(self.clean, total),
            rate_json(self.repaired, total),
            rate_json(self.refused, total),
            rate_json(self.mismatches, total),
            rate_json(self.proposed_mismatches, total),
            rate_json(self.unreviewed, total)
        );
        for (index, row) in self.rows.iter().enumerate() {
            if index > 0 {
                output.push(',');
            }
            append_row_json(&mut output, row);
        }
        output.push_str("],\"authority\":\"retained-population-metric-not-universal-import-rate\",\
                         \"intersection_authority\":\"sampled-diagnostic-not-self-intersection-proof\"}");
        output
    }
}

/// Whole-run refusal. No partial scorecard should be published.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CorpusRunError {
    /// The explicit `Cx` observed cancellation.
    Cancelled {
        /// Case active at observation, or `<scorecard>` at a boundary.
        case_id: String,
        /// Stable stage.
        stage: &'static str,
    },
}

impl std::fmt::Display for CorpusRunError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cancelled { case_id, stage } => {
                write!(
                    formatter,
                    "supplier corpus scorecard cancelled for {case_id} at {stage}"
                )
            }
        }
    }
}

impl std::error::Error for CorpusRunError {}

/// Parse the strict tab-separated corpus manifest.
///
/// Blank lines and `#` comments do not create rows, although their exact bytes
/// remain bound into the manifest identity. The first remaining line must equal
/// [`SUPPLIER_CORPUS_MANIFEST_HEADER`], and data rows must be strictly
/// increasing by `case_id`.
///
/// # Errors
/// [`CorpusManifestError`] for any schema, path, provenance, hash, annotation,
/// ordering, or duplicate refusal.
pub fn parse_supplier_corpus_manifest(input: &str) -> Result<CorpusManifest, CorpusManifestError> {
    if input.len() > MAX_SUPPLIER_CORPUS_MANIFEST_BYTES {
        return Err(manifest_error(
            0,
            "manifest",
            format!(
                "manifest exceeds the admitted {}-byte bound",
                MAX_SUPPLIER_CORPUS_MANIFEST_BYTES
            ),
        ));
    }
    let mut records = input.lines().enumerate().filter(|(_, line)| {
        let trimmed = line.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });
    let Some((header_index, header)) = records.next() else {
        return Err(manifest_error(
            0,
            "header",
            "manifest has no non-comment header",
        ));
    };
    if header != SUPPLIER_CORPUS_MANIFEST_HEADER {
        return Err(manifest_error(
            header_index + 1,
            "header",
            "header does not exactly match the v1 schema",
        ));
    }

    let mut cases = Vec::new();
    let mut case_ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut content_hashes = BTreeSet::new();
    let mut previous_case_id: Option<String> = None;
    for (line_index, line) in records {
        let line_number = line_index + 1;
        if cases.len() == MAX_SUPPLIER_CORPUS_CASES {
            return Err(manifest_error(
                line_number,
                "population",
                format!(
                    "manifest exceeds the admitted {}-case bound",
                    MAX_SUPPLIER_CORPUS_CASES
                ),
            ));
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 19 {
            return Err(manifest_error(
                line_number,
                "row",
                format!("expected 19 tab-separated fields; got {}", fields.len()),
            ));
        }
        let case_id = fields[0];
        validate_case_id(case_id, line_number)?;
        if let Some(previous) = &previous_case_id
            && case_id <= previous.as_str()
        {
            return Err(manifest_error(
                line_number,
                "case_id",
                format!("rows must be strictly increasing; {case_id:?} follows {previous:?}"),
            ));
        }
        if !case_ids.insert(case_id.to_string()) {
            return Err(manifest_error(
                line_number,
                "case_id",
                format!("duplicate case identifier {case_id:?}"),
            ));
        }
        previous_case_id = Some(case_id.to_string());

        let relative_path = fields[1];
        validate_relative_path(relative_path, line_number)?;
        if !paths.insert(relative_path.to_string()) {
            return Err(manifest_error(
                line_number,
                "relative_path",
                format!("duplicate retained path {relative_path:?}"),
            ));
        }

        let format = SupplierCadFormat::parse(fields[2]).ok_or_else(|| {
            manifest_error(
                line_number,
                "format",
                "expected one of stl, obj, ply, or step",
            )
        })?;
        validate_extension(relative_path, format, line_number)?;
        let quality_tier = CorpusQualityTier::parse(fields[3]).ok_or_else(|| {
            manifest_error(line_number, "quality_tier", "unknown source-quality tier")
        })?;
        let source_kind = CorpusSourceKind::parse(fields[4]).ok_or_else(|| {
            manifest_error(line_number, "source_kind", "expected git or http-snapshot")
        })?;
        validate_https_url(fields[5], line_number, "source_origin")?;
        validate_source_provenance(source_kind, fields[6], fields[8], line_number)?;
        validate_source_path(fields[7], line_number)?;
        if fields[9].trim().is_empty() {
            return Err(manifest_error(
                line_number,
                "license_spdx",
                "license or explicit permission token is required",
            ));
        }
        validate_https_url(fields[10], line_number, "license_url")?;
        let content_blake3 = ContentHash::from_hex(fields[11]).ok_or_else(|| {
            manifest_error(
                line_number,
                "content_blake3",
                "expected exactly 64 hexadecimal BLAKE3 characters",
            )
        })?;
        if !content_hashes.insert(content_blake3) {
            return Err(manifest_error(
                line_number,
                "content_blake3",
                format!("duplicate retained byte identity {content_blake3}"),
            ));
        }
        let expected = parse_expected(fields[12], fields[13], line_number)?;
        let annotation_authority =
            parse_authority(fields[14], fields[15], fields[16], fields[17], line_number)?;
        if fields[18].trim().is_empty() {
            return Err(manifest_error(
                line_number,
                "justification",
                "baseline rationale is required",
            ));
        }

        cases.push(SupplierCorpusCase {
            case_id: case_id.to_string(),
            relative_path: relative_path.to_string(),
            format,
            quality_tier,
            source_kind,
            source_origin: fields[5].to_string(),
            source_revision: fields[6].to_string(),
            source_path: fields[7].to_string(),
            source_object_identity: fields[8].to_string(),
            license_spdx: fields[9].to_string(),
            license_url: fields[10].to_string(),
            content_blake3,
            expected,
            annotation_authority,
            justification: fields[18].to_string(),
        });
    }
    if cases.is_empty() {
        return Err(manifest_error(
            0,
            "population",
            "manifest must retain at least one case",
        ));
    }

    Ok(CorpusManifest {
        manifest_identity: hash_domain(SUPPLIER_CORPUS_MANIFEST_IDENTITY_DOMAIN, input.as_bytes()),
        cases,
    })
}

/// Run every retained case with a caller-owned byte loader.
///
/// Loader errors become deterministic `source-unavailable` refusals so a
/// missing retained asset appears as a row and annotation mismatch. `Cx`
/// cancellation instead refuses the whole artifact.
///
/// # Errors
/// [`CorpusRunError::Cancelled`] when cancellation is observed before complete
/// scorecard publication.
pub fn run_supplier_corpus(
    manifest: &CorpusManifest,
    policy: SupplierCorpusPolicy,
    cx: &Cx<'_>,
    mut load: impl FnMut(&SupplierCorpusCase) -> Result<Vec<u8>, String>,
) -> Result<CorpusScorecard, CorpusRunError> {
    checkpoint(cx, "<scorecard>", "entry")?;
    let mut rows = Vec::with_capacity(manifest.cases.len());
    let mut clean = 0usize;
    let mut repaired = 0usize;
    let mut refused = 0usize;
    let mut mismatches = 0usize;
    let mut proposed_mismatches = 0usize;
    let mut unreviewed = 0usize;

    for case in &manifest.cases {
        checkpoint(cx, &case.case_id, "case-entry")?;
        let observed = match load(case) {
            Ok(bytes) => evaluate_case(case, &bytes, policy, cx)?,
            Err(reason) => refusal_observation(
                "source-unavailable",
                format!("retained source could not be loaded: {reason}"),
                None,
            ),
        };
        match observed.verdict {
            CorpusVerdict::Clean => clean += 1,
            CorpusVerdict::Repaired => repaired += 1,
            CorpusVerdict::Refused => refused += 1,
        }
        let agreement = compare_annotation(case, &observed);
        match &agreement {
            AnnotationAgreement::Mismatch { .. } => mismatches += 1,
            AnnotationAgreement::Unreviewed {
                proposed_matches, ..
            } => {
                unreviewed += 1;
                if !proposed_matches {
                    proposed_mismatches += 1;
                }
            }
            AnnotationAgreement::Match => {}
        }
        rows.push(CorpusScorecardRow {
            case_id: case.case_id.clone(),
            relative_path: case.relative_path.clone(),
            content_blake3: case.content_blake3,
            format: case.format,
            source_kind: case.source_kind,
            source_origin: case.source_origin.clone(),
            source_revision: case.source_revision.clone(),
            source_path: case.source_path.clone(),
            source_object_identity: case.source_object_identity.clone(),
            license_spdx: case.license_spdx.clone(),
            license_url: case.license_url.clone(),
            quality_tier: case.quality_tier,
            expected: case.expected.clone(),
            annotation_authority: case.annotation_authority.clone(),
            observed,
            agreement,
        });
    }
    checkpoint(cx, "<scorecard>", "publication")?;

    Ok(CorpusScorecard {
        manifest_identity: manifest.manifest_identity,
        rows,
        clean,
        repaired,
        refused,
        mismatches,
        proposed_mismatches,
        unreviewed,
    })
}

fn evaluate_case(
    case: &SupplierCorpusCase,
    bytes: &[u8],
    policy: SupplierCorpusPolicy,
    cx: &Cx<'_>,
) -> Result<CorpusObservation, CorpusRunError> {
    let actual_hash = hash_bytes(bytes);
    if actual_hash != case.content_blake3 {
        return Ok(refusal_observation(
            "source-content-mismatch",
            format!(
                "retained bytes hash to {actual_hash}, manifest pins {}",
                case.content_blake3
            ),
            None,
        ));
    }
    let (quarantined, materialization_receipt) = match case.format {
        SupplierCadFormat::Stl | SupplierCadFormat::Obj | SupplierCadFormat::Ply => {
            match import_mesh(bytes, case.format.as_str()) {
                Ok(value) => (value, None),
                Err(error) => {
                    return Ok(io_refusal(case.format, error));
                }
            }
        }
        SupplierCadFormat::Step => match step_to_quarantine(bytes, policy.step_limits, cx)? {
            Ok((value, receipt)) => (value, Some(receipt)),
            Err(observation) => return Ok(observation),
        },
    };
    promote_observation(
        quarantined,
        policy,
        cx,
        &case.case_id,
        materialization_receipt.as_deref(),
    )
}

fn step_to_quarantine(
    bytes: &[u8],
    limits: StepFacetedLimits,
    cx: &Cx<'_>,
) -> Result<Result<(Quarantined<fs_rep_mesh::Soup>, String), CorpusObservation>, CorpusRunError> {
    let parsed = match parse_step(bytes) {
        Ok(parsed) => parsed,
        Err(error) => return Ok(Err(io_refusal(SupplierCadFormat::Step, error))),
    };
    let mut roots = Vec::new();
    for instance in &parsed.document().instances {
        if instance
            .components
            .iter()
            .any(|component| component.name == "FACETED_BREP")
        {
            roots.push(instance.id);
        }
    }
    roots.sort_unstable();
    roots.dedup();
    let root_id = match roots.as_slice() {
        [] => {
            return Ok(Err(refusal_observation(
                "step-root-missing",
                "no FACETED_BREP instance is present in the admitted Part-21 document".to_string(),
                Some(parsed.receipt().to_json()),
            )));
        }
        [root] => *root,
        _ => {
            return Ok(Err(refusal_observation(
                "step-root-ambiguous",
                format!(
                    "{} FACETED_BREP roots are present; the corpus row does not guess",
                    roots.len()
                ),
                Some(parsed.receipt().to_json()),
            )));
        }
    };
    match decode_faceted_brep_with_limits(&parsed, root_id, limits, cx) {
        Ok(decoded) => {
            let (soup, decoder_receipt) = decoded.into_parts();
            let mut quarantined = quarantine(soup, "step", bytes);
            quarantined.source_receipt.parser_version =
                crate::step_faceted::STEP_FACETED_DECODER_VERSION;
            let materialization_receipt = format!(
                "{{\"kind\":\"strict-faceted-step-materialization\",\
                 \"syntax\":{},\"decoder\":{}}}",
                parsed.receipt().to_json(),
                decoder_receipt.to_json()
            );
            Ok(Ok((quarantined, materialization_receipt)))
        }
        Err(StepFacetedRefusal::Cancelled { stage, .. }) => Err(CorpusRunError::Cancelled {
            case_id: "<step-decode>".to_string(),
            stage,
        }),
        Err(error) => {
            let code = match error {
                StepFacetedRefusal::Schema { .. } => "step-decode-schema",
                StepFacetedRefusal::Entity { .. } => "step-decode-entity",
                StepFacetedRefusal::Resource { .. } => "step-decode-resource",
                StepFacetedRefusal::Cancelled { stage, .. } => {
                    return Err(CorpusRunError::Cancelled {
                        case_id: "<step-decode>".to_string(),
                        stage,
                    });
                }
            };
            Ok(Err(refusal_observation(
                code,
                error.to_string(),
                Some(parsed.receipt().to_json()),
            )))
        }
    }
}

fn promote_observation(
    quarantined: Quarantined<fs_rep_mesh::Soup>,
    policy: SupplierCorpusPolicy,
    cx: &Cx<'_>,
    case_id: &str,
    materialization_receipt: Option<&str>,
) -> Result<CorpusObservation, CorpusRunError> {
    match promote_with_policy(quarantined, policy.promotion, cx) {
        Ok((_evidence, receipt)) => Ok(promoted_observation(receipt, materialization_receipt)),
        Err(ImportPromotionError::Census(error)) if error.cancelled => {
            Err(CorpusRunError::Cancelled {
                case_id: case_id.to_string(),
                stage: error.stage,
            })
        }
        Err(ImportPromotionError::Census(error)) => Ok(refusal_observation(
            "census-policy-refused",
            error.to_string(),
            None,
        )),
        Err(ImportPromotionError::Refused(refusal)) => Ok(refusal_observation(
            "promotion-residual-refused",
            refusal.blocking.join("; "),
            Some(combine_promotion_receipts(
                materialization_receipt,
                &refusal.receipt_json,
            )),
        )),
    }
}

fn promoted_observation(
    receipt: ImportPromotionReceipt,
    materialization_receipt: Option<&str>,
) -> CorpusObservation {
    let mut counts = BTreeMap::<String, usize>::new();
    for finding in &receipt.before.findings {
        counts.insert(format!("census:{}", finding.class), finding.count);
    }
    for repair in &receipt.repairs {
        *counts
            .entry(format!("repair:{}", repair.defect))
            .or_default() += 1;
    }
    let findings = counts
        .into_iter()
        .map(|(class, count)| CorpusExpectedFinding { class, count })
        .collect::<Vec<_>>();
    let verdict = if findings.is_empty() {
        CorpusVerdict::Clean
    } else {
        CorpusVerdict::Repaired
    };
    let receipt_json = combine_promotion_receipts(materialization_receipt, &receipt.to_json());
    CorpusObservation {
        verdict,
        findings,
        refusal_code: None,
        receipt_json,
        diagnostic: String::new(),
    }
}

fn combine_promotion_receipts(
    materialization_receipt: Option<&str>,
    promotion_receipt: &str,
) -> String {
    materialization_receipt.map_or_else(
        || promotion_receipt.to_string(),
        |materialization| {
            format!(
                "{{\"kind\":\"supplier-corpus-step-promotion\",\
                 \"materialization\":{materialization},\
                 \"promotion\":{promotion_receipt}}}"
            )
        },
    )
}

fn io_refusal(format: SupplierCadFormat, error: IoError) -> CorpusObservation {
    let suffix = match &error {
        IoError::Malformed { .. } => "parse-malformed",
        IoError::Unsupported { .. } => "parse-unsupported",
        IoError::ResourceBound { .. } => "parse-resource-bound",
        IoError::Schema { .. } => "parse-schema",
        IoError::Cancelled { .. } => "parse-cancelled",
    };
    refusal_observation(
        &format!("{}-{suffix}", format.as_str()),
        error.to_string(),
        None,
    )
}

fn refusal_observation(
    code: &str,
    diagnostic: String,
    upstream_receipt: Option<String>,
) -> CorpusObservation {
    let upstream = upstream_receipt.unwrap_or_else(|| "null".to_string());
    let receipt_json = format!(
        "{{\"kind\":\"supplier-corpus-refusal\",\"semantics\":\"{}\",\
         \"code\":\"{}\",\"diagnostic\":\"{}\",\"upstream_receipt\":{}}}",
        SUPPLIER_CORPUS_SEMANTICS_VERSION,
        json_escape(code),
        json_escape(&diagnostic),
        upstream
    );
    CorpusObservation {
        verdict: CorpusVerdict::Refused,
        findings: Vec::new(),
        refusal_code: Some(code.to_string()),
        receipt_json,
        diagnostic,
    }
}

fn compare_annotation(
    case: &SupplierCorpusCase,
    observed: &CorpusObservation,
) -> AnnotationAgreement {
    let comparison = expected_comparison(&case.expected, observed);
    match &case.annotation_authority {
        CorpusAnnotationAuthority::HumanLocked { .. } => match comparison {
            Ok(()) => AnnotationAgreement::Match,
            Err(reason) => AnnotationAgreement::Mismatch { reason },
        },
        CorpusAnnotationAuthority::Proposed => AnnotationAgreement::Unreviewed {
            proposed_matches: comparison.is_ok(),
            reason: comparison
                .err()
                .unwrap_or_else(|| "proposal agrees but lacks human lock authority".to_string()),
        },
    }
}

fn expected_comparison(
    expected: &CorpusExpectedOutcome,
    observed: &CorpusObservation,
) -> Result<(), String> {
    match (expected, observed.verdict) {
        (CorpusExpectedOutcome::ImportClean, CorpusVerdict::Clean) => Ok(()),
        (
            CorpusExpectedOutcome::Repair {
                findings: expected_findings,
            },
            CorpusVerdict::Repaired,
        ) if expected_findings == &observed.findings => Ok(()),
        (CorpusExpectedOutcome::Refuse { code }, CorpusVerdict::Refused)
            if observed.refusal_code.as_deref() == Some(code) =>
        {
            Ok(())
        }
        _ => Err(format!(
            "expected {} detail {:?}; observed {} detail {:?}",
            expected.kind(),
            expected.detail(),
            observed.verdict.as_str(),
            observed.refusal_code.as_ref().cloned().unwrap_or_else(|| {
                observed
                    .findings
                    .iter()
                    .map(|finding| format!("{}={}", finding.class, finding.count))
                    .collect::<Vec<_>>()
                    .join(";")
            })
        )),
    }
}

fn parse_expected(
    kind: &str,
    detail: &str,
    line: usize,
) -> Result<CorpusExpectedOutcome, CorpusManifestError> {
    match kind {
        "clean" => {
            if !detail.is_empty() {
                return Err(manifest_error(
                    line,
                    "expected_detail",
                    "clean expectation must have an empty detail",
                ));
            }
            Ok(CorpusExpectedOutcome::ImportClean)
        }
        "repair" => {
            if detail.is_empty() {
                return Err(manifest_error(
                    line,
                    "expected_detail",
                    "repair expectation requires exact class=count findings",
                ));
            }
            let mut findings = Vec::new();
            let mut classes = BTreeSet::new();
            for token in detail.split(';') {
                let Some((class, count)) = token.split_once('=') else {
                    return Err(manifest_error(
                        line,
                        "expected_detail",
                        format!("repair token {token:?} is not class=count"),
                    ));
                };
                if !(class.starts_with("census:") || class.starts_with("repair:"))
                    || class.len() <= "repair:".len()
                {
                    return Err(manifest_error(
                        line,
                        "expected_detail",
                        format!("finding class {class:?} must use census: or repair: namespace"),
                    ));
                }
                if !classes.insert(class) {
                    return Err(manifest_error(
                        line,
                        "expected_detail",
                        format!("duplicate finding class {class:?}"),
                    ));
                }
                let count = count.parse::<usize>().map_err(|_| {
                    manifest_error(
                        line,
                        "expected_detail",
                        format!("finding count {count:?} is not an unsigned integer"),
                    )
                })?;
                if count == 0 {
                    return Err(manifest_error(
                        line,
                        "expected_detail",
                        "finding counts must be positive",
                    ));
                }
                findings.push(CorpusExpectedFinding {
                    class: class.to_string(),
                    count,
                });
            }
            findings.sort();
            Ok(CorpusExpectedOutcome::Repair { findings })
        }
        "refuse" => {
            validate_stable_token(detail, line, "expected_detail")?;
            Ok(CorpusExpectedOutcome::Refuse {
                code: detail.to_string(),
            })
        }
        _ => Err(manifest_error(
            line,
            "expected_outcome",
            "expected one of clean, repair, or refuse",
        )),
    }
}

fn parse_authority(
    state: &str,
    reviewer: &str,
    reviewed_at: &str,
    revision: &str,
    line: usize,
) -> Result<CorpusAnnotationAuthority, CorpusManifestError> {
    match state {
        "proposed" => {
            if !reviewer.is_empty() || !reviewed_at.is_empty() || revision != "0" {
                return Err(manifest_error(
                    line,
                    "review_state",
                    "proposed rows require blank reviewer/date and revision 0",
                ));
            }
            Ok(CorpusAnnotationAuthority::Proposed)
        }
        "human-locked" => {
            if reviewer.trim().is_empty() {
                return Err(manifest_error(
                    line,
                    "reviewer",
                    "human-locked row requires an identified reviewer",
                ));
            }
            validate_date(reviewed_at, line)?;
            let revision = revision.parse::<u32>().map_err(|_| {
                manifest_error(
                    line,
                    "annotation_revision",
                    "human-locked revision must be a positive u32",
                )
            })?;
            if revision == 0 {
                return Err(manifest_error(
                    line,
                    "annotation_revision",
                    "human-locked revision must be positive",
                ));
            }
            Ok(CorpusAnnotationAuthority::HumanLocked {
                reviewer: reviewer.to_string(),
                reviewed_at: reviewed_at.to_string(),
                revision,
            })
        }
        _ => Err(manifest_error(
            line,
            "review_state",
            "expected proposed or human-locked",
        )),
    }
}

fn validate_case_id(value: &str, line: usize) -> Result<(), CorpusManifestError> {
    if value.is_empty()
        || value.len() > 96
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
        || value.starts_with('-')
        || value.ends_with('-')
    {
        return Err(manifest_error(
            line,
            "case_id",
            "expected 1..=96 lowercase ASCII letters, digits, or interior hyphens",
        ));
    }
    Ok(())
}

fn validate_relative_path(value: &str, line: usize) -> Result<(), CorpusManifestError> {
    let path = Path::new(value);
    if value.is_empty()
        || path.is_absolute()
        || !path
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
        || !value.starts_with("sources/")
    {
        return Err(manifest_error(
            line,
            "relative_path",
            "path must be normalized, relative, and rooted under sources/",
        ));
    }
    Ok(())
}

fn validate_source_path(value: &str, line: usize) -> Result<(), CorpusManifestError> {
    let path = Path::new(value);
    if value.is_empty()
        || path.is_absolute()
        || !path
            .components()
            .all(|component| matches!(component, Component::Normal(_)))
    {
        return Err(manifest_error(
            line,
            "source_path",
            "upstream path must be normalized and relative",
        ));
    }
    Ok(())
}

fn validate_source_provenance(
    kind: CorpusSourceKind,
    revision: &str,
    object_identity: &str,
    line: usize,
) -> Result<(), CorpusManifestError> {
    match kind {
        CorpusSourceKind::Git => {
            validate_hex(revision, 40, line, "source_revision")?;
            validate_hex(object_identity, 40, line, "source_object_identity")
        }
        CorpusSourceKind::HttpSnapshot => {
            if revision.is_empty()
                || revision.len() > 160
                || revision.chars().any(char::is_whitespace)
                || !revision.bytes().all(|byte| {
                    byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':' | b'@')
                })
            {
                return Err(manifest_error(
                    line,
                    "source_revision",
                    "HTTP snapshot revision must be a non-blank stable ASCII identifier",
                ));
            }
            validate_hex(object_identity, 64, line, "source_object_identity")
        }
    }
}

fn validate_extension(
    value: &str,
    format: SupplierCadFormat,
    line: usize,
) -> Result<(), CorpusManifestError> {
    let extension = Path::new(value)
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if extension != format.as_str() {
        return Err(manifest_error(
            line,
            "relative_path",
            format!(
                "extension {extension:?} does not match declared format {}",
                format.as_str()
            ),
        ));
    }
    Ok(())
}

fn validate_https_url(
    value: &str,
    line: usize,
    field: &'static str,
) -> Result<(), CorpusManifestError> {
    if !value.starts_with("https://")
        || value.len() <= "https://".len()
        || value.chars().any(char::is_whitespace)
    {
        return Err(manifest_error(
            line,
            field,
            "expected a non-blank HTTPS URL with no whitespace",
        ));
    }
    Ok(())
}

fn validate_hex(
    value: &str,
    length: usize,
    line: usize,
    field: &'static str,
) -> Result<(), CorpusManifestError> {
    if value.len() != length
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(manifest_error(
            line,
            field,
            format!("expected exactly {length} lowercase hexadecimal characters"),
        ));
    }
    Ok(())
}

fn validate_stable_token(
    value: &str,
    line: usize,
    field: &'static str,
) -> Result<(), CorpusManifestError> {
    if value.is_empty()
        || value.len() > 96
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
        || value.starts_with('-')
        || value.ends_with('-')
    {
        return Err(manifest_error(
            line,
            field,
            "expected a stable lowercase ASCII token",
        ));
    }
    Ok(())
}

fn validate_date(value: &str, line: usize) -> Result<(), CorpusManifestError> {
    let bytes = value.as_bytes();
    if bytes.len() != 10
        || bytes[4] != b'-'
        || bytes[7] != b'-'
        || bytes
            .iter()
            .enumerate()
            .any(|(index, byte)| index != 4 && index != 7 && !byte.is_ascii_digit())
    {
        return Err(manifest_error(line, "reviewed_at", "expected YYYY-MM-DD"));
    }
    let year = value[0..4].parse::<u16>().unwrap_or(0);
    let month = value[5..7].parse::<u8>().unwrap_or(0);
    let day = value[8..10].parse::<u8>().unwrap_or(0);
    let leap_year =
        year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400));
    let days_in_month = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap_year => 29,
        2 => 28,
        _ => 0,
    };
    if year == 0 || day == 0 || day > days_in_month {
        return Err(manifest_error(
            line,
            "reviewed_at",
            "expected a calendar-valid YYYY-MM-DD date",
        ));
    }
    Ok(())
}

fn manifest_error(
    line: usize,
    field: &'static str,
    reason: impl Into<String>,
) -> CorpusManifestError {
    CorpusManifestError {
        line,
        field,
        reason: reason.into(),
    }
}

fn checkpoint(cx: &Cx<'_>, case_id: &str, stage: &'static str) -> Result<(), CorpusRunError> {
    cx.checkpoint().map_err(|_| CorpusRunError::Cancelled {
        case_id: case_id.to_string(),
        stage,
    })
}

fn rate_json(count: usize, total: usize) -> String {
    let basis_points = if total == 0 {
        0
    } else {
        count.saturating_mul(10_000) / total
    };
    format!("{{\"count\":{count},\"total\":{total},\"basis_points\":{basis_points}}}")
}

fn append_row_json(output: &mut String, row: &CorpusScorecardRow) {
    let (review_state, reviewer, reviewed_at, revision) = match &row.annotation_authority {
        CorpusAnnotationAuthority::Proposed => ("proposed", "", "", 0),
        CorpusAnnotationAuthority::HumanLocked {
            reviewer,
            reviewed_at,
            revision,
        } => (
            "human-locked",
            reviewer.as_str(),
            reviewed_at.as_str(),
            *revision,
        ),
    };
    let mut findings = String::from("[");
    for (index, finding) in row.observed.findings.iter().enumerate() {
        if index > 0 {
            findings.push(',');
        }
        let _ = write!(
            findings,
            "{{\"class\":\"{}\",\"count\":{}}}",
            json_escape(&finding.class),
            finding.count
        );
    }
    findings.push(']');
    let refusal_code = row
        .observed
        .refusal_code
        .as_ref()
        .map_or_else(|| "null".to_string(), |code| json_string(code));
    let proposed_matches = match &row.agreement {
        AnnotationAgreement::Unreviewed {
            proposed_matches, ..
        } => proposed_matches.to_string(),
        _ => "null".to_string(),
    };
    let _ = write!(
        output,
        "{{\"case_id\":\"{}\",\"relative_path\":\"{}\",\"content_blake3\":\"{}\",\
         \"format\":\"{}\",\"source\":{{\"kind\":\"{}\",\"origin\":\"{}\",\
         \"revision\":\"{}\",\"path\":\"{}\",\"object_identity\":\"{}\",\
         \"license_spdx\":\"{}\",\"license_url\":\"{}\"}},\
         \"quality_tier\":\"{}\",\
         \"expected\":{{\"outcome\":\"{}\",\"detail\":\"{}\"}},\
         \"annotation\":{{\"state\":\"{}\",\"reviewer\":\"{}\",\"reviewed_at\":\"{}\",\
         \"revision\":{}}},\"observed\":{{\"verdict\":\"{}\",\"findings\":{},\
         \"refusal_code\":{},\"diagnostic\":\"{}\",\"receipt\":{}}},\
         \"agreement\":{{\"status\":\"{}\",\"proposed_matches\":{},\"reason\":\"{}\"}}}}",
        json_escape(&row.case_id),
        json_escape(&row.relative_path),
        row.content_blake3,
        row.format.as_str(),
        row.source_kind.as_str(),
        json_escape(&row.source_origin),
        json_escape(&row.source_revision),
        json_escape(&row.source_path),
        json_escape(&row.source_object_identity),
        json_escape(&row.license_spdx),
        json_escape(&row.license_url),
        row.quality_tier.as_str(),
        row.expected.kind(),
        json_escape(&row.expected.detail()),
        review_state,
        json_escape(reviewer),
        json_escape(reviewed_at),
        revision,
        row.observed.verdict.as_str(),
        findings,
        refusal_code,
        json_escape(&row.observed.diagnostic),
        row.observed.receipt_json,
        row.agreement.as_str(),
        proposed_matches,
        json_escape(row.agreement.reason())
    );
}

fn json_string(value: &str) -> String {
    format!("\"{}\"", json_escape(value))
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
            control if control <= '\u{1f}' => {
                let _ = write!(output, "\\u{:04x}", control as u32);
            }
            character => output.push(character),
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x11_06,
                    kernel_id: 11,
                    tile: 6,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn locked_case(bytes: &[u8], expected: CorpusExpectedOutcome) -> SupplierCorpusCase {
        SupplierCorpusCase {
            case_id: "seeded-parser-fault".to_string(),
            relative_path: "sources/seeded-parser-fault.stl".to_string(),
            format: SupplierCadFormat::Stl,
            quality_tier: CorpusQualityTier::KnownBroken,
            source_kind: CorpusSourceKind::Git,
            source_origin: "https://example.invalid/source".to_string(),
            source_revision: "0".repeat(40),
            source_path: "source.stl".to_string(),
            source_object_identity: "0".repeat(40),
            license_spdx: "LicenseRef-test".to_string(),
            license_url: "https://example.invalid/license".to_string(),
            content_blake3: hash_bytes(bytes),
            expected,
            annotation_authority: CorpusAnnotationAuthority::HumanLocked {
                reviewer: "test-reviewer".to_string(),
                reviewed_at: "2026-07-24".to_string(),
                revision: 1,
            },
            justification: "fault-injection fixture".to_string(),
        }
    }

    fn manifest_row(
        case_id: &str,
        relative_path: &str,
        content_blake3: &str,
        review_state: &str,
        reviewer: &str,
        reviewed_at: &str,
        revision: &str,
    ) -> String {
        [
            case_id.to_string(),
            relative_path.to_string(),
            "stl".to_string(),
            "known-broken".to_string(),
            "git".to_string(),
            "https://example.invalid/repo".to_string(),
            "0".repeat(40),
            format!("{case_id}.stl"),
            "1".repeat(40),
            "MIT".to_string(),
            "https://example.invalid/license".to_string(),
            content_blake3.to_string(),
            "refuse".to_string(),
            "stl-parse-malformed".to_string(),
            review_state.to_string(),
            reviewer.to_string(),
            reviewed_at.to_string(),
            revision.to_string(),
            "review fixture".to_string(),
        ]
        .join("\t")
    }

    #[test]
    fn g0_manifest_refuses_unlocked_rows_that_claim_review_metadata() {
        let row = manifest_row(
            "case-a",
            "sources/case-a.stl",
            &"2".repeat(64),
            "proposed",
            "not-a-reviewer",
            "",
            "0",
        );
        let source = format!("{SUPPLIER_CORPUS_MANIFEST_HEADER}\n{row}\n");
        let error = parse_supplier_corpus_manifest(&source)
            .expect_err("proposed row cannot carry review authority");
        assert_eq!(error.field, "review_state");
    }

    #[test]
    fn g0_manifest_refuses_duplicate_retained_byte_identities() {
        let content_blake3 = "2".repeat(64);
        let first = manifest_row(
            "case-a",
            "sources/case-a.stl",
            &content_blake3,
            "proposed",
            "",
            "",
            "0",
        );
        let second = manifest_row(
            "case-b",
            "sources/case-b.stl",
            &content_blake3,
            "proposed",
            "",
            "",
            "0",
        );
        let source = format!("{SUPPLIER_CORPUS_MANIFEST_HEADER}\n{first}\n{second}\n");
        let error = parse_supplier_corpus_manifest(&source)
            .expect_err("one retained byte sequence cannot inflate the population");
        assert_eq!(error.field, "content_blake3");
    }

    #[test]
    fn g0_manifest_refuses_impossible_human_review_dates() {
        let row = manifest_row(
            "case-a",
            "sources/case-a.stl",
            &"2".repeat(64),
            "human-locked",
            "identified reviewer",
            "2026-13-40",
            "1",
        );
        let source = format!("{SUPPLIER_CORPUS_MANIFEST_HEADER}\n{row}\n");
        let error = parse_supplier_corpus_manifest(&source)
            .expect_err("human review dates must be calendar-valid");
        assert_eq!(error.field, "reviewed_at");
    }

    #[test]
    fn g0_http_snapshot_requires_sha256_object_identity() {
        assert!(
            validate_source_provenance(
                CorpusSourceKind::HttpSnapshot,
                "3d-package:c9105426-6818-4c25-b04c-135e79203b20@2020-04-16",
                &"2".repeat(64),
                1,
            )
            .is_ok()
        );
        let error = validate_source_provenance(
            CorpusSourceKind::HttpSnapshot,
            "3d-package:c9105426-6818-4c25-b04c-135e79203b20@2020-04-16",
            &"2".repeat(40),
            1,
        )
        .expect_err("HTTP snapshots must not disguise Git SHA-1 as object identity");
        assert_eq!(error.field, "source_object_identity");
    }

    #[test]
    fn g0_seeded_parser_fault_becomes_locked_annotation_mismatch() {
        let corrupt = b"this is not an STL resource";
        let case = locked_case(corrupt, CorpusExpectedOutcome::ImportClean);
        let policy = SupplierCorpusPolicy::try_standing_lane().expect("standing policy");
        let observed = with_cx(|cx| evaluate_case(&case, corrupt, policy, cx))
            .expect("scorecard run is not cancelled");
        assert_eq!(observed.verdict, CorpusVerdict::Refused);
        assert_eq!(
            observed.refusal_code.as_deref(),
            Some("stl-parse-malformed")
        );
        assert!(matches!(
            compare_annotation(&case, &observed),
            AnnotationAgreement::Mismatch { .. }
        ));
    }

    #[test]
    fn g0_program_metrics_exclude_proposed_annotations_from_rate_counts() {
        let bytes = b"this is not an STL resource";
        let expected = CorpusExpectedOutcome::Refuse {
            code: "stl-parse-malformed".to_string(),
        };
        let locked = locked_case(bytes, expected.clone());
        let mut proposed = locked_case(bytes, expected);
        proposed.case_id = "proposed-parser-fault".to_string();
        proposed.relative_path = "sources/proposed-parser-fault.stl".to_string();
        proposed.annotation_authority = CorpusAnnotationAuthority::Proposed;
        let manifest = CorpusManifest {
            manifest_identity: hash_bytes(b"authority-aware-metrics"),
            cases: vec![locked, proposed],
        };
        let scorecard = with_cx(|cx| {
            run_supplier_corpus(
                &manifest,
                SupplierCorpusPolicy::try_standing_lane().expect("standing policy"),
                cx,
                |_| Ok(bytes.to_vec()),
            )
        })
        .expect("scorecard run is not cancelled");

        assert_eq!(scorecard.refused, 2, "the diagnostic population sees both");
        let metrics = scorecard.import_metrics();
        assert_eq!(metrics.total(), 2);
        assert_eq!(metrics.reviewed(), 1);
        assert_eq!(metrics.clean(), 0);
        assert_eq!(metrics.repaired(), 0);
        assert_eq!(metrics.refused(), 1);
        assert_eq!(metrics.annotation_mismatches(), 0);
        let summary = scorecard.summary_json();
        assert!(summary.contains("\"population\":{\"total\":2"));
        assert!(summary.contains("\"reviewed\":{\"total\":1"));
        assert!(summary.contains("\"authority\":\"human-locked-only-dashboard-denominator\""));
    }

    #[test]
    fn g4_cancellation_refuses_partial_scorecard_publication() {
        let bytes = b"not needed because cancellation wins";
        let case = locked_case(
            bytes,
            CorpusExpectedOutcome::Refuse {
                code: "stl-parse-malformed".to_string(),
            },
        );
        let manifest = CorpusManifest {
            manifest_identity: hash_bytes(b"cancelled-manifest"),
            cases: vec![case],
        };
        let gate = CancelGate::new_clock_free();
        gate.request();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x11_06,
                    kernel_id: 11,
                    tile: 6,
                    iteration: 1,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let error = run_supplier_corpus(
                &manifest,
                SupplierCorpusPolicy::try_standing_lane().expect("standing policy"),
                &cx,
                |_| Ok(bytes.to_vec()),
            )
            .expect_err("cancelled run must publish no scorecard");
            assert!(matches!(error, CorpusRunError::Cancelled { .. }));
        });
    }
}

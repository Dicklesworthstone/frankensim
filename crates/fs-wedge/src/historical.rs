//! Offline replay for the exact source snapshot behind the measured comparison.
//!
//! Live readiness probes intentionally inspect the current workspace. Historical
//! comparison factors do not: they replay against a retained archive whose
//! revision, manifest, source bytes, and pointer set are all content-bound.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{DomainHasher, hash_domain};

use crate::{
    ComparisonCandidate, DEFAULT_FACTOR_WEIGHTS, EvidenceKind, FactorWeight, ScoringFactor,
};

/// Schema carried by the retained comparison-evidence snapshot descriptor.
pub const COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA: &str =
    "frankensim-wedge-comparison-evidence-snapshot-v3";
/// Parser and admission policy version for the retained source bundle.
pub const COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION: u32 = 2;
/// Domain for the complete retained-evidence descriptor identity.
pub const COMPARISON_EVIDENCE_DESCRIPTOR_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-descriptor.v1";
/// Domain for the exact manifest-byte identity.
pub const COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-manifest.v1";
/// Domain for the exact retained TAR-byte identity.
pub const COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-bundle.v1";
/// Domain for one revision/path/source-byte identity.
pub const COMPARISON_EVIDENCE_SOURCE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-source.v1";
/// Domain for the canonical default weights and complete comparison model.
pub const COMPARISON_MODEL_IDENTITY_DOMAIN: &str = "frankensim.fs-wedge.comparison-model.v1";

/// Maximum retained source-bundle bytes accepted by the replay adapter.
pub const MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES: usize = 1024 * 1024;
/// Maximum manifest bytes accepted by the replay adapter.
pub const MAX_COMPARISON_EVIDENCE_MANIFEST_BYTES: usize = 128 * 1024;
/// Maximum source or pointer rows accepted from one manifest.
pub const MAX_COMPARISON_EVIDENCE_ROWS: usize = 512;
/// Maximum UTF-8 bytes accepted in one manifest scalar.
pub const MAX_COMPARISON_EVIDENCE_FIELD_BYTES: usize = 4096;
/// Maximum conservative KMP comparison-work units across all locator scans.
///
/// One admitted pointer is charged `2 * (source_bytes + locator_bytes)`,
/// where one work unit is one needle/source byte equality in the
/// single-comparison KMP transition. This covers the standard `2N` prefix-table
/// and `2H` source-scan upper bounds. The complete charge is checked before the
/// first marker search begins.
pub const MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK: u64 = 16 * 1024 * 1024;

const TAR_BLOCK_BYTES: usize = 512;
const MANIFEST_BYTES: &str = include_str!("../data/comparison-evidence-b3b5f2c1.tsv");
const BUNDLE_BYTES: &[u8] = include_bytes!("../data/comparison-evidence-b3b5f2c1.tar");

#[cfg(test)]
std::thread_local! {
    static MARKER_SEARCH_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Immutable adapter metadata for one retained comparison source snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HistoricalEvidenceSnapshot {
    /// Stable schema label.
    pub schema: &'static str,
    /// Replay-policy version.
    pub policy_version: u32,
    /// Git revision declared by the retained descriptor and candidate rows.
    ///
    /// Content replay checks consistency with this label; it does not prove
    /// that Git produced the bytes or that an authorized human reviewed them.
    pub declared_inventory_revision: &'static str,
    /// Domain-separated BLAKE3 of every other descriptor field.
    ///
    /// This authenticates descriptor labels, including both artifact paths;
    /// it does not prove that either path exists or names the artifact's
    /// filesystem origin.
    pub descriptor_identity_blake3: &'static str,
    /// Workspace-relative retained-manifest path.
    pub manifest_path: &'static str,
    /// Domain-separated BLAKE3 of the exact manifest bytes.
    pub manifest_identity_blake3: &'static str,
    /// Workspace-relative retained-bundle path.
    pub bundle_path: &'static str,
    /// Domain-separated BLAKE3 of the exact TAR bytes.
    pub bundle_identity_blake3: &'static str,
    /// Domain-separated BLAKE3 of the complete comparison model and weights.
    pub comparison_model_identity_blake3: &'static str,
    /// Required regular-file count.
    pub source_count: usize,
    /// Required historical `WorkspacePath` occurrence count.
    pub pointer_count: usize,
}

impl HistoricalEvidenceSnapshot {
    /// Is every descriptor field structurally admissible?
    #[must_use]
    pub fn is_complete(self) -> bool {
        self.schema == COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA
            && self.policy_version == COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION
            && is_lower_hex(self.declared_inventory_revision, 40)
            && is_lower_hex(self.descriptor_identity_blake3, 64)
            && self
                .descriptor_identity_blake3
                .bytes()
                .any(|byte| byte != b'0')
            && valid_relative_path(self.manifest_path)
            && is_lower_hex(self.manifest_identity_blake3, 64)
            && valid_relative_path(self.bundle_path)
            && is_lower_hex(self.bundle_identity_blake3, 64)
            && is_lower_hex(self.comparison_model_identity_blake3, 64)
            && self
                .comparison_model_identity_blake3
                .bytes()
                .any(|byte| byte != b'0')
            && self.source_count > 0
            && self.source_count <= MAX_COMPARISON_EVIDENCE_ROWS
            && self.pointer_count > 0
            && self.pointer_count <= MAX_COMPARISON_EVIDENCE_ROWS
    }
}

/// The source snapshot bound to the default measured comparison.
///
/// The four BLAKE3 identities are intentionally literal protocol fields. Tests
/// recompute them from the complete descriptor, embedded artifacts, and
/// comparison model before any default recommendation can be consumed.
pub const COMPARISON_EVIDENCE_SNAPSHOT: HistoricalEvidenceSnapshot = HistoricalEvidenceSnapshot {
    schema: COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA,
    policy_version: COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION,
    declared_inventory_revision: "b3b5f2c1c809eec06cde1e40cbc916d6995469b5",
    descriptor_identity_blake3: "87778e8008f2d6ad5b898828dcede63cd2c30b94b4a57332e30446e88589d6e8",
    manifest_path: "crates/fs-wedge/data/comparison-evidence-b3b5f2c1.tsv",
    manifest_identity_blake3: "2ea962dc43b416b22e1591e23b065d20a2c93be4b1440aeff7228e9c37ace6ec",
    bundle_path: "crates/fs-wedge/data/comparison-evidence-b3b5f2c1.tar",
    bundle_identity_blake3: "4a13f162b435a126c979c1da6a8743e4c11388648e6e225ebf577da663f0e9e3",
    comparison_model_identity_blake3: "d9821a2b7cc45bfd635991fc1886e8577cb1138af8ff13c5469693c905cb8b54",
    source_count: 13,
    pointer_count: 31,
};

/// How a successful replay obtained its protocol inputs.
///
/// This is a trust-origin label, not an authority upgrade. Both variants prove
/// at most content integrity and protocol consistency and explicitly carry no
/// current-decision or human-review authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistoricalEvidenceTrustOrigin {
    /// The crate's private default constructor used the embedded descriptor,
    /// manifest, bundle, weights, and comparison model.
    EmbeddedDefault,
    /// A public protocol verifier used caller-supplied descriptor, artifacts,
    /// or comparison records, even if their bytes equal the embedded default.
    CallerSuppliedProtocolConsistency,
}

impl HistoricalEvidenceTrustOrigin {
    /// Stable machine-readable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::EmbeddedDefault => "embedded-default",
            Self::CallerSuppliedProtocolConsistency => "caller-supplied-protocol-consistency",
        }
    }
}

/// Successful offline replay of a comparison source snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HistoricalEvidenceReceipt {
    schema: &'static str,
    policy_version: u32,
    declared_inventory_revision: &'static str,
    descriptor_identity_blake3: String,
    manifest_path: &'static str,
    manifest_identity_blake3: String,
    bundle_path: &'static str,
    bundle_identity_blake3: String,
    comparison_model_identity_blake3: String,
    source_count: usize,
    pointer_count: usize,
    trust_origin: HistoricalEvidenceTrustOrigin,
    current_decision_authority: bool,
    human_review_authority: bool,
}

impl HistoricalEvidenceReceipt {
    /// Snapshot schema.
    #[must_use]
    pub const fn schema(&self) -> &'static str {
        self.schema
    }

    /// Replay policy version.
    #[must_use]
    pub const fn policy_version(&self) -> u32 {
        self.policy_version
    }

    /// Descriptor-declared inventory revision.
    ///
    /// This is a consistency label, not proof of Git provenance or review.
    #[must_use]
    pub const fn declared_inventory_revision(&self) -> &'static str {
        self.declared_inventory_revision
    }

    /// Canonical identity of every snapshot field except the identity itself.
    #[must_use]
    pub fn descriptor_identity_blake3(&self) -> &str {
        &self.descriptor_identity_blake3
    }

    /// Descriptor-authenticated manifest path label.
    ///
    /// This label does not prove filesystem existence or origin.
    #[must_use]
    pub const fn manifest_path(&self) -> &'static str {
        self.manifest_path
    }

    /// Exact domain-separated manifest identity.
    #[must_use]
    pub fn manifest_identity_blake3(&self) -> &str {
        &self.manifest_identity_blake3
    }

    /// Descriptor-authenticated bundle path label.
    ///
    /// This label does not prove filesystem existence or origin.
    #[must_use]
    pub const fn bundle_path(&self) -> &'static str {
        self.bundle_path
    }

    /// Exact domain-separated source-bundle identity.
    #[must_use]
    pub fn bundle_identity_blake3(&self) -> &str {
        &self.bundle_identity_blake3
    }

    /// Canonical identity of the complete comparison model and default weights.
    #[must_use]
    pub fn comparison_model_identity_blake3(&self) -> &str {
        &self.comparison_model_identity_blake3
    }

    /// Number of content-bound source files.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }

    /// Number of replayed historical pointer occurrences.
    #[must_use]
    pub const fn pointer_count(&self) -> usize {
        self.pointer_count
    }

    /// Machine-readable origin of the replayed inputs.
    #[must_use]
    pub const fn trust_origin(&self) -> HistoricalEvidenceTrustOrigin {
        self.trust_origin
    }

    /// Whether this receipt authorizes a current recommendation.
    ///
    /// Historical replay never carries this authority.
    #[must_use]
    pub const fn current_decision_authority(&self) -> bool {
        self.current_decision_authority
    }

    /// Whether this receipt proves authorized human review.
    ///
    /// Historical replay never carries this authority.
    #[must_use]
    pub const fn human_review_authority(&self) -> bool {
        self.human_review_authority
    }
}

/// Fail-closed historical source replay diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HistoricalEvidenceError {
    /// Static snapshot metadata is incomplete or malformed.
    InvalidDescriptor {
        /// Exact invalid field.
        field: &'static str,
        /// Actionable detail.
        detail: String,
    },
    /// The complete descriptor does not match its declared identity.
    DescriptorIdentityMismatch {
        /// Snapshot-declared domain-separated BLAKE3.
        expected: String,
        /// Recomputed domain-separated BLAKE3.
        observed: String,
    },
    /// An embedded artifact exceeds its admitted byte envelope.
    ArtifactTooLarge {
        /// `manifest` or `bundle`.
        artifact: &'static str,
        /// Observed bytes.
        observed: usize,
        /// Maximum admitted bytes.
        maximum: usize,
    },
    /// Exact artifact bytes do not match the declared content identity.
    IdentityMismatch {
        /// `manifest` or `bundle`.
        artifact: &'static str,
        /// Declared domain-separated BLAKE3.
        expected: String,
        /// Recomputed domain-separated BLAKE3.
        observed: String,
    },
    /// The complete comparison model or default weights drifted.
    ComparisonModelIdentityMismatch {
        /// Snapshot-declared domain-separated BLAKE3.
        expected: String,
        /// Recomputed domain-separated BLAKE3.
        observed: String,
    },
    /// The supplied comparison model is structurally outside the protocol.
    InvalidComparisonModel {
        /// Canonical field path.
        field: String,
        /// Actionable detail.
        detail: String,
    },
    /// The canonical TSV manifest is malformed.
    MalformedManifest {
        /// One-based line number, or zero for whole-manifest failures.
        line: usize,
        /// Actionable detail.
        detail: String,
    },
    /// A candidate declares a revision other than the snapshot descriptor.
    CandidateRevisionMismatch {
        /// Candidate slug.
        candidate: &'static str,
        /// Snapshot revision.
        expected: &'static str,
        /// Candidate revision.
        observed: &'static str,
    },
    /// Manifest and comparison pointer sequences differ.
    PointerMismatch {
        /// Descriptor-declared snapshot revision.
        declared_inventory_revision: &'static str,
        /// Authenticated manifest identity.
        manifest_identity_blake3: &'static str,
        /// Zero-based canonical pointer index.
        index: usize,
        /// Expected candidate/factor/path/locator.
        expected: String,
        /// Manifest candidate/factor/path/locator, or `<missing>`.
        observed: String,
    },
    /// A caller-supplied comparison pointer is outside the admitted protocol.
    InvalidCandidatePointer {
        /// Candidate slug.
        candidate: &'static str,
        /// Factor label.
        factor: &'static str,
        /// Workspace-relative path.
        reference: &'static str,
        /// Locator substring.
        locator: &'static str,
        /// Actionable detail.
        detail: String,
    },
    /// Caller-supplied candidate evidence exceeds the admitted row envelope.
    CandidateEvidenceTooLarge {
        /// Observed lower bound when iteration stopped.
        observed: usize,
        /// Maximum admitted pointer rows.
        maximum: usize,
    },
    /// Aggregate marker-search work exceeds the synchronous replay budget.
    MarkerScanBudgetExceeded {
        /// Conservative KMP work units required by the complete pointer table,
        /// saturated to `u64::MAX` if arithmetic cannot represent the bound.
        observed: u64,
        /// Maximum admitted conservative KMP work units.
        maximum: u64,
    },
    /// The manifest retains a source not consumed by any comparison pointer.
    UnexpectedSource {
        /// Unconsumed workspace-relative source path.
        reference: String,
    },
    /// The retained TAR container is malformed.
    MalformedBundle {
        /// Byte offset of the failing TAR header or payload.
        offset: usize,
        /// Actionable detail.
        detail: String,
    },
    /// The TAR does not contain exactly the manifest-declared file.
    SourceMismatch {
        /// Descriptor-declared snapshot revision.
        declared_inventory_revision: &'static str,
        /// Authenticated whole-bundle identity.
        bundle_identity_blake3: &'static str,
        /// Workspace-relative source path.
        reference: String,
        /// Expected byte count, when the source exists in the manifest.
        expected_bytes: Option<usize>,
        /// Observed byte count, or `None` when absent.
        observed_bytes: Option<usize>,
    },
    /// One retained file does not match its revision/path/source identity.
    SourceIdentityMismatch {
        /// Descriptor-declared snapshot revision.
        declared_inventory_revision: &'static str,
        /// Candidate owning the first canonical pointer into this source.
        candidate: Box<str>,
        /// Factor owning the first canonical pointer into this source.
        factor: Box<str>,
        /// Workspace-relative source path.
        reference: Box<str>,
        /// Locator carried by the first canonical pointer into this source.
        locator: Box<str>,
        /// Manifest-declared domain-separated source identity.
        expected_blake3: Box<str>,
        /// Recomputed domain-separated source identity.
        observed_blake3: Box<str>,
    },
    /// Authenticated historical bytes do not contain their recorded locator.
    MarkerMissing {
        /// Descriptor-declared snapshot revision.
        declared_inventory_revision: &'static str,
        /// Verified domain-separated identity of the exact source bytes.
        source_identity_blake3: Box<str>,
        /// Candidate slug.
        candidate: Box<str>,
        /// Factor label.
        factor: Box<str>,
        /// Workspace-relative source path.
        reference: Box<str>,
        /// Missing marker.
        locator: Box<str>,
    },
}

impl fmt::Display for HistoricalEvidenceError {
    #[allow(clippy::too_many_lines)] // One exhaustive table keeps every typed refusal's wording adjacent.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDescriptor { field, detail } => {
                write!(f, "invalid historical snapshot {field}: {detail}")
            }
            Self::DescriptorIdentityMismatch { expected, observed } => write!(
                f,
                "historical descriptor identity mismatch: expected {expected}, observed {observed}"
            ),
            Self::ArtifactTooLarge {
                artifact,
                observed,
                maximum,
            } => write!(
                f,
                "historical {artifact} is {observed} bytes; maximum is {maximum}"
            ),
            Self::IdentityMismatch {
                artifact,
                expected,
                observed,
            } => write!(
                f,
                "historical {artifact} identity mismatch: expected {expected}, observed {observed}"
            ),
            Self::ComparisonModelIdentityMismatch { expected, observed } => write!(
                f,
                "historical comparison-model identity mismatch: expected {expected}, observed {observed}"
            ),
            Self::InvalidComparisonModel { field, detail } => {
                write!(f, "invalid historical comparison model {field}: {detail}")
            }
            Self::MalformedManifest { line, detail } => {
                write!(f, "malformed historical manifest at line {line}: {detail}")
            }
            Self::CandidateRevisionMismatch {
                candidate,
                expected,
                observed,
            } => write!(
                f,
                "candidate {candidate} declares inventory revision {observed}; snapshot descriptor declares {expected}"
            ),
            Self::PointerMismatch {
                declared_inventory_revision,
                manifest_identity_blake3,
                index,
                expected,
                observed,
            } => write!(
                f,
                "historical pointer {index} mismatch at declared revision {declared_inventory_revision} \
                 (manifest BLAKE3 {manifest_identity_blake3}): expected {expected}; observed {observed}"
            ),
            Self::InvalidCandidatePointer {
                candidate,
                factor,
                reference,
                locator,
                detail,
            } => write!(
                f,
                "invalid historical pointer {candidate}/{factor}/{reference}/{locator}: {detail}"
            ),
            Self::CandidateEvidenceTooLarge { observed, maximum } => write!(
                f,
                "historical candidate evidence has at least {observed} pointer rows; maximum is {maximum}"
            ),
            Self::MarkerScanBudgetExceeded { observed, maximum } => write!(
                f,
                "historical marker scan requires {observed} conservative KMP work units; maximum is {maximum}"
            ),
            Self::UnexpectedSource { reference } => {
                write!(
                    f,
                    "historical manifest retains unreferenced source {reference}"
                )
            }
            Self::MalformedBundle { offset, detail } => {
                write!(f, "malformed historical TAR at byte {offset}: {detail}")
            }
            Self::SourceMismatch {
                declared_inventory_revision,
                bundle_identity_blake3,
                reference,
                expected_bytes,
                observed_bytes,
            } => write!(
                f,
                "historical source {reference} set/length mismatch at declared revision \
                 {declared_inventory_revision} (bundle BLAKE3 {bundle_identity_blake3}): \
                 expected {expected_bytes:?} bytes, observed {observed_bytes:?}"
            ),
            Self::SourceIdentityMismatch {
                declared_inventory_revision,
                candidate,
                factor,
                reference,
                locator,
                expected_blake3,
                observed_blake3,
            } => write!(
                f,
                "historical source identity mismatch at declared revision {declared_inventory_revision} for \
                 {candidate}/{factor}/{reference}/{locator}: expected BLAKE3 \
                 {expected_blake3}, observed {observed_blake3}"
            ),
            Self::MarkerMissing {
                declared_inventory_revision,
                source_identity_blake3,
                candidate,
                factor,
                reference,
                locator,
            } => write!(
                f,
                "historical candidate {candidate} factor {factor} lost marker {locator:?} in \
                 {reference} at declared revision {declared_inventory_revision} after verifying source BLAKE3 \
                 {source_identity_blake3}"
            ),
        }
    }
}

impl std::error::Error for HistoricalEvidenceError {}

#[derive(Debug)]
struct SourceSpec<'a> {
    reference: &'a str,
    bytes: usize,
    identity_blake3: &'a str,
}

#[derive(Debug)]
struct PointerSpec<'a> {
    candidate: &'a str,
    factor: &'a str,
    reference: &'a str,
    locator: &'a str,
}

#[derive(Debug)]
struct ParsedManifest<'a> {
    revision: &'a str,
    sources: Vec<SourceSpec<'a>>,
    pointers: Vec<PointerSpec<'a>>,
}

/// Exact embedded manifest bytes used by the default comparison replay.
#[must_use]
pub const fn comparison_evidence_manifest() -> &'static str {
    MANIFEST_BYTES
}

/// Exact embedded TAR bytes used by the default comparison replay.
#[must_use]
pub const fn comparison_evidence_bundle() -> &'static [u8] {
    BUNDLE_BYTES
}

/// Canonical identity of every retained-evidence descriptor field except its
/// own identity.
///
/// The preimage is field-tagged, length-framed, and order-sensitive. It binds
/// the schema, policy version, declared revision, both artifact path labels,
/// both artifact identities, the comparison-model identity, and both admitted
/// row counts. The `descriptor_identity_blake3` field is deliberately excluded
/// to avoid a circular preimage. A valid root authenticates the supplied path
/// labels as protocol data; it does not prove filesystem existence or origin.
#[must_use]
pub fn comparison_evidence_descriptor_identity_blake3(
    snapshot: HistoricalEvidenceSnapshot,
) -> String {
    let mut hasher = DomainHasher::new(COMPARISON_EVIDENCE_DESCRIPTOR_IDENTITY_DOMAIN);
    update_framed_str(&mut hasher, "schema");
    update_framed_str(&mut hasher, snapshot.schema);
    update_framed_str(&mut hasher, "policy_version");
    hasher.update(&snapshot.policy_version.to_le_bytes());
    update_framed_str(&mut hasher, "declared_inventory_revision");
    update_framed_str(&mut hasher, snapshot.declared_inventory_revision);
    update_framed_str(&mut hasher, "manifest_path");
    update_framed_str(&mut hasher, snapshot.manifest_path);
    update_framed_str(&mut hasher, "manifest_identity_blake3");
    update_framed_str(&mut hasher, snapshot.manifest_identity_blake3);
    update_framed_str(&mut hasher, "bundle_path");
    update_framed_str(&mut hasher, snapshot.bundle_path);
    update_framed_str(&mut hasher, "bundle_identity_blake3");
    update_framed_str(&mut hasher, snapshot.bundle_identity_blake3);
    update_framed_str(&mut hasher, "comparison_model_identity_blake3");
    update_framed_str(&mut hasher, snapshot.comparison_model_identity_blake3);
    update_framed_str(&mut hasher, "source_count");
    update_framed_count(&mut hasher, snapshot.source_count);
    update_framed_str(&mut hasher, "pointer_count");
    update_framed_count(&mut hasher, snapshot.pointer_count);
    hasher.finalize().to_hex()
}

/// Canonical identity of supplied weights and the complete comparison model.
///
/// The preimage is length-framed and order-sensitive. It binds every weight;
/// candidate identity/date/revision/minority case; factor/rating; measurement
/// readiness, score, method, finding, and rationale; and every evidence
/// pointer kind, reference, and locator.
pub fn comparison_model_identity_blake3(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
) -> Result<String, HistoricalEvidenceError> {
    validate_comparison_model(weights, candidates)?;
    let mut hasher = DomainHasher::new(COMPARISON_MODEL_IDENTITY_DOMAIN);
    update_framed_str(&mut hasher, "weights");
    update_framed_count(&mut hasher, weights.len());
    for weight in weights {
        update_framed_str(&mut hasher, weight.factor.label());
        hasher.update(&[weight.weight]);
    }
    update_framed_str(&mut hasher, "candidates");
    update_framed_count(&mut hasher, candidates.len());
    for candidate in candidates {
        update_framed_str(&mut hasher, candidate.name);
        update_framed_str(&mut hasher, candidate.display);
        update_framed_str(&mut hasher, candidate.measured_on);
        update_framed_str(&mut hasher, candidate.declared_inventory_revision);
        update_framed_count(&mut hasher, candidate.factors.len());
        for input in candidate.factors {
            update_framed_str(&mut hasher, input.factor.label());
            hasher.update(&[input.rating]);
            update_framed_str(&mut hasher, input.measurement.readiness.label());
            hasher.update(&[input.measurement.score]);
            update_framed_str(&mut hasher, input.measurement.method.label());
            update_framed_str(&mut hasher, input.measurement.finding);
            update_framed_str(&mut hasher, input.rationale);
            update_framed_count(&mut hasher, input.measurement.evidence.len());
            for pointer in input.measurement.evidence {
                update_framed_str(&mut hasher, pointer.kind.label());
                update_framed_str(&mut hasher, pointer.reference);
                update_framed_str(&mut hasher, pointer.locator);
            }
        }
        update_framed_str(&mut hasher, candidate.minority_case);
    }
    Ok(hasher.finalize().to_hex())
}

/// Replay the default comparison against its retained source bytes.
pub fn verify_default_comparison_evidence()
-> Result<HistoricalEvidenceReceipt, HistoricalEvidenceError> {
    verify_comparison_evidence_with_origin(
        COMPARISON_EVIDENCE_SNAPSHOT,
        MANIFEST_BYTES,
        BUNDLE_BYTES,
        crate::comparison_candidates(),
        HistoricalEvidenceTrustOrigin::EmbeddedDefault,
    )
}

/// Check caller-supplied comparison records against one retained source set.
///
/// This adapter reads no filesystem or Git state. Callers may supply alternate
/// bytes for mutation and migration tests, but the descriptor must bind those
/// bytes before any marker can be admitted. Even when every supplied byte
/// equals the embedded default, the returned trust origin remains
/// [`HistoricalEvidenceTrustOrigin::CallerSuppliedProtocolConsistency`].
pub fn verify_comparison_evidence(
    snapshot: HistoricalEvidenceSnapshot,
    manifest: &str,
    bundle: &[u8],
    candidates: &[ComparisonCandidate],
) -> Result<HistoricalEvidenceReceipt, HistoricalEvidenceError> {
    verify_comparison_evidence_with_origin(
        snapshot,
        manifest,
        bundle,
        candidates,
        HistoricalEvidenceTrustOrigin::CallerSuppliedProtocolConsistency,
    )
}

fn verify_comparison_evidence_with_origin(
    snapshot: HistoricalEvidenceSnapshot,
    manifest: &str,
    bundle: &[u8],
    candidates: &[ComparisonCandidate],
    trust_origin: HistoricalEvidenceTrustOrigin,
) -> Result<HistoricalEvidenceReceipt, HistoricalEvidenceError> {
    validate_descriptor(snapshot)?;
    let observed_descriptor = comparison_evidence_descriptor_identity_blake3(snapshot);
    if observed_descriptor != snapshot.descriptor_identity_blake3 {
        return Err(HistoricalEvidenceError::DescriptorIdentityMismatch {
            expected: snapshot.descriptor_identity_blake3.to_string(),
            observed: observed_descriptor,
        });
    }
    if manifest.len() > MAX_COMPARISON_EVIDENCE_MANIFEST_BYTES {
        return Err(HistoricalEvidenceError::ArtifactTooLarge {
            artifact: "manifest",
            observed: manifest.len(),
            maximum: MAX_COMPARISON_EVIDENCE_MANIFEST_BYTES,
        });
    }
    if bundle.len() > MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES {
        return Err(HistoricalEvidenceError::ArtifactTooLarge {
            artifact: "bundle",
            observed: bundle.len(),
            maximum: MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES,
        });
    }

    let observed_manifest = hash_domain(
        COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN,
        manifest.as_bytes(),
    )
    .to_hex();
    if observed_manifest != snapshot.manifest_identity_blake3 {
        return Err(HistoricalEvidenceError::IdentityMismatch {
            artifact: "manifest",
            expected: snapshot.manifest_identity_blake3.to_string(),
            observed: observed_manifest,
        });
    }
    let observed_bundle = hash_domain(COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN, bundle).to_hex();
    if observed_bundle != snapshot.bundle_identity_blake3 {
        return Err(HistoricalEvidenceError::IdentityMismatch {
            artifact: "bundle",
            expected: snapshot.bundle_identity_blake3.to_string(),
            observed: observed_bundle,
        });
    }

    let parsed = parse_manifest(manifest, snapshot)?;
    validate_pointer_table(&parsed, snapshot, candidates)?;
    let observed_model = comparison_model_identity_blake3(&DEFAULT_FACTOR_WEIGHTS, candidates)?;
    if observed_model != snapshot.comparison_model_identity_blake3 {
        return Err(HistoricalEvidenceError::ComparisonModelIdentityMismatch {
            expected: snapshot.comparison_model_identity_blake3.to_string(),
            observed: observed_model,
        });
    }
    let files = parse_tar_bundle(bundle, snapshot.declared_inventory_revision)?;
    validate_sources_and_markers(&parsed, &files, snapshot)?;

    Ok(HistoricalEvidenceReceipt {
        schema: snapshot.schema,
        policy_version: snapshot.policy_version,
        declared_inventory_revision: snapshot.declared_inventory_revision,
        descriptor_identity_blake3: observed_descriptor,
        manifest_path: snapshot.manifest_path,
        manifest_identity_blake3: observed_manifest,
        bundle_path: snapshot.bundle_path,
        bundle_identity_blake3: observed_bundle,
        comparison_model_identity_blake3: observed_model,
        source_count: parsed.sources.len(),
        pointer_count: parsed.pointers.len(),
        trust_origin,
        current_decision_authority: false,
        human_review_authority: false,
    })
}

fn validate_descriptor(
    snapshot: HistoricalEvidenceSnapshot,
) -> Result<(), HistoricalEvidenceError> {
    if snapshot.schema != COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA {
        return Err(HistoricalEvidenceError::InvalidDescriptor {
            field: "schema",
            detail: format!(
                "expected {COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA}, observed {}",
                snapshot.schema
            ),
        });
    }
    if snapshot.policy_version != COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION {
        return Err(HistoricalEvidenceError::InvalidDescriptor {
            field: "policy_version",
            detail: format!(
                "expected {COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION}, observed {}",
                snapshot.policy_version
            ),
        });
    }
    for (field, value, length) in [
        (
            "declared_inventory_revision",
            snapshot.declared_inventory_revision,
            40,
        ),
        (
            "descriptor_identity_blake3",
            snapshot.descriptor_identity_blake3,
            64,
        ),
        (
            "manifest_identity_blake3",
            snapshot.manifest_identity_blake3,
            64,
        ),
        (
            "bundle_identity_blake3",
            snapshot.bundle_identity_blake3,
            64,
        ),
        (
            "comparison_model_identity_blake3",
            snapshot.comparison_model_identity_blake3,
            64,
        ),
    ] {
        if !is_lower_hex(value, length) {
            return Err(HistoricalEvidenceError::InvalidDescriptor {
                field,
                detail: format!("must be exactly {length} lowercase hexadecimal characters"),
            });
        }
    }
    if snapshot
        .descriptor_identity_blake3
        .bytes()
        .all(|byte| byte == b'0')
    {
        return Err(HistoricalEvidenceError::InvalidDescriptor {
            field: "descriptor_identity_blake3",
            detail: "must be a pinned non-placeholder identity".to_string(),
        });
    }
    if snapshot
        .comparison_model_identity_blake3
        .bytes()
        .all(|byte| byte == b'0')
    {
        return Err(HistoricalEvidenceError::InvalidDescriptor {
            field: "comparison_model_identity_blake3",
            detail: "must be a pinned non-placeholder identity".to_string(),
        });
    }
    for (field, path) in [
        ("manifest_path", snapshot.manifest_path),
        ("bundle_path", snapshot.bundle_path),
    ] {
        if !valid_relative_path(path) {
            return Err(HistoricalEvidenceError::InvalidDescriptor {
                field,
                detail: "must be a normalized workspace-relative path".to_string(),
            });
        }
    }
    for (field, count) in [
        ("source_count", snapshot.source_count),
        ("pointer_count", snapshot.pointer_count),
    ] {
        if count == 0 || count > MAX_COMPARISON_EVIDENCE_ROWS {
            return Err(HistoricalEvidenceError::InvalidDescriptor {
                field,
                detail: format!("must be in 1..={MAX_COMPARISON_EVIDENCE_ROWS}"),
            });
        }
    }
    Ok(())
}

#[allow(clippy::too_many_lines)] // One ordered validation pass mirrors the canonical hash preimage.
fn validate_comparison_model(
    weights: &[FactorWeight],
    candidates: &[ComparisonCandidate],
) -> Result<(), HistoricalEvidenceError> {
    if weights.len() != ScoringFactor::ALL.len() {
        return Err(invalid_comparison_model(
            "weights.len",
            format!(
                "expected {}, observed {}",
                ScoringFactor::ALL.len(),
                weights.len()
            ),
        ));
    }
    let mut weight_sum = 0u16;
    for (index, (weight, expected_factor)) in weights.iter().zip(ScoringFactor::ALL).enumerate() {
        if weight.factor != expected_factor {
            return Err(invalid_comparison_model(
                format!("weights[{index}].factor"),
                format!(
                    "expected canonical factor {}, observed {}",
                    expected_factor.label(),
                    weight.factor.label()
                ),
            ));
        }
        weight_sum += u16::from(weight.weight);
    }
    if weight_sum != 100 {
        return Err(invalid_comparison_model(
            "weights.sum",
            format!("expected 100, observed {weight_sum}"),
        ));
    }
    if candidates.is_empty() || candidates.len() > MAX_COMPARISON_EVIDENCE_ROWS {
        return Err(invalid_comparison_model(
            "candidates.len",
            format!(
                "must be in 1..={MAX_COMPARISON_EVIDENCE_ROWS}; observed {}",
                candidates.len()
            ),
        ));
    }

    let mut candidate_names = BTreeSet::new();
    let mut evidence_count = 0usize;
    for (candidate_index, candidate) in candidates.iter().enumerate() {
        for (field, value) in [
            ("name", candidate.name),
            ("display", candidate.display),
            ("measured_on", candidate.measured_on),
            ("minority_case", candidate.minority_case),
        ] {
            if !valid_scalar(value) {
                return Err(invalid_comparison_model(
                    format!("candidates[{candidate_index}].{field}"),
                    "must be a non-empty bounded scalar without control characters",
                ));
            }
        }
        if !candidate_names.insert(candidate.name) {
            return Err(invalid_comparison_model(
                format!("candidates[{candidate_index}].name"),
                format!("duplicate candidate slug {}", candidate.name),
            ));
        }
        if !is_lower_hex(candidate.declared_inventory_revision, 40) {
            return Err(invalid_comparison_model(
                format!("candidates[{candidate_index}].declared_inventory_revision"),
                "must be exactly 40 lowercase hexadecimal characters",
            ));
        }
        if candidate.factors.len() != ScoringFactor::ALL.len() {
            return Err(invalid_comparison_model(
                format!("candidates[{candidate_index}].factors.len"),
                format!(
                    "expected {}, observed {}",
                    ScoringFactor::ALL.len(),
                    candidate.factors.len()
                ),
            ));
        }
        for (factor_index, (input, expected_factor)) in
            candidate.factors.iter().zip(ScoringFactor::ALL).enumerate()
        {
            let factor_path = format!("candidates[{candidate_index}].factors[{factor_index}]");
            if input.factor != expected_factor {
                return Err(invalid_comparison_model(
                    format!("{factor_path}.factor"),
                    format!(
                        "expected canonical factor {}, observed {}",
                        expected_factor.label(),
                        input.factor.label()
                    ),
                ));
            }
            if input.rating > 10 {
                return Err(invalid_comparison_model(
                    format!("{factor_path}.rating"),
                    format!("must be in 0..=10; observed {}", input.rating),
                ));
            }
            if input.measurement.score > input.measurement.readiness.score_ceiling() {
                return Err(invalid_comparison_model(
                    format!("{factor_path}.measurement.score"),
                    format!(
                        "score {} exceeds {} readiness ceiling {}",
                        input.measurement.score,
                        input.measurement.readiness.label(),
                        input.measurement.readiness.score_ceiling()
                    ),
                ));
            }
            for (field, value) in [
                ("measurement.finding", input.measurement.finding),
                ("rationale", input.rationale),
            ] {
                if !valid_scalar(value) {
                    return Err(invalid_comparison_model(
                        format!("{factor_path}.{field}"),
                        "must be a non-empty bounded scalar without control characters",
                    ));
                }
            }
            if input.measurement.evidence.is_empty() {
                return Err(invalid_comparison_model(
                    format!("{factor_path}.measurement.evidence"),
                    "must contain at least one pointer",
                ));
            }
            for (pointer_index, pointer) in input.measurement.evidence.iter().enumerate() {
                evidence_count = evidence_count.checked_add(1).ok_or_else(|| {
                    invalid_comparison_model(
                        "evidence.count",
                        "pointer count overflowed the admitted envelope",
                    )
                })?;
                if evidence_count > MAX_COMPARISON_EVIDENCE_ROWS {
                    return Err(invalid_comparison_model(
                        "evidence.count",
                        format!(
                            "maximum is {MAX_COMPARISON_EVIDENCE_ROWS}; observed at least {evidence_count}"
                        ),
                    ));
                }
                let pointer_path = format!("{factor_path}.measurement.evidence[{pointer_index}]");
                let valid_reference = match pointer.kind {
                    EvidenceKind::WorkspacePath => valid_relative_path(pointer.reference),
                    EvidenceKind::Bead | EvidenceKind::OfficialSource => {
                        valid_scalar(pointer.reference)
                    }
                };
                if !valid_reference {
                    return Err(invalid_comparison_model(
                        format!("{pointer_path}.reference"),
                        "is empty, oversized, contains controls, or violates its pointer-kind syntax",
                    ));
                }
                if !valid_scalar(pointer.locator) {
                    return Err(invalid_comparison_model(
                        format!("{pointer_path}.locator"),
                        "must be a non-empty bounded scalar without control characters",
                    ));
                }
            }
        }
    }
    Ok(())
}

fn invalid_comparison_model(
    field: impl Into<String>,
    detail: impl Into<String>,
) -> HistoricalEvidenceError {
    HistoricalEvidenceError::InvalidComparisonModel {
        field: field.into(),
        detail: detail.into(),
    }
}

fn update_framed_count(hasher: &mut DomainHasher, count: usize) {
    let count = u64::try_from(count).expect("canonical identity protocol counts fit in u64");
    hasher.update(&count.to_le_bytes());
}

fn update_framed_str(hasher: &mut DomainHasher, value: &str) {
    update_framed_count(hasher, value.len());
    hasher.update(value.as_bytes());
}

#[allow(clippy::too_many_lines)] // One ordered pass makes the canonical TSV refusal sequence auditable.
fn parse_manifest(
    manifest: &str,
    snapshot: HistoricalEvidenceSnapshot,
) -> Result<ParsedManifest<'_>, HistoricalEvidenceError> {
    if manifest.as_bytes().contains(&b'\r') {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 0,
            detail: "carriage returns are not canonical".to_string(),
        });
    }
    if !manifest.ends_with('\n') {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 0,
            detail: "manifest must end with one LF".to_string(),
        });
    }
    let mut lines = manifest.lines().enumerate();
    let Some((_, header)) = lines.next() else {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 0,
            detail: "manifest is empty".to_string(),
        });
    };
    if header != "FS-WEDGE-COMPARISON-EVIDENCE-SNAPSHOT\t1" {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 1,
            detail: format!("unexpected schema header {header:?}"),
        });
    }
    let Some((_, revision_line)) = lines.next() else {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 2,
            detail: "missing REVISION row".to_string(),
        });
    };
    let revision_fields: Vec<&str> = revision_line.split('\t').collect();
    if revision_fields.as_slice() != ["REVISION", snapshot.declared_inventory_revision] {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 2,
            detail: format!(
                "expected REVISION {}, observed {revision_line:?}",
                snapshot.declared_inventory_revision
            ),
        });
    }

    let mut sources = Vec::new();
    let mut pointers = Vec::new();
    let mut pointer_section = false;
    let mut previous_source = None;
    for (zero_indexed, line) in lines {
        let line_number = zero_indexed + 1;
        if line.is_empty() {
            return Err(HistoricalEvidenceError::MalformedManifest {
                line: line_number,
                detail: "blank rows are not canonical".to_string(),
            });
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("SOURCE") if fields.len() == 6 && !pointer_section => {
                let reference = fields[1];
                validate_manifest_path(line_number, reference)?;
                validate_scalar(line_number, "sha256", fields[3])?;
                validate_scalar(line_number, "git_blob_oid", fields[4])?;
                validate_scalar(line_number, "source_identity_blake3", fields[5])?;
                if !is_lower_hex(fields[3], 64) {
                    return Err(HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE sha256 must be 64 lowercase hex characters".to_string(),
                    });
                }
                if !is_lower_hex(fields[4], 40) {
                    return Err(HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE Git blob oid must be 40 lowercase hex characters"
                            .to_string(),
                    });
                }
                if !is_lower_hex(fields[5], 64) {
                    return Err(HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE identity must be 64 lowercase hex characters".to_string(),
                    });
                }
                let bytes = fields[2].parse::<usize>().map_err(|_| {
                    HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE byte count is not an unsigned integer".to_string(),
                    }
                })?;
                if bytes == 0 || bytes > MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES {
                    return Err(HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE byte count is outside the admitted bundle envelope"
                            .to_string(),
                    });
                }
                if previous_source.is_some_and(|previous: &str| previous >= reference) {
                    return Err(HistoricalEvidenceError::MalformedManifest {
                        line: line_number,
                        detail: "SOURCE rows must be strictly path-sorted and unique".to_string(),
                    });
                }
                previous_source = Some(reference);
                sources.push(SourceSpec {
                    reference,
                    bytes,
                    identity_blake3: fields[5],
                });
            }
            Some("POINTER") if fields.len() == 5 => {
                pointer_section = true;
                for (name, value) in [
                    ("candidate", fields[1]),
                    ("factor", fields[2]),
                    ("reference", fields[3]),
                    ("locator", fields[4]),
                ] {
                    validate_scalar(line_number, name, value)?;
                }
                validate_manifest_path(line_number, fields[3])?;
                pointers.push(PointerSpec {
                    candidate: fields[1],
                    factor: fields[2],
                    reference: fields[3],
                    locator: fields[4],
                });
            }
            Some("SOURCE") if pointer_section => {
                return Err(HistoricalEvidenceError::MalformedManifest {
                    line: line_number,
                    detail: "SOURCE rows cannot follow POINTER rows".to_string(),
                });
            }
            _ => {
                return Err(HistoricalEvidenceError::MalformedManifest {
                    line: line_number,
                    detail: format!("unexpected row shape {line:?}"),
                });
            }
        }
        if sources.len() > MAX_COMPARISON_EVIDENCE_ROWS
            || pointers.len() > MAX_COMPARISON_EVIDENCE_ROWS
        {
            return Err(HistoricalEvidenceError::MalformedManifest {
                line: line_number,
                detail: format!("row count exceeds {MAX_COMPARISON_EVIDENCE_ROWS}"),
            });
        }
    }
    if sources.len() != snapshot.source_count {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 0,
            detail: format!(
                "expected {} SOURCE rows, observed {}",
                snapshot.source_count,
                sources.len()
            ),
        });
    }
    if pointers.len() != snapshot.pointer_count {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 0,
            detail: format!(
                "expected {} POINTER rows, observed {}",
                snapshot.pointer_count,
                pointers.len()
            ),
        });
    }
    Ok(ParsedManifest {
        revision: revision_fields[1],
        sources,
        pointers,
    })
}

#[allow(clippy::too_many_lines)] // One ordered pass binds the complete caller-to-manifest pointer sequence.
fn validate_pointer_table(
    manifest: &ParsedManifest<'_>,
    snapshot: HistoricalEvidenceSnapshot,
    candidates: &[ComparisonCandidate],
) -> Result<(), HistoricalEvidenceError> {
    if manifest.revision != snapshot.declared_inventory_revision {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 2,
            detail: format!(
                "manifest revision {} does not match descriptor {}",
                manifest.revision, snapshot.declared_inventory_revision
            ),
        });
    }
    if candidates.len() > MAX_COMPARISON_EVIDENCE_ROWS {
        return Err(HistoricalEvidenceError::CandidateEvidenceTooLarge {
            observed: candidates.len(),
            maximum: MAX_COMPARISON_EVIDENCE_ROWS,
        });
    }
    let mut expected = Vec::new();
    let mut visited_evidence = 0usize;
    for candidate in candidates {
        if !valid_scalar(candidate.name) {
            return Err(HistoricalEvidenceError::InvalidCandidatePointer {
                candidate: diagnostic_scalar(candidate.name),
                factor: "<candidate>",
                reference: "<none>",
                locator: "<none>",
                detail: "candidate slug is empty, oversized, or contains a control character"
                    .to_string(),
            });
        }
        if candidate.declared_inventory_revision != snapshot.declared_inventory_revision {
            return Err(HistoricalEvidenceError::CandidateRevisionMismatch {
                candidate: candidate.name,
                expected: snapshot.declared_inventory_revision,
                observed: candidate.declared_inventory_revision,
            });
        }
        if candidate.factors.len() > MAX_COMPARISON_EVIDENCE_ROWS {
            return Err(HistoricalEvidenceError::CandidateEvidenceTooLarge {
                observed: candidate.factors.len(),
                maximum: MAX_COMPARISON_EVIDENCE_ROWS,
            });
        }
        for input in candidate.factors {
            visited_evidence = visited_evidence
                .checked_add(input.measurement.evidence.len())
                .ok_or(HistoricalEvidenceError::CandidateEvidenceTooLarge {
                    observed: usize::MAX,
                    maximum: MAX_COMPARISON_EVIDENCE_ROWS,
                })?;
            if visited_evidence > MAX_COMPARISON_EVIDENCE_ROWS {
                return Err(HistoricalEvidenceError::CandidateEvidenceTooLarge {
                    observed: visited_evidence,
                    maximum: MAX_COMPARISON_EVIDENCE_ROWS,
                });
            }
            for pointer in input
                .measurement
                .evidence
                .iter()
                .filter(|pointer| pointer.kind == EvidenceKind::WorkspacePath)
            {
                if !valid_relative_path(pointer.reference) || !valid_scalar(pointer.locator) {
                    return Err(HistoricalEvidenceError::InvalidCandidatePointer {
                        candidate: candidate.name,
                        factor: input.factor.label(),
                        reference: diagnostic_scalar(pointer.reference),
                        locator: diagnostic_scalar(pointer.locator),
                        detail: "reference must be normalized and locator must be a bounded non-control scalar"
                            .to_string(),
                    });
                }
                if expected.len() == MAX_COMPARISON_EVIDENCE_ROWS {
                    return Err(HistoricalEvidenceError::CandidateEvidenceTooLarge {
                        observed: MAX_COMPARISON_EVIDENCE_ROWS + 1,
                        maximum: MAX_COMPARISON_EVIDENCE_ROWS,
                    });
                }
                expected.push((
                    candidate.name,
                    input.factor.label(),
                    pointer.reference,
                    pointer.locator,
                ));
            }
        }
    }
    let maximum = expected.len().max(manifest.pointers.len());
    for index in 0..maximum {
        match (expected.get(index), manifest.pointers.get(index)) {
            (Some(expected_row), Some(observed_row))
                if expected_row.0 == observed_row.candidate
                    && expected_row.1 == observed_row.factor
                    && expected_row.2 == observed_row.reference
                    && expected_row.3 == observed_row.locator => {}
            (expected_row, observed_row) => {
                return Err(HistoricalEvidenceError::PointerMismatch {
                    declared_inventory_revision: snapshot.declared_inventory_revision,
                    manifest_identity_blake3: snapshot.manifest_identity_blake3,
                    index,
                    expected: expected_row.map_or_else(
                        || "<missing>".to_string(),
                        |row| format_pointer(row.0, row.1, row.2, row.3),
                    ),
                    observed: observed_row.map_or_else(
                        || "<missing>".to_string(),
                        |row| format_pointer(row.candidate, row.factor, row.reference, row.locator),
                    ),
                });
            }
        }
    }

    let referenced: BTreeSet<&str> = expected.iter().map(|row| row.2).collect();
    for source in &manifest.sources {
        if !referenced.contains(source.reference) {
            return Err(HistoricalEvidenceError::UnexpectedSource {
                reference: source.reference.to_string(),
            });
        }
    }
    Ok(())
}

#[allow(clippy::too_many_lines)] // One straight-line TAR admission pass preserves exact refusal offsets.
fn parse_tar_bundle<'a>(
    bundle: &'a [u8],
    declared_inventory_revision: &str,
) -> Result<BTreeMap<String, &'a [u8]>, HistoricalEvidenceError> {
    if bundle.len() < 3 * TAR_BLOCK_BYTES || !bundle.len().is_multiple_of(TAR_BLOCK_BYTES) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset: 0,
            detail: "TAR length must be a multiple of 512 with header and terminator blocks"
                .to_string(),
        });
    }
    let mut files = BTreeMap::new();
    let mut offset = 0usize;
    let mut saw_revision = false;
    let mut saw_terminator = false;
    while offset + TAR_BLOCK_BYTES <= bundle.len() {
        let header = &bundle[offset..offset + TAR_BLOCK_BYTES];
        if header.iter().all(|byte| *byte == 0) {
            if bundle.len() - offset < 2 * TAR_BLOCK_BYTES {
                return Err(HistoricalEvidenceError::MalformedBundle {
                    offset,
                    detail: "TAR termination requires at least two zero blocks".to_string(),
                });
            }
            if !bundle[offset..].iter().all(|byte| *byte == 0) {
                return Err(HistoricalEvidenceError::MalformedBundle {
                    offset,
                    detail: "nonzero bytes follow the first TAR terminator block".to_string(),
                });
            }
            saw_terminator = true;
            break;
        }
        validate_tar_header_profile(header, offset)?;
        validate_tar_checksum(header, offset)?;
        let size = parse_tar_octal(&header[124..136], offset, "size")?;
        let data_start = offset + TAR_BLOCK_BYTES;
        let data_end = data_start.checked_add(size).ok_or_else(|| {
            HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: "TAR payload length overflows usize".to_string(),
            }
        })?;
        if data_end > bundle.len() {
            return Err(HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: format!("TAR payload ends at {data_end}, beyond {}", bundle.len()),
            });
        }
        let padded = size.checked_add(TAR_BLOCK_BYTES - 1).ok_or_else(|| {
            HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: "TAR padded length overflows usize".to_string(),
            }
        })? / TAR_BLOCK_BYTES
            * TAR_BLOCK_BYTES;
        let next = data_start.checked_add(padded).ok_or_else(|| {
            HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: "next TAR header offset overflows usize".to_string(),
            }
        })?;
        if next > bundle.len() {
            return Err(HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: "padded TAR payload exceeds bundle".to_string(),
            });
        }
        if !bundle[data_end..next].iter().all(|byte| *byte == 0) {
            return Err(HistoricalEvidenceError::MalformedBundle {
                offset: data_end,
                detail: "TAR payload padding must be zero".to_string(),
            });
        }
        let data = &bundle[data_start..data_end];
        match header[156] {
            b'g' => {
                let path = tar_path(header, offset)?;
                if path != "pax_global_header" {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!(
                            "leading global PAX path must be \"pax_global_header\", observed {path:?}"
                        ),
                    });
                }
                let payload = std::str::from_utf8(data).map_err(|_| {
                    HistoricalEvidenceError::MalformedBundle {
                        offset: data_start,
                        detail: "global PAX header is not UTF-8".to_string(),
                    }
                })?;
                let expected = format!("52 comment={declared_inventory_revision}\n");
                if payload != expected || saw_revision || offset != 0 {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset: data_start,
                        detail: format!(
                            "expected one leading PAX revision record {expected:?}, observed {payload:?}"
                        ),
                    });
                }
                saw_revision = true;
            }
            b'5' => {
                if !saw_revision {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: "TAR entry precedes the revision PAX record".to_string(),
                    });
                }
                if size != 0 {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: "directory TAR entry has a nonzero payload".to_string(),
                    });
                }
                let path = tar_path(header, offset)?;
                if !path.ends_with('/') || path[..path.len() - 1].ends_with('/') {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!(
                            "retained TAR directory path must end in exactly one slash, observed {path:?}"
                        ),
                    });
                }
                let normalized = path.strip_suffix('/').unwrap_or(&path);
                if !valid_relative_path(normalized) {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!("directory path {path:?} is not normalized"),
                    });
                }
            }
            b'0' => {
                if !saw_revision {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: "TAR entry precedes the revision PAX record".to_string(),
                    });
                }
                let path = tar_path(header, offset)?;
                if !valid_relative_path(&path) {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!("file path {path:?} is not normalized"),
                    });
                }
                if files.insert(path.clone(), data).is_some() {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!("duplicate TAR file {path}"),
                    });
                }
                if files.len() > MAX_COMPARISON_EVIDENCE_ROWS {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!(
                            "regular-file count exceeds {MAX_COMPARISON_EVIDENCE_ROWS}"
                        ),
                    });
                }
            }
            kind => {
                return Err(HistoricalEvidenceError::MalformedBundle {
                    offset,
                    detail: format!("unsupported TAR typeflag 0x{kind:02x}"),
                });
            }
        }
        offset = next;
    }
    if !saw_revision {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset: 0,
            detail: "TAR lacks the Git revision PAX record".to_string(),
        });
    }
    if !saw_terminator {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "TAR lacks zero-block termination".to_string(),
        });
    }
    Ok(files)
}

fn validate_sources_and_markers(
    manifest: &ParsedManifest<'_>,
    files: &BTreeMap<String, &[u8]>,
    snapshot: HistoricalEvidenceSnapshot,
) -> Result<(), HistoricalEvidenceError> {
    for source in &manifest.sources {
        let bytes = match files.get(source.reference) {
            Some(bytes) if bytes.len() == source.bytes => *bytes,
            Some(bytes) => {
                return Err(HistoricalEvidenceError::SourceMismatch {
                    declared_inventory_revision: snapshot.declared_inventory_revision,
                    bundle_identity_blake3: snapshot.bundle_identity_blake3,
                    reference: source.reference.to_string(),
                    expected_bytes: Some(source.bytes),
                    observed_bytes: Some(bytes.len()),
                });
            }
            None => {
                return Err(HistoricalEvidenceError::SourceMismatch {
                    declared_inventory_revision: snapshot.declared_inventory_revision,
                    bundle_identity_blake3: snapshot.bundle_identity_blake3,
                    reference: source.reference.to_string(),
                    expected_bytes: Some(source.bytes),
                    observed_bytes: None,
                });
            }
        };
        let observed = source_identity(
            snapshot.declared_inventory_revision,
            source.reference,
            bytes,
        );
        if observed != source.identity_blake3 {
            let pointer = first_pointer_for_source(manifest, source.reference)?;
            return Err(HistoricalEvidenceError::SourceIdentityMismatch {
                declared_inventory_revision: snapshot.declared_inventory_revision,
                candidate: Box::from(pointer.candidate),
                factor: Box::from(pointer.factor),
                reference: Box::from(source.reference),
                locator: Box::from(pointer.locator),
                expected_blake3: Box::from(source.identity_blake3),
                observed_blake3: observed.into_boxed_str(),
            });
        }
    }
    for reference in files.keys() {
        if !manifest
            .sources
            .iter()
            .any(|source| source.reference == reference)
        {
            return Err(HistoricalEvidenceError::SourceMismatch {
                declared_inventory_revision: snapshot.declared_inventory_revision,
                bundle_identity_blake3: snapshot.bundle_identity_blake3,
                reference: reference.clone(),
                expected_bytes: None,
                observed_bytes: files.get(reference).map(|bytes| bytes.len()),
            });
        }
    }
    let marker_scan_work = marker_scan_work_upper_bound(manifest, files, snapshot)?;
    if marker_scan_work > MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK {
        return Err(HistoricalEvidenceError::MarkerScanBudgetExceeded {
            observed: marker_scan_work,
            maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
        });
    }
    for pointer in &manifest.pointers {
        let Some(bytes) = files.get(pointer.reference) else {
            return Err(HistoricalEvidenceError::SourceMismatch {
                declared_inventory_revision: snapshot.declared_inventory_revision,
                bundle_identity_blake3: snapshot.bundle_identity_blake3,
                reference: pointer.reference.to_string(),
                expected_bytes: None,
                observed_bytes: None,
            });
        };
        if !contains_bytes(bytes, pointer.locator.as_bytes()) {
            let source = manifest
                .sources
                .iter()
                .find(|source| source.reference == pointer.reference)
                .ok_or_else(|| HistoricalEvidenceError::UnexpectedSource {
                    reference: pointer.reference.to_string(),
                })?;
            return Err(HistoricalEvidenceError::MarkerMissing {
                declared_inventory_revision: snapshot.declared_inventory_revision,
                source_identity_blake3: Box::from(source.identity_blake3),
                candidate: Box::from(pointer.candidate),
                factor: Box::from(pointer.factor),
                reference: Box::from(pointer.reference),
                locator: Box::from(pointer.locator),
            });
        }
    }
    Ok(())
}

fn marker_scan_work_upper_bound(
    manifest: &ParsedManifest<'_>,
    files: &BTreeMap<String, &[u8]>,
    snapshot: HistoricalEvidenceSnapshot,
) -> Result<u64, HistoricalEvidenceError> {
    let mut total = 0u64;
    for pointer in &manifest.pointers {
        let Some(bytes) = files.get(pointer.reference) else {
            return Err(HistoricalEvidenceError::SourceMismatch {
                declared_inventory_revision: snapshot.declared_inventory_revision,
                bundle_identity_blake3: snapshot.bundle_identity_blake3,
                reference: pointer.reference.to_string(),
                expected_bytes: None,
                observed_bytes: None,
            });
        };
        let source_bytes = u64::try_from(bytes.len()).map_err(|_| {
            HistoricalEvidenceError::MarkerScanBudgetExceeded {
                observed: u64::MAX,
                maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
            }
        })?;
        let locator_bytes = u64::try_from(pointer.locator.len()).map_err(|_| {
            HistoricalEvidenceError::MarkerScanBudgetExceeded {
                observed: u64::MAX,
                maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
            }
        })?;
        let pointer_work = source_bytes
            .checked_add(locator_bytes)
            .and_then(|work| work.checked_mul(2))
            .ok_or(HistoricalEvidenceError::MarkerScanBudgetExceeded {
                observed: u64::MAX,
                maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
            })?;
        total = total.checked_add(pointer_work).ok_or(
            HistoricalEvidenceError::MarkerScanBudgetExceeded {
                observed: u64::MAX,
                maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
            },
        )?;
    }
    Ok(total)
}

fn source_identity(declared_inventory_revision: &str, reference: &str, bytes: &[u8]) -> String {
    let mut hasher = DomainHasher::new(COMPARISON_EVIDENCE_SOURCE_IDENTITY_DOMAIN);
    hasher.update(declared_inventory_revision.as_bytes());
    hasher.update(&[0]);
    hasher.update(reference.as_bytes());
    hasher.update(&[0]);
    hasher.update(bytes);
    hasher.finalize().to_hex()
}

fn first_pointer_for_source<'manifest, 'data>(
    manifest: &'manifest ParsedManifest<'data>,
    reference: &str,
) -> Result<&'manifest PointerSpec<'data>, HistoricalEvidenceError> {
    manifest
        .pointers
        .iter()
        .find(|pointer| pointer.reference == reference)
        .ok_or_else(|| HistoricalEvidenceError::UnexpectedSource {
            reference: reference.to_string(),
        })
}

fn validate_tar_header_profile(
    header: &[u8],
    offset: usize,
) -> Result<(), HistoricalEvidenceError> {
    if &header[257..263] != b"ustar\0" {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "TAR magic must be the retained USTAR encoding \"ustar\\0\"".to_string(),
        });
    }
    if &header[263..265] != b"00" {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "TAR version must be the retained USTAR encoding \"00\"".to_string(),
        });
    }

    validate_tar_text_padding(&header[..100], offset, "name")?;
    validate_tar_octal_encoding(&header[100..108], offset, "mode")?;
    validate_tar_octal_encoding(&header[108..116], offset, "uid")?;
    validate_tar_octal_encoding(&header[116..124], offset, "gid")?;
    validate_tar_octal_encoding(&header[124..136], offset, "size")?;
    validate_tar_octal_encoding(&header[136..148], offset, "mtime")?;
    validate_tar_octal_encoding(&header[148..156], offset, "checksum")?;
    validate_tar_text_padding(&header[157..257], offset, "linkname")?;
    validate_tar_text_padding(&header[265..297], offset, "uname")?;
    validate_tar_text_padding(&header[297..329], offset, "gname")?;
    validate_tar_octal_encoding(&header[329..337], offset, "devmajor")?;
    validate_tar_octal_encoding(&header[337..345], offset, "devminor")?;
    validate_tar_text_padding(&header[345..500], offset, "prefix")?;

    let expected_mode = match header[156] {
        b'g' => b"0000666\0".as_slice(),
        b'0' => b"0000664\0".as_slice(),
        b'5' => b"0000775\0".as_slice(),
        kind => {
            return Err(HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: format!(
                    "retained TAR typeflag must be 'g', '0', or '5', observed 0x{kind:02x}"
                ),
            });
        }
    };
    if &header[100..108] != expected_mode {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!(
                "retained TAR mode does not match typeflag 0x{:02x}",
                header[156]
            ),
        });
    }
    if !header[157..257].iter().all(|byte| *byte == 0) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "retained TAR linkname field must be empty".to_string(),
        });
    }
    if &header[265..297] != padded_tar_text::<32>(b"root").as_slice() {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "retained TAR uname field must be \"root\" with zero padding".to_string(),
        });
    }
    if &header[297..329] != padded_tar_text::<32>(b"root").as_slice() {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "retained TAR gname field must be \"root\" with zero padding".to_string(),
        });
    }
    for (field, name) in [
        (&header[108..116], "uid"),
        (&header[116..124], "gid"),
        (&header[329..337], "devmajor"),
        (&header[337..345], "devminor"),
    ] {
        if field != b"0000000\0" {
            return Err(HistoricalEvidenceError::MalformedBundle {
                offset,
                detail: format!("retained TAR {name} field must be canonical zero"),
            });
        }
    }
    if !header[345..500].iter().all(|byte| *byte == 0) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "retained TAR prefix field must be empty".to_string(),
        });
    }
    if !header[500..].iter().all(|byte| *byte == 0) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: "retained TAR reserved header padding must be zero".to_string(),
        });
    }
    Ok(())
}

fn padded_tar_text<const N: usize>(text: &[u8]) -> [u8; N] {
    let mut field = [0; N];
    field[..text.len()].copy_from_slice(text);
    field
}

fn validate_tar_octal_encoding(
    field: &[u8],
    offset: usize,
    name: &str,
) -> Result<(), HistoricalEvidenceError> {
    let Some((&terminator, digits)) = field.split_last() else {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR {name} field is empty"),
        });
    };
    if terminator != 0 || !digits.iter().all(|byte| (b'0'..=b'7').contains(byte)) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!(
                "TAR {name} field must be fixed-width octal digits followed by one NUL"
            ),
        });
    }
    Ok(())
}

fn validate_tar_text_padding(
    field: &[u8],
    offset: usize,
    name: &str,
) -> Result<(), HistoricalEvidenceError> {
    let Some(terminator) = field.iter().position(|byte| *byte == 0) else {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR {name} field must be NUL-terminated"),
        });
    };
    if !field[terminator..].iter().all(|byte| *byte == 0) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR {name} field has nonzero bytes after its first NUL"),
        });
    }
    Ok(())
}

fn validate_tar_checksum(header: &[u8], offset: usize) -> Result<(), HistoricalEvidenceError> {
    let expected = parse_tar_octal(&header[148..156], offset, "checksum")?;
    let observed: usize = header
        .iter()
        .enumerate()
        .map(|(index, byte)| {
            if (148..156).contains(&index) {
                usize::from(b' ')
            } else {
                usize::from(*byte)
            }
        })
        .sum();
    if expected != observed {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR checksum expected {expected}, observed {observed}"),
        });
    }
    Ok(())
}

fn parse_tar_octal(
    field: &[u8],
    offset: usize,
    name: &str,
) -> Result<usize, HistoricalEvidenceError> {
    let end = field
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(field.len());
    let text = std::str::from_utf8(&field[..end])
        .map_err(|_| HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR {name} field is not ASCII"),
        })?
        .trim();
    if text.is_empty() || !text.bytes().all(|byte| (b'0'..=b'7').contains(&byte)) {
        return Err(HistoricalEvidenceError::MalformedBundle {
            offset,
            detail: format!("TAR {name} field {text:?} is not canonical octal"),
        });
    }
    usize::from_str_radix(text, 8).map_err(|_| HistoricalEvidenceError::MalformedBundle {
        offset,
        detail: format!("TAR {name} field {text:?} overflows usize"),
    })
}

fn tar_path(header: &[u8], offset: usize) -> Result<String, HistoricalEvidenceError> {
    let name = tar_string(&header[..100], offset, "name")?;
    let prefix = tar_string(&header[345..500], offset, "prefix")?;
    if prefix.is_empty() {
        Ok(name.to_string())
    } else {
        Ok(format!("{prefix}/{name}"))
    }
}

fn tar_string<'a>(
    field: &'a [u8],
    offset: usize,
    name: &str,
) -> Result<&'a str, HistoricalEvidenceError> {
    validate_tar_text_padding(field, offset, name)?;
    let end = field
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(field.len());
    std::str::from_utf8(&field[..end]).map_err(|_| HistoricalEvidenceError::MalformedBundle {
        offset,
        detail: format!("TAR {name} is not UTF-8"),
    })
}

fn validate_manifest_path(line: usize, path: &str) -> Result<(), HistoricalEvidenceError> {
    validate_scalar(line, "reference", path)?;
    if valid_relative_path(path) {
        Ok(())
    } else {
        Err(HistoricalEvidenceError::MalformedManifest {
            line,
            detail: format!("path {path:?} is not normalized workspace-relative"),
        })
    }
}

fn validate_scalar(line: usize, field: &str, value: &str) -> Result<(), HistoricalEvidenceError> {
    if value.is_empty() {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line,
            detail: format!("{field} is empty"),
        });
    }
    if value.len() > MAX_COMPARISON_EVIDENCE_FIELD_BYTES {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line,
            detail: format!(
                "{field} is {} bytes; maximum is {MAX_COMPARISON_EVIDENCE_FIELD_BYTES}",
                value.len()
            ),
        });
    }
    if value.chars().any(char::is_control) {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line,
            detail: format!("{field} contains a control character"),
        });
    }
    Ok(())
}

fn valid_relative_path(path: &str) -> bool {
    !path.is_empty()
        && !path.starts_with('/')
        && !path.contains('\\')
        && path.len() <= MAX_COMPARISON_EVIDENCE_FIELD_BYTES
        && path
            .split('/')
            .all(|component| !component.is_empty() && component != "." && component != "..")
        && !path.chars().any(char::is_control)
}

fn valid_scalar(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= MAX_COMPARISON_EVIDENCE_FIELD_BYTES
        && !value.chars().any(char::is_control)
}

fn diagnostic_scalar(value: &'static str) -> &'static str {
    if valid_scalar(value) {
        value
    } else {
        "<invalid-or-oversized>"
    }
}

fn is_lower_hex(value: &str, length: usize) -> bool {
    value.len() == length
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
    #[cfg(test)]
    MARKER_SEARCH_CALLS.with(|calls| calls.set(calls.get() + 1));
    contains_bytes_counted(haystack, needle, || {})
}

fn contains_bytes_counted(
    haystack: &[u8],
    needle: &[u8],
    mut count_comparison: impl FnMut(),
) -> bool {
    if needle.is_empty() || needle.len() > haystack.len() {
        return false;
    }

    let mut prefix = vec![0usize; needle.len()];
    let mut matched = 0usize;
    for index in 1..needle.len() {
        loop {
            count_comparison();
            if needle[index] == needle[matched] {
                matched += 1;
                prefix[index] = matched;
                break;
            }
            if matched == 0 {
                break;
            }
            matched = prefix[matched - 1];
        }
    }

    matched = 0;
    for byte in haystack {
        loop {
            count_comparison();
            if *byte == needle[matched] {
                matched += 1;
                if matched == needle.len() {
                    return true;
                }
                break;
            }
            if matched == 0 {
                break;
            }
            matched = prefix[matched - 1];
        }
    }
    false
}

fn format_pointer(candidate: &str, factor: &str, reference: &str, locator: &str) -> String {
    format!("{candidate}/{factor}/{reference}/{locator}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_source_identities_match_manifest() {
        let manifest =
            parse_manifest(MANIFEST_BYTES, COMPARISON_EVIDENCE_SNAPSHOT).expect("manifest parses");
        let files = parse_tar_bundle(
            BUNDLE_BYTES,
            COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
        )
        .expect("bundle parses");
        let mismatches: Vec<String> = manifest
            .sources
            .iter()
            .filter_map(|source| {
                let bytes = files.get(source.reference).expect("source exists");
                let observed = source_identity(
                    COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
                    source.reference,
                    bytes,
                );
                (observed != source.identity_blake3)
                    .then(|| format!("{}\t{observed}", source.reference))
            })
            .collect();
        assert!(
            mismatches.is_empty(),
            "per-source BLAKE3 mismatches:\n{}",
            mismatches.join("\n")
        );
    }

    #[test]
    fn embedded_artifact_identities_match_descriptor() {
        let manifest = hash_domain(
            COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN,
            MANIFEST_BYTES.as_bytes(),
        )
        .to_hex();
        let bundle = hash_domain(COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN, BUNDLE_BYTES).to_hex();
        assert_eq!(
            manifest, COMPARISON_EVIDENCE_SNAPSHOT.manifest_identity_blake3,
            "manifest identity drift"
        );
        assert_eq!(
            bundle, COMPARISON_EVIDENCE_SNAPSHOT.bundle_identity_blake3,
            "bundle identity drift"
        );
    }

    #[test]
    fn embedded_descriptor_identity_matches_descriptor() {
        let observed = comparison_evidence_descriptor_identity_blake3(COMPARISON_EVIDENCE_SNAPSHOT);
        assert_ne!(
            COMPARISON_EVIDENCE_SNAPSHOT.descriptor_identity_blake3,
            "0000000000000000000000000000000000000000000000000000000000000000",
            "descriptor identity must be a pinned non-placeholder literal"
        );
        assert_eq!(
            observed, COMPARISON_EVIDENCE_SNAPSHOT.descriptor_identity_blake3,
            "descriptor identity drift"
        );
    }

    #[test]
    fn descriptor_identity_binds_every_other_snapshot_field() {
        let baseline = comparison_evidence_descriptor_identity_blake3(COMPARISON_EVIDENCE_SNAPSHOT);
        let mutations = [
            (
                "schema",
                HistoricalEvidenceSnapshot {
                    schema: "frankensim-wedge-comparison-evidence-snapshot-v4",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "policy_version",
                HistoricalEvidenceSnapshot {
                    policy_version: 3,
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "declared_inventory_revision",
                HistoricalEvidenceSnapshot {
                    declared_inventory_revision: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "manifest_path",
                HistoricalEvidenceSnapshot {
                    manifest_path: "retained/manifest.tsv",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "manifest_identity_blake3",
                HistoricalEvidenceSnapshot {
                    manifest_identity_blake3: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "bundle_path",
                HistoricalEvidenceSnapshot {
                    bundle_path: "retained/bundle.tar",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "bundle_identity_blake3",
                HistoricalEvidenceSnapshot {
                    bundle_identity_blake3: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "comparison_model_identity_blake3",
                HistoricalEvidenceSnapshot {
                    comparison_model_identity_blake3: "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc",
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "source_count",
                HistoricalEvidenceSnapshot {
                    source_count: COMPARISON_EVIDENCE_SNAPSHOT.source_count + 1,
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
            (
                "pointer_count",
                HistoricalEvidenceSnapshot {
                    pointer_count: COMPARISON_EVIDENCE_SNAPSHOT.pointer_count + 1,
                    ..COMPARISON_EVIDENCE_SNAPSHOT
                },
            ),
        ];
        for (field, mutation) in mutations {
            assert_ne!(
                comparison_evidence_descriptor_identity_blake3(mutation),
                baseline,
                "descriptor identity omitted {field}"
            );
        }
        let self_only_mutation = HistoricalEvidenceSnapshot {
            descriptor_identity_blake3: "dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd",
            ..COMPARISON_EVIDENCE_SNAPSHOT
        };
        assert_eq!(
            comparison_evidence_descriptor_identity_blake3(self_only_mutation),
            baseline,
            "descriptor identity must exclude its own root"
        );
    }

    #[test]
    fn embedded_comparison_model_identity_matches_descriptor() {
        let observed = comparison_model_identity_blake3(
            &DEFAULT_FACTOR_WEIGHTS,
            crate::comparison_candidates(),
        )
        .expect("embedded comparison model is structurally valid");
        assert_ne!(
            COMPARISON_EVIDENCE_SNAPSHOT.comparison_model_identity_blake3,
            "0000000000000000000000000000000000000000000000000000000000000000",
            "comparison-model identity must be a pinned non-placeholder literal"
        );
        assert_eq!(
            observed, COMPARISON_EVIDENCE_SNAPSHOT.comparison_model_identity_blake3,
            "comparison-model identity drift"
        );
    }

    #[test]
    fn kmp_marker_search_handles_repeated_prefix_near_matches() {
        let counterexample_needle = b"aab";
        let counterexample_source = b"aaaaaac";
        let mut counterexample_comparisons = 0u64;
        assert!(!contains_bytes_counted(
            counterexample_source,
            counterexample_needle,
            || counterexample_comparisons += 1,
        ));
        assert_eq!(
            counterexample_comparisons, 16,
            "the single-comparison transition shape must remain explicit"
        );
        assert!(
            counterexample_comparisons
                <= 2 * (counterexample_source.len() + counterexample_needle.len()) as u64
        );

        let mut needle = vec![b'a'; MAX_COMPARISON_EVIDENCE_FIELD_BYTES];
        *needle.last_mut().expect("bounded locator is nonempty") = b'b';

        let mut absent = vec![b'a'; 64 * 1024];
        *absent.last_mut().expect("test source is nonempty") = b'c';
        let mut absent_comparisons = 0u64;
        assert!(!contains_bytes_counted(&absent, &needle, || {
            absent_comparisons += 1;
        }));
        assert!(absent_comparisons <= 2 * (absent.len() + needle.len()) as u64);

        let mut present_at_end = vec![b'c'; 64 * 1024];
        present_at_end.extend_from_slice(&needle);
        assert!(contains_bytes(&present_at_end, &needle));
    }

    #[test]
    fn exact_cap_near_matches_stay_within_charged_comparisons() {
        const POINTERS: usize = 8;
        let reference = "bounded/source.txt";
        let mut source =
            vec![b'a'; MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES - MAX_COMPARISON_EVIDENCE_FIELD_BYTES];
        *source.last_mut().expect("cap-boundary source is nonempty") = b'c';
        let mut locator = vec![b'a'; MAX_COMPARISON_EVIDENCE_FIELD_BYTES];
        *locator.last_mut().expect("bounded locator is nonempty") = b'b';
        let locator_text = std::str::from_utf8(&locator).expect("test locator is UTF-8");
        let identity = source_identity(
            COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
            reference,
            &source,
        );
        let sources = vec![SourceSpec {
            reference,
            bytes: source.len(),
            identity_blake3: &identity,
        }];
        let pointers = (0..POINTERS)
            .map(|_| PointerSpec {
                candidate: "adversarial-candidate",
                factor: "adversarial-factor",
                reference,
                locator: locator_text,
            })
            .collect();
        let manifest = ParsedManifest {
            revision: COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
            sources,
            pointers,
        };
        let mut files = BTreeMap::new();
        files.insert(reference.to_string(), source.as_slice());
        assert_eq!(
            marker_scan_work_upper_bound(&manifest, &files, COMPARISON_EVIDENCE_SNAPSHOT)
                .expect("cap-boundary work is representable"),
            MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK
        );

        let mut comparisons = 0u64;
        for _ in 0..POINTERS {
            assert!(!contains_bytes_counted(&source, &locator, || {
                comparisons += 1;
            }));
        }
        assert!(comparisons <= MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK);
    }

    #[test]
    fn maximum_row_preflight_refuses_before_any_marker_search() {
        let reference = "bounded/source.txt";
        let source = vec![b'a'; 32 * 1024];
        let mut locator = "a".repeat(64);
        locator.replace_range(locator.len() - 1.., "b");
        let identity = source_identity(
            COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
            reference,
            &source,
        );
        let sources = vec![SourceSpec {
            reference,
            bytes: source.len(),
            identity_blake3: &identity,
        }];
        let pointers = (0..MAX_COMPARISON_EVIDENCE_ROWS)
            .map(|_| PointerSpec {
                candidate: "adversarial-candidate",
                factor: "adversarial-factor",
                reference,
                locator: &locator,
            })
            .collect();
        let manifest = ParsedManifest {
            revision: COMPARISON_EVIDENCE_SNAPSHOT.declared_inventory_revision,
            sources,
            pointers,
        };
        let mut files = BTreeMap::new();
        files.insert(reference.to_string(), source.as_slice());

        MARKER_SEARCH_CALLS.with(|calls| calls.set(0));
        let error = validate_sources_and_markers(&manifest, &files, COMPARISON_EVIDENCE_SNAPSHOT)
            .expect_err("aggregate work must refuse before any near-match search");
        assert_eq!(
            MARKER_SEARCH_CALLS.with(std::cell::Cell::get),
            0,
            "preflight refusal must precede every marker search"
        );
        let expected = u64::try_from(MAX_COMPARISON_EVIDENCE_ROWS).expect("row cap fits u64")
            * 2
            * (u64::try_from(source.len()).expect("source cap fits u64")
                + u64::try_from(locator.len()).expect("locator cap fits u64"));
        assert_eq!(
            error,
            HistoricalEvidenceError::MarkerScanBudgetExceeded {
                observed: expected,
                maximum: MAX_COMPARISON_EVIDENCE_MARKER_SCAN_WORK,
            }
        );
    }
}

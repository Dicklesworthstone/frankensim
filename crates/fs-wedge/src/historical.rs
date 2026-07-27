//! Offline replay for the exact source snapshot behind the measured comparison.
//!
//! Live readiness probes intentionally inspect the current workspace. Historical
//! comparison factors do not: they replay against a retained archive whose
//! revision, manifest, source bytes, and pointer set are all content-bound.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{DomainHasher, hash_domain};

use crate::{ComparisonCandidate, EvidenceKind};

/// Schema carried by the retained comparison-evidence manifest.
pub const COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA: &str =
    "frankensim-wedge-comparison-evidence-snapshot-v1";
/// Parser and admission policy version for the retained source bundle.
pub const COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION: u32 = 1;
/// Domain for the exact manifest-byte identity.
pub const COMPARISON_EVIDENCE_MANIFEST_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-manifest.v1";
/// Domain for the exact retained TAR-byte identity.
pub const COMPARISON_EVIDENCE_BUNDLE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-bundle.v1";
/// Domain for one revision/path/source-byte identity.
pub const COMPARISON_EVIDENCE_SOURCE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-wedge.comparison-evidence-source.v1";

/// Maximum retained source-bundle bytes accepted by the replay adapter.
pub const MAX_COMPARISON_EVIDENCE_BUNDLE_BYTES: usize = 1024 * 1024;
/// Maximum manifest bytes accepted by the replay adapter.
pub const MAX_COMPARISON_EVIDENCE_MANIFEST_BYTES: usize = 128 * 1024;
/// Maximum source or pointer rows accepted from one manifest.
pub const MAX_COMPARISON_EVIDENCE_ROWS: usize = 512;
/// Maximum UTF-8 bytes accepted in one manifest scalar.
pub const MAX_COMPARISON_EVIDENCE_FIELD_BYTES: usize = 4096;

const TAR_BLOCK_BYTES: usize = 512;
const MANIFEST_BYTES: &str = include_str!("../data/comparison-evidence-b3b5f2c1.tsv");
const BUNDLE_BYTES: &[u8] = include_bytes!("../data/comparison-evidence-b3b5f2c1.tar");

/// Immutable adapter metadata for one retained comparison source snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HistoricalEvidenceSnapshot {
    /// Stable schema label.
    pub schema: &'static str,
    /// Replay-policy version.
    pub policy_version: u32,
    /// Git revision whose exact file bytes are retained.
    pub inventory_revision: &'static str,
    /// Workspace-relative retained-manifest path.
    pub manifest_path: &'static str,
    /// Domain-separated BLAKE3 of the exact manifest bytes.
    pub manifest_identity_blake3: &'static str,
    /// Workspace-relative retained-bundle path.
    pub bundle_path: &'static str,
    /// Domain-separated BLAKE3 of the exact TAR bytes.
    pub bundle_identity_blake3: &'static str,
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
            && is_lower_hex(self.inventory_revision, 40)
            && valid_relative_path(self.manifest_path)
            && is_lower_hex(self.manifest_identity_blake3, 64)
            && valid_relative_path(self.bundle_path)
            && is_lower_hex(self.bundle_identity_blake3, 64)
            && self.source_count > 0
            && self.source_count <= MAX_COMPARISON_EVIDENCE_ROWS
            && self.pointer_count > 0
            && self.pointer_count <= MAX_COMPARISON_EVIDENCE_ROWS
    }
}

/// The source snapshot bound to the default measured comparison.
///
/// The two BLAKE3 identities are intentionally literal protocol fields. Tests
/// recompute them from the embedded artifacts before any source row can be
/// consumed.
pub const COMPARISON_EVIDENCE_SNAPSHOT: HistoricalEvidenceSnapshot = HistoricalEvidenceSnapshot {
    schema: COMPARISON_EVIDENCE_SNAPSHOT_SCHEMA,
    policy_version: COMPARISON_EVIDENCE_SNAPSHOT_POLICY_VERSION,
    inventory_revision: "b3b5f2c1c809eec06cde1e40cbc916d6995469b5",
    manifest_path: "crates/fs-wedge/data/comparison-evidence-b3b5f2c1.tsv",
    manifest_identity_blake3: "2ea962dc43b416b22e1591e23b065d20a2c93be4b1440aeff7228e9c37ace6ec",
    bundle_path: "crates/fs-wedge/data/comparison-evidence-b3b5f2c1.tar",
    bundle_identity_blake3: "4a13f162b435a126c979c1da6a8743e4c11388648e6e225ebf577da663f0e9e3",
    source_count: 13,
    pointer_count: 31,
};

/// Successful offline replay of a comparison source snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HistoricalEvidenceReceipt {
    schema: &'static str,
    policy_version: u32,
    inventory_revision: &'static str,
    manifest_identity_blake3: String,
    bundle_identity_blake3: String,
    source_count: usize,
    pointer_count: usize,
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

    /// Exact reviewed Git revision.
    #[must_use]
    pub const fn inventory_revision(&self) -> &'static str {
        self.inventory_revision
    }

    /// Exact domain-separated manifest identity.
    #[must_use]
    pub fn manifest_identity_blake3(&self) -> &str {
        &self.manifest_identity_blake3
    }

    /// Exact domain-separated source-bundle identity.
    #[must_use]
    pub fn bundle_identity_blake3(&self) -> &str {
        &self.bundle_identity_blake3
    }

    /// Number of authenticated source files.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }

    /// Number of replayed historical pointer occurrences.
    #[must_use]
    pub const fn pointer_count(&self) -> usize {
        self.pointer_count
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
    /// The canonical TSV manifest is malformed.
    MalformedManifest {
        /// One-based line number, or zero for whole-manifest failures.
        line: usize,
        /// Actionable detail.
        detail: String,
    },
    /// A candidate claims a revision other than the authenticated snapshot.
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
        /// Authenticated snapshot revision.
        inventory_revision: &'static str,
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
        /// Authenticated snapshot revision.
        inventory_revision: &'static str,
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
        /// Authenticated snapshot revision.
        inventory_revision: &'static str,
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
        /// Authenticated snapshot revision.
        inventory_revision: &'static str,
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
            Self::MalformedManifest { line, detail } => {
                write!(f, "malformed historical manifest at line {line}: {detail}")
            }
            Self::CandidateRevisionMismatch {
                candidate,
                expected,
                observed,
            } => write!(
                f,
                "candidate {candidate} uses inventory revision {observed}; authenticated snapshot is {expected}"
            ),
            Self::PointerMismatch {
                inventory_revision,
                manifest_identity_blake3,
                index,
                expected,
                observed,
            } => write!(
                f,
                "historical pointer {index} mismatch at revision {inventory_revision} \
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
                inventory_revision,
                bundle_identity_blake3,
                reference,
                expected_bytes,
                observed_bytes,
            } => write!(
                f,
                "historical source {reference} set/length mismatch at revision \
                 {inventory_revision} (bundle BLAKE3 {bundle_identity_blake3}): \
                 expected {expected_bytes:?} bytes, observed {observed_bytes:?}"
            ),
            Self::SourceIdentityMismatch {
                inventory_revision,
                candidate,
                factor,
                reference,
                locator,
                expected_blake3,
                observed_blake3,
            } => write!(
                f,
                "historical source identity mismatch at revision {inventory_revision} for \
                 {candidate}/{factor}/{reference}/{locator}: expected BLAKE3 \
                 {expected_blake3}, observed {observed_blake3}"
            ),
            Self::MarkerMissing {
                inventory_revision,
                source_identity_blake3,
                candidate,
                factor,
                reference,
                locator,
            } => write!(
                f,
                "historical candidate {candidate} factor {factor} lost marker {locator:?} in \
                 {reference} at revision {inventory_revision} after verifying source BLAKE3 \
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

/// Replay the default comparison against its retained source bytes.
pub fn verify_default_comparison_evidence()
-> Result<HistoricalEvidenceReceipt, HistoricalEvidenceError> {
    verify_comparison_evidence(
        COMPARISON_EVIDENCE_SNAPSHOT,
        MANIFEST_BYTES,
        BUNDLE_BYTES,
        crate::comparison_candidates(),
    )
}

/// Replay supplied comparison records against one exact retained source set.
///
/// This adapter reads no filesystem or Git state. Callers may supply alternate
/// bytes for mutation and migration tests, but the descriptor must authenticate
/// those bytes before any marker can be admitted.
pub fn verify_comparison_evidence(
    snapshot: HistoricalEvidenceSnapshot,
    manifest: &str,
    bundle: &[u8],
    candidates: &[ComparisonCandidate],
) -> Result<HistoricalEvidenceReceipt, HistoricalEvidenceError> {
    validate_descriptor(snapshot)?;
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
    let files = parse_tar_bundle(bundle, snapshot.inventory_revision)?;
    validate_sources_and_markers(&parsed, &files, snapshot)?;

    Ok(HistoricalEvidenceReceipt {
        schema: snapshot.schema,
        policy_version: snapshot.policy_version,
        inventory_revision: snapshot.inventory_revision,
        manifest_identity_blake3: observed_manifest,
        bundle_identity_blake3: observed_bundle,
        source_count: parsed.sources.len(),
        pointer_count: parsed.pointers.len(),
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
        ("inventory_revision", snapshot.inventory_revision, 40),
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
    ] {
        if !is_lower_hex(value, length) {
            return Err(HistoricalEvidenceError::InvalidDescriptor {
                field,
                detail: format!("must be exactly {length} lowercase hexadecimal characters"),
            });
        }
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
    if revision_fields.as_slice() != ["REVISION", snapshot.inventory_revision] {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 2,
            detail: format!(
                "expected REVISION {}, observed {revision_line:?}",
                snapshot.inventory_revision
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
    if manifest.revision != snapshot.inventory_revision {
        return Err(HistoricalEvidenceError::MalformedManifest {
            line: 2,
            detail: format!(
                "manifest revision {} does not match descriptor {}",
                manifest.revision, snapshot.inventory_revision
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
        if candidate.inventory_revision != snapshot.inventory_revision {
            return Err(HistoricalEvidenceError::CandidateRevisionMismatch {
                candidate: candidate.name,
                expected: snapshot.inventory_revision,
                observed: candidate.inventory_revision,
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
                    inventory_revision: snapshot.inventory_revision,
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
    inventory_revision: &str,
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
                let payload = std::str::from_utf8(data).map_err(|_| {
                    HistoricalEvidenceError::MalformedBundle {
                        offset: data_start,
                        detail: "global PAX header is not UTF-8".to_string(),
                    }
                })?;
                let expected = format!("52 comment={inventory_revision}\n");
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
                let normalized = path.strip_suffix('/').unwrap_or(&path);
                if !valid_relative_path(normalized) {
                    return Err(HistoricalEvidenceError::MalformedBundle {
                        offset,
                        detail: format!("directory path {path:?} is not normalized"),
                    });
                }
            }
            0 | b'0' => {
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
                    inventory_revision: snapshot.inventory_revision,
                    bundle_identity_blake3: snapshot.bundle_identity_blake3,
                    reference: source.reference.to_string(),
                    expected_bytes: Some(source.bytes),
                    observed_bytes: Some(bytes.len()),
                });
            }
            None => {
                return Err(HistoricalEvidenceError::SourceMismatch {
                    inventory_revision: snapshot.inventory_revision,
                    bundle_identity_blake3: snapshot.bundle_identity_blake3,
                    reference: source.reference.to_string(),
                    expected_bytes: Some(source.bytes),
                    observed_bytes: None,
                });
            }
        };
        let observed = source_identity(snapshot.inventory_revision, source.reference, bytes);
        if observed != source.identity_blake3 {
            let pointer = first_pointer_for_source(manifest, source.reference)?;
            return Err(HistoricalEvidenceError::SourceIdentityMismatch {
                inventory_revision: snapshot.inventory_revision,
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
                inventory_revision: snapshot.inventory_revision,
                bundle_identity_blake3: snapshot.bundle_identity_blake3,
                reference: reference.clone(),
                expected_bytes: None,
                observed_bytes: files.get(reference).map(|bytes| bytes.len()),
            });
        }
    }
    for pointer in &manifest.pointers {
        let Some(bytes) = files.get(pointer.reference) else {
            return Err(HistoricalEvidenceError::SourceMismatch {
                inventory_revision: snapshot.inventory_revision,
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
                inventory_revision: snapshot.inventory_revision,
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

fn source_identity(inventory_revision: &str, reference: &str, bytes: &[u8]) -> String {
    let mut hasher = DomainHasher::new(COMPARISON_EVIDENCE_SOURCE_IDENTITY_DOMAIN);
    hasher.update(inventory_revision.as_bytes());
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
    !needle.is_empty()
        && needle.len() <= haystack.len()
        && haystack
            .windows(needle.len())
            .any(|window| window == needle)
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
            COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision,
        )
        .expect("bundle parses");
        let mismatches: Vec<String> = manifest
            .sources
            .iter()
            .filter_map(|source| {
                let bytes = files.get(source.reference).expect("source exists");
                let observed = source_identity(
                    COMPARISON_EVIDENCE_SNAPSHOT.inventory_revision,
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
}

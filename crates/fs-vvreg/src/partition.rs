//! Purpose-typed corpus access and calibration-taint enforcement (EXTREAL
//! E07, bead `frankensim-extreal-program-f85xj.7.1`).
//!
//! A caller does not obtain a dataset by merely repeating its declared
//! partition. It must state the intended use. Training and calibration rows
//! are admitted only for [`DatasetPurpose::Calibration`] and both taint every
//! derived model; validation rows are admitted only for
//! [`DatasetPurpose::Validation`]. Blind rows additionally require a
//! generation-bound release receipt.
//!
//! All five receipt families share a closed, self-verifying wire envelope for
//! generic content-addressed persistence. This UTIL-layer module still retains
//! only an ordered in-memory partition event log; the wire envelope does not
//! claim dedicated ledger membership, audit-event atomicity, authorization, or
//! partition-state restoration.

use crate::corpus::{
    ContextValue, CorpusDataset, CorpusQueryRefusal, CorpusRegistry, DatasetPartition,
    MAX_CORPUS_TEXT_BYTES,
};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::Evidence;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::Deref;

/// Current purpose/taint/repartition semantic version.
pub const PARTITION_POLICY_SCHEMA_VERSION: u32 = 1;
/// Maximum datasets, direct inputs, parents, or validation inputs in one
/// bounded operation.
pub const MAX_PARTITION_ITEMS: usize = 4_096;
/// Maximum model-artifact depth retained in a taint explanation.
pub const MAX_TAINT_DEPTH: usize = 256;
/// Maximum UTF-8 bytes in a repartition/release justification.
pub const MAX_PARTITION_JUSTIFICATION_BYTES: usize = 4_096;
/// Maximum accepted bytes in one canonical partition receipt wire record.
///
/// The cap covers the worst currently admitted model-taint closure while
/// keeping hostile length prefixes from driving unbounded allocation.
pub const MAX_PARTITION_RECORD_BYTES: usize = 64 * 1024 * 1024;
/// fs-ledger artifact kind for canonical partition receipt wire records.
pub const PARTITION_RECEIPT_ARTIFACT_KIND: &str = "vv-partition-receipt-v1";

const ACCESS_DOMAIN: &str = "org.frankensim.fs-vvreg.dataset-access.v1";
const QUERY_CONTEXT_DOMAIN: &str = "org.frankensim.fs-vvreg.query-context.v1";
const REPARTITION_DOMAIN: &str = "org.frankensim.fs-vvreg.repartition.v1";
const BLIND_RELEASE_DOMAIN: &str = "org.frankensim.fs-vvreg.blind-release.v1";
const MODEL_TAINT_DOMAIN: &str = "org.frankensim.fs-vvreg.model-taint.v1";
const VALIDATION_DOMAIN: &str = "org.frankensim.fs-vvreg.taint-validation.v1";
const PARTITION_RECORD_MAGIC: &[u8; 8] = b"FSVVREC\0";
const PARTITION_RECORD_WIRE_VERSION: u32 = 1;

/// The semantic purpose for which a corpus row is requested.
///
/// `Training` rows deliberately share the `Calibration` purpose: learned
/// training data is calibration data for leakage accounting and therefore
/// enters the same transitive taint closure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DatasetPurpose {
    /// Fit, train, tune, or calibrate a model artifact.
    Calibration,
    /// Evaluate a frozen model on ordinary held-out validation data.
    Validation,
    /// Evaluate a frozen model on an explicitly released blind holdout.
    BlindEvaluation,
}

impl DatasetPurpose {
    /// Stable diagnostic/identity tag.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Calibration => "calibration",
            Self::Validation => "validation",
            Self::BlindEvaluation => "blind-evaluation",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Calibration => 1,
            Self::Validation => 2,
            Self::BlindEvaluation => 3,
        }
    }

    const fn admits(self, partition: DatasetPartition) -> bool {
        matches!(
            (self, partition),
            (
                Self::Calibration,
                DatasetPartition::Training | DatasetPartition::Calibration
            ) | (Self::Validation, DatasetPartition::Validation)
                | (Self::BlindEvaluation, DatasetPartition::BlindHoldout)
        )
    }
}

/// Typed refusal from the purpose, repartition, taint, or freshness boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum PartitionRefusal {
    /// The underlying evidence-bearing corpus query refused.
    Corpus(CorpusQueryRefusal),
    /// No captured partition state exists for the requested dataset.
    UnknownPartitionState {
        /// Requested dataset id.
        dataset_id: String,
    },
    /// The corpus row changed after the partition ledger was captured.
    DatasetRevisionMismatch {
        /// Changed dataset id.
        dataset_id: String,
    },
    /// The intended use is incompatible with the current partition.
    PurposeMismatch {
        /// Requested dataset id.
        dataset_id: String,
        /// Current governed partition.
        declared: DatasetPartition,
        /// Attempted use.
        attempted: DatasetPurpose,
    },
    /// A blind row has not received a generation-bound release.
    BlindReleaseRequired {
        /// Requested dataset id.
        dataset_id: String,
        /// Current partition generation.
        generation: u64,
    },
    /// A blind release was requested for a non-blind row.
    NotBlindHoldout {
        /// Requested dataset id.
        dataset_id: String,
        /// Current partition.
        partition: DatasetPartition,
    },
    /// A second, different release attempted to replace the current release.
    BlindReleaseConflict {
        /// Requested dataset id.
        dataset_id: String,
    },
    /// A previously minted access receipt no longer matches governed state.
    StaleAccess {
        /// Dataset whose receipt is stale.
        dataset_id: String,
        /// Generation bound by the receipt.
        receipt_generation: u64,
        /// Current governed generation.
        current_generation: u64,
    },
    /// A model-training input was not obtained for calibration use.
    WrongModelInputPurpose {
        /// Dataset supplied to model construction.
        dataset_id: String,
        /// Purpose carried by its sealed access receipt.
        purpose: DatasetPurpose,
    },
    /// A validation input was not obtained for validation/blind use.
    WrongEvaluationPurpose {
        /// Dataset supplied to validation.
        dataset_id: String,
        /// Purpose carried by its sealed access receipt.
        purpose: DatasetPurpose,
    },
    /// A model declaration named neither direct calibration data nor a parent
    /// with a non-empty transitive calibration closure.
    EmptyModelLineage,
    /// A model-artifact identity was the all-zero missing-value sentinel.
    ZeroArtifactIdentity,
    /// A parent chain contains the artifact being constructed.
    TaintCycle {
        /// Repeated artifact identity.
        artifact: ContentHash,
    },
    /// A bounded collection or explanation path exceeded its cap.
    ResourceLimit {
        /// Resource name.
        resource: &'static str,
        /// Maximum accepted count/length.
        limit: usize,
        /// Observed count/length.
        observed: usize,
    },
    /// Repartitioning to the already-current partition is not an event.
    RepartitionNoop {
        /// Dataset id.
        dataset_id: String,
        /// Already-current partition.
        partition: DatasetPartition,
    },
    /// Repartition/release justification was empty or oversized.
    InvalidJustification,
    /// A non-zero preregistration or blind-manifest identity was required.
    InvalidBlindIdentity {
        /// Which identity was invalid.
        field: &'static str,
    },
    /// The monotonically increasing generation counter overflowed.
    GenerationOverflow {
        /// Dataset whose generation exhausted.
        dataset_id: String,
    },
    /// Validation reused calibration data, directly or through parent models.
    TaintIntersection {
        /// Model artifact being evaluated.
        model_artifact: ContentHash,
        /// Reused dataset id.
        dataset_id: String,
        /// Reused dataset content identity.
        dataset: ContentHash,
        /// Model-artifact path from evaluated model to the direct consumer.
        model_path: Vec<ContentHash>,
    },
}

impl fmt::Display for PartitionRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Corpus(error) => write!(formatter, "corpus query refused: {error}"),
            Self::UnknownPartitionState { dataset_id } => {
                write!(
                    formatter,
                    "dataset `{dataset_id}` has no governed partition state"
                )
            }
            Self::DatasetRevisionMismatch { dataset_id } => write!(
                formatter,
                "dataset `{dataset_id}` changed after partition state was captured"
            ),
            Self::PurposeMismatch {
                dataset_id,
                declared,
                attempted,
            } => write!(
                formatter,
                "dataset `{dataset_id}` is partitioned as `{}` and refuses `{}` use",
                declared.name(),
                attempted.name()
            ),
            Self::BlindReleaseRequired {
                dataset_id,
                generation,
            } => write!(
                formatter,
                "blind dataset `{dataset_id}` generation {generation} has no release receipt"
            ),
            Self::NotBlindHoldout {
                dataset_id,
                partition,
            } => write!(
                formatter,
                "dataset `{dataset_id}` is `{}`, not a blind holdout",
                partition.name()
            ),
            Self::BlindReleaseConflict { dataset_id } => write!(
                formatter,
                "blind dataset `{dataset_id}` already has a different release receipt"
            ),
            Self::StaleAccess {
                dataset_id,
                receipt_generation,
                current_generation,
            } => write!(
                formatter,
                "dataset `{dataset_id}` access generation {receipt_generation} is stale; current generation is {current_generation}"
            ),
            Self::WrongModelInputPurpose {
                dataset_id,
                purpose,
            } => write!(
                formatter,
                "dataset `{dataset_id}` was opened for `{}` use, not model calibration",
                purpose.name()
            ),
            Self::WrongEvaluationPurpose {
                dataset_id,
                purpose,
            } => write!(
                formatter,
                "dataset `{dataset_id}` was opened for `{}` use, not held-out evaluation",
                purpose.name()
            ),
            Self::EmptyModelLineage => formatter.write_str(
                "model artifact has no calibration/training dataset in its transitive lineage",
            ),
            Self::ZeroArtifactIdentity => {
                formatter.write_str("model artifact identity must be non-zero")
            }
            Self::TaintCycle { artifact } => write!(
                formatter,
                "model taint graph repeats artifact {}",
                artifact.to_hex()
            ),
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "partition resource `{resource}` exceeds limit {limit} (observed {observed})"
            ),
            Self::RepartitionNoop {
                dataset_id,
                partition,
            } => write!(
                formatter,
                "dataset `{dataset_id}` is already partitioned as `{}`",
                partition.name()
            ),
            Self::InvalidJustification
            | Self::InvalidBlindIdentity { .. }
            | Self::GenerationOverflow { .. }
            | Self::TaintIntersection { .. } => self.fmt_blind_and_taint_tail(formatter),
        }
    }
}

impl PartitionRefusal {
    /// Render the justification, blind-release, and taint refusal family.
    fn fmt_blind_and_taint_tail(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidJustification => write!(
                formatter,
                "repartition/release justification must contain 1..={MAX_PARTITION_JUSTIFICATION_BYTES} UTF-8 bytes"
            ),
            Self::InvalidBlindIdentity { field } => {
                write!(
                    formatter,
                    "blind release `{field}` identity must be non-zero"
                )
            }
            Self::GenerationOverflow { dataset_id } => {
                write!(
                    formatter,
                    "dataset `{dataset_id}` exhausted partition generations"
                )
            }
            Self::TaintIntersection {
                model_artifact,
                dataset_id,
                dataset,
                model_path,
            } => write!(
                formatter,
                "model {} refuses validation: dataset `{dataset_id}` ({}) appears in its calibration taint through {} model artifact(s)",
                model_artifact.to_hex(),
                dataset.to_hex(),
                model_path.len()
            ),
            _ => Ok(()),
        }
    }
}

impl std::error::Error for PartitionRefusal {}

impl From<CorpusQueryRefusal> for PartitionRefusal {
    fn from(error: CorpusQueryRefusal) -> Self {
        Self::Corpus(error)
    }
}

/// Fail-closed decoding error for a canonical partition receipt wire record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PartitionRecordError {
    /// The record does not begin with the exact fs-vvreg partition magic.
    InvalidMagic,
    /// The record ended before a declared fixed or bounded field completed.
    Truncated {
        /// Field being decoded.
        field: &'static str,
    },
    /// The wire version is not implemented by this build.
    UnsupportedVersion {
        /// Observed future or historic version.
        observed: u32,
    },
    /// The closed record variant tag is unknown.
    UnknownVariant {
        /// Observed variant tag.
        observed: u8,
    },
    /// A closed enum or boolean field used an unknown tag.
    UnknownTag {
        /// Field being decoded.
        field: &'static str,
        /// Observed tag.
        observed: u8,
    },
    /// A bounded byte or item count exceeded its explicit cap.
    ResourceLimit {
        /// Resource whose prefix exceeded the cap.
        resource: &'static str,
        /// Maximum accepted count or byte length.
        limit: usize,
        /// Count or byte length declared by the record.
        observed: u64,
    },
    /// A length-framed text field was not valid UTF-8.
    InvalidUtf8 {
        /// Field being decoded.
        field: &'static str,
    },
    /// A structurally impossible semantic field was encoded.
    InvalidValue {
        /// Field or invariant that failed.
        field: &'static str,
    },
    /// The embedded semantic identity did not match the decoded fields.
    IdentityMismatch {
        /// Record variant whose identity failed.
        record: &'static str,
    },
    /// Bytes followed the exact closed record.
    TrailingBytes {
        /// Number of unconsumed bytes.
        observed: usize,
    },
    /// The parsed record did not reproduce the exact input bytes.
    NonCanonical,
}

impl fmt::Display for PartitionRecordError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMagic => formatter.write_str("partition receipt has invalid wire magic"),
            Self::Truncated { field } => {
                write!(formatter, "partition receipt is truncated in `{field}`")
            }
            Self::UnsupportedVersion { observed } => write!(
                formatter,
                "partition receipt wire version {observed} is unsupported"
            ),
            Self::UnknownVariant { observed } => {
                write!(
                    formatter,
                    "partition receipt variant tag {observed} is unknown"
                )
            }
            Self::UnknownTag { field, observed } => write!(
                formatter,
                "partition receipt field `{field}` has unknown tag {observed}"
            ),
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "partition receipt resource `{resource}` exceeds limit {limit} (observed {observed})"
            ),
            Self::InvalidUtf8 { field } => {
                write!(formatter, "partition receipt field `{field}` is not UTF-8")
            }
            Self::InvalidValue { field } => {
                write!(formatter, "partition receipt field `{field}` is invalid")
            }
            Self::IdentityMismatch { record } => {
                write!(
                    formatter,
                    "partition `{record}` receipt identity does not match its fields"
                )
            }
            Self::TrailingBytes { observed } => write!(
                formatter,
                "partition receipt has {observed} trailing byte(s)"
            ),
            Self::NonCanonical => {
                formatter.write_str("partition receipt wire bytes are not canonical")
            }
        }
    }
}

impl std::error::Error for PartitionRecordError {}

/// Immutable receipt for one purpose-checked dataset access.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetAccessReceipt {
    dataset_id: String,
    dataset: ContentHash,
    generation: u64,
    partition: DatasetPartition,
    purpose: DatasetPurpose,
    context: ContentHash,
    preceding_event: Option<ContentHash>,
    blind_release: Option<ContentHash>,
    identity: ContentHash,
}

impl DatasetAccessReceipt {
    /// Dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Exact dataset content identity.
    #[must_use]
    pub const fn dataset(&self) -> ContentHash {
        self.dataset
    }

    /// Partition generation observed at access.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Governed partition observed at access.
    #[must_use]
    pub const fn partition(&self) -> DatasetPartition {
        self.partition
    }

    /// Intended use checked at access.
    #[must_use]
    pub const fn purpose(&self) -> DatasetPurpose {
        self.purpose
    }

    /// Canonical identity of the exact, order-independent query context.
    #[must_use]
    pub const fn context(&self) -> ContentHash {
        self.context
    }

    /// Repartition event immediately preceding this access, if any.
    #[must_use]
    pub const fn preceding_event(&self) -> Option<ContentHash> {
        self.preceding_event
    }

    /// Blind release bound to this access, if blind.
    #[must_use]
    pub const fn blind_release(&self) -> Option<ContentHash> {
        self.blind_release
    }

    /// Canonical access identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Evidence-bearing dataset access plus its sealed purpose/freshness receipt.
#[derive(Debug)]
pub struct DatasetAccess<'a> {
    evidence: Evidence<&'a CorpusDataset>,
    receipt: DatasetAccessReceipt,
}

impl<'a> DatasetAccess<'a> {
    /// Underlying non-certifying corpus evidence.
    #[must_use]
    pub const fn evidence(&self) -> &Evidence<&'a CorpusDataset> {
        &self.evidence
    }

    /// Purpose/freshness receipt.
    #[must_use]
    pub const fn receipt(&self) -> &DatasetAccessReceipt {
        &self.receipt
    }
}

impl<'a> Deref for DatasetAccess<'a> {
    type Target = Evidence<&'a CorpusDataset>;

    fn deref(&self) -> &Self::Target {
        &self.evidence
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PartitionState {
    dataset: ContentHash,
    partition: DatasetPartition,
    generation: u64,
    preceding_event: Option<ContentHash>,
    blind_release: Option<BlindReleaseReceipt>,
}

/// Canonical repartition event retained by the partition ledger.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepartitionReceipt {
    dataset_id: String,
    dataset: ContentHash,
    from: DatasetPartition,
    to: DatasetPartition,
    generation: u64,
    preceding_event: Option<ContentHash>,
    justification: String,
    stales_validation_claims: bool,
    identity: ContentHash,
}

impl RepartitionReceipt {
    /// Dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Dataset content identity.
    #[must_use]
    pub const fn dataset(&self) -> ContentHash {
        self.dataset
    }

    /// Previous partition.
    #[must_use]
    pub const fn from(&self) -> DatasetPartition {
        self.from
    }

    /// New partition.
    #[must_use]
    pub const fn to(&self) -> DatasetPartition {
        self.to
    }

    /// New monotonically increasing generation.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Identity of the preceding event for this dataset.
    #[must_use]
    pub const fn preceding_event(&self) -> Option<ContentHash> {
        self.preceding_event
    }

    /// Required human/scientific justification.
    #[must_use]
    pub fn justification(&self) -> &str {
        &self.justification
    }

    /// Whether this change invalidates earlier validation/coverage claims.
    #[must_use]
    pub const fn stales_validation_claims(&self) -> bool {
        self.stales_validation_claims
    }

    /// Canonical event identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Receipt permitting one exact blind dataset generation to be evaluated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlindReleaseReceipt {
    dataset_id: String,
    dataset: ContentHash,
    generation: u64,
    preregistration: ContentHash,
    blind_manifest: ContentHash,
    justification: String,
    identity: ContentHash,
}

impl BlindReleaseReceipt {
    /// Dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Exact dataset content identity released.
    #[must_use]
    pub const fn dataset(&self) -> ContentHash {
        self.dataset
    }

    /// Partition generation released.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Preregistration identity supplied at release.
    #[must_use]
    pub const fn preregistration(&self) -> ContentHash {
        self.preregistration
    }

    /// Exact sealed blind-manifest identity supplied at release.
    #[must_use]
    pub const fn blind_manifest(&self) -> ContentHash {
        self.blind_manifest
    }

    /// Required human/scientific release justification.
    #[must_use]
    pub fn justification(&self) -> &str {
        &self.justification
    }

    /// Canonical release identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Versioned partition state and ordered canonical event receipts.
///
/// This object is deterministic state suitable for replay. It is not itself a
/// durability boundary; callers at HELM must persist returned receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PartitionLedger {
    states: BTreeMap<String, PartitionState>,
    events: Vec<RepartitionReceipt>,
    blind_releases: Vec<BlindReleaseReceipt>,
}

impl PartitionLedger {
    /// Capture the admitted partition and content identity of every corpus row.
    #[must_use]
    pub fn capture(registry: &CorpusRegistry) -> Self {
        let states = registry
            .datasets()
            .iter()
            .map(|dataset| {
                (
                    dataset.id().to_string(),
                    PartitionState {
                        dataset: dataset.digest(),
                        partition: dataset.partition(),
                        generation: 0,
                        preceding_event: None,
                        blind_release: None,
                    },
                )
            })
            .collect();
        Self {
            states,
            events: Vec::new(),
            blind_releases: Vec::new(),
        }
    }

    /// Ordered repartition event log. HELM persists these exact receipts.
    #[must_use]
    pub fn events(&self) -> &[RepartitionReceipt] {
        &self.events
    }

    /// Ordered blind-release event log. HELM persists these exact receipts.
    #[must_use]
    pub fn blind_releases(&self) -> &[BlindReleaseReceipt] {
        &self.blind_releases
    }

    /// Current governed partition and generation for one dataset.
    #[must_use]
    pub fn current(&self, dataset_id: &str) -> Option<(DatasetPartition, u64)> {
        self.states
            .get(dataset_id)
            .map(|state| (state.partition, state.generation))
    }

    /// Record a versioned repartition event.
    pub fn repartition(
        &mut self,
        dataset_id: &str,
        to: DatasetPartition,
        justification: impl Into<String>,
    ) -> Result<RepartitionReceipt, PartitionRefusal> {
        let justification = justification.into();
        validate_justification(&justification)?;
        let state = self.states.get_mut(dataset_id).ok_or_else(|| {
            PartitionRefusal::UnknownPartitionState {
                dataset_id: dataset_id.to_string(),
            }
        })?;
        if state.partition == to {
            return Err(PartitionRefusal::RepartitionNoop {
                dataset_id: dataset_id.to_string(),
                partition: to,
            });
        }
        let generation = state.generation.checked_add(1).ok_or_else(|| {
            PartitionRefusal::GenerationOverflow {
                dataset_id: dataset_id.to_string(),
            }
        })?;
        let from = state.partition;
        let stales_validation_claims = matches!(
            (from, to),
            (
                DatasetPartition::Validation | DatasetPartition::BlindHoldout,
                DatasetPartition::Training | DatasetPartition::Calibration
            )
        );
        let identity = repartition_identity(
            dataset_id,
            state.dataset,
            from,
            to,
            generation,
            state.preceding_event,
            &justification,
            stales_validation_claims,
        );
        let receipt = RepartitionReceipt {
            dataset_id: dataset_id.to_string(),
            dataset: state.dataset,
            from,
            to,
            generation,
            preceding_event: state.preceding_event,
            justification,
            stales_validation_claims,
            identity,
        };
        state.partition = to;
        state.generation = generation;
        state.preceding_event = Some(identity);
        state.blind_release = None;
        self.events.push(receipt.clone());
        Ok(receipt)
    }

    /// Release one exact blind-holdout generation.
    ///
    /// This binds non-zero preregistration and blind-manifest identities. It
    /// does not assert that either source was independently authenticated.
    pub fn release_blind(
        &mut self,
        dataset_id: &str,
        preregistration: ContentHash,
        blind_manifest: ContentHash,
        justification: impl Into<String>,
    ) -> Result<BlindReleaseReceipt, PartitionRefusal> {
        let justification = justification.into();
        validate_justification(&justification)?;
        require_nonzero(preregistration, "preregistration")?;
        require_nonzero(blind_manifest, "blind_manifest")?;
        let state = self.states.get_mut(dataset_id).ok_or_else(|| {
            PartitionRefusal::UnknownPartitionState {
                dataset_id: dataset_id.to_string(),
            }
        })?;
        if state.partition != DatasetPartition::BlindHoldout {
            return Err(PartitionRefusal::NotBlindHoldout {
                dataset_id: dataset_id.to_string(),
                partition: state.partition,
            });
        }
        let identity = blind_release_identity(
            dataset_id,
            state.dataset,
            state.generation,
            preregistration,
            blind_manifest,
            &justification,
        );
        if let Some(existing) = &state.blind_release {
            if existing.identity == identity {
                return Ok(existing.clone());
            }
            return Err(PartitionRefusal::BlindReleaseConflict {
                dataset_id: dataset_id.to_string(),
            });
        }
        let receipt = BlindReleaseReceipt {
            dataset_id: dataset_id.to_string(),
            dataset: state.dataset,
            generation: state.generation,
            preregistration,
            blind_manifest,
            justification,
            identity,
        };
        state.blind_release = Some(receipt.clone());
        self.blind_releases.push(receipt.clone());
        Ok(receipt)
    }

    fn access<'a>(
        &self,
        registry: &'a CorpusRegistry,
        dataset_id: &str,
        purpose: DatasetPurpose,
        context: &[ContextValue],
    ) -> Result<DatasetAccess<'a>, PartitionRefusal> {
        let state =
            self.states
                .get(dataset_id)
                .ok_or_else(|| PartitionRefusal::UnknownPartitionState {
                    dataset_id: dataset_id.to_string(),
                })?;
        let dataset = registry.dataset(dataset_id).ok_or_else(|| {
            PartitionRefusal::DatasetRevisionMismatch {
                dataset_id: dataset_id.to_string(),
            }
        })?;
        if dataset.digest() != state.dataset {
            return Err(PartitionRefusal::DatasetRevisionMismatch {
                dataset_id: dataset_id.to_string(),
            });
        }
        if !purpose.admits(state.partition) {
            return Err(PartitionRefusal::PurposeMismatch {
                dataset_id: dataset_id.to_string(),
                declared: state.partition,
                attempted: purpose,
            });
        }
        let blind_release = if purpose == DatasetPurpose::BlindEvaluation {
            Some(
                state
                    .blind_release
                    .as_ref()
                    .ok_or_else(|| PartitionRefusal::BlindReleaseRequired {
                        dataset_id: dataset_id.to_string(),
                        generation: state.generation,
                    })?
                    .identity,
            )
        } else {
            None
        };
        // The original declared partition is used only to cross the seeded
        // corpus-authority boundary. The governed state above owns current
        // purpose semantics and binds its generation into the returned receipt.
        let evidence =
            registry.query_declared_partition(dataset_id, dataset.partition(), context)?;
        let context = query_context_identity(context);
        let identity = access_identity(
            dataset_id,
            state.dataset,
            state.generation,
            state.partition,
            purpose,
            context,
            state.preceding_event,
            blind_release,
        );
        Ok(DatasetAccess {
            evidence,
            receipt: DatasetAccessReceipt {
                dataset_id: dataset_id.to_string(),
                dataset: state.dataset,
                generation: state.generation,
                partition: state.partition,
                purpose,
                context,
                preceding_event: state.preceding_event,
                blind_release,
                identity,
            },
        })
    }

    fn require_fresh(&self, access: &DatasetAccess<'_>) -> Result<(), PartitionRefusal> {
        let receipt = access.receipt();
        let state = self.states.get(receipt.dataset_id()).ok_or_else(|| {
            PartitionRefusal::UnknownPartitionState {
                dataset_id: receipt.dataset_id().to_string(),
            }
        })?;
        let fresh = state.dataset == receipt.dataset
            && state.generation == receipt.generation
            && state.partition == receipt.partition
            && receipt.purpose.admits(state.partition)
            && (receipt.purpose != DatasetPurpose::BlindEvaluation
                || state.blind_release.as_ref().map(|release| release.identity)
                    == receipt.blind_release);
        if fresh {
            Ok(())
        } else {
            Err(PartitionRefusal::StaleAccess {
                dataset_id: receipt.dataset_id().to_string(),
                receipt_generation: receipt.generation,
                current_generation: state.generation,
            })
        }
    }

    /// Register one model artifact from fresh calibration inputs and complete
    /// parent-model taint closures.
    ///
    /// Direct receipts are checked against this exact ledger before the model
    /// identity is minted, preventing stale or mismatched-generation
    /// calibration accesses from bypassing the taint boundary.
    pub fn register_model(
        &self,
        artifact: ContentHash,
        direct: &[&DatasetAccess<'_>],
        parents: &[&ModelTaint],
    ) -> Result<ModelTaint, PartitionRefusal> {
        for access in direct {
            self.require_fresh(access)?;
        }
        ModelTaint::build(artifact, direct, parents)
    }

    /// Check that evaluation data are fresh, held out, and disjoint from a
    /// model's complete calibration-taint closure.
    pub fn validate_model(
        &self,
        model: &ModelTaint,
        evaluation: &[&DatasetAccess<'_>],
    ) -> Result<TaintValidationReceipt, PartitionRefusal> {
        check_count("evaluation inputs", evaluation.len())?;
        if evaluation.is_empty() {
            return Err(PartitionRefusal::ResourceLimit {
                resource: "evaluation inputs",
                limit: MAX_PARTITION_ITEMS,
                observed: 0,
            });
        }
        let mut accesses = BTreeSet::new();
        for access in evaluation {
            self.require_fresh(access)?;
            let receipt = access.receipt();
            if !matches!(
                receipt.purpose,
                DatasetPurpose::Validation | DatasetPurpose::BlindEvaluation
            ) {
                return Err(PartitionRefusal::WrongEvaluationPurpose {
                    dataset_id: receipt.dataset_id.clone(),
                    purpose: receipt.purpose,
                });
            }
            if let Some(source) = model.sources.get(&receipt.dataset) {
                return Err(PartitionRefusal::TaintIntersection {
                    model_artifact: model.artifact,
                    dataset_id: source.dataset_id.clone(),
                    dataset: receipt.dataset,
                    model_path: source.model_path.clone(),
                });
            }
            accesses.insert(receipt.identity);
        }
        let identity = validation_identity(model.identity, &accesses);
        Ok(TaintValidationReceipt {
            model_taint: model.identity,
            evaluation_accesses: accesses.into_iter().collect(),
            identity,
        })
    }
}

impl CorpusRegistry {
    /// Query a dataset for one semantic purpose under explicit versioned
    /// partition state.
    pub fn query<'a>(
        &'a self,
        partitions: &PartitionLedger,
        dataset_id: &str,
        purpose: DatasetPurpose,
        context: &[ContextValue],
    ) -> Result<DatasetAccess<'a>, PartitionRefusal> {
        partitions.access(self, dataset_id, purpose, context)
    }
}

/// One dataset in a model's transitive calibration-taint closure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaintSource {
    dataset_id: String,
    dataset: ContentHash,
    access_generation: u64,
    model_path: Vec<ContentHash>,
}

impl TaintSource {
    /// Dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Dataset content identity.
    #[must_use]
    pub const fn dataset(&self) -> ContentHash {
        self.dataset
    }

    /// Partition generation at calibration access.
    #[must_use]
    pub const fn access_generation(&self) -> u64 {
        self.access_generation
    }

    /// Path from the outer model to the model that directly consumed data.
    #[must_use]
    pub fn model_path(&self) -> &[ContentHash] {
        &self.model_path
    }
}

/// Exact transitive calibration-taint closure for one model artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelTaint {
    artifact: ContentHash,
    sources: BTreeMap<ContentHash, TaintSource>,
    identity: ContentHash,
}

impl ModelTaint {
    fn build(
        artifact: ContentHash,
        direct: &[&DatasetAccess<'_>],
        parents: &[&Self],
    ) -> Result<Self, PartitionRefusal> {
        require_nonzero_artifact(artifact)?;
        check_count("direct calibration inputs", direct.len())?;
        check_count("parent model inputs", parents.len())?;
        let mut sources = BTreeMap::new();
        for access in direct {
            let receipt = access.receipt();
            if receipt.purpose != DatasetPurpose::Calibration {
                return Err(PartitionRefusal::WrongModelInputPurpose {
                    dataset_id: receipt.dataset_id.clone(),
                    purpose: receipt.purpose,
                });
            }
            insert_source(
                &mut sources,
                TaintSource {
                    dataset_id: receipt.dataset_id.clone(),
                    dataset: receipt.dataset,
                    access_generation: receipt.generation,
                    model_path: vec![artifact],
                },
            );
        }
        for parent in parents {
            if parent.artifact == artifact
                || parent
                    .sources
                    .values()
                    .any(|source| source.model_path.contains(&artifact))
            {
                return Err(PartitionRefusal::TaintCycle { artifact });
            }
            for source in parent.sources.values() {
                let observed = source.model_path.len().saturating_add(1);
                if observed > MAX_TAINT_DEPTH {
                    return Err(PartitionRefusal::ResourceLimit {
                        resource: "taint depth",
                        limit: MAX_TAINT_DEPTH,
                        observed,
                    });
                }
                let mut model_path = Vec::with_capacity(observed);
                model_path.push(artifact);
                model_path.extend_from_slice(&source.model_path);
                insert_source(
                    &mut sources,
                    TaintSource {
                        dataset_id: source.dataset_id.clone(),
                        dataset: source.dataset,
                        access_generation: source.access_generation,
                        model_path,
                    },
                );
            }
        }
        if sources.is_empty() {
            return Err(PartitionRefusal::EmptyModelLineage);
        }
        check_count("transitive taint sources", sources.len())?;
        let identity = model_taint_identity(artifact, &sources);
        Ok(Self {
            artifact,
            sources,
            identity,
        })
    }

    /// Model artifact identity.
    #[must_use]
    pub const fn artifact(&self) -> ContentHash {
        self.artifact
    }

    /// Canonically sorted transitive calibration sources.
    #[must_use]
    pub fn sources(&self) -> impl ExactSizeIterator<Item = &TaintSource> {
        self.sources.values()
    }

    /// Canonical closure identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Receipt proving that one exact model-taint closure was disjoint from a
/// bounded set of fresh evaluation-access receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaintValidationReceipt {
    model_taint: ContentHash,
    evaluation_accesses: Vec<ContentHash>,
    identity: ContentHash,
}

impl TaintValidationReceipt {
    /// Exact model-taint identity checked.
    #[must_use]
    pub const fn model_taint(&self) -> ContentHash {
        self.model_taint
    }

    /// Sorted, deduplicated evaluation-access identities.
    #[must_use]
    pub fn evaluation_accesses(&self) -> &[ContentHash] {
        &self.evaluation_accesses
    }

    /// Canonical validation-check identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Closed wire envelope for one purpose, partition, taint, or validation
/// receipt.
///
/// The embedded receipt identity is re-derived during decode. The envelope is
/// suitable for content-addressed storage, but is not an authorization,
/// ledger-membership, or audit-event receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PartitionReceiptRecord {
    /// One purpose-checked dataset access.
    DatasetAccess(DatasetAccessReceipt),
    /// One governed partition transition.
    Repartition(RepartitionReceipt),
    /// One blind-holdout release.
    BlindRelease(BlindReleaseReceipt),
    /// One model's complete calibration-taint closure.
    ModelTaint(ModelTaint),
    /// One successful disjoint held-out validation check.
    Validation(TaintValidationReceipt),
}

impl PartitionReceiptRecord {
    /// Stable lowercase record variant name for diagnostics and metadata.
    #[must_use]
    pub const fn record_type(&self) -> &'static str {
        match self {
            Self::DatasetAccess(_) => "dataset-access",
            Self::Repartition(_) => "repartition",
            Self::BlindRelease(_) => "blind-release",
            Self::ModelTaint(_) => "model-taint",
            Self::Validation(_) => "validation",
        }
    }

    /// Existing semantic identity embedded in this transport record.
    #[must_use]
    pub const fn semantic_identity(&self) -> ContentHash {
        match self {
            Self::DatasetAccess(receipt) => receipt.identity(),
            Self::Repartition(receipt) => receipt.identity(),
            Self::BlindRelease(receipt) => receipt.identity(),
            Self::ModelTaint(receipt) => receipt.identity(),
            Self::Validation(receipt) => receipt.identity(),
        }
    }

    /// Encode this record in canonical partition receipt wire format v1.
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(PARTITION_RECORD_MAGIC);
        bytes.extend_from_slice(&PARTITION_RECORD_WIRE_VERSION.to_le_bytes());
        match self {
            Self::DatasetAccess(receipt) => {
                bytes.push(1);
                push_record_text(&mut bytes, &receipt.dataset_id);
                push_hash(&mut bytes, receipt.dataset);
                push_u64(&mut bytes, receipt.generation);
                bytes.push(partition_tag(receipt.partition));
                bytes.push(receipt.purpose.tag());
                push_hash(&mut bytes, receipt.context);
                push_optional_hash(&mut bytes, receipt.preceding_event);
                push_optional_hash(&mut bytes, receipt.blind_release);
                push_hash(&mut bytes, receipt.identity);
            }
            Self::Repartition(receipt) => {
                bytes.push(2);
                push_record_text(&mut bytes, &receipt.dataset_id);
                push_hash(&mut bytes, receipt.dataset);
                bytes.push(partition_tag(receipt.from));
                bytes.push(partition_tag(receipt.to));
                push_u64(&mut bytes, receipt.generation);
                push_optional_hash(&mut bytes, receipt.preceding_event);
                push_record_text(&mut bytes, &receipt.justification);
                bytes.push(u8::from(receipt.stales_validation_claims));
                push_hash(&mut bytes, receipt.identity);
            }
            Self::BlindRelease(receipt) => {
                bytes.push(3);
                push_record_text(&mut bytes, &receipt.dataset_id);
                push_hash(&mut bytes, receipt.dataset);
                push_u64(&mut bytes, receipt.generation);
                push_hash(&mut bytes, receipt.preregistration);
                push_hash(&mut bytes, receipt.blind_manifest);
                push_record_text(&mut bytes, &receipt.justification);
                push_hash(&mut bytes, receipt.identity);
            }
            Self::ModelTaint(receipt) => {
                bytes.push(4);
                push_hash(&mut bytes, receipt.artifact);
                push_u64(&mut bytes, receipt.sources.len() as u64);
                for source in receipt.sources.values() {
                    push_record_text(&mut bytes, &source.dataset_id);
                    push_hash(&mut bytes, source.dataset);
                    push_u64(&mut bytes, source.access_generation);
                    push_u64(&mut bytes, source.model_path.len() as u64);
                    for artifact in &source.model_path {
                        push_hash(&mut bytes, *artifact);
                    }
                }
                push_hash(&mut bytes, receipt.identity);
            }
            Self::Validation(receipt) => {
                bytes.push(5);
                push_hash(&mut bytes, receipt.model_taint);
                push_u64(&mut bytes, receipt.evaluation_accesses.len() as u64);
                for access in &receipt.evaluation_accesses {
                    push_hash(&mut bytes, *access);
                }
                push_hash(&mut bytes, receipt.identity);
            }
        }
        bytes
    }

    /// Decode, structurally validate, and independently re-derive one
    /// canonical partition receipt record.
    pub fn decode(bytes: &[u8]) -> Result<Self, PartitionRecordError> {
        let observed = u64::try_from(bytes.len()).unwrap_or(u64::MAX);
        if bytes.len() > MAX_PARTITION_RECORD_BYTES {
            return Err(PartitionRecordError::ResourceLimit {
                resource: "record bytes",
                limit: MAX_PARTITION_RECORD_BYTES,
                observed,
            });
        }
        let mut reader = PartitionRecordReader::new(bytes);
        if reader.take(PARTITION_RECORD_MAGIC.len(), "wire magic")? != PARTITION_RECORD_MAGIC {
            return Err(PartitionRecordError::InvalidMagic);
        }
        let version = reader.u32("wire version")?;
        if version != PARTITION_RECORD_WIRE_VERSION {
            return Err(PartitionRecordError::UnsupportedVersion { observed: version });
        }
        let record = match reader.byte("record variant")? {
            1 => Self::DatasetAccess(decode_access_receipt(&mut reader)?),
            2 => Self::Repartition(decode_repartition_receipt(&mut reader)?),
            3 => Self::BlindRelease(decode_blind_release_receipt(&mut reader)?),
            4 => Self::ModelTaint(decode_model_taint(&mut reader)?),
            5 => Self::Validation(decode_validation_receipt(&mut reader)?),
            observed => return Err(PartitionRecordError::UnknownVariant { observed }),
        };
        let trailing = reader.remaining();
        if trailing != 0 {
            return Err(PartitionRecordError::TrailingBytes { observed: trailing });
        }
        if record.encode() != bytes {
            return Err(PartitionRecordError::NonCanonical);
        }
        Ok(record)
    }
}

impl From<DatasetAccessReceipt> for PartitionReceiptRecord {
    fn from(receipt: DatasetAccessReceipt) -> Self {
        Self::DatasetAccess(receipt)
    }
}

impl From<RepartitionReceipt> for PartitionReceiptRecord {
    fn from(receipt: RepartitionReceipt) -> Self {
        Self::Repartition(receipt)
    }
}

impl From<BlindReleaseReceipt> for PartitionReceiptRecord {
    fn from(receipt: BlindReleaseReceipt) -> Self {
        Self::BlindRelease(receipt)
    }
}

impl From<ModelTaint> for PartitionReceiptRecord {
    fn from(receipt: ModelTaint) -> Self {
        Self::ModelTaint(receipt)
    }
}

impl From<TaintValidationReceipt> for PartitionReceiptRecord {
    fn from(receipt: TaintValidationReceipt) -> Self {
        Self::Validation(receipt)
    }
}

struct PartitionRecordReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> PartitionRecordReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.offset)
    }

    fn take(&mut self, len: usize, field: &'static str) -> Result<&'a [u8], PartitionRecordError> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or(PartitionRecordError::Truncated { field })?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(PartitionRecordError::Truncated { field })?;
        self.offset = end;
        Ok(value)
    }

    fn byte(&mut self, field: &'static str) -> Result<u8, PartitionRecordError> {
        Ok(self.take(1, field)?[0])
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, PartitionRecordError> {
        let bytes: [u8; 4] = self
            .take(4, field)?
            .try_into()
            .map_err(|_| PartitionRecordError::Truncated { field })?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self, field: &'static str) -> Result<u64, PartitionRecordError> {
        let bytes: [u8; 8] = self
            .take(8, field)?
            .try_into()
            .map_err(|_| PartitionRecordError::Truncated { field })?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn count(
        &mut self,
        resource: &'static str,
        limit: usize,
    ) -> Result<usize, PartitionRecordError> {
        let observed = self.u64(resource)?;
        let limit_u64 = u64::try_from(limit).unwrap_or(u64::MAX);
        if observed > limit_u64 {
            return Err(PartitionRecordError::ResourceLimit {
                resource,
                limit,
                observed,
            });
        }
        usize::try_from(observed).map_err(|_| PartitionRecordError::ResourceLimit {
            resource,
            limit,
            observed,
        })
    }

    fn text(&mut self, field: &'static str, limit: usize) -> Result<String, PartitionRecordError> {
        let len = self.count(field, limit)?;
        let bytes = self.take(len, field)?;
        let text =
            std::str::from_utf8(bytes).map_err(|_| PartitionRecordError::InvalidUtf8 { field })?;
        Ok(text.to_owned())
    }

    fn hash(&mut self, field: &'static str) -> Result<ContentHash, PartitionRecordError> {
        let bytes: [u8; 32] = self
            .take(32, field)?
            .try_into()
            .map_err(|_| PartitionRecordError::Truncated { field })?;
        Ok(ContentHash(bytes))
    }

    fn optional_hash(
        &mut self,
        field: &'static str,
    ) -> Result<Option<ContentHash>, PartitionRecordError> {
        match self.byte(field)? {
            0 => Ok(None),
            1 => self.hash(field).map(Some),
            observed => Err(PartitionRecordError::UnknownTag { field, observed }),
        }
    }

    fn partition(&mut self, field: &'static str) -> Result<DatasetPartition, PartitionRecordError> {
        match self.byte(field)? {
            1 => Ok(DatasetPartition::Training),
            2 => Ok(DatasetPartition::Calibration),
            3 => Ok(DatasetPartition::Validation),
            4 => Ok(DatasetPartition::BlindHoldout),
            observed => Err(PartitionRecordError::UnknownTag { field, observed }),
        }
    }

    fn purpose(&mut self, field: &'static str) -> Result<DatasetPurpose, PartitionRecordError> {
        match self.byte(field)? {
            1 => Ok(DatasetPurpose::Calibration),
            2 => Ok(DatasetPurpose::Validation),
            3 => Ok(DatasetPurpose::BlindEvaluation),
            observed => Err(PartitionRecordError::UnknownTag { field, observed }),
        }
    }

    fn boolean(&mut self, field: &'static str) -> Result<bool, PartitionRecordError> {
        match self.byte(field)? {
            0 => Ok(false),
            1 => Ok(true),
            observed => Err(PartitionRecordError::UnknownTag { field, observed }),
        }
    }
}

fn decode_access_receipt(
    reader: &mut PartitionRecordReader<'_>,
) -> Result<DatasetAccessReceipt, PartitionRecordError> {
    let dataset_id = reader.text("dataset id", MAX_CORPUS_TEXT_BYTES)?;
    validate_record_dataset_id(&dataset_id)?;
    let dataset = reader.hash("dataset identity")?;
    let generation = reader.u64("partition generation")?;
    let partition = reader.partition("dataset partition")?;
    let purpose = reader.purpose("dataset purpose")?;
    let context = reader.hash("query context")?;
    let preceding_event = reader.optional_hash("preceding repartition event")?;
    let blind_release = reader.optional_hash("blind release")?;
    let identity = reader.hash("access identity")?;
    if !purpose.admits(partition)
        || (purpose == DatasetPurpose::BlindEvaluation) != blind_release.is_some()
    {
        return Err(PartitionRecordError::InvalidValue {
            field: "purpose/partition/blind-release relation",
        });
    }
    let expected = access_identity(
        &dataset_id,
        dataset,
        generation,
        partition,
        purpose,
        context,
        preceding_event,
        blind_release,
    );
    require_record_identity("dataset-access", identity, expected)?;
    Ok(DatasetAccessReceipt {
        dataset_id,
        dataset,
        generation,
        partition,
        purpose,
        context,
        preceding_event,
        blind_release,
        identity,
    })
}

fn decode_repartition_receipt(
    reader: &mut PartitionRecordReader<'_>,
) -> Result<RepartitionReceipt, PartitionRecordError> {
    let dataset_id = reader.text("dataset id", MAX_CORPUS_TEXT_BYTES)?;
    validate_record_dataset_id(&dataset_id)?;
    let dataset = reader.hash("dataset identity")?;
    let from = reader.partition("source partition")?;
    let to = reader.partition("target partition")?;
    let generation = reader.u64("partition generation")?;
    let preceding_event = reader.optional_hash("preceding repartition event")?;
    let justification = reader.text(
        "repartition justification",
        MAX_PARTITION_JUSTIFICATION_BYTES,
    )?;
    let stales_validation_claims = reader.boolean("stales validation claims")?;
    let identity = reader.hash("repartition identity")?;
    let expected_staleness = matches!(
        (from, to),
        (
            DatasetPartition::Validation | DatasetPartition::BlindHoldout,
            DatasetPartition::Training | DatasetPartition::Calibration
        )
    );
    if from == to
        || generation == 0
        || (generation == 1) != preceding_event.is_none()
        || stales_validation_claims != expected_staleness
        || validate_justification(&justification).is_err()
    {
        return Err(PartitionRecordError::InvalidValue {
            field: "repartition transition",
        });
    }
    let expected = repartition_identity(
        &dataset_id,
        dataset,
        from,
        to,
        generation,
        preceding_event,
        &justification,
        stales_validation_claims,
    );
    require_record_identity("repartition", identity, expected)?;
    Ok(RepartitionReceipt {
        dataset_id,
        dataset,
        from,
        to,
        generation,
        preceding_event,
        justification,
        stales_validation_claims,
        identity,
    })
}

fn decode_blind_release_receipt(
    reader: &mut PartitionRecordReader<'_>,
) -> Result<BlindReleaseReceipt, PartitionRecordError> {
    let dataset_id = reader.text("dataset id", MAX_CORPUS_TEXT_BYTES)?;
    validate_record_dataset_id(&dataset_id)?;
    let dataset = reader.hash("dataset identity")?;
    let generation = reader.u64("partition generation")?;
    let preregistration = reader.hash("preregistration identity")?;
    let blind_manifest = reader.hash("blind manifest identity")?;
    let justification = reader.text(
        "blind release justification",
        MAX_PARTITION_JUSTIFICATION_BYTES,
    )?;
    let identity = reader.hash("blind release identity")?;
    if require_nonzero(preregistration, "preregistration").is_err()
        || require_nonzero(blind_manifest, "blind_manifest").is_err()
        || validate_justification(&justification).is_err()
    {
        return Err(PartitionRecordError::InvalidValue {
            field: "blind release",
        });
    }
    let expected = blind_release_identity(
        &dataset_id,
        dataset,
        generation,
        preregistration,
        blind_manifest,
        &justification,
    );
    require_record_identity("blind-release", identity, expected)?;
    Ok(BlindReleaseReceipt {
        dataset_id,
        dataset,
        generation,
        preregistration,
        blind_manifest,
        justification,
        identity,
    })
}

fn decode_model_taint(
    reader: &mut PartitionRecordReader<'_>,
) -> Result<ModelTaint, PartitionRecordError> {
    let artifact = reader.hash("model artifact")?;
    if require_nonzero_artifact(artifact).is_err() {
        return Err(PartitionRecordError::InvalidValue {
            field: "model artifact",
        });
    }
    let source_count = reader.count("taint sources", MAX_PARTITION_ITEMS)?;
    if source_count == 0 {
        return Err(PartitionRecordError::InvalidValue {
            field: "taint sources",
        });
    }
    let mut sources = BTreeMap::new();
    let mut preceding_dataset = None;
    for _ in 0..source_count {
        let dataset_id = reader.text("taint dataset id", MAX_CORPUS_TEXT_BYTES)?;
        validate_record_dataset_id(&dataset_id)?;
        let dataset = reader.hash("taint dataset identity")?;
        if preceding_dataset.is_some_and(|preceding| preceding >= dataset) {
            return Err(PartitionRecordError::InvalidValue {
                field: "taint source ordering",
            });
        }
        preceding_dataset = Some(dataset);
        let access_generation = reader.u64("taint access generation")?;
        let path_count = reader.count("taint model path", MAX_TAINT_DEPTH)?;
        if path_count == 0 {
            return Err(PartitionRecordError::InvalidValue {
                field: "taint model path",
            });
        }
        let mut model_path = Vec::with_capacity(path_count);
        for _ in 0..path_count {
            model_path.push(reader.hash("taint model path artifact")?);
        }
        if model_path.first() != Some(&artifact) {
            return Err(PartitionRecordError::InvalidValue {
                field: "taint model path root",
            });
        }
        sources.insert(
            dataset,
            TaintSource {
                dataset_id,
                dataset,
                access_generation,
                model_path,
            },
        );
    }
    let identity = reader.hash("model taint identity")?;
    let expected = model_taint_identity(artifact, &sources);
    require_record_identity("model-taint", identity, expected)?;
    Ok(ModelTaint {
        artifact,
        sources,
        identity,
    })
}

fn decode_validation_receipt(
    reader: &mut PartitionRecordReader<'_>,
) -> Result<TaintValidationReceipt, PartitionRecordError> {
    let model_taint = reader.hash("model taint identity")?;
    let access_count = reader.count("validation accesses", MAX_PARTITION_ITEMS)?;
    if access_count == 0 {
        return Err(PartitionRecordError::InvalidValue {
            field: "validation accesses",
        });
    }
    let mut accesses = BTreeSet::new();
    let mut evaluation_accesses = Vec::with_capacity(access_count);
    let mut preceding_access = None;
    for _ in 0..access_count {
        let access = reader.hash("validation access identity")?;
        if preceding_access.is_some_and(|preceding| preceding >= access) {
            return Err(PartitionRecordError::InvalidValue {
                field: "validation access ordering",
            });
        }
        preceding_access = Some(access);
        accesses.insert(access);
        evaluation_accesses.push(access);
    }
    let identity = reader.hash("validation identity")?;
    let expected = validation_identity(model_taint, &accesses);
    require_record_identity("validation", identity, expected)?;
    Ok(TaintValidationReceipt {
        model_taint,
        evaluation_accesses,
        identity,
    })
}

fn validate_record_dataset_id(dataset_id: &str) -> Result<(), PartitionRecordError> {
    let valid = !dataset_id.is_empty()
        && dataset_id.len() <= MAX_CORPUS_TEXT_BYTES
        && dataset_id.as_bytes()[0].is_ascii_lowercase()
        && dataset_id.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-' || byte == b'_'
        });
    if valid {
        Ok(())
    } else {
        Err(PartitionRecordError::InvalidValue {
            field: "dataset id",
        })
    }
}

fn require_record_identity(
    record: &'static str,
    observed: ContentHash,
    expected: ContentHash,
) -> Result<(), PartitionRecordError> {
    if observed == expected {
        Ok(())
    } else {
        Err(PartitionRecordError::IdentityMismatch { record })
    }
}

fn push_record_text(bytes: &mut Vec<u8>, value: &str) {
    push_u64(bytes, value.len() as u64);
    bytes.extend_from_slice(value.as_bytes());
}

fn insert_source(sources: &mut BTreeMap<ContentHash, TaintSource>, candidate: TaintSource) {
    match sources.get(&candidate.dataset) {
        Some(current) if source_order_key(current) <= source_order_key(&candidate) => {}
        _ => {
            sources.insert(candidate.dataset, candidate);
        }
    }
}

fn source_order_key(source: &TaintSource) -> (&str, u64, &[ContentHash]) {
    (
        source.dataset_id.as_str(),
        source.access_generation,
        source.model_path.as_slice(),
    )
}

fn validate_justification(justification: &str) -> Result<(), PartitionRefusal> {
    if justification.trim().is_empty() || justification.len() > MAX_PARTITION_JUSTIFICATION_BYTES {
        Err(PartitionRefusal::InvalidJustification)
    } else {
        Ok(())
    }
}

fn require_nonzero(value: ContentHash, field: &'static str) -> Result<(), PartitionRefusal> {
    if value.as_bytes().iter().all(|byte| *byte == 0) {
        Err(PartitionRefusal::InvalidBlindIdentity { field })
    } else {
        Ok(())
    }
}

fn require_nonzero_artifact(value: ContentHash) -> Result<(), PartitionRefusal> {
    if value.as_bytes().iter().all(|byte| *byte == 0) {
        Err(PartitionRefusal::ZeroArtifactIdentity)
    } else {
        Ok(())
    }
}

fn check_count(resource: &'static str, observed: usize) -> Result<(), PartitionRefusal> {
    if observed > MAX_PARTITION_ITEMS {
        Err(PartitionRefusal::ResourceLimit {
            resource,
            limit: MAX_PARTITION_ITEMS,
            observed,
        })
    } else {
        Ok(())
    }
}

const fn partition_tag(partition: DatasetPartition) -> u8 {
    match partition {
        DatasetPartition::Training => 1,
        DatasetPartition::Calibration => 2,
        DatasetPartition::Validation => 3,
        DatasetPartition::BlindHoldout => 4,
    }
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_u64(bytes, value.len() as u64);
    bytes.extend_from_slice(value.as_bytes());
}

fn push_hash(bytes: &mut Vec<u8>, value: ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}

fn push_optional_hash(bytes: &mut Vec<u8>, value: Option<ContentHash>) {
    match value {
        Some(value) => {
            bytes.push(1);
            push_hash(bytes, value);
        }
        None => bytes.push(0),
    }
}

// The eight ordered fields ARE the canonical preimage; bundling them into
// a struct would hide the domain-separation order this hash commits to.
// Same ruling as repartition_identity below.
#[allow(clippy::too_many_arguments)]
fn access_identity(
    dataset_id: &str,
    dataset: ContentHash,
    generation: u64,
    partition: DatasetPartition,
    purpose: DatasetPurpose,
    context: ContentHash,
    preceding_event: Option<ContentHash>,
    blind_release: Option<ContentHash>,
) -> ContentHash {
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_text(&mut bytes, dataset_id);
    push_hash(&mut bytes, dataset);
    push_u64(&mut bytes, generation);
    bytes.push(partition_tag(partition));
    bytes.push(purpose.tag());
    push_hash(&mut bytes, context);
    push_optional_hash(&mut bytes, preceding_event);
    push_optional_hash(&mut bytes, blind_release);
    hash_domain(ACCESS_DOMAIN, &bytes)
}

fn query_context_identity(context: &[ContextValue]) -> ContentHash {
    let mut ordered = context.iter().collect::<Vec<_>>();
    ordered.sort_by(|left, right| left.name.cmp(&right.name));
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_u64(&mut bytes, ordered.len() as u64);
    for coordinate in ordered {
        push_text(&mut bytes, &coordinate.name);
        bytes.extend_from_slice(&coordinate.value.value.to_bits().to_le_bytes());
        for exponent in coordinate.value.dims.0 {
            bytes.push(exponent as u8);
        }
    }
    hash_domain(QUERY_CONTEXT_DOMAIN, &bytes)
}

#[allow(clippy::too_many_arguments)]
fn repartition_identity(
    dataset_id: &str,
    dataset: ContentHash,
    from: DatasetPartition,
    to: DatasetPartition,
    generation: u64,
    preceding_event: Option<ContentHash>,
    justification: &str,
    stales_validation_claims: bool,
) -> ContentHash {
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_text(&mut bytes, dataset_id);
    push_hash(&mut bytes, dataset);
    bytes.push(partition_tag(from));
    bytes.push(partition_tag(to));
    push_u64(&mut bytes, generation);
    push_optional_hash(&mut bytes, preceding_event);
    push_text(&mut bytes, justification);
    bytes.push(u8::from(stales_validation_claims));
    hash_domain(REPARTITION_DOMAIN, &bytes)
}

fn blind_release_identity(
    dataset_id: &str,
    dataset: ContentHash,
    generation: u64,
    preregistration: ContentHash,
    blind_manifest: ContentHash,
    justification: &str,
) -> ContentHash {
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_text(&mut bytes, dataset_id);
    push_hash(&mut bytes, dataset);
    push_u64(&mut bytes, generation);
    push_hash(&mut bytes, preregistration);
    push_hash(&mut bytes, blind_manifest);
    push_text(&mut bytes, justification);
    hash_domain(BLIND_RELEASE_DOMAIN, &bytes)
}

fn model_taint_identity(
    artifact: ContentHash,
    sources: &BTreeMap<ContentHash, TaintSource>,
) -> ContentHash {
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_hash(&mut bytes, artifact);
    push_u64(&mut bytes, sources.len() as u64);
    for source in sources.values() {
        push_text(&mut bytes, &source.dataset_id);
        push_hash(&mut bytes, source.dataset);
        push_u64(&mut bytes, source.access_generation);
        push_u64(&mut bytes, source.model_path.len() as u64);
        for artifact in &source.model_path {
            push_hash(&mut bytes, *artifact);
        }
    }
    hash_domain(MODEL_TAINT_DOMAIN, &bytes)
}

fn validation_identity(model_taint: ContentHash, accesses: &BTreeSet<ContentHash>) -> ContentHash {
    let mut bytes = Vec::new();
    push_u64(&mut bytes, u64::from(PARTITION_POLICY_SCHEMA_VERSION));
    push_hash(&mut bytes, model_taint);
    push_u64(&mut bytes, accesses.len() as u64);
    for access in accesses {
        push_hash(&mut bytes, *access);
    }
    hash_domain(VALIDATION_DOMAIN, &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn purpose_matrix_is_exact() {
        assert!(DatasetPurpose::Calibration.admits(DatasetPartition::Training));
        assert!(DatasetPurpose::Calibration.admits(DatasetPartition::Calibration));
        assert!(DatasetPurpose::Validation.admits(DatasetPartition::Validation));
        assert!(DatasetPurpose::BlindEvaluation.admits(DatasetPartition::BlindHoldout));
        for (purpose, partition) in [
            (DatasetPurpose::Calibration, DatasetPartition::Validation),
            (DatasetPurpose::Calibration, DatasetPartition::BlindHoldout),
            (DatasetPurpose::Validation, DatasetPartition::Training),
            (DatasetPurpose::Validation, DatasetPartition::Calibration),
            (DatasetPurpose::Validation, DatasetPartition::BlindHoldout),
            (DatasetPurpose::BlindEvaluation, DatasetPartition::Training),
            (
                DatasetPurpose::BlindEvaluation,
                DatasetPartition::Calibration,
            ),
            (
                DatasetPurpose::BlindEvaluation,
                DatasetPartition::Validation,
            ),
        ] {
            assert!(
                !purpose.admits(partition),
                "{purpose:?} admitted {partition:?}"
            );
        }
    }

    #[test]
    fn source_tie_break_is_order_independent() {
        let dataset = hash_domain("test-dataset", b"x");
        let a = hash_domain("test-model", b"a");
        let b = hash_domain("test-model", b"b");
        let mut forward = BTreeMap::new();
        insert_source(
            &mut forward,
            TaintSource {
                dataset_id: "x".to_string(),
                dataset,
                access_generation: 1,
                model_path: vec![b],
            },
        );
        insert_source(
            &mut forward,
            TaintSource {
                dataset_id: "x".to_string(),
                dataset,
                access_generation: 1,
                model_path: vec![a],
            },
        );
        let mut reverse = BTreeMap::new();
        for source in forward.values().cloned().chain([TaintSource {
            dataset_id: "x".to_string(),
            dataset,
            access_generation: 1,
            model_path: vec![b],
        }]) {
            insert_source(&mut reverse, source);
        }
        assert_eq!(forward, reverse);
        assert_eq!(forward[&dataset].model_path, vec![a.min(b)]);
    }

    #[test]
    fn query_context_identity_is_order_independent_and_value_sensitive() {
        let a = ContextValue {
            name: "a".to_string(),
            value: fs_qty::QtyAny::dimensionless(1.0),
        };
        let b = ContextValue {
            name: "b".to_string(),
            value: fs_qty::QtyAny::new(2.0, fs_qty::Dims([1, 0, 0, 0, 0, 0])),
        };
        assert_eq!(
            query_context_identity(&[a.clone(), b.clone()]),
            query_context_identity(&[b.clone(), a.clone()])
        );
        let changed = ContextValue {
            value: fs_qty::QtyAny::new(3.0, b.value.dims),
            ..b.clone()
        };
        assert_ne!(
            query_context_identity(&[a.clone(), changed]),
            query_context_identity(&[a, b])
        );
    }
}

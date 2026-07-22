//! Evidence-bearing validation corpus schema (EXTREAL E04, bead
//! `frankensim-extreal-program-f85xj.4.1`).
//!
//! Dataset admission is fail-closed and serialization is canonical.  The
//! seeded registry can return [`Evidence`] only after an exact partition and
//! context-of-use check.  The wrapper deliberately carries a numerical
//! no-claim: dataset admission and retrieval are not solver verification and
//! cannot manufacture physical validation.

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_evidence::{
    ColorRank, Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain,
};
use fs_qty::{Dims, QtyAny};
use std::collections::BTreeSet;
use std::fmt;
use std::path::Component;
use std::sync::LazyLock;

use crate::thermal_level_a::{
    THERMAL_LEVEL_A_MANIFEST, ThermalLevelAAcceptance, ThermalLevelACase, thermal_level_a_cases,
};

/// Current canonical corpus wire and identity schema.
pub const CORPUS_SCHEMA_VERSION: u32 = 1;
/// Maximum admitted datasets in one caller-built registry.
pub const MAX_CORPUS_DATASETS: usize = 4_096;
/// Maximum sensors on one dataset.
pub const MAX_DATASET_SENSORS: usize = 4_096;
/// Maximum elements in any other dataset collection.
pub const MAX_DATASET_ITEMS: usize = 4_096;
/// Maximum UTF-8 bytes in any schema string.
pub const MAX_CORPUS_TEXT_BYTES: usize = 4_096;
/// Maximum encoded bytes for one dataset.
pub const MAX_DATASET_CANONICAL_BYTES: usize = 16 * 1024 * 1024;

const DATASET_DOMAIN: &str = "org.frankensim.fs-vvreg.corpus-dataset.v1";
const REGISTRY_DOMAIN: &str = "org.frankensim.fs-vvreg.corpus-registry.v1";
const MAGIC: &[u8; 8] = b"FSVVCRP\0";
const RAW_CHT_FIXTURE: &[u8] =
    include_bytes!("../../../data/vv-corpus/fs-benchmark-cht-query-v1/raw-sensors.csv");
const MARTIN_MOYCE_FIXTURE: &[u8] =
    include_bytes!("../../../data/reference/martin-moyce-1952.jsonl");

/// Top-level required field named by a typed missing-field refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DatasetField {
    /// Stable dataset identifier.
    Id,
    /// Human-readable dataset title.
    Title,
    /// Earliest retained payload and its retention state.
    RawPayload,
    /// Sensor or synthetic-channel roster.
    Sensors,
    /// Geometry authority or its declared absence.
    Geometry,
    /// Acquisition environment or its declared absence.
    Environment,
    /// Training, calibration, or validation partition.
    Partition,
    /// Complete or explicitly unreplayable preprocessing history.
    Preprocessing,
    /// Final retained transform output.
    FinalArtifact,
    /// Bounded context in which the row may be queried.
    ContextOfUse,
    /// License authority or its declared absence.
    License,
    /// Acquisition/source provenance.
    Provenance,
    /// Payload and calibration retention policy.
    Retention,
    /// Metric-specific acceptance rules.
    AcceptanceEnvelopes,
    /// Portfolio evidence level A through E.
    EvidenceLevel,
}

impl DatasetField {
    /// Stable dotted schema name used in diagnostics and audit output.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Id => "id",
            Self::Title => "title",
            Self::RawPayload => "raw_payload",
            Self::Sensors => "sensors",
            Self::Geometry => "geometry",
            Self::Environment => "environment",
            Self::Partition => "partition",
            Self::Preprocessing => "preprocessing",
            Self::FinalArtifact => "final_artifact",
            Self::ContextOfUse => "context_of_use",
            Self::License => "license",
            Self::Provenance => "provenance",
            Self::Retention => "retention",
            Self::AcceptanceEnvelopes => "acceptance_envelopes",
            Self::EvidenceLevel => "evidence_level",
        }
    }
}

/// Corpus admission or canonical-codec refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum CorpusError {
    /// A mandatory top-level field was omitted.
    MissingField {
        /// Omitted field.
        field: DatasetField,
    },
    /// A present field failed semantic validation.
    InvalidField {
        /// Invalid field.
        field: DatasetField,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// A bounded collection or payload exceeded its cap.
    ResourceLimit {
        /// Capped resource name.
        resource: &'static str,
        /// Maximum accepted value.
        limit: usize,
        /// Observed value.
        observed: usize,
    },
    /// Two rows declared the same dataset id.
    DuplicateDatasetId {
        /// Conflicting id.
        id: String,
    },
    /// Two sensors declared the same id.
    DuplicateSensorId {
        /// Conflicting id.
        id: String,
    },
    /// A named subcollection contains a duplicate key.
    DuplicateName {
        /// Subcollection name.
        collection: &'static str,
        /// Conflicting key.
        name: String,
    },
    /// Preprocessing hashes or ordinals do not form the declared lineage.
    BrokenLineage {
        /// Zero-based failing transform position.
        step: usize,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// Canonical bytes have the wrong magic prefix.
    BadMagic,
    /// Canonical bytes use an unsupported schema version.
    UnsupportedSchema {
        /// Version read from the payload.
        observed: u32,
    },
    /// Canonical bytes end before a declared field is complete.
    Truncated,
    /// Canonical bytes contain unconsumed suffix data.
    TrailingBytes {
        /// Number of unconsumed bytes.
        count: usize,
    },
    /// A length-framed string is not valid UTF-8.
    InvalidUtf8,
    /// A wire discriminant is unknown.
    InvalidTag {
        /// Discriminated value kind.
        kind: &'static str,
        /// Unknown tag byte.
        tag: u8,
    },
    /// Decoded semantics re-encode to different bytes.
    NonCanonicalEncoding,
}

impl fmt::Display for CorpusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingField { field } => {
                write!(f, "dataset is missing mandatory field `{}`", field.name())
            }
            Self::InvalidField { field, reason } => {
                write!(f, "dataset field `{}` is invalid: {reason}", field.name())
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "corpus resource `{resource}` exceeds limit {limit} (observed {observed})"
            ),
            Self::DuplicateDatasetId { id } => write!(f, "duplicate dataset id `{id}`"),
            Self::DuplicateSensorId { id } => write!(f, "duplicate sensor id `{id}`"),
            Self::DuplicateName { collection, name } => {
                write!(f, "duplicate {collection} name `{name}`")
            }
            Self::BrokenLineage { step, reason } => {
                write!(f, "preprocessing lineage breaks at step {step}: {reason}")
            }
            Self::BadMagic => f.write_str("dataset bytes do not start with FSVVCRP magic"),
            Self::UnsupportedSchema { observed } => {
                write!(f, "unsupported corpus schema {observed}")
            }
            Self::Truncated => f.write_str("dataset bytes are truncated"),
            Self::TrailingBytes { count } => write!(f, "dataset bytes have {count} trailing bytes"),
            Self::InvalidUtf8 => f.write_str("dataset bytes contain invalid UTF-8"),
            Self::InvalidTag { kind, tag } => write!(f, "invalid {kind} tag {tag}"),
            Self::NonCanonicalEncoding => {
                f.write_str("dataset bytes are semantically valid but not canonical")
            }
        }
    }
}

impl std::error::Error for CorpusError {}

/// A content-addressed payload retained by the corpus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusArtifact {
    /// Content digest of the retained bytes.
    pub digest: ContentHash,
    /// Exact retained byte length.
    pub byte_len: u64,
    /// Declared media type.
    pub media_type: String,
    /// Normalized repository-relative locator.
    pub locator: String,
}

/// Presence of a required semantic field. `Unavailable` is not silently
/// defaulted: it carries a reason, remains identity-bearing, and demotes use.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Availability<T> {
    /// The value is present and may be validated.
    Available(T),
    /// The value is absent for the stated identity-bearing reason.
    Unavailable {
        /// Why the value cannot be supplied.
        reason: String,
    },
}

/// Whether the retained source artifact is original raw data or the earliest
/// surviving derived representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PayloadRetention {
    /// The retained artifact is the original raw acquisition.
    OriginalRaw(CorpusArtifact),
    /// Only a derived artifact survives.
    DerivedOnly {
        /// Earliest surviving artifact.
        retained: CorpusArtifact,
        /// Why original raw data are unavailable.
        reason: String,
    },
}

impl PayloadRetention {
    /// Earliest retained artifact, regardless of originality.
    #[must_use]
    pub const fn artifact(&self) -> &CorpusArtifact {
        match self {
            Self::OriginalRaw(artifact) => artifact,
            Self::DerivedOnly { retained, .. } => retained,
        }
    }
}

/// Exact calibration identity and validity dates for one sensor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CalibrationRecord {
    /// Calibration certificate identifier.
    pub certificate_id: String,
    /// Declared certificate content hash.
    pub certificate_hash: ContentHash,
    /// Certificate issue date.
    pub issued_on: String,
    /// Optional last-valid date.
    pub valid_through: Option<String>,
}

/// Sensor placement and componentwise placement uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorPlacement {
    /// Coordinate frame identifier.
    pub frame: String,
    /// Three-dimensional sensor coordinates.
    pub coordinates: [QtyAny; 3],
    /// Componentwise placement uncertainty.
    pub uncertainty: [QtyAny; 3],
}

/// Measurement uncertainty is mandatory as a field. `Unstated` is an
/// explicit admitted no-claim that caps every use at `Estimated`.
#[derive(Debug, Clone, PartialEq)]
pub enum MeasurementUncertainty {
    /// Symmetric absolute half-width.
    Bounded {
        /// Non-negative half-width in the measured quantity dimensions.
        half_width: QtyAny,
    },
    /// One diagonal covariance entry.
    CovarianceDiagonal {
        /// Non-negative variance in squared quantity dimensions.
        variance: QtyAny,
    },
    /// The retained source supplies no quantitative uncertainty.
    Unstated,
}

/// One sensor/channel represented in the retained raw payload.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorRecord {
    /// Stable channel identifier.
    pub id: String,
    /// Instrument identity or explicit absence.
    pub instrument_id: Availability<String>,
    /// Column/channel name in the retained payload.
    pub raw_channel: String,
    /// Dimensions of measured values.
    pub quantity_dims: Dims,
    /// Calibration authority or explicit absence.
    pub calibration: Availability<CalibrationRecord>,
    /// Placement authority or explicit absence.
    pub placement: Availability<SensorPlacement>,
    /// Quantitative uncertainty declaration.
    pub uncertainty: MeasurementUncertainty,
}

/// Nominal and optional as-built geometry artifacts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometryRecord {
    /// Nominal geometry artifact.
    pub nominal: CorpusArtifact,
    /// Optional as-built geometry artifact.
    pub as_built: Option<CorpusArtifact>,
    /// Geometry coordinate frame.
    pub frame: String,
}

/// A measured environmental condition at acquisition time.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvironmentCondition {
    /// Stable condition name.
    pub name: String,
    /// Recorded condition value.
    pub value: QtyAny,
    /// Non-negative condition uncertainty.
    pub uncertainty: QtyAny,
}

/// Declared training/calibration/validation partition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DatasetPartition {
    /// Data available to fitting/training procedures.
    Training,
    /// Data reserved for calibration.
    Calibration,
    /// Data reserved for validation.
    Validation,
}

impl DatasetPartition {
    /// Stable lowercase partition name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Training => "training",
            Self::Calibration => "calibration",
            Self::Validation => "validation",
        }
    }
}

/// Every preprocessing transform forms one exact input/output hash edge.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreprocessingStep {
    /// Contiguous zero-based transform ordinal.
    pub ordinal: u32,
    /// Transform operation identifier.
    pub operation: String,
    /// Transform implementation/version identifier.
    pub version: String,
    /// Exact input artifact hash.
    pub input: ContentHash,
    /// Exact output artifact hash.
    pub output: ContentHash,
}

/// Complete transform chain or an explicit unreplayable historical lineage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreprocessingLineage {
    /// Fully replayable, gap-free transform chain.
    Complete(Vec<PreprocessingStep>),
    /// Historical transform chain cannot be replayed.
    Unreplayable {
        /// Earliest retained input hash.
        retained_input: ContentHash,
        /// Retained output hash.
        retained_output: ContentHash,
        /// Why the transform chain cannot be replayed.
        reason: String,
    },
}

/// Inclusive typed validity range.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextRange {
    /// Stable context-axis name.
    pub name: String,
    /// Inclusive lower bound.
    pub lo: QtyAny,
    /// Inclusive upper bound.
    pub hi: QtyAny,
}

/// One caller-supplied context coordinate for a query.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextValue {
    /// Context-axis name.
    pub name: String,
    /// Requested coordinate.
    pub value: QtyAny,
}

/// License and redistribution policy for the retained artifacts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RedistributionPolicy {
    /// Retained bytes may be redistributed.
    Allowed,
    /// Only metadata may be redistributed.
    MetadataOnly,
    /// Redistribution is prohibited.
    Prohibited,
}

/// Corpus-level license declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusLicense {
    /// License or policy identifier.
    pub identifier: String,
    /// Human-readable terms or policy summary.
    pub terms: String,
    /// Redistribution decision.
    pub redistribution: RedistributionPolicy,
}

/// Acquisition provenance. Instrument identity lives on each sensor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AcquisitionProvenance {
    /// People or source authors responsible for acquisition.
    pub measured_by: String,
    /// Responsible organization or publication.
    pub organization: String,
    /// Acquisition date or explicit absence.
    pub measured_on: Availability<String>,
    /// Exact retained source-record locator.
    pub source_record: String,
}

/// How long artifacts remain available and what must be retained together.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RetentionClass {
    /// Retain without a scheduled expiry.
    Permanent,
    /// Retain for the declared number of years.
    Years(u16),
}

/// Retention policy; raw data and calibration evidence are inseparable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetentionPolicy {
    /// Retention duration class.
    pub class: RetentionClass,
    /// Whether retained payloads must be preserved.
    pub preserve_raw: bool,
    /// Whether calibration evidence must be preserved.
    pub preserve_calibration: bool,
    /// Stable governing policy identifier.
    pub policy_id: String,
}

/// Arithmetic acceptance rule for one validation metric.
#[derive(Debug, Clone, PartialEq)]
pub enum CorpusEnvelope {
    /// Absolute-plus-relative scalar tolerance.
    Tolerance {
        /// Non-negative absolute tolerance.
        atol: f64,
        /// Non-negative relative tolerance.
        rtol: f64,
    },
    /// Inclusive scalar acceptance interval.
    Interval {
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
    /// No defensible scalar envelope is pinned.
    Unpinned {
        /// Scientific or historical basis for the no-claim.
        basis: String,
    },
}

/// Acceptance rule plus the exact regime in which it applies.
#[derive(Debug, Clone, PartialEq)]
pub struct AcceptanceRecord {
    /// Stable validation metric id.
    pub metric: String,
    /// Metric dimensions.
    pub dims: Dims,
    /// Pinned rule or explicit unpinned state.
    pub envelope: CorpusEnvelope,
    /// Context subdomain where the rule applies.
    pub regime: Vec<ContextRange>,
}

/// Portfolio evidence doctrine A through E.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidenceLevel {
    /// Level A analytic or manufactured evidence.
    Analytic,
    /// Level B independent cross-code evidence.
    CrossCode,
    /// Level C published experiment.
    PublishedExperiment,
    /// Level D preregistered blind evidence.
    Blind,
    /// Level E field evidence.
    Field,
}

impl EvidenceLevel {
    /// Stable A-E portfolio code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Analytic => "A",
            Self::CrossCode => "B",
            Self::PublishedExperiment => "C",
            Self::Blind => "D",
            Self::Field => "E",
        }
    }
}

/// Untrusted, optional-field input. Admission names the first missing
/// mandatory field in schema order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DatasetDraft {
    /// Candidate stable id.
    pub id: Option<String>,
    /// Candidate title.
    pub title: Option<String>,
    /// Candidate earliest-retained payload.
    pub raw_payload: Option<PayloadRetention>,
    /// Candidate sensor/channel roster.
    pub sensors: Option<Vec<SensorRecord>>,
    /// Candidate geometry authority.
    pub geometry: Option<Availability<GeometryRecord>>,
    /// Candidate acquisition environment.
    pub environment: Option<Availability<Vec<EnvironmentCondition>>>,
    /// Candidate dataset partition.
    pub partition: Option<DatasetPartition>,
    /// Candidate preprocessing lineage.
    pub preprocessing: Option<PreprocessingLineage>,
    /// Candidate final artifact hash.
    pub final_artifact: Option<ContentHash>,
    /// Candidate context-of-use box.
    pub context_of_use: Option<Vec<ContextRange>>,
    /// Candidate license authority.
    pub license: Option<Availability<CorpusLicense>>,
    /// Candidate acquisition/source provenance.
    pub provenance: Option<AcquisitionProvenance>,
    /// Candidate retention policy.
    pub retention: Option<RetentionPolicy>,
    /// Candidate acceptance records.
    pub acceptance_envelopes: Option<Vec<AcceptanceRecord>>,
    /// Candidate evidence portfolio level.
    pub evidence_level: Option<EvidenceLevel>,
}

/// Admitted, immutable dataset. Fields are private so admission cannot be
/// bypassed by a struct literal.
#[derive(Debug, Clone, PartialEq)]
pub struct CorpusDataset {
    id: String,
    title: String,
    raw_payload: PayloadRetention,
    sensors: Vec<SensorRecord>,
    geometry: Availability<GeometryRecord>,
    environment: Availability<Vec<EnvironmentCondition>>,
    partition: DatasetPartition,
    preprocessing: PreprocessingLineage,
    final_artifact: ContentHash,
    context_of_use: Vec<ContextRange>,
    license: Availability<CorpusLicense>,
    provenance: AcquisitionProvenance,
    retention: RetentionPolicy,
    acceptance_envelopes: Vec<AcceptanceRecord>,
    evidence_level: EvidenceLevel,
}

impl CorpusDataset {
    /// Stable dataset id.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Human-readable dataset title.
    #[must_use]
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Earliest-retained payload declaration.
    #[must_use]
    pub const fn raw_payload(&self) -> &PayloadRetention {
        &self.raw_payload
    }

    /// Sorted sensor/channel roster.
    #[must_use]
    pub fn sensors(&self) -> &[SensorRecord] {
        &self.sensors
    }

    /// Geometry authority or explicit absence.
    #[must_use]
    pub const fn geometry(&self) -> &Availability<GeometryRecord> {
        &self.geometry
    }

    /// Acquisition environment or explicit absence.
    #[must_use]
    pub const fn environment(&self) -> &Availability<Vec<EnvironmentCondition>> {
        &self.environment
    }

    /// Declared dataset partition.
    #[must_use]
    pub const fn partition(&self) -> DatasetPartition {
        self.partition
    }

    /// Complete or unreplayable preprocessing lineage.
    #[must_use]
    pub const fn preprocessing(&self) -> &PreprocessingLineage {
        &self.preprocessing
    }

    /// Final retained artifact hash.
    #[must_use]
    pub const fn final_artifact(&self) -> ContentHash {
        self.final_artifact
    }

    /// Sorted context-of-use axes.
    #[must_use]
    pub fn context_of_use(&self) -> &[ContextRange] {
        &self.context_of_use
    }

    /// License authority or explicit absence.
    #[must_use]
    pub const fn license(&self) -> &Availability<CorpusLicense> {
        &self.license
    }

    /// Acquisition/source provenance.
    #[must_use]
    pub const fn provenance(&self) -> &AcquisitionProvenance {
        &self.provenance
    }

    /// Governing retention policy.
    #[must_use]
    pub const fn retention(&self) -> &RetentionPolicy {
        &self.retention
    }

    /// Sorted metric acceptance records.
    #[must_use]
    pub fn acceptance_envelopes(&self) -> &[AcceptanceRecord] {
        &self.acceptance_envelopes
    }

    /// Portfolio evidence level.
    #[must_use]
    pub const fn evidence_level(&self) -> EvidenceLevel {
        self.evidence_level
    }

    /// Maximum physical claim rank this dataset may support. A/B evidence is
    /// not a physical experiment, and any explicit gap in raw retention,
    /// metrology, geometry, environment, lineage, licensing, acquisition date,
    /// or a pinned acceptance envelope demotes C/D/E to `Estimated`.
    #[must_use]
    pub fn physical_claim_cap(&self) -> ColorRank {
        if matches!(self.raw_payload, PayloadRetention::DerivedOnly { .. })
            || matches!(self.geometry, Availability::Unavailable { .. })
            || matches!(self.environment, Availability::Unavailable { .. })
            || matches!(
                self.preprocessing,
                PreprocessingLineage::Unreplayable { .. }
            )
            || matches!(self.license, Availability::Unavailable { .. })
            || matches!(
                self.provenance.measured_on,
                Availability::Unavailable { .. }
            )
            || self.sensors.iter().any(|sensor| {
                matches!(sensor.instrument_id, Availability::Unavailable { .. })
                    || matches!(sensor.calibration, Availability::Unavailable { .. })
                    || matches!(sensor.placement, Availability::Unavailable { .. })
                    || sensor.uncertainty == MeasurementUncertainty::Unstated
            })
            || self
                .acceptance_envelopes
                .iter()
                .any(|record| matches!(record.envelope, CorpusEnvelope::Unpinned { .. }))
        {
            return ColorRank::Estimated;
        }
        match self.evidence_level {
            EvidenceLevel::PublishedExperiment | EvidenceLevel::Blind | EvidenceLevel::Field => {
                ColorRank::Validated
            }
            EvidenceLevel::Analytic | EvidenceLevel::CrossCode => ColorRank::Estimated,
        }
    }

    /// Canonical dataset identity.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        hash_domain(DATASET_DOMAIN, &self.encode())
    }
}

/// Admit an untrusted dataset draft.
pub fn admit_dataset(draft: DatasetDraft) -> Result<CorpusDataset, CorpusError> {
    let mut dataset = CorpusDataset {
        id: require(draft.id, DatasetField::Id)?,
        title: require(draft.title, DatasetField::Title)?,
        raw_payload: require(draft.raw_payload, DatasetField::RawPayload)?,
        sensors: require(draft.sensors, DatasetField::Sensors)?,
        geometry: require(draft.geometry, DatasetField::Geometry)?,
        environment: require(draft.environment, DatasetField::Environment)?,
        partition: require(draft.partition, DatasetField::Partition)?,
        preprocessing: require(draft.preprocessing, DatasetField::Preprocessing)?,
        final_artifact: require(draft.final_artifact, DatasetField::FinalArtifact)?,
        context_of_use: require(draft.context_of_use, DatasetField::ContextOfUse)?,
        license: require(draft.license, DatasetField::License)?,
        provenance: require(draft.provenance, DatasetField::Provenance)?,
        retention: require(draft.retention, DatasetField::Retention)?,
        acceptance_envelopes: require(
            draft.acceptance_envelopes,
            DatasetField::AcceptanceEnvelopes,
        )?,
        evidence_level: require(draft.evidence_level, DatasetField::EvidenceLevel)?,
    };
    dataset.sensors.sort_by(|a, b| a.id.cmp(&b.id));
    if let Availability::Available(environment) = &mut dataset.environment {
        environment.sort_by(|a, b| a.name.cmp(&b.name));
    }
    if let PreprocessingLineage::Complete(preprocessing) = &mut dataset.preprocessing {
        preprocessing.sort_by_key(|step| step.ordinal);
    }
    dataset.context_of_use.sort_by(|a, b| a.name.cmp(&b.name));
    dataset
        .acceptance_envelopes
        .sort_by(|a, b| a.metric.cmp(&b.metric));
    for acceptance in &mut dataset.acceptance_envelopes {
        acceptance.regime.sort_by(|a, b| a.name.cmp(&b.name));
    }
    validate_dataset(&dataset)?;
    Ok(dataset)
}

fn require<T>(value: Option<T>, field: DatasetField) -> Result<T, CorpusError> {
    value.ok_or(CorpusError::MissingField { field })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CorpusAuthority {
    Seeded,
    Unauthoritative,
}

/// Sorted dataset registry. Caller-built registries can be audited and
/// serialized, but only [`corpus`] has query authority.
#[derive(Debug, Clone, PartialEq)]
pub struct CorpusRegistry {
    datasets: Vec<CorpusDataset>,
    authority: CorpusAuthority,
}

/// Fail-closed dataset query refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum CorpusQueryRefusal {
    /// Caller-built registries cannot return evidence wrappers.
    UnauthoritativeRegistry,
    /// Requested id violates the bounded slug grammar.
    InvalidDatasetId,
    /// No seeded row has the requested id.
    UnknownDataset {
        /// Requested id.
        id: String,
    },
    /// The stored row no longer passes schema validation.
    InvalidDataset(CorpusError),
    /// Requested partition differs from the declared partition.
    PartitionMismatch {
        /// Dataset's declared partition.
        declared: DatasetPartition,
        /// Caller's requested partition.
        requested: DatasetPartition,
    },
    /// A required context axis was omitted.
    MissingContext {
        /// Omitted axis.
        name: String,
    },
    /// Caller supplied an axis not declared by the dataset.
    UnknownContext {
        /// Unknown axis.
        name: String,
    },
    /// Caller supplied the same context axis more than once.
    DuplicateContext {
        /// Duplicated axis.
        name: String,
    },
    /// Context coordinate dimensions differ from the declared axis.
    ContextDimensionMismatch {
        /// Mismatched axis.
        name: String,
        /// Declared dimensions.
        expected: Dims,
        /// Supplied dimensions.
        observed: Dims,
    },
    /// Context coordinate is non-finite or outside the inclusive range.
    OutOfContext {
        /// Out-of-range axis.
        name: String,
        /// Supplied scalar value.
        value: f64,
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
}

impl fmt::Display for CorpusQueryRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnauthoritativeRegistry => f.write_str(
                "caller-built corpus has no evidence-query authority; use the seeded corpus",
            ),
            Self::InvalidDatasetId => f.write_str("dataset id is not a bounded lowercase slug"),
            Self::UnknownDataset { id } => write!(f, "unknown dataset `{id}`"),
            Self::InvalidDataset(error) => write!(f, "dataset failed revalidation: {error}"),
            Self::PartitionMismatch {
                declared,
                requested,
            } => write!(
                f,
                "partition mismatch: dataset is {}, request is {}",
                declared.name(),
                requested.name()
            ),
            Self::MissingContext { name } => write!(f, "missing context coordinate `{name}`"),
            Self::UnknownContext { name } => write!(f, "unknown context coordinate `{name}`"),
            Self::DuplicateContext { name } => write!(f, "duplicate context coordinate `{name}`"),
            Self::ContextDimensionMismatch { name, .. } => {
                write!(f, "context coordinate `{name}` has incompatible dimensions")
            }
            Self::OutOfContext {
                name,
                value,
                lo,
                hi,
            } => write!(
                f,
                "context coordinate `{name}` = {value} is outside inclusive [{lo}, {hi}]"
            ),
        }
    }
}

impl std::error::Error for CorpusQueryRefusal {}

impl CorpusRegistry {
    /// Build a validated but unauthoritative registry.
    pub fn build(drafts: Vec<DatasetDraft>) -> Result<Self, CorpusError> {
        if drafts.len() > MAX_CORPUS_DATASETS {
            return Err(CorpusError::ResourceLimit {
                resource: "datasets",
                limit: MAX_CORPUS_DATASETS,
                observed: drafts.len(),
            });
        }
        let mut datasets = drafts
            .into_iter()
            .map(admit_dataset)
            .collect::<Result<Vec<_>, _>>()?;
        datasets.sort_by(|a, b| a.id.cmp(&b.id));
        for pair in datasets.windows(2) {
            if pair[0].id == pair[1].id {
                return Err(CorpusError::DuplicateDatasetId {
                    id: pair[0].id.clone(),
                });
            }
        }
        Ok(Self {
            datasets,
            authority: CorpusAuthority::Unauthoritative,
        })
    }

    /// Sorted admitted dataset rows.
    #[must_use]
    pub fn datasets(&self) -> &[CorpusDataset] {
        &self.datasets
    }

    /// Find an admitted row by exact id.
    #[must_use]
    pub fn dataset(&self, id: &str) -> Option<&CorpusDataset> {
        self.datasets
            .binary_search_by_key(&id, |dataset| dataset.id.as_str())
            .ok()
            .map(|index| &self.datasets[index])
    }

    /// Query one dataset under an exact partition and complete context.
    /// The returned [`Evidence`] is intentionally non-certifying.
    pub fn query<'a>(
        &'a self,
        id: &str,
        partition: DatasetPartition,
        context: &[ContextValue],
    ) -> Result<Evidence<&'a CorpusDataset>, CorpusQueryRefusal> {
        if self.authority != CorpusAuthority::Seeded {
            return Err(CorpusQueryRefusal::UnauthoritativeRegistry);
        }
        if !valid_slug(id) {
            return Err(CorpusQueryRefusal::InvalidDatasetId);
        }
        let dataset = self
            .dataset(id)
            .ok_or_else(|| CorpusQueryRefusal::UnknownDataset { id: id.to_string() })?;
        validate_dataset(dataset).map_err(CorpusQueryRefusal::InvalidDataset)?;
        if partition != dataset.partition {
            return Err(CorpusQueryRefusal::PartitionMismatch {
                declared: dataset.partition,
                requested: partition,
            });
        }
        validate_query_context(dataset, context)?;
        let validity = dataset
            .context_of_use
            .iter()
            .fold(ValidityDomain::unconstrained(), |domain, range| {
                domain.with(range.name.clone(), range.lo.value, range.hi.value)
            });
        Ok(Evidence {
            value: dataset,
            qoi: 0.0,
            numerical: NumericalCertificate::no_claim(),
            statistical: StatisticalCertificate::HalfWidth {
                half_width: f64::INFINITY,
                confidence: 0.0,
            },
            model: ModelEvidence {
                cards: vec![format!("vv-corpus:{}", dataset.id)],
                assumptions: vec![
                    "corpus registration does not quantify model-form discrepancy".to_string(),
                ],
                validity,
                discrepancy_rel: f64::INFINITY,
                in_domain: true,
            },
            sensitivity: SensitivitySummary::default(),
            provenance: ProvenanceHash::of_bytes(dataset.digest().as_bytes()),
            adjoint_ref: None,
        })
    }

    /// Deterministic registry identity over sorted canonical dataset bytes.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(&CORPUS_SCHEMA_VERSION.to_le_bytes());
        push_u32(&mut payload, self.datasets.len());
        for dataset in &self.datasets {
            let bytes = dataset.encode();
            push_u32(&mut payload, bytes.len());
            payload.extend_from_slice(&bytes);
        }
        hash_domain(REGISTRY_DOMAIN, &payload)
    }
}

static CORPUS: LazyLock<CorpusRegistry> = LazyLock::new(|| {
    let mut datasets = vec![seed_cht_dataset(), seed_martin_moyce_dataset()];
    datasets.extend(
        thermal_level_a_cases()
            .iter()
            .map(seed_thermal_level_a_dataset),
    );
    datasets.sort_by(|a, b| a.id.cmp(&b.id));
    CorpusRegistry {
        datasets,
        authority: CorpusAuthority::Seeded,
    }
});

/// Seeded workspace corpus: reference-only Level-A thermal rows, one explicitly
/// synthetic Level-B CHT fixture, and one Level-C Martin-Moyce retained curve.
/// Every query remains non-certifying; authority gaps remain explicit.
#[must_use]
pub fn corpus() -> &'static CorpusRegistry {
    &CORPUS
}

fn seed_thermal_level_a_dataset(case: &ThermalLevelACase) -> CorpusDataset {
    let retained = hash_bytes(THERMAL_LEVEL_A_MANIFEST);
    let context_of_use = case
        .context
        .iter()
        .map(|axis| ContextRange {
            name: axis.name.to_string(),
            lo: QtyAny::new(axis.lo, axis.dims),
            hi: QtyAny::new(axis.hi, axis.dims),
        })
        .collect::<Vec<_>>();
    let envelope = match case.acceptance {
        ThermalLevelAAcceptance::Tolerance { atol, rtol } => {
            CorpusEnvelope::Tolerance { atol, rtol }
        }
        ThermalLevelAAcceptance::OrderGate {
            theoretical,
            tolerance,
        } => CorpusEnvelope::Interval {
            lo: theoretical - tolerance,
            hi: theoretical + tolerance,
        },
    };
    CorpusDataset {
        id: case.id.to_string(),
        title: case.title.to_string(),
        raw_payload: PayloadRetention::DerivedOnly {
            retained: CorpusArtifact {
                digest: retained,
                byte_len: THERMAL_LEVEL_A_MANIFEST.len() as u64,
                media_type: "text/tab-separated-values".to_string(),
                locator: "data/vv-corpus/thermal-level-a/thermal-level-a-v1.tsv".to_string(),
            },
            reason: "this is an authored analytic reference/target manifest; no sensor acquisition occurred"
                .to_string(),
        },
        sensors: vec![SensorRecord {
            id: "reference-answer".to_string(),
            instrument_id: unavailable(
                "an analytic formula or theoretical order target has no physical instrument",
            ),
            raw_channel: case.metric.to_string(),
            quantity_dims: case.metric_dims,
            calibration: unavailable(
                "an analytic reference has derivation checks, not a calibration certificate",
            ),
            placement: unavailable(
                "the formula manifest declares idealized parameters, not a sensor placement",
            ),
            uncertainty: MeasurementUncertainty::Unstated,
        }],
        geometry: unavailable(
            "idealized geometry is declared in the retained formula row, but no standalone geometry or as-built artifact is bound",
        ),
        environment: unavailable(
            "analytic references and G1 targets are not physical acquisitions with an environment record",
        ),
        partition: DatasetPartition::Validation,
        preprocessing: PreprocessingLineage::Complete(vec![PreprocessingStep {
            ordinal: 0,
            operation: "analytic-manifest-identity-import".to_string(),
            version: "1".to_string(),
            input: retained,
            output: retained,
        }]),
        final_artifact: retained,
        context_of_use: context_of_use.clone(),
        license: Availability::Available(CorpusLicense {
            identifier: "MIT".to_string(),
            terms: "Repository-authored formula and target manifest; redistribution allowed"
                .to_string(),
            redistribution: RedistributionPolicy::Allowed,
        }),
        provenance: AcquisitionProvenance {
            measured_by: "FrankenSim Level-A thermal reference derivation".to_string(),
            organization: "FrankenSim".to_string(),
            measured_on: unavailable(
                "the row is a versioned analytic definition, not a dated physical acquisition",
            ),
            source_record: format!(
                "data/vv-corpus/thermal-level-a/thermal-level-a-v1.tsv:{}",
                case.id
            ),
        },
        retention: RetentionPolicy {
            class: RetentionClass::Permanent,
            preserve_raw: true,
            preserve_calibration: true,
            policy_id: "frankensim-vv-corpus-permanent-v1".to_string(),
        },
        acceptance_envelopes: vec![AcceptanceRecord {
            metric: case.metric.to_string(),
            dims: case.metric_dims,
            envelope,
            regime: context_of_use,
        }],
        evidence_level: EvidenceLevel::Analytic,
    }
}

fn seed_cht_dataset() -> CorpusDataset {
    let raw = hash_bytes(RAW_CHT_FIXTURE);
    let temperature = Dims([0, 0, 0, 1, 0, 0]);
    CorpusDataset {
        id: "fs-benchmark-cht-query-v1".to_string(),
        title: "Synthetic tabulation of the fs-benchmark CHT query fixture".to_string(),
        raw_payload: PayloadRetention::DerivedOnly {
            retained: CorpusArtifact {
                digest: raw,
                byte_len: RAW_CHT_FIXTURE.len() as u64,
                media_type: "text/csv".to_string(),
                locator: "data/vv-corpus/fs-benchmark-cht-query-v1/raw-sensors.csv".to_string(),
            },
            reason: "the retained CSV is an authored tabulation of a hard-coded synthetic query, not original sensor data"
                .to_string(),
        },
        sensors: vec![SensorRecord {
            id: "cht-q3".to_string(),
            instrument_id: unavailable(
                "the synthetic query was not acquired by a physical instrument",
            ),
            raw_channel: "hotspot_thermal_margin".to_string(),
            quantity_dims: temperature,
            calibration: unavailable(
                "the query tolerance is an acceptance rule, not a calibration certificate",
            ),
            placement: unavailable(
                "the synthetic query has no physical sensor placement or placement uncertainty",
            ),
            uncertainty: MeasurementUncertainty::Unstated,
        }],
        geometry: unavailable(
            "the retained nominal-geometry file is an explicit disclaimer, not a geometry artifact",
        ),
        environment: unavailable(
            "the hard-coded synthetic query records no physical acquisition environment",
        ),
        partition: DatasetPartition::Validation,
        preprocessing: PreprocessingLineage::Unreplayable {
            retained_input: raw,
            retained_output: raw,
            reason: "the one-row tabulation was authored manually; no replayable exporter and versioned transform parameters are retained"
                .to_string(),
        },
        final_artifact: raw,
        context_of_use: vec![ContextRange {
            name: "reference_cost_work_units".to_string(),
            lo: QtyAny::dimensionless(250.0),
            hi: QtyAny::dimensionless(250.0),
        }],
        license: Availability::Available(CorpusLicense {
            identifier: "MIT".to_string(),
            terms: "Repository-authored synthetic fixture; redistribution allowed".to_string(),
            redistribution: RedistributionPolicy::Allowed,
        }),
        provenance: AcquisitionProvenance {
            measured_by: "fs-benchmark deterministic query fixture".to_string(),
            organization: "FrankenSim".to_string(),
            measured_on: unavailable(
                "the synthetic query has a source revision, not a physical acquisition date",
            ),
            source_record: "crates/fs-benchmark/src/lib.rs:cht-q3".to_string(),
        },
        retention: RetentionPolicy {
            class: RetentionClass::Permanent,
            preserve_raw: true,
            preserve_calibration: true,
            policy_id: "frankensim-vv-corpus-permanent-v1".to_string(),
        },
        acceptance_envelopes: vec![AcceptanceRecord {
            metric: "hotspot_thermal_margin".to_string(),
            dims: temperature,
            envelope: CorpusEnvelope::Tolerance {
                atol: 1.0,
                rtol: 0.0,
            },
            regime: vec![ContextRange {
                name: "reference_cost_work_units".to_string(),
                lo: QtyAny::dimensionless(250.0),
                hi: QtyAny::dimensionless(250.0),
            }],
        }],
        evidence_level: EvidenceLevel::CrossCode,
    }
}

fn seed_martin_moyce_dataset() -> CorpusDataset {
    let retained = hash_bytes(MARTIN_MOYCE_FIXTURE);
    let dimensionless = Dims::NONE;
    CorpusDataset {
        id: "martin-moyce-1952-square-column".to_string(),
        title: "Digitized Martin-Moyce square-column collapse curve".to_string(),
        raw_payload: PayloadRetention::DerivedOnly {
            retained: CorpusArtifact {
                digest: retained,
                byte_len: MARTIN_MOYCE_FIXTURE.len() as u64,
                media_type: "application/x-ndjson".to_string(),
                locator: "data/reference/martin-moyce-1952.jsonl".to_string(),
            },
            reason: "only digitized figure coordinates survive; no cine frames, raw timing records, or calibration frames are retained"
                .to_string(),
        },
        sensors: vec![SensorRecord {
            id: "surge-front-coordinate".to_string(),
            instrument_id: unavailable(
                "the retained artifact does not identify the original imaging/timing instruments or the digitizer",
            ),
            raw_channel: "t_star,z".to_string(),
            quantity_dims: dimensionless,
            calibration: unavailable(
                "no calibration certificate for the imaging, timing, length-scale, or digitization chain is retained",
            ),
            placement: unavailable(
                "camera geometry, measurement station, and placement uncertainty are unrecorded",
            ),
            uncertainty: MeasurementUncertainty::Unstated,
        }],
        geometry: unavailable(
            "the retained curve comment describes nominal nondimensional geometry but carries no geometry artifact or as-built measurement",
        ),
        environment: unavailable(
            "ambient/fluid temperature, surface condition, and other acquisition environment fields are unrecorded",
        ),
        partition: DatasetPartition::Validation,
        preprocessing: PreprocessingLineage::Unreplayable {
            retained_input: retained,
            retained_output: retained,
            reason: "the published-figure scan, digitizer tool/version, operator, and transform parameters are unrecorded"
                .to_string(),
        },
        final_artifact: retained,
        context_of_use: vec![ContextRange {
            name: "t_star".to_string(),
            lo: QtyAny::dimensionless(0.41),
            hi: QtyAny::dimensionless(2.95),
        }],
        license: unavailable(
            "the underlying figure is Royal Society copyright and redistribution terms for the digitized coordinates are unresolved",
        ),
        provenance: AcquisitionProvenance {
            measured_by: "J. C. Martin and W. J. Moyce; coordinate digitizer unrecorded"
                .to_string(),
            organization: "Phil. Trans. R. Soc. Lond. A 244:312-324".to_string(),
            measured_on: unavailable(
                "1952 is the publication year, not a retained acquisition date",
            ),
            source_record: "crates/fs-lbm/tests/d3q19_freesurface3.rs::lbm3_105_martin_moyce_front"
                .to_string(),
        },
        retention: RetentionPolicy {
            class: RetentionClass::Permanent,
            preserve_raw: true,
            preserve_calibration: true,
            policy_id: "frankensim-vv-corpus-permanent-v1".to_string(),
        },
        acceptance_envelopes: vec![AcceptanceRecord {
            metric: "surge-front-position-z".to_string(),
            dims: dimensionless,
            envelope: CorpusEnvelope::Unpinned {
                basis: "the live consumer applies monotone advance plus z <= 2.2*t_star + 1 for 0.5 < t_star < 2 and compares the digitized curve report-only; no scalar central band is defensible"
                    .to_string(),
            },
            regime: vec![ContextRange {
                name: "t_star".to_string(),
                lo: QtyAny::dimensionless(0.5),
                hi: QtyAny::dimensionless(2.0),
            }],
        }],
        evidence_level: EvidenceLevel::PublishedExperiment,
    }
}

fn unavailable<T>(reason: &str) -> Availability<T> {
    Availability::Unavailable {
        reason: reason.to_string(),
    }
}

fn invalid(field: DatasetField, reason: &'static str) -> CorpusError {
    CorpusError::InvalidField { field, reason }
}

fn validate_dataset(dataset: &CorpusDataset) -> Result<(), CorpusError> {
    if !valid_slug(&dataset.id) {
        return Err(invalid(
            DatasetField::Id,
            "must be a bounded lowercase ASCII slug",
        ));
    }
    validate_text(&dataset.title, DatasetField::Title)?;
    validate_artifact(dataset.raw_payload.artifact(), DatasetField::RawPayload)?;
    if let PayloadRetention::DerivedOnly { reason, .. } = &dataset.raw_payload {
        validate_text(reason, DatasetField::RawPayload)?;
    }

    if dataset.sensors.is_empty() {
        return Err(invalid(
            DatasetField::Sensors,
            "at least one sensor is required",
        ));
    }
    check_count("sensors", dataset.sensors.len(), MAX_DATASET_SENSORS)?;
    let mut sensor_ids = BTreeSet::new();
    for sensor in &dataset.sensors {
        if !valid_slug(&sensor.id) {
            return Err(invalid(
                DatasetField::Sensors,
                "sensor ids must be bounded lowercase ASCII slugs",
            ));
        }
        if !sensor_ids.insert(sensor.id.as_str()) {
            return Err(CorpusError::DuplicateSensorId {
                id: sensor.id.clone(),
            });
        }
        validate_availability_text(&sensor.instrument_id, DatasetField::Sensors)?;
        validate_text(&sensor.raw_channel, DatasetField::Sensors)?;
        match &sensor.calibration {
            Availability::Available(calibration) => validate_calibration(calibration)?,
            Availability::Unavailable { reason } => {
                validate_text(reason, DatasetField::Sensors)?;
            }
        }
        match &sensor.placement {
            Availability::Available(placement) => validate_placement(placement)?,
            Availability::Unavailable { reason } => {
                validate_text(reason, DatasetField::Sensors)?;
            }
        }
        match sensor.uncertainty {
            MeasurementUncertainty::Bounded { half_width } => {
                if half_width.dims != sensor.quantity_dims
                    || !half_width.value.is_finite()
                    || half_width.value < 0.0
                {
                    return Err(invalid(
                        DatasetField::Sensors,
                        "bounded uncertainty must be finite, non-negative, and match quantity dimensions",
                    ));
                }
            }
            MeasurementUncertainty::CovarianceDiagonal { variance } => {
                let expected = sensor
                    .quantity_dims
                    .checked_plus(sensor.quantity_dims)
                    .ok_or_else(|| {
                        invalid(
                            DatasetField::Sensors,
                            "measurement dimensions overflow when squared for covariance",
                        )
                    })?;
                if variance.dims != expected || !variance.value.is_finite() || variance.value < 0.0
                {
                    return Err(invalid(
                        DatasetField::Sensors,
                        "covariance diagonal must be finite, non-negative, and have squared quantity dimensions",
                    ));
                }
            }
            MeasurementUncertainty::Unstated => {}
        }
    }

    match &dataset.geometry {
        Availability::Available(geometry) => validate_geometry(geometry)?,
        Availability::Unavailable { reason } => validate_text(reason, DatasetField::Geometry)?,
    }
    match &dataset.environment {
        Availability::Available(conditions) => {
            if conditions.is_empty() {
                return Err(invalid(
                    DatasetField::Environment,
                    "available acquisition conditions must be nonempty",
                ));
            }
            check_count(
                "environment conditions",
                conditions.len(),
                MAX_DATASET_ITEMS,
            )?;
            validate_conditions(conditions)?;
        }
        Availability::Unavailable { reason } => {
            validate_text(reason, DatasetField::Environment)?;
        }
    }

    match &dataset.preprocessing {
        PreprocessingLineage::Complete(steps) => {
            if steps.is_empty() {
                return Err(invalid(
                    DatasetField::Preprocessing,
                    "a complete lineage needs at least a raw-import or identity transform",
                ));
            }
            check_count("preprocessing steps", steps.len(), MAX_DATASET_ITEMS)?;
            let mut expected_input = dataset.raw_payload.artifact().digest;
            for (index, step) in steps.iter().enumerate() {
                if usize::try_from(step.ordinal).ok() != Some(index) {
                    return Err(CorpusError::BrokenLineage {
                        step: index,
                        reason: "ordinals must be contiguous from zero",
                    });
                }
                validate_text(&step.operation, DatasetField::Preprocessing)?;
                validate_text(&step.version, DatasetField::Preprocessing)?;
                if step.input != expected_input {
                    return Err(CorpusError::BrokenLineage {
                        step: index,
                        reason: "input hash does not equal the preceding retained artifact",
                    });
                }
                if zero_hash(step.output) {
                    return Err(CorpusError::BrokenLineage {
                        step: index,
                        reason: "output hash is zero",
                    });
                }
                expected_input = step.output;
            }
            if dataset.final_artifact != expected_input {
                return Err(CorpusError::BrokenLineage {
                    step: steps.len(),
                    reason: "final artifact does not equal the last transform output",
                });
            }
        }
        PreprocessingLineage::Unreplayable {
            retained_input,
            retained_output,
            reason,
        } => {
            validate_text(reason, DatasetField::Preprocessing)?;
            if *retained_input != dataset.raw_payload.artifact().digest {
                return Err(CorpusError::BrokenLineage {
                    step: 0,
                    reason: "unreplayable lineage input must bind the retained payload",
                });
            }
            if zero_hash(*retained_output) || *retained_output != dataset.final_artifact {
                return Err(CorpusError::BrokenLineage {
                    step: 0,
                    reason: "unreplayable lineage output must bind the final artifact",
                });
            }
        }
    }

    if dataset.context_of_use.is_empty() {
        return Err(invalid(
            DatasetField::ContextOfUse,
            "at least one bounded context coordinate is required",
        ));
    }
    check_count(
        "context ranges",
        dataset.context_of_use.len(),
        MAX_DATASET_ITEMS,
    )?;
    validate_ranges(&dataset.context_of_use, DatasetField::ContextOfUse)?;

    match &dataset.license {
        Availability::Available(license) => {
            validate_text(&license.identifier, DatasetField::License)?;
            validate_text(&license.terms, DatasetField::License)?;
        }
        Availability::Unavailable { reason } => validate_text(reason, DatasetField::License)?,
    }
    validate_text(&dataset.provenance.measured_by, DatasetField::Provenance)?;
    validate_text(&dataset.provenance.organization, DatasetField::Provenance)?;
    validate_text(&dataset.provenance.source_record, DatasetField::Provenance)?;
    match &dataset.provenance.measured_on {
        Availability::Available(date) if valid_date(date) => {}
        Availability::Available(_) => {
            return Err(invalid(
                DatasetField::Provenance,
                "available measured_on must be a real YYYY-MM-DD date",
            ));
        }
        Availability::Unavailable { reason } => {
            validate_text(reason, DatasetField::Provenance)?;
        }
    }

    validate_text(&dataset.retention.policy_id, DatasetField::Retention)?;
    if !dataset.retention.preserve_raw || !dataset.retention.preserve_calibration {
        return Err(invalid(
            DatasetField::Retention,
            "raw payloads and calibration records must be retained together",
        ));
    }
    if matches!(dataset.retention.class, RetentionClass::Years(0)) {
        return Err(invalid(
            DatasetField::Retention,
            "a finite retention period must be at least one year",
        ));
    }

    if dataset.acceptance_envelopes.is_empty() {
        return Err(invalid(
            DatasetField::AcceptanceEnvelopes,
            "at least one metric/envelope/regime record is required",
        ));
    }
    check_count(
        "acceptance envelopes",
        dataset.acceptance_envelopes.len(),
        MAX_DATASET_ITEMS,
    )?;
    let mut metrics = BTreeSet::new();
    for acceptance in &dataset.acceptance_envelopes {
        if !valid_slug(&acceptance.metric) {
            return Err(invalid(
                DatasetField::AcceptanceEnvelopes,
                "metric names must be bounded lowercase ASCII slugs",
            ));
        }
        if !metrics.insert(acceptance.metric.as_str()) {
            return Err(CorpusError::DuplicateName {
                collection: "acceptance metric",
                name: acceptance.metric.clone(),
            });
        }
        match &acceptance.envelope {
            CorpusEnvelope::Tolerance { atol, rtol }
                if atol.is_finite() && rtol.is_finite() && *atol >= 0.0 && *rtol >= 0.0 => {}
            CorpusEnvelope::Interval { lo, hi } if lo.is_finite() && hi.is_finite() && lo <= hi => {
            }
            CorpusEnvelope::Unpinned { basis } => {
                validate_text(basis, DatasetField::AcceptanceEnvelopes)?;
            }
            _ => {
                return Err(invalid(
                    DatasetField::AcceptanceEnvelopes,
                    "envelope bounds must be finite, ordered, and non-negative where applicable",
                ));
            }
        }
        if acceptance.regime.is_empty() {
            return Err(invalid(
                DatasetField::AcceptanceEnvelopes,
                "each envelope needs an explicit nonempty regime",
            ));
        }
        validate_ranges(&acceptance.regime, DatasetField::AcceptanceEnvelopes)?;
        for range in &acceptance.regime {
            let outer = dataset
                .context_of_use
                .iter()
                .find(|candidate| candidate.name == range.name)
                .ok_or_else(|| {
                    invalid(
                        DatasetField::AcceptanceEnvelopes,
                        "envelope regime names must exist in context_of_use",
                    )
                })?;
            if range.lo.dims != outer.lo.dims
                || range.lo.value < outer.lo.value
                || range.hi.value > outer.hi.value
            {
                return Err(invalid(
                    DatasetField::AcceptanceEnvelopes,
                    "envelope regime must be dimension-compatible and contained in context_of_use",
                ));
            }
        }
    }

    if dataset.encode().len() > MAX_DATASET_CANONICAL_BYTES {
        return Err(CorpusError::ResourceLimit {
            resource: "canonical dataset bytes",
            limit: MAX_DATASET_CANONICAL_BYTES,
            observed: dataset.encode().len(),
        });
    }
    Ok(())
}

fn validate_calibration(calibration: &CalibrationRecord) -> Result<(), CorpusError> {
    validate_text(&calibration.certificate_id, DatasetField::Sensors)?;
    if zero_hash(calibration.certificate_hash) {
        return Err(invalid(
            DatasetField::Sensors,
            "calibration certificate hash must be nonzero",
        ));
    }
    if !valid_date(&calibration.issued_on)
        || calibration
            .valid_through
            .as_deref()
            .is_some_and(|date| !valid_date(date) || date < calibration.issued_on.as_str())
    {
        return Err(invalid(
            DatasetField::Sensors,
            "calibration dates must be ordered real YYYY-MM-DD dates",
        ));
    }
    Ok(())
}

fn validate_availability_text(
    value: &Availability<String>,
    field: DatasetField,
) -> Result<(), CorpusError> {
    match value {
        Availability::Available(text) | Availability::Unavailable { reason: text } => {
            validate_text(text, field)
        }
    }
}

fn validate_placement(placement: &SensorPlacement) -> Result<(), CorpusError> {
    validate_text(&placement.frame, DatasetField::Sensors)?;
    let length = Dims([1, 0, 0, 0, 0, 0]);
    for coordinate in placement.coordinates {
        if coordinate.dims != length || !coordinate.value.is_finite() {
            return Err(invalid(
                DatasetField::Sensors,
                "placement coordinates must be finite lengths",
            ));
        }
    }
    for uncertainty in placement.uncertainty {
        if uncertainty.dims != length || !uncertainty.value.is_finite() || uncertainty.value < 0.0 {
            return Err(invalid(
                DatasetField::Sensors,
                "placement uncertainty must be a finite non-negative length",
            ));
        }
    }
    Ok(())
}

fn validate_geometry(geometry: &GeometryRecord) -> Result<(), CorpusError> {
    validate_artifact(&geometry.nominal, DatasetField::Geometry)?;
    if let Some(as_built) = &geometry.as_built {
        validate_artifact(as_built, DatasetField::Geometry)?;
    }
    validate_text(&geometry.frame, DatasetField::Geometry)
}

fn validate_conditions(conditions: &[EnvironmentCondition]) -> Result<(), CorpusError> {
    let mut names = BTreeSet::new();
    for condition in conditions {
        validate_text(&condition.name, DatasetField::Environment)?;
        if !names.insert(condition.name.as_str()) {
            return Err(CorpusError::DuplicateName {
                collection: "environment condition",
                name: condition.name.clone(),
            });
        }
        if condition.value.dims != condition.uncertainty.dims
            || !condition.value.value.is_finite()
            || !condition.uncertainty.value.is_finite()
            || condition.uncertainty.value < 0.0
        {
            return Err(invalid(
                DatasetField::Environment,
                "condition value/uncertainty must be finite, dimension-compatible, and non-negative",
            ));
        }
    }
    Ok(())
}

fn validate_ranges(ranges: &[ContextRange], field: DatasetField) -> Result<(), CorpusError> {
    check_count("ranges", ranges.len(), MAX_DATASET_ITEMS)?;
    let mut names = BTreeSet::new();
    for range in ranges {
        validate_text(&range.name, field)?;
        if !names.insert(range.name.as_str()) {
            return Err(CorpusError::DuplicateName {
                collection: "context range",
                name: range.name.clone(),
            });
        }
        if range.lo.dims != range.hi.dims
            || !range.lo.value.is_finite()
            || !range.hi.value.is_finite()
            || range.lo.value > range.hi.value
        {
            return Err(invalid(
                field,
                "context ranges must have finite ordered bounds with identical dimensions",
            ));
        }
    }
    Ok(())
}

fn validate_artifact(artifact: &CorpusArtifact, field: DatasetField) -> Result<(), CorpusError> {
    if zero_hash(artifact.digest) || artifact.byte_len == 0 {
        return Err(invalid(
            field,
            "artifact hash must be nonzero and byte length must be positive",
        ));
    }
    validate_text(&artifact.media_type, field)?;
    validate_text(&artifact.locator, field)?;
    if artifact.locator.contains('\\')
        || std::path::Path::new(&artifact.locator)
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(invalid(
            field,
            "artifact locator must be a normalized relative path",
        ));
    }
    Ok(())
}

fn validate_text(text: &str, field: DatasetField) -> Result<(), CorpusError> {
    if text.trim().is_empty() {
        return Err(invalid(field, "text must not be blank"));
    }
    if text.len() > MAX_CORPUS_TEXT_BYTES {
        return Err(CorpusError::ResourceLimit {
            resource: field.name(),
            limit: MAX_CORPUS_TEXT_BYTES,
            observed: text.len(),
        });
    }
    if text.chars().any(char::is_control) {
        return Err(invalid(field, "text must not contain control characters"));
    }
    Ok(())
}

fn check_count(resource: &'static str, observed: usize, limit: usize) -> Result<(), CorpusError> {
    if observed > limit {
        Err(CorpusError::ResourceLimit {
            resource,
            limit,
            observed,
        })
    } else {
        Ok(())
    }
}

fn valid_slug(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= MAX_CORPUS_TEXT_BYTES
        && value.as_bytes()[0].is_ascii_lowercase()
        && value.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-' || byte == b'_'
        })
}

fn valid_date(value: &str) -> bool {
    let bytes = value.as_bytes();
    if bytes.len() != 10
        || bytes[4] != b'-'
        || bytes[7] != b'-'
        || bytes
            .iter()
            .enumerate()
            .any(|(index, byte)| index != 4 && index != 7 && !byte.is_ascii_digit())
    {
        return false;
    }
    let parse = |range: std::ops::Range<usize>| value[range].parse::<u32>().ok();
    let (Some(year), Some(month), Some(day)) = (parse(0..4), parse(5..7), parse(8..10)) else {
        return false;
    };
    if year == 0 || !(1..=12).contains(&month) {
        return false;
    }
    let leap = year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400));
    let max_day = match month {
        2 if leap => 29,
        2 => 28,
        4 | 6 | 9 | 11 => 30,
        _ => 31,
    };
    (1..=max_day).contains(&day)
}

fn zero_hash(hash: ContentHash) -> bool {
    hash.as_bytes().iter().all(|byte| *byte == 0)
}

fn validate_query_context(
    dataset: &CorpusDataset,
    context: &[ContextValue],
) -> Result<(), CorpusQueryRefusal> {
    let mut seen = BTreeSet::new();
    for coordinate in context {
        if !seen.insert(coordinate.name.as_str()) {
            return Err(CorpusQueryRefusal::DuplicateContext {
                name: coordinate.name.clone(),
            });
        }
        let Some(range) = dataset
            .context_of_use
            .iter()
            .find(|range| range.name == coordinate.name)
        else {
            return Err(CorpusQueryRefusal::UnknownContext {
                name: coordinate.name.clone(),
            });
        };
        if coordinate.value.dims != range.lo.dims {
            return Err(CorpusQueryRefusal::ContextDimensionMismatch {
                name: coordinate.name.clone(),
                expected: range.lo.dims,
                observed: coordinate.value.dims,
            });
        }
        if !coordinate.value.value.is_finite()
            || coordinate.value.value < range.lo.value
            || coordinate.value.value > range.hi.value
        {
            return Err(CorpusQueryRefusal::OutOfContext {
                name: coordinate.name.clone(),
                value: coordinate.value.value,
                lo: range.lo.value,
                hi: range.hi.value,
            });
        }
    }
    for range in &dataset.context_of_use {
        if !seen.contains(range.name.as_str()) {
            return Err(CorpusQueryRefusal::MissingContext {
                name: range.name.clone(),
            });
        }
    }
    Ok(())
}

impl CorpusDataset {
    /// Canonical, length-framed binary serialization.
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&CORPUS_SCHEMA_VERSION.to_le_bytes());
        push_text(&mut out, &self.id);
        push_text(&mut out, &self.title);
        match &self.raw_payload {
            PayloadRetention::OriginalRaw(artifact) => {
                out.push(1);
                push_artifact(&mut out, artifact);
            }
            PayloadRetention::DerivedOnly { retained, reason } => {
                out.push(2);
                push_artifact(&mut out, retained);
                push_text(&mut out, reason);
            }
        }

        push_u32(&mut out, self.sensors.len());
        for sensor in &self.sensors {
            push_text(&mut out, &sensor.id);
            push_availability_text(&mut out, &sensor.instrument_id);
            push_text(&mut out, &sensor.raw_channel);
            push_dims(&mut out, sensor.quantity_dims);
            push_availability_calibration(&mut out, &sensor.calibration);
            push_availability_placement(&mut out, &sensor.placement);
            match &sensor.uncertainty {
                MeasurementUncertainty::Bounded { half_width } => {
                    out.push(1);
                    push_qty(&mut out, *half_width);
                }
                MeasurementUncertainty::CovarianceDiagonal { variance } => {
                    out.push(2);
                    push_qty(&mut out, *variance);
                }
                MeasurementUncertainty::Unstated => out.push(3),
            }
        }

        match &self.geometry {
            Availability::Available(geometry) => {
                out.push(1);
                push_artifact(&mut out, &geometry.nominal);
                match &geometry.as_built {
                    Some(artifact) => {
                        out.push(1);
                        push_artifact(&mut out, artifact);
                    }
                    None => out.push(0),
                }
                push_text(&mut out, &geometry.frame);
            }
            Availability::Unavailable { reason } => {
                out.push(2);
                push_text(&mut out, reason);
            }
        }

        match &self.environment {
            Availability::Available(conditions) => {
                out.push(1);
                push_u32(&mut out, conditions.len());
                for condition in conditions {
                    push_text(&mut out, &condition.name);
                    push_qty(&mut out, condition.value);
                    push_qty(&mut out, condition.uncertainty);
                }
            }
            Availability::Unavailable { reason } => {
                out.push(2);
                push_text(&mut out, reason);
            }
        }

        out.push(partition_tag(self.partition));
        match &self.preprocessing {
            PreprocessingLineage::Complete(steps) => {
                out.push(1);
                push_u32(&mut out, steps.len());
                for step in steps {
                    out.extend_from_slice(&step.ordinal.to_le_bytes());
                    push_text(&mut out, &step.operation);
                    push_text(&mut out, &step.version);
                    push_hash(&mut out, step.input);
                    push_hash(&mut out, step.output);
                }
            }
            PreprocessingLineage::Unreplayable {
                retained_input,
                retained_output,
                reason,
            } => {
                out.push(2);
                push_hash(&mut out, *retained_input);
                push_hash(&mut out, *retained_output);
                push_text(&mut out, reason);
            }
        }
        push_hash(&mut out, self.final_artifact);

        push_ranges(&mut out, &self.context_of_use);
        match &self.license {
            Availability::Available(license) => {
                out.push(1);
                push_text(&mut out, &license.identifier);
                push_text(&mut out, &license.terms);
                out.push(redistribution_tag(license.redistribution));
            }
            Availability::Unavailable { reason } => {
                out.push(2);
                push_text(&mut out, reason);
            }
        }
        push_text(&mut out, &self.provenance.measured_by);
        push_text(&mut out, &self.provenance.organization);
        push_availability_text(&mut out, &self.provenance.measured_on);
        push_text(&mut out, &self.provenance.source_record);
        match self.retention.class {
            RetentionClass::Permanent => out.push(1),
            RetentionClass::Years(years) => {
                out.push(2);
                out.extend_from_slice(&years.to_le_bytes());
            }
        }
        out.push(u8::from(self.retention.preserve_raw));
        out.push(u8::from(self.retention.preserve_calibration));
        push_text(&mut out, &self.retention.policy_id);

        push_u32(&mut out, self.acceptance_envelopes.len());
        for acceptance in &self.acceptance_envelopes {
            push_text(&mut out, &acceptance.metric);
            push_dims(&mut out, acceptance.dims);
            match &acceptance.envelope {
                CorpusEnvelope::Tolerance { atol, rtol } => {
                    out.push(1);
                    push_f64(&mut out, *atol);
                    push_f64(&mut out, *rtol);
                }
                CorpusEnvelope::Interval { lo, hi } => {
                    out.push(2);
                    push_f64(&mut out, *lo);
                    push_f64(&mut out, *hi);
                }
                CorpusEnvelope::Unpinned { basis } => {
                    out.push(3);
                    push_text(&mut out, basis);
                }
            }
            push_ranges(&mut out, &acceptance.regime);
        }
        out.push(evidence_tag(self.evidence_level));
        out
    }

    /// Decode, validate, canonicalize, and byte-compare one dataset.
    pub fn decode(bytes: &[u8]) -> Result<Self, CorpusError> {
        if bytes.len() > MAX_DATASET_CANONICAL_BYTES {
            return Err(CorpusError::ResourceLimit {
                resource: "canonical dataset bytes",
                limit: MAX_DATASET_CANONICAL_BYTES,
                observed: bytes.len(),
            });
        }
        let mut reader = Reader::new(bytes);
        if reader.take(MAGIC.len())? != MAGIC {
            return Err(CorpusError::BadMagic);
        }
        let schema = reader.u32()?;
        if schema != CORPUS_SCHEMA_VERSION {
            return Err(CorpusError::UnsupportedSchema { observed: schema });
        }

        let id = reader.text()?;
        let title = reader.text()?;
        let raw_payload = match reader.u8()? {
            1 => PayloadRetention::OriginalRaw(reader.artifact()?),
            2 => PayloadRetention::DerivedOnly {
                retained: reader.artifact()?,
                reason: reader.text()?,
            },
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "payload retention",
                    tag,
                });
            }
        };
        let sensor_count = reader.count("sensors", MAX_DATASET_SENSORS)?;
        let mut sensors = Vec::with_capacity(sensor_count);
        for _ in 0..sensor_count {
            let sensor_id = reader.text()?;
            let instrument_id = reader.availability_text("sensor instrument")?;
            let raw_channel = reader.text()?;
            let quantity_dims = reader.dims()?;
            let calibration = match reader.u8()? {
                1 => Availability::Available(CalibrationRecord {
                    certificate_id: reader.text()?,
                    certificate_hash: reader.hash()?,
                    issued_on: reader.text()?,
                    valid_through: reader.option_text()?,
                }),
                2 => Availability::Unavailable {
                    reason: reader.text()?,
                },
                tag => {
                    return Err(CorpusError::InvalidTag {
                        kind: "sensor calibration availability",
                        tag,
                    });
                }
            };
            let placement = match reader.u8()? {
                1 => Availability::Available(SensorPlacement {
                    frame: reader.text()?,
                    coordinates: [reader.qty()?, reader.qty()?, reader.qty()?],
                    uncertainty: [reader.qty()?, reader.qty()?, reader.qty()?],
                }),
                2 => Availability::Unavailable {
                    reason: reader.text()?,
                },
                tag => {
                    return Err(CorpusError::InvalidTag {
                        kind: "sensor placement availability",
                        tag,
                    });
                }
            };
            let uncertainty = match reader.u8()? {
                1 => MeasurementUncertainty::Bounded {
                    half_width: reader.qty()?,
                },
                2 => MeasurementUncertainty::CovarianceDiagonal {
                    variance: reader.qty()?,
                },
                3 => MeasurementUncertainty::Unstated,
                tag => {
                    return Err(CorpusError::InvalidTag {
                        kind: "measurement uncertainty",
                        tag,
                    });
                }
            };
            sensors.push(SensorRecord {
                id: sensor_id,
                instrument_id,
                raw_channel,
                quantity_dims,
                calibration,
                placement,
                uncertainty,
            });
        }

        let geometry = match reader.u8()? {
            1 => {
                let nominal = reader.artifact()?;
                let as_built = match reader.u8()? {
                    0 => None,
                    1 => Some(reader.artifact()?),
                    tag => {
                        return Err(CorpusError::InvalidTag {
                            kind: "optional as-built artifact",
                            tag,
                        });
                    }
                };
                Availability::Available(GeometryRecord {
                    nominal,
                    as_built,
                    frame: reader.text()?,
                })
            }
            2 => Availability::Unavailable {
                reason: reader.text()?,
            },
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "geometry availability",
                    tag,
                });
            }
        };

        let environment = match reader.u8()? {
            1 => {
                let count = reader.count("environment conditions", MAX_DATASET_ITEMS)?;
                let mut conditions = Vec::with_capacity(count);
                for _ in 0..count {
                    conditions.push(EnvironmentCondition {
                        name: reader.text()?,
                        value: reader.qty()?,
                        uncertainty: reader.qty()?,
                    });
                }
                Availability::Available(conditions)
            }
            2 => Availability::Unavailable {
                reason: reader.text()?,
            },
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "environment availability",
                    tag,
                });
            }
        };

        let partition = parse_partition(reader.u8()?)?;
        let preprocessing = match reader.u8()? {
            1 => {
                let count = reader.count("preprocessing steps", MAX_DATASET_ITEMS)?;
                let mut steps = Vec::with_capacity(count);
                for _ in 0..count {
                    steps.push(PreprocessingStep {
                        ordinal: reader.u32()?,
                        operation: reader.text()?,
                        version: reader.text()?,
                        input: reader.hash()?,
                        output: reader.hash()?,
                    });
                }
                PreprocessingLineage::Complete(steps)
            }
            2 => PreprocessingLineage::Unreplayable {
                retained_input: reader.hash()?,
                retained_output: reader.hash()?,
                reason: reader.text()?,
            },
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "preprocessing lineage",
                    tag,
                });
            }
        };
        let final_artifact = reader.hash()?;
        let context_of_use = reader.ranges()?;
        let license = match reader.u8()? {
            1 => Availability::Available(CorpusLicense {
                identifier: reader.text()?,
                terms: reader.text()?,
                redistribution: parse_redistribution(reader.u8()?)?,
            }),
            2 => Availability::Unavailable {
                reason: reader.text()?,
            },
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "license availability",
                    tag,
                });
            }
        };
        let provenance = AcquisitionProvenance {
            measured_by: reader.text()?,
            organization: reader.text()?,
            measured_on: reader.availability_text("acquisition date")?,
            source_record: reader.text()?,
        };
        let retention_class = match reader.u8()? {
            1 => RetentionClass::Permanent,
            2 => RetentionClass::Years(reader.u16()?),
            tag => {
                return Err(CorpusError::InvalidTag {
                    kind: "retention class",
                    tag,
                });
            }
        };
        let preserve_raw = reader.bool("retention preserve_raw")?;
        let preserve_calibration = reader.bool("retention preserve_calibration")?;
        let retention = RetentionPolicy {
            class: retention_class,
            preserve_raw,
            preserve_calibration,
            policy_id: reader.text()?,
        };

        let acceptance_count = reader.count("acceptance envelopes", MAX_DATASET_ITEMS)?;
        let mut acceptance_envelopes = Vec::with_capacity(acceptance_count);
        for _ in 0..acceptance_count {
            let metric = reader.text()?;
            let dims = reader.dims()?;
            let envelope = match reader.u8()? {
                1 => CorpusEnvelope::Tolerance {
                    atol: reader.f64()?,
                    rtol: reader.f64()?,
                },
                2 => CorpusEnvelope::Interval {
                    lo: reader.f64()?,
                    hi: reader.f64()?,
                },
                3 => CorpusEnvelope::Unpinned {
                    basis: reader.text()?,
                },
                tag => {
                    return Err(CorpusError::InvalidTag {
                        kind: "acceptance envelope",
                        tag,
                    });
                }
            };
            acceptance_envelopes.push(AcceptanceRecord {
                metric,
                dims,
                envelope,
                regime: reader.ranges()?,
            });
        }
        let evidence_level = parse_evidence(reader.u8()?)?;
        if reader.remaining() != 0 {
            return Err(CorpusError::TrailingBytes {
                count: reader.remaining(),
            });
        }

        let dataset = admit_dataset(DatasetDraft {
            id: Some(id),
            title: Some(title),
            raw_payload: Some(raw_payload),
            sensors: Some(sensors),
            geometry: Some(geometry),
            environment: Some(environment),
            partition: Some(partition),
            preprocessing: Some(preprocessing),
            final_artifact: Some(final_artifact),
            context_of_use: Some(context_of_use),
            license: Some(license),
            provenance: Some(provenance),
            retention: Some(retention),
            acceptance_envelopes: Some(acceptance_envelopes),
            evidence_level: Some(evidence_level),
        })?;
        if dataset.encode() != bytes {
            return Err(CorpusError::NonCanonicalEncoding);
        }
        Ok(dataset)
    }
}

fn push_u32(out: &mut Vec<u8>, value: usize) {
    out.extend_from_slice(&(value as u32).to_le_bytes());
}

fn push_text(out: &mut Vec<u8>, value: &str) {
    push_u32(out, value.len());
    out.extend_from_slice(value.as_bytes());
}

fn push_option_text(out: &mut Vec<u8>, value: Option<&str>) {
    match value {
        Some(value) => {
            out.push(1);
            push_text(out, value);
        }
        None => out.push(0),
    }
}

fn push_availability_text(out: &mut Vec<u8>, value: &Availability<String>) {
    match value {
        Availability::Available(text) => {
            out.push(1);
            push_text(out, text);
        }
        Availability::Unavailable { reason } => {
            out.push(2);
            push_text(out, reason);
        }
    }
}

fn push_availability_calibration(out: &mut Vec<u8>, value: &Availability<CalibrationRecord>) {
    match value {
        Availability::Available(calibration) => {
            out.push(1);
            push_text(out, &calibration.certificate_id);
            push_hash(out, calibration.certificate_hash);
            push_text(out, &calibration.issued_on);
            push_option_text(out, calibration.valid_through.as_deref());
        }
        Availability::Unavailable { reason } => {
            out.push(2);
            push_text(out, reason);
        }
    }
}

fn push_availability_placement(out: &mut Vec<u8>, value: &Availability<SensorPlacement>) {
    match value {
        Availability::Available(placement) => {
            out.push(1);
            push_text(out, &placement.frame);
            for coordinate in placement.coordinates {
                push_qty(out, coordinate);
            }
            for uncertainty in placement.uncertainty {
                push_qty(out, uncertainty);
            }
        }
        Availability::Unavailable { reason } => {
            out.push(2);
            push_text(out, reason);
        }
    }
}

fn push_hash(out: &mut Vec<u8>, hash: ContentHash) {
    out.extend_from_slice(hash.as_bytes());
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn push_dims(out: &mut Vec<u8>, dims: Dims) {
    for exponent in dims.0 {
        out.push(exponent as u8);
    }
}

fn push_qty(out: &mut Vec<u8>, quantity: QtyAny) {
    push_f64(out, quantity.value);
    push_dims(out, quantity.dims);
}

fn push_artifact(out: &mut Vec<u8>, artifact: &CorpusArtifact) {
    push_hash(out, artifact.digest);
    out.extend_from_slice(&artifact.byte_len.to_le_bytes());
    push_text(out, &artifact.media_type);
    push_text(out, &artifact.locator);
}

fn push_ranges(out: &mut Vec<u8>, ranges: &[ContextRange]) {
    push_u32(out, ranges.len());
    for range in ranges {
        push_text(out, &range.name);
        push_qty(out, range.lo);
        push_qty(out, range.hi);
    }
}

const fn partition_tag(partition: DatasetPartition) -> u8 {
    match partition {
        DatasetPartition::Training => 1,
        DatasetPartition::Calibration => 2,
        DatasetPartition::Validation => 3,
    }
}

fn parse_partition(tag: u8) -> Result<DatasetPartition, CorpusError> {
    match tag {
        1 => Ok(DatasetPartition::Training),
        2 => Ok(DatasetPartition::Calibration),
        3 => Ok(DatasetPartition::Validation),
        tag => Err(CorpusError::InvalidTag {
            kind: "dataset partition",
            tag,
        }),
    }
}

const fn redistribution_tag(policy: RedistributionPolicy) -> u8 {
    match policy {
        RedistributionPolicy::Allowed => 1,
        RedistributionPolicy::MetadataOnly => 2,
        RedistributionPolicy::Prohibited => 3,
    }
}

fn parse_redistribution(tag: u8) -> Result<RedistributionPolicy, CorpusError> {
    match tag {
        1 => Ok(RedistributionPolicy::Allowed),
        2 => Ok(RedistributionPolicy::MetadataOnly),
        3 => Ok(RedistributionPolicy::Prohibited),
        tag => Err(CorpusError::InvalidTag {
            kind: "redistribution policy",
            tag,
        }),
    }
}

const fn evidence_tag(level: EvidenceLevel) -> u8 {
    match level {
        EvidenceLevel::Analytic => 1,
        EvidenceLevel::CrossCode => 2,
        EvidenceLevel::PublishedExperiment => 3,
        EvidenceLevel::Blind => 4,
        EvidenceLevel::Field => 5,
    }
}

fn parse_evidence(tag: u8) -> Result<EvidenceLevel, CorpusError> {
    match tag {
        1 => Ok(EvidenceLevel::Analytic),
        2 => Ok(EvidenceLevel::CrossCode),
        3 => Ok(EvidenceLevel::PublishedExperiment),
        4 => Ok(EvidenceLevel::Blind),
        5 => Ok(EvidenceLevel::Field),
        tag => Err(CorpusError::InvalidTag {
            kind: "evidence level",
            tag,
        }),
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.cursor)
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], CorpusError> {
        let end = self
            .cursor
            .checked_add(count)
            .ok_or(CorpusError::Truncated)?;
        let value = self
            .bytes
            .get(self.cursor..end)
            .ok_or(CorpusError::Truncated)?;
        self.cursor = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, CorpusError> {
        Ok(self.take(1)?[0])
    }

    fn bool(&mut self, kind: &'static str) -> Result<bool, CorpusError> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            tag => Err(CorpusError::InvalidTag { kind, tag }),
        }
    }

    fn u16(&mut self) -> Result<u16, CorpusError> {
        let bytes: [u8; 2] = self
            .take(2)?
            .try_into()
            .map_err(|_| CorpusError::Truncated)?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self) -> Result<u32, CorpusError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| CorpusError::Truncated)?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, CorpusError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| CorpusError::Truncated)?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn f64(&mut self) -> Result<f64, CorpusError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn count(&mut self, resource: &'static str, limit: usize) -> Result<usize, CorpusError> {
        let observed = self.u32()? as usize;
        check_count(resource, observed, limit)?;
        Ok(observed)
    }

    fn text(&mut self) -> Result<String, CorpusError> {
        let len = self.count("string bytes", MAX_CORPUS_TEXT_BYTES)?;
        let bytes = self.take(len)?;
        let value = std::str::from_utf8(bytes).map_err(|_| CorpusError::InvalidUtf8)?;
        Ok(value.to_string())
    }

    fn option_text(&mut self) -> Result<Option<String>, CorpusError> {
        match self.u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.text()?)),
            tag => Err(CorpusError::InvalidTag {
                kind: "optional string",
                tag,
            }),
        }
    }

    fn availability_text(
        &mut self,
        kind: &'static str,
    ) -> Result<Availability<String>, CorpusError> {
        match self.u8()? {
            1 => Ok(Availability::Available(self.text()?)),
            2 => Ok(Availability::Unavailable {
                reason: self.text()?,
            }),
            tag => Err(CorpusError::InvalidTag { kind, tag }),
        }
    }

    fn hash(&mut self) -> Result<ContentHash, CorpusError> {
        let bytes: [u8; 32] = self
            .take(32)?
            .try_into()
            .map_err(|_| CorpusError::Truncated)?;
        Ok(ContentHash(bytes))
    }

    fn dims(&mut self) -> Result<Dims, CorpusError> {
        let mut dims = [0_i8; 6];
        for exponent in &mut dims {
            *exponent = self.u8()? as i8;
        }
        Ok(Dims(dims))
    }

    fn qty(&mut self) -> Result<QtyAny, CorpusError> {
        Ok(QtyAny::new(self.f64()?, self.dims()?))
    }

    fn artifact(&mut self) -> Result<CorpusArtifact, CorpusError> {
        Ok(CorpusArtifact {
            digest: self.hash()?,
            byte_len: self.u64()?,
            media_type: self.text()?,
            locator: self.text()?,
        })
    }

    fn ranges(&mut self) -> Result<Vec<ContextRange>, CorpusError> {
        let count = self.count("ranges", MAX_DATASET_ITEMS)?;
        let mut ranges = Vec::with_capacity(count);
        for _ in 0..count {
            ranges.push(ContextRange {
                name: self.text()?,
                lo: self.qty()?,
                hi: self.qty()?,
            });
        }
        Ok(ranges)
    }
}

/// One deterministic corpus-audit row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusAuditRow {
    dataset_id: String,
    mandatory_present: usize,
    mandatory_total: usize,
    optional_present: usize,
    optional_total: usize,
    partition: DatasetPartition,
    evidence_level: EvidenceLevel,
    physical_cap: ColorRank,
    status: &'static str,
}

impl CorpusAuditRow {
    /// Audited dataset id.
    #[must_use]
    pub fn dataset_id(&self) -> &str {
        &self.dataset_id
    }

    /// Number of mandatory fields present and valid.
    #[must_use]
    pub const fn mandatory_present(&self) -> usize {
        self.mandatory_present
    }

    /// Total mandatory-field count for this schema.
    #[must_use]
    pub const fn mandatory_total(&self) -> usize {
        self.mandatory_total
    }

    /// Stable `OK`, `WARN`, or `ERROR` row status.
    #[must_use]
    pub const fn status(&self) -> &'static str {
        self.status
    }
}

/// Complete deterministic audit, including warn-level optional gaps.
#[derive(Debug, Clone, PartialEq)]
pub struct CorpusAuditReport {
    rows: Vec<CorpusAuditRow>,
    warnings: Vec<String>,
    errors: Vec<String>,
}

impl CorpusAuditReport {
    /// Sorted per-dataset audit rows.
    #[must_use]
    pub fn rows(&self) -> &[CorpusAuditRow] {
        &self.rows
    }

    /// Stable warn-level gap diagnostics.
    #[must_use]
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Stable validation errors.
    #[must_use]
    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    /// Whether the audit contains no structural errors.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.errors.is_empty()
    }

    /// Stable field-completeness table followed by structured WARN/ERROR
    /// lines. Warnings do not make the audit fail.
    #[must_use]
    pub fn render_table(&self) -> String {
        use std::fmt::Write as _;

        let mut out = String::from(
            "dataset_id | mandatory | optional | partition | evidence | physical_cap | status\n",
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "{} | {}/{} | {}/{} | {} | {} | {} | {}",
                row.dataset_id,
                row.mandatory_present,
                row.mandatory_total,
                row.optional_present,
                row.optional_total,
                row.partition.name(),
                row.evidence_level.code(),
                color_name(row.physical_cap),
                row.status
            );
        }
        for warning in &self.warnings {
            let _ = writeln!(out, "level=WARN {warning}");
        }
        for error in &self.errors {
            let _ = writeln!(out, "level=ERROR {error}");
        }
        out
    }
}

impl CorpusRegistry {
    /// Revalidate every dataset and report mandatory completeness plus
    /// optional gaps in stable dataset-id order.
    #[must_use]
    pub fn audit(&self) -> CorpusAuditReport {
        let mut rows = Vec::with_capacity(self.datasets.len());
        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        for dataset in &self.datasets {
            let validation = validate_dataset(dataset);
            let validation_ok = if let Err(error) = validation {
                errors.push(format!("dataset={} error={error}", dataset.id));
                false
            } else {
                true
            };
            let warning_start = warnings.len();

            if let PayloadRetention::DerivedOnly { reason, .. } = &dataset.raw_payload {
                warnings.push(format!(
                    "dataset={} claim_gap=raw_payload.original reason={reason}",
                    dataset.id
                ));
            }
            if let Availability::Unavailable { reason } = &dataset.geometry {
                warnings.push(format!(
                    "dataset={} claim_gap=geometry.nominal reason={reason}",
                    dataset.id
                ));
            }
            for sensor in &dataset.sensors {
                if let Availability::Unavailable { reason } = &sensor.instrument_id {
                    warnings.push(format!(
                        "dataset={} claim_gap=sensors.{}.instrument_id reason={reason}",
                        dataset.id, sensor.id
                    ));
                }
                if let Availability::Unavailable { reason } = &sensor.calibration {
                    warnings.push(format!(
                        "dataset={} claim_gap=sensors.{}.calibration reason={reason}",
                        dataset.id, sensor.id
                    ));
                }
                if let Availability::Unavailable { reason } = &sensor.placement {
                    warnings.push(format!(
                        "dataset={} claim_gap=sensors.{}.placement reason={reason}",
                        dataset.id, sensor.id
                    ));
                }
                if sensor.uncertainty == MeasurementUncertainty::Unstated {
                    warnings.push(format!(
                        "dataset={} claim_gap=sensors.{}.uncertainty",
                        dataset.id, sensor.id
                    ));
                }
            }
            if let Availability::Unavailable { reason } = &dataset.environment {
                warnings.push(format!(
                    "dataset={} claim_gap=environment reason={reason}",
                    dataset.id
                ));
            }
            if let PreprocessingLineage::Unreplayable { reason, .. } = &dataset.preprocessing {
                warnings.push(format!(
                    "dataset={} claim_gap=preprocessing.replay reason={reason}",
                    dataset.id
                ));
            }
            if let Availability::Unavailable { reason } = &dataset.license {
                warnings.push(format!(
                    "dataset={} claim_gap=license reason={reason}",
                    dataset.id
                ));
            }
            if let Availability::Unavailable { reason } = &dataset.provenance.measured_on {
                warnings.push(format!(
                    "dataset={} claim_gap=provenance.measured_on reason={reason}",
                    dataset.id
                ));
            }
            for acceptance in &dataset.acceptance_envelopes {
                if let CorpusEnvelope::Unpinned { basis } = &acceptance.envelope {
                    warnings.push(format!(
                        "dataset={} claim_gap=acceptance.{} reason={basis}",
                        dataset.id, acceptance.metric
                    ));
                }
            }

            let mut optional_present = 0;
            let optional_total = 2;
            if matches!(
                &dataset.geometry,
                Availability::Available(geometry) if geometry.as_built.is_some()
            ) {
                optional_present += 1;
            } else {
                warnings.push(format!(
                    "dataset={} optional_gap=geometry.as_built",
                    dataset.id
                ));
            }
            if dataset.sensors.iter().all(|sensor| {
                matches!(
                    &sensor.calibration,
                    Availability::Available(calibration)
                        if calibration.valid_through.is_some()
                )
            }) {
                optional_present += 1;
            } else {
                warnings.push(format!(
                    "dataset={} optional_gap=sensors.calibration.valid_through",
                    dataset.id
                ));
            }
            let status = if !validation_ok {
                "ERROR"
            } else if warnings.len() > warning_start {
                "WARN"
            } else {
                "OK"
            };
            rows.push(CorpusAuditRow {
                dataset_id: dataset.id.clone(),
                mandatory_present: if validation_ok { 15 } else { 0 },
                mandatory_total: 15,
                optional_present,
                optional_total,
                partition: dataset.partition,
                evidence_level: dataset.evidence_level,
                physical_cap: dataset.physical_claim_cap(),
                status,
            });
        }
        CorpusAuditReport {
            rows,
            warnings,
            errors,
        }
    }
}

const fn color_name(color: ColorRank) -> &'static str {
    match color {
        ColorRank::Estimated => "estimated",
        ColorRank::Validated => "validated",
        ColorRank::Verified => "verified",
    }
}

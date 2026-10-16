//! Stable, bounded experiment-campaign IR compiled from an exact V&V
//! ContextOfUse.
//!
//! The compiler admits campaign structure only. It binds why runs and
//! measurements exist, fixes calibration/validation/blind partitions before
//! observations arrive, and produces a canonical identity. It does not
//! authorize laboratory execution, establish scientific adequacy, authenticate
//! artifacts, or prove instrumentation safety.

use core::fmt;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, Field, FieldSpec, IdentityReceipt,
    NeverCancel, ProblemSemanticId, WireType,
};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::vv::{ContextOfUse, QoiId, UnitId, VV_ARTIFACT_FAMILY};

use crate::IR_VERSION;

/// Version of the first stable experiment-campaign schema.
pub const EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum rows in any one campaign collection.
pub const MAX_CAMPAIGN_ITEMS_V1: usize = 4_096;
/// Maximum bytes in a stable campaign identifier.
pub const MAX_CAMPAIGN_KEY_BYTES_V1: usize = 128;
/// Maximum bytes in a descriptive campaign field.
pub const MAX_CAMPAIGN_TEXT_BYTES_V1: usize = 4_096;
/// Maximum canonical bytes in one admitted campaign.
pub const MAX_CAMPAIGN_CANONICAL_BYTES_V1: usize = 16 * 1_024 * 1_024;

const CAMPAIGN_MAGIC_V1: &[u8; 8] = b"FSCAMP01";
const CAMPAIGN_WIRE_DOMAIN_V1: &str = "org.frankensim.fs-ir.experiment-campaign-wire.v1";
const CAMPAIGN_IDENTITY_LIMITS_V1: CanonicalLimits = CanonicalLimits::new(
    32 * 1_024 * 1_024,
    MAX_CAMPAIGN_CANONICAL_BYTES_V1 as u64,
    4,
    8,
    64,
);

macro_rules! campaign_id {
    ($name:ident, $role:literal) => {
        #[doc = concat!("Stable ", $role, " identifier.")]
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(Box<str>);

        impl $name {
            /// Construct a bounded canonical identifier.
            ///
            /// Segments use the exact grammar [a-z][a-z0-9-]* and may be
            /// separated by forward slashes.
            pub fn try_new(value: impl Into<String>) -> Result<Self, CampaignError> {
                let value = value.into();
                validate_key($role, &value)?;
                Ok(Self(value.into_boxed_str()))
            }

            /// Exact canonical identifier text.
            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(&self.0)
            }
        }
    };
}

campaign_id!(CampaignClaimId, "claim");
campaign_id!(EvidenceGapId, "evidence-gap");
campaign_id!(SpecimenId, "specimen");
campaign_id!(AssemblyId, "assembly");
campaign_id!(FactorId, "factor");
campaign_id!(ResourceId, "resource");
campaign_id!(MeasurementChannelId, "measurement-channel");
campaign_id!(CampaignRunId, "run");
campaign_id!(AnalysisId, "analysis");
campaign_id!(CampaignRuleId, "campaign-rule");

/// Calibration, validation, and blind-holdout data are non-interchangeable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CampaignPartition {
    /// Data that may influence fitted parameters.
    Calibration,
    /// Held-apart data used for declared validation analyses.
    Validation,
    /// Capability-separated validation data whose assignment remains blinded.
    BlindHoldout,
}

impl CampaignPartition {
    const fn tag(self) -> u8 {
        match self {
            Self::Calibration => 1,
            Self::Validation => 2,
            Self::BlindHoldout => 3,
        }
    }

    fn from_tag(tag: u8, offset: usize) -> Result<Self, CampaignError> {
        match tag {
            1 => Ok(Self::Calibration),
            2 => Ok(Self::Validation),
            3 => Ok(Self::BlindHoldout),
            _ => Err(CampaignError::MalformedCanonical {
                offset,
                detail: format!("unknown campaign partition tag {tag}"),
            }),
        }
    }
}

/// Semantic role of a directed claim dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceUse {
    /// The prerequisite may influence calibration of the dependent claim.
    CalibrationInput,
    /// The prerequisite is consumed only during validation of the dependent.
    ValidationInput,
}

impl EvidenceUse {
    const fn tag(self) -> u8 {
        match self {
            Self::CalibrationInput => 1,
            Self::ValidationInput => 2,
        }
    }

    fn from_tag(tag: u8, offset: usize) -> Result<Self, CampaignError> {
        match tag {
            1 => Ok(Self::CalibrationInput),
            2 => Ok(Self::ValidationInput),
            _ => Err(CampaignError::MalformedCanonical {
                offset,
                detail: format!("unknown evidence-use tag {tag}"),
            }),
        }
    }
}

/// A campaign termination rule is either an orderly stop or a safety abort.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CampaignRuleKind {
    /// Stop after already-started work reaches its declared drain boundary.
    StopAfterDrain,
    /// Abort toward the named minimum-risk state.
    AbortToSafeState,
}

impl CampaignRuleKind {
    const fn tag(self) -> u8 {
        match self {
            Self::StopAfterDrain => 1,
            Self::AbortToSafeState => 2,
        }
    }

    fn from_tag(tag: u8, offset: usize) -> Result<Self, CampaignError> {
        match tag {
            1 => Ok(Self::StopAfterDrain),
            2 => Ok(Self::AbortToSafeState),
            _ => Err(CampaignError::MalformedCanonical {
                offset,
                detail: format!("unknown campaign-rule tag {tag}"),
            }),
        }
    }
}

/// One missing evidence item that motivates the campaign.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceGap {
    /// Stable gap identity.
    pub id: EvidenceGapId,
    /// Exact ContextOfUse QoI affected by the gap.
    pub qoi: QoiId,
    /// Expected evidence family or artifact contract.
    pub expected_evidence: String,
    /// Human-auditable explanation of what remains unknown.
    pub description: String,
}

/// One decision claim served by the campaign.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignClaim {
    /// Stable claim identity.
    pub id: CampaignClaimId,
    /// Exact ContextOfUse QoIs this claim can change.
    pub qois: Vec<QoiId>,
    /// Preregistered scientific hypothesis.
    pub hypothesis: String,
    /// Decision consequence if the claim passes or fails.
    pub decision_consequence: String,
    /// Evidence gaps the planned work is intended to reduce.
    pub evidence_gaps: Vec<EvidenceGap>,
}

/// One physical or synthetic specimen occurrence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpecimenSpec {
    /// Stable specimen identity; collection position is never identity.
    pub id: SpecimenId,
    /// Declared specimen kind or configuration.
    pub kind: String,
}

/// One assembly and its exact specimen membership.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblySpec {
    /// Stable assembly identity.
    pub id: AssemblyId,
    /// Nonempty exact specimen membership.
    pub specimens: Vec<SpecimenId>,
}

/// One experimental factor and its finite preregistered levels.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorSpec {
    /// Stable factor identity.
    pub id: FactorId,
    /// Explicit unit for every level.
    pub unit: UnitId,
    /// Finite, duplicate-free levels. Signed zero is normalized.
    pub levels: Vec<f64>,
}

/// One factor assignment in a run.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorSetting {
    /// Referenced factor.
    pub factor: FactorId,
    /// Exact finite level, required to appear in the factor declaration.
    pub level: f64,
}

/// One bounded execution resource class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignResource {
    /// Stable resource identity.
    pub id: ResourceId,
    /// Explicit capabilities required from the eventual session.
    pub capabilities: Vec<String>,
    /// Maximum simultaneous runs the resource declaration permits.
    pub max_concurrent_runs: u32,
}

/// One measurement channel tied to a claim, QoI, and decision consequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MeasurementChannel {
    /// Stable channel identity.
    pub id: MeasurementChannelId,
    /// Claim this channel informs.
    pub claim: CampaignClaimId,
    /// Exact ContextOfUse QoI measured by the channel.
    pub qoi: QoiId,
    /// Explicit unit, required to equal the ContextOfUse QoI unit.
    pub unit: UnitId,
    /// Human-auditable role in the decision.
    pub decision_consequence: String,
}

/// One preregistered campaign run.
#[derive(Debug, Clone, PartialEq)]
pub struct CampaignRun {
    /// Stable run identity.
    pub id: CampaignRunId,
    /// Exact specimen occurrence.
    pub specimen: SpecimenId,
    /// Assembly that must contain the specimen.
    pub assembly: AssemblyId,
    /// Calibration/validation/blind partition fixed before data exists.
    pub partition: CampaignPartition,
    /// Nonempty claim set served by this run.
    pub claims: Vec<CampaignClaimId>,
    /// Nonempty measurement-channel set produced by this run.
    pub channels: Vec<MeasurementChannelId>,
    /// Explicit factor assignments.
    pub factors: Vec<FactorSetting>,
    /// Deterministic slot generated under the campaign randomization plan.
    pub randomization_slot: u64,
    /// Whether operators and analysis paths are blinded to assignment.
    pub blinded: bool,
    /// Required execution resource.
    pub resource: ResourceId,
    /// Declared worst-case wall duration for budget admission.
    pub wall_time_ms: u64,
    /// Declared peak memory for budget admission.
    pub memory_bytes: u64,
}

/// One immutable preregistered analysis commitment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreregisteredAnalysis {
    /// Stable analysis identity.
    pub id: AnalysisId,
    /// Claim adjudicated by this analysis.
    pub claim: CampaignClaimId,
    /// Exact QoIs consumed by the analysis.
    pub qois: Vec<QoiId>,
    /// Partition from which the analysis may read.
    pub partition: CampaignPartition,
    /// Domain-separated commitment to the complete analysis plan.
    pub preregistration_hash: ContentHash,
    /// Bounded method/estimator identifier.
    pub method: String,
}

/// One directed dependency between claim evidence paths.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ClaimDependency {
    /// Claim whose evidence is consumed.
    pub prerequisite: CampaignClaimId,
    /// Claim whose adjudication consumes it.
    pub dependent: CampaignClaimId,
    /// Calibration or validation use.
    pub use_kind: EvidenceUse,
}

/// One preregistered stop or abort rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignRule {
    /// Stable rule identity.
    pub id: CampaignRuleId,
    /// Orderly stop or safety abort semantics.
    pub kind: CampaignRuleKind,
    /// Machine-oriented predicate description.
    pub predicate: String,
    /// Required drain/finalize or minimum-risk-state action.
    pub action: String,
}

/// Explicit campaign-wide resource budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CampaignBudget {
    /// Maximum admitted run count.
    pub max_runs: u64,
    /// Maximum admitted specimen count.
    pub max_specimens: u64,
    /// Maximum conservative sum of declared run wall times.
    pub max_wall_time_ms: u64,
    /// Maximum peak memory of any one run.
    pub max_memory_bytes: u64,
}

/// Reproducible randomization and blind-assignment commitment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RandomizationPlan {
    /// Counter-stream seed fixed before data exists.
    pub seed: u64,
    /// Exact algorithm/version token.
    pub algorithm: String,
    /// Commitment to blind assignments; all-zero is never an identity.
    pub blind_assignment_commitment: ContentHash,
}

/// Optional predecessor binding for an explicit future schema migration.
///
/// The anchor preserves source bytes and declared intent as separate hashes.
/// It does not assert semantic equivalence between schema versions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignHistoryAnchor {
    /// Strictly older source schema version.
    pub source_schema_version: u32,
    /// Exact predecessor canonical-byte coordinate.
    pub source_canonical_hash: ContentHash,
    /// Separate commitment to the predecessor's declared campaign intent.
    pub source_intent_hash: ContentHash,
}

/// Caller-authored campaign draft. Admission canonicalizes collection order.
#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentCampaignDraft {
    /// Optional explicit predecessor binding.
    pub history: Option<CampaignHistoryAnchor>,
    /// Campaign-wide budget.
    pub budget: CampaignBudget,
    /// Randomization/blinding commitment.
    pub randomization: RandomizationPlan,
    /// Target claims.
    pub claims: Vec<CampaignClaim>,
    /// Directed evidence dependencies.
    pub dependencies: Vec<ClaimDependency>,
    /// Specimen occurrences.
    pub specimens: Vec<SpecimenSpec>,
    /// Assemblies and exact membership.
    pub assemblies: Vec<AssemblySpec>,
    /// Experimental factors.
    pub factors: Vec<FactorSpec>,
    /// Execution resource declarations.
    pub resources: Vec<CampaignResource>,
    /// Measurement channels.
    pub channels: Vec<MeasurementChannel>,
    /// Preregistered runs.
    pub runs: Vec<CampaignRun>,
    /// Preregistered analyses.
    pub analyses: Vec<PreregisteredAnalysis>,
    /// Stop and abort rules.
    pub rules: Vec<CampaignRule>,
}

/// Non-fatal admission diagnostics retained with the compiled artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CampaignWarning {
    /// A declared channel is not used by any run.
    UnusedMeasurement {
        /// Exact orphaned channel.
        channel: MeasurementChannelId,
    },
}

/// Typed refusal from campaign construction or canonical decoding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CampaignError {
    /// A stable key violated the bounded canonical grammar.
    InvalidKey {
        /// Semantic role of the key.
        role: &'static str,
        /// Exact reason.
        detail: String,
    },
    /// Descriptive text was empty, unbounded, or contained controls.
    InvalidText {
        /// Stable field path.
        field: &'static str,
        /// Exact reason.
        detail: String,
    },
    /// A collection or transport exceeded a public cap.
    ResourceLimit {
        /// Stable field path.
        field: &'static str,
        /// Submitted amount.
        actual: usize,
        /// Public maximum.
        max: usize,
    },
    /// A mandatory collection or reference set was empty.
    Empty {
        /// Stable field path.
        field: &'static str,
    },
    /// A semantic identity appeared more than once.
    Duplicate {
        /// Stable field path.
        field: &'static str,
        /// Duplicate identity.
        key: String,
    },
    /// A referenced identity does not exist in the campaign or ContextOfUse.
    UnknownReference {
        /// Stable field path.
        field: &'static str,
        /// Missing identity.
        key: String,
    },
    /// One specimen crossed exclusive evidence partitions.
    PartitionLeakage {
        /// Leaking specimen.
        specimen: SpecimenId,
        /// First admitted partition.
        first: CampaignPartition,
        /// Conflicting partition.
        second: CampaignPartition,
    },
    /// The evidence-dependency graph is cyclic.
    DependencyCycle {
        /// Deterministic set of claims retained in the cycle remainder.
        claims: Vec<CampaignClaimId>,
    },
    /// A declared run or collection exceeds an explicit campaign budget.
    BudgetConflict {
        /// Budget axis.
        field: &'static str,
        /// Required amount.
        required: u64,
        /// Declared limit.
        limit: u64,
    },
    /// A numeric factor or per-run resource declaration is invalid.
    InvalidValue {
        /// Stable field path.
        field: &'static str,
        /// Exact reason.
        detail: String,
    },
    /// An all-zero digest was supplied where an identity commitment is required.
    ZeroHash {
        /// Stable field path.
        field: &'static str,
    },
    /// Canonical V&V ContextOfUse transport failed.
    ContextCodec {
        /// Underlying bounded transport diagnostic.
        detail: String,
    },
    /// Canonical campaign bytes were malformed.
    MalformedCanonical {
        /// Byte offset of the first refusal.
        offset: usize,
        /// Exact reason.
        detail: String,
    },
    /// Bytes decoded semantically but were not the unique canonical fixed point.
    NonCanonical,
    /// Strong-identity construction refused.
    Identity {
        /// Underlying canonical-frame diagnostic.
        detail: String,
    },
}

impl CampaignError {
    /// Stable machine diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidKey { .. } => "CampaignInvalidKey",
            Self::InvalidText { .. } => "CampaignInvalidText",
            Self::ResourceLimit { .. } => "CampaignResourceLimit",
            Self::Empty { .. } => "CampaignEmpty",
            Self::Duplicate { .. } => "CampaignDuplicate",
            Self::UnknownReference { .. } => "CampaignUnknownReference",
            Self::PartitionLeakage { .. } => "CampaignPartitionLeakage",
            Self::DependencyCycle { .. } => "CampaignDependencyCycle",
            Self::BudgetConflict { .. } => "CampaignBudgetConflict",
            Self::InvalidValue { .. } => "CampaignInvalidValue",
            Self::ZeroHash { .. } => "CampaignZeroHash",
            Self::ContextCodec { .. } => "CampaignContextCodec",
            Self::MalformedCanonical { .. } => "CampaignMalformedCanonical",
            Self::NonCanonical => "CampaignNonCanonical",
            Self::Identity { .. } => "CampaignIdentity",
        }
    }

    /// Concise repair direction suitable for structured agent diagnostics.
    #[must_use]
    pub const fn hint(&self) -> &'static str {
        match self {
            Self::InvalidKey { .. } => "use a bounded lowercase segmented identity",
            Self::InvalidText { .. } => "supply bounded nonblank control-free text",
            Self::ResourceLimit { .. } => "split the campaign or reduce the declared set",
            Self::Empty { .. } => "declare the mandatory campaign rows explicitly",
            Self::Duplicate { .. } => "assign unique stable identities",
            Self::UnknownReference { .. } => {
                "bind the exact declared campaign or ContextOfUse identity"
            }
            Self::PartitionLeakage { .. } => {
                "use disjoint specimen occurrences for exclusive partitions"
            }
            Self::DependencyCycle { .. } => {
                "remove circular calibration or validation evidence flow"
            }
            Self::BudgetConflict { .. } => "increase the explicit budget or reduce planned work",
            Self::InvalidValue { .. } => "supply a finite in-domain declared value",
            Self::ZeroHash { .. } => "supply a nonzero content-addressed commitment",
            Self::ContextCodec { .. } => "repair the source ContextOfUse canonical artifact",
            Self::MalformedCanonical { .. } | Self::NonCanonical => {
                "re-emit through the current canonical campaign encoder"
            }
            Self::Identity { .. } => {
                "reduce the campaign to the public canonical identity envelope"
            }
        }
    }
}

impl fmt::Display for CampaignError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey { role, detail } => write!(formatter, "invalid {role} key: {detail}"),
            Self::InvalidText { field, detail } => {
                write!(formatter, "invalid {field} text: {detail}")
            }
            Self::ResourceLimit { field, actual, max } => {
                write!(
                    formatter,
                    "{field} has {actual} items/bytes; maximum is {max}"
                )
            }
            Self::Empty { field } => write!(formatter, "{field} must not be empty"),
            Self::Duplicate { field, key } => {
                write!(formatter, "{field} contains duplicate identity {key}")
            }
            Self::UnknownReference { field, key } => {
                write!(formatter, "{field} references unknown identity {key}")
            }
            Self::PartitionLeakage {
                specimen,
                first,
                second,
            } => write!(
                formatter,
                "specimen {specimen} crosses exclusive partitions {first:?} and {second:?}"
            ),
            Self::DependencyCycle { claims } => {
                write!(
                    formatter,
                    "claim evidence dependencies are cyclic: {claims:?}"
                )
            }
            Self::BudgetConflict {
                field,
                required,
                limit,
            } => write!(
                formatter,
                "{field} requires {required}, exceeding declared limit {limit}"
            ),
            Self::InvalidValue { field, detail } => {
                write!(formatter, "invalid {field}: {detail}")
            }
            Self::ZeroHash { field } => write!(formatter, "{field} must not be all zero"),
            Self::ContextCodec { detail } => {
                write!(
                    formatter,
                    "ContextOfUse canonical transport refused: {detail}"
                )
            }
            Self::MalformedCanonical { offset, detail } => {
                write!(
                    formatter,
                    "malformed campaign bytes at offset {offset}: {detail}"
                )
            }
            Self::NonCanonical => formatter
                .write_str("campaign bytes decode but are not the unique canonical fixed point"),
            Self::Identity { detail } => {
                write!(
                    formatter,
                    "campaign identity construction refused: {detail}"
                )
            }
        }
    }
}

impl std::error::Error for CampaignError {}

/// Canonical identity schema for one admitted experiment campaign.
pub enum ExperimentCampaignIdentitySchemaV1 {}

impl CanonicalSchema for ExperimentCampaignIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.experiment-campaign.v1";
    const NAME: &'static str = "experiment-campaign-ir";
    const VERSION: u32 = EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "exact ContextOfUse and canonical preregistered experiment-campaign intent";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("campaign-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("context-of-use-content", WireType::Bytes),
        FieldSpec::required("canonical-campaign", WireType::Bytes),
    ];
}

/// Strong semantic identity of one admitted campaign.
pub type ExperimentCampaignIdV1 = ProblemSemanticId<ExperimentCampaignIdentitySchemaV1>;

/// Immutable admitted experiment-campaign IR.
#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentCampaignIr {
    context: ContextOfUse,
    history: Option<CampaignHistoryAnchor>,
    budget: CampaignBudget,
    randomization: RandomizationPlan,
    claims: Vec<CampaignClaim>,
    dependencies: Vec<ClaimDependency>,
    specimens: Vec<SpecimenSpec>,
    assemblies: Vec<AssemblySpec>,
    factors: Vec<FactorSpec>,
    resources: Vec<CampaignResource>,
    channels: Vec<MeasurementChannel>,
    runs: Vec<CampaignRun>,
    analyses: Vec<PreregisteredAnalysis>,
    rules: Vec<CampaignRule>,
    warnings: Vec<CampaignWarning>,
    canonical_bytes: Vec<u8>,
    wire_hash: ContentHash,
    identity: IdentityReceipt<ExperimentCampaignIdV1>,
}

impl ExperimentCampaignIr {
    /// Compile one exact ContextOfUse and caller draft into canonical campaign IR.
    ///
    /// Admission canonicalizes caller collection order, fixes all identities
    /// and partitions, rejects dangling/circular/leaking structure, and
    /// publishes no identity on refusal.
    #[allow(clippy::too_many_lines)]
    pub fn compile(
        context: ContextOfUse,
        mut draft: ExperimentCampaignDraft,
    ) -> Result<Self, CampaignError> {
        validate_history(draft.history.as_ref())?;
        validate_budget(draft.budget)?;
        validate_randomization(&draft.randomization)?;
        validate_collection_count("campaign.claims", draft.claims.len(), true)?;
        validate_collection_count("campaign.dependencies", draft.dependencies.len(), false)?;
        validate_collection_count("campaign.specimens", draft.specimens.len(), true)?;
        validate_collection_count("campaign.assemblies", draft.assemblies.len(), true)?;
        validate_collection_count("campaign.factors", draft.factors.len(), false)?;
        validate_collection_count("campaign.resources", draft.resources.len(), true)?;
        validate_collection_count("campaign.channels", draft.channels.len(), true)?;
        validate_collection_count("campaign.runs", draft.runs.len(), true)?;
        validate_collection_count("campaign.analyses", draft.analyses.len(), true)?;
        validate_collection_count("campaign.rules", draft.rules.len(), true)?;

        canonicalize_top_level(&mut draft)?;

        let claim_index = validate_claims(&context, &mut draft.claims)?;
        validate_dependencies(&draft.dependencies, &claim_index)?;
        let specimen_index = validate_specimens(&draft.specimens)?;
        let assembly_owner = validate_assemblies(&mut draft.assemblies, &specimen_index)?;
        let factor_index = validate_factors(&mut draft.factors)?;
        let resource_index = validate_resources(&mut draft.resources)?;
        let channel_index = validate_channels(&context, &claim_index, &draft.channels)?;
        let run_summary = validate_runs(
            &mut draft.runs,
            &claim_index,
            &specimen_index,
            &assembly_owner,
            &factor_index,
            &resource_index,
            &channel_index,
            draft.budget,
        )?;
        validate_analyses(
            &mut draft.analyses,
            &context,
            &claim_index,
            &run_summary.claim_partitions,
        )?;
        validate_rules(&draft.rules)?;

        let warnings = draft
            .channels
            .iter()
            .filter(|channel| !run_summary.used_channels.contains(&channel.id))
            .map(|channel| CampaignWarning::UnusedMeasurement {
                channel: channel.id.clone(),
            })
            .collect();

        let canonical_bytes = encode_campaign(&context, &draft)?;
        let context_bytes =
            context
                .canonical_bytes()
                .map_err(|error| CampaignError::ContextCodec {
                    detail: error.to_string(),
                })?;
        let context_hash = hash_domain(VV_ARTIFACT_FAMILY, &context_bytes);
        let identity = CanonicalEncoder::<ExperimentCampaignIdV1, _>::new(
            CAMPAIGN_IDENTITY_LIMITS_V1,
            NeverCancel,
        )
        .map_err(identity_error)?
        .u64(
            Field::new(0, "campaign-schema-version"),
            u64::from(EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1),
        )
        .map_err(identity_error)?
        .u64(
            Field::new(1, "frankenscript-ir-version"),
            u64::from(IR_VERSION),
        )
        .map_err(identity_error)?
        .bytes(
            Field::new(2, "context-of-use-content"),
            context_hash.as_bytes(),
        )
        .map_err(identity_error)?
        .bytes(Field::new(3, "canonical-campaign"), &canonical_bytes)
        .map_err(identity_error)?
        .finish()
        .map_err(identity_error)?;
        let wire_hash = hash_domain(CAMPAIGN_WIRE_DOMAIN_V1, &canonical_bytes);

        Ok(Self {
            context,
            history: draft.history,
            budget: draft.budget,
            randomization: draft.randomization,
            claims: draft.claims,
            dependencies: draft.dependencies,
            specimens: draft.specimens,
            assemblies: draft.assemblies,
            factors: draft.factors,
            resources: draft.resources,
            channels: draft.channels,
            runs: draft.runs,
            analyses: draft.analyses,
            rules: draft.rules,
            warnings,
            canonical_bytes,
            wire_hash,
            identity,
        })
    }

    /// Decode and readmit the unique current canonical representation.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, CampaignError> {
        if bytes.len() > MAX_CAMPAIGN_CANONICAL_BYTES_V1 {
            return Err(CampaignError::ResourceLimit {
                field: "campaign.canonical-bytes",
                actual: bytes.len(),
                max: MAX_CAMPAIGN_CANONICAL_BYTES_V1,
            });
        }
        let mut decoder = Decoder::new(bytes);
        decoder.expect_magic()?;
        let schema_offset = decoder.offset();
        let schema_version = decoder.u32("campaign.schema-version")?;
        if schema_version != EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1 {
            return Err(CampaignError::MalformedCanonical {
                offset: schema_offset,
                detail: format!(
                    "unsupported campaign schema version {schema_version}; current is {}",
                    EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1
                ),
            });
        }
        let ir_offset = decoder.offset();
        let ir_version = decoder.u32("campaign.ir-version")?;
        if ir_version != IR_VERSION {
            return Err(CampaignError::MalformedCanonical {
                offset: ir_offset,
                detail: format!(
                    "unsupported FrankenScript IR version {ir_version}; current is {IR_VERSION}"
                ),
            });
        }
        let context_offset = decoder.offset();
        let context_bytes = decoder.bytes("campaign.context")?;
        let context = ContextOfUse::from_canonical_bytes(context_bytes).map_err(|error| {
            CampaignError::MalformedCanonical {
                offset: context_offset,
                detail: format!("invalid ContextOfUse artifact: {error}"),
            }
        })?;
        let history = decode_history(&mut decoder)?;
        let budget = CampaignBudget {
            max_runs: decoder.u64("campaign.budget.max-runs")?,
            max_specimens: decoder.u64("campaign.budget.max-specimens")?,
            max_wall_time_ms: decoder.u64("campaign.budget.max-wall-time-ms")?,
            max_memory_bytes: decoder.u64("campaign.budget.max-memory-bytes")?,
        };
        let randomization = RandomizationPlan {
            seed: decoder.u64("campaign.randomization.seed")?,
            algorithm: decoder.string("campaign.randomization.algorithm")?,
            blind_assignment_commitment: decoder
                .hash("campaign.randomization.blind-assignment-commitment")?,
        };
        let claims = decode_claims(&mut decoder)?;
        let dependencies = decode_dependencies(&mut decoder)?;
        let specimens = decode_specimens(&mut decoder)?;
        let assemblies = decode_assemblies(&mut decoder)?;
        let factors = decode_factors(&mut decoder)?;
        let resources = decode_resources(&mut decoder)?;
        let channels = decode_channels(&mut decoder)?;
        let runs = decode_runs(&mut decoder)?;
        let analyses = decode_analyses(&mut decoder)?;
        let rules = decode_rules(&mut decoder)?;
        decoder.finish()?;

        let admitted = Self::compile(
            context,
            ExperimentCampaignDraft {
                history,
                budget,
                randomization,
                claims,
                dependencies,
                specimens,
                assemblies,
                factors,
                resources,
                channels,
                runs,
                analyses,
                rules,
            },
        )?;
        if admitted.canonical_bytes.as_slice() != bytes {
            return Err(CampaignError::NonCanonical);
        }
        Ok(admitted)
    }

    /// Strong semantic identity.
    #[must_use]
    pub const fn id(&self) -> ExperimentCampaignIdV1 {
        self.identity.id()
    }

    /// Complete canonical strong-identity receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<ExperimentCampaignIdV1> {
        self.identity
    }

    /// Domain-separated content coordinate for the exact wire bytes.
    #[must_use]
    pub const fn wire_hash(&self) -> ContentHash {
        self.wire_hash
    }

    /// Unique current canonical representation.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Exact source ContextOfUse.
    #[must_use]
    pub const fn context(&self) -> &ContextOfUse {
        &self.context
    }

    /// Optional predecessor intent binding.
    #[must_use]
    pub const fn history(&self) -> Option<&CampaignHistoryAnchor> {
        self.history.as_ref()
    }

    /// Explicit admitted budget.
    #[must_use]
    pub const fn budget(&self) -> CampaignBudget {
        self.budget
    }

    /// Exact randomization/blinding commitment.
    #[must_use]
    pub const fn randomization(&self) -> &RandomizationPlan {
        &self.randomization
    }

    /// Canonically ordered claims.
    #[must_use]
    pub fn claims(&self) -> &[CampaignClaim] {
        &self.claims
    }

    /// Canonically ordered directed evidence dependencies.
    #[must_use]
    pub fn dependencies(&self) -> &[ClaimDependency] {
        &self.dependencies
    }

    /// Canonically ordered specimens.
    #[must_use]
    pub fn specimens(&self) -> &[SpecimenSpec] {
        &self.specimens
    }

    /// Canonically ordered assemblies.
    #[must_use]
    pub fn assemblies(&self) -> &[AssemblySpec] {
        &self.assemblies
    }

    /// Canonically ordered factors.
    #[must_use]
    pub fn factors(&self) -> &[FactorSpec] {
        &self.factors
    }

    /// Canonically ordered resource declarations.
    #[must_use]
    pub fn resources(&self) -> &[CampaignResource] {
        &self.resources
    }

    /// Canonically ordered measurement channels.
    #[must_use]
    pub fn channels(&self) -> &[MeasurementChannel] {
        &self.channels
    }

    /// Canonically ordered runs.
    #[must_use]
    pub fn runs(&self) -> &[CampaignRun] {
        &self.runs
    }

    /// Canonically ordered preregistered analyses.
    #[must_use]
    pub fn analyses(&self) -> &[PreregisteredAnalysis] {
        &self.analyses
    }

    /// Canonically ordered stop/abort rules.
    #[must_use]
    pub fn rules(&self) -> &[CampaignRule] {
        &self.rules
    }

    /// Derived non-fatal diagnostics.
    #[must_use]
    pub fn warnings(&self) -> &[CampaignWarning] {
        &self.warnings
    }
}

fn identity_error(error: impl fmt::Display) -> CampaignError {
    CampaignError::Identity {
        detail: error.to_string(),
    }
}

fn validate_key(role: &'static str, value: &str) -> Result<(), CampaignError> {
    if value.is_empty() {
        return Err(CampaignError::InvalidKey {
            role,
            detail: "identity is empty".to_owned(),
        });
    }
    if value.len() > MAX_CAMPAIGN_KEY_BYTES_V1 {
        return Err(CampaignError::InvalidKey {
            role,
            detail: format!(
                "identity has {} bytes; maximum is {}",
                value.len(),
                MAX_CAMPAIGN_KEY_BYTES_V1
            ),
        });
    }
    for (segment_index, segment) in value.split('/').enumerate() {
        let bytes = segment.as_bytes();
        if bytes.is_empty() || !bytes[0].is_ascii_lowercase() {
            return Err(CampaignError::InvalidKey {
                role,
                detail: format!(
                    "segment {segment_index} must start with an ASCII lowercase letter"
                ),
            });
        }
        if let Some((offset, byte)) = bytes.iter().copied().enumerate().find(|(_, byte)| {
            !(byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'-')
        }) {
            return Err(CampaignError::InvalidKey {
                role,
                detail: format!(
                    "segment {segment_index} byte {offset} ({byte:#04x}) violates [a-z0-9-]"
                ),
            });
        }
    }
    Ok(())
}

fn validate_text(field: &'static str, value: &str) -> Result<(), CampaignError> {
    if value.trim().is_empty() {
        return Err(CampaignError::InvalidText {
            field,
            detail: "value is blank".to_owned(),
        });
    }
    if value.len() > MAX_CAMPAIGN_TEXT_BYTES_V1 {
        return Err(CampaignError::InvalidText {
            field,
            detail: format!(
                "value has {} bytes; maximum is {}",
                value.len(),
                MAX_CAMPAIGN_TEXT_BYTES_V1
            ),
        });
    }
    if let Some((offset, character)) = value.char_indices().find(|(_, value)| value.is_control()) {
        return Err(CampaignError::InvalidText {
            field,
            detail: format!(
                "control character U+{:04X} at byte {offset}",
                character as u32
            ),
        });
    }
    Ok(())
}

fn validate_token(field: &'static str, value: &str) -> Result<(), CampaignError> {
    validate_text(field, value)?;
    if let Some((offset, byte)) = value.bytes().enumerate().find(|(_, byte)| {
        !(byte.is_ascii_alphanumeric() || matches!(*byte, b'-' | b'_' | b'.' | b':' | b'/'))
    }) {
        return Err(CampaignError::InvalidText {
            field,
            detail: format!("token byte {offset} ({byte:#04x}) is not canonical"),
        });
    }
    Ok(())
}

fn validate_collection_count(
    field: &'static str,
    count: usize,
    nonempty: bool,
) -> Result<(), CampaignError> {
    if nonempty && count == 0 {
        return Err(CampaignError::Empty { field });
    }
    if count > MAX_CAMPAIGN_ITEMS_V1 {
        return Err(CampaignError::ResourceLimit {
            field,
            actual: count,
            max: MAX_CAMPAIGN_ITEMS_V1,
        });
    }
    Ok(())
}

fn reject_adjacent_duplicate<T, F>(
    rows: &[T],
    field: &'static str,
    key: F,
) -> Result<(), CampaignError>
where
    F: Fn(&T) -> &str,
{
    if let Some(pair) = rows.windows(2).find(|pair| key(&pair[0]) == key(&pair[1])) {
        return Err(CampaignError::Duplicate {
            field,
            key: key(&pair[0]).to_owned(),
        });
    }
    Ok(())
}

fn canonicalize_refs<T>(values: &mut Vec<T>, field: &'static str) -> Result<(), CampaignError>
where
    T: Ord + fmt::Display,
{
    validate_collection_count(field, values.len(), true)?;
    values.sort();
    if let Some(pair) = values.windows(2).find(|pair| pair[0] == pair[1]) {
        return Err(CampaignError::Duplicate {
            field,
            key: pair[0].to_string(),
        });
    }
    Ok(())
}

fn canonicalize_top_level(draft: &mut ExperimentCampaignDraft) -> Result<(), CampaignError> {
    draft.claims.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.claims, "campaign.claims", |row| row.id.as_str())?;
    draft.dependencies.sort();
    if let Some(pair) = draft
        .dependencies
        .windows(2)
        .find(|pair| pair[0] == pair[1])
    {
        return Err(CampaignError::Duplicate {
            field: "campaign.dependencies",
            key: format!(
                "{}>{}:{:?}",
                pair[0].prerequisite, pair[0].dependent, pair[0].use_kind
            ),
        });
    }
    draft
        .specimens
        .sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.specimens, "campaign.specimens", |row| {
        row.id.as_str()
    })?;
    draft
        .assemblies
        .sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.assemblies, "campaign.assemblies", |row| {
        row.id.as_str()
    })?;
    draft.factors.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.factors, "campaign.factors", |row| row.id.as_str())?;
    draft
        .resources
        .sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.resources, "campaign.resources", |row| {
        row.id.as_str()
    })?;
    draft.channels.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.channels, "campaign.channels", |row| row.id.as_str())?;
    draft.runs.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.runs, "campaign.runs", |row| row.id.as_str())?;
    draft.analyses.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.analyses, "campaign.analyses", |row| row.id.as_str())?;
    draft.rules.sort_by(|left, right| left.id.cmp(&right.id));
    reject_adjacent_duplicate(&draft.rules, "campaign.rules", |row| row.id.as_str())?;
    Ok(())
}

fn validate_history(history: Option<&CampaignHistoryAnchor>) -> Result<(), CampaignError> {
    let Some(history) = history else {
        return Ok(());
    };
    if history.source_schema_version >= EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1 {
        return Err(CampaignError::InvalidValue {
            field: "campaign.history.source-schema-version",
            detail: format!(
                "source version {} must be older than {}",
                history.source_schema_version, EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1
            ),
        });
    }
    require_nonzero_hash(
        "campaign.history.source-canonical-hash",
        history.source_canonical_hash,
    )?;
    require_nonzero_hash(
        "campaign.history.source-intent-hash",
        history.source_intent_hash,
    )
}

fn validate_budget(budget: CampaignBudget) -> Result<(), CampaignError> {
    for (field, value) in [
        ("campaign.budget.max-runs", budget.max_runs),
        ("campaign.budget.max-specimens", budget.max_specimens),
        ("campaign.budget.max-wall-time-ms", budget.max_wall_time_ms),
        ("campaign.budget.max-memory-bytes", budget.max_memory_bytes),
    ] {
        if value == 0 {
            return Err(CampaignError::InvalidValue {
                field,
                detail: "budget must be positive".to_owned(),
            });
        }
    }
    Ok(())
}

fn validate_randomization(plan: &RandomizationPlan) -> Result<(), CampaignError> {
    validate_token("campaign.randomization.algorithm", &plan.algorithm)?;
    require_nonzero_hash(
        "campaign.randomization.blind-assignment-commitment",
        plan.blind_assignment_commitment,
    )
}

fn validate_claims(
    context: &ContextOfUse,
    claims: &mut [CampaignClaim],
) -> Result<BTreeMap<CampaignClaimId, BTreeSet<QoiId>>, CampaignError> {
    let mut index = BTreeMap::new();
    let mut gap_ids = BTreeSet::new();
    for claim in claims.iter_mut() {
        validate_text("campaign.claim.hypothesis", &claim.hypothesis)?;
        validate_text(
            "campaign.claim.decision-consequence",
            &claim.decision_consequence,
        )?;
        canonicalize_refs(&mut claim.qois, "campaign.claim.qois")?;
        for qoi in &claim.qois {
            if !context.qois().contains_key(qoi) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.claim.qois",
                    key: qoi.to_string(),
                });
            }
        }
        validate_collection_count(
            "campaign.claim.evidence-gaps",
            claim.evidence_gaps.len(),
            false,
        )?;
        claim
            .evidence_gaps
            .sort_by(|left, right| left.id.cmp(&right.id));
        reject_adjacent_duplicate(
            &claim.evidence_gaps,
            "campaign.claim.evidence-gaps",
            |row| row.id.as_str(),
        )?;
        for gap in &claim.evidence_gaps {
            if !gap_ids.insert(gap.id.clone()) {
                return Err(CampaignError::Duplicate {
                    field: "campaign.evidence-gap-id",
                    key: gap.id.to_string(),
                });
            }
            if claim.qois.binary_search(&gap.qoi).is_err() {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.claim.evidence-gap.qoi",
                    key: gap.qoi.to_string(),
                });
            }
            validate_token(
                "campaign.claim.evidence-gap.expected-evidence",
                &gap.expected_evidence,
            )?;
            validate_text("campaign.claim.evidence-gap.description", &gap.description)?;
        }
        index.insert(claim.id.clone(), claim.qois.iter().cloned().collect());
    }
    Ok(index)
}

fn validate_dependencies(
    dependencies: &[ClaimDependency],
    claims: &BTreeMap<CampaignClaimId, BTreeSet<QoiId>>,
) -> Result<(), CampaignError> {
    let mut indegree: BTreeMap<CampaignClaimId, usize> =
        claims.keys().cloned().map(|id| (id, 0)).collect();
    let mut outgoing: BTreeMap<CampaignClaimId, Vec<CampaignClaimId>> = BTreeMap::new();
    let mut endpoint_pairs = BTreeSet::new();
    for dependency in dependencies {
        for (field, claim) in [
            ("campaign.dependency.prerequisite", &dependency.prerequisite),
            ("campaign.dependency.dependent", &dependency.dependent),
        ] {
            if !claims.contains_key(claim) {
                return Err(CampaignError::UnknownReference {
                    field,
                    key: claim.to_string(),
                });
            }
        }
        if !endpoint_pairs.insert((
            dependency.prerequisite.clone(),
            dependency.dependent.clone(),
        )) {
            return Err(CampaignError::Duplicate {
                field: "campaign.dependency.endpoints",
                key: format!("{}>{}", dependency.prerequisite, dependency.dependent),
            });
        }
        outgoing
            .entry(dependency.prerequisite.clone())
            .or_default()
            .push(dependency.dependent.clone());
        *indegree
            .get_mut(&dependency.dependent)
            .expect("dependency membership checked") += 1;
    }
    for values in outgoing.values_mut() {
        values.sort();
    }
    let mut ready: BTreeSet<CampaignClaimId> = indegree
        .iter()
        .filter(|(_, degree)| **degree == 0)
        .map(|(claim, _)| claim.clone())
        .collect();
    let mut visited = 0usize;
    while let Some(claim) = ready.pop_first() {
        visited += 1;
        if let Some(dependents) = outgoing.get(&claim) {
            for dependent in dependents {
                let degree = indegree
                    .get_mut(dependent)
                    .expect("dependency membership checked");
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(dependent.clone());
                }
            }
        }
    }
    if visited != claims.len() {
        let claims = indegree
            .into_iter()
            .filter(|(_, degree)| *degree != 0)
            .map(|(claim, _)| claim)
            .collect();
        return Err(CampaignError::DependencyCycle { claims });
    }
    Ok(())
}

fn validate_specimens(specimens: &[SpecimenSpec]) -> Result<BTreeSet<SpecimenId>, CampaignError> {
    for specimen in specimens {
        validate_token("campaign.specimen.kind", &specimen.kind)?;
    }
    Ok(specimens.iter().map(|row| row.id.clone()).collect())
}

fn validate_assemblies(
    assemblies: &mut [AssemblySpec],
    specimens: &BTreeSet<SpecimenId>,
) -> Result<BTreeMap<SpecimenId, AssemblyId>, CampaignError> {
    let mut owner = BTreeMap::new();
    for assembly in assemblies {
        canonicalize_refs(&mut assembly.specimens, "campaign.assembly.specimens")?;
        for specimen in &assembly.specimens {
            if !specimens.contains(specimen) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.assembly.specimens",
                    key: specimen.to_string(),
                });
            }
            if let Some(previous) = owner.insert(specimen.clone(), assembly.id.clone()) {
                return Err(CampaignError::Duplicate {
                    field: "campaign.specimen.assembly-membership",
                    key: format!("{specimen} in {previous} and {}", assembly.id),
                });
            }
        }
    }
    for specimen in specimens {
        if !owner.contains_key(specimen) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.specimen.assembly-membership",
                key: specimen.to_string(),
            });
        }
    }
    Ok(owner)
}

fn normalize_finite(field: &'static str, value: &mut f64) -> Result<(), CampaignError> {
    if !value.is_finite() {
        return Err(CampaignError::InvalidValue {
            field,
            detail: "value must be finite".to_owned(),
        });
    }
    if value.to_bits() == (-0.0_f64).to_bits() {
        *value = 0.0;
    }
    Ok(())
}

fn validate_factors(
    factors: &mut [FactorSpec],
) -> Result<BTreeMap<FactorId, BTreeSet<u64>>, CampaignError> {
    let mut index = BTreeMap::new();
    for factor in factors {
        validate_collection_count("campaign.factor.levels", factor.levels.len(), true)?;
        for level in &mut factor.levels {
            normalize_finite("campaign.factor.level", level)?;
        }
        factor.levels.sort_by(f64::total_cmp);
        if factor
            .levels
            .windows(2)
            .any(|pair| pair[0].to_bits() == pair[1].to_bits())
        {
            return Err(CampaignError::Duplicate {
                field: "campaign.factor.levels",
                key: factor.id.to_string(),
            });
        }
        index.insert(
            factor.id.clone(),
            factor.levels.iter().map(|value| value.to_bits()).collect(),
        );
    }
    Ok(index)
}

fn validate_resources(
    resources: &mut [CampaignResource],
) -> Result<BTreeSet<ResourceId>, CampaignError> {
    for resource in resources.iter_mut() {
        if resource.max_concurrent_runs == 0 {
            return Err(CampaignError::InvalidValue {
                field: "campaign.resource.max-concurrent-runs",
                detail: format!("resource {} declares zero concurrency", resource.id),
            });
        }
        validate_collection_count(
            "campaign.resource.capabilities",
            resource.capabilities.len(),
            true,
        )?;
        resource.capabilities.sort();
        if let Some(pair) = resource
            .capabilities
            .windows(2)
            .find(|pair| pair[0] == pair[1])
        {
            return Err(CampaignError::Duplicate {
                field: "campaign.resource.capabilities",
                key: pair[0].clone(),
            });
        }
        for capability in &resource.capabilities {
            validate_token("campaign.resource.capability", capability)?;
        }
    }
    Ok(resources.iter().map(|row| row.id.clone()).collect())
}

fn validate_channels(
    context: &ContextOfUse,
    claims: &BTreeMap<CampaignClaimId, BTreeSet<QoiId>>,
    channels: &[MeasurementChannel],
) -> Result<BTreeMap<MeasurementChannelId, (CampaignClaimId, QoiId)>, CampaignError> {
    let mut index = BTreeMap::new();
    for channel in channels {
        let Some(claim_qois) = claims.get(&channel.claim) else {
            return Err(CampaignError::UnknownReference {
                field: "campaign.channel.claim",
                key: channel.claim.to_string(),
            });
        };
        let Some(qoi) = context.qois().get(&channel.qoi) else {
            return Err(CampaignError::UnknownReference {
                field: "campaign.channel.qoi",
                key: channel.qoi.to_string(),
            });
        };
        if !claim_qois.contains(&channel.qoi) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.channel.claim-qoi",
                key: format!("{}:{}", channel.claim, channel.qoi),
            });
        }
        if qoi.unit() != &channel.unit {
            return Err(CampaignError::UnknownReference {
                field: "campaign.channel.unit",
                key: format!(
                    "{} expects {}, got {}",
                    channel.qoi,
                    qoi.unit(),
                    channel.unit
                ),
            });
        }
        validate_text(
            "campaign.channel.decision-consequence",
            &channel.decision_consequence,
        )?;
        index.insert(
            channel.id.clone(),
            (channel.claim.clone(), channel.qoi.clone()),
        );
    }
    Ok(index)
}

struct RunSummary {
    used_channels: BTreeSet<MeasurementChannelId>,
    claim_partitions: BTreeSet<(CampaignClaimId, CampaignPartition)>,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn validate_runs(
    runs: &mut [CampaignRun],
    claims: &BTreeMap<CampaignClaimId, BTreeSet<QoiId>>,
    specimens: &BTreeSet<SpecimenId>,
    assembly_owner: &BTreeMap<SpecimenId, AssemblyId>,
    factors: &BTreeMap<FactorId, BTreeSet<u64>>,
    resources: &BTreeSet<ResourceId>,
    channels: &BTreeMap<MeasurementChannelId, (CampaignClaimId, QoiId)>,
    budget: CampaignBudget,
) -> Result<RunSummary, CampaignError> {
    budget_check(
        "campaign.budget.max-runs",
        runs.len() as u64,
        budget.max_runs,
    )?;
    budget_check(
        "campaign.budget.max-specimens",
        specimens.len() as u64,
        budget.max_specimens,
    )?;
    let mut used_channels = BTreeSet::new();
    let mut claim_partitions = BTreeSet::new();
    let mut specimen_partitions = BTreeMap::new();
    let mut randomization_slots = BTreeMap::new();
    let mut wall_time_ms = 0u64;
    let mut peak_memory_bytes = 0u64;

    for run in runs {
        canonicalize_refs(&mut run.claims, "campaign.run.claims")?;
        canonicalize_refs(&mut run.channels, "campaign.run.channels")?;
        if !specimens.contains(&run.specimen) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.run.specimen",
                key: run.specimen.to_string(),
            });
        }
        let expected_assembly = assembly_owner
            .get(&run.specimen)
            .expect("complete specimen ownership checked");
        if expected_assembly != &run.assembly {
            return Err(CampaignError::UnknownReference {
                field: "campaign.run.assembly",
                key: format!(
                    "{} belongs to {}, not {}",
                    run.specimen, expected_assembly, run.assembly
                ),
            });
        }
        if !resources.contains(&run.resource) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.run.resource",
                key: run.resource.to_string(),
            });
        }
        if run.wall_time_ms == 0 || run.memory_bytes == 0 {
            return Err(CampaignError::InvalidValue {
                field: "campaign.run.resource-budget",
                detail: format!("run {} requires positive wall and memory bounds", run.id),
            });
        }
        if run.partition == CampaignPartition::BlindHoldout && !run.blinded {
            return Err(CampaignError::InvalidValue {
                field: "campaign.run.blinded",
                detail: format!("blind-holdout run {} must be blinded", run.id),
            });
        }
        if let Some(previous) = specimen_partitions.insert(run.specimen.clone(), run.partition)
            && previous != run.partition
        {
            return Err(CampaignError::PartitionLeakage {
                specimen: run.specimen.clone(),
                first: previous,
                second: run.partition,
            });
        }
        if let Some(previous) = randomization_slots.insert(run.randomization_slot, run.id.clone()) {
            return Err(CampaignError::Duplicate {
                field: "campaign.run.randomization-slot",
                key: format!(
                    "{} used by {previous} and {}",
                    run.randomization_slot, run.id
                ),
            });
        }
        for claim in &run.claims {
            if !claims.contains_key(claim) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.run.claims",
                    key: claim.to_string(),
                });
            }
            claim_partitions.insert((claim.clone(), run.partition));
        }
        for channel in &run.channels {
            let Some((channel_claim, _)) = channels.get(channel) else {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.run.channels",
                    key: channel.to_string(),
                });
            };
            if run.claims.binary_search(channel_claim).is_err() {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.run.channel-claim",
                    key: format!("{channel}:{channel_claim}"),
                });
            }
            used_channels.insert(channel.clone());
        }
        validate_collection_count("campaign.run.factors", run.factors.len(), false)?;
        for setting in &mut run.factors {
            normalize_finite("campaign.run.factor-level", &mut setting.level)?;
        }
        run.factors
            .sort_by(|left, right| left.factor.cmp(&right.factor));
        if let Some(pair) = run
            .factors
            .windows(2)
            .find(|pair| pair[0].factor == pair[1].factor)
        {
            return Err(CampaignError::Duplicate {
                field: "campaign.run.factors",
                key: pair[0].factor.to_string(),
            });
        }
        for setting in &run.factors {
            let Some(levels) = factors.get(&setting.factor) else {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.run.factor",
                    key: setting.factor.to_string(),
                });
            };
            if !levels.contains(&setting.level.to_bits()) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.run.factor-level",
                    key: format!("{}:{:016x}", setting.factor, setting.level.to_bits()),
                });
            }
        }
        wall_time_ms =
            wall_time_ms
                .checked_add(run.wall_time_ms)
                .ok_or(CampaignError::BudgetConflict {
                    field: "campaign.budget.max-wall-time-ms",
                    required: u64::MAX,
                    limit: budget.max_wall_time_ms,
                })?;
        peak_memory_bytes = peak_memory_bytes.max(run.memory_bytes);
    }
    budget_check(
        "campaign.budget.max-wall-time-ms",
        wall_time_ms,
        budget.max_wall_time_ms,
    )?;
    budget_check(
        "campaign.budget.max-memory-bytes",
        peak_memory_bytes,
        budget.max_memory_bytes,
    )?;
    Ok(RunSummary {
        used_channels,
        claim_partitions,
    })
}

fn validate_analyses(
    analyses: &mut [PreregisteredAnalysis],
    context: &ContextOfUse,
    claims: &BTreeMap<CampaignClaimId, BTreeSet<QoiId>>,
    claim_partitions: &BTreeSet<(CampaignClaimId, CampaignPartition)>,
) -> Result<(), CampaignError> {
    let mut claims_with_analysis = BTreeSet::new();
    for analysis in analyses {
        let Some(claim_qois) = claims.get(&analysis.claim) else {
            return Err(CampaignError::UnknownReference {
                field: "campaign.analysis.claim",
                key: analysis.claim.to_string(),
            });
        };
        canonicalize_refs(&mut analysis.qois, "campaign.analysis.qois")?;
        for qoi in &analysis.qois {
            if !context.qois().contains_key(qoi) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.analysis.qois",
                    key: qoi.to_string(),
                });
            }
            if !claim_qois.contains(qoi) {
                return Err(CampaignError::UnknownReference {
                    field: "campaign.analysis.claim-qoi",
                    key: format!("{}:{qoi}", analysis.claim),
                });
            }
        }
        if !claim_partitions.contains(&(analysis.claim.clone(), analysis.partition)) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.analysis.partition",
                key: format!("{}:{:?}", analysis.claim, analysis.partition),
            });
        }
        require_nonzero_hash(
            "campaign.analysis.preregistration-hash",
            analysis.preregistration_hash,
        )?;
        validate_token("campaign.analysis.method", &analysis.method)?;
        claims_with_analysis.insert(analysis.claim.clone());
    }
    for claim in claims.keys() {
        if !claims_with_analysis.contains(claim) {
            return Err(CampaignError::UnknownReference {
                field: "campaign.claim.preregistered-analysis",
                key: claim.to_string(),
            });
        }
    }
    Ok(())
}

fn validate_rules(rules: &[CampaignRule]) -> Result<(), CampaignError> {
    let mut has_stop = false;
    let mut has_abort = false;
    for rule in rules {
        validate_text("campaign.rule.predicate", &rule.predicate)?;
        validate_text("campaign.rule.action", &rule.action)?;
        match rule.kind {
            CampaignRuleKind::StopAfterDrain => has_stop = true,
            CampaignRuleKind::AbortToSafeState => has_abort = true,
        }
    }
    if !has_stop {
        return Err(CampaignError::Empty {
            field: "campaign.rules.stop-after-drain",
        });
    }
    if !has_abort {
        return Err(CampaignError::Empty {
            field: "campaign.rules.abort-to-safe-state",
        });
    }
    Ok(())
}

fn budget_check(field: &'static str, required: u64, limit: u64) -> Result<(), CampaignError> {
    if required > limit {
        Err(CampaignError::BudgetConflict {
            field,
            required,
            limit,
        })
    } else {
        Ok(())
    }
}

fn require_nonzero_hash(field: &'static str, hash: ContentHash) -> Result<(), CampaignError> {
    if hash.as_bytes() == &[0; 32] {
        Err(CampaignError::ZeroHash { field })
    } else {
        Ok(())
    }
}

struct Encoder {
    bytes: Vec<u8>,
}

impl Encoder {
    fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    fn raw(&mut self, field: &'static str, value: &[u8]) -> Result<(), CampaignError> {
        let actual =
            self.bytes
                .len()
                .checked_add(value.len())
                .ok_or(CampaignError::ResourceLimit {
                    field,
                    actual: usize::MAX,
                    max: MAX_CAMPAIGN_CANONICAL_BYTES_V1,
                })?;
        if actual > MAX_CAMPAIGN_CANONICAL_BYTES_V1 {
            return Err(CampaignError::ResourceLimit {
                field,
                actual,
                max: MAX_CAMPAIGN_CANONICAL_BYTES_V1,
            });
        }
        self.bytes
            .try_reserve(value.len())
            .map_err(|_| CampaignError::ResourceLimit {
                field,
                actual,
                max: MAX_CAMPAIGN_CANONICAL_BYTES_V1,
            })?;
        self.bytes.extend_from_slice(value);
        Ok(())
    }

    fn u8(&mut self, field: &'static str, value: u8) -> Result<(), CampaignError> {
        self.raw(field, &[value])
    }

    fn bool(&mut self, field: &'static str, value: bool) -> Result<(), CampaignError> {
        self.u8(field, u8::from(value))
    }

    fn u32(&mut self, field: &'static str, value: u32) -> Result<(), CampaignError> {
        self.raw(field, &value.to_le_bytes())
    }

    fn u64(&mut self, field: &'static str, value: u64) -> Result<(), CampaignError> {
        self.raw(field, &value.to_le_bytes())
    }

    fn hash(&mut self, field: &'static str, value: ContentHash) -> Result<(), CampaignError> {
        self.raw(field, value.as_bytes())
    }

    fn bytes(&mut self, field: &'static str, value: &[u8]) -> Result<(), CampaignError> {
        let length = u32::try_from(value.len()).map_err(|_| CampaignError::ResourceLimit {
            field,
            actual: value.len(),
            max: u32::MAX as usize,
        })?;
        self.u32(field, length)?;
        self.raw(field, value)
    }

    fn string(&mut self, field: &'static str, value: &str) -> Result<(), CampaignError> {
        self.bytes(field, value.as_bytes())
    }

    fn count(&mut self, field: &'static str, count: usize) -> Result<(), CampaignError> {
        let count = u32::try_from(count).map_err(|_| CampaignError::ResourceLimit {
            field,
            actual: count,
            max: u32::MAX as usize,
        })?;
        self.u32(field, count)
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

fn encode_campaign(
    context: &ContextOfUse,
    draft: &ExperimentCampaignDraft,
) -> Result<Vec<u8>, CampaignError> {
    let mut encoder = Encoder::new();
    encoder.raw("campaign.magic", CAMPAIGN_MAGIC_V1)?;
    encoder.u32(
        "campaign.schema-version",
        EXPERIMENT_CAMPAIGN_SCHEMA_VERSION_V1,
    )?;
    encoder.u32("campaign.ir-version", IR_VERSION)?;
    let context_bytes = context
        .canonical_bytes()
        .map_err(|error| CampaignError::ContextCodec {
            detail: error.to_string(),
        })?;
    encoder.bytes("campaign.context", &context_bytes)?;
    encode_history(&mut encoder, draft.history.as_ref())?;
    encoder.u64("campaign.budget.max-runs", draft.budget.max_runs)?;
    encoder.u64("campaign.budget.max-specimens", draft.budget.max_specimens)?;
    encoder.u64(
        "campaign.budget.max-wall-time-ms",
        draft.budget.max_wall_time_ms,
    )?;
    encoder.u64(
        "campaign.budget.max-memory-bytes",
        draft.budget.max_memory_bytes,
    )?;
    encoder.u64("campaign.randomization.seed", draft.randomization.seed)?;
    encoder.string(
        "campaign.randomization.algorithm",
        &draft.randomization.algorithm,
    )?;
    encoder.hash(
        "campaign.randomization.blind-assignment-commitment",
        draft.randomization.blind_assignment_commitment,
    )?;
    encode_claims(&mut encoder, &draft.claims)?;
    encode_dependencies(&mut encoder, &draft.dependencies)?;
    encode_specimens(&mut encoder, &draft.specimens)?;
    encode_assemblies(&mut encoder, &draft.assemblies)?;
    encode_factors(&mut encoder, &draft.factors)?;
    encode_resources(&mut encoder, &draft.resources)?;
    encode_channels(&mut encoder, &draft.channels)?;
    encode_runs(&mut encoder, &draft.runs)?;
    encode_analyses(&mut encoder, &draft.analyses)?;
    encode_rules(&mut encoder, &draft.rules)?;
    Ok(encoder.finish())
}

fn encode_history(
    encoder: &mut Encoder,
    history: Option<&CampaignHistoryAnchor>,
) -> Result<(), CampaignError> {
    encoder.bool("campaign.history.present", history.is_some())?;
    if let Some(history) = history {
        encoder.u32(
            "campaign.history.source-schema-version",
            history.source_schema_version,
        )?;
        encoder.hash(
            "campaign.history.source-canonical-hash",
            history.source_canonical_hash,
        )?;
        encoder.hash(
            "campaign.history.source-intent-hash",
            history.source_intent_hash,
        )?;
    }
    Ok(())
}

fn encode_claims(encoder: &mut Encoder, claims: &[CampaignClaim]) -> Result<(), CampaignError> {
    encoder.count("campaign.claims", claims.len())?;
    for claim in claims {
        encoder.string("campaign.claim.id", claim.id.as_str())?;
        encoder.string("campaign.claim.hypothesis", &claim.hypothesis)?;
        encoder.string(
            "campaign.claim.decision-consequence",
            &claim.decision_consequence,
        )?;
        encoder.count("campaign.claim.qois", claim.qois.len())?;
        for qoi in &claim.qois {
            encoder.string("campaign.claim.qoi", qoi.as_str())?;
        }
        encoder.count("campaign.claim.evidence-gaps", claim.evidence_gaps.len())?;
        for gap in &claim.evidence_gaps {
            encoder.string("campaign.claim.evidence-gap.id", gap.id.as_str())?;
            encoder.string("campaign.claim.evidence-gap.qoi", gap.qoi.as_str())?;
            encoder.string(
                "campaign.claim.evidence-gap.expected-evidence",
                &gap.expected_evidence,
            )?;
            encoder.string("campaign.claim.evidence-gap.description", &gap.description)?;
        }
    }
    Ok(())
}

fn encode_dependencies(
    encoder: &mut Encoder,
    dependencies: &[ClaimDependency],
) -> Result<(), CampaignError> {
    encoder.count("campaign.dependencies", dependencies.len())?;
    for dependency in dependencies {
        encoder.string(
            "campaign.dependency.prerequisite",
            dependency.prerequisite.as_str(),
        )?;
        encoder.string(
            "campaign.dependency.dependent",
            dependency.dependent.as_str(),
        )?;
        encoder.u8("campaign.dependency.use-kind", dependency.use_kind.tag())?;
    }
    Ok(())
}

fn encode_specimens(
    encoder: &mut Encoder,
    specimens: &[SpecimenSpec],
) -> Result<(), CampaignError> {
    encoder.count("campaign.specimens", specimens.len())?;
    for specimen in specimens {
        encoder.string("campaign.specimen.id", specimen.id.as_str())?;
        encoder.string("campaign.specimen.kind", &specimen.kind)?;
    }
    Ok(())
}

fn encode_assemblies(
    encoder: &mut Encoder,
    assemblies: &[AssemblySpec],
) -> Result<(), CampaignError> {
    encoder.count("campaign.assemblies", assemblies.len())?;
    for assembly in assemblies {
        encoder.string("campaign.assembly.id", assembly.id.as_str())?;
        encoder.count("campaign.assembly.specimens", assembly.specimens.len())?;
        for specimen in &assembly.specimens {
            encoder.string("campaign.assembly.specimen", specimen.as_str())?;
        }
    }
    Ok(())
}

fn encode_factors(encoder: &mut Encoder, factors: &[FactorSpec]) -> Result<(), CampaignError> {
    encoder.count("campaign.factors", factors.len())?;
    for factor in factors {
        encoder.string("campaign.factor.id", factor.id.as_str())?;
        encoder.string("campaign.factor.unit", factor.unit.as_str())?;
        encoder.count("campaign.factor.levels", factor.levels.len())?;
        for level in &factor.levels {
            encoder.u64("campaign.factor.level", level.to_bits())?;
        }
    }
    Ok(())
}

fn encode_resources(
    encoder: &mut Encoder,
    resources: &[CampaignResource],
) -> Result<(), CampaignError> {
    encoder.count("campaign.resources", resources.len())?;
    for resource in resources {
        encoder.string("campaign.resource.id", resource.id.as_str())?;
        encoder.u32(
            "campaign.resource.max-concurrent-runs",
            resource.max_concurrent_runs,
        )?;
        encoder.count(
            "campaign.resource.capabilities",
            resource.capabilities.len(),
        )?;
        for capability in &resource.capabilities {
            encoder.string("campaign.resource.capability", capability)?;
        }
    }
    Ok(())
}

fn encode_channels(
    encoder: &mut Encoder,
    channels: &[MeasurementChannel],
) -> Result<(), CampaignError> {
    encoder.count("campaign.channels", channels.len())?;
    for channel in channels {
        encoder.string("campaign.channel.id", channel.id.as_str())?;
        encoder.string("campaign.channel.claim", channel.claim.as_str())?;
        encoder.string("campaign.channel.qoi", channel.qoi.as_str())?;
        encoder.string("campaign.channel.unit", channel.unit.as_str())?;
        encoder.string(
            "campaign.channel.decision-consequence",
            &channel.decision_consequence,
        )?;
    }
    Ok(())
}

fn encode_runs(encoder: &mut Encoder, runs: &[CampaignRun]) -> Result<(), CampaignError> {
    encoder.count("campaign.runs", runs.len())?;
    for run in runs {
        encoder.string("campaign.run.id", run.id.as_str())?;
        encoder.string("campaign.run.specimen", run.specimen.as_str())?;
        encoder.string("campaign.run.assembly", run.assembly.as_str())?;
        encoder.u8("campaign.run.partition", run.partition.tag())?;
        encoder.count("campaign.run.claims", run.claims.len())?;
        for claim in &run.claims {
            encoder.string("campaign.run.claim", claim.as_str())?;
        }
        encoder.count("campaign.run.channels", run.channels.len())?;
        for channel in &run.channels {
            encoder.string("campaign.run.channel", channel.as_str())?;
        }
        encoder.count("campaign.run.factors", run.factors.len())?;
        for setting in &run.factors {
            encoder.string("campaign.run.factor", setting.factor.as_str())?;
            encoder.u64("campaign.run.factor-level", setting.level.to_bits())?;
        }
        encoder.u64("campaign.run.randomization-slot", run.randomization_slot)?;
        encoder.bool("campaign.run.blinded", run.blinded)?;
        encoder.string("campaign.run.resource", run.resource.as_str())?;
        encoder.u64("campaign.run.wall-time-ms", run.wall_time_ms)?;
        encoder.u64("campaign.run.memory-bytes", run.memory_bytes)?;
    }
    Ok(())
}

fn encode_analyses(
    encoder: &mut Encoder,
    analyses: &[PreregisteredAnalysis],
) -> Result<(), CampaignError> {
    encoder.count("campaign.analyses", analyses.len())?;
    for analysis in analyses {
        encoder.string("campaign.analysis.id", analysis.id.as_str())?;
        encoder.string("campaign.analysis.claim", analysis.claim.as_str())?;
        encoder.u8("campaign.analysis.partition", analysis.partition.tag())?;
        encoder.hash(
            "campaign.analysis.preregistration-hash",
            analysis.preregistration_hash,
        )?;
        encoder.string("campaign.analysis.method", &analysis.method)?;
        encoder.count("campaign.analysis.qois", analysis.qois.len())?;
        for qoi in &analysis.qois {
            encoder.string("campaign.analysis.qoi", qoi.as_str())?;
        }
    }
    Ok(())
}

fn encode_rules(encoder: &mut Encoder, rules: &[CampaignRule]) -> Result<(), CampaignError> {
    encoder.count("campaign.rules", rules.len())?;
    for rule in rules {
        encoder.string("campaign.rule.id", rule.id.as_str())?;
        encoder.u8("campaign.rule.kind", rule.kind.tag())?;
        encoder.string("campaign.rule.predicate", &rule.predicate)?;
        encoder.string("campaign.rule.action", &rule.action)?;
    }
    Ok(())
}

struct Decoder<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Decoder<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    const fn offset(&self) -> usize {
        self.offset
    }

    fn take(&mut self, length: usize, field: &'static str) -> Result<&'a [u8], CampaignError> {
        let end =
            self.offset
                .checked_add(length)
                .ok_or_else(|| CampaignError::MalformedCanonical {
                    offset: self.offset,
                    detail: format!("{field} length overflow"),
                })?;
        let Some(value) = self.bytes.get(self.offset..end) else {
            return Err(CampaignError::MalformedCanonical {
                offset: self.offset,
                detail: format!("{field} is truncated"),
            });
        };
        self.offset = end;
        Ok(value)
    }

    fn expect_magic(&mut self) -> Result<(), CampaignError> {
        let offset = self.offset;
        if self.take(CAMPAIGN_MAGIC_V1.len(), "campaign.magic")? != CAMPAIGN_MAGIC_V1 {
            return Err(CampaignError::MalformedCanonical {
                offset,
                detail: "campaign magic mismatch".to_owned(),
            });
        }
        Ok(())
    }

    fn u8(&mut self, field: &'static str) -> Result<u8, CampaignError> {
        Ok(self.take(1, field)?[0])
    }

    fn bool(&mut self, field: &'static str) -> Result<bool, CampaignError> {
        let offset = self.offset;
        match self.u8(field)? {
            0 => Ok(false),
            1 => Ok(true),
            value => Err(CampaignError::MalformedCanonical {
                offset,
                detail: format!("{field} has non-Boolean tag {value}"),
            }),
        }
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, CampaignError> {
        let bytes: [u8; 4] = self.take(4, field)?.try_into().expect("fixed-size slice");
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self, field: &'static str) -> Result<u64, CampaignError> {
        let bytes: [u8; 8] = self.take(8, field)?.try_into().expect("fixed-size slice");
        Ok(u64::from_le_bytes(bytes))
    }

    fn count(&mut self, field: &'static str) -> Result<usize, CampaignError> {
        let count = self.u32(field)? as usize;
        if count > MAX_CAMPAIGN_ITEMS_V1 {
            return Err(CampaignError::ResourceLimit {
                field,
                actual: count,
                max: MAX_CAMPAIGN_ITEMS_V1,
            });
        }
        Ok(count)
    }

    fn bytes(&mut self, field: &'static str) -> Result<&'a [u8], CampaignError> {
        let length = self.u32(field)? as usize;
        if length > MAX_CAMPAIGN_CANONICAL_BYTES_V1 {
            return Err(CampaignError::ResourceLimit {
                field,
                actual: length,
                max: MAX_CAMPAIGN_CANONICAL_BYTES_V1,
            });
        }
        self.take(length, field)
    }

    fn string(&mut self, field: &'static str) -> Result<String, CampaignError> {
        let offset = self.offset;
        let bytes = self.bytes(field)?;
        let value =
            core::str::from_utf8(bytes).map_err(|error| CampaignError::MalformedCanonical {
                offset,
                detail: format!("{field} is not UTF-8: {error}"),
            })?;
        Ok(value.to_owned())
    }

    fn hash(&mut self, field: &'static str) -> Result<ContentHash, CampaignError> {
        let offset = self.offset;
        ContentHash::from_slice(self.take(32, field)?).ok_or_else(|| {
            CampaignError::MalformedCanonical {
                offset,
                detail: format!("{field} is not a 32-byte hash"),
            }
        })
    }

    fn finish(self) -> Result<(), CampaignError> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err(CampaignError::MalformedCanonical {
                offset: self.offset,
                detail: format!(
                    "{} trailing bytes after the campaign",
                    self.bytes.len() - self.offset
                ),
            })
        }
    }
}

fn decode_history(
    decoder: &mut Decoder<'_>,
) -> Result<Option<CampaignHistoryAnchor>, CampaignError> {
    if !decoder.bool("campaign.history.present")? {
        return Ok(None);
    }
    Ok(Some(CampaignHistoryAnchor {
        source_schema_version: decoder.u32("campaign.history.source-schema-version")?,
        source_canonical_hash: decoder.hash("campaign.history.source-canonical-hash")?,
        source_intent_hash: decoder.hash("campaign.history.source-intent-hash")?,
    }))
}

fn decode_id<T>(
    decoder: &mut Decoder<'_>,
    field: &'static str,
    constructor: impl FnOnce(String) -> Result<T, CampaignError>,
) -> Result<T, CampaignError> {
    constructor(decoder.string(field)?)
}

fn decode_qoi(decoder: &mut Decoder<'_>, field: &'static str) -> Result<QoiId, CampaignError> {
    let value = decoder.string(field)?;
    QoiId::try_new(value).map_err(|error| CampaignError::MalformedCanonical {
        offset: decoder.offset(),
        detail: format!("{field} is invalid: {error}"),
    })
}

fn decode_unit(decoder: &mut Decoder<'_>, field: &'static str) -> Result<UnitId, CampaignError> {
    let value = decoder.string(field)?;
    UnitId::try_new(value).map_err(|error| CampaignError::MalformedCanonical {
        offset: decoder.offset(),
        detail: format!("{field} is invalid: {error}"),
    })
}

fn decode_claims(decoder: &mut Decoder<'_>) -> Result<Vec<CampaignClaim>, CampaignError> {
    let count = decoder.count("campaign.claims")?;
    let mut claims = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.claim.id", CampaignClaimId::try_new)?;
        let hypothesis = decoder.string("campaign.claim.hypothesis")?;
        let decision_consequence = decoder.string("campaign.claim.decision-consequence")?;
        let qoi_count = decoder.count("campaign.claim.qois")?;
        let mut qois = Vec::with_capacity(qoi_count);
        for _ in 0..qoi_count {
            qois.push(decode_qoi(decoder, "campaign.claim.qoi")?);
        }
        let gap_count = decoder.count("campaign.claim.evidence-gaps")?;
        let mut evidence_gaps = Vec::with_capacity(gap_count);
        for _ in 0..gap_count {
            evidence_gaps.push(EvidenceGap {
                id: decode_id(
                    decoder,
                    "campaign.claim.evidence-gap.id",
                    EvidenceGapId::try_new,
                )?,
                qoi: decode_qoi(decoder, "campaign.claim.evidence-gap.qoi")?,
                expected_evidence: decoder
                    .string("campaign.claim.evidence-gap.expected-evidence")?,
                description: decoder.string("campaign.claim.evidence-gap.description")?,
            });
        }
        claims.push(CampaignClaim {
            id,
            qois,
            hypothesis,
            decision_consequence,
            evidence_gaps,
        });
    }
    Ok(claims)
}

fn decode_dependencies(decoder: &mut Decoder<'_>) -> Result<Vec<ClaimDependency>, CampaignError> {
    let count = decoder.count("campaign.dependencies")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let prerequisite = decode_id(
            decoder,
            "campaign.dependency.prerequisite",
            CampaignClaimId::try_new,
        )?;
        let dependent = decode_id(
            decoder,
            "campaign.dependency.dependent",
            CampaignClaimId::try_new,
        )?;
        let offset = decoder.offset();
        let use_kind = EvidenceUse::from_tag(decoder.u8("campaign.dependency.use-kind")?, offset)?;
        rows.push(ClaimDependency {
            prerequisite,
            dependent,
            use_kind,
        });
    }
    Ok(rows)
}

fn decode_specimens(decoder: &mut Decoder<'_>) -> Result<Vec<SpecimenSpec>, CampaignError> {
    let count = decoder.count("campaign.specimens")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        rows.push(SpecimenSpec {
            id: decode_id(decoder, "campaign.specimen.id", SpecimenId::try_new)?,
            kind: decoder.string("campaign.specimen.kind")?,
        });
    }
    Ok(rows)
}

fn decode_assemblies(decoder: &mut Decoder<'_>) -> Result<Vec<AssemblySpec>, CampaignError> {
    let count = decoder.count("campaign.assemblies")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.assembly.id", AssemblyId::try_new)?;
        let specimen_count = decoder.count("campaign.assembly.specimens")?;
        let mut specimens = Vec::with_capacity(specimen_count);
        for _ in 0..specimen_count {
            specimens.push(decode_id(
                decoder,
                "campaign.assembly.specimen",
                SpecimenId::try_new,
            )?);
        }
        rows.push(AssemblySpec { id, specimens });
    }
    Ok(rows)
}

fn decode_factors(decoder: &mut Decoder<'_>) -> Result<Vec<FactorSpec>, CampaignError> {
    let count = decoder.count("campaign.factors")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.factor.id", FactorId::try_new)?;
        let unit = decode_unit(decoder, "campaign.factor.unit")?;
        let level_count = decoder.count("campaign.factor.levels")?;
        let mut levels = Vec::with_capacity(level_count);
        for _ in 0..level_count {
            levels.push(f64::from_bits(decoder.u64("campaign.factor.level")?));
        }
        rows.push(FactorSpec { id, unit, levels });
    }
    Ok(rows)
}

fn decode_resources(decoder: &mut Decoder<'_>) -> Result<Vec<CampaignResource>, CampaignError> {
    let count = decoder.count("campaign.resources")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.resource.id", ResourceId::try_new)?;
        let max_concurrent_runs = decoder.u32("campaign.resource.max-concurrent-runs")?;
        let capability_count = decoder.count("campaign.resource.capabilities")?;
        let mut capabilities = Vec::with_capacity(capability_count);
        for _ in 0..capability_count {
            capabilities.push(decoder.string("campaign.resource.capability")?);
        }
        rows.push(CampaignResource {
            id,
            capabilities,
            max_concurrent_runs,
        });
    }
    Ok(rows)
}

fn decode_channels(decoder: &mut Decoder<'_>) -> Result<Vec<MeasurementChannel>, CampaignError> {
    let count = decoder.count("campaign.channels")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        rows.push(MeasurementChannel {
            id: decode_id(
                decoder,
                "campaign.channel.id",
                MeasurementChannelId::try_new,
            )?,
            claim: decode_id(decoder, "campaign.channel.claim", CampaignClaimId::try_new)?,
            qoi: decode_qoi(decoder, "campaign.channel.qoi")?,
            unit: decode_unit(decoder, "campaign.channel.unit")?,
            decision_consequence: decoder.string("campaign.channel.decision-consequence")?,
        });
    }
    Ok(rows)
}

fn decode_runs(decoder: &mut Decoder<'_>) -> Result<Vec<CampaignRun>, CampaignError> {
    let count = decoder.count("campaign.runs")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.run.id", CampaignRunId::try_new)?;
        let specimen = decode_id(decoder, "campaign.run.specimen", SpecimenId::try_new)?;
        let assembly = decode_id(decoder, "campaign.run.assembly", AssemblyId::try_new)?;
        let offset = decoder.offset();
        let partition = CampaignPartition::from_tag(decoder.u8("campaign.run.partition")?, offset)?;
        let claim_count = decoder.count("campaign.run.claims")?;
        let mut claims = Vec::with_capacity(claim_count);
        for _ in 0..claim_count {
            claims.push(decode_id(
                decoder,
                "campaign.run.claim",
                CampaignClaimId::try_new,
            )?);
        }
        let channel_count = decoder.count("campaign.run.channels")?;
        let mut channels = Vec::with_capacity(channel_count);
        for _ in 0..channel_count {
            channels.push(decode_id(
                decoder,
                "campaign.run.channel",
                MeasurementChannelId::try_new,
            )?);
        }
        let factor_count = decoder.count("campaign.run.factors")?;
        let mut factors = Vec::with_capacity(factor_count);
        for _ in 0..factor_count {
            factors.push(FactorSetting {
                factor: decode_id(decoder, "campaign.run.factor", FactorId::try_new)?,
                level: f64::from_bits(decoder.u64("campaign.run.factor-level")?),
            });
        }
        rows.push(CampaignRun {
            id,
            specimen,
            assembly,
            partition,
            claims,
            channels,
            factors,
            randomization_slot: decoder.u64("campaign.run.randomization-slot")?,
            blinded: decoder.bool("campaign.run.blinded")?,
            resource: decode_id(decoder, "campaign.run.resource", ResourceId::try_new)?,
            wall_time_ms: decoder.u64("campaign.run.wall-time-ms")?,
            memory_bytes: decoder.u64("campaign.run.memory-bytes")?,
        });
    }
    Ok(rows)
}

fn decode_analyses(decoder: &mut Decoder<'_>) -> Result<Vec<PreregisteredAnalysis>, CampaignError> {
    let count = decoder.count("campaign.analyses")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.analysis.id", AnalysisId::try_new)?;
        let claim = decode_id(decoder, "campaign.analysis.claim", CampaignClaimId::try_new)?;
        let offset = decoder.offset();
        let partition =
            CampaignPartition::from_tag(decoder.u8("campaign.analysis.partition")?, offset)?;
        let preregistration_hash = decoder.hash("campaign.analysis.preregistration-hash")?;
        let method = decoder.string("campaign.analysis.method")?;
        let qoi_count = decoder.count("campaign.analysis.qois")?;
        let mut qois = Vec::with_capacity(qoi_count);
        for _ in 0..qoi_count {
            qois.push(decode_qoi(decoder, "campaign.analysis.qoi")?);
        }
        rows.push(PreregisteredAnalysis {
            id,
            claim,
            qois,
            partition,
            preregistration_hash,
            method,
        });
    }
    Ok(rows)
}

fn decode_rules(decoder: &mut Decoder<'_>) -> Result<Vec<CampaignRule>, CampaignError> {
    let count = decoder.count("campaign.rules")?;
    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let id = decode_id(decoder, "campaign.rule.id", CampaignRuleId::try_new)?;
        let offset = decoder.offset();
        let kind = CampaignRuleKind::from_tag(decoder.u8("campaign.rule.kind")?, offset)?;
        rows.push(CampaignRule {
            id,
            kind,
            predicate: decoder.string("campaign.rule.predicate")?,
            action: decoder.string("campaign.rule.action")?,
        });
    }
    Ok(rows)
}

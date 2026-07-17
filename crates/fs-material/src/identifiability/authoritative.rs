//! Authority-separated, multi-case identifiability schema.
//!
//! This module is the public I10.1 contract.  It deliberately separates four
//! monotone stages that answer different questions:
//!
//! 1. [`IdentifiabilityProblemDocument`] is an unresolved, canonical statement
//!    of the physical/statistical question.  Decoding bytes never grants source
//!    authority.
//! 2. [`AdmittedIdentifiabilityProblem`] resolves every source against concrete
//!    artifacts or an explicit authority disposition. Successful admission
//!    gates [`ProblemId`] while the resolution envelope itself is retained in
//!    the distinct [`SourceAdmissionId`].
//! 3. [`IdentifiabilityExecutionPlan`] binds coordinates, algorithms, seeds,
//!    tolerances, budgets, and build semantics and mints [`ExecutionId`].
//! 4. [`IdentifiabilityAssessment`] binds typed claims to evidence and mints
//!    [`AssessmentId`].
//!
//! Consequently, changing a coordinate system cannot change the physical
//! problem identity, and adding evidence cannot silently rewrite either the
//! problem or the execution that generated it.  Multi-case campaigns are a
//! first-class v3 primitive: complementary protocols, specimens, environments,
//! and observation operators can jointly break symmetries that no single case
//! can resolve.

use super::*;
use fs_evidence::vv::{ExperimentOrigin, MAX_VV_MATRIX_DIMENSION, ObservationLocatorIdentity};

/// Umbrella API generation for the authority-separated I10.1 module. Identity
/// preimages use the four stage-specific versions below, so changing one stage
/// never silently rewrites the other three identities.
pub const IDENTIFIABILITY_AUTHORITY_SCHEMA_VERSION: u32 = 3;
pub const IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION: u32 = 3;
pub const IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION: u32 = 3;
pub const IDENTIFIABILITY_EXECUTION_IDENTITY_VERSION: u32 = 3;
pub const IDENTIFIABILITY_ASSESSMENT_IDENTITY_VERSION: u32 = 3;

const PROBLEM_MAGIC: &[u8] = b"fs-material-identifiability-problem\0";
const SOURCE_ADMISSION_MAGIC: &[u8] = b"fs-material-identifiability-source-admission\0";
const EXECUTION_MAGIC: &[u8] = b"fs-material-identifiability-execution\0";
const ASSESSMENT_MAGIC: &[u8] = b"fs-material-identifiability-assessment\0";
pub const IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.identifiability-problem.v3";
pub const IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.identifiability-source-admission.v3";
pub const IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.identifiability-execution.v3";
pub const IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.identifiability-assessment.v3";

/// Exact fs-evidence V&V artifact digest domain accepted by typed admission.
pub const VV_ARTIFACT_SOURCE_DOMAIN: &str = fs_evidence::vv::VV_ARTIFACT_FAMILY;
/// Exact fs-matdb material-card digest domain accepted by typed admission.
pub const MATERIAL_CARD_SOURCE_DOMAIN: &str = "org.frankensim.fs-matdb.material-card.v1";
/// Exact fs-matdb constitutive-model-card digest domain accepted by typed admission.
pub const CONSTITUTIVE_MODEL_CARD_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-matdb.constitutive-model-card.v1";
/// Contract version shared by the narrow case-physics component artifacts.
pub const CASE_PHYSICS_SOURCE_CONTRACT_VERSION: u32 = 1;
/// Exact digest domains for case-physics artifacts that do not yet have a
/// richer typed artifact schema. Keeping the domains distinct prevents two
/// equal byte strings from aliasing different physical roles.
pub const FRAME_TRANSFORM_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.case-frame-transform.v1";
pub const SPECIMEN_GEOMETRY_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.case-specimen-geometry.v1";
pub const SPECIMEN_PROCESS_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.case-specimen-process.v1";
pub const SPECIMEN_PREPARATION_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.case-specimen-preparation.v1";
pub const LOAD_PATH_SOURCE_DOMAIN: &str = "org.frankensim.fs-material.case-load-path.v1";
pub const ENVIRONMENT_PATH_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.case-environment-path.v1";
pub const TIME_GRID_SOURCE_DOMAIN: &str = "org.frankensim.fs-material.case-time-grid.v1";
pub const INITIAL_STATE_SOURCE_DOMAIN: &str = "org.frankensim.fs-material.case-initial-state.v1";
/// Exact source contract for an external claim that one canonical witness
/// belongs to the full conjunction of opaque admissible-domain constraints.
pub const ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_VERSION: u32 = 2;
pub const ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN: &str =
    "org.frankensim.fs-material.admissible-domain-membership.v2";
/// Exact typed receipt binding one producer artifact to the complete
/// domain/version/kind/key/hash identity of a forward-model source.
pub const FORWARD_MODEL_PRODUCTION_BINDING_VERSION: u32 = 1;
pub const FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN: &str =
    "org.frankensim.fs-material.forward-model-production-binding.v1";
const ADMISSIBLE_DOMAIN_WITNESS_BINDING_DOMAIN: &str =
    "org.frankensim.fs-material.admissible-domain-witness-binding.v2";
pub const BLIND_RELEASE_TRUST_RECEIPT_VERSION: u32 = 1;
pub const BLIND_RELEASE_TRUST_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-evidence.blind-release-trust-receipt.v1";

macro_rules! authority_id {
    ($name:ident, $doc:literal) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ContentHash);

        impl $name {
            /// Inspect the domain-separated digest.
            #[must_use]
            pub const fn digest(self) -> ContentHash {
                self.0
            }
        }
    };
}

authority_id!(
    ProblemId,
    "Admission-gated identity of the unresolved physical-question bytes."
);
authority_id!(
    SourceAdmissionId,
    "Identity of the source-resolution and authority envelope."
);
authority_id!(ExecutionId, "Identity of one numerical execution plan.");
authority_id!(
    AssessmentId,
    "Identity of one typed, evidence-bound assessment."
);

macro_rules! authority_token {
    ($name:ident, $field:literal) => {
        #[doc = concat!("Canonical ", $field, " token.")]
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(String);

        impl $name {
            /// Construct a bounded ASCII machine token.
            pub fn try_new(value: impl Into<String>) -> Result<Self, IdentifiabilityError> {
                let value = value.into();
                validate_token(&value, $field)?;
                Ok(Self(value))
            }

            /// Inspect canonical token text.
            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }
    };
}

authority_token!(CaseId, "study case");
authority_token!(SourceKey, "source key");
authority_token!(ConstraintId, "joint constraint");
authority_token!(InfluenceId, "influence declaration");
authority_token!(ClaimId, "identifiability claim");
authority_token!(GaugeCompositionId, "gauge composition");
authority_token!(GaugeReductionId, "gauge reduction");

/// Composite observation identity.  Local channel names are never treated as
/// globally unique across a campaign.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ObservationKey {
    case: CaseId,
    channel: ObservationChannelId,
}

impl ObservationKey {
    /// Construct a case-qualified observation endpoint.
    #[must_use]
    pub const fn new(case: CaseId, channel: ObservationChannelId) -> Self {
        Self { case, channel }
    }

    /// Owning case.
    #[must_use]
    pub const fn case(&self) -> &CaseId {
        &self.case
    }

    /// Case-local channel.
    #[must_use]
    pub const fn channel(&self) -> &ObservationChannelId {
        &self.channel
    }
}

/// Semantic class of an immutable source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SourceKind {
    ContextOfUse,
    MaterialCard,
    ConstitutiveModelCard,
    ConstitutiveGraph,
    ExperimentArtifact,
    CalibrationSplit,
    ForwardModel,
    Geometry,
    Process,
    Protocol,
    ObservationOperator,
    Metrology,
    Parser,
    Preprocessing,
    Likelihood,
    Prior,
    Constraint,
    GaugeAction,
    GaugeSection,
    Discrepancy,
    Assumption,
    Analyzer,
    DerivativeProvider,
    Build,
    EvidenceReceipt,
    ExternalManifold,
    AlgebraicExtension,
    Stratification,
    DerivedFunctional,
    AdmissibleDomainCertificate,
    DimensionlessErrorMetric,
    Nondimensionalization,
    QuantifierRealization,
    ReferenceMeasure,
    ProbabilityMeasure,
    QuantifierDomain,
    GaugeComposition,
    GaugeOrbitTypeProfile,
    GaugeHypothesis,
    GaugeGroupPresentation,
    GaugeOrbitPresentation,
    InfluenceComposition,
    GaugeQuotientProfile,
    FiberCardinalityProfile,
    FiberDimensionProfile,
    FunctionalModelSpace,
    GaugeQuotientMap,
    GaugeInvariantMap,
    GaugeGroupoidPresentation,
    GaugeReductionLaw,
    GaugeSubgroupCertificate,
    GaugeResidualAction,
    GaugeMeasureTransport,
    MeasureTransport,
    ParameterizedLikelihood,
    UnitDefinition,
    ForwardModelProductionBinding,
}

/// Unresolved content reference under an exact digest domain and source
/// contract version. A hash is a binding, not an authentication or scientific
/// correctness claim; authority is supplied only during source admission.
#[derive(Clone, PartialEq, Eq)]
pub struct SourceRef {
    key: SourceKey,
    kind: SourceKind,
    expected_hash: ContentHash,
    content_hash_domain: String,
    contract_version: u32,
}

impl fmt::Debug for SourceRef {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SourceRef")
            .field("kind", &self.kind)
            .field("contract_version", &self.contract_version)
            .finish_non_exhaustive()
    }
}

impl SourceRef {
    /// Construct a source reference with an exact fs-blake3 domain and source
    /// contract version.
    pub fn try_new(
        key: SourceKey,
        kind: SourceKind,
        expected_hash: ContentHash,
        content_hash_domain: impl Into<String>,
        contract_version: u32,
    ) -> Result<Self, IdentifiabilityError> {
        let content_hash_domain = content_hash_domain.into();
        validate_token(&content_hash_domain, "source content-hash domain")?;
        if !hash_is_nonzero(expected_hash) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "source reference",
            });
        }
        if contract_version == 0 {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "source contract",
                expected: 1,
                actual: 0,
            });
        }
        Ok(Self {
            key,
            kind,
            expected_hash,
            content_hash_domain,
            contract_version,
        })
    }

    /// Construct an exact typed reference to a Context-of-Use artifact.
    pub fn context(key: SourceKey, context: &ContextOfUse) -> Result<Self, IdentifiabilityError> {
        let expected_hash = context
            .content_hash()
            .map_err(|error| IdentifiabilityError::Vv {
                detail: error.to_string(),
            })?;
        Self::try_new(
            key,
            SourceKind::ContextOfUse,
            expected_hash,
            VV_ARTIFACT_SOURCE_DOMAIN,
            VV_SCHEMA_VERSION,
        )
    }

    /// Construct an exact typed reference to an experiment artifact.
    pub fn experiment(
        key: SourceKey,
        experiment: &ExperimentArtifact,
    ) -> Result<Self, IdentifiabilityError> {
        let expected_hash =
            experiment
                .content_hash()
                .map_err(|error| IdentifiabilityError::Vv {
                    detail: error.to_string(),
                })?;
        Self::try_new(
            key,
            SourceKind::ExperimentArtifact,
            expected_hash,
            VV_ARTIFACT_SOURCE_DOMAIN,
            VV_SCHEMA_VERSION,
        )
    }

    /// Construct an exact typed reference to a calibration split.
    pub fn calibration_split(
        key: SourceKey,
        split: &CalibrationSplit,
    ) -> Result<Self, IdentifiabilityError> {
        let expected_hash = split
            .content_hash()
            .map_err(|error| IdentifiabilityError::Vv {
                detail: error.to_string(),
            })?;
        Self::try_new(
            key,
            SourceKind::CalibrationSplit,
            expected_hash,
            VV_ARTIFACT_SOURCE_DOMAIN,
            VV_SCHEMA_VERSION,
        )
    }

    /// Construct an exact typed reference to a material card.
    pub fn material_card(
        key: SourceKey,
        material: &MaterialCard,
    ) -> Result<Self, IdentifiabilityError> {
        Self::try_new(
            key,
            SourceKind::MaterialCard,
            material.content_hash(),
            MATERIAL_CARD_SOURCE_DOMAIN,
            MATDB_SCHEMA_VERSION,
        )
    }

    /// Construct an exact typed reference to a constitutive-model card.
    pub fn constitutive_model_card(
        key: SourceKey,
        model: &ConstitutiveModelCard,
    ) -> Result<Self, IdentifiabilityError> {
        Self::try_new(
            key,
            SourceKind::ConstitutiveModelCard,
            model.content_hash(),
            CONSTITUTIVE_MODEL_CARD_SOURCE_DOMAIN,
            MATDB_SCHEMA_VERSION,
        )
    }

    #[must_use]
    pub const fn key(&self) -> &SourceKey {
        &self.key
    }

    #[must_use]
    pub const fn kind(&self) -> SourceKind {
        self.kind
    }

    #[must_use]
    pub const fn expected_hash(&self) -> ContentHash {
        self.expected_hash
    }

    #[must_use]
    pub fn content_hash_domain(&self) -> &str {
        &self.content_hash_domain
    }

    #[must_use]
    pub const fn contract_version(&self) -> u32 {
        self.contract_version
    }
}

/// Canonical preimage for an independently resolvable producer-to-forward-
/// model receipt. The binding covers the producer artifact and every semantic
/// field of the exact [`SourceRef`], not only a display label or bare digest.
pub fn forward_model_production_binding_preimage(
    producer: &ArtifactId,
    forward_model: &SourceRef,
) -> Result<Vec<u8>, IdentifiabilityError> {
    let mut writer = CanonicalWriter::new();
    writer.raw(b"fs-material-forward-model-production-binding\0");
    writer.u32(FORWARD_MODEL_PRODUCTION_BINDING_VERSION);
    encode_artifact_id(&mut writer, producer)?;
    encode_source_ref(&mut writer, forward_model)?;
    writer.finish()
}

/// Optional authentication envelope for a typed external trust receipt.
/// `Unauthenticated` is deliberately representable for legacy authority
/// records such as the current blind-release artifact, but it carries an
/// explicit no-authentication boundary rather than silently omitting issuer
/// and policy semantics.
#[derive(Clone, PartialEq, Eq)]
pub enum TrustAuthentication {
    Unauthenticated,
    IssuerPolicy {
        issuer: ArtifactId,
        trust_policy: SourceRef,
    },
}

impl fmt::Debug for TrustAuthentication {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Unauthenticated => "TrustAuthentication::Unauthenticated",
            Self::IssuerPolicy { .. } => "TrustAuthentication::IssuerPolicy(<redacted>)",
        })
    }
}

/// Domain-, version-, subject-, and optionally issuer/policy-bound trust
/// receipt identity. Equal 32-byte digests in different receipt systems or
/// for different subjects are intentionally distinct.
#[derive(Clone, PartialEq, Eq)]
pub struct TrustReceiptRef {
    receipt: SourceRef,
    subject: SourceRef,
    subject_artifact: Option<ArtifactId>,
    authentication: TrustAuthentication,
}

impl fmt::Debug for TrustReceiptRef {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrustReceiptRef")
            .field("receipt_kind", &self.receipt.kind)
            .field("contract_version", &self.receipt.contract_version)
            .field(
                "issuer_policy_declared",
                &matches!(
                    &self.authentication,
                    TrustAuthentication::IssuerPolicy { .. }
                ),
            )
            .field(
                "subject_artifact_declared",
                &self.subject_artifact.is_some(),
            )
            .finish_non_exhaustive()
    }
}

impl TrustReceiptRef {
    pub fn try_new(
        receipt: SourceRef,
        subject: SourceRef,
        authentication: TrustAuthentication,
    ) -> Result<Self, IdentifiabilityError> {
        Self::try_new_with_subject_artifact(receipt, subject, None, authentication)
    }

    fn try_new_with_subject_artifact(
        receipt: SourceRef,
        subject: SourceRef,
        subject_artifact: Option<ArtifactId>,
        authentication: TrustAuthentication,
    ) -> Result<Self, IdentifiabilityError> {
        if receipt.kind != SourceKind::EvidenceReceipt {
            return Err(IdentifiabilityError::InvalidText {
                field: "trust receipt kind",
                detail: "trust receipt identity must use EvidenceReceipt semantics".to_string(),
            });
        }
        if let TrustAuthentication::IssuerPolicy { trust_policy, .. } = &authentication
            && trust_policy.kind != SourceKind::Assumption
        {
            return Err(IdentifiabilityError::InvalidText {
                field: "trust receipt policy",
                detail: "issuer-policy-bound trust receipts need an exact Assumption policy source; this declaration does not itself verify the issuer"
                    .to_string(),
            });
        }
        let uses_blind_release_domain =
            receipt.content_hash_domain == BLIND_RELEASE_TRUST_RECEIPT_DOMAIN;
        if uses_blind_release_domain
            && receipt.contract_version != BLIND_RELEASE_TRUST_RECEIPT_VERSION
        {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "blind-release trust receipt",
                expected: BLIND_RELEASE_TRUST_RECEIPT_VERSION,
                actual: receipt.contract_version,
            });
        }
        let is_typed_blind_release = uses_blind_release_domain
            && receipt.contract_version == BLIND_RELEASE_TRUST_RECEIPT_VERSION;
        if is_typed_blind_release != subject_artifact.is_some()
            || (is_typed_blind_release && subject.kind != SourceKind::CalibrationSplit)
        {
            return Err(IdentifiabilityError::InvalidText {
                field: "trust receipt subject artifact",
                detail: "the typed blind-release namespace and an exact CalibrationSplit artifact subject must occur together; generic trust receipts may carry neither"
                    .to_string(),
            });
        }
        Ok(Self {
            receipt,
            subject,
            subject_artifact,
            authentication,
        })
    }

    pub fn blind_release(
        subject: &SourceRef,
        split: ArtifactId,
        receipt_hash: ContentHash,
    ) -> Result<Self, IdentifiabilityError> {
        Self::try_new_with_subject_artifact(
            SourceRef::try_new(
                subject.key.clone(),
                SourceKind::EvidenceReceipt,
                receipt_hash,
                BLIND_RELEASE_TRUST_RECEIPT_DOMAIN,
                BLIND_RELEASE_TRUST_RECEIPT_VERSION,
            )?,
            subject.clone(),
            Some(split),
            TrustAuthentication::Unauthenticated,
        )
    }

    #[must_use]
    pub const fn receipt(&self) -> &SourceRef {
        &self.receipt
    }

    #[must_use]
    pub const fn subject(&self) -> &SourceRef {
        &self.subject
    }

    #[must_use]
    pub const fn subject_artifact(&self) -> Option<&ArtifactId> {
        self.subject_artifact.as_ref()
    }

    #[must_use]
    pub const fn authentication(&self) -> &TrustAuthentication {
        &self.authentication
    }
}

/// Authority attached to a resolved source. `ContentVerified` proves only
/// byte identity. `ExternalTrustReceipt` retains typed external authority;
/// neither variant asserts scientific correctness, and an explicitly
/// unauthenticated receipt does not assert issuer or policy authentication.
#[derive(Clone, PartialEq, Eq)]
pub enum AuthorityDisposition {
    ContentVerified,
    ExternalTrustReceipt { trust_receipt: TrustReceiptRef },
    Unverified { reason: String },
}

impl fmt::Debug for AuthorityDisposition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::ContentVerified => "AuthorityDisposition::ContentVerified",
            Self::ExternalTrustReceipt { .. } => {
                "AuthorityDisposition::ExternalTrustReceipt(<redacted>)"
            }
            Self::Unverified { .. } => "AuthorityDisposition::Unverified(<redacted>)",
        })
    }
}

/// Evidence retained for the byte-identity part of source admission.
///
/// This deliberately says nothing about scientific correctness or source
/// authentication. [`HashPreimage`](Self::HashPreimage) proves only that the
/// supplied bytes reproduce the declared domain-separated digest; callers
/// must not infer that those bytes obey a canonical source schema.
#[derive(Clone, PartialEq, Eq)]
pub enum SourceVerification {
    TypedArtifact,
    HashPreimage { byte_len: u64 },
    Unverified,
}

impl fmt::Debug for SourceVerification {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::TypedArtifact => "SourceVerification::TypedArtifact",
            Self::HashPreimage { .. } => "SourceVerification::HashPreimage(<redacted>)",
            Self::Unverified => "SourceVerification::Unverified",
        })
    }
}

/// Resolution supplied for one opaque source key.
#[derive(Clone, PartialEq, Eq)]
pub struct SourceResolution {
    key: SourceKey,
    kind: SourceKind,
    resolved_hash: ContentHash,
    content_hash_domain: String,
    contract_version: u32,
    authority: AuthorityDisposition,
    verification: SourceVerification,
}

impl fmt::Debug for SourceResolution {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SourceResolution")
            .field("kind", &self.kind)
            .field("contract_version", &self.contract_version)
            .field("authority", &self.authority)
            .field("verification", &self.verification)
            .finish_non_exhaustive()
    }
}

impl SourceResolution {
    /// Verify an opaque source from retained hash-preimage bytes. The
    /// `SourceRef.content_hash_domain` is the fs-blake3 domain; callers cannot
    /// self-assert content equality without supplying bytes that reproduce the
    /// expected digest.
    pub fn verify(
        reference: &SourceRef,
        hash_preimage: &[u8],
        authority: AuthorityDisposition,
    ) -> Result<Self, IdentifiabilityError> {
        validate_authority_disposition(&authority)?;
        validate_authority_subject(reference, &authority)?;
        if hash_preimage.len() > MAX_IDENTIFIABILITY_CANONICAL_BYTES {
            return Err(IdentifiabilityError::Cardinality {
                field: "source hash preimage",
                detail: format!(
                    "synchronous verification is bounded to {MAX_IDENTIFIABILITY_CANONICAL_BYTES} bytes; larger artifacts require a typed or streaming artifact-store verifier"
                ),
            });
        }
        if matches!(&authority, AuthorityDisposition::Unverified { .. }) {
            return Err(IdentifiabilityError::InvalidText {
                field: "source verification",
                detail: "verified resolution cannot carry Unverified authority".to_string(),
            });
        }
        let actual = hash_domain(&reference.content_hash_domain, hash_preimage);
        if actual != reference.expected_hash {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "opaque source hash preimage",
            });
        }
        Ok(Self {
            key: reference.key.clone(),
            kind: reference.kind,
            resolved_hash: actual,
            content_hash_domain: reference.content_hash_domain.clone(),
            contract_version: reference.contract_version,
            authority,
            verification: SourceVerification::HashPreimage {
                byte_len: u64::try_from(hash_preimage.len()).map_err(|_| {
                    IdentifiabilityError::Cardinality {
                        field: "source hash preimage",
                        detail: "source preimage length exceeds u64".to_string(),
                    }
                })?,
            },
        })
    }

    /// Retain an explicit unresolved record for diagnostics. Admission rejects
    /// this variant deterministically.
    pub fn unresolved(
        reference: &SourceRef,
        reason: impl Into<String>,
    ) -> Result<Self, IdentifiabilityError> {
        let reason = reason.into();
        validate_reason(&reason, "unverified source reason")?;
        Ok(Self {
            key: reference.key.clone(),
            kind: reference.kind,
            resolved_hash: reference.expected_hash,
            content_hash_domain: reference.content_hash_domain.clone(),
            contract_version: reference.contract_version,
            authority: AuthorityDisposition::Unverified { reason },
            verification: SourceVerification::Unverified,
        })
    }

    #[must_use]
    pub const fn key(&self) -> &SourceKey {
        &self.key
    }

    #[must_use]
    pub const fn kind(&self) -> SourceKind {
        self.kind
    }

    #[must_use]
    pub const fn resolved_hash(&self) -> ContentHash {
        self.resolved_hash
    }

    #[must_use]
    pub fn content_hash_domain(&self) -> &str {
        &self.content_hash_domain
    }

    #[must_use]
    pub const fn contract_version(&self) -> u32 {
        self.contract_version
    }

    #[must_use]
    pub const fn authority(&self) -> &AuthorityDisposition {
        &self.authority
    }

    /// Inspect how byte identity was established without conferring additional
    /// scientific or trust authority.
    #[must_use]
    pub const fn verification(&self) -> &SourceVerification {
        &self.verification
    }
}

/// Exact opaque-source resolutions.  Duplicate keys are refused instead of
/// taking last-writer-wins authority.
#[derive(Clone, PartialEq, Eq, Default)]
pub struct SourceResolutionSet {
    entries: BTreeMap<SourceKey, SourceResolution>,
}

impl fmt::Debug for SourceResolutionSet {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SourceResolutionSet")
            .field("source_count", &self.entries.len())
            .finish()
    }
}

/// Retained preimage of source-authority identity.  This is deliberately
/// distinct from the physical problem document: trust-policy receipts may move
/// this record without rewriting the scientific question.
#[derive(Clone, PartialEq, Eq)]
struct SourceAdmissionRecord {
    schema_version: u32,
    problem_id: ProblemId,
    resolutions: BTreeMap<SourceKey, SourceResolution>,
}

impl fmt::Debug for SourceAdmissionRecord {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SourceAdmissionRecord")
            .field("schema_version", &self.schema_version)
            .field("source_count", &self.resolutions.len())
            .finish_non_exhaustive()
    }
}

impl SourceResolutionSet {
    /// Canonicalize a bounded resolution set.
    pub fn try_new(entries: Vec<SourceResolution>) -> Result<Self, IdentifiabilityError> {
        if entries.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "source resolutions",
                detail: "too many source resolutions".to_string(),
            });
        }
        let mut by_key = BTreeMap::new();
        for entry in entries {
            validate_authority_subject_fields(
                &entry.key,
                entry.kind,
                entry.resolved_hash,
                &entry.content_hash_domain,
                entry.contract_version,
                &entry.authority,
            )?;
            let key = entry.key.clone();
            if by_key.insert(key.clone(), entry).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "source resolution",
                    id: key.to_string(),
                });
            }
        }
        Ok(Self { entries: by_key })
    }

    #[must_use]
    pub const fn entries(&self) -> &BTreeMap<SourceKey, SourceResolution> {
        &self.entries
    }
}

/// Decision-facing role of a physical parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterPurpose {
    Estimand,
    Nuisance,
    Hyperparameter,
    CalibrationControl,
}

/// Exact value and provenance for a parameter conditioned outside inference.
#[derive(Debug, Clone, PartialEq)]
pub struct ConditionedValue {
    value_si: f64,
    source: SourceKey,
}

impl ConditionedValue {
    pub fn try_new(value_si: f64, source: SourceKey) -> Result<Self, IdentifiabilityError> {
        if !value_si.is_finite() {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "conditioned parameter value",
                detail: "value must be finite".to_string(),
            });
        }
        Ok(Self {
            value_si: canonical_f64(value_si),
            source,
        })
    }

    #[must_use]
    pub const fn value_si(&self) -> f64 {
        self.value_si
    }

    #[must_use]
    pub const fn source(&self) -> &SourceKey {
        &self.source
    }
}

/// Inferential treatment is orthogonal to decision purpose.
#[derive(Debug, Clone, PartialEq)]
pub enum ParameterTreatment {
    Estimated,
    Profiled,
    Marginalized,
    Conditioned(ConditionedValue),
    Derived {
        definition: SourceKey,
        parents: BTreeSet<ParameterRoleId>,
    },
}

/// Prior semantics distinguish absence from not-applicable.
#[derive(Debug, Clone, PartialEq)]
pub enum PriorPolicy {
    Distribution(ParameterPrior),
    Absent { reason: String },
    NotApplicable { reason: String },
}

/// Whether schema-level influence connectivity is declared.  This is not an
/// identifiability result and contains no evidence receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InfluenceCoverage {
    Declared,
    IntentionallyAbsent { reason: String },
    NotApplicable { reason: String },
}

/// Semantic owner with an immutable payload binding.  A bare category label is
/// insufficient to distinguish two instruments, discrepancy families, or
/// protocol controls.
#[derive(Clone, PartialEq, Eq)]
pub enum ParameterOwnerBinding {
    ConstitutiveModel,
    InitialState {
        state_path: SourceKey,
    },
    Instrument {
        instrument: ArtifactId,
        acquisition_channel: ArtifactId,
        metrology: SourceKey,
    },
    Discrepancy {
        family: SourceKey,
    },
    ControlledInput {
        protocol: SourceKey,
    },
    Population {
        hierarchy: SourceKey,
    },
}

impl fmt::Debug for ParameterOwnerBinding {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::ConstitutiveModel => "ParameterOwnerBinding::ConstitutiveModel",
            Self::InitialState { .. } => "ParameterOwnerBinding::InitialState(<redacted>)",
            Self::Instrument { .. } => "ParameterOwnerBinding::Instrument(<redacted>)",
            Self::Discrepancy { .. } => "ParameterOwnerBinding::Discrepancy(<redacted>)",
            Self::ControlledInput { .. } => "ParameterOwnerBinding::ControlledInput(<redacted>)",
            Self::Population { .. } => "ParameterOwnerBinding::Population(<redacted>)",
        })
    }
}

/// Population/realization scope, including explicit multi-case scope.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParameterScopeBinding {
    Global,
    Cases(BTreeSet<CaseId>),
    MaterialLot {
        lot: ArtifactId,
        cases: BTreeSet<CaseId>,
    },
    Specimen {
        case: CaseId,
        specimen: ArtifactId,
    },
    Field {
        support: SourceKey,
        cases: BTreeSet<CaseId>,
    },
    Hierarchical {
        population: ArtifactId,
        level: u32,
        hierarchy: SourceKey,
        cases: BTreeSet<CaseId>,
    },
}

/// Coordinate-free physical parameter declaration.
#[derive(Clone, PartialEq)]
pub struct StudyParameter {
    role: ParameterRoleId,
    quantity: QuantitySpec,
    domain: ParameterDomain,
    purpose: ParameterPurpose,
    treatment: ParameterTreatment,
    owner: ParameterOwnerBinding,
    scope: ParameterScopeBinding,
    prior: PriorPolicy,
    influence_coverage: InfluenceCoverage,
}

impl fmt::Debug for StudyParameter {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let treatment = match &self.treatment {
            ParameterTreatment::Estimated => "estimated",
            ParameterTreatment::Profiled => "profiled",
            ParameterTreatment::Marginalized => "marginalized",
            ParameterTreatment::Conditioned(_) => "conditioned",
            ParameterTreatment::Derived { .. } => "derived",
        };
        let scope = match &self.scope {
            ParameterScopeBinding::Global => "global",
            ParameterScopeBinding::Cases(_) => "cases",
            ParameterScopeBinding::MaterialLot { .. } => "material-lot",
            ParameterScopeBinding::Specimen { .. } => "specimen",
            ParameterScopeBinding::Field { .. } => "field",
            ParameterScopeBinding::Hierarchical { .. } => "hierarchical",
        };
        let prior = match &self.prior {
            PriorPolicy::Distribution(ParameterPrior::None { .. }) => "invalid-none",
            PriorPolicy::Distribution(ParameterPrior::Uniform { .. }) => "uniform",
            PriorPolicy::Distribution(ParameterPrior::Gaussian { .. }) => "gaussian",
            PriorPolicy::Distribution(ParameterPrior::LogNormal { .. }) => "log-normal",
            PriorPolicy::Absent { .. } => "absent",
            PriorPolicy::NotApplicable { .. } => "not-applicable",
        };
        let influence = match &self.influence_coverage {
            InfluenceCoverage::Declared => "declared",
            InfluenceCoverage::IntentionallyAbsent { .. } => "intentionally-absent",
            InfluenceCoverage::NotApplicable { .. } => "not-applicable",
        };
        formatter
            .debug_struct("StudyParameter")
            .field("role", &self.role)
            .field("purpose", &self.purpose)
            .field("treatment", &treatment)
            .field("owner", &self.owner)
            .field("scope", &scope)
            .field("prior_family", &prior)
            .field("influence_coverage", &influence)
            .finish_non_exhaustive()
    }
}

impl StudyParameter {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        role: ParameterRoleId,
        quantity: QuantitySpec,
        domain: ParameterDomain,
        purpose: ParameterPurpose,
        treatment: ParameterTreatment,
        owner: ParameterOwnerBinding,
        scope: ParameterScopeBinding,
        mut prior: PriorPolicy,
        influence_coverage: InfluenceCoverage,
    ) -> Result<Self, IdentifiabilityError> {
        if !domain.lo.is_finite() || !domain.hi.is_finite() || domain.lo > domain.hi {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "physical parameter domain",
                detail: format!("parameter {role} has invalid finite bounds"),
            });
        }
        match &mut prior {
            PriorPolicy::Distribution(ParameterPrior::None { .. }) => {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "prior policy",
                    detail: "PriorPolicy::Absent is the sole representation of no prior"
                        .to_string(),
                });
            }
            PriorPolicy::Distribution(distribution) => distribution.validate_against(domain)?,
            PriorPolicy::Absent { reason } => validate_reason(reason, "prior absence reason")?,
            PriorPolicy::NotApplicable { reason } => {
                validate_reason(reason, "prior not-applicable reason")?
            }
        }
        match &treatment {
            ParameterTreatment::Estimated
            | ParameterTreatment::Profiled
            | ParameterTreatment::Marginalized
                if domain.is_degenerate() =>
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "free parameter domain",
                    detail: format!("free parameter {role} needs a non-degenerate domain"),
                });
            }
            ParameterTreatment::Conditioned(value) => {
                if value.value_si < domain.lo || value.value_si > domain.hi {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "conditioned parameter value",
                        detail: format!("conditioned value for {role} lies outside its domain"),
                    });
                }
                if !matches!(&prior, PriorPolicy::NotApplicable { .. }) {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "conditioned parameter prior",
                        detail: format!(
                            "conditioned parameter {role} requires NotApplicable prior"
                        ),
                    });
                }
            }
            ParameterTreatment::Derived { parents, .. } => {
                if parents.is_empty() || parents.contains(&role) {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "derived parameter parent",
                        id: role.to_string(),
                    });
                }
                if !matches!(&prior, PriorPolicy::NotApplicable { .. }) {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "derived parameter prior",
                        detail: format!("derived parameter {role} requires NotApplicable prior"),
                    });
                }
            }
            _ => {}
        }
        if matches!(&treatment, ParameterTreatment::Marginalized)
            && !matches!(&prior, PriorPolicy::Distribution(_))
        {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "marginalized parameter prior",
                detail: "marginalization requires an explicit probability measure".to_string(),
            });
        }
        let free = matches!(
            &treatment,
            ParameterTreatment::Estimated
                | ParameterTreatment::Profiled
                | ParameterTreatment::Marginalized
        );
        match &influence_coverage {
            InfluenceCoverage::IntentionallyAbsent { reason } if free => {
                validate_reason(reason, "intentionally absent influence reason")?;
            }
            InfluenceCoverage::NotApplicable { reason } if !free => {
                validate_reason(reason, "not-applicable influence reason")?;
            }
            InfluenceCoverage::Declared => {}
            InfluenceCoverage::IntentionallyAbsent { .. } => {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "parameter influence coverage",
                    detail: "only free inference parameters may carry an influence no-claim"
                        .to_string(),
                });
            }
            InfluenceCoverage::NotApplicable { .. } => {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "parameter influence coverage",
                    detail: "free inference parameters cannot mark influence not applicable"
                        .to_string(),
                });
            }
        }
        Ok(Self {
            role,
            quantity,
            domain,
            purpose,
            treatment,
            owner,
            scope,
            prior,
            influence_coverage,
        })
    }

    #[must_use]
    pub const fn role(&self) -> &ParameterRoleId {
        &self.role
    }

    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    #[must_use]
    pub const fn domain(&self) -> ParameterDomain {
        self.domain
    }

    #[must_use]
    pub const fn treatment(&self) -> &ParameterTreatment {
        &self.treatment
    }

    #[must_use]
    pub const fn purpose(&self) -> ParameterPurpose {
        self.purpose
    }

    #[must_use]
    pub const fn owner(&self) -> &ParameterOwnerBinding {
        &self.owner
    }

    #[must_use]
    pub const fn scope(&self) -> &ParameterScopeBinding {
        &self.scope
    }

    #[must_use]
    pub const fn prior(&self) -> &PriorPolicy {
        &self.prior
    }

    #[must_use]
    pub const fn influence_coverage(&self) -> &InfluenceCoverage {
        &self.influence_coverage
    }
}

/// A typed scalar coefficient used by a joint affine constraint.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineConstraintTerm {
    parameter: ParameterRoleId,
    coefficient: f64,
    coefficient_quantity: QuantitySpec,
}

impl AffineConstraintTerm {
    pub fn try_new(
        parameter: ParameterRoleId,
        coefficient: f64,
        coefficient_quantity: QuantitySpec,
    ) -> Result<Self, IdentifiabilityError> {
        if !coefficient.is_finite() {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "constraint coefficient",
                detail: "coefficient must be finite".to_string(),
            });
        }
        Ok(Self {
            parameter,
            coefficient: canonical_f64(coefficient),
            coefficient_quantity,
        })
    }

    #[must_use]
    pub const fn parameter(&self) -> &ParameterRoleId {
        &self.parameter
    }

    #[must_use]
    pub const fn coefficient(&self) -> f64 {
        self.coefficient
    }

    #[must_use]
    pub const fn coefficient_quantity(&self) -> QuantitySpec {
        self.coefficient_quantity
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintRelation {
    Equal,
    LessOrEqual,
    GreaterOrEqual,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstraintCodimension {
    Finite {
        codimension: u64,
    },
    /// Functional/continuum codimension whose exact analytic category and
    /// closed-range/transversality semantics are external to this crate.
    InfiniteDimensional {
        profile: SourceKey,
    },
}

/// Cross-parameter admissible-domain constraint.
#[derive(Debug, Clone, PartialEq)]
pub enum JointConstraintKind {
    Affine {
        terms: Vec<AffineConstraintTerm>,
        relation: ConstraintRelation,
        rhs_si: f64,
        residual_quantity: QuantitySpec,
    },
    Simplex {
        members: BTreeSet<ParameterRoleId>,
        total_si: f64,
        quantity: QuantitySpec,
    },
    Ordered {
        members: Vec<ParameterRoleId>,
        strict: bool,
    },
    ExternalManifold {
        members: BTreeSet<ParameterRoleId>,
        definition: SourceKey,
        codimension: ConstraintCodimension,
    },
    StochasticCoupling {
        members: BTreeSet<ParameterRoleId>,
        distribution: SourceKey,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct JointConstraint {
    id: ConstraintId,
    kind: JointConstraintKind,
}

impl JointConstraint {
    #[must_use]
    pub fn new(id: ConstraintId, mut kind: JointConstraintKind) -> Self {
        match &mut kind {
            JointConstraintKind::Affine { rhs_si, .. } => {
                *rhs_si = canonical_f64(*rhs_si);
            }
            JointConstraintKind::Simplex { total_si, .. } => {
                *total_si = canonical_f64(*total_si);
            }
            JointConstraintKind::Ordered { .. }
            | JointConstraintKind::ExternalManifold { .. }
            | JointConstraintKind::StochasticCoupling { .. } => {}
        }
        Self { id, kind }
    }

    #[must_use]
    pub const fn id(&self) -> &ConstraintId {
        &self.id
    }

    #[must_use]
    pub const fn kind(&self) -> &JointConstraintKind {
        &self.kind
    }
}

/// Constructive evidence that the declared parameter domain is not empty.
///
/// `values` supplies one exact finite point for every parameter. Built-in
/// affine, simplex, and ordering constraints are checked directly against that
/// point. Constraints whose semantics are intentionally external cannot be
/// evaluated in this crate, so each one additionally requires an exact,
/// source-admitted membership claim bound to this witness and the *full*
/// canonical constraint conjunction. Those claims are evidence inputs, not
/// theorem tokens; their scientific validity remains a downstream verifier
/// responsibility.
#[derive(Clone, PartialEq, Eq)]
pub struct OpaqueDomainMembershipClaim {
    source: SourceKey,
    witness_binding: Option<ContentHash>,
}

impl fmt::Debug for OpaqueDomainMembershipClaim {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("OpaqueDomainMembershipClaim")
            .field("bound", &self.witness_binding.is_some())
            .finish_non_exhaustive()
    }
}

impl OpaqueDomainMembershipClaim {
    #[must_use]
    pub const fn new(source: SourceKey) -> Self {
        Self {
            source,
            witness_binding: None,
        }
    }

    const fn from_bound_source(source: SourceKey, witness_binding: ContentHash) -> Self {
        Self {
            source,
            witness_binding: Some(witness_binding),
        }
    }

    #[must_use]
    pub const fn source(&self) -> &SourceKey {
        &self.source
    }

    #[must_use]
    pub const fn witness_binding(&self) -> Option<ContentHash> {
        self.witness_binding
    }
}

#[derive(Clone, PartialEq)]
pub struct AdmissibleDomainWitness {
    values: BTreeMap<ParameterRoleId, f64>,
    opaque_membership_claim: Option<OpaqueDomainMembershipClaim>,
}

impl fmt::Debug for AdmissibleDomainWitness {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AdmissibleDomainWitness")
            .field("value_count", &self.values.len())
            .field(
                "opaque_membership_present",
                &self.opaque_membership_claim.is_some(),
            )
            .field(
                "opaque_membership_bound",
                &self
                    .opaque_membership_claim
                    .as_ref()
                    .is_some_and(|claim| claim.witness_binding.is_some()),
            )
            .finish_non_exhaustive()
    }
}

impl AdmissibleDomainWitness {
    pub fn try_new(
        values: Vec<(ParameterRoleId, f64)>,
        opaque_membership_claim: Option<OpaqueDomainMembershipClaim>,
    ) -> Result<Self, IdentifiabilityError> {
        if values.is_empty() || values.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "admissible-domain witness values",
                detail: "the witness needs bounded nonempty parameter values".to_string(),
            });
        }
        let mut canonical_values = BTreeMap::new();
        for (role, value) in values {
            if !value.is_finite() {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "admissible-domain witness value",
                    detail: format!("parameter {role} has a non-finite witness value"),
                });
            }
            if canonical_values
                .insert(role.clone(), canonical_f64(value))
                .is_some()
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "admissible-domain witness parameter",
                    id: role.to_string(),
                });
            }
        }
        if let Some(binding) = opaque_membership_claim
            .as_ref()
            .and_then(|claim| claim.witness_binding)
            && !hash_is_nonzero(binding)
        {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "admissible-domain witness binding",
            });
        }
        Ok(Self {
            values: canonical_values,
            opaque_membership_claim,
        })
    }

    #[must_use]
    pub const fn values(&self) -> &BTreeMap<ParameterRoleId, f64> {
        &self.values
    }

    #[must_use]
    pub const fn opaque_membership_claim(&self) -> Option<&OpaqueDomainMembershipClaim> {
        self.opaque_membership_claim.as_ref()
    }

    fn bind_opaque_membership(
        &mut self,
        expected: ContentHash,
    ) -> Result<(), IdentifiabilityError> {
        let Some(claim) = &mut self.opaque_membership_claim else {
            return Ok(());
        };
        match claim.witness_binding {
            None => claim.witness_binding = Some(expected),
            Some(actual) if actual == expected => {}
            Some(_) => {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "admissible-domain witness/conjunction binding",
                });
            }
        }
        Ok(())
    }
}

/// Why a case participates in the campaign.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CasePurpose {
    Calibration,
    SymmetryBreaking,
    ValidationOnly,
    BlindFalsification,
    ProspectiveDesign,
    Complementary { reason: String },
}

/// Whether observations already exist.  Retrospective lineage is re-derived
/// from concrete V&V artifacts at admission; it is never trusted from bytes.
#[derive(Clone, PartialEq, Eq)]
pub enum CaseDataDeclaration {
    Prospective,
    Retrospective {
        experiment: SourceKey,
        split: SourceKey,
        parser: SourceKey,
        preprocessing: SourceKey,
        parser_version: u32,
        split_grouping: ArtifactId,
    },
}

impl fmt::Debug for CaseDataDeclaration {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Prospective => formatter.write_str("CaseDataDeclaration::Prospective"),
            Self::Retrospective { parser_version, .. } => formatter
                .debug_struct("CaseDataDeclaration::Retrospective")
                .field("parser_version", parser_version)
                .finish_non_exhaustive(),
        }
    }
}

/// Raw-row declaration for one channel.
#[derive(Clone, PartialEq, Eq)]
pub enum ObservationRows {
    Prospective,
    Retrospective(BTreeSet<ObservationId>),
}

impl fmt::Debug for ObservationRows {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (variant, row_count) = match self {
            Self::Prospective => ("prospective", 0),
            Self::Retrospective(rows) => ("retrospective", rows.len()),
        };
        formatter
            .debug_struct("ObservationRows")
            .field("variant", &variant)
            .field("row_count", &row_count)
            .finish()
    }
}

/// Marginal noise family.  Joint dependence is modeled separately so bounded
/// or unknown marginals are never silently assigned a standard deviation.
#[derive(Debug, Clone, PartialEq)]
pub enum MarginalNoiseSpec {
    Gaussian {
        standard_deviation: f64,
    },
    StudentT {
        scale: f64,
        degrees_of_freedom: f64,
    },
    Empirical {
        distribution: SourceKey,
        standard_deviation: f64,
        finite_variance_model: SourceKey,
    },
    Bounded {
        half_width: f64,
    },
    Unknown {
        reason: String,
    },
}

impl MarginalNoiseSpec {
    fn finite_standard_deviation(&self) -> bool {
        match self {
            Self::Gaussian { standard_deviation } => {
                standard_deviation.is_finite() && *standard_deviation > 0.0
            }
            Self::StudentT {
                scale,
                degrees_of_freedom,
            } => scale.is_finite() && *scale > 0.0 && *degrees_of_freedom > 2.0,
            Self::Empirical {
                standard_deviation, ..
            } => standard_deviation.is_finite() && *standard_deviation > 0.0,
            Self::Bounded { .. } | Self::Unknown { .. } => false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MissingnessAssumption {
    Complete { assumption: SourceKey },
    Modeled { mechanism: SourceKey },
    Unknown { reason: String },
}

/// Evidence-free physical observation schema.
#[derive(Clone, PartialEq)]
pub struct StudyObservation {
    id: ObservationChannelId,
    qoi: QoiId,
    unit: UnitId,
    quantity: QuantitySpec,
    unit_definition: SourceKey,
    frame: FrameBinding,
    graph_node: String,
    graph_port: String,
    operator: SourceKey,
    aggregation: SourceKey,
    sensor: SourceKey,
    instrument: ArtifactId,
    acquisition_channel: ArtifactId,
    clock: ArtifactId,
    operator_version: u32,
    noise: MarginalNoiseSpec,
    missingness: MissingnessAssumption,
    saturation: Option<ParameterDomain>,
    protocol_version: u32,
    refinement_version: u32,
    rows: ObservationRows,
}

impl fmt::Debug for StudyObservation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let noise_family = match &self.noise {
            MarginalNoiseSpec::Gaussian { .. } => "gaussian",
            MarginalNoiseSpec::StudentT { .. } => "student-t",
            MarginalNoiseSpec::Empirical { .. } => "empirical",
            MarginalNoiseSpec::Bounded { .. } => "bounded",
            MarginalNoiseSpec::Unknown { .. } => "unknown",
        };
        let missingness_family = match &self.missingness {
            MissingnessAssumption::Complete { .. } => "complete",
            MissingnessAssumption::Modeled { .. } => "modeled",
            MissingnessAssumption::Unknown { .. } => "unknown",
        };
        formatter
            .debug_struct("StudyObservation")
            .field("rows", &self.rows)
            .field("noise_family", &noise_family)
            .field("missingness_family", &missingness_family)
            .field("saturation_declared", &self.saturation.is_some())
            .finish_non_exhaustive()
    }
}

impl StudyObservation {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: ObservationChannelId,
        qoi: QoiId,
        unit: UnitId,
        quantity: QuantitySpec,
        unit_definition: SourceKey,
        frame: FrameBinding,
        graph_node: impl Into<String>,
        graph_port: impl Into<String>,
        operator: SourceKey,
        aggregation: SourceKey,
        sensor: SourceKey,
        instrument: ArtifactId,
        acquisition_channel: ArtifactId,
        clock: ArtifactId,
        operator_version: u32,
        mut noise: MarginalNoiseSpec,
        missingness: MissingnessAssumption,
        saturation: Option<ParameterDomain>,
        protocol_version: u32,
        refinement_version: u32,
        rows: ObservationRows,
    ) -> Result<Self, IdentifiabilityError> {
        let graph_node = graph_node.into();
        let graph_port = graph_port.into();
        validate_token(&graph_node, "observation graph node")?;
        validate_token(&graph_port, "observation graph port")?;
        if operator_version == 0 || protocol_version == 0 || refinement_version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "observation versions",
                detail: "operator, protocol, and refinement versions must be positive".to_string(),
            });
        }
        match &mut noise {
            MarginalNoiseSpec::Gaussian { standard_deviation }
                if standard_deviation.is_finite() && *standard_deviation > 0.0 => {}
            MarginalNoiseSpec::StudentT {
                scale,
                degrees_of_freedom,
            } if scale.is_finite()
                && *scale > 0.0
                && degrees_of_freedom.is_finite()
                && *degrees_of_freedom > 0.0 => {}
            MarginalNoiseSpec::Empirical {
                standard_deviation, ..
            } if standard_deviation.is_finite() && *standard_deviation > 0.0 => {}
            MarginalNoiseSpec::Bounded { half_width }
                if half_width.is_finite() && *half_width >= 0.0 =>
            {
                *half_width = canonical_f64(*half_width);
            }
            MarginalNoiseSpec::Unknown { reason } => {
                validate_reason(reason, "unknown noise reason")?;
            }
            _ => {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "marginal noise",
                    detail: "noise parameters must be finite and physically admissible".to_string(),
                });
            }
        }
        if let MissingnessAssumption::Unknown { reason } = &missingness {
            validate_reason(reason, "unknown missingness reason")?;
        }
        if let ObservationRows::Retrospective(rows) = &rows {
            if rows.is_empty() || rows.len() > MAX_IDENTIFIABILITY_ITEMS {
                return Err(IdentifiabilityError::Cardinality {
                    field: "observation rows",
                    detail: "retrospective channels need bounded nonempty raw-row sets".to_string(),
                });
            }
        }
        Ok(Self {
            id,
            qoi,
            unit,
            quantity,
            unit_definition,
            frame,
            graph_node,
            graph_port,
            operator,
            aggregation,
            sensor,
            instrument,
            acquisition_channel,
            clock,
            operator_version,
            noise,
            missingness,
            saturation,
            protocol_version,
            refinement_version,
            rows,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ObservationChannelId {
        &self.id
    }

    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn unit(&self) -> &UnitId {
        &self.unit
    }

    #[must_use]
    pub const fn unit_definition(&self) -> &SourceKey {
        &self.unit_definition
    }

    #[must_use]
    pub const fn frame(&self) -> &FrameBinding {
        &self.frame
    }

    #[must_use]
    pub fn graph_node(&self) -> &str {
        &self.graph_node
    }

    #[must_use]
    pub fn graph_port(&self) -> &str {
        &self.graph_port
    }

    #[must_use]
    pub const fn operator(&self) -> &SourceKey {
        &self.operator
    }

    #[must_use]
    pub const fn aggregation(&self) -> &SourceKey {
        &self.aggregation
    }

    #[must_use]
    pub const fn sensor(&self) -> &SourceKey {
        &self.sensor
    }

    #[must_use]
    pub const fn instrument(&self) -> &ArtifactId {
        &self.instrument
    }

    /// Physical/data-acquisition channel bound by every consumed manifest row.
    #[must_use]
    pub const fn acquisition_channel(&self) -> &ArtifactId {
        &self.acquisition_channel
    }

    #[must_use]
    pub const fn clock(&self) -> &ArtifactId {
        &self.clock
    }

    #[must_use]
    pub const fn operator_version(&self) -> u32 {
        self.operator_version
    }

    #[must_use]
    pub const fn noise(&self) -> &MarginalNoiseSpec {
        &self.noise
    }

    #[must_use]
    pub const fn missingness(&self) -> &MissingnessAssumption {
        &self.missingness
    }

    #[must_use]
    pub const fn saturation(&self) -> Option<ParameterDomain> {
        self.saturation
    }

    #[must_use]
    pub const fn protocol_version(&self) -> u32 {
        self.protocol_version
    }

    #[must_use]
    pub const fn refinement_version(&self) -> u32 {
        self.refinement_version
    }

    #[must_use]
    pub const fn rows(&self) -> &ObservationRows {
        &self.rows
    }
}

/// Joint noise/correlation semantics over composite observation keys.
#[derive(Debug, Clone, PartialEq)]
pub enum JointNoiseModel {
    Independent {
        assumption: SourceKey,
    },
    DenseCorrelation {
        order: Vec<ObservationKey>,
        correlation: CovarianceMatrix,
        model: SourceKey,
    },
    ExternalKernel {
        model: SourceKey,
    },
    Unknown {
        reason: String,
    },
}

/// Discrepancy is never inferred from absence.  Even an assumed-zero model is
/// an explicit, source-bound assumption rather than evidence of correctness.
#[derive(Clone, PartialEq, Eq)]
pub enum DiscrepancyInapplicability {
    /// A physical retrospective channel may omit a discrepancy term only
    /// under an exact, admitted applicability assumption. The assumption is
    /// not promoted to a theorem by this transport layer.
    PhysicalApplicability { assumption: SourceKey },
    /// The problem *declares* that data came from the same forward-model source
    /// used for inference. The constructor binds producer identity and every
    /// `SourceRef` field through a domain-separated production binding, while
    /// the explicit assumption supplies the conditional no-discrepancy claim.
    /// This is not authenticated run provenance: theorem promotion remains
    /// gated on a future signed execution receipt binding experiment hash,
    /// executable/build, producer, and the exact forward-model source.
    DeclaredSyntheticSelfModel {
        generator: SourceKey,
        producer: ArtifactId,
        production_binding: SourceKey,
        assumption: SourceKey,
    },
    /// No observations exist yet. This branch is restricted to an exact
    /// prospective-design case and cannot support a noisy finite-data or
    /// observed-posterior claim.
    ProspectiveDesign { assumption: SourceKey },
}

impl fmt::Debug for DiscrepancyInapplicability {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::PhysicalApplicability { .. } => {
                "DiscrepancyInapplicability::PhysicalApplicability(<redacted>)"
            }
            Self::DeclaredSyntheticSelfModel { .. } => {
                "DiscrepancyInapplicability::DeclaredSyntheticSelfModel(<redacted>)"
            }
            Self::ProspectiveDesign { .. } => {
                "DiscrepancyInapplicability::ProspectiveDesign(<redacted>)"
            }
        })
    }
}

#[derive(Clone, PartialEq, Eq)]
pub enum StudyDiscrepancy {
    Uncharacterized {
        reason: String,
    },
    NotApplicable {
        basis: DiscrepancyInapplicability,
    },
    AssumedZero {
        assumption: SourceKey,
    },
    Modeled {
        family: SourceKey,
        parameters: BTreeSet<ParameterRoleId>,
        support: SourceKey,
        confounding_guard: SourceKey,
    },
}

impl fmt::Debug for StudyDiscrepancy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Uncharacterized { .. } => {
                formatter.write_str("StudyDiscrepancy::Uncharacterized(<redacted>)")
            }
            Self::NotApplicable { basis } => formatter
                .debug_tuple("StudyDiscrepancy::NotApplicable")
                .field(basis)
                .finish(),
            Self::AssumedZero { .. } => {
                formatter.write_str("StudyDiscrepancy::AssumedZero(<redacted>)")
            }
            Self::Modeled { parameters, .. } => formatter
                .debug_struct("StudyDiscrepancy::Modeled")
                .field("parameter_count", &parameters.len())
                .finish_non_exhaustive(),
        }
    }
}

/// Source-registry keys that close the byte authority of every case-physics
/// component embedded in [`StudyCaseDocument`].
///
/// These references are intentionally separate from the bare component
/// digests in the current physical bindings. The document constructor proves
/// exact digest, role-specific domain, contract-version, and source-kind
/// equality; source admission then requires a typed or hash-preimage-verified
/// resolution for every key before a [`ProblemId`] can exist.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CasePhysicsSources {
    frame_transform: SourceKey,
    specimen_geometry: SourceKey,
    specimen_process: SourceKey,
    specimen_preparation: SourceKey,
    load_path: SourceKey,
    environment_path: SourceKey,
    time_grid: SourceKey,
    initial_state: Option<SourceKey>,
}

impl CasePhysicsSources {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub const fn new(
        frame_transform: SourceKey,
        specimen_geometry: SourceKey,
        specimen_process: SourceKey,
        specimen_preparation: SourceKey,
        load_path: SourceKey,
        environment_path: SourceKey,
        time_grid: SourceKey,
        initial_state: Option<SourceKey>,
    ) -> Self {
        Self {
            frame_transform,
            specimen_geometry,
            specimen_process,
            specimen_preparation,
            load_path,
            environment_path,
            time_grid,
            initial_state,
        }
    }

    #[must_use]
    pub const fn frame_transform(&self) -> &SourceKey {
        &self.frame_transform
    }

    #[must_use]
    pub const fn specimen_geometry(&self) -> &SourceKey {
        &self.specimen_geometry
    }

    #[must_use]
    pub const fn specimen_process(&self) -> &SourceKey {
        &self.specimen_process
    }

    #[must_use]
    pub const fn specimen_preparation(&self) -> &SourceKey {
        &self.specimen_preparation
    }

    #[must_use]
    pub const fn load_path(&self) -> &SourceKey {
        &self.load_path
    }

    #[must_use]
    pub const fn environment_path(&self) -> &SourceKey {
        &self.environment_path
    }

    #[must_use]
    pub const fn time_grid(&self) -> &SourceKey {
        &self.time_grid
    }

    #[must_use]
    pub const fn initial_state(&self) -> Option<&SourceKey> {
        self.initial_state.as_ref()
    }
}

/// Exact within-case declaration for intentionally consuming one raw row in
/// more than one observation channel. The group is semantic: it names every
/// repeated row, its exact consumer-channel set, and the joint-likelihood
/// source that prevents duplicate likelihood factors from being invented.
#[derive(Clone, PartialEq, Eq)]
pub struct ObservationSharingGroup {
    channels: BTreeSet<ObservationChannelId>,
    rows: BTreeSet<ObservationId>,
    joint_likelihood: SourceKey,
    justification: String,
}

impl fmt::Debug for ObservationSharingGroup {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ObservationSharingGroup")
            .field("channel_count", &self.channels.len())
            .field("row_count", &self.rows.len())
            .finish_non_exhaustive()
    }
}

impl ObservationSharingGroup {
    pub fn try_new(
        channels: BTreeSet<ObservationChannelId>,
        rows: BTreeSet<ObservationId>,
        joint_likelihood: SourceKey,
        justification: impl Into<String>,
    ) -> Result<Self, IdentifiabilityError> {
        if channels.len() < 2 || channels.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observation-sharing channels",
                detail: "a sharing group needs at least two bounded channels".to_string(),
            });
        }
        if rows.is_empty() || rows.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observation-sharing rows",
                detail: "a sharing group needs bounded nonempty raw rows".to_string(),
            });
        }
        let justification = justification.into();
        validate_reason(&justification, "observation-sharing justification")?;
        Ok(Self {
            channels,
            rows,
            joint_likelihood,
            justification,
        })
    }

    #[must_use]
    pub const fn channels(&self) -> &BTreeSet<ObservationChannelId> {
        &self.channels
    }

    #[must_use]
    pub const fn rows(&self) -> &BTreeSet<ObservationId> {
        &self.rows
    }

    #[must_use]
    pub const fn joint_likelihood(&self) -> &SourceKey {
        &self.joint_likelihood
    }

    #[must_use]
    pub fn justification(&self) -> &str {
        &self.justification
    }
}

/// One physical or prospective campaign case.
#[derive(Clone, PartialEq)]
pub struct StudyCaseDocument {
    id: CaseId,
    purpose: CasePurpose,
    initial_state: InitialStateBinding,
    specimen: SpecimenBinding,
    protocol: ProtocolBinding,
    physics_sources: CasePhysicsSources,
    forward_model: SourceKey,
    data: CaseDataDeclaration,
    observations: BTreeMap<ObservationChannelId, StudyObservation>,
    discrepancies: BTreeMap<ObservationChannelId, StudyDiscrepancy>,
    observation_sharing: Vec<ObservationSharingGroup>,
}

impl fmt::Debug for StudyCaseDocument {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let purpose = match &self.purpose {
            CasePurpose::Calibration => "calibration",
            CasePurpose::SymmetryBreaking => "symmetry-breaking",
            CasePurpose::ValidationOnly => "validation-only",
            CasePurpose::BlindFalsification => "blind-falsification",
            CasePurpose::ProspectiveDesign => "prospective-design",
            CasePurpose::Complementary { .. } => "complementary",
        };
        let data_phase = match &self.data {
            CaseDataDeclaration::Prospective => "prospective",
            CaseDataDeclaration::Retrospective { .. } => "retrospective",
        };
        formatter
            .debug_struct("StudyCaseDocument")
            .field("purpose", &purpose)
            .field("data_phase", &data_phase)
            .field("observation_count", &self.observations.len())
            .field("discrepancy_count", &self.discrepancies.len())
            .field("sharing_group_count", &self.observation_sharing.len())
            .finish_non_exhaustive()
    }
}

impl StudyCaseDocument {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: CaseId,
        purpose: CasePurpose,
        initial_state: InitialStateBinding,
        specimen: SpecimenBinding,
        protocol: ProtocolBinding,
        physics_sources: CasePhysicsSources,
        forward_model: SourceKey,
        data: CaseDataDeclaration,
        observations: Vec<StudyObservation>,
        discrepancies: Vec<(ObservationChannelId, StudyDiscrepancy)>,
        mut observation_sharing: Vec<ObservationSharingGroup>,
    ) -> Result<Self, IdentifiabilityError> {
        if let CasePurpose::Complementary { reason } = &purpose {
            validate_reason(reason, "complementary case reason")?;
        }
        if observations.is_empty() || observations.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "case observations",
                detail: "each case needs bounded nonempty observations".to_string(),
            });
        }
        let mut observation_map = BTreeMap::new();
        for observation in observations {
            let channel = observation.id.clone();
            if observation_map
                .insert(channel.clone(), observation)
                .is_some()
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "case observation",
                    id: channel.to_string(),
                });
            }
        }
        let mut discrepancy_map = BTreeMap::new();
        for (channel, discrepancy) in discrepancies {
            if !observation_map.contains_key(&channel) {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "discrepancy observation",
                    id: channel.to_string(),
                });
            }
            if discrepancy_map
                .insert(channel.clone(), discrepancy)
                .is_some()
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "discrepancy observation",
                    id: channel.to_string(),
                });
            }
        }
        if discrepancy_map.len() != observation_map.len() {
            return Err(IdentifiabilityError::Cardinality {
                field: "case discrepancies",
                detail: "every observation needs explicit discrepancy semantics".to_string(),
            });
        }
        observation_sharing.sort_by(|left, right| {
            (
                &left.rows,
                &left.channels,
                &left.joint_likelihood,
                &left.justification,
            )
                .cmp(&(
                    &right.rows,
                    &right.channels,
                    &right.joint_likelihood,
                    &right.justification,
                ))
        });
        if observation_sharing.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observation-sharing groups",
                detail: "too many within-case sharing groups".to_string(),
            });
        }
        let mut row_consumers = BTreeMap::<ObservationId, BTreeSet<ObservationChannelId>>::new();
        for (channel, observation) in &observation_map {
            if let ObservationRows::Retrospective(rows) = &observation.rows {
                for row in rows {
                    row_consumers
                        .entry(row.clone())
                        .or_default()
                        .insert(channel.clone());
                }
            }
        }
        let mut declared_rows = BTreeSet::new();
        for group in &observation_sharing {
            for channel in &group.channels {
                if !observation_map.contains_key(channel) {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "observation-sharing channel",
                        id: channel.to_string(),
                    });
                }
            }
            for row in &group.rows {
                if !declared_rows.insert(row.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "observation-sharing row",
                        id: row.as_str().to_string(),
                    });
                }
                if row_consumers.get(row) != Some(&group.channels) {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "observation-sharing consumers",
                        detail: format!(
                            "row {} is not consumed by exactly the declared channel set",
                            row.as_str()
                        ),
                    });
                }
            }
        }
        Ok(Self {
            id,
            purpose,
            initial_state,
            specimen,
            protocol,
            physics_sources,
            forward_model,
            data,
            observations: observation_map,
            discrepancies: discrepancy_map,
            observation_sharing,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &CaseId {
        &self.id
    }

    #[must_use]
    pub const fn observations(&self) -> &BTreeMap<ObservationChannelId, StudyObservation> {
        &self.observations
    }

    #[must_use]
    pub const fn purpose(&self) -> &CasePurpose {
        &self.purpose
    }

    #[must_use]
    pub const fn initial_state(&self) -> InitialStateBinding {
        self.initial_state
    }

    #[must_use]
    pub const fn specimen(&self) -> &SpecimenBinding {
        &self.specimen
    }

    #[must_use]
    pub const fn protocol(&self) -> &ProtocolBinding {
        &self.protocol
    }

    /// Exact source keys closing every embedded case-physics digest.
    #[must_use]
    pub const fn physics_sources(&self) -> &CasePhysicsSources {
        &self.physics_sources
    }

    #[must_use]
    pub const fn forward_model(&self) -> &SourceKey {
        &self.forward_model
    }

    #[must_use]
    pub const fn data(&self) -> &CaseDataDeclaration {
        &self.data
    }

    #[must_use]
    pub const fn discrepancies(&self) -> &BTreeMap<ObservationChannelId, StudyDiscrepancy> {
        &self.discrepancies
    }

    #[must_use]
    pub fn observation_sharing(&self) -> &[ObservationSharingGroup] {
        &self.observation_sharing
    }
}

/// Exact observation-distribution functional whose parameter dependence is
/// part of the physical question.  The derivative quantity is derived from
/// endpoint quantities and therefore cannot be supplied inconsistently.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum DistributionFunctional {
    Location {
        observation: ObservationKey,
    },
    LogScale {
        observation: ObservationKey,
    },
    Correlation {
        left: ObservationKey,
        right: ObservationKey,
    },
    MissingnessLogit {
        observation: ObservationKey,
    },
    CensoringLogit {
        observation: ObservationKey,
    },
}

/// Structural representation of an influence declaration.  Receipts proving
/// nonzero influence belong to an assessment, not here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InfluenceRepresentation {
    Direct,
    StateMediated {
        state_path: SourceKey,
    },
    Composite {
        operator: SourceKey,
        inputs: BTreeSet<InfluenceId>,
    },
    ExternalDefinition {
        definition: SourceKey,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InfluenceDeclaration {
    id: InfluenceId,
    parameter: ParameterRoleId,
    functional: DistributionFunctional,
    representation: InfluenceRepresentation,
}

impl InfluenceDeclaration {
    #[must_use]
    pub const fn new(
        id: InfluenceId,
        parameter: ParameterRoleId,
        functional: DistributionFunctional,
        representation: InfluenceRepresentation,
    ) -> Self {
        Self {
            id,
            parameter,
            functional,
            representation,
        }
    }

    #[must_use]
    pub const fn id(&self) -> &InfluenceId {
        &self.id
    }

    #[must_use]
    pub const fn parameter(&self) -> &ParameterRoleId {
        &self.parameter
    }

    #[must_use]
    pub const fn functional(&self) -> &DistributionFunctional {
        &self.functional
    }

    #[must_use]
    pub const fn representation(&self) -> &InfluenceRepresentation {
        &self.representation
    }
}

/// Algebraic size of the acting symmetry group. Group dimension/order is not
/// orbit dimension/cardinality: stabilizers are represented independently by
/// [`GaugeOrbitGeometry`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeDiscreteSize {
    Finite { order: u64 },
    CountablyInfinite { presentation: SourceKey },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeContinuousDimension {
    Finite {
        dimension: u64,
    },
    /// Infinite-dimensional Lie/model space. The source names exact
    /// Banach/Frechet/Sobolev semantics; it is not a basis-cardinality label.
    InfiniteDimensional {
        model_space: SourceKey,
    },
}

impl GaugeContinuousDimension {
    fn is_zero(&self) -> bool {
        matches!(self, Self::Finite { dimension: 0 })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeAlgebra {
    Continuous {
        group_dimension: GaugeContinuousDimension,
    },
    Discrete {
        size: GaugeDiscreteSize,
    },
    Mixed {
        continuous_group_dimension: GaugeContinuousDimension,
        component_group: GaugeDiscreteSize,
    },
}

/// Principal regular-orbit invariants after stabilizers. Only finite acting
/// and orbit dimensions admit the arithmetic stabilizer dimension
/// `group_dimension - orbit_dimension`; no subtraction theorem is claimed for
/// infinite-dimensional model spaces. The discrete cardinality is the
/// component-group orbit index, not the full stabilizer topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeDiscreteOrbitCardinality {
    Finite { cardinality: u64 },
    CountablyInfinite { presentation: SourceKey },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegularGaugeOrbit {
    continuous_orbit_dimension: GaugeContinuousDimension,
    discrete_orbit_cardinality: GaugeDiscreteOrbitCardinality,
}

impl RegularGaugeOrbit {
    #[must_use]
    pub const fn new(
        continuous_orbit_dimension: GaugeContinuousDimension,
        discrete_orbit_cardinality: GaugeDiscreteOrbitCardinality,
    ) -> Self {
        Self {
            continuous_orbit_dimension,
            discrete_orbit_cardinality,
        }
    }

    #[must_use]
    pub const fn continuous_orbit_dimension(&self) -> &GaugeContinuousDimension {
        &self.continuous_orbit_dimension
    }

    #[must_use]
    pub const fn discrete_orbit_cardinality(&self) -> &GaugeDiscreteOrbitCardinality {
        &self.discrete_orbit_cardinality
    }
}

/// Effective orbit geometry, kept orthogonal to the acting algebra. A
/// stratified action retains typed principal invariants plus an exact
/// orbit-type/stabilizer profile for singular strata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeOrbitGeometry {
    Regular {
        principal: RegularGaugeOrbit,
        stabilizer_profile: Option<SourceKey>,
    },
    /// Exact orbit-type/stabilizer profile across singular strata.
    Stratified {
        principal: RegularGaugeOrbit,
        orbit_type_stabilizer_profile: SourceKey,
    },
}

fn principal_gauge_orbit(geometry: &GaugeOrbitGeometry) -> &RegularGaugeOrbit {
    match geometry {
        GaugeOrbitGeometry::Regular { principal, .. }
        | GaugeOrbitGeometry::Stratified { principal, .. } => principal,
    }
}

fn gauge_algebra_orbit_compatible(
    algebra: &GaugeAlgebra,
    orbit: &RegularGaugeOrbit,
    require_effective: bool,
) -> bool {
    let algebra_is_nontrivial = match algebra {
        GaugeAlgebra::Continuous { group_dimension } => !group_dimension.is_zero(),
        GaugeAlgebra::Discrete { size } => match size {
            GaugeDiscreteSize::Finite { order } => *order > 1,
            GaugeDiscreteSize::CountablyInfinite { .. } => true,
        },
        GaugeAlgebra::Mixed {
            continuous_group_dimension,
            component_group,
        } => {
            !continuous_group_dimension.is_zero()
                && match component_group {
                    GaugeDiscreteSize::Finite { order } => *order > 1,
                    GaugeDiscreteSize::CountablyInfinite { .. } => true,
                }
        }
    };
    let orbit_matches_algebra = match algebra {
        GaugeAlgebra::Continuous { group_dimension } => {
            continuous_orbit_dimension_compatible(
                group_dimension,
                &orbit.continuous_orbit_dimension,
            ) && matches!(
                &orbit.discrete_orbit_cardinality,
                GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 }
            )
        }
        GaugeAlgebra::Discrete { size } => {
            orbit.continuous_orbit_dimension.is_zero()
                && match (size, &orbit.discrete_orbit_cardinality) {
                    (
                        GaugeDiscreteSize::Finite { order },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality },
                    ) => *cardinality >= 1 && *order % *cardinality == 0,
                    (
                        GaugeDiscreteSize::CountablyInfinite { .. },
                        GaugeDiscreteOrbitCardinality::Finite { cardinality },
                    ) => *cardinality >= 1,
                    (
                        GaugeDiscreteSize::CountablyInfinite { .. },
                        GaugeDiscreteOrbitCardinality::CountablyInfinite { .. },
                    ) => true,
                    (
                        GaugeDiscreteSize::Finite { .. },
                        GaugeDiscreteOrbitCardinality::CountablyInfinite { .. },
                    ) => false,
                }
        }
        GaugeAlgebra::Mixed {
            continuous_group_dimension,
            component_group,
        } => {
            continuous_orbit_dimension_compatible(
                continuous_group_dimension,
                &orbit.continuous_orbit_dimension,
            ) && match (component_group, &orbit.discrete_orbit_cardinality) {
                (
                    GaugeDiscreteSize::Finite { order },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality },
                ) => *cardinality >= 1 && *order % *cardinality == 0,
                (
                    GaugeDiscreteSize::CountablyInfinite { .. },
                    GaugeDiscreteOrbitCardinality::Finite { cardinality },
                ) => *cardinality >= 1,
                (
                    GaugeDiscreteSize::CountablyInfinite { .. },
                    GaugeDiscreteOrbitCardinality::CountablyInfinite { .. },
                ) => true,
                (
                    GaugeDiscreteSize::Finite { .. },
                    GaugeDiscreteOrbitCardinality::CountablyInfinite { .. },
                ) => false,
            }
        }
    };
    algebra_is_nontrivial
        && orbit_matches_algebra
        && (!require_effective
            || !orbit.continuous_orbit_dimension.is_zero()
            || !matches!(
                &orbit.discrete_orbit_cardinality,
                GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 }
            ))
}

fn continuous_orbit_dimension_compatible(
    group: &GaugeContinuousDimension,
    orbit: &GaugeContinuousDimension,
) -> bool {
    match (group, orbit) {
        (
            GaugeContinuousDimension::Finite {
                dimension: group_dimension,
            },
            GaugeContinuousDimension::Finite {
                dimension: orbit_dimension,
            },
        ) => orbit_dimension <= group_dimension,
        (
            GaugeContinuousDimension::InfiniteDimensional { .. },
            GaugeContinuousDimension::Finite { .. }
            | GaugeContinuousDimension::InfiniteDimensional { .. },
        ) => true,
        (
            GaugeContinuousDimension::Finite { .. },
            GaugeContinuousDimension::InfiniteDimensional { .. },
        ) => false,
    }
}

fn infinite_dimensional_profile_is_explicit(
    algebra: &GaugeAlgebra,
    geometry: &GaugeOrbitGeometry,
) -> bool {
    let group_is_infinite = matches!(
        algebra,
        GaugeAlgebra::Continuous {
            group_dimension: GaugeContinuousDimension::InfiniteDimensional { .. }
        } | GaugeAlgebra::Mixed {
            continuous_group_dimension: GaugeContinuousDimension::InfiniteDimensional { .. },
            ..
        }
    );
    let orbit_is_infinite = matches!(
        &principal_gauge_orbit(geometry).continuous_orbit_dimension,
        GaugeContinuousDimension::InfiniteDimensional { .. }
    );
    if !(group_is_infinite || orbit_is_infinite) {
        return true;
    }
    match geometry {
        GaugeOrbitGeometry::Regular {
            stabilizer_profile: Some(_),
            ..
        }
        | GaugeOrbitGeometry::Stratified { .. } => true,
        GaugeOrbitGeometry::Regular {
            stabilizer_profile: None,
            ..
        } => false,
    }
}

fn regular_orbit_support_compatible(
    geometry: &GaugeOrbitGeometry,
    support: &GaugeExtentSupport,
) -> bool {
    let GaugeOrbitGeometry::Regular { principal, .. } = geometry else {
        // Singular/fixed/nonproper local behavior is authorized by the exact
        // orbit-type profile and cannot be inferred from principal data.
        return true;
    };
    let has_continuous_orbit = !principal.continuous_orbit_dimension.is_zero();
    let has_nontrivial_discrete_orbit = !matches!(
        &principal.discrete_orbit_cardinality,
        GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 }
    );
    // `Regular` is the proper/local-covering branch for a finite discrete
    // principal orbit. Nonproper accumulation or fixed-stratum behavior must
    // use `Stratified` with its exact orbit-type profile.
    let finite_discrete_only = principal.continuous_orbit_dimension.is_zero()
        && matches!(
            &principal.discrete_orbit_cardinality,
            GaugeDiscreteOrbitCardinality::Finite { .. }
        );
    !((has_continuous_orbit && support.local_obstruction_parameters.is_empty())
        || (finite_discrete_only && !support.local_obstruction_parameters.is_empty())
        || ((has_continuous_orbit || has_nontrivial_discrete_orbit)
            != !support.global_obstruction_parameters.is_empty()))
}

fn gauge_algebra_source_keys(algebra: &GaugeAlgebra) -> BTreeSet<SourceKey> {
    let mut keys = match algebra {
        GaugeAlgebra::Discrete {
            size: GaugeDiscreteSize::CountablyInfinite { presentation },
        }
        | GaugeAlgebra::Mixed {
            component_group: GaugeDiscreteSize::CountablyInfinite { presentation },
            ..
        } => BTreeSet::from([presentation.clone()]),
        GaugeAlgebra::Continuous { .. }
        | GaugeAlgebra::Discrete {
            size: GaugeDiscreteSize::Finite { .. },
        }
        | GaugeAlgebra::Mixed {
            component_group: GaugeDiscreteSize::Finite { .. },
            ..
        } => BTreeSet::new(),
    };
    let continuous = match algebra {
        GaugeAlgebra::Continuous { group_dimension } => Some(group_dimension),
        GaugeAlgebra::Mixed {
            continuous_group_dimension,
            ..
        } => Some(continuous_group_dimension),
        GaugeAlgebra::Discrete { .. } => None,
    };
    if let Some(GaugeContinuousDimension::InfiniteDimensional { model_space }) = continuous {
        keys.insert(model_space.clone());
    }
    keys
}

fn gauge_orbit_source_keys(geometry: &GaugeOrbitGeometry) -> BTreeSet<SourceKey> {
    let mut keys = BTreeSet::new();
    let principal = principal_gauge_orbit(geometry);
    if let GaugeContinuousDimension::InfiniteDimensional { model_space } =
        &principal.continuous_orbit_dimension
    {
        keys.insert(model_space.clone());
    }
    if let GaugeDiscreteOrbitCardinality::CountablyInfinite { presentation } =
        &principal.discrete_orbit_cardinality
    {
        keys.insert(presentation.clone());
    }
    match geometry {
        GaugeOrbitGeometry::Regular {
            stabilizer_profile: Some(profile),
            ..
        } => {
            keys.insert(profile.clone());
        }
        GaugeOrbitGeometry::Stratified {
            orbit_type_stabilizer_profile,
            ..
        } => {
            keys.insert(orbit_type_stabilizer_profile.clone());
        }
        GaugeOrbitGeometry::Regular {
            stabilizer_profile: None,
            ..
        } => {}
    }
    keys
}

/// Epistemic posture is separate from the mathematical action. Neither status
/// is a theorem token; candidate gauges may be tested/refuted by assessment,
/// while assumed gauges are part of the declared physical problem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeStatus {
    Candidate { rationale: SourceKey },
    Assumed { assumption: SourceKey },
}

/// Exact information regime in which a declared gauge remains invariant.
/// Posterior persistence names the prior source in the problem registry; claim
/// admission later requires the claim's complete [`SourceRef`] to equal that
/// registry entry, so a same-key but different prior cannot inherit the gauge.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum GaugeInformationRegime {
    StructuralExactModel,
    ExactInputOutputMap,
    NoisyFiniteData,
    PosteriorUnderDeclaredPrior { joint_prior: SourceKey },
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum GaugeScalarDomain {
    Real,
    Complex { extension: SourceKey },
    MixedDiscreteContinuous { stratification: SourceKey },
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum GaugeLocus {
    WholeDomain,
    Stratum { definition: SourceKey },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GaugeProbabilityThreshold(u64);

impl GaugeProbabilityThreshold {
    pub fn try_new(probability: f64) -> Result<Self, IdentifiabilityError> {
        if !probability.is_finite() || probability <= 0.0 || probability > 1.0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "gauge probability threshold",
                detail: "probability must lie in (0,1]".to_string(),
            });
        }
        Ok(Self(canonical_f64(probability).to_bits()))
    }

    #[must_use]
    pub fn value(self) -> f64 {
        f64::from_bits(self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum GaugeQuantifierScope {
    AtRealization {
        realization: SourceKey,
    },
    AlmostEverywhere {
        measure: SourceKey,
    },
    ForAll {
        domain: SourceKey,
    },
    ProbabilityAtLeast {
        probability: GaugeProbabilityThreshold,
        measure: SourceKey,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct GaugeApplicabilityAxes {
    information: GaugeInformationRegime,
    scalar_domain: GaugeScalarDomain,
    locus: GaugeLocus,
    quantifier: GaugeQuantifierScope,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// Parameters whose projection can obstruct uniqueness in one exact claim
/// cell. These are not pointwise orbit-movement sets: a finite action can have
/// a singleton orbit at a fixed point yet destroy local injectivity nearby.
/// Local germ obstruction is always a subset of global fiber obstruction.
pub struct GaugeExtentSupport {
    local_obstruction_parameters: BTreeSet<ParameterRoleId>,
    global_obstruction_parameters: BTreeSet<ParameterRoleId>,
}

impl GaugeExtentSupport {
    pub fn try_new(
        local_obstruction_parameters: BTreeSet<ParameterRoleId>,
        global_obstruction_parameters: BTreeSet<ParameterRoleId>,
    ) -> Result<Self, IdentifiabilityError> {
        if global_obstruction_parameters.len() > MAX_IDENTIFIABILITY_ITEMS
            || !local_obstruction_parameters.is_subset(&global_obstruction_parameters)
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge extent support",
                detail: "local moved support must be a bounded subset of global moved support"
                    .to_string(),
            });
        }
        Ok(Self {
            local_obstruction_parameters,
            global_obstruction_parameters,
        })
    }

    #[must_use]
    pub const fn local_obstruction_parameters(&self) -> &BTreeSet<ParameterRoleId> {
        &self.local_obstruction_parameters
    }

    #[must_use]
    pub const fn global_obstruction_parameters(&self) -> &BTreeSet<ParameterRoleId> {
        &self.global_obstruction_parameters
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeCellDomain {
    case_obstruction_support: BTreeMap<CaseId, GaugeExtentSupport>,
}

impl GaugeCellDomain {
    pub fn try_new(
        case_obstruction_support: BTreeMap<CaseId, GaugeExtentSupport>,
    ) -> Result<Self, IdentifiabilityError> {
        if case_obstruction_support.is_empty()
            || case_obstruction_support.len() > MAX_IDENTIFIABILITY_ITEMS
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge applicability cell",
                detail: "a gauge cell needs bounded nonempty case-indexed moved support"
                    .to_string(),
            });
        }
        Ok(Self {
            case_obstruction_support,
        })
    }

    #[must_use]
    pub const fn case_obstruction_support(&self) -> &BTreeMap<CaseId, GaugeExtentSupport> {
        &self.case_obstruction_support
    }
}

impl GaugeApplicabilityAxes {
    #[must_use]
    pub const fn new(
        information: GaugeInformationRegime,
        scalar_domain: GaugeScalarDomain,
        locus: GaugeLocus,
        quantifier: GaugeQuantifierScope,
    ) -> Self {
        Self {
            information,
            scalar_domain,
            locus,
            quantifier,
        }
    }

    #[must_use]
    pub const fn information(&self) -> &GaugeInformationRegime {
        &self.information
    }

    #[must_use]
    pub const fn scalar_domain(&self) -> &GaugeScalarDomain {
        &self.scalar_domain
    }

    #[must_use]
    pub const fn locus(&self) -> &GaugeLocus {
        &self.locus
    }

    #[must_use]
    pub const fn quantifier(&self) -> &GaugeQuantifierScope {
        &self.quantifier
    }
}

/// Exact non-Cartesian validity domain of a gauge declaration.
///
/// Each cell keys information, scalar, locus, and quantifier semantics, then
/// records case-specific local/global obstruction support. The scope is
/// structural metadata bound into problem identity; neither case membership
/// nor a [`CasePurpose`] label proves symmetry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeValidityScope {
    cells: BTreeMap<GaugeApplicabilityAxes, GaugeCellDomain>,
}

impl GaugeValidityScope {
    pub fn try_new(
        cells: BTreeMap<GaugeApplicabilityAxes, GaugeCellDomain>,
    ) -> Result<Self, IdentifiabilityError> {
        if cells.is_empty()
            || cells.len() > MAX_IDENTIFIABILITY_ITEMS
            || cells
                .values()
                .any(|cell| cell.case_obstruction_support.is_empty())
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge validity scope",
                detail: "gauge validity needs bounded nonempty information-regime/case cells"
                    .to_string(),
            });
        }
        Ok(Self { cells })
    }

    #[must_use]
    pub const fn cells(&self) -> &BTreeMap<GaugeApplicabilityAxes, GaugeCellDomain> {
        &self.cells
    }
}

fn gauge_applicability_source_keys(axes: &GaugeApplicabilityAxes) -> BTreeSet<SourceKey> {
    let mut keys = BTreeSet::new();
    if let GaugeInformationRegime::PosteriorUnderDeclaredPrior { joint_prior } = &axes.information {
        keys.insert(joint_prior.clone());
    }
    match &axes.scalar_domain {
        GaugeScalarDomain::Real => {}
        GaugeScalarDomain::Complex { extension } => {
            keys.insert(extension.clone());
        }
        GaugeScalarDomain::MixedDiscreteContinuous { stratification } => {
            keys.insert(stratification.clone());
        }
    }
    if let GaugeLocus::Stratum { definition } = &axes.locus {
        keys.insert(definition.clone());
    }
    match &axes.quantifier {
        GaugeQuantifierScope::AtRealization { realization } => {
            keys.insert(realization.clone());
        }
        GaugeQuantifierScope::AlmostEverywhere { measure }
        | GaugeQuantifierScope::ProbabilityAtLeast { measure, .. } => {
            keys.insert(measure.clone());
        }
        GaugeQuantifierScope::ForAll { domain } => {
            keys.insert(domain.clone());
        }
    }
    keys
}

/// Declared physical gauge action. Computational quotient/slice choices live
/// in [`GaugeReductionPlan`] and therefore do not perturb physical ProblemId.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeDeclaration {
    id: GaugeClassId,
    members: BTreeSet<ParameterRoleId>,
    action: SourceKey,
    algebra: GaugeAlgebra,
    orbit_geometry: GaugeOrbitGeometry,
    status: GaugeStatus,
    validity: GaugeValidityScope,
}

impl GaugeDeclaration {
    pub fn try_new(
        id: GaugeClassId,
        members: BTreeSet<ParameterRoleId>,
        action: SourceKey,
        algebra: GaugeAlgebra,
        orbit_geometry: GaugeOrbitGeometry,
        status: GaugeStatus,
        validity: GaugeValidityScope,
    ) -> Result<Self, IdentifiabilityError> {
        if members.is_empty() || members.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::InvalidGauge {
                gauge: id,
                detail: "a gauge needs at least one bounded member".to_string(),
            });
        }
        if !gauge_algebra_orbit_compatible(
            &algebra,
            principal_gauge_orbit(&orbit_geometry),
            matches!(&status, GaugeStatus::Assumed { .. })
                && matches!(&orbit_geometry, GaugeOrbitGeometry::Regular { .. }),
        ) || !infinite_dimensional_profile_is_explicit(&algebra, &orbit_geometry)
        {
            return Err(IdentifiabilityError::InvalidGauge {
                gauge: id,
                detail: "principal orbit geometry must respect group/member bounds and finite component-group orbit-index divisibility; an assumed gauge must have a nontrivial effective principal orbit"
                    .to_string(),
            });
        }
        Ok(Self {
            id,
            members,
            action,
            algebra,
            orbit_geometry,
            status,
            validity,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &GaugeClassId {
        &self.id
    }

    #[must_use]
    pub const fn members(&self) -> &BTreeSet<ParameterRoleId> {
        &self.members
    }

    #[must_use]
    pub const fn action(&self) -> &SourceKey {
        &self.action
    }

    #[must_use]
    pub const fn validity(&self) -> &GaugeValidityScope {
        &self.validity
    }

    #[must_use]
    pub const fn algebra(&self) -> &GaugeAlgebra {
        &self.algebra
    }

    #[must_use]
    pub const fn orbit_geometry(&self) -> &GaugeOrbitGeometry {
        &self.orbit_geometry
    }

    #[must_use]
    pub const fn status(&self) -> &GaugeStatus {
        &self.status
    }
}

/// Exact hyperedge declaration for simultaneously active, overlapping gauge
/// actions. `law` names the action order/commutation/composition semantics;
/// pairwise declarations cannot substitute for a connected 3+-action
/// component. The effective geometry is that of the composed action.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeCompositionKind {
    /// Exact authority asserts a commuting/independent direct-product action.
    IndependentProduct,
    /// Exact authority defines generated, ordered, semidirect, bracket, or
    /// otherwise interacting composition semantics.
    Generated,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeCompositionDeclaration {
    id: GaugeCompositionId,
    members: BTreeSet<GaugeClassId>,
    kind: GaugeCompositionKind,
    law: SourceKey,
    effective_algebra: GaugeAlgebra,
    effective_orbit_geometry: GaugeOrbitGeometry,
    status: GaugeStatus,
    validity: GaugeValidityScope,
}

impl GaugeCompositionDeclaration {
    pub fn try_new(
        id: GaugeCompositionId,
        members: BTreeSet<GaugeClassId>,
        kind: GaugeCompositionKind,
        law: SourceKey,
        effective_algebra: GaugeAlgebra,
        effective_orbit_geometry: GaugeOrbitGeometry,
        status: GaugeStatus,
        validity: GaugeValidityScope,
    ) -> Result<Self, IdentifiabilityError> {
        if members.len() < 2 || members.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge composition members",
                detail: "a gauge composition needs at least two bounded members".to_string(),
            });
        }
        if !infinite_dimensional_profile_is_explicit(&effective_algebra, &effective_orbit_geometry)
        {
            return Err(IdentifiabilityError::InvalidText {
                field: "infinite-dimensional gauge composition profile",
                detail: "infinite-dimensional generated action needs an exact stabilizer/model-space/submersion profile"
                    .to_string(),
            });
        }
        Ok(Self {
            id,
            members,
            kind,
            law,
            effective_algebra,
            effective_orbit_geometry,
            status,
            validity,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &GaugeCompositionId {
        &self.id
    }

    #[must_use]
    pub const fn members(&self) -> &BTreeSet<GaugeClassId> {
        &self.members
    }

    #[must_use]
    pub const fn kind(&self) -> &GaugeCompositionKind {
        &self.kind
    }

    #[must_use]
    pub const fn law(&self) -> &SourceKey {
        &self.law
    }

    #[must_use]
    pub const fn effective_algebra(&self) -> &GaugeAlgebra {
        &self.effective_algebra
    }

    #[must_use]
    pub const fn effective_orbit_geometry(&self) -> &GaugeOrbitGeometry {
        &self.effective_orbit_geometry
    }

    #[must_use]
    pub const fn status(&self) -> &GaugeStatus {
        &self.status
    }

    #[must_use]
    pub const fn validity(&self) -> &GaugeValidityScope {
        &self.validity
    }
}

fn validate_independent_product_invariants(
    composition: &GaugeCompositionDeclaration,
    gauges: &BTreeMap<GaugeClassId, GaugeDeclaration>,
) -> Result<(), IdentifiabilityError> {
    if !matches!(&composition.kind, GaugeCompositionKind::IndependentProduct) {
        return Ok(());
    }

    let mut finite_group_dimension = 0_u64;
    let mut infinite_group_dimension = false;
    let mut finite_group_order = 1_u64;
    let mut infinite_component_group = false;
    let mut has_continuous_group = false;
    let mut has_discrete_group = false;

    let mut finite_orbit_dimension = 0_u64;
    let mut infinite_orbit_dimension = false;
    let mut finite_orbit_cardinality = 1_u64;
    let mut infinite_discrete_orbit = false;
    let mut any_stratified_orbit = false;

    for member in &composition.members {
        let gauge = gauges
            .get(member)
            .expect("composition member existence is checked before product invariants");
        let (continuous, discrete) = match &gauge.algebra {
            GaugeAlgebra::Continuous { group_dimension } => {
                has_continuous_group = true;
                (Some(group_dimension), None)
            }
            GaugeAlgebra::Discrete { size } => {
                has_discrete_group = true;
                (None, Some(size))
            }
            GaugeAlgebra::Mixed {
                continuous_group_dimension,
                component_group,
            } => {
                has_continuous_group = true;
                has_discrete_group = true;
                (Some(continuous_group_dimension), Some(component_group))
            }
        };
        if let Some(continuous) = continuous {
            match continuous {
                GaugeContinuousDimension::Finite { dimension } => {
                    finite_group_dimension = finite_group_dimension.checked_add(*dimension).ok_or(
                        IdentifiabilityError::Cardinality {
                            field: "independent-product group dimension",
                            detail: "finite direct-product dimension exceeds u64".to_string(),
                        },
                    )?;
                }
                GaugeContinuousDimension::InfiniteDimensional { .. } => {
                    infinite_group_dimension = true;
                }
            }
        }
        if let Some(discrete) = discrete {
            match discrete {
                GaugeDiscreteSize::Finite { order } => {
                    finite_group_order = finite_group_order.checked_mul(*order).ok_or(
                        IdentifiabilityError::Cardinality {
                            field: "independent-product group order",
                            detail: "finite direct-product order exceeds u64".to_string(),
                        },
                    )?;
                }
                GaugeDiscreteSize::CountablyInfinite { .. } => {
                    infinite_component_group = true;
                }
            }
        }

        let principal = principal_gauge_orbit(&gauge.orbit_geometry);
        match &principal.continuous_orbit_dimension {
            GaugeContinuousDimension::Finite { dimension } => {
                finite_orbit_dimension = finite_orbit_dimension.checked_add(*dimension).ok_or(
                    IdentifiabilityError::Cardinality {
                        field: "independent-product orbit dimension",
                        detail: "finite direct-product orbit dimension exceeds u64".to_string(),
                    },
                )?;
            }
            GaugeContinuousDimension::InfiniteDimensional { .. } => {
                infinite_orbit_dimension = true;
            }
        }
        match &principal.discrete_orbit_cardinality {
            GaugeDiscreteOrbitCardinality::Finite { cardinality } => {
                finite_orbit_cardinality = finite_orbit_cardinality
                    .checked_mul(*cardinality)
                    .ok_or(IdentifiabilityError::Cardinality {
                        field: "independent-product orbit cardinality",
                        detail: "finite direct-product orbit cardinality exceeds u64".to_string(),
                    })?;
            }
            GaugeDiscreteOrbitCardinality::CountablyInfinite { .. } => {
                infinite_discrete_orbit = true;
            }
        }
        any_stratified_orbit |=
            matches!(&gauge.orbit_geometry, GaugeOrbitGeometry::Stratified { .. });
    }

    let effective_continuous = match &composition.effective_algebra {
        GaugeAlgebra::Continuous { group_dimension }
        | GaugeAlgebra::Mixed {
            continuous_group_dimension: group_dimension,
            ..
        } => Some(group_dimension),
        GaugeAlgebra::Discrete { .. } => None,
    };
    let effective_discrete = match &composition.effective_algebra {
        GaugeAlgebra::Discrete { size }
        | GaugeAlgebra::Mixed {
            component_group: size,
            ..
        } => Some(size),
        GaugeAlgebra::Continuous { .. } => None,
    };
    let continuous_matches = if has_continuous_group {
        effective_continuous.is_some_and(|dimension| {
            if infinite_group_dimension {
                matches!(dimension, GaugeContinuousDimension::InfiniteDimensional { .. })
            } else {
                matches!(dimension, GaugeContinuousDimension::Finite { dimension } if *dimension == finite_group_dimension)
            }
        })
    } else {
        effective_continuous.is_none()
    };
    let discrete_matches = if has_discrete_group {
        effective_discrete.is_some_and(|size| {
            if infinite_component_group {
                matches!(size, GaugeDiscreteSize::CountablyInfinite { .. })
            } else {
                matches!(size, GaugeDiscreteSize::Finite { order } if *order == finite_group_order)
            }
        })
    } else {
        effective_discrete.is_none()
    };

    let effective_principal = principal_gauge_orbit(&composition.effective_orbit_geometry);
    let orbit_dimension_matches = if infinite_orbit_dimension {
        matches!(
            &effective_principal.continuous_orbit_dimension,
            GaugeContinuousDimension::InfiniteDimensional { .. }
        )
    } else {
        matches!(
            &effective_principal.continuous_orbit_dimension,
            GaugeContinuousDimension::Finite { dimension } if *dimension == finite_orbit_dimension
        )
    };
    let orbit_cardinality_matches = if infinite_discrete_orbit {
        matches!(
            &effective_principal.discrete_orbit_cardinality,
            GaugeDiscreteOrbitCardinality::CountablyInfinite { .. }
        )
    } else {
        matches!(
            &effective_principal.discrete_orbit_cardinality,
            GaugeDiscreteOrbitCardinality::Finite { cardinality } if *cardinality == finite_orbit_cardinality
        )
    };
    let stratification_matches = any_stratified_orbit
        == matches!(
            &composition.effective_orbit_geometry,
            GaugeOrbitGeometry::Stratified { .. }
        );

    if !(continuous_matches
        && discrete_matches
        && orbit_dimension_matches
        && orbit_cardinality_matches
        && stratification_matches)
    {
        return Err(IdentifiabilityError::InvalidText {
            field: "independent-product gauge invariants",
            detail: format!(
                "composition {} must exactly add finite continuous dimensions, multiply finite component orders/orbit cardinalities with checked arithmetic, preserve infinite categories, and be stratified iff any factor is stratified",
                composition.id
            ),
        });
    }
    Ok(())
}

/// Expected codimension of an execution-time gauge-fixing slice. Fixed
/// codimension is compared with effective orbit dimension, never group
/// dimension. Singular actions retain an exact profile source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeSliceCodimension {
    FixedFinite {
        codimension: u64,
    },
    FixedInfinite {
        codimension_model: SourceRef,
        compatibility: SourceRef,
    },
    Stratified {
        profile: SourceRef,
    },
}

/// Computational slice choice, deliberately separate from physical joint
/// admissibility constraints and ProblemId.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeSlicePlan {
    support: BTreeSet<ParameterRoleId>,
    constraint: SourceRef,
    expected_codimension: GaugeSliceCodimension,
    coverage: SourceRef,
}

impl GaugeSlicePlan {
    pub fn try_new(
        support: BTreeSet<ParameterRoleId>,
        constraint: SourceRef,
        expected_codimension: GaugeSliceCodimension,
        coverage: SourceRef,
    ) -> Result<Self, IdentifiabilityError> {
        if support.is_empty() || support.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge slice support",
                detail: "a gauge slice needs bounded nonempty parameter support".to_string(),
            });
        }
        Ok(Self {
            support,
            constraint,
            expected_codimension,
            coverage,
        })
    }

    #[must_use]
    pub const fn support(&self) -> &BTreeSet<ParameterRoleId> {
        &self.support
    }

    #[must_use]
    pub const fn constraint(&self) -> &SourceRef {
        &self.constraint
    }

    #[must_use]
    pub const fn expected_codimension(&self) -> &GaugeSliceCodimension {
        &self.expected_codimension
    }

    #[must_use]
    pub const fn coverage(&self) -> &SourceRef {
        &self.coverage
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeQuotientPlan {
    RegularAtlas {
        quotient_map: SourceRef,
        local_section_atlas: SourceRef,
        coverage: SourceRef,
    },
    SingularOrGeneralized {
        quotient_map: SourceRef,
        quotient_profile: SourceRef,
        local_models: Option<SourceRef>,
    },
    InvariantMap {
        invariants: SourceRef,
        completeness_profile: SourceRef,
    },
    GroupoidOrStack {
        presentation: SourceRef,
        quotient_profile: SourceRef,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContinuousGaugeReductionPlan {
    Quotient { quotient: GaugeQuotientPlan },
    Slice { slice: GaugeSlicePlan },
}

/// Exact relation from predecessor reductions to one later stage. This is
/// execution semantics, never inferred from carrier disjointness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeReductionStageRelation {
    NormalSubgroupTower {
        normality: SourceRef,
        induced_residual_action: SourceRef,
    },
    SemidirectOrGenerated {
        extension: SourceRef,
        induced_residual_action: SourceRef,
    },
    TransverseSlices {
        transversality: SourceRef,
    },
    GaugeForGauge {
        reducibility: SourceRef,
        induced_residual_action: SourceRef,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeReductionStage {
    Root,
    After {
        predecessors: BTreeSet<GaugeReductionId>,
        composition_law: SourceRef,
        relation: GaugeReductionStageRelation,
    },
}

/// Measure semantics for a quotient/slice applied to posterior or marginalized
/// coordinates. A generic algorithm name is not a change-of-variables proof.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeMeasureSemantics {
    NotApplicable {
        reason: String,
    },
    Pushforward {
        source_measure: SourceRef,
        reduced_measure: SourceRef,
        transport: SourceRef,
        jacobian_or_disintegration: SourceRef,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeReductionPlan {
    Unreduced {
        reason: String,
    },
    Quotient {
        quotient: GaugeQuotientPlan,
    },
    Slice {
        slice: GaugeSlicePlan,
    },
    ContinuousReductionWithDiscreteResidual {
        reduction: ContinuousGaugeReductionPlan,
        normal_subgroup: SourceRef,
        factor_extension: SourceRef,
        residual_quotient_action: SourceRef,
        compatibility: SourceRef,
    },
}

/// One execution-time reduction over exact preregistered claim cells. The
/// claim IDs bind information/scalar/locus/quantifier/case axes without
/// duplicating physical applicability metadata into execution identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeReductionBinding {
    id: GaugeReductionId,
    action: GaugeActionReference,
    claims: BTreeSet<ClaimId>,
    plan: GaugeReductionPlan,
    stage: GaugeReductionStage,
    measure: GaugeMeasureSemantics,
}

impl GaugeReductionBinding {
    pub fn try_new(
        id: GaugeReductionId,
        action: GaugeActionReference,
        claims: BTreeSet<ClaimId>,
        plan: GaugeReductionPlan,
        stage: GaugeReductionStage,
        measure: GaugeMeasureSemantics,
    ) -> Result<Self, IdentifiabilityError> {
        if claims.is_empty() || claims.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge reduction claim cells",
                detail: "a reduction must bind bounded nonempty exact claim IDs".to_string(),
            });
        }
        if let GaugeReductionPlan::Unreduced { reason } = &plan {
            validate_reason(reason, "unreduced gauge reason")?;
        }
        match &stage {
            GaugeReductionStage::Root => {}
            GaugeReductionStage::After {
                predecessors,
                composition_law,
                ..
            } => {
                if predecessors.is_empty() || predecessors.len() > MAX_IDENTIFIABILITY_ITEMS {
                    return Err(IdentifiabilityError::Cardinality {
                        field: "gauge reduction stage predecessors",
                        detail: "a non-root stage needs bounded nonempty predecessors".to_string(),
                    });
                }
                if predecessors.contains(&id) {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "gauge reduction stage",
                        detail: "a reduction cannot depend on itself".to_string(),
                    });
                }
                if composition_law.kind != SourceKind::GaugeComposition {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "gauge reduction stage composition law",
                        detail: "staged reductions need exact GaugeComposition ordering semantics"
                            .to_string(),
                    });
                }
            }
        }
        if let GaugeMeasureSemantics::NotApplicable { reason } = &measure {
            validate_reason(reason, "gauge reduction measure semantics")?;
        }
        Ok(Self {
            id,
            action,
            claims,
            plan,
            stage,
            measure,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &GaugeReductionId {
        &self.id
    }

    #[must_use]
    pub const fn action(&self) -> &GaugeActionReference {
        &self.action
    }

    #[must_use]
    pub const fn claims(&self) -> &BTreeSet<ClaimId> {
        &self.claims
    }

    #[must_use]
    pub const fn plan(&self) -> &GaugeReductionPlan {
        &self.plan
    }

    #[must_use]
    pub const fn stage(&self) -> &GaugeReductionStage {
        &self.stage
    }

    #[must_use]
    pub const fn measure(&self) -> &GaugeMeasureSemantics {
        &self.measure
    }
}

fn gauge_quotient_sources(
    quotient: &GaugeQuotientPlan,
) -> Result<Vec<&SourceRef>, IdentifiabilityError> {
    let sources = match quotient {
        GaugeQuotientPlan::RegularAtlas {
            quotient_map,
            local_section_atlas,
            coverage,
        } => {
            if quotient_map.kind != SourceKind::GaugeQuotientMap
                || local_section_atlas.kind != SourceKind::GaugeSection
                || coverage.kind != SourceKind::GaugeSection
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "regular gauge quotient plan",
                    detail: "regular atlas needs exact quotient-map, section-atlas, and coverage semantics"
                        .to_string(),
                });
            }
            vec![quotient_map, local_section_atlas, coverage]
        }
        GaugeQuotientPlan::SingularOrGeneralized {
            quotient_map,
            quotient_profile,
            local_models,
        } => {
            if quotient_map.kind != SourceKind::GaugeQuotientMap
                || quotient_profile.kind != SourceKind::GaugeQuotientProfile
                || local_models
                    .as_ref()
                    .is_some_and(|source| source.kind != SourceKind::FunctionalModelSpace)
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "singular/generalized gauge quotient plan",
                    detail: "generalized quotient needs exact map/profile and typed optional local models"
                        .to_string(),
                });
            }
            let mut result = vec![quotient_map, quotient_profile];
            result.extend(local_models.iter());
            result
        }
        GaugeQuotientPlan::InvariantMap {
            invariants,
            completeness_profile,
        } => {
            if invariants.kind != SourceKind::GaugeInvariantMap
                || completeness_profile.kind != SourceKind::GaugeQuotientProfile
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge invariant-map quotient plan",
                    detail: "invariant map and completeness profile need exact dedicated semantics"
                        .to_string(),
                });
            }
            vec![invariants, completeness_profile]
        }
        GaugeQuotientPlan::GroupoidOrStack {
            presentation,
            quotient_profile,
        } => {
            if presentation.kind != SourceKind::GaugeGroupoidPresentation
                || quotient_profile.kind != SourceKind::GaugeQuotientProfile
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge groupoid/stack quotient plan",
                    detail: "groupoid/stack plan needs exact presentation and quotient profile"
                        .to_string(),
                });
            }
            vec![presentation, quotient_profile]
        }
    };
    Ok(sources)
}

fn gauge_reduction_sources(
    plan: &GaugeReductionPlan,
) -> Result<Vec<&SourceRef>, IdentifiabilityError> {
    fn slice_sources(slice: &GaugeSlicePlan) -> Vec<&SourceRef> {
        let mut sources = vec![&slice.constraint, &slice.coverage];
        match &slice.expected_codimension {
            GaugeSliceCodimension::FixedFinite { .. } => {}
            GaugeSliceCodimension::FixedInfinite {
                codimension_model,
                compatibility,
            } => {
                sources.push(codimension_model);
                sources.push(compatibility);
            }
            GaugeSliceCodimension::Stratified { profile } => sources.push(profile),
        }
        sources
    }
    let mut sources = Vec::new();
    match plan {
        GaugeReductionPlan::Unreduced { .. } => {}
        GaugeReductionPlan::Quotient { quotient } => {
            sources.extend(gauge_quotient_sources(quotient)?);
        }
        GaugeReductionPlan::Slice { slice } => sources.extend(slice_sources(slice)),
        GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction,
            normal_subgroup,
            factor_extension,
            residual_quotient_action,
            compatibility,
        } => {
            match reduction {
                ContinuousGaugeReductionPlan::Quotient { quotient } => {
                    sources.extend(gauge_quotient_sources(quotient)?);
                }
                ContinuousGaugeReductionPlan::Slice { slice } => {
                    sources.extend(slice_sources(slice));
                }
            }
            if normal_subgroup.kind != SourceKind::GaugeSubgroupCertificate
                || factor_extension.kind != SourceKind::GaugeReductionLaw
                || residual_quotient_action.kind != SourceKind::GaugeResidualAction
                || compatibility.kind != SourceKind::GaugeReductionLaw
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "continuous reduction with discrete residual",
                    detail: "mixed reduction needs exact normal-subgroup, factor-extension, residual-quotient-action, and compatibility semantics"
                        .to_string(),
                });
            }
            sources.extend([
                normal_subgroup,
                factor_extension,
                residual_quotient_action,
                compatibility,
            ]);
        }
    }
    Ok(sources)
}

fn gauge_reduction_stage_sources(
    stage: &GaugeReductionStage,
) -> Result<Vec<&SourceRef>, IdentifiabilityError> {
    let GaugeReductionStage::After {
        composition_law,
        relation,
        ..
    } = stage
    else {
        return Ok(Vec::new());
    };
    if composition_law.kind != SourceKind::GaugeComposition {
        return Err(IdentifiabilityError::InvalidText {
            field: "gauge reduction stage law",
            detail: "stage composition law needs GaugeComposition semantics".to_string(),
        });
    }
    let mut sources = vec![composition_law];
    match relation {
        GaugeReductionStageRelation::NormalSubgroupTower {
            normality,
            induced_residual_action,
        } => {
            if normality.kind != SourceKind::GaugeSubgroupCertificate
                || induced_residual_action.kind != SourceKind::GaugeResidualAction
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "normal-subgroup reduction stage",
                    detail:
                        "normal tower needs exact subgroup certificate and induced residual action"
                            .to_string(),
                });
            }
            sources.extend([normality, induced_residual_action]);
        }
        GaugeReductionStageRelation::SemidirectOrGenerated {
            extension,
            induced_residual_action,
        } => {
            if extension.kind != SourceKind::GaugeReductionLaw
                || induced_residual_action.kind != SourceKind::GaugeResidualAction
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "generated reduction stage",
                    detail: "generated stage needs exact extension and induced residual action"
                        .to_string(),
                });
            }
            sources.extend([extension, induced_residual_action]);
        }
        GaugeReductionStageRelation::TransverseSlices { transversality } => {
            if transversality.kind != SourceKind::GaugeReductionLaw {
                return Err(IdentifiabilityError::InvalidText {
                    field: "transverse-slice reduction stage",
                    detail: "transverse slices need an exact transversality law".to_string(),
                });
            }
            sources.push(transversality);
        }
        GaugeReductionStageRelation::GaugeForGauge {
            reducibility,
            induced_residual_action,
        } => {
            if reducibility.kind != SourceKind::GaugeReductionLaw
                || induced_residual_action.kind != SourceKind::GaugeResidualAction
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge-for-gauge reduction stage",
                    detail: "gauge-for-gauge stages need exact reducibility and induced-action semantics"
                        .to_string(),
                });
            }
            sources.extend([reducibility, induced_residual_action]);
        }
    }
    Ok(sources)
}

fn gauge_measure_sources(
    measure: &GaugeMeasureSemantics,
) -> Result<Vec<&SourceRef>, IdentifiabilityError> {
    match measure {
        GaugeMeasureSemantics::NotApplicable { .. } => Ok(Vec::new()),
        GaugeMeasureSemantics::Pushforward {
            source_measure,
            reduced_measure,
            transport,
            jacobian_or_disintegration,
        } => {
            if source_measure.kind != SourceKind::ProbabilityMeasure
                || reduced_measure.kind != SourceKind::ProbabilityMeasure
                || transport.kind != SourceKind::GaugeMeasureTransport
                || jacobian_or_disintegration.kind != SourceKind::GaugeMeasureTransport
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge reduction measure transport",
                    detail: "measure-aware reduction needs source/reduced probability measures plus exact pushforward and Jacobian/disintegration semantics"
                        .to_string(),
                });
            }
            Ok(vec![
                source_measure,
                reduced_measure,
                transport,
                jacobian_or_disintegration,
            ])
        }
    }
}

fn gauge_reduction_precedes(
    predecessor: &GaugeReductionId,
    successor: &GaugeReductionId,
    reductions: &BTreeMap<GaugeReductionId, GaugeReductionBinding>,
) -> bool {
    let mut pending = vec![successor.clone()];
    let mut seen = BTreeSet::new();
    while let Some(current) = pending.pop() {
        if !seen.insert(current.clone()) {
            continue;
        }
        let Some(binding) = reductions.get(&current) else {
            continue;
        };
        if let GaugeReductionStage::After { predecessors, .. } = &binding.stage {
            if predecessors.contains(predecessor) {
                return true;
            }
            pending.extend(predecessors.iter().cloned());
        }
    }
    false
}

fn validate_gauge_reduction_dag(
    reductions: &BTreeMap<GaugeReductionId, GaugeReductionBinding>,
) -> Result<(), IdentifiabilityError> {
    let mut incoming = reductions
        .keys()
        .cloned()
        .map(|id| (id, 0_usize))
        .collect::<BTreeMap<_, _>>();
    let mut dependents = BTreeMap::<GaugeReductionId, BTreeSet<GaugeReductionId>>::new();
    for (id, binding) in reductions {
        let GaugeReductionStage::After { predecessors, .. } = &binding.stage else {
            continue;
        };
        if matches!(&binding.plan, GaugeReductionPlan::Unreduced { .. }) {
            return Err(IdentifiabilityError::InvalidText {
                field: "unreduced gauge stage",
                detail: format!("unreduced binding {id} cannot be a later reduction stage"),
            });
        }
        incoming.insert(id.clone(), predecessors.len());
        for predecessor in predecessors {
            let predecessor_binding = reductions.get(predecessor).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "gauge reduction predecessor",
                    id: predecessor.to_string(),
                }
            })?;
            if matches!(
                &predecessor_binding.plan,
                GaugeReductionPlan::Unreduced { .. }
            ) || !binding.claims.is_subset(&predecessor_binding.claims)
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge reduction predecessor coverage",
                    detail: format!(
                        "stage {id} requires a reduced predecessor {predecessor} covering every successor claim cell"
                    ),
                });
            }
            dependents
                .entry(predecessor.clone())
                .or_default()
                .insert(id.clone());
        }
    }
    let mut ready = incoming
        .iter()
        .filter_map(|(id, degree)| (*degree == 0).then_some(id.clone()))
        .collect::<BTreeSet<_>>();
    let mut visited = 0_usize;
    while let Some(id) = ready.pop_first() {
        visited += 1;
        if let Some(children) = dependents.get(&id) {
            for child in children {
                let degree = incoming.get_mut(child).expect("known reduction child");
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(child.clone());
                }
            }
        }
    }
    if visited != reductions.len() {
        return Err(IdentifiabilityError::InvalidText {
            field: "gauge reduction stage graph",
            detail: "gauge reduction stages must form an acyclic dependency graph".to_string(),
        });
    }

    let claim_ids = reductions
        .values()
        .flat_map(|binding| binding.claims.iter().cloned())
        .collect::<BTreeSet<_>>();
    for claim in claim_ids {
        let active = reductions
            .values()
            .filter(|binding| {
                binding.claims.contains(&claim)
                    && !matches!(&binding.plan, GaugeReductionPlan::Unreduced { .. })
            })
            .map(|binding| binding.id.clone())
            .collect::<Vec<_>>();
        for left_index in 0..active.len() {
            for right_index in (left_index + 1)..active.len() {
                let left = &active[left_index];
                let right = &active[right_index];
                if !gauge_reduction_precedes(left, right, reductions)
                    && !gauge_reduction_precedes(right, left, reductions)
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "unordered gauge reductions",
                        detail: format!(
                            "claim {claim} has reductions {left} and {right} without an exact staged law; disjoint carriers do not imply commutation—use one declared IndependentProduct action or an ordered reduction stage"
                        ),
                    });
                }
            }
        }
    }
    Ok(())
}

fn reduction_uses_regular_atlas(plan: &GaugeReductionPlan) -> bool {
    matches!(
        plan,
        GaugeReductionPlan::Quotient {
            quotient: GaugeQuotientPlan::RegularAtlas { .. }
        } | GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction: ContinuousGaugeReductionPlan::Quotient {
                quotient: GaugeQuotientPlan::RegularAtlas { .. }
            },
            ..
        }
    )
}

/// Explicit sharing group for cases that intentionally reuse observations or
/// raw experiment sources under one joint likelihood.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataSharingGroup {
    cases: BTreeSet<CaseId>,
    joint_likelihood: SourceKey,
    justification: String,
}

impl DataSharingGroup {
    pub fn try_new(
        cases: BTreeSet<CaseId>,
        joint_likelihood: SourceKey,
        justification: impl Into<String>,
    ) -> Result<Self, IdentifiabilityError> {
        let justification = justification.into();
        if cases.len() < 2 || cases.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "data sharing group",
                detail: "a sharing group needs at least two bounded cases".to_string(),
            });
        }
        validate_reason(&justification, "data sharing justification")?;
        Ok(Self {
            cases,
            joint_likelihood,
            justification,
        })
    }

    #[must_use]
    pub const fn cases(&self) -> &BTreeSet<CaseId> {
        &self.cases
    }

    #[must_use]
    pub const fn joint_likelihood(&self) -> &SourceKey {
        &self.joint_likelihood
    }

    #[must_use]
    pub fn justification(&self) -> &str {
        &self.justification
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DataReusePolicy {
    Disjoint,
    Shared { groups: Vec<DataSharingGroup> },
}

fn retrospective_experiment(case: &StudyCaseDocument) -> Option<&SourceKey> {
    match &case.data {
        CaseDataDeclaration::Retrospective { experiment, .. } => Some(experiment),
        CaseDataDeclaration::Prospective => None,
    }
}

fn sharing_group_membership(policy: &DataReusePolicy, case: &CaseId) -> Option<usize> {
    match policy {
        DataReusePolicy::Disjoint => None,
        DataReusePolicy::Shared { groups } => {
            groups.iter().position(|group| group.cases.contains(case))
        }
    }
}

fn normalize_joint_noise(noise: JointNoiseModel) -> Result<JointNoiseModel, IdentifiabilityError> {
    let JointNoiseModel::DenseCorrelation {
        order,
        correlation,
        model,
    } = noise
    else {
        return Ok(noise);
    };
    let positions = order
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, key)| (key, index))
        .collect::<BTreeMap<_, _>>();
    if positions.len() != order.len() || correlation.dimension() != order.len() {
        return Err(IdentifiabilityError::Covariance {
            detail: "dense-correlation order must be unique and match matrix dimension".to_string(),
        });
    }
    let sorted = positions.keys().cloned().collect::<Vec<_>>();
    let mut lower = Vec::with_capacity(correlation.lower_triangle().len());
    for row in 0..sorted.len() {
        for column in 0..=row {
            lower.push(matrix_get(
                &correlation,
                positions[&sorted[row]],
                positions[&sorted[column]],
            ));
        }
    }
    let correlation = CovarianceMatrix::try_new(sorted.len(), lower).map_err(|error| {
        IdentifiabilityError::Vv {
            detail: error.to_string(),
        }
    })?;
    Ok(JointNoiseModel::DenseCorrelation {
        order: sorted,
        correlation,
        model,
    })
}

fn problem_source_reachability(
    context_source: &SourceKey,
    material_source: &SourceKey,
    model_source: &SourceKey,
    graph_source: &SourceKey,
    joint_prior: Option<&SourceKey>,
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
    constraints: &BTreeMap<ConstraintId, JointConstraint>,
    admissible_domain: &AdmissibleDomainWitness,
    cases: &BTreeMap<CaseId, StudyCaseDocument>,
    influences: &BTreeMap<InfluenceId, InfluenceDeclaration>,
    gauges: &BTreeMap<GaugeClassId, GaugeDeclaration>,
    gauge_compositions: &BTreeMap<GaugeCompositionId, GaugeCompositionDeclaration>,
    joint_noise: &JointNoiseModel,
    data_reuse: &DataReusePolicy,
) -> BTreeSet<SourceKey> {
    let mut used = BTreeSet::from([
        context_source.clone(),
        material_source.clone(),
        model_source.clone(),
        graph_source.clone(),
    ]);
    used.extend(joint_prior.cloned());
    for parameter in parameters.values() {
        match &parameter.treatment {
            ParameterTreatment::Conditioned(value) => {
                used.insert(value.source.clone());
            }
            ParameterTreatment::Derived { definition, .. } => {
                used.insert(definition.clone());
            }
            _ => {}
        }
        match &parameter.owner {
            ParameterOwnerBinding::ConstitutiveModel => {}
            ParameterOwnerBinding::InitialState { state_path } => {
                used.insert(state_path.clone());
            }
            ParameterOwnerBinding::Instrument { metrology, .. } => {
                used.insert(metrology.clone());
            }
            ParameterOwnerBinding::Discrepancy { family } => {
                used.insert(family.clone());
            }
            ParameterOwnerBinding::ControlledInput { protocol } => {
                used.insert(protocol.clone());
            }
            ParameterOwnerBinding::Population { hierarchy } => {
                used.insert(hierarchy.clone());
            }
        }
        match &parameter.scope {
            ParameterScopeBinding::Field { support, .. } => {
                used.insert(support.clone());
            }
            ParameterScopeBinding::Hierarchical { hierarchy, .. } => {
                used.insert(hierarchy.clone());
            }
            _ => {}
        }
    }
    for constraint in constraints.values() {
        match &constraint.kind {
            JointConstraintKind::ExternalManifold {
                definition,
                codimension,
                ..
            } => {
                used.insert(definition.clone());
                if let ConstraintCodimension::InfiniteDimensional { profile } = codimension {
                    used.insert(profile.clone());
                }
            }
            JointConstraintKind::StochasticCoupling { distribution, .. } => {
                used.insert(distribution.clone());
            }
            _ => {}
        }
    }
    if let Some(claim) = &admissible_domain.opaque_membership_claim {
        used.insert(claim.source.clone());
    }
    for case in cases.values() {
        used.insert(case.forward_model.clone());
        used.extend([
            case.physics_sources.frame_transform.clone(),
            case.physics_sources.specimen_geometry.clone(),
            case.physics_sources.specimen_process.clone(),
            case.physics_sources.specimen_preparation.clone(),
            case.physics_sources.load_path.clone(),
            case.physics_sources.environment_path.clone(),
            case.physics_sources.time_grid.clone(),
        ]);
        if let Some(initial_state) = &case.physics_sources.initial_state {
            used.insert(initial_state.clone());
        }
        if let CaseDataDeclaration::Retrospective {
            experiment,
            split,
            parser,
            preprocessing,
            ..
        } = &case.data
        {
            used.extend([
                experiment.clone(),
                split.clone(),
                parser.clone(),
                preprocessing.clone(),
            ]);
        }
        for observation in case.observations.values() {
            used.extend([
                observation.operator.clone(),
                observation.aggregation.clone(),
                observation.sensor.clone(),
                observation.unit_definition.clone(),
            ]);
            if let MarginalNoiseSpec::Empirical {
                distribution,
                finite_variance_model,
                ..
            } = &observation.noise
            {
                used.extend([distribution.clone(), finite_variance_model.clone()]);
            }
            match &observation.missingness {
                MissingnessAssumption::Complete { assumption } => {
                    used.insert(assumption.clone());
                }
                MissingnessAssumption::Modeled { mechanism } => {
                    used.insert(mechanism.clone());
                }
                MissingnessAssumption::Unknown { .. } => {}
            }
        }
        for discrepancy in case.discrepancies.values() {
            match discrepancy {
                StudyDiscrepancy::NotApplicable { basis } => match basis {
                    DiscrepancyInapplicability::PhysicalApplicability { assumption }
                    | DiscrepancyInapplicability::ProspectiveDesign { assumption } => {
                        used.insert(assumption.clone());
                    }
                    DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                        generator,
                        producer: _,
                        production_binding,
                        assumption,
                    } => {
                        used.extend([
                            generator.clone(),
                            production_binding.clone(),
                            assumption.clone(),
                        ]);
                    }
                },
                StudyDiscrepancy::AssumedZero { assumption } => {
                    used.insert(assumption.clone());
                }
                StudyDiscrepancy::Modeled {
                    family,
                    support,
                    confounding_guard,
                    ..
                } => {
                    used.extend([family.clone(), support.clone(), confounding_guard.clone()]);
                }
                StudyDiscrepancy::Uncharacterized { .. } => {}
            }
        }
        for group in &case.observation_sharing {
            used.insert(group.joint_likelihood.clone());
        }
    }
    for influence in influences.values() {
        match &influence.representation {
            InfluenceRepresentation::StateMediated { state_path } => {
                used.insert(state_path.clone());
            }
            InfluenceRepresentation::Composite { operator, .. } => {
                used.insert(operator.clone());
            }
            InfluenceRepresentation::ExternalDefinition { definition } => {
                used.insert(definition.clone());
            }
            InfluenceRepresentation::Direct => {}
        }
    }
    for gauge in gauges.values() {
        used.insert(gauge.action.clone());
        used.insert(match &gauge.status {
            GaugeStatus::Candidate { rationale } => rationale.clone(),
            GaugeStatus::Assumed { assumption } => assumption.clone(),
        });
        used.extend(gauge_algebra_source_keys(&gauge.algebra));
        used.extend(gauge_orbit_source_keys(&gauge.orbit_geometry));
        for axes in gauge.validity.cells.keys() {
            used.extend(gauge_applicability_source_keys(axes));
        }
    }
    for composition in gauge_compositions.values() {
        used.insert(composition.law.clone());
        used.insert(match &composition.status {
            GaugeStatus::Candidate { rationale } => rationale.clone(),
            GaugeStatus::Assumed { assumption } => assumption.clone(),
        });
        used.extend(gauge_algebra_source_keys(&composition.effective_algebra));
        used.extend(gauge_orbit_source_keys(
            &composition.effective_orbit_geometry,
        ));
        for axes in composition.validity.cells.keys() {
            used.extend(gauge_applicability_source_keys(axes));
        }
    }
    match joint_noise {
        JointNoiseModel::Independent { assumption } => {
            used.insert(assumption.clone());
        }
        JointNoiseModel::DenseCorrelation { model, .. }
        | JointNoiseModel::ExternalKernel { model } => {
            used.insert(model.clone());
        }
        _ => {}
    }
    if let DataReusePolicy::Shared { groups } = data_reuse {
        for group in groups {
            used.insert(group.joint_likelihood.clone());
        }
    }
    used
}

/// Canonical unresolved physical/statistical question.  No coordinate,
/// tolerance, algorithm, random seed, build fingerprint, or result receipt is
/// permitted in this type.
#[derive(Clone, PartialEq)]
pub struct IdentifiabilityProblemDocument {
    schema_version: u32,
    context_source: SourceKey,
    material_source: SourceKey,
    model_source: SourceKey,
    graph_source: SourceKey,
    joint_prior: Option<SourceKey>,
    sources: BTreeMap<SourceKey, SourceRef>,
    parameters: BTreeMap<ParameterRoleId, StudyParameter>,
    constraints: BTreeMap<ConstraintId, JointConstraint>,
    admissible_domain: AdmissibleDomainWitness,
    cases: BTreeMap<CaseId, StudyCaseDocument>,
    influences: BTreeMap<InfluenceId, InfluenceDeclaration>,
    gauges: BTreeMap<GaugeClassId, GaugeDeclaration>,
    gauge_compositions: BTreeMap<GaugeCompositionId, GaugeCompositionDeclaration>,
    joint_noise: JointNoiseModel,
    data_reuse: DataReusePolicy,
}

impl fmt::Debug for IdentifiabilityProblemDocument {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prospective_cases = self
            .cases
            .values()
            .filter(|case| matches!(&case.data, CaseDataDeclaration::Prospective))
            .count();
        formatter
            .debug_struct("IdentifiabilityProblemDocument")
            .field("schema_version", &self.schema_version)
            .field("source_count", &self.sources.len())
            .field("parameter_count", &self.parameters.len())
            .field("constraint_count", &self.constraints.len())
            .field("case_count", &self.cases.len())
            .field("prospective_case_count", &prospective_cases)
            .field(
                "retrospective_case_count",
                &self.cases.len().saturating_sub(prospective_cases),
            )
            .field("influence_count", &self.influences.len())
            .field("gauge_count", &self.gauges.len())
            .field("gauge_composition_count", &self.gauge_compositions.len())
            .finish_non_exhaustive()
    }
}

/// Aggregate bound across nested collections. Per-vector caps alone permit a
/// Cartesian explosion of individually legal 4096-entry maps.
pub const MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS: usize = 65_536;

#[derive(Default)]
struct StructuralItemBudget {
    used: usize,
}

impl StructuralItemBudget {
    fn add(&mut self, count: usize, field: &'static str) -> Result<(), IdentifiabilityError> {
        self.used =
            self.used
                .checked_add(count)
                .ok_or_else(|| IdentifiabilityError::Cardinality {
                    field,
                    detail: "aggregate structural item count overflow".to_string(),
                })?;
        if self.used > MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field,
                detail: format!(
                    "aggregate structural item budget {} exceeds {}",
                    self.used, MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS
                ),
            });
        }
        Ok(())
    }
}

fn validate_problem_structural_budget(
    document: &IdentifiabilityProblemDocument,
) -> Result<(), IdentifiabilityError> {
    let mut budget = StructuralItemBudget::default();
    budget.add(
        4 + usize::from(document.joint_prior.is_some()),
        "problem roots",
    )?;
    budget.add(document.sources.len(), "problem sources")?;
    budget.add(document.parameters.len(), "problem parameters")?;
    for parameter in document.parameters.values() {
        if let ParameterTreatment::Derived { parents, .. } = &parameter.treatment {
            budget.add(parents.len(), "derived parameter parents")?;
        }
        let scoped_cases = match &parameter.scope {
            ParameterScopeBinding::Cases(cases)
            | ParameterScopeBinding::MaterialLot { cases, .. }
            | ParameterScopeBinding::Field { cases, .. }
            | ParameterScopeBinding::Hierarchical { cases, .. } => cases.len(),
            ParameterScopeBinding::Global | ParameterScopeBinding::Specimen { .. } => 0,
        };
        budget.add(scoped_cases, "parameter case scope")?;
    }
    budget.add(document.constraints.len(), "problem constraints")?;
    for constraint in document.constraints.values() {
        let members = match &constraint.kind {
            JointConstraintKind::Affine { terms, .. } => terms.len(),
            JointConstraintKind::Simplex { members, .. }
            | JointConstraintKind::ExternalManifold { members, .. }
            | JointConstraintKind::StochasticCoupling { members, .. } => members.len(),
            JointConstraintKind::Ordered { members, .. } => members.len(),
        };
        budget.add(members, "constraint members")?;
    }
    budget.add(
        document.admissible_domain.values.len(),
        "admissible-domain witness",
    )?;
    budget.add(document.cases.len(), "problem cases")?;
    for case in document.cases.values() {
        budget.add(case.observations.len(), "case observations")?;
        budget.add(case.discrepancies.len(), "case discrepancies")?;
        budget.add(case.observation_sharing.len(), "observation sharing groups")?;
        for observation in case.observations.values() {
            if let ObservationRows::Retrospective(rows) = &observation.rows {
                budget.add(rows.len(), "observation rows")?;
            }
        }
        for discrepancy in case.discrepancies.values() {
            if let StudyDiscrepancy::Modeled { parameters, .. } = discrepancy {
                budget.add(parameters.len(), "modeled discrepancy parameters")?;
            }
        }
        for group in &case.observation_sharing {
            budget.add(group.channels.len(), "observation-sharing channels")?;
            budget.add(group.rows.len(), "observation-sharing rows")?;
        }
    }
    budget.add(document.influences.len(), "problem influences")?;
    for influence in document.influences.values() {
        if let InfluenceRepresentation::Composite { inputs, .. } = &influence.representation {
            budget.add(inputs.len(), "composite influence inputs")?;
        }
    }
    budget.add(document.gauges.len(), "problem gauges")?;
    for gauge in document.gauges.values() {
        budget.add(gauge.members.len(), "gauge members")?;
        budget.add(gauge.validity.cells.len(), "gauge validity cells")?;
        for cell in gauge.validity.cells.values() {
            budget.add(cell.case_obstruction_support.len(), "gauge validity cases")?;
            for support in cell.case_obstruction_support.values() {
                budget.add(
                    support.local_obstruction_parameters.len(),
                    "local gauge support",
                )?;
                budget.add(
                    support.global_obstruction_parameters.len(),
                    "global gauge support",
                )?;
            }
        }
    }
    budget.add(document.gauge_compositions.len(), "gauge compositions")?;
    for composition in document.gauge_compositions.values() {
        budget.add(composition.members.len(), "gauge composition members")?;
        budget.add(
            composition.validity.cells.len(),
            "gauge composition validity cells",
        )?;
        for cell in composition.validity.cells.values() {
            budget.add(
                cell.case_obstruction_support.len(),
                "gauge composition cases",
            )?;
            for support in cell.case_obstruction_support.values() {
                budget.add(
                    support.local_obstruction_parameters.len(),
                    "local composition support",
                )?;
                budget.add(
                    support.global_obstruction_parameters.len(),
                    "global composition support",
                )?;
            }
        }
    }
    if let JointNoiseModel::DenseCorrelation {
        order, correlation, ..
    } = &document.joint_noise
    {
        budget.add(order.len(), "joint-noise order")?;
        budget.add(correlation.lower_triangle().len(), "joint-noise entries")?;
    }
    if let DataReusePolicy::Shared { groups } = &document.data_reuse {
        budget.add(groups.len(), "data-sharing groups")?;
        for group in groups {
            budget.add(group.cases.len(), "data-sharing cases")?;
        }
    }
    Ok(())
}

fn require_source<'a>(
    sources: &'a BTreeMap<SourceKey, SourceRef>,
    key: &SourceKey,
    field: &'static str,
) -> Result<&'a SourceRef, IdentifiabilityError> {
    sources
        .get(key)
        .ok_or_else(|| IdentifiabilityError::UnknownReference {
            field,
            id: key.to_string(),
        })
}

fn require_source_kind(
    sources: &BTreeMap<SourceKey, SourceRef>,
    key: &SourceKey,
    expected: SourceKind,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    let source = require_source(sources, key, field)?;
    if source.kind != expected {
        return Err(IdentifiabilityError::InvalidText {
            field,
            detail: format!(
                "source {} has kind {:?}, expected {:?}",
                key, source.kind, expected
            ),
        });
    }
    Ok(())
}

fn require_source_kind_in(
    sources: &BTreeMap<SourceKey, SourceRef>,
    key: &SourceKey,
    allowed: &[SourceKind],
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    let source = require_source(sources, key, field)?;
    if !allowed.contains(&source.kind) {
        return Err(IdentifiabilityError::InvalidText {
            field,
            detail: format!(
                "source {} has kind {:?}, expected one of {allowed:?}",
                source.key, source.kind
            ),
        });
    }
    Ok(())
}

fn validate_gauge_applicability_sources(
    axes: &GaugeApplicabilityAxes,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<(), IdentifiabilityError> {
    if let GaugeInformationRegime::PosteriorUnderDeclaredPrior { joint_prior } = &axes.information {
        require_source_kind(
            sources,
            joint_prior,
            SourceKind::ProbabilityMeasure,
            "gauge posterior prior",
        )?;
    }
    match &axes.scalar_domain {
        GaugeScalarDomain::Real => {}
        GaugeScalarDomain::Complex { extension } => require_source_kind(
            sources,
            extension,
            SourceKind::AlgebraicExtension,
            "gauge complex scalar extension",
        )?,
        GaugeScalarDomain::MixedDiscreteContinuous { stratification } => require_source_kind(
            sources,
            stratification,
            SourceKind::Stratification,
            "gauge mixed scalar stratification",
        )?,
    }
    if let GaugeLocus::Stratum { definition } = &axes.locus {
        require_source_kind(
            sources,
            definition,
            SourceKind::Stratification,
            "gauge locus stratum",
        )?;
    }
    match &axes.quantifier {
        GaugeQuantifierScope::AtRealization { realization } => require_source_kind(
            sources,
            realization,
            SourceKind::QuantifierRealization,
            "gauge realization scope",
        )?,
        GaugeQuantifierScope::AlmostEverywhere { measure } => require_source_kind_in(
            sources,
            measure,
            &[SourceKind::ReferenceMeasure, SourceKind::ProbabilityMeasure],
            "gauge almost-everywhere measure",
        )?,
        GaugeQuantifierScope::ForAll { domain } => require_source_kind(
            sources,
            domain,
            SourceKind::QuantifierDomain,
            "gauge universal domain",
        )?,
        GaugeQuantifierScope::ProbabilityAtLeast { measure, .. } => require_source_kind(
            sources,
            measure,
            SourceKind::ProbabilityMeasure,
            "gauge probability measure",
        )?,
    }
    Ok(())
}

fn validate_gauge_algebra_orbit_sources(
    algebra: &GaugeAlgebra,
    geometry: &GaugeOrbitGeometry,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<(), IdentifiabilityError> {
    for key in gauge_algebra_source_keys(algebra) {
        require_source_kind(
            sources,
            &key,
            SourceKind::GaugeGroupPresentation,
            "gauge group presentation",
        )?;
    }
    let principal = principal_gauge_orbit(geometry);
    if let GaugeContinuousDimension::InfiniteDimensional { model_space } =
        &principal.continuous_orbit_dimension
    {
        require_source_kind(
            sources,
            model_space,
            SourceKind::GaugeOrbitPresentation,
            "gauge infinite-dimensional orbit model space",
        )?;
    }
    if let GaugeDiscreteOrbitCardinality::CountablyInfinite { presentation } =
        &principal.discrete_orbit_cardinality
    {
        require_source_kind(
            sources,
            presentation,
            SourceKind::GaugeOrbitPresentation,
            "gauge orbit presentation",
        )?;
    }
    let profile = match geometry {
        GaugeOrbitGeometry::Regular {
            stabilizer_profile, ..
        } => stabilizer_profile.as_ref(),
        GaugeOrbitGeometry::Stratified {
            orbit_type_stabilizer_profile,
            ..
        } => Some(orbit_type_stabilizer_profile),
    };
    if let Some(profile) = profile {
        require_source_kind(
            sources,
            profile,
            SourceKind::GaugeOrbitTypeProfile,
            "gauge orbit-type/stabilizer profile",
        )?;
    }
    Ok(())
}

fn require_case_physics_source(
    sources: &BTreeMap<SourceKey, SourceRef>,
    key: &SourceKey,
    expected_kind: SourceKind,
    expected_hash: ContentHash,
    expected_domain: &'static str,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    let source = require_source(sources, key, field)?;
    if source.kind != expected_kind
        || source.expected_hash != expected_hash
        || source.content_hash_domain != expected_domain
        || source.contract_version != CASE_PHYSICS_SOURCE_CONTRACT_VERSION
    {
        return Err(IdentifiabilityError::SourceMismatch { field });
    }
    Ok(())
}

fn insert_unique<K: Ord + Clone + fmt::Display, V>(
    rows: Vec<V>,
    field: &'static str,
    key_of: impl Fn(&V) -> &K,
) -> Result<BTreeMap<K, V>, IdentifiabilityError> {
    if rows.is_empty() || rows.len() > MAX_IDENTIFIABILITY_ITEMS {
        return Err(IdentifiabilityError::Cardinality {
            field,
            detail: "collection must be bounded and nonempty".to_string(),
        });
    }
    let mut result = BTreeMap::new();
    for row in rows {
        let key = key_of(&row).clone();
        if result.insert(key.clone(), row).is_some() {
            return Err(IdentifiabilityError::Duplicate {
                field,
                id: key.to_string(),
            });
        }
    }
    Ok(result)
}

fn validate_source_key(
    sources: &BTreeMap<SourceKey, SourceRef>,
    key: &SourceKey,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    require_source(sources, key, field).map(|_| ())
}

fn validate_derived_parameter_dag(
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
) -> Result<(), IdentifiabilityError> {
    let mut incoming = parameters
        .keys()
        .cloned()
        .map(|role| (role, 0_usize))
        .collect::<BTreeMap<_, _>>();
    let mut children = BTreeMap::<ParameterRoleId, BTreeSet<ParameterRoleId>>::new();
    for (role, parameter) in parameters {
        if let ParameterTreatment::Derived { parents, .. } = &parameter.treatment {
            for parent in parents {
                if !parameters.contains_key(parent) {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "derived parameter parent",
                        id: parent.to_string(),
                    });
                }
                children
                    .entry(parent.clone())
                    .or_default()
                    .insert(role.clone());
            }
            incoming.insert(role.clone(), parents.len());
        }
    }
    let mut ready = incoming
        .iter()
        .filter_map(|(role, degree)| (*degree == 0).then_some(role.clone()))
        .collect::<BTreeSet<_>>();
    let mut visited = 0_usize;
    while let Some(role) = ready.pop_first() {
        visited += 1;
        if let Some(dependents) = children.get(&role) {
            for dependent in dependents {
                let degree = incoming.get_mut(dependent).expect("known derived child");
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(dependent.clone());
                }
            }
        }
    }
    if visited != parameters.len() {
        let role = incoming
            .iter()
            .find_map(|(role, degree)| (*degree > 0).then_some(role))
            .expect("an unvisited node has positive indegree");
        return Err(IdentifiabilityError::InvalidNumeric {
            field: "derived parameter graph",
            detail: format!("cycle reaches parameter {role}"),
        });
    }
    Ok(())
}

fn observation_for<'a>(
    cases: &'a BTreeMap<CaseId, StudyCaseDocument>,
    key: &ObservationKey,
) -> Result<&'a StudyObservation, IdentifiabilityError> {
    cases
        .get(&key.case)
        .and_then(|case| case.observations.get(&key.channel))
        .ok_or_else(|| IdentifiabilityError::UnknownReference {
            field: "composite observation key",
            id: format!("{}:{}", key.case, key.channel),
        })
}

fn initial_state_schema_version(state: InitialStateBinding) -> u32 {
    match state {
        InitialStateBinding::Zero { schema_version }
        | InitialStateBinding::Explicit { schema_version, .. } => schema_version,
    }
}

fn functional_observations(functional: &DistributionFunctional) -> Vec<&ObservationKey> {
    match functional {
        DistributionFunctional::Location { observation }
        | DistributionFunctional::LogScale { observation }
        | DistributionFunctional::MissingnessLogit { observation }
        | DistributionFunctional::CensoringLogit { observation } => vec![observation],
        DistributionFunctional::Correlation { left, right } => vec![left, right],
    }
}

fn transitive_influence_ids(
    root: &InfluenceId,
    influences: &BTreeMap<InfluenceId, InfluenceDeclaration>,
) -> Result<BTreeSet<InfluenceId>, IdentifiabilityError> {
    let mut pending = vec![root.clone()];
    let mut closure = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !closure.insert(id.clone()) {
            continue;
        }
        let influence =
            influences
                .get(&id)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field: "transitive influence input",
                    id: id.to_string(),
                })?;
        if let InfluenceRepresentation::Composite { inputs, .. } = &influence.representation {
            pending.extend(inputs.iter().cloned());
        }
    }
    Ok(closure)
}

fn validate_joint_constraint(
    constraint: &JointConstraint,
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<(), IdentifiabilityError> {
    let require_members = |members: &BTreeSet<ParameterRoleId>, minimum: usize| {
        if members.len() < minimum || members.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "joint constraint members",
                detail: format!(
                    "this joint-constraint variant needs at least {minimum} bounded member{}",
                    if minimum == 1 { "" } else { "s" }
                ),
            });
        }
        for member in members {
            if !parameters.contains_key(member) {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "joint constraint member",
                    id: member.to_string(),
                });
            }
        }
        Ok(())
    };
    match &constraint.kind {
        JointConstraintKind::Affine {
            terms,
            relation,
            rhs_si,
            residual_quantity,
        } => {
            if terms.is_empty() || terms.len() > MAX_IDENTIFIABILITY_ITEMS || !rhs_si.is_finite() {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "affine joint constraint",
                    detail: "requires at least one bounded term and a finite RHS".to_string(),
                });
            }
            if matches!(relation, ConstraintRelation::Equal)
                && terms.iter().all(|term| term.coefficient == 0.0)
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "affine equality rank",
                    detail: "an all-zero affine equality has rank zero and cannot supply the declared codimension"
                        .to_string(),
                });
            }
            let mut seen = BTreeSet::new();
            let mut minimum = 0.0;
            let mut maximum = 0.0;
            for term in terms {
                let parameter = parameters.get(&term.parameter).ok_or_else(|| {
                    IdentifiabilityError::UnknownReference {
                        field: "affine constraint member",
                        id: term.parameter.to_string(),
                    }
                })?;
                if !seen.insert(term.parameter.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "affine constraint member",
                        id: term.parameter.to_string(),
                    });
                }
                let product =
                    checked_add_dims(parameter.quantity.dims(), term.coefficient_quantity.dims())
                        .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                        field: "affine constraint units",
                        detail: "dimension exponent overflow".to_string(),
                    })?;
                if product != residual_quantity.dims() {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "affine constraint units",
                        detail: format!(
                            "coefficient times {} does not have the residual dimensions",
                            term.parameter
                        ),
                    });
                }
                let endpoints = [
                    term.coefficient * parameter.domain.lo,
                    term.coefficient * parameter.domain.hi,
                ];
                minimum += endpoints[0].min(endpoints[1]);
                maximum += endpoints[0].max(endpoints[1]);
            }
            let feasible = match relation {
                ConstraintRelation::Equal => *rhs_si >= minimum && *rhs_si <= maximum,
                ConstraintRelation::LessOrEqual => minimum <= *rhs_si,
                ConstraintRelation::GreaterOrEqual => maximum >= *rhs_si,
            };
            if !minimum.is_finite() || !maximum.is_finite() || !feasible {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "affine constraint feasibility",
                    detail: "affine constraint has no witness in the Cartesian domain enclosure"
                        .to_string(),
                });
            }
        }
        JointConstraintKind::Simplex {
            members,
            total_si,
            quantity,
        } => {
            require_members(members, 1)?;
            if !total_si.is_finite()
                || members
                    .iter()
                    .any(|role| parameters[role].quantity != *quantity)
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "simplex constraint",
                    detail: "members require one exact quantity and a finite total".to_string(),
                });
            }
            let minimum = members
                .iter()
                .map(|role| parameters[role].domain.lo)
                .sum::<f64>();
            let maximum = members
                .iter()
                .map(|role| parameters[role].domain.hi)
                .sum::<f64>();
            if members.iter().any(|role| parameters[role].domain.lo < 0.0)
                || !minimum.is_finite()
                || !maximum.is_finite()
                || *total_si < minimum
                || *total_si > maximum
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "simplex constraint feasibility",
                    detail: "simplex members must be nonnegative and their total attainable"
                        .to_string(),
                });
            }
        }
        JointConstraintKind::Ordered { members, strict } => {
            let member_set = members.iter().cloned().collect::<BTreeSet<_>>();
            require_members(&member_set, 2)?;
            if member_set.len() != members.len() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "ordered constraint member",
                    id: constraint.id.to_string(),
                });
            }
            let first = parameters[&members[0]].quantity;
            if members
                .iter()
                .any(|role| parameters[role].quantity != first)
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "ordered constraint units",
                    detail: "ordered members need one exact quantity".to_string(),
                });
            }
            // Greedily track the infimum of the smallest feasible prefix.
            // For a strict chain, an interval whose upper endpoint equals the
            // prefix infimum has no room for the required positive separation;
            // any positive gap contains enough real points for this finite
            // chain, without inventing a machine-epsilon semantics.
            let mut prefix_infimum = parameters[&members[0]].domain.lo;
            for role in members.iter().skip(1) {
                let domain = parameters[role].domain;
                let feasible = if *strict {
                    prefix_infimum < domain.hi
                } else {
                    prefix_infimum <= domain.hi
                };
                if !feasible {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "ordered constraint feasibility",
                        detail: format!(
                            "ordered member {role} has no {} witness after the preceding domain prefix",
                            if *strict { "strict" } else { "non-strict" }
                        ),
                    });
                }
                prefix_infimum = prefix_infimum.max(domain.lo);
            }
        }
        JointConstraintKind::ExternalManifold {
            members,
            definition,
            codimension,
        } => {
            require_members(members, 1)?;
            require_source_kind(
                sources,
                definition,
                SourceKind::ExternalManifold,
                "external manifold",
            )?;
            match codimension {
                ConstraintCodimension::Finite { codimension: 0 } => {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "external manifold codimension",
                        detail: "finite codimension must be positive".to_string(),
                    });
                }
                ConstraintCodimension::Finite { codimension }
                    if usize::try_from(*codimension)
                        .map_or(true, |codimension| codimension > members.len()) =>
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "external manifold codimension",
                        detail:
                            "finite codimension cannot exceed the finite scalar carrier dimension"
                                .to_string(),
                    });
                }
                ConstraintCodimension::InfiniteDimensional { profile } => {
                    require_source_kind(
                        sources,
                        profile,
                        SourceKind::ExternalManifold,
                        "external manifold infinite-dimensional profile",
                    )?;
                    return Err(IdentifiabilityError::InvalidText {
                        field: "external manifold infinite-dimensional carrier",
                        detail: "infinite-dimensional codimension is reserved for the typed function-space carrier model; finite scalar parameter members cannot realize it"
                            .to_string(),
                    });
                }
                ConstraintCodimension::Finite { .. } => {}
            }
        }
        JointConstraintKind::StochasticCoupling {
            members,
            distribution,
        } => {
            require_members(members, 2)?;
            require_source_kind(
                sources,
                distribution,
                SourceKind::Prior,
                "joint distribution",
            )?;
        }
    }
    Ok(())
}

fn joint_constraint_support(constraint: &JointConstraint) -> BTreeSet<ParameterRoleId> {
    match &constraint.kind {
        JointConstraintKind::Affine { terms, .. } => {
            terms.iter().map(|term| term.parameter.clone()).collect()
        }
        JointConstraintKind::Simplex { members, .. }
        | JointConstraintKind::ExternalManifold { members, .. }
        | JointConstraintKind::StochasticCoupling { members, .. } => members.clone(),
        JointConstraintKind::Ordered { members, .. } => members.iter().cloned().collect(),
    }
}

/// Codimension that follows from an internally understood equality geometry.
/// Inequalities, order cones, and stochastic couplings deliberately return
/// `None`: a [`GaugeSlicePlan`] may use them only through its exact external
/// transversality/coverage declaration.
fn intrinsic_constraint_codimension(constraint: &JointConstraint) -> Option<u64> {
    match &constraint.kind {
        JointConstraintKind::Affine {
            relation: ConstraintRelation::Equal,
            ..
        }
        | JointConstraintKind::Simplex { .. } => Some(1),
        JointConstraintKind::ExternalManifold {
            codimension: ConstraintCodimension::Finite { codimension },
            ..
        } => Some(*codimension),
        JointConstraintKind::ExternalManifold {
            codimension: ConstraintCodimension::InfiniteDimensional { .. },
            ..
        } => None,
        JointConstraintKind::Affine { .. }
        | JointConstraintKind::Ordered { .. }
        | JointConstraintKind::StochasticCoupling { .. } => None,
    }
}

fn validate_gauge_slice(
    action: &GaugeActionReference,
    slice: &GaugeSlicePlan,
    problem: &IdentifiabilityProblemDocument,
) -> Result<(), IdentifiabilityError> {
    let (carrier, geometry) = match action {
        GaugeActionReference::Single(id) => {
            let gauge =
                problem
                    .gauges
                    .get(id)
                    .ok_or_else(|| IdentifiabilityError::UnknownReference {
                        field: "execution gauge reduction action",
                        id: id.to_string(),
                    })?;
            (gauge.members.clone(), &gauge.orbit_geometry)
        }
        GaugeActionReference::Product(id) | GaugeActionReference::Composition(id) => {
            let composition = problem.gauge_compositions.get(id).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "execution gauge reduction composition",
                    id: id.to_string(),
                }
            })?;
            let kind_matches = matches!(
                (action, &composition.kind),
                (
                    GaugeActionReference::Product(_),
                    GaugeCompositionKind::IndependentProduct
                ) | (
                    GaugeActionReference::Composition(_),
                    GaugeCompositionKind::Generated
                )
            );
            if !kind_matches {
                return Err(IdentifiabilityError::InvalidText {
                    field: "execution gauge reduction composition kind",
                    detail: "Product/Composition reference disagrees with the exact declaration"
                        .to_string(),
                });
            }
            (
                composition
                    .members
                    .iter()
                    .flat_map(|member| problem.gauges[member].members.iter().cloned())
                    .collect(),
                &composition.effective_orbit_geometry,
            )
        }
    };
    if !slice.support.is_subset(&carrier) {
        return Err(IdentifiabilityError::InvalidText {
            field: "execution gauge slice support",
            detail: "slice support must be a nonempty subset of the exact action carrier"
                .to_string(),
        });
    }
    if !matches!(
        slice.constraint.kind,
        SourceKind::Constraint | SourceKind::ExternalManifold
    ) || slice.coverage.kind != SourceKind::GaugeSection
    {
        return Err(IdentifiabilityError::InvalidText {
            field: "execution gauge slice sources",
            detail: "slice needs an exact constraint source and GaugeSection coverage source"
                .to_string(),
        });
    }
    match (&slice.expected_codimension, geometry) {
        (
            GaugeSliceCodimension::FixedFinite { codimension },
            GaugeOrbitGeometry::Regular { principal, .. },
        ) if matches!(
            &principal.continuous_orbit_dimension,
            GaugeContinuousDimension::Finite { dimension } if dimension == codimension
        ) => {}
        (
            GaugeSliceCodimension::FixedInfinite {
                codimension_model,
                compatibility,
            },
            GaugeOrbitGeometry::Regular { principal, .. },
        ) if matches!(
            &principal.continuous_orbit_dimension,
            GaugeContinuousDimension::InfiniteDimensional { .. }
                if codimension_model.kind == SourceKind::FunctionalModelSpace
                    && compatibility.kind == SourceKind::GaugeSection
        ) => {}
        (
            GaugeSliceCodimension::Stratified { profile },
            GaugeOrbitGeometry::Stratified {
                orbit_type_stabilizer_profile,
                ..
            },
        ) if profile.kind == SourceKind::GaugeOrbitTypeProfile
            && problem.sources.get(orbit_type_stabilizer_profile) == Some(profile) => {}
        _ => {
            return Err(IdentifiabilityError::InvalidText {
                field: "execution gauge slice codimension",
                detail: "slice codimension/profile does not match effective orbit geometry"
                    .to_string(),
            });
        }
    }
    Ok(())
}

fn parameter_membership_source_keys(parameter: &StudyParameter) -> BTreeSet<SourceKey> {
    let mut keys = BTreeSet::new();
    match &parameter.treatment {
        ParameterTreatment::Conditioned(value) => {
            keys.insert(value.source.clone());
        }
        ParameterTreatment::Derived { definition, .. } => {
            keys.insert(definition.clone());
        }
        ParameterTreatment::Estimated
        | ParameterTreatment::Profiled
        | ParameterTreatment::Marginalized => {}
    }
    match &parameter.owner {
        ParameterOwnerBinding::InitialState { state_path } => {
            keys.insert(state_path.clone());
        }
        ParameterOwnerBinding::Instrument { metrology, .. } => {
            keys.insert(metrology.clone());
        }
        ParameterOwnerBinding::Discrepancy { family } => {
            keys.insert(family.clone());
        }
        ParameterOwnerBinding::ControlledInput { protocol } => {
            keys.insert(protocol.clone());
        }
        ParameterOwnerBinding::Population { hierarchy } => {
            keys.insert(hierarchy.clone());
        }
        ParameterOwnerBinding::ConstitutiveModel => {}
    }
    match &parameter.scope {
        ParameterScopeBinding::Field { support, .. } => {
            keys.insert(support.clone());
        }
        ParameterScopeBinding::Hierarchical { hierarchy, .. } => {
            keys.insert(hierarchy.clone());
        }
        ParameterScopeBinding::Global
        | ParameterScopeBinding::Cases(_)
        | ParameterScopeBinding::MaterialLot { .. }
        | ParameterScopeBinding::Specimen { .. } => {}
    }
    keys
}

fn admissible_domain_witness_binding(
    witness: &AdmissibleDomainWitness,
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
    constraints: &BTreeMap<ConstraintId, JointConstraint>,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<ContentHash, IdentifiabilityError> {
    let mut writer = CanonicalWriter::new();
    writer.count(witness.values.len(), "admissible-domain witness values")?;
    for (role, value) in &witness.values {
        encode_role(&mut writer, role)?;
        writer.f64(*value);
    }
    writer.count(parameters.len(), "admissible-domain parameter semantics")?;
    for parameter in parameters.values() {
        encode_study_parameter(&mut writer, parameter)?;
    }
    let parameter_sources = parameters
        .values()
        .flat_map(parameter_membership_source_keys)
        .collect::<BTreeSet<_>>();
    writer.count(
        parameter_sources.len(),
        "admissible-domain parameter source bindings",
    )?;
    for key in parameter_sources {
        encode_source_ref(
            &mut writer,
            sources
                .get(&key)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field: "admissible-domain parameter source",
                    id: key.to_string(),
                })?,
        )?;
    }
    writer.count(
        constraints.len(),
        "admissible-domain constraint conjunction",
    )?;
    for constraint in constraints.values() {
        encode_constraint(&mut writer, constraint)?;
        let external_sources = match &constraint.kind {
            JointConstraintKind::ExternalManifold {
                definition,
                codimension,
                ..
            } => {
                let mut keys = vec![definition];
                if let ConstraintCodimension::InfiniteDimensional { profile } = codimension {
                    keys.push(profile);
                }
                keys
            }
            JointConstraintKind::StochasticCoupling { distribution, .. } => vec![distribution],
            _ => Vec::new(),
        };
        for source in external_sources {
            encode_source_ref(
                &mut writer,
                sources
                    .get(source)
                    .ok_or_else(|| IdentifiabilityError::UnknownReference {
                        field: "opaque admissible-domain constraint source",
                        id: source.to_string(),
                    })?,
            )?;
        }
    }
    Ok(hash_domain(
        ADMISSIBLE_DOMAIN_WITNESS_BINDING_DOMAIN,
        &writer.finish()?,
    ))
}

/// Build the canonical byte preimage of an admissible-domain membership
/// certificate from the exact witness, parameter semantics, full constraint
/// conjunction, and referenced opaque-constraint sources.
///
/// The returned bytes are a content-addressing envelope, not a scientific
/// proof or an issuer authentication. A certificate source must hash these
/// bytes under [`ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN`]; external trust
/// policy remains represented separately by source admission authority.
pub fn admissible_domain_membership_certificate_preimage(
    witness: &AdmissibleDomainWitness,
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
    constraints: &BTreeMap<ConstraintId, JointConstraint>,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<[u8; 32], IdentifiabilityError> {
    Ok(*admissible_domain_witness_binding(witness, parameters, constraints, sources)?.as_bytes())
}

fn validate_admissible_domain_witness(
    witness: &AdmissibleDomainWitness,
    parameters: &BTreeMap<ParameterRoleId, StudyParameter>,
    constraints: &BTreeMap<ConstraintId, JointConstraint>,
    sources: &BTreeMap<SourceKey, SourceRef>,
) -> Result<(), IdentifiabilityError> {
    let expected_roles = parameters.keys().cloned().collect::<BTreeSet<_>>();
    let witnessed_roles = witness.values.keys().cloned().collect::<BTreeSet<_>>();
    if witnessed_roles != expected_roles {
        let detail = expected_roles
            .difference(&witnessed_roles)
            .next()
            .map_or_else(
                || "witness contains an unknown parameter".to_string(),
                |missing| format!("witness omits parameter {missing}"),
            );
        return Err(IdentifiabilityError::InvalidNumeric {
            field: "admissible-domain witness coverage",
            detail,
        });
    }
    for (role, parameter) in parameters {
        let value = witness.values[role];
        if !value.is_finite() || value < parameter.domain.lo || value > parameter.domain.hi {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "admissible-domain witness bounds",
                detail: format!("witness value for {role} lies outside its parameter domain"),
            });
        }
        if let ParameterTreatment::Conditioned(conditioned) = &parameter.treatment
            && !same_f64(value, conditioned.value_si)
        {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "admissible-domain conditioned value",
                detail: format!(
                    "witness value for {role} differs from its exact conditioned value"
                ),
            });
        }
    }

    let needs_opaque_membership = constraints.values().any(|constraint| {
        matches!(
            &constraint.kind,
            JointConstraintKind::ExternalManifold { .. }
                | JointConstraintKind::StochasticCoupling { .. }
        )
    }) || parameters
        .values()
        .any(|parameter| matches!(&parameter.treatment, ParameterTreatment::Derived { .. }));
    if witness.opaque_membership_claim.is_some() != needs_opaque_membership {
        return Err(IdentifiabilityError::InvalidText {
            field: "admissible-domain opaque membership claim coverage",
            detail: if needs_opaque_membership {
                "derived or opaque constraint semantics require one global witness-membership claim"
                    .to_string()
            } else {
                "an opaque membership claim is forbidden when every witness predicate is locally evaluable"
                    .to_string()
            },
        });
    }
    let expected_binding =
        admissible_domain_witness_binding(witness, parameters, constraints, sources)?;
    if let Some(claim) = &witness.opaque_membership_claim {
        let source = require_source(
            sources,
            &claim.source,
            "admissible-domain opaque membership claim",
        )?;
        if source.kind != SourceKind::AdmissibleDomainCertificate
            || source.content_hash_domain != ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN
            || source.contract_version != ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_VERSION
        {
            return Err(IdentifiabilityError::InvalidText {
                field: "admissible-domain opaque membership source contract",
                detail: format!(
                    "source {} must use the exact membership-certificate kind/domain/version",
                    source.key
                ),
            });
        }
        let expected_certificate_hash = hash_domain(
            ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN,
            expected_binding.as_bytes(),
        );
        if source.expected_hash != expected_certificate_hash {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "admissible-domain membership certificate content",
            });
        }
        if claim.witness_binding != Some(expected_binding) {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "admissible-domain witness/conjunction binding",
            });
        }
    }

    for (id, constraint) in constraints {
        let feasible = match &constraint.kind {
            JointConstraintKind::Affine {
                terms,
                relation,
                rhs_si,
                ..
            } => {
                let lhs = terms.iter().try_fold(0.0, |sum, term| {
                    let next = sum + term.coefficient * witness.values[&term.parameter];
                    next.is_finite().then_some(next)
                });
                lhs.is_some_and(|lhs| match relation {
                    ConstraintRelation::Equal => same_f64(lhs, *rhs_si),
                    ConstraintRelation::LessOrEqual => lhs <= *rhs_si,
                    ConstraintRelation::GreaterOrEqual => lhs >= *rhs_si,
                })
            }
            JointConstraintKind::Simplex {
                members, total_si, ..
            } => {
                let sum = members.iter().try_fold(0.0, |sum, role| {
                    let value = witness.values[role];
                    let domain = parameters[role].domain;
                    let next = sum + value;
                    // `Simplex` is the regular codimension-one composition
                    // chart, not its lower-dimensional boundary strata.  A
                    // boundary solution remains expressible as an affine
                    // equality, but must not inherit the simplex's intrinsic
                    // codimension claim.
                    (value > 0.0 && value > domain.lo && value < domain.hi && next.is_finite())
                        .then_some(next)
                });
                sum.is_some_and(|sum| same_f64(sum, *total_si))
            }
            JointConstraintKind::Ordered { members, strict } => members.windows(2).all(|pair| {
                let left = witness.values[&pair[0]];
                let right = witness.values[&pair[1]];
                if *strict { left < right } else { left <= right }
            }),
            JointConstraintKind::ExternalManifold { .. }
            | JointConstraintKind::StochasticCoupling { .. } => true,
        };
        if !feasible {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "admissible-domain witness feasibility",
                detail: format!("constructive witness violates joint constraint {id}"),
            });
        }
    }
    Ok(())
}

fn declared_parameter_cases(scope: &ParameterScopeBinding) -> Option<&BTreeSet<CaseId>> {
    match scope {
        ParameterScopeBinding::Global | ParameterScopeBinding::Specimen { .. } => None,
        ParameterScopeBinding::Cases(cases)
        | ParameterScopeBinding::MaterialLot { cases, .. }
        | ParameterScopeBinding::Field { cases, .. }
        | ParameterScopeBinding::Hierarchical { cases, .. } => Some(cases),
    }
}

fn validate_declared_parameter_cases(
    parameter: &StudyParameter,
    cases: &BTreeMap<CaseId, StudyCaseDocument>,
) -> Result<(), IdentifiabilityError> {
    let Some(scoped) = declared_parameter_cases(&parameter.scope) else {
        return Ok(());
    };
    if scoped.is_empty() {
        return Err(IdentifiabilityError::Cardinality {
            field: "parameter case scope",
            detail: format!(
                "parameter {} has an empty declared case scope",
                parameter.role
            ),
        });
    }
    for case in scoped {
        if !cases.contains_key(case) {
            return Err(IdentifiabilityError::UnknownReference {
                field: "parameter case scope",
                id: case.to_string(),
            });
        }
    }
    Ok(())
}

fn parameter_applicable_cases(
    parameter: &StudyParameter,
    cases: &BTreeMap<CaseId, StudyCaseDocument>,
) -> BTreeSet<CaseId> {
    match &parameter.scope {
        ParameterScopeBinding::Global => cases.keys().cloned().collect(),
        ParameterScopeBinding::Cases(scoped)
        | ParameterScopeBinding::MaterialLot { cases: scoped, .. }
        | ParameterScopeBinding::Field { cases: scoped, .. }
        | ParameterScopeBinding::Hierarchical { cases: scoped, .. } => scoped.clone(),
        ParameterScopeBinding::Specimen { case, .. } => BTreeSet::from([case.clone()]),
    }
}

impl IdentifiabilityProblemDocument {
    /// Validate and canonicalize a multi-case physical question.  This is
    /// structural admission only; [`Self::from_canonical_bytes`] returns the
    /// same unresolved type and cannot mint [`ProblemId`].
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        context_source: SourceKey,
        material_source: SourceKey,
        model_source: SourceKey,
        graph_source: SourceKey,
        joint_prior: Option<SourceKey>,
        sources: Vec<SourceRef>,
        parameters: Vec<StudyParameter>,
        mut constraints: Vec<JointConstraint>,
        mut admissible_domain: AdmissibleDomainWitness,
        cases: Vec<StudyCaseDocument>,
        mut influences: Vec<InfluenceDeclaration>,
        gauges: Vec<GaugeDeclaration>,
        gauge_compositions: Vec<GaugeCompositionDeclaration>,
        joint_noise: JointNoiseModel,
        mut data_reuse: DataReusePolicy,
    ) -> Result<Self, IdentifiabilityError> {
        let sources = insert_unique(sources, "source registry", |source| &source.key)?;
        require_source_kind(
            &sources,
            &context_source,
            SourceKind::ContextOfUse,
            "context source",
        )?;
        require_source_kind(
            &sources,
            &material_source,
            SourceKind::MaterialCard,
            "material source",
        )?;
        require_source_kind(
            &sources,
            &model_source,
            SourceKind::ConstitutiveModelCard,
            "model source",
        )?;
        require_source_kind(
            &sources,
            &graph_source,
            SourceKind::ConstitutiveGraph,
            "graph source",
        )?;
        if let Some(joint_prior) = &joint_prior {
            require_source_kind(
                &sources,
                joint_prior,
                SourceKind::ProbabilityMeasure,
                "problem joint prior",
            )?;
        }
        let parameters =
            insert_unique(parameters, "study parameters", |parameter| &parameter.role)?;
        if !parameters.values().any(|parameter| {
            matches!(
                &parameter.treatment,
                ParameterTreatment::Estimated
                    | ParameterTreatment::Profiled
                    | ParameterTreatment::Marginalized
            )
        }) {
            return Err(IdentifiabilityError::Cardinality {
                field: "inferential parameter targets",
                detail: "an identifiability problem needs at least one free inferential target"
                    .to_string(),
            });
        }
        validate_derived_parameter_dag(&parameters)?;
        for constraint in &mut constraints {
            if let JointConstraintKind::Affine { terms, .. } = &mut constraint.kind {
                terms.sort_by(|left, right| left.parameter.cmp(&right.parameter));
            }
        }
        let constraints = if constraints.is_empty() {
            BTreeMap::new()
        } else {
            insert_unique(constraints, "joint constraints", |constraint| {
                &constraint.id
            })?
        };
        for constraint in constraints.values() {
            validate_joint_constraint(constraint, &parameters, &sources)?;
        }
        if admissible_domain.opaque_membership_claim.is_some() {
            let binding = admissible_domain_witness_binding(
                &admissible_domain,
                &parameters,
                &constraints,
                &sources,
            )?;
            admissible_domain.bind_opaque_membership(binding)?;
        }
        validate_admissible_domain_witness(
            &admissible_domain,
            &parameters,
            &constraints,
            &sources,
        )?;
        let cases = insert_unique(cases, "study cases", |case| &case.id)?;
        for influence in &mut influences {
            if let DistributionFunctional::Correlation { left, right } = &mut influence.functional
                && right < left
            {
                core::mem::swap(left, right);
            }
        }
        let influences = if influences.is_empty() {
            BTreeMap::new()
        } else {
            insert_unique(influences, "influence declarations", |influence| {
                &influence.id
            })?
        };
        let gauges = if gauges.is_empty() {
            BTreeMap::new()
        } else {
            insert_unique(gauges, "gauge declarations", |gauge| &gauge.id)?
        };
        let gauge_compositions = if gauge_compositions.is_empty() {
            BTreeMap::new()
        } else {
            insert_unique(
                gauge_compositions,
                "gauge composition declarations",
                |composition| &composition.id,
            )?
        };
        let joint_noise = normalize_joint_noise(joint_noise)?;
        if let DataReusePolicy::Shared { groups } = &mut data_reuse {
            groups.sort_by(|left, right| {
                (&left.cases, &left.joint_likelihood, &left.justification).cmp(&(
                    &right.cases,
                    &right.joint_likelihood,
                    &right.justification,
                ))
            });
        }

        for parameter in parameters.values() {
            validate_declared_parameter_cases(parameter, &cases)?;
            match &parameter.owner {
                ParameterOwnerBinding::ConstitutiveModel => {}
                ParameterOwnerBinding::InitialState { state_path } => require_source_kind(
                    &sources,
                    state_path,
                    SourceKind::Assumption,
                    "initial-state owner",
                )?,
                ParameterOwnerBinding::Instrument { metrology, .. } => require_source_kind(
                    &sources,
                    metrology,
                    SourceKind::Metrology,
                    "instrument owner",
                )?,
                ParameterOwnerBinding::Discrepancy { family } => require_source_kind(
                    &sources,
                    family,
                    SourceKind::Discrepancy,
                    "discrepancy owner",
                )?,
                ParameterOwnerBinding::ControlledInput { protocol } => require_source_kind(
                    &sources,
                    protocol,
                    SourceKind::Protocol,
                    "controlled-input owner",
                )?,
                ParameterOwnerBinding::Population { hierarchy } => {
                    require_source_kind(&sources, hierarchy, SourceKind::Prior, "population owner")?
                }
            }
            match &parameter.scope {
                ParameterScopeBinding::Global
                | ParameterScopeBinding::Cases(_)
                | ParameterScopeBinding::MaterialLot { .. } => {}
                ParameterScopeBinding::Specimen { case, specimen } => {
                    let case_doc =
                        cases
                            .get(case)
                            .ok_or_else(|| IdentifiabilityError::UnknownReference {
                                field: "parameter specimen case",
                                id: case.to_string(),
                            })?;
                    if case_doc.specimen.id() != specimen {
                        return Err(IdentifiabilityError::UnknownReference {
                            field: "parameter specimen",
                            id: specimen.as_str().to_string(),
                        });
                    }
                }
                ParameterScopeBinding::Field { support, .. } => {
                    require_source_kind_in(
                        &sources,
                        support,
                        &[SourceKind::Geometry, SourceKind::ExternalManifold],
                        "field support",
                    )?;
                    return Err(IdentifiabilityError::InvalidText {
                        field: "field parameter carrier",
                        detail: format!(
                            "parameter {} declares Field scope, but v3 has no typed function-space carrier/discretization/reconstruction/measure contract; scalar conditioned, derived, and free actions are all rejected rather than silently treating a field as a constant scalar",
                            parameter.role
                        ),
                    });
                }
                ParameterScopeBinding::Hierarchical {
                    level, hierarchy, ..
                } => {
                    if *level == 0 {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "hierarchical level",
                            detail: "level zero is reserved for the global population".to_string(),
                        });
                    }
                    require_source_kind(
                        &sources,
                        hierarchy,
                        SourceKind::Prior,
                        "hierarchy source",
                    )?;
                }
            }
            match &parameter.treatment {
                ParameterTreatment::Conditioned(value) => require_source_kind_in(
                    &sources,
                    &value.source,
                    &[SourceKind::EvidenceReceipt, SourceKind::Metrology],
                    "conditioned value source",
                )?,
                ParameterTreatment::Derived { definition, .. } => require_source_kind(
                    &sources,
                    definition,
                    SourceKind::Constraint,
                    "derived parameter definition",
                )?,
                _ => {}
            }
        }

        let mut all_observations = BTreeMap::new();
        let mut modeled_discrepancy_parameters = BTreeSet::new();
        for (case_id, case) in &cases {
            let mut row_consumers = BTreeMap::<ObservationId, usize>::new();
            for observation in case.observations.values() {
                if let ObservationRows::Retrospective(rows) = &observation.rows {
                    for row in rows {
                        *row_consumers.entry(row.clone()).or_default() += 1;
                    }
                }
            }
            let repeated_rows = row_consumers
                .into_iter()
                .filter_map(|(row, consumers)| (consumers > 1).then_some(row))
                .collect::<BTreeSet<_>>();
            let declared_shared_rows = case
                .observation_sharing
                .iter()
                .flat_map(|group| group.rows.iter().cloned())
                .collect::<BTreeSet<_>>();
            if repeated_rows != declared_shared_rows {
                return Err(IdentifiabilityError::InvalidText {
                    field: "observation-sharing coverage",
                    detail: format!(
                        "case {case_id} must declare every and only multiply consumed raw row"
                    ),
                });
            }
            if case.protocol.state_schema_version
                != initial_state_schema_version(case.initial_state)
            {
                return Err(IdentifiabilityError::VersionMismatch {
                    field: "case initial-state/protocol schema",
                    expected: case.protocol.state_schema_version,
                    actual: initial_state_schema_version(case.initial_state),
                });
            }
            require_case_physics_source(
                &sources,
                &case.physics_sources.frame_transform,
                SourceKind::Geometry,
                case.specimen.frame().transform(),
                FRAME_TRANSFORM_SOURCE_DOMAIN,
                "case frame-transform source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.specimen_geometry,
                SourceKind::Geometry,
                case.specimen.geometry(),
                SPECIMEN_GEOMETRY_SOURCE_DOMAIN,
                "case specimen-geometry source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.specimen_process,
                SourceKind::Process,
                case.specimen.process(),
                SPECIMEN_PROCESS_SOURCE_DOMAIN,
                "case specimen-process source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.specimen_preparation,
                SourceKind::Process,
                case.specimen.preparation(),
                SPECIMEN_PREPARATION_SOURCE_DOMAIN,
                "case specimen-preparation source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.load_path,
                SourceKind::Protocol,
                case.protocol.load_path(),
                LOAD_PATH_SOURCE_DOMAIN,
                "case load-path source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.environment_path,
                SourceKind::Protocol,
                case.protocol.environment_path(),
                ENVIRONMENT_PATH_SOURCE_DOMAIN,
                "case environment-path source",
            )?;
            require_case_physics_source(
                &sources,
                &case.physics_sources.time_grid,
                SourceKind::Protocol,
                case.protocol.time_grid(),
                TIME_GRID_SOURCE_DOMAIN,
                "case time-grid source",
            )?;
            match (case.initial_state, &case.physics_sources.initial_state) {
                (InitialStateBinding::Zero { .. }, None) => {}
                (InitialStateBinding::Explicit { artifact, .. }, Some(source)) => {
                    require_case_physics_source(
                        &sources,
                        source,
                        SourceKind::Assumption,
                        artifact,
                        INITIAL_STATE_SOURCE_DOMAIN,
                        "case initial-state source",
                    )?;
                }
                (InitialStateBinding::Zero { .. }, Some(_)) => {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "zero initial-state source",
                    });
                }
                (InitialStateBinding::Explicit { .. }, None) => {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "case initial-state source",
                        id: case_id.to_string(),
                    });
                }
            }
            require_source_kind(
                &sources,
                &case.forward_model,
                SourceKind::ForwardModel,
                "case forward model",
            )?;
            for group in &case.observation_sharing {
                require_source_kind_in(
                    &sources,
                    &group.joint_likelihood,
                    &[SourceKind::Likelihood, SourceKind::ParameterizedLikelihood],
                    "observation-sharing joint likelihood",
                )?;
            }
            match &case.data {
                CaseDataDeclaration::Prospective => {
                    if matches!(&case.purpose, CasePurpose::BlindFalsification) {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "blind-falsification case data",
                            detail: format!(
                                "case {case_id} must bind retrospective blind rows and a release"
                            ),
                        });
                    }
                    if case.observations.values().any(|observation| {
                        !matches!(&observation.rows, ObservationRows::Prospective)
                    }) {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "prospective observation rows",
                            detail: format!("case {case_id} contains retrospective row IDs"),
                        });
                    }
                }
                CaseDataDeclaration::Retrospective {
                    experiment,
                    split,
                    parser,
                    preprocessing,
                    parser_version,
                    ..
                } => {
                    if matches!(&case.purpose, CasePurpose::ProspectiveDesign) {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "prospective-design case data",
                            detail: format!(
                                "case {case_id} is ProspectiveDesign but binds retrospective data"
                            ),
                        });
                    }
                    require_source_kind(
                        &sources,
                        experiment,
                        SourceKind::ExperimentArtifact,
                        "case experiment",
                    )?;
                    require_source_kind(
                        &sources,
                        split,
                        SourceKind::CalibrationSplit,
                        "case split",
                    )?;
                    require_source_kind(&sources, parser, SourceKind::Parser, "case parser")?;
                    let parser_source = require_source(&sources, parser, "case parser")?;
                    require_source_kind(
                        &sources,
                        preprocessing,
                        SourceKind::Preprocessing,
                        "case preprocessing",
                    )?;
                    if *parser_version == 0
                        || case.observations.values().any(|observation| {
                            !matches!(&observation.rows, ObservationRows::Retrospective(_))
                        })
                    {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "retrospective case",
                            detail: format!(
                                "case {case_id} needs a positive parser version and raw rows"
                            ),
                        });
                    }
                    if parser_source.contract_version != *parser_version {
                        return Err(IdentifiabilityError::VersionMismatch {
                            field: "case parser source contract",
                            expected: *parser_version,
                            actual: parser_source.contract_version,
                        });
                    }
                }
            }
            for (channel, observation) in &case.observations {
                if &observation.frame != case.specimen.frame() {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "observation/specimen frame",
                    });
                }
                if observation.protocol_version != case.protocol.version {
                    return Err(IdentifiabilityError::VersionMismatch {
                        field: "case observation protocol version",
                        expected: case.protocol.version,
                        actual: observation.protocol_version,
                    });
                }
                if observation.refinement_version != case.protocol.refinement_version {
                    return Err(IdentifiabilityError::VersionMismatch {
                        field: "case observation refinement version",
                        expected: case.protocol.refinement_version,
                        actual: observation.refinement_version,
                    });
                }
                if observation.clock != case.protocol.clock {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "case observation protocol clock",
                        detail: format!(
                            "observation {channel} names clock {}, but protocol {} names {}",
                            observation.clock, case.protocol.id, case.protocol.clock
                        ),
                    });
                }
                require_source_kind(
                    &sources,
                    &observation.unit_definition,
                    SourceKind::UnitDefinition,
                    "observation unit definition",
                )?;
                require_source_kind(
                    &sources,
                    &observation.operator,
                    SourceKind::ObservationOperator,
                    "observation operator",
                )?;
                let operator_source =
                    require_source(&sources, &observation.operator, "observation operator")?;
                require_source_kind(
                    &sources,
                    &observation.aggregation,
                    SourceKind::ObservationOperator,
                    "observation aggregation",
                )?;
                let aggregation_source = require_source(
                    &sources,
                    &observation.aggregation,
                    "observation aggregation",
                )?;
                for (field, source) in [
                    ("observation operator source contract", operator_source),
                    (
                        "observation aggregation source contract",
                        aggregation_source,
                    ),
                ] {
                    if source.contract_version != observation.operator_version {
                        return Err(IdentifiabilityError::VersionMismatch {
                            field,
                            expected: observation.operator_version,
                            actual: source.contract_version,
                        });
                    }
                }
                require_source_kind(
                    &sources,
                    &observation.sensor,
                    SourceKind::Metrology,
                    "observation sensor",
                )?;
                if let MarginalNoiseSpec::Empirical {
                    distribution,
                    finite_variance_model,
                    ..
                } = &observation.noise
                {
                    require_source_kind_in(
                        &sources,
                        distribution,
                        &[
                            SourceKind::Likelihood,
                            SourceKind::ParameterizedLikelihood,
                            SourceKind::EvidenceReceipt,
                        ],
                        "empirical noise",
                    )?;
                    require_source_kind(
                        &sources,
                        finite_variance_model,
                        SourceKind::EvidenceReceipt,
                        "empirical finite-variance model",
                    )?;
                }
                match &observation.missingness {
                    MissingnessAssumption::Complete { assumption } => require_source_kind(
                        &sources,
                        assumption,
                        SourceKind::Assumption,
                        "completeness assumption",
                    )?,
                    MissingnessAssumption::Modeled { mechanism } => require_source_kind_in(
                        &sources,
                        mechanism,
                        &[SourceKind::Likelihood, SourceKind::ParameterizedLikelihood],
                        "missingness mechanism",
                    )?,
                    MissingnessAssumption::Unknown { .. } => {}
                }
                let key = ObservationKey::new(case_id.clone(), channel.clone());
                all_observations.insert(key, observation);
            }
            for (channel, discrepancy) in &case.discrepancies {
                match discrepancy {
                    StudyDiscrepancy::Uncharacterized { reason } => {
                        validate_reason(reason, "discrepancy reason")?
                    }
                    StudyDiscrepancy::NotApplicable { basis } => match basis {
                        DiscrepancyInapplicability::PhysicalApplicability { assumption } => {
                            if !matches!(&case.data, CaseDataDeclaration::Retrospective { .. }) {
                                return Err(IdentifiabilityError::InvalidText {
                                    field: "physical discrepancy applicability",
                                    detail: format!(
                                        "case {case_id} channel {channel} is not a retrospective physical observation"
                                    ),
                                });
                            }
                            require_source_kind(
                                &sources,
                                assumption,
                                SourceKind::Assumption,
                                "physical discrepancy applicability assumption",
                            )?;
                        }
                        DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                            generator,
                            producer,
                            production_binding,
                            assumption,
                        } => {
                            if !matches!(&case.data, CaseDataDeclaration::Retrospective { .. })
                                || generator != &case.forward_model
                            {
                                return Err(IdentifiabilityError::SourceMismatch {
                                    field: "declared-synthetic self-model discrepancy",
                                });
                            }
                            require_source_kind(
                                &sources,
                                generator,
                                SourceKind::ForwardModel,
                                "declared-synthetic generator",
                            )?;
                            require_source_kind(
                                &sources,
                                production_binding,
                                SourceKind::ForwardModelProductionBinding,
                                "declared-synthetic producer/forward-model binding",
                            )?;
                            let binding_source = &sources[production_binding];
                            let binding_preimage = forward_model_production_binding_preimage(
                                producer,
                                &sources[generator],
                            )?;
                            if binding_source.content_hash_domain
                                != FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN
                                || binding_source.contract_version
                                    != FORWARD_MODEL_PRODUCTION_BINDING_VERSION
                                || binding_source.expected_hash
                                    != hash_domain(
                                        FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN,
                                        &binding_preimage,
                                    )
                            {
                                return Err(IdentifiabilityError::SourceMismatch {
                                    field: "declared-synthetic producer/forward-model production binding",
                                });
                            }
                            require_source_kind(
                                &sources,
                                assumption,
                                SourceKind::Assumption,
                                "declared-synthetic no-discrepancy assumption",
                            )?;
                        }
                        DiscrepancyInapplicability::ProspectiveDesign { assumption } => {
                            if !matches!(&case.data, CaseDataDeclaration::Prospective)
                                || !matches!(&case.purpose, CasePurpose::ProspectiveDesign)
                            {
                                return Err(IdentifiabilityError::InvalidText {
                                    field: "prospective discrepancy applicability",
                                    detail: format!(
                                        "case {case_id} channel {channel} must be an exact prospective-design case"
                                    ),
                                });
                            }
                            require_source_kind(
                                &sources,
                                assumption,
                                SourceKind::Assumption,
                                "prospective discrepancy applicability assumption",
                            )?;
                        }
                    },
                    StudyDiscrepancy::AssumedZero { assumption } => require_source_kind(
                        &sources,
                        assumption,
                        SourceKind::Assumption,
                        "zero-discrepancy assumption",
                    )?,
                    StudyDiscrepancy::Modeled {
                        family,
                        parameters: discrepancy_parameters,
                        support,
                        confounding_guard,
                    } => {
                        require_source_kind(
                            &sources,
                            family,
                            SourceKind::Discrepancy,
                            "modeled discrepancy family",
                        )?;
                        require_source_kind_in(
                            &sources,
                            support,
                            &[SourceKind::Geometry, SourceKind::ExternalManifold],
                            "modeled discrepancy support",
                        )?;
                        require_source_kind(
                            &sources,
                            confounding_guard,
                            SourceKind::Constraint,
                            "modeled discrepancy confounding guard",
                        )?;
                        if discrepancy_parameters.is_empty() {
                            return Err(IdentifiabilityError::Cardinality {
                                field: "discrepancy parameters",
                                detail: "modeled discrepancy needs explicit parameter roles"
                                    .to_string(),
                            });
                        }
                        for role in discrepancy_parameters {
                            let parameter = parameters.get(role).ok_or_else(|| {
                                IdentifiabilityError::UnknownReference {
                                    field: "discrepancy parameter",
                                    id: role.to_string(),
                                }
                            })?;
                            if !matches!(
                                &parameter.owner,
                                ParameterOwnerBinding::Discrepancy {
                                    family: owner_family
                                } if owner_family == family
                            ) {
                                return Err(IdentifiabilityError::InvalidText {
                                    field: "discrepancy parameter owner",
                                    detail: format!(
                                        "parameter {role} is not owned by modeled family {family}"
                                    ),
                                });
                            }
                            if !parameter_active_in_cases(
                                parameter,
                                &BTreeSet::from([case_id.clone()]),
                            ) {
                                return Err(IdentifiabilityError::InvalidText {
                                    field: "discrepancy parameter/case scope",
                                    detail: format!(
                                        "modeled discrepancy in case {case_id} uses parameter {role} outside its exact applicability"
                                    ),
                                });
                            }
                            modeled_discrepancy_parameters.insert(role.clone());
                        }
                    }
                }
            }
        }
        for parameter in parameters.values() {
            if matches!(&parameter.owner, ParameterOwnerBinding::Discrepancy { .. })
                && !modeled_discrepancy_parameters.contains(&parameter.role)
            {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "modeled discrepancy parameter",
                    id: parameter.role.to_string(),
                });
            }
            let applicable = parameter_applicable_cases(parameter, &cases);
            let owner_is_bound_in_case = |case: &StudyCaseDocument| match &parameter.owner {
                ParameterOwnerBinding::ConstitutiveModel
                | ParameterOwnerBinding::Population { .. } => true,
                ParameterOwnerBinding::InitialState { state_path } => {
                    case.physics_sources.initial_state.as_ref() == Some(state_path)
                }
                ParameterOwnerBinding::Instrument {
                    instrument,
                    acquisition_channel,
                    metrology,
                } => case.observations.values().any(|observation| {
                    &observation.instrument == instrument
                        && &observation.acquisition_channel == acquisition_channel
                        && &observation.sensor == metrology
                }),
                ParameterOwnerBinding::Discrepancy { family } => {
                    case.discrepancies.values().any(|discrepancy| {
                        matches!(
                            discrepancy,
                            StudyDiscrepancy::Modeled {
                                family: candidate,
                                parameters,
                                ..
                            } if candidate == family && parameters.contains(&parameter.role)
                        )
                    })
                }
                ParameterOwnerBinding::ControlledInput { protocol } => [
                    &case.physics_sources.load_path,
                    &case.physics_sources.environment_path,
                    &case.physics_sources.time_grid,
                ]
                .contains(&protocol),
            };
            if let Some(case_id) = applicable
                .iter()
                .find(|case_id| !owner_is_bound_in_case(&cases[*case_id]))
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "parameter owner/case binding",
                    detail: format!(
                        "parameter {} owner is not physically bound in applicable case {case_id}",
                        parameter.role
                    ),
                });
            }
            if let ParameterOwnerBinding::Population { hierarchy } = &parameter.owner
                && let ParameterScopeBinding::Hierarchical {
                    hierarchy: scoped_hierarchy,
                    ..
                } = &parameter.scope
                && hierarchy != scoped_hierarchy
            {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "population owner/hierarchical scope",
                });
            }
        }

        let mut influenced_parameters = BTreeSet::new();
        for influence in influences.values() {
            let parameter = parameters.get(&influence.parameter).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "influence parameter",
                    id: influence.parameter.to_string(),
                }
            })?;
            if let DistributionFunctional::Correlation { left, right } = &influence.functional {
                let left_observation = observation_for(&cases, left)?;
                let right_observation = observation_for(&cases, right)?;
                if !left_observation.noise.finite_standard_deviation()
                    || !right_observation.noise.finite_standard_deviation()
                {
                    return Err(IdentifiabilityError::Covariance {
                        detail: "Pearson-correlation influence requires two finite-second-moment marginals"
                            .to_string(),
                    });
                }
            }
            for key in functional_observations(&influence.functional) {
                let observation = observation_for(&cases, key)?;
                if !parameter_active_in_cases(parameter, &BTreeSet::from([key.case.clone()])) {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "influence parameter/case scope",
                        detail: format!(
                            "influence {} targets case {} outside parameter {} applicability",
                            influence.id, key.case, influence.parameter
                        ),
                    });
                }
                match &influence.functional {
                    DistributionFunctional::Correlation { left, right } if left == right => {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "correlation functional",
                            detail: "self-correlation is a constant, not an identifiability route"
                                .to_string(),
                        });
                    }
                    DistributionFunctional::MissingnessLogit { .. }
                        if !matches!(
                            &observation.missingness,
                            MissingnessAssumption::Modeled { mechanism }
                                if sources[mechanism].kind == SourceKind::ParameterizedLikelihood
                        ) =>
                    {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "missingness functional",
                            detail: "MissingnessLogit influence needs a modeled mechanism with exact ParameterizedLikelihood semantics"
                                .to_string(),
                        });
                    }
                    DistributionFunctional::LogScale { .. }
                        if !matches!(
                            &observation.noise,
                            MarginalNoiseSpec::Empirical { distribution, .. }
                                if sources[distribution].kind == SourceKind::ParameterizedLikelihood
                        ) =>
                    {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "log-scale functional",
                            detail: "LogScale influence needs an exact parameterized marginal likelihood; fixed Gaussian/Student-t/bounded parameters are not influence routes"
                                .to_string(),
                        });
                    }
                    DistributionFunctional::CensoringLogit { .. } => {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "censoring functional",
                            detail: if observation.saturation.is_none() {
                                "CensoringLogit influence requires saturation and an exact parameterized censoring law; saturation is absent"
                            } else {
                                "CensoringLogit influence requires an exact parameterized censoring-law source; a fixed saturation interval alone is insufficient"
                            }
                            .to_string(),
                        });
                    }
                    _ => {}
                }
            }
            match &influence.representation {
                InfluenceRepresentation::Direct => {}
                InfluenceRepresentation::StateMediated { state_path } => require_source_kind(
                    &sources,
                    state_path,
                    SourceKind::Assumption,
                    "state-mediated influence",
                )?,
                InfluenceRepresentation::Composite { operator, inputs } => {
                    require_source_kind(
                        &sources,
                        operator,
                        SourceKind::InfluenceComposition,
                        "composite influence operator",
                    )?;
                    if inputs.is_empty() || inputs.contains(&influence.id) {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "composite influence inputs",
                            detail: format!("influence {} has empty or self input", influence.id),
                        });
                    }
                    for input in inputs {
                        let input_declaration = influences.get(input).ok_or_else(|| {
                            IdentifiabilityError::UnknownReference {
                                field: "composite influence input",
                                id: input.to_string(),
                            }
                        })?;
                        if input_declaration.parameter != influence.parameter {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "composite influence parameter",
                                detail: format!(
                                    "composite influence {} targets parameter {} but input {} targets {}; cross-parameter chain rules require a dedicated typed composition schema",
                                    influence.id,
                                    influence.parameter,
                                    input_declaration.id,
                                    input_declaration.parameter
                                ),
                            });
                        }
                    }
                }
                InfluenceRepresentation::ExternalDefinition { definition } => {
                    require_source_kind(
                        &sources,
                        definition,
                        SourceKind::Constraint,
                        "external influence definition",
                    )?;
                }
            }
            influenced_parameters.insert(influence.parameter.clone());
        }
        for parameter in parameters.values() {
            let free = matches!(
                &parameter.treatment,
                ParameterTreatment::Estimated
                    | ParameterTreatment::Profiled
                    | ParameterTreatment::Marginalized
            );
            match (&parameter.influence_coverage, free) {
                (InfluenceCoverage::Declared, _)
                    if !influenced_parameters.contains(&parameter.role) =>
                {
                    return Err(IdentifiabilityError::DisconnectedEstimatedParameter {
                        parameter: parameter.role.clone(),
                    });
                }
                (InfluenceCoverage::IntentionallyAbsent { .. }, true)
                    if influenced_parameters.contains(&parameter.role) =>
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "parameter influence coverage",
                        detail: format!(
                            "parameter {} both declares and denies influence routes",
                            parameter.role
                        ),
                    });
                }
                (InfluenceCoverage::NotApplicable { .. }, false)
                    if influenced_parameters.contains(&parameter.role) =>
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "parameter influence coverage",
                        detail: format!(
                            "parameter {} marks influence not applicable but owns a route",
                            parameter.role
                        ),
                    });
                }
                _ => {}
            }
        }

        // Composite influence declarations form a DAG; use an iterative
        // topological pass so an adversarial but bounded chain cannot consume
        // the process stack.
        let mut incoming = influences
            .keys()
            .cloned()
            .map(|id| (id, 0_usize))
            .collect::<BTreeMap<_, _>>();
        let mut dependents = BTreeMap::<InfluenceId, BTreeSet<InfluenceId>>::new();
        for (id, influence) in &influences {
            if let InfluenceRepresentation::Composite { inputs, .. } = &influence.representation {
                incoming.insert(id.clone(), inputs.len());
                for input in inputs {
                    dependents
                        .entry(input.clone())
                        .or_default()
                        .insert(id.clone());
                }
            }
        }
        let mut ready = incoming
            .iter()
            .filter_map(|(id, degree)| (*degree == 0).then_some(id.clone()))
            .collect::<BTreeSet<_>>();
        let mut visited = 0_usize;
        while let Some(id) = ready.pop_first() {
            visited += 1;
            if let Some(children) = dependents.get(&id) {
                for child in children {
                    let degree = incoming.get_mut(child).expect("known composite child");
                    *degree -= 1;
                    if *degree == 0 {
                        ready.insert(child.clone());
                    }
                }
            }
        }
        if visited != influences.len() {
            let id = incoming
                .iter()
                .find_map(|(id, degree)| (*degree > 0).then_some(id))
                .expect("an unvisited influence has positive indegree");
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "composite influence graph",
                detail: format!("cycle reaches influence {id}"),
            });
        }

        for gauge in gauges.values() {
            for (axes, cell) in &gauge.validity.cells {
                for (case, extent_support) in &cell.case_obstruction_support {
                    if !cases.contains_key(case) {
                        return Err(IdentifiabilityError::UnknownReference {
                            field: "gauge invariant case",
                            id: case.to_string(),
                        });
                    }
                    if !extent_support
                        .global_obstruction_parameters
                        .is_subset(&gauge.members)
                        || !extent_support
                            .local_obstruction_parameters
                            .is_subset(&extent_support.global_obstruction_parameters)
                    {
                        return Err(IdentifiabilityError::InvalidGauge {
                            gauge: gauge.id.clone(),
                            detail: "local obstruction support must be a subset of global obstruction support, and both must stay inside gauge members"
                                .to_string(),
                        });
                    }
                    if !regular_orbit_support_compatible(&gauge.orbit_geometry, extent_support) {
                        return Err(IdentifiabilityError::InvalidGauge {
                            gauge: gauge.id.clone(),
                            detail: "regular orbit geometry and local/global obstruction support disagree; positive-dimensional regular orbits obstruct locally and globally, finite discrete proper/local-covering regular orbits obstruct only globally, every nontrivial regular orbit obstructs globally, and singular/nonproper exceptions require Stratified profile authority"
                                .to_string(),
                        });
                    }
                    if let Some(member) =
                        extent_support
                            .global_obstruction_parameters
                            .iter()
                            .find(|member| {
                                parameters.get(*member).is_none_or(|parameter| {
                                    !parameter_active_in_cases(
                                        parameter,
                                        &BTreeSet::from([case.clone()]),
                                    )
                                })
                            })
                    {
                        return Err(IdentifiabilityError::InvalidGauge {
                            gauge: gauge.id.clone(),
                            detail: format!(
                                "cell-local obstruction parameter {member} is inactive in case {case}"
                            ),
                        });
                    }
                }
                validate_gauge_applicability_sources(axes, &sources)?;
            }
            let gauge_cases = gauge
                .validity
                .cells
                .values()
                .flat_map(|cell| cell.case_obstruction_support.keys())
                .cloned()
                .collect::<BTreeSet<_>>();
            for member in &gauge.members {
                let parameter = parameters.get(member).ok_or_else(|| {
                    IdentifiabilityError::UnknownReference {
                        field: "gauge member",
                        id: member.to_string(),
                    }
                })?;
                if !matches!(
                    &parameter.treatment,
                    ParameterTreatment::Estimated
                        | ParameterTreatment::Profiled
                        | ParameterTreatment::Marginalized
                ) {
                    return Err(IdentifiabilityError::InvalidGauge {
                        gauge: gauge.id.clone(),
                        detail: format!(
                            "gauge member {member} must be a free inferential coordinate; conditioned and derived coordinates need an explicit induced-action schema"
                        ),
                    });
                }
                if !parameter_active_in_cases(parameter, &gauge_cases) {
                    return Err(IdentifiabilityError::InvalidGauge {
                        gauge: gauge.id.clone(),
                        detail: format!(
                            "gauge member {member} is inactive throughout the declared invariant cases"
                        ),
                    });
                }
            }
            require_source_kind(
                &sources,
                &gauge.action,
                SourceKind::GaugeAction,
                "gauge action",
            )?;
            match &gauge.status {
                GaugeStatus::Candidate { rationale } => require_source_kind(
                    &sources,
                    rationale,
                    SourceKind::GaugeHypothesis,
                    "gauge candidate rationale",
                )?,
                GaugeStatus::Assumed { assumption } => require_source_kind(
                    &sources,
                    assumption,
                    SourceKind::Assumption,
                    "assumed gauge authority",
                )?,
            }
            validate_gauge_algebra_orbit_sources(&gauge.algebra, &gauge.orbit_geometry, &sources)?;
        }

        for composition in gauge_compositions.values() {
            require_source_kind(
                &sources,
                &composition.law,
                SourceKind::GaugeComposition,
                "gauge composition law",
            )?;
            match &composition.status {
                GaugeStatus::Candidate { rationale } => require_source_kind(
                    &sources,
                    rationale,
                    SourceKind::GaugeHypothesis,
                    "gauge composition candidate rationale",
                )?,
                GaugeStatus::Assumed { assumption } => require_source_kind(
                    &sources,
                    assumption,
                    SourceKind::Assumption,
                    "assumed gauge composition authority",
                )?,
            }
            for (axes, cell) in &composition.validity.cells {
                for case in cell.case_obstruction_support.keys() {
                    if !cases.contains_key(case) {
                        return Err(IdentifiabilityError::UnknownReference {
                            field: "gauge composition invariant case",
                            id: case.to_string(),
                        });
                    }
                }
                validate_gauge_applicability_sources(axes, &sources)?;
            }
            let mut support = BTreeSet::new();
            let mut all_assumed = matches!(&composition.status, GaugeStatus::Assumed { .. });
            for member in &composition.members {
                let gauge =
                    gauges
                        .get(member)
                        .ok_or_else(|| IdentifiabilityError::UnknownReference {
                            field: "gauge composition member",
                            id: member.to_string(),
                        })?;
                if composition.validity.cells.iter().any(|(axes, cell)| {
                    gauge.validity.cells.get(axes).is_none_or(|member_cell| {
                        !cell
                            .case_obstruction_support
                            .keys()
                            .all(|case| member_cell.case_obstruction_support.contains_key(case))
                    })
                }) {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "gauge composition validity",
                        detail: format!(
                            "composition {} exceeds member {member} validity scope",
                            composition.law
                        ),
                    });
                }
                support.extend(gauge.members.iter().cloned());
                all_assumed &= matches!(&gauge.status, GaugeStatus::Assumed { .. });
            }
            if matches!(&composition.status, GaugeStatus::Assumed { .. }) && !all_assumed {
                return Err(IdentifiabilityError::InvalidText {
                    field: "assumed gauge composition status",
                    detail: format!(
                        "composition {} cannot be Assumed while any member action remains Candidate",
                        composition.id
                    ),
                });
            }
            for (axes, cell) in &composition.validity.cells {
                for (case, extent_support) in &cell.case_obstruction_support {
                    // A generated/noncommuting action may obstruct coordinates
                    // absent from every member action's point-local support.
                    // Its exact composition law is authoritative for that
                    // effective support; the only structural bound is the
                    // union of all coordinates the member actions may move.
                    let member_local = composition
                        .members
                        .iter()
                        .flat_map(|member| {
                            gauges[member].validity.cells[axes].case_obstruction_support[case]
                                .local_obstruction_parameters
                                .iter()
                                .cloned()
                        })
                        .collect::<BTreeSet<_>>();
                    let member_global = composition
                        .members
                        .iter()
                        .flat_map(|member| {
                            gauges[member].validity.cells[axes].case_obstruction_support[case]
                                .global_obstruction_parameters
                                .iter()
                                .cloned()
                        })
                        .collect::<BTreeSet<_>>();
                    if !extent_support
                        .global_obstruction_parameters
                        .is_subset(&support)
                        || !extent_support
                            .local_obstruction_parameters
                            .is_subset(&extent_support.global_obstruction_parameters)
                        || match &composition.kind {
                            GaugeCompositionKind::IndependentProduct => {
                                extent_support.local_obstruction_parameters != member_local
                                    || extent_support.global_obstruction_parameters != member_global
                            }
                            GaugeCompositionKind::Generated => {
                                !member_local
                                    .is_subset(&extent_support.local_obstruction_parameters)
                                    || !member_global
                                        .is_subset(&extent_support.global_obstruction_parameters)
                            }
                        }
                        || !regular_orbit_support_compatible(
                            &composition.effective_orbit_geometry,
                            extent_support,
                        )
                    {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "gauge composition moved support",
                            detail: format!(
                                "composition {} has invalid effective moved support in case {case}",
                                composition.law
                            ),
                        });
                    }
                    if let Some(member) =
                        extent_support
                            .global_obstruction_parameters
                            .iter()
                            .find(|member| {
                                parameters.get(*member).is_none_or(|parameter| {
                                    !parameter_active_in_cases(
                                        parameter,
                                        &BTreeSet::from([case.clone()]),
                                    )
                                })
                            })
                    {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "gauge composition moved support",
                            detail: format!(
                                "composition {} moves case-inactive parameter {member}",
                                composition.law
                            ),
                        });
                    }
                }
            }
            validate_gauge_algebra_orbit_sources(
                &composition.effective_algebra,
                &composition.effective_orbit_geometry,
                &sources,
            )?;
            let principal = principal_gauge_orbit(&composition.effective_orbit_geometry);
            if !gauge_algebra_orbit_compatible(
                &composition.effective_algebra,
                principal,
                all_assumed
                    && matches!(
                        &composition.effective_orbit_geometry,
                        GaugeOrbitGeometry::Regular { .. }
                    ),
            ) {
                return Err(IdentifiabilityError::InvalidText {
                    field: "gauge composition geometry",
                    detail: format!(
                        "composition {} has incompatible effective principal-orbit invariants",
                        composition.law
                    ),
                });
            }
            validate_independent_product_invariants(composition, &gauges)?;
        }

        for case in cases.keys() {
            let axes_set = gauges
                .values()
                .flat_map(|gauge| {
                    gauge
                        .validity
                        .cells
                        .iter()
                        .filter(|(_, cell)| cell.case_obstruction_support.contains_key(case))
                        .map(|(axes, _)| axes.clone())
                })
                .collect::<BTreeSet<_>>();
            for axes in axes_set {
                let active_assumed = gauges
                    .values()
                    .filter(|gauge| {
                        matches!(&gauge.status, GaugeStatus::Assumed { .. })
                            && gauge.validity.cells.get(&axes).is_some_and(|cell| {
                                cell.case_obstruction_support.contains_key(case)
                            })
                    })
                    .map(|gauge| gauge.id.clone())
                    .collect::<BTreeSet<_>>();
                if active_assumed.len() < 2 {
                    continue;
                }
                let matching = gauge_compositions
                    .values()
                    .filter(|composition| {
                        composition.members == active_assumed
                            && matches!(&composition.status, GaugeStatus::Assumed { .. })
                            && composition.validity.cells.get(&axes).is_some_and(|cell| {
                                cell.case_obstruction_support.contains_key(case)
                            })
                    })
                    .count();
                if matching != 1 {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "assumed gauge composition",
                        detail: format!(
                            "simultaneously active assumed gauge system {:?} needs exactly one exact assumed composition declaration, found {matching}",
                            active_assumed
                        ),
                    });
                }
            }
        }

        match &joint_noise {
            JointNoiseModel::Independent { assumption } => require_source_kind(
                &sources,
                assumption,
                SourceKind::Assumption,
                "independent-noise assumption",
            )?,
            JointNoiseModel::DenseCorrelation {
                order,
                correlation,
                model,
            } => {
                require_source_kind(
                    &sources,
                    model,
                    SourceKind::Likelihood,
                    "dense correlation model",
                )?;
                let unique = order.iter().cloned().collect::<BTreeSet<_>>();
                let all = all_observations.keys().cloned().collect::<BTreeSet<_>>();
                if order.len() != all.len()
                    || unique != all
                    || correlation.dimension() != order.len()
                    || order
                        .iter()
                        .any(|key| !all_observations[key].noise.finite_standard_deviation())
                {
                    return Err(IdentifiabilityError::Covariance {
                        detail: "dense correlation needs every composite channel exactly once and finite marginal standard deviations"
                            .to_string(),
                    });
                }
                for index in 0..order.len() {
                    if !same_f64(matrix_get(correlation, index, index), 1.0) {
                        return Err(IdentifiabilityError::Covariance {
                            detail: format!("correlation diagonal {index} is not exactly one"),
                        });
                    }
                }
            }
            JointNoiseModel::ExternalKernel { model } => require_source_kind_in(
                &sources,
                model,
                &[SourceKind::Likelihood, SourceKind::ParameterizedLikelihood],
                "external noise kernel",
            )?,
            JointNoiseModel::Unknown { reason } => {
                validate_reason(reason, "unknown joint noise reason")?
            }
        }

        if influences.values().any(|influence| {
            matches!(
                &influence.functional,
                DistributionFunctional::Correlation { .. }
            )
        }) {
            let parameterized = matches!(
                &joint_noise,
                JointNoiseModel::ExternalKernel { model }
                    if sources[model].kind == SourceKind::ParameterizedLikelihood
            );
            if !parameterized {
                return Err(IdentifiabilityError::InvalidText {
                    field: "correlation influence/joint likelihood",
                    detail: "correlation influence requires an ExternalKernel whose source has ParameterizedLikelihood semantics; Independent fixes correlation to zero, DenseCorrelation stores a fixed matrix, and a generic/unknown kernel does not declare parameter dependence"
                        .to_string(),
                });
            }
        }

        let declared_sharing_likelihoods = cases
            .values()
            .flat_map(|case| {
                case.observation_sharing
                    .iter()
                    .map(|group| group.joint_likelihood.clone())
            })
            .chain(match &data_reuse {
                DataReusePolicy::Disjoint => Vec::new(),
                DataReusePolicy::Shared { groups } => groups
                    .iter()
                    .map(|group| group.joint_likelihood.clone())
                    .collect(),
            })
            .collect::<BTreeSet<_>>();
        if !declared_sharing_likelihoods.is_empty() {
            let global_model = match &joint_noise {
                JointNoiseModel::DenseCorrelation { model, .. }
                | JointNoiseModel::ExternalKernel { model } => model,
                JointNoiseModel::Independent { .. } | JointNoiseModel::Unknown { .. } => {
                    return Err(IdentifiabilityError::Covariance {
                        detail:
                            "shared raw data requires one explicit global joint-likelihood model"
                                .to_string(),
                    });
                }
            };
            if declared_sharing_likelihoods != BTreeSet::from([global_model.clone()]) {
                return Err(IdentifiabilityError::Covariance {
                    detail: "every sharing declaration must name the exact global joint-likelihood model"
                        .to_string(),
                });
            }
        }

        match &data_reuse {
            DataReusePolicy::Disjoint => {
                let mut seen = BTreeMap::<ContentHash, CaseId>::new();
                for (case_id, case) in &cases {
                    if let Some(experiment) = retrospective_experiment(case) {
                        let hash = sources[experiment].expected_hash;
                        if let Some(other) = seen.insert(hash, case_id.clone()) {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "data reuse policy",
                                detail: format!(
                                    "cases {other} and {case_id} reuse one experiment under Disjoint"
                                ),
                            });
                        }
                    }
                }
            }
            DataReusePolicy::Shared { groups } => {
                if groups.is_empty() || groups.len() > MAX_IDENTIFIABILITY_ITEMS {
                    return Err(IdentifiabilityError::Cardinality {
                        field: "data sharing groups",
                        detail: "Shared policy needs bounded nonempty groups".to_string(),
                    });
                }
                let mut membership = BTreeMap::<CaseId, usize>::new();
                let mut shared_hash_owners = BTreeMap::<ContentHash, usize>::new();
                for (index, group) in groups.iter().enumerate() {
                    require_source_kind_in(
                        &sources,
                        &group.joint_likelihood,
                        &[SourceKind::Likelihood, SourceKind::ParameterizedLikelihood],
                        "sharing-group likelihood",
                    )?;
                    for case_id in &group.cases {
                        if membership.insert(case_id.clone(), index).is_some() {
                            return Err(IdentifiabilityError::Duplicate {
                                field: "data sharing group membership",
                                id: case_id.to_string(),
                            });
                        }
                        let case = cases.get(case_id).ok_or_else(|| {
                            IdentifiabilityError::UnknownReference {
                                field: "data sharing case",
                                id: case_id.to_string(),
                            }
                        })?;
                        let experiment = retrospective_experiment(case).ok_or_else(|| {
                            IdentifiabilityError::InvalidText {
                                field: "data sharing case",
                                detail: format!("prospective case {case_id} cannot share raw data"),
                            }
                        })?;
                        let hash = sources[experiment].expected_hash;
                        if let Some(other) = shared_hash_owners.insert(hash, index)
                            && other != index
                        {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "data reuse policy",
                                detail: format!(
                                    "sharing groups {other} and {index} reuse one experiment"
                                ),
                            });
                        }
                    }
                }
                let mut ungrouped = BTreeMap::<ContentHash, CaseId>::new();
                for (case_id, case) in &cases {
                    if membership.contains_key(case_id) {
                        continue;
                    }
                    if let Some(experiment) = retrospective_experiment(case) {
                        let hash = sources[experiment].expected_hash;
                        if shared_hash_owners.contains_key(&hash) {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "data reuse policy",
                                detail: format!(
                                    "ungrouped case {case_id} reuses an experiment owned by a sharing group"
                                ),
                            });
                        }
                        if let Some(other) = ungrouped.insert(hash, case_id.clone()) {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "data reuse policy",
                                detail: format!(
                                    "ungrouped cases {other} and {case_id} reuse one experiment"
                                ),
                            });
                        }
                    }
                }
            }
        }

        let reachable = problem_source_reachability(
            &context_source,
            &material_source,
            &model_source,
            &graph_source,
            joint_prior.as_ref(),
            &parameters,
            &constraints,
            &admissible_domain,
            &cases,
            &influences,
            &gauges,
            &gauge_compositions,
            &joint_noise,
            &data_reuse,
        );
        let registered = sources.keys().cloned().collect::<BTreeSet<_>>();
        if reachable != registered {
            let detail = registered.difference(&reachable).next().map_or_else(
                || "a referenced source is absent from the registry".to_string(),
                |unused| format!("source {unused} is registered but unreachable"),
            );
            return Err(IdentifiabilityError::InvalidText {
                field: "source registry closure",
                detail,
            });
        }

        let document = Self {
            schema_version: IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION,
            context_source,
            material_source,
            model_source,
            graph_source,
            joint_prior,
            sources,
            parameters,
            constraints,
            admissible_domain,
            cases,
            influences,
            gauges,
            gauge_compositions,
            joint_noise,
            data_reuse,
        };
        validate_problem_structural_budget(&document)?;
        // Construction itself enforces the canonical byte budget; callers do
        // not need to discover an oversized identity only when hashing later.
        let _ = encode_problem(&document)?;
        Ok(document)
    }

    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub const fn context_source(&self) -> &SourceKey {
        &self.context_source
    }

    #[must_use]
    pub const fn material_source(&self) -> &SourceKey {
        &self.material_source
    }

    #[must_use]
    pub const fn model_source(&self) -> &SourceKey {
        &self.model_source
    }

    #[must_use]
    pub const fn graph_source(&self) -> &SourceKey {
        &self.graph_source
    }

    #[must_use]
    pub const fn joint_prior(&self) -> Option<&SourceKey> {
        self.joint_prior.as_ref()
    }

    #[must_use]
    pub const fn sources(&self) -> &BTreeMap<SourceKey, SourceRef> {
        &self.sources
    }

    #[must_use]
    pub const fn parameters(&self) -> &BTreeMap<ParameterRoleId, StudyParameter> {
        &self.parameters
    }

    #[must_use]
    pub const fn cases(&self) -> &BTreeMap<CaseId, StudyCaseDocument> {
        &self.cases
    }

    #[must_use]
    pub const fn constraints(&self) -> &BTreeMap<ConstraintId, JointConstraint> {
        &self.constraints
    }

    #[must_use]
    pub const fn admissible_domain(&self) -> &AdmissibleDomainWitness {
        &self.admissible_domain
    }

    #[must_use]
    pub const fn influences(&self) -> &BTreeMap<InfluenceId, InfluenceDeclaration> {
        &self.influences
    }

    #[must_use]
    pub const fn gauges(&self) -> &BTreeMap<GaugeClassId, GaugeDeclaration> {
        &self.gauges
    }

    #[must_use]
    pub const fn gauge_compositions(
        &self,
    ) -> &BTreeMap<GaugeCompositionId, GaugeCompositionDeclaration> {
        &self.gauge_compositions
    }

    #[must_use]
    pub const fn joint_noise(&self) -> &JointNoiseModel {
        &self.joint_noise
    }

    #[must_use]
    pub const fn data_reuse(&self) -> &DataReusePolicy {
        &self.data_reuse
    }

    /// Derived derivative quantity for an influence functional with respect to
    /// its physical parameter.  Log-scale, correlation, missingness-logit, and
    /// censoring-logit functionals are dimensionless by definition.
    pub fn influence_derivative_quantity(
        &self,
        id: &InfluenceId,
    ) -> Result<QuantitySpec, IdentifiabilityError> {
        let influence =
            self.influences
                .get(id)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field: "influence derivative",
                    id: id.to_string(),
                })?;
        let output_dims = match &influence.functional {
            DistributionFunctional::Location { observation } => {
                observation_for(&self.cases, observation)?.quantity.dims()
            }
            DistributionFunctional::LogScale { .. }
            | DistributionFunctional::Correlation { .. }
            | DistributionFunctional::MissingnessLogit { .. }
            | DistributionFunctional::CensoringLogit { .. } => Dims([0; 6]),
        };
        let input_dims = self.parameters[&influence.parameter].quantity.dims();
        let dims = checked_derivative_dims(output_dims, input_dims).ok_or_else(|| {
            IdentifiabilityError::InvalidNumeric {
                field: "influence derivative quantity",
                detail: "dimension exponent overflow".to_string(),
            }
        })?;
        Ok(QuantitySpec::dimensional(dims))
    }

    /// Canonical unresolved bytes.  These bytes contain no source authority.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, IdentifiabilityError> {
        encode_problem(self)
    }

    /// Decode and fully revalidate an unresolved document.  This method cannot
    /// return an admitted problem or mint a [`ProblemId`].
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, IdentifiabilityError> {
        decode_problem(bytes)
    }
}

/// Concrete V&V sources for one retrospective case.
#[derive(Clone, Copy)]
pub struct CaseSourceBundle<'a> {
    experiment: &'a ExperimentArtifact,
    split: &'a CalibrationSplit,
    blind_release: Option<&'a BlindReleaseReceipt>,
}

impl fmt::Debug for CaseSourceBundle<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CaseSourceBundle")
            .field(
                "observation_count",
                &self.experiment.observation_ids().len(),
            )
            .field("calibration_count", &self.split.calibration_ids().len())
            .field("validation_count", &self.split.validation_ids().len())
            .field("blind_holdout_count", &self.split.blind_holdout_len())
            .field("blind_release_present", &self.blind_release.is_some())
            .finish_non_exhaustive()
    }
}

impl<'a> CaseSourceBundle<'a> {
    #[must_use]
    pub const fn new(experiment: &'a ExperimentArtifact, split: &'a CalibrationSplit) -> Self {
        Self {
            experiment,
            split,
            blind_release: None,
        }
    }

    /// Attach the authority release required to consume sealed blind rows.
    /// Supplying a release to a non-blind case is rejected during admission.
    #[must_use]
    pub const fn with_blind_release(mut self, release: &'a BlindReleaseReceipt) -> Self {
        self.blind_release = Some(release);
        self
    }

    #[must_use]
    pub const fn experiment(&self) -> &'a ExperimentArtifact {
        self.experiment
    }

    #[must_use]
    pub const fn split(&self) -> &'a CalibrationSplit {
        self.split
    }

    #[must_use]
    pub const fn blind_release(&self) -> Option<&'a BlindReleaseReceipt> {
        self.blind_release
    }
}

/// Concrete and opaque artifacts required to resolve a problem document.
/// Extra, missing, unverified, stale-kind, or stale-version resolutions refuse.
pub struct ProblemSourceBundle<'a> {
    context: &'a ContextOfUse,
    material: &'a MaterialCard,
    model: &'a ConstitutiveModelCard,
    cases: BTreeMap<CaseId, CaseSourceBundle<'a>>,
    opaque: SourceResolutionSet,
    concrete_authority: BTreeMap<SourceKey, AuthorityDisposition>,
}

impl fmt::Debug for ProblemSourceBundle<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ProblemSourceBundle")
            .field("case_count", &self.cases.len())
            .field("opaque_source_count", &self.opaque.entries.len())
            .field(
                "concrete_authority_count",
                &self.concrete_authority.len(),
            )
            // Concrete artifacts and resolutions remain intentionally
            // inspectable through explicit APIs. Debug is never that API:
            // recursively formatting them would disclose row/source
            // capabilities, custody links, preprocessing identities, and
            // external trust material into logs and panic reports.
            .finish_non_exhaustive()
    }
}

impl<'a> ProblemSourceBundle<'a> {
    #[must_use]
    pub fn new(
        context: &'a ContextOfUse,
        material: &'a MaterialCard,
        model: &'a ConstitutiveModelCard,
        cases: BTreeMap<CaseId, CaseSourceBundle<'a>>,
        opaque: SourceResolutionSet,
    ) -> Self {
        Self {
            context,
            material,
            model,
            cases,
            opaque,
            concrete_authority: BTreeMap::new(),
        }
    }

    /// Attach external trust-policy dispositions to concrete sources. Missing
    /// entries remain honestly `ContentVerified`; duplicate or malformed
    /// dispositions refuse.
    pub fn with_concrete_authority(
        mut self,
        entries: Vec<(SourceKey, AuthorityDisposition)>,
    ) -> Result<Self, IdentifiabilityError> {
        if self
            .concrete_authority
            .len()
            .checked_add(entries.len())
            .is_none_or(|total| total > MAX_IDENTIFIABILITY_ITEMS)
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "concrete source authority",
                detail: "too many concrete source authority entries".to_string(),
            });
        }
        for (key, authority) in entries {
            validate_authority_disposition(&authority)?;
            if self
                .concrete_authority
                .insert(key.clone(), authority)
                .is_some()
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "concrete source authority",
                    id: key.to_string(),
                });
            }
        }
        Ok(self)
    }
}

/// Source-resolved problem.  Its retained document remains inspectable, while
/// all derived bindings are read-only and recomputable from the source bundle.
#[derive(Clone, PartialEq)]
pub struct AdmittedIdentifiabilityProblem {
    document: IdentifiabilityProblemDocument,
    problem_id: ProblemId,
    source_admission_id: SourceAdmissionId,
    context: ContextBinding,
    model: MaterialModelBinding,
    data: BTreeMap<CaseId, DataLineage>,
    source_admission: SourceAdmissionRecord,
}

impl fmt::Debug for AdmittedIdentifiabilityProblem {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AdmittedIdentifiabilityProblem")
            .field("stage", &"source-admitted")
            .field("case_count", &self.document.cases.len())
            .field("parameter_count", &self.document.parameters.len())
            .field("source_count", &self.source_admission.resolutions.len())
            .finish_non_exhaustive()
    }
}

fn concrete_resolution(
    reference: &SourceRef,
    actual_kind: SourceKind,
    actual_hash: ContentHash,
    subject_artifact: Option<&ArtifactId>,
    authority: AuthorityDisposition,
) -> Result<SourceResolution, IdentifiabilityError> {
    validate_authority_disposition(&authority)?;
    validate_authority_subject_with_artifact(reference, subject_artifact, &authority)?;
    if matches!(&authority, AuthorityDisposition::Unverified { .. }) {
        return Err(IdentifiabilityError::InvalidText {
            field: "source authority",
            detail: format!(
                "concrete source {} is unresolved and cannot mint ProblemId",
                reference.key
            ),
        });
    }
    if reference.kind != actual_kind
        || reference.expected_hash != actual_hash
        || !hash_is_nonzero(actual_hash)
    {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "concrete source",
        });
    }
    let (expected_domain, expected_version) = match actual_kind {
        SourceKind::ContextOfUse
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit => (VV_ARTIFACT_SOURCE_DOMAIN, VV_SCHEMA_VERSION),
        SourceKind::MaterialCard => (MATERIAL_CARD_SOURCE_DOMAIN, MATDB_SCHEMA_VERSION),
        SourceKind::ConstitutiveModelCard => {
            (CONSTITUTIVE_MODEL_CARD_SOURCE_DOMAIN, MATDB_SCHEMA_VERSION)
        }
        _ => {
            return Err(IdentifiabilityError::InvalidText {
                field: "typed source contract",
                detail: format!("source kind {actual_kind:?} has no typed resolver in fs-material"),
            });
        }
    };
    if reference.content_hash_domain != expected_domain
        || reference.contract_version != expected_version
    {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "typed source digest domain/contract version",
        });
    }
    let resolution = SourceResolution {
        key: reference.key.clone(),
        kind: actual_kind,
        resolved_hash: actual_hash,
        content_hash_domain: reference.content_hash_domain.clone(),
        contract_version: reference.contract_version,
        authority,
        verification: SourceVerification::TypedArtifact,
    };
    validate_admitted_resolution(
        reference,
        &resolution,
        subject_artifact,
        "typed source resolution",
    )?;
    Ok(resolution)
}

fn validate_authority_disposition(
    authority: &AuthorityDisposition,
) -> Result<(), IdentifiabilityError> {
    match authority {
        AuthorityDisposition::ContentVerified => Ok(()),
        AuthorityDisposition::ExternalTrustReceipt { trust_receipt } => {
            TrustReceiptRef::try_new_with_subject_artifact(
                trust_receipt.receipt.clone(),
                trust_receipt.subject.clone(),
                trust_receipt.subject_artifact.clone(),
                trust_receipt.authentication.clone(),
            )?;
            Ok(())
        }
        AuthorityDisposition::Unverified { reason } => {
            validate_reason(reason, "unverified source reason")
        }
    }
}

fn validate_authority_subject(
    reference: &SourceRef,
    authority: &AuthorityDisposition,
) -> Result<(), IdentifiabilityError> {
    validate_authority_subject_with_artifact(reference, None, authority)
}

fn validate_authority_subject_with_artifact(
    reference: &SourceRef,
    subject_artifact: Option<&ArtifactId>,
    authority: &AuthorityDisposition,
) -> Result<(), IdentifiabilityError> {
    validate_authority_subject_fields(
        &reference.key,
        reference.kind,
        reference.expected_hash,
        &reference.content_hash_domain,
        reference.contract_version,
        authority,
    )?;
    if let AuthorityDisposition::ExternalTrustReceipt { trust_receipt } = authority {
        if let Some(declared) = trust_receipt.subject_artifact.as_ref()
            && Some(declared) != subject_artifact
        {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "trust receipt subject artifact/source resolution",
            });
        }
    }
    Ok(())
}

fn validate_authority_subject_fields(
    key: &SourceKey,
    kind: SourceKind,
    content_hash: ContentHash,
    content_hash_domain: &str,
    contract_version: u32,
    authority: &AuthorityDisposition,
) -> Result<(), IdentifiabilityError> {
    let AuthorityDisposition::ExternalTrustReceipt { trust_receipt } = authority else {
        return Ok(());
    };
    let subject = &trust_receipt.subject;
    if subject.key != *key
        || subject.kind != kind
        || subject.expected_hash != content_hash
        || subject.content_hash_domain != content_hash_domain
        || subject.contract_version != contract_version
    {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "trust receipt subject/source resolution",
        });
    }
    Ok(())
}

fn concrete_authority_for(
    bundle: &ProblemSourceBundle<'_>,
    key: &SourceKey,
) -> AuthorityDisposition {
    bundle
        .concrete_authority
        .get(key)
        .cloned()
        .unwrap_or(AuthorityDisposition::ContentVerified)
}

fn insert_exact_resolution(
    resolutions: &mut BTreeMap<SourceKey, SourceResolution>,
    key: &SourceKey,
    resolution: SourceResolution,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    if let Some(existing) = resolutions.get(key) {
        if existing != &resolution {
            return Err(IdentifiabilityError::SourceMismatch { field });
        }
    } else {
        resolutions.insert(key.clone(), resolution);
    }
    Ok(())
}

fn validate_admitted_resolution(
    reference: &SourceRef,
    resolution: &SourceResolution,
    subject_artifact: Option<&ArtifactId>,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    validate_authority_disposition(&resolution.authority)?;
    validate_authority_subject_with_artifact(reference, subject_artifact, &resolution.authority)?;
    if resolution.key != reference.key
        || resolution.kind != reference.kind
        || resolution.resolved_hash != reference.expected_hash
        || resolution.content_hash_domain != reference.content_hash_domain
        || resolution.contract_version != reference.contract_version
    {
        return Err(IdentifiabilityError::SourceMismatch { field });
    }
    if matches!(
        &resolution.authority,
        AuthorityDisposition::Unverified { .. }
    ) || matches!(&resolution.verification, SourceVerification::Unverified)
    {
        return Err(IdentifiabilityError::InvalidText {
            field: "source authority",
            detail: format!(
                "source {} is explicitly unresolved and cannot mint ProblemId",
                reference.key
            ),
        });
    }
    Ok(())
}

fn validate_discrepancy_origin(
    case_id: &CaseId,
    case: &StudyCaseDocument,
    experiment: &ExperimentArtifact,
) -> Result<(), IdentifiabilityError> {
    if matches!(&case.purpose, CasePurpose::ValidationOnly)
        && !matches!(experiment.origin(), ExperimentOrigin::Physical { .. })
    {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "validation-only case/physical experiment origin",
        });
    }
    for (channel, discrepancy) in &case.discrepancies {
        let StudyDiscrepancy::NotApplicable { basis } = discrepancy else {
            continue;
        };
        match basis {
            DiscrepancyInapplicability::PhysicalApplicability { .. } => {
                if !matches!(experiment.origin(), ExperimentOrigin::Physical { .. }) {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "physical discrepancy/experiment origin",
                    });
                }
            }
            DiscrepancyInapplicability::DeclaredSyntheticSelfModel { producer, .. } => {
                let actual_producer = match experiment.origin() {
                    ExperimentOrigin::SyntheticHighFidelity { producer } => producer,
                    ExperimentOrigin::SecondImplementation { .. }
                    | ExperimentOrigin::Physical { .. } => {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "declared-synthetic discrepancy/experiment origin",
                        });
                    }
                };
                if actual_producer != producer {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "declared-synthetic producer/forward-model binding",
                    });
                }
            }
            DiscrepancyInapplicability::ProspectiveDesign { .. } => {
                return Err(IdentifiabilityError::InvalidText {
                    field: "prospective discrepancy/experiment origin",
                    detail: format!(
                        "retrospective case {case_id} channel {channel} cannot use a prospective-design inapplicability basis"
                    ),
                });
            }
        }
    }
    Ok(())
}

fn admit_opaque_resolution(
    reference: &SourceRef,
    resolution: &SourceResolution,
) -> Result<(), IdentifiabilityError> {
    validate_admitted_resolution(reference, resolution, None, "opaque source resolution")
}

fn bind_source_reference(
    references: &mut BTreeMap<SourceKey, SourceRef>,
    source: &SourceRef,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    if let Some(prior) = references.insert(source.key.clone(), source.clone())
        && prior != *source
    {
        return Err(IdentifiabilityError::SourceMismatch { field });
    }
    Ok(())
}

fn validate_source_authority_closure(
    references: &BTreeMap<SourceKey, SourceRef>,
    authority: &SourceResolutionSet,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    if authority.entries.len() != references.len() {
        return Err(IdentifiabilityError::Cardinality {
            field,
            detail: "every referenced source needs exactly one locally verified resolution"
                .to_string(),
        });
    }
    for (key, reference) in references {
        let resolution =
            authority
                .entries
                .get(key)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field,
                    id: key.to_string(),
                })?;
        admit_opaque_resolution(reference, resolution)?;
    }
    Ok(())
}

fn encode_source_admission(
    admission: &SourceAdmissionRecord,
) -> Result<Vec<u8>, IdentifiabilityError> {
    check_source_admission_identity_version(admission.schema_version)?;
    let mut writer = CanonicalWriter::new();
    writer.raw(SOURCE_ADMISSION_MAGIC);
    writer.u32(admission.schema_version);
    writer.hash(admission.problem_id.0);
    writer.count(admission.resolutions.len(), "source admission resolutions")?;
    for resolution in admission.resolutions.values() {
        encode_resolution_entry(&mut writer, resolution)?;
    }
    writer.finish()
}

fn source_admission_identity_hash(
    admission: &SourceAdmissionRecord,
) -> Result<SourceAdmissionId, IdentifiabilityError> {
    Ok(SourceAdmissionId(hash_domain(
        IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN,
        &encode_source_admission(admission)?,
    )))
}

fn problem_identity_hash(
    document: &IdentifiabilityProblemDocument,
) -> Result<ProblemId, IdentifiabilityError> {
    Ok(ProblemId(hash_domain(
        IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN,
        &encode_problem(document)?,
    )))
}

impl AdmittedIdentifiabilityProblem {
    /// Resolve exact concrete sources, require a closed authority set for every
    /// opaque reference, re-derive V&V bindings/lineage, and only then mint
    /// problem and source-admission identities.
    pub fn resolve_and_admit(
        document: IdentifiabilityProblemDocument,
        bundle: ProblemSourceBundle<'_>,
    ) -> Result<Self, IdentifiabilityError> {
        if bundle.cases.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "retrospective case source bundles",
                detail: "too many concrete case source bundles".to_string(),
            });
        }
        if bundle.concrete_authority.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "concrete source authority",
                detail: "too many concrete source authority entries".to_string(),
            });
        }
        let context = ContextBinding::from_vv(bundle.context)?;
        let graph_hash = document.sources[&document.graph_source].expected_hash;
        let model = MaterialModelBinding::from_cards(bundle.material, bundle.model, graph_hash)?;
        let mut resolutions = BTreeMap::new();
        let mut concrete_keys = BTreeSet::new();

        let context_ref = &document.sources[&document.context_source];
        let context_resolution = concrete_resolution(
            context_ref,
            SourceKind::ContextOfUse,
            context.reference.hash(),
            None,
            concrete_authority_for(&bundle, &context_ref.key),
        )?;
        concrete_keys.insert(context_ref.key.clone());
        resolutions.insert(context_ref.key.clone(), context_resolution);

        let material_ref = &document.sources[&document.material_source];
        let material_resolution = concrete_resolution(
            material_ref,
            SourceKind::MaterialCard,
            bundle.material.content_hash(),
            None,
            concrete_authority_for(&bundle, &material_ref.key),
        )?;
        concrete_keys.insert(material_ref.key.clone());
        resolutions.insert(material_ref.key.clone(), material_resolution);

        let model_ref = &document.sources[&document.model_source];
        let model_resolution = concrete_resolution(
            model_ref,
            SourceKind::ConstitutiveModelCard,
            bundle.model.content_hash(),
            None,
            concrete_authority_for(&bundle, &model_ref.key),
        )?;
        concrete_keys.insert(model_ref.key.clone());
        resolutions.insert(model_ref.key.clone(), model_resolution);

        // The physical parameter roster is closed against the exact model card.
        for (role, roster) in &model.parameter_roster {
            let parameter = document.parameters.get(role).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "model-card parameter declaration",
                    id: role.to_string(),
                }
            })?;
            if !matches!(&parameter.owner, ParameterOwnerBinding::ConstitutiveModel)
                || parameter.quantity != roster.quantity
                || roster.nominal() < parameter.domain.lo
                || roster.nominal() > parameter.domain.hi
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "model-card parameter binding",
                    detail: format!(
                        "parameter {role} must match owner, exact quantity, and nominal domain"
                    ),
                });
            }
        }
        for parameter in document.parameters.values() {
            if matches!(&parameter.owner, ParameterOwnerBinding::ConstitutiveModel)
                && !model.parameter_roster.contains_key(&parameter.role)
            {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "constitutive-model parameter",
                    id: parameter.role.to_string(),
                });
            }
        }
        if matches!(
            model.initial_state_policy,
            InitialStatePolicy::ZeroInternalState
        ) && document
            .parameters
            .values()
            .any(|parameter| matches!(&parameter.owner, ParameterOwnerBinding::InitialState { .. }))
        {
            return Err(IdentifiabilityError::InitialStatePolicy {
                detail:
                    "zero-internal-state model cannot expose inferential initial-state parameters"
                        .to_string(),
            });
        }

        // A blind release is authority over a split, not an observation-local
        // annotation. Pre-scan by split key so shared uses cannot acquire
        // order-dependent or contradictory authority dispositions.
        let mut blind_releases = BTreeMap::<SourceKey, &BlindReleaseReceipt>::new();
        for (case_id, case) in &document.cases {
            let CaseDataDeclaration::Retrospective { split, .. } = &case.data else {
                continue;
            };
            let case_sources = bundle.cases.get(case_id).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "retrospective case source bundle",
                    id: case_id.to_string(),
                }
            })?;
            match (&case.purpose, case_sources.blind_release) {
                (CasePurpose::BlindFalsification, Some(release)) => {
                    if let Some(existing) = blind_releases.insert(split.clone(), release)
                        && existing != release
                    {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "shared split blind release",
                        });
                    }
                }
                (CasePurpose::BlindFalsification, None) => {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "blind release",
                        detail: format!(
                            "blind-falsification case {case_id} requires an authority release"
                        ),
                    });
                }
                (_, Some(_)) => {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "blind release",
                        detail: format!(
                            "non-blind case {case_id} must not receive blind-release authority"
                        ),
                    });
                }
                (_, None) => {}
            }
        }

        let mut data = BTreeMap::new();
        for (case_id, case) in &document.cases {
            case.initial_state.validate_against(&model)?;
            if case.protocol.state_schema_version != model.state_schema_version {
                return Err(IdentifiabilityError::VersionMismatch {
                    field: "case protocol/model state schema",
                    expected: model.state_schema_version,
                    actual: case.protocol.state_schema_version,
                });
            }
            for observation in case.observations.values() {
                if context.qoi_units.get(&observation.qoi) != Some(&observation.unit) {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "context QoI/unit",
                        id: observation.qoi.as_str().to_string(),
                    });
                }
            }
            match &case.data {
                CaseDataDeclaration::Prospective => {
                    if bundle.cases.contains_key(case_id) {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "prospective case source bundle",
                            detail: format!(
                                "prospective case {case_id} must not receive experiment data"
                            ),
                        });
                    }
                }
                CaseDataDeclaration::Retrospective {
                    experiment,
                    split,
                    parser,
                    preprocessing,
                    parser_version,
                    split_grouping,
                } => {
                    let case_sources = bundle.cases.get(case_id).ok_or_else(|| {
                        IdentifiabilityError::UnknownReference {
                            field: "retrospective case source bundle",
                            id: case_id.to_string(),
                        }
                    })?;
                    validate_discrepancy_origin(case_id, case, case_sources.experiment)?;
                    let experiment_hash =
                        case_sources.experiment.content_hash().map_err(|error| {
                            IdentifiabilityError::Vv {
                                detail: error.to_string(),
                            }
                        })?;
                    let split_hash = case_sources.split.content_hash().map_err(|error| {
                        IdentifiabilityError::Vv {
                            detail: error.to_string(),
                        }
                    })?;
                    let experiment_resolution = concrete_resolution(
                        &document.sources[experiment],
                        SourceKind::ExperimentArtifact,
                        experiment_hash,
                        None,
                        concrete_authority_for(&bundle, experiment),
                    )?;
                    concrete_keys.insert(experiment.clone());
                    insert_exact_resolution(
                        &mut resolutions,
                        experiment,
                        experiment_resolution,
                        "shared experiment source resolution",
                    )?;

                    let release_authority = blind_releases
                        .get(split)
                        .map(|release| {
                            Ok(AuthorityDisposition::ExternalTrustReceipt {
                                trust_receipt: TrustReceiptRef::blind_release(
                                    &document.sources[split],
                                    release.split().id().clone(),
                                    release.authority_receipt_hash(),
                                )?,
                            })
                        })
                        .transpose()?;
                    if let (Some(required), Some(declared)) = (
                        release_authority.as_ref(),
                        bundle.concrete_authority.get(split),
                    ) && required != declared
                    {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "blind release/concrete source authority",
                        });
                    }
                    let split_authority =
                        release_authority.unwrap_or_else(|| concrete_authority_for(&bundle, split));
                    let split_resolution = concrete_resolution(
                        &document.sources[split],
                        SourceKind::CalibrationSplit,
                        split_hash,
                        Some(case_sources.split.id()),
                        split_authority,
                    )?;
                    concrete_keys.insert(split.clone());
                    insert_exact_resolution(
                        &mut resolutions,
                        split,
                        split_resolution,
                        "shared split source resolution",
                    )?;
                    let parser_hash = document.sources[parser].expected_hash;
                    let preprocessing_hash = document.sources[preprocessing].expected_hash;
                    let lineage = DataLineage::from_vv(
                        case_sources.experiment,
                        case_sources.split,
                        parser_hash,
                        *parser_version,
                        preprocessing_hash,
                        split_grouping.clone(),
                    )?;
                    for observation in case.observations.values() {
                        if !lineage.qois().contains(&observation.qoi) {
                            return Err(IdentifiabilityError::UnknownReference {
                                field: "experiment observation QoI",
                                id: observation.qoi.as_str().to_string(),
                            });
                        }
                        let instrument = case_sources
                            .experiment
                            .instrument_calibration(&observation.instrument)
                            .ok_or_else(|| IdentifiabilityError::UnknownReference {
                                field: "experiment observation instrument",
                                id: observation.instrument.as_str().to_string(),
                            })?;
                        if document.sources[&observation.sensor].expected_hash
                            != instrument.certificate_hash()
                        {
                            return Err(IdentifiabilityError::SourceMismatch {
                                field: "observation sensor/instrument calibration",
                            });
                        }
                        if !case_sources.experiment.contains_clock(&observation.clock) {
                            return Err(IdentifiabilityError::UnknownReference {
                                field: "experiment observation clock",
                                id: observation.clock.as_str().to_string(),
                            });
                        }
                    }
                    let allowed_rows: BTreeSet<ObservationId> = match &case.purpose {
                        CasePurpose::Calibration
                        | CasePurpose::SymmetryBreaking
                        | CasePurpose::Complementary { .. } => lineage.calibration_ids.clone(),
                        CasePurpose::ValidationOnly => lineage.validation_ids.clone(),
                        CasePurpose::BlindFalsification => {
                            lineage.blind_sources.keys().cloned().collect()
                        }
                        CasePurpose::ProspectiveDesign => BTreeSet::new(),
                    };
                    let declared_rows = case
                        .observations
                        .values()
                        .filter_map(|observation| match &observation.rows {
                            ObservationRows::Retrospective(rows) => Some(rows.iter().cloned()),
                            ObservationRows::Prospective => None,
                        })
                        .flatten()
                        .collect::<BTreeSet<_>>();
                    for observation in case.observations.values() {
                        let ObservationRows::Retrospective(rows) = &observation.rows else {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "retrospective observation rows",
                                detail: format!(
                                    "case {case_id} contains a prospective observation after structural admission"
                                ),
                            });
                        };
                        if !rows.is_subset(&lineage.observation_ids) {
                            return Err(IdentifiabilityError::UnknownReference {
                                field: "observation raw row",
                                id: observation.id.to_string(),
                            });
                        }
                        if !rows.is_subset(&allowed_rows) {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "case-purpose data partition",
                                detail: format!(
                                    "observation {} consumes rows outside the partition authorized by case {case_id} purpose",
                                    observation.id
                                ),
                            });
                        }
                        for row in rows {
                            let binding =
                                case_sources.experiment.manifest().row(row).ok_or_else(|| {
                                    IdentifiabilityError::UnknownReference {
                                        field: "experiment manifest row binding",
                                        id: row.as_str().to_string(),
                                    }
                                })?;
                            if binding.qoi() != &observation.qoi
                                || binding.instrument() != &observation.instrument
                                || binding.acquisition_channel() != &observation.acquisition_channel
                                || binding.clock() != &observation.clock
                            {
                                return Err(IdentifiabilityError::SourceMismatch {
                                    field: "observation/manifest row binding",
                                });
                            }
                        }
                    }
                    if matches!(&case.purpose, CasePurpose::BlindFalsification) {
                        let release = blind_releases
                            .get(split)
                            .expect("blind release pre-scan established exact presence");
                        let split_reference = ArtifactRef::new(
                            ArtifactKind::CalibrationSplit,
                            case_sources.split.id().clone(),
                            split_hash,
                        );
                        case_sources
                            .split
                            .blind_selection(
                                split_reference,
                                declared_rows.iter().cloned().collect(),
                                (**release).clone(),
                            )
                            .map_err(|error| IdentifiabilityError::Vv {
                                detail: error.to_string(),
                            })?;
                    }
                    data.insert(case_id.clone(), lineage);
                }
            }
        }
        if bundle.cases.len() != data.len() {
            return Err(IdentifiabilityError::Cardinality {
                field: "case source bundles",
                detail: "source bundle contains an unknown or prospective case".to_string(),
            });
        }

        let mut source_bytes_owners = BTreeMap::<ContentHash, BTreeSet<CaseId>>::new();
        let mut manifest_owners = BTreeMap::<ContentHash, BTreeSet<CaseId>>::new();
        let mut row_locator_owners =
            BTreeMap::<ObservationLocatorIdentity, BTreeSet<CaseId>>::new();
        for (case_id, lineage) in &data {
            source_bytes_owners
                .entry(lineage.source_bytes())
                .or_default()
                .insert(case_id.clone());
            manifest_owners
                .entry(lineage.raw_manifest())
                .or_default()
                .insert(case_id.clone());
            for row in lineage.row_bindings().values() {
                row_locator_owners
                    .entry(row.source_ref().locator_identity())
                    .or_default()
                    .insert(case_id.clone());
            }
        }
        fn admit_shared_owner_sets<K: Ord>(
            owners: &BTreeMap<K, BTreeSet<CaseId>>,
            policy: &DataReusePolicy,
            sharing_participation: &mut BTreeSet<CaseId>,
        ) -> Result<(), IdentifiabilityError> {
            for cases in owners.values().filter(|cases| cases.len() > 1) {
                let mut memberships = cases
                    .iter()
                    .map(|case| sharing_group_membership(policy, case));
                let first = memberships.next().expect("owner set is nonempty");
                if first.is_none() || memberships.any(|membership| membership != first) {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "data reuse policy",
                        detail: format!(
                            "cases {} share admitted provenance without one exact joint sharing group",
                            cases
                                .iter()
                                .map(CaseId::as_str)
                                .collect::<Vec<_>>()
                                .join(",")
                        ),
                    });
                }
                sharing_participation.extend(cases.iter().cloned());
            }
            Ok(())
        }
        let mut sharing_participation = BTreeSet::<CaseId>::new();
        admit_shared_owner_sets(
            &source_bytes_owners,
            &document.data_reuse,
            &mut sharing_participation,
        )?;
        admit_shared_owner_sets(
            &manifest_owners,
            &document.data_reuse,
            &mut sharing_participation,
        )?;
        admit_shared_owner_sets(
            &row_locator_owners,
            &document.data_reuse,
            &mut sharing_participation,
        )?;
        if let DataReusePolicy::Shared { groups } = &document.data_reuse {
            for group in groups {
                for case_id in &group.cases {
                    if !sharing_participation.contains(case_id) {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "data sharing group",
                            detail: format!(
                                "case {case_id} declares raw-data sharing but overlaps no peer by admitted bytes, manifest, or row source"
                            ),
                        });
                    }
                }
            }
        }

        for (key, reference) in &document.sources {
            if concrete_keys.contains(key) {
                if bundle.opaque.entries.contains_key(key) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "concrete/opaque source resolution",
                        id: key.to_string(),
                    });
                }
                continue;
            }
            let resolution = bundle.opaque.entries.get(key).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "opaque source resolution",
                    id: key.to_string(),
                }
            })?;
            admit_opaque_resolution(reference, resolution)?;
            resolutions.insert(key.clone(), resolution.clone());
        }
        if bundle.opaque.entries.len() != document.sources.len() - concrete_keys.len()
            || resolutions.len() != document.sources.len()
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "source resolution closure",
                detail: "resolution set has missing or extra source keys".to_string(),
            });
        }
        if let Some(key) = bundle
            .concrete_authority
            .keys()
            .find(|key| !concrete_keys.contains(*key))
        {
            return Err(IdentifiabilityError::UnknownReference {
                field: "concrete source authority",
                id: key.to_string(),
            });
        }

        let problem_id = problem_identity_hash(&document)?;
        let source_admission = SourceAdmissionRecord {
            schema_version: IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION,
            problem_id,
            resolutions,
        };
        Ok(Self {
            problem_id,
            source_admission_id: source_admission_identity_hash(&source_admission)?,
            document,
            context,
            model,
            data,
            source_admission,
        })
    }

    #[must_use]
    pub const fn id(&self) -> ProblemId {
        self.problem_id
    }

    #[must_use]
    pub const fn source_admission_id(&self) -> SourceAdmissionId {
        self.source_admission_id
    }

    #[must_use]
    pub const fn document(&self) -> &IdentifiabilityProblemDocument {
        &self.document
    }

    #[must_use]
    pub const fn data(&self) -> &BTreeMap<CaseId, DataLineage> {
        &self.data
    }

    #[must_use]
    pub const fn context(&self) -> &ContextBinding {
        &self.context
    }

    #[must_use]
    pub const fn model(&self) -> &MaterialModelBinding {
        &self.model
    }

    #[must_use]
    pub const fn source_resolutions(&self) -> &BTreeMap<SourceKey, SourceResolution> {
        &self.source_admission.resolutions
    }

    /// Exact source-admission identity preimage retained for ledger audit.
    pub fn source_admission_canonical_bytes(&self) -> Result<Vec<u8>, IdentifiabilityError> {
        encode_source_admission(&self.source_admission)
    }
}

/// One explicit action for every physical parameter.  Conditioned and derived
/// parameters remain explicit so a plan cannot silently omit them.
#[derive(Debug, Clone, PartialEq)]
pub enum ParameterExecutionAction {
    Optimize {
        coordinate: ParameterCoordinate,
    },
    Profile {
        coordinate: ParameterCoordinate,
    },
    Marginalize {
        coordinate: ParameterCoordinate,
        integrator: SourceRef,
        /// Exact physical-prior to coordinate-measure transport, including
        /// Jacobian, truncation, and normalization semantics.
        measure_transport: SourceRef,
    },
    Conditioned,
    Derived,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RequestedClaimAxis {
    Structural,
    Local,
    Generic,
    Global,
    Practical,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArithmeticPolicy {
    ExactSymbolic,
    CertifiedInterval,
    DeterministicFloatingPoint,
    FastFloatingPoint,
}

/// Exact dimensionless semantics under which a numerical claim error is
/// interpreted. The metric and nondimensionalization are separately sourced:
/// equal scalar bounds under different policies are not interchangeable.
#[derive(Debug, Clone, PartialEq)]
pub struct DimensionlessErrorPolicy {
    metric: SourceRef,
    nondimensionalization: SourceRef,
    maximum_certified_error: f64,
}

impl DimensionlessErrorPolicy {
    pub fn try_new(
        metric: SourceRef,
        nondimensionalization: SourceRef,
        maximum_certified_error: f64,
    ) -> Result<Self, IdentifiabilityError> {
        if metric.kind != SourceKind::DimensionlessErrorMetric {
            return Err(IdentifiabilityError::InvalidText {
                field: "claim error metric",
                detail: "claim error metric must use DimensionlessErrorMetric semantics"
                    .to_string(),
            });
        }
        if nondimensionalization.kind != SourceKind::Nondimensionalization {
            return Err(IdentifiabilityError::InvalidText {
                field: "claim nondimensionalization",
                detail: "claim scale must use Nondimensionalization semantics".to_string(),
            });
        }
        if !maximum_certified_error.is_finite() || maximum_certified_error < 0.0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "maximum certified claim error",
                detail: "bound must be finite and nonnegative".to_string(),
            });
        }
        Ok(Self {
            metric,
            nondimensionalization,
            maximum_certified_error: canonical_f64(maximum_certified_error),
        })
    }

    #[must_use]
    pub const fn metric(&self) -> &SourceRef {
        &self.metric
    }

    #[must_use]
    pub const fn nondimensionalization(&self) -> &SourceRef {
        &self.nondimensionalization
    }

    #[must_use]
    pub const fn maximum_certified_error(&self) -> f64 {
        self.maximum_certified_error
    }
}

/// Confirmatory proposition preregistered before execution together with its
/// exact acceptance metric. Assessments may not substitute another proposition
/// merely because it shares a coarse local/global/structural label.
#[derive(Debug, Clone, PartialEq)]
pub struct ClaimRequest {
    claim: TypedIdentifiabilityClaim,
    error_policy: DimensionlessErrorPolicy,
}

impl ClaimRequest {
    #[must_use]
    pub const fn new(
        claim: TypedIdentifiabilityClaim,
        error_policy: DimensionlessErrorPolicy,
    ) -> Self {
        Self {
            claim,
            error_policy,
        }
    }

    #[must_use]
    pub const fn claim(&self) -> &TypedIdentifiabilityClaim {
        &self.claim
    }

    #[must_use]
    pub const fn error_policy(&self) -> &DimensionlessErrorPolicy {
        &self.error_policy
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IdentifiabilityNumericalPolicy {
    rank_tolerance: f64,
    singular_value_floor: f64,
    maximum_condition_number: f64,
    arithmetic: ArithmeticPolicy,
    nondimensionalization: SourceRef,
}

impl IdentifiabilityNumericalPolicy {
    pub fn try_new(
        rank_tolerance: f64,
        singular_value_floor: f64,
        maximum_condition_number: f64,
        arithmetic: ArithmeticPolicy,
        nondimensionalization: SourceRef,
    ) -> Result<Self, IdentifiabilityError> {
        if !rank_tolerance.is_finite()
            || rank_tolerance <= 0.0
            || !singular_value_floor.is_finite()
            || singular_value_floor < 0.0
            || !maximum_condition_number.is_finite()
            || maximum_condition_number < 1.0
        {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "identifiability numerical policy",
                detail: "tolerances must be finite and physically ordered".to_string(),
            });
        }
        if nondimensionalization.kind != SourceKind::Nondimensionalization {
            return Err(IdentifiabilityError::InvalidText {
                field: "numerical nondimensionalization",
                detail: "rank and singular-value controls require exact dimensionless scaling"
                    .to_string(),
            });
        }
        Ok(Self {
            rank_tolerance: canonical_f64(rank_tolerance),
            singular_value_floor: canonical_f64(singular_value_floor),
            maximum_condition_number: canonical_f64(maximum_condition_number),
            arithmetic,
            nondimensionalization,
        })
    }

    #[must_use]
    pub const fn rank_tolerance(&self) -> f64 {
        self.rank_tolerance
    }

    #[must_use]
    pub const fn singular_value_floor(&self) -> f64 {
        self.singular_value_floor
    }

    #[must_use]
    pub const fn maximum_condition_number(&self) -> f64 {
        self.maximum_condition_number
    }

    #[must_use]
    pub const fn arithmetic(&self) -> ArithmeticPolicy {
        self.arithmetic
    }

    #[must_use]
    pub const fn nondimensionalization(&self) -> &SourceRef {
        &self.nondimensionalization
    }
}

/// Numerical configuration whose identity is deliberately separate from the
/// physical problem.
#[derive(Clone, PartialEq)]
pub struct IdentifiabilityExecutionPlan {
    schema_version: u32,
    header: ArtifactHeader,
    problem_id: ProblemId,
    source_admission_id: SourceAdmissionId,
    analyzer: SourceRef,
    build: SourceRef,
    derivative_provider: Option<SourceRef>,
    claim_requests: BTreeMap<ClaimId, ClaimRequest>,
    actions: BTreeMap<ParameterRoleId, ParameterExecutionAction>,
    gauge_reductions: BTreeMap<GaugeReductionId, GaugeReductionBinding>,
    numerical: IdentifiabilityNumericalPolicy,
    initialization: SourceRef,
    stopping: SourceRef,
    determinism_contract: SourceRef,
    source_authority: SourceResolutionSet,
}

impl fmt::Debug for IdentifiabilityExecutionPlan {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IdentifiabilityExecutionPlan")
            .field("stage", &"execution-plan")
            .field("claim_count", &self.claim_requests.len())
            .field("action_count", &self.actions.len())
            .field("gauge_reduction_count", &self.gauge_reductions.len())
            .field("source_count", &self.source_authority.entries.len())
            .finish_non_exhaustive()
    }
}

fn add_claim_structural_items(
    budget: &mut StructuralItemBudget,
    claim: &TypedIdentifiabilityClaim,
) -> Result<(), IdentifiabilityError> {
    match &claim.subject {
        ClaimSubject::ParameterSet(roles)
        | ClaimSubject::DerivedFunctional {
            parameters: roles, ..
        } => budget.add(roles.len(), "claim parameter support")?,
        ClaimSubject::Parameter(_)
        | ClaimSubject::Influence(_)
        | ClaimSubject::GaugeAction(_)
        | ClaimSubject::WholeProblem => {}
    }
    match &claim.scope {
        ClaimScope::Cases(cases) | ClaimScope::Stratum { cases, .. } => {
            budget.add(cases.len(), "claim case scope")?;
        }
        ClaimScope::WholeCampaign => {}
    }
    Ok(())
}

fn add_reduction_structural_items(
    budget: &mut StructuralItemBudget,
    binding: &GaugeReductionBinding,
) -> Result<(), IdentifiabilityError> {
    budget.add(binding.claims.len(), "gauge reduction claims")?;
    if let GaugeReductionStage::After { predecessors, .. } = &binding.stage {
        budget.add(predecessors.len(), "gauge reduction predecessors")?;
    }
    let slice = match &binding.plan {
        GaugeReductionPlan::Slice { slice }
        | GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction: ContinuousGaugeReductionPlan::Slice { slice },
            ..
        } => Some(slice),
        GaugeReductionPlan::Unreduced { .. }
        | GaugeReductionPlan::Quotient { .. }
        | GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction: ContinuousGaugeReductionPlan::Quotient { .. },
            ..
        } => None,
    };
    if let Some(slice) = slice {
        budget.add(slice.support.len(), "gauge slice support")?;
    }
    Ok(())
}

fn validate_execution_structural_budget(
    plan: &IdentifiabilityExecutionPlan,
) -> Result<(), IdentifiabilityError> {
    let mut budget = StructuralItemBudget::default();
    budget.add(plan.claim_requests.len(), "execution claim requests")?;
    for request in plan.claim_requests.values() {
        add_claim_structural_items(&mut budget, &request.claim)?;
    }
    budget.add(plan.actions.len(), "execution parameter actions")?;
    budget.add(plan.gauge_reductions.len(), "execution gauge reductions")?;
    for binding in plan.gauge_reductions.values() {
        add_reduction_structural_items(&mut budget, binding)?;
    }
    budget.add(
        plan.source_authority.entries.len(),
        "execution source authority",
    )?;
    Ok(())
}

fn validate_coordinate_for_parameter(
    parameter: &StudyParameter,
    coordinate: &ParameterCoordinate,
) -> Result<(), IdentifiabilityError> {
    match coordinate.transform {
        CoordinateTransform::Identity => {
            if coordinate.quantity != parameter.quantity {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "execution coordinate quantity",
                    detail: format!(
                        "identity coordinate for {} must preserve exact QuantitySpec",
                        parameter.role
                    ),
                });
            }
        }
        CoordinateTransform::Affine { scale_quantity, .. } => {
            let mapped = checked_add_dims(coordinate.quantity.dims(), scale_quantity.dims())
                .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                    field: "execution affine coordinate",
                    detail: "dimension exponent overflow".to_string(),
                })?;
            if mapped != parameter.quantity.dims() {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "execution affine coordinate",
                    detail: format!("coordinate for {} has wrong dimensions", parameter.role),
                });
            }
        }
        CoordinateTransform::LogPositive { .. } => {
            if coordinate.quantity.dims() != Dims([0; 6]) || parameter.domain.lo <= 0.0 {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "execution log coordinate",
                    detail: format!(
                        "parameter {} is not positive/dimensionless-charted",
                        parameter.role
                    ),
                });
            }
        }
    }
    let mapped = coordinate.transform.mapped_domain(coordinate.domain)?;
    if !same_f64(mapped.lo, parameter.domain.lo) || !same_f64(mapped.hi, parameter.domain.hi) {
        return Err(IdentifiabilityError::InvalidNumeric {
            field: "execution coordinate domain",
            detail: format!(
                "coordinate does not bijectively cover parameter {}",
                parameter.role
            ),
        });
    }
    Ok(())
}

impl IdentifiabilityExecutionPlan {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        header: ArtifactHeader,
        problem: &AdmittedIdentifiabilityProblem,
        analyzer: SourceRef,
        build: SourceRef,
        derivative_provider: Option<SourceRef>,
        claim_requests: Vec<ClaimRequest>,
        actions: Vec<(ParameterRoleId, ParameterExecutionAction)>,
        gauge_reductions: Vec<GaugeReductionBinding>,
        numerical: IdentifiabilityNumericalPolicy,
        initialization: SourceRef,
        stopping: SourceRef,
        determinism_contract: SourceRef,
        source_authority: SourceResolutionSet,
    ) -> Result<Self, IdentifiabilityError> {
        validate_header_profile(&header)?;
        if !header.capabilities().contains("identifiability.execute") {
            return Err(IdentifiabilityError::InvalidText {
                field: "execution capability",
                detail: "missing identifiability.execute capability".to_string(),
            });
        }
        if claim_requests.is_empty() || claim_requests.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "claim requests",
                detail: "execution needs bounded nonempty exact claim preregistration".to_string(),
            });
        }
        if actions.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "execution parameter actions",
                detail: "parameter action input exceeds the canonical collection bound".to_string(),
            });
        }
        if gauge_reductions.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "execution gauge reductions",
                detail: "gauge reduction input exceeds the canonical collection bound".to_string(),
            });
        }
        let mut execution_sources = BTreeMap::new();
        for (source, kind, field) in [
            (&analyzer, SourceKind::Analyzer, "analyzer"),
            (&build, SourceKind::Build, "build"),
            (&initialization, SourceKind::Assumption, "initialization"),
            (&stopping, SourceKind::Assumption, "stopping policy"),
            (
                &determinism_contract,
                SourceKind::Assumption,
                "determinism contract",
            ),
        ] {
            if source.kind != kind {
                return Err(IdentifiabilityError::InvalidText {
                    field,
                    detail: format!("source {} has wrong kind", source.key),
                });
            }
            bind_source_reference(&mut execution_sources, source, "execution source alias")?;
        }
        if let Some(provider) = &derivative_provider {
            if provider.kind != SourceKind::DerivativeProvider {
                return Err(IdentifiabilityError::InvalidText {
                    field: "derivative provider",
                    detail: "source has wrong kind".to_string(),
                });
            }
            bind_source_reference(&mut execution_sources, provider, "execution source alias")?;
        }
        bind_source_reference(
            &mut execution_sources,
            numerical.nondimensionalization(),
            "execution source alias",
        )?;
        let claim_requests = insert_unique(claim_requests, "claim requests", |request| {
            request.claim.id()
        })?;
        for request in claim_requests.values() {
            validate_claim_compatibility(&request.claim, &problem.document)?;
            for source in validate_claim_sources(&request.claim, &problem.document)? {
                bind_source_reference(
                    &mut execution_sources,
                    source,
                    "execution claim source alias",
                )?;
            }
            for source in [
                request.error_policy.metric(),
                request.error_policy.nondimensionalization(),
            ] {
                bind_source_reference(
                    &mut execution_sources,
                    source,
                    "execution claim policy source alias",
                )?;
            }
        }
        let mut action_map = BTreeMap::new();
        for (role, action) in actions {
            if action_map.insert(role.clone(), action).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "execution parameter action",
                    id: role.to_string(),
                });
            }
        }
        if action_map.len() != problem.document.parameters.len() {
            return Err(IdentifiabilityError::Cardinality {
                field: "execution parameter actions",
                detail: "every physical parameter needs exactly one explicit action".to_string(),
            });
        }
        let mut coordinate_ids = BTreeSet::new();
        for (role, parameter) in &problem.document.parameters {
            let action =
                action_map
                    .get(role)
                    .ok_or_else(|| IdentifiabilityError::UnknownReference {
                        field: "execution parameter action",
                        id: role.to_string(),
                    })?;
            if matches!(&parameter.scope, ParameterScopeBinding::Field { .. })
                && matches!(
                    action,
                    ParameterExecutionAction::Optimize { .. }
                        | ParameterExecutionAction::Profile { .. }
                        | ParameterExecutionAction::Marginalize { .. }
                )
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "field parameter execution carrier",
                    detail: format!(
                        "field-valued parameter {role} cannot be scalarized through ParameterCoordinate; a model-space/discretization/reconstruction-bound function-space coordinate is required"
                    ),
                });
            }
            match (&parameter.treatment, action) {
                (
                    ParameterTreatment::Estimated,
                    ParameterExecutionAction::Optimize { coordinate },
                )
                | (
                    ParameterTreatment::Profiled,
                    ParameterExecutionAction::Profile { coordinate },
                ) => {
                    validate_coordinate_for_parameter(parameter, coordinate)?;
                }
                (
                    ParameterTreatment::Marginalized,
                    ParameterExecutionAction::Marginalize {
                        coordinate,
                        integrator,
                        measure_transport,
                    },
                ) => {
                    validate_coordinate_for_parameter(parameter, coordinate)?;
                    if integrator.kind != SourceKind::Analyzer {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "marginalization integrator",
                            detail: "integrator source must have Analyzer kind".to_string(),
                        });
                    }
                    if measure_transport.kind != SourceKind::MeasureTransport {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "marginalization measure transport",
                            detail: "marginalized coordinates need exact change-of-variables, Jacobian, truncation, and normalization semantics"
                                .to_string(),
                        });
                    }
                    bind_source_reference(
                        &mut execution_sources,
                        integrator,
                        "execution source alias",
                    )?;
                    bind_source_reference(
                        &mut execution_sources,
                        measure_transport,
                        "execution marginal measure source alias",
                    )?;
                }
                (ParameterTreatment::Conditioned(_), ParameterExecutionAction::Conditioned)
                | (ParameterTreatment::Derived { .. }, ParameterExecutionAction::Derived) => {}
                _ => {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "execution parameter treatment",
                        detail: format!("action for {role} contradicts physical treatment"),
                    });
                }
            }
            let coordinate = match action {
                ParameterExecutionAction::Optimize { coordinate }
                | ParameterExecutionAction::Profile { coordinate }
                | ParameterExecutionAction::Marginalize { coordinate, .. } => Some(coordinate),
                ParameterExecutionAction::Conditioned | ParameterExecutionAction::Derived => None,
            };
            if let Some(coordinate) = coordinate
                && !coordinate_ids.insert(coordinate.id().clone())
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "execution scalar coordinate",
                    id: coordinate.id().to_string(),
                });
            }
        }
        for constraint in problem.document.constraints.values() {
            let unsupported = joint_constraint_support(constraint)
                .into_iter()
                .find(|role| {
                    matches!(
                        action_map.get(role),
                        Some(
                            ParameterExecutionAction::Optimize { .. }
                                | ParameterExecutionAction::Profile { .. }
                                | ParameterExecutionAction::Marginalize { .. }
                        )
                    )
                });
            if let Some(role) = unsupported {
                return Err(IdentifiabilityError::InvalidText {
                    field: "joint-constraint execution carrier",
                    detail: format!(
                        "constraint {} contains free parameter {role}, but v3 has no identity-bearing joint chart/retraction/projection/constrained-solver or coupled-measure plan; independent scalar actions are rejected",
                        constraint.id
                    ),
                });
            }
        }
        let mut gauge_reduction_map = BTreeMap::new();
        let mut reduced_action_claims = BTreeSet::new();
        for binding in gauge_reductions {
            if gauge_reduction_map.contains_key(&binding.id) {
                return Err(IdentifiabilityError::Duplicate {
                    field: "execution gauge reduction",
                    id: binding.id.to_string(),
                });
            }
            for claim_id in &binding.claims {
                let request = claim_requests.get(claim_id).ok_or_else(|| {
                    IdentifiabilityError::UnknownReference {
                        field: "gauge reduction claim",
                        id: claim_id.to_string(),
                    }
                })?;
                let cases = claim_case_set(&request.claim, &problem.document)?;
                let view =
                    gauge_action_view(&binding.action, &request.claim, &cases, &problem.document)?;
                let (action_status, action_geometry) =
                    gauge_action_status_and_geometry(&binding.action, &problem.document)?;
                if matches!(action_status, GaugeStatus::Candidate { .. })
                    && !matches!(&binding.plan, GaugeReductionPlan::Unreduced { .. })
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "candidate gauge reduction",
                        detail: format!(
                            "reduction {} cannot quotient or slice a Candidate action as though it were established; retain it Unreduced until hypothesis evidence changes the physical declaration",
                            binding.id
                        ),
                    });
                }
                if matches!(action_geometry, GaugeOrbitGeometry::Stratified { .. })
                    && reduction_uses_regular_atlas(&binding.plan)
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "stratified gauge quotient",
                        detail: format!(
                            "reduction {} cannot use a regular principal-bundle atlas for stratified orbit geometry",
                            binding.id
                        ),
                    });
                }
                if !reduced_action_claims.insert((binding.action.clone(), claim_id.clone())) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "gauge reduction action/claim cell",
                        id: format!("{}:{}", binding.id, claim_id),
                    });
                }
                let measure_required = matches!(
                    &request.claim.information,
                    InformationRegime::PosteriorUnderDeclaredPrior { .. }
                ) || view.carrier.iter().any(|role| {
                    matches!(
                        &problem.document.parameters[role].treatment,
                        ParameterTreatment::Marginalized
                    )
                });
                let reduced = !matches!(&binding.plan, GaugeReductionPlan::Unreduced { .. });
                if reduced
                    && measure_required
                    && !matches!(&binding.measure, GaugeMeasureSemantics::Pushforward { .. })
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "gauge reduction measure semantics",
                        detail: format!(
                            "reduction {} touches a posterior or marginalized claim cell and needs exact pushforward/Jacobian-or-disintegration semantics",
                            binding.id
                        ),
                    });
                }
                if !reduced && matches!(&binding.measure, GaugeMeasureSemantics::Pushforward { .. })
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "unreduced gauge measure semantics",
                        detail: format!(
                            "unreduced binding {} cannot claim a quotient pushforward",
                            binding.id
                        ),
                    });
                }
                if matches!(
                    &binding.plan,
                    GaugeReductionPlan::ContinuousReductionWithDiscreteResidual { .. }
                ) && view.orbit_kind != EffectiveGaugeOrbitKind::Mixed
                {
                    return Err(IdentifiabilityError::InvalidText {
                        field: "continuous reduction with discrete residual",
                        detail:
                            "this reduction is available only for an exact mixed effective orbit"
                                .to_string(),
                    });
                }
            }
            match &binding.plan {
                GaugeReductionPlan::Slice { slice } => {
                    validate_gauge_slice(&binding.action, slice, &problem.document)?;
                }
                GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
                    reduction: ContinuousGaugeReductionPlan::Slice { slice },
                    ..
                } => {
                    validate_gauge_slice(&binding.action, slice, &problem.document)?;
                }
                GaugeReductionPlan::Unreduced { .. }
                | GaugeReductionPlan::Quotient { .. }
                | GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
                    reduction: ContinuousGaugeReductionPlan::Quotient { .. },
                    ..
                } => {}
            }
            for source in gauge_reduction_sources(&binding.plan)? {
                bind_source_reference(
                    &mut execution_sources,
                    source,
                    "execution gauge reduction source alias",
                )?;
            }
            for source in gauge_reduction_stage_sources(&binding.stage)? {
                bind_source_reference(
                    &mut execution_sources,
                    source,
                    "execution gauge stage source alias",
                )?;
            }
            for source in gauge_measure_sources(&binding.measure)? {
                bind_source_reference(
                    &mut execution_sources,
                    source,
                    "execution gauge measure source alias",
                )?;
            }
            gauge_reduction_map.insert(binding.id.clone(), binding);
        }
        validate_gauge_reduction_dag(&gauge_reduction_map)?;
        for (claim_id, request) in &claim_requests {
            if let Some(action) = claim_gauge_action(&request.claim) {
                let coverage = gauge_reduction_map
                    .values()
                    .filter(|binding| {
                        &binding.action == action && binding.claims.contains(claim_id)
                    })
                    .count();
                if coverage != 1 {
                    return Err(IdentifiabilityError::Cardinality {
                        field: "gauge reduction coverage",
                        detail: format!(
                            "claim {claim_id} action target needs exactly one explicit reduced-or-unreduced plan, found {coverage}"
                        ),
                    });
                }
            }
        }
        validate_source_authority_closure(
            &execution_sources,
            &source_authority,
            "execution source authority",
        )?;
        for (key, resolution) in &source_authority.entries {
            if let Some(problem_resolution) = problem.source_admission.resolutions.get(key)
                && problem_resolution != resolution
            {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "execution/problem source authority",
                });
            }
        }
        let plan = Self {
            schema_version: IDENTIFIABILITY_EXECUTION_IDENTITY_VERSION,
            header,
            problem_id: problem.problem_id,
            source_admission_id: problem.source_admission_id,
            analyzer,
            build,
            derivative_provider,
            claim_requests,
            actions: action_map,
            gauge_reductions: gauge_reduction_map,
            numerical,
            initialization,
            stopping,
            determinism_contract,
            source_authority,
        };
        validate_execution_structural_budget(&plan)?;
        let _ = encode_execution(&plan)?;
        Ok(plan)
    }

    pub fn canonical_bytes(&self) -> Result<Vec<u8>, IdentifiabilityError> {
        encode_execution(self)
    }

    pub fn id(&self) -> Result<ExecutionId, IdentifiabilityError> {
        execution_identity_hash(self)
    }

    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn problem_id(&self) -> ProblemId {
        self.problem_id
    }

    #[must_use]
    pub const fn source_admission_id(&self) -> SourceAdmissionId {
        self.source_admission_id
    }

    /// Analyzer identity to which assessment methods must be exactly bound.
    #[must_use]
    pub const fn analyzer(&self) -> &SourceRef {
        &self.analyzer
    }

    #[must_use]
    pub const fn build(&self) -> &SourceRef {
        &self.build
    }

    #[must_use]
    pub const fn derivative_provider(&self) -> Option<&SourceRef> {
        self.derivative_provider.as_ref()
    }

    #[must_use]
    pub const fn claim_requests(&self) -> &BTreeMap<ClaimId, ClaimRequest> {
        &self.claim_requests
    }

    /// Coarse planner projection derived from the exact preregistered
    /// propositions. It is intentionally not independent identity state.
    #[must_use]
    pub fn requested_axes(&self) -> BTreeSet<RequestedClaimAxis> {
        self.claim_requests
            .values()
            .flat_map(|request| required_axes(&request.claim))
            .collect()
    }

    #[must_use]
    pub const fn actions(&self) -> &BTreeMap<ParameterRoleId, ParameterExecutionAction> {
        &self.actions
    }

    #[must_use]
    pub const fn gauge_reductions(&self) -> &BTreeMap<GaugeReductionId, GaugeReductionBinding> {
        &self.gauge_reductions
    }

    /// Analyzer-internal numerical policy. Its rank and singular-value controls
    /// are never compared directly with a claim error bound; claim requests
    /// carry their own exact dimensionless metric and scale.
    #[must_use]
    pub const fn numerical_policy(&self) -> &IdentifiabilityNumericalPolicy {
        &self.numerical
    }

    /// Locally verified authority for every analyzer, build, derivative,
    /// integrator, initialization, stopping, and determinism source.
    #[must_use]
    pub const fn source_authority(&self) -> &SourceResolutionSet {
        &self.source_authority
    }

    #[must_use]
    pub const fn initialization(&self) -> &SourceRef {
        &self.initialization
    }

    #[must_use]
    pub const fn stopping(&self) -> &SourceRef {
        &self.stopping
    }

    #[must_use]
    pub const fn determinism_contract(&self) -> &SourceRef {
        &self.determinism_contract
    }

    pub fn from_canonical_bytes(
        bytes: &[u8],
        problem: &AdmittedIdentifiabilityProblem,
        verified_sources: &SourceResolutionSet,
    ) -> Result<Self, IdentifiabilityError> {
        decode_execution(bytes, problem, verified_sources)
    }
}

/// Information assumed by an identifiability claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InformationRegime {
    StructuralExactModel,
    ExactInputOutputMap,
    NoisyFiniteData,
    PosteriorUnderDeclaredPrior { joint_prior: SourceRef },
}

/// Extent is independent of information regime and quantifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentifiabilityExtent {
    Local,
    Global,
}

/// Scalar/algebraic domain, with every non-real extension bound to exact
/// semantics and independently admitted source bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScalarDomain {
    Real,
    Complex { extension: SourceRef },
    MixedDiscreteContinuous { stratification: SourceRef },
}

/// Geometry and cardinality of the parameter-to-observation fiber, orthogonal
/// to local/global extent. This makes local set-valued, global finite-to-one,
/// orbit-unique, and stratified claims representable without overloading one
/// ordinal axis.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum GaugeActionReference {
    Single(GaugeClassId),
    /// Exact direct product. The referenced declaration must have
    /// [`GaugeCompositionKind::IndependentProduct`]; disjoint support alone is
    /// never treated as a commutation theorem.
    Product(GaugeCompositionId),
    /// Exact generated/interacting action. The referenced declaration must
    /// have [`GaugeCompositionKind::Generated`].
    Composition(GaugeCompositionId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FiberCardinalityBound {
    UniformU64(u64),
    SymbolicProfile(SourceRef),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FiberDimensionLowerBound {
    Finite { minimum_dimension: u64 },
    InfiniteDimensional { model_space: SourceRef },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FiberStructure {
    Unique,
    FiniteToOne {
        maximum_cardinality: Option<FiberCardinalityBound>,
    },
    DiscreteOrbit {
        action: GaugeActionReference,
    },
    /// Orbit with both continuous and discrete effective components. A
    /// computational reduction choice is execution identity and is not needed
    /// merely to state this physical proposition.
    MixedOrbit {
        action: GaugeActionReference,
    },
    OrbitQuotientUnique {
        action: GaugeActionReference,
    },
    GeneralizedQuotientUnique {
        action: GaugeActionReference,
        equivalence: SourceRef,
    },
    PositiveDimensional {
        lower_bound: FiberDimensionLowerBound,
    },
    StratifiedOrbit {
        action: GaugeActionReference,
        orbit_type_profile: SourceRef,
    },
    Stratified {
        strata: SourceRef,
    },
}

fn claim_gauge_action(claim: &TypedIdentifiabilityClaim) -> Option<&GaugeActionReference> {
    match &claim.fiber {
        FiberStructure::DiscreteOrbit { action }
        | FiberStructure::MixedOrbit { action }
        | FiberStructure::OrbitQuotientUnique { action }
        | FiberStructure::GeneralizedQuotientUnique { action, .. }
        | FiberStructure::StratifiedOrbit { action, .. } => Some(action),
        FiberStructure::Unique
        | FiberStructure::FiniteToOne { .. }
        | FiberStructure::PositiveDimensional { .. }
        | FiberStructure::Stratified { .. } => match &claim.subject {
            ClaimSubject::GaugeAction(action) => Some(action),
            _ => None,
        },
    }
}

/// Mathematical quantifier and its exact domain/measure source.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimQuantifier {
    AtRealization {
        realization: SourceRef,
    },
    AlmostEverywhere {
        measure: SourceRef,
    },
    ForAll {
        domain: SourceRef,
    },
    ProbabilityAtLeast {
        probability: f64,
        measure: SourceRef,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ClaimSubject {
    Parameter(ParameterRoleId),
    ParameterSet(BTreeSet<ParameterRoleId>),
    DerivedFunctional {
        definition: SourceRef,
        parameters: BTreeSet<ParameterRoleId>,
    },
    Influence(InfluenceId),
    GaugeAction(GaugeActionReference),
    WholeProblem,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ClaimScope {
    WholeCampaign,
    Cases(BTreeSet<CaseId>),
    Stratum {
        definition: SourceRef,
        cases: BTreeSet<CaseId>,
    },
}

/// Coordinate-free proposition.  Its truth status and receipts live in the
/// paired [`ClaimAssessment`], preserving a product type instead of collapsing
/// “structural/local/generic/global/practical” into one ordinal label.
#[derive(Debug, Clone, PartialEq)]
pub struct TypedIdentifiabilityClaim {
    id: ClaimId,
    information: InformationRegime,
    extent: IdentifiabilityExtent,
    fiber: FiberStructure,
    quantifier: ClaimQuantifier,
    scalar_domain: ScalarDomain,
    subject: ClaimSubject,
    scope: ClaimScope,
}

impl TypedIdentifiabilityClaim {
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub const fn new(
        id: ClaimId,
        information: InformationRegime,
        extent: IdentifiabilityExtent,
        fiber: FiberStructure,
        quantifier: ClaimQuantifier,
        scalar_domain: ScalarDomain,
        subject: ClaimSubject,
        scope: ClaimScope,
    ) -> Self {
        Self {
            id,
            information,
            extent,
            fiber,
            quantifier,
            scalar_domain,
            subject,
            scope,
        }
    }

    #[must_use]
    pub const fn id(&self) -> &ClaimId {
        &self.id
    }

    #[must_use]
    pub const fn information(&self) -> &InformationRegime {
        &self.information
    }

    #[must_use]
    pub const fn extent(&self) -> IdentifiabilityExtent {
        self.extent
    }

    #[must_use]
    pub const fn fiber(&self) -> &FiberStructure {
        &self.fiber
    }

    #[must_use]
    pub const fn quantifier(&self) -> &ClaimQuantifier {
        &self.quantifier
    }

    #[must_use]
    pub const fn scalar_domain(&self) -> &ScalarDomain {
        &self.scalar_domain
    }

    #[must_use]
    pub const fn subject(&self) -> &ClaimSubject {
        &self.subject
    }

    #[must_use]
    pub const fn scope(&self) -> &ClaimScope {
        &self.scope
    }
}

/// Evidence-bound *claim* about a proposition. Positive/refuting variants are
/// deliberately prefixed `Claimed`: content verification and even an external
/// trust receipt do not, by themselves, prove the mathematical proposition.
/// A future method-specific verifier may promote a subject-bound receipt to a
/// sealed theorem token without changing this honest transport layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GaugeResolutionDisposition {
    CandidateRefuted,
    NoProjectionOnSubject,
    SubjectDescendsToQuotient,
    BrokenByJointInformation,
    TrivialResidualIntersection,
    ConsistentWithClaimedFiber,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeResolutionEvidence {
    action: GaugeActionReference,
    disposition: GaugeResolutionDisposition,
    method: SourceRef,
    receipt: SourceRef,
}

impl GaugeResolutionEvidence {
    #[must_use]
    pub const fn new(
        action: GaugeActionReference,
        disposition: GaugeResolutionDisposition,
        method: SourceRef,
        receipt: SourceRef,
    ) -> Self {
        Self {
            action,
            disposition,
            method,
            receipt,
        }
    }

    #[must_use]
    pub const fn action(&self) -> &GaugeActionReference {
        &self.action
    }

    #[must_use]
    pub const fn disposition(&self) -> &GaugeResolutionDisposition {
        &self.disposition
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ClaimAssessment {
    ClaimedEstablished {
        method: SourceRef,
        receipt: SourceRef,
        metric: SourceRef,
        nondimensionalization: SourceRef,
        certified_error_bound: f64,
        /// Canonical, action-keyed evidence. A map is intentional: caller
        /// insertion order must never perturb assessment identity, and the
        /// key/value action agreement is validated before admission.
        gauge_resolutions: BTreeMap<GaugeActionReference, GaugeResolutionEvidence>,
    },
    ClaimedRefuted {
        method: SourceRef,
        receipt: SourceRef,
        metric: SourceRef,
        nondimensionalization: SourceRef,
        certified_error_bound: f64,
    },
    ClaimedInconclusive {
        method: Option<SourceRef>,
        receipt: Option<SourceRef>,
        reason: String,
    },
    NotAssessed {
        reason: String,
    },
}

/// Typed conclusions for one exact execution.
#[derive(Clone, PartialEq)]
pub struct IdentifiabilityAssessment {
    schema_version: u32,
    header: ArtifactHeader,
    problem_id: ProblemId,
    execution_id: ExecutionId,
    claims: BTreeMap<ClaimId, TypedIdentifiabilityClaim>,
    evidence: BTreeMap<ClaimId, ClaimAssessment>,
    source_authority: SourceResolutionSet,
}

impl fmt::Debug for IdentifiabilityAssessment {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IdentifiabilityAssessment")
            .field("stage", &"assessment")
            .field("claim_count", &self.claims.len())
            .field("evidence_count", &self.evidence.len())
            .field("source_count", &self.source_authority.entries.len())
            .finish_non_exhaustive()
    }
}

fn validate_assessment_structural_budget(
    assessment: &IdentifiabilityAssessment,
) -> Result<(), IdentifiabilityError> {
    let mut budget = StructuralItemBudget::default();
    budget.add(assessment.claims.len(), "assessment claims")?;
    for claim in assessment.claims.values() {
        add_claim_structural_items(&mut budget, claim)?;
    }
    budget.add(assessment.evidence.len(), "assessment evidence")?;
    for evidence in assessment.evidence.values() {
        if let ClaimAssessment::ClaimedEstablished {
            gauge_resolutions, ..
        } = evidence
        {
            budget.add(
                gauge_resolutions.len(),
                "assessment gauge resolution evidence",
            )?;
        }
    }
    budget.add(
        assessment.source_authority.entries.len(),
        "assessment source authority",
    )?;
    Ok(())
}

#[allow(
    dead_code,
    unused_variables,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn classify_identifiability_problem_identity_fields(
    document: &IdentifiabilityProblemDocument,
    source_ref: &SourceRef,
    source_kind: &SourceKind,
    parameter_purpose: &ParameterPurpose,
    conditioned_value: &ConditionedValue,
    parameter_treatment: &ParameterTreatment,
    prior_policy: &PriorPolicy,
    influence_coverage: &InfluenceCoverage,
    owner: &ParameterOwnerBinding,
    parameter_scope: &ParameterScopeBinding,
    parameter: &StudyParameter,
    affine_term: &AffineConstraintTerm,
    constraint_relation: &ConstraintRelation,
    constraint_kind: &JointConstraintKind,
    constraint: &JointConstraint,
    opaque_membership: &OpaqueDomainMembershipClaim,
    admissible_domain: &AdmissibleDomainWitness,
    case_purpose: &CasePurpose,
    case_data: &CaseDataDeclaration,
    observation_rows: &ObservationRows,
    marginal_noise: &MarginalNoiseSpec,
    missingness: &MissingnessAssumption,
    observation: &StudyObservation,
    joint_noise: &JointNoiseModel,
    discrepancy_inapplicability: &DiscrepancyInapplicability,
    discrepancy: &StudyDiscrepancy,
    physics_sources: &CasePhysicsSources,
    observation_sharing: &ObservationSharingGroup,
    case: &StudyCaseDocument,
    observation_key: &ObservationKey,
    functional: &DistributionFunctional,
    representation: &InfluenceRepresentation,
    influence: &InfluenceDeclaration,
    constraint_codimension: &ConstraintCodimension,
    gauge_continuous_dimension: &GaugeContinuousDimension,
    gauge_discrete_size: &GaugeDiscreteSize,
    gauge_discrete_orbit_cardinality: &GaugeDiscreteOrbitCardinality,
    regular_gauge_orbit: &RegularGaugeOrbit,
    gauge_algebra: &GaugeAlgebra,
    gauge_orbit: &GaugeOrbitGeometry,
    gauge_status: &GaugeStatus,
    gauge_information: &GaugeInformationRegime,
    gauge_scalar_domain: &GaugeScalarDomain,
    gauge_locus: &GaugeLocus,
    gauge_probability: &GaugeProbabilityThreshold,
    gauge_quantifier: &GaugeQuantifierScope,
    gauge_axes: &GaugeApplicabilityAxes,
    gauge_extent: &GaugeExtentSupport,
    gauge_cell: &GaugeCellDomain,
    gauge_validity: &GaugeValidityScope,
    gauge: &GaugeDeclaration,
    composition_kind: &GaugeCompositionKind,
    composition: &GaugeCompositionDeclaration,
    data_sharing: &DataSharingGroup,
    data_reuse: &DataReusePolicy,
) {
    let IdentifiabilityProblemDocument {
        schema_version,
        context_source,
        material_source,
        model_source,
        graph_source,
        joint_prior,
        sources,
        parameters,
        constraints,
        admissible_domain,
        cases,
        influences,
        gauges,
        gauge_compositions,
        joint_noise,
        data_reuse,
    } = document;
    let SourceRef {
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
    } = source_ref;
    match source_kind {
        SourceKind::ContextOfUse
        | SourceKind::MaterialCard
        | SourceKind::ConstitutiveModelCard
        | SourceKind::ConstitutiveGraph
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit
        | SourceKind::ForwardModel
        | SourceKind::Geometry
        | SourceKind::Process
        | SourceKind::Protocol
        | SourceKind::ObservationOperator
        | SourceKind::Metrology
        | SourceKind::Parser
        | SourceKind::Preprocessing
        | SourceKind::Likelihood
        | SourceKind::Prior
        | SourceKind::Constraint
        | SourceKind::GaugeAction
        | SourceKind::GaugeSection
        | SourceKind::Discrepancy
        | SourceKind::Assumption
        | SourceKind::Analyzer
        | SourceKind::DerivativeProvider
        | SourceKind::Build
        | SourceKind::EvidenceReceipt
        | SourceKind::ExternalManifold
        | SourceKind::AlgebraicExtension
        | SourceKind::Stratification
        | SourceKind::DerivedFunctional
        | SourceKind::AdmissibleDomainCertificate
        | SourceKind::DimensionlessErrorMetric
        | SourceKind::Nondimensionalization
        | SourceKind::QuantifierRealization
        | SourceKind::ReferenceMeasure
        | SourceKind::ProbabilityMeasure
        | SourceKind::QuantifierDomain
        | SourceKind::GaugeComposition
        | SourceKind::GaugeOrbitTypeProfile
        | SourceKind::GaugeHypothesis
        | SourceKind::GaugeGroupPresentation
        | SourceKind::GaugeOrbitPresentation
        | SourceKind::InfluenceComposition
        | SourceKind::GaugeQuotientProfile
        | SourceKind::FiberCardinalityProfile
        | SourceKind::FiberDimensionProfile
        | SourceKind::FunctionalModelSpace
        | SourceKind::GaugeQuotientMap
        | SourceKind::GaugeInvariantMap
        | SourceKind::GaugeGroupoidPresentation
        | SourceKind::GaugeReductionLaw
        | SourceKind::GaugeSubgroupCertificate
        | SourceKind::GaugeResidualAction
        | SourceKind::GaugeMeasureTransport
        | SourceKind::MeasureTransport
        | SourceKind::ParameterizedLikelihood
        | SourceKind::UnitDefinition
        | SourceKind::ForwardModelProductionBinding => {}
    }
    match parameter_purpose {
        ParameterPurpose::Estimand
        | ParameterPurpose::Nuisance
        | ParameterPurpose::Hyperparameter
        | ParameterPurpose::CalibrationControl => {}
    }
    let ConditionedValue { value_si, source } = conditioned_value;
    match parameter_treatment {
        ParameterTreatment::Estimated
        | ParameterTreatment::Profiled
        | ParameterTreatment::Marginalized => {}
        ParameterTreatment::Conditioned(value) => {
            let _ = value;
        }
        ParameterTreatment::Derived {
            definition,
            parents,
        } => {
            let _ = (definition, parents);
        }
    }
    match prior_policy {
        PriorPolicy::Distribution(prior) => {
            let _ = prior;
        }
        PriorPolicy::Absent { reason } | PriorPolicy::NotApplicable { reason } => {
            let _ = reason;
        }
    }
    match influence_coverage {
        InfluenceCoverage::Declared => {}
        InfluenceCoverage::IntentionallyAbsent { reason }
        | InfluenceCoverage::NotApplicable { reason } => {
            let _ = reason;
        }
    }
    match owner {
        ParameterOwnerBinding::ConstitutiveModel => {}
        ParameterOwnerBinding::InitialState { state_path } => {
            let _ = state_path;
        }
        ParameterOwnerBinding::Instrument {
            instrument,
            acquisition_channel,
            metrology,
        } => {
            let _ = (instrument, acquisition_channel, metrology);
        }
        ParameterOwnerBinding::Discrepancy { family } => {
            let _ = family;
        }
        ParameterOwnerBinding::ControlledInput { protocol } => {
            let _ = protocol;
        }
        ParameterOwnerBinding::Population { hierarchy } => {
            let _ = hierarchy;
        }
    }
    match parameter_scope {
        ParameterScopeBinding::Global => {}
        ParameterScopeBinding::Cases(cases) => {
            let _ = cases;
        }
        ParameterScopeBinding::MaterialLot { lot, cases } => {
            let _ = (lot, cases);
        }
        ParameterScopeBinding::Specimen { case, specimen } => {
            let _ = (case, specimen);
        }
        ParameterScopeBinding::Field { support, cases } => {
            let _ = (support, cases);
        }
        ParameterScopeBinding::Hierarchical {
            population,
            level,
            hierarchy,
            cases,
        } => {
            let _ = (population, level, hierarchy, cases);
        }
    }
    let StudyParameter {
        role,
        quantity,
        domain,
        purpose,
        treatment,
        owner,
        scope,
        prior,
        influence_coverage,
    } = parameter;
    let AffineConstraintTerm {
        parameter,
        coefficient,
        coefficient_quantity,
    } = affine_term;
    match constraint_relation {
        ConstraintRelation::Equal
        | ConstraintRelation::LessOrEqual
        | ConstraintRelation::GreaterOrEqual => {}
    }
    match constraint_kind {
        JointConstraintKind::Affine {
            terms,
            relation,
            rhs_si,
            residual_quantity,
        } => {
            let _ = (terms, relation, rhs_si, residual_quantity);
        }
        JointConstraintKind::Simplex {
            members,
            total_si,
            quantity,
        } => {
            let _ = (members, total_si, quantity);
        }
        JointConstraintKind::Ordered { members, strict } => {
            let _ = (members, strict);
        }
        JointConstraintKind::ExternalManifold {
            members,
            definition,
            codimension,
        } => {
            let _ = (members, definition, codimension);
        }
        JointConstraintKind::StochasticCoupling {
            members,
            distribution,
        } => {
            let _ = (members, distribution);
        }
    }
    let JointConstraint { id, kind } = constraint;
    let OpaqueDomainMembershipClaim {
        source,
        witness_binding,
    } = opaque_membership;
    let AdmissibleDomainWitness {
        values,
        opaque_membership_claim,
    } = admissible_domain;
    match case_purpose {
        CasePurpose::Calibration
        | CasePurpose::SymmetryBreaking
        | CasePurpose::ValidationOnly
        | CasePurpose::BlindFalsification
        | CasePurpose::ProspectiveDesign => {}
        CasePurpose::Complementary { reason } => {
            let _ = reason;
        }
    }
    match case_data {
        CaseDataDeclaration::Prospective => {}
        CaseDataDeclaration::Retrospective {
            experiment,
            split,
            parser,
            preprocessing,
            parser_version,
            split_grouping,
        } => {
            let _ = (
                experiment,
                split,
                parser,
                preprocessing,
                parser_version,
                split_grouping,
            );
        }
    }
    match observation_rows {
        ObservationRows::Prospective => {}
        ObservationRows::Retrospective(rows) => {
            let _ = rows;
        }
    }
    match marginal_noise {
        MarginalNoiseSpec::Gaussian { standard_deviation } => {
            let _ = standard_deviation;
        }
        MarginalNoiseSpec::StudentT {
            scale,
            degrees_of_freedom,
        } => {
            let _ = (scale, degrees_of_freedom);
        }
        MarginalNoiseSpec::Empirical {
            distribution,
            standard_deviation,
            finite_variance_model,
        } => {
            let _ = (distribution, standard_deviation, finite_variance_model);
        }
        MarginalNoiseSpec::Bounded { half_width } => {
            let _ = half_width;
        }
        MarginalNoiseSpec::Unknown { reason } => {
            let _ = reason;
        }
    }
    match missingness {
        MissingnessAssumption::Complete { assumption } => {
            let _ = assumption;
        }
        MissingnessAssumption::Modeled { mechanism } => {
            let _ = mechanism;
        }
        MissingnessAssumption::Unknown { reason } => {
            let _ = reason;
        }
    }
    let StudyObservation {
        id,
        qoi,
        unit,
        quantity,
        unit_definition,
        frame,
        graph_node,
        graph_port,
        operator,
        aggregation,
        sensor,
        instrument,
        acquisition_channel,
        clock,
        operator_version,
        noise,
        missingness,
        saturation,
        protocol_version,
        refinement_version,
        rows,
    } = observation;
    match joint_noise {
        JointNoiseModel::Independent { assumption } => {
            let _ = assumption;
        }
        JointNoiseModel::DenseCorrelation {
            order,
            correlation,
            model,
        } => {
            let _ = (order, correlation, model);
        }
        JointNoiseModel::ExternalKernel { model } => {
            let _ = model;
        }
        JointNoiseModel::Unknown { reason } => {
            let _ = reason;
        }
    }
    match discrepancy_inapplicability {
        DiscrepancyInapplicability::PhysicalApplicability { assumption }
        | DiscrepancyInapplicability::ProspectiveDesign { assumption } => {
            let _ = assumption;
        }
        DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
            generator,
            producer,
            production_binding,
            assumption,
        } => {
            let _ = (generator, producer, production_binding, assumption);
        }
    }
    match discrepancy {
        StudyDiscrepancy::Uncharacterized { reason } => {
            let _ = reason;
        }
        StudyDiscrepancy::NotApplicable { basis } => {
            let _ = basis;
        }
        StudyDiscrepancy::AssumedZero { assumption } => {
            let _ = assumption;
        }
        StudyDiscrepancy::Modeled {
            family,
            parameters,
            support,
            confounding_guard,
        } => {
            let _ = (family, parameters, support, confounding_guard);
        }
    }
    let CasePhysicsSources {
        frame_transform,
        specimen_geometry,
        specimen_process,
        specimen_preparation,
        load_path,
        environment_path,
        time_grid,
        initial_state,
    } = physics_sources;
    let ObservationSharingGroup {
        channels,
        rows,
        joint_likelihood,
        justification,
    } = observation_sharing;
    let StudyCaseDocument {
        id,
        purpose,
        initial_state,
        specimen,
        protocol,
        physics_sources,
        forward_model,
        data,
        observations,
        discrepancies,
        observation_sharing,
    } = case;
    let ObservationKey { case, channel } = observation_key;
    match functional {
        DistributionFunctional::Location { observation }
        | DistributionFunctional::LogScale { observation }
        | DistributionFunctional::MissingnessLogit { observation }
        | DistributionFunctional::CensoringLogit { observation } => {
            let _ = observation;
        }
        DistributionFunctional::Correlation { left, right } => {
            let _ = (left, right);
        }
    }
    match representation {
        InfluenceRepresentation::Direct => {}
        InfluenceRepresentation::StateMediated { state_path } => {
            let _ = state_path;
        }
        InfluenceRepresentation::Composite { operator, inputs } => {
            let _ = (operator, inputs);
        }
        InfluenceRepresentation::ExternalDefinition { definition } => {
            let _ = definition;
        }
    }
    let InfluenceDeclaration {
        id,
        parameter,
        functional,
        representation,
    } = influence;
    match constraint_codimension {
        ConstraintCodimension::Finite { codimension } => {
            let _ = codimension;
        }
        ConstraintCodimension::InfiniteDimensional { profile } => {
            let _ = profile;
        }
    }
    match gauge_continuous_dimension {
        GaugeContinuousDimension::Finite { dimension } => {
            let _ = dimension;
        }
        GaugeContinuousDimension::InfiniteDimensional { model_space } => {
            let _ = model_space;
        }
    }
    match gauge_discrete_size {
        GaugeDiscreteSize::Finite { order } => {
            let _ = order;
        }
        GaugeDiscreteSize::CountablyInfinite { presentation } => {
            let _ = presentation;
        }
    }
    match gauge_discrete_orbit_cardinality {
        GaugeDiscreteOrbitCardinality::Finite { cardinality } => {
            let _ = cardinality;
        }
        GaugeDiscreteOrbitCardinality::CountablyInfinite { presentation } => {
            let _ = presentation;
        }
    }
    let RegularGaugeOrbit {
        continuous_orbit_dimension,
        discrete_orbit_cardinality,
    } = regular_gauge_orbit;
    match gauge_algebra {
        GaugeAlgebra::Continuous { group_dimension } => {
            let _ = group_dimension;
        }
        GaugeAlgebra::Discrete { size } => {
            let _ = size;
        }
        GaugeAlgebra::Mixed {
            continuous_group_dimension,
            component_group,
        } => {
            let _ = (continuous_group_dimension, component_group);
        }
    }
    match gauge_orbit {
        GaugeOrbitGeometry::Regular {
            principal,
            stabilizer_profile,
        } => {
            let _ = (principal, stabilizer_profile);
        }
        GaugeOrbitGeometry::Stratified {
            principal,
            orbit_type_stabilizer_profile,
        } => {
            let _ = (principal, orbit_type_stabilizer_profile);
        }
    }
    match gauge_status {
        GaugeStatus::Candidate { rationale } => {
            let _ = rationale;
        }
        GaugeStatus::Assumed { assumption } => {
            let _ = assumption;
        }
    }
    match gauge_information {
        GaugeInformationRegime::StructuralExactModel
        | GaugeInformationRegime::ExactInputOutputMap
        | GaugeInformationRegime::NoisyFiniteData => {}
        GaugeInformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            let _ = joint_prior;
        }
    }
    match gauge_scalar_domain {
        GaugeScalarDomain::Real => {}
        GaugeScalarDomain::Complex { extension } => {
            let _ = extension;
        }
        GaugeScalarDomain::MixedDiscreteContinuous { stratification } => {
            let _ = stratification;
        }
    }
    match gauge_locus {
        GaugeLocus::WholeDomain => {}
        GaugeLocus::Stratum { definition } => {
            let _ = definition;
        }
    }
    let GaugeProbabilityThreshold(probability_bits) = gauge_probability;
    match gauge_quantifier {
        GaugeQuantifierScope::AtRealization { realization } => {
            let _ = realization;
        }
        GaugeQuantifierScope::AlmostEverywhere { measure } => {
            let _ = measure;
        }
        GaugeQuantifierScope::ForAll { domain } => {
            let _ = domain;
        }
        GaugeQuantifierScope::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            let _ = (probability, measure);
        }
    }
    let GaugeApplicabilityAxes {
        information,
        scalar_domain,
        locus,
        quantifier,
    } = gauge_axes;
    let GaugeExtentSupport {
        local_obstruction_parameters,
        global_obstruction_parameters,
    } = gauge_extent;
    let GaugeCellDomain {
        case_obstruction_support,
    } = gauge_cell;
    let GaugeValidityScope { cells } = gauge_validity;
    let GaugeDeclaration {
        id,
        members,
        action,
        algebra,
        orbit_geometry,
        status,
        validity,
    } = gauge;
    match composition_kind {
        GaugeCompositionKind::IndependentProduct | GaugeCompositionKind::Generated => {}
    }
    let GaugeCompositionDeclaration {
        id: composition_id,
        members: composition_members,
        kind: declared_composition_kind,
        law,
        effective_algebra,
        effective_orbit_geometry,
        status: composition_status,
        validity: composition_validity,
    } = composition;
    let DataSharingGroup {
        cases,
        joint_likelihood,
        justification,
    } = data_sharing;
    match data_reuse {
        DataReusePolicy::Disjoint => {}
        DataReusePolicy::Shared { groups } => {
            let _ = groups;
        }
    }
    let _ = (
        schema_version,
        context_source,
        material_source,
        model_source,
        graph_source,
        joint_prior,
        sources,
        parameters,
        constraints,
        admissible_domain,
        cases,
        influences,
        gauges,
        gauge_compositions,
        joint_noise,
        data_reuse,
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
        value_si,
        source,
        role,
        quantity,
        domain,
        purpose,
        treatment,
        owner,
        scope,
        prior,
        influence_coverage,
        parameter,
        coefficient,
        coefficient_quantity,
        id,
        kind,
        source,
        witness_binding,
        values,
        opaque_membership_claim,
        qoi,
        unit,
        unit_definition,
        frame,
        graph_node,
        graph_port,
        operator,
        aggregation,
        sensor,
        instrument,
        acquisition_channel,
        clock,
        operator_version,
        noise,
        missingness,
        saturation,
        protocol_version,
        refinement_version,
        rows,
        frame_transform,
        specimen_geometry,
        specimen_process,
        specimen_preparation,
        load_path,
        environment_path,
        time_grid,
        initial_state,
        channels,
        joint_likelihood,
        justification,
        specimen,
        protocol,
        physics_sources,
        forward_model,
        data,
        observations,
        discrepancies,
        observation_sharing,
        case,
        channel,
        functional,
        representation,
        members,
        action,
        algebra,
        orbit_geometry,
        status,
        validity,
        composition_id,
        composition_members,
        declared_composition_kind,
        law,
        effective_algebra,
        effective_orbit_geometry,
        composition_status,
        composition_validity,
        continuous_orbit_dimension,
        discrete_orbit_cardinality,
        probability_bits,
        cases,
    );
}

/// Tuple-payload enums cannot yet be registered as `sources` by the identity
/// checker. This compiler-exhaustive classifier keeps their variants and
/// payload positions inside the problem schema-function fingerprint until the
/// checker gains variant-qualified tuple-field support.
#[allow(dead_code)]
fn classify_identifiability_problem_tuple_schema(
    treatment: &ParameterTreatment,
    prior: &PriorPolicy,
    parameter_prior: &ParameterPrior,
    scope: &ParameterScopeBinding,
    rows: &ObservationRows,
) {
    match treatment {
        ParameterTreatment::Estimated
        | ParameterTreatment::Profiled
        | ParameterTreatment::Marginalized => {}
        ParameterTreatment::Conditioned(value) => {
            let _ = value;
        }
        ParameterTreatment::Derived {
            definition,
            parents,
        } => {
            let _ = (definition, parents);
        }
    }
    match prior {
        PriorPolicy::Distribution(distribution) => {
            let _ = distribution;
        }
        PriorPolicy::Absent { reason } | PriorPolicy::NotApplicable { reason } => {
            let _ = reason;
        }
    }
    match parameter_prior {
        ParameterPrior::None { reason, version } => {
            let _ = (reason, version);
        }
        ParameterPrior::Uniform { domain, version } => {
            let _ = (domain, version);
        }
        ParameterPrior::Gaussian {
            mean,
            standard_deviation,
            version,
        } => {
            let _ = (mean, standard_deviation, version);
        }
        ParameterPrior::LogNormal {
            log_mean,
            log_standard_deviation,
            reference,
            version,
        } => {
            let _ = (log_mean, log_standard_deviation, reference, version);
        }
    }
    match scope {
        ParameterScopeBinding::Global => {}
        ParameterScopeBinding::Cases(cases) => {
            let _ = cases;
        }
        ParameterScopeBinding::MaterialLot { lot, cases } => {
            let _ = (lot, cases);
        }
        ParameterScopeBinding::Specimen { case, specimen } => {
            let _ = (case, specimen);
        }
        ParameterScopeBinding::Field { support, cases } => {
            let _ = (support, cases);
        }
        ParameterScopeBinding::Hierarchical {
            population,
            level,
            hierarchy,
            cases,
        } => {
            let _ = (population, level, hierarchy, cases);
        }
    }
    match rows {
        ObservationRows::Prospective => {}
        ObservationRows::Retrospective(rows) => {
            let _ = rows;
        }
    }
}

#[allow(
    dead_code,
    unused_variables,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn classify_identifiability_source_admission_identity_fields(
    admission: &SourceAdmissionRecord,
    resolution: &SourceResolution,
    source_kind: &SourceKind,
    authority: &AuthorityDisposition,
    trust_receipt: &TrustReceiptRef,
    trust_authentication: &TrustAuthentication,
    verification: &SourceVerification,
) {
    let SourceAdmissionRecord {
        schema_version,
        problem_id,
        resolutions,
    } = admission;
    let SourceResolution {
        key,
        kind,
        resolved_hash,
        content_hash_domain,
        contract_version,
        authority,
        verification,
    } = resolution;
    match source_kind {
        SourceKind::ContextOfUse
        | SourceKind::MaterialCard
        | SourceKind::ConstitutiveModelCard
        | SourceKind::ConstitutiveGraph
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit
        | SourceKind::ForwardModel
        | SourceKind::Geometry
        | SourceKind::Process
        | SourceKind::Protocol
        | SourceKind::ObservationOperator
        | SourceKind::Metrology
        | SourceKind::Parser
        | SourceKind::Preprocessing
        | SourceKind::Likelihood
        | SourceKind::Prior
        | SourceKind::Constraint
        | SourceKind::GaugeAction
        | SourceKind::GaugeSection
        | SourceKind::Discrepancy
        | SourceKind::Assumption
        | SourceKind::Analyzer
        | SourceKind::DerivativeProvider
        | SourceKind::Build
        | SourceKind::EvidenceReceipt
        | SourceKind::ExternalManifold
        | SourceKind::AlgebraicExtension
        | SourceKind::Stratification
        | SourceKind::DerivedFunctional
        | SourceKind::AdmissibleDomainCertificate
        | SourceKind::DimensionlessErrorMetric
        | SourceKind::Nondimensionalization
        | SourceKind::QuantifierRealization
        | SourceKind::ReferenceMeasure
        | SourceKind::ProbabilityMeasure
        | SourceKind::QuantifierDomain
        | SourceKind::GaugeComposition
        | SourceKind::GaugeOrbitTypeProfile
        | SourceKind::GaugeHypothesis
        | SourceKind::GaugeGroupPresentation
        | SourceKind::GaugeOrbitPresentation
        | SourceKind::InfluenceComposition
        | SourceKind::GaugeQuotientProfile
        | SourceKind::FiberCardinalityProfile
        | SourceKind::FiberDimensionProfile
        | SourceKind::FunctionalModelSpace
        | SourceKind::GaugeQuotientMap
        | SourceKind::GaugeInvariantMap
        | SourceKind::GaugeGroupoidPresentation
        | SourceKind::GaugeReductionLaw
        | SourceKind::GaugeSubgroupCertificate
        | SourceKind::GaugeResidualAction
        | SourceKind::GaugeMeasureTransport
        | SourceKind::MeasureTransport
        | SourceKind::ParameterizedLikelihood
        | SourceKind::UnitDefinition
        | SourceKind::ForwardModelProductionBinding => {}
    }
    match authority {
        AuthorityDisposition::ContentVerified => {}
        AuthorityDisposition::ExternalTrustReceipt { trust_receipt } => {
            let _ = trust_receipt;
        }
        AuthorityDisposition::Unverified { reason } => {
            let _ = reason;
        }
    }
    let TrustReceiptRef {
        receipt,
        subject,
        subject_artifact,
        authentication,
    } = trust_receipt;
    match trust_authentication {
        TrustAuthentication::Unauthenticated => {}
        TrustAuthentication::IssuerPolicy {
            issuer,
            trust_policy,
        } => {
            let _ = (issuer, trust_policy);
        }
    }
    match verification {
        SourceVerification::TypedArtifact | SourceVerification::Unverified => {}
        SourceVerification::HashPreimage { byte_len } => {
            let _ = byte_len;
        }
    }
    let _ = (
        schema_version,
        problem_id,
        resolutions,
        key,
        kind,
        resolved_hash,
        content_hash_domain,
        contract_version,
        authority,
        receipt,
        subject,
        subject_artifact,
        authentication,
        verification,
    );
}

#[allow(
    dead_code,
    unused_variables,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn classify_identifiability_execution_identity_fields(
    plan: &IdentifiabilityExecutionPlan,
    source_ref: &SourceRef,
    source_kind: &SourceKind,
    resolution_set: &SourceResolutionSet,
    resolution: &SourceResolution,
    authority: &AuthorityDisposition,
    trust_receipt: &TrustReceiptRef,
    trust_authentication: &TrustAuthentication,
    verification: &SourceVerification,
    action: &ParameterExecutionAction,
    arithmetic: &ArithmeticPolicy,
    error_policy: &DimensionlessErrorPolicy,
    claim_request: &ClaimRequest,
    numerical_policy: &IdentifiabilityNumericalPolicy,
    information: &InformationRegime,
    extent: &IdentifiabilityExtent,
    scalar_domain: &ScalarDomain,
    gauge_action_reference: &GaugeActionReference,
    fiber_cardinality_bound: &FiberCardinalityBound,
    fiber_dimension_lower_bound: &FiberDimensionLowerBound,
    fiber: &FiberStructure,
    quantifier: &ClaimQuantifier,
    claim_subject: &ClaimSubject,
    claim_scope: &ClaimScope,
    claim: &TypedIdentifiabilityClaim,
    gauge_slice_codimension: &GaugeSliceCodimension,
    gauge_slice: &GaugeSlicePlan,
    gauge_quotient: &GaugeQuotientPlan,
    continuous_reduction: &ContinuousGaugeReductionPlan,
    gauge_reduction_plan: &GaugeReductionPlan,
    gauge_stage_relation: &GaugeReductionStageRelation,
    gauge_stage: &GaugeReductionStage,
    gauge_measure: &GaugeMeasureSemantics,
    gauge_binding: &GaugeReductionBinding,
) {
    let IdentifiabilityExecutionPlan {
        schema_version,
        header,
        problem_id,
        source_admission_id,
        analyzer,
        build,
        derivative_provider,
        claim_requests,
        actions,
        gauge_reductions,
        numerical,
        initialization,
        stopping,
        determinism_contract,
        source_authority,
    } = plan;
    let SourceRef {
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
    } = source_ref;
    match source_kind {
        SourceKind::ContextOfUse
        | SourceKind::MaterialCard
        | SourceKind::ConstitutiveModelCard
        | SourceKind::ConstitutiveGraph
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit
        | SourceKind::ForwardModel
        | SourceKind::Geometry
        | SourceKind::Process
        | SourceKind::Protocol
        | SourceKind::ObservationOperator
        | SourceKind::Metrology
        | SourceKind::Parser
        | SourceKind::Preprocessing
        | SourceKind::Likelihood
        | SourceKind::Prior
        | SourceKind::Constraint
        | SourceKind::GaugeAction
        | SourceKind::GaugeSection
        | SourceKind::Discrepancy
        | SourceKind::Assumption
        | SourceKind::Analyzer
        | SourceKind::DerivativeProvider
        | SourceKind::Build
        | SourceKind::EvidenceReceipt
        | SourceKind::ExternalManifold
        | SourceKind::AlgebraicExtension
        | SourceKind::Stratification
        | SourceKind::DerivedFunctional
        | SourceKind::AdmissibleDomainCertificate
        | SourceKind::DimensionlessErrorMetric
        | SourceKind::Nondimensionalization
        | SourceKind::QuantifierRealization
        | SourceKind::ReferenceMeasure
        | SourceKind::ProbabilityMeasure
        | SourceKind::QuantifierDomain
        | SourceKind::GaugeComposition
        | SourceKind::GaugeOrbitTypeProfile
        | SourceKind::GaugeHypothesis
        | SourceKind::GaugeGroupPresentation
        | SourceKind::GaugeOrbitPresentation
        | SourceKind::InfluenceComposition
        | SourceKind::GaugeQuotientProfile
        | SourceKind::FiberCardinalityProfile
        | SourceKind::FiberDimensionProfile
        | SourceKind::FunctionalModelSpace
        | SourceKind::GaugeQuotientMap
        | SourceKind::GaugeInvariantMap
        | SourceKind::GaugeGroupoidPresentation
        | SourceKind::GaugeReductionLaw
        | SourceKind::GaugeSubgroupCertificate
        | SourceKind::GaugeResidualAction
        | SourceKind::GaugeMeasureTransport
        | SourceKind::MeasureTransport
        | SourceKind::ParameterizedLikelihood
        | SourceKind::UnitDefinition
        | SourceKind::ForwardModelProductionBinding => {}
    }
    let SourceResolutionSet { entries } = resolution_set;
    let SourceResolution {
        key: resolution_key,
        kind: resolution_kind,
        resolved_hash,
        content_hash_domain: resolution_domain,
        contract_version: resolution_version,
        authority: resolution_authority,
        verification: resolution_verification,
    } = resolution;
    match authority {
        AuthorityDisposition::ContentVerified => {}
        AuthorityDisposition::ExternalTrustReceipt { trust_receipt } => {
            let _ = trust_receipt;
        }
        AuthorityDisposition::Unverified { reason } => {
            let _ = reason;
        }
    }
    let TrustReceiptRef {
        receipt: trust_receipt_source,
        subject: trust_receipt_subject,
        subject_artifact: trust_receipt_subject_artifact,
        authentication: receipt_authentication,
    } = trust_receipt;
    match trust_authentication {
        TrustAuthentication::Unauthenticated => {}
        TrustAuthentication::IssuerPolicy {
            issuer,
            trust_policy,
        } => {
            let _ = (issuer, trust_policy);
        }
    }
    match verification {
        SourceVerification::TypedArtifact | SourceVerification::Unverified => {}
        SourceVerification::HashPreimage { byte_len } => {
            let _ = byte_len;
        }
    }
    match action {
        ParameterExecutionAction::Optimize { coordinate }
        | ParameterExecutionAction::Profile { coordinate } => {
            let _ = coordinate;
        }
        ParameterExecutionAction::Marginalize {
            coordinate,
            integrator,
            measure_transport,
        } => {
            let _ = (coordinate, integrator, measure_transport);
        }
        ParameterExecutionAction::Conditioned | ParameterExecutionAction::Derived => {}
    }
    match arithmetic {
        ArithmeticPolicy::ExactSymbolic
        | ArithmeticPolicy::CertifiedInterval
        | ArithmeticPolicy::DeterministicFloatingPoint
        | ArithmeticPolicy::FastFloatingPoint => {}
    }
    let DimensionlessErrorPolicy {
        metric,
        nondimensionalization,
        maximum_certified_error,
    } = error_policy;
    let ClaimRequest {
        claim: requested_claim,
        error_policy: requested_error_policy,
    } = claim_request;
    let IdentifiabilityNumericalPolicy {
        rank_tolerance,
        singular_value_floor,
        maximum_condition_number,
        arithmetic,
        nondimensionalization: numerical_nondimensionalization,
    } = numerical_policy;
    match information {
        InformationRegime::StructuralExactModel
        | InformationRegime::ExactInputOutputMap
        | InformationRegime::NoisyFiniteData => {}
        InformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            let _ = joint_prior;
        }
    }
    match extent {
        IdentifiabilityExtent::Local | IdentifiabilityExtent::Global => {}
    }
    match scalar_domain {
        ScalarDomain::Real => {}
        ScalarDomain::Complex { extension } => {
            let _ = extension;
        }
        ScalarDomain::MixedDiscreteContinuous { stratification } => {
            let _ = stratification;
        }
    }
    match gauge_action_reference {
        GaugeActionReference::Single(gauge) => {
            let _ = gauge;
        }
        GaugeActionReference::Product(composition)
        | GaugeActionReference::Composition(composition) => {
            let _ = composition;
        }
    }
    match fiber_cardinality_bound {
        FiberCardinalityBound::UniformU64(maximum) => {
            let _ = maximum;
        }
        FiberCardinalityBound::SymbolicProfile(profile) => {
            let _ = profile;
        }
    }
    match fiber_dimension_lower_bound {
        FiberDimensionLowerBound::Finite { minimum_dimension } => {
            let _ = minimum_dimension;
        }
        FiberDimensionLowerBound::InfiniteDimensional { model_space } => {
            let _ = model_space;
        }
    }
    match fiber {
        FiberStructure::Unique => {}
        FiberStructure::FiniteToOne {
            maximum_cardinality,
        } => {
            let _ = maximum_cardinality;
        }
        FiberStructure::DiscreteOrbit { action }
        | FiberStructure::MixedOrbit { action }
        | FiberStructure::OrbitQuotientUnique { action } => {
            let _ = action;
        }
        FiberStructure::GeneralizedQuotientUnique {
            action,
            equivalence,
        } => {
            let _ = (action, equivalence);
        }
        FiberStructure::PositiveDimensional { lower_bound } => {
            let _ = lower_bound;
        }
        FiberStructure::StratifiedOrbit {
            action,
            orbit_type_profile,
        } => {
            let _ = (action, orbit_type_profile);
        }
        FiberStructure::Stratified { strata } => {
            let _ = strata;
        }
    }
    match quantifier {
        ClaimQuantifier::AtRealization { realization } => {
            let _ = realization;
        }
        ClaimQuantifier::AlmostEverywhere { measure } => {
            let _ = measure;
        }
        ClaimQuantifier::ForAll { domain } => {
            let _ = domain;
        }
        ClaimQuantifier::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            let _ = (probability, measure);
        }
    }
    classify_identifiability_claim_tuple_schema(claim_subject, claim_scope);
    match gauge_slice_codimension {
        GaugeSliceCodimension::FixedFinite { codimension } => {
            let _ = codimension;
        }
        GaugeSliceCodimension::FixedInfinite {
            codimension_model,
            compatibility,
        } => {
            let _ = (codimension_model, compatibility);
        }
        GaugeSliceCodimension::Stratified { profile } => {
            let _ = profile;
        }
    }
    let GaugeSlicePlan {
        support: slice_support,
        constraint: slice_constraint,
        expected_codimension,
        coverage: slice_coverage,
    } = gauge_slice;
    match gauge_quotient {
        GaugeQuotientPlan::RegularAtlas {
            quotient_map,
            local_section_atlas,
            coverage,
        } => {
            let _ = (quotient_map, local_section_atlas, coverage);
        }
        GaugeQuotientPlan::SingularOrGeneralized {
            quotient_map,
            quotient_profile,
            local_models,
        } => {
            let _ = (quotient_map, quotient_profile, local_models);
        }
        GaugeQuotientPlan::InvariantMap {
            invariants,
            completeness_profile,
        } => {
            let _ = (invariants, completeness_profile);
        }
        GaugeQuotientPlan::GroupoidOrStack {
            presentation,
            quotient_profile,
        } => {
            let _ = (presentation, quotient_profile);
        }
    }
    match continuous_reduction {
        ContinuousGaugeReductionPlan::Quotient { quotient } => {
            let _ = quotient;
        }
        ContinuousGaugeReductionPlan::Slice { slice } => {
            let _ = slice;
        }
    }
    match gauge_reduction_plan {
        GaugeReductionPlan::Unreduced { reason } => {
            let _ = reason;
        }
        GaugeReductionPlan::Quotient { quotient } => {
            let _ = quotient;
        }
        GaugeReductionPlan::Slice { slice } => {
            let _ = slice;
        }
        GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction,
            normal_subgroup,
            factor_extension,
            residual_quotient_action,
            compatibility,
        } => {
            let _ = (
                reduction,
                normal_subgroup,
                factor_extension,
                residual_quotient_action,
                compatibility,
            );
        }
    }
    match gauge_stage_relation {
        GaugeReductionStageRelation::NormalSubgroupTower {
            normality,
            induced_residual_action,
        } => {
            let _ = (normality, induced_residual_action);
        }
        GaugeReductionStageRelation::SemidirectOrGenerated {
            extension,
            induced_residual_action,
        } => {
            let _ = (extension, induced_residual_action);
        }
        GaugeReductionStageRelation::TransverseSlices { transversality } => {
            let _ = transversality;
        }
        GaugeReductionStageRelation::GaugeForGauge {
            reducibility,
            induced_residual_action,
        } => {
            let _ = (reducibility, induced_residual_action);
        }
    }
    match gauge_stage {
        GaugeReductionStage::Root => {}
        GaugeReductionStage::After {
            predecessors,
            composition_law,
            relation,
        } => {
            let _ = (predecessors, composition_law, relation);
        }
    }
    match gauge_measure {
        GaugeMeasureSemantics::NotApplicable { reason } => {
            let _ = reason;
        }
        GaugeMeasureSemantics::Pushforward {
            source_measure,
            reduced_measure,
            transport,
            jacobian_or_disintegration,
        } => {
            let _ = (
                source_measure,
                reduced_measure,
                transport,
                jacobian_or_disintegration,
            );
        }
    }
    let GaugeReductionBinding {
        id: reduction_id,
        action: reduction_action,
        claims: reduction_claims,
        plan: reduction_plan,
        stage: reduction_stage,
        measure: reduction_measure,
    } = gauge_binding;
    let TypedIdentifiabilityClaim {
        id,
        information,
        extent,
        fiber,
        quantifier,
        scalar_domain,
        subject,
        scope,
    } = claim;
    let _ = (
        schema_version,
        header,
        problem_id,
        source_admission_id,
        analyzer,
        build,
        derivative_provider,
        claim_requests,
        actions,
        gauge_reductions,
        numerical,
        initialization,
        stopping,
        determinism_contract,
        source_authority,
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
        entries,
        resolution_key,
        resolution_kind,
        resolved_hash,
        resolution_domain,
        resolution_version,
        resolution_authority,
        resolution_verification,
        trust_receipt_source,
        trust_receipt_subject,
        trust_receipt_subject_artifact,
        receipt_authentication,
        metric,
        nondimensionalization,
        maximum_certified_error,
        requested_claim,
        requested_error_policy,
        rank_tolerance,
        singular_value_floor,
        maximum_condition_number,
        arithmetic,
        numerical_nondimensionalization,
        id,
        information,
        extent,
        fiber,
        quantifier,
        scalar_domain,
        subject,
        scope,
        slice_support,
        slice_constraint,
        expected_codimension,
        slice_coverage,
        reduction_id,
        reduction_action,
        reduction_claims,
        reduction_plan,
        reduction_stage,
        reduction_measure,
    );
}

/// The checker cannot yet register tuple-payload claim enums as structural
/// sources. Both execution and assessment declarations include this
/// compiler-exhaustive function in their schema fingerprints instead.
#[allow(dead_code)]
fn classify_identifiability_claim_tuple_schema(subject: &ClaimSubject, scope: &ClaimScope) {
    match subject {
        ClaimSubject::Parameter(parameter) => {
            let _ = parameter;
        }
        ClaimSubject::ParameterSet(parameters) => {
            let _ = parameters;
        }
        ClaimSubject::DerivedFunctional {
            definition,
            parameters,
        } => {
            let _ = (definition, parameters);
        }
        ClaimSubject::Influence(influence) => {
            let _ = influence;
        }
        ClaimSubject::GaugeAction(action) => {
            let _ = action;
        }
        ClaimSubject::WholeProblem => {}
    }
    match scope {
        ClaimScope::WholeCampaign => {}
        ClaimScope::Cases(cases) => {
            let _ = cases;
        }
        ClaimScope::Stratum { definition, cases } => {
            let _ = (definition, cases);
        }
    }
}

#[allow(
    dead_code,
    unused_variables,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn classify_identifiability_assessment_identity_fields(
    assessment: &IdentifiabilityAssessment,
    source_ref: &SourceRef,
    source_kind: &SourceKind,
    resolution_set: &SourceResolutionSet,
    resolution: &SourceResolution,
    authority: &AuthorityDisposition,
    trust_receipt: &TrustReceiptRef,
    trust_authentication: &TrustAuthentication,
    verification: &SourceVerification,
    information: &InformationRegime,
    extent: &IdentifiabilityExtent,
    scalar_domain: &ScalarDomain,
    gauge_action_reference: &GaugeActionReference,
    fiber_cardinality_bound: &FiberCardinalityBound,
    fiber_dimension_lower_bound: &FiberDimensionLowerBound,
    fiber: &FiberStructure,
    quantifier: &ClaimQuantifier,
    claim_subject: &ClaimSubject,
    claim_scope: &ClaimScope,
    claim: &TypedIdentifiabilityClaim,
    gauge_resolution_disposition: &GaugeResolutionDisposition,
    gauge_resolution_evidence: &GaugeResolutionEvidence,
    claim_assessment: &ClaimAssessment,
) {
    let IdentifiabilityAssessment {
        schema_version,
        header,
        problem_id,
        execution_id,
        claims,
        evidence,
        source_authority,
    } = assessment;
    let SourceRef {
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
    } = source_ref;
    match source_kind {
        SourceKind::ContextOfUse
        | SourceKind::MaterialCard
        | SourceKind::ConstitutiveModelCard
        | SourceKind::ConstitutiveGraph
        | SourceKind::ExperimentArtifact
        | SourceKind::CalibrationSplit
        | SourceKind::ForwardModel
        | SourceKind::Geometry
        | SourceKind::Process
        | SourceKind::Protocol
        | SourceKind::ObservationOperator
        | SourceKind::Metrology
        | SourceKind::Parser
        | SourceKind::Preprocessing
        | SourceKind::Likelihood
        | SourceKind::Prior
        | SourceKind::Constraint
        | SourceKind::GaugeAction
        | SourceKind::GaugeSection
        | SourceKind::Discrepancy
        | SourceKind::Assumption
        | SourceKind::Analyzer
        | SourceKind::DerivativeProvider
        | SourceKind::Build
        | SourceKind::EvidenceReceipt
        | SourceKind::ExternalManifold
        | SourceKind::AlgebraicExtension
        | SourceKind::Stratification
        | SourceKind::DerivedFunctional
        | SourceKind::AdmissibleDomainCertificate
        | SourceKind::DimensionlessErrorMetric
        | SourceKind::Nondimensionalization
        | SourceKind::QuantifierRealization
        | SourceKind::ReferenceMeasure
        | SourceKind::ProbabilityMeasure
        | SourceKind::QuantifierDomain
        | SourceKind::GaugeComposition
        | SourceKind::GaugeOrbitTypeProfile
        | SourceKind::GaugeHypothesis
        | SourceKind::GaugeGroupPresentation
        | SourceKind::GaugeOrbitPresentation
        | SourceKind::InfluenceComposition
        | SourceKind::GaugeQuotientProfile
        | SourceKind::FiberCardinalityProfile
        | SourceKind::FiberDimensionProfile
        | SourceKind::FunctionalModelSpace
        | SourceKind::GaugeQuotientMap
        | SourceKind::GaugeInvariantMap
        | SourceKind::GaugeGroupoidPresentation
        | SourceKind::GaugeReductionLaw
        | SourceKind::GaugeSubgroupCertificate
        | SourceKind::GaugeResidualAction
        | SourceKind::GaugeMeasureTransport
        | SourceKind::MeasureTransport
        | SourceKind::ParameterizedLikelihood
        | SourceKind::UnitDefinition
        | SourceKind::ForwardModelProductionBinding => {}
    }
    let SourceResolutionSet { entries } = resolution_set;
    let SourceResolution {
        key: resolution_key,
        kind: resolution_kind,
        resolved_hash,
        content_hash_domain: resolution_domain,
        contract_version: resolution_version,
        authority: resolution_authority,
        verification: resolution_verification,
    } = resolution;
    match authority {
        AuthorityDisposition::ContentVerified => {}
        AuthorityDisposition::ExternalTrustReceipt { trust_receipt } => {
            let _ = trust_receipt;
        }
        AuthorityDisposition::Unverified { reason } => {
            let _ = reason;
        }
    }
    let TrustReceiptRef {
        receipt: trust_receipt_source,
        subject: trust_receipt_subject,
        subject_artifact: trust_receipt_subject_artifact,
        authentication: receipt_authentication,
    } = trust_receipt;
    match trust_authentication {
        TrustAuthentication::Unauthenticated => {}
        TrustAuthentication::IssuerPolicy {
            issuer,
            trust_policy,
        } => {
            let _ = (issuer, trust_policy);
        }
    }
    match verification {
        SourceVerification::TypedArtifact | SourceVerification::Unverified => {}
        SourceVerification::HashPreimage { byte_len } => {
            let _ = byte_len;
        }
    }
    match information {
        InformationRegime::StructuralExactModel
        | InformationRegime::ExactInputOutputMap
        | InformationRegime::NoisyFiniteData => {}
        InformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            let _ = joint_prior;
        }
    }
    match extent {
        IdentifiabilityExtent::Local | IdentifiabilityExtent::Global => {}
    }
    match scalar_domain {
        ScalarDomain::Real => {}
        ScalarDomain::Complex { extension } => {
            let _ = extension;
        }
        ScalarDomain::MixedDiscreteContinuous { stratification } => {
            let _ = stratification;
        }
    }
    match gauge_action_reference {
        GaugeActionReference::Single(gauge) => {
            let _ = gauge;
        }
        GaugeActionReference::Product(composition)
        | GaugeActionReference::Composition(composition) => {
            let _ = composition;
        }
    }
    match fiber_cardinality_bound {
        FiberCardinalityBound::UniformU64(maximum) => {
            let _ = maximum;
        }
        FiberCardinalityBound::SymbolicProfile(profile) => {
            let _ = profile;
        }
    }
    match fiber_dimension_lower_bound {
        FiberDimensionLowerBound::Finite { minimum_dimension } => {
            let _ = minimum_dimension;
        }
        FiberDimensionLowerBound::InfiniteDimensional { model_space } => {
            let _ = model_space;
        }
    }
    match fiber {
        FiberStructure::Unique => {}
        FiberStructure::FiniteToOne {
            maximum_cardinality,
        } => {
            let _ = maximum_cardinality;
        }
        FiberStructure::DiscreteOrbit { action }
        | FiberStructure::MixedOrbit { action }
        | FiberStructure::OrbitQuotientUnique { action } => {
            let _ = action;
        }
        FiberStructure::GeneralizedQuotientUnique {
            action,
            equivalence,
        } => {
            let _ = (action, equivalence);
        }
        FiberStructure::PositiveDimensional { lower_bound } => {
            let _ = lower_bound;
        }
        FiberStructure::StratifiedOrbit {
            action,
            orbit_type_profile,
        } => {
            let _ = (action, orbit_type_profile);
        }
        FiberStructure::Stratified { strata } => {
            let _ = strata;
        }
    }
    match quantifier {
        ClaimQuantifier::AtRealization { realization } => {
            let _ = realization;
        }
        ClaimQuantifier::AlmostEverywhere { measure } => {
            let _ = measure;
        }
        ClaimQuantifier::ForAll { domain } => {
            let _ = domain;
        }
        ClaimQuantifier::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            let _ = (probability, measure);
        }
    }
    classify_identifiability_claim_tuple_schema(claim_subject, claim_scope);
    let TypedIdentifiabilityClaim {
        id,
        information,
        extent,
        fiber,
        quantifier,
        scalar_domain,
        subject,
        scope,
    } = claim;
    match gauge_resolution_disposition {
        GaugeResolutionDisposition::CandidateRefuted
        | GaugeResolutionDisposition::NoProjectionOnSubject
        | GaugeResolutionDisposition::SubjectDescendsToQuotient
        | GaugeResolutionDisposition::BrokenByJointInformation
        | GaugeResolutionDisposition::TrivialResidualIntersection
        | GaugeResolutionDisposition::ConsistentWithClaimedFiber => {}
    }
    let GaugeResolutionEvidence {
        action: resolution_action,
        disposition: resolution_disposition,
        method: resolution_method,
        receipt: resolution_receipt,
    } = gauge_resolution_evidence;
    match claim_assessment {
        ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound,
            gauge_resolutions,
        } => {
            let _ = (
                method,
                receipt,
                metric,
                nondimensionalization,
                certified_error_bound,
                gauge_resolutions,
            );
        }
        ClaimAssessment::ClaimedRefuted {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound,
        } => {
            let _ = (
                method,
                receipt,
                metric,
                nondimensionalization,
                certified_error_bound,
            );
        }
        ClaimAssessment::ClaimedInconclusive {
            method,
            receipt,
            reason,
        } => {
            let _ = (method, receipt, reason);
        }
        ClaimAssessment::NotAssessed { reason } => {
            let _ = reason;
        }
    }
    let _ = (
        schema_version,
        header,
        problem_id,
        execution_id,
        claims,
        evidence,
        source_authority,
        key,
        kind,
        expected_hash,
        content_hash_domain,
        contract_version,
        entries,
        resolution_key,
        resolution_kind,
        resolved_hash,
        resolution_domain,
        resolution_version,
        resolution_authority,
        resolution_verification,
        trust_receipt_source,
        trust_receipt_subject,
        trust_receipt_subject_artifact,
        receipt_authentication,
        id,
        information,
        extent,
        fiber,
        quantifier,
        scalar_domain,
        subject,
        scope,
        resolution_action,
        resolution_disposition,
        resolution_method,
        resolution_receipt,
    );
}

/// Owner-local declaration for the unresolved physical-question identity.
///
/// The declaration is intentionally literal rather than macro-generated:
/// `xtask check-identities` fingerprints this exact source surface and requires
/// every owner field to remain deliberately classified.
#[allow(dead_code)]
pub const IDENTIFIABILITY_PROBLEM_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-material:identifiability-problem",
    "version_const=IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-material.identifiability-problem.v3",
    "domain_const=IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN",
    "encoder=problem_identity_hash",
    "encoder_helpers=IdentifiabilityProblemDocument::canonical_bytes,encode_problem",
    "schema_constants=IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION,IDENTIFIABILITY_PROBLEM_IDENTITY_DOMAIN,PROBLEM_MAGIC,CASE_PHYSICS_SOURCE_CONTRACT_VERSION,FRAME_TRANSFORM_SOURCE_DOMAIN,SPECIMEN_GEOMETRY_SOURCE_DOMAIN,SPECIMEN_PROCESS_SOURCE_DOMAIN,SPECIMEN_PREPARATION_SOURCE_DOMAIN,LOAD_PATH_SOURCE_DOMAIN,ENVIRONMENT_PATH_SOURCE_DOMAIN,TIME_GRID_SOURCE_DOMAIN,INITIAL_STATE_SOURCE_DOMAIN,ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_VERSION,ADMISSIBLE_DOMAIN_MEMBERSHIP_SOURCE_DOMAIN,ADMISSIBLE_DOMAIN_WITNESS_BINDING_DOMAIN,FORWARD_MODEL_PRODUCTION_BINDING_VERSION,FORWARD_MODEL_PRODUCTION_BINDING_DOMAIN,MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ID_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_TEXT_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_CANONICAL_BYTES,crates/fs-evidence/src/vv/model.rs#MAX_VV_MATRIX_DIMENSION,crates/fs-qty/src/semantic.rs#QUANTITY_SPEC_ENCODING_VERSION,crates/fs-qty/src/semantic.rs#QUANTITY_SPEC_ENCODED_LEN",
    "schema_functions=IdentifiabilityProblemDocument::try_new,classify_identifiability_problem_identity_fields,validate_problem_structural_budget,insert_unique,require_source,require_source_kind,require_source_kind_in,validate_source_key,validate_derived_parameter_dag,validate_joint_constraint,normalize_joint_noise,initial_state_schema_version,observation_for,parameter_applicable_cases,parameter_active_in_cases,retrospective_experiment,forward_model_production_binding_preimage,validate_gauge_applicability_sources,validate_gauge_algebra_orbit_sources,validate_independent_product_invariants,gauge_algebra_orbit_compatible,regular_orbit_support_compatible,principal_gauge_orbit,gauge_algebra_source_keys,gauge_orbit_source_keys,infinite_dimensional_profile_is_explicit,continuous_orbit_dimension_compatible,GaugeDeclaration::try_new,GaugeCompositionDeclaration::try_new,GaugeValidityScope::try_new,GaugeCellDomain::try_new,GaugeExtentSupport::try_new,GaugeProbabilityThreshold::try_new,decode_problem,check_problem_identity_version,classify_identifiability_problem_tuple_schema,problem_source_reachability,require_case_physics_source,parameter_membership_source_keys,admissible_domain_witness_binding,admissible_domain_membership_certificate_preimage,validate_admissible_domain_witness,declared_parameter_cases,validate_declared_parameter_cases,functional_observations,transitive_influence_ids,sharing_group_membership,encode_source_key,decode_source_key,encode_case_id,decode_case_id,encode_role,decode_role,encode_channel,decode_channel,encode_source_kind,decode_source_kind,encode_source_ref,decode_source_ref,encode_observation_key,decode_observation_key,encode_parameter_treatment,decode_parameter_treatment,encode_prior_policy,decode_prior_policy,encode_owner,decode_owner,encode_scope,decode_scope,encode_study_parameter,decode_study_parameter,encode_constraint,decode_constraint,encode_admissible_domain_witness,decode_admissible_domain_witness,encode_marginal_noise,decode_marginal_noise,encode_missingness,decode_missingness,encode_study_observation,decode_study_observation,encode_discrepancy,decode_discrepancy,encode_case_physics_sources,decode_case_physics_sources,encode_observation_sharing_group,decode_observation_sharing_group,encode_case,decode_case,encode_functional,decode_functional,encode_influence,decode_influence,encode_gauge_information_regime,decode_gauge_information_regime,encode_gauge_continuous_dimension,decode_gauge_continuous_dimension,encode_gauge_discrete_size,decode_gauge_discrete_size,encode_gauge_algebra,decode_gauge_algebra,encode_gauge_discrete_orbit,decode_gauge_discrete_orbit,encode_regular_gauge_orbit,decode_regular_gauge_orbit,encode_gauge_orbit_geometry,decode_gauge_orbit_geometry,encode_gauge_status,decode_gauge_status,encode_gauge_axes,decode_gauge_axes,encode_gauge_validity_scope,decode_gauge_validity_scope,encode_gauge,decode_gauge,encode_gauge_composition,decode_gauge_composition,encode_joint_noise,decode_joint_noise,encode_data_reuse,decode_data_reuse,crates/fs-material/src/identifiability.rs#same_f64,crates/fs-material/src/identifiability.rs#matrix_get,crates/fs-material/src/identifiability.rs#checked_add_dims,crates/fs-material/src/identifiability.rs#ParameterDomain::try_new,crates/fs-material/src/identifiability.rs#ParameterDomain::is_degenerate,crates/fs-material/src/identifiability.rs#ParameterPrior::validate_against,crates/fs-material/src/identifiability.rs#encode_parameter_domain,crates/fs-material/src/identifiability.rs#decode_parameter_domain,crates/fs-material/src/identifiability.rs#encode_prior,crates/fs-material/src/identifiability.rs#decode_prior,crates/fs-material/src/identifiability.rs#encode_initial_state,crates/fs-material/src/identifiability.rs#decode_initial_state,crates/fs-material/src/identifiability.rs#encode_frame,crates/fs-material/src/identifiability.rs#decode_frame,crates/fs-material/src/identifiability.rs#encode_specimen,crates/fs-material/src/identifiability.rs#decode_specimen,crates/fs-material/src/identifiability.rs#encode_protocol,crates/fs-material/src/identifiability.rs#decode_protocol,crates/fs-material/src/identifiability.rs#encode_artifact_id,crates/fs-material/src/identifiability.rs#decode_artifact_id,crates/fs-material/src/identifiability.rs#encode_qoi_id,crates/fs-material/src/identifiability.rs#decode_qoi_id,crates/fs-material/src/identifiability.rs#encode_observation_row_id,crates/fs-material/src/identifiability.rs#decode_observation_row_id,crates/fs-material/src/identifiability.rs#FrameBinding::try_new,crates/fs-material/src/identifiability.rs#SpecimenBinding::try_new,crates/fs-material/src/identifiability.rs#ProtocolBinding::try_new,crates/fs-material/src/identifiability.rs#canonical_f64,crates/fs-material/src/identifiability.rs#validate_token,crates/fs-material/src/identifiability.rs#validate_reason,crates/fs-material/src/identifiability.rs#CanonicalWriter::new,crates/fs-material/src/identifiability.rs#CanonicalWriter::raw,crates/fs-material/src/identifiability.rs#CanonicalWriter::byte,crates/fs-material/src/identifiability.rs#CanonicalWriter::u32,crates/fs-material/src/identifiability.rs#CanonicalWriter::u64,crates/fs-material/src/identifiability.rs#CanonicalWriter::f64,crates/fs-material/src/identifiability.rs#CanonicalWriter::count,crates/fs-material/src/identifiability.rs#CanonicalWriter::text,crates/fs-material/src/identifiability.rs#CanonicalWriter::hash,crates/fs-material/src/identifiability.rs#CanonicalWriter::quantity,crates/fs-material/src/identifiability.rs#CanonicalWriter::finish,crates/fs-material/src/identifiability.rs#CanonicalReader::new,crates/fs-material/src/identifiability.rs#CanonicalReader::take,crates/fs-material/src/identifiability.rs#CanonicalReader::byte,crates/fs-material/src/identifiability.rs#CanonicalReader::u32,crates/fs-material/src/identifiability.rs#CanonicalReader::u64,crates/fs-material/src/identifiability.rs#CanonicalReader::f64,crates/fs-material/src/identifiability.rs#CanonicalReader::length,crates/fs-material/src/identifiability.rs#CanonicalReader::count,crates/fs-material/src/identifiability.rs#CanonicalReader::text,crates/fs-material/src/identifiability.rs#CanonicalReader::token,crates/fs-material/src/identifiability.rs#CanonicalReader::reason,crates/fs-material/src/identifiability.rs#CanonicalReader::hash,crates/fs-material/src/identifiability.rs#CanonicalReader::quantity,crates/fs-material/src/identifiability.rs#CanonicalReader::expect_byte,crates/fs-material/src/identifiability.rs#CanonicalReader::finish,crates/fs-evidence/src/vv/model.rs#CovarianceMatrix::try_new,crates/fs-evidence/src/vv/model.rs#CovarianceMatrix::get,crates/fs-evidence/src/vv/model.rs#CovarianceMatrix::is_positive_semidefinite,crates/fs-evidence/src/vv/model.rs#CovarianceMatrix::dimension,crates/fs-evidence/src/vv/model.rs#CovarianceMatrix::lower_triangle,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-qty/src/semantic.rs#QuantitySpec::canonical_bytes,crates/fs-qty/src/semantic.rs#QuantitySpec::from_canonical_bytes",
    "schema_dependencies=none",
    "digest=blake3-256-domain-separated",
    "encoding=canonical-transport-exact-bits",
    "sources=IdentifiabilityProblemDocument,SourceRef,SourceKind,ParameterPurpose,ConditionedValue,InfluenceCoverage,ParameterOwnerBinding,StudyParameter,AffineConstraintTerm,ConstraintRelation,JointConstraintKind,JointConstraint,OpaqueDomainMembershipClaim,AdmissibleDomainWitness,CasePurpose,CaseDataDeclaration,MarginalNoiseSpec,MissingnessAssumption,StudyObservation,JointNoiseModel,DiscrepancyInapplicability,StudyDiscrepancy,CasePhysicsSources,ObservationSharingGroup,StudyCaseDocument,ObservationKey,DistributionFunctional,InfluenceRepresentation,InfluenceDeclaration,ConstraintCodimension,GaugeContinuousDimension,GaugeDiscreteSize,GaugeDiscreteOrbitCardinality,RegularGaugeOrbit,GaugeAlgebra,GaugeOrbitGeometry,GaugeStatus,GaugeInformationRegime,GaugeScalarDomain,GaugeLocus,GaugeQuantifierScope,GaugeApplicabilityAxes,GaugeExtentSupport,GaugeCellDomain,GaugeValidityScope,GaugeDeclaration,GaugeCompositionKind,GaugeCompositionDeclaration,DataSharingGroup,DataReusePolicy",
    "source_fields=IdentifiabilityProblemDocument.schema_version:semantic,IdentifiabilityProblemDocument.context_source:semantic,IdentifiabilityProblemDocument.material_source:semantic,IdentifiabilityProblemDocument.model_source:semantic,IdentifiabilityProblemDocument.graph_source:semantic,IdentifiabilityProblemDocument.joint_prior:semantic,IdentifiabilityProblemDocument.sources:semantic,IdentifiabilityProblemDocument.parameters:semantic,IdentifiabilityProblemDocument.constraints:semantic,IdentifiabilityProblemDocument.admissible_domain:semantic,IdentifiabilityProblemDocument.cases:semantic,IdentifiabilityProblemDocument.influences:semantic,IdentifiabilityProblemDocument.gauges:semantic,IdentifiabilityProblemDocument.gauge_compositions:semantic,IdentifiabilityProblemDocument.joint_noise:semantic,IdentifiabilityProblemDocument.data_reuse:semantic,SourceRef.key:derived:encoded-by-parent,SourceRef.kind:derived:encoded-by-parent,SourceRef.expected_hash:derived:encoded-by-parent,SourceRef.content_hash_domain:derived:encoded-by-parent,SourceRef.contract_version:derived:encoded-by-parent,SourceKind.variant:derived:encoded-by-parent,ParameterPurpose.variant:derived:encoded-by-parent,ConditionedValue.value_si:derived:encoded-by-parent,ConditionedValue.source:derived:encoded-by-parent,InfluenceCoverage.variant:derived:encoded-by-parent,InfluenceCoverage.reason:derived:encoded-by-parent,ParameterOwnerBinding.variant:derived:encoded-by-parent,ParameterOwnerBinding.state_path:derived:encoded-by-parent,ParameterOwnerBinding.instrument:derived:encoded-by-parent,ParameterOwnerBinding.acquisition_channel:derived:encoded-by-parent,ParameterOwnerBinding.metrology:derived:encoded-by-parent,ParameterOwnerBinding.family:derived:encoded-by-parent,ParameterOwnerBinding.protocol:derived:encoded-by-parent,ParameterOwnerBinding.hierarchy:derived:encoded-by-parent,StudyParameter.role:derived:encoded-by-parent,StudyParameter.quantity:derived:encoded-by-parent,StudyParameter.domain:derived:encoded-by-parent,StudyParameter.purpose:derived:encoded-by-parent,StudyParameter.treatment:derived:encoded-by-parent,StudyParameter.owner:derived:encoded-by-parent,StudyParameter.scope:derived:encoded-by-parent,StudyParameter.prior:derived:encoded-by-parent,StudyParameter.influence_coverage:derived:encoded-by-parent,AffineConstraintTerm.parameter:derived:encoded-by-parent,AffineConstraintTerm.coefficient:derived:encoded-by-parent,AffineConstraintTerm.coefficient_quantity:derived:encoded-by-parent,ConstraintRelation.variant:derived:encoded-by-parent,JointConstraintKind.variant:derived:encoded-by-parent,JointConstraintKind.terms:derived:encoded-by-parent,JointConstraintKind.relation:derived:encoded-by-parent,JointConstraintKind.rhs_si:derived:encoded-by-parent,JointConstraintKind.residual_quantity:derived:encoded-by-parent,JointConstraintKind.members:derived:encoded-by-parent,JointConstraintKind.total_si:derived:encoded-by-parent,JointConstraintKind.quantity:derived:encoded-by-parent,JointConstraintKind.strict:derived:encoded-by-parent,JointConstraintKind.definition:derived:encoded-by-parent,JointConstraintKind.codimension:derived:encoded-by-parent,JointConstraintKind.distribution:derived:encoded-by-parent,JointConstraint.id:derived:encoded-by-parent,JointConstraint.kind:derived:encoded-by-parent,OpaqueDomainMembershipClaim.source:derived:encoded-by-parent,OpaqueDomainMembershipClaim.witness_binding:derived:encoded-by-parent,AdmissibleDomainWitness.values:derived:encoded-by-parent,AdmissibleDomainWitness.opaque_membership_claim:derived:encoded-by-parent,CasePurpose.variant:derived:encoded-by-parent,CasePurpose.reason:derived:encoded-by-parent,CaseDataDeclaration.variant:derived:encoded-by-parent,CaseDataDeclaration.experiment:derived:encoded-by-parent,CaseDataDeclaration.split:derived:encoded-by-parent,CaseDataDeclaration.parser:derived:encoded-by-parent,CaseDataDeclaration.preprocessing:derived:encoded-by-parent,CaseDataDeclaration.parser_version:derived:encoded-by-parent,CaseDataDeclaration.split_grouping:derived:encoded-by-parent,MarginalNoiseSpec.variant:derived:encoded-by-parent,MarginalNoiseSpec.standard_deviation:derived:encoded-by-parent,MarginalNoiseSpec.scale:derived:encoded-by-parent,MarginalNoiseSpec.degrees_of_freedom:derived:encoded-by-parent,MarginalNoiseSpec.distribution:derived:encoded-by-parent,MarginalNoiseSpec.finite_variance_model:derived:encoded-by-parent,MarginalNoiseSpec.half_width:derived:encoded-by-parent,MarginalNoiseSpec.reason:derived:encoded-by-parent,MissingnessAssumption.variant:derived:encoded-by-parent,MissingnessAssumption.assumption:derived:encoded-by-parent,MissingnessAssumption.mechanism:derived:encoded-by-parent,MissingnessAssumption.reason:derived:encoded-by-parent,StudyObservation.id:derived:encoded-by-parent,StudyObservation.qoi:derived:encoded-by-parent,StudyObservation.unit:derived:encoded-by-parent,StudyObservation.quantity:derived:encoded-by-parent,StudyObservation.unit_definition:derived:encoded-by-parent,StudyObservation.frame:derived:encoded-by-parent,StudyObservation.graph_node:derived:encoded-by-parent,StudyObservation.graph_port:derived:encoded-by-parent,StudyObservation.operator:derived:encoded-by-parent,StudyObservation.aggregation:derived:encoded-by-parent,StudyObservation.sensor:derived:encoded-by-parent,StudyObservation.instrument:derived:encoded-by-parent,StudyObservation.acquisition_channel:derived:encoded-by-parent,StudyObservation.clock:derived:encoded-by-parent,StudyObservation.operator_version:derived:encoded-by-parent,StudyObservation.noise:derived:encoded-by-parent,StudyObservation.missingness:derived:encoded-by-parent,StudyObservation.saturation:derived:encoded-by-parent,StudyObservation.protocol_version:derived:encoded-by-parent,StudyObservation.refinement_version:derived:encoded-by-parent,StudyObservation.rows:derived:encoded-by-parent,JointNoiseModel.variant:derived:encoded-by-parent,JointNoiseModel.assumption:derived:encoded-by-parent,JointNoiseModel.order:derived:encoded-by-parent,JointNoiseModel.correlation:derived:encoded-by-parent,JointNoiseModel.model:derived:encoded-by-parent,JointNoiseModel.reason:derived:encoded-by-parent,DiscrepancyInapplicability.variant:derived:encoded-by-parent,DiscrepancyInapplicability.assumption:derived:encoded-by-parent,DiscrepancyInapplicability.generator:derived:encoded-by-parent,DiscrepancyInapplicability.producer:derived:encoded-by-parent,DiscrepancyInapplicability.production_binding:derived:encoded-by-parent,StudyDiscrepancy.variant:derived:encoded-by-parent,StudyDiscrepancy.reason:derived:encoded-by-parent,StudyDiscrepancy.basis:derived:encoded-by-parent,StudyDiscrepancy.assumption:derived:encoded-by-parent,StudyDiscrepancy.family:derived:encoded-by-parent,StudyDiscrepancy.parameters:derived:encoded-by-parent,StudyDiscrepancy.support:derived:encoded-by-parent,StudyDiscrepancy.confounding_guard:derived:encoded-by-parent,CasePhysicsSources.frame_transform:derived:encoded-by-parent,CasePhysicsSources.specimen_geometry:derived:encoded-by-parent,CasePhysicsSources.specimen_process:derived:encoded-by-parent,CasePhysicsSources.specimen_preparation:derived:encoded-by-parent,CasePhysicsSources.load_path:derived:encoded-by-parent,CasePhysicsSources.environment_path:derived:encoded-by-parent,CasePhysicsSources.time_grid:derived:encoded-by-parent,CasePhysicsSources.initial_state:derived:encoded-by-parent,ObservationSharingGroup.channels:derived:encoded-by-parent,ObservationSharingGroup.rows:derived:encoded-by-parent,ObservationSharingGroup.joint_likelihood:derived:encoded-by-parent,ObservationSharingGroup.justification:derived:encoded-by-parent,StudyCaseDocument.id:derived:encoded-by-parent,StudyCaseDocument.purpose:derived:encoded-by-parent,StudyCaseDocument.initial_state:derived:encoded-by-parent,StudyCaseDocument.specimen:derived:encoded-by-parent,StudyCaseDocument.protocol:derived:encoded-by-parent,StudyCaseDocument.physics_sources:derived:encoded-by-parent,StudyCaseDocument.forward_model:derived:encoded-by-parent,StudyCaseDocument.data:derived:encoded-by-parent,StudyCaseDocument.observations:derived:encoded-by-parent,StudyCaseDocument.discrepancies:derived:encoded-by-parent,StudyCaseDocument.observation_sharing:derived:encoded-by-parent,ObservationKey.case:derived:encoded-by-parent,ObservationKey.channel:derived:encoded-by-parent,DistributionFunctional.variant:derived:encoded-by-parent,DistributionFunctional.observation:derived:encoded-by-parent,DistributionFunctional.left:derived:encoded-by-parent,DistributionFunctional.right:derived:encoded-by-parent,InfluenceRepresentation.variant:derived:encoded-by-parent,InfluenceRepresentation.state_path:derived:encoded-by-parent,InfluenceRepresentation.operator:derived:encoded-by-parent,InfluenceRepresentation.inputs:derived:encoded-by-parent,InfluenceRepresentation.definition:derived:encoded-by-parent,InfluenceDeclaration.id:derived:encoded-by-parent,InfluenceDeclaration.parameter:derived:encoded-by-parent,InfluenceDeclaration.functional:derived:encoded-by-parent,InfluenceDeclaration.representation:derived:encoded-by-parent,ConstraintCodimension.variant:derived:encoded-by-parent,ConstraintCodimension.codimension:derived:encoded-by-parent,ConstraintCodimension.profile:derived:encoded-by-parent,GaugeContinuousDimension.variant:derived:encoded-by-parent,GaugeContinuousDimension.dimension:derived:encoded-by-parent,GaugeContinuousDimension.model_space:derived:encoded-by-parent,GaugeDiscreteSize.variant:derived:encoded-by-parent,GaugeDiscreteSize.order:derived:encoded-by-parent,GaugeDiscreteSize.presentation:derived:encoded-by-parent,GaugeDiscreteOrbitCardinality.variant:derived:encoded-by-parent,GaugeDiscreteOrbitCardinality.cardinality:derived:encoded-by-parent,GaugeDiscreteOrbitCardinality.presentation:derived:encoded-by-parent,RegularGaugeOrbit.continuous_orbit_dimension:derived:encoded-by-parent,RegularGaugeOrbit.discrete_orbit_cardinality:derived:encoded-by-parent,GaugeAlgebra.variant:derived:encoded-by-parent,GaugeAlgebra.group_dimension:derived:encoded-by-parent,GaugeAlgebra.size:derived:encoded-by-parent,GaugeAlgebra.continuous_group_dimension:derived:encoded-by-parent,GaugeAlgebra.component_group:derived:encoded-by-parent,GaugeOrbitGeometry.variant:derived:encoded-by-parent,GaugeOrbitGeometry.principal:derived:encoded-by-parent,GaugeOrbitGeometry.stabilizer_profile:derived:encoded-by-parent,GaugeOrbitGeometry.orbit_type_stabilizer_profile:derived:encoded-by-parent,GaugeStatus.variant:derived:encoded-by-parent,GaugeStatus.rationale:derived:encoded-by-parent,GaugeStatus.assumption:derived:encoded-by-parent,GaugeInformationRegime.variant:derived:encoded-by-parent,GaugeInformationRegime.joint_prior:derived:encoded-by-parent,GaugeScalarDomain.variant:derived:encoded-by-parent,GaugeScalarDomain.extension:derived:encoded-by-parent,GaugeScalarDomain.stratification:derived:encoded-by-parent,GaugeLocus.variant:derived:encoded-by-parent,GaugeLocus.definition:derived:encoded-by-parent,GaugeQuantifierScope.variant:derived:encoded-by-parent,GaugeQuantifierScope.realization:derived:encoded-by-parent,GaugeQuantifierScope.measure:derived:encoded-by-parent,GaugeQuantifierScope.domain:derived:encoded-by-parent,GaugeQuantifierScope.probability:derived:encoded-by-parent,GaugeApplicabilityAxes.information:derived:encoded-by-parent,GaugeApplicabilityAxes.scalar_domain:derived:encoded-by-parent,GaugeApplicabilityAxes.locus:derived:encoded-by-parent,GaugeApplicabilityAxes.quantifier:derived:encoded-by-parent,GaugeExtentSupport.local_obstruction_parameters:derived:encoded-by-parent,GaugeExtentSupport.global_obstruction_parameters:derived:encoded-by-parent,GaugeCellDomain.case_obstruction_support:derived:encoded-by-parent,GaugeValidityScope.cells:derived:encoded-by-parent,GaugeDeclaration.id:derived:encoded-by-parent,GaugeDeclaration.members:derived:encoded-by-parent,GaugeDeclaration.action:derived:encoded-by-parent,GaugeDeclaration.algebra:derived:encoded-by-parent,GaugeDeclaration.orbit_geometry:derived:encoded-by-parent,GaugeDeclaration.status:derived:encoded-by-parent,GaugeDeclaration.validity:derived:encoded-by-parent,GaugeCompositionKind.variant:derived:encoded-by-parent,GaugeCompositionDeclaration.id:derived:encoded-by-parent,GaugeCompositionDeclaration.members:derived:encoded-by-parent,GaugeCompositionDeclaration.kind:derived:encoded-by-parent,GaugeCompositionDeclaration.law:derived:encoded-by-parent,GaugeCompositionDeclaration.effective_algebra:derived:encoded-by-parent,GaugeCompositionDeclaration.effective_orbit_geometry:derived:encoded-by-parent,GaugeCompositionDeclaration.status:derived:encoded-by-parent,GaugeCompositionDeclaration.validity:derived:encoded-by-parent,DataSharingGroup.cases:derived:encoded-by-parent,DataSharingGroup.joint_likelihood:derived:encoded-by-parent,DataSharingGroup.justification:derived:encoded-by-parent,DataReusePolicy.variant:derived:encoded-by-parent,DataReusePolicy.groups:derived:encoded-by-parent",
    "source_bindings=IdentifiabilityProblemDocument.schema_version>wire-schema-version,IdentifiabilityProblemDocument.context_source>context-source-binding,IdentifiabilityProblemDocument.material_source>material-source-binding,IdentifiabilityProblemDocument.model_source>model-source-binding,IdentifiabilityProblemDocument.graph_source>graph-source-binding,IdentifiabilityProblemDocument.joint_prior>joint-prior-source,IdentifiabilityProblemDocument.sources>source-registry,IdentifiabilityProblemDocument.parameters>parameter-registry,IdentifiabilityProblemDocument.constraints>joint-constraint-registry,IdentifiabilityProblemDocument.admissible_domain>admissible-domain-witness,IdentifiabilityProblemDocument.cases>study-case-registry,IdentifiabilityProblemDocument.influences>influence-registry,IdentifiabilityProblemDocument.gauges>gauge-registry,IdentifiabilityProblemDocument.gauge_compositions>gauge-composition-registry,IdentifiabilityProblemDocument.joint_noise>joint-noise-model,IdentifiabilityProblemDocument.data_reuse>data-reuse-policy",
    "external_semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian",
    "semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian,wire-schema-version,context-source-binding,material-source-binding,model-source-binding,graph-source-binding,joint-prior-source,source-registry,parameter-registry,joint-constraint-registry,admissible-domain-witness,study-case-registry,influence-registry,gauge-registry,gauge-composition-registry,joint-noise-model,data-reuse-policy",
    "excluded_fields=source-authority-envelope:admission-authority-not-physical-question,execution-configuration:belongs-to-execution-identity,assessment-claims-and-evidence:belongs-to-assessment-identity,caller-container-order:canonicalized-before-identity",
    "consumers=AdmittedIdentifiabilityProblem::resolve_and_admit,AdmittedIdentifiabilityProblem::id,IdentifiabilityExecutionPlan::try_new",
    "mutations=identity-domain:crates/fs-material/tests/identifiability_authority.rs#identity_domains_and_wire_magics_are_stage_separated,identity-version:crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact,wire-magic:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,canonical-field-order:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,collection-count-u32-le:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,fixed-numeric-little-endian:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,wire-schema-version:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,context-source-binding:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,material-source-binding:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,model-source-binding:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,graph-source-binding:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,joint-prior-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,source-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,parameter-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,joint-constraint-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,admissible-domain-witness:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,study-case-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,influence-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,gauge-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,gauge-composition-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,joint-noise-model:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence,data-reuse-policy:crates/fs-material/tests/identifiability_authority.rs#identifiability_problem_identity_bindings_have_exact_mutation_evidence",
    "nonsemantic_mutations=source-authority-envelope:crates/fs-material/tests/identifiability_authority.rs#problem_and_source_admission_identities_separate_question_from_trust_envelope,execution-configuration:crates/fs-material/tests/identifiability_authority.rs#coordinates_do_not_move_problem_identity,assessment-claims-and-evidence:crates/fs-material/tests/identifiability_authority.rs#evidence_changes_assessment_not_problem_or_execution,caller-container-order:crates/fs-material/tests/identifiability_authority.rs#case_and_registry_input_order_are_nonsemantic",
    "field_guard=classify_identifiability_problem_identity_fields",
    "transport_guard=IdentifiabilityProblemDocument::from_canonical_bytes",
    "version_guard=crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact",
    "coupling_surface=fs-material:identifiability-problem",
];

/// Owner-local declaration for the exact source-resolution authority envelope.
#[allow(dead_code)]
pub const IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-material:identifiability-source-admission",
    "version_const=IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-material.identifiability-source-admission.v3",
    "domain_const=IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN",
    "encoder=source_admission_identity_hash",
    "encoder_helpers=AdmittedIdentifiabilityProblem::source_admission_canonical_bytes,encode_source_admission,encode_resolution_entry",
    "schema_constants=IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION,IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_DOMAIN,SOURCE_ADMISSION_MAGIC,BLIND_RELEASE_TRUST_RECEIPT_VERSION,BLIND_RELEASE_TRUST_RECEIPT_DOMAIN,VV_ARTIFACT_SOURCE_DOMAIN,MATERIAL_CARD_SOURCE_DOMAIN,CONSTITUTIVE_MODEL_CARD_SOURCE_DOMAIN,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ID_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_TEXT_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_CANONICAL_BYTES,crates/fs-evidence/src/vv/model.rs#VV_SCHEMA_VERSION,crates/fs-evidence/src/vv/model.rs#VV_RULESET_VERSION,crates/fs-evidence/src/vv/model.rs#VV_ARTIFACT_FAMILY,crates/fs-evidence/src/vv/model.rs#MAX_VV_ID_BYTES,crates/fs-evidence/src/vv/model.rs#MAX_VV_TEXT_BYTES,crates/fs-evidence/src/vv/model.rs#MAX_VV_ITEMS,crates/fs-evidence/src/vv/model.rs#MAX_VV_MATRIX_DIMENSION,crates/fs-evidence/src/vv/codec.rs#MAGIC,crates/fs-evidence/src/vv/codec.rs#CANONICAL_RULE,crates/fs-evidence/src/vv/codec.rs#ROOT_ARTIFACT,crates/fs-evidence/src/vv/codec.rs#MAX_VV_CANONICAL_BYTES,crates/fs-evidence/src/vv/codec.rs#MAX_VV_STRING_BYTES,crates/fs-evidence/src/vv/codec.rs#MAX_VV_COLLECTION_ITEMS,crates/fs-evidence/src/vv/codec.rs#MAX_VV_TOTAL_COLLECTION_ITEMS,crates/fs-matdb/src/cards.rs#MATDB_SCHEMA_VERSION,crates/fs-matdb/src/cards.rs#MODEL_HASH_DOMAIN,crates/fs-matdb/src/cards.rs#MATERIAL_HASH_DOMAIN,crates/fs-matdb/src/cards.rs#CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION,crates/fs-matdb/src/cards.rs#CANONICAL_PARAMETER_BLOCK_IDENTITY_DOMAIN",
    "schema_functions=AdmittedIdentifiabilityProblem::resolve_and_admit,SourceResolution::verify,SourceResolution::unresolved,SourceResolutionSet::try_new,ProblemSourceBundle::with_concrete_authority,TrustReceiptRef::try_new,TrustReceiptRef::try_new_with_subject_artifact,TrustReceiptRef::blind_release,concrete_resolution,validate_authority_disposition,validate_authority_subject,validate_authority_subject_with_artifact,validate_authority_subject_fields,concrete_authority_for,insert_exact_resolution,admit_opaque_resolution,bind_source_reference,validate_source_authority_closure,sharing_group_membership,check_source_admission_identity_version,encode_source_key,encode_source_kind,encode_trust_receipt_ref,decode_trust_receipt_ref,encode_resolution_verification,crates/fs-material/src/identifiability.rs#hash_is_nonzero,crates/fs-material/src/identifiability.rs#canonical_f64,crates/fs-material/src/identifiability.rs#validate_token,crates/fs-material/src/identifiability.rs#validate_reason,crates/fs-material/src/identifiability.rs#ContextBinding::from_vv,crates/fs-material/src/identifiability.rs#ContextBinding::validate_structural,crates/fs-material/src/identifiability.rs#MaterialModelBinding::from_cards,crates/fs-material/src/identifiability.rs#MaterialModelBinding::validate_structural,crates/fs-material/src/identifiability.rs#ModelParameterBinding::nominal,crates/fs-material/src/identifiability.rs#InitialStateBinding::validate_against,crates/fs-material/src/identifiability.rs#DataLineage::from_vv,crates/fs-material/src/identifiability.rs#DataLineage::validate_structural,crates/fs-material/src/identifiability.rs#DataLineage::qois,crates/fs-material/src/identifiability.rs#DataLineage::row_bindings,crates/fs-material/src/identifiability.rs#DataLineage::source_bytes,crates/fs-material/src/identifiability.rs#DataLineage::raw_manifest,crates/fs-material/src/identifiability.rs#CanonicalWriter::new,crates/fs-material/src/identifiability.rs#CanonicalWriter::raw,crates/fs-material/src/identifiability.rs#CanonicalWriter::byte,crates/fs-material/src/identifiability.rs#CanonicalWriter::u32,crates/fs-material/src/identifiability.rs#CanonicalWriter::u64,crates/fs-material/src/identifiability.rs#CanonicalWriter::count,crates/fs-material/src/identifiability.rs#CanonicalWriter::text,crates/fs-material/src/identifiability.rs#CanonicalWriter::hash,crates/fs-material/src/identifiability.rs#CanonicalWriter::finish,crates/fs-matdb/src/cards.rs#ConstitutiveModelCard::validate,crates/fs-matdb/src/cards.rs#ConstitutiveModelCard::content_hash,crates/fs-matdb/src/cards.rs#ConstitutiveModelCard::canonical_parameters_hash,crates/fs-matdb/src/cards.rs#ConstitutiveModelCard::canonical_parameters_hash_with_schema,crates/fs-matdb/src/cards.rs#MaterialCard::content_hash,crates/fs-matdb/src/cards.rs#MaterialCard::models,crates/fs-matdb/src/cards.rs#MaterialCard::schema_version,crates/fs-matdb/src/lib.rs#dims_bytes,crates/fs-matdb/src/lib.rs#Provenance::validate,crates/fs-evidence/src/vv/model.rs#ContextOfUse::try_new,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::try_new,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::validate,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::try_new,crates/fs-evidence/src/vv/model.rs#ObservationManifest::try_new,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::dataset_source_bytes_hash,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::locator_domain,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::locator_contract_version,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::locator_hash,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::extraction_receipt_hash,crates/fs-evidence/src/vv/model.rs#ObservationSourceRef::locator_identity,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::source_ref,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::locator_hash,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::qoi,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::instrument,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::acquisition_channel,crates/fs-evidence/src/vv/model.rs#ObservationManifestRow::clock,crates/fs-evidence/src/vv/model.rs#ObservationManifest::row,crates/fs-evidence/src/vv/model.rs#ObservationManifest::rows,crates/fs-evidence/src/vv/model.rs#ObservationManifest::locator_hash_of,crates/fs-evidence/src/vv/model.rs#ObservationManifest::source_ref_of,crates/fs-evidence/src/vv/model.rs#ObservationManifest::canonical_hash,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::try_new,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::id,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::qois,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::observation_ids,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::manifest,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::instrument_calibration,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::contains_clock,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::instruments,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::clocks,crates/fs-evidence/src/vv/model.rs#ExperimentArtifact::authenticity,crates/fs-evidence/src/vv/model.rs#DataAuthenticity::source_bytes_hash,crates/fs-evidence/src/vv/model.rs#DataAuthenticity::custody_receipt_hash,crates/fs-evidence/src/vv/model.rs#InstrumentCalibration::instrument_id,crates/fs-evidence/src/vv/model.rs#InstrumentCalibration::certificate_hash,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::try_new,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::id,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::experiment,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::preregistration_hash,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::calibration_ids,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::validation_ids,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::blind_sources,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::blind_commitment,crates/fs-evidence/src/vv/model.rs#CalibrationSplit::blind_selection,crates/fs-evidence/src/vv/model.rs#BlindReleaseReceipt::authority_receipt_hash,crates/fs-evidence/src/vv/model.rs#ArtifactRef::new,crates/fs-evidence/src/vv/model.rs#ClockSynchronization::contains_clock,crates/fs-evidence/src/vv/codec.rs#canonical_artifact_bytes,crates/fs-evidence/src/vv/codec.rs#content_hash_for,crates/fs-evidence/src/vv/codec.rs#encode_context,crates/fs-evidence/src/vv/codec.rs#decode_context,crates/fs-evidence/src/vv/codec.rs#encode_observation_source_ref,crates/fs-evidence/src/vv/codec.rs#decode_observation_source_ref,crates/fs-evidence/src/vv/codec.rs#encode_observation_manifest,crates/fs-evidence/src/vv/codec.rs#decode_observation_manifest,crates/fs-evidence/src/vv/codec.rs#encode_experiment,crates/fs-evidence/src/vv/codec.rs#decode_experiment,crates/fs-evidence/src/vv/codec.rs#encode_calibration_split,crates/fs-evidence/src/vv/codec.rs#decode_calibration_split,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-evidence:observation-manifest,fs-evidence:vv-artifact,fs-evidence:vv-blind-holdout,fs-material:identifiability-problem,fs-matdb:canonical-parameter-block",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=SourceAdmissionRecord,SourceResolution,SourceKind,AuthorityDisposition,TrustReceiptRef,TrustAuthentication,SourceVerification",
    "source_fields=SourceAdmissionRecord.schema_version:semantic,SourceAdmissionRecord.problem_id:semantic,SourceAdmissionRecord.resolutions:semantic,SourceResolution.key:derived:transitively-bound-by-parent-field,SourceResolution.kind:derived:transitively-bound-by-parent-field,SourceResolution.resolved_hash:derived:transitively-bound-by-parent-field,SourceResolution.content_hash_domain:derived:transitively-bound-by-parent-field,SourceResolution.contract_version:derived:transitively-bound-by-parent-field,SourceResolution.authority:derived:transitively-bound-by-parent-field,SourceResolution.verification:derived:transitively-bound-by-parent-field,SourceKind.variant:derived:transitively-bound-by-parent-field,AuthorityDisposition.variant:derived:transitively-bound-by-parent-field,AuthorityDisposition.trust_receipt:derived:transitively-bound-by-parent-field,AuthorityDisposition.reason:derived:transitively-bound-by-parent-field,TrustReceiptRef.receipt:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject_artifact:derived:transitively-bound-by-parent-field,TrustReceiptRef.authentication:derived:transitively-bound-by-parent-field,TrustAuthentication.variant:derived:transitively-bound-by-parent-field,TrustAuthentication.issuer:derived:transitively-bound-by-parent-field,TrustAuthentication.trust_policy:derived:transitively-bound-by-parent-field,SourceVerification.variant:derived:transitively-bound-by-parent-field,SourceVerification.byte_len:derived:transitively-bound-by-parent-field",
    "source_bindings=SourceAdmissionRecord.schema_version>wire-schema-version,SourceAdmissionRecord.problem_id>problem-id,SourceAdmissionRecord.resolutions>source-resolution-registry",
    "external_semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le",
    "semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,wire-schema-version,problem-id,source-resolution-registry",
    "excluded_fields=execution-configuration:belongs-to-execution-identity,assessment-claims-and-evidence:belongs-to-assessment-identity,caller-container-order:canonicalized-before-identity",
    "consumers=AdmittedIdentifiabilityProblem::resolve_and_admit,AdmittedIdentifiabilityProblem::source_admission_id,IdentifiabilityExecutionPlan::try_new",
    "mutations=identity-domain:crates/fs-material/tests/identifiability_authority.rs#identity_domains_and_wire_magics_are_stage_separated,identity-version:crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact,wire-magic:crates/fs-material/tests/identifiability_authority.rs#identity_domains_and_wire_magics_are_stage_separated,canonical-field-order:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,collection-count-u32-le:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,wire-schema-version:crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact,problem-id:crates/fs-material/tests/identifiability_authority.rs#identifiability_source_admission_identity_bindings_have_exact_mutation_evidence,source-resolution-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_source_admission_identity_bindings_have_exact_mutation_evidence",
    "nonsemantic_mutations=execution-configuration:crates/fs-material/tests/identifiability_authority.rs#source_admission_id_is_stable_across_execution_variants,assessment-claims-and-evidence:crates/fs-material/tests/identifiability_authority.rs#evidence_changes_assessment_not_problem_or_execution,caller-container-order:crates/fs-material/tests/identifiability_authority.rs#source_resolution_input_order_is_nonsemantic",
    "field_guard=classify_identifiability_source_admission_identity_fields",
    "transport_guard=AdmittedIdentifiabilityProblem::source_admission_canonical_bytes",
    "version_guard=crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact",
    "coupling_surface=fs-material:identifiability-source-admission",
];

/// Owner-local declaration for an exact, source-authorized execution identity.
#[allow(dead_code)]
pub const IDENTIFIABILITY_EXECUTION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-material:identifiability-execution",
    "version_const=IDENTIFIABILITY_EXECUTION_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-material.identifiability-execution.v3",
    "domain_const=IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN",
    "encoder=execution_identity_hash",
    "encoder_helpers=encode_execution_identity,encode_execution_with_header_mode",
    "schema_constants=IDENTIFIABILITY_EXECUTION_IDENTITY_VERSION,IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN,EXECUTION_MAGIC,MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ID_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_TEXT_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_CANONICAL_BYTES,crates/fs-qty/src/semantic.rs#QUANTITY_SPEC_ENCODING_VERSION,crates/fs-qty/src/semantic.rs#QUANTITY_SPEC_ENCODED_LEN",
    "schema_functions=IdentifiabilityExecutionPlan::try_new,IdentifiabilityExecutionPlan::canonical_bytes,IdentifiabilityExecutionPlan::from_canonical_bytes,DimensionlessErrorPolicy::try_new,IdentifiabilityNumericalPolicy::try_new,classify_identifiability_execution_identity_fields,validate_execution_structural_budget,validate_gauge_reduction_dag,gauge_reduction_sources,gauge_reduction_stage_sources,gauge_measure_sources,gauge_reduction_precedes,validate_gauge_slice,reduction_uses_regular_atlas,gauge_action_view,gauge_cell_for_claim,gauge_axes_match_claim,claim_gauge_action,gauge_action_status_and_geometry,gauge_action_orbit_profile,classify_identifiability_claim_tuple_schema,required_axes,validate_claim_source_kind,validate_claim_sources,is_free_inferential_parameter,claim_case_set,claim_subject_parameters,validate_claim_compatibility,validate_coordinate_for_parameter,bind_source_reference,validate_source_authority_closure,admit_opaque_resolution,encode_execution,decode_execution,check_execution_identity_version,encode_execution_action,decode_execution_action,encode_gauge_action_reference,decode_gauge_action_reference,encode_gauge_slice,decode_gauge_slice,encode_gauge_quotient,decode_gauge_quotient,encode_continuous_gauge_reduction,decode_continuous_gauge_reduction,encode_gauge_reduction_stage,decode_gauge_reduction_stage,encode_gauge_measure_semantics,decode_gauge_measure_semantics,encode_gauge_reduction_binding,decode_gauge_reduction_binding,encode_dimensionless_error_policy,decode_dimensionless_error_policy,encode_claim_request,decode_claim_request,encode_claim,decode_claim,encode_source_ref,decode_source_ref,encode_source_kind,decode_source_kind,encode_source_key,decode_source_key,encode_case_id,decode_case_id,encode_role,decode_role,encode_channel,decode_channel,encode_observation_key,decode_observation_key,encode_resolution_entry,encode_resolution_set,decode_resolution_set,encode_resolution_verification,decode_resolution_verification,crates/fs-material/src/identifiability.rs#encode_coordinate,crates/fs-material/src/identifiability.rs#decode_coordinate,crates/fs-material/src/identifiability.rs#ParameterCoordinate::try_new,crates/fs-material/src/identifiability.rs#ParameterCoordinate::id,crates/fs-material/src/identifiability.rs#CoordinateTransform::validate,crates/fs-material/src/identifiability.rs#CoordinateTransform::map,crates/fs-material/src/identifiability.rs#CoordinateTransform::mapped_domain,crates/fs-material/src/identifiability.rs#ParameterDomain::try_new,crates/fs-material/src/identifiability.rs#checked_add_dims,crates/fs-material/src/identifiability.rs#same_f64,crates/fs-material/src/identifiability.rs#validate_header_profile,crates/fs-material/src/identifiability.rs#encode_header,crates/fs-material/src/identifiability.rs#decode_header,crates/fs-material/src/identifiability.rs#canonical_f64,crates/fs-material/src/identifiability.rs#validate_token,crates/fs-material/src/identifiability.rs#validate_reason,crates/fs-material/src/identifiability.rs#CanonicalWriter::new,crates/fs-material/src/identifiability.rs#CanonicalWriter::raw,crates/fs-material/src/identifiability.rs#CanonicalWriter::byte,crates/fs-material/src/identifiability.rs#CanonicalWriter::u32,crates/fs-material/src/identifiability.rs#CanonicalWriter::u64,crates/fs-material/src/identifiability.rs#CanonicalWriter::f64,crates/fs-material/src/identifiability.rs#CanonicalWriter::count,crates/fs-material/src/identifiability.rs#CanonicalWriter::text,crates/fs-material/src/identifiability.rs#CanonicalWriter::hash,crates/fs-material/src/identifiability.rs#CanonicalWriter::quantity,crates/fs-material/src/identifiability.rs#CanonicalWriter::finish,crates/fs-material/src/identifiability.rs#CanonicalReader::new,crates/fs-material/src/identifiability.rs#CanonicalReader::take,crates/fs-material/src/identifiability.rs#CanonicalReader::byte,crates/fs-material/src/identifiability.rs#CanonicalReader::u32,crates/fs-material/src/identifiability.rs#CanonicalReader::u64,crates/fs-material/src/identifiability.rs#CanonicalReader::f64,crates/fs-material/src/identifiability.rs#CanonicalReader::length,crates/fs-material/src/identifiability.rs#CanonicalReader::count,crates/fs-material/src/identifiability.rs#CanonicalReader::text,crates/fs-material/src/identifiability.rs#CanonicalReader::token,crates/fs-material/src/identifiability.rs#CanonicalReader::reason,crates/fs-material/src/identifiability.rs#CanonicalReader::hash,crates/fs-material/src/identifiability.rs#CanonicalReader::quantity,crates/fs-material/src/identifiability.rs#CanonicalReader::expect_byte,crates/fs-material/src/identifiability.rs#CanonicalReader::finish,crates/fs-blake3/src/lib.rs#hash_domain,crates/fs-qty/src/semantic.rs#QuantitySpec::canonical_bytes,crates/fs-qty/src/semantic.rs#QuantitySpec::from_canonical_bytes",
    "schema_dependencies=fs-material:identifiability-problem,fs-material:identifiability-source-admission",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=IdentifiabilityExecutionPlan,SourceRef,SourceKind,SourceResolutionSet,SourceResolution,AuthorityDisposition,TrustReceiptRef,TrustAuthentication,SourceVerification,ParameterExecutionAction,ArithmeticPolicy,DimensionlessErrorPolicy,ClaimRequest,IdentifiabilityNumericalPolicy,InformationRegime,IdentifiabilityExtent,ScalarDomain,FiberDimensionLowerBound,FiberStructure,ClaimQuantifier,TypedIdentifiabilityClaim,GaugeSliceCodimension,GaugeSlicePlan,GaugeQuotientPlan,ContinuousGaugeReductionPlan,GaugeReductionPlan,GaugeReductionStageRelation,GaugeReductionStage,GaugeMeasureSemantics,GaugeReductionBinding",
    "source_fields=IdentifiabilityExecutionPlan.schema_version:semantic,IdentifiabilityExecutionPlan.header:derived:identity-projection-excludes-artifact-id,IdentifiabilityExecutionPlan.problem_id:semantic,IdentifiabilityExecutionPlan.source_admission_id:semantic,IdentifiabilityExecutionPlan.analyzer:semantic,IdentifiabilityExecutionPlan.build:semantic,IdentifiabilityExecutionPlan.derivative_provider:semantic,IdentifiabilityExecutionPlan.claim_requests:semantic,IdentifiabilityExecutionPlan.actions:semantic,IdentifiabilityExecutionPlan.gauge_reductions:semantic,IdentifiabilityExecutionPlan.numerical:semantic,IdentifiabilityExecutionPlan.initialization:semantic,IdentifiabilityExecutionPlan.stopping:semantic,IdentifiabilityExecutionPlan.determinism_contract:semantic,IdentifiabilityExecutionPlan.source_authority:semantic,SourceRef.key:derived:transitively-bound-by-parent-field,SourceRef.kind:derived:transitively-bound-by-parent-field,SourceRef.expected_hash:derived:transitively-bound-by-parent-field,SourceRef.content_hash_domain:derived:transitively-bound-by-parent-field,SourceRef.contract_version:derived:transitively-bound-by-parent-field,SourceKind.variant:derived:transitively-bound-by-parent-field,SourceResolutionSet.entries:derived:transitively-bound-by-parent-field,SourceResolution.key:derived:transitively-bound-by-parent-field,SourceResolution.kind:derived:transitively-bound-by-parent-field,SourceResolution.resolved_hash:derived:transitively-bound-by-parent-field,SourceResolution.content_hash_domain:derived:transitively-bound-by-parent-field,SourceResolution.contract_version:derived:transitively-bound-by-parent-field,SourceResolution.authority:derived:transitively-bound-by-parent-field,SourceResolution.verification:derived:transitively-bound-by-parent-field,AuthorityDisposition.variant:derived:transitively-bound-by-parent-field,AuthorityDisposition.trust_receipt:derived:transitively-bound-by-parent-field,AuthorityDisposition.reason:derived:transitively-bound-by-parent-field,TrustReceiptRef.receipt:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject_artifact:derived:transitively-bound-by-parent-field,TrustReceiptRef.authentication:derived:transitively-bound-by-parent-field,TrustAuthentication.variant:derived:transitively-bound-by-parent-field,TrustAuthentication.issuer:derived:transitively-bound-by-parent-field,TrustAuthentication.trust_policy:derived:transitively-bound-by-parent-field,SourceVerification.variant:derived:transitively-bound-by-parent-field,SourceVerification.byte_len:derived:transitively-bound-by-parent-field,ParameterExecutionAction.variant:derived:transitively-bound-by-parent-field,ParameterExecutionAction.coordinate:derived:transitively-bound-by-parent-field,ParameterExecutionAction.integrator:derived:transitively-bound-by-parent-field,ParameterExecutionAction.measure_transport:derived:transitively-bound-by-parent-field,ArithmeticPolicy.variant:derived:transitively-bound-by-parent-field,DimensionlessErrorPolicy.metric:derived:transitively-bound-by-parent-field,DimensionlessErrorPolicy.nondimensionalization:derived:transitively-bound-by-parent-field,DimensionlessErrorPolicy.maximum_certified_error:derived:transitively-bound-by-parent-field,ClaimRequest.claim:derived:transitively-bound-by-parent-field,ClaimRequest.error_policy:derived:transitively-bound-by-parent-field,IdentifiabilityNumericalPolicy.rank_tolerance:derived:transitively-bound-by-parent-field,IdentifiabilityNumericalPolicy.singular_value_floor:derived:transitively-bound-by-parent-field,IdentifiabilityNumericalPolicy.maximum_condition_number:derived:transitively-bound-by-parent-field,IdentifiabilityNumericalPolicy.arithmetic:derived:transitively-bound-by-parent-field,IdentifiabilityNumericalPolicy.nondimensionalization:derived:transitively-bound-by-parent-field,InformationRegime.variant:derived:transitively-bound-by-parent-field,InformationRegime.joint_prior:derived:transitively-bound-by-parent-field,IdentifiabilityExtent.variant:derived:transitively-bound-by-parent-field,ScalarDomain.variant:derived:transitively-bound-by-parent-field,ScalarDomain.extension:derived:transitively-bound-by-parent-field,ScalarDomain.stratification:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.variant:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.minimum_dimension:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.model_space:derived:transitively-bound-by-parent-field,FiberStructure.variant:derived:transitively-bound-by-parent-field,FiberStructure.maximum_cardinality:derived:transitively-bound-by-parent-field,FiberStructure.action:derived:transitively-bound-by-parent-field,FiberStructure.equivalence:derived:transitively-bound-by-parent-field,FiberStructure.lower_bound:derived:transitively-bound-by-parent-field,FiberStructure.orbit_type_profile:derived:transitively-bound-by-parent-field,FiberStructure.strata:derived:transitively-bound-by-parent-field,ClaimQuantifier.variant:derived:transitively-bound-by-parent-field,ClaimQuantifier.realization:derived:transitively-bound-by-parent-field,ClaimQuantifier.measure:derived:transitively-bound-by-parent-field,ClaimQuantifier.domain:derived:transitively-bound-by-parent-field,ClaimQuantifier.probability:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.id:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.information:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.extent:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.fiber:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.quantifier:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.scalar_domain:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.subject:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.scope:derived:transitively-bound-by-parent-field,GaugeSliceCodimension.variant:derived:transitively-bound-by-parent-field,GaugeSliceCodimension.codimension:derived:transitively-bound-by-parent-field,GaugeSliceCodimension.codimension_model:derived:transitively-bound-by-parent-field,GaugeSliceCodimension.compatibility:derived:transitively-bound-by-parent-field,GaugeSliceCodimension.profile:derived:transitively-bound-by-parent-field,GaugeSlicePlan.support:derived:transitively-bound-by-parent-field,GaugeSlicePlan.constraint:derived:transitively-bound-by-parent-field,GaugeSlicePlan.expected_codimension:derived:transitively-bound-by-parent-field,GaugeSlicePlan.coverage:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.variant:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.quotient_map:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.local_section_atlas:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.coverage:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.quotient_profile:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.local_models:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.invariants:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.completeness_profile:derived:transitively-bound-by-parent-field,GaugeQuotientPlan.presentation:derived:transitively-bound-by-parent-field,ContinuousGaugeReductionPlan.variant:derived:transitively-bound-by-parent-field,ContinuousGaugeReductionPlan.quotient:derived:transitively-bound-by-parent-field,ContinuousGaugeReductionPlan.slice:derived:transitively-bound-by-parent-field,GaugeReductionPlan.variant:derived:transitively-bound-by-parent-field,GaugeReductionPlan.reason:derived:transitively-bound-by-parent-field,GaugeReductionPlan.quotient:derived:transitively-bound-by-parent-field,GaugeReductionPlan.slice:derived:transitively-bound-by-parent-field,GaugeReductionPlan.reduction:derived:transitively-bound-by-parent-field,GaugeReductionPlan.normal_subgroup:derived:transitively-bound-by-parent-field,GaugeReductionPlan.factor_extension:derived:transitively-bound-by-parent-field,GaugeReductionPlan.residual_quotient_action:derived:transitively-bound-by-parent-field,GaugeReductionPlan.compatibility:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.variant:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.normality:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.induced_residual_action:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.extension:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.transversality:derived:transitively-bound-by-parent-field,GaugeReductionStageRelation.reducibility:derived:transitively-bound-by-parent-field,GaugeReductionStage.variant:derived:transitively-bound-by-parent-field,GaugeReductionStage.predecessors:derived:transitively-bound-by-parent-field,GaugeReductionStage.composition_law:derived:transitively-bound-by-parent-field,GaugeReductionStage.relation:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.variant:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.reason:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.source_measure:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.reduced_measure:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.transport:derived:transitively-bound-by-parent-field,GaugeMeasureSemantics.jacobian_or_disintegration:derived:transitively-bound-by-parent-field,GaugeReductionBinding.id:derived:transitively-bound-by-parent-field,GaugeReductionBinding.action:derived:transitively-bound-by-parent-field,GaugeReductionBinding.claims:derived:transitively-bound-by-parent-field,GaugeReductionBinding.plan:derived:transitively-bound-by-parent-field,GaugeReductionBinding.stage:derived:transitively-bound-by-parent-field,GaugeReductionBinding.measure:derived:transitively-bound-by-parent-field",
    "source_bindings=IdentifiabilityExecutionPlan.schema_version>wire-schema-version,IdentifiabilityExecutionPlan.problem_id>problem-id,IdentifiabilityExecutionPlan.source_admission_id>source-admission-id,IdentifiabilityExecutionPlan.analyzer>analyzer-source,IdentifiabilityExecutionPlan.build>build-source,IdentifiabilityExecutionPlan.derivative_provider>derivative-provider-source,IdentifiabilityExecutionPlan.claim_requests>claim-requests,IdentifiabilityExecutionPlan.actions>parameter-execution-actions,IdentifiabilityExecutionPlan.gauge_reductions>gauge-reduction-plans,IdentifiabilityExecutionPlan.numerical>numerical-policy,IdentifiabilityExecutionPlan.initialization>initialization-source,IdentifiabilityExecutionPlan.stopping>stopping-source,IdentifiabilityExecutionPlan.determinism_contract>determinism-contract-source,IdentifiabilityExecutionPlan.source_authority>execution-source-authority",
    "external_semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian,identity-header-projection-marker,header-units,header-seed,header-accuracy,header-time-ms,header-memory-bytes,header-versions,header-capabilities",
    "semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian,identity-header-projection-marker,header-units,header-seed,header-accuracy,header-time-ms,header-memory-bytes,header-versions,header-capabilities,wire-schema-version,problem-id,source-admission-id,analyzer-source,build-source,derivative-provider-source,claim-requests,parameter-execution-actions,gauge-reduction-plans,numerical-policy,initialization-source,stopping-source,determinism-contract-source,execution-source-authority",
    "excluded_fields=ArtifactHeader.id:ledger-label-not-scientific-identity,assessment-claims-and-evidence:belongs-to-assessment-identity,caller-container-order:canonicalized-before-identity",
    "consumers=IdentifiabilityExecutionPlan::try_new,IdentifiabilityExecutionPlan::id,IdentifiabilityAssessment::try_new",
    "mutations=identity-domain:crates/fs-material/tests/identifiability_authority.rs#identity_domains_and_wire_magics_are_stage_separated,identity-version:crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact,wire-magic:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,canonical-field-order:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,collection-count-u32-le:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,fixed-numeric-little-endian:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,identity-header-projection-marker:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,header-units:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-seed:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-accuracy:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-time-ms:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-memory-bytes:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-versions:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,header-capabilities:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,wire-schema-version:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,problem-id:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,source-admission-id:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,analyzer-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,build-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,derivative-provider-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,claim-requests:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,parameter-execution-actions:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,gauge-reduction-plans:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,numerical-policy:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,initialization-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,stopping-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,determinism-contract-source:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence,execution-source-authority:crates/fs-material/tests/identifiability_authority.rs#identifiability_execution_identity_bindings_have_exact_mutation_evidence",
    "nonsemantic_mutations=ArtifactHeader.id:crates/fs-material/tests/identifiability_authority.rs#artifact_labels_do_not_move_execution_or_assessment_identity,assessment-claims-and-evidence:crates/fs-material/tests/identifiability_authority.rs#evidence_changes_assessment_not_problem_or_execution,caller-container-order:crates/fs-material/tests/identifiability_authority.rs#execution_action_input_order_is_nonsemantic",
    "field_guard=classify_identifiability_execution_identity_fields",
    "transport_guard=IdentifiabilityExecutionPlan::from_canonical_bytes",
    "version_guard=crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact",
    "coupling_surface=fs-material:identifiability-execution",
];

/// Owner-local declaration for exact claim/evidence assessment identity.
#[allow(dead_code)]
pub const IDENTIFIABILITY_ASSESSMENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-material:identifiability-assessment",
    "version_const=IDENTIFIABILITY_ASSESSMENT_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-material.identifiability-assessment.v3",
    "domain_const=IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN",
    "encoder=assessment_identity_hash",
    "encoder_helpers=encode_assessment_identity,encode_assessment_with_header_mode",
    "schema_constants=IDENTIFIABILITY_ASSESSMENT_IDENTITY_VERSION,IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN,ASSESSMENT_MAGIC,MAX_IDENTIFIABILITY_STRUCTURAL_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ID_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_TEXT_BYTES,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_ITEMS,crates/fs-material/src/identifiability.rs#MAX_IDENTIFIABILITY_CANONICAL_BYTES",
    "schema_functions=IdentifiabilityAssessment::try_new,IdentifiabilityAssessment::canonical_bytes,IdentifiabilityAssessment::from_canonical_bytes,classify_identifiability_assessment_identity_fields,validate_assessment_structural_budget,validate_positive_claim_influence_routes,validate_decisive_practical_claim_closure,relevant_gauge_obstructions,matching_gauge_cell,disposition_is_valid,claim_extent_obstruction,effective_gauge_orbit_kind,gauge_axes_match_claim,claim_gauge_action,classify_identifiability_claim_tuple_schema,validate_claim_source_kind,validate_claim_sources,is_free_inferential_parameter,claim_case_set,claim_subject_parameters,validate_claim_compatibility,functional_observations,require_source_kind_in,admit_opaque_resolution,bind_source_reference,validate_source_authority_closure,encode_assessment,decode_assessment,check_assessment_identity_version,encode_claim,decode_claim,encode_claim_assessment,decode_claim_assessment,decode_optional_source_ref,encode_source_ref,decode_source_ref,encode_source_kind,decode_source_kind,encode_source_key,decode_source_key,encode_role,decode_role,encode_case_id,decode_case_id,encode_channel,decode_channel,encode_observation_key,decode_observation_key,encode_resolution_entry,encode_resolution_set,decode_resolution_set,encode_resolution_verification,decode_resolution_verification,crates/fs-material/src/identifiability.rs#validate_header_profile,crates/fs-material/src/identifiability.rs#encode_header,crates/fs-material/src/identifiability.rs#decode_header,crates/fs-material/src/identifiability.rs#canonical_f64,crates/fs-material/src/identifiability.rs#validate_token,crates/fs-material/src/identifiability.rs#validate_reason,crates/fs-material/src/identifiability.rs#CanonicalWriter::new,crates/fs-material/src/identifiability.rs#CanonicalWriter::raw,crates/fs-material/src/identifiability.rs#CanonicalWriter::byte,crates/fs-material/src/identifiability.rs#CanonicalWriter::u32,crates/fs-material/src/identifiability.rs#CanonicalWriter::u64,crates/fs-material/src/identifiability.rs#CanonicalWriter::f64,crates/fs-material/src/identifiability.rs#CanonicalWriter::count,crates/fs-material/src/identifiability.rs#CanonicalWriter::text,crates/fs-material/src/identifiability.rs#CanonicalWriter::hash,crates/fs-material/src/identifiability.rs#CanonicalWriter::finish,crates/fs-material/src/identifiability.rs#CanonicalReader::new,crates/fs-material/src/identifiability.rs#CanonicalReader::take,crates/fs-material/src/identifiability.rs#CanonicalReader::byte,crates/fs-material/src/identifiability.rs#CanonicalReader::u32,crates/fs-material/src/identifiability.rs#CanonicalReader::u64,crates/fs-material/src/identifiability.rs#CanonicalReader::f64,crates/fs-material/src/identifiability.rs#CanonicalReader::length,crates/fs-material/src/identifiability.rs#CanonicalReader::count,crates/fs-material/src/identifiability.rs#CanonicalReader::text,crates/fs-material/src/identifiability.rs#CanonicalReader::token,crates/fs-material/src/identifiability.rs#CanonicalReader::reason,crates/fs-material/src/identifiability.rs#CanonicalReader::hash,crates/fs-material/src/identifiability.rs#CanonicalReader::finish,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-material:identifiability-problem,fs-material:identifiability-execution",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=IdentifiabilityAssessment,SourceRef,SourceKind,SourceResolutionSet,SourceResolution,AuthorityDisposition,TrustReceiptRef,TrustAuthentication,SourceVerification,InformationRegime,IdentifiabilityExtent,ScalarDomain,FiberDimensionLowerBound,FiberStructure,ClaimQuantifier,TypedIdentifiabilityClaim,GaugeResolutionDisposition,GaugeResolutionEvidence,ClaimAssessment",
    "source_fields=IdentifiabilityAssessment.schema_version:semantic,IdentifiabilityAssessment.header:derived:identity-projection-excludes-artifact-id,IdentifiabilityAssessment.problem_id:semantic,IdentifiabilityAssessment.execution_id:semantic,IdentifiabilityAssessment.claims:semantic,IdentifiabilityAssessment.evidence:semantic,IdentifiabilityAssessment.source_authority:semantic,SourceRef.key:derived:transitively-bound-by-parent-field,SourceRef.kind:derived:transitively-bound-by-parent-field,SourceRef.expected_hash:derived:transitively-bound-by-parent-field,SourceRef.content_hash_domain:derived:transitively-bound-by-parent-field,SourceRef.contract_version:derived:transitively-bound-by-parent-field,SourceKind.variant:derived:transitively-bound-by-parent-field,SourceResolutionSet.entries:derived:transitively-bound-by-parent-field,SourceResolution.key:derived:transitively-bound-by-parent-field,SourceResolution.kind:derived:transitively-bound-by-parent-field,SourceResolution.resolved_hash:derived:transitively-bound-by-parent-field,SourceResolution.content_hash_domain:derived:transitively-bound-by-parent-field,SourceResolution.contract_version:derived:transitively-bound-by-parent-field,SourceResolution.authority:derived:transitively-bound-by-parent-field,SourceResolution.verification:derived:transitively-bound-by-parent-field,AuthorityDisposition.variant:derived:transitively-bound-by-parent-field,AuthorityDisposition.trust_receipt:derived:transitively-bound-by-parent-field,AuthorityDisposition.reason:derived:transitively-bound-by-parent-field,TrustReceiptRef.receipt:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject:derived:transitively-bound-by-parent-field,TrustReceiptRef.subject_artifact:derived:transitively-bound-by-parent-field,TrustReceiptRef.authentication:derived:transitively-bound-by-parent-field,TrustAuthentication.variant:derived:transitively-bound-by-parent-field,TrustAuthentication.issuer:derived:transitively-bound-by-parent-field,TrustAuthentication.trust_policy:derived:transitively-bound-by-parent-field,SourceVerification.variant:derived:transitively-bound-by-parent-field,SourceVerification.byte_len:derived:transitively-bound-by-parent-field,InformationRegime.variant:derived:transitively-bound-by-parent-field,InformationRegime.joint_prior:derived:transitively-bound-by-parent-field,IdentifiabilityExtent.variant:derived:transitively-bound-by-parent-field,ScalarDomain.variant:derived:transitively-bound-by-parent-field,ScalarDomain.extension:derived:transitively-bound-by-parent-field,ScalarDomain.stratification:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.variant:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.minimum_dimension:derived:transitively-bound-by-parent-field,FiberDimensionLowerBound.model_space:derived:transitively-bound-by-parent-field,FiberStructure.variant:derived:transitively-bound-by-parent-field,FiberStructure.maximum_cardinality:derived:transitively-bound-by-parent-field,FiberStructure.action:derived:transitively-bound-by-parent-field,FiberStructure.equivalence:derived:transitively-bound-by-parent-field,FiberStructure.lower_bound:derived:transitively-bound-by-parent-field,FiberStructure.orbit_type_profile:derived:transitively-bound-by-parent-field,FiberStructure.strata:derived:transitively-bound-by-parent-field,ClaimQuantifier.variant:derived:transitively-bound-by-parent-field,ClaimQuantifier.realization:derived:transitively-bound-by-parent-field,ClaimQuantifier.measure:derived:transitively-bound-by-parent-field,ClaimQuantifier.domain:derived:transitively-bound-by-parent-field,ClaimQuantifier.probability:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.id:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.information:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.extent:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.fiber:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.quantifier:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.scalar_domain:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.subject:derived:transitively-bound-by-parent-field,TypedIdentifiabilityClaim.scope:derived:transitively-bound-by-parent-field,GaugeResolutionDisposition.variant:derived:transitively-bound-by-parent-field,GaugeResolutionEvidence.action:derived:transitively-bound-by-parent-field,GaugeResolutionEvidence.disposition:derived:transitively-bound-by-parent-field,GaugeResolutionEvidence.method:derived:transitively-bound-by-parent-field,GaugeResolutionEvidence.receipt:derived:transitively-bound-by-parent-field,ClaimAssessment.variant:derived:transitively-bound-by-parent-field,ClaimAssessment.method:derived:transitively-bound-by-parent-field,ClaimAssessment.receipt:derived:transitively-bound-by-parent-field,ClaimAssessment.metric:derived:transitively-bound-by-parent-field,ClaimAssessment.nondimensionalization:derived:transitively-bound-by-parent-field,ClaimAssessment.certified_error_bound:derived:transitively-bound-by-parent-field,ClaimAssessment.gauge_resolutions:derived:transitively-bound-by-parent-field,ClaimAssessment.reason:derived:transitively-bound-by-parent-field",
    "source_bindings=IdentifiabilityAssessment.schema_version>wire-schema-version,IdentifiabilityAssessment.problem_id>problem-id,IdentifiabilityAssessment.execution_id>execution-id,IdentifiabilityAssessment.claims>typed-claim-registry,IdentifiabilityAssessment.evidence>claim-assessment-registry,IdentifiabilityAssessment.source_authority>assessment-source-authority",
    "external_semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian,identity-header-projection-marker,header-units,header-seed,header-accuracy,header-time-ms,header-memory-bytes,header-versions,header-capabilities",
    "semantic_fields=identity-domain,identity-version,wire-magic,canonical-field-order,collection-count-u32-le,fixed-numeric-little-endian,identity-header-projection-marker,header-units,header-seed,header-accuracy,header-time-ms,header-memory-bytes,header-versions,header-capabilities,wire-schema-version,problem-id,execution-id,typed-claim-registry,claim-assessment-registry,assessment-source-authority",
    "excluded_fields=ArtifactHeader.id:ledger-label-not-scientific-identity,caller-container-order:canonicalized-before-identity",
    "consumers=IdentifiabilityAssessment::try_new,IdentifiabilityAssessment::id",
    "mutations=identity-domain:crates/fs-material/tests/identifiability_authority.rs#identity_domains_and_wire_magics_are_stage_separated,identity-version:crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact,wire-magic:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,canonical-field-order:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,collection-count-u32-le:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,fixed-numeric-little-endian:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,identity-header-projection-marker:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_preimages_have_exact_wire_layout,header-units:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-seed:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-accuracy:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-time-ms:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-memory-bytes:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-versions:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,header-capabilities:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,wire-schema-version:crates/fs-material/tests/identifiability_authority.rs#identifiability_identity_versions_and_transports_fail_closed,problem-id:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,execution-id:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,typed-claim-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,claim-assessment-registry:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence,assessment-source-authority:crates/fs-material/tests/identifiability_authority.rs#identifiability_assessment_identity_bindings_have_exact_mutation_evidence",
    "nonsemantic_mutations=ArtifactHeader.id:crates/fs-material/tests/identifiability_authority.rs#artifact_labels_do_not_move_execution_or_assessment_identity,caller-container-order:crates/fs-material/tests/identifiability_authority.rs#assessment_input_order_is_nonsemantic",
    "field_guard=classify_identifiability_assessment_identity_fields",
    "transport_guard=IdentifiabilityAssessment::from_canonical_bytes",
    "version_guard=crates/fs-material/tests/identifiability_authority.rs#identity_version_guard_is_exact",
    "coupling_surface=fs-material:identifiability-assessment",
];

fn required_axes(claim: &TypedIdentifiabilityClaim) -> BTreeSet<RequestedClaimAxis> {
    let mut required = BTreeSet::new();
    match &claim.information {
        InformationRegime::StructuralExactModel | InformationRegime::ExactInputOutputMap => {
            required.insert(RequestedClaimAxis::Structural);
        }
        InformationRegime::NoisyFiniteData
        | InformationRegime::PosteriorUnderDeclaredPrior { .. } => {
            required.insert(RequestedClaimAxis::Practical);
        }
    }
    required.insert(match claim.extent {
        IdentifiabilityExtent::Local => RequestedClaimAxis::Local,
        IdentifiabilityExtent::Global => RequestedClaimAxis::Global,
    });
    if matches!(claim.quantifier, ClaimQuantifier::AlmostEverywhere { .. }) {
        required.insert(RequestedClaimAxis::Generic);
    }
    required
}

fn validate_claim_source_kind(
    source: &SourceRef,
    allowed: &[SourceKind],
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    if !allowed.contains(&source.kind) {
        return Err(IdentifiabilityError::InvalidText {
            field,
            detail: format!(
                "source {} has kind {:?}, expected one of {allowed:?}",
                source.key, source.kind
            ),
        });
    }
    Ok(())
}

fn validate_claim_sources<'a>(
    claim: &'a TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<Vec<&'a SourceRef>, IdentifiabilityError> {
    let mut sources = Vec::new();
    match &claim.information {
        InformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            validate_claim_source_kind(
                joint_prior,
                &[SourceKind::ProbabilityMeasure],
                "claim joint-prior measure",
            )?;
            let declared_key = problem.joint_prior.as_ref().ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "problem joint prior",
                    id: joint_prior.key.to_string(),
                }
            })?;
            if problem.sources.get(declared_key) != Some(joint_prior) {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "claim joint prior/problem joint prior",
                });
            }
            sources.push(joint_prior);
        }
        InformationRegime::StructuralExactModel
        | InformationRegime::ExactInputOutputMap
        | InformationRegime::NoisyFiniteData => {}
    }
    match &claim.scalar_domain {
        ScalarDomain::Real => {}
        ScalarDomain::Complex { extension } => {
            validate_claim_source_kind(
                extension,
                &[SourceKind::AlgebraicExtension],
                "claim complex scalar extension",
            )?;
            sources.push(extension);
        }
        ScalarDomain::MixedDiscreteContinuous { stratification } => {
            validate_claim_source_kind(
                stratification,
                &[SourceKind::Stratification],
                "claim mixed scalar stratification",
            )?;
            sources.push(stratification);
        }
    }
    match &claim.fiber {
        FiberStructure::Stratified { strata } => {
            validate_claim_source_kind(
                strata,
                &[SourceKind::Stratification],
                "claim fiber stratification",
            )?;
            sources.push(strata);
        }
        FiberStructure::FiniteToOne {
            maximum_cardinality: Some(FiberCardinalityBound::SymbolicProfile(profile)),
        } => {
            validate_claim_source_kind(
                profile,
                &[SourceKind::FiberCardinalityProfile],
                "claim fiber cardinality profile",
            )?;
            sources.push(profile);
        }
        FiberStructure::PositiveDimensional {
            lower_bound: FiberDimensionLowerBound::InfiniteDimensional { model_space },
        } => {
            validate_claim_source_kind(
                model_space,
                &[SourceKind::FiberDimensionProfile],
                "claim infinite-dimensional fiber profile",
            )?;
            sources.push(model_space);
        }
        FiberStructure::GeneralizedQuotientUnique { equivalence, .. } => {
            validate_claim_source_kind(
                equivalence,
                &[SourceKind::GaugeQuotientProfile],
                "claim generalized quotient equivalence",
            )?;
            sources.push(equivalence);
        }
        FiberStructure::StratifiedOrbit {
            orbit_type_profile, ..
        } => {
            validate_claim_source_kind(
                orbit_type_profile,
                &[SourceKind::GaugeOrbitTypeProfile],
                "claim stratified-orbit profile",
            )?;
            sources.push(orbit_type_profile);
        }
        FiberStructure::Unique
        | FiberStructure::FiniteToOne { .. }
        | FiberStructure::DiscreteOrbit { .. }
        | FiberStructure::MixedOrbit { .. }
        | FiberStructure::OrbitQuotientUnique { .. }
        | FiberStructure::PositiveDimensional { .. } => {}
    }
    if let ClaimSubject::DerivedFunctional { definition, .. } = &claim.subject {
        validate_claim_source_kind(
            definition,
            &[SourceKind::DerivedFunctional],
            "claim derived-functional definition",
        )?;
        sources.push(definition);
    }
    match &claim.quantifier {
        ClaimQuantifier::AtRealization { realization } => {
            validate_claim_source_kind(
                realization,
                &[SourceKind::QuantifierRealization],
                "claim realization source",
            )?;
            sources.push(realization);
        }
        ClaimQuantifier::AlmostEverywhere { measure } => {
            validate_claim_source_kind(
                measure,
                &[SourceKind::ReferenceMeasure, SourceKind::ProbabilityMeasure],
                "claim almost-everywhere measure",
            )?;
            sources.push(measure);
        }
        ClaimQuantifier::ForAll { domain } => {
            validate_claim_source_kind(
                domain,
                &[SourceKind::QuantifierDomain],
                "claim universal domain",
            )?;
            sources.push(domain);
        }
        ClaimQuantifier::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            if !probability.is_finite() || *probability <= 0.0 || *probability > 1.0 {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "claim probability",
                    detail: "probability must lie in (0,1]".to_string(),
                });
            }
            validate_claim_source_kind(
                measure,
                &[SourceKind::ProbabilityMeasure],
                "claim probability measure",
            )?;
            sources.push(measure);
        }
    }
    if let ClaimScope::Stratum { definition, .. } = &claim.scope {
        validate_claim_source_kind(definition, &[SourceKind::Stratification], "claim stratum")?;
        if problem.sources.get(definition.key()) != Some(definition) {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "claim stratum/problem source",
            });
        }
        sources.push(definition);
    }
    Ok(sources)
}

fn is_free_inferential_parameter(parameter: &StudyParameter) -> bool {
    matches!(
        &parameter.treatment,
        ParameterTreatment::Estimated
            | ParameterTreatment::Profiled
            | ParameterTreatment::Marginalized
    )
}

fn claim_case_set(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<BTreeSet<CaseId>, IdentifiabilityError> {
    let cases = match &claim.scope {
        ClaimScope::WholeCampaign => problem.cases.keys().cloned().collect(),
        ClaimScope::Cases(cases) => cases.clone(),
        ClaimScope::Stratum { definition, cases } => {
            if definition.kind != SourceKind::Stratification
                || problem.sources.get(definition.key()) != Some(definition)
            {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "claim stratum/problem source",
                });
            }
            cases.clone()
        }
    };
    if cases.is_empty() {
        return Err(IdentifiabilityError::Cardinality {
            field: "claim case scope",
            detail: "claim case scope cannot be empty".to_string(),
        });
    }
    for case in &cases {
        if !problem.cases.contains_key(case) {
            return Err(IdentifiabilityError::UnknownReference {
                field: "claim case scope",
                id: case.to_string(),
            });
        }
    }
    Ok(cases)
}

fn claim_subject_parameters(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<BTreeSet<ParameterRoleId>, IdentifiabilityError> {
    let roles = match &claim.subject {
        ClaimSubject::Parameter(role) => BTreeSet::from([role.clone()]),
        ClaimSubject::ParameterSet(roles) => {
            if roles.len() < 2 {
                return Err(IdentifiabilityError::Cardinality {
                    field: "claim parameter set",
                    detail: "set claims need at least two parameters".to_string(),
                });
            }
            roles.clone()
        }
        ClaimSubject::DerivedFunctional { parameters, .. } => {
            if parameters.is_empty() {
                return Err(IdentifiabilityError::Cardinality {
                    field: "derived-functional parameters",
                    detail: "a derived functional must declare its nonempty parameter support"
                        .to_string(),
                });
            }
            parameters.clone()
        }
        ClaimSubject::Influence(influence) => {
            transitive_influence_ids(influence, &problem.influences)?
                .into_iter()
                .map(|id| problem.influences[&id].parameter.clone())
                .collect()
        }
        ClaimSubject::GaugeAction(action) => match action {
            GaugeActionReference::Single(gauge) => problem
                .gauges
                .get(gauge)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field: "claim gauge action",
                    id: gauge.to_string(),
                })?
                .members
                .clone(),
            GaugeActionReference::Product(composition)
            | GaugeActionReference::Composition(composition) => problem
                .gauge_compositions
                .get(composition)
                .ok_or_else(|| IdentifiabilityError::UnknownReference {
                    field: "claim gauge composition",
                    id: composition.to_string(),
                })?
                .members
                .iter()
                .flat_map(|member| problem.gauges[member].members.iter().cloned())
                .collect(),
        },
        ClaimSubject::WholeProblem => problem
            .parameters
            .iter()
            .filter(|(_, parameter)| is_free_inferential_parameter(parameter))
            .map(|(role, _)| role.clone())
            .collect(),
    };
    for role in &roles {
        if !problem.parameters.contains_key(role) {
            return Err(IdentifiabilityError::UnknownReference {
                field: "claim parameter support",
                id: role.to_string(),
            });
        }
    }
    match &claim.subject {
        ClaimSubject::Parameter(_) | ClaimSubject::ParameterSet(_) => {
            if let Some(nonfree) = roles
                .iter()
                .find(|role| !is_free_inferential_parameter(&problem.parameters[*role]))
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "raw identifiability claim subject",
                    detail: format!(
                        "raw parameter subject {nonfree} is conditioned or derived; known-by-conditioning and derived-functional propositions are distinct from identifiability of free inferential coordinates"
                    ),
                });
            }
        }
        ClaimSubject::DerivedFunctional { .. } => {
            if !roles
                .iter()
                .any(|role| is_free_inferential_parameter(&problem.parameters[role]))
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "derived-functional identifiability subject",
                    detail: "a derived-functional identifiability claim must depend on at least one free inferential parameter"
                        .to_string(),
                });
            }
        }
        ClaimSubject::Influence(_) | ClaimSubject::GaugeAction(_) | ClaimSubject::WholeProblem => {}
    }
    Ok(roles)
}

fn parameter_active_in_cases(parameter: &StudyParameter, cases: &BTreeSet<CaseId>) -> bool {
    match &parameter.scope {
        ParameterScopeBinding::Cases(supported)
        | ParameterScopeBinding::MaterialLot {
            cases: supported, ..
        }
        | ParameterScopeBinding::Field {
            cases: supported, ..
        }
        | ParameterScopeBinding::Hierarchical {
            cases: supported, ..
        } => !cases.is_disjoint(supported),
        ParameterScopeBinding::Specimen { case, .. } => cases.contains(case),
        ParameterScopeBinding::Global => true,
    }
}

fn gauge_axes_match_claim(
    axes: &GaugeApplicabilityAxes,
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> bool {
    let exact_source = |key: &SourceKey, source: &SourceRef| {
        problem
            .sources
            .get(key)
            .is_some_and(|stored| stored == source)
    };
    let information = match (&axes.information, &claim.information) {
        (GaugeInformationRegime::StructuralExactModel, InformationRegime::StructuralExactModel)
        | (GaugeInformationRegime::ExactInputOutputMap, InformationRegime::ExactInputOutputMap)
        | (GaugeInformationRegime::NoisyFiniteData, InformationRegime::NoisyFiniteData) => true,
        (
            GaugeInformationRegime::PosteriorUnderDeclaredPrior { joint_prior },
            InformationRegime::PosteriorUnderDeclaredPrior {
                joint_prior: claim_prior,
            },
        ) => exact_source(joint_prior, claim_prior),
        _ => false,
    };
    let scalar = match (&axes.scalar_domain, &claim.scalar_domain) {
        (GaugeScalarDomain::Real, ScalarDomain::Real) => true,
        (
            GaugeScalarDomain::Complex { extension },
            ScalarDomain::Complex {
                extension: claim_extension,
            },
        ) => exact_source(extension, claim_extension),
        (
            GaugeScalarDomain::MixedDiscreteContinuous { stratification },
            ScalarDomain::MixedDiscreteContinuous {
                stratification: claim_stratification,
            },
        ) => exact_source(stratification, claim_stratification),
        _ => false,
    };
    let locus = match (&axes.locus, &claim.scope) {
        (GaugeLocus::WholeDomain, ClaimScope::WholeCampaign | ClaimScope::Cases(_)) => true,
        (
            GaugeLocus::Stratum { definition },
            ClaimScope::Stratum {
                definition: claim_definition,
                ..
            },
        ) => exact_source(definition, claim_definition),
        _ => false,
    };
    let quantifier = match (&axes.quantifier, &claim.quantifier) {
        (
            GaugeQuantifierScope::AtRealization { realization },
            ClaimQuantifier::AtRealization {
                realization: claim_realization,
            },
        ) => exact_source(realization, claim_realization),
        (
            GaugeQuantifierScope::AlmostEverywhere { measure },
            ClaimQuantifier::AlmostEverywhere {
                measure: claim_measure,
            },
        ) => exact_source(measure, claim_measure),
        (
            GaugeQuantifierScope::ForAll { domain },
            ClaimQuantifier::ForAll {
                domain: claim_domain,
            },
        ) => exact_source(domain, claim_domain),
        (
            GaugeQuantifierScope::ProbabilityAtLeast {
                probability,
                measure,
            },
            ClaimQuantifier::ProbabilityAtLeast {
                probability: claim_probability,
                measure: claim_measure,
            },
        ) => {
            GaugeProbabilityThreshold::try_new(*claim_probability)
                .is_ok_and(|claim_threshold| claim_threshold == *probability)
                && exact_source(measure, claim_measure)
        }
        _ => false,
    };
    information && scalar && locus && quantifier
}

fn gauge_cell_for_claim<'a>(
    validity: &'a GaugeValidityScope,
    claim: &TypedIdentifiabilityClaim,
    cases: &BTreeSet<CaseId>,
    problem: &IdentifiabilityProblemDocument,
) -> Option<&'a GaugeCellDomain> {
    validity.cells.iter().find_map(|(axes, cell)| {
        (gauge_axes_match_claim(axes, claim, problem)
            && cases.is_subset(
                &cell
                    .case_obstruction_support
                    .keys()
                    .cloned()
                    .collect::<BTreeSet<_>>(),
            ))
        .then_some(cell)
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EffectiveGaugeOrbitKind {
    Trivial,
    Continuous,
    Discrete,
    Mixed,
    Stratified,
}

fn effective_gauge_orbit_kind(geometry: &GaugeOrbitGeometry) -> EffectiveGaugeOrbitKind {
    let GaugeOrbitGeometry::Regular { principal, .. } = geometry else {
        return EffectiveGaugeOrbitKind::Stratified;
    };
    let continuous = !principal.continuous_orbit_dimension.is_zero();
    let discrete = !matches!(
        &principal.discrete_orbit_cardinality,
        GaugeDiscreteOrbitCardinality::Finite { cardinality: 1 }
    );
    match (continuous, discrete) {
        (false, false) => EffectiveGaugeOrbitKind::Trivial,
        (true, false) => EffectiveGaugeOrbitKind::Continuous,
        (false, true) => EffectiveGaugeOrbitKind::Discrete,
        (true, true) => EffectiveGaugeOrbitKind::Mixed,
    }
}

struct GaugeActionView {
    carrier: BTreeSet<ParameterRoleId>,
    obstruction: BTreeSet<ParameterRoleId>,
    orbit_kind: EffectiveGaugeOrbitKind,
}

fn claim_extent_obstruction(
    cell: &GaugeCellDomain,
    cases: &BTreeSet<CaseId>,
    extent: IdentifiabilityExtent,
) -> BTreeSet<ParameterRoleId> {
    cases
        .iter()
        .flat_map(|case| {
            let support = &cell.case_obstruction_support[case];
            match extent {
                IdentifiabilityExtent::Local => &support.local_obstruction_parameters,
                IdentifiabilityExtent::Global => &support.global_obstruction_parameters,
            }
            .iter()
            .cloned()
        })
        .collect()
}

fn gauge_action_view(
    action: &GaugeActionReference,
    claim: &TypedIdentifiabilityClaim,
    cases: &BTreeSet<CaseId>,
    problem: &IdentifiabilityProblemDocument,
) -> Result<GaugeActionView, IdentifiabilityError> {
    match action {
        GaugeActionReference::Single(id) => {
            let gauge =
                problem
                    .gauges
                    .get(id)
                    .ok_or_else(|| IdentifiabilityError::UnknownReference {
                        field: "claim fiber gauge action",
                        id: id.to_string(),
                    })?;
            let cell = gauge_cell_for_claim(&gauge.validity, claim, cases, problem)
                .ok_or_else(|| IdentifiabilityError::InvalidGauge {
                    gauge: id.clone(),
                    detail: "gauge action lacks one exact full-axis applicability cell covering the claim cases"
                        .to_string(),
                })?;
            Ok(GaugeActionView {
                carrier: gauge.members.clone(),
                obstruction: claim_extent_obstruction(cell, cases, claim.extent),
                orbit_kind: effective_gauge_orbit_kind(&gauge.orbit_geometry),
            })
        }
        GaugeActionReference::Product(id) | GaugeActionReference::Composition(id) => {
            let composition = problem.gauge_compositions.get(id).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "claim fiber gauge composition",
                    id: id.to_string(),
                }
            })?;
            let kind_matches = matches!(
                (action, &composition.kind),
                (
                    GaugeActionReference::Product(_),
                    GaugeCompositionKind::IndependentProduct
                ) | (
                    GaugeActionReference::Composition(_),
                    GaugeCompositionKind::Generated
                )
            );
            let cell = gauge_cell_for_claim(&composition.validity, claim, cases, problem);
            if !kind_matches || cell.is_none() {
                return Err(IdentifiabilityError::InvalidText {
                    field: "claim fiber gauge composition",
                    detail: "action-reference kind and exact full-axis composition applicability must both match"
                        .to_string(),
                });
            }
            let carrier = composition
                .members
                .iter()
                .flat_map(|member| problem.gauges[member].members.iter().cloned())
                .collect::<BTreeSet<_>>();
            if carrier.is_empty() {
                return Err(IdentifiabilityError::Cardinality {
                    field: "claim fiber gauge composition carrier",
                    detail: "composition carrier cannot be empty".to_string(),
                });
            }
            Ok(GaugeActionView {
                carrier,
                obstruction: claim_extent_obstruction(
                    cell.expect("checked exact composition cell"),
                    cases,
                    claim.extent,
                ),
                orbit_kind: effective_gauge_orbit_kind(&composition.effective_orbit_geometry),
            })
        }
    }
}

fn gauge_action_orbit_profile<'a>(
    action: &GaugeActionReference,
    problem: &'a IdentifiabilityProblemDocument,
) -> Option<&'a SourceKey> {
    let geometry = match action {
        GaugeActionReference::Single(id) => &problem.gauges.get(id)?.orbit_geometry,
        GaugeActionReference::Product(id) | GaugeActionReference::Composition(id) => {
            &problem.gauge_compositions.get(id)?.effective_orbit_geometry
        }
    };
    match geometry {
        GaugeOrbitGeometry::Stratified {
            orbit_type_stabilizer_profile,
            ..
        } => Some(orbit_type_stabilizer_profile),
        GaugeOrbitGeometry::Regular { .. } => None,
    }
}

fn gauge_action_status_and_geometry<'a>(
    action: &GaugeActionReference,
    problem: &'a IdentifiabilityProblemDocument,
) -> Result<(&'a GaugeStatus, &'a GaugeOrbitGeometry), IdentifiabilityError> {
    match action {
        GaugeActionReference::Single(id) => {
            let gauge =
                problem
                    .gauges
                    .get(id)
                    .ok_or_else(|| IdentifiabilityError::UnknownReference {
                        field: "gauge action",
                        id: id.to_string(),
                    })?;
            Ok((&gauge.status, &gauge.orbit_geometry))
        }
        GaugeActionReference::Product(id) | GaugeActionReference::Composition(id) => {
            let composition = problem.gauge_compositions.get(id).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "gauge composition action",
                    id: id.to_string(),
                }
            })?;
            Ok((&composition.status, &composition.effective_orbit_geometry))
        }
    }
}

fn action_supports_claim_subject(
    action: &GaugeActionReference,
    view: &GaugeActionView,
    claim: &TypedIdentifiabilityClaim,
    parameters: &BTreeSet<ParameterRoleId>,
) -> bool {
    match &claim.subject {
        ClaimSubject::Parameter(_) | ClaimSubject::ParameterSet(_) => {
            !view.obstruction.is_empty() && view.obstruction.is_subset(parameters)
        }
        // An induced action or descent on a composite observable is a theorem
        // obligation of positive assessment, not a constructor inference.
        ClaimSubject::DerivedFunctional { .. } | ClaimSubject::Influence(_) => {
            !parameters.is_disjoint(&view.obstruction)
        }
        ClaimSubject::GaugeAction(subject_action) => {
            subject_action == action && !view.obstruction.is_empty()
        }
        ClaimSubject::WholeProblem => {
            !view.obstruction.is_empty() && view.obstruction.is_subset(parameters)
        }
    }
}

fn validate_claim_compatibility(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<(), IdentifiabilityError> {
    let cases = claim_case_set(claim, problem)?;
    let parameters = claim_subject_parameters(claim, problem)?;
    for role in &parameters {
        if !parameter_active_in_cases(&problem.parameters[role], &cases) {
            return Err(IdentifiabilityError::InvalidText {
                field: "claim parameter/case scope",
                detail: format!("parameter {role} is inactive throughout the claimed case scope"),
            });
        }
    }

    if let ClaimSubject::Influence(influence_id) = &claim.subject {
        let support = transitive_influence_ids(influence_id, &problem.influences)?
            .into_iter()
            .flat_map(|id| {
                functional_observations(&problem.influences[&id].functional)
                    .into_iter()
                    .map(|observation| observation.case.clone())
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>();
        if !support.is_subset(&cases) {
            return Err(IdentifiabilityError::InvalidText {
                field: "claim influence/case scope",
                detail: format!(
                    "influence {influence_id} has an observation endpoint outside the claimed scope"
                ),
            });
        }
    }

    match &claim.information {
        InformationRegime::NoisyFiniteData => {
            if cases
                .iter()
                .any(|case| matches!(&problem.cases[case].data, CaseDataDeclaration::Prospective))
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "finite-data claim scope",
                    detail:
                        "NoisyFiniteData claims require retrospective data in every claimed case"
                            .to_string(),
                });
            }
        }
        InformationRegime::PosteriorUnderDeclaredPrior { .. } => {
            if cases
                .iter()
                .any(|case| matches!(&problem.cases[case].data, CaseDataDeclaration::Prospective))
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "posterior claim scope",
                    detail: "PosteriorUnderDeclaredPrior denotes an observed posterior and requires retrospective admitted data in every claimed case; planned prior-predictive questions need a distinct non-decisive regime"
                        .to_string(),
                });
            }
            let free = problem
                .parameters
                .iter()
                .filter(|(_, parameter)| {
                    is_free_inferential_parameter(parameter)
                        && parameter_active_in_cases(parameter, &cases)
                })
                .map(|(role, _)| role)
                .collect::<Vec<_>>();
            if free.is_empty()
                || free.iter().any(|role| {
                    !matches!(
                        &problem.parameters[*role].prior,
                        PriorPolicy::Distribution(_)
                    )
                })
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "posterior claim prior",
                    detail:
                        "every free parameter in a posterior claim needs an explicit declared prior"
                            .to_string(),
                });
            }
        }
        InformationRegime::StructuralExactModel | InformationRegime::ExactInputOutputMap => {}
    }

    match &claim.fiber {
        // Requests are propositions, not theorem assertions. Even an unlikely
        // raw Unique request remains preregisterable; only a positive
        // assessment must dispose every exact applicable obstruction.
        FiberStructure::Unique => {}
        FiberStructure::FiniteToOne {
            maximum_cardinality,
        } => {
            if matches!(
                maximum_cardinality,
                Some(FiberCardinalityBound::UniformU64(0 | 1))
            ) {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "finite-to-one cardinality",
                    detail: "a known maximum below two is canonically represented as Unique"
                        .to_string(),
                });
            }
        }
        FiberStructure::PositiveDimensional {
            lower_bound:
                FiberDimensionLowerBound::Finite {
                    minimum_dimension: 0,
                },
        } => {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "positive-dimensional fiber",
                detail: "minimum dimension must be positive".to_string(),
            });
        }
        FiberStructure::PositiveDimensional { .. } => {}
        FiberStructure::DiscreteOrbit { action } => {
            let view = gauge_action_view(action, claim, &cases, problem)?;
            if !action_supports_claim_subject(action, &view, claim, &parameters)
                || view.orbit_kind != EffectiveGaugeOrbitKind::Discrete
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "discrete-orbit fiber",
                    detail: "claim support must lie in the exact action carrier and the effective regular orbit—not merely the acting group—must be purely discrete"
                        .to_string(),
                });
            }
        }
        FiberStructure::MixedOrbit { action } => {
            let view = gauge_action_view(action, claim, &cases, problem)?;
            if !action_supports_claim_subject(action, &view, claim, &parameters)
                || view.orbit_kind != EffectiveGaugeOrbitKind::Mixed
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "mixed-orbit fiber",
                    detail: "claim support must lie in the exact action carrier and the effective regular orbit must have both continuous and discrete components"
                        .to_string(),
                });
            }
        }
        FiberStructure::OrbitQuotientUnique { action } => {
            let view = gauge_action_view(action, claim, &cases, problem)?;
            if !action_supports_claim_subject(action, &view, claim, &parameters)
                || matches!(view.orbit_kind, EffectiveGaugeOrbitKind::Trivial)
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "orbit-quotient fiber",
                    detail: "claim support must lie in a nontrivial exact action carrier"
                        .to_string(),
                });
            }
        }
        FiberStructure::GeneralizedQuotientUnique {
            action,
            equivalence,
        } => {
            let view = gauge_action_view(action, claim, &cases, problem)?;
            if !action_supports_claim_subject(action, &view, claim, &parameters)
                || equivalence.kind != SourceKind::GaugeQuotientProfile
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "generalized-quotient fiber",
                    detail: "generalized quotient needs an exact action carrier and exact orbit-closure/invariant/groupoid/stack equivalence profile"
                        .to_string(),
                });
            }
        }
        FiberStructure::StratifiedOrbit {
            action,
            orbit_type_profile,
        } => {
            let view = gauge_action_view(action, claim, &cases, problem)?;
            let profile_matches = gauge_action_orbit_profile(action, problem)
                .is_some_and(|key| problem.sources.get(key) == Some(orbit_type_profile));
            if !action_supports_claim_subject(action, &view, claim, &parameters)
                || view.orbit_kind != EffectiveGaugeOrbitKind::Stratified
                || !profile_matches
            {
                return Err(IdentifiabilityError::InvalidText {
                    field: "stratified-orbit fiber",
                    detail: "stratified orbit claims require the exact action, exact applicability cell, and exact problem-bound orbit-type profile"
                        .to_string(),
                });
            }
        }
        FiberStructure::Stratified { strata } => validate_claim_source_kind(
            strata,
            &[SourceKind::Stratification],
            "claim fiber stratification",
        )?,
    }
    Ok(())
}

struct RelevantGaugeObstruction {
    action: GaugeActionReference,
    members: BTreeSet<GaugeClassId>,
    assumed: bool,
    covered_cases: BTreeSet<CaseId>,
    obstruction_support: BTreeSet<ParameterRoleId>,
    fully_covers_claim_cases: bool,
    orbit_kind: EffectiveGaugeOrbitKind,
    orbit_geometry: GaugeOrbitGeometry,
}

fn matching_gauge_cell<'a>(
    validity: &'a GaugeValidityScope,
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Option<&'a GaugeCellDomain> {
    validity
        .cells
        .iter()
        .find_map(|(axes, cell)| gauge_axes_match_claim(axes, claim, problem).then_some(cell))
}

fn relevant_gauge_obstructions(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<Vec<RelevantGaugeObstruction>, IdentifiabilityError> {
    let cases = claim_case_set(claim, problem)?;
    let subject = claim_subject_parameters(claim, problem)?;
    let mut relevant = Vec::new();
    let obstruction_for = |cell: &GaugeCellDomain| {
        let covered = cases
            .iter()
            .filter(|case| cell.case_obstruction_support.contains_key(*case))
            .cloned()
            .collect::<BTreeSet<_>>();
        let full = covered.len() == cases.len();
        let obstruction = claim_extent_obstruction(cell, &covered, claim.extent);
        (obstruction, covered, full)
    };
    for gauge in problem.gauges.values() {
        let Some(cell) = matching_gauge_cell(&gauge.validity, claim, problem) else {
            continue;
        };
        let (obstruction, covered, full) = obstruction_for(cell);
        if obstruction.is_disjoint(&subject) {
            continue;
        }
        relevant.push(RelevantGaugeObstruction {
            action: GaugeActionReference::Single(gauge.id.clone()),
            members: BTreeSet::from([gauge.id.clone()]),
            assumed: matches!(&gauge.status, GaugeStatus::Assumed { .. }),
            covered_cases: covered,
            obstruction_support: obstruction,
            fully_covers_claim_cases: full,
            orbit_kind: effective_gauge_orbit_kind(&gauge.orbit_geometry),
            orbit_geometry: gauge.orbit_geometry.clone(),
        });
    }
    for composition in problem.gauge_compositions.values() {
        let Some(cell) = matching_gauge_cell(&composition.validity, claim, problem) else {
            continue;
        };
        let (obstruction, covered, full) = obstruction_for(cell);
        if obstruction.is_disjoint(&subject) {
            continue;
        }
        let action = match &composition.kind {
            GaugeCompositionKind::IndependentProduct => {
                GaugeActionReference::Product(composition.id.clone())
            }
            GaugeCompositionKind::Generated => {
                GaugeActionReference::Composition(composition.id.clone())
            }
        };
        relevant.push(RelevantGaugeObstruction {
            action,
            members: composition.members.clone(),
            assumed: matches!(&composition.status, GaugeStatus::Assumed { .. })
                && composition.members.iter().all(|member| {
                    matches!(&problem.gauges[member].status, GaugeStatus::Assumed { .. })
                }),
            covered_cases: covered,
            obstruction_support: obstruction,
            fully_covers_claim_cases: full,
            orbit_kind: effective_gauge_orbit_kind(&composition.effective_orbit_geometry),
            orbit_geometry: composition.effective_orbit_geometry.clone(),
        });
    }
    let maximal = relevant
        .into_iter()
        .filter(|candidate| {
            !problem.gauge_compositions.values().any(|composition| {
                let composed_action = match &composition.kind {
                    GaugeCompositionKind::IndependentProduct => {
                        GaugeActionReference::Product(composition.id.clone())
                    }
                    GaugeCompositionKind::Generated => {
                        GaugeActionReference::Composition(composition.id.clone())
                    }
                };
                let effectively_assumed =
                    matches!(&composition.status, GaugeStatus::Assumed { .. })
                        && composition.members.iter().all(|member| {
                            matches!(&problem.gauges[member].status, GaugeStatus::Assumed { .. })
                        });
                candidate.members.is_subset(&composition.members)
                    && candidate.members.len() < composition.members.len()
                    && (effectively_assumed || claim_gauge_action(claim) == Some(&composed_action))
                    && matching_gauge_cell(&composition.validity, claim, problem).is_some_and(
                        |cell| {
                            let (obstruction, covered, _) = obstruction_for(cell);
                            candidate.covered_cases.is_subset(&covered)
                                && candidate.obstruction_support.is_subset(&obstruction)
                        },
                    )
            })
        })
        .collect();
    Ok(maximal)
}

fn disposition_is_valid(
    claim: &TypedIdentifiabilityClaim,
    obstruction: &RelevantGaugeObstruction,
    disposition: &GaugeResolutionDisposition,
) -> bool {
    match disposition {
        GaugeResolutionDisposition::CandidateRefuted => !obstruction.assumed,
        GaugeResolutionDisposition::NoProjectionOnSubject
        | GaugeResolutionDisposition::SubjectDescendsToQuotient
        | GaugeResolutionDisposition::TrivialResidualIntersection => matches!(
            &claim.subject,
            ClaimSubject::DerivedFunctional { .. } | ClaimSubject::Influence(_)
        ),
        GaugeResolutionDisposition::BrokenByJointInformation => {
            !obstruction.fully_covers_claim_cases
        }
        GaugeResolutionDisposition::ConsistentWithClaimedFiber => {
            let exact_action = claim_gauge_action(claim) == Some(&obstruction.action);
            let action_consistent = exact_action
                && match (&claim.fiber, obstruction.orbit_kind) {
                    (FiberStructure::DiscreteOrbit { .. }, EffectiveGaugeOrbitKind::Discrete)
                    | (FiberStructure::MixedOrbit { .. }, EffectiveGaugeOrbitKind::Mixed)
                    | (
                        FiberStructure::StratifiedOrbit { .. },
                        EffectiveGaugeOrbitKind::Stratified,
                    ) => true,
                    (
                        FiberStructure::OrbitQuotientUnique { .. }
                        | FiberStructure::GeneralizedQuotientUnique { .. },
                        EffectiveGaugeOrbitKind::Continuous
                        | EffectiveGaugeOrbitKind::Discrete
                        | EffectiveGaugeOrbitKind::Mixed
                        | EffectiveGaugeOrbitKind::Stratified,
                    ) => true,
                    _ => false,
                };
            if action_consistent {
                return true;
            }
            let GaugeOrbitGeometry::Regular { principal, .. } = &obstruction.orbit_geometry else {
                return matches!(&claim.fiber, FiberStructure::Stratified { .. });
            };
            match &claim.fiber {
                FiberStructure::FiniteToOne {
                    maximum_cardinality,
                } if principal.continuous_orbit_dimension.is_zero() => {
                    match (&principal.discrete_orbit_cardinality, maximum_cardinality) {
                        (
                            GaugeDiscreteOrbitCardinality::Finite { cardinality },
                            Some(FiberCardinalityBound::UniformU64(maximum)),
                        ) => cardinality <= maximum,
                        (GaugeDiscreteOrbitCardinality::Finite { .. }, None) => true,
                        (
                            GaugeDiscreteOrbitCardinality::Finite { .. },
                            Some(FiberCardinalityBound::SymbolicProfile(_)),
                        ) => false,
                        (GaugeDiscreteOrbitCardinality::CountablyInfinite { .. }, _) => false,
                    }
                }
                FiberStructure::PositiveDimensional { lower_bound } => {
                    match (&principal.continuous_orbit_dimension, lower_bound) {
                        (
                            GaugeContinuousDimension::Finite { dimension },
                            FiberDimensionLowerBound::Finite { minimum_dimension },
                        ) => dimension >= minimum_dimension && *dimension > 0,
                        (
                            GaugeContinuousDimension::InfiniteDimensional { .. },
                            FiberDimensionLowerBound::Finite { .. },
                        ) => true,
                        (
                            GaugeContinuousDimension::InfiniteDimensional { .. },
                            FiberDimensionLowerBound::InfiniteDimensional { .. },
                        ) => false,
                        (
                            GaugeContinuousDimension::Finite { .. },
                            FiberDimensionLowerBound::InfiniteDimensional { .. },
                        ) => false,
                    }
                }
                FiberStructure::Stratified { .. } => {
                    obstruction.orbit_kind == EffectiveGaugeOrbitKind::Stratified
                }
                _ => false,
            }
        }
    }
}

fn validate_positive_claim_influence_routes(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<(), IdentifiabilityError> {
    if matches!(&claim.subject, ClaimSubject::Influence(_)) {
        return Ok(());
    }
    let parameters = claim_subject_parameters(claim, problem)?;
    let zero_route = parameters.iter().find(|role| {
        matches!(
            &problem.parameters[*role].influence_coverage,
            InfluenceCoverage::IntentionallyAbsent { .. }
        )
    });
    if let Some(role) = zero_route
        && !matches!(&claim.fiber, FiberStructure::PositiveDimensional { .. })
    {
        return Err(IdentifiabilityError::InvalidText {
            field: "positive claim influence coverage",
            detail: format!(
                "positive claim {} includes free parameter {role} whose schema declares no influence route; only an explicit positive-dimensional unidentifiability proposition may be established until a typed constraint-mediated influence route exists",
                claim.id
            ),
        });
    }
    Ok(())
}

fn validate_decisive_practical_claim_closure(
    claim: &TypedIdentifiabilityClaim,
    problem: &IdentifiabilityProblemDocument,
) -> Result<(), IdentifiabilityError> {
    if !matches!(
        &claim.information,
        InformationRegime::NoisyFiniteData | InformationRegime::PosteriorUnderDeclaredPrior { .. }
    ) {
        return Ok(());
    }
    if matches!(&problem.joint_noise, JointNoiseModel::Unknown { .. }) {
        return Err(IdentifiabilityError::InvalidText {
            field: "positive practical claim joint noise",
            detail: "a positive noisy-data/posterior claim cannot close while the joint noise model is Unknown"
                .to_string(),
        });
    }
    for case_id in claim_case_set(claim, problem)? {
        let case = &problem.cases[&case_id];
        for (channel, observation) in &case.observations {
            if matches!(&observation.noise, MarginalNoiseSpec::Unknown { .. }) {
                return Err(IdentifiabilityError::InvalidText {
                    field: "positive practical claim marginal noise",
                    detail: format!(
                        "case {case_id} channel {channel} retains Unknown marginal noise"
                    ),
                });
            }
            if matches!(
                &observation.missingness,
                MissingnessAssumption::Unknown { .. }
            ) {
                return Err(IdentifiabilityError::InvalidText {
                    field: "positive practical claim missingness",
                    detail: format!("case {case_id} channel {channel} retains Unknown missingness"),
                });
            }
            if matches!(
                &case.discrepancies[channel],
                StudyDiscrepancy::Uncharacterized { .. }
            ) {
                return Err(IdentifiabilityError::InvalidText {
                    field: "positive practical claim discrepancy",
                    detail: format!(
                        "case {case_id} channel {channel} retains Uncharacterized model discrepancy"
                    ),
                });
            }
        }
    }
    Ok(())
}

impl IdentifiabilityAssessment {
    pub fn try_new(
        header: ArtifactHeader,
        problem: &AdmittedIdentifiabilityProblem,
        execution: &IdentifiabilityExecutionPlan,
        claims: Vec<TypedIdentifiabilityClaim>,
        evidence: Vec<(ClaimId, ClaimAssessment)>,
        source_authority: SourceResolutionSet,
    ) -> Result<Self, IdentifiabilityError> {
        validate_header_profile(&header)?;
        if !header.capabilities().contains("identifiability.assess") {
            return Err(IdentifiabilityError::InvalidText {
                field: "assessment capability",
                detail: "missing identifiability.assess capability".to_string(),
            });
        }
        if execution.problem_id != problem.problem_id
            || execution.source_admission_id != problem.source_admission_id
        {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "assessment problem/execution",
            });
        }
        if evidence.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "claim assessments",
                detail: "assessment evidence input exceeds the canonical collection bound"
                    .to_string(),
            });
        }
        if evidence.iter().any(|(_, conclusion)| {
            matches!(
                conclusion,
                ClaimAssessment::ClaimedEstablished {
                    gauge_resolutions,
                    ..
                } if gauge_resolutions.len() > MAX_IDENTIFIABILITY_ITEMS
            )
        }) {
            return Err(IdentifiabilityError::Cardinality {
                field: "positive-claim gauge resolutions",
                detail: "gauge-resolution evidence exceeds the canonical collection bound"
                    .to_string(),
            });
        }
        let claims = insert_unique(claims, "identifiability claims", |claim| &claim.id)?;
        let mut evidence_map = BTreeMap::new();
        for (id, mut conclusion) in evidence {
            match &mut conclusion {
                ClaimAssessment::ClaimedEstablished {
                    certified_error_bound,
                    ..
                }
                | ClaimAssessment::ClaimedRefuted {
                    certified_error_bound,
                    ..
                } => {
                    *certified_error_bound = canonical_f64(*certified_error_bound);
                }
                ClaimAssessment::ClaimedInconclusive { .. }
                | ClaimAssessment::NotAssessed { .. } => {}
            }
            if evidence_map.insert(id.clone(), conclusion).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "claim assessment",
                    id: id.to_string(),
                });
            }
        }
        let claim_ids = claims.keys().cloned().collect::<BTreeSet<_>>();
        let evidence_ids = evidence_map.keys().cloned().collect::<BTreeSet<_>>();
        if evidence_ids != claim_ids {
            return Err(IdentifiabilityError::Cardinality {
                field: "claim assessments",
                detail: "claim and assessment identity sets must match exactly".to_string(),
            });
        }
        let preregistered_claims = execution
            .claim_requests
            .iter()
            .map(|(id, request)| (id.clone(), request.claim.clone()))
            .collect::<BTreeMap<_, _>>();
        if claims != preregistered_claims {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "assessment/execution exact claim preregistration",
            });
        }
        let mut referenced_sources = BTreeMap::<SourceKey, SourceRef>::new();
        let mut bind_source = |source: &SourceRef| -> Result<(), IdentifiabilityError> {
            if let Some(prior) = referenced_sources.insert(source.key.clone(), source.clone()) {
                if &prior != source {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "assessment source alias",
                    });
                }
            }
            Ok(())
        };
        for (id, claim) in &claims {
            validate_claim_compatibility(claim, &problem.document)?;
            for source in validate_claim_sources(claim, &problem.document)? {
                bind_source(source)?;
            }
            let request = &execution.claim_requests[id];
            match &evidence_map[id] {
                ClaimAssessment::ClaimedEstablished {
                    method,
                    receipt,
                    metric,
                    nondimensionalization,
                    certified_error_bound,
                    gauge_resolutions,
                } => {
                    validate_positive_claim_influence_routes(claim, &problem.document)?;
                    validate_decisive_practical_claim_closure(claim, &problem.document)?;
                    if method.kind != SourceKind::Analyzer
                        || receipt.kind != SourceKind::EvidenceReceipt
                        || !certified_error_bound.is_finite()
                        || *certified_error_bound < 0.0
                    {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "claim evidence",
                            detail: "established claims need analyzer, receipt, and finite nonnegative certified error"
                                .to_string(),
                        });
                    }
                    if method != &execution.analyzer {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "claim analyzer/execution analyzer",
                        });
                    }
                    if metric != request.error_policy.metric()
                        || nondimensionalization != request.error_policy.nondimensionalization()
                    {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "claim evidence/error policy",
                        });
                    }
                    if *certified_error_bound > request.error_policy.maximum_certified_error() {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "certified claim error",
                            detail: format!(
                                "certified error {certified_error_bound:e} exceeds preregistered maximum {:e}",
                                request.error_policy.maximum_certified_error()
                            ),
                        });
                    }
                    let relevant = relevant_gauge_obstructions(claim, &problem.document)?;
                    let expected = relevant
                        .iter()
                        .map(|obstruction| obstruction.action.clone())
                        .collect::<BTreeSet<_>>();
                    let actual = gauge_resolutions.keys().cloned().collect::<BTreeSet<_>>();
                    if actual != expected
                        || gauge_resolutions
                            .iter()
                            .any(|(action, resolution)| action != &resolution.action)
                    {
                        return Err(IdentifiabilityError::Cardinality {
                            field: "positive-claim gauge resolutions",
                            detail: "positive assessment must resolve every and only exact applicable single/product/composition obstruction once"
                                .to_string(),
                        });
                    }
                    for resolution in gauge_resolutions.values() {
                        if resolution.method.kind != SourceKind::Analyzer
                            || &resolution.method != execution.analyzer()
                            || resolution.receipt.kind != SourceKind::EvidenceReceipt
                        {
                            return Err(IdentifiabilityError::SourceMismatch {
                                field: "gauge resolution evidence/execution analyzer",
                            });
                        }
                        let obstruction = relevant
                            .iter()
                            .find(|obstruction| obstruction.action == resolution.action)
                            .expect("exact gauge resolution set checked");
                        if !disposition_is_valid(claim, obstruction, &resolution.disposition) {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "gauge resolution disposition",
                                detail: format!(
                                    "disposition {:?} is incompatible with the target status, exact claim cell, subject projection, or claimed fiber",
                                    resolution.disposition
                                ),
                            });
                        }
                        bind_source(&resolution.method)?;
                        bind_source(&resolution.receipt)?;
                    }
                    bind_source(method)?;
                    bind_source(receipt)?;
                    bind_source(metric)?;
                    bind_source(nondimensionalization)?;
                }
                ClaimAssessment::ClaimedRefuted {
                    method,
                    receipt,
                    metric,
                    nondimensionalization,
                    certified_error_bound,
                } => {
                    validate_decisive_practical_claim_closure(claim, &problem.document)?;
                    if method.kind != SourceKind::Analyzer
                        || receipt.kind != SourceKind::EvidenceReceipt
                        || !certified_error_bound.is_finite()
                        || *certified_error_bound < 0.0
                    {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "claim evidence",
                            detail: "refuted claims need analyzer, receipt, and finite nonnegative certified error"
                            .to_string(),
                        });
                    }
                    if method != &execution.analyzer {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "claim analyzer/execution analyzer",
                        });
                    }
                    if metric != request.error_policy.metric()
                        || nondimensionalization != request.error_policy.nondimensionalization()
                    {
                        return Err(IdentifiabilityError::SourceMismatch {
                            field: "claim evidence/error policy",
                        });
                    }
                    if *certified_error_bound > request.error_policy.maximum_certified_error() {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "certified claim error",
                            detail: format!(
                                "certified error {certified_error_bound:e} exceeds preregistered maximum {:e}",
                                request.error_policy.maximum_certified_error()
                            ),
                        });
                    }
                    bind_source(method)?;
                    bind_source(receipt)?;
                    bind_source(metric)?;
                    bind_source(nondimensionalization)?;
                }
                ClaimAssessment::ClaimedInconclusive {
                    method,
                    receipt,
                    reason,
                } => {
                    validate_reason(reason, "inconclusive claim reason")?;
                    if method.is_some() != receipt.is_some() {
                        return Err(IdentifiabilityError::InvalidText {
                            field: "inconclusive evidence pair",
                            detail: "inconclusive evidence must bind method and receipt together, or omit both"
                                .to_string(),
                        });
                    }
                    if let Some(method) = method {
                        if method.kind != SourceKind::Analyzer {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "inconclusive method",
                                detail: "method source must have Analyzer kind".to_string(),
                            });
                        }
                        if method != &execution.analyzer {
                            return Err(IdentifiabilityError::SourceMismatch {
                                field: "inconclusive analyzer/execution analyzer",
                            });
                        }
                        bind_source(method)?;
                    }
                    if let Some(receipt) = receipt {
                        if receipt.kind != SourceKind::EvidenceReceipt {
                            return Err(IdentifiabilityError::InvalidText {
                                field: "inconclusive receipt",
                                detail: "receipt source must have EvidenceReceipt kind".to_string(),
                            });
                        }
                        bind_source(receipt)?;
                    }
                }
                ClaimAssessment::NotAssessed { reason } => {
                    validate_reason(reason, "not-assessed claim reason")?
                }
            }
        }
        if source_authority.entries.len() != referenced_sources.len() {
            return Err(IdentifiabilityError::Cardinality {
                field: "assessment source authority",
                detail: "every claim/method/receipt source needs exactly one resolution"
                    .to_string(),
            });
        }
        for (key, reference) in &referenced_sources {
            let resolution = source_authority.entries.get(key).ok_or_else(|| {
                IdentifiabilityError::UnknownReference {
                    field: "assessment source authority",
                    id: key.to_string(),
                }
            })?;
            admit_opaque_resolution(reference, resolution)?;
            if let Some(execution_resolution) = execution.source_authority.entries.get(key)
                && execution_resolution != resolution
            {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "assessment/execution source authority",
                });
            }
            if let Some(problem_resolution) = problem.source_resolutions().get(key)
                && problem_resolution != resolution
            {
                return Err(IdentifiabilityError::SourceMismatch {
                    field: "assessment/problem source authority",
                });
            }
        }
        let assessment = Self {
            schema_version: IDENTIFIABILITY_ASSESSMENT_IDENTITY_VERSION,
            header,
            problem_id: problem.problem_id,
            execution_id: execution.id()?,
            claims,
            evidence: evidence_map,
            source_authority,
        };
        validate_assessment_structural_budget(&assessment)?;
        let _ = encode_assessment(&assessment)?;
        Ok(assessment)
    }

    pub fn canonical_bytes(&self) -> Result<Vec<u8>, IdentifiabilityError> {
        encode_assessment(self)
    }

    pub fn id(&self) -> Result<AssessmentId, IdentifiabilityError> {
        assessment_identity_hash(self)
    }

    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn problem_id(&self) -> ProblemId {
        self.problem_id
    }

    #[must_use]
    pub const fn execution_id(&self) -> ExecutionId {
        self.execution_id
    }

    #[must_use]
    pub const fn claims(&self) -> &BTreeMap<ClaimId, TypedIdentifiabilityClaim> {
        &self.claims
    }

    #[must_use]
    pub const fn evidence(&self) -> &BTreeMap<ClaimId, ClaimAssessment> {
        &self.evidence
    }

    /// Locally verified source resolutions required to replay this assessment.
    /// Serialized verification markers are never accepted as a substitute.
    #[must_use]
    pub const fn source_authority(&self) -> &SourceResolutionSet {
        &self.source_authority
    }

    pub fn from_canonical_bytes(
        bytes: &[u8],
        problem: &AdmittedIdentifiabilityProblem,
        execution: &IdentifiabilityExecutionPlan,
        verified_sources: &SourceResolutionSet,
    ) -> Result<Self, IdentifiabilityError> {
        decode_assessment(bytes, problem, execution, verified_sources)
    }
}

fn encode_source_key(
    writer: &mut CanonicalWriter,
    key: &SourceKey,
) -> Result<(), IdentifiabilityError> {
    writer.text(key.as_str(), "source key")
}

fn decode_source_key(reader: &mut CanonicalReader<'_>) -> Result<SourceKey, IdentifiabilityError> {
    SourceKey::try_new(reader.token("source key")?)
}

fn encode_case_id(writer: &mut CanonicalWriter, id: &CaseId) -> Result<(), IdentifiabilityError> {
    writer.text(id.as_str(), "case id")
}

fn decode_case_id(reader: &mut CanonicalReader<'_>) -> Result<CaseId, IdentifiabilityError> {
    CaseId::try_new(reader.token("case id")?)
}

fn encode_role(
    writer: &mut CanonicalWriter,
    role: &ParameterRoleId,
) -> Result<(), IdentifiabilityError> {
    writer.text(role.as_str(), "parameter role")
}

fn decode_role(reader: &mut CanonicalReader<'_>) -> Result<ParameterRoleId, IdentifiabilityError> {
    ParameterRoleId::try_new(reader.token("parameter role")?)
}

fn encode_channel(
    writer: &mut CanonicalWriter,
    channel: &ObservationChannelId,
) -> Result<(), IdentifiabilityError> {
    writer.text(channel.as_str(), "observation channel")
}

fn decode_channel(
    reader: &mut CanonicalReader<'_>,
) -> Result<ObservationChannelId, IdentifiabilityError> {
    ObservationChannelId::try_new(reader.token("observation channel")?)
}

fn encode_source_kind(writer: &mut CanonicalWriter, kind: SourceKind) {
    writer.byte(match kind {
        SourceKind::ContextOfUse => 0,
        SourceKind::MaterialCard => 1,
        SourceKind::ConstitutiveModelCard => 2,
        SourceKind::ConstitutiveGraph => 3,
        SourceKind::ExperimentArtifact => 4,
        SourceKind::CalibrationSplit => 5,
        SourceKind::ForwardModel => 6,
        SourceKind::Geometry => 7,
        SourceKind::Process => 8,
        SourceKind::Protocol => 9,
        SourceKind::ObservationOperator => 10,
        SourceKind::Metrology => 11,
        SourceKind::Parser => 12,
        SourceKind::Preprocessing => 13,
        SourceKind::Likelihood => 14,
        SourceKind::Prior => 15,
        SourceKind::Constraint => 16,
        SourceKind::GaugeAction => 17,
        SourceKind::GaugeSection => 18,
        SourceKind::Discrepancy => 19,
        SourceKind::Assumption => 20,
        SourceKind::Analyzer => 21,
        SourceKind::DerivativeProvider => 22,
        SourceKind::Build => 23,
        SourceKind::EvidenceReceipt => 24,
        SourceKind::ExternalManifold => 25,
        SourceKind::AlgebraicExtension => 26,
        SourceKind::Stratification => 27,
        SourceKind::DerivedFunctional => 28,
        SourceKind::AdmissibleDomainCertificate => 29,
        SourceKind::DimensionlessErrorMetric => 30,
        SourceKind::Nondimensionalization => 31,
        SourceKind::QuantifierRealization => 32,
        SourceKind::ReferenceMeasure => 33,
        SourceKind::ProbabilityMeasure => 34,
        SourceKind::QuantifierDomain => 35,
        SourceKind::GaugeComposition => 36,
        SourceKind::GaugeOrbitTypeProfile => 37,
        SourceKind::GaugeHypothesis => 38,
        SourceKind::GaugeGroupPresentation => 39,
        SourceKind::GaugeOrbitPresentation => 40,
        SourceKind::InfluenceComposition => 41,
        SourceKind::GaugeQuotientProfile => 42,
        SourceKind::FiberCardinalityProfile => 43,
        SourceKind::FiberDimensionProfile => 44,
        SourceKind::FunctionalModelSpace => 45,
        SourceKind::GaugeQuotientMap => 46,
        SourceKind::GaugeInvariantMap => 47,
        SourceKind::GaugeGroupoidPresentation => 48,
        SourceKind::GaugeReductionLaw => 49,
        SourceKind::GaugeSubgroupCertificate => 50,
        SourceKind::GaugeResidualAction => 51,
        SourceKind::GaugeMeasureTransport => 52,
        SourceKind::MeasureTransport => 53,
        SourceKind::ParameterizedLikelihood => 54,
        SourceKind::UnitDefinition => 55,
        SourceKind::ForwardModelProductionBinding => 56,
    });
}

fn decode_source_kind(
    reader: &mut CanonicalReader<'_>,
) -> Result<SourceKind, IdentifiabilityError> {
    match reader.byte("source kind")? {
        0 => Ok(SourceKind::ContextOfUse),
        1 => Ok(SourceKind::MaterialCard),
        2 => Ok(SourceKind::ConstitutiveModelCard),
        3 => Ok(SourceKind::ConstitutiveGraph),
        4 => Ok(SourceKind::ExperimentArtifact),
        5 => Ok(SourceKind::CalibrationSplit),
        6 => Ok(SourceKind::ForwardModel),
        7 => Ok(SourceKind::Geometry),
        8 => Ok(SourceKind::Process),
        9 => Ok(SourceKind::Protocol),
        10 => Ok(SourceKind::ObservationOperator),
        11 => Ok(SourceKind::Metrology),
        12 => Ok(SourceKind::Parser),
        13 => Ok(SourceKind::Preprocessing),
        14 => Ok(SourceKind::Likelihood),
        15 => Ok(SourceKind::Prior),
        16 => Ok(SourceKind::Constraint),
        17 => Ok(SourceKind::GaugeAction),
        18 => Ok(SourceKind::GaugeSection),
        19 => Ok(SourceKind::Discrepancy),
        20 => Ok(SourceKind::Assumption),
        21 => Ok(SourceKind::Analyzer),
        22 => Ok(SourceKind::DerivativeProvider),
        23 => Ok(SourceKind::Build),
        24 => Ok(SourceKind::EvidenceReceipt),
        25 => Ok(SourceKind::ExternalManifold),
        26 => Ok(SourceKind::AlgebraicExtension),
        27 => Ok(SourceKind::Stratification),
        28 => Ok(SourceKind::DerivedFunctional),
        29 => Ok(SourceKind::AdmissibleDomainCertificate),
        30 => Ok(SourceKind::DimensionlessErrorMetric),
        31 => Ok(SourceKind::Nondimensionalization),
        32 => Ok(SourceKind::QuantifierRealization),
        33 => Ok(SourceKind::ReferenceMeasure),
        34 => Ok(SourceKind::ProbabilityMeasure),
        35 => Ok(SourceKind::QuantifierDomain),
        36 => Ok(SourceKind::GaugeComposition),
        37 => Ok(SourceKind::GaugeOrbitTypeProfile),
        38 => Ok(SourceKind::GaugeHypothesis),
        39 => Ok(SourceKind::GaugeGroupPresentation),
        40 => Ok(SourceKind::GaugeOrbitPresentation),
        41 => Ok(SourceKind::InfluenceComposition),
        42 => Ok(SourceKind::GaugeQuotientProfile),
        43 => Ok(SourceKind::FiberCardinalityProfile),
        44 => Ok(SourceKind::FiberDimensionProfile),
        45 => Ok(SourceKind::FunctionalModelSpace),
        46 => Ok(SourceKind::GaugeQuotientMap),
        47 => Ok(SourceKind::GaugeInvariantMap),
        48 => Ok(SourceKind::GaugeGroupoidPresentation),
        49 => Ok(SourceKind::GaugeReductionLaw),
        50 => Ok(SourceKind::GaugeSubgroupCertificate),
        51 => Ok(SourceKind::GaugeResidualAction),
        52 => Ok(SourceKind::GaugeMeasureTransport),
        53 => Ok(SourceKind::MeasureTransport),
        54 => Ok(SourceKind::ParameterizedLikelihood),
        55 => Ok(SourceKind::UnitDefinition),
        56 => Ok(SourceKind::ForwardModelProductionBinding),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown source kind tag {tag}"),
        }),
    }
}

fn preflight_collection_bytes(
    reader: &CanonicalReader<'_>,
    count: usize,
    minimum_item_bytes: usize,
    field: &'static str,
) -> Result<(), IdentifiabilityError> {
    let required =
        count
            .checked_mul(minimum_item_bytes)
            .ok_or_else(|| IdentifiabilityError::Cardinality {
                field,
                detail: "minimum encoded collection size overflows address space".to_string(),
            })?;
    let remaining = reader.bytes.len().saturating_sub(reader.at);
    if remaining < required {
        return Err(IdentifiabilityError::Canonical {
            at: reader.at,
            detail: format!(
                "truncated {field}: {count} entries require at least {required} bytes, only {remaining} remain"
            ),
        });
    }
    Ok(())
}

fn encode_source_ref(
    writer: &mut CanonicalWriter,
    source: &SourceRef,
) -> Result<(), IdentifiabilityError> {
    encode_source_key(writer, &source.key)?;
    encode_source_kind(writer, source.kind);
    writer.hash(source.expected_hash);
    writer.text(&source.content_hash_domain, "source content-hash domain")?;
    writer.u32(source.contract_version);
    Ok(())
}

fn decode_source_ref(reader: &mut CanonicalReader<'_>) -> Result<SourceRef, IdentifiabilityError> {
    SourceRef::try_new(
        decode_source_key(reader)?,
        decode_source_kind(reader)?,
        reader.hash("source hash")?,
        reader.token("source content-hash domain")?,
        reader.u32("source contract version")?,
    )
}

fn encode_trust_receipt_ref(
    writer: &mut CanonicalWriter,
    trust: &TrustReceiptRef,
) -> Result<(), IdentifiabilityError> {
    encode_source_ref(writer, &trust.receipt)?;
    encode_source_ref(writer, &trust.subject)?;
    match &trust.subject_artifact {
        Some(subject_artifact) => {
            writer.byte(1);
            encode_artifact_id(writer, subject_artifact)?;
        }
        None => writer.byte(0),
    }
    match &trust.authentication {
        TrustAuthentication::Unauthenticated => writer.byte(0),
        TrustAuthentication::IssuerPolicy {
            issuer,
            trust_policy,
        } => {
            writer.byte(1);
            encode_artifact_id(writer, issuer)?;
            encode_source_ref(writer, trust_policy)?;
        }
    }
    Ok(())
}

fn decode_trust_receipt_ref(
    reader: &mut CanonicalReader<'_>,
) -> Result<TrustReceiptRef, IdentifiabilityError> {
    let receipt = decode_source_ref(reader)?;
    let subject = decode_source_ref(reader)?;
    let subject_artifact = match reader.byte("trust receipt subject artifact")? {
        0 => None,
        1 => Some(decode_artifact_id(reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown trust-receipt-subject-artifact tag {tag}"),
            });
        }
    };
    let authentication = match reader.byte("trust receipt authentication")? {
        0 => TrustAuthentication::Unauthenticated,
        1 => TrustAuthentication::IssuerPolicy {
            issuer: decode_artifact_id(reader)?,
            trust_policy: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown trust-receipt-authentication tag {tag}"),
            });
        }
    };
    TrustReceiptRef::try_new_with_subject_artifact(
        receipt,
        subject,
        subject_artifact,
        authentication,
    )
}

fn encode_resolution_verification(writer: &mut CanonicalWriter, verification: &SourceVerification) {
    match verification {
        SourceVerification::TypedArtifact => writer.byte(0),
        SourceVerification::HashPreimage { byte_len } => {
            writer.byte(1);
            writer.u64(*byte_len);
        }
        SourceVerification::Unverified => writer.byte(2),
    }
}

fn decode_resolution_verification(
    reader: &mut CanonicalReader<'_>,
) -> Result<SourceVerification, IdentifiabilityError> {
    match reader.byte("source verification")? {
        0 => Ok(SourceVerification::TypedArtifact),
        1 => {
            let byte_len = reader.u64("verified source byte length")?;
            Ok(SourceVerification::HashPreimage { byte_len })
        }
        2 => Ok(SourceVerification::Unverified),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown source-verification tag {tag}"),
        }),
    }
}

fn encode_resolution_entry(
    writer: &mut CanonicalWriter,
    resolution: &SourceResolution,
) -> Result<(), IdentifiabilityError> {
    encode_source_key(writer, &resolution.key)?;
    encode_source_kind(writer, resolution.kind);
    writer.hash(resolution.resolved_hash);
    writer.text(
        &resolution.content_hash_domain,
        "resolved source content-hash domain",
    )?;
    writer.u32(resolution.contract_version);
    encode_resolution_verification(writer, &resolution.verification);
    match &resolution.authority {
        AuthorityDisposition::ContentVerified => writer.byte(0),
        AuthorityDisposition::ExternalTrustReceipt { trust_receipt } => {
            writer.byte(1);
            encode_trust_receipt_ref(writer, trust_receipt)?;
        }
        AuthorityDisposition::Unverified { reason } => {
            writer.byte(2);
            writer.text(reason, "unverified source reason")?;
        }
    }
    Ok(())
}

fn encode_resolution_set(
    writer: &mut CanonicalWriter,
    resolutions: &SourceResolutionSet,
) -> Result<(), IdentifiabilityError> {
    writer.count(resolutions.entries.len(), "source resolutions")?;
    for resolution in resolutions.entries.values() {
        encode_resolution_entry(writer, resolution)?;
    }
    Ok(())
}

fn decode_resolution_set(
    reader: &mut CanonicalReader<'_>,
) -> Result<SourceResolutionSet, IdentifiabilityError> {
    let count = reader.count("source resolutions")?;
    preflight_collection_bytes(reader, count, 1, "source resolutions")?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        let key = decode_source_key(reader)?;
        let kind = decode_source_kind(reader)?;
        let resolved_hash = reader.hash("resolved source hash")?;
        let content_hash_domain = reader.token("resolved source content-hash domain")?;
        let contract_version = reader.u32("resolved source contract version")?;
        let verification = decode_resolution_verification(reader)?;
        let authority = match reader.byte("source authority")? {
            0 => AuthorityDisposition::ContentVerified,
            1 => AuthorityDisposition::ExternalTrustReceipt {
                trust_receipt: decode_trust_receipt_ref(reader)?,
            },
            2 => AuthorityDisposition::Unverified {
                reason: reader.reason("unverified source reason")?,
            },
            tag => {
                return Err(IdentifiabilityError::Canonical {
                    at: reader.at.saturating_sub(1),
                    detail: format!("unknown source-authority tag {tag}"),
                });
            }
        };
        validate_authority_disposition(&authority)?;
        if !hash_is_nonzero(resolved_hash) || contract_version == 0 {
            return Err(IdentifiabilityError::InvalidText {
                field: "transported source resolution",
                detail: "resolved hash and source contract version must be nonzero".to_string(),
            });
        }
        entries.push(SourceResolution {
            key,
            kind,
            resolved_hash,
            content_hash_domain,
            contract_version,
            authority,
            verification,
        });
    }
    SourceResolutionSet::try_new(entries)
}

fn encode_observation_key(
    writer: &mut CanonicalWriter,
    key: &ObservationKey,
) -> Result<(), IdentifiabilityError> {
    encode_case_id(writer, &key.case)?;
    encode_channel(writer, &key.channel)
}

fn decode_observation_key(
    reader: &mut CanonicalReader<'_>,
) -> Result<ObservationKey, IdentifiabilityError> {
    Ok(ObservationKey::new(
        decode_case_id(reader)?,
        decode_channel(reader)?,
    ))
}

fn encode_parameter_treatment(
    writer: &mut CanonicalWriter,
    treatment: &ParameterTreatment,
) -> Result<(), IdentifiabilityError> {
    match treatment {
        ParameterTreatment::Estimated => writer.byte(0),
        ParameterTreatment::Profiled => writer.byte(1),
        ParameterTreatment::Marginalized => writer.byte(2),
        ParameterTreatment::Conditioned(value) => {
            writer.byte(3);
            writer.f64(value.value_si);
            encode_source_key(writer, &value.source)?;
        }
        ParameterTreatment::Derived {
            definition,
            parents,
        } => {
            writer.byte(4);
            encode_source_key(writer, definition)?;
            writer.count(parents.len(), "derived parameter parents")?;
            for parent in parents {
                encode_role(writer, parent)?;
            }
        }
    }
    Ok(())
}

fn decode_parameter_treatment(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterTreatment, IdentifiabilityError> {
    match reader.byte("parameter treatment")? {
        0 => Ok(ParameterTreatment::Estimated),
        1 => Ok(ParameterTreatment::Profiled),
        2 => Ok(ParameterTreatment::Marginalized),
        3 => Ok(ParameterTreatment::Conditioned(ConditionedValue::try_new(
            reader.f64("conditioned value")?,
            decode_source_key(reader)?,
        )?)),
        4 => {
            let definition = decode_source_key(reader)?;
            let count = reader.count("derived parameter parents")?;
            let mut parents = BTreeSet::new();
            for _ in 0..count {
                let parent = decode_role(reader)?;
                if !parents.insert(parent.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "derived parameter parent",
                        id: parent.to_string(),
                    });
                }
            }
            Ok(ParameterTreatment::Derived {
                definition,
                parents,
            })
        }
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown parameter treatment tag {tag}"),
        }),
    }
}

fn encode_prior_policy(
    writer: &mut CanonicalWriter,
    prior: &PriorPolicy,
) -> Result<(), IdentifiabilityError> {
    match prior {
        PriorPolicy::Distribution(distribution) => {
            writer.byte(0);
            encode_prior(writer, distribution)?;
        }
        PriorPolicy::Absent { reason } => {
            writer.byte(1);
            writer.text(reason, "prior absence reason")?;
        }
        PriorPolicy::NotApplicable { reason } => {
            writer.byte(2);
            writer.text(reason, "prior not-applicable reason")?;
        }
    }
    Ok(())
}

fn decode_prior_policy(
    reader: &mut CanonicalReader<'_>,
) -> Result<PriorPolicy, IdentifiabilityError> {
    match reader.byte("prior policy")? {
        0 => Ok(PriorPolicy::Distribution(decode_prior(reader)?)),
        1 => Ok(PriorPolicy::Absent {
            reason: reader.reason("prior absence reason")?,
        }),
        2 => Ok(PriorPolicy::NotApplicable {
            reason: reader.reason("prior not-applicable reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown prior policy tag {tag}"),
        }),
    }
}

fn encode_owner(
    writer: &mut CanonicalWriter,
    owner: &ParameterOwnerBinding,
) -> Result<(), IdentifiabilityError> {
    match owner {
        ParameterOwnerBinding::ConstitutiveModel => writer.byte(0),
        ParameterOwnerBinding::InitialState { state_path } => {
            writer.byte(1);
            encode_source_key(writer, state_path)?;
        }
        ParameterOwnerBinding::Instrument {
            instrument,
            acquisition_channel,
            metrology,
        } => {
            writer.byte(2);
            encode_artifact_id(writer, instrument)?;
            encode_artifact_id(writer, acquisition_channel)?;
            encode_source_key(writer, metrology)?;
        }
        ParameterOwnerBinding::Discrepancy { family } => {
            writer.byte(3);
            encode_source_key(writer, family)?;
        }
        ParameterOwnerBinding::ControlledInput { protocol } => {
            writer.byte(4);
            encode_source_key(writer, protocol)?;
        }
        ParameterOwnerBinding::Population { hierarchy } => {
            writer.byte(5);
            encode_source_key(writer, hierarchy)?;
        }
    }
    Ok(())
}

fn decode_owner(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterOwnerBinding, IdentifiabilityError> {
    match reader.byte("parameter owner")? {
        0 => Ok(ParameterOwnerBinding::ConstitutiveModel),
        1 => Ok(ParameterOwnerBinding::InitialState {
            state_path: decode_source_key(reader)?,
        }),
        2 => Ok(ParameterOwnerBinding::Instrument {
            instrument: decode_artifact_id(reader)?,
            acquisition_channel: decode_artifact_id(reader)?,
            metrology: decode_source_key(reader)?,
        }),
        3 => Ok(ParameterOwnerBinding::Discrepancy {
            family: decode_source_key(reader)?,
        }),
        4 => Ok(ParameterOwnerBinding::ControlledInput {
            protocol: decode_source_key(reader)?,
        }),
        5 => Ok(ParameterOwnerBinding::Population {
            hierarchy: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown parameter owner tag {tag}"),
        }),
    }
}

fn encode_scope(
    writer: &mut CanonicalWriter,
    scope: &ParameterScopeBinding,
) -> Result<(), IdentifiabilityError> {
    match scope {
        ParameterScopeBinding::Global => writer.byte(0),
        ParameterScopeBinding::Cases(cases) => {
            writer.byte(1);
            writer.count(cases.len(), "parameter case scope")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
        ParameterScopeBinding::MaterialLot { lot, cases } => {
            writer.byte(2);
            encode_artifact_id(writer, lot)?;
            writer.count(cases.len(), "material-lot case scope")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
        ParameterScopeBinding::Specimen { case, specimen } => {
            writer.byte(3);
            encode_case_id(writer, case)?;
            encode_artifact_id(writer, specimen)?;
        }
        ParameterScopeBinding::Field { support, cases } => {
            writer.byte(4);
            encode_source_key(writer, support)?;
            writer.count(cases.len(), "field case scope")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
        ParameterScopeBinding::Hierarchical {
            population,
            level,
            hierarchy,
            cases,
        } => {
            writer.byte(5);
            encode_artifact_id(writer, population)?;
            writer.u32(*level);
            encode_source_key(writer, hierarchy)?;
            writer.count(cases.len(), "hierarchical case scope")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
    }
    Ok(())
}

fn decode_scope(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterScopeBinding, IdentifiabilityError> {
    match reader.byte("parameter scope")? {
        0 => Ok(ParameterScopeBinding::Global),
        1 => {
            let count = reader.count("parameter case scope")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "parameter case scope",
                        id: case.to_string(),
                    });
                }
            }
            Ok(ParameterScopeBinding::Cases(cases))
        }
        2 => {
            let lot = decode_artifact_id(reader)?;
            let count = reader.count("material-lot case scope")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "material-lot case scope",
                        id: case.to_string(),
                    });
                }
            }
            Ok(ParameterScopeBinding::MaterialLot { lot, cases })
        }
        3 => Ok(ParameterScopeBinding::Specimen {
            case: decode_case_id(reader)?,
            specimen: decode_artifact_id(reader)?,
        }),
        4 => {
            let support = decode_source_key(reader)?;
            let count = reader.count("field case scope")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "field case scope",
                        id: case.to_string(),
                    });
                }
            }
            Ok(ParameterScopeBinding::Field { support, cases })
        }
        5 => {
            let population = decode_artifact_id(reader)?;
            let level = reader.u32("hierarchical level")?;
            let hierarchy = decode_source_key(reader)?;
            let count = reader.count("hierarchical case scope")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "hierarchical case scope",
                        id: case.to_string(),
                    });
                }
            }
            Ok(ParameterScopeBinding::Hierarchical {
                population,
                level,
                hierarchy,
                cases,
            })
        }
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown parameter scope tag {tag}"),
        }),
    }
}

fn encode_study_parameter(
    writer: &mut CanonicalWriter,
    parameter: &StudyParameter,
) -> Result<(), IdentifiabilityError> {
    encode_role(writer, &parameter.role)?;
    writer.quantity(parameter.quantity);
    encode_parameter_domain(writer, parameter.domain);
    writer.byte(match parameter.purpose {
        ParameterPurpose::Estimand => 0,
        ParameterPurpose::Nuisance => 1,
        ParameterPurpose::Hyperparameter => 2,
        ParameterPurpose::CalibrationControl => 3,
    });
    encode_parameter_treatment(writer, &parameter.treatment)?;
    encode_owner(writer, &parameter.owner)?;
    encode_scope(writer, &parameter.scope)?;
    encode_prior_policy(writer, &parameter.prior)?;
    match &parameter.influence_coverage {
        InfluenceCoverage::Declared => writer.byte(0),
        InfluenceCoverage::IntentionallyAbsent { reason } => {
            writer.byte(1);
            writer.text(reason, "influence absence reason")?;
        }
        InfluenceCoverage::NotApplicable { reason } => {
            writer.byte(2);
            writer.text(reason, "influence not-applicable reason")?;
        }
    }
    Ok(())
}

fn decode_study_parameter(
    reader: &mut CanonicalReader<'_>,
) -> Result<StudyParameter, IdentifiabilityError> {
    let role = decode_role(reader)?;
    let quantity = reader.quantity("physical parameter quantity")?;
    let domain = decode_parameter_domain(reader)?;
    let purpose = match reader.byte("parameter purpose")? {
        0 => ParameterPurpose::Estimand,
        1 => ParameterPurpose::Nuisance,
        2 => ParameterPurpose::Hyperparameter,
        3 => ParameterPurpose::CalibrationControl,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown parameter purpose tag {tag}"),
            });
        }
    };
    let treatment = decode_parameter_treatment(reader)?;
    let owner = decode_owner(reader)?;
    let scope = decode_scope(reader)?;
    let prior = decode_prior_policy(reader)?;
    let influence_coverage = match reader.byte("influence coverage")? {
        0 => InfluenceCoverage::Declared,
        1 => InfluenceCoverage::IntentionallyAbsent {
            reason: reader.reason("influence absence reason")?,
        },
        2 => InfluenceCoverage::NotApplicable {
            reason: reader.reason("influence not-applicable reason")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown influence coverage tag {tag}"),
            });
        }
    };
    StudyParameter::try_new(
        role,
        quantity,
        domain,
        purpose,
        treatment,
        owner,
        scope,
        prior,
        influence_coverage,
    )
}

fn encode_constraint(
    writer: &mut CanonicalWriter,
    constraint: &JointConstraint,
) -> Result<(), IdentifiabilityError> {
    writer.text(constraint.id.as_str(), "constraint id")?;
    match &constraint.kind {
        JointConstraintKind::Affine {
            terms,
            relation,
            rhs_si,
            residual_quantity,
        } => {
            writer.byte(0);
            writer.count(terms.len(), "affine terms")?;
            for term in terms {
                encode_role(writer, &term.parameter)?;
                writer.f64(term.coefficient);
                writer.quantity(term.coefficient_quantity);
            }
            writer.byte(match relation {
                ConstraintRelation::Equal => 0,
                ConstraintRelation::LessOrEqual => 1,
                ConstraintRelation::GreaterOrEqual => 2,
            });
            writer.f64(*rhs_si);
            writer.quantity(*residual_quantity);
        }
        JointConstraintKind::Simplex {
            members,
            total_si,
            quantity,
        } => {
            writer.byte(1);
            writer.count(members.len(), "simplex members")?;
            for member in members {
                encode_role(writer, member)?;
            }
            writer.f64(*total_si);
            writer.quantity(*quantity);
        }
        JointConstraintKind::Ordered { members, strict } => {
            writer.byte(2);
            writer.count(members.len(), "ordered members")?;
            for member in members {
                encode_role(writer, member)?;
            }
            writer.byte(u8::from(*strict));
        }
        JointConstraintKind::ExternalManifold {
            members,
            definition,
            codimension,
        } => {
            writer.byte(3);
            writer.count(members.len(), "manifold members")?;
            for member in members {
                encode_role(writer, member)?;
            }
            encode_source_key(writer, definition)?;
            match codimension {
                ConstraintCodimension::Finite { codimension } => {
                    writer.byte(0);
                    writer.u64(*codimension);
                }
                ConstraintCodimension::InfiniteDimensional { profile } => {
                    writer.byte(1);
                    encode_source_key(writer, profile)?;
                }
            }
        }
        JointConstraintKind::StochasticCoupling {
            members,
            distribution,
        } => {
            writer.byte(4);
            writer.count(members.len(), "stochastic members")?;
            for member in members {
                encode_role(writer, member)?;
            }
            encode_source_key(writer, distribution)?;
        }
    }
    Ok(())
}

fn decode_role_set(
    reader: &mut CanonicalReader<'_>,
    field: &'static str,
) -> Result<BTreeSet<ParameterRoleId>, IdentifiabilityError> {
    let count = reader.count(field)?;
    let mut result = BTreeSet::new();
    for _ in 0..count {
        let role = decode_role(reader)?;
        if !result.insert(role.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field,
                id: role.to_string(),
            });
        }
    }
    Ok(result)
}

fn decode_constraint(
    reader: &mut CanonicalReader<'_>,
) -> Result<JointConstraint, IdentifiabilityError> {
    let id = ConstraintId::try_new(reader.token("constraint id")?)?;
    let kind = match reader.byte("constraint kind")? {
        0 => {
            let count = reader.count("affine terms")?;
            preflight_collection_bytes(reader, count, 1, "affine terms")?;
            let mut terms = Vec::with_capacity(count);
            for _ in 0..count {
                terms.push(AffineConstraintTerm::try_new(
                    decode_role(reader)?,
                    reader.f64("affine coefficient")?,
                    reader.quantity("affine coefficient quantity")?,
                )?);
            }
            let relation = match reader.byte("constraint relation")? {
                0 => ConstraintRelation::Equal,
                1 => ConstraintRelation::LessOrEqual,
                2 => ConstraintRelation::GreaterOrEqual,
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown constraint relation tag {tag}"),
                    });
                }
            };
            JointConstraintKind::Affine {
                terms,
                relation,
                rhs_si: reader.f64("constraint RHS")?,
                residual_quantity: reader.quantity("constraint residual quantity")?,
            }
        }
        1 => JointConstraintKind::Simplex {
            members: decode_role_set(reader, "simplex members")?,
            total_si: reader.f64("simplex total")?,
            quantity: reader.quantity("simplex quantity")?,
        },
        2 => {
            let count = reader.count("ordered members")?;
            preflight_collection_bytes(reader, count, 4, "ordered members")?;
            let mut members = Vec::with_capacity(count);
            for _ in 0..count {
                members.push(decode_role(reader)?);
            }
            JointConstraintKind::Ordered {
                members,
                strict: match reader.byte("strict ordering")? {
                    0 => false,
                    1 => true,
                    tag => {
                        return Err(IdentifiabilityError::Canonical {
                            at: reader.at.saturating_sub(1),
                            detail: format!("invalid strict-ordering tag {tag}"),
                        });
                    }
                },
            }
        }
        3 => {
            let members = decode_role_set(reader, "manifold members")?;
            let definition = decode_source_key(reader)?;
            let codimension = match reader.byte("manifold codimension")? {
                0 => ConstraintCodimension::Finite {
                    codimension: reader.u64("finite manifold codimension")?,
                },
                1 => ConstraintCodimension::InfiniteDimensional {
                    profile: decode_source_key(reader)?,
                },
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown manifold-codimension tag {tag}"),
                    });
                }
            };
            JointConstraintKind::ExternalManifold {
                members,
                definition,
                codimension,
            }
        }
        4 => JointConstraintKind::StochasticCoupling {
            members: decode_role_set(reader, "stochastic members")?,
            distribution: decode_source_key(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown constraint kind tag {tag}"),
            });
        }
    };
    Ok(JointConstraint::new(id, kind))
}

fn encode_admissible_domain_witness(
    writer: &mut CanonicalWriter,
    witness: &AdmissibleDomainWitness,
) -> Result<(), IdentifiabilityError> {
    writer.count(witness.values.len(), "admissible-domain witness values")?;
    for (role, value) in &witness.values {
        encode_role(writer, role)?;
        writer.f64(*value);
    }
    match &witness.opaque_membership_claim {
        None => writer.byte(0),
        Some(claim) => {
            writer.byte(1);
            encode_source_key(writer, &claim.source)?;
            writer.hash(claim.witness_binding.ok_or_else(|| {
                IdentifiabilityError::InvalidText {
                    field: "admissible-domain witness binding",
                    detail: "an unbound opaque membership claim cannot be encoded".to_string(),
                }
            })?);
        }
    }
    Ok(())
}

fn decode_admissible_domain_witness(
    reader: &mut CanonicalReader<'_>,
) -> Result<AdmissibleDomainWitness, IdentifiabilityError> {
    let value_count = reader.count("admissible-domain witness values")?;
    preflight_collection_bytes(reader, value_count, 12, "admissible-domain witness values")?;
    let mut values = Vec::with_capacity(value_count);
    for _ in 0..value_count {
        values.push((
            decode_role(reader)?,
            reader.f64("admissible-domain witness value")?,
        ));
    }
    let claim = match reader.byte("admissible-domain opaque membership claim")? {
        0 => None,
        1 => Some(OpaqueDomainMembershipClaim::from_bound_source(
            decode_source_key(reader)?,
            reader.hash("admissible-domain witness binding")?,
        )),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown admissible-domain membership-claim tag {tag}"),
            });
        }
    };
    AdmissibleDomainWitness::try_new(values, claim)
}

fn encode_marginal_noise(
    writer: &mut CanonicalWriter,
    noise: &MarginalNoiseSpec,
) -> Result<(), IdentifiabilityError> {
    match noise {
        MarginalNoiseSpec::Gaussian { standard_deviation } => {
            writer.byte(0);
            writer.f64(*standard_deviation);
        }
        MarginalNoiseSpec::StudentT {
            scale,
            degrees_of_freedom,
        } => {
            writer.byte(1);
            writer.f64(*scale);
            writer.f64(*degrees_of_freedom);
        }
        MarginalNoiseSpec::Empirical {
            distribution,
            standard_deviation,
            finite_variance_model,
        } => {
            writer.byte(2);
            encode_source_key(writer, distribution)?;
            writer.f64(*standard_deviation);
            encode_source_key(writer, finite_variance_model)?;
        }
        MarginalNoiseSpec::Bounded { half_width } => {
            writer.byte(3);
            writer.f64(*half_width);
        }
        MarginalNoiseSpec::Unknown { reason } => {
            writer.byte(4);
            writer.text(reason, "unknown noise reason")?;
        }
    }
    Ok(())
}

fn decode_marginal_noise(
    reader: &mut CanonicalReader<'_>,
) -> Result<MarginalNoiseSpec, IdentifiabilityError> {
    match reader.byte("marginal noise")? {
        0 => Ok(MarginalNoiseSpec::Gaussian {
            standard_deviation: reader.f64("Gaussian standard deviation")?,
        }),
        1 => Ok(MarginalNoiseSpec::StudentT {
            scale: reader.f64("Student-t scale")?,
            degrees_of_freedom: reader.f64("Student-t degrees of freedom")?,
        }),
        2 => Ok(MarginalNoiseSpec::Empirical {
            distribution: decode_source_key(reader)?,
            standard_deviation: reader.f64("empirical standard deviation")?,
            finite_variance_model: decode_source_key(reader)?,
        }),
        3 => Ok(MarginalNoiseSpec::Bounded {
            half_width: reader.f64("bounded-noise half width")?,
        }),
        4 => Ok(MarginalNoiseSpec::Unknown {
            reason: reader.reason("unknown noise reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown marginal-noise tag {tag}"),
        }),
    }
}

fn encode_missingness(
    writer: &mut CanonicalWriter,
    missingness: &MissingnessAssumption,
) -> Result<(), IdentifiabilityError> {
    match missingness {
        MissingnessAssumption::Complete { assumption } => {
            writer.byte(0);
            encode_source_key(writer, assumption)?;
        }
        MissingnessAssumption::Modeled { mechanism } => {
            writer.byte(1);
            encode_source_key(writer, mechanism)?;
        }
        MissingnessAssumption::Unknown { reason } => {
            writer.byte(2);
            writer.text(reason, "unknown missingness reason")?;
        }
    }
    Ok(())
}

fn decode_missingness(
    reader: &mut CanonicalReader<'_>,
) -> Result<MissingnessAssumption, IdentifiabilityError> {
    match reader.byte("missingness")? {
        0 => Ok(MissingnessAssumption::Complete {
            assumption: decode_source_key(reader)?,
        }),
        1 => Ok(MissingnessAssumption::Modeled {
            mechanism: decode_source_key(reader)?,
        }),
        2 => Ok(MissingnessAssumption::Unknown {
            reason: reader.reason("unknown missingness reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown missingness tag {tag}"),
        }),
    }
}

fn encode_study_observation(
    writer: &mut CanonicalWriter,
    observation: &StudyObservation,
) -> Result<(), IdentifiabilityError> {
    encode_channel(writer, &observation.id)?;
    encode_qoi_id(writer, &observation.qoi)?;
    writer.text(observation.unit.as_str(), "observation unit")?;
    writer.quantity(observation.quantity);
    encode_source_key(writer, &observation.unit_definition)?;
    encode_frame(writer, &observation.frame)?;
    writer.text(&observation.graph_node, "observation graph node")?;
    writer.text(&observation.graph_port, "observation graph port")?;
    encode_source_key(writer, &observation.operator)?;
    encode_source_key(writer, &observation.aggregation)?;
    encode_source_key(writer, &observation.sensor)?;
    encode_artifact_id(writer, &observation.instrument)?;
    encode_artifact_id(writer, &observation.acquisition_channel)?;
    encode_artifact_id(writer, &observation.clock)?;
    writer.u32(observation.operator_version);
    encode_marginal_noise(writer, &observation.noise)?;
    encode_missingness(writer, &observation.missingness)?;
    match observation.saturation {
        Some(domain) => {
            writer.byte(1);
            encode_parameter_domain(writer, domain);
        }
        None => writer.byte(0),
    }
    writer.u32(observation.protocol_version);
    writer.u32(observation.refinement_version);
    match &observation.rows {
        ObservationRows::Prospective => writer.byte(0),
        ObservationRows::Retrospective(rows) => {
            writer.byte(1);
            writer.count(rows.len(), "observation rows")?;
            for row in rows {
                encode_observation_row_id(writer, row)?;
            }
        }
    }
    Ok(())
}

fn decode_study_observation(
    reader: &mut CanonicalReader<'_>,
) -> Result<StudyObservation, IdentifiabilityError> {
    let id = decode_channel(reader)?;
    let qoi = decode_qoi_id(reader)?;
    let unit = UnitId::try_new(reader.token("observation unit")?).map_err(|error| {
        IdentifiabilityError::Vv {
            detail: error.to_string(),
        }
    })?;
    let quantity = reader.quantity("observation quantity")?;
    let unit_definition = decode_source_key(reader)?;
    let frame = decode_frame(reader)?;
    let graph_node = reader.token("observation graph node")?;
    let graph_port = reader.token("observation graph port")?;
    let operator = decode_source_key(reader)?;
    let aggregation = decode_source_key(reader)?;
    let sensor = decode_source_key(reader)?;
    let instrument = decode_artifact_id(reader)?;
    let acquisition_channel = decode_artifact_id(reader)?;
    let clock = decode_artifact_id(reader)?;
    let operator_version = reader.u32("observation operator version")?;
    let noise = decode_marginal_noise(reader)?;
    let missingness = decode_missingness(reader)?;
    let saturation = match reader.byte("saturation option")? {
        0 => None,
        1 => Some(decode_parameter_domain(reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("invalid saturation option tag {tag}"),
            });
        }
    };
    let protocol_version = reader.u32("observation protocol version")?;
    let refinement_version = reader.u32("observation refinement version")?;
    let rows = match reader.byte("observation rows")? {
        0 => ObservationRows::Prospective,
        1 => {
            let count = reader.count("observation rows")?;
            let mut rows = BTreeSet::new();
            for _ in 0..count {
                let row = decode_observation_row_id(reader)?;
                if !rows.insert(row.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "observation row",
                        id: row.as_str().to_string(),
                    });
                }
            }
            ObservationRows::Retrospective(rows)
        }
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("invalid observation-row tag {tag}"),
            });
        }
    };
    StudyObservation::try_new(
        id,
        qoi,
        unit,
        quantity,
        unit_definition,
        frame,
        graph_node,
        graph_port,
        operator,
        aggregation,
        sensor,
        instrument,
        acquisition_channel,
        clock,
        operator_version,
        noise,
        missingness,
        saturation,
        protocol_version,
        refinement_version,
        rows,
    )
}

fn encode_discrepancy(
    writer: &mut CanonicalWriter,
    discrepancy: &StudyDiscrepancy,
) -> Result<(), IdentifiabilityError> {
    match discrepancy {
        StudyDiscrepancy::Uncharacterized { reason } => {
            writer.byte(0);
            writer.text(reason, "uncharacterized discrepancy reason")?;
        }
        StudyDiscrepancy::NotApplicable { basis } => {
            writer.byte(1);
            match basis {
                DiscrepancyInapplicability::PhysicalApplicability { assumption } => {
                    writer.byte(0);
                    encode_source_key(writer, assumption)?;
                }
                DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                    generator,
                    producer,
                    production_binding,
                    assumption,
                } => {
                    writer.byte(1);
                    encode_source_key(writer, generator)?;
                    encode_artifact_id(writer, producer)?;
                    encode_source_key(writer, production_binding)?;
                    encode_source_key(writer, assumption)?;
                }
                DiscrepancyInapplicability::ProspectiveDesign { assumption } => {
                    writer.byte(2);
                    encode_source_key(writer, assumption)?;
                }
            }
        }
        StudyDiscrepancy::AssumedZero { assumption } => {
            writer.byte(2);
            encode_source_key(writer, assumption)?;
        }
        StudyDiscrepancy::Modeled {
            family,
            parameters,
            support,
            confounding_guard,
        } => {
            writer.byte(3);
            encode_source_key(writer, family)?;
            writer.count(parameters.len(), "discrepancy parameters")?;
            for parameter in parameters {
                encode_role(writer, parameter)?;
            }
            encode_source_key(writer, support)?;
            encode_source_key(writer, confounding_guard)?;
        }
    }
    Ok(())
}

fn decode_discrepancy(
    reader: &mut CanonicalReader<'_>,
) -> Result<StudyDiscrepancy, IdentifiabilityError> {
    match reader.byte("discrepancy")? {
        0 => Ok(StudyDiscrepancy::Uncharacterized {
            reason: reader.reason("uncharacterized discrepancy reason")?,
        }),
        1 => {
            let basis = match reader.byte("discrepancy inapplicability basis")? {
                0 => DiscrepancyInapplicability::PhysicalApplicability {
                    assumption: decode_source_key(reader)?,
                },
                1 => DiscrepancyInapplicability::DeclaredSyntheticSelfModel {
                    generator: decode_source_key(reader)?,
                    producer: decode_artifact_id(reader)?,
                    production_binding: decode_source_key(reader)?,
                    assumption: decode_source_key(reader)?,
                },
                2 => DiscrepancyInapplicability::ProspectiveDesign {
                    assumption: decode_source_key(reader)?,
                },
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown discrepancy-inapplicability tag {tag}"),
                    });
                }
            };
            Ok(StudyDiscrepancy::NotApplicable { basis })
        }
        2 => Ok(StudyDiscrepancy::AssumedZero {
            assumption: decode_source_key(reader)?,
        }),
        3 => Ok(StudyDiscrepancy::Modeled {
            family: decode_source_key(reader)?,
            parameters: decode_role_set(reader, "discrepancy parameters")?,
            support: decode_source_key(reader)?,
            confounding_guard: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown discrepancy tag {tag}"),
        }),
    }
}

fn encode_case_physics_sources(
    writer: &mut CanonicalWriter,
    sources: &CasePhysicsSources,
) -> Result<(), IdentifiabilityError> {
    for key in [
        &sources.frame_transform,
        &sources.specimen_geometry,
        &sources.specimen_process,
        &sources.specimen_preparation,
        &sources.load_path,
        &sources.environment_path,
        &sources.time_grid,
    ] {
        encode_source_key(writer, key)?;
    }
    match &sources.initial_state {
        Some(key) => {
            writer.byte(1);
            encode_source_key(writer, key)?;
        }
        None => writer.byte(0),
    }
    Ok(())
}

fn decode_case_physics_sources(
    reader: &mut CanonicalReader<'_>,
) -> Result<CasePhysicsSources, IdentifiabilityError> {
    let frame_transform = decode_source_key(reader)?;
    let specimen_geometry = decode_source_key(reader)?;
    let specimen_process = decode_source_key(reader)?;
    let specimen_preparation = decode_source_key(reader)?;
    let load_path = decode_source_key(reader)?;
    let environment_path = decode_source_key(reader)?;
    let time_grid = decode_source_key(reader)?;
    let initial_state = match reader.byte("case initial-state source option")? {
        0 => None,
        1 => Some(decode_source_key(reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown case initial-state source option tag {tag}"),
            });
        }
    };
    Ok(CasePhysicsSources::new(
        frame_transform,
        specimen_geometry,
        specimen_process,
        specimen_preparation,
        load_path,
        environment_path,
        time_grid,
        initial_state,
    ))
}

fn encode_observation_sharing_group(
    writer: &mut CanonicalWriter,
    group: &ObservationSharingGroup,
) -> Result<(), IdentifiabilityError> {
    writer.count(group.channels.len(), "observation-sharing channels")?;
    for channel in &group.channels {
        encode_channel(writer, channel)?;
    }
    writer.count(group.rows.len(), "observation-sharing rows")?;
    for row in &group.rows {
        encode_observation_row_id(writer, row)?;
    }
    encode_source_key(writer, &group.joint_likelihood)?;
    writer.text(&group.justification, "observation-sharing justification")
}

fn decode_observation_sharing_group(
    reader: &mut CanonicalReader<'_>,
) -> Result<ObservationSharingGroup, IdentifiabilityError> {
    let channel_count = reader.count("observation-sharing channels")?;
    let mut channels = BTreeSet::new();
    for _ in 0..channel_count {
        let channel = decode_channel(reader)?;
        if !channels.insert(channel.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "observation-sharing channel",
                id: channel.to_string(),
            });
        }
    }
    let row_count = reader.count("observation-sharing rows")?;
    let mut rows = BTreeSet::new();
    for _ in 0..row_count {
        let row = decode_observation_row_id(reader)?;
        if !rows.insert(row.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "observation-sharing row",
                id: row.as_str().to_string(),
            });
        }
    }
    ObservationSharingGroup::try_new(
        channels,
        rows,
        decode_source_key(reader)?,
        reader.reason("observation-sharing justification")?,
    )
}

fn encode_case(
    writer: &mut CanonicalWriter,
    case: &StudyCaseDocument,
) -> Result<(), IdentifiabilityError> {
    encode_case_id(writer, &case.id)?;
    match &case.purpose {
        CasePurpose::Calibration => writer.byte(0),
        CasePurpose::SymmetryBreaking => writer.byte(1),
        CasePurpose::ValidationOnly => writer.byte(2),
        CasePurpose::BlindFalsification => writer.byte(3),
        CasePurpose::ProspectiveDesign => writer.byte(4),
        CasePurpose::Complementary { reason } => {
            writer.byte(5);
            writer.text(reason, "complementary case reason")?;
        }
    }
    encode_initial_state(writer, case.initial_state);
    encode_specimen(writer, &case.specimen)?;
    encode_protocol(writer, &case.protocol)?;
    encode_case_physics_sources(writer, &case.physics_sources)?;
    encode_source_key(writer, &case.forward_model)?;
    match &case.data {
        CaseDataDeclaration::Prospective => writer.byte(0),
        CaseDataDeclaration::Retrospective {
            experiment,
            split,
            parser,
            preprocessing,
            parser_version,
            split_grouping,
        } => {
            writer.byte(1);
            encode_source_key(writer, experiment)?;
            encode_source_key(writer, split)?;
            encode_source_key(writer, parser)?;
            encode_source_key(writer, preprocessing)?;
            writer.u32(*parser_version);
            encode_artifact_id(writer, split_grouping)?;
        }
    }
    writer.count(case.observations.len(), "case observations")?;
    for observation in case.observations.values() {
        encode_study_observation(writer, observation)?;
    }
    writer.count(case.discrepancies.len(), "case discrepancies")?;
    for (channel, discrepancy) in &case.discrepancies {
        encode_channel(writer, channel)?;
        encode_discrepancy(writer, discrepancy)?;
    }
    writer.count(case.observation_sharing.len(), "observation-sharing groups")?;
    for group in &case.observation_sharing {
        encode_observation_sharing_group(writer, group)?;
    }
    Ok(())
}

fn decode_case(
    reader: &mut CanonicalReader<'_>,
) -> Result<StudyCaseDocument, IdentifiabilityError> {
    let id = decode_case_id(reader)?;
    let purpose = match reader.byte("case purpose")? {
        0 => CasePurpose::Calibration,
        1 => CasePurpose::SymmetryBreaking,
        2 => CasePurpose::ValidationOnly,
        3 => CasePurpose::BlindFalsification,
        4 => CasePurpose::ProspectiveDesign,
        5 => CasePurpose::Complementary {
            reason: reader.reason("complementary case reason")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown case-purpose tag {tag}"),
            });
        }
    };
    let initial_state = decode_initial_state(reader)?;
    let specimen = decode_specimen(reader)?;
    let protocol = decode_protocol(reader)?;
    let physics_sources = decode_case_physics_sources(reader)?;
    let forward_model = decode_source_key(reader)?;
    let data = match reader.byte("case data")? {
        0 => CaseDataDeclaration::Prospective,
        1 => CaseDataDeclaration::Retrospective {
            experiment: decode_source_key(reader)?,
            split: decode_source_key(reader)?,
            parser: decode_source_key(reader)?,
            preprocessing: decode_source_key(reader)?,
            parser_version: reader.u32("parser version")?,
            split_grouping: decode_artifact_id(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown case-data tag {tag}"),
            });
        }
    };
    let observation_count = reader.count("case observations")?;
    preflight_collection_bytes(reader, observation_count, 1, "case observations")?;
    let mut observations = Vec::with_capacity(observation_count);
    for _ in 0..observation_count {
        observations.push(decode_study_observation(reader)?);
    }
    let discrepancy_count = reader.count("case discrepancies")?;
    preflight_collection_bytes(reader, discrepancy_count, 1, "case discrepancies")?;
    let mut discrepancies = Vec::with_capacity(discrepancy_count);
    for _ in 0..discrepancy_count {
        discrepancies.push((decode_channel(reader)?, decode_discrepancy(reader)?));
    }
    let sharing_count = reader.count("observation-sharing groups")?;
    preflight_collection_bytes(reader, sharing_count, 1, "observation-sharing groups")?;
    let mut observation_sharing = Vec::with_capacity(sharing_count);
    for _ in 0..sharing_count {
        observation_sharing.push(decode_observation_sharing_group(reader)?);
    }
    StudyCaseDocument::try_new(
        id,
        purpose,
        initial_state,
        specimen,
        protocol,
        physics_sources,
        forward_model,
        data,
        observations,
        discrepancies,
        observation_sharing,
    )
}

fn encode_functional(
    writer: &mut CanonicalWriter,
    functional: &DistributionFunctional,
) -> Result<(), IdentifiabilityError> {
    match functional {
        DistributionFunctional::Location { observation } => {
            writer.byte(0);
            encode_observation_key(writer, observation)?;
        }
        DistributionFunctional::LogScale { observation } => {
            writer.byte(1);
            encode_observation_key(writer, observation)?;
        }
        DistributionFunctional::Correlation { left, right } => {
            writer.byte(2);
            encode_observation_key(writer, left)?;
            encode_observation_key(writer, right)?;
        }
        DistributionFunctional::MissingnessLogit { observation } => {
            writer.byte(3);
            encode_observation_key(writer, observation)?;
        }
        DistributionFunctional::CensoringLogit { observation } => {
            writer.byte(4);
            encode_observation_key(writer, observation)?;
        }
    }
    Ok(())
}

fn decode_functional(
    reader: &mut CanonicalReader<'_>,
) -> Result<DistributionFunctional, IdentifiabilityError> {
    match reader.byte("distribution functional")? {
        0 => Ok(DistributionFunctional::Location {
            observation: decode_observation_key(reader)?,
        }),
        1 => Ok(DistributionFunctional::LogScale {
            observation: decode_observation_key(reader)?,
        }),
        2 => Ok(DistributionFunctional::Correlation {
            left: decode_observation_key(reader)?,
            right: decode_observation_key(reader)?,
        }),
        3 => Ok(DistributionFunctional::MissingnessLogit {
            observation: decode_observation_key(reader)?,
        }),
        4 => Ok(DistributionFunctional::CensoringLogit {
            observation: decode_observation_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown distribution-functional tag {tag}"),
        }),
    }
}

fn encode_influence(
    writer: &mut CanonicalWriter,
    influence: &InfluenceDeclaration,
) -> Result<(), IdentifiabilityError> {
    writer.text(influence.id.as_str(), "influence id")?;
    encode_role(writer, &influence.parameter)?;
    encode_functional(writer, &influence.functional)?;
    match &influence.representation {
        InfluenceRepresentation::Direct => writer.byte(0),
        InfluenceRepresentation::StateMediated { state_path } => {
            writer.byte(1);
            encode_source_key(writer, state_path)?;
        }
        InfluenceRepresentation::Composite { operator, inputs } => {
            writer.byte(2);
            encode_source_key(writer, operator)?;
            writer.count(inputs.len(), "composite influence inputs")?;
            for input in inputs {
                writer.text(input.as_str(), "influence id")?;
            }
        }
        InfluenceRepresentation::ExternalDefinition { definition } => {
            writer.byte(3);
            encode_source_key(writer, definition)?;
        }
    }
    Ok(())
}

fn decode_influence(
    reader: &mut CanonicalReader<'_>,
) -> Result<InfluenceDeclaration, IdentifiabilityError> {
    let id = InfluenceId::try_new(reader.token("influence id")?)?;
    let parameter = decode_role(reader)?;
    let functional = decode_functional(reader)?;
    let representation = match reader.byte("influence representation")? {
        0 => InfluenceRepresentation::Direct,
        1 => InfluenceRepresentation::StateMediated {
            state_path: decode_source_key(reader)?,
        },
        2 => {
            let operator = decode_source_key(reader)?;
            let count = reader.count("composite influence inputs")?;
            let mut inputs = BTreeSet::new();
            for _ in 0..count {
                let input = InfluenceId::try_new(reader.token("influence id")?)?;
                if !inputs.insert(input.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "composite influence input",
                        id: input.to_string(),
                    });
                }
            }
            InfluenceRepresentation::Composite { operator, inputs }
        }
        3 => InfluenceRepresentation::ExternalDefinition {
            definition: decode_source_key(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown influence-representation tag {tag}"),
            });
        }
    };
    Ok(InfluenceDeclaration::new(
        id,
        parameter,
        functional,
        representation,
    ))
}

fn encode_gauge_information_regime(
    writer: &mut CanonicalWriter,
    regime: &GaugeInformationRegime,
) -> Result<(), IdentifiabilityError> {
    match regime {
        GaugeInformationRegime::StructuralExactModel => writer.byte(0),
        GaugeInformationRegime::ExactInputOutputMap => writer.byte(1),
        GaugeInformationRegime::NoisyFiniteData => writer.byte(2),
        GaugeInformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            writer.byte(3);
            encode_source_key(writer, joint_prior)?;
        }
    }
    Ok(())
}

fn decode_gauge_information_regime(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeInformationRegime, IdentifiabilityError> {
    match reader.byte("gauge information regime")? {
        0 => Ok(GaugeInformationRegime::StructuralExactModel),
        1 => Ok(GaugeInformationRegime::ExactInputOutputMap),
        2 => Ok(GaugeInformationRegime::NoisyFiniteData),
        3 => Ok(GaugeInformationRegime::PosteriorUnderDeclaredPrior {
            joint_prior: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-information-regime tag {tag}"),
        }),
    }
}

fn encode_gauge_continuous_dimension(
    writer: &mut CanonicalWriter,
    dimension: &GaugeContinuousDimension,
) -> Result<(), IdentifiabilityError> {
    match dimension {
        GaugeContinuousDimension::Finite { dimension } => {
            writer.byte(0);
            writer.u64(*dimension);
        }
        GaugeContinuousDimension::InfiniteDimensional { model_space } => {
            writer.byte(1);
            encode_source_key(writer, model_space)?;
        }
    }
    Ok(())
}

fn decode_gauge_continuous_dimension(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeContinuousDimension, IdentifiabilityError> {
    match reader.byte("gauge continuous dimension")? {
        0 => Ok(GaugeContinuousDimension::Finite {
            dimension: reader.u64("finite gauge dimension")?,
        }),
        1 => Ok(GaugeContinuousDimension::InfiniteDimensional {
            model_space: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-continuous-dimension tag {tag}"),
        }),
    }
}

fn encode_gauge_discrete_size(
    writer: &mut CanonicalWriter,
    size: &GaugeDiscreteSize,
) -> Result<(), IdentifiabilityError> {
    match size {
        GaugeDiscreteSize::Finite { order } => {
            writer.byte(0);
            writer.u64(*order);
        }
        GaugeDiscreteSize::CountablyInfinite { presentation } => {
            writer.byte(1);
            encode_source_key(writer, presentation)?;
        }
    }
    Ok(())
}

fn decode_gauge_discrete_size(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeDiscreteSize, IdentifiabilityError> {
    match reader.byte("gauge discrete size")? {
        0 => Ok(GaugeDiscreteSize::Finite {
            order: reader.u64("finite gauge order")?,
        }),
        1 => Ok(GaugeDiscreteSize::CountablyInfinite {
            presentation: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-discrete-size tag {tag}"),
        }),
    }
}

fn encode_gauge_algebra(
    writer: &mut CanonicalWriter,
    algebra: &GaugeAlgebra,
) -> Result<(), IdentifiabilityError> {
    match algebra {
        GaugeAlgebra::Continuous { group_dimension } => {
            writer.byte(0);
            encode_gauge_continuous_dimension(writer, group_dimension)?;
        }
        GaugeAlgebra::Discrete { size } => {
            writer.byte(1);
            encode_gauge_discrete_size(writer, size)?;
        }
        GaugeAlgebra::Mixed {
            continuous_group_dimension,
            component_group,
        } => {
            writer.byte(2);
            encode_gauge_continuous_dimension(writer, continuous_group_dimension)?;
            encode_gauge_discrete_size(writer, component_group)?;
        }
    }
    Ok(())
}

fn decode_gauge_algebra(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeAlgebra, IdentifiabilityError> {
    match reader.byte("gauge algebra")? {
        0 => Ok(GaugeAlgebra::Continuous {
            group_dimension: decode_gauge_continuous_dimension(reader)?,
        }),
        1 => Ok(GaugeAlgebra::Discrete {
            size: decode_gauge_discrete_size(reader)?,
        }),
        2 => Ok(GaugeAlgebra::Mixed {
            continuous_group_dimension: decode_gauge_continuous_dimension(reader)?,
            component_group: decode_gauge_discrete_size(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-algebra tag {tag}"),
        }),
    }
}

fn encode_gauge_discrete_orbit(
    writer: &mut CanonicalWriter,
    orbit: &GaugeDiscreteOrbitCardinality,
) -> Result<(), IdentifiabilityError> {
    match orbit {
        GaugeDiscreteOrbitCardinality::Finite { cardinality } => {
            writer.byte(0);
            writer.u64(*cardinality);
        }
        GaugeDiscreteOrbitCardinality::CountablyInfinite { presentation } => {
            writer.byte(1);
            encode_source_key(writer, presentation)?;
        }
    }
    Ok(())
}

fn decode_gauge_discrete_orbit(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeDiscreteOrbitCardinality, IdentifiabilityError> {
    match reader.byte("gauge discrete orbit")? {
        0 => Ok(GaugeDiscreteOrbitCardinality::Finite {
            cardinality: reader.u64("finite gauge orbit cardinality")?,
        }),
        1 => Ok(GaugeDiscreteOrbitCardinality::CountablyInfinite {
            presentation: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-discrete-orbit tag {tag}"),
        }),
    }
}

fn encode_regular_gauge_orbit(
    writer: &mut CanonicalWriter,
    orbit: &RegularGaugeOrbit,
) -> Result<(), IdentifiabilityError> {
    encode_gauge_continuous_dimension(writer, &orbit.continuous_orbit_dimension)?;
    encode_gauge_discrete_orbit(writer, &orbit.discrete_orbit_cardinality)
}

fn decode_regular_gauge_orbit(
    reader: &mut CanonicalReader<'_>,
) -> Result<RegularGaugeOrbit, IdentifiabilityError> {
    Ok(RegularGaugeOrbit::new(
        decode_gauge_continuous_dimension(reader)?,
        decode_gauge_discrete_orbit(reader)?,
    ))
}

fn encode_gauge_orbit_geometry(
    writer: &mut CanonicalWriter,
    geometry: &GaugeOrbitGeometry,
) -> Result<(), IdentifiabilityError> {
    match geometry {
        GaugeOrbitGeometry::Regular {
            principal,
            stabilizer_profile,
        } => {
            writer.byte(0);
            encode_regular_gauge_orbit(writer, principal)?;
            writer.byte(u8::from(stabilizer_profile.is_some()));
            if let Some(profile) = stabilizer_profile {
                encode_source_key(writer, profile)?;
            }
        }
        GaugeOrbitGeometry::Stratified {
            principal,
            orbit_type_stabilizer_profile,
        } => {
            writer.byte(1);
            encode_regular_gauge_orbit(writer, principal)?;
            encode_source_key(writer, orbit_type_stabilizer_profile)?;
        }
    }
    Ok(())
}

fn decode_gauge_orbit_geometry(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeOrbitGeometry, IdentifiabilityError> {
    match reader.byte("gauge orbit geometry")? {
        0 => {
            let principal = decode_regular_gauge_orbit(reader)?;
            let stabilizer_profile = match reader.byte("gauge stabilizer profile presence")? {
                0 => None,
                1 => Some(decode_source_key(reader)?),
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown stabilizer-profile-presence tag {tag}"),
                    });
                }
            };
            Ok(GaugeOrbitGeometry::Regular {
                principal,
                stabilizer_profile,
            })
        }
        1 => Ok(GaugeOrbitGeometry::Stratified {
            principal: decode_regular_gauge_orbit(reader)?,
            orbit_type_stabilizer_profile: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-orbit-geometry tag {tag}"),
        }),
    }
}

fn encode_gauge_status(
    writer: &mut CanonicalWriter,
    status: &GaugeStatus,
) -> Result<(), IdentifiabilityError> {
    match status {
        GaugeStatus::Candidate { rationale } => {
            writer.byte(0);
            encode_source_key(writer, rationale)?;
        }
        GaugeStatus::Assumed { assumption } => {
            writer.byte(1);
            encode_source_key(writer, assumption)?;
        }
    }
    Ok(())
}

fn decode_gauge_status(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeStatus, IdentifiabilityError> {
    match reader.byte("gauge status")? {
        0 => Ok(GaugeStatus::Candidate {
            rationale: decode_source_key(reader)?,
        }),
        1 => Ok(GaugeStatus::Assumed {
            assumption: decode_source_key(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-status tag {tag}"),
        }),
    }
}

fn encode_gauge_axes(
    writer: &mut CanonicalWriter,
    axes: &GaugeApplicabilityAxes,
) -> Result<(), IdentifiabilityError> {
    encode_gauge_information_regime(writer, &axes.information)?;
    match &axes.scalar_domain {
        GaugeScalarDomain::Real => writer.byte(0),
        GaugeScalarDomain::Complex { extension } => {
            writer.byte(1);
            encode_source_key(writer, extension)?;
        }
        GaugeScalarDomain::MixedDiscreteContinuous { stratification } => {
            writer.byte(2);
            encode_source_key(writer, stratification)?;
        }
    }
    match &axes.locus {
        GaugeLocus::WholeDomain => writer.byte(0),
        GaugeLocus::Stratum { definition } => {
            writer.byte(1);
            encode_source_key(writer, definition)?;
        }
    }
    match &axes.quantifier {
        GaugeQuantifierScope::AtRealization { realization } => {
            writer.byte(0);
            encode_source_key(writer, realization)?;
        }
        GaugeQuantifierScope::AlmostEverywhere { measure } => {
            writer.byte(1);
            encode_source_key(writer, measure)?;
        }
        GaugeQuantifierScope::ForAll { domain } => {
            writer.byte(2);
            encode_source_key(writer, domain)?;
        }
        GaugeQuantifierScope::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            writer.byte(3);
            writer.u64(probability.0);
            encode_source_key(writer, measure)?;
        }
    }
    Ok(())
}

fn decode_gauge_axes(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeApplicabilityAxes, IdentifiabilityError> {
    let information = decode_gauge_information_regime(reader)?;
    let scalar_domain = match reader.byte("gauge scalar domain")? {
        0 => GaugeScalarDomain::Real,
        1 => GaugeScalarDomain::Complex {
            extension: decode_source_key(reader)?,
        },
        2 => GaugeScalarDomain::MixedDiscreteContinuous {
            stratification: decode_source_key(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-scalar-domain tag {tag}"),
            });
        }
    };
    let locus = match reader.byte("gauge locus")? {
        0 => GaugeLocus::WholeDomain,
        1 => GaugeLocus::Stratum {
            definition: decode_source_key(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-locus tag {tag}"),
            });
        }
    };
    let quantifier = match reader.byte("gauge quantifier")? {
        0 => GaugeQuantifierScope::AtRealization {
            realization: decode_source_key(reader)?,
        },
        1 => GaugeQuantifierScope::AlmostEverywhere {
            measure: decode_source_key(reader)?,
        },
        2 => GaugeQuantifierScope::ForAll {
            domain: decode_source_key(reader)?,
        },
        3 => GaugeQuantifierScope::ProbabilityAtLeast {
            probability: GaugeProbabilityThreshold::try_new(f64::from_bits(
                reader.u64("gauge probability threshold")?,
            ))?,
            measure: decode_source_key(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-quantifier tag {tag}"),
            });
        }
    };
    Ok(GaugeApplicabilityAxes::new(
        information,
        scalar_domain,
        locus,
        quantifier,
    ))
}

fn encode_gauge_validity_scope(
    writer: &mut CanonicalWriter,
    validity: &GaugeValidityScope,
) -> Result<(), IdentifiabilityError> {
    writer.count(validity.cells.len(), "gauge applicability cells")?;
    for (axes, cell) in &validity.cells {
        encode_gauge_axes(writer, axes)?;
        writer.count(
            cell.case_obstruction_support.len(),
            "gauge cell case support",
        )?;
        for (case, support) in &cell.case_obstruction_support {
            encode_case_id(writer, case)?;
            writer.count(
                support.local_obstruction_parameters.len(),
                "local gauge obstruction support",
            )?;
            for role in &support.local_obstruction_parameters {
                encode_role(writer, role)?;
            }
            writer.count(
                support.global_obstruction_parameters.len(),
                "global gauge obstruction support",
            )?;
            for role in &support.global_obstruction_parameters {
                encode_role(writer, role)?;
            }
        }
    }
    Ok(())
}

fn decode_gauge_validity_scope(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeValidityScope, IdentifiabilityError> {
    let count = reader.count("gauge applicability cells")?;
    let mut cells = BTreeMap::new();
    for _ in 0..count {
        let axes = decode_gauge_axes(reader)?;
        let case_count = reader.count("gauge cell case support")?;
        let mut case_support = BTreeMap::new();
        for _ in 0..case_count {
            let case = decode_case_id(reader)?;
            let local = decode_role_set(reader, "local gauge obstruction support")?;
            let global = decode_role_set(reader, "global gauge obstruction support")?;
            let support = GaugeExtentSupport::try_new(local, global)?;
            if case_support.insert(case.clone(), support).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "gauge cell case",
                    id: case.to_string(),
                });
            }
        }
        let cell = GaugeCellDomain::try_new(case_support)?;
        if cells.insert(axes, cell).is_some() {
            return Err(IdentifiabilityError::Duplicate {
                field: "gauge applicability axes",
                id: "duplicate canonical axes".to_string(),
            });
        }
    }
    GaugeValidityScope::try_new(cells)
}

fn encode_gauge(
    writer: &mut CanonicalWriter,
    gauge: &GaugeDeclaration,
) -> Result<(), IdentifiabilityError> {
    writer.text(gauge.id.as_str(), "gauge id")?;
    writer.count(gauge.members.len(), "gauge members")?;
    for member in &gauge.members {
        encode_role(writer, member)?;
    }
    encode_source_key(writer, &gauge.action)?;
    encode_gauge_algebra(writer, &gauge.algebra)?;
    encode_gauge_orbit_geometry(writer, &gauge.orbit_geometry)?;
    encode_gauge_status(writer, &gauge.status)?;
    encode_gauge_validity_scope(writer, &gauge.validity)
}

fn decode_gauge(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeDeclaration, IdentifiabilityError> {
    GaugeDeclaration::try_new(
        GaugeClassId::try_new(reader.token("gauge id")?)?,
        decode_role_set(reader, "gauge members")?,
        decode_source_key(reader)?,
        decode_gauge_algebra(reader)?,
        decode_gauge_orbit_geometry(reader)?,
        decode_gauge_status(reader)?,
        decode_gauge_validity_scope(reader)?,
    )
}

fn encode_gauge_composition(
    writer: &mut CanonicalWriter,
    composition: &GaugeCompositionDeclaration,
) -> Result<(), IdentifiabilityError> {
    writer.text(composition.id.as_str(), "gauge composition id")?;
    writer.count(composition.members.len(), "gauge composition members")?;
    for member in &composition.members {
        writer.text(member.as_str(), "gauge id")?;
    }
    writer.byte(match &composition.kind {
        GaugeCompositionKind::IndependentProduct => 0,
        GaugeCompositionKind::Generated => 1,
    });
    encode_source_key(writer, &composition.law)?;
    encode_gauge_algebra(writer, &composition.effective_algebra)?;
    encode_gauge_orbit_geometry(writer, &composition.effective_orbit_geometry)?;
    encode_gauge_status(writer, &composition.status)?;
    encode_gauge_validity_scope(writer, &composition.validity)
}

fn decode_gauge_composition(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeCompositionDeclaration, IdentifiabilityError> {
    let id = GaugeCompositionId::try_new(reader.token("gauge composition id")?)?;
    let count = reader.count("gauge composition members")?;
    let mut members = BTreeSet::new();
    for _ in 0..count {
        let member = GaugeClassId::try_new(reader.token("gauge id")?)?;
        if !members.insert(member.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "gauge composition member",
                id: member.to_string(),
            });
        }
    }
    let kind = match reader.byte("gauge composition kind")? {
        0 => GaugeCompositionKind::IndependentProduct,
        1 => GaugeCompositionKind::Generated,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-composition-kind tag {tag}"),
            });
        }
    };
    GaugeCompositionDeclaration::try_new(
        id,
        members,
        kind,
        decode_source_key(reader)?,
        decode_gauge_algebra(reader)?,
        decode_gauge_orbit_geometry(reader)?,
        decode_gauge_status(reader)?,
        decode_gauge_validity_scope(reader)?,
    )
}

fn encode_gauge_action_reference(
    writer: &mut CanonicalWriter,
    action: &GaugeActionReference,
) -> Result<(), IdentifiabilityError> {
    match action {
        GaugeActionReference::Single(id) => {
            writer.byte(0);
            writer.text(id.as_str(), "gauge id")?;
        }
        GaugeActionReference::Product(id) => {
            writer.byte(1);
            writer.text(id.as_str(), "gauge composition id")?;
        }
        GaugeActionReference::Composition(id) => {
            writer.byte(2);
            writer.text(id.as_str(), "gauge composition id")?;
        }
    }
    Ok(())
}

fn decode_gauge_action_reference(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeActionReference, IdentifiabilityError> {
    match reader.byte("gauge action reference")? {
        0 => Ok(GaugeActionReference::Single(GaugeClassId::try_new(
            reader.token("gauge id")?,
        )?)),
        1 => Ok(GaugeActionReference::Product(GaugeCompositionId::try_new(
            reader.token("gauge composition id")?,
        )?)),
        2 => Ok(GaugeActionReference::Composition(
            GaugeCompositionId::try_new(reader.token("gauge composition id")?)?,
        )),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-action-reference tag {tag}"),
        }),
    }
}

fn encode_gauge_slice(
    writer: &mut CanonicalWriter,
    slice: &GaugeSlicePlan,
) -> Result<(), IdentifiabilityError> {
    writer.count(slice.support.len(), "gauge slice support")?;
    for role in &slice.support {
        encode_role(writer, role)?;
    }
    encode_source_ref(writer, &slice.constraint)?;
    match &slice.expected_codimension {
        GaugeSliceCodimension::FixedFinite { codimension } => {
            writer.byte(0);
            writer.u64(*codimension);
        }
        GaugeSliceCodimension::FixedInfinite {
            codimension_model,
            compatibility,
        } => {
            writer.byte(1);
            encode_source_ref(writer, codimension_model)?;
            encode_source_ref(writer, compatibility)?;
        }
        GaugeSliceCodimension::Stratified { profile } => {
            writer.byte(2);
            encode_source_ref(writer, profile)?;
        }
    }
    encode_source_ref(writer, &slice.coverage)
}

fn decode_gauge_slice(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeSlicePlan, IdentifiabilityError> {
    let support = decode_role_set(reader, "gauge slice support")?;
    let constraint = decode_source_ref(reader)?;
    let expected_codimension = match reader.byte("gauge slice codimension")? {
        0 => GaugeSliceCodimension::FixedFinite {
            codimension: reader.u64("finite gauge slice codimension")?,
        },
        1 => GaugeSliceCodimension::FixedInfinite {
            codimension_model: decode_source_ref(reader)?,
            compatibility: decode_source_ref(reader)?,
        },
        2 => GaugeSliceCodimension::Stratified {
            profile: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-slice-codimension tag {tag}"),
            });
        }
    };
    GaugeSlicePlan::try_new(
        support,
        constraint,
        expected_codimension,
        decode_source_ref(reader)?,
    )
}

fn encode_gauge_quotient(
    writer: &mut CanonicalWriter,
    quotient: &GaugeQuotientPlan,
) -> Result<(), IdentifiabilityError> {
    match quotient {
        GaugeQuotientPlan::RegularAtlas {
            quotient_map,
            local_section_atlas,
            coverage,
        } => {
            writer.byte(0);
            encode_source_ref(writer, quotient_map)?;
            encode_source_ref(writer, local_section_atlas)?;
            encode_source_ref(writer, coverage)?;
        }
        GaugeQuotientPlan::SingularOrGeneralized {
            quotient_map,
            quotient_profile,
            local_models,
        } => {
            writer.byte(1);
            encode_source_ref(writer, quotient_map)?;
            encode_source_ref(writer, quotient_profile)?;
            writer.byte(u8::from(local_models.is_some()));
            if let Some(local_models) = local_models {
                encode_source_ref(writer, local_models)?;
            }
        }
        GaugeQuotientPlan::InvariantMap {
            invariants,
            completeness_profile,
        } => {
            writer.byte(2);
            encode_source_ref(writer, invariants)?;
            encode_source_ref(writer, completeness_profile)?;
        }
        GaugeQuotientPlan::GroupoidOrStack {
            presentation,
            quotient_profile,
        } => {
            writer.byte(3);
            encode_source_ref(writer, presentation)?;
            encode_source_ref(writer, quotient_profile)?;
        }
    }
    Ok(())
}

fn decode_gauge_quotient(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeQuotientPlan, IdentifiabilityError> {
    match reader.byte("gauge quotient plan")? {
        0 => Ok(GaugeQuotientPlan::RegularAtlas {
            quotient_map: decode_source_ref(reader)?,
            local_section_atlas: decode_source_ref(reader)?,
            coverage: decode_source_ref(reader)?,
        }),
        1 => {
            let quotient_map = decode_source_ref(reader)?;
            let quotient_profile = decode_source_ref(reader)?;
            let local_models = match reader.byte("gauge local models presence")? {
                0 => None,
                1 => Some(decode_source_ref(reader)?),
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown gauge-local-model-presence tag {tag}"),
                    });
                }
            };
            Ok(GaugeQuotientPlan::SingularOrGeneralized {
                quotient_map,
                quotient_profile,
                local_models,
            })
        }
        2 => Ok(GaugeQuotientPlan::InvariantMap {
            invariants: decode_source_ref(reader)?,
            completeness_profile: decode_source_ref(reader)?,
        }),
        3 => Ok(GaugeQuotientPlan::GroupoidOrStack {
            presentation: decode_source_ref(reader)?,
            quotient_profile: decode_source_ref(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-quotient-plan tag {tag}"),
        }),
    }
}

fn encode_continuous_gauge_reduction(
    writer: &mut CanonicalWriter,
    reduction: &ContinuousGaugeReductionPlan,
) -> Result<(), IdentifiabilityError> {
    match reduction {
        ContinuousGaugeReductionPlan::Quotient { quotient } => {
            writer.byte(0);
            encode_gauge_quotient(writer, quotient)?;
        }
        ContinuousGaugeReductionPlan::Slice { slice } => {
            writer.byte(1);
            encode_gauge_slice(writer, slice)?;
        }
    }
    Ok(())
}

fn decode_continuous_gauge_reduction(
    reader: &mut CanonicalReader<'_>,
) -> Result<ContinuousGaugeReductionPlan, IdentifiabilityError> {
    match reader.byte("continuous gauge reduction")? {
        0 => Ok(ContinuousGaugeReductionPlan::Quotient {
            quotient: decode_gauge_quotient(reader)?,
        }),
        1 => Ok(ContinuousGaugeReductionPlan::Slice {
            slice: decode_gauge_slice(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown continuous-gauge-reduction tag {tag}"),
        }),
    }
}

fn encode_gauge_reduction_stage(
    writer: &mut CanonicalWriter,
    stage: &GaugeReductionStage,
) -> Result<(), IdentifiabilityError> {
    match stage {
        GaugeReductionStage::Root => writer.byte(0),
        GaugeReductionStage::After {
            predecessors,
            composition_law,
            relation,
        } => {
            writer.byte(1);
            writer.count(predecessors.len(), "gauge reduction stage predecessors")?;
            for predecessor in predecessors {
                writer.text(predecessor.as_str(), "gauge reduction id")?;
            }
            encode_source_ref(writer, composition_law)?;
            match relation {
                GaugeReductionStageRelation::NormalSubgroupTower {
                    normality,
                    induced_residual_action,
                } => {
                    writer.byte(0);
                    encode_source_ref(writer, normality)?;
                    encode_source_ref(writer, induced_residual_action)?;
                }
                GaugeReductionStageRelation::SemidirectOrGenerated {
                    extension,
                    induced_residual_action,
                } => {
                    writer.byte(1);
                    encode_source_ref(writer, extension)?;
                    encode_source_ref(writer, induced_residual_action)?;
                }
                GaugeReductionStageRelation::TransverseSlices { transversality } => {
                    writer.byte(2);
                    encode_source_ref(writer, transversality)?;
                }
                GaugeReductionStageRelation::GaugeForGauge {
                    reducibility,
                    induced_residual_action,
                } => {
                    writer.byte(3);
                    encode_source_ref(writer, reducibility)?;
                    encode_source_ref(writer, induced_residual_action)?;
                }
            }
        }
    }
    Ok(())
}

fn decode_gauge_reduction_stage(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeReductionStage, IdentifiabilityError> {
    match reader.byte("gauge reduction stage")? {
        0 => Ok(GaugeReductionStage::Root),
        1 => {
            let count = reader.count("gauge reduction stage predecessors")?;
            let mut predecessors = BTreeSet::new();
            for _ in 0..count {
                let predecessor = GaugeReductionId::try_new(reader.token("gauge reduction id")?)?;
                if !predecessors.insert(predecessor.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "gauge reduction stage predecessor",
                        id: predecessor.to_string(),
                    });
                }
            }
            let composition_law = decode_source_ref(reader)?;
            let relation = match reader.byte("gauge reduction stage relation")? {
                0 => GaugeReductionStageRelation::NormalSubgroupTower {
                    normality: decode_source_ref(reader)?,
                    induced_residual_action: decode_source_ref(reader)?,
                },
                1 => GaugeReductionStageRelation::SemidirectOrGenerated {
                    extension: decode_source_ref(reader)?,
                    induced_residual_action: decode_source_ref(reader)?,
                },
                2 => GaugeReductionStageRelation::TransverseSlices {
                    transversality: decode_source_ref(reader)?,
                },
                3 => GaugeReductionStageRelation::GaugeForGauge {
                    reducibility: decode_source_ref(reader)?,
                    induced_residual_action: decode_source_ref(reader)?,
                },
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown gauge-reduction-stage-relation tag {tag}"),
                    });
                }
            };
            Ok(GaugeReductionStage::After {
                predecessors,
                composition_law,
                relation,
            })
        }
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-reduction-stage tag {tag}"),
        }),
    }
}

fn encode_gauge_measure_semantics(
    writer: &mut CanonicalWriter,
    measure: &GaugeMeasureSemantics,
) -> Result<(), IdentifiabilityError> {
    match measure {
        GaugeMeasureSemantics::NotApplicable { reason } => {
            writer.byte(0);
            writer.text(reason, "gauge reduction measure semantics")?;
        }
        GaugeMeasureSemantics::Pushforward {
            source_measure,
            reduced_measure,
            transport,
            jacobian_or_disintegration,
        } => {
            writer.byte(1);
            encode_source_ref(writer, source_measure)?;
            encode_source_ref(writer, reduced_measure)?;
            encode_source_ref(writer, transport)?;
            encode_source_ref(writer, jacobian_or_disintegration)?;
        }
    }
    Ok(())
}

fn decode_gauge_measure_semantics(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeMeasureSemantics, IdentifiabilityError> {
    match reader.byte("gauge reduction measure semantics")? {
        0 => Ok(GaugeMeasureSemantics::NotApplicable {
            reason: reader.reason("gauge reduction measure semantics")?,
        }),
        1 => Ok(GaugeMeasureSemantics::Pushforward {
            source_measure: decode_source_ref(reader)?,
            reduced_measure: decode_source_ref(reader)?,
            transport: decode_source_ref(reader)?,
            jacobian_or_disintegration: decode_source_ref(reader)?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown gauge-reduction-measure-semantics tag {tag}"),
        }),
    }
}

fn encode_gauge_reduction_binding(
    writer: &mut CanonicalWriter,
    binding: &GaugeReductionBinding,
) -> Result<(), IdentifiabilityError> {
    writer.text(binding.id.as_str(), "gauge reduction id")?;
    encode_gauge_action_reference(writer, &binding.action)?;
    writer.count(binding.claims.len(), "gauge reduction claims")?;
    for claim in &binding.claims {
        writer.text(claim.as_str(), "claim id")?;
    }
    match &binding.plan {
        GaugeReductionPlan::Unreduced { reason } => {
            writer.byte(0);
            writer.text(reason, "unreduced gauge reason")?;
        }
        GaugeReductionPlan::Quotient { quotient } => {
            writer.byte(1);
            encode_gauge_quotient(writer, quotient)?;
        }
        GaugeReductionPlan::Slice { slice } => {
            writer.byte(2);
            encode_gauge_slice(writer, slice)?;
        }
        GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction,
            normal_subgroup,
            factor_extension,
            residual_quotient_action,
            compatibility,
        } => {
            writer.byte(3);
            encode_continuous_gauge_reduction(writer, reduction)?;
            encode_source_ref(writer, normal_subgroup)?;
            encode_source_ref(writer, factor_extension)?;
            encode_source_ref(writer, residual_quotient_action)?;
            encode_source_ref(writer, compatibility)?;
        }
    }
    encode_gauge_reduction_stage(writer, &binding.stage)?;
    encode_gauge_measure_semantics(writer, &binding.measure)?;
    Ok(())
}

fn decode_gauge_reduction_binding(
    reader: &mut CanonicalReader<'_>,
) -> Result<GaugeReductionBinding, IdentifiabilityError> {
    let id = GaugeReductionId::try_new(reader.token("gauge reduction id")?)?;
    let action = decode_gauge_action_reference(reader)?;
    let count = reader.count("gauge reduction claims")?;
    let mut claims = BTreeSet::new();
    for _ in 0..count {
        let claim = ClaimId::try_new(reader.token("claim id")?)?;
        if !claims.insert(claim.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "gauge reduction claim",
                id: claim.to_string(),
            });
        }
    }
    let plan = match reader.byte("gauge reduction plan")? {
        0 => GaugeReductionPlan::Unreduced {
            reason: reader.reason("unreduced gauge reason")?,
        },
        1 => GaugeReductionPlan::Quotient {
            quotient: decode_gauge_quotient(reader)?,
        },
        2 => GaugeReductionPlan::Slice {
            slice: decode_gauge_slice(reader)?,
        },
        3 => GaugeReductionPlan::ContinuousReductionWithDiscreteResidual {
            reduction: decode_continuous_gauge_reduction(reader)?,
            normal_subgroup: decode_source_ref(reader)?,
            factor_extension: decode_source_ref(reader)?,
            residual_quotient_action: decode_source_ref(reader)?,
            compatibility: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown gauge-reduction-plan tag {tag}"),
            });
        }
    };
    let stage = decode_gauge_reduction_stage(reader)?;
    let measure = decode_gauge_measure_semantics(reader)?;
    GaugeReductionBinding::try_new(id, action, claims, plan, stage, measure)
}

fn encode_joint_noise(
    writer: &mut CanonicalWriter,
    noise: &JointNoiseModel,
) -> Result<(), IdentifiabilityError> {
    match noise {
        JointNoiseModel::Independent { assumption } => {
            writer.byte(0);
            encode_source_key(writer, assumption)?;
        }
        JointNoiseModel::DenseCorrelation {
            order,
            correlation,
            model,
        } => {
            writer.byte(1);
            writer.count(order.len(), "joint-noise order")?;
            for key in order {
                encode_observation_key(writer, key)?;
            }
            writer.count(correlation.dimension(), "correlation dimension")?;
            writer.count(correlation.lower_triangle().len(), "correlation entries")?;
            for value in correlation.lower_triangle() {
                writer.f64(*value);
            }
            encode_source_key(writer, model)?;
        }
        JointNoiseModel::ExternalKernel { model } => {
            writer.byte(2);
            encode_source_key(writer, model)?;
        }
        JointNoiseModel::Unknown { reason } => {
            writer.byte(3);
            writer.text(reason, "unknown joint noise reason")?;
        }
    }
    Ok(())
}

fn decode_joint_noise(
    reader: &mut CanonicalReader<'_>,
) -> Result<JointNoiseModel, IdentifiabilityError> {
    match reader.byte("joint noise")? {
        0 => Ok(JointNoiseModel::Independent {
            assumption: decode_source_key(reader)?,
        }),
        1 => {
            let count = reader.count("joint-noise order")?;
            preflight_collection_bytes(reader, count, 8, "joint-noise order")?;
            let mut order = Vec::with_capacity(count);
            for _ in 0..count {
                order.push(decode_observation_key(reader)?);
            }
            let dimension = reader.count("correlation dimension")?;
            if dimension == 0 || dimension > MAX_VV_MATRIX_DIMENSION {
                return Err(IdentifiabilityError::Cardinality {
                    field: "correlation dimension",
                    detail: format!(
                        "matrix dimension {dimension} lies outside 1..={MAX_VV_MATRIX_DIMENSION}"
                    ),
                });
            }
            let expected_entries = dimension
                .checked_mul(dimension.saturating_add(1))
                .and_then(|value| value.checked_div(2))
                .ok_or_else(|| IdentifiabilityError::Cardinality {
                    field: "correlation entries",
                    detail: "matrix entry count overflows address space".to_string(),
                })?;
            let encoded_entries =
                usize::try_from(reader.u32("correlation entries")?).map_err(|_| {
                    IdentifiabilityError::Cardinality {
                        field: "correlation entries",
                        detail: "matrix entry count exceeds address space".to_string(),
                    }
                })?;
            if encoded_entries != expected_entries {
                return Err(IdentifiabilityError::Cardinality {
                    field: "correlation entries",
                    detail: format!(
                        "dimension {dimension} requires exactly {expected_entries} entries, found {encoded_entries}"
                    ),
                });
            }
            let required_bytes = expected_entries
                .checked_mul(core::mem::size_of::<f64>())
                .ok_or_else(|| IdentifiabilityError::Cardinality {
                    field: "correlation entries",
                    detail: "matrix byte count overflows address space".to_string(),
                })?;
            if reader.bytes.len().saturating_sub(reader.at) < required_bytes {
                return Err(IdentifiabilityError::Canonical {
                    at: reader.at,
                    detail: format!(
                        "truncated correlation entries: need {required_bytes} bytes before trailing fields"
                    ),
                });
            }
            let mut lower = Vec::with_capacity(expected_entries);
            for _ in 0..expected_entries {
                lower.push(reader.f64("correlation entry")?);
            }
            let correlation = CovarianceMatrix::try_new(dimension, lower).map_err(|error| {
                IdentifiabilityError::Vv {
                    detail: error.to_string(),
                }
            })?;
            Ok(JointNoiseModel::DenseCorrelation {
                order,
                correlation,
                model: decode_source_key(reader)?,
            })
        }
        2 => Ok(JointNoiseModel::ExternalKernel {
            model: decode_source_key(reader)?,
        }),
        3 => Ok(JointNoiseModel::Unknown {
            reason: reader.reason("unknown joint noise reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown joint-noise tag {tag}"),
        }),
    }
}

fn encode_data_reuse(
    writer: &mut CanonicalWriter,
    policy: &DataReusePolicy,
) -> Result<(), IdentifiabilityError> {
    match policy {
        DataReusePolicy::Disjoint => writer.byte(0),
        DataReusePolicy::Shared { groups } => {
            writer.byte(1);
            writer.count(groups.len(), "data sharing groups")?;
            for group in groups {
                writer.count(group.cases.len(), "sharing-group cases")?;
                for case in &group.cases {
                    encode_case_id(writer, case)?;
                }
                encode_source_key(writer, &group.joint_likelihood)?;
                writer.text(&group.justification, "data sharing justification")?;
            }
        }
    }
    Ok(())
}

fn decode_data_reuse(
    reader: &mut CanonicalReader<'_>,
) -> Result<DataReusePolicy, IdentifiabilityError> {
    match reader.byte("data reuse policy")? {
        0 => Ok(DataReusePolicy::Disjoint),
        1 => {
            let count = reader.count("data sharing groups")?;
            preflight_collection_bytes(reader, count, 1, "data sharing groups")?;
            let mut groups = Vec::with_capacity(count);
            for _ in 0..count {
                let case_count = reader.count("sharing-group cases")?;
                let mut cases = BTreeSet::new();
                for _ in 0..case_count {
                    let case = decode_case_id(reader)?;
                    if !cases.insert(case.clone()) {
                        return Err(IdentifiabilityError::Duplicate {
                            field: "sharing-group case",
                            id: case.to_string(),
                        });
                    }
                }
                groups.push(DataSharingGroup::try_new(
                    cases,
                    decode_source_key(reader)?,
                    reader.reason("data sharing justification")?,
                )?);
            }
            Ok(DataReusePolicy::Shared { groups })
        }
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown data-reuse tag {tag}"),
        }),
    }
}

fn encode_problem(
    document: &IdentifiabilityProblemDocument,
) -> Result<Vec<u8>, IdentifiabilityError> {
    check_problem_identity_version(document.schema_version)?;
    let mut writer = CanonicalWriter::new();
    writer.raw(PROBLEM_MAGIC);
    writer.u32(document.schema_version);
    encode_source_key(&mut writer, &document.context_source)?;
    encode_source_key(&mut writer, &document.material_source)?;
    encode_source_key(&mut writer, &document.model_source)?;
    encode_source_key(&mut writer, &document.graph_source)?;
    match &document.joint_prior {
        Some(joint_prior) => {
            writer.byte(1);
            encode_source_key(&mut writer, joint_prior)?;
        }
        None => writer.byte(0),
    }
    writer.count(document.sources.len(), "source registry")?;
    for source in document.sources.values() {
        encode_source_ref(&mut writer, source)?;
    }
    writer.count(document.parameters.len(), "study parameters")?;
    for parameter in document.parameters.values() {
        encode_study_parameter(&mut writer, parameter)?;
    }
    writer.count(document.constraints.len(), "joint constraints")?;
    for constraint in document.constraints.values() {
        encode_constraint(&mut writer, constraint)?;
    }
    encode_admissible_domain_witness(&mut writer, &document.admissible_domain)?;
    writer.count(document.cases.len(), "study cases")?;
    for case in document.cases.values() {
        encode_case(&mut writer, case)?;
    }
    writer.count(document.influences.len(), "influence declarations")?;
    for influence in document.influences.values() {
        encode_influence(&mut writer, influence)?;
    }
    writer.count(document.gauges.len(), "gauge declarations")?;
    for gauge in document.gauges.values() {
        encode_gauge(&mut writer, gauge)?;
    }
    writer.count(
        document.gauge_compositions.len(),
        "gauge composition declarations",
    )?;
    for composition in document.gauge_compositions.values() {
        encode_gauge_composition(&mut writer, composition)?;
    }
    encode_joint_noise(&mut writer, &document.joint_noise)?;
    encode_data_reuse(&mut writer, &document.data_reuse)?;
    writer.finish()
}

fn decode_problem(bytes: &[u8]) -> Result<IdentifiabilityProblemDocument, IdentifiabilityError> {
    let mut reader = CanonicalReader::new(bytes)?;
    let magic = reader.take(PROBLEM_MAGIC.len(), "problem magic")?;
    if magic != PROBLEM_MAGIC {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "wrong identifiability-problem magic".to_string(),
        });
    }
    let version = reader.u32("problem schema version")?;
    check_problem_identity_version(version)?;
    let context_source = decode_source_key(&mut reader)?;
    let material_source = decode_source_key(&mut reader)?;
    let model_source = decode_source_key(&mut reader)?;
    let graph_source = decode_source_key(&mut reader)?;
    let joint_prior = match reader.byte("problem joint-prior presence")? {
        0 => None,
        1 => Some(decode_source_key(&mut reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown problem joint-prior presence tag {tag}"),
            });
        }
    };
    let source_count = reader.count("source registry")?;
    preflight_collection_bytes(&reader, source_count, 45, "source registry")?;
    let mut sources = Vec::with_capacity(source_count);
    for _ in 0..source_count {
        sources.push(decode_source_ref(&mut reader)?);
    }
    let parameter_count = reader.count("study parameters")?;
    preflight_collection_bytes(&reader, parameter_count, 1, "study parameters")?;
    let mut parameters = Vec::with_capacity(parameter_count);
    for _ in 0..parameter_count {
        parameters.push(decode_study_parameter(&mut reader)?);
    }
    let constraint_count = reader.count("joint constraints")?;
    preflight_collection_bytes(&reader, constraint_count, 1, "joint constraints")?;
    let mut constraints = Vec::with_capacity(constraint_count);
    for _ in 0..constraint_count {
        constraints.push(decode_constraint(&mut reader)?);
    }
    let admissible_domain = decode_admissible_domain_witness(&mut reader)?;
    let case_count = reader.count("study cases")?;
    preflight_collection_bytes(&reader, case_count, 1, "study cases")?;
    let mut cases = Vec::with_capacity(case_count);
    for _ in 0..case_count {
        cases.push(decode_case(&mut reader)?);
    }
    let influence_count = reader.count("influence declarations")?;
    preflight_collection_bytes(&reader, influence_count, 1, "influence declarations")?;
    let mut influences = Vec::with_capacity(influence_count);
    for _ in 0..influence_count {
        influences.push(decode_influence(&mut reader)?);
    }
    let gauge_count = reader.count("gauge declarations")?;
    preflight_collection_bytes(&reader, gauge_count, 1, "gauge declarations")?;
    let mut gauges = Vec::with_capacity(gauge_count);
    for _ in 0..gauge_count {
        gauges.push(decode_gauge(&mut reader)?);
    }
    let composition_count = reader.count("gauge composition declarations")?;
    preflight_collection_bytes(
        &reader,
        composition_count,
        1,
        "gauge composition declarations",
    )?;
    let mut gauge_compositions = Vec::with_capacity(composition_count);
    for _ in 0..composition_count {
        gauge_compositions.push(decode_gauge_composition(&mut reader)?);
    }
    let joint_noise = decode_joint_noise(&mut reader)?;
    let data_reuse = decode_data_reuse(&mut reader)?;
    reader.finish()?;
    let document = IdentifiabilityProblemDocument::try_new(
        context_source,
        material_source,
        model_source,
        graph_source,
        joint_prior,
        sources,
        parameters,
        constraints,
        admissible_domain,
        cases,
        influences,
        gauges,
        gauge_compositions,
        joint_noise,
        data_reuse,
    )?;
    if document.canonical_bytes()?.as_slice() != bytes {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "non-canonical problem encoding".to_string(),
        });
    }
    Ok(document)
}

fn check_identity_version(declared: u32, supported: u32) -> Result<(), IdentifiabilityError> {
    if declared == supported {
        Ok(())
    } else {
        Err(IdentifiabilityError::UnsupportedSchemaVersion {
            declared,
            supported,
        })
    }
}

/// Fail closed on a stale/future umbrella API generation. Identity transports
/// must use their stage-specific checkers below.
pub fn check_authority_schema_version(declared: u32) -> Result<(), IdentifiabilityError> {
    check_identity_version(declared, IDENTIFIABILITY_AUTHORITY_SCHEMA_VERSION)
}

/// Fail closed on a stale/future physical-problem identity version.
pub fn check_problem_identity_version(declared: u32) -> Result<(), IdentifiabilityError> {
    check_identity_version(declared, IDENTIFIABILITY_PROBLEM_IDENTITY_VERSION)
}

/// Fail closed on a stale/future source-admission identity version.
pub fn check_source_admission_identity_version(declared: u32) -> Result<(), IdentifiabilityError> {
    check_identity_version(declared, IDENTIFIABILITY_SOURCE_ADMISSION_IDENTITY_VERSION)
}

/// Fail closed on a stale/future execution identity version.
pub fn check_execution_identity_version(declared: u32) -> Result<(), IdentifiabilityError> {
    check_identity_version(declared, IDENTIFIABILITY_EXECUTION_IDENTITY_VERSION)
}

/// Fail closed on a stale/future assessment identity version.
pub fn check_assessment_identity_version(declared: u32) -> Result<(), IdentifiabilityError> {
    check_identity_version(declared, IDENTIFIABILITY_ASSESSMENT_IDENTITY_VERSION)
}

fn encode_execution_action(
    writer: &mut CanonicalWriter,
    action: &ParameterExecutionAction,
) -> Result<(), IdentifiabilityError> {
    match action {
        ParameterExecutionAction::Optimize { coordinate } => {
            writer.byte(0);
            encode_coordinate(writer, coordinate)?;
        }
        ParameterExecutionAction::Profile { coordinate } => {
            writer.byte(1);
            encode_coordinate(writer, coordinate)?;
        }
        ParameterExecutionAction::Marginalize {
            coordinate,
            integrator,
            measure_transport,
        } => {
            writer.byte(2);
            encode_coordinate(writer, coordinate)?;
            encode_source_ref(writer, integrator)?;
            encode_source_ref(writer, measure_transport)?;
        }
        ParameterExecutionAction::Conditioned => writer.byte(3),
        ParameterExecutionAction::Derived => writer.byte(4),
    }
    Ok(())
}

fn decode_execution_action(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterExecutionAction, IdentifiabilityError> {
    match reader.byte("execution action")? {
        0 => Ok(ParameterExecutionAction::Optimize {
            coordinate: decode_coordinate(reader)?,
        }),
        1 => Ok(ParameterExecutionAction::Profile {
            coordinate: decode_coordinate(reader)?,
        }),
        2 => Ok(ParameterExecutionAction::Marginalize {
            coordinate: decode_coordinate(reader)?,
            integrator: decode_source_ref(reader)?,
            measure_transport: decode_source_ref(reader)?,
        }),
        3 => Ok(ParameterExecutionAction::Conditioned),
        4 => Ok(ParameterExecutionAction::Derived),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown execution-action tag {tag}"),
        }),
    }
}

fn encode_dimensionless_error_policy(
    writer: &mut CanonicalWriter,
    policy: &DimensionlessErrorPolicy,
) -> Result<(), IdentifiabilityError> {
    encode_source_ref(writer, &policy.metric)?;
    encode_source_ref(writer, &policy.nondimensionalization)?;
    writer.f64(policy.maximum_certified_error);
    Ok(())
}

fn decode_dimensionless_error_policy(
    reader: &mut CanonicalReader<'_>,
) -> Result<DimensionlessErrorPolicy, IdentifiabilityError> {
    DimensionlessErrorPolicy::try_new(
        decode_source_ref(reader)?,
        decode_source_ref(reader)?,
        reader.f64("maximum certified claim error")?,
    )
}

fn encode_claim_request(
    writer: &mut CanonicalWriter,
    request: &ClaimRequest,
) -> Result<(), IdentifiabilityError> {
    encode_claim(writer, &request.claim)?;
    encode_dimensionless_error_policy(writer, &request.error_policy)
}

fn decode_claim_request(
    reader: &mut CanonicalReader<'_>,
) -> Result<ClaimRequest, IdentifiabilityError> {
    Ok(ClaimRequest::new(
        decode_claim(reader)?,
        decode_dimensionless_error_policy(reader)?,
    ))
}

fn encode_execution_with_header_mode(
    plan: &IdentifiabilityExecutionPlan,
    exact_header: bool,
) -> Result<Vec<u8>, IdentifiabilityError> {
    check_execution_identity_version(plan.schema_version)?;
    let mut writer = CanonicalWriter::new();
    writer.raw(EXECUTION_MAGIC);
    writer.u32(plan.schema_version);
    encode_header(&mut writer, &plan.header, exact_header)?;
    writer.hash(plan.problem_id.0);
    writer.hash(plan.source_admission_id.0);
    encode_source_ref(&mut writer, &plan.analyzer)?;
    encode_source_ref(&mut writer, &plan.build)?;
    match &plan.derivative_provider {
        Some(source) => {
            writer.byte(1);
            encode_source_ref(&mut writer, source)?;
        }
        None => writer.byte(0),
    }
    writer.count(plan.claim_requests.len(), "claim requests")?;
    for request in plan.claim_requests.values() {
        encode_claim_request(&mut writer, request)?;
    }
    writer.count(plan.actions.len(), "parameter actions")?;
    for (role, action) in &plan.actions {
        encode_role(&mut writer, role)?;
        encode_execution_action(&mut writer, action)?;
    }
    writer.count(plan.gauge_reductions.len(), "gauge reductions")?;
    for binding in plan.gauge_reductions.values() {
        encode_gauge_reduction_binding(&mut writer, binding)?;
    }
    writer.f64(plan.numerical.rank_tolerance);
    writer.f64(plan.numerical.singular_value_floor);
    writer.f64(plan.numerical.maximum_condition_number);
    writer.byte(match plan.numerical.arithmetic {
        ArithmeticPolicy::ExactSymbolic => 0,
        ArithmeticPolicy::CertifiedInterval => 1,
        ArithmeticPolicy::DeterministicFloatingPoint => 2,
        ArithmeticPolicy::FastFloatingPoint => 3,
    });
    encode_source_ref(&mut writer, &plan.numerical.nondimensionalization)?;
    encode_source_ref(&mut writer, &plan.initialization)?;
    encode_source_ref(&mut writer, &plan.stopping)?;
    encode_source_ref(&mut writer, &plan.determinism_contract)?;
    encode_resolution_set(&mut writer, &plan.source_authority)?;
    writer.finish()
}

fn encode_execution(plan: &IdentifiabilityExecutionPlan) -> Result<Vec<u8>, IdentifiabilityError> {
    encode_execution_with_header_mode(plan, true)
}

fn encode_execution_identity(
    plan: &IdentifiabilityExecutionPlan,
) -> Result<Vec<u8>, IdentifiabilityError> {
    encode_execution_with_header_mode(plan, false)
}

fn execution_identity_hash(
    plan: &IdentifiabilityExecutionPlan,
) -> Result<ExecutionId, IdentifiabilityError> {
    Ok(ExecutionId(hash_domain(
        IDENTIFIABILITY_EXECUTION_IDENTITY_DOMAIN,
        &encode_execution_identity(plan)?,
    )))
}

fn decode_execution(
    bytes: &[u8],
    problem: &AdmittedIdentifiabilityProblem,
    verified_sources: &SourceResolutionSet,
) -> Result<IdentifiabilityExecutionPlan, IdentifiabilityError> {
    let mut reader = CanonicalReader::new(bytes)?;
    if reader.take(EXECUTION_MAGIC.len(), "execution magic")? != EXECUTION_MAGIC {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "wrong identifiability-execution magic".to_string(),
        });
    }
    check_execution_identity_version(reader.u32("execution schema version")?)?;
    let header = decode_header(&mut reader)?;
    let problem_id = ProblemId(reader.hash("execution problem id")?);
    let source_admission_id = SourceAdmissionId(reader.hash("execution source-admission id")?);
    if problem_id != problem.problem_id || source_admission_id != problem.source_admission_id {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "execution problem/source admission",
        });
    }
    let analyzer = decode_source_ref(&mut reader)?;
    let build = decode_source_ref(&mut reader)?;
    let derivative_provider = match reader.byte("derivative-provider option")? {
        0 => None,
        1 => Some(decode_source_ref(&mut reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("invalid derivative-provider option tag {tag}"),
            });
        }
    };
    let request_count = reader.count("claim requests")?;
    preflight_collection_bytes(&reader, request_count, 1, "claim requests")?;
    let mut claim_requests = Vec::with_capacity(request_count);
    for _ in 0..request_count {
        claim_requests.push(decode_claim_request(&mut reader)?);
    }
    let action_count = reader.count("parameter actions")?;
    preflight_collection_bytes(&reader, action_count, 1, "parameter actions")?;
    let mut actions = Vec::with_capacity(action_count);
    for _ in 0..action_count {
        actions.push((
            decode_role(&mut reader)?,
            decode_execution_action(&mut reader)?,
        ));
    }
    let gauge_reduction_count = reader.count("gauge reductions")?;
    preflight_collection_bytes(&reader, gauge_reduction_count, 1, "gauge reductions")?;
    let mut gauge_reductions = Vec::with_capacity(gauge_reduction_count);
    for _ in 0..gauge_reduction_count {
        gauge_reductions.push(decode_gauge_reduction_binding(&mut reader)?);
    }
    let rank_tolerance = reader.f64("rank tolerance")?;
    let singular_value_floor = reader.f64("singular-value floor")?;
    let maximum_condition_number = reader.f64("maximum condition number")?;
    let arithmetic = match reader.byte("arithmetic policy")? {
        0 => ArithmeticPolicy::ExactSymbolic,
        1 => ArithmeticPolicy::CertifiedInterval,
        2 => ArithmeticPolicy::DeterministicFloatingPoint,
        3 => ArithmeticPolicy::FastFloatingPoint,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown arithmetic-policy tag {tag}"),
            });
        }
    };
    let numerical_nondimensionalization = decode_source_ref(&mut reader)?;
    let numerical = IdentifiabilityNumericalPolicy::try_new(
        rank_tolerance,
        singular_value_floor,
        maximum_condition_number,
        arithmetic,
        numerical_nondimensionalization,
    )?;
    let initialization = decode_source_ref(&mut reader)?;
    let stopping = decode_source_ref(&mut reader)?;
    let determinism_contract = decode_source_ref(&mut reader)?;
    let transported_sources = decode_resolution_set(&mut reader)?;
    reader.finish()?;
    if transported_sources != *verified_sources {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "execution source-resolution replay",
        });
    }
    let plan = IdentifiabilityExecutionPlan::try_new(
        header,
        problem,
        analyzer,
        build,
        derivative_provider,
        claim_requests,
        actions,
        gauge_reductions,
        numerical,
        initialization,
        stopping,
        determinism_contract,
        verified_sources.clone(),
    )?;
    if plan.canonical_bytes()?.as_slice() != bytes {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "non-canonical execution encoding".to_string(),
        });
    }
    Ok(plan)
}

fn encode_claim(
    writer: &mut CanonicalWriter,
    claim: &TypedIdentifiabilityClaim,
) -> Result<(), IdentifiabilityError> {
    writer.text(claim.id.as_str(), "claim id")?;
    match &claim.information {
        InformationRegime::StructuralExactModel => writer.byte(0),
        InformationRegime::ExactInputOutputMap => writer.byte(1),
        InformationRegime::NoisyFiniteData => writer.byte(2),
        InformationRegime::PosteriorUnderDeclaredPrior { joint_prior } => {
            writer.byte(3);
            encode_source_ref(writer, joint_prior)?;
        }
    }
    writer.byte(match claim.extent {
        IdentifiabilityExtent::Local => 0,
        IdentifiabilityExtent::Global => 1,
    });
    match &claim.fiber {
        FiberStructure::Unique => writer.byte(0),
        FiberStructure::FiniteToOne {
            maximum_cardinality,
        } => {
            writer.byte(1);
            match maximum_cardinality {
                Some(FiberCardinalityBound::UniformU64(maximum)) => {
                    writer.byte(1);
                    writer.u64(*maximum);
                }
                Some(FiberCardinalityBound::SymbolicProfile(profile)) => {
                    writer.byte(2);
                    encode_source_ref(writer, profile)?;
                }
                None => writer.byte(0),
            }
        }
        FiberStructure::DiscreteOrbit { action } => {
            writer.byte(2);
            encode_gauge_action_reference(writer, action)?;
        }
        FiberStructure::OrbitQuotientUnique { action } => {
            writer.byte(3);
            encode_gauge_action_reference(writer, action)?;
        }
        FiberStructure::PositiveDimensional { lower_bound } => {
            writer.byte(4);
            match lower_bound {
                FiberDimensionLowerBound::Finite { minimum_dimension } => {
                    writer.byte(0);
                    writer.u64(*minimum_dimension);
                }
                FiberDimensionLowerBound::InfiniteDimensional { model_space } => {
                    writer.byte(1);
                    encode_source_ref(writer, model_space)?;
                }
            }
        }
        FiberStructure::Stratified { strata } => {
            writer.byte(5);
            encode_source_ref(writer, strata)?;
        }
        FiberStructure::MixedOrbit { action } => {
            writer.byte(6);
            encode_gauge_action_reference(writer, action)?;
        }
        FiberStructure::GeneralizedQuotientUnique {
            action,
            equivalence,
        } => {
            writer.byte(7);
            encode_gauge_action_reference(writer, action)?;
            encode_source_ref(writer, equivalence)?;
        }
        FiberStructure::StratifiedOrbit {
            action,
            orbit_type_profile,
        } => {
            writer.byte(8);
            encode_gauge_action_reference(writer, action)?;
            encode_source_ref(writer, orbit_type_profile)?;
        }
    }
    match &claim.quantifier {
        ClaimQuantifier::AtRealization { realization } => {
            writer.byte(0);
            encode_source_ref(writer, realization)?;
        }
        ClaimQuantifier::AlmostEverywhere { measure } => {
            writer.byte(1);
            encode_source_ref(writer, measure)?;
        }
        ClaimQuantifier::ForAll { domain } => {
            writer.byte(2);
            encode_source_ref(writer, domain)?;
        }
        ClaimQuantifier::ProbabilityAtLeast {
            probability,
            measure,
        } => {
            writer.byte(3);
            writer.f64(*probability);
            encode_source_ref(writer, measure)?;
        }
    }
    match &claim.scalar_domain {
        ScalarDomain::Real => writer.byte(0),
        ScalarDomain::Complex { extension } => {
            writer.byte(1);
            encode_source_ref(writer, extension)?;
        }
        ScalarDomain::MixedDiscreteContinuous { stratification } => {
            writer.byte(2);
            encode_source_ref(writer, stratification)?;
        }
    }
    match &claim.subject {
        ClaimSubject::Parameter(role) => {
            writer.byte(0);
            encode_role(writer, role)?;
        }
        ClaimSubject::ParameterSet(roles) => {
            writer.byte(1);
            writer.count(roles.len(), "claim parameter set")?;
            for role in roles {
                encode_role(writer, role)?;
            }
        }
        ClaimSubject::DerivedFunctional {
            definition,
            parameters,
        } => {
            writer.byte(2);
            encode_source_ref(writer, definition)?;
            writer.count(parameters.len(), "derived-functional parameters")?;
            for role in parameters {
                encode_role(writer, role)?;
            }
        }
        ClaimSubject::Influence(influence) => {
            writer.byte(3);
            writer.text(influence.as_str(), "influence id")?;
        }
        ClaimSubject::GaugeAction(action) => {
            writer.byte(4);
            encode_gauge_action_reference(writer, action)?;
        }
        ClaimSubject::WholeProblem => writer.byte(5),
    }
    match &claim.scope {
        ClaimScope::WholeCampaign => writer.byte(0),
        ClaimScope::Cases(cases) => {
            writer.byte(1);
            writer.count(cases.len(), "claim case scope")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
        ClaimScope::Stratum { definition, cases } => {
            writer.byte(2);
            encode_source_ref(writer, definition)?;
            writer.count(cases.len(), "claim stratum cases")?;
            for case in cases {
                encode_case_id(writer, case)?;
            }
        }
    }
    Ok(())
}

fn decode_claim(
    reader: &mut CanonicalReader<'_>,
) -> Result<TypedIdentifiabilityClaim, IdentifiabilityError> {
    let id = ClaimId::try_new(reader.token("claim id")?)?;
    let information = match reader.byte("claim information regime")? {
        0 => InformationRegime::StructuralExactModel,
        1 => InformationRegime::ExactInputOutputMap,
        2 => InformationRegime::NoisyFiniteData,
        3 => InformationRegime::PosteriorUnderDeclaredPrior {
            joint_prior: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown information-regime tag {tag}"),
            });
        }
    };
    let extent = match reader.byte("claim extent")? {
        0 => IdentifiabilityExtent::Local,
        1 => IdentifiabilityExtent::Global,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown claim-extent tag {tag}"),
            });
        }
    };
    let fiber = match reader.byte("claim fiber structure")? {
        0 => FiberStructure::Unique,
        1 => {
            let maximum_cardinality = match reader.byte("finite-to-one maximum option")? {
                0 => None,
                1 => Some(FiberCardinalityBound::UniformU64(
                    reader.u64("finite-to-one maximum")?,
                )),
                2 => Some(FiberCardinalityBound::SymbolicProfile(decode_source_ref(
                    reader,
                )?)),
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown finite-to-one maximum option tag {tag}"),
                    });
                }
            };
            FiberStructure::FiniteToOne {
                maximum_cardinality,
            }
        }
        2 => FiberStructure::DiscreteOrbit {
            action: decode_gauge_action_reference(reader)?,
        },
        3 => FiberStructure::OrbitQuotientUnique {
            action: decode_gauge_action_reference(reader)?,
        },
        4 => FiberStructure::PositiveDimensional {
            lower_bound: match reader.byte("fiber-dimension lower-bound kind")? {
                0 => FiberDimensionLowerBound::Finite {
                    minimum_dimension: reader.u64("positive-dimensional fiber minimum")?,
                },
                1 => FiberDimensionLowerBound::InfiniteDimensional {
                    model_space: decode_source_ref(reader)?,
                },
                tag => {
                    return Err(IdentifiabilityError::Canonical {
                        at: reader.at.saturating_sub(1),
                        detail: format!("unknown fiber-dimension lower-bound tag {tag}"),
                    });
                }
            },
        },
        5 => FiberStructure::Stratified {
            strata: decode_source_ref(reader)?,
        },
        6 => FiberStructure::MixedOrbit {
            action: decode_gauge_action_reference(reader)?,
        },
        7 => FiberStructure::GeneralizedQuotientUnique {
            action: decode_gauge_action_reference(reader)?,
            equivalence: decode_source_ref(reader)?,
        },
        8 => FiberStructure::StratifiedOrbit {
            action: decode_gauge_action_reference(reader)?,
            orbit_type_profile: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown claim-fiber tag {tag}"),
            });
        }
    };
    let quantifier = match reader.byte("claim quantifier")? {
        0 => ClaimQuantifier::AtRealization {
            realization: decode_source_ref(reader)?,
        },
        1 => ClaimQuantifier::AlmostEverywhere {
            measure: decode_source_ref(reader)?,
        },
        2 => ClaimQuantifier::ForAll {
            domain: decode_source_ref(reader)?,
        },
        3 => ClaimQuantifier::ProbabilityAtLeast {
            probability: reader.f64("claim probability")?,
            measure: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown claim-quantifier tag {tag}"),
            });
        }
    };
    let scalar_domain = match reader.byte("claim scalar domain")? {
        0 => ScalarDomain::Real,
        1 => ScalarDomain::Complex {
            extension: decode_source_ref(reader)?,
        },
        2 => ScalarDomain::MixedDiscreteContinuous {
            stratification: decode_source_ref(reader)?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown scalar-domain tag {tag}"),
            });
        }
    };
    let subject = match reader.byte("claim subject")? {
        0 => ClaimSubject::Parameter(decode_role(reader)?),
        1 => ClaimSubject::ParameterSet(decode_role_set(reader, "claim parameter set")?),
        2 => ClaimSubject::DerivedFunctional {
            definition: decode_source_ref(reader)?,
            parameters: decode_role_set(reader, "derived-functional parameters")?,
        },
        3 => ClaimSubject::Influence(InfluenceId::try_new(reader.token("influence id")?)?),
        4 => ClaimSubject::GaugeAction(decode_gauge_action_reference(reader)?),
        5 => ClaimSubject::WholeProblem,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown claim-subject tag {tag}"),
            });
        }
    };
    let scope = match reader.byte("claim scope")? {
        0 => ClaimScope::WholeCampaign,
        1 => {
            let count = reader.count("claim case scope")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "claim case scope",
                        id: case.to_string(),
                    });
                }
            }
            ClaimScope::Cases(cases)
        }
        2 => {
            let definition = decode_source_ref(reader)?;
            let count = reader.count("claim stratum cases")?;
            let mut cases = BTreeSet::new();
            for _ in 0..count {
                let case = decode_case_id(reader)?;
                if !cases.insert(case.clone()) {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "claim stratum case",
                        id: case.to_string(),
                    });
                }
            }
            ClaimScope::Stratum { definition, cases }
        }
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown claim-scope tag {tag}"),
            });
        }
    };
    Ok(TypedIdentifiabilityClaim::new(
        id,
        information,
        extent,
        fiber,
        quantifier,
        scalar_domain,
        subject,
        scope,
    ))
}

fn encode_claim_assessment(
    writer: &mut CanonicalWriter,
    assessment: &ClaimAssessment,
) -> Result<(), IdentifiabilityError> {
    match assessment {
        ClaimAssessment::ClaimedEstablished {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound,
            gauge_resolutions,
        } => {
            writer.byte(0);
            encode_source_ref(writer, method)?;
            encode_source_ref(writer, receipt)?;
            encode_source_ref(writer, metric)?;
            encode_source_ref(writer, nondimensionalization)?;
            writer.f64(*certified_error_bound);
            writer.count(gauge_resolutions.len(), "gauge resolution evidence")?;
            for (action, evidence) in gauge_resolutions {
                if action != &evidence.action {
                    return Err(IdentifiabilityError::SourceMismatch {
                        field: "gauge resolution map key",
                    });
                }
                encode_gauge_action_reference(writer, action)?;
                writer.byte(match &evidence.disposition {
                    GaugeResolutionDisposition::CandidateRefuted => 0,
                    GaugeResolutionDisposition::NoProjectionOnSubject => 1,
                    GaugeResolutionDisposition::SubjectDescendsToQuotient => 2,
                    GaugeResolutionDisposition::BrokenByJointInformation => 3,
                    GaugeResolutionDisposition::TrivialResidualIntersection => 4,
                    GaugeResolutionDisposition::ConsistentWithClaimedFiber => 5,
                });
                encode_source_ref(writer, &evidence.method)?;
                encode_source_ref(writer, &evidence.receipt)?;
            }
        }
        ClaimAssessment::ClaimedRefuted {
            method,
            receipt,
            metric,
            nondimensionalization,
            certified_error_bound,
        } => {
            writer.byte(1);
            encode_source_ref(writer, method)?;
            encode_source_ref(writer, receipt)?;
            encode_source_ref(writer, metric)?;
            encode_source_ref(writer, nondimensionalization)?;
            writer.f64(*certified_error_bound);
        }
        ClaimAssessment::ClaimedInconclusive {
            method,
            receipt,
            reason,
        } => {
            writer.byte(2);
            match method {
                Some(source) => {
                    writer.byte(1);
                    encode_source_ref(writer, source)?;
                }
                None => writer.byte(0),
            }
            match receipt {
                Some(source) => {
                    writer.byte(1);
                    encode_source_ref(writer, source)?;
                }
                None => writer.byte(0),
            }
            writer.text(reason, "inconclusive reason")?;
        }
        ClaimAssessment::NotAssessed { reason } => {
            writer.byte(3);
            writer.text(reason, "not-assessed reason")?;
        }
    }
    Ok(())
}

fn decode_optional_source_ref(
    reader: &mut CanonicalReader<'_>,
    field: &'static str,
) -> Result<Option<SourceRef>, IdentifiabilityError> {
    match reader.byte(field)? {
        0 => Ok(None),
        1 => Ok(Some(decode_source_ref(reader)?)),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("invalid {field} tag {tag}"),
        }),
    }
}

fn decode_claim_assessment(
    reader: &mut CanonicalReader<'_>,
) -> Result<ClaimAssessment, IdentifiabilityError> {
    match reader.byte("claim assessment")? {
        0 => {
            let method = decode_source_ref(reader)?;
            let receipt = decode_source_ref(reader)?;
            let metric = decode_source_ref(reader)?;
            let nondimensionalization = decode_source_ref(reader)?;
            let certified_error_bound = reader.f64("certified claim error")?;
            let resolution_count = reader.count("gauge resolution evidence")?;
            let mut gauge_resolutions = BTreeMap::new();
            for _ in 0..resolution_count {
                let action = decode_gauge_action_reference(reader)?;
                let disposition = match reader.byte("gauge resolution disposition")? {
                    0 => GaugeResolutionDisposition::CandidateRefuted,
                    1 => GaugeResolutionDisposition::NoProjectionOnSubject,
                    2 => GaugeResolutionDisposition::SubjectDescendsToQuotient,
                    3 => GaugeResolutionDisposition::BrokenByJointInformation,
                    4 => GaugeResolutionDisposition::TrivialResidualIntersection,
                    5 => GaugeResolutionDisposition::ConsistentWithClaimedFiber,
                    tag => {
                        return Err(IdentifiabilityError::Canonical {
                            at: reader.at.saturating_sub(1),
                            detail: format!("unknown gauge-resolution disposition tag {tag}"),
                        });
                    }
                };
                let evidence = GaugeResolutionEvidence::new(
                    action.clone(),
                    disposition,
                    decode_source_ref(reader)?,
                    decode_source_ref(reader)?,
                );
                if gauge_resolutions.insert(action.clone(), evidence).is_some() {
                    return Err(IdentifiabilityError::Duplicate {
                        field: "gauge resolution evidence",
                        id: format!("{action:?}"),
                    });
                }
            }
            Ok(ClaimAssessment::ClaimedEstablished {
                method,
                receipt,
                metric,
                nondimensionalization,
                certified_error_bound,
                gauge_resolutions,
            })
        }
        1 => Ok(ClaimAssessment::ClaimedRefuted {
            method: decode_source_ref(reader)?,
            receipt: decode_source_ref(reader)?,
            metric: decode_source_ref(reader)?,
            nondimensionalization: decode_source_ref(reader)?,
            certified_error_bound: reader.f64("certified claim error")?,
        }),
        2 => Ok(ClaimAssessment::ClaimedInconclusive {
            method: decode_optional_source_ref(reader, "inconclusive method option")?,
            receipt: decode_optional_source_ref(reader, "inconclusive receipt option")?,
            reason: reader.reason("inconclusive reason")?,
        }),
        3 => Ok(ClaimAssessment::NotAssessed {
            reason: reader.reason("not-assessed reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown claim-assessment tag {tag}"),
        }),
    }
}

fn encode_assessment_with_header_mode(
    assessment: &IdentifiabilityAssessment,
    exact_header: bool,
) -> Result<Vec<u8>, IdentifiabilityError> {
    check_assessment_identity_version(assessment.schema_version)?;
    let mut writer = CanonicalWriter::new();
    writer.raw(ASSESSMENT_MAGIC);
    writer.u32(assessment.schema_version);
    encode_header(&mut writer, &assessment.header, exact_header)?;
    writer.hash(assessment.problem_id.0);
    writer.hash(assessment.execution_id.0);
    writer.count(assessment.claims.len(), "identifiability claims")?;
    for claim in assessment.claims.values() {
        encode_claim(&mut writer, claim)?;
    }
    writer.count(assessment.evidence.len(), "claim assessments")?;
    for (id, conclusion) in &assessment.evidence {
        writer.text(id.as_str(), "claim id")?;
        encode_claim_assessment(&mut writer, conclusion)?;
    }
    encode_resolution_set(&mut writer, &assessment.source_authority)?;
    writer.finish()
}

fn encode_assessment(
    assessment: &IdentifiabilityAssessment,
) -> Result<Vec<u8>, IdentifiabilityError> {
    encode_assessment_with_header_mode(assessment, true)
}

fn encode_assessment_identity(
    assessment: &IdentifiabilityAssessment,
) -> Result<Vec<u8>, IdentifiabilityError> {
    encode_assessment_with_header_mode(assessment, false)
}

fn assessment_identity_hash(
    assessment: &IdentifiabilityAssessment,
) -> Result<AssessmentId, IdentifiabilityError> {
    Ok(AssessmentId(hash_domain(
        IDENTIFIABILITY_ASSESSMENT_IDENTITY_DOMAIN,
        &encode_assessment_identity(assessment)?,
    )))
}

fn decode_assessment(
    bytes: &[u8],
    problem: &AdmittedIdentifiabilityProblem,
    execution: &IdentifiabilityExecutionPlan,
    verified_sources: &SourceResolutionSet,
) -> Result<IdentifiabilityAssessment, IdentifiabilityError> {
    let mut reader = CanonicalReader::new(bytes)?;
    if reader.take(ASSESSMENT_MAGIC.len(), "assessment magic")? != ASSESSMENT_MAGIC {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "wrong identifiability-assessment magic".to_string(),
        });
    }
    check_assessment_identity_version(reader.u32("assessment schema version")?)?;
    let header = decode_header(&mut reader)?;
    let problem_id = ProblemId(reader.hash("assessment problem id")?);
    let execution_id = ExecutionId(reader.hash("assessment execution id")?);
    if problem_id != problem.problem_id || execution_id != execution.id()? {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "assessment problem/execution identity",
        });
    }
    let claim_count = reader.count("identifiability claims")?;
    preflight_collection_bytes(&reader, claim_count, 1, "identifiability claims")?;
    let mut claims = Vec::with_capacity(claim_count);
    for _ in 0..claim_count {
        claims.push(decode_claim(&mut reader)?);
    }
    let evidence_count = reader.count("claim assessments")?;
    preflight_collection_bytes(&reader, evidence_count, 1, "claim assessments")?;
    let mut evidence = Vec::with_capacity(evidence_count);
    for _ in 0..evidence_count {
        evidence.push((
            ClaimId::try_new(reader.token("claim id")?)?,
            decode_claim_assessment(&mut reader)?,
        ));
    }
    // Resolution evidence is retained in transport for identity/replay, but it
    // is not itself proof. Compare it to a caller-held, locally verified set;
    // never pass deserialized verification markers into the admitting
    // constructor.
    let transported_sources = decode_resolution_set(&mut reader)?;
    reader.finish()?;
    if transported_sources != *verified_sources {
        return Err(IdentifiabilityError::SourceMismatch {
            field: "assessment source-resolution replay",
        });
    }
    let assessment = IdentifiabilityAssessment::try_new(
        header,
        problem,
        execution,
        claims,
        evidence,
        verified_sources.clone(),
    )?;
    if assessment.canonical_bytes()?.as_slice() != bytes {
        return Err(IdentifiabilityError::Canonical {
            at: 0,
            detail: "non-canonical assessment encoding".to_string(),
        });
    }
    Ok(assessment)
}

#[cfg(test)]
mod decoder_resource_tests {
    use super::*;

    #[test]
    fn dense_correlation_decoder_rejects_oversize_before_entry_allocation() {
        let mut bytes = vec![1];
        bytes.extend_from_slice(&0_u32.to_le_bytes());
        bytes.extend_from_slice(
            &u32::try_from(MAX_VV_MATRIX_DIMENSION + 1)
                .expect("matrix limit fits u32")
                .to_le_bytes(),
        );
        let mut reader = CanonicalReader::new(&bytes).expect("bounded fixture");
        assert!(matches!(
            decode_joint_noise(&mut reader),
            Err(IdentifiabilityError::Cardinality {
                field: "correlation dimension",
                ..
            })
        ));
    }

    #[test]
    fn dense_correlation_decoder_rejects_truncation_before_entry_allocation() {
        let dimension = MAX_VV_MATRIX_DIMENSION;
        let entries = dimension * (dimension + 1) / 2;
        let mut bytes = vec![1];
        bytes.extend_from_slice(&0_u32.to_le_bytes());
        bytes.extend_from_slice(
            &u32::try_from(dimension)
                .expect("matrix limit fits u32")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(
            &u32::try_from(entries)
                .expect("entry count fits u32")
                .to_le_bytes(),
        );
        let mut reader = CanonicalReader::new(&bytes).expect("bounded fixture");
        assert!(matches!(
            decode_joint_noise(&mut reader),
            Err(IdentifiabilityError::Canonical { .. })
        ));
    }
}

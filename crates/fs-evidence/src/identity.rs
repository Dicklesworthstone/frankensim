//! Typed, canonical identities for evidence semantics.
//!
//! This module covers exact color-evidence graph replay, normalized validity
//! domains, standalone numerical/statistical certificate declarations, exact
//! two-fidelity observations and discrepancy-band declarations, model-form
//! evidence slices, model-card declarations with exact calibration sources,
//! an opaque strong-identity projection of locally certified scalar evidence,
//! its exact-source binding, its exact ordered provenance-lineage graph, and
//! its recomputed local decision assessment through separate schemas. It never
//! promotes [`crate::ProvenanceHash`] into strong identity: exact-source and
//! lineage helpers use retained values only as legacy FNV compatibility
//! checks and exclude them from strong frames. This module
//! publishes only unanchored [`IdentityReceipt`] values. Origin verification,
//! policy admission, structural [`crate::Certified`] consistency, and
//! scientific color rank remain separate axes.

use core::fmt;

pub use fs_blake3::identity::{
    CancellationProbe as EvidenceIdentityCancellationProbe,
    CanonicalLimits as EvidenceIdentityLimits, TrustState as EvidenceIdentityTrustState,
};
use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalSchema, ChildSpec, EvidenceNodeId, Field, FieldSpec,
    IdentityAdjudication, IdentityReceipt, LimitKind, ModelId, ObservedIdentity,
    OrderedBytesStreamError, SchemaId, SemanticId, SourceByteId, SourceId, StrongIdentity,
    WireType, adjudicate,
};
use fs_obs::ident::{BoundedIdentityBuilder, IdentityBuildError};

use crate::{
    Ambition, COLOR_ALGEBRA_VERSION, Certified, Color, ColorPayloadError, DecisionStatus,
    DiscrepancyBand, EscalationAdvice, FidelityPair, IntervalOp, ModelBracket, ModelCard,
    ModelEvidence, NumericalKind, Op, ProvenanceHash, StatisticalCertificate, UncertaintyBreakdown,
    UncertaintySource, ValidityDomain, color_identity_reason, color_leaf_identity_reason, compose,
    validate_color_payload,
};

/// Identity schema version for exact retained color-evidence sources.
pub const COLOR_EVIDENCE_SOURCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for color-evidence graph nodes.
pub const COLOR_EVIDENCE_NODE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for normalized evidence-validity domains.
pub const VALIDITY_DOMAIN_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one model-evidence semantic slice.
pub const MODEL_EVIDENCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one exact model-evidence/card declaration join.
pub const CARD_BOUND_MODEL_EVIDENCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one standalone numerical-certificate declaration.
pub const NUMERICAL_CERTIFICATE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one standalone statistical-certificate declaration.
pub const STATISTICAL_CERTIFICATE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one standalone two-fidelity observation.
pub const FIDELITY_PAIR_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one standalone discrepancy-band declaration.
pub const DISCREPANCY_BAND_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one exact admitted model-bracket declaration.
pub const MODEL_BRACKET_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one certified-f64 decision assessment.
pub const CERTIFIED_F64_DECISION_ASSESSMENT_IDENTITY_VERSION_V1: u32 = 1;
/// Semantic version of the bound breakdown, tie-break, status, and advice law.
pub const DECISION_ASSESSMENT_ALGORITHM_VERSION_V1: u32 = 1;
/// Identity schema version for one locally certified scalar-evidence projection.
pub const CERTIFIED_F64_EVIDENCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for an exact certified-f64 source artifact.
pub const CERTIFIED_F64_SOURCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one exact-source-bound certified-f64 projection.
pub const SOURCED_CERTIFIED_F64_EVIDENCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one exact ordered certified-f64 provenance node.
pub const CERTIFIED_F64_LINEAGE_NODE_IDENTITY_VERSION_V1: u32 = 1;
/// Semantic version of the exact-source/ordered-parent lineage law.
pub const CERTIFIED_F64_LINEAGE_ALGORITHM_VERSION_V1: u32 = 1;
/// Identity schema version for exact model-card calibration source bytes.
pub const MODEL_CARD_CALIBRATION_SOURCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one helper-validated model-card declaration.
pub const MODEL_CARD_IDENTITY_VERSION_V1: u32 = 1;
/// Hard allocation ceiling for one canonical output-color payload.
pub const MAX_COLOR_EVIDENCE_NODE_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for the ordered axes field of one validity domain.
pub const MAX_VALIDITY_DOMAIN_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each variable model-evidence identity field.
pub const MAX_MODEL_EVIDENCE_IDENTITY_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for the ordered typed-card field of one bound model-evidence join.
pub const MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each variable certified-evidence field.
pub const MAX_CERTIFIED_F64_EVIDENCE_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each exact certified-f64 source field.
pub const MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for one certified-f64 lineage parent field.
pub const MAX_CERTIFIED_F64_LINEAGE_PARENT_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each variable model-card declaration field.
pub const MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard retained-byte ceiling for one model-card calibration source.
pub const MAX_MODEL_CARD_CALIBRATION_SOURCE_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for the ordered parameter field of one fidelity pair.
pub const MAX_FIDELITY_PAIR_PARAMETERS_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard parameter-count ceiling matching discrepancy-fit v1 admission.
pub const MAX_FIDELITY_PAIR_PARAMETERS_V1: u64 = 1_024;
/// Hard payload ceiling for one exact model-bracket member field.
pub const MAX_MODEL_BRACKET_MEMBERS_FIELD_BYTES_V1: u64 = 1 << 20;
/// Minimum member count needed to express model-choice spread.
pub const MIN_MODEL_BRACKET_MEMBERS_V1: u64 = 2;
/// Hard member-count ceiling matching model-bracket v1 admission.
pub const MAX_MODEL_BRACKET_MEMBERS_V1: u64 = 1_024;
const _: () =
    assert!(MIN_MODEL_BRACKET_MEMBERS_V1 as usize == crate::discrepancy::MIN_BRACKET_MEMBERS);
const _: () =
    assert!(MAX_MODEL_BRACKET_MEMBERS_V1 as usize == crate::discrepancy::MAX_BRACKET_MEMBERS);
/// Non-semantic scatter/gather writes emitted for each streamed axis row.
const VALIDITY_DOMAIN_STREAM_CHUNKS_PER_AXIS_V1: u64 = 4;
/// Non-semantic scatter/gather writes emitted for each sensitivity row.
const CERTIFIED_F64_SENSITIVITY_STREAM_CHUNKS_PER_ROW_V1: u64 = 3;
/// Non-semantic scatter/gather writes emitted for each fidelity parameter row.
const FIDELITY_PAIR_STREAM_CHUNKS_PER_PARAMETER_V1: u64 = 3;
/// Non-semantic scatter/gather writes emitted for each model-bracket row.
const MODEL_BRACKET_STREAM_CHUNKS_PER_MEMBER_V1: u64 = 3;

/// Canonical identity schema for one retained source that may root a color
/// evidence graph. The resulting identity is content-bound but untrusted.
pub enum ColorEvidenceSourceIdentitySchemaV1 {}

impl CanonicalSchema for ColorEvidenceSourceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.color-evidence-source.v1";
    const NAME: &'static str = "color-evidence-source";
    const VERSION: u32 = COLOR_EVIDENCE_SOURCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact retained source schema domain, source schema version, and canonical source bytes; no origin or scientific authority";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source-domain", WireType::Utf8),
        FieldSpec::required("source-schema-version", WireType::U64),
        FieldSpec::required("canonical-source", WireType::Bytes),
    ];
}

/// Low-level canonical-frame identity for one retained color-evidence source.
///
/// Direct generic encoder output proves only schema-shaped framing. The
/// helper-enforced source-domain and version invariants belong to
/// [`ColorEvidenceSourceV1`].
pub type ColorEvidenceSourceIdV1 = SourceId<ColorEvidenceSourceIdentitySchemaV1>;

/// Canonical identity schema for one normalized evidence-validity domain.
pub enum ValidityDomainIdentitySchemaV1 {}

impl CanonicalSchema for ValidityDomainIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.validity-domain.v1";
    const NAME: &'static str = "validity-domain";
    const VERSION: u32 = VALIDITY_DOMAIN_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "sorted exact validity-axis UTF-8 bytes and finite IEEE-754 bounds; no empirical membership or model authority";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("axes", WireType::OrderedBytes)];
}

/// Low-level canonical-frame identity for one normalized validity domain.
///
/// Direct generic encoder output proves only schema-shaped framing. Exact axis
/// decoding, finite ordered bounds, normalization, and resource admission
/// belong to [`IdentifiedValidityDomainV1`].
pub type ValidityDomainIdV1 = SemanticId<ValidityDomainIdentitySchemaV1>;

/// Low-level producer receipt underlying an opaque validated domain.
pub type ValidityDomainReceiptV1 = IdentityReceipt<ValidityDomainIdV1>;

/// A normalized validity domain kept attached to its unanchored identity.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentifiedValidityDomainV1 {
    domain: ValidityDomain,
    receipt: ValidityDomainReceiptV1,
}

impl IdentifiedValidityDomainV1 {
    /// Read-only normalized domain committed by this identity.
    #[must_use]
    pub const fn domain(&self) -> &ValidityDomain {
        &self.domain
    }

    /// Typed semantic identity.
    #[must_use]
    pub const fn id(&self) -> ValidityDomainIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> ValidityDomainReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the plain domain.
    #[must_use]
    pub fn into_domain(self) -> ValidityDomain {
        self.domain
    }
}

static MODEL_EVIDENCE_VALIDITY_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<ValidityDomainIdV1>();

/// Canonical semantic schema for one model-form evidence slice.
pub enum ModelEvidenceIdentitySchemaV1 {}

impl CanonicalSchema for ModelEvidenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.model-evidence.v1";
    const NAME: &'static str = "model-evidence";
    const VERSION: u32 = MODEL_EVIDENCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact model-card name and assumption sets, typed declared validity, discrepancy bits, and in-domain claim state; no card content, units, evaluation point, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("model-card-names", WireType::CanonicalSet),
        FieldSpec::required("assumptions", WireType::CanonicalSet),
        FieldSpec::child_of("validity", &MODEL_EVIDENCE_VALIDITY_CHILD_V1),
        FieldSpec::required("discrepancy-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("in-domain", WireType::Bool),
    ];
}

/// Low-level schema-shaped identity for one model-evidence semantic frame.
///
/// Only [`IdentifiedModelEvidenceV1`] proves correspondence with an attached
/// public [`ModelEvidence`] and helper-validated validity child.
pub type ModelEvidenceIdV1 = SemanticId<ModelEvidenceIdentitySchemaV1>;

/// Low-level producer receipt for one model-evidence semantic frame.
pub type ModelEvidenceReceiptV1 = IdentityReceipt<ModelEvidenceIdV1>;

/// A model-form evidence slice kept attached to its unanchored semantic
/// identity and helper-built validity receipt.
#[derive(Debug, Clone)]
pub struct IdentifiedModelEvidenceV1 {
    model_evidence: ModelEvidence,
    validity_receipt: ValidityDomainReceiptV1,
    receipt: ModelEvidenceReceiptV1,
}

impl IdentifiedModelEvidenceV1 {
    /// Read-only model evidence committed by this identity.
    #[must_use]
    pub const fn model_evidence(&self) -> &ModelEvidence {
        &self.model_evidence
    }

    /// Typed semantic identity.
    #[must_use]
    pub const fn id(&self) -> ModelEvidenceIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored semantic receipt.
    #[must_use]
    pub const fn receipt(&self) -> ModelEvidenceReceiptV1 {
        self.receipt
    }

    /// Typed normalized validity identity bound as a child.
    #[must_use]
    pub const fn validity_id(&self) -> ValidityDomainIdV1 {
        self.validity_receipt.id()
    }

    /// Complete helper-built validity receipt.
    #[must_use]
    pub const fn validity_receipt(&self) -> ValidityDomainReceiptV1 {
        self.validity_receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the model evidence.
    #[must_use]
    pub fn into_model_evidence(self) -> ModelEvidence {
        self.model_evidence
    }
}

/// Canonical source schema for exact calibration artifact bytes supplied while
/// migrating a model card away from its legacy FNV correlation token.
pub enum ModelCardCalibrationSourceIdentitySchemaV1 {}

impl CanonicalSchema for ModelCardCalibrationSourceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.model-card-calibration-source.v1";
    const NAME: &'static str = "model-card-calibration-source";
    const VERSION: u32 = MODEL_CARD_CALIBRATION_SOURCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact caller-retained calibration artifact bytes supplied for legacy crosswalk; no format, origin, custody, currentness, efficacy, or authority";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "canonical-calibration-artifact",
        WireType::Bytes,
    )];
}

/// Low-level schema-shaped source-byte identity for exact calibration bytes.
///
/// Only [`IdentifiedModelCardV1`] proves that a present source was cross-checked
/// against the attached card's legacy correlation token. For an uncalibrated
/// card, the helper binds this schema's empty-byte root as an absence sentinel.
pub type ModelCardCalibrationSourceIdV1 = SourceByteId<ModelCardCalibrationSourceIdentitySchemaV1>;

/// Low-level receipt for one calibration source-byte frame.
pub type ModelCardCalibrationSourceReceiptV1 = IdentityReceipt<ModelCardCalibrationSourceIdV1>;

static MODEL_CARD_VALIDITY_CHILD_V1: ChildSpec = ChildSpec::for_identity::<ValidityDomainIdV1>();
static MODEL_CARD_CALIBRATION_SOURCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<ModelCardCalibrationSourceIdV1>();

/// Canonical model identity for one declaration plus its exact calibration
/// source binding. The legacy FNV value is checked but never framed.
pub enum ModelCardIdentitySchemaV1 {}

impl CanonicalSchema for ModelCardIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.model-card.v1";
    const NAME: &'static str = "model-card";
    const VERSION: u32 = MODEL_CARD_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "model declaration frame with exact names, ambition, canonical assumptions and failures, typed validity, discrepancy bits, and exact calibration source bytes; legacy FNV excluded; no model correctness, registry admission, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("name", WireType::Utf8),
        FieldSpec::required("version", WireType::Utf8),
        FieldSpec::required("ambition", WireType::Variant),
        FieldSpec::required("assumptions", WireType::CanonicalSet),
        FieldSpec::child_of("validity", &MODEL_CARD_VALIDITY_CHILD_V1),
        FieldSpec::required("known-failures", WireType::CanonicalSet),
        FieldSpec::required("calibration-present", WireType::Bool),
        FieldSpec::child_of(
            "calibration-source",
            &MODEL_CARD_CALIBRATION_SOURCE_CHILD_V1,
        ),
        FieldSpec::required("discrepancy-rel-ieee754-bits", WireType::U64),
    ];
}

/// Low-level schema-shaped model-card identity.
///
/// Direct encoder output does not prove consistency with a retained
/// [`ModelCard`], normalized validity, or exact calibration bytes.
pub type ModelCardIdV1 = ModelId<ModelCardIdentitySchemaV1>;

/// Low-level producer receipt for one model-card frame.
pub type ModelCardReceiptV1 = IdentityReceipt<ModelCardIdV1>;

/// A model declaration kept attached to its unanchored model identity, exact
/// optional calibration bytes, and helper-built child receipts.
#[derive(Debug, Clone)]
pub struct IdentifiedModelCardV1 {
    card: ModelCard,
    calibration_bytes: Option<Vec<u8>>,
    validity_receipt: ValidityDomainReceiptV1,
    calibration_source_receipt: ModelCardCalibrationSourceReceiptV1,
    receipt: ModelCardReceiptV1,
}

impl IdentifiedModelCardV1 {
    /// Read-only model card committed by this identity.
    #[must_use]
    pub const fn card(&self) -> &ModelCard {
        &self.card
    }

    /// Exact supplied calibration bytes, if the card declares calibration.
    #[must_use]
    pub fn calibration_bytes(&self) -> Option<&[u8]> {
        self.calibration_bytes.as_deref()
    }

    /// Typed model identity.
    #[must_use]
    pub const fn id(&self) -> ModelCardIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored model receipt.
    #[must_use]
    pub const fn receipt(&self) -> ModelCardReceiptV1 {
        self.receipt
    }

    /// Typed normalized validity identity bound as a child.
    #[must_use]
    pub const fn validity_id(&self) -> ValidityDomainIdV1 {
        self.validity_receipt.id()
    }

    /// Complete helper-built validity receipt.
    #[must_use]
    pub const fn validity_receipt(&self) -> ValidityDomainReceiptV1 {
        self.validity_receipt
    }

    /// Exact calibration source-byte identity for a calibrated card.
    ///
    /// Returns `None` for an uncalibrated card so the internal empty-byte
    /// absence sentinel cannot be mistaken for a declared artifact.
    #[must_use]
    pub fn calibration_source_id(&self) -> Option<ModelCardCalibrationSourceIdV1> {
        self.calibration_bytes
            .as_ref()
            .map(|_| self.calibration_source_receipt.id())
    }

    /// Complete calibration source-byte receipt for a calibrated card.
    #[must_use]
    pub fn calibration_source_receipt(&self) -> Option<ModelCardCalibrationSourceReceiptV1> {
        self.calibration_bytes
            .as_ref()
            .map(|_| self.calibration_source_receipt)
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment without discarding retained source
    /// bytes.
    #[must_use]
    pub fn into_parts(self) -> (ModelCard, Option<Vec<u8>>) {
        (self.card, self.calibration_bytes)
    }
}

static CARD_BOUND_MODEL_EVIDENCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<ModelEvidenceIdV1>();
static CARD_BOUND_MODEL_CARD_CHILD_V1: ChildSpec = ChildSpec::for_identity::<ModelCardIdV1>();

/// Canonical identity schema joining one exact model-evidence declaration to
/// the ordered typed model-card declarations named by that evidence.
pub enum CardBoundModelEvidenceIdentitySchemaV1 {}

impl CanonicalSchema for CardBoundModelEvidenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.card-bound-model-evidence.v1";
    const NAME: &'static str = "card-bound-model-evidence";
    const VERSION: u32 = CARD_BOUND_MODEL_EVIDENCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "typed model-evidence child plus exact ordered typed model-card declarations after byte-exact positional name matching; no derivation, execution, registry membership, occupancy, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("model-evidence", &CARD_BOUND_MODEL_EVIDENCE_CHILD_V1),
        FieldSpec::ordered_children_of("model-card-declarations", &CARD_BOUND_MODEL_CARD_CHILD_V1),
    ];
}

/// Low-level schema-shaped identity for one model-evidence/card declaration
/// join.
///
/// Only [`IdentifiedCardBoundModelEvidenceV1`] proves exact positional name
/// correspondence with retained opaque children.
pub type CardBoundModelEvidenceIdV1 = SemanticId<CardBoundModelEvidenceIdentitySchemaV1>;

/// Low-level producer receipt for one card-bound model-evidence frame.
pub type CardBoundModelEvidenceReceiptV1 = IdentityReceipt<CardBoundModelEvidenceIdV1>;

/// One opaque model-evidence child kept attached to the exact ordered opaque
/// model-card declarations named by its canonical card-name set.
#[derive(Debug, Clone)]
pub struct IdentifiedCardBoundModelEvidenceV1 {
    model_evidence: IdentifiedModelEvidenceV1,
    model_cards: Vec<IdentifiedModelCardV1>,
    receipt: CardBoundModelEvidenceReceiptV1,
}

impl IdentifiedCardBoundModelEvidenceV1 {
    /// Opaque model-evidence child committed by this join.
    #[must_use]
    pub const fn model_evidence(&self) -> &IdentifiedModelEvidenceV1 {
        &self.model_evidence
    }

    /// Exact ordered opaque model-card declarations committed by this join.
    #[must_use]
    pub fn model_cards(&self) -> &[IdentifiedModelCardV1] {
        &self.model_cards
    }

    /// Typed card-bound model-evidence identity.
    #[must_use]
    pub const fn id(&self) -> CardBoundModelEvidenceIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> CardBoundModelEvidenceReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the join while preserving every retained opaque child.
    #[must_use]
    pub fn into_parts(self) -> (IdentifiedModelEvidenceV1, Vec<IdentifiedModelCardV1>) {
        (self.model_evidence, self.model_cards)
    }
}

/// Canonical semantic schema for the exact structural state of one standalone
/// numerical-certificate declaration.
pub enum NumericalCertificateIdentitySchemaV1 {}

impl CanonicalSchema for NumericalCertificateIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.numerical-certificate.v1";
    const NAME: &'static str = "numerical-certificate";
    const VERSION: u32 = NUMERICAL_CERTIFICATE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact admitted numerical-certificate kind and IEEE-754 bound bits; no QoI, value, units, containment, rigor, derivation, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("kind", WireType::Variant),
        FieldSpec::required("lo-ieee754-bits", WireType::U64),
        FieldSpec::required("hi-ieee754-bits", WireType::U64),
    ];
}

/// Low-level schema-shaped identity for one numerical-certificate frame.
///
/// Only [`IdentifiedNumericalCertificateV1`] proves correspondence with a
/// retained certificate that passed this module's structural admission.
pub type NumericalCertificateIdV1 = SemanticId<NumericalCertificateIdentitySchemaV1>;

/// Low-level producer receipt for one numerical-certificate frame.
pub type NumericalCertificateReceiptV1 = IdentityReceipt<NumericalCertificateIdV1>;

/// A structurally admitted numerical-certificate declaration kept attached to
/// its unanchored semantic identity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdentifiedNumericalCertificateV1 {
    certificate: crate::NumericalCertificate,
    receipt: NumericalCertificateReceiptV1,
}

impl IdentifiedNumericalCertificateV1 {
    /// Read-only certificate declaration committed by this identity.
    #[must_use]
    pub const fn certificate(&self) -> &crate::NumericalCertificate {
        &self.certificate
    }

    /// Typed structural identity.
    #[must_use]
    pub const fn id(&self) -> NumericalCertificateIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> NumericalCertificateReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the plain certificate.
    #[must_use]
    pub const fn into_certificate(self) -> crate::NumericalCertificate {
        self.certificate
    }
}

/// Canonical semantic schema for the exact structural state of one standalone
/// statistical-certificate declaration.
pub enum StatisticalCertificateIdentitySchemaV1 {}

impl CanonicalSchema for StatisticalCertificateIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.statistical-certificate.v1";
    const NAME: &'static str = "statistical-certificate";
    const VERSION: u32 = STATISTICAL_CERTIFICATE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact admitted statistical-certificate variant and IEEE-754 payload bits; no hypothesis, estimand, method, sample, stopping or dependence context, coverage, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("certificate", WireType::Variant)];
}

/// Low-level schema-shaped identity for one statistical-certificate frame.
///
/// Only [`IdentifiedStatisticalCertificateV1`] proves correspondence with a
/// retained certificate whose numeric parameters passed local shape checks.
pub type StatisticalCertificateIdV1 = SemanticId<StatisticalCertificateIdentitySchemaV1>;

/// Low-level producer receipt for one statistical-certificate frame.
pub type StatisticalCertificateReceiptV1 = IdentityReceipt<StatisticalCertificateIdV1>;

/// A structurally admitted statistical-certificate declaration kept attached
/// to its unanchored semantic identity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdentifiedStatisticalCertificateV1 {
    certificate: StatisticalCertificate,
    receipt: StatisticalCertificateReceiptV1,
}

impl IdentifiedStatisticalCertificateV1 {
    /// Read-only certificate declaration committed by this identity.
    #[must_use]
    pub const fn certificate(&self) -> &StatisticalCertificate {
        &self.certificate
    }

    /// Typed structural identity.
    #[must_use]
    pub const fn id(&self) -> StatisticalCertificateIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> StatisticalCertificateReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the plain certificate.
    #[must_use]
    pub const fn into_certificate(self) -> StatisticalCertificate {
        self.certificate
    }
}

/// Canonical semantic schema for one exact two-fidelity observation.
pub enum FidelityPairIdentitySchemaV1 {}

impl CanonicalSchema for FidelityPairIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.fidelity-pair.v1";
    const NAME: &'static str = "fidelity-pair";
    const VERSION: u32 = FIDELITY_PAIR_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "sorted exact parameter-name and finite IEEE-754 coordinate rows plus finite low/high QoI bits; no units, model identity, run, source, pairing authenticity, reference truth, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("parameters", WireType::OrderedBytes),
        FieldSpec::required("lo-fi-qoi-ieee754-bits", WireType::U64),
        FieldSpec::required("hi-fi-qoi-ieee754-bits", WireType::U64),
    ];
}

/// Low-level schema-shaped identity for one two-fidelity observation frame.
///
/// Only [`IdentifiedFidelityPairV1`] proves correspondence with a retained pair
/// whose names and numeric fields passed local structural admission.
pub type FidelityPairIdV1 = SemanticId<FidelityPairIdentitySchemaV1>;

/// Low-level producer receipt for one two-fidelity observation frame.
pub type FidelityPairReceiptV1 = IdentityReceipt<FidelityPairIdV1>;

/// A structurally admitted two-fidelity observation kept attached to its
/// unanchored semantic identity.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentifiedFidelityPairV1 {
    pair: FidelityPair,
    receipt: FidelityPairReceiptV1,
}

impl IdentifiedFidelityPairV1 {
    /// Read-only fidelity pair committed by this identity.
    #[must_use]
    pub const fn pair(&self) -> &FidelityPair {
        &self.pair
    }

    /// Typed structural identity.
    #[must_use]
    pub const fn id(&self) -> FidelityPairIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> FidelityPairReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the plain fidelity pair.
    #[must_use]
    pub fn into_pair(self) -> FidelityPair {
        self.pair
    }
}

/// Canonical semantic schema for one standalone discrepancy-band declaration.
pub enum DiscrepancyBandIdentitySchemaV1 {}

impl CanonicalSchema for DiscrepancyBandIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.discrepancy-band.v1";
    const NAME: &'static str = "discrepancy-band";
    const VERSION: u32 = DISCREPANCY_BAND_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact admitted mean and maximum relative-discrepancy IEEE-754 bits; no corpus, domain, pair count, query point, metric, denominator, aggregation, confidence, derivation, rigor, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("mean-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("max-rel-ieee754-bits", WireType::U64),
    ];
}

/// Low-level schema-shaped identity for one discrepancy-band frame.
///
/// Only [`IdentifiedDiscrepancyBandV1`] proves correspondence with a retained
/// band that passed this module's structural admission.
pub type DiscrepancyBandIdV1 = SemanticId<DiscrepancyBandIdentitySchemaV1>;

/// Low-level producer receipt for one discrepancy-band frame.
pub type DiscrepancyBandReceiptV1 = IdentityReceipt<DiscrepancyBandIdV1>;

/// A structurally admitted discrepancy band kept attached to its unanchored
/// semantic identity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdentifiedDiscrepancyBandV1 {
    band: DiscrepancyBand,
    receipt: DiscrepancyBandReceiptV1,
}

impl IdentifiedDiscrepancyBandV1 {
    /// Read-only discrepancy band committed by this identity.
    #[must_use]
    pub const fn band(&self) -> &DiscrepancyBand {
        &self.band
    }

    /// Typed structural identity.
    #[must_use]
    pub const fn id(&self) -> DiscrepancyBandIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> DiscrepancyBandReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the plain band.
    #[must_use]
    pub const fn into_band(self) -> DiscrepancyBand {
        self.band
    }
}

/// Canonical semantic schema for one exact admitted model-bracket mapping.
pub enum ModelBracketIdentitySchemaV1 {}

impl CanonicalSchema for ModelBracketIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.model-bracket.v1";
    const NAME: &'static str = "model-bracket";
    const VERSION: u32 = MODEL_BRACKET_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "model-name-ordered exact admitted member name and finite QoI IEEE-754 rows; no model-card contents, model execution, units, quantity kind, bracket evidence, source, provenance, scientific authority, origin, or trust";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("members", WireType::OrderedBytes)];
}

/// Low-level schema-shaped identity for one model-bracket member mapping.
///
/// Only [`IdentifiedModelBracketV1`] proves correspondence with a retained
/// bracket whose canonical members passed local structural admission.
pub type ModelBracketIdV1 = SemanticId<ModelBracketIdentitySchemaV1>;

/// Low-level producer receipt for one model-bracket member mapping.
pub type ModelBracketReceiptV1 = IdentityReceipt<ModelBracketIdV1>;

/// An admitted model bracket kept attached to its unanchored semantic identity.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentifiedModelBracketV1 {
    bracket: ModelBracket,
    receipt: ModelBracketReceiptV1,
}

impl IdentifiedModelBracketV1 {
    /// Read-only canonical bracket committed by this identity.
    #[must_use]
    pub const fn bracket(&self) -> &ModelBracket {
        &self.bracket
    }

    /// Typed structural identity.
    #[must_use]
    pub const fn id(&self) -> ModelBracketIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> ModelBracketReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the canonical bracket.
    #[must_use]
    pub fn into_bracket(self) -> ModelBracket {
        self.bracket
    }
}

static CERTIFIED_F64_VALIDITY_CHILD_V1: ChildSpec = ChildSpec::for_identity::<ValidityDomainIdV1>();

/// Canonical identity schema for the strong semantic projection of one locally
/// certified scalar. This schema is intentionally unqualified by units or
/// quantity kind because those concepts are absent from `Evidence<f64>` today.
pub enum CertifiedF64EvidenceIdentitySchemaV1 {}

impl CanonicalSchema for CertifiedF64EvidenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.certified-f64-evidence.v1";
    const NAME: &'static str = "certified-f64-evidence";
    const VERSION: u32 = CERTIFIED_F64_EVIDENCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "locally certified unqualified scalar evidence semantics with typed validity and adjoint-claim presence; legacy FNV values excluded; no units, quantity kind, origin, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("value", WireType::FiniteF64),
        FieldSpec::required("qoi", WireType::FiniteF64),
        FieldSpec::required("numerical-kind", WireType::Variant),
        FieldSpec::required("numerical-lo", WireType::FiniteF64),
        FieldSpec::required("numerical-hi", WireType::FiniteF64),
        FieldSpec::required("statistical", WireType::Variant),
        FieldSpec::required("model-cards", WireType::CanonicalSet),
        FieldSpec::required("model-assumptions", WireType::CanonicalSet),
        FieldSpec::child_of("model-validity", &CERTIFIED_F64_VALIDITY_CHILD_V1),
        FieldSpec::required("model-discrepancy-ieee754-bits", WireType::U64),
        FieldSpec::required("model-in-domain", WireType::Bool),
        FieldSpec::required("sensitivity", WireType::OrderedBytes),
        FieldSpec::required("legacy-adjoint-correlation-present", WireType::Bool),
    ];
}

/// Low-level schema-shaped identity for one certified-f64 semantic frame.
///
/// Only [`IdentifiedCertifiedF64EvidenceV1`] proves that the frame was built
/// from an opaque [`Certified<f64>`] and helper-validated validity child.
pub type CertifiedF64EvidenceIdV1 = SemanticId<CertifiedF64EvidenceIdentitySchemaV1>;

/// Low-level producer receipt for a certified-f64 semantic frame.
pub type CertifiedF64EvidenceReceiptV1 = IdentityReceipt<CertifiedF64EvidenceIdV1>;

/// A locally certified scalar record kept attached to its unanchored semantic
/// identity and helper-built validity child.
#[derive(Debug, Clone)]
pub struct IdentifiedCertifiedF64EvidenceV1 {
    certified: Certified<f64>,
    validity_receipt: ValidityDomainReceiptV1,
    receipt: CertifiedF64EvidenceReceiptV1,
}

impl IdentifiedCertifiedF64EvidenceV1 {
    /// Read-only certified scalar evidence committed by this projection.
    #[must_use]
    pub const fn certified(&self) -> &Certified<f64> {
        &self.certified
    }

    /// Typed semantic-projection identity.
    #[must_use]
    pub const fn id(&self) -> CertifiedF64EvidenceIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored semantic receipt.
    #[must_use]
    pub const fn receipt(&self) -> CertifiedF64EvidenceReceiptV1 {
        self.receipt
    }

    /// Typed normalized validity identity bound as the semantic child.
    #[must_use]
    pub const fn validity_id(&self) -> ValidityDomainIdV1 {
        self.validity_receipt.id()
    }

    /// Complete helper-built validity receipt.
    #[must_use]
    pub const fn validity_receipt(&self) -> ValidityDomainReceiptV1 {
        self.validity_receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover the certified record.
    #[must_use]
    pub fn into_certified(self) -> Certified<f64> {
        self.certified
    }
}

/// Canonical identity schema for one exact source artifact supplied while
/// migrating certified scalar evidence away from a legacy FNV correlation.
pub enum CertifiedF64SourceIdentitySchemaV1 {}

impl CanonicalSchema for CertifiedF64SourceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.certified-f64-source.v1";
    const NAME: &'static str = "certified-f64-source";
    const VERSION: u32 = CERTIFIED_F64_SOURCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact caller-retained source domain, source schema version, and canonical bytes supplied for a certified-f64 legacy-FNV crosswalk; no format truth, origin, custody, execution, adjoint validity, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source-domain", WireType::Utf8),
        FieldSpec::required("source-schema-version", WireType::U64),
        FieldSpec::required("canonical-source", WireType::Bytes),
    ];
}

/// Low-level schema-shaped identity for an exact certified-f64 source.
///
/// Direct encoder output does not prove the source-domain/version admission or
/// correspondence with retained bytes. Those invariants belong to
/// [`CertifiedF64SourceV1`].
pub type CertifiedF64SourceIdV1 = SourceId<CertifiedF64SourceIdentitySchemaV1>;

/// Low-level producer receipt for one exact certified-f64 source frame.
pub type CertifiedF64SourceReceiptV1 = IdentityReceipt<CertifiedF64SourceIdV1>;

/// Exact source metadata and bytes kept attached to their unanchored identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CertifiedF64SourceV1 {
    domain: String,
    schema_version: u32,
    canonical_bytes: Vec<u8>,
    receipt: CertifiedF64SourceReceiptV1,
}

impl CertifiedF64SourceV1 {
    /// Exact source-domain UTF-8 committed by this identity.
    #[must_use]
    pub fn domain(&self) -> &str {
        &self.domain
    }

    /// Nonzero caller-declared source-schema version committed by this identity.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact canonical source bytes committed by this identity.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Typed exact-source identity.
    #[must_use]
    pub const fn id(&self) -> CertifiedF64SourceIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored source receipt.
    #[must_use]
    pub const fn receipt(&self) -> CertifiedF64SourceReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the identity attachment and recover exact source inputs.
    #[must_use]
    pub fn into_parts(self) -> (String, u32, Vec<u8>) {
        (self.domain, self.schema_version, self.canonical_bytes)
    }
}

static SOURCED_CERTIFIED_F64_EVIDENCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<CertifiedF64EvidenceIdV1>();
static SOURCED_CERTIFIED_F64_SOURCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<CertifiedF64SourceIdV1>();

/// Canonical parent schema attaching exact producer and optional adjoint
/// sources to one opaque certified-f64 semantic child.
pub enum SourcedCertifiedF64EvidenceIdentitySchemaV1 {}

impl CanonicalSchema for SourcedCertifiedF64EvidenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.sourced-certified-f64-evidence.v1";
    const NAME: &'static str = "sourced-certified-f64-evidence";
    const VERSION: u32 = SOURCED_CERTIFIED_F64_EVIDENCE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "typed certified-f64 semantic child plus exact typed producer source and optional exact typed adjoint source after pinned legacy-FNV-1a-64 correlation checks; weak FNV values excluded; no origin, custody, execution, units, adjoint correctness, authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of(
            "certified-f64-evidence",
            &SOURCED_CERTIFIED_F64_EVIDENCE_CHILD_V1,
        ),
        FieldSpec::child_of("producer-source", &SOURCED_CERTIFIED_F64_SOURCE_CHILD_V1),
        FieldSpec::required("adjoint-source-present", WireType::Bool),
        FieldSpec::ordered_children_of("adjoint-source", &SOURCED_CERTIFIED_F64_SOURCE_CHILD_V1),
    ];
}

/// Low-level schema-shaped identity for one exact-source-bound certified child.
///
/// Only [`IdentifiedSourcedCertifiedF64EvidenceV1`] proves that retained exact
/// source bytes passed both legacy-correlation crosswalks.
pub type SourcedCertifiedF64EvidenceIdV1 = SemanticId<SourcedCertifiedF64EvidenceIdentitySchemaV1>;

/// Low-level producer receipt for one sourced certified-f64 parent frame.
pub type SourcedCertifiedF64EvidenceReceiptV1 = IdentityReceipt<SourcedCertifiedF64EvidenceIdV1>;

/// A certified-f64 semantic child kept attached to exact producer and optional
/// adjoint sources whose bytes reproduce its legacy correlation tokens.
#[derive(Debug, Clone)]
pub struct IdentifiedSourcedCertifiedF64EvidenceV1 {
    certified_evidence: IdentifiedCertifiedF64EvidenceV1,
    producer_source: CertifiedF64SourceV1,
    adjoint_source: Option<CertifiedF64SourceV1>,
    receipt: SourcedCertifiedF64EvidenceReceiptV1,
}

impl IdentifiedSourcedCertifiedF64EvidenceV1 {
    /// Opaque certified semantic child committed by this parent.
    #[must_use]
    pub const fn certified_evidence(&self) -> &IdentifiedCertifiedF64EvidenceV1 {
        &self.certified_evidence
    }

    /// Exact producer source attached by the successful legacy crosswalk.
    #[must_use]
    pub const fn producer_source(&self) -> &CertifiedF64SourceV1 {
        &self.producer_source
    }

    /// Exact adjoint source attached when the certified child carries an
    /// adjoint correlation.
    #[must_use]
    pub const fn adjoint_source(&self) -> Option<&CertifiedF64SourceV1> {
        self.adjoint_source.as_ref()
    }

    /// Typed exact-source-bound parent identity.
    #[must_use]
    pub const fn id(&self) -> SourcedCertifiedF64EvidenceIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored parent receipt.
    #[must_use]
    pub const fn receipt(&self) -> SourcedCertifiedF64EvidenceReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the parent attachment without discarding either exact source.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        IdentifiedCertifiedF64EvidenceV1,
        CertifiedF64SourceV1,
        Option<CertifiedF64SourceV1>,
    ) {
        (
            self.certified_evidence,
            self.producer_source,
            self.adjoint_source,
        )
    }
}

static CERTIFIED_F64_LINEAGE_SOURCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<CertifiedF64SourceIdV1>();

/// Canonical graph schema for exact certified-f64 provenance lineage.
///
/// The graph binds exact source identities or an ordered binary operation over
/// typed parent nodes. The legacy FNV compatibility value is deliberately not
/// a field: it cannot become authority merely by being carried beside this
/// strong graph.
pub enum CertifiedF64LineageNodeIdentitySchemaV1 {}

impl CanonicalSchema for CertifiedF64LineageNodeIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.certified-f64-lineage-node.v1";
    const NAME: &'static str = "certified-f64-lineage-node";
    const VERSION: u32 = CERTIFIED_F64_LINEAGE_NODE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "exact typed source or versioned operation over exactly ordered typed parent lineage nodes; legacy FNV compatibility values excluded; no evidence-output semantics, execution, units, origin, custody, scientific authority, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("node-kind", WireType::Variant),
        FieldSpec::required("operation", WireType::Variant),
        FieldSpec::required("parent-semantics", WireType::Variant),
        FieldSpec::required("lineage-algorithm-version", WireType::U64),
        FieldSpec::ordered_children_of("source", &CERTIFIED_F64_LINEAGE_SOURCE_CHILD_V1),
        // A self-recursive ChildSpec would make this static schema recursive.
        // The opaque builder accepts only CertifiedF64LineageNodeIdV1 values
        // and frames their complete role/schema/root descriptors here.
        FieldSpec::required("parents", WireType::OrderedBytes),
    ];
}

/// Low-level schema-shaped identity for one certified-f64 lineage node.
///
/// Direct encoder output proves only canonical framing. Source/parent arity,
/// exact-source retention, operation mapping, ordered-parent preservation, and
/// compatibility crosswalk invariants belong to [`CertifiedF64LineageNodeV1`].
pub type CertifiedF64LineageNodeIdV1 = EvidenceNodeId<CertifiedF64LineageNodeIdentitySchemaV1>;

/// Low-level producer receipt for one certified-f64 lineage frame.
pub type CertifiedF64LineageNodeReceiptV1 = IdentityReceipt<CertifiedF64LineageNodeIdV1>;

/// Whether a lineage node is an exact source or an ordered composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CertifiedF64LineageNodeKindV1 {
    /// One retained exact-source child and no parents.
    Source,
    /// No source child and exactly two ordered typed parents.
    Composition,
}

impl CertifiedF64LineageNodeKindV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Source => 1,
            Self::Composition => 2,
        }
    }
}

/// Stable operation vocabulary committed by lineage schema v1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CertifiedF64LineageOperationV1 {
    /// Root a node at one exact retained source.
    Source,
    /// Ordered addition derivation.
    Add,
    /// Ordered subtraction derivation.
    Sub,
    /// Ordered multiplication derivation.
    Mul,
    /// Ordered minimum derivation.
    Min,
    /// Ordered maximum derivation.
    Max,
}

impl CertifiedF64LineageOperationV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Source => 1,
            Self::Add => 2,
            Self::Sub => 3,
            Self::Mul => 4,
            Self::Min => 5,
            Self::Max => 6,
        }
    }

    const fn kind(self) -> CertifiedF64LineageNodeKindV1 {
        match self {
            Self::Source => CertifiedF64LineageNodeKindV1::Source,
            Self::Add | Self::Sub | Self::Mul | Self::Min | Self::Max => {
                CertifiedF64LineageNodeKindV1::Composition
            }
        }
    }

    const fn parent_semantics(self) -> CertifiedF64LineageParentSemanticsV1 {
        CertifiedF64LineageParentSemanticsV1::Ordered
    }
}

/// The five operations supported by [`crate::Evidence::combine`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CertifiedF64LineageCompositionOpV1 {
    /// Addition.
    Add,
    /// Subtraction.
    Sub,
    /// Multiplication.
    Mul,
    /// Pointwise minimum.
    Min,
    /// Pointwise maximum.
    Max,
}

impl CertifiedF64LineageCompositionOpV1 {
    const fn node_operation(self) -> CertifiedF64LineageOperationV1 {
        match self {
            Self::Add => CertifiedF64LineageOperationV1::Add,
            Self::Sub => CertifiedF64LineageOperationV1::Sub,
            Self::Mul => CertifiedF64LineageOperationV1::Mul,
            Self::Min => CertifiedF64LineageOperationV1::Min,
            Self::Max => CertifiedF64LineageOperationV1::Max,
        }
    }

    const fn legacy_operation(self) -> Op {
        match self {
            Self::Add => Op::Add,
            Self::Sub => Op::Sub,
            Self::Mul => Op::Mul,
            Self::Min => Op::Min,
            Self::Max => Op::Max,
        }
    }
}

/// Parent ordering law for lineage schema v1.
///
/// All operations are ordered, including mathematically commutative ones,
/// because this schema mirrors the current ordered legacy chain and preserves
/// evaluation order, multiplicity, and association exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CertifiedF64LineageParentSemanticsV1 {
    /// Preserve the caller's left/right order without sorting or deduplication.
    Ordered,
}

impl CertifiedF64LineageParentSemanticsV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Ordered => 1,
        }
    }
}

/// An exact source or ordered derivation kept attached to its strong receipt.
///
/// The legacy provenance is retained only to interoperate with current
/// `Evidence` constructors. It is recomputed by the helper and excluded from
/// the strong root. Fields remain private so callers cannot detach exact source
/// bytes or typed parent IDs from the receipt that commits them.
#[derive(Debug, Clone)]
pub struct CertifiedF64LineageNodeV1 {
    source: Option<CertifiedF64SourceV1>,
    parents: Option<[CertifiedF64LineageNodeIdV1; 2]>,
    legacy_provenance: ProvenanceHash,
    receipt: CertifiedF64LineageNodeReceiptV1,
    operation: CertifiedF64LineageOperationV1,
}

impl CertifiedF64LineageNodeV1 {
    /// Exact source retained by a source node; absent for compositions.
    #[must_use]
    pub const fn source(&self) -> Option<&CertifiedF64SourceV1> {
        self.source.as_ref()
    }

    /// Ordered typed parent IDs; empty for a source node.
    #[must_use]
    pub fn parent_ids(&self) -> &[CertifiedF64LineageNodeIdV1] {
        self.parents
            .as_ref()
            .map_or(&[], |parents| parents.as_slice())
    }

    /// Weak compatibility correlation recomputed from exact source bytes or
    /// the current ordered `ProvenanceHash::chain` law. Never use as authority.
    #[must_use]
    pub const fn legacy_provenance(&self) -> ProvenanceHash {
        self.legacy_provenance
    }

    /// Typed strong lineage-node identity.
    #[must_use]
    pub const fn id(&self) -> CertifiedF64LineageNodeIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> CertifiedF64LineageNodeReceiptV1 {
        self.receipt
    }

    /// Stable operation committed by this node.
    #[must_use]
    pub const fn operation(&self) -> CertifiedF64LineageOperationV1 {
        self.operation
    }

    /// Source or composition node kind.
    #[must_use]
    pub const fn kind(&self) -> CertifiedF64LineageNodeKindV1 {
        self.operation.kind()
    }

    /// Exact ordered-parent law committed by schema v1.
    #[must_use]
    pub const fn parent_semantics(&self) -> CertifiedF64LineageParentSemanticsV1 {
        self.operation.parent_semantics()
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }
}

static DECISION_ASSESSMENT_CERTIFIED_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<CertifiedF64EvidenceIdV1>();

/// Canonical semantic schema for one decision assessment derived from an
/// opaque certified-f64 identity under the current local uncertainty law.
pub enum CertifiedF64DecisionAssessmentIdentitySchemaV1 {}

impl CanonicalSchema for CertifiedF64DecisionAssessmentIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.certified-f64-decision-assessment.v1";
    const NAME: &'static str = "certified-f64-decision-assessment";
    const VERSION: u32 = CERTIFIED_F64_DECISION_ASSESSMENT_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "typed certified-f64 child, assessment-algorithm version, exact threshold and derived uncertainty bits, status, and advice; no units, loss or decision context, policy authority, governor execution, scientific authority, origin, or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of(
            "certified-f64-evidence",
            &DECISION_ASSESSMENT_CERTIFIED_CHILD_V1,
        ),
        FieldSpec::required("assessment-algorithm-version", WireType::U64),
        FieldSpec::required("threshold-rel", WireType::FiniteF64),
        FieldSpec::required("numerical-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("statistical-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("model-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("total-rel-ieee754-bits", WireType::U64),
        FieldSpec::required("status", WireType::Variant),
        FieldSpec::required("advice", WireType::Variant),
    ];
}

/// Low-level schema-shaped identity for one certified-f64 decision assessment.
///
/// Only [`IdentifiedCertifiedF64DecisionAssessmentV1`] proves that all derived
/// fields were recomputed from the retained opaque child and threshold.
pub type CertifiedF64DecisionAssessmentIdV1 =
    SemanticId<CertifiedF64DecisionAssessmentIdentitySchemaV1>;

/// Low-level producer receipt for one certified-f64 decision assessment.
pub type CertifiedF64DecisionAssessmentReceiptV1 =
    IdentityReceipt<CertifiedF64DecisionAssessmentIdV1>;

/// A certified-f64 child and threshold kept attached to their recomputed local
/// uncertainty assessment and unanchored semantic identity.
#[derive(Debug, Clone)]
pub struct IdentifiedCertifiedF64DecisionAssessmentV1 {
    certified_evidence: IdentifiedCertifiedF64EvidenceV1,
    threshold_rel: f64,
    breakdown: UncertaintyBreakdown,
    total_rel: f64,
    status: DecisionStatus,
    advice: EscalationAdvice,
    receipt: CertifiedF64DecisionAssessmentReceiptV1,
}

impl IdentifiedCertifiedF64DecisionAssessmentV1 {
    /// Opaque certified evidence attachment used to recompute this assessment.
    #[must_use]
    pub const fn certified_evidence(&self) -> &IdentifiedCertifiedF64EvidenceV1 {
        &self.certified_evidence
    }

    /// Typed child identity bound by the assessment frame.
    #[must_use]
    pub const fn certified_evidence_id(&self) -> CertifiedF64EvidenceIdV1 {
        self.certified_evidence.id()
    }

    /// Exact caller-supplied relative decision threshold.
    #[must_use]
    pub const fn threshold_rel(&self) -> f64 {
        self.threshold_rel
    }

    /// Recomputed per-source relative uncertainty bands.
    #[must_use]
    pub const fn breakdown(&self) -> UncertaintyBreakdown {
        self.breakdown
    }

    /// Recomputed first-order total relative band.
    #[must_use]
    pub const fn total_rel(&self) -> f64 {
        self.total_rel
    }

    /// Recomputed local decision status, including presentation detail.
    ///
    /// The detail string is derived presentation and is not framed separately.
    #[must_use]
    pub const fn status(&self) -> &DecisionStatus {
        &self.status
    }

    /// Recomputed local escalation advice.
    #[must_use]
    pub const fn advice(&self) -> EscalationAdvice {
        self.advice
    }

    /// Typed assessment identity.
    #[must_use]
    pub const fn id(&self) -> CertifiedF64DecisionAssessmentIdV1 {
        self.receipt.id()
    }

    /// Complete unanchored producer receipt.
    #[must_use]
    pub const fn receipt(&self) -> CertifiedF64DecisionAssessmentReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }

    /// Surrender the assessment attachment and recover its semantic inputs.
    #[must_use]
    pub fn into_inputs(self) -> (IdentifiedCertifiedF64EvidenceV1, f64) {
        (self.certified_evidence, self.threshold_rel)
    }
}

static COLOR_EVIDENCE_SOURCE_CHILD_V1: ChildSpec =
    ChildSpec::for_identity::<ColorEvidenceSourceIdV1>();

/// Canonical identity schema for one color-evidence graph node.
pub enum ColorEvidenceNodeIdentitySchemaV1 {}

impl CanonicalSchema for ColorEvidenceNodeIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-evidence.color-evidence-node.v1";
    const NAME: &'static str = "color-evidence-node";
    const VERSION: u32 = COLOR_EVIDENCE_NODE_IDENTITY_VERSION_V1;
    const CONTEXT: &'static str = "node kind, operation law, color algebra, typed source, exact output color, and typed parent multiset or sequence";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("node-kind", WireType::Variant),
        FieldSpec::required("operation", WireType::Variant),
        FieldSpec::required("parent-semantics", WireType::Variant),
        FieldSpec::required("color-algebra-version", WireType::U64),
        FieldSpec::ordered_children_of("source", &COLOR_EVIDENCE_SOURCE_CHILD_V1),
        FieldSpec::required("output-color", WireType::Bytes),
        // A self-recursive ChildSpec would make the static schema recursive.
        // The public builder accepts only ColorEvidenceNodeIdV1 values and
        // frames their exact 32-byte roots here.
        FieldSpec::required("parents", WireType::OrderedBytes),
    ];
}

/// Low-level canonical-frame identity for one color-evidence graph node.
///
/// Direct generic encoder output proves only schema-shaped framing. Operation,
/// arity, source, parent-row, and recomputation invariants belong to the opaque
/// [`ColorEvidenceNodeV1`] returned by this module's helpers.
pub type ColorEvidenceNodeIdV1 = EvidenceNodeId<ColorEvidenceNodeIdentitySchemaV1>;

/// Low-level producer receipt underlying an opaque validated source.
pub type ColorEvidenceSourceReceiptV1 = IdentityReceipt<ColorEvidenceSourceIdV1>;
/// Low-level producer receipt underlying an opaque validated graph node.
pub type ColorEvidenceNodeReceiptV1 = IdentityReceipt<ColorEvidenceNodeIdV1>;

/// Unanchored canonical receipt for one exact retained source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColorEvidenceSourceV1 {
    receipt: ColorEvidenceSourceReceiptV1,
}

impl ColorEvidenceSourceV1 {
    /// Typed source identity.
    #[must_use]
    pub const fn id(&self) -> ColorEvidenceSourceIdV1 {
        self.receipt.id()
    }

    /// Complete producer receipt for downstream observation or authority work.
    #[must_use]
    pub const fn receipt(&self) -> ColorEvidenceSourceReceiptV1 {
        self.receipt
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }
}

/// A color plus its exact typed graph-node receipt.
///
/// Fields are private so a parent ID cannot be detached from the color whose
/// canonical bytes it commits. Construction is source-rooted or recomputed by
/// the v2 color algebra; neither route adds external trust.
#[derive(Debug, Clone, PartialEq)]
pub struct ColorEvidenceNodeV1 {
    color: Color,
    receipt: ColorEvidenceNodeReceiptV1,
    operation: ColorEvidenceOperationV1,
}

impl ColorEvidenceNodeV1 {
    /// Exact epistemic color committed by this node.
    #[must_use]
    pub const fn color(&self) -> &Color {
        &self.color
    }

    /// Typed graph-node identity.
    #[must_use]
    pub const fn id(&self) -> ColorEvidenceNodeIdV1 {
        self.receipt.id()
    }

    /// Complete producer receipt for downstream observation or authority work.
    #[must_use]
    pub const fn receipt(&self) -> ColorEvidenceNodeReceiptV1 {
        self.receipt
    }

    /// Stable operation committed by the node.
    #[must_use]
    pub const fn operation(&self) -> ColorEvidenceOperationV1 {
        self.operation
    }

    /// Source or derived node kind.
    #[must_use]
    pub const fn kind(&self) -> ColorEvidenceNodeKindV1 {
        self.operation.kind()
    }

    /// Ordered or commutative-multiset parent law.
    #[must_use]
    pub const fn parent_semantics(&self) -> ColorEvidenceParentSemanticsV1 {
        self.operation.parent_semantics()
    }

    /// Fixed-size typed digest bytes.
    #[must_use]
    pub fn id_bytes(&self) -> [u8; 32] {
        *self.id().as_bytes()
    }

    /// Identity state of a producer receipt. This is always unanchored.
    #[must_use]
    pub fn trust_state(&self) -> EvidenceIdentityTrustState {
        self.receipt.audit_record().trust()
    }
}

/// Whether the node is an independently retained source or a derivation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorEvidenceNodeKindV1 {
    /// A source node with one typed retained-source identity and no parents.
    Source,
    /// A derived node with typed parent-node identities and no source slot.
    Composition,
}

impl ColorEvidenceNodeKindV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Source => 1,
            Self::Composition => 2,
        }
    }
}

/// Stable operation vocabulary for color-evidence graph identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorEvidenceOperationV1 {
    /// Root a node at exact retained source bytes.
    Source,
    /// Addition.
    Add,
    /// Multiplication.
    Mul,
    /// Conservative interval hull.
    Hull,
}

impl ColorEvidenceOperationV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Source => 1,
            Self::Add => 2,
            Self::Mul => 3,
            Self::Hull => 4,
        }
    }

    const fn kind(self) -> ColorEvidenceNodeKindV1 {
        match self {
            Self::Source => ColorEvidenceNodeKindV1::Source,
            Self::Add | Self::Mul | Self::Hull => ColorEvidenceNodeKindV1::Composition,
        }
    }

    const fn parent_semantics(self) -> ColorEvidenceParentSemanticsV1 {
        match self {
            Self::Source => ColorEvidenceParentSemanticsV1::Ordered,
            Self::Add | Self::Mul | Self::Hull => {
                ColorEvidenceParentSemanticsV1::CommutativeMultiset
            }
        }
    }
}

/// The three operations implemented by the current versioned color algebra.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorEvidenceCompositionOpV1 {
    /// Addition.
    Add,
    /// Multiplication.
    Mul,
    /// Conservative interval hull.
    Hull,
}

impl ColorEvidenceCompositionOpV1 {
    const fn node_operation(self) -> ColorEvidenceOperationV1 {
        match self {
            Self::Add => ColorEvidenceOperationV1::Add,
            Self::Mul => ColorEvidenceOperationV1::Mul,
            Self::Hull => ColorEvidenceOperationV1::Hull,
        }
    }

    const fn interval_operation(self) -> IntervalOp {
        match self {
            Self::Add => IntervalOp::Add,
            Self::Mul => IntervalOp::Mul,
            Self::Hull => IntervalOp::Hull,
        }
    }
}

/// Whether parent order is semantic for this operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ColorEvidenceParentSemanticsV1 {
    /// Preserve caller order exactly.
    Ordered,
    /// Sort full typed parent roots lexicographically while preserving
    /// duplicates. This is a multiset, not a set.
    CommutativeMultiset,
}

impl ColorEvidenceParentSemanticsV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Ordered => 1,
            Self::CommutativeMultiset => 2,
        }
    }
}

/// Fail-closed refusal from color-evidence identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColorEvidenceIdentityError {
    /// A source schema domain must be explicit and non-empty.
    EmptySourceDomain,
    /// Source schema version zero is reserved for unknown/legacy data.
    ZeroSourceSchemaVersion,
    /// The output color is structurally malformed.
    MalformedColor(ColorPayloadError),
    /// Two parents presented the same typed ID with different retained-byte
    /// observations. Neither observation wins.
    ParentObservationConflict,
    /// The bounded canonical color buffer could not reserve its exact size.
    ColorBufferAllocationFailed {
        /// Exact preflighted payload bytes requested from the allocator.
        requested_bytes: u64,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

/// Fail-closed refusal from normalized validity-domain identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidityDomainIdentityError {
    /// A sorted axis has unusable interval bounds.
    InvalidBounds {
        /// Zero-based axis position in canonical `BTreeMap` order.
        axis_index: u64,
        /// Finite-ordering refusal.
        reason: &'static str,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

/// Fail-closed refusal from model-evidence semantic identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelEvidenceIdentityError {
    /// The relative discrepancy is NaN or negative.
    InvalidDiscrepancy {
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// The typed validity child refused normalization, limits, or cancellation.
    Validity(ValidityDomainIdentityError),
    /// Canonical framing, set admission, resources, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for ModelEvidenceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDiscrepancy { bits, reason } => write!(
                formatter,
                "model-evidence identity refused discrepancy bits 0x{bits:016x}: {reason}"
            ),
            Self::Validity(error) => {
                write!(
                    formatter,
                    "model-evidence identity refused validity: {error}"
                )
            }
            Self::Canonical(error) => write!(formatter, "model-evidence identity refused: {error}"),
        }
    }
}

impl std::error::Error for ModelEvidenceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Validity(error) => Some(error),
            Self::Canonical(error) => Some(error),
            Self::InvalidDiscrepancy { .. } => None,
        }
    }
}

impl From<ValidityDomainIdentityError> for ModelEvidenceIdentityError {
    fn from(error: ValidityDomainIdentityError) -> Self {
        Self::Validity(error)
    }
}

impl From<CanonicalError> for ModelEvidenceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from standalone numerical-certificate identity
/// construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NumericalCertificateIdentityError {
    /// The public certificate does not have the canonical structural shape for
    /// its declared kind.
    InvalidShape {
        /// Caller-carried certificate kind.
        kind: NumericalKind,
        /// Exact refused lower-bound bits.
        lo_bits: u64,
        /// Exact refused upper-bound bits.
        hi_bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for NumericalCertificateIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidShape {
                kind,
                lo_bits,
                hi_bits,
                reason,
            } => write!(
                formatter,
                "numerical-certificate identity refused {kind:?} bounds 0x{lo_bits:016x}..0x{hi_bits:016x}: {reason}"
            ),
            Self::Canonical(error) => {
                write!(formatter, "numerical-certificate identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for NumericalCertificateIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::InvalidShape { .. } => None,
        }
    }
}

impl From<CanonicalError> for NumericalCertificateIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from standalone statistical-certificate identity
/// construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StatisticalCertificateIdentityError {
    /// One public numeric parameter falls outside the local structural domain.
    InvalidParameter {
        /// Stable parameter name from `StatisticalCertificate` validation.
        field: &'static str,
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for StatisticalCertificateIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter {
                field,
                bits,
                reason,
            } => write!(
                formatter,
                "statistical-certificate identity refused {field} bits 0x{bits:016x}: {reason}"
            ),
            Self::Canonical(error) => {
                write!(
                    formatter,
                    "statistical-certificate identity refused: {error}"
                )
            }
        }
    }
}

impl std::error::Error for StatisticalCertificateIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::InvalidParameter { .. } => None,
        }
    }
}

impl From<CanonicalError> for StatisticalCertificateIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from standalone fidelity-pair identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FidelityPairIdentityError {
    /// A fit observation must locate the paired QoIs at a nonempty point.
    EmptyParameterPoint,
    /// One parameter name cannot enter the discrepancy-fit identity grammar.
    InvalidParameterName {
        /// Zero-based position in deterministic `BTreeMap` order.
        parameter_index: u64,
        /// Shared identity-grammar rejection reason.
        reason: &'static str,
    },
    /// One parameter coordinate is NaN or infinite.
    NonFiniteParameter {
        /// Zero-based position in deterministic `BTreeMap` order.
        parameter_index: u64,
        /// Exact refused IEEE-754 bits.
        bits: u64,
    },
    /// A low- or high-fidelity QoI is NaN or infinite.
    NonFiniteQoi {
        /// Stable field name: `lo_fi` or `hi_fi`.
        field: &'static str,
        /// Exact refused IEEE-754 bits.
        bits: u64,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for FidelityPairIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyParameterPoint => {
                formatter.write_str("fidelity-pair identity refused an empty parameter point")
            }
            Self::InvalidParameterName {
                parameter_index,
                reason,
            } => write!(
                formatter,
                "fidelity-pair identity refused parameter {parameter_index}: invalid name ({reason})"
            ),
            Self::NonFiniteParameter {
                parameter_index,
                bits,
            } => write!(
                formatter,
                "fidelity-pair identity refused parameter {parameter_index} value bits 0x{bits:016x}: coordinate must be finite"
            ),
            Self::NonFiniteQoi { field, bits } => write!(
                formatter,
                "fidelity-pair identity refused {field} QoI bits 0x{bits:016x}: QoI must be finite"
            ),
            Self::Canonical(error) => {
                write!(formatter, "fidelity-pair identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for FidelityPairIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::EmptyParameterPoint
            | Self::InvalidParameterName { .. }
            | Self::NonFiniteParameter { .. }
            | Self::NonFiniteQoi { .. } => None,
        }
    }
}

impl From<CanonicalError> for FidelityPairIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from standalone discrepancy-band identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiscrepancyBandIdentityError {
    /// A public band field is NaN, negative, or violates mean <= maximum.
    InvalidBand {
        /// Stable public field name.
        field: &'static str,
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for DiscrepancyBandIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBand {
                field,
                bits,
                reason,
            } => write!(
                formatter,
                "discrepancy-band identity refused {field} bits 0x{bits:016x}: {reason}"
            ),
            Self::Canonical(error) => {
                write!(formatter, "discrepancy-band identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for DiscrepancyBandIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::InvalidBand { .. } => None,
        }
    }
}

impl From<CanonicalError> for DiscrepancyBandIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from model-bracket semantic identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelBracketIdentityError {
    /// A bracket needs at least two distinct models to describe model choice.
    TooFewMembers {
        /// Admitted member count.
        count: u64,
        /// Minimum meaningful bracket size.
        minimum: u64,
    },
    /// One member name cannot enter the model-leaf identity grammar.
    InvalidMemberName {
        /// Zero-based position in canonical model-name order.
        member_index: u64,
        /// Shared identity-grammar rejection reason.
        reason: &'static str,
    },
    /// Internal bracket rows were not in strict model-name order.
    NonCanonicalMemberOrder {
        /// First row that was not greater than its predecessor.
        member_index: u64,
    },
    /// One member QoI is NaN or infinite.
    NonFiniteMemberQoi {
        /// Zero-based position in canonical model-name order.
        member_index: u64,
        /// Exact refused IEEE-754 bits.
        bits: u64,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for ModelBracketIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewMembers { count, minimum } => write!(
                formatter,
                "model-bracket identity refused {count} member(s): at least {minimum} are required"
            ),
            Self::InvalidMemberName {
                member_index,
                reason,
            } => write!(
                formatter,
                "model-bracket identity refused member {member_index}: invalid model name ({reason})"
            ),
            Self::NonCanonicalMemberOrder { member_index } => write!(
                formatter,
                "model-bracket identity refused member {member_index}: names must be strictly increasing"
            ),
            Self::NonFiniteMemberQoi { member_index, bits } => write!(
                formatter,
                "model-bracket identity refused member {member_index} QoI bits 0x{bits:016x}: QoI must be finite"
            ),
            Self::Canonical(error) => {
                write!(formatter, "model-bracket identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for ModelBracketIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::TooFewMembers { .. }
            | Self::InvalidMemberName { .. }
            | Self::NonCanonicalMemberOrder { .. }
            | Self::NonFiniteMemberQoi { .. } => None,
        }
    }
}

impl From<CanonicalError> for ModelBracketIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from model-card identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelCardIdentityError {
    /// The headline discrepancy is NaN or negative.
    InvalidDiscrepancy {
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// The legacy card and supplied exact-byte source disagree on whether a
    /// calibration artifact exists.
    CalibrationPresenceMismatch {
        /// Whether the legacy card carries a calibration correlation token.
        declared: bool,
        /// Whether exact calibration bytes were supplied.
        supplied: bool,
    },
    /// Exact supplied bytes do not reproduce the retained legacy correlation.
    CalibrationCorrelationMismatch {
        /// Legacy FNV value retained by the card.
        declared: u64,
        /// Incrementally recomputed FNV value over the supplied bytes.
        computed: u64,
    },
    /// Cancellation was observed during the incremental legacy crosswalk.
    CalibrationCrosswalkCancelled {
        /// Exact calibration-source bytes processed before observation.
        processed_bytes: u64,
    },
    /// The typed validity child refused normalization, limits, or cancellation.
    Validity(ValidityDomainIdentityError),
    /// Canonical framing, set admission, resources, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for ModelCardIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDiscrepancy { bits, reason } => write!(
                formatter,
                "model-card identity refused discrepancy bits 0x{bits:016x}: {reason}"
            ),
            Self::CalibrationPresenceMismatch { declared, supplied } => write!(
                formatter,
                "model-card identity refused calibration shape: legacy card presence is {declared}, exact-byte source presence is {supplied}"
            ),
            Self::CalibrationCorrelationMismatch { declared, computed } => write!(
                formatter,
                "model-card identity refused calibration crosswalk: legacy FNV 0x{declared:016x} does not match exact-byte FNV 0x{computed:016x}"
            ),
            Self::CalibrationCrosswalkCancelled { processed_bytes } => write!(
                formatter,
                "model-card identity cancelled during calibration crosswalk after {processed_bytes} source bytes"
            ),
            Self::Validity(error) => {
                write!(formatter, "model-card identity refused validity: {error}")
            }
            Self::Canonical(error) => write!(formatter, "model-card identity refused: {error}"),
        }
    }
}

impl std::error::Error for ModelCardIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Validity(error) => Some(error),
            Self::Canonical(error) => Some(error),
            Self::InvalidDiscrepancy { .. }
            | Self::CalibrationPresenceMismatch { .. }
            | Self::CalibrationCorrelationMismatch { .. }
            | Self::CalibrationCrosswalkCancelled { .. } => None,
        }
    }
}

impl From<ValidityDomainIdentityError> for ModelCardIdentityError {
    fn from(error: ValidityDomainIdentityError) -> Self {
        Self::Validity(error)
    }
}

impl From<CanonicalError> for ModelCardIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from card-bound model-evidence construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CardBoundModelEvidenceIdentityError {
    /// The canonical card-name set and supplied opaque declarations differ in
    /// cardinality.
    CardCountMismatch {
        /// Number of names committed by the model-evidence child.
        declared: u64,
        /// Number of opaque card declarations supplied for binding.
        supplied: u64,
    },
    /// One supplied card name differs byte-for-byte from its declared name.
    CardNameMismatch {
        /// Zero-based position in canonical card-name order.
        card_index: u64,
        /// Declared model-evidence name length in bytes.
        declared_bytes: u64,
        /// Supplied model-card name length in bytes.
        supplied_bytes: u64,
        /// First differing byte, or the common-prefix length when lengths differ.
        first_different_byte: u64,
    },
    /// Cancellation was observed while comparing one potentially long name.
    CardNameComparisonCancelled {
        /// Zero-based position in canonical card-name order.
        card_index: u64,
        /// Common-prefix bytes compared before cancellation observation.
        compared_bytes: u64,
    },
    /// Parent framing, resources, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for CardBoundModelEvidenceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CardCountMismatch { declared, supplied } => write!(
                formatter,
                "card-bound model-evidence identity refused {supplied} declaration(s): model evidence names {declared}"
            ),
            Self::CardNameMismatch {
                card_index,
                declared_bytes,
                supplied_bytes,
                first_different_byte,
            } => write!(
                formatter,
                "card-bound model-evidence identity refused card {card_index}: declared {declared_bytes} name bytes and supplied {supplied_bytes} name bytes first differ at byte {first_different_byte}"
            ),
            Self::CardNameComparisonCancelled {
                card_index,
                compared_bytes,
            } => write!(
                formatter,
                "card-bound model-evidence identity cancelled while comparing card {card_index} after {compared_bytes} name bytes"
            ),
            Self::Canonical(error) => {
                write!(
                    formatter,
                    "card-bound model-evidence identity refused: {error}"
                )
            }
        }
    }
}

impl std::error::Error for CardBoundModelEvidenceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::CardCountMismatch { .. }
            | Self::CardNameMismatch { .. }
            | Self::CardNameComparisonCancelled { .. } => None,
        }
    }
}

impl From<CanonicalError> for CardBoundModelEvidenceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from certified-f64 semantic-identity construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CertifiedF64EvidenceIdentityError {
    /// The typed validity child refused normalization, limits, or cancellation.
    Validity(ValidityDomainIdentityError),
    /// Outer semantic framing, set admission, resources, or cancellation
    /// refused.
    Canonical(CanonicalError),
}

impl fmt::Display for CertifiedF64EvidenceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Validity(error) => {
                write!(
                    formatter,
                    "certified-f64 identity refused validity: {error}"
                )
            }
            Self::Canonical(error) => {
                write!(formatter, "certified-f64 identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for CertifiedF64EvidenceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Validity(error) => Some(error),
            Self::Canonical(error) => Some(error),
        }
    }
}

impl From<ValidityDomainIdentityError> for CertifiedF64EvidenceIdentityError {
    fn from(error: ValidityDomainIdentityError) -> Self {
        Self::Validity(error)
    }
}

impl From<CanonicalError> for CertifiedF64EvidenceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from exact certified-f64 source construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CertifiedF64SourceIdentityError {
    /// Source domains are required identity-bearing interpretation metadata.
    EmptySourceDomain,
    /// Source schema version zero is reserved for invalid/absent metadata.
    ZeroSourceSchemaVersion,
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for CertifiedF64SourceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySourceDomain => {
                write!(
                    formatter,
                    "certified-f64 source identity refused an empty domain"
                )
            }
            Self::ZeroSourceSchemaVersion => write!(
                formatter,
                "certified-f64 source identity refused source schema version zero"
            ),
            Self::Canonical(error) => {
                write!(formatter, "certified-f64 source identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for CertifiedF64SourceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::EmptySourceDomain | Self::ZeroSourceSchemaVersion => None,
        }
    }
}

impl From<CanonicalError> for CertifiedF64SourceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from exact-source-bound certified-f64 construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourcedCertifiedF64EvidenceIdentityError {
    /// The carried adjoint correlation and supplied exact adjoint source
    /// disagree on presence.
    AdjointSourcePresenceMismatch {
        /// Whether the certified record carries a legacy adjoint correlation.
        declared: bool,
        /// Whether an exact adjoint source was supplied.
        supplied: bool,
    },
    /// Exact supplied bytes do not reproduce one retained legacy correlation.
    SourceCorrelationMismatch {
        /// Stable source role: `producer` or `adjoint`.
        source: &'static str,
        /// Legacy FNV value retained by the certified record.
        declared: u64,
        /// Incrementally recomputed FNV value over exact source bytes.
        computed: u64,
    },
    /// Cancellation was observed during one incremental legacy crosswalk.
    SourceCrosswalkCancelled {
        /// Stable source role: `producer` or `adjoint`.
        source: &'static str,
        /// Exact source bytes processed before cancellation observation.
        processed_bytes: u64,
    },
    /// Parent framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for SourcedCertifiedF64EvidenceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AdjointSourcePresenceMismatch { declared, supplied } => write!(
                formatter,
                "sourced certified-f64 identity refused adjoint shape: legacy correlation presence is {declared}, exact-source presence is {supplied}"
            ),
            Self::SourceCorrelationMismatch {
                source,
                declared,
                computed,
            } => write!(
                formatter,
                "sourced certified-f64 identity refused {source} crosswalk: legacy FNV 0x{declared:016x} does not match exact-byte FNV 0x{computed:016x}"
            ),
            Self::SourceCrosswalkCancelled {
                source,
                processed_bytes,
            } => write!(
                formatter,
                "sourced certified-f64 identity cancelled during {source} crosswalk after {processed_bytes} source bytes"
            ),
            Self::Canonical(error) => {
                write!(formatter, "sourced certified-f64 identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for SourcedCertifiedF64EvidenceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::AdjointSourcePresenceMismatch { .. }
            | Self::SourceCorrelationMismatch { .. }
            | Self::SourceCrosswalkCancelled { .. } => None,
        }
    }
}

impl From<CanonicalError> for SourcedCertifiedF64EvidenceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

/// Fail-closed refusal from certified-f64 lineage construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CertifiedF64LineageIdentityError {
    /// Two parents presented the same typed identity with conflicting receipt
    /// observations. Neither observation wins.
    ParentObservationConflict,
    /// Cancellation was observed while deriving the weak compatibility token
    /// from exact source bytes.
    SourceCrosswalkCancelled {
        /// Exact source bytes processed before cancellation observation.
        processed_bytes: u64,
    },
    /// The fixed-size fallible compatibility-chain builder refused its byte
    /// budget, framing arithmetic, or allocation.
    LegacyChain(IdentityBuildError),
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for CertifiedF64LineageIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParentObservationConflict => write!(
                formatter,
                "certified-f64 lineage refused conflicting observations for one parent identity"
            ),
            Self::SourceCrosswalkCancelled { processed_bytes } => write!(
                formatter,
                "certified-f64 lineage cancelled during exact-source legacy crosswalk after {processed_bytes} bytes"
            ),
            Self::LegacyChain(error) => write!(
                formatter,
                "certified-f64 lineage legacy compatibility chain refused: {error}"
            ),
            Self::Canonical(error) => {
                write!(formatter, "certified-f64 lineage identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for CertifiedF64LineageIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::LegacyChain(error) => Some(error),
            Self::ParentObservationConflict | Self::SourceCrosswalkCancelled { .. } => None,
        }
    }
}

impl From<CanonicalError> for CertifiedF64LineageIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

impl From<IdentityBuildError> for CertifiedF64LineageIdentityError {
    fn from(error: IdentityBuildError) -> Self {
        Self::LegacyChain(error)
    }
}

/// Fail-closed refusal from certified-f64 decision-assessment identity
/// construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CertifiedF64DecisionAssessmentIdentityError {
    /// The caller-supplied decision threshold is NaN, infinite, or negative.
    InvalidThreshold {
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// A recomputed relative band is NaN or negative, indicating algorithm or
    /// invariant drift that must not receive an identity.
    InvalidDerivedBand {
        /// Stable derived-field name.
        field: &'static str,
        /// Exact refused IEEE-754 bits.
        bits: u64,
        /// Structural requirement that was violated.
        reason: &'static str,
    },
    /// Canonical framing, resource admission, or cancellation refused.
    Canonical(CanonicalError),
}

impl fmt::Display for CertifiedF64DecisionAssessmentIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidThreshold { bits, reason } => write!(
                formatter,
                "certified-f64 decision assessment refused threshold bits 0x{bits:016x}: {reason}"
            ),
            Self::InvalidDerivedBand {
                field,
                bits,
                reason,
            } => write!(
                formatter,
                "certified-f64 decision assessment refused derived {field} bits 0x{bits:016x}: {reason}"
            ),
            Self::Canonical(error) => write!(
                formatter,
                "certified-f64 decision assessment identity refused: {error}"
            ),
        }
    }
}

impl std::error::Error for CertifiedF64DecisionAssessmentIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::InvalidThreshold { .. } | Self::InvalidDerivedBand { .. } => None,
        }
    }
}

impl From<CanonicalError> for CertifiedF64DecisionAssessmentIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

impl fmt::Display for ValidityDomainIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBounds { axis_index, reason } => write!(
                formatter,
                "validity-domain identity refused bounds for axis {axis_index}: {reason}"
            ),
            Self::Canonical(error) => {
                write!(formatter, "validity-domain identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for ValidityDomainIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Canonical(error) => Some(error),
            Self::InvalidBounds { .. } => None,
        }
    }
}

impl From<CanonicalError> for ValidityDomainIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

impl fmt::Display for ColorEvidenceIdentityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySourceDomain => {
                formatter.write_str("color-evidence source schema domain must be non-empty")
            }
            Self::ZeroSourceSchemaVersion => formatter
                .write_str("color-evidence source schema version zero is reserved for legacy data"),
            Self::MalformedColor(error) => {
                write!(
                    formatter,
                    "color-evidence identity refused malformed output: {error}"
                )
            }
            Self::ParentObservationConflict => formatter.write_str(
                "color-evidence composition refused one typed parent ID backed by different byte observations",
            ),
            Self::ColorBufferAllocationFailed { requested_bytes } => write!(
                formatter,
                "color-evidence identity could not reserve its {requested_bytes}-byte canonical color buffer"
            ),
            Self::Canonical(error) => write!(formatter, "color-evidence identity refused: {error}"),
        }
    }
}

impl std::error::Error for ColorEvidenceIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::MalformedColor(error) => Some(error),
            Self::Canonical(error) => Some(error),
            Self::EmptySourceDomain
            | Self::ZeroSourceSchemaVersion
            | Self::ParentObservationConflict
            | Self::ColorBufferAllocationFailed { .. } => None,
        }
    }
}

impl From<CanonicalError> for ColorEvidenceIdentityError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

fn poll_identity_cancellation<C>(cancellation: &mut C) -> Result<(), CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if cancellation.is_cancelled() {
        Err(CanonicalError::Cancelled { absorbed_bytes: 0 })
    } else {
        Ok(())
    }
}

fn add_bounded_color_bytes(
    length: &mut u64,
    additional: u64,
    limit: u64,
) -> Result<(), ColorEvidenceIdentityError> {
    let requested = length
        .checked_add(additional)
        .ok_or(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LengthOverflow,
        ))?;
    if requested > limit {
        return Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested,
                limit,
            },
        ));
    }
    *length = requested;
    Ok(())
}

fn bounded_len(value: usize) -> Result<u64, CanonicalError> {
    u64::try_from(value).map_err(|_| CanonicalError::LengthOverflow)
}

fn poll_color_buffer_cancellation<C>(
    output: &[u8],
    cancellation: &mut C,
) -> Result<(), CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if cancellation.is_cancelled() {
        Err(CanonicalError::Cancelled {
            absorbed_bytes: bounded_len(output.len())?,
        })
    } else {
        Ok(())
    }
}

fn append_color_bytes<C>(
    output: &mut Vec<u8>,
    bytes: &[u8],
    cancellation_poll_bytes: usize,
    cancellation: &mut C,
) -> Result<(), CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    for chunk in bytes.chunks(cancellation_poll_bytes) {
        poll_color_buffer_cancellation(output, cancellation)?;
        output.extend_from_slice(chunk);
    }
    Ok(())
}

fn push_color_len<C>(
    output: &mut Vec<u8>,
    length: usize,
    cancellation_poll_bytes: usize,
    cancellation: &mut C,
) -> Result<(), CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    append_color_bytes(
        output,
        &bounded_len(length)?.to_le_bytes(),
        cancellation_poll_bytes,
        cancellation,
    )
}

fn push_color_field<C>(
    output: &mut Vec<u8>,
    bytes: &[u8],
    cancellation_poll_bytes: usize,
    cancellation: &mut C,
) -> Result<(), CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    push_color_len(output, bytes.len(), cancellation_poll_bytes, cancellation)?;
    append_color_bytes(output, bytes, cancellation_poll_bytes, cancellation)
}

/// Normalize and identify one evidence-validity domain.
///
/// The owned input is retained inside the opaque result, preventing the
/// admitted semantic identity from being detached from different bounds. Axis
/// rows use `BTreeMap` order and bind the axis byte length, exact UTF-8 bytes
/// without normalization, and both IEEE-754 endpoint bit patterns. An
/// unconstrained domain is the canonical empty row sequence; it is not an
/// invalid domain.
///
/// # Errors
/// Refuses non-finite or inverted bounds, invalid limits, resource overflow, or
/// cancellation. No partial identity is published.
pub fn identify_validity_domain_v1<C>(
    domain: ValidityDomain,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedValidityDomainV1, ValidityDomainIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let receipt = identify_validity_domain_receipt_v1(&domain, limits, &mut cancellation)?;
    Ok(IdentifiedValidityDomainV1 { domain, receipt })
}

fn identify_validity_domain_receipt_v1<C>(
    domain: &ValidityDomain,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<ValidityDomainReceiptV1, ValidityDomainIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive"),
        ));
    }
    poll_identity_cancellation(cancellation)?;
    let axis_count = bounded_len(domain.bounds().len())?;
    if axis_count > limits.max_collection_items() {
        return Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: LimitKind::CollectionItems,
                requested: axis_count,
                limit: limits.max_collection_items(),
            },
        ));
    }

    let field_limit = limits
        .max_field_bytes()
        .min(MAX_VALIDITY_DOMAIN_FIELD_BYTES_V1);
    let mut field_payload_bytes = u64::from(u64::BITS / 8);
    if field_payload_bytes > field_limit {
        return Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: field_payload_bytes,
                limit: field_limit,
            },
        ));
    }
    for (axis_index, (axis, (lo, hi))) in domain.bounds().iter().enumerate() {
        poll_identity_cancellation(cancellation)?;
        let axis_index = bounded_len(axis_index)?;
        if !lo.is_finite() || !hi.is_finite() {
            return Err(ValidityDomainIdentityError::InvalidBounds {
                axis_index,
                reason: "bounds must be finite",
            });
        }
        if lo > hi {
            return Err(ValidityDomainIdentityError::InvalidBounds {
                axis_index,
                reason: "lower bound exceeds upper bound",
            });
        }
        let row_bytes = 24_u64
            .checked_add(bounded_len(axis.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
        let framed_row_bytes = u64::from(u64::BITS / 8)
            .checked_add(row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        field_payload_bytes = field_payload_bytes
            .checked_add(framed_row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        if field_payload_bytes > field_limit {
            return Err(ValidityDomainIdentityError::Canonical(
                CanonicalError::LimitExceeded {
                    kind: LimitKind::FieldBytes,
                    requested: field_payload_bytes,
                    limit: field_limit,
                },
            ));
        }
    }
    let required_stream_chunks = axis_count
        .checked_mul(VALIDITY_DOMAIN_STREAM_CHUNKS_PER_AXIS_V1)
        .ok_or(CanonicalError::LengthOverflow)?;
    if required_stream_chunks > limits.max_collection_items() {
        return Err(ValidityDomainIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: LimitKind::StreamChunks,
                requested: required_stream_chunks,
                limit: limits.max_collection_items(),
            },
        ));
    }

    let receipt = {
        let row_lengths = domain.bounds().keys().map(|axis| {
            bounded_len(axis.len()).and_then(|axis_bytes| {
                24_u64
                    .checked_add(axis_bytes)
                    .ok_or(CanonicalError::LengthOverflow)
            })
        });
        let mut rows = domain.bounds().iter();
        CanonicalEncoder::<ValidityDomainIdV1, _>::new(limits, || cancellation.is_cancelled())?
            .ordered_bytes_stream(
                Field::new(0, "axes"),
                axis_count,
                row_lengths,
                |row_index, mut sink| -> Result<(), CanonicalError> {
                    let Some((axis, (lo, hi))) = rows.next() else {
                        return Err(CanonicalError::DeclaredLengthMismatch {
                            declared: axis_count,
                            observed: row_index,
                        });
                    };
                    sink.write(&bounded_len(axis.len())?.to_le_bytes())?;
                    sink.write(axis.as_bytes())?;
                    sink.write(&lo.to_bits().to_le_bytes())?;
                    sink.write(&hi.to_bits().to_le_bytes())?;
                    Ok(())
                },
            )
            .map_err(|error| match error {
                OrderedBytesStreamError::Canonical { source, .. }
                | OrderedBytesStreamError::Producer { source, .. } => {
                    ValidityDomainIdentityError::Canonical(source)
                }
            })?
            .finish()?
    };
    Ok(receipt)
}

const fn numerical_certificate_kind_tag_v1(kind: NumericalKind) -> u32 {
    match kind {
        NumericalKind::Exact => 1,
        NumericalKind::Enclosure => 2,
        NumericalKind::Estimate => 3,
        NumericalKind::NoClaim => 4,
    }
}

fn validate_numerical_certificate_shape_v1(
    certificate: &crate::NumericalCertificate,
) -> Result<(), NumericalCertificateIdentityError> {
    let invalid = |reason| NumericalCertificateIdentityError::InvalidShape {
        kind: certificate.kind,
        lo_bits: certificate.lo.to_bits(),
        hi_bits: certificate.hi.to_bits(),
        reason,
    };
    if certificate.lo.is_nan() || certificate.hi.is_nan() {
        return Err(invalid("bounds must not be NaN"));
    }
    match certificate.kind {
        NumericalKind::Exact if certificate.lo.to_bits() != certificate.hi.to_bits() => {
            Err(invalid("Exact bounds must be bit-identical"))
        }
        NumericalKind::Enclosure | NumericalKind::Estimate if certificate.lo > certificate.hi => {
            Err(invalid("lower bound must not exceed upper bound"))
        }
        NumericalKind::NoClaim
            if certificate.lo.to_bits() != f64::NEG_INFINITY.to_bits()
                || certificate.hi.to_bits() != f64::INFINITY.to_bits() =>
        {
            Err(invalid(
                "NoClaim must use the canonical negative-infinity to positive-infinity bounds",
            ))
        }
        NumericalKind::Exact
        | NumericalKind::Enclosure
        | NumericalKind::Estimate
        | NumericalKind::NoClaim => Ok(()),
    }
}

/// Identify the exact admitted structural state of one standalone numerical
/// certificate.
///
/// The helper consumes and retains the public/mutable certificate without
/// normalization. It binds the stable kind tag and exact endpoint bits.
/// Enclosure and estimate declarations may use ordered infinite bounds, and an
/// exact declaration may use matching infinite bits; those are structural
/// states only, not local certification. `NoClaim` has one canonical
/// negative-infinity to positive-infinity representation.
///
/// # Errors
/// Refuses NaN, inverted bounds, non-identical exact endpoints, a forged
/// `NoClaim` shape, invalid limits, resource overflow, or cancellation. No
/// partial identity is published.
pub fn identify_numerical_certificate_v1<C>(
    certificate: crate::NumericalCertificate,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedNumericalCertificateV1, NumericalCertificateIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    validate_numerical_certificate_shape_v1(&certificate)?;
    let receipt = CanonicalEncoder::<NumericalCertificateIdV1, _>::new(limits, cancellation)?
        .variant(
            Field::new(0, "kind"),
            numerical_certificate_kind_tag_v1(certificate.kind),
            &[],
        )?
        .u64(Field::new(1, "lo-ieee754-bits"), certificate.lo.to_bits())?
        .u64(Field::new(2, "hi-ieee754-bits"), certificate.hi.to_bits())?
        .finish()?;
    Ok(IdentifiedNumericalCertificateV1 {
        certificate,
        receipt,
    })
}

/// Identify the exact admitted structural state of one standalone statistical
/// certificate.
///
/// The stable variant tag and exact accepted numeric bits are retained without
/// normalization. This is a local shape check only: the frame carries no null,
/// estimand, method, sample, stopping rule, dependence context, or coverage
/// proof. `None` binds only the caller's local no-stochastic-component state.
///
/// # Errors
/// Refuses a malformed public numeric parameter, invalid limits, resource
/// overflow, or cancellation. No partial identity is published.
pub fn identify_statistical_certificate_v1<C>(
    certificate: StatisticalCertificate,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedStatisticalCertificateV1, StatisticalCertificateIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    if let Some((field, value, reason)) = certificate.validation_issue() {
        return Err(StatisticalCertificateIdentityError::InvalidParameter {
            field,
            bits: value.to_bits(),
            reason,
        });
    }

    let mut payload = [0_u8; 16];
    let (tag, payload_len) = match certificate {
        StatisticalCertificate::None => (1, 0),
        StatisticalCertificate::EValue { e, alpha } => {
            payload[..8].copy_from_slice(&e.to_bits().to_le_bytes());
            payload[8..].copy_from_slice(&alpha.to_bits().to_le_bytes());
            (2, payload.len())
        }
        StatisticalCertificate::HalfWidth {
            half_width,
            confidence,
        } => {
            payload[..8].copy_from_slice(&half_width.to_bits().to_le_bytes());
            payload[8..].copy_from_slice(&confidence.to_bits().to_le_bytes());
            (3, payload.len())
        }
    };
    let receipt = CanonicalEncoder::<StatisticalCertificateIdV1, _>::new(limits, cancellation)?
        .variant(Field::new(0, "certificate"), tag, &payload[..payload_len])?
        .finish()?;
    Ok(IdentifiedStatisticalCertificateV1 {
        certificate,
        receipt,
    })
}

fn preflight_fidelity_pair_parameters_v1<C>(
    pair: &FidelityPair,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<u64, FidelityPairIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let parameter_count = bounded_len(pair.params.len())?;
    if parameter_count == 0 {
        return Err(FidelityPairIdentityError::EmptyParameterPoint);
    }
    let parameter_limit = limits
        .max_collection_items()
        .min(MAX_FIDELITY_PAIR_PARAMETERS_V1);
    if parameter_count > parameter_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: parameter_count,
            limit: parameter_limit,
        }
        .into());
    }

    let field_limit = limits
        .max_field_bytes()
        .min(MAX_FIDELITY_PAIR_PARAMETERS_FIELD_BYTES_V1);
    let mut field_payload_bytes = u64::from(u64::BITS / 8);
    if field_payload_bytes > field_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: field_payload_bytes,
            limit: field_limit,
        }
        .into());
    }
    for (parameter_index, (name, value)) in pair.params.iter().enumerate() {
        poll_identity_cancellation(cancellation)?;
        let parameter_index = bounded_len(parameter_index)?;
        if let Some(reason) = color_identity_reason(name) {
            return Err(FidelityPairIdentityError::InvalidParameterName {
                parameter_index,
                reason,
            });
        }
        if !value.is_finite() {
            return Err(FidelityPairIdentityError::NonFiniteParameter {
                parameter_index,
                bits: value.to_bits(),
            });
        }
        let row_bytes = 16_u64
            .checked_add(bounded_len(name.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
        let framed_row_bytes = u64::from(u64::BITS / 8)
            .checked_add(row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        field_payload_bytes = field_payload_bytes
            .checked_add(framed_row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        if field_payload_bytes > field_limit {
            return Err(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: field_payload_bytes,
                limit: field_limit,
            }
            .into());
        }
    }
    let required_stream_chunks = parameter_count
        .checked_mul(FIDELITY_PAIR_STREAM_CHUNKS_PER_PARAMETER_V1)
        .ok_or(CanonicalError::LengthOverflow)?;
    if required_stream_chunks > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::StreamChunks,
            requested: required_stream_chunks,
            limit: limits.max_collection_items(),
        }
        .into());
    }
    Ok(parameter_count)
}

/// Identify one exact, structurally admitted two-fidelity observation.
///
/// The helper consumes and retains the public pair while binding its
/// nonempty, `BTreeMap`-ordered parameter point and exact low/high QoI bits.
/// Parameter names follow the discrepancy fit's machine-readable identity
/// grammar. Every coordinate and QoI must be finite; accepted signed zero is
/// preserved without normalization.
///
/// # Errors
/// Refuses an empty or oversized point, an invalid parameter name, non-finite
/// coordinates or QoIs, invalid limits, resource overflow, or cancellation.
/// No partial identity is published.
pub fn identify_fidelity_pair_v1<C>(
    pair: FidelityPair,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedFidelityPairV1, FidelityPairIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    if !pair.lo_fi.is_finite() {
        return Err(FidelityPairIdentityError::NonFiniteQoi {
            field: "lo_fi",
            bits: pair.lo_fi.to_bits(),
        });
    }
    if !pair.hi_fi.is_finite() {
        return Err(FidelityPairIdentityError::NonFiniteQoi {
            field: "hi_fi",
            bits: pair.hi_fi.to_bits(),
        });
    }
    let parameter_count = preflight_fidelity_pair_parameters_v1(&pair, limits, &mut cancellation)?;
    let receipt = {
        let row_lengths = pair.params.keys().map(|name| {
            bounded_len(name.len()).and_then(|name_bytes| {
                16_u64
                    .checked_add(name_bytes)
                    .ok_or(CanonicalError::LengthOverflow)
            })
        });
        let mut rows = pair.params.iter();
        CanonicalEncoder::<FidelityPairIdV1, _>::new(limits, cancellation)?
            .ordered_bytes_stream(
                Field::new(0, "parameters"),
                parameter_count,
                row_lengths,
                |row_index, mut sink| -> Result<(), CanonicalError> {
                    let Some((name, value)) = rows.next() else {
                        return Err(CanonicalError::DeclaredLengthMismatch {
                            declared: parameter_count,
                            observed: row_index,
                        });
                    };
                    sink.write(&bounded_len(name.len())?.to_le_bytes())?;
                    sink.write(name.as_bytes())?;
                    sink.write(&value.to_bits().to_le_bytes())?;
                    Ok(())
                },
            )
            .map_err(|error| match error {
                OrderedBytesStreamError::Canonical { source, .. }
                | OrderedBytesStreamError::Producer { source, .. } => {
                    FidelityPairIdentityError::Canonical(source)
                }
            })?
            .u64(
                Field::new(1, "lo-fi-qoi-ieee754-bits"),
                pair.lo_fi.to_bits(),
            )?
            .u64(
                Field::new(2, "hi-fi-qoi-ieee754-bits"),
                pair.hi_fi.to_bits(),
            )?
            .finish()?
    };
    Ok(IdentifiedFidelityPairV1 { pair, receipt })
}

fn validate_discrepancy_band_v1(band: DiscrepancyBand) -> Result<(), DiscrepancyBandIdentityError> {
    let invalid = |field, value: f64, reason| DiscrepancyBandIdentityError::InvalidBand {
        field,
        bits: value.to_bits(),
        reason,
    };
    if band.mean_rel.is_nan() {
        return Err(invalid("mean_rel", band.mean_rel, "value must not be NaN"));
    }
    if band.max_rel.is_nan() {
        return Err(invalid("max_rel", band.max_rel, "value must not be NaN"));
    }
    if band.mean_rel < 0.0 {
        return Err(invalid(
            "mean_rel",
            band.mean_rel,
            "value must be non-negative",
        ));
    }
    if band.max_rel < 0.0 {
        return Err(invalid(
            "max_rel",
            band.max_rel,
            "value must be non-negative",
        ));
    }
    if band.mean_rel > band.max_rel {
        return Err(invalid(
            "mean_rel",
            band.mean_rel,
            "mean_rel must not exceed max_rel",
        ));
    }
    Ok(())
}

/// Identify the exact admitted structural state of one discrepancy band.
///
/// The helper consumes and retains both raw relative-discrepancy values without
/// normalization. Positive infinity is an explicit unbounded state and signed
/// zero remains bit-distinct. This frame carries no training corpus, validity
/// domain, pair count, query point, discrepancy definition, or derivation.
///
/// # Errors
/// Refuses NaN, negative values, mean greater than maximum, invalid limits,
/// resource overflow, or cancellation. No partial identity is published.
pub fn identify_discrepancy_band_v1<C>(
    band: DiscrepancyBand,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedDiscrepancyBandV1, DiscrepancyBandIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    validate_discrepancy_band_v1(band)?;
    let receipt = CanonicalEncoder::<DiscrepancyBandIdV1, _>::new(limits, cancellation)?
        .u64(
            Field::new(0, "mean-rel-ieee754-bits"),
            band.mean_rel.to_bits(),
        )?
        .u64(
            Field::new(1, "max-rel-ieee754-bits"),
            band.max_rel.to_bits(),
        )?
        .finish()?;
    Ok(IdentifiedDiscrepancyBandV1 { band, receipt })
}

fn preflight_model_bracket_members_v1<C>(
    bracket: &ModelBracket,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<u64, ModelBracketIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let members = bracket.identity_members();
    let member_count = bounded_len(members.len())?;
    if member_count < MIN_MODEL_BRACKET_MEMBERS_V1 {
        return Err(ModelBracketIdentityError::TooFewMembers {
            count: member_count,
            minimum: MIN_MODEL_BRACKET_MEMBERS_V1,
        });
    }
    let member_limit = limits
        .max_collection_items()
        .min(MAX_MODEL_BRACKET_MEMBERS_V1);
    if member_count > member_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: member_count,
            limit: member_limit,
        }
        .into());
    }

    let field_limit = limits
        .max_field_bytes()
        .min(MAX_MODEL_BRACKET_MEMBERS_FIELD_BYTES_V1);
    let mut field_payload_bytes = u64::from(u64::BITS / 8);
    if field_payload_bytes > field_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: field_payload_bytes,
            limit: field_limit,
        }
        .into());
    }
    let mut previous_name: Option<&str> = None;
    for (member_index, (name, qoi)) in members.iter().enumerate() {
        poll_identity_cancellation(cancellation)?;
        let member_index = bounded_len(member_index)?;
        if let Some(reason) = color_leaf_identity_reason(name) {
            return Err(ModelBracketIdentityError::InvalidMemberName {
                member_index,
                reason,
            });
        }
        if previous_name.is_some_and(|previous| previous >= name.as_str()) {
            return Err(ModelBracketIdentityError::NonCanonicalMemberOrder { member_index });
        }
        if !qoi.is_finite() {
            return Err(ModelBracketIdentityError::NonFiniteMemberQoi {
                member_index,
                bits: qoi.to_bits(),
            });
        }
        let row_bytes = 16_u64
            .checked_add(bounded_len(name.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
        let framed_row_bytes = u64::from(u64::BITS / 8)
            .checked_add(row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        field_payload_bytes = field_payload_bytes
            .checked_add(framed_row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        if field_payload_bytes > field_limit {
            return Err(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: field_payload_bytes,
                limit: field_limit,
            }
            .into());
        }
        previous_name = Some(name);
    }
    let required_stream_chunks = member_count
        .checked_mul(MODEL_BRACKET_STREAM_CHUNKS_PER_MEMBER_V1)
        .ok_or(CanonicalError::LengthOverflow)?;
    if required_stream_chunks > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::StreamChunks,
            requested: required_stream_chunks,
            limit: limits.max_collection_items(),
        }
        .into());
    }
    Ok(member_count)
}

/// Identify the exact admitted member mapping of one model bracket.
///
/// [`ModelBracket`] canonicalizes construction order by exact model name. This
/// helper consumes and retains that bracket and binds every canonical member
/// name and finite QoI bit without deriving bracket evidence or accepting a
/// legacy provenance token. Accepted signed zero remains bit-distinct.
///
/// # Errors
/// Refuses fewer than two members, an oversized member mapping, invalid names,
/// non-canonical order, non-finite QoIs, invalid limits, resource overflow, or
/// cancellation. No partial identity is published.
pub fn identify_model_bracket_v1<C>(
    bracket: ModelBracket,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedModelBracketV1, ModelBracketIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    let member_count = preflight_model_bracket_members_v1(&bracket, limits, &mut cancellation)?;
    let receipt = {
        let row_lengths = bracket.identity_members().iter().map(|(name, _)| {
            16_u64
                .checked_add(bounded_len(name.len())?)
                .ok_or(CanonicalError::LengthOverflow)
        });
        let mut rows = bracket.identity_members().iter();
        CanonicalEncoder::<ModelBracketIdV1, _>::new(limits, cancellation)?
            .ordered_bytes_stream(
                Field::new(0, "members"),
                member_count,
                row_lengths,
                |row_index, mut sink| -> Result<(), CanonicalError> {
                    let Some((name, qoi)) = rows.next() else {
                        return Err(CanonicalError::DeclaredLengthMismatch {
                            declared: member_count,
                            observed: row_index,
                        });
                    };
                    sink.write(&bounded_len(name.len())?.to_le_bytes())?;
                    sink.write(name.as_bytes())?;
                    sink.write(&qoi.to_bits().to_le_bytes())?;
                    Ok(())
                },
            )
            .map_err(|error| match error {
                OrderedBytesStreamError::Canonical { source, .. }
                | OrderedBytesStreamError::Producer { source, .. } => {
                    ModelBracketIdentityError::Canonical(source)
                }
            })?
            .finish()?
    };
    Ok(IdentifiedModelBracketV1 { bracket, receipt })
}

fn preflight_bounded_field_bytes_v1(
    length: usize,
    limits: EvidenceIdentityLimits,
    hard_limit: u64,
) -> Result<u64, CanonicalError> {
    let requested = bounded_len(length)?;
    let limit = limits.max_field_bytes().min(hard_limit);
    if requested > limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested,
            limit,
        });
    }
    Ok(requested)
}

fn preflight_canonical_string_set_v1<C>(
    values: &[String],
    limits: EvidenceIdentityLimits,
    hard_limit: u64,
    cancellation: &mut C,
) -> Result<u64, CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let count = bounded_len(values.len())?;
    if count > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: count,
            limit: limits.max_collection_items(),
        });
    }
    let field_limit = limits.max_field_bytes().min(hard_limit);
    let mut field_payload_bytes = u64::from(u64::BITS / 8);
    if field_payload_bytes > field_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: field_payload_bytes,
            limit: field_limit,
        });
    }
    for value in values {
        poll_identity_cancellation(cancellation)?;
        let framed_bytes = u64::from(u64::BITS / 8)
            .checked_add(bounded_len(value.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
        field_payload_bytes = field_payload_bytes
            .checked_add(framed_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        if field_payload_bytes > field_limit {
            return Err(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: field_payload_bytes,
                limit: field_limit,
            });
        }
    }
    Ok(count)
}

/// Identify one model-form evidence slice without promoting its declarations
/// into scientific authority.
///
/// The helper consumes and retains the public/mutable [`ModelEvidence`] while
/// binding exact canonical model-card name and assumption sets, a typed
/// normalized validity child, raw discrepancy bits, and the exact `in_domain`
/// claim bit. Card names are identifiers only; this frame does not bind
/// [`ModelCard`] contents, an evaluation point, units, or occupancy evidence.
/// Empty card/assumption sets are preserved as explicit state rather than
/// interpreted as model absence by this helper.
///
/// Cards and assumptions must already be strictly byte-sorted and
/// duplicate-free. Discrepancy refuses NaN and negative values, accepts
/// positive infinity as explicit unbounded state, and preserves signed-zero
/// bits.
///
/// # Errors
/// Refuses malformed discrepancy or validity, non-canonical sets, invalid
/// limits, resource overflow, or cancellation. No partial identity is
/// published.
pub fn identify_model_evidence_v1<C>(
    model_evidence: ModelEvidence,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedModelEvidenceV1, ModelEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    if model_evidence.discrepancy_rel.is_nan() {
        return Err(ModelEvidenceIdentityError::InvalidDiscrepancy {
            bits: model_evidence.discrepancy_rel.to_bits(),
            reason: "discrepancy must not be NaN",
        });
    }
    if model_evidence.discrepancy_rel < 0.0 {
        return Err(ModelEvidenceIdentityError::InvalidDiscrepancy {
            bits: model_evidence.discrepancy_rel.to_bits(),
            reason: "discrepancy must be non-negative; positive infinity is explicit unbounded state",
        });
    }

    let (receipt, validity_receipt) = {
        let validity_receipt = identify_validity_domain_receipt_v1(
            &model_evidence.validity,
            limits,
            &mut cancellation,
        )?;
        let card_count = preflight_canonical_string_set_v1(
            &model_evidence.cards,
            limits,
            MAX_MODEL_EVIDENCE_IDENTITY_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let assumption_count = preflight_canonical_string_set_v1(
            &model_evidence.assumptions,
            limits,
            MAX_MODEL_EVIDENCE_IDENTITY_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let receipt = CanonicalEncoder::<ModelEvidenceIdV1, _>::new(limits, cancellation)?
            .canonical_set(
                Field::new(0, "model-card-names"),
                card_count,
                model_evidence.cards.iter().map(|card| card.as_bytes()),
            )?
            .canonical_set(
                Field::new(1, "assumptions"),
                assumption_count,
                model_evidence
                    .assumptions
                    .iter()
                    .map(|assumption| assumption.as_bytes()),
            )?
            .child(Field::new(2, "validity"), validity_receipt.id())?
            .u64(
                Field::new(3, "discrepancy-rel-ieee754-bits"),
                model_evidence.discrepancy_rel.to_bits(),
            )?
            .flag(Field::new(4, "in-domain"), model_evidence.in_domain)?
            .finish()?;
        (receipt, validity_receipt)
    };
    Ok(IdentifiedModelEvidenceV1 {
        model_evidence,
        validity_receipt,
        receipt,
    })
}

const MODEL_CARD_LEGACY_FNV_OFFSET_BASIS_V1: u64 = 0xcbf2_9ce4_8422_2325;
const MODEL_CARD_LEGACY_FNV_PRIME_V1: u64 = 0x0000_0100_0000_01b3;

fn model_card_legacy_provenance_v1<C>(
    bytes: &[u8],
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<ProvenanceHash, ModelCardIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    let stride = usize::try_from(limits.cancellation_poll_bytes())
        .map_err(|_| CanonicalError::LengthOverflow)?;
    let mut processed_bytes = 0_u64;
    let mut hash = MODEL_CARD_LEGACY_FNV_OFFSET_BASIS_V1;
    if cancellation.is_cancelled() {
        return Err(ModelCardIdentityError::CalibrationCrosswalkCancelled { processed_bytes });
    }
    for chunk in bytes.chunks(stride) {
        if processed_bytes != 0 && cancellation.is_cancelled() {
            return Err(ModelCardIdentityError::CalibrationCrosswalkCancelled { processed_bytes });
        }
        for byte in chunk {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(MODEL_CARD_LEGACY_FNV_PRIME_V1);
        }
        processed_bytes = processed_bytes
            .checked_add(bounded_len(chunk.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
    }
    if cancellation.is_cancelled() {
        return Err(ModelCardIdentityError::CalibrationCrosswalkCancelled { processed_bytes });
    }
    Ok(ProvenanceHash(hash))
}

const fn model_card_ambition_tag_v1(ambition: Ambition) -> u32 {
    match ambition {
        Ambition::Solid => 1,
        Ambition::Frontier => 2,
        Ambition::Moonshot => 3,
    }
}

/// Identify one model-card declaration together with exact calibration bytes.
///
/// The helper consumes and retains the public/mutable card plus its optional
/// calibration bytes. A calibrated card must supply exact bytes whose
/// incrementally recomputed FNV matches the legacy correlation token; the weak
/// token is then excluded from the strong parent frame. An uncalibrated card
/// binds a false presence bit plus the calibration schema's empty-byte child as
/// an internal absence sentinel.
///
/// Assumptions and known failures must already be strictly byte-sorted and
/// duplicate-free. Names and versions bind exact UTF-8 without normalization
/// or semantic-version parsing. Discrepancy refuses NaN and negative values,
/// accepts positive infinity as explicit unbounded state, and preserves exact
/// signed-zero bits.
///
/// # Errors
/// Refuses calibration shape/crosswalk mismatch, malformed discrepancy or
/// validity, non-canonical sets, invalid limits, resource overflow, or
/// cancellation. No partial identity is published.
#[allow(
    clippy::too_many_lines,
    reason = "one linear helper keeps the legacy crosswalk and both typed child bindings auditable"
)]
pub fn identify_model_card_v1<C>(
    card: ModelCard,
    calibration_bytes: Option<Vec<u8>>,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedModelCardV1, ModelCardIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    let (receipt, validity_receipt, calibration_source_receipt) = {
        let declared_calibration = card.calibration;
        let supplied_calibration = calibration_bytes.is_some();
        if declared_calibration.is_some() != supplied_calibration {
            return Err(ModelCardIdentityError::CalibrationPresenceMismatch {
                declared: declared_calibration.is_some(),
                supplied: supplied_calibration,
            });
        }

        preflight_bounded_field_bytes_v1(
            card.name.len(),
            limits,
            MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1,
        )?;
        preflight_bounded_field_bytes_v1(
            card.version.len(),
            limits,
            MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1,
        )?;
        let assumption_count = preflight_canonical_string_set_v1(
            &card.assumptions,
            limits,
            MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let known_failure_count = preflight_canonical_string_set_v1(
            &card.known_failures,
            limits,
            MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let calibration_slice = calibration_bytes.as_deref().unwrap_or(&[]);
        preflight_bounded_field_bytes_v1(
            calibration_slice.len(),
            limits,
            MAX_MODEL_CARD_CALIBRATION_SOURCE_BYTES_V1,
        )?;
        if card.discrepancy_rel.is_nan() {
            return Err(ModelCardIdentityError::InvalidDiscrepancy {
                bits: card.discrepancy_rel.to_bits(),
                reason: "discrepancy must not be NaN",
            });
        }
        if card.discrepancy_rel < 0.0 {
            return Err(ModelCardIdentityError::InvalidDiscrepancy {
                bits: card.discrepancy_rel.to_bits(),
                reason: "discrepancy must be non-negative; positive infinity is explicit unbounded state",
            });
        }
        if let Some(declared) = declared_calibration {
            let computed =
                model_card_legacy_provenance_v1(calibration_slice, limits, &mut cancellation)?;
            if computed != declared {
                return Err(ModelCardIdentityError::CalibrationCorrelationMismatch {
                    declared: declared.0,
                    computed: computed.0,
                });
            }
        }

        let validity_receipt =
            identify_validity_domain_receipt_v1(&card.validity, limits, &mut cancellation)?;
        let calibration_source_receipt =
            CanonicalEncoder::<ModelCardCalibrationSourceIdV1, _>::new(limits, || {
                cancellation.is_cancelled()
            })?
            .bytes(
                Field::new(0, "canonical-calibration-artifact"),
                calibration_slice,
            )?
            .finish()?;

        let receipt = CanonicalEncoder::<ModelCardIdV1, _>::new(limits, cancellation)?
            .utf8(Field::new(0, "name"), &card.name)?
            .utf8(Field::new(1, "version"), &card.version)?
            .variant(
                Field::new(2, "ambition"),
                model_card_ambition_tag_v1(card.ambition),
                &[],
            )?
            .canonical_set(
                Field::new(3, "assumptions"),
                assumption_count,
                card.assumptions.iter().map(|value| value.as_bytes()),
            )?
            .child(Field::new(4, "validity"), validity_receipt.id())?
            .canonical_set(
                Field::new(5, "known-failures"),
                known_failure_count,
                card.known_failures.iter().map(|value| value.as_bytes()),
            )?
            .flag(Field::new(6, "calibration-present"), supplied_calibration)?
            .child(
                Field::new(7, "calibration-source"),
                calibration_source_receipt.id(),
            )?
            .u64(
                Field::new(8, "discrepancy-rel-ieee754-bits"),
                card.discrepancy_rel.to_bits(),
            )?
            .finish()?;
        (receipt, validity_receipt, calibration_source_receipt)
    };
    Ok(IdentifiedModelCardV1 {
        card,
        calibration_bytes,
        validity_receipt,
        calibration_source_receipt,
        receipt,
    })
}

fn card_bound_model_card_field_bytes_v1(card_count: u64) -> Result<u64, CanonicalError> {
    let descriptor_bytes = [
        1_u64,
        u64::from(u64::BITS / 8),
        bounded_len(ModelCardIdentitySchemaV1::DOMAIN.len())?,
        u64::from(u64::BITS / 8),
        bounded_len(ModelCardIdentitySchemaV1::NAME.len())?,
        32,
        u64::from(u32::BITS / 8),
        u64::from(u64::BITS / 8),
        bounded_len(ModelCardIdentitySchemaV1::CONTEXT.len())?,
    ]
    .into_iter()
    .try_fold(0_u64, |sum, value| {
        sum.checked_add(value).ok_or(CanonicalError::LengthOverflow)
    })?;
    u64::from(u64::BITS / 8)
        .checked_add(descriptor_bytes)
        .and_then(|bytes| {
            card_count
                .checked_mul(32)
                .and_then(|roots| bytes.checked_add(roots))
        })
        .ok_or(CanonicalError::LengthOverflow)
}

fn preflight_card_bound_model_card_field_v1(
    card_count: u64,
    limits: EvidenceIdentityLimits,
) -> Result<u64, CanonicalError> {
    let requested_field_bytes = card_bound_model_card_field_bytes_v1(card_count)?;
    let field_limit = limits
        .max_field_bytes()
        .min(MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1);
    if requested_field_bytes > field_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: requested_field_bytes,
            limit: field_limit,
        });
    }
    Ok(requested_field_bytes)
}

fn compare_card_bound_model_name_v1<C>(
    card_index: u64,
    declared: &str,
    supplied: &str,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<(), CardBoundModelEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let stride = usize::try_from(limits.cancellation_poll_bytes())
        .map_err(|_| CanonicalError::LengthOverflow)?;
    let declared = declared.as_bytes();
    let supplied = supplied.as_bytes();
    let common_bytes = declared.len().min(supplied.len());
    let mut compared_bytes = 0_usize;
    if cancellation.is_cancelled() {
        return Err(
            CardBoundModelEvidenceIdentityError::CardNameComparisonCancelled {
                card_index,
                compared_bytes: 0,
            },
        );
    }
    while compared_bytes < common_bytes {
        if compared_bytes != 0 && cancellation.is_cancelled() {
            return Err(
                CardBoundModelEvidenceIdentityError::CardNameComparisonCancelled {
                    card_index,
                    compared_bytes: bounded_len(compared_bytes)?,
                },
            );
        }
        let end = compared_bytes.saturating_add(stride).min(common_bytes);
        let declared_chunk = &declared[compared_bytes..end];
        let supplied_chunk = &supplied[compared_bytes..end];
        if let Some(different_in_chunk) = declared_chunk
            .iter()
            .zip(supplied_chunk)
            .position(|(left, right)| left != right)
        {
            return Err(CardBoundModelEvidenceIdentityError::CardNameMismatch {
                card_index,
                declared_bytes: bounded_len(declared.len())?,
                supplied_bytes: bounded_len(supplied.len())?,
                first_different_byte: bounded_len(
                    compared_bytes
                        .checked_add(different_in_chunk)
                        .ok_or(CanonicalError::LengthOverflow)?,
                )?,
            });
        }
        compared_bytes = end;
    }
    if cancellation.is_cancelled() {
        return Err(
            CardBoundModelEvidenceIdentityError::CardNameComparisonCancelled {
                card_index,
                compared_bytes: bounded_len(compared_bytes)?,
            },
        );
    }
    if declared.len() != supplied.len() {
        return Err(CardBoundModelEvidenceIdentityError::CardNameMismatch {
            card_index,
            declared_bytes: bounded_len(declared.len())?,
            supplied_bytes: bounded_len(supplied.len())?,
            first_different_byte: bounded_len(common_bytes)?,
        });
    }
    Ok(())
}

/// Bind one exact model-evidence declaration to the ordered opaque model-card
/// declarations named by its canonical card-name set.
///
/// The helper consumes and retains both the opaque model-evidence child and all
/// opaque model-card children. Counts must agree, and every card name must
/// match the corresponding model-evidence name byte-for-byte in the caller's
/// existing canonical order. No sorting, normalization, aliasing, or implicit
/// lookup occurs.
///
/// This join proves only positional name-to-declaration attachment. It does not
/// prove that the model evidence was derived from or agrees with the cards,
/// that any model ran, or that assumptions, validity, discrepancy, occupancy,
/// calibration, registry state, or scientific claims are correct.
///
/// # Errors
/// Refuses count/name mismatch, invalid limits, exact ordered-child field or
/// collection overflow, or cancellation during comparison/framing. No partial
/// join is published.
pub fn identify_card_bound_model_evidence_v1<C>(
    model_evidence: IdentifiedModelEvidenceV1,
    model_cards: Vec<IdentifiedModelCardV1>,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedCardBoundModelEvidenceV1, CardBoundModelEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;

    let declared_names = &model_evidence.model_evidence().cards;
    let declared_count = bounded_len(declared_names.len())?;
    let supplied_count = bounded_len(model_cards.len())?;
    if declared_count != supplied_count {
        return Err(CardBoundModelEvidenceIdentityError::CardCountMismatch {
            declared: declared_count,
            supplied: supplied_count,
        });
    }
    if supplied_count > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: supplied_count,
            limit: limits.max_collection_items(),
        }
        .into());
    }
    preflight_card_bound_model_card_field_v1(supplied_count, limits)?;

    for (card_index, (declared_name, supplied_card)) in
        declared_names.iter().zip(&model_cards).enumerate()
    {
        compare_card_bound_model_name_v1(
            bounded_len(card_index)?,
            declared_name,
            &supplied_card.card().name,
            limits,
            &mut cancellation,
        )?;
    }
    poll_identity_cancellation(&mut cancellation)?;

    let receipt = CanonicalEncoder::<CardBoundModelEvidenceIdV1, _>::new(limits, cancellation)?
        .child(Field::new(0, "model-evidence"), model_evidence.id())?
        .ordered_children(
            Field::new(1, "model-card-declarations"),
            supplied_count,
            model_cards.iter().map(IdentifiedModelCardV1::id),
        )?
        .finish()?;
    Ok(IdentifiedCardBoundModelEvidenceV1 {
        model_evidence,
        model_cards,
        receipt,
    })
}

fn preflight_certified_f64_sensitivity_v1<C>(
    certified: &Certified<f64>,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<u64, CanonicalError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let sensitivity = &certified.evidence().sensitivity.d_qoi;
    let count = bounded_len(sensitivity.len())?;
    if count > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: count,
            limit: limits.max_collection_items(),
        });
    }
    let field_limit = limits
        .max_field_bytes()
        .min(MAX_CERTIFIED_F64_EVIDENCE_FIELD_BYTES_V1);
    let mut field_payload_bytes = u64::from(u64::BITS / 8);
    if field_payload_bytes > field_limit {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::FieldBytes,
            requested: field_payload_bytes,
            limit: field_limit,
        });
    }
    for parameter in sensitivity.keys() {
        poll_identity_cancellation(cancellation)?;
        let row_bytes = 16_u64
            .checked_add(bounded_len(parameter.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
        let framed_row_bytes = u64::from(u64::BITS / 8)
            .checked_add(row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        field_payload_bytes = field_payload_bytes
            .checked_add(framed_row_bytes)
            .ok_or(CanonicalError::LengthOverflow)?;
        if field_payload_bytes > field_limit {
            return Err(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: field_payload_bytes,
                limit: field_limit,
            });
        }
    }
    let required_stream_chunks = count
        .checked_mul(CERTIFIED_F64_SENSITIVITY_STREAM_CHUNKS_PER_ROW_V1)
        .ok_or(CanonicalError::LengthOverflow)?;
    if required_stream_chunks > limits.max_collection_items() {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::StreamChunks,
            requested: required_stream_chunks,
            limit: limits.max_collection_items(),
        });
    }
    Ok(count)
}

/// Identify the strong semantic projection carried by one locally certified
/// scalar value.
///
/// The projection binds the carried scalar, QoI, numerical and statistical
/// certificates, canonical model-card and assumption sets, a typed validity
/// child, discrepancy and in-domain state, every sensitivity row, and exact
/// adjoint-hook presence. The existing `Certified<f64>` is consumed and
/// retained without changing its layout or trust state.
///
/// Cards and assumptions must already satisfy their documented sorted,
/// duplicate-free set representation. Legacy FNV provenance and adjoint values
/// are deliberately excluded from the strong root rather than rehashed into
/// apparent authority; `None` versus `Some` adjoint remains semantic claim
/// state, while the original correlation tokens remain inspectable through the
/// attached certified record.
///
/// # Errors
/// Refuses a non-canonical set, validity-child refusal, invalid limits,
/// resource overflow, or cancellation. No partial identity is published.
#[allow(
    clippy::too_many_lines,
    reason = "one linear canonical frame keeps all field-order and ownership invariants visible"
)]
pub fn identify_certified_f64_evidence_v1<C>(
    certified: Certified<f64>,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedCertifiedF64EvidenceV1, CertifiedF64EvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let (receipt, validity_receipt) = {
        let evidence = certified.evidence();
        let validity_receipt = identify_validity_domain_receipt_v1(
            &evidence.model.validity,
            limits,
            &mut cancellation,
        )?;
        let card_count = preflight_canonical_string_set_v1(
            &evidence.model.cards,
            limits,
            MAX_CERTIFIED_F64_EVIDENCE_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let assumption_count = preflight_canonical_string_set_v1(
            &evidence.model.assumptions,
            limits,
            MAX_CERTIFIED_F64_EVIDENCE_FIELD_BYTES_V1,
            &mut cancellation,
        )?;
        let sensitivity_count =
            preflight_certified_f64_sensitivity_v1(&certified, limits, &mut cancellation)?;

        let mut statistical_payload = [0_u8; 16];
        let (statistical_tag, statistical_payload_len) = match evidence.statistical {
            StatisticalCertificate::None => (1, 0),
            StatisticalCertificate::EValue { e, alpha } => {
                statistical_payload[..8].copy_from_slice(&e.to_bits().to_le_bytes());
                statistical_payload[8..].copy_from_slice(&alpha.to_bits().to_le_bytes());
                (2, statistical_payload.len())
            }
            StatisticalCertificate::HalfWidth {
                half_width,
                confidence,
            } => {
                statistical_payload[..8].copy_from_slice(&half_width.to_bits().to_le_bytes());
                statistical_payload[8..].copy_from_slice(&confidence.to_bits().to_le_bytes());
                (3, statistical_payload.len())
            }
        };
        let encoder = CanonicalEncoder::<CertifiedF64EvidenceIdV1, _>::new(limits, cancellation)?
            .finite_f64(Field::new(0, "value"), evidence.value)?
            .finite_f64(Field::new(1, "qoi"), evidence.qoi)?
            .variant(
                Field::new(2, "numerical-kind"),
                numerical_certificate_kind_tag_v1(evidence.numerical.kind),
                &[],
            )?
            .finite_f64(Field::new(3, "numerical-lo"), evidence.numerical.lo)?
            .finite_f64(Field::new(4, "numerical-hi"), evidence.numerical.hi)?
            .variant(
                Field::new(5, "statistical"),
                statistical_tag,
                &statistical_payload[..statistical_payload_len],
            )?
            .canonical_set(
                Field::new(6, "model-cards"),
                card_count,
                evidence.model.cards.iter().map(|card| card.as_bytes()),
            )?
            .canonical_set(
                Field::new(7, "model-assumptions"),
                assumption_count,
                evidence
                    .model
                    .assumptions
                    .iter()
                    .map(|assumption| assumption.as_bytes()),
            )?
            .child(Field::new(8, "model-validity"), validity_receipt.id())?
            .u64(
                Field::new(9, "model-discrepancy-ieee754-bits"),
                evidence.model.discrepancy_rel.to_bits(),
            )?
            .flag(Field::new(10, "model-in-domain"), evidence.model.in_domain)?;

        let sensitivity = &evidence.sensitivity.d_qoi;
        let row_lengths = sensitivity.keys().map(|parameter| {
            bounded_len(parameter.len()).and_then(|parameter_bytes| {
                16_u64
                    .checked_add(parameter_bytes)
                    .ok_or(CanonicalError::LengthOverflow)
            })
        });
        let mut rows = sensitivity.iter();
        let encoder = encoder
            .ordered_bytes_stream(
                Field::new(11, "sensitivity"),
                sensitivity_count,
                row_lengths,
                |row_index, mut sink| -> Result<(), CanonicalError> {
                    let Some((parameter, derivative)) = rows.next() else {
                        return Err(CanonicalError::DeclaredLengthMismatch {
                            declared: sensitivity_count,
                            observed: row_index,
                        });
                    };
                    sink.write(&bounded_len(parameter.len())?.to_le_bytes())?;
                    sink.write(parameter.as_bytes())?;
                    sink.write(&derivative.to_bits().to_le_bytes())?;
                    Ok(())
                },
            )
            .map_err(|error| match error {
                OrderedBytesStreamError::Canonical { source, .. }
                | OrderedBytesStreamError::Producer { source, .. } => {
                    CertifiedF64EvidenceIdentityError::Canonical(source)
                }
            })?
            .flag(
                Field::new(12, "legacy-adjoint-correlation-present"),
                evidence.adjoint_ref.is_some(),
            )?
            .finish()?;
        (encoder, validity_receipt)
    };
    Ok(IdentifiedCertifiedF64EvidenceV1 {
        certified,
        validity_receipt,
        receipt,
    })
}

/// Identify one exact source artifact for a certified-f64 legacy crosswalk.
///
/// The helper consumes and retains the exact domain, nonzero source-schema
/// version, and canonical bytes. Empty byte artifacts are valid and remain
/// distinct from source absence at the parent layer. The receipt is
/// content-bound but unanchored.
///
/// # Errors
/// Refuses an empty domain, schema version zero, invalid limits, resource
/// overflow, or cancellation. No partial source identity is published.
pub fn identify_certified_f64_source_v1<C>(
    domain: String,
    schema_version: u32,
    canonical_bytes: Vec<u8>,
    limits: EvidenceIdentityLimits,
    cancellation: C,
) -> Result<CertifiedF64SourceV1, CertifiedF64SourceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if domain.is_empty() {
        return Err(CertifiedF64SourceIdentityError::EmptySourceDomain);
    }
    if schema_version == 0 {
        return Err(CertifiedF64SourceIdentityError::ZeroSourceSchemaVersion);
    }
    preflight_bounded_field_bytes_v1(
        domain.len(),
        limits,
        MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1,
    )?;
    preflight_bounded_field_bytes_v1(
        canonical_bytes.len(),
        limits,
        MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1,
    )?;
    let receipt = CanonicalEncoder::<CertifiedF64SourceIdV1, _>::new(limits, cancellation)?
        .utf8(Field::new(0, "source-domain"), &domain)?
        .u64(
            Field::new(1, "source-schema-version"),
            u64::from(schema_version),
        )?
        .bytes(Field::new(2, "canonical-source"), &canonical_bytes)?
        .finish()?;
    Ok(CertifiedF64SourceV1 {
        domain,
        schema_version,
        canonical_bytes,
        receipt,
    })
}

const CERTIFIED_F64_LEGACY_FNV_OFFSET_BASIS_V1: u64 = 0xcbf2_9ce4_8422_2325;
const CERTIFIED_F64_LEGACY_FNV_PRIME_V1: u64 = 0x0000_0100_0000_01b3;
// Exact fs-obs replay-identity v1 bytes for kind `provenance-chain`, one
// three-byte operation name, and two eight-byte `operand` child roots.
const CERTIFIED_F64_LEGACY_CHAIN_CANONICAL_BYTES_V1: usize = 118;
const _: () =
    assert!(CERTIFIED_F64_LEGACY_FNV_OFFSET_BASIS_V1 == MODEL_CARD_LEGACY_FNV_OFFSET_BASIS_V1);
const _: () = assert!(CERTIFIED_F64_LEGACY_FNV_PRIME_V1 == MODEL_CARD_LEGACY_FNV_PRIME_V1);

#[derive(Debug, Clone, PartialEq, Eq)]
enum CertifiedF64LegacyProvenanceError {
    Cancelled { processed_bytes: u64 },
    Canonical(CanonicalError),
}

impl From<CanonicalError> for CertifiedF64LegacyProvenanceError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

fn certified_f64_legacy_provenance_v1<C>(
    bytes: &[u8],
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<ProvenanceHash, CertifiedF64LegacyProvenanceError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    let stride = usize::try_from(limits.cancellation_poll_bytes())
        .map_err(|_| CanonicalError::LengthOverflow)?;
    let mut processed_bytes = 0_u64;
    let mut hash = CERTIFIED_F64_LEGACY_FNV_OFFSET_BASIS_V1;
    if cancellation.is_cancelled() {
        return Err(CertifiedF64LegacyProvenanceError::Cancelled { processed_bytes });
    }
    for chunk in bytes.chunks(stride) {
        if processed_bytes != 0 && cancellation.is_cancelled() {
            return Err(CertifiedF64LegacyProvenanceError::Cancelled { processed_bytes });
        }
        for byte in chunk {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(CERTIFIED_F64_LEGACY_FNV_PRIME_V1);
        }
        processed_bytes = processed_bytes
            .checked_add(bounded_len(chunk.len())?)
            .ok_or(CanonicalError::LengthOverflow)?;
    }
    if cancellation.is_cancelled() {
        return Err(CertifiedF64LegacyProvenanceError::Cancelled { processed_bytes });
    }
    Ok(ProvenanceHash(hash))
}

fn sourced_certified_f64_legacy_provenance_v1<C>(
    source: &'static str,
    bytes: &[u8],
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<ProvenanceHash, SourcedCertifiedF64EvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    certified_f64_legacy_provenance_v1(bytes, limits, cancellation).map_err(|error| match error {
        CertifiedF64LegacyProvenanceError::Cancelled { processed_bytes } => {
            SourcedCertifiedF64EvidenceIdentityError::SourceCrosswalkCancelled {
                source,
                processed_bytes,
            }
        }
        CertifiedF64LegacyProvenanceError::Canonical(error) => {
            SourcedCertifiedF64EvidenceIdentityError::Canonical(error)
        }
    })
}

fn certified_f64_bounded_legacy_chain_v1(
    operation: Op,
    left: ProvenanceHash,
    right: ProvenanceHash,
) -> Result<ProvenanceHash, IdentityBuildError> {
    let identity = BoundedIdentityBuilder::new(
        "provenance-chain",
        CERTIFIED_F64_LEGACY_CHAIN_CANONICAL_BYTES_V1,
    )?
    .str("op", operation.name())?
    .child_root64("operand", left.0)?
    .child_root64("operand", right.0)?
    .finish();
    debug_assert_eq!(
        identity.canonical_bytes().len(),
        CERTIFIED_F64_LEGACY_CHAIN_CANONICAL_BYTES_V1
    );
    Ok(ProvenanceHash(identity.root()))
}

/// Attach exact producer and optional adjoint sources to one opaque certified
/// scalar semantic identity.
///
/// Exact producer bytes must reproduce the certified record's legacy FNV-1a-64
/// provenance. Exact adjoint-source presence must match `adjoint_ref`, and a
/// present source must reproduce that correlation too. Both weak `u64` values
/// are admission crosswalks only and never enter the strong parent frame. The
/// source domain, source-schema version, and exact bytes remain identity-bearing
/// through typed children and are retained for inspection and replay.
///
/// This v1 helper provides only a leaf-source binding. It cannot recover how a
/// bare legacy `u64` was formed or distinguish an accidental/adversarial FNV
/// collision from the intended source. Callers must not represent
/// `ProvenanceHash::chain` lineage with this leaf API; use
/// [`identify_certified_f64_lineage_source_v1`] and
/// [`compose_certified_f64_lineage_nodes_v1`] for chained derivations.
///
/// # Errors
/// Refuses adjoint-source presence mismatch, either legacy-FNV correlation
/// mismatch, invalid limits, resource overflow, or cancellation. No partial
/// parent identity is published.
pub fn identify_sourced_certified_f64_evidence_v1<C>(
    certified_evidence: IdentifiedCertifiedF64EvidenceV1,
    producer_source: CertifiedF64SourceV1,
    adjoint_source: Option<CertifiedF64SourceV1>,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedSourcedCertifiedF64EvidenceV1, SourcedCertifiedF64EvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;

    let (declared_producer, declared_adjoint) = {
        let evidence = certified_evidence.certified().evidence();
        (evidence.provenance, evidence.adjoint_ref)
    };
    if declared_adjoint.is_some() != adjoint_source.is_some() {
        return Err(
            SourcedCertifiedF64EvidenceIdentityError::AdjointSourcePresenceMismatch {
                declared: declared_adjoint.is_some(),
                supplied: adjoint_source.is_some(),
            },
        );
    }

    let computed_producer = sourced_certified_f64_legacy_provenance_v1(
        "producer",
        producer_source.canonical_bytes(),
        limits,
        &mut cancellation,
    )?;
    if computed_producer != declared_producer {
        return Err(
            SourcedCertifiedF64EvidenceIdentityError::SourceCorrelationMismatch {
                source: "producer",
                declared: declared_producer.0,
                computed: computed_producer.0,
            },
        );
    }

    if let (Some(declared), Some(source)) = (declared_adjoint, adjoint_source.as_ref()) {
        let computed = sourced_certified_f64_legacy_provenance_v1(
            "adjoint",
            source.canonical_bytes(),
            limits,
            &mut cancellation,
        )?;
        if computed != declared {
            return Err(
                SourcedCertifiedF64EvidenceIdentityError::SourceCorrelationMismatch {
                    source: "adjoint",
                    declared: declared.0,
                    computed: computed.0,
                },
            );
        }
    }

    let producer_source_id = producer_source.id();
    let adjoint_source_id = adjoint_source.as_ref().map(CertifiedF64SourceV1::id);
    let adjoint_source_count = u64::from(adjoint_source_id.is_some());
    let receipt =
        CanonicalEncoder::<SourcedCertifiedF64EvidenceIdV1, _>::new(limits, cancellation)?
            .child(
                Field::new(0, "certified-f64-evidence"),
                certified_evidence.id(),
            )?
            .child(Field::new(1, "producer-source"), producer_source_id)?
            .flag(
                Field::new(2, "adjoint-source-present"),
                adjoint_source_id.is_some(),
            )?
            .ordered_children(
                Field::new(3, "adjoint-source"),
                adjoint_source_count,
                adjoint_source_id,
            )?
            .finish()?;

    Ok(IdentifiedSourcedCertifiedF64EvidenceV1 {
        certified_evidence,
        producer_source,
        adjoint_source,
        receipt,
    })
}

const CERTIFIED_F64_LINEAGE_PARENT_REFERENCE_BYTES_V1: usize = 65;

fn certified_f64_lineage_parent_reference_bytes_v1(
    parent: CertifiedF64LineageNodeIdV1,
) -> [u8; CERTIFIED_F64_LINEAGE_PARENT_REFERENCE_BYTES_V1] {
    let mut output = [0_u8; CERTIFIED_F64_LINEAGE_PARENT_REFERENCE_BYTES_V1];
    output[0] = CertifiedF64LineageNodeIdV1::ROLE.tag();
    output[1..33].copy_from_slice(
        SchemaId::<CertifiedF64LineageNodeIdentitySchemaV1>::for_schema().as_bytes(),
    );
    output[33..].copy_from_slice(parent.as_bytes());
    output
}

fn build_certified_f64_lineage_node_v1<C>(
    operation: CertifiedF64LineageOperationV1,
    source: Option<CertifiedF64SourceIdV1>,
    parents: Option<[CertifiedF64LineageNodeIdV1; 2]>,
    limits: EvidenceIdentityLimits,
    cancellation: C,
) -> Result<CertifiedF64LineageNodeReceiptV1, CertifiedF64LineageIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let source_count = u64::from(source.is_some());
    let parent_count = if parents.is_some() { 2_u64 } else { 0 };
    let parent_count_usize = if parents.is_some() { 2 } else { 0 };
    let parent_field_bytes = core::mem::size_of::<u64>()
        + parent_count_usize
            * (core::mem::size_of::<u64>() + CERTIFIED_F64_LINEAGE_PARENT_REFERENCE_BYTES_V1);
    preflight_bounded_field_bytes_v1(
        parent_field_bytes,
        limits,
        MAX_CERTIFIED_F64_LINEAGE_PARENT_FIELD_BYTES_V1,
    )?;
    let parent_rows =
        parents.map(|parents| parents.map(certified_f64_lineage_parent_reference_bytes_v1));
    let kind = operation.kind();
    let parent_semantics = operation.parent_semantics();

    Ok(
        CanonicalEncoder::<CertifiedF64LineageNodeIdV1, _>::new(limits, cancellation)?
            .variant(Field::new(0, "node-kind"), kind.tag(), &[])?
            .variant(Field::new(1, "operation"), operation.tag(), &[])?
            .variant(
                Field::new(2, "parent-semantics"),
                parent_semantics.tag(),
                &[],
            )?
            .u64(
                Field::new(3, "lineage-algorithm-version"),
                u64::from(CERTIFIED_F64_LINEAGE_ALGORITHM_VERSION_V1),
            )?
            .ordered_children(Field::new(4, "source"), source_count, source)?
            .ordered_bytes(
                Field::new(5, "parents"),
                parent_count,
                parent_rows
                    .iter()
                    .flat_map(|rows| rows.iter())
                    .map(|row| row.as_slice()),
            )?
            .finish()?,
    )
}

/// Root one strong lineage node at a retained exact certified-f64 source.
///
/// The exact source domain, version, and bytes remain authoritative through
/// the typed source child. A pinned, cancellable FNV-1a-64 pass also derives the
/// legacy compatibility token, but that weak value is excluded from the strong
/// frame and proves no source origin or authenticity.
///
/// # Errors
/// Refuses invalid limits, source-work or frame budget overflow, and
/// cancellation. No partial lineage node is published.
pub fn identify_certified_f64_lineage_source_v1<C>(
    source: CertifiedF64SourceV1,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<CertifiedF64LineageNodeV1, CertifiedF64LineageIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    poll_identity_cancellation(&mut cancellation)?;
    preflight_bounded_field_bytes_v1(
        source.canonical_bytes().len(),
        limits,
        MAX_CERTIFIED_F64_SOURCE_FIELD_BYTES_V1,
    )?;
    let legacy_provenance =
        certified_f64_legacy_provenance_v1(source.canonical_bytes(), limits, &mut cancellation)
            .map_err(|error| match error {
                CertifiedF64LegacyProvenanceError::Cancelled { processed_bytes } => {
                    CertifiedF64LineageIdentityError::SourceCrosswalkCancelled { processed_bytes }
                }
                CertifiedF64LegacyProvenanceError::Canonical(error) => {
                    CertifiedF64LineageIdentityError::Canonical(error)
                }
            })?;
    let receipt = build_certified_f64_lineage_node_v1(
        CertifiedF64LineageOperationV1::Source,
        Some(source.id()),
        None,
        limits,
        cancellation,
    )?;
    Ok(CertifiedF64LineageNodeV1 {
        source: Some(source),
        parents: None,
        legacy_provenance,
        receipt,
        operation: CertifiedF64LineageOperationV1::Source,
    })
}

/// Identify one ordered binary provenance derivation over typed parent nodes.
///
/// Parent order, duplicates, and association-tree shape are all semantic for
/// every operation, including Add/Mul/Min/Max. This exactly preserves the
/// current order-sensitive legacy chain and avoids silently claiming float or
/// evidence-metadata commutativity. The helper derives a weak compatibility
/// token through a fixed 118-byte, allocation-fallible builder that is
/// schema-equivalent to `ProvenanceHash::chain`, but only the exact operation
/// and typed ordered parents enter the strong root.
///
/// This graph does not bind an output `Evidence`, scalar value, certificate,
/// units, quantity kind, or execution. Those require a separate typed semantic
/// attachment rather than inference from lineage alone.
///
/// # Errors
/// Refuses conflicting parent observations, bounded compatibility-chain
/// allocation/framing failure, invalid limits, resource overflow, or
/// cancellation. No partial lineage node is published.
pub fn compose_certified_f64_lineage_nodes_v1<C>(
    operation: CertifiedF64LineageCompositionOpV1,
    left: &CertifiedF64LineageNodeV1,
    right: &CertifiedF64LineageNodeV1,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<CertifiedF64LineageNodeV1, CertifiedF64LineageIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    poll_identity_cancellation(&mut cancellation)?;
    if matches!(
        adjudicate(
            ObservedIdentity::from_receipt(left.receipt()),
            ObservedIdentity::from_receipt(right.receipt()),
        ),
        IdentityAdjudication::Refused(_)
    ) {
        return Err(CertifiedF64LineageIdentityError::ParentObservationConflict);
    }

    let node_operation = operation.node_operation();
    let parents = [left.id(), right.id()];
    poll_identity_cancellation(&mut cancellation)?;
    let legacy_operation = operation.legacy_operation();
    let legacy_provenance = certified_f64_bounded_legacy_chain_v1(
        legacy_operation,
        left.legacy_provenance(),
        right.legacy_provenance(),
    )?;
    poll_identity_cancellation(&mut cancellation)?;
    let receipt = build_certified_f64_lineage_node_v1(
        node_operation,
        None,
        Some(parents),
        limits,
        cancellation,
    )?;
    Ok(CertifiedF64LineageNodeV1 {
        source: None,
        parents: Some(parents),
        legacy_provenance,
        receipt,
        operation: node_operation,
    })
}

const fn decision_assessment_uncertainty_source_tag_v1(source: UncertaintySource) -> u32 {
    match source {
        UncertaintySource::ModelForm => 1,
        UncertaintySource::Statistical => 2,
        UncertaintySource::Numerical => 3,
    }
}

const fn decision_assessment_advice_tag_v1(advice: EscalationAdvice) -> u32 {
    match advice {
        EscalationAdvice::NoneNeeded => 1,
        EscalationAdvice::RefineNumerics => 2,
        EscalationAdvice::GatherMoreSamples => 3,
        EscalationAdvice::EscalateModelFidelity => 4,
    }
}

fn validate_decision_assessment_band_v1(
    field: &'static str,
    value: f64,
) -> Result<(), CertifiedF64DecisionAssessmentIdentityError> {
    if value.is_nan() {
        return Err(
            CertifiedF64DecisionAssessmentIdentityError::InvalidDerivedBand {
                field,
                bits: value.to_bits(),
                reason: "derived relative band must not be NaN",
            },
        );
    }
    if value < 0.0 {
        return Err(
            CertifiedF64DecisionAssessmentIdentityError::InvalidDerivedBand {
                field,
                bits: value.to_bits(),
                reason: "derived relative band must be non-negative",
            },
        );
    }
    Ok(())
}

/// Recompute and identify one local decision assessment over an opaque
/// certified-f64 semantic child.
///
/// The frame binds the assessment-algorithm version, the complete typed child,
/// exact threshold bits, all recomputed relative-band bits, the local status
/// and deterministic dominant-source tag, and the resulting advice. The
/// presentation-only status detail string is excluded. Positive infinity is an
/// honest non-decision-grade band; accepted signed zero remains bit-distinct.
///
/// # Errors
/// Refuses a non-finite or negative threshold, invalid derived-band state,
/// invalid limits, resource overflow, or cancellation. No partial assessment
/// identity is published.
#[allow(
    clippy::too_many_lines,
    reason = "one linear frame keeps recomputation, status payload, and field order auditable"
)]
pub fn identify_certified_f64_decision_assessment_v1<C>(
    certified_evidence: IdentifiedCertifiedF64EvidenceV1,
    threshold_rel: f64,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<IdentifiedCertifiedF64DecisionAssessmentV1, CertifiedF64DecisionAssessmentIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive").into(),
        );
    }
    poll_identity_cancellation(&mut cancellation)?;
    if !threshold_rel.is_finite() {
        return Err(
            CertifiedF64DecisionAssessmentIdentityError::InvalidThreshold {
                bits: threshold_rel.to_bits(),
                reason: "threshold must be finite",
            },
        );
    }
    if threshold_rel < 0.0 {
        return Err(
            CertifiedF64DecisionAssessmentIdentityError::InvalidThreshold {
                bits: threshold_rel.to_bits(),
                reason: "threshold must be non-negative",
            },
        );
    }

    let certified = certified_evidence.certified();
    let breakdown = certified.breakdown();
    let total_rel = breakdown.total_rel();
    for (field, value) in [
        ("numerical_rel", breakdown.numerical_rel),
        ("statistical_rel", breakdown.statistical_rel),
        ("model_rel", breakdown.model_rel),
        ("total_rel", total_rel),
    ] {
        validate_decision_assessment_band_v1(field, value)?;
    }
    let status = certified.assess(threshold_rel);
    let advice = certified.escalation_advice(threshold_rel);
    let mut status_payload = [0_u8; 4];
    let (status_tag, status_payload_len) = match &status {
        DecisionStatus::DecisionGrade => (1, 0),
        DecisionStatus::NotDecisionGrade { dominant, .. } => {
            status_payload.copy_from_slice(
                &decision_assessment_uncertainty_source_tag_v1(*dominant).to_le_bytes(),
            );
            (2, status_payload.len())
        }
    };

    let receipt =
        CanonicalEncoder::<CertifiedF64DecisionAssessmentIdV1, _>::new(limits, cancellation)?
            .child(
                Field::new(0, "certified-f64-evidence"),
                certified_evidence.id(),
            )?
            .u64(
                Field::new(1, "assessment-algorithm-version"),
                u64::from(DECISION_ASSESSMENT_ALGORITHM_VERSION_V1),
            )?
            .finite_f64(Field::new(2, "threshold-rel"), threshold_rel)?
            .u64(
                Field::new(3, "numerical-rel-ieee754-bits"),
                breakdown.numerical_rel.to_bits(),
            )?
            .u64(
                Field::new(4, "statistical-rel-ieee754-bits"),
                breakdown.statistical_rel.to_bits(),
            )?
            .u64(
                Field::new(5, "model-rel-ieee754-bits"),
                breakdown.model_rel.to_bits(),
            )?
            .u64(Field::new(6, "total-rel-ieee754-bits"), total_rel.to_bits())?
            .variant(
                Field::new(7, "status"),
                status_tag,
                &status_payload[..status_payload_len],
            )?
            .variant(
                Field::new(8, "advice"),
                decision_assessment_advice_tag_v1(advice),
                &[],
            )?
            .finish()?;
    Ok(IdentifiedCertifiedF64DecisionAssessmentV1 {
        certified_evidence,
        threshold_rel,
        breakdown,
        total_rel,
        status,
        advice,
        receipt,
    })
}

/// Reproduce `Color::canonical_bytes` under a hard allocation ceiling, a
/// fallible exact reservation, and byte-stride cancellation checks.
fn bounded_color_bytes<C>(
    color: &Color,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
) -> Result<Vec<u8>, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    bounded_color_bytes_with_reservation(color, limits, cancellation, |output, capacity| {
        output
            .try_reserve_exact(capacity)
            .map_err(|_| ColorBufferReservationError)
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ColorBufferReservationError;

fn bounded_color_bytes_with_reservation<C, R>(
    color: &Color,
    limits: EvidenceIdentityLimits,
    cancellation: &mut C,
    reserve: R,
) -> Result<Vec<u8>, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
    R: FnOnce(&mut Vec<u8>, usize) -> Result<(), ColorBufferReservationError>,
{
    if limits.cancellation_poll_bytes() == 0 {
        return Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::InvalidLimits("cancellation_poll_bytes must be positive"),
        ));
    }
    let cancellation_poll_bytes = usize::try_from(limits.cancellation_poll_bytes())
        .map_err(|_| ColorEvidenceIdentityError::Canonical(CanonicalError::LengthOverflow))?;
    poll_identity_cancellation(cancellation)?;
    let limit = limits
        .max_field_bytes()
        .min(limits.max_canonical_bytes())
        .min(MAX_COLOR_EVIDENCE_NODE_BYTES_V1);
    let mut length = 0;
    add_bounded_color_bytes(&mut length, 2, limit)?;
    match color {
        Color::Verified { .. } => {
            add_bounded_color_bytes(&mut length, 32, limit)?;
        }
        Color::Validated { regime, dataset } => {
            let axis_count = bounded_len(regime.bounds().len())?;
            if axis_count > limits.max_collection_items() {
                return Err(ColorEvidenceIdentityError::Canonical(
                    CanonicalError::LimitExceeded {
                        kind: LimitKind::CollectionItems,
                        requested: axis_count,
                        limit: limits.max_collection_items(),
                    },
                ));
            }
            add_bounded_color_bytes(
                &mut length,
                8_u64.checked_add(bounded_len(dataset.len())?).ok_or(
                    ColorEvidenceIdentityError::Canonical(CanonicalError::LengthOverflow),
                )?,
                limit,
            )?;
            add_bounded_color_bytes(&mut length, 8, limit)?;
            for axis in regime.bounds().keys() {
                poll_identity_cancellation(cancellation)?;
                let row_bytes = 40_u64.checked_add(bounded_len(axis.len())?).ok_or(
                    ColorEvidenceIdentityError::Canonical(CanonicalError::LengthOverflow),
                )?;
                add_bounded_color_bytes(&mut length, row_bytes, limit)?;
            }
        }
        Color::Estimated { estimator, .. } => {
            let payload = 24_u64.checked_add(bounded_len(estimator.len())?).ok_or(
                ColorEvidenceIdentityError::Canonical(CanonicalError::LengthOverflow),
            )?;
            add_bounded_color_bytes(&mut length, payload, limit)?;
        }
    }
    validate_color_payload(color).map_err(ColorEvidenceIdentityError::MalformedColor)?;
    poll_identity_cancellation(cancellation)?;

    let capacity = usize::try_from(length)
        .map_err(|_| ColorEvidenceIdentityError::Canonical(CanonicalError::LengthOverflow))?;
    let mut output = Vec::new();
    reserve(&mut output, capacity).map_err(|ColorBufferReservationError| {
        ColorEvidenceIdentityError::ColorBufferAllocationFailed {
            requested_bytes: length,
        }
    })?;
    match color {
        Color::Verified { lo, hi } => {
            append_color_bytes(
                &mut output,
                &[COLOR_ALGEBRA_VERSION as u8, 0],
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_field(
                &mut output,
                &lo.to_bits().to_le_bytes(),
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_field(
                &mut output,
                &hi.to_bits().to_le_bytes(),
                cancellation_poll_bytes,
                cancellation,
            )?;
        }
        Color::Validated { regime, dataset } => {
            append_color_bytes(
                &mut output,
                &[COLOR_ALGEBRA_VERSION as u8, 1],
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_field(
                &mut output,
                dataset.as_bytes(),
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_len(
                &mut output,
                regime.bounds().len(),
                cancellation_poll_bytes,
                cancellation,
            )?;
            for (axis, (lo, hi)) in regime.bounds() {
                push_color_field(
                    &mut output,
                    axis.as_bytes(),
                    cancellation_poll_bytes,
                    cancellation,
                )?;
                push_color_field(
                    &mut output,
                    &lo.to_bits().to_le_bytes(),
                    cancellation_poll_bytes,
                    cancellation,
                )?;
                push_color_field(
                    &mut output,
                    &hi.to_bits().to_le_bytes(),
                    cancellation_poll_bytes,
                    cancellation,
                )?;
            }
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            append_color_bytes(
                &mut output,
                &[COLOR_ALGEBRA_VERSION as u8, 2],
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_field(
                &mut output,
                estimator.as_bytes(),
                cancellation_poll_bytes,
                cancellation,
            )?;
            push_color_field(
                &mut output,
                &dispersion.to_bits().to_le_bytes(),
                cancellation_poll_bytes,
                cancellation,
            )?;
        }
    }
    debug_assert_eq!(output.len(), capacity);
    Ok(output)
}

fn parent_reference_bytes(parent: ColorEvidenceNodeIdV1) -> [u8; 65] {
    let mut output = [0_u8; 65];
    output[0] = ColorEvidenceNodeIdV1::ROLE.tag();
    output[1..33]
        .copy_from_slice(SchemaId::<ColorEvidenceNodeIdentitySchemaV1>::for_schema().as_bytes());
    output[33..].copy_from_slice(parent.as_bytes());
    output
}

fn build_color_evidence_node_v1<C>(
    operation: ColorEvidenceOperationV1,
    source: Option<ColorEvidenceSourceIdV1>,
    output: &Color,
    parents: Option<[ColorEvidenceNodeIdV1; 2]>,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<ColorEvidenceNodeReceiptV1, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let output_bytes = bounded_color_bytes(output, limits, &mut cancellation)?;
    let parent_count = if parents.is_some() { 2_u64 } else { 0 };
    let parent_rows = parents.map(|parents| parents.map(parent_reference_bytes));
    let source_count: u64 = if source.is_some() { 1 } else { 0 };
    let kind = operation.kind();
    let parent_semantics = operation.parent_semantics();

    Ok(
        CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(limits, cancellation)?
            .variant(Field::new(0, "node-kind"), kind.tag(), &[])?
            .variant(Field::new(1, "operation"), operation.tag(), &[])?
            .variant(
                Field::new(2, "parent-semantics"),
                parent_semantics.tag(),
                &[],
            )?
            .u64(
                Field::new(3, "color-algebra-version"),
                u64::from(COLOR_ALGEBRA_VERSION),
            )?
            .ordered_children(Field::new(4, "source"), source_count, source)?
            .bytes(Field::new(5, "output-color"), &output_bytes)?
            .ordered_bytes(
                Field::new(6, "parents"),
                parent_count,
                parent_rows
                    .iter()
                    .flat_map(|rows| rows.iter())
                    .map(|row| row.as_slice()),
            )?
            .finish()?,
    )
}

/// Identify exact retained source bytes in the color-evidence source role.
///
/// The source schema domain and nonzero version describe the meaning of
/// `canonical_source`; they are identity-bearing rather than naming
/// conventions. The returned value is content-bound and explicitly unanchored.
///
/// # Errors
/// Refuses an empty domain, schema version zero, invalid resource limits,
/// budget overflow, or cancellation. No partial identity is published.
pub fn identify_color_evidence_source_v1<C>(
    source_domain: &str,
    source_schema_version: u32,
    canonical_source: &[u8],
    limits: EvidenceIdentityLimits,
    cancellation: C,
) -> Result<ColorEvidenceSourceV1, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    if source_domain.is_empty() {
        return Err(ColorEvidenceIdentityError::EmptySourceDomain);
    }
    if source_schema_version == 0 {
        return Err(ColorEvidenceIdentityError::ZeroSourceSchemaVersion);
    }
    let receipt = CanonicalEncoder::<ColorEvidenceSourceIdV1, _>::new(limits, cancellation)?
        .utf8(Field::new(0, "source-domain"), source_domain)?
        .u64(
            Field::new(1, "source-schema-version"),
            u64::from(source_schema_version),
        )?
        .bytes(Field::new(2, "canonical-source"), canonical_source)?
        .finish()?;
    Ok(ColorEvidenceSourceV1 { receipt })
}

/// Root one typed graph node at an exact retained source.
///
/// # Errors
/// Refuses malformed colors, resource overflow, invalid limits, or
/// cancellation. The result remains unanchored.
pub fn identify_color_evidence_source_node_v1<C>(
    source: &ColorEvidenceSourceV1,
    color: Color,
    limits: EvidenceIdentityLimits,
    cancellation: C,
) -> Result<ColorEvidenceNodeV1, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    let receipt = build_color_evidence_node_v1(
        ColorEvidenceOperationV1::Source,
        Some(source.id()),
        &color,
        None,
        limits,
        cancellation,
    )?;
    Ok(ColorEvidenceNodeV1 {
        color,
        receipt,
        operation: ColorEvidenceOperationV1::Source,
    })
}

/// Recompute one Add/Mul/Hull color result and identify the exact derivation.
///
/// Parent order is canonicalized by full typed ID before both color composition
/// and identity encoding, so commutative construction paths agree even where
/// legacy display-lineage strings were input-order-sensitive. Multiplicity is
/// retained. The opaque parent values prevent detaching an ID from its color.
///
/// # Errors
/// Refuses conflicting observations for one parent ID, malformed recomputed
/// output, resource overflow, invalid limits, or cancellation. No authority is
/// added.
pub fn compose_color_evidence_nodes_v1<C>(
    operation: ColorEvidenceCompositionOpV1,
    left: &ColorEvidenceNodeV1,
    right: &ColorEvidenceNodeV1,
    limits: EvidenceIdentityLimits,
    mut cancellation: C,
) -> Result<ColorEvidenceNodeV1, ColorEvidenceIdentityError>
where
    C: EvidenceIdentityCancellationProbe,
{
    poll_identity_cancellation(&mut cancellation)?;
    if matches!(
        adjudicate(
            ObservedIdentity::from_receipt(left.receipt()),
            ObservedIdentity::from_receipt(right.receipt()),
        ),
        IdentityAdjudication::Refused(_)
    ) {
        return Err(ColorEvidenceIdentityError::ParentObservationConflict);
    }

    let (first, second) = if left.id().as_bytes() <= right.id().as_bytes() {
        (left, right)
    } else {
        (right, left)
    };
    poll_identity_cancellation(&mut cancellation)?;
    let color = compose(
        first.color(),
        second.color(),
        operation.interval_operation(),
    );
    poll_identity_cancellation(&mut cancellation)?;
    let node_operation = operation.node_operation();
    let parents = [first.id(), second.id()];
    let receipt = build_color_evidence_node_v1(
        node_operation,
        None,
        &color,
        Some(parents),
        limits,
        cancellation,
    )?;
    Ok(ColorEvidenceNodeV1 {
        color,
        receipt,
        operation: node_operation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }

    impl EvidenceIdentityCancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }

    #[test]
    fn model_card_legacy_crosswalk_matches_public_fnv_and_cancels_by_stride() {
        let limits = EvidenceIdentityLimits::new(4096, 1024, 32, 64, 4);
        for bytes in [
            Vec::new(),
            vec![0x01, 0x02, 0x03],
            vec![0x01, 0x02, 0x03, 0x04],
            vec![0x01, 0x02, 0x03, 0x04, 0x05],
            vec![0x00, 0xff, 0x80, 0x7f, 0x42, 0x24, 0x11, 0xee, 0x99],
        ] {
            let mut never_cancel = || false;
            assert_eq!(
                model_card_legacy_provenance_v1(&bytes, limits, &mut never_cancel),
                Ok(ProvenanceHash::of_bytes(&bytes))
            );
        }

        let bytes = vec![0x5a; 12];
        let mut cancel = CancelAfter {
            successful_polls: 2,
        };
        assert_eq!(
            model_card_legacy_provenance_v1(&bytes, limits, &mut cancel),
            Err(ModelCardIdentityError::CalibrationCrosswalkCancelled { processed_bytes: 8 })
        );
    }

    #[test]
    fn ordered_card_field_arithmetic_straddles_the_hard_cap_exactly() {
        let empty_field = card_bound_model_card_field_bytes_v1(0).expect("empty field size");
        let root_bytes = 32_u64;
        let maximum_cards =
            (MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1 - empty_field) / root_bytes;
        let maximum_field =
            card_bound_model_card_field_bytes_v1(maximum_cards).expect("maximum field size");
        let first_oversized = card_bound_model_card_field_bytes_v1(maximum_cards + 1)
            .expect("first oversized field arithmetic");

        assert!(maximum_field <= MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1);
        assert!(first_oversized > MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1);
        assert_eq!(first_oversized - maximum_field, root_bytes);

        let caller_above_hard_cap =
            EvidenceIdentityLimits::new(1 << 22, 1 << 21, 32, maximum_cards + 1, 256);
        assert_eq!(
            preflight_card_bound_model_card_field_v1(maximum_cards, caller_above_hard_cap),
            Ok(maximum_field)
        );
        assert_eq!(
            preflight_card_bound_model_card_field_v1(maximum_cards + 1, caller_above_hard_cap),
            Err(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: first_oversized,
                limit: MAX_CARD_BOUND_MODEL_EVIDENCE_CARD_FIELD_BYTES_V1,
            })
        );
    }

    #[test]
    fn color_buffer_allocation_refusal_is_typed() {
        let mut cancellation = || false;
        let result = bounded_color_bytes_with_reservation(
            &Color::Verified { lo: 0.0, hi: 1.0 },
            EvidenceIdentityLimits::new(4096, 1024, 32, 64, 16),
            &mut cancellation,
            |output, capacity| {
                assert!(output.is_empty());
                assert_eq!(capacity, 34);
                Err(ColorBufferReservationError)
            },
        );

        assert_eq!(
            result,
            Err(ColorEvidenceIdentityError::ColorBufferAllocationFailed {
                requested_bytes: 34,
            })
        );
    }
}

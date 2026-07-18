//! Typed, canonical identities for evidence semantics.
//!
//! This module covers exact color-evidence graph replay, normalized validity
//! domains, model-card declarations with exact calibration sources, and an
//! opaque strong-identity projection of locally certified scalar evidence
//! through separate schemas. It does not reinterpret
//! [`crate::ProvenanceHash`], and it publishes only unanchored
//! [`IdentityReceipt`] values. Origin verification, policy admission,
//! structural [`crate::Certified`] consistency, and scientific color rank
//! remain separate axes.

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

use crate::{
    Ambition, COLOR_ALGEBRA_VERSION, Certified, Color, ColorPayloadError, IntervalOp, ModelCard,
    NumericalKind, ProvenanceHash, StatisticalCertificate, ValidityDomain, compose,
    validate_color_payload,
};

/// Identity schema version for exact retained color-evidence sources.
pub const COLOR_EVIDENCE_SOURCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for color-evidence graph nodes.
pub const COLOR_EVIDENCE_NODE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for normalized evidence-validity domains.
pub const VALIDITY_DOMAIN_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one locally certified scalar-evidence projection.
pub const CERTIFIED_F64_EVIDENCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for exact model-card calibration source bytes.
pub const MODEL_CARD_CALIBRATION_SOURCE_IDENTITY_VERSION_V1: u32 = 1;
/// Identity schema version for one helper-validated model-card declaration.
pub const MODEL_CARD_IDENTITY_VERSION_V1: u32 = 1;
/// Hard allocation ceiling for one canonical output-color payload.
pub const MAX_COLOR_EVIDENCE_NODE_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for the ordered axes field of one validity domain.
pub const MAX_VALIDITY_DOMAIN_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each variable certified-evidence field.
pub const MAX_CERTIFIED_F64_EVIDENCE_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard payload ceiling for each variable model-card declaration field.
pub const MAX_MODEL_CARD_IDENTITY_FIELD_BYTES_V1: u64 = 1 << 20;
/// Hard retained-byte ceiling for one model-card calibration source.
pub const MAX_MODEL_CARD_CALIBRATION_SOURCE_BYTES_V1: u64 = 1 << 20;
/// Non-semantic scatter/gather writes emitted for each streamed axis row.
const VALIDITY_DOMAIN_STREAM_CHUNKS_PER_AXIS_V1: u64 = 4;
/// Non-semantic scatter/gather writes emitted for each sensitivity row.
const CERTIFIED_F64_SENSITIVITY_STREAM_CHUNKS_PER_ROW_V1: u64 = 3;

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

const fn certified_f64_numerical_kind_tag_v1(kind: NumericalKind) -> u32 {
    match kind {
        NumericalKind::Exact => 1,
        NumericalKind::Enclosure => 2,
        NumericalKind::Estimate => 3,
        NumericalKind::NoClaim => 4,
    }
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
                certified_f64_numerical_kind_tag_v1(evidence.numerical.kind),
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

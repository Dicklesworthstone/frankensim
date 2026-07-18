//! Immutable, exact-byte crosswalks from quarantined legacy identities to
//! typed strong identities.
//!
//! This module keeps five concepts deliberately non-confusable:
//!
//! - exact legacy bytes and their plain BLAKE3 content ID;
//! - the historical FNV `u64`, retained only through its quarantine type;
//! - exact canonical payload bytes and their distinct plain content ID;
//! - a schema-typed semantic identity plus its canonical-frame audit record;
//! - authority state, which is never inferred from any digest or row.
//!
//! Multiple receipts may name the same legacy content ID. Bounded candidate
//! lookup exposes that ambiguity without selecting a winner; callers obtain a
//! typed semantic ID only by naming the exact expected Rust schema.

use std::{fmt, str};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, ContentId, Field, FieldSpec,
    IdentityAuditRecord, IdentityReceipt, IdentityRole, NoClaimState, SchemaId, SemanticId,
    StrongIdentity, TrustState, WireType, legacy::LegacyProvenanceV1,
};
use fsqlite::SqliteValue;

use crate::{
    ContentHash, EdgeRole, FiveExplicits, Ledger, LedgerError, MAX_TUNE_KERNEL_BYTES,
    MAX_TUNE_MACHINE_BYTES, MAX_TUNE_MEASURED_BYTES, MAX_TUNE_PARAMS_BYTES,
    MAX_TUNE_SHAPE_CLASS_BYTES, TuneRow, blob_param, now_wall_ns, row_text, sql_err, text_param,
    tune_corrupt,
};

/// Schema version of one artifact compatibility-hash to typed-content-ID row.
pub const ARTIFACT_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Schema version of one lineage-edge typed-content-ID companion row.
pub const EDGE_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Schema version of one frozen-operation-field typed-content-ID sidecar.
pub const OP_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Schema version of one autotuner-cache typed-content-ID sidecar.
pub const TUNE_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Maximum receipt IDs exposed by one artifact-semantic candidate lookup.
pub const MAX_ARTIFACT_SEMANTIC_BINDING_CANDIDATES: usize = 256;
/// Maximum UTF-8 bytes in an evidence name admitted to a semantic binding.
pub const MAX_EVIDENCE_SEMANTIC_BINDING_NAME_BYTES: usize = 64 * 1024;
/// Maximum receipt IDs exposed by one evidence-semantic candidate lookup.
pub const MAX_EVIDENCE_SEMANTIC_BINDING_CANDIDATES: usize = 256;

/// Identity schema version for immutable migration receipts.
pub const IDENTITY_MIGRATION_RECEIPT_VERSION: u32 = 1;
/// Canonical domain for immutable migration receipt identities.
pub const IDENTITY_MIGRATION_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-ledger.identity-migration-receipt.v1";
/// Maximum exact legacy or canonical payload retained by one receipt.
pub const MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES: usize = 1024 * 1024;
/// Maximum visible-ASCII semantic-rule bytes.
pub const MAX_IDENTITY_MIGRATION_RULE_BYTES: usize = 256;
/// Maximum static identity-domain bytes persisted by this schema.
pub const MAX_IDENTITY_MIGRATION_DOMAIN_BYTES: usize = 256;
/// Maximum static identity-schema-name bytes persisted by this schema.
pub const MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES: usize = 256;
/// Maximum static identity-context bytes persisted by this schema.
pub const MAX_IDENTITY_MIGRATION_CONTEXT_BYTES: usize = 4096;
/// Maximum receipt IDs exposed by one legacy-candidate lookup.
pub const MAX_IDENTITY_MIGRATION_CANDIDATES: usize = 256;
/// Version of the complete migration-receipt package/process transport.
pub const IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION: u32 = 1;
/// Maximum encoded bytes accepted by the complete migration-receipt transport.
pub const MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES: usize =
    IDENTITY_MIGRATION_RECEIPT_WIRE_BASE_BYTES
        + (2 * MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES)
        + MAX_IDENTITY_MIGRATION_RULE_BYTES
        + MAX_IDENTITY_MIGRATION_DOMAIN_BYTES
        + MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES
        + MAX_IDENTITY_MIGRATION_CONTEXT_BYTES
        + (3 * 32);

const IDENTITY_MIGRATION_RECEIPT_WIRE_MAGIC: &[u8; 8] = b"FSMIGR01";
const IDENTITY_MIGRATION_RECEIPT_WIRE_BASE_BYTES: usize = 290;

const RECEIPT_ID_LIMITS: CanonicalLimits = CanonicalLimits::new(64 * 1024, 8 * 1024, 64, 64, 4096);

/// Canonical schema for the immutable receipt over one legacy-to-strong-ID
/// crosswalk. Exact payloads are retained separately and bound here by their
/// plain content IDs and byte counts.
pub enum IdentityMigrationReceiptSchemaV1 {}

impl CanonicalSchema for IdentityMigrationReceiptSchemaV1 {
    const DOMAIN: &'static str = IDENTITY_MIGRATION_RECEIPT_DOMAIN;
    const NAME: &'static str = "identity-migration-receipt";
    const VERSION: u32 = IDENTITY_MIGRATION_RECEIPT_VERSION;
    const CONTEXT: &'static str = "exact legacy and canonical content roots, quarantined FNV, typed semantic schema, producer audit metadata, and explicit authority state; no inferred equivalence or trust";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("legacy-content-id", WireType::Bytes),
        FieldSpec::required("legacy-byte-count", WireType::U64),
        FieldSpec::required("legacy-fnv-le-u64", WireType::Bytes),
        FieldSpec::required("canonical-content-id", WireType::Bytes),
        FieldSpec::required("canonical-byte-count", WireType::U64),
        FieldSpec::required("semantic-rule", WireType::Utf8),
        FieldSpec::required("identity-role", WireType::Variant),
        FieldSpec::required("semantic-id", WireType::Bytes),
        FieldSpec::required("identity-domain", WireType::Utf8),
        FieldSpec::required("identity-schema-name", WireType::Utf8),
        FieldSpec::required("identity-schema-id", WireType::Bytes),
        FieldSpec::required("identity-schema-version", WireType::U64),
        FieldSpec::required("identity-context", WireType::Utf8),
        FieldSpec::required("canonical-preimage-content-id", WireType::Bytes),
        FieldSpec::required("canonical-frame-bytes", WireType::U64),
        FieldSpec::required("field-count", WireType::U64),
        FieldSpec::required("collection-items", WireType::U64),
        FieldSpec::required("max-canonical-bytes", WireType::U64),
        FieldSpec::required("max-field-bytes", WireType::U64),
        FieldSpec::required("max-fields", WireType::U64),
        FieldSpec::required("max-collection-items", WireType::U64),
        FieldSpec::required("cancellation-poll-bytes", WireType::U64),
        FieldSpec::required("trust-state", WireType::Variant),
        FieldSpec::optional_bytes("anchor-content-id"),
        FieldSpec::optional_bytes("verifier-id"),
        FieldSpec::optional_bytes("key-policy-id"),
        FieldSpec::required("no-claim-state", WireType::Variant),
    ];
}

/// Typed identity of one immutable migration receipt.
pub type IdentityMigrationReceiptId = SemanticId<IdentityMigrationReceiptSchemaV1>;

/// Structured refusal from the bounded migration-receipt wire codec.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdentityMigrationWireError {
    /// The whole transport exceeds its public allocation envelope.
    TransportLength { found: usize, max: usize },
    /// The fixed transport domain marker is absent or changed.
    Magic,
    /// The transport version is unknown to this decoder.
    UnsupportedVersion { found: u32 },
    /// A fixed or declared-length field runs past the available bytes.
    Truncated {
        field: &'static str,
        needed: usize,
        remaining: usize,
    },
    /// A declared variable field exceeds its field-specific bound.
    FieldLength {
        field: &'static str,
        found: usize,
        max: usize,
    },
    /// A required UTF-8 field contains invalid bytes.
    Utf8 { field: &'static str },
    /// A closed enum or optional-field tag is outside its declared universe.
    InvalidTag { field: &'static str, found: u8 },
    /// Extra bytes follow the complete v1 field sequence.
    TrailingBytes { remaining: usize },
    /// Complete decoded fields violate the receipt contract.
    InvalidReceipt { detail: String },
    /// The carried primary identity differs from independent reconstruction.
    ReceiptIdMismatch {
        stored: IdentityMigrationReceiptId,
        derived: IdentityMigrationReceiptId,
    },
}

impl fmt::Display for IdentityMigrationWireError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TransportLength { found, max } => {
                write!(formatter, "transport length {found} exceeds {max} bytes")
            }
            Self::Magic => formatter.write_str("migration receipt wire magic differs from v1"),
            Self::UnsupportedVersion { found } => {
                write!(
                    formatter,
                    "unsupported migration receipt wire version {found}"
                )
            }
            Self::Truncated {
                field,
                needed,
                remaining,
            } => write!(
                formatter,
                "wire field {field} needs {needed} bytes but only {remaining} remain"
            ),
            Self::FieldLength { field, found, max } => write!(
                formatter,
                "wire field {field} declares {found} bytes, exceeding {max}"
            ),
            Self::Utf8 { field } => write!(formatter, "wire field {field} is not UTF-8"),
            Self::InvalidTag { field, found } => {
                write!(formatter, "wire field {field} has unknown tag {found}")
            }
            Self::TrailingBytes { remaining } => {
                write!(formatter, "wire transport has {remaining} trailing bytes")
            }
            Self::InvalidReceipt { detail } => {
                write!(formatter, "decoded migration receipt is invalid: {detail}")
            }
            Self::ReceiptIdMismatch { stored, derived } => write!(
                formatter,
                "wire receipt ID {stored} differs from independently derived {derived}"
            ),
        }
    }
}

impl std::error::Error for IdentityMigrationWireError {}

/// Owner-local declaration consumed by `xtask check-identities`.
pub const IDENTITY_MIGRATION_RECEIPT_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:identity-migration-receipt",
    "version_const=IDENTITY_MIGRATION_RECEIPT_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.identity-migration-receipt.v1",
    "domain_const=IDENTITY_MIGRATION_RECEIPT_DOMAIN",
    "encoder=derive_receipt_id",
    "encoder_helpers=receipt_body_from_claim,validate_receipt_body",
    "schema_constants=IDENTITY_MIGRATION_RECEIPT_VERSION,IDENTITY_MIGRATION_RECEIPT_DOMAIN,IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION,MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES,MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,MAX_IDENTITY_MIGRATION_RULE_BYTES,MAX_IDENTITY_MIGRATION_DOMAIN_BYTES,MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES,MAX_IDENTITY_MIGRATION_CONTEXT_BYTES,crates/fs-ledger/src/schema.rs#V13",
    "schema_functions=derive_receipt_id,receipt_body_from_claim,validate_receipt_body,IdentityMigrationReceipt::typed_semantic_id,IdentityMigrationReceipt::to_wire_bytes,IdentityMigrationReceipt::from_wire_bytes,Ledger::record_identity_migration,Ledger::identity_migration_receipt,Ledger::identity_migration_candidates,crates/fs-blake3/src/identity.rs#CanonicalEncoder::finish",
    "schema_dependencies=fs-blake3:canonical-identity-frame,fs-blake3:schema-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=IdentityMigrationBody",
    "source_fields=IdentityMigrationBody.legacy_content_id:semantic,IdentityMigrationBody.legacy_fnv:semantic,IdentityMigrationBody.canonical_content_id:semantic,IdentityMigrationBody.semantic_rule:semantic,IdentityMigrationBody.semantic_id:semantic,IdentityMigrationBody.identity_role:semantic,IdentityMigrationBody.identity_domain:semantic,IdentityMigrationBody.identity_schema_name:semantic,IdentityMigrationBody.identity_schema_id:semantic,IdentityMigrationBody.identity_schema_version:semantic,IdentityMigrationBody.identity_context:semantic,IdentityMigrationBody.canonical_preimage_id:semantic,IdentityMigrationBody.canonical_frame_bytes:semantic,IdentityMigrationBody.field_count:semantic,IdentityMigrationBody.collection_items:semantic,IdentityMigrationBody.limits:semantic,IdentityMigrationBody.trust_state:semantic,IdentityMigrationBody.anchor_content_id:semantic,IdentityMigrationBody.verifier_id:semantic,IdentityMigrationBody.key_policy_id:semantic,IdentityMigrationBody.no_claim_state:semantic,IdentityMigrationBody.legacy_bytes:derived:bound-by-content-id-and-byte-count,IdentityMigrationBody.canonical_bytes:derived:bound-by-content-id-and-byte-count",
    "source_bindings=IdentityMigrationBody.legacy_bytes>legacy-content-id+legacy-byte-count,IdentityMigrationBody.legacy_fnv>legacy-fnv-le-u64,IdentityMigrationBody.canonical_bytes>canonical-content-id+canonical-byte-count,IdentityMigrationBody.semantic_rule>semantic-rule,IdentityMigrationBody.semantic_id>semantic-id,IdentityMigrationBody.identity_role>identity-role,IdentityMigrationBody.identity_domain>identity-domain,IdentityMigrationBody.identity_schema_name>identity-schema-name,IdentityMigrationBody.identity_schema_id>identity-schema-id,IdentityMigrationBody.identity_schema_version>identity-schema-version,IdentityMigrationBody.identity_context>identity-context,IdentityMigrationBody.canonical_preimage_id>canonical-preimage-content-id,IdentityMigrationBody.canonical_frame_bytes>canonical-frame-bytes,IdentityMigrationBody.field_count>field-count,IdentityMigrationBody.collection_items>collection-items,IdentityMigrationBody.limits>max-canonical-bytes+max-field-bytes+max-fields+max-collection-items+cancellation-poll-bytes,IdentityMigrationBody.trust_state>trust-state,IdentityMigrationBody.anchor_content_id>anchor-content-id,IdentityMigrationBody.verifier_id>verifier-id,IdentityMigrationBody.key_policy_id>key-policy-id,IdentityMigrationBody.no_claim_state>no-claim-state",
    "external_semantic_fields=receipt-schema-domain,receipt-schema-version,canonical-field-order",
    "semantic_fields=receipt-schema-domain,receipt-schema-version,canonical-field-order,legacy-content-id,legacy-byte-count,legacy-fnv-le-u64,canonical-content-id,canonical-byte-count,semantic-rule,identity-role,semantic-id,identity-domain,identity-schema-name,identity-schema-id,identity-schema-version,identity-context,canonical-preimage-content-id,canonical-frame-bytes,field-count,collection-items,max-canonical-bytes,max-field-bytes,max-fields,max-collection-items,cancellation-poll-bytes,trust-state,anchor-content-id,verifier-id,key-policy-id,no-claim-state",
    "excluded_fields=exact-legacy-bytes:bound-by-content-id,exact-canonical-bytes:bound-by-content-id,created-at:provenance-envelope-only",
    "consumers=Ledger::record_identity_migration,Ledger::identity_migration_receipt,Ledger::identity_migration_candidates,IdentityMigrationReceipt::typed_semantic_id,IdentityMigrationReceipt::to_wire_bytes,IdentityMigrationReceipt::from_wire_bytes,Ledger::bind_artifact_semantic_identity,Ledger::artifact_semantic_binding,Ledger::bind_evidence_semantic_identity,Ledger::evidence_semantic_binding",
    "mutations=receipt-schema-domain:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,receipt-schema-version:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-field-order:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-byte-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-fnv-le-u64:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-byte-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,semantic-rule:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-role:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,semantic-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-domain:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-name:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-version:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-context:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-preimage-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-frame-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,field-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,collection-items:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-canonical-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-field-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-fields:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-collection-items:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,cancellation-poll-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,trust-state:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,anchor-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,verifier-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,key-policy-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,no-claim-state:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state",
    "nonsemantic_mutations=created-at:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state",
    "field_guard=classify_identity_migration_receipt_fields",
    "transport_guard=IdentityMigrationReceipt::from_wire_bytes",
    "version_guard=crates/fs-ledger/tests/identity_migration.rs#wire_transport_refuses_truncation_extension_future_version_and_forged_id",
    "coupling_surface=fs-ledger:identity-migration-receipt",
];

/// Exact inputs for one immutable legacy-to-strong-ID crosswalk.
#[derive(Debug, Clone, Copy)]
pub struct IdentityMigrationClaim<'a, I: StrongIdentity> {
    /// Exact historical bytes. Empty is valid and remains distinguishable.
    pub legacy_bytes: &'a [u8],
    /// Exact historical FNV value, retained only in its quarantine type.
    pub legacy_fnv: LegacyProvenanceV1,
    /// Exact canonical payload bytes for the owner schema.
    pub canonical_bytes: &'a [u8],
    /// Stable visible-ASCII rule identifying the owner-defined transformation.
    pub semantic_rule: &'a str,
    /// Typed semantic identity producer receipt.
    pub receipt: IdentityReceipt<I>,
    /// Audit record derived from the receipt or an explicit authority ref.
    pub audit: IdentityAuditRecord,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct IdentityMigrationBody {
    legacy_bytes: Box<[u8]>,
    legacy_content_id: ContentId,
    legacy_fnv: LegacyProvenanceV1,
    canonical_bytes: Box<[u8]>,
    canonical_content_id: ContentId,
    semantic_rule: Box<str>,
    semantic_id: [u8; 32],
    identity_role: IdentityRole,
    identity_domain: Box<str>,
    identity_schema_name: Box<str>,
    identity_schema_id: [u8; 32],
    identity_schema_version: u32,
    identity_context: Box<str>,
    canonical_preimage_id: ContentId,
    canonical_frame_bytes: u64,
    field_count: u32,
    collection_items: u64,
    limits: CanonicalLimits,
    trust_state: TrustState,
    anchor_content_id: Option<ContentId>,
    verifier_id: Option<[u8; 32]>,
    key_policy_id: Option<[u8; 32]>,
    no_claim_state: NoClaimState,
}

/// Exhaustive owner-type classifier for identity governance. Adding a stored
/// receipt field must break this destructure until its identity role is
/// classified deliberately in the declaration above.
#[allow(dead_code)]
fn classify_identity_migration_receipt_fields(source: &IdentityMigrationBody) {
    let IdentityMigrationBody {
        legacy_bytes,
        legacy_content_id,
        legacy_fnv,
        canonical_bytes,
        canonical_content_id,
        semantic_rule,
        semantic_id,
        identity_role,
        identity_domain,
        identity_schema_name,
        identity_schema_id,
        identity_schema_version,
        identity_context,
        canonical_preimage_id,
        canonical_frame_bytes,
        field_count,
        collection_items,
        limits,
        trust_state,
        anchor_content_id,
        verifier_id,
        key_policy_id,
        no_claim_state,
    } = source;
    let _ = (
        legacy_bytes,
        legacy_content_id,
        legacy_fnv,
        canonical_bytes,
        canonical_content_id,
        semantic_rule,
        semantic_id,
        identity_role,
        identity_domain,
        identity_schema_name,
        identity_schema_id,
        identity_schema_version,
        identity_context,
        canonical_preimage_id,
        canonical_frame_bytes,
        field_count,
        collection_items,
        limits,
        trust_state,
        anchor_content_id,
        verifier_id,
        key_policy_id,
        no_claim_state,
    );
}

/// Independently reverified immutable migration receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityMigrationReceipt {
    receipt_id: IdentityMigrationReceiptId,
    body: IdentityMigrationBody,
}

impl IdentityMigrationReceipt {
    /// Typed identity of this complete crosswalk receipt.
    #[must_use]
    pub const fn receipt_id(&self) -> IdentityMigrationReceiptId {
        self.receipt_id
    }

    /// Exact historical bytes retained for replay and inspection.
    #[must_use]
    pub fn legacy_bytes(&self) -> &[u8] {
        &self.body.legacy_bytes
    }

    /// Plain BLAKE3 content ID of the exact historical bytes.
    #[must_use]
    pub const fn legacy_content_id(&self) -> ContentId {
        self.body.legacy_content_id
    }

    /// Quarantined historical FNV value; never widened into strong identity.
    #[must_use]
    pub const fn legacy_fnv(&self) -> LegacyProvenanceV1 {
        self.body.legacy_fnv
    }

    /// Exact canonical owner payload retained for replay and inspection.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.body.canonical_bytes
    }

    /// Plain BLAKE3 content ID of the exact canonical owner payload.
    #[must_use]
    pub const fn canonical_content_id(&self) -> ContentId {
        self.body.canonical_content_id
    }

    /// Exact owner-defined migration rule.
    #[must_use]
    pub fn semantic_rule(&self) -> &str {
        &self.body.semantic_rule
    }

    /// Stored semantic digest bytes. Use [`Self::typed_semantic_id`] before
    /// treating them as a semantic identity.
    #[must_use]
    pub const fn semantic_id_bytes(&self) -> [u8; 32] {
        self.body.semantic_id
    }

    /// Non-interchangeable identity role recorded by the producer.
    #[must_use]
    pub const fn identity_role(&self) -> IdentityRole {
        self.body.identity_role
    }

    /// Exact static identity domain.
    #[must_use]
    pub fn identity_domain(&self) -> &str {
        &self.body.identity_domain
    }

    /// Exact static identity schema name.
    #[must_use]
    pub fn identity_schema_name(&self) -> &str {
        &self.body.identity_schema_name
    }

    /// Exact schema descriptor digest.
    #[must_use]
    pub const fn identity_schema_id(&self) -> [u8; 32] {
        self.body.identity_schema_id
    }

    /// Exact semantic schema version.
    #[must_use]
    pub const fn identity_schema_version(&self) -> u32 {
        self.body.identity_schema_version
    }

    /// Exact static identity context.
    #[must_use]
    pub fn identity_context(&self) -> &str {
        &self.body.identity_context
    }

    /// Plain root of the strong identity's complete canonical frame.
    #[must_use]
    pub const fn canonical_preimage_id(&self) -> ContentId {
        self.body.canonical_preimage_id
    }

    /// Exact canonical-frame byte count from the producer receipt.
    #[must_use]
    pub const fn canonical_frame_bytes(&self) -> u64 {
        self.body.canonical_frame_bytes
    }

    /// Exact encoded top-level field count from the producer receipt.
    #[must_use]
    pub const fn field_count(&self) -> u32 {
        self.body.field_count
    }

    /// Exact encoded collection-item count from the producer receipt.
    #[must_use]
    pub const fn collection_items(&self) -> u64 {
        self.body.collection_items
    }

    /// Explicit resource envelope used by the semantic identity producer.
    #[must_use]
    pub const fn canonical_limits(&self) -> CanonicalLimits {
        self.body.limits
    }

    /// Trust state retained from the exact producer audit record.
    #[must_use]
    pub const fn trust_state(&self) -> TrustState {
        self.body.trust_state
    }

    /// Presented external anchor, when authority data exists.
    #[must_use]
    pub const fn anchor_content_id(&self) -> Option<ContentId> {
        self.body.anchor_content_id
    }

    /// Exact verifier identity bytes, when authority data exists.
    #[must_use]
    pub const fn verifier_id(&self) -> Option<[u8; 32]> {
        self.body.verifier_id
    }

    /// Exact key-policy identity bytes, when authority data exists.
    #[must_use]
    pub const fn key_policy_id(&self) -> Option<[u8; 32]> {
        self.body.key_policy_id
    }

    /// Explicit no-claim boundary retained from the producer audit record.
    #[must_use]
    pub const fn no_claim_state(&self) -> NoClaimState {
        self.body.no_claim_state
    }

    /// Parse the semantic digest only when every nominal role/schema field
    /// equals the caller's exact expected strong-identity type.
    #[must_use]
    pub fn typed_semantic_id<I: StrongIdentity>(&self) -> Option<I> {
        if self.body.identity_role != I::ROLE
            || self.body.identity_domain.as_ref() != I::Schema::DOMAIN
            || self.body.identity_schema_name.as_ref() != I::Schema::NAME
            || self.body.identity_schema_version != I::Schema::VERSION
            || self.body.identity_context.as_ref() != I::Schema::CONTEXT
            || self.body.identity_schema_id != *SchemaId::<I::Schema>::for_schema().as_bytes()
        {
            return None;
        }
        I::parse_slice(&self.body.semantic_id)
    }

    /// Encode every retained receipt field into the exact bounded v1
    /// package/process transport. Fixed-width integers are little-endian;
    /// variable fields carry explicit lengths and optional identities carry
    /// closed 0/1 tags.
    ///
    /// # Errors
    /// Refuses an internally inconsistent receipt or any field outside the
    /// published v1 transport bounds.
    pub fn to_wire_bytes(&self) -> Result<Vec<u8>, IdentityMigrationWireError> {
        encode_identity_migration_receipt(self)
    }

    /// Decode and independently reconstruct one complete transport candidate.
    ///
    /// Successful decoding proves only self-consistency. It does not grant
    /// ledger membership, semantic authority, verifier trust, or promotion.
    ///
    /// # Errors
    /// Refuses oversized, truncated, extended, future-version, invalid-UTF-8,
    /// unknown-tag, incoherent, content-mismatched, or forged-ID transports.
    pub fn from_wire_bytes(bytes: &[u8]) -> Result<Self, IdentityMigrationWireError> {
        decode_identity_migration_receipt(bytes)
    }
}

/// Result of an idempotent identity-migration write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IdentityMigrationWrite {
    receipt_id: IdentityMigrationReceiptId,
    legacy_content_id: ContentId,
    canonical_content_id: ContentId,
    deduped: bool,
}

impl IdentityMigrationWrite {
    /// Complete typed receipt identity.
    #[must_use]
    pub const fn receipt_id(self) -> IdentityMigrationReceiptId {
        self.receipt_id
    }

    /// Plain content ID of the retained historical bytes.
    #[must_use]
    pub const fn legacy_content_id(self) -> ContentId {
        self.legacy_content_id
    }

    /// Plain content ID of the retained canonical owner bytes.
    #[must_use]
    pub const fn canonical_content_id(self) -> ContentId {
        self.canonical_content_id
    }

    /// Whether an identical independently verified row already existed.
    #[must_use]
    pub const fn deduped(self) -> bool {
        self.deduped
    }
}

/// Bounded non-authoritative receipt IDs associated with one legacy content
/// ID. More than one candidate is expected and never resolved implicitly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityMigrationCandidates {
    receipt_ids: Vec<IdentityMigrationReceiptId>,
    truncated: bool,
}

/// Independently verified typed content identity for one retained artifact.
///
/// `artifact_hash` is the schema-v1 compatibility key. `content_id` is the
/// non-confusable raw-byte identity type introduced by the v14 sidecar. Their
/// digest bytes must be exactly equal; neither field carries semantic meaning
/// or authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArtifactContentIdentity {
    artifact_hash: ContentHash,
    content_id: ContentId,
    row_schema_version: u32,
}

impl ArtifactContentIdentity {
    /// Compatibility hash still used by legacy artifact and edge rows.
    #[must_use]
    pub const fn artifact_hash(self) -> ContentHash {
        self.artifact_hash
    }

    /// Typed plain-BLAKE3 identity of the exact retained artifact bytes.
    #[must_use]
    pub const fn content_id(self) -> ContentId {
        self.content_id
    }

    /// Exact sidecar row schema version.
    #[must_use]
    pub const fn row_schema_version(self) -> u32 {
        self.row_schema_version
    }
}

/// Independently verified typed artifact identity carried by one lineage edge.
///
/// Operation identity and edge role remain separate from the artifact's raw
/// byte identity. This row makes no claim about the artifact's semantic type
/// or the authority of the producing operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeContentIdentity {
    op: i64,
    role: EdgeRole,
    artifact_hash: ContentHash,
    content_id: ContentId,
    row_schema_version: u32,
}

/// Independently verified typed raw-content identities for one operation's
/// frozen input fields.
///
/// The ledger-local row ID, branch/mode, clocks, terminal outcome, diagnostic,
/// lineage, semantic meaning, and authority are intentionally not identities
/// in this sidecar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpContentIdentity {
    op: i64,
    session_content_id: Option<ContentId>,
    ir_content_id: ContentId,
    seed_content_id: ContentId,
    versions_content_id: ContentId,
    budget_content_id: ContentId,
    capability_content_id: ContentId,
    row_schema_version: u32,
}

impl OpContentIdentity {
    /// Ledger-local operation row owning these exact content identities.
    #[must_use]
    pub const fn op(self) -> i64 {
        self.op
    }

    /// Typed raw-content identity of the optional exact session bytes.
    #[must_use]
    pub const fn session_content_id(self) -> Option<ContentId> {
        self.session_content_id
    }

    /// Typed raw-content identity of the exact frozen IR JSON bytes.
    #[must_use]
    pub const fn ir_content_id(self) -> ContentId {
        self.ir_content_id
    }

    /// Typed raw-content identity of the exact frozen RNG seed bytes.
    #[must_use]
    pub const fn seed_content_id(self) -> ContentId {
        self.seed_content_id
    }

    /// Typed raw-content identity of the exact versions JSON bytes.
    #[must_use]
    pub const fn versions_content_id(self) -> ContentId {
        self.versions_content_id
    }

    /// Typed raw-content identity of the exact budget JSON bytes.
    #[must_use]
    pub const fn budget_content_id(self) -> ContentId {
        self.budget_content_id
    }

    /// Typed raw-content identity of the exact capability JSON bytes.
    #[must_use]
    pub const fn capability_content_id(self) -> ContentId {
        self.capability_content_id
    }

    /// Exact sidecar row schema version.
    #[must_use]
    pub const fn row_schema_version(self) -> u32 {
        self.row_schema_version
    }
}

/// Independently verified typed raw-content identities for one autotuner row.
///
/// Each exact cache-key and cache-value field receives its own plain content
/// identity. The sidecar does not claim that a kernel, shape, machine, params,
/// or measurement payload has an owner-defined semantic schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TuneContentIdentity {
    kernel_content_id: ContentId,
    shape_class_content_id: ContentId,
    machine_content_id: ContentId,
    params_content_id: ContentId,
    measured_content_id: ContentId,
    row_schema_version: u32,
}

impl TuneContentIdentity {
    /// Typed raw-content identity of the exact kernel bytes.
    #[must_use]
    pub const fn kernel_content_id(self) -> ContentId {
        self.kernel_content_id
    }

    /// Typed raw-content identity of the exact shape-class bytes.
    #[must_use]
    pub const fn shape_class_content_id(self) -> ContentId {
        self.shape_class_content_id
    }

    /// Typed raw-content identity of the exact machine-fingerprint bytes.
    #[must_use]
    pub const fn machine_content_id(self) -> ContentId {
        self.machine_content_id
    }

    /// Typed raw-content identity of the exact params JSON bytes.
    #[must_use]
    pub const fn params_content_id(self) -> ContentId {
        self.params_content_id
    }

    /// Typed raw-content identity of the exact measured JSON bytes.
    #[must_use]
    pub const fn measured_content_id(self) -> ContentId {
        self.measured_content_id
    }

    /// Exact sidecar row schema version.
    #[must_use]
    pub const fn row_schema_version(self) -> u32 {
        self.row_schema_version
    }
}

impl EdgeContentIdentity {
    /// Operation owning this role-qualified lineage edge.
    #[must_use]
    pub const fn op(self) -> i64 {
        self.op
    }

    /// Whether the operation consumes or produces the artifact.
    #[must_use]
    pub const fn role(self) -> EdgeRole {
        self.role
    }

    /// Schema-v1 artifact compatibility hash stored by the edge.
    #[must_use]
    pub const fn artifact_hash(self) -> ContentHash {
        self.artifact_hash
    }

    /// Typed raw-byte identity of the linked artifact.
    #[must_use]
    pub const fn content_id(self) -> ContentId {
        self.content_id
    }

    /// Exact sidecar row schema version.
    #[must_use]
    pub const fn row_schema_version(self) -> u32 {
        self.row_schema_version
    }
}

/// One immutable, independently reverified artifact-to-semantic binding.
///
/// The embedded receipt retains the exact nominal schema and authority state;
/// callers still choose the expected Rust type through
/// [`IdentityMigrationReceipt::typed_semantic_id`]. Row presence alone is
/// non-authoritative.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactSemanticBinding {
    artifact_hash: ContentHash,
    receipt: IdentityMigrationReceipt,
}

impl ArtifactSemanticBinding {
    /// Exact retained artifact named by this binding.
    #[must_use]
    pub const fn artifact_hash(&self) -> ContentHash {
        self.artifact_hash
    }

    /// Complete independently reverified migration receipt.
    #[must_use]
    pub const fn receipt(&self) -> &IdentityMigrationReceipt {
        &self.receipt
    }

    /// Project the semantic digest only for the caller's exact nominal type.
    #[must_use]
    pub fn typed_semantic_id<I: StrongIdentity>(&self) -> Option<I> {
        self.receipt.typed_semantic_id::<I>()
    }
}

/// Result of one idempotent artifact-semantic binding write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArtifactSemanticBindingWrite {
    artifact_hash: ContentHash,
    receipt_id: IdentityMigrationReceiptId,
    deduped: bool,
}

impl ArtifactSemanticBindingWrite {
    /// Exact retained artifact named by the binding.
    #[must_use]
    pub const fn artifact_hash(self) -> ContentHash {
        self.artifact_hash
    }

    /// Exact immutable migration receipt authorizing the semantic tuple.
    #[must_use]
    pub const fn receipt_id(self) -> IdentityMigrationReceiptId {
        self.receipt_id
    }

    /// Whether an identical independently verified binding already existed.
    #[must_use]
    pub const fn deduped(self) -> bool {
        self.deduped
    }
}

/// Bounded, deterministic, non-authoritative receipt candidates for one
/// artifact. Multiple candidates are preserved and never ranked implicitly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactSemanticBindingCandidates {
    receipt_ids: Vec<IdentityMigrationReceiptId>,
    truncated: bool,
}

impl ArtifactSemanticBindingCandidates {
    /// Deterministic receipt-ID prefix in bytewise order.
    #[must_use]
    pub fn receipt_ids(&self) -> &[IdentityMigrationReceiptId] {
        &self.receipt_ids
    }

    /// Whether at least one additional candidate exists beyond the caller cap.
    #[must_use]
    pub const fn truncated(&self) -> bool {
        self.truncated
    }
}

/// One immutable, independently reverified evidence-to-semantic binding.
///
/// The exact retained evidence JSON bytes must equal the receipt canonical
/// bytes. JSON equivalence, evidence names, and row presence never establish
/// semantic equivalence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceSemanticBinding {
    evidence_name: Box<str>,
    receipt: IdentityMigrationReceipt,
}

impl EvidenceSemanticBinding {
    /// Exact retained evidence record name.
    #[must_use]
    pub fn evidence_name(&self) -> &str {
        &self.evidence_name
    }

    /// Typed raw-byte identity of the exact retained evidence JSON.
    #[must_use]
    pub const fn content_id(&self) -> ContentId {
        self.receipt.canonical_content_id()
    }

    /// Exact immutable migration receipt authorizing the semantic tuple.
    #[must_use]
    pub const fn receipt_id(&self) -> IdentityMigrationReceiptId {
        self.receipt.receipt_id()
    }

    /// Complete independently reverified migration receipt.
    #[must_use]
    pub const fn receipt(&self) -> &IdentityMigrationReceipt {
        &self.receipt
    }

    /// Project the semantic digest only for the caller's exact nominal type.
    #[must_use]
    pub fn typed_semantic_id<I: StrongIdentity>(&self) -> Option<I> {
        self.receipt.typed_semantic_id::<I>()
    }
}

/// Result of one idempotent evidence-semantic binding write.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceSemanticBindingWrite {
    evidence_name: Box<str>,
    receipt_id: IdentityMigrationReceiptId,
    content_id: ContentId,
    deduped: bool,
}

impl EvidenceSemanticBindingWrite {
    /// Exact retained evidence record name.
    #[must_use]
    pub fn evidence_name(&self) -> &str {
        &self.evidence_name
    }

    /// Exact immutable migration receipt authorizing the semantic tuple.
    #[must_use]
    pub const fn receipt_id(&self) -> IdentityMigrationReceiptId {
        self.receipt_id
    }

    /// Typed raw-byte identity of the exact retained evidence JSON.
    #[must_use]
    pub const fn content_id(&self) -> ContentId {
        self.content_id
    }

    /// Whether an identical independently verified binding already existed.
    #[must_use]
    pub const fn deduped(&self) -> bool {
        self.deduped
    }
}

/// Bounded, deterministic, non-authoritative receipt candidates for one
/// evidence record. Multiple candidates are preserved and never ranked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceSemanticBindingCandidates {
    receipt_ids: Vec<IdentityMigrationReceiptId>,
    truncated: bool,
}

impl EvidenceSemanticBindingCandidates {
    /// Deterministic receipt-ID prefix in bytewise order.
    #[must_use]
    pub fn receipt_ids(&self) -> &[IdentityMigrationReceiptId] {
        &self.receipt_ids
    }

    /// Whether at least one additional candidate exists beyond the caller cap.
    #[must_use]
    pub const fn truncated(&self) -> bool {
        self.truncated
    }
}

impl IdentityMigrationCandidates {
    /// Deterministic receipt-ID prefix in bytewise order.
    #[must_use]
    pub fn receipt_ids(&self) -> &[IdentityMigrationReceiptId] {
        &self.receipt_ids
    }

    /// Whether at least one further candidate exists beyond the caller cap.
    #[must_use]
    pub const fn truncated(&self) -> bool {
        self.truncated
    }
}

fn invalid(field: &str, problem: impl Into<String>) -> LedgerError {
    LedgerError::Invalid {
        field: field.to_string(),
        problem: problem.into(),
    }
}

fn artifact_identity_corrupt(
    artifact_hash: &ContentHash,
    detail: impl Into<String>,
) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: artifact_hash.to_hex(),
        detail: detail.into(),
    }
}

fn edge_identity_corrupt(op: i64, detail: impl Into<String>) -> LedgerError {
    LedgerError::OpCorrupt {
        op,
        detail: detail.into(),
    }
}

fn op_identity_corrupt(op: i64, detail: impl Into<String>) -> LedgerError {
    LedgerError::OpCorrupt {
        op,
        detail: detail.into(),
    }
}

fn op_content_id(
    value: Option<&SqliteValue>,
    op: i64,
    field: &'static str,
) -> Result<ContentId, LedgerError> {
    match value {
        Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes),
        _ => None,
    }
    .ok_or_else(|| {
        op_identity_corrupt(
            op,
            format!("operation {field} is not an exact 32-byte typed content identity"),
        )
    })
}

fn optional_op_content_id(
    value: Option<&SqliteValue>,
    op: i64,
    field: &'static str,
) -> Result<Option<ContentId>, LedgerError> {
    match value {
        Some(SqliteValue::Null) => Ok(None),
        other => op_content_id(other, op, field).map(Some),
    }
}

fn derive_op_content_identity(
    op: i64,
    session: Option<&[u8]>,
    ir: &str,
    explicits: &FiveExplicits<'_>,
) -> OpContentIdentity {
    OpContentIdentity {
        op,
        session_content_id: session.map(ContentId::of_bytes),
        ir_content_id: ContentId::of_bytes(ir.as_bytes()),
        seed_content_id: ContentId::of_bytes(explicits.seed),
        versions_content_id: ContentId::of_bytes(explicits.versions.as_bytes()),
        budget_content_id: ContentId::of_bytes(explicits.budget.as_bytes()),
        capability_content_id: ContentId::of_bytes(explicits.capability.as_bytes()),
        row_schema_version: OP_CONTENT_IDENTITY_ROW_VERSION,
    }
}

fn tune_content_id(
    value: Option<&SqliteValue>,
    kernel: &str,
    field: &'static str,
) -> Result<ContentId, LedgerError> {
    match value {
        Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes),
        _ => None,
    }
    .ok_or_else(|| {
        tune_corrupt(
            kernel,
            format!("{field} is not an exact 32-byte typed content identity"),
        )
    })
}

fn derive_tune_content_identity(
    kernel: &str,
    shape_class: &str,
    machine: &[u8],
    params: &str,
    measured: &str,
) -> TuneContentIdentity {
    TuneContentIdentity {
        kernel_content_id: ContentId::of_bytes(kernel.as_bytes()),
        shape_class_content_id: ContentId::of_bytes(shape_class.as_bytes()),
        machine_content_id: ContentId::of_bytes(machine),
        params_content_id: ContentId::of_bytes(params.as_bytes()),
        measured_content_id: ContentId::of_bytes(measured.as_bytes()),
        row_schema_version: TUNE_CONTENT_IDENTITY_ROW_VERSION,
    }
}

fn artifact_semantic_binding_corrupt(
    artifact_hash: &ContentHash,
    receipt_id: IdentityMigrationReceiptId,
    detail: impl Into<String>,
) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: artifact_hash.to_hex(),
        detail: format!("artifact semantic binding {receipt_id}: {}", detail.into()),
    }
}

fn evidence_semantic_binding_corrupt(
    evidence_name: &str,
    receipt_id: IdentityMigrationReceiptId,
    detail: impl Into<String>,
) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: ContentId::of_bytes(evidence_name.as_bytes()).to_hex(),
        detail: format!("evidence semantic binding {receipt_id}: {}", detail.into()),
    }
}

fn validate_evidence_semantic_binding_name(evidence_name: &str) -> Result<(), LedgerError> {
    let len = evidence_name.len();
    if len == 0 || len > MAX_EVIDENCE_SEMANTIC_BINDING_NAME_BYTES {
        return Err(invalid(
            "evidence_semantic_binding.evidence_name",
            format!(
                "evidence name must contain 1..={MAX_EVIDENCE_SEMANTIC_BINDING_NAME_BYTES} UTF-8 bytes, found {len}"
            ),
        ));
    }
    Ok(())
}

fn stored_corrupt(
    receipt_id: IdentityMigrationReceiptId,
    detail: impl Into<String>,
) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: receipt_id.to_hex(),
        detail: format!("identity migration receipt is corrupt: {}", detail.into()),
    }
}

fn validate_rule(rule: &str) -> Result<(), &'static str> {
    if rule.is_empty() {
        return Err("semantic rule must be nonempty");
    }
    if rule.len() > MAX_IDENTITY_MIGRATION_RULE_BYTES {
        return Err("semantic rule exceeds 256 bytes");
    }
    if !rule.bytes().all(|byte| matches!(byte, b'!'..=b'~')) {
        return Err("semantic rule must contain visible ASCII without whitespace");
    }
    Ok(())
}

fn authority_is_coherent(body: &IdentityMigrationBody) -> bool {
    let complete = body.anchor_content_id.is_some()
        && body.verifier_id.is_some()
        && body.key_policy_id.is_some();
    let absent = body.anchor_content_id.is_none()
        && body.verifier_id.is_none()
        && body.key_policy_id.is_none();
    match (body.trust_state, body.no_claim_state) {
        (TrustState::Unanchored, NoClaimState::ExternalTrustRequired) => absent,
        (TrustState::Presented | TrustState::Verified, NoClaimState::ExternalTrustRequired) => {
            complete
        }
        (TrustState::Admitted, NoClaimState::ScientificCorrectnessNotProven) => complete,
        _ => false,
    }
}

fn validate_receipt_body(body: &IdentityMigrationBody) -> Result<(), &'static str> {
    if body.legacy_bytes.len() > MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES {
        return Err("legacy bytes exceed the 1 MiB receipt envelope");
    }
    if body.canonical_bytes.len() > MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES {
        return Err("canonical bytes exceed the 1 MiB receipt envelope");
    }
    validate_rule(&body.semantic_rule)?;
    if body.identity_domain.is_empty()
        || body.identity_domain.len() > MAX_IDENTITY_MIGRATION_DOMAIN_BYTES
    {
        return Err("identity domain is empty or exceeds 256 bytes");
    }
    if body.identity_schema_name.is_empty()
        || body.identity_schema_name.len() > MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES
    {
        return Err("identity schema name is empty or exceeds 256 bytes");
    }
    if body.identity_context.is_empty()
        || body.identity_context.len() > MAX_IDENTITY_MIGRATION_CONTEXT_BYTES
    {
        return Err("identity context is empty or exceeds 4096 bytes");
    }
    if body.identity_schema_version == 0 {
        return Err("identity schema version zero is invalid");
    }
    if ContentId::of_bytes(&body.legacy_bytes) != body.legacy_content_id {
        return Err("legacy content ID does not match the exact retained bytes");
    }
    if ContentId::of_bytes(&body.canonical_bytes) != body.canonical_content_id {
        return Err("canonical content ID does not match the exact retained bytes");
    }
    if body.canonical_frame_bytes > body.limits.max_canonical_bytes() {
        return Err("canonical frame byte count exceeds producer limit");
    }
    if body.field_count > body.limits.max_fields() {
        return Err("encoded field count exceeds producer limit");
    }
    if body.collection_items > body.limits.max_collection_items() {
        return Err("encoded collection count exceeds producer limit");
    }
    if body.limits.max_canonical_bytes() == 0
        || body.limits.max_field_bytes() == 0
        || body.limits.max_fields() == 0
        || body.limits.cancellation_poll_bytes() == 0
    {
        return Err("producer canonical limits contain a zero resource bound");
    }
    if !authority_is_coherent(body) {
        return Err("trust, authority, and no-claim fields are incoherent");
    }
    Ok(())
}

fn receipt_body_from_claim<I: StrongIdentity>(
    claim: IdentityMigrationClaim<'_, I>,
) -> Result<IdentityMigrationBody, LedgerError> {
    let base = claim.receipt.audit_record();
    let audit = claim.audit;
    if base.id() != audit.id()
        || base.canonical_preimage() != audit.canonical_preimage()
        || base.role() != audit.role()
        || base.domain() != audit.domain()
        || base.schema_name() != audit.schema_name()
        || base.schema_id() != audit.schema_id()
        || base.version() != audit.version()
        || base.context() != audit.context()
        || base.canonical_bytes() != audit.canonical_bytes()
        || base.field_count() != audit.field_count()
        || base.collection_items() != audit.collection_items()
        || base.limits() != audit.limits()
    {
        return Err(invalid(
            "identity_migration.audit",
            "audit record does not describe the exact offered typed receipt",
        ));
    }
    let body = IdentityMigrationBody {
        legacy_bytes: claim.legacy_bytes.into(),
        legacy_content_id: ContentId::of_bytes(claim.legacy_bytes),
        legacy_fnv: claim.legacy_fnv,
        canonical_bytes: claim.canonical_bytes.into(),
        canonical_content_id: ContentId::of_bytes(claim.canonical_bytes),
        semantic_rule: claim.semantic_rule.into(),
        semantic_id: audit.id(),
        identity_role: audit.role(),
        identity_domain: audit.domain().into(),
        identity_schema_name: audit.schema_name().into(),
        identity_schema_id: audit.schema_id(),
        identity_schema_version: audit.version(),
        identity_context: audit.context().into(),
        canonical_preimage_id: audit.canonical_preimage(),
        canonical_frame_bytes: audit.canonical_bytes(),
        field_count: audit.field_count(),
        collection_items: audit.collection_items(),
        limits: audit.limits(),
        trust_state: audit.trust(),
        anchor_content_id: audit.anchor(),
        verifier_id: audit.verifier(),
        key_policy_id: audit.key_policy(),
        no_claim_state: audit.no_claim(),
    };
    validate_receipt_body(&body).map_err(|problem| invalid("identity_migration.claim", problem))?;
    Ok(body)
}

fn derive_receipt_id(
    body: &IdentityMigrationBody,
) -> Result<IdentityMigrationReceiptId, fs_blake3::identity::CanonicalError> {
    let legacy_fnv = body.legacy_fnv.value().to_le_bytes();
    let anchor = body.anchor_content_id.map(|value| *value.as_bytes());
    let legacy_byte_count = u64::try_from(body.legacy_bytes.len())
        .map_err(|_| fs_blake3::identity::CanonicalError::LengthOverflow)?;
    let canonical_byte_count = u64::try_from(body.canonical_bytes.len())
        .map_err(|_| fs_blake3::identity::CanonicalError::LengthOverflow)?;
    CanonicalEncoder::<IdentityMigrationReceiptId, _>::new(RECEIPT_ID_LIMITS, || false)?
        .bytes(
            Field::new(0, "legacy-content-id"),
            body.legacy_content_id.as_bytes(),
        )?
        .u64(Field::new(1, "legacy-byte-count"), legacy_byte_count)?
        .bytes(Field::new(2, "legacy-fnv-le-u64"), &legacy_fnv)?
        .bytes(
            Field::new(3, "canonical-content-id"),
            body.canonical_content_id.as_bytes(),
        )?
        .u64(Field::new(4, "canonical-byte-count"), canonical_byte_count)?
        .utf8(Field::new(5, "semantic-rule"), &body.semantic_rule)?
        .variant(
            Field::new(6, "identity-role"),
            u32::from(body.identity_role.tag()),
            &[],
        )?
        .bytes(Field::new(7, "semantic-id"), &body.semantic_id)?
        .utf8(Field::new(8, "identity-domain"), &body.identity_domain)?
        .utf8(
            Field::new(9, "identity-schema-name"),
            &body.identity_schema_name,
        )?
        .bytes(
            Field::new(10, "identity-schema-id"),
            &body.identity_schema_id,
        )?
        .u64(
            Field::new(11, "identity-schema-version"),
            u64::from(body.identity_schema_version),
        )?
        .utf8(Field::new(12, "identity-context"), &body.identity_context)?
        .bytes(
            Field::new(13, "canonical-preimage-content-id"),
            body.canonical_preimage_id.as_bytes(),
        )?
        .u64(
            Field::new(14, "canonical-frame-bytes"),
            body.canonical_frame_bytes,
        )?
        .u64(Field::new(15, "field-count"), u64::from(body.field_count))?
        .u64(Field::new(16, "collection-items"), body.collection_items)?
        .u64(
            Field::new(17, "max-canonical-bytes"),
            body.limits.max_canonical_bytes(),
        )?
        .u64(
            Field::new(18, "max-field-bytes"),
            body.limits.max_field_bytes(),
        )?
        .u64(
            Field::new(19, "max-fields"),
            u64::from(body.limits.max_fields()),
        )?
        .u64(
            Field::new(20, "max-collection-items"),
            body.limits.max_collection_items(),
        )?
        .u64(
            Field::new(21, "cancellation-poll-bytes"),
            u64::from(body.limits.cancellation_poll_bytes()),
        )?
        .variant(
            Field::new(22, "trust-state"),
            trust_state_tag(body.trust_state),
            &[],
        )?
        .optional_bytes(
            Field::new(23, "anchor-content-id"),
            anchor.as_ref().map(<[u8; 32]>::as_slice),
        )?
        .optional_bytes(
            Field::new(24, "verifier-id"),
            body.verifier_id.as_ref().map(<[u8; 32]>::as_slice),
        )?
        .optional_bytes(
            Field::new(25, "key-policy-id"),
            body.key_policy_id.as_ref().map(<[u8; 32]>::as_slice),
        )?
        .variant(
            Field::new(26, "no-claim-state"),
            no_claim_state_tag(body.no_claim_state),
            &[],
        )?
        .finish()
        .map(|receipt| receipt.id())
}

const fn trust_state_tag(state: TrustState) -> u32 {
    match state {
        TrustState::Unanchored => 0,
        TrustState::Presented => 1,
        TrustState::Verified => 2,
        TrustState::Admitted => 3,
    }
}

const fn no_claim_state_tag(state: NoClaimState) -> u32 {
    match state {
        NoClaimState::ExternalTrustRequired => 0,
        NoClaimState::ScientificCorrectnessNotProven => 1,
    }
}

fn identity_role_from_tag(tag: i64) -> Option<IdentityRole> {
    match tag {
        1 => Some(IdentityRole::Semantic),
        2 => Some(IdentityRole::WireContent),
        3 => Some(IdentityRole::EvidenceNode),
        4 => Some(IdentityRole::Entity),
        5 => Some(IdentityRole::SourceBytes),
        6 => Some(IdentityRole::Source),
        7 => Some(IdentityRole::Model),
        8 => Some(IdentityRole::Checker),
        9 => Some(IdentityRole::Schema),
        10 => Some(IdentityRole::Verifier),
        11 => Some(IdentityRole::KeyPolicy),
        12 => Some(IdentityRole::ProblemSemantic),
        _ => None,
    }
}

fn trust_state_from_tag(tag: i64) -> Option<TrustState> {
    match tag {
        0 => Some(TrustState::Unanchored),
        1 => Some(TrustState::Presented),
        2 => Some(TrustState::Verified),
        3 => Some(TrustState::Admitted),
        _ => None,
    }
}

fn no_claim_state_from_tag(tag: i64) -> Option<NoClaimState> {
    match tag {
        0 => Some(NoClaimState::ExternalTrustRequired),
        1 => Some(NoClaimState::ScientificCorrectnessNotProven),
        _ => None,
    }
}

fn wire_invalid(detail: impl Into<String>) -> IdentityMigrationWireError {
    IdentityMigrationWireError::InvalidReceipt {
        detail: detail.into(),
    }
}

fn push_wire_u16_text(
    output: &mut Vec<u8>,
    field: &'static str,
    value: &str,
    max: usize,
) -> Result<(), IdentityMigrationWireError> {
    let len = value.len();
    if len > max {
        return Err(IdentityMigrationWireError::FieldLength {
            field,
            found: len,
            max,
        });
    }
    let len = u16::try_from(len).map_err(|_| IdentityMigrationWireError::FieldLength {
        field,
        found: len,
        max,
    })?;
    output.extend_from_slice(&len.to_le_bytes());
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn push_wire_u32_bytes(
    output: &mut Vec<u8>,
    field: &'static str,
    value: &[u8],
    max: usize,
) -> Result<(), IdentityMigrationWireError> {
    let len = value.len();
    if len > max {
        return Err(IdentityMigrationWireError::FieldLength {
            field,
            found: len,
            max,
        });
    }
    let len = u32::try_from(len).map_err(|_| IdentityMigrationWireError::FieldLength {
        field,
        found: len,
        max,
    })?;
    output.extend_from_slice(&len.to_le_bytes());
    output.extend_from_slice(value);
    Ok(())
}

fn push_wire_optional_fixed(output: &mut Vec<u8>, value: Option<&[u8; 32]>) {
    match value {
        Some(value) => {
            output.push(1);
            output.extend_from_slice(value);
        }
        None => output.push(0),
    }
}

fn encode_identity_migration_receipt(
    receipt: &IdentityMigrationReceipt,
) -> Result<Vec<u8>, IdentityMigrationWireError> {
    validate_receipt_body(&receipt.body).map_err(wire_invalid)?;
    let derived = derive_receipt_id(&receipt.body)
        .map_err(|error| wire_invalid(format!("receipt ID reconstruction refused: {error}")))?;
    if derived != receipt.receipt_id {
        return Err(IdentityMigrationWireError::ReceiptIdMismatch {
            stored: receipt.receipt_id,
            derived,
        });
    }

    let body = &receipt.body;
    let capacity = IDENTITY_MIGRATION_RECEIPT_WIRE_BASE_BYTES
        .saturating_add(body.legacy_bytes.len())
        .saturating_add(body.canonical_bytes.len())
        .saturating_add(body.semantic_rule.len())
        .saturating_add(body.identity_domain.len())
        .saturating_add(body.identity_schema_name.len())
        .saturating_add(body.identity_context.len())
        .saturating_add(body.anchor_content_id.map_or(0, |_| 32))
        .saturating_add(body.verifier_id.map_or(0, |_| 32))
        .saturating_add(body.key_policy_id.map_or(0, |_| 32));
    if capacity > MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES {
        return Err(IdentityMigrationWireError::TransportLength {
            found: capacity,
            max: MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES,
        });
    }

    let mut output = Vec::with_capacity(capacity);
    output.extend_from_slice(IDENTITY_MIGRATION_RECEIPT_WIRE_MAGIC);
    output.extend_from_slice(&IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION.to_le_bytes());
    output.extend_from_slice(receipt.receipt_id.as_bytes());
    push_wire_u32_bytes(
        &mut output,
        "legacy_bytes",
        &body.legacy_bytes,
        MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
    )?;
    output.extend_from_slice(body.legacy_content_id.as_bytes());
    output.extend_from_slice(&body.legacy_fnv.value().to_le_bytes());
    push_wire_u32_bytes(
        &mut output,
        "canonical_bytes",
        &body.canonical_bytes,
        MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
    )?;
    output.extend_from_slice(body.canonical_content_id.as_bytes());
    push_wire_u16_text(
        &mut output,
        "semantic_rule",
        &body.semantic_rule,
        MAX_IDENTITY_MIGRATION_RULE_BYTES,
    )?;
    output.extend_from_slice(&body.semantic_id);
    output.push(body.identity_role.tag());
    push_wire_u16_text(
        &mut output,
        "identity_domain",
        &body.identity_domain,
        MAX_IDENTITY_MIGRATION_DOMAIN_BYTES,
    )?;
    push_wire_u16_text(
        &mut output,
        "identity_schema_name",
        &body.identity_schema_name,
        MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES,
    )?;
    output.extend_from_slice(&body.identity_schema_id);
    output.extend_from_slice(&body.identity_schema_version.to_le_bytes());
    push_wire_u16_text(
        &mut output,
        "identity_context",
        &body.identity_context,
        MAX_IDENTITY_MIGRATION_CONTEXT_BYTES,
    )?;
    output.extend_from_slice(body.canonical_preimage_id.as_bytes());
    output.extend_from_slice(&body.canonical_frame_bytes.to_le_bytes());
    output.extend_from_slice(&body.field_count.to_le_bytes());
    output.extend_from_slice(&body.collection_items.to_le_bytes());
    output.extend_from_slice(&body.limits.max_canonical_bytes().to_le_bytes());
    output.extend_from_slice(&body.limits.max_field_bytes().to_le_bytes());
    output.extend_from_slice(&body.limits.max_fields().to_le_bytes());
    output.extend_from_slice(&body.limits.max_collection_items().to_le_bytes());
    output.extend_from_slice(&body.limits.cancellation_poll_bytes().to_le_bytes());
    output.push(
        u8::try_from(trust_state_tag(body.trust_state))
            .map_err(|_| wire_invalid("trust-state tag exceeds one byte"))?,
    );
    let anchor = body.anchor_content_id.as_ref().map(ContentId::as_bytes);
    push_wire_optional_fixed(&mut output, anchor);
    push_wire_optional_fixed(&mut output, body.verifier_id.as_ref());
    push_wire_optional_fixed(&mut output, body.key_policy_id.as_ref());
    output.push(
        u8::try_from(no_claim_state_tag(body.no_claim_state))
            .map_err(|_| wire_invalid("no-claim-state tag exceeds one byte"))?,
    );
    if output.len() != capacity {
        return Err(wire_invalid(format!(
            "wire length accounting expected {capacity} bytes but encoded {}",
            output.len()
        )));
    }
    Ok(output)
}

struct IdentityMigrationWireCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> IdentityMigrationWireCursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(
        &mut self,
        len: usize,
        field: &'static str,
    ) -> Result<&'a [u8], IdentityMigrationWireError> {
        let remaining = self.bytes.len().saturating_sub(self.offset);
        let Some(end) = self.offset.checked_add(len) else {
            return Err(IdentityMigrationWireError::Truncated {
                field,
                needed: len,
                remaining,
            });
        };
        let Some(value) = self.bytes.get(self.offset..end) else {
            return Err(IdentityMigrationWireError::Truncated {
                field,
                needed: len,
                remaining,
            });
        };
        self.offset = end;
        Ok(value)
    }

    fn fixed<const N: usize>(
        &mut self,
        field: &'static str,
    ) -> Result<[u8; N], IdentityMigrationWireError> {
        let bytes = self.take(N, field)?;
        let mut value = [0_u8; N];
        value.copy_from_slice(bytes);
        Ok(value)
    }

    fn u8(&mut self, field: &'static str) -> Result<u8, IdentityMigrationWireError> {
        self.take(1, field)?
            .first()
            .copied()
            .ok_or(IdentityMigrationWireError::Truncated {
                field,
                needed: 1,
                remaining: 0,
            })
    }

    fn u16(&mut self, field: &'static str) -> Result<u16, IdentityMigrationWireError> {
        Ok(u16::from_le_bytes(self.fixed(field)?))
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, IdentityMigrationWireError> {
        Ok(u32::from_le_bytes(self.fixed(field)?))
    }

    fn u64(&mut self, field: &'static str) -> Result<u64, IdentityMigrationWireError> {
        Ok(u64::from_le_bytes(self.fixed(field)?))
    }

    fn u32_bytes(
        &mut self,
        field: &'static str,
        max: usize,
    ) -> Result<Box<[u8]>, IdentityMigrationWireError> {
        let found = usize::try_from(self.u32(field)?).map_err(|_| {
            IdentityMigrationWireError::FieldLength {
                field,
                found: usize::MAX,
                max,
            }
        })?;
        if found > max {
            return Err(IdentityMigrationWireError::FieldLength { field, found, max });
        }
        Ok(self.take(found, field)?.into())
    }

    fn u16_text(
        &mut self,
        field: &'static str,
        max: usize,
    ) -> Result<Box<str>, IdentityMigrationWireError> {
        let found = usize::from(self.u16(field)?);
        if found > max {
            return Err(IdentityMigrationWireError::FieldLength { field, found, max });
        }
        let bytes = self.take(found, field)?;
        let value =
            str::from_utf8(bytes).map_err(|_| IdentityMigrationWireError::Utf8 { field })?;
        Ok(value.into())
    }

    fn optional_fixed(
        &mut self,
        field: &'static str,
    ) -> Result<Option<[u8; 32]>, IdentityMigrationWireError> {
        match self.u8(field)? {
            0 => Ok(None),
            1 => self.fixed(field).map(Some),
            found => Err(IdentityMigrationWireError::InvalidTag { field, found }),
        }
    }

    fn finish(self) -> Result<(), IdentityMigrationWireError> {
        let remaining = self.bytes.len().saturating_sub(self.offset);
        if remaining == 0 {
            Ok(())
        } else {
            Err(IdentityMigrationWireError::TrailingBytes { remaining })
        }
    }
}

fn wire_content_id(
    bytes: [u8; 32],
    field: &'static str,
) -> Result<ContentId, IdentityMigrationWireError> {
    ContentId::parse_slice(&bytes)
        .ok_or_else(|| wire_invalid(format!("{field} is not a typed 32-byte content ID")))
}

fn decode_identity_migration_receipt(
    bytes: &[u8],
) -> Result<IdentityMigrationReceipt, IdentityMigrationWireError> {
    if bytes.len() > MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES {
        return Err(IdentityMigrationWireError::TransportLength {
            found: bytes.len(),
            max: MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES,
        });
    }
    let mut cursor = IdentityMigrationWireCursor::new(bytes);
    if cursor.fixed::<8>("magic")? != *IDENTITY_MIGRATION_RECEIPT_WIRE_MAGIC {
        return Err(IdentityMigrationWireError::Magic);
    }
    let version = cursor.u32("wire_version")?;
    if version != IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION {
        return Err(IdentityMigrationWireError::UnsupportedVersion { found: version });
    }
    let stored_id_bytes = cursor.fixed::<32>("receipt_id")?;
    let stored_id = IdentityMigrationReceiptId::parse_slice(&stored_id_bytes)
        .ok_or_else(|| wire_invalid("receipt_id is not a typed 32-byte identity"))?;
    let legacy_bytes = cursor.u32_bytes("legacy_bytes", MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES)?;
    let legacy_content_id =
        wire_content_id(cursor.fixed("legacy_content_id")?, "legacy_content_id")?;
    let legacy_fnv = LegacyProvenanceV1::new(cursor.u64("legacy_fnv")?);
    let canonical_bytes =
        cursor.u32_bytes("canonical_bytes", MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES)?;
    let canonical_content_id = wire_content_id(
        cursor.fixed("canonical_content_id")?,
        "canonical_content_id",
    )?;
    let semantic_rule = cursor.u16_text("semantic_rule", MAX_IDENTITY_MIGRATION_RULE_BYTES)?;
    let semantic_id = cursor.fixed("semantic_id")?;
    let identity_role_tag = cursor.u8("identity_role")?;
    let identity_role = identity_role_from_tag(i64::from(identity_role_tag)).ok_or(
        IdentityMigrationWireError::InvalidTag {
            field: "identity_role",
            found: identity_role_tag,
        },
    )?;
    let identity_domain =
        cursor.u16_text("identity_domain", MAX_IDENTITY_MIGRATION_DOMAIN_BYTES)?;
    let identity_schema_name = cursor.u16_text(
        "identity_schema_name",
        MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES,
    )?;
    let identity_schema_id = cursor.fixed("identity_schema_id")?;
    let identity_schema_version = cursor.u32("identity_schema_version")?;
    let identity_context =
        cursor.u16_text("identity_context", MAX_IDENTITY_MIGRATION_CONTEXT_BYTES)?;
    let canonical_preimage_id = wire_content_id(
        cursor.fixed("canonical_preimage_id")?,
        "canonical_preimage_id",
    )?;
    let canonical_frame_bytes = cursor.u64("canonical_frame_bytes")?;
    let field_count = cursor.u32("field_count")?;
    let collection_items = cursor.u64("collection_items")?;
    let limits = CanonicalLimits::new(
        cursor.u64("max_canonical_bytes")?,
        cursor.u64("max_field_bytes")?,
        cursor.u32("max_fields")?,
        cursor.u64("max_collection_items")?,
        cursor.u32("cancellation_poll_bytes")?,
    );
    let trust_tag = cursor.u8("trust_state")?;
    let trust_state = trust_state_from_tag(i64::from(trust_tag)).ok_or(
        IdentityMigrationWireError::InvalidTag {
            field: "trust_state",
            found: trust_tag,
        },
    )?;
    let anchor_content_id = cursor
        .optional_fixed("anchor_content_id")?
        .map(|value| wire_content_id(value, "anchor_content_id"))
        .transpose()?;
    let verifier_id = cursor.optional_fixed("verifier_id")?;
    let key_policy_id = cursor.optional_fixed("key_policy_id")?;
    let no_claim_tag = cursor.u8("no_claim_state")?;
    let no_claim_state = no_claim_state_from_tag(i64::from(no_claim_tag)).ok_or(
        IdentityMigrationWireError::InvalidTag {
            field: "no_claim_state",
            found: no_claim_tag,
        },
    )?;
    cursor.finish()?;

    let body = IdentityMigrationBody {
        legacy_bytes,
        legacy_content_id,
        legacy_fnv,
        canonical_bytes,
        canonical_content_id,
        semantic_rule,
        semantic_id,
        identity_role,
        identity_domain,
        identity_schema_name,
        identity_schema_id,
        identity_schema_version,
        identity_context,
        canonical_preimage_id,
        canonical_frame_bytes,
        field_count,
        collection_items,
        limits,
        trust_state,
        anchor_content_id,
        verifier_id,
        key_policy_id,
        no_claim_state,
    };
    validate_receipt_body(&body).map_err(wire_invalid)?;
    let derived = derive_receipt_id(&body)
        .map_err(|error| wire_invalid(format!("receipt ID reconstruction refused: {error}")))?;
    if stored_id != derived {
        return Err(IdentityMigrationWireError::ReceiptIdMismatch {
            stored: stored_id,
            derived,
        });
    }
    Ok(IdentityMigrationReceipt {
        receipt_id: stored_id,
        body,
    })
}

fn fixed_bytes<const N: usize>(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<[u8; N], LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = value else {
        return Err(stored_corrupt(receipt_id, format!("{field} is not a BLOB")));
    };
    bytes.as_slice().try_into().map_err(|_| {
        stored_corrupt(
            receipt_id,
            format!(
                "{field} must contain exactly {N} bytes, found {}",
                bytes.len()
            ),
        )
    })
}

fn bounded_bytes(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
    max: usize,
) -> Result<Box<[u8]>, LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = value else {
        return Err(stored_corrupt(receipt_id, format!("{field} is not a BLOB")));
    };
    if bytes.len() > max {
        return Err(stored_corrupt(
            receipt_id,
            format!("{field} exceeds its {max}-byte storage envelope"),
        ));
    }
    Ok(bytes.clone().into_boxed_slice())
}

fn bounded_utf8(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
    max: usize,
) -> Result<Box<str>, LedgerError> {
    let bytes = bounded_bytes(value, receipt_id, field, max)?;
    let text = str::from_utf8(&bytes)
        .map_err(|_| stored_corrupt(receipt_id, format!("{field} is not UTF-8")))?;
    Ok(text.to_string().into_boxed_str())
}

fn integer(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<i64, LedgerError> {
    let Some(SqliteValue::Integer(value)) = value else {
        return Err(stored_corrupt(
            receipt_id,
            format!("{field} is not an INTEGER"),
        ));
    };
    Ok(*value)
}

fn u32_integer(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<u32, LedgerError> {
    u32::try_from(integer(value, receipt_id, field)?)
        .map_err(|_| stored_corrupt(receipt_id, format!("{field} is outside the u32 domain")))
}

fn u64_blob(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<u64, LedgerError> {
    Ok(u64::from_le_bytes(fixed_bytes::<8>(
        value, receipt_id, field,
    )?))
}

fn content_id(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<ContentId, LedgerError> {
    let bytes = fixed_bytes::<32>(value, receipt_id, field)?;
    ContentId::parse_slice(&bytes).ok_or_else(|| {
        stored_corrupt(
            receipt_id,
            format!("{field} is not a 32-byte content identity"),
        )
    })
}

fn optional_content_id(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<Option<ContentId>, LedgerError> {
    match value {
        Some(SqliteValue::Null) => Ok(None),
        other => content_id(other, receipt_id, field).map(Some),
    }
}

fn optional_fixed_32(
    value: Option<&SqliteValue>,
    receipt_id: IdentityMigrationReceiptId,
    field: &'static str,
) -> Result<Option<[u8; 32]>, LedgerError> {
    match value {
        Some(SqliteValue::Null) => Ok(None),
        other => fixed_bytes::<32>(other, receipt_id, field).map(Some),
    }
}

fn optional_content_param(value: Option<ContentId>) -> SqliteValue {
    value.map_or(SqliteValue::Null, |id| blob_param(id.as_bytes()))
}

fn optional_fixed_param(value: Option<[u8; 32]>) -> SqliteValue {
    value.map_or(SqliteValue::Null, |bytes| blob_param(&bytes))
}

impl Ledger {
    /// Return the typed raw-byte identity sidecar for a retained artifact.
    ///
    /// The read fails closed unless the compatibility hash, typed content ID,
    /// row version, and independently re-hashed artifact bytes all agree.
    /// The result deliberately contains no semantic identity or authority.
    ///
    /// # Errors
    /// [`LedgerError::Corrupt`] when a retained artifact lacks an exact v14
    /// sidecar or either identity disagrees with its bytes; storage failures
    /// otherwise.
    pub fn artifact_content_identity(
        &self,
        artifact_hash: &ContentHash,
    ) -> Result<Option<ArtifactContentIdentity>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT i.artifact_hash, i.content_id, i.row_schema_version \
                 FROM artifacts AS a \
                 LEFT JOIN artifact_content_identities AS i \
                   ON i.artifact_hash = a.hash \
                 WHERE a.hash = ?1 LIMIT 2",
                &[blob_param(artifact_hash.as_bytes())],
            )
            .map_err(|error| sql_err("artifact content identity read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(artifact_identity_corrupt(
                artifact_hash,
                "one artifact compatibility hash selected multiple identity sidecars",
            ));
        }
        let row = rows.first().expect("non-empty row set checked above");
        let stored_hash = match row.first() {
            Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            artifact_identity_corrupt(
                artifact_hash,
                "retained artifact has no exact 32-byte content-identity sidecar",
            )
        })?;
        if stored_hash != *artifact_hash {
            return Err(artifact_identity_corrupt(
                artifact_hash,
                "sidecar compatibility hash differs from the requested artifact hash",
            ));
        }
        let content_id = match row.get(1) {
            Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            artifact_identity_corrupt(
                artifact_hash,
                "artifact content_id is not an exact 32-byte typed content identity",
            )
        })?;
        if content_id.as_bytes() != artifact_hash.as_bytes() {
            return Err(artifact_identity_corrupt(
                artifact_hash,
                "typed content ID differs from the compatibility hash",
            ));
        }
        let row_schema_version = match row.get(2) {
            Some(SqliteValue::Integer(value)) => u32::try_from(*value).ok(),
            _ => None,
        }
        .ok_or_else(|| {
            artifact_identity_corrupt(
                artifact_hash,
                "artifact content-identity row version is not a u32 integer",
            )
        })?;
        if row_schema_version != ARTIFACT_CONTENT_IDENTITY_ROW_VERSION {
            return Err(artifact_identity_corrupt(
                artifact_hash,
                format!(
                    "artifact content-identity row version {row_schema_version} differs from supported {}",
                    ARTIFACT_CONTENT_IDENTITY_ROW_VERSION
                ),
            ));
        }
        if self
            .read_artifact_chunks(artifact_hash, &mut |_| {})?
            .is_none()
        {
            return Err(artifact_identity_corrupt(
                artifact_hash,
                "artifact disappeared while its typed content identity was being verified",
            ));
        }
        Ok(Some(ArtifactContentIdentity {
            artifact_hash: stored_hash,
            content_id,
            row_schema_version,
        }))
    }

    /// Authenticate the complete v14 artifact backfill before its schema
    /// marker commits. The migration transaction rolls every sidecar row and
    /// trigger back if any source artifact is corrupt, missing, or divergent.
    pub(crate) fn verify_artifact_content_identity_backfill(&self) -> Result<(), LedgerError> {
        let integrity = self.verify_artifact_integrity()?;
        if let Some(corrupt) = integrity.corrupted.first() {
            return Err(LedgerError::Corrupt {
                hash_hex: corrupt.clone(),
                detail: "v14 artifact identity backfill source failed independent content re-hash"
                    .to_string(),
            });
        }

        let invalid_rows = self
            .conn
            .query(
                "SELECT a.hash \
                 FROM artifacts AS a \
                 LEFT JOIN artifact_content_identities AS i \
                   ON i.artifact_hash = a.hash \
                 WHERE i.artifact_hash IS NULL \
                    OR typeof(i.artifact_hash) != 'blob' \
                    OR length(i.artifact_hash) != 32 \
                    OR i.artifact_hash != a.hash \
                    OR typeof(i.content_id) != 'blob' \
                    OR length(i.content_id) != 32 \
                    OR i.content_id != a.hash \
                    OR typeof(i.row_schema_version) != 'integer' \
                    OR i.row_schema_version != 1 \
                 LIMIT 1",
            )
            .map_err(|error| sql_err("verify artifact identity backfill", &error))?;
        if let Some(row) = invalid_rows.first() {
            let hash_hex = match row.first() {
                Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes)
                    .map_or_else(|| "<malformed>".to_string(), |hash| hash.to_hex()),
                _ => "<malformed>".to_string(),
            };
            return Err(LedgerError::Corrupt {
                hash_hex,
                detail: "v14 backfill did not produce one exact typed content identity for every artifact"
                    .to_string(),
            });
        }

        let orphan_rows = self
            .conn
            .query(
                "SELECT i.artifact_hash \
                 FROM artifact_content_identities AS i \
                 LEFT JOIN artifacts AS a ON a.hash = i.artifact_hash \
                 WHERE a.hash IS NULL LIMIT 1",
            )
            .map_err(|error| sql_err("verify artifact identity orphan", &error))?;
        if let Some(row) = orphan_rows.first() {
            let hash_hex = match row.first() {
                Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes)
                    .map_or_else(|| "<malformed>".to_string(), |hash| hash.to_hex()),
                _ => "<malformed>".to_string(),
            };
            return Err(LedgerError::Corrupt {
                hash_hex,
                detail: "v14 artifact content-identity sidecar has no retained artifact"
                    .to_string(),
            });
        }
        Ok(())
    }

    /// Return the typed artifact content identity carried by one exact lineage
    /// edge. Absence means the compatibility edge itself does not exist.
    ///
    /// The read verifies the edge sidecar, the v14 artifact sidecar, and the
    /// retained artifact bytes. It never projects a semantic identity or
    /// operation authority from edge presence.
    ///
    /// # Errors
    /// [`LedgerError::OpCorrupt`] for a missing or malformed v15 edge sidecar;
    /// [`LedgerError::Corrupt`] when the linked artifact identity or bytes are
    /// invalid; storage failures otherwise.
    pub fn edge_content_identity(
        &self,
        op: i64,
        artifact_hash: &ContentHash,
        role: EdgeRole,
    ) -> Result<Option<EdgeContentIdentity>, LedgerError> {
        let role_text = match role {
            EdgeRole::In => "in",
            EdgeRole::Out => "out",
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT i.op, i.artifact_hash, i.role, i.content_id, i.row_schema_version \
                 FROM edges AS e \
                 LEFT JOIN edge_content_identities AS i \
                   ON i.op = e.op AND i.artifact_hash = e.artifact AND i.role = e.role \
                 WHERE e.op = ?1 AND e.artifact = ?2 AND e.role = ?3 LIMIT 2",
                &[
                    SqliteValue::Integer(op),
                    blob_param(artifact_hash.as_bytes()),
                    SqliteValue::Text(role_text.into()),
                ],
            )
            .map_err(|error| sql_err("edge content identity read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(edge_identity_corrupt(
                op,
                "one role-qualified edge selected multiple typed identity sidecars",
            ));
        }
        let row = rows.first().expect("non-empty row set checked above");
        if !matches!(row.first(), Some(SqliteValue::Integer(stored)) if *stored == op) {
            return Err(edge_identity_corrupt(
                op,
                "edge identity sidecar has a missing or divergent operation identity",
            ));
        }
        let stored_hash = match row.get(1) {
            Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            edge_identity_corrupt(
                op,
                "edge identity sidecar has no exact 32-byte artifact compatibility hash",
            )
        })?;
        if stored_hash != *artifact_hash {
            return Err(edge_identity_corrupt(
                op,
                "edge identity sidecar names a different artifact compatibility hash",
            ));
        }
        if !matches!(row.get(2), Some(SqliteValue::Text(stored)) if stored.as_str() == role_text) {
            return Err(edge_identity_corrupt(
                op,
                "edge identity sidecar has a missing or divergent in/out role",
            ));
        }
        let content_id = match row.get(3) {
            Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            edge_identity_corrupt(
                op,
                "edge content_id is not an exact 32-byte typed content identity",
            )
        })?;
        if content_id.as_bytes() != artifact_hash.as_bytes() {
            return Err(edge_identity_corrupt(
                op,
                "edge typed content ID differs from its artifact compatibility hash",
            ));
        }
        let row_schema_version = match row.get(4) {
            Some(SqliteValue::Integer(value)) => u32::try_from(*value).ok(),
            _ => None,
        }
        .ok_or_else(|| {
            edge_identity_corrupt(op, "edge content-identity row version is not a u32 integer")
        })?;
        if row_schema_version != EDGE_CONTENT_IDENTITY_ROW_VERSION {
            return Err(edge_identity_corrupt(
                op,
                format!(
                    "edge content-identity row version {row_schema_version} differs from supported {}",
                    EDGE_CONTENT_IDENTITY_ROW_VERSION
                ),
            ));
        }
        let artifact_identity =
            self.artifact_content_identity(artifact_hash)?
                .ok_or_else(|| {
                    edge_identity_corrupt(op, "edge identity names an artifact that disappeared")
                })?;
        if artifact_identity.content_id() != content_id {
            return Err(edge_identity_corrupt(
                op,
                "edge and artifact typed content identities disagree",
            ));
        }
        Ok(Some(EdgeContentIdentity {
            op,
            role,
            artifact_hash: stored_hash,
            content_id,
            row_schema_version,
        }))
    }

    /// Authenticate the complete v15 edge backfill before its marker commits.
    pub(crate) fn verify_edge_content_identity_backfill(&self) -> Result<(), LedgerError> {
        self.verify_artifact_content_identity_backfill()?;
        let invalid_rows = self
            .conn
            .query(
                "SELECT e.op, e.artifact \
                 FROM edges AS e \
                 LEFT JOIN edge_content_identities AS i \
                   ON i.op = e.op AND i.artifact_hash = e.artifact AND i.role = e.role \
                 WHERE i.op IS NULL \
                    OR typeof(i.op) != 'integer' \
                    OR i.op != e.op \
                    OR typeof(i.artifact_hash) != 'blob' \
                    OR length(i.artifact_hash) != 32 \
                    OR i.artifact_hash != e.artifact \
                    OR typeof(i.role) != 'text' \
                    OR i.role != e.role \
                    OR typeof(i.content_id) != 'blob' \
                    OR length(i.content_id) != 32 \
                    OR i.content_id != e.artifact \
                    OR typeof(i.row_schema_version) != 'integer' \
                    OR i.row_schema_version != 1 \
                 LIMIT 1",
            )
            .map_err(|error| sql_err("verify edge identity backfill", &error))?;
        if let Some(row) = invalid_rows.first() {
            let op = match row.first() {
                Some(SqliteValue::Integer(op)) => *op,
                _ => -1,
            };
            return Err(edge_identity_corrupt(
                op,
                "v15 backfill did not produce one exact typed content identity for every edge",
            ));
        }

        let orphan_rows = self
            .conn
            .query(
                "SELECT i.op \
                 FROM edge_content_identities AS i \
                 LEFT JOIN edges AS e \
                   ON e.op = i.op AND e.artifact = i.artifact_hash AND e.role = i.role \
                 WHERE e.op IS NULL LIMIT 1",
            )
            .map_err(|error| sql_err("verify edge identity orphan", &error))?;
        if let Some(row) = orphan_rows.first() {
            let op = match row.first() {
                Some(SqliteValue::Integer(op)) => *op,
                _ => -1,
            };
            return Err(edge_identity_corrupt(
                op,
                "v15 edge content-identity sidecar has no compatibility edge",
            ));
        }
        Ok(())
    }

    fn write_op_content_identity(
        &self,
        op: i64,
        session: Option<&[u8]>,
        ir: &str,
        explicits: &FiveExplicits<'_>,
        if_absent: bool,
    ) -> Result<OpContentIdentity, LedgerError> {
        let identity = derive_op_content_identity(op, session, ir, explicits);
        let sql = if if_absent {
            "INSERT OR IGNORE INTO op_content_identities(
                 op, session_content_id, ir_content_id, seed_content_id,
                 versions_content_id, budget_content_id, capability_content_id,
                 row_schema_version
             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)"
        } else {
            "INSERT INTO op_content_identities(
                 op, session_content_id, ir_content_id, seed_content_id,
                 versions_content_id, budget_content_id, capability_content_id,
                 row_schema_version
             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)"
        };
        let affected = self
            .conn
            .prepare(sql)
            .map_err(|error| sql_err("operation content identity prepare", &error))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                optional_content_param(identity.session_content_id),
                blob_param(identity.ir_content_id.as_bytes()),
                blob_param(identity.seed_content_id.as_bytes()),
                blob_param(identity.versions_content_id.as_bytes()),
                blob_param(identity.budget_content_id.as_bytes()),
                blob_param(identity.capability_content_id.as_bytes()),
                SqliteValue::Integer(i64::from(identity.row_schema_version)),
            ])
            .map_err(|error| sql_err("operation content identity write", &error))?;
        if affected == 1 || (if_absent && affected == 0) {
            Ok(identity)
        } else {
            Err(op_identity_corrupt(
                op,
                format!("one operation identity write changed {affected} rows"),
            ))
        }
    }

    /// Persist exact typed content IDs for the operation fields admitted by
    /// [`Ledger::begin_op_on`](crate::Ledger::begin_op_on). The caller keeps
    /// this insert in the same transaction as the compatibility operation row.
    pub(crate) fn insert_op_content_identity(
        &self,
        op: i64,
        session: Option<&[u8]>,
        ir: &str,
        explicits: &FiveExplicits<'_>,
    ) -> Result<(), LedgerError> {
        self.write_op_content_identity(op, session, ir, explicits, false)
            .map(|_| ())
    }

    /// Return independently re-hashed typed identities for one operation's
    /// exact frozen session, IR, seed, versions, budget, and capability bytes.
    ///
    /// Row IDs, branches, execution mode, clocks, outcomes, diagnostics,
    /// lineage, semantic meaning, and authority remain separate. Absence means
    /// the compatibility operation itself does not exist.
    ///
    /// # Errors
    /// [`LedgerError::OpCorrupt`] when the operation envelope or its v18
    /// sidecar is missing, malformed, future-versioned, or content-divergent;
    /// storage failures otherwise.
    pub fn op_content_identity(&self, op: i64) -> Result<Option<OpContentIdentity>, LedgerError> {
        let Some(source) = self.op(op)? else {
            return Ok(None);
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT op, session_content_id, ir_content_id, seed_content_id,
                        versions_content_id, budget_content_id, capability_content_id,
                        row_schema_version
                 FROM op_content_identities WHERE op = ?1 LIMIT 2",
                &[SqliteValue::Integer(op)],
            )
            .map_err(|error| sql_err("operation content identity read", &error))?;
        if rows.len() != 1 {
            return Err(op_identity_corrupt(
                op,
                if rows.is_empty() {
                    "retained operation has no typed content-identity sidecar"
                } else {
                    "one operation row selected multiple typed content-identity sidecars"
                },
            ));
        }
        let row = rows
            .first()
            .expect("single operation identity row checked above");
        if !matches!(row.first(), Some(SqliteValue::Integer(stored)) if *stored == op) {
            return Err(op_identity_corrupt(
                op,
                "operation identity sidecar has a missing or divergent row ID",
            ));
        }
        let stored = OpContentIdentity {
            op,
            session_content_id: optional_op_content_id(row.get(1), op, "session_content_id")?,
            ir_content_id: op_content_id(row.get(2), op, "ir_content_id")?,
            seed_content_id: op_content_id(row.get(3), op, "seed_content_id")?,
            versions_content_id: op_content_id(row.get(4), op, "versions_content_id")?,
            budget_content_id: op_content_id(row.get(5), op, "budget_content_id")?,
            capability_content_id: op_content_id(row.get(6), op, "capability_content_id")?,
            row_schema_version: match row.get(7) {
                Some(SqliteValue::Integer(value)) => u32::try_from(*value).map_err(|_| {
                    op_identity_corrupt(
                        op,
                        "operation content-identity row version is outside the u32 domain",
                    )
                })?,
                _ => {
                    return Err(op_identity_corrupt(
                        op,
                        "operation content-identity row version is not an INTEGER",
                    ));
                }
            },
        };
        if stored.row_schema_version != OP_CONTENT_IDENTITY_ROW_VERSION {
            return Err(op_identity_corrupt(
                op,
                format!(
                    "operation content-identity row version {} differs from supported {}",
                    stored.row_schema_version, OP_CONTENT_IDENTITY_ROW_VERSION
                ),
            ));
        }
        let explicits = FiveExplicits {
            seed: &source.seed,
            versions: &source.versions,
            budget: &source.budget,
            capability: &source.capability,
        };
        let expected =
            derive_op_content_identity(op, source.session.as_deref(), &source.ir, &explicits);
        for (field, found, required) in [
            (
                "ir_content_id",
                stored.ir_content_id,
                expected.ir_content_id,
            ),
            (
                "seed_content_id",
                stored.seed_content_id,
                expected.seed_content_id,
            ),
            (
                "versions_content_id",
                stored.versions_content_id,
                expected.versions_content_id,
            ),
            (
                "budget_content_id",
                stored.budget_content_id,
                expected.budget_content_id,
            ),
            (
                "capability_content_id",
                stored.capability_content_id,
                expected.capability_content_id,
            ),
        ] {
            if found != required {
                return Err(op_identity_corrupt(
                    op,
                    format!("{field} differs from the independently re-hashed source bytes"),
                ));
            }
        }
        if stored.session_content_id != expected.session_content_id {
            return Err(op_identity_corrupt(
                op,
                "session_content_id differs from the optional exact session bytes",
            ));
        }
        Ok(Some(stored))
    }

    /// Reconcile one operation inserted by a compatible pre-v18 writer with
    /// its exact typed raw-content sidecar.
    ///
    /// The operation is re-read through the bounded storage envelope and each
    /// frozen field is independently hashed. A missing sidecar is inserted in
    /// the same transaction; an existing sidecar must already agree exactly
    /// and is never rewritten. This assigns no IR schema, semantic identity,
    /// or authority. Absence means the compatibility operation does not exist.
    ///
    /// When the caller has an open transaction, reconciliation participates in
    /// it. Otherwise this method owns one transaction and rolls it back on any
    /// refusal, so a partial sidecar is never published.
    ///
    /// # Errors
    /// [`LedgerError::OpCorrupt`] when the operation is malformed or an
    /// existing sidecar is partial, future-versioned, or divergent; storage
    /// failures otherwise.
    pub fn reconcile_op_content_identity(
        &self,
        op: i64,
    ) -> Result<Option<OpContentIdentity>, LedgerError> {
        let owns_transaction = !self.in_transaction();
        if owns_transaction {
            self.begin()?;
        }
        let reconciled = (|| {
            let Some(source) = self.op(op)? else {
                return Ok(None);
            };
            let explicits = FiveExplicits {
                seed: &source.seed,
                versions: &source.versions,
                budget: &source.budget,
                capability: &source.capability,
            };
            self.write_op_content_identity(
                op,
                source.session.as_deref(),
                &source.ir,
                &explicits,
                true,
            )?;
            self.op_content_identity(op)
        })();
        match (&reconciled, owns_transaction) {
            (Ok(_), true) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
            }
            (Err(_), true) => {
                let _ = self.rollback();
            }
            _ => {}
        }
        reconciled
    }

    /// Backfill and authenticate v18 operation content identities before the
    /// schema marker commits. Work is paged by 64 fixed-size row IDs and only
    /// one bounded operation payload is retained at a time.
    pub(crate) fn backfill_and_verify_op_content_identities(&self) -> Result<(), LedgerError> {
        let mut after = None;
        loop {
            let rows = match after {
                Some(last) => self.conn.query_with_params(
                    "SELECT id FROM ops WHERE id > ?1 ORDER BY id LIMIT 64",
                    &[SqliteValue::Integer(last)],
                ),
                None => self.conn.query("SELECT id FROM ops ORDER BY id LIMIT 64"),
            }
            .map_err(|error| sql_err("operation identity backfill page", &error))?;
            if rows.is_empty() {
                break;
            }
            for row in &rows {
                let op = match row.first() {
                    Some(SqliteValue::Integer(op)) => *op,
                    _ => {
                        return Err(op_identity_corrupt(
                            -1,
                            "operation backfill selected a non-integer row ID",
                        ));
                    }
                };
                let source = self.op(op)?.ok_or_else(|| {
                    op_identity_corrupt(op, "operation disappeared during v18 backfill")
                })?;
                let explicits = FiveExplicits {
                    seed: &source.seed,
                    versions: &source.versions,
                    budget: &source.budget,
                    capability: &source.capability,
                };
                self.write_op_content_identity(
                    op,
                    source.session.as_deref(),
                    &source.ir,
                    &explicits,
                    true,
                )?;
                let _ = self.op_content_identity(op)?.ok_or_else(|| {
                    op_identity_corrupt(op, "v18 backfill lost its compatibility operation")
                })?;
                after = Some(op);
            }
            if rows.len() < 64 {
                break;
            }
        }

        let orphans = self
            .conn
            .query(
                "SELECT i.op FROM op_content_identities AS i
                 LEFT JOIN ops AS o ON o.id = i.op
                 WHERE o.id IS NULL LIMIT 1",
            )
            .map_err(|error| sql_err("verify operation identity orphan", &error))?;
        if let Some(row) = orphans.first() {
            let op = match row.first() {
                Some(SqliteValue::Integer(op)) => *op,
                _ => -1,
            };
            return Err(op_identity_corrupt(
                op,
                "v18 operation content-identity sidecar has no compatibility operation",
            ));
        }
        Ok(())
    }

    fn bounded_tune_key_at_rowid(
        &self,
        rowid: i64,
    ) -> Result<(String, String, Vec<u8>), LedgerError> {
        let sql = format!(
            "SELECT kernel, shape_class, machine FROM tune
             WHERE rowid = ?1 AND
                   typeof(kernel) = 'text' AND
                   length(CAST(kernel AS BLOB)) BETWEEN 1 AND {MAX_TUNE_KERNEL_BYTES} AND
                   length(CAST(kernel AS BLOB)) = length(kernel) AND
                   kernel NOT GLOB '*[^!-~]*' AND
                   typeof(shape_class) = 'text' AND
                   length(CAST(shape_class AS BLOB)) BETWEEN 1 AND {MAX_TUNE_SHAPE_CLASS_BYTES} AND
                   length(CAST(shape_class AS BLOB)) = length(shape_class) AND
                   shape_class NOT GLOB '*[^!-~]*' AND
                   typeof(machine) = 'blob' AND
                   length(machine) BETWEEN 1 AND {MAX_TUNE_MACHINE_BYTES} AND
                   CASE WHEN typeof(params) = 'text' THEN
                       CASE WHEN length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES}
                           THEN json_valid(params) ELSE 0 END
                       ELSE 0 END = 1 AND
                   CASE WHEN typeof(measured) = 'text' THEN
                       CASE WHEN length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES}
                           THEN json_valid(measured) ELSE 0 END
                       ELSE 0 END = 1 LIMIT 2"
        );
        let rows = self
            .conn
            .query_with_params(&sql, &[SqliteValue::Integer(rowid)])
            .map_err(|error| sql_err("tune identity bounded-key read", &error))?;
        let row_label = format!("<rowid:{rowid}>");
        if rows.len() != 1 {
            let exists = self
                .conn
                .query_with_params(
                    "SELECT 1 FROM tune WHERE rowid = ?1 LIMIT 1",
                    &[SqliteValue::Integer(rowid)],
                )
                .map_err(|error| sql_err("tune identity rowid existence", &error))?;
            return Err(tune_corrupt(
                &row_label,
                if exists.is_empty() {
                    "cache row disappeared during the v19 migration transaction"
                } else {
                    "cache row is outside the bounded canonical tune storage contract"
                },
            ));
        }
        let row = rows
            .first()
            .expect("single bounded tune key row checked above");
        let kernel = row_text(row, 0, "tune identity kernel")?;
        let shape_class = row_text(row, 1, "tune identity shape_class")?;
        let machine = match row.get(2) {
            Some(SqliteValue::Blob(bytes)) => bytes.to_vec(),
            other => {
                return Err(tune_corrupt(
                    &kernel,
                    format!("machine is not a bounded BLOB: {other:?}"),
                ));
            }
        };
        Ok((kernel, shape_class, machine))
    }

    fn write_tune_content_identity(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
        if_absent: bool,
    ) -> Result<TuneContentIdentity, LedgerError> {
        let identity = derive_tune_content_identity(kernel, shape_class, machine, params, measured);
        let sql = if if_absent {
            "INSERT OR IGNORE INTO tune_content_identities(
                 kernel, shape_class, machine, kernel_content_id,
                 shape_class_content_id, machine_content_id, params_content_id,
                 measured_content_id, row_schema_version
             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)"
        } else {
            "INSERT INTO tune_content_identities(
                 kernel, shape_class, machine, kernel_content_id,
                 shape_class_content_id, machine_content_id, params_content_id,
                 measured_content_id, row_schema_version
             ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
             ON CONFLICT(kernel, shape_class, machine) DO UPDATE SET
                 params_content_id = excluded.params_content_id,
                 measured_content_id = excluded.measured_content_id"
        };
        let affected = self
            .conn
            .prepare(sql)
            .map_err(|error| sql_err("tune content identity prepare", &error))?
            .execute_with_params(&[
                text_param(kernel),
                text_param(shape_class),
                blob_param(machine),
                blob_param(identity.kernel_content_id.as_bytes()),
                blob_param(identity.shape_class_content_id.as_bytes()),
                blob_param(identity.machine_content_id.as_bytes()),
                blob_param(identity.params_content_id.as_bytes()),
                blob_param(identity.measured_content_id.as_bytes()),
                SqliteValue::Integer(i64::from(identity.row_schema_version)),
            ])
            .map_err(|error| sql_err("tune content identity write", &error))?;
        if affected != 1 && !(if_absent && affected == 0) {
            return Err(tune_corrupt(
                kernel,
                format!("one tune identity write changed {affected} rows"),
            ));
        }
        let stored = self
            .tune_content_identity_inner(kernel, shape_class, machine)?
            .ok_or_else(|| tune_corrupt(kernel, "cache row disappeared after sidecar write"))?;
        if stored != identity {
            return Err(tune_corrupt(
                kernel,
                "stored tune content identity differs from the requested exact row bytes",
            ));
        }
        Ok(stored)
    }

    pub(crate) fn upsert_tune_content_identity(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        self.write_tune_content_identity(kernel, shape_class, machine, params, measured, false)
            .map(|_| ())
    }

    /// Return independently re-hashed typed identities for one autotuner
    /// cache row's exact kernel, shape, machine, params, and measured bytes.
    ///
    /// Cache-key semantics, JSON schemas, scientific validity, freshness, and
    /// authority remain separate. Absence means the compatibility cache row
    /// itself does not exist.
    ///
    /// # Errors
    /// [`LedgerError::TuneCorrupt`] when a retained row or its v19 sidecar is
    /// missing, malformed, future-versioned, or content-divergent; invalid
    /// lookup keys and storage failures otherwise.
    pub fn tune_content_identity(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneContentIdentity>, LedgerError> {
        self.note_read_query();
        self.tune_content_identity_inner(kernel, shape_class, machine)
    }

    /// Reconcile one autotuner row written by a compatible pre-v19 writer with
    /// its exact typed raw-content sidecar.
    ///
    /// The complete bounded source row is re-read and independently hashed. A
    /// missing sidecar is inserted; for an existing sidecar only the mutable
    /// params/measured content IDs may move. The immutable kernel, shape,
    /// machine, and row-schema projection must still verify exactly. This does
    /// not assign cache-key semantics, freshness, scientific validity, or
    /// authority. Absence means the compatibility cache row does not exist.
    ///
    /// When the caller has an open transaction, reconciliation participates in
    /// it. Otherwise this method owns one transaction and rolls it back on any
    /// refusal, so source and sidecar cannot be partially published by this
    /// path.
    ///
    /// # Errors
    /// [`LedgerError::TuneCorrupt`] when the source row or immutable sidecar
    /// projection is malformed, future-versioned, or divergent; invalid lookup
    /// keys and storage failures otherwise.
    pub fn reconcile_tune_content_identity(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneContentIdentity>, LedgerError> {
        let owns_transaction = !self.in_transaction();
        if owns_transaction {
            self.begin()?;
        }
        let reconciled = (|| {
            let Some(source) = self.tune_get_inner(kernel, shape_class, machine)? else {
                return Ok(None);
            };
            self.write_tune_content_identity(
                &source.kernel,
                &source.shape_class,
                &source.machine,
                &source.params,
                &source.measured,
                false,
            )
            .map(Some)
        })();
        match (&reconciled, owns_transaction) {
            (Ok(_), true) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
            }
            (Err(_), true) => {
                let _ = self.rollback();
            }
            _ => {}
        }
        reconciled
    }

    fn tune_content_identity_inner(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneContentIdentity>, LedgerError> {
        let Some(source) = self.tune_get_inner(kernel, shape_class, machine)? else {
            return Ok(None);
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT kernel_content_id, shape_class_content_id,
                        machine_content_id, params_content_id,
                        measured_content_id, row_schema_version
                 FROM tune_content_identities
                 WHERE kernel = ?1 AND shape_class = ?2 AND machine = ?3 LIMIT 2",
                &[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                ],
            )
            .map_err(|error| sql_err("tune content identity read", &error))?;
        if rows.len() != 1 {
            return Err(tune_corrupt(
                kernel,
                if rows.is_empty() {
                    "retained cache row has no typed content-identity sidecar"
                } else {
                    "one cache key selected multiple typed content-identity sidecars"
                },
            ));
        }
        let row = rows
            .first()
            .expect("single tune content identity row checked above");
        let stored = TuneContentIdentity {
            kernel_content_id: tune_content_id(row.first(), kernel, "kernel_content_id")?,
            shape_class_content_id: tune_content_id(row.get(1), kernel, "shape_class_content_id")?,
            machine_content_id: tune_content_id(row.get(2), kernel, "machine_content_id")?,
            params_content_id: tune_content_id(row.get(3), kernel, "params_content_id")?,
            measured_content_id: tune_content_id(row.get(4), kernel, "measured_content_id")?,
            row_schema_version: match row.get(5) {
                Some(SqliteValue::Integer(value)) => u32::try_from(*value).map_err(|_| {
                    tune_corrupt(
                        kernel,
                        "tune content-identity row version is outside the u32 domain",
                    )
                })?,
                _ => {
                    return Err(tune_corrupt(
                        kernel,
                        "tune content-identity row version is not an INTEGER",
                    ));
                }
            },
        };
        if stored.row_schema_version != TUNE_CONTENT_IDENTITY_ROW_VERSION {
            return Err(tune_corrupt(
                kernel,
                format!(
                    "tune content-identity row version {} differs from supported {}",
                    stored.row_schema_version, TUNE_CONTENT_IDENTITY_ROW_VERSION
                ),
            ));
        }
        let expected = derive_tune_content_identity(
            &source.kernel,
            &source.shape_class,
            &source.machine,
            &source.params,
            &source.measured,
        );
        for (field, found, required) in [
            (
                "kernel_content_id",
                stored.kernel_content_id,
                expected.kernel_content_id,
            ),
            (
                "shape_class_content_id",
                stored.shape_class_content_id,
                expected.shape_class_content_id,
            ),
            (
                "machine_content_id",
                stored.machine_content_id,
                expected.machine_content_id,
            ),
            (
                "params_content_id",
                stored.params_content_id,
                expected.params_content_id,
            ),
            (
                "measured_content_id",
                stored.measured_content_id,
                expected.measured_content_id,
            ),
        ] {
            if found != required {
                return Err(tune_corrupt(
                    kernel,
                    format!("{field} differs from the independently re-hashed source bytes"),
                ));
            }
        }
        Ok(Some(stored))
    }

    pub(crate) fn verify_tune_content_identity(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<(), LedgerError> {
        self.tune_content_identity_inner(kernel, shape_class, machine)?
            .map(|_| ())
            .ok_or_else(|| tune_corrupt(kernel, "cache row has no compatibility source"))
    }

    /// Backfill and authenticate v19 autotuner content identities before the
    /// schema marker commits. Only fixed-size row IDs are paged; each bounded
    /// cache row is materialized and re-hashed independently.
    pub(crate) fn backfill_and_verify_tune_content_identities(&self) -> Result<(), LedgerError> {
        let mut after = None;
        loop {
            let rows = match after {
                Some(last) => self.conn.query_with_params(
                    "SELECT rowid FROM tune WHERE rowid > ?1 ORDER BY rowid LIMIT 64",
                    &[SqliteValue::Integer(last)],
                ),
                None => self
                    .conn
                    .query("SELECT rowid FROM tune ORDER BY rowid LIMIT 64"),
            }
            .map_err(|error| sql_err("tune identity backfill page", &error))?;
            if rows.is_empty() {
                break;
            }
            for row in &rows {
                let rowid = match row.first() {
                    Some(SqliteValue::Integer(rowid)) => *rowid,
                    _ => {
                        return Err(tune_corrupt(
                            "<migration>",
                            "tune identity backfill selected a non-integer row ID",
                        ));
                    }
                };
                let (kernel, shape_class, machine) = self.bounded_tune_key_at_rowid(rowid)?;
                let source = self
                    .tune_get_inner(&kernel, &shape_class, &machine)?
                    .ok_or_else(|| {
                        tune_corrupt(&kernel, "cache row disappeared during v19 backfill")
                    })?;
                self.write_tune_content_identity(
                    &source.kernel,
                    &source.shape_class,
                    &source.machine,
                    &source.params,
                    &source.measured,
                    true,
                )?;
                after = Some(rowid);
            }
            if rows.len() < 64 {
                break;
            }
        }

        let orphans = self
            .conn
            .query(
                "SELECT i.kernel_content_id FROM tune_content_identities AS i
                 LEFT JOIN tune AS t
                   ON t.kernel = i.kernel
                  AND t.shape_class = i.shape_class
                  AND t.machine = i.machine
                 WHERE t.kernel IS NULL LIMIT 1",
            )
            .map_err(|error| sql_err("verify tune identity orphan", &error))?;
        if let Some(row) = orphans.first() {
            let kernel = match row.first() {
                Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes)
                    .map(|id| format!("<content:{}>", id.to_hex()))
                    .unwrap_or_else(|| "<malformed-content-id>".to_string()),
                _ => "<malformed-content-id>".to_string(),
            };
            return Err(tune_corrupt(
                &kernel,
                "v19 tune content-identity sidecar has no compatibility cache row",
            ));
        }
        Ok(())
    }

    fn stored_artifact_semantic_binding(
        &self,
        artifact_hash: &ContentHash,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<ArtifactSemanticBinding>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT artifact_hash, receipt_id, semantic_id, identity_role, \
                        identity_schema_id, identity_schema_version, trust_state, no_claim_state \
                 FROM artifact_semantic_bindings \
                 WHERE artifact_hash = ?1 AND receipt_id = ?2 LIMIT 2",
                &[
                    blob_param(artifact_hash.as_bytes()),
                    blob_param(receipt_id.as_bytes()),
                ],
            )
            .map_err(|error| sql_err("artifact semantic binding read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "one artifact/receipt key selected multiple rows",
            ));
        }
        let row = rows.first().expect("non-empty row set checked above");
        let stored_artifact = match row.first() {
            Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "artifact_hash is not an exact 32-byte compatibility hash",
            )
        })?;
        if stored_artifact != *artifact_hash {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "stored artifact hash differs from the requested key",
            ));
        }
        let stored_receipt_bytes = fixed_bytes::<32>(row.get(1), receipt_id, "receipt_id")?;
        let stored_receipt = IdentityMigrationReceiptId::parse_slice(&stored_receipt_bytes)
            .ok_or_else(|| {
                artifact_semantic_binding_corrupt(
                    artifact_hash,
                    receipt_id,
                    "receipt_id is not a typed 32-byte identity",
                )
            })?;
        if stored_receipt != receipt_id {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "stored receipt ID differs from the requested key",
            ));
        }

        let receipt = self.stored_identity_migration(receipt_id)?.ok_or_else(|| {
            artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "referenced migration receipt is missing",
            )
        })?;
        let semantic_id = fixed_bytes::<32>(row.get(2), receipt_id, "semantic_id")?;
        let role_tag = integer(row.get(3), receipt_id, "identity_role")?;
        let schema_id = fixed_bytes::<32>(row.get(4), receipt_id, "identity_schema_id")?;
        let schema_version = u32_integer(row.get(5), receipt_id, "identity_schema_version")?;
        let trust_tag = integer(row.get(6), receipt_id, "trust_state")?;
        let no_claim_tag = integer(row.get(7), receipt_id, "no_claim_state")?;
        if semantic_id != receipt.semantic_id_bytes()
            || role_tag != i64::from(receipt.identity_role().tag())
            || schema_id != receipt.identity_schema_id()
            || schema_version != receipt.identity_schema_version()
            || trust_tag != i64::from(trust_state_tag(receipt.trust_state()))
            || no_claim_tag != i64::from(no_claim_state_tag(receipt.no_claim_state()))
        {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "stored semantic/schema/authority projection differs from the exact receipt",
            ));
        }
        if receipt.canonical_content_id().as_bytes() != artifact_hash.as_bytes() {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "receipt canonical content ID does not name the bound artifact",
            ));
        }
        let artifact_identity =
            self.artifact_content_identity(artifact_hash)?
                .ok_or_else(|| {
                    artifact_semantic_binding_corrupt(
                        artifact_hash,
                        receipt_id,
                        "bound artifact is not retained",
                    )
                })?;
        if artifact_identity.content_id() != receipt.canonical_content_id() {
            return Err(artifact_semantic_binding_corrupt(
                artifact_hash,
                receipt_id,
                "artifact typed content identity differs from the receipt canonical content ID",
            ));
        }
        Ok(Some(ArtifactSemanticBinding {
            artifact_hash: stored_artifact,
            receipt,
        }))
    }

    /// Bind one exact retained artifact to one independently verified semantic
    /// migration receipt. Exact retries dedupe; multiple distinct receipts are
    /// retained without ranking or replacement.
    ///
    /// # Errors
    /// Refuses caller-owned transactions, missing receipts or artifacts,
    /// canonical-content disagreement, malformed stored rows, and database
    /// failures.
    pub fn bind_artifact_semantic_identity(
        &self,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<ArtifactSemanticBindingWrite, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "artifact_semantic_binding.transaction",
                "artifact semantic binding must own its transaction",
            ));
        }
        self.begin()?;
        let result = (|| {
            let receipt = self.stored_identity_migration(receipt_id)?.ok_or_else(|| {
                LedgerError::NotFound {
                    what: format!("identity migration receipt {receipt_id}"),
                }
            })?;
            let artifact_hash = ContentHash(*receipt.canonical_content_id().as_bytes());
            let artifact_identity = self.artifact_content_identity(&artifact_hash)?.ok_or_else(
                || LedgerError::NotFound {
                    what: format!(
                        "canonical artifact {} required by identity migration receipt {receipt_id}",
                        artifact_hash.to_hex()
                    ),
                },
            )?;
            if artifact_identity.content_id() != receipt.canonical_content_id() {
                return Err(artifact_semantic_binding_corrupt(
                    &artifact_hash,
                    receipt_id,
                    "retained artifact content ID differs from receipt canonical content ID",
                ));
            }
            if self
                .stored_artifact_semantic_binding(&artifact_hash, receipt_id)?
                .is_some()
            {
                return Ok(ArtifactSemanticBindingWrite {
                    artifact_hash,
                    receipt_id,
                    deduped: true,
                });
            }
            self.conn
                .prepare(
                    "INSERT INTO artifact_semantic_bindings(\
                        artifact_hash, receipt_id, semantic_id, identity_role, \
                        identity_schema_id, identity_schema_version, trust_state, \
                        no_claim_state, created_at\
                     ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
                )
                .map_err(|error| sql_err("artifact semantic binding insert prepare", &error))?
                .execute_with_params(&[
                    blob_param(artifact_hash.as_bytes()),
                    blob_param(receipt_id.as_bytes()),
                    blob_param(&receipt.semantic_id_bytes()),
                    SqliteValue::Integer(i64::from(receipt.identity_role().tag())),
                    blob_param(&receipt.identity_schema_id()),
                    SqliteValue::Integer(i64::from(receipt.identity_schema_version())),
                    SqliteValue::Integer(i64::from(trust_state_tag(receipt.trust_state()))),
                    SqliteValue::Integer(i64::from(no_claim_state_tag(receipt.no_claim_state()))),
                    SqliteValue::Integer(now_wall_ns()),
                ])
                .map_err(|error| sql_err("artifact semantic binding insert", &error))?;
            Ok(ArtifactSemanticBindingWrite {
                artifact_hash,
                receipt_id,
                deduped: false,
            })
        })();
        match result {
            Ok(write) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(write)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    /// Read and independently reverify one exact artifact-semantic binding.
    ///
    /// # Errors
    /// Refuses malformed projections, missing retained prerequisites,
    /// content-identity disagreement, and database failures.
    pub fn artifact_semantic_binding(
        &self,
        artifact_hash: &ContentHash,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<ArtifactSemanticBinding>, LedgerError> {
        self.stored_artifact_semantic_binding(artifact_hash, receipt_id)
    }

    /// Return a bounded deterministic receipt-ID prefix for one artifact.
    /// This lookup never chooses or promotes a semantic interpretation.
    ///
    /// # Errors
    /// Refuses caps outside the public bound, malformed candidate keys, and
    /// database failures.
    pub fn artifact_semantic_binding_candidates(
        &self,
        artifact_hash: &ContentHash,
        cap: usize,
    ) -> Result<ArtifactSemanticBindingCandidates, LedgerError> {
        if cap > MAX_ARTIFACT_SEMANTIC_BINDING_CANDIDATES {
            return Err(invalid(
                "artifact_semantic_binding.cap",
                format!("candidate cap {cap} exceeds {MAX_ARTIFACT_SEMANTIC_BINDING_CANDIDATES}"),
            ));
        }
        let probe = cap.checked_add(1).ok_or_else(|| {
            invalid(
                "artifact_semantic_binding.cap",
                "candidate cap overflows the bounded probe",
            )
        })?;
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE \
                            WHEN typeof(receipt_id) = 'blob' AND length(receipt_id) = 32 \
                            THEN receipt_id ELSE NULL \
                        END \
                 FROM artifact_semantic_bindings \
                 WHERE artifact_hash = ?1 \
                 ORDER BY receipt_id LIMIT ?2",
                &[
                    blob_param(artifact_hash.as_bytes()),
                    SqliteValue::Integer(i64::try_from(probe).unwrap_or(i64::MAX)),
                ],
            )
            .map_err(|error| sql_err("artifact semantic candidate scan", &error))?;
        let truncated = rows.len() > cap;
        let mut receipt_ids = Vec::with_capacity(rows.len().min(cap));
        for row in rows.iter().take(cap) {
            let Some(SqliteValue::Blob(bytes)) = row.first() else {
                return Err(LedgerError::Corrupt {
                    hash_hex: artifact_hash.to_hex(),
                    detail: "artifact semantic candidate has malformed receipt ID".to_string(),
                });
            };
            let receipt_id = IdentityMigrationReceiptId::parse_slice(bytes).ok_or_else(|| {
                LedgerError::Corrupt {
                    hash_hex: artifact_hash.to_hex(),
                    detail: "artifact semantic candidate is not a typed 32-byte receipt ID"
                        .to_string(),
                }
            })?;
            receipt_ids.push(receipt_id);
        }
        Ok(ArtifactSemanticBindingCandidates {
            receipt_ids,
            truncated,
        })
    }

    /// Authenticate every pre-marker v16 binding plus all prior identity
    /// layers before initialization or migration commits.
    pub(crate) fn verify_artifact_semantic_bindings(&self) -> Result<(), LedgerError> {
        self.verify_edge_content_identity_backfill()?;
        let rows = self
            .conn
            .query(
                "SELECT artifact_hash, receipt_id \
                 FROM artifact_semantic_bindings ORDER BY artifact_hash, receipt_id",
            )
            .map_err(|error| sql_err("verify artifact semantic bindings", &error))?;
        for row in &rows {
            let artifact_hash = match row.first() {
                Some(SqliteValue::Blob(bytes)) => ContentHash::from_slice(bytes),
                _ => None,
            }
            .ok_or_else(|| LedgerError::Corrupt {
                hash_hex: "<malformed>".to_string(),
                detail: "artifact semantic binding has a malformed artifact hash".to_string(),
            })?;
            let receipt_bytes = match row.get(1) {
                Some(SqliteValue::Blob(bytes)) => bytes.as_slice(),
                _ => &[],
            };
            let receipt_id =
                IdentityMigrationReceiptId::parse_slice(receipt_bytes).ok_or_else(|| {
                    LedgerError::Corrupt {
                        hash_hex: artifact_hash.to_hex(),
                        detail: "artifact semantic binding has a malformed receipt ID".to_string(),
                    }
                })?;
            if self
                .stored_artifact_semantic_binding(&artifact_hash, receipt_id)?
                .is_none()
            {
                return Err(artifact_semantic_binding_corrupt(
                    &artifact_hash,
                    receipt_id,
                    "binding disappeared during migration verification",
                ));
            }
        }
        Ok(())
    }

    fn evidence_body_for_semantic_binding(
        &self,
        evidence_name: &str,
        receipt_id: IdentityMigrationReceiptId,
        binding_exists: bool,
    ) -> Result<Option<String>, LedgerError> {
        if let Err(error) = validate_evidence_semantic_binding_name(evidence_name) {
            return if binding_exists {
                Err(evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    format!("stored evidence name violates the binding envelope: {error}"),
                ))
            } else {
                Err(error)
            };
        }
        let refuse = |detail: String| {
            if binding_exists {
                evidence_semantic_binding_corrupt(evidence_name, receipt_id, detail)
            } else {
                invalid("evidence_semantic_binding.body", detail)
            }
        };
        let metadata_sql = format!(
            "SELECT typeof(body), length(CAST(body AS BLOB)), \
                    CASE WHEN typeof(body) = 'text' THEN \
                        CASE WHEN length(CAST(body AS BLOB)) BETWEEN 1 AND {MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES} \
                            THEN json_valid(body) ELSE 0 END \
                    ELSE 0 END \
             FROM evidence WHERE name = ?1 LIMIT 2"
        );
        let metadata = self
            .conn
            .query_with_params(&metadata_sql, &[SqliteValue::Text(evidence_name.into())])
            .map_err(|error| sql_err("evidence semantic binding metadata read", &error))?;
        if metadata.is_empty() {
            return Ok(None);
        }
        if metadata.len() != 1 {
            return Err(refuse(
                "one evidence name selected multiple source rows".to_string(),
            ));
        }
        let row = metadata
            .first()
            .ok_or_else(|| refuse("evidence metadata row disappeared".to_string()))?;
        if !matches!(row.first(), Some(SqliteValue::Text(kind)) if kind.as_str() == "text") {
            return Err(refuse("evidence body is not stored as TEXT".to_string()));
        }
        let body_len = match row.get(1) {
            Some(SqliteValue::Integer(len)) => usize::try_from(*len).ok(),
            _ => None,
        }
        .ok_or_else(|| refuse("evidence body length is not a nonnegative usize".to_string()))?;
        if body_len == 0 || body_len > MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES {
            return Err(refuse(format!(
                "evidence body must contain 1..={MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES} bytes for exact receipt binding, found {body_len}"
            )));
        }
        if !matches!(row.get(2), Some(SqliteValue::Integer(1))) {
            return Err(refuse(
                "evidence body is not valid JSON inside the binding envelope".to_string(),
            ));
        }

        let guarded_sql = format!(
            "SELECT body FROM evidence \
             WHERE name = ?1 \
               AND typeof(body) = 'text' \
               AND length(CAST(body AS BLOB)) BETWEEN 1 AND {MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES} \
               AND CASE WHEN typeof(body) = 'text' THEN \
                       CASE WHEN length(CAST(body AS BLOB)) BETWEEN 1 AND {MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES} \
                           THEN json_valid(body) ELSE 0 END \
                   ELSE 0 END = 1 \
             LIMIT 2"
        );
        let guarded = self
            .conn
            .query_with_params(&guarded_sql, &[SqliteValue::Text(evidence_name.into())])
            .map_err(|error| sql_err("evidence semantic binding guarded read", &error))?;
        if guarded.len() != 1 {
            return Err(LedgerError::Busy {
                context: "evidence semantic binding guarded read".to_string(),
                detail: "evidence row changed after bounded metadata preflight".to_string(),
            });
        }
        match guarded.first().and_then(|row| row.first()) {
            Some(SqliteValue::Text(body)) => Ok(Some(body.as_str().to_string())),
            _ => Err(refuse(
                "guarded evidence body is not stored as TEXT".to_string(),
            )),
        }
    }

    fn stored_evidence_semantic_binding(
        &self,
        evidence_name: &str,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<EvidenceSemanticBinding>, LedgerError> {
        validate_evidence_semantic_binding_name(evidence_name)?;
        let rows = self
            .conn
            .query_with_params(
                "SELECT evidence_name, receipt_id, content_id, semantic_id, identity_role, \
                        identity_schema_id, identity_schema_version, trust_state, no_claim_state \
                 FROM evidence_semantic_bindings \
                 WHERE evidence_name = ?1 AND receipt_id = ?2 LIMIT 2",
                &[
                    SqliteValue::Text(evidence_name.into()),
                    blob_param(receipt_id.as_bytes()),
                ],
            )
            .map_err(|error| sql_err("evidence semantic binding read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "one evidence/receipt key selected multiple rows",
            ));
        }
        let row = rows.first().ok_or_else(|| {
            evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "binding row disappeared after selection",
            )
        })?;
        let stored_name = match row.first() {
            Some(SqliteValue::Text(name)) => name.as_str(),
            _ => {
                return Err(evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    "evidence_name is not stored as TEXT",
                ));
            }
        };
        if stored_name != evidence_name {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "stored evidence name differs from the requested key",
            ));
        }
        let stored_receipt_bytes = fixed_bytes::<32>(row.get(1), receipt_id, "receipt_id")?;
        let stored_receipt = IdentityMigrationReceiptId::parse_slice(&stored_receipt_bytes)
            .ok_or_else(|| {
                evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    "receipt_id is not a typed 32-byte identity",
                )
            })?;
        if stored_receipt != receipt_id {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "stored receipt ID differs from the requested key",
            ));
        }
        let content_id = match row.get(2) {
            Some(SqliteValue::Blob(bytes)) => ContentId::parse_slice(bytes),
            _ => None,
        }
        .ok_or_else(|| {
            evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "content_id is not a typed 32-byte raw content identity",
            )
        })?;
        let receipt = self.stored_identity_migration(receipt_id)?.ok_or_else(|| {
            evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "referenced migration receipt is missing",
            )
        })?;
        let semantic_id = fixed_bytes::<32>(row.get(3), receipt_id, "semantic_id")?;
        let role_tag = integer(row.get(4), receipt_id, "identity_role")?;
        let schema_id = fixed_bytes::<32>(row.get(5), receipt_id, "identity_schema_id")?;
        let schema_version = u32_integer(row.get(6), receipt_id, "identity_schema_version")?;
        let trust_tag = integer(row.get(7), receipt_id, "trust_state")?;
        let no_claim_tag = integer(row.get(8), receipt_id, "no_claim_state")?;
        if content_id != receipt.canonical_content_id()
            || semantic_id != receipt.semantic_id_bytes()
            || role_tag != i64::from(receipt.identity_role().tag())
            || schema_id != receipt.identity_schema_id()
            || schema_version != receipt.identity_schema_version()
            || trust_tag != i64::from(trust_state_tag(receipt.trust_state()))
            || no_claim_tag != i64::from(no_claim_state_tag(receipt.no_claim_state()))
        {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "stored content/semantic/schema/authority projection differs from the exact receipt",
            ));
        }
        let body = self
            .evidence_body_for_semantic_binding(evidence_name, receipt_id, true)?
            .ok_or_else(|| {
                evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    "bound evidence row is missing",
                )
            })?;
        if body.as_bytes() != receipt.canonical_bytes() {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "retained evidence JSON bytes differ from the exact receipt canonical bytes",
            ));
        }
        if ContentId::of_bytes(body.as_bytes()) != content_id {
            return Err(evidence_semantic_binding_corrupt(
                evidence_name,
                receipt_id,
                "retained evidence JSON does not re-hash to the stored typed content ID",
            ));
        }
        Ok(Some(EvidenceSemanticBinding {
            evidence_name: evidence_name.into(),
            receipt,
        }))
    }

    /// Bind one exact retained evidence JSON body to one independently
    /// verified migration receipt. Exact retries dedupe; distinct receipts
    /// remain visible without ranking. A successful binding makes the exact
    /// evidence name/body pair immutable.
    ///
    /// # Errors
    /// Refuses caller-owned transactions, missing or oversized evidence,
    /// non-identical canonical bytes, malformed stored rows, and database
    /// failures.
    pub fn bind_evidence_semantic_identity(
        &self,
        evidence_name: &str,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<EvidenceSemanticBindingWrite, LedgerError> {
        validate_evidence_semantic_binding_name(evidence_name)?;
        if self.in_transaction() {
            return Err(invalid(
                "evidence_semantic_binding.transaction",
                "evidence semantic binding must own its transaction",
            ));
        }
        self.begin()?;
        let result = (|| {
            let receipt = self.stored_identity_migration(receipt_id)?.ok_or_else(|| {
                LedgerError::NotFound {
                    what: format!("identity migration receipt {receipt_id}"),
                }
            })?;
            let body = self
                .evidence_body_for_semantic_binding(evidence_name, receipt_id, false)?
                .ok_or_else(|| LedgerError::NotFound {
                    what: format!("evidence record {evidence_name:?}"),
                })?;
            if body.as_bytes() != receipt.canonical_bytes() {
                return Err(invalid(
                    "evidence_semantic_binding.canonical_bytes",
                    "retained evidence JSON must equal the receipt canonical bytes exactly; semantic JSON equivalence is insufficient",
                ));
            }
            let content_id = ContentId::of_bytes(body.as_bytes());
            if content_id != receipt.canonical_content_id() {
                return Err(evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    "exact evidence bytes do not re-hash to the receipt canonical content ID",
                ));
            }
            if self
                .stored_evidence_semantic_binding(evidence_name, receipt_id)?
                .is_some()
            {
                return Ok(EvidenceSemanticBindingWrite {
                    evidence_name: evidence_name.into(),
                    receipt_id,
                    content_id,
                    deduped: true,
                });
            }
            self.conn
                .prepare(
                    "INSERT INTO evidence_semantic_bindings(\
                        evidence_name, receipt_id, content_id, semantic_id, identity_role, \
                        identity_schema_id, identity_schema_version, trust_state, \
                        no_claim_state, created_at\
                     ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
                )
                .map_err(|error| sql_err("evidence semantic binding insert prepare", &error))?
                .execute_with_params(&[
                    SqliteValue::Text(evidence_name.into()),
                    blob_param(receipt_id.as_bytes()),
                    blob_param(content_id.as_bytes()),
                    blob_param(&receipt.semantic_id_bytes()),
                    SqliteValue::Integer(i64::from(receipt.identity_role().tag())),
                    blob_param(&receipt.identity_schema_id()),
                    SqliteValue::Integer(i64::from(receipt.identity_schema_version())),
                    SqliteValue::Integer(i64::from(trust_state_tag(receipt.trust_state()))),
                    SqliteValue::Integer(i64::from(no_claim_state_tag(receipt.no_claim_state()))),
                    SqliteValue::Integer(now_wall_ns()),
                ])
                .map_err(|error| sql_err("evidence semantic binding insert", &error))?;
            Ok(EvidenceSemanticBindingWrite {
                evidence_name: evidence_name.into(),
                receipt_id,
                content_id,
                deduped: false,
            })
        })();
        match result {
            Ok(write) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(write)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    /// Read and independently reverify one exact evidence-semantic binding.
    ///
    /// # Errors
    /// Refuses malformed projections, missing or changed source evidence,
    /// content-identity disagreement, and database failures.
    pub fn evidence_semantic_binding(
        &self,
        evidence_name: &str,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<EvidenceSemanticBinding>, LedgerError> {
        self.stored_evidence_semantic_binding(evidence_name, receipt_id)
    }

    /// Return a bounded deterministic receipt-ID prefix for one evidence row.
    /// This lookup never chooses or promotes a semantic interpretation.
    ///
    /// # Errors
    /// Refuses names or caps outside their public bounds, malformed candidate
    /// keys, and database failures.
    pub fn evidence_semantic_binding_candidates(
        &self,
        evidence_name: &str,
        cap: usize,
    ) -> Result<EvidenceSemanticBindingCandidates, LedgerError> {
        validate_evidence_semantic_binding_name(evidence_name)?;
        if cap > MAX_EVIDENCE_SEMANTIC_BINDING_CANDIDATES {
            return Err(invalid(
                "evidence_semantic_binding.cap",
                format!("candidate cap {cap} exceeds {MAX_EVIDENCE_SEMANTIC_BINDING_CANDIDATES}"),
            ));
        }
        let probe = cap.checked_add(1).ok_or_else(|| {
            invalid(
                "evidence_semantic_binding.cap",
                "candidate cap overflows the bounded probe",
            )
        })?;
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE \
                            WHEN typeof(receipt_id) = 'blob' AND length(receipt_id) = 32 \
                            THEN receipt_id ELSE NULL \
                        END \
                 FROM evidence_semantic_bindings \
                 WHERE evidence_name = ?1 \
                 ORDER BY receipt_id LIMIT ?2",
                &[
                    SqliteValue::Text(evidence_name.into()),
                    SqliteValue::Integer(i64::try_from(probe).unwrap_or(i64::MAX)),
                ],
            )
            .map_err(|error| sql_err("evidence semantic candidate scan", &error))?;
        let truncated = rows.len() > cap;
        let mut receipt_ids = Vec::with_capacity(rows.len().min(cap));
        for row in rows.iter().take(cap) {
            let Some(SqliteValue::Blob(bytes)) = row.first() else {
                return Err(LedgerError::Corrupt {
                    hash_hex: ContentId::of_bytes(evidence_name.as_bytes()).to_hex(),
                    detail: "evidence semantic candidate has malformed receipt ID".to_string(),
                });
            };
            let receipt_id = IdentityMigrationReceiptId::parse_slice(bytes).ok_or_else(|| {
                LedgerError::Corrupt {
                    hash_hex: ContentId::of_bytes(evidence_name.as_bytes()).to_hex(),
                    detail: "evidence semantic candidate is not a typed 32-byte receipt ID"
                        .to_string(),
                }
            })?;
            receipt_ids.push(receipt_id);
        }
        Ok(EvidenceSemanticBindingCandidates {
            receipt_ids,
            truncated,
        })
    }

    /// Authenticate every pre-marker v17 evidence binding plus all prior
    /// identity layers before initialization or migration commits.
    pub(crate) fn verify_evidence_semantic_bindings(&self) -> Result<(), LedgerError> {
        self.verify_artifact_semantic_bindings()?;
        let rows = self
            .conn
            .query(
                "SELECT evidence_name, receipt_id \
                 FROM evidence_semantic_bindings ORDER BY evidence_name, receipt_id",
            )
            .map_err(|error| sql_err("verify evidence semantic bindings", &error))?;
        for row in &rows {
            let evidence_name = match row.first() {
                Some(SqliteValue::Text(name)) => name.as_str(),
                _ => {
                    return Err(LedgerError::Corrupt {
                        hash_hex: "evidence:<malformed>".to_string(),
                        detail: "evidence semantic binding has a malformed evidence name"
                            .to_string(),
                    });
                }
            };
            let receipt_bytes = match row.get(1) {
                Some(SqliteValue::Blob(bytes)) => bytes.as_slice(),
                _ => &[],
            };
            let receipt_id =
                IdentityMigrationReceiptId::parse_slice(receipt_bytes).ok_or_else(|| {
                    LedgerError::Corrupt {
                        hash_hex: ContentId::of_bytes(evidence_name.as_bytes()).to_hex(),
                        detail: "evidence semantic binding has a malformed receipt ID".to_string(),
                    }
                })?;
            if self
                .stored_evidence_semantic_binding(evidence_name, receipt_id)?
                .is_none()
            {
                return Err(evidence_semantic_binding_corrupt(
                    evidence_name,
                    receipt_id,
                    "binding disappeared during migration verification",
                ));
            }
        }
        Ok(())
    }

    fn stored_identity_migration(
        &self,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<IdentityMigrationReceipt>, LedgerError> {
        let present = self
            .conn
            .query_with_params(
                "SELECT receipt_id FROM identity_migration_receipts \
                 WHERE receipt_id = ?1 LIMIT 2",
                &[blob_param(receipt_id.as_bytes())],
            )
            .map_err(|error| sql_err("identity migration existence read", &error))?;
        if present.is_empty() {
            return Ok(None);
        }
        if present.len() != 1 {
            return Err(stored_corrupt(
                receipt_id,
                "one receipt identity names multiple rows",
            ));
        }

        const GUARD: &str = "typeof(receipt_id) = 'blob' AND length(receipt_id) = 32 \
            AND typeof(legacy_bytes) = 'blob' AND length(legacy_bytes) <= 1048576 \
            AND typeof(legacy_content_id) = 'blob' AND length(legacy_content_id) = 32 \
            AND typeof(legacy_fnv) = 'blob' AND length(legacy_fnv) = 8 \
            AND typeof(canonical_bytes) = 'blob' AND length(canonical_bytes) <= 1048576 \
            AND typeof(canonical_content_id) = 'blob' AND length(canonical_content_id) = 32 \
            AND typeof(semantic_rule) = 'blob' AND length(semantic_rule) BETWEEN 1 AND 256 \
            AND typeof(semantic_id) = 'blob' AND length(semantic_id) = 32 \
            AND typeof(identity_role) = 'integer' AND identity_role BETWEEN 1 AND 12 \
            AND typeof(identity_domain) = 'blob' AND length(identity_domain) BETWEEN 1 AND 256 \
            AND typeof(identity_schema_name) = 'blob' \
                AND length(identity_schema_name) BETWEEN 1 AND 256 \
            AND typeof(identity_schema_id) = 'blob' AND length(identity_schema_id) = 32 \
            AND typeof(identity_schema_version) = 'integer' \
                AND identity_schema_version BETWEEN 1 AND 4294967295 \
            AND typeof(identity_context) = 'blob' \
                AND length(identity_context) BETWEEN 1 AND 4096 \
            AND typeof(canonical_preimage_id) = 'blob' AND length(canonical_preimage_id) = 32 \
            AND typeof(canonical_frame_bytes) = 'blob' AND length(canonical_frame_bytes) = 8 \
            AND typeof(field_count) = 'integer' AND field_count BETWEEN 0 AND 4294967295 \
            AND typeof(collection_items) = 'blob' AND length(collection_items) = 8 \
            AND typeof(max_canonical_bytes) = 'blob' AND length(max_canonical_bytes) = 8 \
            AND typeof(max_field_bytes) = 'blob' AND length(max_field_bytes) = 8 \
            AND typeof(max_fields) = 'integer' AND max_fields BETWEEN 1 AND 4294967295 \
            AND typeof(max_collection_items) = 'blob' AND length(max_collection_items) = 8 \
            AND typeof(cancellation_poll_bytes) = 'integer' \
                AND cancellation_poll_bytes BETWEEN 1 AND 4294967295 \
            AND typeof(trust_state) = 'integer' AND trust_state BETWEEN 0 AND 3 \
            AND (anchor_content_id IS NULL OR \
                (typeof(anchor_content_id) = 'blob' AND length(anchor_content_id) = 32)) \
            AND (verifier_id IS NULL OR \
                (typeof(verifier_id) = 'blob' AND length(verifier_id) = 32)) \
            AND (key_policy_id IS NULL OR \
                (typeof(key_policy_id) = 'blob' AND length(key_policy_id) = 32)) \
            AND typeof(no_claim_state) = 'integer' AND no_claim_state BETWEEN 0 AND 1 \
            AND typeof(created_at) = 'integer'";
        let query = format!(
            "SELECT receipt_id, legacy_bytes, legacy_content_id, legacy_fnv, \
                    canonical_bytes, canonical_content_id, semantic_rule, semantic_id, \
                    identity_role, identity_domain, identity_schema_name, identity_schema_id, \
                    identity_schema_version, identity_context, canonical_preimage_id, \
                    canonical_frame_bytes, field_count, collection_items, \
                    max_canonical_bytes, max_field_bytes, max_fields, max_collection_items, \
                    cancellation_poll_bytes, trust_state, anchor_content_id, verifier_id, \
                    key_policy_id, no_claim_state \
             FROM identity_migration_receipts \
             WHERE receipt_id = ?1 AND {GUARD} LIMIT 2"
        );
        let rows = self
            .conn
            .query_with_params(&query, &[blob_param(receipt_id.as_bytes())])
            .map_err(|error| sql_err("identity migration guarded read", &error))?;
        if rows.len() != 1 {
            return Err(stored_corrupt(
                receipt_id,
                "row violates the bounded v13 storage envelope",
            ));
        }
        let Some(row) = rows.first() else {
            return Err(stored_corrupt(
                receipt_id,
                "row disappeared after guarded selection",
            ));
        };
        let stored_id_bytes = fixed_bytes::<32>(row.first(), receipt_id, "receipt_id")?;
        let stored_id = IdentityMigrationReceiptId::parse_slice(&stored_id_bytes)
            .ok_or_else(|| stored_corrupt(receipt_id, "receipt_id is not a typed digest"))?;
        if stored_id != receipt_id {
            return Err(stored_corrupt(
                receipt_id,
                "indexed receipt ID disagrees with the selected row",
            ));
        }
        let identity_role_tag = integer(row.get(8), receipt_id, "identity_role")?;
        let trust_tag = integer(row.get(23), receipt_id, "trust_state")?;
        let no_claim_tag = integer(row.get(27), receipt_id, "no_claim_state")?;
        let body = IdentityMigrationBody {
            legacy_bytes: bounded_bytes(
                row.get(1),
                receipt_id,
                "legacy_bytes",
                MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
            )?,
            legacy_content_id: content_id(row.get(2), receipt_id, "legacy_content_id")?,
            legacy_fnv: LegacyProvenanceV1::new(u64_blob(row.get(3), receipt_id, "legacy_fnv")?),
            canonical_bytes: bounded_bytes(
                row.get(4),
                receipt_id,
                "canonical_bytes",
                MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
            )?,
            canonical_content_id: content_id(row.get(5), receipt_id, "canonical_content_id")?,
            semantic_rule: bounded_utf8(
                row.get(6),
                receipt_id,
                "semantic_rule",
                MAX_IDENTITY_MIGRATION_RULE_BYTES,
            )?,
            semantic_id: fixed_bytes::<32>(row.get(7), receipt_id, "semantic_id")?,
            identity_role: identity_role_from_tag(identity_role_tag).ok_or_else(|| {
                stored_corrupt(
                    receipt_id,
                    "identity_role is outside the closed role universe",
                )
            })?,
            identity_domain: bounded_utf8(
                row.get(9),
                receipt_id,
                "identity_domain",
                MAX_IDENTITY_MIGRATION_DOMAIN_BYTES,
            )?,
            identity_schema_name: bounded_utf8(
                row.get(10),
                receipt_id,
                "identity_schema_name",
                MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES,
            )?,
            identity_schema_id: fixed_bytes::<32>(row.get(11), receipt_id, "identity_schema_id")?,
            identity_schema_version: u32_integer(
                row.get(12),
                receipt_id,
                "identity_schema_version",
            )?,
            identity_context: bounded_utf8(
                row.get(13),
                receipt_id,
                "identity_context",
                MAX_IDENTITY_MIGRATION_CONTEXT_BYTES,
            )?,
            canonical_preimage_id: content_id(row.get(14), receipt_id, "canonical_preimage_id")?,
            canonical_frame_bytes: u64_blob(row.get(15), receipt_id, "canonical_frame_bytes")?,
            field_count: u32_integer(row.get(16), receipt_id, "field_count")?,
            collection_items: u64_blob(row.get(17), receipt_id, "collection_items")?,
            limits: CanonicalLimits::new(
                u64_blob(row.get(18), receipt_id, "max_canonical_bytes")?,
                u64_blob(row.get(19), receipt_id, "max_field_bytes")?,
                u32_integer(row.get(20), receipt_id, "max_fields")?,
                u64_blob(row.get(21), receipt_id, "max_collection_items")?,
                u32_integer(row.get(22), receipt_id, "cancellation_poll_bytes")?,
            ),
            trust_state: trust_state_from_tag(trust_tag).ok_or_else(|| {
                stored_corrupt(
                    receipt_id,
                    "trust_state is outside the closed state universe",
                )
            })?,
            anchor_content_id: optional_content_id(row.get(24), receipt_id, "anchor_content_id")?,
            verifier_id: optional_fixed_32(row.get(25), receipt_id, "verifier_id")?,
            key_policy_id: optional_fixed_32(row.get(26), receipt_id, "key_policy_id")?,
            no_claim_state: no_claim_state_from_tag(no_claim_tag).ok_or_else(|| {
                stored_corrupt(
                    receipt_id,
                    "no_claim_state is outside the closed state universe",
                )
            })?,
        };
        validate_receipt_body(&body).map_err(|detail| stored_corrupt(receipt_id, detail))?;
        let recomputed = derive_receipt_id(&body).map_err(|error| {
            stored_corrupt(
                receipt_id,
                format!("receipt identity reconstruction refused: {error}"),
            )
        })?;
        if recomputed != receipt_id {
            return Err(stored_corrupt(
                receipt_id,
                "complete stored receipt preimage does not reproduce receipt_id",
            ));
        }
        Ok(Some(IdentityMigrationReceipt { receipt_id, body }))
    }

    /// Persist one exact legacy-to-strong-ID crosswalk in an owned transaction.
    ///
    /// Exact retry is idempotent. The ledger independently hashes both byte
    /// payloads, checks that the audit record describes the offered typed
    /// receipt, validates authority/no-claim coherence, derives a typed receipt
    /// identity, and only then inserts the immutable row.
    ///
    /// # Errors
    /// Refuses caller-owned transactions, oversized or malformed inputs,
    /// receipt/audit disagreement, incoherent authority state, stored
    /// corruption, and database failures.
    pub fn record_identity_migration<I: StrongIdentity>(
        &self,
        claim: IdentityMigrationClaim<'_, I>,
    ) -> Result<IdentityMigrationWrite, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "identity_migration.transaction",
                "identity migration recording must own its transaction",
            ));
        }
        let body = receipt_body_from_claim(claim)?;
        let receipt_id = derive_receipt_id(&body).map_err(|error| {
            invalid(
                "identity_migration.receipt_id",
                format!("canonical receipt construction refused: {error}"),
            )
        })?;
        let offered = IdentityMigrationReceipt { receipt_id, body };

        self.begin()?;
        let result = (|| {
            if let Some(stored) = self.stored_identity_migration(receipt_id)? {
                if stored != offered {
                    return Err(stored_corrupt(
                        receipt_id,
                        "same receipt ID resolved to different exact fields",
                    ));
                }
                return Ok(IdentityMigrationWrite {
                    receipt_id,
                    legacy_content_id: stored.legacy_content_id(),
                    canonical_content_id: stored.canonical_content_id(),
                    deduped: true,
                });
            }

            let body = &offered.body;
            let legacy_fnv = body.legacy_fnv.value().to_le_bytes();
            let canonical_frame_bytes = body.canonical_frame_bytes.to_le_bytes();
            let collection_items = body.collection_items.to_le_bytes();
            let max_canonical_bytes = body.limits.max_canonical_bytes().to_le_bytes();
            let max_field_bytes = body.limits.max_field_bytes().to_le_bytes();
            let max_collection_items = body.limits.max_collection_items().to_le_bytes();
            let params = [
                blob_param(receipt_id.as_bytes()),
                blob_param(&body.legacy_bytes),
                blob_param(body.legacy_content_id.as_bytes()),
                blob_param(&legacy_fnv),
                blob_param(&body.canonical_bytes),
                blob_param(body.canonical_content_id.as_bytes()),
                blob_param(body.semantic_rule.as_bytes()),
                blob_param(&body.semantic_id),
                SqliteValue::Integer(i64::from(body.identity_role.tag())),
                blob_param(body.identity_domain.as_bytes()),
                blob_param(body.identity_schema_name.as_bytes()),
                blob_param(&body.identity_schema_id),
                SqliteValue::Integer(i64::from(body.identity_schema_version)),
                blob_param(body.identity_context.as_bytes()),
                blob_param(body.canonical_preimage_id.as_bytes()),
                blob_param(&canonical_frame_bytes),
                SqliteValue::Integer(i64::from(body.field_count)),
                blob_param(&collection_items),
                blob_param(&max_canonical_bytes),
                blob_param(&max_field_bytes),
                SqliteValue::Integer(i64::from(body.limits.max_fields())),
                blob_param(&max_collection_items),
                SqliteValue::Integer(i64::from(body.limits.cancellation_poll_bytes())),
                SqliteValue::Integer(i64::from(trust_state_tag(body.trust_state))),
                optional_content_param(body.anchor_content_id),
                optional_fixed_param(body.verifier_id),
                optional_fixed_param(body.key_policy_id),
                SqliteValue::Integer(i64::from(no_claim_state_tag(body.no_claim_state))),
                SqliteValue::Integer(now_wall_ns()),
            ];
            self.conn
                .prepare(
                    "INSERT INTO identity_migration_receipts(\
                        receipt_id, legacy_bytes, legacy_content_id, legacy_fnv, \
                        canonical_bytes, canonical_content_id, semantic_rule, semantic_id, \
                        identity_role, identity_domain, identity_schema_name, identity_schema_id, \
                        identity_schema_version, identity_context, canonical_preimage_id, \
                        canonical_frame_bytes, field_count, collection_items, \
                        max_canonical_bytes, max_field_bytes, max_fields, max_collection_items, \
                        cancellation_poll_bytes, trust_state, anchor_content_id, verifier_id, \
                        key_policy_id, no_claim_state, created_at\
                     ) VALUES (\
                        ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, \
                        ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, \
                        ?21, ?22, ?23, ?24, ?25, ?26, ?27, ?28, ?29\
                     )",
                )
                .map_err(|error| sql_err("identity migration insert prepare", &error))?
                .execute_with_params(&params)
                .map_err(|error| sql_err("identity migration insert", &error))?;
            Ok(IdentityMigrationWrite {
                receipt_id,
                legacy_content_id: body.legacy_content_id,
                canonical_content_id: body.canonical_content_id,
                deduped: false,
            })
        })();
        match result {
            Ok(write) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(write)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    /// Read and independently reverify one immutable migration receipt.
    ///
    /// Absence is `Ok(None)`. Presence proves exact storage and receipt-ID
    /// consistency only; callers must use `typed_semantic_id::<I>()` and apply
    /// their own authority policy before consuming the semantic identity.
    ///
    /// # Errors
    /// Refuses malformed storage, content-root drift, receipt-ID drift, or
    /// database failures.
    pub fn identity_migration_receipt(
        &self,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<IdentityMigrationReceipt>, LedgerError> {
        self.stored_identity_migration(receipt_id)
    }

    /// Return a bounded, deterministic, non-authoritative candidate list for
    /// one exact legacy content ID.
    ///
    /// A zero cap is a bounded existence probe. `truncated` means at least one
    /// additional receipt exists. No candidate is selected or promoted.
    ///
    /// # Errors
    /// Refuses caps above [`MAX_IDENTITY_MIGRATION_CANDIDATES`], malformed
    /// fixed-size receipt IDs, or database failures.
    pub fn identity_migration_candidates(
        &self,
        legacy_content_id: ContentId,
        cap: usize,
    ) -> Result<IdentityMigrationCandidates, LedgerError> {
        if cap > MAX_IDENTITY_MIGRATION_CANDIDATES {
            return Err(invalid(
                "identity_migration.cap",
                format!("candidate cap {cap} exceeds {MAX_IDENTITY_MIGRATION_CANDIDATES}"),
            ));
        }
        let probe = cap.checked_add(1).ok_or_else(|| {
            invalid(
                "identity_migration.cap",
                "candidate cap overflows the bounded probe",
            )
        })?;
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE \
                            WHEN typeof(receipt_id) = 'blob' AND length(receipt_id) = 32 \
                            THEN receipt_id ELSE NULL \
                        END \
                 FROM identity_migration_receipts INDEXED BY idx_identity_migration_legacy \
                 WHERE legacy_content_id = ?1 \
                 ORDER BY receipt_id LIMIT ?2",
                &[
                    blob_param(legacy_content_id.as_bytes()),
                    SqliteValue::Integer(i64::try_from(probe).unwrap_or(i64::MAX)),
                ],
            )
            .map_err(|error| sql_err("identity migration candidate scan", &error))?;
        let truncated = rows.len() > cap;
        let mut receipt_ids = Vec::with_capacity(rows.len().min(cap));
        for row in rows.iter().take(cap) {
            let Some(SqliteValue::Blob(bytes)) = row.first() else {
                return Err(LedgerError::Corrupt {
                    hash_hex: legacy_content_id.to_hex(),
                    detail: "identity migration candidate has malformed receipt ID".to_string(),
                });
            };
            let receipt_id = IdentityMigrationReceiptId::parse_slice(bytes).ok_or_else(|| {
                LedgerError::Corrupt {
                    hash_hex: legacy_content_id.to_hex(),
                    detail: "identity migration candidate is not a typed 32-byte receipt ID"
                        .to_string(),
                }
            })?;
            receipt_ids.push(receipt_id);
        }
        Ok(IdentityMigrationCandidates {
            receipt_ids,
            truncated,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SCHEMA_VERSION, schema};

    fn drop_v13_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_identity_migration_receipts_immutable_reinsert",
            "DROP TRIGGER IF EXISTS trg_identity_migration_receipts_immutable_delete",
            "DROP TRIGGER IF EXISTS trg_identity_migration_receipts_immutable_update",
            "DROP INDEX IF EXISTS idx_identity_migration_semantic",
            "DROP INDEX IF EXISTS idx_identity_migration_canonical",
            "DROP INDEX IF EXISTS idx_identity_migration_legacy",
            "DROP TABLE IF EXISTS identity_migration_receipts",
        ] {
            ledger.conn.execute(ddl).expect("remove v13 fixture object");
        }
    }

    fn drop_v14_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_artifact_content_identity_guard_delete",
            "DROP TRIGGER IF EXISTS trg_artifact_content_identity_immutable_update",
            "DROP TRIGGER IF EXISTS trg_artifact_content_identity_dual_write",
            "DROP INDEX IF EXISTS idx_artifact_content_identity_content",
            "DROP TABLE IF EXISTS artifact_content_identities",
        ] {
            ledger.conn.execute(ddl).expect("remove v14 fixture object");
        }
    }

    fn drop_v15_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_edge_content_identity_guard_delete",
            "DROP TRIGGER IF EXISTS trg_edge_content_identity_immutable_update",
            "DROP TRIGGER IF EXISTS trg_edge_content_identity_dual_write",
            "DROP INDEX IF EXISTS idx_edge_content_identity_content",
            "DROP TABLE IF EXISTS edge_content_identities",
        ] {
            ledger.conn.execute(ddl).expect("remove v15 fixture object");
        }
    }

    fn drop_v16_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_artifact_semantic_binding_immutable_reinsert",
            "DROP TRIGGER IF EXISTS trg_artifact_semantic_binding_immutable_delete",
            "DROP TRIGGER IF EXISTS trg_artifact_semantic_binding_immutable_update",
            "DROP INDEX IF EXISTS idx_artifact_semantic_binding_semantic",
            "DROP INDEX IF EXISTS idx_artifact_semantic_binding_receipt",
            "DROP TABLE IF EXISTS artifact_semantic_bindings",
        ] {
            ledger.conn.execute(ddl).expect("remove v16 fixture object");
        }
    }

    fn drop_v17_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_evidence_semantic_binding_guard_source_delete",
            "DROP TRIGGER IF EXISTS trg_evidence_semantic_binding_guard_source_update",
            "DROP TRIGGER IF EXISTS trg_evidence_semantic_binding_immutable_reinsert",
            "DROP TRIGGER IF EXISTS trg_evidence_semantic_binding_immutable_delete",
            "DROP TRIGGER IF EXISTS trg_evidence_semantic_binding_immutable_update",
            "DROP INDEX IF EXISTS idx_evidence_semantic_binding_semantic",
            "DROP INDEX IF EXISTS idx_evidence_semantic_binding_receipt",
            "DROP INDEX IF EXISTS idx_evidence_semantic_binding_content",
            "DROP TABLE IF EXISTS evidence_semantic_bindings",
        ] {
            ledger.conn.execute(ddl).expect("remove v17 fixture object");
        }
    }

    fn drop_v18_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_op_content_identity_guard_delete",
            "DROP TRIGGER IF EXISTS trg_op_content_identity_immutable_update",
            "DROP INDEX IF EXISTS idx_op_content_identity_ir",
            "DROP INDEX IF EXISTS idx_op_content_identity_session",
            "DROP TABLE IF EXISTS op_content_identities",
        ] {
            ledger.conn.execute(ddl).expect("remove v18 fixture object");
        }
    }

    fn drop_v19_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_tune_content_identity_guard_delete",
            "DROP TRIGGER IF EXISTS trg_tune_content_identity_key_immutable",
            "DROP INDEX IF EXISTS idx_tune_content_identity_measured",
            "DROP INDEX IF EXISTS idx_tune_content_identity_params",
            "DROP INDEX IF EXISTS idx_tune_content_identity_key",
            "DROP TABLE IF EXISTS tune_content_identities",
        ] {
            ledger.conn.execute(ddl).expect("remove v19 fixture object");
        }
    }

    fn v14_edge_fixture(ledger: &Ledger) -> (i64, ContentHash) {
        let artifact = ledger
            .put_artifact("v14-edge-artifact", b"exact v14 edge bytes", None)
            .expect("store v14 edge artifact");
        let explicits = crate::FiveExplicits {
            seed: b"v14-edge-seed",
            versions: "{}",
            budget: "{}",
            capability: "{}",
        };
        let op = ledger
            .begin_op(None, "{}", &explicits, 1)
            .expect("begin v14 edge operation");
        ledger
            .link(op, &artifact.hash, EdgeRole::Out)
            .expect("store pre-v15 edge");
        (op, artifact.hash)
    }

    #[test]
    fn genuine_v12_and_stale_v13_markers_migrate_through_v19() {
        let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
        drop_v13_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 12")
            .expect("mark genuine v12 fixture");
        ledger
            .migrate_from_observed_version(12)
            .expect("migrate genuine v12 fixture");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(
            ledger.table_count("identity_migration_receipts").unwrap(),
            0
        );

        ledger
            .conn
            .execute("PRAGMA user_version = 12")
            .expect("install stale v12 marker over exact v13 objects");
        ledger
            .migrate_from_observed_version(12)
            .expect("heal exact pre-applied v13 objects");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(
            ledger.table_count("identity_migration_receipts").unwrap(),
            0
        );
    }

    #[test]
    fn divergent_early_v13_object_refuses_before_marker_advances() {
        let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
        drop_v13_objects(&ledger);
        ledger
            .conn
            .execute("CREATE TABLE identity_migration_receipts(alien INTEGER) STRICT")
            .expect("install divergent early object");
        ledger
            .conn
            .execute("PRAGMA user_version = 12")
            .expect("mark v12 fixture");
        assert!(matches!(
            ledger.migrate_from_observed_version(12),
            Err(LedgerError::SchemaMismatch {
                claimed_version: 12,
                ..
            })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 12);
    }

    #[test]
    fn v14_backfills_v13_artifacts_and_replays_a_stale_marker() {
        for preapply_v14 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
            drop_v14_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 13")
                .expect("mark genuine v13 fixture");
            let artifact = ledger
                .put_artifact("v13-artifact", b"exact v13 artifact bytes", None)
                .expect("store pre-v14 artifact");
            if preapply_v14 {
                for ddl in schema::V14 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v14 migration batch");
                }
            }

            ledger
                .migrate_from_observed_version(13)
                .expect("authenticate and backfill exact v13 artifact");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            let identity = ledger
                .artifact_content_identity(&artifact.hash)
                .unwrap()
                .expect("backfilled identity exists");
            assert_eq!(identity.artifact_hash(), artifact.hash);
            assert_eq!(
                identity.content_id(),
                ContentId::of_bytes(b"exact v13 artifact bytes")
            );
            assert_eq!(
                ledger.table_count("artifact_content_identities").unwrap(),
                1
            );
        }
    }

    #[test]
    fn v14_corrupt_source_rolls_back_rows_objects_and_marker() {
        let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
        drop_v14_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 13")
            .expect("mark genuine v13 fixture");
        let artifact = ledger
            .put_artifact("corrupt-v13-artifact", b"pre-migration bytes", None)
            .expect("store pre-v14 artifact");
        ledger
            .corrupt_artifact_for_test(&artifact.hash)
            .expect("inject source corruption");

        assert!(matches!(
            ledger.migrate_from_observed_version(13),
            Err(LedgerError::Corrupt { .. })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 13);
        let objects = ledger
            .conn
            .query(
                "SELECT name FROM sqlite_master \
                 WHERE name = 'artifact_content_identities' LIMIT 1",
            )
            .expect("inspect rolled-back v14 schema");
        assert!(objects.is_empty());
    }

    #[test]
    fn v15_backfills_v14_edges_and_replays_a_stale_marker() {
        for preapply_v15 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
            drop_v15_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 14")
                .expect("mark genuine v14 fixture");
            let (op, artifact) = v14_edge_fixture(&ledger);
            if preapply_v15 {
                for ddl in schema::V15 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v15 migration batch");
                }
            }

            ledger
                .migrate_from_observed_version(14)
                .expect("authenticate and backfill exact v14 edge");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            let identity = ledger
                .edge_content_identity(op, &artifact, EdgeRole::Out)
                .unwrap()
                .expect("backfilled edge identity exists");
            assert_eq!(identity.op(), op);
            assert_eq!(identity.artifact_hash(), artifact);
            assert_eq!(
                identity.content_id(),
                ContentId::of_bytes(b"exact v14 edge bytes")
            );
            assert_eq!(ledger.table_count("edge_content_identities").unwrap(), 1);
        }
    }

    #[test]
    fn v15_missing_artifact_sidecar_rolls_back_objects_rows_and_marker() {
        let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
        drop_v15_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 14")
            .expect("mark genuine v14 fixture");
        let (_op, artifact) = v14_edge_fixture(&ledger);
        ledger
            .conn
            .execute("DROP TRIGGER trg_artifact_content_identity_guard_delete")
            .expect("open corruption-injection window");
        ledger
            .conn
            .prepare("DELETE FROM artifact_content_identities WHERE artifact_hash = ?1")
            .expect("prepare sidecar corruption")
            .execute_with_params(&[blob_param(artifact.as_bytes())])
            .expect("remove exact artifact sidecar");
        ledger
            .conn
            .execute(schema::V14.last().expect("v14 delete-guard DDL"))
            .expect("restore exact v14 guard");

        assert!(matches!(
            ledger.migrate_from_observed_version(14),
            Err(LedgerError::Corrupt { .. })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 14);
        let objects = ledger
            .conn
            .query(
                "SELECT name FROM sqlite_master \
                 WHERE name = 'edge_content_identities' LIMIT 1",
            )
            .expect("inspect rolled-back v15 schema");
        assert!(objects.is_empty());
    }

    #[test]
    fn v16_empty_migration_and_exact_stale_marker_replay_infer_no_meaning() {
        for preapply_v16 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
            drop_v16_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 15")
                .expect("mark genuine v15 fixture");
            if preapply_v16 {
                for ddl in schema::V16 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v16 migration batch");
                }
            }
            ledger
                .migrate_from_observed_version(15)
                .expect("migrate v15 without inferring semantic bindings");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            assert_eq!(ledger.table_count("artifact_semantic_bindings").unwrap(), 0);
        }
    }

    #[test]
    fn v17_empty_migration_and_exact_stale_marker_replay_infer_no_evidence_meaning() {
        for preapply_v17 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v17 ledger");
            drop_v17_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 16")
                .expect("mark genuine v16 fixture");
            if preapply_v17 {
                for ddl in schema::V17 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v17 migration batch");
                }
            }
            ledger
                .migrate_from_observed_version(16)
                .expect("migrate v16 without inferring evidence semantic bindings");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            assert_eq!(ledger.table_count("evidence_semantic_bindings").unwrap(), 0);
        }
    }

    #[test]
    fn pre_v18_operation_writer_can_be_reconciled_without_rewriting_existing_identity() {
        let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
        ledger
            .conn
            .prepare(
                "INSERT INTO ops(
                     session, ir, seed, versions, budget, capability, t_start, branch, exec_mode
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, 1, 'deterministic')",
            )
            .expect("prepare old-writer operation")
            .execute_with_params(&[
                blob_param(b"old-writer-session"),
                text_param(r#"{"op":"old-writer"}"#),
                blob_param(b"old-writer-seed"),
                text_param(r#"{"version":17}"#),
                text_param(r#"{"memory":1024}"#),
                text_param(r#"{"cores":1}"#),
                SqliteValue::Integer(1),
            ])
            .expect("insert compatibility operation without v18 sidecar");
        let op = ledger.conn.last_insert_rowid();
        assert!(matches!(
            ledger.op_content_identity(op),
            Err(LedgerError::OpCorrupt { op: rejected, detail })
                if rejected == op && detail.contains("no typed content-identity sidecar")
        ));

        let reconciled = ledger
            .reconcile_op_content_identity(op)
            .expect("reconcile exact bounded operation")
            .expect("compatibility operation exists");
        assert_eq!(
            reconciled.session_content_id(),
            Some(ContentId::of_bytes(b"old-writer-session"))
        );
        assert_eq!(
            reconciled.ir_content_id(),
            ContentId::of_bytes(br#"{"op":"old-writer"}"#)
        );
        assert_eq!(
            ledger
                .reconcile_op_content_identity(op)
                .expect("exact retry")
                .expect("operation remains present"),
            reconciled,
            "an exact retry must retain the immutable sidecar"
        );

        ledger
            .begin()
            .expect("caller-owned reconciliation transaction");
        ledger
            .conn
            .prepare(
                "INSERT INTO ops(
                     session, ir, seed, versions, budget, capability, t_start, branch, exec_mode
                 ) VALUES (NULL, '{}', ?1, '{}', '{}', '{}', 2, 1, 'deterministic')",
            )
            .expect("prepare transactional old-writer operation")
            .execute_with_params(&[blob_param(b"transactional-old-writer-seed")])
            .expect("insert transactional compatibility operation");
        let rolled_back = ledger.conn.last_insert_rowid();
        assert!(
            ledger
                .reconcile_op_content_identity(rolled_back)
                .expect("reconcile inside caller transaction")
                .is_some()
        );
        ledger
            .rollback()
            .expect("roll back source and sidecar together");
        assert_eq!(ledger.op(rolled_back).unwrap(), None);
        assert_eq!(ledger.op_content_identity(rolled_back).unwrap(), None);
    }

    #[test]
    fn pre_v19_tune_writer_insert_and_value_update_reconcile_atomically() {
        let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
        let kernel = "old-writer-kernel";
        let shape = "old-writer-shape";
        let machine = b"old-writer-machine";
        ledger
            .conn
            .prepare(
                "INSERT INTO tune(kernel, shape_class, machine, params, measured)
                 VALUES (?1, ?2, ?3, ?4, ?5)",
            )
            .expect("prepare old-writer cache insert")
            .execute_with_params(&[
                text_param(kernel),
                text_param(shape),
                blob_param(machine),
                text_param(r#"{"tile":64}"#),
                text_param(r#"{"gflops":90}"#),
            ])
            .expect("insert compatibility cache row without v19 sidecar");
        assert!(matches!(
            ledger.tune_content_identity(kernel, shape, machine),
            Err(LedgerError::TuneCorrupt { detail, .. })
                if detail.contains("no typed content-identity sidecar")
        ));

        let inserted = ledger
            .reconcile_tune_content_identity(kernel, shape, machine)
            .expect("reconcile old-writer cache insert")
            .expect("compatibility cache row exists");
        assert_eq!(
            inserted.params_content_id(),
            ContentId::of_bytes(br#"{"tile":64}"#)
        );

        ledger
            .conn
            .prepare(
                "UPDATE tune SET params = ?1, measured = ?2
                 WHERE kernel = ?3 AND shape_class = ?4 AND machine = ?5",
            )
            .expect("prepare old-writer cache value update")
            .execute_with_params(&[
                text_param(r#"{"tile":96}"#),
                text_param(r#"{"gflops":103}"#),
                text_param(kernel),
                text_param(shape),
                blob_param(machine),
            ])
            .expect("update compatibility cache values without v19 sidecar update");
        assert!(matches!(
            ledger.tune_content_identity(kernel, shape, machine),
            Err(LedgerError::TuneCorrupt { detail, .. })
                if detail.contains("independently re-hashed source bytes")
        ));
        let updated = ledger
            .reconcile_tune_content_identity(kernel, shape, machine)
            .expect("reconcile old-writer cache value update")
            .expect("compatibility cache row remains present");
        assert_eq!(updated.kernel_content_id(), inserted.kernel_content_id());
        assert_eq!(updated.machine_content_id(), inserted.machine_content_id());
        assert_eq!(
            updated.params_content_id(),
            ContentId::of_bytes(br#"{"tile":96}"#)
        );
        assert_ne!(updated.params_content_id(), inserted.params_content_id());

        ledger
            .begin()
            .expect("caller-owned cache reconciliation transaction");
        ledger
            .conn
            .prepare(
                "UPDATE tune SET params = ?1, measured = ?2
                 WHERE kernel = ?3 AND shape_class = ?4 AND machine = ?5",
            )
            .expect("prepare transactional old-writer cache update")
            .execute_with_params(&[
                text_param(r#"{"tile":128}"#),
                text_param(r#"{"gflops":111}"#),
                text_param(kernel),
                text_param(shape),
                blob_param(machine),
            ])
            .expect("update compatibility values inside caller transaction");
        assert_ne!(
            ledger
                .reconcile_tune_content_identity(kernel, shape, machine)
                .expect("reconcile inside caller transaction")
                .expect("cache row remains present")
                .params_content_id(),
            updated.params_content_id()
        );
        ledger
            .rollback()
            .expect("roll back source and sidecar update together");
        assert_eq!(
            ledger
                .tune_content_identity(kernel, shape, machine)
                .unwrap(),
            Some(updated)
        );
        assert_eq!(
            ledger
                .reconcile_tune_content_identity(kernel, shape, machine)
                .expect("exact cache reconciliation retry"),
            Some(updated)
        );
    }

    #[test]
    fn v18_backfills_frozen_operation_fields_and_replays_a_stale_marker() {
        for preapply_v18 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
            let explicits = crate::FiveExplicits {
                seed: b"v17-operation-seed",
                versions: r#"{"version":17}"#,
                budget: r#"{"memory":2048}"#,
                capability: r#"{"cores":1}"#,
            };
            let op = ledger
                .begin_op(
                    Some(b"v17-operation-session"),
                    r#"{"op":"v17-fixture"}"#,
                    &explicits,
                    1,
                )
                .expect("record operation before recreating v17 fixture");
            drop_v18_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 17")
                .expect("mark genuine v17 fixture");
            if preapply_v18 {
                for ddl in schema::V18 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v18 migration batch");
                }
            }

            ledger
                .migrate_from_observed_version(17)
                .expect("backfill and authenticate exact v17 operation fields");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            let identity = ledger
                .op_content_identity(op)
                .unwrap()
                .expect("backfilled operation identity exists");
            assert_eq!(
                identity.session_content_id(),
                Some(ContentId::of_bytes(b"v17-operation-session"))
            );
            assert_eq!(
                identity.seed_content_id(),
                ContentId::of_bytes(explicits.seed)
            );
            assert_eq!(ledger.table_count("op_content_identities").unwrap(), 1);
        }
    }

    #[test]
    fn v18_malformed_source_rolls_back_objects_rows_and_marker() {
        let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
        let explicits = crate::FiveExplicits {
            seed: b"bounded-seed",
            versions: "{}",
            budget: "{}",
            capability: "{}",
        };
        let op = ledger
            .begin_op(None, "{}", &explicits, 1)
            .expect("record operation before recreating v17 fixture");
        drop_v18_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 17")
            .expect("mark genuine v17 fixture");
        ledger
            .conn
            .prepare("UPDATE ops SET seed = ?1 WHERE id = ?2")
            .expect("prepare oversized frozen-field corruption")
            .execute_with_params(&[
                SqliteValue::Blob(vec![0xA5; crate::MAX_OP_SEED_BYTES + 1].into()),
                SqliteValue::Integer(op),
            ])
            .expect("inject oversized frozen-field corruption");

        assert!(matches!(
            ledger.migrate_from_observed_version(17),
            Err(LedgerError::OpCorrupt { op: rejected, .. }) if rejected == op
        ));
        assert_eq!(ledger.schema_version().unwrap(), 17);
        let objects = ledger
            .conn
            .query(
                "SELECT name FROM sqlite_master
                 WHERE name = 'op_content_identities' LIMIT 1",
            )
            .expect("inspect rolled-back v18 schema");
        assert!(objects.is_empty());
    }

    #[test]
    fn v19_backfills_tune_fields_and_replays_a_stale_marker() {
        for preapply_v19 in [false, true] {
            let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
            ledger
                .tune_put(
                    "v18-kernel",
                    "v18-shape",
                    b"v18-machine",
                    r#"{"tile":64}"#,
                    r#"{"gflops":91}"#,
                )
                .expect("record tune row before recreating v18 fixture");
            drop_v19_objects(&ledger);
            ledger
                .conn
                .execute("PRAGMA user_version = 18")
                .expect("mark genuine v18 fixture");
            if preapply_v19 {
                for ddl in schema::V19 {
                    ledger
                        .conn
                        .execute(ddl)
                        .expect("preapply exact v19 migration batch");
                }
            }

            ledger
                .migrate_from_observed_version(18)
                .expect("backfill and authenticate exact v18 tune fields");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            let identity = ledger
                .tune_content_identity("v18-kernel", "v18-shape", b"v18-machine")
                .unwrap()
                .expect("backfilled tune identity exists");
            assert_eq!(
                identity.kernel_content_id(),
                ContentId::of_bytes(b"v18-kernel")
            );
            assert_eq!(
                identity.params_content_id(),
                ContentId::of_bytes(br#"{"tile":64}"#)
            );
            assert_eq!(ledger.table_count("tune_content_identities").unwrap(), 1);
        }
    }

    #[test]
    fn v19_malformed_tune_source_rolls_back_objects_rows_and_marker() {
        let ledger = Ledger::open(":memory:").expect("fresh v19 ledger");
        ledger
            .tune_put("v18-kernel", "v18-shape", b"v18-machine", "{}", "{}")
            .expect("record tune row before recreating v18 fixture");
        drop_v19_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 18")
            .expect("mark genuine v18 fixture");
        let oversized = format!("\"{}\"", "x".repeat(crate::MAX_TUNE_PARAMS_BYTES));
        ledger
            .conn
            .prepare(
                "UPDATE tune SET params = ?1
                 WHERE kernel = ?2 AND shape_class = ?3 AND machine = ?4",
            )
            .expect("prepare oversized cache-source corruption")
            .execute_with_params(&[
                text_param(&oversized),
                text_param("v18-kernel"),
                text_param("v18-shape"),
                blob_param(b"v18-machine"),
            ])
            .expect("inject oversized cache-source corruption");

        assert!(matches!(
            ledger.migrate_from_observed_version(18),
            Err(LedgerError::TuneCorrupt { .. })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 18);
        let objects = ledger
            .conn
            .query(
                "SELECT name FROM sqlite_master
                 WHERE name = 'tune_content_identities' LIMIT 1",
            )
            .expect("inspect rolled-back v19 schema");
        assert!(objects.is_empty());
    }

    #[test]
    fn migration_ladder_ends_with_v19() {
        assert_eq!(
            schema::MIGRATIONS.len(),
            usize::try_from(SCHEMA_VERSION).unwrap()
        );
        assert_eq!(schema::MIGRATIONS.get(12), Some(&schema::V13));
        assert_eq!(schema::MIGRATIONS.get(13), Some(&schema::V14));
        assert_eq!(schema::MIGRATIONS.get(14), Some(&schema::V15));
        assert_eq!(schema::MIGRATIONS.get(15), Some(&schema::V16));
        assert_eq!(schema::MIGRATIONS.get(16), Some(&schema::V17));
        assert_eq!(schema::MIGRATIONS.get(17), Some(&schema::V18));
        assert_eq!(schema::MIGRATIONS.last(), Some(&schema::V19));
    }
}

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

use std::str;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, ContentId, Field, FieldSpec,
    IdentityAuditRecord, IdentityReceipt, IdentityRole, NoClaimState, SchemaId, SemanticId,
    StrongIdentity, TrustState, WireType, legacy::LegacyProvenanceV1,
};
use fsqlite::SqliteValue;

use crate::{ContentHash, EdgeRole, Ledger, LedgerError, blob_param, now_wall_ns, sql_err};

/// Schema version of one artifact compatibility-hash to typed-content-ID row.
pub const ARTIFACT_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Schema version of one lineage-edge typed-content-ID companion row.
pub const EDGE_CONTENT_IDENTITY_ROW_VERSION: u32 = 1;
/// Maximum receipt IDs exposed by one artifact-semantic candidate lookup.
pub const MAX_ARTIFACT_SEMANTIC_BINDING_CANDIDATES: usize = 256;

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
    "schema_constants=IDENTITY_MIGRATION_RECEIPT_VERSION,IDENTITY_MIGRATION_RECEIPT_DOMAIN,MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,MAX_IDENTITY_MIGRATION_RULE_BYTES,MAX_IDENTITY_MIGRATION_DOMAIN_BYTES,MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES,MAX_IDENTITY_MIGRATION_CONTEXT_BYTES,crates/fs-ledger/src/schema.rs#V13",
    "schema_functions=derive_receipt_id,receipt_body_from_claim,validate_receipt_body,IdentityMigrationReceipt::typed_semantic_id,Ledger::record_identity_migration,Ledger::identity_migration_receipt,Ledger::identity_migration_candidates,crates/fs-blake3/src/identity.rs#CanonicalEncoder::finish",
    "schema_dependencies=fs-blake3:canonical-identity-frame,fs-blake3:schema-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=IdentityMigrationBody",
    "source_fields=IdentityMigrationBody.legacy_content_id:semantic,IdentityMigrationBody.legacy_fnv:semantic,IdentityMigrationBody.canonical_content_id:semantic,IdentityMigrationBody.semantic_rule:semantic,IdentityMigrationBody.semantic_id:semantic,IdentityMigrationBody.identity_role:semantic,IdentityMigrationBody.identity_domain:semantic,IdentityMigrationBody.identity_schema_name:semantic,IdentityMigrationBody.identity_schema_id:semantic,IdentityMigrationBody.identity_schema_version:semantic,IdentityMigrationBody.identity_context:semantic,IdentityMigrationBody.canonical_preimage_id:semantic,IdentityMigrationBody.canonical_frame_bytes:semantic,IdentityMigrationBody.field_count:semantic,IdentityMigrationBody.collection_items:semantic,IdentityMigrationBody.limits:semantic,IdentityMigrationBody.trust_state:semantic,IdentityMigrationBody.anchor_content_id:semantic,IdentityMigrationBody.verifier_id:semantic,IdentityMigrationBody.key_policy_id:semantic,IdentityMigrationBody.no_claim_state:semantic,IdentityMigrationBody.legacy_bytes:derived:bound-by-content-id-and-byte-count,IdentityMigrationBody.canonical_bytes:derived:bound-by-content-id-and-byte-count",
    "source_bindings=IdentityMigrationBody.legacy_bytes>legacy-content-id+legacy-byte-count,IdentityMigrationBody.legacy_fnv>legacy-fnv-le-u64,IdentityMigrationBody.canonical_bytes>canonical-content-id+canonical-byte-count,IdentityMigrationBody.semantic_rule>semantic-rule,IdentityMigrationBody.semantic_id>semantic-id,IdentityMigrationBody.identity_role>identity-role,IdentityMigrationBody.identity_domain>identity-domain,IdentityMigrationBody.identity_schema_name>identity-schema-name,IdentityMigrationBody.identity_schema_id>identity-schema-id,IdentityMigrationBody.identity_schema_version>identity-schema-version,IdentityMigrationBody.identity_context>identity-context,IdentityMigrationBody.canonical_preimage_id>canonical-preimage-content-id,IdentityMigrationBody.canonical_frame_bytes>canonical-frame-bytes,IdentityMigrationBody.field_count>field-count,IdentityMigrationBody.collection_items>collection-items,IdentityMigrationBody.limits>max-canonical-bytes+max-field-bytes+max-fields+max-collection-items+cancellation-poll-bytes,IdentityMigrationBody.trust_state>trust-state,IdentityMigrationBody.anchor_content_id>anchor-content-id,IdentityMigrationBody.verifier_id>verifier-id,IdentityMigrationBody.key_policy_id>key-policy-id,IdentityMigrationBody.no_claim_state>no-claim-state",
    "external_semantic_fields=receipt-schema-domain,receipt-schema-version,canonical-field-order",
    "semantic_fields=receipt-schema-domain,receipt-schema-version,canonical-field-order,legacy-content-id,legacy-byte-count,legacy-fnv-le-u64,canonical-content-id,canonical-byte-count,semantic-rule,identity-role,semantic-id,identity-domain,identity-schema-name,identity-schema-id,identity-schema-version,identity-context,canonical-preimage-content-id,canonical-frame-bytes,field-count,collection-items,max-canonical-bytes,max-field-bytes,max-fields,max-collection-items,cancellation-poll-bytes,trust-state,anchor-content-id,verifier-id,key-policy-id,no-claim-state",
    "excluded_fields=exact-legacy-bytes:bound-by-content-id,exact-canonical-bytes:bound-by-content-id,created-at:provenance-envelope-only",
    "consumers=Ledger::record_identity_migration,Ledger::identity_migration_receipt,Ledger::identity_migration_candidates,IdentityMigrationReceipt::typed_semantic_id",
    "mutations=receipt-schema-domain:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,receipt-schema-version:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-field-order:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-byte-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,legacy-fnv-le-u64:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-byte-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,semantic-rule:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-role:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,semantic-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-domain:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-name:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-schema-version:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,identity-context:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-preimage-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,canonical-frame-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,field-count:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,collection-items:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-canonical-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-field-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-fields:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,max-collection-items:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,cancellation-poll-bytes:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,trust-state:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,anchor-content-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,verifier-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,key-policy-id:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state,no-claim-state:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state",
    "nonsemantic_mutations=created-at:crates/fs-ledger/tests/identity_migration.rs#receipt_identity_binds_exact_bytes_schema_and_audit_state",
    "field_guard=classify_identity_migration_receipt_fields",
    "transport_guard=Ledger::identity_migration_receipt",
    "version_guard=crates/fs-ledger/tests/identity_migration.rs#typed_projection_refuses_a_different_schema",
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
    pub fn artifact_semantic_binding(
        &self,
        artifact_hash: &ContentHash,
        receipt_id: IdentityMigrationReceiptId,
    ) -> Result<Option<ArtifactSemanticBinding>, LedgerError> {
        self.stored_artifact_semantic_binding(artifact_hash, receipt_id)
    }

    /// Return a bounded deterministic receipt-ID prefix for one artifact.
    /// This lookup never chooses or promotes a semantic interpretation.
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
    fn genuine_v12_and_stale_v13_markers_migrate_through_v16() {
        let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
        let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
            let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
        let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
            let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
        let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
            let ledger = Ledger::open(":memory:").expect("fresh v16 ledger");
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
    fn migration_ladder_ends_with_v16() {
        assert_eq!(
            schema::MIGRATIONS.len(),
            usize::try_from(SCHEMA_VERSION).unwrap()
        );
        assert_eq!(schema::MIGRATIONS.get(12), Some(&schema::V13));
        assert_eq!(schema::MIGRATIONS.get(13), Some(&schema::V14));
        assert_eq!(schema::MIGRATIONS.get(14), Some(&schema::V15));
        assert_eq!(schema::MIGRATIONS.last(), Some(&schema::V16));
    }
}
